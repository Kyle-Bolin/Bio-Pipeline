"""REST API: the same tools as the MCP server, over HTTP for the web front end and scripts.

Run ``cellfactory-api`` (needs the ``[api]`` extra) and open http://127.0.0.1:8000/docs
for interactive documentation. Every function in ``tools.TOOLS`` becomes an endpoint,
so a tool added there shows up here and in the MCP server with no extra code.

    GET  /api/health                health check and package version
    GET  /api/tools                 every tool: name, description, JSON schema of its parameters
    POST /api/tools/{name}          call a tool; JSON body = its keyword arguments
    POST /api/jobs                  {"tool": name, "args": {...}}: run a slow tool in the background
    GET  /api/jobs                  recent jobs, newest first (without results)
    GET  /api/jobs/{id}             one job, with its result once done

Errors come back as ``{"error": "<message>"}``: 400 for input the caller can fix
(a ``ToolInputError``, or a missing, unknown or mistyped argument), 404 for an
unknown tool or job, 500 for a bug. If ``web/static/index.html`` exists in the
package it is served at ``/``. Set ``CELLFACTORY_CORS_ORIGINS`` to a comma-separated
list of origins to allow a front end served from somewhere else.
"""

from __future__ import annotations

import inspect
import json
import logging
import os
import threading
import typing
import uuid
from collections import OrderedDict
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Any, Callable

from fastapi import FastAPI, Request
from fastapi.concurrency import run_in_threadpool
from fastapi.exceptions import RequestValidationError
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel, ConfigDict, Field, ValidationError, create_model

from .tools import TOOLS, ToolInputError

log = logging.getLogger(__name__)

STATIC_DIR = Path(__file__).parent / "web" / "static"
MAX_JOBS = 100  # finished jobs beyond this are forgotten, oldest first
JOB_WORKERS = int(os.environ.get("CELLFACTORY_JOB_WORKERS", "2"))

MAX_PENDING = int(os.environ.get("CELLFACTORY_MAX_PENDING_JOBS", "20"))  # queued + running; more get a 429
# Tools that can take minutes on a genome-scale model: call them through /api/jobs, not synchronously.
SLOW_TOOLS = {"search_knockouts", "search_knockouts_moma", "find_precursors"}

try:
    __version__ = version("cellfactory")
except PackageNotFoundError:  # running from a source tree that was never installed
    __version__ = "unknown"


class APIError(Exception):
    def __init__(self, status: int, message: str):
        super().__init__(message)
        self.status = status
        self.message = message


def _params_model(fn: Callable) -> type[BaseModel]:
    """A pydantic model of ``fn``'s keyword arguments, so FastAPI validates and documents them."""
    try:
        hints = typing.get_type_hints(fn)  # tools.py uses `from __future__ import annotations`
    except Exception as err:
        raise TypeError(f"Cannot resolve the argument types of tool {fn.__name__!r}: {err}") from err
    fields: dict[str, Any] = {}
    for p in inspect.signature(fn).parameters.values():
        if p.name not in hints:
            raise TypeError(f"Tool {fn.__name__!r} argument {p.name!r} needs a type annotation.")
        annotation = hints[p.name]
        default = ... if p.default is inspect.Parameter.empty else p.default
        fields[p.name] = (annotation, Field(default))
    name = "".join(part.title() for part in fn.__name__.split("_")) + "Args"
    return create_model(name, __config__=ConfigDict(extra="forbid"), **fields)


def _description(fn: Callable) -> str:
    return inspect.cleandoc(fn.__doc__ or "")


def _format_validation(err: ValidationError | RequestValidationError) -> str:
    parts = []
    for e in err.errors():
        loc = [str(x) for x in e.get("loc", ()) if x not in ("body",)]
        where = ".".join(loc)
        msg = "unknown argument" if e.get("type") == "extra_forbidden" else e.get("msg", "invalid")
        parts.append(f"{where}: {msg}" if where else msg)
    return "Invalid arguments: " + "; ".join(parts) + "."


class ToolRegistry:
    def __init__(self, tools: list[Callable]):
        self.fns = {fn.__name__: fn for fn in tools}
        self.models = {name: _params_model(fn) for name, fn in self.fns.items()}

    def get(self, name: str) -> Callable:
        if name not in self.fns:
            raise APIError(404, f"Unknown tool {name!r}. GET /api/tools lists them.")
        return self.fns[name]

    def validate(self, name: str, args: dict | None) -> dict:
        self.get(name)
        if args is not None and not isinstance(args, dict):
            raise APIError(400, "Arguments must be a JSON object.")
        try:
            return self.models[name].model_validate(args or {}).model_dump()
        except ValidationError as err:
            raise APIError(400, _format_validation(err)) from err

    def describe(self) -> list[dict]:
        out = []
        for name, fn in self.fns.items():
            schema = self.models[name].model_json_schema()
            schema.pop("title", None)
            schema.pop("additionalProperties", None)
            schema.setdefault("required", [])
            out.append({"name": name, "description": _description(fn), "parameters": schema,
                        "slow": name in SLOW_TOOLS})
        return out


def _call(fn: Callable, kwargs: dict) -> dict:
    try:
        return fn(**kwargs)
    except ToolInputError as err:
        raise APIError(400, str(err)) from err


def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


class JobStore:
    """Background runs of slow tools, kept in memory. One store per app."""

    def __init__(self, workers: int = JOB_WORKERS, keep: int = MAX_JOBS, max_pending: int = MAX_PENDING):
        self._max_pending = max_pending
        self._jobs: OrderedDict[str, dict] = OrderedDict()
        self._lock = threading.Lock()
        self._pool = ThreadPoolExecutor(max_workers=max(1, workers), thread_name_prefix="cellfactory-job")
        self._keep = keep

    def submit(self, tool: str, fn: Callable, kwargs: dict) -> dict:
        job = {"id": uuid.uuid4().hex, "tool": tool, "args": kwargs, "status": "queued",
               "result": None, "error": None, "created_at": _now(), "started_at": None, "finished_at": None}
        with self._lock:
            pending = sum(j["status"] in ("queued", "running") for j in self._jobs.values())
            if pending >= self._max_pending:
                raise APIError(429, f"{pending} jobs are already queued or running; wait for some to finish.")
            self._jobs[job["id"]] = job
            self._prune()
            snapshot = dict(job)
        self._pool.submit(self._run, job, fn, kwargs)
        return snapshot

    def _run(self, job: dict, fn: Callable, kwargs: dict) -> None:
        with self._lock:
            job.update(status="running", started_at=_now())
        try:
            result, status, error = _call(fn, kwargs), "done", None
        except APIError as err:
            result, status, error = None, "failed", err.message
        except Exception as err:  # a bug in a tool; report it rather than leave the job running forever
            log.exception("job %s (%s) crashed", job["id"], job["tool"])
            result, status, error = None, "failed", f"Internal error: {type(err).__name__}: {err}"
        with self._lock:
            job.update(status=status, result=result, error=error, finished_at=_now())

    def _prune(self) -> None:
        finished = [k for k, j in self._jobs.items() if j["status"] in ("done", "failed")]
        while len(self._jobs) > self._keep and finished:
            del self._jobs[finished.pop(0)]

    def get(self, job_id: str) -> dict:
        with self._lock:
            if job_id not in self._jobs:
                raise APIError(404, f"No job {job_id!r}; it may have expired.")
            return dict(self._jobs[job_id])

    def list(self) -> list[dict]:
        with self._lock:
            return [{k: v for k, v in j.items() if k != "result"} for j in reversed(self._jobs.values())]

    def shutdown(self) -> None:
        self._pool.shutdown(wait=False, cancel_futures=True)


class JobRequest(BaseModel):
    model_config = ConfigDict(extra="forbid")
    tool: str = Field(description="Tool name, as listed by GET /api/tools")
    args: dict[str, Any] = Field(default_factory=dict, description="The tool's keyword arguments")


def _add_tool_route(app: FastAPI, registry: ToolRegistry, name: str) -> None:
    fn, model = registry.fns[name], registry.models[name]

    async def endpoint(request: Request) -> dict:
        # Read the body ourselves so an empty one is fine for tools whose arguments all have
        # defaults; the schema is still published to /docs through openapi_extra below.
        raw = await request.body()
        try:
            args = json.loads(raw) if raw.strip() else {}
        except json.JSONDecodeError as err:
            raise APIError(400, f"Body is not valid JSON: {err}.") from err
        kwargs = registry.validate(name, args)
        return await run_in_threadpool(_call, fn, kwargs)

    summary = _description(fn).split("\n\n")[0].replace("\n", " ")
    description = _description(fn)
    if name in SLOW_TOOLS:
        description = ("**Slow** (minutes on a genome-scale model): prefer `POST /api/jobs`.\n\n" + description)
    app.post(f"/api/tools/{name}", name=name, summary=summary, description=description,
             tags=["tools"], openapi_extra={"requestBody": {
                 "required": bool(model.model_json_schema().get("required")),
                 "content": {"application/json": {"schema": model.model_json_schema()}}}})(endpoint)


def create_app(tools: list[Callable] | None = None, static_dir: Path | None = STATIC_DIR) -> FastAPI:
    registry = ToolRegistry(TOOLS if tools is None else tools)
    jobs = JobStore()
    app = FastAPI(title="cellfactory", version=__version__,
                  description="Flux balance analysis of microbes turning a substrate into a product. "
                              "The same tools as the cellfactory MCP server.")
    app.state.registry = registry
    app.state.jobs = jobs
    app.router.add_event_handler("shutdown", jobs.shutdown)

    origins = [o.strip() for o in os.environ.get("CELLFACTORY_CORS_ORIGINS", "").split(",") if o.strip()]
    if origins:
        app.add_middleware(CORSMiddleware, allow_origins=origins, allow_methods=["GET", "POST"],
                           allow_headers=["*"])

    @app.exception_handler(APIError)
    async def _api_error(request: Request, err: APIError):
        return JSONResponse({"error": err.message}, status_code=err.status)

    @app.exception_handler(RequestValidationError)
    async def _validation_error(request: Request, err: RequestValidationError):
        return JSONResponse({"error": _format_validation(err)}, status_code=400)

    @app.exception_handler(Exception)
    async def _crash(request: Request, err: Exception):
        log.exception("unhandled error on %s", request.url.path)
        return JSONResponse({"error": f"Internal error: {type(err).__name__}: {err}"}, status_code=500)

    @app.get("/api/health", tags=["meta"])
    async def health() -> dict:
        return {"status": "ok", "version": __version__}

    @app.get("/api/tools", tags=["meta"], summary="Every tool, with a JSON schema of its parameters")
    async def list_tools() -> list[dict]:
        return registry.describe()

    for name in registry.fns:
        _add_tool_route(app, registry, name)

    @app.post("/api/tools/{name}", include_in_schema=False)
    def unknown_tool(name: str):
        registry.get(name)  # always raises: every known tool has its own route above

    @app.post("/api/jobs", status_code=202, tags=["jobs"],
              summary="Run a tool in the background; poll GET /api/jobs/{id} for the result")
    async def submit_job(req: JobRequest) -> dict:
        kwargs = registry.validate(req.tool, req.args)
        return jobs.submit(req.tool, registry.get(req.tool), kwargs)

    @app.get("/api/jobs", tags=["jobs"], summary="Recent jobs, newest first, without results")
    async def list_jobs() -> dict:
        return {"jobs": jobs.list()}

    @app.get("/api/jobs/{job_id}", tags=["jobs"], summary="One job, with its result once done")
    async def get_job(job_id: str) -> dict:
        return jobs.get(job_id)

    if static_dir is not None and (Path(static_dir) / "index.html").is_file():
        app.mount("/", StaticFiles(directory=static_dir, html=True), name="web")

    return app


def main(argv: list[str] | None = None) -> None:
    import argparse

    import uvicorn

    parser = argparse.ArgumentParser(prog="cellfactory-api", description="Serve the cellfactory REST API.")
    parser.add_argument("--host", default="127.0.0.1", help="interface to bind (default 127.0.0.1)")
    parser.add_argument("--port", type=int, default=8000)
    parser.add_argument("--reload", action="store_true", help="restart on code changes (development)")
    args = parser.parse_args(argv)
    logging.getLogger("cobra").setLevel(logging.WARNING)
    if args.reload:
        uvicorn.run("cellfactory.api:create_app", factory=True, host=args.host, port=args.port, reload=True)
    else:
        uvicorn.run(create_app(), host=args.host, port=args.port)


if __name__ == "__main__":
    main()
