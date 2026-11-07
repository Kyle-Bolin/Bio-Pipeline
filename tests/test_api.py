"""The REST API (api.py): every tool over HTTP, errors as {"error": ...}, background jobs."""

import threading
import time

import pytest

pytest.importorskip("fastapi")
pytest.importorskip("httpx")

from fastapi.testclient import TestClient  # noqa: E402

from cellfactory import api, tools  # noqa: E402


@pytest.fixture(scope="module")
def client():
    with TestClient(api.create_app(static_dir=None)) as c:
        yield c


def _wait(client, job_id, timeout=60):
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        job = client.get(f"/api/jobs/{job_id}").json()
        if job["status"] in ("done", "failed"):
            return job
        time.sleep(0.1)
    raise AssertionError(f"job {job_id} did not finish")


def test_health(client):
    assert client.get("/api/health").json()["status"] == "ok"


def test_lists_every_tool_with_a_parameter_schema(client):
    listed = {t["name"]: t for t in client.get("/api/tools").json()}
    assert set(listed) == {fn.__name__ for fn in tools.TOOLS}
    sim = listed["simulate_production"]
    assert sim["description"].startswith("Feed ``substrate``")
    assert sim["parameters"]["required"] == ["product"]
    props = sim["parameters"]["properties"]
    assert props["aerobic"] == {"default": True, "title": "Aerobic", "type": "boolean"}
    assert props["uptake"]["type"] == "number"
    assert listed["list_models"]["parameters"]["properties"] == {}
    assert listed["search_knockouts"]["slow"] and not listed["list_models"]["slow"]
    assert api.SLOW_TOOLS <= set(listed)


def test_calls_a_tool(client):
    r = client.post("/api/tools/simulate_production", json={"product": "ethanol", "aerobic": False})
    assert r.status_code == 200, r.text
    out = r.json()
    assert out["product"] == "EX_etoh_e"
    assert out["production"]["max_theoretical_yield"] == pytest.approx(2.0, abs=1e-3)


def test_tool_without_body_uses_defaults(client):
    r = client.post("/api/tools/list_models")
    assert r.status_code == 200
    assert {m["id"] for m in r.json()["bundled"]} >= {"e_coli_core"}
    assert client.post("/api/tools/describe_model").json()["model"] == "e_coli_core"


@pytest.mark.parametrize("name, body, status, message", [
    ("simulate_production", {"product": "unobtainium"}, 400, "not a metabolite"),
    ("simulate_production", {}, 400, "product: Field required"),
    ("simulate_production", {"product": "ethanol", "colour": "red"}, 400, "colour: unknown argument"),
    ("simulate_production", {"product": "ethanol", "uptake": "lots"}, 400, "uptake:"),
    ("search_knockouts", {"product": "ethanol", "max_knockouts": 0}, 400, "between 1 and 10"),
    ("no_such_tool", {}, 404, "Unknown tool"),
])
def test_errors_come_back_as_json(client, name, body, status, message):
    r = client.post(f"/api/tools/{name}", json=body)
    assert r.status_code == status, r.text
    assert message in r.json()["error"]


def test_job_runs_a_tool_in_the_background(client):
    r = client.post("/api/jobs", json={"tool": "simulate_production",
                                       "args": {"product": "ethanol", "aerobic": False}})
    assert r.status_code == 202, r.text
    job = r.json()
    assert job["status"] == "queued" and job["result"] is None
    done = _wait(client, job["id"])
    assert done["status"] == "done" and done["error"] is None
    assert done["result"]["product"] == "EX_etoh_e"
    assert done["args"]["uptake"] == 10.0  # defaults are filled in
    listed = client.get("/api/jobs").json()["jobs"]
    assert listed[0]["id"] == job["id"] and "result" not in listed[0]


def test_job_reports_tool_errors(client):
    job = client.post("/api/jobs", json={"tool": "simulate_production", "args": {"product": "unobtainium"}}).json()
    done = _wait(client, job["id"])
    assert done["status"] == "failed" and "not a metabolite" in done["error"]


def test_job_input_is_checked_before_queueing(client):
    assert client.post("/api/jobs", json={"tool": "nope"}).status_code == 404
    r = client.post("/api/jobs", json={"tool": "simulate_production", "args": {"uptake": 5}})
    assert r.status_code == 400 and "product" in r.json()["error"]
    assert client.get("/api/jobs/doesnotexist").status_code == 404


def test_tool_crash_is_a_500_not_a_hang():
    def explode() -> dict:
        """Always fails."""
        raise RuntimeError("boom")

    with TestClient(api.create_app(tools=[explode], static_dir=None), raise_server_exceptions=False) as c:
        r = c.post("/api/tools/explode")
        assert r.status_code == 500 and "boom" in r.json()["error"]
        job = c.post("/api/jobs", json={"tool": "explode"}).json()
        done = _wait(c, job["id"])
        assert done["status"] == "failed" and "boom" in done["error"]


def test_serves_the_front_end_when_present(tmp_path):
    (tmp_path / "index.html").write_text("<h1>cellfactory</h1>")
    with TestClient(api.create_app(tools=[tools.list_models], static_dir=tmp_path)) as c:
        assert "cellfactory" in c.get("/").text
        assert c.get("/api/health").status_code == 200


def test_cors_is_opt_in(monkeypatch):
    monkeypatch.setenv("CELLFACTORY_CORS_ORIGINS", "http://localhost:5173")
    with TestClient(api.create_app(tools=[tools.list_models], static_dir=None)) as c:
        r = c.get("/api/health", headers={"Origin": "http://localhost:5173"})
        assert r.headers["access-control-allow-origin"] == "http://localhost:5173"


def test_tools_with_richer_annotations_still_load():
    def optional_arg(label: "str | None" = None, count: int = 1) -> dict:
        """Echo."""
        return {"label": label, "count": count}

    with TestClient(api.create_app(tools=[optional_arg], static_dir=None)) as c:
        assert c.post("/api/tools/optional_arg", json={"count": 2}).json() == {"label": None, "count": 2}


def test_unannotated_tool_fails_naming_the_argument():
    def bare(x):
        """No types."""
        return {}

    with pytest.raises(TypeError, match="'bare' argument 'x'"):
        api.create_app(tools=[bare], static_dir=None)


def test_job_queue_is_capped():
    gate = threading.Event()

    def wait() -> dict:
        """Blocks until released."""
        gate.wait(10)
        return {}

    app = api.create_app(tools=[wait], static_dir=None)
    app.state.jobs._max_pending = 2
    with TestClient(app) as c:
        try:
            assert c.post("/api/jobs", json={"tool": "wait"}).status_code == 202
            assert c.post("/api/jobs", json={"tool": "wait"}).status_code == 202
            r = c.post("/api/jobs", json={"tool": "wait"})
            assert r.status_code == 429 and "already queued or running" in r.json()["error"]
        finally:
            gate.set()
