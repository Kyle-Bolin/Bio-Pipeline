// Client for the cellfactory REST API (src/cellfactory/api.py):
//   GET  /api/health, GET /api/tools
//   POST /api/tools/<name>      body = the tool's kwargs, sync
//   POST /api/jobs              body = {tool, args}, then poll GET /api/jobs/<id>
// Tool names are typed as ToolName so tests/test_web.py can check each one still exists.

export const API = "/api";

/** Every tool the page calls by name. Keep in sync with TOOLS_USED in src/cellfactory/web/__init__.py. */
export const TOOL_NAMES = [
  "list_models", "describe_model", "list_pathways", "search_metabolites",
  "simulate_production", "production_envelope", "search_knockouts",
  "search_knockouts_moma", "find_precursors", "check_product_fate", "check_model",
  "compare_conditions", "flux_distribution", "check_thermodynamics", "qa_model",
] as const;
export type ToolName = (typeof TOOL_NAMES)[number];

/** Models small enough to answer synchronously; anything else goes through the job queue. */
const BUNDLED_FAST = new Set(["e_coli_core"]);
/** Tools that take seconds to minutes. The API's `slow` flag adds to this at start-up. */
export const SLOW = new Set<string>(["search_knockouts", "search_knockouts_moma", "find_precursors", "check_model",
  "check_product_fate", "check_thermodynamics", "qa_model"]);

export class ApiError extends Error {
  status?: number;
}

export type Args = Record<string, unknown>;
export type JobStatus = "queued" | "running" | "done" | "failed";

export async function request<T = any>(method: string, path: string, body?: unknown): Promise<T> {
  let res: Response;
  try {
    res = await fetch(API + path, {
      method,
      headers: body === undefined ? {} : { "Content-Type": "application/json" },
      body: body === undefined ? undefined : JSON.stringify(body),
    });
  } catch {
    throw new ApiError("Can't reach the API. Is `cellfactory-api` running?");
  }
  let data: any = null;
  try { data = await res.json(); } catch { /* empty or non-JSON body */ }
  if (!res.ok) {
    const msg = (data && (data.error || data.detail)) || `${res.status} ${res.statusText}`;
    const err = new ApiError(typeof msg === "string" ? msg : JSON.stringify(msg));
    err.status = res.status;
    throw err;
  }
  return data as T;
}

const JOB_TIMEOUT_MS = 20 * 60 * 1000; // a genome-scale MOMA search can take several minutes
const sleep = (ms: number) => new Promise((r) => setTimeout(r, ms));

async function runJob<T>(tool: string, args: Args, onStatus?: (s: JobStatus) => void): Promise<T> {
  const job = await request<{ id: string }>("POST", "/jobs", { tool, args });
  let delay = 400;
  const deadline = Date.now() + JOB_TIMEOUT_MS;
  for (;;) {
    if (Date.now() > deadline) {
      throw new ApiError(`Still not finished after ${JOB_TIMEOUT_MS / 60000} minutes, so the page stopped waiting. ` +
        `The job may still be running: check ${API}/jobs/${job.id}.`);
    }
    await sleep(delay);
    delay = Math.min(delay * 1.5, 2000);
    const j = await request<{ status: JobStatus; result?: T; error?: string }>("GET", `/jobs/${encodeURIComponent(job.id)}`);
    onStatus?.(j.status);
    if (j.status === "done") return j.result as T;
    if (j.status === "failed") throw new ApiError(j.error || "The job failed.");
  }
}

let jobsSupported = false;
export async function probeJobs(): Promise<boolean> {
  try { await request("GET", "/jobs"); jobsSupported = true; } catch { jobsSupported = false; }
  return jobsSupported;
}

/** Call a tool; slow tools (or anything on a genome-scale model) go through /api/jobs. */
export async function callTool<T = any>(tool: string, args: Args = {}, onStatus?: (s: JobStatus) => void): Promise<T> {
  const model = typeof args.model === "string" ? args.model : "";
  const heavy = SLOW.has(tool) || (model !== "" && !BUNDLED_FAST.has(model));
  if (heavy && jobsSupported) return runJob<T>(tool, args, onStatus);
  return request<T>("POST", `/tools/${encodeURIComponent(tool)}`, args);
}

/** Typed entry point for the named views: the name must be in TOOL_NAMES. */
export function call<T = any>(tool: ToolName, args: Args = {}, onStatus?: (s: JobStatus) => void): Promise<T> {
  return callTool<T>(tool, args, onStatus);
}

export const isBundledFast = (model: string) => BUNDLED_FAST.has(model);

// ---------- tool specs from GET /api/tools ----------
export interface ParamSchema {
  type?: string;
  anyOf?: { type?: string }[];
  default?: unknown;
  enum?: unknown[];
  description?: string;
}
export interface ToolSpec {
  name: string;
  description?: string;
  slow?: boolean;
  parameters?: { properties?: Record<string, ParamSchema>; required?: string[] };
}
