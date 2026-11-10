import { createContext, useCallback, useContext, useRef, useState, type ReactNode } from "react";
import { callTool, request, type Args, type RunStatus } from "./api";
import { numberProblem } from "./format";

// ---------- shared setup: every tab reads its conditions from here ----------
export interface Setup {
  model: string;
  substrate: string;
  product: string;
  uptake: number;
  aerobic: boolean;
  pathway: string;
}

export const DEFAULT_SETUP: Setup = {
  model: "e_coli_core", substrate: "glucose", product: "ethanol", uptake: 10, aerobic: true, pathway: "",
};

/** Why the shared setup can't be run yet, or null. */
export const setupProblem = (s: Setup) => numberProblem("Uptake", s.uptake, { min: 0 });

/** The condition kwargs most tools share. */
export function conditions(s: Setup) {
  return { model: s.model.trim() || "e_coli_core", substrate: s.substrate.trim() || "glucose", uptake: s.uptake, aerobic: s.aerobic };
}

interface SetupCtx { setup: Setup; update: (patch: Partial<Setup>) => void }
const SetupContext = createContext<SetupCtx | null>(null);

export function SetupProvider({ children }: { children: ReactNode }) {
  const [setup, setSetup] = useState<Setup>(DEFAULT_SETUP);
  const update = useCallback((patch: Partial<Setup>) => setSetup((s) => ({ ...s, ...patch })), []);
  return <SetupContext.Provider value={{ setup, update }}>{children}</SetupContext.Provider>;
}

export function useSetup() {
  const ctx = useContext(SetupContext);
  if (!ctx) throw new Error("useSetup outside SetupProvider");
  return ctx;
}

// ---------- run history ----------
export interface HistoryEntry { tool: string; args: Args; ok: boolean; ms: number; at: Date }
interface HistoryCtx { history: HistoryEntry[]; timed: <T>(tool: string, args: Args, onStatus?: (s: RunStatus) => void) => Promise<T> }
const HistoryContext = createContext<HistoryCtx | null>(null);

export function HistoryProvider({ children }: { children: ReactNode }) {
  const [history, setHistory] = useState<HistoryEntry[]>([]);
  const timed = useCallback(async <T,>(tool: string, args: Args, onStatus?: (s: RunStatus) => void): Promise<T> => {
    const t0 = performance.now();
    const log = (ok: boolean) => setHistory((h) => [{ tool, args, ok, ms: performance.now() - t0, at: new Date() }, ...h].slice(0, 12));
    try {
      const r = await callTool<T>(tool, args, onStatus);
      log(true);
      return r;
    } catch (err) {
      log(false);
      throw err;
    }
  }, []);
  return <HistoryContext.Provider value={{ history, timed }}>{children}</HistoryContext.Provider>;
}

export function useHistory() {
  const ctx = useContext(HistoryContext);
  if (!ctx) throw new Error("useHistory outside HistoryProvider");
  return ctx;
}

// ---------- running an analysis ----------
export type RunState<T> =
  | { kind: "idle" }
  | { kind: "running"; status: RunStatus | null }
  | { kind: "done"; data: T }
  | { kind: "error"; message: string };

/**
 * State for one Run button. `run(work)` passes `work` a status callback for job
 * polling and stores its result. A newer run supersedes an older one still in flight.
 */
export function useRun<T>() {
  const [state, setState] = useState<RunState<T>>({ kind: "idle" });
  const seq = useRef(0);
  const run = useCallback(async (work: (onStatus: (s: RunStatus) => void) => Promise<T>) => {
    const mine = ++seq.current;
    setState({ kind: "running", status: null });
    try {
      const data = await work((status) => { if (seq.current === mine) setState({ kind: "running", status }); });
      if (seq.current === mine) setState({ kind: "done", data });
    } catch (err) {
      if (seq.current === mine) setState({ kind: "error", message: err instanceof Error ? err.message : String(err) });
    }
  }, []);
  return { state, run, busy: state.kind === "running" };
}

// ---------- metabolite suggestions shared by every compound input ----------
export interface MetOption { id: string; name: string; compartment?: string }
interface MetCtx { options: MetOption[]; suggest: (q: string) => void }
const MetContext = createContext<MetCtx | null>(null);

export function MetaboliteProvider({ children }: { children: ReactNode }) {
  const { setup } = useSetup();
  const [options, setOptions] = useState<MetOption[]>([]);
  const timer = useRef<ReturnType<typeof setTimeout>>(undefined);
  const last = useRef("");
  const model = setup.model.trim() || "e_coli_core";
  const suggest = useCallback((q: string) => {
    clearTimeout(timer.current);
    const key = `${model}|${q}`;
    if (q.length < 3 || key === last.current) return;
    timer.current = setTimeout(async () => {
      last.current = key;
      try {
        const r = await request<{ metabolites: MetOption[] }>("POST", "/tools/search_metabolites", { query: q, model, limit: 25 });
        setOptions(r.metabolites);
      } catch { /* suggestions are optional */ }
    }, 300);
  }, [model]);
  return <MetContext.Provider value={{ options, suggest }}>{children}</MetContext.Provider>;
}

export function useMetaboliteSuggest() {
  const ctx = useContext(MetContext);
  if (!ctx) throw new Error("useMetaboliteSuggest outside MetaboliteProvider");
  return ctx;
}

// ---------- catalog: models and pathways the server offers ----------
export interface ModelInfo { id: string; description: string }
export interface PathwayInfo { name: string; description: string; product: string; steps: unknown[] }
export interface Catalog { models: ModelInfo[]; pathways: PathwayInfo[] }
export const CatalogContext = createContext<Catalog>({ models: [], pathways: [] });
export const useCatalog = () => useContext(CatalogContext);
