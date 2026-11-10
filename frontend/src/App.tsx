import { useEffect, useState, type ReactNode } from "react";
import { call, isBundledFast, probeJobs, request, SLOW, type ToolSpec } from "./api";
import { CatalogContext, HistoryProvider, MetaboliteProvider, SetupProvider, useCatalog, useHistory,
  useMetaboliteSuggest, useSetup, type Catalog, type ModelInfo, type PathwayInfo } from "./state";
import { MetaboliteInput, NumberInput } from "./components/run";
import { Pill } from "./components/ui";
import { fmt } from "./format";
import { YieldTab } from "./tabs/Yield";
import { KnockoutsTab } from "./tabs/Knockouts";
import { DesignTab } from "./tabs/Design";
import { FateTab } from "./tabs/Fate";
import { CheckTab } from "./tabs/Check";
import { CompareTab } from "./tabs/Compare";
import { FluxTab } from "./tabs/Flux";
import { ThermoTab } from "./tabs/Thermo";
import { KineticsTab } from "./tabs/Kinetics";
import { ToolsTab } from "./tabs/Tools";

type TabId = "yield" | "compare" | "flux" | "knockouts" | "design" | "fate" | "thermo" | "kinetics" | "check" | "tools";
const TABS: { id: TabId; label: string }[] = [
  { id: "yield", label: "Yield" },
  { id: "compare", label: "Compare substrates" },
  { id: "flux", label: "Fluxes" },
  { id: "knockouts", label: "Knockouts" },
  { id: "design", label: "Pathway design" },
  { id: "fate", label: "Product fate" },
  { id: "thermo", label: "Thermodynamics" },
  { id: "kinetics", label: "Enzyme kinetics" },
  { id: "check", label: "Model check" },
  { id: "tools", label: "All tools" },
];
const TAB_KEY = "cellfactory.tab";

function savedTab(): TabId {
  try {
    const t = localStorage.getItem(TAB_KEY);
    if (TABS.some((x) => x.id === t)) return t as TabId;
  } catch { /* storage blocked */ }
  return "yield";
}

export function App() {
  const [catalog, setCatalog] = useState<Catalog>({ models: [], pathways: [] });
  const [specs, setSpecs] = useState<ToolSpec[]>([]);
  const [specError, setSpecError] = useState<string | null>(null);
  const [api, setApi] = useState<{ ok: boolean | null; text: string; title?: string }>({ ok: null, text: "Connecting…" });

  useEffect(() => {
    (async () => {
      try {
        const health = await request<{ version?: string }>("GET", "/health");
        setApi({ ok: true, text: `API ${health.version ? "v" + health.version : "connected"}` });
      } catch (err) {
        setApi({ ok: false, text: "API unreachable", title: err instanceof Error ? err.message : String(err) });
      }
      await probeJobs();
      const [models, pathways, tools] = await Promise.allSettled([
        call<{ bundled: ModelInfo[] }>("list_models"),
        call<{ pathways: PathwayInfo[] }>("list_pathways"),
        request<ToolSpec[] | { tools: ToolSpec[] }>("GET", "/tools"),
      ]);
      setCatalog({
        models: models.status === "fulfilled" ? models.value.bundled : [],
        pathways: pathways.status === "fulfilled" ? pathways.value.pathways : [],
      });
      if (tools.status === "fulfilled") {
        const list = Array.isArray(tools.value) ? tools.value : tools.value.tools ?? [];
        for (const t of list) if (t.slow) SLOW.add(t.name); // the API marks its slow tools
        setSpecs(list);
      } else {
        setSpecError("Couldn't load the tool list from the API.");
      }
    })();
  }, []);

  return (
    <CatalogContext.Provider value={catalog}>
      <SetupProvider>
        <MetaboliteProvider>
          <HistoryProvider>
            <header className="top">
              <div className="brand">
                <span className="logo" aria-hidden="true" />
                <div>
                  <h1>cellfactory</h1>
                  <p className="tagline">Feed substrate X to organism Y, ask how much product Z it can make.</p>
                </div>
              </div>
              <div className={`status ${api.ok === true ? "ok" : api.ok === false ? "bad" : ""}`} role="status" title={api.title}>{api.text}</div>
            </header>
            <main>
              <SetupCard />
              <Tabs specs={specs} specError={specError} />
              <History />
            </main>
            <footer className="muted small">Model predictions from flux balance analysis, not lab results: they ignore regulation, kinetics and toxicity.</footer>
            <Datalists />
          </HistoryProvider>
        </MetaboliteProvider>
      </SetupProvider>
    </CatalogContext.Provider>
  );
}

function SetupCard() {
  const { setup, update } = useSetup();
  const { models, pathways } = useCatalog();
  const [info, setInfo] = useState("Any BiGG id or a model file path also works.");
  const chosen = pathways.find((p) => p.name === setup.pathway);

  // Describe the model once the user settles on one (debounced, and only for models the server has).
  useEffect(() => {
    const id = setup.model.trim();
    if (!id) return;
    if (!isBundledFast(id) && !models.some((m) => m.id === id)) {
      setInfo("Not bundled: downloaded from BiGG on first use, if the server can reach it.");
      return;
    }
    let live = true;
    const t = setTimeout(async () => {
      setInfo("Loading model…");
      try {
        const d = await call<{ name?: string; model: string; reactions: number; metabolites: number; genes: number; growth_rate_default_medium: number }>("describe_model", { model: id });
        if (live) setInfo(`${d.name || d.model}: ${d.reactions} reactions, ${d.metabolites} metabolites, ${d.genes} genes. Grows at ${fmt(d.growth_rate_default_medium)} 1/h on its default medium.`);
      } catch (err) {
        if (live) setInfo(err instanceof Error ? err.message : String(err));
      }
    }, 400);
    return () => { live = false; clearTimeout(t); };
  }, [setup.model, models]);

  return (
    <section className="setup card" aria-labelledby="setup-h">
      <h2 id="setup-h">Setup</h2>
      <div className="grid">
        <label>Organism model
          <input list="model-list" value={setup.model} autoComplete="off" onChange={(e) => update({ model: e.target.value })} />
          <span className="hint">{info}</span>
        </label>
        <label>Substrate <MetaboliteInput value={setup.substrate} onChange={(v) => update({ substrate: v })} /></label>
        <label>Product <MetaboliteInput value={setup.product} onChange={(v) => update({ product: v })} /></label>
        <label><span>Uptake <span className="unit">mmol/gDW/h</span></span>
          <NumberInput value={setup.uptake} onChange={(v) => update({ uptake: v })} min={0} step={0.5} />
        </label>
        <label>Non-native pathway
          <select value={setup.pathway} onChange={(e) => update({ pathway: e.target.value })}>
            <option value="">None (native metabolism)</option>
            {pathways.map((p) => <option key={p.name} value={p.name}>{p.name} ({p.steps.length} enzymes)</option>)}
          </select>
          {chosen ? <span className="hint">{chosen.description} Needs a genome-scale model such as iJO1366; set Product to "{chosen.product}".</span> : null}
        </label>
        <fieldset className="oxygen">
          <legend>Oxygen</legend>
          <label className="radio"><input type="radio" name="aerobic" checked={setup.aerobic} onChange={() => update({ aerobic: true })} /> Aerobic</label>
          <label className="radio"><input type="radio" name="aerobic" checked={!setup.aerobic} onChange={() => update({ aerobic: false })} /> Anaerobic</label>
        </fieldset>
      </div>
    </section>
  );
}

function Tabs({ specs, specError }: { specs: ToolSpec[]; specError: string | null }) {
  const [tab, setTab] = useState<TabId>(savedTab);
  const select = (id: TabId) => {
    setTab(id);
    try { localStorage.setItem(TAB_KEY, id); } catch { /* storage blocked */ }
  };
  // Panels stay mounted (just hidden) so a result survives switching tabs.
  const panels: Record<TabId, ReactNode> = {
    yield: <YieldTab />, compare: <CompareTab />, flux: <FluxTab />, knockouts: <KnockoutsTab />,
    design: <DesignTab />, fate: <FateTab />, thermo: <ThermoTab />, kinetics: <KineticsTab />, check: <CheckTab />,
    tools: <ToolsTab specs={specs} error={specError} />,
  };
  return (
    <>
      <nav className="tabs" role="tablist" aria-label="Analyses">
        {TABS.map((t) => (
          <button key={t.id} role="tab" id={`tab-${t.id}`} aria-controls={`panel-${t.id}`} aria-selected={tab === t.id}
            onClick={() => select(t.id)}>{t.label}</button>
        ))}
      </nav>
      {TABS.map((t) => (
        <section key={t.id} id={`panel-${t.id}`} className="panel card" role="tabpanel" aria-labelledby={`tab-${t.id}`} hidden={tab !== t.id}>
          {panels[t.id]}
        </section>
      ))}
    </>
  );
}

function History() {
  const { history } = useHistory();
  return (
    <section className="card jobs" aria-labelledby="jobs-h">
      <h2 id="jobs-h">Recent runs</h2>
      <ol className="history">
        {history.length ? history.map((e, i) => (
          <li key={i}>
            <span className="when">{e.at.toLocaleTimeString()}</span>
            <code>{e.tool}</code>
            {e.ok ? <Pill kind="good">{(e.ms / 1000).toFixed(1)} s</Pill> : <Pill kind="bad">error</Pill>}
            <span className="muted">{Object.entries(e.args).filter(([, v]) => v !== "" && v != null).map(([k, v]) => `${k}=${v}`).join(" ")}</span>
          </li>
        )) : <li className="muted">Nothing run yet.</li>}
      </ol>
    </section>
  );
}

function Datalists() {
  const { models } = useCatalog();
  const { options } = useMetaboliteSuggest();
  return (
    <>
      <datalist id="model-list">{models.map((m) => <option key={m.id} value={m.id}>{m.description}</option>)}</datalist>
      <datalist id="met-list">{options.map((m) => <option key={m.id} value={m.id}>{m.name}{m.compartment ? ` [${m.compartment}]` : ""}</option>)}</datalist>
    </>
  );
}
