import { useState } from "react";
import { useHistory, useRun, useSetup, conditions, setupProblem } from "../state";
import { NumberInput, RunButton, RunOutput } from "../components/run";
import { Caveats, Notice, PanelHead, Pill, RawJson, Stat, Stats, Table, countLabel, type Column, type RowBlock } from "../components/ui";
import { FluxBar } from "../components/charts";
import { firstProblem, fmt, numberProblem } from "../format";

interface KineticsRow {
  reaction: string; flux: number | null; ec: string; ec_source: string; value_type: string;
  kcat_per_s: number | null; kcat_min: number | null; kcat_max: number | null; n_kcat: number;
  kcat_organism: string | null; organism_match: string | null; substrate_matched: boolean | null;
  source: string | null; km_mM: number | null; km_substrate: string | null; n_km: number; note: string | null;
}
interface Result {
  model: string; substrate: string; route_for: string; organism: string | null; sources: string[];
  summary: { reactions: number; enzyme_steps: number; measured: number; predicted: number; missing: number;
    measured_in_host: number };
  reactions: RowBlock<KineticsRow>; caveats?: string[];
}

/** A measured kcat from the host is the one to trust; a predicted one can be ten times off. */
const valuePill = (t: string) => (t === "measured" ? "good" : t === "predicted" ? "new" : undefined);

export function KineticsTab() {
  const { setup } = useSetup();
  const { timed } = useHistory();
  const [goal, setGoal] = useState<"growth" | "product">("product");
  const [minGrowth, setMinGrowth] = useState(0);
  const [organism, setOrganism] = useState("");
  const [online, setOnline] = useState(true);
  const { state, run, busy } = useRun<Result>();

  const problem = firstProblem(setupProblem(setup),
    goal === "product" ? numberProblem("Min growth", minGrowth, { min: 0 }) : null);

  const go = () => run((onStatus) => timed<Result>("check_enzyme_kinetics", {
    ...conditions(setup), pathway: setup.pathway, organism: organism.trim(), online,
    ...(goal === "product" ? { product: setup.product.trim(), min_growth: minGrowth } : {}),
  }, onStatus));

  const maxFlux = state.kind === "done"
    ? Math.max(1e-9, ...state.data.reactions.rows.map((r) => Math.abs(r.flux ?? 0))) : 1;

  const COLS: Column<KineticsRow>[] = [
    { key: "reaction", label: "Reaction", render: (r) => <><strong>{r.reaction}</strong>
      <div className="muted small">{r.ec ? `EC ${r.ec}` : "no EC number"}{r.ec_source ? ` (${r.ec_source})` : ""}</div></> },
    { key: "flux", label: "Flux", render: (r) => <><span className="num">{fmt(r.flux)}</span>
      <FluxBar value={r.flux ?? 0} max={maxFlux} /></> },
    { key: "value_type", label: "kcat", render: (r) => <Pill kind={valuePill(r.value_type)}>{r.value_type}</Pill> },
    { key: "kcat_per_s", label: "kcat (1/s)", render: (r) => r.kcat_per_s == null ? "–"
      : <>{fmt(r.kcat_per_s)}{r.n_kcat > 1 ? <div className="muted small">{fmt(r.kcat_min)}–{fmt(r.kcat_max)}, n={r.n_kcat}</div> : null}</> },
    { key: "organism_match", label: "Measured in", render: (r) => r.kcat_organism
      ? <>{r.kcat_organism}<div className="muted small">{r.organism_match}{r.substrate_matched === false ? ", other substrate" : ""}</div></> : "–" },
    { key: "km_mM", label: "Km (mM)", render: (r) => r.km_mM == null ? "–"
      : <>{fmt(r.km_mM)}<div className="muted small">{r.km_substrate}</div></> },
    { key: "source", label: "Source", render: (r) => r.source ?? "–" },
    { key: "note", label: "Note", render: (r) => r.note || "–" },
  ];

  return (
    <>
      <PanelHead title="Enzyme kinetics"
        blurb="How fast each enzyme on the route turns over: measured kcat and Km from BRENDA and SABIO-RK, matched by EC number. Flux divided by kcat is the least enzyme a step needs, so a slow enzyme on a busy step is where a real strain stalls."
        action={<RunButton busy={busy} onClick={go} problem={problem}>Look up kinetics</RunButton>} />
      <div className="grid small-grid">
        <label>Route to look up
          <select value={goal} onChange={(e) => setGoal(e.target.value as "growth" | "product")}>
            <option value="product">Making the product (Setup)</option>
            <option value="growth">Growing</option>
          </select>
        </label>
        {goal === "product" ? (
          <label><span>Min growth <span className="unit">1/h</span></span>
            <NumberInput value={minGrowth} onChange={setMinGrowth} min={0} step={0.01} />
            <span className="hint">0 is the stoichiometric best case; about 0.05 for a growing strain.</span></label>
        ) : null}
        <label>Host species
          <input value={organism} placeholder="from the model id" autoComplete="off" onChange={(e) => setOrganism(e.target.value)} />
          <span className="hint">E.g. "Escherichia coli". Values measured in the host are the ones to trust.</span>
        </label>
        <label className="check"><input type="checkbox" checked={online} onChange={(e) => setOnline(e.target.checked)} />
          Query SABIO-RK <span className="hint">one request per EC number the first time, then cached on the server</span></label>
      </div>
      <RunOutput state={state} render={(res) => {
        const s = res.summary;
        return (
          <>
            {s.enzyme_steps && !s.measured && !s.predicted
              ? <Notice>No kcat was found for any step. BRENDA's file has to be downloaded by the server's owner (its licence forbids bundling it), and SABIO-RK has to be reachable.</Notice>
              : null}
            <Stats>
              <Stat k="Enzyme steps" v={String(s.enzyme_steps)} s={`of ${s.reactions} reactions on the route`} />
              <Stat k="Measured kcat" v={String(s.measured)} s={`${s.measured_in_host} in the host species`} />
              <Stat k="Predicted kcat" v={String(s.predicted)} s="can be ten times off" />
              <Stat k="No kcat at all" v={String(s.missing)} />
            </Stats>
            <p className="muted small">{res.model}, {res.substrate}, route for {res.route_for}.
              Host {res.organism || "unknown"}. Sources: {res.sources.length ? res.sources.join(", ") : "none available"}.</p>
            <h3>Route steps, busiest first {countLabel(res.reactions)}</h3>
            <Table rows={res.reactions.rows} columns={COLS} />
            <Caveats list={res.caveats} />
            <RawJson data={res} />
          </>
        );
      }} />
    </>
  );
}
