import { useState } from "react";
import { useHistory, useRun, useSetup, conditions, setupProblem } from "../state";
import { NumberInput, RunButton, RunOutput } from "../components/run";
import { Caveats, ErrorBox, Notice, PanelHead, Pill, RawJson, Stat, Stats, Table, countLabel, type Column, type RowBlock } from "../components/ui";
import { firstProblem, fmt, numberProblem } from "../format";

interface ThermoRow {
  reaction: string; name: string; flux: number; dg0_prime: number | null; dg0_sigma: number | null;
  dg_prime_min: number | null; dg_prime_max: number | null; thermo_direction: string | null; model_direction: string; method: string;
}
interface Result {
  model: string; substrate: string; route_for: string; conditions: string;
  summary: { reactions_checked: number; estimated: number; active_reactions: number; against_thermodynamics: number;
    against_thermodynamics_group_contribution: number; conflicts: number; can_tighten: number; route_feasible: boolean };
  against_thermodynamics: RowBlock<ThermoRow>; conflicts: RowBlock<ThermoRow>; can_tighten: RowBlock<ThermoRow>;
  unestimated: RowBlock; caveats?: string[]; note?: string;
}

const range = (r: ThermoRow) => r.dg_prime_min == null ? "–" : `${fmt(r.dg_prime_min, 1)} to ${fmt(r.dg_prime_max, 1)}`;
const COLS: Column<ThermoRow>[] = [
  { key: "reaction", label: "Reaction", render: (r) => <><strong>{r.reaction}</strong><div className="muted small">{r.name}</div></> },
  { key: "flux", label: "Flux" },
  { key: "dg0_prime", label: "ΔG'° (kJ/mol)", render: (r) => r.dg0_prime == null ? "–" : `${fmt(r.dg0_prime, 1)} ± ${fmt(r.dg0_sigma, 1)}` },
  { key: "range", label: "ΔG' range", render: range },
  { key: "thermo_direction", label: "Physics allows" },
  { key: "model_direction", label: "Model allows" },
  { key: "method", label: "Estimate" },
];

interface MdfStep {
  reaction: string; name: string; flux: number; dg0_prime: number | null; dg0_sigma: number | null;
  method: string | null; driving_force: number | null; bottleneck: boolean; note: string;
}
interface MdfResult {
  model: string; substrate: string; route_for: string; conditions: string;
  mdf: number | null; mdf_optimistic: number | null; verdict: string; reason: string; bottlenecks: string[];
  summary: { steps: number; estimated: number; excluded: number; bottlenecks: number; route_feasible: boolean };
  steps: RowBlock<MdfStep>; concentrations: RowBlock<{ metabolite: string; concentration: number }>;
  caveats?: string[];
}

export function ThermoTab() {
  const { setup } = useSetup();
  const { timed } = useHistory();
  const [scope, setScope] = useState<"reactions" | "route">("reactions");
  const [goal, setGoal] = useState<"growth" | "product">("product");
  const [minGrowth, setMinGrowth] = useState(0);
  const [all, setAll] = useState(false);
  const { state, run, busy } = useRun<Result | MdfResult>();

  const problem = firstProblem(setupProblem(setup),
    goal === "product" ? numberProblem("Min growth", minGrowth, { min: 0 }) : null);

  const go = () => run((onStatus) => {
    const base = {
      ...conditions(setup), pathway: setup.pathway,
      ...(goal === "product" ? { product: setup.product.trim(), min_growth: minGrowth } : {}),
    };
    return scope === "route"
      ? timed<MdfResult>("check_driving_force", base, onStatus)
      : timed<Result>("check_thermodynamics", { ...base, all_reactions: all }, onStatus);
  });

  return (
    <>
      <PanelHead title="Thermodynamic feasibility"
        blurb="Does the route ask any reaction to run in a direction physics rules out? ΔG' estimates from eQuilibrator, across 1 µM to 10 mM metabolite levels."
        action={<RunButton busy={busy} onClick={go} problem={problem}>Check route</RunButton>} />
      <div className="grid small-grid">
        <label>Question
          <select value={scope} onChange={(e) => setScope(e.target.value as "reactions" | "route")}>
            <option value="reactions">Each reaction on its own (ΔG' directions)</option>
            <option value="route">Every step at once (max-min driving force)</option>
          </select>
          <span className="hint">{scope === "route"
            ? "One concentration per metabolite for the whole route, so neighbouring steps can't each assume the shared intermediate sits at a different level."
            : "Whether some concentrations in range let each reaction run the way the route uses it."}</span>
        </label>
        <label>Route to check
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
        {scope === "reactions" ? (
          <label className="check"><input type="checkbox" checked={all} onChange={(e) => setAll(e.target.checked)} />
            Check every reaction <span className="hint">about a minute on genome scale</span></label>
        ) : null}
      </div>
      <RunOutput state={state} render={(data) => {
        if ("mdf" in data) return <Mdf res={data} />;
        const res = data;
        const s = res.summary;
        return (
          <>
            {s.route_feasible
              ? <Notice>No active reaction runs in a direction the estimates rule out. That rules nothing in: a feasible direction can still be slow.</Notice>
              : <ErrorBox><strong>{s.against_thermodynamics} reaction{s.against_thermodynamics === 1 ? "" : "s"} run against thermodynamics. </strong>
                  {s.against_thermodynamics_group_contribution
                    ? `${s.against_thermodynamics_group_contribution} rest on the weaker group-contribution estimate, so may be the estimate's fault.`
                    : "All rest on the stronger estimate."}</ErrorBox>}
            <Stats>
              <Stat k="Reactions checked" v={String(s.reactions_checked)} s={`${s.estimated} with an estimate`} />
              <Stat k="Against thermodynamics" v={String(s.against_thermodynamics)} s="active, impossible direction" />
              <Stat k="Model conflicts" v={String(s.conflicts)} s="only allowed the impossible way" />
              <Stat k="Could be one-way" v={String(s.can_tighten)} s="reversible in the model" />
            </Stats>
            <p className="muted small">{res.model}, {res.substrate}, route for {res.route_for}. {res.conditions}</p>
            {s.against_thermodynamics ? <><h3>Against thermodynamics {countLabel(res.against_thermodynamics)}</h3><Table rows={res.against_thermodynamics.rows} columns={COLS} /></> : null}
            <h3>Conflicts with the model's bounds {countLabel(res.conflicts)}</h3>
            <Table rows={res.conflicts.rows} columns={COLS} />
            <details className="raw"><summary>Reversible reactions that could be made one-way {countLabel(res.can_tighten)}</summary>
              <Table rows={res.can_tighten.rows} columns={COLS} /></details>
            {res.unestimated.total ? <details className="raw"><summary>No estimate {countLabel(res.unestimated, 20)}</summary><Table rows={res.unestimated.rows} /></details> : null}
            <Caveats list={res.caveats} />
            {res.note ? <p className="muted small">{res.note}</p> : null}
            <RawJson data={res} />
          </>
        );
      }} />
    </>
  );
}

const MDF_COLS: Column<MdfStep>[] = [
  { key: "reaction", label: "Step", render: (r) => <>
      <strong>{r.reaction}</strong>{r.bottleneck ? <> <Pill kind="bad">bottleneck</Pill></> : null}
      <div className="muted small">{r.name}</div></> },
  { key: "flux", label: "Flux" },
  { key: "driving_force", label: "Driving force (kJ/mol)", render: (r) => r.driving_force == null ? "–" : fmt(r.driving_force, 1) },
  { key: "dg0_prime", label: "ΔG'° (kJ/mol)", render: (r) => r.dg0_prime == null ? "–" : `${fmt(r.dg0_prime, 1)} ± ${fmt(r.dg0_sigma, 1)}` },
  { key: "method", label: "Estimate", render: (r) => r.method ?? "–" },
  { key: "note", label: "Left out because", render: (r) => r.note || "–" },
];

/** Max-min driving force: one verdict for the whole route, plus the steps holding it down. */
function Mdf({ res }: { res: MdfResult }) {
  const s = res.summary;
  const bad = res.verdict === "infeasible" || res.verdict === "uncertain";
  return (
    <>
      {bad
        ? <ErrorBox><strong>The route is {res.verdict}: </strong>{res.reason}.</ErrorBox>
        : <Notice><strong>The route is {res.verdict}: </strong>{res.reason}. A route that can run is not a route that runs fast.</Notice>}
      <Stats>
        <Stat k="Max-min driving force" v={fmt(res.mdf, 1)} s="kJ/mol, best estimates" />
        <Stat k="Allowing for estimate error" v={fmt(res.mdf_optimistic, 1)} s="kJ/mol, optimistic" />
        <Stat k="Steps" v={String(s.steps)} s={`${s.estimated} with a ΔG'° estimate`} />
        <Stat k="Bottleneck steps" v={String(s.bottlenecks)} s={res.bottlenecks.slice(0, 3).join(", ") || "none"} />
      </Stats>
      <p className="muted small">{res.model}, {res.substrate}, route for {res.route_for}. {res.conditions} Below about 5.7 kJ/mol a step runs backwards more than once every ten turnovers.</p>
      <h3>Steps, least downhill first {countLabel(res.steps, 200)}</h3>
      <Table rows={res.steps.rows} columns={MDF_COLS} />
      <details className="raw"><summary>Metabolite concentrations at the optimum {countLabel(res.concentrations, 200)}</summary>
        <Table rows={res.concentrations.rows} columns={[
          { key: "metabolite", label: "Metabolite" },
          { key: "concentration", label: "Concentration (M)", render: (r) => r.concentration.toExponential(2) },
        ]} /></details>
      <Caveats list={res.caveats} />
      <RawJson data={res} />
    </>
  );
}
