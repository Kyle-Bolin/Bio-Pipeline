import { useState } from "react";
import { useHistory, useRun, useSetup, conditions } from "../state";
import { NumberInput, RunButton, RunOutput } from "../components/run";
import { Caveats, ErrorBox, Notice, PanelHead, RawJson, Stat, Stats, Table, countLabel, type Column, type RowBlock } from "../components/ui";
import { fmt } from "../format";

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

export function ThermoTab() {
  const { setup } = useSetup();
  const { timed } = useHistory();
  const [goal, setGoal] = useState<"growth" | "product">("product");
  const [minGrowth, setMinGrowth] = useState(0);
  const [all, setAll] = useState(false);
  const { state, run, busy } = useRun<Result>();

  const go = () => run((onStatus) => timed<Result>("check_thermodynamics", {
    ...conditions(setup), pathway: setup.pathway, all_reactions: all,
    ...(goal === "product" ? { product: setup.product.trim(), min_growth: minGrowth } : {}),
  }, onStatus));

  return (
    <>
      <PanelHead title="Thermodynamic feasibility"
        blurb="Does the route ask any reaction to run in a direction physics rules out? ΔG' estimates from eQuilibrator, across 1 µM to 10 mM metabolite levels."
        action={<RunButton busy={busy} onClick={go}>Check route</RunButton>} />
      <div className="grid small-grid">
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
        <label className="check"><input type="checkbox" checked={all} onChange={(e) => setAll(e.target.checked)} />
          Check every reaction <span className="hint">about a minute on genome scale</span></label>
      </div>
      <RunOutput state={state} render={(res) => {
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
