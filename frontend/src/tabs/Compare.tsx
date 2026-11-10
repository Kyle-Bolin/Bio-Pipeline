import { useState } from "react";
import { useHistory, useRun, useSetup, conditions, setupProblem } from "../state";
import { BarChart } from "../components/charts";
import { RunButton, RunOutput } from "../components/run";
import { Notice, PanelHead, RawJson, Table } from "../components/ui";
import { firstProblem, fmt, pct, splitList } from "../format";

interface Row {
  substrate: string; aerobic: boolean; error: string | null; growth?: number;
  guaranteed_yield?: number; max_theoretical_yield?: number; guaranteed_carbon_yield?: number; max_carbon_yield?: number;
}
interface Result { model: string; product: string; rows: Row[]; best: Row | null }

const condition = (r: Row) => `${r.substrate}, ${r.aerobic ? "aerobic" : "anaerobic"}`;

export function CompareTab() {
  const { setup } = useSetup();
  const { timed } = useHistory();
  const [substrates, setSubstrates] = useState("glucose, fructose, pyruvate");
  const [oxygen, setOxygen] = useState("both");
  const { state, run, busy } = useRun<Result>();
  const count = splitList(substrates).length;

  const problem = firstProblem(setupProblem(setup),
    count < 1 || count > 8 ? "List between 1 and 8 substrates." : null);

  const go = () => run((onStatus) => {
    const { model, uptake } = conditions(setup);
    return timed<Result>("compare_conditions",
      { model, uptake, product: setup.product.trim(), substrates, oxygen, pathway: setup.pathway }, onStatus);
  });

  return (
    <>
      <PanelHead title="Compare substrates"
        blurb="Which feedstock, with or without oxygen, turns the most of its carbon into the product? Carbon yield is the fair comparison between substrates of different size."
        action={<RunButton busy={busy} onClick={go} problem={problem}>Compare</RunButton>} />
      <div className="grid small-grid">
        <label className="wide">Substrates
          <input value={substrates} onChange={(e) => setSubstrates(e.target.value)} />
          <span className="hint">Comma-separated, up to 8. {count > 8 ? "Too many." : `${count} listed.`} Product and uptake come from Setup.</span>
        </label>
        <label>Oxygen
          <select value={oxygen} onChange={(e) => setOxygen(e.target.value)}>
            <option value="both">Both</option>
            <option value="aerobic">Aerobic only</option>
            <option value="anaerobic">Anaerobic only</option>
          </select>
        </label>
      </div>
      <RunOutput state={state} render={(res) => (
        <>
          {res.best
            ? <p>Best: <strong>{condition(res.best)}</strong>, with {pct(res.best.guaranteed_carbon_yield)} of the substrate's carbon guaranteed to end up in {res.product}.</p>
            : <Notice>The cell couldn't grow under any of these conditions.</Notice>}
          <BarChart ariaLabel="Guaranteed and maximum carbon yield for each condition"
            legend={["guaranteed at max growth", "theoretical max"]} format={(x) => pct(x)}
            rows={res.rows.map((r) => ({
              label: condition(r), value: r.error ? null : r.guaranteed_carbon_yield ?? 0,
              ceiling: r.error ? null : r.max_carbon_yield, note: r.error ?? undefined,
            }))} />
          <Table rows={res.rows} columns={[
            { key: "substrate", label: "Substrate" },
            { key: "aerobic", label: "Oxygen", render: (r) => r.aerobic ? "aerobic" : "anaerobic" },
            { key: "growth", label: "Growth (1/h)", render: (r) => r.error ? "" : fmt(r.growth) },
            { key: "guaranteed_yield", label: "Guaranteed (mol/mol)", render: (r) => r.error ? "" : fmt(r.guaranteed_yield) },
            { key: "max_theoretical_yield", label: "Max (mol/mol)", render: (r) => r.error ? "" : fmt(r.max_theoretical_yield) },
            { key: "guaranteed_carbon_yield", label: "Carbon yield", render: (r) => r.error ? "" : pct(r.guaranteed_carbon_yield) },
            { key: "error", label: "Problem", render: (r) => r.error ?? "" },
          ]} />
          <RawJson data={res} />
        </>
      )} />
    </>
  );
}
