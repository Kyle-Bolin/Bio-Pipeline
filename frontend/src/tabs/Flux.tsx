import { useState } from "react";
import { useHistory, useRun, useSetup, conditions, setupProblem } from "../state";
import { FluxBar } from "../components/charts";
import { NumberInput, RunButton, RunOutput } from "../components/run";
import { PanelHead, RawJson, Stat, Stats, Table } from "../components/ui";
import { firstProblem, fmt, numberProblem } from "../format";

interface Exchange { reaction: string; metabolite: string; name: string; flux: number }
interface Reaction { reaction: string; name: string; flux: number; equation: string; genes: string[] }
interface Result {
  model: string; substrate: string; optimised_for: string; growth: number; product_flux: number | null;
  active_reactions: number; uptake: Exchange[]; secretion: Exchange[]; reactions: Reaction[]; truncated: boolean;
}

export function FluxTab() {
  const { setup } = useSetup();
  const { timed } = useHistory();
  const [goal, setGoal] = useState<"growth" | "product">("growth");
  const [minGrowth, setMinGrowth] = useState(0.05);
  const [limit, setLimit] = useState(30);
  const { state, run, busy } = useRun<Result>();

  const problem = firstProblem(setupProblem(setup),
    goal === "product" ? numberProblem("Min growth", minGrowth, { min: 0 }) : null,
    numberProblem("Reactions to list", limit, { min: 1, max: 200, integer: true }));

  const go = () => run((onStatus) => timed<Result>("flux_distribution", {
    ...conditions(setup), pathway: setup.pathway, limit,
    ...(goal === "product" ? { product: setup.product.trim(), min_growth: minGrowth } : {}),
  }, onStatus));

  return (
    <>
      <PanelHead title="Where the carbon goes"
        blurb="One flux distribution (parsimonious FBA): what the cell takes up, what it secretes, and its busiest reactions."
        action={<RunButton busy={busy} onClick={go} problem={problem}>Show fluxes</RunButton>} />
      <div className="grid small-grid">
        <label>Optimise for
          <select value={goal} onChange={(e) => setGoal(e.target.value as "growth" | "product")}>
            <option value="growth">Growth (wild type)</option>
            <option value="product">Product, from Setup</option>
          </select>
        </label>
        {goal === "product" ? (
          <label><span>Min growth <span className="unit">1/h</span></span>
            <NumberInput value={minGrowth} onChange={setMinGrowth} min={0} step={0.01} />
            <span className="hint">0 is a non-growing cell; about 0.05 is a growing production strain.</span></label>
        ) : null}
        <label>Reactions to list <NumberInput value={limit} onChange={setLimit} min={1} max={200} /></label>
      </div>
      <RunOutput state={state} render={(res) => {
        const exMax = Math.max(1e-9, ...[...res.uptake, ...res.secretion].map((e) => Math.abs(e.flux)));
        const rxMax = Math.max(1e-9, ...res.reactions.map((r) => Math.abs(r.flux)));
        const exCols = [
          { key: "name", label: "Compound", render: (r: Exchange) => <>{r.name || r.metabolite} <span className="muted mono">{r.metabolite}</span></> },
          { key: "flux", label: "Flux", render: (r: Exchange) => <span className="num-cell">{fmt(r.flux)}</span> },
          { key: "bar", label: "", render: (r: Exchange) => <FluxBar value={r.flux} max={exMax} /> },
        ];
        return (
          <>
            <Stats>
              <Stat k="Growth rate" v={fmt(res.growth)} s="1/h" />
              {res.product_flux != null ? <Stat k="Product flux" v={fmt(res.product_flux)} s={`${res.optimised_for}, mmol/gDW/h`} /> : null}
              <Stat k="Active reactions" v={String(res.active_reactions)} s="internal, carrying flux" />
              <Stat k="Uptakes / secretions" v={`${res.uptake.length} / ${res.secretion.length}`} />
            </Stats>
            <div className="two-col">
              <div><h3>Taken up</h3><Table rows={res.uptake} columns={exCols} /></div>
              <div><h3>Secreted</h3><Table rows={res.secretion} columns={exCols} /></div>
            </div>
            <h3>Busiest reactions{res.truncated ? ` (top ${res.reactions.length})` : ""}</h3>
            <p className="muted small">Bars to the left run the reaction right to left. Units mmol/gDW/h.</p>
            <Table rows={res.reactions} columns={[
              { key: "reaction", label: "Reaction", render: (r) => <><strong>{r.reaction}</strong><div className="muted small">{r.name}</div></> },
              { key: "flux", label: "Flux", render: (r) => <span className="num-cell">{fmt(r.flux)}</span> },
              { key: "bar", label: "", render: (r) => <FluxBar value={r.flux} max={rxMax} /> },
              { key: "equation", label: "Equation", render: (r) => <span className="mono">{r.equation}</span> },
              { key: "genes", label: "Genes" },
            ]} />
            <p className="muted small">The optimum is often not unique, so other equally good routes are possible.</p>
            <RawJson data={res} />
          </>
        );
      }} />
    </>
  );
}
