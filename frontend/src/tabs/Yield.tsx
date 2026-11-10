import { useHistory, useRun, useSetup, conditions } from "../state";
import { EnvelopeChart, type EnvelopePoint } from "../components/charts";
import { RunButton, RunOutput } from "../components/run";
import { Caveats, PanelHead, RawJson, Stat, Stats, Table } from "../components/ui";
import { fmt, pct } from "../format";

interface Sim {
  model_id: string; substrate: string; product: string; aerobic: boolean; uptake: number;
  production: { growth: number; guaranteed_yield: number; best_case_yield: number; max_theoretical_yield: number;
    carbon_yield: number; target_flux_min: number };
  pathway?: { name: string; genes: string[] } | null;
  caveats?: string[];
  checks?: Record<string, Record<string, unknown>[]>;
}
interface Env { points: EnvelopePoint[] }

export function YieldTab() {
  const { setup } = useSetup();
  const { timed } = useHistory();
  const { state, run, busy } = useRun<{ sim: Sim; env: Env | null }>();

  const go = () => run(async (onStatus) => {
    const args = { ...conditions(setup), product: setup.product.trim(), pathway: setup.pathway };
    const [sim, env] = await Promise.all([
      timed<Sim>("simulate_production", args, onStatus),
      timed<Env>("production_envelope", { ...args, points: 25 }, onStatus).catch(() => null),
    ]);
    return { sim, env };
  });

  return (
    <>
      <PanelHead title="Growth and yield"
        blurb="How fast the cell grows and how much product it makes per mole of substrate, plus the growth vs. production trade-off."
        action={<RunButton busy={busy} onClick={go}>Simulate</RunButton>} />
      <RunOutput state={state} render={({ sim, env }) => {
        const p = sim.production;
        return (
          <>
            <Stats>
              <Stat k="Growth rate" v={fmt(p.growth)} s="1/h" />
              <Stat k="Guaranteed yield" v={fmt(p.guaranteed_yield)} s="mol product / mol substrate, at max growth" />
              <Stat k="Best-case yield" v={fmt(p.best_case_yield)} s="at max growth" />
              <Stat k="Theoretical max" v={fmt(p.max_theoretical_yield)} s="ignoring growth" />
              <Stat k="Carbon yield" v={pct(p.carbon_yield)} s="of substrate carbon, guaranteed" />
            </Stats>
            <p className="muted small">
              {sim.model_id}: {sim.substrate} → {sim.product}, {sim.aerobic ? "aerobic" : "anaerobic"}, uptake {fmt(sim.uptake)} mmol/gDW/h.
            </p>
            {sim.pathway ? (
              <div>
                <h3>Added pathway: {sim.pathway.name}</h3>
                <p className="small">Genes cloned in: {sim.pathway.genes.join(", ")}</p>
              </div>
            ) : null}
            {env?.points?.length ? <EnvelopeChart points={env.points} marker={{ growth: p.growth, flux: p.target_flux_min }} /> : null}
            <Caveats list={sim.caveats} />
            {Object.entries(sim.checks ?? {}).filter(([, rows]) => rows && rows.length).map(([title, rows]) => (
              <details className="raw" key={title}><summary>{title}</summary><Table rows={rows} /></details>
            ))}
            <RawJson data={{ simulate_production: sim, production_envelope: env }} />
          </>
        );
      }} />
    </>
  );
}
