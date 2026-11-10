import { useState } from "react";
import { useHistory, useRun, useSetup, conditions, setupProblem } from "../state";
import { MetaboliteInput, NumberInput, RunButton, RunOutput } from "../components/run";
import { PanelHead, Pill, RawJson } from "../components/ui";
import { firstProblem, fmt, numberProblem } from "../format";

interface Step { reaction: string; name?: string; native: boolean; equation: string; genes?: string[]; also_needs?: string[] }
interface Route {
  precursor: string; precursor_name?: string; target: string; target_name?: string;
  heterologous_steps: number; total_steps: number; precursor_max_flux?: number | null; steps: Step[];
}
interface Result { model: string; target: string; donor: string; routes?: Route[]; summary?: string; note?: string }

export function DesignTab() {
  const { setup } = useSetup();
  const { timed } = useHistory();
  const [target, setTarget] = useState("chor_c");
  const [donor, setDonor] = useState("iJO1366");
  const [maxNew, setMaxNew] = useState(8);
  const [maxSteps, setMaxSteps] = useState(10);
  const { state, run, busy } = useRun<Result>();

  const problem = firstProblem(setupProblem(setup),
    numberProblem("Max new enzymes", maxNew, { min: 0, max: 15, integer: true }),
    numberProblem("Max total steps", maxSteps, { min: 1, max: 20, integer: true }));

  const go = () => run((onStatus) => timed<Result>("find_precursors", {
    ...conditions(setup), target: target.trim(), donor: donor.trim(), max_new_enzymes: maxNew, max_total_steps: maxSteps,
  }, onStatus));

  return (
    <>
      <PanelHead title="Nearest native precursors"
        blurb="Know the compound but not the pathway? Search backwards from it to what the host already makes, and list the enzymes to clone in."
        action={<RunButton busy={busy} onClick={go} problem={problem}>Find routes</RunButton>} />
      <div className="grid small-grid">
        <label>Target metabolite <MetaboliteInput value={target} onChange={setTarget} /></label>
        <label>Donor model
          <input list="model-list" value={donor} autoComplete="off" onChange={(e) => setDonor(e.target.value)} />
          <span className="hint">Leave empty to use only the host's own chemistry.</span>
        </label>
        <label>Max new enzymes <NumberInput value={maxNew} onChange={setMaxNew} min={0} max={15} /></label>
        <label>Max total steps <NumberInput value={maxSteps} onChange={setMaxSteps} min={1} max={20} /></label>
      </div>
      <RunOutput state={state} render={(res) => {
        const routes = res.routes ?? [];
        return (
          <>
            <p>{routes.length
              ? `${routes.length} route${routes.length > 1 ? "s" : ""} from ${res.model} metabolism to ${res.target} (donor: ${res.donor}).`
              : res.summary}</p>
            {routes.map((r, i) => (
              <div className="route" key={i}>
                <h3 style={{ marginTop: 0 }}>Route {i + 1}: {r.precursor_name || r.precursor} → {r.target_name || r.target}</h3>
                <p className="small muted">
                  {r.heterologous_steps} new enzyme{r.heterologous_steps === 1 ? "" : "s"}, {r.total_steps} steps in total.
                  {r.precursor_max_flux != null ? ` The host can make up to ${fmt(r.precursor_max_flux)} mmol/gDW/h of the starting metabolite.` : ""}
                </p>
                <ol>
                  {r.steps.map((st, j) => (
                    <li key={j}>
                      <strong>{st.name || st.reaction}</strong>{" "}
                      {st.native ? <Pill kind="good">host has it</Pill> : <Pill kind="new">clone in</Pill>}
                      <div className="mono">{st.equation}</div>
                      {st.genes?.length ? <div className="muted">Genes: {st.genes.join(", ")}</div> : null}
                      {st.also_needs?.length ? <div className="muted">Also needs: {st.also_needs.join(", ")}</div> : null}
                    </li>
                  ))}
                </ol>
              </div>
            ))}
            {res.note ? <p className="muted small">{res.note}</p> : null}
            <RawJson data={res} />
          </>
        );
      }} />
    </>
  );
}
