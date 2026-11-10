import { useState } from "react";
import { useHistory, useRun, useSetup, conditions, setupProblem } from "../state";
import { NumberInput, RunButton, RunOutput } from "../components/run";
import { Caveats, ErrorBox, Notice, PanelHead, Pill, RawJson, Stat, Stats, Table, countLabel, type RowBlock } from "../components/ui";
import { firstProblem, fmt, numberProblem } from "../format";

interface FbaResult {
  summary: { wild_type_guaranteed_yield: number; best_guaranteed_yield: number; knockouts: string[]; genes: string[]; improved: boolean };
  knockouts: { knocked_out: string[]; genes: string[]; growth: number; guaranteed_yield: number }[];
  caveats?: string[];
}
interface MomaResult {
  summary: { parent_yield: number; best_yield: number; genes: string[]; improved: boolean };
  steps: { knocked_out: string[]; gene_names: string[]; reactions: string[]; growth: number; product_flux: number; product_yield: number }[];
  note?: string;
}
export interface DesignRow {
  rank: number; knockouts: string; n_knockouts: number; growth: number; guaranteed_yield: number;
  best_case_yield: number; coupling: string; genes: string; also_disabled: string;
  gene_level_growth: number | null; gene_level_guaranteed_yield: number | null;
}
interface DesignResult {
  method: string; solver: string | null; status: string; wild_type_guaranteed_yield: number;
  min_growth: number; min_yield: number | null;
  summary: { designs: number; best_guaranteed_yield: number | null; best_knockouts: string | null;
    best_genes: string | null; best_coupling: string | null; improved: boolean };
  designs: RowBlock<DesignRow>; notes?: string[]; caveats?: string[];
}
type Method = "fba" | "moma" | "design";
type Result = FbaResult | MomaResult | DesignResult;

const FORMULATIONS: { id: string; label: string; hint: string }[] = [
  { id: "optknock", label: "OptKnock (optimistic)", hint: "Most product at the mutant's fastest growth, assuming the cell picks the flux state we like among the equally fast ones." },
  { id: "robustknock", label: "RobustKnock (guaranteed)", hint: "Most guaranteed product at fastest growth: the yield reported everywhere else here." },
  { id: "mcs", label: "Minimal cut sets (strong coupling)", hint: "Fewest deletions after which every growing state makes the product, not just the fastest one." },
];

/** "strong" coupling is the design worth building; "none" means the MILP's design didn't survive the FBA re-check. */
const couplingPill = (c: string) => (c === "strong" ? "good" : c === "none" || c === "too slow" ? "bad" : undefined);

export function KnockoutsTab() {
  const { setup } = useSetup();
  const { timed } = useHistory();
  const [method, setMethod] = useState<Method>("fba");
  const [max, setMax] = useState(3);
  const [minGrowth, setMinGrowth] = useState(0.05);
  const [formulation, setFormulation] = useState("optknock");
  const [minYield, setMinYield] = useState(0);
  const [timeLimit, setTimeLimit] = useState(120);
  const { state, run, busy } = useRun<Result>();

  const maxAllowed = method === "moma" ? 5 : 10;
  const problem = firstProblem(setupProblem(setup),
    numberProblem("Max deletions", max, { min: 1, max: maxAllowed, integer: true }),
    method === "moma" ? null : numberProblem("Min growth", minGrowth, { min: 0 }),
    method === "design" && formulation === "mcs" ? numberProblem("Min yield", minYield, { min: 0 }) : null,
    method === "design" ? numberProblem("Solver time limit", timeLimit, { min: 1, max: 3600 }) : null);

  const go = () => run((onStatus) => {
    const base = { ...conditions(setup), product: setup.product.trim(), pathway: setup.pathway, max_knockouts: max };
    if (method === "moma") return timed<MomaResult>("search_knockouts_moma", base, onStatus);
    if (method === "design") {
      return timed<DesignResult>("design_knockouts", {
        ...base, method: formulation, min_growth: minGrowth, time_limit: timeLimit,
        ...(formulation === "mcs" ? { min_yield: minYield } : {}),
      }, onStatus);
    }
    return timed<FbaResult>("search_knockouts", { ...base, min_growth: minGrowth }, onStatus);
  });

  const chosen = FORMULATIONS.find((f) => f.id === formulation);
  return (
    <>
      <PanelHead title="Gene knockouts"
        blurb="Search for deletions that force the cell to make more product. One at a time (FBA or MOMA), or a whole combination at once with a MILP solver."
        action={<RunButton busy={busy} onClick={go} problem={problem}>{method === "design" ? "Design" : "Search"}</RunButton>} />
      <div className="grid small-grid">
        <label>Method
          <select value={method} onChange={(e) => {
            const m = e.target.value as Method;
            setMethod(m);
            if (m === "moma" && max > 5) setMax(2);
          }}>
            <option value="fba">FBA, one deletion at a time (growth-coupled)</option>
            <option value="moma">MOMA, one deletion at a time (non-native products)</option>
            <option value="design">Designed knockouts, a combination at once (StrainDesign MILP)</option>
          </select>
          <span className="hint">{method === "design"
            ? "No single deletion has to help on its own: the solver searches every combination of up to this many."
            : method === "moma"
              ? "Assumes a fresh mutant keeps its parent's fluxes, which is how a deletion can push carbon to a product the cell has no reason to make."
              : "Keeps the single deletion that raises the guaranteed yield most, then repeats."}</span>
        </label>
        {method === "design" ? (
          <label>Formulation
            <select value={formulation} onChange={(e) => setFormulation(e.target.value)}>
              {FORMULATIONS.map((f) => <option key={f.id} value={f.id}>{f.label}</option>)}
            </select>
            <span className="hint">{chosen?.hint}</span>
          </label>
        ) : null}
        <label>Max deletions <NumberInput value={max} onChange={setMax} min={1} max={maxAllowed} /></label>
        {method === "moma" ? null : (
          <label><span>Min growth <span className="unit">1/h</span></span>
            <NumberInput value={minGrowth} onChange={setMinGrowth} min={0} step={0.01} /></label>
        )}
        {method === "design" && formulation === "mcs" ? (
          <label><span>Min yield <span className="unit">mol/mol</span></span>
            <NumberInput value={minYield} onChange={setMinYield} min={0} step={0.1} />
            <span className="hint">0 asks for 30% of the theoretical maximum.</span></label>
        ) : null}
        {method === "design" ? (
          <label><span>Solver time limit <span className="unit">s</span></span>
            <NumberInput value={timeLimit} onChange={setTimeLimit} min={1} max={3600} step={10} />
            <span className="hint">At most an hour: the solve holds the model while it runs. A design that hits the limit is a good one, not provably the best.</span></label>
        ) : null}
      </div>
      <RunOutput state={state} render={(res) => "designs" in res ? <Design res={res} />
        : "steps" in res ? <Moma res={res} /> : <Fba res={res} />} />
    </>
  );
}

function Design({ res }: { res: DesignResult }) {
  const s = res.summary;
  return (
    <>
      {res.status === "infeasible" && !s.designs
        ? <Notice>No combination of deletions forces the product out {res.min_yield ? `at ${fmt(res.min_yield)} mol/mol or more ` : ""}within these limits. Try more deletions, a lower minimum growth, or a less strict formulation.</Notice>
        : res.status === "error" ? <ErrorBox><strong>The solver failed. </strong>{(res.notes ?? []).join(" ")}</ErrorBox>
        : res.status === "time limit" ? <Notice>The solver hit its time limit, so better designs may exist.</Notice>
        : null}
      <Stats>
        <Stat k="Wild-type guaranteed yield" v={fmt(res.wild_type_guaranteed_yield)} s="mol/mol" />
        <Stat k="Best design yield" v={fmt(s.best_guaranteed_yield)} s="mol/mol, FBA re-checked" />
        <Stat k="Designs found" v={String(s.designs)} s={`${res.method}${res.solver ? `, ${res.solver}` : ""}`} />
        <Stat k="Best coupling" v={s.best_coupling ? <Pill kind={couplingPill(s.best_coupling)}>{s.best_coupling}</Pill> : "–"}
          s="strong = forced at every growth rate" />
      </Stats>
      {s.best_genes ? <p className="muted small">Best design deletes <strong>{s.best_knockouts}</strong> (genes {s.best_genes}).</p> : null}
      <h3>Designs {countLabel(res.designs)}</h3>
      <Table rows={res.designs.rows} columns={[
        { key: "rank", label: "#" },
        { key: "knockouts", label: "Reactions to delete", render: (r) => r.knockouts || "(none)" },
        { key: "genes", label: "Genes to delete", render: (r) => r.genes || "(no genes)" },
        { key: "growth", label: "Growth (1/h)" },
        { key: "guaranteed_yield", label: "Guaranteed yield" },
        { key: "coupling", label: "Coupling", render: (r) => <Pill kind={couplingPill(r.coupling)}>{r.coupling}</Pill> },
        { key: "also_disabled", label: "Genes also disable", render: (r) => r.also_disabled || "–" },
        { key: "gene_level_guaranteed_yield", label: "Yield at gene level", render: (r) => fmt(r.gene_level_guaranteed_yield) },
      ]} />
      <Caveats list={res.caveats} />
      <RawJson data={res} />
    </>
  );
}

function Moma({ res }: { res: MomaResult }) {
  const s = res.summary;
  return (
    <>
      <Stats>
        <Stat k="Parent yield" v={fmt(s.parent_yield)} s="mol/mol" />
        <Stat k="Best yield" v={fmt(s.best_yield)} s="mol/mol" />
        <Stat k="Genes deleted" v={String(s.genes.length)} />
      </Stats>
      {s.improved ? null : <Notice>No deletion improved the predicted yield.</Notice>}
      <Table rows={res.steps.map((r, i) => ({ ...r, step: i }))} columns={[
        { key: "step", label: "Step" },
        { key: "knocked_out", label: "Genes deleted", render: (r) => r.knocked_out.length
          ? r.knocked_out.map((g, i) => `${r.gene_names[i] || g} (${g})`).join(", ") : "(parent)" },
        { key: "reactions", label: "Reactions disabled" },
        { key: "growth", label: "Growth (1/h)" },
        { key: "product_flux", label: "Product flux" },
        { key: "product_yield", label: "Yield (mol/mol)" },
      ]} />
      {res.note ? <p className="muted small">{res.note}</p> : null}
      <RawJson data={res} />
    </>
  );
}

function Fba({ res }: { res: FbaResult }) {
  const s = res.summary;
  return (
    <>
      <Stats>
        <Stat k="Wild-type guaranteed yield" v={fmt(s.wild_type_guaranteed_yield)} s="mol/mol" />
        <Stat k="Best guaranteed yield" v={fmt(s.best_guaranteed_yield)} s="mol/mol" />
        <Stat k="Reactions deleted" v={String(s.knockouts.length)} />
        <Stat k="Genes to delete" v={String(s.genes.length)} />
      </Stats>
      {s.improved ? null : <Notice>No deletion raised the guaranteed yield. For a product the cell has no reason to make, try the MOMA method; for one where only a combination helps, try the designed knockouts.</Notice>}
      <Table rows={res.knockouts.map((r, i) => ({ ...r, step: i }))} columns={[
        { key: "step", label: "Step" },
        { key: "knocked_out", label: "Reactions deleted", render: (r) => r.knocked_out.length ? r.knocked_out.join(", ") : "(wild type)" },
        { key: "genes", label: "Genes" },
        { key: "growth", label: "Growth (1/h)" },
        { key: "guaranteed_yield", label: "Guaranteed yield" },
      ]} />
      <Caveats list={res.caveats} />
      <RawJson data={res} />
    </>
  );
}
