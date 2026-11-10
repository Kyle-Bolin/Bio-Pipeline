import { useState } from "react";
import { useHistory, useRun, useSetup, conditions } from "../state";
import { NumberInput, RunButton, RunOutput } from "../components/run";
import { Caveats, Notice, PanelHead, RawJson, Stat, Stats, Table } from "../components/ui";
import { fmt } from "../format";

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
type Result = FbaResult | MomaResult;

export function KnockoutsTab() {
  const { setup } = useSetup();
  const { timed } = useHistory();
  const [method, setMethod] = useState<"fba" | "moma">("fba");
  const [max, setMax] = useState(3);
  const [minGrowth, setMinGrowth] = useState(0.05);
  const { state, run, busy } = useRun<Result>();

  const go = () => run((onStatus) => {
    const base = { ...conditions(setup), product: setup.product.trim(), pathway: setup.pathway, max_knockouts: max };
    return method === "moma"
      ? timed<MomaResult>("search_knockouts_moma", base, onStatus)
      : timed<FbaResult>("search_knockouts", { ...base, min_growth: minGrowth }, onStatus);
  });

  return (
    <>
      <PanelHead title="Gene knockouts"
        blurb="Search for deletions that force the cell to make more product. FBA finds growth-coupled designs; MOMA suits products the cell has no reason to make."
        action={<RunButton busy={busy} onClick={go}>Search</RunButton>} />
      <div className="grid small-grid">
        <label>Method
          <select value={method} onChange={(e) => {
            const m = e.target.value as "fba" | "moma";
            setMethod(m);
            if (m === "moma" && max > 5) setMax(2);
          }}>
            <option value="fba">FBA (growth-coupled)</option>
            <option value="moma">MOMA (non-native products)</option>
          </select>
        </label>
        <label>Max deletions <NumberInput value={max} onChange={setMax} min={1} max={method === "moma" ? 5 : 10} /></label>
        {method === "fba" ? (
          <label><span>Min growth <span className="unit">1/h</span></span>
            <NumberInput value={minGrowth} onChange={setMinGrowth} min={0} step={0.01} /></label>
        ) : null}
      </div>
      <RunOutput state={state} render={(res) => "steps" in res ? <Moma res={res} /> : <Fba res={res} />} />
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
      {s.improved ? null : <Notice>No deletion raised the guaranteed yield. For a product the cell has no reason to make, try the MOMA method.</Notice>}
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
