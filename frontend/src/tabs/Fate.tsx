import { useState } from "react";
import { useHistory, useRun, useSetup, conditions } from "../state";
import { MetaboliteInput, RunButton, RunOutput } from "../components/run";
import { Caveats, PanelHead, Pill, RawJson, Stat, Stats, Table } from "../components/ui";

interface Consumer { reaction: string; name: string; equation: string; genes: string[]; flux_at_max_growth: number; active: boolean; safe_to_delete: boolean }
interface Result {
  target: string; target_name?: string; formula?: string; charge?: number | null; charge_note?: string;
  in_host_model: boolean; protective_knockouts?: string[]; consumed_by: Consumer[];
  reactive_pools_the_host_makes: Record<string, unknown>[]; functional_group_hints?: Record<string, unknown>[]; caveats?: string[];
}

export function FateTab() {
  const { setup } = useSetup();
  const { timed } = useHistory();
  const [target, setTarget] = useState("succ_c");
  const [donor, setDonor] = useState("");
  const { state, run, busy } = useRun<Result>();

  const go = () => run((onStatus) => timed<Result>("check_product_fate",
    { ...conditions(setup), target: target.trim(), donor: donor.trim() }, onStatus));

  return (
    <>
      <PanelHead title="What the cell does to the product"
        blurb="Which host reactions consume it, which deletions protect it, and whether its charge means it needs a transporter."
        action={<RunButton busy={busy} onClick={go}>Check fate</RunButton>} />
      <div className="grid small-grid">
        <label>Target metabolite <MetaboliteInput value={target} onChange={setTarget} />
          <span className="hint">A metabolite id inside the cell, e.g. succ_c.</span></label>
        <label>Donor model
          <input list="model-list" value={donor} autoComplete="off" onChange={(e) => setDonor(e.target.value)} />
          <span className="hint">Only needed for a product the host can't make yet.</span></label>
      </div>
      <RunOutput state={state} render={(res) => (
        <>
          <Stats>
            <Stat k="Product" v={res.target_name || res.target} s={res.formula || ""} />
            <Stat k="Charge at cytosolic pH" v={res.charge == null ? "–" : String(res.charge)} />
            <Stat k="In host model" v={res.in_host_model ? "yes" : "no"} s={res.in_host_model ? "" : "no native enzyme touches it"} />
            <Stat k="Protective knockouts" v={String((res.protective_knockouts ?? []).length)} />
          </Stats>
          {res.charge_note ? <p className="small">{res.charge_note}</p> : null}
          <h3>Host reactions that consume it</h3>
          <Table rows={res.consumed_by} columns={[
            { key: "reaction", label: "Reaction" },
            { key: "name", label: "Name" },
            { key: "equation", label: "Equation", render: (r) => <span className="mono">{r.equation}</span> },
            { key: "genes", label: "Genes" },
            { key: "flux_at_max_growth", label: "Flux at max growth" },
            { key: "active", label: "Active", render: (r) => r.active ? <Pill kind="bad">yes</Pill> : <Pill>no</Pill> },
            { key: "safe_to_delete", label: "Safe to delete", render: (r) => r.safe_to_delete ? <Pill kind="good">yes</Pill> : <Pill>no</Pill> },
          ]} />
          <h3>Reactive pools the host makes</h3>
          <Table rows={res.reactive_pools_the_host_makes} columns={[
            { key: "metabolite", label: "Metabolite" }, { key: "kind", label: "Kind" },
            { key: "max_flux", label: "Max flux" }, { key: "why_it_matters", label: "Why it matters" },
          ]} />
          {res.functional_group_hints?.length ? <><h3>Functional group hints (guessed from the name)</h3><Table rows={res.functional_group_hints} /></> : null}
          <Caveats list={res.caveats} />
          <RawJson data={res} />
        </>
      )} />
    </>
  );
}
