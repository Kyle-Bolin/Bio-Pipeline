import { useState } from "react";
import { useHistory, useRun, useSetup, conditions } from "../state";
import { RunButton, RunOutput } from "../components/run";
import { Caveats, ErrorBox, Notice, PanelHead, Pill, RawJson, Stat, Stats, Table, countLabel, type Column, type RowBlock } from "../components/ui";

interface CheckResult {
  ok: boolean;
  free_energy_checks: { check: string; value: unknown; ok: boolean }[];
  mass_imbalanced_reactions: RowBlock;
  predicted_essential_genes: RowBlock;
}
interface QaRow { test: string; parameter: string; title: string; memote_result: string; metric: unknown; gate: string; note: string }
interface QaResult {
  model: string; ok: boolean; engine: string; cached: boolean; counts: Record<string, number>;
  hard_failures: RowBlock<QaRow>; warnings: RowBlock<QaRow>; info: RowBlock<QaRow>; caveats?: string[];
}

const QA_COLUMNS: Column<QaRow>[] = [
  { key: "title", label: "Check" },
  { key: "parameter", label: "On" },
  { key: "memote_result", label: "Result" },
  { key: "metric", label: "Metric" },
  { key: "note", label: "Note" },
];

export function CheckTab() {
  const { setup } = useSetup();
  const { timed } = useHistory();
  const checks = useRun<CheckResult>();
  const qa = useRun<QaResult>();
  const [full, setFull] = useState(false);

  const runChecks = () => checks.run((onStatus) => {
    const { model, substrate, aerobic } = conditions(setup);
    return timed<CheckResult>("check_model", { model, substrate, aerobic }, onStatus);
  });
  const runQa = () => qa.run((onStatus) => timed<QaResult>("qa_model", { model: conditions(setup).model, full }, onStatus));

  return (
    <>
      <PanelHead title="Model sanity checks"
        blurb="Before trusting a striking result: does the model make energy from nothing, are reactions mass-balanced, which genes are essential?"
        action={<RunButton busy={checks.busy} onClick={runChecks}>Run checks</RunButton>} />
      <RunOutput state={checks.state} render={(res) => (
        <>
          {res.ok
            ? <Notice>Passes the free-energy checks: the model can't grow or make ATP from nothing.</Notice>
            : <ErrorBox><strong>Fails a free-energy check. </strong>The model makes something from nothing, which inflates every yield.</ErrorBox>}
          <h3>Free-energy checks</h3>
          <Table rows={res.free_energy_checks} columns={[
            { key: "check", label: "Check" }, { key: "value", label: "Value" },
            { key: "ok", label: "Result", render: (r) => r.ok ? <Pill kind="good">pass</Pill> : <Pill kind="bad">fail</Pill> },
          ]} />
          <h3>Mass-imbalanced reactions {countLabel(res.mass_imbalanced_reactions)}</h3>
          <Table rows={res.mass_imbalanced_reactions.rows} />
          <h3>Predicted essential genes {countLabel(res.predicted_essential_genes)}</h3>
          <Table rows={res.predicted_essential_genes.rows} />
          <RawJson data={res} />
        </>
      )} />

      <hr className="rule" />
      <PanelHead title="MEMOTE quality report"
        blurb="The community test suite for metabolic models. Hard failures (mass or energy from nothing) make every yield suspect, and the yield tools refuse such a model."
        action={<RunButton busy={qa.busy} onClick={runQa}>Run MEMOTE</RunButton>} />
      <div className="grid small-grid">
        <label className="check"><input type="checkbox" checked={full} onChange={(e) => setFull(e.target.checked)} />
          Full suite <span className="hint">minutes on a genome-scale model</span></label>
      </div>
      <RunOutput state={qa.state} render={(res) => (
        <>
          {res.ok
            ? <Notice>No hard failures: {res.model} passes the gate the yield tools apply.</Notice>
            : <ErrorBox><strong>Hard failures. </strong>The yield tools will refuse {res.model} until these are fixed.</ErrorBox>}
          <Stats>
            <Stat k="Hard failures" v={String(res.hard_failures.total)} />
            <Stat k="Warnings" v={String(res.warnings.total)} />
            <Stat k="Passed" v={String(res.counts.pass ?? 0)} />
            <Stat k="Engine" v={res.engine} s={res.cached ? "cached result" : "fresh run"} />
          </Stats>
          {res.hard_failures.total ? <><h3>Hard failures {countLabel(res.hard_failures)}</h3><Table rows={res.hard_failures.rows} columns={QA_COLUMNS} /></> : null}
          <h3>Warnings {countLabel(res.warnings)}</h3>
          <Table rows={res.warnings.rows} columns={QA_COLUMNS} />
          {res.info.total ? (
            <details className="raw"><summary>Info {countLabel(res.info, 20)}</summary><Table rows={res.info.rows} columns={QA_COLUMNS} /></details>
          ) : null}
          <Caveats list={res.caveats} />
          <RawJson data={res} />
        </>
      )} />
    </>
  );
}
