import { useMemo, useState } from "react";
import type { ToolSpec } from "../api";
import { useHistory, useRun } from "../state";
import { RunButton, RunOutput } from "../components/run";
import { PanelHead } from "../components/ui";
import { coerceArgs, paramType } from "../format";

/** Every tool the API exposes, with a form built from its JSON schema. New tools appear here with no UI change. */
export function ToolsTab({ specs, error }: { specs: ToolSpec[]; error: string | null }) {
  const [name, setName] = useState("");
  const spec = specs.find((t) => t.name === name) ?? specs[0];
  return (
    <>
      <PanelHead title="All tools" blurb="Every function the API exposes, with a form built from its parameters. Results are shown as JSON." />
      {error ? <p className="muted">{error}</p> : null}
      <div className="grid small-grid">
        <label>Tool
          <select value={spec?.name ?? ""} onChange={(e) => setName(e.target.value)}>
            {specs.map((t) => <option key={t.name} value={t.name}>{t.name}{t.slow ? " (slow)" : ""}</option>)}
          </select>
        </label>
      </div>
      {spec ? <ToolForm key={spec.name} spec={spec} /> : null}
    </>
  );
}

function initialValues(spec: ToolSpec) {
  const props = spec.parameters?.properties ?? {};
  return Object.fromEntries(Object.entries(props).map(([k, p]) => [k,
    paramType(p) === "boolean" ? p.default === true : p.default == null ? "" : String(p.default)])) as Record<string, string | boolean>;
}

function ToolForm({ spec }: { spec: ToolSpec }) {
  const { timed } = useHistory();
  const props = spec.parameters?.properties ?? {};
  const required = useMemo(() => new Set(spec.parameters?.required ?? []), [spec]);
  const [values, setValues] = useState(() => initialValues(spec));
  const { state, run, busy } = useRun<unknown>();
  const set = (k: string, v: string | boolean) => setValues((s) => ({ ...s, [k]: v }));

  return (
    <>
      <p className="muted doc">{spec.description}</p>
      <form className="grid small-grid" onSubmit={(e) => { e.preventDefault(); run((onStatus) => timed(spec.name, coerceArgs(values, props), onStatus)); }}>
        {Object.entries(props).map(([k, p]) => {
          const type = paramType(p);
          if (type === "boolean") {
            return <label className="check" key={k}><input type="checkbox" checked={values[k] === true} onChange={(e) => set(k, e.target.checked)} />{k}</label>;
          }
          const list = /model|donor/.test(k) ? "model-list" : /product|substrate|target|query|metabolite|compound/.test(k) ? "met-list" : undefined;
          return (
            <label key={k}>
              <span>{k}{required.has(k) ? <span className="unit"> required</span> : null}</span>
              <input value={String(values[k] ?? "")} list={list} required={required.has(k)} autoComplete="off"
                type={type === "integer" || type === "number" ? "number" : "text"} step={type === "number" ? "any" : undefined}
                placeholder={p.default === undefined ? (required.has(k) ? "required" : "") : String(p.default)}
                onChange={(e) => set(k, e.target.value)} />
            </label>
          );
        })}
        <div className="form-actions"><RunButton busy={busy}>Run tool</RunButton></div>
      </form>
      <RunOutput state={state} render={(data) => <pre>{JSON.stringify(data, null, 2)}</pre>} />
    </>
  );
}
