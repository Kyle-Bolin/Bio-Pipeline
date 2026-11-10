import type { ReactNode } from "react";
import type { RunState } from "../state";
import { useMetaboliteSuggest } from "../state";
import { ErrorBox } from "./ui";

/**
 * The panel's Run button. Without onClick it submits the surrounding form. With a
 * `problem` (an empty or out-of-range input) it is disabled and says what to fix.
 */
export function RunButton({ busy, onClick, problem, children }: {
  busy: boolean; onClick?: () => void; problem?: string | null; children: ReactNode;
}) {
  return (
    <div className="run-wrap">
      <button className="run" type={onClick ? "button" : "submit"} disabled={busy || !!problem} onClick={onClick}
>{busy ? "Running…" : children}</button>
      {problem ? <span className="hint problem" role="alert">{problem}</span> : null}
    </div>
  );
}

/** Spinner while running, the error, or `render(data)` once done. */
export function RunOutput<T>({ state, render }: { state: RunState<T>; render: (data: T) => ReactNode }) {
  if (state.kind === "idle") return null;
  return (
    <div className="out" aria-live="polite">
      {state.kind === "running" ? (
        <div className="loading"><span className="spinner" />{state.status === "queued" ? "Queued…" : state.status === "reconnecting" ? "Lost contact with the server, retrying…" : "Running…"}</div>
      ) : state.kind === "error" ? (
        <ErrorBox><strong>Couldn't run that. </strong>{state.message}</ErrorBox>
      ) : render(state.data)}
    </div>
  );
}

/** A text input with metabolite suggestions from the current model. */
export function MetaboliteInput({ value, onChange, id }: { value: string; onChange: (v: string) => void; id?: string }) {
  const { suggest } = useMetaboliteSuggest();
  return (
    <input id={id} list="met-list" value={value} autoComplete="off"
      onChange={(e) => { onChange(e.target.value); suggest(e.target.value.trim()); }} />
  );
}

export function NumberInput({ value, onChange, min, max, step }: {
  value: number; onChange: (v: number) => void; min?: number; max?: number; step?: number | "any";
}) {
  return <input type="number" value={Number.isNaN(value) ? "" : value} min={min} max={max} step={step}
    onChange={(e) => onChange(e.target.value === "" ? Number.NaN : Number(e.target.value))} />;
}
