import type { ReactNode } from "react";
import { fmt } from "../format";

export function Stat({ k, v, s }: { k: string; v: ReactNode; s?: ReactNode }) {
  return (
    <div className="stat">
      <div className="k">{k}</div>
      <div className="v">{v}</div>
      {s ? <div className="s">{s}</div> : null}
    </div>
  );
}

export const Stats = ({ children }: { children: ReactNode }) => <div className="stats">{children}</div>;

export function Pill({ kind, children }: { kind?: "good" | "bad" | "new"; children: ReactNode }) {
  return <span className={`pill ${kind ?? ""}`}>{children}</span>;
}

export function Notice({ children }: { children: ReactNode }) {
  return <div className="notice">{children}</div>;
}

export function ErrorBox({ children }: { children: ReactNode }) {
  return <div className="error">{children}</div>;
}

export function Caveats({ list }: { list?: string[] | null }) {
  if (!list || !list.length) return null;
  return (
    <Notice>
      <strong>Caveats</strong>
      <ul>{list.map((c, i) => <li key={i}>{c}</li>)}</ul>
    </Notice>
  );
}

export function RawJson({ data }: { data: unknown }) {
  return (
    <details className="raw">
      <summary>Raw JSON</summary>
      <pre>{JSON.stringify(data, null, 2)}</pre>
    </details>
  );
}

export interface Column<R> {
  key: string;
  label: string;
  render?: (row: R) => ReactNode;
}

type Row = Record<string, any>;

function cell(v: unknown): ReactNode {
  if (typeof v === "number") return fmt(v);
  if (Array.isArray(v)) return v.join(", ");
  if (v == null) return "";
  if (typeof v === "boolean") return v ? "yes" : "no";
  return String(v);
}

/** A table of API rows. Without columns, shows every key of the first row. */
export function Table<R extends Row>({ rows, columns }: { rows?: R[] | null; columns?: Column<R>[] }) {
  if (!rows || !rows.length) return <p className="muted small">None.</p>;
  const cols: Column<R>[] = columns ?? Object.keys(rows[0]).map((k) => ({ key: k, label: k }));
  return (
    <div className="table-wrap">
      <table>
        <thead><tr>{cols.map((c) => <th key={c.key}>{c.label}</th>)}</tr></thead>
        <tbody>
          {rows.map((r, i) => (
            <tr key={i}>
              {cols.map((c) => {
                const isNum = !c.render && typeof r[c.key] === "number";
                return <td key={c.key} className={isNum ? "num" : undefined}>{c.render ? c.render(r) : cell(r[c.key])}</td>;
              })}
            </tr>
          ))}
        </tbody>
      </table>
    </div>
  );
}

/** A {rows, total, truncated} block from the API's _rows helper. */
export interface RowBlock<R = Row> { rows: R[]; total: number; truncated: boolean }

export function countLabel(block: RowBlock | undefined, shown = 50) {
  if (!block) return "";
  return `(${block.total}${block.truncated ? `, first ${shown} shown` : ""})`;
}

export function PanelHead({ title, blurb, action }: { title: string; blurb: ReactNode; action?: ReactNode }) {
  return (
    <div className="panel-head">
      <div>
        <h2>{title}</h2>
        <p className="muted">{blurb}</p>
      </div>
      {action}
    </div>
  );
}
