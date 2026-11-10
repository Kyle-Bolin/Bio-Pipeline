import { fmt } from "../format";

export interface EnvelopePoint { growth: number; flux_minimum: number; flux_maximum: number }

/** Production envelope: shaded band between min and max product flux at each growth rate. */
export function EnvelopeChart({ points, marker }: { points: EnvelopePoint[]; marker?: { growth: number; flux: number } }) {
  const W = 640, H = 300, L = 52, R = 14, T = 14, B = 40;
  const ymax = Math.max(...points.map((p) => p.flux_maximum), 1e-9);
  const xmax = Math.max(...points.map((p) => p.growth), 1e-9);
  const X = (x: number) => L + (x / xmax) * (W - L - R);
  const Y = (y: number) => H - B - (y / ymax) * (H - T - B);
  const line = (key: "flux_minimum" | "flux_maximum") =>
    points.map((p, i) => `${i ? "L" : "M"}${X(p.growth).toFixed(1)},${Y(p[key]).toFixed(1)}`).join("");
  const band = line("flux_maximum") +
    [...points].reverse().map((p) => `L${X(p.growth).toFixed(1)},${Y(p.flux_minimum).toFixed(1)}`).join("") + "Z";
  const ticks = [0, 1, 2, 3, 4];
  const midY = (T + H - B) / 2;
  return (
    <figure style={{ margin: 0 }}>
      <svg viewBox={`0 0 ${W} ${H}`} className="chart" role="img" aria-label="Production envelope: product flux range at each growth rate">
        {ticks.map((i) => {
          const yv = (ymax * i) / 4, xv = (xmax * i) / 4;
          return (
            <g key={i}>
              <line x1={L} x2={W - R} y1={Y(yv)} y2={Y(yv)} className="axis" strokeOpacity={i ? 0.5 : 1} />
              <text x={L - 6} y={Y(yv) + 4} textAnchor="end">{fmt(yv, 2)}</text>
              <text x={X(xv)} y={H - B + 16} textAnchor="middle">{fmt(xv, 2)}</text>
            </g>
          );
        })}
        <line x1={L} x2={L} y1={T} y2={H - B} className="axis" />
        <path d={band} className="band" />
        <path d={line("flux_maximum")} className="max" />
        <path d={line("flux_minimum")} className="min" />
        {marker ? <circle cx={X(marker.growth)} cy={Y(marker.flux)} r={4.5} className="dot" /> : null}
        <text x={(L + W - R) / 2} y={H - 6} textAnchor="middle">growth rate (1/h)</text>
        <text x={12} y={midY} textAnchor="middle" transform={`rotate(-90 12 ${midY})`}>product flux (mmol/gDW/h)</text>
      </svg>
      <figcaption className="muted small">
        Shaded: every product flux the cell can have at that growth rate. Solid line is the best case, dashed the
        guaranteed minimum. {marker ? "The dot is the simulated cell growing as fast as it can." : ""}
      </figcaption>
    </figure>
  );
}

export interface BarRow { label: string; value: number | null; ceiling?: number | null; note?: string }

/**
 * Horizontal bars on a shared 0..max scale. `ceiling` draws a faint bar behind the
 * solid one (e.g. theoretical max behind guaranteed). Rows with a null value show their note.
 */
export function BarChart({ rows, format, ariaLabel, legend }: {
  rows: BarRow[]; format: (x: number) => string; ariaLabel: string; legend?: [string, string];
}) {
  const max = Math.max(1e-9, ...rows.flatMap((r) => [r.value ?? 0, r.ceiling ?? 0]));
  return (
    <figure style={{ margin: 0 }}>
      <div className="bars" role="img" aria-label={ariaLabel}>
        {rows.map((r, i) => (
          <div className="bar-row" key={i}>
            <div className="bar-label">{r.label}</div>
            <div className="bar-track">
              {r.value == null ? <span className="muted small">{r.note ?? "–"}</span> : (
                <>
                  {r.ceiling != null ? <div className="bar ceiling" style={{ width: `${(100 * r.ceiling) / max}%` }} /> : null}
                  <div className="bar" style={{ width: `${(100 * r.value) / max}%` }} />
                </>
              )}
            </div>
            <div className="bar-value">{r.value == null ? "" : format(r.value)}</div>
          </div>
        ))}
      </div>
      {legend ? (
        <figcaption className="muted small legend">
          <span><i className="swatch" /> {legend[0]}</span>
          <span><i className="swatch ceiling" /> {legend[1]}</span>
        </figcaption>
      ) : null}
    </figure>
  );
}

/** A signed flux as an inline bar: right of centre is forward/secretion, left is reverse/uptake. */
export function FluxBar({ value, max }: { value: number; max: number }) {
  const w = Math.min(50, (50 * Math.abs(value)) / Math.max(max, 1e-9));
  return (
    <div className="fluxbar" aria-hidden="true">
      <div className={`fill ${value < 0 ? "neg" : "pos"}`} style={value < 0 ? { right: "50%", width: `${w}%` } : { left: "50%", width: `${w}%` }} />
    </div>
  );
}
