import type { ParamSchema } from "./api";

/** A number for display: 3 significant decimals, exponent form for tiny or huge values, an en dash for missing. */
export function fmt(x: unknown, digits = 3): string {
  if (x == null || (typeof x === "number" && Number.isNaN(x))) return "–";
  if (typeof x !== "number") return String(x);
  if (x !== 0 && (Math.abs(x) < 1e-3 || Math.abs(x) >= 1e5)) return x.toExponential(2);
  return Number(x.toFixed(digits)).toString();
}

export const pct = (x: number | null | undefined) => (x == null ? "–" : `${fmt(100 * x, 1)}%`);

/** The JSON-schema type of a tool parameter, looking through Optional[...] (anyOf with null). */
export function paramType(p: ParamSchema): string {
  return p.type || p.anyOf?.map((a) => a.type).find((t) => t && t !== "null") || "string";
}

/**
 * Turn raw form values into tool kwargs. Empty fields are left out so the server
 * uses its default; numbers are parsed by schema type.
 */
export function coerceArgs(values: Record<string, string | boolean>, props: Record<string, ParamSchema>): Record<string, unknown> {
  const args: Record<string, unknown> = {};
  for (const [name, raw] of Object.entries(values)) {
    const type = paramType(props[name] || {});
    if (type === "boolean") { args[name] = raw === true; continue; }
    const text = String(raw).trim();
    if (text === "") continue;
    if (type === "integer") args[name] = parseInt(text, 10);
    else if (type === "number") args[name] = Number(text);
    else args[name] = text;
  }
  return args;
}

/** "glucose, glycerol" -> ["glucose", "glycerol"], matching the server's _split. */
export const splitList = (text: string) => text.split(",").map((t) => t.trim()).filter(Boolean);
