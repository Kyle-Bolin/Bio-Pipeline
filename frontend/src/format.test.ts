import { describe, expect, it } from "vitest";
import { coerceArgs, fmt, paramType, pct, splitList } from "./format";

describe("fmt", () => {
  it("rounds, switches to exponent form at the extremes and shows missing as a dash", () => {
    expect(fmt(0.123456)).toBe("0.123");
    expect(fmt(2)).toBe("2");
    expect(fmt(0)).toBe("0");
    expect(fmt(0.00001234)).toBe("1.23e-5");
    expect(fmt(123456)).toBe("1.23e+5");
    expect(fmt(null)).toBe("–");
    expect(fmt(Number.NaN)).toBe("–");
    expect(fmt("abc")).toBe("abc");
  });
  it("formats fractions as percentages", () => {
    expect(pct(0.5)).toBe("50%");
    expect(pct(null)).toBe("–");
  });
});

describe("tool forms", () => {
  const props = {
    product: { type: "string" },
    limit: { type: "integer", default: 30 },
    uptake: { type: "number", default: 10 },
    aerobic: { type: "boolean", default: true },
    donor: { anyOf: [{ type: "string" }, { type: "null" }] },
  };
  it("reads Optional types", () => {
    expect(paramType(props.donor)).toBe("string");
    expect(paramType({})).toBe("string");
  });
  it("parses by type and leaves empty fields to the server's default", () => {
    expect(coerceArgs({ product: " etoh ", limit: "5", uptake: "", aerobic: false, donor: "" }, props))
      .toEqual({ product: "etoh", limit: 5, aerobic: false });
  });
  it("splits comma lists like the server", () => {
    expect(splitList("glucose, glycerol,,xylose ")).toEqual(["glucose", "glycerol", "xylose"]);
  });
});
