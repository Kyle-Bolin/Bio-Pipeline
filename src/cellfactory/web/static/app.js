// cellfactory front end. Talks to the REST API:
//   GET  /api/health, GET /api/tools
//   POST /api/tools/<name>      body = the tool's kwargs, sync
//   POST /api/jobs              body = {tool, args}, then poll GET /api/jobs/<id>
// Every tool in cellfactory.tools.TOOLS is reachable from the "All tools" tab;
// the other tabs are friendlier views over the main workflow.
"use strict";

const API = "/api";
const $ = (sel) => document.querySelector(sel);
const BUNDLED_FAST = new Set(["e_coli_core"]);
// Tools that take seconds to minutes; run as background jobs when the API supports them.
const SLOW = new Set(["search_knockouts", "search_knockouts_moma", "find_precursors", "check_model", "check_product_fate"]);

let jobsSupported = false; // set by a probe of GET /api/jobs at start-up
let toolSpecs = [];

// ---------- DOM helper: never inject HTML strings built from data ----------
function h(tag, attrs = {}, ...kids) {
  const el = document.createElement(tag);
  for (const [k, v] of Object.entries(attrs || {})) {
    if (v == null || v === false) continue;
    if (k === "class") el.className = v;
    else if (k.startsWith("on")) el.addEventListener(k.slice(2), v);
    else el.setAttribute(k, v === true ? "" : v);
  }
  for (const kid of kids.flat(Infinity)) {
    if (kid == null || kid === false) continue;
    el.append(kid instanceof Node ? kid : String(kid));
  }
  return el;
}

// ---------- API ----------
class ApiError extends Error {}

async function request(method, path, body) {
  let res;
  try {
    res = await fetch(API + path, {
      method,
      headers: body === undefined ? {} : { "Content-Type": "application/json" },
      body: body === undefined ? undefined : JSON.stringify(body),
    });
  } catch (err) {
    throw new ApiError("Can't reach the API. Is `cellfactory-api` running?");
  }
  let data = null;
  try { data = await res.json(); } catch { /* empty or non-JSON body */ }
  if (!res.ok) {
    const msg = (data && (data.error || data.detail)) || `${res.status} ${res.statusText}`;
    const err = new ApiError(typeof msg === "string" ? msg : JSON.stringify(msg));
    err.status = res.status;
    throw err;
  }
  return data;
}

const JOB_TIMEOUT_MS = 20 * 60 * 1000; // a genome-scale MOMA search can take several minutes
const sleep = (ms) => new Promise((r) => setTimeout(r, ms));

async function runJob(tool, args, onStatus) {
  const job = await request("POST", "/jobs", { tool, args });
  let delay = 400;
  const deadline = Date.now() + JOB_TIMEOUT_MS;
  for (;;) {
    if (Date.now() > deadline) {
      throw new ApiError(`Still not finished after ${JOB_TIMEOUT_MS / 60000} minutes, so the page stopped waiting. ` +
        `The job may still be running: check ${API}/jobs/${job.id}.`);
    }
    await sleep(delay);
    delay = Math.min(delay * 1.5, 2000);
    const j = await request("GET", `/jobs/${encodeURIComponent(job.id)}`);
    if (onStatus) onStatus(j.status);
    if (j.status === "done") return j.result;
    if (j.status === "failed") throw new ApiError(j.error || "The job failed.");
  }
}

// Call a tool; slow tools (or anything on a genome-scale model) go through /api/jobs.
async function call(tool, args = {}, onStatus) {
  const heavy = SLOW.has(tool) || (args.model && !BUNDLED_FAST.has(args.model));
  if (heavy && jobsSupported) return runJob(tool, args, onStatus);
  return request("POST", `/tools/${encodeURIComponent(tool)}`, args);
}

// ---------- shared setup ----------
const val = (id) => $(id).value.trim();
const num = (id) => Number($(id).value);
function conditions() {
  const args = {
    model: val("#model") || "e_coli_core",
    substrate: val("#substrate") || "glucose",
    uptake: num("#uptake"),
    aerobic: document.querySelector('input[name="aerobic"]:checked').value === "true",
  };
  return args;
}
const pathway = () => $("#pathway").value;

// ---------- formatting ----------
function fmt(x, digits = 3) {
  if (x == null || Number.isNaN(x)) return "–";
  if (typeof x !== "number") return String(x);
  if (x !== 0 && (Math.abs(x) < 1e-3 || Math.abs(x) >= 1e5)) return x.toExponential(2);
  return Number(x.toFixed(digits)).toString();
}
const pct = (x) => (x == null ? "–" : `${fmt(100 * x, 1)}%`);

function stat(k, v, s) {
  return h("div", { class: "stat" }, h("div", { class: "k" }, k), h("div", { class: "v" }, v), s ? h("div", { class: "s" }, s) : null);
}

function table(rows, columns) {
  if (!rows || !rows.length) return h("p", { class: "muted small" }, "None.");
  columns = columns || Object.keys(rows[0]).map((k) => ({ key: k, label: k }));
  return h("div", { class: "table-wrap" },
    h("table", {},
      h("thead", {}, h("tr", {}, columns.map((c) => h("th", {}, c.label)))),
      h("tbody", {}, rows.map((r) => h("tr", {}, columns.map((c) => {
        const v = c.render ? c.render(r) : r[c.key];
        const isNum = typeof r[c.key] === "number" && !c.render;
        return h("td", { class: isNum ? "num" : null }, v instanceof Node ? v : isNum ? fmt(v) : Array.isArray(v) ? v.join(", ") : v == null ? "" : String(v));
      }))))));
}

function rawJson(data) {
  return h("details", { class: "raw" }, h("summary", {}, "Raw JSON"), h("pre", {}, JSON.stringify(data, null, 2)));
}

function caveats(list) {
  if (!list || !list.length) return null;
  return h("div", { class: "notice" }, h("strong", {}, "Caveats"), h("ul", {}, list.map((c) => h("li", {}, c))));
}

function pill(text, kind) { return h("span", { class: `pill ${kind || ""}` }, text); }

// Production envelope: shaded band between min and max product flux at each growth rate.
function envelopeChart(points, marker) {
  const W = 640, H = 300, L = 52, R = 14, T = 14, B = 40;
  const xs = points.map((p) => p.growth), ymax = Math.max(...points.map((p) => p.flux_maximum), 1e-9);
  const xmax = Math.max(...xs, 1e-9);
  const X = (x) => L + (x / xmax) * (W - L - R);
  const Y = (y) => H - B - (y / ymax) * (H - T - B);
  const line = (key) => points.map((p, i) => `${i ? "L" : "M"}${X(p.growth).toFixed(1)},${Y(p[key]).toFixed(1)}`).join("");
  const band = line("flux_maximum") + [...points].reverse().map((p) => `L${X(p.growth).toFixed(1)},${Y(p.flux_minimum).toFixed(1)}`).join("") + "Z";
  const ns = "http://www.w3.org/2000/svg";
  const s = (tag, attrs, text) => {
    const el = document.createElementNS(ns, tag);
    for (const [k, v] of Object.entries(attrs)) el.setAttribute(k, v);
    if (text != null) el.textContent = text;
    return el;
  };
  const svg = s("svg", { viewBox: `0 0 ${W} ${H}`, class: "chart", role: "img",
    "aria-label": "Production envelope: product flux range at each growth rate" });
  for (let i = 0; i <= 4; i++) {
    const yv = (ymax * i) / 4, xv = (xmax * i) / 4;
    svg.append(s("line", { x1: L, x2: W - R, y1: Y(yv), y2: Y(yv), class: "axis", "stroke-opacity": i ? 0.5 : 1 }));
    svg.append(s("text", { x: L - 6, y: Y(yv) + 4, "text-anchor": "end" }, fmt(yv, 2)));
    svg.append(s("text", { x: X(xv), y: H - B + 16, "text-anchor": "middle" }, fmt(xv, 2)));
  }
  svg.append(s("line", { x1: L, x2: L, y1: T, y2: H - B, class: "axis" }));
  svg.append(s("path", { d: band, class: "band" }));
  svg.append(s("path", { d: line("flux_maximum"), class: "max" }));
  svg.append(s("path", { d: line("flux_minimum"), class: "min" }));
  if (marker) svg.append(s("circle", { cx: X(marker.growth), cy: Y(marker.flux), r: 4.5, class: "dot" }));
  svg.append(s("text", { x: (L + W - R) / 2, y: H - 6, "text-anchor": "middle" }, "growth rate (1/h)"));
  svg.append(s("text", { x: 12, y: (T + H - B) / 2, "text-anchor": "middle", transform: `rotate(-90 12 ${(T + H - B) / 2})` }, "product flux (mmol/gDW/h)"));
  return h("figure", { style: "margin:0" }, svg,
    h("figcaption", { class: "muted small" },
      "Shaded: every product flux the cell can have at that growth rate. Solid line is the best case, dashed the guaranteed minimum. ",
      marker ? "The dot is the simulated cell growing as fast as it can." : ""));
}

// ---------- renderers ----------
const render = {
  yield([sim, env]) {
    const p = sim.production;
    const sub = sim.substrate, prod = sim.product;
    return [
      h("div", { class: "stats" },
        stat("Growth rate", fmt(p.growth), "1/h"),
        stat("Guaranteed yield", fmt(p.guaranteed_yield), "mol product / mol substrate, at max growth"),
        stat("Best-case yield", fmt(p.best_case_yield), "at max growth"),
        stat("Theoretical max", fmt(p.max_theoretical_yield), "ignoring growth"),
        stat("Carbon yield", pct(p.carbon_yield), "of substrate carbon, guaranteed")),
      h("p", { class: "muted small" }, `${sim.model_id}: ${sub} → ${prod}, ${sim.aerobic ? "aerobic" : "anaerobic"}, uptake ${fmt(sim.uptake)} mmol/gDW/h.`),
      sim.pathway ? h("div", {}, h("h3", {}, `Added pathway: ${sim.pathway.name}`),
        h("p", { class: "small" }, "Genes cloned in: ", sim.pathway.genes.join(", "))) : null,
      env && env.points && env.points.length ? envelopeChart(env.points, { growth: p.growth, flux: p.target_flux_min }) : null,
      caveats(sim.caveats),
      ...checksTables(sim.checks),
      rawJson({ simulate_production: sim, production_envelope: env }),
    ];
  },

  knockouts(res) {
    if (res.steps) { // MOMA
      const s = res.summary;
      return [
        h("div", { class: "stats" },
          stat("Parent yield", fmt(s.parent_yield), "mol/mol"),
          stat("Best yield", fmt(s.best_yield), "mol/mol"),
          stat("Genes deleted", String(s.genes.length))),
        s.improved ? null : h("div", { class: "notice" }, "No deletion improved the predicted yield."),
        table(res.steps.map((r, i) => ({ ...r, step: i })), [
          { key: "step", label: "Step" },
          { key: "knocked_out", label: "Genes deleted", render: (r) => r.knocked_out.length ? r.knocked_out.map((g, i) => `${r.gene_names[i] || g} (${g})`).join(", ") : "(parent)" },
          { key: "reactions", label: "Reactions disabled" },
          { key: "growth", label: "Growth (1/h)" },
          { key: "product_flux", label: "Product flux" },
          { key: "product_yield", label: "Yield (mol/mol)" },
        ]),
        h("p", { class: "muted small" }, res.note),
        rawJson(res),
      ];
    }
    const s = res.summary;
    return [
      h("div", { class: "stats" },
        stat("Wild-type guaranteed yield", fmt(s.wild_type_guaranteed_yield), "mol/mol"),
        stat("Best guaranteed yield", fmt(s.best_guaranteed_yield), "mol/mol"),
        stat("Reactions deleted", String(s.knockouts.length)),
        stat("Genes to delete", String(s.genes.length))),
      s.improved ? null : h("div", { class: "notice" },
        "No deletion raised the guaranteed yield. For a product the cell has no reason to make, try the MOMA method."),
      table(res.knockouts.map((r, i) => ({ ...r, step: i })), [
        { key: "step", label: "Step" },
        { key: "knocked_out", label: "Reactions deleted", render: (r) => r.knocked_out.length ? r.knocked_out.join(", ") : "(wild type)" },
        { key: "genes", label: "Genes" },
        { key: "growth", label: "Growth (1/h)" },
        { key: "guaranteed_yield", label: "Guaranteed yield" },
      ]),
      caveats(res.caveats),
      rawJson(res),
    ];
  },

  design(res) {
    const routes = res.routes || [];
    return [
      h("p", {}, routes.length ? `${routes.length} route${routes.length > 1 ? "s" : ""} from ${res.model} metabolism to ${res.target} (donor: ${res.donor}).` : res.summary),
      ...routes.map((r, i) => h("div", { class: "route" },
        h("h3", { style: "margin-top:0" }, `Route ${i + 1}: ${r.precursor_name || r.precursor} → ${r.target_name || r.target}`),
        h("p", { class: "small muted" }, `${r.heterologous_steps} new enzyme${r.heterologous_steps === 1 ? "" : "s"}, ${r.total_steps} steps in total.`,
          r.precursor_max_flux != null ? ` The host can make up to ${fmt(r.precursor_max_flux)} mmol/gDW/h of the starting metabolite.` : ""),
        h("ol", {}, r.steps.map((st) => h("li", {},
          h("strong", {}, st.name || st.reaction), " ", st.native ? pill("host has it", "good") : pill("clone in", "new"),
          h("div", { class: "mono" }, st.equation),
          st.genes && st.genes.length ? h("div", { class: "muted" }, "Genes: ", st.genes.join(", ")) : null,
          st.also_needs && st.also_needs.length ? h("div", { class: "muted" }, "Also needs: ", st.also_needs.join(", ")) : null))))),
      h("p", { class: "muted small" }, res.note),
      rawJson(res),
    ];
  },

  fate(res) {
    return [
      h("div", { class: "stats" },
        stat("Product", res.target_name || res.target, res.formula || ""),
        stat("Charge at cytosolic pH", res.charge == null ? "–" : String(res.charge)),
        stat("In host model", res.in_host_model ? "yes" : "no", res.in_host_model ? "" : "no native enzyme touches it"),
        stat("Protective knockouts", String((res.protective_knockouts || []).length))),
      res.charge_note ? h("p", { class: "small" }, res.charge_note) : null,
      h("h3", {}, "Host reactions that consume it"),
      table(res.consumed_by, [
        { key: "reaction", label: "Reaction" },
        { key: "name", label: "Name" },
        { key: "equation", label: "Equation", render: (r) => h("span", { class: "mono" }, r.equation) },
        { key: "genes", label: "Genes" },
        { key: "flux_at_max_growth", label: "Flux at max growth" },
        { key: "active", label: "Active", render: (r) => r.active ? pill("yes", "bad") : pill("no") },
        { key: "safe_to_delete", label: "Safe to delete", render: (r) => r.safe_to_delete ? pill("yes", "good") : pill("no") },
      ]),
      h("h3", {}, "Reactive pools the host makes"),
      table(res.reactive_pools_the_host_makes, [
        { key: "metabolite", label: "Metabolite" }, { key: "kind", label: "Kind" },
        { key: "max_flux", label: "Max flux" }, { key: "why_it_matters", label: "Why it matters" },
      ]),
      res.functional_group_hints && res.functional_group_hints.length ? [h("h3", {}, "Functional group hints (guessed from the name)"), table(res.functional_group_hints)] : null,
      caveats(res.caveats),
      rawJson(res),
    ];
  },

  check(res) {
    const mb = res.mass_imbalanced_reactions, eg = res.predicted_essential_genes;
    return [
      h("div", { class: res.ok ? "notice" : "error" },
        res.ok ? "Passes the free-energy checks: the model can't grow or make ATP from nothing."
          : h("span", {}, h("strong", {}, "Fails a free-energy check. "), "The model makes something from nothing, which inflates every yield.")),
      h("h3", {}, "Free-energy checks"),
      table(res.free_energy_checks, [
        { key: "check", label: "Check" }, { key: "value", label: "Value" },
        { key: "ok", label: "Result", render: (r) => r.ok ? pill("pass", "good") : pill("fail", "bad") },
      ]),
      h("h3", {}, `Mass-imbalanced reactions (${mb.total}${mb.truncated ? ", first 50 shown" : ""})`),
      table(mb.rows),
      h("h3", {}, `Predicted essential genes (${eg.total}${eg.truncated ? ", first 50 shown" : ""})`),
      table(eg.rows),
      rawJson(res),
    ];
  },
};

function checksTables(checks) {
  if (!checks) return [];
  return Object.entries(checks).filter(([, rows]) => rows && rows.length).map(([title, rows]) =>
    h("details", { class: "raw" }, h("summary", {}, title), table(rows)));
}

// ---------- running an analysis ----------
const history = [];
function logRun(tool, args, ok, ms) {
  history.unshift({ tool, args, ok, ms, at: new Date() });
  history.length = Math.min(history.length, 12);
  $("#history").replaceChildren(...history.map((e) => h("li", {},
    h("span", { class: "when" }, e.at.toLocaleTimeString()),
    h("code", {}, e.tool), e.ok ? pill(`${(e.ms / 1000).toFixed(1)} s`, "good") : pill("error", "bad"),
    h("span", { class: "muted" }, Object.entries(e.args).filter(([, v]) => v !== "" && v != null).map(([k, v]) => `${k}=${v}`).join(" ")))));
}

async function withOutput(outSel, button, work) {
  const out = $(outSel);
  const label = h("span", {}, "Running…");
  out.replaceChildren(h("div", { class: "loading" }, h("span", { class: "spinner" }), label));
  if (button) button.disabled = true;
  try {
    const nodes = await work((status) => { label.textContent = status === "queued" ? "Queued…" : "Running…"; });
    out.replaceChildren(...[nodes].flat(Infinity).filter(Boolean));
  } catch (err) {
    out.replaceChildren(h("div", { class: "error" }, h("strong", {}, "Couldn't run that. "), err.message));
  } finally {
    if (button) button.disabled = false;
  }
}

async function timed(tool, args, onStatus) {
  const t0 = performance.now();
  try {
    const r = await call(tool, args, onStatus);
    logRun(tool, args, true, performance.now() - t0);
    return r;
  } catch (err) {
    logRun(tool, args, false, performance.now() - t0);
    throw err;
  }
}

const runners = {
  async yield(onStatus) {
    const args = { ...conditions(), product: val("#product"), pathway: pathway() };
    const [sim, env] = await Promise.all([
      timed("simulate_production", args, onStatus),
      timed("production_envelope", { ...args, points: 25 }, onStatus).catch(() => null),
    ]);
    return render.yield([sim, env]);
  },
  async knockouts(onStatus) {
    const c = conditions(), max = num("#ko-max");
    if ($("#ko-method").value === "moma") {
      return render.knockouts(await timed("search_knockouts_moma",
        { ...c, product: val("#product"), pathway: pathway(), max_knockouts: max }, onStatus));
    }
    return render.knockouts(await timed("search_knockouts",
      { ...c, product: val("#product"), pathway: pathway(), max_knockouts: max, min_growth: num("#ko-min-growth") }, onStatus));
  },
  async design(onStatus) {
    return render.design(await timed("find_precursors", {
      ...conditions(), target: val("#design-target"), donor: val("#design-donor"),
      max_new_enzymes: num("#design-max-new"), max_total_steps: num("#design-max-steps"),
    }, onStatus));
  },
  async fate(onStatus) {
    return render.fate(await timed("check_product_fate",
      { ...conditions(), target: val("#fate-target"), donor: val("#fate-donor") }, onStatus));
  },
  async check(onStatus) {
    const { model, substrate, aerobic } = conditions();
    return render.check(await timed("check_model", { model, substrate, aerobic }, onStatus));
  },
};

document.querySelectorAll("button[data-run]").forEach((btn) => {
  const name = btn.dataset.run;
  btn.addEventListener("click", () => withOutput(`#out-${name}`, btn, runners[name]));
});

// ---------- tabs ----------
function showTab(name) {
  document.querySelectorAll('[role="tab"]').forEach((t) => t.setAttribute("aria-selected", String(t.dataset.tab === name)));
  document.querySelectorAll("[data-panel]").forEach((p) => { p.hidden = p.dataset.panel !== name; });
  try { localStorage.setItem("cellfactory.tab", name); } catch { /* storage blocked */ }
}
document.querySelectorAll('[role="tab"]').forEach((t) => t.addEventListener("click", () => showTab(t.dataset.tab)));
$("#ko-method").addEventListener("change", () => {
  const moma = $("#ko-method").value === "moma";
  $("#ko-min-growth-l").hidden = moma;
  $("#ko-max").max = moma ? 5 : 10;
  if (moma && num("#ko-max") > 5) $("#ko-max").value = 2;
});

// ---------- metabolite suggestions ----------
let metTimer = null, lastQuery = "";
function suggestMetabolites(q) {
  clearTimeout(metTimer);
  if (q.length < 3 || q === lastQuery) return;
  metTimer = setTimeout(async () => {
    lastQuery = q;
    try {
      const r = await request("POST", "/tools/search_metabolites", { query: q, model: val("#model") || "e_coli_core", limit: 25 });
      $("#met-list").replaceChildren(...r.metabolites.map((m) => h("option", { value: m.id }, `${m.name}${m.compartment ? ` [${m.compartment}]` : ""}`)));
    } catch { /* suggestions are optional */ }
  }, 300);
}
["#substrate", "#product", "#design-target", "#fate-target"].forEach((sel) =>
  $(sel).addEventListener("input", (e) => suggestMetabolites(e.target.value.trim())));

// ---------- model info ----------
let modelTimer = null;
async function describeModel() {
  const id = val("#model");
  const info = $("#model-info");
  if (!id) return;
  if (!BUNDLED_FAST.has(id) && !document.querySelector(`#model-list option[value="${CSS.escape(id)}"]`)) {
    info.textContent = "Not bundled: downloaded from BiGG on first use, if the server can reach it.";
    return;
  }
  info.textContent = "Loading model…";
  try {
    const d = await call("describe_model", { model: id });
    info.textContent = `${d.name || d.model}: ${d.reactions} reactions, ${d.metabolites} metabolites, ${d.genes} genes. Grows at ${fmt(d.growth_rate_default_medium)} 1/h on its default medium.`;
  } catch (err) {
    info.textContent = err.message;
  }
}
$("#model").addEventListener("change", () => { clearTimeout(modelTimer); modelTimer = setTimeout(describeModel, 200); lastQuery = ""; });

// ---------- "All tools": forms generated from the API's JSON schemas ----------
function buildToolForm(spec) {
  $("#tool-doc").textContent = spec.description || "";
  const props = (spec.parameters && spec.parameters.properties) || {};
  const required = new Set((spec.parameters && spec.parameters.required) || []);
  $("#tool-form").replaceChildren(...Object.entries(props).map(([name, p]) => {
    const id = `tf-${name}`;
    const type = p.type || (p.anyOf && p.anyOf.map((a) => a.type).find((t) => t !== "null")) || "string";
    if (type === "boolean") {
      return h("label", { class: "check" }, h("input", { type: "checkbox", id, "data-name": name, "data-type": type, checked: p.default === true }), name);
    }
    const input = h("input", {
      id, "data-name": name, "data-type": type, type: type === "integer" || type === "number" ? "number" : "text",
      step: type === "number" ? "any" : null, required: required.has(name),
      list: /model|donor/.test(name) ? "model-list" : /product|substrate|target|query|metabolite|compound/.test(name) ? "met-list" : null,
      placeholder: p.default === undefined ? (required.has(name) ? "required" : "") : String(p.default),
    });
    if (p.default !== undefined && p.default !== "") input.value = String(p.default);
    return h("label", { for: id }, h("span", {}, name, required.has(name) ? h("span", { class: "unit" }, " required") : null), input);
  }));
}

function toolArgs() {
  const args = {};
  $("#tool-form").querySelectorAll("[data-name]").forEach((el) => {
    const { name, type } = el.dataset;
    if (type === "boolean") args[name] = el.checked;
    else if (el.value.trim() === "") return; // let the server use its default
    else if (type === "integer") args[name] = parseInt(el.value, 10);
    else if (type === "number") args[name] = Number(el.value);
    else args[name] = el.value.trim();
  });
  return args;
}

$("#tool-pick").addEventListener("change", () => {
  const spec = toolSpecs.find((t) => t.name === $("#tool-pick").value);
  if (spec) buildToolForm(spec);
  $("#out-tools").replaceChildren();
});
$("#tool-run").addEventListener("click", () => {
  const tool = $("#tool-pick").value;
  if (!tool) return;
  withOutput("#out-tools", $("#tool-run"), async (onStatus) => {
    const r = await timed(tool, toolArgs(), onStatus);
    return h("pre", {}, JSON.stringify(r, null, 2));
  });
});

// ---------- start-up ----------
async function init() {
  let saved = null;
  try { saved = localStorage.getItem("cellfactory.tab"); } catch { /* storage blocked */ }
  if (saved && document.querySelector(`[data-tab="${CSS.escape(saved)}"]`)) showTab(saved);

  const status = $("#api-status");
  try {
    const health = await request("GET", "/health");
    status.textContent = `API ${health.version ? "v" + health.version : "connected"}`;
    status.className = "status ok";
  } catch (err) {
    status.textContent = "API unreachable";
    status.className = "status bad";
    status.title = err.message;
  }

  try { await request("GET", "/jobs"); jobsSupported = true; } catch { jobsSupported = false; }

  const [models, pathways, specs] = await Promise.allSettled([
    request("POST", "/tools/list_models", {}),
    request("POST", "/tools/list_pathways", {}),
    request("GET", "/tools"),
  ]);
  if (models.status === "fulfilled") {
    $("#model-list").replaceChildren(...models.value.bundled.map((m) => h("option", { value: m.id }, m.description)));
  }
  if (pathways.status === "fulfilled") {
    const byName = {};
    for (const p of pathways.value.pathways) {
      byName[p.name] = p;
      $("#pathway").append(h("option", { value: p.name }, `${p.name} (${p.steps.length} enzymes)`));
    }
    $("#pathway").addEventListener("change", () => {
      const p = byName[pathway()];
      $("#pathway-info").textContent = p ? `${p.description} Needs a genome-scale model such as iJO1366; set Product to "${p.product}".` : "";
    });
  }
  if (specs.status === "fulfilled") {
    toolSpecs = Array.isArray(specs.value) ? specs.value : specs.value.tools || [];
    for (const t of toolSpecs) if (t.slow) SLOW.add(t.name); // the API marks its slow tools
    $("#tool-pick").replaceChildren(...toolSpecs.map((t) => h("option", { value: t.name }, t.name)));
    if (toolSpecs.length) buildToolForm(toolSpecs[0]);
  } else {
    $("#tool-doc").textContent = "Couldn't load the tool list from the API.";
  }
  describeModel();
}
init();
