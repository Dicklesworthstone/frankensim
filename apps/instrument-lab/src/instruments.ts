// The #instruments section: one card per instrument family, built from public/data/instruments/ (tools/build_instruments.py).
// Every number on a card comes from that data; this file only lays it out. A family with no render says so.
import "./instruments.css";
import { createNativeWindPanel } from "./native-wind";
import { GroupPlayer } from "./audio";
import { blindLine, blindOrder, type Order } from "./blind";
import { LineChart, type Series } from "./chart";
import { notes as sharedNotes, saveNotes } from "./notes";
import { asset, el, fmt, json, onVisible } from "./util";

type Badge = "none" | "blocked" | "ungraded" | "listening" | "close" | "mixed" | "off";
type Value = number | (number | null)[] | null;
interface Metric {
  id: string; label: string; unit: string; gap_unit: string;
  real: Value; ours: Value; gap: number | null;
  ruler: { label: string; value: number }[]; ruler_max: number | null; ruler_min: number | null;
  verdict: "within" | "outside" | "no-ruler" | "no-gap"; ratio: number | null; note: string;
  // The strict ruler: real recordings of the SAME note only. null when the pair predates it ("not computed").
  ruler_strict: { label: string; value: number }[] | null; ruler_strict_max: number | null; ruler_strict_min: number | null;
  verdict_strict: "within" | "outside" | "no-ruler" | "no-gap" | null; ratio_strict: number | null;
  inside_loose: boolean | null; inside_strict: boolean | null;
}
interface Tally { inside: number; outside: number; no_ruler: number }
type RowVerdict = Metric["verdict"] | "loose-only";
interface FigureSet { real: string | null; ours: string | null; diff: string | null; f_lo_hz: number; f_hi_hz: number; span_db: number; diff_range_db: number | null; size?: [number, number] }
interface Side { level_db: (number | null)[]; cents: (number | null)[] }
interface Track { n: number; ref_hz: number | null; real: Side | null; ours: Side | null }
interface Pair {
  note: string; dynamic: string; duration_s: number;
  real: { source: { name: string; credit: string; licence: string; licence_url: string | null; url: string | null; external: boolean }; reference_id: string; file: string; segment_s: [number, number] | null; url: string | null; licence: string; mic: string | null; pattern: string | null; distance_m: number | null; distance_status: string | null; room: string | null; format: string | null } | null;
  no_reference?: string | null; // why a render has no recording to compare with; real is null then
  ours: { wav: string; git_sha: string; observer: string; sample_rate: number; recipe: string } | null;
  audio: { real: string | null; ours: string | null }; level_match: string;
  alignment: { onset_real_s: number; onset_ours_s: number | null; highpass_hz: number; level_offset_db: number | null } | null;
  figures: { spectrogram: FigureSet; scalogram: FigureSet };
  partials: { t_s: number[]; f0_hz: Record<string, number | null>; tracks: Track[] } | null;
  metrics: Metric[];
  verdicts: { measures: number; loose: Tally; strict: Tally | null } | null; // null: no render, or a listening-only note
}
// What real recordings from another library earn against the same rulers (P2's cross-library-baseline.json, summed and checked by the generator).
interface Share { inside: number; outside: number; no_ruler: number; pct: number | null }
interface Spread { median: number; min: number; max: number }
type Baseline = { shown: false; why: string } | {
  shown: true; source: { name: string; url: string; licence: string; licence_url: string | null; credit: string }; independence: "other_library" | "other_singer_same_library";
  compared_with: string; grader: string[]; method: string[]; caveats: string[]; generated_utc: string | null;
  recordings: { graded: number; refused: { slot: string; reason: string }[]; retuned: number; notes: string[]; nothing_outside: { loose: number; strict: number; both: number } };
  real: { loose: Share; strict: Share };
  render: { loose: Share | null; strict: Share | null; recordings: number; missing_notes: string[] };
  ruler_rows: { loose: Spread; strict: Spread };
  by_measure: { id: string; label: string; loose: { outside: number; ruled: number }; strict: { outside: number; ruled: number } }[];
};
interface Card {
  family: string; title: string; built_utc: string; baseline?: Baseline | null;
  reference: { instrument: string; files: number; by_dynamic: Record<string, number>; chains: string[]; licence: string; licence_url: string | null; credit: string | null } | null;
  status: { badge: Badge; label: string; line: string; summary: string; blocked: string | null; physics_present: string[]; physics_missing: string[]; beads: string[]; described: boolean };
  pairs: Pair[]; refused: { note: string; kind?: "grader" | "site"; problems: string[]; reason?: string; render_sha?: string | null; time?: string | null }[];
  withheld?: string[];
}
interface IndexFamily { id: string; title: string; reference: string | null; badge: Badge; label: string; line: string; counts: { notes_with_render: number; notes_reference_only: number; notes_no_reference?: number; notes_refused_by_grader: number; within: number; outside: number; unruled: number }; refused: number }
interface Index { built_utc: string; families: IndexFamily[] }

const data = (path: string): string => asset(`data/instruments/${path}`);
const MARK: Record<Badge, string> = { close: "✓", mixed: "◐", off: "✕", ungraded: "◌", listening: "♪", blocked: "■", none: "○" };
const PARTIAL_COLOURS = ["#d9a94f", "#63d2c0", "#8fb9d6", "#e86a3c", "#c79ad8", "#a6d46f", "#f1eadb", "#e0708f"];
const save = (_notes: Record<string, string>): void => saveNotes(); // the one shared store (notes.ts), written whole

const num = (v: number | null | undefined, unit = "", digits = 3): string => (v === null || v === undefined || !Number.isFinite(v) ? "—" : `${fmt(v, digits)}${unit ? ` ${unit}` : ""}`);
// Four significant digits: a real 146.8 Hz and a render at 147.2 Hz must not both read "147".
const show = (v: Value, unit: string): string => (Array.isArray(v) ? v.map((x) => (x === null ? "—" : fmt(x, 2))).join(" · ") + (unit ? ` ${unit}` : "") : num(v, unit, 4));
const pill = (badge: Badge, label: string): HTMLElement => el("span", { class: `ir-pill ir-${badge}` }, el("span", { "aria-hidden": "true" }, MARK[badge]), label);
const hz = (f: number): string => (f >= 1000 ? `${f / 1000}k` : String(f));

/** Frequency labels down one side of a log-frequency map, and time labels along its foot. */
function axes(fig: FigureSet, duration: number): HTMLElement {
  const box = el("div", { class: "ir-axis", "aria-hidden": "true" });
  for (const f of [50, 100, 200, 500, 1000, 2000, 5000, 10000]) {
    if (f <= fig.f_lo_hz * 1.05 || f >= fig.f_hi_hz) continue;
    const y = 1 - Math.log(f / fig.f_lo_hz) / Math.log(fig.f_hi_hz / fig.f_lo_hz);
    box.append(el("span", { class: "ir-f", style: `top: calc(${(y * 100).toFixed(2)}% - 7px)` }, hz(f)));
  }
  const step = duration <= 1.5 ? 0.25 : duration <= 3.5 ? 0.5 : 1;
  for (let t = 0; t < duration - step * 0.3; t += step) box.append(el("span", { class: "ir-t", style: `left: ${((t / duration) * 100).toFixed(2)}%` }, `${t % 1 ? t.toFixed(2).replace(/0$/, "") : t} s`));
  return box;
}

function figure(label: string, src: string | null, fig: FigureSet, duration: number, alt: string, absent: string): HTMLElement {
  const [w, h] = fig.size ?? [700, 260];
  if (!src) return el("figure", { class: "ir-fig empty" }, el("figcaption", {}, label), el("div", { class: "ir-img-box", style: `aspect-ratio: ${w} / ${h}` }, el("p", {}, absent)));
  const img = el("img", { class: "ir-img", src: data(src), alt, width: String(w), height: String(h), decoding: "async" });
  return el("figure", { class: "ir-fig" }, el("figcaption", {}, label), el("div", { class: "ir-img-box", style: `aspect-ratio: ${w} / ${h}` }, img, axes(fig, duration)));
}

/** One verdict for a row, from both rulers: a measure inside the neighbouring-note spread but outside the same-note one is its own case. */
const rowVerdict = (m: Metric): RowVerdict => (m.verdict === "within" && m.verdict_strict === "outside" ? "loose-only" : m.verdict);

/** The headline of a row, plainly: which spreads the gap is inside and which it is outside. */
function headline(m: Metric): string {
  const s = m.verdict_strict;
  if (m.verdict === "within") return s === "within" ? "inside both spreads" : s === "outside" ? "inside the neighbouring-note spread, OUTSIDE the same-note spread" : "inside the neighbouring-note spread";
  if (m.verdict === "outside") return s === "outside" ? "outside both spreads" : "outside the neighbouring-note spread";
  return "";
}

/** A measure's gap against the spread between real recordings: pale band = this and neighbouring notes, dark band = the same note only, mark = this render. */
function rulerBar(m: Metric): HTMLElement {
  const top = Math.max(m.gap ?? 0, m.ruler_max ?? 0) * 1.15 || 1;
  const pct = (v: number): string => `${Math.min(100, (v / top) * 100).toFixed(2)}%`;
  const bar = el("div", { class: "ir-ruler", role: "img", "aria-label": m.gap === null ? "No render to measure" : `${num(m.gap, m.gap_unit)} from the recording; ${headline(m) || m.verdict}` });
  const band = (lo: number, hi: number, cls: string): HTMLElement => el("span", { class: cls, style: `left: ${pct(lo)}; width: max(3px, calc(${pct(hi)} - ${pct(lo)}))` });
  if (m.ruler_max !== null && m.ruler_min !== null) bar.append(el("span", { class: "ir-zero" }), band(m.ruler_min, m.ruler_max, "ir-band"));
  if (m.ruler_strict_max !== null && m.ruler_strict_min !== null) bar.append(band(m.ruler_strict_min, m.ruler_strict_max, "ir-band strict"));
  if (m.gap !== null) bar.append(el("span", { class: `ir-mark ${rowVerdict(m)}`, style: `left: ${pct(m.gap)}` }));
  return bar;
}

const distances = (n: number): string => `${n} real-vs-real distance${n === 1 ? "" : "s"}`;

function metricsTable(p: Pair): HTMLElement {
  const rows = p.metrics.map((m) => {
    const list = (rs: { label: string; value: number }[]): string => rs.map((r) => `${r.label} ${num(r.value, "", 3)}`).join(" · ");
    const range = (lo: number | null, hi: number | null): string => (lo === hi ? num(hi, m.gap_unit, 3) : `${num(lo, "", 3)}–${num(hi, m.gap_unit, 3)}`);
    const wide = m.ruler_max === null ? "no ruler for this measure" : m.verdict === "outside" ? `outside: ${num(m.ratio, "", 3)}× the largest such distance (${num(m.ruler_max, m.gap_unit, 3)})` : `${m.verdict === "within" ? "inside: " : ""}real recordings of these notes differ by ${range(m.ruler_min, m.ruler_max)}`;
    const strict = m.verdict_strict === null ? "not computed for this note"
      : m.verdict_strict === "no-ruler" ? "no comparison: there is no other recording of this note"
      : m.verdict_strict === "outside" ? `outside: ${num(m.ratio_strict, "", 3)}× the largest such distance (${num(m.ruler_strict_max, m.gap_unit, 3)})`
      : `${m.verdict_strict === "within" ? "inside: " : ""}real recordings of this note differ by ${range(m.ruler_strict_min, m.ruler_strict_max)}`;
    return el("tr", { class: `ir-row ${rowVerdict(m)}` },
      el("td", { class: "vs-measure" }, el("div", {}, m.label), m.note ? el("div", { class: "small muted" }, m.note) : null,
        m.ruler.length ? el("details", { class: "small muted ir-rulers" }, el("summary", {}, `neighbouring-note ruler: ${distances(m.ruler.length)}`), el("div", {}, list(m.ruler))) : null,
        m.ruler_strict?.length ? el("details", { class: "small muted ir-rulers strict" }, el("summary", {}, `same-note ruler: ${distances(m.ruler_strict.length)}`), el("div", {}, list(m.ruler_strict))) : null),
      el("td", { class: "num", "data-label": "Real" }, show(m.real, m.unit)),
      el("td", { class: "num", "data-label": "FrankenSim" }, p.ours ? show(m.ours, m.unit) : "—"),
      el("td", { class: "ir-gap", "data-label": "Gap" },
        m.gap === null ? null : el("div", { class: "num" }, num(m.gap, m.gap_unit, 3), el("span", { class: `ir-verdict ${rowVerdict(m)}` }, ` ${headline(m)}`)),
        rulerBar(m),
        m.gap === null ? null : el("div", { class: "small muted ir-two" }, el("div", {}, el("b", {}, "Neighbouring notes included: "), wide), el("div", {}, el("b", {}, "Same note only: "), strict))));
  });
  return el("div", { class: "table-scroll", tabindex: "0", "aria-label": "Measures against the two real-vs-real rulers" },
    el("table", { class: "vs-table ir-table" }, el("thead", {}, el("tr", {}, ...["Measure", "Real", "FrankenSim", "Gap against how much real recordings differ"].map((h) => el("th", {}, h)))), el("tbody", {}, ...rows)));
}

/** The note's own counts, both rulers, in plain words. */
function verdictLine(p: Pair): HTMLElement | null {
  const v = p.verdicts;
  if (!v) return null;
  const s = v.strict;
  const wide = `Against real recordings of this and neighbouring notes: ${v.loose.inside} of ${v.loose.inside + v.loose.outside} measures inside${v.loose.no_ruler ? ` (${v.loose.no_ruler} with no ruler)` : ""}.`;
  const same = s ? `Against real recordings of the same note only: ${s.inside} of ${s.inside + s.outside} inside${s.no_ruler ? ` (${s.no_ruler} ${s.no_ruler === 1 ? "has" : "have"} no same-note comparison)` : ""}.` : "Against real recordings of the same note only: not computed for this note.";
  return el("div", { class: "ir-counts" },
    el("p", {}, el("b", {}, `This note, ${v.measures} measures. `), wide, " ", same),
    el("p", { class: "small muted" }, "How to read a row: the pale band is how far real recordings of this and neighbouring notes (other dynamics and takes too) differ from each other, the darker band is the same note only, and the mark is how far this render is from the recording. A row marked ‘OUTSIDE the same-note spread’ is inside the pale band but outside the dark one."));
}

/** A share as the owner reads it: the whole-number percent of the measures that had a ruler, or why there is none. */
const pctText = (s: Share | null): string => (s === null ? "not computed" : s.pct === null ? "no measure had a ruler" : `${s.pct}%`);
const ofText = (s: Share | null): string => (s === null ? "not computed" : `${s.inside} of ${s.inside + s.outside}`);
const plural = (n: number, one: string, many = `${one}s`): string => `${n} ${n === 1 ? one : many}`;

/** Beside the badge: what a real instrument from another library earns against the same rulers, next to this render's own share.
 *  A family with a graded badge and no baseline says so; the badge rule is not changed by any of this. */
function baselineBlock(c: Card): HTMLElement | null {
  const b = c.baseline ?? null;
  if (b === null) {
    if (!["close", "mixed", "off"].includes(c.status.badge)) return null;
    return el("div", { class: "ir-baseline none" }, el("p", {}, el("b", {}, "No cross-library baseline yet. "), "Real recordings of this instrument from another library have not been graded against these rulers, so there is no figure for what a real instrument earns here."));
  }
  if (!b.shown) return el("div", { class: "ir-baseline refused", role: "status" }, el("p", {}, el("b", {}, "Cross-library baseline not shown: "), b.why));
  const number = (s: Share | null): HTMLElement => el("b", {}, pctText(s));
  // Another singer from the same recording set is NOT another library: same room, same microphone. The data says which, the page says it plainly.
  const lead = b.independence === "other_library" ? "Real recordings from another library" : "Another singer from the same recording set (the same library, room and microphone: not another library, so a weaker check)";
  const render = b.render.loose === null ? el("span", {}, "no graded render of these notes") : el("span", {}, number(b.render.loose), " / ", number(b.render.strict));
  const covered = b.render.recordings < b.recordings.graded ? `, ${plural(b.render.recordings, "recording")} of ${b.recordings.graded} have a graded render (none for ${b.render.missing_notes.join(", ")})` : "";
  const r = b.recordings;
  const rows = b.ruler_rows;
  const worst = b.by_measure.filter((m) => m.loose.outside > 0).slice(0, 6);
  return el("div", { class: "ir-baseline" },
    el("p", { class: "ir-bl-head" }, `${lead}: `, number(b.real.loose), " of measures inside (loose) / ", number(b.real.strict), " (strict). This render: ", render, "."),
    el("p", { class: "small muted" }, "A real instrument does not land inside every ruler either, so a measure outside one is not by itself a sign that a render is wrong. Loose: against real recordings of this and neighbouring notes. Strict: the same note only."),
    el("p", { class: "small muted" }, `Real: ${ofText(b.real.loose)} measures loose, ${ofText(b.real.strict)} strict, over ${plural(r.graded, "recording")} of ${plural(r.notes.length, "note")} (${r.notes.join(", ")}). This render, the same notes: ${ofText(b.render.loose)} loose, ${ofText(b.render.strict)} strict${covered}. Source: ${b.source.name}, ${b.source.licence}, ${b.source.credit} (`, el("a", { href: b.source.url, rel: "noopener" }, b.source.url), ")."),
    el("details", { class: "small muted ir-bl-more" }, el("summary", {}, "How this was measured, and its limits"),
      ...b.method.map((t) => el("p", {}, t)),
      el("p", {}, `Compared with: ${b.compared_with}. ${r.nothing_outside.both} of ${r.graded} recordings had no measure outside either ruler${r.nothing_outside.both === 0 ? ": by the badge rule none of these real recordings would earn ‘Inside the real spread’" : ""}.`),
      r.refused.length ? el("p", {}, `Not counted, refused by the grader: ${r.refused.map((x) => `${x.slot} (${x.reason})`).join("; ")}.`) : null,
      r.retuned ? el("p", {}, `${plural(r.retuned, "recording")} of ${r.graded} ${r.retuned === 1 ? "was" : "were"} digitally retuned by the library's authors.`) : null,
      el("p", {}, `The rulers are built from a median of ${rows.loose.median} real-vs-real distances per measure (loose; ${rows.loose.min}–${rows.loose.max}) and ${rows.strict.median} (strict; ${rows.strict.min}–${rows.strict.max}); a measure is inside when its gap is no larger than the largest.`),
      worst.length ? el("p", {}, `Measures real recordings land outside (loose, strict): ${worst.map((m) => `${m.label} ${m.loose.outside}/${m.loose.ruled}, ${m.strict.outside}/${m.strict.ruled}`).join("; ")}.`) : null,
      el("ul", {}, ...b.caveats.map((t) => el("li", {}, t))),
      el("p", {}, `Graded with ${b.grader.join(", ") || "an unstated tool"}${b.generated_utc ? `; file written ${b.generated_utc.replace("T", " ").replace("Z", " UTC")}` : ""}.`)));
}

function physicsLists(c: Card): HTMLElement {
  if (!c.status.described) {
    return el("p", { class: "small muted" }, `No description data has been written for this family yet (data/instruments/descriptions/${c.family}/status.json is missing), so the physics that runs and the physics that is missing are not listed.`);
  }
  const list = (cls: string, title: string, items: string[], none: string): HTMLElement =>
    el("div", { class: `ir-phys ${cls}` }, el("h4", {}, title), items.length ? el("ul", {}, ...items.map((t) => el("li", {}, t))) : el("p", { class: "small muted" }, none));
  return el("div", { class: "grid two" },
    list("present", "Physics that runs", c.status.physics_present, "None listed."),
    list("missing", "Physics that is missing", c.status.physics_missing, "None listed."));
}

function referenceLine(c: Card): string {
  if (!c.reference) {
    if (c.pairs.length > 0 && c.pairs.every((p) => p.real === null)) return "There is no reference recording for this instrument, so these renders are not compared with anything and not graded: listen to them and look at their spectrogram and scalogram.";
    const s = c.pairs.find((p) => p.real)?.real?.source;
    if (s) return `Reference: ${s.name}. Credit: ${s.credit}. Licence: ${s.licence}`;
    return "The reference library has no recording of this instrument, so there is nothing to grade it against yet.";
  }
  const dyn = Object.entries(c.reference.by_dynamic).map(([k, v]) => `${v} ${k}`).join(", ");
  return `Reference: the University of Iowa Musical Instrument Samples, ${c.reference.instrument} (${c.reference.files} recordings: ${dyn}).${c.reference.credit ? ` Credit: ${c.reference.credit}.` : ""} The library's licence, quoted: “${c.reference.licence}”`;
}

/** The attribution that sits beside the clips: who recorded this note, and on what terms. */
function creditLine(p: Pair): HTMLElement {
  if (!p.real) return el("p", { class: "small muted ir-credit" }, "No recording was used for this note, so there is no recording credit.");
  const s = p.real.source;
  const quoted = s.licence.length > 30 ? `“${s.licence}”` : s.licence;
  return el("p", { class: "small muted ir-credit" }, `Recording credit: ${s.credit}. Licence: ${quoted}`, s.licence_url ? " (" : null, s.licence_url ? el("a", { href: s.licence_url, rel: "noopener" }, "terms") : null, s.licence_url ? ")" : null, /[.”]$/.test(quoted) ? null : ".");
}

function sourceLines(p: Pair): HTMLElement {
  const r = p.real;
  const level = p.level_match.trim().replace(/\.$/, ""); // a sentence from the data may or may not end in a full stop
  const levelSentence = `${level.charAt(0).toUpperCase()}${level.slice(1)}.`;
  if (!r) {
    return el("div", { class: "small muted ir-source" },
      p.ours ? el("p", {}, `Render: commit ${p.ours.git_sha.slice(0, 12)}, observed at ${p.ours.observer}, ${fmt(p.ours.sample_rate / 1000)} kHz; recipe ${p.ours.recipe}.`) : null,
      el("p", {}, levelSentence));
  }
  const mic = [r.mic, r.pattern === "unstated" ? null : r.pattern, r.distance_m !== null ? `at ${num(r.distance_m, "m")}` : r.distance_status === "absent-on-page" ? "distance not stated on the source page" : null, r.room].filter(Boolean).join(", ");
  const seg = r.segment_s ? `, ${num(r.segment_s[0], "")}–${num(r.segment_s[1], "s")} of the file` : "";
  const a = p.alignment;
  return el("div", { class: "small muted ir-source" },
    el("p", {}, "Recording: ", r.url ? el("a", { href: r.url, rel: "noopener" }, r.file) : r.file, `${seg}. ${mic}${r.format ? `; ${r.format}` : ""}.`),
    p.ours ? el("p", {}, `Render: commit ${p.ours.git_sha.slice(0, 12)}, observed at ${p.ours.observer}, ${fmt(p.ours.sample_rate / 1000)} kHz; recipe ${p.ours.recipe}.`) : el("p", {}, "Render: none for this note."),
    el("p", {}, `${levelSentence}${a ?` ${p.ours ? "Both signals are" : "The recording is"} high-passed at ${num(a.highpass_hz, "Hz")}; onsets ${num(a.onset_real_s, "s")} (recording)${a.onset_ours_s !== null ? ` and ${num(a.onset_ours_s, "s")} (render)` : ""} aligned.` : ""}`));
}

function buildCard(c: Card, notes: Record<string, string>): HTMLElement {
  const root = el("article", { class: "panel ir-card", "data-family": c.family, "data-blind": "false" });
  const status = el("div", { class: `vs-verdict ir-status ir-${c.status.badge}` },
    el("div", { class: "vs-verdict-head" }, el("span", { class: "vs-mark", "aria-hidden": "true" }, MARK[c.status.badge]), el("span", { class: "vs-word" }, c.status.label), el("h4", {}, "Status")),
    el("p", {}, c.status.line),
    baselineBlock(c),
    c.status.summary ? el("p", { class: "small muted" }, c.status.summary) : null);
  root.append(
    el("div", { class: "panel-head" }, el("h3", {}, c.title), el("div", { class: "tags" }, pill(c.status.badge, c.status.label), c.status.beads.length ? el("span", { class: "tag plain" }, c.status.beads.join(" · ")) : null)),
    status,
    el("p", { class: "small muted", style: "margin: 12px 0 0" }, referenceLine(c)));
  const byGrader = c.refused.filter((r) => r.kind === "grader");
  const bySite = c.refused.filter((r) => r.kind !== "grader");
  if (byGrader.length) {
    root.append(el("div", { class: "ir-refused small", role: "status" }, el("b", {}, `${byGrader.length} note${byGrader.length === 1 ? "" : "s"} refused by the grader, in no badge: `),
      ...byGrader.map((r) => el("div", {}, `${r.note}: refused: ${r.reason ?? r.problems.join("; ")}${r.render_sha ? ` (render ${r.render_sha.slice(0, 12)}${r.time ? `, ${r.time}` : ""})` : ""}`))));
  }
  if (bySite.length) {
    root.append(el("div", { class: "ir-refused small", role: "status" }, el("b", {}, `${bySite.length} analysis pair${bySite.length === 1 ? "" : "s"} refused, not shown: `),
      ...bySite.map((r) => el("div", {}, `${r.note}: ${r.problems.join("; ")}`))));
  }
  if (c.withheld?.length) {
    root.append(el("div", { class: "ir-refused small", role: "status" }, el("b", {}, `Showing ${c.pairs.length} of ${c.pairs.length + c.withheld.length} notes: `),
      `this review build leaves out ${c.withheld.join(", ")} to fit the page host's file limit. They are not shown here and are in no badge.`));
  }
  if (c.pairs.length === 0) {
    const kind = c.status.badge === "blocked" ? "blocked" : c.status.badge === "ungraded" ? "ungraded" : "none";
    root.append(el("div", { class: "ir-empty" }, el("h4", {}, kind === "blocked" ? "Blocked" : kind === "ungraded" ? "Rendered, nothing to compare" : "No render yet"),
      el("p", {}, kind === "blocked" ? c.status.blocked ?? "" : kind === "ungraded" ? c.status.line : "Nothing has been rendered for this family, so there is nothing to hear or measure. The comparison with a real recording appears here when there is."),
      c.reference && !c.withheld?.length ? el("p", { class: "small muted" }, "The real recordings are in the reference library; their spectrograms, scalograms and the real-vs-real ruler appear here once they have been analysed.") : null),
      el("h4", { class: "ir-h" }, "What the description says"), physicsLists(c));
    return root;
  }

  let pair = c.pairs[0]!;
  let view: "spectrogram" | "scalogram" = "spectrogram";
  let blind = false;
  let order: Order = ["real", "ours"]; // [which clip plays as X, which as Y]: seeded per note (blind.ts), so a note reopens in the same order
  const shuffles = new Map<string, number>(); // how many times each note was re-shuffled in this session
  let allPartials = false;
  let player: GroupPlayer | null = null;
  const listening = (): boolean => pair.real === null; // a render with no reference recording: only our sound and our pictures

  const seg = el("div", { class: "seg ir-notes", role: "group", "aria-label": "Choose a note" });
  const playA = el("button", { class: "btn sm", type: "button" });
  const playB = el("button", { class: "btn sm", type: "button" });
  const stop = el("button", { class: "btn ghost sm", type: "button" }, "Stop");
  const blindBox = el("label", { class: "toggle" }, el("input", { type: "checkbox" }), "Blind test: hide which is which");
  const reveal = el("div", { class: "small vs-reveal", "aria-live": "polite" });
  const listenNote = el("p", { class: "small muted" });
  const credit = el("div", {});
  const tabs = el("div", { class: "seg vs-tabs", role: "group", "aria-label": "Decomposition to show" });
  const figs = el("div", { class: "ir-figs" });
  const figLegend = el("p", { class: "small muted" });
  const levelCanvas = el("canvas", { "aria-label": "Partial levels over time, recording and render" });
  const centsCanvas = el("canvas", { "aria-label": "Partial frequencies over time, in cents from the recording's median" });
  const partialsBox = el("div", { class: "ir-partials" });
  const partialLegend = el("div", { class: "legend" });
  const morePartials = el("button", { class: "btn ghost sm", type: "button" });
  const tableBox = el("div", {});
  const source = el("div", {});
  const lookHead = el("h4", { class: "ir-h" });
  const measuredHead = el("h4", { class: "ir-h" });
  const blindHidden = el("p", { class: "small muted ir-blind-note" }, "The analysis is hidden during the blind test so that it cannot give the answer away. Switch the blind test off to see it.");
  const feedback = el("textarea", { class: "ir-feedback", rows: "3", placeholder: "What do you hear? What is wrong, and where? Kept in this browser only; use “Copy all my notes” to send them.", "aria-label": "Your notes on this note" });
  const copied = el("span", { class: "small muted", "aria-live": "polite" });
  const copyAll = el("button", { class: "btn sm", type: "button" }, "Copy all my notes");

  // LineChart hangs its tooltip on the canvas's parent, so the wrappers exist before the charts.
  const levelWrap = el("div", { class: "chart" }, levelCanvas);
  const centsWrap = el("div", { class: "chart" }, centsCanvas);
  const levelChart = new LineChart(levelCanvas,{ xLabel: "time [s]", yLabel: "level [dB]", aspect: 0.5, xFormat: (v) => v.toFixed(1), yFormat: (v) => v.toFixed(0) });
  const centsChart = new LineChart(centsCanvas, { xLabel: "time [s]", yLabel: "frequency [cents vs real median]", aspect: 0.5, xFormat: (v) => v.toFixed(1), yFormat: (v) => v.toFixed(0), zeroLine: true });

  const which = (slot: 0 | 1): "real" | "ours" => (listening() ? "ours" : blind ? order[slot]! : slot === 0 ? "real" : "ours");
  const labelButtons = (): void => {
    const playing = player?.state === "playing" ? player.current : null;
    const name = (slot: 0 | 1): string => listening() ? "FrankenSim" : blind ? `Sound ${slot === 0 ? "X" : "Y"}` : slot === 0 ? `Real ${c.title.toLowerCase()}` : "FrankenSim";
    playA.textContent = `${playing === which(0) ? "■" : "▶"} ${name(0)}`;
    playB.textContent = `${playing === which(1) ? "■" : "▶"} ${name(1)}`;
    playB.hidden = pair.ours === null || listening();
    blindBox.hidden = pair.ours === null || listening();
  };
  const makePlayer = (): void => {
    player?.pause();
    const urls: Record<string, string> = {};
    if (pair.audio.real) urls["real"] = data(pair.audio.real);
    if (pair.audio.ours) urls["ours"] = data(pair.audio.ours);
    player = new GroupPlayer(urls, listening() ? "ours" : "real");
    player.onState = labelButtons;
  };
  const playSlot = (slot: 0 | 1): void => {
    if (!player) return;
    const id = which(slot);
    if (player.state === "playing" && player.current === id) { player.pause(); return; }
    player.select(id);
    if (player.state !== "playing") void player.play(0).catch(() => undefined);
  };
  playA.addEventListener("click", () => playSlot(0));
  playB.addEventListener("click", () => playSlot(1));
  stop.addEventListener("click", () => player?.pause());

  // Blind A/B: X and Y play the level-matched recording and render in an order seeded by the note; he picks the real one, the
  // page reveals it and adds "blind: picked X, real was Y" to his notes for the note. It records his guess and scores nothing.
  const askBlind = (): void => {
    order = blindOrder(key(), shuffles.get(key()) ?? 0);
    reveal.replaceChildren();
    if (!blind || !pair.ours || listening()) return;
    const ask = el("div", { class: "btn-row" }, el("span", {}, `Which one is the real ${c.title.toLowerCase()}?`));
    const guesses: HTMLButtonElement[] = [];
    for (const slot of [0, 1] as const) {
      const b = el("button", { class: "btn sm", type: "button" }, `Sound ${slot === 0 ? "X" : "Y"}`);
      b.addEventListener("click", () => {
        for (const g of guesses) g.disabled = true;
        const line = blindLine(order, slot);
        notes[key()] = `${(notes[key()] ?? "").trimEnd()}${(notes[key()] ?? "").trim() ? "\n" : ""}${line}`;
        save(notes);
        feedback.value = notes[key()]!;
        const is = (clip: "real" | "ours"): string => (clip === "real" ? "the recording" : "FrankenSim");
        reveal.append(el("p", { style: "margin-top: 8px" }, `You picked ${slot === 0 ? "X" : "Y"}. X was ${is(order[0])}, Y was ${is(order[1])}. Added to your notes for this note: “${line}”.`));
      });
      guesses.push(b);
      ask.append(b);
    }
    const again = el("button", { class: "btn ghost sm", type: "button" }, "Shuffle again");
    again.addEventListener("click", () => { shuffles.set(key(), (shuffles.get(key()) ?? 0) + 1); player?.pause(); askBlind(); labelButtons(); });
    ask.append(again);
    reveal.append(ask);
  };
  blindBox.querySelector("input")!.addEventListener("change", (e) => {
    blind = (e.target as HTMLInputElement).checked;
    root.dataset["blind"] = String(blind);
    player?.pause();
    askBlind();
    labelButtons();
  });

  const drawFigures = (): void => {
    const fig = pair.figures[view];
    const d = pair.duration_s;
    const kind = view === "spectrogram" ? "Spectrogram" : "Wavelet scalogram";
    const alt = (who: string): string => `${kind} of ${who}, ${pair.note} ${pair.dynamic}, ${fmt(fig.f_lo_hz)} Hz to ${fmt(fig.f_hi_hz)} Hz on a log axis`;
    figs.className = listening() ? "ir-figs one" : "ir-figs";
    figs.replaceChildren(
      ...(listening() ? [] : [figure(`Real ${c.title.toLowerCase()}`, fig.real, fig, d, alt("the real recording"), "")]),
      figure("FrankenSim", fig.ours, fig, d, alt("the render"), "No render yet for this note."),
      ...(listening() ? [] : [figure("Difference", fig.diff, fig, d, alt("the difference"), "Needs a render.")]));
    figLegend.textContent = listening()
      ? `${fmt(d)} s from the start of the render, ${fmt(fig.f_lo_hz)} Hz–${fmt(fig.f_hi_hz)} Hz (log), ${fmt(fig.span_db)} dB colour range. There is no recording to compare with, so there is no real map and no difference map.`
      : `${fmt(d)} s from the aligned onset, ${fmt(fig.f_lo_hz)} Hz–${fmt(fig.f_hi_hz)} Hz (log), ${fmt(fig.span_db)} dB colour range` +
        (pair.ours ? ` shared by the recording and the render. Difference: FrankenSim minus the recording, ±${fmt(fig.diff_range_db ?? 0)} dB (red: FrankenSim louder, blue: quieter, white: the same).` : ".");
    for (const [k, b] of Object.entries(tabButtons)) b.setAttribute("aria-pressed", String(k === view));
  };
  const tabButtons: Record<string, HTMLButtonElement> = {};
  for (const [k, label] of [["spectrogram", "Spectrogram (STFT)"], ["scalogram", "Wavelet scalogram"]] as const) {
    const b = el("button", { type: "button", "aria-pressed": String(k === view) }, label);
    b.addEventListener("click", () => { view = k; drawFigures(); });
    tabButtons[k] = b;
    tabs.append(b);
  }

  const drawPartials = (): void => {
    const pt = pair.partials;
    partialsBox.hidden = pt === null;
    if (!pt) return;
    const tracks = pt.tracks.slice(0, allPartials ? 8 : 4);
    const series = (key: keyof Side): Series[] => tracks.flatMap((t): Series[] => {
      const colour = PARTIAL_COLOURS[(t.n - 1) % PARTIAL_COLOURS.length]!;
      const one = (side: Side | null, label: string, dash?: number[]): Series[] => (side ? [{ label: `n=${t.n} ${label}`, color: colour, x: pt.t_s, y: side[key].map((v) => v ?? NaN), width: 1.6, ...(dash ? { dash } : {}) }] : []);
      return [...one(t.real, "recording"), ...one(t.ours, "FrankenSim", [5, 3])];
    });
    // A chart with nothing finite to draw is hidden: LineChart would leave the previous note's drawing on the canvas.
    const drawable = (s: Series[]): boolean => s.some((x) => Array.from(x.y).some((v) => Number.isFinite(v)));
    const levels = series("level_db");
    const cents = series("cents");
    levelWrap.hidden = !drawable(levels);
    centsWrap.hidden = !drawable(cents);
    if (drawable(levels)) levelChart.set(levels);
    if (drawable(cents)) centsChart.set(cents);
    partialLegend.replaceChildren(...tracks.map((t) => {
      const i = el("i");
      i.style.background = PARTIAL_COLOURS[(t.n - 1) % PARTIAL_COLOURS.length]!;
      return el("span", {}, i, `partial ${t.n}${t.ref_hz ? ` (${num(t.ref_hz, "Hz", 4)})` : ""}`);
    }), el("span", {}, tracks.some((t) => t.ours) ? "solid: recording · dashed: FrankenSim" : pair.ours ? "solid: recording; none of the render's partials matched these" : "solid: recording"));
    morePartials.textContent = allPartials ? "Show fewer partials" : "Show 8 partials";
    morePartials.hidden = pt.tracks.length <= 4;
  };
  morePartials.addEventListener("click", () => { allPartials = !allPartials; drawPartials(); });

  const key = (): string => `${c.family}/${pair.note}`;
  const render = (): void => {
    for (const b of seg.querySelectorAll("button")) b.setAttribute("aria-pressed", String(b.dataset["note"] === pair.note));
    credit.replaceChildren(creditLine(pair));
    listenNote.textContent = listening() ? `${pair.note}, ${pair.dynamic}. No reference recording: ${pair.no_reference ?? ""}`
      : pair.ours ? `${pair.note}, ${pair.dynamic}. Switching keeps the playhead, so the two are heard at the same instant of the note.` : `${pair.note}, ${pair.dynamic}: only the real recording exists for this note. No render yet.`;
    lookHead.textContent = listening() ? "Look: the render, decomposed" : "Look: the same note, decomposed";
    measuredHead.textContent = listening() ? "How this render was made" : "Measured against how much real recordings differ";
    drawFigures();
    drawPartials();
    tableBox.replaceChildren(...(listening() ? [] : [verdictLine(pair), metricsTable(pair)].filter((x): x is HTMLElement => x !== null)));
    source.replaceChildren(sourceLines(pair));
    feedback.value = notes[key()] ?? "";
    labelButtons();
  };
  for (const p of c.pairs) {
    const b = el("button", { type: "button", "data-note": p.note, "aria-pressed": String(p === pair) }, `${p.note} · ${p.dynamic}`, p.real === null ? el("small", {}, "no reference") : p.ours ? null : el("small", {}, "recording only"));
    b.addEventListener("click", () => {
      pair = p;
      if (listening() && blind) { blind = false; root.dataset["blind"] = "false"; (blindBox.querySelector("input") as HTMLInputElement).checked = false; }
      makePlayer(); askBlind(); render();
    });
    seg.append(b);
  }
  feedback.addEventListener("input", () => {
    const v = feedback.value.trim();
    if (v) notes[key()] = feedback.value; else delete notes[key()];
    save(notes);
  });
  copyAll.addEventListener("click", () => {
    const entries = Object.entries(notes).filter(([, v]) => v.trim());
    const text = entries.length ? entries.map(([k, v]) => `## ${k}\n${v.trim()}`).join("\n\n") : "";
    if (!text) { copied.textContent = "Nothing written yet."; return; }
    void navigator.clipboard.writeText(`# Instrument listening notes\n\n${text}\n`).then(() => { copied.textContent = `Copied ${entries.length} note${entries.length === 1 ? "" : "s"}.`; }, () => { copied.textContent = "The browser would not copy; select the text of each box instead."; });
  });

  root.append(
    el("h4", { class: "ir-h" }, "Listen"), seg, listenNote,
    el("div", { class: "btn-row" }, playA, playB, stop, blindBox), credit, reveal,
    blindHidden,
    el("div", { class: "ir-analysis stack" },
      el("div", { class: "stack" }, lookHead, tabs, figs, figLegend),
      partialsBox,
      el("div", { class: "stack" }, measuredHead, tableBox, source)),
    el("h4", { class: "ir-h" }, "What the description says"), physicsLists(c),
    el("h4", { class: "ir-h" }, "Your notes"), feedback, el("div", { class: "btn-row" }, copyAll, copied));
  partialsBox.append(el("h4", { class: "ir-h" }, "Partial trajectories"),
    el("p", { class: "small muted" }, "The strongest partials tracked through the note. Level is in dB on the same scale for both; frequency is in cents from each partial's median in the recording, so a mistuned render sits off zero and vibrato shows as motion."),
    el("div", { class: "grid two" }, levelWrap, centsWrap), partialLegend, el("div", { class: "btn-row" }, morePartials));
  makePlayer();
  render();
  return root;
}

export async function mountInstruments(root: HTMLElement): Promise<void> {
  const nativeWind = createNativeWindPanel();
  root.replaceChildren(nativeWind, el("p", { class: "small muted" }, "Loading instrument cards…"));
  let index: Index;
  try { index = await json<Index>(data("index.json")); } catch (e) {
    root.style.minHeight = "0";
    root.replaceChildren(nativeWind, el("p", { class: "panel small", role: "status" }, `The instrument cards are not in this build (${(e as Error).message}). They are generated by tools/build_instruments.py from the analysis of the renders; run npm run data:instruments and build again.`));
    return;
  }
  const rank = (f: IndexFamily): number => (["close", "mixed", "off", "ungraded"].includes(f.badge) ? 0 : f.badge === "listening" ? 1 : f.counts.notes_reference_only > 0 ? 2 : f.badge === "blocked" ? 3 : 4);
  const families = index.families.map((f, i) => ({ f, i })).sort((a, b) => rank(a.f) - rank(b.f) || a.i - b.i).map((x) => x.f);
  const rendered = families.filter((f) => f.counts.notes_with_render > 0).length;
  const graded = families.filter((f) => ["close", "mixed", "off"].includes(f.badge)).length;
  const stat = (k: string, v: string): HTMLElement => el("div", { class: "stat" }, el("dt", {}, k), el("dd", {}, v));
  const tiles = el("div", { class: "ir-tiles", role: "group", "aria-label": "Choose an instrument family" });
  const host = el("div", { class: "ir-cards" });
  const notes = sharedNotes;
  const cards = new Map<string, HTMLElement>();
  const pending = new Map<string, Promise<HTMLElement>>();
  let current = "";

  const select = async (id: string): Promise<void> => {
    current = id;
    for (const t of tiles.querySelectorAll("button")) t.setAttribute("aria-pressed", String(t.dataset["family"] === id));
    let p = pending.get(id);
    if (!p) {
      p = json<Card>(data(`${id}/card.json`)).then((c) => buildCard(c, notes));
      pending.set(id, p);
      p.catch(() => pending.delete(id));
    }
    try {
      const card = await p;
      cards.set(id, card);
      if (current !== id) return;
      host.replaceChildren(card);
    } catch (e) {
      if (current === id) host.replaceChildren(el("p", { class: "panel small", role: "status" }, `This card could not be loaded (${(e as Error).message}). Reload to try again.`));
    }
  };
  for (const f of families) {
    const sub = f.badge === "listening" ? `${f.counts.notes_with_render} note${f.counts.notes_with_render === 1 ? "" : "s"} to listen to, nothing to grade against` : f.counts.notes_with_render ?`${f.counts.notes_with_render} note${f.counts.notes_with_render === 1 ? "" : "s"} rendered` : f.counts.notes_refused_by_grader ? `${f.counts.notes_refused_by_grader} refused by the grader` : f.counts.notes_reference_only ? "recordings analysed" : f.reference ? "recordings only in library" : "no recordings in library";
    const b = el("button", { type: "button", class: `ir-tile ir-${f.badge}`, "data-family": f.id, "aria-pressed": "false" }, el("b", {}, f.title), el("span", { class: "ir-tile-badge" }, el("span", { "aria-hidden": "true" }, MARK[f.badge]), f.label), el("small", {}, sub));
    b.addEventListener("click", () => void select(f.id));
    tiles.append(b);
  }
  root.replaceChildren(
    nativeWind,
    el("dl", { class: "stats" }, stat("families", String(families.length)), stat("with a render", String(rendered)), stat("graded", String(graded)), stat("no render yet", String(families.length - rendered))),
    tiles, host,
    el("p", { class: "small muted", style: "margin-top: 18px" }, `Cards built ${index.built_utc.replace("T", " ").replace("Z", " UTC")}. The badge only reports where the measured gaps fall against the spread between real recordings; it is not a claim that anything sounds like the real instrument.`));
  root.style.minHeight = "0";
  const first = families.find((f) => f.counts.notes_with_render > 0) ?? families.find((f) => f.counts.notes_reference_only > 0) ?? families[0];
  if (first) void select(first.id);
}

const root = document.getElementById("instruments-root");
if (root) onVisible(root, () => void mountInstruments(root), "800px");
