// The #instruments section: one card per instrument family, built from public/data/instruments/ (tools/build_instruments.py).
// Every number on a card comes from that data; this file only lays it out. A family with no render says so.
import "./instruments.css";
import { GroupPlayer } from "./audio";
import { LineChart, type Series } from "./chart";
import { asset, el, fmt, json, onVisible } from "./util";

type Badge = "none" | "blocked" | "ungraded" | "close" | "mixed" | "off";
type Value = number | (number | null)[] | null;
interface Metric {
  id: string; label: string; unit: string; gap_unit: string;
  real: Value; ours: Value; gap: number | null;
  ruler: { label: string; value: number }[]; ruler_max: number | null; ruler_min: number | null;
  verdict: "within" | "outside" | "no-ruler" | "no-gap"; ratio: number | null; note: string;
}
interface FigureSet { real: string; ours: string | null; diff: string | null; f_lo_hz: number; f_hi_hz: number; span_db: number; diff_range_db: number; size?: [number, number] }
interface Side { level_db: (number | null)[]; cents: (number | null)[] }
interface Track { n: number; ref_hz: number | null; real: Side | null; ours: Side | null }
interface Pair {
  note: string; dynamic: string; duration_s: number;
  real: { reference_id: string; file: string; segment_s: [number, number] | null; url: string; licence: string; mic: string | null; pattern: string | null; distance_m: number | null; distance_status: string | null; room: string | null; format: string | null };
  ours: { wav: string; git_sha: string; observer: string; sample_rate: number; recipe: string } | null;
  audio: { real: string; ours: string | null }; level_match: string;
  alignment: { onset_real_s: number; onset_ours_s: number | null; highpass_hz: number; level_offset_db: number | null } | null;
  figures: { spectrogram: FigureSet; scalogram: FigureSet };
  partials: { t_s: number[]; f0_hz: Record<string, number | null>; tracks: Track[] } | null;
  metrics: Metric[];
}
interface Card {
  family: string; title: string; built_utc: string;
  reference: { instrument: string; files: number; by_dynamic: Record<string, number>; chains: string[]; licence: string } | null;
  status: { badge: Badge; label: string; line: string; summary: string; blocked: string | null; physics_present: string[]; physics_missing: string[]; beads: string[]; described: boolean };
  pairs: Pair[]; refused: { note: string; problems: string[] }[];
}
interface IndexFamily { id: string; title: string; reference: string | null; badge: Badge; label: string; line: string; counts: { notes_with_render: number; notes_reference_only: number; within: number; outside: number; unruled: number }; refused: number }
interface Index { built_utc: string; families: IndexFamily[] }

const data = (path: string): string => asset(`data/instruments/${path}`);
const MARK: Record<Badge, string> = { close: "✓", mixed: "◐", off: "✕", ungraded: "◌", blocked: "■", none: "○" };
const PARTIAL_COLOURS = ["#d9a94f", "#63d2c0", "#8fb9d6", "#e86a3c", "#c79ad8", "#a6d46f", "#f1eadb", "#e0708f"];
const NOTES_KEY = "struck:instruments:notes";

function load(): Record<string, string> {
  try { return JSON.parse(localStorage.getItem(NOTES_KEY) ?? "{}") as Record<string, string>; } catch { return {}; }
}
function save(notes: Record<string, string>): void {
  try { localStorage.setItem(NOTES_KEY, JSON.stringify(notes)); } catch { /* storage unavailable: notes stay on the page only */ }
}

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

/** A measure's gap against the spread between real recordings: shaded band = real-vs-real, mark = this render. */
function rulerBar(m: Metric): HTMLElement {
  const top = Math.max(m.gap ?? 0, m.ruler_max ?? 0) * 1.15 || 1;
  const pct = (v: number): string => `${Math.min(100, (v / top) * 100).toFixed(2)}%`;
  const bar = el("div", { class: "ir-ruler", role: "img", "aria-label": m.gap === null ? "No render to measure" : `${num(m.gap, m.gap_unit)} from the recording; ${m.verdict}` });
  if (m.ruler_max !== null && m.ruler_min !== null) {
    const band = el("span", { class: "ir-band", style: `left: ${pct(m.ruler_min)}; width: max(3px, calc(${pct(m.ruler_max)} - ${pct(m.ruler_min)}))` });
    bar.append(el("span", { class: "ir-zero" }), band);
  }
  if (m.gap !== null) bar.append(el("span", { class: `ir-mark ${m.verdict}`, style: `left: ${pct(m.gap)}` }));
  return bar;
}

function metricsTable(p: Pair): HTMLElement {
  const rows = p.metrics.map((m) => {
    const spread = m.ruler_max === null ? "no real-vs-real ruler for this measure" : `real recordings differ by ${num(m.ruler_min, "", 3)}–${num(m.ruler_max, m.gap_unit, 3)}`;
    const verdict = m.verdict === "within" ? "inside the real spread" : m.verdict === "outside" ? `outside it (${num(m.ratio, "", 2)}× the largest real-vs-real distance)` : "";
    const ruler = m.ruler.length ? m.ruler.map((r) => `${r.label} ${num(r.value, "", 3)}`).join(" · ") : "";
    return el("tr", { class: `ir-row ${m.verdict}` },
      el("td", { class: "vs-measure" }, el("div", {}, m.label), m.note ? el("div", { class: "small muted" }, m.note) : null,
        ruler ? el("details", { class: "small muted ir-rulers" }, el("summary", {}, `ruler: ${m.ruler.length} real-vs-real distance${m.ruler.length === 1 ? "" : "s"}`), el("div", {}, ruler)) : null),
      el("td", { class: "num", "data-label": "Real" }, show(m.real, m.unit)),
      el("td", { class: "num", "data-label": "FrankenSim" }, p.ours ? show(m.ours, m.unit) : "—"),
      el("td", { class: "ir-gap", "data-label": "Gap" },
        m.gap === null ? null : el("div", { class: "num" }, num(m.gap, m.gap_unit, 3), el("span", { class: `ir-verdict ${m.verdict}` }, ` ${verdict}`)),
        rulerBar(m),
        el("div", { class: "small muted" }, spread)));
  });
  return el("div", { class: "table-scroll", tabindex: "0", "aria-label": "Measures against the real-vs-real ruler" },
    el("table", { class: "vs-table ir-table" }, el("thead", {}, el("tr", {}, ...["Measure", "Real", "FrankenSim", "Gap vs real-vs-real spread"].map((h) => el("th", {}, h)))), el("tbody", {}, ...rows)));
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
  if (!c.reference) return "The reference library has no recording of this instrument, so there is nothing to grade it against yet.";
  const dyn = Object.entries(c.reference.by_dynamic).map(([k, v]) => `${v} ${k}`).join(", ");
  return `Reference: the University of Iowa Musical Instrument Samples, ${c.reference.instrument} (${c.reference.files} recordings: ${dyn}). The library's licence, quoted: “${c.reference.licence}”`;
}

function sourceLines(p: Pair): HTMLElement {
  const r = p.real;
  const mic = [r.mic, r.pattern === "unstated" ? null : r.pattern, r.distance_m !== null ? `at ${num(r.distance_m, "m")}` : r.distance_status === "absent-on-page" ? "distance not stated on the source page" : null, r.room].filter(Boolean).join(", ");
  const seg = r.segment_s ? `, ${num(r.segment_s[0], "")}–${num(r.segment_s[1], "s")} of the file` : "";
  const a = p.alignment;
  return el("div", { class: "small muted ir-source" },
    el("p", {}, "Recording: ", el("a", { href: r.url, rel: "noopener" }, r.file), `${seg}. ${mic}${r.format ? `; ${r.format}` : ""}.`),
    p.ours ? el("p", {}, `Render: commit ${p.ours.git_sha.slice(0, 12)}, observed at ${p.ours.observer}, ${fmt(p.ours.sample_rate / 1000)} kHz; recipe ${p.ours.recipe}.`) : el("p", {}, "Render: none for this note."),
    el("p", {}, `${p.level_match.charAt(0).toUpperCase()}${p.level_match.slice(1)}.${a ? ` ${p.ours ? "Both signals are" : "The recording is"} high-passed at ${num(a.highpass_hz, "Hz")}; onsets ${num(a.onset_real_s, "s")} (recording)${a.onset_ours_s !== null ? ` and ${num(a.onset_ours_s, "s")} (render)` : ""} aligned.` : ""}`));
}

function buildCard(c: Card, notes: Record<string, string>): HTMLElement {
  const root = el("article", { class: "panel ir-card", "data-family": c.family, "data-blind": "false" });
  const status = el("div", { class: `vs-verdict ir-status ir-${c.status.badge}` },
    el("div", { class: "vs-verdict-head" }, el("span", { class: "vs-mark", "aria-hidden": "true" }, MARK[c.status.badge]), el("span", { class: "vs-word" }, c.status.label), el("h4", {}, "Status")),
    el("p", {}, c.status.line),
    c.status.summary ? el("p", { class: "small muted" }, c.status.summary) : null);
  root.append(
    el("div", { class: "panel-head" }, el("h3", {}, c.title), el("div", { class: "tags" }, pill(c.status.badge, c.status.label), c.status.beads.length ? el("span", { class: "tag plain" }, c.status.beads.join(" · ")) : null)),
    status,
    el("p", { class: "small muted", style: "margin: 12px 0 0" }, referenceLine(c)));
  if (c.refused.length) {
    root.append(el("div", { class: "ir-refused small", role: "status" }, el("b", {}, `${c.refused.length} analysis pair${c.refused.length === 1 ? "" : "s"} refused, not shown: `),
      ...c.refused.map((r) => el("div", {}, `${r.note}: ${r.problems.join("; ")}`))));
  }
  if (c.pairs.length === 0) {
    root.append(el("div", { class: "ir-empty" }, el("h4", {}, c.status.badge === "blocked" ? "Blocked" : "No render yet"),
      el("p", {}, c.status.badge === "blocked" ? c.status.blocked ?? "" : "Nothing has been rendered for this family, so there is nothing to hear or measure. The comparison with a real recording appears here when there is."),
      c.reference ? el("p", { class: "small muted" }, "The real recordings are in the reference library; their spectrograms, scalograms and the real-vs-real ruler appear here once they have been analysed.") : null),
      el("h4", { class: "ir-h" }, "What the description says"), physicsLists(c));
    return root;
  }

  let pair = c.pairs[0]!;
  let view: "spectrogram" | "scalogram" = "spectrogram";
  let blind = false;
  let order: ("real" | "ours")[] = ["real", "ours"];
  let allPartials = false;
  let player: GroupPlayer | null = null;

  const seg = el("div", { class: "seg ir-notes", role: "group", "aria-label": "Choose a note" });
  const playA = el("button", { class: "btn sm", type: "button" });
  const playB = el("button", { class: "btn sm", type: "button" });
  const stop = el("button", { class: "btn ghost sm", type: "button" }, "Stop");
  const blindBox = el("label", { class: "toggle" }, el("input", { type: "checkbox" }), "Blind test: hide which is which");
  const reveal = el("div", { class: "small vs-reveal", "aria-live": "polite" });
  const listenNote = el("p", { class: "small muted" });
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
  const blindHidden = el("p", { class: "small muted ir-blind-note" }, "The analysis is hidden during the blind test so that it cannot give the answer away. Switch the blind test off to see it.");
  const feedback = el("textarea", { class: "ir-feedback", rows: "3", placeholder: "What do you hear? What is wrong, and where? Kept in this browser only; use “Copy all my notes” to send them.", "aria-label": "Your notes on this note" });
  const copied = el("span", { class: "small muted", "aria-live": "polite" });
  const copyAll = el("button", { class: "btn sm", type: "button" }, "Copy all my notes");

  // LineChart hangs its tooltip on the canvas's parent, so the wrappers exist before the charts.
  const levelWrap = el("div", { class: "chart" }, levelCanvas);
  const centsWrap = el("div", { class: "chart" }, centsCanvas);
  const levelChart = new LineChart(levelCanvas,{ xLabel: "time [s]", yLabel: "level [dB]", aspect: 0.5, xFormat: (v) => v.toFixed(1), yFormat: (v) => v.toFixed(0) });
  const centsChart = new LineChart(centsCanvas, { xLabel: "time [s]", yLabel: "frequency [cents vs real median]", aspect: 0.5, xFormat: (v) => v.toFixed(1), yFormat: (v) => v.toFixed(0), zeroLine: true });

  const which = (slot: 0 | 1): "real" | "ours" => (blind ? order[slot]! : slot === 0 ? "real" : "ours");
  const labelButtons = (): void => {
    const playing = player?.state === "playing" ? player.current : null;
    const name = (slot: 0 | 1): string => blind ? `Sound ${slot === 0 ? "X" : "Y"}` : slot === 0 ? `Real ${c.title.toLowerCase()}` : "FrankenSim";
    playA.textContent = `${playing === which(0) ? "■" : "▶"} ${name(0)}`;
    playB.textContent = `${playing === which(1) ? "■" : "▶"} ${name(1)}`;
    playB.hidden = pair.ours === null;
    blindBox.hidden = pair.ours === null;
  };
  const makePlayer = (): void => {
    player?.pause();
    const urls: Record<string, string> = { real: data(pair.audio.real) };
    if (pair.audio.ours) urls["ours"] = data(pair.audio.ours);
    player = new GroupPlayer(urls, "real");
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

  const askBlind = (): void => {
    order = Math.random() < 0.5 ? ["real", "ours"] : ["ours", "real"];
    reveal.replaceChildren();
    if (!blind || !pair.ours) return;
    const ask = el("div", { class: "btn-row" }, el("span", {}, `Which one is the real ${c.title.toLowerCase()}?`));
    const guesses: HTMLButtonElement[] = [];
    for (const slot of [0, 1] as const) {
      const b = el("button", { class: "btn sm", type: "button" }, `Sound ${slot === 0 ? "X" : "Y"}`);
      b.addEventListener("click", () => {
        for (const g of guesses) g.disabled = true;
        reveal.append(el("p", { style: "margin-top: 8px" }, which(slot) === "real" ? "Correct: that was the recording. " : "That one was FrankenSim's render. ",
          `X was ${which(0) === "real" ? "the recording" : "FrankenSim"}, Y was ${which(1) === "real" ? "the recording" : "FrankenSim"}. Your answer stays in this browser. Pick another note to try again.`));
      });
      guesses.push(b);
      ask.append(b);
    }
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
    figs.replaceChildren(
      figure(`Real ${c.title.toLowerCase()}`, fig.real, fig, d, alt("the real recording"), ""),
      figure("FrankenSim", fig.ours, fig, d, alt("the render"), "No render yet for this note."),
      figure("Difference", fig.diff, fig, d, alt("the difference"), "Needs a render."));
    figLegend.textContent = `${fmt(d)} s from the aligned onset, ${fmt(fig.f_lo_hz)} Hz–${fmt(fig.f_hi_hz)} Hz (log), ${fmt(fig.span_db)} dB colour range` +
      (pair.ours ? ` shared by the recording and the render. Difference: FrankenSim minus the recording, ±${fmt(fig.diff_range_db)} dB (red: FrankenSim louder, blue: quieter, white: the same).` : ".");
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
    listenNote.textContent = pair.ours ? `${pair.note}, ${pair.dynamic}. Switching keeps the playhead, so the two are heard at the same instant of the note.` : `${pair.note}, ${pair.dynamic}: only the real recording exists for this note. No render yet.`;
    drawFigures();
    drawPartials();
    tableBox.replaceChildren(metricsTable(pair));
    source.replaceChildren(sourceLines(pair));
    feedback.value = notes[key()] ?? "";
    labelButtons();
  };
  for (const p of c.pairs) {
    const b = el("button", { type: "button", "data-note": p.note, "aria-pressed": String(p === pair) }, `${p.note} · ${p.dynamic}`, p.ours ? null : el("small", {}, "recording only"));
    b.addEventListener("click", () => { pair = p; makePlayer(); askBlind(); render(); });
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
    el("div", { class: "btn-row" }, playA, playB, stop, blindBox), reveal,
    blindHidden,
    el("div", { class: "ir-analysis stack" },
      el("div", { class: "stack" }, el("h4", { class: "ir-h" }, "Look: the same note, decomposed"), tabs, figs, figLegend),
      partialsBox,
      el("div", { class: "stack" }, el("h4", { class: "ir-h" }, "Measured against how much real recordings differ"), tableBox, source)),
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
  let index: Index;
  try { index = await json<Index>(data("index.json")); } catch (e) {
    root.style.minHeight = "0";
    root.replaceChildren(el("p", { class: "panel small", role: "status" }, `The instrument cards are not in this build (${(e as Error).message}). They are generated by tools/build_instruments.py from the analysis of the renders; run npm run data:instruments and build again.`));
    return;
  }
  const rank = (f: IndexFamily): number => (["close", "mixed", "off", "ungraded"].includes(f.badge) ? 0 : f.counts.notes_reference_only > 0 ? 1 : f.badge === "blocked" ? 2 : 3);
  const families = index.families.map((f, i) => ({ f, i })).sort((a, b) => rank(a.f) - rank(b.f) || a.i - b.i).map((x) => x.f);
  const rendered = families.filter((f) => f.counts.notes_with_render > 0).length;
  const graded = families.filter((f) => ["close", "mixed", "off"].includes(f.badge)).length;
  const stat = (k: string, v: string): HTMLElement => el("div", { class: "stat" }, el("dt", {}, k), el("dd", {}, v));
  const tiles = el("div", { class: "ir-tiles", role: "group", "aria-label": "Choose an instrument family" });
  const host = el("div", { class: "ir-cards" });
  const notes = load();
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
    const sub = f.counts.notes_with_render ? `${f.counts.notes_with_render} note${f.counts.notes_with_render === 1 ? "" : "s"} rendered` : f.counts.notes_reference_only ? "recordings analysed" : f.reference ? "recordings only in library" : "no recordings in library";
    const b = el("button", { type: "button", class: `ir-tile ir-${f.badge}`, "data-family": f.id, "aria-pressed": "false" }, el("b", {}, f.title), el("span", { class: "ir-tile-badge" }, el("span", { "aria-hidden": "true" }, MARK[f.badge]), f.label), el("small", {}, sub));
    b.addEventListener("click", () => void select(f.id));
    tiles.append(b);
  }
  root.replaceChildren(
    el("dl", { class: "stats" }, stat("families", String(families.length)), stat("with a render", String(rendered)), stat("graded", String(graded)), stat("no render yet", String(families.length - rendered))),
    tiles, host,
    el("p", { class: "small muted", style: "margin-top: 18px" }, `Cards built ${index.built_utc.replace("T", " ").replace("Z", " UTC")}. The badge only reports where the measured gaps fall against the spread between real recordings; it is not a claim that anything sounds like the real instrument.`));
  root.style.minHeight = "0";
  void select((families.find((f) => f.counts.notes_with_render > 0) ?? families.find((f) => f.counts.notes_reference_only > 0) ?? families[0]!).id);
}

const root = document.getElementById("instruments-root");
if (root) onVisible(root, () => void mountInstruments(root), "800px");
