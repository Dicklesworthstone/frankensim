import { GroupPlayer } from "./audio";
import { dpr, el, fmt, hidpi, sciNode, asset } from "./util";

interface Descriptor { real?: number | number[] | null; ours?: number | number[] | null; note?: string; rms_gap_db?: number | null; n?: number[]; bands_hz?: number[]; partials_used?: number[] }
interface Pair {
  id: string;
  kind: "piano" | "cymbal";
  title: string;
  duration_s: number;
  real: { label: string; url: string; sha256_wav48: string; rights: string; recording: string };
  ours: { clip: string; label: string };
  audio: { real: string; ours: string; level_match: string };
  images: { real: string; ours: string; diff: string; f_lo: number; f_hi: number; span_db: number; diff_range_db: number };
  descriptors: Record<string, Descriptor>;
}
interface KitSummary { pairs: { ruler: { ours_db: number; real_other_dynamic_db: number } }[] }
interface Verisim {
  pairs: Pair[];
  piano_ruler: { real_pp_vs_real_mf_db: number; real_mf_vs_real_ff_db: number; real_bass_mic_vs_treble_mic_db?: number; note: string };
  piano_dynamics: { real_db: Record<string, number>; ours_db: Record<string, number>; real_high_band_db: Record<string, number>; ours_high_band_db: Record<string, number>; note: string };
  cymbal_ruler: { reference: string; rows: { label: string; rms_gap_db: number; real: boolean; grid_mm?: number }[]; note: string };
  sources: { iowa_piano: string; iowa_cymbals: string; rights: string };
}

const LABELS: Record<string, string> = {
  f1_hz: "First partial",
  inharmonicity_B: "Inharmonicity B",
  partial_envelope_db: "Partial levels vs partial 1 (RMS gap)",
  partial_decay_db_s: "Partial decay rates 1–6",
  high_band_db: "Brightness: 2–8 kHz vs 0.2–2 kHz",
  attack_rise_ms: "Attack rise, 10→90 %",
  octave_early_db: "Octave-band balance, first 250 ms (RMS gap)",
  octave_decay_db_s: "Octave-band decay rates",
  low_knock_early_db: "Low knock, first half second",
  low_knock_late_db: "Low component, 1.5–2.9 s",
};

function show(v: Descriptor[keyof Descriptor], unit: string, digits = 3): string {
  if (v === null || v === undefined) return "—";
  if (Array.isArray(v)) return v.map((x) => (typeof x === "number" ? fmt(x, 2) : "—")).join(" · ") + (unit ? ` ${unit}` : "");
  if (typeof v !== "number") return String(v);
  return `${fmt(v, digits)}${unit ? ` ${unit}` : ""}`;
}

function cell(v: Descriptor[keyof Descriptor], unit: string): Node {
  if (typeof v === "number" && v !== 0 && Math.abs(v) < 1e-2) return sciNode(v, 3);
  return document.createTextNode(show(v, unit));
}

const UNITS: Record<string, string> = { f1_hz: "Hz", high_band_db: "dB", attack_rise_ms: "ms", partial_decay_db_s: "dB/s", octave_decay_db_s: "dB/s", low_knock_early_db: "dB", low_knock_late_db: "dB" };

/** Bars against a ruler made only of real recordings. */
function rulerChart(canvas: HTMLCanvasElement, rows: { label: string; value: number; real: boolean }[], unit: string): void {
  const draw = (): void => {
    // Size once: resizing a canvas resets its 2D context, font included.
    const narrow = (canvas.clientWidth || canvas.parentElement?.clientWidth || 600) < 520;
    hidpi(canvas, 0, rows.length * (narrow ? 44 : 34) + 40);
    const g = canvas.getContext("2d")!;
    const d = dpr();
    const { width: W, height: H } = canvas;
    g.fillStyle = "#07080b";
    g.fillRect(0, 0, W, H);
    g.font = `${(narrow ? 10.5 : 11.5) * d}px "JetBrains Mono", monospace`;
    const labelW = Math.max(...rows.map((r) => g.measureText(r.label).width));
    const left = narrow ? 10 * d : Math.min(labelW + 22 * d, W * 0.5);
    const rowH = (narrow ? 44 : 34) * d;
    const max = Math.max(...rows.map((r) => r.value)) * 1.15;
    const span = W - left - 74 * d;
    const reals = rows.filter((r) => r.real).map((r) => r.value);
    if (reals.length) {
      // The spread of real-vs-real distances: the ruler itself.
      const x0 = left + (span * Math.min(...reals)) / max;
      const x1 = left + (span * Math.max(...reals)) / max;
      g.fillStyle = "rgba(143,185,214,0.10)";
      g.fillRect(x0, 0, Math.max(2 * d, x1 - x0), canvas.height);
      g.strokeStyle = "rgba(143,185,214,0.45)";
      g.setLineDash([4 * d, 4 * d]);
      g.beginPath(); g.moveTo(x1, 0); g.lineTo(x1, canvas.height); g.stroke();
      g.setLineDash([]);
    }
    g.textBaseline = "middle";
    rows.forEach((r, i) => {
      const y = 14 * d + i * rowH;
      const barY = narrow ? y + 18 * d : y + 3 * d;
      g.fillStyle = "rgba(201,194,179,0.9)";
      g.textAlign = narrow ? "left" : "right";
      g.fillText(r.label, narrow ? left : left - 10 * d, narrow ? y + 6 * d : y + 11 * d);
      const w = (span * r.value) / max;
      g.fillStyle = r.real ? "rgba(143,185,214,0.85)" : "#d9a94f";
      g.fillRect(left, barY, w, 14 * d);
      g.textAlign = "left";
      g.fillStyle = "#f1eadb";
      g.fillText(`${r.value.toFixed(1)} ${unit}`, left + w + 8 * d, barY + 7 * d);
    });
  };
  new ResizeObserver(draw).observe(canvas);
  draw();
}

export type VerdictKind = "close" | "mixed" | "off";

export function pianoToneVerdict(
  dull: number,
  softInside: boolean,
  matched: boolean,
): { kind: VerdictKind; text: string } {
  const absDull = Math.abs(dull);
  const brightInside = Number.isFinite(dull) && absDull <= 3;
  const kind: VerdictKind =
    softInside && brightInside ? "close" : softInside || brightInside ? "mixed" : "off";
  const softClause = softInside
    ? "pp and mf partial levels sit inside the real piano's own change of touch."
    : "pp and mf partial levels sit outside the real piano's own change of touch.";
  const matchedSuffix = matched ? ", even at matched loudness" : "";
  const ffClause = !Number.isFinite(dull)
    ? "ff brightness could not be measured."
    : brightInside
      ? `ff brightness is within ${absDull.toFixed(0)} dB of the recording${matched ? " at matched loudness" : ""}.`
      : `ff is ${absDull.toFixed(0)} dB ${dull >= 0 ? "too dull" : "too bright"}${matchedSuffix}.`;
  return {
    kind,
    text: `${softClause} ${ffClause} Several sounds the real piano lacks are listed below.`,
  };
}

/** Three-line verdict up front, computed like the findings: both the words and the grade follow the data. */
function verdicts(v: Verisim, kit: KitSummary | null): HTMLElement {
  const d = (id: string, k: string): Descriptor | undefined => v.pairs.find((x) => x.id === id)?.descriptors[k];
  const n = (id: string, k: string, side: "real" | "ours"): number => (d(id, k)?.[side] as number | null | undefined) ?? NaN;
  const relErr = (id: string, realRef: number): number => Math.abs(n(id, "inharmonicity_B", "ours") - realRef) / realRef;
  const realA3 = ["a3-pp", "a3-mf", "a3-ff"].map((id) => n(id, "inharmonicity_B", "real"));
  const meanA3 = realA3.reduce((a, b) => a + b, 0) / realA3.length;
  const bErr = Math.max(relErr("a3-mf", meanA3), relErr("c6-mf", n("c6-mf", "inharmonicity_B", "real")));
  const matched = v.pairs.some((p) => p.id === "a3-ff8");
  const ff = matched ? "a3-ff8" : "a3-ff";
  const dull = n(ff, "high_band_db", "real") - n(ff, "high_band_db", "ours");
  const gap = (id: string): number => (d(id, "partial_envelope_db")?.rms_gap_db as number | null | undefined) ?? NaN;
  const softInside = gap("a3-pp") <= v.piano_ruler.real_pp_vs_real_mf_db && gap("a3-mf") <= v.piano_ruler.real_pp_vs_real_mf_db;
  const grade = (k: VerdictKind): [VerdictKind, string, string] => (k === "close" ? [k, "✓", "Close"] : k === "mixed" ? [k, "◐", "Mixed"] : [k, "✕", "Not yet"]);
  const card = ([kind, mark, word]: [VerdictKind, string, string], title: string, text: string): HTMLElement =>
    el("div", { class: `vs-verdict ${kind}` },
      el("div", { class: "vs-verdict-head" }, el("span", { class: "vs-mark", "aria-hidden": "true" }, mark), el("span", { class: "vs-word" }, word), el("h4", {}, title)),
      el("p", {}, text));
  const tone = pianoToneVerdict(dull, softInside, matched);
  // The drum kit is graded piece by piece in its own section against a real kit (Salamander overheads).
  const pieces = kit?.pairs ?? [];
  const within = pieces.filter((p) => p.ruler.ours_db <= p.ruler.real_other_dynamic_db).length;
  const ratio = pieces.length ? pieces.map((p) => p.ruler.ours_db / Math.max(p.ruler.real_other_dynamic_db, 0.1)).sort((a, b) => a - b)[Math.floor(pieces.length / 2)]! : NaN;
  const kitKind: VerdictKind = pieces.length && within === pieces.length ? "close" : within > 0 || ratio <= 3 ? "mixed" : "off";
  return el("div", { class: "vs-verdicts" },
    card(grade(bErr <= 0.15 ? "close" : bErr <= 0.3 ? "mixed" : "off"), "Piano strings", `Inharmonicity within ${Math.ceil(bErr * 100)}% of a real Steinway at A3 and C6, measured from the partials of both.`),
    card(grade(tone.kind), "Piano tone", tone.text),
    ...(pieces.length ? [card(grade(kitKind), "Drum kit", `${within} of ${pieces.length} pieces sit as close to a real stroke as the same real drum played one dynamic away; the median piece is ${ratio.toFixed(1)}× that distance. Graded in the kit section above.`)] : []),
  );
}

/** Findings are derived from the data so the words cannot drift from the numbers. */
function findings(v: Verisim): HTMLElement[] {
  const p = (id: string): Pair | undefined => v.pairs.find((x) => x.id === id);
  const num = (id: string, k: string, side: "real" | "ours"): number => (p(id)?.descriptors[k]?.[side] as number | null | undefined) ?? NaN;
  const gap = (id: string): number => (p(id)?.descriptors["partial_envelope_db"]?.rms_gap_db as number | null | undefined) ?? NaN;
  const hz = (id: string, k: string, side: "real" | "ours"): number => ((p(id)?.descriptors[k] as { hz?: Record<string, number> } | undefined)?.hz?.[side]) ?? NaN;
  const sci = (x: number): string => {
    const [m, e] = x.toExponential(2).split("e");
    return `${m}×10${String(Number(e)).replace("-", "⁻").replace(/\d/g, (c) => "⁰¹²³⁴⁵⁶⁷⁸⁹"[Number(c)]!)}`;
  };
  const signed = (x: number): string => `${x >= 0 ? "+" : "−"}${Math.abs(x).toFixed(1)}`;
  const rB = ["a3-pp", "a3-mf", "a3-ff"].map((id) => num(id, "inharmonicity_B", "real"));
  const dyn = v.piano_dynamics;
  const brightGain = (side: "real_high_band_db" | "ours_high_band_db", hi: string): number => (dyn[side][hi] ?? NaN) - (dyn[side]["mf"] ?? NaN);
  const items: [string, string][] = [];
  items.push(["The string physics lands on a real Steinway.",
    `Inharmonicity fitted from the partials of A3 is ${sci(num("a3-mf", "inharmonicity_B", "ours"))} for FrankenSim and ${sci(Math.min(...rB))}–${sci(Math.max(...rB))} for the real piano across three dynamics. At C6 it is ${sci(num("c6-mf", "inharmonicity_B", "ours"))} against ${sci(num("c6-mf", "inharmonicity_B", "real"))}. Our scale describes a Model D and the recording is a Model B. Both are tuned near A = 440 Hz (A3 first partials ${num("a3-mf", "f1_hz", "ours").toFixed(2)} and ${num("a3-mf", "f1_hz", "real").toFixed(2)} Hz), so like is compared with like.`]);
  const ruler = v.piano_ruler;
  const within = gap("a3-pp") <= ruler.real_pp_vs_real_mf_db && gap("a3-mf") <= ruler.real_pp_vs_real_mf_db && gap("a3-ff") <= ruler.real_mf_vs_real_ff_db;
  items.push([within ? "Partial levels stay within one real dynamic step." : "Partial levels drift beyond one real dynamic step.",
    `The partial-level gap to the recording is ${gap("a3-pp").toFixed(1)} dB at pp, ${gap("a3-mf").toFixed(1)} dB at mf and ${gap("a3-ff").toFixed(1)} dB at ff. The real piano changes by ${ruler.real_pp_vs_real_mf_db.toFixed(1)} dB from pp to mf and ${ruler.real_mf_vs_real_ff_db.toFixed(1)} dB from mf to ff. At mf, 2–8 kHz brightness is within ${Math.abs(num("a3-mf", "high_band_db", "real") - num("a3-mf", "high_band_db", "ours")).toFixed(1)} dB. This is one coarse measure; the sounds below are what set the two apart by ear.`]);
  if (dyn.ours_db["ff8"] !== undefined) {
    const ff8Dull = num("a3-ff8", "high_band_db", "real") - num("a3-ff8", "high_band_db", "ours");
    const ff8Abs = Math.abs(ff8Dull);
    const ff8Word = ff8Dull >= 0 ? "dull" : "bright";
    const ff8Comp = ff8Dull >= 0 ? "duller" : "brighter";
    items.push([
      ff8Abs <= 3 ? "ff brightness lands within 3 dB at matched loudness." : `ff is too ${ff8Word}, even at matched loudness.`,
      `From mf to ff the real piano gains ${brightGain("real_high_band_db", "ff").toFixed(1)} dB of 2–8 kHz brightness and ${dyn.real_db["ff"]!.toFixed(1)} dB of level. FrankenSim at 4.5 m/s gains ${brightGain("ours_high_band_db", "ff").toFixed(1)} dB of brightness and ${dyn.ours_db["ff"]!.toFixed(1)} dB of level. At 8 m/s, the hardest stroke rendered, it reaches ${dyn.ours_db["ff8"]!.toFixed(1)} dB of level, close to the real ff, ${ff8Abs <= 3 ? `and within ${ff8Abs.toFixed(1)} dB of its brightness.` : `but is still ${ff8Abs.toFixed(1)} dB ${ff8Comp}. Stroke speed alone does not explain it.`}`,
    ]);
  }
  const kA = hz("a3-mf", "low_knock_early_db", "ours"), kC = hz("c6-mf", "low_knock_early_db", "ours");
  const kEarly = num("a3-mf", "low_knock_early_db", "ours") - num("a3-mf", "low_knock_early_db", "real");
  const kLateO = num("a3-mf", "low_knock_late_db", "ours"), kLateR = num("a3-mf", "low_knock_late_db", "real");
  items.push([kEarly > 3 ? "A low resonance rings too loudly at the start." : kEarly < -3 ? "The low band at the start is quieter than the real piano's." : "The low band at the start matches the real piano within 3 dB.",
    `In the first half second FrankenSim has a component near ${kA.toFixed(0)} Hz at ${num("a3-mf", "low_knock_early_db", "ours").toFixed(1)} dB relative to the fundamental. The real piano's strongest component in that band sits at ${num("a3-mf", "low_knock_early_db", "real").toFixed(1)} dB. ${Math.abs(kA - kC) < 2 ? `It appears at the same frequency when C6 is struck (${kC.toFixed(1)} Hz), so it belongs to the modelled instrument, not to the A3 string. ` : ""}${Math.abs(kLateO - kLateR) <= 3 ? "Later the two agree" : "Later they still differ"} (${kLateO.toFixed(1)} dB and ${kLateR.toFixed(1)} dB).`]);
  const dr = p("a3-mf")?.descriptors["partial_decay_db_s"];
  const at = (arr: unknown, k: number): number => (Array.isArray(arr) && typeof arr[k - 1] === "number" ? (arr[k - 1] as number) : NaN);
  if (dr) {
    // A slower fall (less negative dB/s) than the recording means the partial rings too long.
    const way = (k: number): string => (at(dr.ours, k) > at(dr.real, k) ? "rings too long" : "dies away too fast");
    items.push([way(4) !== way(2) ? "Some partials ring too long, others too short." : `Partials 2 and 4 both ${way(4).replace("rings", "ring").replace("dies", "die")}.`,
      `At mf, partial 4 decays at ${at(dr.ours, 4).toFixed(1)} dB/s in FrankenSim against ${at(dr.real, 4).toFixed(1)} dB/s on the real piano: it ${way(4)}. Partial 2 decays at ${at(dr.ours, 2).toFixed(1)} dB/s against ${at(dr.real, 2).toFixed(1)} dB/s: it ${way(2)}. Decay is set by losses at the bridge and in the string, and by how the three strings of a unison beat against each other, so this points at the coupling, the damping and the unison tuning.`]);
  }
  const riseO = num("a3-mf", "attack_rise_ms", "ours"), riseR = num("a3-mf", "attack_rise_ms", "real");
  const atIowa = p("a3-mf")?.ours.label.includes("where the Iowa microphones sat") ?? false;
  items.push([riseO < riseR - 1.5 ? "The attack is sharper." : riseO > riseR + 1.5 ? "The attack is slower." : "The attack rises like the recording's.",
    `FrankenSim rises in about ${riseO.toFixed(1)} ms, the recording in about ${riseR.toFixed(1)} ms. ` + (atIowa
      ? "Both are heard 8 inches above the strings, where the Iowa microphones sat; the recording's room and cardioid microphones remain differences."
      : "The Iowa microphones sat 8 inches above the strings in a normal room; ours is a single point 1 m above the board with no room, so part of this gap may belong to the recording rather than the piano.")]);
  const sp = p("splash-mf")?.descriptors["octave_early_db"]?.rms_gap_db as number | null | undefined;
  const spSoft = p("splash-pp")?.descriptors["octave_early_db"]?.rms_gap_db as number | null | undefined;
  const crashes = v.cymbal_ruler.rows.filter((row) => row.real && row.label.includes("crash")).map((row) => row.rms_gap_db);
  // The same stroke on other grids shows how far the graded render's figure has converged.
  const pairGrid = /on a ([\d.]+) mm grid/.exec(p("splash-mf")?.ours.label ?? "")?.[1];
  const grids = v.cymbal_ruler.rows.filter((row) => !row.real && typeof row.grid_mm === "number" && String(row.grid_mm) !== pairGrid).sort((a, b) => b.grid_mm! - a.grid_mm!);
  if (typeof sp === "number" && p("splash-mf")?.ours.label.includes("fs-vkplate") && crashes.length) {
    const lo = Math.min(...crashes), hi = Math.max(...crashes);
    // The octave furthest below the recording, from the same early octave balance the gap is measured on.
    const oct = p("splash-mf")?.descriptors["octave_early_db"] as { bands_hz?: number[]; real?: number[]; ours?: number[] } | undefined;
    let weak = "";
    if (oct?.bands_hz && oct.real && oct.ours) {
      const { real, ours } = oct;
      const short = oct.bands_hz.map((hz, i) => ({ hz, d: (ours[i] ?? NaN) - (real[i] ?? NaN) })).filter((b) => b.hz >= 1000 && Number.isFinite(b.d)).sort((a, b) => a.d - b.d)[0];
      if (short && short.d < -3) weak = ` Its ${short.hz / 1000} kHz octave sits ${(-short.d).toFixed(1)} dB below the recording's.`;
    }
    items.push([sp <= lo ? "The splash's octave balance is closer to the real splash than real crashes are." : sp <= hi ? "The splash's octave balance sits among the real crashes." : "The splash's octave balance sits beyond the real crashes.",
      `Struck at 2 m/s${pairGrid ? ` on a ${pairGrid} mm grid` : ""}, the 8-inch splash is ${sp.toFixed(1)} dB from the real splash in octave-band balance over the first 250 ms${typeof spSoft === "number" ? ` (${spSoft.toFixed(1)} dB at 0.6 m/s against the real pp)` : ""}.${grids.length ? ` The same stroke on other grids: ${grids.map((g) => `${g.rms_gap_db.toFixed(1)} dB (${g.grid_mm} mm)`).join(", ")}.` : ""} Three real crashes, struck the same way, sit ${lo.toFixed(1)}–${hi.toFixed(1)} dB from it. Its 2–8 kHz brightness is ${num("splash-mf", "high_band_db", "ours").toFixed(1)} dB against the recording's ${num("splash-mf", "high_band_db", "real").toFixed(1)} dB.${weak}`]);
  }
  // A missing number must not reach the page as "NaN dB".
  return items.filter(([, t]) => !t.includes("NaN")).map(([h, t]) => el("li", {}, el("b", {}, h), " ", t));
}

export async function mountVerisim(root: HTMLElement): Promise<void> {
  // The diagnosis has its own data file: it must not vanish when the scorecard fails to load.
  void mountDiagnosis();
  let data: Verisim;
  try {
    const res = await fetch(asset("data/verisim.json"));
    if (!res.ok) throw new Error(`HTTP ${res.status}`);
    data = (await res.json()) as Verisim;
    // A cymbal pair is shown only when it comes from the kit's own solver; the retired modal splash stays off the page.
    data.pairs = data.pairs.filter((p) => p.kind === "piano" || p.ours.label.includes("fs-vkplate"));
  } catch (e) {
    root.style.minHeight = "0";
    root.replaceChildren(el("p", { class: "panel small", role: "status" }, `The comparison data could not be loaded (${(e as Error).message}). Reload to try again.`));
    return;
  }
  const pairs = data.pairs;
  let current = pairs[0]!;
  let blind = false;
  let blindOrder: ("real" | "ours")[] = ["real", "ours"];

  const seg = el("div", { class: "seg vs-pairs", role: "group", "aria-label": "Choose a comparison" });
  const playA = el("button", { class: "btn sm", type: "button" });
  const playB = el("button", { class: "btn sm", type: "button" });
  const stop = el("button", { class: "btn ghost sm", type: "button" }, "Stop");
  const blindBox = el("label", { class: "toggle" }, el("input", { type: "checkbox" }), "Blind test: hide which is which");
  const reveal = el("div", { class: "small vs-reveal", "aria-live": "polite" });
  const subtitle = el("p", { class: "small muted", style: "margin-top: 8px" });
  const tabs = el("div", { class: "seg vs-tabs", role: "group", "aria-label": "Spectrogram to show" });
  const img = el("img", { alt: "", class: "vs-img", width: "700", height: "260", decoding: "async" });
  const axis = el("div", { class: "vs-axis", "aria-hidden": "true" });
  const imgBox = el("div", { class: "clip-view", style: "position:relative" }, img, axis);
  const legend = el("p", { class: "small muted" });
  const table = el("table", { class: "vs-table" });
  const tableBox = el("div", { class: "table-scroll", tabindex: "0", "aria-label": "Scorecard" }, table);
  const source = el("p", { class: "small muted" });

  let view: "real" | "ours" | "diff" = "real";
  const viewBtns: Record<string, HTMLButtonElement> = {};
  for (const [k, label] of [["real", "Real recording"], ["ours", "FrankenSim"], ["diff", "Difference"]] as const) {
    const b = el("button", { type: "button", "aria-pressed": String(k === view) }, label);
    b.addEventListener("click", () => { view = k; render(); });
    viewBtns[k] = b;
    tabs.append(b);
  }
  for (const f of [100, 250, 500, 1000, 2000, 4000, 8000]) {
    const y = 1 - Math.log(f / 40) / Math.log(16000 / 40);
    axis.append(el("span", { style: `position:absolute; right:6px; top:calc(${(y * 100).toFixed(2)}% - 7px); font: 500 10.5px var(--font-mono); color: rgba(241,234,219,0.8); text-shadow: 0 1px 2px #000` }, f >= 1000 ? `${f / 1000}k` : String(f)));
  }

  let player: GroupPlayer | null = null;
  const makePlayer = (): void => {
    player?.pause();
    player = new GroupPlayer({ real: asset(current.audio.real), ours: asset(current.audio.ours) }, "real");
    player.onState = () => labelButtons();
  };
  const which = (slot: 0 | 1): "real" | "ours" => (blind ? blindOrder[slot]! : slot === 0 ? "real" : "ours");
  const labelButtons = (): void => {
    const playing = player?.state === "playing" ? player.current : null;
    const name = (slot: 0 | 1): string => blind ? `Sound ${slot === 0 ? "X" : "Y"}` : slot === 0 ? (current.kind === "piano" ? "Real Steinway" : "Real splash") : "FrankenSim";
    playA.textContent = `${playing === which(0) ? "■" : "▶"} ${name(0)}`;
    playB.textContent = `${playing === which(1) ? "■" : "▶"} ${name(1)}`;
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
    blindOrder = Math.random() < 0.5 ? ["real", "ours"] : ["ours", "real"];
    reveal.replaceChildren();
    if (!blind) return;
    const ask = el("div", { class: "btn-row" }, el("span", {}, `Which one is the real ${current.kind === "piano" ? "piano" : "cymbal"}?`));
    const guesses: HTMLButtonElement[] = [];
    for (const slot of [0, 1] as const) {
      const b = el("button", { class: "btn sm", type: "button" }, `Sound ${slot === 0 ? "X" : "Y"}`);
      b.addEventListener("click", () => {
        for (const g of guesses) g.disabled = true;
        const right = which(slot) === "real";
        reveal.append(el("p", { style: "margin-top: 8px" }, right ? "Correct: that was the recording. " : "That one was FrankenSim's render. ",
          `X was ${which(0) === "real" ? "the recording" : "FrankenSim"}, Y was ${which(1) === "real" ? "the recording" : "FrankenSim"}. Your answer stays in this browser; nothing is collected. Pick another pair to try again.`));
      });
      guesses.push(b);
      ask.append(b);
    }
    reveal.append(ask);
  };
  blindBox.querySelector("input")!.addEventListener("change", (e) => {
    blind = (e.target as HTMLInputElement).checked;
    player?.pause();
    if (blind) view = "diff";
    askBlind();
    render();
  });

  const render = (): void => {
    // In a blind test the labelled spectrograms would give the answer away.
    if (blind && view !== "diff") view = "diff";
    for (const [k, b] of Object.entries(viewBtns)) b.setAttribute("aria-pressed", String(k === view));
    const sub = current.title.split(": ").slice(1).join(": ");
    subtitle.textContent = sub.charAt(0).toUpperCase() + sub.slice(1);
    viewBtns["real"]!.disabled = blind;
    viewBtns["ours"]!.disabled = blind;
    img.src = asset(current.images[view]);
    img.alt = `${view === "diff" ? "Difference" : view === "real" ? "Real recording" : "FrankenSim"} spectrogram, ${current.title}, 40 Hz to 16 kHz on a log axis`;
    legend.textContent = view === "diff"
      ? `Difference in dB, FrankenSim minus the recording, after matching early loudness: orange = FrankenSim louder, blue = FrankenSim quieter (±${current.images.diff_range_db} dB). Near-black = within a few dB of each other, or both below the shared floor.`
      : `${current.duration_s} s from the onset, 40 Hz–16 kHz (log), ${current.images.span_db} dB colour range shared by both sounds.`;
    const rows: HTMLElement[] = [el("tr", {}, el("th", {}, "Measure"), el("th", {}, current.kind === "piano" ? "Real Steinway" : "Real splash"), el("th", {}, "FrankenSim"))];
    for (const [k, dsc] of Object.entries(current.descriptors)) {
      const unit = UNITS[k] ?? "";
      const gap = dsc.rms_gap_db;
      const measure = el("td", { class: "vs-measure" }, el("div", {}, LABELS[k] ?? k), el("div", { class: "small muted", style: "margin-top:3px" }, dsc.note ?? ""));
      if (gap !== undefined && gap !== null) {
        rows.push(el("tr", {}, measure, el("td", { class: "num", colspan: "2", "data-label": "Real ↔ FrankenSim" }, `${gap.toFixed(1)} dB apart`, el("span", { class: "vs-wide" }, " (real ↔ FrankenSim)"))));
        continue;
      }
      const cls = Array.isArray(dsc.real) ? "num small" : "num";
      const realName = current.kind === "piano" ? "Real Steinway" : "Real splash";
      rows.push(el("tr", {}, measure, el("td", { class: cls, "data-label": realName }, cell(dsc.real, unit)), el("td", { class: cls, "data-label": "FrankenSim" }, cell(dsc.ours, unit))));
    }
    table.replaceChildren(el("thead", {}, rows[0]!), el("tbody", {}, ...rows.slice(1)));
    source.replaceChildren(
      el("span", {}, `Recording: ${current.real.label}, `), el("a", { href: current.real.url, rel: "noopener" }, "source file"),
      el("span", {}, `. ${current.real.recording}. ${current.audio.level_match.charAt(0).toUpperCase()}${current.audio.level_match.slice(1)}. Both signals high-passed at 40–70 Hz to remove room rumble. Rights: ${current.real.rights}.`),
    );
    labelButtons();
  };

  for (const p of pairs) {
    const b = el("button", { type: "button", "aria-pressed": String(p === current) }, p.title.split(":")[0]!.replace(" · ", " "));
    b.addEventListener("click", () => {
      current = p;
      for (const x of seg.querySelectorAll("button")) x.setAttribute("aria-pressed", String(x === b));
      makePlayer();
      askBlind();
      render();
    });
    seg.append(b);
  }

  const pianoRuler = el("canvas", { "aria-label": "Partial-envelope gaps: FrankenSim versus the real piano, and the real piano versus itself across dynamics" });
  const gap = (id: string): number => (pairs.find((p) => p.id === id)?.descriptors["partial_envelope_db"]?.rms_gap_db ?? NaN) as number;

  const kit = await fetch(asset("data/kit.json")).then((r) => (r.ok ? (r.json() as Promise<KitSummary>) : null)).catch(() => null);
  root.append(
    verdicts(data, kit),
    el("div", { class: "panel" },
      el("div", { class: "panel-head" }, el("h3", {}, "Listen and look, side by side"), el("div", { class: "tags" }, el("span", { class: "tag offline" }, "Offline reference"), el("span", { class: "tag plain" }, "Real recording"))),
      seg, subtitle,
      el("div", { class: "stack", style: "margin-top: 16px" },
        el("div", { class: "stack" },
          el("div", { class: "btn-row" }, playA, playB, stop, blindBox),
          reveal, tabs, imgBox, legend),
        el("div", { class: "stack" }, el("h4", {}, "Scorecard"), tableBox, source)),
    ),
    el("div", { style: "margin-top: 22px" },
      el("figure", { class: "panel" },
        el("h3", {}, "Piano: a ruler made of the real piano"),
        el("p", { class: "muted small", style: "margin: 8px 0 14px" }, "How far apart are the partial-level envelopes (partials 2–16 relative to partial 1)? Blue bars are the real Steinway against itself: one dynamic apart (what a change of touch alone does), and the same note heard by its two microphones (what the listening position alone does). Brass bars are FrankenSim against the real piano at the same dynamic. A brass bar left of the dashed line differs from the recording by less than the real piano differs from itself between mf and ff."),
        el("div", { class: "chart" }, pianoRuler)),
    ),
  );
  rulerChart(pianoRuler, [
    { label: "real pp vs real mf", value: data.piano_ruler.real_pp_vs_real_mf_db, real: true },
    { label: "real mf vs real ff", value: data.piano_ruler.real_mf_vs_real_ff_db, real: true },
    { label: "real mf, bass mic vs treble mic", value: data.piano_ruler.real_bass_mic_vs_treble_mic_db ?? NaN, real: true },
    { label: "FrankenSim vs real, pp", value: gap("a3-pp"), real: false },
    { label: "FrankenSim vs real, mf", value: gap("a3-mf"), real: false },
    { label: "FrankenSim vs real, ff", value: gap("a3-ff"), real: false },
    { label: "FrankenSim vs real, ff at 8 m/s", value: gap("a3-ff8"), real: false },
    { label: "FrankenSim vs real, C6 mf", value: gap("c6-mf"), real: false },
  ].filter((r) => Number.isFinite(r.value)), "dB");
  const host = document.getElementById("verisim-findings");
  if (host) host.replaceChildren(...findings(data));
  makePlayer();
  render();
}

interface DiagRow { id: string; symptom: string; evidence: string; suspect: string; test: string; status: "open" | "testing" | "confirmed" | "refuted" | "improved" | "fixed"; result?: string; fix?: string }

/** The differential diagnosis written by tools/diagnose.py: words and numbers come from the same run. */
async function mountDiagnosis(): Promise<void> {
  const root = document.getElementById("diagnosis-root");
  if (!root) return;
  let rows: DiagRow[];
  try {
    const res = await fetch(asset("data/diagnosis.json"));
    if (!res.ok) throw new Error(`HTTP ${res.status}`);
    rows = ((await res.json()) as { rows: DiagRow[] }).rows;
  } catch (e) {
    root.replaceChildren(el("p", { class: "small", role: "status" }, `The diagnosis could not be loaded (${(e as Error).message}).`));
    return;
  }
  const word: Record<DiagRow["status"], string> = { open: "Open", testing: "Partly explained", confirmed: "Cause confirmed by experiment", refuted: "First suspect ruled out", improved: "Partly fixed in the renders above", fixed: "Fixed in the renders above" };
  root.replaceChildren(...rows.map((r) =>
    el("article", { class: `diag ${r.status}` },
      el("div", { class: "diag-head" }, el("span", { class: "diag-status" }, word[r.status]), el("h4", {}, r.symptom)),
      el("dl", {},
        el("dt", {}, "Evidence"), el("dd", {}, r.evidence),
        el("dt", {}, "Suspected cause"), el("dd", {}, r.suspect),
        ...(r.result ? [el("dt", {}, "Experiment"), el("dd", {}, r.result)] : []),
        ...(r.fix ? [el("dt", {}, "Fix"), el("dd", {}, r.fix)] : []),
        ...(r.status === "confirmed" || r.status === "fixed" ? [] : [el("dt", {}, "Next test"), el("dd", {}, r.test)])))));
}
