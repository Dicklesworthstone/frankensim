import { GroupPlayer } from "./audio";
import { dpr, el, fmt, hidpi } from "./util";

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
interface Verisim {
  pairs: Pair[];
  piano_ruler: { real_pp_vs_real_mf_db: number; real_mf_vs_real_ff_db: number; note: string };
  piano_dynamics: { real_db: Record<string, number>; ours_db: Record<string, number>; real_high_band_db: Record<string, number>; ours_high_band_db: Record<string, number>; note: string };
  cymbal_ruler: { reference: string; rows: { label: string; rms_gap_db: number; real: boolean }[]; note: string };
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
};

function show(v: Descriptor[keyof Descriptor], unit: string, digits = 3): string {
  if (v === null || v === undefined) return "—";
  if (Array.isArray(v)) return v.map((x) => (typeof x === "number" ? fmt(x, 2) : "—")).join(" · ") + (unit ? ` ${unit}` : "");
  if (typeof v !== "number") return String(v);
  return `${Math.abs(v) < 1e-2 && v !== 0 ? v.toExponential(2) : fmt(v, digits)}${unit ? ` ${unit}` : ""}`;
}

const UNITS: Record<string, string> = { f1_hz: "Hz", high_band_db: "dB", attack_rise_ms: "ms", partial_decay_db_s: "dB/s", octave_decay_db_s: "dB/s" };

/** Bars against a ruler made only of real recordings. */
function rulerChart(canvas: HTMLCanvasElement, rows: { label: string; value: number; real: boolean }[], unit: string): void {
  const draw = (): void => {
    hidpi(canvas, 0, rows.length * 34 + 36);
    const g = canvas.getContext("2d")!;
    const d = dpr();
    const { width: W, height: H } = canvas;
    g.fillStyle = "#07080b";
    g.fillRect(0, 0, W, H);
    const max = Math.max(...rows.map((r) => r.value)) * 1.12;
    const left = Math.min(W * 0.42, 260 * d);
    g.font = `${11.5 * d}px "JetBrains Mono", monospace`;
    g.textBaseline = "middle";
    rows.forEach((r, i) => {
      const y = 16 * d + i * 34 * d;
      g.fillStyle = "rgba(201,194,179,0.85)";
      g.textAlign = "right";
      g.fillText(r.label, left - 10 * d, y + 11 * d);
      const w = ((W - left - 70 * d) * r.value) / max;
      g.fillStyle = r.real ? "rgba(143,185,214,0.85)" : "#d9a94f";
      g.fillRect(left, y + 3 * d, w, 16 * d);
      g.textAlign = "left";
      g.fillStyle = "#f1eadb";
      g.fillText(`${r.value.toFixed(1)} ${unit}`, left + w + 8 * d, y + 11 * d);
    });
  };
  new ResizeObserver(draw).observe(canvas);
  draw();
}

export async function mountVerisim(root: HTMLElement): Promise<void> {
  const data = (await (await fetch("/data/verisim.json")).json()) as Verisim;
  const pairs = data.pairs;
  let current = pairs[0]!;
  let blind = false;
  let blindOrder: ("real" | "ours")[] = ["real", "ours"];

  const seg = el("div", { class: "seg", role: "group", "aria-label": "Choose a comparison" });
  const playA = el("button", { class: "btn sm", type: "button" });
  const playB = el("button", { class: "btn sm", type: "button" });
  const stop = el("button", { class: "btn ghost sm", type: "button" }, "Stop");
  const blindBox = el("label", { class: "toggle" }, el("input", { type: "checkbox" }), "Blind test: hide which is which");
  const reveal = el("div", { class: "small", "aria-live": "polite" });
  const tabs = el("div", { class: "seg", role: "tablist", "aria-label": "Spectrogram" });
  const img = el("img", { alt: "", style: "width:100%; height:auto; display:block; image-rendering:auto" });
  const axis = el("div", { class: "vs-axis", "aria-hidden": "true" });
  const imgBox = el("div", { class: "clip-view", style: "position:relative" }, img, axis);
  const legend = el("p", { class: "small muted" });
  const table = el("table");
  const tableBox = el("div", { class: "table-scroll", tabindex: "0", "aria-label": "Scorecard" }, table);
  const source = el("p", { class: "small muted" });

  let view: "real" | "ours" | "diff" = "real";
  const viewBtns: Record<string, HTMLButtonElement> = {};
  for (const [k, label] of [["real", "Real recording"], ["ours", "FrankenSim"], ["diff", "Difference"]] as const) {
    const b = el("button", { type: "button", role: "tab", "aria-pressed": String(k === view) }, label);
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
    player = new GroupPlayer({ real: "/" + current.audio.real, ours: "/" + current.audio.ours }, "real");
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
  blindBox.querySelector("input")!.addEventListener("change", (e) => {
    blind = (e.target as HTMLInputElement).checked;
    blindOrder = Math.random() < 0.5 ? ["real", "ours"] : ["ours", "real"];
    reveal.replaceChildren();
    if (blind) {
      const ask = el("div", { class: "btn-row" }, el("span", {}, "Which one is the real piano or cymbal?"));
      for (const slot of [0, 1] as const) {
        const b = el("button", { class: "btn sm", type: "button" }, `Sound ${slot === 0 ? "X" : "Y"}`);
        b.addEventListener("click", () => {
          const right = which(slot) === "real";
          reveal.replaceChildren(el("p", {}, right ? "Correct: that was the recording." : "That one was FrankenSim's render.", ` X was ${which(0) === "real" ? "the recording" : "FrankenSim"}, Y was ${which(1) === "real" ? "the recording" : "FrankenSim"}. Your answer stays in this browser; nothing is collected.`));
        });
        ask.append(b);
      }
      reveal.append(ask);
      view = "diff";
    }
    render();
  });

  const render = (): void => {
    for (const [k, b] of Object.entries(viewBtns)) b.setAttribute("aria-pressed", String(k === view));
    // In a blind test the labelled spectrograms would give the answer away.
    viewBtns["real"]!.disabled = blind;
    viewBtns["ours"]!.disabled = blind;
    if (blind && view !== "diff") view = "diff";
    img.src = "/" + current.images[view];
    img.alt = `${view === "diff" ? "Difference" : view === "real" ? "Real recording" : "FrankenSim"} spectrogram, ${current.title}, 40 Hz to 16 kHz on a log axis`;
    legend.textContent = view === "diff"
      ? `Difference in dB, FrankenSim minus the recording, after matching early loudness: orange = FrankenSim louder, blue = FrankenSim quieter (±${current.images.diff_range_db} dB). Black = both below the shared floor.`
      : `${current.duration_s} s from the onset, 40 Hz–16 kHz (log), ${current.images.span_db} dB colour range shared by both sounds.`;
    const rows: HTMLElement[] = [el("tr", {}, el("th", {}, "Measure"), el("th", {}, current.kind === "piano" ? "Real Steinway" : "Real splash"), el("th", {}, "FrankenSim"))];
    for (const [k, dsc] of Object.entries(current.descriptors)) {
      const unit = UNITS[k] ?? "";
      const gap = dsc.rms_gap_db;
      const measure = el("td", {}, el("div", {}, LABELS[k] ?? k), el("div", { class: "small muted", style: "margin-top:3px" }, dsc.note ?? ""));
      if (gap !== undefined && gap !== null) {
        rows.push(el("tr", {}, measure, el("td", { class: "num", colspan: "2" }, `${gap.toFixed(1)} dB apart`)));
        continue;
      }
      const cls = Array.isArray(dsc.real) ? "num small" : "num";
      rows.push(el("tr", {}, measure, el("td", { class: cls }, show(dsc.real, unit)), el("td", { class: cls }, show(dsc.ours, unit))));
    }
    table.replaceChildren(el("thead", {}, rows[0]!), el("tbody", {}, ...rows.slice(1)));
    source.replaceChildren(
      el("span", {}, `Recording: ${current.real.label}, `), el("a", { href: current.real.url, rel: "noopener" }, "source file"),
      el("span", {}, `. ${current.real.recording}. ${current.audio.level_match}. Rights: ${current.real.rights}.`),
    );
    labelButtons();
  };

  for (const p of pairs) {
    const b = el("button", { type: "button", "aria-pressed": String(p === current) }, p.title.split(":")[0]!);
    b.addEventListener("click", () => {
      current = p;
      for (const x of seg.querySelectorAll("button")) x.setAttribute("aria-pressed", String(x === b));
      makePlayer();
      render();
    });
    seg.append(b);
  }

  const pianoRuler = el("canvas", { "aria-label": "Partial-envelope gaps: FrankenSim versus the real piano, and the real piano versus itself across dynamics" });
  const cymRuler = el("canvas", { "aria-label": "Distance from the real splash cymbal: real crashes, the same splash at other dynamics, and FrankenSim" });
  const gap = (id: string): number => (pairs.find((p) => p.id === id)?.descriptors["partial_envelope_db"]?.rms_gap_db ?? NaN) as number;

  root.append(
    el("div", { class: "panel" },
      el("div", { class: "panel-head" }, el("h3", {}, "Listen and look, side by side"), el("div", { class: "tags" }, el("span", { class: "tag offline" }, "Offline reference"), el("span", { class: "tag plain" }, "Real recording"))),
      seg,
      el("div", { class: "stack", style: "margin-top: 16px" },
        el("div", { class: "stack" },
          el("div", { class: "btn-row" }, playA, playB, stop, blindBox),
          reveal, tabs, imgBox, legend),
        el("div", { class: "stack" }, el("h4", {}, "Scorecard"), tableBox, source)),
    ),
    el("div", { class: "grid two", style: "margin-top: 22px" },
      el("figure", { class: "panel" },
        el("h3", {}, "Piano: a ruler made of the real piano"),
        el("p", { class: "muted small", style: "margin: 8px 0 14px" }, "How far apart are the partial-level envelopes (partials 2–16 relative to partial 1)? The blue bars are the real Steinway against itself, one dynamic apart: what a change of touch does. The brass bars are FrankenSim against the real piano at the same dynamic."),
        el("div", { class: "chart" }, pianoRuler)),
      el("figure", { class: "panel" },
        el("h3", {}, "Cymbal: a ruler made of real cymbals"),
        el("p", { class: "muted small", style: "margin: 8px 0 14px" }, "How far is each sound from the real splash struck mf, in octave-band balance over the first 250 ms? Blue bars are real cymbals from the same anechoic session, brass bars are FrankenSim. A FrankenSim bar inside the blue spread would sound as close to this splash as another real cymbal does."),
        el("div", { class: "chart" }, cymRuler)),
    ),
  );
  rulerChart(pianoRuler, [
    { label: "real pp vs real mf", value: data.piano_ruler.real_pp_vs_real_mf_db, real: true },
    { label: "real mf vs real ff", value: data.piano_ruler.real_mf_vs_real_ff_db, real: true },
    { label: "FrankenSim vs real, pp", value: gap("a3-pp"), real: false },
    { label: "FrankenSim vs real, mf", value: gap("a3-mf"), real: false },
    { label: "FrankenSim vs real, ff", value: gap("a3-ff"), real: false },
    { label: "FrankenSim vs real, C6 mf", value: gap("c6-mf"), real: false },
  ].filter((r) => Number.isFinite(r.value)), "dB");
  rulerChart(cymRuler, data.cymbal_ruler.rows.map((r) => ({ label: r.label, value: r.rms_gap_db, real: r.real })), "dB");
  makePlayer();
  render();
}
