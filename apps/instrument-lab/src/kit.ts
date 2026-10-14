import { GroupPlayer } from "./audio";
import { el, fmt } from "./util";

interface Pair {
  id: string;
  title: string;
  duration_s: number;
  real: { file: string; url: string; rights: string; recording: string };
  ours: { clip: string; label: string };
  audio: { real: string; ours: string; level_match: string };
  images: { real: string; ours: string; diff: string; span_db: number; diff_range_db: number };
  descriptors: Record<string, { real: number | (number | null)[] | null; ours: number | (number | null)[] | null; note: string; bands_hz?: number[]; rms_gap_db?: number }>;
  ruler: { real_same_stroke_db: number; real_other_dynamic_db: number; ours_db: number; decay_same_stroke_oct: number | null; decay_ours_oct: number | null; strokes: string[] };
}
interface KitData {
  pairs: Pair[];
  rights: string;
  source: string;
  groove?: { real: string; ours: string; pattern: string; level_match: string };
}

const SCALAR: [string, string, (v: number) => string][] = [
  ["fundamental_hz", "Strongest low component", (v) => `${fmt(v, 3)} Hz`],
  ["centroid_hz", "Spectral centroid, first 250 ms", (v) => `${fmt(v, 3)} Hz`],
  ["above_1k_pct", "Energy above 1 kHz", (v) => `${v.toFixed(1)} %`],
  ["above_4k_pct", "Energy above 4 kHz", (v) => `${v.toFixed(1)} %`],
];

function bandLabel(f: number): string {
  return f >= 1000 ? `${f / 1000}k` : String(f);
}

/** Real kit vs FrankenSim's kit, piece by piece, plus the same groove on both. */
export async function mountKit(root: HTMLElement): Promise<void> {
  let data: KitData;
  try {
    const res = await fetch("/data/kit.json");
    if (!res.ok) throw new Error(`HTTP ${res.status}`);
    data = (await res.json()) as KitData;
  } catch (e) {
    root.style.minHeight = "0";
    root.replaceChildren(el("p", { class: "panel small", role: "status" }, `The kit data could not be loaded (${(e as Error).message}). Reload to try again.`));
    return;
  }
  const pairs = data.pairs;
  if (!pairs.length) return;
  let current = pairs[0]!;
  let view: "real" | "ours" | "diff" = "ours";

  const seg = el("div", { class: "seg kit-pieces", role: "group", "aria-label": "Choose a piece of the kit" });
  const playOurs = el("button", { class: "btn sm", type: "button" });
  const playReal = el("button", { class: "btn sm", type: "button" });
  const stop = el("button", { class: "btn ghost sm", type: "button" }, "Stop");
  const tabs = el("div", { class: "seg vs-tabs", role: "group", "aria-label": "Spectrogram to show" });
  const img = el("img", { alt: "", class: "vs-img", width: "700", height: "260", decoding: "async" });
  const axis = el("div", { class: "vs-axis", "aria-hidden": "true" });
  for (const f of [100, 250, 500, 1000, 2000, 4000, 8000]) {
    const y = 1 - Math.log(f / 40) / Math.log(16000 / 40);
    axis.append(el("span", { style: `position:absolute; right:6px; top:calc(${(y * 100).toFixed(2)}% - 7px); font: 500 10.5px var(--font-mono); color: rgba(241,234,219,0.8); text-shadow: 0 1px 2px #000` }, bandLabel(f)));
  }
  const legend = el("p", { class: "small muted" });
  const table = el("table", { class: "vs-table" });
  const ruler = el("div", { class: "kit-ruler" });
  const source = el("p", { class: "small muted" });
  const viewBtns: Record<string, HTMLButtonElement> = {};
  for (const [k, label] of [["ours", "FrankenSim"], ["real", "Real kit"], ["diff", "Difference"]] as const) {
    const b = el("button", { type: "button", "aria-pressed": String(k === view) }, label);
    b.addEventListener("click", () => { view = k; render(); });
    viewBtns[k] = b;
    tabs.append(b);
  }

  let player: GroupPlayer | null = null;
  const label = (): void => {
    const playing = player?.state === "playing" ? player.current : null;
    playOurs.textContent = `${playing === "ours" ? "■" : "▶"} FrankenSim`;
    playReal.textContent = `${playing === "real" ? "■" : "▶"} Real kit`;
  };
  const makePlayer = (): void => {
    player?.pause();
    player = new GroupPlayer({ real: "/" + current.audio.real, ours: "/" + current.audio.ours }, "ours");
    player.onState = label;
  };
  const play = (id: "real" | "ours"): void => {
    if (!player) return;
    if (player.state === "playing" && player.current === id) { player.pause(); return; }
    player.select(id);
    if (player.state !== "playing") void player.play(0).catch(() => undefined);
  };
  playOurs.addEventListener("click", () => play("ours"));
  playReal.addEventListener("click", () => play("real"));
  stop.addEventListener("click", () => player?.pause());

  const render = (): void => {
    for (const [k, b] of Object.entries(viewBtns)) b.setAttribute("aria-pressed", String(k === view));
    img.src = "/" + current.images[view];
    img.alt = `${view === "diff" ? "Difference" : view === "real" ? "Real kit" : "FrankenSim"} spectrogram of the ${current.title.toLowerCase()}, 40 Hz to 16 kHz on a log axis`;
    legend.textContent = view === "diff"
      ? `FrankenSim minus the recording, in dB after matching early loudness: orange = FrankenSim louder, blue = quieter (±${current.images.diff_range_db} dB).`
      : `${current.duration_s} s from the strike, 40 Hz–16 kHz (log), ${current.images.span_db} dB colour range shared by both.`;
    const d = current.descriptors;
    const rows: HTMLElement[] = [];
    for (const [k, name, f] of SCALAR) {
      const v = d[k];
      if (!v) continue;
      rows.push(el("tr", {}, el("td", { class: "vs-measure" }, el("div", {}, name)),
        el("td", { class: "num", "data-label": "Real kit" }, typeof v.real === "number" ? f(v.real) : "—"),
        el("td", { class: "num", "data-label": "FrankenSim" }, typeof v.ours === "number" ? f(v.ours) : "—")));
    }
    const dec = d["decay_ms"];
    if (dec && Array.isArray(dec.real) && Array.isArray(dec.ours)) {
      const cells = (arr: (number | null)[]): string => arr.map((x) => (x === null ? "—" : String(x))).join(" · ");
      rows.push(el("tr", {}, el("td", { class: "vs-measure" }, el("div", {}, "Time to fall 20 dB, ms"),
        el("div", { class: "small muted", style: "margin-top:3px" }, `octave bands ${(dec.bands_hz ?? []).map(bandLabel).join(" · ")} Hz; a dash never falls that far`)),
        el("td", { class: "num small", "data-label": "Real kit" }, cells(dec.real as (number | null)[])),
        el("td", { class: "num small", "data-label": "FrankenSim" }, cells(dec.ours as (number | null)[]))));
    }
    table.replaceChildren(el("thead", {}, el("tr", {}, el("th", {}, "Measure"), el("th", {}, "Real kit"), el("th", {}, "FrankenSim"))), el("tbody", {}, ...rows));
    const r = current.ruler;
    const max = Math.max(r.real_same_stroke_db, r.real_other_dynamic_db, r.ours_db) * 1.1 || 1;
    const bar = (name: string, v: number, real: boolean): HTMLElement =>
      el("div", { class: "kit-bar" }, el("span", { class: "kit-bar-label" }, name),
        el("span", { class: `kit-bar-fill ${real ? "real" : "ours"}`, style: `width:${((100 * v) / max).toFixed(1)}%` }),
        el("span", { class: "kit-bar-value" }, `${v.toFixed(1)} dB`));
    ruler.replaceChildren(
      el("h4", {}, "Octave-band balance, first 250 ms: distance from this real stroke"),
      bar("Another real stroke, same dynamic", r.real_same_stroke_db, true),
      bar("The same real drum, one dynamic away", r.real_other_dynamic_db, true),
      bar("FrankenSim", r.ours_db, false),
      el("p", { class: "small muted" }, r.ours_db <= r.real_other_dynamic_db
        ? "FrankenSim is closer to this stroke than the same real drum played at another dynamic."
        : `FrankenSim is ${(r.ours_db / Math.max(r.real_other_dynamic_db, 0.1)).toFixed(1)}× as far from this stroke as the same real drum at another dynamic.`));
    source.replaceChildren(el("span", {}, `Recording: ${current.real.file}, `), el("a", { href: current.real.url, rel: "noopener" }, "Salamander Drumkit"),
      el("span", {}, `. ${current.real.recording}. ${current.audio.level_match.charAt(0).toUpperCase()}${current.audio.level_match.slice(1)}. ${current.real.rights}.`));
    label();
  };

  for (const p of pairs) {
    const f = p.descriptors["fundamental_hz"]?.ours;
    const b = el("button", { type: "button", "aria-pressed": String(p === current) }, p.title, el("small", {}, typeof f === "number" ? `${fmt(f, 3)} Hz` : ""));
    b.addEventListener("click", () => {
      current = p;
      for (const x of seg.querySelectorAll("button")) x.setAttribute("aria-pressed", String(x === b));
      makePlayer();
      render();
    });
    seg.append(b);
  }

  const groove = data.groove;
  const grooveBox = groove ? (() => {
    const gp = new GroupPlayer({ real: "/" + groove.real, ours: "/" + groove.ours }, "ours");
    const gOurs = el("button", { class: "btn primary", type: "button" });
    const gReal = el("button", { class: "btn", type: "button" });
    const gl = (): void => {
      const playing = gp.state === "playing" ? gp.current : null;
      gOurs.textContent = `${playing === "ours" ? "■" : "▶"} FrankenSim's kit`;
      gReal.textContent = `${playing === "real" ? "■" : "▶"} The real kit`;
    };
    gp.onState = gl;
    const go = (id: "real" | "ours"): void => {
      if (gp.state === "playing" && gp.current === id) { gp.pause(); return; }
      gp.select(id);
      if (gp.state !== "playing") void gp.play(0).catch(() => undefined);
    };
    gOurs.addEventListener("click", () => go("ours"));
    gReal.addEventListener("click", () => go("real"));
    gl();
    return el("div", { class: "panel kit-groove" },
      el("div", { class: "panel-head" }, el("h3", {}, "The same groove on both kits"), el("div", { class: "tags" }, el("span", { class: "tag offline" }, "Offline reference"), el("span", { class: "tag plain" }, "Real recording"))),
      el("p", { class: "muted small", style: "margin: 0 0 14px" }, groove.pattern),
      el("div", { class: "btn-row" }, gOurs, gReal),
      el("p", { class: "small muted", style: "margin-top: 12px" }, groove.level_match));
  })() : null;

  root.replaceChildren(
    ...(grooveBox ? [grooveBox] : []),
    el("div", { class: "panel", style: "margin-top: 22px" },
      el("div", { class: "panel-head" }, el("h3", {}, "Piece by piece, against the real kit"), el("div", { class: "tags" }, el("span", { class: "tag offline" }, "Offline reference"), el("span", { class: "tag plain" }, "Real recording"))),
      seg,
      el("div", { class: "kit-grid", style: "margin-top: 16px" },
        el("div", { class: "stack" }, el("div", { class: "btn-row" }, playOurs, playReal, stop), tabs,
          el("div", { class: "clip-view", style: "position:relative" }, img, axis), legend),
        el("div", { class: "stack" }, el("div", { class: "table-scroll", tabindex: "0", "aria-label": "Measurements" }, table), ruler, source))),
  );
  makePlayer();
  render();
}
