import "@fontsource-variable/fraunces/full.css";
import "@fontsource-variable/inter";
import "@fontsource/jetbrains-mono/400.css";
import "@fontsource/jetbrains-mono/500.css";
import "./styles.css";
import drumSpec from "../../../crates/fs-couple/examples/percussion/estimated_drum.fsd?raw";
import cymbalSpec from "../../../crates/fs-couple/examples/percussion/estimated-splash.profile?raw";
import { GroupPlayer, claimPlayback } from "./audio";
import { ClipGroupView, type Manifest } from "./clipgroup";
import { LineChart, type Series } from "./chart";
import { BoardView, buildKeyboard, type PianoModel } from "./board3d";
import { ModalViewer } from "./modal3d";
import { SoundboardMixer } from "./mixer";
import { dpr, el, fmt, fmtSci, hidpi, json, noteName, onVisible, sciNode } from "./util";

type SeriesData = Record<string, Record<string, number[] | number>>;

const $ = <T extends Element = HTMLElement>(sel: string, root: ParentNode = document): T | null => root.querySelector<T>(sel);
const $$ = <T extends Element = HTMLElement>(sel: string, root: ParentNode = document): T[] => [...root.querySelectorAll<T>(sel)];

const manifestP = json<Manifest & { performances?: Record<string, [number, string, number, number][]> }>("/data/manifest.json");
let seriesCache: Promise<SeriesData> | null = null;
const seriesP = (): Promise<SeriesData> => (seriesCache ??= json<SeriesData>("/data/series.json"));
const idle = (cb: () => void): void => {
  const ric = (window as unknown as { requestIdleCallback?: (f: () => void, o?: { timeout: number }) => number }).requestIdleCallback;
  if (ric) ric(cb, { timeout: 1500 });
  else setTimeout(cb, 300);
};
const pianoP = json<PianoModel>("/data/piano-model.json");

// ------------------------------------------------------------------- navigation
function initNav(): void {
  const links = $$<HTMLAnchorElement>(".topnav a");
  const byId = new Map(links.map((a) => [a.getAttribute("href")!.slice(1), a]));
  const io = new IntersectionObserver((entries) => {
    for (const e of entries) {
      if (!e.isIntersecting) continue;
      for (const a of links) a.removeAttribute("aria-current");
      const a = byId.get(e.target.id);
      if (a) {
        a.setAttribute("aria-current", "true");
        // Scroll only the nav strip sideways; never the page (it would cancel a jump in progress).
        const nav = a.parentElement!;
        const left = a.offsetLeft - nav.clientWidth / 2 + a.clientWidth / 2;
        nav.scrollTo({ left: Math.max(0, left), behavior: "instant" as ScrollBehavior });
      }
    }
  }, { rootMargin: "-45% 0px -50% 0px" });
  for (const s of $$("main > section[id]")) io.observe(s);
}

// ------------------------------------------------------------------------ hero
async function initHero(): Promise<void> {
  const [manifest, model] = await Promise.all([manifestP, pianoP]);
  const clip = manifest.clips["hero"];
  if (!clip) return;
  const ev = clip.evidence;
  const facts = $("#hero-facts")!;
  const stat = (k: string, v: string | Node, unit = ""): HTMLElement => el("div", { class: "stat" }, el("dt", {}, k), el("dd", {}, v, unit ? el("small", {}, unit) : null));
  facts.append(
    stat("string partials", fmt(ev.string_modes ?? 0)),
    stat("board modes", String(ev.board_modes ?? "—")),
    stat("render", fmt((ev.wall_s ?? 0) / 60, 2), "min"),
    stat("energy closure", sciNode(ev.energy_closure_j ?? 0), "J"),
  );
  $("#hero-render-note")!.textContent = `rendered in ${fmt((ev.wall_s ?? 0) / 60, 2)} min for ${clip.duration_s.toFixed(1)} s of sound (${fmt(ev.wall_audio_ratio ?? 0)}× slower than real time)`;

  const player = new GroupPlayer({ hero: "/" + clip.audio }, "hero");
  const btn = $<HTMLButtonElement>("#hero-play")!;
  const label = btn.querySelector("span")!;
  player.onState = (s) => {
    label.textContent = s === "playing" ? "Pause" : s === "loading" ? "Loading…" : s === "paused" ? "Resume the offline grand" : "Hear the offline grand";
    btn.querySelector("svg")!.innerHTML = s === "playing" ? '<path d="M4 2.5h3v11H4zM9 2.5h3v11H9z" fill="currentColor"/>' : '<path d="M4 2.5v11l9-5.5z" fill="currentColor"/>';
  };
  btn.addEventListener("click", () => player.toggle());

  const stage = $("#hero-stage")!;
  const board = new BoardView(stage, model);
  onVisible(stage, () => board.mount(), "100px");
  const keys = buildKeyboard($("#hero-kbd")!);
  const perf = manifest.performances?.["hero"] ?? [];
  const scope = $<HTMLCanvasElement>("#hero-scope")!;
  const drawScope = (t: number): void => {
    hidpi(scope, 0.12, 64);
    const g = scope.getContext("2d")!;
    const { width: W, height: H } = scope;
    g.fillStyle = "#07080b";
    g.fillRect(0, 0, W, H);
    const env = clip.envelope_pa;
    const n = env.hi.length;
    const peak = clip.peak_pa;
    const px = (t / clip.duration_s) * W;
    for (let i = 0; i < n; i++) {
      const x = (i / n) * W;
      const top = H / 2 - (env.hi[i]! / peak) * (H / 2 - 3);
      const bot = H / 2 - (env.lo[i]! / peak) * (H / 2 - 3);
      g.fillStyle = x <= px ? "rgba(217,169,79,0.95)" : "rgba(201,194,179,0.35)";
      g.fillRect(x, top, Math.max(1, W / n + 0.5), Math.max(1, bot - top));
    }
    g.fillStyle = "rgba(143,137,125,0.95)";
    g.font = `${10 * dpr()}px "JetBrains Mono", monospace`;
    g.fillText(`left microphone · peak ${fmt(peak)} Pa`, 8 * dpr(), 13 * dpr());
  };
  new ResizeObserver(() => drawScope(player.time)).observe(scope);
  const clock = $("#hero-clock")!;
  const hud = $("#hero-hud-right")!;
  player.onTime = (t) => {
    drawScope(t);
    clock.textContent = `${t.toFixed(2)} s`;
    let sustain = 0;
    const down = new Map<number, { on: number; off: number; v: number }>();
    for (const [time, ev2, key, val] of perf) {
      if (time > t) break;
      if (ev2 === "sustain") sustain = val;
      else if (ev2 === "note_on") down.set(key, { on: time, off: Infinity, v: val });
      else if (ev2 === "note_off") {
        const d = down.get(key);
        if (d) d.off = time;
      }
    }
    const act = new Map<number, number>();
    const held: string[] = [];
    for (const [key, d] of down) {
      const age = t - d.on;
      let a = Math.exp(-age / 3.5) * Math.min(1, d.v / 1.7);
      if (d.off < t) a *= sustain > 0.5 ? 1 : Math.exp(-(t - d.off) / 0.12);
      else held.push(noteName(key));
      if (a > 0.02) act.set(key, a);
    }
    for (const [m, k] of keys) k.classList.toggle("on", (down.get(m)?.off ?? -1) > t);
    board.setActivity(act);
    hud.textContent = `${sustain > 0.5 ? "sustain pedal down" : "pedal up"}${held.length ? " · keys " + held.join(" ") : ""}`;
  };
  drawScope(0);
}

// ---------------------------------------------------------------- clip groups
async function initGroups(): Promise<void> {
  const manifest = await manifestP;
  const views = new Map<string, ClipGroupView>();
  for (const article of $$("[data-group]")) {
    const id = article.dataset.group!;
    const host = $(".clip-host", article);
    const pending = $(`[data-group-pending="${id}"]`);
    if (!manifest.groups[id]) {
      if (host && !pending) host.append(el("p", { class: "small muted" }, "This comparison is not available in this build."));
      continue;
    }
    pending?.remove();
    if (!host) continue;
    views.set(id, new ClipGroupView(host, manifest, id, id === "piano-silent" ? { listenCloserDb: 18 } : {}));
  }
}

// ----------------------------------------------------------------------- charts
const SPEED = ["#8fb9d6", "#d9a94f", "#e86a3c"];
async function initCharts(): Promise<void> {
  const [manifest, model] = await Promise.all([manifestP, pianoP]);
  const seriesHost = $('[data-chart="microscope"]');
  const series = await new Promise<SeriesData>((resolve) => {
    if (!seriesHost) resolve({});
    else onVisible(seriesHost, () => void seriesP().then(resolve), "1200px");
  });

  const pitch = $<HTMLCanvasElement>('[data-chart="pitch"]');
  const tracks = manifest.groups["hard"]?.pitch_tracks;
  if (pitch && tracks) onVisible(pitch, () => {
    const s: Series[] = [];
    const speeds = ["08", "25", "50"];
    const label = (v: string): string => `${v === "08" ? "0.8" : v === "25" ? "2.5" : "5"} m/s`;
    speeds.forEach((v, i) => {
      const nl = tracks[`hard-nl-${v}`];
      const ln = tracks[`hard-lin-${v}`];
      if (ln) s.push({ label: `linear ${label(v)}`, color: SPEED[i]!, x: ln.t, y: ln.f, dash: [5, 4], width: 1.4 });
      if (nl) s.push({ label: `stretching ${label(v)}`, color: SPEED[i]!, x: nl.t, y: nl.f, width: 2.2 });
    });
    new LineChart(pitch, { xLabel: "time [s]", yLabel: "fundamental [Hz]", aspect: 0.58, xFormat: (v) => v.toFixed(2), yFormat: (v) => v.toFixed(1) }, $('[data-legend="pitch"]')).set(s);
  });

  const micro = $<HTMLCanvasElement>('[data-chart="microscope"]');
  const d = series["drum_stretch_4ms"];
  if (micro && d) onVisible(micro, () => {
    const t = (d["time_s"] as number[]).map((v) => v * 1000);
    const stored = d["total_energy_j"] as number[];
    const lost = d["loss_cumulative_j"] as number[];
    const sum = stored.map((v, i) => v + lost[i]!);
    const stretch = d["head_stretching_energy_j"] as number[] | undefined;
    const s: Series[] = [
      { label: "stored in the instrument", color: "#63d2c0", x: t, y: stored.map((v) => v * 1000), width: 2 },
      { label: "dissipated (cumulative)", color: "#e86a3c", x: t, y: lost.map((v) => v * 1000), width: 2 },
      { label: "stored + dissipated", color: "#d9a94f", x: t, y: sum.map((v) => v * 1000), dash: [6, 4], width: 1.6 },
    ];
    if (stretch) s.push({ label: "head stretching energy", color: "#8fb9d6", x: t, y: stretch.map((v) => v * 1000), width: 1.4 });
    new LineChart(micro, { xLabel: "time [ms]", yLabel: "energy [mJ]", aspect: 0.62, xFormat: (v) => v.toFixed(0) }, $('[data-legend="microscope"]')).set(s);
    const worst = d["balance_abs_max_j"] as number | undefined;
    const cap = micro.closest("figure")?.querySelector("figcaption");
    if (cap && worst !== undefined) cap.replaceChildren(`The first 50 ms of a 4 m/s strike on the stretching drum, 25,000 steps of 2 µs. Stored energy (stick, heads, air, contact) and cumulative loss always add up to the stick's initial kinetic energy of ${fmt((stored[0] ?? 0) * 1000)} mJ. The largest single-step balance residual the integrator reported was `, sciNode(worst), " J.");
  });

  const snare = $<HTMLCanvasElement>('[data-chart="snare-ledger"]');
  const on = series["snare_on"], off = series["snare_off"];
  if (snare && on && off) onVisible(snare, () => {
    const mk = (src: Record<string, number[] | number>, label: string, color: string, dash?: number[]): Series => ({ label, color, dash, x: (src["time_s"] as number[]).map((v) => v * 1000), y: (src["point_b_displacement_m"] as number[]).map((v) => -v * 1e6) });
    new LineChart(snare, { xLabel: "time [ms]", yLabel: "resonant head, up [µm]", aspect: 0.56, zeroLine: true, xFormat: (v) => v.toFixed(0) }, $('[data-legend="snare-ledger"]')).set([mk(on, "wires on", "#d9a94f"), mk(off, "wires off", "#8fb9d6", [5, 4])]);
  });

  const drag = $<HTMLCanvasElement>('[data-chart="drag"]');
  const a = series["drag_ordinary"], b = series["drag_counterfactual"];
  if (drag && a && b) onVisible(drag, () => {
    const mk = (src: Record<string, number[] | number>, label: string, color: string, dash?: number[]): Series => {
      const pa = src["cavity_point_a_pa"] as number[], pb = src["cavity_point_b_pa"] as number[];
      return { label, color, dash, x: (src["time_s"] as number[]).map((v) => v * 1000), y: pa.map((v, i) => v - pb[i]!) };
    };
    new LineChart(drag, { xLabel: "time [ms]", yLabel: "p(A) − p(B) [Pa]", aspect: 0.56, zeroLine: true, xFormat: (v) => v.toFixed(0) }, $('[data-legend="drag"]')).set([mk(a, "ordinary air", "#8fb9d6"), mk(b, "10⁵ s⁻¹ momentum drag", "#e86a3c", [5, 4])]);
  });

  const bridge = $<HTMLCanvasElement>('[data-chart="bridge"]');
  if (bridge) onVisible(bridge, () => {
    const s: Series[] = [];
    (["a3-pp", "a3-mf", "a3-ff"] as const).forEach((id, i) => {
      const bv = manifest.clips[id]?.evidence["bridge_velocity_first_60ms"] as { t: number[]; v: number[] } | undefined;
      if (!bv) return;
      const peak = Math.max(...bv.v.map(Math.abs)) || 1;
      const n = bv.t.findIndex((v) => v > 0.03);
      const end = n < 0 ? bv.t.length : n;
      s.push({ label: `${["pp 0.6", "mf 2.0", "ff 4.5"][i]} m/s (peak ${fmt(peak * 1000)} mm/s)`, color: SPEED[i]!, x: bv.t.slice(0, end).map((v) => v * 1000), y: bv.v.slice(0, end).map((v) => v / peak + (2 - i) * 2.2), width: 1.4 });
    });
    new LineChart(bridge, { xLabel: "time [ms]", yLabel: "bridge velocity / peak (offset)", aspect: 0.62, xFormat: (v) => v.toFixed(0), yFormat: () => "" }, $('[data-legend="bridge"]')).set(s);
  });

  const beats = $<HTMLCanvasElement>('[data-chart="beats"]');
  if (beats) onVisible(beats, () => {
    const s: Series[] = [];
    const partial2 = (r: PianoModel["scale_et"][number]): number => 2 * r.f1 * Math.sqrt(1 + 4 * r.B) / Math.sqrt(1 + r.B);
    for (const [id, label, color, dash, sc] of [["octhi2-et", "equal-tempered", "#e86a3c", undefined, model.scale_et], ["octhi3-stretch", "stretch-tuned", "#63d2c0", [6, 4], model.scale_stretch]] as const) {
      const z = manifest.clips[id]?.evidence["zoom_1760"] as { f: number[]; db: number[] } | undefined;
      if (!z) continue;
      s.push({ label: `${label}: rendered spectrum`, color, dash: dash ? [...dash] : undefined, x: z.f, y: z.db, width: 1.6 });
      const a5 = sc.find((r) => r.midi === 81), a6 = sc.find((r) => r.midi === 93);
      if (a5 && a6) s.push({ label: `${label}: scale predicts A6 f₁ ${a6.f1.toFixed(2)} Hz, A5 partial 2 ${partial2(a5).toFixed(2)} Hz`, color, x: [a6.f1, a6.f1, NaN, partial2(a5), partial2(a5)], y: [-40, 3, NaN, -40, 3], width: 1, dash: [2, 3] });
    }
    if (s.length) new LineChart(beats, { xLabel: "frequency [Hz]", yLabel: "level [dB]", aspect: 0.5, xRange: [1752, 1774], yRange: [-40, 3], xFormat: (v) => v.toFixed(0), yFormat: (v) => v.toFixed(0) }, $('[data-legend="beats"]')).set(s);
    else beats.closest("figure")?.setAttribute("hidden", "");
  });

  const stretch = $<HTMLCanvasElement>('[data-chart="stretch"]');
  if (stretch) onVisible(stretch, () => {
    const et = model.scale_et, st = model.scale_stretch;
    const midi = et.map((r) => r.midi);
    const cents = et.map((r, i) => 1200 * Math.log2(st[i]!.f1 / r.f1));
    new LineChart(stretch, { xLabel: "key (MIDI number; 21 = A0, 60 = C4, 108 = C8)", yLabel: "stretch [cents]", y2Label: "inharmonicity B", y2Log: true, y2Range: [1e-5, 1e-2], aspect: 0.5, xFormat: (v) => v.toFixed(0), yFormat: (v) => v.toFixed(0), y2Format: (v) => fmtSci(v, 1), zeroLine: true }).set([
      { label: "stretch", color: "#63d2c0", x: midi, y: cents, width: 2 },
      { label: "B", color: "#d9a94f", x: midi, y: et.map((r) => r.B), axis: "right", width: 1.4, dash: [4, 3] },
    ]);
  });
}

// --------------------------------------------------------------------- viewers
function initViewers(): void {
  for (const stage of $$("[data-viz]")) {
    const panel = stage.parentElement!;
    const viewer = new ModalViewer(stage, {
      play: $<HTMLButtonElement>("[data-viz-play]", panel),
      scrub: $<HTMLInputElement>("[data-viz-scrub]", panel),
      clock: $("[data-viz-clock]", stage),
      scaleLabel: $("[data-viz-scale]", stage),
      label: $("[data-viz-label]", stage),
    });
    const labels: Record<string, string> = {
      "drum-slowmo": "Both heads · 4 m/s stick at (60, 10) mm · stretching",
      "cymbal-slowmo": "Splash shell · 2 m/s stick at 68 mm",
      "cymbal-ring": "Splash shell · 2 m/s stick · first second",
    };
    onVisible(stage, () => void viewer.load(stage.dataset.viz!, labels[stage.dataset.viz!]).catch(() => {
      stage.insertAdjacentHTML("beforeend", '<div class="stage-note">This solver-state replay is still being computed for this build.</div>');
    }));
    const pick = $("[data-viz-pick]", panel);
    pick?.addEventListener("click", (e) => {
      const b = (e.target as HTMLElement).closest<HTMLButtonElement>("button[data-viz-file]");
      if (!b) return;
      for (const x of $$("button", pick)) x.setAttribute("aria-pressed", String(x === b));
      void viewer.load(b.dataset.vizFile!, labels[b.dataset.vizFile!]);
    });
  }
}

// ------------------------------------------------------------------------ live
function initLive(): void {
  // The sections only build DOM here; no wasm is fetched until a Start button is pressed.
  const drum = $("#live-drum-root");
  const piano = $("#live-piano-root");
  const opened = (root: HTMLElement): void => {
    const d = root.querySelector<HTMLDetailsElement>("details.evidence");
    if (d) d.open = true;
  };
  idle(() => {
    const jobs: Promise<void>[] = [];
    if (drum) jobs.push(import("./live/drumSection").then((m) => { m.mountLiveDrum(drum); opened(drum); }));
    if (piano) jobs.push(import("./live/pianoSection").then((m) => { m.mountLivePiano(piano); opened(piano); }));
    void Promise.all(jobs).then(() => {
      if (location.hash) document.querySelector(location.hash)?.scrollIntoView();
    });
  });
  // An engine start is a deliberate act: offline clips yield to it.
  for (const root of [drum, piano]) root?.addEventListener("click", (e) => {
    const b = (e.target as HTMLElement).closest("button");
    if (b && /start/i.test(b.textContent ?? "")) claimPlayback(() => {});
  }, true);
}

// ---------------------------------------------------------------------- method
function highlightSpec(text: string): string {
  const esc = (s: string): string => s.replace(/&/g, "&amp;").replace(/</g, "&lt;");
  return text.trimEnd().split("\n").map((line) => {
    if (line.trim().startsWith("#")) return `<span class="c">${esc(line)}</span>`;
    const [k, ...rest] = line.split(",");
    if (!rest.length) return `<span class="k">${esc(line)}</span>`;
    return `<span class="k">${esc(k!)}</span>,${rest.map((v) => `<span class="n">${esc(v)}</span>`).join(",")}`;
  }).join("\n");
}

const REFUSALS: [string, string][] = [
  ["percussion reference refused: BEM reports invalid diagnostics or negative radiation power beyond roundoff at 1460 Hz, input 26, formulation PlainCbie: power_interval=(-4.9e-22, -4.9e-22) … refine the acoustic solve",
    "Every drum and cymbal render at HEAD stopped here. The diagnosis: a cymbal mode that radiates almost nothing produced a radiated power of −4.9×10⁻²² W, which is solve noise 10⁻¹⁴ times smaller than its reactive power, but the gate compared it only with floating-point summation roundoff. The fix keeps the refusal and scales it to the dense solve's own accuracy (n·ε·κ of the apparent power). Genuine sign errors, like the 1.3% Burton–Miller artifact on the drum's quadrupole modes, still refuse; the drum renders here use a 40–1000 Hz band, where the drum qualifies for plain CBIE."],
  ["percussion reference refused: observer order budget exhausted at 24; last candidate: observer holdout refused: max=0.146, rms=0.0743",
    "The cymbal's close-microphone transfer could not be fitted by a stable filter that passes its held-out audit (5% RMS) even at order 24. The cymbal clips on this page therefore use the far-field microphone, whose fit passed its audit at order 4."],
  ["grand_piano: audio sample 713760 after 32 frames in block: sample 713760, event 83: hammer/key has not returned for another strike",
    "A mistake in our own score: Satie's melody holds F♯4 while the left-hand chord also asks for F♯4. One physical key cannot be struck again while it is held, so the performance was refused before the 41-minute render. The chord now leaves out the held key."],
  ["grand_piano: adaptive bridge response work estimate 9667910336 exceeds 8000000000; use fewer source modes, retained modes or admitted bridge ports",
    "An adaptive soundboard reduction to 3.2 kHz was refused against its declared work budget before it started. The renders here use an explicit 128-mode Ritz reduction of the 0–2.4 kHz board instead."],
  ["percussion reference refused: percussion mechanics: Owner(\"FS-PLATE-BAD-SECTION: shell state exceeds its declared linearized rigid-rotation limit\")",
    "We let the splash rock on its felts. A 2 m/s off-centre stroke exceeded the declared small-angle rigid-rotation limit at 0.03 rad, and again at the model's largest admitted limit of 0.1 rad. The linearization was refused rather than extrapolated. The rocking clip on this page is a gentler 0.6 m/s stroke, which stays inside 0.1 rad."],
];

async function initMethod(): Promise<void> {
  $("#spec-drum")!.innerHTML = highlightSpec(drumSpec);
  $("#spec-cymbal")!.innerHTML = highlightSpec(cymbalSpec);
  const pm = await pianoP;
  const table = $("#scale-table")!;
  const cols: [string, (r: PianoModel["scale_et"][number]) => string][] = [
    ["key", (r) => noteName(r.midi)], ["strings", (r) => String(r.unison)], ["length", (r) => `${r.L.toFixed(3)} m`],
    ["mass/length", (r) => `${(r.mu * 1000).toFixed(2)} g/m`], ["tension", (r) => `${r.T.toFixed(0)} N`],
    ["stiffness EI", (r) => r.EI.toExponential(2) + " N·m²"], ["B", (r) => r.B.toExponential(2)], ["f₁", (r) => `${r.f1.toFixed(2)} Hz`],
  ];
  table.append(el("thead", {}, el("tr", {}, ...cols.map(([h]) => el("th", {}, h)))));
  const body = el("tbody");
  for (const midi of [21, 33, 48, 60, 69, 96]) {
    const r = pm.scale_et.find((x) => x.midi === midi);
    if (r) body.append(el("tr", {}, ...cols.map(([, f]) => el("td", { class: "num" }, f(r)))));
  }
  table.append(body);
  const list = $("#refusals")!;
  for (const [msg, why] of REFUSALS) list.append(el("div", { class: "refusal" }, el("p", {}, why), el("pre", {}, msg)));
  const manifest = await manifestP;
  $("#commit")!.textContent = manifest.commit;
  const clips = Object.values(manifest.clips);
  const cpu = clips.reduce((a, c) => a + (c.evidence.wall_s ?? 0), 0);
  const audio = clips.reduce((a, c) => a + c.duration_s, 0);
  const closures = clips.map((c) => Math.abs(c.evidence.energy_closure_j ?? 0)).filter((v) => v > 0);
  const stats = $("#method-stats")!;
  const stat = (k: string, v: string | Node, unit = ""): HTMLElement => el("div", { class: "stat" }, el("dt", {}, k), el("dd", {}, v, unit ? el("small", {}, unit) : null));
  stats.append(
    stat("offline clips", String(clips.length)),
    stat("seconds of sound", audio.toFixed(1), "s"),
    stat("CPU time to render", (cpu / 3600).toFixed(1), "h"),
    stat("worst piano energy closure", sciNode(Math.max(...closures)), "J"),
    stat("comparison groups", String(Object.keys(manifest.groups).length)),
  );
}

// --------------------------------------------------------------------- mixer
function initMixer(): void {
  const host = $("[data-mixer]");
  if (host) onVisible(host, () => new SoundboardMixer(host));
}

async function initRegisterNote(): Promise<void> {
  const model = await pianoP;
  const a0 = model.scale_et.find((r) => r.midi === 21);
  const p = $('[data-group="piano-register"] p.muted');
  if (!a0 || !p) return;
  const cents = (n: number): number => 1200 * Math.log2(Math.sqrt(1 + a0.B * n * n) / Math.sqrt(1 + a0.B));
  p.textContent = `A0 is a single ${a0.L.toFixed(2)} m wound string with B = ${fmtSci(a0.B)}, so its 20th partial lies ${cents(20).toFixed(0)} cents and its 40th ${cents(40).toFixed(0)} cents sharp of the harmonic series. C6 has three short plain strings and a much shorter hammer contact.`;
}

function initVerisim(): void {
  const root = $("#verisim-root");
  if (root) onVisible(root, () => void import("./verisim").then((m) => m.mountVerisim(root)), "800px");
}

initNav();
void initHero();
initVerisim();
void initGroups();
void initCharts();
initViewers();
initLive();
void initMethod();
initMixer();
void initRegisterNote();
