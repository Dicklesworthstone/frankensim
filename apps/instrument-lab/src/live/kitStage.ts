// A playable drum kit drawn like a real one. Every piece is the fs-vkplate solver running in its
// own Worker; touching a head or a cymbal strikes it at that point, and the piece is painted with
// its own computed velocity field, so the wave you see spreading is the one you hear.

import { claimPlayback, registerPlayer } from "../audio";
import { el } from "../util";
import { createLiveAudioContext } from "./liveAudio";
import { createLiveMix, KitLane, type KitBurst, type KitField, type KitReady, type LiveMix } from "./mix";

/** Slow-motion factor for the replay of each stroke's first milliseconds. */
const SLOW = 60;

type PieceId = "kick" | "snare" | "tom-12" | "tom-14" | "hihat-closed" | "hihat-open" | "crash" | "ride";

interface PieceDef {
  id: PieceId;
  name: string;
  key: string;
  /** Default stick or beater speed, m/s. */
  speed: number;
  /** Fixed output gain (Pa to full scale), see the evidence note. */
  gain: number;
  kind: "drum" | "kick" | "cymbal" | "hat";
  /** Finish of a drum shell, or bronze tone of a cymbal. */
  finish?: string;
  /**
   * How much of one stroke the Worker times each grid over (default 112 ms). The snare's wires touch
   * the head from ~85 ms after a stroke, and its steps then cost ~1.6x as much, so a 112 ms test
   * promised ~1.5x the speed it kept under rapid tapping, and it dropped audio. Over 450 ms the
   * test reads within ~5 % of its speed with a stroke every 0.56 s.
   */
  benchMs?: number;
  /**
   * Grid coarsenings to try, finest first (default 2.5-8). The snare stops at 6x: 5x and 4x sound no
   * closer to the 1x grid than 6x for 1.2-1.9x the work, and timing 5x only to fail cost 0.4 s of start-up.
   */
  scales?: number[];
}

/**
 * The hats also try 4.5x, between 4x and 5x. On 5x a hat rings a tonal line near 558 Hz (closed
 * standard stroke 33 dB above the local spectrum, against 23 dB at 550 Hz on the 1x grid). 4.5x keeps
 * a weaker one (519 Hz, 28.5 dB) and is closer to the 1x grid than 5x in every measure: octave balance
 * 9.8 against 14.9 dB (soft stroke 16.8 / 21.3, open hat 7.2 / 10.4) and 16.5 against 20.5 dB from
 * the real hat. It costs ~1.3x 5x; 4x stays first (9.3 dB, no line) where it fits.
 * They skip 5x altogether: with the rims rubbing (fs-vkplate Friction), a closed stroke there keeps a
 * 250 Hz component ringing ~4 s, where 2.5-4.5x, 6x and 7x put 80-99 % of its energy in its first 50 ms
 * (8x rings ~6 s too, and stays only as the last resort of a slow device).
 */
const HAT_SCALES = [2.5, 3, 3.5, 4, 4.5, 6, 7, 8];

/**
 * The kick and toms stop at the grids past which they sound no closer to the 1x grid, so the CPU
 * goes to the pieces that do gain. The kick's octave balance against 1x is one plateau from 3x to 8x
 * (14.1 / 14.2 / 13.9 dB at 4x / 5x / 6x), and 6x is ~0.45x the work of 4x. The rack tom is as close
 * at 5x as at 3-4x (1.3 against 0.8-2.0 dB). The floor tom's 63-500 Hz bands stay within 0.9 dB of
 * 3.5x down to 4x but move 1.9 dB at 5x, so it stops at 4x.
 */
const DONOR = { kick: [6, 7, 8], rack: [5, 6], floor: [4, 5, 6] };

// Gains: each piece's 50 ms loudness matched to the real kit's (Salamander overheads), capped so a
// standard stroke peaks at 0.55 of full scale (the live cymbals' attacks are too sharp to match).
const PIECES: PieceDef[] = [
  { id: "kick", name: "Kick", key: "B", speed: 3.0, gain: 0.0319, kind: "kick", finish: "#5c1d1a", scales: DONOR.kick },
  { id: "snare", name: "Snare", key: "S", speed: 4.5, gain: 0.0129, kind: "drum", finish: "#c9cdd3", benchMs: 450, scales: [6, 7, 8] },
  { id: "tom-12", name: "Rack tom", key: "T", speed: 3.5, gain: 0.0276, kind: "drum", finish: "#5c1d1a", scales: DONOR.rack },
  { id: "tom-14", name: "Floor tom", key: "F", speed: 3.5, gain: 0.0267, kind: "drum", finish: "#5c1d1a", scales: DONOR.floor },
  { id: "hihat-closed", name: "Hi-hat", key: "H", speed: 2.5, gain: 0.134, kind: "hat", scales: HAT_SCALES },
  { id: "hihat-open", name: "Open hi-hat", key: "O", speed: 2.5, gain: 0.134, kind: "hat", scales: HAT_SCALES },
  { id: "crash", name: "Crash", key: "C", speed: 6.0, gain: 0.0496, kind: "cymbal" },
  { id: "ride", name: "Ride", key: "R", speed: 2.5, gain: 0.182, kind: "cymbal" },
];

/** Where each piece sits, in a 1000-wide design space: centre, radii (perspective ellipse), draw order. */
interface Spot { cx: number; cy: number; rx: number; ry: number; z: number }
const WIDE: Record<string, Spot> = {
  kick: { cx: 500, cy: 470, rx: 165, ry: 165, z: 0 },
  "tom-12": { cx: 430, cy: 268, rx: 92, ry: 54, z: 1 },
  "tom-14": { cx: 712, cy: 432, rx: 122, ry: 70, z: 2 },
  snare: { cx: 300, cy: 440, rx: 112, ry: 64, z: 3 },
  hat: { cx: 132, cy: 312, rx: 122, ry: 38, z: 4 },
  crash: { cx: 245, cy: 128, rx: 158, ry: 50, z: 5 },
  ride: { cx: 772, cy: 190, rx: 176, ry: 56, z: 6 },
};
const TALL: Record<string, Spot> = {
  kick: { cx: 500, cy: 900, rx: 210, ry: 210, z: 0 },
  "tom-12": { cx: 500, cy: 470, rx: 150, ry: 88, z: 1 },
  "tom-14": { cx: 760, cy: 700, rx: 190, ry: 110, z: 2 },
  snare: { cx: 245, cy: 700, rx: 185, ry: 106, z: 3 },
  hat: { cx: 210, cy: 410, rx: 190, ry: 62, z: 4 },
  crash: { cx: 280, cy: 170, rx: 250, ry: 80, z: 5 },
  ride: { cx: 735, cy: 290, rx: 245, ry: 80, z: 6 },
};

interface Live {
  def: PieceDef;
  lane: KitLane | null;
  ready: KitReady | null;
  /** grid cell of every surface node, and the image the field is painted into */
  cells: Int32Array | null;
  m: number;
  img: ImageData | null;
  off: HTMLCanvasElement | null;
  field: KitField | null;
  vref: number;
  kinetic: number;
  kpeak: number;
  force: number;
  lastStrike: { x: number; y: number; at: number } | null;
  /** a stroke's first milliseconds, replayed SLOW times slower from `at` */
  replay: { b: KitBurst; at: number; n: number } | null;
  /** why the piece's engine stopped after it was ready (it has left the mix), or null */
  failure: string | null;
}

export function mountLiveKit(root: HTMLElement): void {
  const canvas = el("canvas", { class: "kit-canvas", "aria-label": "A drum kit: press a head or a cymbal to strike it there" });
  const start = el("button", { class: "btn primary kit-start", type: "button" }, "▶ Start the kit");
  const overlay = el("div", { class: "kit-overlay" }, start, el("p", { class: "small muted" }, "Eight pieces, each solved in its own Worker as you play."));
  const slowLabel = el("div", { class: "kit-slow", hidden: "" }, `first ${(32 * 0.25).toFixed(0)} ms after the stroke, ${SLOW}× slower`);
  const stage = el("div", { class: "kit-stage" }, canvas, overlay, slowLabel);
  const hud = el("div", { class: "kit-hud", role: "status" });
  const hatMode = el("div", { class: "seg", role: "group", "aria-label": "Hi-hat pedal" });
  let hatOpen = false;
  for (const [open, label] of [[false, "Pedal down (closed)"], [true, "Pedal up (open)"]] as const) {
    const b = el("button", { type: "button", "aria-pressed": String(open === hatOpen) }, label);
    b.addEventListener("click", () => {
      hatOpen = open;
      for (const x of hatMode.querySelectorAll("button")) x.setAttribute("aria-pressed", String(x === b));
      paintStatic = true;
    });
    hatMode.append(b);
  }
  const dyn = el("input", { type: "range", min: "0.4", max: "1.6", step: "0.05", value: "1", "aria-label": "Stroke strength" }) as HTMLInputElement;
  const dynLabel = el("span", { class: "num small" }, "mf");
  dyn.addEventListener("input", () => { dynLabel.textContent = dynName(Number(dyn.value)); });
  const buttons = el("div", { class: "kit-keys" });
  const controls = el("div", { class: "kit-controls" },
    el("label", { class: "kit-dyn" }, el("span", { class: "small muted" }, "Stroke"), dyn, dynLabel),
    hatMode);
  const numberInput = (value: number, min: number, max: number): HTMLInputElement =>
    el("input", { type: "number", value: String(value), min: String(min), max: String(max), step: "any", required: "" });
  const temperature = numberInput(20, 0, 50);
  const pressure = numberInput(101.325, 10, 200);
  const humidity = numberInput(0, 0, 100);
  const medium = el("select", {}, el("option", { value: "air" }, "Air"),
    el("option", { value: "mixture" }, "Custom gas mixture"));
  const basis = el("select", {}, el("option", { value: "mole" }, "Mole percentages"),
    el("option", { value: "mass" }, "Mass percentages"));
  const components = ([
    ["N2", "Nitrogen", 70], ["O2", "Oxygen", 20], ["Ar", "Argon", 0],
    ["CO2", "Carbon dioxide", 10], ["H2O", "Water vapor", 0],
  ] as const).map(([symbol, name, value]) => ({ symbol, name, input: numberInput(value, 0, 100) }));
  const labelled = (name: string, input: HTMLElement): HTMLElement => el("label", {}, el("span", {}, name), input);
  const mixtureFields = el("div", { class: "kit-atmosphere-grid", hidden: "" },
    labelled("Composition basis", basis), ...components.map(c => labelled(`${c.name} · %`, c.input)));
  const atmosphereStatus = el("p", { class: "small muted", role: "status" },
    "One atmosphere supplies every piece's cavity, air loading and microphone. Applying changes restarts the kit from rest.");
  const applyAtmosphere = el("button", { class: "btn sm", type: "submit" }, "Apply atmosphere");
  const atmosphereForm = el("form", {},
    el("div", { class: "kit-atmosphere-grid" }, labelled("Gas", medium), labelled("Temperature · °C", temperature),
      labelled("Pressure · kPa", pressure), labelled("Relative humidity in air · %", humidity)),
    mixtureFields, el("p", { class: "small muted" }, "Mixture percentages must total 100. Water vapor is a mixture component; relative humidity applies to air."),
    applyAtmosphere, atmosphereStatus);
  const atmosphere = el("details", { class: "kit-atmosphere" }, el("summary", {}, "Atmosphere"), atmosphereForm);
  const updateMedium = (): void => {
    const isMixture = medium.value === "mixture";
    humidity.disabled = isMixture;
    mixtureFields.hidden = !isMixture;
    basis.disabled = !isMixture;
    for (const component of components) component.input.disabled = !isMixture;
  };
  medium.addEventListener("change", updateMedium);
  updateMedium();

  function readAtmosphere(): { record: string; label: string } {
    if (!atmosphereForm.reportValidity()) throw new Error("Enter temperature, pressure and composition within the shown ranges.");
    const t = temperature.valueAsNumber + 273.15, p = pressure.valueAsNumber * 1000;
    if (!Number.isFinite(t) || !Number.isFinite(p)) throw new Error("Temperature and pressure must be finite.");
    if (medium.value === "air") {
      const rh = humidity.valueAsNumber / 100;
      if (!Number.isFinite(rh)) throw new Error("Relative humidity must be finite.");
      return { record: `${t} ${p} ${rh}`, label: `Air · ${temperature.value} °C · ${pressure.value} kPa · ${humidity.value}% RH` };
    }
    const selected = components.filter(c => c.input.valueAsNumber > 0);
    const total = components.reduce((sum, c) => sum + c.input.valueAsNumber, 0);
    if (!Number.isFinite(total) || Math.abs(total - 100) > 1e-6) throw new Error("Gas mixture percentages must total 100.");
    const fields = selected.flatMap(c => [c.symbol, String(c.input.valueAsNumber / 100)]);
    return { record: `mixture ${t} ${p} ${basis.value} ${selected.length} ${fields.join(" ")}`,
      label: `${selected.map(c => `${c.name} ${c.input.value}%`).join(" · ")} (${basis.value}) · ${temperature.value} °C · ${pressure.value} kPa` };
  }
  const evidence = el("details", { class: "evidence" }, el("summary", {}, "What is computed"),
    el("p", { class: "small" }, "Each piece is the same finite-difference solver as the drum kit above: tensioned Mylar heads with the air between them, von Kármán bronze shells, Hertz stick contacts, a felt beater, snare wires and a pedal-clamped hi-hat whose rims rub with dry friction. Here it runs on a grid several times coarser (picked per piece for this device when you start), so it is cheaper and loses detail above a few kHz. A piece that has gone quiet is not stepped until it is struck again."),
    el("p", { class: "small muted" }, "The colours are the computed transverse velocity of the surface (warm up, cool down), drawn from the solver state about 20 ms ahead of what you hear. Energy is the kinetic energy of the piece; force is the peak stick force of the last stroke. Levels: one gain per piece, matched on one grid to a real kit's recorded balance, except that the cymbals are capped so their very sharp attacks do not clip. A stroke's loudness changes with the grid (the crash is about 5 dB louder on a 3× grid than on 5×), so each piece's gain is rescaled for the grid it runs on, from a table of that stroke measured on every grid; it is raised only while the stroke stays under the same peak cap."));
  root.replaceChildren(el("div", { class: "panel kit-live" }, stage, el("div", { class: "kit-bar" }, controls, buttons), atmosphere, hud, evidence));
  root.style.minHeight = "0";

  const g = canvas.getContext("2d")!;
  const lives = new Map<PieceId, Live>();
  for (const def of PIECES) lives.set(def.id, { def, lane: null, ready: null, cells: null, m: 0, img: null, off: null, field: null, vref: 1e-3, kinetic: 0, kpeak: 0, force: 0, lastStrike: null, replay: null, failure: null });
  for (const def of PIECES) {
    const b = el("button", { type: "button", class: "kit-key" }, el("kbd", {}, def.key), def.name);
    b.addEventListener("click", () => hitCentre(def.id));
    buttons.append(b);
  }

  let ctx: AudioContext | null = null;
  let mix: LiveMix | null = null;
  let running = false;
  let W = 0, H = 0, S = 1, layout = WIDE, ox = 0, oy = 0;
  let paintStatic = true;
  const base = document.createElement("canvas");
  let raf = 0;
  let visible = true; // the canvas is on screen (an IntersectionObserver below keeps this)

  const spotOf = (id: PieceId): Spot => layout[id.startsWith("hihat") ? "hat" : id]!;
  const px = (s: Spot): { cx: number; cy: number; rx: number; ry: number } => ({ cx: ox + s.cx * S, cy: oy + s.cy * S, rx: s.rx * S, ry: s.ry * S });

  function resize(): void {
    const r = canvas.getBoundingClientRect();
    const d = Math.min(2, window.devicePixelRatio || 1);
    W = Math.max(2, Math.round(r.width * d));
    H = Math.max(2, Math.round(r.height * d));
    canvas.width = W; canvas.height = H;
    base.width = W; base.height = H;
    layout = r.width / r.height < 0.95 ? TALL : WIDE;
    const dh = layout === TALL ? 1150 : 660;
    S = Math.min(W / 1000, H / dh);
    ox = (W - 1000 * S) / 2;
    oy = (H - dh * S) / 2;
    paintStatic = true;
    // draw() schedules its own next frame: the one already pending is cancelled, or every resize
    // while the kit runs would start one more loop drawing the whole kit each frame
    cancelAnimationFrame(raf);
    draw();
  }
  const ro = new ResizeObserver(resize);
  ro.observe(canvas);

  // ------------------------------------------------------------- drawing the kit
  function drawStatic(): void {
    const b = base.getContext("2d")!;
    b.clearRect(0, 0, W, H);
    // floor and a warm spot from above
    const floor = b.createRadialGradient(W / 2, H * 0.62, 10, W / 2, H * 0.62, Math.max(W, H) * 0.75);
    floor.addColorStop(0, "#2a211b");
    floor.addColorStop(0.45, "#16120f");
    floor.addColorStop(1, "#08090c");
    b.fillStyle = floor;
    b.fillRect(0, 0, W, H);
    const rug = px({ cx: 500, cy: layout === TALL ? 760 : 470, rx: 470, ry: layout === TALL ? 400 : 200, z: 0 });
    b.fillStyle = "rgba(122, 40, 32, 0.22)";
    ellipse(b, rug.cx, rug.cy, rug.rx, rug.ry);
    b.fill();
    const order = [...PIECES.filter((p) => p.id !== "hihat-open")].sort((a, c) => spotOf(a.id).z - spotOf(c.id).z);
    // stands first, behind everything
    for (const p of order) if (p.kind === "cymbal" || p.kind === "hat") stand(b, px(spotOf(p.id)), p.kind === "hat");
    for (const p of order) {
      const s = px(spotOf(p.id));
      if (p.kind === "kick") kick(b, s, p.finish!);
      else if (p.kind === "drum") drum(b, s, p.finish!, p.id === "snare");
      else cymbal(b, s, p.kind === "hat");
    }
    paintStatic = false;
  }

  function draw(): void {
    if (paintStatic) drawStatic();
    g.clearRect(0, 0, W, H);
    g.drawImage(base, 0, 0);
    const now = performance.now();
    let busy = false;
    let slowTag = false;
    for (const p of PIECES) {
      const L = lives.get(p.id)!;
      if (p.id === "hihat-open" && !hatOpen && !L.kinetic) continue;
      const s = px(spotOf(p.id));
      let fade = L.kpeak > 0 ? Math.max(0, Math.min(1, (Math.log10(Math.max(L.kinetic, 1e-30) / L.kpeak) + 4.5) / 4.5)) : 0;
      let src: Float32Array | null = L.field ? L.field.v : null;
      let fresh = false;
      if (L.replay && L.off) {
        const r = L.replay, nodes = r.b.v.length / r.b.count;
        const k = Math.floor((now - r.at) / (r.b.stepMs * SLOW));
        if (k < r.b.count) {
          src = r.b.v.subarray(k * nodes, (k + 1) * nodes);
          fade = 1;
          fresh = true;
          slowTag = true;
        } else L.replay = null;
      }
      if (src && L.off && fade > 0.003) {
        paintField(L, src, fade, fresh);
        g.save();
        ellipse(g, s.cx, s.cy, s.rx * 0.985, s.ry * 0.985);
        g.clip();
        // light coated heads take the field as paint; bronze and the dark kick head let it glow
        g.globalCompositeOperation = p.kind === "drum" ? "source-over" : "screen";
        g.imageSmoothingEnabled = true;
        const R = L.ready!.radius, h = L.ready!.h;
        const cw = (s.rx / R) * h, ch = (s.ry / R) * h;
        if ("filter" in g) g.filter = `blur(${(0.45 * Math.min(cw, ch)).toFixed(1)}px)`;
        g.drawImage(L.off, s.cx - (L.m / 2) * cw, s.cy - (L.m / 2) * ch, L.m * cw, L.m * ch);
        g.restore();
        busy = true;
      }
      // an energy halo: brightness follows the piece's kinetic energy on a log scale
      if (fade > 0.01) {
        g.save();
        g.globalCompositeOperation = "screen";
        g.strokeStyle = p.kind === "cymbal" || p.kind === "hat" ? `rgba(255, 196, 92, ${0.55 * fade})` : `rgba(99, 210, 192, ${0.55 * fade})`;
        g.lineWidth = Math.max(1.5, 3 * S) * (0.5 + fade);
        ellipse(g, s.cx, s.cy, s.rx * 1.02 + 6 * S * fade, s.ry * 1.02 + 6 * S * fade);
        g.stroke();
        g.restore();
        busy = true;
      }
      // where the stick landed (fades in 250 ms)
      if (L.lastStrike && now - L.lastStrike.at < 250) {
        const a = 1 - (now - L.lastStrike.at) / 250;
        const R = L.ready ? L.ready.radius : 0.2;
        const x = s.cx + (L.lastStrike.x / R) * s.rx, y = s.cy + (L.lastStrike.y / R) * s.ry;
        g.fillStyle = `rgba(255,255,255,${0.7 * a})`;
        ellipse(g, x, y, 5 * S * (1 + 1.5 * (1 - a)), 3 * S * (1 + 1.5 * (1 - a)));
        g.fill();
        busy = true;
      }
    }
    slowLabel.hidden = !slowTag;
    if ((busy || running) && visible) raf = requestAnimationFrame(draw);
    else raf = 0;
  }
  const kickDraw = (): void => { if (!raf && visible) raf = requestAnimationFrame(draw); };
  // a running kit scrolled out of view is not drawn every frame; it draws again as it comes back
  new IntersectionObserver((entries) => {
    visible = entries.some((e) => e.isIntersecting);
    if (visible) kickDraw();
  }).observe(canvas);

  function paintField(L: Live, v: Float32Array, fade: number, fresh: boolean): void {
    const img = L.img!, cells = L.cells!;
    const paint = L.def.kind === "drum";
    let vmax = 0;
    for (let k = 0; k < v.length; k++) vmax = Math.max(vmax, Math.abs(v[k]!));
    // the colour scale follows the motion down as it decays, so the shape stays visible; in the
    // slow-motion replay each frame is scaled on its own so the front shows as it spreads
    L.vref = fresh ? Math.max(vmax, 1e-9) : Math.max(vmax, L.vref * 0.92, 1e-7);
    const d = img.data;
    d.fill(0);
    for (let k = 0; k < v.length; k++) {
      const u = v[k]! / L.vref;
      const a = Math.min(1, Math.abs(u) * 1.6) * fade;
      const o = cells[k]! * 4;
      if (paint) {
        // ember up, deep verdigris down, opaque enough to read on a white head
        if (u >= 0) { d[o] = 232; d[o + 1] = 96 - 40 * Math.min(1, u); d[o + 2] = 44; }
        else { d[o] = 24; d[o + 1] = 128; d[o + 2] = 120; }
        d[o + 3] = 235 * a;
      } else {
        if (u >= 0) { d[o] = 255; d[o + 1] = 120 + 120 * Math.min(1, u); d[o + 2] = 40 + 160 * Math.max(0, u - 0.6); }
        else { d[o] = 60; d[o + 1] = 200 + 40 * Math.min(1, -u); d[o + 2] = 210 + 45 * Math.min(1, -u); }
        d[o + 3] = 255 * a;
      }
    }
    L.off!.getContext("2d")!.putImageData(img, 0, 0);
  }

  // ------------------------------------------------------------- playing
  function hitAt(id: PieceId, x: number, y: number, speed: number): void {
    const L = lives.get(id)!;
    if (!L.lane) return;
    L.lane.strike(x, y, speed);
    L.force = 0;
    L.lastStrike = { x, y, at: performance.now() };
    kickDraw();
  }
  function hitCentre(id: PieceId): void {
    if (!running) { void begin(); return; }
    const L = lives.get(id)!;
    const R = L.ready?.radius ?? 0.2;
    const at = id === "kick" ? [0, 0.02] : id.startsWith("hihat") || id === "crash" ? [0.85 * R, 0] : id === "ride" ? [0.6 * R, 0] : [0.35 * R, 0.1 * R];
    hitAt(id, at[0]!, at[1]!, L.def.speed * Number(dyn.value));
  }
  function pick(cx: number, cy: number): { id: PieceId; x: number; y: number } | null {
    const order = [...PIECES.filter((p) => p.id !== "hihat-open")].sort((a, b) => spotOf(b.id).z - spotOf(a.id).z);
    for (const p of order) {
      const s = px(spotOf(p.id));
      const u = (cx - s.cx) / s.rx, v = (cy - s.cy) / s.ry;
      if (u * u + v * v <= 1.04) {
        const id: PieceId = p.id === "hihat-closed" && hatOpen ? "hihat-open" : p.id;
        const R = lives.get(id)!.ready?.radius ?? 0.2;
        return { id, x: u * R, y: v * R };
      }
    }
    return null;
  }
  canvas.addEventListener("pointerdown", (e) => {
    if (!running) return;
    const r = canvas.getBoundingClientRect();
    const d = W / r.width;
    const hit = pick((e.clientX - r.left) * d, (e.clientY - r.top) * d);
    if (!hit) return;
    e.preventDefault();
    // a pen or a pressure-sensing screen sets the stroke; otherwise the Stroke slider does
    const pressure = e.pointerType === "pen" || (e.pressure > 0 && e.pressure !== 0.5) ? 0.5 + e.pressure : 1;
    const def = lives.get(hit.id)!.def;
    hitAt(hit.id, hit.x, hit.y, def.speed * Number(dyn.value) * pressure);
  });
  window.addEventListener("keydown", (e) => {
    if (!running || e.metaKey || e.ctrlKey || e.altKey || (e.target as HTMLElement).closest("input, textarea, select")) return;
    const def = PIECES.find((p) => p.key === e.key.toUpperCase());
    if (!def || !root.isConnected) return;
    const r = root.getBoundingClientRect();
    if (r.bottom < 0 || r.top > innerHeight) return;
    e.preventDefault();
    hitCentre(def.id);
  });

  // every start is a session; stop() ends it, and a start still building finds out at its next await
  let session = 0;
  const stop = (): void => {
    session++;
    running = false;
    for (const L of lives.values()) { L.lane?.dispose(); L.lane = null; L.field = null; L.kinetic = 0; L.replay = null; L.failure = null; }
    mix?.dispose();
    mix = null;
    // nothing plays: the audio device rests until the next start resumes it
    if (ctx?.state === "running") void ctx.suspend();
    overlay.hidden = false;
    start.textContent = "▶ Start the kit";
    start.disabled = false;
    hud.textContent = "";
    readout.hidden = true;
  };
  // another player starting stops the kit, also while it is still being built
  const onClaimed = (): void => { if (running || start.disabled) stop(); };
  registerPlayer(onClaimed);

  async function begin(): Promise<void> {
    if (running || start.disabled) return;
    let environment: { record: string; label: string };
    try { environment = readAtmosphere(); }
    catch (e) { atmosphereStatus.textContent = (e as Error).message; return; }
    // Every lane, including grids rebuilt after failed probes, receives this
    // same immutable record even if controls are edited while startup awaits.
    claimPlayback(onClaimed);
    const my = ++session;
    start.disabled = true;
    start.textContent = "Building the kit…";
    atmosphereStatus.textContent = `Building in ${environment.label}…`;
    try {
      ctx ??= createLiveAudioContext();
      await ctx.resume();
      if (my !== session) return;
      const m = await createLiveMix(ctx, ctx.destination);
      if (my !== session) { m.dispose(); return; }
      mix = m;
      let built = 0;
      await Promise.all(PIECES.map(async (def) => {
        const L = lives.get(def.id)!;
        const lane = new KitLane(m, def.id, { gain: def.gain, gasRecord: environment.record, aheadFrames: 1024, minRealtime: 1.7, benchMs: def.benchMs, scales: def.scales, fieldHz: 45 });
        L.lane = lane;
        lane.onField = (f) => {
          L.field = f;
          L.kinetic = f.kinetic;
          L.kpeak = Math.max(L.kpeak * 0.999, f.kinetic);
          L.force = Math.max(L.force, f.force);
          kickDraw();
        };
        lane.onStats = (s) => { if (s.asleep) { L.kinetic = 0; L.kpeak = 0; } };
        lane.onBurst = (b) => { L.replay = { b, at: performance.now(), n: b.count }; L.kpeak = Math.max(L.kpeak, 1e-30); kickDraw(); };
        // an engine that fails once ready has left the mix (mix.ts LaneWorker): its surface goes dark (no
        // more fields arrive to fade it) and the readout says so, also if the others are still building
        lane.onFailure = (e) => {
          if (my !== session) return;
          L.failure = e.message;
          L.field = null; L.kinetic = 0; L.kpeak = 0; L.replay = null;
          if (running) paintHud();
          kickDraw();
        };
        const r = await lane.ready;
        if (my !== session) { lane.dispose(); return; }
        if (r.gasRecord !== environment.record) throw new Error("The kit engine did not retain the requested atmosphere. Reload after rebuilding the live kit.");
        L.ready = r;
        // map every surface node to a cell of a square image the size of the grid
        const m2 = 2 * Math.ceil(r.radius / r.h) + 3;
        L.m = m2;
        L.cells = new Int32Array(r.xy.length / 2);
        for (let k = 0; k < L.cells.length; k++) {
          const i = Math.round(r.xy[2 * k]! / r.h) + (m2 - 1) / 2, j = Math.round(r.xy[2 * k + 1]! / r.h) + (m2 - 1) / 2;
          L.cells[k] = j * m2 + i;
        }
        L.off = document.createElement("canvas");
        L.off.width = L.off.height = m2;
        L.img = new ImageData(m2, m2);
        built += 1;
        if (my === session) start.textContent = `Building the kit… ${built}/${PIECES.length}`;
      }));
      if (my !== session) return;
      running = true;
      const gas = lives.get("kick")?.ready?.info.gas;
      atmosphereStatus.textContent = `${environment.label}. ${gas
        ? `Density ${gas.densityKgM3.toFixed(3)} kg/m³ · sound speed ${gas.soundSpeedMps.toFixed(1)} m/s. ` : ""}All eight pieces use this atmosphere.`;
      overlay.hidden = true;
      kickDraw();
      paintHud();
      energyLoop(my);
    } catch (e) {
      if (my !== session) return; // stopped while building: its lanes were disposed on purpose
      stop();
      hud.textContent = `Could not start: ${(e as Error).message}`;
      atmosphereStatus.textContent = `Atmosphere was not applied: ${(e as Error).message}`;
    }
  }
  start.addEventListener("click", () => void begin());
  atmosphereForm.addEventListener("submit", e => {
    e.preventDefault();
    try {
      const selected = readAtmosphere();
      const active = running || start.disabled;
      if (active) stop();
      atmosphereStatus.textContent = `${selected.label}. ${active ? "Restarting from rest…" : "Press Start the kit to apply."}`;
      if (active) void begin();
    } catch (error) { atmosphereStatus.textContent = (error as Error).message; }
  });

  function paintHud(): void {
    hud.replaceChildren(...PIECES.map((p) => {
      const L = lives.get(p.id)!, r = L.ready!;
      return el("span", { class: "kit-hud-item" }, el("b", {}, p.name),
        L.failure ? ` stopped: ${L.failure}` : ` ${r.info.surfaceNodes} nodes · grid ${r.scale}× · ${r.realtime.toFixed(1)}× real time`);
    }));
  }

  // a compact live readout of the piece played last
  const readout = el("div", { class: "kit-readout", "aria-live": "off" });
  stage.append(readout);
  function energyLoop(my: number): void {
    if (!running || my !== session) return; // one loop per start, even after a quick stop and restart
    let best: Live | null = null;
    for (const L of lives.values()) if (L.lastStrike && (!best || L.lastStrike.at > best.lastStrike!.at)) best = L;
    if (best && best.kinetic > 0) {
      readout.hidden = false;
      readout.replaceChildren(el("b", {}, best.def.name),
        el("span", {}, `energy ${fmtJ(best.kinetic)}`),
        el("span", {}, `stick force ${best.force > 0 ? best.force.toFixed(0) + " N" : "—"}`));
    } else readout.hidden = true;
    setTimeout(() => energyLoop(my), 120);
  }

  resize();
}

function dynName(v: number): string {
  return v < 0.55 ? "pp" : v < 0.75 ? "p" : v < 0.95 ? "mp" : v < 1.15 ? "mf" : v < 1.35 ? "f" : "ff";
}
function fmtJ(j: number): string {
  return j >= 1 ? `${j.toFixed(2)} J` : j >= 1e-3 ? `${(j * 1e3).toFixed(1)} mJ` : j >= 1e-6 ? `${(j * 1e6).toFixed(1)} µJ` : `${(j * 1e9).toPrecision(2)} nJ`;
}

// ---------------------------------------------------------------- painters
type Box = { cx: number; cy: number; rx: number; ry: number };
function ellipse(g: CanvasRenderingContext2D, cx: number, cy: number, rx: number, ry: number): void {
  g.beginPath();
  g.ellipse(cx, cy, Math.max(0.5, rx), Math.max(0.5, ry), 0, 0, Math.PI * 2);
}

function shell(g: CanvasRenderingContext2D, s: Box, depth: number, finish: string, chromeShell: boolean): void {
  // the shell's side: a cylinder seen from above and in front
  const grad = g.createLinearGradient(s.cx - s.rx, 0, s.cx + s.rx, 0);
  if (chromeShell) {
    grad.addColorStop(0, "#5b6068"); grad.addColorStop(0.18, "#e9edf2"); grad.addColorStop(0.35, "#8c929b");
    grad.addColorStop(0.62, "#f4f6f8"); grad.addColorStop(0.85, "#6d737c"); grad.addColorStop(1, "#3b3f45");
  } else {
    grad.addColorStop(0, shade(finish, -0.55)); grad.addColorStop(0.22, shade(finish, 0.25)); grad.addColorStop(0.5, finish);
    grad.addColorStop(0.7, shade(finish, 0.35)); grad.addColorStop(1, shade(finish, -0.6));
  }
  g.fillStyle = grad;
  g.beginPath();
  g.ellipse(s.cx, s.cy, s.rx, s.ry, 0, 0, Math.PI);
  g.lineTo(s.cx - s.rx, s.cy + depth);
  g.ellipse(s.cx, s.cy + depth, s.rx, s.ry, 0, Math.PI, 0, true);
  g.closePath();
  g.fill();
  // lugs
  for (let k = 0; k < 8; k++) {
    const t = Math.PI * (0.06 + (0.88 * k) / 7);
    const x = s.cx + Math.cos(t) * s.rx * 0.98, y = s.cy + Math.sin(t) * s.ry + depth * 0.32;
    const w = Math.max(2, s.rx * 0.045 * Math.abs(Math.sin(t)) + 1.5), hgt = depth * 0.36;
    const lg = g.createLinearGradient(x - w, 0, x + w, 0);
    lg.addColorStop(0, "#4a4f57"); lg.addColorStop(0.5, "#f1f4f7"); lg.addColorStop(1, "#5b6068");
    g.fillStyle = lg;
    g.beginPath();
    g.roundRect(x - w / 2, y, w, hgt, w / 2);
    g.fill();
  }
  // bottom hoop
  g.strokeStyle = "rgba(220, 225, 232, 0.65)";
  g.lineWidth = Math.max(1, s.ry * 0.05);
  g.beginPath();
  g.ellipse(s.cx, s.cy + depth, s.rx, s.ry, 0, 0.05, Math.PI - 0.05);
  g.stroke();
}

function head(g: CanvasRenderingContext2D, s: Box): void {
  // coated white head with a soft spot of light and a darker edge where it meets the hoop
  const hg = g.createRadialGradient(s.cx - s.rx * 0.25, s.cy - s.ry * 0.35, 1, s.cx, s.cy, Math.max(s.rx, s.ry));
  hg.addColorStop(0, "#fbf8f1"); hg.addColorStop(0.6, "#e9e3d6"); hg.addColorStop(1, "#bdb5a6");
  g.fillStyle = hg;
  ellipse(g, s.cx, s.cy, s.rx, s.ry);
  g.fill();
  // chrome hoop
  const ring = g.createLinearGradient(s.cx - s.rx, s.cy - s.ry, s.cx + s.rx, s.cy + s.ry);
  ring.addColorStop(0, "#f5f7fa"); ring.addColorStop(0.3, "#7c828b"); ring.addColorStop(0.55, "#ffffff"); ring.addColorStop(0.8, "#60666e"); ring.addColorStop(1, "#d9dde3");
  g.strokeStyle = ring;
  g.lineWidth = Math.max(2, s.rx * 0.055);
  ellipse(g, s.cx, s.cy, s.rx, s.ry);
  g.stroke();
}

function drum(g: CanvasRenderingContext2D, s: Box, finish: string, snare: boolean): void {
  const depth = s.ry * (snare ? 0.75 : 1.25);
  g.fillStyle = "rgba(0,0,0,0.45)";
  ellipse(g, s.cx + s.rx * 0.05, s.cy + depth + s.ry * 0.25, s.rx * 1.02, s.ry * 0.7);
  g.fill();
  shell(g, s, depth, finish, snare);
  head(g, s);
}

function kick(g: CanvasRenderingContext2D, s: Box, finish: string): void {
  // seen from the front: the resonant head with a port-less logo, hoops and spurs
  g.fillStyle = "rgba(0,0,0,0.5)";
  ellipse(g, s.cx, s.cy + s.ry * 0.96, s.rx * 1.05, s.ry * 0.12);
  g.fill();
  const rim = g.createRadialGradient(s.cx, s.cy, s.rx * 0.9, s.cx, s.cy, s.rx * 1.08);
  rim.addColorStop(0, shade(finish, -0.2)); rim.addColorStop(0.5, finish); rim.addColorStop(1, shade(finish, -0.7));
  g.fillStyle = rim;
  ellipse(g, s.cx, s.cy, s.rx * 1.07, s.ry * 1.07);
  g.fill();
  const hg = g.createRadialGradient(s.cx - s.rx * 0.3, s.cy - s.ry * 0.35, 2, s.cx, s.cy, s.rx);
  hg.addColorStop(0, "#2a2a30"); hg.addColorStop(1, "#111116");
  g.fillStyle = hg;
  ellipse(g, s.cx, s.cy, s.rx, s.ry);
  g.fill();
  g.font = `600 ${Math.round(s.rx * 0.16)}px Fraunces Variable, Georgia, serif`;
  g.textAlign = "center";
  g.textBaseline = "middle";
  g.fillStyle = "rgba(241, 234, 219, 0.82)";
  g.fillText("FrankenSim", s.cx, s.cy + s.ry * 0.1);
  g.font = `500 ${Math.round(s.rx * 0.07)}px JetBrains Mono, monospace`;
  g.fillStyle = "rgba(217, 169, 79, 0.85)";
  g.fillText("18 × 18", s.cx, s.cy + s.ry * 0.32);
  const ring = g.createLinearGradient(s.cx - s.rx, s.cy - s.ry, s.cx + s.rx, s.cy + s.ry);
  ring.addColorStop(0, "#f5f7fa"); ring.addColorStop(0.4, "#6f757e"); ring.addColorStop(0.6, "#ffffff"); ring.addColorStop(1, "#555b63");
  g.strokeStyle = ring;
  g.lineWidth = Math.max(2, s.rx * 0.035);
  ellipse(g, s.cx, s.cy, s.rx * 1.075, s.ry * 1.075);
  g.stroke();
  // spurs
  g.strokeStyle = "#c9cdd3";
  g.lineWidth = Math.max(2, s.rx * 0.025);
  for (const sg of [-1, 1]) {
    g.beginPath();
    g.moveTo(s.cx + sg * s.rx * 0.8, s.cy + s.ry * 0.6);
    g.lineTo(s.cx + sg * s.rx * 1.25, s.cy + s.ry * 1.0);
    g.stroke();
  }
}

function stand(g: CanvasRenderingContext2D, s: Box, hat: boolean): void {
  g.strokeStyle = "#8d939b";
  g.lineWidth = Math.max(2, s.rx * 0.03);
  const foot = s.cy + s.ry + s.rx * (hat ? 1.35 : 1.15);
  g.beginPath();
  g.moveTo(s.cx, s.cy);
  g.lineTo(s.cx, foot);
  g.stroke();
  g.lineWidth = Math.max(1.5, s.rx * 0.02);
  g.beginPath();
  for (const sg of [-1, 0, 1]) { g.moveTo(s.cx, foot - s.rx * 0.18); g.lineTo(s.cx + sg * s.rx * 0.38, foot + (sg === 0 ? s.rx * 0.05 : 0)); }
  g.stroke();
}

function cymbal(g: CanvasRenderingContext2D, s: Box, hat: boolean): void {
  g.fillStyle = "rgba(0,0,0,0.35)";
  ellipse(g, s.cx + s.rx * 0.04, s.cy + s.ry * 0.5, s.rx, s.ry);
  g.fill();
  if (hat) {
    // the bottom shell's edge peeking out under the top one
    g.fillStyle = "#6e5326";
    ellipse(g, s.cx, s.cy + s.ry * 0.16, s.rx * 0.995, s.ry);
    g.fill();
  }
  const br = g.createLinearGradient(s.cx - s.rx, s.cy - s.ry, s.cx + s.rx, s.cy + s.ry);
  br.addColorStop(0, "#8a6a2e"); br.addColorStop(0.28, "#e8c27a"); br.addColorStop(0.45, "#b58a3e");
  br.addColorStop(0.62, "#f3d796"); br.addColorStop(0.8, "#9c7533"); br.addColorStop(1, "#5e4520");
  g.fillStyle = br;
  ellipse(g, s.cx, s.cy, s.rx, s.ry);
  g.fill();
  // lathing: fine concentric grooves
  g.lineWidth = 1;
  for (let k = 1; k < 26; k++) {
    const f = k / 26;
    g.strokeStyle = k % 2 ? "rgba(70, 45, 12, 0.18)" : "rgba(255, 238, 190, 0.14)";
    ellipse(g, s.cx, s.cy, s.rx * f, s.ry * f);
    g.stroke();
  }
  // bell
  const bell = g.createRadialGradient(s.cx - s.rx * 0.05, s.cy - s.ry * 0.25, 1, s.cx, s.cy, s.rx * 0.2);
  bell.addColorStop(0, "#fff0c4"); bell.addColorStop(0.6, "#c99a48"); bell.addColorStop(1, "#7d5c26");
  g.fillStyle = bell;
  ellipse(g, s.cx, s.cy - s.ry * 0.08, s.rx * 0.2, s.ry * 0.32);
  g.fill();
  g.fillStyle = "#3a3d42";
  ellipse(g, s.cx, s.cy - s.ry * 0.14, s.rx * 0.03, s.ry * 0.07);
  g.fill();
}

function shade(hex: string, k: number): string {
  const n = parseInt(hex.slice(1), 16);
  const f = (c: number): number => Math.round(k >= 0 ? c + (255 - c) * k : c * (1 + k));
  return `rgb(${f(n >> 16)}, ${f((n >> 8) & 255)}, ${f(n & 255)})`;
}
