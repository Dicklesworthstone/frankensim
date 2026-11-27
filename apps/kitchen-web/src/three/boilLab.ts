/**
 * Put the water on (chapter 02, second figure): the same water four ways on two six-burner ranges. One 40-quart pot
 * without a lid and one with, on the left range; the water split between two 20-quart pots on two burners, without
 * lids and with, on the right. Each pot heats as sim/potBoil.ts computes it: its water rises in temperature, steams
 * off the open pots near the boil, rolls at the boil, and the lids rattle.
 *
 * Beside the stage, each way's time to the boil, when to light it for 6:30, and the gas and water it costs; under it,
 * each pot's temperature over time, drawn as the stage plays.
 */
import * as THREE from "three";
import { BOIL_C, BURNER_W, GAS_INTO_POT, MAX_LITRES, SETUPS, TAP_C, boilAll, lossesW, type Setup, type SetupRun } from "../sim/potBoil";
import { bindCanvasResize, clipFitter, createScrollSafeOrbit, createStudioRenderer, observeInView, shaderGate } from "./kit";
import { BURNER_X, BURNER_Z, GRATE_Y, RANGE } from "./morningKitchen";
import { buildFixture, type Model, type Piece } from "./models/all";
import { makeStaff } from "./models/people";
import { studioEnvironment, studioLights } from "./models/stage";
import { backdrop, place, tilePlane } from "./weekLab";

const POT_40 = { id: "NSF-STOCKPOT-40QT", w: 0.38, d: 0.46, h: 0.375 };
const POT_20 = { id: "NSF-STOCKPOT-20QT", w: 0.31, d: 0.38, h: 0.295 };
/** The two ranges, side by side, their fronts at z = 0; a hair apart, so no side panel lies on the other's. */
const GAP = 0.004;
const RANGE_X0 = [-RANGE.w - GAP / 2, GAP / 2] as const;

/** Where each way's pots stand: range, burner (0..2 front row left to right, 3..5 back row), and whether turned so
 *  its handles point across. The 20-quart pots in front keep theirs fore and aft, clear of the big pot on the next
 *  range; the ones behind are turned, or their front handles would meet the back handles of the ones in front. */
const SPOTS: Record<Setup["key"], { range: 0 | 1; burner: number; turned: boolean }[]> = {
  "one-open": [{ range: 0, burner: 0, turned: false }],
  "one-lid": [{ range: 0, burner: 2, turned: false }],
  "two-open": [
    { range: 1, burner: 0, turned: false },
    { range: 1, burner: 2, turned: false },
  ],
  "two-lid": [
    { range: 1, burner: 3, turned: true },
    { range: 1, burner: 5, turned: true },
  ],
};
/** What the steam drawn over a pot scales by: watts going off as steam that the model draws as a hard boil. */
const STEAM_W = 1500;
/** Each way's colour, on the stage's labels, the tiles and the chart. */
export const SETUP_COLOR: Record<Setup["key"], string> = { "one-open": "#f97316", "one-lid": "#facc15", "two-open": "#22d3ee", "two-lid": "#a5b4fc" };

/** The ranges and the pots, for main.ts to build ahead. */
export const BOIL_PIECES: readonly Piece[] = [
  [RANGE.id, RANGE, 2],
  [RANGE.id, RANGE, 2],
  [POT_40.id, POT_40, 2],
  [POT_40.id, POT_40, 2],
  ...Array.from({ length: 4 }, (): Piece => [POT_20.id, POT_20, 2]),
];

/** The flame settings offered: the valve's opening, which the gas burned follows. */
const FLAMES = [
  { key: "full", name: "Full", v: 1 },
  { key: "three", name: "Three-quarters", v: 0.75 },
  { key: "half", name: "Half", v: 0.5 },
] as const;
/** When the water is wanted at the boil: the macaroni goes in at 6:30 (minutes after midnight). */
const NEEDED = 6 * 60 + 30;
/** How long to wait on a pot that has not boiled (minutes). */
const MAX_MIN = 180;
/** The loop: the burners lighting, the heating, and a hold on the end (seconds). */
const LIGHT_S = 1.4;
const HEAT_S = 30;
const END_S = 4;

const clock = (m: number) => {
  const r = Math.round(m);
  const h = Math.floor(r / 60);
  return `${h}:${String(r - h * 60).padStart(2, "0")}`;
};
/** "1 h 05" or "25 min". */
const span = (m: number) => {
  const r = Math.round(m);
  return r >= 60 ? `${Math.floor(r / 60)} h${r % 60 ? ` ${String(r % 60).padStart(2, "0")}` : ""}` : `${r} min`;
};
/** Minutes to the boil, whole minutes up: the water is at the boil by then. */
const upTo = (m: number) => span(Math.ceil(m - 1e-6));
/** A temperature as every readout shows it: whole degrees, tenths (down) near the boil. */
const degC = (c: number) => (c >= 95 ? (Math.floor(c * 10 + 1e-6) / 10).toFixed(1) : String(Math.floor(c + 1e-6)));
const smooth = (x: number) => {
  const t = Math.min(1, Math.max(0, x));
  return t * t * (3 - 2 * t);
};

/** A run's temperature and litres at minute `t` (straight between samples; at the boil from its minute on). */
function at(r: SetupRun, t: number): { c: number; litres: number; boiled: boolean } {
  if (r.boilMin !== null && t >= r.boilMin) return { c: BOIL_C, litres: r.litres[r.litres.length - 1]!, boiled: true };
  const n = r.t.length;
  if (t >= r.t[n - 1]!) return { c: r.c[n - 1]!, litres: r.litres[n - 1]!, boiled: false };
  let lo = 0;
  let hi = n - 1;
  while (hi - lo > 1) {
    const m = (lo + hi) >> 1;
    if (r.t[m]! <= t) lo = m;
    else hi = m;
  }
  const f = (t - r.t[lo]!) / Math.max(1e-9, r.t[hi]! - r.t[lo]!);
  return { c: r.c[lo]! + f * (r.c[hi]! - r.c[lo]!), litres: r.litres[lo]! + f * (r.litres[hi]! - r.litres[lo]!), boiled: false };
}

export function mountBoilLab(root: HTMLElement): void {
  const mountT0 = performance.now();
  const stageEl = root.querySelector<HTMLElement>("[data-bo-stage]");
  const canvas = root.querySelector<HTMLCanvasElement>("[data-bo-canvas]");
  if (!stageEl || !canvas) return;
  const maybeRenderer = createStudioRenderer(canvas, false);
  if (!maybeRenderer) return;
  const renderer = maybeRenderer;
  renderer.toneMappingExposure = 1.0;
  const phone = window.matchMedia("(max-width: 720px)").matches;

  /* The scene: two ranges against the back, the pots on them, a cook beside them. */
  const scene = new THREE.Scene();
  scene.background = backdrop();
  scene.fog = new THREE.Fog(0x151510, 8, 16);
  scene.environment = studioEnvironment(renderer);
  scene.environmentIntensity = 0.5;
  const lights = studioLights(scene, new THREE.Vector3(0, 1.0, -0.3), 3.6, phone ? 1024 : 2048);
  lights.key.position.set(-1.6, 5.6, 3.8);
  lights.key.intensity = 1.5;
  lights.key.shadow.normalBias = 0.02;
  lights.rim.position.set(3.0, 3.2, -2.6);
  lights.rim.intensity = 0.75;
  lights.hemi.intensity = 0.5;
  {
    const floor = tilePlane(10, 7);
    floor.rotation.x = -Math.PI / 2;
    floor.position.set(0.3, 0, 1.4);
    scene.add(floor);
  }
  const ranges = RANGE_X0.map((x0) => {
    const m = buildFixture(RANGE.id, RANGE, 2);
    scene.add(place(m, x0, 0, true));
    const flames = m.hooks.flames;
    flames?.setAll(0);
    return { m, flames, flow: [0, 0, 0, 0, 0, 0], x0 };
  });
  /** A burner's centre on a range. */
  const burnerAt = (r: number, b: number) => new THREE.Vector3(RANGE_X0[r]! + BURNER_X[b % 3]!, GRATE_Y, -BURNER_Z[b >= 3 ? 1 : 0]!);

  interface Pot {
    setup: Setup;
    model: Model;
    ext: typeof POT_40;
    range: number;
    burner: number;
    /** The rim's centre, for the labels. */
    top: THREE.Vector3;
    phase: number;
    level: number;
    boil: number;
    steam: number;
  }
  const pots: Pot[] = [];
  for (const setup of SETUPS) {
    const ext = setup.pot.key === "40qt" ? POT_40 : POT_20;
    for (const s of SPOTS[setup.key]) {
      const model = buildFixture(ext.id, ext, 2);
      const b = burnerAt(s.range, s.burner);
      // The fixture's frame has the pot's axis at (w/2, 0, -d/2): set it on the burner, turned if need be.
      const holder = new THREE.Group();
      model.group.position.set(-ext.w / 2, 0, ext.d / 2);
      holder.add(model.group);
      holder.position.set(b.x, GRATE_Y + 0.0006, b.z);
      if (s.turned) holder.rotation.y = Math.PI / 2;
      scene.add(holder);
      model.hooks.set?.("water", 1);
      if (!setup.lid) {
        // No lid: put away, and the steam rises freely.
        const lid = model.joints.get("pot_lid_lift");
        if (lid) lid.visible = false;
        model.hooks.open?.(1);
      }
      pots.push({ setup, model, ext, range: s.range, burner: s.burner, top: new THREE.Vector3(b.x, GRATE_Y + ext.h, b.z), phase: pots.length * 1.7, level: -1, boil: -1, steam: -1 });
    }
  }
  // The cook: beside the right-hand range, an eye on the pots.
  const cook = makeStaff("p50", 2, { seed: 41, apron: 0x26344d, sleeves: "cuffed" });
  {
    const x = RANGE_X0[1] + RANGE.w + 0.42;
    const z = 0.32;
    const look = new THREE.Vector3(RANGE_X0[1] + RANGE.w * 0.55, 0, -0.35);
    cook.place(x, z, Math.atan2(look.x - x, look.z - z));
  }
  scene.add(cook.group);

  /* The camera: from in front, raised, both ranges and the cook in view, looking into the pots. */
  const camera = new THREE.PerspectiveCamera(32, 16 / 9, 0.05, 40);
  const TARGET = new THREE.Vector3(0.15, 1.0, -0.4);
  const VIEW = new THREE.Vector3();
  const keep: THREE.Vector3[] = [];
  const aim = () => {
    const wide = stageEl.clientWidth >= 1.3 * stageEl.clientHeight;
    VIEW.set(0.12, wide ? 0.78 : 0.95, 1).normalize();
    const L = RANGE_X0[0];
    const R = RANGE_X0[1] + RANGE.w;
    keep.length = 0;
    // The ranges' fronts at the floor and at the grates, the pots' rims, the back shelf, the cook's head.
    keep.push(new THREE.Vector3(L, 0.02, 0), new THREE.Vector3(R, 0.02, 0), new THREE.Vector3(L, RANGE.h, -RANGE.d), new THREE.Vector3(R, RANGE.h, -RANGE.d));
    for (const p of pots) for (const s of [-1, 1]) keep.push(p.top.clone().add(new THREE.Vector3((s * p.ext.w) / 2, 0.08, 0)));
    keep.push(new THREE.Vector3(R + 0.42, 1.85, 0.32), new THREE.Vector3(R + 0.62, 0, 0.32));
  };
  const controls = createScrollSafeOrbit(camera, canvas, [0, 1, 0], { minDistance: 1.0, maxDistance: 12, maxPolarAngle: Math.PI * 0.47 });
  let userMoved = false;
  controls.addEventListener("start", () => {
    userMoved = true;
  });
  const hudEl = stageEl.querySelector<HTMLElement>(".bo-hud");
  /** The HUD's bottom (px from the stage's top) with the longest words it ever shows. */
  const hudRoom = (): number => {
    if (!hudEl) return 0;
    const parts = ["[data-bo-hud-clock]", "[data-bo-hud-doing]", "[data-bo-hud-note]"].map((k) => hudEl.querySelector<HTMLElement>(k));
    const was = parts.map((e) => e?.textContent ?? "");
    const doing = [`Burners lit, three-quarters flame`, `${MAX_LITRES} litres heating from ${TAP_C} °C`, "All four at the boil"];
    for (const s of SETUPS) doing.push(`${s.name}: at the boil`, `${s.name}: simmering at 99.9 °C`);
    let most = 0;
    for (const d of doing) {
      if (parts[0]) parts[0].textContent = "2 h 59";
      if (parts[1]) parts[1].textContent = d;
      if (parts[2]) parts[2].textContent = "3 of 4 boiling";
      most = Math.max(most, hudEl.offsetTop + hudEl.offsetHeight);
    }
    parts.forEach((e, i) => {
      if (e) e.textContent = was[i]!;
    });
    return most;
  };
  /** Room over the pots for their labels, and one more stacked over another (px). */
  const LABEL_ROOM = 26;
  const ndc = new THREE.Vector3();
  let fitW = 0;
  let fitH = 0;
  const frame = () => {
    if (userMoved) return;
    aim();
    const w = Math.max(1, stageEl.clientWidth);
    const h = Math.max(1, stageEl.clientHeight);
    fitW = w;
    fitH = h;
    camera.aspect = Math.max(0.2, w / h);
    camera.clearViewOffset();
    camera.updateProjectionMatrix();
    // Everything between the HUD (and the labels' room under it) and the bottom edge, centred across.
    const top = Math.min(0.88, 1 - (2 * (hudRoom() + 8 + LABEL_ROOM)) / h);
    const bot = -0.94;
    let dist = 2;
    let lo = 0;
    let hi = 0;
    let left = 0;
    let right = 0;
    for (; dist < 16; dist += 0.04) {
      camera.position.copy(TARGET).addScaledVector(VIEW, dist);
      camera.lookAt(TARGET);
      camera.updateMatrixWorld(true);
      lo = left = Infinity;
      hi = right = -Infinity;
      for (const p of keep) {
        ndc.copy(p).project(camera);
        lo = Math.min(lo, ndc.y);
        hi = Math.max(hi, ndc.y);
        left = Math.min(left, ndc.x);
        right = Math.max(right, ndc.x);
      }
      if (right - left < 1.88 && hi - lo < top - bot) break;
    }
    const shiftY = (top + bot) / 2 - (hi + lo) / 2;
    const shiftX = (left + right) / 2;
    if (Math.abs(shiftY) > 0.005 || Math.abs(shiftX) > 0.005) {
      camera.setViewOffset(w, h, (shiftX * w) / 2, (shiftY * h) / 2, w, h);
      camera.updateProjectionMatrix();
    }
    const fog = scene.fog as THREE.Fog;
    fog.near = dist * 1.4;
    fog.far = dist * 3.2;
    controls.target.copy(TARGET);
    controls.update();
  };
  const gate = shaderGate(renderer, camera);
  bindCanvasResize(stageEl, renderer, camera, gate);
  frame();
  new ResizeObserver(() => {
    if (stageEl.clientWidth === fitW && stageEl.clientHeight === fitH) return;
    tagSize.clear();
    frame();
  }).observe(stageEl);
  // The HUD's room is measured in its font: once that has loaded, again.
  void document.fonts?.ready.then(() => frame());

  /* ---------------------------------------------------------------------------------------------------------- */
  /* The page.                                                                                                   */
  /* ---------------------------------------------------------------------------------------------------------- */
  const q = <T extends Element = HTMLElement>(s: string) => root.querySelector<T>(s);
  const litresIn = q<HTMLInputElement>("[data-bo-litres]");
  const litresVal = q("[data-bo-litres-val]");
  const flameBtns = Array.from(root.querySelectorAll<HTMLButtonElement>("[data-bo-flame]"));
  const statusEl = q("[data-bo-status]");
  const chartEl = q("[data-bo-chart]");
  const hudClock = q("[data-bo-hud-clock]");
  const hudDoing = q("[data-bo-hud-doing]");
  const hudNote = q("[data-bo-hud-note]");
  const tagEls = new Map(SETUPS.map((s) => [s.key, q<HTMLElement>(`[data-bo-tag="${s.key}"]`)] as const));
  const tileEls = new Map(SETUPS.map((s) => [s.key, { v: q(`[data-bo-tile="${s.key}"]`), s: q(`[data-bo-tile-s="${s.key}"]`) }] as const));

  /** The loop's clock (seconds), a minute held by the hook, and the minute on show. */
  let loopT = 0;
  let frozen: number | null = null;
  let playMin = 0;
  let litres = Math.min(MAX_LITRES, Number(litresIn?.value ?? 30));
  let flame: (typeof FLAMES)[number] = FLAMES[0];
  let runs: SetupRun[] = boilAll(litres, flame.v, MAX_MIN);
  const runOf = (k: Setup["key"]) => runs.find((r) => r.setup.key === k)!;
  /** The minutes the stage plays to: a little past the last pot to boil, or, if one never does, all the time given. */
  let endMin = 0;
  const endOf = () => (runs.every((r) => r.boilMin !== null) ? Math.max(...runs.map((r) => r.boilMin!)) + 2 : MAX_MIN);

  function renderControls(): void {
    if (litresIn) litresIn.value = String(litres);
    if (litresVal) litresVal.textContent = `${litres} litres`;
    for (const b of flameBtns) b.setAttribute("aria-pressed", String(b.dataset.boFlame === flame.key));
  }
  function renderTiles(): void {
    for (const r of runs) {
      const t = tileEls.get(r.setup.key);
      if (!t) continue;
      const lines: string[] = [];
      const burners = r.setup.pots === 2 ? "2 burners" : "1 burner";
      if (r.boilMin === null) {
        if (t.v) t.v.textContent = `Not in ${span(MAX_MIN)}`;
        lines.push(`levels off at ${degC(r.c[r.c.length - 1]!)} °C`, `${burners} × ${span(MAX_MIN)}`);
      } else {
        if (t.v) t.v.textContent = upTo(r.boilMin);
        lines.push(`light it by ${clock(NEEDED - Math.ceil(r.boilMin - 1e-6))}`, `${burners} × ${upTo(r.boilMin)}`);
      }
      const steamed = r.steamL * r.setup.pots;
      lines.push(steamed >= 0.05 ? `${steamed.toFixed(1)} L to steam` : "under 0.1 L to steam");
      if (t.s) t.s.textContent = lines.join("\n");
    }
  }
  function renderStatus(): void {
    if (!statusEl) return;
    const [oo, ol, tl] = (["one-open", "one-lid", "two-lid"] as const).map(runOf) as [SetupRun, SetupRun, SetupRun];
    const parts: string[] = [];
    if (tl.boilMin !== null && ol.boilMin !== null) {
      const x = ol.boilMin / tl.boilMin;
      const gas = tl.burnerMin / ol.burnerMin - 1;
      const how = x > 1.85 && x < 2.15 ? `half the ${upTo(ol.boilMin)} of one big pot with its lid` : `${x.toFixed(1)} times as fast as one big pot with its lid (${upTo(ol.boilMin)})`;
      const gasWords = Math.abs(gas) < 0.05 ? "for the same gas" : gas > 0 ? `for ${Math.round(gas * 100)}% more gas` : `for ${Math.round(-gas * 100)}% less gas`;
      parts.push(`Split between two pots on two burners, ${litres} litres boil in ${upTo(tl.boilMin)} with lids on: ${how}, ${gasWords}, as twice the heat goes into the same water.`);
    }
    const kept = Math.max(0, oo.steamL - ol.steamL);
    if (oo.boilMin !== null && ol.boilMin !== null) {
      const saves = Math.ceil(oo.boilMin - 1e-6) - Math.ceil(ol.boilMin - 1e-6);
      parts.push(`On the big pot the lid saves ${saves >= 1 ? span(saves) : "under a minute"} and ${kept.toFixed(1)} L of water: near the boil, most of the heat an open pot loses goes off as steam.`);
    } else if (oo.boilMin === null) {
      parts.push(`With no lid the big pot does not boil in ${span(MAX_MIN)}: it levels off at ${degC(oo.c[oo.c.length - 1]!)} °C, a simmer, with the steam off its open top carrying away most of the burner's heat and ${oo.steamL.toFixed(1)} L of the water.`);
    }
    if (tl.boilMin !== null) parts.push(`For the macaroni at ${clock(NEEDED)}, light the two pots with lids by ${clock(NEEDED - Math.ceil(tl.boilMin - 1e-6))}${ol.boilMin !== null ? `, or the big pot with its lid by ${clock(NEEDED - Math.ceil(ol.boilMin - 1e-6))}` : ""}.`);
    statusEl.textContent = parts.join(" ");
    statusEl.dataset.state = oo.boilMin === null ? "warn" : "ok";
  }

  /* The chart: each pot's temperature over the minutes, drawn as the stage plays. */
  let playX: ((m: number) => number) | null = null;
  let playLine: SVGLineElement | null = null;
  let playClip: SVGRectElement | null = null;
  let chartL = 0;
  function renderChart(): void {
    if (!chartEl) return;
    const W = Math.max(280, Math.round(chartEl.clientWidth));
    const narrow = W < 520;
    const H = narrow ? 190 : 220;
    const L = narrow ? 34 : 46;
    const R = W - (narrow ? 10 : 16);
    const T = 14;
    const B = H - 24;
    const c0 = Math.min(TAP_C, 15) - 3;
    const c1 = BOIL_C + 12;
    const x = (m: number) => L + (m / endMin) * (R - L);
    const y = (c: number) => B - ((c - c0) / (c1 - c0)) * (B - T);
    chartL = L;
    let svg = `<svg viewBox="0 0 ${W} ${H}" width="${W}" height="${H}" class="rs-svg bo-svg"><defs><clipPath id="bo-clip"><rect x="0" y="0" width="${L}" height="${H}" data-bo-clip /></clipPath></defs>`;
    for (const c of [20, 40, 60, 80, 100]) {
      svg += `<line class="grid${c === BOIL_C ? " boil" : ""}" x1="${L}" x2="${R}" y1="${y(c).toFixed(1)}" y2="${y(c).toFixed(1)}" />`;
      svg += `<text class="ya" x="${L - 6}" y="${(y(c) + 4).toFixed(1)}" text-anchor="end">${c}°</text>`;
    }
    const step = endMin > 120 ? 30 : endMin > 50 ? 15 : endMin > 24 ? 10 : 5;
    const every = x(step) - x(0) < 44 ? 2 : 1;
    for (let m = 0; m <= endMin + 1e-9; m += step) {
      if (m > 0 && (m / step) % every !== 0) continue;
      svg += `<text class="xa${m === 0 ? " st" : x(m) > R - 16 ? " en" : ""}" x="${x(m).toFixed(1)}" y="${H - 6}">${m === 0 ? "lit" : span(m)}</text>`;
    }
    const pathOf = (r: SetupRun) => {
      const pts: string[] = [];
      for (let i = 0; i < r.t.length && r.t[i]! <= endMin + 1e-9; i++) pts.push(`${x(r.t[i]!).toFixed(1)},${y(r.c[i]!).toFixed(1)}`);
      if (r.boilMin !== null && r.boilMin < endMin) pts.push(`${x(endMin).toFixed(1)},${y(BOIL_C).toFixed(1)}`);
      return `M${pts.join("L")}`;
    };
    for (const r of runs) svg += `<path class="ln faint" d="${pathOf(r)}" style="--c:${SETUP_COLOR[r.setup.key]}" />`;
    svg += `<g clip-path="url(#bo-clip)">`;
    for (const r of runs) svg += `<path class="ln" d="${pathOf(r)}" style="--c:${SETUP_COLOR[r.setup.key]}" />`;
    svg += `</g>`;
    // Where each boils, with its minutes where there is room; a pot that never does, its last temperature.
    const marks = runs
      .filter((r) => r.boilMin !== null && r.boilMin <= endMin)
      .map((r) => ({ r, px: x(r.boilMin!) }))
      .sort((a, b) => a.px - b.px);
    // Over the line, a second row up where two would touch.
    const tier = [-Infinity, -Infinity];
    for (const { r, px } of marks) {
      svg += `<circle class="bp" cx="${px.toFixed(1)}" cy="${y(BOIL_C).toFixed(1)}" r="3.5" style="--c:${SETUP_COLOR[r.setup.key]}" />`;
      const k = px - tier[0]! >= 50 ? 0 : px - tier[1]! >= 50 ? 1 : -1;
      if (k < 0) continue;
      tier[k] = px;
      const tx = Math.min(R - 22, Math.max(L + 22, px));
      svg += `<text class="bpl" x="${tx.toFixed(1)}" y="${(y(BOIL_C) - 8 - 12 * k).toFixed(1)}" text-anchor="middle" style="--c:${SETUP_COLOR[r.setup.key]}">${upTo(r.boilMin!)}</text>`;
    }
    for (const r of runs)
      if (r.boilMin === null) {
        const last = at(r, endMin);
        svg += `<text class="bpl" x="${(R - 2).toFixed(1)}" y="${(y(last.c) + 16).toFixed(1)}" text-anchor="end" style="--c:${SETUP_COLOR[r.setup.key]}">${degC(last.c)} °C</text>`;
      }
    svg += `<line class="play" x1="${L}" x2="${L}" y1="${T - 4}" y2="${B}" />`;
    svg += `</svg>`;
    chartEl.innerHTML = svg;
    playX = x;
    playLine = chartEl.querySelector<SVGLineElement>("line.play");
    playClip = chartEl.querySelector<SVGRectElement>("[data-bo-clip]");
  }

  function apply(): void {
    runs = boilAll(litres, flame.v, MAX_MIN);
    endMin = endOf();
    loopT = 0;
    renderControls();
    renderTiles();
    renderStatus();
    renderChart();
    lastHud = "";
    for (const t of tagEls.values()) if (t) t.dataset.last = "";
  }
  let pending = 0;
  litresIn?.addEventListener("input", () => {
    litres = Math.min(MAX_LITRES, Number(litresIn.value));
    renderControls();
    window.clearTimeout(pending);
    pending = window.setTimeout(apply, 120);
  });
  for (const b of flameBtns)
    b.addEventListener("click", () => {
      const f = FLAMES.find((x) => x.key === b.dataset.boFlame);
      if (!f || f === flame) return;
      flame = f;
      apply();
    });

  /* The readouts. */
  let lastHud = "";
  function showHud(m: number, lighting: boolean): void {
    const boiled = runs.filter((r) => r.boilMin !== null && r.boilMin <= m + 1e-9).sort((a, b) => b.boilMin! - a.boilMin!);
    // Waiting on: the hottest pot not yet boiling.
    const waiting = runs
      .filter((r) => !(r.boilMin !== null && r.boilMin <= m + 1e-9))
      .map((r) => ({ r, c: at(r, m).c }))
      .sort((a, b) => b.c - a.c)[0];
    const levelled = waiting && waiting.r.boilMin === null && waiting.c > waiting.r.c[waiting.r.c.length - 1]! - 0.3;
    const doing = lighting
      ? `Burners lit, ${flame.name.toLowerCase()} flame`
      : boiled.length === runs.length
        ? "All four at the boil"
        : boiled.length && m - boiled[0]!.boilMin! < 4
          ? `${boiled[0]!.setup.name}: at the boil`
          : boiled.length && waiting
            ? levelled
              ? `${waiting.r.setup.name}: simmering at ${degC(waiting.c)} °C`
              : `${waiting.r.setup.name}: ${degC(waiting.c)} °C`
            : `${litres} litres heating from ${TAP_C} °C`;
    const note = boiled.length ? `${boiled.length} of ${runs.length} boiling` : "Heating";
    const sig = `${Math.floor(m)}|${doing}|${note}`;
    if (sig === lastHud) return;
    lastHud = sig;
    if (hudClock) hudClock.textContent = span(Math.floor(m));
    if (hudDoing) hudDoing.textContent = doing;
    if (hudNote) {
      hudNote.textContent = note;
      hudNote.dataset.kind = boiled.length ? "boil" : "heat";
    }
  }
  /** The labels over each way's pots: its temperature, or that it is boiling. The lowest on screen sits just over
   *  its pots; one that would cover another is lifted clear, on a longer stem. */
  const tagAt = new THREE.Vector3();
  const tagSize = new Map<HTMLElement, { w: number; h: number }>();
  const tagStem = new Map<HTMLElement, string>();
  const tagItems: { el: HTMLElement; x: number; y: number; w: number; h: number }[] = [];
  const tagPlaced: { l: number; r: number; t: number; b: number }[] = [];
  function showTags(m: number): void {
    const W = stageEl!.clientWidth;
    const H = stageEl!.clientHeight;
    tagItems.length = 0;
    for (const s of SETUPS) {
      const el = tagEls.get(s.key);
      if (!el) continue;
      const ps = pots.filter((p) => p.setup === s);
      tagAt.set(0, 0, 0);
      for (const p of ps) tagAt.add(p.top);
      tagAt.divideScalar(ps.length);
      tagAt.y += 0.03;
      tagAt.project(camera);
      const on = tagAt.z < 1 && Math.abs(tagAt.x) < 0.98 && Math.abs(tagAt.y) < 0.98;
      el.hidden = !on;
      if (!on) continue;
      const a = at(runOf(s.key), m);
      const text = a.boiled ? "Boiling" : `${degC(a.c)} °C`;
      if (el.dataset.last !== text || !tagSize.has(el)) {
        el.dataset.last = text;
        const v = el.querySelector("b");
        if (v) v.textContent = text;
        el.toggleAttribute("data-boiled", a.boiled);
        tagSize.set(el, { w: el.offsetWidth, h: el.offsetHeight });
      }
      const sz = tagSize.get(el)!;
      tagItems.push({ el, x: ((tagAt.x + 1) / 2) * W, y: ((1 - tagAt.y) / 2) * H, w: sz.w, h: sz.h });
    }
    tagItems.sort((a, b) => b.y - a.y);
    tagPlaced.length = 0;
    for (const it of tagItems) {
      const l = it.x - it.w / 2 - 3;
      const r = it.x + it.w / 2 + 3;
      let b = it.y - 6;
      for (let moved = true, n = 0; moved && n < 8; n++) {
        moved = false;
        for (const p of tagPlaced)
          if (l < p.r && r > p.l && b > p.t - 3 && b - it.h < p.b + 3) {
            b = p.t - 3;
            moved = true;
          }
      }
      tagPlaced.push({ l, r, t: b - it.h, b });
      // Kept inside the stage's sides; the stem still drops to the pots.
      const x = Math.min(W - it.w / 2 - 4, Math.max(it.w / 2 + 4, it.x));
      const stem = `${Math.round(it.y - b)}|${Math.round(it.x - x)}`;
      it.el.style.transform = `translate(${x.toFixed(1)}px, ${b.toFixed(1)}px)`;
      if (tagStem.get(it.el) !== stem) {
        tagStem.set(it.el, stem);
        it.el.style.setProperty("--stem", `${Math.round(it.y - b)}px`);
        it.el.style.setProperty("--sx", `${Math.round(it.x - x)}px`);
      }
    }
  }
  /** Each pot's water, flame and lid at minute `m` (`lit`: how far the burners are lit, 0..1). */
  function showPots(m: number, lit: number, t: number): void {
    ranges.forEach((r, k) => {
      for (let b = 0; b < 6; b++) {
        const v = pots.some((p) => p.range === k && p.burner === b) ? flame.v * lit : 0;
        if (Math.abs(v - r.flow[b]!) < 1e-3) continue;
        r.flow[b] = v;
        r.flames?.setFlow(b, v);
      }
    });
    for (const p of pots) {
      const a = at(runOf(p.setup.key), m);
      const area = Math.PI * (p.setup.pot.d / 2) ** 2;
      const body = p.ext.h - 0.036;
      const level = a.litres / 1000 / area / (body - 0.03 - (0.0085 + 0.0014));
      // The water stirs as it nears the boil: bubbles on the floor from about 60 °C, a simmer, then a rolling boil.
      const boil = a.boiled ? 1 : 0.8 * Math.min(1, Math.max(0, (a.c - 60) / 40)) ** 1.6;
      if (Math.abs(level - p.level) > 2e-3) {
        p.level = level;
        p.model.hooks.set?.("level", level);
      }
      if (Math.abs(boil - p.boil) > 0.01 || (boil === 1) !== (p.boil === 1)) {
        p.boil = boil;
        p.model.hooks.set?.("boil", boil);
      }
      // The steam off it, as the model has it: off the open water, or past the lid, short of the boil; at the boil,
      // all the heat the burner gives that the pot does not lose otherwise.
      const wet = a.litres / 1000 / area;
      const lost = lossesW(a.c, p.setup.pot.d, wet, p.setup.lid);
      const steamW = a.boiled ? Math.max(0, BURNER_W * GAS_INTO_POT * flame.v - lost.radiation - lost.air) : lost.steam;
      const steam = steamW < 0.02 * STEAM_W ? 0 : Math.min(3, steamW / STEAM_W);
      if (Math.abs(steam - p.steam) > 0.02 || (steam === 0) !== (p.steam === 0)) {
        p.steam = steam;
        p.model.hooks.set?.("steam", steam);
      }
      // A lid on a pot at the boil rattles as the steam lifts it.
      if (p.setup.lid) p.model.hooks.open?.(a.boiled ? 0.02 + 0.02 * Math.max(0, Math.sin(t * 19 + p.phase)) : 0);
      p.model.hooks.tick?.(t, 0, 0);
    }
  }

  apply();
  if (chartEl) {
    let lastW = chartEl.clientWidth;
    new ResizeObserver(() => {
      if (Math.abs(chartEl.clientWidth - lastW) < 4) return;
      lastW = chartEl.clientWidth;
      renderChart();
    }).observe(chartEl);
  }
  void cook.ready.then(() => void gate.refresh(scene));

  let inView = false;
  observeInView(stageEl, (v) => {
    inView = v;
  });
  let last = performance.now();
  const fitClip = clipFitter(camera, scene, 0.05);
  function animate(now: number): void {
    requestAnimationFrame(animate);
    const dt = Math.max(0, Math.min(0.05, (now - last) / 1000));
    last = now;
    if (!inView || document.hidden || !gate.ready()) return;
    const loopS = LIGHT_S + HEAT_S + END_S;
    loopT = (loopT + dt) % loopS;
    const t = now / 1000;
    const lighting = frozen === null && loopT < LIGHT_S;
    playMin = frozen !== null ? Math.min(endMin, frozen) : lighting ? 0 : Math.min(endMin, ((loopT - LIGHT_S) / HEAT_S) * endMin);
    const lit = frozen !== null ? 1 : smooth(loopT / 0.8);
    showPots(playMin, lit, t);
    for (const r of ranges) r.m.hooks.tick?.(t, dt, 1);
    cook.stand(t);
    showHud(playMin, lighting);
    if (playX && playLine) {
      const xx = playX(playMin).toFixed(1);
      playLine.setAttribute("x1", xx);
      playLine.setAttribute("x2", xx);
      playClip?.setAttribute("width", String(Math.max(chartL, +xx)));
    }
    controls.update();
    fitClip();
    renderer.render(scene, camera);
    showTags(playMin);
  }
  requestAnimationFrame(animate);
  void gate.refresh(scene);

  // A hook for checking the page by script.
  Object.assign(stageEl, {
    boil: {
      info: () => ({
        mountMs: Math.round(mountMs),
        litres,
        flame: flame.key,
        endMin: +endMin.toFixed(1),
        t: +playMin.toFixed(1),
        runs: runs.map((r) => ({ key: r.setup.key, boilMin: r.boilMin === null ? null : +r.boilMin.toFixed(1), steamL: +(r.steamL * r.setup.pots).toFixed(2), peak: +r.c[r.c.length - 1]!.toFixed(1), burnerMin: +r.burnerMin.toFixed(1) })),
        tags: SETUPS.map((s) => tagEls.get(s.key)?.textContent ?? null),
        status: statusEl?.textContent ?? "",
      }),
      hold(min: number | null) {
        frozen = min;
      },
      set(s: { litres?: number; flame?: "full" | "three" | "half" }) {
        if (s.litres !== undefined) litres = Math.max(8, Math.min(MAX_LITRES, Math.round(s.litres)));
        if (s.flame !== undefined) flame = FLAMES.find((f) => f.key === s.flame) ?? flame;
        apply();
      },
    },
  });
  const mountMs = performance.now() - mountT0;
}
