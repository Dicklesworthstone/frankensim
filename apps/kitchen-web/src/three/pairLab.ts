/**
 * Two kitchens, one order (chapter 07, third figure): the morning's kitchen twice, one above the other on the same
 * clock: the hub, whose van leaves at 10:30, and a second kitchen across town, whose van leaves at 10:00.
 *
 * Everything shown follows the best split for the settings on the page (bestPair, sim/twoKitchens.ts, run in a
 * worker): each kitchen acts out its own planned morning (MorningKitchen), its carriers rolling out of the door when
 * its van leaves; with the hub alone, the kitchen across town stands dark. The chart shows why the split falls where
 * it does: for each dish in each kitchen, the stretch in which a batch must come off the heat to make that kitchen's
 * van and still be at its best at lunch, and when each batch did.
 */
import * as THREE from "three";
import { DOOR, SERVE, SLOTS, clockOf, menuDishes, type Crew, type Morning, type MorningDish } from "../sim/cateringMorning";
import type { PairAsk, PairDone } from "../sim/pairRuns.worker";
import { FAR_DOOR, limitKey, planKey, travels, type PairPlan, type PairSettings } from "../sim/twoKitchens";
import { bindCanvasResize, clipFitter, createScrollSafeOrbit, createStudioRenderer, observeInView, shaderGate } from "./kit";
import { CARRIER, CARRIER_COLS, CARRIER_X0, END_HOLD_S, LINE_Z, MORNING_PIECES, MorningKitchen, OVEN_X0, SEC_PER_MIN, SLOT_COLOR, TABLE, TABLE_Z, fmtLate, shortName, type CookSpec, type KitchenMoment } from "./morningKitchen";
import type { Piece } from "./models/all";
import { studioEnvironment, studioLights } from "./models/stage";
import { backdrop, tilePlane } from "./weekLab";

/** Two kitchens' pieces, for main.ts to build ahead. */
export const PAIR_PIECES: readonly Piece[] = [...MORNING_PIECES, ...MORNING_PIECES];

/** The hub's three cooks, and up to three across town. */
const HUB_COOKS: CookSpec[] = [
  { cohort: "p50", seed: 41, apron: 0x26344d },
  { cohort: "p95", seed: 57, apron: 0x3b3d43 },
  { cohort: "p05", seed: 63, apron: 0x45563f },
];
const FAR_COOKS: CookSpec[] = [
  { cohort: "p95", seed: 83, apron: 0x5a2f35 },
  { cohort: "p50", seed: 91, apron: 0x2f4a4a },
  { cohort: "p05", seed: 97, apron: 0x3e3552 },
];
const HUB_CREW: Crew = 3;
/** The kitchen across town's lights when it is not drafted. */
const DARK = 0.3;

type Mode = "alone" | "both";
type Side = "hub" | "far";

/** "a", "a and b", "a, b and c". */
const listOf = (xs: string[]) => (xs.length < 2 ? (xs[0] ?? "") : `${xs.slice(0, -1).join(", ")} and ${xs[xs.length - 1]}`);
const num = (v: number) => Math.round(v).toLocaleString("en-US");
const holdText = (min: number) => (!Number.isFinite(min) ? "all morning, cold" : min >= 120 ? `${min / 60} h` : `${min} min`);

export function mountPairLab(root: HTMLElement): void {
  const mountT0 = performance.now();
  const stageEl = root.querySelector<HTMLElement>("[data-pair-stage]");
  const halfEl = { hub: root.querySelector<HTMLElement>('[data-pr-half="hub"]'), far: root.querySelector<HTMLElement>('[data-pr-half="far"]') };
  const canvasEl = { hub: root.querySelector<HTMLCanvasElement>('[data-pair-canvas="hub"]'), far: root.querySelector<HTMLCanvasElement>('[data-pair-canvas="far"]') };
  if (!stageEl || !halfEl.hub || !halfEl.far || !canvasEl.hub || !canvasEl.far) return;
  const rHub = createStudioRenderer(canvasEl.hub, false);
  if (!rHub) return;
  const rFar = createStudioRenderer(canvasEl.far, false);
  if (!rFar) return;
  rHub.toneMappingExposure = 1.0;
  rFar.toneMappingExposure = 1.0;
  const phone = window.matchMedia("(max-width: 720px)").matches;

  /* Two scenes, one kitchen in each, lit alike. */
  const env = studioEnvironment(rHub);
  const makeScene = (kitchen: MorningKitchen) => {
    const scene = new THREE.Scene();
    scene.background = backdrop();
    scene.fog = new THREE.Fog(0x151510, 8, 16);
    scene.environment = env;
    scene.environmentIntensity = 0.6;
    const lights = studioLights(scene, new THREE.Vector3(-0.4, 0.9, 0.1), 4.6, phone ? 1024 : 2048);
    lights.key.position.set(-1.2, 6.4, 4.6);
    lights.key.intensity = 1.9;
    lights.key.shadow.normalBias = 0.02;
    lights.rim.position.set(3.4, 3.2, -3.0);
    lights.rim.intensity = 0.75;
    lights.hemi.intensity = 0.5;
    // Wider than the morning's floor: a half's view is squarer, and its far edges would show.
    const floor = tilePlane(24, 16);
    floor.rotation.x = -Math.PI / 2;
    floor.position.set(0, 0, -1.5);
    scene.add(floor);
    kitchen.rollOut = true;
    scene.add(kitchen.group);
    const base = { key: lights.key.intensity, rim: lights.rim.intensity, hemi: lights.hemi.intensity, env: scene.environmentIntensity };
    /** Lights up (1) or down (DARK). */
    const lit = (f: number) => {
      lights.key.intensity = base.key * f;
      lights.rim.intensity = base.rim * f;
      lights.hemi.intensity = base.hemi * f;
      scene.environmentIntensity = base.env * f;
    };
    return { scene, lit };
  };
  const kitchen = { hub: new MorningKitchen(HUB_COOKS), far: new MorningKitchen(FAR_COOKS) };
  const view = { hub: makeScene(kitchen.hub), far: makeScene(kitchen.far) };
  const renderer = { hub: rHub, far: rFar };

  /* One camera for both: the kitchens stand alike in their scenes, so a drag on either turns both. */
  const camera = new THREE.PerspectiveCamera(32, 16 / 6, 0.05, 60);
  const TARGET = new THREE.Vector3();
  const VIEW = new THREE.Vector3();
  const keep: THREE.Vector3[] = [];
  const aim = () => {
    keep.length = 0;
    TARGET.set(-0.1, 0.75, 0.0);
    VIEW.set(0.1, 1.15, 1).normalize();
    for (const x of [OVEN_X0[0], CARRIER_X0 + CARRIER_COLS * (CARRIER.w + 0.03)]) keep.push(new THREE.Vector3(x, 1.85, LINE_Z - 0.4), new THREE.Vector3(x, 0.2, TABLE_Z + TABLE.d));
    // The top of the carriers' stack, nearer the camera than the line's top: on a narrow half it is the widest point.
    keep.push(new THREE.Vector3(CARRIER_X0 + CARRIER_COLS * (CARRIER.w + 0.03), 3 * CARRIER.h, TABLE_Z));
  };
  const controls = (["hub", "far"] as const).map((k) => createScrollSafeOrbit(camera, canvasEl[k]!, [TARGET.x, TARGET.y, TARGET.z], { minDistance: 1.2, maxDistance: 14, maxPolarAngle: Math.PI * 0.47 }));
  controls[1]!.target = controls[0]!.target;
  let userMoved = false;
  for (const c of controls)
    c.addEventListener("start", () => {
      userMoved = true;
    });
  const ndc = new THREE.Vector3();
  const hudEls = Array.from(stageEl.querySelectorAll<HTMLElement>(".pr-hud"));
  const frame = () => {
    if (userMoved) return;
    aim();
    const half = halfEl.hub!;
    const w = Math.max(1, half.clientWidth);
    const h = Math.max(1, half.clientHeight);
    camera.aspect = Math.max(0.2, w / h);
    camera.clearViewOffset();
    camera.updateProjectionMatrix();
    // The kitchen fits under the HUD (the taller of the two halves'): fitted to the band's height and the half's
    // width, then shifted into the band and across to the middle.
    const hudBottom = Math.max(0, ...hudEls.map((e) => e.offsetTop + e.offsetHeight));
    const top = Math.min(0.86, 1 - (2 * (hudBottom + 6)) / h);
    const bot = -0.94;
    let dist = 2.5;
    let lo = 0;
    let hi = 0;
    let left = 0;
    let right = 0;
    for (; dist < 20; dist += 0.06) {
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
      if (right - left < 1.9 && hi - lo < top - bot) break;
    }
    const shiftY = (top + bot) / 2 - (hi + lo) / 2;
    const shiftX = (left + right) / 2;
    if (Math.abs(shiftY) > 0.005 || Math.abs(shiftX) > 0.005) {
      camera.setViewOffset(w, h, (shiftX * w) / 2, (shiftY * h) / 2, w, h);
      camera.updateProjectionMatrix();
    }
    for (const v of [view.hub, view.far]) {
      const fog = v.scene.fog as THREE.Fog;
      fog.near = dist * 1.1;
      fog.far = dist * 2.6;
    }
    controls[0]!.target.copy(TARGET);
    for (const c of controls) c.update();
  };
  const gate = { hub: shaderGate(rHub, camera), far: shaderGate(rFar, camera) };
  bindCanvasResize(halfEl.hub, rHub, camera, gate.hub);
  bindCanvasResize(halfEl.far, rFar, camera, gate.far);
  frame();
  // The HUDs too: they wrap to a second line once their text is in, and the view must clear them.
  const ro = new ResizeObserver(() => frame());
  for (const el of [halfEl.hub, ...hudEls]) ro.observe(el);

  /* ---------------------------------------------------------------------------------------------------------- */
  /* The page.                                                                                                   */
  /* ---------------------------------------------------------------------------------------------------------- */

  const q = <T extends Element = HTMLElement>(s: string) => root.querySelector<T>(s);
  const modeBtns = Array.from(root.querySelectorAll<HTMLButtonElement>("[data-pr-mode]"));
  const crewBtns = Array.from(root.querySelectorAll<HTMLButtonElement>("[data-pr-crew]"));
  const menuEl = q("[data-pr-menu]");
  const mealsIn = q<HTMLInputElement>("[data-pr-meals]");
  const mealsVal = q("[data-pr-meals-val]");
  const statusEl = q("[data-pr-status]");
  const chartEl = q("[data-pair-chart]");
  const hudClock = q("[data-pr-hud-clock]");
  const hudNote = { hub: q('[data-pr-hud-note="hub"]'), far: q('[data-pr-hud-note="far"]') };
  const tile = (k: string) => ({ v: q(`[data-pr-${k}]`), d: q(`[data-pr-${k}-d]`), s: q(`[data-pr-${k}-s]`) });
  const tiles = { lunch: tile("lunch"), past: tile("past"), far: tile("far"), limit: tile("limit") };
  const setTile = (t: (typeof tiles)[keyof typeof tiles], v: string, s: string, d: string | null = null, good = true) => {
    if (t.v) t.v.textContent = v;
    if (t.s) t.s.textContent = s;
    if (t.d) {
      t.d.hidden = d === null;
      if (d !== null) {
        t.d.textContent = d;
        t.d.dataset.good = String(good);
      }
    }
  };

  const settings: PairSettings = { meals: Number(mealsIn?.value ?? 1400), menu: SLOTS.length === 3 ? 7 : 0, hubCooks: HUB_CREW, farCooks: 2 };
  let mode: Mode = "both";
  const plans = new Map<string, Extract<PairDone, { kind: "plan" }>>();
  const limits = new Map<string, Extract<PairDone, { kind: "limit" }>>();
  /** The plan on screen and the settings it is for (it stays while the next is planned). */
  let shown: PairPlan | null = null;
  let shownKey = "";

  /* The menu card: each place's dish, how long it holds and where it is cooked; a press swaps it. */
  const slotBtns: HTMLButtonElement[] = [];
  let sidesLi: HTMLLIElement | null = null;
  if (menuEl) {
    menuEl.innerHTML = "";
    SLOTS.forEach((_, k) => {
      const li = document.createElement("li");
      const b = document.createElement("button");
      b.type = "button";
      b.className = "wk-dish mo-dish";
      b.style.setProperty("--dish", SLOT_COLOR[k]!);
      b.addEventListener("click", () => {
        settings.menu ^= 1 << k;
        changed(true);
      });
      li.appendChild(b);
      menuEl.appendChild(li);
      slotBtns.push(b);
    });
    sidesLi = document.createElement("li");
    sidesLi.className = "wk-dish mo-dish pr-sides";
    sidesLi.style.setProperty("--dish", SLOT_COLOR[SLOTS.length]!);
    menuEl.appendChild(sidesLi);
  }
  const holdBars = (min: number) => {
    const n = min >= 210 ? 4 : min >= 140 ? 3 : min >= 85 ? 2 : 1;
    return `<span class="wk-keeps" aria-hidden="true">${[1, 2, 3, 4].map((k) => `<i${k <= n ? ' class="on"' : ""}></i>`).join("")}<small>${Number.isFinite(min) ? (min >= 120 ? `${min / 60} h` : `${min} min`) : "cold"}</small></span>`;
  };
  const ready = () => shown !== null && shownKey === planKey(settings);
  /** Where the plan on screen cooks place k. */
  const whereOf = (k: number, d: MorningDish) => {
    if (!ready()) return "…";
    const part = shown!.parts[k] ?? 0;
    if (part >= 1) return "Across town";
    if (part > 0) return "Half and half";
    return mode === "both" && !travels(d.holdMin) ? "Hub only" : "At the hub";
  };
  function renderMenu(): void {
    const dishes = menuDishes(settings.menu);
    slotBtns.forEach((b, k) => {
      const d = dishes[k]!;
      const swapped = (settings.menu & (1 << k)) !== 0;
      const other = swapped ? SLOTS[k]!.written : SLOTS[k]!.swap;
      const where = whereOf(k, d);
      b.innerHTML = `<b>${d.name}</b>${holdBars(d.holdMin)}<em data-where="${where === "Across town" || where === "Half and half" ? "far" : where === "Hub only" ? "stuck" : "hub"}">${where}</em>`;
      b.setAttribute("aria-pressed", String(swapped));
      b.title = `Swap for ${other.name.toLowerCase()}`;
    });
    if (sidesLi) {
      const d = dishes[SLOTS.length]!;
      const where = whereOf(SLOTS.length, d);
      sidesLi.innerHTML = `<b>${d.name}</b>${holdBars(d.holdMin)}<em data-where="${where === "Across town" || where === "Half and half" ? "far" : "hub"}">${where}</em>`;
    }
    for (const b of modeBtns) b.setAttribute("aria-pressed", String(b.dataset.prMode === mode));
    for (const b of crewBtns) b.setAttribute("aria-pressed", String(Number(b.dataset.prCrew) === settings.farCooks));
    if (mealsVal) mealsVal.textContent = `${num(settings.meals)} meals`;
  }

  /* The chart: for each kitchen and each dish it cooks, the stretch a batch must come off the heat in (at its best at
   * lunch, and in time for that kitchen's van), and when each batch did; the vans, the drive and lunch. */
  let X0 = 150;
  let X1 = 600;
  let T1 = 360;
  const xOf = (m: number) => X0 + (Math.max(0, Math.min(m, T1)) / T1) * (X1 - X0);
  let cursor: SVGLineElement | null = null;
  const ROW_H = 20;
  const HEAD_H = 20;
  function renderChart(): void {
    if (!chartEl) return;
    const p = ready() ? shown : null;
    const W = Math.max(300, Math.round(chartEl.clientWidth - 2 * 12));
    const narrow = W < 480;
    X0 = narrow ? 100 : 168;
    X1 = W - 8;
    const lunch = SERVE + (p?.lunchLate ?? 0);
    T1 = Math.max(360, Math.ceil((lunch + 10) / 60) * 60);
    const dishes = menuDishes(settings.menu);
    const groups: { side: Side; door: number; m: Morning; portions: number[] }[] = [];
    if (p) {
      groups.push({ side: "hub", door: DOOR, m: p.hub, portions: p.here });
      if (p.far) groups.push({ side: "far", door: FAR_DOOR, m: p.far, portions: p.there });
    }
    let y = 8;
    let body = "";
    const r = (a: number, b: number, yy: number, h: number, attrs: string) => (b - a > 0.05 ? `<rect x="${xOf(a).toFixed(1)}" y="${yy.toFixed(1)}" width="${Math.max(0.8, xOf(b) - xOf(a)).toFixed(1)}" height="${h}" ${attrs}/>` : "");
    for (const g of groups) {
      const late = g.m.late > 0.5;
      const name = g.side === "hub" ? "The hub" : "Across town";
      // The group's head: its name, the van, and its road to lunch.
      body += `<text class="gh" x="${X0 - 8}" y="${y + 13}">${name}</text>`;
      body += `<line class="road" x1="${xOf(g.m.depart)}" x2="${xOf(lunch)}" y1="${y + 9}" y2="${y + 9}"/>`;
      const roadW = xOf(lunch) - xOf(g.m.depart);
      if (roadW > 34) body += `<text class="roadl" x="${((xOf(g.m.depart) + xOf(lunch)) / 2).toFixed(1)}" y="${y + 6}">${roadW > 64 ? "on the road" : "road"}</text>`;
      const top = y + HEAD_H - 4;
      let rowsH = 0;
      dishes.forEach((d, k) => {
        if (!(g.portions[k]! > 0)) return;
        const ry = y + HEAD_H + rowsH;
        rowsH += ROW_H;
        const label = narrow ? (k === SLOTS.length ? "Sides" : shortName(d.name)) : `${shortName(d.name)} · ${num(g.portions[k]!)}`;
        body += `<text class="ya" x="${X0 - 8}" y="${ry + ROW_H / 2 + 3.5}">${label}</text>`;
        body += `<rect class="lane" x="${X0}" y="${ry + 2}" width="${X1 - X0}" height="${ROW_H - 4}"/>`;
        // The window: off the heat no sooner than lunch less its holding time, and in time for the van.
        const from = Number.isFinite(d.holdMin) ? SERVE - d.holdMin : 0;
        if (from < g.door) body += r(Math.max(0, from), g.door, ry + 3, ROW_H - 6, `class="win" fill="${SLOT_COLOR[k]}"`);
        // Each batch: when it came off the heat (the sides, cold, when they were boxed).
        for (const b of g.m.batches) {
          if (b.slot !== k) continue;
          const at = b.cooked >= 0 ? b.cooked : b.packed;
          if (at < 0) continue;
          const bad = (Number.isFinite(d.holdMin) && b.cooked >= 0 && lunch - b.cooked > d.holdMin + 1e-6) || b.packed > g.m.depart + 1e-6 || at > g.door + 0.5;
          const x = xOf(at).toFixed(1);
          body += `<line class="tick${bad ? " bad" : ""}" x1="${x}" x2="${x}" y1="${ry + 4}" y2="${ry + ROW_H - 4}"${bad ? "" : ` stroke="${SLOT_COLOR[k]}"`}/>`;
        }
      });
      // The van: due, and when it left if later.
      const bottom = y + HEAD_H + rowsH;
      body += `<line class="van" x1="${xOf(g.door)}" x2="${xOf(g.door)}" y1="${top}" y2="${bottom}"/>`;
      if (late) body += `<rect class="late" x="${xOf(g.door)}" y="${top}" width="${(xOf(g.m.depart) - xOf(g.door)).toFixed(1)}" height="${bottom - top}"/><line class="van late" x1="${xOf(g.m.depart)}" x2="${xOf(g.m.depart)}" y1="${top}" y2="${bottom}"/>`;
      const vx = xOf(late ? g.m.depart : g.door);
      body += `<text class="vanl r${late ? " late" : ""}" x="${(vx - 4).toFixed(1)}" y="${y + 13}">van ${clockOf(late ? g.m.depart : g.door)}</text>`;
      y = bottom + 10;
    }
    const axisY = y + 2;
    const H = axisY + 34;
    let svg = `<svg viewBox="0 0 ${W} ${H}" width="${W}" height="${H}" class="rs-svg mo-svg pr-svg" aria-hidden="true">`;
    // An hour's gridline; a label every hour, or every other on a narrow chart.
    const every = (X1 - X0) / (T1 / 60) < 40 ? 120 : 60;
    for (let t = 0; t <= T1; t += 60) svg += `<line class="g v" x1="${xOf(t)}" x2="${xOf(t)}" y1="4" y2="${axisY}"/>${xOf(t) > X1 - 14 || t % every ? "" : `<text class="xa" x="${xOf(t)}" y="${axisY + 14}">${clockOf(t)}</text>`}`;
    svg += body;
    // Lunch, at the client.
    svg += `<line class="lunch${p && p.lunchLate > 0.5 ? " late" : ""}" x1="${xOf(lunch)}" x2="${xOf(lunch)}" y1="4" y2="${axisY}"/><text class="lunchl${p && p.lunchLate > 0.5 ? " late" : ""}${xOf(lunch) > X1 - 70 ? " r" : ""}" x="${(xOf(lunch) + (xOf(lunch) > X1 - 70 ? -4 : 4)).toFixed(1)}" y="${axisY + 28}">lunch ${clockOf(lunch)}</text>`;
    if (!p) svg += `<text class="wait" x="${((X0 + X1) / 2).toFixed(1)}" y="${axisY / 2 + 4}">Planning…</text>`;
    svg += `<line class="cur" x1="${X0}" x2="${X0}" y1="4" y2="${axisY + 2}"/></svg>`;
    chartEl.innerHTML = svg;
    cursor = chartEl.querySelector<SVGLineElement>(".cur");
  }

  /* Tiles and status: the plan on screen, against the other mode's. */
  function renderTiles(): void {
    if (!ready()) {
      for (const t of Object.values(tiles)) setTile(t, "…", "");
    } else {
      const r = plans.get(shownKey)!;
      const p = shown!;
      const other = mode === "both" ? r.alone : null;
      // Against the hub alone, when both kitchens are on screen and the hub alone does differently.
      const lateText = (x: PairPlan) => (x.lunchLate > 0.5 ? fmtLate(x.lunchLate) : "On time");
      const vans = p.far ? `vans at ${clockOf(p.far.depart)} and ${clockOf(p.hub.depart)}` : `the hub's van at ${clockOf(p.hub.depart)}`;
      setTile(tiles.lunch, lateText(p), other && Math.abs(other.lunchLate - p.lunchLate) > 0.5 ? `the hub alone: ${lateText(other).toLowerCase()}` : vans);
      setTile(tiles.past, p.pastBest ? num(p.pastBest) : "None", other && other.pastBest !== p.pastBest ? `the hub alone: ${other.pastBest ? num(other.pastBest) : "none"}` : p.pastBest ? "portions, by lunch" : "every portion at its best");
      const there = p.there.reduce((a, b) => a + b, 0);
      if (mode === "alone") setTile(tiles.far, "Nothing", "the kitchen stands dark");
      else if (!p.far) setTile(tiles.far, "Nothing", menuDishes(settings.menu).some((d, k) => k < SLOTS.length && travels(d.holdMin)) ? "no split does better here" : "only the sides keep long enough");
      else setTile(tiles.far, `${Math.round((100 * p.farHands) / Math.max(1, p.farHands + p.hubHands))}%`, `of the cooks' work: ${num(there)} portions`);
    }
    const l = limits.get(limitKey(settings.menu, HUB_CREW, settings.farCooks));
    if (!l) setTile(tiles.limit, "…", "");
    else if (mode === "both") setTile(tiles.limit, num(l.both), `the hub alone: ${num(l.alone)}`, l.both > l.alone ? `+${num(l.both - l.alone)}` : null, true);
    else setTile(tiles.limit, num(l.alone), `with the kitchen across town: ${num(l.both)}`);
  }

  function renderStatus(): void {
    if (!statusEl) return;
    if (!ready()) {
      statusEl.textContent = mode === "both" ? "Planning both kitchens…" : "Planning the hub's morning…";
      statusEl.dataset.state = "warn";
      return;
    }
    const r = plans.get(shownKey)!;
    const p = shown!;
    const dishes = menuDishes(settings.menu);
    const outcome = (x: PairPlan) => {
      const bits: string[] = [];
      if (x.lunchLate > 0.5) bits.push(`lunch ${fmtLate(x.lunchLate).toLowerCase()}`);
      if (x.pastBest) bits.push(`${num(x.pastBest)} portions past their best`);
      return bits.length ? listOf(bits) : "lunch on time with every dish at its best";
    };
    /** What a plan has cooked across town: "the braise and half the salad". */
    const movedOf = (x: PairPlan) => dishes.flatMap((d, k) => (x.there[k]! > 0 ? [`${x.parts[k]! >= 1 ? "the" : "half the"} ${d.name.toLowerCase()}`] : []));
    const moved = movedOf(p);
    // The dishes on the menu that cannot travel: they keep too short a time for the van across town.
    const stuck = dishes.filter((d, k) => k < SLOTS.length && !travels(d.holdMin));
    const why =
      stuck.length > 0
        ? ` The ${listOf(stuck.map((d) => d.name.toLowerCase()))} ${stuck.length > 1 ? "keep" : "keeps"} ${listOf(stuck.map((d) => `${d.holdMin} min`))}, so ${stuck.length > 1 ? "they have" : "it has"} to come off the heat after ${listOf(stuck.map((d) => clockOf(SERVE - d.holdMin)))}: too late for a van at ${clockOf(FAR_DOOR)}.`
        : "";
    let text: string;
    if (mode === "alone") {
      text = `The hub alone: ${outcome(p)}, its van away at ${clockOf(p.hub.depart)}.`;
      if (!p.made && r.both.made) text += ` With the kitchen across town making ${listOf(movedOf(r.both))}, lunch is on time with every dish at its best.`;
    } else if (p.made) {
      text = moved.length
        ? `Across town they make ${listOf(moved)}, and their van leaves at ${clockOf(p.far!.depart)}; the hub's leaves at ${clockOf(p.hub.depart)}. Lunch is on time with every dish at its best.${r.alone.made ? "" : ` The hub alone: ${outcome(r.alone)}.`}`
        : `The hub manages this order on its own, lunch on time with every dish at its best, so nothing goes across town.`;
    } else if (moved.length) {
      text = `Even with ${listOf(moved)} made across town: ${outcome(p)}.${why}${stuck.length ? ` Swap ${stuck.length > 1 ? "them" : "it"} for a dish that keeps, or lower the order.` : settings.farCooks < 3 ? " Add a cook across town, or lower the order." : " Lower the order."}`;
    } else {
      text = `A second kitchen can't help this order: ${outcome(p)}.${why}${stuck.length ? ` Swap ${stuck.length > 1 ? "them" : "it"} for a dish that keeps.` : ""}`;
    }
    statusEl.textContent = text;
    statusEl.dataset.state = p.made ? "ok" : p.lunchLate <= 30 && p.pastBest === 0 ? "warn" : "danger";
  }

  /* ---------------------------------------------------------------------------------------------------------- */
  /* The plans, from the worker.                                                                                 */
  /* ---------------------------------------------------------------------------------------------------------- */

  let worker: Worker | null = null;
  function ask(): void {
    if (!worker) {
      try {
        worker = new Worker(new URL("../sim/pairRuns.worker.ts", import.meta.url), { type: "module" });
      } catch {
        if (statusEl) statusEl.textContent = "This browser can't plan the kitchens off the page.";
        return;
      }
      worker.onmessage = (e: MessageEvent<PairDone>) => {
        const d = e.data;
        if (d.kind === "plan") {
          // Where each dish went, for the hub-alone status line.
          plans.set(d.key, d);
          if (d.key === planKey(settings)) apply();
        } else {
          limits.set(d.key, d);
          renderTiles();
        }
      };
    }
    const key = planKey(settings);
    if (!plans.has(key)) worker.postMessage({ kind: "plan", s: { ...settings } } satisfies PairAsk);
    if (!limits.has(limitKey(settings.menu, HUB_CREW, settings.farCooks))) worker.postMessage({ kind: "limit", menu: settings.menu, hubCooks: HUB_CREW, farCooks: settings.farCooks } satisfies PairAsk);
  }

  let endMin = DOOR + 6;
  let loopS = endMin * SEC_PER_MIN + END_HOLD_S;
  let tMin = 0;
  let loopT = 0;
  /** Put the plan for the settings on screen, if it has come, and everything that shows it. */
  function apply(): void {
    const key = planKey(settings);
    const r = plans.get(key);
    if (r) {
      const p = mode === "both" ? r.both : r.alone;
      const frac = endMin > 0 ? tMin / endMin : 0;
      shown = p;
      shownKey = key;
      let last = 0;
      for (const m of p.far ? [p.hub, p.far] : [p.hub]) for (const j of m.jobs) last = Math.max(last, j.t1);
      endMin = Math.max(DOOR, p.hub.depart, p.far?.depart ?? 0, last) + 6;
      let loop = endMin * SEC_PER_MIN + END_HOLD_S;
      loop = kitchen.hub.plan(p.hub, HUB_CREW, endMin, loop);
      const farLoop = kitchen.far.plan(p.far, settings.farCooks, endMin, loop);
      if (farLoop > loop + 1e-6) loop = kitchen.hub.plan(p.hub, HUB_CREW, endMin, farLoop);
      loopS = Math.max(loop, farLoop);
      view.far.lit(p.far ? 1 : DARK);
      tMin = frac * endMin;
      loopT = tMin * SEC_PER_MIN;
    }
    renderMenu();
    renderChart();
    renderTiles();
    renderStatus();
    lastHud = "";
  }
  let pending = 0;
  /** The settings changed: ask for their plan (a slider's steps only once it rests a moment). */
  function changed(now: boolean): void {
    window.clearTimeout(pending);
    const go = () => {
      ask();
      apply();
    };
    if (now) go();
    else {
      renderMenu();
      pending = window.setTimeout(go, 150);
    }
  }
  for (const b of modeBtns)
    b.addEventListener("click", () => {
      mode = b.dataset.prMode === "alone" ? "alone" : "both";
      changed(true);
    });
  for (const b of crewBtns)
    b.addEventListener("click", () => {
      settings.farCooks = Number(b.dataset.prCrew) as Crew;
      changed(true);
    });
  mealsIn?.addEventListener("input", () => {
    settings.meals = Number(mealsIn.value);
    changed(false);
  });

  /* The readouts on the stage, each minute. */
  let lastHud = "";
  function showHud(t: number, a: KitchenMoment, b: KitchenMoment): void {
    const p = shown;
    const note = (side: Side, m: Morning | null, door: number, k: KitchenMoment): [string, string] => {
      if (!m) return ["Not drafted", "idle"];
      if (t >= m.depart) return [`Van gone at ${clockOf(m.depart)}`, m.late > 0.5 ? "late" : "held"];
      if (t >= door) return ["The van is waiting", "late"];
      return [`${k.packed} of ${k.hot} hot pans packed`, side === "hub" || k.packed > 0 ? "held" : "idle"];
    };
    const [hn, hk] = note("hub", p?.hub ?? null, DOOR, a);
    const [fn, fk] = note("far", p?.far ?? null, FAR_DOOR, b);
    const sig = `${Math.floor(t)}|${hn}|${fn}`;
    if (sig === lastHud) return;
    lastHud = sig;
    if (hudClock) hudClock.textContent = clockOf(t);
    if (hudNote.hub) {
      hudNote.hub.textContent = hn;
      hudNote.hub.dataset.kind = hk;
    }
    if (hudNote.far) {
      hudNote.far.textContent = fn;
      hudNote.far.dataset.kind = fk;
    }
  }

  renderMenu();
  apply();
  ask();
  if (chartEl) {
    let lastW = chartEl.clientWidth;
    new ResizeObserver(() => {
      if (Math.abs(chartEl.clientWidth - lastW) < 4) return;
      lastW = chartEl.clientWidth;
      renderChart();
    }).observe(chartEl);
  }
  for (const side of ["hub", "far"] as const) for (const k of kitchen[side].cooks) void k.st.ready.then(() => void gate[side].refresh(view[side].scene));

  let inView = false;
  observeInView(stageEl, (v) => {
    inView = v;
  });
  let frozen: number | null = null;
  let last = performance.now();
  const fitClip = clipFitter(camera, view.hub.scene, 0.05);
  function animate(now: number): void {
    requestAnimationFrame(animate);
    const dt = Math.max(0, Math.min(0.05, (now - last) / 1000));
    last = now;
    if (!inView || document.hidden || !gate.hub.ready() || !gate.far.ready()) return;
    if (frozen === null) loopT = (loopT + dt) % loopS;
    else loopT = frozen * SEC_PER_MIN;
    tMin = Math.min(endMin, loopT / SEC_PER_MIN);
    kitchen.hub.pose(loopT);
    kitchen.far.pose(loopT);
    const a = kitchen.hub.show(tMin, dt);
    const b = kitchen.far.show(tMin, dt);
    if (a.minute || b.minute || !lastHud) showHud(tMin, a, b);
    if (cursor) {
      const x = xOf(tMin).toFixed(1);
      cursor.setAttribute("x1", x);
      cursor.setAttribute("x2", x);
    }
    for (const c of controls) c.update();
    fitClip();
    renderer.hub.render(view.hub.scene, camera);
    renderer.far.render(view.far.scene, camera);
  }
  requestAnimationFrame(animate);
  void gate.hub.refresh(view.hub.scene);
  void gate.far.refresh(view.far.scene);

  // A hook for checking the page by script.
  const brief = (m: Morning | null) => (m ? { depart: clockOf(m.depart), late: Math.round(m.late), pastBest: m.pastBest, jobs: m.jobs.length } : null);
  Object.assign(stageEl, {
    pair: {
      info: () => ({
        mountMs: Math.round(mountMs),
        key: planKey(settings),
        mode,
        ready: ready(),
        made: shown?.made ?? null,
        parts: shown?.parts ?? null,
        lunchLate: shown ? Math.round(shown.lunchLate) : null,
        pastBest: shown?.pastBest ?? null,
        hub: brief(shown?.hub ?? null),
        far: brief(shown?.far ?? null),
        alone: plans.get(planKey(settings)) ? { made: plans.get(planKey(settings))!.alone.made, lunchLate: Math.round(plans.get(planKey(settings))!.alone.lunchLate) } : null,
        limit: limits.get(limitKey(settings.menu, HUB_CREW, settings.farCooks)) ?? null,
        loopS: +loopS.toFixed(2),
        t: +tMin.toFixed(1),
        cooks: { hub: kitchen.hub.cooks.filter((k) => k.st.group.visible).length, far: kitchen.far.cooks.filter((k) => k.st.group.visible).length },
        holdText: holdText(menuDishes(settings.menu)[0]!.holdMin),
      }),
      hold(min: number | null) {
        frozen = min;
      },
      view(eye: [number, number, number], at: [number, number, number]) {
        userMoved = true;
        camera.position.set(...eye);
        controls[0]!.target.set(...at);
        camera.clearViewOffset();
        camera.updateProjectionMatrix();
        for (const c of controls) c.update();
      },
      set(s: Partial<PairSettings & { mode: Mode }>) {
        if (s.mode) mode = s.mode;
        if (s.meals !== undefined) settings.meals = s.meals;
        if (s.menu !== undefined) settings.menu = s.menu;
        if (s.farCooks !== undefined) settings.farCooks = s.farCooks;
        if (mealsIn) mealsIn.value = String(settings.meals);
        changed(true);
      },
    },
  });
  const mountMs = performance.now() - mountT0;
}
