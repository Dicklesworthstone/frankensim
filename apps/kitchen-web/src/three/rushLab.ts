/**
 * The lunch rush (chapter 06, third figure): a fast-casual hot line seen from the guests' side of the pass.
 *
 * Everything shown follows `runRush` (sim/serviceRush.ts) for the settings on the page, on one day: the median day
 * by orders of those simulated, the same day for every setting, so a change of plan replays the same guests.
 *   - the tickets on the rail over the pass are the orders waiting at the clock's time: paper white when new, amber
 *     once a guest has waited 8 minutes, red past 15;
 *   - the holding cabinet's rails carry the held portions, a 2.5-inch half pan per 12 (8 for burgers);
 *   - the fryer's baskets drop and the griddle and the combi run by the batches each has cooking;
 *   - the bags on the pass are the orders handed out in the last two minutes;
 *   - the three cooks work their stations and plate at the pass for the minutes the model has them doing each in
 *     every ten-minute stretch, and stand when it has them idle: a kinematic illustration (people.ts Routines).
 */
import * as THREE from "three";
import { KEEPERS, RUSH_MENU, SERVICE, modesFor, runRush, rushMenu, type DayRun, type Mode, type RushDay, type RushRun, type RushSettings } from "../sim/serviceRush";
import { bindCanvasResize, clipFitter, createScrollSafeOrbit, createStudioRenderer, observeInView, shaderGate } from "./kit";
import { buildFixture, type Model } from "./models/all";
import { palette } from "./models/materials";
import { Routine, makeStaff, type RoutineStep, type Staff, type WorkKind } from "./models/people";
import { studioEnvironment, studioLights } from "./models/stage";
import { Instances, backdrop, place, tilePlane } from "./weekLab";

/* ------------------------------------------------------------------------------------------------------------ */
/* The line.                                                                                                     */
/* ------------------------------------------------------------------------------------------------------------ */

/** Animation seconds per minute of the service, minutes per stretch the cooks' work is laid out in. */
const SEC_PER_MIN = 0.22;
const BIN = 10;
/** The pause on the finished day before the clock starts again. */
const END_HOLD_S = 2.5;
/** Simulated days per setting. */
const DAYS_N = 120;

/** The hot line's fronts face +z at LINE_Z; the cooks stand at it, backs to the expo pass behind them. */
const LINE_Z = -0.55;
const STATION_Z = LINE_Z + 0.42;
const PASS_Z = STATION_Z + 0.62;
const CAB = { id: "NSF-HOLDING-CABINET-EMPTY", w: 0.66, d: 0.813, h: 1.778, x0: -2.16 };
const FRYER = { id: "NSF-FRYER-TWIN-50LB", w: 0.394, d: 0.762, h: 1.143, x0: -1.42 };
const GRIDDLE = { id: "NSF-GRIDDLE-36", w: 0.914, d: 0.813, h: 0.381, x0: -0.96 };
const STAND_H = 0.6;
/** The combi: a 6-pan on a 48-inch worktable, the table's free end its cook's landing. */
const TABLE = { id: "NSF-PREPTABLE-48", w: 1.219, d: 0.762, h: 0.914, x0: 0.02 };
const COMBI = { id: "NSF-COMBI-6PAN", w: 0.847, d: 0.776, h: 0.782, x0: 0.04 };
const PASS = { id: "NSF-PREPTABLE-96", w: 2.438, d: 0.762, h: 0.914, x0: -1.62 };
const LAMP = { id: "NSF-HEATLAMP-STRIP-48", w: 1.219, d: 0.0636, h: 0.152 };
const PAN = { id: "NSF-HOTEL-PAN-HALF-4IN", w: 0.325, d: 0.265, h: 0.065 };
/** The overshelf above the pass: its underside, its depth, centred over the table. */
const SHELF_Y = 1.5;
const SHELF_D = 0.3;
const SHELF_Z = PASS_Z + PASS.d / 2;

/** The cooks, one per station in STATION_KEYS order: what their hands do at it and where, and where they plate. */
const COOKS = [
  { cohort: "p50" as const, x: FRYER.x0 + FRYER.w / 2, seed: 71, apron: 0x26344d, work: "stir" as WorkKind, at: [0, 0.95, LINE_Z - 0.16], surf: 0.9, land: [0.12, 0.905, LINE_Z - 0.04] },
  { cohort: "p95" as const, x: GRIDDLE.x0 + GRIDDLE.w / 2, seed: 83, apron: 0x3b3d43, work: "saute" as WorkKind, at: [0.06, STAND_H + 0.31, LINE_Z - 0.24], surf: STAND_H + 0.3, land: [-0.26, STAND_H + 0.29, LINE_Z - 0.05] },
  { cohort: "p05" as const, x: TABLE.x0 + 0.86, seed: 97, apron: 0x45563f, work: "stir" as WorkKind, at: [-0.3, TABLE.h + 0.3, LINE_Z - 0.06], surf: TABLE.h + 0.12, land: [0.2, TABLE.h + 0.004, LINE_Z - 0.3] },
];

const TICKET = { w: 0.07, h: 0.1, step: 0.078, perRow: 26 };
const AGE_AMBER = 8;
const AGE_RED = 15;
const COLOR = { paper: new THREE.Color(0xf1ece2), amber: new THREE.Color(0xf2b134), red: new THREE.Color(0xf0525f) };

const clock = (m: number) => {
  const t = Math.max(0, Math.round(m));
  const h = 10 + Math.floor(t / 60);
  return `${h > 12 ? h - 12 : h}:${String(t % 60).padStart(2, "0")}`;
};
const fmtMin = (m: number) => (m < 0.95 ? "<1 min" : `${m >= 9.95 ? m.toFixed(0) : m.toFixed(1)} min`);

/* ------------------------------------------------------------------------------------------------------------ */
/* Small builds: the griddle's stand, the overshelf, the ticket rail, a takeout bag.                             */
/* ------------------------------------------------------------------------------------------------------------ */

/** An open stainless equipment stand: top, four legs on bullet feet, an undershelf. Front-left corner at the origin. */
function equipmentStand(w: number, d: number, h: number): THREE.Group {
  const p = palette();
  const g = new THREE.Group();
  const add = (geo: THREE.BufferGeometry, mat: THREE.Material, x: number, y: number, z: number) => {
    const m = new THREE.Mesh(geo, mat);
    m.position.set(x, y, z);
    m.castShadow = m.receiveShadow = true;
    g.add(m);
  };
  add(new THREE.BoxGeometry(w, 0.032, d), p.ss, w / 2, h - 0.016, -d / 2);
  const leg = new THREE.BoxGeometry(0.038, h - 0.05, 0.038);
  const foot = new THREE.CylinderGeometry(0.016, 0.022, 0.03, 16);
  for (const x of [0.045, w - 0.045])
    for (const z of [-0.045, -d + 0.045]) {
      // Each leg runs from inside its foot to inside the top: no face of one lies on a face of the other.
      add(leg, p.ss, x, 0.025 + (h - 0.05) / 2, z);
      add(foot, p.ssDull, x, 0.015, z);
    }
  add(new THREE.BoxGeometry(w - 0.06, 0.018, d - 0.06), p.ss, w / 2, 0.24, -d / 2);
  return g;
}

/** A kraft takeout bag, folded over once, its base centre at the origin. */
function bagTemplate(): THREE.Group {
  const kraft = new THREE.MeshStandardMaterial({ color: 0xb48c5c, roughness: 0.88 });
  const fold = new THREE.MeshStandardMaterial({ color: 0x9c774b, roughness: 0.9 });
  const g = new THREE.Group();
  const body = new THREE.Mesh(new THREE.BoxGeometry(0.16, 0.19, 0.105), kraft);
  body.position.y = 0.095;
  const top = new THREE.Mesh(new THREE.BoxGeometry(0.162, 0.034, 0.03), fold);
  top.position.set(0, 0.203, 0.012);
  const seal = new THREE.Mesh(new THREE.BoxGeometry(0.05, 0.03, 0.003), new THREE.MeshStandardMaterial({ color: 0xf2b134, roughness: 0.6 }));
  seal.position.set(0, 0.2, 0.0285);
  for (const m of [body, top, seal]) {
    m.castShadow = m.receiveShadow = true;
    g.add(m);
  }
  return g;
}

/** A 2.5-inch half pan of a food at a fill, its base centre at the origin, its length along z (as a cabinet holds it). */
function panTemplate(food: number, fill: number): THREE.Group {
  const m = buildFixture(PAN.id, PAN, 1);
  m.hooks.set?.("food", food);
  m.hooks.set?.("fill", fill);
  m.group.position.set(-PAN.w / 2, 0, PAN.d / 2);
  const turn = new THREE.Group();
  turn.rotation.y = Math.PI / 2;
  turn.add(m.group);
  const g = new THREE.Group();
  g.add(turn);
  return g;
}

/* ------------------------------------------------------------------------------------------------------------ */
/* The cooks' day: per ten-minute stretch, the computed minutes loading the station and plating at its landing.  */
/* ------------------------------------------------------------------------------------------------------------ */

interface CookRig {
  c: (typeof COOKS)[number];
  st: Staff;
  routine: Routine | null;
}

/** Each cook faces their station all service: loading it, plating at its landing, or stepped back from it, idle. */
function cookSteps(day: DayRun, i: number, rig: CookRig, stepS: number, endMin: number): RoutineStep[] {
  const { c } = rig;
  type Seg = { kind: "load" | "plate" | "idle"; dur: number; busy: number };
  const segs: Seg[] = [];
  const push = (kind: Seg["kind"], dur: number, busy: number) => {
    if (dur <= 0) return;
    const prev = segs[segs.length - 1];
    if (prev && prev.kind === kind) {
      prev.busy = (prev.busy * prev.dur + busy * dur) / (prev.dur + dur);
      prev.dur += dur;
    } else segs.push({ kind, dur, busy });
  };
  for (let b = 0; b * BIN < endMin; b++) {
    let busy = 0;
    let plate = 0;
    const m1 = Math.min(endMin, (b + 1) * BIN);
    for (let m = b * BIN; m < m1; m++) {
      busy += day.busy[i]![m] ?? 0;
      plate += day.plating[i]![m] ?? 0;
    }
    const share = busy / (m1 - b * BIN);
    push("load", (busy - plate) * SEC_PER_MIN, share);
    push("plate", plate * SEC_PER_MIN, share);
    push("idle", (m1 - b * BIN - busy) * SEC_PER_MIN, share);
  }
  // Slivers go to the stretch before them, so no motion flickers in and out.
  for (let k = segs.length - 1; k > 0; k--)
    if (segs[k]!.dur < 0.35) {
      segs[k - 1]!.dur += segs[k]!.dur;
      segs.splice(k, 1);
    }
  const station = new THREE.Vector3(c.x + c.at[0]!, c.at[1]!, c.at[2]!);
  const landing = new THREE.Vector3(c.x + c.land[0]!, c.land[1]!, c.land[2]!);
  const steps: RoutineStep[] = [];
  for (const g of segs) {
    const tempo = 0.9 + 0.45 * Math.min(1, g.busy);
    if (g.kind === "load") steps.push({ kind: "work", work: c.work, duration: g.dur, opts: { at: station, surfaceY: c.surf, tempo } });
    else if (g.kind === "plate") steps.push({ kind: "work", work: "plate", duration: g.dur, opts: { at: landing, surfaceY: c.land[1]! - 0.004, tempo } });
    else if (g.dur > 2 * stepS + 1.2) {
      // A long wait: a step back from the station, and back up to it when the work comes.
      steps.push({ kind: "shuffle", to: [c.x, STATION_Z + 0.16] }, { kind: "stand", duration: g.dur - 2 * stepS }, { kind: "shuffle", to: [c.x, STATION_Z] });
    } else steps.push({ kind: "stand", duration: g.dur });
  }
  steps.push({ kind: "stand", duration: END_HOLD_S });
  return steps;
}
/* ------------------------------------------------------------------------------------------------------------ */

export function mountRushLab(root: HTMLElement): void {
  const mountT0 = performance.now();
  const stageEl = root.querySelector<HTMLElement>("[data-rush-stage]");
  const canvas = root.querySelector<HTMLCanvasElement>("[data-rush-canvas]");
  if (!stageEl || !canvas) return;
  const maybeRenderer = createStudioRenderer(canvas, false);
  if (!maybeRenderer) return;
  const renderer = maybeRenderer;
  renderer.toneMappingExposure = 1.0;
  const phone = window.matchMedia("(max-width: 720px)").matches;

  const scene = new THREE.Scene();
  scene.background = backdrop();
  scene.fog = new THREE.Fog(0x151510, 8, 16);
  scene.environment = studioEnvironment(renderer);
  scene.environmentIntensity = 0.6;
  const lights = studioLights(scene, new THREE.Vector3(-0.6, 0.9, -0.2), 4.2, phone ? 1024 : 2048);
  lights.key.position.set(-1.8, 6.2, 4.4);
  lights.key.intensity = 1.9;
  lights.key.shadow.normalBias = 0.02;
  lights.rim.position.set(3.6, 3.2, -3.0);
  lights.rim.intensity = 0.75;
  lights.hemi.intensity = 0.5;

  const floor = tilePlane(9, 6);
  floor.rotation.x = -Math.PI / 2;
  floor.position.set(-0.6, 0, -0.3);
  scene.add(floor);

  /* The line: holding cabinet, fryer, griddle on its stand, combi on its worktable. */
  const cab = buildFixture(CAB.id, CAB, 2);
  scene.add(place(cab, CAB.x0, LINE_Z, true));
  cab.hooks.open?.(0.92);
  const fryer = buildFixture(FRYER.id, FRYER, 2);
  scene.add(place(fryer, FRYER.x0, LINE_Z, true));
  const stand = equipmentStand(GRIDDLE.w, GRIDDLE.d - 0.04, STAND_H);
  stand.position.set(GRIDDLE.x0, 0, LINE_Z - 0.02);
  scene.add(stand);
  const griddle = buildFixture(GRIDDLE.id, GRIDDLE, 2);
  scene.add(place(griddle, GRIDDLE.x0, LINE_Z, true));
  griddle.group.position.y = STAND_H + 0.0005;
  const table = buildFixture(TABLE.id, TABLE, 2);
  scene.add(place(table, TABLE.x0, LINE_Z, true));
  const combi = buildFixture(COMBI.id, COMBI, 2);
  scene.add(place(combi, COMBI.x0, LINE_Z - 0.0, true));
  combi.group.position.y = TABLE.h + 0.0005;

  /* The expo pass: the table, its overshelf on two posts, the strip warmers under it, the ticket rail on its front. */
  const lamps: Model[] = [];
  const pass = buildFixture(PASS.id, PASS, 2);
  scene.add(place(pass, PASS.x0, PASS_Z, false));
  {
    const p = palette();
    const shelf = new THREE.Mesh(new THREE.BoxGeometry(PASS.w - 0.12, 0.014, SHELF_D), p.ss);
    shelf.position.set(PASS.x0 + PASS.w / 2, SHELF_Y + 0.007, SHELF_Z);
    shelf.castShadow = shelf.receiveShadow = true;
    scene.add(shelf);
    const post = new THREE.CylinderGeometry(0.0127, 0.0127, SHELF_Y + 0.01 - PASS.h + 0.006, 16);
    for (const x of [PASS.x0 + 0.1, PASS.x0 + PASS.w - 0.1]) {
      const m = new THREE.Mesh(post, p.ss);
      m.position.set(x, (SHELF_Y + 0.01 + PASS.h - 0.006) / 2, SHELF_Z);
      m.castShadow = true;
      scene.add(m);
    }
    for (const x0 of [PASS.x0, PASS.x0 + PASS.w - LAMP.w]) {
      const lamp = buildFixture(LAMP.id, LAMP, 1);
      scene.add(place(lamp, x0, SHELF_Z + LAMP.d / 2, true));
      lamp.group.position.y = SHELF_Y - LAMP.h - 0.003;
      lamp.hooks.set?.("heat", 0.85);
      lamps.push(lamp);
    }
    const rail = new THREE.Mesh(new THREE.BoxGeometry(PASS.w - 0.26, 0.03, 0.016), p.alu);
    rail.position.set(PASS.x0 + PASS.w / 2, SHELF_Y - 0.016, SHELF_Z + SHELF_D / 2 - 0.008);
    rail.castShadow = true;
    scene.add(rail);
  }
  const railX0 = PASS.x0 + 0.17 + TICKET.w / 2;
  const railZ = SHELF_Z + SHELF_D / 2 + 0.002;

  /* The tickets: one instanced card each, coloured by how long its guest has waited. */
  const ticketMat = new THREE.MeshStandardMaterial({ color: 0xffffff, roughness: 0.92, side: THREE.DoubleSide });
  const tickets = new THREE.InstancedMesh(new THREE.PlaneGeometry(TICKET.w, TICKET.h), ticketMat, TICKET.perRow * 2);
  tickets.count = 0;
  tickets.castShadow = true;
  tickets.frustumCulled = false;
  tickets.setColorAt(0, COLOR.paper);
  scene.add(tickets);

  /* The cabinet's pans: per food, three fills (full, part, last scrapings). */
  const panGroup = new THREE.Group();
  scene.add(panGroup);
  const FILLS = [0.85, 0.55, 0.28];
  const panSets = [0, 1, 2].map((food) => FILLS.map((f) => new Instances(panTemplate(food, f), SERVICE.cabinetPans, panGroup)));
  // Rails as the cabinet builds them (models/ovens.ts holdingCabinet), less the lowest, over its water pan.
  const cavY0 = 0.17 + 0.065;
  const rails: number[] = [];
  for (let y = cavY0 + 0.1; y < CAB.h - 0.178 - 0.03 - 0.06; y += 0.0762) rails.push(y);
  const panSpots: THREE.Vector3[] = [];
  // Filled from eye height outwards, as a cook loads the rails easiest to reach.
  const order = rails.slice(1).sort((a, b) => Math.abs(a - 1.15) - Math.abs(b - 1.15));
  for (const y of order) for (const sx of [-1, 1]) panSpots.push(new THREE.Vector3(CAB.x0 + CAB.w / 2 + sx * 0.134, y - 0.06, LINE_Z - 0.3875));

  /* The bags on the pass. */
  const bags = new Instances(bagTemplate(), 22, scene);
  const bagSpot = (k: number) => new THREE.Vector3(PASS.x0 + 0.26 + (k % 11) * 0.2, PASS.h + 0.0005, PASS_Z + 0.2 + Math.floor(k / 11) * 0.16);

  /* The cooks. */
  const cooks: CookRig[] = COOKS.map((c, i) => {
    const st = makeStaff(c.cohort, 2, { seed: c.seed, apron: c.apron, sleeves: i === 1 ? "rolled" : "cuffed" });
    st.place(c.x, STATION_Z, Math.PI);
    scene.add(st.group);
    return { c, st, routine: null };
  });

  /* Camera: from the guests' side, raised over the pass. */
  const camera = new THREE.PerspectiveCamera(32, 16 / 10, 0.05, 60);
  const TARGET = new THREE.Vector3();
  const VIEW = new THREE.Vector3();
  const keep: THREE.Vector3[] = [];
  const aim = () => {
    keep.length = 0;
    const wide = stageEl.clientWidth >= 1.3 * stageEl.clientHeight;
    if (wide) {
      TARGET.set(-0.5, 0.92, -0.15);
      VIEW.set(0.3, 0.66, 1).normalize();
      for (const x of [CAB.x0, TABLE.x0 + TABLE.w]) keep.push(new THREE.Vector3(x, 1.8, LINE_Z - 0.3), new THREE.Vector3(x, 0.75, PASS_Z + PASS.d));
    } else {
      // From the combi's end: the line runs away to the cabinet.
      TARGET.set(-0.6, 0.95, -0.2);
      VIEW.set(0.62, 0.7, 1).normalize();
      for (const x of [CAB.x0 + 0.1, TABLE.x0 + TABLE.w]) for (const y of [0.1, 1.85]) for (const z of [LINE_Z, PASS_Z + PASS.d]) keep.push(new THREE.Vector3(x, y, z));
    }
  };
  const controls = createScrollSafeOrbit(camera, canvas, [TARGET.x, TARGET.y, TARGET.z], { minDistance: 1.2, maxDistance: 14, maxPolarAngle: Math.PI * 0.47 });
  let userMoved = false;
  controls.addEventListener("start", () => {
    userMoved = true;
  });
  const ndc = new THREE.Vector3();
  const frame = () => {
    if (userMoved) return;
    aim();
    camera.aspect = Math.max(0.2, stageEl.clientWidth / Math.max(1, stageEl.clientHeight));
    camera.updateProjectionMatrix();
    let dist = 2.5;
    for (; dist < 20; dist += 0.06) {
      camera.position.copy(TARGET).addScaledVector(VIEW, dist);
      camera.lookAt(TARGET);
      camera.updateMatrixWorld(true);
      if (keep.every((p) => (ndc.copy(p).project(camera), Math.abs(ndc.x) < 0.95 && Math.abs(ndc.y) < 0.86))) break;
    }
    const fog = scene.fog as THREE.Fog;
    fog.near = dist * 1.1;
    fog.far = dist * 2.6;
    controls.target.copy(TARGET);
    controls.update();
  };
  const gate = shaderGate(renderer, camera);
  bindCanvasResize(stageEl, renderer, camera, gate);
  frame();
  new ResizeObserver(() => frame()).observe(stageEl);

  /* ---------------------------------------------------------------------------------------------------------- */
  /* The page: day, menu (what is held), crew, the thighs' station, forecast error; the chart; the tiles.        */
  /* ---------------------------------------------------------------------------------------------------------- */

  const q = <T extends Element = HTMLElement>(s: string) => root.querySelector<T>(s);
  const planBtns = Array.from(root.querySelectorAll<HTMLButtonElement>("[data-rs-plan]"));
  const dayBtns = Array.from(root.querySelectorAll<HTMLButtonElement>("[data-rs-day]"));
  const crewBtns = Array.from(root.querySelectorAll<HTMLButtonElement>("[data-rs-pool]"));
  const swapBtns = Array.from(root.querySelectorAll<HTMLButtonElement>("[data-rs-menu]"));
  const menuEl = q("[data-rs-dishes]");
  const errIn = q<HTMLInputElement>("[data-rs-err]");
  const errVal = q("[data-rs-err-val]");
  const statusEl = q("[data-rs-status]");
  const chartEl = q("[data-rush-chart]");
  const hudClock = q("[data-rs-hud-clock]");
  const hudTickets = q("[data-rs-hud-tickets]");
  const hudNote = q("[data-rs-hud-note]");
  const tile = (k: string) => ({ v: q(`[data-rs-${k}]`), bad: q(`[data-rs-${k}-p90]`), delta: q(`[data-rs-${k}-d]`) });
  const tiles = { wait: tile("wait"), long: tile("long"), rush: tile("rush"), waste: tile("waste") };

  const settings: RushSettings = { modes: modesFor([]), day: "busy", sigma: Number(errIn?.value ?? 10) / 100, menu: "written", pool: false };

  /* The menu card: each dish, how long it holds, and whether it is held or cooked to order. */
  const dishBtns = new Map<string, HTMLButtonElement>();
  if (menuEl) {
    menuEl.innerHTML = "";
    for (const dish of RUSH_MENU) {
      const li = document.createElement("li");
      const b = document.createElement("button");
      b.type = "button";
      b.className = "wk-dish rs-dish";
      b.dataset.rsDish = dish.key;
      b.addEventListener("click", () => {
        settings.modes = { ...settings.modes, [dish.key]: settings.modes[dish.key] === "hold" ? "order" : "hold" };
        update();
      });
      li.appendChild(b);
      menuEl.appendChild(li);
      dishBtns.set(dish.key, b);
    }
  }
  const holdBars = (min: number) => {
    const n = min >= 90 ? 4 : min >= 45 ? 3 : min >= 20 ? 2 : 1;
    return `<span class="wk-keeps" aria-hidden="true">${[1, 2, 3, 4].map((k) => `<i${k <= n ? ' class="on"' : ""}></i>`).join("")}<small>${min >= 60 ? `${min / 60} h` : `${min} min`}</small></span>`;
  };
  function renderMenu(): void {
    const menu = rushMenu(settings.menu ?? "written");
    for (const dish of menu) {
      const b = dishBtns.get(dish.key);
      if (!b) continue;
      const held = settings.modes[dish.key] === "hold";
      b.innerHTML = `<b>${dish.name}</b>${holdBars(dish.holdMin)}<em>${held ? "Held" : "To order"}</em>`;
      b.setAttribute("aria-pressed", String(held));
      b.title = held ? `Cook ${dish.name.toLowerCase()} to order` : `Cook ${dish.name.toLowerCase()} ahead and hold it (it keeps ${dish.holdMin} minutes)`;
    }
    const held = RUSH_MENU.filter((d) => settings.modes[d.key] === "hold").map((d) => d.key);
    const preset = held.length === 0 ? "order" : held.length === KEEPERS.length && KEEPERS.every((k) => held.includes(k)) ? "keepers" : "";
    for (const b of planBtns) b.setAttribute("aria-pressed", String(b.dataset.rsPlan === preset));
    for (const b of dayBtns) b.setAttribute("aria-pressed", String(b.dataset.rsDay === settings.day));
    for (const b of crewBtns) b.setAttribute("aria-pressed", String((b.dataset.rsPool === "yes") === !!settings.pool));
    for (const b of swapBtns) b.setAttribute("aria-pressed", String(b.dataset.rsMenu === settings.menu));
  }

  /* The chart: orders in, tickets waiting and portions held over the service; each cook's busy share beneath. */
  const X0 = 52;
  let X1 = 592;
  let W = 600;
  const TOP0 = 10;
  const TOP1 = 128;
  const ROW0 = 146;
  const ROW_H = 13;
  let domain = 300;
  const xOf = (m: number) => X0 + (Math.min(m, domain) / domain) * (X1 - X0);
  let cursor: SVGLineElement | null = null;
  function renderChart(d: DayRun, endMin: number): void {
    if (!chartEl) return;
    // Drawn at the panel's own pixel width, so the labels keep their size on a phone.
    W = Math.max(300, Math.round(chartEl.clientWidth - 2 * 12));
    X1 = W - 8;
    domain = Math.max(300, Math.ceil(endMin / 30) * 30);
    const n = Math.ceil(domain);
    const arrivals = new Float32Array(Math.ceil(domain / 5));
    const diff = new Float32Array(n + 2);
    for (const o of d.orders) {
      const a = Math.min(n, Math.ceil(o.arrive));
      arrivals[Math.min(arrivals.length - 1, Math.floor(o.arrive / 5))]!++;
      diff[a]!++;
      diff[Math.min(n + 1, Math.ceil(o.plated))]!--;
    }
    const waiting = new Float32Array(n);
    let w = 0;
    for (let m = 0; m < n; m++) waiting[m] = w += diff[m]!;
    const held = new Float32Array(n);
    for (const s of d.stock) for (let m = 0; m < Math.min(n, s.length); m++) held[m]! += s[m]!;
    const top = Math.max(10, ...waiting, ...held);
    const nice = [20, 40, 60, 80, 100, 120, 160, 200, 240, 300, 400, 500].find((v) => v >= top) ?? Math.ceil(top / 100) * 100;
    const yOf = (v: number) => TOP1 - (v / nice) * (TOP1 - TOP0);
    // The orders coming in are a rate, not a count on the axis: drawn as the day's shape, its peak labelled.
    const inTop = Math.max(1, ...arrivals);
    const yIn = (v: number) => TOP1 - (v / inTop) * 0.62 * (TOP1 - TOP0);
    const line = (a: Float32Array) => {
      let s = "";
      for (let m = 0; m < a.length; m++) s += `${m ? "L" : "M"}${xOf(m + 0.5).toFixed(1)} ${yOf(a[m]!).toFixed(1)}`;
      return s;
    };
    const H = ROW0 + 3 * ROW_H + 22;
    let svg = `<svg viewBox="0 0 ${W} ${H}" width="${W}" height="${H}" class="rs-svg" aria-hidden="true">`;
    // Grid and axis.
    for (const v of [0, nice / 2, nice]) svg += `<line class="g" x1="${X0}" x2="${X1}" y1="${yOf(v)}" y2="${yOf(v)}"/><text class="ya" x="${X0 - 6}" y="${yOf(v) + 3.5}">${v}</text>`;
    const step = W < 420 ? 120 : 60;
    for (let m = 0; m <= domain; m += step) svg += `<line class="g v" x1="${xOf(m)}" x2="${xOf(m)}" y1="${TOP0}" y2="${ROW0 + 3 * ROW_H}"/>${xOf(m) > X1 - 16 ? "" : `<text class="xa" x="${xOf(m)}" y="${ROW0 + 3 * ROW_H + 16}">${clock(m)}</text>`}`;
    // Doors open to last order.
    svg += `<rect class="open" x="${xOf(SERVICE.open)}" y="${TOP0}" width="${xOf(SERVICE.lastOrder) - xOf(SERVICE.open)}" height="${TOP1 - TOP0}"/>`;
    // Orders in, per 5 minutes.
    const bw = ((X1 - X0) * 5) / domain;
    arrivals.forEach((v, k) => {
      if (v > 0) svg += `<rect class="in" x="${(xOf(k * 5) + 0.6).toFixed(1)}" y="${yIn(v).toFixed(1)}" width="${(bw - 1.2).toFixed(1)}" height="${(TOP1 - yIn(v)).toFixed(1)}"/>`;
    });
    svg += `<path class="held" d="${line(held)}L${xOf(n - 0.5).toFixed(1)} ${TOP1}L${xOf(0.5).toFixed(1)} ${TOP1}Z"/>`;
    svg += `<path class="wait" d="${line(waiting)}"/>`;
    // The crew: one row per cook, a cell per 5 minutes, brighter the busier.
    const names = ["Fryer", "Griddle", "Combi"];
    d.busy.forEach((b, i) => {
      const y = ROW0 + i * ROW_H;
      svg += `<text class="ya" x="${X0 - 6}" y="${y + ROW_H / 2 + 3}">${names[i]}</text>`;
      for (let k = 0; k * 5 < domain; k++) {
        let s = 0;
        for (let m = k * 5; m < k * 5 + 5; m++) s += b[m] ?? 0;
        const f = s / 5;
        svg += `<rect class="cell" x="${(xOf(k * 5) + 0.4).toFixed(1)}" y="${y + 1}" width="${(bw - 0.8).toFixed(1)}" height="${ROW_H - 2}" style="opacity:${(0.06 + 0.94 * f).toFixed(2)}"/>`;
      }
    });
    svg += `<line class="cur" x1="${X0}" x2="${X0}" y1="${TOP0 - 4}" y2="${ROW0 + 3 * ROW_H + 2}"/></svg>`;
    chartEl.innerHTML =
      svg +
      `<p class="wk-keys rs-keys"><span class="in">Orders in (peak ${inTop} in 5 min)</span><span class="wait">Tickets on the rail</span><span class="held">Portions in the cabinet</span><span class="busy">Cook busy</span></p>`;
    cursor = chartEl.querySelector<SVGLineElement>(".cur");
  }

  function renderTiles(run: RushRun, base: RushRun | null): void {
    const set = (t: (typeof tiles)[keyof typeof tiles], v: string, bad: string, delta: number | null, unit: string, digits = 0, upGood = false) => {
      if (t.v) t.v.textContent = v;
      if (t.bad) t.bad.textContent = bad;
      if (t.delta) {
        t.delta.hidden = delta === null || Math.abs(delta) < 0.5 * 10 ** -digits;
        if (delta !== null) {
          t.delta.textContent = `${delta > 0 ? "+" : "−"}${Math.abs(delta).toFixed(Math.abs(delta) >= 10 ? 0 : digits)}${unit}`;
          t.delta.dataset.good = String(upGood ? delta > 0 : delta < 0);
        }
      }
    };
    set(tiles.wait, fmtMin(run.waitP50.p50), fmtMin(run.waitP50.p90), base ? run.waitP50.p50 - base.waitP50.p50 : null, " min", 1);
    set(tiles.long, fmtMin(run.waitP90.p50), fmtMin(run.waitP90.p90), base ? run.waitP90.p50 - base.waitP90.p50 : null, " min", 1);
    const pct = (x: number) => `${Math.round(x * 100)}%`;
    set(tiles.rush, pct(run.prepLoad.p50), pct(run.rushLoad.p50), base ? Math.round((run.prepLoad.p50 - base.prepLoad.p50) * 100) : null, " pts", 0, true);
    set(tiles.waste, `${Math.round(run.waste.p50)}`, `${Math.round(run.waste.p90)}`, base ? Math.round(run.waste.p50 - base.waste.p50) : null, "");
  }

  function renderStatus(run: RushRun): void {
    if (!statusEl) return;
    const held = RUSH_MENU.filter((d) => settings.modes[d.key] === "hold");
    const p90 = run.waitP90.p50;
    // The bottleneck: the cook busiest over the shown day's rush hour.
    const names = ["fryer", "griddle", "combi"];
    const load = day.busy.map((b) => {
      let best = 0;
      let sum = 0;
      for (let m = 0; m < b.length; m++) {
        sum += b[m]! - (m >= 60 ? b[m - 60]! : 0);
        best = Math.max(best, sum);
      }
      return best;
    });
    const hot = load.indexOf(Math.max(...load));
    const hint = !settings.pool ? "Let whoever is free plate." : settings.menu === "written" && hot === 1 ? "Roast the thighs in the combi." : "";
    const short = held.filter((d) => d.holdMin < 30);
    let text: string;
    if (!held.length) text = `Nothing gets cooked before 11:00, then the line falls behind: one guest in ten waits ${fmtMin(p90)}. Hold the dishes that keep.`;
    else if (p90 > 12) text = `The cabinet covers ${held.length === 1 ? held[0]!.name.toLowerCase() : `${held.length} dishes`}, but the ${names[hot]} still holds the line up: one guest in ten waits ${fmtMin(p90)}. ${hint}`;
    else text = `One guest in ten waits ${fmtMin(p90)}, for ${Math.round(run.waste.p50)} portions thrown out.`;
    if (short.length) text += ` ${short.map((d) => d.name).join(" and ")} ${short.length > 1 ? "don't" : "doesn't"} hold long enough: most of the waste is ${short.length > 1 ? "theirs" : "its"}.`;
    statusEl.textContent = text.trim();
    statusEl.dataset.state = p90 <= 12 ? "ok" : p90 <= 25 ? "warn" : "danger";
  }

  /* ---------------------------------------------------------------------------------------------------------- */
  /* Running it.                                                                                                 */
  /* ---------------------------------------------------------------------------------------------------------- */

  const cache = new Map<string, RushRun>();
  const runFor = (s: RushSettings) => {
    const k = JSON.stringify([s.day, s.sigma, s.menu, s.pool, RUSH_MENU.map((d) => s.modes[d.key] ?? "order")]);
    let r = cache.get(k);
    if (!r) {
      r = runRush(s, DAYS_N);
      cache.set(k, r);
    }
    return r;
  };
  let run: RushRun = runFor(settings);
  let day: DayRun = run.shown;
  let endMin: number = SERVICE.close;
  let loopS = 1;
  let tMin = 0;
  let loopT = 0;
  let stepS = 0;

  function buildRoutines(): void {
    // How long a step back from the station takes, measured once.
    if (cooks[0] && stepS === 0) {
      const c = cooks[0];
      stepS = new Routine(c.st, { x: c.c.x, z: STATION_Z, yaw: Math.PI }, [{ kind: "shuffle", to: [c.c.x, STATION_Z + 0.16] }, { kind: "shuffle", to: [c.c.x, STATION_Z] }], false).duration / 2;
    }
    cooks.forEach((k, i) => {
      const steps = cookSteps(day, i, k, stepS, endMin);
      let r = new Routine(k.st, { x: k.c.x, z: STATION_Z, yaw: Math.PI }, steps, true);
      // Square the loop with the clock: pad (or trim) the closing pause.
      const err = loopS - r.duration;
      const last = steps[steps.length - 1];
      if (Math.abs(err) > 0.01 && last && last.kind === "stand") {
        last.duration = Math.max(0.1, last.duration + err);
        r = new Routine(k.st, { x: k.c.x, z: STATION_Z, yaw: Math.PI }, steps, true);
      }
      k.routine = r;
    });
  }

  let pending = 0;
  function update(): void {
    const frac = loopS > 0 ? tMin / Math.max(1, endMin) : 0;
    run = runFor(settings);
    const anyHeld = Object.values(settings.modes).some((m: Mode) => m === "hold");
    const base = anyHeld || settings.pool || settings.menu === "balanced" ? runFor({ ...settings, modes: modesFor([]), pool: false, menu: "written" }) : null;
    day = run.shown;
    let last: number = SERVICE.close;
    for (const o of day.orders) last = Math.max(last, o.plated);
    endMin = Math.min(SERVICE.end, Math.ceil(last) + 4);
    loopS = endMin * SEC_PER_MIN + END_HOLD_S;
    buildRoutines();
    tMin = frac * endMin;
    loopT = tMin * SEC_PER_MIN;
    renderMenu();
    renderChart(day, endMin);
    renderTiles(run, base);
    renderStatus(run);
    if (errVal && settings.sigma !== undefined) errVal.textContent = settings.sigma > 0 ? `±${Math.round(settings.sigma * 100)}%` : "Exact";
    lastMinute = -1;
  }
  const later = () => {
    if (pending) return;
    pending = window.setTimeout(() => {
      pending = 0;
      update();
    }, 90);
  };

  for (const b of planBtns)
    b.addEventListener("click", () => {
      settings.modes = modesFor(b.dataset.rsPlan === "keepers" ? [...KEEPERS] : []);
      update();
    });
  for (const b of dayBtns)
    b.addEventListener("click", () => {
      settings.day = b.dataset.rsDay as RushDay;
      update();
    });
  for (const b of crewBtns)
    b.addEventListener("click", () => {
      settings.pool = b.dataset.rsPool === "yes";
      update();
    });
  for (const b of swapBtns)
    b.addEventListener("click", () => {
      settings.menu = b.dataset.rsMenu === "balanced" ? "balanced" : "written";
      update();
    });
  errIn?.addEventListener("input", () => {
    settings.sigma = Number(errIn.value) / 100;
    if (errVal) errVal.textContent = settings.sigma > 0 ? `±${Math.round(settings.sigma * 100)}%` : "Exact";
    later();
  });

  /* One moment of the service. */
  const M = new THREE.Matrix4();
  const Q = new THREE.Quaternion();
  const ONE = new THREE.Vector3(1, 1, 1);
  const P = new THREE.Vector3();
  let lastMinute = -1;
  const lift = [1, 1];
  function show(t: number, dt: number): void {
    const m = Math.min(endMin - 1, Math.floor(t));
    // Stations: by the batches each has cooking this minute.
    const cookingAt = (s: number) => day.cooking[s]![Math.max(0, m)] ?? 0;
    const fr = cookingAt(0);
    const target = [fr > 0.5 ? 0 : 1, fr > 1.5 ? 0 : 1];
    for (let k = 0; k < 2; k++) lift[k]! += (target[k]! - lift[k]!) * Math.min(1, dt * 3);
    fryer.hooks.set?.("lift_l", lift[0]!);
    fryer.hooks.set?.("lift_r", lift[1]!);
    fryer.hooks.tick?.(t, dt, Math.min(1, fr / 2));
    const gr = Math.min(1, cookingAt(1) / 2);
    griddle.hooks.set?.("heat", 0.4 + 0.6 * gr);
    griddle.hooks.tick?.(t, dt, gr);
    const co = Math.min(1, cookingAt(2) / 2);
    combi.hooks.set?.("heat", co);
    combi.hooks.set?.("fan", co > 0 ? 0.4 + 0.6 * co : 0);
    combi.hooks.set?.("steam", 0.5 * co);
    combi.hooks.tick?.(t, dt, co);
    cab.hooks.set?.("heat", 0.8);
    cab.hooks.set?.("humidity", 0.5);
    cab.hooks.tick?.(t, dt, 0.6);
    for (const l of lamps) l.hooks.tick?.(t, dt, 0.85);

    if (m !== lastMinute) {
      lastMinute = m;
      // Tickets on the rail: orders in and not yet out, oldest at the left.
      let n = 0;
      let waiting = 0;
      for (const o of day.orders) {
        if (o.arrive > t) break;
        if (o.plated <= t) continue;
        waiting++;
        if (n >= TICKET.perRow * 2) continue;
        const row = Math.floor(n / TICKET.perRow);
        const col = n % TICKET.perRow;
        P.set(railX0 + col * TICKET.step, SHELF_Y - 0.026 - TICKET.h / 2 - row * (TICKET.h + 0.012), railZ + (col % 2) * 0.0025);
        M.compose(P, Q.identity(), ONE);
        tickets.setMatrixAt(n, M);
        const age = t - o.arrive;
        tickets.setColorAt(n, age >= AGE_RED ? COLOR.red : age >= AGE_AMBER ? COLOR.amber : COLOR.paper);
        n++;
      }
      tickets.count = n;
      tickets.instanceMatrix.needsUpdate = true;
      if (tickets.instanceColor) tickets.instanceColor.needsUpdate = true;

      // The cabinet: each held dish's portions as half pans, full ones first.
      for (const set of panSets) for (const s of set) s.begin();
      let spot = 0;
      let heldNow = 0;
      const menu = rushMenu(settings.menu ?? "written");
      menu.forEach((dish, d) => {
        const s = Math.round(day.stock[d]![Math.max(0, m)] ?? 0);
        heldNow += s;
        let left = s;
        while (left > 0 && spot < panSpots.length) {
          const f = Math.min(1, left / dish.perPan);
          const fi = f > 0.7 ? 0 : f > 0.4 ? 1 : 2;
          M.compose(panSpots[spot++]!, Q.identity(), ONE);
          panSets[dish.food === 3 ? 0 : Math.min(2, dish.food)]![fi]!.push(M);
          left -= dish.perPan;
        }
      });
      for (const set of panSets) for (const s of set) s.end();

      // Bags out on the pass: orders handed over in the last two minutes.
      bags.begin();
      let nb = 0;
      for (const o of day.orders) {
        if (o.arrive > t) break;
        if (o.plated <= t && o.plated > t - 2 && nb < 22) {
          M.compose(bagSpot(nb), Q.setFromAxisAngle(THREE.Object3D.DEFAULT_UP, ((nb * 37) % 11) * 0.02 - 0.1), ONE);
          bags.push(M);
          nb++;
        }
      }
      bags.end();

      if (hudClock) hudClock.textContent = clock(t);
      if (hudTickets) {
        hudTickets.textContent = waiting ? `${waiting} ${waiting === 1 ? "ticket" : "tickets"}` : "Rail clear";
        hudTickets.dataset.kind = waiting >= 25 ? "late" : "";
      }
      if (hudNote) {
        hudNote.textContent = heldNow ? `${heldNow} portions held` : t < SERVICE.open ? "Doors at 11:00" : t >= SERVICE.lastOrder ? "Last orders in" : "Nothing held";
        hudNote.dataset.kind = heldNow ? "held" : "idle";
      }
    }
    if (cursor) {
      const x = xOf(t).toFixed(1);
      cursor.setAttribute("x1", x);
      cursor.setAttribute("x2", x);
    }
  }

  update();
  show(0, 0);
  if (chartEl) {
    let lastW = chartEl.clientWidth;
    new ResizeObserver(() => {
      if (Math.abs(chartEl.clientWidth - lastW) < 4) return;
      lastW = chartEl.clientWidth;
      renderChart(day, endMin);
    }).observe(chartEl);
  }
  for (const k of cooks) void k.st.ready.then(() => void gate.refresh(scene));

  let inView = false;
  observeInView(stageEl, (v) => {
    inView = v;
  });
  let frozen: number | null = null;
  let last = performance.now();
  const fitClip = clipFitter(camera, scene, 0.05);
  function animate(now: number): void {
    requestAnimationFrame(animate);
    const dt = Math.max(0, Math.min(0.05, (now - last) / 1000));
    last = now;
    if (!inView || document.hidden || !gate.ready()) return;
    if (frozen === null) loopT = (loopT + dt) % loopS;
    else loopT = frozen * SEC_PER_MIN;
    tMin = Math.min(endMin, loopT / SEC_PER_MIN);
    for (const k of cooks) k.routine?.apply(loopT);
    show(tMin, dt);
    controls.update();
    fitClip();
    renderer.render(scene, camera);
  }
  requestAnimationFrame(animate);
  void gate.refresh(scene);

  // A hook for checking the page by script.
  Object.assign(stageEl, {
    rush: {
      info: () => ({
        mountMs: Math.round(mountMs),
        tris: renderer.info.render.triangles,
        calls: renderer.info.render.calls,
        loopS: +loopS.toFixed(2),
        stepS: +stepS.toFixed(2),
        routineS: cooks.map((k) => +(k.routine?.duration ?? 0).toFixed(2)),
        t: +tMin.toFixed(1),
        tickets: tickets.count,
        acts: cooks.map((k) => {
          const st = k.routine?.apply(loopT);
          return st ? `${st.activity}${st.work ? `:${st.work}` : ""}` : "none";
        }),
        busyNow: day.busy.map((b) => +(b[Math.floor(tMin)] ?? 0).toFixed(2)),
        waitP90: +run.waitP90.p50.toFixed(1),
        waste: run.waste.p50,
      }),
      hold(min: number | null) {
        frozen = min;
      },
      view(eye: [number, number, number], at: [number, number, number]) {
        userMoved = true;
        camera.position.set(...eye);
        controls.target.set(...at);
        controls.update();
      },
    },
  });
  const mountMs = performance.now() - mountT0;
}
