/**
 * Every step at the station (chapter 05, second figure). From the burner to the plate should be a pivot, not a
 * pivot, step, step, drop. Two stations side by side, the same cook and the same plate at each: on the left the
 * plating table as built, across a 60-inch aisle and along it; on the right where you move it, straight behind across a
 * 36-inch aisle or beside the range. Each cook sautes, turns to the table, plates and goes back, on a routine the gait
 * planner turns into footsteps (Routine); the figure counts those footsteps, lays them on the floor, and adds them up
 * over a service and a year.
 */
import * as THREE from "three";
import { bindCanvasResize, clipFitter, createScrollSafeOrbit, createStudioRenderer, observeInView, shaderGate } from "./kit";
import { makePlate, makeTool, type ToolKind } from "./heroKitchen";
import { buildFixture, type Model, type Piece } from "./models/all";
import { Routine, makeStaff, type RoutineStep, type Staff, type WorkKind } from "./models/people";
import { studioEnvironment, studioLights } from "./models/stage";
import { backdrop, place, tilePlane } from "./weekLab";

const RANGE = { id: "NSF-RANGE-GAS-6B", w: 0.914, d: 0.864, h: 1.473 };
const TABLE = { id: "NSF-PREPTABLE-48", w: 1.219, d: 0.762, h: 0.914 };
const SKILLET = { id: "NSF-SKILLET-10IN", w: 0.254, d: 0.42, h: 0.05 };
const GRATE_Y = 0.914;
/** The two stations' range centres, as built on the left and moved on the right. */
const STATION_X = [-1.5, 1.5] as const;
/** The front middle burner, from the range's centre and front. */
const BURNER = { x: 0, z: -0.22 };
/** Seconds of sauteing and of plating a plate, the same at both stations (shortened; the moves between are real). */
const SAUTE_S = 4.2;
const PLATE_S = 3.2;
/** Services a year, for the year's walking: six a week, fifty weeks. */
const SERVICES_A_YEAR = 300;

export type PivotMove = "behind" | "beside";
export const PIVOT_PIECES: readonly Piece[] = [
  [RANGE.id, RANGE, 2],
  [RANGE.id, RANGE, 2],
  [TABLE.id, TABLE, 2],
  [TABLE.id, TABLE, 2],
  [SKILLET.id, SKILLET, 1],
  [SKILLET.id, SKILLET, 1],
];

interface Spot {
  x: number;
  z: number;
  yaw: number;
}
interface Layout {
  name: string;
  /** Where the cook sautes and where they plate, the plate's place on the table, and how the table stands. */
  burner: Spot;
  plating: Spot;
  plate: THREE.Vector3;
  table: { x: number; z: number; yaw: number };
}

/** A station's layout, from its range centre `c`. */
function layoutFor(c: number, move: PivotMove | "built"): Layout {
  if (move === "built") {
    // Across a 60-inch aisle and down the line, its near end 0.55 m along from the range's centre: turn, walk, turn.
    const front = 1.52;
    const x0 = c + 0.55;
    return {
      name: "As built",
      burner: { x: c, z: 0.3, yaw: Math.PI },
      plating: { x: x0 + 0.36, z: front - 0.3, yaw: 0 },
      plate: new THREE.Vector3(x0 + 0.36, TABLE.h, front + 0.24),
      table: { x: x0, z: front, yaw: Math.PI },
    };
  }
  if (move === "behind") {
    // Straight behind across a 36-inch aisle: turn round and a short step.
    const front = 0.914;
    return {
      name: "Behind, 36 in",
      burner: { x: c, z: 0.3, yaw: Math.PI },
      plating: { x: c, z: front - 0.3, yaw: 0 },
      plate: new THREE.Vector3(c, TABLE.h, front + 0.24),
      table: { x: c - TABLE.w / 2, z: front, yaw: Math.PI },
    };
  }
  // Beside the range, its front facing the cook's right shoulder: a quarter turn and back.
  const fx = c + RANGE.w / 2 + 0.03;
  return {
    name: "Beside the range",
    burner: { x: c + 0.14, z: 0.3, yaw: Math.PI },
    plating: { x: c + 0.14, z: 0.3, yaw: Math.PI / 2 },
    plate: new THREE.Vector3(fx + 0.26, TABLE.h, 0.3),
    table: { x: fx, z: -0.3, yaw: -Math.PI / 2 },
  };
}

/** How the cook gets from one place to the next: a turn on the spot, a turn and a short step, or a turn, a walk and a
 * turn. The walk's first turn, to face its way, is its own step (the walk would turn on the spot there itself), so its
 * footfalls count as a pivot's. */
function moveSteps(a: Spot, b: Spot): RoutineStep[] {
  const d = Math.hypot(b.x - a.x, b.z - a.z);
  if (d < 0.02) return [{ kind: "turn", yaw: b.yaw }];
  if (d < 0.45) return [{ kind: "turn", yaw: b.yaw }, { kind: "shuffle", to: [b.x, b.z] }];
  return [
    { kind: "turn", yaw: Math.atan2(b.x - a.x, b.z - a.z) },
    { kind: "walk", path: [[b.x, b.z]], speed: 1.0 },
    { kind: "turn", yaw: b.yaw },
  ];
}

/** A footstep of the loop, and whether it walked the cook somewhere or turned them on the spot. */
interface Footfall {
  side: "L" | "R";
  x: number;
  z: number;
  yaw: number;
  tOn: number;
  pivot: boolean;
}

/** What one plate costs at a layout: walking steps, pivots (turns on the spot), seconds on the move, metres walked. */
export interface PlateCost {
  steps: number;
  pivots: number;
  moveS: number;
  metres: number;
}

interface Station {
  c: number;
  st: Staff;
  routine: Routine | null;
  steps: Footfall[];
  cost: PlateCost;
  prints: THREE.InstancedMesh;
  path: THREE.Mesh;
  table: Model;
  plate: THREE.Object3D;
  skillet: THREE.Group;
  range: Model;
  tools: Map<ToolKind, THREE.Object3D>;
  held: ToolKind | null;
  layout: Layout;
}

const TOOL_FOR: Partial<Record<WorkKind, ToolKind>> = { saute: "skillet", plate: "spoon" };

/** A flat strip `width` wide along the floor through points (x, z), at height `y`: the path the body takes. */
function ribbon(pts: readonly THREE.Vector2[], width: number, y: number): THREE.BufferGeometry {
  const pos: number[] = [];
  const idx: number[] = [];
  const n = pts.length;
  for (let i = 0; i < n; i++) {
    const a = pts[Math.max(0, i - 1)]!;
    const b = pts[Math.min(n - 1, i + 1)]!;
    let dx = b.x - a.x;
    let dz = b.y - a.y;
    const l = Math.hypot(dx, dz) || 1;
    dx /= l;
    dz /= l;
    const p = pts[i]!;
    pos.push(p.x - dz * width * 0.5, y, p.y + dx * width * 0.5, p.x + dz * width * 0.5, y, p.y - dx * width * 0.5);
    if (i + 1 < n) idx.push(2 * i, 2 * i + 2, 2 * i + 1, 2 * i + 1, 2 * i + 2, 2 * i + 3);
  }
  const g = new THREE.BufferGeometry();
  g.setAttribute("position", new THREE.Float32BufferAttribute(pos, 3));
  g.setIndex(idx);
  g.computeVertexNormals();
  return g;
}
const MAX_PRINTS = 48;

export function mountPivotLab(root: HTMLElement): void {
  const mountT0 = performance.now();
  const stageEl = root.querySelector<HTMLElement>("[data-pivot-stage]");
  const canvas = root.querySelector<HTMLCanvasElement>("[data-pivot-canvas]");
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
  const lights = studioLights(scene, new THREE.Vector3(0, 0.9, 0.3), 4.8, phone ? 1024 : 2048);
  lights.key.position.set(-1.0, 6.4, 4.8);
  lights.key.intensity = 1.9;
  lights.key.shadow.normalBias = 0.02;
  lights.rim.position.set(3.2, 3.2, -3.0);
  lights.rim.intensity = 0.75;
  lights.hemi.intensity = 0.5;

  const floor = tilePlane(9, 5);
  floor.rotation.x = -Math.PI / 2;
  floor.position.set(0, 0, 0.6);
  scene.add(floor);

  // A footprint: a sole 25 cm long, toe forward (+z once laid flat), lifted clear of the floor's tiles.
  const sole = new THREE.Shape();
  sole.absellipse(0, 0, 0.05, 0.135, 0, Math.PI * 2, false, 0);
  const soleGeo = new THREE.ShapeGeometry(sole, 16).rotateX(-Math.PI / 2);
  const soleMat = new THREE.MeshStandardMaterial({ roughness: 0.85, metalness: 0, polygonOffset: true, polygonOffsetFactor: -2, polygonOffsetUnits: -2 });
  const DONE = new THREE.Color(0xfbbf24);
  const NOW = new THREE.Color(0xfff1c2);
  const TO_COME = new THREE.Color(0x8f867c);
  const PIVOT = new THREE.Color(0x7dd3fc);
  const PIVOT_NEXT = new THREE.Color(0x4f6f80);
  const pathMat = new THREE.MeshBasicMaterial({ color: 0xf5efe6, transparent: true, opacity: 0.5, depthWrite: false, polygonOffset: true, polygonOffsetFactor: -3, polygonOffsetUnits: -3 });

  let move: PivotMove = "behind";
  let plates = 150;

  const stations: Station[] = STATION_X.map((c, i) => {
    const range = buildFixture(RANGE.id, RANGE, 2);
    scene.add(place(range, c - RANGE.w / 2, 0, true));
    range.hooks.flames?.setAll?.(0);
    range.hooks.flames?.setFlow(1, 0.75);
    const table = buildFixture(TABLE.id, TABLE, 2);
    scene.add(table.group);
    const plate = makePlate();
    scene.add(plate);
    const sk = buildFixture(SKILLET.id, SKILLET, 1);
    const skillet = new THREE.Group();
    // The fixture frame puts the pan's centre at (w/2, 0, -d + w/2), its handle towards +z.
    sk.group.position.set(-SKILLET.w / 2, 0, SKILLET.d - SKILLET.w / 2);
    skillet.add(sk.group);
    skillet.position.set(c + BURNER.x, GRATE_Y + 0.0006, BURNER.z);
    scene.add(skillet);
    const prints = new THREE.InstancedMesh(soleGeo, soleMat, MAX_PRINTS);
    prints.count = 0;
    prints.frustumCulled = false;
    prints.receiveShadow = true;
    scene.add(prints);
    const path = new THREE.Mesh(new THREE.BufferGeometry(), pathMat);
    path.frustumCulled = false;
    path.renderOrder = 1;
    scene.add(path);
    // The same cook at both stations.
    const st = makeStaff("p50", 2, { seed: 41, apron: 0x26344d, sleeves: "cuffed" });
    scene.add(st.group);
    return { c, st, routine: null, steps: [], cost: { steps: 0, pivots: 0, moveS: 0, metres: 0 }, prints, path, table, plate, skillet, range, tools: new Map(), held: null, layout: layoutFor(c, i === 0 ? "built" : move) };
  });

  /** Lays a station out and plans its cook's loop: saute, to the table, plate, back. */
  function build(s: Station, lay: Layout): void {
    s.layout = lay;
    const g = s.table.group;
    g.rotation.set(0, lay.table.yaw, 0);
    if (lay.table.yaw === Math.PI) g.position.set(lay.table.x + TABLE.w, 0, lay.table.z);
    else g.position.set(lay.table.x, 0, lay.table.z);
    s.plate.position.copy(lay.plate).setY(TABLE.h + 0.0008);
    const at = new THREE.Vector3(s.c + BURNER.x, GRATE_Y + 0.03, BURNER.z);
    const steps: RoutineStep[] = [
      { kind: "work", work: "saute", duration: SAUTE_S, opts: { at, surfaceY: GRATE_Y } },
      ...moveSteps(lay.burner, lay.plating),
      { kind: "work", work: "plate", duration: PLATE_S, opts: { at: lay.plate.clone(), surfaceY: TABLE.h } },
      ...moveSteps(lay.plating, lay.burner),
    ];
    const r = new Routine(s.st, lay.burner, steps, true);
    s.routine = r;
    // Each footfall belongs to the routine step it lands in: a turn's are a pivot's, a walk's or a shuffle's are steps.
    const kindAt = (t: number) => {
      let i = 0;
      while (i + 1 < steps.length && r.stepTimes[i + 1]! <= t) i++;
      return steps[i]!.kind;
    };
    s.steps = r.footfalls().map((f) => ({ ...f, pivot: kindAt(f.tOn - 1e-4) === "turn" }));
    const pivots = steps.filter((x, i) => x.kind === "turn" && s.steps.some((f) => f.pivot && f.tOn > r.stepTimes[i]! && f.tOn <= r.stepTimes[i + 1]! + 1e-4)).length;
    // Seconds on the move: everything but the two stretches of work.
    const workS = SAUTE_S + PLATE_S;
    const moveS = Math.max(0, r.duration - workS);
    // Metres walked: each foot's travel, the body's being the mean of the two.
    let metres = 0;
    for (const side of ["L", "R"] as const) {
      const xs = s.steps.filter((f) => f.side === side);
      let prev: { x: number; z: number } | null = null;
      let m = 0;
      for (const f of xs) {
        if (prev) m += Math.hypot(f.x - prev.x, f.z - prev.z);
        prev = f;
      }
      metres += m / 2;
    }
    s.cost = { steps: s.steps.filter((f) => !f.pivot).length, pivots, moveS, metres };
    // The body's path, foot to foot: a line on the floor (a spaghetti diagram), the same at every plate.
    const mids: THREE.Vector2[] = [new THREE.Vector2(lay.burner.x, lay.burner.z)];
    for (let k = 1; k < s.steps.length; k++) {
      const a = s.steps[k - 1]!;
      const b = s.steps[k]!;
      if (a.side !== b.side) mids.push(new THREE.Vector2((a.x + b.x) / 2, (a.z + b.z) / 2));
    }
    mids.push(new THREE.Vector2(lay.burner.x, lay.burner.z));
    s.path.geometry.dispose();
    s.path.geometry = ribbon(mids, 0.035, 0.0025);
    // The footprints, where they fall.
    const n = Math.min(MAX_PRINTS, s.steps.length);
    const m4 = new THREE.Matrix4();
    const q = new THREE.Quaternion();
    const one = new THREE.Vector3(1, 1, 1);
    const pos = new THREE.Vector3();
    for (let k = 0; k < n; k++) {
      const f = s.steps[k]!;
      q.setFromAxisAngle(THREE.Object3D.DEFAULT_UP, f.yaw);
      m4.compose(pos.set(f.x, 0.004, f.z), q, one);
      s.prints.setMatrixAt(k, m4);
      s.prints.setColorAt(k, f.pivot ? PIVOT_NEXT : TO_COME);
    }
    s.prints.count = n;
    s.prints.instanceMatrix.needsUpdate = true;
    if (s.prints.instanceColor) s.prints.instanceColor.needsUpdate = true;
  }
  function buildAll(): void {
    build(stations[0]!, layoutFor(stations[0]!.c, "built"));
    build(stations[1]!, layoutFor(stations[1]!.c, move));
    renderPage();
  }

  /* Camera: from in front and above, both stations across the view. */
  const camera = new THREE.PerspectiveCamera(32, 16 / 10, 0.05, 60);
  const TARGET = new THREE.Vector3(0, 0.6, 0.35);
  const VIEW = new THREE.Vector3(0, 1.05, 1).normalize();
  const keep = [-1, 1].flatMap((sx) => [new THREE.Vector3(sx * 2.4, 1.5, -0.9), new THREE.Vector3(sx * 2.4, 0.1, 2.35), new THREE.Vector3(sx * 0.6, 0.0, 2.35)]);
  const controls = createScrollSafeOrbit(camera, canvas, [TARGET.x, TARGET.y, TARGET.z], { minDistance: 1.2, maxDistance: 14, maxPolarAngle: Math.PI * 0.47 });
  let userMoved = false;
  controls.addEventListener("start", () => {
    userMoved = true;
  });
  const hudEls = Array.from(stageEl.querySelectorAll<HTMLElement>(".pv-hud"));
  const legendEl = stageEl.querySelector<HTMLElement>(".pv-legend");
  const ndc = new THREE.Vector3();
  const frame = () => {
    if (userMoved) return;
    const w = Math.max(1, stageEl.clientWidth);
    const h = Math.max(1, stageEl.clientHeight);
    camera.aspect = Math.max(0.2, w / h);
    camera.clearViewOffset();
    camera.updateProjectionMatrix();
    // The stations fit between the HUDs and the legend, then the view shifts to sit in that band.
    const hudB = hudEls.reduce((m, e) => Math.max(m, e.offsetTop + e.offsetHeight), 0);
    const top = Math.min(0.88, 1 - (2 * (hudB + 6)) / h);
    const bot = -Math.min(0.9, 1 - (2 * ((legendEl ? h - legendEl.offsetTop : 0) + 4)) / h);
    let dist = 2.5;
    let lo = 0;
    let hi = 0;
    for (; dist < 24; dist += 0.06) {
      camera.position.copy(TARGET).addScaledVector(VIEW, dist);
      camera.lookAt(TARGET);
      camera.updateMatrixWorld(true);
      lo = Infinity;
      hi = -Infinity;
      let wide = false;
      for (const p of keep) {
        ndc.copy(p).project(camera);
        lo = Math.min(lo, ndc.y);
        hi = Math.max(hi, ndc.y);
        if (Math.abs(ndc.x) >= 0.96) wide = true;
      }
      if (!wide && hi - lo < top - bot) break;
    }
    const shift = (top + bot) / 2 - (hi + lo) / 2;
    if (Math.abs(shift) > 0.005) {
      camera.setViewOffset(w, h, 0, (shift * h) / 2, w, h);
      camera.updateProjectionMatrix();
    }
    const fog = scene.fog as THREE.Fog;
    fog.near = dist * 1.1;
    fog.far = dist * 2.6;
    controls.target.copy(TARGET);
    controls.update();
  };
  const gate = shaderGate(renderer, camera);
  bindCanvasResize(stageEl, renderer, camera, gate);

  /* The page. */
  const q = <T extends Element = HTMLElement>(sel: string) => root.querySelector<T>(sel);
  const moveBtns = Array.from(root.querySelectorAll<HTMLButtonElement>("[data-pv-move]"));
  const platesIn = q<HTMLInputElement>("[data-pv-plates]");
  const platesVal = q("[data-pv-plates-val]");
  const hudA = q("[data-pv-hud-a]");
  const hudB = q("[data-pv-hud-b]");
  const nameB = q("[data-pv-name-b]");
  const statusEl = q("[data-pv-status]");
  const tile = (k: string) => ({ v: q(`[data-pv-${k}]`), d: q(`[data-pv-${k}-d]`), s: q(`[data-pv-${k}-s]`) });
  const tiles = { steps: tile("steps"), service: tile("service"), time: tile("time"), year: tile("year") };
  const num = (v: number) => Math.round(v).toLocaleString("en-US");
  function setTile(t: { v: HTMLElement | null; d: HTMLElement | null; s: HTMLElement | null }, value: string, sub: string, delta: string | null): void {
    if (t.v) t.v.textContent = value;
    if (t.s) t.s.textContent = sub;
    if (t.d) {
      t.d.hidden = delta === null;
      t.d.textContent = delta ?? "";
      t.d.dataset.good = "true";
    }
  }
  function renderPage(): void {
    const a = stations[0]!.cost;
    const b = stations[1]!.cost;
    if (nameB) nameB.textContent = stations[1]!.layout.name;
    const pv = a.pivots === b.pivots ? `${a.pivots} pivots at both` : `${a.pivots} → ${b.pivots} pivots`;
    setTile(tiles.steps, `${a.steps} → ${b.steps}`, `as built → moved, and ${pv}`, a.steps > b.steps ? `−${a.steps - b.steps}` : null);
    setTile(tiles.service, num((a.steps - b.steps) * plates), `fewer, ${num(a.steps * plates)} → ${num(b.steps * plates)}`, null);
    const minA = (a.moveS * plates) / 60;
    const minB = (b.moveS * plates) / 60;
    setTile(tiles.time, `${num(minA - minB)} min`, `${num(minA)} → ${num(minB)} min on the move`, null);
    const km = ((a.metres - b.metres) * plates * SERVICES_A_YEAR) / 1000;
    setTile(tiles.year, `${num(km)} km`, `less walking, ${SERVICES_A_YEAR} services`, null);
    if (statusEl) {
      const pivots = b.pivots === 1 ? "a pivot" : `${b.pivots} pivots`;
      const how = b.steps === 0 ? `${pivots} and no steps` : `${pivots} and ${b.steps === 1 ? "a step" : `${b.steps} steps`}`;
      statusEl.textContent = `Moved ${move === "beside" ? "beside the range" : "straight behind"}, plating is ${how}: ${a.steps - b.steps} steps a plate fewer than as built, ${num((a.steps - b.steps) * plates)} steps and ${num(minA - minB)} minutes on the move in a service of ${num(plates)} plates, ${num(km)} km a year.`;
      statusEl.dataset.state = "ok";
    }
    moveBtns.forEach((x) => x.setAttribute("aria-pressed", String(x.dataset.pvMove === move)));
    if (platesVal) platesVal.textContent = `${num(plates)} plates`;
  }
  for (const btn of moveBtns)
    btn.addEventListener("click", () => {
      move = btn.dataset.pvMove === "beside" ? "beside" : "behind";
      build(stations[1]!, layoutFor(stations[1]!.c, move));
      t0 = clock;
      renderPage();
    });
  platesIn?.addEventListener("input", () => {
    plates = Number(platesIn.value) || 150;
    renderPage();
  });

  buildAll();
  frame();
  new ResizeObserver(() => frame()).observe(stageEl);
  const mountMs = performance.now() - mountT0;
  for (const s of stations) void s.st.ready.then(() => void gate.refresh(scene));

  /** Pose a station's cook at their own loop time, the tool for the work in hand and the skillet off the grate. */
  function show(s: Station, t: number): void {
    const r = s.routine;
    if (!r) return;
    const lt = t % r.duration;
    const st = r.apply(lt);
    const want = st.activity === "work" && st.work ? (TOOL_FOR[st.work] ?? null) : null;
    if (want !== s.held) {
      if (s.held) s.tools.get(s.held)?.removeFromParent();
      s.held = want;
      if (want) {
        let tool = s.tools.get(want);
        if (!tool) s.tools.set(want, (tool = makeTool(want)));
        s.st.hand("R").add(tool);
        tool.position.set(0, 0, 0);
        tool.quaternion.identity();
      }
    }
    s.skillet.visible = s.held !== "skillet";
    // Footprints: this plate's so far bright, the last one brightest, the rest to come dim.
    let walked = 0;
    for (let k = 0; k < s.prints.count; k++) {
      const f = s.steps[k]!;
      const on = f.tOn <= lt;
      const latest = on && !(k + 1 < s.steps.length && s.steps[k + 1]!.tOn <= lt);
      s.prints.setColorAt(k, f.pivot ? (on ? (latest ? NOW : PIVOT) : PIVOT_NEXT) : on ? (latest ? NOW : DONE) : TO_COME);
      if (on && !f.pivot) walked++;
    }
    if (s.prints.instanceColor) s.prints.instanceColor.needsUpdate = true;
    const el = s === stations[0] ? hudA : hudB;
    const piv = `${s.cost.pivots} pivot${s.cost.pivots === 1 ? "" : "s"} a plate`;
    if (el) el.textContent = s.cost.steps ? `${walked} of ${s.cost.steps} steps · ${piv}` : `No steps · ${piv}`;
  }

  let inView = false;
  observeInView(stageEl, (v) => {
    inView = v;
  });
  let clock = 0;
  let t0 = 0;
  let frozen: number | null = null;
  let last = performance.now();
  const fitClip = clipFitter(camera, scene, 0.05);
  function animate(now: number): void {
    requestAnimationFrame(animate);
    const dt = Math.max(0, Math.min(0.05, (now - last) / 1000));
    last = now;
    if (!inView || document.hidden || !gate.ready()) return;
    clock += dt;
    const t = frozen ?? clock - t0;
    show(stations[0]!, frozen ?? clock);
    show(stations[1]!, t);
    stations.forEach((s) => s.range.hooks.tick?.(clock, dt, 0.5));
    controls.update();
    fitClip();
    renderer.render(scene, camera);
  }
  requestAnimationFrame(animate);
  void gate.refresh(scene);

  // A hook for checking the page by script.
  Object.assign(stageEl, {
    pivot: {
      info: () => ({
        mountMs: Math.round(mountMs),
        tris: renderer.info.render.triangles,
        calls: renderer.info.render.calls,
        move,
        plates,
        built: { ...stations[0]!.cost, duration: +(stations[0]!.routine?.duration ?? 0).toFixed(2) },
        moved: { ...stations[1]!.cost, duration: +(stations[1]!.routine?.duration ?? 0).toFixed(2) },
      }),
      hold(t: number | null) {
        frozen = t;
      },
      view(eye: [number, number, number], at: [number, number, number]) {
        userMoved = true;
        camera.clearViewOffset();
        camera.position.set(...eye);
        controls.target.set(...at);
        controls.update();
      },
    },
  });
}
