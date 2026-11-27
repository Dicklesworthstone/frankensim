/**
 * A combi trips (chapter 07, fifth figure): the morning's kitchen, planned as usual, and the same morning with one of
 * its two combis out of use for a while, the kitchen keeping to its plan (sim/combiTrip.ts, run in a worker).
 *
 * The kitchen on the stage is the morning with the trip: the combi that is out stands dark, its screen off and its
 * program stopped, the floor in front of it red, and a cook with free hands while work waits for a combi wears a pause
 * badge. A camera on a stand at the front left watches the line. Under the stage are two screens from its place:
 * the twin's picture of the planned morning, and the morning with the trip as the camera shows it; they part the
 * minute the combi goes dark, and the twin marks what differs. Under them, a map of every trip on the grid and what
 * it costs the van.
 */
import * as THREE from "three";
import { RoundedBoxGeometry } from "three/examples/jsm/geometries/RoundedBoxGeometry.js";
import { DOOR, MENUS, OVENS, SLOTS, clockOf, madeIt, menuDishes, type Morning } from "../sim/cateringMorning";
import { AT_N, FOR_N, TRIP_AT, TRIP_FOR, costOf, free, tripKey, tripMorning, type Trip, type TripMap, type TripSettings } from "../sim/combiTrip";
import type { TripAsk, TripDone } from "../sim/tripRuns.worker";
import { bindCanvasResize, clipFitter, createScrollSafeOrbit, createStudioRenderer, observeInView, shaderGate } from "./kit";
import { CARRIER, CARRIER_COLS, CARRIER_X0, CUTOFF, END_HOLD_S, LINE_Z, MORNING_PIECES, MorningKitchen, OVEN, OVEN_X0, RANGE, ROLL_MIN, SEC_PER_MIN, SLOT_COLOR, TABLE, TABLE_Z, fmtLate, ovenX, ovenX0, pauseBadgeMaterial, rangeX0, shortName, softPadTexture, type CookSpec } from "./morningKitchen";
import type { Piece } from "./models/all";
import { studioEnvironment, studioLights } from "./models/stage";
import { backdrop, tilePlane } from "./weekLab";

/** Two kitchens' pieces (the twin's and the morning's), for main.ts to build ahead. */
export const TRIP_PIECES: readonly Piece[] = [...MORNING_PIECES, ...MORNING_PIECES];

/** The morning's three cooks, the same in the twin. */
const COOKS: CookSpec[] = [
  { cohort: "p50", seed: 41, apron: 0x26344d },
  { cohort: "p95", seed: 57, apron: 0x3b3d43 },
  { cohort: "p05", seed: 63, apron: 0x45563f },
];
const OVEN_NAME = ["the left combi", "the right combi"] as const;
const num = (v: number) => Math.round(v).toLocaleString("en-US");
const cap = (s: string) => s.replace(/^./, (c) => c.toUpperCase());
/** "1 h 10" or "40 min". */
const span = (m: number) => {
  const r = Math.round(m);
  return r >= 60 ? `${Math.floor(r / 60)} h${r % 60 ? ` ${String(r % 60).padStart(2, "0")}` : ""}` : `${r} min`;
};
/** Plans kept on the page for settings seen before. */
const PLANS_KEPT = 24;
const RED = "#ef4444";
const AMBER = "#fbbf24";

/**
 * The camera's stand on the floor, out past the front left corner of the table (clear of the places idle cooks wait
 * by its end), its head 2.6 m up, and where it looks: across the combis and down the line, the carriers at the far
 * side of its view.
 */
const STAND = new THREE.Vector3(-3.3, 0, 2.55);
const HEAD_Y = 2.6;
const LOOK = new THREE.Vector3(-0.4, 0.45, -0.2);
/** The camera's field of view, degrees up and down. */
const CCTV_FOV = 48;
/** The stand's legs: angles on the floor, one leg straight back, away from the kitchen, and how far they spread. */
const LEG_A = Array.from({ length: 3 }, (_, i) => Math.atan2(STAND.z - LOOK.z, STAND.x - LOOK.x) + (i * 2 * Math.PI) / 3);
const SPREAD = 0.4;
/** Seconds between frames on the camera's screens (each screen draws at half that rate). */
const SCREEN_S = 1 / 16;

/**
 * A camera on a light stand, built for the stage: three legs spread from a collar, a two-piece column with its lock,
 * and a white bullet camera under its sun shield on a tilt bracket, its cable run down the column to the floor. The
 * head looks at `look`; `lens` is where its view starts.
 */
function cameraOnStand(look: THREE.Vector3): { group: THREE.Group; lens: THREE.Vector3; led: THREE.MeshBasicMaterial } {
  const group = new THREE.Group();
  const metal = new THREE.MeshStandardMaterial({ color: 0x1b1c1f, metalness: 0.62, roughness: 0.38 });
  const collar = new THREE.MeshStandardMaterial({ color: 0x2a2b2f, metalness: 0.35, roughness: 0.52 });
  const rubber = new THREE.MeshStandardMaterial({ color: 0x111112, metalness: 0, roughness: 0.92 });
  const shell = new THREE.MeshPhysicalMaterial({ color: 0xe8eaec, metalness: 0, roughness: 0.42, clearcoat: 0.25, clearcoatRoughness: 0.4 });
  const glass = new THREE.MeshPhysicalMaterial({ color: 0x07080a, metalness: 0, roughness: 0.07, clearcoat: 1, clearcoatRoughness: 0.05 });
  const lensGlass = new THREE.MeshPhysicalMaterial({ color: 0x0c1220, metalness: 0.15, roughness: 0.03, clearcoat: 1, clearcoatRoughness: 0.02 });
  const ir = new THREE.MeshStandardMaterial({ color: 0x2b0d0f, metalness: 0, roughness: 0.3 });
  const led = new THREE.MeshBasicMaterial({ color: 0xff2b2b, toneMapped: false });
  const mesh = (g: THREE.BufferGeometry, m: THREE.Material, shadow = true) => {
    const x = new THREE.Mesh(g, m);
    x.castShadow = shadow;
    x.receiveShadow = shadow;
    return x;
  };
  /** A tube of radius r from a to b. */
  const rod = (a: THREE.Vector3, b: THREE.Vector3, r: number, m: THREE.Material, sides = 16) => {
    const d = b.clone().sub(a);
    const x = mesh(new THREE.CylinderGeometry(r, r, d.length(), sides), m);
    x.position.copy(a).addScaledVector(d, 0.5);
    x.quaternion.setFromUnitVectors(new THREE.Vector3(0, 1, 0), d.normalize());
    return x;
  };
  const at = (x: number, y: number, z: number) => new THREE.Vector3(STAND.x + x, y, STAND.z + z);

  /* The stand: legs spread 0.4 m from a collar at 0.66 m, braced to the column; rubber feet. */
  const HUB = 0.66;
  for (const a of LEG_A) {
    const foot = at(Math.cos(a) * SPREAD, 0.03, Math.sin(a) * SPREAD);
    group.add(rod(at(Math.cos(a) * 0.03, HUB, Math.sin(a) * 0.03), foot, 0.0105, metal));
    const f = mesh(new THREE.CylinderGeometry(0.017, 0.019, 0.03, 18), rubber);
    f.position.set(foot.x, 0.015, foot.z);
    group.add(f);
    // The brace: from a collar low on the column to the leg's middle.
    const mid = at(Math.cos(a) * SPREAD * 0.55, HUB - (HUB - 0.03) * 0.55, Math.sin(a) * SPREAD * 0.55);
    group.add(rod(at(Math.cos(a) * 0.02, 0.34, Math.sin(a) * 0.02), mid, 0.0055, metal, 10));
  }
  const ring = (y: number, r: number, h: number) => {
    const x = mesh(new THREE.CylinderGeometry(r, r, h, 24), collar);
    x.position.set(STAND.x, y, STAND.z);
    group.add(x);
  };
  ring(HUB, 0.027, 0.075);
  ring(0.34, 0.022, 0.035);
  group.add(rod(at(0, 0.3, 0), at(0, 1.33, 0), 0.0175, metal, 20));
  // The height lock: a collar, and its T-knob out to one side.
  ring(1.33, 0.024, 0.05);
  group.add(rod(at(0.024, 1.33, 0), at(0.052, 1.33, 0), 0.0045, collar, 10));
  const knob = mesh(new THREE.CylinderGeometry(0.004, 0.004, 0.04, 10), rubber);
  knob.position.copy(at(0.056, 1.33, 0));
  knob.rotation.x = Math.PI / 2;
  group.add(knob);
  group.add(rod(at(0, 1.33, 0), at(0, HEAD_Y - 0.06, 0), 0.0128, metal, 20));
  // The spigot and the tilt bracket's foot.
  group.add(rod(at(0, HEAD_Y - 0.07, 0), at(0, HEAD_Y - 0.02, 0), 0.0085, collar, 14));

  /* The head: a bullet camera, its lens along its +z, turned to look at `look`. */
  const head = new THREE.Group();
  head.position.copy(at(0, HEAD_Y, 0));
  head.lookAt(look);
  group.add(head);
  // The bracket: a block on the spigot, the camera on a short arm above it.
  const block = mesh(new RoundedBoxGeometry(0.034, 0.03, 0.04, 2, 0.004), collar);
  block.position.set(0, -0.005, 0);
  head.add(block);
  const arm = mesh(new THREE.CylinderGeometry(0.007, 0.007, 0.03, 12), collar);
  arm.position.set(0, 0.022, -0.02);
  head.add(arm);
  const BODY = { w: 0.08, h: 0.074, d: 0.2, y: 0.074 };
  const body = mesh(new RoundedBoxGeometry(BODY.w, BODY.h, BODY.d, 3, 0.012), shell);
  body.position.set(0, BODY.y, 0);
  head.add(body);
  // The sun shield: a curved hood over the top, clear of the body all along, a little past its front.
  const SH = { r: 0.062, arc: 1.6, len: 0.25 };
  const shield = mesh(new THREE.CylinderGeometry(SH.r, SH.r, SH.len, 28, 1, true, Math.PI - SH.arc / 2, SH.arc).rotateX(Math.PI / 2), new THREE.MeshPhysicalMaterial({ color: 0xe8eaec, roughness: 0.42, clearcoat: 0.25, clearcoatRoughness: 0.4, side: THREE.DoubleSide }));
  shield.position.set(0, BODY.y, 0.014);
  head.add(shield);
  // The front: black glass in the face's flat middle, the lens proud of it, a ring of infrared lamps, the red lamp.
  const front = BODY.d / 2;
  const bezel = mesh(new THREE.CircleGeometry(0.024, 40), glass, false);
  bezel.position.set(0, BODY.y, front + 0.0012);
  head.add(bezel);
  const barrel = mesh(new THREE.CylinderGeometry(0.0125, 0.0135, 0.012, 32).rotateX(Math.PI / 2), metal, false);
  barrel.position.set(0, BODY.y, front + 0.0012 + 0.006);
  head.add(barrel);
  const lensDisc = mesh(new THREE.CircleGeometry(0.0108, 32), lensGlass, false);
  lensDisc.position.set(0, BODY.y, front + 0.0012 + 0.0122);
  head.add(lensDisc);
  for (let i = 0; i < 6; i++) {
    const a = (i / 6) * Math.PI * 2 + Math.PI / 6;
    const d = mesh(new THREE.CircleGeometry(0.0021, 14), ir, false);
    d.position.set(Math.cos(a) * 0.0185, BODY.y + Math.sin(a) * 0.0185, front + 0.0012 + 0.0004);
    head.add(d);
  }
  const lamp = mesh(new THREE.CircleGeometry(0.0019, 14), led, false);
  lamp.position.set(0.026, BODY.y + 0.024, front + 0.0004);
  head.add(lamp);
  // The cable gland at the back.
  const gland = mesh(new THREE.CylinderGeometry(0.0065, 0.0075, 0.016, 14).rotateX(Math.PI / 2), rubber);
  gland.position.set(0, BODY.y - 0.012, -front - 0.008);
  head.add(gland);
  head.updateMatrixWorld(true);
  const lens = head.localToWorld(new THREE.Vector3(0, BODY.y, front + 0.02));

  /* The cable: out of the gland, down behind the head, along the column a few centimetres off it, to the floor. */
  const out = head.localToWorld(new THREE.Vector3(0, BODY.y - 0.012, -front - 0.016));
  const back = head.localToWorld(new THREE.Vector3(0, BODY.y - 0.06, -front - 0.07));
  // Down the column's back, and off along the floor away from the kitchen, beside the back leg.
  const bx = Math.cos(LEG_A[0]!);
  const bz = Math.sin(LEG_A[0]!);
  const sx = -bz * 0.012;
  const sz = bx * 0.012;
  const off = (k: number, y: number, s = 1) => at(bx * 0.03 * k + sx * s, y, bz * 0.03 * k + sz * s);
  const floor = (d: number, s: number) => at(bx * d - bz * s, 0.0036, bz * d + bx * s);
  const pts = [out, back, off(1, HEAD_Y - 0.22), off(1, 1.6), off(1.2, 1.0), off(1.5, 0.72, 1.5), off(3.6, 0.2, 2), floor(0.3, 0.1), floor(0.75, 0.24)];
  const cable = mesh(new THREE.TubeGeometry(new THREE.CatmullRomCurve3(pts, false, "centripetal"), 140, 0.0032, 8, false), rubber);
  group.add(cable);
  return { group, lens, led };
}

export function mountTripLab(root: HTMLElement): void {
  const mountT0 = performance.now();
  const stageEl = root.querySelector<HTMLElement>("[data-trip-stage]");
  const canvas = root.querySelector<HTMLCanvasElement>("[data-trip-canvas]");
  const camCanvas = { twin: root.querySelector<HTMLCanvasElement>('[data-trip-cam="twin"]'), seen: root.querySelector<HTMLCanvasElement>('[data-trip-cam="seen"]') };
  const monsEl = root.querySelector<HTMLElement>(".tp-mons");
  if (!stageEl || !canvas || !camCanvas.twin || !camCanvas.seen || !monsEl) return;
  const maybeRenderer = createStudioRenderer(canvas, false);
  if (!maybeRenderer) return;
  const renderer = maybeRenderer;
  renderer.toneMappingExposure = 1.0;
  const rTwin = createStudioRenderer(camCanvas.twin, false);
  const rSeen = createStudioRenderer(camCanvas.seen, false);
  if (!rTwin || !rSeen) return;
  // A camera's picture, not a photograph: one pixel to a pixel of the page.
  for (const r of [rTwin, rSeen]) {
    r.setPixelRatio(1);
    r.toneMappingExposure = 1.0;
  }
  const phone = window.matchMedia("(max-width: 720px)").matches;

  /* Two scenes lit alike: the morning with the trip, on the stage and the camera's right-hand screen, and the twin's
   * planned morning, on its left-hand screen. */
  const env = studioEnvironment(renderer);
  const makeScene = (kitchen: MorningKitchen, shadowMap: number) => {
    const scene = new THREE.Scene();
    scene.background = backdrop();
    scene.fog = new THREE.Fog(0x151510, 8, 16);
    scene.environment = env;
    scene.environmentIntensity = 0.6;
    const lights = studioLights(scene, new THREE.Vector3(-0.4, 0.9, 0.1), 5.2, shadowMap);
    lights.key.position.set(-1.4, 6.4, 4.6);
    lights.key.intensity = 1.9;
    lights.key.shadow.normalBias = 0.02;
    lights.rim.position.set(3.4, 3.2, -3.0);
    lights.rim.intensity = 0.75;
    lights.hemi.intensity = 0.5;
    const floor = tilePlane(16, 10);
    floor.rotation.x = -Math.PI / 2;
    floor.position.set(0.2, 0, 0.6);
    scene.add(floor);
    kitchen.rollOut = true;
    scene.add(kitchen.group);
    return scene;
  };
  const kitchen = new MorningKitchen(COOKS);
  const twin = new MorningKitchen(COOKS);
  const scene = makeScene(kitchen, phone ? 1024 : 2048);
  const twinScene = makeScene(twin, phone ? 1024 : 2048);

  /* What only the stage shows: the camera on its stand and its view's edges, the red floor, the pause badges. */
  const marks = new THREE.Group();
  scene.add(marks);
  const rig = cameraOnStand(LOOK);
  marks.add(rig.group);
  const outPads = Array.from({ length: OVENS }, (_, o) => {
    const mat = new THREE.MeshBasicMaterial({ color: RED, map: softPadTexture(), toneMapped: false, transparent: true, opacity: 0, depthWrite: false, polygonOffset: true, polygonOffsetFactor: -2, polygonOffsetUnits: -2 });
    const mesh = new THREE.Mesh(new THREE.PlaneGeometry(OVEN.w * 0.96, 0.62).rotateX(-Math.PI / 2), mat);
    mesh.position.set(ovenX(o), 0.004, LINE_Z + 0.36);
    mesh.renderOrder = 1;
    mesh.visible = false;
    marks.add(mesh);
    return { mesh, mat, on: 0 };
  });
  const badgeMat = pauseBadgeMaterial(AMBER);
  const badges = kitchen.cooks.map(() => {
    const b = new THREE.Sprite(badgeMat);
    b.scale.setScalar(0.26);
    b.renderOrder = 5;
    b.visible = false;
    marks.add(b);
    return b;
  });

  /* The camera: from its lens, a wide view across the table to the line. */
  const cctv = new THREE.PerspectiveCamera(CCTV_FOV, 16 / 9, 0.05, 30);
  cctv.position.copy(rig.lens);
  cctv.lookAt(LOOK);
  cctv.updateMatrixWorld(true);
  // Its view's edges, from the lens, fading out over 2.4 m.
  {
    const pos: number[] = [];
    const col: number[] = [];
    const c = cctv.position;
    cctv.updateProjectionMatrix();
    for (const [sx, sy] of [
      [-1, -1],
      [1, -1],
      [1, 1],
      [-1, 1],
    ] as const) {
      const p = new THREE.Vector3(sx, sy, 0.5).unproject(cctv).sub(c).normalize().multiplyScalar(2.4).add(c);
      pos.push(c.x, c.y, c.z, p.x, p.y, p.z);
      col.push(1, 1, 1, 0.5, 1, 1, 1, 0);
    }
    const g = new THREE.BufferGeometry();
    g.setAttribute("position", new THREE.Float32BufferAttribute(pos, 3));
    g.setAttribute("color", new THREE.Float32BufferAttribute(col, 4));
    marks.add(new THREE.LineSegments(g, new THREE.LineBasicMaterial({ vertexColors: true, transparent: true, depthWrite: false, toneMapped: false })));
  }

  /* The stage's camera: from behind the worktable, raised, the whole line and the camera's stand in view. */
  const camera = new THREE.PerspectiveCamera(32, 16 / 9, 0.05, 60);
  const TARGET = new THREE.Vector3();
  const VIEW = new THREE.Vector3();
  const keep: THREE.Vector3[] = [];
  const LEFT = OVEN_X0[0];
  const RIGHT = Math.max(rangeX0(0) + RANGE.w, CARRIER_X0 + CARRIER_COLS * (CARRIER.w + 0.03));
  const aim = () => {
    keep.length = 0;
    const wide = stageEl.clientWidth >= 1.3 * stageEl.clientHeight;
    TARGET.set((STAND.x + RIGHT) / 2, 0.75, 0.5);
    VIEW.set(0.06, wide ? 1.15 : 1.35, 1).normalize();
    for (const x of [LEFT, RIGHT]) keep.push(new THREE.Vector3(x, OVEN.h + 0.04, LINE_Z - OVEN.d), new THREE.Vector3(x, 0.2, TABLE_Z + TABLE.d));
    keep.push(new THREE.Vector3(RIGHT, 3 * CARRIER.h, TABLE_Z));
    // The camera's head, and its stand's feet.
    keep.push(new THREE.Vector3(STAND.x, HEAD_Y + 0.12, STAND.z));
    for (const a of LEG_A) keep.push(new THREE.Vector3(STAND.x + Math.cos(a) * SPREAD, 0, STAND.z + Math.sin(a) * SPREAD));
  };
  const controls = createScrollSafeOrbit(camera, canvas, [0, 0.75, 0], { minDistance: 1.2, maxDistance: 16, maxPolarAngle: Math.PI * 0.47 });
  let userMoved = false;
  controls.addEventListener("start", () => {
    userMoved = true;
  });
  const ndc = new THREE.Vector3();
  const hudEl = stageEl.querySelector<HTMLElement>(".tp-hud");
  const legendEl = stageEl.querySelector<HTMLElement>(".tp-legend");
  /** The HUD's lowest edge so far at this size: a HUD that wraps for a moment and back must not refit the view. */
  let hudMax = 0;
  let fitW = 0;
  let fitH = 0;
  const frame = () => {
    if (userMoved) return;
    aim();
    const w = Math.max(1, stageEl.clientWidth);
    const h = Math.max(1, stageEl.clientHeight);
    if (w !== fitW || h !== fitH) hudMax = 0;
    fitW = w;
    fitH = h;
    hudMax = Math.max(hudMax, hudEl ? hudEl.offsetTop + hudEl.offsetHeight : 0);
    camera.aspect = Math.max(0.2, w / h);
    camera.clearViewOffset();
    camera.updateProjectionMatrix();
    // The kitchen fits between the HUD and the legend: fitted to that band's height and the stage's width, then
    // shifted into the band and across to the middle.
    const top = Math.min(0.86, 1 - (2 * (hudMax + 6)) / h);
    const bot = -Math.min(0.92, 1 - (2 * ((legendEl ? h - legendEl.offsetTop : 0) + 4)) / h);
    let dist = 2.5;
    let lo = 0;
    let hi = 0;
    let left = 0;
    let right = 0;
    for (; dist < 22; dist += 0.06) {
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
    const fog = scene.fog as THREE.Fog;
    fog.near = dist * 1.1;
    fog.far = dist * 2.6;
    controls.target.copy(TARGET);
    controls.update();
  };
  const gate = shaderGate(renderer, camera);
  const twinGate = shaderGate(rTwin, cctv);
  const seenGate = shaderGate(rSeen, cctv);
  bindCanvasResize(stageEl, renderer, camera, gate);
  const screenEl = { twin: camCanvas.twin.parentElement!, seen: camCanvas.seen.parentElement! };
  bindCanvasResize(screenEl.twin, rTwin, cctv, twinGate);
  bindCanvasResize(screenEl.seen, rSeen, cctv, seenGate);
  frame();
  new ResizeObserver(() => {
    const w = stageEl.clientWidth;
    const h = stageEl.clientHeight;
    if (w === fitW && h === fitH && (hudEl ? hudEl.offsetTop + hudEl.offsetHeight : 0) <= hudMax) return;
    frame();
  }).observe(stageEl);
  if (hudEl)
    new ResizeObserver(() => {
      if (hudEl.offsetTop + hudEl.offsetHeight > hudMax) frame();
    }).observe(hudEl);

  /* The twin's mark on the camera's picture: a box round the combi that is out, where the camera shows it. */
  const boxEl = root.querySelector<HTMLElement>("[data-tp-box]");
  const placeBox = (o: number) => {
    if (!boxEl) return;
    cctv.updateMatrixWorld(true);
    let x0 = Infinity;
    let x1 = -Infinity;
    let y0 = Infinity;
    let y1 = -Infinity;
    for (let i = 0; i < 8; i++) {
      ndc.set(i & 1 ? ovenX0(o) + OVEN.w : ovenX0(o), i & 2 ? OVEN.h : 0, i & 4 ? LINE_Z : LINE_Z - OVEN.d).project(cctv);
      x0 = Math.min(x0, ndc.x);
      x1 = Math.max(x1, ndc.x);
      y0 = Math.min(y0, ndc.y);
      y1 = Math.max(y1, ndc.y);
    }
    const pct = (v: number) => `${(Math.max(0, Math.min(1, v)) * 100).toFixed(2)}%`;
    boxEl.style.left = pct((x0 + 1) / 2);
    boxEl.style.width = pct((x1 - x0) / 2);
    boxEl.style.top = pct((1 - y1) / 2);
    boxEl.style.height = pct((y1 - y0) / 2);
  };

  /* ---------------------------------------------------------------------------------------------------------- */
  /* The page.                                                                                                   */
  /* ---------------------------------------------------------------------------------------------------------- */

  const q = <T extends Element = HTMLElement>(s: string) => root.querySelector<T>(s);
  const ovenBtns = Array.from(root.querySelectorAll<HTMLButtonElement>("[data-tp-oven]"));
  const atIn = q<HTMLInputElement>("[data-tp-at]");
  const atVal = q("[data-tp-at-val]");
  const forIn = q<HTMLInputElement>("[data-tp-for]");
  const forVal = q("[data-tp-for-val]");
  const menuEl = q("[data-tp-menu]");
  const mealsIn = q<HTMLInputElement>("[data-tp-meals]");
  const mealsVal = q("[data-tp-meals-val]");
  const statusEl = q("[data-tp-status]");
  const mapEl = q("[data-trip-map]");
  const mapNote = q("[data-trip-map-note]");
  const hudClock = q("[data-tp-hud-clock]");
  const hudVan = q("[data-tp-hud-van]");
  const hudNote = q("[data-tp-hud-note]");
  const flagEl = q("[data-tp-flag]");
  const tsEl = { twin: q('[data-tp-ts="twin"]'), seen: q('[data-tp-ts="seen"]') };
  const seenMon = q('[data-tp-mon="seen"]');
  const tile = (k: string) => ({ v: q(`[data-tp-${k}]`), d: q(`[data-tp-${k}-d]`), s: q(`[data-tp-${k}-s]`) });
  const tiles = { van: tile("van"), known: tile("known"), free: tile("free"), best: tile("best") };
  const tilesEl = q(".tp-lab .wk-tiles");
  q(".tp-keys .load")?.style.setProperty("--k", `linear-gradient(90deg, ${SLOT_COLOR.map((c, i) => `${c} ${(i * 100) / SLOT_COLOR.length}% ${((i + 1) * 100) / SLOT_COLOR.length}%`).join(", ")})`);
  const setTile = (t: (typeof tiles)[keyof typeof tiles], v: string, s: string, d: string | null = null, good = true) => {
    if (t.v) t.v.textContent = v;
    if (t.s) t.s.textContent = s;
    if (t.d) {
      t.d.hidden = d === null;
      t.d.textContent = d ?? "";
      if (d !== null) t.d.dataset.good = String(good);
    }
  };

  // Every dish swapped: the menu that leans on the combis.
  const settings: TripSettings = { menu: MENUS - 1, meals: Number(mealsIn?.value ?? 800) };
  const trip: Trip = { oven: 0, at: Number(atIn?.value ?? 200), mins: Number(forIn?.value ?? 60) };
  const runs = new Map<string, TripDone>();
  /** The plan and map on screen, and the morning with the trip (they stay while the next are made). */
  let shown: TripDone | null = null;
  let shownKey = "";
  let actual: Morning | null = null;
  /** The trip of the morning on the stage. */
  let shownTrip: Trip = { ...trip };
  const ready = () => shown !== null && shownKey === tripKey(settings);

  /* The menu card: each place's dish, where it cooks and how long it holds; a press swaps it. */
  const slotBtns: HTMLButtonElement[] = [];
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
  }
  const holdBars = (min: number) => {
    const n = min >= 210 ? 4 : min >= 140 ? 3 : min >= 85 ? 2 : 1;
    return `<span class="wk-keeps" aria-hidden="true">${[1, 2, 3, 4].map((k) => `<i${k <= n ? ' class="on"' : ""}></i>`).join("")}<small>${min >= 120 ? `${min / 60} h` : `${min} min`}</small></span>`;
  };
  /** Where a dish cooks: whether a combi going out can hold it up. */
  const whereOf = (k: number) => {
    const st = menuDishes(settings.menu)[k]!.stages;
    const range = st.some((x) => x.at === "range");
    const oven = st.some((x) => x.at === "oven");
    if (range && oven) return ["Range, then combi", "both"];
    if (range) return ["Range only", "range"];
    return ["Combi", "oven"];
  };
  function renderControls(): void {
    const dishes = menuDishes(settings.menu);
    slotBtns.forEach((b, k) => {
      const d = dishes[k]!;
      const swapped = (settings.menu & (1 << k)) !== 0;
      const other = swapped ? SLOTS[k]!.written : SLOTS[k]!.swap;
      const [where, at] = whereOf(k);
      b.innerHTML = `<b>${d.name}</b>${holdBars(d.holdMin)}<em data-at="${at}">${where}</em>`;
      b.setAttribute("aria-pressed", String(swapped));
      b.title = `Swap for ${other.name.toLowerCase()}`;
    });
    for (const b of ovenBtns) b.setAttribute("aria-pressed", String(Number(b.dataset.tpOven) === trip.oven));
    if (atIn) atIn.value = String(trip.at);
    if (forIn) forIn.value = String(trip.mins);
    if (atVal) atVal.textContent = clockOf(trip.at);
    if (forVal) forVal.textContent = `${trip.mins} min, back at ${clockOf(trip.at + trip.mins)}`;
    if (mealsVal) mealsVal.textContent = `${num(settings.meals)} meals`;
  }

  /* The map: every trip on the grid, for each combi, and what it costs; above each, its loads as planned. */
  const costAt = (map: TripMap, o: number, k: number) => ({ late: map.late[o]![k]!, pastBest: map.pastBest[o]![k]! });
  const freeShare = (map: TripMap) => {
    let n = 0;
    for (let o = 0; o < OVENS; o++) for (let k = 0; k < AT_N * FOR_N; k++) if (free(costAt(map, o, k))) n++;
    return { n, of: OVENS * AT_N * FOR_N };
  };
  const cellClass = (c: { late: number; pastBest: number }) => (c.late > 0.5 ? (c.late >= 59.5 ? "l4" : c.late >= 29.5 ? "l3" : c.late >= 9.5 ? "l2" : "l1") : c.pastBest > 0 ? "pb" : "ok");
  const costText = (c: { late: number; pastBest: number }) =>
    !Number.isFinite(c.late) ? `the van still waiting at ${CUTOFF}` : c.late > 0.5 ? `the van ${span(c.late)} later than planned${c.pastBest ? `, ${num(c.pastBest)} more portions past their best` : ""}` : c.pastBest ? `the van as planned, ${num(c.pastBest)} more portions past their best` : "costs nothing";
  function renderMap(): void {
    if (!mapEl) return;
    const stale = !ready();
    const W = Math.max(260, Math.round(mapEl.clientWidth));
    const narrow = W < 520;
    const L = narrow ? 52 : 112;
    const R = W - 6;
    const cw = (R - L) / AT_N;
    const ch = narrow ? 11 : 13;
    const lane = narrow ? 9 : 12;
    const gap = narrow ? 30 : 34;
    const block = lane + 6 + FOR_N * ch;
    // On a narrow map each combi's name stands over its grid.
    const T = narrow ? 16 : 6;
    const H = T + OVENS * block + (OVENS - 1) * gap + 22;
    const xOf = (m: number) => L + ((m - TRIP_AT.from) / TRIP_AT.step) * cw;
    const plan = shown?.plan ?? null;
    const map = shown?.map ?? null;
    let svg = `<svg viewBox="0 0 ${W} ${H}" width="${W}" height="${H}" class="rs-svg tp-svg${stale ? " stale" : ""}">`;
    for (let o = 0; o < OVENS; o++) {
      const y0 = T + o * (block + gap);
      const gy = y0 + lane + 6;
      // The combi's loads in the plan: when it holds food, by dish.
      svg += `<text class="nm${narrow ? " st" : ""}" x="${narrow ? L : L - 8}" y="${narrow ? y0 - 1 : y0 + lane - 1}">${cap(OVEN_NAME[o]!)}</text>`;
      if (!narrow) svg += `<text class="of" x="4" y="${gy + 4 * ch + ch / 2 + 3.5}">out for</text>`;
      if (plan)
        for (const l of plan.loads) {
          if (l.oven !== o) continue;
          const a = Math.max(TRIP_AT.from, l.t0);
          const b = Math.min(TRIP_AT.to + TRIP_AT.step, l.t3);
          if (b <= a) continue;
          svg += `<rect class="load" x="${xOf(a).toFixed(1)}" y="${(narrow ? y0 + 2 : y0).toFixed(1)}" width="${Math.max(1, xOf(b) - xOf(a)).toFixed(1)}" height="${narrow ? lane - 2 : lane}" style="fill:${SLOT_COLOR[l.slot]}"><title>${l.verb}, ${clockOf(l.t0)} to ${clockOf(l.t3)}, as planned</title></rect>`;
        }
      for (let r = 0; r < FOR_N; r++) {
        const dur = TRIP_FOR.from + r * TRIP_FOR.step;
        if (r % (narrow ? 4 : 2) === 0) svg += `<text class="ya" x="${L - 5}" y="${gy + r * ch + ch / 2 + 3.5}">${dur} min</text>`;
        for (let c = 0; c < AT_N; c++) {
          const at = TRIP_AT.from + c * TRIP_AT.step;
          const k = r * AT_N + c;
          const cost = map && plan?.finished ? costAt(map, o, k) : null;
          const on = trip.oven === o && trip.at === at && trip.mins === dur;
          svg += `<rect class="cell ${cost ? cellClass(cost) : plan && !plan.finished ? "na" : "pend"}${on ? " on" : ""}" data-o="${o}" data-at="${at}" data-mins="${dur}" x="${(xOf(at) + 0.5).toFixed(1)}" y="${(gy + r * ch + 0.5).toFixed(1)}" width="${(cw - 1).toFixed(1)}" height="${ch - 1}"><title>${cap(OVEN_NAME[o]!)} out ${clockOf(at)} to ${clockOf(at + dur)}: ${cost ? costText(cost) : plan && !plan.finished ? `no van to set it against: the plan is not packed by ${CUTOFF}` : "…"}</title></rect>`;
        }
      }
    }
    // The clock under the bottom grid, and the van at its right edge.
    const axisY = T + OVENS * block + (OVENS - 1) * gap + 14;
    // An hour's label stands clear of the van's, which ends at the grid's right edge.
    const vanW = narrow ? 34 : 60;
    for (let m = TRIP_AT.from; m < DOOR; m += 60) if (xOf(m) + 17 + 6 < xOf(DOOR) - vanW) svg += `<text class="xa" x="${xOf(m).toFixed(1)}" y="${axisY}">${clockOf(m)}</text>`;
    svg += `<text class="xa vanl" x="${xOf(DOOR).toFixed(1)}" y="${axisY}" text-anchor="end">${narrow ? "" : "van "}10:30</text>`;
    if (stale) svg += `<text class="planning" x="${((L + R) / 2).toFixed(1)}" y="${(H / 2).toFixed(1)}">Mapping every trip…</text>`;
    svg += `</svg>`;
    mapEl.innerHTML = svg;
    if (mapNote) {
      if (!map || !plan) mapNote.textContent = "Mapping every trip…";
      else if (!plan.finished) mapNote.textContent = `Across, when the combi trips; down, how long it is out. This order is not packed by ${CUTOFF} even with both combis, so there is no van to set a trip against.`;
      else {
        const { n, of } = freeShare(map);
        const braise = plan.loads.filter((l) => l.verb === "Braise");
        // The shortest trip that costs anything.
        let shortest = Infinity;
        for (let o = 0; o < OVENS; o++) for (let k = 0; k < AT_N * FOR_N; k++) if (!free(costAt(map, o, k))) shortest = Math.min(shortest, TRIP_FOR.from + Math.floor(k / AT_N) * TRIP_FOR.step);
        let text = `Across, when the combi trips; down, how long it is out. ${num(n)} of the ${num(of)} trips cost nothing at ${num(settings.meals)} meals.`;
        if (n === of) text += plan.loads.length === 0 ? " No dish on this menu goes in a combi." : " The combis have time to spare all morning.";
        else if (n / of < 0.3 && new Set(braise.map((l) => l.oven)).size > 1) text += ` The braise takes both combis${braise.every((l) => l.t0 < 30) ? " from the start" : ""}, and leaves little room to lose.`;
        else if (n / of > 0.8 && shortest > TRIP_FOR.from) text += ` The combis have time to spare: no trip shorter than ${shortest} min costs anything.`;
        mapNote.textContent = `${text} Pick a cell to trip it.`;
      }
    }
  }
  mapEl?.addEventListener("click", (e) => {
    const r = (e.target as Element).closest<SVGRectElement>("rect.cell");
    if (!r) return;
    trip.oven = Number(r.dataset.o);
    trip.at = Number(r.dataset.at);
    trip.mins = Number(r.dataset.mins);
    changed(true);
  });

  /* Tiles and status: the morning with the trip, against the plan. */
  /** What the combi held when it tripped, as planned. A load the plan starts that very minute is the one due in: with
   * the combi out it never goes in. */
  const heldAt = (plan: Morning, o: number, t: number) => {
    const l = plan.loads.find((x) => x.oven === o && x.t0 < t && x.t3 > t);
    const dishes = menuDishes(plan.settings.menu);
    if (!l) {
      const next = plan.loads.filter((x) => x.oven === o && x.t0 >= t).sort((a, b) => a.t0 - b.t0)[0];
      return { load: null, text: next ? `while it stands empty, its next load (the ${shortName(dishes[next.slot]!.name).toLowerCase()}) due in at ${clockOf(next.t0)}` : "while it stands empty for the rest of the morning" };
    }
    const name = shortName(dishes[l.slot]!.name).toLowerCase();
    if (t < l.t1) return { load: l, text: `while a cook loads it with the ${name}` };
    if (t < l.t2) return { load: l, text: t - l.t1 < 0.5 ? `with the ${name} in it, its cooking just started` : `with the ${name} in it, ${span(t - l.t1)} into ${span(l.t2 - l.t1)} of cooking` };
    if (t < l.tu) return { load: l, text: `with the ${name} done in it, waiting for a cook to take it out` };
    return { load: l, text: `while a cook takes the ${name} out` };
  };
  function renderTiles(): void {
    tilesEl?.toggleAttribute("data-stale", !ready());
    if (!shown || !actual) {
      for (const t of [tiles.van, tiles.known, tiles.best]) setTile(t, "…", "");
    } else {
      const { plan } = shown;
      const c = costOf(actual, plan);
      // Between the clocks as shown, so the difference reads off them.
      const d = Math.round(actual.depart) - Math.round(plan.depart);
      if (!actual.finished) setTile(tiles.van, "Not yet", `still waiting at ${CUTOFF}, some dishes not packed`);
      else setTile(tiles.van, clockOf(actual.depart), actual.late > 0.5 ? fmtLate(actual.late).toLowerCase() : "on time", d > 0 && plan.finished ? `+${d} min` : null, false);
      const at = shownTrip.at;
      if (free(c)) setTile(tiles.known, clockOf(at), "and the van still leaves as planned");
      else setTile(tiles.known, clockOf(at), at < DOOR ? `${span(DOOR - at)} before the van is due` : "as the van is due");
      if (!plan.finished || !actual.finished) setTile(tiles.best, "–", `the van not away by ${CUTOFF}`);
      else setTile(tiles.best, c.pastBest ? `+${num(c.pastBest)}` : "+0", c.pastBest ? "portions past their best at lunch" : "no portion more past its best", null);
    }
    if (!shown) setTile(tiles.free, "…", "");
    else if (!shown.plan.finished) setTile(tiles.free, "–", `not packed by ${CUTOFF} even with both combis`);
    else {
      const { n, of } = freeShare(shown.map);
      setTile(tiles.free, `${Math.round((100 * n) / of)}%`, `${num(n)} of ${num(of)} trips mapped, both combis`);
    }
  }

  function renderStatus(): void {
    if (!statusEl) return;
    if (failed) {
      statusEl.textContent = failed;
      statusEl.dataset.state = "warn";
      return;
    }
    if (!ready() || !actual) {
      statusEl.textContent = "Planning the morning and mapping every trip…";
      statusEl.dataset.state = "warn";
      return;
    }
    const { plan } = shown!;
    const c = costOf(actual, plan);
    const tr = shownTrip;
    const held = heldAt(plan, tr.oven, tr.at);
    let text = `${cap(OVEN_NAME[tr.oven]!)} trips at ${clockOf(tr.at)} ${held.text}, and is back at ${clockOf(tr.at + tr.mins)}.`;
    if (!plan.finished) text = `This order does not make the van even without a trip: the plan is still waiting at ${CUTOFF}.`;
    else if (!madeIt(plan)) text += ` The plan itself ${plan.late > 0.5 ? `leaves at ${clockOf(plan.depart)}, ${fmtLate(plan.late).toLowerCase()}` : `has ${num(plan.pastBest)} portions past their best`} for this order; against it,`;
    else text += " Against the plan,";
    if (plan.finished) {
      if (free(c)) text += ` it costs nothing: the van leaves at ${clockOf(actual.depart)}${madeIt(actual) ? " with every portion at its best" : ""}.`;
      else if (!actual.finished) text += ` the van is still waiting at ${CUTOFF}.`;
      else {
        const late = c.late > 0.5;
        // The load the trip held up, and when it came out instead.
        const after = held.load ? actual.loads.find((x) => x.oven === held.load!.oven && x.t0 === held.load!.t0) : undefined;
        const dish = after ? shortName(menuDishes(plan.settings.menu)[after.slot]!.name).toLowerCase() : "";
        const out = after && held.load && after.t2 > held.load.t2 + 0.5 ? ` The ${dish} come${/s$/.test(dish) ? "" : "s"} out at ${clockOf(after.t2)}, not ${clockOf(held.load.t2)}.` : "";
        // Against a plan that is itself late, the trip's minutes are minutes later, not minutes late.
        const planLate = plan.late > 0.5;
        text += `${late ? ` the van leaves at ${clockOf(actual.depart)}, ${planLate ? `${span(c.late)} later` : fmtLate(c.late).toLowerCase()}` : planLate ? " the van leaves as planned" : " the van leaves on time"}${c.pastBest ? `${late ? " and" : " but"} ${num(c.pastBest)} more portions are past their best at lunch` : ""}.${out}`;
        if (tr.at < DOOR) text += ` The camera shows it dark at ${clockOf(tr.at)}, ${span(DOOR - tr.at)} before the van is due, and the map's ${clockOf(tr.at)} column shows what it costs the van however long it stays out.`;
      }
    }
    statusEl.textContent = text;
    statusEl.dataset.state = !plan.finished || !actual.finished ? "danger" : free(c) ? "ok" : c.late > 30 ? "danger" : "warn";
  }

  /* ---------------------------------------------------------------------------------------------------------- */
  /* The runs, from the worker.                                                                                  */
  /* ---------------------------------------------------------------------------------------------------------- */

  let worker: Worker | null = null;
  /** Why the runs cannot come, once the worker could not start or has failed. */
  let failed = "";
  function ask(): void {
    if (failed) return;
    if (!worker) {
      try {
        worker = new Worker(new URL("../sim/tripRuns.worker.ts", import.meta.url), { type: "module" });
      } catch {
        failed = "This browser could not start the planner for the kitchen.";
        renderStatus();
        return;
      }
      worker.onerror = () => {
        failed = "The planner for the kitchen stopped; reload the page to start it again.";
        renderStatus();
      };
      worker.onmessage = (e: MessageEvent<TripDone>) => {
        const d = e.data;
        runs.set(d.key, d);
        // The oldest go first, never the one on screen or the one asked for.
        for (const k of runs.keys()) {
          if (runs.size <= PLANS_KEPT) break;
          if (k !== shownKey && k !== tripKey(settings)) runs.delete(k);
        }
        if (d.key === tripKey(settings)) apply();
      };
    }
    if (!runs.has(tripKey(settings))) worker.postMessage({ s: { ...settings } } satisfies TripAsk);
  }

  let endMin = DOOR + 6;
  let loopS = endMin * SEC_PER_MIN + END_HOLD_S;
  /** What the twin's kitchen was last planned with. */
  let twinFor: { plan: Morning; endMin: number; loop: number } | null = null;
  let tMin = 0;
  let loopT = 0;
  /** Put the plan for the settings on screen, if it has come, with the morning the trip makes of it. */
  function apply(): void {
    const key = tripKey(settings);
    const r = runs.get(key);
    if (r) {
      // Kept as the newest.
      runs.delete(key);
      runs.set(key, r);
      const frac = endMin > 0 ? tMin / endMin : 0;
      shown = r;
      shownKey = key;
      actual = tripMorning(r.plan, trip);
      shownTrip = { ...trip };
      let last = 0;
      for (const m of [r.plan, actual]) for (const j of m.jobs) last = Math.max(last, j.t1);
      // On past the last work, and past each van until its carriers are out of the door.
      endMin = Math.max(DOOR + 6, last + 6, r.plan.depart + ROLL_MIN + 2, actual.depart + ROLL_MIN + 2);
      let loop = kitchen.plan(actual, actual.settings.cooks, endMin, endMin * SEC_PER_MIN + END_HOLD_S);
      let twinLoop = twinFor && twinFor.plan === r.plan && twinFor.endMin === endMin && Math.abs(twinFor.loop - loop) < 1e-6 ? twinFor.loop : twin.plan(r.plan, r.plan.settings.cooks, endMin, loop);
      if (twinLoop > loop + 1e-6) loop = kitchen.plan(actual, actual.settings.cooks, endMin, twinLoop);
      // The twin's walks home can grow its loop as well: then both run to the longer.
      if (loop > twinLoop + 1e-6) twinLoop = twin.plan(r.plan, r.plan.settings.cooks, endMin, loop);
      twinFor = { plan: r.plan, endMin, loop: twinLoop };
      loopS = Math.max(loop, twinLoop);
      tMin = frac * endMin;
      loopT = tMin * SEC_PER_MIN;
      placeBox(shownTrip.oven);
    }
    renderControls();
    renderMap();
    renderTiles();
    renderStatus();
    lastHud = "";
  }
  let pending = 0;
  /** The settings changed: ask for their runs (a slider's steps only once it rests a moment). */
  function changed(now: boolean): void {
    window.clearTimeout(pending);
    const go = () => {
      ask();
      apply();
    };
    if (now) go();
    else {
      renderControls();
      renderMap();
      pending = window.setTimeout(go, 150);
    }
  }
  for (const b of ovenBtns)
    b.addEventListener("click", () => {
      trip.oven = Number(b.dataset.tpOven);
      changed(true);
    });
  atIn?.addEventListener("input", () => {
    trip.at = Number(atIn.value);
    changed(false);
  });
  forIn?.addEventListener("input", () => {
    trip.mins = Number(forIn.value);
    changed(false);
  });
  mealsIn?.addEventListener("input", () => {
    settings.meals = Number(mealsIn.value);
    changed(false);
  });

  /* The readouts on the stage and the screens, each minute. */
  let lastHud = "";
  const cctvTime = (t: number) => {
    const m = Math.max(0, Math.floor(t));
    return `CAM 1 · ${String(6 + Math.floor(m / 60)).padStart(2, "0")}:${String(m % 60).padStart(2, "0")}`;
  };
  function showHud(t: number): void {
    const plan = shown?.plan ?? null;
    const m = actual;
    const k = Math.floor(t);
    const tr = shownTrip;
    const out = t >= tr.at && t < tr.at + tr.mins;
    const seen = t >= tr.at;
    let van = "…";
    let note = "";
    let kind = "idle";
    if (plan && m) {
      const c = costOf(m, plan);
      // From the minute the camera shows the combi dark, the twin's forecast is the morning with the trip: while it
      // is out, if it is back when picked.
      const forecast = seen ? m : plan;
      van = !forecast.finished ? `Van: not by ${CUTOFF}` : `Van ${clockOf(forecast.depart)}${!seen || free(c) ? " as planned" : out ? ` if back at ${clockOf(tr.at + tr.mins)}` : ""}`;
      if (t >= m.depart && m.finished) {
        note = `Van gone at ${clockOf(m.depart)}`;
        kind = m.late > 0.5 ? "late" : "held";
      } else if (out) {
        note = `${cap(OVEN_NAME[tr.oven]!)} out until ${clockOf(tr.at + tr.mins)}`;
        kind = "full";
      } else if (t >= DOOR) {
        note = "The van is waiting";
        kind = "late";
      } else {
        // Waits for a combi only: a cook waiting for a range or the table is the menu's, not the trip's.
        const w = m.waits.find((x) => x.at === "oven" && x.t0 <= t && x.t1 > t);
        if (w) {
          note = `${w.cooks.length === 1 ? "A cook" : `${w.cooks.length} cooks`} waiting for a combi`;
          kind = "wait";
        }
      }
    }
    const sig = `${k}|${van}|${note}`;
    if (sig === lastHud) return;
    lastHud = sig;
    if (hudClock) hudClock.textContent = clockOf(t);
    if (hudVan) hudVan.textContent = van;
    if (hudNote) {
      hudNote.textContent = note;
      hudNote.hidden = !note;
      hudNote.dataset.kind = kind;
    }
    for (const el of [tsEl.twin, tsEl.seen]) if (el) el.textContent = cctvTime(t);
    // The camera's screen: as planned until the combi goes dark, then the twin's mark on what differs.
    const differs = !!plan && seen;
    seenMon?.toggleAttribute("data-differs", differs);
    if (boxEl) boxEl.hidden = !out;
    if (flagEl) {
      flagEl.hidden = !differs;
      if (differs && plan && m) {
        const c = costOf(m, plan);
        flagEl.textContent = out ? `Not as planned: ${OVEN_NAME[tr.oven]} is dark` : free(c) ? "Combi back: the van as planned" : `Combi back: van ${m.finished ? clockOf(m.depart) : `not by ${CUTOFF}`}`;
      }
    }
  }

  /** The red floor and the badges at minute t, `dt` seconds after the last frame. */
  function showMarks(t: number, dt: number): void {
    const m = actual;
    outPads.forEach((p, o) => {
      const want = m && o === shownTrip.oven && t >= shownTrip.at && t < shownTrip.at + shownTrip.mins ? 1 : 0;
      p.on += (want - p.on) * Math.min(1, dt * 6);
      if (p.on < 0.01 && want === 0) p.on = 0;
      p.mesh.visible = p.on > 0;
      p.mat.opacity = 0.62 * p.on;
    });
    // A cook free while work waits for a combi: their badge, where the routine has them now.
    const waits = m ? m.waits.filter((w) => w.at === "oven" && w.t0 <= t && w.t1 > t) : [];
    kitchen.cooks.forEach((c, i) => {
      const b = badges[i]!;
      const on = c.st.group.visible && !!c.routine && waits.some((w) => w.cooks.includes(kitchen.whoAt(i, loopT)));
      b.visible = on;
      if (on) {
        const at = c.routine!.rootAt(loopT);
        b.position.set(at.x, c.st.dims.standing_height_m + 0.3, at.z);
      }
    });
    // The camera's red lamp: on while it records, blinking on its own clock, once a second.
    rig.led.color.setHex(Math.floor(performance.now() / 500) % 2 === 0 ? 0xff2b2b : 0x3a0a0a);
  }

  renderControls();
  apply();
  ask();
  if (mapEl) {
    let lastW = mapEl.clientWidth;
    new ResizeObserver(() => {
      if (Math.abs(mapEl.clientWidth - lastW) < 4) return;
      lastW = mapEl.clientWidth;
      renderMap();
    }).observe(mapEl);
  }
  new ResizeObserver(() => placeBox(shownTrip.oven)).observe(screenEl.seen);
  for (const k of kitchen.cooks) void k.st.ready.then(() => void gate.refresh(scene).then(() => void seenGate.refresh(scene)));
  for (const k of twin.cooks) void k.st.ready.then(() => void twinGate.refresh(twinScene));

  let inView = false;
  observeInView(stageEl, (v) => {
    inView = v;
  });
  let monsInView = false;
  observeInView(monsEl, (v) => {
    monsInView = v;
  });
  let frozen: number | null = null;
  let last = performance.now();
  let screenT = 0;
  let screenNext: "twin" | "seen" = "twin";
  const fitClip = clipFitter(camera, scene, 0.05);
  const fitCctv = clipFitter(cctv, scene, 0.05);
  function animate(now: number): void {
    requestAnimationFrame(animate);
    const dt = Math.max(0, Math.min(0.05, (now - last) / 1000));
    last = now;
    if ((!inView && !monsInView) || document.hidden || !gate.ready()) return;
    if (frozen === null) loopT = (loopT + dt) % loopS;
    else loopT = frozen * SEC_PER_MIN;
    tMin = Math.min(endMin, loopT / SEC_PER_MIN);
    kitchen.pose(loopT);
    const a = kitchen.show(tMin, dt);
    showMarks(tMin, dt);
    if (a.minute || !lastHud) showHud(tMin);
    if (inView) {
      controls.update();
      fitClip();
      renderer.render(scene, camera);
    }
    // The camera's screens, one then the other, a few times a second; neither shows what only the stage marks.
    screenT += dt;
    if (monsInView && (screenT >= SCREEN_S || frozen !== null)) {
      screenT = 0;
      fitCctv();
      if (screenNext === "twin" || frozen !== null) {
        twin.pose(loopT);
        twin.show(tMin, SCREEN_S * 2);
        if (twinGate.ready()) rTwin!.render(twinScene, cctv);
      }
      if ((screenNext === "seen" || frozen !== null) && seenGate.ready()) {
        const fog = scene.fog as THREE.Fog;
        const near = fog.near;
        const far = fog.far;
        const twinFog = twinScene.fog as THREE.Fog;
        fog.near = twinFog.near;
        fog.far = twinFog.far;
        marks.visible = false;
        rSeen!.render(scene, cctv);
        marks.visible = true;
        fog.near = near;
        fog.far = far;
      }
      screenNext = screenNext === "twin" ? "seen" : "twin";
    }
  }
  requestAnimationFrame(animate);
  void gate.refresh(scene).then(() => void twinGate.refresh(twinScene).then(() => void seenGate.refresh(scene)));

  // A hook for checking the page by script.
  Object.assign(stageEl, {
    trip: {
      info: () => {
        const plan = shown?.plan ?? null;
        const free0 = shown ? freeShare(shown.map) : null;
        return {
          mountMs: Math.round(mountMs),
          key: tripKey(settings),
          ready: ready(),
          trip: { ...trip },
          plan: plan ? { depart: clockOf(plan.depart), made: madeIt(plan) } : null,
          depart: actual ? (actual.finished ? clockOf(actual.depart) : null) : null,
          late: actual && plan ? Math.round(costOf(actual, plan).late) : null,
          pastBest: actual && plan ? costOf(actual, plan).pastBest : null,
          free: free0 ? free0.n : null,
          cells: free0 ? free0.of : null,
          cooks: kitchen.cooks.filter((k) => k.st.group.visible).length,
          twinCooks: twin.cooks.filter((k) => k.st.group.visible).length,
          pads: outPads.filter((p) => p.mesh.visible).length,
          badges: badges.filter((b) => b.visible).length,
          box: boxEl ? !boxEl.hidden : null,
          flag: flagEl && !flagEl.hidden ? flagEl.textContent : null,
          loopS: +loopS.toFixed(2),
          t: +tMin.toFixed(1),
        };
      },
      hold(min: number | null) {
        frozen = min;
      },
      view(eye: [number, number, number], at: [number, number, number]) {
        userMoved = true;
        camera.position.set(...eye);
        controls.target.set(...at);
        camera.clearViewOffset();
        camera.updateProjectionMatrix();
        controls.update();
      },
      set(s: Partial<TripSettings & { oven: number; at: number; mins: number }>) {
        if (s.menu !== undefined) settings.menu = s.menu;
        if (s.meals !== undefined) settings.meals = s.meals;
        if (s.oven !== undefined) trip.oven = s.oven;
        if (s.at !== undefined) trip.at = s.at;
        if (s.mins !== undefined) trip.mins = s.mins;
        if (mealsIn) mealsIn.value = String(settings.meals);
        changed(true);
      },
    },
  });
  const mountMs = performance.now() - mountT0;
}
