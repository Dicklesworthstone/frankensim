/**
 * One cook, two dishes (chapter 07, first figure): a cook wearing camera glasses makes two dishes in the morning's
 * kitchen. The twin keeps each dish as a card (its steps, where each happens, how long, and how much of that time it
 * needs the cook's hands) and plans the two together for one pair of hands (sim/oneCook.ts).
 *
 * The stage acts out the plan picked: both at once, or one then the other. Beside it, the glasses' own picture, from
 * the camera at the corner of the frame. Under them, the two cards, each step lit the first time the glasses see it;
 * and a chart of the cook's hands over both plans, the work of each dish stacked as much of the cook as it takes.
 */
import * as THREE from "three";
import { RoundedBoxGeometry } from "three/examples/jsm/geometries/RoundedBoxGeometry.js";
import { clockOf, type Morning, type MorningDish } from "../sim/cateringMorning";
import { COOK_DISHES, cardOf, oneCook, type CardStep, type OneCook } from "../sim/oneCook";
import { bindCanvasResize, clipFitter, createScrollSafeOrbit, createStudioRenderer, observeInView, shaderGate } from "./kit";
import { CARRIER, CARRIER_X0, END_HOLD_S, LINE_Z, MORNING_PIECES, MorningKitchen, OVEN, OVEN_X0, RANGE, SEC_PER_MIN, SLOT_COLOR, TABLE, TABLE_Z, rangeX0, shortName, type CookSpec } from "./morningKitchen";
import type { Piece } from "./models/all";
import { studioEnvironment, studioLights } from "./models/stage";
import { backdrop, tilePlane } from "./weekLab";

/** One kitchen's pieces, for main.ts to build ahead. */
export const GLASSES_PIECES: readonly Piece[] = MORNING_PIECES;

/** The cook with the glasses. */
const COOK: CookSpec[] = [{ cohort: "p50", seed: 41, apron: 0x26344d }];
const num = (v: number) => Math.round(v).toLocaleString("en-US");
const cap = (s: string) => s.replace(/^./, (c) => c.toUpperCase());
/** "1 h 10" or "40 min". */
const span = (m: number) => {
  const r = Math.round(m);
  return r >= 60 ? `${Math.floor(r / 60)} h${r % 60 ? ` ${String(r % 60).padStart(2, "0")}` : ""}` : `${r} min`;
};
/** A dish's short name, lower case: "thighs", "mac and cheese". */
const dishName = (d: MorningDish) => shortName(d.name).toLowerCase();
/** Whether a dish's short name takes a plural verb ("the thighs roast", "the braise bakes"). */
const plural = (d: MorningDish) => /s$/.test(dishName(d));
/** The glasses' camera: its field of view up and down (a portrait picture, 3 by 4), and how far it looks down. */
const POV_FOV = 76;
const POV_TILT = THREE.MathUtils.degToRad(6);
/** Seconds between frames of the glasses' picture. */
const POV_S = 1 / 30;

/**
 * The stage's clock. At the morning's pace a step of a minute or two would be over before the cook got to it, so each
 * step is slowed to at least STEP_MIN minutes of the stage's clock; and the plan's clock stands still while the cook
 * walks from one station to the next (the plan counts no time for it): the stage's cook sets off for a station as its
 * work starts (routineFor), so each walk is put at the start of the stretch it walks to. The rest runs at the
 * morning's pace, or faster where the loop would run past LOOP_S. `at` takes a minute of the plan to the stage's (a
 * station's first minute to the cook setting off for it), `back` the stage's to the plan's, and `walk` says where the
 * cook is walking to at a minute of the stage's.
 */
const STEP_MIN = 10;
const LOOP_S = 75;
/** A walk's seconds, until the stage has walked it (routineFor's walks are the same every time). */
const WALK_GUESS_S = 5;
interface Warp {
  at(t: number): number;
  back(u: number): number;
  walk(u: number): string | null;
  /** The walks it took, by their routine's key ("table:0>oven:0"), and the seconds it gave each. */
  walked: Map<string, number>;
}
interface Place {
  key: string;
  name: string;
}
/** Where the plan has the cook at minute `t`, as the stage puts them (rawActs in morningKitchen): at the table or a
 * combi for work that needs the whole cook, else at the range. */
function placeAt(m: Morning, t: number): Place | null {
  const on = m.jobs.filter((j) => j.t0 <= t && j.t1 > t);
  if (!on.length) return null;
  const whole = on.find((j) => j.hands >= 1 - 1e-9 && j.at !== "range");
  if (whole?.at === "table") return { key: `table:${whole.spot}`, name: "the table" };
  if (whole) return { key: `oven:${whole.spot}`, name: whole.spot === 0 ? "the left combi" : "the right combi" };
  return { key: `range:${on[0]!.side ?? 0}`, name: "the range" };
}
/** The stage's clock for a plan to minute `end`, each walk given `walkS(key)` seconds. */
function warpOf(m: Morning, end: number, walkS: (key: string) => number): Warp {
  const raw = [0, end];
  for (const j of m.jobs) raw.push(j.t0, j.t1);
  for (const l of m.loads) raw.push(l.t0, l.t1, l.t2, l.tu, l.t3);
  const ts: number[] = [];
  for (const t of raw.filter((t) => t >= 0 && t <= end).sort((a, b) => a - b)) if (!ts.length || t - ts[ts.length - 1]! > 1e-6) ts.push(t);
  // Each stretch between two cuts: a walk (stage minutes), or work slowed or not (minutes of the plan, and the factor).
  const parts: { a: number; b: number; f: number; walk?: string; g?: number }[] = [];
  const walked = new Map<string, number>();
  let last: Place | null = null;
  for (let i = 0; i + 1 < ts.length; i++) {
    const [a, b] = [ts[i]!, ts[i + 1]!];
    const mid = (a + b) / 2;
    const here = placeAt(m, mid);
    if (here && last && here.key !== last.key) {
      const key = `${last.key}>${here.key}`;
      const g = walkS(key) / SEC_PER_MIN;
      walked.set(key, g * SEC_PER_MIN);
      parts.push({ a, b: a, f: 0, walk: here.name, g });
    }
    if (here) last = here;
    let f = 1;
    for (const j of m.jobs) if (j.t0 <= mid && j.t1 > mid) f = Math.max(f, STEP_MIN / Math.max(1e-6, j.t1 - j.t0));
    parts.push({ a, b, f });
  }
  // Past LOOP_S, the plain stretches go faster (down to a quarter of the morning's pace).
  const fixed = parts.reduce((s, p) => s + (p.walk ? p.g! : p.f > 1 ? (p.b - p.a) * p.f : 0), 0);
  const plain = parts.reduce((s, p) => s + (!p.walk && p.f <= 1 ? p.b - p.a : 0), 0);
  const room = LOOP_S / SEC_PER_MIN - fixed;
  const c = plain > 0 && fixed + plain > LOOP_S / SEC_PER_MIN ? Math.max(0.25, room / plain) : 1;
  const xs = [ts[0]!];
  const us = [0];
  const walks: [number, number, string][] = [];
  for (const p of parts) {
    const u = us[us.length - 1]!;
    if (p.walk) walks.push([u, u + p.g!, p.walk]);
    xs.push(p.b);
    us.push(u + (p.walk ? p.g! : (p.b - p.a) * (p.f > 1 ? p.f : c)));
  }
  /** Piecewise straight from xs (rising, two equal at a walk) to ys, the first of equal xs, and at the morning's pace
   * past either end. */
  const map = (xs: number[], ys: number[], v: number) => {
    if (v <= xs[0]!) return ys[0]! + (v - xs[0]!);
    const n = xs.length - 1;
    if (v > xs[n]!) return ys[n]! + (v - xs[n]!);
    // The first k with xs[k] >= v, then straight from the one before.
    let lo = 0;
    let hi = n;
    while (lo < hi) {
      const k = (lo + hi) >> 1;
      if (xs[k]! < v) lo = k + 1;
      else hi = k;
    }
    if (xs[lo] === v) return ys[lo]!;
    const d = xs[lo]! - xs[lo - 1]!;
    return ys[lo - 1]! + ((v - xs[lo - 1]!) / d) * (ys[lo]! - ys[lo - 1]!);
  };
  return { at: (t) => map(xs, us, t), back: (u) => map(us, xs, u), walk: (u) => walks.find(([a, b]) => u >= a && u < b)?.[2] ?? null, walked };
}
/** A plan on the stage's clock: its every time taken through `at` (-1, never, stays). */
function onStage(m: Morning, at: (t: number) => number): Morning {
  const w = (x: number) => (x < 0 ? x : at(x));
  return {
    ...m,
    jobs: m.jobs.map((j) => ({ ...j, t0: w(j.t0), t1: w(j.t1) })),
    burns: m.burns.map((x) => ({ ...x, t0: w(x.t0), t1: w(x.t1) })),
    loads: m.loads.map((x) => ({ ...x, t0: w(x.t0), t1: w(x.t1), t2: w(x.t2), tu: w(x.tu), t3: w(x.t3) })),
    waits: m.waits.map((x) => ({ ...x, t0: w(x.t0), t1: w(x.t1) })),
    batches: m.batches.map((x) => ({ ...x, release: w(x.release), started: w(x.started), cooked: w(x.cooked), packed: w(x.packed) })),
    lastPacked: m.lastPacked.map(w),
    depart: w(m.depart),
  };
}

/**
 * Camera glasses, built in the head bone's frame (+z forward, +y up, +x the wearer's left): a black acetate front with
 * two rounded lenses, temples back over the ears, a camera in the left-hand corner piece and a capture lamp in the
 * right. `eye`: the middle of the wearer's left eye's front, in that frame; `half`: half the head's width at the temples.
 */
function cameraGlasses(eye: THREE.Vector3, half: number): { group: THREE.Group; lens: THREE.Object3D; led: THREE.MeshBasicMaterial } {
  const group = new THREE.Group();
  const acetate = new THREE.MeshPhysicalMaterial({ color: 0x0b0b0c, roughness: 0.28, clearcoat: 0.8, clearcoatRoughness: 0.18 });
  const tint = new THREE.MeshPhysicalMaterial({ color: 0x2c3138, roughness: 0.05, metalness: 0, transparent: true, opacity: 0.32, clearcoat: 1, clearcoatRoughness: 0.03, depthWrite: false });
  const glass = new THREE.MeshPhysicalMaterial({ color: 0x050608, roughness: 0.06, clearcoat: 1, clearcoatRoughness: 0.02 });
  const led = new THREE.MeshBasicMaterial({ color: 0xf5f7ff, toneMapped: false });
  const mesh = (g: THREE.BufferGeometry, m: THREE.Material, shadow = true) => {
    const x = new THREE.Mesh(g, m);
    x.castShadow = shadow;
    return x;
  };
  /** A rounded rectangle, centred. */
  const rr = <P extends THREE.Path>(w: number, h: number, r: number, path: P): P => {
    path.moveTo(-w / 2 + r, -h / 2);
    path.lineTo(w / 2 - r, -h / 2);
    path.quadraticCurveTo(w / 2, -h / 2, w / 2, -h / 2 + r);
    path.lineTo(w / 2, h / 2 - r);
    path.quadraticCurveTo(w / 2, h / 2, w / 2 - r, h / 2);
    path.lineTo(-w / 2 + r, h / 2);
    path.quadraticCurveTo(-w / 2, h / 2, -w / 2, h / 2 - r);
    path.lineTo(-w / 2, -h / 2 + r);
    path.quadraticCurveTo(-w / 2, -h / 2, -w / 2 + r, -h / 2);
    return path;
  };
  // The front sits a little ahead of the eyes, its lenses' middles a touch above them.
  const LENS = { w: 0.05, h: 0.041, r: 0.014 };
  const RIM = 0.0055;
  const DEPTH = 0.0065;
  const z0 = eye.z + 0.013;
  const y0 = eye.y + 0.002;
  const cx = Math.max(eye.x, LENS.w / 2 + 0.009);
  for (const s of [1, -1]) {
    const shape = rr(LENS.w + 2 * RIM, LENS.h + 2 * RIM, LENS.r + RIM, new THREE.Shape());
    shape.holes.push(rr(LENS.w, LENS.h, LENS.r, new THREE.Path()));
    const rim = mesh(new THREE.ExtrudeGeometry(shape, { depth: DEPTH, bevelEnabled: true, bevelThickness: 0.0012, bevelSize: 0.0009, bevelSegments: 2, curveSegments: 10 }), acetate);
    rim.position.set(s * cx, y0, z0 - DEPTH / 2);
    group.add(rim);
    const lens = mesh(new THREE.ShapeGeometry(rr(LENS.w + 0.001, LENS.h + 0.001, LENS.r, new THREE.Shape()), 10), tint, false);
    lens.position.set(s * cx, y0, z0 + 0.0004);
    lens.renderOrder = 2;
    group.add(lens);
  }
  // The bridge, set back a little within the rims' depth so no face lies on theirs.
  const gap = 2 * cx - LENS.w - 2 * RIM;
  const bridge = mesh(new RoundedBoxGeometry(gap + 2 * RIM, 0.007, DEPTH * 0.7, 2, 0.002), acetate);
  bridge.position.set(0, y0 + LENS.h * 0.22, z0 - 0.0004);
  group.add(bridge);
  // The corner pieces: the frame's ends, thicker, where the temples hinge; the camera in the left one (clear of a
  // right hand at work), the lamp in the right.
  const endX = cx + LENS.w / 2 + RIM;
  const ends = [1, -1].map((s) => {
    const e = mesh(new RoundedBoxGeometry(0.011, 0.016, 0.016, 2, 0.0035), acetate);
    e.position.set(s * (endX + 0.002), y0 + LENS.h * 0.18, z0 - 0.004);
    group.add(e);
    return e;
  });
  const lensObj = new THREE.Object3D();
  {
    const e = ends[0]!;
    const barrel = mesh(new THREE.CylinderGeometry(0.0034, 0.0034, 0.002, 20).rotateX(Math.PI / 2), glass, false);
    barrel.position.set(e.position.x, e.position.y, e.position.z + 0.0082);
    group.add(barrel);
    lensObj.position.set(e.position.x, e.position.y, e.position.z + 0.0095);
    group.add(lensObj);
    const l = ends[1]!;
    const lamp = mesh(new THREE.CircleGeometry(0.0013, 14), led, false);
    lamp.position.set(l.position.x, l.position.y, l.position.z + 0.0081);
    group.add(lamp);
  }
  // The temples: from each corner piece straight back along the side of the head, then down behind the ear.
  for (const s of [1, -1]) {
    const x = Math.max(endX + 0.002, half + 0.004);
    const a = new THREE.Vector3(s * (endX + 0.004), y0 + LENS.h * 0.18, z0 - 0.012);
    const b = new THREE.Vector3(s * x, y0 + LENS.h * 0.1, eye.z - 0.085);
    const c = new THREE.Vector3(s * (x - 0.006), y0 - 0.03, eye.z - 0.112);
    const curve = new THREE.CatmullRomCurve3([a, a.clone().lerp(b, 0.5), b, c], false, "centripetal");
    group.add(mesh(new THREE.TubeGeometry(curve, 40, 0.0026, 8, false), acetate));
  }
  return { group, lens: lensObj, led };
}

export function mountGlassesLab(root: HTMLElement): void {
  const mountT0 = performance.now();
  const stageEl = root.querySelector<HTMLElement>("[data-gl-stage]");
  const canvas = root.querySelector<HTMLCanvasElement>("[data-gl-canvas]");
  const povCanvas = root.querySelector<HTMLCanvasElement>("[data-gl-pov]");
  const povEl = root.querySelector<HTMLElement>(".gl-pov");
  if (!stageEl || !canvas || !povCanvas || !povEl) return;
  const maybeRenderer = createStudioRenderer(canvas, false);
  if (!maybeRenderer) return;
  const renderer = maybeRenderer;
  renderer.toneMappingExposure = 1.0;
  const rPov = createStudioRenderer(povCanvas, false);
  if (!rPov) return;
  rPov.toneMappingExposure = 1.0;
  const phone = window.matchMedia("(max-width: 720px)").matches;

  /* The scene: the morning's kitchen, one cook in it. */
  const scene = new THREE.Scene();
  scene.background = backdrop();
  scene.fog = new THREE.Fog(0x151510, 8, 16);
  scene.environment = studioEnvironment(renderer);
  scene.environmentIntensity = 0.6;
  const lights = studioLights(scene, new THREE.Vector3(-0.4, 0.9, 0.1), 5.2, phone ? 1024 : 2048);
  lights.key.position.set(-1.4, 6.4, 4.6);
  lights.key.intensity = 1.9;
  lights.key.shadow.normalBias = 0.02;
  lights.rim.position.set(3.4, 3.2, -3.0);
  lights.rim.intensity = 0.75;
  lights.hemi.intensity = 0.5;
  {
    const floor = tilePlane(16, 10);
    floor.rotation.x = -Math.PI / 2;
    floor.position.set(0.2, 0, 0.6);
    scene.add(floor);
  }
  const kitchen = new MorningKitchen(COOK);
  kitchen.rollOut = false;
  scene.add(kitchen.group);
  const cook = kitchen.cooks[0]!;

  /* The glasses, on the cook's head: the eyes' place in the head's frame, measured on the drawn body. */
  const head = cook.st.bone("head");
  /** Where the eyes are in the head bone's frame, and half the head's width at the temples. */
  const face = { eye: new THREE.Vector3(0.032, 0.07, 0.1), half: 0.071, measured: false };
  let glasses = cameraGlasses(face.eye, face.half);
  head.add(glasses.group);
  /** Measure the face on the drawn body (its bind pose is the body's rest pose: identity binds), and refit. */
  const measureFace = () => {
    const skel = cook.st.skeleton;
    const hi = skel.bones.indexOf(head);
    if (hi < 0) return;
    // The meshes are bound at the body's rest pose with identity binds: a vertex is at rest in the body's frame, and
    // the head bone's inverse takes it into the head's.
    const toHead = skel.boneInverses[hi]!;
    const eyeY = cook.st.dims.eye_height_m + cook.st.dims.shoe_sole_m;
    let eyeZ = -Infinity;
    let half = 0;
    const p = new THREE.Vector3();
    for (const m of cook.st.meshes) {
      const pos = m.geometry.getAttribute("position") as THREE.BufferAttribute | undefined;
      if (!pos) continue;
      for (let i = 0; i < pos.count; i++) {
        p.fromBufferAttribute(pos, i);
        if (Math.abs(p.y - eyeY) < 0.008 && Math.abs(Math.abs(p.x) - 0.032) < 0.008) eyeZ = Math.max(eyeZ, p.z);
        if (p.y > eyeY - 0.005 && p.y < eyeY + 0.03 && p.z < 0.02 && p.z > -0.06) half = Math.max(half, Math.abs(p.x));
      }
    }
    if (!Number.isFinite(eyeZ) || half <= 0) return;
    face.eye.set(0.032, eyeY, eyeZ).applyMatrix4(toHead);
    face.eye.x = Math.abs(face.eye.x);
    face.half = half;
    face.measured = true;
    head.remove(glasses.group);
    glasses = cameraGlasses(face.eye, face.half);
    head.add(glasses.group);
  };
  void cook.st.ready.then(() => {
    measureFace();
    void gate.refresh(scene).then(() => void povGate.refresh(scene));
  });

  /* The stage's camera: from in front of the worktable, raised, the line and the table in view. */
  const camera = new THREE.PerspectiveCamera(32, 16 / 9, 0.05, 60);
  const TARGET = new THREE.Vector3();
  const VIEW = new THREE.Vector3();
  const keep: THREE.Vector3[] = [];
  const LEFT = OVEN_X0[0];
  const RIGHT = rangeX0(0) + RANGE.w;
  const aim = () => {
    keep.length = 0;
    const wide = stageEl.clientWidth >= 1.3 * stageEl.clientHeight;
    TARGET.set((LEFT + RIGHT) / 2, 0.75, 0.1);
    VIEW.set(0.38, wide ? 0.52 : 0.8, 1).normalize();
    // The combis' tops, the range's back, the table's front corners, and the first carriers by the door.
    keep.push(new THREE.Vector3(LEFT, OVEN.h, LINE_Z), new THREE.Vector3(LEFT, OVEN.h, LINE_Z - OVEN.d), new THREE.Vector3(RIGHT, RANGE.h, LINE_Z - RANGE.d));
    keep.push(new THREE.Vector3(TABLE.x0, 0, TABLE_Z + TABLE.d / 2), new THREE.Vector3(TABLE.x0 + TABLE.w, 0, TABLE_Z + TABLE.d / 2));
    keep.push(new THREE.Vector3(CARRIER_X0 + CARRIER.w, 0, TABLE_Z + CARRIER.d / 2), new THREE.Vector3(CARRIER_X0 + CARRIER.w, 2 * CARRIER.h, TABLE_Z));
  };
  const controls = createScrollSafeOrbit(camera, canvas, [0, 0.8, 0], { minDistance: 1.2, maxDistance: 16, maxPolarAngle: Math.PI * 0.47 });
  let userMoved = false;
  controls.addEventListener("start", () => {
    userMoved = true;
  });
  const ndc = new THREE.Vector3();
  const hudEl = stageEl.querySelector<HTMLElement>(".gl-hud");
  const legendEl = stageEl.querySelector<HTMLElement>(".gl-legend");
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
    // The kitchen fits between the HUD and the legend, centred across.
    const top = Math.min(0.86, 1 - (2 * (hudMax + 6)) / h);
    const bot = -Math.min(0.92, 1 - (2 * ((legendEl ? h - legendEl.offsetTop : 0) + 4)) / h);
    let dist = 2.5;
    let lo = 0;
    let hi = 0;
    let left = 0;
    let right = 0;
    for (; dist < 22; dist += 0.05) {
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
  /* The glasses' camera: a portrait picture from the corner of the frame, looking where the head looks. */
  const pov = new THREE.PerspectiveCamera(POV_FOV, 3 / 4, 0.012, 30);
  const povGate = shaderGate(rPov, pov);
  bindCanvasResize(stageEl, renderer, camera, gate);
  const povScreen = povCanvas.parentElement!;
  bindCanvasResize(povScreen, rPov, pov, povGate);
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
  const FLIP = new THREE.Quaternion().setFromAxisAngle(new THREE.Vector3(0, 1, 0), Math.PI).multiply(new THREE.Quaternion().setFromAxisAngle(new THREE.Vector3(1, 0, 0), -POV_TILT));
  const placePov = () => {
    glasses.lens.updateWorldMatrix(true, false);
    glasses.lens.getWorldPosition(pov.position);
    glasses.lens.getWorldQuaternion(pov.quaternion).multiply(FLIP);
    pov.updateMatrixWorld(true);
  };

  /* ---------------------------------------------------------------------------------------------------------- */
  /* The page.                                                                                                   */
  /* ---------------------------------------------------------------------------------------------------------- */

  const q = <T extends Element = HTMLElement>(s: string) => root.querySelector<T>(s);
  const dishGroups = [0, 1].map((k) => q(`[data-gl-dishes="${k}"]`));
  const portionsIn = q<HTMLInputElement>("[data-gl-portions]");
  const portionsVal = q("[data-gl-portions-val]");
  const planBtns = Array.from(root.querySelectorAll<HTMLButtonElement>("[data-gl-plan]"));
  const statusEl = q("[data-gl-status]");
  const cardsEl = q("[data-gl-cards]");
  const chartEl = q("[data-gl-chart]");
  const chartNote = q("[data-gl-chart-note]");
  const keyEls = [q("[data-gl-key-a]"), q("[data-gl-key-b]")];
  const hudClock = q("[data-gl-hud-clock]");
  const hudDoing = q("[data-gl-hud-doing]");
  const hudNote = q("[data-gl-hud-note]");
  const tsEl = q("[data-gl-ts]");
  const seenEl = q("[data-gl-seen]");
  const tagEl = q("[data-gl-tag]");
  const tile = (k: string) => ({ v: q(`[data-gl-${k}]`), s: q(`[data-gl-${k}-s]`) });
  const tiles = { both: tile("both"), apart: tile("apart"), saved: tile("saved"), free: tile("free") };
  const setTile = (t: (typeof tiles)[keyof typeof tiles], v: string, s: string) => {
    if (t.v) t.v.textContent = v;
    if (t.s) t.s.textContent = s;
  };

  /** The two dishes (indices into COOK_DISHES), the lunches of each, and the plan the stage plays. */
  const pick: [number, number] = [1, 2];
  let portions = Number(portionsIn?.value ?? 60);
  let playing: "together" | "apart" = "together";
  let run: OneCook = oneCook(COOK_DISHES[pick[0]]!, COOK_DISHES[pick[1]]!, portions);
  const shown = () => (playing === "together" ? run.together : run.apart);

  /* The dish buttons: one group for each of the two dishes; a dish picked in one is not offered in the other. */
  const dishBtns = dishGroups.map((g, k) => {
    if (!g) return [] as HTMLButtonElement[];
    g.innerHTML = "";
    return COOK_DISHES.map((d, i) => {
      const b = document.createElement("button");
      b.type = "button";
      b.textContent = shortName(d.name);
      b.addEventListener("click", () => {
        if (pick[1 - k] === i) return;
        pick[k] = i;
        changed(true);
      });
      g.appendChild(b);
      return b;
    });
  });
  function renderControls(): void {
    dishBtns.forEach((bs, k) =>
      bs.forEach((b, i) => {
        b.setAttribute("aria-pressed", String(pick[k] === i));
        b.disabled = pick[1 - k] === i;
      }),
    );
    for (const b of planBtns) b.setAttribute("aria-pressed", String(b.dataset.glPlan === playing));
    if (portionsIn) portionsIn.value = String(portions);
    if (portionsVal) portionsVal.textContent = `${num(portions)} lunches`;
    keyEls.forEach((el, k) => {
      if (el) el.textContent = shortName(run.dishes[k]!.name);
    });
  }

  /* The cards: each dish's steps as the twin keeps them, lit as the glasses see them. */
  /** "3.6 s a portion", "8 min a skillet of 8", "35 min in a combi". */
  const perText = (d: MorningDish, c: CardStep) => {
    if (c.per === "portion") return `${c.min < 1 ? `${+(c.min * 60).toFixed(1)} s` : `${+c.min.toFixed(1)} min`} a portion`;
    if (c.per === "batch") {
      const st = d.stages.find((s) => s.verb === c.verb);
      return `${c.min} min a ${st?.vessel ?? "pot"} of ${d.batch}`;
    }
    return `${span(c.min)} in a combi`;
  };
  const handsText = (c: CardStep) =>
    c.at === "oven" ? `Only to load and unload, ${+((c.perPan ?? 0) * 60).toFixed(0)} s a pan each way` : c.hands >= 1 - 1e-9 ? "All the time" : "Half the time: one cook minds two";
  const WHERE = { table: "Table", range: "Range", oven: "Combi" } as const;
  const holdText = (m: number) => (m >= 120 ? `${+(m / 60).toFixed(1)} h` : `${m} min`);
  /** Each card's rows, and each step's stretches in the plan on the stage: [dish][step] -> [t0, t1][]. */
  let cardRows: HTMLElement[][] = [];
  let cardSeen: (HTMLElement | null)[] = [];
  let stretches: [number, number][][][] = [];
  function renderCards(): void {
    if (!cardsEl) return;
    cardsEl.innerHTML = run.dishes
      .map((d, k) => {
        const card = cardOf(d);
        const rows = card
          .map(
            (c, i) =>
              `<li data-step="${i}"><span class="v">${c.verb}</span><span class="w" data-at="${c.at}">${WHERE[c.at]}</span><span class="m">${perText(d, c)}</span><span class="h" aria-hidden="true"><i style="width:${Math.max(3, Math.round(c.hands * 100))}%"></i></span><em>${handsText(c)}</em></li>`,
          )
          .join("");
        return `<article class="gl-card" style="--dish:${SLOT_COLOR[k]}"><header><b>${d.name}</b><small>${d.batch} to a batch · ${d.pan} to a pan · at its best ${holdText(d.holdMin)}</small></header><ol>${rows}</ol><footer><span data-gl-card-seen>Not seen yet</span></footer></article>`;
      })
      .join("");
    cardRows = Array.from(cardsEl.querySelectorAll<HTMLElement>(".gl-card")).map((c) => Array.from(c.querySelectorAll<HTMLElement>("li")));
    cardSeen = Array.from(cardsEl.querySelectorAll<HTMLElement>("[data-gl-card-seen]"));
    lastCardSig = "";
  }
  /** Each step's stretches in a plan: the table and range jobs with its verb, the combi loads from loading to out. */
  const stretchesOf = (m: Morning): [number, number][][][] =>
    run.dishes.map((d, k) =>
      d.stages.map((st) =>
        st.at === "oven"
          ? m.loads.filter((l) => l.slot === k && l.verb === st.verb).map((l): [number, number] => [l.t0, l.t3])
          : m.jobs.filter((j) => j.slot === k && j.at === st.at && j.verb === st.verb).map((j): [number, number] => [j.t0, j.t1]),
      ),
    );
  let lastCardSig = "";
  function showCards(t: number): void {
    const sig = `${Math.floor(t * 2)}`;
    if (sig === lastCardSig) return;
    lastCardSig = sig;
    stretches.forEach((steps, k) => {
      let seen = 0;
      steps.forEach((xs, i) => {
        const li = cardRows[k]?.[i];
        if (!li) return;
        const was = xs.some(([a]) => a <= t);
        const now = xs.some(([a, b]) => a <= t && b > t);
        if (was) seen++;
        li.toggleAttribute("data-seen", was);
        li.toggleAttribute("data-now", now);
      });
      const el = cardSeen[k];
      if (el) el.textContent = seen === 0 ? "Not seen yet" : seen === steps.length ? `Every step seen: ${steps.length} of ${steps.length}` : `Seen so far: ${seen} of ${steps.length} steps`;
    });
  }

  /* The chart: the cook's hands over each plan, each dish's work stacked as much of the cook as it takes; above, a
   * dish's time in a combi, needing nobody. */
  let playX: ((t: number) => number) | null = null;
  let playY: [number, number] = [0, 0];
  let playLine: SVGLineElement | null = null;
  function renderChart(): void {
    if (!chartEl) return;
    const W = Math.max(280, Math.round(chartEl.clientWidth));
    const narrow = W < 520;
    const L = narrow ? 8 : 132;
    const R = W - (narrow ? 8 : 14);
    const H = narrow ? 24 : 30;
    const OV = 5;
    const max = Math.max(run.endApart, run.end) * 1.02;
    const x = (t: number) => L + (t / max) * (R - L);
    const plans: { key: "together" | "apart"; name: string; m: Morning; end: number }[] = [
      { key: "together", name: "Both at once", m: run.together, end: run.end },
      { key: "apart", name: "One, then the other", m: run.apart, end: run.endApart },
    ];
    const lab = narrow ? 15 : 0;
    const block = lab + 2 * OV + 3 + H + 26;
    const T = 4;
    const total = T + plans.length * block + 20;
    let svg = `<svg viewBox="0 0 ${W} ${total}" width="${W}" height="${total}" class="rs-svg gl-svg">`;
    plans.forEach((p, n) => {
      const y0 = T + n * block;
      const yo = y0 + lab;
      const base = yo + 2 * OV + 3 + H;
      svg += `<text class="nm${p.key === playing ? " on" : ""}" x="${narrow ? L : L - 10}" y="${narrow ? y0 + 10 : base - H / 2 + 4}" text-anchor="${narrow ? "start" : "end"}">${p.name}</text>`;
      svg += `<rect class="lane" x="${L}" y="${(base - H).toFixed(1)}" width="${(R - L).toFixed(1)}" height="${H}" />`;
      // In a combi: from loaded to done, one row per combi; done and waiting for the cook's hands, dashed.
      for (const l of p.m.loads) {
        const yy = yo + l.oven * OV;
        svg += `<rect class="ov" x="${x(l.t1).toFixed(1)}" y="${yy.toFixed(1)}" width="${Math.max(1, x(l.t2) - x(l.t1)).toFixed(1)}" height="${OV - 1}" style="--c:${SLOT_COLOR[l.slot]}"><title>${cap(dishName(run.dishes[l.slot]!))} in a combi, ${span(l.t2 - l.t1)}: nobody needed</title></rect>`;
        if (l.tu > l.t2 + 0.05) svg += `<rect class="ov wait" x="${x(l.t2).toFixed(1)}" y="${yy.toFixed(1)}" width="${Math.max(1, x(l.tu) - x(l.t2)).toFixed(1)}" height="${OV - 1}" style="--c:${SLOT_COLOR[l.slot]}"><title>Done, waiting ${span(l.tu - l.t2)} for the cook</title></rect>`;
      }
      // The cook's hands: each job in the lowest halves of the lane free for its whole stretch.
      const placed: { a: number; b: number; k: number; n: number }[] = [];
      for (const j of [...p.m.jobs].sort((a, b) => a.t0 - b.t0 || b.hands - a.hands)) {
        const nh = Math.max(1, Math.round(j.hands * 2));
        let k = 0;
        for (; k + nh <= 2; k++) if (!placed.some((u) => u.a < j.t1 - 1e-9 && u.b > j.t0 + 1e-9 && u.k < k + nh && u.k + u.n > k)) break;
        if (k + nh > 2) k = 2 - nh;
        placed.push({ a: j.t0, b: j.t1, k, n: nh });
        const yTop = base - ((k + nh) * H) / 2;
        const xa = x(j.t0);
        const w = Math.max(0.8, x(j.t1) - xa - 0.5);
        svg += `<rect class="job" x="${xa.toFixed(1)}" y="${(yTop + 0.5).toFixed(1)}" width="${w.toFixed(1)}" height="${((nh * H) / 2 - 1).toFixed(1)}" style="fill:${SLOT_COLOR[j.slot]}"><title>${j.verb}, ${dishName(run.dishes[j.slot]!)}: ${clockOf(j.t0)} to ${clockOf(j.t1)}${j.hands < 1 ? ", half the cook" : ""}</title></rect>`;
      }
      // When the last of it is done.
      const xe = x(p.end);
      svg += `<line class="end" x1="${xe.toFixed(1)}" x2="${xe.toFixed(1)}" y1="${(yo - 1).toFixed(1)}" y2="${(base + 3).toFixed(1)}" />`;
      svg += `<text class="endl" x="${(xe - 4).toFixed(1)}" y="${(base + 15).toFixed(1)}" text-anchor="end">${span(p.end)}</text>`;
      if (p.key === "together" && run.endApart - run.end > 0.5) {
        const xa = x(run.endApart);
        const yy = base - H / 2;
        svg += `<line class="saved" x1="${(xe + 3).toFixed(1)}" x2="${(xa - 1).toFixed(1)}" y1="${yy.toFixed(1)}" y2="${yy.toFixed(1)}" />`;
        if (xa - xe > 70) svg += `<text class="savedl" x="${((xe + xa) / 2).toFixed(1)}" y="${(yy - 5).toFixed(1)}" text-anchor="middle">${span(run.endApart - run.end)} saved</text>`;
      }
      if (p.key === playing) playY = [yo - 2, base + 2];
    });
    // The clock: minutes from the start.
    const axisY = total - 4;
    const step = max > 150 ? 30 : 15;
    for (let m = 0; m <= max; m += step) {
      if (narrow && m > 0 && m % (2 * step) !== 0) continue;
      svg += `<text class="xa${m === 0 ? " st" : ""}" x="${x(m).toFixed(1)}" y="${axisY}">${m === 0 ? "start" : span(m)}</text>`;
    }
    svg += `<line class="play" x1="0" x2="0" y1="${playY[0]}" y2="${playY[1]}" />`;
    svg += `</svg>`;
    chartEl.innerHTML = svg;
    playX = x;
    playLine = chartEl.querySelector<SVGLineElement>("line.play");
    if (chartNote)
      chartNote.textContent = `A bar's height is how much of the cook the work takes: full height, the whole cook; half, a pan only minded. Above it, outlined, a dish in a combi needing nobody. ${playing === "together" ? "The stage plays both at once." : "The stage plays one, then the other."}`;
  }

  /* Tiles and status. */
  function renderTiles(): void {
    const [a, b] = run.dishes;
    const saved = run.endApart - run.end;
    setTile(tiles.both, span(run.end), `both done at ${clockOf(run.end)}`);
    setTile(tiles.apart, span(run.endApart), `${dishName(a)} ${span(run.endA)}, then ${dishName(b)} ${span(run.endB)}`);
    setTile(tiles.saved, saved > 0.5 ? span(saved) : "None", saved > 0.5 ? `${Math.round((100 * saved) / run.endApart)}% of one after the other` : "both need the cook's hands to cook");
    setTile(tiles.free, run.free > 0.5 ? span(run.free) : "None", run.free > 0.5 ? "of the cook's hands, both at once: room for more work" : "the cook busy from start to finish");
  }
  /** A dish's longest program in a combi, minutes (0: none). */
  const combiMin = (d: MorningDish) => d.stages.reduce((m, s) => (s.at === "oven" ? Math.max(m, s.min) : m), 0);
  function renderStatus(): void {
    if (!statusEl) return;
    const [a, b] = run.dishes;
    const saved = run.endApart - run.end;
    let text: string;
    if (saved <= 0.5) {
      const both = !combiMin(a) && !combiMin(b);
      text = `${both ? `The ${dishName(a)} and the ${dishName(b)} both cook at the range with the cook's hands on them, so neither leaves the cook free while it cooks. ` : ""}Together they take ${span(run.end)}, no less than one after the other.${both ? " Give the cook one dish that cooks in a combi instead." : ""}`;
      statusEl.dataset.state = "warn";
    } else {
      const w = combiMin(a) >= combiMin(b) ? a : b;
      const o = w === a ? b : a;
      text = `The ${dishName(w)} spend${plural(w) ? "" : "s"} ${span(combiMin(w))} in a combi needing nobody, and the ${dishName(o)}${plural(o) ? "'" : "'s"} work fits in while ${plural(w) ? "they" : "it"} cook${plural(w) ? "" : "s"}: together ${span(run.end)}, ${span(saved)} less than one after the other.`;
      if (run.free >= 5) text += ` The cook's hands are still free for ${span(run.free)} of it, room for a third dish.`;
      statusEl.dataset.state = "ok";
    }
    statusEl.textContent = text;
  }

  /* The plan on the stage. */
  /** The plan's end on the stage (a few minutes past its last pan), and the loop: on the plan's clock (minutes), on
   * the stage's (`warp`), and in seconds. */
  let endMin = 0;
  let endStage = 0;
  let warp: Warp = { at: (t) => t, back: (u) => u, walk: () => null, walked: new Map() };
  let loopS = 0;
  let tMin = 0;
  let loopT = 0;
  function stage(): void {
    const m = shown();
    const frac = endMin > 0 ? tMin / endMin : 0;
    const end = playing === "together" ? run.end : run.endApart;
    endMin = end + 4;
    // Laid out once with the walks guessed, then again with them as long as the stage's cook takes, if any differ.
    const walkS = (key: string) => cook.walks.get(key)?.duration ?? WALK_GUESS_S;
    for (let pass = 0; pass < 2; pass++) {
      warp = warpOf(m, endMin, walkS);
      endStage = warp.at(endMin);
      loopS = kitchen.plan(onStage(m, warp.at), 1, endStage, endStage * SEC_PER_MIN + END_HOLD_S);
      if (![...warp.walked].some(([k, g]) => Math.abs(walkS(k) - g) > 0.05)) break;
    }
    tMin = frac * endMin;
    loopT = warp.at(tMin) * SEC_PER_MIN;
    stretches = stretchesOf(m);
    lastCardSig = "";
    lastHud = "";
  }
  function apply(): void {
    run = oneCook(COOK_DISHES[pick[0]]!, COOK_DISHES[pick[1]]!, portions);
    stage();
    renderControls();
    renderCards();
    renderChart();
    renderTiles();
    renderStatus();
  }
  let pending = 0;
  function changed(now: boolean): void {
    window.clearTimeout(pending);
    if (now) apply();
    else {
      renderControls();
      pending = window.setTimeout(apply, 120);
    }
  }
  portionsIn?.addEventListener("input", () => {
    portions = Number(portionsIn.value);
    changed(false);
  });
  for (const b of planBtns)
    b.addEventListener("click", () => {
      const p = b.dataset.glPlan === "apart" ? "apart" : "together";
      if (p === playing) return;
      playing = p;
      tMin = 0;
      endMin = 0;
      stage();
      renderControls();
      renderChart();
    });

  /* The readouts: the clock, what the cook is doing, and what the glasses are seeing. */
  let lastHud = "";
  const doingAt = (m: Morning, t: number) => m.jobs.filter((j) => j.t0 <= t && j.t1 > t);
  /** What the cook is doing at minute `t` of the plan (`walking`: to where, while the plan's clock stands still). */
  function showHud(t: number, walking: string | null): void {
    const m = shown();
    const end = playing === "together" ? run.end : run.endApart;
    const jobs = walking ? [] : doingAt(m, t);
    const hands = jobs.reduce((a, j) => a + j.hands, 0);
    const words = [...new Set(jobs.map((j) => `${j.verb} · ${dishName(run.dishes[j.slot]!)}`))];
    const inOven = m.loads.filter((l) => l.t1 <= t && l.t2 > t);
    const ovenWords = inOven.map((l) => `${cap(dishName(run.dishes[l.slot]!))} in a combi, ${span(l.t2 - t)} to go`);
    let doing = walking ? `Walking to ${walking}` : words.length ? words.join(" + ") : t >= end ? "Both dishes done" : inOven.length ? "Waiting on the combi" : "Waiting";
    if (ovenWords.length && t < end) doing += ` · ${ovenWords.join(" · ")}`;
    let note: string;
    let kind: string;
    if (t >= end) [note, kind] = [`Done at ${clockOf(end)}`, "done"];
    else if (walking) [note, kind] = ["On the move", "walk"];
    else if (hands >= 1 - 1e-9) [note, kind] = ["Hands full", "full"];
    else if (hands > 0) [note, kind] = ["Hands half free", "half"];
    else [note, kind] = ["Hands free", "free"];
    const where = { table: "table", range: "range", oven: "combi" } as const;
    const j0 = jobs[0];
    const seen = walking ? `Walking to ${walking}` : j0 ? `${j0.verb} · ${dishName(run.dishes[j0.slot]!)} · ${where[j0.at]}` : t >= end ? "Done" : inOven.length ? "Combi running · hands free" : "Waiting";
    const sig = `${Math.floor(t)}|${doing}|${note}|${seen}`;
    if (sig === lastHud) return;
    lastHud = sig;
    if (hudClock) hudClock.textContent = clockOf(t);
    if (hudDoing) hudDoing.textContent = doing;
    if (hudNote) {
      hudNote.textContent = note;
      hudNote.hidden = false;
      hudNote.dataset.kind = kind;
    }
    const mm = Math.max(0, Math.floor(t));
    if (tsEl) tsEl.textContent = `GLASSES · ${String(6 + Math.floor(mm / 60)).padStart(2, "0")}:${String(mm % 60).padStart(2, "0")}`;
    if (seenEl) seenEl.textContent = seen;
  }
  /** The label on the stage that points at the glasses. */
  const tagAt = new THREE.Vector3();
  function showTag(): void {
    if (!tagEl) return;
    glasses.group.getWorldPosition(tagAt);
    tagAt.y += 0.1;
    tagAt.project(camera);
    const on = cook.st.group.visible && tagAt.z < 1 && Math.abs(tagAt.x) < 0.95 && Math.abs(tagAt.y) < 0.95;
    tagEl.hidden = !on;
    if (on) tagEl.style.transform = `translate(${(((tagAt.x + 1) / 2) * stageEl!.clientWidth).toFixed(1)}px, ${(((1 - tagAt.y) / 2) * stageEl!.clientHeight).toFixed(1)}px)`;
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
  for (const k of kitchen.cooks) void k.st.ready.then(() => void gate.refresh(scene));

  let inView = false;
  observeInView(stageEl, (v) => {
    inView = v;
  });
  let povInView = false;
  observeInView(povEl, (v) => {
    povInView = v;
  });
  let frozen: number | null = null;
  let last = performance.now();
  let povT = 0;
  const fitClip = clipFitter(camera, scene, 0.05);
  function animate(now: number): void {
    requestAnimationFrame(animate);
    const dt = Math.max(0, Math.min(0.05, (now - last) / 1000));
    last = now;
    if ((!inView && !povInView) || document.hidden || !gate.ready()) return;
    if (frozen === null) loopT = (loopT + dt) % loopS;
    else loopT = warp.at(frozen) * SEC_PER_MIN;
    const u = Math.min(endStage, loopT / SEC_PER_MIN);
    tMin = Math.min(endMin, warp.back(u));
    const walking = warp.walk(u);
    kitchen.pose(loopT);
    kitchen.show(u, dt);
    showHud(tMin, walking);
    // A step lights up as the cook gets to it.
    showCards(walking ? tMin - 1e-4 : tMin);
    if (playX && playLine) {
      const xx = playX(Math.min(tMin, playing === "together" ? run.end : run.endApart)).toFixed(1);
      playLine.setAttribute("x1", xx);
      playLine.setAttribute("x2", xx);
    }
    glasses.led.color.setHex(Math.floor(now / 500) % 2 === 0 ? 0xf5f7ff : 0x3a3d48);
    if (inView) {
      controls.update();
      fitClip();
      renderer.render(scene, camera);
      showTag();
    }
    povT += dt;
    if (povInView && (povT >= POV_S || frozen !== null) && povGate.ready()) {
      povT = 0;
      placePov();
      glasses.group.visible = false;
      rPov!.render(scene, pov);
      glasses.group.visible = true;
      if (!povCanvas!.dataset.drawn) povCanvas!.dataset.drawn = "1";
    }
  }
  requestAnimationFrame(animate);
  void gate.refresh(scene).then(() => void povGate.refresh(scene));

  // A hook for checking the page by script.
  Object.assign(stageEl, {
    glasses: {
      info: () => ({
        mountMs: Math.round(mountMs),
        dishes: run.dishes.map((d) => d.key),
        portions,
        playing,
        end: +run.end.toFixed(1),
        endA: +run.endA.toFixed(1),
        endB: +run.endB.toFixed(1),
        endApart: +run.endApart.toFixed(1),
        saved: +(run.endApart - run.end).toFixed(1),
        free: +run.free.toFixed(1),
        cooks: kitchen.cooks.filter((k) => k.st.group.visible).length,
        seen: cardRows.map((rows) => rows.filter((r) => r.hasAttribute("data-seen")).length),
        steps: cardRows.map((rows) => rows.length),
        tag: tagEl ? !tagEl.hidden : null,
        face: { measured: face.measured, eye: face.eye.toArray().map((v) => +v.toFixed(4)), half: +face.half.toFixed(4) },
        glassesAt: glasses.group.getWorldPosition(new THREE.Vector3()).toArray().map((v) => +v.toFixed(3)),
        headAt: head.getWorldPosition(new THREE.Vector3()).toArray().map((v) => +v.toFixed(3)),
        loopS: +loopS.toFixed(2),
        t: +tMin.toFixed(1),
      }),
      /** Where the cook's head is at each of some minutes of the plan: [minute, the stage's minute, x, z, walking to]. */
      track: (ts: number[]) =>
        ts.map((t) => {
          const u = warp.at(t);
          kitchen.pose(u * SEC_PER_MIN);
          const r = head.getWorldPosition(new THREE.Vector3());
          return [t, +u.toFixed(2), +r.x.toFixed(2), +r.z.toFixed(2), warp.walk(u)];
        }),
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
      set(s: { a?: string; b?: string; portions?: number; plan?: "together" | "apart" }) {
        const ix = (key: string) => COOK_DISHES.findIndex((d) => d.key === key);
        const was: [number, number] = [pick[0], pick[1]];
        if (s.a !== undefined && ix(s.a) >= 0) pick[0] = ix(s.a);
        if (s.b !== undefined && ix(s.b) >= 0) pick[1] = ix(s.b);
        if (pick[0] === pick[1]) [pick[0], pick[1]] = was;
        if (s.portions !== undefined) portions = s.portions;
        if (s.plan !== undefined && s.plan !== playing) {
          playing = s.plan;
          tMin = 0;
          endMin = 0;
        }
        changed(true);
      },
    },
  });
  const mountMs = performance.now() - mountT0;
}
