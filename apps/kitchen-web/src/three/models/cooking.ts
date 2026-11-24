/**
 * The rest of the hot line: griddles, the charbroiler, the fryer, the salamander, induction, the tilting skillet, the
 * steam kettle and the strip warmer. Dimensions are the archetypes' (assets/kitchens/specs); the construction is the
 * equipment class's: welded and folded 304 sheet on adjustable legs, gas valves behind knob panels, cast iron where it
 * takes the heat, and every part that moves in life a joint that moves here.
 *
 * Live names, `Hooks.set(name, v)`:
 *   heat         0..1  glow of whatever radiates: charbroiler radiants, salamander ceramics, strip-warmer calrod,
 *                      kettle and skillet "heating" lamps.
 *   lift_l/lift_r 0..1 fryer baskets: 0 down in the oil, 1 hung on the hanger rod to drain.
 *   rack         0..1  salamander broiling grid height, lowest to highest notch.
 *   zone0..zone3 0..1  induction zones' power: lights each zone's indicator and its bar of LEDs.
 *   tilt         0..1  tilting skillet pan, level to fully tipped (about 85°); the handwheel turns with it.
 *   lid          0..1  skillet and kettle lids, shut to open.
 *   level        0..1  kettle contents, empty to the brim line.
 * `open(f)` works the main door, drawer or lid: griddle grease drawer, charbroiler grease pan, fryer cabinet door,
 * salamander grid (rolls out), induction-suite doors, skillet and kettle lids. Gas burners follow `flames.setFlow`, and
 * their knobs turn with their valves.
 */
import * as THREE from "three";
import { mergeGeometries } from "three/examples/jsm/utils/BufferGeometryUtils.js";
import { FlameSet, linePorts, type Port } from "./flames";
import { valveAngle } from "./hotline";
import { register, type Builder, type Ext, type Hooks } from "./index";
import { glowMaterial, palette, setGlow } from "./materials";
import {
  ball,
  box,
  circleHole,
  cyl,
  knob,
  lathe,
  plate,
  rbox,
  rectPath,
  screw,
  slideJoint,
  torus,
  turnJoint,
  type Built,
  type Shop,
  type V3,
} from "./shop";

/* ------------------------------------------------------------------------------------------------------------ */
/* Parts                                                                                                          */
/* ------------------------------------------------------------------------------------------------------------ */

/** A sheet-metal part: its edge radius drawn close up, a plain box (12 triangles) from further off. */
export function sheet(s: Shop, w: number, h: number, d: number, r = 0.0008): THREE.BufferGeometry {
  return s.fine(2) ? rbox(s, w, h, d, r) : box(w, h, d);
}

/**
 * A round bar bent along a polyline, as `bar` makes it but drawn economically: each straight run one cylinder, each
 * bend (of radius `bend`) a few rings, the ends capped. Wire frames, rods, pipes and handles are mostly straight.
 */
export function pipe(s: Shop, pts: readonly V3[], r: number, bend = r * 3, radial = 12): THREE.BufferGeometry {
  const v = pts.map((q) => new THREE.Vector3(q[0], q[1], q[2]));
  const rad = s.seg(radial, 5);
  const parts: THREE.BufferGeometry[] = [];
  const up = new THREE.Vector3(0, 1, 0);
  const straight = (a: THREE.Vector3, b: THREE.Vector3) => {
    const len = a.distanceTo(b);
    if (len < 1e-5) return;
    const g = new THREE.CylinderGeometry(r, r, len, rad, 1, true);
    g.applyQuaternion(new THREE.Quaternion().setFromUnitVectors(up, b.clone().sub(a).normalize()));
    g.translate((a.x + b.x) / 2, (a.y + b.y) / 2, (a.z + b.z) / 2);
    parts.push(g);
  };
  let from = v[0]!.clone();
  for (let i = 1; i < v.length - 1; i++) {
    const [p0, p1, p2] = [v[i - 1]!, v[i]!, v[i + 1]!];
    const cut = Math.min(bend, p0.distanceTo(p1) * 0.49, p1.distanceTo(p2) * 0.49);
    const a = p1.clone().addScaledVector(p1.clone().sub(p0).normalize(), -cut);
    const b = p1.clone().addScaledVector(p2.clone().sub(p1).normalize(), cut);
    straight(from, a);
    if (cut > 1e-5) parts.push(new THREE.TubeGeometry(new THREE.QuadraticBezierCurve3(a, p1.clone(), b), s.seg(6, 2), r, rad, false));
    from = b;
  }
  straight(from, v[v.length - 1]!);
  if (s.fine(1))
    for (const [e0, e1] of [[v[0]!, v[1]!], [v[v.length - 1]!, v[v.length - 2]!]] as const) {
      const c = new THREE.CircleGeometry(r, rad);
      c.applyQuaternion(new THREE.Quaternion().setFromUnitVectors(new THREE.Vector3(0, 0, 1), e0.clone().sub(e1).normalize()));
      c.translate(e0.x, e0.y, e0.z);
      parts.push(c);
    }
  for (const g of parts) for (const k of Object.keys(g.attributes)) if (k !== "position" && k !== "normal" && k !== "uv") g.deleteAttribute(k);
  return mergeGeometries(parts, false) ?? parts[0]!;
}

/** An adjustable stainless leg `h` tall on its flanged foot, with the lock nut on its thread. */
export function flangedLeg(s: Shop, x: number, z: number, h: number, r = 0.019): void {
  const p = palette();
  const fh = Math.min(0.03, h * 0.3);
  s.add(p.ss, cyl(s, r, h - fh, "y", 20), [x, fh + (h - fh) / 2, z]);
  s.add(p.ssBright, lathe(s, [
    [0, 0],
    [r * 1.5, 0],
    [r * 1.58, 0.005],
    [r * 1.18, fh * 0.45],
    [r * 0.8, fh],
    [0, fh],
  ], 24), [x, 0, z]);
  if (s.fine(2)) s.add(p.ssBright, new THREE.CylinderGeometry(r * 1.2, r * 1.2, 0.008, 6), [x, fh + 0.004, z]);
}

/** The rating plate every listed appliance carries: a riveted foil plate with ruled fields and no maker's name. */
export function dataPlate(s: Shop, at: V3, w = 0.08, h = 0.05, rot: V3 = [0, 0, 0]): void {
  if (!s.fine(1)) return;
  const p = palette();
  s.at(at, rot, () => {
    s.add(p.ssBright, box(w, h, 0.0006), [0, 0, 0.0003]);
    if (!s.fine(2)) return;
    for (let i = 0; i < 5; i++) s.add(p.charcoal, box(w * (i % 2 ? 0.46 : 0.7), 0.0007, 0.0003), [-(w * 0.1) * (i % 2) - w * 0.06, h * 0.26 - i * h * 0.13, 0.00075]);
    for (const sy of [-1, 1]) s.add(p.charcoal, box(w - 0.014, 0.0006, 0.0003), [0, sy * (h / 2 - 0.0035), 0.00075]);
    for (const [sx, sy] of [[-1, -1], [1, -1], [-1, 1], [1, 1]] as const) s.add(p.ssBright, ball(s, 0.0016, 8, Math.PI / 2), [sx * (w / 2 - 0.004), sy * (h / 2 - 0.004), 0.0006], [Math.PI / 2, 0, 0]);
  });
}

/** A pressed louvre panel on a face looking +z: `n` lanced slats proud of the sheet over dark slots, `w` by `h`. */
export function louvres(s: Shop, at: V3, w: number, h: number, n: number, mat?: THREE.Material, rot: V3 = [0, 0, 0]): void {
  const p = palette();
  const m = mat ?? p.ss;
  s.at(at, rot, () => {
    s.add(p.enamel, box(w, h, 0.0005), [0, 0, 0.00065]);
    if (!s.fine(1)) return;
    const pitch = h / n;
    for (let i = 0; i < n; i++) s.add(m, sheet(s, w + 0.004, pitch * 0.62, 0.0009, 0.0003), [0, -h / 2 + pitch * (i + 0.5), 0.0024], [-0.75, 0, 0]);
  });
}

/** A valve or thermostat knob on its own joint (so it turns with its valve), facing +z at `at`. */
function valveKnob(s: Shop, name: string, at: V3, r = 0.021): string {
  const p = palette();
  if (!s.fine(1)) {
    // From across the room a knob is a dark stub merged into the body (it no longer turns).
    s.add(p.bakelite, cyl(s, r * 0.8, 0.026, "z", 10), [at[0], at[1], at[2] + 0.013]);
    return name;
  }
  const j = s.joint(name, at);
  if (s.fine(2)) knob(j, p.bakelite, p.ssBright, [0, 0, 0], r);
  else if (s.fine(1)) {
    // From the aisle: the same skirted body and pointer, in fewer facets.
    j.add(p.ssBright, cyl(j, r * 1.18, 0.002, "z", 16), [0, 0, 0.001]);
    j.add(p.bakelite, lathe(j, [[0, 0.002], [r, 0.002], [r, 0.006], [r * 0.72, 0.012], [r * 0.7, 0.024], [0, 0.029]], 16).rotateX(Math.PI / 2));
    j.add(p.bakelite, box(r * 0.28, r * 1.3, 0.012), [0, r * 0.2, 0.022]);
  }
  return name;
}

/** The printed dial round a knob: tick marks from off (up) counter-clockwise to low (down). */
function dialTicks(s: Shop, at: V3, r: number, n = 9): void {
  if (!s.fine(2)) return;
  const p = palette();
  for (let k = 0; k < n; k++) {
    const a = Math.PI / 2 + (k / (n - 1)) * Math.PI;
    const big = k === 0 || k === n - 1 || k === (n - 1) / 2;
    s.add(p.enamel, box(0.0013, big ? 0.006 : 0.0035, 0.0004), [at[0] + Math.cos(a) * r, at[1] + Math.sin(a) * r, at[2] + 0.0003], [0, 0, a - Math.PI / 2]);
  }
}

/** Ports along a burner pipe running front to back at (x, y), from z0 to z1: `linePorts` turned through 90°. */
export function zLinePorts(burner: number, y: number, x: number, z0: number, z1: number, n: number, tilt = 0.2): Port[] {
  return linePorts(burner, y, x, z0, z1, n, tilt).map((q) => ({ at: [q.at[2], q.at[1], q.at[0]] as V3, dir: [q.dir[2], q.dir[1], q.dir[0]] as V3, burner: q.burner }));
}

/** How many of a burner's `n` ports carry a flame at this detail: each flame is a few hundred triangles. */
function portCount(s: Shop, n: number): number {
  return s.fine(2) ? n : s.fine(1) ? Math.ceil(n * 0.45) : Math.ceil(n * 0.25);
}

/** Gas hooks: the flame set added to the model, pilots (numbered after the valves) lit, and knobs following valves. */
function gasHooks(b: Built, ports: Port[], valves: number, size: number, knobs: string[]): Hooks & { flames: FlameSet } {
  const flames = new FlameSet(ports, size);
  b.group.add(flames.mesh);
  flames.setAll(0);
  if (flames.burners > valves) flames.setFlow(valves, 0.12);
  const set = flames.setFlow.bind(flames);
  flames.setFlow = (i: number, f: number) => {
    set(i, f);
    if (i < valves) turnJoint(b.joints.get(knobs[i]!), "z", valveAngle(f));
  };
  return { flames, burners: valves, tick: (t) => flames.tick(t) };
}

/** A small emissive part whose brightness is live (an indicator lamp, a display): one per model, as glows are. */
function lampMaterial(color: number): THREE.MeshStandardMaterial {
  const m = new THREE.MeshStandardMaterial({ color: 0x1a0d0b, metalness: 0, roughness: 0.35, emissive: color, emissiveIntensity: 0 });
  m.userData.uv = "own";
  return m;
}

/**
 * Bubbles breaking on a liquid surface, all in one instanced draw: each grows, sits a moment and bursts, then rises
 * again somewhere else. Lay it out in the frame it is added to: `at` the surface centre, `w` by `d` its extent.
 */
export class Bubbles {
  readonly mesh: THREE.InstancedMesh;
  private life: Float32Array;
  private m = new THREE.Matrix4();
  private q = new THREE.Quaternion();
  constructor(
    mat: THREE.Material,
    n: number,
    public at: V3,
    public w: number,
    public d: number,
    private r = 0.004,
    round = false,
  ) {
    const g = new THREE.SphereGeometry(1, 8, 4, 0, Math.PI * 2, 0, Math.PI / 2);
    this.mesh = new THREE.InstancedMesh(g, mat, n);
    this.mesh.frustumCulled = false;
    this.mesh.castShadow = false;
    this.mesh.count = 0;
    this.life = new Float32Array(n);
    for (let i = 0; i < n; i++) this.life[i] = ((i * 0.618034) % 1) * (0.5 + (i % 3) * 0.2);
    this.round = round;
  }
  private round: boolean;
  /** `rate` 0 is still, 1 a hard boil. */
  tick(t: number, rate: number): void {
    const n = this.life.length;
    const k = Math.round(n * Math.min(1, Math.max(0, rate)));
    this.mesh.count = k;
    if (!k) return;
    const sc = new THREE.Vector3();
    const pos = new THREE.Vector3();
    for (let i = 0; i < k; i++) {
      const period = 0.5 + ((i * 0.37) % 1) * 0.7;
      const ph = (t / period + this.life[i]!) % 1;
      const cycle = Math.floor(t / period + this.life[i]!);
      const h1 = Math.sin((i + 1) * 12.9898 + cycle * 78.233) * 43758.5453;
      const h2 = Math.sin((i + 1) * 39.3468 + cycle * 11.135) * 24634.6345;
      let u = h1 - Math.floor(h1) - 0.5;
      let v = h2 - Math.floor(h2) - 0.5;
      if (this.round) {
        const a = u * Math.PI * 2;
        const rr = Math.sqrt(v + 0.5) * 0.5;
        u = Math.cos(a) * rr;
        v = Math.sin(a) * rr;
      }
      const size = this.r * (0.5 + ((i * 0.7548) % 1)) * Math.sin(Math.PI * Math.min(1, ph * 1.15));
      pos.set(this.at[0] + u * this.w, this.at[1], this.at[2] + v * this.d);
      sc.set(Math.max(1e-5, size), Math.max(1e-5, size * 0.7), Math.max(1e-5, size));
      this.m.compose(pos, this.q, sc);
      this.mesh.setMatrixAt(i, this.m);
    }
    this.mesh.instanceMatrix.needsUpdate = true;
  }
}

/** A grid of round wires `w` by `h` in its own xy plane, `pitch` apart: fry-basket mesh, broiler grids, screens. */
export function wireGrid(s: Shop, mat: THREE.Material, w: number, h: number, pitch: number, r: number, at: V3, rot: V3 = [0, 0, 0]): void {
  const nx = Math.max(1, Math.round(w / pitch));
  const ny = Math.max(1, Math.round(h / pitch));
  const wire = (len: number, axis: "x" | "y") =>
    s.fine(2) ? new THREE.CylinderGeometry(r, r, len, 5, 1, true).rotateZ(axis === "x" ? Math.PI / 2 : 0) : axis === "x" ? box(len, r * 2, r * 2) : box(r * 2, len, r * 2);
  s.at(at, rot, () => {
    for (let i = 1; i < nx; i++) s.add(mat, wire(h, "y"), [-w / 2 + (i * w) / nx, 0, r]);
    for (let j = 1; j < ny; j++) s.add(mat, wire(w, "x"), [0, -h / 2 + (j * h) / ny, -r * 0.2]);
  });
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Countertop gas appliances: griddles and the charbroiler share a chassis.                                       */
/* ------------------------------------------------------------------------------------------------------------ */

interface Shell {
  /** Leg height to the body's underside (4 in). */
  legH: number;
  /** The control panel's face. */
  zf: number;
  /** The body sides' top. */
  sideTop: number;
  /** The splashes' foot, and the side splashes' inside face. */
  splashBase: number;
  splashX: number;
  /** The cooking surface's front edge (the trough's back wall) and back edge (the back splash). */
  zTop0: number;
  zTop1: number;
  /** The splashes' top. */
  yS: number;
  /** The grease trough's floor. */
  troughY: number;
  /** Openings cut in the control panel, traced in its (x, y) the hole way round. */
  holes: THREE.Path[];
}

/** The bull-nose front rail's centre line. */
const NOSE: [number, number, number] = [0.289, -0.031, 0.013];

/**
 * A countertop chassis: legs, body, back and side splashes with the flue outlet behind, the grease trough along the
 * front of the cooking surface, the bull-nose rail and the control panel under it.
 */
function counterShell(s: Shop, e: Ext, o: Shell): void {
  const p = palette();
  const { w, d } = e;
  const back = -d;
  const { legH, zf, yS } = o;
  const bodyD = zf - back;
  // 4 in adjustable legs, each on a welded pad.
  for (const lx of [0.045, w - 0.045])
    for (const lz of [zf - 0.05, back + 0.055]) {
      flangedLeg(s, lx, lz, legH + 0.002, 0.016);
      s.add(p.ssDull, sheet(s, 0.05, 0.005, 0.05, 0.001), [lx, legH - 0.0025, lz]);
    }
  // Body: bottom, sides, back. The sides rise in front of the cooking surface to close the trough's ends.
  s.add(p.ssDull, sheet(s, w - 0.0032, 0.0016, bodyD), [w / 2, legH + 0.0008, (zf + back) / 2]);
  const yN = NOSE[0] + NOSE[2] - 0.002;
  for (const sx of [0.0008, w - 0.0008]) {
    s.add(p.ss, sheet(s, 0.0016, o.sideTop - legH, bodyD, 0.0006), [sx, (legH + o.sideTop) / 2, (zf + back) / 2]);
    s.add(p.ss, sheet(s, 0.0016, yN - o.sideTop + 0.001, o.zTop0 - zf + 0.004, 0.0006), [sx, (o.sideTop + yN) / 2 - 0.0005, (o.zTop0 + zf) / 2 - 0.002]);
    // Behind the back splash: the flue chamber's sides.
    s.add(p.ss, sheet(s, 0.0016, yS - o.sideTop + 0.001, o.zTop1 - back, 0.0006), [sx, (yS + o.sideTop) / 2, (o.zTop1 + back) / 2]);
  }
  s.add(p.ssDull, sheet(s, w - 0.0032, yS - legH, 0.0016), [w / 2, (yS + legH) / 2, back + 0.0008]);
  // The back splash, its top folded back into the flue cap, slotted where the flue gases leave.
  s.add(p.ss, sheet(s, w - 0.0032, yS - o.splashBase, 0.0016), [w / 2, (yS + o.splashBase) / 2, o.zTop1 - 0.0008]);
  s.add(p.ss, sheet(s, w, 0.0016, o.zTop1 - back + 0.0016), [w / 2, yS + 0.0008, (o.zTop1 + back) / 2 - 0.0008]);
  if (s.fine(1)) {
    const slots = Math.max(2, Math.round((w - 0.1) / 0.11));
    const sw = (w - 0.1) / slots;
    for (let i = 0; i < slots; i++)
      for (const zz of [0.35, 0.65]) s.add(p.enamel, box(sw - 0.014, 0.0004, 0.01), [0.05 + sw * (i + 0.5), yS + 0.0018, o.zTop1 + (back - o.zTop1) * zz]);
  }
  // Side splashes: tapering from the back splash's height to a low front, the top edge rolled.
  const hb = yS - o.splashBase;
  const hf = 0.034;
  const zFront = o.zTop0 - 0.01;
  const L = zFront - o.zTop1;
  const sh = new THREE.Shape();
  sh.moveTo(0, 0);
  sh.lineTo(L, 0);
  sh.lineTo(L, hb);
  sh.lineTo(0.02, hf);
  sh.quadraticCurveTo(0, hf, 0, hf - 0.02);
  sh.lineTo(0, 0);
  for (const [sx, flip] of [[o.splashX, 0], [w - o.splashX - 0.0016, 1]] as const) {
    const g = plate(s, sh, 0.0016, 0, 6);
    g.rotateY(Math.PI / 2);
    s.add(p.ss, g, [sx + (flip ? 0 : 0), o.splashBase, zFront]);
    if (s.fine(2)) {
      const x = sx + 0.0008;
      s.add(p.ss, pipe(s, [[x, o.splashBase + hf - 0.004, zFront - 0.006], [x, o.splashBase + hf, zFront - 0.02], [x, yS, o.zTop1]], 0.0019, 0.01, 8));
    }
  }
  // The grease trough: a formed channel between the cooking surface's front edge and the bull-nose.
  const zt1 = NOSE[1] - 0.011;
  s.add(p.ss, sheet(s, w - 0.0032, 0.0012, zt1 - o.zTop0 + 0.002), [w / 2, o.troughY, (o.zTop0 + zt1) / 2]);
  s.add(p.ss, sheet(s, w - 0.0032, o.sideTop - o.troughY + 0.002, 0.0012), [w / 2, (o.sideTop + o.troughY) / 2, o.zTop0 - 0.0008]);
  s.add(p.ss, sheet(s, w - 0.0032, NOSE[0] - o.troughY, 0.0012), [w / 2, (NOSE[0] + o.troughY) / 2, zt1]);
  // The bull-nose rail across the front.
  s.add(p.ss, cyl(s, NOSE[2], w - 0.0032, "x", 28), [w / 2, NOSE[0], NOSE[1]]);
  // The control panel, cut for the knob stems' bezels, the sight holes and the drawer.
  const py1 = NOSE[0] - 0.008;
  if (s.fine(1) && o.holes.length) {
    const sp = new THREE.Shape();
    rectPath(sp, w - 0.0032, py1 - legH, 0, w / 2, (py1 + legH) / 2);
    sp.holes.push(...o.holes);
    s.add(p.ss, plate(s, sp, 0.0012, 0, 16), [0, 0, zf - 0.0012]);
  } else s.add(p.ss, box(w - 0.0032, py1 - legH, 0.0012), [w / 2, (py1 + legH) / 2, zf - 0.0006]);
  if (s.fine(2)) for (const sx of [0.018, w - 0.018]) for (const sy of [legH + 0.018, py1 - 0.016]) screw(s, p.ssBright, [sx, sy, zf], 0.0026);
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Griddles                                                                                                       */
/* ------------------------------------------------------------------------------------------------------------ */

/**
 * A countertop gas griddle: a 1 in steel plate over one U-shaped tube burner per 12 in section, each on its own
 * thermostat. Grease is scraped into the front trough and drains through a hole at the left end into a drawer below.
 * Each section has a sight hole beside its knob for the standing pilot.
 */
function griddle(): Builder {
  return (s, e) => {
    const p = palette();
    const { w } = e;
    const legH = 0.102;
    const top = 0.305; // the cooking surface (the archetype's platen)
    const bot = top - 0.0254;
    const zf = -0.032;
    const zp0 = -0.084;
    const zp1 = -0.7;
    const n = Math.max(1, Math.round(w / 0.3048));
    const sec = w / n;
    const R = sec * 0.22; // half the U's span
    const by = 0.226; // burner tube axis
    const tr = 0.0125;
    const ky = 0.228;
    const dw = Math.min(0.3, w * 0.4);
    const dx0 = 0.028;
    const dy0 = 0.114;
    const dh = 0.05;
    const holes: THREE.Path[] = [];
    holes.push(rectPath(new THREE.Path(), dw - 0.01, dh - 0.008, 0.002, dx0 + dw / 2, dy0 + dh / 2, true));
    const sight = (i: number) => sec * (i + 0.5) + R;
    for (let i = 0; i < n; i++) holes.push(circleHole(sight(i), ky, 0.0075));
    counterShell(s, e, { legH, zf, sideTop: top - 0.0008, splashBase: top - 0.001, splashX: 0.0022, zTop0: zp0, zTop1: zp1, yS: 0.3775, troughY: 0.2585, holes });

    // The plate: 1 in steel, seasoned black, its edges eased.
    s.add(p.griddle, rbox(s, w - 0.0036, top - bot, zp0 - zp1, 0.0016), [w / 2, (top + bot) / 2, (zp0 + zp1) / 2]);
    // The trough drains at the left end, down a chute into the drawer.
    const drainX = dx0 + 0.04;
    const drainZ = (zp0 + NOSE[1] - 0.011) / 2;
    s.add(p.enamel, cyl(s, 0.009, 0.0008, "y", 20), [drainX, 0.2585 + 0.001, drainZ]);
    if (s.fine(1)) {
      s.add(p.ss, torus(s, 0.0095, 0.0012, "y", Math.PI * 2, 20, 6), [drainX, 0.2595, drainZ]);
      s.add(p.ssDull, cyl(s, 0.0095, 0.085, "y", 14, 0.0095, true), [drainX, 0.2585 - 0.0425, drainZ]);
    }

    // Knobs: one thermostat per section, with its dial; a chrome grommet round each sight hole.
    const knobs: string[] = [];
    for (let i = 0; i < n; i++) {
      const kx = sec * (i + 0.5) - R;
      knobs.push(valveKnob(s, `knob_${i}`, [kx, ky, zf]));
      dialTicks(s, [kx, ky, zf], 0.031);
      if (s.fine(1)) s.add(p.ssBright, torus(s, 0.0088, 0.0015, "z", Math.PI * 2, 20, 6), [sight(i), ky, zf + 0.0006]);
    }
    dataPlate(s, [w - 0.075, dy0 + 0.028, zf], 0.075, 0.042);

    // Under the plate: per section a U-tube burner, its venturi at the front of the left leg on the valve's orifice,
    // the right leg closed and its front end beside the standing pilot behind the sight hole.
    const ports: Port[] = [];
    const zA = -0.118;
    const zB = -0.6;
    for (let i = 0; i < n; i++) {
      const cx = sec * (i + 0.5);
      if (s.fine(1)) {
        s.add(p.ssHeat, cyl(s, tr, zA - zB, "z", 14), [cx - R, by, (zA + zB) / 2]);
        s.add(p.ssHeat, cyl(s, tr, zA + 0.02 - zB, "z", 14), [cx + R, by, (zA + 0.02 + zB) / 2]);
        s.add(p.ssHeat, torus(s, R, tr, "y", Math.PI, 24, 10), [cx - 0, by, zB]);
        // Venturi bell and air shutter, the valve's orifice hood in its mouth.
        s.add(p.ssHeat, cyl(s, tr, 0.04, "z", 14, 0.019), [cx - R, by, zA + 0.02]);
        s.add(p.ssDull, cyl(s, 0.021, 0.0015, "z", 16), [cx - R, by, zA + 0.041]);
        s.add(p.brass, rbox(s, 0.026, 0.028, 0.026, 0.002), [cx - R, ky, zf - 0.016]);
        s.add(p.brass, cyl(s, 0.004, 0.022, "z", 8), [cx - R, by, zA + 0.05]);
        // The pilot: a capillary up to a hooded tip just behind the sight hole.
        if (s.fine(2)) {
          s.add(p.ssBright, pipe(s, [[cx - R + 0.012, 0.16, zf - 0.03], [cx + R, 0.16, zf - 0.03], [cx + R, by - 0.004, zf - 0.03]], 0.0018, 0.008, 6));
          s.add(p.ssBright, cyl(s, 0.0035, 0.012, "y", 10), [cx + R, by, zf - 0.03]);
        }
      }
      ports.push(...zLinePorts(i, by + tr, cx - R, zA + 0.01, zB, portCount(s, 14), 0.25), ...zLinePorts(i, by + tr, cx + R, zA + 0.03, zB, portCount(s, 13), 0.25));
      const nb = portCount(s, 6);
      for (let k = 1; k < nb; k++) {
        const a = (k / nb) * Math.PI;
        ports.push({ at: [cx + R * Math.cos(a), by + tr, zB - R * Math.sin(a)], dir: [0, 1, 0], burner: i });
      }
    }
    // Pilots, numbered after the valves: a small flame each, licking back toward its burner.
    for (let i = 0; i < n; i++)
      for (const dx of [-0.002, 0.002]) ports.push({ at: [sight(i) + dx, by + 0.006, zf - 0.03], dir: [dx * 40, 0.5, -1], burner: n });
    if (s.fine(2)) s.add(p.enamel, cyl(s, 0.0105, w - 0.1, "x", 12), [w / 2, 0.165, zf - 0.03]);
    // A heat shield under the burners, dark from the heat.
    if (s.fine(1)) s.add(p.ssHeat, box(w - 0.01, 0.0012, zA - zp1 + 0.03), [w / 2, 0.186, (zA + zp1) / 2 - 0.015]);

    // The grease drawer: a stainless pan on slides, its front lapping the panel's opening.
    const dr = s.joint("grease_drawer", [dx0, dy0, zf]);
    dr.add(p.ss, sheet(dr, dw, dh + 0.002, 0.0016, 0.0008), [dw / 2, dh / 2, 0.0008]);
    pullHandle0(dr, [dw / 2, dh / 2, 0.0016], Math.min(0.1, dw * 0.5));
    if (s.fine(1)) {
      const dd = 0.3;
      dr.add(p.ssDull, box(dw - 0.016, 0.001, dd), [dw / 2, 0.002, -dd / 2]);
      for (const sx of [0.008, dw - 0.008]) dr.add(p.ssDull, box(0.001, dh - 0.01, dd), [sx, dh / 2 - 0.003, -dd / 2]);
      dr.add(p.ssDull, box(dw - 0.016, dh - 0.01, 0.001), [dw / 2, dh / 2 - 0.003, -dd]);
      // Last service's grease, gone hard and dark.
      dr.add(p.griddle, box(dw - 0.02, 0.004, dd - 0.01), [dw / 2, 0.0045, -dd / 2]);
    }

    return (b) => {
      const g = gasHooks(b, ports, n, 0.02, knobs);
      const rest: V3 = [dx0, dy0, zf];
      g.open = (f) => slideJoint(b.joints.get("grease_drawer"), "z", f * 0.24, rest);
      return g;
    };
  };
}

/** A small bar pull on two posts, along x, for drawers and pans: stands 16 mm proud. */
function pullHandle0(s: Shop, at: V3, len: number): void {
  const p = palette();
  const h = len / 2;
  s.at(at, [0, 0, 0], () => {
    s.add(p.ssBright, pipe(s, [[-h, 0, 0], [-h, 0, 0.016], [h, 0, 0.016], [h, 0, 0]], 0.0055, 0.008, 10));
    if (s.fine(2)) for (const x of [-h, h]) s.add(p.ssBright, cyl(s, 0.008, 0.003, "z", 14), [x, 0, 0.0015]);
  });
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Charbroiler                                                                                                    */
/* ------------------------------------------------------------------------------------------------------------ */

/**
 * A radiant charbroiler: one straight pipe burner per 6 in section, each under a stainless radiant (an inverted V
 * slotted along its ridge, so the flame tips lick through) that glows and vaporises the drippings, under cast-iron
 * grates. The grates' bars carry a grease channel along their tops and slope down to the front trough; what falls
 * through lands in a water pan that slides out under the knobs.
 */
const charbroiler: Builder = (s, e) => {
  const p = palette();
  const { w } = e;
  const legH = 0.102;
  const zf = -0.032;
  const zg0 = -0.072;
  const zg1 = -0.698;
  const yf = 0.318; // the grates' top at the front and at the back: they fall 18 mm toward the trough
  const yb = 0.336;
  const n = Math.max(1, Math.round(w / 0.152));
  const sec = (w - 0.004) / n;
  const by = 0.254;
  const br = 0.013;
  const ky = 0.236;
  const panY0 = 0.112;
  const panH = 0.05;
  const holes = [rectPath(new THREE.Path(), w - 0.07, panH - 0.008, 0.002, w / 2, panY0 + panH / 2, true)];
  counterShell(s, e, { legH, zf, sideTop: 0.305, splashBase: 0.304, splashX: 0, zTop0: -0.078, zTop1: zg1 - 0.002, yS: 0.3775, troughY: 0.262, holes });
  const glow = glowMaterial(0x5d5853);

  // Grate rails across the front and back, angle iron welded to the sides.
  const L = zg0 - zg1;
  const slope = Math.atan2(yb - yf, L);
  for (const [zz, yy] of [[zg0 - 0.012, yf - 0.03], [zg1 + 0.012, yb - 0.03]] as const) {
    s.add(p.ssHeat, sheet(s, w - 0.0032, 0.003, 0.025), [w / 2, yy - 0.0015, zz]);
    s.add(p.ssHeat, sheet(s, w - 0.0032, 0.025, 0.003), [w / 2, yy - 0.014, zz + (zz > -0.3 ? 0.011 : -0.011)]);
  }

  const knobs: string[] = [];
  const ports: Port[] = [];
  const nb = 5;
  const bw = 0.019;
  for (let i = 0; i < n; i++) {
    const cx = 0.002 + sec * (i + 0.5);
    const gw = sec - 0.004;
    // The grate: a casting of five bars, each with a grease channel along its top between two lips, on three ribs.
    s.at([cx, yf, zg0], [-slope, 0, 0], () => {
      const pitch = (gw - bw) / (nb - 1);
      for (let k = 0; k < nb; k++) {
        const x = -gw / 2 + bw / 2 + k * pitch;
        if (s.fine(1)) {
          s.add(p.castIron, sheet(s, bw, 0.016, L, 0.002), [x, -0.0115, -L / 2]);
          for (const sx of [-1, 1]) s.add(p.castIron, sheet(s, 0.0045, 0.006, L, 0.0015), [x + sx * (bw / 2 - 0.00225), -0.003, -L / 2]);
          // The bar's front end is notched where its channel spills into the trough.
          if (s.fine(2)) s.add(p.castIron, rbox(s, bw, 0.012, 0.006, 0.002), [x, -0.017, -0.003]);
        } else s.add(p.castIron, box(bw, 0.02, L), [x, -0.01, -L / 2]);
      }
      for (const zz of [0.012, L / 2, L - 0.012]) s.add(p.castIron, sheet(s, gw, 0.012, 0.014, 0.002), [0, -0.025, -zz]);
      if (s.fine(2)) s.add(p.castIron, rbox(s, 0.03, 0.008, 0.02, 0.002), [0, -0.004, -L - 0.004], [0.5, 0, 0]);
    });
    // The radiant: two stainless leaves leaning together, a slot along the ridge held by bridges.
    const ry = 0.284;
    const hw = 0.058;
    const gap = 0.0045;
    const rl = L - 0.03;
    const rz = (zg0 + zg1) / 2;
    const ang = Math.atan2(0.026, hw - gap);
    const leaf = Math.hypot(hw - gap, 0.026);
    for (const sx of [-1, 1]) s.add(glow, sheet(s, leaf, 0.0014, rl, 0.0005), [cx + sx * (gap + (hw - gap) / 2), ry - 0.013, rz], [0, 0, -sx * ang]);
    if (s.fine(2)) for (let k = 0; k < 9; k++) s.add(glow, box(gap * 2 + 0.002, 0.003, 0.008), [cx, ry - 0.0005, rz - rl / 2 + ((k + 0.5) * rl) / 9]);
    if (s.fine(1)) for (const zz of [rz + rl / 2, rz - rl / 2]) s.add(p.ssHeat, box(hw * 2 + 0.01, 0.004, 0.003), [cx, ry - 0.026, zz]);
    // The burner: a pipe front to back, ports along its top, its venturi on the valve at the front.
    if (s.fine(1)) {
      s.add(p.castIron, cyl(s, br, 0.57, "z", 14), [cx, by, -0.385]);
      s.add(p.castIron, cyl(s, br, 0.03, "z", 14, 0.018), [cx, by, -0.085]);
      s.add(p.ssDull, cyl(s, 0.02, 0.0015, "z", 16), [cx, by, -0.069]);
      s.add(p.brass, rbox(s, 0.026, 0.03, 0.026, 0.002), [cx, ky, zf - 0.016]);
      s.add(p.brass, pipe(s, [[cx, ky, zf - 0.028], [cx, ky, -0.064], [cx, by, -0.072]], 0.004, 0.006, 8));
    }
    ports.push(...zLinePorts(i, by + br - 0.002, cx, -0.1, -0.66, portCount(s, 18), 0.05));
    knobs.push(valveKnob(s, `knob_${i}`, [cx, ky, zf]));
    dialTicks(s, [cx, ky, zf], 0.03, 7);
  }
  if (s.fine(1)) s.add(p.enamel, cyl(s, 0.0105, w - 0.06, "x", 12), [w / 2, 0.178, zf - 0.03]);
  dataPlate(s, [w - 0.06, ky + 0.002, zf], 0.05, 0.034);

  // The water pan: full width under the burners, its front lapping the panel's opening.
  const pw = w - 0.06;
  const pan = s.joint("grease_pan", [0.03, panY0, zf]);
  pan.add(p.ss, sheet(pan, pw, panH + 0.002, 0.0016, 0.0008), [pw / 2, panH / 2, 0.0008]);
  for (const hx of [0.14, pw - 0.14]) pullHandle0(pan, [hx, panH / 2, 0.0016], 0.09);
  if (s.fine(1)) {
    const dd = 0.62;
    pan.add(p.ssDull, box(pw - 0.016, 0.001, dd), [pw / 2, 0.002, -dd / 2]);
    for (const sx of [0.008, pw - 0.008]) pan.add(p.ssDull, box(0.001, panH - 0.01, dd), [sx, panH / 2 - 0.003, -dd / 2]);
    pan.add(p.ssDull, box(pw - 0.016, panH - 0.01, 0.001), [pw / 2, panH / 2 - 0.003, -dd]);
    pan.add(p.water, box(pw - 0.02, 0.012, dd - 0.01), [pw / 2, 0.009, -dd / 2]);
  }

  return (b) => {
    const g = gasHooks(b, ports, n, 0.045, knobs);
    const rest: V3 = [0.03, panY0, zf];
    g.glows = [glow];
    g.open = (f) => slideJoint(b.joints.get("grease_pan"), "z", f * 0.4, rest);
    g.set = (name, v) => {
      if (name === "heat") setGlow(glow, v);
    };
    return g;
  };
};

/* ------------------------------------------------------------------------------------------------------------ */
/* Fryer                                                                                                          */
/* ------------------------------------------------------------------------------------------------------------ */

const smooth = (a: number, b: number, x: number) => {
  const t = Math.min(1, Math.max(0, (x - a) / (b - a)));
  return t * t * (3 - 2 * t);
};

/**
 * A fry basket in its joint's frame: the origin at the middle of its front top edge, the basket hanging back (-z)
 * and down from it. Nickel-plated wire mesh on a heavier rod frame, a long handle with an insulated grip forward, a
 * pair of front hooks that rest on the fryer's front ledge when it is down, and a pair of rear hooks that catch the
 * hanger rod when it is up. `ledge` is the ledge's height above the basket's rim when down.
 */
function fryBasket(s: Shop, bw: number, bd: number, bh: number, ledge: number): void {
  const p = palette();
  const m = p.chrome;
  const [hx, r] = [bw / 2, 0.0028];
  const rim = (y: number, rr: number) => s.add(m, pipe(s, [[-hx, y, 0], [hx, y, 0], [hx, y, -bd], [-hx, y, -bd], [-hx, y, 0], [hx * 0.2, y, 0]], rr, 0.012, 8));
  rim(0, r);
  rim(-bh, r * 0.9);
  for (const [x, z] of [[-hx, 0], [hx, 0], [-hx, -bd], [hx, -bd]] as const) s.add(m, cyl(s, r * 0.85, bh, "y", 8), [x, -bh / 2, z]);
  if (s.fine(1)) {
    const pitch = s.fine(2) ? 0.009 : 0.02;
    const wr = s.fine(2) ? 0.0009 : 0.0012;
    wireGrid(s, m, bw, bh, pitch, wr, [0, -bh / 2, 0]);
    wireGrid(s, m, bw, bh, pitch, wr, [0, -bh / 2, -bd]);
    for (const sx of [-1, 1]) wireGrid(s, m, bd, bh, pitch, wr, [sx * hx, -bh / 2, -bd / 2], [0, Math.PI / 2, 0]);
    wireGrid(s, m, bw, bd, pitch, wr, [0, -bh, -bd / 2], [-Math.PI / 2, 0, 0]);
    // Two skids under the floor, so the mesh sits clear of the tube screen.
    for (const sx of [-0.5, 0.5]) s.add(m, cyl(s, r * 0.8, bd, "z", 8), [sx * hx, -bh - r, -bd / 2]);
  }
  // The handle: two rods up from the front rim, over the ledge and forward into the grip.
  const gy = ledge + 0.04;
  for (const sx of [-1, 1]) s.add(m, pipe(s, [[sx * 0.024, -0.02, 0.001], [sx * 0.024, 0, 0.002], [sx * 0.014, ledge + 0.012, 0.002], [sx * 0.007, gy, 0.045], [sx * 0.007, gy, 0.07]], 0.0032, 0.02, 10));
  s.add(p.ssBright, cyl(s, 0.0135, 0.009, "z", 16), [0, gy, 0.066]);
  s.add(p.bakelite, lathe(s, [
    [0, 0],
    [0.0125, 0],
    [0.0135, 0.02],
    [0.0128, 0.05],
    [0.0118, 0.058],
    [0.006, 0.062],
    [0, 0.062],
  ], 24).rotateX(Math.PI / 2), [0, gy, 0.07]);
  if (s.fine(2)) for (let k = 0; k < 5; k++) s.add(p.bakelite, torus(s, 0.0136, 0.0012, "z", Math.PI * 2, 16, 4), [0, gy, 0.078 + k * 0.007]);
  // Front hooks: up past the rim and out onto the ledge.
  for (const sx of [-1, 1]) s.add(m, pipe(s, [[sx * hx * 0.62, -0.03, 0.001], [sx * hx * 0.62, ledge + 0.0025, 0.002], [sx * hx * 0.62, ledge + 0.0025, 0.03]], 0.0024, 0.006, 8));
  // Rear hooks: up off the back rim, over and down, to hang on the rod.
  for (const sx of [-1, 1]) s.add(m, pipe(s, [[sx * hx * 0.55, -0.025, -bd], [sx * hx * 0.55, 0.02, -bd], [sx * hx * 0.55, 0.022, -bd - 0.016], [sx * hx * 0.55, 0.004, -bd - 0.02]], 0.0025, 0.006, 8));
}

/**
 * A 50 lb tube-fired floor fryer, 15.5 in wide: a stainless cabinet on 6 in legs, a door hinged on its left, the
 * frypot sunk in the top deck with four fire tubes crossing it below the oil line (the cold zone under them keeps
 * crumbs from scorching), a wire screen over the tubes for the baskets to sit on, the flue riser at the back with
 * its deflector, and the hanger rod in front of it. Behind the door: the inshot burners firing into the tubes, the
 * gas valve, the drain valve.
 */
const fryer: Builder = (s, e) => {
  const p = palette();
  const { w, d, h } = e;
  const legH = 0.152;
  const deck = 0.9;
  const zf = -0.03;
  const back = -d;
  const [px0, px1, pz0, pz1, pyB] = [0.037, w - 0.037, -0.125, -0.515, 0.575];
  const oilY = 0.78;
  const tubeY = 0.605;
  const doorY0 = 0.17;
  const doorY1 = 0.7;

  for (const lx of [0.045, w - 0.045])
    for (const lz of [zf - 0.05, back + 0.05]) {
      flangedLeg(s, lx, lz, legH + 0.002, 0.02);
      s.add(p.ssDull, sheet(s, 0.055, 0.005, 0.055, 0.001), [lx, legH - 0.0025, lz]);
    }
  // Cabinet: sides, back, bottom frame, the front frame round the door and the control panel above it.
  const bodyD = zf - back;
  for (const sx of [0.0008, w - 0.0008]) s.add(p.ss, sheet(s, 0.0016, deck - legH, bodyD, 0.0006), [sx, (deck + legH) / 2, (zf + back) / 2]);
  s.add(p.ssDull, sheet(s, w - 0.0032, deck - legH, 0.0016), [w / 2, (deck + legH) / 2, back + 0.0008]);
  s.add(p.ssDull, sheet(s, w - 0.0032, 0.0016, bodyD), [w / 2, legH + 0.0008, (zf + back) / 2]);
  for (const sx of [0.006, w - 0.006]) s.add(p.ss, sheet(s, 0.012, doorY1 - legH, 0.0016), [sx, (doorY1 + legH) / 2, zf - 0.0008]);
  s.add(p.ss, sheet(s, w, doorY0 - legH, 0.0016), [w / 2, (doorY0 + legH) / 2, zf - 0.0008]);
  const panY1 = deck - 0.014;
  s.add(p.ss, sheet(s, w, panY1 - doorY1, 0.0016), [w / 2, (panY1 + doorY1) / 2, zf - 0.0008]);
  // The deck's front nosing, rolled from the top sheet.
  s.add(p.ss, cyl(s, 0.012, w, "x", 24), [w / 2, deck - 0.012, zf - 0.006]);
  // The top deck round the frypot opening, and the frypot's radiused lip.
  s.add(p.ss, sheet(s, w, 0.0016, pz0 - zf + 0.006), [w / 2, deck - 0.0008, (pz0 + zf - 0.006) / 2]);
  s.add(p.ss, sheet(s, w, 0.0016, pz1 - back), [w / 2, deck - 0.0008, (pz1 + back) / 2]);
  for (const [x0, x1] of [[0, px0], [px1, w]] as const) s.add(p.ss, sheet(s, x1 - x0, 0.0016, pz0 - pz1), [(x0 + x1) / 2, deck - 0.0008, (pz0 + pz1) / 2]);
  if (s.fine(1)) {
    for (const z of [pz0, pz1]) s.add(p.ss, cyl(s, 0.003, px1 - px0, "x", 8), [w / 2, deck - 0.003, z]);
    for (const x of [px0, px1]) s.add(p.ss, cyl(s, 0.003, pz0 - pz1, "z", 8), [x, deck - 0.003, (pz0 + pz1) / 2]);
  }
  // The frypot: welded stainless, its floor a shallow V falling to the drain at the front.
  for (const x of [px0 - 0.001, px1 + 0.001]) s.add(p.ss, box(0.002, deck - pyB, pz0 - pz1), [x, (deck + pyB) / 2, (pz0 + pz1) / 2]);
  for (const z of [pz0 + 0.001, pz1 - 0.001]) s.add(p.ss, box(px1 - px0 + 0.004, deck - pyB + 0.02, 0.002), [w / 2, (deck + pyB - 0.02) / 2, z]);
  const vw = (px1 - px0) / 2;
  const va = Math.atan2(0.02, vw);
  for (const sx of [-1, 1]) s.add(p.ss, box(Math.hypot(vw, 0.02), 0.002, pz0 - pz1), [w / 2 + (sx * vw) / 2, pyB - 0.01, (pz0 + pz1) / 2], [0, 0, sx * va]);
  // The oil: 50 lb of it, up to the line.
  s.add(p.oil, box(px1 - px0 - 0.001, 0.002, pz0 - pz1 - 0.001), [w / 2, oilY, (pz0 + pz1) / 2]);
  // Fire tubes across the frypot below the oil, flattened ovals, flanged where they pass the front wall.
  const tubes = 4;
  const tx = (k: number) => px0 + ((px1 - px0) * (k + 0.5)) / tubes;
  for (let k = 0; k < tubes; k++) {
    // Hollow, so that the burner's flame shows in its mouth.
    const tl = pz0 - pz1 + 0.045;
    const tube = lathe(s, [[0.026, 0], [0.03, 0], [0.03, tl], [0.026, tl], [0.026, 0]], 20).rotateX(Math.PI / 2);
    s.add(p.ssHeat, tube, [tx(k), tubeY, pz1 - 0.015], undefined, [0.82, 1, 1]);
    if (s.fine(1)) s.add(p.ssHeat, torus(s, 0.03, 0.004, "z", Math.PI * 2, 20, 6), [tx(k), tubeY, pz0 + 0.03], undefined, [0.82, 1, 1]);
  }
  // The tube screen the baskets rest on.
  if (s.fine(1)) {
    wireGrid(s, p.ssDull, px1 - px0 - 0.004, pz0 - pz1 - 0.01, 0.03, 0.0018, [w / 2, 0.641, (pz0 + pz1) / 2], [-Math.PI / 2, 0, 0]);
    for (const z of [pz0 - 0.01, pz1 + 0.01]) s.add(p.ssDull, box(px1 - px0 - 0.004, 0.012, 0.003), [w / 2, 0.636, z]);
  }

  // Flue riser at the back, slotted at its top, and the deflector that turns the hot gases up and back.
  const [fx0, fx1, fz0, fz1] = [0.06, w - 0.06, -0.635, back + 0.012];
  const fTop = 1.105;
  s.add(p.ss, sheet(s, fx1 - fx0, fTop - deck, 0.0016), [(fx0 + fx1) / 2, (fTop + deck) / 2, fz0]);
  s.add(p.ss, sheet(s, fx1 - fx0, fTop - deck, 0.0016), [(fx0 + fx1) / 2, (fTop + deck) / 2, fz1]);
  for (const x of [fx0, fx1]) s.add(p.ss, sheet(s, 0.0016, fTop - deck, fz0 - fz1), [x, (fTop + deck) / 2, (fz0 + fz1) / 2]);
  s.add(p.enamel, box(fx1 - fx0 - 0.004, 0.002, fz0 - fz1 - 0.004), [(fx0 + fx1) / 2, fTop - 0.012, (fz0 + fz1) / 2]);
  if (s.fine(1)) for (let i = 0; i < 6; i++) s.add(p.ss, sheet(s, 0.006, 0.003, fz0 - fz1), [fx0 + ((fx1 - fx0) * (i + 0.5)) / 6, fTop - 0.0015, (fz0 + fz1) / 2]);
  s.add(p.ss, sheet(s, fx1 - fx0 + 0.02, 0.0016, 0.07), [(fx0 + fx1) / 2, h - 0.026, (fz0 + fz1) / 2 - 0.006], [-0.62, 0, 0]);
  for (const x of [fx0 - 0.008, fx1 + 0.008]) s.add(p.ss, sheet(s, 0.0016, 0.05, 0.05), [x, fTop + 0.006, (fz0 + fz1) / 2 + 0.006]);
  // The basket hanger: a rod on two brackets off the flue's face, its ends turned up.
  const [rodY, rodZ] = [0.98, -0.5];
  s.add(p.ssBright, pipe(s, [[0.02, rodY + 0.02, rodZ], [0.02, rodY, rodZ], [w - 0.02, rodY, rodZ], [w - 0.02, rodY + 0.02, rodZ]], 0.0055, 0.008, 10));
  for (const x of [0.1, w - 0.1]) s.add(p.ss, sheet(s, 0.004, 0.022, fz0 - rodZ + 0.012, 0.001), [x, rodY + 0.004, (fz0 + rodZ) / 2 - 0.004], [-0.12, 0, 0]);

  // The thermostat dial on the control panel, and the frypot's drain behind the door.
  const ky = (doorY1 + panY1) / 2;
  knob(s, p.bakelite, p.ssBright, [w / 2, ky, zf], 0.022);
  dialTicks(s, [w / 2, ky, zf], 0.032, 11);
  if (s.fine(2)) for (const sx of [0.02, w - 0.02]) screw(s, p.ssBright, [sx, ky, zf], 0.0026);

  // Behind the door: inshot burners at the tube mouths on their manifold, the gas valve, the drain valve.
  const ports: Port[] = [];
  for (let k = 0; k < tubes; k++) {
    if (s.fine(1)) {
      s.add(p.castIron, cyl(s, 0.0095, 0.026, "z", 14), [tx(k), tubeY, pz0 + 0.046]);
      s.add(p.castIron, pipe(s, [[tx(k), tubeY - 0.065, pz0 + 0.05], [tx(k), tubeY, pz0 + 0.05], [tx(k), tubeY, pz0 + 0.056]], 0.006, 0.012, 8));
    }
    const np = portCount(s, 6);
    for (let i = 0; i < np; i++) {
      const a = (i / np) * Math.PI * 2;
      ports.push({ at: [tx(k) + Math.cos(a) * 0.0065, tubeY + Math.sin(a) * 0.0065, pz0 + 0.032], dir: [Math.cos(a) * 0.3, Math.sin(a) * 0.3, -1], burner: 0 });
    }
  }
  ports.push({ at: [tx(0) - 0.022, tubeY - 0.012, pz0 + 0.04], dir: [0.4, 0.3, -1], burner: 1 });
  let gasKnob = "gas_valve_knob";
  if (s.fine(1)) {
    s.add(p.enamel, cyl(s, 0.011, px1 - px0, "x", 12), [w / 2, tubeY - 0.065, pz0 + 0.05]);
    // The combination gas valve, its knob on top, the supply in from the back and up to the manifold.
    const [gx, gy, gz] = [0.1, 0.46, -0.16];
    s.add(p.alu, rbox(s, 0.08, 0.05, 0.05, 0.004), [gx, gy, gz]);
    s.add(p.alu, rbox(s, 0.03, 0.03, 0.03, 0.003), [gx, gy + 0.035, gz]);
    s.add(p.enamel, pipe(s, [[gx - 0.04, gy, gz], [0.03, gy, gz], [0.03, gy, back + 0.05], [0.03, 0.2, back + 0.05]], 0.011, 0.03, 12));
    s.add(p.enamel, pipe(s, [[gx + 0.04, gy, gz], [gx + 0.07, gy, gz], [gx + 0.07, tubeY - 0.065, gz], [gx + 0.07, tubeY - 0.065, pz0 + 0.05]], 0.009, 0.025, 12));
    if (s.fine(2)) s.add(p.copper, pipe(s, [[gx, gy + 0.025, gz + 0.025], [gx, gy + 0.08, gz + 0.06], [tx(0) - 0.022, tubeY - 0.03, pz0 + 0.06], [tx(0) - 0.022, tubeY - 0.014, pz0 + 0.045]], 0.0022, 0.02, 6));
    const gk = s.joint("gas_valve_knob", [gx, gy + 0.05, gz], [-Math.PI / 2, 0, 0]);
    gk.add(p.bakelite, cyl(gk, 0.018, 0.016, "z", 20), [0, 0, 0.008]);
    gk.add(p.bakelite, rbox(gk, 0.008, 0.034, 0.014, 0.002), [0, 0, 0.022]);
    dataPlate(s, [w - 0.0016, 0.4, -0.25], 0.08, 0.05, [0, -Math.PI / 2, 0]);
    // The drain: 1.25 in from the frypot's front floor, a ball valve, its lever handle forward.
    s.add(p.ss, pipe(s, [[w / 2, pyB - 0.012, pz0 - 0.02], [w / 2, pyB - 0.04, pz0 - 0.02], [w / 2, pyB - 0.04, -0.07], [w / 2, 0.42, -0.07]], 0.016, 0.03, 14));
    s.add(p.brass, rbox(s, 0.05, 0.045, 0.05, 0.006), [w / 2, 0.4, -0.07]);
    s.add(p.ssBright, cyl(s, 0.016, 0.03, "y", 14), [w / 2, 0.365, -0.07]);
    s.add(p.bakelite, rbox(s, 0.016, 0.01, 0.13, 0.004), [w / 2 + 0.03, 0.405, -0.07 + 0.01], [0, -0.9, 0]);
  } else gasKnob = "";

  // The door: double-skinned, hinged on its left edge, a bar pull up its right.
  const dw = w - 0.008;
  const dh = doorY1 - doorY0 - 0.004;
  const door = s.joint("front_door", [0.004, doorY0 + 0.002, zf]);
  door.add(p.ss, rbox(door, dw, dh, 0.018, 0.003), [dw / 2, dh / 2, 0.009]);
  if (s.fine(1)) {
    door.add(p.ss, sheet(door, dw - 0.05, dh - 0.06, 0.0012, 0.003), [dw / 2, dh / 2, 0.0184]);
    door.add(p.ssDull, sheet(door, dw - 0.03, dh - 0.03, 0.002), [dw / 2, dh / 2, -0.0008]);
    for (const y of [0.04, dh - 0.04]) door.add(p.ssBright, cyl(door, 0.005, 0.05, "y", 12), [-0.001, y, 0.004]);
  }
  door.add(p.ssBright, pipe(door, [[dw - 0.028, dh * 0.3, 0.018], [dw - 0.028, dh * 0.3, 0.0265], [dw - 0.028, dh * 0.62, 0.0265], [dw - 0.028, dh * 0.62, 0.018]], 0.0055, 0.006, 10));

  // The baskets, each a joint posed by its lift.
  const bw = (px1 - px0 - 0.012) / 2;
  const [bd, bh] = [0.33, 0.135];
  const downY = 0.79;
  const downZ = pz0 - 0.01;
  const bx = [px0 + 0.004 + bw / 2, px1 - 0.004 - bw / 2];
  const names = ["basket_left_lift", "basket_right_lift"];
  for (let i = 0; i < 2; i++) fryBasket(s.joint(names[i]!, [bx[i]!, downY, downZ]), bw, bd, bh, deck - downY);
  const hookZ = bd + 0.0105;
  const hung: [number, number] = [rodY + 0.0055 + 0.0025 - 0.021, rodZ + hookZ];

  return (b) => {
    const lift = [0, 0];
    const flames = new FlameSet(ports, 0.03);
    b.group.add(flames.mesh);
    flames.setAll(0);
    flames.setFlow(1, 0.12);
    if (gasKnob) {
      const set = flames.setFlow.bind(flames);
      flames.setFlow = (i: number, f: number) => {
        set(i, f);
        if (i === 0) turnJoint(b.joints.get(gasKnob), "z", valveAngle(f));
      };
    }
    const pose = (i: number) => {
      const j = b.joints.get(names[i]!);
      if (!j) return;
      const t = lift[i]!;
      const a = smooth(0, 0.65, t);
      const c = smooth(0.85, 1, t);
      const y = downY + (rodY + 0.03 - downY) * a - (rodY + 0.03 - hung[0]) * c;
      const z = downZ + (hung[1] - downZ) * smooth(0.65, 0.85, t);
      j.position.set(bx[i]!, y, z);
      turnJoint(j, "x", -0.12 * Math.sin(Math.PI * smooth(0.05, 0.75, t)));
    };
    const bub = bx.map((x) => new Bubbles(p.hdpe, s.fine(2) ? 70 : 30, [x, oilY + 0.001, downZ - bd / 2], bw - 0.012, bd - 0.02, 0.0045));
    for (const q of bub) b.group.add(q.mesh);
    return {
      flames,
      burners: 1,
      open: (f) => turnJoint(b.joints.get("front_door"), "y", -f * 1.9),
      set: (name, v) => {
        const k = name === "lift_l" ? 0 : name === "lift_r" ? 1 : -1;
        if (k < 0) return;
        lift[k] = Math.min(1, Math.max(0, v));
        pose(k);
      },
      tick: (t, _dt, load) => {
        flames.tick(t);
        for (let i = 0; i < 2; i++) bub[i]!.tick(t, (1 - smooth(0.05, 0.3, lift[i]!)) * (0.35 + 0.65 * load));
      },
    };
  };
};

/* ------------------------------------------------------------------------------------------------------------ */
/* Salamander                                                                                                     */
/* ------------------------------------------------------------------------------------------------------------ */

/**
 * A 36 in infrared salamander, hung on a range's riser by the brackets on its back: an open cavity under two
 * ceramic infrared burners (perforated plaques the gas burns on the face of, so they glow rather than flame), a
 * spring-balanced carriage that rides up and down in the side walls on a lever through a notched slot, and on it a
 * chrome broiling grid that rolls out to load.
 */
const salamander: Builder = (s, e) => {
  const p = palette();
  const { w, d, h } = e;
  const zf = -0.03;
  const back = -d;
  const top = 0.4;
  const [cx0, cx1, cy0, cy1, cz1] = [0.032, w - 0.124, 0.035, 0.36, back + 0.03];
  const colX = (cx1 + w) / 2;
  const yLow = 0.08;
  const travel = 0.2;
  const glow = glowMaterial(0xbdb2a2);

  // The shell: top, sides, back and bottom, each double-skinned, the front a frame round the opening.
  s.add(p.ss, sheet(s, w, 0.0016, d + zf), [w / 2, top - 0.0008, (zf + back) / 2]);
  for (const sx of [0.0008, w - 0.0008]) s.add(p.ss, sheet(s, 0.0016, top, d + zf, 0.0006), [sx, top / 2, (zf + back) / 2]);
  s.add(p.ssDull, sheet(s, w - 0.0032, top, 0.0016), [w / 2, top / 2, back + 0.0008]);
  s.add(p.ss, sheet(s, w, 0.0016, d + zf), [w / 2, 0.0008, (zf + back) / 2]);
  s.add(p.ss, sheet(s, w, top - cy1, 0.0016), [w / 2, (top + cy1) / 2, zf - 0.0008]);
  s.add(p.ss, sheet(s, cx0, cy1, 0.0016), [cx0 / 2, cy1 / 2, zf - 0.0008]);
  s.add(p.ss, sheet(s, w - cx1, cy1, 0.0016), [(w + cx1) / 2, cy1 / 2, zf - 0.0008]);
  s.add(p.ss, sheet(s, cx1 - cx0, cy0, 0.0016), [(cx0 + cx1) / 2, cy0 / 2, zf - 0.0008]);
  // A rolled lip along the top of the opening, where the heat comes out.
  s.add(p.ss, cyl(s, 0.006, cx1 - cx0, "x", 16), [(cx0 + cx1) / 2, cy1 + 0.002, zf + 0.003]);
  // The cavity, heat-tinted: walls with the carriage's guide channels, floor and back.
  for (const x of [cx0, cx1]) {
    s.add(p.ssHeat, box(0.0016, cy1 - cy0, cz1 - zf), [x, (cy0 + cy1) / 2, (zf + cz1) / 2]);
    if (s.fine(1)) for (const z of [zf - 0.03, cz1 + 0.03]) s.add(p.ssHeat, box(0.01, cy1 - cy0 - 0.03, 0.014), [x + (x === cx0 ? 0.005 : -0.005), (cy0 + cy1) / 2, z]);
  }
  s.add(p.ssHeat, box(cx1 - cx0, cy1 - cy0, 0.0016), [(cx0 + cx1) / 2, (cy0 + cy1) / 2, cz1]);
  s.add(p.ssHeat, box(cx1 - cx0, 0.0016, cz1 - zf), [(cx0 + cx1) / 2, cy0 + 0.0008, (zf + cz1) / 2]);
  // The drip pan on the floor, its front turned up into a lip with a pull.
  s.add(p.ssDull, sheet(s, cx1 - cx0 - 0.01, 0.022, 0.0012), [(cx0 + cx1) / 2, cy0 + 0.012, zf - 0.006]);
  s.add(p.ssDull, box(cx1 - cx0 - 0.01, 0.0012, cz1 - zf - 0.01), [(cx0 + cx1) / 2, cy0 + 0.003, (zf + cz1) / 2]);
  pullHandle0(s, [(cx0 + cx1) / 2, cy0 + 0.012, zf - 0.0054], 0.12);

  // The burners: two ceramic plaques in stainless frames across the ceiling, tilted to throw their heat forward.
  const bwid = (cx1 - cx0 - 0.03) / 2;
  const bdep = 0.28;
  for (let i = 0; i < 2; i++) {
    const bx = cx0 + 0.01 + bwid / 2 + i * (bwid + 0.01);
    s.at([bx, cy1 - 0.03, (zf + cz1) / 2 - 0.02], [0.1, 0, 0], () => {
      s.add(p.ssHeat, sheet(s, bwid, 0.03, bdep, 0.001), [0, 0.012, 0]);
      s.add(p.ss, sheet(s, bwid + 0.006, 0.008, bdep + 0.006, 0.001), [0, -0.0015, 0]);
      if (s.fine(2)) {
        const [nx, nz] = [5, 4];
        const tw = (bwid - 0.02) / nx;
        const td = (bdep - 0.02) / nz;
        for (let a = 0; a < nx; a++)
          for (let c = 0; c < nz; c++) s.add(glow, box(tw - 0.003, 0.006, td - 0.003), [-bwid / 2 + 0.01 + tw * (a + 0.5), -0.003, -bdep / 2 + 0.01 + td * (c + 0.5)]);
      } else s.add(glow, box(bwid - 0.02, 0.006, bdep - 0.02), [0, -0.003, 0]);
    });
  }
  if (s.fine(1)) s.add(p.ssHeat, box(cx1 - cx0, 0.0016, cz1 - zf), [(cx0 + cx1) / 2, cy1 - 0.0008, (zf + cz1) / 2]);

  // The control column: two burner valves, and the carriage lever's slot with its four notches.
  const knobs = [valveKnob(s, "knob_0", [colX, 0.325, zf]), valveKnob(s, "knob_1", [colX, 0.262, zf])];
  for (const y of [0.325, 0.262]) dialTicks(s, [colX, y, zf], 0.029, 7);
  const slotX = colX;
  const sy0 = yLow - 0.012 - 0.012;
  const sy1 = yLow + travel - 0.012 + 0.012;
  s.add(p.enamel, box(0.012, sy1 - sy0, 0.0006), [slotX, (sy0 + sy1) / 2, zf + 0.0003]);
  for (let k = 0; k < 4; k++) s.add(p.enamel, box(0.014, 0.012, 0.0006), [slotX - 0.012, yLow - 0.012 + (k * travel) / 3, zf + 0.0003]);
  if (s.fine(2)) for (const y of [0.04, 0.375]) for (const x of [0.016, w - 0.016]) screw(s, p.ssBright, [x, y, zf], 0.0026);

  // Flue on the top at the back, louvred; louvres for the cooling air down each side.
  s.add(p.ss, sheet(s, w - 0.2, h - top, 0.11, 0.001), [w / 2, (top + h) / 2 - 0.0008, back + 0.07]);
  louvres(s, [w / 2, (top + h) / 2, back + 0.125], w - 0.26, h - top - 0.012, 3);
  for (const [x, ry] of [[0, -Math.PI / 2], [w, Math.PI / 2]] as const) louvres(s, [x, 0.2, back + 0.15], 0.14, 0.12, 6, undefined, [0, ry, 0]);
  dataPlate(s, [w, 0.3, back + 0.15], 0.07, 0.04, [0, Math.PI / 2, 0]);
  // The range-mount brackets on the back: hanger channels with their slots.
  for (const x of [0.16, w - 0.16]) {
    s.add(p.ssDull, sheet(s, 0.05, 0.3, 0.004, 0.001), [x, 0.2, back + 0.002]);
    if (s.fine(1)) s.add(p.enamel, box(0.012, 0.05, 0.0006), [x, 0.31, back - 0.0003]);
  }

  // The carriage, at its lowest when its joint is at rest: side rails in the guides, the lever through the slot.
  const car = s.joint("broiler_rack_lift", [0, yLow, 0]);
  const rz0 = zf - 0.012;
  const rl = rz0 - (cz1 + 0.01);
  for (const x of [cx0 + 0.012, cx1 - 0.012]) {
    car.add(p.ssDull, sheet(car, 0.014, 0.003, rl, 0.001), [x, -0.0135, rz0 - rl / 2]);
    car.add(p.ssDull, sheet(car, 0.0025, 0.03, rl, 0.001), [x + (x < w / 2 ? -0.0075 : 0.0075), -0.025, rz0 - rl / 2]);
  }
  car.add(p.ssDull, sheet(car, cx1 - cx0 - 0.03, 0.02, 0.003, 0.001), [(cx0 + cx1) / 2, -0.022, cz1 + 0.012]);
  car.add(p.ssBright, pipe(car, [[cx1 - 0.01, -0.02, -0.12], [slotX, -0.02, -0.12], [slotX, -0.012, zf - 0.02], [slotX, -0.012, zf + 0.004]], 0.004, 0.012, 10));
  car.add(p.bakelite, ball(car, 0.012, 18), [slotX, -0.012, zf + 0.014], undefined, [1, 1, 0.8]);
  // The grid: chrome rods front to back on a heavier frame, three cross bars under, a pull at the front.
  const grid = car.joint("broiler_grid", [0, 0, 0]);
  const [gx0, gx1, gz0, gz1] = [cx0 + 0.02, cx1 - 0.02, zf - 0.016, cz1 + 0.02];
  grid.add(p.chrome, pipe(grid, [[gx0, -0.006, gz0], [gx1, -0.006, gz0], [gx1, -0.006, gz1], [gx0, -0.006, gz1], [gx0, -0.006, gz0], [w / 2, -0.006, gz0]], 0.004, 0.01, 8));
  const rods = s.fine(1) ? Math.round((gx1 - gx0) / 0.02) : 8;
  for (let i = 1; i < rods; i++) grid.add(p.chrome, cyl(grid, 0.0025, gz0 - gz1, "z", 8), [gx0 + ((gx1 - gx0) * i) / rods, -0.0025, (gz0 + gz1) / 2]);
  for (const f of [0.25, 0.5, 0.75]) grid.add(p.chrome, cyl(grid, 0.0035, gx1 - gx0, "x", 8), [w / 2, -0.0085, gz0 + (gz1 - gz0) * f]);
  grid.add(p.chrome, pipe(grid, [[w / 2 - 0.07, -0.006, gz0], [w / 2 - 0.07, -0.004, gz0 + 0.022], [w / 2 + 0.07, -0.004, gz0 + 0.022], [w / 2 + 0.07, -0.006, gz0]], 0.0045, 0.01, 10));

  return (b) => {
    let rack = 0.5;
    let out = 0;
    const pose = () => {
      const j = b.joints.get("broiler_rack_lift");
      if (j) j.position.set(0, yLow + rack * travel, 0);
      const g = b.joints.get("broiler_grid");
      if (g) g.position.set(0, 0, out * 0.22);
    };
    pose();
    return {
      glows: [glow],
      open: (f) => {
        out = Math.min(1, Math.max(0, f));
        pose();
      },
      set: (name, v) => {
        if (name === "heat") {
          setGlow(glow, v);
          for (const k of knobs) turnJoint(b.joints.get(k), "z", valveAngle(v));
        } else if (name === "rack") {
          rack = Math.round(Math.min(1, Math.max(0, v)) * 3) / 3;
          pose();
        }
      },
    };
  };
};

/* ------------------------------------------------------------------------------------------------------------ */
/* Induction                                                                                                      */
/* ------------------------------------------------------------------------------------------------------------ */

/** The marks printed on the underside of a ceramic glass top for one zone: its rim circle, an inner ring, centre ticks. */
function zoneRings(s: Shop, x: number, y: number, z: number, r: number): void {
  const p = palette();
  const n = s.seg(96, 24);
  const flat = (r0: number, r1: number) => new THREE.RingGeometry(r0, r1, n).rotateX(-Math.PI / 2);
  s.add(p.gasket, flat(r - 0.0022, r), [x, y + 0.0006, z]);
  if (!s.fine(2)) return;
  s.add(p.gasket, flat(r * 0.52 - 0.0012, r * 0.52), [x, y + 0.0006, z]);
  for (let k = 0; k < 4; k++) {
    const a = (k / 4) * Math.PI * 2 + Math.PI / 4;
    s.add(p.gasket, box(0.0016, 0.0002, 0.012), [x + Math.cos(a) * (r + 0.012), y + 0.0006, z + Math.sin(a) * (r + 0.012)], [0, -a + Math.PI / 2, 0]);
  }
}

/**
 * An induction control: a power knob with its dial, over it a row of LEDs that light one by one with the zone's
 * power and the zone's lamp; the lit LEDs are one instanced draw whose count follows the level. Faces +z at `at`.
 */
interface ZoneCtl {
  lamp: THREE.MeshStandardMaterial;
  leds: THREE.InstancedMesh;
  knob: string;
}
function zoneControl(s: Shop, i: number, at: V3, ledsAbove = 0.036): ZoneCtl {
  const p = palette();
  const [x, y, z] = at;
  const knobName = valveKnob(s, `knob_${i}`, at, 0.019);
  dialTicks(s, at, 0.028, 9);
  const n = 8;
  const pitch = 0.007;
  const ly = y + ledsAbove;
  const leds = new THREE.InstancedMesh(new THREE.BoxGeometry(0.0045, 0.0028, 0.0012), p.ledRed, n);
  const m = new THREE.Matrix4();
  for (let k = 0; k < n; k++) {
    const lx = x - ((n - 1) * pitch) / 2 + k * pitch;
    m.makeTranslation(lx, ly, z + 0.0009);
    leds.setMatrixAt(k, m);
    if (s.fine(2)) s.add(p.charcoal, box(0.0045, 0.0028, 0.0006), [lx, ly, z + 0.0003]);
  }
  leds.count = 0;
  leds.castShadow = false;
  s.add(p.screen, box(n * pitch + 0.008, 0.009, 0.0004), [x, ly, z + 0.0001]);
  const lamp = lampMaterial(0xff3a1a);
  s.add(lamp, cyl(s, 0.0035, 0.002, "z", 12), [x + (n * pitch) / 2 + 0.01, ly, z + 0.001]);
  if (s.fine(2)) s.add(p.ssBright, torus(s, 0.0042, 0.0008, "z", Math.PI * 2, 12, 4), [x + (n * pitch) / 2 + 0.01, ly, z + 0.0008]);
  return { lamp, leds, knob: knobName };
}
function setZone(c: ZoneCtl, b: Built, v: number): void {
  const k = Math.min(1, Math.max(0, v));
  c.leds.count = Math.ceil(k * 8 - 1e-6);
  c.lamp.emissiveIntensity = k > 0.001 ? 2.2 : 0;
  turnJoint(b.joints.get(c.knob), "z", k <= 0.001 ? 0 : -(0.25 + k * 1.1 * Math.PI));
}

/** A ceramic glass top `w` by `d` centred at (x, z) on y, in a stainless bezel: black, glossy, its edges eased. */
function glassTop(s: Shop, x: number, z: number, w: number, d: number, y: number, t = 0.006): void {
  const p = palette();
  s.add(p.bakelite, rbox(s, w, t, d, 0.0025), [x, y - t / 2, z]);
  const b = 0.004;
  for (const sz of [-1, 1]) s.add(p.ss, sheet(s, w + 2 * b, t - 0.0015, b, 0.0008), [x, y - t / 2 - 0.00075, z + sz * (d / 2 + b / 2)]);
  for (const sx of [-1, 1]) s.add(p.ss, sheet(s, b, t - 0.0015, d, 0.0008), [x + sx * (w / 2 + b / 2), y - t / 2 - 0.00075, z]);
}

/**
 * A countertop dual-zone induction range, the zones one behind the other: a ceramic glass top on a stainless case
 * on rubber feet, the generators' cooling air drawn in through the sides and blown out of the back.
 */
const inductionHob: Builder = (s, e) => {
  const p = palette();
  const { w, d, h } = e;
  const zf = -0.03;
  const back = -d;
  const footH = 0.016;
  const caseTop = h - 0.006;
  for (const x of [0.04, w - 0.04])
    for (const z of [zf - 0.04, back + 0.04]) {
      s.add(p.rubber, cyl(s, 0.014, footH - 0.002, "y", 16, 0.012), [x, (footH - 0.002) / 2, z]);
      s.add(p.ssDull, cyl(s, 0.009, 0.003, "y", 12), [x, footH - 0.0015, z]);
    }
  // The case: a wrap of stainless round a black base pan, the front a control fascia.
  s.add(p.ss, rbox(s, w, caseTop - footH, d + zf, 0.004), [w / 2, (caseTop + footH) / 2, (zf + back) / 2]);
  s.add(p.charcoal, box(w - 0.01, 0.003, d + zf - 0.01), [w / 2, footH - 0.0005, (zf + back) / 2]);
  glassTop(s, w / 2, (zf + back) / 2, w - 0.012, d + zf - 0.012, h);
  const zones: [number, number][] = [[w / 2, -0.2], [w / 2, -0.51]];
  for (const [x, z] of zones) zoneRings(s, x, h, z, 0.135);
  // Intake louvres down both sides, the exhaust grille across the back.
  for (const [x, ry] of [[0, -Math.PI / 2], [w, Math.PI / 2]] as const) louvres(s, [x, (caseTop + footH) / 2, -0.25], 0.26, 0.06, 6, undefined, [0, ry, 0]);
  louvres(s, [w / 2, (caseTop + footH) / 2, back], w - 0.08, 0.07, 7, undefined, [0, Math.PI, 0]);
  dataPlate(s, [w * 0.25, (caseTop + footH) / 2, back], 0.06, 0.035, [0, Math.PI, 0]);
  // The controls: a knob, LEDs and lamp per zone; a rocker for the mains between them.
  const ky = 0.058;
  const ctl = [zoneControl(s, 0, [0.1, ky, zf], 0.04), zoneControl(s, 1, [w - 0.1, ky, zf], 0.04)];
  s.add(p.plastic, rbox(s, 0.022, 0.03, 0.006, 0.002), [w / 2, ky, zf + 0.003]);
  s.add(p.ledGreen, rbox(s, 0.016, 0.012, 0.01, 0.002), [w / 2, ky + 0.006, zf + 0.004], [-0.15, 0, 0]);
  if (s.fine(2)) {
    // Which zone each knob runs: two rings, front or back one filled.
    for (const [i, x] of [[0, 0.1], [1, w - 0.1]] as const)
      for (const [k, dy] of [[0, -0.006], [1, 0.006]] as const) {
        const cx = x + 0.045;
        s.add(p.enamel, (k === i ? new THREE.CircleGeometry(0.004, 16) : new THREE.RingGeometry(0.0032, 0.0042, 16)), [cx, ky + dy, zf + 0.0003]);
      }
  }

  return (b) => {
    for (const c of ctl) b.group.add(c.leds);
    return {
      set: (name, v) => {
        const m = /^zone(\d)$/.exec(name);
        const c = m ? ctl[Number(m[1])] : undefined;
        if (c) setZone(c, b, v);
      },
    };
  };
};

/**
 * A 36 in floor induction suite: four zones under one glass top, the generators in the cabinet under it on a
 * shelf, the cabinet's two doors hinged on their outer edges, the control fascia across the top front, 6 in legs.
 */
const inductionSuite: Builder = (s, e) => {
  const p = palette();
  const { w, d, h } = e;
  const legH = 0.152;
  const zf = -0.03;
  const back = -d;
  const topY = h;
  const deckY = h - 0.008;
  const panY0 = 0.77;
  const doorY0 = legH + 0.07;
  const doorY1 = panY0 - 0.006;
  for (const lx of [0.05, w - 0.05])
    for (const lz of [zf - 0.05, back + 0.05]) {
      flangedLeg(s, lx, lz, legH + 0.002, 0.02);
      s.add(p.ssDull, sheet(s, 0.055, 0.005, 0.055, 0.001), [lx, legH - 0.0025, lz]);
    }
  // Body: sides, back, bottom, top deck; the fascia proud of the doors, a kick plate under them.
  const bodyD = zf - back;
  for (const sx of [0.0008, w - 0.0008]) s.add(p.ss, sheet(s, 0.0016, deckY - legH, bodyD, 0.0006), [sx, (deckY + legH) / 2, (zf + back) / 2]);
  s.add(p.ssDull, sheet(s, w - 0.0032, deckY - legH, 0.0016), [w / 2, (deckY + legH) / 2, back + 0.0008]);
  s.add(p.ssDull, sheet(s, w - 0.0032, 0.0016, bodyD), [w / 2, legH + 0.0008, (zf + back) / 2]);
  s.add(p.ss, sheet(s, w, 0.0016, bodyD + 0.002), [w / 2, deckY - 0.0008, (zf + back) / 2 + 0.001]);
  s.add(p.ss, rbox(s, w, deckY - panY0, 0.012, 0.003), [w / 2, (deckY + panY0) / 2, zf - 0.006]);
  s.add(p.ss, sheet(s, w, doorY0 - legH, 0.0016), [w / 2, (doorY0 + legH) / 2, zf - 0.008]);
  louvres(s, [w / 2, (doorY0 + legH) / 2, zf - 0.0072], w - 0.2, 0.04, 4);
  for (const sx of [0.008, w - 0.008]) s.add(p.ss, sheet(s, 0.016, doorY1 - doorY0 + 0.006, 0.0016), [sx, (doorY1 + doorY0) / 2, zf - 0.0008]);
  s.add(p.ss, sheet(s, 0.016, doorY1 - doorY0 + 0.006, 0.0016), [w / 2, (doorY1 + doorY0) / 2, zf - 0.0028]);
  // The glass, set in the deck, its four zones' marks; the exhaust grille behind it.
  const gz1 = back + 0.05;
  glassTop(s, w / 2, (zf - 0.012 + gz1) / 2, w - 0.03, zf - 0.012 - gz1, topY);
  const zones: [number, number][] = [[0.24, -0.24], [0.674, -0.24], [0.24, -0.573], [0.674, -0.573]];
  for (const [x, z] of zones) zoneRings(s, x, topY, z, 0.14);
  s.add(p.enamel, box(w - 0.1, 0.001, 0.026), [w / 2, deckY + 0.0005, back + 0.027]);
  if (s.fine(1)) for (let i = 0; i < 24; i++) s.add(p.ss, sheet(s, 0.003, 0.004, 0.026, 0.0008), [0.05 + ((w - 0.1) * (i + 0.5)) / 24, deckY + 0.002, back + 0.027]);
  // The fascia: one control per zone, left pair for the left zones; the mains rocker and a data plate.
  const ky = 0.805;
  // Zones 0 and 1 are the front pair, 2 and 3 behind them: the back zones' knobs are the outer ones.
  const kx = [0.33, 0.584, 0.12, 0.794];
  const ctl = [0, 1, 2, 3].map((i) => zoneControl(s, i, [kx[i]!, ky, zf], 0.044));
  if (s.fine(2))
    for (let i = 0; i < 4; i++)
      for (let k = 0; k < 4; k++) {
        const [cx, cy] = [kx[i]! + 0.04 + (k % 2) * 0.009, ky + (k < 2 ? -0.005 : 0.005)];
        s.add(p.enamel, k === i ? new THREE.CircleGeometry(0.0036, 14) : new THREE.RingGeometry(0.0028, 0.0037, 14), [cx, cy, zf + 0.0003]);
      }
  s.add(p.plastic, rbox(s, 0.024, 0.034, 0.006, 0.002), [w / 2, ky + 0.01, zf + 0.003]);
  s.add(p.ledGreen, rbox(s, 0.018, 0.014, 0.01, 0.002), [w / 2, ky + 0.016, zf + 0.004], [-0.15, 0, 0]);
  dataPlate(s, [w - 0.0016, 0.6, -0.3], 0.08, 0.05, [0, -Math.PI / 2, 0]);
  // Inside: an undershelf, and the generator housings under the glass with their fans.
  if (s.fine(1)) {
    s.add(p.ssDull, sheet(s, w - 0.01, 0.0016, bodyD - 0.04), [w / 2, 0.38, (zf + back) / 2 - 0.01]);
    for (const x of [0.25, w - 0.25]) {
      s.add(p.charcoal, rbox(s, 0.34, 0.09, 0.6, 0.004), [x, deckY - 0.06, -0.42]);
      s.add(p.enamel, cyl(s, 0.045, 0.004, "z", 24), [x, deckY - 0.06, -0.118]);
      if (s.fine(2)) for (let k = 0; k < 5; k++) s.add(p.charcoal, torus(s, 0.012 + k * 0.008, 0.0012, "z", Math.PI * 2, 24, 4), [x, deckY - 0.06, -0.115]);
    }
  }
  // Doors: hinged on their outer edges, pulls at the meeting stiles.
  const dw = w / 2 - 0.012;
  const dh = doorY1 - doorY0;
  for (const [name, hx, sx] of [["door_left", 0.008, 1], ["door_right", w - 0.008, -1]] as const) {
    const j = s.joint(name, [hx, doorY0, zf]);
    j.add(p.ss, rbox(j, dw, dh, 0.018, 0.003), [(sx * dw) / 2, dh / 2, 0.009]);
    if (s.fine(1)) {
      j.add(p.ssDull, sheet(j, dw - 0.03, dh - 0.03, 0.002), [(sx * dw) / 2, dh / 2, -0.0008]);
      for (const y of [0.05, dh - 0.05]) j.add(p.ssBright, cyl(j, 0.005, 0.05, "y", 12), [0, y, 0.004]);
    }
    const px = sx * (dw - 0.03);
    j.add(p.ssBright, pipe(j, [[px, dh * 0.62, 0.018], [px, dh * 0.62, 0.027], [px, dh * 0.86, 0.027], [px, dh * 0.86, 0.018]], 0.0055, 0.006, 10));
  }

  return (b) => {
    for (const c of ctl) b.group.add(c.leds);
    return {
      open: (f) => {
        turnJoint(b.joints.get("door_left"), "y", -f * 1.75);
        turnJoint(b.joints.get("door_right"), "y", f * 1.75);
      },
      set: (name, v) => {
        const m = /^zone(\d)$/.exec(name);
        const c = m ? ctl[Number(m[1])] : undefined;
        if (c) setZone(c, b, v);
      },
    };
  };
};

/* ------------------------------------------------------------------------------------------------------------ */
/* Tilting skillet                                                                                                */
/* ------------------------------------------------------------------------------------------------------------ */

/** Run `f` in a joint with its parts placed in the fixture frame rather than about the pivot. */
function inFixtureFrame(j: Shop, pivot: V3, f: () => void): void {
  j.push([-pivot[0], -pivot[1], -pivot[2]]);
  try {
    f();
  } finally {
    j.pop();
  }
}

/** A flanged floor foot: a round plate with two anchor-bolt holes and a tapered socket for the leg. */
function flangeFoot(s: Shop, x: number, z: number, r: number): void {
  const p = palette();
  s.add(p.ss, cyl(s, r * 2.3, 0.006, "y", 24), [x, 0.003, z]);
  if (s.fine(1)) s.add(p.ss, lathe(s, [[0, 0.006], [r * 1.5, 0.006], [r * 1.05, 0.04], [0, 0.04]], 20), [x, 0, z]);
  if (s.fine(2)) for (const sx of [-1, 1]) s.add(p.enamel, cyl(s, 0.0045, 0.0008, "y", 10), [x + sx * r * 1.8, 0.0062, z]);
}

/**
 * A 40 gal tilting braising pan: a rectangular stainless pan, double-walled round its heated bed, with a pouring
 * lip at the front, carried on trunnions near the front so that the lip stays put as it tips. A handwheel on the
 * console drives a worm gearbox on the right trunnion; the left rides a pillow block on the open tubular stand. The
 * lid is hinged to the pan's back and counterweighted, so it stays where it is put; a swing faucet fills the pan.
 */
const tiltSkillet: Builder = (s, e) => {
  const p = palette();
  const { w, d } = e;
  const back = -d;
  const piv: V3 = [0.6095, 0.76, -0.24];
  const [ix0, ix1, iz1, bed, rim] = [0.09, 1.02, -0.76, 0.68, 0.905];
  const [ox0, ox1, oz1, oBot] = [0.065, 1.045, -0.785, 0.62];
  const frameY = 0.56;
  const [cx0, cx1, cz0, cz1, cy0, cy1] = [1.07, w - 0.008, -0.06, back + 0.09, 0.1, 1.0];
  const lamp = lampMaterial(0xffa030);

  // The stand: square-tube frame under the pan on four round legs, flanged feet bolted to the floor.
  const tr = 0.024;
  for (const x of [0.05, 1.04])
    for (const z of [-0.12, back + 0.08]) {
      s.add(p.ss, cyl(s, tr, frameY - 0.04, "y", 18), [x, (frameY + 0.04) / 2, z]);
      flangeFoot(s, x, z, tr);
    }
  const sq = 0.04;
  s.add(p.ss, rbox(s, 1.04 - 0.05 + sq, sq, sq, 0.004), [(0.05 + 1.04) / 2, frameY, -0.12]);
  s.add(p.ss, rbox(s, 1.04 - 0.05 + sq, sq, sq, 0.004), [(0.05 + 1.04) / 2, frameY, back + 0.08]);
  for (const x of [0.05, 1.04]) s.add(p.ss, rbox(s, sq, sq, -0.12 - back - 0.08, 0.004), [x, frameY, (-0.12 + back + 0.08) / 2]);
  s.add(p.ss, cyl(s, 0.016, -0.12 - back - 0.08, "z", 14), [0.05, 0.2, (-0.12 + back + 0.08) / 2]);
  // The left trunnion's pillow block on its upright.
  s.add(p.ss, rbox(s, sq, piv[1] - frameY - 0.04, sq, 0.004), [0.035, (piv[1] + frameY) / 2 - 0.02, piv[2]]);
  s.add(p.castIron, rbox(s, 0.04, 0.05, 0.11, 0.008), [0.035, piv[1] - 0.012, piv[2]]);
  s.add(p.castIron, cyl(s, 0.042, 0.04, "x", 24), [0.035, piv[1], piv[2]]);
  if (s.fine(2)) s.add(p.brass, cyl(s, 0.004, 0.012, "y", 8), [0.035, piv[1] + 0.046, piv[2]]);
  if (s.fine(2)) for (const sz of [-1, 1]) s.add(p.ssBright, new THREE.CylinderGeometry(0.008, 0.008, 0.008, 6), [0.035, piv[1] - 0.033, piv[2] + sz * 0.042]);

  // The console: a stainless pedestal on its own two legs, the gearbox on its flank.
  const cw = cx1 - cx0;
  for (const z of [cz0 - 0.06, cz1 + 0.06]) {
    s.add(p.ss, cyl(s, tr, cy0 - 0.04, "y", 18), [(cx0 + cx1) / 2, (cy0 + 0.04) / 2, z]);
    flangeFoot(s, (cx0 + cx1) / 2, z, tr);
  }
  s.add(p.ss, rbox(s, cw, cy1 - cy0, cz0 - cz1, 0.006), [(cx0 + cx1) / 2, (cy0 + cy1) / 2, (cz0 + cz1) / 2]);
  louvres(s, [(cx0 + cx1) / 2, 0.24, cz0], cw - 0.05, 0.12, 8);
  s.add(p.charcoal, rbox(s, 0.03, 0.15, 0.16, 0.012), [cx0 - 0.012, piv[1], piv[2]]);
  s.add(p.charcoal, lathe(s, [[0, 0], [0.05, 0], [0.05, 0.008], [0.04, 0.016], [0, 0.016]], 28).rotateZ(Math.PI / 2), [cx0 - 0.027, piv[1], piv[2]]);
  if (s.fine(2)) for (let k = 0; k < 6; k++) {
    const a = (k / 6) * Math.PI * 2;
    s.add(p.ssBright, new THREE.CylinderGeometry(0.004, 0.004, 0.006, 6).rotateZ(Math.PI / 2), [cx0 - 0.03, piv[1] + Math.sin(a) * 0.04, piv[2] + Math.cos(a) * 0.04]);
  }
  // Controls on the console's face: the thermostat, power switch, the heating lamp.
  const fx = (cx0 + cx1) / 2;
  const tk = valveKnob(s, "knob_temp", [fx, 0.925, cz0]);
  dialTicks(s, [fx, 0.925, cz0], 0.031, 11);
  s.add(lamp, cyl(s, 0.006, 0.006, "z", 14), [fx - 0.035, 0.86, cz0 + 0.003]);
  s.add(p.ssBright, torus(s, 0.007, 0.0014, "z", Math.PI * 2, 14, 5), [fx - 0.035, 0.86, cz0 + 0.001]);
  s.add(p.ssBright, cyl(s, 0.011, 0.006, "z", 16), [fx + 0.03, 0.86, cz0 + 0.003]);
  s.add(p.bakelite, cyl(s, 0.003, 0.02, "z", 8), [fx + 0.03, 0.865, cz0 + 0.014], [0.45, 0, 0]);
  dataPlate(s, [cx1, 0.5, (cz0 + cz1) / 2], 0.08, 0.05, [0, Math.PI / 2, 0]);
  // The handwheel: rim, three spokes, a spinner, on its shaft to the worm.
  const wheelAt: V3 = [fx, 0.72, cz0];
  const wheel = s.joint("tilt_wheel", wheelAt);
  const WR = Math.min(0.068, cw / 2 - 0.006);
  wheel.add(p.ssBright, cyl(wheel, 0.008, 0.03, "z", 12), [0, 0, 0.015]);
  wheel.add(p.ssBright, cyl(wheel, 0.017, 0.016, "z", 18), [0, 0, 0.03]);
  wheel.add(p.ssBright, torus(wheel, WR, 0.0065, "z", Math.PI * 2, 48, 10), [0, 0, 0.03]);
  for (let k = 0; k < 3; k++) {
    const a = (k / 3) * Math.PI * 2 + Math.PI / 2;
    wheel.add(p.ssBright, pipe(wheel, [[Math.cos(a) * 0.015, Math.sin(a) * 0.015, 0.03], [Math.cos(a) * WR, Math.sin(a) * WR, 0.03]], 0.0045, 0.005, 8));
  }
  wheel.add(p.bakelite, lathe(wheel, [[0, 0], [0.009, 0], [0.011, 0.012], [0.0095, 0.024], [0, 0.026]], 16).rotateX(Math.PI / 2), [0, -WR, 0.034]);

  // The pan, on its trunnions.
  const pan = s.joint("tilt_pan", piv);
  inFixtureFrame(pan, piv, () => {
    const L = ix1 - ix0;
    // The bed: the heated plate the food sits on; the inner walls, the front one raked forward into the lip.
    pan.add(p.ss, sheet(pan, L, 0.003, -0.14 - iz1, 0.001), [(ix0 + ix1) / 2, bed - 0.0015, (iz1 - 0.14) / 2]);
    for (const x of [ix0 - 0.001, ix1 + 0.001]) pan.add(p.ss, box(0.002, rim - bed + 0.003, -0.1 - iz1), [x, (rim + bed) / 2, (iz1 - 0.1) / 2]);
    pan.add(p.ss, box(L + 0.004, rim - bed + 0.003, 0.002), [(ix0 + ix1) / 2, (rim + bed) / 2, iz1 - 0.001]);
    const rake = Math.atan2(0.035, rim - bed);
    const fl = Math.hypot(0.035, rim - bed);
    pan.add(p.ss, box(L + 0.004, fl + 0.002, 0.002), [(ix0 + ix1) / 2, (rim + bed) / 2, -0.1225], [rake, 0, 0]);
    // Gallon marks etched on the back wall.
    if (pan.fine(2)) for (let k = 1; k <= 4; k++) pan.add(p.charcoal, box(k % 2 ? 0.05 : 0.08, 0.0012, 0.0004), [ix0 + 0.06, bed + k * 0.045, iz1 + 0.0002]);
    // The outer jacket and its insulated floor.
    for (const x of [ox0 + 0.0008, ox1 - 0.0008]) pan.add(p.ss, sheet(pan, 0.0016, rim - oBot, -0.085 - oz1, 0.0006), [x, (rim + oBot) / 2, (oz1 - 0.085) / 2]);
    pan.add(p.ss, sheet(pan, ox1 - ox0, rim - oBot, 0.0016), [(ox0 + ox1) / 2, (rim + oBot) / 2, oz1 + 0.0008]);
    pan.add(p.ss, box(ox1 - ox0, Math.hypot(0.035, rim - oBot), 0.0016), [(ox0 + ox1) / 2, (rim + oBot) / 2, -0.1025], [rake, 0, 0]);
    pan.add(p.ssDull, sheet(pan, ox1 - ox0, 0.0016, -0.115 - oz1), [(ox0 + ox1) / 2, oBot + 0.0008, (oz1 - 0.115) / 2]);
    // The rim: a rolled cap over both walls all round.
    pan.add(p.ss, rbox(pan, 0.03, 0.012, -0.08 - oz1, 0.005), [ox0 + 0.014, rim + 0.003, (oz1 - 0.08) / 2]);
    pan.add(p.ss, rbox(pan, 0.03, 0.012, -0.08 - oz1, 0.005), [ox1 - 0.014, rim + 0.003, (oz1 - 0.08) / 2]);
    pan.add(p.ss, rbox(pan, ox1 - ox0, 0.012, 0.03, 0.005), [(ox0 + ox1) / 2, rim + 0.003, oz1 + 0.014]);
    pan.add(p.ss, rbox(pan, ox1 - ox0, 0.012, 0.032, 0.005), [(ox0 + ox1) / 2, rim + 0.003, -0.094]);
    // The pouring lip: a tapered spout out of the front rim, its sides turned up.
    const spout = taperLip(pan, 0.05, 0.26, 0.12);
    pan.add(p.ss, spout, [(ox0 + ox1) / 2, rim + 0.0035, -0.088], [0, -Math.PI / 2, 0]);
    // Trunnion bosses and stubs.
    pan.add(p.ss, rbox(pan, 0.012, 0.1, 0.1, 0.005), [ox0 - 0.006, piv[1], piv[2]]);
    pan.add(p.ss, rbox(pan, 0.012, 0.1, 0.1, 0.005), [ox1 + 0.006, piv[1], piv[2]]);
    pan.add(p.ssBright, cyl(pan, 0.026, ox0 - 0.012 - 0.022, "x", 20), [(ox0 - 0.012 + 0.022) / 2, piv[1], piv[2]]);
    pan.add(p.ssBright, cyl(pan, 0.026, cx0 - 0.03 - ox1 - 0.012, "x", 20), [(ox1 + 0.012 + cx0 - 0.03) / 2, piv[1], piv[2]]);
    // The lid's hinge brackets on the pan's back.
    for (const x of [0.2, 0.91]) pan.add(p.ss, rbox(pan, 0.012, 0.05, 0.03, 0.003), [x, rim + 0.02, oz1 - 0.012]);
  });

  // The lid, hinged at the pan's back, its counterweight on arms behind the hinge.
  const hinge: V3 = [0.555, rim + 0.031, oz1 - 0.022];
  let lid = pan;
  inFixtureFrame(pan, piv, () => (lid = pan.joint("lid_hinge", hinge)));
  inFixtureFrame(lid, hinge, () => {
    const [lz0, lz1] = [-0.07, oz1 - 0.006];
    const lw = ox1 - ox0 + 0.01;
    const ly = rim + 0.009 + 0.0006;
    lid.add(p.ss, rbox(lid, lw, 0.018, lz0 - lz1, 0.005), [(ox0 + ox1) / 2, ly + 0.009, (lz0 + lz1) / 2]);
    lid.add(p.ss, sheet(lid, lw - 0.12, 0.004, lz0 - lz1 - 0.12, 0.002), [(ox0 + ox1) / 2, ly + 0.019, (lz0 + lz1) / 2]);
    lid.add(p.ss, sheet(lid, lw, 0.03, 0.0016), [(ox0 + ox1) / 2, ly - 0.006, lz0 + 0.0008]);
    // The handle across the front, and the steam vent near the back.
    lid.add(p.ssBright, pipe(lid, [[0.25, ly + 0.018, -0.1], [0.25, ly + 0.048, -0.1], [0.86, ly + 0.048, -0.1], [0.86, ly + 0.018, -0.1]], 0.011, 0.025, 14));
    for (const x of [0.25, 0.86]) lid.add(p.ssBright, cyl(lid, 0.017, 0.004, "y", 16), [x, ly + 0.02, -0.1]);
    lid.add(p.ss, lathe(lid, [[0, 0], [0.04, 0], [0.038, 0.012], [0.03, 0.022], [0, 0.026]], 28), [0.555, ly + 0.018, -0.66]);
    if (lid.fine(2)) for (let k = 0; k < 8; k++) {
      const a = (k / 8) * Math.PI * 2;
      lid.add(p.enamel, box(0.012, 0.002, 0.003), [0.555 + Math.cos(a) * 0.03, ly + 0.034, -0.66 + Math.sin(a) * 0.03], [0, -a, -0.5]);
    }
    // Hinge knuckles and pin; counterweight arms and their weight.
    for (const x of [0.2, 0.91]) {
      lid.add(p.ssBright, cyl(lid, 0.011, 0.04, "x", 16), [x - 0.026, hinge[1], hinge[2]]);
      lid.add(p.ss, rbox(lid, 0.01, 0.03, 0.15, 0.003), [x + 0.02, hinge[1], hinge[2] - 0.05]);
    }
    lid.add(p.ssBright, cyl(lid, 0.004, 0.75, "x", 8), [0.555, hinge[1], hinge[2]]);
    lid.add(p.ss, cyl(lid, 0.028, 0.75, "x", 24), [0.555, hinge[1], hinge[2] - 0.12]);
  });

  // The pan-fill faucet on the console's top: deck plate, two lever valves, the swing spout parked over the console.
  const fz = cz1 + 0.12;
  s.add(p.ssBright, rbox(s, 0.05, 0.008, 0.2, 0.003), [fx, cy1 + 0.004, fz]);
  for (const dz of [-0.07, 0.07]) {
    s.add(p.ssBright, cyl(s, 0.014, 0.03, "y", 16), [fx, cy1 + 0.023, fz + dz]);
    s.add(p.ssBright, pipe(s, [[fx, cy1 + 0.038, fz + dz], [fx, cy1 + 0.044, fz + dz], [fx + 0.008, cy1 + 0.046, fz + dz + Math.sign(dz) * 0.05]], 0.004, 0.008, 8));
  }
  s.add(p.ssBright, cyl(s, 0.016, 0.04, "y", 18), [fx, cy1 + 0.028, fz]);
  const sp = s.joint("faucet_spout", [fx, cy1 + 0.048, fz]);
  sp.add(p.ssBright, pipe(sp, [[0, 0, 0], [0, 0.012, 0], [0, 0.012, 0.16], [0, -0.008, 0.18]], 0.0075, 0.02, 12));
  // Supply lines from the floor, and the conduit for the elements.
  if (s.fine(1)) {
    s.add(p.ssDull, rbox(s, 0.1, 0.1, 0.06, 0.004), [fx, 0.4, cz1 - 0.03]);
  }

  return (b) => {
    const tilt = (v: number) => {
      const k = Math.min(1, Math.max(0, v));
      turnJoint(b.joints.get("tilt_pan"), "x", k * 1.4835);
      turnJoint(b.joints.get("tilt_wheel"), "z", -k * 12 * Math.PI);
    };
    const lidTo = (v: number) => turnJoint(b.joints.get("lid_hinge"), "x", -Math.min(1, Math.max(0, v)) * 1.4835);
    return {
      glows: [lamp],
      open: lidTo,
      set: (name, v) => {
        if (name === "tilt") tilt(v);
        else if (name === "lid") lidTo(v);
        else if (name === "heat") {
          lamp.emissiveIntensity = v > 0.05 ? 2 : 0;
          turnJoint(b.joints.get(tk), "z", -v * 1.5 * Math.PI);
        }
      },
    };
  };
};

/** A pouring lip in plan: `len` long out of the rim from a root `w0` wide to a tip `w1` wide, its sides turned up. */
function taperLip(s: Shop, len: number, w0: number, w1: number): THREE.BufferGeometry {
  const shp = new THREE.Shape([new THREE.Vector2(0, -w0 / 2), new THREE.Vector2(len, -w1 / 2), new THREE.Vector2(len + 0.004, 0), new THREE.Vector2(len, w1 / 2), new THREE.Vector2(0, w0 / 2)]);
  const g = plate(s, shp, 0.003, 0.0008, 6);
  g.rotateX(-Math.PI / 2);
  return g;
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Steam kettle                                                                                                   */
/* ------------------------------------------------------------------------------------------------------------ */

/** A dial gauge facing +z at its base: chrome case, white face with a scale and a needle, a glass. `r` its radius. */
function gauge(s: Shop, at: V3, rot: V3, r: number, needle = 0.3): void {
  const p = palette();
  s.at(at, rot, () => {
    s.add(p.ssBright, cyl(s, r, 0.022, "z", 28), [0, 0, 0.011]);
    s.add(p.hdpe, cyl(s, r * 0.88, 0.002, "z", 28), [0, 0, 0.0215]);
    if (!s.fine(1)) return;
    s.add(p.ssBright, torus(s, r * 0.93, r * 0.07, "z", Math.PI * 2, 28, 6), [0, 0, 0.022]);
    if (s.fine(2)) {
      for (let k = 0; k <= 10; k++) {
        const a = Math.PI * 1.25 - (k / 10) * Math.PI * 1.5;
        s.add(p.enamel, box(0.0009, k % 5 ? r * 0.1 : r * 0.18, 0.0004), [Math.cos(a) * r * 0.72, Math.sin(a) * r * 0.72, 0.0228], [0, 0, a - Math.PI / 2]);
      }
      // The green band of normal working pressure is printed dark; the needle and its hub.
      s.add(p.enamel, new THREE.RingGeometry(r * 0.5, r * 0.56, 16, 1, Math.PI * 0.55, Math.PI * 0.3), [0, 0, 0.0227]);
      const a = Math.PI * 1.25 - needle * Math.PI * 1.5;
      s.add(p.enamel, box(0.0014, r * 0.7, 0.0006), [Math.cos(a) * r * 0.3, Math.sin(a) * r * 0.3, 0.0232], [0, 0, a - Math.PI / 2]);
      s.add(p.enamel, cyl(s, r * 0.08, 0.002, "z", 12), [0, 0, 0.0235]);
    }
    s.add(p.glass, cyl(s, r * 0.9, 0.001, "z", 28), [0, 0, 0.0255]);
  });
}

/**
 * A 40 gal floor kettle: a hemispherical bowl with a short straight wall above, jacketed over its lower two-thirds
 * (the steam space between, raised by elements in the jacket's floor), on three tubular legs. A 2 in tangent
 * draw-off leaves the bowl's lowest point and runs forward to its valve; the jacket carries the pressure gauge and
 * the safety valve; the lid is hinged at the back on a spring; the controls are in a console off the right leg.
 */
const kettle: Builder = (s, e) => {
  const p = palette();
  const { w } = e;
  const [cx, cz] = [w / 2, -0.482];
  const R = 0.33; // the bowl's inside
  const yc = 0.72; // the hemisphere's centre
  const rimY = 0.935;
  const Ro = 0.37; // the jacket's outside
  const jTop = 0.77;
  const lamp = lampMaterial(0xffa030);
  const n = s.fine(1) ? 64 : 36;
  const arc = (r: number, a0: number, a1: number, k: number, down: boolean) => {
    const pts: [number, number][] = [];
    for (let i = 0; i <= k; i++) {
      const a = a0 + ((a1 - a0) * i) / k;
      pts.push([r * Math.sin(a), yc - r * Math.cos(a)]);
    }
    return down ? pts.reverse() : pts;
  };
  const k = s.fine(1) ? s.seg(18, 6) : 4;
  // The bowl's inside, seen from within (traced rim to bottom so its faces look in).
  s.add(p.ss, lathe(s, [[R, rimY - 0.006], [R, yc], ...arc(R, Math.PI / 2, 0, k, false).slice(1)], n), [cx, 0, cz]);
  // The jacket, its band welded to the bowl's wall, the wall above it, the rolled rim.
  s.add(p.ss, lathe(s, [...arc(Ro, 0, Math.PI / 2, k, false), [Ro, jTop - 0.012], [R + 0.03, jTop], [R + 0.003, jTop + 0.004]], n), [cx, 0, cz]);
  s.add(p.ss, lathe(s, [[R + 0.0025, jTop], [R + 0.0025, rimY - 0.006]], n), [cx, 0, cz]);
  if (s.fine(1)) {
    s.add(p.ss, torus(s, R + 0.009, 0.0095, "y", Math.PI * 2, n, 10), [cx, rimY - 0.0055, cz]);
    s.add(p.ss, torus(s, Ro - 0.004, 0.004, "y", Math.PI * 2, n, 6), [cx, jTop - 0.012, cz]);
  }

  // Three legs: two forward, one behind, each welded to a pad on the jacket, on flanged feet.
  const legR = 0.025;
  const legs: [number, number][] = [[cx - 0.31, cz + 0.18], [cx + 0.31, cz + 0.18], [cx, cz - 0.36]];
  for (const [lx, lz] of legs) {
    const rr = Math.hypot(lx - cx, lz - cz);
    const top = yc - Math.sqrt(Math.max(0, Ro * Ro - rr * rr)) + 0.003;
    s.add(p.ss, cyl(s, legR, top - 0.04, "y", 20), [lx, (top + 0.04) / 2, lz]);
    flangeFoot(s, lx, lz, legR);
    if (s.fine(1)) s.add(p.ss, lathe(s, [[0, 0], [legR * 1.05, 0], [legR * 1.35, 0.045], [0, 0.045]], 20), [lx, top - 0.045, lz]);
  }
  // The draw-off: out of the bowl's lowest point, forward under the jacket to the valve, its outlet turned down.
  const dy = 0.31;
  const vz = -0.1;
  s.add(p.ss, pipe(s, [[cx, yc - Ro + 0.01, cz], [cx, dy, cz], [cx, dy, vz - 0.04]], 0.03, 0.05, 16));
  if (s.fine(1)) s.add(p.ss, cyl(s, 0.036, 0.012, "z", 20), [cx, dy, vz - 0.042]);
  s.add(p.ssBright, rbox(s, 0.075, 0.085, 0.075, 0.012), [cx, dy, vz]);
  s.add(p.ssBright, cyl(s, 0.026, 0.06, "y", 20), [cx, dy - 0.07, vz]);
  if (s.fine(1)) s.add(p.ssBright, lathe(s, [[0.026, 0], [0.031, 0], [0.031, 0.008], [0.026, 0.008]], 20), [cx, dy - 0.106, vz]);
  s.add(p.ssBright, cyl(s, 0.012, 0.03, "y", 14), [cx, dy + 0.055, vz]);
  // The valve's T-handle, swung forward.
  s.add(p.ssBright, pipe(s, [[cx, dy + 0.07, vz], [cx, dy + 0.075, vz + 0.03], [cx, dy + 0.06, vz + 0.075]], 0.008, 0.02, 10));
  s.add(p.bakelite, cyl(s, 0.012, 0.11, "x", 16), [cx, dy + 0.06, vz + 0.078]);

  // The pressure gauge on a stub off the jacket's left flank, the safety valve on the back, the vacuum breaker.
  const ga = -Math.PI * 0.72;
  const gr = Ro + 0.02;
  s.add(p.ss, cyl(s, 0.008, 0.05, "z", 10), [cx + Math.sin(ga) * (Ro + 0.005), 0.66, cz + Math.cos(ga) * (Ro + 0.005)], [0, ga, 0]);
  gauge(s, [cx + Math.sin(ga) * gr, 0.66, cz + Math.cos(ga) * gr], [0, ga, 0], 0.045, 0.42);
  const sa = Math.PI * 0.8;
  const [sx, sz] = [cx + Math.sin(sa) * (Ro - 0.02), cz + Math.cos(sa) * (Ro - 0.02)];
  s.add(p.brass, cyl(s, 0.014, 0.05, "y", 16), [sx, jTop + 0.0, sz]);
  s.add(p.brass, lathe(s, [[0, 0], [0.022, 0], [0.022, 0.03], [0.012, 0.05], [0, 0.052]], 18), [sx, jTop + 0.02, sz]);
  if (s.fine(1)) {
    s.add(p.brass, torus(s, 0.012, 0.0022, "x", Math.PI * 2, 14, 5), [sx, jTop + 0.085, sz]);
    s.add(p.brass, cyl(s, 0.0025, 0.02, "y", 6), [sx, jTop + 0.066, sz]);
  }
  const va = -Math.PI * 0.82;
  s.add(p.ssBright, cyl(s, 0.01, 0.045, "y", 12), [cx + Math.sin(va) * (Ro - 0.02), jTop + 0.02, cz + Math.cos(va) * (Ro - 0.02)]);
  // The elements' wiring box under the jacket.
  if (s.fine(1)) s.add(p.ssDull, rbox(s, 0.2, 0.08, 0.14, 0.005), [cx, yc - Ro - 0.02, cz - 0.16]);

  // The console: a stainless box off the right front leg, the thermostat, the switch and lamps.
  const [k0, k1, kz0, kz1, ky0, ky1] = [0.745, w - 0.01, -0.045, -0.2, 0.46, 0.86];
  const kxm = (k0 + k1) / 2;
  s.add(p.ss, rbox(s, k1 - k0, ky1 - ky0, kz0 - kz1, 0.005), [kxm, (ky0 + ky1) / 2, (kz0 + kz1) / 2]);
  for (const y of [ky0 + 0.06, ky1 - 0.06]) s.add(p.ss, rbox(s, 0.03, 0.025, legs[1]![1] - kz1 + 0.02, 0.003), [legs[1]![0] + 0.01, y, (legs[1]![1] + kz1) / 2]);
  const tk = valveKnob(s, "knob_temp", [kxm, ky1 - 0.08, kz0]);
  dialTicks(s, [kxm, ky1 - 0.08, kz0], 0.031, 11);
  for (const [i, lx] of [[0, kxm - 0.04], [1, kxm + 0.04]] as const) {
    s.add(i ? p.ledGreen : lamp, cyl(s, 0.006, 0.006, "z", 14), [lx, ky1 - 0.165, kz0 + 0.003]);
    if (s.fine(1)) s.add(p.ssBright, torus(s, 0.007, 0.0014, "z", Math.PI * 2, 14, 5), [lx, ky1 - 0.165, kz0 + 0.001]);
  }
  s.add(p.ssBright, cyl(s, 0.012, 0.006, "z", 16), [kxm, ky1 - 0.23, kz0 + 0.003]);
  s.add(p.bakelite, cyl(s, 0.003, 0.02, "z", 8), [kxm, ky1 - 0.225, kz0 + 0.014], [0.45, 0, 0]);
  dataPlate(s, [kxm, ky0 + 0.07, kz0], 0.08, 0.05);
  // Its supply comes up from the floor in a conduit.
  if (s.fine(1)) {
    s.add(p.ssDull, cyl(s, 0.012, ky0, "y", 10), [kxm, ky0 / 2, (kz0 + kz1) / 2]);
    s.add(p.ssDull, cyl(s, 0.024, 0.006, "y", 14), [kxm, 0.003, (kz0 + kz1) / 2]);
  }

  // The contents: a level surface whose radius follows the bowl, and its simmer.
  const contents = s.joint("contents", [cx, 0.6, cz]);
  contents.add(p.oil, new THREE.CircleGeometry(1, s.seg(64, 16)).rotateX(-Math.PI / 2));

  // The lid: a shallow dome hinged on a sprung bracket behind the rim, a loop handle at its front.
  const hinge: V3 = [cx, rimY + 0.012, cz - R - 0.06];
  const lid = s.joint("kettle_lid", hinge);
  inFixtureFrame(lid, hinge, () => {
    const lr = R + 0.022;
    const ly = rimY + 0.004;
    lid.add(p.ss, lathe(lid, [[lr, -0.008], [lr + 0.002, 0], [lr * 0.75, 0.022], [lr * 0.4, 0.034], [0, 0.038]], n), [cx, ly, cz]);
    if (lid.fine(1)) lid.add(p.ss, lathe(lid, [[0, 0], [lr - 0.004, -0.006], [lr - 0.004, -0.0065], [0, -0.0005]], n), [cx, ly - 0.002, cz]);
    lid.add(p.ssBright, pipe(lid, [[cx - 0.07, ly + 0.012, cz + lr - 0.03], [cx - 0.07, ly + 0.04, cz + lr + 0.01], [cx + 0.07, ly + 0.04, cz + lr + 0.01], [cx + 0.07, ly + 0.012, cz + lr - 0.03]], 0.008, 0.02, 12));
    for (const sxx of [-1, 1]) lid.add(p.ss, rbox(lid, 0.014, 0.012, 0.09, 0.003), [cx + sxx * 0.06, ly + 0.012, cz - lr + 0.02], [0.2, 0, 0]);
    lid.add(p.ssBright, cyl(lid, 0.01, 0.15, "x", 14), [cx, hinge[1], hinge[2]]);
  });
  // The hinge bracket on the rim's back and its counterbalance spring.
  s.add(p.ss, rbox(s, 0.2, 0.05, 0.05, 0.005), [cx, rimY - 0.012, cz - R - 0.035]);
  if (s.fine(1)) s.add(p.ssBright, cyl(s, 0.018, 0.2, "x", 18), [cx, rimY + 0.012, cz - R - 0.06]);
  if (s.fine(2)) for (let i = 0; i < 9; i++) s.add(p.ssDull, torus(s, 0.0205, 0.0025, "x", Math.PI * 2, 14, 5), [cx - 0.09 + i * 0.0225, rimY + 0.012, cz - R - 0.06]);

  return (b) => {
    const bub = new Bubbles(p.hdpe, s.fine(2) ? 50 : 20, [cx, 0.6, cz], 0.5, 0.5, 0.005, true);
    b.group.add(bub.mesh);
    let level = 0.6;
    let surf = 0.7;
    let rad = R;
    const setLevel = (v: number) => {
      level = Math.min(1, Math.max(0, v));
      const bottom = yc - R;
      surf = bottom + 0.004 + level * (0.895 - bottom);
      rad = surf < yc ? Math.sqrt(Math.max(0, R * R - (yc - surf) * (yc - surf))) : R;
      const j = b.joints.get("contents");
      if (j) {
        j.visible = level > 0.005;
        j.position.set(cx, surf, cz);
        j.scale.set(rad, 1, rad);
      }
      bub.at = [cx, surf + 0.0005, cz];
      bub.w = bub.d = rad * 1.6;
    };
    setLevel(level);
    return {
      glows: [lamp],
      open: (f) => turnJoint(b.joints.get("kettle_lid"), "x", -Math.min(1, Math.max(0, f)) * 1.396),
      set: (name, v) => {
        if (name === "level") setLevel(v);
        else if (name === "lid") turnJoint(b.joints.get("kettle_lid"), "x", -Math.min(1, Math.max(0, v)) * 1.396);
        else if (name === "heat") {
          lamp.emissiveIntensity = v > 0.05 ? 2 : 0;
          turnJoint(b.joints.get(tk), "z", -v * 1.5 * Math.PI);
        }
      },
      tick: (t, _dt, load) => bub.tick(t, level > 0.02 ? 0.25 + 0.5 * load : 0),
    };
  };
};

/* ------------------------------------------------------------------------------------------------------------ */
/* Strip warmer                                                                                                   */
/* ------------------------------------------------------------------------------------------------------------ */

/**
 * A 48 in infrared strip warmer, hung under a shelf over the pass: an extruded aluminium channel, finned on top, with
 * a polished reflector inside throwing a calrod element's heat down; end caps carry the element's ceramic terminals,
 * the toggle switch and the cord; two slotted brackets on top bolt it up.
 */
const stripWarmer: Builder = (s, e) => {
  const p = palette();
  const { w, d, h } = e;
  const bodyH = 0.05;
  const cap = 0.018;
  const [x0, x1] = [0.012 + cap, w - 0.012 - cap];
  const glow = glowMaterial(0x3a3633);
  // The extrusion's section, traced in (depth back from the front, height): a C open below, fins along its top.
  const sec = new THREE.Shape();
  const t = 0.0028;
  const lip = 0.02;
  sec.moveTo(d - lip, 0);
  sec.lineTo(d, 0);
  sec.lineTo(d, bodyH - 0.004);
  sec.quadraticCurveTo(d, bodyH, d - 0.004, bodyH);
  const fins = s.fine(1) ? 6 : 0;
  for (let k = fins - 1; k >= 0; k--) {
    const u = 0.025 + ((d - 0.05) * (k + 0.5)) / fins;
    sec.lineTo(u + 0.004, bodyH);
    sec.lineTo(u + 0.004, bodyH - 0.0025);
    sec.lineTo(u - 0.004, bodyH - 0.0025);
    sec.lineTo(u - 0.004, bodyH);
  }
  sec.lineTo(0.004, bodyH);
  sec.quadraticCurveTo(0, bodyH, 0, bodyH - 0.004);
  sec.lineTo(0, 0);
  sec.lineTo(lip, 0);
  sec.lineTo(lip, t);
  sec.lineTo(t, t);
  sec.lineTo(t, bodyH - t);
  sec.lineTo(d - t, bodyH - t);
  sec.lineTo(d - t, t);
  sec.lineTo(d - lip, t);
  sec.lineTo(d - lip, 0);
  const g = plate(s, sec, x1 - x0, 0, 4);
  g.rotateY(Math.PI / 2);
  s.add(p.alu, g, [x0, 0, 0]);
  // The reflector: a polished trough arched over the element.
  const ref = new THREE.Shape();
  const [ra, rb] = [d / 2 - 0.006, bodyH - t - 0.005];
  const m = s.seg(20, 8);
  for (let i = 0; i <= m; i++) {
    const a = Math.PI * (i / m);
    const [u, v] = [d / 2 - Math.cos(a) * ra, 0.012 + Math.sin(a) * (rb - 0.012)];
    if (i === 0) ref.moveTo(u, v);
    else ref.lineTo(u, v);
  }
  for (let i = m; i >= 0; i--) {
    const a = Math.PI * (i / m);
    ref.lineTo(d / 2 - Math.cos(a) * (ra - 0.0012), 0.012 + Math.sin(a) * (rb - 0.0132));
  }
  const rg = plate(s, ref, x1 - x0 - 0.002, 0, 4);
  rg.rotateY(Math.PI / 2);
  s.add(p.ssBright, rg, [x0 + 0.001, 0, 0]);
  // The element: a calrod tube end to end, its terminals in ceramic bushings, clipped to the reflector twice.
  const [ey, ez] = [0.02, -d / 2];
  s.add(glow, cyl(s, 0.0055, x1 - x0 + 0.01, "x", 14), [w / 2, ey, ez]);
  for (const x of [x0 + 0.004, x1 - 0.004]) s.add(p.ceramic, cyl(s, 0.0095, 0.014, "x", 16), [x, ey, ez]);
  if (s.fine(1))
    for (const f of [0.33, 0.67]) {
      const x = x0 + (x1 - x0) * f;
      s.add(p.ssBright, pipe(s, [[x, rb + 0.008, ez], [x, ey + 0.009, ez], [x, ey, ez + 0.008]], 0.0012, 0.004, 6));
    }
  // End caps; the toggle switch on the right one, the cord's strain relief through the left.
  for (const x of [0.012 + cap / 2, w - 0.012 - cap / 2]) s.add(p.ss, rbox(s, cap, bodyH + 0.002, d, 0.004), [x, bodyH / 2, -d / 2]);
  s.add(p.ss, rbox(s, 0.003, 0.03, 0.05, 0.002), [w - 0.012 + 0.0015, bodyH / 2, -d / 2]);
  s.add(p.ssBright, cyl(s, 0.008, 0.004, "x", 16), [w - 0.008, bodyH / 2, -d / 2]);
  s.add(p.ssBright, cyl(s, 0.0025, 0.006, "x", 8), [w - 0.005, bodyH / 2, -d / 2]);
  const sw = s.joint("toggle", [w - 0.003, bodyH / 2, -d / 2]);
  sw.add(p.ssBright, cyl(sw, 0.0016, 0.01, "x", 8, 0.0012), [0.004, 0, 0], [0, 0, 0.35]);
  s.add(p.rubber, lathe(s, [[0, 0], [0.008, 0], [0.006, 0.01], [0.0042, 0.012], [0, 0.012]], 14).rotateZ(Math.PI / 2), [0.012, bodyH / 2, -d * 0.7]);
  s.add(p.rubber, cyl(s, 0.004, 0.012, "x", 10), [0.006 - 0.006, bodyH / 2, -d * 0.7]);
  dataPlate(s, [w * 0.7, bodyH * 0.5, -d], 0.06, 0.03, [0, Math.PI, 0]);
  // Mounting brackets: hat sections on top, slotted for the bolts.
  for (const x of [0.22, w - 0.22]) {
    s.add(p.ss, sheet(s, 0.05, 0.002, d - 0.02, 0.0006), [x, h - 0.001, -d / 2]);
    for (const z of [-0.016, -d + 0.016]) s.add(p.ss, sheet(s, 0.05, h - bodyH, 0.002, 0.0006), [x, (h + bodyH) / 2, z]);
    if (s.fine(1)) for (const z of [-0.045, -d + 0.045]) s.add(p.enamel, box(0.008, 0.0004, 0.022), [x, h + 0.0001, z]);
  }

  return (b) => ({
    glows: [glow],
    set: (name, v) => {
      if (name === "heat") {
        setGlow(glow, v);
        turnJoint(b.joints.get("toggle"), "y", 0);
        turnJoint(b.joints.get("toggle"), "z", v > 0.05 ? -0.7 : 0);
      }
    },
  });
};

register({
  "NSF-HEATLAMP-STRIP-48": stripWarmer,
  "NSF-KETTLE-40GAL": kettle,
  "NSF-TILT-SKILLET-40GAL": tiltSkillet,
  "NSF-INDUCTION-HOB-2Z": inductionHob,
  "NSF-INDUCTION-SUITE-4Z": inductionSuite,
  "NSF-SALAMANDER-36": salamander,
  "NSF-FRYER-TWIN-50LB": fryer,
  "NSF-CHARBROILER-36": charbroiler,
  "NSF-GRIDDLE-24": griddle(),
  "NSF-GRIDDLE-36": griddle(),
  "NSF-GRIDDLE-48": griddle(),
});
