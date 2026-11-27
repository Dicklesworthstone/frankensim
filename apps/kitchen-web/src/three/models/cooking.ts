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
import { mergeGeometries, toCreasedNormals } from "three/examples/jsm/utils/BufferGeometryUtils.js";
import { FlameSet, linePorts, type Port } from "./flames";
import { valveAngle } from "./hotline";
import { register, type Builder, type Ext, type Hooks } from "./index";
import { palette } from "./materials";
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
  roundRect,
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
  // Close up the leg starts inside its nut, so that its end does not lie flush with the nut's underside.
  const y0 = s.fine(2) ? fh + 0.004 : fh;
  s.add(p.ss, cyl(s, r, h - y0, "y", 20), [x, y0 + (h - y0) / 2, z]);
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

/**
 * A 6 in adjustable floor leg: a stainless tube from a welded socket under the cabinet, a hex adjuster and a
 * polished bullet foot that screws up and down in it.
 */
function bulletLeg(s: Shop, x: number, z: number, h: number, r = 0.019): void {
  const p = palette();
  const fh = Math.min(0.055, h * 0.36);
  s.add(p.ss, cyl(s, r, h - fh + 0.002, "y", 20), [x, fh + (h - fh) / 2 - 0.001, z]);
  s.add(p.ssBright, lathe(s, [[0, 0], [r * 0.5, 0], [r * 0.78, 0.004], [r * 1.0, fh * 0.32], [r * 1.12, fh * 0.78], [r * 1.1, fh], [0, fh]], 24), [x, 0, z]);
  if (s.fine(1)) s.add(p.ssBright, new THREE.CylinderGeometry(r * 1.22, r * 1.22, 0.009, 6), [x, fh + 0.0045, z]);
  if (s.fine(2)) {
    // The thread showing between the adjuster and the foot, and the socket's collar under the cabinet.
    for (let k = 0; k < 4; k++) s.add(p.ssDull, torus(s, r * 0.82, 0.0009, "y", Math.PI * 2, 18, 4), [x, fh - 0.002 - k * 0.0022, z]);
    s.add(p.ss, cyl(s, r * 1.3, 0.012, "y", 22), [x, h - 0.006, z]);
  }
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
  private sc = new THREE.Vector3();
  private pos = new THREE.Vector3();
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
    const { sc, pos } = this;
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
/* Formed sheet: smooth sweeps, drawn pans, rolled beads                                                          */
/* ------------------------------------------------------------------------------------------------------------ */

type P2 = readonly [number, number];
const clamp01 = (v: number) => Math.min(1, Math.max(0, v));
const smooth = (a: number, b: number, x: number) => {
  const t = clamp01((x - a) / (b - a));
  return t * t * (3 - 2 * t);
};

/**
 * A smooth surface through a grid of points, rows `i` along u and columns `j` along v. Its face is du × dv (the
 * back of it with `flip`); rows or columns that close on themselves share their seam, so the shading runs on round
 * it. UVs run 0..1 along both.
 */
function gridMesh(rows: readonly (readonly THREE.Vector3[])[], closeU = false, closeV = false, flip = false): THREE.BufferGeometry {
  const nu = rows.length;
  const nv = rows[0]!.length;
  const pos = new Float32Array(nu * nv * 3);
  const uv = new Float32Array(nu * nv * 2);
  for (let i = 0; i < nu; i++)
    for (let j = 0; j < nv; j++) {
      rows[i]![j]!.toArray(pos, (i * nv + j) * 3);
      uv[(i * nv + j) * 2] = i / Math.max(1, nu - 1);
      uv[(i * nv + j) * 2 + 1] = j / Math.max(1, nv - 1);
    }
  const idx: number[] = [];
  for (let i = 0; i < (closeU ? nu : nu - 1); i++)
    for (let j = 0; j < (closeV ? nv : nv - 1); j++) {
      const a = i * nv + j;
      const b = ((i + 1) % nu) * nv + j;
      const c = i * nv + ((j + 1) % nv);
      const d = ((i + 1) % nu) * nv + ((j + 1) % nv);
      if (flip) idx.push(a, c, b, b, c, d);
      else idx.push(a, b, c, b, d, c);
    }
  const g = new THREE.BufferGeometry();
  g.setAttribute("position", new THREE.BufferAttribute(pos, 3));
  g.setAttribute("uv", new THREE.BufferAttribute(uv, 2));
  g.setIndex(idx);
  g.computeVertexNormals();
  return g;
}

/**
 * A section swept straight along x from x0 to x1, its points [z, y]. Drawn with z to the right and y up, the face
 * is on the right of the way the section is traced (`flip` for the left). A point given twice makes a sharp fold.
 */
function sweepX(prof: readonly P2[], x0: number, x1: number, flip = false): THREE.BufferGeometry {
  return gridMesh([x0, x1].map((x) => prof.map(([z, y]) => new THREE.Vector3(x, y, z))), false, false, flip);
}

/**
 * A rounded rectangle in plan, x0..x1 by z0..z1 (z1 the front), its corners of radius `r` in `n` steps each, traced
 * from the front-left along the front: a section swept round it from the rim downward then faces inward.
 */
function loopRR(x0: number, x1: number, z0: number, z1: number, r: number, n: number): P2[] {
  const rr = Math.max(1e-4, Math.min(r, (x1 - x0) / 2 - 1e-4, (z1 - z0) / 2 - 1e-4));
  const pts: P2[] = [];
  const corners: [number, number, number][] = [
    [x1 - rr, z1 - rr, Math.PI / 2],
    [x1 - rr, z0 + rr, 0],
    [x0 + rr, z0 + rr, -Math.PI / 2],
    [x0 + rr, z1 - rr, -Math.PI],
  ];
  for (const [cx, cz, a0] of corners)
    for (let k = 0; k <= n; k++) {
      const a = a0 - (k / n) * (Math.PI / 2);
      pts.push([cx + Math.cos(a) * rr, cz + Math.sin(a) * rr]);
    }
  return pts;
}

/** Inward unit normals of a convex plan loop, at each of its points. */
function loopNormals(loop: readonly P2[]): P2[] {
  const n = loop.length;
  let cx = 0;
  let cz = 0;
  for (const [x, z] of loop) {
    cx += x / n;
    cz += z / n;
  }
  return loop.map(([x, z], i) => {
    const [xa, za] = loop[(i + n - 1) % n]!;
    const [xb, zb] = loop[(i + 1) % n]!;
    let nx = -(zb - za);
    let nz = xb - xa;
    const l = Math.hypot(nx, nz) || 1;
    nx /= l;
    nz /= l;
    if (nx * (cx - x) + nz * (cz - z) < 0) [nx, nz] = [-nx, -nz];
    return [nx, nz] as P2;
  });
}

/**
 * A section swept round a closed plan loop: `prof` is [inset, y] points, the inset measured inward from the loop.
 * Traced from the rim downward it faces in (`flip` for out); `closeV` closes the section on itself (a bead).
 */
function sweepLoop(loop: readonly P2[], prof: readonly P2[], flip = false, closeV = false): THREE.BufferGeometry {
  const nrm = loopNormals(loop);
  const rows = loop.map(([x, z], i) => prof.map(([ins, y]) => new THREE.Vector3(x + nrm[i]![0] * ins, y, z + nrm[i]![1] * ins)));
  return gridMesh(rows, true, closeV, flip);
}

/** A flat fan filling a closed plan loop at height y (inset by `ins`), facing up, or down with `down`. */
function loopCap(loop: readonly P2[], y: number, ins = 0, down = false): THREE.BufferGeometry {
  const nrm = loopNormals(loop);
  const n = loop.length;
  const pos = new Float32Array((n + 1) * 3);
  let cx = 0;
  let cz = 0;
  loop.forEach(([x, z], i) => {
    const px = x + nrm[i]![0] * ins;
    const pz = z + nrm[i]![1] * ins;
    pos.set([px, y, pz], i * 3);
    cx += px / n;
    cz += pz / n;
  });
  pos.set([cx, y, cz], n * 3);
  // Planar UVs in metres, so a map laid on the cap keeps its scale.
  const uv = new Float32Array((n + 1) * 2);
  for (let i = 0; i <= n; i++) uv.set([pos[i * 3]!, -pos[i * 3 + 2]!], i * 2);
  const idx: number[] = [];
  // The loop runs anticlockwise seen from above, so (centre, i, i+1) faces up.
  for (let i = 0; i < n; i++) down ? idx.push(n, (i + 1) % n, i) : idx.push(n, i, (i + 1) % n);
  const g = new THREE.BufferGeometry();
  g.setAttribute("position", new THREE.BufferAttribute(pos, 3));
  g.setAttribute("uv", new THREE.BufferAttribute(uv, 2));
  g.setIndex(idx);
  g.computeVertexNormals();
  return g;
}

/** Points round a circle in a section, [inset, y], starting at angle `a0` and closing on itself. */
function ring2(ci: number, cy: number, r: number, n: number, a0 = 0): P2[] {
  const pts: P2[] = [];
  for (let k = 0; k < n; k++) {
    const a = a0 + (k / n) * Math.PI * 2;
    pts.push([ci + Math.cos(a) * r, cy + Math.sin(a) * r]);
  }
  return pts;
}

/** A quarter-round (cove or radius) from a to b about centre c, `n` steps, [u, v] points. */
function arc2(c: P2, r: number, a0: number, a1: number, n: number): P2[] {
  const pts: P2[] = [];
  for (let k = 0; k <= n; k++) {
    const a = a0 + ((a1 - a0) * k) / n;
    pts.push([c[0] + Math.cos(a) * r, c[1] + Math.sin(a) * r]);
  }
  return pts;
}

/**
 * A pan drawn from one sheet, in its own frame: plan `w` by `d` centred on the origin, the rim at y = 0 and the
 * floor `h` below, plan corners of radius `rc`, the floor coved by `rb` inside, a rolled bead of radius `bead` round
 * the rim (0 for a plain turned edge). Drawer pans, water pans, drip pans: one merged solid.
 */
function drawnPan(s: Shop, w: number, d: number, h: number, rc: number, rb: number, t = 0.0012, bead = 0.003): THREE.BufferGeometry {
  const cs = s.fine(2) ? 4 : 2;
  const loop = loopRR(-w / 2, w / 2, -d / 2, d / 2, rc, cs);
  const cove = s.fine(2) ? 4 : 2;
  const rbi = Math.max(0.001, rb - t);
  // Inside: down the wall, round the cove, out onto the floor; the floor's flat middle.
  const inner: P2[] = [[t, 0], [t, -h + t + rbi], ...arc2([t + rbi, -h + t + rbi], rbi, Math.PI, Math.PI * 1.5, cove).slice(1)];
  const outer: P2[] = [[0, 0], [0, -h + rb], ...arc2([rb, -h + rb], rb, Math.PI, Math.PI * 1.5, cove).slice(1)];
  const parts = [sweepLoop(loop, inner), loopCap(loop, -h + t, t + rbi), sweepLoop(loop, outer, true), loopCap(loop, -h, rb, true)];
  if (bead > 0) parts.push(sweepLoop(loop, ring2(-bead + t, -bead * 0.9, bead, s.seg(10, 6), 0), true, true));
  else parts.push(sweepLoop(loop, [[0, 0], [t, 0]], true));
  return mergeGeometries(parts, false)!;
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Finishes only these fixtures use, each made once                                                               */
/* ------------------------------------------------------------------------------------------------------------ */

function rand(seed: number): () => number {
  let s = seed >>> 0;
  return () => {
    s = (s * 1664525 + 1013904223) >>> 0;
    return s / 4294967296;
  };
}
function canvas2d(w: number, h: number): [HTMLCanvasElement, CanvasRenderingContext2D] {
  const c = document.createElement("canvas");
  c.width = w;
  c.height = h;
  return [c, c.getContext("2d")!];
}
function canvasTex(c: HTMLCanvasElement, srgb: boolean): THREE.CanvasTexture {
  const t = new THREE.CanvasTexture(c);
  t.colorSpace = srgb ? THREE.SRGBColorSpace : THREE.NoColorSpace;
  t.anisotropy = 8;
  t.needsUpdate = true;
  return t;
}
/** Set a part's UVs from its own coordinates. */
function mapUV(g: THREE.BufferGeometry, f: (x: number, y: number, z: number, nx: number, ny: number, nz: number) => [number, number]): THREE.BufferGeometry {
  if (!g.attributes.normal) g.computeVertexNormals();
  const p = g.attributes.position as THREE.BufferAttribute;
  const n = g.attributes.normal as THREE.BufferAttribute;
  const uv = new Float32Array(p.count * 2);
  for (let i = 0; i < p.count; i++) {
    const [u, v] = f(p.getX(i), p.getY(i), p.getZ(i), n.getX(i), n.getY(i), n.getZ(i));
    uv[i * 2] = u;
    uv[i * 2 + 1] = v;
  }
  g.setAttribute("uv", new THREE.BufferAttribute(uv, 2));
  return g;
}

/**
 * A seasoned griddle plate of `n` 12 in sections, mapped u across its width and v from its front (0) to its back
 * (1): brown-black polymerised oil over the steel, each section's working area scraped back to burnished metal in
 * front-to-back strokes, the edges along the splashes less used, lighter and heat-tinted, carbon specks. The rough
 * and metal map carries roughness in green and metalness in blue: the burnished metal smooth and bare, the
 * seasoning glossy but dielectric, the carbon dry.
 */
const plateCache = new Map<number, THREE.MeshStandardMaterial>();
function seasonedPlate(n: number): THREE.MeshStandardMaterial {
  const hit = plateCache.get(n);
  if (hit) return hit;
  const r = rand(4211 + n * 97);
  const W = 400 * n;
  const H = 800;
  const [c, x] = canvas2d(W, H);
  const [cm, xm] = canvas2d(W / 2, H / 2);
  // Seasoned ground: near-black brown, its coats uneven.
  x.fillStyle = "#1e1813";
  x.fillRect(0, 0, W, H);
  xm.fillStyle = "rgb(0,120,70)";
  xm.fillRect(0, 0, W / 2, H / 2);
  const blob = (ctx: CanvasRenderingContext2D, cx: number, cy: number, rad: number, col: string, a: number) => {
    const g = ctx.createRadialGradient(cx, cy, 0, cx, cy, rad);
    g.addColorStop(0, col.replace("A", a.toFixed(3)));
    g.addColorStop(1, col.replace("A", "0"));
    ctx.fillStyle = g;
    ctx.fillRect(cx - rad, cy - rad, rad * 2, rad * 2);
  };
  for (let i = 0; i < 220 * n; i++) {
    const v = r();
    blob(x, r() * W, r() * H, 30 + r() * 140, v < 0.5 ? "rgba(10,8,6,A)" : v < 0.85 ? "rgba(52,38,26,A)" : "rgba(70,50,30,A)", 0.12 + r() * 0.22);
  }
  // The edges, outside the cooks' reach: lighter, browner, tinted straw and blue by the heat.
  const edge = 0.055 * H;
  for (const [x0, y0, w0, h0, gx0, gy0, gx1, gy1] of [
    [0, 0, edge, H, 0, 0, edge, 0],
    [W - edge, 0, edge, H, W, 0, W - edge, 0],
    [0, 0, W, edge * 1.3, 0, 0, 0, edge * 1.3],
  ] as const) {
    const g = x.createLinearGradient(gx0, gy0, gx1, gy1);
    g.addColorStop(0, "rgba(104,86,62,0.85)");
    g.addColorStop(0.45, "rgba(78,66,58,0.55)");
    g.addColorStop(1, "rgba(40,32,26,0)");
    x.fillStyle = g;
    x.fillRect(x0, y0, w0, h0);
    const gm = xm.createLinearGradient(gx0 / 2, gy0 / 2, gx1 / 2, gy1 / 2);
    gm.addColorStop(0, "rgba(0,150,175,0.9)");
    gm.addColorStop(1, "rgba(0,150,175,0)");
    xm.fillStyle = gm;
    xm.fillRect(x0 / 2, y0 / 2, w0 / 2, h0 / 2);
  }
  for (let i = 0; i < 90 * n; i++) {
    const side = r() < 0.5;
    const px = side ? (r() < 0.5 ? r() * edge : W - r() * edge) : r() * W;
    const py = side ? r() * H : r() * edge;
    blob(x, px, py, 6 + r() * 26, r() < 0.5 ? "rgba(70,78,104,A)" : "rgba(128,98,52,A)", 0.12 + r() * 0.2);
  }
  // Each section's working area, scraped back to the steel: a soft-edged burnish biased to the front.
  for (let k = 0; k < n; k++) {
    const cx = (k + 0.5) * (W / n) + (r() - 0.5) * 30;
    const cy = H * (0.6 + r() * 0.06);
    const rx = (W / n) * (0.42 + r() * 0.05);
    const ry = H * (0.36 + r() * 0.04);
    for (let j = 0; j < 4; j++) {
      x.save();
      x.translate(cx + (r() - 0.5) * 40, cy + (r() - 0.5) * 40);
      x.scale(1, ry / rx);
      const g = x.createRadialGradient(0, 0, 0, 0, 0, rx * (0.75 + j * 0.12));
      g.addColorStop(0, "rgba(104,98,90,0.5)");
      g.addColorStop(0.55, "rgba(82,76,68,0.32)");
      g.addColorStop(1, "rgba(60,52,44,0)");
      x.fillStyle = g;
      x.fillRect(-rx * 1.2, -rx * 1.2, rx * 2.4, rx * 2.4);
      x.restore();
      xm.save();
      xm.translate(cx / 2, cy / 2);
      xm.scale(1, ry / rx);
      const gm = xm.createRadialGradient(0, 0, 0, 0, 0, (rx / 2) * (0.8 + j * 0.1));
      gm.addColorStop(0, "rgba(0,62,235,0.5)");
      gm.addColorStop(1, "rgba(0,62,235,0)");
      xm.fillStyle = gm;
      xm.fillRect(-rx, -rx, rx * 2, rx * 2);
      xm.restore();
    }
    // Scraper strokes: short, fine and many, mostly front to back toward the trough, the rest every way.
    for (let i = 0; i < 2600; i++) {
      const a = Math.sqrt(r()) * 0.98;
      const t = r() * Math.PI * 2;
      const px = cx + Math.cos(t) * rx * a;
      const py = cy + Math.sin(t) * ry * a;
      const ang = r() < 0.6 ? Math.PI / 2 + (r() - 0.5) * 0.5 : r() * Math.PI;
      const len = 5 + r() * r() * 70;
      const bright = r() < 0.55;
      const fade = 1 - a * 0.6;
      x.strokeStyle = bright ? `rgba(${160 + ((r() * 50) | 0)},${152 + ((r() * 46) | 0)},${140 + ((r() * 40) | 0)},${((0.025 + r() * 0.05) * fade).toFixed(3)})` : `rgba(12,9,7,${(0.03 + r() * 0.07).toFixed(3)})`;
      x.lineWidth = 0.5 + r() * 1.1;
      x.beginPath();
      x.moveTo(px, py);
      x.lineTo(px + Math.cos(ang) * len, py + Math.sin(ang) * len);
      x.stroke();
      if (bright && i % 4 === 0) {
        xm.strokeStyle = "rgba(0,40,255,0.12)";
        xm.lineWidth = 0.6;
        xm.beginPath();
        xm.moveTo(px / 2, py / 2);
        xm.lineTo((px + Math.cos(ang) * len) / 2, (py + Math.sin(ang) * len) / 2);
        xm.stroke();
      }
    }
    // Burnished smears where the spatula's heel rides, and dark islands of fresh seasoning in the working area.
    for (let i = 0; i < 40; i++) {
      const t = r() * Math.PI * 2;
      const a = Math.sqrt(r()) * 0.8;
      const bright = r() < 0.6;
      blob(x, cx + Math.cos(t) * rx * a, cy + Math.sin(t) * ry * a, 10 + r() * 40, bright ? "rgba(128,122,112,A)" : "rgba(22,16,12,A)", 0.1 + r() * 0.18);
    }
    // A brick's swirl here and there.
    for (let i = 0; i < 6; i++) {
      const sx = cx + (r() - 0.5) * rx;
      const sy = cy + (r() - 0.5) * ry;
      for (let q = 0; q < 26; q++) {
        x.strokeStyle = `rgba(170,162,150,${(0.03 + r() * 0.05).toFixed(3)})`;
        x.lineWidth = 0.7;
        x.beginPath();
        x.arc(sx + (r() - 0.5) * 8, sy + (r() - 0.5) * 8, 10 + r() * 40, r() * 6, r() * 6 + 1 + r() * 2);
        x.stroke();
      }
    }
  }
  // The front band the scraper rides off into the trough: worn bright.
  {
    const g = x.createLinearGradient(0, H, 0, H - 26);
    g.addColorStop(0, "rgba(120,114,104,0.7)");
    g.addColorStop(1, "rgba(120,114,104,0)");
    x.fillStyle = g;
    x.fillRect(0, H - 26, W, 26);
  }
  // Carbon: specks and the odd burnt patch, dry and dull.
  for (let i = 0; i < 450 * n; i++) {
    const px = r() * W;
    const py = r() * H;
    const sz = 0.5 + r() * (r() < 0.04 ? 4 : 1.2);
    x.fillStyle = `rgba(6,5,4,${(0.15 + r() * 0.35).toFixed(2)})`;
    x.beginPath();
    x.arc(px, py, sz, 0, Math.PI * 2);
    x.fill();
    if (sz > 2.5) {
      xm.fillStyle = "rgba(0,230,40,0.7)";
      xm.fillRect(px / 2, py / 2, sz / 2, sz / 2);
    }
  }
  // The bottom rows: the plate's flame-cut edges, plain dark steel (its sides are mapped here).
  x.fillStyle = "#3a3734";
  x.fillRect(0, H - 6, W, 6);
  xm.fillStyle = "rgb(0,150,200)";
  xm.fillRect(0, H / 2 - 3, W / 2, 3);
  const map = canvasTex(c, true);
  const rm = canvasTex(cm, false);
  const m = new THREE.MeshStandardMaterial({ color: 0xffffff, map, roughnessMap: rm, metalnessMap: rm, roughness: 1, metalness: 1, envMapIntensity: 1.1 });
  m.userData.uv = "own";
  plateCache.set(n, m);
  return m;
}

/**
 * The printed parts, one atlas: dial plates in the top half (eight 256 px slots, 4 by 2) and rating plates and
 * labels in the bottom half (four 512 by 256 slots). Dials are screen-printed black on satin stainless.
 */
type DialKind = "gas" | "griddle" | "power" | "temp450" | "fryer" | "steps";
const DIAL_SLOT: Record<DialKind, number> = { gas: 0, griddle: 1, power: 2, temp450: 3, fryer: 4, steps: 5 };
type LabelKind = "gas" | "elec" | "caution" | "elec120";
const LABEL_SLOT: Record<LabelKind, number> = { gas: 0, elec: 1, caution: 2, elec120: 3 };
let printCache: THREE.MeshStandardMaterial | null = null;
function printMat(): THREE.MeshStandardMaterial {
  if (printCache) return printCache;
  const S = 1024;
  const [c, x] = canvas2d(S, S);
  x.fillStyle = "#c9ccd0";
  x.fillRect(0, 0, S, S);
  // Dials: a ring of graduations and figures; OFF at the top. `a` is the knob turn from OFF (positive is anticlockwise).
  const dial = (slot: number, marks: { a: number; text?: string; major?: boolean }[], arcs: [number, number][] = [], centre?: string) => {
    const ox = (slot % 4) * 256 + 128;
    const oy = Math.floor(slot / 4) * 256 + 128;
    const g = x.createRadialGradient(ox - 30, oy - 40, 10, ox, oy, 128);
    g.addColorStop(0, "#e3e5e8");
    g.addColorStop(1, "#b8bcc1");
    x.fillStyle = g;
    x.beginPath();
    x.arc(ox, oy, 127, 0, Math.PI * 2);
    x.fill();
    for (let i = 0; i < 900; i++) {
      const yy = oy - 127 + Math.random() * 254;
      x.fillStyle = `rgba(255,255,255,${(Math.random() * 0.08).toFixed(3)})`;
      x.fillRect(ox - 127, yy, 254, 1);
    }
    x.strokeStyle = "#141414";
    x.fillStyle = "#141414";
    x.lineWidth = 3;
    x.beginPath();
    x.arc(ox, oy, 123, 0, Math.PI * 2);
    x.stroke();
    for (const [a0, a1] of arcs) {
      x.lineWidth = 5;
      x.beginPath();
      x.arc(ox, oy, 76, -Math.PI / 2 - a1, -Math.PI / 2 - a0);
      x.stroke();
    }
    x.textAlign = "center";
    x.textBaseline = "middle";
    for (const m of marks) {
      const ang = -Math.PI / 2 - m.a; // canvas angles run clockwise from +x
      const [cx, cy] = [Math.cos(ang), Math.sin(ang)];
      x.lineWidth = m.major ? 5 : 3;
      x.beginPath();
      x.moveTo(ox + cx * 80, oy + cy * 80);
      x.lineTo(ox + cx * (m.major ? 96 : 90), oy + cy * (m.major ? 96 : 90));
      x.stroke();
      if (m.text) {
        x.font = `bold ${m.text.length > 3 ? 19 : 23}px sans-serif`;
        x.fillText(m.text, ox + cx * 110, oy + cy * 110);
      }
    }
    if (centre) {
      x.font = "bold 15px sans-serif";
      x.fillText(centre, ox, oy + 66);
    }
  };
  const sweep = (a0: number, a1: number, labels: string[], minors = 1) => {
    const out: { a: number; text?: string; major?: boolean }[] = [];
    const n = labels.length - 1;
    for (let k = 0; k <= n * (minors + 1); k++) {
      const major = k % (minors + 1) === 0;
      out.push({ a: a0 + ((a1 - a0) * k) / (n * (minors + 1)), major, text: major ? labels[k / (minors + 1)] : undefined });
    }
    return out;
  };
  // Gas valve: OFF up, HI at a quarter turn anticlockwise, down to LO at a half turn.
  dial(DIAL_SLOT.gas, [{ a: 0, text: "OFF", major: true }, { a: Math.PI / 2, text: "HI", major: true }, ...sweep(Math.PI / 2, Math.PI, ["", "", "", "", ""]).slice(1, -1), { a: Math.PI, text: "LO", major: true }], [[Math.PI / 2, Math.PI]]);
  // Griddle: the valve's sweep marked in °F, hottest at full flow.
  dial(DIAL_SLOT.griddle, [{ a: 0, text: "OFF", major: true }, ...sweep(Math.PI / 2, Math.PI, ["550", "450", "350", "250"], 1)], [], "°F");
  // Induction power: 0 up, 1 to 10 clockwise to 212°.
  dial(DIAL_SLOT.power, [{ a: 0, text: "0", major: true }, ...sweep(-0.25, -(0.25 + 1.1 * Math.PI), ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10"], 1)], [], "POWER");
  // Braising pan: OFF up, 100 to 450 °F clockwise over three quarters of a turn.
  dial(DIAL_SLOT.temp450, [{ a: 0, text: "OFF", major: true }, ...sweep(-0.3, -1.5 * Math.PI, ["100", "150", "200", "250", "300", "350", "400", "450"], 1)], [], "°F");
  // Fryer thermostat: OFF, then 200 to 400 °F clockwise.
  dial(DIAL_SLOT.fryer, [{ a: 0, text: "OFF", major: true }, ...sweep(-0.4, -1.5 * Math.PI, ["200", "250", "300", "350", "400"], 4)], [], "°F");
  // Kettle steam: OFF, then 1 to 10 clockwise.
  dial(DIAL_SLOT.steps, [{ a: 0, text: "OFF", major: true }, ...sweep(-0.3, -1.5 * Math.PI, ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10"], 1)], [], "HEAT");
  // Rating plates and a caution label in the bottom half.
  const label = (slot: number, draw: (ox: number, oy: number) => void) => draw((slot % 2) * 512, 512 + Math.floor(slot / 2) * 256);
  const rating = (ox: number, oy: number, lines: [string, string][], title: string) => {
    x.fillStyle = "#d4d6d8";
    x.fillRect(ox, oy, 512, 256);
    for (let i = 0; i < 400; i++) {
      x.fillStyle = `rgba(255,255,255,${(Math.random() * 0.1).toFixed(3)})`;
      x.fillRect(ox, oy + Math.random() * 256, 512, 1);
    }
    x.strokeStyle = "#1a1a1a";
    x.lineWidth = 4;
    x.strokeRect(ox + 14, oy + 14, 484, 228);
    x.fillStyle = "#1a1a1a";
    x.fillRect(ox + 14, oy + 14, 484, 34);
    x.fillStyle = "#d4d6d8";
    x.font = "bold 22px sans-serif";
    x.textAlign = "center";
    x.textBaseline = "middle";
    x.fillText(title, ox + 256, oy + 32);
    x.textAlign = "left";
    x.fillStyle = "#1a1a1a";
    lines.forEach(([k, v], i) => {
      const y = oy + 66 + i * 27;
      x.font = "bold 15px sans-serif";
      x.fillText(k, ox + 26, y);
      x.font = "17px monospace";
      x.fillText(v, ox + 250, y);
      x.fillRect(ox + 246, y + 11, 236, 1);
    });
  };
  label(LABEL_SLOT.gas, (ox, oy) =>
    rating(ox, oy, [["MODEL NO.", "CG-36-T"], ["SERIAL NO.", "4471 0926 18"], ["GAS", "NATURAL"], ["TOTAL INPUT", "90,000 BTU/HR"], ["MANIFOLD PRESS.", '4.0" W.C.'], ["MIN. CLEARANCE", 'SIDES 6"  BACK 6"']], "COMMERCIAL COOKING APPLIANCE"),
  );
  label(LABEL_SLOT.elec, (ox, oy) =>
    rating(ox, oy, [["MODEL NO.", "EK-40-3"], ["SERIAL NO.", "5528 1034 07"], ["VOLTS", "208 / 3 PH"], ["KW", "18.0"], ["AMPS", "50.0"], ["HZ", "60"]], "COMMERCIAL COOKING APPLIANCE"),
  );
  label(LABEL_SLOT.elec120, (ox, oy) =>
    rating(ox, oy, [["MODEL NO.", "IR-48-S"], ["SERIAL NO.", "6093 2210 45"], ["VOLTS", "120 / 1 PH"], ["WATTS", "1100"], ["AMPS", "9.2"], ["HZ", "60"]], "INFRARED FOOD WARMER"),
  );
  label(LABEL_SLOT.caution, (ox, oy) => {
    x.fillStyle = "#f2c31b";
    x.fillRect(ox, oy, 512, 256);
    x.fillStyle = "#111";
    x.fillRect(ox + 12, oy + 12, 488, 232);
    x.fillStyle = "#f2c31b";
    x.fillRect(ox + 20, oy + 20, 472, 216);
    x.fillStyle = "#111";
    x.beginPath();
    x.moveTo(ox + 90, oy + 52);
    x.lineTo(ox + 150, oy + 160);
    x.lineTo(ox + 30, oy + 160);
    x.closePath();
    x.fill();
    x.fillStyle = "#f2c31b";
    x.font = "bold 70px sans-serif";
    x.textAlign = "center";
    x.textBaseline = "middle";
    x.fillText("!", ox + 90, oy + 122);
    x.fillStyle = "#111";
    x.font = "bold 46px sans-serif";
    x.fillText("CAUTION", ox + 330, oy + 80);
    x.font = "bold 26px sans-serif";
    x.fillText("HOT SURFACE", ox + 330, oy + 136);
    x.font = "18px sans-serif";
    x.fillText("DO NOT TOUCH WHEN IN USE", ox + 256, oy + 200);
  });
  const m = new THREE.MeshStandardMaterial({ color: 0xffffff, map: canvasTex(c, true), metalness: 0.55, roughness: 0.38 });
  m.userData.uv = "own";
  printCache = m;
  return m;
}

/** A dial plate under a knob: a printed disc of radius `r`, facing +z, standing 0.4 mm proud of the face at `at`. */
function dialPlate(s: Shop, kind: DialKind, at: V3, r: number): void {
  if (!s.fine(2)) return;
  const slot = DIAL_SLOT[kind];
  const [col, row] = [slot % 4, Math.floor(slot / 4)];
  const g = new THREE.CylinderGeometry(r, r, 0.0005, s.seg(48, 16), 1).rotateX(Math.PI / 2);
  mapUV(g, (x, y, _z, _nx, _ny, nz) => (nz > 0.5 ? [(col + 0.5 + (x / r) * 0.496) / 4, 1 - (row + 0.5 - (y / r) * 0.496) / 4] : [(col + 0.5) / 4, 1 - (row + 0.03) / 4]));
  s.add(printMat(), g, [at[0], at[1], at[2] + 0.00035]);
}

/** A rating plate or label `w` by `h`, facing +z at `at` (turned by `rot`), on four rivets, 0.3 mm proud. */
function ratingPlate(s: Shop, kind: LabelKind, at: V3, w: number, h: number, rot: V3 = [0, 0, 0]): void {
  if (!s.fine(1)) return;
  const p = palette();
  s.at(at, rot, () => {
    if (!s.fine(2)) {
      s.add(p.ssBright, box(w, h, 0.0006), [0, 0, 0.0003]);
      return;
    }
    const slot = LABEL_SLOT[kind];
    const [col, row] = [slot % 2, Math.floor(slot / 2)];
    const g = rbox(s, w, h, 0.0006, 0.0002);
    mapUV(g, (x, y, _z, _nx, _ny, nz) => (nz > 0.5 ? [(col + 0.5 + x / w) / 2, 1 - (2 + row + 0.5 - y / h) / 4] : [(col + 0.02) / 2, 1 - (2 + row + 0.02) / 4]));
    s.add(printMat(), g, [0, 0, 0.0003]);
    if (kind !== "caution") for (const [sx, sy] of [[-1, -1], [1, -1], [-1, 1], [1, 1]] as const) s.add(p.ssBright, ball(s, 0.0017, 10, Math.PI / 2), [sx * (w / 2 - 0.0045), sy * (h / 2 - 0.0045), 0.0005], [Math.PI / 2, 0, 0]);
  });
}

/**
 * Red heat: steel from dark through dull cherry to orange as `heat` rises, the emission held below what tone
 * mapping would bleach to cream. Radiants, sheathed elements and their glow maps.
 */
function setRedHeat(m: THREE.MeshStandardMaterial, heat: number, peak = 1.5): void {
  const h = clamp01(heat);
  const k = smooth(0.12, 1, h);
  m.emissive.setRGB(1, 0.07 + 0.25 * k * k, 0.008 + 0.03 * k * k * k);
  m.emissiveIntensity = h <= 0.12 ? 0 : peak * Math.pow(k, 1.2);
}

/**
 * An infrared plaque's heat: the ceramic face burning from dull red to orange, its own cream giving way to the glow
 * as it heats (the base colour it was made with is kept in userData.cold).
 */
function setCeramicHeat(m: THREE.MeshStandardMaterial, heat: number, peak = 1.6): void {
  const h = clamp01(heat);
  const k = smooth(0.12, 1, h);
  m.userData.cold ??= m.color.getHex();
  m.color.setHex(m.userData.cold as number).multiplyScalar(1 - 0.85 * k);
  m.emissive.setRGB(1, 0.1 + 0.2 * k * k, 0.01 + 0.02 * k * k);
  m.emissiveIntensity = h <= 0.12 ? 0 : peak * Math.pow(k, 1.1);
}

/** A glowing part whose own map shades its glow (hotter where `glowMap` is bright): own UVs. */
function heatMaterial(base: number, glowMap: THREE.Texture | null, map: THREE.Texture | null = null, rough = 0.6, metal = 0.4): THREE.MeshStandardMaterial {
  const m = new THREE.MeshStandardMaterial({ color: base, metalness: metal, roughness: rough, emissive: 0xff3000, emissiveIntensity: 0, emissiveMap: glowMap, map });
  m.userData.uv = "own";
  return m;
}

/** A radiant's glow along its leaf: brightest along the ridge (v = 1), mottled, dark toward its foot (v = 0). */
let radiantTex: [THREE.CanvasTexture, THREE.CanvasTexture] | null = null;
function radiantMaps(): [THREE.CanvasTexture, THREE.CanvasTexture] {
  if (radiantTex) return radiantTex;
  const r = rand(733);
  const [c, x] = canvas2d(256, 128);
  const [cb, xb] = canvas2d(256, 128);
  const g = x.createLinearGradient(0, 128, 0, 0);
  g.addColorStop(0, "#050302");
  g.addColorStop(0.35, "#2a1408");
  g.addColorStop(0.75, "#a8a8a8");
  g.addColorStop(1, "#ffffff");
  x.fillStyle = g;
  x.fillRect(0, 0, 256, 128);
  // Heat-tinted steel: bronze, straw and blue under the soot.
  xb.fillStyle = "#57483a";
  xb.fillRect(0, 0, 256, 128);
  for (let i = 0; i < 500; i++) {
    const [px, py, s] = [r() * 256, r() * 128, 4 + r() * 22];
    x.fillStyle = `rgba(0,0,0,${(r() * 0.18).toFixed(3)})`;
    x.fillRect(px, py, s * 2, s * 0.5);
    const v = r();
    xb.fillStyle = v < 0.3 ? "rgba(52,58,92,0.25)" : v < 0.6 ? "rgba(140,104,58,0.22)" : "rgba(20,16,14,0.3)";
    xb.fillRect(px, py, s * 2, s * 0.6);
  }
  const a = canvasTex(c, true);
  const b = canvasTex(cb, true);
  for (const t of [a, b]) t.wrapS = THREE.RepeatWrapping;
  radiantTex = [a, b];
  return radiantTex;
}

/**
 * An infrared burner's face: 2 by 2 cemented ceramic tiles of fine ports, cream cold. Its glow map burns hottest
 * mid-tile, a little unevenly from tile to tile, cooler at the cemented joints, the ports brighter than the web.
 */
let ceramicTex: [THREE.CanvasTexture, THREE.CanvasTexture] | null = null;
function ceramicMaps(): [THREE.CanvasTexture, THREE.CanvasTexture] {
  if (ceramicTex) return ceramicTex;
  const S = 1024;
  const T = S / 2;
  const r = rand(4411);
  const [c, x] = canvas2d(S, S);
  const [ce, xe] = canvas2d(S, S);
  x.fillStyle = "#d4c9b3";
  x.fillRect(0, 0, S, S);
  xe.fillStyle = "#000000";
  xe.fillRect(0, 0, S, S);
  for (let ty = 0; ty < 2; ty++)
    for (let tx = 0; tx < 2; tx++) {
      const [mx, my] = [tx * T + T / 2 + (r() - 0.5) * 60, ty * T + T / 2 + (r() - 0.5) * 60];
      const v = 0.82 + r() * 0.18;
      const g = xe.createRadialGradient(mx, my, T * 0.05, mx, my, T * 0.75);
      g.addColorStop(0, `rgb(${(255 * v) | 0},${(255 * v) | 0},${(255 * v) | 0})`);
      g.addColorStop(0.7, `rgb(${(200 * v) | 0},${(200 * v) | 0},${(200 * v) | 0})`);
      g.addColorStop(1, `rgb(${(120 * v) | 0},${(120 * v) | 0},${(120 * v) | 0})`);
      xe.fillStyle = g;
      xe.fillRect(tx * T, ty * T, T, T);
      // Soot and wear on the cold face, a few darker clouds per tile.
      for (let k = 0; k < 6; k++) {
        const [px, py, rr] = [tx * T + r() * T, ty * T + r() * T, 20 + r() * 80];
        const gs = x.createRadialGradient(px, py, 0, px, py, rr);
        gs.addColorStop(0, `rgba(90,78,64,${(0.08 + r() * 0.12).toFixed(3)})`);
        gs.addColorStop(1, "rgba(90,78,64,0)");
        x.fillStyle = gs;
        x.fillRect(px - rr, py - rr, rr * 2, rr * 2);
      }
    }
  // Brighter patches where the burning settles, darker where a port row runs lean.
  for (let k = 0; k < 160; k++) {
    const [px, py, rr] = [r() * S, r() * S, 14 + r() * 60];
    const hot = r() < 0.55;
    const gs = xe.createRadialGradient(px, py, 0, px, py, rr);
    gs.addColorStop(0, hot ? `rgba(255,255,255,${(0.08 + r() * 0.16).toFixed(3)})` : `rgba(0,0,0,${(0.05 + r() * 0.12).toFixed(3)})`);
    gs.addColorStop(1, hot ? "rgba(255,255,255,0)" : "rgba(0,0,0,0)");
    xe.fillStyle = gs;
    xe.fillRect(px - rr, py - rr, rr * 2, rr * 2);
  }
  // The ports: round, staggered row to row; dark holes in the cold face, darker than the glowing web when lit.
  const pitch = 8;
  x.fillStyle = "rgba(70,60,50,0.9)";
  xe.fillStyle = "rgba(0,0,0,0.55)";
  for (let j = 0; j < S / pitch; j++)
    for (let i = 0; i < S / pitch; i++) {
      const px = i * pitch + (j % 2 ? pitch / 2 : 0) + pitch / 2;
      const py = j * pitch + pitch / 2;
      x.beginPath();
      x.arc(px, py, 1.9, 0, Math.PI * 2);
      x.fill();
      xe.beginPath();
      xe.arc(px, py, 2.1, 0, Math.PI * 2);
      xe.fill();
    }
  // The cemented joints between tiles, and round the edge.
  for (const k of [0, 1, 2]) {
    for (const [ctx, col] of [[x, "#8a7f6d"], [xe, "#141414"]] as const) {
      ctx.fillStyle = col;
      const at = Math.min(S - 5, Math.max(0, k * T - 5));
      ctx.fillRect(at, 0, k === 1 ? 10 : 5, S);
      ctx.fillRect(0, at, S, k === 1 ? 10 : 5);
    }
  }
  const a = canvasTex(c, true);
  const b = canvasTex(ce, true);
  for (const tx of [a, b]) tx.anisotropy = 8;
  return (ceramicTex = [a, b]);
}

/**
 * A tileable normal map of gentle surface movement: convection rolling hot oil, a simmer's swell. Built from a
 * height field of whole-period waves, so it repeats without a seam.
 */
let rippleCache: THREE.DataTexture | null = null;
function rippleNormal(): THREE.DataTexture {
  if (rippleCache) return rippleCache;
  const S = 256;
  const r = rand(5150);
  const waves = Array.from({ length: 14 }, () => ({ kx: Math.round((r() - 0.5) * 12), ky: Math.round((r() - 0.5) * 12), ph: r() * 6.283, a: 0.4 + r() }));
  const amp = waves.map((w) => w.a / (1 + Math.hypot(w.kx, w.ky)));
  const hgt = (x: number, y: number) => {
    let acc = 0;
    for (let k = 0; k < waves.length; k++) {
      const w = waves[k]!;
      acc = acc + amp[k]! * Math.sin(((w.kx * x + w.ky * y) / S) * 6.283 + w.ph);
    }
    return acc;
  };
  // The height at every whole texel from -1 to S, each read by up to four neighbours' slopes.
  const W = S + 2;
  const hg = new Float64Array(W * W);
  for (let y = -1; y <= S; y++) for (let x = -1; x <= S; x++) hg[(y + 1) * W + x + 1] = hgt(x, y);
  const data = new Uint8Array(S * S * 4);
  for (let y = 0; y < S; y++)
    for (let x = 0; x < S; x++) {
      const o = (y + 1) * W + x + 1;
      const dx = (hg[o + 1]! - hg[o - 1]!) * 9;
      const dy = (hg[o + W]! - hg[o - W]!) * 9;
      const l = Math.hypot(dx, dy, 1);
      const i = (y * S + x) * 4;
      data[i] = ((-dx / l) * 0.5 + 0.5) * 255;
      data[i + 1] = ((-dy / l) * 0.5 + 0.5) * 255;
      data[i + 2] = ((1 / l) * 0.5 + 0.5) * 255;
      data[i + 3] = 255;
    }
  const t = new THREE.DataTexture(data, S, S, THREE.RGBAFormat);
  t.wrapS = t.wrapT = THREE.RepeatWrapping;
  t.generateMipmaps = true;
  t.minFilter = THREE.LinearMipmapLinearFilter;
  t.magFilter = THREE.LinearFilter;
  t.needsUpdate = true;
  return (rippleCache = t);
}

/** Frying oil a few days into service: dark amber, glassy, rolling gently with the heat. UVs in metres. */
let oilCache: THREE.MeshPhysicalMaterial | null = null;
function fryOil(): THREE.MeshPhysicalMaterial {
  if (oilCache) return oilCache;
  const n = rippleNormal().clone();
  n.repeat.set(5, 5);
  n.needsUpdate = true;
  const m = new THREE.MeshPhysicalMaterial({
    color: 0x4a2604,
    metalness: 0,
    roughness: 0.015,
    transparent: true,
    opacity: 0.93,
    ior: 1.47,
    specularIntensity: 1,
    envMapIntensity: 1.25,
    normalMap: n,
    normalScale: new THREE.Vector2(0.35, 0.35),
  });
  m.userData.uv = "own";
  return (oilCache = m);
}
/** Fryer foam: the bubbles that break on hot oil, pale gold and wet. */
let foamCache: THREE.MeshPhysicalMaterial | null = null;
function oilFoam(): THREE.MeshPhysicalMaterial {
  if (foamCache) return foamCache;
  const m = new THREE.MeshPhysicalMaterial({ color: 0xd9a54a, metalness: 0, roughness: 0.08, transparent: true, opacity: 0.55, clearcoat: 1, depthWrite: false });
  m.userData.uv = "own";
  return (foamCache = m);
}

/** Red moulded rubber: igniter caps, reset buttons. */
let redCache: THREE.MeshStandardMaterial | null = null;
function redCap(): THREE.MeshStandardMaterial {
  return (redCache ??= new THREE.MeshStandardMaterial({ color: 0xa3120b, metalness: 0, roughness: 0.42 }));
}

/** Polished anodised reflector aluminium: near mirror. */
let reflectorCache: THREE.MeshStandardMaterial | null = null;
function reflectorMat(): THREE.MeshStandardMaterial {
  if (reflectorCache) return reflectorCache;
  const m = new THREE.MeshStandardMaterial({ color: 0xe9ebee, metalness: 1, roughness: 0.07 });
  m.userData.uv = "own";
  return (reflectorCache = m);
}

/**
 * A palette material taking a part's own UVs (in metres) instead of box-projected ones: for formed surfaces that
 * curve round from one face to the next, where a box projection would break the grain at every change of facing.
 */
const ownCache = new Map<THREE.Material, THREE.Material>();
function ownGrain(m: THREE.Material): THREE.Material {
  let o = ownCache.get(m);
  if (!o) {
    o = m.clone();
    o.userData = { ...m.userData, uv: "own" };
    ownCache.set(m, o);
  }
  return o;
}

/** A simmering stock: golden brown, clear enough to show a little depth, glossy, rolling gently. */
let stockCache: THREE.MeshPhysicalMaterial | null = null;
function stock(): THREE.MeshPhysicalMaterial {
  if (stockCache) return stockCache;
  const n = rippleNormal().clone();
  n.repeat.set(2, 2);
  n.needsUpdate = true;
  const m = new THREE.MeshPhysicalMaterial({
    color: 0x7a4a14,
    metalness: 0,
    roughness: 0.04,
    transparent: true,
    opacity: 0.9,
    ior: 1.34,
    specularIntensity: 1,
    envMapIntensity: 1.2,
    normalMap: n,
    normalScale: new THREE.Vector2(0.1, 0.1),
  });
  m.userData.uv = "own";
  return (stockCache = m);
}
/** The simmer's bubbles and the scum they raise: the stock's colour, paler. */
let stockFoamCache: THREE.MeshStandardMaterial | null = null;
function stockFoam(): THREE.MeshStandardMaterial {
  return (stockFoamCache ??= new THREE.MeshStandardMaterial({ color: 0xc89a5a, roughness: 0.25, metalness: 0, transparent: true, opacity: 0.85 }));
}

/** Cast iron as it comes off a charbroiler: black, sand-cast pitted, the tops of the bars rubbed to a dull shine. */
let castCache: THREE.MeshStandardMaterial | null = null;
function castGrate(): THREE.MeshStandardMaterial {
  if (castCache) return castCache;
  const r = rand(9127);
  const S = 256;
  const [c, x] = canvas2d(S, S);
  const [cb, xb] = canvas2d(S, S);
  x.fillStyle = "#2b2a2a";
  x.fillRect(0, 0, S, S);
  xb.fillStyle = "#808080";
  xb.fillRect(0, 0, S, S);
  for (let i = 0; i < 9000; i++) {
    const v = (r() * 70) | 0;
    const [px, py, s] = [r() * S, r() * S, 0.6 + r() * 2.2];
    x.fillStyle = `rgba(${v},${v * 0.95},${v * 0.9},${(0.15 + r() * 0.35).toFixed(2)})`;
    x.fillRect(px, py, s, s);
    const b = r() < 0.5 ? 40 : 200;
    xb.fillStyle = `rgba(${b},${b},${b},${(0.2 + r() * 0.4).toFixed(2)})`;
    xb.beginPath();
    xb.arc(px, py, s * 0.7, 0, Math.PI * 2);
    xb.fill();
  }
  const map = canvasTex(c, true);
  const bump = canvasTex(cb, false);
  for (const t of [map, bump]) {
    t.wrapS = t.wrapT = THREE.RepeatWrapping;
  }
  const m = new THREE.MeshStandardMaterial({ color: 0xffffff, map, bumpMap: bump, bumpScale: 1.2, metalness: 0.45, roughness: 0.72 });
  m.userData.uvScale = 9;
  return (castCache = m);
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
  /** The splashes' foot (the cooking surface), and the side splashes' outside face. */
  splashBase: number;
  splashX: number;
  /** The cooking surface's front edge (the trough's back wall) and back edge (the back splash). */
  zTop0: number;
  zTop1: number;
  /** The splashes' top. */
  yS: number;
  /** The grease trough's floor, the top of its back wall, and the x of the drain hole in its floor. */
  troughY: number;
  troughTop: number;
  drainX: number;
  /** Openings cut in the control panel, traced in its (x, y) the hole way round. */
  holes: THREE.Path[];
}

/** The bull-nose front rail's centre line: height, depth, radius. */
const NOSE: [number, number, number] = [0.289, -0.019, 0.013];

/** A sheet cut to an outline drawn in the (z, y) plane, `t` thick from x = 0 toward +x, with holes likewise. */
function zyPlate(s: Shop, pts: readonly P2[], t: number, holes: readonly (readonly P2[])[] = []): THREE.BufferGeometry {
  const sh = new THREE.Shape(pts.map(([z, y]) => new THREE.Vector2(-z, y)));
  for (const h of holes) sh.holes.push(new THREE.Path(h.map(([z, y]) => new THREE.Vector2(-z, y))));
  return plate(s, sh, t, 0, 8).rotateY(Math.PI / 2);
}

/** A sheet cut to an outline drawn in plan, (x, z), lying `t` thick from y = 0 upward, with holes likewise. */
function xzPlate(s: Shop, outline: THREE.Shape, t: number): THREE.BufferGeometry {
  return plate(s, outline, t, 0, 12).rotateX(-Math.PI / 2);
}

/** A weld bead or fillet: a fine round bar along a straight run, ground smooth. Only drawn close up. */
function bead(s: Shop, mat: THREE.Material, a: V3, b: V3, r = 0.0011): void {
  if (!s.fine(2)) return;
  s.add(mat, pipe(s, [a, b], r, r, 6));
}

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
  // 4 in adjustable legs, each in a socket on a welded pad.
  for (const lx of [0.045, w - 0.045])
    for (const lz of [zf - 0.05, back + 0.055]) {
      flangedLeg(s, lx, lz, legH + 0.002, 0.016);
      s.add(p.ssDull, sheet(s, 0.05, 0.005, 0.05, 0.001), [lx, legH - 0.0025, lz]);
    }
  // The bottom pan; the sides, each one sheet cut to the body and rising behind the back splash round the flue
  // chamber; the back.
  // (It stops behind the control panel, whose lower edge would otherwise lie flush with its own.)
  s.add(p.ssDull, sheet(s, w - 0.0032, 0.0016, bodyD - 0.0012), [w / 2, legH + 0.0008, (zf - 0.0012 + back) / 2]);
  const rr = 0.005;
  const side: P2[] = [[back, legH], [zf, legH], [zf, o.sideTop - rr], ...arc2([zf - rr, o.sideTop - rr], rr, 0, Math.PI / 2, s.seg(5, 2)).slice(1), [o.zTop1, o.sideTop], [o.zTop1, yS], [back, yS]];
  for (const x of [0, w - 0.0016]) s.add(p.ss, zyPlate(s, side, 0.0016), [x, 0, 0]);
  s.add(p.ssDull, sheet(s, w - 0.0032, yS - legH, 0.0016), [w / 2, (yS + legH) / 2, back + 0.0008]);
  // The back splash, its top folded back into the flue cap, slotted where the flue gases leave.
  s.add(p.ss, sheet(s, w - 0.0032, yS - o.splashBase, 0.0016), [w / 2, (yS + o.splashBase) / 2, o.zTop1 - 0.0008]);
  const capD = o.zTop1 - back + 0.0016;
  const capZ = (o.zTop1 + back) / 2 - 0.0008;
  const slots = Math.max(2, Math.round((w - 0.1) / 0.11));
  const sw = (w - 0.1) / slots;
  const slotZ = (zz: number) => o.zTop1 + (back - o.zTop1) * zz;
  if (s.fine(2)) {
    // Close up the slots are cut through, the flue chamber sooted dark below them.
    const cap = roundRect(w, capD, 0, w / 2, -capZ);
    for (let i = 0; i < slots; i++) for (const zz of [0.33, 0.67]) cap.holes.push(rectPath(new THREE.Path(), sw - 0.016, 0.009, 0.0045, 0.05 + sw * (i + 0.5), -slotZ(zz), true));
    s.add(p.ss, plate(s, cap, 0.0016, 0, 5).rotateX(-Math.PI / 2), [0, yS, 0]);
    s.add(p.enamel, box(w - 0.006, 0.001, capD - 0.006), [w / 2, yS - 0.03, capZ]);
  } else {
    s.add(p.ss, box(w, 0.0016, capD), [w / 2, yS + 0.0008, capZ]);
    if (s.fine(1)) for (let i = 0; i < slots; i++) for (const zz of [0.33, 0.67]) s.add(p.enamel, box(sw - 0.016, 0.0004, 0.009), [0.05 + sw * (i + 0.5), yS + 0.0018, slotZ(zz)]);
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
  for (const [sx, inner] of [[o.splashX, o.splashX + 0.0016], [w - o.splashX - 0.0016, w - o.splashX - 0.0016]] as const) {
    const g = plate(s, sh, 0.0016, 0, 6);
    g.rotateY(Math.PI / 2);
    s.add(p.ss, g, [sx, o.splashBase, zFront]);
    if (s.fine(2)) {
      const x = sx + 0.0008;
      s.add(p.ss, pipe(s, [[x, o.splashBase + hf - 0.004, zFront - 0.006], [x, o.splashBase + hf, zFront - 0.02], [x, yS, o.zTop1]], 0.0019, 0.01, 8));
      // The fillet weld along its foot, inside.
      const xi = inner + (sx < w / 2 ? 0.0004 : -0.0004);
      bead(s, p.ss, [xi, o.splashBase + 0.0003, zFront], [xi, o.splashBase + 0.0003, o.zTop1 + 0.0006]);
    }
  }
  bead(s, p.ss, [o.splashX + 0.0016, o.splashBase + 0.0003, o.zTop1 + 0.0004], [w - o.splashX - 0.0016, o.splashBase + 0.0003, o.zTop1 + 0.0004]);
  // The grease trough: one formed channel from the cooking surface's front edge to the rail, coved at its foot,
  // a hole at one end of its floor draining down a chute.
  const zt0 = o.zTop0 + 0.0006;
  const zt1 = NOSE[1] - NOSE[2] - 0.0012;
  const ty = o.troughY;
  const [x0, x1] = [0.0016, w - 0.0016];
  if (s.fine(1)) {
    const rc = 0.008;
    const cove = s.seg(8, 3);
    const backWall: P2[] = [[zt0 + rc, ty], ...arc2([zt0 + rc, ty + rc], rc, -Math.PI / 2, -Math.PI, cove).slice(1), [zt0, o.troughTop]];
    const frontWall: P2[] = [[zt1 + 0.004, NOSE[0] + 0.006], [zt1, NOSE[0] + 0.002], [zt1, ty + rc], ...arc2([zt1 - rc, ty + rc], rc, 0, -Math.PI / 2, cove).slice(1)];
    s.add(p.ss, sweepX(backWall, x0, x1));
    s.add(p.ss, sweepX(frontWall, x0, x1));
    const [fz0, fz1] = [zt0 + rc, zt1 - rc];
    const zc = (fz0 + fz1) / 2;
    const fl = roundRect(x1 - x0, fz1 - fz0, 0, w / 2, -zc);
    fl.holes.push(circleHole(o.drainX, -zc, 0.011));
    s.add(p.ss, xzPlate(s, fl, 0.0012), [0, ty - 0.0012, 0]);
    // The chute, seen down the hole; a pressed lip round it.
    s.add(p.ssHeat, lathe(s, [[0.011, 0], [0.011, -0.09]], 18), [o.drainX, ty - 0.0012, zc]);
    if (s.fine(2)) s.add(p.ss, torus(s, 0.0118, 0.0011, "y", Math.PI * 2, 24, 6), [o.drainX, ty + 0.0001, zc]);
    // The seam where the trough's back wall is welded under the cooking surface's edge.
    bead(s, p.ssHeat, [x0, o.troughTop - 0.0004, zt0 + 0.0002], [x1, o.troughTop - 0.0004, zt0 + 0.0002], 0.0009);
  } else {
    s.add(p.ss, box(x1 - x0, 0.0012, zt1 - zt0), [w / 2, ty - 0.0006, (zt0 + zt1) / 2]);
    s.add(p.ss, box(x1 - x0, o.troughTop - ty, 0.0012), [w / 2, (o.troughTop + ty) / 2, zt0 + 0.0006]);
    s.add(p.ss, box(x1 - x0, NOSE[0] - ty, 0.0012), [w / 2, (NOSE[0] + ty) / 2, zt1 - 0.0006]);
  }
  // The bull-nose rail across the front, its ends closed by pressed caps.
  s.add(p.ss, cyl(s, NOSE[2], w - 0.004, "x", 28), [w / 2, NOSE[0], NOSE[1]]);
  if (s.fine(2)) for (const [x, sx] of [[0.002, -1], [w - 0.002, 1]] as const) s.add(p.ss, ball(s, NOSE[2], 24, Math.PI / 2), [x, NOSE[0], NOSE[1]], [0, 0, -sx * Math.PI / 2], [1, 0.15, 1]);
  // The control panel, cut for the sight holes and the drawer, screwed at its corners.
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
 * thermostatic valve. Grease is scraped into the front trough and drains through a hole at the left end into a
 * drawer below. Each section has a sight hole beside its knob for the standing pilot.
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
    const drainX = dx0 + 0.04;
    counterShell(s, e, { legH, zf, sideTop: top - 0.0008, splashBase: top, splashX: 0.0022, zTop0: zp0, zTop1: zp1, yS: 0.3775, troughY: 0.2585, troughTop: bot + 0.003, drainX, holes });

    // The plate: 1 in steel, seasoned, its edges eased; the burnish and the patina are in its own map.
    const pw = w - 0.0036;
    const pd = zp0 - zp1;
    const pg = s.fine(1) ? rbox(s, pw, top - bot, pd, 0.0016) : box(pw, top - bot, pd);
    mapUV(pg, (x, _y, z, _nx, ny) => (ny > 0.7 ? [x / pw + 0.5, Math.max(0.012, 0.5 - z / pd)] : [x / pw + 0.5 - z * 0.2, 0.002]));
    s.add(seasonedPlate(n), pg, [w / 2, (top + bot) / 2, (zp0 + zp1) / 2]);

    // Knobs: one thermostat per section, on its printed dial; a chrome grommet round each sight hole.
    const knobs: string[] = [];
    for (let i = 0; i < n; i++) {
      const kx = sec * (i + 0.5) - R;
      knobs.push(valveKnob(s, `knob_${i}`, [kx, ky, zf]));
      dialPlate(s, "griddle", [kx, ky, zf], 0.035);
      if (s.fine(1)) s.add(p.ssBright, torus(s, 0.0088, 0.0015, "z", Math.PI * 2, 20, 6), [sight(i), ky, zf + 0.0006]);
    }
    ratingPlate(s, "gas", [w - 0.075, dy0 + 0.026, zf], 0.08, 0.05);

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

    // The grease drawer: a drawn stainless pan on slides, its front lapping the panel's opening.
    const dr = s.joint("grease_drawer", [dx0, dy0, zf]);
    dr.add(p.ss, sheet(dr, dw, dh + 0.002, 0.0016, 0.0008), [dw / 2, dh / 2, 0.0008]);
    pullHandle0(dr, [dw / 2, dh / 2, 0.0016], Math.min(0.1, dw * 0.5));
    const dd = 0.3;
    if (s.fine(2)) {
      dr.add(p.ssDull, drawnPan(dr, dw - 0.016, dd, dh - 0.008, 0.008, 0.008, 0.001, 0.0025), [dw / 2, dh - 0.005, -dd / 2 - 0.003]);
      // Last service's grease, gone hard and dark in the bottom.
      dr.add(p.griddle, rbox(dr, dw * 0.45, 0.002, dd * 0.35, 0.0009), [dw * 0.4, 0.0042, -dd * 0.7]);
    } else if (s.fine(1)) {
      dr.add(p.ssDull, box(dw - 0.016, 0.001, dd), [dw / 2, 0.002, -dd / 2]);
      for (const sx of [0.008, dw - 0.008]) dr.add(p.ssDull, box(0.001, dh - 0.01, dd), [sx, dh / 2 - 0.003, -dd / 2]);
      dr.add(p.ssDull, box(dw - 0.016, dh - 0.01, 0.001), [dw / 2, dh / 2 - 0.003, -dd]);
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
 * A cast-iron grate bar's section, counter-clockwise in (x, y) from its crown: the narrow sear ridge on top, the
 * flanks swelling to the bar's width and tapering into a keel. `bw` wide and `bh` deep.
 */
function barSection(s: Shop, bw: number, bh: number): P2[] {
  const hw = bw / 2;
  const half: P2[] = s.fine(2)
    ? [[0, 0], [hw * 0.2, -bh * 0.015], [hw * 0.36, -bh * 0.06], [hw * 0.5, -bh * 0.12], [hw * 0.7, -bh * 0.2], [hw * 0.9, -bh * 0.27], [hw, -bh * 0.34], [hw * 0.94, -bh * 0.44], [hw * 0.66, -bh * 0.68], [hw * 0.4, -bh * 0.9], [hw * 0.18, -bh * 0.98], [0, -bh]]
    : [[0, 0], [hw * 0.42, -bh * 0.08], [hw, -bh * 0.34], [hw * 0.45, -bh * 0.9], [0, -bh]];
  const left = half.slice(1).map(([x, y]) => [-x, y] as P2);
  const right = half.slice(1, -1).reverse();
  return [half[0]!, ...left, ...right];
}

/** A bar of section `sec` (counter-clockwise in x, y) running from z = 0 back to z = -len, its ends closed. */
function barAlongZ(sec: readonly P2[], len: number): THREE.BufferGeometry {
  const rows = [0, -len].map((z) => sec.map(([x, y]) => new THREE.Vector3(x, y, z)));
  const side = gridMesh(rows, false, true);
  const n = sec.length;
  let cy = 0;
  for (const [, y] of sec) cy += y / n;
  const cap = (z: number, front: boolean) => {
    const pos = new Float32Array((n + 1) * 3);
    sec.forEach(([x, y], i) => pos.set([x, y, z], i * 3));
    pos.set([0, cy, z], n * 3);
    const idx: number[] = [];
    for (let i = 0; i < n; i++) front ? idx.push(n, i, (i + 1) % n) : idx.push(n, (i + 1) % n, i);
    const g = new THREE.BufferGeometry();
    g.setAttribute("position", new THREE.BufferAttribute(pos, 3));
    g.setAttribute("uv", new THREE.BufferAttribute(new Float32Array((n + 1) * 2), 2));
    g.setIndex(idx);
    g.computeVertexNormals();
    return g;
  };
  return mergeGeometries([side, cap(0, true), cap(-len, false)], false)!;
}

/**
 * A radiant charbroiler: one straight pipe burner per 6 in section, each under a stainless radiant (an inverted V
 * slotted along its ridge, so the flame tips lick through) that glows and vaporises the drippings, under cast-iron
 * grates. The grates' bars are cast with a sear ridge along the crown and slope down to the front trough; what falls
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
  const panY0 = 0.114;
  const panH = 0.04;
  const pw = Math.min(w - 0.12, w * 0.86);
  const px0 = (w - pw) / 2;
  const holes = [rectPath(new THREE.Path(), pw - 0.01, panH - 0.006, 0.002, w / 2, panY0 + panH / 2, true)];
  counterShell(s, e, { legH, zf, sideTop: 0.305, splashBase: 0.305, splashX: 0.0022, zTop0: -0.078, zTop1: zg1 - 0.002, yS: 0.3775, troughY: 0.262, troughTop: 0.303, drainX: 0.06, holes });
  const [radMap, radTint] = radiantMaps();
  const glow = heatMaterial(0xffffff, radMap, radTint, 0.55, 0.55);
  const iron = castGrate();

  // Grate rails across the front and back, angle iron welded to the sides.
  const L = zg0 - zg1;
  const slope = Math.atan2(yb - yf, L);
  for (const [zz, yy, fwd] of [[zg0 - 0.02, yf - 0.03, 1], [zg1 + 0.012, yb - 0.03, -1]] as const) {
    s.add(p.ssHeat, sheet(s, w - 0.0032, 0.003, 0.025), [w / 2, yy - 0.0015, zz]);
    s.add(p.ssHeat, sheet(s, w - 0.0032, 0.025, 0.003), [w / 2, yy - 0.014, zz + fwd * 0.011]);
  }

  const knobs: string[] = [];
  const ports: Port[] = [];
  const nb = 5;
  const bw = 0.019;
  const bh = 0.021;
  const bar = s.fine(1) ? barSection(s, bw, bh) : null;
  for (let i = 0; i < n; i++) {
    const cx = 0.002 + sec * (i + 0.5);
    const gw = sec - 0.004;
    // The grate: one casting of five bars on three ribs, a front lip the bars run over into the trough and a back
    // lug that sits on the rear rail; sloped, so the fat runs forward along the bars.
    s.at([cx, yf, zg0], [slope, 0, 0], () => {
      const pitch = (gw - bw) / (nb - 1);
      for (let k = 0; k < nb; k++) {
        const x = -gw / 2 + bw / 2 + k * pitch;
        if (bar) s.add(iron, barAlongZ(bar, L), [x, 0, 0]);
        else s.add(iron, box(bw, 0.02, L), [x, -0.01, -L / 2]);
      }
      for (const zz of [0.006, L / 2, L - 0.006]) s.add(iron, sheet(s, gw, 0.011, zz === L / 2 ? 0.016 : 0.012, 0.003), [0, zz === L / 2 ? -0.024 : -0.0105, -zz]);
      if (s.fine(2)) {
        s.add(iron, rbox(s, 0.03, 0.008, 0.02, 0.002), [0, -0.004, -L - 0.004], [0.5, 0, 0]);
        // The sprue's witness on the middle rib, the casting's parting line along the outer bars.
        for (const sx of [-1, 1]) s.add(iron, rbox(s, 0.004, 0.012, 0.016, 0.0015), [sx * (gw / 2 - 0.002), -0.026, -L / 2]);
      }
    });
    // The radiant: one sheet folded into an inverted V, flanged at its feet, slotted along its ridge and bridged.
    const ry = 0.284;
    const hw = 0.058;
    const gap = 0.0045;
    const rl = L - 0.03;
    const rz = (zg0 + zg1) / 2;
    if (s.fine(1)) {
      s.at([cx, ry, rz], [0, Math.PI / 2, 0], () => {
        const legL: P2[] = [[-hw - 0.008, -0.0262], [-hw, -0.026], [-hw, -0.026], [-gap - 0.0012, -0.0008], [-gap, 0]];
        const legR = legL.map(([x, y]) => [-x, y] as P2);
        s.add(glow, sweepX(legL, -rl / 2, rl / 2, true));
        s.add(glow, sweepX(legR, -rl / 2, rl / 2));
        if (s.fine(2)) {
          const under = (q: P2[]) => q.map(([x, y]) => [x * 0.985, y - 0.0012] as P2);
          s.add(glow, sweepX(under(legL), -rl / 2, rl / 2));
          s.add(glow, sweepX(under(legR), -rl / 2, rl / 2, true));
        }
      });
      if (s.fine(2)) for (let k = 0; k < 9; k++) s.add(glow, mapUV(box(gap * 2 + 0.002, 0.0014, 0.008), () => [0.5, 0.98]), [cx, ry - 0.0008, rz - rl / 2 + ((k + 0.5) * rl) / 9]);
    } else for (const sx of [-1, 1]) s.add(glow, box(Math.hypot(hw - gap, 0.026), 0.0014, rl), [cx + sx * (gap + (hw - gap) / 2), ry - 0.013, rz], [0, 0, -sx * Math.atan2(0.026, hw - gap)]);
    if (s.fine(1)) for (const zz of [rz + rl / 2, rz - rl / 2]) s.add(p.ssHeat, box(hw * 2 + 0.01, 0.004, 0.003), [cx, ry - 0.028, zz]);
    // The burner: a pipe front to back, ports along its top, its venturi on the valve at the front.
    if (s.fine(1)) {
      s.add(iron, cyl(s, br, 0.57, "z", 14), [cx, by, -0.385]);
      s.add(iron, cyl(s, br, 0.03, "z", 14, 0.018), [cx, by, -0.085]);
      s.add(p.ssDull, cyl(s, 0.02, 0.0015, "z", 16), [cx, by, -0.069]);
      s.add(p.brass, rbox(s, 0.026, 0.03, 0.026, 0.002), [cx, ky, zf - 0.016]);
      s.add(p.brass, pipe(s, [[cx, ky, zf - 0.028], [cx, ky, -0.064], [cx, by, -0.072]], 0.004, 0.006, 8));
    }
    ports.push(...zLinePorts(i, by + br - 0.002, cx, -0.1, -0.66, portCount(s, 18), 0.05));
    knobs.push(valveKnob(s, `knob_${i}`, [cx, ky, zf]));
    dialPlate(s, "gas", [cx, ky, zf], 0.034);
  }
  if (s.fine(1)) s.add(p.enamel, cyl(s, 0.0105, w - 0.06, "x", 12), [w / 2, 0.178, zf - 0.03]);
  ratingPlate(s, "gas", [w, 0.2, -0.36], 0.08, 0.05, [0, Math.PI / 2, 0]);

  // The water pan: under the burners, narrower than the body, its front lapping the panel's opening.
  const dd = 0.56;
  const pan = s.joint("grease_pan", [px0, panY0, zf]);
  pan.add(p.ss, sheet(pan, pw, panH + 0.002, 0.0016, 0.0008), [pw / 2, panH / 2, 0.0008]);
  pullHandle0(pan, [pw / 2, panH / 2, 0.0016], 0.12);
  if (s.fine(2)) {
    pan.add(p.ssDull, drawnPan(pan, pw - 0.012, dd, panH - 0.006, 0.012, 0.01, 0.001, 0.003), [pw / 2, panH - 0.004, -dd / 2 - 0.003]);
    pan.add(p.water, new THREE.PlaneGeometry(pw - 0.034, dd - 0.024).rotateX(-Math.PI / 2), [pw / 2, 0.018, -dd / 2 - 0.003]);
  } else if (s.fine(1)) {
    pan.add(p.ssDull, box(pw - 0.016, 0.001, dd), [pw / 2, 0.002, -dd / 2]);
    for (const sx of [0.008, pw - 0.008]) pan.add(p.ssDull, box(0.001, panH - 0.01, dd), [sx, panH / 2 - 0.003, -dd / 2]);
    pan.add(p.ssDull, box(pw - 0.016, panH - 0.01, 0.001), [pw / 2, panH / 2 - 0.003, -dd]);
    pan.add(p.water, box(pw - 0.02, 0.012, dd - 0.01), [pw / 2, 0.009, -dd / 2]);
  }

  return (b) => {
    const g = gasHooks(b, ports, n, 0.045, knobs);
    const rest: V3 = [px0, panY0, zf];
    g.glows = [glow];
    g.open = (f) => slideJoint(b.joints.get("grease_pan"), "z", clamp01(f) * 0.6 * dd, rest);
    g.set = (name, v) => {
      if (name === "heat") setRedHeat(glow, v, 1.7);
    };
    return g;
  };
};

/* ------------------------------------------------------------------------------------------------------------ */
/* Fryer                                                                                                          */
/* ------------------------------------------------------------------------------------------------------------ */

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
      bulletLeg(s, lx, lz, legH + 0.002, 0.02);
      s.add(p.ssDull, sheet(s, 0.055, 0.005, 0.055, 0.001), [lx, legH - 0.0025, lz]);
    }
  // Cabinet: sides, back, bottom frame, the front frame round the door and the control panel above it. Every
  // skin butts against the next, none laps another's face: the sides stop under the deck, the back, bottom and
  // front frame sit between the sides, the bottom between the back and the front frame.
  const bodyD = zf - back;
  const sk = 0.0016;
  for (const sx of [sk / 2, w - sk / 2]) s.add(p.ss, sheet(s, sk, deck - sk - legH, bodyD, 0.0006), [sx, (deck - sk + legH) / 2, (zf + back) / 2]);
  s.add(p.ssDull, box(w - 2 * sk, deck - legH - 2 * sk, sk), [w / 2, (deck + legH) / 2, back + sk / 2]);
  s.add(p.ssDull, box(w - 2 * sk, sk, bodyD - 2 * sk), [w / 2, legH + sk / 2, (zf + back) / 2]);
  for (const sx of [sk + 0.006, w - sk - 0.006]) s.add(p.ss, sheet(s, 0.012, doorY1 - doorY0, sk, 0.0005), [sx, (doorY1 + doorY0) / 2, zf - sk / 2]);
  s.add(p.ss, box(w - 2 * sk, doorY0 - legH, sk), [w / 2, (doorY0 + legH) / 2, zf - sk / 2]);
  const panY1 = deck - 0.014;
  s.add(p.ss, box(w - 2 * sk, panY1 - doorY1, sk), [w / 2, (panY1 + doorY1) / 2, zf - sk / 2]);
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
  // Its walls stop under the deck's 1.6 mm, so the deck's edge, not a second face, finishes the opening.
  const potTop = deck - 0.0016;
  for (const x of [px0 - 0.001, px1 + 0.001]) s.add(p.ss, box(0.002, potTop - pyB, pz0 - pz1), [x, (potTop + pyB) / 2, (pz0 + pz1) / 2]);
  for (const z of [pz0 + 0.001, pz1 - 0.001]) s.add(p.ss, box(px1 - px0 + 0.004, potTop - pyB + 0.02, 0.002), [w / 2, (potTop + pyB - 0.02) / 2, z]);
  const vw = (px1 - px0) / 2;
  const va = Math.atan2(0.02, vw);
  for (const sx of [-1, 1]) s.add(p.ss, box(Math.hypot(vw, 0.02), 0.002, pz0 - pz1), [w / 2 + (sx * vw) / 2, pyB - 0.01, (pz0 + pz1) / 2], [0, 0, sx * va]);
  // The oil: 50 lb of it, up to the line, dark amber after a few days' frying; it climbs the pot's walls a
  // little in a meniscus.
  {
    const loop = loopRR(px0, px1, pz1, pz0, 0.004, s.fine(2) ? 3 : 1);
    if (s.fine(1)) {
      const men: P2[] = s.fine(2) ? [[0, 0.0026], [0.0006, 0.0016], [0.0018, 0.0008], [0.0036, 0.0003], [0.0065, 0]] : [[0, 0.0022], [0.003, 0.0004], [0.0065, 0]];
      s.add(fryOil(), sweepLoop(loop, men.map(([i, y]) => [i, oilY + y] as P2), true));
      s.add(fryOil(), loopCap(loop, oilY, 0.0065));
    } else s.add(fryOil(), loopCap(loop, oilY, 0));
  }
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
  for (const x of [fx0, fx1]) s.add(p.ss, sheet(s, 0.0016, fTop - deck, fz0 - fz1 - 0.0016), [x, (fTop + deck) / 2, (fz0 + fz1) / 2]);
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
  dialPlate(s, "fryer", [w / 2, ky, zf], 0.036);
  if (s.fine(2)) for (const sx of [0.02, w - 0.02]) screw(s, p.ssBright, [sx, ky, zf], 0.0026);
  // The piezo spark igniter's red cap on the left; the high-limit reset's black button on the right.
  if (s.fine(2)) {
    s.add(p.ssBright, cyl(s, 0.009, 0.003, "z", 18), [0.06, ky, zf + 0.0015]);
    s.add(redCap(), lathe(s, [[0, 0], [0.0072, 0], [0.0072, 0.006], [0.006, 0.0095], [0, 0.0105]], 18).rotateX(Math.PI / 2), [0.06, ky, zf + 0.003]);
    s.add(p.ssBright, cyl(s, 0.0065, 0.003, "z", 16), [w - 0.06, ky, zf + 0.0015]);
    s.add(p.bakelite, cyl(s, 0.0042, 0.005, "z", 14), [w - 0.06, ky, zf + 0.004]);
  }

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
    ratingPlate(s, "gas", [w - 0.0016, 0.42, -0.25], 0.08, 0.05, [0, -Math.PI / 2, 0]);
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
    const bub = bx.map((x) => new Bubbles(oilFoam(), s.fine(2) ? 70 : 30, [x, oilY + 0.0005, downZ - bd / 2], bw - 0.012, bd - 0.02, 0.0045));
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

/**
 * The light a glowing burner throws on the cavity round it, one additive map: the top half a wash fading down a wall
 * from its top edge (soft at both ends), the bottom half a soft pool for the floor below.
 */
let spillTex: THREE.CanvasTexture | null = null;
function spillTexture(): THREE.CanvasTexture {
  if (spillTex) return spillTex;
  const [c, x] = canvas2d(128, 256);
  const gv = x.createLinearGradient(0, 0, 0, 126);
  gv.addColorStop(0, "rgba(255,255,255,1)");
  gv.addColorStop(0.3, "rgba(255,255,255,0.5)");
  gv.addColorStop(1, "rgba(255,255,255,0)");
  x.fillStyle = gv;
  x.fillRect(0, 0, 128, 126);
  x.globalCompositeOperation = "destination-in";
  const gh = x.createLinearGradient(0, 0, 128, 0);
  gh.addColorStop(0, "rgba(0,0,0,0)");
  gh.addColorStop(0.18, "rgba(0,0,0,1)");
  gh.addColorStop(0.82, "rgba(0,0,0,1)");
  gh.addColorStop(1, "rgba(0,0,0,0)");
  x.fillStyle = gh;
  x.fillRect(0, 0, 128, 126);
  x.globalCompositeOperation = "source-over";
  const gr = x.createRadialGradient(64, 192, 0, 64, 192, 62);
  gr.addColorStop(0, "rgba(255,255,255,1)");
  gr.addColorStop(0.45, "rgba(255,255,255,0.55)");
  gr.addColorStop(1, "rgba(255,255,255,0)");
  x.fillStyle = gr;
  x.fillRect(0, 130, 128, 126);
  return (spillTex = canvasTex(c, true));
}

/** A burner's light on its surroundings: additive, unlit, its strength set with the heat (`setSpill`). */
function spillMaterial(): THREE.MeshBasicMaterial {
  const m = new THREE.MeshBasicMaterial({ color: 0x000000, map: spillTexture(), transparent: true, depthWrite: false, blending: THREE.AdditiveBlending });
  m.userData.uv = "own";
  return m;
}
function setSpill(m: THREE.MeshBasicMaterial, heat: number, peak = 0.5): void {
  const k = smooth(0.12, 1, clamp01(heat));
  const a = peak * Math.pow(k, 1.3);
  m.color.setRGB(a, a * (0.2 + 0.14 * k), a * 0.035);
}

/** A plane for a spill map: `wall` uses the wash (bright along its top edge), otherwise the pool. */
function spillPlane(w: number, h: number, wall: boolean): THREE.BufferGeometry {
  const g = new THREE.PlaneGeometry(w, h);
  const uv = g.attributes.uv as THREE.BufferAttribute;
  for (let i = 0; i < uv.count; i++) uv.setY(i, wall ? 0.5 + uv.getY(i) * 0.5 : uv.getY(i) * 0.5);
  return g;
}

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
  // 18 gauge skin; the top's corners formed over a radius.
  const t = 0.0012;
  const rT = 0.008;
  const [cx0, cx1, cy0, cy1, cz1] = [0.032, w - 0.124, 0.035, 0.36, back + 0.03];
  // The valves at the right of the column, the rack lever's slot at its left beside them.
  const colX = cx1 + 0.074;
  const slotX = cx1 + 0.018;
  const yLow = 0.08;
  const travel = 0.2;
  const [cerMap, cerGlow] = ceramicMaps();
  const glow = heatMaterial(0xc4b9a6, cerGlow, cerMap, 0.9, 0);
  const spill = s.fine(1) ? spillMaterial() : null;

  // The shell: top and sides one wrapper formed over the corners, a bottom and back between its sides, and the
  // front one panel cut round the broiler's mouth, the lever's notched slot and the valves.
  const panel = new THREE.Shape();
  panel.moveTo(0, 0);
  panel.lineTo(w, 0);
  panel.lineTo(w, top - rT);
  panel.quadraticCurveTo(w, top, w - rT, top);
  panel.lineTo(rT, top);
  panel.quadraticCurveTo(0, top, 0, top - rT);
  panel.lineTo(0, 0);
  const sy0 = yLow - 0.012 - 0.007;
  const sy1 = yLow + travel - 0.012 + 0.007;
  if (s.fine(1)) {
    const u = new THREE.Shape();
    u.moveTo(0, 0);
    u.lineTo(t, 0);
    u.lineTo(t, top - rT);
    u.quadraticCurveTo(t, top - t, rT, top - t);
    u.lineTo(w - rT, top - t);
    u.quadraticCurveTo(w - t, top - t, w - t, top - rT);
    u.lineTo(w - t, 0);
    u.lineTo(w, 0);
    u.lineTo(w, top - rT);
    u.quadraticCurveTo(w, top, w - rT, top);
    u.lineTo(rT, top);
    u.quadraticCurveTo(0, top, 0, top - rT);
    u.lineTo(0, 0);
    s.add(p.ss, toCreasedNormals(plate(s, u, zf - t - back, 0, 6), 0.5), [0, 0, back]);
    panel.holes.push(rectPath(new THREE.Path(), cx1 - cx0, cy1 - cy0, 0.004, (cx0 + cx1) / 2, (cy0 + cy1) / 2, true));
    // The slot, the lever's four rests notched off it toward the valves.
    const sl: THREE.Vector2[] = [new THREE.Vector2(slotX - 0.0055, sy0), new THREE.Vector2(slotX + 0.0055, sy0)];
    for (let k = 0; k < 4; k++) {
      const ny = yLow - 0.012 + (k * travel) / 3;
      sl.push(new THREE.Vector2(slotX + 0.0055, ny - 0.0055), new THREE.Vector2(slotX + 0.017, ny - 0.0055), new THREE.Vector2(slotX + 0.017, ny + 0.0055), new THREE.Vector2(slotX + 0.0055, ny + 0.0055));
    }
    sl.push(new THREE.Vector2(slotX + 0.0055, sy1), new THREE.Vector2(slotX - 0.0055, sy1));
    panel.holes.push(new THREE.Path(sl));
    s.add(p.ss, toCreasedNormals(plate(s, panel, t, 0, 6), 0.5), [0, 0, zf - t]);
    // Dark behind the slot: the column's inside.
    s.add(p.enamel, box(0.04, sy1 - sy0 + 0.02, 0.001), [slotX + 0.006, (sy0 + sy1) / 2, zf - 0.03]);
  } else {
    for (const x of [t / 2, w - t / 2]) s.add(p.ss, box(t, top, zf - t - back), [x, top / 2, (zf - t + back) / 2]);
    s.add(p.ss, box(w - 2 * t, t, zf - t - back), [w / 2, top - t / 2, (zf - t + back) / 2]);
    s.add(p.ss, box(w, top - cy1, t), [w / 2, (top + cy1) / 2, zf - t / 2]);
    s.add(p.ss, box(cx0, cy1, t), [cx0 / 2, cy1 / 2, zf - t / 2]);
    s.add(p.ss, box(w - cx1, cy1, t), [(w + cx1) / 2, cy1 / 2, zf - t / 2]);
    s.add(p.ss, box(cx1 - cx0, cy0, t), [(cx0 + cx1) / 2, cy0 / 2, zf - t / 2]);
    s.add(p.enamel, box(0.011, sy1 - sy0, 0.0006), [slotX, (sy0 + sy1) / 2, zf + 0.0003]);
  }
  s.add(p.ss, box(w - 2 * t, t, zf - t - back), [w / 2, t / 2, (zf - t + back) / 2]);
  // The back, set in behind the hanger channels.
  s.add(p.ssDull, box(w - 2 * t, top - 2 * t, t), [w / 2, top / 2, back + 0.0041 + t / 2]);
  // The mouth's top edge hemmed back under itself.
  s.add(p.ss, cyl(s, 0.005, cx1 - cx0 - 0.002, "x", 16), [(cx0 + cx1) / 2, cy1 + 0.001, zf + 0.0035]);
  if (s.fine(2)) for (const x of [cx0 + 0.001, cx1 - 0.001]) s.add(p.ss, ball(s, 0.005, 12), [x, cy1 + 0.001, zf + 0.0035], undefined, [0.2, 1, 1]);

  // The cavity, heat-tinted: walls, floor, ceiling and back, butted at their corners; the right wall slotted for
  // the lever.
  const cd = zf - t - cz1;
  const czm = (zf - t + cz1) / 2;
  s.add(p.ssHeat, box(t, cy1 - cy0, cd), [cx0 - t / 2, (cy0 + cy1) / 2, czm]);
  s.add(p.ssHeat, box(t, cy1 - cy0, cd), [cx1 + t / 2, (cy0 + cy1) / 2, czm]);
  s.add(p.ssHeat, box(cx1 - cx0, t, cd), [(cx0 + cx1) / 2, cy0 - t / 2, czm]);
  s.add(p.ssHeat, box(cx1 - cx0, cy1 - cy0, t), [(cx0 + cx1) / 2, (cy0 + cy1) / 2, cz1 - t / 2]);
  if (s.fine(1)) {
    s.add(p.ssHeat, box(cx1 - cx0, t, cd), [(cx0 + cx1) / 2, cy1 + t / 2, czm]);
    s.add(p.enamel, box(0.0006, sy1 - sy0, 0.011), [cx1 - 0.0003, (sy0 + sy1) / 2, -0.12]);
    // The carriage's guide channels down each side.
    for (const x of [cx0, cx1])
      for (const z of [zf - 0.03, cz1 + 0.03]) s.add(p.ssHeat, s.fine(2) ? rbox(s, 0.008, cy1 - cy0 - 0.03, 0.014, 0.001) : box(0.008, cy1 - cy0 - 0.03, 0.014), [x + (x === cx0 ? 0.004 : -0.004), (cy0 + cy1) / 2, z]);
  }

  // The drip pan on the floor: a drawn pan, its front under the mouth with a pull.
  const pz0 = zf - t - 0.006;
  const pz1 = cz1 + 0.004;
  const pw = cx1 - cx0 - 0.012;
  const ph = 0.021;
  const panY = cy0 + 0.0015;
  if (s.fine(2)) s.add(p.ssDull, drawnPan(s, pw, pz0 - pz1, ph, 0.01, 0.006, 0.0012, 0.0025), [(cx0 + cx1) / 2, panY + ph, (pz0 + pz1) / 2]);
  else {
    s.add(p.ssDull, box(pw, ph, 0.0012), [(cx0 + cx1) / 2, panY + ph / 2, pz0 - 0.0006]);
    s.add(p.ssDull, box(pw, 0.0012, pz0 - pz1), [(cx0 + cx1) / 2, panY + 0.0006, (pz0 + pz1) / 2]);
  }
  pullHandle0(s, [(cx0 + cx1) / 2, panY + ph * 0.55, pz0], 0.12);

  // The burners: two infrared plaques in stainless boxes across the ceiling, pitched a little to face the mouth.
  // Each face is ceramic tiles of fine ports in a frame turned down round it.
  const bwid = (cx1 - cx0 - 0.03) / 2;
  const bdep = 0.28;
  const bz = (zf + cz1) / 2 - 0.02;
  const pitch = -0.1;
  for (let i = 0; i < 2; i++) {
    const bx = cx0 + 0.01 + bwid / 2 + i * (bwid + 0.01);
    s.at([bx, cy1 - 0.028, bz], [pitch, 0, 0], () => {
      s.add(p.ssHeat, box(bwid, 0.026, bdep), [0, 0.013, 0]);
      if (s.fine(1)) {
        const fr = roundRect(bwid + 0.008, bdep + 0.008, 0.002);
        fr.holes.push(rectPath(new THREE.Path(), bwid - 0.012, bdep - 0.012, 0.0015, 0, 0, true));
        s.add(p.ss, plate(s, fr, 0.008, 0, 3).rotateX(Math.PI / 2), [0, 0.003, 0]);
      }
      s.add(glow, box(bwid - 0.012, 0.004, bdep - 0.012), [0, -0.0005, 0]);
    });
  }

  // The burners' light on the cavity: down the walls from the ceiling, and a pool on the drip pan.
  if (spill) {
    const ch = cy1 - cy0 - 0.004;
    s.add(spill, spillPlane(cx1 - cx0 - 0.004, ch, true), [(cx0 + cx1) / 2, cy0 + 0.002 + ch / 2, cz1 + 0.002]);
    for (const [x, ry] of [[cx0 + 0.002, Math.PI / 2], [cx1 - 0.002, -Math.PI / 2]] as const) s.add(spill, spillPlane(cd - 0.01, ch, true), [x, cy0 + 0.002 + ch / 2, czm - 0.003], [0, ry, 0]);
    s.add(spill, spillPlane(pw - 0.03, pz0 - pz1 - 0.03, false), [(cx0 + cx1) / 2, panY + 0.0032, (pz0 + pz1) / 2], [-Math.PI / 2, 0, 0]);
  }

  // The control column: two burner valves over their printed dials, and the carriage lever's slot.
  const ky = [0.335, 0.27] as const;
  const knobs = ky.map((y, i) => valveKnob(s, `knob_${i}`, [colX, y, zf]));
  for (const y of ky) dialPlate(s, "gas", [colX, y, zf], 0.03);
  if (s.fine(2)) for (const y of [0.018, 0.382]) for (const x of [0.016, w - 0.016]) screw(s, p.ssBright, [x, y, zf], 0.0026);
  ratingPlate(s, "caution", [cx0 + 0.06, 0.3835, zf], 0.05, 0.025);

  // Flue on the top at the back, louvred; louvres for the cooling air down each side; the rating plate.
  s.add(p.ss, sheet(s, w - 0.2, h - top, 0.11, 0.001), [w / 2, (top + h) / 2 - 0.0008, back + 0.07]);
  louvres(s, [w / 2, (top + h) / 2, back + 0.125], w - 0.26, h - top - 0.012, 3);
  for (const [x, ry] of [[0, -Math.PI / 2], [w, Math.PI / 2]] as const) louvres(s, [x, 0.2, back + 0.15], 0.14, 0.12, 6, undefined, [0, ry, 0]);
  ratingPlate(s, "gas", [w, 0.31, back + 0.15], 0.07, 0.044, [0, Math.PI / 2, 0]);
  // The range-mount brackets on the back: hanger channels with their slots.
  for (const x of [0.16, w - 0.16]) {
    s.add(p.ssDull, s.fine(2) ? rbox(s, 0.05, 0.3, 0.0035, 0.0008) : box(0.05, 0.3, 0.0035), [x, 0.2, back + 0.0006 + 0.00175]);
    if (s.fine(1)) s.add(p.enamel, box(0.012, 0.05, 0.0006), [x, 0.31, back + 0.0003]);
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

  // The grid draws out to its stops at 175 mm, about half its depth. heroKitchen.ts sala() moves the cook's hand on the
  // pull by the same 0.175 m, so change both together.
  const draw = 0.175;
  return (b) => {
    let rack = 0.5;
    let out = 0;
    const pose = () => {
      const j = b.joints.get("broiler_rack_lift");
      if (j) j.position.set(0, yLow + rack * travel, 0);
      const g = b.joints.get("broiler_grid");
      if (g) g.position.set(0, 0, out * draw);
    };
    pose();
    return {
      glows: [glow],
      open: (f) => {
        out = clamp01(f);
        pose();
      },
      set: (name, v) => {
        if (name === "heat") {
          setCeramicHeat(glow, v);
          if (spill) setSpill(spill, v, 0.42);
          for (const k of knobs) turnJoint(b.joints.get(k), "z", valveAngle(v));
        } else if (name === "rack") {
          rack = Math.round(clamp01(v) * 3) / 3;
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
  dialPlate(s, "power", at, 0.032);
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
  ratingPlate(s, "elec", [w * 0.25, (caseTop + footH) / 2, back], 0.06, 0.03, [0, Math.PI, 0]);
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
  // The sides stop under the deck, the back and bottom between them, so no two skins share a face.
  const sk = 0.0016;
  const sideTop = deckY - sk;
  for (const sx of [sk / 2, w - sk / 2]) s.add(p.ss, sheet(s, sk, sideTop - legH, bodyD, 0.0006), [sx, (sideTop + legH) / 2, (zf + back) / 2]);
  s.add(p.ssDull, box(w - 2 * sk, sideTop - legH - sk, sk), [w / 2, (sideTop + legH + sk) / 2, back + sk / 2]);
  s.add(p.ssDull, box(w - 2 * sk, sk, bodyD), [w / 2, legH + sk / 2, (zf + back) / 2]);
  s.add(p.ss, sheet(s, w, sk, bodyD + 0.002, 0.0006), [w / 2, deckY - sk / 2, (zf + back) / 2 + 0.001]);
  s.add(p.ss, rbox(s, w - 2 * sk, sideTop - panY0, 0.012, 0.003), [w / 2, (sideTop + panY0) / 2, zf - 0.006]);
  s.add(p.ss, box(w - 2 * sk, doorY0 - legH - sk, sk), [w / 2, (doorY0 + legH + sk) / 2, zf - 0.008]);
  louvres(s, [w / 2, (doorY0 + legH) / 2, zf - 0.0072], w - 0.2, 0.04, 4);
  for (const sx of [sk + 0.008, w - sk - 0.008]) s.add(p.ss, sheet(s, 0.016, doorY1 - doorY0 + 0.006, sk, 0.0006), [sx, (doorY1 + doorY0) / 2, zf - sk / 2]);
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
  ratingPlate(s, "elec", [w, 0.6, -0.3], 0.08, 0.04, [0, Math.PI / 2, 0]);
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
    if (s.fine(1)) {
      // A double-pan door: the face drawn from one sheet, its edges radiused and returned 18 mm, the return's raw
      // edge hemmed inward; behind it a liner pan inset, its flange on the hem.
      j.add(p.ss, drawnPan(j, dw, dh, 0.018, 0.003, 0.0025, 0.001, 0).rotateX(-Math.PI / 2), [(sx * dw) / 2, dh / 2, 0]);
      if (s.fine(2)) {
        const hem = roundRect(dw, dh, 0.003, 0, 0);
        hem.holes.push(rectPath(new THREE.Path(), dw - 0.024, dh - 0.024, 0.003, 0, 0, true));
        j.add(p.ss, plate(j, hem, 0.001, 0, 4), [(sx * dw) / 2, dh / 2, -0.0004]);
        j.add(p.ssDull, drawnPan(j, dw - 0.03, dh - 0.03, 0.008, 0.004, 0.003, 0.001, 0).rotateX(Math.PI / 2), [(sx * dw) / 2, dh / 2, -0.0004]);
      } else j.add(p.ssDull, box(dw - 0.006, dh - 0.006, 0.006), [(sx * dw) / 2, dh / 2, -0.0026]);
      for (const y of [0.05, dh - 0.05]) j.add(p.ssBright, cyl(j, 0.005, 0.05, "y", 12), [0, y, 0.004]);
    } else j.add(p.ss, box(dw, dh, 0.018), [(sx * dw) / 2, dh / 2, 0.009]);
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
  // The side rails butt against the front and back ones' inner faces.
  for (const x of [0.05, 1.04]) s.add(p.ss, rbox(s, sq, sq, -0.12 - back - 0.08 - sq, 0.004), [x, frameY, (-0.12 + back + 0.08) / 2]);
  s.add(p.ss, cyl(s, 0.016, -0.12 - back - 0.08, "z", 14), [0.05, 0.2, (-0.12 + back + 0.08) / 2]);
  // The left trunnion's pillow block on its upright.
  s.add(p.ss, rbox(s, sq, piv[1] - frameY - 0.04, sq, 0.004), [0.035, (piv[1] + frameY) / 2 - 0.02, piv[2]]);
  s.add(p.castIron, rbox(s, 0.04, 0.05, 0.11, 0.008), [0.035, piv[1] - 0.012, piv[2]]);
  s.add(p.castIron, cyl(s, 0.042, 0.04, "x", 24), [0.035, piv[1], piv[2]]);
  if (s.fine(2)) s.add(p.brass, cyl(s, 0.004, 0.012, "y", 8), [0.035, piv[1] + 0.046, piv[2]]);
  // Its bolts' nuts, under the block.
  if (s.fine(2)) for (const sz of [-1, 1]) s.add(p.ssBright, new THREE.CylinderGeometry(0.008, 0.008, 0.008, 6), [0.035, piv[1] - 0.0405, piv[2] + sz * 0.042]);

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
  dialPlate(s, "temp450", [fx, 0.925, cz0], 0.034);
  s.add(lamp, cyl(s, 0.006, 0.006, "z", 14), [fx - 0.035, 0.86, cz0 + 0.003]);
  s.add(p.ssBright, torus(s, 0.007, 0.0014, "z", Math.PI * 2, 14, 5), [fx - 0.035, 0.86, cz0 + 0.001]);
  s.add(p.ssBright, cyl(s, 0.011, 0.006, "z", 16), [fx + 0.03, 0.86, cz0 + 0.003]);
  s.add(p.bakelite, cyl(s, 0.003, 0.02, "z", 8), [fx + 0.03, 0.865, cz0 + 0.014], [0.45, 0, 0]);
  ratingPlate(s, "elec", [cx1, 0.5, (cz0 + cz1) / 2], 0.08, 0.04, [0, Math.PI / 2, 0]);
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

  // The pan, on its trunnions: one drawn tub, its plan corners and its floor coved, the front wall raked forward and
  // pressed out at the middle of the rim into the pouring lip, which dips as it leaves; a jacket round it, and a rolled
  // cap over both.
  const pan = s.joint("tilt_pan", piv);
  inFixtureFrame(pan, piv, () => {
    const g = 0.025;
    const rcp = 0.04;
    const rb = 0.025;
    const rake = 0.035;
    const zFb = -0.14;
    const xm = (ix0 + ix1) / 2;
    const zMid = (iz1 + zFb) / 2;
    const [hw, lipOut, lipH, lipHo, dip] = [0.15, 0.055, 0.1, 0.035, 0.012];
    const cn = pan.seg(6, 1);
    const nf = pan.fine(2) ? 26 : pan.fine(1) ? 9 : 5;
    const frontZ = (y: number) => zFb + ((y - bed) / (rim - bed)) * rake;
    const shape = (x: number) => 1 - smooth(0.25, 1, Math.abs(x - xm) / hw);
    const level = (x0: number, x1: number, z0: number, z1: number, r: number, out: number): P2[] => {
      const rr = Math.max(1e-4, Math.min(r, (x1 - x0) / 2 - 1e-4, (z1 - z0) / 2 - 1e-4));
      const front: P2[] = [];
      for (let i = 1; i < nf; i++) {
        const x = x0 + rr + ((x1 - x0 - 2 * rr) * i) / nf;
        front.push([x, z1 + out * shape(x)]);
      }
      return [...front, ...loopRR(x0, x1, z0, z1, r, cn)];
    };
    // How far a point dips below its level: only the lip's, most at the rim and none a lip's height down.
    const dipAt = (x: number, z: number, dy: number) => (z > zMid ? dip * shape(x) * Math.max(0, 1 - dy / lipH) : 0);
    type Level = { loop: P2[]; y: number; dy: number };
    // A loft through the levels, rim first, its UVs in metres (round the plan, and down) so the grain runs on
    // unbroken round the corners; the first column repeated at the end to carry the seam.
    const loft = (levels: readonly Level[], flip: boolean) => {
      const n = levels[0]!.loop.length;
      const rows: THREE.Vector3[][] = [];
      for (let i = 0; i <= n; i++) rows.push(levels.map((l) => {
        const [x, z] = l.loop[i % n]!;
        return new THREE.Vector3(x, l.y - dipAt(x, z, l.dy), z);
      }));
      const geo = gridMesh(rows, false, false, flip);
      const uv = geo.attributes.uv as THREE.BufferAttribute;
      const run = levels.map(() => 0);
      const sc = (p.ss.userData.uvScale as number | undefined) ?? 1;
      for (let i = 0; i <= n; i++)
        for (let j = 0; j < levels.length; j++) {
          if (i > 0) run[j] += rows[i]![j]!.distanceTo(rows[i - 1]![j]!);
          uv.setXY(i * levels.length + j, run[j]! * sc, rows[i]![j]!.y * sc);
        }
      return geo;
    };
    const wallYs = pan.fine(2) ? [0, 0.012, 0.025, 0.04, 0.055, 0.07, 0.085, 0.1, 0.13] : pan.fine(1) ? [0, 0.04, 0.1] : [0, 0.1];
    // Inside: down the walls from the rim, round the cove onto the bed.
    const inner: Level[] = wallYs.map((dy) => {
      const y = rim - dy;
      return { y, dy, loop: level(ix0, ix1, iz1, frontZ(y), rcp, lipOut * smooth(rim - lipH, rim, y)) };
    });
    const cs = pan.seg(6, 1);
    for (let k = 0; k <= cs; k++) {
      const a = (k / cs) * (Math.PI / 2);
      const y = bed + rb - rb * Math.sin(a);
      const ins = rb * (1 - Math.cos(a));
      inner.push({ y, dy: rim - y, loop: level(ix0 + ins, ix1 - ins, iz1 + ins, frontZ(y) - ins, rcp - ins, 0) });
    }
    const ssOwn = ownGrain(p.ss);
    pan.add(ssOwn, loft(inner, false));
    pan.add(p.ss, loopCap(inner[inner.length - 1]!.loop, bed));
    // Outside: the jacket, the same plan grown by the gap, down to its floor round a radius; under the lip it turns
    // out sharply, so the lip stands proud of the jacket's face.
    const ob = 0.02;
    const outer: Level[] = [...wallYs.filter((dy) => dy < 0.06), 0.06, rim - oBot - ob].map((dy) => {
      const y = rim - dy;
      return { y, dy, loop: level(ix0 - g, ix1 + g, iz1 - g, frontZ(y) + g, rcp + g, lipOut * smooth(rim - lipHo, rim, y)) };
    });
    for (let k = 1; k <= cs; k++) {
      const a = (k / cs) * (Math.PI / 2);
      const y = oBot + ob - ob * Math.sin(a);
      const ins = ob * (1 - Math.cos(a));
      outer.push({ y, dy: rim - y, loop: level(ix0 - g + ins, ix1 + g - ins, iz1 - g + ins, frontZ(y) + g - ins, rcp + g - ins, 0) });
    }
    pan.add(ssOwn, loft(outer, true));
    pan.add(p.ssDull, loopCap(outer[outer.length - 1]!.loop, oBot, 0, true));
    // The rim: a flat cap over both walls, its edges rolled down over each, dipping with the lip.
    const rimLoop = inner[0]!.loop;
    const ra = pan.seg(5, 1);
    const cap: P2[] = [[-g - 0.001, -0.006], ...arc2([-g + 0.004, 0.003], 0.005, Math.PI, Math.PI / 2, ra), ...arc2([-0.0032, 0.004], 0.004, Math.PI / 2, 0, ra), [0.0008, -0.006]];
    const nrm = loopNormals(rimLoop);
    const capRows = rimLoop.map(([x, z], i) => cap.map(([ins, y]) => new THREE.Vector3(x + nrm[i]![0] * ins, rim + y - dipAt(x, z, 0), z + nrm[i]![1] * ins)));
    pan.add(p.ss, gridMesh(capRows, true, false, false));
    // Gallon marks etched on the back wall.
    if (pan.fine(2)) for (let k = 1; k <= 4; k++) pan.add(p.charcoal, box(k % 2 ? 0.05 : 0.08, 0.0012, 0.0004), [ix0 + rcp + 0.05, bed + k * 0.045, iz1 + 0.0003]);
    // Trunnion bosses and stubs.
    pan.add(p.ss, rbox(pan, 0.012, 0.1, 0.1, 0.005), [ox0 - 0.006, piv[1], piv[2]]);
    pan.add(p.ss, rbox(pan, 0.012, 0.1, 0.1, 0.005), [ox1 + 0.006, piv[1], piv[2]]);
    pan.add(p.ssBright, cyl(pan, 0.026, ox0 - 0.012 - 0.022, "x", 20), [(ox0 - 0.012 + 0.022) / 2, piv[1], piv[2]]);
    pan.add(p.ssBright, cyl(pan, 0.026, cx0 - 0.03 - ox1 - 0.012, "x", 20), [(ox1 + 0.012 + cx0 - 0.03) / 2, piv[1], piv[2]]);
    // The lid's hinge brackets on the pan's back.
    for (const x of [0.2, 0.91]) pan.add(p.ss, rbox(pan, 0.012, 0.05, 0.03, 0.003), [x, rim + 0.02, oz1 - 0.012]);
  });

  // The lid, hinged at the pan's back: a pressed cover with a raised panel, its edge turned down and beaded, a
  // tubular handle across the front, a steam vent; its counterweight on arms behind the hinge.
  const hinge: V3 = [0.555, rim + 0.031, oz1 - 0.022];
  let lid = pan;
  inFixtureFrame(pan, piv, () => (lid = pan.joint("lid_hinge", hinge)));
  inFixtureFrame(lid, hinge, () => {
    const [lz0, lz1] = [-0.07, oz1 - 0.006];
    const lw = ox1 - ox0 + 0.01;
    const ld = lz0 - lz1;
    const ly = rim + 0.0085;
    const lh = 0.026;
    const top = ly + lh;
    const lx = (ox0 + ox1) / 2;
    const lzm = (lz0 + lz1) / 2;
    if (lid.fine(1)) {
      lid.add(p.ss, drawnPan(lid, lw, ld, lh, 0.03, 0.01, 0.0012, 0.004).rotateX(Math.PI), [lx, ly, lzm]);
      if (lid.fine(2)) lid.add(p.ss, drawnPan(lid, lw - 0.16, ld - 0.16, 0.006, 0.03, 0.006, 0.001, 0).rotateX(Math.PI), [lx, top + 0.0004, lzm]);
    } else lid.add(p.ss, box(lw, lh, ld), [lx, ly + lh / 2, lzm]);
    const hy = 0.958;
    const hz = -0.086;
    lid.add(p.ssBright, pipe(lid, [[0.25, top, -0.104], [0.25, hy, hz], [0.86, hy, hz], [0.86, top, -0.104]], 0.011, 0.022, 14));
    for (const x of [0.25, 0.86]) lid.add(p.ssBright, cyl(lid, 0.017, 0.004, "y", 16), [x, top + 0.002, -0.104]);
    lid.add(p.ss, lathe(lid, [[0.04, 0], [0.038, 0.012], [0.03, 0.022], [0, 0.026]], 28), [0.555, top, -0.66]);
    if (lid.fine(2))
      for (let k = 0; k < 8; k++) {
        const a = (k / 8) * Math.PI * 2;
        lid.add(p.enamel, box(0.012, 0.002, 0.003), [0.555 + Math.cos(a) * 0.03, top + 0.016, -0.66 + Math.sin(a) * 0.03], [0, -a, -0.5]);
      }
    // Hinge knuckles and pin; counterweight arms and their weight, a bar with domed ends.
    for (const x of [0.2, 0.91]) {
      lid.add(p.ssBright, cyl(lid, 0.011, 0.04, "x", 16), [x - 0.026, hinge[1], hinge[2]]);
      lid.add(p.ss, rbox(lid, 0.01, 0.03, 0.15, 0.003), [x + 0.02, hinge[1], hinge[2] - 0.05]);
    }
    lid.add(p.ssBright, cyl(lid, 0.004, 0.75, "x", 8), [0.555, hinge[1], hinge[2]]);
    lid.add(p.ss, cyl(lid, 0.028, 0.75, "x", 24), [0.555, hinge[1], hinge[2] - 0.12]);
    if (lid.fine(1)) for (const [x, sx] of [[0.555 - 0.375, -1], [0.555 + 0.375, 1]] as const) lid.add(p.ss, ball(lid, 0.028, 18), [x, hinge[1], hinge[2] - 0.12], undefined, [sx * 0.25, 1, 1]);
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
  const ssOwn = ownGrain(p.ss);
  const uvScale = (p.ss.userData.uvScale as number | undefined) ?? 1;
  // The bowl and the wall above the jacket, one drawn shell: round in plan but at the front, where the wall's top is
  // pressed out and down into the pouring lip. Points round the plan crowd toward the lip so it stays smooth.
  const nr = s.fine(2) ? 96 : s.fine(1) ? 32 : 12;
  const warp = s.fine(2) ? 1.35 : s.fine(1) ? 1.45 : 1;
  const phis = Array.from({ length: nr }, (_, i) => {
    const u = -1 + (2 * i) / nr;
    return Math.PI * Math.sign(u) * Math.pow(Math.abs(u), warp);
  });
  const [lipOut, lipH, lipW, dip] = [0.06, 0.09, 0.36, 0.01];
  // Across the lip: a rounded V, pointed at the middle, easing into the round at its sides.
  const lipShape = (phi: number, wide = 1) => {
    // The tip rounded over a few degrees so the V carries no crease down the wall.
    const u = Math.min(1, Math.hypot(phi / (lipW * wide), 0.07) - 0.07);
    return Math.pow(1 - u, 1.6) * (1 - smooth(0.7, 1, u));
  };
  const wallTop = rimY - 0.006;
  // A ring of the shell at height y, radius r, swelling into the lip by how far up the wall it is; dipping too.
  const ring = (r: number, y: number, lip: boolean): THREE.Vector3[] =>
    phis.map((phi) => {
      const up = lip ? smooth(rimY - lipH, rimY, y) : 0;
      const rr = r + lipOut * up * lipShape(phi);
      const yy = y - (lip ? dip * lipShape(phi, 1.6) * smooth(rimY - lipH, rimY, y) : 0);
      return new THREE.Vector3(cx + Math.sin(phi) * rr, yy, cz + Math.cos(phi) * rr);
    });
  // A shell through rings (top first), UVs in metres round and down it, the first column repeated for the seam.
  const shell = (rings: THREE.Vector3[][], flip: boolean): THREE.BufferGeometry => {
    const rows = phis.map((_, i) => rings.map((rg) => rg[i]!));
    rows.push(rings.map((rg) => rg[0]!.clone()));
    const g = gridMesh(rows, false, false, flip);
    const uv = g.attributes.uv as THREE.BufferAttribute;
    const nv = rings.length;
    const down = rows.map(() => 0);
    for (let i = 0; i < rows.length; i++)
      for (let j = 0; j < nv; j++) {
        if (j > 0) down[i] += rows[i]![j]!.distanceTo(rows[i]![j - 1]!);
        const around = (i / nr) * Math.PI * 2 * R;
        uv.setXY(i * nv + j, around * uvScale, -down[i]! * uvScale);
      }
    return g;
  };
  const wallYs = s.fine(2) ? [0, 0.012, 0.025, 0.04, 0.055, 0.07, 0.09, 0.12] : s.fine(1) ? [0, 0.03, 0.06, 0.09] : [0, 0.09];
  // The inside, seen from within: down the wall, round the hemisphere to the bottom.
  const inner: THREE.Vector3[][] = [...wallYs.map((dy) => ring(R, wallTop - dy, true)), ...arc(R, Math.PI / 2, 0, k, false).slice(1).map(([r, y]) => ring(Math.max(r, 1e-4), y, false))];
  s.add(ssOwn, shell(inner, false));
  // The wall's outside above the jacket, and the jacket: its band welded to the bowl's wall, round under the bowl.
  const outerYs = wallYs.map((dy) => wallTop - dy).filter((y) => y > jTop + 0.006);
  s.add(ssOwn, shell([...outerYs, jTop + 0.004].map((y) => ring(R + 0.0025, y, true)), true));
  const jp = [...arc(Ro, 0, Math.PI / 2, k, false), [Ro, jTop - 0.012], [R + 0.03, jTop], [R + 0.003, jTop + 0.004]] as [number, number][];
  const jg = lathe(s, jp, n);
  {
    // Own UVs in metres: round the jacket, and down its profile.
    const uv = jg.attributes.uv as THREE.BufferAttribute;
    const run = [0];
    for (let i = 1; i < jp.length; i++) run.push(run[i - 1]! + Math.hypot(jp[i]![0] - jp[i - 1]![0], jp[i]![1] - jp[i - 1]![1]));
    for (let i = 0; i < uv.count; i++) {
      uv.setXY(i, uv.getX(i) * Math.PI * 2 * Ro * uvScale, run[i % jp.length]! * uvScale);
    }
  }
  s.add(ssOwn, jg, [cx, 0, cz]);
  // The rim rolled over into a bead all round, following the lip down and out.
  if (s.fine(1)) {
    const sec = ring2(-0.009, -0.0055, 0.0095, s.seg(10, 6), 0);
    const top = ring(R, rimY, true);
    const rows = top.map((pt, i) => {
      const phi = phis[i]!;
      const [nx, nz] = [-Math.sin(phi), -Math.cos(phi)];
      return sec.map(([ins, y]) => new THREE.Vector3(pt.x + nx * ins, pt.y + y, pt.z + nz * ins));
    });
    rows.push(rows[0]!.map((v) => v.clone()));
    const bg = gridMesh(rows, false, true, true);
    const buv = bg.attributes.uv as THREE.BufferAttribute;
    for (let i = 0; i < buv.count; i++) buv.setXY(i, buv.getX(i) * Math.PI * 2 * (R + 0.009) * uvScale, buv.getY(i) * Math.PI * 2 * 0.0095 * uvScale);
    s.add(ssOwn, bg);
    s.add(p.ssBright, torus(s, Ro - 0.004, 0.004, "y", Math.PI * 2, n, 6), [cx, jTop - 0.012, cz]);
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
  const [sx, sz] = [cx + Math.sin(sa) * (Ro - 0.012), cz + Math.cos(sa) * (Ro - 0.012)];
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
  dialPlate(s, "steps", [kxm, ky1 - 0.08, kz0], 0.034);
  for (const [i, lx] of [[0, kxm - 0.04], [1, kxm + 0.04]] as const) {
    s.add(i ? p.ledGreen : lamp, cyl(s, 0.006, 0.006, "z", 14), [lx, ky1 - 0.165, kz0 + 0.003]);
    if (s.fine(1)) s.add(p.ssBright, torus(s, 0.007, 0.0014, "z", Math.PI * 2, 14, 5), [lx, ky1 - 0.165, kz0 + 0.001]);
  }
  s.add(p.ssBright, cyl(s, 0.012, 0.006, "z", 16), [kxm, ky1 - 0.23, kz0 + 0.003]);
  s.add(p.bakelite, cyl(s, 0.003, 0.02, "z", 8), [kxm, ky1 - 0.225, kz0 + 0.014], [0.45, 0, 0]);
  ratingPlate(s, "elec", [kxm, ky0 + 0.07, kz0], 0.08, 0.04);
  // Its supply comes up from the floor in a conduit.
  if (s.fine(1)) {
    s.add(p.ssDull, cyl(s, 0.012, ky0, "y", 10), [kxm, ky0 / 2, (kz0 + kz1) / 2]);
    s.add(p.ssDull, cyl(s, 0.024, 0.006, "y", 14), [kxm, 0.003, (kz0 + kz1) / 2]);
  }

  // The contents: a level surface whose radius follows the bowl, and its simmer.
  const contents = s.joint("contents", [cx, 0.6, cz]);
  contents.add(stock(), new THREE.CircleGeometry(0.985, s.seg(64, 16)).rotateX(-Math.PI / 2));
  // Its meniscus, climbing the wall (in the joint's unit plan, scaled to the bowl's radius at the level).
  if (s.fine(1)) contents.add(stock(), lathe(contents, [[1, 0.003], [0.996, 0.0016], [0.99, 0.0006], [0.985, 0]], s.seg(64, 16)));

  // The lid: a shallow dome hinged on a sprung bracket behind the rim, a loop handle at its front.
  const hinge: V3 = [cx, rimY + 0.012, cz - R - 0.06];
  const lid = s.joint("kettle_lid", hinge);
  inFixtureFrame(lid, hinge, () => {
    // A shallow dome, its edge rolled into a bead that rests on the rim's, so the lip runs out from under it.
    const lr = R + 0.022;
    const ly = rimY + 0.007;
    lid.add(p.ss, lathe(lid, [[lr - 0.004, 0.004], [lr * 0.75, 0.024], [lr * 0.4, 0.035], [0, 0.038]], n), [cx, ly, cz]);
    lid.add(p.ss, torus(lid, lr - 0.004, 0.004, "y", Math.PI * 2, n, 8), [cx, ly + 0.001, cz]);
    // Its underside, seen when it stands open.
    if (lid.fine(1)) lid.add(p.ss, lathe(lid, [[0, 0.034], [lr * 0.4, 0.031], [lr * 0.75, 0.02], [lr - 0.006, 0.0005]], n), [cx, ly, cz]);
    lid.add(p.ssBright, pipe(lid, [[cx - 0.07, ly + 0.012, cz + lr - 0.03], [cx - 0.07, ly + 0.04, cz + lr + 0.01], [cx + 0.07, ly + 0.04, cz + lr + 0.01], [cx + 0.07, ly + 0.012, cz + lr - 0.03]], 0.008, 0.02, 12));
    for (const sxx of [-1, 1]) lid.add(p.ss, rbox(lid, 0.014, 0.012, 0.09, 0.003), [cx + sxx * 0.06, ly + 0.012, cz - lr + 0.02], [0.2, 0, 0]);
    lid.add(p.ssBright, cyl(lid, 0.01, 0.15, "x", 14), [cx, hinge[1], hinge[2]]);
  });
  // The hinge bracket on the rim's back and its counterbalance spring.
  s.add(p.ss, rbox(s, 0.2, 0.05, 0.05, 0.005), [cx, rimY - 0.012, cz - R - 0.035]);
  if (s.fine(1)) s.add(p.ssBright, cyl(s, 0.018, 0.19, "x", 18), [cx, rimY + 0.012, cz - R - 0.06]);
  if (s.fine(2)) for (let i = 0; i < 9; i++) s.add(p.ssDull, torus(s, 0.0205, 0.0025, "x", Math.PI * 2, 14, 5), [cx - 0.09 + i * 0.0225, rimY + 0.012, cz - R - 0.06]);

  return (b) => {
    const bub = new Bubbles(stockFoam(), s.fine(2) ? 50 : 20, [cx, 0.6, cz], 0.5, 0.5, 0.005, true);
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

/** A sheathed element's glow along its length: cold for its last few centimetres into each terminal. */
let elementTex: THREE.CanvasTexture | null = null;
function elementGlowMap(): THREE.CanvasTexture {
  if (elementTex) return elementTex;
  const [c, x] = canvas2d(8, 512);
  const g = x.createLinearGradient(0, 0, 0, 512);
  g.addColorStop(0, "#000000");
  g.addColorStop(0.012, "#000000");
  g.addColorStop(0.04, "#ffffff");
  g.addColorStop(0.96, "#ffffff");
  g.addColorStop(0.988, "#000000");
  g.addColorStop(1, "#000000");
  x.fillStyle = g;
  x.fillRect(0, 0, 8, 512);
  return (elementTex = canvasTex(c, true));
}

/**
 * A 48 in infrared strip warmer, hung under a shelf over the pass: an extruded aluminium housing, ribbed on top, its
 * lower edges rolled under into lips that hold a polished parabolic reflector; a sheathed element along the
 * reflector's focus on two stirrups, its terminals in ceramic bushings; formed end caps, the right one carrying the
 * toggle switch and the amber pilot, the cord leaving the left one through a strain relief and running up into the
 * shelf; two hat-section brackets on top bolt it up.
 */
const stripWarmer: Builder = (s, e) => {
  const p = palette();
  const { w, d, h } = e;
  const bodyH = 0.05;
  const cap = 0.018;
  const xc0 = 0.014;
  const [x0, x1] = [xc0 + cap, w - xc0 - cap];
  const glow = heatMaterial(0x34302d, elementGlowMap(), null, 0.5, 0.65);
  const spill = s.fine(1) ? spillMaterial() : null;
  const pilot = s.fine(2) ? lampMaterial(0xffa21a) : null;
  const refl = s.fine(1) ? reflectorMat() : p.ssBright;

  // The housing's section, [z, y]: from the front lip's inner edge, under it and round the rolled front edge, up the
  // front, over the ribbed top, down the back and under to the back lip's edge; its face outward throughout.
  const ee = 0.0015;
  const yb = 0.0012;
  const [zF, zB] = [-ee, -d + ee];
  const [rc, rr, lipW, tl] = [0.006, 0.003, 0.02, 0.0025];
  const ca = s.seg(6, 2);
  const sec: P2[] = [
    [zF - lipW, yb + tl],
    [zF - lipW, yb],
    [zF - lipW, yb],
    ...arc2([zF - rr, yb + rr], rr, -Math.PI / 2, 0, ca),
    ...arc2([zF - rc, bodyH - rc], rc, 0, Math.PI / 2, ca),
  ];
  if (s.fine(2)) {
    // Seven low ribs along the top, each a rounded bead off the flat.
    const n = 7;
    for (let k = 0; k < n; k++) {
      const zc = zF - rc - 0.012 - ((d - 2 * ee - 2 * rc - 0.024) * k) / (n - 1);
      sec.push([zc + 0.0034, bodyH], [zc + 0.0022, bodyH + 0.0016], [zc, bodyH + 0.0022], [zc - 0.0022, bodyH + 0.0016], [zc - 0.0034, bodyH]);
    }
  }
  sec.push(...arc2([zB + rc, bodyH - rc], rc, Math.PI / 2, Math.PI, ca), ...arc2([zB + rr, yb + rr], rr, Math.PI, Math.PI * 1.5, ca), [zB + lipW, yb], [zB + lipW, yb], [zB + lipW, yb + tl]);
  s.add(p.alu, sweepX(sec, x0, x1));

  // The reflector: a parabola from lip to lip, its focus on the element; flanged out over the lips that hold it.
  const zc = -d / 2;
  const yv = bodyH - 0.004;
  const ye = yb + tl + 0.0004;
  const a = d / 2 - ee - lipW;
  const f = (a * a) / (4 * (yv - ye));
  const [ey, ez] = [yv - f, zc];
  const pn = s.seg(28, 8);
  const para: P2[] = [];
  for (let i = 0; i <= pn; i++) {
    const z = zc - a + (2 * a * i) / pn;
    para.push([z, yv - ((z - zc) * (z - zc)) / (4 * f)]);
  }
  s.add(refl, sweepX([[zc - a - 0.004, ye], [zc - a, ye], ...para, [zc + a, ye], [zc + a + 0.004, ye]], x0, x1));
  // The element's image in the polished reflector: a warm wash strongest at the crown, fading toward the lips and
  // the ends.
  if (spill) {
    const g = sweepX(para.map(([z, y]) => [z, y - 0.0005] as P2), x0, x1);
    const uv = g.attributes.uv as THREE.BufferAttribute;
    for (let i = 0; i < uv.count; i++) uv.setXY(i, 0.14 + 0.72 * uv.getX(i), 1 - Math.min(0.5, Math.abs(2 * uv.getY(i) - 1) * 0.8));
    s.add(spill, g);
  }

  // The element: a sheathed tube from terminal to terminal, the ends in ceramic bushings, on two stirrups.
  s.add(glow, cyl(s, 0.0042, x1 - x0 - 0.004, "x", 14), [w / 2, ey, ez]);
  for (const x of [x0 + 0.007, x1 - 0.007]) {
    s.add(p.ceramic, s.fine(1) ? lathe(s, [[0.0046, 0], [0.0085, 0], [0.0095, 0.002], [0.0095, 0.012], [0.0085, 0.014], [0.0046, 0.014]], 16).rotateZ(x < w / 2 ? -Math.PI / 2 : Math.PI / 2) : cyl(s, 0.0095, 0.014, "x", 8), [x < w / 2 ? x - 0.007 : x + 0.007, ey, ez]);
  }
  if (s.fine(1))
    for (const fx of [0.33, 0.67]) {
      const x = x0 + (x1 - x0) * fx;
      s.add(p.ssBright, pipe(s, [[x, yv - 0.0005, ez], [x, ey + 0.0058, ez]], 0.0011, 0.0011, 6));
      s.add(p.ssBright, torus(s, 0.0052, 0.0011, "x", Math.PI * 2, 18, 6), [x, ey, ez]);
    }

  // End caps: formed to the section, standing a little proud of it all round, screwed to the extrusion.
  for (const [cx, sx] of [[xc0, -1], [w - xc0 - cap, 1]] as const) {
    const ch = bodyH + 0.0028;
    if (s.fine(1)) {
      const g = plate(s, roundRect(d - 0.0024, ch - 0.0024, 0.0065, d / 2, ch / 2), cap, 0.0012, 6).rotateY(Math.PI / 2);
      s.add(p.ss, g, [cx, 0, 0]);
    } else s.add(p.ss, box(cap, ch, d), [cx + cap / 2, ch / 2, -d / 2]);
    if (s.fine(2))
      for (const [zz, yy] of [[-0.022, 0.012], [-d + 0.022, 0.012], [-0.022, ch - 0.012], [-d + 0.022, ch - 0.012]] as const)
        screw(s, p.ssBright, [sx < 0 ? cx : cx + cap, yy, zz], 0.0026, [0, (sx * Math.PI) / 2, 0]);
  }

  // The right cap: a switch plate with the toggle in its threaded bushing and hex nut, the pilot beside it.
  const xo = w - xc0;
  const [tz, lz, my] = [-d / 2 + 0.016, -d / 2 - 0.018, bodyH / 2];
  s.add(p.ss, rbox(s, 0.0025, 0.032, 0.07, 0.0012), [xo + 0.00125, my, -d / 2]);
  s.add(p.ssBright, cyl(s, 0.0042, 0.004, "x", 14), [xo + 0.0045, my, tz]);
  if (s.fine(1)) s.add(p.ssBright, new THREE.CylinderGeometry(0.0068, 0.0068, 0.0022, 6).rotateZ(Math.PI / 2), [xo + 0.0036, my, tz]);
  const sw = s.joint("toggle", [xo + 0.0065, my, tz]);
  sw.add(p.ssBright, cyl(sw, 0.0017, 0.0062, "x", 8, 0.0012).translate(0.0031, 0, 0), [0, 0, 0], [0, 0, 0.35]);
  if (s.fine(2)) sw.add(p.ssBright, ball(sw, 0.0017, 10).translate(0.0062, 0, 0), [0, 0, 0], [0, 0, 0.35]);
  if (pilot) {
    s.add(p.ssBright, torus(s, 0.0058, 0.0012, "x", Math.PI * 2, 20, 6), [xo + 0.0028, my, lz]);
    s.add(pilot, lathe(s, [[0.0052, 0], [0.0052, 0.0008], [0.0042, 0.0028], [0.0026, 0.0042], [0, 0.0045]], 18).rotateZ(-Math.PI / 2), [xo + 0.0025, my, lz]);
  }

  // The cord: out of the left cap through a moulded strain relief, turning up into the shelf above.
  const cz = -d * 0.72;
  s.add(p.rubber, lathe(s, [[0, 0], [0.0072, 0], [0.0068, 0.0025], [0.0052, 0.006], [0.0042, 0.0085], [0, 0.0085]], 14).rotateZ(Math.PI / 2), [xc0, my, cz]);
  s.add(p.rubber, pipe(s, [[xc0 - 0.004, my, cz], [0.0041, my, cz], [0.0041, h, cz]], 0.0033, 0.0045, 10));
  ratingPlate(s, "elec120", [w * 0.7, my, zB], 0.06, 0.03, [0, Math.PI, 0]);

  // Mounting brackets: hat sections on top, their crowns slotted for the bolts.
  for (const x of [0.22, w - 0.22]) {
    if (s.fine(2)) {
      const crown = roundRect(0.05, d - 0.02, 0.002, 0, 0);
      for (const z of [0.045, d - 0.045]) crown.holes.push(rectPath(new THREE.Path(), 0.0085, 0.022, 0.0042, 0, z - d / 2, true));
      s.add(p.ss, xzPlate(s, crown, 0.002), [x, h - 0.002, -d / 2]);
    } else s.add(p.ss, sheet(s, 0.05, 0.002, d - 0.02, 0.0006), [x, h - 0.001, -d / 2]);
    for (const z of [-0.016, -d + 0.016]) s.add(p.ss, sheet(s, 0.05, h - bodyH - 0.002, 0.002, 0.0006), [x, (h - 0.002 + bodyH) / 2, z]);
  }

  return (b) => ({
    glows: [glow],
    set: (name, v) => {
      if (name === "heat") {
        setRedHeat(glow, v, 1.05);
        if (spill) setSpill(spill, v, 0.22);
        if (pilot) pilot.emissiveIntensity = v > 0.05 ? 1.3 : 0;
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
