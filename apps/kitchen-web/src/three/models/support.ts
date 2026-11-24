/**
 * Model family: ventilation, sanitation and prep. Hoods, sinks, work tables, dishmachines, the insulated carrier,
 * the food processor, the slicer, the portion scale, the box grater and the cutting board. Proportions are those of the
 * equipment class (NSF/ANSI 2, 3, 8; UL 710 and UL 300 for hoods and their fire systems), sized to the archetype's
 * extents.
 *
 * Construction is drawn as it is made: hoods are 18-gauge double-wall canopies with welded skirts, a sloped bank of
 * interlocking S-baffle filters over a grease trough, vapour-proof globes, a wet-chemical system's piping, nozzles and
 * fusible links; sinks are drawn bowls with coved corners on 1-5/8 in legs and gussets; tables are 16-gauge tops with
 * bull-nose or marine edges on hat channels. Anything that moves in life is a joint.
 *
 * Joints, open() and set() names, by archetype:
 *   NSF-HOOD-TYPE1-*, NSF-HOOD-TYPE2-CONDENSATE: no joints. set("lights", 0..1) dims the globes (default lit).
 *   NSF-SINK-3COMP, NSF-VEG-PREP-SINK-48IN: joints "waste_<i>" (lever wastes, detail 1 up), "spout_<i>" (swing
 *     spouts), "water". set("water", 0..1) fills the bowls; set("drain", 0..1) turns the waste levers from shut to open.
 *   NSF-HANDWASH-SINK-WALL-17IN: joint "stream". set("flow", 0..1) runs the faucet.
 *   NSF-PREPTABLE-48/72/96: none.
 *   NSF-DISHWASHER-DOOR: joints "hood_lift" (slides up 0.44 m; open(frac)), "wash_arm_upper", "wash_arm_lower",
 *     "needle_wash", "needle_rinse", "needle_pressure". tick spins the arms while the hood is down;
 *     set("wash_temp" | "rinse_temp", °C, or K when above 200) moves the gauge needles.
 *   NSF-DISHWASHER-CONVEYOR: joints "inspect_door" (open(frac) swings it 90°), "conveyor_pawl", "curtain_<i>" (detail
 *     1 up), "needle_wash", "needle_rinse". tick strokes the pawl bar and sways the curtains with load.
 *   NSF-CARRIER-INSULATED: joints "carrier_door", "latch_<i>", "vent_cap", "pan_<i>". open(frac) throws the two latches
 *     over, then swings the door round to 180°; set("vent", 0..1) turns the cap from shut to open, its windows over
 *     the slots; set("pans", 0..6) loads that many lidded GN 1/1 pans, bottom up (default 3).
 *   NSF-PROCESSOR-SBLADE: joints "s_blade_spin", "pusher". tick spins the blade with load; open(frac) lifts the pusher.
 *   NSF-SLICER-ROTARY: joints "rotary_blade_spin", "carriage_slide", "gauge_plate", "index_knob". tick spins the blade
 *     and strokes the carriage with load; set("carriage", -1..1) places the carriage; set("thickness", 0..1) turns the
 *     index knob and opens the gauge plate.
 *   NSF-PORTION-SCALE-DIGI: set("weight", kg) shows the reading.
 *   NSF-GRATER-BOX, NSF-CUTTING-BOARD-HDPE-24X18: none.
 */
import * as THREE from "three";
import { RoundedBoxGeometry } from "three/examples/jsm/geometries/RoundedBoxGeometry.js";
import { register, type Builder, type Hooks } from "./index";
import { palette, perforation } from "./materials";
import { alongAxis, ball, cyl, knob, lathe, plate, pullHandle, rbox, rectPath, roundRect, screw, Shop, slideJoint, torus, turnJoint, type Axis, type Built, type V3 } from "./shop";

type Mat = THREE.Material;
type P2 = readonly [number, number];

/* ------------------------------------------------------------------------------------------------------------ */
/* Workshop helpers                                                                                               */
/* ------------------------------------------------------------------------------------------------------------ */

/**
 * A block from corner `a` to corner `b`, its edges rounded by `r` (0 for a plain box, 12 triangles). Roundings under
 * 4 mm are drawn only at full detail: from where detail 1 is seen they are below a pixel.
 */
export function blk(s: Shop, m: Mat, a: V3, b: V3, r = 0.0008): void {
  const w = Math.abs(b[0] - a[0]);
  const h = Math.abs(b[1] - a[1]);
  const d = Math.abs(b[2] - a[2]);
  const round = r > 0 && (s.detail === 2 || (s.detail === 1 && r >= 0.004));
  const g = round ? rbox(s, w, h, d, r) : new THREE.BoxGeometry(w, h, d);
  s.add(m, g, [(a[0] + b[0]) / 2, (a[1] + b[1]) / 2, (a[2] + b[2]) / 2]);
}

/** A plain `n`-sided rod for the far detail levels, `h` long along `axis` and centred at `at`; open-ended by default. */
export function rod(s: Shop, m: Mat, r: number, h: number, axis: Axis, at: V3, n = 6, open = true): void {
  s.add(m, alongAxis(new THREE.CylinderGeometry(r, r, h, n, 1, open), axis), at);
}

/**
 * A round bar or pipe through `pts`, each corner rounded at `bend`: shop's `bar`, but with its straight runs as single
 * spans and only its bends subdivided, the section carried round them by parallel transport so that the rings line up.
 */
export function pipe(s: Shop, pts: readonly V3[], r: number, bend = r * 3, radial = 12): THREE.BufferGeometry {
  const v = pts.map((q) => new THREE.Vector3(q[0], q[1], q[2]));
  const line: THREE.Vector3[] = [v[0]!.clone()];
  const nb = s.seg(8, 3);
  for (let i = 1; i < v.length - 1; i++) {
    const [a0, c0, b0] = [v[i - 1]!, v[i]!, v[i + 1]!];
    const cut = Math.min(bend, c0.distanceTo(a0) * 0.49, b0.distanceTo(c0) * 0.49);
    const a = c0.clone().addScaledVector(c0.clone().sub(a0).normalize(), -cut);
    const b = c0.clone().addScaledVector(b0.clone().sub(c0).normalize(), cut);
    const curve = new THREE.QuadraticBezierCurve3(a, c0.clone(), b);
    for (let k = 0; k <= nb; k++) line.push(curve.getPoint(k / nb));
  }
  line.push(v[v.length - 1]!.clone());
  const c = line.filter((q, i) => i === 0 || q.distanceToSquared(line[i - 1]!) > 1e-12);
  const m = c.length;
  const n = s.seg(radial, 5);
  const T = c.map((_, i) => c[Math.min(m - 1, i + 1)]!.clone().sub(c[Math.max(0, i - 1)]!).normalize());
  const N = new THREE.Vector3().crossVectors(T[0]!, Math.abs(T[0]!.x) < 0.9 ? new THREE.Vector3(1, 0, 0) : new THREE.Vector3(0, 1, 0)).normalize();
  const B = new THREE.Vector3();
  const q = new THREE.Quaternion();
  const pos: number[] = [];
  const nor: number[] = [];
  const uv: number[] = [];
  const idx: number[] = [];
  const rings: [THREE.Vector3, THREE.Vector3][] = [];
  let along = 0;
  for (let i = 0; i < m; i++) {
    if (i > 0) {
      N.applyQuaternion(q.setFromUnitVectors(T[i - 1]!, T[i]!));
      along += c[i]!.distanceTo(c[i - 1]!);
    }
    N.addScaledVector(T[i]!, -N.dot(T[i]!)).normalize();
    B.crossVectors(T[i]!, N);
    rings.push([N.clone(), B.clone()]);
    for (let k = 0; k <= n; k++) {
      const a = (k / n) * Math.PI * 2;
      const [cx, sy] = [Math.cos(a), Math.sin(a)];
      const nx = N.x * cx + B.x * sy;
      const ny = N.y * cx + B.y * sy;
      const nz = N.z * cx + B.z * sy;
      pos.push(c[i]!.x + r * nx, c[i]!.y + r * ny, c[i]!.z + r * nz);
      nor.push(nx, ny, nz);
      uv.push(along, k / n);
    }
  }
  for (let i = 0; i < m - 1; i++)
    for (let k = 0; k < n; k++) {
      const a = i * (n + 1) + k;
      const b = a + n + 1;
      idx.push(a, a + 1, b, b, a + 1, b + 1);
    }
  // End caps, as bar has them, from detail 1.
  if (s.detail > 0)
    for (const [i, sign] of [
      [0, -1],
      [m - 1, 1],
    ] as const) {
      const base = pos.length / 3;
      const [Ni, Bi] = rings[i]!;
      const t = T[i]!;
      pos.push(c[i]!.x, c[i]!.y, c[i]!.z);
      nor.push(t.x * sign, t.y * sign, t.z * sign);
      uv.push(0.5, 0.5);
      for (let k = 0; k <= n; k++) {
        const a = (k / n) * Math.PI * 2;
        const [cx, sy] = [Math.cos(a), Math.sin(a)];
        pos.push(c[i]!.x + r * (Ni.x * cx + Bi.x * sy), c[i]!.y + r * (Ni.y * cx + Bi.y * sy), c[i]!.z + r * (Ni.z * cx + Bi.z * sy));
        nor.push(t.x * sign, t.y * sign, t.z * sign);
        uv.push(0.5 + cx / 2, 0.5 + sy / 2);
      }
      for (let k = 0; k < n; k++) idx.push(...(sign > 0 ? [base, base + 1 + k, base + 2 + k] : [base, base + 2 + k, base + 1 + k]));
    }
  const g = new THREE.BufferGeometry();
  g.setAttribute("position", new THREE.Float32BufferAttribute(pos, 3));
  g.setAttribute("normal", new THREE.Float32BufferAttribute(nor, 3));
  g.setAttribute("uv", new THREE.Float32BufferAttribute(uv, 2));
  g.setIndex(idx);
  return g;
}

/** Points on an arc about (cx, cy), from angle a0 to a1 inclusive, in n steps. */
export function arc(cx: number, cy: number, r: number, a0: number, a1: number, n: number): [number, number][] {
  const out: [number, number][] = [];
  for (let i = 0; i <= n; i++) {
    const a = a0 + ((a1 - a0) * i) / n;
    out.push([cx + r * Math.cos(a), cy + r * Math.sin(a)]);
  }
  return out;
}

/** A sheet's section: the closed outline of a centreline polyline thickened to `t`, mitred at the bends. */
export function sheetProfile(pts: readonly P2[], t: number): THREE.Shape {
  const n = pts.length;
  const L: THREE.Vector2[] = [];
  const R: THREE.Vector2[] = [];
  for (let i = 0; i < n; i++) {
    const p = pts[i]!;
    const a = pts[Math.max(0, i - 1)]!;
    const b = pts[Math.min(n - 1, i + 1)]!;
    let tx = b[0] - a[0];
    let ty = b[1] - a[1];
    const l = Math.hypot(tx, ty) || 1;
    tx /= l;
    ty /= l;
    let k = 1;
    if (i > 0 && i < n - 1) {
      const u = [p[0] - a[0], p[1] - a[1]];
      const v = [b[0] - p[0], b[1] - p[1]];
      const lu = Math.hypot(u[0]!, u[1]!) || 1;
      const lv = Math.hypot(v[0]!, v[1]!) || 1;
      const c = (u[0]! * v[0]! + u[1]! * v[1]!) / (lu * lv);
      k = 1 / Math.max(0.4, Math.sqrt((1 + c) / 2));
    }
    const [nx, ny] = [-ty * (t / 2) * k, tx * (t / 2) * k];
    L.push(new THREE.Vector2(p[0] + nx, p[1] + ny));
    R.push(new THREE.Vector2(p[0] - nx, p[1] - ny));
  }
  return new THREE.Shape([...L, ...R.reverse()]);
}

/**
 * Smooth normals across shallow bends only, matching vertices to 0.01 mm (the stock helper matches to 1 cm, which
 * fuses the two faces of a sheet). Extruded sections come out faceted without this.
 */
export function creased(g: THREE.BufferGeometry, angle = 0.62): THREE.BufferGeometry {
  const src = g.index ? g.toNonIndexed() : g;
  const pos = src.attributes.position as THREE.BufferAttribute;
  const n = pos.count;
  const fn = new Float32Array(n * 3);
  const map = new Map<string, number[]>();
  const a = new THREE.Vector3();
  const b = new THREE.Vector3();
  const c = new THREE.Vector3();
  const key = (i: number) => `${Math.round(pos.getX(i) * 1e5)},${Math.round(pos.getY(i) * 1e5)},${Math.round(pos.getZ(i) * 1e5)}`;
  for (let i = 0; i < n; i += 3) {
    a.fromBufferAttribute(pos, i);
    b.fromBufferAttribute(pos, i + 1);
    c.fromBufferAttribute(pos, i + 2);
    const nrm = c.sub(b).cross(a.sub(b)).normalize();
    for (let k = 0; k < 3; k++) {
      fn[(i + k) * 3] = nrm.x;
      fn[(i + k) * 3 + 1] = nrm.y;
      fn[(i + k) * 3 + 2] = nrm.z;
      const kk = key(i + k);
      let list = map.get(kk);
      if (!list) map.set(kk, (list = []));
      list.push(i + k);
    }
  }
  const out = new Float32Array(n * 3);
  const cos = Math.cos(angle);
  for (let i = 0; i < n; i++) {
    let [x, y, z] = [0, 0, 0];
    for (const j of map.get(key(i))!) {
      const d = fn[i * 3]! * fn[j * 3]! + fn[i * 3 + 1]! * fn[j * 3 + 1]! + fn[i * 3 + 2]! * fn[j * 3 + 2]!;
      if (d >= cos) {
        x += fn[j * 3]!;
        y += fn[j * 3 + 1]!;
        z += fn[j * 3 + 2]!;
      }
    }
    const l = Math.hypot(x, y, z) || 1;
    out[i * 3] = x / l;
    out[i * 3 + 1] = y / l;
    out[i * 3 + 2] = z / l;
  }
  src.setAttribute("normal", new THREE.BufferAttribute(out, 3));
  return src;
}

/**
 * A folded-sheet section `t` thick extruded `len` along `axis` from 0. The profile is its centreline in the plane
 * across the axis: (z, y) for "x", (x, z) for "y", (x, y) for "z".
 */
export function sweep(s: Shop, pts: readonly P2[], t: number, len: number, axis: Axis = "x"): THREE.BufferGeometry {
  return section(s, sheetProfile(axis === "z" ? pts : pts.map(([u, v]) => (axis === "x" ? [-u, v] : [u, -v]) as P2), t), len, axis, true);
}
/** A solid section (a shape in the same plane as `sweep`'s profile, already mapped) extruded along `axis`. */
function section(s: Shop, shape: THREE.Shape, len: number, axis: Axis, mapped = false): THREE.BufferGeometry {
  let sh = shape;
  if (!mapped && axis !== "z") sh = new THREE.Shape(shape.getPoints().map((q) => (axis === "x" ? new THREE.Vector2(-q.x, q.y) : new THREE.Vector2(q.x, -q.y))));
  let g = plate(s, sh, len);
  if (s.detail > 0) g = creased(g);
  if (axis === "x") g.rotateY(Math.PI / 2);
  else if (axis === "y") g.rotateX(-Math.PI / 2);
  return g;
}
/** Add a folded-sheet section run along x from x0 to x1, its (z, y) profile offset by (y, z). */
export function sweepX(s: Shop, m: Mat, zy: readonly P2[], t: number, x0: number, x1: number, y = 0, z = 0): void {
  s.add(m, sweep(s, zy, t, x1 - x0, "x"), [x0, y, z]);
}
/** Add a folded-sheet section run along z from z0 to z1 (z0 < z1), its (x, y) profile offset by (x, y). */
export function sweepZ(s: Shop, m: Mat, xy: readonly P2[], t: number, z0: number, z1: number, x = 0, y = 0): void {
  s.add(m, sweep(s, xy, t, z1 - z0, "z"), [x, y, z0]);
}
/** Add a folded-sheet section run up y from y0 to y1, its (x, z) profile offset by (x, z). */
export function sweepY(s: Shop, m: Mat, xz: readonly P2[], t: number, y0: number, y1: number, x = 0, z = 0): void {
  s.add(m, sweep(s, xz, t, y1 - y0, "y"), [x, y0, z]);
}

/** A lathe whose profile may be given top-down (hanging parts): ordered bottom-up so its normals face out. */
export function hang(s: Shop, prof: readonly (readonly [number, number])[], n = 48): THREE.BufferGeometry {
  const pr = prof[0]![1] > prof[prof.length - 1]![1] ? [...prof].reverse() : prof;
  return lathe(s, pr, n);
}

/** Drop every triangle whose centroid lies above `y` (in the geometry's own frame): opens a rounded box's top. */
function openTop(g: THREE.BufferGeometry, y: number): THREE.BufferGeometry {
  return cutAbove(g, "y", y);
}
/** Drop every triangle whose centroid lies beyond `v` along `axis` (in the geometry's own frame). */
function cutAbove(g: THREE.BufferGeometry, axis: Axis, v: number): THREE.BufferGeometry {
  const src = g.index ? g.toNonIndexed() : g;
  const p = src.attributes.position.array as ArrayLike<number>;
  const n = src.attributes.normal.array as ArrayLike<number>;
  const k0 = axis === "x" ? 0 : axis === "y" ? 1 : 2;
  const kp: number[] = [];
  const kn: number[] = [];
  for (let i = 0; i < p.length; i += 9) {
    if ((p[i + k0]! + p[i + 3 + k0]! + p[i + 6 + k0]!) / 3 > v) continue;
    for (let k = 0; k < 9; k++) {
      kp.push(p[i + k]!);
      kn.push(n[i + k]!);
    }
  }
  const out = new THREE.BufferGeometry();
  out.setAttribute("position", new THREE.Float32BufferAttribute(kp, 3));
  out.setAttribute("normal", new THREE.Float32BufferAttribute(kn, 3));
  return out;
}
/** Turn a non-indexed surface inside out: winding reversed, normals negated. */
function inward(g: THREE.BufferGeometry): THREE.BufferGeometry {
  for (const name of ["position", "normal"]) {
    const a = g.attributes[name]!.array as Float32Array;
    for (let i = 0; i < a.length; i += 9)
      for (let k = 0; k < 3; k++) {
        const t = a[i + 3 + k]!;
        a[i + 3 + k] = a[i + 6 + k]!;
        a[i + 6 + k] = t;
      }
    if (name === "normal") for (let i = 0; i < a.length; i++) a[i] = -a[i]!;
  }
  return g;
}

/**
 * A drawn sink bowl, open at its rim: the inner and outer skins of a rounded box, its corners coved at `r`. Inner size
 * `w` across, `dz` front to back, `depth` deep, centred at (cx, cz), its rim at `rimY`.
 */
export function basin(s: Shop, m: Mat, cx: number, cz: number, w: number, dz: number, depth: number, rimY: number, r = 0.025, t = 0.0016): void {
  const seg = s.detail === 2 ? 3 : 1;
  for (const [ww, dd, hh, rr, flip] of [
    [w, dz, depth, r, true],
    [w + 2 * t, dz + 2 * t, depth + t, r + t, false],
  ] as const) {
    if (s.detail === 0) {
      const g = openTop(new THREE.BoxGeometry(ww, hh, dd), hh / 2 - 1e-4);
      s.add(m, flip ? inward(g) : g, [cx, rimY - hh / 2, cz]);
      continue;
    }
    const H = hh + rr;
    const g = openTop(new RoundedBoxGeometry(ww, H, dd, seg, rr), H / 2 - rr + 1e-5);
    s.add(m, flip ? inward(g) : g, [cx, rimY - hh + H / 2, cz]);
  }
}

/** A flat sheet over [x0,x1]×[z0,z1] between heights y0 and y1, cut round rectangular openings [x0, x1, z0, z1]. */
export function sheetWithHoles(s: Shop, m: Mat, x0: number, x1: number, z0: number, z1: number, y0: number, y1: number, holes: readonly (readonly [number, number, number, number])[]): void {
  const xs = [...new Set([x0, x1, ...holes.flatMap((q) => [q[0], q[1]])])].filter((x) => x >= x0 && x <= x1).sort((a, b) => a - b);
  for (let i = 0; i < xs.length - 1; i++) {
    const [a, b] = [xs[i]!, xs[i + 1]!];
    if (b - a < 1e-6) continue;
    const mid = (a + b) / 2;
    const cuts = holes.filter((q) => q[0] < mid && q[1] > mid).sort((p, q) => p[2] - q[2]);
    let z = z0;
    for (const c of cuts) {
      if (c[2] > z + 1e-6) blk(s, m, [a, y0, z], [b, y1, c[2]], 0);
      z = Math.max(z, c[3]);
    }
    if (z1 > z + 1e-6) blk(s, m, [a, y0, z], [b, y1, z1], 0);
  }
}

/** Parts built once and stamped many times (filters, plates, pans): cloning is far cheaper than building. */
interface Part {
  m: Mat;
  g: THREE.BufferGeometry;
  at?: V3;
  rot?: V3;
}
function stamp(s: Shop, parts: readonly Part[], at: V3, rot: V3 = [0, 0, 0]): void {
  s.at(at, rot, () => {
    for (const q of parts) s.add(q.m, q.g.clone(), q.at, q.rot);
  });
}

/** A hex body (nut, union, nozzle) along `axis`. */
function hex(r: number, h: number, axis: Axis = "y"): THREE.BufferGeometry {
  const g = new THREE.CylinderGeometry(r, r, h, 6);
  if (axis === "x") g.rotateZ(-Math.PI / 2);
  else if (axis === "z") g.rotateX(Math.PI / 2);
  return g;
}

/**
 * A blank data plate (no maker's name): a satin plate with rivets and rows of engraved lines, `w` by `h`, its face
 * toward +z of the frame at `at` turned by `rot`.
 */
export function dataPlate(s: Shop, at: V3, rot: V3, w = 0.09, h = 0.055): void {
  if (!s.fine(1)) return;
  const p = palette();
  s.at(at, rot, () => {
    s.add(p.ssBright, rbox(s, w, h, 0.0008, 0.0003), [0, 0, 0.0004]);
    if (!s.fine(2)) return;
    for (let i = 0; i < 5; i++) {
      const lw = (i === 0 ? 0.55 : 0.8 - (i % 2) * 0.25) * w;
      s.add(p.enamel, new THREE.BoxGeometry(lw, h * 0.05, 0.0003), [-w / 2 + 0.008 + lw / 2, h / 2 - 0.012 - i * h * 0.16, 0.00092]);
    }
    for (const sx of [-1, 1]) for (const sy of [-1, 1]) s.add(p.ssBright, ball(s, 0.0016, 8, Math.PI / 2), [sx * (w / 2 - 0.004), sy * (h / 2 - 0.004), 0.0008], [Math.PI / 2, 0, 0]);
  });
}

/** Perforated stainless sheet, alpha-tested round holes on `pitch` centres. Cached: one material per pattern. */
const perfCache = new Map<string, THREE.MeshStandardMaterial>();
export function perforatedSheet(pitch: number, frac = 0.55, stagger = true): THREE.MeshStandardMaterial {
  const key = `${pitch}|${frac}|${stagger}`;
  const hit = perfCache.get(key);
  if (hit) return hit;
  const holes = 32;
  const m = new THREE.MeshStandardMaterial({ color: 0xc4c8cd, metalness: 1, roughness: 0.38, alphaMap: perforation(holes, frac, stagger), alphaTest: 0.5, side: THREE.DoubleSide });
  m.alphaToCoverage = true;
  m.userData.uv = "box";
  m.userData.uvScale = 1 / (pitch * holes);
  perfCache.set(key, m);
  return m;
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Legs, feet and edges shared by sinks, tables and machines                                                       */
/* ------------------------------------------------------------------------------------------------------------ */

/** 1-5/8 in stainless tube: the leg every NSF table and sink stands on. */
const LEG_R = 0.0206;

/**
 * An adjustable bullet foot: a stainless bullet threaded into the leg's foot, an inch of thread showing, its domed
 * pad on the floor at y = 0.
 */
export function bulletFoot(s: Shop, x: number, z: number, r = LEG_R): number {
  const p = palette();
  const top = 0.05;
  if (s.detail === 0) return top;
  const prof: [number, number][] = s.fine(2)
    ? [[0, 0], [r * 0.55, 0.0005], [r * 0.9, 0.006], [r * 1.02, 0.018], [r * 1.02, 0.034], [r * 0.78, 0.04], [0, 0.04]]
    : [[0, 0], [r * 0.95, 0.008], [r * 1.02, 0.034], [0, 0.04]];
  s.add(p.ssBright, lathe(s, prof, s.fine(2) ? 22 : 16), [x, 0, z]);
  // The thread between bullet and leg, then the leg's swaged end.
  s.add(p.ssDull, cyl(s, r * 0.42, top - 0.038, "y", 12, r * 0.42, true), [x, (0.038 + top) / 2, z]);
  return top;
}

/**
 * A table or sink leg at (x, z): bullet foot, 1-5/8 in tube up to `y1`, and its gusset, the flared socket welded under
 * the top that the leg slides into and is pinned by a set screw. Returns nothing; the gusset's flange sits at y1.
 */
export function gussetLeg(s: Shop, x: number, z: number, y1: number, gusset = true): void {
  const p = palette();
  if (s.detail === 0) {
    rod(s, p.ss, LEG_R, y1, "y", [x, y1 / 2, z]);
    return;
  }
  const foot = bulletFoot(s, x, z);
  const gh = gusset ? 0.075 : 0;
  const n = s.fine(2) ? 22 : 16;
  s.add(p.ss, cyl(s, LEG_R, y1 - gh - foot + 0.004, "y", n, LEG_R, true), [x, (foot + y1 - gh) / 2, z]);
  if (!gusset) return;
  const prof: [number, number][] = s.fine(2)
    ? [[LEG_R + 0.0012, 0], [LEG_R + 0.0022, 0.004], [LEG_R + 0.003, 0.04], [0.031, 0.066], [0.036, 0.0715], [0.036, gh], [0, gh]]
    : [[LEG_R + 0.002, 0], [LEG_R + 0.003, 0.04], [0.036, 0.07], [0, gh]];
  s.add(p.ss, lathe(s, prof, n), [x, y1 - gh, z]);
  if (s.fine(2)) {
    s.add(p.ssBright, hex(0.0042, 0.006, "z"), [x, y1 - gh + 0.022, z + LEG_R + 0.004]);
    s.add(p.ssBright, cyl(s, 0.002, 0.004, "z", 8), [x, y1 - gh + 0.022, z + LEG_R + 0.008]);
  }
}

/** A 1 in round cross-rail between two leg positions at height y, with its welded collars. */
function crossRail(s: Shop, a: V3, b: V3, r = 0.0127): void {
  const p = palette();
  const len = Math.hypot(b[0] - a[0], b[2] - a[2]);
  const axis: Axis = Math.abs(b[0] - a[0]) > Math.abs(b[2] - a[2]) ? "x" : "z";
  const mid: V3 = [(a[0] + b[0]) / 2, a[1], (a[2] + b[2]) / 2];
  if (s.detail === 0) return rod(s, p.ss, r, len - 2 * LEG_R, axis, mid);
  s.add(p.ss, cyl(s, r, len - 2 * LEG_R, axis, s.fine(2) ? 18 : 12, r, true), mid);
  if (s.fine(2)) for (const q of [a, b]) s.add(p.ss, cyl(s, r * 1.35, 0.012, axis, 18), [q[0] + (axis === "x" ? Math.sign(b[0] - a[0]) * (q === a ? 1 : -1) * (LEG_R + 0.004) : 0), q[1], q[2] + (axis === "z" ? Math.sign(b[2] - a[2]) * (q === a ? 1 : -1) * (LEG_R + 0.004) : 0)]);
}

/**
 * Edge profiles of 16-gauge tops, as (z, y) centrelines with the top's front edge at (0, 0) and the sheet running back
 * toward -z: a boxed bull-nose (1-5/8 in, rolled down and returned under), and a marine edge (a raised rolled lip that
 * keeps spills on the top, then down and under).
 */
function bullNose(s: Shop, back = 0.02): P2[] {
  const n = s.fine(2) ? s.seg(8, 3) : 2;
  return [[-back, 0], ...arc(-0.01, -0.01, 0.01, Math.PI / 2, 0, n), ...arc(-0.01, -0.031, 0.01, 0, -Math.PI / 2, n), [-0.03, -0.041], [-0.034, -0.037]];
}
function marineEdge(s: Shop, back = 0.03): P2[] {
  const n = s.fine(2) ? s.seg(10, 3) : 2;
  return [[-back, 0], [-0.022, 0], ...arc(-0.022, 0.012, 0.012, -Math.PI / 2, -0.25, n).slice(1), ...arc(-0.006, 0.006, 0.006, Math.PI * 0.72, 0, n), [0, -0.034], ...arc(-0.006, -0.034, 0.006, 0, -Math.PI / 2, 3), [-0.026, -0.04]];
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Formed sections swept round a plan                                                                             */
/* ------------------------------------------------------------------------------------------------------------ */

/**
 * A plan path of (x, z) points with each corner rounded on a true arc of radius `r` (one for all, or one per point) in
 * `n` steps; an open path keeps its two ends.
 */
export function filletPath(pts: readonly P2[], r: number | readonly number[], n: number, closed = false): P2[] {
  const m = pts.length;
  const out: P2[] = [];
  for (let i = 0; i < m; i++) {
    const q = pts[i]!;
    const ri = typeof r === "number" ? r : (r[i] ?? 0);
    if (ri <= 0 || n < 1 || (!closed && (i === 0 || i === m - 1))) {
      out.push(q);
      continue;
    }
    const a = pts[(i - 1 + m) % m]!;
    const b = pts[(i + 1) % m]!;
    const la = Math.hypot(q[0] - a[0], q[1] - a[1]);
    const lb = Math.hypot(b[0] - q[0], b[1] - q[1]);
    const da: P2 = [(q[0] - a[0]) / la, (q[1] - a[1]) / la];
    const db: P2 = [(b[0] - q[0]) / lb, (b[1] - q[1]) / lb];
    const cos = Math.max(-1, Math.min(1, da[0] * db[0] + da[1] * db[1]));
    const th = Math.acos(cos);
    if (th < 1e-4) {
      out.push(q);
      continue;
    }
    let rr = ri;
    let tl = rr * Math.tan(th / 2);
    const lim = Math.min(la, lb) * 0.49;
    if (tl > lim) {
      tl = lim;
      rr = tl / Math.tan(th / 2);
    }
    const s0: P2 = [q[0] - da[0] * tl, q[1] - da[1] * tl];
    const s1: P2 = [q[0] + db[0] * tl, q[1] + db[1] * tl];
    // The centre lies off s0 toward the turn.
    let [px, pz] = [db[0] - da[0] * cos, db[1] - da[1] * cos];
    const pl = Math.hypot(px, pz) || 1;
    px /= pl;
    pz /= pl;
    const c: P2 = [s0[0] + px * rr, s0[1] + pz * rr];
    const a0 = Math.atan2(s0[1] - c[1], s0[0] - c[0]);
    let a1 = Math.atan2(s1[1] - c[1], s1[0] - c[0]);
    while (a1 - a0 > Math.PI) a1 -= 2 * Math.PI;
    while (a1 - a0 < -Math.PI) a1 += 2 * Math.PI;
    for (let k = 0; k <= n; k++) {
      const t = a0 + ((a1 - a0) * k) / n;
      out.push([c[0] + rr * Math.cos(t), c[1] + rr * Math.sin(t)]);
    }
  }
  return out;
}

/** The closed outline of a sheet `t` thick drawn along a centreline, as (u, v) points for `pathSweep`. */
export function sheetOutline(pts: readonly P2[], t: number): P2[] {
  return sheetProfile(pts, t)
    .getPoints()
    .map((q) => [q.x, q.y] as P2);
}

/** Each point of a plan path with its mitred inward offset direction, scaled so that offsets stay parallel. */
function pathFrames(path: readonly P2[], closed: boolean): [number, number, number, number][] {
  const m = path.length;
  const segN = (i: number): P2 => {
    const a = path[i]!;
    const b = path[(i + 1) % m]!;
    const l = Math.hypot(b[0] - a[0], b[1] - a[1]) || 1;
    return [(b[1] - a[1]) / l, -(b[0] - a[0]) / l];
  };
  const out: [number, number, number, number][] = [];
  for (let i = 0; i < m; i++) {
    const n1 = closed || i > 0 ? segN((i - 1 + m) % m) : segN(i);
    const n2 = closed || i < m - 1 ? segN(i) : n1;
    let [nx, nz] = [n1[0] + n2[0], n1[1] + n2[1]];
    const l = Math.hypot(nx, nz) || 1;
    nx /= l;
    nz /= l;
    const sc = 1 / Math.max(0.3, nx * n1[0] + nz * n1[1]);
    out.push([path[i]![0], path[i]![1], nx * sc, nz * sc]);
  }
  return out;
}

/** A plan path moved `u` across itself, the way `pathSweep` carries its section. */
export function offsetPath(path: readonly P2[], u: number, closed = false): P2[] {
  return pathFrames(path, closed).map(([x, z, nx, nz]) => [x + nx * u, z + nz * u] as P2);
}

/**
 * A section swept along a plan path: the closed `outline` is (u, v), u measured across the path toward (t.z, -t.x)
 * of its direction of travel t, v up from `y0`. Corners are mitred, so a filleted path gives a formed corner and a sharp
 * one a welded mitre, with nothing lapped or doubled; an open path is capped at both ends. Normals are smoothed across
 * shallow bends only.
 */
export function pathSweep(outline: readonly P2[], path: readonly P2[], y0: number, closed = false): THREE.BufferGeometry {
  const m = path.length;
  const k = outline.length;
  const pos: number[] = [];
  for (const [x, z, nx, nz] of pathFrames(path, closed)) for (const [u, v] of outline) pos.push(x + nx * u, y0 + v, z + nz * u);
  let area = 0;
  for (let j = 0; j < k; j++) {
    const a = outline[j]!;
    const b = outline[(j + 1) % k]!;
    area += a[0] * b[1] - b[0] * a[1];
  }
  const ccw = area > 0;
  const idx: number[] = [];
  const rings = closed ? m : m - 1;
  for (let i = 0; i < rings; i++) {
    const i1 = (i + 1) % m;
    for (let j = 0; j < k; j++) {
      const j1 = (j + 1) % k;
      const [a0, b0, a1, b1] = [i * k + j, i * k + j1, i1 * k + j, i1 * k + j1];
      if (ccw) idx.push(a0, b0, b1, a0, b1, a1);
      else idx.push(a0, b1, b0, a0, a1, b1);
    }
  }
  if (!closed) {
    const tri = THREE.ShapeUtils.triangulateShape(outline.map(([u, v]) => new THREE.Vector2(u, v)), []);
    for (const [ring, facing] of [
      [0, -1],
      [m - 1, 1],
    ] as const)
      for (const [a, b, c] of tri) {
        const [pa, pb, pc] = [outline[a]!, outline[b]!, outline[c]!];
        const sgn = (pb[0] - pa[0]) * (pc[1] - pa[1]) - (pc[0] - pa[0]) * (pb[1] - pa[1]);
        // Counter-clockwise in (u, v) faces along the travel.
        if (sgn * facing > 0) idx.push(ring * k + a, ring * k + b, ring * k + c);
        else idx.push(ring * k + a, ring * k + c, ring * k + b);
      }
  }
  const g = new THREE.BufferGeometry();
  g.setAttribute("position", new THREE.Float32BufferAttribute(pos, 3));
  g.setIndex(idx);
  return creased(g);
}

/** A plan outline (x, z) as a shape for `plate` laid flat by rotateX(-π/2): shape y is -z. */
function planShape(pts: readonly P2[]): THREE.Shape {
  return new THREE.Shape(pts.map(([x, z]) => new THREE.Vector2(x, -z)));
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Hoods                                                                                                          */
/* ------------------------------------------------------------------------------------------------------------ */

/** The baffle bank's slope from the horizontal: steep enough to drain grease, shallow enough to catch it. */
const SLOPE = (55 * Math.PI) / 180;
/** Nominal 16 × 20 in baffle filters are 15-1/2 in across and 19-1/2 in up the slope, 2 in deep. */
const FW = 0.394;
const FH = 0.495;
const FT = 0.048;

/**
 * One baffle grease filter, built once in its own frame and stamped along the bank: x across (centred), y up its face
 * from the bottom edge, z out of the face (its back at -FT). A U-channel frame, its stiles butted between the rails,
 * holds two interlocking rows of chevron baffles: the front row's ridges toward the room, the rear row's toward the
 * plenum, each with its edges turned back into hooks, so the air turns twice through every gap and throws its grease
 * on the walls; it runs down the channels to weep holes in the bottom rail. Two wire handles lift it out.
 */
function baffleFilter(s: Shop): Part[] {
  const p = palette();
  const parts: Part[] = [];
  const fw = FW - 0.003;
  const hw = fw / 2;
  const fr = 0.025;
  const ft = -0.0012;
  if (s.detail === 0) {
    parts.push({ m: p.ss, g: new THREE.BoxGeometry(fw, FH - 0.004, FT * 0.6), at: [0, FH / 2, -FT * 0.3] });
    return parts;
  }
  const fine = s.fine(2);
  const rc = fine ? 0.003 : 0;
  // Frame: the rails run the full width, the stiles between them; flanges flush, nothing lapped at the corners.
  const stile: P2[] = fine ? [[fr, ft], [rc, ft], [0, ft - rc], [0, -FT + rc], [rc, -FT], [fr, -FT]] : [[fr, ft], [0, ft], [0, -FT], [fr, -FT]];
  for (const sx of [-1, 1]) {
    const g = sweep(s, stile.map(([x, z]) => [x * sx, z] as P2), 0.0012, FH - 2 * fr, "y");
    parts.push({ m: p.ss, g, at: [-sx * hw, fr, 0] });
  }
  const rail = (y0: number, dir: number): P2[] => (fine ? [[ft, y0 + dir * fr], [ft, y0 + dir * rc], [ft - rc, y0], [-FT + rc, y0], [-FT, y0 + dir * rc], [-FT, y0 + dir * fr]] : [[ft, y0 + dir * fr], [ft, y0], [-FT, y0], [-FT, y0 + dir * fr]]);
  for (const [y0, dir] of [
    [0, 1],
    [FH, -1],
  ] as const) {
    const g = sweep(s, rail(y0, dir), 0.0012, fw, "x");
    parts.push({ m: p.ss, g, at: [-hw, 0, 0] });
  }
  // The baffles, their ends in the rails' channels.
  const inner = fw - 2 * fr;
  const nf = 9;
  const pitch = inner / nf;
  const bl = FH - 2 * fr + 0.01;
  const by = fr - 0.005;
  if (fine) {
    const rw = pitch * 0.42;
    const front: P2[] = [[-rw + 0.0035, -0.031], [-rw, -0.0285], [-rw, -0.0135], [-0.0024, -0.0049], [0, -0.0042], [0.0024, -0.0049], [rw, -0.0135], [rw, -0.0285], [rw - 0.0035, -0.031]];
    const back: P2[] = front.map(([x, z]) => [x, -FT - z] as P2);
    const fg = sweep(s, front, 0.0008, bl, "y");
    const bg = sweep(s, back, 0.0008, bl, "y");
    for (let i = 0; i < nf; i++) parts.push({ m: p.ss, g: fg, at: [-inner / 2 + pitch * (i + 0.5), by, 0] });
    for (let i = 0; i <= nf; i++) parts.push({ m: p.ssHeat, g: bg, at: [-inner / 2 + pitch * i, by, 0] });
  } else {
    const fg = new THREE.BoxGeometry(pitch * 0.78, bl, 0.012);
    const bg = new THREE.BoxGeometry(pitch * 0.78, bl, 0.012);
    for (let i = 0; i < nf; i++) parts.push({ m: p.ss, g: fg, at: [-inner / 2 + pitch * (i + 0.5), by + bl / 2, -0.012] });
    for (let i = 0; i <= nf; i++) parts.push({ m: p.ssHeat, g: bg, at: [-inner / 2 + pitch * i, by + bl / 2, -0.034] });
  }
  if (fine) {
    // Wire lift handles on the bottom rail, weep holes in its floor, and the corner rivets.
    const hg = pipe(s, [[-0.035, 0.013, -0.001], [-0.035, 0.013, 0.02], [0.035, 0.013, 0.02], [0.035, 0.013, -0.001]], 0.003, 0.008, 8);
    for (const sx of [-1, 1]) {
      parts.push({ m: p.ssBright, g: hg, at: [sx * fw * 0.26, 0, 0] });
      for (const ex of [-0.035, 0.035]) parts.push({ m: p.ssBright, g: new THREE.CylinderGeometry(0.0055, 0.0055, 0.002, 12).rotateX(Math.PI / 2), at: [sx * fw * 0.26 + ex, 0.013, -0.0002] });
    }
    const weep = new THREE.CylinderGeometry(0.0028, 0.0028, 0.0024, 10);
    for (let i = 0; i < 6; i++) parts.push({ m: p.enamel, g: weep, at: [-inner / 2 + inner * ((i + 0.5) / 6), 0, -FT * 0.5] });
    const rivet = ball(s, 0.0022, 8, Math.PI / 2);
    rivet.rotateX(Math.PI / 2);
    for (const sx of [-1, 1]) for (const y of [0.011, FH - 0.011]) parts.push({ m: p.ssBright, g: rivet, at: [sx * (hw - 0.011), y, ft] });
  }
  return parts;
}

/** A soft radial falloff, full at the centre and clear at the rim: the pool a lamp throws on the faces round it. */
let poolTex: THREE.Texture | null = null;
function poolTexture(): THREE.Texture {
  if (poolTex) return poolTex;
  const c = document.createElement("canvas");
  c.width = c.height = 128;
  const x = c.getContext("2d")!;
  const g = x.createRadialGradient(64, 64, 0, 64, 64, 64);
  for (const [t, a] of [
    [0, 1],
    [0.16, 0.8],
    [0.38, 0.4],
    [0.62, 0.13],
    [0.84, 0.03],
    [1, 0],
  ] as const)
    g.addColorStop(t, `rgba(255,255,255,${a})`);
  x.fillStyle = g;
  x.fillRect(0, 0, 128, 128);
  poolTex = new THREE.CanvasTexture(c);
  return poolTex;
}

interface Lamps {
  lamp: THREE.MeshStandardMaterial;
  globe: THREE.MeshStandardMaterial;
  pool: THREE.MeshBasicMaterial;
  set(v: number): void;
}
/**
 * A hood's lights: the lamps; the frosted globes, lit from within and brightest where they face the eye, as frosted
 * glass is; and the warm pools they throw on the faces round them. Every hood builds with the one lit set, so a kitchen
 * of hoods draws its lights in three calls; `lampHooks` gives a hood a set of its own only when it is dimmed.
 */
let litLamps: Lamps | null = null;
function hoodLamps(): Lamps {
  return (litLamps ??= makeLamps());
}

/** The "lights" hook: lit (1) keeps the shared set; anything else moves this hood onto its own set and dims that. */
function lampHooks(b: Built): Hooks {
  let own: Lamps | null = null;
  return {
    set: (name, v) => {
      if (name !== "lights") return;
      if (!own) {
        if (v >= 0.999) return;
        const shared = hoodLamps();
        const mine = makeLamps();
        const swap = new Map<THREE.Material, THREE.Material>([
          [shared.lamp, mine.lamp],
          [shared.globe, mine.globe],
          [shared.pool, mine.pool],
        ]);
        b.group.traverse((o) => {
          const m = o as THREE.Mesh;
          if (!m.isMesh) return;
          const to = swap.get(m.material as THREE.Material);
          if (to) m.material = to;
        });
        own = mine;
      }
      own.set(v);
    },
  };
}

function makeLamps(): Lamps {
  const lamp = new THREE.MeshStandardMaterial({ color: 0xfff6ea, roughness: 0.5, metalness: 0, emissive: 0xffe8c4, emissiveIntensity: 3.2 });
  lamp.userData.uv = "own";
  const globe = new THREE.MeshStandardMaterial({ color: 0xe6e3dc, roughness: 0.34, metalness: 0, emissive: 0xffd8a4, emissiveIntensity: 1.15, transparent: true, opacity: 0.86 });
  globe.userData.uv = "own";
  globe.onBeforeCompile = (sh) => {
    sh.fragmentShader = sh.fragmentShader.replace("#include <emissivemap_fragment>", "#include <emissivemap_fragment>\n\ttotalEmissiveRadiance *= 0.28 + 0.72 * pow(abs(dot(normal, normalize(vViewPosition))), 1.6);");
  };
  globe.customProgramCacheKey = () => "hood-globe-frost";
  const pool = new THREE.MeshBasicMaterial({ color: 0xffd09a, map: poolTexture(), transparent: true, opacity: 0.42, blending: THREE.AdditiveBlending, depthWrite: false, polygonOffset: true, polygonOffsetFactor: -2, polygonOffsetUnits: -4 });
  pool.userData.uv = "own";
  const set = (v: number) => {
    const k = Math.max(0, Math.min(1, v));
    lamp.emissiveIntensity = 3.2 * k;
    globe.emissiveIntensity = 1.15 * k;
    pool.opacity = 0.42 * k;
    pool.visible = k > 0.01;
  };
  return { lamp, globe, pool, set };
}

/** A pool of lamplight on a face: a quad `w` by `h` centred at `at`, its front toward `rot` turned +z. */
function lightPool(s: Shop, L: Lamps, at: V3, w: number, h: number, rot: V3): void {
  if (!s.fine(1) || w <= 0.02 || h <= 0.02) return;
  s.add(L.pool, new THREE.PlaneGeometry(w, h), at, rot);
}

/**
 * A vapour-proof hood light hung from a ceiling whose underside is at y: the cast housing on its flange, the threaded
 * collar, a frosted jelly-jar globe with the lamp glowing in it, and the cast guard round the globe.
 */
function hoodLight(s: Shop, L: Lamps, x: number, y: number, z: number): void {
  const p = palette();
  s.add(p.alu, hang(s, [
    [0, 0],
    [0.064, 0],
    [0.064, -0.005],
    [0.061, -0.007],
    [0.058, -0.02],
    [0.047, -0.026],
    [0, -0.026],
  ], s.fine(2) ? 40 : 32), [x, y - 0.0005, z]);
  if (s.detail === 0) {
    s.add(L.lamp, cyl(s, 0.04, 0.09, "y", 16), [x, y - 0.07, z]);
    return;
  }
  if (s.fine(2)) {
    for (let i = 0; i < 4; i++) {
      const a = (i * Math.PI) / 2 + Math.PI / 4;
      screw(s, p.ssBright, [x + Math.cos(a) * 0.057, y - 0.0055, z + Math.sin(a) * 0.057], 0.0028, [Math.PI / 2, 0, 0]);
    }
    // The conduit hub on the housing's side.
    s.add(p.alu, cyl(s, 0.011, 0.02, "x", 16), [x + 0.06, y - 0.014, z]);
  }
  // The threaded collar the globe screws into, ribbed close up.
  s.add(p.alu, cyl(s, 0.045, 0.012, "y", 32), [x, y - 0.031, z]);
  if (s.fine(2)) for (let i = 0; i < 3; i++) s.add(p.alu, torus(s, 0.0452, 0.0012, "y", Math.PI * 2, 32, 4), [x, y - 0.027 - i * 0.004, z]);
  // The globe: a heavy jelly jar of frosted heat-resistant glass, and the lamp inside it.
  s.add(L.globe, hang(s, [
    [0.04, 0],
    [0.043, -0.008],
    [0.0455, -0.03],
    [0.046, -0.05],
    [0.0445, -0.075],
    [0.042, -0.092],
    [0.036, -0.103],
    [0.03, -0.108],
    [0.016, -0.1105],
    [0, -0.111],
  ], s.fine(2) ? 40 : 28), [x, y - 0.037, z]);
  s.add(L.lamp, hang(s, [
    [0.008, 0],
    [0.011, -0.012],
    [0.018, -0.026],
    [0.022, -0.035],
    [0.024, -0.05],
    [0.021, -0.062],
    [0.012, -0.069],
    [0, -0.07],
  ], s.fine(2) ? 24 : 16), [x, y - 0.04, z]);
  // The cast guard: four bent wires from the collar to a cap under the globe, and a ring at the globe's waist.
  for (let i = 0; i < 4; i++) {
    const a = (i * Math.PI) / 2 + Math.PI / 4;
    const [c, sn] = [Math.cos(a), Math.sin(a)];
    s.add(p.alu, pipe(s, [[x + c * 0.047, y - 0.03, z + sn * 0.047], [x + c * 0.052, y - 0.07, z + sn * 0.052], [x + c * 0.048, y - 0.12, z + sn * 0.048], [x + c * 0.012, y - 0.154, z + sn * 0.012]], 0.0028, 0.02, 6));
  }
  s.add(p.alu, torus(s, 0.051, 0.003, "y", Math.PI * 2, 32, 6), [x, y - 0.085, z]);
  s.add(p.alu, cyl(s, 0.016, 0.008, "y", 16), [x, y - 0.153, z]);
}

/**
 * A wet-chemical nozzle on its drop below a branch at (x, y, z): a tee on the branch, a 3/8 in drop, the hex body,
 * its orifice tip and the rubber blow-off cap that keeps grease out until the system fires. `aim` tilts it.
 */
function nozzleDrop(s: Shop, x: number, y: number, z: number, drop: number, aim: V3 = [0, 0, 0], tee = true): void {
  const p = palette();
  if (tee && s.fine(1)) s.add(p.chrome, cyl(s, 0.0145, 0.03, "x", 16), [x, y, z]);
  s.at([x, y, z], aim, () => {
    s.add(p.chrome, cyl(s, 0.0086, drop, "y", 14), [0, -drop / 2, 0]);
    if (!s.fine(1)) return;
    s.add(p.ssBright, hex(0.0115, 0.014), [0, -drop - 0.006, 0]);
    s.add(p.ssBright, cyl(s, 0.0075, 0.012, "y", 14, 0.0095), [0, -drop - 0.019, 0]);
    s.add(p.rubber, hang(s, [[0, 0], [0.0088, 0], [0.0088, -0.007], [0.006, -0.0095], [0, -0.0098]], 14), [0, -drop - 0.024, 0]);
  });
}

/**
 * A fusible-link detector inline on the 1/2 in EMT detection line at (x, y, z), the line running along x: the bracket
 * the conduit connects to at both ends, the cable through it and the soldered two-plate link that parts at its rating
 * and lets the cable go slack, which trips the release.
 */
function linkDetector(s: Shop, x: number, y: number, z: number): void {
  const p = palette();
  const L = 0.13;
  // The bracket: a U channel open below.
  blk(s, p.ssBright, [x - L / 2, y + 0.008, z - 0.013], [x + L / 2, y + 0.0105, z + 0.013], 0.0004);
  for (const sz of [-1, 1]) blk(s, p.ssBright, [x - L / 2 + 0.0013, y - 0.02, z + sz * 0.0118], [x + L / 2 - 0.0013, y + 0.008, z + sz * 0.0131], 0.0004);
  for (const sx of [-1, 1]) {
    blk(s, p.ssBright, [x + sx * (L / 2 - 0.0013), y - 0.02, z - 0.013], [x + sx * (L / 2 + 0.0002), y + 0.008, z + 0.013], 0.0004);
    // Compression connectors where the conduit joins.
    s.add(p.ssDull, cyl(s, 0.0115, 0.014, "x", 16), [x + sx * (L / 2 + 0.007), y, z]);
    if (s.fine(2)) s.add(p.ssDull, hex(0.0125, 0.006, "x"), [x + sx * (L / 2 + 0.003), y, z]);
  }
  if (!s.fine(1)) return;
  // The cable through it, and the link: two brass plates soldered at their overlap, a loop at each end.
  for (const sx of [-1, 1]) s.add(p.ssBright, cyl(s, 0.0008, L / 2 - 0.022, "x", 6), [x + sx * (L / 4 + 0.011), y - 0.004, z]);
  s.add(p.brass, rbox(s, 0.022, 0.009, 0.0007, 0.0002), [x - 0.006, y - 0.004, z + 0.0004]);
  s.add(p.copper, rbox(s, 0.022, 0.009, 0.0007, 0.0002), [x + 0.006, y - 0.004, z - 0.0004]);
  if (s.fine(2)) for (const sx of [-1, 1]) s.add(p.ssBright, torus(s, 0.003, 0.0007, "z", Math.PI * 2, 12, 4), [x + sx * 0.019, y - 0.004, z]);
}

let rodMat: THREE.MeshStandardMaterial | null = null;
/** Zinc-plated all-thread: its 13-to-the-inch thread in the bump and a little in the colour. UVs: v in threads/16. */
function allThread(): THREE.MeshStandardMaterial {
  if (rodMat) return rodMat;
  const c = document.createElement("canvas");
  c.width = 8;
  c.height = 128;
  const x = c.getContext("2d")!;
  for (let i = 0; i < 16; i++) {
    const g = x.createLinearGradient(0, i * 8, 0, i * 8 + 8);
    g.addColorStop(0, "#202020");
    g.addColorStop(0.5, "#f0f0f0");
    g.addColorStop(1, "#202020");
    x.fillStyle = g;
    x.fillRect(0, i * 8, 8, 8);
  }
  const t = new THREE.CanvasTexture(c);
  t.wrapS = t.wrapT = THREE.RepeatWrapping;
  rodMat = new THREE.MeshStandardMaterial({ color: 0xb9bec4, metalness: 0.85, roughness: 0.42, bumpMap: t, bumpScale: 1.6, roughnessMap: t });
  rodMat.userData.uv = "own";
  return rodMat;
}
/** A length of 1/2 in all-thread up y from `y0` to `y1` at (x, z), its thread drawn close up. */
function threadedRod(s: Shop, x: number, z: number, y0: number, y1: number): void {
  const len = y1 - y0;
  if (len <= 0) return;
  if (!s.fine(1)) return rod(s, palette().ssDull, 0.0064, len, "y", [x, (y0 + y1) / 2, z]);
  const g = cyl(s, 0.0064, len, "y", 12, 0.0064, true);
  if (!s.fine(2)) return void s.add(palette().ssDull, g, [x, (y0 + y1) / 2, z]);
  const uv = g.attributes.uv as THREE.BufferAttribute;
  for (let i = 0; i < uv.count; i++) uv.setY(i, (uv.getY(i) * len) / (16 * 0.00195));
  s.add(allThread(), g, [x, (y0 + y1) / 2, z]);
}

/**
 * A hanger lug on a hood top at height y, upright along x, and the hanger on it: a forged eye nut on a bolt through
 * the lug, the 1/2 in all-thread rising from it to `rodTop`, locked by a jam nut, and a coupling nut where the
 * building's rod joins it.
 */
function hangerClip(s: Shop, x: number, y: number, z: number, rodTop: number): void {
  const p = palette();
  blk(s, p.ss, [x - 0.04, y, z - 0.028], [x + 0.04, y + 0.003, z + 0.028], 0.0006);
  const sh = roundRect(0.056, 0.07, 0.012, 0, 0.035);
  if (s.fine(2)) sh.holes.push(new THREE.Path().absarc(0, 0.045, 0.009, 0, Math.PI * 2, true));
  const g = plate(s, sh, 0.003);
  s.add(p.ss, g, [x - 0.0015, y, z], [0, Math.PI / 2, 0]);
  if (s.fine(2)) for (const sz of [-1, 1]) screw(s, p.ssBright, [x + 0.025, y + 0.003, z + sz * 0.017], 0.003, [-Math.PI / 2, 0, 0]);
  const ex = x + 0.0068;
  const ye = y + 0.045;
  if (s.fine(2)) {
    s.add(p.ssDull, torus(s, 0.0125, 0.0045, "x", Math.PI * 2, 24, 8), [ex, ye, z]);
    s.add(p.ssDull, cyl(s, 0.0085, 0.024, "y", 16), [ex, ye + 0.028, z]);
    s.add(p.ssDull, cyl(s, 0.0062, 0.026, "x", 12), [x + 0.004, ye, z]);
    s.add(p.ssDull, hex(0.0105, 0.006, "x"), [x - 0.0045, ye, z]);
    s.add(p.ssDull, hex(0.0105, 0.007, "x"), [ex + 0.008, ye, z]);
    s.add(p.ssDull, hex(0.011, 0.008), [ex, ye + 0.044, z]);
    s.add(p.ssDull, hex(0.011, 0.04), [ex, rodTop - 0.02, z]);
  }
  threadedRod(s, ex, z, s.fine(2) ? ye + 0.04 : ye + 0.012, s.fine(2) ? rodTop - 0.035 : rodTop);
}

/**
 * A duct collar standing on a hood top at `y0` up to `y1`: four walls of 16-gauge, a welded flange on the top, the
 * upper edge turned out where the duct is welded on. [x0, x1] × [z0, z1] is its opening.
 */
function ductCollar(s: Shop, x0: number, x1: number, z0: number, z1: number, y0: number, y1: number): void {
  const p = palette();
  const t = 0.0016;
  blk(s, p.ss, [x0, y0, z0], [x1, y1, z0 + t], 0.0005);
  blk(s, p.ss, [x0, y0, z1 - t], [x1, y1, z1], 0.0005);
  blk(s, p.ss, [x0, y0, z0 + t], [x0 + t, y1, z1 - t], 0.0005);
  blk(s, p.ss, [x1 - t, y0, z0 + t], [x1, y1, z1 - t], 0.0005);
  // The flange round its foot, and its top edge flared out.
  const f = 0.025;
  sheetWithHoles(s, p.ss, x0 - f, x1 + f, z0 - f, z1 + f, y0, y0 + 0.002, [[x0, x1, z0, z1]]);
  if (s.fine(1)) {
    for (const [a, b] of [
      [z0 - 0.004, z0],
      [z1, z1 + 0.004],
    ] as const) blk(s, p.ss, [x0 - 0.004, y1 - 0.006, a], [x1 + 0.004, y1, b], 0.0005);
    for (const [a, b] of [
      [x0 - 0.004, x0],
      [x1, x1 + 0.004],
    ] as const) blk(s, p.ss, [a, y1 - 0.006, z0], [b, y1, z1], 0.0005);
  }
  if (s.fine(2)) {
    // Weld bead round the foot, and the flange's tack screws.
    for (const [a, b] of [
      [z0 - 0.0025, z0],
      [z1, z1 + 0.0025],
    ] as const) blk(s, p.ssDull, [x0 - 0.002, y0 + 0.002, a], [x1 + 0.002, y0 + 0.0045, b], 0.001);
    for (const [a, b] of [
      [x0 - 0.0025, x0],
      [x1, x1 + 0.0025],
    ] as const) blk(s, p.ssDull, [a, y0 + 0.002, z0], [b, y0 + 0.0045, z1], 0.001);
    const n = Math.max(2, Math.round((x1 - x0) / 0.15));
    for (let i = 0; i <= n; i++) for (const z of [z0 - f / 2, z1 + f / 2]) screw(s, p.ssBright, [x0 + ((x1 - x0) * i) / n, y0 + 0.002, z], 0.003, [-Math.PI / 2, 0, 0]);
  }
}

/**
 * A hood skin's lower edge, as (u, v) for `pathSweep` with u inward from the outer face: the skin creased out into a
 * 3/4 in hem, turned under round a tight bend, the return across the double wall `wl` thick, its inner corner
 * radiused up into the liner, closed at `top` where the panels stand on it.
 */
function hemSection(s: Shop, wl: number, top: number): P2[] {
  if (s.detail === 0) return [[0, top], [0, 0], [wl, 0], [wl, top]];
  if (s.detail === 1) return [[0, top], [0, 0.0055], [0.0016, 0.0016], [0.0055, 0], [wl, 0], [wl, top]];
  const n = s.seg(6, 2);
  return [[0, top], [0, 0.0205], [-0.0013, 0.019], [-0.0013, 0.0055], ...arc(0.0042, 0.0055, 0.0055, Math.PI, Math.PI * 1.5, n).slice(1), [wl - 0.003, 0], ...arc(wl - 0.003, 0.003, 0.003, -Math.PI / 2, 0, Math.max(2, n - 2)).slice(1), [wl, top]];
}

/**
 * A formed outside corner of a hood, radius `r` about (x, z), from y0 to y1: the quarter round that shows, facing out
 * of a w × d plan.
 */
function cornerPost(s: Shop, r: number, x: number, z: number, w: number, d: number, y0: number, y1: number): void {
  if (r <= 0) return;
  const t0 = Math.atan2(x < w / 2 ? -1 : 1, z > -d / 2 ? 1 : -1) - Math.PI / 4;
  s.add(palette().ss, new THREE.CylinderGeometry(r, r, y1 - y0, s.seg(8, 4), 1, true, t0, Math.PI / 2), [x, (y0 + y1) / 2, z]);
}

/** A standing seam on a hood top, as an (x, y) centreline: both sheets turned up 22 mm and folded over together. */
function standingSeam(s: Shop): P2[] {
  return [[-0.013, 0.0006], [-0.0024, 0.0006], [-0.0024, 0.019], ...arc(0, 0.019, 0.0024, Math.PI, 0, s.seg(6, 3)).slice(1, -1), [0.0024, 0.019], [0.0024, 0.0006], [0.013, 0.0006]];
}

/** Positions along a hood top for its standing seams, about every `pitch`, moved off anything in their way. */
function seamPositions(w: number, pitch: number, blocked: readonly [number, number][]): number[] {
  const n = Math.round(w / pitch);
  const out: number[] = [];
  for (let i = 1; i < n; i++) {
    let x = (w * i) / n;
    for (let k = 0; k < 40; k++) {
      const hit = blocked.find(([a, b]) => x > a && x < b);
      if (!hit) break;
      const want = (w * i) / n;
      x = Math.abs(hit[0] - want) < Math.abs(hit[1] - want) ? hit[0] - 0.005 : hit[1] + 0.005;
    }
    if (x > 0.15 && x < w - 0.15 && !blocked.some(([a, b]) => x > a && x < b) && !out.some((q) => Math.abs(q - x) < 0.3)) out.push(x);
  }
  return out;
}

/**
 * A Type I wall canopy, or the double island when `island`: an 18-gauge box canopy with formed corners and a hemmed
 * drip edge, double-walled at the ends, standing seams across its top, hung on all-thread from welded lugs. A front
 * make-up air plenum throws supply air forward through perforated panels; under it the capture reservoir, lit by
 * vapour-proof globes, with a perimeter gutter at its foot. A 55° bank of baffle filters stands over a grease trough
 * at the back that drains to a grease cup (on the island, a V of two banks back to back over one centre trough), the
 * exhaust plenum behind them and its duct collars on top. The wet-chemical system's piping hangs in the reservoir,
 * nozzles aimed down at the line, and its detection line runs fusible links in front of the filters.
 *
 * y = 0 is the hood's lower edge; the canopy is h - 0.07 tall and the collars stand on it to h.
 */
function typeOneHood(island: boolean): Builder {
  return (s, e) => {
    const p = palette();
    const { w, d, h } = e;
    const hb = h - 0.07;
    const T = 0.0012;
    const wall = 0.025;
    const D = island ? d / 2 : d;
    const mua = 0.2;
    const perf = perforatedSheet(0.0065, 0.55);
    const r = s.fine(1) ? 0.0012 : 0;
    const L = hoodLamps();
    const yH = 0.035;
    const R = s.detail ? wall : 0;
    const ceil = hb - T;

    // The hemmed lower edge round the outside, the formed corners, the skins between them.
    const corners: P2[] = island ? [[0, 0], [w, 0], [w, -d], [0, -d]] : [[0, -d], [0, 0], [w, 0], [w, -d]];
    s.add(p.ss, pathSweep(hemSection(s, wall, yH), filletPath(corners, R, s.seg(8, 2), island), 0, island));
    const posts: P2[] = island ? [[R, -R], [w - R, -R], [w - R, -d + R], [R, -d + R]] : [[R, -R], [w - R, -R]];
    for (const [x, z] of posts) cornerPost(s, R, x, z, w, d, yH, ceil);
    for (const [x0, sg] of [
      [0, 1],
      [w, -1],
    ] as const) {
      blk(s, p.ss, [x0, yH, island ? -d + R : -d + T], [x0 + sg * T, ceil, -R], r);
      blk(s, p.ss, [x0 + sg * (wall - T), yH, island ? -d + mua : -d + T], [x0 + sg * wall, ceil, -mua], r);
      if (!island) blk(s, p.ss, [x0, yH, -d], [x0 + sg * wall, ceil, -d + T], 0);
    }
    // The back panel against the wall.
    if (!island) blk(s, p.ss, [wall, 0, -d], [w - wall, ceil, -d + T], r);

    // Collars: exhaust over the plenum, make-up air over the front plenum(s).
    const nEx = Math.max(1, Math.round(w / 3.6));
    const exW = island ? 0.6 : 0.5;
    const exZ: [number, number] = island ? [-d / 2 - 0.2, -d / 2 + 0.2] : [-d + 0.05, -d + 0.3];
    const exhausts: [number, number, number, number][] = [];
    for (let i = 0; i < nEx; i++) {
      const xc = (w * (i + 0.5)) / nEx;
      exhausts.push([xc - exW / 2, xc + exW / 2, exZ[0], exZ[1]]);
    }
    const nMua = Math.max(1, Math.round(w / 2.4));
    const muas: [number, number, number, number][] = [];
    for (let i = 0; i < nMua; i++) {
      const xc = (w * (i + 0.5)) / nMua;
      muas.push([xc - 0.25, xc + 0.25, -mua + 0.04, -0.04]);
      if (island) muas.push([xc - 0.25, xc + 0.25, -d + 0.04, -d + mua - 0.04]);
    }
    // The top: one sheet, its corners formed where the skins' are, cut for the collars.
    const Rb = island ? R : 0;
    const top = planShape(filletPath([[0, 0], [w, 0], [w, -d], [0, -d]], [R, R, Rb, Rb], s.seg(8, 2), true));
    for (const [x0, x1, z0, z1] of [...exhausts, ...muas]) top.holes.push(rectPath(new THREE.Path(), x1 - x0, z1 - z0, 0, (x0 + x1) / 2, -(z0 + z1) / 2, true) as THREE.Path);
    const tg = plate(s, top, T, 0, 8);
    tg.rotateX(-Math.PI / 2);
    s.add(p.ss, tg, [0, ceil, 0]);
    for (const [x0, x1, z0, z1] of exhausts) ductCollar(s, x0, x1, z0, z1, hb, h);
    for (const [x0, x1, z0, z1] of muas) {
      ductCollar(s, x0, x1, z0, z1, hb, h - 0.03);
      // The supply collar's fire damper: a blade held open by its own fusible link.
      if (s.fine(2)) {
        s.add(p.ssDull, rbox(s, x1 - x0 - 0.006, 0.0015, (z1 - z0) * 0.9, 0.0004), [(x0 + x1) / 2, hb + 0.03, (z0 + z1) / 2], [0.5, 0, 0]);
        s.add(p.ssDull, cyl(s, 0.003, x1 - x0, "x", 8), [(x0 + x1) / 2, hb + 0.03, (z0 + z1) / 2]);
        s.add(p.brass, rbox(s, 0.02, 0.008, 0.0008, 0.0002), [(x0 + x1) / 2, hb + 0.045, z0 + 0.012]);
      }
    }

    // Hanger lugs at the corners and every 2.4 m between, each on its rod.
    const nh = Math.max(1, Math.round(w / 2.4));
    const lugX: number[] = [];
    for (let i = 0; i <= nh; i++) lugX.push(0.1 + ((w - 0.2) * i) / nh);
    for (const x of lugX) for (const z of [-0.09, -d + 0.09]) hangerClip(s, x, hb, z, h + 0.14);

    // A field joint on hoods too long to ship in one piece: bolted flanges inside, a cover strip on the face, and a
    // standing seam over the top.
    const xj = w > 4.9 ? w / 2 : -1;
    if (xj > 0) {
      blk(s, p.ss, [xj - 0.02, yH, 0], [xj + 0.02, ceil, 0.0022], 0.0008);
      for (const sx of [-1, 1]) blk(s, p.ss, [xj + sx * 0.0015, ceil - 0.03, -d + 0.002], [xj + sx * 0.0035, ceil, -0.002], 0);
      if (s.fine(2))
        for (let z = -0.25; z > -d + 0.1; z -= 0.15) {
          s.add(p.ssBright, hex(0.006, 0.012, "x"), [xj, ceil - 0.015, z]);
        }
    }

    // One side of the hood in its own frame: front at z = 0, its "wall" (the back wall, or the island's centre) at -D.
    const filter = baffleFilter(s);
    const phi = Math.PI / 2 - SLOPE;
    const cphi = Math.cos(phi);
    const yb = 0.1;
    const zb = island ? -D + 0.062 : -D + 0.08;
    const zpf = zb + 0.302;
    const zl = -mua - 0.17;
    const zf = zpf + 0.1;
    const yf = hb - 0.06;
    const zc = zpf + 0.035;
    const yc = hb - 0.035;
    const nl = Math.max(1, Math.round((w - 2 * wall) / 1.2));
    const lx: number[] = [];
    for (let i = 0; i < nl; i++) lx.push(wall + ((w - 2 * wall) * (i + 0.5)) / nl);
    const usable = w - 2 * wall - 0.01;
    const nfl = Math.max(1, Math.floor(usable / FW));
    const fx0 = (w - nfl * FW) / 2;
    const cupX = w - wall - 0.08;
    const side = () => {
      // The make-up air plenum: inner wall, floor, the front face over its perforated supply panels and diffuser.
      blk(s, p.ss, [wall, 0, -mua], [w - wall, ceil, -mua + T], r);
      blk(s, p.ss, [wall, 0, -mua + T], [w - wall, T, -wall], 0);
      blk(s, p.ss, [R, 0.305, -T], [w - R, ceil, 0], r);
      const nSec = Math.max(1, Math.round((w - 2 * wall) / 1.22));
      const secW = (w - 2 * wall - 0.03) / nSec;
      const mull = (a: number, b: number) => blk(s, p.ss, [a, yH, -T], [b, 0.305, 0], r);
      mull(R, wall + 0.03);
      mull(w - wall - 0.03, w - R);
      for (let i = 1; i < nSec; i++) mull(wall + secW * i, wall + 0.03 + secW * i);
      for (let i = 0; i < nSec; i++) {
        const xa = wall + 0.03 + secW * i;
        const xb = xa + secW - 0.03;
        if (s.fine(1)) blk(s, perf, [xa - 0.006, yH - 0.004, -0.0058], [xb + 0.006, 0.31, -0.0049], 0);
        else blk(s, p.ssDull, [xa, yH, -0.003], [xb, 0.305, -0.002], 0);
      }
      blk(s, p.charcoal, [wall, 0.028, -0.075], [w - wall, 0.312, -0.074], 0);
      if (s.fine(2)) for (let i = 0; i <= nSec; i++) for (const y of [yH + 0.012, 0.32]) screw(s, p.ssBright, [wall + 0.015 + secW * i, y, 0], 0.0028);

      // The capture reservoir's lights, wired through boxes on the top, and the pools they throw round them.
      for (const x of lx) {
        hoodLight(s, L, x, ceil, zl);
        lightPool(s, L, [x, ceil - 0.003, zl], 1.15, 2 * Math.min(-mua - zl, zl - zpf) - 0.02, [Math.PI / 2, 0, 0]);
        lightPool(s, L, [x, ceil / 2 + 0.03, -mua - 0.003], 1.05, ceil - 0.1, [0, Math.PI, 0]);
        if (s.fine(1)) {
          blk(s, p.ssDull, [x - 0.05, hb, zl - 0.05], [x + 0.05, hb + 0.045, zl + 0.05], 0.002);
          if (s.fine(2)) for (const sx of [-1, 1]) for (const sz of [-1, 1]) screw(s, p.ssBright, [x + sx * 0.04, hb + 0.045, zl + sz * 0.04], 0.0022, [-Math.PI / 2, 0, 0]);
        }
      }
      if (s.fine(1)) {
        s.add(p.ssDull, cyl(s, 0.009, lx[lx.length - 1]! - 0.12, "x", 12), [(0.12 + lx[lx.length - 1]!) / 2, hb + 0.034, zl]);
        blk(s, p.ssDull, [0.07, hb, zl - 0.055], [0.17, hb + 0.055, zl + 0.055], 0.002);
      }

      // The filter bank: filters stamped across, blank fillers at the ends, the top channel and the plenum's front.
      s.at([0, yb, zb], [phi, 0, 0], () => {
        for (let i = 0; i < nfl; i++) stamp(s, filter, [fx0 + FW * (i + 0.5), 0, 0]);
        for (const [a, b] of [
          [wall, fx0],
          [fx0 + nfl * FW, w - wall],
        ] as const) {
          if (b - a < 0.005) continue;
          blk(s, p.ss, [a, 0.004, -FT], [b, FH - 0.004, -0.0028], 0);
        }
        // The top channel the filters lift into.
        sweepX(s, p.ss, [[0.012, FH - 0.024], [0.012, FH + 0.006], [-FT - 0.004, FH + 0.006], [-FT - 0.004, FH - 0.012]], 0.0012, wall, w - wall);
        for (const x of lx) lightPool(s, L, [x, FH * 0.5, 0.012], 1.3, FH * 0.96, [0, 0, 0]);
      });
      blk(s, p.ss, [wall, yb + FH * cphi - 0.012, zpf], [w - wall, ceil, zpf + T], r);
      // The ledge inside the trough that the filters stand on.
      blk(s, p.ss, [wall, yb - 0.006, zb - 0.045], [w - wall, yb - 0.003, zb + 0.004], 0);

      // Fire suppression: the 1/2 in branch line with its nozzles aimed down at the line, hung from the ceiling.
      const xa = 0.16;
      const xb = w - 0.16;
      s.add(p.chrome, pipe(s, [[xb, yf, zf], [xa, yf, zf], [xa, h - 0.006, zf]], 0.0107, 0.035, 12));
      if (s.fine(1)) {
        s.add(p.ssBright, hex(0.017, 0.024), [xa, h - 0.034, zf]);
        s.add(p.ssBright, hex(0.02, 0.012), [xa, hb + 0.006, zf]);
        if (s.fine(2)) s.add(p.ssBright, hex(0.02, 0.01), [xa, ceil - 0.005, zf]);
        s.add(p.chrome, cyl(s, 0.0118, 0.012, "x", 12), [xb + 0.006, yf, zf]);
      }
      const nn = Math.max(2, Math.round((w - 0.4) / 0.6));
      for (let i = 0; i < nn; i++) {
        const x = 0.3 + ((w - 0.6) * (i + 0.5)) / nn;
        nozzleDrop(s, x, yf, zf, 0.065, [i % 2 ? 0.12 : -0.12, 0, 0]);
      }
      if (s.fine(2)) {
        const nh2 = Math.max(2, Math.round(w / 1.1));
        for (let i = 0; i < nh2; i++) {
          const x = 0.25 + ((w - 0.5) * (i + 0.5)) / nh2;
          s.add(p.ssBright, cyl(s, 0.003, hb - yf - 0.012, "y", 8), [x, (hb + yf + 0.012) / 2, zf]);
          s.add(p.ssBright, torus(s, 0.0125, 0.0022, "x", Math.PI * 2, 16, 5), [x, yf, zf]);
        }
      }

      // The detection line: EMT conduit along the bank's top, a link detector every metre or so, a pulley elbow at the
      // end where it turns up through the top to the release.
      const nd = Math.max(2, Math.round((w - 0.4) / 1.3));
      const dx: number[] = [];
      for (let i = 0; i < nd; i++) dx.push(0.35 + ((w - 0.7) * (i + 0.5)) / nd);
      let from = 0.12;
      for (const x of dx) {
        s.add(p.ssDull, cyl(s, 0.0089, x - 0.08 - from, "x", 12), [(from + x - 0.08) / 2, yc, zc]);
        linkDetector(s, x, yc, zc);
        from = x + 0.08;
      }
      s.add(p.alu, cyl(s, 0.028, 0.022, "z", 24), [0.1, yc + 0.012, zc]);
      if (s.fine(1)) s.add(p.alu, cyl(s, 0.022, 0.004, "z", 24), [0.1, yc + 0.012, zc + 0.013]);
      s.add(p.ssDull, cyl(s, 0.0089, h - 0.004 - yc - 0.03, "y", 12), [0.1, (h - 0.004 + yc + 0.03) / 2, zc]);
      if (s.fine(1)) s.add(p.ssDull, hex(0.016, 0.014), [0.1, hb + 0.007, zc]);
    };
    side();
    if (island) s.at([w, 0, -d], [0, Math.PI, 0], side);

    // The reservoir's perimeter gutter round the foot of its walls: a turned-up lip, its top rolled.
    const gut = sheetOutline([[0, 0.0006], [0.018, 0.0006], [0.0194, 0.002], [0.0194, 0.014], ...arc(0.0172, 0.014, 0.0022, 0, Math.PI, s.seg(5, 2)).slice(1)], 0.0012);
    const gPath: P2[] = island ? [[wall, -mua], [w - wall, -mua], [w - wall, -d + mua], [wall, -d + mua]] : [[wall, -d + T], [wall, -mua], [w - wall, -mua], [w - wall, -d + T]];
    if (s.detail) s.add(p.ss, pathSweep(gut, gPath, 0, island));

    // The grease trough: on a wall canopy a U against the wall, its front lip rolled in; on the island one trough on
    // the centreline under both banks. It falls to the drain at its right end over the grease cup.
    const nArc = s.seg(6, 2);
    let zcup: number;
    if (island) {
      const zm = -d / 2;
      const b = 0.082;
      sweepX(s, p.ss, [[zm - b + 0.007, 0.116], [zm - b, 0.11], [zm - b, 0.075], ...arc(zm - b + 0.015, 0.075, 0.015, Math.PI, Math.PI * 1.5, nArc).slice(1), [zm + b - 0.015, 0.06], ...arc(zm + b - 0.015, 0.075, 0.015, -Math.PI / 2, 0, nArc).slice(1), [zm + b, 0.11], [zm + b - 0.007, 0.116]], 0.0012, wall, w - wall);
      zcup = zm;
    } else {
      const zt = -D + 0.004;
      sweepX(s, p.ss, [[zt, 0.135], [zt, 0.075], ...arc(zt + 0.015, 0.075, 0.015, Math.PI, Math.PI * 1.5, nArc).slice(1), [zt + 0.08, 0.06], ...arc(zt + 0.08, 0.075, 0.015, -Math.PI / 2, 0, nArc).slice(1), [zt + 0.095, 0.11], [zt + 0.088, 0.116]], 0.0012, wall, w - wall);
      zcup = zt + 0.0475;
    }
    // End caps closing the trough, and the drain spud into the cup.
    if (s.fine(1)) for (const x of [wall + 0.0008, w - wall - 0.0008]) blk(s, p.ss, [x - 0.0006, 0.062, zcup - (island ? 0.08 : 0.045)], [x + 0.0006, 0.108, zcup + (island ? 0.08 : 0.045)], 0);
    s.add(p.ss, cyl(s, 0.011, 0.02, "y", 16), [cupX, 0.052, zcup]);
    // The grease cup: a drawn cup with a rolled rim on its bail, hooked under the drain, part full of old grease.
    s.add(p.ss, lathe(s, [
      [0, 0.004],
      [0.032, 0.004],
      [0.036, 0.0065],
      [0.0385, 0.012],
      [0.0395, 0.05],
      [0.0412, 0.0515],
      [0.0412, 0.0535],
      [0.0392, 0.0535],
      [0.0375, 0.05],
      [0.0365, 0.0125],
      [0.0345, 0.0095],
      [0, 0.0095],
    ], 32), [cupX, 0, zcup]);
    if (s.fine(1)) {
      s.add(p.ssBright, torus(s, 0.0425, 0.002, "y", Math.PI * 2, 32, 5), [cupX, 0.04, zcup]);
      for (const sz of [-1, 1]) s.add(p.ssBright, pipe(s, [[cupX, 0.04, zcup + sz * 0.0425], [cupX, 0.062, zcup + sz * 0.0425]], 0.0022, 0.004, 6));
      if (s.fine(2)) {
        s.add(p.ssBright, pipe(s, [[cupX - 0.04, 0.03, zcup], [cupX - 0.042, 0.06, zcup], [cupX, 0.075, zcup], [cupX + 0.042, 0.06, zcup], [cupX + 0.04, 0.03, zcup]], 0.0014, 0.01, 6));
        s.add(greaseMat(), cyl(s, 0.0368, 0.002, "y", 28), [cupX, 0.034, zcup]);
      }
    }

    // The duct nozzle in each exhaust collar, fed down through the top: the pipe drops into the plenum, turns under the
    // collar's opening and points up the duct.
    for (const [x0, x1, z0, z1] of exhausts) {
      const xc = (x0 + x1) / 2;
      const zm = (z0 + z1) / 2;
      const xp = x1 + 0.06;
      s.add(p.chrome, pipe(s, [[xp, h - 0.006, zm], [xp, hb - 0.05, zm], [xc, hb - 0.05, zm], [xc, hb - 0.005, zm]], 0.0086, 0.03, 12));
      if (s.fine(1)) {
        s.add(p.ssBright, hex(0.017, 0.012), [xp, hb + 0.006, zm]);
        s.add(p.ssBright, hex(0.0115, 0.014), [xc, hb + 0.002, zm]);
        s.add(p.rubber, cyl(s, 0.0088, 0.01, "y", 14), [xc, hb + 0.014, zm]);
      }
    }

    // Standing seams across the top where its sheets join, clear of the collars, boxes and lugs.
    if (s.fine(1)) {
      const blocked: [number, number][] = [
        ...muas.map(([x0, x1]) => [x0 - 0.045, x1 + 0.045] as [number, number]),
        ...exhausts.map(([x0, x1]) => [x0 - 0.045, x1 + 0.095] as [number, number]),
        ...lx.map((x) => [x - 0.065, x + 0.065] as [number, number]),
        ...lugX.map((x) => [x - 0.05, x + 0.05] as [number, number]),
        [0.05, 0.19],
        [w - 0.19, w - 0.05],
      ];
      if (xj > 0) blocked.push([xj - 0.3, xj + 0.3]);
      const seams = seamPositions(w, 1.22, blocked);
      if (xj > 0) seams.push(xj);
      for (const x of seams) sweepZ(s, p.ss, standingSeam(s), 0.0012, -d + 0.004, -0.004, x, hb);
    }

    // A plate on the left end's inner skin: the listing label, blank.
    dataPlate(s, [wall + 0.0002, hb * 0.55, -mua - 0.12], [0, Math.PI / 2, 0], 0.1, 0.06);

    return lampHooks;
  };
}

let greaseM: THREE.MeshStandardMaterial | null = null;
/** Old fryer grease in a cup: dark amber, glossy, opaque at any depth that matters. */
function greaseMat(): THREE.MeshStandardMaterial {
  if (!greaseM) {
    greaseM = new THREE.MeshStandardMaterial({ color: 0x4a2c0c, roughness: 0.18, metalness: 0 });
    greaseM.userData.uv = "own";
  }
  return greaseM;
}

/**
 * A Type II condensate hood, over a dishmachine or steamer: no grease filters, a smooth welded liner, and a perimeter
 * gutter inside the skirts that catches what condenses on the liner and drains it from one corner. An internal baffle
 * under the collar makes the vapour turn and shed its water; its low edge drips into the back gutter. One light, and
 * the exhaust collar on top. The skirts have formed corners and a hemmed drip edge. y = 0 is the skirts' lower edge.
 */
const condensateHood: Builder = (s, e) => {
  const p = palette();
  const { w, d, h } = e;
  const hb = h - 0.065;
  const T = 0.0012;
  const wall = 0.035;
  const gy = 0.02;
  const yH = 0.035;
  const ceil = hb - T;
  const r = s.fine(1) ? 0.0012 : 0;
  const R = s.detail ? wall : 0;
  const L = hoodLamps();
  // The hemmed edge all round, the formed corners and the skins between them; the liner inside.
  s.add(p.ss, pathSweep(hemSection(s, wall, yH), filletPath([[0, 0], [w, 0], [w, -d], [0, -d]], R, s.seg(8, 2), true), 0, true));
  for (const [x, z] of [[R, -R], [w - R, -R], [w - R, -d + R], [R, -d + R]] as const) cornerPost(s, R, x, z, w, d, yH, ceil);
  blk(s, p.ss, [R, yH, -T], [w - R, ceil, 0], r);
  blk(s, p.ss, [R, yH, -d], [w - R, ceil, -d + T], r);
  for (const x of [0, w - T]) blk(s, p.ss, [x, yH, -d + R], [x + T, ceil, -R], r);
  blk(s, p.ss, [wall, yH, -wall], [w - wall, ceil, -wall + T], r);
  blk(s, p.ss, [wall, yH, -d + wall - T], [w - wall, ceil, -d + wall], r);
  for (const x of [wall - T, w - wall]) blk(s, p.ss, [x, yH, -d + wall], [x + T, ceil, -wall], r);
  // The top, its corners formed, cut for the collar.
  const cx0 = w / 2 - 0.225;
  const cx1 = w / 2 + 0.225;
  const cz0 = -d + 0.06;
  const cz1 = -d + 0.34;
  const top = planShape(filletPath([[0, 0], [w, 0], [w, -d], [0, -d]], R, s.seg(8, 2), true));
  top.holes.push(rectPath(new THREE.Path(), cx1 - cx0, cz1 - cz0, 0, (cx0 + cx1) / 2, -(cz0 + cz1) / 2, true) as THREE.Path);
  const tg = plate(s, top, T, 0, 8);
  tg.rotateX(-Math.PI / 2);
  s.add(p.ss, tg, [0, ceil, 0]);
  ductCollar(s, cx0, cx1, cz0, cz1, hb, h);

  // The perimeter gutter round the liner's foot, one welded piece: a U with a rolled lip, its floor 20 mm up.
  const n = s.fine(2) ? 6 : 2;
  const gut: P2[] = [[0.0006, 0.065], [0.0006, gy + 0.008], ...arc(0.0086, gy + 0.008, 0.008, Math.PI, Math.PI * 1.5, n).slice(1), [0.037, gy], ...arc(0.037, gy + 0.008, 0.008, -Math.PI / 2, 0, n).slice(1), [0.045, gy + 0.035], ...arc(0.041, gy + 0.035, 0.004, 0, Math.PI, 3).slice(1)];
  s.add(p.ss, pathSweep(sheetOutline(gut, 0.0012), [[wall, -wall], [w - wall, -wall], [w - wall, -d + wall], [wall, -d + wall]], 0, true));
  // The drain: a 1 in half coupling under the low corner, a close nipple and its cap.
  const dx = w - wall - 0.022;
  const dz = -wall - 0.022;
  s.add(p.ss, cyl(s, 0.017, 0.012, "y", 20), [dx, gy - 0.0066, dz]);
  s.add(p.ssDull, cyl(s, 0.0167, 0.008, "y", 20), [dx, 0.006, dz]);
  if (s.fine(1)) s.add(p.ssDull, hex(0.02, 0.006), [dx, 0.0035, dz]);

  // The baffle: a sloped plate under the collar, its upper edge hemmed, its low edge turned down into a drip lip over
  // the back gutter; hung from the top by straps at each end.
  const [za, ya] = [-d + 0.55, hb - 0.095];
  const [zb, yb] = [-d + 0.058, 0.16];
  const Lb = Math.hypot(za - zb, ya - yb);
  const ang = Math.atan2(ya - yb, za - zb);
  s.at([w / 2, yb, zb], [-ang, 0, 0], () => {
    blk(s, p.ss, [-w / 2 + wall + 0.05, -0.0008, 0], [w / 2 - wall - 0.05, 0.0008, Lb], r);
    s.add(p.ss, cyl(s, 0.004, w - 2 * wall - 0.1, "x", 12), [0, 0.003, Lb]);
    blk(s, p.ss, [-w / 2 + wall + 0.05, -0.018, -0.0008], [w / 2 - wall - 0.05, -0.0008, 0.0008], 0.0004);
  });
  for (const x of [wall + 0.08, w - wall - 0.08]) {
    s.add(p.ss, pipe(s, [[x, ceil, za - 0.03], [x, ya + 0.004, za - 0.03]], 0.0025, 0.01, 8));
    s.add(p.ss, pipe(s, [[x, ceil, zb + 0.12], [x, yb + 0.12 * Math.tan(ang) + 0.004, zb + 0.12]], 0.0025, 0.01, 8));
  }

  // The light and the pools it throws, its box on top, hanger lugs and rods, the label.
  hoodLight(s, L, w / 2, ceil, -0.3);
  lightPool(s, L, [w / 2, ceil - 0.003, -0.3], Math.min(1.1, w - 2 * wall - 0.1), 0.5, [Math.PI / 2, 0, 0]);
  lightPool(s, L, [w / 2, (ceil + 0.07) / 2, -wall - 0.003], Math.min(1.0, w - 2 * wall - 0.1), ceil - 0.1, [0, Math.PI, 0]);
  if (s.fine(1)) blk(s, p.ssDull, [w / 2 - 0.05, hb, -0.35], [w / 2 + 0.05, hb + 0.045, -0.25], 0.002);
  for (const x of [0.1, w - 0.1]) for (const z of [-0.09, -d + 0.09]) hangerClip(s, x, hb, z, h + 0.14);
  dataPlate(s, [wall + 0.0002, hb * 0.55, -0.25], [0, Math.PI / 2, 0], 0.09, 0.055);
  return lampHooks;
};

register({
  "NSF-HOOD-TYPE2-CONDENSATE": condensateHood,
  "NSF-HOOD-TYPE2-CONDENSATE-60IN": condensateHood,
  "NSF-HOOD-TYPE1-CANOPY": typeOneHood(false),
  "NSF-HOOD-TYPE1-CANOPY-96IN": typeOneHood(false),
  "NSF-HOOD-TYPE1-CANOPY-144IN": typeOneHood(false),
  "NSF-HOOD-TYPE1-CANOPY-168IN": typeOneHood(false),
  "NSF-HOOD-TYPE1-CANOPY-288IN": typeOneHood(false),
  "NSF-HOOD-TYPE1-ISLAND-144IN": typeOneHood(true),
  "NSF-HOOD-TYPE1-ISLAND-120IN": typeOneHood(true),
});

/* ------------------------------------------------------------------------------------------------------------ */
/* Sinks                                                                                                          */
/* ------------------------------------------------------------------------------------------------------------ */

/**
 * The raised rolled rim round a sink top's front and ends, as a (u, v) centreline for `pathSweep`, u in from its outer
 * face and v up from the deck: an upstand rising out of the deck (and out of the drainboards, which pitch up toward it
 * under a rim that stays level), rolled over into a bead, down the apron and hemmed back under.
 */
function sinkRim(s: Shop): P2[] {
  const n = s.fine(2) ? s.seg(10, 3) : 3;
  return [[0.018, -0.012], [0.018, 0.02], ...arc(0.009, 0.02, 0.009, 0, Math.PI, n).slice(1), [0, -0.04], ...arc(0.006, -0.04, 0.006, Math.PI, Math.PI * 1.5, s.fine(2) ? 3 : 1).slice(1), [0.018, -0.046]];
}

/**
 * A splash-mounted faucet body on 8 in centres at (x, y), its shanks through the splash face at z: two valve bonnets
 * with wrist blades (hot on the left, its index red), and between them the outlet, a swing spout of `spout` reach on a
 * swivel joint, or nothing when `spout` is 0 (the pre-rinse riser takes the outlet).
 */
function splashFaucet(s: Shop, x: number, y: number, z: number, spout: number, name: string): void {
  const p = palette();
  const c = 0.1016;
  const zb = z + 0.034;
  if (s.detail === 0) {
    s.add(p.chrome, new THREE.BoxGeometry(2 * c + 0.04, 0.03, 0.034), [x, y, z + 0.02]);
    if (!spout) return;
    const j = s.joint(name, [x, y + 0.03, zb]);
    j.add(p.chrome, new THREE.BoxGeometry(0.02, 0.02, spout), [0, 0.03, spout / 2]);
    return;
  }
  for (const sx of [-1, 1]) {
    s.add(p.chrome, cyl(s, 0.011, 0.036, "z", 16, 0.011, true), [x + sx * c, y, z + 0.018]);
    if (s.fine(2)) s.add(p.chrome, cyl(s, 0.019, 0.004, "z", 20), [x + sx * c, y, z + 0.002]);
    // The bonnet, then a wrist blade on its stem, swept out and forward so a forearm can work it.
    const bonnet: [number, number][] = s.fine(2) ? [[0.016, -0.012], [0.017, 0.012], [0.013, 0.02], [0.0085, 0.024], [0.0085, 0.04], [0, 0.04]] : [[0.017, -0.012], [0.016, 0.018], [0.0085, 0.024], [0, 0.04]];
    s.add(p.chrome, lathe(s, bonnet, s.fine(2) ? 20 : 16), [x + sx * c, y, zb]);
    s.at([x + sx * c, y + 0.042, zb], [0, sx * 0.55, 0], () => {
      s.add(p.chrome, cyl(s, 0.011, 0.012, "y", 16), [0, 0, 0]);
      s.add(p.chrome, taperBlade(s, 0.105, 0.016, 0.022, 0.004), [0, 0.004, 0], [0, -Math.PI / 2, 0]);
      if (s.fine(2)) s.add(sx < 0 ? p.ledRed : p.enamel, cyl(s, 0.0065, 0.0012, "y", 14), [0, 0.0066, 0]);
    });
  }
  s.add(p.chrome, cyl(s, 0.0135, 2 * c, "x", 20, 0.0135, true), [x, y, zb]);
  if (!spout) return;
  s.add(p.chrome, cyl(s, 0.015, 0.03, "y", 20), [x, y + 0.016, zb]);
  const j = s.joint(name, [x, y + 0.03, zb]);
  j.add(p.chrome, cyl(j, 0.0125, 0.03, "y", 18), [0, 0.012, 0]);
  j.add(p.chrome, pipe(j, [[0, 0.024, 0], [0, 0.03, 0.03], [0, 0.03, spout - 0.03], [0, -0.005, spout]], 0.0095, 0.025, 14));
  j.add(p.chrome, cyl(j, 0.012, 0.02, "y", 16), [0, -0.012, spout]);
  if (j.fine(2)) j.add(p.enamel, cyl(j, 0.008, 0.001, "y", 14), [0, -0.0225, spout]);
}

/** A wrist blade: a flat paddle tapering from `w0` at its hub to `w1`, `len` long along +x, `t` thick. */
function taperBlade(s: Shop, len: number, w0: number, w1: number, t: number): THREE.BufferGeometry {
  const sh = new THREE.Shape();
  sh.moveTo(0, -w0 / 2);
  sh.lineTo(len - w1 / 2, -w1 / 2);
  sh.quadraticCurveTo(len, -w1 / 2, len, 0);
  sh.quadraticCurveTo(len, w1 / 2, len - w1 / 2, w1 / 2);
  sh.lineTo(0, w0 / 2);
  sh.lineTo(0, -w0 / 2);
  const g = plate(s, sh, t, s.fine(2) ? 0.0012 : 0, 8);
  g.rotateX(-Math.PI / 2);
  return g;
}

/**
 * A pre-rinse unit on the splash at (x, y, z): its mixing body, the riser up to a wall bracket, an add-on swing
 * faucet, and at the top the spring gooseneck, a long coil spring over the hose that keeps it arched over the sink and
 * pulls the spray valve back up out of the way. The spray valve hangs on the hose over the bowl.
 */
function preRinse(s: Shop, x: number, y: number, z: number, top: number, reach: number, wallZ: number): void {
  const p = palette();
  splashFaucet(s, x, y, z, 0, "");
  const zr = z + 0.034;
  const yr = top - reach / 2 - 0.052;
  if (s.detail === 0) {
    // Far off: the riser, the arch of the spring as a tube, the hose down to the spray valve.
    rod(s, p.chrome, 0.011, yr - y, "y", [x, (y + yr) / 2, zr]);
    const R0 = reach / 2;
    const pts: THREE.Vector3[] = [];
    for (let i = 0; i <= 6; i++) {
      const a = Math.PI - (Math.PI * i) / 6;
      pts.push(new THREE.Vector3(x, yr + 0.03 + R0 * Math.sin(a), zr + R0 + R0 * Math.cos(a)));
    }
    s.add(p.ssBright, new THREE.TubeGeometry(new THREE.CatmullRomCurve3(pts), 6, 0.019, 4, false));
    rod(s, p.rubber, 0.009, yr + 0.03 - y - 0.1, "y", [x, (yr + 0.03 + y + 0.1) / 2, zr + reach]);
    s.add(p.chrome, new THREE.BoxGeometry(0.035, 0.12, 0.035), [x, y + 0.06, zr + reach]);
    return;
  }
  s.add(p.chrome, cyl(s, 0.0105, yr - y, "y", 18), [x, (y + yr) / 2, zr]);
  s.add(p.chrome, cyl(s, 0.015, 0.04, "y", 18), [x, y + 0.02, zr]);
  // The add-on faucet's tee and swing nozzle.
  const yt = y + 0.11;
  s.add(p.chrome, cyl(s, 0.0145, 0.04, "y", 18), [x, yt, zr]);
  s.add(p.chrome, pipe(s, [[x, yt, zr], [x, yt, zr + 0.25], [x, yt - 0.03, zr + 0.27]], 0.0085, 0.02, 12));
  // The wall bracket: a clamp round the riser on an arm back to a flange on the wall.
  const yw = yr - 0.11;
  s.add(p.chrome, torus(s, 0.014, 0.0035, "y", Math.PI * 2, 24, 6), [x, yw, zr]);
  s.add(p.chrome, cyl(s, 0.006, zr - wallZ - 0.01, "z", 12), [x, yw, (zr + wallZ) / 2 - 0.005]);
  s.add(p.chrome, cyl(s, 0.026, 0.006, "z", 24), [x, yw, wallZ + 0.003]);
  if (s.fine(2)) for (const sx of [-1, 1]) screw(s, p.ssBright, [x + sx * 0.016, yw, wallZ + 0.006], 0.0025);
  // The spring-checked head at the riser's top, where the hose leaves.
  s.add(p.chrome, cyl(s, 0.014, 0.03, "y", 18), [x, yr + 0.012, zr]);
  // The gooseneck: a semicircle from the riser forward over to the hanging side, the spring coiled round it.
  const R = reach / 2;
  const cz = zr + R;
  const cy = yr + 0.03;
  const pts: THREE.Vector3[] = [];
  const nArc = s.seg(40, 12);
  for (let i = 0; i <= nArc; i++) {
    const a = Math.PI - (Math.PI * i) / nArc;
    pts.push(new THREE.Vector3(x, cy + R * Math.sin(a), cz + R * Math.cos(a)));
  }
  const arcCurve = new THREE.CatmullRomCurve3(pts);
  const coilR = 0.019;
  if (!s.fine(2)) s.add(p.ssBright, new THREE.TubeGeometry(arcCurve, 14, coilR, 6, false));
  else s.add(p.rubber, new THREE.TubeGeometry(arcCurve, 40, 0.0085, 8, false));
  // The spring, coiled at its working pitch round the hose; from where detail 1 is seen, a bright sleeve over it.
  const pitch = 0.0135;
  const turns = s.fine(2) ? Math.round(arcCurve.getLength() / pitch) : 0;
  const per = 7;
  const hp: THREE.Vector3[] = [];
  const side = new THREE.Vector3(1, 0, 0);
  const tv = new THREE.Vector3();
  const nv = new THREE.Vector3();
  for (let i = 0; turns && i <= turns * per; i++) {
    const u = i / (turns * per);
    const c = arcCurve.getPointAt(u);
    arcCurve.getTangentAt(u, tv);
    nv.crossVectors(tv, side).normalize();
    const a = (i / per) * Math.PI * 2;
    hp.push(c.clone().addScaledVector(side, Math.cos(a) * coilR).addScaledVector(nv, Math.sin(a) * coilR));
  }
  if (turns) s.add(p.ssBright, new THREE.TubeGeometry(new THREE.CatmullRomCurve3(hp), turns * per, 0.0028, 4, false));
  // The hose drops from the spring's end to the spray valve: a pistol-grip body, its squeeze lever and spray face.
  const zh = cz + R;
  const yh = y + 0.1;
  s.add(p.rubber, cyl(s, 0.0085, cy - yh - 0.13, "y", 12, 0.0085, true), [x, (cy + yh + 0.13) / 2, zh]);
  s.add(p.chrome, cyl(s, 0.011, 0.02, "y", 16, 0.011, true), [x, yh + 0.125, zh]);
  const head: [number, number][] = s.fine(2) ? [[0, 0], [0.024, 0], [0.025, 0.006], [0.019, 0.02], [0.016, 0.05], [0.016, 0.115], [0.012, 0.12], [0, 0.12]] : [[0, 0], [0.025, 0.004], [0.017, 0.03], [0.016, 0.118], [0, 0.12]];
  s.add(p.chrome, lathe(s, head, s.fine(2) ? 24 : 16), [x, yh, zh]);
  s.add(p.plastic, cyl(s, 0.0175, 0.06, "y", 20, 0.0175, true), [x, yh + 0.07, zh]);
  s.add(p.chrome, pipe(s, [[x, yh + 0.035, zh - 0.019], [x, yh + 0.06, zh - 0.034], [x, yh + 0.105, zh - 0.03], [x, yh + 0.112, zh - 0.017]], 0.0045, 0.012, s.fine(2) ? 8 : 5));
  if (s.fine(2)) {
    s.add(p.enamel, cyl(s, 0.019, 0.001, "y", 20), [x, yh - 0.0004, zh]);
    for (let i = 0; i < 8; i++) {
      const a = (i / 8) * Math.PI * 2;
      s.add(p.chrome, cyl(s, 0.0018, 0.0016, "y", 6), [x + Math.cos(a) * 0.011, yh - 0.001, zh + Math.sin(a) * 0.011]);
    }
  }
}

/** A lever waste under a bowl whose floor is at y, centred at (x, z): strainer, valve body, overflow, the lever at the front. */
function leverWaste(s: Shop, x: number, y: number, z: number, zBack: number, zFront: number, name: string): void {
  const p = palette();
  if (s.detail === 0) return;
  const f2 = s.fine(2);
  // The strainer in the bowl floor: a drawn flange round the outlet, a perforated basket down in it over the dark of
  // the valve, its lift knob in the middle.
  if (f2) {
    s.add(p.ss, lathe(s, [[0.0508, 0.0001], [0.0503, 0.0008], [0.047, 0.0024], [0.0385, 0.0024], [0.0356, 0.0012], [0.0352, -0.006], [0.033, -0.0075]], 40), [x, y, z]);
    s.add(perforatedSheet(0.0042, 0.5), new THREE.CircleGeometry(0.0335, 32).rotateX(-Math.PI / 2), [x, y - 0.0068, z]);
    s.add(p.enamel, new THREE.CircleGeometry(0.034, 24).rotateX(-Math.PI / 2), [x, y - 0.018, z]);
    s.add(p.ss, lathe(s, [[0.0065, -0.007], [0.0068, -0.002], [0.0058, 0.0042], [0.0035, 0.0058], [0, 0.006]], 16), [x, y, z]);
  } else {
    s.add(p.ss, cyl(s, 0.05, 0.0025, "y", 20), [x, y + 0.0013, z]);
    s.add(p.enamel, cyl(s, 0.036, 0.001, "y", 16), [x, y + 0.003, z]);
  }
  // The valve body under it: cast, chrome-plated, a 2 in tailpiece below.
  const yv = y - 0.045;
  const body: [number, number][] = f2 ? [[0, -0.05], [0.026, -0.05], [0.03, -0.04], [0.034, -0.02], [0.036, 0.0], [0.034, 0.025], [0.028, 0.04], [0.026, 0.045], [0, 0.045]] : [[0, -0.05], [0.03, -0.042], [0.036, 0], [0.027, 0.042], [0, 0.045]];
  s.add(p.chrome, lathe(s, body, f2 ? 24 : 16), [x, yv, z]);
  s.add(p.chrome, cyl(s, 0.0245, yv - 0.05 - 0.4, "y", f2 ? 24 : 16, 0.0245, true), [x, (yv - 0.05 + 0.4) / 2, z]);
  // The overflow, from a slotted cover high on the bowl's back wall down into the body (behind the bowl: close views).
  const yo = y + 0.2;
  if (f2) {
    // The overflow's cover: a domed stadium plate with slots, standing off the dark opening behind it.
    const cover = roundRect(0.058, 0.03, 0.015);
    for (let i = -2; i <= 2; i++) cover.holes.push(rectPath(new THREE.Path(), 0.0042, i % 2 ? 0.016 : 0.019, 0.0021, i * 0.0085, 0, true) as THREE.Path);
    s.add(p.enamel, rbox(s, 0.05, 0.024, 0.0012, 0.0006), [x, yo, zBack + 0.0006]);
    s.add(p.ss, plate(s, cover, 0.0014, 0.0006, 6), [x, yo, zBack + 0.0016]);
    for (const sx of [-1, 1]) screw(s, p.ssBright, [x + sx * 0.022, yo, zBack + 0.0036], 0.0022);
    s.add(p.chrome, pipe(s, [[x, yo, zBack - 0.012], [x, yo, zBack - 0.03], [x, yv + 0.02, zBack - 0.03], [x, yv + 0.02, z - 0.03]], 0.0125, 0.03, 12));
  }
  // The lever's rod forward to just under the rim, and the turn handle on its own joint.
  s.add(p.chrome, cyl(s, 0.006, zFront - z - 0.03, "z", 12, 0.006, true), [x, yv, (z + zFront - 0.03) / 2]);
  if (f2) s.add(p.chrome, cyl(s, 0.012, 0.03, "z", 16), [x, yv, zFront - 0.035]);
  const j = s.joint(name, [x, yv, zFront - 0.02]);
  j.add(p.chrome, cyl(j, 0.01, 0.014, "z", 16), [0, 0, 0]);
  j.add(p.chrome, f2 ? rbox(j, 0.012, 0.095, 0.012, 0.005) : new THREE.BoxGeometry(0.012, 0.095, 0.012), [0, 0.035, 0.006]);
  j.add(p.chrome, ball(j, 0.013, f2 ? 16 : 8), [0, 0.08, 0.006]);
}

interface SinkSpec {
  /** Bowl centres across. */
  bowls: number[];
  bw: number;
  bdz: number;
  depth: number;
  /** Drainboards, [x0, x1]. */
  boards: [number, number][];
  /** Rolled end rims at the left and right ends of the top. */
  faucets: { x: number; spout: number }[];
  preRinse?: number;
  legs: { x: number; z0: number; z1: number; top: number }[];
}

/**
 * A commercial sink top on legs: drawn bowls with coved corners welded into a deck, raised rolled rims along the front
 * and round the drainboards' ends, drainboards pitched toward the bowls, a coved backsplash with a 45° top, splash
 * faucets, and under each bowl its lever waste into a common 2 in drain.
 */
function sinkTop(sp: (e: { w: number; d: number; h: number }) => SinkSpec): Builder {
  return (s, e) => {
    const p = palette();
    const { d, h } = e;
    const Y = 0.914;
    const t = 0.0016;
    const S = sp(e);
    const zb0 = -0.075;
    const zb1 = zb0 - S.bdz;
    const zs = -d + 0.03;
    const bx0 = S.bowls[0]! - S.bw / 2 - 0.025;
    const bx1 = S.bowls[S.bowls.length - 1]! + S.bw / 2 + 0.025;
    const r = 0.025;
    const zc = (zb0 + zb1) / 2;
    const x0 = Math.min(bx0, ...S.boards.map((q) => q[0]));
    const x1 = Math.max(bx1, ...S.boards.map((q) => q[1]));
    const sh = h > 1.3 ? 0.254 : Math.min(0.254, h - Y - 0.002);
    // Water in the bowls: one joint at the floors, scaled to the level; and the hooks.
    const finish = (): ((b: Built) => Hooks) => {
      const water = s.joint("water", [0, Y - t - S.depth, 0]);
      for (const x of S.bowls) {
        const g = s.fine(2) ? new RoundedBoxGeometry(S.bw - 0.002, S.depth - 0.035, S.bdz - 0.002, 1, r - 0.001) : new THREE.BoxGeometry(S.bw - 0.002, S.depth - 0.035, S.bdz - 0.002);
        g.translate(0, (S.depth - 0.035) / 2, 0);
        water.add(p.water, g, [x, 0.001, zc]);
      }
      return (b) => {
        const wj = b.joints.get("water");
        if (wj) wj.visible = false;
        return {
          set: (name, v) => {
            if (name === "water" && wj) {
              wj.visible = v > 0.01;
              wj.scale.y = Math.max(0.01, Math.min(1, v));
            } else if (name === "drain") S.bowls.forEach((_, i) => turnJoint(b.joints.get(`waste_${i}`), "z", -Math.max(0, Math.min(1, v)) * (Math.PI / 2)));
          },
        };
      };
    };
    if (s.detail === 0) {
      // Far off: the bowls in a framed deck, boards, splash, legs and rails, the faucets and the pre-rinse's outline.
      for (const x of S.bowls) basin(s, p.ss, x, zc, S.bw, S.bdz, S.depth, Y - t, r);
      blk(s, p.ss, [bx0, Y - 0.04, zb0], [bx1, Y, 0], 0);
      blk(s, p.ss, [bx0, Y - 0.04, zs], [bx1, Y, zb1], 0);
      const edges = [bx0, ...S.bowls.flatMap((x) => [x - S.bw / 2, x + S.bw / 2]), bx1];
      for (let i = 0; i < edges.length; i += 2) blk(s, p.ss, [edges[i]!, Y - 0.04, zb1], [edges[i + 1]!, Y, zb0], 0);
      for (const [a, b] of S.boards) blk(s, p.ss, [a, Y - 0.04, zs], [b, Y + 0.012, 0], 0);
      blk(s, p.ss, [x0, Y, -d], [x1, Y + sh, zs], 0);
      for (const L of S.legs) {
        for (const z of [L.z0, L.z1]) gussetLeg(s, L.x, z, L.top);
        crossRail(s, [L.x, 0.2, L.z1], [L.x, 0.2, L.z0]);
      }
      S.faucets.forEach((f, i) => splashFaucet(s, f.x, Y + 0.13, zs, f.spout, `spout_${i}`));
      if (S.preRinse !== undefined) preRinse(s, S.preRinse, Y + 0.13, zs, h, Math.max(0.2, zc - zs - 0.034 + 0.06), -d);
      return finish();
    }
    // Bowls.
    for (const x of S.bowls) basin(s, p.ss, x, (zb0 + zb1) / 2, S.bw, S.bdz, S.depth, Y - t, r);
    // The rim: one formed piece along the front and both ends, its corners radiused, level all round.
    const zr = -0.001;
    const zu = zr - 0.018;
    const R = 0.025;
    s.add(p.ss, pathSweep(sheetOutline(sinkRim(s), t), filletPath([[x0, zs], [x0, zr], [x1, zr], [x1, zs]], R, s.fine(2) ? 7 : 2), Y - t / 2));
    // The deck over the bowl block, cut for the bowls, run into the rim's upstand and back to the splash's cove.
    const boardAt = (x: number) => S.boards.some(([a, b]) => Math.abs(a - x) < 1e-3 || Math.abs(b - x) < 1e-3);
    const dx0 = bx0 + (boardAt(bx0) ? 0 : 0.008);
    const dx1 = bx1 - (boardAt(bx1) ? 0 : 0.008);
    const deck = new THREE.Shape();
    rectPath(deck, dx1 - dx0, zu - (zs + 0.008), 0, (dx0 + dx1) / 2, -(zu + zs + 0.008) / 2);
    for (const x of S.bowls) deck.holes.push(rectPath(new THREE.Path(), S.bw, S.bdz, r, x, -(zb0 + zb1) / 2, true) as THREE.Path);
    const dg = plate(s, deck, t, 0, s.fine(2) ? 16 : 8);
    dg.rotateX(-Math.PI / 2);
    s.add(p.ss, dg, [0, Y - t, 0]);
    // Drainboards: pitched 15 mm up toward their outer ends, butted to the deck, run under the rim and into the splash.
    const rise = 0.015;
    const boardTop = (x: number): number | null => {
      for (const [a, b] of S.boards) {
        if (x < Math.min(a, b) || x > Math.max(a, b)) continue;
        const outer = Math.abs(a - bx0) < Math.abs(b - bx1) && b <= bx0 + 1e-6 ? a : b;
        const inner = outer === a ? b : a;
        return Y + (rise * Math.abs(x - inner)) / Math.abs(b - a);
      }
      return null;
    };
    for (const [a, b] of S.boards) {
      const outer = Math.abs(a - bx0) < Math.abs(b - bx1) && b <= bx0 + 1e-6 ? a : b;
      const inner = outer === a ? b : a;
      const sx = outer > inner ? 1 : -1;
      const L0 = Math.abs(b - a);
      const xo = outer - sx * 0.008;
      const len = Math.abs(xo - inner);
      const xm = (xo + inner) / 2;
      const ang = Math.atan2(rise, L0) * sx;
      const cz = (zs + zu) / 2;
      s.at([xm, Y - t / 2 + (rise * Math.abs(xm - inner)) / L0, cz], [0, 0, ang], () => blk(s, p.ss, [-len / 2, -t / 2, zs - cz], [len / 2, t / 2, zu - cz], 0));
      // A hat channel welded under the board, front to back, the gussets of its legs welded to it.
      const xh = outer - sx * 0.06;
      sweepZ(s, p.ss, [[-0.035, 0], [-0.018, 0], [-0.018, -0.026], [0.018, -0.026], [0.018, 0], [0.035, 0]], 0.0012, zs + 0.01, -0.05, xh, boardTop(xh)! - t - 0.0006);
    }
    // The backsplash: coved up from the deck, a 45° return to the wall at its top, its ends closed by welded caps.
    const n = s.seg(6, 2);
    const splash: P2[] = [[zs + 0.008, 0], ...arc(zs + 0.008, 0.008, 0.008, -Math.PI / 2, -Math.PI, n).slice(1), [zs, sh - 0.03], [zs - 0.026, sh], [-d + 0.0008, sh], [-d + 0.0008, sh - 0.03]];
    sweepX(s, p.ss, splash, t, x0 + t, x1 - t, Y - t / 2);
    for (const x of [x0, x1]) {
      const cap = new THREE.Shape([new THREE.Vector2(-zs, 0), new THREE.Vector2(-zs, sh - 0.03), new THREE.Vector2(d - 0.004, sh), new THREE.Vector2(d, sh), new THREE.Vector2(d, 0)]);
      const g = plate(s, cap, t);
      g.rotateY(Math.PI / 2);
      s.add(p.ss, g, [x - (x === x0 ? 0 : t), Y, 0]);
    }
    // Lever wastes and the common drain behind them.
    S.bowls.forEach((x, i) => leverWaste(s, x, Y - t - S.depth, zc, zb1, -0.001, `waste_${i}`));
    const xa = S.bowls[0]!;
    const xb = S.bowls[S.bowls.length - 1]!;
    s.add(p.chrome, pipe(s, [[xa, 0.4, zc], [xa, 0.37, zc], [xb + 0.12, 0.37, zc], [xb + 0.12, 0.37, -d + 0.003]], 0.0245, 0.06, s.fine(2) ? 20 : 8));
    for (const x of S.bowls.slice(1)) s.add(p.chrome, cyl(s, 0.03, 0.06, "y", 20), [x, 0.37, zc]);

    // Faucets and the pre-rinse.
    const yf = Y + 0.13;
    S.faucets.forEach((f, i) => splashFaucet(s, f.x, yf, zs, f.spout, `spout_${i}`));
    if (S.preRinse !== undefined) preRinse(s, S.preRinse, yf, zs, h, Math.max(0.2, zc - zs - 0.034 + 0.06), -d);

    // Legs: under the bowls' floors and under the boards' hat channels, rails between front and back.
    for (const L of S.legs) {
      const bt = boardTop(L.x);
      const top = bt === null ? L.top : bt - t - 0.0012 - 0.026 - 0.0002;
      for (const z of [L.z0, L.z1]) gussetLeg(s, L.x, z, top);
      crossRail(s, [L.x, 0.2, L.z1], [L.x, 0.2, L.z0]);
    }
    dataPlate(s, [x1 - 0.09, Y - 0.021, zr + 0.0008 + 0.0002], [0, 0, 0], 0.06, 0.022);
    return finish();
  };
}

const threeCompartment = sinkTop((e) => {
  const board = 0.26;
  const bw = (e.w - 2 * board - 2 * 0.025 - 2 * 0.03) / 3;
  const bowls = [0, 1, 2].map((i) => board + 0.025 + bw / 2 + i * (bw + 0.03));
  const zc = -0.075 - 0.26;
  return {
    bowls,
    bw,
    bdz: 0.52,
    depth: 0.305,
    boards: [
      [0, board],
      [e.w - board, e.w],
    ],
    faucets: [{ x: (bowls[1]! + bowls[2]!) / 2, spout: 0.3 }],
    preRinse: bowls[0]!,
    legs: [
      { x: 0.05, z0: -0.06, z1: -e.d + 0.08, top: 0.914 - 0.03 },
      { x: e.w - 0.05, z0: -0.06, z1: -e.d + 0.08, top: 0.914 - 0.03 },
      { x: board + 0.075, z0: zc + 0.2, z1: zc - 0.2, top: 0.914 - 0.0016 - 0.305 - 0.0016 },
      { x: e.w - board - 0.075, z0: zc + 0.2, z1: zc - 0.2, top: 0.914 - 0.0016 - 0.305 - 0.0016 },
    ],
  };
});

const vegPrepSink = sinkTop((e) => {
  const bw = 0.56;
  const bowls = [0.03 + 0.025 + bw / 2];
  const zc = -0.075 - 0.26;
  return {
    bowls,
    bw,
    bdz: 0.52,
    depth: 0.33,
    boards: [[0.03 + 0.05 + bw, e.w]],
    faucets: [{ x: bowls[0]!, spout: 0.25 }],
    legs: [
      { x: 0.1, z0: zc + 0.2, z1: zc - 0.2, top: 0.914 - 0.0032 - 0.33 },
      { x: 0.03 + 0.05 + bw - 0.07, z0: zc + 0.2, z1: zc - 0.2, top: 0.914 - 0.0032 - 0.33 },
      { x: e.w - 0.05, z0: -0.06, z1: -e.d + 0.08, top: 0.914 - 0.03 },
    ],
  };
});

/* ------------------------------------------------------------------------------------------------------------ */
/* Work tables                                                                                                    */
/* ------------------------------------------------------------------------------------------------------------ */

/**
 * A 16-gauge 304 work table: the top, its edge a boxed bull-nose or a raised marine lip formed all round in one piece
 * with radiused corners (or round the front and sides, with a coved backsplash and a 45° return closed by welded
 * end caps), hat channels under it front to back that the leg gussets weld to, 1-5/8 in legs on bullet feet, and an
 * undershelf, its edges turned down and hemmed, welded to flared gussets that clamp the legs.
 */
function workTable(opts: { legs: 4 | 6; edge: "bull" | "marine"; splash: number; shelf: "galv" | "ss" }): Builder {
  return (s, e) => {
    const p = palette();
    const { w, d } = e;
    const Y = 0.914;
    const t = 0.0016;
    const zs = -d + 0.025;
    const legX = opts.legs === 6 ? [0.055, w / 2, w - 0.055] : [0.055, w - 0.055];
    const sm = opts.shelf === "galv" ? p.ssDull : p.ss;
    const ys = 0.254;
    const ts = 0.0012;
    const sx0 = 0.024;
    const sx1 = w - 0.024;
    const sz0 = -d + 0.028;
    const sz1 = -0.028;
    if (s.detail === 0) {
      // Far off: top and edge as one slab, the splash, legs to the floor and the undershelf.
      blk(s, p.ss, [0, Y - 0.04, -d], [w, Y, 0], 0);
      if (opts.splash) blk(s, p.ss, [0, Y, -d], [w, Y + opts.splash, -d + 0.025], 0);
      for (const x of legX) for (const z of [-0.06, -d + 0.06]) gussetLeg(s, x, z, Y - 0.04);
      blk(s, sm, [sx0, ys - 0.032, sz0], [sx1, ys, sz1], 0);
      return;
    }
    // The edge, one formed piece round the plan, and the top sheet butted inside it.
    const fine = s.fine(2);
    const back = opts.edge === "bull" ? 0.02 : 0.03;
    const R = 0.045;
    const fn = fine ? 7 : 2;
    const edge = sheetOutline((opts.edge === "bull" ? bullNose(s, back) : marineEdge(s, back)).map(([z, y]) => [-z, y] as P2), t);
    const plan: P2[] = opts.splash ? [[0, -d], [0, 0], [w, 0], [w, -d]] : [[0, 0], [w, 0], [w, -d], [0, -d]];
    const edgePath = filletPath(plan, R, fn, !opts.splash);
    s.add(p.ss, pathSweep(edge, edgePath, Y - t / 2, !opts.splash));
    const topLine = opts.splash ? offsetPath(filletPath([[0, zs + 0.008], [0, 0], [w, 0], [w, zs + 0.008]], R, fn), back) : offsetPath(edgePath, back, true);
    const tg = plate(s, planShape(topLine), t);
    tg.rotateX(-Math.PI / 2);
    s.add(p.ss, tg, [0, Y - t, 0]);
    if (opts.splash) {
      // The splash: coved up from the top, a 45° return to the wall, its ends closed by welded caps.
      const sh = opts.splash;
      const n = s.seg(6, 2);
      sweepX(s, p.ss, [[zs + 0.008, 0], ...arc(zs + 0.008, 0.008, 0.008, -Math.PI / 2, -Math.PI, n).slice(1), [zs, sh - 0.025], [zs - 0.024, sh - 0.001], [-d + 0.0008, sh - 0.001], [-d + 0.0008, sh - 0.025]], t, t, w - t, Y - t / 2);
      for (const x of [0, w - t]) {
        const cap = new THREE.Shape([[-zs, 0], [-zs, sh - 0.025], [d - 0.001, sh], [d, sh], [d, 0]].map(([u, v]) => new THREE.Vector2(u!, v!)));
        const g = plate(s, cap, t);
        g.rotateY(Math.PI / 2);
        s.add(p.ss, g, [x, Y - t, 0]);
      }
    }

    // Hat channels under the top, front to back, at every pair of legs and between.
    const hats = [...legX];
    for (let i = 0; i < legX.length - 1; i++) {
      const span = legX[i + 1]! - legX[i]!;
      const k = Math.round(span / 0.6);
      for (let j = 1; j < k; j++) hats.push(legX[i]! + (span * j) / k);
    }
    const hatY = Y - t;
    const hz0 = -d + 0.045;
    for (const x of hats) sweepZ(s, p.ss, [[-0.04, 0], [-0.02, 0], [-0.02, -0.03], [0.02, -0.03], [0.02, 0], [0.04, 0]], 0.0012, hz0, -0.045, x, hatY - 0.0006);
    // Legs, gussets welded under the hat channels.
    const zf = -0.06;
    const zr = -d + 0.06;
    for (const x of legX) for (const z of [zf, zr]) gussetLeg(s, x, z, hatY - 0.031);

    // The undershelf: 18-gauge, its edges turned down and hemmed under, mitred at the corners, on the legs' gussets.
    const nb = fine ? 4 : 1;
    const lip = sheetOutline([[0.02, 0], [0.004, 0], ...arc(0.004, -0.004, 0.004, Math.PI / 2, Math.PI, nb).slice(1), [0, -0.029], ...arc(0.003, -0.029, 0.003, Math.PI, Math.PI * 1.5, nb).slice(1), [0.012, -0.032]], ts);
    const shelfPath: P2[] = [[sx0, sz1], [sx1, sz1], [sx1, sz0], [sx0, sz0]];
    s.add(sm, pathSweep(lip, shelfPath, ys - ts / 2, true));
    const sg = plate(s, planShape(offsetPath(shelfPath, 0.02, true)), ts);
    sg.rotateX(-Math.PI / 2);
    s.add(sm, sg, [0, ys - ts, 0]);
    // Each gusset: a drawn socket flared up to a flange welded under the shelf, pinned to the leg by a set screw.
    const gp: [number, number][] = fine
      ? [[LEG_R + 0.0012, 0], [LEG_R + 0.002, 0.003], [LEG_R + 0.0025, 0.022], [0.0245, 0.031], [0.0272, 0.0345], [0.028, 0.0365], [0.028, 0.0378], [0, 0.0378]]
      : [[LEG_R + 0.0015, 0], [LEG_R + 0.0025, 0.022], [0.028, 0.035], [0, 0.0378]];
    for (const x of legX)
      for (const z of [zf, zr]) {
        s.add(sm, lathe(s, gp, fine ? 24 : 14), [x, ys - ts - 0.0378, z]);
        if (fine) {
          const sz = z > -d / 2 ? -1 : 1;
          s.add(p.ssBright, hex(0.0042, 0.006, "z"), [x, ys - ts - 0.026, z + sz * (LEG_R + 0.0045)]);
          s.add(p.ssBright, cyl(s, 0.0021, 0.004, "z", 8), [x, ys - ts - 0.026, z + sz * (LEG_R + 0.009)]);
        }
      }
    if (w > 1.5) sweepX(s, sm, [[-0.03, 0], [-0.015, 0], [-0.015, -0.022], [0.015, -0.022], [0.015, 0], [0.03, 0]], 0.0012, sx0 + 0.03, sx1 - 0.03, ys - ts - 0.0006, -d / 2);
    // The listing plate on the web of the right-hand hat channel.
    const xh = Math.max(...hats);
    dataPlate(s, [xh + 0.02 + 0.0008, hatY - 0.016, -d / 2], [0, Math.PI / 2, 0], 0.06, 0.022);
  };
}


let paperM: THREE.MeshStandardMaterial | null = null;
/** Unbleached paper towel: warm white, matt, faintly embossed, thin enough to light through. */
function towelMat(): THREE.MeshStandardMaterial {
  if (paperM) return paperM;
  const c = document.createElement("canvas");
  c.width = c.height = 64;
  const x = c.getContext("2d")!;
  x.fillStyle = "#808080";
  x.fillRect(0, 0, 64, 64);
  for (let i = 0; i < 64; i += 4) for (let j = (i / 4) % 2 ? 2 : 0; j < 64; j += 4) {
    x.fillStyle = "#9a9a9a";
    x.fillRect(j, i, 2, 2);
  }
  const t = new THREE.CanvasTexture(c);
  t.wrapS = t.wrapT = THREE.RepeatWrapping;
  t.repeat.set(6, 6);
  paperM = new THREE.MeshStandardMaterial({ color: 0xeeeae0, roughness: 0.96, metalness: 0, bumpMap: t, bumpScale: 0.5, side: THREE.DoubleSide });
  paperM.userData.uv = "own";
  return paperM;
}
let stackM: THREE.MeshStandardMaterial | null = null;
/** The edges of a stack of folded towels seen through a dispenser's window: fine uneven layers, warm white. */
function towelStack(): THREE.MeshStandardMaterial {
  if (stackM) return stackM;
  const c = document.createElement("canvas");
  c.width = 16;
  c.height = 64;
  const x = c.getContext("2d")!;
  // Each folded towel is a band, lit at its rounded fold and shadowed where the next one presses on it.
  let y = 0;
  let k = 7;
  while (y < 64) {
    const hgt = 4 + ((k * 37) % 4) * 0.5;
    const g = x.createLinearGradient(0, y, 0, y + hgt);
    g.addColorStop(0, "#fbf9f3");
    g.addColorStop(0.55, "#e9e4d8");
    g.addColorStop(1, "#8f877a");
    x.fillStyle = g;
    x.fillRect(0, y, 16, hgt);
    y += hgt;
    k = (k * 13 + 5) % 97;
  }
  const t = new THREE.CanvasTexture(c);
  t.colorSpace = THREE.SRGBColorSpace;
  stackM = new THREE.MeshStandardMaterial({ map: t, roughness: 0.95, metalness: 0 });
  stackM.userData.uv = "own";
  return stackM;
}
let paneM: THREE.MeshPhysicalMaterial | null = null;
/** A dispenser's level window: a lightly smoked acrylic pane, glossy, mostly clear. */
function windowPane(): THREE.MeshPhysicalMaterial {
  if (!paneM) {
    paneM = new THREE.MeshPhysicalMaterial({ color: 0x4a4f55, roughness: 0.06, metalness: 0, transparent: true, opacity: 0.32, clearcoat: 1, depthWrite: false });
    paneM.userData.uv = "own";
  }
  return paneM;
}
/** A paper towel `w` wide hanging `len` from a dispenser's slot at the origin, curling forward as it falls. */
function towelTail(s: Shop, w: number, len: number): THREE.BufferGeometry {
  const ny = s.fine(2) ? 10 : 3;
  const g = new THREE.PlaneGeometry(w, len, s.fine(2) ? 6 : 1, ny);
  const pos = g.attributes.position as THREE.BufferAttribute;
  for (let i = 0; i < pos.count; i++) {
    const u = (len / 2 - pos.getY(i)) / len;
    const xx = pos.getX(i) / (w / 2);
    // Down from the slot, its lower edge kicking forward; a soft ripple across where the folds were.
    // A multifold towel: its panel creased down the middle and its lower edge curling out.
    const crease = 0.006 * (1 - Math.abs(xx)) * Math.min(1, u * 3);
    pos.setXYZ(i, pos.getX(i) * (1 - 0.08 * u), -u * len * 0.95, 0.032 * u * u * u + crease + 0.0035 * Math.sin(xx * Math.PI * 2.5 + 0.7) * u);
  }
  g.computeVertexNormals();
  return g;
}

/**
 * A falling column of water `len` long down from y = 0: it leaves the aerator 13 mm across and necks as it speeds up,
 * then surface tension pinches it into a string of growing beads, and the last of it falls as separate drops.
 */
function waterStream(s: Shop, m: Mat, len: number): void {
  const r0 = 0.0062;
  const v0 = 0.65;
  const neck = (y: number) => r0 * Math.pow((v0 * v0) / (v0 * v0 + 2 * 9.81 * y), 0.25);
  const n = s.fine(2) ? 64 : 5;
  const intact = len * 0.86;
  const prof: [number, number][] = [[0, 0.0005], [r0 * 1.02, 0.0005]];
  for (let i = 1; i <= n; i++) {
    const y = (intact * i) / n;
    const grow = Math.max(0, (y - len * 0.45) / (intact - len * 0.45));
    const ripple = s.fine(2) ? 1 + 0.34 * grow * grow * Math.sin((y / (9 * neck(y))) * Math.PI * 2) : 1;
    prof.push([neck(y) * ripple, -y]);
  }
  prof.push([neck(intact) * 0.5, -intact - 0.002], [0, -intact - 0.0035]);
  s.add(m, hang(s, prof, s.fine(2) ? 20 : 10));
  // The drops it has shed, between the column's end and the strainer.
  const drops = s.fine(2) ? 4 : 1;
  for (let i = 0; i < drops; i++) {
    const y = intact + 0.008 + ((len - intact - 0.012) * (i + 0.5)) / drops;
    const rd = neck(y) * (1.05 - 0.12 * i);
    const g = ball(s, rd, s.fine(2) ? 12 : 6);
    g.scale(1, 1.18, 1);
    s.add(m, g, [0.0004 * (i % 2 ? 1 : -1), -y, 0]);
  }
}

/**
 * A wall-hung hand sink: a drawn bowl in a box skirt with a rolled rim, side splash guards, a backsplash with a 45° top,
 * a gooseneck faucet on 4 in centres with wrist blades, a soap dispenser and a towel dispenser on the splash, and the
 * hanger bracket on the wall it hooks over.
 */
const handSink: Builder = (s, e) => {
  const p = palette();
  const { w, d, h } = e;
  const Y = 0.165;
  const t = 0.0012;
  const zs = -d + 0.03;
  const bw = 0.33;
  const bdz = 0.245;
  const zb0 = -0.045;
  const zb1 = zb0 - bdz;
  const bowlDepth = 0.125;
  const r = 0.03;
  const hooks = (b: Built): Hooks => {
    const st = b.joints.get("stream");
    if (st) st.visible = false;
    return {
      set: (name, v) => {
        if (name === "flow" && st) {
          st.visible = v > 0.02;
          const k = 0.5 + Math.min(1, v);
          st.scale.set(k, 1, k);
        }
      },
    };
  };
  const xf = w / 2 - 0.02;
  const zc = (zb0 + zb1) / 2;
  const yf = Y - t - bowlDepth;
  const R = 0.07;
  const top = 0.395;
  const yo = top - R - 0.06;
  const zb = zs + 0.03;
  if (s.detail === 0) {
    // Far off: the bowl in its box, splash and side guards, the faucet's outline and its stream, the dispensers.
    basin(s, p.ss, w / 2, zc, bw, bdz, bowlDepth, Y - t, r, 0.0012);
    blk(s, p.ss, [0, Y - 0.115, zb0], [w, Y, 0], 0);
    blk(s, p.ss, [0, Y - 0.115, zs], [w, Y, zb1], 0);
    for (const [a, b] of [[0, w / 2 - bw / 2], [w / 2 + bw / 2, w]]) blk(s, p.ss, [a!, Y - 0.115, zb1], [b!, Y, zb0], 0);
    blk(s, p.ss, [0, 0, -d], [w, h, zs], 0);
    for (const x of [0.0006, w - 0.0018]) blk(s, p.ss, [x, Y, zs], [x + 0.0012, Y + 0.12, -0.03], 0);
    s.add(p.chrome, new THREE.BoxGeometry(0.13, 0.025, 0.03), [xf, Y + 0.07, zb]);
    const pts = [new THREE.Vector3(xf, Y + 0.07, zb), new THREE.Vector3(xf, top - 0.01, zb + 0.03), new THREE.Vector3(xf, top - 0.03, zb + 2 * R), new THREE.Vector3(xf, yo, zb + 2 * R)];
    s.add(p.chrome, new THREE.TubeGeometry(new THREE.CatmullRomCurve3(pts), 6, 0.01, 4, false));
    const stream = s.joint("stream", [xf, yo, zb + 2 * R]);
    const sg = new THREE.CylinderGeometry(0.0065, 0.0055, yo - yf, 6, 1, true);
    sg.translate(0, -(yo - yf) / 2, 0);
    stream.add(p.water, sg);
    blk(s, p.ss, [0.035, 0.3, zs + 0.001], [0.13, 0.448, zs + 0.08], 0);
    blk(s, p.ss, [0.262, 0.3, zs + 0.001], [0.424, 0.452, zs + 0.09], 0);
    return hooks;
  }
  basin(s, p.ss, w / 2, zc, bw, bdz, bowlDepth, Y - t, r, 0.0012);
  // The rolled rim along the front, its corners formed round to short returns where the side guards begin.
  const zr = -0.001;
  const zg = -0.03;
  const rimPath = filletPath([[0, zg], [0, zr], [w, zr], [w, zg]], 0.02, s.fine(2) ? 7 : 2);
  s.add(p.ss, pathSweep(sheetOutline(sinkRim(s), t), rimPath, Y - t / 2));
  // The deck round the bowl: run into the rim's upstand in front, out to the side guards, back to the splash's cove.
  const xg = 0.0022;
  const deck = planShape([[xg, zs + 0.008], [xg, zg], ...offsetPath(rimPath, 0.018).slice(1, -1), [w - xg, zg], [w - xg, zs + 0.008]]);
  deck.holes.push(rectPath(new THREE.Path(), bw, bdz, r, w / 2, -(zb0 + zb1) / 2, true) as THREE.Path);
  const dg = plate(s, deck, t, 0, s.fine(2) ? 16 : 8);
  dg.rotateX(-Math.PI / 2);
  s.add(p.ss, dg, [0, Y - t, 0]);
  // The box skirt under it, set in behind the rim's apron, and its turned-in bottom flange.
  const skirt = 0.115;
  const zk = zr - 0.002;
  blk(s, p.ss, [0.001, Y - skirt, zs], [xg, Y, zk], 0);
  blk(s, p.ss, [w - xg, Y - skirt, zs], [w - 0.001, Y, zk], 0);
  blk(s, p.ss, [xg, Y - skirt, zk - 0.0012], [w - xg, Y - 0.044, zk], 0.0004);
  blk(s, p.ss, [xg, Y - skirt - 0.004, zk - 0.012], [w - xg, Y - skirt, zk], 0.0004);
  blk(s, p.ss, [0.001, Y - skirt - 0.004, zs], [0.012, Y - skirt, zk - 0.012], 0.0004);
  blk(s, p.ss, [w - 0.012, Y - skirt - 0.004, zs], [w - 0.001, Y - skirt, zk - 0.012], 0.0004);
  // Backsplash to the top of the box, coved up from the deck, with its return to the wall.
  const sh = h - Y;
  const n = s.seg(6, 2);
  sweepX(s, p.ss, [[zs + 0.008, 0], ...arc(zs + 0.008, 0.008, 0.008, -Math.PI / 2, -Math.PI, n).slice(1), [zs, sh - 0.022], [zs - 0.02, sh - 0.001], [-d + 0.0008, sh - 0.001], [-d + 0.0008, sh - 0.025]], t, xg, w - xg, Y - t / 2);
  blk(s, p.ss, [0.004, 0.0, -d], [w - 0.004, Y - t, -d + 0.0012], 0.0004);
  // Side splash guards: welded trapezoids rising from the rim's returns to the splash, their top edges hemmed round.
  for (const x of [0.001, w - xg]) {
    const g = plate(s, new THREE.Shape([new THREE.Vector2(-zg, 0), new THREE.Vector2(-zg, 0.04), new THREE.Vector2(-zs, 0.21), new THREE.Vector2(d, 0.21), new THREE.Vector2(d, 0)]), 0.0012);
    g.rotateY(Math.PI / 2);
    s.add(p.ss, g, [x, Y, 0]);
    s.add(p.ss, pipe(s, [[x + 0.0006, Y + 0.04, zg], [x + 0.0006, Y + 0.21, zs]], 0.003, 0.01, 8));
  }
  // The drain: basket strainer, tailpiece and slip nut.
  s.add(p.ssBright, cyl(s, 0.032, 0.0022, "y", 28), [w / 2, yf + 0.0011, zc]);
  s.add(p.enamel, cyl(s, 0.022, 0.001, "y", 24), [w / 2, yf + 0.0025, zc]);
  if (s.fine(2)) for (let i = 0; i < 6; i++) s.add(p.ssBright, rbox(s, 0.04, 0.0015, 0.0035, 0.0006), [w / 2, yf + 0.003, zc], [0, (i * Math.PI) / 6, 0]);
  s.add(p.chrome, cyl(s, 0.019, yf - 0.004, "y", 20, 0.019, true), [w / 2, yf / 2 + 0.002, zc]);
  s.add(p.chrome, cyl(s, 0.024, 0.022, "y", 20), [w / 2, yf - 0.018, zc]);
  s.add(p.chrome, hex(0.026, 0.014), [w / 2, 0.016, zc]);
  // The hanger bracket on the wall, the sink's back hooked over its top lip, and its lag bolts.
  blk(s, p.ssDull, [0.07, Y - 0.115, -d + 0.0012], [w - 0.07, Y - 0.04, -d + 0.0042], 0.0006);
  blk(s, p.ssDull, [0.07, Y - 0.045, -d + 0.004], [w - 0.07, Y - 0.04, -d + 0.016], 0.0006);
  if (s.fine(1)) for (const x of [0.1, w / 2, w - 0.1]) s.add(p.ssBright, hex(0.008, 0.006, "z"), [x, Y - 0.085, -d + 0.0072]);

  // The faucet: 4 in centres on the splash, wrist blades, a gooseneck over the bowl.
  const yv = Y + 0.07;
  for (const sx of [-1, 1]) {
    const x = xf + sx * 0.0508;
    s.add(p.chrome, cyl(s, 0.009, 0.03, "z", 14, 0.009, true), [x, yv, zs + 0.015]);
    const bonnet: [number, number][] = s.fine(2) ? [[0.014, -0.01], [0.015, 0.01], [0.011, 0.017], [0.007, 0.02], [0.007, 0.034], [0, 0.034]] : [[0.015, -0.01], [0.014, 0.015], [0.007, 0.02], [0, 0.034]];
    s.add(p.chrome, lathe(s, bonnet, 18), [x, yv, zb]);
    s.at([x, yv + 0.036, zb], [0, sx * 0.6, 0], () => {
      s.add(p.chrome, cyl(s, 0.0095, 0.01, "y", 14), [0, 0, 0]);
      s.add(p.chrome, taperBlade(s, 0.085, 0.014, 0.019, 0.0035), [0, 0.003, 0], [0, -Math.PI / 2, 0]);
      if (s.fine(2)) s.add(sx < 0 ? p.ledRed : p.enamel, cyl(s, 0.0055, 0.001, "y", 12), [0, 0.0055, 0]);
    });
  }
  s.add(p.chrome, cyl(s, 0.012, 0.1016, "x", 18), [xf, yv, zb]);
  s.add(p.chrome, cyl(s, 0.013, 0.03, "y", 18), [xf, yv + 0.02, zb]);
  const goose: V3[] = [[xf, yv + 0.03, zb], [xf, top - R, zb]];
  for (let i = 1; i <= 8; i++) {
    const a = Math.PI - (Math.PI * i) / 8;
    goose.push([xf, top - R + R * Math.sin(a), zb + R + R * Math.cos(a)]);
  }
  goose.push([xf, top - R - 0.045, zb + 2 * R]);
  s.add(p.chrome, pipe(s, goose, 0.009, 0.012, s.fine(2) ? 14 : 8));
  s.add(p.chrome, cyl(s, 0.0115, 0.02, "y", 16), [xf, top - R - 0.05, zb + 2 * R]);
  if (s.fine(2)) s.add(p.enamel, cyl(s, 0.0075, 0.001, "y", 14), [xf, yo - 0.0003, zb + 2 * R]);
  // The running water, a joint scaled by the flow: a column leaving the aerator, necking as it falls faster, then
  // rippling into beads as surface tension pinches it, and breaking into drops just above the strainer.
  const stream = s.joint("stream", [xf, yo, zb + 2 * R]);
  waterStream(stream, p.water, yo - yf);

  // Soap dispenser, upper left: a stainless cartridge housing, a level window, the push plate and its nozzle.
  const zd = zs;
  {
    const [x0, x1, y0, y1] = [0.035, 0.13, 0.3, 0.448];
    blk(s, p.ss, [x0, y0, zd + 0.001], [x1, y1, zd + 0.08], 0.008);
    blk(s, p.glassSmoke, [x0 + 0.02, y0 + 0.07, zd + 0.079], [x1 - 0.02, y1 - 0.025, zd + 0.0815], 0.002);
    blk(s, p.plastic, [x0 + 0.01, y0 + 0.008, zd + 0.06], [x1 - 0.01, y0 + 0.055, zd + 0.088], 0.006);
    s.add(p.plastic, cyl(s, 0.0045, 0.012, "y", 10), [(x0 + x1) / 2, y0 - 0.002, zd + 0.05]);
    if (s.fine(2)) s.add(p.ssBright, cyl(s, 0.004, 0.003, "z", 12), [(x0 + x1) / 2, y1 - 0.012, zd + 0.081]);
  }
  // Towel dispenser, upper right: a stainless cabinet with a key lock and a towel at its slot.
  {
    const [x0, x1, y0, y1] = [0.262, 0.424, 0.3, 0.452];
    blk(s, p.ss, [x0, y0, zd + 0.001], [x1, y1, zd + 0.09], 0.01);
    blk(s, p.ssDull, [x0 + 0.03, y0 - 0.0005, zd + 0.03], [x1 - 0.03, y0 + 0.002, zd + 0.06], 0.001);
    // The next towel hanging from the slot, its fold curling toward the user.
    s.add(towelMat(), towelTail(s, x1 - x0 - 0.075, 0.055), [(x0 + x1) / 2, y0 - 0.0005, zd + 0.045]);
    // The level window: the edges of the stacked towels behind a lightly tinted pane in a pressed bezel.
    const [wx0, wx1, wy0, wy1] = [x0 + 0.055, x1 - 0.055, y0 + 0.028, y0 + 0.072];
    if (s.fine(2)) {
    s.add(towelStack(), new THREE.PlaneGeometry(wx1 - wx0, wy1 - wy0), [(wx0 + wx1) / 2, (wy0 + wy1) / 2, zd + 0.0906]);
    const bez = roundRect(wx1 - wx0 + 0.012, wy1 - wy0 + 0.012, 0.006);
    bez.holes.push(rectPath(new THREE.Path(), wx1 - wx0, wy1 - wy0, 0.003, 0, 0, true) as THREE.Path);
    s.add(p.ss, plate(s, bez, 0.0012, s.fine(2) ? 0.0004 : 0, 4), [(wx0 + wx1) / 2, (wy0 + wy1) / 2, zd + 0.09]);
    s.add(windowPane(), rbox(s, wx1 - wx0 + 0.002, wy1 - wy0 + 0.002, 0.001, 0.0004), [(wx0 + wx1) / 2, (wy0 + wy1) / 2, zd + 0.0914]);
    }
    s.add(p.ssBright, cyl(s, 0.006, 0.004, "z", 14), [(x0 + x1) / 2, y1 - 0.02, zd + 0.091]);
    if (s.fine(2)) s.add(p.enamel, rbox(s, 0.0015, 0.006, 0.001, 0.0003), [(x0 + x1) / 2, y1 - 0.02, zd + 0.0932]);
  }
  dataPlate(s, [w - 0.06, Y - 0.075, zk + 0.0002], [0, 0, 0], 0.05, 0.022);
  return hooks;
};


/* ------------------------------------------------------------------------------------------------------------ */
/* Dishmachines                                                                                                   */
/* ------------------------------------------------------------------------------------------------------------ */

/**
 * A dial thermometer on a panel facing +z at `at`: chrome bezel, white face, a scale of ticks over 270°, and its
 * needle on a joint `name` (turned by `needleAngle`).
 */
function dialGauge(s: Shop, at: V3, r: number, name: string): void {
  const p = palette();
  s.at(at, [0, 0, 0], () => {
    s.add(p.chrome, lathe(s, [[r * 0.92, 0], [r * 1.12, 0], [r * 1.14, 0.006], [r * 1.05, 0.011], [r * 0.98, 0.011]], 32).rotateX(Math.PI / 2));
    s.add(p.hdpe, cyl(s, r, 0.004, "z", 32), [0, 0, 0.003]);
    if (s.fine(1))
      for (let i = 0; i <= 12; i++) {
        const a = Math.PI * 1.25 - (i / 12) * Math.PI * 1.5;
        const len = i % 3 ? r * 0.1 : r * 0.2;
        s.add(p.enamel, new THREE.BoxGeometry(r * 0.035, len, 0.0006), [Math.cos(a) * (r * 0.82 - len / 2), Math.sin(a) * (r * 0.82 - len / 2), 0.0053], [0, 0, a - Math.PI / 2]);
      }
    if (s.fine(2)) s.add(p.ledRed, rbox(s, r * 0.5, r * 0.06, 0.0005, 0.0002), [0, -r * 0.45, 0.0053]);
    s.add(p.glass, cyl(s, r * 0.99, 0.0015, "z", 32), [0, 0, 0.0095]);
    const j = s.joint(name, [0, 0, 0.006]);
    j.add(p.enamel, taperBlade(j, r * 0.78, r * 0.08, r * 0.025, 0.0008), [0, 0, 0.0005], [Math.PI / 2, 0, Math.PI / 2]);
    j.add(p.enamel, cyl(j, r * 0.09, 0.002, "z", 12), [0, 0, 0.001]);
  });
}
/** Needle angle (from straight up, counter-clockwise positive) for a reading `v` on a 0..`max` °C dial over 270°. */
function needleAngle(v: number, max: number): number {
  const c = v > 200 ? v - 273.15 : v;
  return Math.PI * 0.75 - Math.max(0, Math.min(1, c / max)) * Math.PI * 1.5;
}

/**
 * A wash arm on its hub: two flattened tubes, their jets slotted on the side facing the rack, their ends canted so the
 * jets' reaction spins the arm; a thinner rinse arm crosses it. Built in a joint at the hub, jets toward `dir` (±1 y).
 */
function washArms(j: Shop, span: number, dir: number): void {
  const p = palette();
  j.add(p.ss, cyl(j, 0.028, 0.03, "y", 24), [0, 0, 0]);
  for (const sx of [-1, 1]) {
    j.add(p.ss, rbox(j, span / 2 - 0.02, 0.022, 0.04, 0.009), [sx * (span / 4 + 0.01), 0, 0]);
    j.add(p.ss, rbox(j, 0.02, 0.022, 0.032, 0.009), [sx * (span / 2 - 0.004), 0, 0], [0, sx * 0.5, 0]);
    if (j.fine(1))
      for (let i = 0; i < 6; i++) j.add(p.enamel, new THREE.BoxGeometry(0.012, 0.001, 0.003), [sx * (0.05 + i * ((span / 2 - 0.08) / 5)), dir * 0.0112, (i % 2 ? 1 : -1) * 0.008], [0, 0.4 * sx, 0]);
    j.add(p.ssBright, cyl(j, 0.008, span / 2 - 0.05, "z", 12), [0, -dir * 0.002, sx * (span / 4 + 0.012)]);
    if (j.fine(2)) for (let i = 0; i < 4; i++) j.add(p.ssBright, cyl(j, 0.003, 0.006, "y", 8), [0, dir * 0.009, sx * (0.06 + i * 0.05)]);
  }
}

let china: THREE.MeshPhysicalMaterial | null = null;
/** Vitrified hotel china: a warm white body under a hard clear glaze that holds sharp reflections. */
function chinaMat(): THREE.MeshPhysicalMaterial {
  if (!china) {
    china = new THREE.MeshPhysicalMaterial({ color: 0xf2efe7, metalness: 0, roughness: 0.32, clearcoat: 1, clearcoatRoughness: 0.05 });
    china.userData.uv = "box";
  }
  return china;
}

/** A 20 × 20 in peg rack, its floor at y, centred at (x, z), loaded with `plates` dinner plates standing in its pegs. */
function dishRack(s: Shop, x: number, y: number, z: number, plates: number): void {
  const p = palette();
  const S = 0.5;
  const H = 0.1;
  const c = 0.018;
  // The moulded polypropylene frame: corner posts, rails top and bottom butted between them, an open lattice of
  // uprights in each wall, a floor of deep ribs one way and shallow ones the other, and two rows of pegs.
  for (const sx of [-1, 1]) for (const sz of [-1, 1]) blk(s, p.poly, [x + sx * (S / 2 - c), y, z + sz * (S / 2 - c)], [x + (sx * S) / 2, y + H, z + (sz * S) / 2], 0.004);
  const rt = 0.012;
  for (const yy of [y, y + H - rt])
    for (const sz of [-1, 1]) {
      blk(s, p.poly, [x - S / 2 + c, yy, z + sz * (S / 2 - rt)], [x + S / 2 - c, yy + rt, z + (sz * S) / 2], 0.003);
      blk(s, p.poly, [x + sz * (S / 2 - rt), yy, z - S / 2 + c], [x + (sz * S) / 2, yy + rt, z + S / 2 - c], 0.003);
    }
  const nu = s.fine(1) ? 7 : 0;
  for (let i = 1; i <= nu; i++) {
    const q = -S / 2 + c + ((S - 2 * c) * i) / (nu + 1);
    for (const sz of [-1, 1]) {
      blk(s, p.poly, [x + q - 0.004, y + rt, z + sz * (S / 2 - 0.008)], [x + q + 0.004, y + H - rt, z + (sz * S) / 2 - sz * 0.001], 0);
      blk(s, p.poly, [x + sz * (S / 2 - 0.008), y + rt, z + q - 0.004], [x + (sz * S) / 2 - sz * 0.001, y + H - rt, z + q + 0.004], 0);
    }
  }
  const nl = s.fine(1) ? 9 : 4;
  for (let i = 1; i < nl; i++) {
    const q = -S / 2 + (S * i) / nl;
    blk(s, p.poly, [x - S / 2 + rt, y + 0.001, z + q - 0.003], [x + S / 2 - rt, y + 0.009, z + q + 0.003], 0);
    blk(s, p.poly, [x + q - 0.002, y + 0.002, z - S / 2 + rt], [x + q + 0.002, y + 0.0065, z + S / 2 - rt], 0);
  }
  // Pegs in two rows across, and the plates standing between them, facing front.
  const pitch = 0.045;
  const n = Math.min(plates, 10);
  if (s.fine(1))
    for (let i = 0; i <= 10; i++)
      for (const px of [-0.09, 0.09]) s.add(p.poly, cyl(s, 0.0035, 0.06, "y", 8, 0.002), [x + px, y + 0.039, z - S / 2 + 0.03 + i * pitch]);
  // A 10-1/2 in dinner plate: a foot ring under it, the well, the marly rising to a rolled rim.
  const prof: [number, number][] = s.fine(2)
    ? [[0, 0.0035], [0.058, 0.0035], [0.061, 0.0012], [0.064, 0], [0.068, 0], [0.071, 0.0022], [0.085, 0.0068], [0.105, 0.0125], [0.124, 0.0172], [0.1305, 0.0196], [0.1325, 0.0222], [0.1312, 0.0243], [0.1272, 0.0246], [0.112, 0.0198], [0.096, 0.0142], [0.088, 0.0102], [0.08, 0.0082], [0.06, 0.0072], [0, 0.0069]]
    : [[0, 0.0035], [0.064, 0], [0.071, 0.0022], [0.124, 0.0172], [0.1325, 0.0222], [0.1272, 0.0246], [0.096, 0.0142], [0.08, 0.0082], [0, 0.0069]];
  const plateG = lathe(s, prof, s.fine(2) ? 56 : 24);
  plateG.rotateX(Math.PI / 2);
  for (let i = 0; i < n; i++) s.add(chinaMat(), plateG.clone(), [x, y + 0.009 + 0.1325, z - S / 2 + 0.03 + (i + 0.5) * pitch], [-0.12, 0, 0]);
}

const doorDishwasher: Builder = (s, e) => {
  const p = palette();
  const { w, d, h } = e;
  const t = 0.0015;
  const legH = 0.15;
  const deck = 0.8;
  const rearD = 0.065;
  const r = 0.0015;
  for (const x of [0.05, w - 0.05]) for (const z of [-0.05, -d + 0.05]) gussetLeg(s, x, z, legH + 0.01, false);
  // The tank cabinet.
  blk(s, p.ss, [t, legH, -d + t], [w - t, legH + t, -t], 0);
  for (const x of [0, w - t]) blk(s, p.ss, [x, legH, -d], [x + t, deck, 0], r);
  blk(s, p.ss, [t, legH, -d], [w - t, deck, -d + t], r);
  // Front: the lower access panel, louvred, and the control fascia above it.
  blk(s, p.ss, [t, legH, -t], [w - t, 0.6, 0], r);
  const nv = s.fine(1) ? 8 : 0;
  for (let i = 0; i < nv; i++) blk(s, p.enamel, [0.12 + i * 0.05, legH + 0.03, 0.0001], [0.12 + i * 0.05 + 0.008, legH + 0.1, 0.0006], 0);
  if (s.fine(2)) for (const x of [0.03, w - 0.03]) for (const y of [legH + 0.03, 0.57]) screw(s, p.ssBright, [x, y, 0], 0.003);
  blk(s, p.ss, [t, 0.602, -0.03], [w - t, deck - 0.002, -0.028], r);
  blk(s, p.ss, [t, 0.6, -0.03], [w - t, 0.602, 0], 0);
  blk(s, p.ss, [t, deck - 0.002, -0.03], [w - t, deck, 0], 0);
  dialGauge(s, [0.12, 0.7, -0.028], 0.034, "needle_wash");
  dialGauge(s, [0.225, 0.7, -0.028], 0.034, "needle_rinse");
  blk(s, p.screen, [0.3, 0.68, -0.028], [0.42, 0.73, -0.026], 0.002);
  if (s.fine(1)) {
    s.add(p.ledGreen, cyl(s, 0.014, 0.008, "z", 20), [0.46, 0.7, -0.024]);
    s.add(p.rubber, cyl(s, 0.014, 0.008, "z", 20), [0.505, 0.7, -0.024]);
    s.add(p.ledAmber, cyl(s, 0.005, 0.003, "z", 12), [0.55, 0.72, -0.026]);
    s.add(p.bakelite, rbox(s, 0.025, 0.036, 0.012, 0.003), [0.575, 0.7, -0.022]);
  }
  // The tank top: a scrap screen of perforated sheet over the tank, and the rack track.
  const perf = perforatedSheet(0.004, 0.5);
  blk(s, p.ss, [t, deck - t, -d + rearD], [w - t, deck, -d + rearD + 0.04], 0);
  blk(s, p.ss, [t, deck - t, -0.06], [w - t, deck, -0.03], 0);
  blk(s, s.fine(1) ? perf : p.ssDull, [0.06, deck - 0.02, -d + rearD + 0.04], [w - 0.06, deck - 0.019, -0.06], 0);
  blk(s, p.ssHeat, [0.04, deck - 0.12, -d + rearD + 0.02], [w - 0.04, deck - 0.118, -0.04], 0);
  for (const z of [-0.09, -d + rearD + 0.07]) sweepX(s, p.ss, [[0.014, 0.03], [0, 0.03], [0, 0.06], [0.018, 0.06]].map(([u, v]) => [u * (z > -0.3 ? 1 : -1), v] as P2), 0.002, 0, w, deck, z);
  // Lower wash arm on its riser.
  const cx = w / 2;
  const cz = (-d + rearD - 0.03) / 2 - 0.0;
  s.add(p.ss, cyl(s, 0.022, 0.03, "y", 20), [cx, deck + 0.0, cz]);
  const lower = s.joint("wash_arm_lower", [cx, deck + 0.033, cz]);
  washArms(lower, 0.5, 1);
  // The rack, in the chamber on its track.
  dishRack(s, cx, deck + 0.063, cz, s.fine(1) ? 9 : 5);

  // The rear column: guides for the hood, springs in its head, the upper manifold forward to the upper arms.
  const top = h;
  blk(s, p.ss, [0, deck, -d], [w, top - 0.01, -d + rearD], r);
  s.add(p.ss, rbox(s, w, 0.02, rearD + 0.01, 0.008), [w / 2, top - 0.01, -d + rearD / 2]);
  for (const x of [0.018, w - 0.018]) blk(s, p.ssBright, [x - 0.01, deck, -d + rearD], [x + 0.01, top - 0.04, -d + rearD + 0.012], 0.002);
  const yu = 1.3;
  s.add(p.ss, pipe(s, [[cx, yu + 0.04, -d + rearD], [cx, yu + 0.04, cz], [cx, yu + 0.01, cz]], 0.017, 0.03, 16));
  const upper = s.joint("wash_arm_upper", [cx, yu, cz]);
  washArms(upper, 0.5, -1);
  // The rinse line up the back: copper, a pressure gauge and the vacuum breaker above the flood rim.
  if (s.fine(1)) {
    s.add(p.copper, pipe(s, [[w - 0.08, deck - 0.1, -d - 0.0], [w - 0.08, top - 0.05, -d + 0.0], [w - 0.08, top - 0.05, -d + 0.03]], 0.008, 0.03, 10).translate(0, 0, 0.012));
    s.add(p.brass, lathe(s, [[0, 0], [0.018, 0], [0.02, 0.02], [0.014, 0.05], [0.01, 0.06], [0, 0.06]], 20), [w - 0.08, top - 0.12, -d + 0.012]);
    dialGauge(s, [w - 0.16, top - 0.08, -d + rearD + 0.012], 0.022, "needle_pressure");
  }
  dataPlate(s, [w + 0.0002, 0.5, -0.2], [0, Math.PI / 2, 0], 0.08, 0.05);

  // The hood: front and side double walls and the top, lifted straight up by its handle; its back slides on the guides.
  const hy0 = deck - 0.004;
  const hy1 = 1.42;
  const pivot: V3 = [w / 2, 0.864, -d / 2];
  const hood = s.joint("hood_lift", pivot);
  const H = (q: V3): V3 => [q[0] - pivot[0], q[1] - pivot[1], q[2] - pivot[2]];
  const hz0 = -d + rearD + 0.001;
  const wt = 0.025;
  // One double-walled shell formed round the front and sides: the outer skin radiused over into the top, the bottom
  // closed by a rounded return, the inner skin up to the inner top; both tops are sheets laid into it.
  {
    const Rc = wt + 0.004;
    const nb = s.fine(2) ? 4 : 1;
    const sec: P2[] = [[wt, hy1 - wt], [wt, hy0 + 0.004], ...arc(wt - 0.004, hy0 + 0.004, 0.004, 0, -Math.PI / 2, nb).slice(1), [0.004, hy0], ...arc(0.004, hy0 + 0.004, 0.004, -Math.PI / 2, -Math.PI, nb).slice(1), [0, hy1 - 0.01], ...arc(0.01, hy1 - 0.01, 0.01, Math.PI, Math.PI / 2, s.fine(2) ? 6 : 2).slice(1), [Rc, hy1]];
    const path = filletPath([[0.004, hz0], [0.004, 0], [w - 0.004, 0], [w - 0.004, hz0]], Rc, s.fine(2) ? 8 : 2);
    hood.add(p.ss, pathSweep(sheetOutline(sec, 0.0012), path, 0), H([0, 0, 0]));
    for (const [u, y] of [
      [Rc, hy1],
      [wt, hy1 - wt],
    ] as const) {
      const g = plate(hood, planShape(offsetPath(path, u)), 0.0012);
      g.rotateX(-Math.PI / 2);
      hood.add(p.ss, g, H([0, y - 0.0006, 0]));
    }
  }
  // The lift handle along the front, low, on two standoffs; and guide shoes at the back edges.
  const hh = 0.959;
  hood.add(p.ssBright, pipe(hood, [H([0.12, hh, 0]), H([0.12, hh, 0.05]), H([w - 0.12, hh, 0.05]), H([w - 0.12, hh, 0])], 0.0115, 0.03, 16));
  for (const x of [0.12, w - 0.12]) hood.add(p.ssBright, cyl(hood, 0.018, 0.006, "z", 20), H([x, hh, 0.003]));
  for (const x of [0.018, w - 0.018]) for (const y of [hy0 + 0.05, hy1 - 0.05]) blk(hood, p.plastic, H([x - 0.016, y - 0.03, hz0 - 0.001 + 0.012]), H([x + 0.016, y + 0.03, hz0 + 0.03]), 0.003);
  if (s.fine(2)) for (const x of [0.06, w - 0.06]) for (const y of [hy0 + 0.03, hy1 - 0.03]) screw(hood, p.ssBright, H([x, y, 0]), 0.0028);

  return (b) => {
    const hj = b.joints.get("hood_lift");
    const rest: V3 = hj ? [hj.position.x, hj.position.y, hj.position.z] : pivot;
    let lift = 0;
    let aL = 0;
    let aU = 0;
    const temps = { wash_temp: 66, rinse_temp: 85 };
    const show = () => {
      turnJoint(b.joints.get("needle_wash"), "z", needleAngle(temps.wash_temp, 120));
      turnJoint(b.joints.get("needle_rinse"), "z", needleAngle(temps.rinse_temp, 120));
    };
    show();
    turnJoint(b.joints.get("needle_pressure"), "z", needleAngle(55, 120));
    return {
      open: (f) => {
        lift = Math.max(0, Math.min(1, f));
        slideJoint(hj, "y", lift * 0.44, rest);
      },
      tick: (_t, dt, load) => {
        if (lift > 0.03) return;
        const wv = (2.5 + 5 * Math.max(0, load)) * dt;
        aL += wv;
        aU -= wv * 1.1;
        turnJoint(b.joints.get("wash_arm_lower"), "y", aL);
        turnJoint(b.joints.get("wash_arm_upper"), "y", aU);
      },
      set: (name, v) => {
        if (name === "wash_temp" || name === "rinse_temp") {
          temps[name] = v;
          show();
        }
      },
    };
  };
};

let pvcM: THREE.MeshPhysicalMaterial | null = null;
/** Clear flexible PVC strip curtain: faintly blue-grey, glossy, mostly see-through, its edges catching the light. */
function pvcCurtain(): THREE.MeshPhysicalMaterial {
  if (!pvcM) {
    pvcM = new THREE.MeshPhysicalMaterial({ color: 0xc9d6da, roughness: 0.22, metalness: 0, transparent: true, opacity: 0.38, clearcoat: 0.6, side: THREE.DoubleSide, depthWrite: false });
    pvcM.userData.uv = "own";
  }
  return pvcM;
}

/** A sheet-metal transition, open top and bottom: a `w0` × `d0` base at y = 0 to a `w1` × `d1` top at `h`, centred. */
function frustum(w0: number, d0: number, w1: number, d1: number, h: number): THREE.BufferGeometry {
  const b = [[-w0 / 2, 0, d0 / 2], [w0 / 2, 0, d0 / 2], [w0 / 2, 0, -d0 / 2], [-w0 / 2, 0, -d0 / 2]];
  const u = [[-w1 / 2, h, d1 / 2], [w1 / 2, h, d1 / 2], [w1 / 2, h, -d1 / 2], [-w1 / 2, h, -d1 / 2]];
  const pos: number[] = [];
  for (let i = 0; i < 4; i++) {
    const j = (i + 1) % 4;
    pos.push(...b[i]!, ...b[j]!, ...u[j]!, ...b[i]!, ...u[j]!, ...u[i]!);
  }
  const g = new THREE.BufferGeometry();
  g.setAttribute("position", new THREE.Float32BufferAttribute(pos, 3));
  g.computeVertexNormals();
  return g;
}

/**
 * A single-tank rack conveyor: racks ride the track left to right through a tunnel, pushed by a reciprocating pawl bar;
 * vinyl curtains hold the spray in at the entry, between zones and at the exit; vent cowls at each end draw the vapour
 * off. The wash pump hangs under the tank, its motor beside it. An inspection door on the front, controls on top.
 */
const conveyorDishwasher: Builder = (s, e) => {
  const p = palette();
  const { w, d, h } = e;
  const t = 0.0015;
  const r = 0.0015;
  const tank0 = 0.36;
  const deck = 0.8;
  const zf = -0.055;
  const zb = -d + 0.03;
  const roof = 1.42;
  for (const x of [0.05, w / 2, w - 0.05]) for (const z of [-0.08, -d + 0.06]) gussetLeg(s, x, z, tank0 + 0.01, false);
  // The tank: a welded box under the tunnel, its front panel with the access cover.
  blk(s, p.ss, [t, tank0, zb], [w - t, tank0 + t, zf - 0.0015], 0);
  blk(s, p.ss, [0, tank0, zf - 0.0015], [w, deck, zf], r);
  blk(s, p.ss, [0, tank0, zb - t], [w, deck, zb], r);
  for (const x of [0, w - t]) blk(s, p.ss, [x, tank0, zb], [x + t, deck, zf - 0.0015], r);
  blk(s, p.ss, [0.25, tank0 + 0.08, zf], [0.85, deck - 0.06, zf + 0.002], 0.003);
  if (s.fine(2)) for (const x of [0.27, 0.83]) for (const y of [tank0 + 0.1, deck - 0.08]) screw(s, p.ssBright, [x, y, zf + 0.002], 0.0035);
  // The wash pump and its motor under the tank: a cast volute, the motor's finned frame, its conduit box.
  const py = 0.25;
  const pz = -0.33;
  s.add(p.ssDull, lathe(s, [[0, -0.04], [0.06, -0.04], [0.075, -0.025], [0.075, 0.025], [0.06, 0.04], [0, 0.04]], 32).rotateZ(Math.PI / 2), [0.48, py, pz]);
  s.add(p.ssDull, cyl(s, 0.028, tank0 - py, "y", 18), [0.48, (tank0 + py) / 2 + 0.03, pz]);
  s.add(p.ssDull, pipe(s, [[0.48, py + 0.05, pz - 0.07], [0.48, py + 0.05, zb + 0.03], [0.48, deck + 0.1, zb + 0.03]], 0.024, 0.05, 14));
  s.add(p.charcoal, cyl(s, 0.07, 0.2, "x", 32), [0.64, py, pz]);
  if (s.fine(1)) for (let i = 0; i < 9; i++) s.add(p.charcoal, cyl(s, 0.075, 0.006, "x", 32), [0.565 + i * 0.019, py, pz]);
  s.add(p.charcoal, lathe(s, [[0, 0], [0.072, 0], [0.068, 0.02], [0.03, 0.04], [0, 0.04]], 32).rotateZ(-Math.PI / 2), [0.74, py, pz]);
  blk(s, p.charcoal, [0.6, py + 0.07, pz - 0.04], [0.68, py + 0.12, pz + 0.04], 0.004);
  s.add(p.ssDull, cyl(s, 0.009, deck - py - 0.12, "y", 10), [0.64, (deck + py + 0.12) / 2, pz]);
  dataPlate(s, [0.66, py + 0.03, pz + 0.073], [0, 0, 0], 0.05, 0.03);

  // The tunnel: front with the door opening, back, roof, and end walls with the rack openings.
  const o: [number, number, number, number] = [0.864, 1.31, -d + 0.13, -0.14];
  const [dx0, dx1, dy0, dy1] = [0.28, 0.85, 0.864, 1.364];
  blk(s, p.ss, [0, deck, zf - t], [dx0, roof, zf], r);
  blk(s, p.ss, [dx1, deck, zf - t], [w, roof, zf], r);
  blk(s, p.ss, [dx0, deck, zf - t], [dx1, dy0, zf], r);
  blk(s, p.ss, [dx0, dy1, zf - t], [dx1, roof, zf], r);
  blk(s, p.ss, [0, deck, zb - t], [w, roof, zb], r);
  blk(s, p.ss, [0, roof - t, zb], [w, roof, zf - t], r);
  for (const x of [0, w - t]) {
    blk(s, p.ss, [x, deck, zb], [x + t, o[0], zf - t], 0);
    blk(s, p.ss, [x, o[1], zb], [x + t, roof - t, zf - t], 0);
    blk(s, p.ss, [x, o[0], zb], [x + t, o[1], o[2]], 0);
    blk(s, p.ss, [x, o[0], o[3]], [x + t, o[1], zf - t], 0);
    // A rolled lip round the opening.
    if (s.fine(1)) {
      s.add(p.ss, cyl(s, 0.005, o[1] - o[0], "y", 10), [x + 0.0007, (o[0] + o[1]) / 2, o[2]]);
      s.add(p.ss, cyl(s, 0.005, o[1] - o[0], "y", 10), [x + 0.0007, (o[0] + o[1]) / 2, o[3]]);
      s.add(p.ss, cyl(s, 0.005, o[3] - o[2], "z", 10), [x + 0.0007, o[1], (o[2] + o[3]) / 2]);
    }
  }
  // The track through the tunnel and the pawl bar between its rails.
  for (const z of [o[2] + 0.015, o[3] - 0.015]) sweepX(s, p.ss, [[0.014, 0.0], [0, 0.0], [0, 0.03], [0.018, 0.03]].map(([u, v]) => [u * (z < -0.4 ? 1 : -1), v] as P2), 0.002, 0.002, w - 0.002, deck + 0.034, z);
  const zp = (o[2] + o[3]) / 2;
  const pawl = s.joint("conveyor_pawl", [w / 2, deck + 0.05, zp]);
  pawl.add(p.ssBright, rbox(pawl, w - 0.12, 0.012, 0.03, 0.002), [0, 0, 0]);
  for (let i = 0; i < 8; i++) pawl.add(p.ssBright, taperBlade(pawl, 0.04, 0.018, 0.012, 0.012), [-w / 2 + 0.1 + i * 0.13, 0.006, 0], [0, 0, 0.5]);
  blk(s, p.ssHeat, [0.01, deck - t, zb], [w - 0.01, deck, zf - 0.01], 0);
  blk(s, s.fine(1) ? perforatedSheet(0.004, 0.5) : p.ssDull, [0.05, deck + 0.006, o[2]], [w - 0.05, deck + 0.007, o[3]], 0);
  // Spray manifolds over and under the rack line, and the curtains.
  for (const x of s.fine(1) ? [0.38, 0.62] : []) {
    s.add(p.ss, cyl(s, 0.016, o[3] - o[2], "z", 16), [x, 1.26, zp]);
    s.add(p.ss, cyl(s, 0.016, o[3] - o[2], "z", 16), [x, deck + 0.02, zp]);
    if (s.fine(2)) for (let i = 0; i < 7; i++) for (const y of [1.245, deck + 0.035]) s.add(p.ssBright, cyl(s, 0.004, 0.008, "y", 8), [x, y, o[2] + 0.04 + i * 0.07]);
  }
  const curtains = [0.02, 0.27, 0.86, w - 0.02];
  if (s.fine(1)) curtains.forEach((x, i) => {
    s.add(p.ssBright, cyl(s, 0.005, o[3] - o[2] + 0.04, "z", 10), [x, o[1] - 0.004, zp]);
    const c = s.joint(`curtain_${i}`, [x, o[1] - 0.01, zp]);
    // Clear PVC strips, each lapping the next, hung on hooks from the rod; they stop short of the racks' path.
    const n = 6;
    const sw = (o[3] - o[2]) / n;
    const sl = o[1] - o[0] - 0.06;
    const strip = rbox(c, 0.002, sl, sw + 0.012, 0.0009);
    for (let k = 0; k < n; k++) {
      const zk = -(o[3] - o[2]) / 2 + sw * (k + 0.5);
      c.add(pvcCurtain(), strip.clone(), [(k % 2 ? 1 : -1) * 0.0016, -sl / 2 - 0.004, zk]);
      if (c.fine(2)) c.add(p.ssBright, torus(c, 0.007, 0.0012, "z", Math.PI * 2, 12, 4), [0, 0.006, zk]);
    }
  });

  // Vent cowls over each end, their dampers' handles, and the control box on the roof.
  for (const xc of [0.16, w - 0.16]) {
    s.add(p.ss, frustum(0.38, 0.22, 0.22, 0.12, 0.09), [xc, roof, -d * 0.45]);
    ductCollar(s, xc - 0.11, xc + 0.11, -d * 0.45 - 0.06, -d * 0.45 + 0.06, roof + 0.09, h);
    s.add(p.ssDull, cyl(s, 0.004, 0.24, "x", 8), [xc, roof + 0.2, -d * 0.45]);
    s.add(p.bakelite, rbox(s, 0.012, 0.06, 0.012, 0.004), [xc + 0.13, roof + 0.22, -d * 0.45]);
  }
  const [cx0, cx1, cy0, cy1, cz0, cz1] = [0.38, 0.74, roof, roof + 0.22, -0.32, -0.065];
  blk(s, p.ss, [cx0, cy0, cz0], [cx1, cy1, cz1 - 0.07], 0.004);
  s.at([0, cy0, cz1 - 0.07], [-0.5, 0, 0], () => {
    blk(s, p.ss, [cx0, 0, -0.001], [cx1, 0.2, 0.075], 0.004);
  });
  s.at([(cx0 + cx1) / 2, cy0 + 0.1, cz1 - 0.035], [-0.5, 0, 0], () => {
    s.add(p.screen, rbox(s, 0.12, 0.05, 0.002, 0.002), [-0.05, 0.035, 0.0385]);
    dialGauge(s, [0.1, 0.035, 0.0375], 0.026, "needle_wash");
    dialGauge(s, [0.1, -0.035, 0.0375], 0.026, "needle_rinse");
    for (let i = 0; i < 3; i++) s.add(i ? p.rubber : p.ledGreen, cyl(s, 0.011, 0.006, "z", 16), [-0.09 + i * 0.04, -0.03, 0.04]);
    s.add(p.ledRed, lathe(s, [[0, 0], [0.018, 0], [0.019, 0.008], [0.014, 0.014], [0, 0.015]], 20).rotateX(Math.PI / 2), [0.03, -0.03, 0.038]);
  });
  // The final rinse: brass solenoid valve and vacuum breaker on the roof at the exit end, the line into the tunnel.
  if (s.fine(1)) {
    s.add(p.copper, pipe(s, [[w - 0.06, roof + 0.003, zb + 0.06], [w - 0.06, roof + 0.12, zb + 0.06], [w - 0.32, roof + 0.12, zb + 0.06]], 0.0095, 0.03, 10));
    s.add(p.brass, rbox(s, 0.05, 0.04, 0.04, 0.005), [w - 0.2, roof + 0.12, zb + 0.06]);
    s.add(p.charcoal, cyl(s, 0.016, 0.04, "y", 16), [w - 0.2, roof + 0.16, zb + 0.06]);
    s.add(p.brass, lathe(s, [[0, 0], [0.016, 0], [0.018, 0.02], [0.012, 0.045], [0, 0.05]], 18), [w - 0.06, roof + 0.13, zb + 0.06]);
  }

  // The inspection door on its hinge at the left of the opening, a pull handle at its right.
  const hinge: V3 = [dx0, 1.124, zf];
  const door = s.joint("inspect_door", hinge);
  const dw = dx1 - dx0;
  blk(door, p.ss, [0.003, dy0 - hinge[1] + 0.003, -0.0015], [dw - 0.003, dy1 - hinge[1] - 0.003, 0.001], r);
  blk(door, p.ss, [0.03, dy0 - hinge[1] + 0.03, -0.03], [dw - 0.03, dy1 - hinge[1] - 0.03, -0.0015], r);
  blk(door, p.gasket, [0.02, dy0 - hinge[1] + 0.02, -0.008], [dw - 0.02, dy1 - hinge[1] - 0.02, -0.0025], 0.003);
  pullHandle(door, p.ssBright, [0.76 - dx0, 0, 0.001], 0.2, 0.045, 0.01, "y");
  for (const y of [-0.18, 0.18]) door.add(p.ssBright, cyl(door, 0.008, 0.06, "y", 14), [-0.004, y, 0.004]);
  if (s.fine(2)) for (const x of [0.04, dw - 0.04]) for (const y of [-0.21, 0.21]) screw(door, p.ssBright, [x, y, 0.001], 0.003);
  for (const y of [-0.18, 0.18]) s.add(p.ssBright, cyl(s, 0.0085, 0.03, "y", 14), [dx0 - 0.004, 1.124 + y + 0.045, zf + 0.004]);

  return (b) => {
    turnJoint(b.joints.get("needle_wash"), "z", needleAngle(66, 120));
    turnJoint(b.joints.get("needle_rinse"), "z", needleAngle(85, 120));
    const pj = b.joints.get("conveyor_pawl");
    const rest: V3 = pj ? [pj.position.x, pj.position.y, pj.position.z] : [w / 2, deck + 0.05, zp];
    let phase = 0;
    return {
      open: (f) => turnJoint(b.joints.get("inspect_door"), "y", -Math.max(0, Math.min(1, f)) * (Math.PI / 2)),
      tick: (tt, dt, load) => {
        const l = Math.max(0, Math.min(1, load));
        phase += dt * 0.5 * l;
        slideJoint(pj, "x", 0.06 * Math.sin(phase * Math.PI * 2), rest);
        curtains.forEach((_, i) => turnJoint(b.joints.get(`curtain_${i}`), "z", l * 0.05 * Math.sin(tt * 1.7 + i * 1.3)));
      },
    };
  };
};


/* ------------------------------------------------------------------------------------------------------------ */
/* The insulated carrier                                                                                          */
/* ------------------------------------------------------------------------------------------------------------ */

let peMat: THREE.MeshStandardMaterial | null = null;
/** Rotomoulded polyethylene in deep charcoal, as catering carriers are moulded. One material for every carrier. */
export function carrierPE(): THREE.MeshStandardMaterial {
  if (!peMat) {
    // Rotomoulded polyethylene: the mould's fine orange-peel stipple in the bump and a little in the sheen.
    const c = document.createElement("canvas");
    c.width = c.height = 128;
    const x = c.getContext("2d")!;
    x.fillStyle = "#808080";
    x.fillRect(0, 0, 128, 128);
    let k = 12345;
    const rnd = () => ((k = (k * 1103515245 + 12345) & 0x7fffffff) / 0x7fffffff);
    for (let i = 0; i < 2600; i++) {
      const [px, py, r] = [rnd() * 128, rnd() * 128, 0.8 + rnd() * 1.8];
      const v = Math.round(110 + rnd() * 50);
      x.fillStyle = `rgba(${v},${v},${v},0.55)`;
      for (const [ox, oy] of [[0, 0], [128, 0], [0, 128], [-128, 0], [0, -128]]) {
        x.beginPath();
        x.arc(px + ox!, py + oy!, r, 0, Math.PI * 2);
        x.fill();
      }
    }
    const t = new THREE.CanvasTexture(c);
    t.wrapS = t.wrapT = THREE.RepeatWrapping;
    peMat = new THREE.MeshStandardMaterial({ color: 0x2d2f33, roughness: 0.64, metalness: 0, bumpMap: t, bumpScale: 0.3 });
    peMat.userData.uv = "box";
    peMat.userData.uvScale = 1 / 0.06;
  }
  return peMat;
}

/** A rounded box `w`×`h`×`d` centred on the origin with its +z face and the rounding round it cut away. */
function openFrontBox(s: Shop, w: number, h: number, d: number, r: number): THREE.BufferGeometry {
  if (s.detail === 0) return cutAbove(new THREE.BoxGeometry(w, h, d), "z", d / 2 - 1e-4);
  return cutAbove(new RoundedBoxGeometry(w, h, d, s.detail === 2 ? 3 : 1, r), "z", d / 2 - r + 1e-5);
}

/** A GN 1/1 pan, 65 mm deep, with its flat lid, in its own frame: rim top at y = 0, long side along z. */
function gnPan(s: Shop): Part[] {
  const p = palette();
  const parts: Part[] = [];
  const [W, L] = [0.325, 0.53];
  const rim = roundRect(W, L, 0.022);
  rim.holes.push(rectPath(new THREE.Path(), W - 0.026, L - 0.026, 0.016, 0, 0, true) as THREE.Path);
  const rg = plate(s, rim, 0.0025, s.fine(2) ? 0.0008 : 0, 12);
  rg.rotateX(-Math.PI / 2);
  parts.push({ m: p.ss, g: rg, at: [0, -0.0025, 0] });
  // The body, drawn in one piece with coved corners, tapering a little to its floor.
  const tmp = new Shop(s.detail, "pan");
  basin(tmp, p.ss, 0, 0, W - 0.03, L - 0.03, 0.063, -0.0025, 0.022, 0.0008);
  for (const g of shopGeoms(tmp)) parts.push({ m: p.ss, g });
  // The lid: a pressed sheet with a raised centre field, notched handle at its middle.
  parts.push({ m: p.ss, g: rbox(s, W + 0.004, 0.0018, L + 0.004, 0.0008), at: [0, 0.0009, 0] });
  if (s.fine(2)) {
    parts.push({ m: p.ss, g: rbox(s, W - 0.05, 0.004, L - 0.06, 0.003), at: [0, 0.003, 0] });
    parts.push({ m: p.ssBright, g: pipe(s, [[0, 0.005, -0.05], [0, 0.022, -0.035], [0, 0.022, 0.035], [0, 0.005, 0.05]], 0.0035, 0.01, 8) });
  }
  return parts;
}
/** The geometries a scratch shop collected (its parts are all one material, the caller knows which). */
function shopGeoms(tmp: Shop): THREE.BufferGeometry[] {
  const b = tmp.build(false);
  const out: THREE.BufferGeometry[] = [];
  b.group.traverse((o) => {
    const m = o as THREE.Mesh;
    if (m.isMesh) out.push(m.geometry);
  });
  return out;
}

/**
 * A front-loading insulated food-pan carrier: a rotomoulded double-wall polyethylene shell foamed between its skins,
 * moulded feet and side hand-holds, rails inside for six GN 1/1 pans; its door hinged at the left on stainless pins,
 * sealed by a gasket, held by two over-centre latches, a vent cap in its face that turns to let steam out.
 */
const insulatedCarrier: Builder = (s, e) => {
  const p = palette();
  const pe = carrierPE();
  const { w, d, h } = e;
  const [x0, x1] = [0.016, w - 0.014];
  const [y0, y1] = [0.02, h];
  const zF = -0.058;
  const R = 0.034;
  const bw = x1 - x0;
  const bh = y1 - y0;
  const bd = zF + d;
  const bx = (x0 + x1) / 2;
  const by = (y0 + y1) / 2;
  // The outer skin, open at the front where the moulded front frame closes the double wall.
  s.add(pe, openFrontBox(s, bw, bh, bd, R), [bx, by, -d + bd / 2]);
  const [cw, ch] = [0.36, 0.525];
  const [cx, cy] = [bx, 0.06 + ch / 2];
  const rr = 0.026;
  const bev = s.fine(1) ? 0.007 : 0;
  const frame = roundRect(bw - 2 * bev, bh - 2 * bev, R - bev, 0, 0);
  frame.holes.push(rectPath(new THREE.Path(), cw + 2 * bev + 0.004, ch + 2 * bev + 0.004, rr + bev + 0.002, cx - bx, cy - by, true) as THREE.Path);
  s.add(pe, plate(s, frame, R, bev, 16), [bx, by, zF - R]);
  // The inner skin: the cavity, its corners radiused, run forward through the frame to its face.
  const back = -d + 0.034;
  const Lc = zF - R - back;
  s.add(pe, inward(openFrontBox(s, cw, ch, Lc + R + rr, rr)), [cx, cy, back + (Lc + R + rr) / 2]);
  // The moulded nosing where liner and frame meet: a rounded lip round the opening, standing a little proud of the
  // frame for the door's gasket to seal on.
  if (bev) {
    const nb = s.fine(2) ? 5 : 2;
    const sec: P2[] = [[0.0015, zF - 0.02], [0.0015, zF - 0.0025], ...arc(-0.0025, zF - 0.0025, 0.004, 0, Math.PI / 2, nb).slice(1), [-(bev + 0.006), zF + 0.0015], ...arc(-(bev + 0.006), zF, 0.0015, Math.PI / 2, Math.PI, 2).slice(1), [-(bev + 0.0075), zF - 0.003], [-0.0005, zF - 0.003], [-0.0005, zF - 0.02]];
    const [xl, xr, yb, yt] = [cx - cw / 2, cx + cw / 2, cy - ch / 2, cy + ch / 2];
    const g = pathSweep(sec, filletPath([[xl, -yb], [xr, -yb], [xr, -yt], [xl, -yt]], rr, s.fine(2) ? 8 : 2, true), 0, true);
    g.rotateX(Math.PI / 2);
    s.add(pe, g);
  }
  // Moulded rails for six pans down both side walls.
  const pitch = 0.085;
  const rims: number[] = [];
  for (let i = 0; i < 6; i++) rims.push(0.135 + i * pitch);
  if (s.fine(1))
    for (const y of rims)
      for (const sx of [-1, 1]) s.add(pe, rbox(s, 0.022, 0.012, Lc - 0.03, 0.005), [cx + sx * (cw / 2 - 0.009), y - 0.006, back + (Lc - 0.03) / 2 + 0.02]);
  // Feet, stacking lugs and ribs on top, and the side hand-holds: a raised rim round a deep pocket with its grip.
  for (const sx of [-1, 1]) for (const sz of [-1, 1]) s.add(pe, rbox(s, 0.09, 0.03, 0.11, 0.01), [bx + sx * (bw / 2 - 0.06), y0 + 0.005, -d / 2 - 0.03 + sz * (bd / 2 - 0.07)]);
  if (s.fine(1)) {
    for (const sx of [-1, 1]) s.add(pe, rbox(s, 0.04, 0.012, bd - 0.14, 0.005), [bx + sx * 0.12, y1 - 0.002, -d + bd / 2]);
    for (const sx of [-1, 1]) for (const sz of [-1, 1]) s.add(pe, cyl(s, 0.018, 0.01, "y", 20), [bx + sx * (bw / 2 - 0.05), y1 + 0.0, -d + bd / 2 + sz * (bd / 2 - 0.05)]);
  }
  const hy = 0.495;
  const hz = -d + bd / 2 + 0.02;
  for (const [xs, sx] of [
    [x0, -1],
    [x1, 1],
  ] as const) {
    const ring = roundRect(0.2, 0.075, 0.03);
    ring.holes.push(rectPath(new THREE.Path(), 0.16, 0.045, 0.018, 0, 0, true) as THREE.Path);
    const g = plate(s, ring, 0.012, s.fine(1) ? 0.004 : 0, 12);
    g.rotateY(Math.PI / 2);
    s.add(pe, g, [xs + sx * 0.004 - 0.006, hy, hz]);
    s.add(p.rubber, rbox(s, 0.003, 0.045, 0.16, 0.002), [xs + sx * 0.0, hy, hz]);
    s.add(pe, cyl(s, 0.011, 0.165, "z", 16), [xs + sx * 0.0, hy + 0.012, hz]);
  }
  // Latch keepers moulded on the right side, hinge leaves on the left.
  const lY = [0.2, 0.46];
  for (const y of lY) s.add(pe, rbox(s, 0.008, 0.04, 0.02, 0.003), [x1 + 0.003, y, zF - 0.022]);
  const ax: V3 = [0.01, h / 2, -0.004];
  const hingeY = [0.09, 0.31, 0.53];
  for (const y of hingeY) {
    s.add(pe, rbox(s, 0.01, 0.04, 0.06, 0.004), [x0 + 0.002, y + 0.038, zF - 0.012]);
    s.add(pe, cyl(s, 0.011, 0.04, "y", 16), [ax[0], y + 0.038, ax[2]]);
    s.add(p.ssBright, cyl(s, 0.0035, 0.082, "y", 10), [ax[0], y + 0.0205, ax[2]]);
    if (s.fine(1)) for (const yy of [y - 0.0205, y + 0.0615]) s.add(p.ssBright, cyl(s, 0.0055, 0.003, "y", 12), [ax[0], yy, ax[2]]);
  }
  dataPlate(s, [bx, 0.08, -d + 0.0002], [0, Math.PI, 0], 0.08, 0.05);

  // Pans, each its own joint so that a load shows.
  const pan = gnPan(s);
  rims.forEach((y, i) => stamp(s.joint(`pan_${i}`, [cx, y, back + 0.004 + 0.265]), pan, [0, 0, 0]));

  // The door: a foamed double-wall panel, its gasket, moulded face, vent cap and latches.
  const door = s.joint("carrier_door", ax);
  const Dq = (q: V3): V3 => [q[0] - ax[0], q[1] - ax[1], q[2] - ax[2]];
  const [dx0, dx1, dy0, dy1, dz0, dz1] = [0.02, w - 0.016, y0 + 0.008, y1 - 0.008, -0.052, -0.006];
  const dc = Dq([(dx0 + dx1) / 2, (dy0 + dy1) / 2, (dz0 + dz1) / 2]);
  door.add(pe, s.detail ? new RoundedBoxGeometry(dx1 - dx0, dy1 - dy0, dz1 - dz0, s.detail === 2 ? 3 : 1, 0.02) : new THREE.BoxGeometry(dx1 - dx0, dy1 - dy0, dz1 - dz0), dc);
  const gk = roundRect(cw + 0.05, ch + 0.05, rr + 0.02);
  gk.holes.push(rectPath(new THREE.Path(), cw + 0.022, ch + 0.022, rr + 0.008, 0, 0, true) as THREE.Path);
  door.add(p.rubber, plate(door, gk, 0.0055, s.fine(1) ? 0.0015 : 0, 12), Dq([cx, cy, dz0 - 0.0055]));
  if (s.fine(1)) {
    const field = roundRect(dx1 - dx0 - 0.05, dy1 - dy0 - 0.05, 0.03);
    field.holes.push(rectPath(new THREE.Path(), dx1 - dx0 - 0.08, dy1 - dy0 - 0.08, 0.02, 0, 0, true) as THREE.Path);
    door.add(pe, plate(door, field, 0.004, 0.0015, 12), Dq([(dx0 + dx1) / 2, (dy0 + dy1) / 2, dz1 - 0.0015]));
    // The label holder: a satin frame round a white card.
    door.add(p.ssBright, rbox(door, 0.13, 0.08, 0.002, 0.002), Dq([w / 2, 0.36, dz1 + 0.0008]));
    door.add(p.hdpe, rbox(door, 0.12, 0.068, 0.0015, 0.0008), Dq([w / 2, 0.36, dz1 + 0.0021]));
  }
  // The vent seat: a boss with four slots through the door; the cap turns over it.
  const vc: V3 = [w / 2, 0.52, dz1];
  door.add(pe, cyl(door, 0.04, 0.006, "z", 32), Dq([vc[0], vc[1], vc[2] + 0.002]));
  const sector = (r0: number, r1: number, a0: number, a1: number) => {
    const sh = new THREE.Shape();
    sh.absarc(0, 0, r1, a0, a1, false);
    sh.absarc(0, 0, r0, a1, a0, true);
    return sh;
  };
  for (let i = 0; i < 4; i++) {
    const a = (i * Math.PI) / 2;
    door.add(p.rubber, plate(door, sector(0.012, 0.03, a - 0.32, a + 0.32), 0.0008, 0, 8), Dq([vc[0], vc[1], vc[2] + 0.0052]));
  }
  const cap = door.joint("vent_cap", Dq([vc[0], vc[1], vc[2] + 0.006]));
  const disc = new THREE.Shape();
  disc.absarc(0, 0, 0.036, 0, Math.PI * 2, false);
  for (let i = 0; i < 4; i++) {
    const a = (i * Math.PI) / 2 + Math.PI / 4;
    const hole = new THREE.Path();
    hole.absarc(0, 0, 0.0305, a + 0.3, a - 0.3, true);
    hole.absarc(0, 0, 0.0115, a - 0.3, a + 0.3, false);
    disc.holes.push(hole);
  }
  cap.add(pe, plate(cap, disc, 0.006, s.fine(1) ? 0.0015 : 0, 32), [0, 0, 0]);
  cap.add(pe, rbox(cap, 0.012, 0.064, 0.012, 0.004), [0, 0, 0.009]);
  cap.add(p.ssBright, rbox(cap, 0.003, 0.01, 0.002, 0.001), [0, 0.04, 0.003]);
  // Index marks on the boss: shut where the pointer rests, open a quarter-turn's half clockwise.
  for (const a of [Math.PI / 2, Math.PI / 4]) door.add(a === Math.PI / 2 ? p.ssDull : p.hdpe, cyl(door, 0.0022, 0.001, "z", 10), Dq([vc[0] + Math.cos(a) * 0.046, vc[1] + Math.sin(a) * 0.046, vc[2] + 0.0003]));
  // Latches: over-centre levers on the door's face, their hooks round the right edge onto the keepers.
  lY.forEach((y, i) => {
    const lj = door.joint(`latch_${i}`, Dq([x1 - 0.022, y, dz1 + 0.004]));
    lj.add(p.plastic, rbox(lj, 0.04, 0.034, 0.012, 0.004), [0.012, 0, 0.002]);
    lj.add(p.plastic, rbox(lj, 0.008, 0.03, 0.03, 0.003), [0.034, 0, -0.012]);
    lj.add(p.plastic, rbox(lj, 0.009, 0.026, 0.012, 0.003), [0.031, 0, -0.03]);
    lj.add(p.ssBright, cyl(lj, 0.003, 0.04, "y", 8), [-0.004, 0, 0]);
    door.add(p.plastic, rbox(door, 0.016, 0.044, 0.008, 0.003), Dq([x1 - 0.026, y, dz1 + 0.002]));
  });
  for (const y of hingeY) {
    door.add(pe, cyl(door, 0.011, 0.034, "y", 16), Dq([ax[0], y - 0.0, ax[2]]));
    door.add(pe, rbox(door, 0.02, 0.034, 0.016, 0.004), Dq([0.022, y, ax[2] - 0.004]));
  }

  return (b) => {
    let pans = 3;
    const showPans = () => {
      for (let i = 0; i < 6; i++) {
        const j = b.joints.get(`pan_${i}`);
        if (j) j.visible = i < pans;
      }
    };
    showPans();
    const vent = (v: number) => turnJoint(b.joints.get("vent_cap"), "z", -Math.max(0, Math.min(1, v)) * (Math.PI / 4));
    vent(0);
    return {
      open: (f) => {
        const k = Math.max(0, Math.min(1, f));
        for (let i = 0; i < 2; i++) turnJoint(b.joints.get(`latch_${i}`), "z", Math.min(1, k * 8) * 1.3);
        turnJoint(b.joints.get("carrier_door"), "y", -Math.max(0, (k - 0.1) / 0.9) * Math.PI);
      },
      set: (name, v) => {
        if (name === "vent") vent(v);
        else if (name === "pans") {
          pans = Math.max(0, Math.min(6, Math.round(v)));
          showPans();
        }
      },
    };
  };
};


/* ------------------------------------------------------------------------------------------------------------ */
/* Prep machines and smallwares                                                                                   */
/* ------------------------------------------------------------------------------------------------------------ */

/**
 * A commercial batch-bowl food processor: a stainless motor base with its paddle switches, the drive shaft through a
 * clear polycarbonate bowl's centre column, the S-blade on its hub (two scimitar blades at two heights, so that
 * everything in the bowl meets an edge), and a lid with its feed tube and pusher, the lid's tab working the safety rod.
 */
const foodProcessor: Builder = (s, e) => {
  const p = palette();
  const { w, d } = e;
  const cx = w / 2;
  const cz = -0.155;
  const r = 0.02;
  if (s.detail === 0) {
    // Far off: the base, the bowl and lid, the feed tube; the blade and pusher still on their joints.
    blk(s, p.ss, [0.012, 0, -d + 0.008], [w - 0.012, 0.172, -0.02], 0);
    const bowl = new THREE.CylinderGeometry(0.098, 0.09, 0.133, 10, 1, true);
    s.add(p.glass, bowl, [cx, 0.2385, cz]);
    s.add(p.glass, new THREE.CylinderGeometry(0.1, 0.1, 0.02, 10), [cx, 0.315, cz]);
    s.add(p.glass, new THREE.CylinderGeometry(0.033, 0.033, 0.07, 8, 1, true), [cx - 0.035, 0.358, cz + 0.03]);
    s.joint("s_blade_spin", [cx, 0.2, cz]).add(p.ssBright, new THREE.BoxGeometry(0.15, 0.004, 0.02), [0, 0.01, 0]);
    s.joint("pusher", [cx - 0.035, 0.393, cz + 0.03]).add(p.hdpe, new THREE.CylinderGeometry(0.03, 0.03, 0.08, 8), [0, -0.03, 0]);
    return (b) => {
      const pj = b.joints.get("pusher");
      const rest: V3 = pj ? [pj.position.x, pj.position.y, pj.position.z] : [cx - 0.035, 0.393, cz + 0.03];
      let a = 0;
      return {
        tick: (_t, dt, load) => {
          a += dt * 45 * Math.max(0, Math.min(1, load));
          turnJoint(b.joints.get("s_blade_spin"), "y", a);
        },
        // Far off the lid stays put; the pusher draws up its tube, never out of it.
        open: (f) => slideJoint(pj, "y", Math.max(0, Math.min(1, f)) * 0.025, rest),
      };
    };
  }
  // The base.
  for (const sx of [-1, 1]) for (const sz of [-1, 1]) s.add(p.rubber, cyl(s, 0.012, 0.012, "y", 16), [cx + sx * 0.085, 0.006, -d / 2 + sz * 0.11]);
  blk(s, p.ss, [0.012, 0.012, -d + 0.008], [w - 0.012, 0.165, -0.02], r);
  blk(s, p.plastic, [0.016, 0.163, -d + 0.012], [w - 0.016, 0.172, -0.024], 0.006);
  // The front panel and its two paddles: on, and off / pulse below it.
  blk(s, p.plastic, [0.04, 0.035, -0.0215], [w - 0.04, 0.14, -0.019], 0.006);
  s.add(p.bakelite, rbox(s, 0.11, 0.032, 0.012, 0.005), [cx, 0.11, -0.014]);
  s.add(p.bakelite, rbox(s, 0.11, 0.032, 0.014, 0.005), [cx, 0.068, -0.013], [0.12, 0, 0]);
  s.add(p.ledGreen, rbox(s, 0.03, 0.006, 0.0015, 0.0006), [cx, 0.11, -0.0076]);
  s.add(p.ledRed, rbox(s, 0.03, 0.006, 0.0015, 0.0006), [cx, 0.068, -0.0058]);
  if (s.fine(2)) for (const x of [0.0118, w - 0.0128]) for (let i = 0; i < 6; i++) blk(s, p.enamel, [x, 0.04, -0.24 + i * 0.022], [x + 0.001, 0.12, -0.232 + i * 0.022], 0);
  dataPlate(s, [cx, 0.09, -d + 0.0078], [0, Math.PI, 0], 0.07, 0.04);
  // The drive: the shaft up through the bowl's column; the safety rod at the back right.
  const yb = 0.172;
  s.add(p.ssBright, cyl(s, 0.006, 0.1, "y", 12), [cx, yb + 0.05, cz]);
  s.add(p.plastic, cyl(s, 0.03, 0.01, "y", 24), [cx, yb + 0.005, cz]);
  s.add(p.ssBright, cyl(s, 0.004, 0.02, "y", 10), [cx + 0.07, yb + 0.01, cz - 0.07]);
  // The bowl: a thick-walled clear lathe with its centre column, a handle front-right, bayonet lugs at its foot.
  const yt = 0.305;
  s.add(p.glass, lathe(s, [[0.012, yb + 0.012], [0.086, yb + 0.012], [0.09, yb + 0.02], [0.095, yt - 0.004], [0.098, yt], [0.092, yt], [0.089, yb + 0.022], [0.084, yb + 0.016], [0.014, yb + 0.016], [0.012, 0.29], [0.009, 0.29], [0.009, yb + 0.012]], 48), [cx, 0, cz]);
  s.add(p.glass, pipe(s, [[cx + 0.067, yt - 0.02, cz + 0.067], [cx + 0.098, yt - 0.03, cz + 0.098], [cx + 0.098, yb + 0.03, cz + 0.098], [cx + 0.067, yb + 0.03, cz + 0.067]], 0.009, 0.02, 10));
  if (s.fine(2)) for (let i = 0; i < 3; i++) {
    const a = (i / 3) * Math.PI * 2 + 0.3;
    s.add(p.glass, rbox(s, 0.02, 0.008, 0.012, 0.003), [cx + Math.cos(a) * 0.09, yb + 0.016, cz + Math.sin(a) * 0.09], [0, -a, 0]);
  }
  // The S-blade.
  const blade = s.joint("s_blade_spin", [cx, 0.2, cz]);
  blade.add(p.plastic, lathe(blade, [[0.0095, -0.02], [0.017, -0.02], [0.016, 0.06], [0.012, 0.075], [0.008, 0.082], [0, 0.083]], 24), [0, 0, 0]);
  const scim = new THREE.Shape();
  scim.moveTo(0, -0.011);
  scim.quadraticCurveTo(0.05, -0.016, 0.074, -0.002);
  scim.lineTo(0.074, 0.004);
  scim.quadraticCurveTo(0.045, 0.002, 0, 0.011);
  scim.lineTo(0, -0.011);
  for (const [y, a, tilt] of [
    [-0.004, 0, 0.1],
    [0.022, Math.PI, -0.1],
  ] as const) {
    const g = plate(blade, scim, 0.0016, 0, 10);
    g.rotateX(-Math.PI / 2);
    blade.add(p.ssBright, g, [0, y, 0], [0, a, tilt]);
  }
  // The lid on its hinge at the bowl's back: moulded clear with its feed tube and the tab that works the safety rod,
  // the pusher riding in the tube. Everything on the lid is built about the hinge pin.
  const H: V3 = [cx, yt + 0.004, cz - 0.104];
  const L = (q: V3): V3 => [q[0] - H[0], q[1] - H[1], q[2] - H[2]];
  // The bowl's two hinge knuckles either side of the lid's, moulded on its rim.
  for (const sx of [-1, 1]) {
    s.add(p.plastic, cyl(s, 0.0055, 0.014, "x", 16), [cx + sx * 0.0175, H[1], H[2]]);
    s.add(p.plastic, rbox(s, 0.014, 0.012, 0.01, 0.003), [cx + sx * 0.0175, yt - 0.004, cz - 0.098]);
  }
  if (s.fine(2)) s.add(p.ssBright, cyl(s, 0.0018, 0.05, "x", 10), [H[0], H[1], H[2]]);
  const lid = s.joint("lid", H);
  lid.add(p.glass, lathe(lid, [[0.0, yt + 0.024], [0.07, yt + 0.02], [0.094, yt + 0.008], [0.1, yt], [0.1, yt - 0.006], [0.097, yt - 0.006], [0.097, yt + 0.002], [0.09, yt + 0.006], [0.068, yt + 0.016], [0, yt + 0.02]], 48), L([cx, 0, cz]));
  lid.add(p.glass, cyl(lid, 0.0052, 0.02, "x", 16), [0, 0, 0]);
  lid.add(p.glass, rbox(lid, 0.02, 0.008, 0.012, 0.003), [0, -0.0008, 0.006]);
  const fx = cx - 0.035;
  const fz = cz + 0.03;
  lid.add(p.glass, lathe(lid, [[0.03, yt + 0.01], [0.031, 0.385], [0.036, 0.39], [0.036, 0.393], [0.027, 0.393], [0.026, yt + 0.012]], 32), L([fx, 0, fz]));
  // The safety rod's channel moulded down the bowl's back right, and the short tab on the lid that presses into its
  // top when the lid is shut: a tab reaching down the channel would swing through the bowl as the lid lifts.
  s.add(p.glass, rbox(s, 0.016, 0.115, 0.012, 0.004), [cx + 0.07, yt - 0.0575, cz - 0.07], [0, Math.PI / 4, 0]);
  lid.add(p.glass, new THREE.BoxGeometry(0.016, 0.016, 0.012), L([cx + 0.07, yt + 0.007, cz - 0.07]), [0, Math.PI / 4, 0]);
  const pusher = lid.joint("pusher", L([fx, 0.393, fz]));
  pusher.add(p.hdpe, cyl(pusher, 0.024, 0.075, "y", 24), [0, -0.035, 0]);
  pusher.add(p.hdpe, lathe(pusher, [[0, 0], [0.036, 0], [0.037, 0.006], [0.03, 0.014], [0, 0.016]], 28), [0, 0.0005, 0]);
  return (b) => {
    let a = 0;
    const pj = b.joints.get("pusher");
    const rest: V3 = pj ? [pj.position.x, pj.position.y, pj.position.z] : L([fx, 0.393, fz]);
    return {
      tick: (_t, dt, load) => {
        a += dt * 45 * Math.max(0, Math.min(1, load));
        turnJoint(b.joints.get("s_blade_spin"), "y", a);
      },
      // Open: the lid swings up on its hinge, front first, the pusher seated in its tube and drawn up a little.
      open: (f) => {
        const k = Math.max(0, Math.min(1, f));
        turnJoint(b.joints.get("lid"), "x", -0.62 * k);
        slideJoint(pj, "y", k * 0.022, rest);
      },
    };
  };
};

/**
 * A 12 in gravity-feed slicer: a cast base on feet, the motor housing behind the blade, the hollow-ground knife in its
 * ring guard under a centre plate, the gauge plate before it set by the index knob, the carriage tray tilted against the
 * knife on its slide with the food pusher on its rod, and the two-stone sharpener stowed on top of the guard.
 */
const rotarySlicer: Builder = (s, e) => {
  const p = palette();
  const { w, d } = e;
  const R = 0.152;
  // A gravity-feed slicer: the knife leans 30° so that the product table, square to it, falls toward the gauge plate.
  const tilt = 0.52;
  const [ct, st] = [Math.cos(tilt), Math.sin(tilt)];
  const bc: V3 = [0.3, 0.27, -0.33];
  const f2 = s.fine(2);
  // Rounded boxes at full detail; plain ones from further off, where 1 mm rounds are lost.
  const rb = (sh: Shop, x: number, y: number, z: number, r: number) => (f2 ? rbox(sh, x, y, z, r) : new THREE.BoxGeometry(x, y, z));
  // The product table's centre, in the knife's frame 11 cm out from its plane and 13 cm under its axis.
  const piv: V3 = [bc[0] + 0.112 * ct + 0.134 * st, bc[1] + 0.112 * st - 0.134 * ct, -0.24];
  const C = (q: V3): V3 => [q[0] - piv[0], q[1] - piv[1], q[2] - piv[2]];
  const slicerHooks = (b: Built): Hooks => {
    let a = 0;
    let pos = 0;
    let thick = 0.3;
    const cj = b.joints.get("carriage_slide");
    const rest: V3 = cj ? [cj.position.x, cj.position.y, cj.position.z] : piv;
    const gj = b.joints.get("gauge_plate");
    const grest: V3 = gj ? [gj.position.x, gj.position.y, gj.position.z] : [0, 0, 0];
    const place = () => slideJoint(cj, "z", -pos * 0.15, rest);
    const gauge = () => {
      slideJoint(gj, "x", -thick * 0.02, grest);
      turnJoint(b.joints.get("index_knob"), "z", -thick * Math.PI * 1.2);
    };
    place();
    gauge();
    let tt = 0;
    return {
      tick: (_t, dt, load) => {
        const l = Math.max(0, Math.min(1, load));
        a += dt * 30 * l;
        turnJoint(b.joints.get("rotary_blade_spin"), "x", a);
        if (l > 0) {
          tt += dt * l;
          pos = Math.sin(tt * Math.PI * 1.2);
          place();
        }
      },
      set: (name, v) => {
        if (name === "carriage") {
          pos = Math.max(-1, Math.min(1, v));
          place();
        } else if (name === "thickness") {
          thick = Math.max(0, Math.min(1, v));
          gauge();
        }
      },
    };
  };
  // The motor housing: one casting from the base up behind the knife, its face inclined parallel to the knife's
  // back, its top and left side rounded over in long smooth curves.
  const faceX = (y: number) => bc[0] + (-0.035 - (y - bc[1]) * st) / ct;
  const housing = new THREE.Shape();
  housing.moveTo(0.06, 0.1);
  housing.lineTo(faceX(0.1), 0.1);
  housing.lineTo(faceX(0.392), 0.392);
  housing.quadraticCurveTo(faceX(0.392) - 0.022, 0.438, 0.125, 0.437);
  housing.quadraticCurveTo(0.064, 0.436, 0.06, 0.37);
  housing.lineTo(0.06, 0.1);
  const hz0 = -0.565;
  const hz1 = -0.125;
  // Close up its edges are rounded over in six steps, so the big radii read as cast, not as facets.
  const hb = f2 ? 0.014 : 0;
  const housingGeo = (sh: Shop) => {
    if (!f2) return creased(plate(sh, housing, hz1 - hz0, 0, 10), 0.86);
    const g = new THREE.ExtrudeGeometry(housing, { depth: hz1 - hz0 - 2 * hb, bevelEnabled: true, bevelThickness: hb, bevelSize: hb, bevelSegments: 6, curveSegments: 28 });
    g.translate(0, 0, hb);
    return creased(g, 0.5);
  };
  // The housing's outer side, where the bevel leaves it.
  const wx = 0.06 - hb;
  if (s.detail === 0) {
    // Far off: base and housing, the knife disc and its guard ring, gauge plate, carriage tray and handle, the knob.
    blk(s, p.alu, [0.02, 0, -d + 0.02], [w - 0.02, 0.1, -0.02], 0);
    s.add(p.alu, housingGeo(s), [0, 0, hz0]);
    s.at(bc, [0, 0, tilt], () => {
      s.joint("rotary_blade_spin", [0, 0, 0]).add(p.ssBright, alongAxis(new THREE.CylinderGeometry(R + 0.009, R + 0.009, 0.008, 14), "x"), [0.002, 0, 0]);
      s.joint("gauge_plate", [0, 0, 0]).add(p.alu, new THREE.BoxGeometry(0.012, 0.14, 0.19), [-0.011, -0.07, -0.095]);
    });
    const car = s.joint("carriage_slide", piv);
    car.at([0, 0, 0], [0, 0, tilt], () => car.add(p.alu, new THREE.BoxGeometry(0.21, 0.03, 0.26), [0, 0, 0]));
    car.add(p.alu, new THREE.BoxGeometry(0.02, 0.2, 0.02), C([0.47, 0.24, -0.08]));
    s.joint("index_knob", [0.12, 0.065, -0.02]).add(p.bakelite, new THREE.CylinderGeometry(0.03, 0.03, 0.025, 8).rotateX(Math.PI / 2), [0, 0, 0.0125]);
    return slicerHooks;
  }
  // Feet and the cast base.
  for (const x of [0.06, w - 0.06]) for (const z of [-0.06, -d + 0.06]) s.add(p.rubber, lathe(s, f2 ? [[0, 0], [0.026, 0], [0.024, 0.02], [0.018, 0.03], [0, 0.03]] : [[0, 0], [0.026, 0.002], [0.02, 0.03], [0, 0.03]], 20), [x, 0, z]);
  blk(s, p.alu, [0.02, 0.03, -d + 0.02], [w - 0.02, 0.1, -0.02], 0.03);
  s.add(p.alu, housingGeo(s), [0, 0, hz0]);
  // The motor's end bell on the housing's outer side, its vents and the cord gland.
  s.add(p.alu, lathe(s, [[0, 0], [0.085, 0], [0.085, 0.012], [0.07, 0.03], [0.03, 0.04], [0, 0.04]], f2 ? 48 : 40).rotateZ(Math.PI / 2), [wx + 0.0005, 0.27, -0.36]);
  // Cooling slots in the side under the end bell.
  if (f2) for (let i = 0; i < 6; i++) s.add(p.enamel, rb(s, 0.002, 0.0045, 0.07, 0.0015), [wx - 0.0009, 0.125 + i * 0.009, -0.36]);
  s.add(p.rubber, cyl(s, 0.009, 0.03, "z", 12), [0.12, 0.13, hz0 - 0.012]);
  // The receiving tray on the base's front left, where the slices land: pressed stainless, its edges turned up.
  {
    const [rx0, rx1, rz0, rz1] = [0.045, 0.3, -0.115, -0.032];
    const path = filletPath([[rx0, rz1], [rx1, rz1], [rx1, rz0], [rx0, rz0]], 0.012, f2 ? 5 : 1, true);
    s.add(p.ssDull, pathSweep(sheetOutline([[0.008, 0], [0.003, 0], ...arc(0.003, 0.003, 0.003, -Math.PI / 2, -Math.PI, f2 ? 3 : 1).slice(1), [0, 0.008], ...arc(0.0015, 0.008, 0.0015, Math.PI, 0, f2 ? 4 : 1).slice(1)], 0.001), path, 0.1006, true));
    const tg = plate(s, planShape(offsetPath(path, 0.008, true)), 0.001);
    tg.rotateX(-Math.PI / 2);
    s.add(p.ssDull, tg, [0, 0.1001, 0]);
  }
  dataPlate(s, [0.3, 0.065, -d + 0.0198], [0, Math.PI, 0], 0.07, 0.035);
  // Start and stop on the base's front right.
  s.add(p.ledGreen, cyl(s, 0.011, 0.008, "z", 18), [0.44, 0.065, -0.016]);
  s.add(p.rubber, cyl(s, 0.011, 0.008, "z", 18), [0.48, 0.065, -0.016]);
  // The index knob: a big fluted knob over a graduated ring on the base's front left.
  if (s.fine(1)) for (let i = 0; i < 9; i++) {
    const a = Math.PI * 1.1 - (i / 8) * Math.PI * 1.2;
    s.add(p.enamel, new THREE.BoxGeometry(0.0015, 0.007, 0.0008), [0.12 + Math.cos(a) * 0.038, 0.065 + Math.sin(a) * 0.038, -0.0196], [0, 0, a - Math.PI / 2]);
  }
  const ik = s.joint("index_knob", [0.12, 0.065, -0.02]);
  if (f2) knob(ik, p.bakelite, p.ssBright, [0, 0, 0], 0.03);
  else ik.add(p.bakelite, knobGeo(ik, 0.03).rotateX(Math.PI / 2), [0, 0, 0]);

  // The knife and everything in its plane, in a frame turned so its axis (local x) leans over the carriage.
  s.at(bc, [0, 0, tilt], () => {
    // Hub and bearing boss from the housing.
    s.add(p.alu, cyl(s, 0.05, 0.034, "x", f2 ? 40 : 32), [-0.022, 0, 0]);
    const knife = s.joint("rotary_blade_spin", [0, 0, 0]);
    const kp: [number, number][] = f2 ? [[0, -0.003], [0.11, -0.003], [R - 0.004, -0.0007], [R, 0], [R - 0.004, 0.0007], [0.11, 0.0022], [0.05, 0.004], [0, 0.004]] : [[0, -0.003], [R - 0.004, -0.001], [R, 0], [R - 0.004, 0.0012], [0, 0.004]];
    knife.add(p.ssBright, lathe(knife, kp, f2 ? 96 : 72).rotateZ(-Math.PI / 2), [0.001, 0, 0]);
    if (knife.fine(1)) for (let i = 0; i < 4; i++) {
      const a = (i / 4) * Math.PI * 2;
      knife.add(p.ssDull, cyl(knife, 0.006, 0.003, "x", 12), [0.004, Math.cos(a) * 0.085, Math.sin(a) * 0.085]);
    }
    // The centre plate over the knife's middle, its removal knob.
    s.add(p.alu, lathe(s, [[0, 0.016], [0.06, 0.013], [0.1, 0.006], [0.116, 0.0], [0, 0]], f2 ? 64 : 48).rotateZ(-Math.PI / 2), [0.007, 0, 0]);
    s.add(p.bakelite, lathe(s, [[0, 0], [0.022, 0], [0.024, 0.008], [0.016, 0.022], [0, 0.024]], 24).rotateZ(-Math.PI / 2), [0.022, 0, 0]);
    // The ring guard round the edge from the front top over to the back, open where the carriage passes, on its
    // mounting block from the housing.
    const a0 = -0.35;
    const a1 = Math.PI + 0.9;
    const ring = torus(s, R + 0.009, 0.007, "x", a1 - a0, f2 ? 96 : 64, f2 ? 10 : 8);
    ring.rotateX(-a0);
    s.add(p.alu, ring, [0.004, 0, 0]);
    s.add(p.alu, rb(s, 0.046, 0.04, 0.05, 0.01), [-0.02, R + 0.016, 0]);
    // The gauge plate before the knife, sliding on its axis as the knob sets the slice, on a boss from the housing.
    const gp = s.joint("gauge_plate", [0.0, 0, 0]);
    const gs = new THREE.Shape();
    gs.moveTo(0.0, -0.14);
    gs.lineTo(0.19, -0.14);
    gs.lineTo(0.2, 0.01);
    gs.absarc(0, 0, 0.157, 0.05, -0.9, true);
    gs.lineTo(0, -0.14);
    const gpg = plate(gp, gs, 0.012, f2 ? 0.002 : 0, 24);
    gpg.rotateY(Math.PI / 2);
    gp.add(p.alu, gpg, [-0.011, 0, 0]);
    gp.add(p.alu, cyl(gp, 0.013, 0.036, "x", 20), [-0.029, -0.085, -0.11]);
    // The sharpener stowed on top: a cast housing with its two stones and knob.
    s.at([0.01, R + 0.05, -0.02], [0, 0, 0], () => {
      s.add(p.alu, rb(s, 0.07, 0.04, 0.09, 0.012), [0, 0, 0]);
      for (const z of [-0.025, 0.025]) s.add(p.ceramic, cyl(s, 0.014, 0.012, "x", 16), [0.03, -0.012, z]);
      s.add(p.bakelite, knobGeo(s, 0.016), [0, 0.02, 0]);
    });
  });

  // The carriage: an arm from the slide on the base's right up to the product table, inclined with the knife so that
  // the product slides down it onto the gauge plate; its end weight on a rod, and the handle at its front.
  s.add(p.ssBright, cyl(s, 0.012, d - 0.12, "z", f2 ? 20 : 16), [0.45, 0.115, -d / 2]);
  for (const z of [-0.07, -d + 0.07]) s.add(p.alu, rb(s, 0.04, 0.03, 0.03, 0.006), [0.45, 0.112, z]);
  const car = s.joint("carriage_slide", piv);
  car.add(p.alu, rb(car, 0.05, 0.04, 0.08, 0.008), C([0.45, 0.13, -0.24]));
  car.add(p.alu, rb(car, 0.032, 0.075, 0.05, 0.008), C([0.45, 0.18, -0.24]));
  car.at([0, 0, 0], [0, 0, tilt], () => {
    car.add(p.alu, rb(car, 0.09, 0.014, 0.12, 0.005), [0, -0.011, 0]);
    car.add(p.alu, rb(car, 0.21, 0.008, 0.26, 0.003), [0.0, 0, 0]);
    for (const sz of [-1, 1]) car.add(p.alu, rb(car, 0.21, 0.03, 0.006, 0.002), [0, 0.014, sz * 0.128]);
    car.add(p.alu, rb(car, 0.006, 0.05, 0.26, 0.002), [0.102, 0.022, 0]);
    // Ribs pressed in the table's face, running downhill, so the product slides.
    if (car.fine(2)) for (let i = 0; i < 7; i++) car.add(p.alu, rb(car, 0.19, 0.0016, 0.004, 0.0008), [0, 0.0046, -0.09 + i * 0.03]);
    // The end weight: a sleeve on a rod over the table, a spiked plate that bears the product toward the knife, a grip.
    car.add(p.ssBright, cyl(car, 0.006, 0.2, "x", 12), [0.0, 0.13, -0.11]);
    car.add(p.alu, rb(car, 0.012, 0.12, 0.012, 0.004), [0.1, 0.065, -0.11]);
    car.add(p.alu, rb(car, 0.04, 0.025, 0.025, 0.006), [0.02, 0.13, -0.11]);
    car.add(p.alu, rb(car, 0.008, 0.08, 0.12, 0.003), [0.0, 0.08, -0.05]);
    car.add(p.alu, rb(car, 0.02, 0.012, 0.07, 0.004), [0.015, 0.125, -0.08]);
    car.add(p.bakelite, knobGeo(car, 0.016), [0.02, 0.15, -0.11]);
    if (car.fine(2)) for (let i = 0; i < 4; i++) for (let k = 0; k < 3; k++) car.add(p.ssBright, cyl(car, 0.0016, 0.01, "x", 6, 0.0003), [-0.008, 0.055 + k * 0.02, -0.09 + i * 0.026]);
  });
  // The handle: forward from the slide block and up in front of the table, a black grip on top.
  car.add(p.alu, pipe(car, [C([0.45, 0.13, -0.2]), C([0.45, 0.13, -0.085]), C([0.47, 0.34, -0.08])], 0.008, 0.03, 12));
  car.add(p.bakelite, lathe(car, [[0, 0], [0.016, 0], [0.019, 0.02], [0.016, 0.04], [0, 0.043]], 20), C([0.47, 0.335, -0.08]));

  return slicerHooks;
};
/** A small domed knob, its base at y = 0. */
function knobGeo(s: Shop, r: number): THREE.BufferGeometry {
  return lathe(s, [[0, 0], [r * 0.7, 0], [r, r * 0.4], [r * 0.95, r * 0.9], [r * 0.6, r * 1.25], [0, r * 1.35]], 20);
}

/** A seven-segment readout drawn on a canvas: the digits of a weight, green on black. */
function drawReadout(c: HTMLCanvasElement, kg: number): void {
  const x = c.getContext("2d")!;
  x.fillStyle = "#050807";
  x.fillRect(0, 0, c.width, c.height);
  const text = Math.max(0, Math.min(99.99, kg)).toFixed(2).padStart(5, " ");
  const seg: Record<string, string> = { "0": "abcdef", "1": "bc", "2": "abged", "3": "abgcd", "4": "fgbc", "5": "afgcd", "6": "afgedc", "7": "abc", "8": "abcdefg", "9": "abcdfg", " ": "" };
  const H = c.height * 0.7;
  const W = H * 0.5;
  const t = H * 0.11;
  const y0 = (c.height - H) / 2;
  let px = c.width - 12 - 5 * (W + 14) + 6;
  const draw = (on: boolean, rx: number, ry: number, rw: number, rh: number) => {
    x.fillStyle = on ? "#5dff8a" : "#0c1a12";
    x.fillRect(rx, ry, rw, rh);
  };
  for (const ch of text) {
    if (ch === ".") {
      x.fillStyle = "#5dff8a";
      x.fillRect(px - 10, y0 + H - t, t, t);
      continue;
    }
    const on = seg[ch] ?? "";
    draw(on.includes("a"), px + t, y0, W - 2 * t, t);
    draw(on.includes("g"), px + t, y0 + H / 2 - t / 2, W - 2 * t, t);
    draw(on.includes("d"), px + t, y0 + H - t, W - 2 * t, t);
    draw(on.includes("f"), px, y0 + t * 0.6, t, H / 2 - t);
    draw(on.includes("b"), px + W - t, y0 + t * 0.6, t, H / 2 - t);
    draw(on.includes("e"), px, y0 + H / 2 + t * 0.4, t, H / 2 - t);
    draw(on.includes("c"), px + W - t, y0 + H / 2 + t * 0.4, t, H / 2 - t);
    px += W + 14;
  }
  x.fillStyle = "#5dff8a";
  x.font = `${Math.round(c.height * 0.22)}px sans-serif`;
  x.fillText("kg", c.width - 34, c.height - 10);
}

/**
 * A washdown portion scale: a stainless housing on levelling feet with a bubble level, a removable platter on its
 * spider over the load cell, and a sloped display panel, its seven-segment readout live, three membrane keys under it.
 */
const portionScale: Builder = (s, e) => {
  const p = palette();
  const { w, d, h } = e;
  const far = s.detail === 0;
  if (!far) for (const x of [0.04, w - 0.04]) for (const z of [-0.06, -d + 0.04]) {
    s.add(p.rubber, cyl(s, 0.012, 0.006, "y", 16), [x, 0.003, z]);
    s.add(p.ssBright, cyl(s, 0.004, 0.01, "y", 10), [x, 0.01, z]);
  }
  blk(s, p.ss, [0.01, 0.012, -d + 0.008], [w - 0.01, 0.088, -0.045], 0.008);
  // The display wedge at the front: its face sloping back, the readout and keys.
  const wedge = new THREE.Shape([new THREE.Vector2(0, 0), new THREE.Vector2(-0.04, 0), new THREE.Vector2(-0.032, 0.058), new THREE.Vector2(0, 0.076)]);
  const wg = plate(s, wedge, w - 0.02, s.fine(1) ? 0.0015 : 0, 4);
  if (far) blk(s, p.ss, [0.01, 0, -d + 0.008], [w - 0.01, 0.012, -0.005], 0);
  wg.rotateY(Math.PI / 2);
  s.add(p.ss, wg, [0.01, 0.012, -0.045]);
  const slope = Math.atan2(0.008, 0.058);
  const canvas = document.createElement("canvas");
  canvas.width = 256;
  canvas.height = 72;
  drawReadout(canvas, 0);
  const tex = new THREE.CanvasTexture(canvas);
  tex.colorSpace = THREE.SRGBColorSpace;
  const lcd = new THREE.MeshBasicMaterial({ map: tex });
  lcd.userData.uv = "own";
  s.at([w / 2 - 0.04, 0.048, -0.0086], [-slope, 0, 0], () => {
    s.add(p.screen, rbox(s, 0.13, 0.034, 0.002, 0.001), [0, 0.004, 0.0015]);
    s.add(lcd, new THREE.PlaneGeometry(0.12, 0.028), [0, 0.004, 0.0028]);
    // Three membrane keys (zero, tare, units) right of the readout, a power lamp over the first.
    if (!far) for (let i = 0; i < 3; i++) s.add(p.rubber, rbox(s, 0.016, 0.011, 0.003, 0.002), [0.1 + i * 0.022, -0.002, 0.0012]);
    if (s.fine(1)) s.add(p.ledGreen, cyl(s, 0.0025, 0.002, "z", 10), [0.1, 0.012, 0.0012]);
  });
  // The load cell's sealed boss through the housing top, its spider and the rubber pads the platter sits on.
  const pz = -0.04 - 0.145;
  const pw = 0.29;
  if (far) {
    blk(s, p.ss, [w / 2 - pw / 2, h - 0.01, pz - pw / 2], [w / 2 + pw / 2, h, pz + pw / 2], 0);
    s.add(p.ssDull, new THREE.CylinderGeometry(0.03, 0.03, h - 0.01 - 0.088, 6, 1, true), [w / 2, (h - 0.01 + 0.088) / 2, pz]);
    return () => ({
      set: (name, v) => {
        if (name !== "weight") return;
        drawReadout(canvas, v);
        tex.needsUpdate = true;
      },
    });
  }
  s.add(p.ssDull, cyl(s, 0.03, 0.013, "y", 24), [w / 2, 0.0945, pz]);
  if (s.fine(1)) s.add(p.rubber, torus(s, 0.031, 0.002, "y", Math.PI * 2, 24, 6), [w / 2, 0.089, pz]);
  for (const a of [0, Math.PI / 2]) s.add(p.ssDull, rbox(s, 0.2, 0.006, 0.014, 0.002), [w / 2, 0.1045, pz], [0, a + Math.PI / 4, 0]);
  for (let i = 0; i < 4; i++) {
    const a = (i / 4) * Math.PI * 2;
    s.add(p.rubber, cyl(s, 0.006, 0.004, "y", 10), [w / 2 + Math.cos(a) * 0.095, 0.1095, pz + Math.sin(a) * 0.095]);
  }
  // The platter: one pressed sheet, its edges rolled over and turned down all round, the corners drawn round.
  {
    const [x0, x1, z0, z1] = [w / 2 - pw / 2, w / 2 + pw / 2, pz - pw / 2, pz + pw / 2];
    const f2 = s.fine(2);
    const path = filletPath([[x0, z1], [x1, z1], [x1, z0], [x0, z0]], 0.007, f2 ? 5 : 2, true);
    const rr = 0.0025;
    s.add(p.ss, pathSweep(sheetOutline([[rr + 0.0006, 0], ...arc(rr + 0.0006, -rr, rr, Math.PI / 2, Math.PI, f2 ? 5 : 2).slice(1), [0.0006, -0.011]], 0.0012), path, h - 0.0006, true));
    const top = plate(s, planShape(offsetPath(path, rr + 0.0006, true)), 0.0012);
    top.rotateX(-Math.PI / 2);
    s.add(p.ss, top, [0, h - 0.0012, 0]);
  }
  return () => ({
    set: (name, v) => {
      if (name !== "weight") return;
      drawReadout(canvas, v);
      tex.needsUpdate = true;
    },
  });
};

/** The box grater's top width over its foot's, which its face patterns are drawn to. */
const GRATER_TAPER = 0.71;
/** Grater faces: alpha-tested punched sheet with raised cutting lips in the bump map, one material per pattern. */
const graterCache = new Map<string, THREE.MeshStandardMaterial>();
function graterFace(kind: "coarse" | "medium" | "fine" | "slicer"): THREE.MeshStandardMaterial {
  const hit = graterCache.get(kind);
  if (hit) return hit;
  const W = 256;
  const H = 512;
  const mk = () => {
    const c = document.createElement("canvas");
    c.width = W;
    c.height = H;
    return [c, c.getContext("2d")!] as const;
  };
  const [ca, xa] = mk();
  const [cb, xb] = mk();
  xa.fillStyle = "#fff";
  xa.fillRect(0, 0, W, H);
  xb.fillStyle = "#808080";
  xb.fillRect(0, 0, W, H);
  const hole = (cx: number, cy: number, rx: number, ry: number) => {
    xa.fillStyle = "#000";
    xa.beginPath();
    xa.ellipse(cx, cy, rx, ry, 0, 0, Math.PI * 2);
    xa.fill();
    // Where the lip was punched up out of the sheet, a soft swell just under the hole (the lip itself is modelled).
    const g = xb.createRadialGradient(cx, cy + ry * 0.6, ry * 0.4, cx, cy + ry * 0.6, ry * 1.6);
    g.addColorStop(0, "rgba(200,200,200,0.9)");
    g.addColorStop(1, "rgba(128,128,128,0)");
    xb.fillStyle = g;
    xb.beginPath();
    xb.ellipse(cx, cy + ry * 0.6, rx * 1.4, ry * 1.4, 0, 0, Math.PI);
    xb.fill();
  };
  if (kind === "slicer") {
    const margin = 18;
    const half = (cy: number) => ((GRATER_TAPER + (1 - GRATER_TAPER) * (cy / H)) * W) / 2 - margin;
    for (const cy of [H * 0.38, H * 0.62]) {
      const hw = half(cy);
      xa.fillStyle = "#000";
      xa.fillRect(W / 2 - hw + 20, cy - 10, 2 * hw - 40, 20);
      xb.fillStyle = "#f0f0f0";
      xb.fillRect(W / 2 - hw + 16, cy + 10, 2 * hw - 32, 8);
    }
  } else
    for (const [cx, cy, rx, ry] of graterHoles(kind)) {
      hole(cx, cy, rx, ry);
      if (kind === "fine") {
        xb.fillStyle = "#f2f2f2";
        xb.fillRect(cx - 4, cy + 3, 8, 3);
      }
    }
  const ta = new THREE.CanvasTexture(ca);
  const tb = new THREE.CanvasTexture(cb);
  for (const t of [ta, tb]) t.anisotropy = 8;
  const m = new THREE.MeshStandardMaterial({ color: 0xd3d7dc, metalness: 1, roughness: 0.3, alphaMap: ta, alphaTest: 0.5, bumpMap: tb, bumpScale: 2.2, side: THREE.DoubleSide });
  m.alphaToCoverage = true;
  m.userData.uv = "own";
  graterCache.set(kind, m);
  return m;
}

/** The punched holes of a grater face, [cx, cy, rx, ry] in its 256 × 512 texture's pixels (y down from the top). */
function graterHoles(kind: "coarse" | "medium" | "fine"): [number, number, number, number][] {
  const [W, H, margin] = [256, 512, 18];
  // The face tapers to GRATER_TAPER of its foot's width at the top; the pattern keeps its pitch and loses holes instead.
  const half = (cy: number) => ((GRATER_TAPER + (1 - GRATER_TAPER) * (cy / H)) * W) / 2 - margin;
  const [cols, rows, rx, ry] = kind === "coarse" ? [5, 11, 13, 9] : kind === "medium" ? [8, 18, 7, 5] : [16, 36, 2.6, 2.6];
  const out: [number, number, number, number][] = [];
  for (let j = 0; j < rows; j++)
    for (let i = 0; i < cols - (j % 2); i++) {
      const cx = margin + ((i + 0.5 + (j % 2) * 0.5) * (W - 2 * margin)) / cols;
      const cy = margin + ((j + 0.5) * (H - 2 * margin)) / rows;
      if (Math.abs(cx - W / 2) + rx * 1.2 > half(cy)) continue;
      out.push([cx, cy, rx, ry]);
    }
  return out;
}

/**
 * The cutting teeth of a shredding face whose quad runs a (foot left), b (foot right), c (top right), d (top left):
 * each the lower rim of its hole raised out of the sheet into a curved scoop, highest at its middle, its edge turned up
 * the face, open behind to the hole so that the shred passes through into the box.
 */
function graterTeeth(kind: "coarse" | "medium", a: V3, b: V3, c: V3, d: V3, n: number): THREE.BufferGeometry {
  const v3 = (q: V3) => new THREE.Vector3(q[0], q[1], q[2]);
  const [A, B, Cc, D] = [v3(a), v3(b), v3(c), v3(d)];
  const across = B.clone().sub(A);
  const up = Cc.clone().add(D).multiplyScalar(0.5).sub(A.clone().add(B).multiplyScalar(0.5));
  const ex = across.clone().normalize();
  const ey = up.clone().normalize();
  const nn = ex.clone().cross(ey).normalize();
  const hgt = kind === "coarse" ? 0.0034 : 0.0016;
  const pos: number[] = [];
  const idx: number[] = [];
  for (const [cx, cy, rx, ry] of graterHoles(kind)) {
    const ctr = A.clone().addScaledVector(across, cx / 256).addScaledVector(up, 1 - cy / 512);
    const [Rx, Ry] = [(rx / 256) * across.length(), (ry / 512) * up.length()];
    const base = pos.length / 3;
    for (let k = 0; k <= n; k++) {
      const th = Math.PI + (Math.PI * k) / n;
      const f = -Math.sin(th);
      const p0 = ctr.clone().addScaledVector(ex, Rx * Math.cos(th)).addScaledVector(ey, Ry * Math.sin(th));
      const p1 = p0.clone().addScaledVector(nn, hgt * (0.15 + 0.85 * f)).addScaledVector(ey, hgt * 0.55 * f).addScaledVector(ex, -Rx * 0.08 * Math.cos(th));
      pos.push(p0.x, p0.y, p0.z, p1.x, p1.y, p1.z);
      if (k < n) {
        const q = base + 2 * k;
        idx.push(q, q + 2, q + 1, q + 1, q + 2, q + 3);
      }
    }
  }
  const g = new THREE.BufferGeometry();
  g.setAttribute("position", new THREE.Float32BufferAttribute(pos, 3));
  // Every vertex samples the solid margin of the face's alpha map, so the teeth share the face's material and draw.
  g.setAttribute("uv", new THREE.Float32BufferAttribute(new Array((pos.length / 3) * 2).fill(0.02), 2));
  g.setIndex(idx);
  g.computeVertexNormals();
  return g;
}

/**
 * A quad from four corners (bottom-left, bottom-right, top-right, top-left). A trapezoid's UVs keep its true shape:
 * the foot spans u 0..1 and the top its share of that, centred, so a pattern is not stretched by the taper.
 */
function quad(a: V3, b: V3, c: V3, d: V3): THREE.BufferGeometry {
  const g = new THREE.BufferGeometry();
  g.setAttribute("position", new THREE.Float32BufferAttribute([...a, ...b, ...c, ...d], 3));
  const dist = (p: V3, q: V3) => Math.hypot(p[0] - q[0], p[1] - q[1], p[2] - q[2]);
  const k = dist(c, d) / dist(a, b);
  g.setAttribute("uv", new THREE.Float32BufferAttribute([0, 0, 1, 0, (1 + k) / 2, 1, (1 - k) / 2, 1], 2));
  g.setIndex([0, 1, 2, 0, 2, 3]);
  g.computeVertexNormals();
  return g;
}

/**
 * A four-sided box grater: a tapered stainless box with a different cutting face on each side, coarse and medium
 * shredders, a fine rasp and a slicing slot, its corners folded and seamed, a rolled foot with a rubber grip, a top cap
 * and a loop handle with its grip.
 */
const boxGrater: Builder = (s, e) => {
  const p = palette();
  const { w, d } = e;
  const [bw, bd, tw, td] = [w - 0.006, d - 0.006, 0.078, 0.062];
  const [y0, y1] = [0.006, 0.19];
  const cx = w / 2;
  const cz = -d / 2;
  const B = (sx: number, sz: number): V3 => [cx + (sx * bw) / 2, y0, cz + (sz * bd) / 2];
  const T = (sx: number, sz: number): V3 => [cx + (sx * tw) / 2, y1, cz + (sz * td) / 2];
  const faces: ["coarse" | "medium" | "fine" | "slicer", V3, V3, V3, V3][] = [
    ["coarse", B(-1, 1), B(1, 1), T(1, 1), T(-1, 1)],
    ["medium", B(1, 1), B(1, -1), T(1, -1), T(1, 1)],
    ["fine", B(1, -1), B(-1, -1), T(-1, -1), T(1, -1)],
    ["slicer", B(-1, -1), B(-1, 1), T(-1, 1), T(-1, -1)],
  ];
  for (const [k, a, b, c, dd] of faces) {
    s.add(s.fine(1) ? graterFace(k) : p.ss, quad(a, b, c, dd));
    if (s.fine(2) && (k === "coarse" || k === "medium")) s.add(graterFace(k), graterTeeth(k, a, b, c, dd, k === "coarse" ? 10 : 6));
  }
  const hh = e.h - 0.006;
  if (s.detail === 0) {
    // Far off: the four faces, the cap, the loop of the handle.
    blk(s, p.ss, [cx - tw / 2, y1 - 0.001, cz - td / 2], [cx + tw / 2, y1 + 0.001, cz + td / 2], 0);
    s.add(p.ssBright, new THREE.TubeGeometry(new THREE.CatmullRomCurve3([new THREE.Vector3(cx - 0.03, y1, cz), new THREE.Vector3(cx - 0.028, hh - 0.012, cz), new THREE.Vector3(cx + 0.028, hh - 0.012, cz), new THREE.Vector3(cx + 0.03, y1, cz)]), 6, 0.004, 4, false));
    return;
  }
  // Corner seams up the four edges, the rolled foot with its rubber, the top cap and its rolled edge.
  for (const [sx, sz] of [
    [-1, 1],
    [1, 1],
    [1, -1],
    [-1, -1],
  ] as const)
    s.add(p.ss, pipe(s, [B(sx, sz), T(sx, sz)], 0.0022, 0.004, 6));
  const ring = (y: number, ww: number, dd: number, r: number, m: THREE.Material) => {
    for (const sz of [-1, 1]) s.add(m, cyl(s, r, ww, "x", 10), [cx, y, cz + (sz * dd) / 2]);
    for (const sx of [-1, 1]) s.add(m, cyl(s, r, dd, "z", 10), [cx + (sx * ww) / 2, y, cz]);
  };
  ring(y0, bw, bd, 0.0032, p.ss);
  ring(0.0025, bw - 0.002, bd - 0.002, 0.0026, p.rubber);
  ring(y1, tw, td, 0.0028, p.ss);
  blk(s, p.ss, [cx - tw / 2, y1 - 0.0008, cz - td / 2], [cx + tw / 2, y1 + 0.0008, cz + td / 2], 0.0004);
  // The handle: a stainless loop riveted to the cap, a black grip over its top.
  s.add(p.ssBright, pipe(s, [[cx - 0.03, y1, cz], [cx - 0.03, hh - 0.012, cz], [cx + 0.03, hh - 0.012, cz], [cx + 0.03, y1, cz]], 0.0035, 0.012, 8));
  s.add(p.plastic, lathe(s, [[0, -0.022], [0.006, -0.022], [0.0072, -0.016], [0.0072, 0.016], [0.006, 0.022], [0, 0.022]], 16).rotateZ(Math.PI / 2), [cx, hh - 0.012, cz]);
  for (const sx of [-1, 1]) s.add(p.ssBright, ball(s, 0.0035, 8, Math.PI / 2), [cx + sx * 0.03, y1 + 0.0008, cz]);
};

let boardMat: THREE.MeshStandardMaterial | null = null;
/** White HDPE with a season's faint knife scoring in its bump and roughness: one material for every board. */
function scoredHDPE(): THREE.MeshStandardMaterial {
  if (boardMat) return boardMat;
  const S = 512;
  const mk = (bg: string) => {
    const c = document.createElement("canvas");
    c.width = c.height = S;
    const x = c.getContext("2d")!;
    x.fillStyle = bg;
    x.fillRect(0, 0, S, S);
    return [c, x] as const;
  };
  // The same cuts drawn three ways: sunk in the bump, greyed a little in the colour, rougher in the roughness.
  const [cb, xb] = mk("#808080");
  const [cc, xc] = mk("#ffffff");
  const [cr, xr] = mk("#c8c8c8");
  let seed = 90127;
  const rnd = () => (seed = (seed * 1664525 + 1013904223) >>> 0) / 4294967296;
  for (let i = 0; i < 1100; i++) {
    const a = (rnd() < 0.6 ? 0 : Math.PI / 2) + (rnd() - 0.5) * 0.7;
    const len = 10 + rnd() * 90;
    const px = rnd() * S;
    const py = rnd() * S;
    const k = 0.25 + rnd() * 0.5;
    const lw = 0.7 + rnd() * 0.9;
    for (const [x, col] of [
      [xb, `rgba(20,20,20,${k.toFixed(2)})`],
      [xc, `rgba(170,168,160,${(k * 0.55).toFixed(2)})`],
      [xr, `rgba(255,255,255,${k.toFixed(2)})`],
    ] as const) {
      x.strokeStyle = col;
      x.lineWidth = lw;
      x.beginPath();
      x.moveTo(px, py);
      x.lineTo(px + Math.cos(a) * len, py + Math.sin(a) * len);
      x.stroke();
    }
  }
  const tex = (c: HTMLCanvasElement, srgb: boolean) => {
    const t = new THREE.CanvasTexture(c);
    t.wrapS = t.wrapT = THREE.RepeatWrapping;
    t.anisotropy = 8;
    if (srgb) t.colorSpace = THREE.SRGBColorSpace;
    return t;
  };
  boardMat = new THREE.MeshStandardMaterial({ color: 0xf1f0eb, map: tex(cc, true), roughness: 0.72, roughnessMap: tex(cr, false), metalness: 0, bumpMap: tex(cb, false), bumpScale: 2.5 });
  boardMat.userData.uv = "box";
  boardMat.userData.uvScale = 1 / 0.3;
  return boardMat;
}

/** A 24 × 18 × 3/4 in HDPE board: eased edges, rounded corners, a juice groove round its working face. */
const cuttingBoard: Builder = (s, e) => {
  const m = s.fine(1) ? scoredHDPE() : palette().hdpe;
  const { w, d, h } = e;
  const g0 = 0.003;
  const flat = (shape: THREE.Shape, t: number, bev: number, y: number) => {
    const g = plate(s, shape, t, s.fine(1) ? bev : 0, 10);
    g.rotateX(-Math.PI / 2);
    s.add(m, g, [w / 2, y, -d / 2]);
  };
  if (!s.fine(1)) {
    flat(roundRect(w - 0.003, d - 0.003, 0.012), h, 0.0015, 0);
    return;
  }
  // One solid board with the groove routed into it: the border full height round the groove, eased on both of its
  // edges; the field inside, full height and eased; the groove's floor between them, buried at its edges in both.
  const rim = roundRect(w - 0.003, d - 0.003, 0.012);
  rim.holes.push(rectPath(new THREE.Path(), w - 0.057, d - 0.057, 0.028, 0, 0, true) as THREE.Path);
  flat(rim, h, 0.0015, 0);
  flat(roundRect(w - 0.075, d - 0.075, 0.02), h, 0.0012, 0);
  const floor = plate(s, roundRect(w - 0.059, d - 0.059, 0.027), h - g0 - 0.002);
  floor.rotateX(-Math.PI / 2);
  s.add(m, floor, [w / 2, 0.002, -d / 2]);
};

register({
  "NSF-PROCESSOR-SBLADE": foodProcessor,
  "NSF-SLICER-ROTARY": rotarySlicer,
  "NSF-PORTION-SCALE-DIGI": portionScale,
  "NSF-GRATER-BOX": boxGrater,
  "NSF-CUTTING-BOARD-HDPE-24X18": cuttingBoard,
  "NSF-CARRIER-INSULATED": insulatedCarrier,
  "NSF-INSULATED-TRANSIT-CARRIER": insulatedCarrier,
  "NSF-DISHWASHER-DOOR": doorDishwasher,
  "NSF-DISHWASHER-CONVEYOR": conveyorDishwasher,
  "NSF-HANDWASH-SINK-WALL-17IN": handSink,
  "NSF-SINK-3COMP": threeCompartment,
  "NSF-3COMP-SINK-72IN": threeCompartment,
  "NSF-VEG-PREP-SINK-48IN": vegPrepSink,
  "NSF-PREPTABLE-48": workTable({ legs: 4, edge: "bull", splash: 0, shelf: "galv" }),
  "NSF-PREPTABLE-72": workTable({ legs: 4, edge: "bull", splash: 0.102, shelf: "ss" }),
  "NSF-PREPTABLE-96": workTable({ legs: 6, edge: "marine", splash: 0, shelf: "galv" }),
});

export type { V3 };
