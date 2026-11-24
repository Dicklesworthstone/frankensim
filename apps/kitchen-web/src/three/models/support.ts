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
  const n = s.seg(8, 3);
  return [[-back, 0], ...arc(-0.01, -0.01, 0.01, Math.PI / 2, 0, n), ...arc(-0.01, -0.031, 0.01, 0, -Math.PI / 2, n), [-0.03, -0.041], [-0.034, -0.037]];
}
function marineEdge(s: Shop, back = 0.03): P2[] {
  const n = s.seg(10, 3);
  return [[-back, 0], [-0.022, 0], ...arc(-0.022, 0.012, 0.012, -Math.PI / 2, -0.25, n).slice(1), ...arc(-0.006, 0.006, 0.006, Math.PI * 0.72, 0, n), [0, -0.034], ...arc(-0.006, -0.034, 0.006, 0, -Math.PI / 2, 3), [-0.026, -0.04]];
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
 * from the bottom edge, z out of the face (its back at -FT). A U-channel frame holds two interlocking rows of S-baffles,
 * front arches opening back and rear arches opening forward, so the air turns twice through every gap and throws its
 * grease out on the walls; it runs down the channels to weep holes in the bottom rail. Two wire handles lift it out.
 */
function baffleFilter(s: Shop): Part[] {
  const p = palette();
  const parts: Part[] = [];
  const fw = FW - 0.003;
  const hw = fw / 2;
  const fr = 0.025;
  if (s.detail === 0) {
    parts.push({ m: p.ss, g: new THREE.BoxGeometry(fw, FH - 0.004, FT * 0.6), at: [0, FH / 2, -FT * 0.3] });
    return parts;
  }
  const fine = s.fine(2);
  // Frame: side stiles and top and bottom rails, U-channels opening inward.
  const rc = fine ? 0.003 : 0;
  const stile: P2[] = fine ? [[fr, -0.0012], [rc, -0.0012], [0, -0.0012 - rc], [0, -FT + rc], [rc, -FT], [fr, -FT]] : [[fr, 0], [0, 0], [0, -FT], [fr, -FT]];
  for (const sx of [-1, 1]) {
    const g = sweep(s, stile.map(([x, z]) => [x * sx, z] as P2), 0.0012, FH, "y");
    parts.push({ m: p.ss, g, at: [-sx * hw, 0, 0] });
  }
  const rail = (y0: number, dir: number): P2[] => (fine ? [[-0.0012, y0 + dir * fr], [-0.0012, y0 + dir * rc], [-0.0012 - rc, y0], [-FT + rc, y0], [-FT, y0 + dir * rc], [-FT, y0 + dir * fr]] : [[0, y0 + dir * fr], [0, y0], [-FT, y0], [-FT, y0 + dir * fr]]);
  for (const [y0, dir] of [
    [0, 1],
    [FH, -1],
  ] as const) {
    const g = sweep(s, rail(y0, dir), 0.0012, fw - 0.004, "x");
    parts.push({ m: p.ss, g, at: [-hw + 0.002, 0, 0] });
  }
  // The baffles, between the rails.
  const inner = fw - 2 * fr;
  const nf = 9;
  const pitch = inner / nf;
  const bl = FH - 2 * fr + 0.01;
  const by = fr - 0.005;
  if (fine) {
    const rw = pitch * 0.4;
    const n = s.seg(8, 4);
    const front: P2[] = [[-rw + 0.003, -0.038], [-rw, -0.035], [-rw, -0.0215], ...arc(0, -0.0215, rw, Math.PI, 0, n).slice(1, -1), [rw, -0.0215], [rw, -0.035], [rw - 0.003, -0.038]];
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
    const weep = new THREE.CylinderGeometry(0.0028, 0.0028, 0.0016, 10);
    for (let i = 0; i < 6; i++) parts.push({ m: p.enamel, g: weep, at: [-inner / 2 + inner * ((i + 0.5) / 6), 0.0003, -FT * 0.5] });
    const rivet = ball(s, 0.0022, 8, Math.PI / 2);
    rivet.rotateX(Math.PI / 2);
    for (const sx of [-1, 1]) for (const y of [0.011, FH - 0.011]) parts.push({ m: p.ssBright, g: rivet, at: [sx * (hw - 0.011), y, -0.0012] });
  }
  return parts;
}

/** The lamps in a hood's globes: a live material, so that the lights can be dimmed. */
function lampMaterial(): THREE.MeshStandardMaterial {
  const m = new THREE.MeshStandardMaterial({ color: 0xfffaf2, roughness: 0.6, metalness: 0, emissive: 0xffeed2, emissiveIntensity: 1.6 });
  m.userData.uv = "own";
  return m;
}
const setLamp = (m: THREE.MeshStandardMaterial, v: number) => {
  m.emissiveIntensity = 1.6 * Math.max(0, Math.min(1, v));
};

/** A vapour-proof hood light hung from a ceiling whose underside is at y: cast base, threaded collar, globe, guard, lamp. */
function hoodLight(s: Shop, lamp: Mat, x: number, y: number, z: number): void {
  const p = palette();
  s.add(p.alu, hang(s, [
    [0, 0],
    [0.064, 0],
    [0.064, -0.006],
    [0.058, -0.02],
    [0.047, -0.026],
    [0, -0.026],
  ], 40), [x, y - 0.0005, z]);
  if (s.detail === 0) {
    s.add(lamp, cyl(s, 0.04, 0.09, "y", 16), [x, y - 0.07, z]);
    return;
  }
  // The threaded collar the globe screws into, ribbed close up.
  s.add(p.alu, cyl(s, 0.045, 0.012, "y", 32), [x, y - 0.031, z]);
  if (s.fine(2)) for (let i = 0; i < 3; i++) s.add(p.alu, torus(s, 0.0452, 0.0012, "y", Math.PI * 2, 32, 4), [x, y - 0.027 - i * 0.004, z]);
  // The globe: a heavy jelly-jar of heat-resistant glass, and the frosted lamp inside it.
  s.add(p.glass, hang(s, [
    [0.04, 0],
    [0.043, -0.008],
    [0.046, -0.05],
    [0.042, -0.092],
    [0.03, -0.108],
    [0, -0.111],
  ], 32), [x, y - 0.037, z]);
  s.add(lamp, hang(s, [
    [0.008, 0],
    [0.011, -0.012],
    [0.022, -0.035],
    [0.024, -0.05],
    [0.017, -0.066],
    [0, -0.07],
  ], 20), [x, y - 0.04, z]);
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
  for (const sz of [-1, 1]) blk(s, p.ssBright, [x - L / 2, y - 0.02, z + sz * 0.0118], [x + L / 2, y + 0.0105, z + sz * 0.0131], 0.0004);
  for (const sx of [-1, 1]) {
    blk(s, p.ssBright, [x + sx * (L / 2 - 0.0013), y - 0.02, z - 0.013], [x + sx * (L / 2 + 0.0002), y + 0.0105, z + 0.013], 0.0004);
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

/** A hanger bracket on a hood top at height y: a welded angle clip with a hole for the threaded rod. Upright along x. */
function hangerClip(s: Shop, x: number, y: number, z: number): void {
  const p = palette();
  blk(s, p.ss, [x - 0.04, y, z - 0.028], [x + 0.04, y + 0.003, z + 0.028], 0.0006);
  const sh = roundRect(0.056, 0.07, 0.012, 0, 0.035);
  if (s.fine(1)) sh.holes.push(new THREE.Path().absarc(0, 0.045, 0.009, 0, Math.PI * 2, true));
  const g = plate(s, sh, 0.003);
  s.add(p.ss, g, [x - 0.0015, y, z], [0, Math.PI / 2, 0]);
  if (s.fine(2)) for (const sz of [-1, 1]) screw(s, p.ssBright, [x + 0.025, y + 0.003, z + sz * 0.017], 0.003, [-Math.PI / 2, 0, 0]);
}

/**
 * A duct collar standing on a hood top at `y0` up to `y1`: four walls of 16-gauge, a welded flange on the top, the
 * upper edge turned out where the duct is welded on. [x0, x1] × [z0, z1] is its opening.
 */
function ductCollar(s: Shop, x0: number, x1: number, z0: number, z1: number, y0: number, y1: number): void {
  const p = palette();
  const t = 0.0016;
  const hgt = y1 - y0;
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
    const n = Math.max(2, Math.round((x1 - x0) / 0.15));
    for (let i = 0; i <= n; i++) for (const z of [z0 - f / 2, z1 + f / 2]) screw(s, p.ssBright, [x0 + ((x1 - x0) * i) / n, y0 + 0.002, z], 0.003, [-Math.PI / 2, 0, 0]);
  }
  void hgt;
}

/**
 * A Type I wall canopy, or the double island when `island`: an 18-gauge box canopy, double-walled at the ends, with a
 * front make-up air plenum whose perforated face throws supply air forward, a lit capture reservoir, and a 55° bank of
 * baffle filters over a grease trough at the back (on the island, two banks in a V back to back), the exhaust plenum
 * behind them and its duct collars on top. The wet-chemical system's piping hangs in the reservoir, nozzles aimed down
 * at the line, and its detection line runs fusible links in front of the filters.
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
    const lamp = lampMaterial();

    // The ends: double-wall panels, closed below by a return with a rolled safety hem.
    for (const [x0, sg] of [
      [0, 1],
      [w, -1],
    ] as const) {
      blk(s, p.ss, [x0, 0, -d], [x0 + sg * T, hb, 0], r);
      blk(s, p.ss, [x0 + sg * (wall - T), 0, -d + (island ? 0 : T)], [x0 + sg * wall, hb - T, 0], r);
      blk(s, p.ss, [x0, 0, -d], [x0 + sg * wall, T, 0], 0);
      if (s.fine(2)) s.add(p.ss, cyl(s, 0.0035, d, "z", 12), [x0 + sg * 0.0035, 0.0035, -d / 2]);
    }
    // The back panel against the wall.
    if (!island) blk(s, p.ss, [wall, 0, -d], [w - wall, hb - T, -d + T], r);

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
    const holes = [...exhausts, ...muas];
    sheetWithHoles(s, p.ss, 0, w, -d, 0, hb - T, hb, holes);
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

    // Hanger clips at the corners and every 2.4 m between.
    const nh = Math.max(1, Math.round(w / 2.4));
    for (let i = 0; i <= nh; i++) {
      const x = 0.1 + ((w - 0.2) * i) / nh;
      for (const z of [-0.09, -d + 0.09]) hangerClip(s, x, hb, z);
    }

    // A field joint on hoods too long to ship in one piece: bolted flanges inside, cover strips outside.
    if (w > 4.9) {
      const xj = w / 2;
      blk(s, p.ss, [xj - 0.02, 0.0, -0.0001], [xj + 0.02, hb, 0.0014], 0.0008);
      blk(s, p.ss, [xj - 0.02, hb, -d], [xj + 0.02, hb + 0.0016, 0], 0.0006);
      for (const sx of [-1, 1]) blk(s, p.ss, [xj + sx * 0.0015, hb - 0.03, -d + 0.002], [xj + sx * 0.0035, hb - T, -0.002], 0);
      if (s.fine(2))
        for (let z = -0.25; z > -d + 0.1; z -= 0.15) {
          s.add(p.ssBright, hex(0.006, 0.012, "x"), [xj, hb - 0.015, z]);
        }
    }

    // A plate on the left end's inner skin: the listing label, blank.
    dataPlate(s, [wall + 0.0002, hb * 0.55, -mua - 0.12], [0, Math.PI / 2, 0], 0.1, 0.06);

    // One side of the hood in its own frame: front at z = 0, its "wall" (the back wall, or the island's centre) at -D.
    const filter = baffleFilter(s);
    const phi = Math.PI / 2 - SLOPE;
    const cphi = Math.cos(phi);
    const yb = 0.1;
    const zb = -D + 0.08;
    const zpf = zb + 0.302;
    const zl = -mua - 0.17;
    const zf = zpf + 0.1;
    const yf = hb - 0.06;
    const zc = zpf + 0.035;
    const yc = hb - 0.035;
    const side = (k: number) => {
      // The make-up air plenum: inner wall, floor, front face with its perforated supply panels and diffuser behind.
      blk(s, p.ss, [wall, 0, -mua], [w - wall, hb - T, -mua + T], r);
      blk(s, p.ss, [wall, 0, -mua], [w - wall, T, -T], 0);
      blk(s, p.ss, [0, 0.305, -T], [w, hb, 0], r);
      blk(s, p.ss, [0, 0, -T], [w, 0.035, 0], r);
      if (s.fine(2)) s.add(p.ss, cyl(s, 0.0035, w, "x", 12), [w / 2, 0.0035, -0.0035]);
      const nSec = Math.max(1, Math.round((w - 2 * wall) / 1.22));
      const secW = (w - 2 * wall - 0.03) / nSec;
      const mull = (x: number, hw: number) => blk(s, p.ss, [x - hw, 0.034, -T], [x + hw, 0.306, 0], r);
      mull(wall / 2 + 0.0075, wall / 2 + 0.0075);
      mull(w - wall / 2 - 0.0075, wall / 2 + 0.0075);
      for (let i = 1; i < nSec; i++) mull(wall + 0.015 + secW * i, 0.015);
      for (let i = 0; i < nSec; i++) {
        const xa = wall + 0.03 + secW * i;
        const xb = xa + secW - 0.03;
        if (s.fine(1)) blk(s, perf, [xa - 0.006, 0.03, -0.0058], [xb + 0.006, 0.31, -0.0049], 0);
        else blk(s, p.ssDull, [xa, 0.035, -0.003], [xb, 0.305, -0.002], 0);
      }
      blk(s, p.charcoal, [wall, 0.028, -0.075], [w - wall, 0.312, -0.074], 0);
      if (s.fine(2)) for (let i = 0; i <= nSec; i++) for (const y of [0.018, 0.32]) screw(s, p.ssBright, [wall + 0.015 + secW * i, y, 0], 0.0028);

      // The capture reservoir's lights, wired through boxes on the top.
      const nl = Math.max(1, Math.round((w - 2 * wall) / 1.2));
      const lx: number[] = [];
      for (let i = 0; i < nl; i++) lx.push(wall + ((w - 2 * wall) * (i + 0.5)) / nl);
      for (const x of lx) {
        hoodLight(s, lamp, x, hb - T, zl);
        if (s.fine(1)) {
          blk(s, p.ssDull, [x - 0.05, hb, zl - 0.05], [x + 0.05, hb + 0.045, zl + 0.05], 0.002);
          if (s.fine(2)) for (const sx of [-1, 1]) for (const sz of [-1, 1]) screw(s, p.ssBright, [x + sx * 0.04, hb + 0.045, zl + sz * 0.04], 0.0022, [-Math.PI / 2, 0, 0]);
        }
      }
      if (s.fine(1)) {
        s.add(p.ssDull, cyl(s, 0.009, lx[lx.length - 1]! - 0.12, "x", 12), [(0.12 + lx[lx.length - 1]!) / 2, hb + 0.022, zl]);
        blk(s, p.ssDull, [0.07, hb, zl - 0.055], [0.17, hb + 0.055, zl + 0.055], 0.002);
      }

      // The filter bank: filters stamped across, blank fillers at the ends, the top channel and the plenum's front.
      const usable = w - 2 * wall - 0.01;
      const nf = Math.max(1, Math.floor(usable / FW));
      const x0 = (w - nf * FW) / 2;
      s.at([0, yb, zb], [phi, 0, 0], () => {
        for (let i = 0; i < nf; i++) stamp(s, filter, [x0 + FW * (i + 0.5), 0, 0]);
        for (const [a, b] of [
          [wall, x0],
          [x0 + nf * FW, w - wall],
        ] as const) {
          if (b - a < 0.005) continue;
          blk(s, p.ss, [a, 0.004, -0.004], [b, FH - 0.004, -0.0028], 0);
          blk(s, p.ss, [a, 0.004, -FT], [b, FH - 0.004, -0.004], 0);
        }
        // The top channel the filters lift into.
        sweepX(s, p.ss, [[0.012, FH - 0.024], [0.012, FH + 0.006], [-FT - 0.004, FH + 0.006], [-FT - 0.004, FH - 0.012]], 0.0012, wall, w - wall);
      });
      blk(s, p.ss, [wall, yb + FH * cphi - 0.012, zpf], [w - wall, hb - T, zpf + T], r);

      // The grease trough: a U channel under the bank, its front lip rolled in; a ledge inside the filters stand on.
      const zt = -D + 0.004;
      const nArc = s.seg(6, 2);
      sweepX(s, p.ss, [[zt, 0.135], [zt, 0.075], ...arc(zt + 0.015, 0.075, 0.015, Math.PI, Math.PI * 1.5, nArc).slice(1), [zt + 0.08, 0.06], ...arc(zt + 0.08, 0.075, 0.015, -Math.PI / 2, 0, nArc).slice(1), [zt + 0.095, 0.11], [zt + 0.088, 0.116]], 0.0012, wall, w - wall);
      blk(s, p.ss, [wall, yb - 0.006, zb - 0.045], [w - wall, yb - 0.003, zb + 0.004], 0);
      // The drain at the low end and the grease cup on its bracket under it.
      const xcup = w - wall - 0.08;
      const zcup = zt + 0.0475;
      s.add(p.ss, cyl(s, 0.011, 0.02, "y", 16), [xcup, 0.052, zcup]);
      s.add(p.ss, lathe(s, [
        [0, 0.004],
        [0.034, 0.004],
        [0.037, 0.008],
        [0.039, 0.05],
        [0.0405, 0.052],
        [0.0395, 0.0535],
        [0.0375, 0.05],
        [0.0355, 0.0095],
        [0, 0.0095],
      ], 32), [xcup, 0, zcup]);
      if (s.fine(1)) {
        s.add(p.ssBright, torus(s, 0.0425, 0.002, "y", Math.PI * 2, 32, 5), [xcup, 0.04, zcup]);
        for (const sz of [-1, 1]) s.add(p.ssBright, pipe(s, [[xcup, 0.04, zcup + sz * 0.0425], [xcup, 0.062, zcup + sz * 0.0425]], 0.0022, 0.004, 6));
        if (s.fine(2)) s.add(p.ssBright, pipe(s, [[xcup - 0.04, 0.03, zcup], [xcup - 0.042, 0.06, zcup], [xcup, 0.075, zcup], [xcup + 0.042, 0.06, zcup], [xcup + 0.04, 0.03, zcup]], 0.0014, 0.01, 6));
      }

      // Fire suppression: the 1/2 in branch line with its nozzles aimed down at the line, hung from the ceiling.
      const xa = 0.16;
      const xb = w - 0.16;
      s.add(p.chrome, pipe(s, [[xb, yf, zf], [xa, yf, zf], [xa, h - 0.006, zf]], 0.0107, 0.035, 12));
      if (s.fine(1)) {
        s.add(p.ssBright, hex(0.017, 0.024), [xa, h - 0.034, zf]);
        s.add(p.ssBright, hex(0.02, 0.012), [xa, hb + 0.006, zf]);
        if (s.fine(2)) s.add(p.ssBright, hex(0.02, 0.01), [xa, hb - T - 0.005, zf]);
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
      void k;
    };
    side(0);
    if (island) s.at([w, 0, -d], [0, Math.PI, 0], () => side(1));

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

    return (): Hooks => ({ set: (name, v) => name === "lights" && setLamp(lamp, v) });
  };
}


/**
 * A Type II condensate hood, over a dishmachine or steamer: no grease filters, a smooth welded liner, and a perimeter
 * gutter inside the skirts that catches what condenses on the liner and drains it from one corner. An internal baffle
 * under the collar makes the vapour turn and shed its water; its low edge drips into the back gutter. One light, and
 * the exhaust collar on top. y = 0 is the skirts' lower edge.
 */
const condensateHood: Builder = (s, e) => {
  const p = palette();
  const { w, d, h } = e;
  const hb = h - 0.065;
  const T = 0.0012;
  const wall = 0.035;
  const gy = 0.02;
  const r = s.fine(1) ? 0.0012 : 0;
  const lamp = lampMaterial();
  // The skins: outer, the inner liner, and the return closing the double wall, hemmed below.
  blk(s, p.ss, [0, 0, -T], [w, hb, 0], r);
  blk(s, p.ss, [0, 0, -d], [w, hb, -d + T], r);
  for (const x of [0, w - T]) blk(s, p.ss, [x, 0, -d + T], [x + T, hb, -T], r);
  blk(s, p.ss, [wall, 0, -wall], [w - wall, hb - T, -wall + T], r);
  blk(s, p.ss, [wall, 0, -d + wall - T], [w - wall, hb - T, -d + wall], r);
  for (const x of [wall - T, w - wall]) blk(s, p.ss, [x, 0, -d + wall], [x + T, hb - T, -wall], r);
  sheetWithHoles(s, p.ss, 0, w, -d, 0, 0, T, [[wall, w - wall, -d + wall, -wall]]);
  if (s.fine(2)) {
    for (const z of [-0.0035, -d + 0.0035]) s.add(p.ss, cyl(s, 0.0035, w, "x", 12), [w / 2, 0.0035, z]);
    for (const x of [0.0035, w - 0.0035]) s.add(p.ss, cyl(s, 0.0035, d, "z", 12), [x, 0.0035, -d / 2]);
  }
  // The top, cut for the collar.
  const cx0 = w / 2 - 0.225;
  const cx1 = w / 2 + 0.225;
  const cz0 = -d + 0.06;
  const cz1 = -d + 0.34;
  sheetWithHoles(s, p.ss, 0, w, -d, 0, hb - T, hb, [[cx0, cx1, cz0, cz1]]);
  ductCollar(s, cx0, cx1, cz0, cz1, hb, h);

  // The perimeter gutter round the liner's foot: a U with a rolled lip, its floor 20 mm above the skirt hem.
  const n = s.seg(6, 2);
  const gut = (o: number): P2[] => [[o, 0.065], [o, gy + 0.008], ...arc(o - 0.008, gy + 0.008, 0.008, 0, -Math.PI / 2, n).slice(1), [o - 0.037, gy], ...arc(o - 0.037, gy + 0.008, 0.008, -Math.PI / 2, -Math.PI, n).slice(1), [o - 0.045, gy + 0.035], ...arc(o - 0.041, gy + 0.035, 0.004, Math.PI, 0, 3).slice(1)];
  // Profiles are (z, y) for the front and back runs, (x, y) for the ends; each runs inward from its liner.
  sweepX(s, p.ss, gut(-wall), 0.0012, wall, w - wall);
  sweepX(s, p.ss, gut(-wall).map(([z, y]) => [-d - z, y] as P2), 0.0012, wall, w - wall);
  sweepZ(s, p.ss, gut(0).map(([x, y]) => [wall - x, y] as P2), 0.0012, -d + wall, -wall);
  sweepZ(s, p.ss, gut(0).map(([x, y]) => [w - wall + x, y] as P2), 0.0012, -d + wall, -wall);
  // The drain: a 1 in half coupling under the low corner, a close nipple and its cap.
  const dx = w - wall - 0.022;
  const dz = -wall - 0.022;
  s.add(p.ss, cyl(s, 0.017, 0.012, "y", 20), [dx, gy - 0.006, dz]);
  s.add(p.ssDull, cyl(s, 0.0167, 0.008, "y", 20), [dx, 0.006, dz]);
  if (s.fine(1)) s.add(p.ssDull, hex(0.02, 0.006), [dx, 0.0035, dz]);

  // The baffle: a sloped plate under the collar, its upper edge hemmed, its low edge turned down into a drip lip over
  // the back gutter; hung from the top by straps at each end.
  const [za, ya] = [-d + 0.55, hb - 0.095];
  const [zb, yb] = [-d + 0.058, 0.16];
  const L = Math.hypot(za - zb, ya - yb);
  const ang = Math.atan2(ya - yb, za - zb);
  s.at([w / 2, yb, zb], [-ang, 0, 0], () => {
    blk(s, p.ss, [-w / 2 + wall + 0.05, -0.0008, 0], [w / 2 - wall - 0.05, 0.0008, L], r);
    s.add(p.ss, cyl(s, 0.004, w - 2 * wall - 0.1, "x", 12), [0, 0.003, L]);
    blk(s, p.ss, [-w / 2 + wall + 0.05, -0.018, -0.0008], [w / 2 - wall - 0.05, 0, 0.0008], 0.0004);
  });
  for (const x of [wall + 0.08, w - wall - 0.08]) {
    s.add(p.ss, pipe(s, [[x, hb - T, za - 0.03], [x, ya + 0.004, za - 0.03]], 0.0025, 0.01, 8));
    s.add(p.ss, pipe(s, [[x, hb - T, zb + 0.12], [x, yb + 0.12 * Math.tan(ang) + 0.004, zb + 0.12]], 0.0025, 0.01, 8));
  }

  // The light, its box on top, hanger clips, the label.
  hoodLight(s, lamp, w / 2, hb - T, -0.3);
  if (s.fine(1)) blk(s, p.ssDull, [w / 2 - 0.05, hb, -0.35], [w / 2 + 0.05, hb + 0.045, -0.25], 0.002);
  for (const x of [0.1, w - 0.1]) for (const z of [-0.09, -d + 0.09]) hangerClip(s, x, hb, z);
  dataPlate(s, [wall + 0.0002, hb * 0.55, -0.25], [0, Math.PI / 2, 0], 0.09, 0.055);
  return () => ({ set: (name, v) => name === "lights" && setLamp(lamp, v) });
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

/** The raised rolled rim along a sink's front, as (z, y) from the deck's front edge: up into a bead, down an apron. */
function sinkRim(s: Shop): P2[] {
  const n = s.seg(10, 3);
  return [[-0.034, 0], [-0.018, 0], [-0.018, 0.012], ...arc(-0.009, 0.012, 0.009, Math.PI, 0, n).slice(1), [0, -0.04], ...arc(-0.006, -0.04, 0.006, 0, -Math.PI / 2, 3).slice(1), [-0.018, -0.046]];
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
  // The strainer in the bowl floor.
  s.add(p.ssBright, cyl(s, 0.05, 0.0025, "y", f2 ? 32 : 20), [x, y + 0.0013, z]);
  s.add(p.enamel, cyl(s, 0.036, 0.001, "y", f2 ? 28 : 20), [x, y + 0.003, z]);
  if (f2) {
    for (const a of [0, Math.PI / 2]) s.add(p.ssBright, rbox(s, 0.072, 0.003, 0.006, 0.001), [x, y + 0.0035, z], [0, a, 0]);
    s.add(p.ssBright, cyl(s, 0.009, 0.006, "y", 16), [x, y + 0.005, z]);
  }
  // The valve body under it: cast, chrome-plated, a 2 in tailpiece below.
  const yv = y - 0.045;
  const body: [number, number][] = f2 ? [[0, -0.05], [0.026, -0.05], [0.03, -0.04], [0.034, -0.02], [0.036, 0.0], [0.034, 0.025], [0.028, 0.04], [0.026, 0.045], [0, 0.045]] : [[0, -0.05], [0.03, -0.042], [0.036, 0], [0.027, 0.042], [0, 0.045]];
  s.add(p.chrome, lathe(s, body, f2 ? 24 : 16), [x, yv, z]);
  s.add(p.chrome, cyl(s, 0.0245, yv - 0.05 - 0.4, "y", f2 ? 24 : 16, 0.0245, true), [x, (yv - 0.05 + 0.4) / 2, z]);
  // The overflow, from a slotted cover high on the bowl's back wall down into the body (behind the bowl: close views).
  const yo = y + 0.2;
  if (f2) {
    s.add(p.ssBright, rbox(s, 0.05, 0.03, 0.004, 0.003), [x, yo, zBack + 0.002]);
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
    // The deck over the bowl block, cut for the bowls.
    const deck = new THREE.Shape();
    rectPath(deck, bx1 - bx0, -0.03 - zs, 0, (bx0 + bx1) / 2, -(-0.03 + zs) / 2);
    for (const x of S.bowls) deck.holes.push(rectPath(new THREE.Path(), S.bw, S.bdz, r, x, -(zb0 + zb1) / 2, true) as THREE.Path);
    const dg = plate(s, deck, t, 0, s.fine(2) ? 16 : 8);
    dg.rotateX(-Math.PI / 2);
    s.add(p.ss, dg, [0, Y - t, 0]);
    // Drainboards: pitched 15 mm toward the bowls, their outer ends turned up into rolled rims.
    for (const [a, b] of S.boards) {
      const outer = Math.abs(a - bx0) < Math.abs(b - bx1) && b <= bx0 + 1e-6 ? a : b;
      const inner = outer === a ? b : a;
      const len = Math.abs(b - a);
      const rise = 0.015;
      const ang = Math.atan2(rise, len) * (outer > inner ? 1 : -1);
      const cz = (zs - 0.03) / 2;
      s.at([(a + b) / 2, Y + rise / 2 - t / 2, cz], [0, 0, ang], () => blk(s, p.ss, [-len / 2 - 0.002, -t / 2, zs - cz], [len / 2 + 0.002, t / 2, -0.03 - cz], 0));
      const sx = outer > inner ? 1 : -1;
      sweepZ(s, p.ss, sinkRim(s).map(([z, y]) => [sx * z, y + rise] as P2), t, zs, -0.03, outer, Y - t / 2);
      // A hat channel under the board, front to back, the gussets of its legs welded to it.
      if (s.fine(2)) sweepZ(s, p.ss, [[-0.035, 0], [-0.018, 0], [-0.018, -0.026], [0.018, -0.026], [0.018, 0], [0.035, 0]], 0.0012, zs + 0.01, -0.05, outer - sx * 0.06, Y - 0.004);
    }
    // The ends of the bowl block where no board continues it.
    for (const [x, sx] of [
      [bx0, -1],
      [bx1, 1],
    ] as const)
      if (!S.boards.some(([a, b]) => Math.abs((sx < 0 ? b : a) - x) < 1e-3)) sweepZ(s, p.ss, sinkRim(s).map(([z, y]) => [sx * z, y] as P2), t, zs, -0.03, x, Y - t / 2);
    // The front rim, full length.
    sweepX(s, p.ss, sinkRim(s), t, x0, x1, Y - t / 2, -0.001);
    // The backsplash: coved up from the deck, a 45° return to the wall at its top, its ends closed.
    const n = s.seg(6, 2);
    const splash: P2[] = [[zs + 0.03, 0], ...arc(zs + 0.008, 0.008, 0.008, -Math.PI / 2, -Math.PI, n).slice(1), [zs, sh - 0.03], [zs - 0.026, sh], [-d + 0.0008, sh], [-d + 0.0008, sh - 0.03]];
    sweepX(s, p.ss, splash, t, x0, x1, Y - t / 2);
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
      for (const z of [L.z0, L.z1]) gussetLeg(s, L.x, z, L.top);
      crossRail(s, [L.x, 0.2, L.z1], [L.x, 0.2, L.z0]);
    }
    dataPlate(s, [x1 - 0.09, Y - 0.021, -0.0001], [0, 0, 0], 0.06, 0.022);
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
 * A 16-gauge 304 work table: the top, its front edge a boxed bull-nose or a marine lip, sides and back turned down (or a
 * coved backsplash with a 45° return), hat channels under it front to back that the leg gussets weld to, 1-5/8 in legs
 * on bullet feet, and an undershelf with turned-down edges on shelf collars.
 */
function workTable(opts: { legs: 4 | 6; edge: "bull" | "marine"; splash: number; shelf: "galv" | "ss" }): Builder {
  return (s, e) => {
    const p = palette();
    const { w, d } = e;
    const Y = 0.914;
    const t = 0.0016;
    const zBack = opts.splash ? -d + 0.025 : -d;
    const legX = opts.legs === 6 ? [0.055, w / 2, w - 0.055] : [0.055, w - 0.055];
    if (s.detail === 0) {
      // Far off: top and edge as one slab, the splash, legs to the floor and the undershelf.
      blk(s, p.ss, [0, Y - 0.04, -d], [w, Y, 0], 0);
      if (opts.splash) blk(s, p.ss, [0, Y, -d], [w, Y + opts.splash, -d + 0.025], 0);
      for (const x of legX) for (const z of [-0.06, -d + 0.06]) gussetLeg(s, x, z, Y - 0.04);
      blk(s, opts.shelf === "galv" ? p.ssDull : p.ss, [0.032, 0.222, -d + 0.036], [w - 0.032, 0.254, -0.036], 0);
      return;
    }
    // The top sheet and its edges.
    blk(s, p.ss, [0.0, Y - t, zBack + (opts.splash ? 0.012 : 0.02)], [w, Y, -0.02], 0);
    sweepX(s, p.ss, opts.edge === "bull" ? bullNose(s) : marineEdge(s), t, 0, w, Y - t / 2, 0);
    const n = s.seg(6, 2);
    const down: P2[] = [[0.024, 0], ...arc(0.006, -0.006, 0.006, Math.PI / 2, Math.PI, n).slice(1), [0, -0.038], [0.012, -0.038]];
    for (const [x, sx] of [
      [0, 1],
      [w, -1],
    ] as const) {
      sweepZ(s, p.ss, down.map(([u, v]) => [sx * u, v] as P2), t, zBack + (opts.splash ? 0 : 0.012), -0.02, x + sx * (t / 2), Y - t / 2);
      // The corner closed where side meets front edge.
      const edge = opts.edge === "bull" ? bullNose(s) : marineEdge(s);
      const cap = new THREE.Shape(edge.map(([z, y]) => new THREE.Vector2(-z, y)));
      const g = plate(s, cap, t);
      g.rotateY(Math.PI / 2);
      s.add(p.ss, g, [x - (sx < 0 ? t : 0), Y - t / 2, 0]);
    }
    if (opts.splash) {
      const sh = opts.splash;
      const zs = -d + 0.025;
      sweepX(s, p.ss, [[zs + 0.03, 0], ...arc(zs + 0.008, 0.008, 0.008, -Math.PI / 2, -Math.PI, n).slice(1), [zs, sh - 0.025], [zs - 0.024, sh - 0.001], [-d + 0.0008, sh - 0.001], [-d + 0.0008, sh - 0.025]], t, 0, w, Y - t / 2);
      for (const x of [0, w - t]) {
        const cap = new THREE.Shape([[-zs, 0], [-zs, sh - 0.025], [d - 0.001, sh], [d, sh], [d, 0]].map(([u, v]) => new THREE.Vector2(u!, v!)));
        const g = plate(s, cap, t);
        g.rotateY(Math.PI / 2);
        s.add(p.ss, g, [x, Y - t, 0]);
      }
    } else sweepX(s, p.ss, down.map(([u, v]) => [-d + u, v] as P2), t, 0.002, w - 0.002, Y - t / 2);

    // Hat channels under the top, front to back, at every pair of legs and between.
    const hats = [...legX];
    for (let i = 0; i < legX.length - 1; i++) {
      const span = legX[i + 1]! - legX[i]!;
      const k = Math.round(span / 0.6);
      for (let j = 1; j < k; j++) hats.push(legX[i]! + (span * j) / k);
    }
    const hatY = Y - t;
    for (const x of hats) sweepZ(s, p.ss, [[-0.04, 0], [-0.02, 0], [-0.02, -0.03], [0.02, -0.03], [0.02, 0], [0.04, 0]], 0.0012, zBack + 0.03, -0.045, x, hatY - 0.0006);
    // Legs, gussets welded under the hat channels.
    const zf = -0.06;
    const zr = -d + 0.06;
    for (const x of legX) for (const z of [zf, zr]) gussetLeg(s, x, z, hatY - 0.031);
    // The undershelf: 18-gauge, edges turned down, on collars round the legs.
    const ys = 0.254;
    const sm = opts.shelf === "galv" ? p.ssDull : p.ss;
    const sx0 = 0.032;
    const sx1 = w - 0.032;
    const sz0 = -d + 0.036;
    const sz1 = -0.036;
    blk(s, sm, [sx0, ys - 0.0012, sz0], [sx1, ys, sz1], 0);
    const lip: P2[] = [[-0.02, 0], [0, 0], [0, -0.032], [-0.01, -0.032]];
    sweepX(s, sm, lip, 0.0012, sx0, sx1, ys - 0.0006, sz1);
    sweepX(s, sm, lip.map(([z, y]) => [-z, y] as P2), 0.0012, sx0, sx1, ys - 0.0006, sz0);
    sweepZ(s, sm, lip, 0.0012, sz0, sz1, sx0, ys - 0.0006);
    sweepZ(s, sm, lip.map(([x, y]) => [-x, y] as P2), 0.0012, sz0, sz1, sx1, ys - 0.0006);
    for (const x of legX)
      for (const z of [zf, zr]) {
        s.add(sm, cyl(s, LEG_R + 0.0035, 0.045, "y", 24), [x, ys - 0.018, z]);
        if (s.fine(2)) s.add(p.ssBright, hex(0.0042, 0.006, "z"), [x, ys - 0.025, z + LEG_R + 0.006]);
      }
    if (w > 1.5) sweepX(s, sm, [[-0.03, 0], [-0.015, 0], [-0.015, -0.022], [0.015, -0.022], [0.015, 0], [0.03, 0]], 0.0012, sx0 + 0.03, sx1 - 0.03, ys - 0.0018, -d / 2);
    dataPlate(s, [w - 0.0016 - 0.0002, Y - 0.02, -d / 2], [0, Math.PI / 2, 0], 0.06, 0.025);
  };
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
  // Deck and the skirt round the bowl, front and sides.
  const deck = new THREE.Shape();
  rectPath(deck, w, -0.03 - zs, 0, w / 2, -(-0.03 + zs) / 2);
  deck.holes.push(rectPath(new THREE.Path(), bw, bdz, r, w / 2, -(zb0 + zb1) / 2, true) as THREE.Path);
  const dg = plate(s, deck, t, 0, s.fine(2) ? 16 : 8);
  dg.rotateX(-Math.PI / 2);
  s.add(p.ss, dg, [0, Y - t, 0]);
  sweepX(s, p.ss, sinkRim(s), t, 0, w, Y - t / 2, -0.001);
  const skirt = 0.115;
  blk(s, p.ss, [0.0006, Y - skirt, zs], [0.0018, Y - 0.002, -0.0006], 0.0004);
  blk(s, p.ss, [w - 0.0018, Y - skirt, zs], [w - 0.0006, Y - 0.002, -0.0006], 0.0004);
  blk(s, p.ss, [0.0006, Y - skirt, -0.0012], [w - 0.0006, Y - 0.04, 0], 0.0004);
  blk(s, p.ss, [0.0006, Y - skirt - 0.004, -0.012], [w - 0.0006, Y - skirt, 0], 0.0004);
  blk(s, p.ss, [0.0006, Y - skirt - 0.004, zs], [0.012, Y - skirt, 0], 0.0004);
  blk(s, p.ss, [w - 0.012, Y - skirt - 0.004, zs], [w - 0.0006, Y - skirt, 0], 0.0004);
  // Backsplash to the top of the box, with its return to the wall.
  const sh = h - Y;
  const n = s.seg(6, 2);
  sweepX(s, p.ss, [[zs + 0.03, 0], ...arc(zs + 0.008, 0.008, 0.008, -Math.PI / 2, -Math.PI, n).slice(1), [zs, sh - 0.022], [zs - 0.02, sh - 0.001], [-d + 0.0008, sh - 0.001], [-d + 0.0008, sh - 0.025]], t, 0, w, Y - t / 2);
  blk(s, p.ss, [0.004, 0.0, -d], [w - 0.004, Y - t, -d + 0.0012], 0.0004);
  // Side splash guards: welded trapezoids rising from the rim to the splash.
  for (const x of [0.0006, w - 0.0018]) {
    const g = plate(s, new THREE.Shape([new THREE.Vector2(0.03, 0), new THREE.Vector2(0.03, 0.04), new THREE.Vector2(-zs, 0.21), new THREE.Vector2(-zs, 0)]), 0.0012);
    g.rotateY(Math.PI / 2);
    s.add(p.ss, g, [x, Y, 0]);
    s.add(p.ss, pipe(s, [[x + 0.0006, Y + 0.04, -0.03], [x + 0.0006, Y + 0.21, zs]], 0.003, 0.01, 8));
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
  // The running water, a joint scaled by the flow.
  const stream = s.joint("stream", [xf, yo, zb + 2 * R]);
  const sg = cyl(stream, 0.0055, yo - yf, "y", 12, 0.0065);
  sg.translate(0, -(yo - yf) / 2, 0);
  stream.add(p.water, sg);

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
    s.add(p.hdpe, rbox(s, x1 - x0 - 0.07, 0.045, 0.0012, 0.0005), [(x0 + x1) / 2, y0 - 0.02, zd + 0.045], [0.2, 0, 0]);
    blk(s, p.glassSmoke, [x0 + 0.06, y0 + 0.03, zd + 0.089], [x1 - 0.06, y0 + 0.07, zd + 0.0915], 0.001);
    s.add(p.ssBright, cyl(s, 0.006, 0.004, "z", 14), [(x0 + x1) / 2, y1 - 0.02, zd + 0.091]);
    if (s.fine(2)) s.add(p.enamel, rbox(s, 0.0015, 0.006, 0.001, 0.0003), [(x0 + x1) / 2, y1 - 0.02, zd + 0.0932]);
  }
  dataPlate(s, [w - 0.06, Y - 0.07, 0.0002], [0, 0, 0], 0.05, 0.022);
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

/** A 20 × 20 in peg rack, its floor at y, centred at (x, z), loaded with `plates` dinner plates standing in its pegs. */
function dishRack(s: Shop, x: number, y: number, z: number, plates: number): void {
  const p = palette();
  const S = 0.5;
  const H = 0.1;
  // The open-sided polypropylene frame: corner posts, rails top and bottom, a lattice floor and the peg rows.
  for (const sx of [-1, 1]) for (const sz of [-1, 1]) blk(s, p.poly, [x + sx * (S / 2 - 0.02), y, z + sz * (S / 2 - 0.02)], [x + (sx * S) / 2, y + H, z + (sz * S) / 2], 0.004);
  for (const yy of [y + 0.006, y + H - 0.006])
    for (const sz of [-1, 1]) {
      blk(s, p.poly, [x - S / 2, yy - 0.006, z + sz * (S / 2 - 0.006) - 0.006], [x + S / 2, yy + 0.006, z + sz * (S / 2 - 0.006) + 0.006], 0.003);
      blk(s, p.poly, [x + sz * (S / 2 - 0.006) - 0.006, yy - 0.006, z - S / 2], [x + sz * (S / 2 - 0.006) + 0.006, yy + 0.006, z + S / 2], 0.003);
    }
  const nl = s.fine(1) ? 9 : 4;
  for (let i = 1; i < nl; i++) {
    blk(s, p.poly, [x - S / 2, y, z - S / 2 + (S * i) / nl - 0.003], [x + S / 2, y + 0.008, z - S / 2 + (S * i) / nl + 0.003], 0);
    blk(s, p.poly, [x - S / 2 + (S * i) / nl - 0.003, y, z - S / 2], [x - S / 2 + (S * i) / nl + 0.003, y + 0.008, z + S / 2], 0);
  }
  // Pegs in two rows across, and the plates standing between them, facing front.
  const pitch = 0.045;
  const n = Math.min(plates, 10);
  if (s.fine(1))
    for (let i = 0; i <= 10; i++)
      for (const px of [-0.09, 0.09]) s.add(p.poly, cyl(s, 0.0035, 0.06, "y", 8, 0.002), [x + px, y + 0.038, z - S / 2 + 0.03 + i * pitch]);
  const plateG = lathe(s, [[0, 0], [0.085, 0], [0.09, 0.004], [0.118, 0.016], [0.124, 0.02], [0.12, 0.022], [0.088, 0.008], [0, 0.006]], 40);
  plateG.rotateX(Math.PI / 2);
  for (let i = 0; i < n; i++) s.add(p.ceramic, plateG.clone(), [x, y + 0.008 + 0.124, z - S / 2 + 0.03 + (i + 0.5) * pitch], [-0.12, 0, 0]);
}

/**
 * A high-temperature door-type dishmachine: the wash tank on legs with its control fascia and gauges, a chamber with
 * rotating upper and lower wash arms and rinse arms, the rack track, and the three-sided hood that lifts straight up
 * on guides in the rear column to load a rack from either side. The rinse line, its vacuum breaker and pressure gauge
 * climb the back.
 */
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
  blk(s, p.ss, [0, legH, -d], [w, legH + t, 0], 0);
  for (const x of [0, w - t]) blk(s, p.ss, [x, legH, -d], [x + t, deck, -t], r);
  blk(s, p.ss, [t, legH, -d], [w - t, deck, -d + t], r);
  // Front: the lower access panel, louvred, and the control fascia above it.
  blk(s, p.ss, [0, legH, -t], [w, 0.6, 0], r);
  const nv = s.fine(1) ? 8 : 0;
  for (let i = 0; i < nv; i++) blk(s, p.enamel, [0.12 + i * 0.05, legH + 0.03, 0.0001], [0.12 + i * 0.05 + 0.008, legH + 0.1, 0.0006], 0);
  if (s.fine(2)) for (const x of [0.03, w - 0.03]) for (const y of [legH + 0.03, 0.57]) screw(s, p.ssBright, [x, y, 0], 0.003);
  blk(s, p.ss, [0, 0.6, -0.03], [w, deck, -0.028], r);
  blk(s, p.ss, [0, 0.6, -0.03], [w, 0.602, 0], 0);
  blk(s, p.ss, [0, deck - 0.002, -0.03], [w, deck, 0], 0);
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
  blk(s, p.ss, [0, deck - t, -d + rearD], [w, deck, -d + rearD + 0.04], 0);
  blk(s, p.ss, [0, deck - t, -0.06], [w, deck, -0.03], 0);
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
  blk(hood, p.ss, H([0.004, hy0, -0.0015]), H([w - 0.004, hy1, 0]), r);
  blk(hood, p.ss, H([0.004 + wt, hy0, -wt]), H([w - 0.004 - wt, hy1 - wt, -wt + 0.0012]), r);
  for (const [x0, sx] of [
    [0.004, 1],
    [w - 0.004, -1],
  ] as const) {
    blk(hood, p.ss, H([x0, hy0, hz0]), H([x0 + sx * 0.0015, hy1, 0]), r);
    blk(hood, p.ss, H([x0 + sx * (wt - 0.0012), hy0, hz0]), H([x0 + sx * wt, hy1 - wt, -wt]), r);
    blk(hood, p.ss, H([x0, hy0, hz0]), H([x0 + sx * wt, hy0 + 0.0015, 0]), 0);
  }
  blk(hood, p.ss, H([0.004, hy0, -wt]), H([w - 0.004, hy0 + 0.0015, 0]), 0);
  blk(hood, p.ss, H([0.004, hy1 - 0.0015, hz0]), H([w - 0.004, hy1, 0]), r);
  blk(hood, p.ss, H([0.004 + wt, hy1 - wt, hz0]), H([w - 0.004 - wt, hy1 - wt + 0.0012, -wt]), r);
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
  blk(s, p.ss, [0, tank0, zf - 0.01], [w, tank0 + t, zb], 0);
  blk(s, p.ss, [0, tank0, zf - 0.0015], [w, deck, zf], r);
  blk(s, p.ss, [0, tank0, zb - t], [w, deck, zb], r);
  for (const x of [0, w - t]) blk(s, p.ss, [x, tank0, zb], [x + t, deck, zf], r);
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
  blk(s, p.ss, [0, roof - t, zb], [w, roof, zf], r);
  for (const x of [0, w - t]) {
    blk(s, p.ss, [x, deck, zb], [x + t, o[0], zf], 0);
    blk(s, p.ss, [x, o[1], zb], [x + t, roof, zf], 0);
    blk(s, p.ss, [x, o[0], zb], [x + t, o[1], o[2]], 0);
    blk(s, p.ss, [x, o[0], o[3]], [x + t, o[1], zf], 0);
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
    const n = 6;
    const sw = (o[3] - o[2]) / n;
    for (let k = 0; k < n; k++) c.add(p.gasket, rbox(c, 0.003, o[1] - o[0] - 0.06, sw - 0.004, 0.001), [0, -(o[1] - o[0] - 0.06) / 2, -(o[3] - o[2]) / 2 + sw * (k + 0.5)]);
  });

  // Vent cowls over each end, their dampers' handles, and the control box on the roof.
  for (const xc of [0.16, w - 0.16]) {
    const g = cutAbove(new THREE.CylinderGeometry(0.11, 0.21, 0.09, 4, 1, false), "y", 0.0449);
    g.rotateY(Math.PI / 4);
    g.scale(1, 1, 0.55);
    s.add(p.ss, g, [xc, roof + 0.045, -d * 0.45]);
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
    peMat = new THREE.MeshStandardMaterial({ color: 0x2d2f33, roughness: 0.62, metalness: 0 });
    peMat.userData.uv = "box";
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
  frame.holes.push(rectPath(new THREE.Path(), cw + 2 * bev, ch + 2 * bev, rr + bev, cx - bx, cy - by, true) as THREE.Path);
  s.add(pe, plate(s, frame, R, bev, 16), [bx, by, zF - R]);
  // The inner skin: the cavity, its corners radiused, open to the front frame.
  const back = -d + 0.034;
  const Lc = zF - R - back;
  s.add(pe, inward(openFrontBox(s, cw, ch, Lc + rr, rr)), [cx, cy, back + (Lc + rr) / 2]);
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
        open: (f) => slideJoint(pj, "y", Math.max(0, Math.min(1, f)) * 0.09, rest),
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
  // The lid, its feed tube and the pusher on its own joint.
  s.add(p.glass, lathe(s, [[0.0, yt + 0.024], [0.07, yt + 0.02], [0.094, yt + 0.008], [0.1, yt], [0.1, yt - 0.006], [0.097, yt - 0.006], [0.097, yt + 0.002], [0.09, yt + 0.006], [0.068, yt + 0.016], [0, yt + 0.02]], 48), [cx, 0, cz]);
  const fx = cx - 0.035;
  const fz = cz + 0.03;
  s.add(p.glass, lathe(s, [[0.03, yt + 0.01], [0.031, 0.385], [0.036, 0.39], [0.036, 0.393], [0.027, 0.393], [0.026, yt + 0.012]], 32), [fx, 0, fz]);
  s.add(p.glass, rbox(s, 0.016, 0.13, 0.012, 0.004), [cx + 0.07, yt - 0.05, cz - 0.07], [0, Math.PI / 4, 0]);
  const pusher = s.joint("pusher", [fx, 0.393, fz]);
  pusher.add(p.hdpe, cyl(pusher, 0.024, 0.075, "y", 24), [0, -0.035, 0]);
  pusher.add(p.hdpe, lathe(pusher, [[0, 0], [0.036, 0], [0.037, 0.006], [0.03, 0.014], [0, 0.016]], 28), [0, 0.0005, 0]);
  return (b) => {
    let a = 0;
    const pj = b.joints.get("pusher");
    const rest: V3 = pj ? [pj.position.x, pj.position.y, pj.position.z] : [fx, 0.393, fz];
    return {
      tick: (_t, dt, load) => {
        a += dt * 45 * Math.max(0, Math.min(1, load));
        turnJoint(b.joints.get("s_blade_spin"), "y", a);
      },
      open: (f) => slideJoint(pj, "y", Math.max(0, Math.min(1, f)) * 0.09, rest),
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
  const tilt = 0.22;
  const bc: V3 = [0.29, 0.29, -0.33];
  const f2 = s.fine(2);
  // Rounded boxes at full detail; plain ones from further off, where 1 mm rounds are lost.
  const rb = (sh: Shop, x: number, y: number, z: number, r: number) => (f2 ? rbox(sh, x, y, z, r) : new THREE.BoxGeometry(x, y, z));
  const piv: V3 = [0.42, 0.2, -0.24];
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
  if (s.detail === 0) {
    // Far off: base and column, the knife disc and its guard ring, gauge plate, carriage tray and handle, the knob.
    blk(s, p.alu, [0.02, 0, -d + 0.02], [w - 0.02, 0.1, -0.02], 0);
    blk(s, p.alu, [0.065, 0.1, -0.58], [0.262, 0.41, -0.14], 0);
    s.at(bc, [0, 0, tilt], () => {
      s.joint("rotary_blade_spin", [0, 0, 0]).add(p.ssBright, alongAxis(new THREE.CylinderGeometry(R + 0.009, R + 0.009, 0.008, 14), "x"), [0.002, 0, 0]);
      s.joint("gauge_plate", [0, 0, 0]).add(p.alu, new THREE.BoxGeometry(0.012, 0.14, 0.19), [-0.011, -0.07, 0.095]);
    });
    const car = s.joint("carriage_slide", piv);
    car.at([0.01, 0, 0], [0, 0, tilt], () => car.add(p.alu, new THREE.BoxGeometry(0.21, 0.03, 0.26), [0, 0, 0]));
    car.add(p.alu, new THREE.BoxGeometry(0.02, 0.13, 0.02), [0.06, 0.07, 0.13]);
    s.joint("index_knob", [0.12, 0.065, -0.02]).add(p.bakelite, new THREE.CylinderGeometry(0.03, 0.03, 0.025, 8).rotateX(Math.PI / 2), [0, 0, 0.0125]);
    return slicerHooks;
  }
  // Feet and the cast base, its receiving area to the front left.
  for (const x of [0.06, w - 0.06]) for (const z of [-0.06, -d + 0.06]) s.add(p.rubber, lathe(s, f2 ? [[0, 0], [0.026, 0], [0.024, 0.02], [0.018, 0.03], [0, 0.03]] : [[0, 0], [0.026, 0.002], [0.02, 0.03], [0, 0.03]], 20), [x, 0, z]);
  blk(s, p.alu, [0.02, 0.03, -d + 0.02], [w - 0.02, 0.1, -0.02], 0.03);
  // The cast column behind the knife: one casting, its back rounded over, its face flat to the knife's hub.
  const col = new THREE.Shape();
  col.moveTo(0.07, 0.1);
  col.lineTo(0.262, 0.1);
  col.lineTo(0.262, 0.3);
  col.quadraticCurveTo(0.262, 0.4, 0.2, 0.41);
  col.lineTo(0.13, 0.4);
  col.quadraticCurveTo(0.07, 0.39, 0.065, 0.3);
  col.lineTo(0.07, 0.1);
  const cg = plate(s, col, 0.44, f2 ? 0.012 : 0, f2 ? 16 : 8);
  s.add(p.alu, cg, [0, 0, -0.58]);
  // The motor's end bell on the column's outer side, its vents and the cord gland.
  s.add(p.alu, lathe(s, [[0, 0], [0.085, 0], [0.085, 0.012], [0.07, 0.03], [0.03, 0.04], [0, 0.04]], 40).rotateZ(Math.PI / 2), [0.062, 0.27, -0.4]);
  if (f2) for (let i = 0; i < 7; i++) s.add(p.enamel, rb(s, 0.002, 0.006, 0.07, 0.002), [0.052, 0.24 + i * 0.01, -0.4]);
  s.add(p.rubber, cyl(s, 0.009, 0.03, "z", 12), [0.12, 0.13, -0.595]);
  blk(s, p.ss, [0.04, 0.1, -0.15], [0.25, 0.104, -0.03], 0.002);
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
    s.add(p.alu, cyl(s, 0.05, 0.03, "x", 32), [-0.02, 0, 0]);
    const knife = s.joint("rotary_blade_spin", [0, 0, 0]);
    const kp: [number, number][] = f2 ? [[0, -0.003], [0.11, -0.003], [R - 0.004, -0.0007], [R, 0], [R - 0.004, 0.0007], [0.11, 0.0022], [0.05, 0.004], [0, 0.004]] : [[0, -0.003], [R - 0.004, -0.001], [R, 0], [R - 0.004, 0.0012], [0, 0.004]];
    knife.add(p.ssBright, lathe(knife, kp, 72).rotateZ(-Math.PI / 2), [0.001, 0, 0]);
    if (knife.fine(1)) for (let i = 0; i < 4; i++) {
      const a = (i / 4) * Math.PI * 2;
      knife.add(p.ssDull, cyl(knife, 0.006, 0.003, "x", 12), [0.004, Math.cos(a) * 0.085, Math.sin(a) * 0.085]);
    }
    // The centre plate over the knife's middle, its removal knob.
    s.add(p.alu, lathe(s, [[0, 0.016], [0.06, 0.013], [0.1, 0.006], [0.116, 0.0], [0, 0]], 48).rotateZ(-Math.PI / 2), [0.007, 0, 0]);
    s.add(p.bakelite, lathe(s, [[0, 0], [0.022, 0], [0.024, 0.008], [0.016, 0.022], [0, 0.024]], 24).rotateZ(-Math.PI / 2), [0.022, 0, 0]);
    // The ring guard round the edge from the front top over to the back, open where the carriage passes.
    const a0 = -0.35;
    const a1 = Math.PI + 0.9;
    const ring = torus(s, R + 0.009, 0.007, "x", a1 - a0, 64, 8);
    ring.rotateX(-a0);
    s.add(p.alu, ring, [0.004, 0, 0]);
    s.add(p.alu, rb(s, 0.022, 0.05, 0.05, 0.008), [-0.012, Math.sin(Math.PI / 2) * (R + 0.012), 0]);
    // The gauge plate before the knife, sliding on its axis as the knob sets the slice.
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
    // The sharpener stowed on top: a cast housing with its two stones and knob.
    s.at([0.01, R + 0.035, -0.02], [0, 0, 0], () => {
      s.add(p.alu, rb(s, 0.07, 0.04, 0.09, 0.012), [0, 0, 0]);
      for (const z of [-0.025, 0.025]) s.add(p.ceramic, cyl(s, 0.014, 0.012, "x", 16), [0.03, -0.012, z]);
      s.add(p.bakelite, knobGeo(s, 0.016), [0, 0.02, 0]);
    });
  });

  // The carriage: an arm from the slide under the base's right, the tray tilted to the knife, its handle and pusher.
  s.add(p.ssBright, cyl(s, 0.012, d - 0.12, "z", 16), [0.45, 0.115, -d / 2]);
  for (const z of [-0.07, -d + 0.07]) s.add(p.alu, rb(s, 0.04, 0.03, 0.03, 0.006), [0.45, 0.112, z]);
  const car = s.joint("carriage_slide", piv);
  const C = (q: V3): V3 => [q[0] - piv[0], q[1] - piv[1], q[2] - piv[2]];
  car.add(p.alu, rb(car, 0.05, 0.04, 0.08, 0.008), C([0.45, 0.13, -0.24]));
  car.add(p.alu, rb(car, 0.03, 0.06, 0.05, 0.008), C([0.43, 0.165, -0.24]), [0, 0, tilt]);
  car.at(C([0.43, 0.2, -0.24]), [0, 0, tilt], () => {
    car.add(p.alu, rb(car, 0.21, 0.008, 0.26, 0.003), [0.0, 0, 0]);
    for (const sz of [-1, 1]) car.add(p.alu, rb(car, 0.21, 0.03, 0.006, 0.002), [0, 0.014, sz * 0.128]);
    car.add(p.alu, rb(car, 0.006, 0.05, 0.26, 0.002), [0.102, 0.022, 0]);
    // The pusher: a sleeve on a rod over the tray, a spiked plate that bears the product toward the knife, a grip.
    car.add(p.ssBright, cyl(car, 0.006, 0.2, "x", 12), [0.0, 0.13, -0.11]);
    car.add(p.alu, rb(car, 0.012, 0.12, 0.012, 0.004), [0.1, 0.065, -0.11]);
    car.add(p.alu, rb(car, 0.04, 0.025, 0.025, 0.006), [0.02, 0.13, -0.11]);
    car.add(p.alu, rb(car, 0.008, 0.08, 0.12, 0.003), [0.0, 0.08, -0.05]);
    car.add(p.alu, rb(car, 0.02, 0.012, 0.07, 0.004), [0.015, 0.125, -0.08]);
    car.add(p.bakelite, knobGeo(car, 0.016), [0.02, 0.15, -0.11]);
    if (car.fine(2)) for (let i = 0; i < 4; i++) for (let k = 0; k < 3; k++) car.add(p.ssBright, cyl(car, 0.0016, 0.01, "x", 6, 0.0003), [-0.008, 0.055 + k * 0.02, -0.09 + i * 0.026]);
  });
  // The handle at the carriage's front right.
  car.add(p.alu, pipe(car, [C([0.48, 0.21, -0.12]), C([0.48, 0.33, -0.11])], 0.008, 0.02, 12));
  car.add(p.bakelite, lathe(car, [[0, 0], [0.016, 0], [0.019, 0.02], [0.016, 0.04], [0, 0.043]], 20), C([0.48, 0.325, -0.11]));

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
  // The platter: a pressed sheet with its edges turned down.
  blk(s, p.ss, [w / 2 - pw / 2, h - 0.0012, pz - pw / 2], [w / 2 + pw / 2, h, pz + pw / 2], 0.0005);
  const lip: P2[] = [[-0.006, 0], ...arc(-0.004, -0.002, 0.002, Math.PI / 2, 0, 3), [-0.002, -0.01]];
  sweepX(s, p.ss, lip.map(([z, y]) => [z, y] as P2), 0.0012, w / 2 - pw / 2, w / 2 + pw / 2, h - 0.0006, pz + pw / 2);
  sweepX(s, p.ss, lip.map(([z, y]) => [-z, y] as P2), 0.0012, w / 2 - pw / 2, w / 2 + pw / 2, h - 0.0006, pz - pw / 2);
  sweepZ(s, p.ss, lip.map(([x, y]) => [-x, y] as P2), 0.0012, pz - pw / 2, pz + pw / 2, w / 2 + pw / 2, h - 0.0006);
  sweepZ(s, p.ss, lip.map(([x, y]) => [x, y] as P2), 0.0012, pz - pw / 2, pz + pw / 2, w / 2 - pw / 2, h - 0.0006);
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
    // The lip: punched outward along the hole's lower edge, so food pushed down the face meets the cut.
    xb.fillStyle = "#f2f2f2";
    xb.beginPath();
    xb.ellipse(cx, cy + ry * 0.55, rx * 1.2, ry * 1.1, 0, 0, Math.PI);
    xb.fill();
    xb.fillStyle = "#3a3a3a";
    xb.beginPath();
    xb.ellipse(cx, cy - ry * 0.2, rx * 1.05, ry * 0.7, 0, Math.PI, Math.PI * 2);
    xb.fill();
  };
  const margin = 18;
  // The face tapers to GRATER_TAPER of its foot's width at the top; the pattern keeps its pitch and loses holes instead.
  const half = (cy: number) => ((GRATER_TAPER + (1 - GRATER_TAPER) * (cy / H)) * W) / 2 - margin;
  if (kind === "slicer") {
    for (const cy of [H * 0.38, H * 0.62]) {
      const hw = half(cy);
      xa.fillStyle = "#000";
      xa.fillRect(W / 2 - hw + 20, cy - 10, 2 * hw - 40, 20);
      xb.fillStyle = "#f0f0f0";
      xb.fillRect(W / 2 - hw + 16, cy + 10, 2 * hw - 32, 8);
    }
  } else {
    const [cols, rows, rx, ry] = kind === "coarse" ? [5, 11, 13, 9] : kind === "medium" ? [8, 18, 7, 5] : [16, 36, 2.6, 2.6];
    for (let j = 0; j < rows; j++)
      for (let i = 0; i < cols - (j % 2); i++) {
        const cx = margin + ((i + 0.5 + (j % 2) * 0.5) * (W - 2 * margin)) / cols;
        const cy = margin + ((j + 0.5) * (H - 2 * margin)) / rows;
        if (Math.abs(cx - W / 2) + rx * 1.2 > half(cy)) continue;
        hole(cx, cy, rx, ry);
        if (kind === "fine") {
          xb.fillStyle = "#f2f2f2";
          xb.fillRect(cx - 4, cy + 3, 8, 3);
        }
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
  for (const [k, a, b, c, dd] of faces) s.add(s.fine(1) ? graterFace(k) : p.ss, quad(a, b, c, dd));
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
  flat(roundRect(w - 0.003, d - 0.003, 0.012), h - g0, 0.0015, 0);
  if (!s.fine(1)) return;
  flat(roundRect(w - 0.075, d - 0.075, 0.02), g0 + 0.0002, 0.0012, h - g0 - 0.0002);
  const rim = roundRect(w - 0.003, d - 0.003, 0.012);
  rim.holes.push(rectPath(new THREE.Path(), w - 0.057, d - 0.057, 0.028, 0, 0, true) as THREE.Path);
  flat(rim, g0 + 0.0002, 0.0012, h - g0 - 0.0002);
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
