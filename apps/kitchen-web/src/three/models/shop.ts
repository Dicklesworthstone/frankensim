/**
 * The workshop every kitchen model is built in: solid primitives placed in a fixture's own frame, gathered by
 * material and merged, so that a whole kitchen draws in a few dozen calls however many screws it has.
 *
 * Fixture frame (metres): x runs across the width from the left end (0) to the right (w) as seen from the front;
 * y is up from the floor or the fixture's base (0) to its top (h); z runs from the back (-d) to the front face (0),
 * which faces +z. A manifest pose maps onto it exactly: position (x_m, z_m, -y_m), turned about +y by yaw_rad
 * (see `placeFixture`). Moving sub-assemblies (doors, baskets, fans, lids) are `joint`s: groups pivoting about a
 * point in the frame, their own parts merged separately so that they can turn.
 *
 * Detail: 0 is a whole floor seen from across the room, 1 a kitchen seen from its aisle, 2 a close-up. Segment
 * counts and the small hardware (screws, seams, perforations) scale with it.
 */
import * as THREE from "three";
import { RoundedBoxGeometry } from "three/examples/jsm/geometries/RoundedBoxGeometry.js";
import { mergeGeometries } from "three/examples/jsm/utils/BufferGeometryUtils.js";

export type V3 = readonly [number, number, number];
export type Detail = 0 | 1 | 2;
export type Axis = "x" | "y" | "z";

const ID = new THREE.Matrix4();
const tmpM = new THREE.Matrix4();
const tmpQ = new THREE.Quaternion();
const tmpE = new THREE.Euler();
const tmpV = new THREE.Vector3();
const ONE = new THREE.Vector3(1, 1, 1);

/** Turn a solid built about +y so that its axis runs along `axis`. */
export function alongAxis(g: THREE.BufferGeometry, axis: Axis): THREE.BufferGeometry {
  if (axis === "x") g.rotateZ(-Math.PI / 2);
  else if (axis === "z") g.rotateX(Math.PI / 2);
  return g;
}

/** A float attribute stored one value per component, `size` to an item: the loops below read and write its array. */
function plain(a: THREE.BufferAttribute | THREE.InterleavedBufferAttribute | undefined, size: number): a is THREE.BufferAttribute {
  return a !== undefined && !(a as THREE.InterleavedBufferAttribute).isInterleavedBufferAttribute && !a.normalized && a.itemSize === size && a.array instanceof Float32Array;
}

const normalM = new THREE.Matrix3();

/**
 * An attribute's first `size` components, item after item, as three.js's getX/getY/getZ/getW give them. A plain
 * array is returned as it is, the attribute's own, so read it and do not write to it; anything else (interleaved,
 * normalised integers) is copied out with three.js's arithmetic: a normalised integer over its type's largest value,
 * no lower than -1 when signed. Half floats go through three.js's accessors.
 */
export function attributeValues(a: THREE.BufferAttribute | THREE.InterleavedBufferAttribute, size: number): ArrayLike<number> {
  const inter = (a as THREE.InterleavedBufferAttribute).isInterleavedBufferAttribute === true;
  // three.js sets this flag on Float16BufferAttribute; its type declarations omit it.
  const half = (a as { isFloat16BufferAttribute?: boolean }).isFloat16BufferAttribute === true;
  if (!inter && !half && !a.normalized && a.itemSize === size) return a.array;
  const n = a.count;
  const out = new Float64Array(n * size);
  if (half) {
    const get = [a.getX, a.getY, a.getZ, a.getW];
    for (let i = 0; i < n; i++) for (let c = 0; c < size; c++) out[i * size + c] = get[c]!.call(a, i);
    return out;
  }
  const arr = a.array;
  const stride = inter ? (a as THREE.InterleavedBufferAttribute).data.stride : a.itemSize;
  const off = inter ? (a as THREE.InterleavedBufferAttribute).offset : 0;
  // three.js's denormalize, by the array's type: 0 leaves the value as it is.
  let div = 0;
  let signed = false;
  if (a.normalized) {
    const C = arr.constructor;
    if (C === Uint32Array) div = 4294967295.0;
    else if (C === Uint16Array) div = 65535.0;
    else if (C === Uint8Array) div = 255.0;
    else if (C === Int32Array) [div, signed] = [2147483647.0, true];
    else if (C === Int16Array) [div, signed] = [32767.0, true];
    else if (C === Int8Array) [div, signed] = [127.0, true];
    else if (C !== Float32Array) throw new Error("Invalid component type.");
  }
  for (let i = 0, o = 0; i < n; i++) {
    const b = i * stride + off;
    for (let c = 0; c < size; c++, o++) {
      const x = arr[b + c]!;
      out[o] = div === 0 ? x : signed ? Math.max(x / div, -1.0) : x / div;
    }
  }
  return out;
}

/** Points `src` (x, y, z after each other) moved by `m` into `dst`, as Vector3.applyMatrix4 does; `dst` may be `src`. */
export function applyMatrixInto(src: ArrayLike<number>, dst: Float32Array, m: THREE.Matrix4): void {
  const e = m.elements;
  for (let i = 0, l = dst.length; i < l; i += 3) {
    const x = src[i]!;
    const y = src[i + 1]!;
    const z = src[i + 2]!;
    const w = 1 / (e[3]! * x + e[7]! * y + e[11]! * z + e[15]!);
    dst[i] = (e[0]! * x + e[4]! * y + e[8]! * z + e[12]!) * w;
    dst[i + 1] = (e[1]! * x + e[5]! * y + e[9]! * z + e[13]!) * w;
    dst[i + 2] = (e[2]! * x + e[6]! * y + e[10]! * z + e[14]!) * w;
  }
}

/**
 * Normals `src` turned by `k3` (a move's normal matrix, from Matrix3.getNormalMatrix) and made unit length into `dst`,
 * as Vector3.applyMatrix3 then normalize do; `dst` may be `src`.
 */
export function applyNormalInto(src: ArrayLike<number>, dst: Float32Array, k3: THREE.Matrix3): void {
  const k = k3.elements;
  for (let i = 0, l = dst.length; i < l; i += 3) {
    const x = src[i]!;
    const y = src[i + 1]!;
    const z = src[i + 2]!;
    const nx = k[0]! * x + k[3]! * y + k[6]! * z;
    const ny = k[1]! * x + k[4]! * y + k[7]! * z;
    const nz = k[2]! * x + k[5]! * y + k[8]! * z;
    const s = 1 / (Math.sqrt(nx * nx + ny * ny + nz * nz) || 1);
    dst[i] = nx * s;
    dst[i + 1] = ny * s;
    dst[i + 2] = nz * s;
  }
}

/**
 * `g.applyMatrix4(m)` with three.js's per-vertex loops written out over the arrays: Vector3.applyMatrix4 for
 * positions, and applyMatrix3 then normalize for normals, in the same order, so the floats are the same. three.js
 * reads and writes every vertex through getX/getY/getZ and setXYZ instead. Tangents, and data that is not plain
 * floats, go through three.js.
 */
export function transformGeometry(g: THREE.BufferGeometry, m: THREE.Matrix4): THREE.BufferGeometry {
  const p = g.attributes.position;
  const n = g.attributes.normal;
  if (!plain(p, 3) || (n !== undefined && !plain(n, 3)) || g.attributes.tangent !== undefined) return g.applyMatrix4(m);
  const a = p.array as Float32Array;
  applyMatrixInto(a, a, m);
  p.needsUpdate = true;
  if (n !== undefined) {
    const b = n.array as Float32Array;
    applyNormalInto(b, b, normalM.getNormalMatrix(m));
    n.needsUpdate = true;
  }
  if (g.boundingBox !== null) g.computeBoundingBox();
  if (g.boundingSphere !== null) boundSphere(g);
  return g;
}

/**
 * `g.computeBoundingSphere()` over a plain float position array: the centre of the bounding box, then the farthest
 * vertex from it, with three.js's arithmetic, so the sphere is the same. Morphed or other data goes through three.js.
 */
export function boundSphere(g: THREE.BufferGeometry): void {
  const p = g.attributes.position;
  if (!plain(p, 3) || g.morphAttributes.position !== undefined) {
    g.computeBoundingSphere();
    return;
  }
  const a = p.array as Float32Array;
  const l = p.count * 3;
  let x0 = Infinity;
  let y0 = Infinity;
  let z0 = Infinity;
  let x1 = -Infinity;
  let y1 = -Infinity;
  let z1 = -Infinity;
  for (let i = 0; i < l; i += 3) {
    x0 = Math.min(x0, a[i]!);
    y0 = Math.min(y0, a[i + 1]!);
    z0 = Math.min(z0, a[i + 2]!);
    x1 = Math.max(x1, a[i]!);
    y1 = Math.max(y1, a[i + 1]!);
    z1 = Math.max(z1, a[i + 2]!);
  }
  const s = (g.boundingSphere ??= new THREE.Sphere());
  const empty = x1 < x0 || y1 < y0 || z1 < z0;
  const cx = empty ? 0 : (x0 + x1) * 0.5;
  const cy = empty ? 0 : (y0 + y1) * 0.5;
  const cz = empty ? 0 : (z0 + z1) * 0.5;
  s.center.set(cx, cy, cz);
  let r2 = 0;
  for (let i = 0; i < l; i += 3) {
    const dx = cx - a[i]!;
    const dy = cy - a[i + 1]!;
    const dz = cz - a[i + 2]!;
    r2 = Math.max(r2, dx * dx + dy * dy + dz * dz);
  }
  s.radius = Math.sqrt(r2);
  if (isNaN(s.radius)) console.error('THREE.BufferGeometry.computeBoundingSphere(): Computed radius is NaN. The "position" attribute is likely to have NaN values.', g);
}

/** Every geometry merged must carry the same attributes and an index. */
function normalise(g: THREE.BufferGeometry): THREE.BufferGeometry {
  for (const name of Object.keys(g.attributes)) if (name !== "position" && name !== "normal" && name !== "uv") g.deleteAttribute(name);
  if (!g.attributes.normal) g.computeVertexNormals();
  if (!g.attributes.uv) g.setAttribute("uv", new THREE.BufferAttribute(new Float32Array(g.attributes.position.count * 2), 2));
  if (!g.index) {
    const n = g.attributes.position.count;
    const idx = n > 65535 ? new Uint32Array(n) : new Uint16Array(n);
    for (let i = 0; i < n; i++) idx[i] = i;
    g.setIndex(new THREE.BufferAttribute(idx, 1));
  }
  g.clearGroups();
  return g;
}

/**
 * Box-projected UVs in metres, `scale` repeats per metre: each face takes the two axes its normal is not along,
 * so brushed grain runs level on every vertical face and front to back on every top.
 */
export function boxUVs(g: THREE.BufferGeometry, scale = 1): void {
  const p = g.attributes.position;
  const n = g.attributes.normal;
  const uv = g.attributes.uv as THREE.BufferAttribute;
  if (plain(p, 3) && plain(n, 3) && plain(uv, 2)) {
    // The accessor loop below, reading and writing the arrays directly: the same branches and products.
    const P = p.array as Float32Array;
    const N = n.array as Float32Array;
    const U = uv.array as Float32Array;
    for (let i = 0, c = p.count; i < c; i++) {
      const j = i * 3;
      const ax = Math.abs(N[j]!);
      const ay = Math.abs(N[j + 1]!);
      const az = Math.abs(N[j + 2]!);
      const k = i * 2;
      if (ay >= ax && ay >= az) {
        U[k] = P[j]! * scale;
        U[k + 1] = P[j + 2]! * scale;
      } else if (ax >= az) {
        U[k] = P[j + 2]! * scale;
        U[k + 1] = P[j + 1]! * scale;
      } else {
        U[k] = P[j]! * scale;
        U[k + 1] = P[j + 1]! * scale;
      }
    }
    uv.needsUpdate = true;
    return;
  }
  for (let i = 0; i < p.count; i++) {
    const [ax, ay, az] = [Math.abs(n.getX(i)), Math.abs(n.getY(i)), Math.abs(n.getZ(i))];
    if (ay >= ax && ay >= az) uv.setXY(i, p.getX(i) * scale, p.getZ(i) * scale);
    else if (ax >= az) uv.setXY(i, p.getZ(i) * scale, p.getY(i) * scale);
    else uv.setXY(i, p.getX(i) * scale, p.getY(i) * scale);
  }
  uv.needsUpdate = true;
}

/** How a material's parts are textured once merged: projected in metres, or keeping their own UVs. */
export function uvMode(m: THREE.Material): "box" | "own" {
  return (m.userData.uv as "box" | "own" | undefined) ?? "box";
}

interface Bin {
  mat: THREE.Material;
  geos: THREE.BufferGeometry[];
}

/** One rigid body's parts, gathered by material. */
class Bins {
  private map = new Map<THREE.Material, Bin>();
  add(mat: THREE.Material, g: THREE.BufferGeometry) {
    let bin = this.map.get(mat);
    if (!bin) this.map.set(mat, (bin = { mat, geos: [] }));
    bin.geos.push(g);
  }
  /** Merged meshes, one per material. */
  meshes(shadows: boolean): THREE.Mesh[] {
    const out: THREE.Mesh[] = [];
    for (const { mat, geos } of this.map.values()) {
      const merged = geos.length === 1 ? geos[0]! : mergeGeometries(geos, false);
      if (!merged) continue;
      if (geos.length > 1) for (const g of geos) g.dispose();
      boundSphere(merged);
      const mesh = new THREE.Mesh(merged, mat);
      const lit = !(mat as THREE.MeshBasicMaterial).isMeshBasicMaterial && !mat.transparent;
      mesh.castShadow = shadows && lit && mat.userData.castShadow !== false;
      mesh.receiveShadow = shadows && lit;
      if (mat.transparent) mesh.renderOrder = 2;
      out.push(mesh);
    }
    this.map.clear();
    return out;
  }
  get size() {
    return this.map.size;
  }
}

export interface JointSpec {
  /** Pivot in the parent's frame. */
  at: V3;
  /** Fixed turn of the joint frame (Euler XYZ, radians). */
  rot?: V3;
}

/**
 * A model under construction. Parts go into the current frame (a stack of placements); `joint` opens a child that
 * becomes its own group, pivoting about a point, for anything that moves.
 */
export class Shop {
  readonly detail: Detail;
  readonly name: string;
  private bins = new Bins();
  private stack: THREE.Matrix4[] = [ID.clone()];
  private children: { shop: Shop; spec: JointSpec }[] = [];
  /** Anything that must sit in the result as is (lights, instanced flames, sprites). */
  private extras: THREE.Object3D[] = [];
  /** Counts of solids, for the size readouts. */
  tris = 0;

  constructor(detail: Detail, name = "model") {
    this.detail = detail;
    this.name = name;
  }

  /** Segment count for a curve of `n` segments at full detail. */
  seg(n: number, min = 6): number {
    const k = this.detail === 2 ? 1 : this.detail === 1 ? 0.5 : 0.28;
    return Math.max(min, Math.round(n * k));
  }
  /** True when parts this small are worth drawing. */
  fine(level: Detail = 2): boolean {
    return this.detail >= level;
  }

  /** Place the following parts at `at`, turned by `rot`, until the matching `pop`. */
  push(at: V3, rot: V3 = [0, 0, 0]): this {
    tmpQ.setFromEuler(tmpE.set(rot[0], rot[1], rot[2]));
    const m = new THREE.Matrix4().compose(tmpV.set(at[0], at[1], at[2]), tmpQ, ONE);
    this.stack.push(this.stack[this.stack.length - 1]!.clone().multiply(m));
    return this;
  }
  pop(): this {
    if (this.stack.length > 1) this.stack.pop();
    return this;
  }
  /** Run `f` with the parts it adds placed at `at`, turned by `rot`. */
  at(at: V3, rot: V3, f: () => void): this {
    this.push(at, rot);
    try {
      f();
    } finally {
      this.pop();
    }
    return this;
  }

  /** Add a solid, built about its own origin, at `at` turned by `rot` (Euler XYZ) within the current frame. */
  add(mat: THREE.Material, g: THREE.BufferGeometry, at: V3 = [0, 0, 0], rot?: V3, scale?: V3): this {
    tmpQ.setFromEuler(tmpE.set(rot?.[0] ?? 0, rot?.[1] ?? 0, rot?.[2] ?? 0));
    tmpM.compose(tmpV.set(at[0], at[1], at[2]), tmpQ, scale ? new THREE.Vector3(scale[0], scale[1], scale[2]) : ONE);
    transformGeometry(g, tmpM.premultiply(this.stack[this.stack.length - 1]!));
    normalise(g);
    if (uvMode(mat) === "box") boxUVs(g, (mat.userData.uvScale as number | undefined) ?? 1);
    this.tris += (g.index ? g.index.count : g.attributes.position.count) / 3;
    this.bins.add(mat, g);
    return this;
  }

  /** An object kept as it is (lights, instanced meshes), placed in the current frame. */
  extra(o: THREE.Object3D, at: V3 = [0, 0, 0], rot?: V3): this {
    const holder = new THREE.Group();
    holder.matrixAutoUpdate = false;
    tmpQ.setFromEuler(tmpE.set(rot?.[0] ?? 0, rot?.[1] ?? 0, rot?.[2] ?? 0));
    holder.matrix.compose(tmpV.set(at[0], at[1], at[2]), tmpQ, ONE).premultiply(this.stack[this.stack.length - 1]!);
    holder.add(o);
    this.extras.push(holder);
    return this;
  }

  /** A moving sub-assembly pivoting about `at` (in the current frame). Its parts are placed relative to the pivot. */
  joint(name: string, at: V3, rot: V3 = [0, 0, 0]): Shop {
    const shop = new Shop(this.detail, name);
    // The joint's pivot, carried through the current frame.
    const p = new THREE.Vector3(at[0], at[1], at[2]).applyMatrix4(this.stack[this.stack.length - 1]!);
    const q = new THREE.Quaternion().setFromRotationMatrix(this.stack[this.stack.length - 1]!).multiply(new THREE.Quaternion().setFromEuler(new THREE.Euler(rot[0], rot[1], rot[2])));
    const e = new THREE.Euler().setFromQuaternion(q);
    this.children.push({ shop, spec: { at: [p.x, p.y, p.z], rot: [e.x, e.y, e.z] } });
    return shop;
  }

  /** Merge everything into a group. Joints become named child groups; `joints` maps their names to them. */
  build(shadows = true): Built {
    const group = new THREE.Group();
    group.name = this.name;
    const joints = new Map<string, THREE.Group>();
    for (const m of this.bins.meshes(shadows)) group.add(m);
    for (const x of this.extras) group.add(x);
    let tris = this.tris;
    for (const { shop, spec } of this.children) {
      const child = shop.build(shadows);
      const pivot = child.group;
      pivot.position.set(spec.at[0], spec.at[1], spec.at[2]);
      pivot.rotation.set(spec.rot?.[0] ?? 0, spec.rot?.[1] ?? 0, spec.rot?.[2] ?? 0);
      pivot.userData.rest = pivot.quaternion.clone();
      group.add(pivot);
      joints.set(shop.name, pivot);
      for (const [k, v] of child.joints) joints.set(k, v);
      tris += child.tris;
    }
    return { group, joints, tris };
  }
}

export interface Built {
  group: THREE.Group;
  joints: Map<string, THREE.Group>;
  tris: number;
}

/** Turn a joint from its rest about its own local axis. */
export function turnJoint(j: THREE.Object3D | undefined, axis: Axis, angle: number): void {
  if (!j) return;
  const rest = j.userData.rest as THREE.Quaternion | undefined;
  tmpV.set(axis === "x" ? 1 : 0, axis === "y" ? 1 : 0, axis === "z" ? 1 : 0);
  tmpQ.setFromAxisAngle(tmpV, angle);
  if (rest) j.quaternion.copy(rest).multiply(tmpQ);
  else j.quaternion.copy(tmpQ);
}

/** Slide a joint from its rest along its own local axis. */
export function slideJoint(j: THREE.Object3D | undefined, axis: Axis, dist: number, restAt: V3): void {
  if (!j) return;
  tmpV.set(axis === "x" ? dist : 0, axis === "y" ? dist : 0, axis === "z" ? dist : 0).applyQuaternion(j.quaternion);
  j.position.set(restAt[0] + tmpV.x, restAt[1] + tmpV.y, restAt[2] + tmpV.z);
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Solids. Each is built about its own origin; boxes and cylinders are centred, unless the name says otherwise.   */
/* ------------------------------------------------------------------------------------------------------------ */

const rboxMade = new Map<string, THREE.BufferGeometry>();
/**
 * A box, its edges rounded by `r` when there is detail to show it. One segment (108 triangles) rounds a small part
 * convincingly with smooth normals; only big panels seen close up get two (300).
 */
export function rbox(shop: Shop, w: number, h: number, d: number, r = 0.0015): THREE.BufferGeometry {
  const rr = Math.min(r, w / 2 - 1e-5, h / 2 - 1e-5, d / 2 - 1e-5);
  if (shop.detail === 0 || rr <= 1e-4) return new THREE.BoxGeometry(w, h, d);
  const big = rr >= 0.004 && Math.max(w, h, d) > 0.2;
  const segs = shop.detail === 2 && big ? 2 : 1;
  // Most rounded boxes on the page repeat one already made (1,190 made of 332 sizes): each size is made once, and a
  // copy handed out, as the shop moves and merges what it is given.
  const key = `${w}|${h}|${d}|${segs}|${rr}`;
  let g = rboxMade.get(key);
  if (!g) rboxMade.set(key, (g = new RoundedBoxGeometry(w, h, d, segs, rr)));
  return new THREE.BufferGeometry().copy(g);
}

/**
 * A bar whose plan is a trapezoid: `len` long along x from its root (width `w0`) at x = 0 to its tip (width `w1`),
 * `h` tall from y = 0, its edges rounded by `r`. Grate fingers, brackets, ribs.
 */
export function taper(shop: Shop, len: number, w0: number, w1: number, h: number, r = 0.002): THREE.BufferGeometry {
  const s = new THREE.Shape([new THREE.Vector2(0, -w0 / 2), new THREE.Vector2(len, -w1 / 2), new THREE.Vector2(len, w1 / 2), new THREE.Vector2(0, w0 / 2)]);
  const bev = shop.detail > 0 ? Math.min(r, w1 / 2 - 1e-4, h / 2 - 1e-4) : 0;
  const g = new THREE.ExtrudeGeometry(s, {
    depth: Math.max(1e-5, h - 2 * bev),
    bevelEnabled: bev > 0,
    bevelThickness: bev,
    bevelSize: bev,
    bevelOffset: -bev,
    bevelSegments: shop.detail === 2 ? 2 : 1,
  });
  g.translate(0, 0, bev);
  // Extruded along z; stand the height up along y.
  g.rotateX(-Math.PI / 2);
  return g;
}
export const box = (w: number, h: number, d: number) => new THREE.BoxGeometry(w, h, d);

/** A cylinder along `axis`, radius `r` (and `r2` at its +axis end), length `h`. */
export function cyl(shop: Shop, r: number, h: number, axis: Axis = "y", n = 32, r2 = r, open = false): THREE.BufferGeometry {
  return alongAxis(new THREE.CylinderGeometry(r2, r, h, shop.seg(n, 8), 1, open), axis);
}

/** A tube: a ring of radius R with a round section of radius r, about `axis`, sweeping `arc` radians. */
export function torus(shop: Shop, R: number, r: number, axis: Axis = "y", arc = Math.PI * 2, n = 48, m = 10): THREE.BufferGeometry {
  const g = new THREE.TorusGeometry(R, r, shop.seg(m, 4), shop.seg(n, 8), arc);
  // TorusGeometry lies in xy about z; stand it about y by default.
  if (axis === "y") g.rotateX(-Math.PI / 2);
  else if (axis === "x") g.rotateY(Math.PI / 2);
  return g;
}

/** A solid of revolution about +y through the profile points [radius, y]. */
export function lathe(shop: Shop, profile: readonly (readonly [number, number])[], n = 48, phi0 = 0, phi = Math.PI * 2): THREE.BufferGeometry {
  const pts = profile.map(([r, y]) => new THREE.Vector2(Math.max(0, r), y));
  return new THREE.LatheGeometry(pts, shop.seg(n, 8), phi0, phi);
}

/** A sphere, or a cap of one when `theta` is less than π. */
export function ball(shop: Shop, r: number, n = 24, theta = Math.PI): THREE.BufferGeometry {
  return new THREE.SphereGeometry(r, shop.seg(n, 8), shop.seg(Math.max(6, n * 0.6), 4), 0, Math.PI * 2, 0, theta);
}

/**
 * A round bar bent along a polyline: straight runs joined by bends of radius `bend` (a handle, a rail, a pipe).
 * The bar is closed at both ends.
 */
export function bar(shop: Shop, pts: readonly V3[], r: number, bend = r * 3, radial = 12): THREE.BufferGeometry {
  const v = pts.map((p) => new THREE.Vector3(p[0], p[1], p[2]));
  const path = new THREE.CurvePath<THREE.Vector3>();
  let from = v[0]!.clone();
  for (let i = 1; i < v.length; i++) {
    const p = v[i]!;
    if (i === v.length - 1) {
      path.add(new THREE.LineCurve3(from, p.clone()));
      break;
    }
    const next = v[i + 1]!;
    const inDir = p.clone().sub(v[i - 1]!).normalize();
    const outDir = next.clone().sub(p).normalize();
    const lenIn = p.distanceTo(v[i - 1]!);
    const lenOut = next.distanceTo(p);
    const cut = Math.min(bend, lenIn * 0.49, lenOut * 0.49);
    const a = p.clone().addScaledVector(inDir, -cut);
    const b = p.clone().addScaledVector(outDir, cut);
    path.add(new THREE.LineCurve3(from, a));
    path.add(new THREE.QuadraticBezierCurve3(a, p.clone(), b));
    from = b;
  }
  const len = path.getLength();
  const steps = Math.max(2, Math.round(len / Math.max(0.004, bend * 0.35)));
  const tube = new THREE.TubeGeometry(path, Math.min(steps, shop.seg(240, 8)), r, shop.seg(radial, 5), false);
  const caps: THREE.BufferGeometry[] = [tube];
  if (shop.detail > 0) {
    for (const [p, q] of [
      [v[0]!, v[1]!],
      [v[v.length - 1]!, v[v.length - 2]!],
    ] as const) {
      const cap = new THREE.CircleGeometry(r, shop.seg(radial, 5));
      const dir = p.clone().sub(q).normalize();
      cap.applyQuaternion(new THREE.Quaternion().setFromUnitVectors(new THREE.Vector3(0, 0, 1), dir));
      cap.translate(p.x, p.y, p.z);
      caps.push(cap);
    }
  }
  return mergeGeometries(caps.map(normalise), false) ?? tube;
}

/**
 * A round wire or tube bent along a polyline, the same shape as `bar` (straight runs, each corner rounded by a curve
 * that starts `bend` before it), but sampled where it bends rather than evenly along its length: a straight run is a
 * single segment, a bend is split by its angle. Its frames are carried round the bends without twisting, and its
 * normals are the tube's own. Oven racks, wire shelving, ladders, long handles and rails cost a tenth of `bar`'s
 * triangles this way. The wire is closed at both ends when there is detail to show it.
 */
export function wire(shop: Shop, pts: readonly V3[], r: number, bend = r * 3, radial = 12): THREE.BufferGeometry {
  const v = pts.map((p) => new THREE.Vector3(p[0], p[1], p[2]));
  const step = shop.detail === 2 ? Math.PI / 14 : shop.detail === 1 ? Math.PI / 7 : Math.PI / 3.5;
  const line: THREE.Vector3[] = [v[0]!.clone()];
  const push = (p: THREE.Vector3) => {
    if (p.distanceToSquared(line[line.length - 1]!) > 1e-12) line.push(p);
  };
  for (let i = 1; i < v.length - 1; i++) {
    const p = v[i]!;
    const inDir = p.clone().sub(v[i - 1]!);
    const outDir = v[i + 1]!.clone().sub(p);
    const [lenIn, lenOut] = [inDir.length(), outDir.length()];
    inDir.normalize();
    outDir.normalize();
    const ang = Math.acos(Math.max(-1, Math.min(1, inDir.dot(outDir))));
    const cut = Math.min(bend, lenIn * 0.49, lenOut * 0.49);
    if (ang < 1e-3 || cut < 1e-6) {
      push(p.clone());
      continue;
    }
    const a = p.clone().addScaledVector(inDir, -cut);
    const b = p.clone().addScaledVector(outDir, cut);
    const k = Math.max(2, Math.ceil(ang / step));
    for (let j = 0; j <= k; j++) {
      const t = j / k;
      push(a.clone().multiplyScalar((1 - t) * (1 - t)).addScaledVector(p, 2 * (1 - t) * t).addScaledVector(b, t * t));
    }
  }
  push(v[v.length - 1]!.clone());
  const n = line.length;
  const m = shop.seg(radial, 5);
  // Tangents, and normals carried along by the rotation from each tangent to the next.
  const T = line.map((_, i) => line[Math.min(n - 1, i + 1)]!.clone().sub(line[Math.max(0, i - 1)]!).normalize());
  const N: THREE.Vector3[] = [];
  const t0 = T[0]!;
  const seed = Math.abs(t0.y) < 0.9 ? new THREE.Vector3(0, 1, 0) : new THREE.Vector3(1, 0, 0);
  N.push(seed.sub(t0.clone().multiplyScalar(seed.dot(t0))).normalize());
  const q = new THREE.Quaternion();
  for (let i = 1; i < n; i++) {
    q.setFromUnitVectors(T[i - 1]!, T[i]!);
    const nn = N[i - 1]!.clone().applyQuaternion(q);
    N.push(nn.sub(T[i]!.clone().multiplyScalar(nn.dot(T[i]!))).normalize());
  }
  const pos: number[] = [];
  const nrm: number[] = [];
  const uv: number[] = [];
  const idx: number[] = [];
  let run = 0;
  const d = new THREE.Vector3();
  for (let i = 0; i < n; i++) {
    if (i > 0) run += line[i]!.distanceTo(line[i - 1]!);
    const B = new THREE.Vector3().crossVectors(T[i]!, N[i]!);
    for (let j = 0; j <= m; j++) {
      const th = (j / m) * Math.PI * 2;
      d.copy(N[i]!).multiplyScalar(Math.cos(th)).addScaledVector(B, Math.sin(th));
      pos.push(line[i]!.x + d.x * r, line[i]!.y + d.y * r, line[i]!.z + d.z * r);
      nrm.push(d.x, d.y, d.z);
      uv.push(run, j / m);
    }
  }
  const w = m + 1;
  for (let i = 0; i < n - 1; i++)
    for (let j = 0; j < m; j++) {
      const [a, b, c, e] = [i * w + j, (i + 1) * w + j, i * w + j + 1, (i + 1) * w + j + 1];
      idx.push(a, c, b, c, e, b);
    }
  if (shop.detail > 0)
    for (const [end, sign] of [
      [0, -1],
      [n - 1, 1],
    ] as const) {
      const base = pos.length / 3;
      const t = T[end]!;
      pos.push(line[end]!.x, line[end]!.y, line[end]!.z);
      nrm.push(t.x * sign, t.y * sign, t.z * sign);
      uv.push(0, 0);
      for (let j = 0; j < m; j++) {
        const o = (end * w + j) * 3;
        pos.push(pos[o]!, pos[o + 1]!, pos[o + 2]!);
        nrm.push(t.x * sign, t.y * sign, t.z * sign);
        uv.push(0, 0);
      }
      for (let j = 0; j < m; j++) {
        const [a, b] = [base + 1 + j, base + 1 + ((j + 1) % m)];
        sign > 0 ? idx.push(base, a, b) : idx.push(base, b, a);
      }
    }
  const g = new THREE.BufferGeometry();
  g.setAttribute("position", new THREE.Float32BufferAttribute(pos, 3));
  g.setAttribute("normal", new THREE.Float32BufferAttribute(nrm, 3));
  g.setAttribute("uv", new THREE.Float32BufferAttribute(uv, 2));
  g.setIndex(idx);
  return g;
}

/**
 * A 2D section swept along a path of points: the section's x runs across the path (to its right, looking along it with
 * `up` overhead) and its y along `up`, turned at every point to the path's direction there. The section is a closed
 * outline, counter-clockwise; it is not repeated at its end. `scale(i)` stretches it across and up at point i (a
 * tapering finger, a flattened coil); `closed` joins the path's end to its start (a frame, a ring); an open path is
 * capped at both ends unless `caps` is false. Normals are smooth round the section and along the path, so a section
 * with sharp corners should repeat each corner point. UVs run along the path (u, metres) and round the section (v).
 */
export function sweep(
  path: readonly V3[],
  section: readonly (readonly [number, number])[],
  opts: { closed?: boolean; caps?: boolean; up?: V3; scale?: (i: number) => readonly [number, number] } = {},
): THREE.BufferGeometry {
  const closed = opts.closed ?? false;
  const P = path.map((p) => new THREE.Vector3(p[0], p[1], p[2]));
  const n = P.length;
  const m = section.length;
  const up0 = new THREE.Vector3(...(opts.up ?? [0, 1, 0]));
  const pos: number[] = [];
  const uv: number[] = [];
  let run = 0;
  const t = new THREE.Vector3();
  const sd = new THREE.Vector3();
  const u = new THREE.Vector3();
  for (let i = 0; i < n; i++) {
    const a = P[closed ? (i - 1 + n) % n : Math.max(0, i - 1)]!;
    const b = P[closed ? (i + 1) % n : Math.min(n - 1, i + 1)]!;
    t.subVectors(b, a).normalize();
    sd.crossVectors(t, up0);
    if (sd.lengthSq() < 1e-10) sd.set(1, 0, 0);
    sd.normalize();
    u.crossVectors(sd, t).normalize();
    // Keep the section's width through a corner: stretch it across by 1/cos of half the turn.
    let k = 1;
    if (closed || (i > 0 && i < n - 1)) {
      const t0 = new THREE.Vector3().subVectors(P[i]!, a).normalize();
      const t1 = new THREE.Vector3().subVectors(b, P[i]!).normalize();
      const c = t0.dot(t1);
      if (t0.lengthSq() > 0 && t1.lengthSq() > 0) k = 1 / Math.max(0.5, Math.sqrt((1 + c) / 2));
    }
    const [sx, sy] = opts.scale ? opts.scale(i) : [1, 1];
    if (i > 0) run += P[i]!.distanceTo(P[i - 1]!);
    for (let j = 0; j < m; j++) {
      const [x, y] = section[j]!;
      const q = P[i]!.clone().addScaledVector(sd, x * sx * k).addScaledVector(u, y * sy);
      pos.push(q.x, q.y, q.z);
      uv.push(run, j / m);
    }
  }
  const idx: number[] = [];
  const segs = closed ? n : n - 1;
  for (let i = 0; i < segs; i++) {
    const i1 = (i + 1) % n;
    for (let j = 0; j < m; j++) {
      const j1 = (j + 1) % m;
      const [A, B, C, D] = [i * m + j, i * m + j1, i1 * m + j1, i1 * m + j];
      idx.push(A, C, B, A, D, C);
    }
  }
  if (!closed && (opts.caps ?? true)) {
    const tri = THREE.ShapeUtils.triangulateShape(section.map(([x, y]) => new THREE.Vector2(x, y)), []);
    // A counter-clockwise section's own triangles face back down the path: right for the start cap, flipped for the end.
    for (const [end, flip] of [
      [0, false],
      [n - 1, true],
    ] as const) {
      const base = pos.length / 3;
      for (let j = 0; j < m; j++) {
        const o = (end * m + j) * 3;
        pos.push(pos[o]!, pos[o + 1]!, pos[o + 2]!);
        uv.push(0, 0);
      }
      for (const [a, b, c] of tri) flip ? idx.push(base + a!, base + c!, base + b!) : idx.push(base + a!, base + b!, base + c!);
    }
  }
  const g = new THREE.BufferGeometry();
  g.setAttribute("position", new THREE.Float32BufferAttribute(pos, 3));
  g.setAttribute("uv", new THREE.Float32BufferAttribute(uv, 2));
  g.setIndex(idx);
  g.computeVertexNormals();
  return g;
}

/**
 * A cast bar's section, top at y = 0 and `h` deep, `w` across its crown: drafted sides flaring to `w · flare` at the
 * foot, a shallow elliptical crown `crown` high, the foot's corners eased. Counter-clockwise, for `sweep`. A crowned top
 * meets a crossing bar's crown along a line, never over an area, so cast grids never show coplanar flicker.
 */
export function castSection(w: number, h: number, crown = Math.min(0.004, w * 0.32), flare = 1.22, arc = 6): [number, number][] {
  const out: [number, number][] = [];
  const wb = (w * flare) / 2;
  const ease = Math.min(0.0015, h * 0.12);
  out.push([-wb + ease, -h], [wb - ease, -h], [wb, -h + ease]);
  for (let i = 0; i <= arc; i++) {
    const a = (i / arc) * Math.PI;
    out.push([(w / 2) * Math.cos(a), -crown + crown * Math.sin(a)]);
  }
  out.push([-wb, -h + ease]);
  return out;
}

/** A flat plate cut from `shape` (in its own xy), `t` thick, standing along +z from 0, with an edge chamfer `bevel`. */
export function plate(shop: Shop, shape: THREE.Shape, t: number, bevel = 0, curve = 12): THREE.BufferGeometry {
  const g = new THREE.ExtrudeGeometry(shape, {
    depth: Math.max(1e-5, t - 2 * bevel),
    bevelEnabled: bevel > 0 && shop.detail > 0,
    bevelThickness: bevel,
    bevelSize: bevel,
    bevelSegments: shop.detail === 2 ? 2 : 1,
    curveSegments: shop.seg(curve, 3),
  });
  if (bevel > 0 && shop.detail > 0) g.translate(0, 0, bevel);
  return g;
}

/** A rectangle shape centred on the origin, its corners rounded by `r`. */
export function roundRect(w: number, h: number, r: number, cx = 0, cy = 0): THREE.Shape {
  const s = new THREE.Shape();
  rectPath(s, w, h, r, cx, cy);
  return s;
}
/** Trace a rounded rectangle into a path (a shape, or a hole when traced the other way with `reverse`). */
export function rectPath(p: THREE.Path, w: number, h: number, r: number, cx = 0, cy = 0, reverse = false): THREE.Path {
  const rr = Math.max(0, Math.min(r, w / 2, h / 2));
  const [x0, x1, y0, y1] = [cx - w / 2, cx + w / 2, cy - h / 2, cy + h / 2];
  if (!reverse) {
    p.moveTo(x0 + rr, y0);
    p.lineTo(x1 - rr, y0);
    if (rr) p.quadraticCurveTo(x1, y0, x1, y0 + rr);
    p.lineTo(x1, y1 - rr);
    if (rr) p.quadraticCurveTo(x1, y1, x1 - rr, y1);
    p.lineTo(x0 + rr, y1);
    if (rr) p.quadraticCurveTo(x0, y1, x0, y1 - rr);
    p.lineTo(x0, y0 + rr);
    if (rr) p.quadraticCurveTo(x0, y0, x0 + rr, y0);
  } else {
    p.moveTo(x0 + rr, y0);
    if (rr) p.quadraticCurveTo(x0, y0, x0, y0 + rr);
    p.lineTo(x0, y1 - rr);
    if (rr) p.quadraticCurveTo(x0, y1, x0 + rr, y1);
    p.lineTo(x1 - rr, y1);
    if (rr) p.quadraticCurveTo(x1, y1, x1, y1 - rr);
    p.lineTo(x1, y0 + rr);
    if (rr) p.quadraticCurveTo(x1, y0, x1 - rr, y0);
    p.lineTo(x0 + rr, y0);
  }
  return p;
}
/** A round hole for a shape. */
export function circleHole(cx: number, cy: number, r: number): THREE.Path {
  return new THREE.Path().absarc(cx, cy, r, 0, Math.PI * 2, true);
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Hardware every fixture uses.                                                                                  */
/* ------------------------------------------------------------------------------------------------------------ */

/** A socket-head or pan-head screw's head, facing +z, base at z = 0. Only drawn close up. */
export function screw(shop: Shop, mat: THREE.Material, at: V3, r = 0.003, rot?: V3): void {
  if (!shop.fine(2)) return;
  const g = lathe(shop, [
    [0, 0],
    [r, 0],
    [r, r * 0.35],
    [r * 0.8, r * 0.62],
    [r * 0.42, r * 0.66],
    [r * 0.42, r * 0.4],
    [0, r * 0.4],
  ], 14);
  g.rotateX(Math.PI / 2);
  shop.add(mat, g, at, rot);
}

/**
 * A commercial range or appliance knob, facing +z from a base at z = 0: a skirted bakelite body with grip flutes and
 * a pointer ridge, on a chrome bezel ring.
 */
export function knob(shop: Shop, body: THREE.Material, bezel: THREE.Material, at: V3, r = 0.022, turn = 0): void {
  shop.at(at, [0, 0, turn], () => {
    shop.add(bezel, cyl(shop, r * 1.18, 0.002, "z", 32), [0, 0, 0.001]);
    const prof: [number, number][] = [
      [0, 0.002],
      [r, 0.002],
      [r, 0.006],
      [r * 0.82, 0.009],
      [r * 0.72, 0.012],
      [r * 0.7, 0.024],
      [r * 0.6, 0.028],
      [0, 0.029],
    ];
    const g = lathe(shop, prof, 36);
    g.rotateX(Math.PI / 2);
    shop.add(body, g);
    // Pointer ridge across the grip, and close up the flutes round it that a thumb finds.
    shop.add(body, rbox(shop, r * 0.28, r * 1.3, 0.012, 0.002), [0, r * 0.2, 0.022]);
    if (shop.fine(2))
      for (let i = 0; i < 12; i++) {
        const a = (i / 12) * Math.PI * 2;
        shop.add(body, rbox(shop, r * 0.16, r * 0.16, 0.011, 0.0012), [Math.cos(a) * r * 0.7, Math.sin(a) * r * 0.7, 0.018]);
      }
  });
}

/** A bar pull handle on two standoffs, along x, its bar `len` long and standing `off` proud of the face at z = 0. */
export function pullHandle(shop: Shop, mat: THREE.Material, at: V3, len: number, off = 0.04, r = 0.0095, axis: "x" | "y" = "x"): void {
  shop.at(at, axis === "x" ? [0, 0, 0] : [0, 0, Math.PI / 2], () => {
    const h = len / 2;
    shop.add(mat, bar(shop, [[-h, 0, 0], [-h, 0, off], [h, 0, off], [h, 0, 0]], r, off * 0.45, 14));
    for (const x of [-h, h]) shop.add(mat, cyl(shop, r * 1.5, 0.004, "z", 16), [x, 0, 0.002]);
  });
}

/** A swivel caster, its wheel touching y = 0 under the plate centred at `at` (the plate's underside). */
export function caster(shop: Shop, wheel: THREE.Material, metal: THREE.Material, at: V3, r = 0.05, brake = false): void {
  const [x, y, z] = at;
  shop.add(metal, rbox(shop, 0.07, 0.006, 0.07, 0.002), [x, y - 0.003, z]);
  shop.add(metal, cyl(shop, 0.012, 0.014, "y", 16), [x, y - 0.013, z]);
  const cy = r;
  // Fork legs either side of the wheel.
  for (const s of [-1, 1]) shop.add(metal, rbox(shop, 0.004, y - 0.02 - cy + 0.012, r * 1.1, 0.0015), [x + s * 0.022, (y - 0.02 + cy - 0.012) / 2 + 0.006, z - r * 0.25]);
  shop.add(metal, rbox(shop, 0.048, 0.006, r * 1.1, 0.002), [x, y - 0.022, z - r * 0.25]);
  const tread = lathe(shop, [
    [r * 0.55, -0.016],
    [r * 0.92, -0.016],
    [r, -0.01],
    [r, 0.01],
    [r * 0.92, 0.016],
    [r * 0.55, 0.016],
  ], 32);
  tread.rotateZ(Math.PI / 2);
  shop.add(wheel, tread, [x, cy, z - r * 0.4]);
  shop.add(metal, cyl(shop, r * 0.56, 0.03, "x", 20), [x, cy, z - r * 0.4]);
  shop.add(metal, cyl(shop, 0.006, 0.05, "x", 10), [x, cy, z - r * 0.4]);
  if (brake && shop.fine(1)) shop.add(wheel, rbox(shop, 0.03, 0.006, 0.03, 0.002), [x, cy + r + 0.004, z + r * 0.45], [-0.35, 0, 0]);
}

/** An adjustable bullet foot under a stand leg: a tapered sleeve and a domed pad. */
export function bulletFoot(shop: Shop, mat: THREE.Material, at: V3, r = 0.02, h = 0.1): void {
  shop.at(at, [0, 0, 0], () => {
    shop.add(mat, lathe(shop, [
      [0, 0],
      [r * 0.9, 0],
      [r, 0.006],
      [r * 0.95, h * 0.35],
      [r * 0.72, h * 0.8],
      [r * 0.72, h],
      [0, h],
    ], 24));
  });
}

/** Square tube legs of a stand or table: a leg of `s` section, rounded, from y0 to y1 at (x, z). */
export function squareLeg(shop: Shop, mat: THREE.Material, x: number, z: number, y0: number, y1: number, s = 0.04): void {
  shop.add(mat, rbox(shop, s, y1 - y0, s, s * 0.18), [x, (y0 + y1) / 2, z]);
}

/** Round tube legs, the stainless 41 mm kind under tables, with a bullet foot. */
export function roundLeg(shop: Shop, mat: THREE.Material, x: number, z: number, y0: number, y1: number, r = 0.0205): void {
  shop.add(mat, cyl(shop, r, y1 - y0 - 0.07, "y", 20), [x, (y0 + 0.07 + y1) / 2, z]);
  bulletFoot(shop, mat, [x, y0, z], r * 1.05, 0.075);
}

/** Place a built fixture as a manifest pose says: back-left corner at (x, y) on the plan, base at z, turned by yaw. */
export function placeFixture(o: THREE.Object3D, pose: { x_m: number; y_m: number; z_m: number; yaw_rad: number }): void {
  o.position.set(pose.x_m, pose.z_m, -pose.y_m);
  o.rotation.set(0, pose.yaw_rad, 0);
}
