/**
 * Kitchen staff: procedural, animated people in chef whites, built to the anthropometric cohorts of
 * assets/kitchens/specs/anthropometrics.tsv, and a motion layer that walks them along floor paths on planted feet and
 * works them at stations with two-bone arm and leg IK inside the table's ranges of motion.
 *
 * Every person is synthetic. Heads are sculpted like an artist's mannequin (skull, jaw, ears, neck, hair under a
 * cap) with no eyes, nose or mouth, so nothing biometric is drawn.
 *
 * Body frame (metres): the root sits on the floor between the feet, y up, the person facing +z, their left at +x.
 * The table's heights are barefoot; shoes lift the whole body by the row's shoe sole (`shoe_sole_m`), so a shod
 * person stands that much taller, as they do at work.
 *
 * Meshes: one SkinnedMesh per material per person (whites, apron, trousers, skin, shoes, hair), six draw calls,
 * every one bound to the same skeleton. Cloth is built as lofts shaped round the body; vertices are rigid (weight 1
 * to a bone) except where cloth or skin crosses a joint, where two bones blend so that joints stay closed.
 *
 * Skeleton: the 29 bones of the old rigs (staff_rig_manifest.json), names unchanged, plus 15 finger bones a hand:
 * thumb_01..03, index_01..03, middle_01..03, ring_01..03, pinky_01..03, suffixed _L and _R (59 bones in all). Every
 * deforming bone's bind rotation is the identity, so a bone's local rotation is a rotation in the body frame.
 *
 * Grip frame (`hand(side)`, the grip_L/grip_R bones): the handle of a held tool runs along +z through the fist, from
 * the heel of the palm out past the thumb and forefinger towards the tool's working end; +y is the thumb side (the
 * top of a pan, the spine of a knife). Model a tool with its handle's middle at the origin, the handle along +z and
 * its top along +y, and add it as a child of `hand("R")`.
 */
import * as THREE from "three";
import { GLTFLoader, type GLTF } from "three/examples/jsm/loaders/GLTFLoader.js";
import { MeshoptDecoder } from "three/examples/jsm/libs/meshopt_decoder.module.js";
import type { Detail } from "./shop";

export type Cohort = "p05" | "p50" | "p95";
export type Side = "L" | "R";
export type WorkKind = "stand" | "saute" | "stir" | "plate" | "chop" | "carry" | "push" | "lift";
export type RomKey =
  | "neck_pitch"
  | "neck_yaw"
  | "shoulder_flex"
  | "shoulder_abd"
  | "elbow_flex"
  | "wrist_flex"
  | "torso_flex"
  | "torso_yaw"
  | "hip_flex"
  | "knee_flex"
  | "ankle_pitch";

/** One cohort's row of anthropometrics.tsv (metres, degrees), with the body mass the brigade chapter uses. */
export interface Anthro {
  readonly percentile_id: string;
  readonly cohort_label: string;
  readonly standing_height_m: number;
  readonly eye_height_m: number;
  readonly shoulder_height_m: number;
  readonly elbow_height_m: number;
  readonly knuckle_height_m: number;
  readonly hip_height_m: number;
  readonly knee_height_m: number;
  readonly ankle_height_m: number;
  readonly shoulder_breadth_m: number;
  readonly chest_depth_m: number;
  readonly upper_arm_m: number;
  readonly forearm_m: number;
  readonly hand_length_m: number;
  readonly forward_reach_m: number;
  readonly overhead_reach_m: number;
  readonly shoe_length_m: number;
  readonly shoe_width_m: number;
  readonly shoe_sole_m: number;
  readonly torso_capsule_radius_m: number;
  readonly joint_rom_deg: Readonly<Record<RomKey, readonly [number, number]>>;
  readonly source_citation: string;
  /** Not in the table: the masses the brigade chapter feeds the L5/S1 torque (52, 74, 98 kg). */
  readonly body_mass_kg: number;
}

const ROM: Record<RomKey, readonly [number, number]> = {
  neck_pitch: [-45, 45],
  neck_yaw: [-65, 65],
  shoulder_flex: [-45, 165],
  shoulder_abd: [0, 150],
  elbow_flex: [0, 142],
  wrist_flex: [-65, 70],
  torso_flex: [-20, 55],
  torso_yaw: [-40, 40],
  hip_flex: [-15, 115],
  knee_flex: [0, 135],
  ankle_pitch: [-35, 20],
};

function row(id: string, label: string, n: readonly number[], mass: number, cite: string): Anthro {
  const [h, eye, sh, el, kn, hip, knee, ank, sb, cd, ua, fa, hl, fr, or, sl, sw, ss, tr] = n as [
    number, number, number, number, number, number, number, number, number, number,
    number, number, number, number, number, number, number, number, number,
  ];
  return {
    percentile_id: id,
    cohort_label: label,
    standing_height_m: h,
    eye_height_m: eye,
    shoulder_height_m: sh,
    elbow_height_m: el,
    knuckle_height_m: kn,
    hip_height_m: hip,
    knee_height_m: knee,
    ankle_height_m: ank,
    shoulder_breadth_m: sb,
    chest_depth_m: cd,
    upper_arm_m: ua,
    forearm_m: fa,
    hand_length_m: hl,
    forward_reach_m: fr,
    overhead_reach_m: or,
    shoe_length_m: sl,
    shoe_width_m: sw,
    shoe_sole_m: ss,
    torso_capsule_radius_m: tr,
    joint_rom_deg: ROM,
    source_citation: cite,
    body_mass_kg: mass,
  };
}

/** The three cohorts, transcribed from assets/kitchens/specs/anthropometrics.tsv. */
export const ANTHROPOMETRICS: Readonly<Record<Cohort, Anthro>> = {
  p05: row("P05_FEMALE", "5th Percentile Adult Female Chef",
    [1.511, 1.402, 1.234, 0.942, 0.678, 0.785, 0.422, 0.068, 0.382, 0.208, 0.305, 0.232, 0.165, 0.676, 1.854, 0.242, 0.092, 0.026, 0.175], 52,
    "ISO 7250-1:2017; NASA-STD-3000 Vol. I sec. 3.3.1 (5th percentile female anthropometry & joint ROM)"),
  p50: row("P50_ADULT", "50th Percentile Combined Adult Chef",
    [1.715, 1.6, 1.406, 1.074, 0.758, 0.892, 0.49, 0.076, 0.448, 0.238, 0.345, 0.264, 0.186, 0.772, 2.115, 0.276, 0.104, 0.028, 0.205], 74,
    "ISO 7250-1:2017; NASA-STD-3000 Vol. I sec. 3.3.1 (50th percentile adult midpoint anthropometry & joint ROM)"),
  p95: row("P95_MALE", "95th Percentile Adult Male Chef",
    [1.878, 1.76, 1.552, 1.186, 0.832, 0.985, 0.548, 0.084, 0.506, 0.272, 0.38, 0.292, 0.206, 0.862, 2.34, 0.308, 0.116, 0.03, 0.235], 98,
    "ISO 7250-1:2017; NASA-STD-3000 Vol. I sec. 3.3.1 (95th percentile male anthropometry & joint ROM; ADA 36 in clearance fit check)"),
};

/* ------------------------------------------------------------------------------------------------------------ */
/* Small maths                                                                                                    */
/* ------------------------------------------------------------------------------------------------------------ */

const DEG = Math.PI / 180;
const clamp = (x: number, a: number, b: number) => (x < a ? a : x > b ? b : x);
const lerp = (a: number, b: number, t: number) => a + (b - a) * t;
/** Hermite smoothstep of x between e0 and e1 (either order). */
function smooth(e0: number, e1: number, x: number): number {
  const t = clamp((x - e0) / (e1 - e0), 0, 1);
  return t * t * (3 - 2 * t);
}
/** Minimum-jerk position profile on [0, 1]. */
function mj(u: number): number {
  const t = clamp(u, 0, 1);
  return t * t * t * (10 + t * (-15 + 6 * t));
}
/** Ease in and out on [0, 1] (cosine). */
function ease(u: number): number {
  return 0.5 - 0.5 * Math.cos(Math.PI * clamp(u, 0, 1));
}
function wrapPi(a: number): number {
  let x = (a + Math.PI) % (Math.PI * 2);
  if (x < 0) x += Math.PI * 2;
  return x - Math.PI;
}
/** Sign-preserving power, for superellipses. */
function spow(x: number, p: number): number {
  return Math.sign(x) * Math.pow(Math.abs(x), p);
}
/** Soft minimum of a and b with blending width k. */
function smin(a: number, b: number, k: number): number {
  const h = clamp(0.5 + (0.5 * (b - a)) / k, 0, 1);
  return lerp(b, a, h) - k * h * (1 - h);
}

function rng(seed: number): () => number {
  let s = (seed * 2654435761) >>> 0 || 1;
  return () => {
    s ^= s << 13;
    s >>>= 0;
    s ^= s >>> 17;
    s ^= s << 5;
    s >>>= 0;
    return s / 4294967296;
  };
}

const V = (x = 0, y = 0, z = 0) => new THREE.Vector3(x, y, z);
const X = V(1, 0, 0);
const Y = V(0, 1, 0);
const Z = V(0, 0, 1);
const ZERO = V(0, 0, 0);

/** Monotone cubic interpolation (Fritsch-Carlson) through (x_i, y_i), x increasing; clamps outside. */
class Mono {
  private m: number[];
  constructor(readonly xs: readonly number[], readonly ys: readonly number[]) {
    const n = xs.length;
    const d: number[] = [];
    for (let i = 0; i < n - 1; i++) d.push((ys[i + 1]! - ys[i]!) / Math.max(1e-9, xs[i + 1]! - xs[i]!));
    const m: number[] = new Array<number>(n).fill(0);
    if (n > 1) {
      m[0] = d[0]!;
      m[n - 1] = d[n - 2]!;
      for (let i = 1; i < n - 1; i++) m[i] = d[i - 1]! * d[i]! <= 0 ? 0 : (d[i - 1]! + d[i]!) / 2;
      for (let i = 0; i < n - 1; i++) {
        if (d[i] === 0) {
          m[i] = 0;
          m[i + 1] = 0;
          continue;
        }
        const a = m[i]! / d[i]!;
        const b = m[i + 1]! / d[i]!;
        const h = a * a + b * b;
        if (h > 9) {
          const t = 3 / Math.sqrt(h);
          m[i] = t * a * d[i]!;
          m[i + 1] = t * b * d[i]!;
        }
      }
    }
    this.m = m;
  }
  /** Force the end slopes (for starting and stopping at rest). */
  ends(m0: number, m1: number): this {
    this.m[0] = m0;
    this.m[this.m.length - 1] = m1;
    return this;
  }
  at(x: number): number {
    const xs = this.xs;
    const n = xs.length;
    if (n === 1 || x <= xs[0]!) return this.ys[0]!;
    if (x >= xs[n - 1]!) return this.ys[n - 1]!;
    let lo = 0;
    let hi = n - 1;
    while (hi - lo > 1) {
      const mid = (lo + hi) >> 1;
      if (xs[mid]! <= x) lo = mid;
      else hi = mid;
    }
    const h = xs[hi]! - xs[lo]!;
    const t = (x - xs[lo]!) / h;
    const t2 = t * t;
    const t3 = t2 * t;
    return (2 * t3 - 3 * t2 + 1) * this.ys[lo]! + (t3 - 2 * t2 + t) * h * this.m[lo]! + (-2 * t3 + 3 * t2) * this.ys[hi]! + (t3 - t2) * h * this.m[hi]!;
  }
}

const _fm0 = new THREE.Matrix4();
const _fm1 = new THREE.Matrix4();
const _fc0 = V();
const _fc1 = V();
/** The rotation taking the orthonormal pair (a0, b0) onto (a, b). */
function frameQuat(a0: THREE.Vector3, b0: THREE.Vector3, a: THREE.Vector3, b: THREE.Vector3, out: THREE.Quaternion): THREE.Quaternion {
  _fm0.makeBasis(a0, b0, _fc0.crossVectors(a0, b0));
  _fm1.makeBasis(a, b, _fc1.crossVectors(a, b));
  _fm0.transpose();
  _fm1.multiply(_fm0);
  return out.setFromRotationMatrix(_fm1);
}
/** Make b perpendicular to unit a and unit length (falling back to `alt`). */
function orthoTo(b: THREE.Vector3, a: THREE.Vector3, alt: THREE.Vector3 = X): THREE.Vector3 {
  b.addScaledVector(a, -b.dot(a));
  if (b.lengthSq() < 1e-10) b.copy(alt).addScaledVector(a, -alt.dot(a));
  return b.normalize();
}
/** Swing-twist: the twist angle of q about unit axis a. */
function twistAngle(q: THREE.Quaternion, a: THREE.Vector3): number {
  const p = q.x * a.x + q.y * a.y + q.z * a.z;
  return 2 * Math.atan2(p, q.w);
}
const qAxis = (axis: THREE.Vector3, ang: number, out = new THREE.Quaternion()) => out.setFromAxisAngle(axis, ang);


/* ------------------------------------------------------------------------------------------------------------ */
/* Cloth, skin and hair. Drawn once on canvases from fixed seeds; shared by every person (colours ride on vertex   */
/* colours, so one material serves every apron colour and skin tone).                                           */
/* ------------------------------------------------------------------------------------------------------------ */

function canvas2d(w: number, h: number): [HTMLCanvasElement, CanvasRenderingContext2D] {
  const c = document.createElement("canvas");
  c.width = w;
  c.height = h;
  return [c, c.getContext("2d")!];
}
/** A texture whose UVs are metres, repeating every `tile` metres. */
function canvasTex(c: HTMLCanvasElement, srgb: boolean, tile: number): THREE.CanvasTexture {
  const t = new THREE.CanvasTexture(c);
  t.colorSpace = srgb ? THREE.SRGBColorSpace : THREE.NoColorSpace;
  t.wrapS = t.wrapT = THREE.RepeatWrapping;
  t.repeat.set(1 / tile, 1 / tile);
  t.anisotropy = 8;
  t.needsUpdate = true;
  return t;
}
/** A tangent-space normal map from a height function on the unit square (periodic), `k` its strength. */
function normalMap(S: number, h: (u: number, v: number) => number, k: number): HTMLCanvasElement {
  const [c, x] = canvas2d(S, S);
  const hs = new Float32Array(S * S);
  for (let j = 0; j < S; j++) for (let i = 0; i < S; i++) hs[j * S + i] = h(i / S, j / S);
  const img = x.createImageData(S, S);
  for (let j = 0; j < S; j++)
    for (let i = 0; i < S; i++) {
      const dx = (hs[j * S + ((i + 1) % S)]! - hs[j * S + ((i + S - 1) % S)]!) * k;
      const dy = (hs[((j + 1) % S) * S + i]! - hs[((j + S - 1) % S) * S + i]!) * k;
      const l = Math.hypot(dx, dy, 1);
      const o = (j * S + i) * 4;
      img.data[o] = Math.round((0.5 - (0.5 * dx) / l) * 255);
      img.data[o + 1] = Math.round((0.5 + (0.5 * dy) / l) * 255);
      img.data[o + 2] = Math.round((0.5 + 0.5 / l) * 255);
      img.data[o + 3] = 255;
    }
  x.putImageData(img, 0, 0);
  return c;
}

/** 2/1 twill: diagonal wales with a little yarn noise. One tile is `n` wales. */
function twill(seed: number, n: number): HTMLCanvasElement {
  const r = rng(seed);
  const noise = new Float32Array(64 * 64).map(() => r());
  const nz = (u: number, v: number) => noise[(Math.floor(v * 64) % 64) * 64 + (Math.floor(u * 64) % 64)]!;
  return normalMap(256, (u, v) => {
    const w = (u + v) * n;
    const f = w - Math.floor(w);
    const ridge = Math.sin(Math.PI * f) ** 0.6;
    const yarn = 0.25 * Math.sin(Math.PI * 2 * (u - v) * n * 3);
    return ridge + yarn * ridge + 0.35 * nz(u, v);
  }, 1.4);
}

/** Houndstooth as weavers make it: a 2/2 twill with four dark and four light threads each way. */
function houndstooth(): HTMLCanvasElement {
  const S = 256;
  const k = 4;
  const [c, x] = canvas2d(S, S);
  const r = rng(77);
  for (let j = 0; j < S / k; j++)
    for (let i = 0; i < S / k; i++) {
      const warpUp = (i + j) % 4 < 2;
      const dark = warpUp ? i % 8 < 4 : j % 8 < 4;
      const v = dark ? 26 + r() * 14 : 222 + r() * 22;
      x.fillStyle = `rgb(${v | 0},${v | 0},${(v + (dark ? 2 : -4)) | 0})`;
      x.fillRect(i * k, j * k, k, k);
      // A thread's shading: lighter along its middle.
      x.fillStyle = dark ? "rgba(255,255,255,0.05)" : "rgba(0,0,0,0.07)";
      if (warpUp) x.fillRect(i * k, j * k, 1, k);
      else x.fillRect(i * k, j * k, k, 1);
    }
  return c;
}

/** Fine skin: pores and a soft mottling, as a bump map. */
function skinBump(): HTMLCanvasElement {
  const S = 256;
  const [c, x] = canvas2d(S, S);
  const r = rng(4242);
  x.fillStyle = "#808080";
  x.fillRect(0, 0, S, S);
  for (let i = 0; i < 900; i++) {
    const v = 110 + r() * 40;
    x.fillStyle = `rgba(${v | 0},${v | 0},${v | 0},0.25)`;
    x.beginPath();
    x.arc(r() * S, r() * S, 4 + r() * 14, 0, Math.PI * 2);
    x.fill();
  }
  for (let i = 0; i < 5200; i++) {
    x.fillStyle = `rgba(40,40,40,${(0.25 + r() * 0.35).toFixed(2)})`;
    x.fillRect(r() * S, r() * S, 1, 1);
  }
  return c;
}

/** Hair: fine strands running down the texture's v, light and dark. */
function hairStrands(): HTMLCanvasElement {
  const S = 256;
  const [c, x] = canvas2d(S, S);
  const r = rng(9091);
  x.fillStyle = "#9a9a9a";
  x.fillRect(0, 0, S, S);
  for (let i = 0; i < 2600; i++) {
    const v = r() < 0.5 ? 40 + r() * 60 : 170 + r() * 85;
    x.strokeStyle = `rgba(${v | 0},${v | 0},${v | 0},${(0.25 + r() * 0.4).toFixed(2)})`;
    x.lineWidth = 0.6 + r() * 0.9;
    const x0 = r() * S;
    const y0 = r() * S;
    const len = 20 + r() * 70;
    const bend = (r() - 0.5) * 8;
    x.beginPath();
    for (const dx of [-S, 0]) {
      x.moveTo(x0 + dx, y0);
      x.quadraticCurveTo(x0 + dx + bend, y0 + len / 2, x0 + dx + bend * 0.3, y0 + len);
    }
    x.stroke();
  }
  return c;
}

export interface PeopleMaterials {
  /** Chef whites: jacket, collar, knotted buttons, caps. Cotton twill with a cloth sheen. */
  whites: THREE.MeshPhysicalMaterial;
  /** Bib aprons: heavier canvas, coloured by the person's vertex colour. */
  apron: THREE.MeshPhysicalMaterial;
  /** Check (houndstooth) trousers. */
  check: THREE.MeshPhysicalMaterial;
  /** Plain trousers (black or dark), coloured by vertex colour. */
  plain: THREE.MeshPhysicalMaterial;
  /** Skin: hands, forearms, neck and the faceless head; tone by vertex colour. */
  skin: THREE.MeshPhysicalMaterial;
  /** Shoes and clogs: glazed leather uppers and rubber soles, by vertex colour. */
  shoes: THREE.MeshPhysicalMaterial;
  /** Hair, with an anisotropic sheen along the strands. */
  hair: THREE.MeshPhysicalMaterial;
}

let matCache: PeopleMaterials | null = null;

/** The shared materials, made the first time they are asked for. */
export function peopleMaterials(): PeopleMaterials {
  if (matCache) return matCache;
  const fine = canvasTex(twill(31, 24), false, 0.03);
  const canvasWeave = canvasTex(twill(57, 14), false, 0.03);
  const hound = canvasTex(houndstooth(), true, 0.056);
  const houndN = canvasTex(twill(73, 28), false, 0.028);
  const skinB = canvasTex(skinBump(), false, 0.05);
  const strands = canvasTex(hairStrands(), true, 0.06);
  const cloth = (p: THREE.MeshPhysicalMaterialParameters) =>
    new THREE.MeshPhysicalMaterial({ vertexColors: true, side: THREE.DoubleSide, metalness: 0, ...p });
  matCache = {
    whites: cloth({ color: 0xffffff, roughness: 0.84, normalMap: fine, normalScale: new THREE.Vector2(0.45, 0.45), sheen: 0.55, sheenRoughness: 0.75, sheenColor: new THREE.Color(0xffffff) }),
    apron: cloth({ color: 0xffffff, roughness: 0.88, normalMap: canvasWeave, normalScale: new THREE.Vector2(0.6, 0.6), sheen: 0.35, sheenRoughness: 0.8, sheenColor: new THREE.Color(0xd8dce4) }),
    check: cloth({ color: 0xffffff, roughness: 0.86, map: hound, normalMap: houndN, normalScale: new THREE.Vector2(0.4, 0.4), sheen: 0.3, sheenRoughness: 0.8, sheenColor: new THREE.Color(0xcccccc) }),
    plain: cloth({ color: 0xffffff, roughness: 0.82, normalMap: houndN, normalScale: new THREE.Vector2(0.45, 0.45), sheen: 0.45, sheenRoughness: 0.7, sheenColor: new THREE.Color(0x8a8f99) }),
    skin: new THREE.MeshPhysicalMaterial({
      vertexColors: true,
      color: 0xffffff,
      roughness: 0.56,
      metalness: 0,
      bumpMap: skinB,
      bumpScale: 0.12,
      sheen: 0.35,
      sheenRoughness: 0.55,
      sheenColor: new THREE.Color(0x9c4a3a),
      clearcoat: 0.06,
      clearcoatRoughness: 0.55,
    }),
    shoes: new THREE.MeshPhysicalMaterial({ vertexColors: true, color: 0xffffff, roughness: 0.42, metalness: 0, clearcoat: 0.55, clearcoatRoughness: 0.32 }),
    hair: new THREE.MeshPhysicalMaterial({ vertexColors: true, color: 0xffffff, roughness: 0.5, metalness: 0, map: strands, anisotropy: 0.7, anisotropyRotation: Math.PI / 2, sheen: 0.2, sheenRoughness: 0.4 }),
  };
  return matCache;
}


/* ------------------------------------------------------------------------------------------------------------ */
/* Skinned geometry: a vertex soup per material, filled by lofts (rings swept along a body part) and primitives.  */
/* ------------------------------------------------------------------------------------------------------------ */

/** Bone influences: [bone index, weight] pairs (normalised on the way in, four at most kept). */
type Wts = readonly (readonly [number, number])[];
type WFn = (p: THREE.Vector3) => Wts;

class Soup {
  pos: number[] = [];
  nrm: number[] = [];
  uv: number[] = [];
  col: number[] = [];
  si: number[] = [];
  sw: number[] = [];
  idx: number[] = [];
  get verts(): number {
    return this.pos.length / 3;
  }
  get tris(): number {
    return this.idx.length / 3;
  }
  vert(p: THREE.Vector3, n: THREE.Vector3, u: number, v: number, c: THREE.Color, w: Wts): number {
    const i = this.verts;
    this.pos.push(p.x, p.y, p.z);
    this.nrm.push(n.x, n.y, n.z);
    this.uv.push(u, v);
    this.col.push(c.r, c.g, c.b);
    const acc: [number, number][] = [];
    for (const [b, x] of w) {
      if (x <= 1e-4) continue;
      const e = acc.find((a) => a[0] === b);
      if (e) e[1] += x;
      else acc.push([b, x]);
    }
    acc.sort((a, b) => b[1] - a[1]);
    acc.length = Math.min(acc.length, 4);
    let sum = 0;
    for (const a of acc) sum += a[1];
    for (let k = 0; k < 4; k++) {
      const a = acc[k];
      this.si.push(a ? a[0] : 0);
      this.sw.push(a && sum > 0 ? a[1] / sum : k === 0 && !a ? 1 : 0);
    }
    return i;
  }
  tri(a: number, b: number, c: number): void {
    this.idx.push(a, b, c);
  }
  geometry(): THREE.BufferGeometry {
    const g = new THREE.BufferGeometry();
    g.setAttribute("position", new THREE.Float32BufferAttribute(this.pos, 3));
    g.setAttribute("normal", new THREE.Float32BufferAttribute(this.nrm, 3));
    g.setAttribute("uv", new THREE.Float32BufferAttribute(this.uv, 2));
    g.setAttribute("color", new THREE.Float32BufferAttribute(this.col, 3));
    g.setAttribute("skinIndex", new THREE.Uint16BufferAttribute(this.si, 4));
    g.setAttribute("skinWeight", new THREE.Float32BufferAttribute(this.sw, 4));
    g.setIndex(this.verts > 65535 ? new THREE.Uint32BufferAttribute(this.idx, 1) : new THREE.Uint16BufferAttribute(this.idx, 1));
    return g;
  }
}

type ColorFn = THREE.Color | ((i: number, j: number, p: THREE.Vector3) => THREE.Color);
interface LoftOpts {
  /** Rings are closed loops (default) or open strips. */
  closed?: boolean;
  /** Influences of the vertex at ring i, point j. */
  w: (i: number, j: number, p: THREE.Vector3) => Wts;
  color: ColorFn;
  /** A direction the surface should face at a point (default: away from its ring's centre). */
  out?: (i: number, j: number, p: THREE.Vector3) => THREE.Vector3;
  /** Override UVs (default: metres round the ring and along the loft). */
  uv?: (i: number, j: number, p: THREE.Vector3, u: number, v: number) => readonly [number, number];
}

const _la = V();
const _lb = V();
const _lc = V();
/**
 * Sweep rings into a skinned surface. Every ring has the same number of points. Normals come from the grid's own
 * derivatives (so seams of closed rings stay smooth), and the winding follows whichever way `out` says is outside.
 */
function loft(soup: Soup, rings: readonly THREE.Vector3[][], o: LoftOpts): void {
  const R = rings.length;
  const n = rings[0]!.length;
  const closed = o.closed !== false;
  const cen = rings.map((r) => r.reduce((a, p) => a.add(p), V()).multiplyScalar(1 / r.length));
  const P = (i: number, j: number) => rings[clamp(i, 0, R - 1)]![closed ? (j + n) % n : clamp(j, 0, n - 1)]!;
  const degenerate = rings.map((r) => r.every((p) => p.distanceToSquared(r[0]!) < 1e-12));
  const N: THREE.Vector3[][] = [];
  for (let i = 0; i < R; i++) {
    const row: THREE.Vector3[] = [];
    for (let j = 0; j < n; j++) {
      if (degenerate[i]) {
        const nb = i === 0 ? cen[1]! : cen[R - 2]!;
        row.push(V().subVectors(rings[i]![0]!, nb).normalize());
        continue;
      }
      _la.subVectors(P(i, j + 1), P(i, j - 1));
      let i0 = i - 1;
      let i1 = i + 1;
      i0 = Math.max(0, i0);
      i1 = Math.min(R - 1, i1);
      _lb.subVectors(P(i1, j), P(i0, j));
      row.push(V().crossVectors(_la, _lb).normalize());
    }
    N.push(row);
  }
  // Orientation, decided by a majority vote over a middle ring.
  const mid = Math.floor(R / 2);
  let vote = 0;
  for (let j = 0; j < n; j++) {
    const p = rings[mid]![j]!;
    const ref = o.out ? o.out(mid, j, p) : _lc.subVectors(p, cen[mid]!);
    vote += Math.sign(N[mid]![j]!.dot(ref));
  }
  const flip = vote < 0;
  if (flip) for (const row of N) for (const nn of row) if (nn.lengthSq() > 0) nn.negate();
  for (let i = 0; i < R; i++)
    if (degenerate[i]) {
      const nb = i === 0 ? cen[1]! : cen[R - 2]!;
      for (const nn of N[i]!) nn.subVectors(rings[i]![0]!, nb).normalize();
    }
  // UVs in metres.
  const vs: number[] = [0];
  for (let i = 1; i < R; i++) {
    let d = 0;
    for (let j = 0; j < n; j++) d += rings[i]![j]!.distanceTo(rings[i - 1]![j]!);
    vs.push(vs[i - 1]! + d / n);
  }
  const cols = n + (closed ? 1 : 0);
  const base = soup.verts;
  for (let i = 0; i < R; i++) {
    let u = 0;
    for (let j = 0; j < cols; j++) {
      const jj = j % n;
      if (j > 0) u += rings[i]![jj]!.distanceTo(rings[i]![(j - 1) % n]!);
      const p = rings[i]![jj]!;
      const c = typeof o.color === "function" ? o.color(i, jj, p) : o.color;
      const [uu, vv] = o.uv ? o.uv(i, jj, p, u, vs[i]!) : [u, vs[i]!];
      soup.vert(p, N[i]![jj]!, uu, vv, c, o.w(i, jj, p));
    }
  }
  for (let i = 0; i < R - 1; i++)
    for (let j = 0; j < cols - 1; j++) {
      const a = base + i * cols + j;
      const b = a + 1;
      const c = a + cols + 1;
      const d = a + cols;
      if (flip) {
        soup.tri(a, c, b);
        soup.tri(a, d, c);
      } else {
        soup.tri(a, b, c);
        soup.tri(a, c, d);
      }
    }
}

const _gm = new THREE.Matrix3();
const _gp = V();
const _gn = V();
/** Add a three.js primitive, transformed by `m`, rigid or with per-point influences. */
function addGeo(soup: Soup, g: THREE.BufferGeometry, m: THREE.Matrix4, w: Wts | WFn, color: THREE.Color | ((p: THREE.Vector3) => THREE.Color), uvScale = 1): void {
  const pos = g.attributes.position!;
  const nrm = g.attributes.normal;
  const uv = g.attributes.uv;
  _gm.getNormalMatrix(m);
  const base = soup.verts;
  for (let i = 0; i < pos.count; i++) {
    _gp.fromBufferAttribute(pos, i).applyMatrix4(m);
    if (nrm) _gn.fromBufferAttribute(nrm, i).applyMatrix3(_gm).normalize();
    else _gn.set(0, 1, 0);
    const c = typeof color === "function" ? color(_gp) : color;
    soup.vert(_gp, _gn, uv ? uv.getX(i) * uvScale : 0, uv ? uv.getY(i) * uvScale : 0, c, typeof w === "function" ? w(_gp) : w);
  }
  const mirrored = m.determinant() < 0;
  const ix = g.index;
  const cnt = ix ? ix.count : pos.count;
  for (let k = 0; k < cnt; k += 3) {
    const a = base + (ix ? ix.getX(k) : k);
    const b = base + (ix ? ix.getX(k + 1) : k + 1);
    const c = base + (ix ? ix.getX(k + 2) : k + 2);
    if (mirrored) soup.tri(a, c, b);
    else soup.tri(a, b, c);
  }
  g.dispose();
}
const mat4 = (p: THREE.Vector3, q?: THREE.Quaternion, s?: THREE.Vector3) => new THREE.Matrix4().compose(p, q ?? new THREE.Quaternion(), s ?? V(1, 1, 1));

/* ------------------------------------------------------------------------------------------------------------ */
/* Skeleton                                                                                                      */
/* ------------------------------------------------------------------------------------------------------------ */

const FINGERS = ["thumb", "index", "middle", "ring", "pinky"] as const;
type Finger = (typeof FINGERS)[number];
const SIDES: readonly Side[] = ["L", "R"];
/** +1 for the left side (at +x), -1 for the right. */
const sgn = (s: Side) => (s === "L" ? 1 : -1);

/** The 29 bone names of the staff rig manifest, kept. */
export const REQUIRED_BONES: readonly string[] = [
  "root", "pelvis", "spine_01", "spine_02", "spine_03", "neck", "head", "gaze_target", "carry_chest",
  "clavicle_L", "clavicle_R", "upper_arm_L", "upper_arm_R", "forearm_L", "forearm_R", "hand_L", "hand_R",
  "grip_L", "grip_R", "claw_grip_L", "pinch_grip_R", "thigh_L", "thigh_R", "shin_L", "shin_R",
  "foot_L", "foot_R", "toe_L", "toe_R",
];
/** The finger bones added to them: three a digit, proximal to distal. */
export const FINGER_BONES: readonly string[] = SIDES.flatMap((s) => FINGERS.flatMap((f) => [1, 2, 3].map((k) => `${f}_0${k}_${s}`)));
/**
 * Two cloth bones: the apron's lower edge hangs from each knee, plumb unless the shin swings forward and pushes
 * it (posed each frame, never animated by hand).
 */
export const CLOTH_BONES: readonly string[] = ["apron_L", "apron_R"];
const ALL_BONES = [...REQUIRED_BONES, ...FINGER_BONES, ...CLOTH_BONES];

function parentOf(b: string): string | null {
  if (b === "root") return null;
  const m = /^(.*)_([LR])$/.exec(b);
  const s = m ? m[2]! : "";
  const base = m ? m[1]! : b;
  const fm = /^(thumb|index|middle|ring|pinky)_0(\d)$/.exec(base);
  if (fm) return fm[2] === "1" ? `hand_${s}` : `${fm[1]}_0${Number(fm[2]) - 1}_${s}`;
  const table: Record<string, string> = {
    pelvis: "root",
    spine_01: "pelvis",
    spine_02: "spine_01",
    spine_03: "spine_02",
    neck: "spine_03",
    head: "neck",
    gaze_target: "head",
    carry_chest: "spine_02",
    clavicle: "spine_03",
    upper_arm: "clavicle",
    forearm: "upper_arm",
    hand: "forearm",
    grip: "hand",
    claw_grip: "hand",
    pinch_grip: "hand",
    thigh: "pelvis",
    shin: "thigh",
    foot: "shin",
    toe: "foot",
    apron: "thigh",
  };
  const p = table[base]!;
  return ["pelvis", "spine_01", "spine_02", "spine_03", "neck", "head", "root"].includes(p) ? p : `${p}_${s}`;
}

interface DigitDef {
  /** Base joint (MCP, or CMC for the thumb), bind body space. */
  base: THREE.Vector3;
  /** Straight direction of the digit at bind. */
  dir: THREE.Vector3;
  /** Flexion axis: a positive turn about it curls the digit towards the palm. */
  k: THREE.Vector3;
  /** Back of the digit (where the nail is). */
  back: THREE.Vector3;
  len: readonly [number, number, number];
  rad: number;
}
interface HandDef {
  side: Side;
  sd: number;
  wrist: THREE.Vector3;
  /** Distal along the hand, palm normal (facing in, towards the body at bind), thumb side. */
  da: THREE.Vector3;
  dn: THREE.Vector3;
  dt: THREE.Vector3;
  HL: number;
  breadth: number;
  sh: number;
  digits: Record<Finger, DigitDef>;
  /** Grip frames relative to the hand bone (bind orientation is the identity): tool (oblique), bar (power), palm. */
  toolPos: THREE.Vector3;
  toolQ: THREE.Quaternion;
  barPos: THREE.Vector3;
  barAxis: THREE.Vector3;
  hookPos: THREE.Vector3;
  palmPos: THREE.Vector3;
}

/** Everything a body's build and motion needs, derived from its row. */
interface Body {
  d: Anthro;
  cohort: Cohort;
  s: number;
  lift: number;
  fem: number;
  J: Record<string, THREE.Vector3>;
  thigh: number;
  shin: number;
  upper: number;
  fore: number;
  footLen: number;
  heelBack: number;
  mtpFwd: number;
  toeFwd: number;
  shoeW: number;
  hipY: number;
  kneeY: number;
  ankY: number;
  shY: number;
  eyeY: number;
  topY: number;
  waistY: number;
  chestY: number;
  armpitY: number;
  neckBaseY: number;
  hemY: number;
  hands: Record<Side, HandDef>;
  /** Stance half-width (ankle centres) standing at work, and walking. */
  standHalf: number;
  walkHalf: number;
}

function digit(h: { wrist: THREE.Vector3; da: THREE.Vector3; dn: THREE.Vector3; dt: THREE.Vector3; sd: number }, a: number, n: number, t: number, fan: number, len: [number, number, number], rad: number): DigitDef {
  const base = h.wrist.clone().addScaledVector(h.da, a).addScaledVector(h.dn, n).addScaledVector(h.dt, t);
  const dir = h.da.clone().multiplyScalar(Math.cos(fan)).addScaledVector(h.dt, Math.sin(fan)).normalize();
  const k = V().crossVectors(dir, h.dn).normalize();
  return { base, dir, k, back: h.dn.clone().negate(), len, rad };
}

function layout(cohort: Cohort): Body {
  const d = ANTHROPOMETRICS[cohort];
  const H = d.standing_height_m;
  const s = H / 1.715;
  const lift = d.shoe_sole_m;
  const fem = cohort === "p05" ? 1 : cohort === "p50" ? 0.4 : 0;
  const hipY = d.hip_height_m + lift;
  const kneeY = d.knee_height_m + lift;
  const ankY = d.ankle_height_m + lift;
  const shY = d.shoulder_height_m + lift;
  const eyeY = d.eye_height_m + lift;
  const topY = H + lift;
  const SB = d.shoulder_breadth_m;
  const hx = 0.215 * SB * (1 + 0.07 * fem);
  const c7y = shY + 0.035 * s;
  const l5y = hipY + 0.085 * s;
  const span = c7y - l5y;
  const J: Record<string, THREE.Vector3> = {};
  J.root = V();
  J.pelvis = V(0, hipY, 0);
  J.spine_01 = V(0, l5y, -0.035 * s);
  J.spine_02 = V(0, l5y + 0.37 * span, -0.045 * s);
  J.spine_03 = V(0, l5y + 0.66 * span, -0.058 * s);
  J.neck = V(0, c7y, -0.066 * s);
  J.head = V(0, eyeY - 0.068 * s, -0.026 * s);
  J.gaze_target = V(0, eyeY, 0.25);
  J.carry_chest = V(0, J.spine_02.y + 0.06 * s, 0.16 * s);
  const footLen = d.shoe_length_m - 0.022;
  const mtpFwd = 0.5 * footLen;
  const hands = {} as Record<Side, HandDef>;
  for (const side of SIDES) {
    const sd = sgn(side);
    J[`clavicle_${side}`] = V(sd * 0.02 * s, shY - 0.03 * s, 0.05 * s);
    const gh = V(sd * (SB / 2 - 0.046 * s), shY - 0.035 * s, -0.012 * s);
    J[`upper_arm_${side}`] = gh;
    const el = gh.clone().add(V(sd * Math.sin(6 * DEG), -Math.cos(6 * DEG), 0).multiplyScalar(d.upper_arm_m));
    J[`forearm_${side}`] = el;
    const da = V(sd * Math.sin(10 * DEG), -Math.cos(10 * DEG), 0);
    const wr = el.clone().addScaledVector(da, d.forearm_m);
    J[`hand_${side}`] = wr;
    J[`thigh_${side}`] = V(sd * hx, hipY, 0);
    J[`shin_${side}`] = V(sd * hx * 0.965, kneeY, 0);
    J[`apron_${side}`] = V(sd * hx * 0.965, kneeY, 0);
    J[`foot_${side}`] = V(sd * hx * 0.93, ankY, 0);
    J[`toe_${side}`] = V(sd * hx * 0.93, lift + 0.022 * s, mtpFwd);
    // The hand, in its own terms: distal, palm normal (medial), thumb side (forward at bind).
    const HL = d.hand_length_m;
    const sh = HL / 0.186;
    const B = 0.45 * HL;
    const dn = V(-sd * Math.cos(10 * DEG), -Math.sin(10 * DEG), 0);
    const dt = V(0, 0, 1);
    const h = { wrist: wr, da, dn, dt, sd };
    const digits: Record<Finger, DigitDef> = {
      index: digit(h, 0.505 * HL, -0.002 * sh, 0.37 * B, 4 * DEG, [0.225 * HL, 0.135 * HL, 0.112 * HL], 0.0093 * sh),
      middle: digit(h, 0.52 * HL, -0.002 * sh, 0.12 * B, 0, [0.245 * HL, 0.155 * HL, 0.116 * HL], 0.0096 * sh),
      ring: digit(h, 0.505 * HL, -0.002 * sh, -0.12 * B, -4 * DEG, [0.23 * HL, 0.148 * HL, 0.114 * HL], 0.009 * sh),
      pinky: digit(h, 0.462 * HL, 0, -0.35 * B, -9 * DEG, [0.182 * HL, 0.108 * HL, 0.104 * HL], 0.0079 * sh),
      thumb: digit(h, 0.13 * HL, 0.006 * sh, 0.36 * B, 0, [0.235 * HL, 0.168 * HL, 0.148 * HL], 0.0112 * sh),
    };
    // The thumb points distally, palmwards and out; it curls across the palm towards the little finger.
    const th = digits.thumb;
    th.dir.copy(da).multiplyScalar(0.8).addScaledVector(dn, 0.38).addScaledVector(dt, 0.46).normalize();
    const curl = orthoTo(dn.clone().multiplyScalar(0.55).addScaledVector(dt, -0.83), th.dir);
    th.k.crossVectors(th.dir, curl).normalize();
    th.back.copy(curl).negate();
    // Oblique grip: the handle crosses the palm from the heel of the hand on the little-finger side to the web
    // of the thumb, 55 degrees off the hand's long axis, so the fingers wrap round it and the thumb lies along it.
    const OBL = 55 * DEG;
    const zg = da.clone().multiplyScalar(Math.cos(OBL)).addScaledVector(dt, Math.sin(OBL)).normalize();
    const yg = orthoTo(da.clone().multiplyScalar(-Math.sin(OBL)).addScaledVector(dt, Math.cos(OBL)), zg);
    const xg = V().crossVectors(yg, zg);
    const toolQ = new THREE.Quaternion().setFromRotationMatrix(new THREE.Matrix4().makeBasis(xg, yg, zg));
    hands[side] = {
      side,
      sd,
      wrist: wr,
      da,
      dn,
      dt,
      HL,
      breadth: B,
      sh,
      digits,
      toolPos: da.clone().multiplyScalar(0.36 * HL).addScaledVector(dn, 0.026 * sh).addScaledVector(dt, 0.004 * sh),
      toolQ,
      barPos: da.clone().multiplyScalar(0.47 * HL).addScaledVector(dn, 0.024 * sh),
      barAxis: da.clone().multiplyScalar(0.22).addScaledVector(dt, 0.975).normalize(),
      hookPos: da.clone().multiplyScalar(0.6 * HL).addScaledVector(dn, 0.018 * sh),
      palmPos: da.clone().multiplyScalar(0.36 * HL).addScaledVector(dn, 0.016 * sh),
    };
    for (const f of FINGERS) {
      const g = digits[f];
      J[`${f}_01_${side}`] = g.base.clone();
      J[`${f}_02_${side}`] = g.base.clone().addScaledVector(g.dir, g.len[0]);
      J[`${f}_03_${side}`] = g.base.clone().addScaledVector(g.dir, g.len[0] + g.len[1]);
    }
    J[`grip_${side}`] = wr.clone().add(hands[side].toolPos);
  }
  J.claw_grip_L = J.hand_L!.clone().addScaledVector(hands.L.da, 0.62 * hands.L.HL).addScaledVector(hands.L.dn, 0.012);
  J.pinch_grip_R = J.hand_R!.clone().addScaledVector(hands.R.da, 0.78 * hands.R.HL).addScaledVector(hands.R.dn, 0.02).addScaledVector(hands.R.dt, 0.03);
  return {
    d,
    cohort,
    s,
    lift,
    fem,
    J,
    thigh: hipY - kneeY,
    shin: kneeY - ankY,
    upper: d.upper_arm_m,
    fore: d.forearm_m,
    footLen,
    heelBack: 0.21 * footLen,
    mtpFwd,
    toeFwd: 0.79 * footLen,
    shoeW: d.shoe_width_m,
    hipY,
    kneeY,
    ankY,
    shY,
    eyeY,
    topY,
    waistY: d.elbow_height_m + lift - 0.01 * s,
    chestY: shY - 0.135 * s,
    armpitY: shY - 0.09 * s,
    neckBaseY: c7y + 0.012 * s,
    hemY: hipY - 0.085 * s,
    hands,
    standHalf: hx * 1.02,
    walkHalf: Math.max(0.5 * d.shoe_width_m + 0.016, hx * 0.62),
  };
}


const BI: Record<string, number> = Object.fromEntries(ALL_BONES.map((b, i) => [b, i]));
const rigid = (bone: string): Wts => [[BI[bone]!, 1]];

/* ------------------------------------------------------------------------------------------------------------ */
/* The torso: the jacket's outer surface as horizontal superellipse sections, and how it follows the spine.       */
/* ------------------------------------------------------------------------------------------------------------ */

interface Sect {
  a: number;
  f: number;
  b: number;
  c: number;
  n: number;
}
interface Torso {
  sect(y: number): Sect;
  /** Point on the jacket surface at height y, angle th (0 front, +pi/2 the person's left). */
  at(y: number, th: number, out?: THREE.Vector3): THREE.Vector3;
  /** Outward horizontal normal at that point (approximate). */
  normal(y: number, th: number, out?: THREE.Vector3): THREE.Vector3;
  /** Angle at which the surface reaches x on the front, at height y. */
  angleAtX(y: number, x: number): number;
  hipHalf: number;
}

function torsoModel(B: Body): Torso {
  const s = B.s;
  const sb = B.d.shoulder_breadth_m / 0.448;
  const cd = B.d.chest_depth_m / 0.238;
  const fem = B.fem;
  const hipHalf = 0.188 * s * (1 + 0.06 * fem);
  const waistHalf = 0.158 * s * (1 - 0.04 * fem);
  const keys: [number, number, number, number, number, number][] = [
    [B.hemY - 0.01, hipHalf * 1.04, 0.124 * cd, (0.138 + 0.006 * fem) * cd, 0, 2.2],
    [B.hipY, hipHalf, 0.116 * cd, (0.134 + 0.008 * fem) * cd, 0, 2.3],
    [(B.hipY + B.waistY) / 2, lerp(hipHalf, waistHalf, 0.55), 0.112 * cd, 0.118 * cd, 0.002, 2.35],
    [B.waistY, waistHalf, 0.111 * cd, 0.103 * cd, 0.004, 2.4],
    [B.chestY - 0.075 * s, 0.158 * sb, 0.12 * cd, 0.104 * cd, 0.002, 2.5],
    [B.chestY, 0.163 * sb, (0.132 + 0.01 * fem) * cd, 0.109 * cd, 0, 2.55],
    [B.armpitY, 0.16 * sb, 0.127 * cd, 0.111 * cd, -0.004, 2.6],
    [B.shY - 0.035 * s, 0.166 * sb, 0.106 * cd, 0.098 * cd, -0.012, 2.7],
    [B.shY + 0.004 * s, 0.158 * sb, 0.084 * cd, 0.086 * cd, -0.02, 2.6],
    [B.shY + 0.026 * s, 0.112 * sb, 0.07 * s, 0.074 * s, -0.03 * s, 2.3],
    [B.neckBaseY, 0.066 * s, 0.056 * s, 0.062 * s, -0.035 * s, 2.1],
  ];
  const ys = keys.map((k) => k[0]);
  const sp = [1, 2, 3, 4, 5].map((c) => new Mono(ys, keys.map((k) => k[c]!)));
  const sect = (y: number): Sect => ({ a: sp[0]!.at(y), f: sp[1]!.at(y), b: sp[2]!.at(y), c: sp[3]!.at(y), n: sp[4]!.at(y) });
  const at = (y: number, th: number, out = V()) => {
    const S = sect(y);
    const e = 2 / S.n;
    const cs = Math.cos(th);
    return out.set(S.a * spow(Math.sin(th), e), y, S.c + (cs >= 0 ? S.f : S.b) * spow(cs, e));
  };
  const _p = V();
  const _q = V();
  const normal = (y: number, th: number, out = V()) => {
    at(y, th - 0.01, _p);
    at(y, th + 0.01, _q);
    _q.sub(_p);
    return out.set(-_q.z, 0, _q.x).normalize();
  };
  const angleAtX = (y: number, x: number) => {
    const S = sect(y);
    const r = clamp(Math.abs(x) / S.a, 0, 0.999);
    return Math.sign(x) * Math.asin(Math.pow(r, S.n / 2));
  };
  return { sect, at, normal, angleAtX, hipHalf };
}

/** How the trunk's skin and cloth follow the spine, the clavicles and (below the hips) the thighs. */
/**
 * How much of a garment over the hips the thighs carry at height y: the same rule for the trousers, the jacket's
 * hem and the apron, so the layers move together and the inner never comes through the outer.
 */
function legShare(B: Body, y: number): number {
  return smooth(-0.05 * B.s, 0.12 * B.s, B.hipY - y);
}

function trunkW(B: Body, p: THREE.Vector3, thighShare = 0.5): Wts {
  const J = B.J;
  const s = B.s;
  const y = p.y;
  const band = 0.045 * s;
  const w = new Map<number, number>();
  const add = (b: string, x: number) => w.set(BI[b]!, (w.get(BI[b]!) ?? 0) + x);
  const cuts = [J.spine_01!.y - 0.015 * s, J.spine_02!.y, J.spine_03!.y];
  const chain = ["pelvis", "spine_01", "spine_02", "spine_03"];
  // Piecewise blend along the chain.
  let k = 0;
  while (k < 3 && y > cuts[k]! + band) k++;
  let spine: [string, number][];
  if (k < 3 && y > cuts[k]! - band) {
    const t = smooth(cuts[k]! - band, cuts[k]! + band, y);
    spine = [[chain[k]!, 1 - t], [chain[k + 1]!, t]];
  } else spine = [[chain[k]!, 1]];
  // Shoulders ride on the clavicles.
  const ax = Math.abs(p.x);
  const cw = 0.85 * smooth(0.06 * s, 0.15 * s, ax) * smooth(B.armpitY - 0.03 * s, B.shY - 0.015 * s, y);
  // Below the hip joints the hem is carried partly by the thighs.
  const tw = Math.min(1, 2 * thighShare) * legShare(B, y);
  const rest = (1 - cw) * (1 - tw);
  for (const [b, x] of spine) add(b, x * rest);
  if (cw > 0) add(p.x > 0 ? "clavicle_L" : "clavicle_R", cw * (1 - tw));
  if (tw > 0) {
    const wl = smooth(-0.05 * s, 0.05 * s, p.x);
    add("thigh_L", tw * wl);
    add("thigh_R", tw * (1 - wl));
  }
  return [...w.entries()];
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Chef's jacket                                                                                                  */
/* ------------------------------------------------------------------------------------------------------------ */

interface Look {
  skin: THREE.Color;
  skinDark: THREE.Color;
  nail: THREE.Color;
  hair: THREE.Color;
  white: THREE.Color;
  button: THREE.Color;
  apron: THREE.Color;
  apronTie: THREE.Color;
  trousers: THREE.Color;
  upper: THREE.Color;
  sole: THREE.Color;
  cap: THREE.Color;
  sleeves: "rolled" | "cuffed";
  shoes: "clog" | "shoe";
  capKind: "skull" | "toque";
  check: boolean;
  bun: boolean;
}

/** The flap edge of the double-breasted front, x on the person's right side, as a function of height. */
function flapEdgeX(B: Body, y: number): number {
  return -lerp(0.088 * B.s, 0.03 * B.s, smooth(B.chestY - 0.02 * B.s, B.neckBaseY, y));
}

function buildJacket(soup: Soup, B: Body, T: Torso, L: Look, detail: Detail): void {
  const s = B.s;
  const nAround = detail === 2 ? 40 : detail === 1 ? 24 : 14;
  const y0 = B.hemY;
  const y1 = B.neckBaseY;
  // Ring heights: even down the body, closer over the shoulders where the section changes fast.
  const ys: number[] = [];
  const nBody = detail === 2 ? 18 : detail === 1 ? 14 : 7;
  const ySh = B.armpitY;
  for (let i = 0; i <= nBody; i++) ys.push(lerp(y0, ySh, i / nBody));
  const nTop = detail === 2 ? 8 : detail === 1 ? 6 : 3;
  for (let i = 1; i <= nTop; i++) ys.push(lerp(ySh, y1, Math.sin((i / nTop) * Math.PI * 0.5)));
  const lift = 0.0035 * s;
  const rings: THREE.Vector3[][] = [];
  for (const y of ys) {
    const xe = flapEdgeX(B, y);
    const te = T.angleAtX(y, xe);
    const d = 0.0025;
    const ring: THREE.Vector3[] = [];
    for (let j = 0; j < nAround; j++) {
      const th = te + d + ((Math.PI * 2 - 2 * d) * j) / (nAround - 1);
      const p = T.at(y, th);
      const tw = wrapPi(th);
      const flap = j === nAround - 1 ? 0 : tw > te - 1e-6 && tw < 1.8 ? 1 - smooth(1.15, 1.8, tw) : 0;
      if (flap > 0) p.add(T.normal(y, th).multiplyScalar(lift * flap));
      ring.push(p);
    }
    rings.push(ring);
  }
  loft(soup, rings, { w: (_i, _j, p) => trunkW(B, p), color: L.white });

  // Stand collar: a band round the neck base, open a little at the front, leaning in towards the neck.
  const cH = 0.042 * s;
  const cR0 = 0.07 * s;
  const cz = -0.03 * s;
  const nC = detail === 2 ? 24 : detail === 1 ? 18 : 8;
  const cRings: THREE.Vector3[][] = [];
  const open = 0.22;
  const prof = detail > 0 ? [[0, 0], [0.003, 0.15], [0.004, 0.85], [0.0015, 1], [-0.002, 1], [-0.0035, 0.2]] : [[0, 0], [0.003, 1], [-0.003, 1]];
  for (let k = 0; k <= nC; k++) {
    const th = open + ((Math.PI * 2 - 2 * open) * k) / nC;
    const ring: THREE.Vector3[] = [];
    for (const [dr, h] of prof as [number, number][]) {
      const y = B.neckBaseY - 0.006 * s + h * cH;
      const r = cR0 * (1 - 0.1 * h) + dr;
      const front = Math.max(0, Math.cos(th));
      ring.push(V(Math.sin(th) * r * 1.02, y - front * 0.012 * s, cz + Math.cos(th) * r * (1 + 0.12 * front)));
    }
    cRings.push(ring);
  }
  const neckW = (p: THREE.Vector3): Wts => {
    const t = smooth(B.neckBaseY, B.neckBaseY + cH, p.y) * 0.45;
    return [[BI.spine_03!, 1 - t], [BI.neck!, t]];
  };
  loft(soup, cRings, { closed: false, w: (_i, _j, p) => neckW(p), color: L.white, out: (_i, _j, p) => V(p.x, 0, p.z - cz) });

  // Knotted cloth buttons: two columns down the double-breasted front.
  const rows = detail === 0 ? 0 : 5;
  // The top rows sit between the collar and the apron's bib, the rest under it.
  const top = B.neckBaseY - 0.058 * s;
  for (let r = 0; r < rows; r++) {
    const y = top - r * 0.072 * s;
    // Rows under the apron's bib would only push through it; the cloth hides them anyway.
    if (y < B.chestY + 0.05 * s + 0.012 * s) continue;
    for (const cx of [-0.058 * s, 0.058 * s]) {
      const th = T.angleAtX(y, cx);
      const p = T.at(y, th);
      const n = T.normal(y, th);
      p.addScaledVector(n, lift + 0.0015);
      knotButton(soup, p, n, 0.0068 * s, trunkW(B, p), L.button, detail);
    }
  }

  // Sleeves.
  for (const side of SIDES) sleeve(soup, B, side, L, detail);
}

/** A Chinese knot button: a ball of braided cord with the loops of its knot round it. */
function knotButton(soup: Soup, p: THREE.Vector3, n: THREE.Vector3, r: number, w: Wts, c: THREE.Color, detail: Detail): void {
  const q = new THREE.Quaternion().setFromUnitVectors(Y, n);
  const seg = detail === 2 ? 8 : 6;
  addGeo(soup, new THREE.SphereGeometry(r, seg, Math.max(4, seg - 4)), mat4(p.clone().addScaledVector(n, r * 0.75), q, V(1, 0.82, 1)), w, c);
  if (detail === 2) {
    // The cord crossing over the ball and the loops of the knot.
    for (const [a, tilt] of [[0, 0.5], [Math.PI / 2, -0.5], [Math.PI / 4, 0]] as [number, number][]) {
      const g = new THREE.TorusGeometry(r * 0.78, r * 0.2, 5, 14);
      const qq = q.clone().multiply(new THREE.Quaternion().setFromEuler(new THREE.Euler(Math.PI / 2 + tilt, a, 0)));
      addGeo(soup, g, mat4(p.clone().addScaledVector(n, r * 0.7), qq, V(1, 1, 0.9)), w, c);
    }
  }
}

/** A sleeve from its cap in the shoulder to a roll below the elbow or a turned-back cuff at the wrist. */
function sleeve(soup: Soup, B: Body, side: Side, L: Look, detail: Detail): void {
  const s = B.s;
  const sd = sgn(side);
  const gh = B.J[`upper_arm_${side}`]!;
  const el = B.J[`forearm_${side}`]!;
  const wr = B.J[`hand_${side}`]!;
  const up = V().subVectors(el, gh).normalize();
  const fo = V().subVectors(wr, el).normalize();
  const lu = B.upper;
  const lf = B.fore;
  const rolled = L.sleeves === "rolled";
  const tEnd = rolled ? lu + 0.07 * s : lu + lf - 0.004 * s;
  // The cloth over the deltoid is fullest just below the shoulder point, easing to the elbow.
  const rAt = (t: number) => {
    if (t < lu) return lerp(0.0565, 0.048, smooth(0.03 * s, lu, t)) * s;
    return lerp(0.049, rolled ? 0.046 : 0.041, clamp((t - lu) / lf, 0, 1)) * s;
  };
  // Stations along the arm (t from the shoulder joint; negative is the cap above it).
  const ts: number[] = [];
  const capTop = -0.05 * s;
  const capN = detail === 2 ? 7 : detail === 1 ? 4 : 2;
  for (let i = 0; i < capN; i++) ts.push(capTop * Math.cos(((i / capN) * Math.PI) / 2));
  const nU = detail === 2 ? 8 : detail === 1 ? 5 : 3;
  const tRollStart = tEnd - (rolled ? 0.05 : 0.07) * s;
  for (let i = 0; i <= nU; i++) ts.push(lerp(0, tRollStart, i / nU));
  // The roll (or the cuff): a fat fold that tucks in round the forearm.
  const tail: [number, number][] = rolled
    ? [[0.008, 0.004], [0.016, 0.011], [0.026, 0.015], [0.036, 0.014], [0.045, 0.009], [0.05, 0.001], [0.051, -0.006]]
    : [[0.0, 0.0045], [0.035, 0.005], [0.066, 0.004], [0.07, -0.002], [0.071, -0.007]];
  const extra: number[] = [];
  for (const [dt, dr] of tail) {
    ts.push(tRollStart + dt * s);
    extra.push(dr * s);
  }
  const nA = detail === 2 ? 14 : detail === 1 ? 10 : 6;
  const rings: THREE.Vector3[][] = [];
  const nTail = tail.length;
  for (let i = 0; i < ts.length; i++) {
    const t = ts[i]!;
    const iTail = i - (ts.length - nTail);
    const c = t <= lu ? gh.clone().addScaledVector(up, t) : el.clone().addScaledVector(fo, t - lu);
    const dir = t < lu - 0.04 * s ? up : t > lu + 0.04 * s ? fo : up.clone().lerp(fo, smooth(lu - 0.04 * s, lu + 0.04 * s, t)).normalize();
    const u = orthoTo(V(sd, 0, 0), dir);
    const v = V().crossVectors(dir, u).multiplyScalar(sd);
    let r = rAt(Math.max(0, t));
    let rz = 1.06;
    if (t < 0) {
      // The cap: the deltoid's dome, its crown drawn in towards the neck so that the jacket's shoulder line runs
      // on over it into the sleeve as one curve (its inner half lies inside the jacket's shoulder, unseen).
      const k = Math.sqrt(Math.max(0, 1 - (t / capTop) ** 2));
      r = rAt(0) * (0.25 + 0.75 * k);
      rz = 1.1;
      c.add(V(-sd * 0.022 * s * (1 - k), 0, 0));
    }
    if (iTail >= 0) r += extra[iTail]!;
    const ring: THREE.Vector3[] = [];
    for (let j = 0; j < nA; j++) {
      const a = (j / nA) * Math.PI * 2;
      ring.push(c.clone().addScaledVector(u, Math.cos(a) * r).addScaledVector(v, Math.sin(a) * r * rz));
    }
    rings.push(ring);
  }
  const iClav = BI[`clavicle_${side}`]!;
  const iUp = BI[`upper_arm_${side}`]!;
  const iFo = BI[`forearm_${side}`]!;
  loft(soup, rings, {
    w: (i, _j, p) => {
      const t = ts[i]!;
      if (t < 0.03 * s) {
        // The cap's inner half (towards the neck) rides more with the shoulder girdle than with the arm.
        const k = smooth(capTop, 0.03 * s, t);
        const inner = smooth(0, -0.03 * s, sd * (p.x - gh.x));
        const wc = (0.3 + 0.35 * inner) * (1 - k);
        return [[iClav, wc], [iUp, 1 - wc]];
      }
      const e = smooth(lu - 0.035 * s, lu + 0.035 * s, t);
      return [[iUp, 1 - e], [iFo, e]];
    },
    color: L.white,
    out: (i, _j, p) => {
      const t = ts[i]!;
      const c = t <= lu ? gh.clone().addScaledVector(up, t) : el.clone().addScaledVector(fo, t - lu);
      return p.clone().sub(c);
    },
  });
}


/* ------------------------------------------------------------------------------------------------------------ */
/* Trousers                                                                                                       */
/* ------------------------------------------------------------------------------------------------------------ */

function buildTrousers(soup: Soup, B: Body, T: Torso, L: Look, detail: Detail): void {
  const s = B.s;
  const nA = detail === 2 ? 18 : detail === 1 ? 14 : 8;
  // The seat: from the waist (under the jacket) down to the crotch, closing between the legs.
  const yTop = B.waistY - 0.03 * s;
  const yC = B.hipY - 0.078 * s;
  const seatYs: number[] = [];
  const nS = detail === 2 ? 8 : 4;
  for (let i = 0; i <= nS; i++) seatYs.push(lerp(yTop, yC + 0.02 * s, i / nS));
  const close: [number, number, number][] = [[yC + 0.008 * s, 0.86, 0.8], [yC, 0.62, 0.55], [yC - 0.004 * s, 0.3, 0.25]];
  const seat: THREE.Vector3[][] = [];
  const ring = (y: number, kx: number, kz: number) => {
    const r: THREE.Vector3[] = [];
    for (let j = 0; j < nA; j++) {
      const p = T.at(Math.max(y, B.hemY), (j / nA) * Math.PI * 2 + Math.PI);
      r.push(V(p.x * 0.955 * kx, y, p.z * 0.95 * kz));
    }
    return r;
  };
  for (const y of seatYs) seat.push(ring(y, 1, 1));
  for (const [y, kx, kz] of close) seat.push(ring(y, kx, kz));
  const seatW = (p: THREE.Vector3): Wts => {
    const tw = legShare(B, p.y);
    const wl = smooth(-0.05 * s, 0.05 * s, p.x);
    return [[BI.pelvis!, 1 - tw], [BI.thigh_L!, tw * wl], [BI.thigh_R!, tw * (1 - wl)]];
  };
  loft(soup, seat, { w: (_i, _j, p) => seatW(p), color: L.trousers });

  for (const side of SIDES) {
    const sd = sgn(side);
    const hip = B.J[`thigh_${side}`]!;
    const knee = B.J[`shin_${side}`]!;
    const ank = B.J[`foot_${side}`]!;
    const th = V().subVectors(knee, hip).normalize();
    const sh = V().subVectors(ank, knee).normalize();
    const lt = B.thigh;
    const ls = B.shin;
    const tHem = lt + ls - (B.ankY - B.lift - 0.034 * s);
    const rAt = (t: number) => {
      const k: [number, number][] = [[-0.08, 0.07], [0.0, 0.085], [0.06, 0.097], [0.16, 0.09], [lt - 0.06, 0.071], [lt, 0.066], [lt + 0.12, 0.066], [lt + ls * 0.75, 0.064], [tHem, 0.069]];
      for (let i = 1; i < k.length; i++) if (t <= k[i]![0]) return lerp(k[i - 1]![1], k[i]![1], smooth(k[i - 1]![0], k[i]![0], t)) * s;
      return k[k.length - 1]![1] * s;
    };
    const ts: number[] = [];
    const nT = detail === 2 ? 16 : detail === 1 ? 12 : 6;
    for (let i = 0; i <= nT; i++) ts.push(lerp(-0.08 * s, tHem, i / nT));
    // Extra rings round the knee so the bend stays round.
    if (detail > 0) for (const dt of [-0.03, -0.012, 0.012, 0.03]) ts.push(lt + dt * s);
    ts.sort((a, b) => a - b);
    const rings: THREE.Vector3[][] = [];
    const lip = detail > 0 ? 1 : 0;
    const all = [...ts];
    if (lip) all.push(tHem + 0.0005);
    for (let i = 0; i < all.length; i++) {
      const t = all[i]!;
      const isLip = lip && i === all.length - 1;
      const c = t <= lt ? hip.clone().addScaledVector(th, t) : knee.clone().addScaledVector(sh, t - lt);
      const dir = t < lt - 0.05 * s ? th : t > lt + 0.05 * s ? sh : th.clone().lerp(sh, smooth(lt - 0.05 * s, lt + 0.05 * s, t)).normalize();
      const u = orthoTo(V(sd, 0, 0), dir);
      const v = V().crossVectors(u, dir).multiplyScalar(-sd);
      const fwd = orthoTo(Z.clone(), dir);
      const r = rAt(Math.min(t, tHem)) - (isLip ? 0.004 * s : 0);
      const near = smooth(tHem - 0.08 * s, tHem, t);
      const rr: THREE.Vector3[] = [];
      for (let j = 0; j < nA; j++) {
        // Start on the inside seam, so the pattern's join is hidden between the legs.
        const a = (j / nA) * Math.PI * 2 + Math.PI;
        const p = c.clone().addScaledVector(u, Math.cos(a) * r * 0.97).addScaledVector(v, Math.sin(a) * r * 1.04);
        // The hem breaks over the instep: higher at the front.
        const front = Math.max(0, V().subVectors(p, c).normalize().dot(fwd));
        p.addScaledVector(dir, -near * front * front * 0.03 * s);
        rr.push(p);
      }
      rings.push(rr);
    }
    const iTh = BI[`thigh_${side}`]!;
    const iSh = BI[`shin_${side}`]!;
    const iFt = BI[`foot_${side}`]!;
    loft(soup, rings, {
      w: (i, _j, p) => {
        const t = all[i]!;
        if (t < 0.12 * s) {
          // legShare's rule, measured along the thigh.
          const k = smooth(-0.05 * s, 0.12 * s, t);
          return [[BI.pelvis!, 1 - k], [iTh, k]];
        }
        const e = smooth(lt - 0.055 * s, lt + 0.055 * s, t);
        // The hem's front rides on the shoe.
        const front = clamp((p.z - ank.z) / (0.06 * s), 0, 1) * smooth(tHem - 0.05 * s, tHem, t) * 0.45;
        return [[iTh, 1 - e], [iSh, e * (1 - front)], [iFt, e * front]];
      },
      color: L.trousers,
      out: (i, _j, p) => {
        const t = all[i]!;
        const c = t <= lt ? hip.clone().addScaledVector(th, t) : knee.clone().addScaledVector(sh, t - lt);
        return p.clone().sub(c);
      },
    });
  }
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Bib apron, neck strap and waist ties                                                                           */
/* ------------------------------------------------------------------------------------------------------------ */

/** A flat tape (an apron tie) along a path, lying on surfaces whose normals are `nrm`. */
function tape(soup: Soup, pts: THREE.Vector3[], nrm: THREE.Vector3[], width: number, thick: number, w: WFn, c: THREE.Color, detail: Detail): void {
  const rings: THREE.Vector3[][] = [];
  const n = detail === 2 ? 6 : 4;
  for (let i = 0; i < pts.length; i++) {
    const p = pts[i]!;
    const tan = V().subVectors(pts[Math.min(i + 1, pts.length - 1)]!, pts[Math.max(i - 1, 0)]!).normalize();
    const nn = orthoTo(nrm[i]!.clone(), tan);
    const side = V().crossVectors(tan, nn);
    const r: THREE.Vector3[] = [];
    for (let j = 0; j < n; j++) {
      const a = (j / n) * Math.PI * 2;
      r.push(p.clone().addScaledVector(side, (Math.cos(a) * width) / 2).addScaledVector(nn, (spow(Math.sin(a), 0.5) * thick) / 2));
    }
    rings.push(r);
  }
  // Close the ends.
  rings.unshift(rings[0]!.map(() => pts[0]!.clone()));
  rings.push(rings[rings.length - 1]!.map(() => pts[pts.length - 1]!.clone()));
  loft(soup, rings, { w: (_i, _j, p) => w(p), color: c });
}

function buildApron(soup: Soup, B: Body, T: Torso, L: Look, detail: Detail): void {
  const s = B.s;
  const off = 0.0068 * s;
  const yTop = B.chestY + 0.05 * s;
  const yBot = B.kneeY - 0.07 * s;
  const yHang = B.hipY - 0.01 * s;
  const bib = 0.122 * s;
  const wrap = 1.72;
  const thHalf = (y: number) => {
    if (y <= B.waistY) return wrap;
    const tb = T.angleAtX(y, bib);
    return lerp(wrap, tb, smooth(B.waistY, B.waistY + 0.17 * s, y));
  };
  const nV = detail === 2 ? 20 : detail === 1 ? 16 : 7;
  const nU = detail === 2 ? 20 : detail === 1 ? 16 : 8;
  const rows: number[] = [];
  for (let i = 0; i <= nV; i++) rows.push(lerp(yTop, yBot, i / nV));
  // A little more ease over the hips, where the thighs pull the layers apart in a stride.
  const offAt = (y: number) => off + 0.003 * s * smooth(B.waistY, B.hipY, y);
  const pt = (y: number, th: number) => {
    if (y >= yHang) return T.at(y, th).addScaledVector(T.normal(y, th), offAt(y));
    const p = T.at(yHang, th).addScaledVector(T.normal(yHang, th), offAt(yHang));
    // Over the jacket's hem the skirt must stay outside the jacket, which widens a little below the hips.
    if (y >= B.hemY - 0.01 * s) {
      const q = T.at(Math.max(y, B.hemY), th).addScaledVector(T.normal(Math.max(y, B.hemY), th), offAt(B.hipY));
      if (q.x * q.x + q.z * q.z > p.x * p.x + p.z * p.z) p.set(q.x, p.y, q.z);
    }
    const flare = 1 + 0.06 * ((yHang - y) / (yHang - yBot));
    return V(p.x * flare, y, p.z * (1 + 0.5 * (flare - 1)) + 0.004 * s * Math.cos(th));
  };
  const grid: THREE.Vector3[][] = [];
  for (const y of rows) {
    const h = thHalf(y);
    const r: THREE.Vector3[] = [];
    for (let j = 0; j <= nU; j++) r.push(pt(y, h * ((2 * j) / nU - 1)));
    grid.push(r);
  }
  const apronW = (p: THREE.Vector3): Wts => {
    if (p.y >= B.hemY) return trunkW(B, p, 0.5);
    // Below the hips the skirt is carried by the thighs, so a forward lean or a crouch lays it over them
    // instead of swinging it back into them; the middle is shared, so a stride stretches it between the legs as
    // cloth would. Below the knee it hangs from the cloth bones.
    const tw = lerp(legShare(B, B.hemY), 0.97, smooth(B.hemY, B.hemY - 0.1 * s, p.y));
    const ks = smooth(B.kneeY + 0.05 * s, B.kneeY - 0.045 * s, p.y);
    const wl = smooth(-0.05 * s, 0.05 * s, p.x);
    return [
      [BI.pelvis!, 1 - tw],
      [BI.thigh_L!, tw * wl * (1 - ks)],
      [BI.thigh_R!, tw * (1 - wl) * (1 - ks)],
      [BI.apron_L!, tw * wl * ks],
      [BI.apron_R!, tw * (1 - wl) * ks],
    ];
  };
  loft(soup, grid, { closed: false, w: (_i, _j, p) => apronW(p), color: L.apron, out: (_i, _j, p) => V(p.x, 0, p.z) });
  if (detail === 0) return;

  // Neck strap: from the bib's top corners up over the shoulders and round the back of the collar.
  const tw = 0.022 * s;
  const tt = 0.0022 * s;
  const strapPts: THREE.Vector3[] = [];
  const strapN: THREE.Vector3[] = [];
  const thC = T.angleAtX(yTop, bib);
  const ySide = B.neckBaseY + 0.012 * s;
  const nUp = detail === 2 ? 12 : 6;
  for (let k = 0; k <= nUp; k++) {
    const u = k / nUp;
    const y = lerp(yTop - 0.006 * s, ySide, ease(u));
    const th = lerp(thC, 1.12, u);
    strapPts.push(T.at(y, th).addScaledVector(T.normal(y, th), off + 0.0025 * s));
    strapN.push(T.normal(y, th).add(V(0, 0.4 * u, 0)).normalize());
  }
  const rC = 0.079 * s;
  const cz = -0.03 * s;
  const nBack = detail === 2 ? 18 : 8;
  const a0 = Math.atan2(strapPts[nUp]!.x, strapPts[nUp]!.z - cz);
  for (let k = 1; k <= nBack; k++) {
    const a = lerp(a0, Math.PI * 2 - a0, k / nBack);
    strapPts.push(V(Math.sin(a) * rC, ySide + 0.006 * s * Math.sin((k / nBack) * Math.PI), cz + Math.cos(a) * rC));
    strapN.push(V(Math.sin(a), 0.3, Math.cos(a)).normalize());
  }
  for (let k = nUp; k >= 0; k--) {
    const p = strapPts[k]!;
    strapPts.push(V(-p.x, p.y, p.z));
    const n = strapN[k]!;
    strapN.push(V(-n.x, n.y, n.z));
  }
  tape(soup, strapPts, strapN, tw, tt, (p) => trunkW(B, p), L.apronTie, detail);

  // Waist ties round to the back, tied in a bow over the spine.
  const yT = B.waistY - 0.012 * s;
  for (const sd of [1, -1]) {
    const pts: THREE.Vector3[] = [];
    const ns: THREE.Vector3[] = [];
    const n = detail === 2 ? 16 : 8;
    for (let k = 0; k <= n; k++) {
      const th = sd * lerp(wrap - 0.05, Math.PI - 0.06, k / n);
      pts.push(T.at(yT, th).addScaledVector(T.normal(yT, th), off + 0.002 * s));
      ns.push(T.normal(yT, th));
    }
    tape(soup, pts, ns, tw * 1.05, tt, (p) => trunkW(B, p), L.apronTie, detail);
  }
  const knot = T.at(yT, Math.PI).addScaledVector(T.normal(yT, Math.PI), off + 0.006 * s);
  const back = T.normal(yT, Math.PI);
  const kw = trunkW(B, knot);
  addGeo(soup, new THREE.SphereGeometry(0.013 * s, 10, 8), mat4(knot, undefined, V(1.2, 0.9, 0.6)), kw, L.apronTie);
  // Two loops of the bow and two hanging tails.
  for (const sd of [1, -1]) {
    const loop: THREE.Vector3[] = [];
    const ln: THREE.Vector3[] = [];
    const n = detail === 2 ? 14 : 7;
    for (let k = 0; k <= n; k++) {
      const a = (k / n) * Math.PI * 2;
      const lx = sd * (0.045 * s) * (1 - Math.cos(a)) * 0.5;
      const ly = 0.018 * s * Math.sin(a);
      loop.push(knot.clone().add(V(lx, ly + 0.004 * s, 0)).addScaledVector(back, 0.004 * s * Math.sin(a / 2)));
      ln.push(back.clone());
    }
    tape(soup, loop, ln, tw * 0.95, tt, (p) => trunkW(B, p), L.apronTie, detail);
    const tail: THREE.Vector3[] = [];
    const tn: THREE.Vector3[] = [];
    for (let k = 0; k <= 8; k++) {
      const u = k / 8;
      tail.push(knot.clone().add(V(sd * (0.012 + 0.03 * u) * s, -0.2 * s * u, 0)).addScaledVector(back, 0.003 * s + 0.01 * s * u * u));
      tn.push(back.clone());
    }
    tape(soup, tail, tn, tw, tt, (p) => trunkW(B, p, 0.5), L.apronTie, detail);
  }
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Shoes: kitchen clogs or lace-up work shoes on non-slip soles                                                   */
/* ------------------------------------------------------------------------------------------------------------ */

function buildShoes(soup: Soup, B: Body, L: Look, detail: Detail): void {
  const s = B.s;
  const W = B.shoeW;
  const clog = L.shoes === "clog";
  const z0 = -B.heelBack - 0.013 * s;
  const z1 = B.toeFwd + 0.015 * s;
  const sole = B.lift;
  const nZ = detail === 2 ? 16 : detail === 1 ? 9 : 4;
  const nA = detail === 2 ? 16 : detail === 1 ? 9 : 6;
  const width = (z: number) => {
    const u = (z - z0) / (z1 - z0);
    const k: [number, number][] = [[0, 0.6], [0.12, 0.7], [0.4, 0.84], [0.66, 1], [0.86, clog ? 0.96 : 0.9], [1, clog ? 0.62 : 0.48]];
    for (let i = 1; i < k.length; i++) if (u <= k[i]![0]) return lerp(k[i - 1]![1], k[i]![1], smooth(k[i - 1]![0], k[i]![0], u)) * W;
    return k[k.length - 1]![1] * W;
  };
  const height = (z: number) => {
    const u = (z - z0) / (z1 - z0);
    const k: [number, number][] = clog
      ? [[0, 0.05], [0.18, 0.062], [0.38, 0.074], [0.55, 0.06], [0.68, 0.05], [0.88, 0.042], [1, 0.022]]
      : [[0, 0.052], [0.18, 0.058], [0.38, 0.066], [0.55, 0.052], [0.68, 0.042], [0.88, 0.032], [1, 0.016]];
    for (let i = 1; i < k.length; i++) if (u <= k[i]![0]) return lerp(k[i - 1]![1], k[i]![1], smooth(k[i - 1]![0], k[i]![0], u)) * s;
    return k[k.length - 1]![1] * s;
  };
  // The toe spring: the sole's front curves up off the floor.
  const spring = (z: number) => 0.009 * s * smooth(B.mtpFwd + 0.01 * s, z1, z) ** 1.6;
  const zs: number[] = [];
  for (let i = 0; i <= nZ; i++) zs.push(lerp(z0, z1, 0.5 - 0.5 * Math.cos((i / nZ) * Math.PI)));
  for (const side of SIDES) {
    const ank = B.J[`foot_${side}`]!;
    const iF = BI[`foot_${side}`]!;
    const iT = BI[`toe_${side}`]!;
    const wt = (p: THREE.Vector3): Wts => {
      const e = smooth(B.mtpFwd - 0.018 * s, B.mtpFwd + 0.018 * s, p.z);
      return [[iF, 1 - e], [iT, e]];
    };
    const at = (x: number, y: number, z: number) => V(ank.x + x, y, z);
    // Sole: a rounded slab, heel and toe domed.
    const soleRings: THREE.Vector3[][] = [];
    const domes = detail > 0 ? [0.0, 0.35, 0.7] : [0];
    const soleRing = (z: number, k: number) => {
      const hw = width(z) / 2 + 0.004 * s;
      const y0 = spring(z);
      const r: THREE.Vector3[] = [];
      const ns = detail === 2 ? 12 : detail === 1 ? 7 : 6;
      for (let j = 0; j < ns; j++) {
        const a = (j / ns) * Math.PI * 2;
        const cx = Math.cos(a);
        const sy = Math.sin(a);
        const yy = y0 + (sole * (1 + spow(sy, 0.25))) / 2;
        r.push(at(hw * spow(cx, 0.3) * k, lerp(y0 + sole / 2, yy, k), z));
      }
      return r;
    };
    for (const k of domes) soleRings.push(soleRing(z0 - 0.004 * s * (1 - k), k * 0.9 + 0.05));
    for (const z of zs) soleRings.push(soleRing(z, 1));
    for (const k of [...domes].reverse()) soleRings.push(soleRing(z1 + 0.004 * s * (1 - k), k * 0.9 + 0.05));
    loft(soup, soleRings, {
      w: (_i, _j, p) => wt(p),
      color: (_i, _j, p) => (p.y - spring(p.z) < sole * 0.3 ? L.sole.clone().multiplyScalar(0.55) : L.sole),
    });
    // Upper: a dome over the foot standing on the sole.
    const upRings: THREE.Vector3[][] = [];
    const upRing = (z: number, k: number) => {
      const hw = (width(z) / 2) * k;
      const h = height(z) * (0.25 + 0.75 * k);
      const y0 = spring(z) + sole - 0.002 * s;
      const r: THREE.Vector3[] = [];
      for (let j = 0; j < nA; j++) {
        const a = (j / nA) * Math.PI * 2;
        const sy = Math.sin(a);
        r.push(at(hw * spow(Math.cos(a), 0.55), y0 + (sy > 0 ? h * spow(sy, 0.8) : 0.002 * s * sy), z));
      }
      return r;
    };
    upRings.push(upRing(z0 + 0.004 * s, 0.35));
    upRings.push(upRing(z0 + 0.01 * s, 0.8));
    for (const z of zs) if (z > z0 + 0.014 * s && z < z1 - 0.012 * s) upRings.push(upRing(z, 1));
    upRings.push(upRing(z1 - 0.007 * s, 0.8));
    upRings.push(upRing(z1 - 0.002 * s, 0.35));
    loft(soup, upRings, { w: (_i, _j, p) => wt(p), color: L.upper });
    // The ankle in its sock, from inside the shoe to well up under the hem.
    const iS = BI[`shin_${side}`]!;
    const sockC = new THREE.Color("#1c1c1f");
    const sockRings: THREE.Vector3[][] = [];
    const nSock = detail === 2 ? 14 : detail === 1 ? 7 : 5;
    const sockYs = [sole + 0.02 * s, B.ankY - 0.02 * s, B.ankY + 0.01 * s, B.ankY + 0.05 * s, B.ankY + 0.085 * s];
    const sockR = [0.034, 0.036, 0.034, 0.033, 0.035];
    sockYs.forEach((y, k) => {
      const r: THREE.Vector3[] = [];
      for (let j = 0; j < nSock; j++) {
        const a = (j / nSock) * Math.PI * 2;
        r.push(V(ank.x + Math.sin(a) * sockR[k]! * 0.92 * s, y, ank.z - 0.012 * s + Math.cos(a) * sockR[k]! * s));
      }
      sockRings.push(r);
    });
    loft(soup, sockRings, {
      w: (_i, _j, p) => {
        const e = smooth(B.ankY - 0.015 * s, B.ankY + 0.03 * s, p.y);
        return [[iF, 1 - e], [iS, e]];
      },
      color: sockC,
    });
    if (detail === 2) {
      if (!clog) {
        // Laces across the instep, and the tongue's edge.
        for (let k = 0; k < 4; k++) {
          const z = lerp(0.01, 0.085, k / 3) * s;
          const y = spring(z) + sole + height(z) - 0.002 * s;
          addGeo(soup, new THREE.CylinderGeometry(0.0018 * s, 0.0018 * s, width(z) * 0.42, 6), mat4(at(0, y, z), new THREE.Quaternion().setFromEuler(new THREE.Euler(0.2, 0, Math.PI / 2))), wt(at(0, y, z)), L.sole);
        }
      } else {
        // A clog's welt seam: a thin bead where the upper meets the sole.
        const bead: THREE.Vector3[] = [];
        const bn: THREE.Vector3[] = [];
        for (let k = 0; k <= 40; k++) {
          const a = (k / 40) * Math.PI * 2;
          const z = lerp(z0 + 0.006, z1 - 0.006, 0.5 - 0.5 * Math.cos(a));
          const hw = width(z) / 2 + 0.001;
          bead.push(at(Math.sin(a) * hw, spring(z) + sole + 0.001, z));
          bn.push(V(Math.sin(a), 0.2, 0).normalize());
        }
        tape(soup, bead, bn, 0.003 * s, 0.002 * s, wt, L.upper.clone().multiplyScalar(0.7), detail);
      }
    }
  }
}


/* ------------------------------------------------------------------------------------------------------------ */
/* Forearms and hands                                                                                             */
/* ------------------------------------------------------------------------------------------------------------ */

function buildArmSkin(soup: Soup, B: Body, side: Side, L: Look, detail: Detail): void {
  const s = B.s;
  const H = B.hands[side];
  const sh = H.sh;
  const HL = H.HL;
  const el = B.J[`forearm_${side}`]!;
  const fo = B.fore;
  const rolled = L.sleeves === "rolled";
  const t0 = rolled ? 0.03 * s : fo - 0.055 * s;
  // Stations from the elbow along the forearm and on through the palm to the knuckles: [t, r_thumbside, r_palm, n].
  const prof: [number, number, number, number][] = [
    [0.1 * fo, 0.04 * s, 0.037 * s, 2.0],
    [0.5 * fo, 0.034 * s, 0.029 * s, 2.0],
    [fo - 0.035 * s, 0.028 * s, 0.02 * s, 2.1],
    [fo, 0.0272 * sh, 0.0182 * sh, 2.3],
    [fo + 0.14 * HL, 0.0335 * sh, 0.0168 * sh, 2.8],
    [fo + 0.34 * HL, 0.0385 * sh, 0.0152 * sh, 3.0],
    [fo + 0.49 * HL, 0.0405 * sh, 0.0135 * sh, 3.0],
    [fo + 0.535 * HL, 0.036 * sh, 0.0105 * sh, 2.6],
    [fo + 0.56 * HL, 0.022 * sh, 0.0065 * sh, 2.2],
    [fo + 0.568 * HL, 0.006 * sh, 0.002 * sh, 2.0],
  ];
  const pts = prof.filter((p) => p[0] >= t0 - 1e-6);
  pts.unshift([t0, ...(() => {
    const i = prof.findIndex((p) => p[0] > t0);
    const a = prof[Math.max(0, i - 1)]!;
    const b = prof[i]!;
    const u = clamp((t0 - a[0]) / Math.max(1e-6, b[0] - a[0]), 0, 1);
    return [lerp(a[1], b[1], u), lerp(a[2], b[2], u), lerp(a[3], b[3], u)] as [number, number, number];
  })()]);
  const ts: number[] = [];
  const steps = detail === 2 ? 2 : detail === 1 ? 1 : 1;
  for (let i = 0; i < pts.length - 1; i++) for (let k = 0; k < steps; k++) ts.push(lerp(pts[i]![0], pts[i + 1]![0], k / steps));
  ts.push(pts[pts.length - 1]![0]);
  const tt = pts.map((p) => p[0]);
  const fR = new Mono(tt, pts.map((p) => p[1]));
  const fN = new Mono(tt, pts.map((p) => p[2]));
  const fE = new Mono(tt, pts.map((p) => p[3]));
  const nA = detail === 2 ? 14 : detail === 1 ? 10 : 5;
  const rings: THREE.Vector3[][] = [];
  for (const t of ts) {
    const c = el.clone().addScaledVector(H.da, t);
    const rt = fR.at(t);
    const rn = fN.at(t);
    const e = 2 / fE.at(t);
    const palm = smooth(fo, fo + 0.1 * HL, t);
    c.addScaledVector(H.dn, 0.0025 * sh * palm);
    const ring: THREE.Vector3[] = [];
    for (let j = 0; j < nA; j++) {
      const a = (j / nA) * Math.PI * 2;
      const sn = Math.sin(a);
      // The back of the hand is flatter than the palm.
      const kn = sn < 0 ? 1 - 0.18 * palm : 1;
      ring.push(c.clone().addScaledVector(H.dt, spow(Math.cos(a), e) * rt).addScaledVector(H.dn, spow(sn, e) * rn * kn));
    }
    rings.push(ring);
  }
  const iFo = BI[`forearm_${side}`]!;
  const iHa = BI[`hand_${side}`]!;
  const fingerAt = (p: THREE.Vector3) => {
    const t = V().subVectors(p, H.wrist).dot(H.dt);
    let best: Finger = "middle";
    let bd = 1e9;
    for (const f of ["index", "middle", "ring", "pinky"] as const) {
      const d = Math.abs(V().subVectors(H.digits[f].base, H.wrist).dot(H.dt) - t);
      if (d < bd) {
        bd = d;
        best = f;
      }
    }
    return BI[`${best}_01_${side}`]!;
  };
  loft(soup, rings, {
    w: (i, _j, p) => {
      const t = ts[i]!;
      const e = smooth(fo - 0.018 * s, fo + 0.012 * s, t);
      const k = 0.5 * smooth(fo + 0.42 * HL, fo + 0.545 * HL, t);
      return [[iFo, 1 - e], [iHa, e * (1 - k)], [fingerAt(p), e * k]];
    },
    color: L.skin,
    out: (i, _j, p) => p.clone().sub(el.clone().addScaledVector(H.da, ts[i]!)),
  });
  for (const f of FINGERS) buildDigit(soup, H, f, L, detail);
}

function buildDigit(soup: Soup, H: HandDef, f: Finger, L: Look, detail: Detail): void {
  const g = H.digits[f];
  const side = H.side;
  const [l1, l2, l3] = g.len;
  const r0 = g.rad;
  const thumb = f === "thumb";
  const total = l1 + l2 + l3;
  const start = thumb ? -0.028 * H.sh : -0.013 * H.sh;
  // Radius along the digit: knuckles a little proud, the pulp of the tip full, the tip domed.
  const prof: [number, number][] = thumb
    ? [[start, 0.6], [0, 1.45], [0.4 * l1, 1.32], [l1, 1.0], [l1 + 0.5 * l2, 0.95], [l1 + l2, 0.92], [l1 + l2 + 0.5 * l3, 0.92], [total - 0.45 * r0, 0.72], [total - 0.12 * r0, 0.35], [total, 0.02]]
    : [[start, 1.02], [0.12 * l1, 1.02], [0.6 * l1, 0.93], [l1, 0.97], [l1 + 0.5 * l2, 0.86], [l1 + l2, 0.86], [l1 + l2 + 0.55 * l3, 0.84], [total - 0.5 * r0, 0.7], [total - 0.15 * r0, 0.36], [total, 0.02]];
  const fr = new Mono(prof.map((p) => p[0]), prof.map((p) => p[1]));
  const n = detail === 2 ? 13 : detail === 1 ? 6 : 3;
  const ss: number[] = [];
  for (let i = 0; i <= n; i++) ss.push(lerp(start, total, 1 - (1 - i / n) ** 1.25));
  const nA = detail === 2 ? 8 : detail === 1 ? 6 : 3;
  const across = g.k;
  const palmward = V().crossVectors(g.k, g.dir).normalize();
  if (palmward.dot(g.back) > 0) palmward.negate();
  const rings: THREE.Vector3[][] = [];
  for (const sv of ss) {
    const c = g.base.clone().addScaledVector(g.dir, sv);
    const r = fr.at(sv) * r0;
    const pulp = smooth(l1 + l2, l1 + l2 + 0.4 * l3, sv) * (1 - smooth(total - 0.6 * r0, total, sv));
    c.addScaledVector(palmward, 0.12 * r0 * pulp);
    const ring: THREE.Vector3[] = [];
    for (let j = 0; j < nA; j++) {
      const a = (j / nA) * Math.PI * 2;
      ring.push(c.clone().addScaledVector(across, Math.cos(a) * r * 1.07).addScaledVector(palmward, Math.sin(a) * r * 0.92));
    }
    rings.push(ring);
  }
  const b1 = BI[`${f}_01_${side}`]!;
  const b2 = BI[`${f}_02_${side}`]!;
  const b3 = BI[`${f}_03_${side}`]!;
  const iHa = BI[`hand_${side}`]!;
  const wAt = (sv: number): Wts => {
    const h = thumb ? smooth(-0.012 * H.sh, 0.012 * H.sh, sv) : smooth(-0.004 * H.sh, 0.006 * H.sh, sv);
    const e2 = smooth(l1 - 0.0038, l1 + 0.0038, sv);
    const e3 = smooth(l1 + l2 - 0.003, l1 + l2 + 0.003, sv);
    return [[iHa, 1 - h], [b1, h * (1 - e2)], [b2, h * e2 * (1 - e3)], [b3, h * e3]];
  };
  loft(soup, rings, { w: (i) => wAt(ss[i]!), color: L.skin, out: (i, _j, p) => p.clone().sub(g.base.clone().addScaledVector(g.dir, ss[i]!)) });
  if (detail === 2) {
    // The nail: a thin curved plate on the back of the distal phalanx.
    const s0 = l1 + l2 + 0.38 * l3;
    const s1 = total - 0.2 * r0;
    const rows: THREE.Vector3[][] = [];
    for (let k = 0; k <= 4; k++) {
      const sv = lerp(s0, s1, k / 4);
      const c = g.base.clone().addScaledVector(g.dir, sv);
      const pulp = smooth(l1 + l2, l1 + l2 + 0.4 * l3, sv) * (1 - smooth(total - 0.6 * r0, total, sv));
      c.addScaledVector(palmward, 0.12 * r0 * pulp);
      const r = fr.at(sv) * r0 + 0.0004;
      const row: THREE.Vector3[] = [];
      const half = 0.78 * (k === 0 ? 0.85 : 1);
      for (let j = 0; j <= 4; j++) {
        const a = -Math.PI / 2 + lerp(-half, half, j / 4);
        row.push(c.clone().addScaledVector(across, Math.cos(a) * r * 1.07).addScaledVector(palmward, Math.sin(a) * r * 0.92));
      }
      rows.push(row);
    }
    loft(soup, rows, { closed: false, w: () => [[b3, 1]], color: L.nail, out: () => palmward.clone().negate() });
  }
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Neck, head (sculpted, faceless), ears, hair and caps                                                           */
/* ------------------------------------------------------------------------------------------------------------ */

interface HeadShape {
  c: THREE.Vector3;
  /** Distance from the centre to the surface along unit direction d. */
  r(d: THREE.Vector3): number;
  /** The cap's lower edge: the plane d . n = k (d unit). */
  capN: THREE.Vector3;
  capK: number;
}

function headShape(B: Body): HeadShape {
  const s = B.s;
  const fem = B.fem;
  const c = V(0, B.topY - 0.098 * s, 0.008 * s);
  const ax = (0.077 - 0.004 * fem) * s;
  const ayU = 0.098 * s;
  const ayD = (0.134 - 0.007 * fem) * s;
  const azF = 0.096 * s;
  const azB = (0.103 - 0.003 * fem) * s;
  const bump = (d: THREE.Vector3, x: number, y: number, z: number, wx: number, wy: number, a: number) => {
    if (d.z * z < -0.2 && z !== 0) return 0;
    return a * Math.exp(-(((d.x - x) / wx) ** 2) - ((d.y - y) / wy) ** 2);
  };
  const r = (d: THREE.Vector3) => {
    const low = smooth(0.05, -0.85, d.y);
    const fwd = smooth(-0.2, 0.6, d.z);
    const axE = ax * (1 - 0.15 * low * (0.5 + 0.5 * fwd));
    const azE = d.z >= 0 ? azF * (1 - 0.05 * low) : azB * (1 - 0.42 * smooth(-0.2, -0.85, d.y) * smooth(0.05, -0.55, d.z));
    const ay = d.y >= 0 ? ayU : ayD * (1 - 0.2 * smooth(0.1, -0.7, d.z));
    // Below the cheekbones the front view squares off (jaw), then rounds at the chin.
    const pe = 2 + 0.95 * low * (0.4 + 0.6 * fwd);
    const xy = Math.pow(Math.pow(Math.abs(d.x / axE), pe) + Math.pow(Math.abs(d.y / ay), pe), 2 / pe);
    let rr = 1 / Math.sqrt(xy + (d.z / azE) ** 2);
    const sc = s;
    rr += bump(d, 0, 0.24, 1, 0.55, 0.11, 0.0028 * sc);
    rr -= bump(d, 0.33, 0.07, 1, 0.17, 0.11, 0.0024 * sc) + bump(d, -0.33, 0.07, 1, 0.17, 0.11, 0.0024 * sc);
    rr += bump(d, 0.6, -0.14, 1, 0.2, 0.15, 0.0042 * sc) + bump(d, -0.6, -0.14, 1, 0.2, 0.15, 0.0042 * sc);
    rr += bump(d, 0, -0.7, 1, 0.28, 0.16, 0.0055 * sc * (1 - 0.4 * fem));
    rr += bump(d, 0.62, -0.6, 0, 0.18, 0.18, 0.0036 * sc * (1 - 0.5 * fem)) + bump(d, -0.62, -0.6, 0, 0.18, 0.18, 0.0036 * sc * (1 - 0.5 * fem));
    rr += d.z < -0.5 ? 0.0045 * sc * Math.exp(-(d.x ** 2) / 0.2 - ((d.y - 0.05) ** 2) / 0.12) : 0;
    rr -= bump(d, 0.95, 0.25, 0, 0.2, 0.2, 0.0028 * sc) + bump(d, -0.95, 0.25, 0, 0.2, 0.2, 0.0028 * sc);
    return rr;
  };
  // The cap sits higher at the forehead than at the nape: a plane through mid-forehead and the occiput.
  const capN = V(0, 1, -0.24).normalize();
  const capK = 0.3 * capN.y + 0.954 * capN.z;
  return { c, r, capN, capK };
}

function dirOf(lat: number, lon: number, out = V()): THREE.Vector3 {
  // lon 0 faces +z (front), increasing towards +x (the left); lat from -pi/2 (down) to pi/2 (up).
  return out.set(Math.cos(lat) * Math.sin(lon), Math.sin(lat), Math.cos(lat) * Math.cos(lon));
}

function buildNeck(soup: Soup, B: Body, L: Look, detail: Detail): void {
  const s = B.s;
  const J = B.J;
  const r = (0.052 + 0.007 * (1 - B.fem)) * s;
  const bot = V(0, B.neckBaseY - 0.035 * s, -0.026 * s);
  const top = V(0, J.head!.y + 0.04 * s, J.head!.z + 0.022 * s);
  const n = detail === 2 ? 6 : 5;
  const nA = detail === 2 ? 14 : detail === 1 ? 12 : 7;
  const rings: THREE.Vector3[][] = [];
  for (let i = 0; i <= n; i++) {
    const u = i / n;
    const c = bot.clone().lerp(top, u);
    const k = 1 + 0.12 * smooth(0.6, 1, u) + 0.06 * smooth(0.25, 0, u);
    const ring: THREE.Vector3[] = [];
    for (let j = 0; j < nA; j++) {
      const a = (j / nA) * Math.PI * 2;
      const cz = Math.cos(a);
      ring.push(c.clone().add(V(Math.sin(a) * r * 1.06 * k, 0, cz * r * (cz > 0 ? 0.9 : 1.02) * k)));
    }
    rings.push(ring);
  }
  loft(soup, rings, {
    w: (_i, _j, p) => {
      const a = smooth(B.neckBaseY - 0.01 * s, B.neckBaseY + 0.02 * s, p.y);
      const b = smooth(J.head!.y - 0.015 * s, J.head!.y + 0.03 * s, p.y);
      return [[BI.spine_03!, 1 - a], [BI.neck!, a * (1 - b)], [BI.head!, a * b]];
    },
    color: L.skin,
  });
}

function buildHead(skin: Soup, whites: Soup, hair: Soup, B: Body, L: Look, detail: Detail): number {
  const s = B.s;
  const Hs = headShape(B);
  const head = rigid("head");
  const nLat = detail === 2 ? 22 : detail === 1 ? 14 : 8;
  const nLon = detail === 2 ? 30 : detail === 1 ? 18 : 8;
  const d = V();
  const surf = (lat: number, lon: number, off = 0) => {
    dirOf(lat, lon, d);
    return Hs.c.clone().addScaledVector(d, Hs.r(d) + off);
  };
  const rings: THREE.Vector3[][] = [];
  for (let i = 0; i <= nLat; i++) {
    const lat = -Math.PI / 2 + (Math.PI * i) / nLat;
    const ring: THREE.Vector3[] = [];
    for (let j = 0; j < nLon; j++) ring.push(surf(lat, (j / nLon) * Math.PI * 2));
    rings.push(ring);
  }
  loft(skin, rings, { w: () => head, color: L.skin, out: (_i, _j, p) => p.clone().sub(Hs.c) });

  // Ears: a rolled helix round a hollowed shell, set out from the side of the head and tilted back.
  for (const sd of [1, -1]) {
    const ec = Hs.c.clone().add(V(sd * (0.074 - 0.004 * B.fem) * s, -0.03 * s, -0.018 * s));
    const out = V(sd, 0, -0.32).normalize();
    const upA = V(0, 1, 0.22).normalize();
    const back = V().crossVectors(upA, out).multiplyScalar(sd);
    const eh = 0.031 * s;
    const ew = 0.018 * s;
    const nK = detail === 2 ? 5 : 3;
    const nE = detail === 2 ? 14 : 10;
    const rimH = [0.0016, 0.0058, 0.0044, 0.0034, 0.0052, 0.0024, 0.0006, 0.0];
    const front: THREE.Vector3[][] = [];
    const backS: THREE.Vector3[][] = [];
    for (let k = 0; k <= nK; k++) {
      const u = k / nK;
      const scale = 1 - u;
      const hh = rimH[Math.round(u * (rimH.length - 1))]! * s;
      const fr: THREE.Vector3[] = [];
      const br: THREE.Vector3[] = [];
      for (let j = 0; j < nE; j++) {
        const a = (j / nE) * Math.PI * 2;
        // Egg outline: broad at the top, a soft lobe below; the front edge straighter where it meets the head.
        const yy = Math.sin(a) * eh * (Math.sin(a) > 0 ? 1 : 0.92);
        const xx = Math.cos(a) * ew * (Math.sin(a) > 0 ? 1.08 : 0.8) * (Math.cos(a) > 0 ? 1 : 0.72);
        const base = ec.clone().addScaledVector(upA, yy * scale).addScaledVector(back, xx * scale);
        const lobe = Math.sin(a) < -0.6 ? 0.0012 * s : 0;
        fr.push(base.clone().addScaledVector(out, hh + 0.002 * s + lobe));
        br.push(base.clone().addScaledVector(out, -0.004 * s * u + 0.0005 * s));
      }
      front.push(fr);
      backS.push(br);
    }
    loft(skin, front, { w: () => head, color: L.skin, out: () => out });
    loft(skin, [front[0]!, ...backS], { w: () => head, color: L.skinDark, out: () => out.clone().negate() });
  }

  // The cap line, as a latitude for each longitude.
  const capLat = (lon: number) => {
    let lo = -Math.PI / 2;
    let hi = Math.PI / 2;
    for (let k = 0; k < 24; k++) {
      const m = (lo + hi) / 2;
      dirOf(m, lon, d);
      if (d.dot(Hs.capN) > Hs.capK) hi = m;
      else lo = m;
    }
    return (lo + hi) / 2;
  };
  // Hair: from the hairline (under the cap at the front, round above the ears, down to the nape) up under the cap.
  const hairLow = (lon: number) => {
    const psi = Math.abs(wrapPi(lon));
    const k: [number, number][] = [[0.85, 0.02], [1.12, -0.36], [1.32, -0.3], [1.48, 0.06], [1.86, 0.06], [2.1, -0.44], [2.6, -0.6], [Math.PI, -0.64]];
    if (psi <= k[0]![0]) return Math.asin(k[0]![1]);
    for (let i = 1; i < k.length; i++) if (psi <= k[i]![0]) return Math.asin(lerp(k[i - 1]![1], k[i]![1], smooth(k[i - 1]![0], k[i]![0], psi)));
    return Math.asin(-0.64);
  };
  const nHu = detail === 2 ? 36 : detail === 1 ? 24 : 12;
  const nHv = detail === 2 ? 6 : 4;
  const lon0 = 0.8;
  const hRows: THREE.Vector3[][] = [];
  for (let i = 0; i <= nHv; i++) {
    const row: THREE.Vector3[] = [];
    for (let j = 0; j <= nHu; j++) {
      const lon = lon0 + ((Math.PI * 2 - 2 * lon0) * j) / nHu;
      const a = hairLow(lon);
      const b = capLat(lon) + 0.12;
      const lat = lerp(a, Math.max(a, b), i / nHv);
      const th = (0.0016 + 0.0042 * smooth(0, 0.35, i / nHv)) * s * (L.bun ? 0.8 : 1);
      row.push(surf(lat, lon, th));
    }
    hRows.push(row);
  }
  loft(hair, hRows, { closed: false, w: () => head, color: L.hair, out: (_i, _j, p) => p.clone().sub(Hs.c), uv: (_i, _j, _p, u, v) => [u, -v] });
  // A bun at the nape for longer hair.
  if (L.bun && detail > 0) {
    const bc = Hs.c.clone().add(V(0, -0.045 * s, -0.1 * s));
    addGeo(hair, new THREE.SphereGeometry(0.036 * s, detail === 2 ? 18 : 9, detail === 2 ? 12 : 6), mat4(bc, new THREE.Quaternion().setFromEuler(new THREE.Euler(0.5, 0, 0)), V(1.15, 0.85, 0.9)), head, L.hair, 6);
  }

  // The cap.
  let capTris = whites.tris;
  if (L.capKind === "skull") {
    const nV = detail === 2 ? 8 : detail === 1 ? 6 : 3;
    const nU = detail === 2 ? 32 : detail === 1 ? 22 : 10;
    const rows: THREE.Vector3[][] = [];
    for (let i = 0; i <= nV; i++) {
      const row: THREE.Vector3[] = [];
      for (let j = 0; j < nU; j++) {
        const lon = (j / nU) * Math.PI * 2;
        const lat = lerp(capLat(lon) - 0.02, Math.PI / 2, (i / nV) ** 1.1);
        row.push(surf(lat, lon, 0.0072 * s + (i === 0 ? 0.0 : 0.0014 * s)));
      }
      rows.push(row);
    }
    loft(whites, rows, { w: () => head, color: L.cap, out: (_i, _j, p) => p.clone().sub(Hs.c) });
    // The turned-up band.
    const band: THREE.Vector3[][] = [];
    for (const [dl, off] of [[-0.03, 0.006], [-0.034, 0.0105], [0.08, 0.0108], [0.095, 0.0082]] as [number, number][]) {
      const row: THREE.Vector3[] = [];
      for (let j = 0; j < nU; j++) {
        const lon = (j / nU) * Math.PI * 2;
        row.push(surf(capLat(lon) + dl, lon, off * s));
      }
      band.push(row);
    }
    loft(whites, band, { w: () => head, color: L.cap, out: (_i, _j, p) => p.clone().sub(Hs.c) });
  } else {
    // A pleated paper toque: a headband, a column of pleats flaring out, a puffed crown.
    const nU = detail === 2 ? 48 : detail === 1 ? 32 : 10;
    const pleats = detail === 0 ? 0 : 24;
    const bandRows: THREE.Vector3[][] = [];
    const topBand: THREE.Vector3[] = [];
    for (const [dl, off] of [[-0.02, 0.007], [0.12, 0.007]] as [number, number][]) {
      const row: THREE.Vector3[] = [];
      for (let j = 0; j < nU; j++) {
        const lon = (j / nU) * Math.PI * 2;
        row.push(surf(capLat(lon) + dl, lon, off * s));
      }
      bandRows.push(row);
    }
    for (const p of bandRows[1]!) topBand.push(p.clone());
    loft(whites, bandRows, { w: () => head, color: L.cap, out: (_i, _j, p) => p.clone().sub(Hs.c) });
    const cy = topBand.reduce((a, p) => a + p.y, 0) / nU;
    const ctr = topBand.reduce((a, p) => a.add(p), V()).multiplyScalar(1 / nU);
    const hT = 0.2 * s;
    const col: THREE.Vector3[][] = [];
    const nV = detail === 2 ? 10 : 4;
    for (let i = 0; i <= nV; i++) {
      const u = i / nV;
      const row: THREE.Vector3[] = [];
      for (let j = 0; j < nU; j++) {
        const lon = (j / nU) * Math.PI * 2;
        const p0 = topBand[j]!;
        const rad = V(p0.x - ctr.x, 0, p0.z - ctr.z);
        const r0 = rad.length();
        const flare = 1 + 0.2 * Math.sin(u * Math.PI * 0.6) + 0.05 * u;
        const pleat = pleats ? 0.006 * s * Math.abs(Math.sin((pleats * lon) / 2)) * smooth(0, 0.2, u) : 0;
        rad.setLength(r0 * flare + pleat);
        row.push(V(ctr.x + rad.x, lerp(p0.y, cy + hT, u), ctr.z + rad.z));
      }
      col.push(row);
    }
    // The crown puffs over the top of the column.
    const last = col[col.length - 1]!;
    for (const [k, dy] of [[1.03, 0.012], [0.85, 0.028], [0.5, 0.036], [0.0, 0.038]] as [number, number][]) {
      col.push(last.map((p) => V(ctr.x + (p.x - ctr.x) * k, cy + hT + dy * s, ctr.z + (p.z - ctr.z) * k)));
    }
    loft(whites, col, { w: () => head, color: L.cap, out: (_i, _j, p) => V(p.x - ctr.x, p.y > cy + hT ? 1 : 0, p.z - ctr.z) });
  }
  capTris = whites.tris - capTris;
  return capTris;
}


/* ------------------------------------------------------------------------------------------------------------ */
/* Poses                                                                                                         */
/* ------------------------------------------------------------------------------------------------------------ */

/**
 * How a hand holds: a tool's handle in the oblique grip with the working end past the thumb (skillet, knife,
 * spoon), the same grip with the working end below the little finger (`dagger`: a ladle or paddle stirring a deep
 * pot), a bar or pan rim (power grip, hook), or flat on a surface.
 */
export type GripKind = "tool" | "dagger" | "bar" | "hook" | "palm";
/**
 * Finger shapes. `close` blends from relaxed to the shape. `grasp` closes every digit, joint by joint, until it
 * touches the held item's surface (`Staff.itemSdf`, in the item socket's frame): fingers flat on a wall, a thumb
 * over a lip.
 */
export type HandShape = "relax" | "tool" | "bar" | "hook" | "claw" | "flat" | "grasp";

export interface FootPose {
  /** Ankle joint, body space. */
  pos: THREE.Vector3;
  /** Foot orientation, body space (identity: flat, pointing +z). */
  quat: THREE.Quaternion;
  /** Toe (MTP) extension relative to the foot, radians. */
  toe: number;
}
export interface ArmPose {
  /** Free arm (used where `ik` < 1): shoulder flexion, abduction, elbow flexion, humeral twist, radians. */
  flex: number;
  abd: number;
  elbow: number;
  twist: number;
  /** Weight of the IK target, 0..1. */
  ik: number;
  grip: GripKind;
  /** The grip point to reach, body space. */
  pos: THREE.Vector3;
  /** Target frame: tool: handle axis (+z of the grip frame) and its up (+y); bar/hook: the bar's axis towards the
   *  thumb side and the palm's facing; palm: the fingers' direction and the palm's facing. */
  a: THREE.Vector3;
  b: THREE.Vector3;
  /** Where the elbow should point, body space. */
  pole: THREE.Vector3;
  shape: HandShape;
  close: number;
  /** Handle radius for the tool and bar shapes, metres. */
  radius: number;
  /** How far the hand may roll round the handle axis to stay comfortable (round handles), radians. */
  roll: number;
  /**
   * A hook on a lip may sit anywhere from fingers-forward under it to fingers-down over its edge: the target
   * (frame and point) may turn by 0..`spinMax` about `spinAxis` through `spinPivot`, shifting by `spinShift` times
   * the sine of the turn, whichever suits the arm. `spinMax` 0 for none.
   */
  spinMax: number;
  spinAxis: THREE.Vector3;
  spinPivot: THREE.Vector3;
  spinShift: THREE.Vector3;
}
export interface BodyPose {
  /** Pelvis (between the hip joints), body space. */
  pelvis: THREE.Vector3;
  /** Pelvis forward tilt, roll (left hip up positive), yaw, radians. */
  tilt: number;
  roll: number;
  yaw: number;
  /** Trunk on the pelvis: flexion (forward positive), side bend (to the right positive), twist (to the left positive). */
  flex: number;
  side: number;
  twist: number;
  /** Clavicles: elevation and protraction, left then right. */
  clav: [number, number, number, number];
  feet: [FootPose, FootPose];
  arms: [ArmPose, ArmPose];
  /** Where the eyes would look, body space. */
  gaze: THREE.Vector3;
  /** The two-handed item's frame, body space (when `itemOn`). */
  itemOn: boolean;
  itemPos: THREE.Vector3;
  itemQuat: THREE.Quaternion;
}

function armPose(): ArmPose {
  return { flex: 0, abd: 0, elbow: 0.2, twist: 0, ik: 0, grip: "tool", pos: V(), a: Z.clone(), b: Y.clone(), pole: V(0, -1, -0.5).normalize(), shape: "relax", close: 0, radius: 0.013, roll: 0, spinMax: 0, spinAxis: X.clone(), spinPivot: V(), spinShift: V() };
}
function footPose(): FootPose {
  return { pos: V(), quat: new THREE.Quaternion(), toe: 0 };
}
export function newPose(): BodyPose {
  return {
    pelvis: V(),
    tilt: 0,
    roll: 0,
    yaw: 0,
    flex: 0,
    side: 0,
    twist: 0,
    clav: [0, 0, 0, 0],
    feet: [footPose(), footPose()],
    arms: [armPose(), armPose()],
    gaze: V(0, 1.6, 3),
    itemOn: false,
    itemPos: V(),
    itemQuat: new THREE.Quaternion(),
  };
}
function copyArm(o: ArmPose, a: ArmPose): void {
  o.flex = a.flex;
  o.abd = a.abd;
  o.elbow = a.elbow;
  o.twist = a.twist;
  o.ik = a.ik;
  o.grip = a.grip;
  o.pos.copy(a.pos);
  o.a.copy(a.a);
  o.b.copy(a.b);
  o.pole.copy(a.pole);
  o.shape = a.shape;
  o.close = a.close;
  o.radius = a.radius;
  o.roll = a.roll;
  o.spinMax = a.spinMax;
  o.spinAxis.copy(a.spinAxis);
  o.spinPivot.copy(a.spinPivot);
  o.spinShift.copy(a.spinShift);
}
export function copyPose(o: BodyPose, p: BodyPose): BodyPose {
  o.pelvis.copy(p.pelvis);
  o.tilt = p.tilt;
  o.roll = p.roll;
  o.yaw = p.yaw;
  o.flex = p.flex;
  o.side = p.side;
  o.twist = p.twist;
  o.clav = [...p.clav];
  for (let i = 0; i < 2; i++) {
    o.feet[i]!.pos.copy(p.feet[i]!.pos);
    o.feet[i]!.quat.copy(p.feet[i]!.quat);
    o.feet[i]!.toe = p.feet[i]!.toe;
    copyArm(o.arms[i]!, p.arms[i]!);
  }
  o.gaze.copy(p.gaze);
  o.itemOn = p.itemOn;
  o.itemPos.copy(p.itemPos);
  o.itemQuat.copy(p.itemQuat);
  return o;
}
/**
 * Blend two poses into `o` (which may be `a`). Feet are expected to agree (planted); arms blend their free swing
 * and their targets, and where only one side reaches, the reach weight ramps instead.
 */
export function blendPose(o: BodyPose, a: BodyPose, b: BodyPose, w: number): BodyPose {
  const k = clamp(w, 0, 1);
  o.pelvis.lerpVectors(a.pelvis, b.pelvis, k);
  o.tilt = lerp(a.tilt, b.tilt, k);
  o.roll = lerp(a.roll, b.roll, k);
  o.yaw = lerp(a.yaw, b.yaw, k);
  o.flex = lerp(a.flex, b.flex, k);
  o.side = lerp(a.side, b.side, k);
  o.twist = lerp(a.twist, b.twist, k);
  o.clav = [0, 1, 2, 3].map((i) => lerp(a.clav[i]!, b.clav[i]!, k)) as [number, number, number, number];
  for (let i = 0; i < 2; i++) {
    const fa = a.feet[i]!;
    const fb = b.feet[i]!;
    const fo = o.feet[i]!;
    fo.pos.lerpVectors(fa.pos, fb.pos, k);
    fo.quat.slerpQuaternions(fa.quat, fb.quat, k);
    fo.toe = lerp(fa.toe, fb.toe, k);
    const aa = a.arms[i]!;
    const ab = b.arms[i]!;
    const ao = o.arms[i]!;
    const ikA = aa.ik;
    const ikB = ab.ik;
    ao.flex = lerp(aa.flex, ab.flex, k);
    ao.abd = lerp(aa.abd, ab.abd, k);
    ao.elbow = lerp(aa.elbow, ab.elbow, k);
    ao.twist = lerp(aa.twist, ab.twist, k);
    if (ikA > 0 && ikB > 0) {
      ao.pos.lerpVectors(aa.pos, ab.pos, k);
      ao.a.lerpVectors(aa.a, ab.a, k).normalize();
      ao.b.lerpVectors(aa.b, ab.b, k).normalize();
      ao.pole.lerpVectors(aa.pole, ab.pole, k).normalize();
      ao.grip = k < 0.5 ? aa.grip : ab.grip;
    } else {
      const src = ikB > 0 ? ab : aa;
      ao.pos.copy(src.pos);
      ao.a.copy(src.a);
      ao.b.copy(src.b);
      ao.pole.copy(src.pole);
      ao.grip = src.grip;
    }
    ao.ik = lerp(ikA, ikB, k);
    const sameShape = aa.shape === ab.shape;
    if (sameShape) {
      ao.shape = aa.shape;
      ao.close = lerp(aa.close, ab.close, k);
    } else {
      // Open the first shape, then close the second.
      ao.shape = k < 0.5 ? aa.shape : ab.shape;
      ao.close = k < 0.5 ? aa.close * (1 - 2 * k) : ab.close * (2 * k - 1);
    }
    ao.radius = lerp(aa.radius, ab.radius, k);
    ao.roll = Math.max(aa.roll, ab.roll);
    const sa = aa.spinMax >= ab.spinMax ? aa : ab;
    const both = aa.spinMax > 0 && ab.spinMax > 0;
    const piv = both ? V().lerpVectors(aa.spinPivot, ab.spinPivot, k) : sa.spinPivot.clone();
    const ax = sa.spinAxis.clone();
    const sh = both ? V().lerpVectors(aa.spinShift, ab.spinShift, k) : sa.spinShift.clone();
    ao.spinMax = sa.spinMax;
    ao.spinAxis.copy(ax);
    ao.spinPivot.copy(piv);
    ao.spinShift.copy(sh);
  }
  o.gaze.lerpVectors(a.gaze, b.gaze, k);
  const aOn = a.itemOn;
  const bOn = b.itemOn;
  if (aOn && bOn) {
    o.itemPos.lerpVectors(a.itemPos, b.itemPos, k);
    o.itemQuat.slerpQuaternions(a.itemQuat, b.itemQuat, k);
  } else if (bOn || aOn) {
    o.itemPos.copy(bOn ? b.itemPos : a.itemPos);
    o.itemQuat.copy(bOn ? b.itemQuat : a.itemQuat);
  }
  o.itemOn = k < 0.5 ? aOn : bOn;
  return o;
}

/* ------------------------------------------------------------------------------------------------------------ */
/* A person                                                                                                       */
/* ------------------------------------------------------------------------------------------------------------ */

export interface StaffOptions {
  /** Apron colour (default: chosen from the seed). */
  apron?: THREE.ColorRepresentation;
  cap?: "skull" | "toque";
  seed?: number;
  trousers?: "check" | "black";
  sleeves?: "rolled" | "cuffed";
  shoes?: "clog" | "shoe";
  /** Cast shadows (default true). */
  shadows?: boolean;
  /**
   * "detailed" (the default): the cohort's MakeHuman-built body from public/staff/, swapped in once loaded (see
   * `ready`), with the procedural body drawn until then and kept if the file cannot be had. "procedural": that alone.
   */
  body?: "detailed" | "procedural";
}

const SKIN = ["#f2cdb0", "#e8b896", "#d29a72", "#b97d56", "#9a6141", "#7a4a30", "#5c3622", "#f0c6a0"];
const HAIR = ["#1a1411", "#2b1c13", "#3f2a1c", "#5d4029", "#8a6440", "#b89766", "#7b7874", "#2e2520"];
const APRONS = ["#26344d", "#e8e5dc", "#3b3d43", "#4a5a6e", "#5a4632", "#45563f", "#7a2f2a"];

function makeLook(cohort: Cohort, o: StaffOptions): Look {
  const r = rng((o.seed ?? 1) * 7919 + (cohort === "p05" ? 11 : cohort === "p50" ? 23 : 37));
  const pick = <T>(a: readonly T[]) => a[Math.floor(r() * a.length) % a.length]!;
  const skin = new THREE.Color(pick(SKIN));
  const hair = new THREE.Color(pick(HAIR));
  const white = new THREE.Color("#f3f1ea");
  const apron = new THREE.Color(o.apron ?? pick(APRONS));
  const tie = apron.clone().multiplyScalar(0.82);
  const check = o.trousers ? o.trousers === "check" : r() < 0.6;
  const shoes = o.shoes ?? (r() < 0.55 ? "clog" : "shoe");
  const upper = new THREE.Color(r() < 0.82 ? "#151517" : "#e9e7e2");
  const capKind = o.cap ?? (r() < 0.72 ? "skull" : "toque");
  return {
    skin,
    skinDark: skin.clone().multiplyScalar(0.8),
    nail: skin.clone().lerp(new THREE.Color("#f4d9d4"), 0.55),
    hair,
    white,
    button: r() < 0.75 ? white.clone().multiplyScalar(0.97) : new THREE.Color("#1d1d1f"),
    apron,
    apronTie: tie,
    trousers: check ? new THREE.Color("#ffffff") : new THREE.Color("#1c1d21"),
    upper,
    sole: new THREE.Color(shoes === "clog" ? "#232325" : "#3a3b3e"),
    cap: capKind === "skull" && r() < 0.4 ? new THREE.Color("#1e1f22") : white.clone(),
    sleeves: o.sleeves ?? (r() < 0.6 ? "rolled" : "cuffed"),
    shoes,
    capKind,
    check,
    bun: r() < (cohort === "p05" ? 0.7 : cohort === "p50" ? 0.3 : 0.05),
  };
}

const RELAX: Record<Finger, [number, number, number]> = { thumb: [0.12, 0.14, 0.12], index: [0.22, 0.32, 0.16], middle: [0.28, 0.4, 0.2], ring: [0.33, 0.45, 0.22], pinky: [0.4, 0.5, 0.25] };
const CLAW: Record<Finger, [number, number, number]> = { thumb: [0.55, 0.5, 0.4], index: [0.42, 1.5, 0.8], middle: [0.36, 1.55, 0.85], ring: [0.4, 1.5, 0.8], pinky: [0.5, 1.45, 0.75] };
const FLAT: Record<Finger, [number, number, number]> = { thumb: [-0.05, 0.05, 0.05], index: [0.04, 0.08, 0.04], middle: [0.05, 0.1, 0.05], ring: [0.06, 0.1, 0.05], pinky: [0.08, 0.12, 0.06] };

export interface WorkOpts {
  /** World point the work centres on: the pan's centre, the pot's rim centre, the board or plate centre, the
   *  rack bar's centre, the hotel pan's rim centre on the shelf. Defaults are set from the body and `surfaceY`. */
  at?: THREE.Vector3;
  /** World height of the work surface (grate, counter, board, pass shelf, rack bar, low shelf). */
  surfaceY?: number;
  /** Horizontal distance from the grip to the pan centre (saute), metres. */
  reach?: number;
  /** Half the distance between a carried item's grip points along its x axis (hotel pan: 0.27). */
  halfWidth?: number;
  /** Tempo multiplier. */
  tempo?: number;
}

const _v1 = V();
const _v2 = V();
const _v3 = V();
const _v4 = V();
const _q1 = new THREE.Quaternion();
const _q2 = new THREE.Quaternion();
const _q3 = new THREE.Quaternion();

/* ------------------------------------------------------------------------------------------------------------ */
/* Detailed bodies. Each cohort has a .glb (public/staff/<cohort>.glb) built in Blender from MakeHuman's CC0 assets */
/* with MPFB 2: a synthetic body shaped to the cohort's row by parametric targets (no scan, no real person), skin    */
/* and eye textures, brows, lashes, short hair under a skull cap, a double-breasted jacket, trousers, a bib apron    */
/* and shoes, rigged to this file's skeleton (the same bone names at the same joints, identity binds). A person is   */
/* drawn with the procedural body until the file is in, then swaps to it; the motion layer drives both alike.       */
/* ------------------------------------------------------------------------------------------------------------ */

let staffBase: string | null = null;

/** Where the people's .glb files are served from (default: the site's base URL + "staff/"). */
export function setStaffAssetBase(url: string): void {
  staffBase = url.endsWith("/") ? url : `${url}/`;
}

function staffUrl(cohort: Cohort): string {
  const base = staffBase ?? `${import.meta.env?.BASE_URL ?? "/"}staff/`;
  return `${base}${cohort}.glb`;
}

/** One mesh of a detailed body: geometry in body space at the bind pose, its bones indexed as `ALL_BONES`. */
interface DetailedPart {
  role: string;
  geometry: THREE.BufferGeometry;
  source: THREE.MeshStandardMaterial;
  /** Metres of surface per unit of texture coordinate (for patterns drawn to a real size). */
  mpu: number;
}
interface DetailedAsset {
  parts: DetailedPart[];
}

const assetCache = new Map<Cohort, Promise<DetailedAsset>>();

/** Load a cohort's detailed body once; every person of that cohort shares its geometry and textures. */
export function loadDetailedBody(cohort: Cohort): Promise<DetailedAsset> {
  let p = assetCache.get(cohort);
  if (!p) {
    const loader = new GLTFLoader();
    loader.setMeshoptDecoder(MeshoptDecoder);
    p = loader.loadAsync(staffUrl(cohort)).then((g) => prepareDetailed(cohort, g));
    assetCache.set(cohort, p);
  }
  return p;
}

const _dm = new THREE.Matrix4();
const _dn = new THREE.Matrix3();
const _dp = new THREE.Vector3();
const _dq = new THREE.Vector3();

function metresPerUv(g: THREE.BufferGeometry): number {
  const pos = g.getAttribute("position");
  const uv = g.getAttribute("uv");
  if (!uv) return 1;
  const ix = g.index;
  const n = ix ? ix.count : pos.count;
  const at = (k: number) => (ix ? ix.getX(k) : k);
  let a3 = 0;
  let a2 = 0;
  const p = [V(), V(), V()];
  const u = [new THREE.Vector2(), new THREE.Vector2(), new THREE.Vector2()];
  for (let t = 0; t + 2 < n; t += 3) {
    for (let k = 0; k < 3; k++) {
      p[k]!.fromBufferAttribute(pos as THREE.BufferAttribute, at(t + k));
      u[k]!.fromBufferAttribute(uv as THREE.BufferAttribute, at(t + k));
    }
    a3 += _dp.subVectors(p[1]!, p[0]!).cross(_dq.subVectors(p[2]!, p[0]!)).length() / 2;
    a2 += Math.abs((u[1]!.x - u[0]!.x) * (u[2]!.y - u[0]!.y) - (u[2]!.x - u[0]!.x) * (u[1]!.y - u[0]!.y)) / 2;
  }
  return a2 > 0 ? Math.sqrt(a3 / a2) : 1;
}

/**
 * The file's meshes, made ready to bind to a Staff's skeleton: positions (which the file may store quantized, undone
 * by its skin) carried to body space at the bind pose, normals to plain floats, bone indices renumbered by name.
 */
function prepareDetailed(cohort: Cohort, gltf: GLTF): DetailedAsset {
  gltf.scene.updateMatrixWorld(true);
  const J = layout(cohort).J;
  const parts: DetailedPart[] = [];
  let spread = 0;
  let jointErr = 0;
  gltf.scene.traverse((o) => {
    const m = o as THREE.SkinnedMesh;
    if (!m.isSkinnedMesh) return;
    const sk = m.skeleton;
    const map = sk.bones.map((b) => BI[b.name]);
    const lacking = sk.bones.filter((_, k) => map[k] === undefined).map((b) => b.name);
    if (lacking.length) throw new Error(`staff ${cohort}: ${m.name} uses bones the rig lacks: ${lacking.join(", ")}`);
    const D = new THREE.Matrix4().multiplyMatrices(sk.bones[0]!.matrixWorld, sk.boneInverses[0]!).multiply(m.bindMatrix);
    for (let k = 0; k < sk.bones.length; k++) {
      _dm.multiplyMatrices(sk.bones[k]!.matrixWorld, sk.boneInverses[k]!).multiply(m.bindMatrix);
      for (let e = 0; e < 16; e++) spread = Math.max(spread, Math.abs(_dm.elements[e]! - D.elements[e]!));
      const j = J[sk.bones[k]!.name];
      if (j) jointErr = Math.max(jointErr, _dp.setFromMatrixPosition(sk.bones[k]!.matrixWorld).distanceTo(j));
    }
    const g = m.geometry;
    const pos = g.getAttribute("position");
    const P = new Float32Array(pos.count * 3);
    for (let i = 0; i < pos.count; i++) _dp.fromBufferAttribute(pos, i).applyMatrix4(D).toArray(P, i * 3);
    g.setAttribute("position", new THREE.BufferAttribute(P, 3));
    const nrm = g.getAttribute("normal");
    if (nrm) {
      _dn.getNormalMatrix(D);
      const N = new Float32Array(nrm.count * 3);
      for (let i = 0; i < nrm.count; i++) _dp.fromBufferAttribute(nrm, i).applyMatrix3(_dn).normalize().toArray(N, i * 3);
      g.setAttribute("normal", new THREE.BufferAttribute(N, 3));
    }
    const si = g.getAttribute("skinIndex");
    const S = new Uint8Array(si.count * 4);
    for (let i = 0; i < si.count; i++) {
      S[i * 4] = map[si.getX(i)]!;
      S[i * 4 + 1] = map[si.getY(i)]!;
      S[i * 4 + 2] = map[si.getZ(i)]!;
      S[i * 4 + 3] = map[si.getW(i)]!;
    }
    g.setAttribute("skinIndex", new THREE.BufferAttribute(S, 4));
    g.boundingBox = null;
    g.boundingSphere = null;
    parts.push({ role: m.name, geometry: g, source: m.material as THREE.MeshStandardMaterial, mpu: metresPerUv(g) });
  });
  if (spread > 1e-4 || jointErr > 0.002)
    console.warn(`staff ${cohort}: the file's bind pose is off this rig's (matrix spread ${spread.toExponential(1)}, joints ${(jointErr * 1000).toFixed(1)} mm)`);
  return { parts };
}

/** Skin: light wrapping past the terminator, reddened, for the scatter under the surface. */
function skinScatter(m: THREE.MeshPhysicalMaterial): THREE.MeshPhysicalMaterial {
  m.onBeforeCompile = (sh) => {
    sh.fragmentShader = sh.fragmentShader.replace(
      "#include <lights_physical_pars_fragment>",
      THREE.ShaderChunk.lights_physical_pars_fragment.replace(
        "reflectedLight.directDiffuse += irradiance * BRDF_Lambert( material.diffuseColor );",
        `reflectedLight.directDiffuse += irradiance * BRDF_Lambert( material.diffuseColor );
	{
		float nl = dot( geometryNormal, directLight.direction );
		float wrapped = saturate( ( nl + 0.5 ) / 1.5 ) - saturate( nl );
		reflectedLight.directDiffuse += wrapped * directLight.color * vec3( 0.55, 0.2, 0.12 ) * BRDF_Lambert( material.diffuseColor );
	}`,
      ),
    );
  };
  m.customProgramCacheKey = () => "staff-skin-scatter-1";
  return m;
}

const detailedMats = new Map<string, THREE.Material>();
function cachedMat<T extends THREE.Material>(key: string, make: () => T): T {
  let m = detailedMats.get(key) as T | undefined;
  if (!m) {
    m = make();
    detailedMats.set(key, m);
  }
  return m;
}
function sharpen(t: THREE.Texture | null): THREE.Texture | null {
  if (t) t.anisotropy = 8;
  return t;
}
/** A copy of a shared pattern texture repeating every `tile` metres on a part whose UVs run `mpu` metres a unit. */
function patternFor(t: THREE.Texture | null, mpu: number, tile: number): THREE.Texture | null {
  if (!t) return null;
  const c = t.clone();
  c.repeat.setScalar(mpu / tile);
  c.needsUpdate = true;
  return c;
}

/** The materials of one person's detailed body: shared per cohort where nothing personal varies. */
function detailedMaterials(cohort: Cohort, a: DetailedAsset, L: Look): Record<string, THREE.Material> {
  const M = peopleMaterials();
  const part = (r: string) => a.parts.find((p) => p.role === r);
  const src = (r: string) => part(r)?.source;
  const mpu = (r: string) => part(r)?.mpu ?? 1;
  const out: Record<string, THREE.Material> = {};
  out.body = cachedMat(`${cohort}:skin`, () =>
    skinScatter(
      new THREE.MeshPhysicalMaterial({
        map: sharpen(src("body")?.map ?? null),
        color: cohort === "p50" ? 0xc8c2bc : 0xf2ece6,
        roughness: 0.5,
        metalness: 0,
        specularIntensity: 0.55,
        ior: 1.4,
        bumpMap: patternFor(M.skin.bumpMap, mpu("body"), 0.05),
        bumpScale: 0.1,
        sheen: 0.25,
        sheenRoughness: 0.5,
        sheenColor: new THREE.Color(0x8a3a2c),
      }),
    ),
  );
  out.eyes = cachedMat(`${cohort}:eyes`, () => new THREE.MeshPhysicalMaterial({ map: sharpen(src("eyes")?.map ?? null), color: 0xeae6e2, alphaTest: 0.5, roughness: 0.3, metalness: 0, clearcoat: 1, clearcoatRoughness: 0.03 }));
  // One hair texture serves every cohort; a tint per cohort makes it black, dark brown or brown.
  const hairTint = cohort === "p50" ? 0x6e6a68 : cohort === "p05" ? 0xa49a94 : 0xffffff;
  const strands = (r: string, rough: number) =>
    cachedMat(`${cohort}:${r}`, () =>
      new THREE.MeshPhysicalMaterial({
        map: sharpen(src(r)?.map ?? null),
        color: r === "hair" ? hairTint : 0xffffff,
        normalMap: src(r)?.normalMap ?? null,
        roughness: rough,
        metalness: 0,
        alphaTest: 0.35,
        alphaToCoverage: true,
        side: THREE.DoubleSide,
        sheen: 0.3,
        sheenRoughness: 0.45,
      }),
    );
  out.brows = strands("brows", 0.6);
  out.lashes = strands("lashes", 0.6);
  out.hair = strands("hair", 0.5);
  const suitAo = sharpen(src("whites")?.map ?? null);
  const suitN = sharpen(src("whites")?.normalMap ?? null);
  const twillSheen = { sheen: 0.55, sheenRoughness: 0.75, sheenColor: new THREE.Color(0xffffff) };
  out.whites = cachedMat(`${cohort}:whites`, () =>
    new THREE.MeshPhysicalMaterial({ map: suitAo, color: L.white, normalMap: suitN, roughness: 0.86, metalness: 0, side: THREE.DoubleSide, ...twillSheen }),
  );
  const capWhite = L.cap.getHex() === L.white.getHex();
  out.cap = capWhite
    ? out.whites
    : cachedMat(`${cohort}:cap:${L.cap.getHexString()}`, () =>
        new THREE.MeshPhysicalMaterial({ map: suitAo, color: L.cap, normalMap: suitN, roughness: 0.8, metalness: 0, side: THREE.DoubleSide, sheen: 0.4, sheenRoughness: 0.7, sheenColor: new THREE.Color(0x888c94) }),
      );
  out.buttons = cachedMat(`buttons:${L.button.getHexString()}`, () => new THREE.MeshPhysicalMaterial({ color: L.button, roughness: 0.55, metalness: 0, clearcoat: 0.3, clearcoatRoughness: 0.4 }));
  out.trousers = L.check
    ? cachedMat(`${cohort}:check`, () =>
        new THREE.MeshPhysicalMaterial({
          map: patternFor(M.check.map, mpu("trousers"), 0.056),
          aoMap: suitAo,
          normalMap: suitN,
          roughness: 0.86,
          metalness: 0,
          side: THREE.DoubleSide,
          sheen: 0.3,
          sheenRoughness: 0.8,
          sheenColor: new THREE.Color(0xcccccc),
        }),
      )
    : cachedMat(`${cohort}:plain`, () =>
        new THREE.MeshPhysicalMaterial({ map: suitAo, color: 0x1c1d21, normalMap: suitN, roughness: 0.8, metalness: 0, side: THREE.DoubleSide, sheen: 0.45, sheenRoughness: 0.7, sheenColor: new THREE.Color(0x8a8f99) }),
      );
  out.apron = cachedMat(`apron:${L.apron.getHexString()}`, () =>
    new THREE.MeshPhysicalMaterial({
      color: L.apron,
      normalMap: patternFor(M.apron.normalMap, mpu("apron"), 0.03),
      normalScale: new THREE.Vector2(0.6, 0.6),
      roughness: 0.88,
      metalness: 0,
      side: THREE.DoubleSide,
      sheen: 0.35,
      sheenRoughness: 0.8,
      sheenColor: new THREE.Color(0xd8dce4),
    }),
  );
  out.shoes = cachedMat(`${cohort}:shoes`, () =>
    new THREE.MeshPhysicalMaterial({ map: sharpen(src("shoes")?.map ?? null), normalMap: src("shoes")?.normalMap ?? null, roughness: 0.4, metalness: 0, clearcoat: 0.5, clearcoatRoughness: 0.3 }),
  );
  return out;
}

/** Mesh names by part: the procedural body's suffixes where they mean the same thing. */
const PART_NAME: Record<string, string> = { body: "skin" };
const NO_SHADOW = new Set(["eyes", "brows", "lashes", "buttons"]);

export class Staff {
  readonly group = new THREE.Group();
  readonly bones = new Map<string, THREE.Bone>();
  readonly dims: Anthro;
  readonly cohort: Cohort;
  readonly detail: Detail;
  readonly skeleton: THREE.Skeleton;
  readonly meshes: THREE.SkinnedMesh[] = [];
  /** Socket for a two-handed item (hotel pan, rack): parent the item to it; carry, push and lift move it. */
  readonly item = new THREE.Object3D();
  /**
   * The held item's solid as a signed distance in the socket's frame (negative inside), which `grasp` hands close
   * on. A hotel pan's by default; set it when the socket holds something else, or null for none.
   */
  itemSdf: ((p: THREE.Vector3) => number) | null = hotelPanSdf;
  tris: number;
  /** Resolves once the detailed body is drawn (or straight away when it is off, or the file could not be had). */
  readonly ready: Promise<void>;
  /** Which body is drawn now. */
  body: "procedural" | "detailed" = "procedural";
  /** Triangles by material. */
  readonly trisBy: Record<string, number> = {};
  /** The last pose applied. */
  readonly current = newPose();
  /** Largest range-of-motion excursions in the last pose, degrees beyond the table's limits (0 when inside). */
  readonly rom: Record<string, number> = {};
  /** The measured angle behind each `rom` entry (the side nearer its limit), degrees. */
  readonly romAngle: Record<string, number> = {};
  /** @internal */
  readonly B: Body;
  /** @internal */
  readonly torso: Torso;
  private bonesArr: THREE.Bone[] = [];
  private Pw: THREE.Vector3[] = [];
  private Qw: THREE.Quaternion[] = [];
  private wrapCache = new Map<string, Record<Finger, [number, number, number]>>();
  private placed = { x: 0, z: 0, yaw: 0 };
  private scratch = newPose();
  /** @internal Last arm solve: iteration, forearm twist, wrist flexion, deviation, residual twist (degrees). */
  armDebug: Record<string, number[]> = {};
  /** The grip point each arm's last solve aimed at (its target, turned by any spin), body space. */
  private reachTarget: Record<Side, THREE.Vector3> = { L: V(), R: V() };
  private outfit: Look;
  private shadows: boolean;
  private shared = new Set<THREE.BufferGeometry>();
  private disposed = false;

  constructor(cohort: Cohort, detail: Detail, opts: StaffOptions = {}) {
    this.cohort = cohort;
    this.detail = detail;
    const B = (this.B = layout(cohort));
    this.dims = B.d;
    const T = (this.torso = torsoModel(B));
    const L = (this.outfit = makeLook(cohort, opts));
    // Bones, all at identity bind rotations.
    for (const name of ALL_BONES) {
      const b = new THREE.Bone();
      b.name = name;
      this.bones.set(name, b);
      this.bonesArr.push(b);
      this.Pw.push(V());
      this.Qw.push(new THREE.Quaternion());
    }
    for (const name of ALL_BONES) {
      const b = this.bones.get(name)!;
      const p = parentOf(name);
      const at = B.J[name]!;
      if (p) {
        this.bones.get(p)!.add(b);
        b.position.copy(at).sub(B.J[p]!);
      } else b.position.copy(at);
    }
    const root = this.bones.get("root")!;
    this.group.add(root);
    this.group.add(this.item);
    this.group.name = `staff-${cohort}`;
    this.group.updateMatrixWorld(true);
    this.skeleton = new THREE.Skeleton(this.bonesArr);
    // Sockets: the grip frames carry the grip's orientation.
    for (const side of SIDES) this.bones.get(`grip_${side}`)!.quaternion.copy(B.hands[side].toolQ);
    this.bones.get("claw_grip_L")!.quaternion.copy(B.hands.L.toolQ);
    this.bones.get("pinch_grip_R")!.quaternion.copy(B.hands.R.toolQ);

    const M = peopleMaterials();
    const soups = { whites: new Soup(), apron: new Soup(), trousers: new Soup(), skin: new Soup(), shoes: new Soup(), hair: new Soup() };
    buildJacket(soups.whites, B, T, L, detail);
    buildApron(soups.apron, B, T, L, detail);
    buildTrousers(soups.trousers, B, T, L, detail);
    for (const side of SIDES) buildArmSkin(soups.skin, B, side, L, detail);
    buildNeck(soups.skin, B, L, detail);
    buildHead(soups.skin, soups.whites, soups.hair, B, L, detail);
    buildShoes(soups.shoes, B, L, detail);
    const mats: Record<keyof typeof soups, THREE.Material> = {
      whites: M.whites,
      apron: M.apron,
      trousers: L.check ? M.check : M.plain,
      skin: M.skin,
      shoes: M.shoes,
      hair: M.hair,
    };
    let tris = 0;
    const shadows = (this.shadows = opts.shadows !== false);
    const H = B.topY;
    for (const key of Object.keys(soups) as (keyof typeof soups)[]) {
      const soup = soups[key];
      if (!soup.tris) continue;
      const mesh = new THREE.SkinnedMesh(soup.geometry(), mats[key]);
      mesh.name = `${this.group.name}-${key}`;
      mesh.castShadow = shadows;
      mesh.receiveShadow = shadows;
      mesh.bind(this.skeleton, new THREE.Matrix4());
      // A generous fixed bound, so culling never trims a limb mid-stride.
      mesh.boundingSphere = new THREE.Sphere(V(0, H * 0.5, 0), H * 0.95);
      mesh.frustumCulled = true;
      this.group.add(mesh);
      this.meshes.push(mesh);
      this.trisBy[key] = soup.tris;
      tris += soup.tris;
    }
    this.tris = tris;
    this.place(0, 0, 0);
    this.stand(0);
    this.ready =
      opts.body !== "detailed"
        ? Promise.resolve()
        : loadDetailedBody(cohort).then(
            (a) => this.useDetailed(a),
            (e: unknown) => console.warn(`staff ${cohort}: keeping the procedural body`, e),
          );
  }

  /** Swap the procedural meshes for the detailed body's, bound to the same skeleton. */
  private useDetailed(a: DetailedAsset): void {
    if (this.disposed) return;
    const mats = detailedMaterials(this.cohort, a, this.outfit);
    for (const m of this.meshes) {
      this.group.remove(m);
      m.geometry.dispose();
    }
    this.meshes.length = 0;
    for (const k of Object.keys(this.trisBy)) delete this.trisBy[k];
    const H = this.B.topY;
    let tris = 0;
    for (const p of a.parts) {
      const mat = mats[p.role];
      if (!mat) continue;
      const key = PART_NAME[p.role] ?? p.role;
      const mesh = new THREE.SkinnedMesh(p.geometry, mat);
      mesh.name = `${this.group.name}-${key}`;
      mesh.castShadow = this.shadows && !NO_SHADOW.has(p.role);
      mesh.receiveShadow = this.shadows;
      mesh.bind(this.skeleton, new THREE.Matrix4());
      mesh.boundingSphere = new THREE.Sphere(V(0, H * 0.5, 0), H * 0.95);
      this.group.add(mesh);
      this.meshes.push(mesh);
      this.shared.add(p.geometry);
      const t = (p.geometry.index ? p.geometry.index.count : p.geometry.getAttribute("position").count) / 3;
      this.trisBy[key] = (this.trisBy[key] ?? 0) + t;
      tris += t;
    }
    this.tris = tris;
    this.body = "detailed";
  }

  /** The grip frame where a held tool's handle goes (the grip_L / grip_R bone). */
  hand(side: Side): THREE.Object3D {
    return this.bones.get(`grip_${side}`)!;
  }
  bone(name: string): THREE.Bone {
    return this.bones.get(name)!;
  }
  get massKg(): number {
    return this.dims.body_mass_kg;
  }
  /** Stand at (x, z) on the floor facing yaw (0 faces +z; positive turns to the person's left). */
  place(x: number, z: number, yaw: number): void {
    this.placed = { x, z, yaw };
    this.group.position.set(x, 0, z);
    this.group.rotation.set(0, yaw, 0);
    this.group.updateMatrix();
  }
  get placement(): { x: number; z: number; yaw: number } {
    return { ...this.placed };
  }
  /** World point to body space (relative to the current placement). */
  toBody(p: THREE.Vector3, out = V()): THREE.Vector3 {
    const c = Math.cos(this.placed.yaw);
    const s = Math.sin(this.placed.yaw);
    const dx = p.x - this.placed.x;
    const dz = p.z - this.placed.z;
    return out.set(c * dx - s * dz, p.y, s * dx + c * dz);
  }
  /** Body-space point to world. */
  toWorld(p: THREE.Vector3, out = V()): THREE.Vector3 {
    const c = Math.cos(this.placed.yaw);
    const s = Math.sin(this.placed.yaw);
    return out.set(this.placed.x + c * p.x + s * p.z, p.y, this.placed.z - s * p.x + c * p.z);
  }

  /* -------------------------------------------------------------------------------------------------------- */
  /* Forward kinematics bookkeeping in body space                                                              */

  private idx(name: string): number {
    return BI[name]!;
  }
  /** Set a bone's local rotation and update its body-space frame from its parent's. */
  private setLocal(name: string, q: THREE.Quaternion): void {
    const i = this.idx(name);
    const b = this.bonesArr[i]!;
    b.quaternion.copy(q);
    const p = parentOf(name);
    if (!p) {
      this.Pw[i]!.copy(b.position);
      this.Qw[i]!.copy(q);
      return;
    }
    const pi = this.idx(p);
    this.Pw[i]!.copy(b.position).applyQuaternion(this.Qw[pi]!).add(this.Pw[pi]!);
    this.Qw[i]!.copy(this.Qw[pi]!).multiply(q);
  }
  /** Set a bone's body-space rotation (its parent's frame must be current). */
  private setWorld(name: string, q: THREE.Quaternion): void {
    const p = parentOf(name)!;
    _q3.copy(this.Qw[this.idx(p)]!).invert().multiply(q);
    this.setLocal(name, _q3);
  }
  /** Body-space position of a bone's joint. */
  jointPos(name: string, out = V()): THREE.Vector3 {
    return out.copy(this.Pw[this.idx(name)]!);
  }
  jointQuat(name: string, out = new THREE.Quaternion()): THREE.Quaternion {
    return out.copy(this.Qw[this.idx(name)]!);
  }

  /* -------------------------------------------------------------------------------------------------------- */
  /* Applying a pose                                                                                           */

  pose(p: BodyPose): void {
    if (p !== this.current) copyPose(this.current, p);
    for (const k of Object.keys(this.rom)) this.rom[k] = 0;
    // Root and pelvis.
    this.setLocal("root", _q1.identity());
    const pel = this.bones.get("pelvis")!;
    pel.position.copy(p.pelvis);
    _q1.setFromAxisAngle(Y, p.yaw).multiply(qAxis(X, p.tilt, _q2)).multiply(qAxis(Z, p.roll, _q3));
    this.setLocal("pelvis", _q1);
    // Spine: flexion mostly lumbar, twist mostly thoracic.
    const fl = [0.48, 0.32, 0.2];
    const tw = [0.2, 0.35, 0.45];
    const names = ["spine_01", "spine_02", "spine_03"];
    for (let k = 0; k < 3; k++) {
      _q1.setFromAxisAngle(Y, p.twist * tw[k]!).multiply(qAxis(X, p.flex * fl[k]!, _q2)).multiply(qAxis(Z, p.side / 3, _q3));
      this.setLocal(names[k]!, _q1);
    }
    this.romNote("torso_flex", (p.flex + 0) / DEG);
    this.romNote("torso_yaw", p.twist / DEG);
    this.setLocal("carry_chest", _q1.identity());
    // Neck and head turn towards the gaze point.
    this.look(p.gaze);
    // Arms.
    for (let s = 0; s < 2; s++) this.arm(SIDES[s]!, p.arms[s]!, p.clav[s * 2]!, p.clav[s * 2 + 1]!);
    // Legs.
    for (let s = 0; s < 2; s++) this.leg(SIDES[s]!, p.feet[s]!);
    // Item socket.
    this.item.visible = true;
    if (p.itemOn) {
      this.item.position.copy(p.itemPos);
      this.item.quaternion.copy(p.itemQuat);
    }
  }

  private romNote(key: RomKey, deg: number): void {
    const [lo, hi] = this.dims.joint_rom_deg[key];
    const over = deg < lo ? lo - deg : deg > hi ? deg - hi : 0;
    if (over >= (this.rom[key] ?? 0)) this.romAngle[key] = deg;
    this.rom[key] = Math.max(this.rom[key] ?? 0, over);
  }

  private look(g: THREE.Vector3): void {
    // Gaze direction in the upper thorax's frame.
    const i3 = this.idx("spine_03");
    const eye = _v1.set(0, this.B.eyeY, 0.08).sub(this.B.J.spine_03!).applyQuaternion(this.Qw[i3]!).add(this.Pw[i3]!);
    const dir = _v2.subVectors(g, eye).applyQuaternion(_q1.copy(this.Qw[i3]!).invert());
    const yaw = Math.atan2(dir.x, dir.z);
    const pitch = Math.atan2(-dir.y, Math.hypot(dir.x, dir.z));
    // The eyes take a share of the turn; the neck and head do the rest.
    const y = clamp(yaw * 0.85, -60 * DEG, 60 * DEG);
    const pt = clamp(pitch * 0.72, -35 * DEG, 50 * DEG);
    this.setLocal("neck", _q1.setFromAxisAngle(Y, y * 0.45).multiply(qAxis(X, pt * 0.55, _q2)));
    this.setLocal("head", _q1.setFromAxisAngle(Y, y * 0.55).multiply(qAxis(X, pt * 0.45, _q2)));
    this.setLocal("gaze_target", _q1.identity());
    this.romNote("neck_yaw", (y * 0.45) / DEG);
    this.romNote("neck_pitch", (pt * 0.55) / DEG);
  }

  /** Grip offset and canonical frame for a grip kind. */
  private gripFrame(side: Side, kind: GripKind): { off: THREE.Vector3; a0: THREE.Vector3; b0: THREE.Vector3 } {
    const H = this.B.hands[side];
    if (kind === "tool") return { off: H.toolPos, a0: V(0, 0, 1).applyQuaternion(H.toolQ), b0: V(0, 1, 0).applyQuaternion(H.toolQ) };
    if (kind === "dagger") return { off: H.toolPos, a0: V(0, 0, -1).applyQuaternion(H.toolQ), b0: V(0, 1, 0).applyQuaternion(H.toolQ) };
    if (kind === "bar") return { off: H.barPos, a0: H.barAxis, b0: H.dn };
    if (kind === "hook") return { off: H.hookPos, a0: H.barAxis, b0: H.dn };
    return { off: H.palmPos, a0: H.da, b0: H.dn };
  }

  private arm(side: Side, a: ArmPose, elev: number, prot: number): void {
    const B = this.B;
    const H = B.hands[side];
    const sd = H.sd;
    // Clavicle.
    _q1.setFromAxisAngle(Z, sd * elev).multiply(qAxis(Y, -sd * prot, _q2));
    this.setLocal(`clavicle_${side}`, _q1);
    const J = B.J;
    const u0 = _v1.subVectors(J[`forearm_${side}`]!, J[`upper_arm_${side}`]!).normalize().clone();
    const f0 = _v1.subVectors(J[`hand_${side}`]!, J[`forearm_${side}`]!).normalize().clone();
    const n0u = V().crossVectors(u0, Z).normalize();
    const n0f = V().crossVectors(f0, Z).normalize();
    const kw = V().crossVectors(H.da, H.dn).normalize();
    // Free swing (FK), as local rotations.
    const fkUp = new THREE.Quaternion().setFromAxisAngle(X, -a.flex).multiply(qAxis(Z, sd * a.abd, _q2)).multiply(qAxis(u0, a.twist, _q3));
    const fkFo = new THREE.Quaternion().setFromAxisAngle(n0f, a.elbow);
    const fkHa = new THREE.Quaternion().setFromAxisAngle(kw, 0.12);
    let upL = fkUp;
    let foL = fkFo;
    let haL = fkHa;
    if (a.ik > 0) {
      const ik = this.solveArm(side, a, u0, f0, n0u, n0f, kw);
      if (a.ik >= 1) {
        upL = ik[0];
        foL = ik[1];
        haL = ik[2];
      } else {
        upL = fkUp.clone().slerp(ik[0], a.ik);
        foL = fkFo.clone().slerp(ik[1], a.ik);
        haL = fkHa.clone().slerp(ik[2], a.ik);
      }
    }
    this.setLocal(`upper_arm_${side}`, upL);
    this.setLocal(`forearm_${side}`, foL);
    this.setLocal(`hand_${side}`, haL);
    // ROM bookkeeping: shoulder flexion/abduction of the upper arm in the thorax, elbow flexion.
    const i3 = this.idx("spine_03");
    const ua = _v3.subVectors(this.Pw[this.idx(`forearm_${side}`)]!, this.Pw[this.idx(`upper_arm_${side}`)]!).normalize().applyQuaternion(_q1.copy(this.Qw[i3]!).invert());
    this.romNote("shoulder_flex", Math.atan2(ua.z, -ua.y) / DEG);
    this.romNote("shoulder_abd", Math.max(0, Math.atan2(sd * ua.x, -ua.y) / DEG));
    const fa = _v4.subVectors(this.Pw[this.idx(`hand_${side}`)]!, this.Pw[this.idx(`forearm_${side}`)]!).normalize();
    const ub = _v3.subVectors(this.Pw[this.idx(`forearm_${side}`)]!, this.Pw[this.idx(`upper_arm_${side}`)]!).normalize();
    this.romNote("elbow_flex", Math.acos(clamp(ub.dot(fa), -1, 1)) / DEG);
    // Fingers.
    this.fingers(side, a.shape, a.close, a.radius);
    // The grip socket carries the tool: its +z is the handle's axis as the fist holds it (towards the working end),
    // its +y the tool's up as asked; a hand that rolled round a round handle does not roll the tool with it.
    if (a.ik > 0 && (a.grip === "tool" || a.grip === "dagger")) {
      const Qhand = this.Qw[this.idx(`hand_${side}`)]!;
      const zh = V(0, 0, a.grip === "dagger" ? -1 : 1).applyQuaternion(H.toolQ).applyQuaternion(Qhand);
      const yh = orthoTo(a.b.clone(), zh);
      const Qt = frameQuat(Z, Y, zh, yh, new THREE.Quaternion());
      const local = Qhand.clone().invert().multiply(Qt);
      this.setLocal(`grip_${side}`, a.ik >= 1 ? local : H.toolQ.clone().slerp(local, a.ik));
    } else this.setLocal(`grip_${side}`, H.toolQ);
    if (side === "L") this.setLocal("claw_grip_L", H.toolQ);
    else this.setLocal("pinch_grip_R", H.toolQ);
  }

  /**
   * Two-bone arm IK to the grip target, inside the joint limits; returns local rotations (upper, fore, hand).
   * The elbow's swivel about the shoulder-wrist line is chosen, near the pose's pole, to keep the forearm's twist
   * and the wrist's flexion and deviation inside their ranges and the elbow out of the trunk, as people do. What
   * still cannot be met is given up in the hand's orientation, never in the joints.
   */
  private solveArm(side: Side, a: ArmPose, u0: THREE.Vector3, f0: THREE.Vector3, n0u: THREE.Vector3, n0f: THREE.Vector3, kw: THREE.Vector3): [THREE.Quaternion, THREE.Quaternion, THREE.Quaternion] {
    const B = this.B;
    const H = B.hands[side];
    const sd = H.sd;
    const gf = this.gripFrame(side, a.grip);
    // The shoulder joint as the clavicle now carries it (its own stored frame is from the last pose).
    const ic = this.idx(`clavicle_${side}`);
    const Qcl = this.Qw[ic]!.clone();
    const S = this.bonesArr[this.idx(`upper_arm_${side}`)]!.position.clone().applyQuaternion(Qcl).add(this.Pw[ic]!);
    const lu = B.upper;
    const lf = B.fore;
    const aT = a.a.clone().normalize();
    const bT = orthoTo(a.b.clone(), aT);
    const a0 = gf.a0.clone().normalize();
    const b0 = orthoTo(gf.b0.clone(), a0);
    const dMin = Math.sqrt(lu * lu + lf * lf - 2 * lu * lf * Math.cos(Math.PI - 140 * DEG));
    const FL: [number, number] = [-62 * DEG, 68 * DEG];
    const DV: [number, number] = [-32 * DEG, 18 * DEG];
    const TW = 85 * DEG;
    const over = (x: number, [lo, hi]: [number, number]) => (x < lo ? lo - x : x > hi ? x - hi : 0);
    // One solve for a hand orientation and an elbow pole.
    let P = a.pos.clone();
    const solve = (Qh: THREE.Quaternion, poleIn: THREE.Vector3) => {
      const W = P.clone().sub(gf.off.clone().applyQuaternion(Qh));
      const dirV = V().subVectors(W, S);
      const short = Math.max(0, dirV.length() - (lu + lf) * 0.9995);
      const d = clamp(dirV.length(), dMin, (lu + lf) * 0.9995);
      const dir = dirV.normalize();
      const alpha = Math.acos(clamp((lu * lu + d * d - lf * lf) / (2 * lu * d), -1, 1));
      let pole = orthoTo(poleIn.clone(), dir, V(sd, 0, 0));
      let E = S.clone().addScaledVector(dir, lu * Math.cos(alpha)).addScaledVector(pole, lu * Math.sin(alpha));
      let inside = 0;
      for (let k = 0; k < 8 && this.insideTrunk(E, 0.052 * B.s); k++) {
        inside++;
        pole = orthoTo(pole.lerp(orthoTo(V(sd, 0, 0), dir), 0.3), dir);
        E = S.clone().addScaledVector(dir, lu * Math.cos(alpha)).addScaledVector(pole, lu * Math.sin(alpha));
      }
      const Wc = S.clone().addScaledVector(dir, d);
      const u = V().subVectors(E, S).normalize();
      const f = V().subVectors(Wc, E).normalize();
      let n = V().crossVectors(u, f);
      if (n.lengthSq() < 1e-6) n = V().crossVectors(u, pole.clone().negate());
      n.normalize();
      const Qup = frameQuat(u0, orthoTo(n0u.clone(), u0), u, n, new THREE.Quaternion());
      const Qfo0 = frameQuat(f0, orthoTo(n0f.clone(), f0), f, n, new THREE.Quaternion());
      // The forearm takes the twist (pronation, supination) the hand needs; the wrist takes the rest.
      const want = wrapPi(twistAngle(Qfo0.clone().invert().multiply(Qh), f0));
      const tau = clamp(want, -TW, TW);
      const Qfo = Qfo0.clone().multiply(qAxis(f0, tau));
      const hl = Qfo.clone().invert().multiply(Qh);
      const v = H.da.clone().applyQuaternion(hl);
      const flex = Math.atan2(v.dot(H.dn), v.dot(H.da));
      // Deviation towards the thumb (radial) positive.
      const dev = Math.atan2(v.dot(H.dt), v.dot(H.da));
      const rest = wrapPi(twistAngle(hl, H.da));
      const cost = over(flex, FL) + over(dev, DV) + Math.abs(want - tau) + Math.max(0, Math.abs(rest) - 10 * DEG) + inside * 0.3 + short * 12;
      return { Qup, Qfo, hl, flex, dev, rest, cost, pole, tau, want, short };
    };
    // Swivel the elbow round the shoulder-wrist line, preferring the pose's own pole; on a round handle the hand
    // may also roll round it.
    let best: (ReturnType<typeof solve> & { Qh: THREE.Quaternion; P: THREE.Vector3 }) | null = null;
    const rolls = a.roll > 0 ? [0, 30, -30, 60, -60, 90, -90, 120, -120, 150, -150, 180].filter((r) => Math.abs(r) * DEG <= a.roll + 1e-6) : [0];
    const swivels = [0, 12, -12, 25, -25, 40, -40, 55, -55, 70, -70, 90, -90, 110];
    const trial = (rd: number, sp: number, deg: number) => {
      const qs = qAxis(a.spinAxis, sp * DEG);
      const aS = aT.clone().applyQuaternion(qs);
      const bS = bT.clone().applyQuaternion(qs).applyAxisAngle(aS, rd * DEG);
      P = a.pos.clone().sub(a.spinPivot).applyQuaternion(qs).add(a.spinPivot).addScaledVector(a.spinShift, Math.sin(sp * DEG));
      const Qh0 = frameQuat(a0, b0, aS, bS, new THREE.Quaternion());
      const W0 = P.clone().sub(gf.off.clone().applyQuaternion(Qh0));
      const axis = V().subVectors(W0, S).normalize();
      const r = solve(Qh0, a.pole.clone().applyAxisAngle(axis, sd * deg * DEG));
      r.cost += Math.abs(deg) * DEG * 0.08 + Math.abs(rd) * DEG * 0.03 + sp * DEG * 0.04;
      return { ...r, Qh: Qh0, P: P.clone() };
    };
    // The spin, where allowed, lays the fingers along the forearm (a straight wrist): a smooth function of the
    // target, found by a few fixed-point passes, so the hand never flicks between samples.
    let spin = 0;
    if (a.spinMax > 0) {
      const ax = a.spinAxis.clone().normalize();
      const e1 = orthoTo(aT.clone(), ax);
      const e2 = V().crossVectors(ax, e1);
      for (let it = 0; it < 4; it++) {
        const r = trial(0, spin, 0);
        const f = f0.clone().applyQuaternion(r.Qfo);
        const fp = f.addScaledVector(ax, -f.dot(ax));
        spin = clamp(Math.atan2(fp.dot(e2), fp.dot(e1)) / DEG, 0, a.spinMax / DEG);
      }
    }
    const spins = [spin];
    // Every roll and swivel sample is scored; the choice is their soft minimum (cost-weighted mean), so the elbow
    // and the hand move continuously as the target moves, never jumping between samples.
    const samples: [number, number, number][] = [];
    let cMin = Infinity;
    for (const rd of rolls)
      for (const sp of spins)
        for (const deg of swivels) {
          const c = trial(rd, sp, deg).cost;
          samples.push([rd, deg, c]);
          cMin = Math.min(cMin, c);
        }
    const TEMP = 0.04;
    let wS = 0;
    let rS = 0;
    let dS = 0;
    for (const [rd, deg, c] of samples) {
      const w = Math.exp(-(c - cMin) / TEMP);
      wS += w;
      rS += w * rd;
      dS += w * deg;
    }
    const bestAt: [number, number, number] = [rS / wS, spin, dS / wS];
    best = trial(bestAt[0], spin, bestAt[2]);
    P = best!.P;
    this.reachTarget[side].copy(P);
    const Qh0 = best!.Qh;
    let { Qup, Qfo, hl } = best!;
    const pole = best!.pole;
    // If the wrist is still out of range, bring it inside and let the hand follow; re-solve the reach.
    let Qh = Qh0;
    for (let iter = 0; iter < 6; iter++) {
      const r = iter === 0 ? best! : solve(Qh, pole);
      Qup = r.Qup;
      Qfo = r.Qfo;
      hl = r.hl;
      const flexC = clamp(r.flex, FL[0], FL[1]);
      const devC = clamp(r.dev, DV[0], DV[1]);
      const restC = clamp(r.rest, -10 * DEG, 10 * DEG);
      if (flexC === r.flex && devC === r.dev && restC === r.rest) break;
      hl = qAxis(kw, flexC).multiply(qAxis(H.dn, sd * devC)).multiply(qAxis(H.da, restC));
      Qh = Qfo.clone().multiply(hl);
    }
    this.armDebug[side] = [best!.cost / DEG, best!.flex / DEG, best!.dev / DEG, best!.rest / DEG, best!.tau / DEG, best!.want / DEG, best!.short * 1000, bestAt[1], bestAt[2]];
    const vf = H.da.clone().applyQuaternion(hl);
    this.romNote("wrist_flex", Math.atan2(vf.dot(H.dn), vf.dot(H.da)) / DEG);
    const upL = Qcl.clone().invert().multiply(Qup);
    const foL = Qup.clone().invert().multiply(Qfo);
    return [upL, foL, hl];
  }

  /** Is a body-space point within `r` of the trunk's surface or inside it? */
  insideTrunk(p: THREE.Vector3, r: number): boolean {
    // Bring the point into the bind frame of the nearest spine segment.
    const seg = p.y > this.Pw[this.idx("spine_03")]!.y - 0.02 ? "spine_03" : p.y > this.Pw[this.idx("spine_02")]!.y ? "spine_02" : "spine_01";
    const i = this.idx(seg);
    const q = _q2.copy(this.Qw[i]!).invert();
    const b = _v4.subVectors(p, this.Pw[i]!).applyQuaternion(q).add(this.B.J[seg]!);
    if (b.y < this.B.hemY || b.y > this.B.neckBaseY) return false;
    const S = this.torso.sect(b.y);
    const zz = b.z - S.c;
    const dz = zz >= 0 ? S.f + r : S.b + r;
    return Math.abs(b.x) ** S.n / (S.a + r) ** S.n + Math.abs(zz) ** S.n / dz ** S.n < 1;
  }

  private leg(side: Side, f0: FootPose): void {
    // The ankle's dorsiflexion (the shin's lean over the foot) has a limit (20 degrees in the table): where a pose
    // would pass it, the heel comes up and the foot turns on its ball joint, the toes staying down, as people do in
    // a deep crouch or a pivot.
    const lim = this.dims.joint_rom_deg.ankle_pitch[1] - 2;
    const df = this.legSolve(side, f0, false);
    if (df <= lim) {
      this.legSolve(side, f0, true);
      return;
    }
    const B = this.B;
    const f = this.liftHeel;
    // Raising the heel by an angle takes back about 0.6 of it at the ankle (the knee straightens a little too).
    const raise = ((df - lim) / 0.6) * DEG;
    // Ankle about the ball: in the foot's frame the ankle sits mtpFwd behind and (ankY - toe height) above it.
    const ball = V(0, B.J[`toe_${side}`]!.y - B.ankY, B.mtpFwd).applyQuaternion(f0.quat).add(f0.pos);
    f.quat.copy(f0.quat).multiply(qAxis(X, raise, _q1));
    f.pos.copy(V(0, B.ankY - B.J[`toe_${side}`]!.y, -B.mtpFwd).applyQuaternion(f.quat)).add(ball);
    f.toe = f0.toe + raise;
    this.legSolve(side, f, true);
  }
  private liftHeel: FootPose = { pos: V(), quat: new THREE.Quaternion(), toe: 0 };

  /** Two-bone leg to the foot; with `apply` false only measures the ankle's dorsiflexion it would need (degrees). */
  private legSolve(side: Side, f: FootPose, apply: boolean): number {
    const B = this.B;
    const J = B.J;
    const lt = B.thigh;
    const ls = B.shin;
    const ip = this.idx("pelvis");
    const Hp = _v1.subVectors(J[`thigh_${side}`]!, J.pelvis!).applyQuaternion(this.Qw[ip]!).add(this.Pw[ip]!).clone();
    const A = f.pos;
    const dirV = V().subVectors(A, Hp);
    let d = dirV.length();
    const dir = dirV.normalize();
    d = clamp(d, Math.abs(lt - ls) + 0.05, (lt + ls) * 0.99995);
    const cosA = clamp((lt * lt + d * d - ls * ls) / (2 * lt * d), -1, 1);
    const alpha = Math.acos(cosA);
    // The knee tracks over the foot (a touch outward).
    const fwd = Z.clone().applyQuaternion(f.quat);
    const pole = orthoTo(fwd.addScaledVector(X, sgn(side) * 0.12), dir);
    const K = Hp.clone().addScaledVector(dir, lt * Math.cos(alpha)).addScaledVector(pole, lt * Math.sin(alpha));
    const Ac = Hp.clone().addScaledVector(dir, d);
    const t = V().subVectors(K, Hp).normalize();
    const sv = V().subVectors(Ac, K).normalize();
    let n = V().crossVectors(t, sv);
    if (n.lengthSq() < 1e-6) n = V().crossVectors(t, pole.clone().negate());
    n.normalize();
    const t0 = V().subVectors(J[`shin_${side}`]!, J[`thigh_${side}`]!).normalize();
    const s0 = V().subVectors(J[`foot_${side}`]!, J[`shin_${side}`]!).normalize();
    const nt0 = orthoTo(V().crossVectors(t0, Z.clone().negate()), t0);
    const ns0 = orthoTo(V().crossVectors(s0, Z.clone().negate()), s0);
    if (!apply) {
      const qs = frameQuat(s0, ns0, sv, n, new THREE.Quaternion()).invert();
      const fl = Z.clone().applyQuaternion(f.quat).applyQuaternion(qs);
      return Math.atan2(fl.y, fl.z) / DEG;
    }
    this.setWorld(`thigh_${side}`, frameQuat(t0, nt0, t, n, _q1));
    this.setWorld(`shin_${side}`, frameQuat(s0, ns0, sv, n, _q1));
    this.setWorld(`foot_${side}`, f.quat);
    this.setLocal(`toe_${side}`, qAxis(X, -f.toe, _q1));
    // The apron's cloth bone: plumb from the knee, unless the shin swings forward of plumb and pushes the cloth.
    const fwdT = Z.clone().applyQuaternion(this.Qw[this.idx(`thigh_${side}`)]!);
    const plumb = V(0, -1, 0);
    const push = smooth(-0.04, 0.04, sv.dot(fwdT) - plumb.dot(fwdT));
    const hang = plumb.lerp(sv, push).normalize();
    _q1.setFromUnitVectors(t, hang).multiply(this.Qw[this.idx(`thigh_${side}`)]!);
    this.setWorld(`apron_${side}`, _q1);
    // ROM: hip flexion against the pelvis, knee flexion, ankle pitch against the shin.
    const tp = t.clone().applyQuaternion(_q2.copy(this.Qw[ip]!).invert());
    this.romNote("hip_flex", Math.atan2(tp.z, -tp.y) / DEG);
    this.romNote("knee_flex", Math.acos(clamp(t.dot(sv), -1, 1)) / DEG);
    const fl = Z.clone().applyQuaternion(f.quat).applyQuaternion(_q2.copy(this.Qw[this.idx(`shin_${side}`)]!).invert());
    const df = Math.atan2(fl.y, fl.z) / DEG;
    this.romNote("ankle_pitch", df);
    return df;
  }

  /** Finger angles for a shape, solved against the handle for tool, bar and hook grips. */
  private fingers(side: Side, shape: HandShape, close: number, radius: number): void {
    const H = this.B.hands[side];
    let target: Record<Finger, [number, number, number]> = RELAX;
    if (shape === "claw") target = CLAW;
    else if (shape === "flat") target = FLAT;
    else if (shape === "tool" || shape === "bar" || shape === "hook") target = this.wrap(side, shape, radius);
    else if (shape === "grasp") target = this.grasp(side);
    const c = clamp(close, 0, 1);
    for (const f of FINGERS) {
      const g = H.digits[f];
      const t = target[f];
      // A grasp opens from its contact towards a flat hand, never bending past the contact (into the item).
      const r = shape === "grasp" ? (FLAT[f].map((x, k) => Math.min(x, t[k]!)) as [number, number, number]) : RELAX[f];
      for (let k = 0; k < 3; k++) this.setLocal(`${f}_0${k + 1}_${side}`, qAxis(g.k, lerp(r[k]!, t[k]!, c), _q1));
    }
  }

  /** Close each digit round a cylinder until it touches, joint by joint (computed once per handle size). */
  private wrap(side: Side, kind: "tool" | "bar" | "hook", radius: number): Record<Finger, [number, number, number]> {
    const key = `${side}${kind}${radius.toFixed(4)}`;
    const hit = this.wrapCache.get(key);
    if (hit) return hit;
    const H = this.B.hands[side];
    const gf = this.gripFrame(side, kind);
    const O = H.wrist.clone().add(gf.off);
    const D = kind === "tool" ? V(0, 0, 1).applyQuaternion(H.toolQ) : H.barAxis.clone();
    const r = kind === "hook" ? Math.min(radius, 0.006) : radius;
    const dist = (p: THREE.Vector3) => {
      const w = _v4.subVectors(p, O);
      return w.addScaledVector(D, -w.dot(D)).length();
    };
    const out = {} as Record<Finger, [number, number, number]>;
    for (const f of FINGERS) {
      const g = H.digits[f];
      const thumb = f === "thumb";
      const pts = (a: [number, number, number]) => {
        const P: THREE.Vector3[] = [g.base.clone()];
        let ang = 0;
        for (let k = 0; k < 3; k++) {
          ang += a[k]!;
          P.push(P[k]!.clone().addScaledVector(g.dir.clone().applyAxisAngle(g.k, ang), g.len[k]!));
        }
        return P;
      };
      const touches = (a: [number, number, number], from: number) => {
        const P = pts(a);
        for (let k = from; k < 3; k++)
          for (let m = 1; m <= 4; m++) {
            const q = P[k]!.clone().lerp(P[k + 1]!, m / 4);
            if (dist(q) < r + g.rad * (m === 4 && k === 2 ? 0.75 : 0.92)) return true;
          }
        return false;
      };
      const max: [number, number, number] = thumb ? [0.75, 0.8, 0.95] : kind === "hook" ? [0.95, 1.5, 0.9] : [1.45, 1.75, 1.25];
      const a: [number, number, number] = [0, 0, 0];
      {
        // Stage 1: all joints together until something touches.
        let u = 0;
        if (!touches(a, 0)) {
          for (; u <= 1; u += 0.02) if (touches([max[0] * u, max[1] * u, max[2] * u], 0)) break;
          u = Math.max(0, u - 0.02);
        }
        a[0] = max[0] * u;
        a[1] = max[1] * u;
        a[2] = max[2] * u;
        // Stage 2: the middle and end joints wrap on.
        for (let v = 0; v < 60 && a[1] < max[1] * 1.08 && !touches([a[0], a[1] + 0.03, a[2] + 0.02], 1); v++) {
          a[1] += 0.03;
          a[2] += 0.02;
        }
        // Stage 3: the tip curls in.
        for (let v = 0; v < 40 && a[2] < max[2] * 1.1 && !touches([a[0], a[1], a[2] + 0.03], 2); v++) a[2] += 0.03;
      }
      out[f] = a;
    }
    this.wrapCache.set(key, out);
    return out;
  }

  /**
   * Close each digit on the held item, as the hand now lies against it: a finger bends at its knuckle until its pad
   * meets the surface, then its middle and end joints until they do (curling round an edge where there is one); a
   * thumb first opens until it clears the item (over a lip, round an edge), then curls down onto it. Computed from
   * the item's signed distance every pose, so the fingers follow the hand's real contact, never a stored shape.
   */
  private grasp(side: Side): Record<Finger, [number, number, number]> {
    const sdf = this.itemSdf;
    const c = this.current;
    if (!sdf || !c.itemOn) return RELAX;
    const H = this.B.hands[side];
    const ih = this.idx(`hand_${side}`);
    const qi = new THREE.Quaternion().copy(c.itemQuat).invert();
    const R = qi.clone().multiply(this.Qw[ih]!);
    const T = this.Pw[ih]!.clone().sub(c.itemPos).applyQuaternion(qi);
    const Jh = this.B.J[`hand_${side}`]!;
    const tmp = V();
    const dist = (p: THREE.Vector3) => sdf(tmp.copy(p).sub(Jh).applyQuaternion(R).add(T));
    const out = {} as Record<Finger, [number, number, number]>;
    const seg = [V(), V(), V(), V()];
    const dir = V();
    const q = V();
    for (const f of FINGERS) {
      const g = H.digits[f];
      const thumb = f === "thumb";
      const place = (a: readonly number[]) => {
        seg[0]!.copy(g.base);
        let ang = 0;
        for (let k = 0; k < 3; k++) {
          ang += a[k]!;
          dir.copy(g.dir).applyAxisAngle(g.k, ang);
          seg[k + 1]!.copy(seg[k]!).addScaledVector(dir, g.len[k]!);
        }
      };
      // Does any part of the digit from segment `from` on come within its own radius of the surface?
      const touches = (a: readonly number[], from: number) => {
        place(a);
        for (let k = from; k < 3; k++)
          for (let m = 1; m <= 6; m++) {
            q.copy(seg[k]!).lerp(seg[k + 1]!, m / 6);
            // The skin stands about 1.05 radii off the bone line on the pad side (less at the domed tip).
            if (dist(q) < g.rad * (k === 2 && m >= 5 ? (m === 6 ? 0.55 : 1.0) : 1.15)) return true;
          }
        return false;
      };
      const max: [number, number, number] = thumb ? [0.75, 0.85, 0.95] : [1.45, 1.75, 1.25];
      const a: [number, number, number] = thumb ? [0.3, 0.12, 0.1] : [0, 0.12, 0.08];
      if (thumb) {
        // Open the thumb (towards the back of the hand) until it is clear of the item.
        while (a[0] > -0.95 && touches(a, 0)) a[0] -= 0.05;
      } else {
        // A finger pressed into the surface at rest straightens back out of it.
        while (a[0] > -0.35 && touches(a, 0)) a[0] -= 0.04;
        // Knuckle first (the middle and end joints follow a little) until the finger meets the surface.
        for (let v = 0; v < 60 && a[0] < max[0] && !touches([a[0] + 0.025, a[1] + 0.008, a[2] + 0.005], 0); v++) {
          a[0] += 0.025;
          a[1] += 0.008;
          a[2] += 0.005;
        }
      }
      // The middle and end joints wrap on, then the tip curls in, each until it touches.
      for (let v = 0; v < 60 && a[1] < max[1] && !touches([a[0], a[1] + 0.03, a[2] + 0.02], 1); v++) {
        a[1] += 0.03;
        a[2] += 0.02;
      }
      for (let v = 0; v < 40 && a[2] < max[2] && !touches([a[0], a[1], a[2] + 0.03], 2); v++) a[2] += 0.03;
      out[f] = a;
    }
    return out;
  }

  /* -------------------------------------------------------------------------------------------------------- */
  /* Readouts                                                                                                  */

  /** Trunk flexion the body shows: the L5/S1-to-C7 line's lean from vertical, forward positive, degrees. */
  lumbarFlexDeg(): number {
    const l5 = this.Pw[this.idx("spine_01")]!;
    const c7 = this.Pw[this.idx("neck")]!;
    const v = _v1.subVectors(c7, l5);
    // Forward is the pelvis's facing.
    const fwd = _v2.set(0, 0, 1).applyQuaternion(this.Qw[this.idx("pelvis")]!);
    fwd.y = 0;
    fwd.normalize();
    return Math.atan2(v.dot(fwd), v.y) / DEG;
  }
  /**
   * What the L5/S1 estimate needs from the pose: trunk flexion, the hands' horizontal reach from L5/S1, and the held
   * item's: the horizontal distance from L5/S1 to `itemCom` (a point in the item socket's frame, its centre of mass;
   * default the socket's origin), or the hands' reach when nothing is held.
   */
  liftMetrics(itemCom?: THREE.Vector3): { torsoFlexDeg: number; reachM: number; handsY: number; loadReachM: number } {
    const l5 = this.Pw[this.idx("spine_01")]!;
    const gl = this.gripPoint("L", _v3);
    const gr = this.gripPoint("R", _v4);
    const mid = gl.add(gr).multiplyScalar(0.5);
    const reachM = Math.hypot(mid.x - l5.x, mid.z - l5.z);
    let loadReachM = reachM;
    const c = this.current;
    if (c.itemOn) {
      const com = _v1.copy(itemCom ?? ZERO).applyQuaternion(c.itemQuat).add(c.itemPos);
      loadReachM = Math.hypot(com.x - l5.x, com.z - l5.z);
    }
    return { torsoFlexDeg: this.lumbarFlexDeg(), reachM, handsY: mid.y, loadReachM };
  }
  /** How far (metres) a reaching hand's grip point is from where its last solve aimed; 0 for a free arm. */
  gripError(side: Side): number {
    const a = this.current.arms[side === "L" ? 0 : 1]!;
    if (a.ik < 1) return 0;
    const i = this.idx(`hand_${side}`);
    const off = this.gripFrame(side, a.grip).off;
    return off.clone().applyQuaternion(this.Qw[i]!).add(this.Pw[i]!).distanceTo(this.reachTarget[side]);
  }
  /** Body-space grip point of a hand (tool grip). */
  gripPoint(side: Side, out = V()): THREE.Vector3 {
    const i = this.idx(`hand_${side}`);
    return out.copy(this.B.hands[side].toolPos).applyQuaternion(this.Qw[i]!).add(this.Pw[i]!);
  }

  /* -------------------------------------------------------------------------------------------------------- */
  /* Motions at the current placement                                                                          */

  /** Stand at ease: slow weight shifts and breathing (t in seconds). */
  stand(t = 0): void {
    const p = standPose(this, t, this.scratch);
    fitHips(this, p);
    this.pose(p);
  }
  /**
   * Work at a station at the current placement, t seconds into the work. Targets in `opts.at` are world points
   * (the pan's centre for saute, the pot rim's centre for stir, the plate, board, rack bar or shelved pan's centre).
   */
  work(kind: WorkKind, t: number, opts: WorkOpts = {}): void {
    const p = workPose(this, kind, t, opts, this.scratch);
    fitHips(this, p);
    this.pose(p);
  }

  dispose(): void {
    this.disposed = true;
    for (const m of this.meshes) if (!this.shared.has(m.geometry)) m.geometry.dispose();
    this.skeleton.dispose();
  }
}


/* ------------------------------------------------------------------------------------------------------------ */
/* Standing and working at a station. Poses are in body space: the root between the feet, facing +z.             */
/* ------------------------------------------------------------------------------------------------------------ */

const TOE_OUT = 7 * DEG;

/**
 * The tools the work motions are fitted to, in their grip frames (metres): +z runs along the handle towards the
 * working end, +y is up (the thumb side). Build props to these numbers and parent them to `hand(side)`.
 *   Skillet: the handle runs from z = -handleBack to the rim at z = gripFromRim, rising at `rise` from the pan;
 *     the pan's base (radius `base`) is `depth` below its rim (radius `rim`).
 *   Knife: handle from z = -handleBack to the bolster at z = bolster; the edge runs straight at y = edgeY from the
 *     heel (z = heelZ) to z = bellyZ, then curves up on a circle of radius bellyR to the tip at (tipY, tipZ).
 *   Ladle: a straight handle from z = -handleBack to the bowl, whose centre is `reach` along z. Stirring holds
 *     it in the dagger grip (working end below the little finger); the socket's +z still points to the bowl.
 *   Spoon (plating): the bowl's centre is `reach` along z.
 */
export const TOOLS = {
  skillet: { rim: 0.14, base: 0.11, depth: 0.05, rise: 20 * DEG, gripFromRim: 0.11, handleBack: 0.075 },
  knife: { handleBack: 0.065, bolster: 0.06, heelZ: 0.075, bellyZ: 0.17, tipZ: 0.28, edgeY: -0.04, tipY: -0.02, bellyR: 0.3125, spineY: 0.008 },
  ladle: { reach: 0.5, bowlR: 0.045, handleBack: 0.06 },
  spoon: { reach: 0.17, handleBack: 0.05 },
} as const;


/** Both feet flat at the standing stance. */
function standFeet(st: Staff, p: BodyPose): void {
  const B = st.B;
  for (let i = 0; i < 2; i++) {
    const sd = i === 0 ? 1 : -1;
    const f = p.feet[i]!;
    f.pos.set(sd * B.standHalf, B.ankY, 0);
    f.quat.setFromAxisAngle(Y, sd * TOE_OUT);
    f.toe = 0;
  }
}

const _fq = new THREE.Quaternion();
const _fv = V();
/**
 * Keep the hips low enough for both legs to reach their feet with the knees just short of locked: the pelvis
 * height becomes a soft minimum of the wanted height and each leg's reach. Returns how far it moved (<= 0).
 */
function fitHips(st: Staff, p: BodyPose, k = 0.012, weight: readonly [number, number] = [1, 1]): number {
  const B = st.B;
  const J = B.J;
  const Lmax = (B.thigh + B.shin) * 0.9985;
  _fq.setFromAxisAngle(Y, p.yaw).multiply(qAxis(X, p.tilt, _q1)).multiply(qAxis(Z, p.roll, _q2));
  let py = p.pelvis.y;
  for (let i = 0; i < 2; i++) {
    const side = SIDES[i]!;
    const off = _fv.subVectors(J[`thigh_${side}`]!, J.pelvis!).applyQuaternion(_fq);
    const f = p.feet[i]!.pos;
    const d2 = (p.pelvis.x + off.x - f.x) ** 2 + (p.pelvis.z + off.z - f.z) ** 2;
    const reach = f.y + Math.sqrt(Math.max(0, Lmax * Lmax - d2)) - off.y;
    py = smin(py, reach + (1 - clamp(weight[i]!, 0, 1)) * 0.3, k * B.s);
  }
  const dy = py - p.pelvis.y;
  p.pelvis.y = py;
  return dy;
}

function relaxArms(p: BodyPose): void {
  for (const a of p.arms) {
    a.flex = 0.06;
    a.abd = 0.05;
    a.elbow = 0.24;
    a.twist = 0;
    a.ik = 0;
    a.roll = 0;
    a.spinMax = 0;
    a.shape = "relax";
    a.close = 1;
  }
}

/** Standing at ease: weight drifting slowly from foot to foot, breathing. */
export function standPose(st: Staff, t: number, p: BodyPose = newPose()): BodyPose {
  const B = st.B;
  standFeet(st, p);
  const w = Math.sin((t * 2 * Math.PI) / 7.3 + 0.7 * (st.dims.standing_height_m * 10));
  p.pelvis.set(0.012 * w * B.s, B.hipY - (0.0015 + 0.0015 * Math.abs(w)) * B.s, 0.006 * B.s);
  p.tilt = 0.04;
  // Contrapposto: the loaded side's hip rides higher, the shoulders tilt the other way.
  p.roll = 0.012 * w;
  p.yaw = 0.01 * Math.sin(t * 0.5);
  p.flex = 0.02 + 0.006 * Math.sin((t * 2 * Math.PI) / 4.2);
  p.side = -0.02 * w;
  p.twist = 0;
  p.clav = [0.01, 0, 0.01, 0];
  relaxArms(p);
  p.gaze.set(0.3 * Math.sin(t * 0.23), B.eyeY - 0.15, 3);
  p.itemOn = false;
  return p;
}

function setTool(a: ArmPose, pos: THREE.Vector3, axis: THREE.Vector3, up: THREE.Vector3, pole: THREE.Vector3, radius = 0.012, grip: GripKind = "tool", roll = 0): void {
  a.ik = 1;
  a.grip = grip;
  a.roll = roll;
  a.spinMax = 0;
  a.pos.copy(pos);
  a.a.copy(axis).normalize();
  a.b.copy(up).normalize();
  a.pole.copy(pole).normalize();
  a.shape = "tool";
  a.close = 1;
  a.radius = radius;
}
function setPalm(a: ArmPose, pos: THREE.Vector3, fingers: THREE.Vector3, facing: THREE.Vector3, pole: THREE.Vector3, shape: HandShape = "flat", close = 1): void {
  a.ik = 1;
  a.roll = 0;
    a.spinMax = 0;
  a.grip = "palm";
  a.pos.copy(pos);
  a.a.copy(fingers).normalize();
  a.b.copy(facing).normalize();
  a.pole.copy(pole).normalize();
  a.shape = shape;
  a.close = close;
}

/** World `at` (if given) to body space, else the default. */
function target(st: Staff, o: WorkOpts, def: THREE.Vector3): THREE.Vector3 {
  return o.at ? st.toBody(o.at) : def;
}

/**
 * A full-size hotel pan (GN 1/1: 530 x 325 mm over the flange, 100 mm deep), its frame at the rim's centre, long side
 * across the body (x). The walls stand `flange` inside the flange's edge at the top and lean in by `taper` (radians)
 * going down; the rolled flange is `lip` thick. (smallwares.ts builds NSF-HOTEL-PAN-FULL-4IN to these numbers.)
 */
export const HOTEL_PAN = { length: 0.53, width: 0.325, depth: 0.1, flange: 0.0125, gripHalf: 0.27, taper: 0.155, lip: 0.0045 } as const;

/**
 * Signed distance (metres) from a point in the item frame to a hotel pan's solid: its tapered body (the food space
 * counted as inside, so nothing reaches into it) and its flange. Negative inside.
 */
export function hotelPanSdf(p: THREE.Vector3): number {
  const K = HOTEL_PAN;
  const ct = Math.cos(K.taper);
  const inset = Math.max(0, -p.y) * Math.tan(K.taper);
  const dx = Math.abs(p.x) - (K.length / 2 - K.flange - inset);
  const dz = Math.abs(p.z) - (K.width / 2 - K.flange - inset);
  const dy = Math.max(-K.depth - p.y, p.y);
  const ex = Math.max(dx * ct, 0);
  const ez = Math.max(dz * ct, 0);
  const ey = Math.max(dy, 0);
  const body = Math.hypot(ex, ez, ey) + Math.min(Math.max(dx * ct, dz * ct, dy), 0);
  const fx = Math.abs(p.x) - K.length / 2;
  const fz = Math.abs(p.z) - K.width / 2;
  const fy = Math.max(-K.lip - p.y, p.y);
  const flange = Math.hypot(Math.max(fx, 0), Math.max(fz, 0), Math.max(fy, 0)) + Math.min(Math.max(fx, fz, fy), 0);
  return Math.min(body, flange);
}

function carryFrame(st: Staff, bob: number, pos: THREE.Vector3, q: THREE.Quaternion): void {
  const B = st.B;
  const y = B.waistY - 0.03 * B.s;
  // The pan's near wall clears the apron at the belly by a couple of centimetres.
  const belly = Math.max(st.torso.at(y - 0.05, 0).z, st.torso.at(y - 0.1, 0).z);
  pos.set(0, y + bob, belly + 0.012 + HOTEL_PAN.width / 2 + 0.025);
  q.setFromAxisAngle(X, -0.03);
}

/**
 * Hands at the pan's ends as cooks hold a full one. Carried at the waist (forearms level), each palm lies flat on
 * the outside of an end wall just under the rolled lip, fingers along the wall towards the far side, thumb over the
 * lip and down on it. Reaching down to it (forearms steep), the fingers tip down the wall and the palm rises and
 * stands out onto the lip's edge, the fingers curling in under the lip to the wall. Between the two the arm chooses
 * (a turn about the wall's normal that keeps the wrist straight). `half` moves the walls for a pan of another
 * length (the default is HOTEL_PAN.gripHalf, the full pan). The fingers close on the pan by contact (`grasp`).
 */
function panGrips(st: Staff, p: BodyPose, itemPos: THREE.Vector3, itemQ: THREE.Quaternion, half: number, close: number, ik: number): void {
  const K = HOTEL_PAN;
  const H = st.B.hands.L;
  const sh = H.sh;
  const tilt0 = 12 * DEG;
  // How far the hand stands above its palm point at a tilt: the thumb side of the palm and wrist (the thumb itself
  // goes over the lip).
  const top = (phi: number) => Math.max(0.0335 * sh * Math.cos(phi) + 0.22 * H.HL * Math.sin(phi), 0.0272 * sh * Math.cos(phi) + 0.36 * H.HL * Math.sin(phi));
  const y = -(K.lip + 0.0035) - top(tilt0);
  const tanT = Math.tan(K.taper);
  const wallX = K.length / 2 - K.flange - -y * tanT + (half - K.gripHalf);
  const ct = Math.cos(K.taper);
  const stT = Math.sin(K.taper);
  for (let i = 0; i < 2; i++) {
    const sd = i === 0 ? 1 : -1;
    const a = p.arms[i]!;
    const n = V(sd * ct, -stT, 0);
    const down = V(-sd * stT, -ct, 0);
    const fingers = V(0, 0, Math.cos(tilt0)).addScaledVector(down, Math.sin(tilt0));
    // The palm point a few millimetres off the wall: the palm's skin and the thenar's swell meet the steel.
    const g = V(sd * wallX, y, -0.055).addScaledVector(n, 0.0045);
    setPalm(
      a,
      g.clone().applyQuaternion(itemQ).add(itemPos),
      fingers.applyQuaternion(itemQ),
      n.clone().negate().applyQuaternion(itemQ),
      V(sd * 0.55, -0.8, -0.25),
      "grasp",
      close,
    );
    a.ik = ik;
    // Fingers from level along the wall to straight down it; the palm rises to the lip and stands out past its
    // edge (by the lip's overhang and a little), so that it never rides into the lip on the way.
    const tip = 78 * DEG;
    a.spinMax = tip;
    // A positive turn about this axis tips the fingers down the wall, on either hand.
    a.spinAxis.copy(n).multiplyScalar(sd).applyQuaternion(itemQ);
    a.spinPivot.copy(g).applyQuaternion(itemQ).add(itemPos);
    const rise = -y - 0.004;
    const out = K.flange + -y * tanT + 0.006;
    a.spinShift.set(sd * out, rise, 0).multiplyScalar(1 / Math.sin(tip)).applyQuaternion(itemQ);
  }
}

/** Hands on a push bar (a speed rack's handle), palms forward and down. */
function barGrips(st: Staff, p: BodyPose, barPos: THREE.Vector3, barQ: THREE.Quaternion): void {
  const half = 0.2 * st.B.s;
  for (let i = 0; i < 2; i++) {
    const sd = i === 0 ? 1 : -1;
    const a = p.arms[i]!;
    a.ik = 1;
    // A round bar: the fist may roll round it to keep the wrist straight.
    a.roll = (75 * Math.PI) / 180;
    a.spinMax = 0;
    a.grip = "bar";
    a.pos.copy(V(sd * half, 0, 0).applyQuaternion(barQ).add(barPos));
    a.a.copy(V(-sd, 0, 0).applyQuaternion(barQ));
    a.b.copy(V(0, -0.45, 0.89).applyQuaternion(barQ));
    a.pole.set(sd * 0.55, -0.8, -0.2).normalize();
    a.shape = "bar";
    a.close = 1;
    a.radius = 0.0125;
  }
}

/**
 * The crouch that brings the hands to a low target: trunk inclination, knee bend and hips back, from c in [0, 1].
 * The shin leans at most 18 degrees (ankle dorsiflexion 20), the knee bends to 85, the trunk leans 50 from vertical.
 */
/**
 * Lower the body for a reach down: c = 0 stands, c = 1 is a half squat with the knees at 85 degrees, the ankles at
 * 18 degrees of dorsiflexion (near their limit) and the trunk leaning 50 degrees; beyond 1 (to 1.3) the legs stay
 * and the back stoops further, as tall people do over a low shelf.
 */
function crouch(st: Staff, p: BodyPose, c: number): void {
  const B = st.B;
  const cl = Math.min(c, 1);
  const cx = Math.max(0, c - 1);
  const phi = 18 * DEG * cl;
  const kappa = 85 * DEG * cl;
  const ank = V(0, B.ankY, 0);
  const knee = ank.clone().add(V(0, Math.cos(phi), Math.sin(phi)).multiplyScalar(B.shin));
  const th = kappa - phi;
  const hip = knee.clone().add(V(0, Math.cos(th), -Math.sin(th)).multiplyScalar(B.thigh));
  // At c = 0 this is the straight standing leg; keep the knees just unlocked.
  p.pelvis.set(p.pelvis.x, Math.min(hip.y, B.hipY - 0.0015 * B.s), hip.z + 0.006 * B.s);
  const lean = 50 * DEG * cl;
  p.tilt = 0.04 + 0.55 * lean + 0.3 * 50 * DEG * cx;
  p.flex = 0.02 + 0.45 * lean + 0.7 * 50 * DEG * cx;
}

const crouchCache = new WeakMap<Staff, Map<string, number>>();

/**
 * The least crouch that brings both grip points (body space, left then right) within 95% of the arms' reach,
 * found by posing the body itself (so the spine's bend and the hips' fit are as they will be), then restoring it.
 * Cached per person and goal.
 */
function crouchFor(st: Staff, goals: readonly [THREE.Vector3, THREE.Vector3]): number {
  const key = goals.map((g) => g.toArray().map((x) => x.toFixed(3)).join(",")).join(";");
  let m = crouchCache.get(st);
  if (!m) crouchCache.set(st, (m = new Map()));
  const hit = m.get(key);
  if (hit !== undefined) return hit;
  const B = st.B;
  // The wrist within 93% of the arm's length, the palm point (0.36 of the hand's length on) beyond it.
  const reach = (B.upper + B.fore) * 0.93 + 0.36 * B.d.hand_length_m;
  const saved = copyPose(newPose(), st.current);
  const tmp = newPose();
  const ok = (c: number) => {
    standPose(st, 0, tmp);
    tmp.pelvis.x *= 0.4;
    crouch(st, tmp, c);
    fitHips(st, tmp);
    st.pose(tmp);
    return st.jointPos("upper_arm_L").distanceTo(goals[0]) <= reach && st.jointPos("upper_arm_R").distanceTo(goals[1]) <= reach;
  };
  let hi = 0;
  if (!ok(0)) {
    let lo = 0;
    hi = 1.3;
    for (let k = 0; k < 16; k++) {
      const mid = (lo + hi) / 2;
      if (ok(mid)) hi = mid;
      else lo = mid;
    }
  }
  st.armDebug.crouch = [hi, st.jointPos("upper_arm_L").distanceTo(goals[0]), reach];
  st.pose(saved);
  m.set(key, hi);
  return hi;
}

/** A work pose at time t (seconds into the work) at the current placement. */
export function workPose(st: Staff, kind: WorkKind, t: number, o: WorkOpts = {}, p: BodyPose = newPose()): BodyPose {
  const B = st.B;
  const s = B.s;
  standPose(st, t, p);
  const tempo = o.tempo ?? 1;
  const tt = t * tempo;
  const R = p.arms[1]!;
  const Lh = p.arms[0]!;
  p.pelvis.x *= 0.4;
  switch (kind) {
    case "stand":
      break;
    case "saute": {
      // The pan's base sits on the grate; the handle comes back towards the right hand, rising from the rim.
      const K = TOOLS.skillet;
      const y0 = o.surfaceY ?? 0.92;
      const base = target(st, o, V(-0.03 * s, y0, 0.5 * s));
      const hh = V(0.38, 0, 0.92).normalize();
      // A toss every few seconds: slide out, snap back and up with the far lip rising; small shakes between.
      const P = 2.7;
      const c = tt % P;
      const toss = c < 0.8 ? c / 0.8 : -1;
      let slide = 0.008 * Math.sin(tt * 2 * Math.PI * 1.7) * smooth(0.95, 1.4, c) * (1 - smooth(2.3, 2.6, c));
      let lift = 0;
      let pitch = 0;
      if (toss >= 0) {
        const out = Math.sin(Math.PI * smooth(0, 0.38, toss));
        const snap = Math.sin(Math.PI * smooth(0.28, 0.95, toss));
        slide = 0.04 * out - 0.05 * snap;
        lift = 0.03 * snap;
        pitch = 15 * DEG * snap;
      }
      const hd = hh.clone().multiplyScalar(Math.cos(K.rise)).addScaledVector(Y, -Math.sin(K.rise));
      const up = orthoTo(Y.clone(), hd);
      const rimNear = base.clone().addScaledVector(Y, K.depth).addScaledVector(hh, -K.rim);
      const grip = rimNear.clone().addScaledVector(hd, -K.gripFromRim);
      // Tilt the whole pan about its near rim (the far lip rises), then slide and lift it.
      const lat = V().crossVectors(hh, Y).normalize();
      const rel = grip.sub(rimNear).applyAxisAngle(lat, pitch);
      hd.applyAxisAngle(lat, pitch);
      up.applyAxisAngle(lat, pitch);
      const g = rimNear.add(rel).addScaledVector(hh, slide).addScaledVector(Y, lift + 0.004 * Math.sin(Math.PI * smooth(0, 1, Math.abs(pitch) / (15 * DEG))));
      setTool(R, g, hd, up, V(-0.6, -0.7, -0.35), 0.011, "tool", (70 * Math.PI) / 180);
      p.itemOn = false;
      // The other hand rests on the range's front rail.
      setPalm(Lh, V(0.24 * s, y0 - 0.01, Math.min(base.z - 0.2, 0.3 * s)), V(-0.3, -0.1, 1), V(0, -1, 0), V(0.6, -0.6, -0.5));
      p.flex += 8 * DEG;
      p.twist = -0.05;
      p.pelvis.z += 0.012 * s;
      p.clav[3] = 0.05;
      p.gaze.copy(base).add(V(0, 0.03, 0));
      break;
    }
    case "stir": {
      // A ladle round a stockpot: the hand circles above the rim, the bowl circles near the bottom.
      const K = TOOLS.ladle;
      const rim = target(st, o, V(0.02 * s, (o.surfaceY ?? 0.6) + 0.33, 0.5 * s));
      const w = (tt * 2 * Math.PI) / 2.3;
      // The hand circles just above the near rim, the bowl near the bottom; the handle crosses the rim's plane
      // inside the pot, so it never rides over the wall.
      // A fist over the pot, the ladle running down out of its little-finger side to the bowl near the bottom,
      // leaning back a little towards the cook; the hand may roll round the round handle to keep the wrist easy.
      const grip = rim.clone().add(V(-0.03 + 0.03 * Math.cos(w + 0.5), 0.27 + 0.006 * Math.sin(2 * w), -0.03 + 0.03 * Math.sin(w + 0.5)));
      const aim = rim.clone().add(V(0.055 * Math.cos(w), -0.2, 0.055 * Math.sin(w) - 0.05));
      const axis = aim.sub(grip).normalize();
      setTool(R, grip, axis, orthoTo(V(0, 0.3, -1), axis), V(-0.85, -0.45, -0.2), 0.01, "dagger", (120 * Math.PI) / 180);
      // The other hand on the hip, fingers forward, elbow out.
      const hipSide = st.torso.at(B.waistY - 0.03 * s, Math.PI / 2 - 0.45).add(V(0.016, 0, 0));
      setPalm(Lh, hipSide, V(0.05, -0.45, 0.89), V(-1, 0, 0), V(1, -0.2, -0.25), "relax", 1);
      p.flex += 10 * DEG;
      p.twist = -0.03;
      p.gaze.copy(rim).add(V(0, -0.1, 0));
      void K;
      break;
    }
    case "plate": {
      const y0 = o.surfaceY ?? 0.9;
      const plate = target(st, o, V(0, y0 + 0.015, 0.4 * s));
      // Spoon placements round the plate, one every 2.2 s: travel, lower, place, lift.
      const P = 2.2;
      const k = Math.floor(tt / P);
      const u = (tt % P) / P;
      const spot = (n: number) => plate.clone().add(V(0.055 * Math.cos(n * 2.4) - 0.02, 0, 0.055 * Math.sin(n * 2.4)));
      const a = spot(k);
      const b = spot(k + 1);
      const move = mj(smooth(0.55, 0.95, u));
      const bowl = a.clone().lerp(b, move);
      bowl.y += 0.02 + 0.05 * (1 - Math.sin(Math.PI * smooth(0.05, 0.5, u))) * (u < 0.5 ? 1 : 0) + 0.05 * Math.sin(Math.PI * smooth(0.5, 1, u)) * (u >= 0.5 ? 1 : 0);
      // The spoon held upright in the fist, its bowl below the little finger.
      const axis = V(-0.1, -0.98, -0.05).normalize();
      const grip = bowl.clone().addScaledVector(axis, -TOOLS.spoon.reach);
      setTool(R, grip, axis, orthoTo(V(0, 0.3, -1), axis), V(-0.85, -0.45, -0.2), 0.006, "dagger", (120 * Math.PI) / 180);
      R.shape = "tool";
      // Left hand steadies the plate at its rim.
      const rimL = plate.clone().add(V(0.115, 0.006, -0.05));
      setPalm(Lh, rimL.add(V(0.03, 0.012, -0.02)), V(-0.55, -0.15, 0.82), V(0.15, -1, 0.1), V(0.6, -0.65, -0.45), "flat");
      p.flex += 15 * DEG;
      p.pelvis.z += 0.018 * s;
      p.gaze.copy(plate);
      break;
    }
    case "chop": {
      // The rocking chop: the curved belly rolls on the board while the handle rises and falls; the edge never
      // goes below the board. The left hand's knuckles guide the blade's flat and step back along the food.
      const K = TOOLS.knife;
      const y0 = o.surfaceY ?? 0.9;
      const board = target(st, o, V(-0.03 * s, y0 + 0.02, 0.43 * s));
      const f = 2.2;
      const ph = tt * f;
      const cut = Math.floor(ph);
      const fr = ph - cut;
      const rock = 0.5 + 0.5 * Math.cos(2 * Math.PI * fr);
      const kd = V(0.3, 0, 0.954).normalize();
      const left = V(kd.z, 0, -kd.x);
      const n9 = cut % 9;
      const feed = n9 < 8 ? 0.005 * (n9 + mj(clamp(fr * 2.5, 0, 1))) : 0.04 * (1 - mj(fr));
      const pitch = lerp(2, 22, rock) * DEG;
      const axis = kd.clone().multiplyScalar(Math.cos(pitch)).addScaledVector(Y, -Math.sin(pitch));
      const up = orthoTo(Y.clone(), axis);
      // The belly's circle centre sits bellyR above the board, rolling forward a little as the handle drops.
      const roll = K.bellyR * (pitch - 12 * DEG) * 0.35;
      const O = board.clone().addScaledVector(Y, K.bellyR).addScaledVector(kd, 0.02 + roll).addScaledVector(left, feed - 0.02);
      const grip = O.addScaledVector(axis, -K.bellyZ).addScaledVector(up, -(K.edgeY + K.bellyR));
      setTool(R, grip, axis, up, V(-0.55, -0.6, -0.55), 0.012, "tool", (30 * Math.PI) / 180);
      // Claw: fingertips on the food just left of the blade, knuckles against its flat, wrist out to the left.
      const contact = board.clone().addScaledVector(kd, 0.02 + roll + K.bellyZ - K.bellyZ * 0.55).addScaledVector(left, feed - 0.02);
      const da = V().addScaledVector(left, -0.62).addScaledVector(kd, 0.45).addScaledVector(Y, -0.5).normalize();
      const knuck = contact.addScaledVector(left, 0.014).addScaledVector(Y, 0.034);
      const palm = knuck.addScaledVector(da, -0.075).addScaledVector(Y, 0.012);
      setPalm(Lh, palm, da, V(0.25, -1, -0.1), V(0.85, -0.45, -0.3), "claw", 1);
      p.flex += 13 * DEG;
      p.pelvis.z += 0.015 * s;
      p.twist = -0.03;
      p.gaze.copy(board).add(V(0.03, 0, 0.05));
      break;
    }
    case "carry": {
      const bob = 0.004 * Math.sin(tt * 2.1);
      carryFrame(st, bob, p.itemPos, p.itemQuat);
      p.itemOn = true;
      panGrips(st, p, p.itemPos, p.itemQuat, o.halfWidth ?? HOTEL_PAN.gripHalf, 1, 1);
      p.flex -= 2 * DEG;
      p.gaze.set(0, B.eyeY - 0.3, 2.5);
      break;
    }
    case "push": {
      const h = o.surfaceY ?? 1.02;
      const bar = target(st, o, V(0, h, 0.5 * s));
      p.itemOn = true;
      p.itemPos.copy(bar);
      p.itemQuat.identity();
      barGrips(st, p, bar, p.itemQuat);
      p.flex += 10 * DEG;
      p.tilt += 0.06;
      p.pelvis.z += 0.03 * s;
      p.gaze.set(0, B.eyeY - 0.2, 3);
      break;
    }
    case "lift": {
      // Lift a pan from a low shelf to the carry, hold it, put it back: an 8.6 s cycle.
      const shelf = target(st, o, V(0, (o.surfaceY ?? 0.32) + 0.1, 0.43 * s));
      const half = o.halfWidth ?? HOTEL_PAN.gripHalf;
      const cyc = tt % 8.6;
      const shelfQ = new THREE.Quaternion();
      const carryP = V();
      const carryQ = new THREE.Quaternion();
      carryFrame(st, 0, carryP, carryQ);
      // The palms on the end walls at the shelf: the lowest the hands go.
      const lip = (sd: number) => shelf.clone().add(V(sd * (HOTEL_PAN.length / 2 - HOTEL_PAN.flange - 0.003), -0.05, -0.055));
      const cLow = crouchFor(st, [lip(1), lip(-1)]);
      let c = 0;
      let reachW = 0;
      let close = 0;
      let carryU = 0;
      if (cyc < 1) {
        c = 0;
      } else if (cyc < 2.6) {
        const u = ease((cyc - 1) / 1.6);
        c = cLow * u;
        // The hands go to the pan once the body is most of the way down.
        reachW = smooth(0.45, 0.95, u);
      } else if (cyc < 3) {
        c = cLow;
        reachW = 1;
        close = smooth(2.6, 2.95, cyc);
      } else if (cyc < 4.6) {
        const u = (cyc - 3) / 1.6;
        carryU = mj(u);
        reachW = 1;
        close = 1;
        c = cLow * (1 - smooth(0.12, 1, u));
      } else if (cyc < 5.6) {
        carryU = 1;
        reachW = 1;
        close = 1;
      } else if (cyc < 7.2) {
        const u = (cyc - 5.6) / 1.6;
        carryU = 1 - mj(u);
        reachW = 1;
        close = 1;
        c = cLow * smooth(0, 0.88, u);
      } else if (cyc < 7.6) {
        c = cLow;
        reachW = 1;
        close = 1 - smooth(7.2, 7.55, cyc);
      } else {
        // Hands off the pan first, then up.
        const u = (cyc - 7.6) / 1.0;
        c = cLow * (1 - ease(clamp((u - 0.2) / 0.8, 0, 1)));
        reachW = 1 - smooth(0, 0.3, u);
      }
      // The pan rises straight off the shelf first, then comes in close to the body.
      const itemP = shelf.clone().lerp(carryP, carryU);
      itemP.y += 0.03 * Math.sin(Math.PI * clamp(carryU * 2.5, 0, 1)) * (carryU < 0.4 ? 1 : 0);
      itemP.z = lerp(shelf.z, carryP.z, smooth(0.15, 1, carryU));
      const itemQ = shelfQ.clone().slerp(carryQ, carryU);
      p.itemOn = true;
      p.itemPos.copy(itemP);
      p.itemQuat.copy(itemQ);
      crouch(st, p, c);
      panGrips(st, p, itemP, itemQ, half, close, reachW);
      // Arms hang down towards the pan as the body lowers, before the hands take it.
      for (const a of p.arms) {
        a.flex = 0.06 + 0.9 * c;
        a.elbow = 0.25;
      }
      // Eyes on the pan, but standing tall the eyes do more of the looking down than the neck.
      const tall = cLow > 1e-3 ? 1 - clamp(c / cLow, 0, 1) : 1;
      p.gaze.copy(itemP).lerp(V(0, B.eyeY - 0.45 * s, 1.8 * s), 0.5 * tall * (1 - carryU)).lerp(V(0, B.eyeY, 3), carryU * 0.6);
      break;
    }
  }
  return p;
}


/* ------------------------------------------------------------------------------------------------------------ */
/* Walking: a footstep plan along a smoothed path, feet planted where they land, the body carried over them.     */
/* ------------------------------------------------------------------------------------------------------------ */

/** A floor path in the scene's xz plane, its corners rounded, measured by arc length. Yaw 0 faces +z. */
export class FloorPath {
  readonly xs: number[] = [];
  readonly zs: number[] = [];
  readonly ss: number[] = [];
  readonly length: number;
  constructor(points: readonly (readonly [number, number])[], fillet = 0.5) {
    const P = points.filter((p, i) => i === 0 || Math.hypot(p[0] - points[i - 1]![0], p[1] - points[i - 1]![1]) > 1e-4);
    const push = (x: number, z: number) => {
      const n = this.xs.length;
      if (n && Math.hypot(x - this.xs[n - 1]!, z - this.zs[n - 1]!) < 1e-5) return;
      this.xs.push(x);
      this.zs.push(z);
    };
    push(P[0]![0], P[0]![1]);
    for (let i = 1; i < P.length - 1; i++) {
      const [ax, az] = P[i - 1]!;
      const [bx, bz] = P[i]!;
      const [cx, cz] = P[i + 1]!;
      const l0 = Math.hypot(bx - ax, bz - az);
      const l1 = Math.hypot(cx - bx, cz - bz);
      const d0x = (bx - ax) / l0;
      const d0z = (bz - az) / l0;
      const d1x = (cx - bx) / l1;
      const d1z = (cz - bz) / l1;
      const turn = Math.acos(clamp(d0x * d1x + d0z * d1z, -1, 1));
      if (turn < 0.02) {
        push(bx, bz);
        continue;
      }
      const tl = Math.min(fillet * Math.tan(turn / 2), 0.45 * Math.min(l0, l1));
      const r = tl / Math.tan(turn / 2);
      const sx = bx - d0x * tl;
      const sz = bz - d0z * tl;
      // Centre of the arc, on the inside of the turn.
      const cross = d0x * d1z - d0z * d1x;
      const nx = cross > 0 ? -d0z : d0z;
      const nz = cross > 0 ? d0x : -d0x;
      const ox = sx + nx * r;
      const oz = sz + nz * r;
      const a0 = Math.atan2(sx - ox, sz - oz);
      const sgnT = cross > 0 ? -1 : 1;
      const n = Math.max(2, Math.ceil(turn / (5 * DEG)));
      for (let k = 0; k <= n; k++) {
        const a = a0 + (sgnT * turn * k) / n;
        push(ox + Math.sin(a) * r, oz + Math.cos(a) * r);
      }
    }
    push(P[P.length - 1]![0], P[P.length - 1]![1]);
    let s = 0;
    this.ss.push(0);
    for (let i = 1; i < this.xs.length; i++) {
      s += Math.hypot(this.xs[i]! - this.xs[i - 1]!, this.zs[i]! - this.zs[i - 1]!);
      this.ss.push(s);
    }
    this.length = s;
  }
  private seg(s: number): number {
    const ss = this.ss;
    let lo = 0;
    let hi = ss.length - 1;
    while (hi - lo > 1) {
      const m = (lo + hi) >> 1;
      if (ss[m]! <= s) lo = m;
      else hi = m;
    }
    return lo;
  }
  private dirYaw(i: number): number {
    const j = clamp(i, 0, this.xs.length - 2);
    return Math.atan2(this.xs[j + 1]! - this.xs[j]!, this.zs[j + 1]! - this.zs[j]!);
  }
  /** Position and heading at arc length s. */
  at(s: number): { x: number; z: number; yaw: number } {
    if (this.xs.length < 2) return { x: this.xs[0] ?? 0, z: this.zs[0] ?? 0, yaw: 0 };
    const c = clamp(s, 0, this.length);
    const i = this.seg(c);
    const u = (c - this.ss[i]!) / Math.max(1e-9, this.ss[i + 1]! - this.ss[i]!);
    // Heading blends between neighbouring segments' directions so it turns smoothly along arcs.
    const y0 = this.dirYaw(i);
    const yPrev = this.dirYaw(i - 1);
    const yNext = this.dirYaw(i + 1);
    const ya = y0 + wrapPi(yPrev - y0) * 0.5 * (1 - u) * (i > 0 ? 1 : 0);
    const yb = y0 + wrapPi(yNext - y0) * 0.5 * u * (i < this.xs.length - 2 ? 1 : 0);
    return { x: lerp(this.xs[i]!, this.xs[i + 1]!, u), z: lerp(this.zs[i]!, this.zs[i + 1]!, u), yaw: ya + (yb - y0) };
  }
}

interface Plant {
  x: number;
  z: number;
  yaw: number;
  tOn: number;
  tOff: number;
  arrive: "walk" | "step" | "init";
  leave: "walk" | "step" | "end";
  hsAng: number;
  hsDur: number;
  toAng: number;
  toDur: number;
  /** Seconds over which the foot takes up the body's weight after landing, and gives it up before lifting. */
  rUp: number;
  rDn: number;
}

/** What the hands are doing while walking or turning: swinging free, carrying a pan, pushing a rack. */
export type Upper = "free" | "carry" | "push";

export type RoutineStep =
  | { kind: "walk"; path: readonly (readonly [number, number])[]; speed?: number; upper?: Upper; opts?: WorkOpts }
  | { kind: "turn"; yaw: number; upper?: Upper; opts?: WorkOpts }
  | { kind: "work"; work: WorkKind; duration: number; opts?: WorkOpts }
  | { kind: "stand"; duration: number };

interface Seg {
  t0: number;
  t1: number;
  type: "walk" | "turn" | "work" | "stand";
  x: number;
  z: number;
  yaw0: number;
  yaw1: number;
  path?: FloorPath;
  prog?: Mono;
  speed?: number;
  stepLen?: number;
  upper: Upper;
  work?: WorkKind;
  opts?: WorkOpts;
}

export interface RoutineState {
  activity: "walk" | "turn" | "work" | "stand";
  work?: WorkKind;
  /** Seconds into the current activity. */
  local: number;
  /** Whether the item socket is in use (a pan carried or lifted, a rack pushed). */
  holding: boolean;
}

const BIG = 1e9;

/**
 * A person's timeline: walks, turns in place and work at stations, planned once into footsteps so that every pose
 * is a pure function of time. `apply(t)` places and poses the person; with `loop` the timeline repeats seamlessly
 * (the routine must end where it started, as `routeRoutine` arranges).
 */
export class Routine {
  readonly duration: number;
  private plants: Record<Side, Plant[]> = { L: [], R: [] };
  private segs: Seg[] = [];
  private cx: number;
  private cz: number;
  private cyaw: number;
  private t = 0;
  private lastStep: Side = "L";
  private pose = newPose();
  private aux = newPose();

  constructor(readonly staff: Staff, start: { x: number; z: number; yaw: number }, steps: readonly RoutineStep[], readonly loop = true) {
    this.cx = start.x;
    this.cz = start.z;
    this.cyaw = start.yaw;
    for (const side of SIDES) {
      const [x, z] = this.stanceAt(this.cx, this.cz, this.cyaw, side);
      this.plants[side].push(this.plant(x, z, this.cyaw + sgn(side) * TOE_OUT, -BIG, "init"));
    }
    for (const st of steps) {
      if (st.kind === "turn") this.planTurn(st.yaw, st.upper ?? "free", st.opts);
      else if (st.kind === "walk") this.planWalk(st.path, st.speed ?? 1.0, st.upper ?? "free", st.opts);
      else {
        const dur = st.duration;
        if (dur <= 0) continue;
        const work = st.kind === "work" ? st.work : undefined;
        const upper: Upper = work === "carry" || work === "push" ? work : "free";
        this.segs.push({ t0: this.t, t1: this.t + dur, type: st.kind === "work" ? "work" : "stand", x: this.cx, z: this.cz, yaw0: this.cyaw, yaw1: this.cyaw, upper, work, opts: st.kind === "work" ? st.opts : undefined });
        this.t += dur;
      }
    }
    if (!this.segs.length) this.segs.push({ t0: 0, t1: 1, type: "stand", x: this.cx, z: this.cz, yaw0: this.cyaw, yaw1: this.cyaw, upper: "free" });
    this.duration = Math.max(this.t, 1e-3);
    // Weight transfer: a foot takes the load while the other is still down (double support) and hands it over
    // before it lifts, so the ramps are the overlaps of the two feet's plants.
    for (const side of SIDES) {
      const other = this.plants[side === "L" ? "R" : "L"];
      for (const p of this.plants[side]) {
        if (p.tOn > -BIG / 2) {
          const q = other.find((o) => o.tOn <= p.tOn && p.tOn <= o.tOff);
          p.rUp = q && q.tOff < BIG / 2 ? clamp(q.tOff - p.tOn, 0.06, 0.3) : 0.3;
        }
        if (p.tOff < BIG / 2) {
          let q: Plant | undefined;
          for (const o of other) if (o.tOn <= p.tOff) q = o;
          p.rDn = q && q.tOn > -BIG / 2 ? clamp(p.tOff - q.tOn, 0.06, 0.3) : 0.24;
        }
      }
    }
  }

  private stanceAt(x: number, z: number, yaw: number, side: Side, half = this.staff.B.standHalf): [number, number] {
    const sd = sgn(side);
    return [x + Math.cos(yaw) * sd * half, z - Math.sin(yaw) * sd * half];
  }
  private plant(x: number, z: number, yaw: number, tOn: number, arrive: Plant["arrive"]): Plant {
    return { x, z, yaw, tOn, tOff: BIG, arrive, leave: "end", hsAng: 0, hsDur: 0, toAng: 0, toDur: 0, rUp: 0.3, rDn: 0.3 };
  }
  private last(side: Side): Plant {
    const a = this.plants[side];
    return a[a.length - 1]!;
  }
  /** Lift a foot at tLift and set it down at tLand on a new plant. */
  private step(side: Side, tLift: number, leave: Plant["leave"], np: Plant): void {
    const lp = this.last(side);
    lp.tOff = tLift;
    lp.leave = leave;
    this.plants[side].push(np);
    this.lastStep = side;
  }

  private planTurn(yaw: number, upper: Upper = "free", opts?: WorkOpts): void {
    const d = wrapPi(yaw - this.cyaw);
    if (Math.abs(d) < 0.01) {
      this.cyaw = yaw;
      return;
    }
    const m = Math.max(1, Math.ceil(Math.abs(d) / (62 * DEG)));
    const lead: Side = d > 0 ? "L" : "R";
    const trail: Side = lead === "L" ? "R" : "L";
    const t0 = this.t;
    const sw = 0.34;
    const dsT = 0.14;
    const pre = 0.2;
    for (let j = 1; j <= m; j++) {
      const yj = this.cyaw + (d * j) / m;
      for (const [k, side] of [[0, lead], [1, trail]] as [number, Side][]) {
        const lift = t0 + pre + (j - 1) * 2 * (sw + dsT) + k * (sw + dsT);
        const [x, z] = this.stanceAt(this.cx, this.cz, yj, side);
        this.step(side, lift, "step", this.plant(x, z, yj + sgn(side) * TOE_OUT, lift + sw, "step"));
      }
    }
    const t1 = t0 + pre + m * 2 * (sw + dsT) + 0.2;
    this.segs.push({ t0, t1, type: "turn", x: this.cx, z: this.cz, yaw0: this.cyaw, yaw1: this.cyaw + d, upper, opts });
    this.t = t1;
    this.cyaw = this.cyaw + d;
  }

  private planWalk(pts: readonly (readonly [number, number])[], v: number, upper: Upper, opts?: WorkOpts): void {
    const B = this.staff.B;
    const s = B.s;
    const path = new FloorPath([[this.cx, this.cz], ...pts], 0.55 * s);
    const D = path.length;
    if (D < 0.04) return;
    const h0 = path.at(0).yaw;
    if (Math.abs(wrapPi(h0 - this.cyaw)) > 0.03) this.planTurn(h0, upper, opts);
    const end = path.at(D);
    if (D < 0.32 * s) {
      // A shuffle: two short steps to the new spot.
      const t0 = this.t;
      const lead: Side = this.lastStep === "L" ? "R" : "L";
      for (const [k, side] of [[0, lead], [1, lead === "L" ? "R" : "L"]] as [number, Side][]) {
        const [x, z] = this.stanceAt(end.x, end.z, this.cyaw, side);
        const lift = t0 + 0.2 + k * 0.5;
        this.step(side, lift, "step", this.plant(x, z, this.cyaw + sgn(side) * TOE_OUT, lift + 0.36, "step"));
      }
      // The root moves between the stances over the shuffle, carried by the feet.
      this.segs.push({ t0, t1: t0 + 1.3, type: "turn", x: this.cx, z: this.cz, yaw0: this.cyaw, yaw1: this.cyaw, path, upper, opts });
      this.t = t0 + 1.3;
      this.cx = end.x;
      this.cz = end.z;
      return;
    }
    // Step length from the walk ratio (step length over cadence, about 0.0066 m per step/min for adults).
    const L0 = clamp(Math.sqrt(v * 60 * 0.006 * s), 0.3 * s, 0.78 * s);
    const n = Math.max(2, Math.round(D / L0 + 0.5));
    const L = D / (n - 0.5);
    const T = L / v;
    const ds = 0.2 * T;
    const sw = T - ds;
    const vk = clamp(v / 1.3, 0.2, 1.2);
    const hsAng = (7 + 7 * vk) * DEG;
    // The foot rolls over the ball until it stands steep at toe-off (about 50 degrees at a brisk walk).
    const toAng = (32 + 20 * vk) * DEG;
    const lead: Side = this.lastStep === "L" ? "R" : "L";
    const other: Side = lead === "L" ? "R" : "L";
    const t0 = this.t;
    const ws = 0.24;
    const side = (k: number): Side => (k % 2 === 1 ? lead : other);
    const hs: number[] = [0];
    const rock = (p: Plant, leaving: boolean, amp = 1) => {
      if (leaving) {
        p.toAng = toAng * amp;
        p.toDur = 0.55 * T;
      }
    };
    // The lead foot leaves its standing plant after the weight shifts onto the other.
    const first = this.last(lead);
    first.leave = "walk";
    first.tOff = t0 + ws;
    rock(first, true, 0.55);
    for (let k = 1; k <= n; k++) {
      const hk = k === 1 ? t0 + ws + sw * 0.9 : hs[k - 1]! + T;
      hs.push(hk);
      const a = (k - 0.5) * L;
      const at = path.at(a);
      const sd = sgn(side(k));
      const half = k === n ? B.standHalf : B.walkHalf;
      const x = at.x + Math.cos(at.yaw) * sd * half;
      const z = at.z - Math.sin(at.yaw) * sd * half;
      const np = this.plant(x, z, at.yaw + sd * TOE_OUT * (k === n ? 1 : 0.5), hk, "walk");
      np.hsAng = hsAng * (k === 1 ? 0.7 : 1);
      np.hsDur = 0.2 * T;
      // The foot this one replaces lifted ds after the previous heel strike.
      const prev = this.last(side(k));
      if (k >= 2) {
        prev.tOff = hs[k - 1]! + ds;
        prev.leave = "walk";
        rock(prev, true, k === 2 ? 0.8 : 1);
      }
      this.plants[side(k)].push(np);
    }
    // The closing step brings the trailing foot alongside.
    const closer = side(n - 1);
    const hn = hs[n]!;
    const tl = hn + ds;
    const hc = tl + sw * 0.85;
    const [cx2, cz2] = this.stanceAt(end.x, end.z, end.yaw, closer);
    const cp = this.last(closer);
    cp.tOff = tl;
    cp.leave = "walk";
    rock(cp, true, 0.7);
    const np = this.plant(cx2, cz2, end.yaw + sgn(closer) * TOE_OUT, hc, "walk");
    np.hsAng = hsAng * 0.4;
    np.hsDur = 0.16 * T;
    this.plants[closer].push(np);
    this.lastStep = closer;
    // Body progress along the path: midway between the feet at each heel strike, at rest at both ends.
    const tk = [t0, t0 + ws * 0.8];
    const sk = [0, 0.004 * s];
    for (let k = 1; k <= n; k++) {
      tk.push(hs[k]!);
      // At each heel strike the body is a little nearer the trailing foot than the leading one.
      sk.push(k === 1 ? 0.22 * L : (k - 1.05) * L);
    }
    tk.push(hc, hc + 0.32);
    sk.push(D - 0.05 * L, D);
    const prog = new Mono(tk, sk).ends(0, 0);
    const t1 = hc + 0.32;
    this.segs.push({ t0, t1, type: "walk", x: this.cx, z: this.cz, yaw0: this.cyaw, yaw1: end.yaw, path, prog, speed: v, stepLen: L, upper, opts });
    this.t = t1;
    this.cx = end.x;
    this.cz = end.z;
    this.cyaw = end.yaw;
  }

  /* -------------------------------------------------------------------------------------------------------- */

  private segIndex(t: number): number {
    for (let i = 0; i < this.segs.length; i++) if (t < this.segs[i]!.t1) return i;
    return this.segs.length - 1;
  }

  private foot(side: Side, t: number): { x: number; y: number; z: number; yaw: number; pitch: number; toe: number; swing: number; base: number; clear: number } {
    const arr = this.plants[side];
    let lo = 0;
    let hi = arr.length - 1;
    while (lo < hi) {
      const m = (lo + hi + 1) >> 1;
      if (arr[m]!.tOn <= t) lo = m;
      else hi = m - 1;
    }
    const pl = arr[lo]!;
    if (t <= pl.tOff || lo === arr.length - 1) {
      const f = this.stance(pl, Math.min(t, pl.tOff));
      return { ...f, swing: -1, base: f.y, clear: 0 };
    }
    const nx = arr[lo + 1]!;
    const a = this.stance(pl, pl.tOff);
    const b = this.stance(nx, nx.tOn);
    const u = clamp((t - pl.tOff) / Math.max(1e-6, nx.tOn - pl.tOff), 0, 1);
    const B = this.staff.B;
    const s = B.s;
    const walking = nx.arrive === "walk";
    // The swing foot leaves with some speed (the toe slides off the floor) and settles to land: a Hermite ramp
    // with a starting slope, rather than minimum jerk, which would leave it trailing behind the hip.
    const v0 = walking ? 0.45 : 0;
    const m = walking ? v0 * u + (3 - 2 * v0) * u * u + (v0 - 2) * u * u * u : mj(u);
    const yaw = a.yaw + wrapPi(b.yaw - a.yaw) * mj(clamp(u * 1.15, 0, 1));
    const pitch = walking ? a.pitch * (1 - smooth(0, 0.45, u)) + b.pitch * smooth(0.45, 1, u) : 5 * DEG * Math.sin(Math.PI * u);
    const toe = a.toe * (1 - smooth(0, 0.35, u));
    // The knee folds early in swing and lifts the ankle; it comes down to the floor as the shank swings through.
    const fold = walking ? (0.055 * s * Math.sin(Math.PI * u) * Math.exp(-3 * u)) / 0.334 : 0;
    const base = lerp(a.y, b.y, m) + (walking ? 0.012 : 0.032) * s * Math.pow(Math.sin(Math.PI * u), 0.8) + fold;
    // Keep heel, ball and toe tip clear of the floor in mid swing.
    const clear = 0.012 * s * Math.sin(Math.PI * u);
    const y = this.clearY(base, pitch, toe, clear);
    return { x: lerp(a.x, b.x, m), y, z: lerp(a.z, b.z, m), yaw, pitch, toe, swing: u, base, clear };
  }

  /** The lowest the ankle may be for the heel, ball and toe tip all to clear the floor by `clear`. */
  private clearY(y: number, pitch: number, toe: number, clear: number): number {
    const B = this.staff.B;
    const s = B.s;
    const heel = -B.ankY * Math.cos(pitch) - (B.heelBack + 0.013 * s) * Math.sin(pitch);
    const mtpY = B.J.toe_R!.y - B.ankY;
    const ball = mtpY * Math.cos(pitch) + B.mtpFwd * Math.sin(pitch);
    const sole = -B.ankY * Math.cos(pitch) + B.mtpFwd * Math.sin(pitch);
    const tip = ball + (B.toeFwd + 0.015 * s - B.mtpFwd) * Math.sin(pitch + toe) - (B.J.toe_R!.y - 0.009 * s) * Math.cos(pitch + toe);
    return Math.max(y, clear - Math.min(heel, sole, tip));
  }

  /** The shin's forward lean for a leg reaching from its hip to a body-space ankle (radians, forward positive). */
  private shinLean(p: BodyPose, i: number, fwd: THREE.Vector3): number {
    const B = this.staff.B;
    const side = SIDES[i]!;
    _q2.setFromAxisAngle(Y, p.yaw).multiply(qAxis(X, p.tilt, _q3)).multiply(qAxis(Z, p.roll, _q1));
    const H = _v2.subVectors(B.J[`thigh_${side}`]!, B.J.pelvis!).applyQuaternion(_q2).add(p.pelvis);
    const A = p.feet[i]!.pos;
    const dirV = _v3.subVectors(A, H);
    const d = clamp(dirV.length(), Math.abs(B.thigh - B.shin) + 0.05, (B.thigh + B.shin) * 0.99995);
    const dir = dirV.normalize();
    const ca = clamp((B.thigh * B.thigh + d * d - B.shin * B.shin) / (2 * B.thigh * d), -1, 1);
    const pole = orthoTo(fwd.clone().addScaledVector(X, sgn(side) * 0.12), dir);
    const K = H.clone().addScaledVector(dir, B.thigh * ca).addScaledVector(pole, B.thigh * Math.sqrt(1 - ca * ca));
    const sv = K.sub(H.clone().addScaledVector(dir, d));
    return Math.atan2(sv.dot(fwd), sv.y);
  }

  /** The foot on a plant at time t: heel rocker after a heel strike, rising on the ball before toe-off. */
  private stance(pl: Plant, t: number): { x: number; y: number; z: number; yaw: number; pitch: number; toe: number } {
    const B = this.staff.B;
    const s = B.s;
    let pitch = 0;
    let toe = 0;
    let pv: [number, number] | null = null;
    if (pl.arrive === "walk" && pl.hsDur > 0 && t < pl.tOn + pl.hsDur) {
      const u = clamp((t - pl.tOn) / pl.hsDur, 0, 1);
      pitch = pl.hsAng * (1 - u) * (1 - u * 0.4);
      pv = [-(B.heelBack + 0.013 * s), -B.ankY];
    }
    if (pl.leave === "walk" && pl.toDur > 0 && t > pl.tOff - pl.toDur) {
      const u = clamp((t - (pl.tOff - pl.toDur)) / pl.toDur, 0, 1);
      pitch = -pl.toAng * Math.pow(u, 1.25);
      toe = -pitch;
      pv = [B.mtpFwd, B.J.toe_R!.y - B.ankY];
    }
    let x = pl.x;
    let y = B.ankY;
    let z = pl.z;
    if (pv) {
      // Rotate the ankle about the pivot (heel edge or ball joint) in the foot's own vertical plane.
      const [pz, py] = pv;
      const c = Math.cos(-pitch);
      const sn = Math.sin(-pitch);
      // Foot frame: rotation about x by -pitch takes (y, z) to (y c - z s, y s + z c).
      const ry = py * c - pz * sn;
      const rz = py * sn + pz * c;
      const dzf = pz - rz;
      y = B.ankY + (py - ry);
      x += Math.sin(pl.yaw) * dzf;
      z += Math.cos(pl.yaw) * dzf;
    }
    return { x, y, z, yaw: pl.yaw, pitch, toe };
  }

  /** The share of the body's weight on the left foot at t (unsmoothed). */
  private loadL(t: number): number {
    const load = (side: Side) => {
      const arr = this.plants[side];
      let lo = 0;
      let hi = arr.length - 1;
      while (lo < hi) {
        const m = (lo + hi + 1) >> 1;
        if (arr[m]!.tOn <= t) lo = m;
        else hi = m - 1;
      }
      const p = arr[lo]!;
      if (t > p.tOff || t < p.tOn) return 0;
      return smooth(p.tOn, p.tOn + p.rUp, t) * (1 - smooth(p.tOff - p.rDn, p.tOff, t));
    };
    const l = load("L");
    const r = load("R");
    return l + r > 1e-6 ? l / (l + r) : 0.5;
  }
  /** The left foot's share, smoothed over a short window (the body's mass cannot jump sideways). */
  private support(t: number): number {
    let w = 0;
    for (let k = -3; k <= 3; k++) w += ((4 - Math.abs(k)) / 16) * this.loadL(t + k * 0.045);
    return w;
  }

  /** The upper-body pose a segment holds at time t (standing, working, or carrying / pushing). */
  private upperPose(g: Seg, tm: number, p: BodyPose): BodyPose {
    const st = this.staff;
    if (g.type === "work") return workPose(st, g.work!, tm - g.t0, g.opts ?? {}, p);
    if (g.upper !== "free") return workPose(st, g.upper, tm, g.opts ?? {}, p);
    return standPose(st, tm, p);
  }

  /** Place and pose the person at time t; returns what they are doing. */
  apply(t: number): RoutineState {
    const st = this.staff;
    const B = st.B;
    const s = B.s;
    const D = this.duration;
    const tm = this.loop ? ((t % D) + D) % D : clamp(t, 0, D);
    const gi = this.segIndex(tm);
    const g = this.segs[gi]!;
    const n = this.segs.length;
    const prev = this.segs[gi > 0 ? gi - 1 : this.loop ? n - 1 : 0]!;
    const next = this.segs[gi < n - 1 ? gi + 1 : this.loop ? 0 : n - 1]!;
    const local = tm - g.t0;
    // Root placement.
    let rx = g.x;
    let rz = g.z;
    let ryaw = g.yaw0;
    let along = 0;
    if (g.type === "walk") {
      along = g.prog!.at(tm);
      const a = g.path!.at(along);
      rx = a.x;
      rz = a.z;
      ryaw = a.yaw;
    } else if (g.type === "turn") {
      const u = ease(clamp((tm - g.t0 - 0.1) / Math.max(0.1, g.t1 - g.t0 - 0.3), 0, 1));
      ryaw = g.yaw0 + wrapPi(g.yaw1 - g.yaw0) * u;
      if (g.path) {
        const a = g.path.at(g.path.length * u);
        rx = a.x;
        rz = a.z;
      }
    }
    st.place(rx, rz, ryaw);
    // Upper body: the segment's own, eased from and to standing where the neighbours hold something else.
    const p = this.upperPose(g, tm, this.pose);
    const ramp = g.type === "work" ? 0.9 : 0.6;
    const span = g.t1 - g.t0;
    const rIn = g.upper !== "free" && prev.upper === g.upper ? 1 : smooth(g.t0, g.t0 + Math.min(ramp, span / 2), tm);
    const rOut = g.upper !== "free" && next.upper === g.upper ? 1 : 1 - smooth(g.t1 - Math.min(ramp, span / 2), g.t1, tm);
    const wU = g.type === "work" || g.upper !== "free" ? rIn * rOut : 1;
    if (wU < 1) blendPose(p, standPose(st, tm, this.aux), p, wU);
    // Feet from the plan.
    const reachW: [number, number] = [1, 1];
    const fs: ReturnType<Routine["foot"]>[] = [];
    for (let i = 0; i < 2; i++) {
      const f = this.foot(SIDES[i]!, tm);
      fs.push(f);
      // A foot in mid swing does not hold the hips down; one about to land or just lifting does.
      reachW[i] = f.swing < 0 ? 1 : (1 - smooth(0, 0.3, f.swing)) + smooth(0.65, 1, f.swing);
      const fo = p.feet[i]!;
      st.toBody(_v1.set(f.x, f.y, f.z), fo.pos);
      fo.quat.setFromAxisAngle(Y, f.yaw - ryaw).multiply(qAxis(X, -f.pitch, _q1));
      fo.toe = f.toe;
    }
    // The gait on top: eased in and out at the segment's ends so it meets standing without a jump.
    const moving = g.type === "walk" || g.type === "turn";
    const restY = p.pelvis.y;
    if (moving) {
      const m = smooth(g.t0, g.t0 + 0.45, tm) * (1 - smooth(g.t1 - 0.45, g.t1, tm));
      const mx = smooth(g.t0, g.t0 + 0.15, tm) * (1 - smooth(g.t1 - 0.15, g.t1, tm));
      const wl = this.support(tm);
      const fl = p.feet[0]!.pos;
      const fr = p.feet[1]!.pos;
      const L = g.stepLen ?? 0.5 * s;
      const fd = clamp((fl.z - fr.z) / L, -1.3, 1.3);
      const vk = clamp((g.speed ?? 0.6) / 1.3, 0.2, 1.2);
      const walk = g.type === "walk" ? 1 : 0.5;
      // Over the loaded foot, a third of the way out.
      p.pelvis.x = lerp(p.pelvis.x, 0.32 * (wl * fl.x + (1 - wl) * fr.x), mx);
      p.pelvis.z = lerp(p.pelvis.z, 0.012 * s, m);
      p.pelvis.y = lerp(p.pelvis.y, B.hipY - 0.004 * s, m);
      p.yaw = lerp(p.yaw, -0.07 * fd, m);
      p.roll = lerp(p.roll, 0.06 * (2 * wl - 1) * walk, mx);
      p.tilt = lerp(p.tilt, 0.03 + 0.01 * Math.abs(fd), m);
      if (g.upper === "free") {
        p.flex = lerp(p.flex, 0.035 + 0.035 * vk * walk, m);
        p.twist = lerp(p.twist, 0.11 * fd, m);
        p.side = lerp(p.side, -0.35 * p.roll, m);
        // Arms swing against the legs, more the faster the walk; the forward arm's elbow bends.
        const A = (0.16 + 0.22 * vk) * walk;
        for (let i = 0; i < 2; i++) {
          const a = p.arms[i]!;
          const sw = (i === 0 ? -1 : 1) * A * fd + 0.04;
          a.flex = lerp(a.flex, sw, m);
          a.elbow = lerp(a.elbow, 0.24 + (0.32 * Math.max(0, sw)) / (A + 0.04), m);
          a.abd = lerp(a.abd, 0.07, m);
        }
        // Eyes on the way ahead, or round the turn.
        if (g.type === "walk") {
          const ahead = g.path!.at(along + 2.6);
          st.toBody(_v2.set(ahead.x, B.eyeY - 0.55, ahead.z), _v3);
        } else {
          // The turn still to go (unwrapped: a half turn must not flip the eyes to the other side), led a little.
          const dir = clamp((g.yaw1 - ryaw) * 1.4, -1.4, 1.4);
          _v3.set(Math.sin(dir) * 3, B.eyeY - 0.25, Math.cos(dir) * 3);
        }
        p.gaze.lerp(_v3, m);
      } else p.twist = lerp(p.twist, p.twist + 0.04 * fd, m);
    }
    // Hips as high as wanted and no higher than both legs reach; a carried pan rides down with the shoulders.
    fitHips(st, p, 0.012, reachW);
    // In swing the foot hangs from the shin: pointed as it leaves, level by mid swing, toes up to land.
    if (g.type === "walk")
      for (let i = 0; i < 2; i++) {
        const f = fs[i]!;
        if (f.swing < 0) continue;
        const fwd = V(Math.sin(f.yaw - ryaw), 0, Math.cos(f.yaw - ryaw));
        const lean = this.shinLean(p, i, fwd);
        const rel = lerp(-16, 2, smooth(0.1, 0.65, f.swing)) * DEG;
        const w = smooth(0, 0.22, f.swing) * (1 - smooth(0.72, 1, f.swing));
        const pitch = lerp(f.pitch, rel - lean, w);
        st.toBody(_v1.set(f.x, this.clearY(f.base, pitch, f.toe, f.clear), f.z), p.feet[i]!.pos);
        p.feet[i]!.quat.setFromAxisAngle(Y, f.yaw - ryaw).multiply(qAxis(X, -pitch, _q1));
      }
    const dy = p.pelvis.y - restY;
    if (g.upper === "carry" || (g.type === "work" && g.work === "carry")) {
      p.itemPos.y += dy;
      for (const a of p.arms) if (a.ik > 0) a.pos.y += dy;
    }
    st.pose(p);
    return { activity: g.type, work: g.work, local, holding: p.itemOn };
  }
}

export interface RouteOptions {
  /** Walking speed, m/s (default 0.8). */
  speed?: number;
  /** Seconds spent at each waypoint (one value, or one per waypoint). Default 4. */
  dwell?: number | readonly number[];
  /** What is done at each waypoint (one kind, or one per waypoint). Default "stand". */
  work?: WorkKind | readonly WorkKind[];
  /** Which way to face at each waypoint (yaw), or null to face the way they arrived. */
  facing?: readonly (number | null | undefined)[];
  /** Work options per waypoint (world targets). */
  workOpts?: readonly (WorkOpts | undefined)[];
  /** Upper body while walking each leg (leg i runs from waypoint i to i + 1). */
  upper?: readonly (Upper | undefined)[];
}

/**
 * A looping route through floor waypoints (the scene's x, z): at each waypoint the person turns to face the
 * station and works there, then turns and walks to the next, round and back to the first.
 */
export function routeRoutine(staff: Staff, waypoints: readonly (readonly [number, number])[], o: RouteOptions = {}): Routine {
  const n = waypoints.length;
  const at = <T>(v: T | readonly T[] | undefined, i: number, d: T): T => (Array.isArray(v) ? ((v as readonly T[])[i] ?? d) : ((v as T | undefined) ?? d));
  const heading = (a: readonly [number, number], b: readonly [number, number]) => Math.atan2(b[0] - a[0], b[1] - a[1]);
  if (n < 2) {
    const w = waypoints[0] ?? [0, 0];
    return new Routine(staff, { x: w[0], z: w[1], yaw: o.facing?.[0] ?? 0 }, [{ kind: "work", work: at(o.work, 0, "stand"), duration: at(o.dwell, 0, 4), opts: o.workOpts?.[0] }]);
  }
  const face = (i: number) => o.facing?.[i] ?? heading(waypoints[(i - 1 + n) % n]!, waypoints[i]!);
  const steps: RoutineStep[] = [];
  for (let i = 0; i < n; i++) {
    const j = (i + 1) % n;
    steps.push({ kind: "work", work: at(o.work, i, "stand"), duration: at(o.dwell, i, 4), opts: o.workOpts?.[i] });
    const upper = o.upper?.[i] ?? "free";
    steps.push({ kind: "walk", path: [waypoints[j]!], speed: o.speed ?? 0.8, upper });
    steps.push({ kind: "turn", yaw: face(j), upper });
  }
  return new Routine(staff, { x: waypoints[0]![0], z: waypoints[0]![1], yaw: face(0) }, steps, true);
}

/** Build a person: chef whites, apron, trousers and shoes on a body built to the cohort's row. */
export function makeStaff(cohort: Cohort, detail: Detail = 1, opts: StaffOptions = {}): Staff {
  return new Staff(cohort, detail, opts);
}
