/**
 * CHAPTER 5 · THE BRIGADE & THE FLOOR.
 *
 * A slice of a cookline on a quarry-tile floor. Two 72 in prep tables with backsplashes form the back line, with a low
 * wire shelf in front of each cook holding a full 4 in hotel pan. Across the aisle are two 36 in six-burner ranges and a 20-pan speed rack. Three cooks
 * built to the anthropometric cohorts (P05 female, P50, P95 male) lift their pans to the carry, hold them, and set
 * them back, in step.
 *
 * Everything shown is computed:
 *   - The aisle slider moves the range line. The strip on the floor is the aisle; its bright band is the floor left
 *     free behind the cooks, measured from their skinned bodies' rearmost points over the whole lift. Where the aisle
 *     is too narrow, a red block marks the depth by which a body would pass into the range or rack behind it.
 *   - The load slider is the loaded pan's mass. The food's depth is that mass less the pan's steel, at the food's
 *     density, in the pan's tapered volume.
 *   - Each card's L5/S1 moment is estimateLumbarTorqueNm, fed the trunk flexion and the load's horizontal reach that
 *     the body shows this frame (people.ts liftMetrics). The pan's weight counts only once it is off the shelf.
 *   - The reach shells are each cohort's thumb-tip forward reach (anthropometrics.tsv), centred on the shoulder joint:
 *     the reach is measured from the back of the shoulder blades, so the shell's radius is that less the joint's
 *     depth in front of the back.
 * The motion is a kinematic illustration: the poses are solved for reach and joint range, not from forces.
 */

import * as THREE from "three";
import { estimateLumbarTorqueNm } from "../sim/physicsKernels";
import { bindCanvasResize, createScrollSafeOrbit, createStudioRenderer, observeInView, shaderGate } from "./kit";
import { buildFixture, type Model } from "./models/all";
import { wireShelving } from "./models/cold";
import { palette } from "./models/materials";
import { Shop } from "./models/shop";
import { studioEnvironment, studioLights } from "./models/stage";
import { HOTEL_PAN, makeStaff, type Cohort, type Staff, type WorkOpts } from "./models/people";

const IN = 0.0254;
/** The lift's cycle in people.ts: stand, stoop, take, rise to the carry, hold, lower, release, stand. */
const LIFT_CYCLE = 8.6;

interface CohortSpec {
  readonly key: Cohort;
  readonly name: string;
  readonly accent: number;
  readonly apron: string;
  readonly seed: number;
  /** Food in the pan (smallwares: 0 braise, 3 soup) and its density, kg/m³. */
  readonly food: 0 | 1 | 2 | 3;
  readonly density: number;
  readonly out: string;
}

const COHORTS: readonly CohortSpec[] = [
  { key: "p05", name: "short", accent: 0x38bdf8, apron: "#26344d", seed: 4, food: 0, density: 1050, out: "[data-br-p05]" },
  { key: "p50", name: "average", accent: 0x34d399, apron: "#4a5a6e", seed: 2, food: 3, density: 1020, out: "[data-br-p50]" },
  { key: "p95", name: "tall", accent: 0xc084fc, apron: "#3b3d43", seed: 7, food: 0, density: 1050, out: "[data-br-p95]" },
];

/* Layout: x along the line, y up, z across the aisle. The prep tables' front edge is z = 0 and they face +z; the
   range line faces -z with its front at z = aisle. The cooks stand between, facing the tables. */
const SPACING = 1.2;
const TABLE_EXT = { w: 1.829, d: 0.762, h: 1.016 };
const RANGE_EXT = { w: 0.914, d: 0.864, h: 1.473 };
const RACK_EXT = { w: 0.521, d: 0.66, h: 1.753 };
const PAN_EXT = { w: 0.53, d: 0.325, h: 0.1 };
/** 36 x 14 in chrome wire shelving on 14 in posts, two shelves: the top one the pan's. */
const SHELF = { L: 36 * IN, D: 14 * IN, H: 14 * IN, tiers: [0.17, 0.33] };
/** The shelf's back edge, tucked under the table top's 60 mm front overhang, clear of its legs. */
const SHELF_BACK = -0.03;
/** The pan's centre on the shelf, a centimetre in from its front edge. */
const PAN_Z = SHELF_BACK + SHELF.D - PAN_EXT.d / 2 - 0.012;
const PAN_Y = SHELF.tiers[1]! + PAN_EXT.h;
/** The pan's centre of mass in the item frame (origin at the rim's centre): about mid-depth. */
const PAN_COM = new THREE.Vector3(0, -0.045, 0);
/** What stands across the aisle, by its centre along the line. */
const LINE_B = [
  { kind: "range", x: -0.914, name: "the range" },
  { kind: "range", x: 0, name: "the range" },
  { kind: "rack", x: 1.2, name: "the speed rack" },
] as const;

/* ------------------------------------------------------------------------------------------------------------ */
/* The pan: its steel and its food                                                                                */
/* ------------------------------------------------------------------------------------------------------------ */

const PAN_TAN = Math.tan(HOTEL_PAN.taper);
/** Wall half-sizes at the rim and at the base. */
const PAN_A1 = HOTEL_PAN.length / 2 - HOTEL_PAN.flange;
const PAN_B1 = HOTEL_PAN.width / 2 - HOTEL_PAN.flange;
const PAN_A0 = PAN_A1 - HOTEL_PAN.depth * PAN_TAN;
const PAN_B0 = PAN_B1 - HOTEL_PAN.depth * PAN_TAN;

/**
 * The empty pan's mass: 22 gauge (0.79 mm) type 304 at 7930 kg/m³ over its drawn area (the base, the tapered walls,
 * the flange and its rolled bead).
 */
const PAN_STEEL_KG = (() => {
  const base = 4 * PAN_A0 * PAN_B0;
  const walls = 2 * (PAN_A0 + PAN_A1 + PAN_B0 + PAN_B1) * (HOTEL_PAN.depth / Math.cos(HOTEL_PAN.taper));
  const flange = HOTEL_PAN.length * HOTEL_PAN.width - 4 * PAN_A1 * PAN_B1;
  const bead = 2 * (HOTEL_PAN.length + HOTEL_PAN.width) * Math.PI * 0.0022;
  return (base + walls + flange + bead) * 0.000794 * 7930;
})();

/** Food volume (m³) up to height y above the pan's floor: the walls lean out by the taper. */
function panVolume(y: number): number {
  return 4 * (PAN_A0 * PAN_B0 * y + ((PAN_A0 + PAN_B0) * PAN_TAN * y * y) / 2 + (PAN_TAN * PAN_TAN * y * y * y) / 3);
}

/** The pan's `fill` (smallwares: 0 at its floor, 1 at 8 mm under the rim) for a loaded mass, and whether it fits. */
function panFill(loadKg: number, density: number): { fill: number; overfull: boolean } {
  const y0 = 0.0008;
  const yMax = PAN_EXT.h - 0.008;
  const v = Math.max(0, loadKg - PAN_STEEL_KG) / density;
  if (v <= 0) return { fill: 0, overfull: false };
  if (v >= panVolume(yMax)) return { fill: 1, overfull: true };
  let lo = 0;
  let hi = yMax;
  for (let k = 0; k < 40; k++) {
    const m = (lo + hi) / 2;
    if (panVolume(m) < v) lo = m;
    else hi = m;
  }
  return { fill: Math.max(0, (lo - y0) / (yMax - y0)), overfull: false };
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Scene pieces                                                                                                   */
/* ------------------------------------------------------------------------------------------------------------ */

function backdrop(): THREE.CanvasTexture {
  const c = document.createElement("canvas");
  c.width = 4;
  c.height = 256;
  const x = c.getContext("2d")!;
  const g = x.createLinearGradient(0, 0, 0, 256);
  g.addColorStop(0, "#0a0807");
  g.addColorStop(0.5, "#17110d");
  g.addColorStop(0.78, "#231811");
  g.addColorStop(1, "#120d0a");
  x.fillStyle = g;
  x.fillRect(0, 0, 4, 256);
  const t = new THREE.CanvasTexture(c);
  t.colorSpace = THREE.SRGBColorSpace;
  return t;
}

/** A plane of quarry tile `w` by `h`, the texture at its true 152 mm module. */
function tilePlane(w: number, h: number): THREE.Mesh {
  const P = palette();
  const g = new THREE.PlaneGeometry(w, h);
  const uv = g.attributes.uv as THREE.BufferAttribute;
  for (let i = 0; i < uv.count; i++) uv.setXY(i, (uv.getX(i) * w) / 0.608, (uv.getY(i) * h) / 0.608);
  const m = new THREE.Mesh(g, P.tile);
  m.receiveShadow = true;
  return m;
}

/** The low wire shelf, from the shelving family's construction. */
function lowShelf(): THREE.Group {
  const s = new Shop(2, "wire-shelf-36x14x14");
  wireShelving(s, SHELF.L, SHELF.D, SHELF.H, { tiers: [...SHELF.tiers] });
  return s.build(true).group;
}

/** A fixture placed facing +z (`flip` false) or -z, its front-face centre at (x, z). */
function place(m: Model, x: number, z: number, flip: boolean): THREE.Group {
  const g = m.group;
  if (flip) {
    g.rotation.y = Math.PI;
    g.position.set(x + m.ext.w / 2, 0, z);
  } else g.position.set(x - m.ext.w / 2, 0, z);
  return g;
}

/**
 * A translucent spherical shell of radius 1 facing +z: ±78° round the vertical, from 28° under the zenith to 50°
 * below the horizontal, with a sparse graticule and its rim drawn on it.
 */
function reachShell(color: number): { group: THREE.Group; mats: THREE.Material[] } {
  const group = new THREE.Group();
  const phi0 = Math.PI / 2 - (78 * Math.PI) / 180;
  const phiL = (156 * Math.PI) / 180;
  const th0 = (28 * Math.PI) / 180;
  const thL = ((90 - 28 + 50) * Math.PI) / 180;
  // three's sphere: x = -cos(phi) sin(theta), y = cos(theta), z = sin(phi) sin(theta).
  const at = (phi: number, th: number) => new THREE.Vector3(-Math.cos(phi) * Math.sin(th), Math.cos(th), Math.sin(phi) * Math.sin(th));
  const skin = new THREE.MeshPhysicalMaterial({
    color,
    emissive: color,
    emissiveIntensity: 0.22,
    roughness: 0.3,
    metalness: 0,
    transparent: true,
    opacity: 0.1,
    side: THREE.DoubleSide,
    depthWrite: false,
  });
  const shell = new THREE.Mesh(new THREE.SphereGeometry(1, 56, 28, phi0, phiL, th0, thL), skin);
  shell.renderOrder = 2;
  group.add(shell);
  const pts: THREE.Vector3[] = [];
  const line = (f: (u: number) => THREE.Vector3, n: number) => {
    for (let i = 0; i < n; i++) pts.push(f(i / n), f((i + 1) / n));
  };
  for (const thd of [45, 90, 120])
    line((u) => at(phi0 + phiL * u, (thd * Math.PI) / 180), 48);
  for (let k = 0; k <= 6; k++) line((u) => at(phi0 + (phiL * k) / 6, th0 + thL * u), 32);
  line((u) => at(phi0 + phiL * u, th0), 48);
  line((u) => at(phi0 + phiL * u, th0 + thL), 48);
  const lines = new THREE.LineSegments(
    new THREE.BufferGeometry().setFromPoints(pts),
    new THREE.LineBasicMaterial({ color, transparent: true, opacity: 0.42, depthWrite: false }),
  );
  lines.renderOrder = 3;
  group.add(lines);
  return { group, mats: [skin, lines.material as THREE.Material] };
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Cooks                                                                                                          */
/* ------------------------------------------------------------------------------------------------------------ */

interface Cook {
  spec: CohortSpec;
  st: Staff;
  opts: WorkOpts;
  pan: Model;
  shell: THREE.Group;
  /** Rest-pose sample points on the back of the trunk (with the cloth over it), by the bone they ride on. */
  back: { bone: string; off: THREE.Vector3 }[];
  skinBack: { mesh: THREE.SkinnedMesh; idx: number[] }[];
  /** The body's rearmost reach into the aisle over the whole lift (world z), and that point's height. */
  rearZ: number;
  rearY: number;
  /** Where the body would pass into what stands behind it, when the aisle is too narrow. */
  clash: THREE.Mesh;
  shellR: number;
  /** The table cell's three parts (back load and peak, lean, reach), each kept on one line, or null without a cell. */
  out: [HTMLElement, HTMLElement, HTMLElement] | null;
  fill: number;
  peakNm: number;
}

/** Split a readout cell into three spans joined by " · ", so a narrow screen wraps between the parts, never inside one. */
function readoutParts(cell: HTMLElement | null): [HTMLElement, HTMLElement, HTMLElement] | null {
  if (!cell) return null;
  const parts: [HTMLElement, HTMLElement, HTMLElement] = [document.createElement("span"), document.createElement("span"), document.createElement("span")];
  cell.replaceChildren(parts[0], " · ", parts[1], " · ", parts[2]);
  return parts;
}

/** Points over the back of the trunk, from the gluteal fold to the shoulders, each tied to its spine bone. */
function backPoints(st: Staff): { bone: string; off: THREE.Vector3 }[] {
  const B = st.B;
  const J = B.J;
  const out: { bone: string; off: THREE.Vector3 }[] = [];
  const y0 = J.pelvis!.y - 0.12 * B.s;
  const y1 = J.neck!.y - 0.03 * B.s;
  for (let i = 0; i <= 8; i++) {
    const y = y0 + ((y1 - y0) * i) / 8;
    const bone = y < J.spine_01!.y ? "pelvis" : y < J.spine_02!.y ? "spine_01" : y < J.spine_03!.y ? "spine_02" : "spine_03";
    for (const da of [-0.9, -0.45, 0, 0.45, 0.9]) {
      const th = Math.PI + da;
      const p = st.torso.at(y, th).addScaledVector(st.torso.normal(y, th), 0.014);
      out.push({ bone, off: p.sub(J[bone]!) });
    }
  }
  return out;
}

const _p = new THREE.Vector3();
const _q = new THREE.Quaternion();
const _w = new THREE.Vector3();

/** The skin and cloth vertices over the back of the body from the knees up (bind pose, body space), every other one. */
function backVertices(st: Staff): { mesh: THREE.SkinnedMesh; idx: number[] }[] {
  const J = st.B.J;
  const y0 = J.shin_L!.y;
  const y1 = J.neck!.y;
  const z0 = J.pelvis!.z;
  return st.meshes.map((mesh) => {
    const pos = mesh.geometry.attributes.position as THREE.BufferAttribute;
    const idx: number[] = [];
    for (let i = 0; i < pos.count; i += 2) {
      const y = pos.getY(i);
      if (y > y0 && y < y1 && pos.getZ(i) < z0) idx.push(i);
    }
    return { mesh, idx };
  });
}

/** How far into the aisle (world z) the body's surface reaches in the current pose. */
function rearSkin(c: Cook): { z: number; y: number } {
  c.st.group.updateMatrixWorld(true);
  let z = -Infinity;
  let y = 0;
  for (const { mesh, idx } of c.skinBack) {
    for (const i of idx) {
      mesh.getVertexPosition(i, _p).applyMatrix4(mesh.matrixWorld);
      if (_p.z > z) {
        z = _p.z;
        y = _p.y;
      }
    }
  }
  return { z, y };
}

/** How far into the aisle (world z) the back of the trunk reaches in the current pose: quick, from the joints. */
function rearNow(c: Cook): number {
  let z = -Infinity;
  for (const b of c.back) {
    c.st.jointQuat(b.bone, _q);
    c.st.jointPos(b.bone, _p).add(_w.copy(b.off).applyQuaternion(_q));
    z = Math.max(z, c.st.toWorld(_p, _w).z);
  }
  return z;
}

/* ------------------------------------------------------------------------------------------------------------ */

export function mountBrigadeStage(root: HTMLElement): void {
  const mountT0 = performance.now();
  const stageEl = root.querySelector<HTMLElement>("[data-brigade-stage]")!;
  const canvas = root.querySelector<HTMLCanvasElement>("[data-brigade-canvas]")!;
  const maybeRenderer = createStudioRenderer(canvas, false);
  if (!maybeRenderer) return;
  const renderer = maybeRenderer;
  renderer.toneMappingExposure = 1.0;

  const scene = new THREE.Scene();
  scene.background = backdrop();
  scene.fog = new THREE.Fog(0x17110d, 6.5, 12.5);
  scene.environment = studioEnvironment(renderer);
  scene.environmentIntensity = 0.62;
  const lights = studioLights(scene, new THREE.Vector3(0, 0.7, 0.55), 2.7, 2048);
  lights.key.position.set(-1.9, 4.6, 3.6);
  lights.key.intensity = 2.1;
  lights.key.shadow.normalBias = 0.02;
  lights.rim.position.set(2.8, 2.2, -2.6);
  lights.rim.intensity = 0.9;
  lights.hemi.intensity = 0.5;

  /* The room: a quarry-tile floor from the line's back edge out across the aisle, falling away into the dark. */
  const wallZ = -TABLE_EXT.d - 0.012;
  // The floor runs from the wall out across the aisle and on past the range line, into the dark.
  const floor = tilePlane(9.6, 6.2);
  floor.rotation.x = -Math.PI / 2;
  floor.position.set(0, 0, wallZ + 3.1);
  scene.add(floor);

  /* Line A: the prep tables, a cutting board, a scale and a half pan on them; the low shelves in front. */
  for (const sx of [-1, 1]) {
    const t = buildFixture("NSF-PREPTABLE-72", TABLE_EXT, 1);
    scene.add(place(t, (sx * TABLE_EXT.w) / 2, 0, false));
  }
  const prop = (id: string, ext: { w: number; d: number; h: number }, x: number, y: number, z: number, yaw = 0) => {
    const m = buildFixture(id, ext, 1);
    const holder = new THREE.Group();
    m.group.position.set(-ext.w / 2, 0, ext.d / 2);
    holder.add(m.group);
    holder.position.set(x, y, z);
    holder.rotation.y = yaw;
    scene.add(holder);
    return m;
  };
  const TOP = 0.914;
  prop("NSF-CUTTING-BOARD-HDPE-24X18", { w: 0.61, d: 0.457, h: 0.019 }, -1.25, TOP, -0.36, 0.04);
  prop("NSF-CHEF-KNIFE-10IN", { w: 0.052, d: 0.385, h: 0.024 }, -1.08, TOP + 0.019, -0.36, 1.2);
  prop("NSF-PORTION-SCALE-DIGI", { w: 0.33, d: 0.356, h: 0.114 }, 0.34, TOP, -0.42, -0.1);
  const half = prop("NSF-HOTEL-PAN-HALF-4IN", { w: 0.325, d: 0.265, h: 0.1 }, -0.32, TOP, -0.4, 0.05);
  half.hooks.set?.("food", 2);
  half.hooks.set?.("fill", 0.7);
  const half2 = prop("NSF-HOTEL-PAN-HALF-4IN", { w: 0.325, d: 0.265, h: 0.1 }, 1.25, TOP, -0.42, -0.03);
  half2.hooks.set?.("food", 1);
  half2.hooks.set?.("fill", 0.55);
  const shelfTpl = lowShelf();
  for (let i = 0; i < COHORTS.length; i++) {
    const s = i === 0 ? shelfTpl : shelfTpl.clone();
    s.position.set((i - 1) * SPACING - SHELF.L / 2, 0, SHELF_BACK + SHELF.D);
    scene.add(s);
    // Spare 2½ in pans nested on the lower shelf.
    const spare = prop("NSF-HOTEL-PAN-FULL-2IN", { w: 0.53, d: 0.325, h: 0.065 }, (i - 1) * SPACING + 0.04, SHELF.tiers[0]!, PAN_Z + 0.005, 0.02);
    spare.hooks.set?.("fill", 0);
  }

  /* Line B: two six-burner ranges (a stockpot on one) and the speed rack, moved by the aisle slider. */
  const lineB = new THREE.Group();
  scene.add(lineB);
  const ranges: Model[] = [];
  for (const b of LINE_B) {
    if (b.kind === "range") {
      const r = buildFixture("NSF-RANGE-GAS-6B", RANGE_EXT, 1);
      lineB.add(place(r, b.x, 0, true));
      ranges.push(r);
    } else {
      const r = buildFixture("NSF-SPEEDRACK-20PAN", RACK_EXT, 1);
      r.hooks.set?.("pans", 9);
      lineB.add(place(r, b.x, 0, true));
    }
  }
  {
    // The stockpot on the second range's back-left burner (burner 3), simmering.
    const r = ranges[1]!;
    const bx = LINE_B[1].x + RANGE_EXT.w / 2 - 0.18;
    const pot = buildFixture("NSF-STOCKPOT-20QT", { w: 0.31, d: 0.38, h: 0.295 }, 1);
    pot.group.position.set(bx - 0.31 / 2, 0.914, 0.58 + 0.38 / 2);
    pot.group.rotation.y = 0;
    lineB.add(pot.group);
    r.hooks.flames?.setFlow(3, 0.35);
    ranges[0]!.hooks.flames?.setFlow(1, 0.2);
  }

  /* The aisle strip: the whole aisle faint, the floor free behind the cooks bright, the state's colour. */
  const stripMat = new THREE.MeshBasicMaterial({ color: 0x34d399, transparent: true, opacity: 0.12, depthWrite: false });
  const freeMat = new THREE.MeshBasicMaterial({ color: 0x34d399, transparent: true, opacity: 0.3, depthWrite: false });
  const edgeMat = new THREE.LineBasicMaterial({ color: 0x34d399, transparent: true, opacity: 0.85 });
  const LINE_LEN = TABLE_EXT.w * 2;
  const unitPlane = new THREE.PlaneGeometry(1, 1).rotateX(-Math.PI / 2).translate(0, 0, 0.5);
  const strip = new THREE.Mesh(unitPlane, stripMat);
  strip.scale.x = LINE_LEN;
  strip.position.y = 0.003;
  strip.renderOrder = 1;
  const free = new THREE.Mesh(unitPlane, freeMat);
  free.scale.x = LINE_LEN;
  free.position.y = 0.004;
  free.renderOrder = 1;
  const edgeGeo = new THREE.BufferGeometry().setFromPoints([new THREE.Vector3(-LINE_LEN / 2, 0, 0), new THREE.Vector3(LINE_LEN / 2, 0, 0)]);
  const edgeA = new THREE.Line(edgeGeo, edgeMat);
  const edgeB = new THREE.Line(edgeGeo, edgeMat);
  const edgeR = new THREE.Line(edgeGeo, edgeMat);
  for (const e of [edgeA, edgeB, edgeR]) e.position.y = 0.005;
  scene.add(strip, free, edgeA, edgeB, edgeR);

  /* The cooks. A red block marks where a body would pass into the range line in too narrow an aisle. */
  const clashGeo = new THREE.BoxGeometry(1, 1, 1);
  const clashMat = new THREE.MeshBasicMaterial({ color: 0xfb7185, transparent: true, opacity: 0.38, depthWrite: false });
  const cooks: Cook[] = COHORTS.map((spec, i) => {
    const st = makeStaff(spec.key, 2, { seed: spec.seed, apron: spec.apron, sleeves: i === 1 ? "rolled" : "cuffed" });
    const x = (i - 1) * SPACING;
    st.place(x, PAN_Z + 0.43 * st.B.s, Math.PI);
    scene.add(st.group);
    const pan = buildFixture("NSF-HOTEL-PAN-FULL-4IN", PAN_EXT, 2);
    pan.hooks.set?.("food", spec.food);
    // Fixture frame (x 0..w, z -d..0, y 0..h) to the item socket's (origin at the rim's centre).
    pan.group.position.set(-PAN_EXT.w / 2, -PAN_EXT.h, PAN_EXT.d / 2);
    st.item.add(pan.group);
    const shell = reachShell(spec.accent);
    scene.add(shell.group);
    // Thumb-tip reach is measured from the back of the shoulder blades; the shell centres on the shoulder joint.
    const J = st.B.J;
    const sh = J.upper_arm_R!;
    const backZ = st.torso.at(sh.y, Math.PI).z;
    shell.group.scale.setScalar(Math.max(0.3, st.dims.forward_reach_m - (sh.z - backZ)));
    const c: Cook = {
      spec,
      st,
      opts: { at: new THREE.Vector3(x, PAN_Y, PAN_Z), surfaceY: SHELF.tiers[1] },
      pan,
      shell: shell.group,
      back: backPoints(st),
      skinBack: backVertices(st),
      rearZ: 0,
      rearY: 0,
      clash: new THREE.Mesh(clashGeo, clashMat),
      shellR: shell.group.scale.x,
      out: readoutParts(root.querySelector<HTMLElement>(spec.out)),
      fill: -1,
      /** The largest back load seen since the pan weight was last set: the bottom of the lift. */
      peakNm: 0,
    };
    measureRear(c);
    c.clash.visible = false;
    c.clash.renderOrder = 4;
    scene.add(c.clash);
    return c;
  });

  /** The rearmost reach over the lift, from the poses themselves: the trunk's joints find the deepest stoops on the
   * way down and back, the skinned surface says how far back the body reaches there. */
  function measureRear(c: Cook): void {
    const st = c.st;
    const N = 43;
    const quick: number[] = [];
    for (let k = 0; k < N; k++) {
      st.work("lift", (k / N) * LIFT_CYCLE, c.opts);
      quick.push(rearNow(c));
    }
    const peaks = quick.map((z, k) => ({ z, k })).filter(({ z, k }) => z >= quick[(k + N - 1) % N]! && z >= quick[(k + 1) % N]!);
    peaks.sort((a, b) => b.z - a.z);
    c.rearZ = -Infinity;
    for (const { k } of peaks.slice(0, 2))
      for (const dk of [-1, -0.5, 0, 0.5, 1]) {
        st.work("lift", ((k + dk) / N) * LIFT_CYCLE, c.opts);
        const r = rearSkin(c);
        if (r.z > c.rearZ) {
          c.rearZ = r.z;
          c.rearY = r.y;
        }
      }
  }

  /* Camera: down the aisle from the P05 end, raised, a little to the tables' side so the hands and pans show. */
  const camera = new THREE.PerspectiveCamera(30, 16 / 10, 0.05, 40);
  const TARGET = new THREE.Vector3(0.1, 0.7, 0.62);
  const VIEW = new THREE.Vector3(-3.4, 2.25, -0.62).normalize();
  const controls = createScrollSafeOrbit(camera, canvas, [TARGET.x, TARGET.y, TARGET.z], {
    minDistance: 0.55,
    maxDistance: 10,
    maxPolarAngle: Math.PI * 0.47,
  });
  let userMoved = false;
  controls.addEventListener("start", () => {
    userMoved = true;
  });
  /** The boxes the view must hold: each cook standing, the pan's shelf in front, the floor behind at full stoop. */
  const keep: THREE.Vector3[] = [];
  const keepBoxes = () => {
    keep.length = 0;
    for (const c of cooks) {
    const x = c.st.placement.x;
    for (const dx of [-0.42, 0.42]) for (const y of [0, c.st.dims.standing_height_m + 0.04]) for (const z of [SHELF_BACK, c.rearZ + 0.06]) keep.push(new THREE.Vector3(x + dx, y, z));
    }
  };
  keepBoxes();
  const ndc = new THREE.Vector3();
  /** Fit the three cooks to the stage, whatever its shape: the nearest distance that holds every box corner. */
  const frame = () => {
    if (userMoved) return;
    let dist = 2;
    for (; dist < 14; dist += 0.05) {
      camera.position.copy(TARGET).addScaledVector(VIEW, dist);
      camera.lookAt(TARGET);
      camera.updateMatrixWorld(true);
      if (keep.every((p) => (ndc.copy(p).project(camera), Math.abs(ndc.x) < 0.94 && Math.abs(ndc.y) < 0.92))) break;
    }
    controls.target.copy(TARGET);
    controls.update();
  };
  const gate = shaderGate(renderer, camera);
  bindCanvasResize(stageEl, renderer, camera, gate);
  frame();
  new ResizeObserver(() => frame()).observe(stageEl);

  /* Controls and readouts. */
  const aisleInput = root.querySelector<HTMLInputElement>("[data-br-aisle]");
  const aisleVal = root.querySelector<HTMLElement>("[data-br-aisle-val]");
  const loadInput = root.querySelector<HTMLInputElement>("[data-br-load]");
  const loadVal = root.querySelector<HTMLElement>("[data-br-load-val]");
  const reachToggle = root.querySelector<HTMLButtonElement>("[data-br-reach]");
  const presetBtns = Array.from(root.querySelectorAll<HTMLButtonElement>("[data-br-preset]"));
  const statusBadge = root.querySelector<HTMLElement>("[data-br-status]");

  let aisleIn = Number(aisleInput?.value ?? 44);
  let loadKg = Number(loadInput?.value ?? 6.5);
  let showReach = reachToggle ? reachToggle.getAttribute("aria-pressed") !== "false" : true;
  // A P95 slipping past sideways needs his chest depth and a hand's breadth either side.
  const passNeed = cooks[2]!.st.dims.chest_depth_m + 2 * 0.05;

  const STATE_COLOR = { ok: 0x34d399, warn: 0xf59e0b, danger: 0xfb7185 } as const;

  function syncAisle(): void {
    const aisleM = aisleIn * IN;
    lineB.position.z = aisleM;
    strip.scale.z = aisleM;
    edgeA.position.z = 0;
    edgeB.position.z = aisleM;
    // The cook reaching furthest back, against what stands behind them.
    let worst = cooks[0]!;
    for (const c of cooks) if (aisleM - c.rearZ < aisleM - worst.rearZ) worst = c;
    const gap = aisleM - worst.rearZ;
    free.position.z = Math.min(aisleM, worst.rearZ);
    free.scale.z = Math.max(0.001, gap);
    free.visible = gap > 0.005;
    edgeR.position.z = Math.min(aisleM, worst.rearZ);
    for (const c of cooks) {
      const over = c.rearZ - aisleM;
      c.clash.visible = over > 0.002;
      if (c.clash.visible) {
        c.clash.scale.set(0.44 * c.st.B.s, 0.42 * c.st.B.s, over);
        c.clash.position.set(c.st.placement.x, c.rearY, aisleM + over / 2);
      }
    }
    const behind = LINE_B.reduce((a, b) => (Math.abs(b.x - worst.st.placement.x) < Math.abs(a.x - worst.st.placement.x) ? b : a));
    const head = `${aisleIn.toFixed(0)}" (${aisleM.toFixed(2)} m)`;
    let state: keyof typeof STATE_COLOR = "ok";
    let msg: string;
    if (gap < 0) {
      state = "danger";
      const overMm = -gap * 1000;
      const by = overMm < 0.5 ? "just touching" : overMm < 10 ? `${overMm.toFixed(0)} mm too close` : `${(overMm / 10).toFixed(0)} cm too close`;
      msg = `Contact at ${head}: at the bottom of the lift the ${worst.spec.name} cook's back meets ${behind.name}, ${by}.`;
    } else if (gap < passNeed) {
      state = "warn";
      const free = gap < 0.01 ? `${(gap * 1000).toFixed(0)} mm` : `${(gap * 100).toFixed(0)} cm`;
      msg = `No passing at ${head}: ${free} free behind the stooping ${worst.spec.name} cook, and a tall cook passing sideways needs ${(passNeed * 100).toFixed(0)} cm.`;
    } else msg = `Clear to pass at ${head}: ${(gap * 100).toFixed(0)} cm free behind the cooks at full stoop.`;
    const col = STATE_COLOR[state];
    stripMat.color.setHex(col);
    freeMat.color.setHex(col);
    edgeMat.color.setHex(col);
    if (statusBadge) {
      statusBadge.dataset.state = state;
      statusBadge.textContent = msg;
    }
    if (aisleVal) aisleVal.textContent = head;
  }

  function syncLoad(): void {
    if (loadVal) loadVal.textContent = `${loadKg.toFixed(1)} kg`;
    for (const c of cooks) {
      c.peakNm = 0;
      const f = panFill(loadKg, c.spec.density).fill;
      if (Math.abs(f - c.fill) > 1e-4) {
        c.fill = f;
        c.pan.hooks.set?.("fill", f);
      }
    }
  }

  function syncReach(): void {
    for (const c of cooks) c.shell.visible = showReach;
    reachToggle?.setAttribute("aria-pressed", String(showReach));
  }

  aisleInput?.addEventListener("input", () => {
    aisleIn = Number(aisleInput.value);
    syncAisle();
  });
  loadInput?.addEventListener("input", () => {
    loadKg = Number(loadInput.value);
    syncLoad();
    readouts();
  });
  reachToggle?.addEventListener("click", () => {
    showReach = !showReach;
    syncReach();
  });
  for (const b of presetBtns) {
    b.addEventListener("click", () => {
      const v = Number(b.dataset.brPreset);
      if (!Number.isFinite(v)) return;
      aisleIn = v;
      if (aisleInput) {
        aisleInput.value = String(v);
        aisleInput.dispatchEvent(new Event("input", { bubbles: true }));
      }
      syncAisle();
    });
  }

  /** L5/S1 moment from what each body shows now: its trunk flexion, and the pan's reach once the pan is lifted. */
  function readouts(): void {
    for (const c of cooks) {
      const m = c.st.liftMetrics(PAN_COM);
      const cur = c.st.current;
      const off = cur.itemOn ? THREE.MathUtils.smoothstep(cur.itemPos.y, PAN_Y + 0.001, PAN_Y + 0.012) : 0;
      const reach = off > 0 ? m.loadReachM : m.reachM;
      const tau = estimateLumbarTorqueNm(c.st.massKg, loadKg * off, reach, m.torsoFlexDeg);
      c.peakNm = Math.max(c.peakNm, tau);
      if (c.out) {
        c.out[0].textContent = `${tau.toFixed(0)} N·m (peak ${c.peakNm.toFixed(0)})`;
        c.out[1].textContent = `${Math.max(0, m.torsoFlexDeg).toFixed(0)}°`;
        c.out[2].textContent = `reach ${reach.toFixed(2)} m`;
      }
    }
  }

  /* Animation. */
  const mountMs = performance.now() - mountT0;
  const reduced = window.matchMedia?.("(prefers-reduced-motion: reduce)").matches ?? false;
  let t = reduced ? 3.6 : 0;
  let frozen: number | null = reduced ? 3.6 : null;
  const shoulder = new THREE.Vector3();
  let poseMs = 0;
  function pose(): void {
    const t0 = performance.now();
    for (const c of cooks) {
      c.st.work("lift", t, c.opts);
      c.st.toWorld(c.st.jointPos("upper_arm_R", shoulder), c.shell.position);
      c.shell.rotation.y = c.st.placement.yaw;
    }
    poseMs = 0.95 * poseMs + 0.05 * (performance.now() - t0);
  }

  // A hook for checking the page by script: set a view, or hold the lift at a moment.
  Object.assign(stageEl, {
    brigade: {
      view(eye: [number, number, number], at: [number, number, number]) {
        userMoved = true;
        camera.position.set(...eye);
        controls.target.set(...at);
        controls.update();
      },
      /** The trunk samples' rear reach against the skinned surface's, this pose (world z, metres). */
      measure() {
        scene.updateMatrixWorld(true);
        const v = new THREE.Vector3();
        return cooks.map((c) => {
          let skin = -Infinity;
          for (const m of c.st.meshes) {
            const pos = m.geometry.attributes.position as THREE.BufferAttribute;
            for (let i = 0; i < pos.count; i += 3) {
              m.getVertexPosition(i, v);
              skin = Math.max(skin, v.applyMatrix4(m.matrixWorld).z);
            }
          }
          return { cook: c.spec.name, sampled: +rearNow(c).toFixed(3), back: +rearSkin(c).z.toFixed(3), skin: +skin.toFixed(3), rearZ: +c.rearZ.toFixed(3) };
        });
      },
      poseMs: () => poseMs,
      info: () => ({
        mountMs: Math.round(mountMs),
        tris: renderer.info.render.triangles,
        calls: renderer.info.render.calls,
        panSteelKg: +PAN_STEEL_KG.toFixed(3),
        panCapacityL: +(panVolume(PAN_EXT.h - 0.008) * 1000).toFixed(2),
        cooks: cooks.map((c) => ({ cook: c.spec.name, shellR: +c.shellR.toFixed(3), rearZ: +c.rearZ.toFixed(3), rearY: +c.rearY.toFixed(3), fill: +c.fill.toFixed(3), root: +c.st.placement.z.toFixed(3) })),
      }),
      hold(at: number | null) {
        frozen = at;
        if (at !== null) t = at;
      },
    },
  });

  syncAisle();
  syncLoad();
  syncReach();
  pose();
  readouts();
  void gate.refresh(scene);
  // When a cook's detailed body arrives, measure the clearance again on its surface (and compile its shaders before
  // the next draw).
  for (const c of cooks)
    void c.st.ready.then(() => {
      if (c.st.body !== "detailed") return;
      void gate.refresh(scene);
      c.skinBack = backVertices(c.st);
      measureRear(c);
      keepBoxes();
      frame();
      syncAisle();
      pose();
    });

  // The stage draws only while it is on screen. On a phone the crew table sits below the stage, so the lift keeps
  // going (posed, not drawn) while the table is read, and its numbers stay live.
  let inView = false;
  let tableInView = false;
  observeInView(stageEl, (v) => {
    inView = v;
  });
  const tableEl = root.querySelector("table.cohorts");
  if (tableEl)
    observeInView(tableEl, (v) => {
      tableInView = v;
    });
  let last = performance.now();
  let sinceReadout = 0;
  function animate(now: number): void {
    requestAnimationFrame(animate);
    const dt = Math.max(0, Math.min(0.05, (now - last) / 1000));
    last = now;
    if (!(inView || tableInView) || document.hidden) return;
    if (frozen === null) t += dt;
    else t = frozen;
    pose();
    sinceReadout += dt;
    if (sinceReadout > 0.12) {
      sinceReadout = 0;
      readouts();
    }
    if (!inView || !gate.ready()) return;
    for (const r of ranges) r.hooks.tick?.(now / 1000, dt, 1);
    controls.update();
    renderer.render(scene, camera);
  }
  requestAnimationFrame(animate);
}
