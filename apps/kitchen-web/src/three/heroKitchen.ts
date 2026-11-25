/**
 * HERO STAGE (#hero): a commercial kitchen built live from its manifest.
 *
 * Every fixture in the exemplar's scene3d manifest is built by the procedural model library at its manifest extents
 * and placed at its manifest pose, inside a room shell of the manifest's size: quarry tile on the floor, FRP on the
 * walls (the walls facing the camera cut away), the site's floor drains, anti-fatigue mats at the stations. Static
 * geometry is merged per material across fixtures; anything that moves (burner valves, fryer baskets, doors people
 * open, the speed rack people push) stays live.
 *
 * Staff work the kitchen's stations. Where they may stand and walk, and when, is crowd.ts: a clearance field from
 * the room and the fixtures' footprints, A* paths kept to the right of the aisle, and space-time legs committed only
 * when they keep every body clear of every other and of the room (see crowd.ts for the invariants).
 *
 * Overlays are computed from the telemetry kernel (computeKitchenTelemetry): fixture surface temperatures from the
 * gas and electric loads over the fixtures' ratings, plumes from each appliance's convective heat rising as a point
 * source plume into its hood (and spilling at the uncaptured fraction), and the aisle clearance field itself.
 */

import * as THREE from "three";
import { mergeGeometries } from "three/examples/jsm/utils/BufferGeometryUtils.js";
import { EXEMPLARS, FLAGSHIP_EQUIPMENT, computeKitchenTelemetry, type CameraPresetId, type ExemplarId, type ExemplarSpec, type HeroMetrics } from "../sim/exemplars";
import { ASSET_BASE, createScrollSafeOrbit, createStudioRenderer, fitClipPlanes, getSoftParticleTexture, observeInView, shaderGate, thermalColor, type Scene3dFixture } from "./kit";
import { buildFixture, type Model } from "./models/all";
import { palette, perforation } from "./models/materials";
import { placeFixture } from "./models/shop";
import { studioEnvironment, studioLights } from "./models/stage";
import { Crowd, NavGrid, PAIR_MARGIN, STATIC_MARGIN, rectDistance, type Footprint, type LegCandidate, type LegVerdict, type Motion, type Rect, type XZ } from "./models/crowd";
import { HOTEL_PAN, PICK_S, PLACE_S, Routine, TOOLS, copyPose, makeStaff, newPose, panHoldZ, panSeconds, pickShortfall, standPose, type BodyPose, type Cohort, type RoutineStep, type Staff, type WorkKind, type WorkOpts } from "./models/people";
import { combiCavity } from "./models/ovens";

type OverlayMode = "pbr" | "thermal" | "plumes" | "aisles";

export interface HeroController {
  /** Settles once the first kitchen is built (or could not be). */
  ready: Promise<void>;
  setExemplar(id: ExemplarId): Promise<void>;
  getLoadedExemplarId(): ExemplarId;
  getDrawCalls(): number;
}

/* ------------------------------------------------------------------------------------------------------------ */
/* The manifest                                                                                                   */
/* ------------------------------------------------------------------------------------------------------------ */

interface Socket {
  socket_id: string;
  kind: string;
  local_pos_m: [number, number, number];
  forward_dir: [number, number, number];
}
interface Fixture extends Scene3dFixture {
  readonly sockets?: readonly Socket[];
}
interface Manifest {
  room_shell: { width_m: number; depth_m: number; ceiling_height_m: number };
  fixtures: Fixture[];
}
interface SiteMeta {
  walk_mesh_corridors?: { corridor_id: string; min_m: number[]; max_m: number[] }[];
  floor_drains_m?: { x_m: number; y_m: number }[];
  aisle_clearances_m?: Record<string, unknown>;
  mua_diffusers?: { diffuser_id: string; cfm: number; pos_m: number[] }[];
}

async function fetchJson<T>(url: string): Promise<T | null> {
  try {
    const r = await fetch(url);
    return r.ok ? ((await r.json()) as T) : null;
  } catch {
    return null;
  }
}

/* ------------------------------------------------------------------------------------------------------------ */
/* What the fixtures draw and hold, from the specs                                                               */
/* ------------------------------------------------------------------------------------------------------------ */

/**
 * Utility ratings: gas and electric input (W) and hood exhaust (CFM), transcribed from the `utility_ports` column of
 * assets/kitchens/specs/equipment_dimensions.tsv.
 */
const RATING: Record<string, { gas?: number; elec?: number; cfm?: number }> = {
  "NSF-RANGE-GAS-6B": { gas: 61500 },
  "NSF-SALAMANDER-36": { gas: 10300 },
  "NSF-GRIDDLE-24": { gas: 17600 },
  "NSF-GRIDDLE-36": { gas: 26400 },
  "NSF-GRIDDLE-48": { gas: 35200 },
  "NSF-CHARBROILER-36": { gas: 35200 },
  "NSF-FRYER-TWIN-50LB": { gas: 35200 },
  "NSF-INDUCTION-SUITE-4Z": { elec: 14000 },
  "NSF-COMBI-6PAN": { elec: 10800 },
  "NSF-COMBI-10PAN": { elec: 18900 },
  "NSF-COMBI-20PAN-ROLLIN": { elec: 37200 },
  "NSF-KETTLE-40GAL": { elec: 18000 },
  "NSF-TILT-SKILLET-40GAL": { elec: 18000 },
  "NSF-MICROWAVE-COMM": { elec: 3200 },
  "NSF-TOASTER-CONVEYOR": { elec: 2700 },
  "NSF-SOUSVIDE-BATH": { elec: 1100 },
  "NSF-HEATLAMP-STRIP-48": { elec: 1100 },
  "NSF-HOLDING-CABINET": { elec: 2000 },
  "NSF-STEAMTABLE-HOT-60": { elec: 3000 },
  "NSF-DISHWASHER-DOOR": { elec: 12500 },
  "NSF-DISHWASHER-CONVEYOR": { elec: 26000 },
  "NSF-BOOSTER-HEATER-36KW": { elec: 36000 },
  "NSF-WALKIN-COOLER": { elec: 2200 },
  "NSF-WALKIN-COOLER-14X16": { elec: 4200 },
  "NSF-WALKIN-FREEZER-8X10": { elec: 3400 },
  "NSF-WALKIN-FREEZER-10X14": { elec: 3800 },
  "NSF-BLAST-CHILLER-15PAN": { elec: 5200 },
  "NSF-REACHIN-2D": { elec: 680 },
  "NSF-UNDERCOUNTER-2D": { elec: 420 },
  "NSF-MAKELINE-COLD-72": { elec: 740 },
  "NSF-ICEMACHINE-AIRCOOLED": { elec: 1450 },
  "NSF-SLICER-ROTARY": { elec: 420 },
  "NSF-PROCESSOR-SBLADE": { elec: 1100 },
  "NSF-HOOD-TYPE1-CANOPY": { cfm: 2800 },
  "NSF-HOOD-TYPE1-CANOPY-144IN": { cfm: 3600 },
  "NSF-HOOD-TYPE1-CANOPY-168IN": { cfm: 4200 },
  "NSF-HOOD-TYPE1-CANOPY-288IN": { cfm: 9800 },
  "NSF-HOOD-TYPE1-ISLAND-144IN": { cfm: 5200 },
  "NSF-HOOD-TYPE2-CONDENSATE": { cfm: 1100 },
};

/** Heat-producing electric fixtures (refrigeration rejects its heat at the condensers, drawn as warm only). */
const isElectricHeat = (id: string) => /COMBI|KETTLE|TILT|MICROWAVE|TOASTER|SOUSVIDE|HEATLAMP|HOLDING|STEAMTABLE|DISHWASHER|BOOSTER|INDUCTION/.test(id);
const isOpenFlame = (id: string) => /RANGE|CHARBROILER|SALAMANDER/.test(id);
const isCold = (id: string) => /WALKIN|REACHIN|UNDERCOUNTER|MAKELINE|BLAST|ICEMACHINE/.test(id);
const isHood = (id: string) => /HOOD/.test(id);

/**
 * Each fixture's operating temperature, °C: the flagship table's nominal where it lists the archetype or its
 * family, else the setpoint the class is built to (cold.ts thermometer defaults; NSF/ANSI 3 final rinse 82 °C;
 * FDA Food Code hot holding).
 */
function setPointC(id: string): number | null {
  const flag = (key: string) => FLAGSHIP_EQUIPMENT.find((e) => e.id === key)?.nominalTempC ?? null;
  const exact = flag(id);
  if (exact !== null) return exact;
  if (/GRIDDLE/.test(id)) return flag("NSF-GRIDDLE-48");
  if (/COMBI-20/.test(id)) return flag("NSF-COMBI-20PAN-ROLLIN");
  if (/COMBI/.test(id)) return flag("NSF-COMBI-10PAN");
  if (/WALKIN-COOLER/.test(id)) return flag("NSF-WALKIN-COOLER");
  if (/WALKIN-FREEZER/.test(id)) return -20;
  if (/HOOD-TYPE1/.test(id)) return flag("NSF-HOOD-TYPE1-CANOPY");
  if (/HOLDING|CARRIER/.test(id)) return flag("NSF-CARRIER-INSULATED");
  if (/STOCKPOT/.test(id)) return flag("NSF-STOCKPOT-40QT");
  if (/SKILLET-\d/.test(id)) return flag("NSF-SKILLET-12IN");
  if (/SALAMANDER/.test(id)) return flag("NSF-CHARBROILER-36");
  if (/REACHIN|UNDERCOUNTER/.test(id)) return 3;
  if (/MAKELINE/.test(id)) return 2;
  if (/ICEMACHINE/.test(id)) return 0;
  if (/DISHWASHER|BOOSTER/.test(id)) return 82;
  if (/STEAMTABLE|HEATLAMP/.test(id)) return flag("NSF-CARRIER-INSULATED");
  if (/SOUSVIDE/.test(id)) return 60;
  if (/TOASTER/.test(id)) return flag("NSF-GRIDDLE-48");
  if (/MICROWAVE/.test(id)) return 60;
  return null;
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Thermal IR: one shader for every fixture surface, temperatures by fixture index                              */
/* ------------------------------------------------------------------------------------------------------------ */

const MAX_FIX = 128;
const T_MIN = -25;
const T_MAX = 340;

function thermalMaterial(): THREE.ShaderMaterial {
  return new THREE.ShaderMaterial({
    uniforms: { uT: { value: new Float32Array(MAX_FIX) }, uMin: { value: T_MIN }, uMax: { value: T_MAX } },
    vertexShader: /* glsl */ `
      attribute float aFix;
      attribute float aHeat;
      uniform float uT[${MAX_FIX}];
      varying float vT;
      varying vec3 vN;
      void main() {
        float room = uT[0];
        vT = room + (uT[int(aFix + 0.5)] - room) * aHeat;
        vN = normalize(normalMatrix * normal);
        gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0);
      }`,
    fragmentShader: /* glsl */ `
      uniform float uMin;
      uniform float uMax;
      varying float vT;
      varying vec3 vN;
      vec3 ironbow(float t) {
        if (t < 0.25) { float u = t / 0.25; return vec3(0.06 + u * 0.06, 0.14 + u * 0.24, 0.38 + u * 0.34); }
        if (t < 0.52) { float u = (t - 0.25) / 0.27; return vec3(0.12 + u * 0.68, 0.38 - u * 0.22, 0.72 - u * 0.46); }
        if (t < 0.8) { float u = (t - 0.52) / 0.28; return vec3(0.8 + u * 0.18, 0.16 + u * 0.48, 0.26 - u * 0.22); }
        float u = (t - 0.8) / 0.2; return vec3(0.98 + u * 0.02, 0.64 + u * 0.32, 0.04 + u * 0.72);
      }
      void main() {
        float t = clamp((vT - uMin) / (uMax - uMin), 0.0, 1.0);
        float shade = 0.72 + 0.28 * abs(normalize(vN).z);
        gl_FragColor = vec4(ironbow(t) * shade, 1.0);
        #include <colorspace_fragment>
      }`,
  });
}

/** Tag a geometry with its fixture index and its surfaces' share of the fixture's excess temperature. */
function tagThermal(g: THREE.BufferGeometry, fix: number, heat: (x: number, y: number, z: number) => number): void {
  const pos = g.getAttribute("position");
  const n = pos.count;
  const f = new Float32Array(n).fill(fix);
  const h = new Float32Array(n);
  for (let i = 0; i < n; i++) h[i] = heat(pos.getX(i), pos.getY(i), pos.getZ(i));
  g.setAttribute("aFix", new THREE.BufferAttribute(f, 1));
  g.setAttribute("aHeat", new THREE.BufferAttribute(h, 1));
}

/* ------------------------------------------------------------------------------------------------------------ */
/* The kitchen                                                                                                    */
/* ------------------------------------------------------------------------------------------------------------ */

/** A placed fixture. Frame conversions use its group's matrix (manifest local [x across, y back, z up]). */
export interface Fx {
  i: number;
  f: Fixture;
  id: string;
  model: Model;
  /** Floor footprint (three xz), or null for fixtures above head height. */
  rect: Rect | null;
  setC: number | null;
  /** Rated gas and electric input, W. */
  gasW: number;
  elecW: number;
  /** Computed each telemetry update. */
  tempC: number;
  firing: number;
  /** The convective heat it gives off as a plume, W. */
  plumeW: number;
  /** For shelving people take pans from: the bounding boxes of its triangles (fixture frame), six floats each. */
  tris?: Float32Array;
}

/** A manifest-local point (x across, y from the front face back, z up) in world coordinates. */
export function fxLocal(fx: Fx, lx: number, ly: number, lz: number, out = new THREE.Vector3()): THREE.Vector3 {
  return out.set(lx, lz, -ly).applyMatrix4(fx.model.group.matrixWorld);
}
/** A fixture's socket by id. */
export function fxSocket(fx: Fx, id: string): Socket | undefined {
  return fx.f.sockets?.find((s) => s.socket_id === id);
}

interface Wall {
  mesh: THREE.Mesh;
  /** Outward normal in xz and a point on the inside face. */
  nx: number;
  nz: number;
  px: number;
  pz: number;
}

interface Kitchen {
  spec: ExemplarSpec;
  manifest: Manifest;
  meta: SiteMeta | null;
  root: THREE.Group;
  W: number;
  D: number;
  H: number;
  fixtures: Fx[];
  walls: Wall[];
  grid: NavGrid;
  /** The room with only its walls and the fixtures taller than a pan carried at the waist. */
  tall: NavGrid;
  merged: THREE.Mesh[];
  /** Every mesh whose material the thermal view swaps (fixtures and shell). */
  thermalMeshes: THREE.Mesh[];
  tris: number;
  hoodLights: THREE.PointLight[];
  dispose(): void;
}

const tileUV = 1 / 0.608;

function shellMaterials() {
  const p = palette();
  return { floor: p.tile, wall: p.frp, cove: p.tile, rubber: p.rubber, steel: p.ssDull, dark: p.charcoal };
}

/** The floor, walls with a tiled cove base, the site's floor drains. */
function buildShell(K: { W: number; D: number; H: number; meta: SiteMeta | null }, group: THREE.Group, thermal: THREE.Mesh[]): Wall[] {
  const M = shellMaterials();
  const { W, D, H } = K;
  const floorGeo = new THREE.PlaneGeometry(W, D, 1, 1).rotateX(-Math.PI / 2).translate(W / 2, 0, -D / 2);
  const uv = floorGeo.getAttribute("uv");
  const pos = floorGeo.getAttribute("position");
  for (let i = 0; i < pos.count; i++) uv.setXY(i, pos.getX(i) * tileUV, -pos.getZ(i) * tileUV);
  tagThermal(floorGeo, 0, () => 0);
  const floor = new THREE.Mesh(floorGeo, M.floor);
  floor.receiveShadow = true;
  floor.name = "floor";
  group.add(floor);
  thermal.push(floor);

  // Walls: FRP above a 150 mm quarry-tile cove, 100 mm thick, standing outside the room's inside faces.
  const T = 0.1;
  const walls: Wall[] = [];
  const sides = [
    { nx: 0, nz: 1, px: W / 2, pz: 0, len: W + 2 * T }, // front (z = 0)
    { nx: 0, nz: -1, px: W / 2, pz: -D, len: W + 2 * T }, // back
    { nx: -1, nz: 0, px: 0, pz: -D / 2, len: D }, // west (x = 0)
    { nx: 1, nz: 0, px: W, pz: -D / 2, len: D }, // east
  ];
  for (const s of sides) {
    const along = s.nx === 0;
    const g = new THREE.BoxGeometry(along ? s.len : T, H, along ? T : s.len);
    g.translate(s.px + s.nx * (T / 2), H / 2, s.pz + s.nz * (T / 2));
    tagThermal(g, 0, () => 0);
    const mesh = new THREE.Mesh(g, M.wall);
    mesh.receiveShadow = true;
    mesh.castShadow = false;
    group.add(mesh);
    thermal.push(mesh);
    const cove = new THREE.BoxGeometry(along ? s.len - 2 * T : 0.012, 0.15, along ? 0.012 : s.len);
    cove.translate(s.px - s.nx * 0.006, 0.075, s.pz - s.nz * 0.006);
    const cuv = cove.getAttribute("uv");
    const cp = cove.getAttribute("position");
    for (let i = 0; i < cp.count; i++) cuv.setXY(i, (along ? cp.getX(i) : cp.getZ(i)) * tileUV, cp.getY(i) * tileUV);
    tagThermal(cove, 0, () => 0);
    const cm = new THREE.Mesh(cove, M.cove);
    cm.receiveShadow = true;
    group.add(cm);
    thermal.push(cm);
    walls.push({ mesh, nx: s.nx, nz: s.nz, px: s.px, pz: s.pz });
  }

  // Floor drains: a 300 mm square nickel-bronze grate over a dark sump, flush with the tile.
  const drains = K.meta?.floor_drains_m ?? [];
  if (drains.length) {
    const bars: THREE.BufferGeometry[] = [];
    const sumps: THREE.BufferGeometry[] = [];
    for (const d of drains) {
      const x = d.x_m;
      const z = -d.y_m;
      sumps.push(new THREE.PlaneGeometry(0.3, 0.3).rotateX(-Math.PI / 2).translate(x, 0.0012, z));
      const frame = (w: number, dd: number, ox: number, oz: number) => bars.push(new THREE.BoxGeometry(w, 0.004, dd).translate(x + ox, 0.002, z + oz));
      frame(0.3, 0.018, 0, 0.141);
      frame(0.3, 0.018, 0, -0.141);
      frame(0.018, 0.3, 0.141, 0);
      frame(0.018, 0.3, -0.141, 0);
      for (let k = -5; k <= 5; k++) frame(0.006, 0.27, k * 0.024, 0);
    }
    const sump = new THREE.Mesh(mergeGeometries(sumps)!, new THREE.MeshStandardMaterial({ color: 0x0c0b0a, roughness: 0.9 }));
    const grate = new THREE.Mesh(mergeGeometries(bars)!, palette().ssBright);
    grate.receiveShadow = true;
    group.add(sump, grate);
  }
  return walls;
}

/** Where the walls of the shell hide: any wall the camera is outside of. */
function cutWalls(walls: Wall[], cam: THREE.Vector3): void {
  for (const w of walls) w.mesh.visible = (cam.x - w.px) * w.nx + (cam.z - w.pz) * w.nz < 0;
}

/** Archetypes whose moving parts people work while the scene runs: their joints stay live, the rest is merged. */
const LIVE = /RANGE|FRYER|REACHIN|UNDERCOUNTER|WALKIN|SPEEDRACK|SALAMANDER/;
/** Fixtures people move bodily (speed racks on their casters): merged into themselves, not into the room. */
const MOVABLE = /SPEEDRACK/;
/** Shelving people take pans from and put them back on: its triangles are kept to find the free places. */
const SHELVES = /WALKIN-COOLER|WIRE-SHELVING/;
/** Fixtures at counter height or with only a back part standing higher: a pan carried at the waist passes over them. */
const LOW = /RANGE|GRIDDLE|FRYER|CHARBROILER|INDUCTION|STEAMTABLE|MAKELINE|PREPTABLE|SINK|KETTLE|TILT|HEATLAMP|TOASTER|UNDERCOUNTER|MICROWAVE|SOUSVIDE|PORTION|BOARD/;

/** The walk-in coolers' doors and panels (cold.ts walkIn): the door's left edge and width along the front wall. */
const WALKINS: Record<string, { door0: number; doorW: number }> = {
  "NSF-WALKIN-COOLER": { door0: 0.76, doorW: 0.914 },
  "NSF-WALKIN-COOLER-14X16": { door0: 1.6, doorW: 1.067 },
};
const PANEL = 0.1016;
const WALL_FRONT = -0.087;

/** A point in a fixture's own frame (x across, y up, z from its front face at 0 back to -d) in world coordinates. */
function fxThree(fx: Fx, x: number, y: number, z: number, out = new THREE.Vector3()): THREE.Vector3 {
  return out.set(x, y, z).applyMatrix4(fx.model.group.matrixWorld);
}
/** A box on the floor in a fixture's frame, as a world rectangle. */
function fxRect(fx: Fx, x0: number, x1: number, z0: number, z1: number): Rect {
  const r: Rect = { minX: Infinity, minZ: Infinity, maxX: -Infinity, maxZ: -Infinity };
  const v = new THREE.Vector3();
  for (const x of [x0, x1])
    for (const z of [z0, z1]) {
      fxThree(fx, x, 0, z, v);
      r.minX = Math.min(r.minX, v.x);
      r.maxX = Math.max(r.maxX, v.x);
      r.minZ = Math.min(r.minZ, v.z);
      r.maxZ = Math.max(r.maxZ, v.z);
    }
  return r;
}
/** The units of shelving down each side of a walk-in (cold.ts walkIn): their lengths, front to back. */
function walkInUnits(w: number, d: number): { Ds: number; lens: number[]; z0: number } {
  const T = PANEL;
  const zi = WALL_FRONT - T;
  const zb = -d + T;
  const Ds = w < 2.7 ? 0.457 : 0.61;
  let left = zi - 0.06 - (zb + 0.03);
  const lens: number[] = [];
  for (const L of [1.219, 1.219, 1.219, 0.914, 0.762, 0.61])
    if (L <= left + 0.001 && lens.length < (w > 2.7 ? 2 : 3)) {
      lens.push(L);
      left -= L + 0.01;
    }
  return { Ds, lens, z0: zi - 0.06 };
}
/**
 * A walk-in cooler as people meet it: its panels, the shelving down its sides and the unit coolers' corner at the
 * back are obstacles; the floor between, and the doorway, are free (and gated: see Crew).
 */
function walkInRects(fx: Fx): Rect[] {
  const [w, d] = fx.f.local_extents_m;
  const k = WALKINS[fx.id]!;
  const T = PANEL;
  const zi = WALL_FRONT - T;
  const zb = -d + T;
  const { Ds, lens, z0 } = walkInUnits(w, d);
  const run = lens.reduce((a, b) => a + b + 0.01, 0);
  return [
    fxRect(fx, 0, T, -d, 0),
    fxRect(fx, w - T, w, -d, 0),
    fxRect(fx, 0, w, -d, zb + 0.85),
    fxRect(fx, 0, k.door0, zi, 0),
    fxRect(fx, k.door0 + k.doorW, w, zi, 0),
    fxRect(fx, T, T + 0.012 + Ds, z0 - run, z0),
    fxRect(fx, w - T - 0.012 - Ds, w - T, z0 - run, z0),
  ];
}
/** Small things on counters that people may pick up: drawn but not obstacles to walking. */
const SMALLWARE = /SKILLET-\d|STOCKPOT|HOTEL-PAN|CUTTING-BOARD|CHEF-KNIFE|GRATER|PORTION-SCALE|PROCESSOR|SLICER|MICROWAVE|TOASTER|SOUSVIDE/;

/** The share of a fixture's excess temperature at a point of its surface (fixture frame: y up, z from -d to 0). */
function heatProfile(id: string, ext: { w: number; d: number; h: number }): (x: number, y: number, z: number) => number {
  if (isCold(id)) return () => 0.3;
  if (isHood(id)) return () => 0.65;
  if (/SALAMANDER|HEATLAMP/.test(id)) return () => 1;
  if (/COMBI|HOLDING|MICROWAVE|DISHWASHER|BOOSTER|TOASTER/.test(id)) return (_x, y) => (y > ext.h * 0.12 ? 0.5 : 0.2);
  if (/RANGE|GRIDDLE|CHARBROILER|FRYER|INDUCTION|KETTLE|TILT|STEAMTABLE|SOUSVIDE|SKILLET-\d|STOCKPOT/.test(id)) {
    const top = ext.h > 1.0 ? 0.914 : ext.h;
    return (_x, y, z) => (y > top - 0.05 && y < top + 0.25 && z < -0.035 ? 1 : 0.24);
  }
  return () => 0;
}

/** Anti-fatigue mats: perforated rubber in front of the stations people stand at, kept clear of other fixtures. */
const MATTED = /RANGE|GRIDDLE|FRYER|CHARBROILER|INDUCTION|PREPTABLE|MAKELINE|STEAMTABLE|SINK-3COMP|VEG-PREP|KETTLE|TILT|DISHWASHER/;
function buildMats(fixtures: Fx[], grid: NavGrid): THREE.Object3D {
  const geos: THREE.BufferGeometry[] = [];
  const v = new THREE.Vector3();
  for (const fx of fixtures) {
    if (!MATTED.test(fx.id) || !fxSocket(fx, "stance_front")) continue;
    const [w] = fx.f.local_extents_m;
    const mw = Math.max(0.6, w - 0.04);
    const md = 0.6;
    // Trim the mat's depth until all of it lies on free floor.
    let depth = md;
    const ok = (dd: number) => {
      for (const [lx, ly] of [[0.02, -0.02], [mw + 0.02, -0.02], [0.02, -dd - 0.02], [mw + 0.02, -dd - 0.02], [mw / 2 + 0.02, -dd - 0.02]] as const) {
        fxLocal(fx, lx + (w - mw) / 2 - 0.02, ly, 0, v);
        if (grid.clearance(v.x, v.z) < 0.005) return false;
      }
      return true;
    };
    while (depth > 0.3 && !ok(depth)) depth -= 0.05;
    if (!ok(depth)) continue;
    const g = new THREE.BoxGeometry(mw, 0.004, depth);
    g.translate(w / 2, 0.002, depth / 2 + 0.03);
    g.applyMatrix4(fx.model.group.matrixWorld);
    const pos = g.getAttribute("position");
    const uv = g.getAttribute("uv");
    for (let i = 0; i < pos.count; i++) uv.setXY(i, pos.getX(i) / 0.42, pos.getZ(i) / 0.42);
    geos.push(g);
  }
  if (!geos.length) return new THREE.Group();
  const alpha = perforation(14, 0.5);
  alpha.wrapS = alpha.wrapT = THREE.RepeatWrapping;
  const mat = new THREE.MeshStandardMaterial({ color: 0x1a1817, roughness: 0.93, metalness: 0, alphaMap: alpha, alphaTest: 0.5 });
  const mesh = new THREE.Mesh(mergeGeometries(geos)!, mat);
  mesh.receiveShadow = true;
  mesh.name = "mats";
  return mesh;
}

interface BuildOpts {
  phone: boolean;
  onProgress?: (done: number, total: number) => void;
  /** Called with each fixture before it is merged, to set its state and say which joints stay live. */
  prepare?: (fx: Fx) => void;
}

async function buildKitchen(spec: ExemplarSpec, manifest: Manifest, meta: SiteMeta | null, o: BuildOpts, alive: () => boolean): Promise<Kitchen | null> {
  const root = new THREE.Group();
  root.name = `kitchen-${spec.id}`;
  const W = manifest.room_shell.width_m;
  const D = manifest.room_shell.depth_m;
  const H = manifest.room_shell.ceiling_height_m;
  const thermalMeshes: THREE.Mesh[] = [];
  const walls = buildShell({ W, D, H, meta }, root, thermalMeshes);
  const fixtures: Fx[] = [];
  const staticRoot = new THREE.Group();
  root.add(staticRoot);
  const bins = new Map<string, { mat: THREE.Material; geos: THREE.BufferGeometry[]; cast: boolean; recv: boolean; order: number }>();
  let tris = 0;
  const big = manifest.fixtures.length > 40;
  let k = 0;
  let slice = performance.now();
  for (const f of manifest.fixtures) {
    if (!alive()) return null;
    const [w, d, h] = f.local_extents_m;
    const ext = { w, d, h };
    const id = f.archetype_id;
    const live = LIVE.test(id);
    // Hot line and the fixtures people work at are drawn finest on desktop when the room is not the commissary.
    const detail = o.phone ? (live ? 1 : 0) : live && !big ? 2 : 1;
    const model = buildFixture(id, ext, detail as 0 | 1 | 2, true);
    placeFixture(model.group, f.pose);
    root.add(model.group);
    model.group.updateMatrixWorld(true);
    const box = f.world_aabb_m;
    const rect: Rect | null =
      box && box.min_m[2] < 1.95 && !SMALLWARE.test(id) ? { minX: box.min_m[0], maxX: box.max_m[0], minZ: -box.max_m[1], maxZ: -box.min_m[1] } : null;
    const rating = RATING[id] ?? {};
    const fx: Fx = { i: fixtures.length + 1, f, id, model, rect, setC: setPointC(id), gasW: rating.gas ?? 0, elecW: isElectricHeat(id) ? (rating.elec ?? 0) : 0, tempC: 20, firing: 0, plumeW: 0 };
    fixtures.push(fx);
    model.hooks.set?.("lights", 1);
    o.prepare?.(fx);
    const heat = heatProfile(id, ext);
    // A combi's door swings when a pan goes in or comes out: that joint alone stays live (its fans and dial merge).
    const doorJ = /COMBI/.test(id) ? model.joints.get("cavity_door") : undefined;
    const liveJoints = new Set<THREE.Object3D>(live ? [...model.joints.values()] : doorJ ? [doorJ] : []);
    const toMerge: THREE.Mesh[] = [];
    const ownMerge: THREE.Mesh[] = [];
    const movable = MOVABLE.test(id);
    // A walk-in's roof is left off so that the people in it can be seen.
    if (/WALKIN-COOLER/.test(id)) model.hooks.set?.("roof", 1);
    // The cold store's shelving in the room is a third stocked: the rest of it holds the pans the cooks bring.
    if (/WIRE-SHELVING/.test(id) && /cold/i.test(fx.f.fixture_id)) model.hooks.set?.("stock", 0.34);
    const boxes: number[] = [];
    const fxInv = new THREE.Matrix4().copy(model.group.matrixWorld).invert();
    model.group.traverse((obj) => {
      const m = obj as THREE.Mesh;
      if (!m.isMesh || (m as THREE.InstancedMesh).isInstancedMesh || (m as THREE.SkinnedMesh).isSkinnedMesh) return;
      const g = m.geometry as THREE.BufferGeometry;
      if (!g.getAttribute("position")) return;
      let shown = true;
      for (let q: THREE.Object3D | null = m; q && q !== model.group; q = q.parent) if (!q.visible) shown = false;
      // Hidden parts (stock left off a shelf) are neither merged nor in the way.
      if (!shown && !live) return;
      if (SHELVES.test(id) && shown) {
        // Every triangle's box in the fixture's frame (stock on the shelves included), to find room for a pan.
        const M = new THREE.Matrix4().copy(fxInv).multiply(m.matrixWorld);
        const pos = g.getAttribute("position");
        const idx = g.index;
        const n = idx ? idx.count : pos.count;
        const a = new THREE.Vector3();
        for (let t = 0; t + 2 < n; t += 3) {
          let x0 = Infinity, y0 = Infinity, z0 = Infinity, x1 = -Infinity, y1 = -Infinity, z1 = -Infinity;
          for (let q = 0; q < 3; q++) {
            a.fromBufferAttribute(pos, idx ? idx.getX(t + q) : t + q).applyMatrix4(M);
            x0 = Math.min(x0, a.x); y0 = Math.min(y0, a.y); z0 = Math.min(z0, a.z);
            x1 = Math.max(x1, a.x); y1 = Math.max(y1, a.y); z1 = Math.max(z1, a.z);
          }
          if (y1 > 0.15 && y0 < 2.0) boxes.push(x0, y0, z0, x1, y1, z1);
        }
      }
      // Thermal tags in the fixture's own frame.
      const inv = new THREE.Matrix4().copy(model.group.matrixWorld).invert().multiply(m.matrixWorld);
      const lp = new THREE.Vector3();
      tagThermal(g, fx.i, (x, y, z) => {
        lp.set(x, y, z).applyMatrix4(inv);
        return heat(lp.x, lp.y, lp.z);
      });
      tris += (g.index ? g.index.count : g.getAttribute("position").count) / 3;
      let p: THREE.Object3D | null = m.parent;
      let underLive = false;
      while (p && p !== model.group) {
        if (liveJoints.has(p)) underLive = true;
        p = p.parent;
      }
      const mat = m.material as THREE.Material;
      if (/CHEF-KNIFE/.test(id) || underLive || Array.isArray(m.material) || (mat as THREE.ShaderMaterial).isShaderMaterial || !(m.parent === model.group || !live)) {
        thermalMeshes.push(m);
        return;
      }
      (movable ? ownMerge : toMerge).push(m);
    });
    if (boxes.length) fx.tris = Float32Array.from(boxes);
    if (ownMerge.length) {
      // A rack's frame merged per material in its own frame, so that it can roll away as one.
      const own = new Map<THREE.Material, THREE.BufferGeometry[]>();
      for (const m of ownMerge) {
        const mat = m.material as THREE.Material;
        const g = (m.geometry as THREE.BufferGeometry).clone().applyMatrix4(new THREE.Matrix4().copy(fxInv).multiply(m.matrixWorld));
        for (const name of Object.keys(g.attributes)) if (!["position", "normal", "uv", "aFix", "aHeat"].includes(name)) g.deleteAttribute(name);
        g.morphAttributes = {};
        if (!g.getAttribute("normal")) g.computeVertexNormals();
        if (!g.getAttribute("uv")) g.setAttribute("uv", new THREE.BufferAttribute(new Float32Array(g.getAttribute("position").count * 2), 2));
        if (!g.index) g.setIndex([...Array(g.getAttribute("position").count).keys()]);
        g.clearGroups();
        if (!own.has(mat)) own.set(mat, []);
        own.get(mat)!.push(g);
        m.parent?.remove(m);
        (m.geometry as THREE.BufferGeometry).dispose();
      }
      for (const [mat, geos] of own) {
        const g = geos.length === 1 ? geos[0]! : mergeGeometries(geos, false);
        if (!g) continue;
        const mesh = new THREE.Mesh(g, mat);
        mesh.castShadow = mesh.receiveShadow = true;
        mesh.name = "own-merged";
        model.group.add(mesh);
        thermalMeshes.push(mesh);
      }
    }
    for (const m of toMerge) {
      const mat = m.material as THREE.Material;
      const g = (m.geometry as THREE.BufferGeometry).clone().applyMatrix4(m.matrixWorld);
      for (const name of Object.keys(g.attributes)) if (!["position", "normal", "uv", "aFix", "aHeat", ...(mat as THREE.MeshStandardMaterial).vertexColors ? ["color"] : []].includes(name)) g.deleteAttribute(name);
      g.morphAttributes = {};
      if (!g.getAttribute("normal")) g.computeVertexNormals();
      if (!g.getAttribute("uv")) g.setAttribute("uv", new THREE.BufferAttribute(new Float32Array(g.getAttribute("position").count * 2), 2));
      if (!g.index) g.setIndex([...Array(g.getAttribute("position").count).keys()]);
      g.clearGroups();
      const key = `${mat.uuid}|${m.castShadow ? 1 : 0}${m.receiveShadow ? 1 : 0}|${Object.keys(g.attributes).sort().join(",")}`;
      let bin = bins.get(key);
      if (!bin) bins.set(key, (bin = { mat, geos: [], cast: m.castShadow, recv: m.receiveShadow, order: m.renderOrder }));
      bin.geos.push(g);
      m.parent?.remove(m);
      (m.geometry as THREE.BufferGeometry).dispose();
    }
    k++;
    o.onProgress?.(k, manifest.fixtures.length);
    // Let the page breathe between fixtures: hand the thread back once a slice has run 40 ms (one fixture can take
    // that alone the first time its finishes are painted).
    if (performance.now() - slice > 40) {
      await new Promise((r) => setTimeout(r, 0));
      slice = performance.now();
    }
  }
  const merged: THREE.Mesh[] = [];
  for (const bin of bins.values()) {
    const g = bin.geos.length === 1 ? bin.geos[0]! : mergeGeometries(bin.geos, false);
    if (!g) continue;
    if (bin.geos.length > 1) for (const x of bin.geos) x.dispose();
    g.computeBoundingSphere();
    const mesh = new THREE.Mesh(g, bin.mat);
    mesh.castShadow = bin.cast;
    mesh.receiveShadow = bin.recv;
    mesh.renderOrder = bin.order;
    mesh.matrixAutoUpdate = false;
    staticRoot.add(mesh);
    merged.push(mesh);
    thermalMeshes.push(mesh);
  }
  // Obstacles to walking: the fixtures below head height; a walk-in cooler's panels and shelving, not its floor. The
  // tall ones (not counters) also stop a pan carried at the waist.
  const rects: Rect[] = [];
  const tallRects: Rect[] = [];
  for (const x of fixtures) {
    if (!x.rect) continue;
    const parts = WALKINS[x.id] ? walkInRects(x) : [x.rect];
    rects.push(...parts);
    if (!LOW.test(x.id) && x.f.world_aabb_m && x.f.world_aabb_m.max_m[2] > 1.05) tallRects.push(...parts);
  }
  const bounds = { minX: 0, minZ: -D, maxX: W, maxZ: 0 };
  const cell = W * D > 150 ? 0.08 : W * D > 60 ? 0.05 : 0.04;
  const grid = new NavGrid(bounds, rects, cell);
  const tall = new NavGrid(bounds, tallRects, cell);
  root.add(buildMats(fixtures, grid));
  const hoodLights: THREE.PointLight[] = [];
  for (const fx of fixtures) {
    if (!/HOOD-TYPE1/.test(fx.id)) continue;
    const [w, d] = fx.f.local_extents_m;
    // A warm lamp under the hood, standing in for its vapour-proof globes.
    const L = new THREE.PointLight(0xffc48a, 0, Math.max(3.5, w * 1.2), 2);
    fxLocal(fx, w / 2, d * 0.55, -0.12, L.position);
    root.add(L);
    hoodLights.push(L);
  }
  return {
    spec,
    manifest,
    meta,
    root,
    W,
    D,
    H,
    fixtures,
    walls,
    grid,
    tall,
    merged,
    thermalMeshes,
    tris,
    hoodLights,
    dispose() {
      root.traverse((obj) => {
        const m = obj as THREE.Mesh;
        if (m.isMesh) m.geometry?.dispose();
      });
    },
  };
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Heat: what every fixture is doing at this rush, electrification and ambient                                  */
/* ------------------------------------------------------------------------------------------------------------ */

/** The share of input heat a fixture spills to the room as convection: the telemetry kernel's own coefficients. */
const SPILL = { gas: 0.58, induction: 0.14, elec: 0.28 };

function applyHeat(K: Kitchen, m: HeroMetrics): void {
  const gasIn = Math.max(0, m.thermalKw - m.electricalKw) * 1000;
  const elecIn = m.electricalKw * 1000;
  const gasRated = K.fixtures.reduce((s, f) => s + f.gasW, 0);
  const elecRated = K.fixtures.reduce((s, f) => s + f.elecW, 0);
  const fg = gasRated > 0 ? Math.min(1, gasIn / gasRated) : 0;
  const fe = elecRated > 0 ? Math.min(1, elecIn / elecRated) : 0;
  const Troom = m.thermalComfortC;
  for (const fx of K.fixtures) {
    const set = fx.setC;
    if (fx.gasW > 0) {
      fx.firing = fg;
      fx.plumeW = fx.gasW * fg * SPILL.gas;
    } else if (fx.elecW > 0) {
      fx.firing = fe;
      fx.plumeW = fx.elecW * fe * (/INDUCTION/.test(fx.id) ? SPILL.induction : SPILL.elec);
    } else {
      fx.firing = 0;
      fx.plumeW = 0;
    }
    if (set === null) fx.tempC = Troom;
    else if (isCold(fx.id)) fx.tempC = set;
    else if (isHood(fx.id)) fx.tempC = Troom + (set - Troom) * fg;
    else if (isOpenFlame(fx.id)) fx.tempC = Troom + (set - Troom) * Math.sqrt(fx.firing);
    else if (fx.gasW > 0 || fx.elecW > 0) fx.tempC = fx.firing > 0.02 ? set : Troom;
    else fx.tempC = Math.max(Troom, Math.min(set, Troom + (set - Troom) * fg));
  }
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Plumes: each appliance's convective heat rising into its hood                                                 */
/* ------------------------------------------------------------------------------------------------------------ */

interface PlumeSource {
  fx: Fx;
  /** Source centre (world) and its half-width. */
  x: number;
  y: number;
  z: number;
  r: number;
  /** The hood over it, as its world box (null: rises to the ceiling). */
  hood: { minX: number; maxX: number; minZ: number; maxZ: number; y: number; top: number; bx: number; bz: number } | null;
}

const N_PLUME = 2400;

class Plumes {
  readonly points: THREE.Points;
  private pos = new Float32Array(N_PLUME * 3);
  private alpha = new Float32Array(N_PLUME);
  private size = new Float32Array(N_PLUME);
  private age = new Float32Array(N_PLUME);
  private life = new Float32Array(N_PLUME);
  private src = new Int16Array(N_PLUME).fill(-1);
  private esc = new Uint8Array(N_PLUME);
  private off = new Float32Array(N_PLUME * 2);
  private acc: number[] = [];
  private sources: PlumeSource[] = [];
  capture = 0.95;
  ceiling = 3;

  constructor() {
    const g = new THREE.BufferGeometry();
    g.setAttribute("position", new THREE.BufferAttribute(this.pos, 3).setUsage(THREE.DynamicDrawUsage));
    g.setAttribute("aAlpha", new THREE.BufferAttribute(this.alpha, 1).setUsage(THREE.DynamicDrawUsage));
    g.setAttribute("aSize", new THREE.BufferAttribute(this.size, 1).setUsage(THREE.DynamicDrawUsage));
    const mat = new THREE.ShaderMaterial({
      uniforms: { uTex: { value: getSoftParticleTexture() }, uScale: { value: 600 }, uColor: { value: new THREE.Color(0xffb070) } },
      vertexShader: /* glsl */ `
        attribute float aAlpha;
        attribute float aSize;
        uniform float uScale;
        varying float vA;
        void main() {
          vA = aAlpha;
          vec4 mv = modelViewMatrix * vec4(position, 1.0);
          gl_PointSize = aSize * uScale / max(0.3, -mv.z);
          gl_Position = projectionMatrix * mv;
        }`,
      fragmentShader: /* glsl */ `
        uniform sampler2D uTex;
        uniform vec3 uColor;
        varying float vA;
        void main() {
          float a = texture2D(uTex, gl_PointCoord).a * vA;
          if (a < 0.004) discard;
          gl_FragColor = vec4(uColor, a);
        }`,
      transparent: true,
      depthWrite: false,
      blending: THREE.NormalBlending,
    });
    this.points = new THREE.Points(g, mat);
    this.points.frustumCulled = false;
    this.points.renderOrder = 5;
  }

  setSources(K: Kitchen): void {
    const hoods = K.fixtures.filter((f) => isHood(f.id) && f.f.world_aabb_m);
    this.ceiling = K.H;
    this.sources = [];
    for (const fx of K.fixtures) {
      if (!(fx.gasW > 0 || /COMBI|KETTLE|TILT|STEAMTABLE|DISHWASHER|INDUCTION|SOUSVIDE|TOASTER/.test(fx.id))) continue;
      const [w, d, h] = fx.f.local_extents_m;
      const top = h > 1.0 ? 0.914 : h;
      const c = fxLocal(fx, w / 2, d * 0.5, top);
      let hood: PlumeSource["hood"] = null;
      for (const hf of hoods) {
        const b = hf.f.world_aabb_m!;
        const minX = b.min_m[0];
        const maxX = b.max_m[0];
        const minZ = -b.max_m[1];
        const maxZ = -b.min_m[1];
        if (c.x >= minX - 0.05 && c.x <= maxX + 0.05 && c.z >= minZ - 0.05 && c.z <= maxZ + 0.05) {
          // The filter bank sits at the back of the hood, high up.
          const back = fxLocal(hf, hf.f.local_extents_m[0] / 2, hf.f.local_extents_m[1] * 0.85, hf.f.local_extents_m[2] * 0.6);
          hood = { minX, maxX, minZ, maxZ, y: b.min_m[2], top: b.max_m[2], bx: back.x, bz: back.z };
          break;
        }
      }
      this.sources.push({ fx, x: c.x, y: c.y + 0.02, z: c.z, r: Math.min(w, d) * 0.35, hood });
    }
    this.acc = this.sources.map(() => 0);
    this.src.fill(-1);
    this.alpha.fill(0);
  }

  /**
   * Advance the particles. A point-source plume's centreline speed is w = 1.03 (g Q / (rho cp T))^(1/3) z^(-1/3)
   * (Heskestad) and it widens as b = 0.12 z; under a hood, what reaches the skirt is drawn up into the filters
   * except the uncaptured fraction, which rolls out past the skirt and drifts up to the ceiling.
   */
  tick(dt: number): void {
    const S = this.sources;
    const g = 9.81;
    const rhoCpT = 1.2 * 1005 * 300;
    // Particles per second per kW of convective heat.
    const rate = 9;
    for (let s = 0; s < S.length; s++) {
      const q = S[s]!.fx.plumeW / 1000;
      this.acc[s]! += q * rate * dt;
    }
    let free = 0;
    for (let i = 0; i < N_PLUME; i++) {
      if (this.src[i]! < 0) {
        // Spawn into this slot from whichever source has particles owed.
        for (let s = free; s < S.length; s++) {
          if (this.acc[s]! >= 1) {
            this.acc[s]! -= 1;
            const src = S[s]!;
            const a = Math.random() * Math.PI * 2;
            const rr = Math.sqrt(Math.random()) * src.r;
            this.pos[i * 3] = src.x + Math.cos(a) * rr;
            this.pos[i * 3 + 1] = src.y;
            this.pos[i * 3 + 2] = src.z + Math.sin(a) * rr;
            this.off[i * 2] = Math.cos(a) * (rr / Math.max(0.01, src.r));
            this.off[i * 2 + 1] = Math.sin(a) * (rr / Math.max(0.01, src.r));
            this.src[i] = s;
            this.age[i] = 0;
            this.life[i] = 4 + Math.random() * 2;
            this.esc[i] = Math.random() > this.capture ? 1 : 0;
            free = s;
            break;
          }
        }
        continue;
      }
      const src = S[this.src[i]!]!;
      const Q = Math.max(200, src.fx.plumeW);
      const z = Math.max(0.05, this.pos[i * 3 + 1]! - src.y) + 0.35;
      const w = 1.03 * Math.cbrt((g * Q) / rhoCpT) * Math.pow(z, -1 / 3);
      const b = src.r + 0.12 * (z - 0.35);
      let x = this.pos[i * 3]!;
      let y = this.pos[i * 3 + 1]!;
      let zz = this.pos[i * 3 + 2]!;
      const hd = src.hood;
      this.age[i]! += dt;
      let fade = Math.min(1, this.age[i]! / 0.4);
      if (hd && y > hd.y - 0.05 && !this.esc[i]) {
        // Captured: drawn to the filter bank and gone.
        const dx = hd.bx - x;
        const dy = hd.top - 0.15 - y;
        const dz = hd.bz - zz;
        const L = Math.hypot(dx, dy, dz) || 1;
        const v = Math.max(w, 1.2);
        x += (dx / L) * v * dt;
        y += (dy / L) * v * dt;
        zz += (dz / L) * v * dt;
        fade *= Math.min(1, L / 0.5);
        if (L < 0.12) this.src[i] = -1;
      } else if (hd && y > hd.y - 0.1 && this.esc[i]) {
        // Spilled: rolls out under the skirt to the nearest edge and rises outside it.
        const ex = Math.min(x - hd.minX, hd.maxX - x);
        const ez = Math.min(zz - hd.minZ, hd.maxZ - zz);
        if (ex > -0.05 && ez > -0.05) {
          if (ex < ez) x += (x - hd.minX < hd.maxX - x ? -1 : 1) * 0.6 * dt;
          else zz += (zz - hd.minZ < hd.maxZ - zz ? -1 : 1) * 0.6 * dt;
          y = Math.min(y, hd.y - 0.02);
        } else y += w * 0.5 * dt;
        fade *= 0.8;
      } else {
        y += w * dt;
        // Widen with height.
        x = src.x + this.off[i * 2]! * b + (Math.random() - 0.5) * 0.01;
        zz = src.z + this.off[i * 2 + 1]! * b + (Math.random() - 0.5) * 0.01;
      }
      if (y > this.ceiling - 0.05 || this.age[i]! > this.life[i]!) this.src[i] = -1;
      fade *= Math.min(1, (this.life[i]! - this.age[i]!) / 0.8);
      this.pos[i * 3] = x;
      this.pos[i * 3 + 1] = y;
      this.pos[i * 3 + 2] = zz;
      this.alpha[i] = this.src[i]! < 0 ? 0 : 0.085 * fade * (this.esc[i] ? 1.6 : 1);
      this.size[i] = 0.1 + 0.14 * Math.min(1, z / 1.2);
    }
    const geo = this.points.geometry;
    geo.getAttribute("position").needsUpdate = true;
    geo.getAttribute("aAlpha").needsUpdate = true;
    geo.getAttribute("aSize").needsUpdate = true;
  }
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Aisles: the walking clearance field, the site's corridors, door swings                                        */
/* ------------------------------------------------------------------------------------------------------------ */

function buildAisles(K: Kitchen, rWalk: number): THREE.Group {
  const g = new THREE.Group();
  g.name = "aisles";
  // The clearance field as a floor texture: green where a walking body fits with room, amber where it just
  // fits, red where it cannot pass.
  const { nx, nz, clr } = K.grid;
  const data = new Uint8Array(nx * nz * 4);
  for (let j = 0; j < nz; j++)
    for (let i = 0; i < nx; i++) {
      const c = clr[j * nx + i]!;
      const o = ((nz - 1 - j) * nx + i) * 4;
      if (c < 0) {
        data[o + 3] = 0;
        continue;
      }
      const fit = c - rWalk;
      const col = fit < 0 ? [239, 68, 68] : fit < 0.15 ? [245, 158, 11] : [52, 211, 153];
      data[o] = col[0]!;
      data[o + 1] = col[1]!;
      data[o + 2] = col[2]!;
      data[o + 3] = fit < 0 ? 150 : fit < 0.15 ? 175 : Math.max(70, 150 - fit * 110);
    }
  const tex = new THREE.DataTexture(data, nx, nz, THREE.RGBAFormat);
  tex.magFilter = THREE.LinearFilter;
  tex.minFilter = THREE.LinearFilter;
  tex.needsUpdate = true;
  const field = new THREE.Mesh(
    new THREE.PlaneGeometry(nx * K.grid.cell, nz * K.grid.cell).rotateX(-Math.PI / 2).translate(K.grid.x0 + (nx * K.grid.cell) / 2, 0.006, K.grid.z0 + (nz * K.grid.cell) / 2),
    new THREE.MeshBasicMaterial({ map: tex, transparent: true, depthWrite: false, polygonOffset: true, polygonOffsetFactor: -2, polygonOffsetUnits: -2 }),
  );
  field.renderOrder = 3;
  g.add(field);
  // The site's walk corridors, outlined.
  const lineMat = new THREE.LineBasicMaterial({ color: 0x6ee7b7, transparent: true, opacity: 0.9 });
  for (const c of K.meta?.walk_mesh_corridors ?? []) {
    const [x0, y0] = c.min_m;
    const [x1, y1] = c.max_m;
    const pts = [
      [x0!, -y0!],
      [x1!, -y0!],
      [x1!, -y1!],
      [x0!, -y1!],
      [x0!, -y0!],
    ].map(([x, z]) => new THREE.Vector3(x, 0.012, z));
    const l = new THREE.Line(new THREE.BufferGeometry().setFromPoints(pts), lineMat);
    l.renderOrder = 4;
    g.add(l);
  }
  // Door swings and other swept volumes, as amber footprints.
  const sweptMat = new THREE.MeshBasicMaterial({ color: 0xf59e0b, transparent: true, opacity: 0.28, depthWrite: false });
  const edgeMat = new THREE.LineBasicMaterial({ color: 0xfbbf24 });
  for (const fx of K.fixtures)
    for (const v of fx.f.swept_volumes ?? []) {
      const a = v.aabb_world_m;
      if (a.min_m[2] > 1.9) continue;
      const w = a.max_m[0] - a.min_m[0];
      const d = a.max_m[1] - a.min_m[1];
      const plane = new THREE.Mesh(new THREE.PlaneGeometry(w, d).rotateX(-Math.PI / 2), sweptMat);
      plane.position.set((a.min_m[0] + a.max_m[0]) / 2, 0.01, -(a.min_m[1] + a.max_m[1]) / 2);
      plane.renderOrder = 4;
      const e = new THREE.LineSegments(new THREE.EdgesGeometry(new THREE.PlaneGeometry(w, d).rotateX(-Math.PI / 2)), edgeMat);
      e.position.copy(plane.position).setY(0.013);
      g.add(plane, e);
    }
  return g;
}

/* ------------------------------------------------------------------------------------------------------------ */
/* The crew: stations from the manifest's sockets, people walking between them through the crowd planner         */
/* ------------------------------------------------------------------------------------------------------------ */

/** Sampling step of the people's committed tracks, seconds. */
const CROWD_DT = 0.1;
/** How long before a person's work ends their next leg is planned, seconds. */
const LOOKAHEAD = 2.5;
/** Planning time allowed each frame, ms. */
const PLAN_BUDGET_MS = 4;

type ToolKind = "skillet" | "knife" | "ladle" | "brush" | "spoon";
/** Work at a station; `flip`: steaks turned on a charbroiler with tongs; `sala`: a salamander's grid pulled out to look. */
type TaskKind = "saute" | "stir" | "chop" | "plate" | "scrub" | "fry" | "flip" | "sala";
const WORK_OF: Record<TaskKind, WorkKind> = { saute: "saute", stir: "stir", chop: "chop", plate: "plate", scrub: "stir", fry: "stand", flip: "stand", sala: "stand" };
const TOOL_OF: Record<TaskKind, ToolKind | null> = { saute: "skillet", stir: "ladle", chop: "knife", plate: "spoon", scrub: "brush", fry: null, flip: null, sala: null };
/** How far each kind of work may reach from where the body stands to its target, metres at p50 scale. */
const REACH: Record<TaskKind, number> = { saute: 0.66, stir: 0.72, chop: 0.6, plate: 0.6, scrub: 0.72, fry: 0.6, flip: 0.7, sala: 0.85 };
/** How long a person stays at each kind of work, seconds (min, max). */
const DWELL: Record<TaskKind, [number, number]> = { saute: [9, 16], stir: [7, 12], chop: [9, 16], plate: [6, 10], scrub: [8, 14], fry: [7, 7], flip: [10.5, 10.5], sala: [8.5, 8.5] };
/** When a tool comes into the hand and goes back down, seconds into the work, and how long the hand takes between. */
const TOOL_IN = 0.95;
const TOOL_BLEND = 0.7;

interface Basket {
  down: boolean;
  /** When it last went down or up, seconds. */
  since: number;
}

/** A tool's resting frame, world: its grip point, the handle's axis towards the working end, and its top. */
interface Rest {
  p: THREE.Vector3;
  z: THREE.Vector3;
  y: THREE.Vector3;
}

export interface Station {
  key: string;
  kind: TaskKind;
  fx: Fx;
  /** The work's target point (pan base, pot rim, board top, plate) and the normal from the fixture to the person. */
  at: THREE.Vector3;
  nx: number;
  nz: number;
  /** Where the floor meets the fixture in front of the target. */
  faceX: number;
  faceZ: number;
  holder: number;
  /** The station's own tool and where it rests when nobody holds it. */
  tool: THREE.Object3D | null;
  rest: Rest | null;
  /** A pickup with the hand turned over (a knife lying flat). */
  turnOver?: boolean;
  baskets?: Basket[];
  /** Portioning from a hotel pan rather than plating: no plate is set out. */
  noPlate: boolean;
  /** Which basket the person working here is lifting or dropping, and which way. */
  fryJob?: { k: number; lift: boolean };
  /** A charbroiler's steaks on the grates (each turned over in turn) and the tongs, with where they rest. */
  steaks?: { obj: THREE.Object3D; at: THREE.Vector3; up: number; done: number }[];
  tongs?: { obj: THREE.Object3D; model: Model; rest: Rest };
}

interface Agent {
  id: number;
  st: Staff;
  ghost: Staff;
  cohort: Cohort;
  /** Walking disc, path clearance and the distance kept from a counter face when working. */
  rWalk: number;
  rPath: number;
  dFace: number;
  shapes: Footprint[];
  /** How far ahead of the root a rack's near face is when pushed, and where the hands grip its uprights. */
  rackFace: number;
  gripY: number;
  /** How far the person leans in to reach the work (a salamander's pull over a range). */
  lean: number;
  speed: number;
  routine: Routine;
  steps: RoutineStep[];
  start: { x: number; z: number; yaw: number };
  t0: number;
  /** The committed end of the current routine (its track's last step), seconds. */
  endT: number;
  /** The station worked at the end of the current leg, or null when the leg is a run (pans, a rack). */
  station: Station | null;
  run: Run | null;
  workDur: number;
  planned: { steps: RoutineStep[]; start: { x: number; z: number; yaw: number }; t0: number; station: Station | null; run: Run | null; workDur: number; path: XZ[] } | null;
  fails: number;
  weights: Map<Station, number>;
  path: XZ[] | null;
  line: THREE.Line | null;
  /** How far the hand has gone from the tool's rest to the work (0..1), and whether it holds the tool. */
  toolK: number;
  holding: boolean;
  accent: THREE.Color;
}

const COHORT_SPEED: Record<Cohort, number> = { p05: 1.02, p50: 1.08, p95: 1.12 };

function rngSeed(seed: number): () => number {
  let s = seed >>> 0 || 1;
  return () => {
    s ^= s << 13;
    s ^= s >>> 17;
    s ^= s << 5;
    return (s >>> 0) / 4294967296;
  };
}

const smoothstep = (e0: number, e1: number, x: number) => {
  const q = Math.min(1, Math.max(0, (x - e0) / (e1 - e0)));
  return q * q * (3 - 2 * q);
};
/** A body-frame direction (x to the left, z forward) in world for a person facing yaw. */
function bodyDir(yaw: number, x: number, y: number, z: number, out = new THREE.Vector3()): THREE.Vector3 {
  const c = Math.cos(yaw);
  const s = Math.sin(yaw);
  return out.set(c * x + s * z, y, -s * x + c * z);
}
function restMatrix(r: Rest, out = new THREE.Matrix4()): THREE.Matrix4 {
  const z = r.z.clone().normalize();
  const y = r.y.clone().addScaledVector(z, -r.y.dot(z)).normalize();
  const x = new THREE.Vector3().crossVectors(y, z);
  return out.makeBasis(x, y, z).setPosition(r.p);
}

/** The knife's grip frame in its model's frame (fixture x across the blade to the spine, z back to the tip at -d). */
function knifeHolder(ext: { w: number; d: number; h: number }): THREE.Matrix4 {
  const xh = ext.w - 0.0012 - TOOLS.knife.spineY;
  const basis = new THREE.Matrix4().makeBasis(new THREE.Vector3(0, 1, 0), new THREE.Vector3(1, 0, 0), new THREE.Vector3(0, 0, -1));
  return basis.multiply(new THREE.Matrix4().makeTranslation(-xh, -ext.h / 2, TOOLS.knife.handleBack));
}

/** The tools people hold, built in their grip frame (handle along +z towards the working end, +y its top). */
function makeTool(kind: ToolKind, knife?: THREE.Object3D, knifeExt?: { w: number; d: number; h: number }): THREE.Object3D {
  const p = palette();
  const grip = new THREE.Group();
  grip.name = `tool-${kind}`;
  if (kind === "skillet") {
    // An 11 in tri-ply skillet whose handle rises at the angle the saute pose holds it (TOOLS.skillet.rise): the
    // model has its base centre at (w/2, 0, -d + D/2) and the handle running to +z, so turn it half round and tip
    // it by the rise to lay the handle along the grip axis, the grip point gripFromRim up the handle from the rim.
    const K = TOOLS.skillet;
    const D = 2 * K.rim;
    const L = 0.1 + 0.27 * D;
    const rimY = 0.026 + 0.075 * D;
    const yRoot = rimY - 0.011;
    const total = yRoot + Math.tan(K.rise) * L + 0.0075;
    const ext = { w: D + 0.001, d: D + L + 0.001, h: total };
    const m = buildFixture("NSF-SKILLET-10IN", ext, 1, true).group;
    const gy = yRoot + K.gripFromRim * Math.sin(K.rise);
    const gz = -ext.d + D / 2 + K.rim + K.gripFromRim * Math.cos(K.rise);
    const holder = new THREE.Group();
    holder.quaternion.setFromEuler(new THREE.Euler(-K.rise, Math.PI, 0, "XYZ"));
    m.position.set(-ext.w / 2, -gy, -gz);
    holder.add(m);
    grip.add(holder);
  } else if (kind === "knife") {
    const ext = knifeExt ?? { w: 0.052, d: TOOLS.knife.tipZ + TOOLS.knife.handleBack, h: 0.019 };
    const m = knife ?? buildFixture("NSF-CHEF-KNIFE-10IN", ext, 1, true).group;
    const holder = new THREE.Group();
    holder.matrixAutoUpdate = false;
    holder.matrix.copy(knifeHolder(ext));
    m.position.set(0, 0, 0);
    m.quaternion.identity();
    m.updateMatrix();
    holder.add(m);
    grip.add(holder);
  } else if (kind === "ladle" || kind === "brush") {
    const K = TOOLS.ladle;
    const len = K.reach + K.handleBack - 0.04;
    const handle = new THREE.Mesh(new THREE.CylinderGeometry(0.0065, 0.0065, len, 10).rotateX(Math.PI / 2).translate(0, 0, (K.reach - 0.04 - K.handleBack) / 2), kind === "ladle" ? p.ssBright : p.poly);
    grip.add(handle);
    if (kind === "ladle") {
      const bowl = new THREE.Mesh(new THREE.SphereGeometry(K.bowlR, 18, 10, 0, Math.PI * 2, Math.PI / 2, Math.PI / 2), p.ssBright);
      bowl.position.set(0, 0, K.reach);
      bowl.rotation.x = -Math.PI / 2;
      grip.add(bowl);
    } else {
      const head = new THREE.Mesh(new THREE.BoxGeometry(0.05, 0.03, 0.08), p.hdpe);
      head.position.set(0, 0, K.reach - 0.03);
      const bristle = new THREE.Mesh(new THREE.BoxGeometry(0.046, 0.026, 0.05), new THREE.MeshStandardMaterial({ color: 0xc9a24a, roughness: 0.95 }));
      bristle.position.set(0, 0, K.reach + 0.035);
      grip.add(head, bristle);
    }
  } else {
    const K = TOOLS.spoon;
    const len = K.reach + K.handleBack - 0.03;
    const handle = new THREE.Mesh(new THREE.CylinderGeometry(0.0045, 0.006, len, 8).rotateX(Math.PI / 2).translate(0, 0, (K.reach - 0.03 - K.handleBack) / 2), p.ssBright);
    const bowl = new THREE.Mesh(new THREE.SphereGeometry(0.024, 14, 8, 0, Math.PI * 2, Math.PI / 2, Math.PI / 2), p.ssBright);
    bowl.scale.set(1, 1.5, 0.45);
    bowl.rotation.x = Math.PI / 2;
    bowl.position.set(0, 0, K.reach);
    grip.add(handle, bowl);
  }
  grip.traverse((o) => {
    if ((o as THREE.Mesh).isMesh) o.castShadow = true;
  });
  return grip;
}

/** A white plate with food on it, for the pass and the makeline. */
function makePlate(): THREE.Object3D {
  const g = new THREE.Group();
  const china = new THREE.MeshPhysicalMaterial({ color: 0xf4f1ea, roughness: 0.25, clearcoat: 0.6 });
  const prof = [[0, 0], [0.075, 0], [0.085, 0.004], [0.125, 0.016], [0.13, 0.019], [0.128, 0.021], [0.08, 0.009], [0, 0.009]].map(([x, y]) => new THREE.Vector2(x, y));
  const plate = new THREE.Mesh(new THREE.LatheGeometry(prof, 40), china);
  plate.castShadow = plate.receiveShadow = true;
  const food = new THREE.Mesh(new THREE.SphereGeometry(0.05, 18, 10, 0, Math.PI * 2, 0, Math.PI / 2), new THREE.MeshStandardMaterial({ color: 0x9a5a2c, roughness: 0.7 }));
  food.scale.set(1.2, 0.45, 0.9);
  food.position.y = 0.009;
  const greens = new THREE.Mesh(new THREE.SphereGeometry(0.025, 12, 8, 0, Math.PI * 2, 0, Math.PI / 2), new THREE.MeshStandardMaterial({ color: 0x4b7a2e, roughness: 0.8 }));
  greens.position.set(0.055, 0.009, 0.02);
  greens.scale.set(1, 0.6, 1);
  g.add(plate, food, greens);
  return g;
}

/** A steel ninth-pan of water the plating spoon stands in. */
function makeBain(): THREE.Object3D {
  const m = new THREE.Mesh(new THREE.CylinderGeometry(0.055, 0.05, 0.1, 24, 1, true), palette().ssBright);
  m.position.y = 0.05;
  const water = new THREE.Mesh(new THREE.CircleGeometry(0.052, 24).rotateX(-Math.PI / 2), new THREE.MeshStandardMaterial({ color: 0x5a6a70, roughness: 0.08, metalness: 0.2, transparent: true, opacity: 0.7 }));
  water.position.y = 0.075;
  const g = new THREE.Group();
  g.add(m, water);
  m.castShadow = true;
  return g;
}

/** Footprint ids (the crowd's shapes): standing, walking; with a hotel pan held; behind a rack. */
const SH_STAND = 0;
const SH_WALK = 1;
const SH_PAN = 2;
const SH_PAN_WALK = 3;
const SH_RACK = 4;
const SH_RACK_WALK = 5;
const SH_TURN = 6;
const SH_TURN_WALK = 7;
/** A speed rack's footprint (NSF-SPEEDRACK-20PAN); a loaded rack's pushing speed, m/s. */
const RACK_DIMS = { w: 0.521, d: 0.66 };
const RACK_SPEED = 0.62;

/**
 * Footprints, body frame (x to the left, z forward): standing (three discs across the shoulders), walking (one); each
 * also with a full hotel pan held at the waist (four discs over the pan and the hands at its ends; `high`, as the pan
 * passes over counters), and with a speed rack ahead in the hands (four discs over its 0.52 m x 0.66 m).
 */
function footprints(st: Staff): { shapes: Footprint[]; rWalk: number; rackFace: number; gripY: number } {
  const sb = st.dims.shoulder_breadth_m;
  const cd = st.dims.chest_depth_m;
  const rr = cd / 2 + 0.02;
  const half = sb / 2;
  const rWalk = half + 0.015;
  const stand: Footprint = [
    { x: -(half - rr), z: 0, r: rr },
    { x: 0, z: 0, r: rr },
    { x: half - rr, z: 0, r: rr },
  ];
  const walk: Footprint = [{ x: 0, z: 0, r: rWalk }];
  // The pan's near wall a couple of centimetres off the apron at the belly (people.ts carryFrame).
  const pz = cd * 0.55 + 0.012 + HOTEL_PAN.width / 2 + 0.025;
  const pan: Footprint = [-1, 1].flatMap((sx) => [-1, 1].map((sz) => ({ x: sx * 0.15, z: pz + sz * 0.081, r: 0.165, high: true })));
  // The rack's near uprights in the hands a forearm and a half ahead, gripped a hand below the shoulders: elbows
  // bent, the push through the arms; the rack's casters well clear of the feet.
  const s = st.B.s;
  const rackFace = 0.45 * s;
  const gripY = st.dims.shoulder_height_m - 0.25 * s;
  const hw = RACK_DIMS.w / 2;
  const hd = RACK_DIMS.d / 2;
  const rack: Footprint = [-1, 1].flatMap((sx) => [-1, 1].map((sz) => ({ x: (sx * hw) / 2, z: rackFace + hd + (sz * hd) / 2, r: Math.hypot(hw / 2, hd / 2) + 0.005 })));
  // A pan held turned (to go onto a combi's runners lengthwise): long side ahead, its near end clear of the apron.
  const tz = panHoldZ(st);
  const turned: Footprint = [-1, 1].flatMap((sx) => [-1, 1].map((sz) => ({ x: sx * 0.081, z: tz + sz * 0.15, r: 0.165, high: true })));
  return { shapes: [stand, walk, [...stand, ...pan], [...walk, ...pan], [...stand, ...rack], [...walk, ...rack], [...stand, ...turned], [...walk, ...turned]], rWalk, rackFace, gripY };
}


/** Where the floor meets the fixture in front of a target: walk from it along the normal until clear. */
function faceOf(grid: NavGrid, at: THREE.Vector3, nx: number, nz: number): [number, number] | null {
  let x = at.x;
  let z = at.z;
  if (grid.clearance(x, z) >= 0) return null;
  for (let k = 0; k < 200; k++) {
    x += nx * 0.01;
    z += nz * 0.01;
    if (grid.clearance(x, z) >= 0) return [x, z];
  }
  return null;
}

/** The fixture's outward front normal on the floor. */
function frontNormal(fx: Fx): [number, number] {
  const y = fx.f.pose.yaw_rad;
  return [Math.sin(y), Math.cos(y)];
}

function makeStations(K: Kitchen, extras: THREE.Group): Station[] {
  const out: Station[] = [];
  const add = (fx: Fx, kind: TaskKind, at: THREE.Vector3, key: string, n?: [number, number], noPlate = false) => {
    const [nx, nz] = n ?? frontNormal(fx);
    const face = faceOf(K.grid, at, nx, nz);
    if (!face) return null;
    const deep = (at.x - face[0]) * -nx + (at.z - face[1]) * -nz;
    if (kind === "saute" && deep > 0.34) at = new THREE.Vector3(face[0] - nx * 0.3, at.y, face[1] - nz * 0.3);
    const s: Station = { key: `${fx.f.fixture_id}:${key}`, kind, fx, at: at.clone(), nx, nz, faceX: face[0], faceZ: face[1], holder: -1, tool: null, rest: null, noPlate };
    out.push(s);
    return s;
  };
  const smalls = K.fixtures.filter((f) => SMALLWARE.test(f.id) && f.f.world_aabb_m);
  const centre = (f: Fx) => {
    const b = f.f.world_aabb_m!;
    return new THREE.Vector3((b.min_m[0] + b.max_m[0]) / 2, (b.min_m[2] + b.max_m[2]) / 2, -(b.min_m[1] + b.max_m[1]) / 2);
  };
  /** Whether a footprint (three xz box) on a counter is free of the smallwares standing on it. */
  const freeOnTop = (x0: number, x1: number, z0: number, z1: number) =>
    !smalls.some((f) => {
      const b = f.f.world_aabb_m!;
      return b.min_m[0] < x1 && b.max_m[0] > x0 && -b.max_m[1] < z1 && -b.min_m[1] > z0;
    });
  const onTop = (fx: Fx, f: Fx) => {
    const c = centre(f).applyMatrix4(new THREE.Matrix4().copy(fx.model.group.matrixWorld).invert());
    const [w, d] = fx.f.local_extents_m;
    return c.x > 0 && c.x < w && c.z < 0 && c.z > -d;
  };
  for (const fx of K.fixtures) {
    const id = fx.id;
    const [w, d, h] = fx.f.local_extents_m;
    const seat = (sid: string) => {
      const s = fxSocket(fx, sid);
      return s ? fxLocal(fx, s.local_pos_m[0], s.local_pos_m[1], s.local_pos_m[2]) : null;
    };
    if (/RANGE-GAS/.test(id)) {
      const at = seat("pan_seat_front_center");
      if (at) add(fx, "saute", at, "pan");
    } else if (/GRIDDLE/.test(id)) {
      const at = seat("pan_seat_platen");
      if (at) add(fx, "saute", at, "platen");
    } else if (/INDUCTION/.test(id)) {
      const at = seat("pan_seat_coil_0");
      if (at) add(fx, "saute", at, "coil");
    } else if (/FRYER/.test(id)) {
      const s = add(fx, "fry", fxLocal(fx, w / 2, 0.3, 0.9), "baskets");
      if (s) s.baskets = [{ down: true, since: -20 }, { down: false, since: -5 }];
    } else if (/STOCKPOT/.test(id)) {
      add(fx, "stir", fxLocal(fx, w / 2, d / 2, h), "pot");
    } else if (/KETTLE/.test(id)) {
      add(fx, "stir", fxLocal(fx, w / 2, 0.36, Math.min(h, 1.0)), "kettle");
    } else if (/TILT-SKILLET/.test(id)) {
      add(fx, "stir", fxLocal(fx, w / 2, 0.38, 0.95), "pan");
    } else if (/SINK-3COMP|VEG-PREP-SINK/.test(id)) {
      add(fx, "scrub", fxLocal(fx, /3COMP/.test(id) ? w / 2 : 0.34, 0.36, 0.914), "basin");
    } else if (/PREPTABLE/.test(id)) {
      // A board on the table: the manifest's own if one lies on it (and its knife), else one of ours where the top
      // is clear.
      const board = smalls.find((f) => /CUTTING-BOARD/.test(f.id) && onTop(fx, f));
      const boards: THREE.Vector3[] = [];
      if (board) {
        const bs = fxSocket(board, "pan_seat_board");
        const s = bs && add(fx, "chop", fxLocal(board, bs.local_pos_m[0], bs.local_pos_m[1], bs.local_pos_m[2]), "board");
        const knife = smalls.find((f) => /CHEF-KNIFE/.test(f.id) && onTop(fx, f) && centre(f).distanceTo(centre(board)) < 0.6);
        if (s) boards.push(s.at);
        if (s && knife) {
          const [kw, kd, kh] = knife.f.local_extents_m;
          const ext = { w: kw, d: kd, h: kh };
          const M = knife.model.group.matrixWorld.clone().multiply(new THREE.Matrix4().copy(knifeHolder(ext)).invert());
          const p = new THREE.Vector3();
          const q = new THREE.Quaternion();
          M.decompose(p, q, new THREE.Vector3());
          s.rest = { p, z: new THREE.Vector3(0, 0, 1).applyQuaternion(q), y: new THREE.Vector3(0, 1, 0).applyQuaternion(q) };
          s.tool = makeTool("knife", knife.model.group, ext);
          s.turnOver = true;
        }
      }
      const pans = smalls.filter((f) => /HOTEL-PAN/.test(f.id) && onTop(fx, f));
      if (pans.length) {
        // Portion from the pans, a cook to every two.
        for (let k = 0; k + 1 < pans.length; k += 2) {
          const c = centre(pans[k]!).add(centre(pans[k + 1]!)).multiplyScalar(0.5);
          c.y = pans[k]!.f.world_aabb_m!.max_m[2];
          add(fx, "plate", c, `pans${k}`, undefined, true);
        }
        continue;
      }
      // A table under a heat lamp is the pass: it is for plates, not boards.
      const lamps = K.fixtures.filter((f) => /HEATLAMP/.test(f.id) && f.f.world_aabb_m && onTop(fx, f));
      if (lamps.length) continue;
      const bw = 0.61;
      const bd = 0.457;
      // A table keeps some of its top as a staging area (for the pans brought out): two boards only on a long one
      // that has none of its own.
      const want = w > 2.2 && !board ? 2 : 1;
      const placedX: number[] = [];
      for (const ox of w > 1.8 ? [0.25, 0.75, 0.5, 0.15, 0.85] : [0.5, 0.3, 0.7, 0.2, 0.8]) {
        if (placedX.length + (board ? 1 : 0) >= want) break;
        const lx = Math.max(bw / 2 + 0.03, Math.min(w - bw / 2 - 0.03, w * ox));
        if (placedX.some((x) => Math.abs(x - lx) < 0.95)) continue;
        const c0 = fxLocal(fx, lx - bw / 2, 0.06, 0.92);
        const c1 = fxLocal(fx, lx + bw / 2, 0.06 + bd, 0.92);
        if (!freeOnTop(Math.min(c0.x, c1.x), Math.max(c0.x, c1.x), Math.min(c0.z, c1.z), Math.max(c0.z, c1.z))) continue;
        const m = buildFixture("NSF-CUTTING-BOARD-HDPE-24X18", { w: bw, d: bd, h: 0.019 }, 1, true);
        // The board's front edge 6 cm behind the table's front, its top at 0.933 m.
        m.group.position.set(lx - bw / 2, h > 0.95 ? h : 0.914, -0.06);
        fx.model.group.add(m.group);
        m.group.updateMatrixWorld(true);
        extras.userData.boards = ((extras.userData.boards as number | undefined) ?? 0) + 1;
        const at = fxLocal(fx, lx, 0.06 + bd / 2, (h > 0.95 ? h : 0.914) + 0.019);
        // Not beside the manifest's own board, where its cook stands.
        if (boards.some((b) => b.distanceTo(at) < 0.95)) {
          fx.model.group.remove(m.group);
          continue;
        }
        placedX.push(lx);
        add(fx, "chop", at, `board${placedX.length}`);
      }
    } else if (/CHARBROILER/.test(id)) {
      // Steaks on the front half of the grates (their top 0.32 m up the unit), within an easy reach with tongs.
      add(fx, "flip", fxLocal(fx, w / 2, 0.27, 0.3235), "grates");
    } else if (/SALAMANDER/.test(id)) {
      // The grid's pull, its carriage at the middle of its travel (cooking.ts salamander).
      add(fx, "sala", fxThree(fx, w / 2, 0.176, -0.024), "grid");
    } else if (/MAKELINE/.test(id)) {
      add(fx, "plate", fxLocal(fx, w / 2, 0.17, 0.936), "board");
    } else if (/HEATLAMP/.test(id)) {
      const at = seat("pan_seat_pass");
      if (!at) continue;
      // The pass under the lamp: worked from whichever side has floor to stand on, nearer.
      const [nx, nz] = frontNormal(fx);
      const a = faceOf(K.grid, at, nx, nz);
      const b = faceOf(K.grid, at, -nx, -nz);
      const da = a ? Math.hypot(a[0] - at.x, a[1] - at.z) : Infinity;
      const db = b ? Math.hypot(b[0] - at.x, b[1] - at.z) : Infinity;
      if (da === Infinity && db === Infinity) continue;
      add(fx, "plate", at, "pass", da <= db ? [nx, nz] : [-nx, -nz]);
    }
  }
  return out;
}

/** Where each station's tool rests (when it is not the manifest's own, already placed). */
function restFor(s: Station): Rest | null {
  const yaw = Math.atan2(-s.nx, -s.nz);
  const n = new THREE.Vector3(s.nx, 0, s.nz);
  const Y = new THREE.Vector3(0, 1, 0);
  switch (TOOL_OF[s.kind]) {
    case "skillet": {
      // On the burner as the saute pose sets it down: the base on the grate, the handle back to the right hand.
      const K = TOOLS.skillet;
      const hh = bodyDir(yaw, 0.38, 0, 0.92).normalize();
      const hd = hh.clone().multiplyScalar(Math.cos(K.rise)).addScaledVector(Y, -Math.sin(K.rise));
      const rimNear = s.at.clone().addScaledVector(Y, K.depth).addScaledVector(hh, -K.rim);
      const p = rimNear.addScaledVector(hd, -K.gripFromRim);
      return { p, z: hd, y: Y.clone() };
    }
    case "ladle":
    case "brush": {
      // Left in the pot (or the sink), the bowl on the bottom at the far side, the handle leaning on the near rim.
      const K = TOOLS.ladle;
      const bowl = s.at.clone().addScaledVector(n, -0.07).addScaledVector(Y, -0.24);
      const rim = s.at.clone().addScaledVector(n, 0.11).addScaledVector(Y, 0.012);
      const dir = rim.sub(bowl).normalize();
      const p = bowl.clone().addScaledVector(dir, K.reach);
      return { p, z: dir.clone().negate(), y: n.clone().addScaledVector(Y, 0.3) };
    }
    case "spoon": {
      // Standing in its bain of water to the right of the plate, bowl down.
      const right = bodyDir(yaw, -1, 0, 0);
      const base = s.at.clone().addScaledVector(right, 0.2).addScaledVector(n, 0.04);
      return { p: base.addScaledVector(Y, 0.012 + TOOLS.spoon.reach), z: Y.clone().negate(), y: n.clone().addScaledVector(Y, 0.0001) };
    }
    case "knife": {
      // Lying flat on the board's right front corner, the edge away, the tip to the left.
      const left = bodyDir(yaw, 1, 0, 0);
      const p = s.at.clone().addScaledVector(left, -0.2).addScaledVector(n, 0.15).addScaledVector(Y, 0.026);
      return { p, z: left, y: n.clone() };
    }
  }
  return null;
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Runs: hotel pans fetched from the cold stores to the line and taken back, speed racks pushed                   */
/* ------------------------------------------------------------------------------------------------------------ */

/**
 * A hinged door people open: a reach-in's, a walk-in's. Its leaf turns about the hinge from `e0` (shut, hinge to
 * free edge) towards `e1` (the outward normal of the shut door) by up to `max` radians.
 */
interface Door {
  fx: Fx;
  /** The fixture hook's name for it: set(name, f) opens it to f in 0..1 ("open": the model's own open(f)). */
  name: string;
  hinge: THREE.Vector3;
  e0: THREE.Vector3;
  e1: THREE.Vector3;
  len: number;
  /** How far it opens here (radians: short of a neighbour in its swing), and the model's own stop (its hook's 1). */
  max: number;
  hookMax: number;
  /** The pull's grip: along the leaf from the hinge, out from its face; its height (world). */
  hAlong: number;
  hOut: number;
  hY: number;
  /** How far the leaf's inner face (towards the cabinet when shut) is from the hinge line. */
  back: number;
  /**
   * Where along the leaf (from the hinge) the pull begins to stand out, and how far the leaf's outer face stands out
   * from the hinge line short of it. Unset: the pull's stand-out is taken along the whole leaf.
   */
  pullFrom?: number;
  face?: number;
  /** How far open the hand takes it before letting it swing on to its stop. */
  pullTo: number;
  /** Seconds the hand pulls; the leaf then swings on, slowing, to its stop (`swing` seconds); closing takes `close`. */
  pull: number;
  swing: number;
  close: number;
  /** The floor its leaf sweeps (world). */
  rect: Rect;
}

const _dv1 = new THREE.Vector3();
const _dv2 = new THREE.Vector3();
/** Swing a door's model to f of its opening here (0 shut, 1 at its `max`). */
function setDoor(d: Door, f: number): void {
  const v = (f * d.max) / d.hookMax;
  if (d.name === "open") d.fx.model.hooks.open?.(v);
  else d.fx.model.hooks.set?.(d.name, v);
}
function leafDir(d: Door, f: number, out = new THREE.Vector3()): THREE.Vector3 {
  const th = f * d.max;
  return out.copy(d.e0).multiplyScalar(Math.cos(th)).addScaledVector(d.e1, Math.sin(th));
}
function leafNormal(d: Door, f: number, out = new THREE.Vector3()): THREE.Vector3 {
  const th = f * d.max;
  return out.copy(d.e0).multiplyScalar(-Math.sin(th)).addScaledVector(d.e1, Math.cos(th));
}
/** How far the leaf's pull side stands out from the hinge line, `u` along the leaf from the hinge. */
function outAt(d: Door, u: number): number {
  return d.pullFrom !== undefined && u < d.pullFrom ? (d.face ?? 0) : d.hOut;
}
/** The pull's grip point (world) with the door open to f. */
function handleAt(d: Door, f: number, out = new THREE.Vector3()): THREE.Vector3 {
  leafDir(d, f, _dv1);
  leafNormal(d, f, _dv2);
  return out.copy(d.hinge).addScaledVector(_dv1, d.hAlong).addScaledVector(_dv2, d.hOut).setY(d.hY);
}
/** How much of the leaf's free edge a hand takes hold of, as a bar of this radius. */
const EDGE_R = 0.02;
/** Where a hand takes the leaf by its free edge (world), at the pull's height, with the door open to f. */
function edgeAt(d: Door, f: number, out = new THREE.Vector3()): THREE.Vector3 {
  leafDir(d, f, _dv1);
  leafNormal(d, f, _dv2);
  return out.copy(d.hinge).addScaledVector(_dv1, d.len - EDGE_R).addScaledVector(_dv2, -d.back / 2).setY(d.hY);
}
/**
 * How far a floor point is from the leaf (its pull included) with the door open to f. The first few centimetres at
 * the hinge hardly move and are the cabinet's own frame (kept clear as a fixture), so the leaf is taken from there.
 */
function leafGap(d: Door, f: number, x: number, z: number): number {
  leafDir(d, f, _dv1);
  const px = x - d.hinge.x;
  const pz = z - d.hinge.z;
  const u = Math.max(0.06, Math.min(d.len, px * _dv1.x + pz * _dv1.z));
  const ex = px - _dv1.x * u;
  const ez = pz - _dv1.z * u;
  leafNormal(d, f, _dv2);
  // On the pull's side the pull stands out from the face; on the other, the leaf's own face.
  return Math.hypot(ex, ez) - (ex * _dv2.x + ez * _dv2.z > 0 ? outAt(d, u) + 0.02 : d.back);
}
/** Whether the leaf (its pull side included) open to f meets an obstacle or a wall. */
function leafHits(d: Door, f: number, rects: readonly Rect[], b: Rect): boolean {
  leafDir(d, f, _dv1);
  leafNormal(d, f, _dv2);
  for (let k = 0; k <= 10; k++) {
    const u = (d.len * k) / 10;
    for (const v of [-d.back, outAt(d, u) + 0.01]) {
      const x = d.hinge.x + _dv1.x * u + _dv2.x * v;
      const z = d.hinge.z + _dv1.z * u + _dv2.z * v;
      if (x < b.minX + 0.01 || x > b.maxX - 0.01 || z < b.minZ + 0.01 || z > b.maxZ - 0.01) return true;
      for (const r of rects) if (rectDistance(r, x, z) < 0.012) return true;
    }
  }
  return false;
}
/**
 * A door: it opens as far as its model lets it, or until its leaf would meet a neighbour standing proud of the
 * cabinet's front or a wall (the obstacles inside the fixture's own footprint are its own, and do not count).
 */
function makeDoor(fx: Fx, o: Omit<Door, "fx" | "rect" | "swing" | "hookMax"> & { swing?: number }, rects: readonly Rect[], bounds: Rect): Door {
  // The free swing continues the pull's speed and slows to the stop: its time follows from the pull's.
  const swing = o.swing ?? (o.pull * (1 - o.pullTo)) / o.pullTo;
  const d: Door = { ...o, fx, swing, hookMax: o.max, rect: { minX: Infinity, minZ: Infinity, maxX: -Infinity, maxZ: -Infinity } };
  const own = fx.rect!;
  const others = rects.filter((q) => !(q.minX >= own.minX - 1e-3 && q.maxX <= own.maxX + 1e-3 && q.minZ >= own.minZ - 1e-3 && q.maxZ <= own.maxZ + 1e-3));
  let fMax = 1;
  for (let i = 1; i <= 60; i++)
    if (leafHits(d, i / 60, others, bounds)) {
      fMax = (i - 1) / 60;
      break;
    }
  d.max = o.max * fMax;
  const r = d.rect;
  const take = (x: number, z: number) => {
    r.minX = Math.min(r.minX, x);
    r.maxX = Math.max(r.maxX, x);
    r.minZ = Math.min(r.minZ, z);
    r.maxZ = Math.max(r.maxZ, z);
  };
  for (let i = 0; i <= 16; i++) {
    const f = i / 16;
    leafDir(d, f, _dv1);
    leafNormal(d, f, _dv2);
    for (const u of [0, d.len])
      for (const v of [-d.back - 0.01, d.hOut + 0.03]) take(d.hinge.x + _dv1.x * u + _dv2.x * v, d.hinge.z + _dv1.z * u + _dv2.z * v);
  }
  r.minX -= 0.02;
  r.minZ -= 0.02;
  r.maxX += 0.02;
  r.maxZ += 0.02;
  return d;
}

/** When a door is worked on a run, seconds from the leg's start: the hand on the pull, letting go, open, shutting, shut. */
interface DoorPlan {
  d: Door;
  t0: number;
  t1: number;
  t2: number;
  c0: number;
  c1: number;
  /** The hand that pulls: 0 left, 1 right. */
  arm: 0 | 1;
  /** A door with no closer (a combi's), shut by hand: the hand on its pull from c0 brings it to `fRel` by `rel`,
   *  and lets go; it swings on and latches at c1. The hand that does it, and where the closer stands. */
  rel?: number;
  fRel?: number;
  closeArm?: 0 | 1;
  /** Shut from its inner side, the hand on its free edge (not its pull). */
  edge?: boolean;
  /** Whether the run drives the door now (from the hand on its pull to its latch), and whether it has latched: the
   *  door is another run's before and after (two runs to one walk-in or reach-in, one after the other). */
  driving?: boolean;
  latched?: boolean;
}
/** The door's opening at time t of its plan: pulled (speeding up), swinging on (slowing to its stop), held, shutting. */
function doorOpening(p: DoorPlan, t: number): number {
  const d = p.d;
  if (t <= p.t0 || t >= p.c1) return 0;
  if (t < p.t1) {
    const u = (t - p.t0) / (p.t1 - p.t0);
    return d.pullTo * u * u;
  }
  if (t < p.t2) {
    const u = (t - p.t1) / (p.t2 - p.t1);
    return d.pullTo + (1 - d.pullTo) * (1 - (1 - u) * (1 - u));
  }
  if (t < p.c0) return 1;
  if (p.rel !== undefined && p.fRel !== undefined) {
    // By hand: drawn round from the stop, faster as it comes, let go, and on to the latch at the speed it had.
    if (t < p.rel) {
      const u = (t - p.c0) / (p.rel - p.c0);
      return 1 - (1 - p.fRel) * 0.5 * (u * u + u);
    }
    return p.fRel * Math.max(0, 1 - (t - p.rel) / (p.c1 - p.rel));
  }
  // A closer: slow from the stop, quicker through the middle, checked before the latch.
  return 1 - smoothstep(p.c0, p.c1, t);
}
/**
 * A door shut by hand from `fRel` of the way round: seconds the hand is on its pull (from the stop to fRel), and
 * seconds it then swings on to the latch at the speed the hand gave it.
 */
function handShut(fRel: number): { hand: number; swing: number } {
  const hand = 0.35 + 1.2 * (1 - fRel);
  const swing = Math.max(0.5, Math.min(2.2, (fRel * hand) / (1.5 * Math.max(0.05, 1 - fRel))));
  return { hand, swing };
}

/** A full hotel pan of food (NSF-HOTEL-PAN-FULL-4IN), its frame at the rim's centre, its long side along x. */
function makePan(food: number): THREE.Group {
  const g = new THREE.Group();
  g.name = "hotel-pan";
  const m = buildFixture("NSF-HOTEL-PAN-FULL-4IN", { w: HOTEL_PAN.length, d: HOTEL_PAN.width, h: HOTEL_PAN.depth }, 1, true);
  m.hooks.set?.("food", food);
  m.hooks.set?.("fill", 0.7);
  g.userData.setFill = (f: number) => m.hooks.set?.("fill", f);
  m.group.position.set(-HOTEL_PAN.length / 2, -HOTEL_PAN.depth, HOTEL_PAN.width / 2);
  g.add(m.group);
  g.traverse((o) => {
    if ((o as THREE.Mesh).isMesh) o.castShadow = true;
  });
  return g;
}

/** A place a full hotel pan stands: on a shelf of a cold store, or on the line (a table, the makeline's board). */
interface Spot {
  key: string;
  cold: boolean;
  fx: Fx;
  /** The pan's rim centre (world), the surface under it, and the heading of whoever picks or places it. */
  at: THREE.Vector3;
  surfaceY: number;
  yaw: number;
  /** The shelf's or counter's edge in front of the pan (floor point) and its outward normal. */
  face: XZ;
  n: XZ;
  pan: THREE.Object3D | null;
  /** The person whose run will take a pan from or put one on it, or -1. */
  booked: number;
  /** The pot sink's drainboard (where emptied pans go); when the pan on the line is used up; when the one at the
   *  sink is washed and put away. */
  dish?: boolean;
  usedAt?: number;
  washAt?: number;
  /** A reach-in's door in front of it, or the walk-in it is in. */
  door?: Door;
  walkin?: WalkIn;
  /** The pan stands here turned (its long side away from whoever takes it): a combi's rail, the counter its pans
   *  wait on (`stage`); the combi it belongs to. */
  turn?: boolean;
  stage?: boolean;
  combi?: Combi;
}
/** The heading of a pan standing on a spot (its long side across the spot's front, or along its depth if turned). */
function panYaw(sp: Spot): number {
  return sp.yaw - (sp.turn ? Math.PI / 2 : 0);
}

/**
 * A combi oven people load: its door (no closer: opened and shut by hand), the rail a pan goes onto lengthwise, the
 * counter its pans wait on, and when it is next due to be loaded, or its pan done.
 */
interface Combi {
  fx: Fx;
  key: string;
  door: Door;
  spot: Spot;
  stage: Spot;
  /** The rim's height on the rail (world); the last stretch the pan is pushed home by its near end. */
  holdY: number;
  push: number;
  /** The cavity's centre line on the floor (its front face) and the opening's and the pan's depths behind it. */
  face: XZ;
  z1: number;
  zc: number;
  /** Where the door is shut from: clear of its sweep, facing the way its pull (or, from its inner side, its free
   *  edge) comes round. */
  close: XZ;
  closeYaw: number;
  closeArm: 0 | 1;
  closeEdge: boolean;
  /** How far round the hand brings the door before letting it swing shut. */
  fRel: number;
  loadAt: number;
  doneAt: number;
  /** Runs found but each with its door left open too long (the counter is too far off by the ways there): after a
   *  few, it is left alone. */
  far: number;
}

/** A walk-in cooler people go into: its door, the gate over its doorway and floor, its strip curtain. */
interface WalkIn {
  fx: Fx;
  key: string;
  door: Door;
  gate: Rect;
  /**
   * The strips (fixture frame x, hanging in the plane z = rz), how far each is swung now (about x, about z), and the
   * model's own sway under it (what its tick set) with what was drawn on top of it last frame.
   */
  strips: { j: THREE.Object3D; x: number; ax: number; az: number; base: THREE.Quaternion; mine: THREE.Quaternion }[];
  rz: number;
  inv: THREE.Matrix4;
  /** How far its door is open now (0..1). */
  open: number;
}

/** A speed rack people push: between its place (home) and a dock along its own axis, pushed from either end. */
interface Rack {
  fx: Fx;
  key: string;
  home: THREE.Vector3;
  away: THREE.Vector3;
  /** Home to away, a unit floor vector along the rack's depth; the fixture's front normal. */
  d: XZ;
  n: XZ;
  out: boolean;
  busy: boolean;
  /** The casters' swivel now (people.ts: 0 trailing to the back, pi to the front). */
  swivel: number;
  half: XZ;
  base: THREE.Vector3;
  /** Where the kitchen's plan has it (its model's own place); the dock it goes to next after its first push, if its
   *  home is in the way and it shuttles between two others; whether it stays where it is once out. */
  origin: THREE.Vector3;
  next: THREE.Vector3 | null;
  stay: boolean;
}
function rackRect(r: Rack, c: THREE.Vector3): Rect {
  return { minX: c.x - r.half[0], maxX: c.x + r.half[0], minZ: c.z - r.half[1], maxZ: c.z + r.half[1] };
}

type RunKind = "reachin" | "walkin" | "shelf" | "table" | "rack" | "combi";

/** What a leg that is a run does besides walking, timed from the leg's start (seconds). */
interface Run {
  kind: RunKind;
  times: number[];
  door?: DoorPlan;
  pan?: { obj: THREE.Object3D; from: Spot; to: Spot; t0: number; t1: number; held: boolean; placed: boolean };
  walkin?: WalkIn;
  rack?: { r: Rack; from: THREE.Vector3; to: THREE.Vector3; A: XZ; d: XZ; t0: number; t1: number; s0: number; s1: number; parked: boolean };
  /** When the run's bookings are settled (the pan placed, the rack let go). */
  end: number;
  settled: boolean;
}

/** World direction to body space for a person facing yaw. */
function dirToBody(v: THREE.Vector3, yaw: number, out = new THREE.Vector3()): THREE.Vector3 {
  const c = Math.cos(yaw);
  const s = Math.sin(yaw);
  return out.set(c * v.x - s * v.z, v.y, s * v.x + c * v.z);
}

/** A hand on a bar or handle (world point, the bar's axis towards the thumb, the palm's facing), weight w. */
function gripArm(st: Staff, p: BodyPose, i: 0 | 1, w: number, at: THREE.Vector3, axis: THREE.Vector3, palm: THREE.Vector3, close: number, radius = 0.011): void {
  const A = p.arms[i];
  const yaw = st.placement.yaw;
  A.ik = w;
  A.grip = "bar";
  A.roll = 0.5;
  A.spinMax = 0;
  st.toBody(at, A.pos);
  dirToBody(axis, yaw, A.a).normalize();
  dirToBody(palm, yaw, A.b).normalize();
  A.pole.set(i === 0 ? 0.55 : -0.55, -0.8, -0.25).normalize();
  A.shape = "bar";
  A.close = close;
  A.radius = radius;
}

const UP = new THREE.Vector3(0, 1, 0);
const _rv1 = new THREE.Vector3();
const _rv2 = new THREE.Vector3();
const _rv3 = new THREE.Vector3();
const _rq = new THREE.Quaternion();
const _re = new THREE.Euler();
const _fl1 = new THREE.Vector3();
const _fl2 = new THREE.Vector3();
const _fl3 = new THREE.Vector3();
const _fl4 = new THREE.Vector3();
const _fl5 = new THREE.Vector3();
const _fl6 = new THREE.Vector3();
/** The outward normal of a shut door, on the floor. */
function doorNormal(d: Door): XZ {
  return [d.e1.x, d.e1.z];
}
function rectsTouch(a: Rect, b: Rect, m: number): boolean {
  return a.minX - m < b.maxX && b.minX - m < a.maxX && a.minZ - m < b.maxZ && b.minZ - m < a.maxZ;
}

/** Where a routine puts the body's root at time t, and whether it is walking, without posing the body. */
function rootSampler(r: Routine): (t: number) => { x: number; z: number; yaw: number; walking: boolean } {
  const out = { x: 0, z: 0, yaw: 0, walking: false };
  return (t: number) => {
    const p = r.rootAt(t);
    out.x = p.x;
    out.z = p.z;
    out.yaw = p.yaw;
    out.walking = r.activityAt(t) === "walk";
    return out;
  };
}

/** A candidate leg with what it is: its steps and path, and the run it makes (null for a station). */
type Cand = LegCandidate & { steps: RoutineStep[]; path: XZ[]; run: Run | null; closing: number };

/** How a run's footprint goes: the steps between which a pan is held (after the pick, to the place), a rack pushed. */
interface RunShape {
  pan?: [number, number];
  rack?: [number, number];
  /** The pan held turned; the steps where it is in or over a fixture (a combi's opening), its footprint left off
   *  there unless it is clear. */
  turn?: boolean;
  inside?: readonly number[];
}

/** A run's steps with its marks: the pick, the place, the door pulled, the rack taken and let go. */
interface RunSeq {
  kind: RunKind;
  steps: RoutineStep[];
  /** Walks still to be given their paths: step index, from, to, carrying, through a gate. */
  walks: { i: number; from: XZ; to: XZ; carry: boolean; arrive?: XZ }[];
  pick: number;
  place: number;
  /** The stand step at which the door is pulled open, and the step after which (once clear of its swing) it shuts. */
  open: number;
  closeAfter: number;
  /** The first step by which the door is shut again (the walks from it on are planned without its leaf). */
  shut?: number;
  door?: Door;
  arm: 0 | 1;
  walkin?: WalkIn;
  from?: Spot;
  to?: Spot;
  rack?: { r: Rack; from: THREE.Vector3; to: THREE.Vector3; A: XZ; d: XZ; grab: number; walk: number; release: number; front: boolean };
  /** A combi loaded or unloaded: the pan held turned; the steps where it is in the oven's opening; the stand step
   *  at which the hand shuts the door. */
  combi?: Combi;
  inside?: number[];
  closeAt?: number;
  /** A candidate was turned down for leaving the door open too long. */
  long?: boolean;
}

class Crew {
  readonly group = new THREE.Group();
  readonly agents: Agent[] = [];
  readonly crowd: Crowd;
  readonly stations: Station[];
  readonly lines = new THREE.Group();
  /** Planning results, for the console and for checks. */
  readonly log = { legs: 0, waits: 0, delays: 0, planMs: 0, planMax: 0, searchMax: 0, runs: {} as Record<string, number>, runFails: {} as Record<string, number>, doorHits: 0, whyNot: {} as Record<string, number> };
  now = 0;
  /** A run to make next (scripted checks): its kind, or null to choose. */
  force: RunKind | null = null;
  readonly spots: Spot[] = [];
  readonly walkins: WalkIn[] = [];
  readonly racks: Rack[] = [];
  readonly combis: Combi[] = [];
  private nextPan = 5;
  /** Pan runs (from > to) that could not be made, and until when they are not tried again. */
  private readonly cool = new Map<string, number>();
  private nextRack = 16;
  private ghosts = new Map<Cohort, Staff>();
  private rnd: () => number;
  private scratch = newPose();
  private lastPrune = 0;
  private standCache = new Map<string, XZ | null>();

  /** The room's stations, runs and tools; its people come with `populate`. */
  constructor(readonly K: Kitchen) {
    this.group.name = "crew";
    this.lines.name = "paths";
    this.rnd = rngSeed(K.fixtures.length * 7919 + Math.round(K.W * 100));
    this.stations = makeStations(K, this.group);
    // The cold stores, the line's pan places and the racks that roll: the racks people push leave the static floor
    // (they are parked rectangles in the crowd instead), the walk-ins' floors are gated.
    const [grid, tall] = this.setupRuns();
    this.crowd = new Crowd(grid, CROWD_DT, tall);
    for (const w of this.walkins) this.crowd.addGate(w.key, w.gate);
    for (const r of this.racks) this.crowd.addTimedRect(rackRect(r, r.home), 0, Infinity, -1, r.key);
    for (const s of this.stations) {
      const tk = TOOL_OF[s.kind];
      if (tk && !s.tool) {
        s.tool = makeTool(tk);
        s.rest = restFor(s);
        if (tk !== "skillet" && tk !== "knife") this.own(s.tool);
      }
      if (s.tool && s.rest) {
        restMatrix(s.rest, s.tool.matrix).decompose(s.tool.position, s.tool.quaternion, s.tool.scale);
        this.group.add(s.tool);
      }
      if (s.kind === "plate" && !s.noPlate) {
        const pl = makePlate();
        pl.position.copy(s.at);
        this.group.add(pl);
        this.own(pl);
        if (s.rest) {
          const bain = makeBain();
          bain.position.copy(s.rest.p).setY(s.at.y);
          this.group.add(bain);
          this.own(bain);
        }
      }
      if (s.kind === "flip") this.setupGrill(s);
      if (s.baskets) s.baskets.forEach((b, k) => s.fx.model.hooks.set?.(k === 0 ? "lift_l" : "lift_r", b.down ? 0 : 1));
    }
  }

  /**
   * Who works here: the spec's head count (3 to 10), its routes' cohorts and areas as hints. Each person is built and
   * placed in a task of their own; false once `alive` turns false between them (the caller disposes the crew).
   */
  async populate(detailed: boolean, phone: boolean, alive: () => boolean): Promise<boolean> {
    const spec = this.K.spec;
    const n = Math.max(Math.min(3, this.stations.length), Math.min(10, spec.staffCount, Math.floor(this.stations.length * 0.7)));
    const routes = spec.staffRoutes;
    const cohorts: Cohort[] = ["p50", "p05", "p95"];
    const taken = new Set<Station>();
    for (let i = 0; i < n; i++) {
      await new Promise((r) => setTimeout(r, 0));
      if (!alive()) return false;
      const route = i < routes.length ? routes[i] : undefined;
      const cohort: Cohort = route ? route.cohort : cohorts[i % 3]!;
      const st = makeStaff(cohort, phone ? 1 : 2, { seed: 11 + i * 17, body: detailed ? "detailed" : "procedural", shadows: !phone });
      let ghost = this.ghosts.get(cohort);
      if (!ghost) {
        ghost = makeStaff(cohort, 1, { body: "procedural", shadows: false });
        this.ghosts.set(cohort, ghost);
      }
      const { shapes, rWalk, rackFace, gripY } = footprints(st);
      const a: Agent = {
        id: -1,
        st,
        ghost,
        cohort,
        rWalk,
        rPath: rWalk + STATIC_MARGIN + 0.012,
        dFace: rWalk + STATIC_MARGIN + 0.016,
        shapes,
        rackFace,
        gripY,
        lean: 0,
        speed: COHORT_SPEED[cohort],
        routine: null as unknown as Routine,
        steps: [],
        start: { x: 0, z: 0, yaw: 0 },
        t0: 0,
        endT: 0,
        station: null,
        run: null,
        workDur: 0,
        planned: null,
        fails: 0,
        weights: new Map(),
        path: null,
        line: null,
        toolK: 0,
        holding: false,
        accent: new THREE.Color(route?.accentHex ?? "#f59e0b"),
      };
      // Prefer the stations near the route's floor area.
      const pts = (route?.waypoints ?? []).map(([x, y]) => [x, -y] as const);
      for (const s of this.stations) {
        const d = pts.length ? Math.min(...pts.map(([x, z]) => Math.hypot(s.faceX - x, s.faceZ - z))) : 3;
        a.weights.set(s, Math.exp(-d / 2.2) + 0.06);
      }
      // Start at the best free station whose standing spot is clear of the room and of everyone placed so far.
      const order = [...this.stations].sort((p, q) => a.weights.get(q)! - a.weights.get(p)!);
      let placed = false;
      for (const s of order) {
        if (taken.has(s)) continue;
        const sp = this.stance(a, s);
        if (!sp) continue;
        if (this.agents.some((b) => Math.hypot(b.start.x - sp.x, b.start.z - sp.z) < a.rWalk + b.rWalk + PAIR_MARGIN + 0.05)) continue;
        const id = this.crowd.addAgent({ shapes }, sp.x, sp.z, sp.yaw, 0, 0);
        const free = this.crowd.standingObstacles(id, 0, 0)(sp.x, sp.z);
        if (this.crowd.staticRoom(id, sp.x, sp.z, sp.yaw, 0) < STATIC_MARGIN || free < a.rWalk) {
          this.crowd.agents.pop();
          this.crowd.tracks.pop();
          continue;
        }
        a.id = id;
        taken.add(s);
        s.holder = id;
        const dur = 2 + this.rnd() * 6;
        const steps = [this.workStep(s, dur)];
        this.begin(a, sp, steps, CROWD_DT, s, null, dur, null);
        this.crowd.commit(id, 1, 0, this.realize(a, sp, steps, dur));
        a.endT = this.crowd.end(id) * CROWD_DT;
        placed = true;
        break;
      }
      if (!placed) {
        st.dispose();
        continue;
      }
      this.agents.push(a);
      this.group.add(st.group);
    }
    return true;
  }

  /* -------------------------------------------------------------------------------------------------------- */
  /* Setting up: cold stores, the line's pan places, racks                                                     */

  /** Find the walk-ins, the reach-ins, the shelving and tables pans go on, and the racks that can roll. */
  private setupRuns(): [NavGrid, NavGrid] {
    const K = this.K;
    const toWorldDir = (fx: Fx, x: number, z: number) => new THREE.Vector3(x, 0, z).transformDirection(fx.model.group.matrixWorld).setY(0).normalize();
    const xz = (v: THREE.Vector3): XZ => [v.x, v.z];
    /** Whether a pan's box (fixture frame: centre x, z, half extents, from the shelf up) is clear of every triangle there. */
    const clearOn = (fx: Fx, cx: number, cz: number, hx: number, hz: number, y: number): boolean => {
      const b = fx.tris;
      if (!b) return true;
      const [x0, x1, z0, z1, y0, y1] = [cx - hx, cx + hx, cz - hz, cz + hz, y + 0.012, y + 0.14];
      for (let i = 0; i < b.length; i += 6) if (b[i]! < x1 && b[i + 3]! > x0 && b[i + 1]! < y1 && b[i + 4]! > y0 && b[i + 2]! < z1 && b[i + 5]! > z0) return false;
      return true;
    };
    const addSpot = (sp: Omit<Spot, "pan" | "booked">) => this.spots.push({ ...sp, pan: null, booked: -1 });

    // Walk-in coolers: the door, the gate, the strips, pans on the shelving down both sides.
    for (const fx of K.fixtures) {
      const k = WALKINS[fx.id];
      if (!k || !fx.rect) continue;
      const [w, d] = fx.f.local_extents_m;
      const T = PANEL;
      const zi = WALL_FRONT - T;
      const hx = k.door0 - 0.049;
      const hz = -0.039;
      // The hand closes on the latch's pull grip (its handle_pull socket): 41 mm in from the slab's free edge, its
      // bar's axis 34 mm proud of the hinge line. Its height is read off the model, the grip being the leaf's outermost
      // part: a freezer's door stands on the insulated floor, a cooler's on the screed.
      const gripY = (() => {
        const leaf = fx.model.joints.get("cam_lift_door");
        if (!leaf) return null;
        leaf.updateWorldMatrix(true, true);
        const inv = new THREE.Matrix4().copy(leaf.matrixWorld).invert();
        const M = new THREE.Matrix4();
        const v = new THREE.Vector3();
        const pts: number[] = [];
        leaf.traverse((o) => {
          const m = o as THREE.Mesh;
          const pos = m.isMesh ? (m.geometry as THREE.BufferGeometry).getAttribute("position") : undefined;
          if (!pos) return;
          M.copy(inv).multiply(m.matrixWorld);
          for (let i = 0; i < pos.count; i++) {
            v.fromBufferAttribute(pos, i).applyMatrix4(M);
            pts.push(v.y, v.z);
          }
        });
        let zMax = -Infinity;
        for (let i = 1; i < pts.length; i += 2) zMax = Math.max(zMax, pts[i]!);
        let y0 = Infinity;
        let y1 = -Infinity;
        for (let i = 0; i < pts.length; i += 2)
          if (pts[i + 1]! > zMax - 0.003) {
            y0 = Math.min(y0, pts[i]!);
            y1 = Math.max(y1, pts[i]!);
          }
        return y1 >= y0 ? v.set(0, (y0 + y1) / 2, 0).applyMatrix4(leaf.matrixWorld).y : null;
      })();
      const door = makeDoor(fx, {
        name: "door",
        hinge: fxThree(fx, hx, 0, hz),
        e0: toWorldDir(fx, 1, 0),
        e1: toWorldDir(fx, 0, 1),
        len: k.doorW + 0.094,
        max: 1.92,
        hAlong: k.doorW - 0.041,
        hOut: 0.034,
        hY: gripY ?? fxThree(fx, 0, 1.08, 0).y,
        back: 0.07,
        pullTo: 0.3,
        pull: 1.05,
        close: 2.4,
      }, K.grid.rects, K.grid.bounds);
      const inv = new THREE.Matrix4().copy(fx.model.group.matrixWorld).invert();
      const strips: WalkIn["strips"] = [];
      for (const [name, j] of fx.model.joints)
        if (name.startsWith("pvc_strip_")) strips.push({ j, x: j.getWorldPosition(new THREE.Vector3()).applyMatrix4(inv).x, ax: 0, az: 0, base: new THREE.Quaternion(), mine: new THREE.Quaternion() });
      const wi: WalkIn = {
        fx,
        key: `walkin:${fx.f.fixture_id}`,
        door,
        gate: fxRect(fx, T, w - T, -d + T, 0),
        strips,
        rz: zi - 0.022,
        inv,
        open: 0,
      };
      this.walkins.push(wi);
      const { Ds, lens, z0 } = walkInUnits(w, d);
      for (const side of [1, -1]) {
        const xf = side > 0 ? T + 0.012 + Ds : w - T - 0.012 - Ds;
        const xc = xf - side * (0.02 + HOTEL_PAN.width / 2);
        const n = toWorldDir(fx, side, 0);
        let zu = z0;
        for (const L of lens) {
          for (const tier of [0.2, 0.7167, 1.2333, 1.75]) {
            // The first two free places along the unit, a pan's length apart.
            let taken = -Infinity;
            let count = 0;
            for (let zc = zu - 0.32; zc > zu - L + 0.32 && count < 2; zc -= 0.04) {
              if (Math.abs(zc - taken) < 0.6) continue;
              if (!clearOn(fx, xc, zc, HOTEL_PAN.width / 2 + 0.01, HOTEL_PAN.length / 2 + 0.015, tier)) continue;
              taken = zc;
              count++;
              const at = fxThree(fx, xc, tier + HOTEL_PAN.depth, zc);
              addSpot({ key: `${fx.f.fixture_id}:${side}:${tier}:${zc.toFixed(2)}`, cold: true, fx, at, surfaceY: at.y - HOTEL_PAN.depth, yaw: Math.atan2(-n.x, -n.z), face: xz(fxThree(fx, xf, 0, zc)), n: [n.x, n.z], walkin: wi });
            }
          }
          zu -= L + 0.01;
        }
      }
    }
    // Reach-ins: a door to each section, pans on its low and middle shelves.
    for (const fx of K.fixtures) {
      if (!/REACHIN-2D/.test(fx.id) || !fx.rect) continue;
      const [w, , h] = fx.f.local_extents_m;
      const n = toWorldDir(fx, 0, 1);
      const dw = w / 2 - 0.003 - 0.023;
      const yb = 0.152 + 0.071;
      const top = h - 0.195 - 0.056 - 0.22;
      const hY = (0.192 + (h - 0.195 - 0.025)) / 2;
      const mid = w / 2;
      for (const side of [1, -1] as const) {
        const px = side > 0 ? 0.035 : w - 0.035;
        const door = makeDoor(fx, {
          name: side > 0 ? "door_left" : "door_right",
          hinge: fxThree(fx, px, 0, -0.065),
          e0: toWorldDir(fx, side, 0),
          e1: n.clone(),
          len: dw - 0.012,
          max: 2.094,
          hAlong: dw - 0.012 - 0.032,
          hOut: 0.054,
          hY,
          back: 0.035,
          pullTo: 0.38,
          pull: 0.85,
          close: 1.4,
        }, K.grid.rects, K.grid.bounds);
        // A door that cannot open past a right angle here leaves no room to take a pan out.
        if (door.max < 1.5) continue;
        const x0 = side > 0 ? 0.054 : mid + 0.012;
        const x1 = side > 0 ? mid - 0.012 : w - 0.054;
        const cx = (x0 + x1) / 2;
        for (const k of [0.27, 0.53]) {
          const tier = yb + (top - yb) * k + 0.006;
          const at = fxThree(fx, cx, tier + HOTEL_PAN.depth, -0.2 - HOTEL_PAN.width / 2);
          addSpot({ key: `${fx.f.fixture_id}:${side}:${k}`, cold: true, fx, at, surfaceY: at.y - HOTEL_PAN.depth, yaw: Math.atan2(-n.x, -n.z), face: xz(fxThree(fx, cx, 0, 0)), n: [n.x, n.z], door });
        }
      }
    }
    // Wire shelving in the cold store's room (a commissary's): pans on its two middle shelves.
    for (const fx of K.fixtures) {
      if (!/WIRE-SHELVING/.test(fx.id) || !/cold/i.test(fx.f.fixture_id) || !fx.rect) continue;
      const [w, , h] = fx.f.local_extents_m;
      const top = h - 0.13;
      const n = toWorldDir(fx, 0, 1);
      const tiers = [0.2, 0.2 + (top - 0.2) / 3, 0.2 + ((top - 0.2) * 2) / 3, top];
      for (const tier of tiers) {
        let taken = -Infinity;
        for (let cx = 0.3; cx < w - 0.3; cx += 0.04) {
          if (Math.abs(cx - taken) < 0.58) continue;
          if (!clearOn(fx, cx, -0.02 - HOTEL_PAN.width / 2, HOTEL_PAN.length / 2 + 0.015, HOTEL_PAN.width / 2 + 0.01, tier)) continue;
          taken = cx;
          const at = fxThree(fx, cx, tier + HOTEL_PAN.depth, -0.02 - HOTEL_PAN.width / 2);
          addSpot({ key: `${fx.f.fixture_id}:${tier.toFixed(2)}:${cx.toFixed(2)}`, cold: true, fx, at, surfaceY: at.y - HOTEL_PAN.depth, yaw: Math.atan2(-n.x, -n.z), face: xz(fxThree(fx, cx, 0, 0)), n: [n.x, n.z] });
        }
      }
    }
    // What stands on the counters (smallwares, heat lamps over a pass).
    const smalls = K.fixtures.filter((f) => SMALLWARE.test(f.id) && f.f.world_aabb_m);
    const lamps = K.fixtures.filter((f) => /HEATLAMP/.test(f.id) && f.f.world_aabb_m);
    const overlaps = (r: Rect, f: Fx) => {
      const b = f.f.world_aabb_m!;
      return b.min_m[0] < r.maxX && b.max_m[0] > r.minX && -b.max_m[1] < r.maxZ && -b.min_m[1] > r.minZ;
    };
    // Combi ovens: a pan goes onto the rail nearest a cook's waist, lengthwise, through a door that has no closer;
    // the pans wait for it, turned, on the nearest free stretch of counter. (Why one is left out, for checks.)
    const skip = (why: string) => {
      this.log.whyNot[`combi: ${why}`] = (this.log.whyNot[`combi: ${why}`] ?? 0) + 1;
    };
    for (const fx of K.fixtures) {
      if (!fx.rect) continue;
      const [w, d, h] = fx.f.local_extents_m;
      const cav = combiCavity(fx.id, { w, d, h });
      if (!cav) continue;
      const doorJ = fx.model.joints.get("cavity_door");
      if (!doorJ) {
        skip("no door joint");
        continue;
      }
      const lip = HOTEL_PAN.lip;
      const railY = (y: number) => fxThree(fx, 0, y + lip, 0).y;
      let rail = cav.rails[0]!;
      for (const y of cav.rails) if (Math.abs(railY(y) - 1.0) < Math.abs(railY(rail) - 1.0)) rail = y;
      const n = toWorldDir(fx, 0, 1);
      const yaw = Math.atan2(-n.x, -n.z);
      const at = fxThree(fx, cav.cx, rail + lip, cav.zc);
      // The door (ovens.ts combi): hinged at the left on pins in the plane of its outer face, so that it swings
      // forward; the grip on standoffs near its free edge is all that stands out from that face (the knuckles turn in
      // place about the pins).
      const big = cav.rails.length >= 20;
      const D = cav.door;
      const leaf = new THREE.Box3().setFromObject(doorJ);
      const ym = (leaf.min.y + leaf.max.y) / 2;
      const door = makeDoor(
        fx,
        {
          name: "open",
          hinge: doorJ.getWorldPosition(new THREE.Vector3()),
          e0: toWorldDir(fx, 1, 0),
          e1: n.clone(),
          len: D.len,
          max: 2.356,
          hAlong: D.hAlong,
          hOut: D.hOut,
          hY: big ? Math.min(ym, fxThree(fx, 0, 1.0, 0).y) : ym,
          back: D.back,
          pullFrom: D.pullFrom,
          face: 0.002,
          pullTo: 0.38,
          pull: 0.9,
          close: 1.6,
        },
        K.grid.rects,
        K.grid.bounds,
      );
      // Opened at least square to the front (a pan goes in past its leaf).
      if (door.max < 1.6) {
        skip("its door is blocked");
        continue;
      }
      const face = xz(fxThree(fx, cav.cx, 0, 0));
      // The counter: free of what stands on it, of the stations and of the line's pan places; nearest first.
      let stage: Spot | null = null;
      // (A pass under heat lamps too, its pan well under them, but as if it were further off.) Either side of a
      // table, or of an undercounter's worktop, that has room to stand at it, the side nearer the oven first.
      const sides: { t: Fx; sg: 1 | -1; dist: number }[] = [];
      for (const t of K.fixtures) {
        if (!/PREPTABLE|UNDERCOUNTER/.test(t.id) || !t.rect || t.f.local_extents_m[1] < 0.66) continue;
        const [tw, td] = t.f.local_extents_m;
        const lamp = lamps.some((l) => overlaps(t.rect!, l)) ? 4 : 0;
        for (const sg of [1, -1] as const) {
          const z = sg > 0 ? 0.45 : -td - 0.45;
          if (K.grid.clearance(fxThree(t, tw / 2, 0, z).x, fxThree(t, tw / 2, 0, z).z) < 0.25) continue;
          const c0 = fxThree(t, tw / 2, 0, sg > 0 ? 0 : -td);
          const dist = Math.hypot(c0.x - face[0], c0.z - face[1]) + lamp;
          if (dist < 9.5) sides.push({ t, sg, dist });
        }
      }
      sides.sort((p, q) => p.dist - q.dist);
      for (const { t, sg } of sides) {
        const [tw, td, th] = t.f.local_extents_m;
        const topY = th > 0.95 ? th : 0.914;
        const tn = toWorldDir(t, 0, sg);
        // (Its edge, and the pan's centre, in the table's own frame.)
        const z0 = sg > 0 ? 0 : -td;
        const zc = z0 - sg * (0.06 + HOTEL_PAN.length / 2);
        const mine = this.stations.filter((q) => q.fx === t);
        // The plates at a pass on it (a heat lamp's station over it): the pan and its way in keep off them.
        const passes = this.stations.filter((q) => q.fx !== t && /HEATLAMP/.test(q.fx.id) && overlaps(t.rect!, q.fx));
        const far = (lx: number) => {
          const v = fxThree(t, lx, 0, zc);
          return Math.hypot(v.x - face[0], v.z - face[1]);
        };
        // (Every 5 cm along it, so that the pan can wait as near the oven as the plates and the stations leave room.)
        const xs: number[] = [];
        for (let lx = 0.24; lx <= tw - 0.24 + 1e-6; lx += 0.05) xs.push(lx);
        xs.sort((p, q) => far(p) - far(q));
        for (const lx of xs) {
          const pa = fxThree(t, lx, topY + HOTEL_PAN.depth, zc);
          // Its place, and the counter between it and the edge (it is drawn across that, held just over it).
          const r = fxRect(t, lx - 0.19, lx + 0.19, Math.min(z0, zc - sg * 0.29), Math.max(z0, zc - sg * 0.29));
          if (smalls.some((f) => overlaps(r, f))) continue;
          if (mine.some((q) => Math.hypot(q.at.x - pa.x, q.at.z - pa.z) < 0.55)) continue;
          if (passes.some((q) => rectDistance(r, q.at.x, q.at.z) < 0.16)) continue;
          if (this.spots.some((q) => Math.hypot(q.at.x - pa.x, q.at.z - pa.z) < (q.turn ? 0.45 : 0.55))) continue;
          let board = false;
          t.model.group.traverse((o) => {
            if (o !== t.model.group && o.parent === t.model.group && o.children.length && Math.abs(o.position.x + 0.305 - lx) < 0.6 && o.position.y > 0.9) board = true;
          });
          if (board) continue;
          stage = { key: `${t.f.fixture_id}:stage:${sg > 0 ? "" : "back:"}${lx.toFixed(2)}`, cold: false, fx: t, at: pa, surfaceY: topY, yaw: Math.atan2(-tn.x, -tn.z), face: xz(fxThree(t, lx, 0, z0)), n: [tn.x, tn.z], pan: null, booked: -1, turn: true, stage: true };
          break;
        }
        if (stage) break;
      }
      if (!stage) {
        skip("no counter near");
        continue;
      }
      // Where the door is shut from: within reach of its pull at the stop, just clear of the leaf's whole sweep (it
      // passes in front of them); the hand draws the pull round as far as it reaches and lets go, and the door swings
      // on and latches. The one that brings it furthest round, then the nearest. (The planner checks every body
      // against the leaf; here, a chest's depth.) Where the open leaf leaves no way round its free end (it stands
      // across a narrow aisle), it is shut from the side the oven is worked from, its inner face: the hand takes it by
      // its free edge and draws it round towards them. They face it as it comes, so half a chest's depth and a hand's
      // breadth keep them clear; the room for that is narrow, so it is looked for more finely.
      const tip = leafDir(door, 1, new THREE.Vector3()).multiplyScalar(door.len + 0.33).add(door.hinge);
      const closeEdge = K.grid.clearance(tip.x, tip.z) < 0.27;
      const nIn = leafNormal(door, 1, new THREE.Vector3());
      const grip = (f: number, out: THREE.Vector3) => (closeEdge ? edgeAt(door, f, out) : handleAt(door, f, out));
      const h1 = grip(1, new THREE.Vector3());
      const clearOf = closeEdge ? 0.24 : 0.3;
      const radii = closeEdge ? [0.34, 0.37, 0.4, 0.43, 0.46, 0.49, 0.52, 0.55] : [0.34, 0.4, 0.46, 0.52];
      const turns = closeEdge ? 96 : 32;
      let best: { p: XZ; fRel: number; r: number } | null = null;
      for (const r of radii)
        for (let k = 0; k < turns; k++) {
          const x = h1.x + r * Math.cos((k / turns) * 2 * Math.PI);
          const z = h1.z + r * Math.sin((k / turns) * 2 * Math.PI);
          if (K.grid.clearance(x, z) < 0.2) continue;
          if (closeEdge && (x - door.hinge.x) * nIn.x + (z - door.hinge.z) * nIn.z > 0) continue;
          let clear = true;
          for (let i = 0; i <= 30 && clear; i++) if (leafGap(door, i / 30, x, z) < clearOf) clear = false;
          if (!clear) continue;
          let f = 1;
          while (f > 0.2) {
            const hp = grip(f - 0.025, _rv1);
            if (Math.hypot(hp.x - x, hp.z - z) > 0.56) break;
            f -= 0.025;
          }
          if (f > 0.9) continue;
          if (!best || f < best.fRel - 0.02 || (f < best.fRel + 0.02 && r < best.r)) best = { p: [x, z], fRel: f, r };
        }
      if (!best) {
        skip("nowhere to shut its door from");
        continue;
      }
      const hm = grip((1 + best.fRel) / 2, _rv1);
      const closeYaw = Math.atan2(hm.x - best.p[0], hm.z - best.p[1]);
      const closeArm: 0 | 1 = (h1.x - best.p[0]) * Math.cos(closeYaw) - (h1.z - best.p[1]) * Math.sin(closeYaw) > 0 ? 0 : 1;
      const spot: Spot = { key: `${fx.f.fixture_id}:rail`, cold: false, fx, at, surfaceY: fxThree(fx, 0, rail, 0).y, yaw, face, n: [n.x, n.z], pan: null, booked: -1, turn: true };
      const c: Combi = {
        fx,
        key: `combi:${fx.f.fixture_id}`,
        door,
        spot,
        stage,
        holdY: at.y,
        push: Math.max(0.06, cav.box.z1 - cav.zc - 0.135),
        face,
        z1: cav.box.z1,
        zc: cav.zc,
        close: best.p,
        closeYaw,
        closeArm,
        closeEdge,
        fRel: best.fRel,
        loadAt: 6 + 10 * this.rnd(),
        doneAt: Infinity,
        far: 0,
      };
      spot.combi = stage.combi = c;
      this.spots.push(spot, stage);
      this.combis.push(c);
    }

    // The line: free places on the prep tables (clear of what stands on them, of the stations there and of the
    // combis' counters), else the makeline's board.
    for (const fx of K.fixtures) {
      if (!/PREPTABLE/.test(fx.id) || !fx.rect) continue;
      if (lamps.some((l) => overlaps(fx.rect!, l))) continue;
      const [w, , h] = fx.f.local_extents_m;
      const topY = h > 0.95 ? h : 0.914;
      const n = toWorldDir(fx, 0, 1);
      const mine = this.stations.filter((s) => s.fx === fx);
      let count = 0;
      for (const ox of [0.18, 0.82, 0.5, 0.32, 0.68]) {
        if (count >= 2) break;
        const lx = Math.max(0.3, Math.min(w - 0.3, w * ox));
        const zc = -0.06 - HOTEL_PAN.width / 2;
        const at = fxThree(fx, lx, topY + HOTEL_PAN.depth, zc);
        const r = fxRect(fx, lx - 0.29, lx + 0.29, zc - 0.18, zc + 0.18);
        if (smalls.some((f) => overlaps(r, f))) continue;
        if (mine.some((s) => Math.hypot(s.at.x - at.x, s.at.z - at.z) < 0.62)) continue;
        if (this.spots.some((q) => !q.cold && q.at.distanceTo(at) < (q.turn ? 0.55 : 0.6))) continue;
        // A board of ours may stand there.
        let board = false;
        fx.model.group.traverse((o) => {
          if (o !== fx.model.group && o.parent === fx.model.group && o.children.length && Math.abs(o.position.x + 0.305 - lx) < 0.6 && o.position.y > 0.9) board = true;
        });
        if (board) continue;
        count++;
        addSpot({ key: `${fx.f.fixture_id}:${lx.toFixed(2)}`, cold: false, fx, at, surfaceY: topY, yaw: Math.atan2(-n.x, -n.z), face: xz(fxThree(fx, lx, 0, 0)), n: [n.x, n.z] });
      }
    }
    if (!this.spots.some((q) => !q.cold && !q.turn))
      for (const fx of K.fixtures) {
        if (!/MAKELINE/.test(fx.id) || !fx.rect) continue;
        const [w] = fx.f.local_extents_m;
        const n = toWorldDir(fx, 0, 1);
        for (const lx of [0.36, w - 0.36]) {
          const at = fxThree(fx, lx, 0.9305 + HOTEL_PAN.depth, -0.166);
          addSpot({ key: `${fx.f.fixture_id}:${lx.toFixed(2)}`, cold: false, fx, at, surfaceY: 0.9305, yaw: Math.atan2(-n.x, -n.z), face: xz(fxThree(fx, lx, 0, 0)), n: [n.x, n.z] });
        }
      }
    // Over the first bowl of a pot sink, beside whoever washes: where a pan used up on the line is taken.
    for (const fx of K.fixtures) {
      if (!/3COMP/.test(fx.id) || !fx.rect) continue;
      const [w] = fx.f.local_extents_m;
      const n = toWorldDir(fx, 0, 1);
      for (const lx of [0.36, w - 0.36]) {
        const zc = -0.06 - HOTEL_PAN.width / 2;
        const at = fxThree(fx, lx, 0.914 + HOTEL_PAN.depth, zc);
        if (this.stations.some((q) => q.fx === fx && Math.hypot(q.at.x - at.x, q.at.z - at.z) < 0.5)) continue;
        addSpot({ key: `${fx.f.fixture_id}:dish:${lx.toFixed(2)}`, cold: false, dish: true, fx, at, surfaceY: 0.914, yaw: Math.atan2(-n.x, -n.z), face: xz(fxThree(fx, lx, 0, 0)), n: [n.x, n.z] });
        break;
      }
    }
    // Pans to begin with: most cold places full, one place on the line.
    let food = 0;
    const put = (sp: Spot) => {
      const pan = makePan(food++ % 4);
      pan.position.copy(sp.at);
      pan.rotation.set(0, panYaw(sp), 0);
      this.group.add(pan);
      sp.pan = pan;
    };
    // (On the shelves a cook can reach without a step: most can take a pan from knee to head height.)
    // (A walk-in's are always stocked.)
    const handy = this.spots.filter((sp) => sp.cold && sp.surfaceY > 0.45 && sp.surfaceY < 1.45);
    for (const sp of handy) if (sp.walkin || this.rnd() < 0.62) put(sp);
    // (At least two to fetch wherever there are places for them: a store with few places is not left empty by chance.)
    for (const sp of handy) if (!sp.pan && handy.filter((q) => q.pan).length < 2) put(sp);
    const line = this.spots.filter((q) => !q.cold && !q.dish && !q.turn);
    if (line.length > 1) {
      const sp = line[Math.floor(this.rnd() * line.length)]!;
      put(sp);
      sp.usedAt = 20 + 25 * this.rnd();
    }
    // A pan waiting by each combi.
    for (const c of this.combis) put(c.stage);

    // Racks that can roll: along their own axis, from either end, out to a dock and back, clear of everything.
    const rects = K.grid.rects;
    const roomOK = (r: Rect, m: number, skip: Rect) => {
      const b = K.grid.bounds;
      if (r.minX - m < b.minX || r.maxX + m > b.maxX || r.minZ - m < b.minZ || r.maxZ + m > b.maxZ) return false;
      for (const o of rects) if (o !== skip && o.minX - m < r.maxX && r.minX - m < o.maxX && o.minZ - m < r.maxZ && r.minZ - m < o.maxZ) return false;
      for (const w of this.walkins) if (rectsTouch(w.gate, r, m) || rectsTouch(w.door.rect, r, m)) return false;
      for (const sp of this.spots) if (sp.door && rectsTouch(sp.door.rect, r, m)) return false;
      for (const c of this.combis) if (rectsTouch(c.door.rect, r, m)) return false;
      return true;
    };
    const clearAt = (x: number, z: number, skip: Rect) => {
      const b = K.grid.bounds;
      let q = Math.min(x - b.minX, b.maxX - x, z - b.minZ, b.maxZ - z);
      for (const o of rects) if (o !== skip) q = Math.min(q, rectDistance(o, x, z));
      return q;
    };
    const stationNear = (x: number, z: number) => this.stations.some((s) => Math.hypot(s.faceX + s.nx * 0.3 - x, s.faceZ + s.nz * 0.3 - z) < 0.9);
    // Where people need to get to (the stations, the pan places), and which of them a body can walk between, and a
    // body with a pan held out ahead carry one between: a rack is never left where it would cut one off from another.
    // Judged without the racks that roll (each placed in turn).
    const rolling = new Set(K.fixtures.filter((fx) => MOVABLE.test(fx.id) && fx.rect).map((fx) => fx.rect!));
    const W0 = new NavGrid(K.grid.bounds, rects.filter((r) => !rolling.has(r)), K.grid.cell);
    const T0 = new NavGrid(K.tall.bounds, K.tall.rects.filter((r) => !rolling.has(r)), K.tall.cell);
    const spotPts = this.spots.map((q) => [q.face[0] + q.n[0] * (q.cold ? 0.62 : 0.36), q.face[1] + q.n[1] * (q.cold ? 0.62 : 0.36)] as XZ);
    const walkPts: XZ[] = [...this.stations.map((q) => [q.faceX + q.nx * 0.3, q.faceZ + q.nz * 0.3] as XZ), ...spotPts];
    const RW = 0.28;
    const RC = 0.36;
    const { nx: NX, nz: NZ, x0: GX, z0: GZ, cell: GC } = W0;
    const queue = new Int32Array(NX * NZ);
    const labels = (blocks: readonly Rect[], carry: boolean): Int32Array => {
      const lab = new Int32Array(NX * NZ).fill(-1);
      const R = carry ? RC : RW;
      const free = (c: number) => {
        if (W0.clr[c]! < RW || (carry && T0.clr[c]! < RC)) return false;
        const x = GX + ((c % NX) + 0.5) * GC;
        const z = GZ + (((c / NX) | 0) + 0.5) * GC;
        for (const b of blocks) if (rectDistance(b, x, z) < R) return false;
        return true;
      };
      let next = 0;
      for (let c0 = 0; c0 < NX * NZ; c0++) {
        if (lab[c0] !== -1 || !free(c0)) continue;
        let h = 0;
        let t = 0;
        queue[t++] = c0;
        lab[c0] = next;
        while (h < t) {
          const c = queue[h++]!;
          const i = c % NX;
          const j = (c / NX) | 0;
          for (let q = 0; q < 4; q++) {
            const ii = i + (q === 0 ? 1 : q === 1 ? -1 : 0);
            const jj = j + (q === 2 ? 1 : q === 3 ? -1 : 0);
            if (ii < 0 || jj < 0 || ii >= NX || jj >= NZ) continue;
            const n = jj * NX + ii;
            if (lab[n] !== -1 || !free(n)) continue;
            lab[n] = next;
            queue[t++] = n;
          }
        }
        next++;
      }
      return lab;
    };
    const cellOf = ([x, z]: XZ) => {
      const i = Math.floor((x - GX) / GC);
      const j = Math.floor((z - GZ) / GC);
      return i < 0 || j < 0 || i >= NX || j >= NZ ? -1 : j * NX + i;
    };
    const baseWalk = labels([], false);
    const baseCarry = labels([], true);
    const splits = (pts: readonly XZ[], base: Int32Array, lab: Int32Array): boolean => {
      const seen = new Map<number, number>();
      for (const p of pts) {
        const c = cellOf(p);
        if (c < 0 || base[c]! < 0) continue;
        const now = lab[c]!;
        if (now < 0) return true;
        const prev = seen.get(base[c]!);
        if (prev === undefined) seen.set(base[c]!, now);
        else if (prev !== now) return true;
      }
      return false;
    };
    const cutsOff = (blocks: readonly Rect[]): boolean => splits(walkPts, baseWalk, labels(blocks, false)) || splits(spotPts, baseCarry, labels(blocks, true));
    // The docks of the racks already taken (each may stand at any of them).
    const docks: Rect[] = [];
    const pushed: Fx[] = [];
    const LENS = [2.2, 1.8, 1.4, 2.6, 3.2];
    for (const fx of K.fixtures) {
      if (!MOVABLE.test(fx.id) || !fx.rect || pushed.length >= 2) continue;
      const rect = fx.rect;
      const home = new THREE.Vector3((rect.minX + rect.maxX) / 2, 0, (rect.minZ + rect.maxZ) / 2);
      const nf = frontNormal(fx);
      const half: XZ = [(rect.maxX - rect.minX) / 2, (rect.maxZ - rect.minZ) / 2];
      const hd = Math.abs(nf[0]) > 0.5 ? half[0] : half[1];
      const at = (x: number, z: number): Rect => ({ minX: x - half[0], maxX: x + half[0], minZ: z - half[1], maxZ: z + half[1] });
      // Every dock it can be rolled to along its axis (either way), clear of everything, with room for the pusher.
      const cands: { d: XZ; L: number; c: THREE.Vector3; blocks: boolean }[] = [];
      for (const sg of [-1, 1]) {
        const d: XZ = [nf[0] * sg, nf[1] * sg];
        let reach = -1;
        for (let s = 0; s <= 3.2 + 1e-6; s += 0.05) {
          const c = { x: home.x + d[0] * s, z: home.z + d[1] * s };
          // The pusher behind it, and room to walk round to push it back from the far end.
          if (!roomOK(at(c.x, c.z), s < 0.3 ? 0.02 : 0.12, rect) || clearAt(c.x - d[0] * (hd + 0.5), c.z - d[1] * (hd + 0.5), rect) < 0.34) break;
          reach = s;
        }
        for (const L of LENS) {
          if (L > reach + 1e-6) continue;
          if (clearAt(home.x + d[0] * (L + hd + 0.5), home.z + d[1] * (L + hd + 0.5), rect) < 0.34) continue;
          const c = home.clone().add(new THREE.Vector3(d[0] * L, 0, d[1] * L));
          if (stationNear(c.x, c.z) || this.spots.some((sp) => Math.hypot(sp.face[0] - c.x, sp.face[1] - c.z) < 1.0)) continue;
          cands.push({ d, L, c, blocks: cutsOff([...docks, at(c.x, c.z)]) });
        }
      }
      const open = cands.filter((q) => !q.blocks);
      if (!open.length) continue;
      const homeBlocks = cutsOff([...docks, rect]);
      const away = open[0]!;
      // Home in the way: it is pushed out first, then between two docks out of the way (or left out).
      let next: THREE.Vector3 | null = null;
      if (homeBlocks) {
        const other = open.find((q) => q.d[0] === away.d[0] && q.d[1] === away.d[1] && Math.abs(q.L - away.L) >= 0.6);
        next = other ? other.c : null;
      }
      pushed.push(fx);
      if (!homeBlocks) docks.push(rect);
      docks.push(at(away.c.x, away.c.z));
      if (next) docks.push(at(next.x, next.z));
      fx.model.hooks.set?.("pans", 9 + this.racks.length * 4);
      this.racks.push({ fx, key: `rack:${fx.f.fixture_id}`, home, away: away.c, d: away.d, n: nf, out: false, busy: false, swivel: 0, half, base: fx.model.group.position.clone(), origin: home.clone(), next, stay: homeBlocks && !next });
      // (One in the way is moved soon.)
      if (homeBlocks) this.nextRack = Math.min(this.nextRack, 3);
    }
    if (!pushed.length) return [K.grid, K.tall];
    // The racks that roll are parked rectangles in the crowd, not fixed obstacles, for walking and for a pan.
    const moving = (r: Rect) => pushed.some((fx) => fx.rect === r);
    return [new NavGrid(K.grid.bounds, rects.filter((r) => !moving(r)), K.grid.cell), new NavGrid(K.tall.bounds, K.tall.rects.filter((r) => !moving(r)), K.tall.cell)];
  }

  /** A charbroiler's steaks on its grates and the tongs on its front ledge. */
  private setupGrill(s: Station): void {
    const yaw = Math.atan2(-s.nx, -s.nz);
    const right = bodyDir(yaw, -1, 0, 0);
    const n = new THREE.Vector3(s.nx, 0, s.nz);
    const seared = new THREE.MeshStandardMaterial({ color: 0x4f2616, roughness: 0.62 });
    const raw = new THREE.MeshStandardMaterial({ color: 0x9a3c2c, roughness: 0.5 });
    this.ownMats.add(seared).add(raw);
    s.steaks = [-0.2, 0, 0.2].map((ox, i) => {
      // A ribeye 165 x 115 x 32 mm: its seared face down on the grate, the raw one up.
      const g = new THREE.Group();
      const geo = new THREE.SphereGeometry(0.5, 20, 10);
      const top = new THREE.Mesh(geo, raw);
      top.scale.set(0.165, 0.032, 0.115);
      const bot = new THREE.Mesh(geo, seared);
      bot.scale.set(0.163, 0.03, 0.113);
      bot.position.y = -0.002;
      top.castShadow = true;
      g.add(top, bot);
      const at = s.at.clone().addScaledVector(right, ox).addScaledVector(n, -0.03 + 0.05 * (i % 2)).setY(s.at.y + 0.016);
      g.position.copy(at);
      g.rotation.y = yaw + 0.25 * (i - 1);
      g.userData.seared = seared;
      this.group.add(g);
      return { obj: g, at, up: 0, done: 0 };
    });
    const m = buildFixture("NSF-TONGS-12IN", { w: 0.055, d: 0.305, h: 0.038 }, 1, true);
    const grip = new THREE.Group();
    grip.name = "tool-tongs";
    // The tongs' hinge end in the hand, the cups 0.185 m on along the grip axis (+z).
    m.group.rotation.y = Math.PI;
    m.group.position.set(0.0275, -0.019, -0.09);
    grip.add(m.group);
    grip.traverse((o) => {
      if ((o as THREE.Mesh).isMesh) o.castShadow = true;
    });
    // At rest on the front ledge to the right, the cups towards the back.
    const rest: Rest = { p: s.at.clone().addScaledVector(right, 0.36).addScaledVector(n, 0.2).setY(s.at.y + 0.021), z: n.clone().negate(), y: UP.clone() };
    restMatrix(rest, grip.matrix).decompose(grip.position, grip.quaternion, grip.scale);
    this.group.add(grip);
    s.tongs = { obj: grip, model: m, rest };
  }

  /** Where a person stands to pick from or place on a spot: as far back as the pan stays in reach. */
  private standFor(a: Agent, sp: Spot): XZ | null {
    const key = `${a.id}|${sp.key}`;
    if (!this.standCache.has(key)) {
      const it = this.standSteps(a, sp);
      while (!it.next().done);
    }
    return this.standCache.get(key)!;
  }
  /** `standFor` worked out a trial stance at a time (each poses the body), yielding between them; cached. */
  private *standSteps(a: Agent, sp: Spot): Generator<void, void> {
    const key = `${a.id}|${sp.key}`;
    if (this.standCache.has(key)) return;
    let out: XZ | null = null;
    const g = a.ghost;
    for (let dp = 0.44; dp >= 0.165; dp -= 0.025) {
      const x = sp.face[0] + sp.n[0] * dp;
      const z = sp.face[1] + sp.n[1] * dp;
      if (this.crowd.staticRoom(a.id, x, z, sp.yaw, SH_STAND) < STATIC_MARGIN + 0.004) continue;
      g.place(x, z, sp.yaw);
      const short = pickShortfall(g, { at: sp.at, surfaceY: sp.surfaceY, panTurn: sp.turn });
      if (short <= 0.004) {
        out = [x, z];
        break;
      }
      yield;
    }
    this.standCache.set(key, out);
  }

  /** Where a person stands to work a station, facing it; null if they cannot stand or reach there. */
  stance(a: Agent, s: Station): { x: number; z: number; yaw: number } | null {
    const x = s.faceX + s.nx * a.dFace;
    const z = s.faceZ + s.nz * a.dFace;
    if (this.crowd.grid.clearance(x, z) < a.rPath) return null;
    const reach = Math.hypot(s.at.x - x, s.at.z - z);
    if (reach > REACH[s.kind] * a.st.B.s) return null;
    return { x, z, yaw: Math.atan2(-s.nx, -s.nz) };
  }

  private workStep(s: Station, dur: number): RoutineStep {
    const opts: WorkOpts = { at: s.at.clone(), surfaceY: s.at.y };
    if (s.kind === "stir" || s.kind === "scrub") opts.tempo = s.kind === "scrub" ? 1.5 : 0.9 + 0.2 * this.rnd();
    return { kind: "work", work: WORK_OF[s.kind], duration: dur, opts };
  }

  private begin(a: Agent, start: { x: number; z: number; yaw: number }, steps: RoutineStep[], t0: number, s: Station | null, run: Run | null, workDur: number, path: XZ[] | null): void {
    a.start = start;
    a.steps = steps;
    a.t0 = t0;
    a.routine = new Routine(a.st, start, steps, false);
    a.station = s;
    a.run = run;
    a.workDur = workDur;
    a.path = path;
    a.toolK = 0;
    a.holding = false;
    if (s?.kind === "sala") {
      const sp = this.stance(a, s);
      a.lean = sp ? this.leanFor(a, sp, s) : 0;
    }
    if (s?.kind === "fry" && s.baskets) {
      // Drop the basket that is up, or lift the one that has cooked longest.
      const up = s.baskets.findIndex((b) => !b.down);
      s.fryJob = up >= 0 ? { k: up, lift: false } : { k: s.baskets[0]!.since < s.baskets[1]!.since ? 0 : 1, lift: true };
    }
  }

  /** Sample a routine's root on the ghost body every CROWD_DT from its start, with the footprint each sample shows. */
  private realize(a: Agent, start: { x: number; z: number; yaw: number }, steps: RoutineStep[], workDur: number): Motion {
    const r = new Routine(a.ghost, start, steps, false);
    const root = rootSampler(r);
    const n = Math.ceil(r.duration / CROWD_DT) + 1;
    const xs = new Float32Array(n);
    const zs = new Float32Array(n);
    const ys = new Float32Array(n);
    const sh = new Uint8Array(n);
    for (let i = 0; i < n; i++) {
      const p = root(Math.min(i * CROWD_DT, r.duration));
      xs[i] = p.x;
      zs[i] = p.z;
      ys[i] = p.yaw;
      sh[i] = p.walking ? SH_WALK : SH_STAND;
    }
    void workDur;
    return { n, x: xs, z: zs, yaw: ys, shape: sh };
  }

  /** The time (seconds from the start) each step of a routine begins, and its end. */
  private stepTimes(a: Agent, start: { x: number; z: number; yaw: number }, steps: readonly RoutineStep[]): number[] {
    return new Routine(a.ghost, start, steps, false).stepTimes;
  }

  /**
   * A run's motion: the root sampled every CROWD_DT (only where the body moves: through work and standing it stays
   * where the step began), with its footprint: a pan's discs from the pick to the place wherever they are clear of
   * the tall fixtures (inside a cabinet's opening they are over its shelf, in nobody's way), a rack's while pushing.
   */
  private realizeRun(a: Agent, start: { x: number; z: number; yaw: number }, steps: readonly RoutineStep[], times: readonly number[], o: RunShape): Motion {
    const r = new Routine(a.ghost, start, steps, false);
    const root = rootSampler(r);
    const D = r.duration;
    const n = Math.ceil(D / CROWD_DT) + 1;
    const xs = new Float32Array(n);
    const zs = new Float32Array(n);
    const ys = new Float32Array(n);
    const sh = new Uint8Array(n);
    const last = steps.length - 1;
    let si = 0;
    for (let i = 0; i < n; i++) {
      const t = Math.min(i * CROWD_DT, D);
      while (si < last && t >= times[si + 1]!) si++;
      const p = root(t);
      xs[i] = p.x;
      zs[i] = p.z;
      ys[i] = p.yaw;
      const moving = p.walking;
      let shape = moving ? SH_WALK : SH_STAND;
      if (o.pan && si > o.pan[0] && si <= o.pan[1]) {
        // The pan's own footprint, but for the steps that take it off its shelf or put it on one (the shuffle
        // after the pick, the one before the place, the place itself): there it is over the shelf, in its opening.
        // (Held turned, it is in a combi's opening for the steps that slide it in or out.)
        const inside = o.inside?.includes(si) ?? false;
        const atShelf = si === o.pan[0] + 1 || si >= o.pan[1] - 1 || inside;
        const [still, walking] = o.turn ? [SH_TURN, SH_TURN_WALK] : [SH_PAN, SH_PAN_WALK];
        if (moving && !inside) shape = walking;
        else if (!atShelf || this.crowd.staticRoom(a.id, xs[i]!, zs[i]!, ys[i]!, still) >= STATIC_MARGIN) shape = still;
      }
      if (o.rack && si >= o.rack[0] && si <= o.rack[1]) shape = moving ? SH_RACK_WALK : SH_RACK;
      sh[i] = shape;
    }
    return { n, x: xs, z: zs, yaw: ys, shape: sh };
  }

  /* -------------------------------------------------------------------------------------------------------- */
  /* Choosing what to do next                                                                                  */

  private choose(a: Agent): Station | null {
    const cands: { s: Station; w: number }[] = [];
    for (const s of this.stations) {
      if (s === a.station || s.holder >= 0) continue;
      const sp = this.stance(a, s);
      if (!sp) continue;
      // Not where someone stands or is going to stand.
      if (
        this.agents.some((b) => {
          if (b === a) return false;
          const p = this.crowd.poseAt(b.id, this.crowd.end(b.id));
          return Math.hypot(p.x - sp.x, p.z - sp.z) < a.rWalk + b.rWalk + PAIR_MARGIN + 0.04;
        })
      )
        continue;
      // Nor where a rack is parked.
      if (this.racks.some((r) => rectDistance(rackRect(r, r.out ? r.away : r.home), sp.x, sp.z) < a.rPath + PAIR_MARGIN)) continue;
      let w = a.weights.get(s) ?? 0.05;
      // The fryer calls when a basket is up or has been down long enough.
      if (s.kind === "fry" && s.baskets) w *= s.baskets.some((b) => !b.down || b.since < this.now - 25) ? 2.5 : 0.3;
      if (s.kind === "flip" || s.kind === "sala") w *= 1.6;
      cands.push({ s, w });
    }
    if (!cands.length) return null;
    let r = this.rnd() * cands.reduce((q, c) => q + c.w, 0);
    for (const c of cands) if ((r -= c.w) <= 0) return c.s;
    return cands[cands.length - 1]!.s;
  }

  /** Runs under way or planned, besides person a's. */
  private runsBusy(a: Agent): number {
    return this.agents.filter((b) => b !== a && ((b.run && !b.run.settled) || b.planned?.run)).length;
  }

  /** Whether a run is due for person a (a forced one, a rack to push, a pan to move) and others leave room for it. */
  private runDue(a: Agent): boolean {
    const force = this.force;
    if (force) return true;
    if (this.runsBusy(a) >= (this.agents.length >= 6 ? 2 : 1)) return false;
    return (
      (this.now >= this.nextRack && this.racks.some((r) => !r.busy && !(r.stay && r.out))) ||
      (this.now >= this.nextPan && this.spots.length > 0) ||
      this.combis.some((c) => this.combiDue(c) !== null)
    );
  }

  /** The runs person a might make now, best first (a generator: each pan run's stances cost some posing). */
  private *runSeqs(a: Agent, here: { x: number; z: number }): Generator<void, RunSeq[]> {
    const force = this.force;
    const out: RunSeq[] = [];
    if (force === "rack" || (!force && this.now >= this.nextRack)) {
      const R = this.racks.find((r) => !r.busy && !(r.stay && r.out));
      if (R) {
        const seq = this.rackSeq(a, R);
        if (seq) out.push(seq);
      }
      if (force === "rack" || out.length) return out;
    }
    // A combi to load, or its pan done.
    if (!force || force === "combi") {
      for (const c of this.combis) {
        const due = this.combiDue(c, force === "combi");
        if (due === null) continue;
        yield* this.standSteps(a, c.stage);
        yield;
        const seq = this.combiSeq(a, c, due);
        if (seq) out.push(seq);
      }
      // (A pan run due as well comes after it: tried if the combi's cannot be made, so that a combi that cannot be
      // served now does not keep every pan run waiting.)
      if (force === "combi") return out;
    }
    if (!force && this.now < this.nextPan) return out;
    const dist = (sp: Spot) => Math.hypot(sp.face[0] - here.x, sp.face[1] - here.z);
    const kindOf = (sp: Spot): RunKind => (sp.walkin ? "walkin" : sp.door ? "reachin" : "shelf");
    const full = (sp: Spot) => sp.pan && sp.booked < 0;
    const empty = (sp: Spot) => !sp.pan && sp.booked < 0;
    const cooling = (from: Spot, to: Spot) => !force && (this.cool.get(`${from.key}>${to.key}`) ?? 0) > this.now;
    const pairs: { from: Spot; to: Spot; cost: number }[] = [];
    const add = (from: Spot, to: Spot) => {
      if (!cooling(from, to)) pairs.push({ from, to, cost: dist(from) + Math.hypot(from.face[0] - to.face[0], from.face[1] - to.face[1]) + 2 * this.rnd() });
    };
    // A pan used up on the line goes to the pot sink first; then a full one is fetched from the cold to the line.
    if (!force || force === "table")
      for (const from of this.spots)
        if (!from.cold && !from.dish && !from.turn && full(from) && (force === "table" || (from.usedAt ?? Infinity) <= this.now)) for (const to of this.spots) if (to.dish && empty(to)) add(from, to);
    if (!pairs.length && force !== "table")
      for (const from of this.spots) {
        if (!from.cold || !full(from) || (force && kindOf(from) !== force)) continue;
        for (const to of this.spots) if (!to.cold && !to.dish && !to.turn && empty(to)) add(from, to);
      }
    pairs.sort((p, q) => p.cost - q.cost);
    for (const { from, to } of pairs.slice(0, 6)) {
      if (out.length >= 3) break;
      yield* this.standSteps(a, from);
      yield* this.standSteps(a, to);
      yield;
      const seq = this.panSeq(a, from, to);
      if (seq) out.push(seq);
    }
    return out;
  }

  /** The steps of a pan run: to the cold store (its door opened first), the pan taken, carried, put down. */
  private panSeq(a: Agent, from: Spot, to: Spot): RunSeq | null {
    const steps: RoutineStep[] = [];
    const walks: RunSeq["walks"] = [];
    const cold = from.cold ? from : to;
    const seq: RunSeq = { kind: cold.walkin ? "walkin" : cold.door ? "reachin" : from.cold ? "shelf" : "table", steps, walks, pick: -1, place: -1, open: -1, closeAfter: -1, arm: 0, from, to };
    let cur: XZ = [0, 0];
    const panOpts = (sp: Spot): WorkOpts => ({ at: sp.at.clone(), surfaceY: sp.surfaceY });
    // (A walk to a counter or a cabinet comes in square to it, its last stretch along `arrive` reversed: a pan held
    // out ahead comes over the counter, not round past the cabinet beside it.)
    const walk = (to: XZ, carry: boolean, arrive?: XZ) => {
      walks.push({ i: steps.length, from: cur, to, carry, arrive });
      steps.push({ kind: "walk", path: [], speed: carry ? a.speed * 0.92 : a.speed, upper: carry ? "carry" : "free" });
      cur = to;
    };
    const turn = (yaw: number, carry: boolean) => steps.push({ kind: "turn", yaw, upper: carry ? "carry" : "free" });
    const shuffle = (to: XZ, carry: boolean) => {
      steps.push({ kind: "shuffle", to: [to[0], to[1]], upper: carry ? "carry" : "free" });
      cur = to;
    };
    const off = (p: XZ, n: XZ, d: number): XZ => [p[0] + n[0] * d, p[1] + n[1] * d];
    const yawOf = (n: XZ) => Math.atan2(-n[0], -n[1]);
    // Where the walk to a shelf or a counter ends, before the turn to it: in a walk-in's aisle, its middle; at the
    // line, close in, out of the way of the cooks' backs across the aisle.
    const appr = (sp: Spot) => (sp.cold ? 0.62 : 0.36);
    const P = (sp: Spot) => this.standFor(a, sp);
    // The door: stand beside its pull (clear of the leaf's swing), pull it open; the stance in front of the pull.
    const doorStance = (d: Door, n: XZ, out: number, along: number): XZ => {
      const h = handleAt(d, 0, _rv1);
      return [h.x + d.e0.x * along + n[0] * out, h.z + d.e0.z * along + n[1] * out];
    };
    const openDoor = (d: Door, n: XZ, along: number, out: number) => {
      const H = doorStance(d, n, out, along);
      walk(H, false);
      turn(yawOf(n), false);
      seq.open = steps.length;
      // The hand to the pull, the pull, and most of the swing on to the stop before stepping in.
      steps.push({ kind: "stand", duration: 0.45 + d.pull + 0.55 * d.swing + 0.2 });
      seq.door = d;
      // The hand nearer the pull.
      const h = handleAt(d, 0, _rv1);
      const yaw = yawOf(n);
      seq.arm = (h.x - H[0]) * Math.cos(yaw) - (h.z - H[1]) * Math.sin(yaw) > 0 ? 0 : 1;
    };
    // Taking the pan from an open place (or a cabinet whose door is open, from where the opener stands).
    const take = (sp: Spot, already: boolean) => {
      const p = P(sp);
      if (!p) return false;
      if (sp.walkin) {
        const B = off(sp.face, sp.n, appr(sp));
        walk(B, false);
        turn(sp.yaw, false);
      } else if (!already) {
        walk(off(sp.face, sp.n, appr(sp)), false, sp.n);
        turn(sp.yaw, false);
      }
      shuffle(p, false);
      seq.pick = steps.length;
      steps.push({ kind: "work", work: "pick", duration: PICK_S, opts: panOpts(sp) });
      return true;
    };
    const backOff = (sp: Spot) => {
      if (sp.door) shuffle(doorStance(sp.door, sp.n, 0.75, 0.36), true);
      else {
        // From shelves, far enough out that the pan held ahead clears them as they turn away (if the aisle allows).
        const far = sp.cold ? off(sp.face, sp.n, 0.72) : null;
        shuffle(far && this.crowd.grid.clearance(far[0], far[1]) >= a.rWalk + 0.06 ? far : off(sp.face, sp.n, appr(sp)), true);
      }
    };
    const put = (sp: Spot, already: boolean) => {
      const p = P(sp);
      if (!p) return false;
      if (sp.door && already) {
        walk(doorStance(sp.door, sp.n, 0.75, 0.36), true, sp.n);
        turn(sp.yaw, true);
      } else {
        walk(off(sp.face, sp.n, appr(sp)), true, sp.walkin ? undefined : sp.n);
        turn(sp.yaw, true);
      }
      shuffle(p, true);
      seq.place = steps.length;
      steps.push({ kind: "work", work: "place", duration: PLACE_S, opts: panOpts(sp) });
      return true;
    };
    if (from.cold) {
      if (from.walkin) {
        openDoor(from.walkin.door, doorNormal(from.walkin.door), 0.3, 0.5);
        seq.walkin = from.walkin;
        if (!take(from, true)) return null;
      } else if (from.door) {
        openDoor(from.door, from.n, 0.27, 0.45);
        if (!take(from, true)) return null;
      } else if (!take(from, false)) return null;
      seq.closeAfter = seq.pick;
      backOff(from);
      // A reach-in's door swings shut on its closer behind them (they wait a moment out of its way, the pan held)
      // before they carry the pan off past where its leaf stood.
      if (from.door) {
        steps.push({ kind: "work", work: "carry", duration: from.door.close + 0.35 });
        seq.shut = steps.length;
      }
      if (!put(to, false)) return null;
      steps.push({ kind: "stand", duration: 1.2 });
    } else {
      // Taking a pan back: a door is opened first (both hands will be on the pan), then the pan fetched.
      if (to.walkin) {
        openDoor(to.walkin.door, doorNormal(to.walkin.door), 0.3, 0.5);
        seq.walkin = to.walkin;
      } else if (to.door) openDoor(to.door, to.n, 0.27, 0.45);
      if (!take(from, false)) return null;
      backOff(from);
      if (!put(to, !!to.door)) return null;
      seq.closeAfter = seq.place;
      if (to.walkin) {
        shuffle(off(to.face, to.n, appr(to)), false);
        walk(doorStance(to.walkin.door, doorNormal(to.walkin.door), 0.95, 0.5), false);
      } else if (to.door) shuffle(doorStance(to.door, to.n, 0.75, 0.36), false);
      steps.push({ kind: "stand", duration: 1.2 });
    }
    return seq;
  }

  /** Whether a combi is due: its pan to go in ("load"), or done and to come out ("unload"); null if not. */
  private combiDue(c: Combi, force = false): "load" | "unload" | null {
    if (c.spot.booked >= 0 || c.stage.booked >= 0) return null;
    if (c.stage.pan && !c.spot.pan && (force || this.now >= c.loadAt)) return "load";
    if (c.spot.pan && !c.stage.pan && (force || this.now >= c.doneAt)) return "unload";
    return null;
  }

  /**
   * The steps of loading a combi, or unloading it: its door pulled open first (both hands will be on the pan); the
   * pan taken off the counter turned, carried square up to the oven, raised to the rail, slid in over the runners
   * and pushed home by its near end (or drawn out by it, and carried back); then round the open door to its pull,
   * and the door swung shut by hand.
   */
  private combiSeq(a: Agent, c: Combi, due: "load" | "unload"): RunSeq | null {
    const load = due === "load";
    const from = load ? c.stage : c.spot;
    const to = load ? c.spot : c.stage;
    const steps: RoutineStep[] = [];
    const walks: RunSeq["walks"] = [];
    const inside: number[] = [];
    const seq: RunSeq = { kind: "combi", steps, walks, pick: -1, place: -1, open: -1, closeAfter: -1, arm: 0, from, to, combi: c, inside };
    const st = a.st;
    const B = st.B;
    const grid = this.crowd.grid;
    let cur: XZ = [0, 0];
    const n = c.spot.n;
    const yawO = c.spot.yaw;
    const off = (p: XZ, nn: XZ, d: number): XZ => [p[0] + nn[0] * d, p[1] + nn[1] * d];
    // The rim's height: carried (a little below the waist, or just clear of the counter), and slid in or out (just
    // over the rail's runners).
    const carryY = Math.max(c.stage.at.y + 0.02, B.waistY - 0.03 * B.s);
    const railIn = c.holdY + 0.012;
    const carry = (y: number): WorkOpts => ({ panTurn: true, holdY: y });
    const walk = (p: XZ, o: WorkOpts | null, arrive?: XZ) => {
      walks.push({ i: steps.length, from: cur, to: p, carry: !!o, arrive });
      steps.push(o ? { kind: "walk", path: [], speed: a.speed * 0.9, upper: "carry", opts: o } : { kind: "walk", path: [], speed: a.speed, upper: "free" });
      cur = p;
    };
    const turn = (yaw: number, o: WorkOpts | null) => steps.push(o ? { kind: "turn", yaw, upper: "carry", opts: o } : { kind: "turn", yaw, upper: "free" });
    const shuffle = (p: XZ, o: WorkOpts | null) => {
      steps.push(o ? { kind: "shuffle", to: [p[0], p[1]], upper: "carry", opts: o } : { kind: "shuffle", to: [p[0], p[1]], upper: "free" });
      cur = p;
    };
    // Before the oven: where the pan's far end is clear of its front (F), and where it is slid home from (N): the
    // arms take it 12 cm on from the hold, and the push the rest.
    const hz = panHoldZ(st, railIn);
    let dFar = hz + HOTEL_PAN.length / 2 + 0.05 + STATIC_MARGIN;
    const dNear = hz + 0.12 + c.push + c.zc;
    // The pan is brought to F from the side away from the hinge (the open leaf stands out on the other): turned to
    // the oven a little way along the front (S), then sidestepped across, clear of the leaf. In a narrow aisle the two
    // stand up to three centimetres further out or further in: wherever whoever holds the pan, finishing the turn to
    // the oven, has most room both from the oven's front and from the counter at their back.
    const e0: XZ = [c.door.e0.x, c.door.e0.z];
    const slack = (d: number) => {
      const P = off(c.face, n, d);
      let q = Infinity;
      for (const p of [P, off(P, e0, 0.55)]) {
        q = Math.min(q, grid.clearance(p[0], p[1]) - a.rPath);
        for (const dy of [-0.12, 0, 0.12]) q = Math.min(q, this.crowd.staticRoom(a.id, p[0], p[1], yawO + dy, SH_TURN) - STATIC_MARGIN);
      }
      return q;
    };
    if (slack(dFar) < 0.004) {
      let best = dFar;
      let most = slack(dFar);
      for (let dd = -0.03; dd <= 0.0301; dd += 0.005) {
        const q = slack(dFar + dd);
        if (q > most + 1e-4) [best, most] = [dFar + dd, q];
      }
      dFar = best;
    }
    if (dNear < 0.22 || dFar - dNear < 0.05) return null;
    const F = off(c.face, n, dFar);
    const N = off(c.face, n, dNear);
    const S = off(F, e0, 0.55);
    // At the counter: square to it, out far enough that the pan held ahead clears its edge.
    const Pst = this.standFor(a, c.stage);
    if (!Pst) return null;
    const A = off(c.stage.face, c.stage.n, panHoldZ(st, carryY) + HOTEL_PAN.length / 2 + 0.05 + STATIC_MARGIN);
    // The door: beside its pull, the pull's hand on it.
    const h0 = handleAt(c.door, 0, _rv1);
    const O: XZ = [h0.x + c.door.e0.x * 0.27 + n[0] * 0.45, h0.z + c.door.e0.z * 0.27 + n[1] * 0.45];
    for (const q of [F, S, A, O]) if (grid.clearance(q[0], q[1]) < a.rWalk + 0.02) return null;
    // Held out ahead, the pan sweeps round with whoever turns on the spot: the way round (split into turns of under
    // half a revolution) that keeps it and them furthest from the open leaf and the room; none if both come too near.
    const room = (x: number, z: number, yaw: number) => {
      let q = this.crowd.staticRoom(a.id, x, z, yaw, SH_TURN) - STATIC_MARGIN;
      const sn = Math.sin(yaw);
      const co = Math.cos(yaw);
      for (const ci of a.shapes[SH_TURN] ?? []) q = Math.min(q, leafGap(c.door, 1, x + co * ci.x + sn * ci.z, z - sn * ci.x + co * ci.z) - ci.r - 0.015);
      return q;
    };
    const turnRound = (p: XZ, y0: number, y1: number, o: WorkOpts): boolean => {
      const d0 = Math.atan2(Math.sin(y1 - y0), Math.cos(y1 - y0));
      let best = d0;
      let most = -Infinity;
      for (const w of [d0, d0 > 0 ? d0 - 2 * Math.PI : d0 + 2 * Math.PI]) {
        let q = Infinity;
        const m = Math.max(2, Math.ceil(Math.abs(w) / 0.1));
        for (let j = 0; j <= m; j++) q = Math.min(q, room(p[0], p[1], y0 + (w * j) / m));
        if (q > most) [best, most] = [w, q];
      }
      if (most < 0.005) return false;
      const k = Math.ceil(Math.abs(best) / (0.9 * Math.PI));
      for (let j = 1; j <= k; j++) turn(y0 + (best * j) / k, o);
      return true;
    };
    walk(O, null);
    turn(yawO, null);
    seq.open = steps.length;
    steps.push({ kind: "stand", duration: 0.45 + c.door.pull + 0.55 * c.door.swing + 0.2 });
    seq.door = c.door;
    seq.arm = (h0.x - O[0]) * Math.cos(yawO) - (h0.z - O[1]) * Math.sin(yawO) > 0 ? 0 : 1;
    const level = Math.abs(railIn - carryY) < 0.008;
    if (load) {
      walk(A, null, c.stage.n);
      turn(c.stage.yaw, null);
      shuffle(Pst, null);
      seq.pick = steps.length;
      const pk: WorkOpts = { panTurn: true, at: c.stage.at.clone(), surfaceY: c.stage.surfaceY, holdY: carryY };
      steps.push({ kind: "work", work: "pick", duration: panSeconds("pick", pk), opts: pk });
      shuffle(A, carry(carryY));
      walk(S, carry(carryY));
      turn(yawO, carry(carryY));
      if (!level) steps.push({ kind: "work", work: "carry", duration: 1.1, opts: { panTurn: true, holdY: railIn, fromY: carryY, raise: "up" } });
      shuffle(F, carry(railIn));
      inside.push(steps.length);
      shuffle(N, carry(railIn));
      seq.place = steps.length;
      inside.push(steps.length);
      const pl: WorkOpts = { panTurn: true, at: c.spot.at.clone(), surfaceY: c.spot.surfaceY, push: c.push, holdY: railIn };
      steps.push({ kind: "work", work: "place", duration: panSeconds("place", pl), opts: pl });
      seq.closeAfter = seq.place;
      shuffle(off(c.face, n, Math.max(dNear + 0.25, 0.6)), null);
    } else {
      shuffle(N, null);
      seq.pick = steps.length;
      const pk: WorkOpts = { panTurn: true, at: c.spot.at.clone(), surfaceY: c.spot.surfaceY, push: c.push, holdY: railIn };
      steps.push({ kind: "work", work: "pick", duration: panSeconds("pick", pk), opts: pk });
      inside.push(steps.length);
      shuffle(F, carry(railIn));
      shuffle(S, carry(railIn));
      if (!level) steps.push({ kind: "work", work: "carry", duration: 1.1, opts: { panTurn: true, holdY: railIn, fromY: carryY, raise: "down" } });
      // With the counter a step or two away across the aisle: turned round to it there, away from the open leaf, and
      // stepped across to it. Otherwise walked to it and turned to it before it.
      const across = Math.hypot(Pst[0] - S[0], Pst[1] - S[1]);
      let clear = across < 0.7;
      for (let k = 0; clear && k <= 8; k++) clear = room(S[0] + ((Pst[0] - S[0]) * k) / 8, S[1] + ((Pst[1] - S[1]) * k) / 8, c.stage.yaw) >= 0.005;
      if (!clear || !turnRound(S, yawO, c.stage.yaw, carry(carryY))) {
        walk(A, carry(carryY), c.stage.n);
        turn(c.stage.yaw, carry(carryY));
      }
      shuffle(Pst, carry(carryY));
      seq.place = steps.length;
      const pl: WorkOpts = { panTurn: true, at: c.stage.at.clone(), surfaceY: c.stage.surfaceY, holdY: carryY };
      steps.push({ kind: "work", work: "place", duration: panSeconds("place", pl), opts: pl });
      seq.closeAfter = seq.place;
    }
    // Round the open door to its pull, and shut it: the hand takes it round as far as it reaches and lets it swing on.
    walk(c.close, null);
    turn(c.closeYaw, null);
    seq.closeAt = steps.length;
    const sh = handShut(c.fRel);
    steps.push({ kind: "stand", duration: 0.5 + sh.hand + sh.swing + 0.3 });
    seq.shut = steps.length;
    return seq;
  }

  /** The steps of a rack push: behind it, hands on its uprights, pushed straight to its other place, let go. */
  private rackSeq(a: Agent, R: Rack): RunSeq | null {
    const from = R.out ? R.away : R.home;
    const to = R.out ? R.home : R.away;
    const L = Math.hypot(to.x - from.x, to.z - from.z);
    if (L < 0.3) return null;
    const d: XZ = [(to.x - from.x) / L, (to.z - from.z) / L];
    const hd = Math.abs(R.n[0]) > 0.5 ? R.half[0] : R.half[1];
    const A: XZ = [from.x - d[0] * (hd + a.rackFace), from.z - d[1] * (hd + a.rackFace)];
    const yaw = Math.atan2(d[0], d[1]);
    const steps: RoutineStep[] = [];
    const walks: RunSeq["walks"] = [{ i: 0, from: [0, 0], to: A, carry: false }];
    steps.push({ kind: "walk", path: [], speed: a.speed, upper: "free" });
    steps.push({ kind: "turn", yaw, upper: "free" });
    const grab = steps.length;
    // Hands round the near face's two corner posts (1 in tube, their centres 12.7 mm in from the sides and 15.9 mm
    // from the face), thumbs up.
    const opts: WorkOpts = { surfaceY: a.gripY, uprights: true, halfWidth: RACK_DIMS.w / 2 - 0.0127, reach: a.rackFace + 0.0159 };
    steps.push({ kind: "work", work: "push", duration: 1.1, opts });
    const walkI = steps.length;
    steps.push({ kind: "walk", path: [[A[0] + d[0] * L, A[1] + d[1] * L]], speed: RACK_SPEED, upper: "push", opts });
    const release = steps.length;
    steps.push({ kind: "work", work: "push", duration: 0.9, opts });
    steps.push({ kind: "stand", duration: 1.2 });
    // The near face is the front when the pusher stands in front of it (pushing it backwards).
    const front = -(d[0] * R.n[0] + d[1] * R.n[1]) > 0;
    return { kind: "rack", steps, walks, pick: -1, place: -1, open: -1, closeAfter: -1, arm: 0, rack: { r: R, from, to, A, d, grab, walk: walkI, release, front } };
  }

  job: {
    a: Agent;
    s: Station | null;
    seq: RunSeq | null;
    k0: number;
    here: { x: number; z: number; yaw: number };
    dur: number;
    cands: Cand[];
    /** Finding a run and its paths (spread over frames), then checking its candidates against the crowd. */
    prep: Generator<void, { seq: RunSeq | null; cands: Cand[] } | null> | null;
    it: Generator<void, LegVerdict> | null;
  } | null = null;

  /** A run for person a: the first of the due runs whose walks all have paths, with its candidates. */
  private *prepRun(a: Agent, here: { x: number; z: number; yaw: number }, k0: number): Generator<void, { seq: RunSeq; cands: Cand[] } | null> {
    const t0ms = performance.now();
    const seqs = yield* this.runSeqs(a, here);
    for (const seq of seqs) {
      yield;
      const cands = yield* this.runCands(a, seq, here, k0);
      if (cands.length) return { seq, cands };
      this.runFailed(seq);
    }
    void t0ms;
    return null;
  }

  /**
   * Plan a person's next leg: a run if one is due and can be made, else the next station, by the first path
   * variant, speed and wait that keep clear. The path search runs now; making the candidates' motions and checking
   * them against everyone's committed tracks runs in slices (`update`).
   */
  private plan(a: Agent, runs = true): void {
    const t0ms = performance.now();
    const id = a.id;
    const k0 = Math.max(this.crowd.end(id) + 1, this.crowd.step(this.now + 1));
    const here = this.crowd.poseAt(id, this.crowd.end(id));
    if (runs && this.runDue(a)) {
      this.job = { a, s: null, seq: null, k0, here, dur: 1.2, cands: [], prep: this.prepRun(a, here, k0), it: null };
      return;
    }
    const s = this.choose(a);
    if (!s) {
      this.extend(a, 2);
      return;
    }
    const dur = DWELL[s.kind][0] + this.rnd() * (DWELL[s.kind][1] - DWELL[s.kind][0]);
    this.job = { a, s, seq: null, k0, here, dur, cands: [], prep: this.prepStation(a, s, here, k0, dur), it: null };
    void t0ms;
  }

  /** A leg to a station: the paths there (keeping right, or nearer the middle), each at two speeds. */
  private *prepStation(a: Agent, s: Station, here: { x: number; z: number; yaw: number }, k0: number, dur: number): Generator<void, { seq: null; cands: Cand[] } | null> {
    const goal = this.stance(a, s)!;
    const dyn = this.crowd.standingObstacles(a.id, k0);
    const cands: Cand[] = [];
    const paths: XZ[][] = [];
    let variant = 0;
    const t0ms = performance.now();
    const found = yield* this.crowd.grid.findPathsSteps([here.x, here.z], [goal.x, goal.z], a.rPath, [0.75, 1, 0.4], dyn);
    this.log.searchMax = Math.max(this.log.searchMax, performance.now() - t0ms);
    for (const path of found) {
      variant++;
      if (!path || path.length < 2) continue;
      if (paths.some((q) => q.length === path.length && q.every((pt, i) => Math.abs(pt[0] - path[i]![0]) < 1e-3 && Math.abs(pt[1] - path[i]![1]) < 1e-3))) continue;
      paths.push(path);
      for (const v of [a.speed, a.speed * 0.78]) {
        const steps: RoutineStep[] = [{ kind: "walk", path: path.slice(1).map(([x, z]) => [x, z] as [number, number]), speed: v }, { kind: "turn", yaw: goal.yaw }, this.workStep(s, dur)];
        const start = { x: here.x, z: here.z, yaw: here.yaw };
        cands.push({ speed: v, variant, steps, path, run: null, closing: dur, realize: () => this.realize(a, start, steps, dur) });
      }
    }
    return cands.length ? { seq: null, cands } : null;
  }

  /** A run that could not be made: not tried again for a while (longer when no way was found than when it met someone). */
  private runFailed(seq: RunSeq, blocked = true): void {
    this.log.runFails[seq.kind] = (this.log.runFails[seq.kind] ?? 0) + 1;
    if (seq.kind === "rack") this.nextRack = this.now + 8;
    else if (seq.combi) {
      // Tried again a little later, but not for ever if the counter is too far off for its door.
      const c = seq.combi;
      if (seq.long) c.far++;
      const later = c.far >= 3 ? Infinity : this.now + (6 + 6 * this.rnd()) * (1 + c.far);
      if (c.far >= 3) this.log.whyNot["combi: its counter too far for its door"] = (this.log.whyNot["combi: its counter too far for its door"] ?? 0) + 1;
      if (seq.to === c.spot) c.loadAt = later;
      else c.doneAt = later;
    } else this.nextPan = this.now + 5;
    if (seq.from && seq.to) this.cool.set(`${seq.from.key}>${seq.to.key}`, this.now + (blocked ? 30 + 30 * this.rnd() : 8));
  }

  /**
   * A run's candidates: its walks given paths (keeping to the right, or nearer the middle; a carried pan kept off
   * the tall fixtures), at two speeds; each made into a motion when the planner gets to it, with its door's timing,
   * its gate window, the rack it takes over and leaves parked, all read off the motion itself.
   */
  private *runCands(a: Agent, seq: RunSeq, here: { x: number; z: number; yaw: number }, k0: number): Generator<void, Cand[]> {
    const id = a.id;
    const tally = (why: string) => {
      const k = `${seq.kind}: ${why}`;
      this.log.whyNot[k] = (this.log.whyNot[k] ?? 0) + 1;
    };
    const ignore = [...(seq.walkin ? [seq.walkin.key] : []), ...(seq.rack ? [seq.rack.r.key] : [])];
    const dyn = this.crowd.standingObstacles(id, k0, 10, ignore);
    // A pan carried at the waist stands out to the sides past the body: kept further off the tall fixtures and the
    // parked racks (it passes over counters).
    const tall = this.crowd.tall;
    const parked = this.racks.filter((q) => q !== seq.rack?.r).map((q) => rackRect(q, q.out ? q.away : q.home));
    const grid = this.crowd.grid;
    const fracs = [0.75, 0.3];
    const door = seq.door;
    const start = { x: here.x, z: here.z, yaw: here.yaw };
    // Paths that keep the clearance the footprint needs on straight runs can still bring it too near a fixture where
    // the walk rounds a corner (a pan out ahead swings into it): if every candidate does, search again keeping
    // further off, twice.
    for (const bump of [0, 0.07, 0.15]) {
      // The pan's discs reach 0.315 m to either side of the body's line, and keep the static margin.
      const extra = Math.max(0, 0.315 + STATIC_MARGIN + 0.005 - a.rPath) + bump;
      const dynCarry = (x: number, z: number) => {
        let q = Math.min(dyn(x, z), tall.clearanceLB(x, z) - extra);
        for (const r of parked) q = Math.min(q, rectDistance(r, x, z) - extra - PAIR_MARGIN);
        return q;
      };
      const dynWalk = (x: number, z: number) => Math.min(dyn(x, z), grid.clearance(x, z) - bump);
      // (The open leaf kept a little further off: a walk rounds its corners inside the path's.)
      const dynOpen = door ? (x: number, z: number) => Math.min(dynWalk(x, z), leafGap(door, 1, x, z) - 0.06 - bump) : dynWalk;
      const dynCarryOpen = door ? (x: number, z: number) => Math.min(dynCarry(x, z), leafGap(door, 1, x, z) - 0.06 - extra) : dynCarry;
      const found: (XZ[] | null)[][] = [];
      let from: XZ = [here.x, here.z];
      for (const w of seq.walks) {
        const f = w.i === 0 ? from : w.from;
        const open = door && w.i > seq.open && w.i < (seq.shut ?? Infinity);
        const t0ms = performance.now();
        // The walk's ends stay where they are (a stance, a door): only its way between keeps further off.
        const dd = w.carry ? (open ? dynCarryOpen : dynCarry) : open ? dynOpen : dynWalk;
        const plain = (x: number, z: number) => Math.min(dyn(x, z), grid.clearance(x, z), open && door ? leafGap(door, 1, x, z) - 0.06 : Infinity);
        // Coming in square: the search ends a little out from the stance, and the last stretch is straight in.
        let goal = w.to;
        let tail: XZ | null = null;
        if (w.arrive) {
          for (const L of [0.5, 0.36, 0.24, 0.15]) {
            const P: XZ = [w.to[0] + w.arrive[0] * L, w.to[1] + w.arrive[1] * L];
            if (plain(P[0], P[1]) >= a.rPath && grid.segmentClear(P[0], P[1], w.to[0], w.to[1], a.rPath, plain)) {
              goal = P;
              tail = w.to;
              break;
            }
          }
        }
        const ends = (x: number, z: number) => (Math.hypot(x - f[0], z - f[1]) < 0.3 || Math.hypot(x - goal[0], z - goal[1]) < 0.3 ? plain(x, z) : dd(x, z));
        let res = yield* grid.findPathsSteps(f, goal, a.rPath, fracs, ends);
        if (tail) {
          const T = tail;
          res = res.map((p) => (p ? [...p, T] : p));
        }
        this.log.searchMax = Math.max(this.log.searchMax, performance.now() - t0ms);
        if (res.every((p) => !p)) {
          tally(`no way ${w.carry ? "with the pan" : "there"}${open ? " past the open door" : ""}`);
          return [];
        }
        found.push(res);
        from = w.to;
        yield;
      }
      const cands: Cand[] = [];
      for (let v = 0; v < fracs.length; v++) {
        const steps = seq.steps.map((s) => ({ ...s })) as RoutineStep[];
        let path: XZ[] = [];
        let same = v > 0;
        seq.walks.forEach((w, j) => {
          const p = found[j]![v] ?? found[j]![0] ?? found[j]![1]!;
          if (found[j]![v] && v > 0 && found[j]![0] && JSON.stringify(found[j]![v]) !== JSON.stringify(found[j]![0])) same = false;
          const st = steps[w.i] as Extract<RoutineStep, { kind: "walk" }>;
          steps[w.i] = { ...st, path: p.slice(1).map(([x, z]) => [x, z] as [number, number]) };
          path = path.length ? [...path, ...p.slice(1)] : p.slice();
        });
        if (same) continue;
        for (const sp of v === 0 ? [1, 0.82] : [1]) {
          const st2 = sp === 1 ? steps : (steps.map((s) => (s.kind === "walk" && s.upper !== "push" ? { ...s, speed: (s.speed ?? 1) * sp } : s)) as RoutineStep[]);
          const c: Cand = { speed: a.speed * sp, variant: v + 1, steps: st2, path, run: null, closing: 1.2 };
          c.realize = () => this.makeRun(a, seq, start, st2, c);
          cands.push(c);
        }
      }
      // Keep the candidates whose motion clears the room (each made once more by the planner: it is cheap).
      const ok: Cand[] = [];
      for (const c of cands) {
        yield;
        const m = this.makeRun(a, seq, start, c.steps, c);
        if (m && this.crowd.motionClearOfRoom(id, m)) {
          c.motion = m;
          ok.push(c);
        }
      }
      if (ok.length) return ok;
      tally(cands.every((c) => c.run === null) ? "the door would meet them" : "too near a fixture");
    }
    return [];
  }

  /** Make a run candidate's motion and its run (null if the door would meet its opener). */
  private makeRun(a: Agent, seq: RunSeq, start: { x: number; z: number; yaw: number }, steps: RoutineStep[], c: Cand): Motion | null {
    const times = this.stepTimes(a, start, steps);
    const shape: RunShape = {};
    if (seq.pick >= 0 || seq.place >= 0) shape.pan = [seq.pick >= 0 ? seq.pick : -1, seq.place >= 0 ? seq.place : steps.length];
    if (seq.rack) shape.rack = [seq.rack.grab, seq.rack.release];
    if (seq.combi) {
      shape.turn = true;
      shape.inside = seq.inside;
    }
    const m = this.realizeRun(a, start, steps, times, shape);
    const dt = CROWD_DT;
    const idx = (t: number) => Math.max(0, Math.min(m.n - 1, Math.round(t / dt)));
    const run: Run = { kind: seq.kind, times, end: times[steps.length - 1]!, settled: false };
    const circles = (i: number, f: (x: number, z: number, r: number) => boolean) => {
      for (const ci of a.shapes[m.shape[i]!] ?? []) {
        const s = Math.sin(m.yaw[i]!);
        const co = Math.cos(m.yaw[i]!);
        if (f(m.x[i]! + co * ci.x + s * ci.z, m.z[i]! - s * ci.x + co * ci.z, ci.r)) return true;
      }
      return false;
    };
    m.timed = [];
    if (seq.walkin) {
      // The gate's window: from the first sample that touches its floor to the last.
      const g = seq.walkin.gate;
      let k0 = -1;
      let k1 = -1;
      for (let i = 0; i < m.n; i++)
        if (circles(i, (x, z, r) => rectDistance(g, x, z) < r + PAIR_MARGIN + 0.01)) {
          if (k0 < 0) k0 = i;
          k1 = i;
        }
      if (k0 >= 0) m.gates = [{ key: seq.walkin.key, k0: Math.max(0, k0 - 2), k1: k1 + 2 }];
      run.walkin = seq.walkin;
    }
    if (seq.door) {
      const d = seq.door;
      const t0 = times[seq.open]! + 0.45;
      const t1 = t0 + d.pull;
      const t2 = t1 + d.swing;
      // The leaf against the opener: their footprint, and the pan in their hands from the pick to the place.
      // The pan as it is (with the hands at its ends), not the crowd's discs over it: its outline against the
      // leaf, and the leaf's free edge against its inside.
      // (Held turned, a combi's pan: its long sides along the arms.)
      const panC = a.shapes[seq.combi ? SH_TURN : SH_PAN]!.filter((q) => q.high);
      const pz = panC.reduce((q, c) => q + c.z, 0) / Math.max(1, panC.length);
      const PX = seq.combi ? HOTEL_PAN.width / 2 + 0.045 : HOTEL_PAN.gripHalf + 0.035;
      const PZ = seq.combi ? HOTEL_PAN.length / 2 + 0.01 : HOTEL_PAN.width / 2 + 0.01;
      const outline: XZ[] = [[-PX, -PZ], [0, -PZ], [PX, -PZ], [PX, 0], [PX, PZ], [0, PZ], [-PX, PZ], [-PX, 0]];
      const panFrom = seq.pick >= 0 ? times[seq.pick + 1]! : 0;
      const panTo = seq.place >= 0 ? times[seq.place]! : Infinity;
      const hit = (i: number, f: number) => {
        if (f < 0.01) return false;
        const k = Math.min(i, m.n - 1);
        if (circles(k, (x, z, r) => leafGap(d, f, x, z) < r + 0.015)) return true;
        const t = k * dt;
        if (t < panFrom || t >= panTo) return false;
        const sn = Math.sin(m.yaw[k]!);
        const co = Math.cos(m.yaw[k]!);
        const x0 = m.x[k]!;
        const z0 = m.z[k]!;
        for (const [ox, oz] of outline) if (leafGap(d, f, x0 + co * ox + sn * (pz + oz), z0 - sn * ox + co * (pz + oz)) < 0.02) return true;
        leafDir(d, f, _rv3);
        const ex = d.hinge.x + _rv3.x * d.len - x0;
        const ez = d.hinge.z + _rv3.z * d.len - z0;
        const lx = co * ex - sn * ez;
        const lz = sn * ex + co * ez - pz;
        return Math.abs(lx) < PX + 0.02 && Math.abs(lz) < PZ + 0.02;
      };
      if (seq.combi && seq.closeAt !== undefined) {
        // Shut by hand, when they come round to it: clear of them (and of the pan) all the way, open to its close.
        const sh = handShut(seq.combi.fRel);
        const c0 = times[seq.closeAt]! + 0.5;
        const dp: DoorPlan = { d, t0, t1, t2, c0, c1: c0 + sh.hand + sh.swing, arm: seq.arm, rel: c0 + sh.hand, fRel: seq.combi.fRel, closeArm: seq.combi.closeArm, edge: seq.combi.closeEdge };
        for (let i = idx(t0); i <= Math.ceil(dp.c1 / dt); i++)
          if (hit(i, doorOpening(dp, i * dt))) {
            this.log.doorHits++;
            const k = `combi: its door would meet them (step ${times.findIndex((q) => q > i * dt) - 1}, open ${doorOpening(dp, i * dt).toFixed(2)})`;
            this.log.whyNot[k] = (this.log.whyNot[k] ?? 0) + 1;
            return null;
          }
        // Not left open while a pan is fetched from far off.
        if (dp.c1 - t0 > 45) {
          this.log.doorHits++;
          seq.long = true;
          const k = `combi: its door open ${10 * Math.round((dp.c1 - t0) / 10)} s`;
          this.log.whyNot[k] = (this.log.whyNot[k] ?? 0) + 1;
          return null;
        }
        run.door = dp;
        m.timed.push({ rect: d.rect, k0: idx(t0) - 2, k1: Math.ceil(dp.c1 / dt) + 2 });
        run.end = Math.max(run.end, dp.c1);
      } else {
        const opening: DoorPlan = { d, t0, t1, t2, c0: Infinity, c1: Infinity, arm: seq.arm };
        const nClose = Math.ceil(d.close / dt) + 1;
        // The first time the opener would meet the leaf if it stayed open.
        let firstHit = Infinity;
        for (let i = idx(t0); i < m.n + nClose; i++)
          if (hit(i, doorOpening(opening, i * dt))) {
            firstHit = i;
            break;
          }
        if (firstHit * dt < t2 + 0.05) {
          this.log.doorHits++;
          return null;
        }
        // It shuts behind them: as soon as they have done what it was opened for (and are out of a walk-in's
        // doorway), and its closing sweep stays clear of them all the way.
        let earliest = Math.max(t2, times[seq.closeAfter + 1]!);
        if (seq.walkin) {
          const g = seq.walkin.gate;
          let last = -1;
          for (let i = 0; i < m.n; i++) if (circles(i, (x, z, r) => rectDistance(g, x, z) < r + 0.05)) last = i;
          earliest = Math.max(earliest, last * dt + 0.3);
        }
        let c0 = -1;
        for (let c = earliest; idx(c) < firstHit && c < m.n * dt + 3; c += dt) {
          let ok = true;
          for (let i = idx(c); i <= idx(c + d.close) && ok; i++) if (hit(i, 1 - smoothstep(c, c + d.close, i * dt))) ok = false;
          if (ok) {
            c0 = c;
            break;
          }
        }
        if (c0 < 0) {
          this.log.doorHits++;
          return null;
        }
        const dp: DoorPlan = { d, t0, t1, t2, c0, c1: c0 + d.close, arm: seq.arm };
        // A reach-in is not left standing open while a pan is fetched from far off.
        if (!seq.walkin && seq.from && !seq.from.cold && dp.c1 - t0 > 16) {
          this.log.doorHits++;
          return null;
        }
        run.door = dp;
        m.timed.push({ rect: d.rect, k0: idx(t0) - 2, k1: idx(dp.c1) + 2 });
        run.end = Math.max(run.end, dp.c1);
      }
    }
    if (seq.from && seq.to && (seq.pick >= 0 || seq.place >= 0)) {
      run.pan = { obj: seq.from.pan!, from: seq.from, to: seq.to, t0: times[seq.pick]!, t1: times[seq.place + 1]!, held: false, placed: false };
      run.end = Math.max(run.end, run.pan.t1);
    }
    if (seq.rack) {
      const R = seq.rack;
      m.release = [{ key: R.r.key, k: idx(times[R.grab]!) }];
      m.timed.push({ rect: rackRect(R.r, R.to), k0: idx(times[R.release + 1]!), k1: Infinity, owner: -1, key: R.r.key });
      run.rack = {
        r: R.r,
        from: R.from,
        to: R.to,
        A: R.A,
        d: R.d,
        t0: times[R.grab]!,
        t1: times[R.release + 1]!,
        s0: R.r.swivel,
        s1: R.d[0] * R.r.n[0] + R.d[1] * R.r.n[1] > 0 ? 0 : Math.PI,
        parked: false,
      };
      run.end = Math.max(run.end, run.rack.t1);
    }
    c.run = run;
    return m;
  }

  /** The planning job's verdict: start the leg when its time comes, or keep working and try again shortly. */
  private finish(verdict: LegVerdict): void {
    const { a, s, seq, k0, here, dur, cands } = this.job!;
    this.job = null;
    if (!verdict.ok) {
      this.log.waits++;
      a.fails++;
      if (seq) {
        this.runFailed(seq, false);
        const k = `${seq.kind}: ${verdict.reason}`;
        this.log.whyNot[k] = (this.log.whyNot[k] ?? 0) + 1;
      }
      this.extend(a, 1.5);
      return;
    }
    const c = cands.find((q) => q.variant === verdict.variant && q.speed === verdict.speed)!;
    const delaySteps = Math.round(verdict.delay / CROWD_DT);
    const t0 = (k0 + delaySteps) * CROWD_DT;
    // The current work runs on until the leg starts.
    this.stretch(a, t0);
    a.planned = { steps: c.steps, start: { x: here.x, z: here.z, yaw: here.yaw }, t0, station: s, run: c.run, workDur: seq ? c.closing : dur, path: c.path };
    if (s) s.holder = a.id;
    if (seq && c.run) {
      this.log.runs[seq.kind] = (this.log.runs[seq.kind] ?? 0) + 1;
      if (c.run.pan) c.run.pan.from.booked = c.run.pan.to.booked = a.id;
      if (c.run.rack) c.run.rack.r.busy = true;
      if (seq.kind === "rack") this.nextRack = this.now + 40 + 25 * this.rnd();
      else if (seq.kind !== "combi") this.nextPan = this.now + 10 + 12 * this.rnd();
      if (this.force === seq.kind || (this.force && seq.kind !== "rack" && this.force !== "rack")) this.force = null;
    }
    a.fails = 0;
    this.log.legs++;
    if (delaySteps > 0) this.log.delays++;
  }

  /** Lengthen the current routine's closing work so that it ends at time t. */
  private stretch(a: Agent, t: number): void {
    const last = a.steps[a.steps.length - 1]!;
    if (last.kind !== "work" && last.kind !== "stand") return;
    const extra = t - (a.t0 + a.routine.duration);
    if (extra <= 1e-6) return;
    const steps = a.steps.slice();
    steps[steps.length - 1] = { ...last, duration: last.duration + extra };
    if (last.kind === "work") a.workDur += extra;
    a.steps = steps;
    a.routine = new Routine(a.st, a.start, steps, false);
  }

  /** Keep working `s` seconds longer (by I3 always safe: everyone else's legs keep clear of this spot), plan later. */
  private extend(a: Agent, s: number): void {
    const k = this.crowd.end(a.id) + Math.round(s / CROWD_DT);
    this.crowd.wait(a.id, k);
    a.endT = k * CROWD_DT;
    this.stretch(a, a.endT);
  }

  /** One slice of the planning job: a step of its preparation (finding the run, its paths), or of the check. */
  private slice(j: NonNullable<Crew["job"]>): void {
    if (j.prep) {
      const r = j.prep.next();
      if (!r.done) return;
      j.prep = null;
      if (!r.value) {
        this.job = null;
        if (j.s) {
          // No way to the station now: keep working a little longer and try again.
          this.log.waits++;
          this.extend(j.a, 1.5);
        } else this.plan(j.a, false);
        return;
      }
      // The preparation took some frames: start no sooner than a moment from now.
      j.k0 = Math.max(j.k0, this.crowd.step(this.now + 0.4));
      j.seq = r.value.seq;
      j.cands = r.value.cands;
      j.it = this.crowd.planLegSteps(j.a.id, j.k0, j.cands, 6);
      return;
    }
    const r = j.it!.next();
    if (r.done) this.finish(r.value);
  }

  /** Advance everyone to time t: plan the legs that are due, start them when their time comes, pose the bodies. */
  update(t: number, dt: number): void {
    this.now = t;
    const t0ms = performance.now();
    // One plan at a time, the most pressing first, given a few milliseconds a frame.
    if (!this.job) {
      let due: Agent | null = null;
      for (const a of this.agents) if (!a.planned && t >= a.endT - LOOKAHEAD && (!due || a.endT < due.endT)) due = a;
      if (due) this.plan(due);
    }
    // Slices until the budget is spent, not starting one that the last one's cost says would overrun it.
    let last = 0;
    while (this.job && (last === 0 || performance.now() - t0ms + last < PLAN_BUDGET_MS)) {
      const ts = performance.now();
      this.slice(this.job);
      last = Math.max(0.05, performance.now() - ts);
    }
    const ms = performance.now() - t0ms;
    this.log.planMs += ms;
    this.log.planMax = Math.max(this.log.planMax, ms);
    for (const a of this.agents) {
      if (a.planned && t >= a.planned.t0) {
        const prev = a.station;
        this.putDown(a);
        if (a.run && !a.run.settled) this.settle(a, a.run);
        if (prev && prev.holder === a.id) prev.holder = -1;
        const p = a.planned;
        a.planned = null;
        this.begin(a, p.start, p.steps, p.t0, p.station, p.run, p.workDur, p.path);
        a.endT = this.crowd.end(a.id) * CROWD_DT;
        this.drawPath(a);
      }
      const local = t - a.t0;
      const state = a.routine.apply(local);
      if (a.run) this.runFrame(a, local);
      else if (a.station) {
        const u = state.activity === "work" ? local - (a.routine.duration - a.workDur) : -1;
        if (a.station.kind === "fry") this.fry(a, u);
        else if (a.station.kind === "flip") this.flip(a, u);
        else if (a.station.kind === "sala") this.sala(a, u);
        else this.tool(a, u, dt);
      }
    }
    this.stripFrame(dt);
    this.tickPans();
    if (t - this.lastPrune > 5) {
      this.crowd.prune(this.crowd.step(t));
      this.lastPrune = t;
    }
  }

  /* -------------------------------------------------------------------------------------------------------- */
  /* Drawing a run: the door, the pan, the rack, the hands on them                                             */

  private runFrame(a: Agent, local: number): void {
    const r = a.run!;
    const st = a.st;
    let p: BodyPose | null = null;
    const pose = () => (p ??= copyPose(this.scratch, st.current));
    if (r.door) {
      const dp = r.door;
      // The door is this run's only from the hand on its pull to its latch (the crowd keeps two runs' holds of one
      // door apart); before and after, it is left as whoever else works it sets it.
      if (local >= dp.t0 && !dp.latched) {
        const f = doorOpening(dp, local);
        setDoor(dp.d, f);
        if (r.walkin) r.walkin.open = f;
        dp.driving = true;
        if (local >= dp.c1) {
          dp.latched = true;
          dp.driving = false;
        }
      }
      // The hand to the pull, on it while it pulls, then off as the door swings on.
      const w = smoothstep(dp.t0 - 0.6, dp.t0 - 0.05, local) * (1 - smoothstep(dp.t1, dp.t1 + 0.5, local));
      if (w > 1e-3) {
        const q = pose();
        const fh = doorOpening(dp, Math.min(local, dp.t1));
        const h = handleAt(dp.d, fh, _rv1);
        const palm = leafNormal(dp.d, fh, _rv2).negate();
        const close = smoothstep(dp.t0 - 0.3, dp.t0 - 0.02, local) * (1 - smoothstep(dp.t1, dp.t1 + 0.18, local));
        gripArm(st, q, dp.arm, w, h, UP, palm, close, 0.0095);
        q.gaze.lerp(st.toBody(h, _rv3), 0.6 * w);
      }
      // Shut by hand: the hand to the pull at the stop, on it as it brings the door round, off as it swings on.
      if (dp.rel !== undefined && dp.closeArm !== undefined) {
        const wc = smoothstep(dp.c0 - 0.5, dp.c0 - 0.05, local) * (1 - smoothstep(dp.rel, dp.rel + 0.45, local));
        if (wc > 1e-3) {
          const q = pose();
          const fh = doorOpening(dp, Math.min(local, dp.rel));
          // (By its free edge from the inner side: the palm on the edge, the fingers round onto the outer face.)
          const h = dp.edge ? edgeAt(dp.d, fh, _rv1) : handleAt(dp.d, fh, _rv1);
          const palm = dp.edge ? leafDir(dp.d, fh, _rv2).negate() : leafNormal(dp.d, fh, _rv2).negate();
          const close = smoothstep(dp.c0 - 0.3, dp.c0 - 0.02, local) * (1 - smoothstep(dp.rel, dp.rel + 0.15, local));
          gripArm(st, q, dp.closeArm, wc, h, UP, palm, close, dp.edge ? EDGE_R : 0.0095);
          q.gaze.lerp(st.toBody(h, _rv3), 0.6 * wc);
        }
        // Watching it swing to.
        if (local > dp.rel && local < dp.c1 + 0.3) {
          const q = pose();
          const g = handleAt(dp.d, doorOpening(dp, local), _rv1);
          q.gaze.lerp(st.toBody(g, _rv3), 0.5 * (1 - smoothstep(dp.c1, dp.c1 + 0.3, local)));
        }
      }
    }
    if (r.pan) {
      const P = r.pan;
      if (local >= P.t0 && local < P.t1 && !P.held && P.obj) {
        st.item.add(P.obj);
        P.obj.position.set(0, 0, 0);
        P.obj.quaternion.identity();
        P.held = true;
      }
      if (local >= P.t1 && !P.placed && P.obj) {
        this.group.add(P.obj);
        P.obj.position.copy(P.to.at);
        P.obj.rotation.set(0, panYaw(P.to), 0);
        P.placed = true;
      }
    }
    if (r.rack) {
      const R = r.rack;
      const rk = R.r;
      const pl = st.placement;
      const L = Math.hypot(R.to.x - R.from.x, R.to.z - R.from.z);
      const travel = local < R.t0 ? 0 : local >= R.t1 ? L : Math.max(0, Math.min(L, (pl.x - R.A[0]) * R.d[0] + (pl.z - R.A[1]) * R.d[1]));
      const c = _rv1.set(R.from.x + R.d[0] * travel, 0, R.from.z + R.d[1] * travel);
      rk.fx.model.group.position.copy(rk.base).add(c).sub(rk.origin);
      rk.fx.model.group.updateMatrixWorld(true);
      // The casters swing round to trail the way it goes over its first few centimetres.
      const sw = R.s0 + (R.s1 - R.s0) * smoothstep(0, 0.16, travel);
      rk.fx.model.hooks.set?.("swivel", sw);
      // The hands are on its corner posts by the push pose itself (people.ts: a push with `uprights`).
    }
    if (p) st.pose(p);
    if (!r.settled && local >= r.end) this.settle(a, r);
  }

  /** A run done: the pan where it was put, the rack parked where it was left, the bookings let go. */
  private settle(a: Agent, r: Run): void {
    void a;
    r.settled = true;
    // A door still swinging to (the run cut short) is shut; one this run never opened, or has latched, is left alone.
    if (r.door?.driving) {
      setDoor(r.door.d, 0);
      if (r.walkin) r.walkin.open = 0;
      r.door.driving = false;
      r.door.latched = true;
    }
    if (r.pan) {
      const P = r.pan;
      if (!P.placed && P.obj) {
        this.group.add(P.obj);
        P.obj.position.copy(P.to.at);
        P.obj.rotation.set(0, panYaw(P.to), 0);
        P.placed = true;
      }
      P.from.pan = null;
      P.to.pan = P.obj;
      P.from.usedAt = undefined;
      // A combi's pan cooks, and comes out; the next waits a while on the counter before it goes in.
      const C = P.to.combi;
      if (C) {
        if (P.to === C.spot) C.doneAt = this.now + 40 + 30 * this.rnd();
        else C.loadAt = this.now + 15 + 20 * this.rnd();
      } else if (P.to.dish) P.to.washAt = this.now + 10 + 10 * this.rnd();
      else if (!P.to.cold) P.to.usedAt = this.now + 35 + 35 * this.rnd();
      P.from.booked = -1;
      P.to.booked = -1;
    }
    if (r.rack) {
      const R = r.rack;
      R.r.out = !R.r.out;
      R.r.busy = false;
      R.r.swivel = R.s1;
      R.r.fx.model.group.position.copy(R.r.base).add(R.to).sub(R.r.origin);
      // Out of the way for good: it shuttles between its two docks from now on (or stays).
      if (R.r.next) {
        R.r.home = R.r.next;
        R.r.next = null;
      }
      R.r.fx.model.hooks.set?.("swivel", R.s1);
    }
  }

  /**
   * The strip curtains part round whoever goes through: the strips over the body (and over a pan carried ahead of
   * it) swing aside, away from its middle, far enough to clear the shoulders, and lean a little the way it is going;
   * they fall back slowly once it has passed. A strip is only touched from outside while the door is open.
   */
  /**
   * Pans on the line used up (their food gone), and those at the pot sink washed and, refilled, put back in a cold
   * store while no one can see in (a reach-in whose door no one is working, a walk-in with no one in it).
   */
  /** Whether a point is in the viewer's frame, set by the stage: no pan is restocked where it would be seen to appear. */
  inView: ((p: THREE.Vector3) => boolean) | null = null;

  private tickPans(): void {
    for (const sp of this.spots) {
      const pan = sp.pan;
      if (!pan) continue;
      if (!sp.cold && !sp.dish && sp.usedAt !== undefined && this.now >= sp.usedAt && !pan.userData.emptied) {
        (pan.userData.setFill as ((f: number) => void) | undefined)?.(0.05);
        pan.userData.emptied = true;
      }
      if (!sp.dish || sp.booked >= 0 || sp.washAt === undefined || this.now < sp.washAt) continue;
      const runs = this.agents.flatMap((b) => [b.run && !b.run.settled ? b.run : null, b.planned?.run ?? null]).filter((r): r is Run => !!r);
      const unseen = (q: Spot) => {
        if (q.door) return !runs.some((r) => r.door?.d === q.door);
        // A walk-in (its roof left off) and open shelving show what is on them: nor may the viewer see the pan appear.
        if (this.inView?.(q.at)) return false;
        if (q.walkin) {
          const R = q.walkin.fx.rect!;
          return !runs.some((r) => r.walkin === q.walkin) && !this.agents.some((b) => rectDistance(R, b.st.placement.x, b.st.placement.z) < 0.5);
        }
        // Open shelving: restocked by the prep crew off the floor plan, while no one is about.
        return !this.agents.some((b) => Math.hypot(q.face[0] - b.st.placement.x, q.face[1] - b.st.placement.z) < 4);
      };
      const to = this.spots.find((q) => q.cold && !q.pan && q.booked < 0 && q.surfaceY > 0.45 && q.surfaceY < 1.45 && unseen(q));
      if (!to) {
        sp.washAt = this.now + 4;
        continue;
      }
      pan.position.copy(to.at);
      pan.rotation.set(0, to.yaw, 0);
      (pan.userData.setFill as ((f: number) => void) | undefined)?.(0.7);
      pan.userData.emptied = false;
      to.pan = pan;
      sp.pan = null;
      sp.washAt = undefined;
    }
  }

  private stripFrame(dt: number): void {
    const v = _rv1;
    for (const w of this.walkins) {
      if (!w.strips.length) continue;
      const kIn = 1 - Math.exp(-dt / 0.09);
      const kOut = 1 - Math.exp(-dt / 0.5);
      const tx = this.stripT;
      const tz = this.stripZ;
      tx.length = tz.length = w.strips.length;
      tx.fill(0);
      tz.fill(0);
      for (const a of this.agents) {
        const pl = a.st.placement;
        v.set(pl.x, 0, pl.z).applyMatrix4(w.inv);
        const dz = v.z - w.rz;
        if (Math.abs(dz) > 0.75) continue;
        if (dz > 0.05 && w.open < 0.25) continue;
        // Which way the body faces through the doorway (fixture frame): in (-z) or out (+z).
        _rv2.set(Math.sin(pl.yaw), 0, Math.cos(pl.yaw)).transformDirection(w.inv);
        const dir = _rv2.z < 0 ? 1 : -1;
        const pan = !!(a.run?.pan?.held && !a.run.pan.placed);
        // The body (and the pan ahead of it, the way it faces) as a span along z about the root.
        const ahead = pan ? 0.55 : 0.16;
        const lead = dir > 0 ? -dz : dz;
        const deep = smoothstep(-0.35, -0.05, lead + ahead) * (1 - smoothstep(0.16, 0.42, lead));
        if (deep <= 0) continue;
        const half = pan && lead + 0.12 > -0.5 ? 0.33 : 0.27;
        w.strips.forEach((s, i) => {
          const lat = s.x - v.x;
          const need = half + 0.1 - Math.abs(lat);
          if (need <= 0) return;
          const side = Math.abs(lat) > 0.02 ? Math.sign(lat) : i % 2 ? 1 : -1;
          // Clear the shoulders 0.7 m below the rail.
          const b = Math.atan(need / 0.7) * deep;
          if (b > Math.abs(tz[i]!)) tz[i] = side * b;
          const f = dir * 0.2 * deep * Math.min(1, need / 0.2);
          if (Math.abs(f) > Math.abs(tx[i]!)) tx[i] = f;
        });
      }
      w.strips.forEach((s, i) => {
        // The model's tick sets each strip every frame; if it has not (paused), take back what was drawn last frame.
        if (s.j.quaternion.equals(s.mine)) s.j.quaternion.copy(s.base);
        s.base.copy(s.j.quaternion);
        const x = tx[i]!;
        const z = tz[i]!;
        s.ax += (x - s.ax) * (Math.abs(x) > Math.abs(s.ax) ? kIn : kOut);
        s.az += (z - s.az) * (Math.abs(z) > Math.abs(s.az) ? kIn : kOut);
        if (Math.abs(s.ax) + Math.abs(s.az) > 1e-4) s.j.quaternion.multiply(_rq.setFromEuler(_re.set(s.ax, 0, s.az)));
        s.mine.copy(s.j.quaternion);
      });
    }
  }
  private stripT: number[] = [];
  private stripZ: number[] = [];

  private putDown(a: Agent): void {
    const s = a.station;
    if (!s || !s.tool || !s.rest || s.tool.parent === this.group) return;
    this.group.add(s.tool);
    restMatrix(s.rest, s.tool.matrix).decompose(s.tool.position, s.tool.quaternion, s.tool.scale);
    a.holding = false;
    a.toolK = 0;
  }

  /**
   * A station's tool: the hand goes to it where it rests, takes it, works it, brings it back to rest and lets go.
   * Between rest and work the hand's target eases from one to the other, the tool in the hand.
   */
  private tool(a: Agent, u: number, dt: number): void {
    const s = a.station!;
    const tool = s.tool;
    const rest = s.rest;
    if (!tool || !rest || u < 0) return;
    const D = a.workDur;
    const want = u >= TOOL_IN && u < D - 1.0 - TOOL_BLEND ? 1 : 0;
    if (want && !a.holding && u >= TOOL_IN) {
      a.st.hand("R").add(tool);
      tool.position.set(0, 0, 0);
      tool.quaternion.identity();
      a.holding = true;
    }
    a.toolK = want ? Math.min(1, a.toolK + dt / TOOL_BLEND) : Math.max(0, a.toolK - dt / TOOL_BLEND);
    if (!want && a.holding && a.toolK <= 0) this.putDown(a);
    const k = smoothstep(0, 1, a.toolK);
    // A griddle's whole top and front trough are hot: the free hand stays off it.
    const hot = /GRIDDLE/.test(s.fx.id);
    if (k >= 1 && !hot) return;
    // Hold the arm towards the tool's rest, blended with the work's own target.
    const st = a.st;
    const yaw = st.placement.yaw;
    const p = copyPose(this.scratch, st.current);
    if (hot) p.arms[0].ik = 0;
    const R = p.arms[1];
    const rp = st.toBody(rest.p);
    const rz = dirToBody(rest.z, yaw);
    const ry = dirToBody(rest.y, yaw);
    R.pos.lerpVectors(rp, R.pos, k);
    R.a.lerpVectors(rz, R.a, k).normalize();
    R.b.lerpVectors(ry, R.b, k).normalize();
    if (s.turnOver) R.roll = Math.max(R.roll, Math.PI * 0.6 * (1 - k));
    if (a.holding || u < TOOL_IN) {
      R.shape = "tool";
      R.close = a.holding ? 1 : smoothstep(0.5, TOOL_IN, u);
    }
    p.gaze.lerpVectors(rp, p.gaze, k);
    st.pose(p);
  }

  /**
   * Fryer baskets. The right hand takes the basket's handle and the basket follows the lift the fryer's model
   * defines (out of the oil and onto the hanger rod; or off the rod and down), the hand on the handle throughout.
   */
  private fry(a: Agent, u: number): void {
    const s = a.station!;
    const job = s.fryJob;
    if (!job || !s.baskets || u < 0) return;
    const joint = s.fx.model.joints.get(job.k === 0 ? "basket_left_lift" : "basket_right_lift");
    if (!joint) return;
    const name = job.k === 0 ? "lift_l" : "lift_r";
    const b = s.baskets[job.k]!;
    let lift: number;
    let w: number;
    if (job.lift) {
      lift = 0.65 * smoothstep(1.3, 3.1, u) + 0.35 * smoothstep(3.9, 5.0, u);
      // A shake over the oil before it is hung.
      if (u > 3.1 && u < 3.9) lift += 0.02 * Math.sin((u - 3.1) * 24) * Math.sin(((u - 3.1) / 0.8) * Math.PI);
      w = smoothstep(0.25, 1.2, u) * (1 - smoothstep(5.1, 5.9, u));
      if (u > 5.0 && b.down) {
        b.down = false;
        b.since = this.now;
      }
    } else {
      lift = 1 - 0.35 * smoothstep(1.3, 2.0, u) - 0.65 * smoothstep(2.1, 3.5, u);
      w = smoothstep(0.25, 1.2, u) * (1 - smoothstep(3.6, 4.4, u));
      if (u > 3.5 && !b.down) {
        b.down = true;
        b.since = this.now;
      }
    }
    s.fx.model.hooks.set?.(name, lift);
    if (w <= 0.001) return;
    joint.updateMatrixWorld(true);
    const st = a.st;
    // The basket handle's grip in the joint's frame (cooking.ts fryBasket), its axis running into the fryer.
    const g = st.toBody(joint.localToWorld(new THREE.Vector3(0, 0.15, 0.1)));
    const q = new THREE.Quaternion();
    joint.getWorldQuaternion(q);
    const yawQ = new THREE.Quaternion().setFromAxisAngle(UP, -st.placement.yaw);
    const axis = new THREE.Vector3(0, 0, -1).applyQuaternion(q).applyQuaternion(yawQ);
    const up = new THREE.Vector3(0, 1, 0).applyQuaternion(q).applyQuaternion(yawQ);
    const p = copyPose(this.scratch, st.current);
    const R = p.arms[1];
    R.ik = w;
    R.grip = "tool";
    R.pos.copy(g);
    R.a.copy(axis);
    R.b.copy(up);
    R.pole.set(-0.6, -0.7, -0.35).normalize();
    R.shape = "tool";
    R.close = w;
    R.radius = 0.0135;
    R.roll = (70 * Math.PI) / 180;
    R.spinMax = 0;
    p.gaze.copy(g);
    st.pose(p);
  }

  /**
   * Steaks on a charbroiler, turned with tongs: the right hand takes the tongs off the front of the grates, and for
   * each steak in turn pinches its near edge, lifts it so it hangs from the tongs, carries the edge forward and down
   * so the steak lies over onto its other face where it was, lets go; then puts the tongs back.
   */
  private flip(a: Agent, u: number): void {
    const s = a.station!;
    const T = s.tongs;
    const steaks = s.steaks;
    if (!T || !steaks || u < 0) return;
    const st = a.st;
    const D = a.workDur;
    const hand = st.hand("R");
    const TAKE = 0.95;
    const PER = 2.6;
    const T0 = 1.25;
    const BACK = T0 + steaks.length * PER + 0.45;
    const w = smoothstep(0.05, 0.85, u) * (1 - smoothstep(D - 0.65, D - 0.05, u));
    const holding = u >= TAKE && u < BACK;
    if (holding && T.obj.parent !== hand) {
      hand.add(T.obj);
      T.obj.position.set(0, 0, 0);
      T.obj.quaternion.identity();
    }
    if (!holding && T.obj.parent !== this.group) {
      this.group.add(T.obj);
      restMatrix(T.rest, T.obj.matrix).decompose(T.obj.position, T.obj.quaternion, T.obj.scale);
      T.model.hooks.open?.(0);
    }
    if (w <= 1e-3) return;
    const yaw = st.placement.yaw;
    const L = bodyDir(yaw, 1, 0, 0, _fl1);
    const h = bodyDir(yaw, 0, 0, 1, _fl2);
    const grate = s.at.y;
    const GRIP = 0.185;
    const HALF = 0.042;
    // The cups' point E (world), the tongs' pitch down from level, how far open; the steak in them and its turn.
    const E = _fl3;
    let pitch = 0.5;
    let open = 0;
    let held = -1;
    let psi = 0;
    const edge = (i: number, off: number, lift: number, out = new THREE.Vector3()) => out.copy(steaks[i]!.at).addScaledVector(h, off).setY(grate + 0.016 + lift);
    const restE = (out: THREE.Vector3) => out.copy(T.rest.p).addScaledVector(T.rest.z, GRIP);
    if (u < T0) {
      // To the tongs at rest, take them, lift them clear.
      restE(E);
      E.y += 0.06 * smoothstep(TAKE, T0, u);
      pitch = 0.42 * smoothstep(TAKE, T0, u);
    } else if (u >= BACK - 0.45) {
      const v = smoothstep(BACK - 0.45, BACK, u);
      edge(steaks.length - 1, -HALF - 0.04, 0.06, _fl4);
      restE(E).lerp(_fl4, 1 - v);
      pitch = 0.42 * (1 - v);
    } else {
      const i = Math.min(steaks.length - 1, Math.floor((u - T0) / PER));
      const tau = u - T0 - i * PER;
      const from = i === 0 ? restE(_fl4).setY(_fl4.y + 0.06) : edge(i - 1, -HALF - 0.04, 0.06, _fl4);
      const app = edge(i, -HALF - 0.045, 0.06, _fl5);
      const grab = edge(i, -HALF, 0, _fl6);
      if (tau < 0.65) {
        E.lerpVectors(from, app, smoothstep(0, 0.65, tau));
        open = 0.5 * smoothstep(0.1, 0.5, tau);
        pitch = 0.42 + 0.13 * smoothstep(0.2, 0.65, tau);
      } else if (tau < 0.95) {
        E.lerpVectors(app, grab, smoothstep(0.65, 0.95, tau));
        open = 0.5;
        pitch = 0.55;
      } else if (tau < 1.15) {
        E.copy(grab);
        open = 0.5 - 0.27 * smoothstep(0.95, 1.12, tau);
        pitch = 0.55;
        if (tau >= 1.08) held = i;
      } else if (tau < 1.6) {
        // Lifted by the edge until it hangs.
        const v = smoothstep(1.15, 1.6, tau);
        psi = (Math.PI / 2) * v;
        E.copy(grab).setY(grate + 0.016 + 0.114 * Math.sin(psi) + 0.012 * v);
        open = 0.23;
        pitch = 0.55 - 0.2 * v;
        held = i;
      } else if (tau < 2.1) {
        // Over: the edge forward and down as the steak falls onto its other face.
        const v = smoothstep(1.6, 2.1, tau);
        psi = Math.PI / 2 + (Math.PI / 2) * v;
        E.copy(steaks[i]!.at).addScaledVector(h, -HALF + 2 * HALF * v).setY(grate + 0.016 + 0.114 * Math.sin(psi) + 0.012 * (1 - v));
        open = 0.23;
        pitch = 0.35 + 0.15 * v;
        held = i;
      } else {
        // Let go and away.
        const v = smoothstep(2.1, 2.3, tau);
        const back = smoothstep(2.25, PER, tau);
        E.copy(steaks[i]!.at).addScaledVector(h, HALF - (2 * HALF + 0.04) * back).setY(grate + 0.016 + 0.06 * back);
        open = 0.23 + 0.27 * v;
        pitch = 0.5 - 0.08 * back;
        psi = Math.PI;
        if (tau < 2.15) held = i;
        else this.flipped(steaks[i]!, L);
      }
    }
    // The steak in the tongs: turned about the hand's left-right axis through its pinched edge.
    if (held >= 0) {
      const sk = steaks[held]!;
      if (!sk.obj.userData.q0) sk.obj.userData.q0 = sk.obj.quaternion.clone();
      const q0 = sk.obj.userData.q0 as THREE.Quaternion;
      _rq.setFromAxisAngle(L, psi);
      sk.obj.quaternion.copy(_rq).multiply(q0);
      sk.obj.position.copy(E).addScaledVector(_fl4.copy(h).applyAxisAngle(L, psi), HALF);
      sk.up = 1;
    }
    T.model.hooks.open?.(open);
    // The hand: the tongs' axis from the grip to the cups, pitched down; their top up.
    const ax = _fl5.copy(h).multiplyScalar(Math.cos(pitch)).addScaledVector(UP, -Math.sin(pitch));
    const up = _fl6.copy(UP).multiplyScalar(Math.cos(pitch)).addScaledVector(h, Math.sin(pitch));
    if (u < TAKE || u >= BACK) {
      ax.copy(T.rest.z);
      up.copy(T.rest.y);
    }
    const g = _fl4.copy(E).addScaledVector(ax, -GRIP);
    if (u < TAKE || u >= BACK) g.copy(T.rest.p);
    const p = copyPose(this.scratch, st.current);
    const R = p.arms[1];
    R.ik = w;
    R.grip = "tool";
    st.toBody(g, R.pos);
    dirToBody(ax, yaw, R.a);
    dirToBody(up, yaw, R.b);
    R.pole.set(-0.6, -0.75, -0.3).normalize();
    R.shape = "tool";
    R.close = holding ? 1 : smoothstep(0.55, TAKE, u) * (1 - smoothstep(BACK, BACK + 0.3, u));
    R.radius = 0.0095;
    R.roll = 0;
    R.spinMax = 0;
    st.toBody(E, _fl3);
    p.gaze.lerp(_fl3, w);
    st.pose(p);
  }

  /** A steak let go after its turn: lying on its other face where it was, seared on both once it has been turned. */
  private flipped(sk: NonNullable<Station["steaks"]>[number], L: THREE.Vector3): void {
    if (!sk.up) return;
    const q0 = sk.obj.userData.q0 as THREE.Quaternion;
    sk.obj.quaternion.copy(_rq.setFromAxisAngle(L, Math.PI)).multiply(q0);
    sk.obj.position.copy(sk.at);
    sk.obj.userData.q0 = null;
    sk.up = 0;
    sk.done++;
    const seared = sk.obj.userData.seared as THREE.Material | undefined;
    if (seared) for (const c of sk.obj.children) (c as THREE.Mesh).material = seared;
  }

  /**
   * A salamander: the right hand up to the grid's pull (leaning in over the range below as far as the reach needs),
   * the grid drawn out to look at what is under the broiler, then slid back in.
   */
  private sala(a: Agent, u: number): void {
    const s = a.station!;
    if (u < 0) return;
    const st = a.st;
    const out = smoothstep(1.1, 2.6, u) * (1 - smoothstep(5.6, 7.0, u));
    s.fx.model.hooks.open?.(out);
    const w = smoothstep(0.1, 1.05, u) * (1 - smoothstep(7.1, 8.1, u));
    if (w <= 1e-3) return;
    const n = _fl1.set(s.nx, 0, s.nz);
    const pull = _fl2.copy(s.at).addScaledVector(n, 0.175 * out);
    const yaw = st.placement.yaw;
    const left = bodyDir(yaw, 1, 0, 0, _fl3);
    const palm = _fl4.copy(n).multiplyScalar(-0.6).addScaledVector(UP, -0.8).normalize();
    const p = copyPose(this.scratch, st.current);
    gripArm(st, p, 1, w, pull, left, palm, smoothstep(0.65, 1.05, u) * (1 - smoothstep(7.0, 7.4, u)), 0.006);
    const lean = a.lean * w;
    p.tilt += 0.4 * lean;
    p.flex += 0.6 * lean;
    const look = _fl5.copy(pull).addScaledVector(n, -0.2).setY(pull.y - 0.02);
    p.gaze.lerp(st.toBody(look, _fl6), w);
    st.pose(p);
  }

  /** How far a person standing at `at` must lean in for the right hand to reach a bar there (a salamander's pull). */
  private leanFor(a: Agent, at: { x: number; z: number; yaw: number }, s: Station): number {
    const g = a.ghost;
    g.place(at.x, at.z, at.yaw);
    const n = new THREE.Vector3(s.nx, 0, s.nz);
    const left = bodyDir(at.yaw, 1, 0, 0);
    const palm = n.clone().multiplyScalar(-0.6).addScaledVector(UP, -0.8).normalize();
    const p = newPose();
    for (let lean = 0; lean <= 0.6; lean += 0.05) {
      standPose(g, 0, p);
      p.tilt += 0.4 * lean;
      p.flex += 0.6 * lean;
      gripArm(g, p, 1, 1, s.at, left, palm, 1, 0.006);
      g.pose(p);
      if (g.gripError("R") < 0.008) return lean;
    }
    return 0.6;
  }

  private drawPath(a: Agent): void {
    if (a.line) {
      this.lines.remove(a.line);
      a.line.geometry.dispose();
      (a.line.material as THREE.Material).dispose();
      a.line = null;
    }
    if (!a.path || a.path.length < 2) return;
    const pts = a.path.map(([x, z]) => new THREE.Vector3(x, 0.02, z));
    const line = new THREE.Line(new THREE.BufferGeometry().setFromPoints(pts), new THREE.LineDashedMaterial({ color: a.accent, dashSize: 0.12, gapSize: 0.08 }));
    line.computeLineDistances();
    line.renderOrder = 6;
    a.line = line;
    this.lines.add(line);
  }

  /** People in the thermal view: clothed skin near 33 °C. */
  thermal(on: boolean, mat: THREE.Material): void {
    for (const a of this.agents)
      for (const m of a.st.meshes) {
        if (!m.userData.pbr) m.userData.pbr = m.material;
        m.material = on ? mat : (m.userData.pbr as THREE.Material);
      }
  }

  tris(): number {
    return this.agents.reduce((s, a) => s + a.st.tris, 0);
  }

  /** The runs under way (scripted checks): who, where, how far into the leg, and its marks (seconds into the leg). */
  runInfo(kind?: RunKind) {
    return this.agents
      .filter((a) => a.run && !a.run.settled && (!kind || a.run.kind === kind))
      .map((a) => {
        const r = a.run!;
        const p = a.st.placement;
        return {
          id: a.id,
          kind: r.kind,
          local: this.now - a.t0,
          x: p.x,
          z: p.z,
          yaw: p.yaw,
          times: r.times,
          end: r.end,
          door: r.door ? [r.door.t0, r.door.t1, r.door.t2, r.door.c0, r.door.c1] : null,
          pan: r.pan ? [r.pan.t0, r.pan.t1] : null,
          rack: r.rack ? [r.rack.t0, r.rack.t1] : null,
          gate: !!r.walkin,
          grip: [a.st.gripError("L"), a.st.gripError("R")],
        };
      });
  }

  dispose(): void {
    // What the crew drew besides the people (pans, tools, plates, steaks, tongs, paths): its geometry, and the
    // materials made for it alone (the shared palette's stay).
    const geos = new Set<THREE.BufferGeometry>();
    const take = (o: THREE.Object3D) =>
      o.traverse((q) => {
        const m = q as THREE.Mesh;
        if (m.geometry) geos.add(m.geometry as THREE.BufferGeometry);
      });
    const bodies = new Set<THREE.Object3D>(this.agents.map((a) => a.st.group));
    for (const c of this.group.children) if (!bodies.has(c)) take(c);
    for (const a of this.agents) for (const c of a.st.item.children) take(c);
    for (const s of this.stations) {
      if (s.tool) take(s.tool);
      if (s.tongs) take(s.tongs.obj);
    }
    take(this.lines);
    for (const g of geos) g.dispose();
    for (const l of this.lines.children) ((l as THREE.Line).material as THREE.Material).dispose();
    const shared = new Set<THREE.Material>(Object.values(palette()) as THREE.Material[]);
    for (const m of this.ownMats) if (!shared.has(m)) m.dispose();
    for (const a of this.agents) a.st.dispose();
    for (const g of this.ghosts.values()) g.dispose();
  }
  /** Materials made for the crew's own props (plates, bains, a brush, steaks), disposed with it. */
  private readonly ownMats = new Set<THREE.Material>();
  private own(o: THREE.Object3D): void {
    o.traverse((q) => {
      const m = (q as THREE.Mesh).material;
      if (m) for (const x of Array.isArray(m) ? m : [m]) this.ownMats.add(x);
    });
  }
}


/* ------------------------------------------------------------------------------------------------------------ */
/* Camera presets, framed on this kitchen's own fixtures and stations                                            */
/* ------------------------------------------------------------------------------------------------------------ */

type Pose = { pos: [number, number, number]; target: [number, number, number] };

function framePresets(K: Kitchen, crew: Crew, spec: ExemplarSpec): Record<CameraPresetId, Pose> {
  const out: Record<CameraPresetId, Pose> = { ...spec.cameras } as unknown as Record<CameraPresetId, Pose>;
  const inside = (x: number, y: number, z: number): [number, number, number] => [
    Math.min(K.W - 0.25, Math.max(0.25, x)),
    Math.min(K.H - 0.3, y),
    Math.min(-0.25, Math.max(-K.D + 0.25, z)),
  ];
  /** A view from in front of a spot (normal n), off to one side, looking back at it. */
  const view = (cx: number, cz: number, nx: number, nz: number, back: number, side: number, y: number, ty: number, cy = 0): Pose => {
    const sx = -nz;
    const sz = nx;
    return { pos: inside(cx + nx * back + sx * side, y, cz + nz * back + sz * side), target: [cx + nx * cy, ty, cz + nz * cy] };
  };
  const hot = K.fixtures.filter((f) => f.rect && /RANGE|GRIDDLE|FRYER|CHARBROILER|INDUCTION|TILT|KETTLE/.test(f.id));
  if (hot.length) {
    const mid = (f: Fx) => [(f.rect!.minX + f.rect!.maxX) / 2, (f.rect!.minZ + f.rect!.maxZ) / 2] as const;
    // The line: the hot fixture with the most hot neighbours, and those neighbours.
    let best = hot[0]!;
    let bestN = -1;
    for (const f of hot) {
      const [x, z] = mid(f);
      const n = hot.filter((g) => Math.hypot(mid(g)[0] - x, mid(g)[1] - z) < 2.6).length;
      if (n > bestN) [best, bestN] = [f, n];
    }
    const [bx, bz] = mid(best);
    const line = hot.filter((g) => Math.hypot(mid(g)[0] - bx, mid(g)[1] - bz) < 2.6);
    const cx = line.reduce((q, f) => q + mid(f)[0], 0) / line.length;
    const cz = line.reduce((q, f) => q + mid(f)[1], 0) / line.length;
    const [nx, nz] = frontNormal(best);
    const span = Math.max(...line.map((f) => Math.abs((mid(f)[0] - cx) * -nz + (mid(f)[1] - cz) * nx))) + 0.6;
    out.cookline = view(cx, cz, nx, nz, 2.2, Math.min(1.6, span * 0.8), 2.3, 0.95, 0.15);
    out.hood = view(cx, cz, nx, nz, 3.2, -0.6, 2.05, 1.55, 0.1);
  }
  const pass = crew.stations.find((st) => st.kind === "plate") ?? crew.stations.find((st) => st.kind === "chop");
  if (pass) out.pass = view(pass.at.x, pass.at.z, pass.nx, pass.nz, 1.9, 0.9, 2.25, pass.at.y, 0.1);
  const cold = K.fixtures.find((f) => /WALKIN/.test(f.id)) ?? K.fixtures.find((f) => /REACHIN/.test(f.id));
  if (cold) {
    const st = fxSocket(cold, "stance_front");
    const [w] = cold.f.local_extents_m;
    const at = st ? fxLocal(cold, st.local_pos_m[0], st.local_pos_m[1], 1.1) : fxLocal(cold, w / 2, -0.3, 1.1);
    const [nx, nz] = frontNormal(cold);
    out.walkin = view(at.x, at.z, nx, nz, 2.7, 1.1, 2.3, 1.1);
  }
  return out;
}

/* ------------------------------------------------------------------------------------------------------------ */
/* The stage                                                                                                      */
/* ------------------------------------------------------------------------------------------------------------ */

export function mountHeroKitchen(root: HTMLElement): HeroController {
  const stageEl = root.querySelector<HTMLElement>("[data-hero-stage]")!;
  const canvas = root.querySelector<HTMLCanvasElement>("[data-hero-canvas]")!;
  const statusEl = root.querySelector<HTMLElement>("[data-hero-status]");

  const maybeRenderer = createStudioRenderer(canvas, true);
  if (!maybeRenderer) {
    // The poster stays as a still of the kitchen.
    if (statusEl) statusEl.textContent = "A still: the live kitchen needs WebGL";
    return { ready: Promise.resolve(), setExemplar: async () => {}, getLoadedExemplarId: () => "mall_qsr_400", getDrawCalls: () => 0 };
  }
  const renderer = maybeRenderer;
  renderer.toneMappingExposure = 1.0;
  const scene = new THREE.Scene();
  scene.environment = studioEnvironment(renderer);
  scene.environmentIntensity = 0.62;
  scene.fog = new THREE.Fog(0x0b0807, 18, 60);
  const camera = new THREE.PerspectiveCamera(38, 16 / 9, 0.05, 160);
  const initialSpec = EXEMPLARS[0]!;
  camera.position.set(...initialSpec.cameras.isometric.pos);
  const controls = createScrollSafeOrbit(camera, canvas, initialSpec.cameras.isometric.target, { minDistance: 0.6, maxDistance: 46, maxPolarAngle: Math.PI * 0.485 });
  const lights = studioLights(scene, new THREE.Vector3(2.5, 0, -4), 6, 2048);
  const gate = shaderGate(renderer, camera);
  lights.key.intensity = 2.3;
  lights.hemi.intensity = 0.42;
  lights.rim.intensity = 0.55;

  // The stage fills the hero behind its copy: on wide screens the copy takes the left and the console the bottom,
  // so the picture's centre moves right and up.
  let wide = false;
  // Resizing the canvas waits on the GPU, so behind a shader compile it is held until the compile is done.
  const fit = () => gate.defer("fit", fitNow);
  const fitNow = () => {
    const w = Math.max(1, Math.floor(stageEl.clientWidth));
    const h = Math.max(1, Math.floor(stageEl.clientHeight));
    renderer.setSize(w, h, false);
    camera.aspect = w / h;
    wide = w >= 1000;
    if (wide) camera.setViewOffset(w, h, -w * 0.2, h * 0.06, w, h);
    else camera.clearViewOffset();
    camera.updateProjectionMatrix();
  };
  fit();
  new ResizeObserver(fit).observe(stageEl);
  const phone = Math.min(window.innerWidth, window.innerHeight) < 600;

  let currentSpec: ExemplarSpec = initialSpec;
  let kitchen: Kitchen | null = null;
  let activeCameraPreset: CameraPresetId = "isometric";
  const activeOverlays = new Set<OverlayMode>(["pbr", "plumes"]);
  let paused = false;
  let simT = 0;
  let inView = true;
  let loadSeq = 0;
  let lastDraws = 0;
  /** "<kitchen> · <tris> tris" once a kitchen is built, null while one is loading. */
  let statusBase: string | null = null;
  const thermalMat = thermalMaterial();
  const plumes = new Plumes();
  scene.add(plumes.points);
  let aisles: THREE.Group | null = null;
  let crew: Crew | null = null;
  let presets: Record<CameraPresetId, Pose> | null = null;
  const peopleThermal = new THREE.MeshBasicMaterial({ color: thermalColor(33, T_MIN, T_MAX) });
  let metrics: HeroMetrics = computeKitchenTelemetry(initialSpec, initialSpec.nominalOrdersPerHr, 35, 26);

  let camAnimT = 1;
  const camStartPos = new THREE.Vector3();
  const camGoalPos = new THREE.Vector3(...initialSpec.cameras.isometric.pos);
  const tgtStartPos = new THREE.Vector3();
  const tgtGoalPos = new THREE.Vector3(...initialSpec.cameras.isometric.target);

  let ordersPerHr = initialSpec.nominalOrdersPerHr;
  let electrificationPct = 35;
  let ambientC = 26;

  const archetypeBtns = Array.from(root.querySelectorAll<HTMLButtonElement>("[data-archetype]"));
  const overlayBtns = Array.from(root.querySelectorAll<HTMLButtonElement>("[data-overlay]"));
  const camBtns = Array.from(root.querySelectorAll<HTMLButtonElement>("[data-cam]"));
  const pauseBtn = root.querySelector<HTMLButtonElement>("[data-hero-pause]");
  const rushInput = root.querySelector<HTMLInputElement>("[data-hero-rush]");
  const elecInput = root.querySelector<HTMLInputElement>("[data-hero-elec]");
  const ambInput = root.querySelector<HTMLInputElement>("[data-hero-amb]");
  const rushVal = root.querySelector<HTMLElement>("[data-hero-rush-val]");
  const elecVal = root.querySelector<HTMLElement>("[data-hero-elec-val]");
  const ambVal = root.querySelector<HTMLElement>("[data-hero-amb-val]");
  const outKw = root.querySelector<HTMLElement>("[data-out-kw]");
  const outCapture = root.querySelector<HTMLElement>("[data-out-capture]");
  const outAisle = root.querySelector<HTMLElement>("[data-out-aisle]");
  const outMeta = root.querySelector<HTMLElement>("[data-out-meta]");
  const barGas = root.querySelector<HTMLElement>("[data-bar='gas']");
  const barElec = root.querySelector<HTMLElement>("[data-bar='elec']");
  const barHvac = root.querySelector<HTMLElement>("[data-bar='hvac']");
  const keyGas = root.querySelector<HTMLElement>("[data-key='gas']");
  const keyElec = root.querySelector<HTMLElement>("[data-key='elec']");
  const keyHvac = root.querySelector<HTMLElement>("[data-key='hvac']");

  function lerpToPreset(preset: CameraPresetId): void {
    activeCameraPreset = preset;
    for (const b of camBtns) b.setAttribute("aria-pressed", String(b.dataset.cam === preset));
    const pose = presets ?? currentSpec.cameras;
    const pz = pose[preset];
    camStartPos.copy(camera.position);
    tgtStartPos.copy(controls.target);
    camGoalPos.set(pz.pos[0], pz.pos[1], pz.pos[2]);
    tgtGoalPos.set(pz.target[0], pz.target[1], pz.target[2]);
    camAnimT = 0;
  }

  function applyOverlays(): void {
    const K = kitchen;
    if (!K) return;
    const thermal = activeOverlays.has("thermal");
    const T = thermalMat.uniforms.uT!.value as Float32Array;
    T[0] = metrics.thermalComfortC;
    for (const fx of K.fixtures) if (fx.i < MAX_FIX) T[fx.i] = fx.tempC;
    for (const m of K.thermalMeshes) {
      if (!m.userData.pbr) m.userData.pbr = m.material;
      m.material = thermal ? thermalMat : (m.userData.pbr as THREE.Material);
    }
    plumes.points.visible = activeOverlays.has("plumes");
    if (aisles) aisles.visible = activeOverlays.has("aisles");
    if (crew) {
      crew.lines.visible = activeOverlays.has("aisles");
      crew.thermal(thermal, peopleThermal);
    }
    scene.environmentIntensity = thermal ? 0 : 0.62;
    for (const L of K.hoodLights) L.visible = !thermal;
  }

  function applyFixtures(): void {
    const K = kitchen;
    if (!K) return;
    applyHeat(K, metrics);
    plumes.capture = metrics.hoodCapturePct / 100;
    // Lit burners: the share of the range's valves the gas input calls for, the rest off.
    for (const fx of K.fixtures) {
      const fl = fx.model.hooks.flames;
      if (!fl || !/RANGE/.test(fx.id)) continue;
      const n = fx.model.hooks.burners ?? 6;
      const lit = Math.round(fx.firing * n);
      for (let b = 0; b < n; b++) fl.setFlow(b, b < lit ? 0.35 + 0.65 * fx.firing : 0);
    }
    for (const L of K.hoodLights) L.intensity = 2.2;
  }

  function syncTelemetryUI(): void {
    metrics = computeKitchenTelemetry(currentSpec, ordersPerHr, electrificationPct, ambientC);
    const m = metrics;
    if (rushVal) rushVal.textContent = `${Math.round(ordersPerHr)} ord/hr`;
    if (elecVal) elecVal.textContent = `${Math.round(electrificationPct)}%`;
    if (ambVal) ambVal.textContent = `${Math.round(ambientC)} °C`;
    if (outKw) outKw.textContent = `${m.thermalKw.toFixed(1)} kW`;
    if (outCapture) outCapture.textContent = `${m.hoodCapturePct.toFixed(1)}%`;
    if (outAisle) outAisle.textContent = `${m.minAisleIn.toFixed(0)}" · ${m.minAisleM.toFixed(2)}m`;
    if (outMeta) {
      const draws = kitchen ? `${lastDraws} draws` : "loading";
      outMeta.textContent = `${kitchen ? kitchen.fixtures.length : currentSpec.fixtures} fixtures · ${crew ? crew.agents.length : currentSpec.staffCount} staff · ${draws}`;
    }
    if (barGas) barGas.style.flexGrow = String(Math.max(1, m.gasSharePct));
    if (barElec) barElec.style.flexGrow = String(Math.max(1, m.elecSharePct));
    if (barHvac) barHvac.style.flexGrow = String(Math.max(1, m.hvacSharePct));
    if (keyGas) keyGas.textContent = `${(m.thermalKw - m.electricalKw).toFixed(1)} kW`;
    if (keyElec) keyElec.textContent = `${m.electricalKw.toFixed(1)} kW`;
    if (keyHvac) keyHvac.textContent = `${m.hvacCoolingKw.toFixed(1)} kW`;
    applyFixtures();
    applyOverlays();
  }

  async function loadExemplar(id: ExemplarId): Promise<void> {
    const seq = ++loadSeq;
    const spec = EXEMPLARS.find((e) => e.id === id) ?? EXEMPLARS[0]!;
    currentSpec = spec;
    for (const b of archetypeBtns) b.setAttribute("aria-pressed", String(b.dataset.archetype === id));
    if (rushInput) {
      // On the slider's step of 5, as the nominal and top rates are, so the slider can sit on the nominal rate.
      rushInput.min = String(Math.round((spec.nominalOrdersPerHr * 0.35) / 5) * 5);
      rushInput.max = String(spec.maxOrdersPerHr);
      ordersPerHr = spec.nominalOrdersPerHr;
      rushInput.value = String(ordersPerHr);
    }
    presets = null;
    statusBase = null;
    if (statusEl) statusEl.textContent = `Building ${spec.shortName}…`;
    const base = `${ASSET_BASE}/exemplars/${spec.id}/${spec.id}`;
    const [manifest, meta] = await Promise.all([fetchJson<Manifest>(`${base}.scene3d.v1.json`), fetchJson<SiteMeta>(`${base}.site_metadata.json`)]);
    if (seq !== loadSeq) return;
    if (!manifest) {
      if (statusEl) statusEl.textContent = `${spec.shortName} · manifest unavailable`;
      return;
    }
    const K = await buildKitchen(spec, manifest, meta, {
      phone,
      onProgress: (d, n) => {
        if (statusEl && seq === loadSeq) statusEl.textContent = `Building ${spec.shortName} · ${d}/${n} fixtures`;
      },
    }, () => seq === loadSeq);
    if (!K || seq !== loadSeq) {
      K?.dispose();
      return;
    }
    // The crew is built before its room replaces the one on screen: the stations and runs in a task apart from the
    // room's last fixtures and its merge, then a person to a task.
    await new Promise((r) => setTimeout(r, 0));
    if (seq !== loadSeq) {
      K.dispose();
      return;
    }
    const C = new Crew(K);
    if (!(await C.populate(!phone, phone, () => seq === loadSeq))) {
      C.dispose();
      K.dispose();
      return;
    }
    if (kitchen) {
      scene.remove(kitchen.root);
      kitchen.dispose();
    }
    if (aisles) {
      scene.remove(aisles);
      aisles = null;
    }
    if (crew) {
      scene.remove(crew.group, crew.lines);
      crew.dispose();
      crew = null;
    }
    kitchen = K;
    scene.add(K.root);
    crew = C;
    crew.inView = panInView;
    scene.add(crew.group, crew.lines);
    presets = framePresets(K, crew, spec);
    simT = 0;
    // A detailed body swapped in brings its own shaders: compile them before the next draw.
    for (const a of crew.agents)
      void a.st.ready.then(() => {
        applyOverlays();
        void gate.refresh(scene);
      });
    // Light the room: the key over its middle, its shadow covering it.
    const span = Math.max(K.W, K.D) * 0.62;
    const centre = new THREE.Vector3(K.W / 2, 0, -K.D / 2);
    lights.key.position.copy(centre).add(new THREE.Vector3(span * 0.55, span * 1.1, span * 0.75));
    lights.key.target.position.copy(centre);
    const sc = lights.key.shadow.camera;
    sc.left = -span;
    sc.right = span;
    sc.top = span;
    sc.bottom = -span;
    sc.far = span * 4;
    sc.updateProjectionMatrix();
    lights.key.shadow.normalBias = 0.012 * Math.max(0.2, span / 4);
    lights.rim.position.copy(centre).add(new THREE.Vector3(-span * 0.9, span * 0.6, -span * 0.8));
    scene.fog = new THREE.Fog(0x0b0807, span * 3, span * 7);
    plumes.setSources(K);
    aisles = buildAisles(K, 0.26);
    scene.add(aisles);
    lerpToPreset(activeCameraPreset);
    syncTelemetryUI();
    // The room's shaders compile in the background; until they are linked the frame loop holds its draw (the poster,
    // or the last kitchen, stays on screen).
    await gate.refresh(scene);
    if (seq !== loadSeq) return;
    renderer.render(scene, camera);
    lastDraws = renderer.info.render.calls;
    // The draw-call count is refreshed every 60 frames from then on (it changes with the view, and counts the shadow
    // pass), so the line shows the frame on screen rather than this first one.
    statusBase = `${spec.shortName} · ${Math.round(K.tris + (crew?.tris() ?? 0)).toLocaleString()} tris`;
    if (statusEl) statusEl.textContent = `${statusBase} · ${lastDraws} draw calls`;
    syncTelemetryUI();
    stageEl.dataset.live = "1";
  }

  for (const btn of archetypeBtns) btn.addEventListener("click", () => void loadExemplar(btn.dataset.archetype as ExemplarId));
  for (const btn of camBtns) btn.addEventListener("click", () => lerpToPreset(btn.dataset.cam as CameraPresetId));
  for (const btn of overlayBtns) {
    btn.addEventListener("click", () => {
      const mode = btn.dataset.overlay as OverlayMode;
      if (mode === "pbr") {
        activeOverlays.delete("thermal");
        activeOverlays.add("pbr");
      } else if (mode === "thermal") {
        if (activeOverlays.has("thermal")) {
          activeOverlays.delete("thermal");
          activeOverlays.add("pbr");
        } else {
          activeOverlays.delete("pbr");
          activeOverlays.add("thermal");
        }
      } else if (activeOverlays.has(mode)) activeOverlays.delete(mode);
      else activeOverlays.add(mode);
      for (const ob of overlayBtns) ob.setAttribute("aria-pressed", String(activeOverlays.has(ob.dataset.overlay as OverlayMode)));
      applyOverlays();
    });
  }
  if (pauseBtn) {
    pauseBtn.addEventListener("click", () => {
      paused = !paused;
      pauseBtn.setAttribute("aria-pressed", String(paused));
    });
  }
  rushInput?.addEventListener("input", () => {
    ordersPerHr = Number(rushInput.value);
    syncTelemetryUI();
  });
  elecInput?.addEventListener("input", () => {
    electrificationPct = Number(elecInput.value);
    syncTelemetryUI();
  });
  ambInput?.addEventListener("input", () => {
    ambientC = Number(ambInput.value);
    syncTelemetryUI();
  });
  window.addEventListener("keydown", (e) => {
    const tag = (e.target as HTMLElement | null)?.tagName;
    if (tag === "INPUT" || tag === "TEXTAREA" || tag === "SELECT") return;
    if (e.key >= "1" && e.key <= "5") {
      const spec = EXEMPLARS[Number(e.key) - 1];
      if (spec) void loadExemplar(spec.id);
    } else if (e.code === "Space" && inView) {
      e.preventDefault();
      paused = !paused;
      pauseBtn?.setAttribute("aria-pressed", String(paused));
    }
  });
  observeInView(stageEl, (visible) => {
    inView = visible;
  });

  let lastTs = performance.now();
  let frames = 0;
  let frameMs = 0;
  // Frame timing for scripted checks: mean ms per frame over the last second.
  const perf = { ms: 0, cpu: 0, cpuMax: 0 };
  let cpuSum = 0;
  let cpuMax = 0;
  (stageEl as unknown as { heroPerf: typeof perf }).heroPerf = perf;
  (stageEl as unknown as { heroCrew: () => Crew | null }).heroCrew = () => crew;
  // Run the kitchen ahead by s seconds in 50 ms steps; ask for a kind of run next; see the runs under way (checks).
  const hooksEl = stageEl as unknown as {
    heroJump: (s: number) => number;
    heroForce: (k: RunKind | null) => void;
    heroRun: (k?: RunKind) => ReturnType<Crew["runInfo"]>;
  };
  hooksEl.heroJump = (s) => {
    const K = kitchen;
    if (!K || !crew) return simT;
    for (let t = 0; t < s - 1e-9; t += 0.05) {
      simT += 0.05;
      for (const fx of K.fixtures) fx.model.hooks.tick?.(simT, 0.05, fx.firing);
      crew.update(simT, 0.05);
    }
    return simT;
  };
  hooksEl.heroForce = (k) => {
    if (crew) crew.force = k;
  };
  hooksEl.heroRun = (k) => (crew ? crew.runInfo(k) : []);
  // Point the camera somewhere, for scripted checks.
  (stageEl as unknown as { heroLook: (p: number[], t: number[]) => void }).heroLook = (p, t) => {
    camStartPos.copy(camera.position);
    tgtStartPos.copy(controls.target);
    camGoalPos.set(p[0]!, p[1]!, p[2]!);
    tgtGoalPos.set(t[0]!, t[1]!, t[2]!);
    camAnimT = 0;
  };

  /** Whether a pan at p would show in the frame (a pan's half-length round it). */
  const viewFrustum = new THREE.Frustum();
  const viewPV = new THREE.Matrix4();
  const viewBall = new THREE.Sphere(new THREE.Vector3(), 0.3);
  function panInView(p: THREE.Vector3): boolean {
    camera.updateMatrixWorld();
    viewPV.multiplyMatrices(camera.projectionMatrix, camera.matrixWorldInverse);
    viewFrustum.setFromProjectionMatrix(viewPV);
    viewBall.center.copy(p);
    return viewFrustum.intersectsSphere(viewBall);
  }

  /**
   * Clip planes that follow the camera: everything drawn lies inside the room's shell, so the near plane may stand
   * out to 0.6 of the distance to the shell's nearest point (a plane, it cuts nothing nearer than that within 53° of
   * the view axis; the frustum's widest corner, the wide stage's view offset included, is 40° off it on a 16:10
   * screen, under 50° on one twice as wide as the stage is tall), and the far
   * plane just past its farthest corner. From outside the room this keeps the depth buffer fine enough for the
   * millimetre steps between surfaces at the far side of the biggest kitchen.
   */
  const clipBox = new THREE.Box3();
  function clipToRoom(K: Kitchen): void {
    const m = 0.15;
    clipBox.min.set(-m, -m, -K.D - m);
    clipBox.max.set(K.W + m, K.H + m, m);
    fitClipPlanes(camera, clipBox, 0.05);
  }

  function animate(now: number): void {
    requestAnimationFrame(animate);
    const dt = Math.max(0, Math.min(0.05, (now - lastTs) / 1000));
    frameMs += now - lastTs;
    lastTs = now;
    if (!inView || document.hidden || !gate.ready()) return;
    const cpu0 = performance.now();
    if (camAnimT < 1) {
      camAnimT = Math.min(1, camAnimT + dt * 1.9);
      const ease = 1 - Math.pow(1 - camAnimT, 3);
      camera.position.lerpVectors(camStartPos, camGoalPos, ease);
      controls.target.lerpVectors(tgtStartPos, tgtGoalPos, ease);
    }
    controls.update();
    const K = kitchen;
    if (K) {
      cutWalls(K.walls, camera.position);
      clipToRoom(K);
      if (!paused) {
        simT += dt;
        for (const fx of K.fixtures) fx.model.hooks.tick?.(simT, dt, fx.firing);
        if (plumes.points.visible) plumes.tick(dt);
        crew?.update(simT, dt);
      }
    }
    renderer.render(scene, camera);
    lastDraws = renderer.info.render.calls;
    const cpu = performance.now() - cpu0;
    cpuSum += cpu;
    cpuMax = Math.max(cpuMax, cpu);
    if (++frames >= 60) {
      if (statusBase && statusEl) statusEl.textContent = `${statusBase} · ${lastDraws} draw calls`;
      perf.cpu = cpuSum / frames;
      perf.cpuMax = cpuMax;
      cpuSum = 0;
      cpuMax = 0;
      perf.ms = frameMs / frames;
      frames = 0;
      frameMs = 0;
    }
  }
  requestAnimationFrame(animate);
  // Settles either way, for what waits on it; a failure is still reported, as it was before anything waited.
  const ready = loadExemplar(initialSpec.id).catch((e: unknown) => console.error(e));
  void activeCameraPreset;
  void wide;

  return {
    ready,
    setExemplar: loadExemplar,
    getLoadedExemplarId: () => currentSpec.id,
    getDrawCalls: () => lastDraws,
  };
}
