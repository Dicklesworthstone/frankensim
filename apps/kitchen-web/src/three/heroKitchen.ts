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
import { ASSET_BASE, createScrollSafeOrbit, createStudioRenderer, getSoftParticleTexture, observeInView, thermalColor, type Scene3dFixture } from "./kit";
import { buildFixture, type Model } from "./models/all";
import { palette, perforation } from "./models/materials";
import { placeFixture } from "./models/shop";
import { studioEnvironment, studioLights } from "./models/stage";
import { Crowd, NavGrid, PAIR_MARGIN, STATIC_MARGIN, type Footprint, type LegCandidate, type LegVerdict, type Motion, type Rect, type XZ } from "./models/crowd";
import { Routine, TOOLS, copyPose, makeStaff, newPose, type Cohort, type RoutineStep, type Staff, type WorkKind, type WorkOpts } from "./models/people";

type OverlayMode = "pbr" | "thermal" | "plumes" | "aisles";

export interface HeroController {
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
const LIVE = /RANGE|FRYER|REACHIN|UNDERCOUNTER|WALKIN|SPEEDRACK/;
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
    const liveJoints = new Set<THREE.Object3D>(live ? [...model.joints.values()] : []);
    const toMerge: THREE.Mesh[] = [];
    model.group.traverse((obj) => {
      const m = obj as THREE.Mesh;
      if (!m.isMesh || (m as THREE.InstancedMesh).isInstancedMesh || (m as THREE.SkinnedMesh).isSkinnedMesh) return;
      const g = m.geometry as THREE.BufferGeometry;
      if (!g.getAttribute("position")) return;
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
      toMerge.push(m);
    });
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
    // Let the page breathe between fixtures.
    if (k % 4 === 0) await new Promise((r) => setTimeout(r, 0));
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
  // Obstacles to walking: the fixtures below head height.
  const grid = new NavGrid(
    { minX: 0, minZ: -D, maxX: W, maxZ: 0 },
    fixtures.filter((x) => x.rect).map((x) => x.rect!),
    W * D > 150 ? 0.08 : W * D > 60 ? 0.05 : 0.04,
  );
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
type TaskKind = "saute" | "stir" | "chop" | "plate" | "scrub" | "fry";
const WORK_OF: Record<TaskKind, WorkKind> = { saute: "saute", stir: "stir", chop: "chop", plate: "plate", scrub: "stir", fry: "stand" };
const TOOL_OF: Record<TaskKind, ToolKind | null> = { saute: "skillet", stir: "ladle", chop: "knife", plate: "spoon", scrub: "brush", fry: null };
/** How far each kind of work may reach from where the body stands to its target, metres at p50 scale. */
const REACH: Record<TaskKind, number> = { saute: 0.66, stir: 0.72, chop: 0.6, plate: 0.6, scrub: 0.72, fry: 0.6 };
/** How long a person stays at each kind of work, seconds (min, max). */
const DWELL: Record<TaskKind, [number, number]> = { saute: [9, 16], stir: [7, 12], chop: [9, 16], plate: [6, 10], scrub: [8, 14], fry: [7, 7] };
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
  speed: number;
  routine: Routine;
  steps: RoutineStep[];
  start: { x: number; z: number; yaw: number };
  t0: number;
  /** The committed end of the current routine (its track's last step), seconds. */
  endT: number;
  station: Station;
  workDur: number;
  planned: { steps: RoutineStep[]; start: { x: number; z: number; yaw: number }; t0: number; station: Station; workDur: number; path: XZ[] } | null;
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

/** Footprints, body frame (x to the left, z forward): standing (three discs across the shoulders), walking (one). */
function footprints(st: Staff): { shapes: Footprint[]; rWalk: number } {
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
  return { shapes: [stand, walk], rWalk };
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
      const want = w > 1.8 ? 2 : 1;
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

class Crew {
  readonly group = new THREE.Group();
  readonly agents: Agent[] = [];
  readonly crowd: Crowd;
  readonly stations: Station[];
  readonly lines = new THREE.Group();
  /** Planning results, for the console and for checks. */
  readonly log = { legs: 0, waits: 0, delays: 0, planMs: 0, planMax: 0, searchMax: 0 };
  now = 0;
  private ghosts = new Map<Cohort, Staff>();
  private rnd: () => number;
  private scratch = newPose();
  private lastPrune = 0;

  constructor(
    readonly K: Kitchen,
    detailed: boolean,
    phone: boolean,
  ) {
    this.group.name = "crew";
    this.lines.name = "paths";
    this.crowd = new Crowd(K.grid, CROWD_DT);
    this.rnd = rngSeed(K.fixtures.length * 7919 + Math.round(K.W * 100));
    this.stations = makeStations(K, this.group);
    for (const s of this.stations) {
      const tk = TOOL_OF[s.kind];
      if (tk && !s.tool) {
        s.tool = makeTool(tk);
        s.rest = restFor(s);
      }
      if (s.tool && s.rest) {
        restMatrix(s.rest, s.tool.matrix).decompose(s.tool.position, s.tool.quaternion, s.tool.scale);
        this.group.add(s.tool);
      }
      if (s.kind === "plate" && !s.noPlate) {
        const pl = makePlate();
        pl.position.copy(s.at);
        this.group.add(pl);
        if (s.rest) {
          const bain = makeBain();
          bain.position.copy(s.rest.p).setY(s.at.y);
          this.group.add(bain);
        }
      }
      if (s.baskets) s.baskets.forEach((b, k) => s.fx.model.hooks.set?.(k === 0 ? "lift_l" : "lift_r", b.down ? 0 : 1));
    }
    // Who works here: the spec's head count (3 to 10), its routes' cohorts and areas as hints.
    const spec = K.spec;
    const n = Math.max(Math.min(3, this.stations.length), Math.min(10, spec.staffCount, Math.floor(this.stations.length * 0.7)));
    const routes = spec.staffRoutes;
    const cohorts: Cohort[] = ["p50", "p05", "p95"];
    const taken = new Set<Station>();
    for (let i = 0; i < n; i++) {
      const route = i < routes.length ? routes[i] : undefined;
      const cohort: Cohort = route ? route.cohort : cohorts[i % 3]!;
      const st = makeStaff(cohort, phone ? 1 : 2, { seed: 11 + i * 17, body: detailed ? "detailed" : "procedural", shadows: !phone });
      let ghost = this.ghosts.get(cohort);
      if (!ghost) {
        ghost = makeStaff(cohort, 1, { body: "procedural", shadows: false });
        this.ghosts.set(cohort, ghost);
      }
      const { shapes, rWalk } = footprints(st);
      const a: Agent = {
        id: -1,
        st,
        ghost,
        cohort,
        rWalk,
        rPath: rWalk + STATIC_MARGIN + 0.012,
        dFace: rWalk + STATIC_MARGIN + 0.016,
        shapes,
        speed: COHORT_SPEED[cohort],
        routine: null as unknown as Routine,
        steps: [],
        start: { x: 0, z: 0, yaw: 0 },
        t0: 0,
        endT: 0,
        station: null as unknown as Station,
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
        if (this.crowd.staticRoom(id, sp.x, sp.z, sp.yaw, 0) < STATIC_MARGIN) {
          this.crowd.agents.pop();
          this.crowd.tracks.pop();
          continue;
        }
        a.id = id;
        taken.add(s);
        s.holder = id;
        const dur = 2 + this.rnd() * 6;
        const steps = [this.workStep(s, dur)];
        this.begin(a, sp, steps, CROWD_DT, s, dur, null);
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
  }

  /** Where a person stands to work a station, facing it; null if they cannot stand or reach there. */
  stance(a: Agent, s: Station): { x: number; z: number; yaw: number } | null {
    const x = s.faceX + s.nx * a.dFace;
    const z = s.faceZ + s.nz * a.dFace;
    if (this.K.grid.clearance(x, z) < a.rPath) return null;
    const reach = Math.hypot(s.at.x - x, s.at.z - z);
    if (reach > REACH[s.kind] * a.st.B.s) return null;
    return { x, z, yaw: Math.atan2(-s.nx, -s.nz) };
  }

  private workStep(s: Station, dur: number): RoutineStep {
    const opts: WorkOpts = { at: s.at.clone(), surfaceY: s.at.y };
    if (s.kind === "stir" || s.kind === "scrub") opts.tempo = s.kind === "scrub" ? 1.5 : 0.9 + 0.2 * this.rnd();
    return { kind: "work", work: WORK_OF[s.kind], duration: dur, opts };
  }

  private begin(a: Agent, start: { x: number; z: number; yaw: number }, steps: RoutineStep[], t0: number, s: Station, workDur: number, path: XZ[] | null): void {
    a.start = start;
    a.steps = steps;
    a.t0 = t0;
    a.routine = new Routine(a.st, start, steps, false);
    a.station = s;
    a.workDur = workDur;
    a.path = path;
    a.toolK = 0;
    a.holding = false;
    if (s.kind === "fry" && s.baskets) {
      // Drop the basket that is up, or lift the one that has cooked longest.
      const up = s.baskets.findIndex((b) => !b.down);
      s.fryJob = up >= 0 ? { k: up, lift: false } : { k: s.baskets[0]!.since < s.baskets[1]!.since ? 0 : 1, lift: true };
    }
  }

  /** Sample a routine's root on the ghost body every CROWD_DT from its start, with the footprint each sample shows. */
  private realize(a: Agent, start: { x: number; z: number; yaw: number }, steps: RoutineStep[], workDur: number): Motion {
    const r = new Routine(a.ghost, start, steps, false);
    const n = Math.ceil(r.duration / CROWD_DT) + 1;
    const xs = new Float32Array(n);
    const zs = new Float32Array(n);
    const ys = new Float32Array(n);
    const sh = new Uint8Array(n);
    let still = -1;
    for (let i = 0; i < n; i++) {
      if (still >= 0) {
        xs[i] = xs[still]!;
        zs[i] = zs[still]!;
        ys[i] = ys[still]!;
        continue;
      }
      const state = r.apply(Math.min(i * CROWD_DT, r.duration));
      const p = a.ghost.placement;
      xs[i] = p.x;
      zs[i] = p.z;
      ys[i] = p.yaw;
      sh[i] = state.activity === "walk" ? 1 : 0;
      // Once the closing work has begun the body stands where it is.
      if (state.activity === "work" && i * CROWD_DT > r.duration - workDur + 1e-6) still = i;
    }
    return { n, x: xs, z: zs, yaw: ys, shape: sh };
  }

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
      let w = a.weights.get(s) ?? 0.05;
      // The fryer calls when a basket is up or has been down long enough.
      if (s.kind === "fry" && s.baskets) w *= s.baskets.some((b) => !b.down || b.since < this.now - 25) ? 2.5 : 0.3;
      cands.push({ s, w });
    }
    if (!cands.length) return null;
    let r = this.rnd() * cands.reduce((q, c) => q + c.w, 0);
    for (const c of cands) if ((r -= c.w) <= 0) return c.s;
    return cands[cands.length - 1]!.s;
  }

  private job: {
    a: Agent;
    s: Station;
    k0: number;
    here: { x: number; z: number; yaw: number };
    dur: number;
    cands: (LegCandidate & { steps: RoutineStep[]; path: XZ[] })[];
    it: Generator<void, LegVerdict>;
  } | null = null;

  /**
   * Plan a person's next leg: to the next station by the first path variant, speed and wait that keep clear. The
   * path search runs now; checking the candidates against everyone's committed tracks runs in slices (`update`).
   */
  private plan(a: Agent): void {
    const t0ms = performance.now();
    const s = this.choose(a);
    if (!s) {
      this.extend(a, 2);
      return;
    }
    const goal = this.stance(a, s)!;
    const id = a.id;
    const k0 = Math.max(this.crowd.end(id) + 1, this.crowd.step(this.now + 1));
    const here = this.crowd.poseAt(id, this.crowd.end(id));
    const dyn = this.crowd.standingObstacles(id, k0);
    const dur = DWELL[s.kind][0] + this.rnd() * (DWELL[s.kind][1] - DWELL[s.kind][0]);
    const cands: (LegCandidate & { steps: RoutineStep[]; path: XZ[] })[] = [];
    const paths: XZ[][] = [];
    let variant = 0;
    const found = this.K.grid.findPaths([here.x, here.z], [goal.x, goal.z], a.rPath, [0.75, 1, 0.4], dyn);
    for (const path of found) {
      variant++;
      if (!path || path.length < 2) continue;
      if (paths.some((q) => q.length === path.length && q.every((pt, i) => Math.abs(pt[0] - path[i]![0]) < 1e-3 && Math.abs(pt[1] - path[i]![1]) < 1e-3))) continue;
      paths.push(path);
      for (const v of [a.speed, a.speed * 0.78]) {
        const steps: RoutineStep[] = [{ kind: "walk", path: path.slice(1).map(([x, z]) => [x, z] as [number, number]), speed: v }, { kind: "turn", yaw: goal.yaw }, this.workStep(s, dur)];
        const start = { x: here.x, z: here.z, yaw: here.yaw };
        cands.push({ speed: v, variant, steps, path, realize: () => this.realize(a, start, steps, dur) });
      }
    }
    this.log.searchMax = Math.max(this.log.searchMax, performance.now() - t0ms);
    if (!cands.length) {
      this.log.waits++;
      this.extend(a, 1.5);
      return;
    }
    this.job = { a, s, k0, here, dur, cands, it: this.crowd.planLegSteps(id, k0, cands, 6) };
  }

  /** The planning job's verdict: start the leg when its time comes, or keep working and try again shortly. */
  private finish(verdict: LegVerdict): void {
    const { a, s, k0, here, dur, cands } = this.job!;
    this.job = null;
    if (!verdict.ok) {
      this.log.waits++;
      a.fails++;
      this.extend(a, 1.5);
      return;
    }
    const c = cands.find((q) => q.variant === verdict.variant && q.speed === verdict.speed)!;
    const delaySteps = Math.round(verdict.delay / CROWD_DT);
    const t0 = (k0 + delaySteps) * CROWD_DT;
    // The current work runs on until the leg starts.
    this.stretch(a, t0);
    a.planned = { steps: c.steps, start: { x: here.x, z: here.z, yaw: here.yaw }, t0, station: s, workDur: dur, path: c.path };
    s.holder = a.id;
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
    while (this.job && performance.now() - t0ms < PLAN_BUDGET_MS) {
      const r = this.job.it.next();
      if (r.done) this.finish(r.value);
    }
    const ms = performance.now() - t0ms;
    this.log.planMs += ms;
    this.log.planMax = Math.max(this.log.planMax, ms);
    for (const a of this.agents) {
      if (a.planned && t >= a.planned.t0) {
        const prev = a.station;
        this.putDown(a);
        if (prev.holder === a.id) prev.holder = -1;
        const p = a.planned;
        a.planned = null;
        this.begin(a, p.start, p.steps, p.t0, p.station, p.workDur, p.path);
        a.endT = this.crowd.end(a.id) * CROWD_DT;
        this.drawPath(a);
      }
      const local = t - a.t0;
      const state = a.routine.apply(local);
      const u = state.activity === "work" ? local - (a.routine.duration - a.workDur) : -1;
      if (a.station.kind === "fry") this.fry(a, u);
      else this.tool(a, u, dt);
    }
    if (t - this.lastPrune > 5) {
      this.crowd.prune(this.crowd.step(t));
      this.lastPrune = t;
    }
  }

  private putDown(a: Agent): void {
    const s = a.station;
    if (!s.tool || !s.rest || s.tool.parent === this.group) return;
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
    const s = a.station;
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
    const toB = (v: THREE.Vector3) => {
      const c = Math.cos(yaw);
      const sn = Math.sin(yaw);
      return new THREE.Vector3(c * v.x - sn * v.z, v.y, sn * v.x + c * v.z);
    };
    const rz = toB(rest.z);
    const ry = toB(rest.y);
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
    const s = a.station;
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
    const yawQ = new THREE.Quaternion().setFromAxisAngle(new THREE.Vector3(0, 1, 0), -st.placement.yaw);
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

  private drawPath(a: Agent): void {
    if (a.line) {
      this.lines.remove(a.line);
      a.line.geometry.dispose();
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

  dispose(): void {
    for (const a of this.agents) a.st.dispose();
    for (const g of this.ghosts.values()) g.dispose();
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
    return { setExemplar: async () => {}, getLoadedExemplarId: () => "mall_qsr_400", getDrawCalls: () => 0 };
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
  lights.key.intensity = 2.3;
  lights.hemi.intensity = 0.42;
  lights.rim.intensity = 0.55;

  // The stage fills the hero behind its copy: on wide screens the copy takes the left and the console the bottom,
  // so the picture's centre moves right and up.
  let wide = false;
  const fit = () => {
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
      rushInput.min = String(Math.round(spec.nominalOrdersPerHr * 0.35));
      rushInput.max = String(spec.maxOrdersPerHr);
      ordersPerHr = spec.nominalOrdersPerHr;
      rushInput.value = String(ordersPerHr);
    }
    presets = null;
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
    crew = new Crew(K, !phone, phone);
    scene.add(crew.group, crew.lines);
    presets = framePresets(K, crew, spec);
    simT = 0;
    for (const a of crew.agents) void a.st.ready.then(() => applyOverlays());
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
    renderer.render(scene, camera);
    lastDraws = renderer.info.render.calls;
    if (statusEl) statusEl.textContent = `${spec.shortName} · ${Math.round(K.tris + (crew?.tris() ?? 0)).toLocaleString()} tris · ${lastDraws} draw calls`;
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
  // Point the camera somewhere, for scripted checks.
  (stageEl as unknown as { heroLook: (p: number[], t: number[]) => void }).heroLook = (p, t) => {
    camStartPos.copy(camera.position);
    tgtStartPos.copy(controls.target);
    camGoalPos.set(p[0]!, p[1]!, p[2]!);
    tgtGoalPos.set(t[0]!, t[1]!, t[2]!);
    camAnimT = 0;
  };

  function animate(now: number): void {
    requestAnimationFrame(animate);
    const dt = Math.max(0, Math.min(0.05, (now - lastTs) / 1000));
    frameMs += now - lastTs;
    lastTs = now;
    if (!inView || document.hidden) return;
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
  void loadExemplar(initialSpec.id);
  void activeCameraPreset;
  void wide;

  return {
    setExemplar: loadExemplar,
    getLoadedExemplarId: () => currentSpec.id,
    getDrawCalls: () => lastDraws,
  };
}
