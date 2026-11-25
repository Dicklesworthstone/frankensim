/**
 * CHAPTER 4 - THE COLD CHAIN (#cold). Two live modes, built on the procedural model library (./models). Every moving or
 * glowing thing on screen is read from a kernel's state each frame; nothing here integrates or invents behaviour.
 *
 * Walk-in (`stepWalkinKernel`): NSF-WALKIN-COOLER with its right wall, front-right panel and roof cut away, receiving
 * hot GN 1/1 hotel pans on NSF-SPEEDRACK-20PAN racks.
 *   door swing      <- cfg.doorOpenFrac, set("door"); the strip curtain sways while it is open (the model's tick)
 *   pans and racks  <- cfg.panCount and cfg.panDepthIn. A 2 in pan takes one 3 in slide, a 4 in pan every second
 *                      slide, a 6 in pan every third, so deep pans need more racks, as they would in a real cooler.
 *                      Each pan is half full: the kernel loads 1.9 kg of food per inch of depth.
 *   pan steam       <- state.pan_core_c over state.box_air_c: the water-vapour pressure excess of the food surface
 *                      over the box air (Magnus), 1 at the kernel's initial 85 °C pans in a 3.3 °C box
 *   evaporator fans <- state.evap_cooling_w / cfg.ratedEvapW (the model's tick load)
 *   coil frost      <- state.frost_mass_kg / 4.5 kg, the kernel's frost cap, set("frost")
 *   door thermometer<- state.box_air_c, set("temp")
 *   doorway fog     <- state.infiltration_w: the air exchanged, q / (ρ_box·c_p·ΔT), sets how many puffs leave
 *                      (one per 1/160 m³) and, over the half-doorway each stream uses, how fast they go; its
 *                      density follows the box-kitchen temperature contrast. Cold box air slumps out under the
 *                      neutral plane and fogs where it meets the kitchen; kitchen air comes in over it and fogs inside.
 *   food tint       <- state.pan_core_c, a faint warm cast on the food while it is hot
 * Carrier (`stepCarrierKernel`; hold times from `fullModelHoldTimes`, or from the measured fast-estimate table in
 * sim/carrierFastEstimate.ts when that is on): NSF-CARRIER-INSULATED on a work table, its right wall cut
 * away to show the load and the foamed double wall.
 *   vent cap        <- cfg.ventFraction, set("vent")
 *   pans            <- cfg.foodMassKg, one lidded GN 1/1 pan per 2 kg, set("pans")
 *   vent plume      <- ventFraction × headspace_rh × p_sat(food_temp_c), 1 for a wide-open vent on a saturated
 *                      headspace over food at its initial temperature; rise speed from food-ambient buoyancy
 *   headspace haze  <- state.headspace_rh
 */

import * as THREE from "three";
import { HOLD_HORIZON_MIN, buildFastEstimate, fullModelHoldTimes } from "../sim/carrierFastEstimate";
import {
  HACCP_COLD_LIMIT_C,
  HOT_HOLD_FLOOR_C,
  createInitialCarrierState,
  createInitialWalkinState,
  stepCarrierKernel,
  stepWalkinKernel,
  type CarrierConfig,
  type CarrierState,
  type WalkinConfig,
  type WalkinState,
} from "../sim/physicsKernels";
import { bindCanvasResize, clipFitter, createScrollSafeOrbit, createStudioRenderer, drawStripChart, getSoftParticleTexture, observeInView, shaderGate } from "./kit";
import { buildFixture, type Ext, type Model, type Piece } from "./models/all";
import { speedRackSlides } from "./models/cold";
import { palette } from "./models/materials";
import { vessel } from "./models/ovens";
import { Shop, type Detail } from "./models/shop";
import { hotelPan, type PanFood } from "./models/smallwares";
import { studioEnvironment, studioLights, type Lights } from "./models/stage";

type ColdMode = "walkin" | "carrier";

const clamp = (v: number, lo: number, hi: number) => Math.max(lo, Math.min(hi, v));
const clamp01 = (v: number) => clamp(v, 0, 1);
const smooth = (e0: number, e1: number, x: number) => {
  const t = clamp01((x - e0) / (e1 - e0));
  return t * t * (3 - 2 * t);
};
/** Saturation vapour pressure of water over a wet surface, Pa (Magnus, Alduchov-Eskridge coefficients). */
const psatPa = (tC: number) => 610.94 * Math.exp((17.625 * tC) / (tC + 243.04));
const G = 9.80665;

/* ------------------------------------------------------------------------------------------------------------ */
/* Soft particles: steam, fog, haze. Screen-facing sprites, each with its own size, opacity, colour and spin.      */
/* ------------------------------------------------------------------------------------------------------------ */

const PUFF_VERT = /* glsl */ `
attribute float aSize;
attribute float aAlpha;
attribute vec3 aColor;
attribute float aSpin;
uniform float uScale;
varying float vAlpha;
varying vec3 vColor;
varying float vSpin;
void main() {
  vec4 mv = modelViewMatrix * vec4(position, 1.0);
  gl_Position = projectionMatrix * mv;
  gl_PointSize = clamp(aSize * uScale / max(0.05, -mv.z), 0.0, 480.0);
  vAlpha = aAlpha;
  vColor = aColor;
  vSpin = aSpin;
}`;
const PUFF_FRAG = /* glsl */ `
uniform sampler2D uMap;
varying float vAlpha;
varying vec3 vColor;
varying float vSpin;
void main() {
  vec2 p = gl_PointCoord - 0.5;
  float c = cos(vSpin);
  float s = sin(vSpin);
  p = mat2(c, -s, s, c) * p + 0.5;
  float a = min(1.0, texture2D(uMap, p).a * vAlpha);
  if (a < 0.002) discard;
  gl_FragColor = vec4(vColor, a);
  #include <tonemapping_fragment>
  #include <colorspace_fragment>
}`;

/** Pixels per metre at unit depth, shared by every sprite system: drawing-buffer height × projection[1][1] / 2. */
const puffScale = { value: 400 };

class Puffs {
  readonly points: THREE.Points;
  readonly pos: Float32Array;
  readonly size: Float32Array;
  readonly alpha: Float32Array;
  readonly color: Float32Array;
  readonly spin: Float32Array;
  private readonly geo: THREE.BufferGeometry;
  count = 0;

  constructor(readonly n: number, renderOrder = 4) {
    this.pos = new Float32Array(n * 3);
    this.size = new Float32Array(n);
    this.alpha = new Float32Array(n);
    this.color = new Float32Array(n * 3);
    this.spin = new Float32Array(n);
    const g = new THREE.BufferGeometry();
    g.setAttribute("position", new THREE.BufferAttribute(this.pos, 3).setUsage(THREE.DynamicDrawUsage));
    g.setAttribute("aSize", new THREE.BufferAttribute(this.size, 1).setUsage(THREE.DynamicDrawUsage));
    g.setAttribute("aAlpha", new THREE.BufferAttribute(this.alpha, 1).setUsage(THREE.DynamicDrawUsage));
    g.setAttribute("aColor", new THREE.BufferAttribute(this.color, 3).setUsage(THREE.DynamicDrawUsage));
    g.setAttribute("aSpin", new THREE.BufferAttribute(this.spin, 1).setUsage(THREE.DynamicDrawUsage));
    g.setDrawRange(0, 0);
    this.geo = g;
    const mat = new THREE.ShaderMaterial({
      uniforms: { uScale: puffScale, uMap: { value: getSoftParticleTexture() } },
      vertexShader: PUFF_VERT,
      fragmentShader: PUFF_FRAG,
      transparent: true,
      depthWrite: false,
    });
    this.points = new THREE.Points(g, mat);
    this.points.frustumCulled = false;
    this.points.renderOrder = renderOrder;
  }
  /** Start a frame's writes. */
  begin(): void {
    this.count = 0;
  }
  push(x: number, y: number, z: number, size: number, alpha: number, r: number, g: number, b: number, spin: number): void {
    if (this.count >= this.n || alpha <= 0.002) return;
    const i = this.count++;
    this.pos[i * 3] = x;
    this.pos[i * 3 + 1] = y;
    this.pos[i * 3 + 2] = z;
    this.size[i] = size;
    this.alpha[i] = alpha;
    this.color[i * 3] = r;
    this.color[i * 3 + 1] = g;
    this.color[i * 3 + 2] = b;
    this.spin[i] = spin;
  }
  end(): void {
    this.geo.setDrawRange(0, this.count);
    for (const k of ["position", "aSize", "aAlpha", "aColor", "aSpin"]) {
      const a = this.geo.getAttribute(k) as THREE.BufferAttribute;
      a.clearUpdateRanges();
      a.addUpdateRange(0, this.count * a.itemSize);
      a.needsUpdate = true;
    }
  }
}

/**
 * A pool of particles with a life: each is born somewhere, travels a distance `s` along a path its system defines and
 * dies at the end of it. Births come from a rate (per second) carried over between frames.
 */
class Flock {
  readonly live: Uint8Array;
  readonly s: Float32Array;
  readonly len: Float32Array;
  readonly a: Float32Array;
  readonly b: Float32Array;
  readonly c: Float32Array;
  readonly d: Float32Array;
  readonly src: Int16Array;
  private debt = 0;
  private seed: number;
  constructor(readonly n: number, seed: number) {
    this.live = new Uint8Array(n);
    this.s = new Float32Array(n);
    this.len = new Float32Array(n);
    this.a = new Float32Array(n);
    this.b = new Float32Array(n);
    this.c = new Float32Array(n);
    this.d = new Float32Array(n);
    this.src = new Int16Array(n);
    this.seed = seed >>> 0;
  }
  rand(): number {
    this.seed = (this.seed * 1664525 + 1013904223) >>> 0;
    return this.seed / 4294967296;
  }
  /** Birth `rate × dt` particles (fractions carried over), calling `born(i)` to set each one up. */
  spawn(rate: number, dt: number, born: (i: number) => void): void {
    this.debt = Math.min(this.debt + Math.max(0, rate) * dt, this.n);
    let i = 0;
    while (this.debt >= 1) {
      while (i < this.n && this.live[i]) i++;
      if (i >= this.n) {
        this.debt = 0;
        return;
      }
      this.live[i] = 1;
      this.s[i] = 0;
      born(i);
      this.debt -= 1;
    }
  }
  kill(i: number): void {
    this.live[i] = 0;
  }
  clear(): void {
    this.live.fill(0);
    this.debt = 0;
  }
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Studio set: a quarry-tile floor that fades out into the page, lights, environment.                             */
/* ------------------------------------------------------------------------------------------------------------ */

function fadeFloor(size: number, r0: number, r1: number, tint = 0x8c8480): THREE.Mesh {
  const p = palette();
  const mat = p.tile.clone();
  mat.color.setHex(tint);
  const c = document.createElement("canvas");
  c.width = c.height = 256;
  const x = c.getContext("2d")!;
  const g = x.createRadialGradient(128, 128, (128 * r0) / (size / 2), 128, 128, (128 * r1) / (size / 2));
  g.addColorStop(0, "#fff");
  g.addColorStop(1, "#000");
  x.fillStyle = g;
  x.fillRect(0, 0, 256, 256);
  const am = new THREE.CanvasTexture(c);
  am.channel = 1;
  mat.alphaMap = am;
  mat.transparent = true;
  const geo = new THREE.PlaneGeometry(size, size);
  geo.rotateX(-Math.PI / 2);
  const uv = geo.attributes.uv as THREE.BufferAttribute;
  geo.setAttribute("uv1", uv.clone());
  const pos = geo.attributes.position as THREE.BufferAttribute;
  const k = (p.tile.userData.uvScale as number | undefined) ?? 1;
  for (let i = 0; i < uv.count; i++) uv.setXY(i, pos.getX(i) * k, -pos.getZ(i) * k);
  const m = new THREE.Mesh(geo, mat);
  m.receiveShadow = true;
  m.renderOrder = -1;
  return m;
}

function studioScene(renderer: THREE.WebGLRenderer, centre: THREE.Vector3, span: number, envIntensity: number): { scene: THREE.Scene; lights: Lights } {
  const scene = new THREE.Scene();
  scene.environment = studioEnvironment(renderer);
  scene.environmentIntensity = envIntensity;
  const lights = studioLights(scene, centre, span);
  return { scene, lights };
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Speed-rack loads: an 18 × 26 in sheet pan on the slide carrying a GN 1/1 hotel pan of hot food.                */
/* ------------------------------------------------------------------------------------------------------------ */

/** An 18 × 26 in aluminium sheet pan, its beaded rim on y = 0 and its floor 25 mm below. */
function sheetPan(s: Shop): void {
  const prof: (readonly [number, number])[] = s.fine(2)
    ? [[0.006, -0.025], [0.003, -0.0235], [0.0012, -0.005], [-0.0012, -0.0038], [-0.0026, -0.0018], [-0.0016, 0.0004], [0.0012, 0.0004], [0.0026, -0.002], [0.0036, -0.0236], [0.0066, -0.0246]]
    : [[0.006, -0.025], [0.0012, -0.005], [-0.0026, -0.0018], [0.0012, 0.0004], [0.0036, -0.0236], [0.0066, -0.0246]];
  s.add(palette().alu, vessel(s, 0.457, 0.66, 0.012, prof));
}

/** Hotel pan depth for the kernel's depth in inches: GN 1/1 65, 100 and 150 mm (the 2, 4 and 6 in pans). */
const PAN_DEPTH_M: Record<number, number> = { 2: 0.065, 4: 0.1, 6: 0.15 };
/** The kernel loads 1.9 kg of food per inch of pan depth; a GN 1/1 pan holds about 3.8 L per inch: half full. */
const PAN_FILL = 0.5;
/** Rack slides are 3 in apart; a pan with its sheet pan needs 9 mm clear above it. */
function slidesPerPan(depthM: number): number {
  return Math.ceil((depthM + 0.0015 + 0.008) / 0.076 - 1e-6);
}
const FOODS: PanFood[] = [0, 1, 3, 2];
const RACK_EXT: Ext = { w: 0.521, d: 0.66, h: 1.753 };
const WI_EXT: Ext = { w: 2.438, d: 3.048, h: 2.438 };
/** The carrier's worktable. */
const TABLE_EXT: Ext = { w: 1.219, d: 0.762, h: 0.914 };
/** The room close up on a laptop, a step coarser on a phone. */
const roomDetailFor = (): Detail => (Math.min(window.innerWidth, window.innerHeight) < 600 ? 1 : 2);
const RACK_SLIDES = speedRackSlides(RACK_EXT.h);
/**
 * Where the racks stand once rolled in, world (x, z) of each rack's centre, in the order they fill: two down the right
 * side (clear of the coil's drain line), one mid-room, the last just in from the door. Staggered so the view sees each.
 */
const RACK_SPOTS: readonly (readonly [number, number])[] = [
  [0.78, -0.3],
  [0.78, -1.02],
  [-0.3, -0.45],
  [0.2, 0.62],
];
const MAX_PANS = 25;

/** One food in one depth of pan, instanced across every slide that carries it. */
interface LoadKind {
  meshes: THREE.InstancedMesh[];
  foodMats: THREE.MeshStandardMaterial[];
}

/* ------------------------------------------------------------------------------------------------------------ */
/* The carrier's cut: where its right wall is cut away, and the section the cut shows.                            */
/* ------------------------------------------------------------------------------------------------------------ */

const CAR_EXT: Ext = { w: 0.46, d: 0.635, h: 0.625 };
/** The pieces mountColdChainLab builds, for main.ts to build ahead: the walk-in, a rack for each spot, the carrier's
 * table and the carrier. */
export function coldPieces(): Piece[] {
  return [
    ["NSF-WALKIN-COOLER", WI_EXT, roomDetailFor()],
    ...RACK_SPOTS.map((): Piece => ["NSF-SPEEDRACK-20PAN", RACK_EXT, 1]),
    ["NSF-PREPTABLE-48", TABLE_EXT, 1],
    ["NSF-CARRIER-INSULATED", CAR_EXT, 2],
  ];
}
/** The section plane, in the carrier's frame: just inside the cavity's right wall (x 0.05..0.41). */
const CAR_CUT_X = 0.4;

/**
 * The carrier shell's section at the cut, in (u = -z, v = y) of its frame: the double wall round the cavity, open at
 * the front where the door closes it. `inset` shrinks it from both skins (to leave the skins' own thickness showing).
 */
function shellSection(inset: number): THREE.Shape {
  const [u0, u1, v0, v1, R] = [0.058, 0.635, 0.02, 0.625, 0.034];
  const [c1, cv0, cv1, rr] = [0.601, 0.06, 0.585, 0.026];
  const s = new THREE.Shape();
  const a = u0;
  const b = u1 - inset;
  const lo = v0 + inset;
  const hi = v1 - inset;
  const ci = c1 + inset;
  const ilo = cv0 - inset;
  const ihi = cv1 + inset;
  const Ro = R - inset;
  const ri = rr + inset;
  s.moveTo(a, lo);
  s.lineTo(b - Ro, lo);
  s.quadraticCurveTo(b, lo, b, lo + Ro);
  s.lineTo(b, hi - Ro);
  s.quadraticCurveTo(b, hi, b - Ro, hi);
  s.lineTo(a, hi);
  s.lineTo(a, ihi);
  s.lineTo(ci - ri, ihi);
  s.quadraticCurveTo(ci, ihi, ci, ihi - ri);
  s.lineTo(ci, ilo + ri);
  s.quadraticCurveTo(ci, ilo, ci - ri, ilo);
  s.lineTo(a, ilo);
  s.lineTo(a, lo);
  return s;
}

export function mountColdChainLab(root: HTMLElement): void {
  const stageEl = root.querySelector<HTMLElement>("[data-cold-stage]")!;
  const canvas = root.querySelector<HTMLCanvasElement>("[data-cold-canvas]")!;
  const chartCanvas = root.querySelector<HTMLCanvasElement>("[data-cold-chart]");

  const maybeRenderer = createStudioRenderer(canvas, true);
  if (!maybeRenderer) return;
  const renderer = maybeRenderer;
  renderer.toneMappingExposure = 0.95;
  renderer.localClippingEnabled = true;
  const phone = Math.min(window.innerWidth, window.innerHeight) < 600;
  // The room close up on a laptop, a step coarser on a phone; racks and their loads, many of them, a step coarser still.
  const roomDetail = roomDetailFor();
  const detail: Detail = 1;

  const camera = new THREE.PerspectiveCamera(36, 16 / 10, 0.03, 40);
  const controls = createScrollSafeOrbit(camera, canvas, [0, 1, 0], { minDistance: 0.6, maxDistance: 9, maxPolarAngle: Math.PI * 0.48 });
  const gate = shaderGate(renderer, camera);
  bindCanvasResize(stageEl, renderer, camera, gate);
  let userMoved = false;
  controls.addEventListener("start", () => {
    userMoved = true;
  });

  // =========================================================================
  // WALK-IN COOLER
  // =========================================================================
  const wiSet = studioScene(renderer, new THREE.Vector3(0.2, 1.0, 0), 3.4, 0.5);
  const wiScene = wiSet.scene;
  // A warm kitchen outside, a cold room inside: the key a little warm and lower, the sky fill and the rim cool.
  wiSet.lights.key.intensity = 1.25;
  wiSet.lights.key.color.setHex(0xffeedd);
  wiSet.lights.hemi.color.setHex(0xd4e2ff);
  wiSet.lights.hemi.groundColor.setHex(0x2b2421);
  wiSet.lights.hemi.intensity = 0.32;
  wiSet.lights.rim.intensity = 0.9;
  wiScene.add(fadeFloor(16, 3.0, 7.0, 0x9c6c52));
  // The box's own light: the cool vapour-proof fixture's worth of light inside (the fixture is on the cut-away roof).
  const boxLight = new THREE.PointLight(0xcfe0ff, 7, 0, 2);
  boxLight.position.set(0.1, 2.2, -0.15);
  wiScene.add(boxLight);

  const walkin: Model = buildFixture("NSF-WALKIN-COOLER", WI_EXT, roomDetail);
  // Box centred on the origin, front face toward +z.
  walkin.group.position.set(-WI_EXT.w / 2, 0, WI_EXT.d / 2);
  walkin.hooks.set?.("cutaway", 1);
  {
    // The box's cream panel skin, its own copy here, a shade cooler and greyer: a cold room under its own light.
    const panel = palette().panel;
    const cool = panel.clone();
    cool.color.setHex(0xc4ccd4);
    walkin.group.traverse((o) => {
      const mesh = o as THREE.Mesh;
      if (mesh.isMesh && mesh.material === panel) mesh.material = cool;
    });
  }
  wiScene.add(walkin.group);
  /** Door clear opening, world: x from hinge to latch, the wall's outer and inner faces. */
  const DOOR = { x0: 0.76 - WI_EXT.w / 2, x1: 0.76 + 0.914 - WI_EXT.w / 2, h: 1.981, zOut: WI_EXT.d / 2 - 0.087, zIn: WI_EXT.d / 2 - 0.087 - 0.1016 };
  const DOOR_AREA = 0.914 * 1.981;

  const racks: Model[] = RACK_SPOTS.map(([x, z]) => {
    const m = buildFixture("NSF-SPEEDRACK-20PAN", RACK_EXT, detail);
    m.group.position.set(x - RACK_EXT.w / 2, 0, z + RACK_EXT.d / 2);
    m.hooks.set?.("pans", 0);
    m.hooks.set?.("swivel", 0.25);
    m.group.visible = false;
    wiScene.add(m.group);
    return m;
  });

  // Loads, instanced: one set per (depth, food), built when first needed.
  const loadKinds = new Map<string, LoadKind>();
  const loadGroup = new THREE.Group();
  wiScene.add(loadGroup);
  const paletteMats = new Set<THREE.Material>(Object.values(palette()) as THREE.Material[]);
  function loadKind(depthIn: number, food: PanFood): LoadKind {
    const key = `${depthIn}:${food}`;
    const hit = loadKinds.get(key);
    if (hit) return hit;
    const s = new Shop(detail, `load_${key}`);
    sheetPan(s);
    s.at([0, -0.0236, 0], [0, Math.PI / 2, 0], () => hotelPan(s, "1/1", PAN_DEPTH_M[depthIn] ?? 0.1, PAN_FILL, food));
    const built = s.build(true);
    const cap = Math.ceil(MAX_PANS / FOODS.length);
    const meshes: THREE.InstancedMesh[] = [];
    const foodMats: THREE.MeshStandardMaterial[] = [];
    built.group.traverse((o) => {
      const mesh = o as THREE.Mesh;
      if (!mesh.isMesh) return;
      let mat = mesh.material as THREE.Material;
      const std = mat as THREE.MeshStandardMaterial;
      if (!paletteMats.has(mat) && std.isMeshStandardMaterial && std.metalness < 0.3) {
        // Food: its own copy, so that its warmth can show without touching any other chapter's food.
        const own = (mat as THREE.MeshStandardMaterial).clone();
        foodMats.push(own);
        mat = own;
      }
      mesh.updateWorldMatrix(true, false);
      const geo = mesh.geometry.clone().applyMatrix4(mesh.matrixWorld);
      const inst = new THREE.InstancedMesh(geo, mat, cap);
      inst.count = 0;
      inst.castShadow = mesh.castShadow;
      inst.receiveShadow = mesh.receiveShadow;
      inst.frustumCulled = false;
      meshes.push(inst);
      loadGroup.add(inst);
    });
    const kind = { meshes, foodMats };
    loadKinds.set(key, kind);
    return kind;
  }

  /** Each loaded pan, world: food surface centre and the tier gap above it (0 for the top pan of its rack). */
  interface PanSpot {
    x: number;
    y: number;
    z: number;
    gap: number;
  }
  const panSpots: PanSpot[] = [];
  const tmpM = new THREE.Matrix4();

  function layoutLoads(): void {
    const n = clamp(Math.round(wiCfg.panCount), 1, MAX_PANS);
    const depthIn = wiCfg.panDepthIn <= 2.5 ? 2 : wiCfg.panDepthIn <= 4.5 ? 4 : 6;
    const depthM = PAN_DEPTH_M[depthIn]!;
    const k = slidesPerPan(depthM);
    const perRack = Math.floor((RACK_SLIDES.length - 1) / k) + 1;
    for (const kind of loadKinds.values()) for (const m of kind.meshes) m.count = 0;
    const kinds = FOODS.map((f) => loadKind(depthIn, f));
    const counts = FOODS.map(() => 0);
    panSpots.length = 0;
    racks.forEach((r) => (r.group.visible = false));
    for (let i = 0; i < n; i++) {
      const ri = Math.floor(i / perRack);
      const rack = racks[Math.min(ri, racks.length - 1)]!;
      const spot = RACK_SPOTS[Math.min(ri, RACK_SPOTS.length - 1)]!;
      rack.group.visible = true;
      const tier = (i % perRack) * k;
      const y = RACK_SLIDES[tier]!;
      const f = i % FOODS.length;
      const kind = kinds[f]!;
      tmpM.makeTranslation(spot[0], y, spot[1]);
      for (const m of kind.meshes) m.setMatrixAt(counts[f]!, tmpM);
      counts[f] = (counts[f] ?? 0) + 1;
      const top = i % perRack === perRack - 1 || i === n - 1;
      panSpots.push({ x: spot[0], y: y - 0.0236 + depthM * PAN_FILL, z: spot[1], gap: top ? 0 : k * 0.076 - depthM - 0.002 });
    }
    kinds.forEach((kind, f) => {
      for (const m of kind.meshes) {
        m.count = counts[f]!;
        m.instanceMatrix.needsUpdate = true;
      }
    });
  }

  // Steam off the pans, and the doorway's two streams.
  const steam = new Puffs(phone ? 300 : 460);
  const steamFlock = new Flock(steam.n, 7177);
  wiScene.add(steam.points);
  const doorFog = new Puffs(phone ? 220 : 320);
  const outFlock = new Flock(Math.floor(doorFog.n * 0.6), 3109);
  const inFlock = new Flock(doorFog.n - outFlock.n, 5501);
  wiScene.add(doorFog.points);

  // =========================================================================
  // INSULATED CATERING CARRIER
  // =========================================================================
  const carSet = studioScene(renderer, new THREE.Vector3(0, 1.1, 0), 1.2, 0.7);
  const carScene = carSet.scene;
  carSet.lights.key.intensity = 1.5;
  carSet.lights.hemi.intensity = 0.35;
  carScene.add(fadeFloor(10, 1.4, 3.8, 0x9c6c52));
  const table = buildFixture("NSF-PREPTABLE-48", TABLE_EXT, 1);
  table.group.position.set(-TABLE_EXT.w / 2 - 0.12, 0, TABLE_EXT.d / 2 - 0.05);
  carScene.add(table.group);
  // A fill from the open side so that the cut shows its load.
  const sideFill = new THREE.DirectionalLight(0xfff1e2, 0.9);
  sideFill.position.set(3, 1.6, 0.6);
  carScene.add(sideFill);

  const carrier: Model = buildFixture("NSF-CARRIER-INSULATED", CAR_EXT, 2);
  const CAR_POS = new THREE.Vector3(-CAR_EXT.w / 2, 0.914, CAR_EXT.d / 2 - 0.04);
  carrier.group.position.copy(CAR_POS);
  carScene.add(carrier.group);
  // Cut the shell (not the door, the pans or the vent cap) at the section plane.
  const cutPlane = new THREE.Plane(new THREE.Vector3(-1, 0, 0), CAR_POS.x + CAR_CUT_X);
  for (const o of carrier.group.children) {
    const mesh = o as THREE.Mesh;
    if (!mesh.isMesh) continue;
    const m = (mesh.material as THREE.Material).clone();
    m.clippingPlanes = [cutPlane];
    m.clipShadows = true;
    mesh.material = m;
  }
  {
    // The section: polyethylene skins round the foamed-in-place polyurethane core.
    const skinMat = new THREE.MeshStandardMaterial({ color: 0x26282c, roughness: 0.6, metalness: 0 });
    const foamMat = new THREE.MeshStandardMaterial({ color: 0xd8c690, roughness: 0.92, metalness: 0 });
    const skin = new THREE.Mesh(new THREE.ShapeGeometry(shellSection(0), 12), skinMat);
    const foam = new THREE.Mesh(new THREE.ShapeGeometry(shellSection(0.0035), 12), foamMat);
    for (const [m, dx] of [[skin, 0.0002], [foam, 0.0005]] as const) {
      m.rotation.y = Math.PI / 2;
      m.position.set(CAR_CUT_X + dx, 0, 0);
      m.receiveShadow = true;
      carrier.group.add(m);
    }
  }
  const plume = new Puffs(phone ? 160 : 240);
  const plumeFlock = new Flock(plume.n, 9133);
  carScene.add(plume.points);
  const haze = new Puffs(48, 3);
  carScene.add(haze.points);
  /** Vent centre, world: the boss on the door's face. */
  const VENT = new THREE.Vector3(CAR_POS.x + CAR_EXT.w / 2, CAR_POS.y + 0.52, CAR_POS.z + 0.004);

  // =========================================================================
  // PHYSICS STATE & UI WIRING
  // =========================================================================
  let mode: ColdMode = "walkin";
  let inView = true;

  const wiCfg: WalkinConfig = {
    panCount: 16,
    panDepthIn: 4,
    doorOpenFrac: 0.3,
    kitchenAmbientC: 32,
    ratedEvapW: 3800,
  };
  let wiState: WalkinState = createInitialWalkinState();
  const wiHistBox: number[] = [];
  const wiHistPan: number[] = [];
  const wiInit = createInitialWalkinState();
  const STEAM_REF = psatPa(wiInit.pan_core_c) - psatPa(wiInit.box_air_c);

  function prewarmWalkin(): void {
    wiHistBox.length = 0;
    wiHistPan.length = 0;
    wiHistBox.push(wiState.box_air_c);
    wiHistPan.push(wiState.pan_core_c);
    for (let s = 0; s < 14; s += 1) {
      stepWalkinKernel(wiState, wiCfg, 30.0);
      wiHistBox.push(wiState.box_air_c);
      wiHistPan.push(wiState.pan_core_c);
    }
  }
  prewarmWalkin();

  const carCfg: CarrierConfig = {
    ventFraction: 0.15,
    foodMassKg: 4.5,
    ambientC: 22,
    initialFoodC: 84,
  };
  let carState: CarrierState = createInitialCarrierState(carCfg);
  const carHistFood: number[] = [];
  const carHistTg: number[] = [];

  function prewarmCarrier(): void {
    carHistFood.length = 0;
    carHistTg.length = 0;
    carHistFood.push(carState.food_temp_c);
    carHistTg.push(carState.crust_tg_c);
    for (let s = 0; s < 12; s += 1) {
      stepCarrierKernel(carState, carCfg, 20.0);
      carHistFood.push(carState.food_temp_c);
      carHistTg.push(carState.crust_tg_c);
    }
  }
  prewarmCarrier();

  let surrogateEnabled = false;
  let chartAccum = 0;
  // The fast estimate's table: built from the full model in idle time after the chapter mounts (about 70 ms on a
  // laptop), or at once if the fast estimate is switched on before that.
  let fast: ReturnType<typeof buildFastEstimate> | null = null;
  const fastTable = () => (fast ??= buildFastEstimate());
  const idle = (window as { requestIdleCallback?: (cb: () => void) => number }).requestIdleCallback;
  if (idle) idle(() => void fastTable());
  else setTimeout(() => void fastTable(), 2500);

  // DOM controls
  const modeBtns = Array.from(root.querySelectorAll<HTMLButtonElement>("[data-cold-mode]"));
  const wiControlsEl = root.querySelector<HTMLElement>("[data-walkin-controls]")!;
  const carControlsEl = root.querySelector<HTMLElement>("[data-carrier-controls]")!;

  const wiPansIn = root.querySelector<HTMLInputElement>("[data-wi-pans]");
  const wiDepthBtns = Array.from(root.querySelectorAll<HTMLButtonElement>("[data-wi-depth]"));
  const wiDoorIn = root.querySelector<HTMLInputElement>("[data-wi-door]");
  const wiAmbIn = root.querySelector<HTMLInputElement>("[data-wi-amb]");
  const wiResetBtn = root.querySelector<HTMLButtonElement>("[data-wi-reset]");
  const wiPansVal = root.querySelector<HTMLElement>("[data-wi-pans-val]");
  const wiDoorVal = root.querySelector<HTMLElement>("[data-wi-door-val]");
  const wiAmbVal = root.querySelector<HTMLElement>("[data-wi-amb-val]");

  const carVentIn = root.querySelector<HTMLInputElement>("[data-car-vent]");
  const carMassIn = root.querySelector<HTMLInputElement>("[data-car-mass]");
  const carAmbIn = root.querySelector<HTMLInputElement>("[data-car-amb]");
  const carSurrBtn = root.querySelector<HTMLButtonElement>("[data-car-surrogate]");
  const carOodPresetBtn = root.querySelector<HTMLButtonElement>("[data-car-ood-preset]");
  const carResetBtn = root.querySelector<HTMLButtonElement>("[data-car-reset]");
  const carVentVal = root.querySelector<HTMLElement>("[data-car-vent-val]");
  const carMassVal = root.querySelector<HTMLElement>("[data-car-mass-val]");
  const carAmbVal = root.querySelector<HTMLElement>("[data-car-amb-val]");

  const outGrid = root.querySelector<HTMLElement>("[data-cold-readouts]")!;
  const alertBanner = root.querySelector<HTMLElement>("[data-cold-alert]")!;

  // Framing: each mode's subject, seen from its own side, fitted to the stage's shape: half its width `rx` and half
  // its height `ry` as seen (the walk-in's is wide: its door swings out to the left).
  const VIEWS: Record<ColdMode, { target: THREE.Vector3; dir: THREE.Vector3; rx: number; ry: number; min: number; max: number }> = {
    walkin: { target: new THREE.Vector3(-0.02, 0.92, 0.22), dir: new THREE.Vector3(1.0, 0.62, 0.95), rx: 2.74, ry: 2.1, min: 0.6, max: 9 },
    carrier: { target: new THREE.Vector3(0.0, CAR_POS.y + 0.36, CAR_POS.z - CAR_EXT.d / 2 + 0.05), dir: new THREE.Vector3(0.9, 0.42, 1.0), rx: 0.62, ry: 0.59, min: 0.3, max: 3.5 },
  };
  let framedAspect = 0;
  function frame(): void {
    const v = VIEWS[mode];
    const ty = Math.tan(THREE.MathUtils.degToRad(camera.fov) / 2);
    const dist = Math.max(v.rx / (ty * camera.aspect), v.ry / ty);
    controls.target.copy(v.target);
    camera.position.copy(v.target).addScaledVector(v.dir.clone().normalize(), dist);
    controls.minDistance = v.min;
    controls.maxDistance = Math.max(v.max, dist * 1.2);
    controls.update();
    framedAspect = camera.aspect;
  }

  function setColdMode(next: ColdMode): void {
    mode = next;
    for (const b of modeBtns) b.setAttribute("aria-pressed", String(b.dataset.coldMode === mode));
    wiControlsEl.hidden = mode !== "walkin";
    carControlsEl.hidden = mode !== "carrier";
    userMoved = false;
    frame();
    syncColdReadouts();
  }

  function syncColdReadouts(): void {
    if (mode === "walkin") {
      if (wiPansVal) wiPansVal.textContent = `${wiCfg.panCount} pans`;
      if (wiDoorVal) wiDoorVal.textContent = `${Math.round(wiCfg.doorOpenFrac * 100)}% open`;
      if (wiAmbVal) wiAmbVal.textContent = `${Math.round(wiCfg.kitchenAmbientC)} °C`;

      const boxF = (wiState.box_air_c * 9) / 5 + 32;
      outGrid.innerHTML = `
        <div class="stat-cell"><span>Air in the box</span><b class="num">${wiState.box_air_c.toFixed(1)} °C (${boxF.toFixed(1)} °F)</b></div>
        <div class="stat-cell"><span>Centre of the pans</span><b class="num">${wiState.pan_core_c.toFixed(1)} °C</b></div>
        <div class="stat-cell"><span>Heat from the pans</span><b class="num">${(wiState.pan_load_w / 1000).toFixed(2)} kW</b></div>
        <div class="stat-cell"><span>Warm air at the door</span><b class="num">${(wiState.infiltration_w / 1000).toFixed(2)} kW</b></div>
        <div class="stat-cell"><span>Coil cooling</span><b class="num">${(wiState.evap_cooling_w / 1000).toFixed(2)} kW</b></div>
        <div class="stat-cell"><span>Frost on the coil</span><b class="num">${wiState.frost_mass_kg.toFixed(2)} kg (coil at ${Math.round(wiState.frost_ua_factor * 100)}%)</b></div>
      `;
      if (wiState.box_excursion_41f) {
        alertBanner.dataset.state = "danger";
        alertBanner.textContent = `Above 41 °F (5 °C): the air in the box reached ${wiState.peak_box_air_c.toFixed(1)} °C. The hot pans and the open door brought in more heat than the coil could take out.`;
      } else {
        alertBanner.dataset.state = "ok";
        alertBanner.textContent = `Holding at or below 41 °F (5 °C), with the warm air through the doorway and the frost on the coil both counted.`;
      }
      if (chartCanvas) {
        drawStripChart(
          chartCanvas,
          [
            { label: "Box °C", color: "#38bdf8", values: wiHistBox },
            { label: "Pan °C", color: "#f97316", values: wiHistPan },
          ],
          {
            minY: 0,
            maxY: 90,
            unit: "°C",
            threshold: { value: HACCP_COLD_LIMIT_C, color: "#fb7185", label: "41 °F (5 °C) limit" },
          },
        );
      }
    } else {
      if (carVentVal) carVentVal.textContent = `${Math.round(carCfg.ventFraction * 100)}% ${carCfg.ventFraction < 0.2 ? "(Sealed)" : "(Vented)"}`;
      if (carMassVal) carMassVal.textContent = `${(carCfg.foodMassKg ?? 4.5).toFixed(1)} kg`;
      if (carAmbVal) carAmbVal.textContent = `${Math.round(carCfg.ambientC ?? 22)} °C`;

      // The two hold times: from the fast estimate when it is on and will answer, otherwise from the full model.
      const est = surrogateEnabled ? fastTable().estimate(carCfg) : null;
      const times = est?.inRange ? est.times : fullModelHoldTimes(carCfg);
      const crispPct = Math.round(carState.crispness_index * 100);
      const foodF = (carState.food_temp_c * 9) / 5 + 32;
      const hold = (min: number) => (min >= HOLD_HORIZON_MIN ? `over ${HOLD_HORIZON_MIN / 60} h` : `${min.toFixed(0)} min`);

      outGrid.innerHTML = `
        <div class="stat-cell"><span>Food</span><b class="num">${carState.food_temp_c.toFixed(1)} °C (${foodF.toFixed(0)} °F)</b></div>
        <div class="stat-cell"><span>Crust softening point</span><b class="num">${carState.crust_tg_c.toFixed(1)} °C</b></div>
        <div class="stat-cell"><span>Crispness</span><b class="num">${crispPct}%</b></div>
        <div class="stat-cell"><span>Humidity in the box</span><b class="num">${Math.round(carState.headspace_rh * 100)}%</b></div>
        <div class="stat-cell"><span>Soggy after packing</span><b class="num">${hold(times.soggyMin)}</b></div>
        <div class="stat-cell"><span>Below 135 °F after packing</span><b class="num">${hold(times.below135Min)}</b></div>
      `;

      // First the state of the food now (a safety warning is never hidden), then where the hold times came from.
      let state: "ok" | "warn" = "ok";
      let now: string;
      if (carState.food_temp_c < HOT_HOLD_FLOOR_C) {
        state = "warn";
        now = `Below 135 °F (57 °C): the food is no longer safely hot-held.`;
      } else if (carState.crispness_index <= 0.35) {
        state = "warn";
        now = `Soggy: steam trapped in the box has softened the crust. Open the vent to trade heat for crispness.`;
      } else if (times.soggyMin < HOLD_HORIZON_MIN) {
        now = `Crisp and above 135 °F, ${(carState.time_s / 60).toFixed(0)} min after packing. The steam in the box is softening the crust: soggy ${hold(times.soggyMin)} after packing.`;
      } else {
        now = `Crisp and above 135 °F, ${(carState.time_s / 60).toFixed(0)} min after packing, and at this vent setting the crust stays crisp for over ${HOLD_HORIZON_MIN / 60} h.`;
      }
      let source = "";
      if (est && !est.inRange) {
        state = "warn";
        source =
          est.reason === "range"
            ? ` Out of range for the fast estimate, which covers 15–32 °C outside, 2–8 kg of food and food packed at 75–92 °C, so the full model gives the times.`
            : ` Too close to the vent setting where the box stops going soggy for the fast estimate to answer, so the full model gives the times.`;
      } else if (est?.inRange) {
        source = ` Times from the fast estimate: checked against the full model at ${est.checkedPoints.toLocaleString("en-US")} points across its range, its worst error was ${est.boundMin.soggyMin.toFixed(1)} min for soggy and ${est.boundMin.below135Min.toFixed(1)} min for 135 °F.`;
      }
      alertBanner.dataset.state = state;
      alertBanner.textContent = now + source;

      if (chartCanvas) {
        drawStripChart(
          chartCanvas,
          [
            { label: "Food °C", color: "#f59e0b", values: carHistFood },
            { label: "Crust softening °C", color: "#38bdf8", values: carHistTg },
          ],
          {
            minY: 15,
            maxY: 105,
            unit: "°C",
            threshold: { value: HOT_HOLD_FLOOR_C, color: "#fb7185", label: "135 °F (57 °C)" },
          },
        );
      }
    }
  }

  let loadsDirty = true;
  for (const b of modeBtns) {
    b.addEventListener("click", () => setColdMode(b.dataset.coldMode as ColdMode));
  }
  wiPansIn?.addEventListener("input", () => {
    wiCfg.panCount = Number(wiPansIn.value);
    loadsDirty = true;
    syncColdReadouts();
  });
  for (const b of wiDepthBtns) {
    b.addEventListener("click", () => {
      wiCfg.panDepthIn = Number(b.dataset.wiDepth);
      for (const x of wiDepthBtns) x.setAttribute("aria-pressed", String(Number(x.dataset.wiDepth) === wiCfg.panDepthIn));
      loadsDirty = true;
      syncColdReadouts();
    });
  }
  wiDoorIn?.addEventListener("input", () => {
    wiCfg.doorOpenFrac = Number(wiDoorIn.value) / 100;
    syncColdReadouts();
  });
  wiAmbIn?.addEventListener("input", () => {
    wiCfg.kitchenAmbientC = Number(wiAmbIn.value);
    syncColdReadouts();
  });
  wiResetBtn?.addEventListener("click", () => {
    wiState = createInitialWalkinState();
    wiHistBox.length = 0;
    wiHistPan.length = 0;
    wiHistBox.push(wiState.box_air_c);
    wiHistPan.push(wiState.pan_core_c);
    syncColdReadouts();
  });

  carVentIn?.addEventListener("input", () => {
    carCfg.ventFraction = Number(carVentIn.value) / 100;
    syncColdReadouts();
  });
  carMassIn?.addEventListener("input", () => {
    carCfg.foodMassKg = Number(carMassIn.value);
    syncColdReadouts();
  });
  carAmbIn?.addEventListener("input", () => {
    carCfg.ambientC = Number(carAmbIn.value);
    syncColdReadouts();
  });
  carSurrBtn?.addEventListener("click", () => {
    surrogateEnabled = !surrogateEnabled;
    carSurrBtn.setAttribute("aria-pressed", String(surrogateEnabled));
    syncColdReadouts();
  });
  carOodPresetBtn?.addEventListener("click", () => {
    surrogateEnabled = true;
    carSurrBtn?.setAttribute("aria-pressed", "true");
    carCfg.ambientC = 42;
    carCfg.foodMassKg = 11.5;
    if (carAmbIn) {
      carAmbIn.value = "42";
      carAmbIn.dispatchEvent(new Event("input", { bubbles: true }));
    }
    if (carMassIn) {
      carMassIn.value = "11.5";
      carMassIn.dispatchEvent(new Event("input", { bubbles: true }));
    }
    syncColdReadouts();
  });
  carResetBtn?.addEventListener("click", () => {
    carState = createInitialCarrierState(carCfg);
    carHistFood.length = 0;
    carHistTg.length = 0;
    carHistFood.push(carState.food_temp_c);
    carHistTg.push(carState.crust_tg_c);
    syncColdReadouts();
  });

  setColdMode("walkin");
  observeInView(stageEl, (visible) => {
    inView = visible;
  });

  // =========================================================================
  // PER-FRAME VISUALS, each read from the kernel state
  // =========================================================================
  const steamCol = [1.0, 1.0, 1.0] as const;
  const coldCol = [0.84, 0.92, 1.0] as const;
  const warmCol = [0.98, 0.95, 0.9] as const;
  const warmTint = new THREE.Color();
  /** Doorway fog: one puff per 1/160 m³ of air exchanged; a puff mixes away within a few seconds. */
  const PUFFS_PER_M3 = 160;
  const FOG_AGE = 4.5;

  function drawWalkin(t: number, dt: number): void {
    const s = wiState;
    walkin.hooks.set?.("door", wiCfg.doorOpenFrac);
    walkin.hooks.set?.("frost", clamp01(s.frost_mass_kg / 4.5));
    walkin.hooks.set?.("temp", s.box_air_c);
    walkin.hooks.tick?.(t, dt, clamp01(s.evap_cooling_w / (wiCfg.ratedEvapW ?? 3800)));

    // Steam off each pan's food: vapour pressure excess at the food surface over the box air.
    const drive = clamp01((psatPa(s.pan_core_c) - psatPa(s.box_air_c)) / STEAM_REF);
    const rise = 0.12 + 0.28 * Math.sqrt(clamp01((s.pan_core_c - s.box_air_c) / 80));
    const nP = panSpots.length;
    const F = steamFlock;
    F.spawn(nP * 9 * drive, dt, (i) => {
      const src = Math.floor(F.rand() * nP);
      F.src[i] = src;
      // Start on the food; leave by the nearest open side of the rack under the pan above, then rise.
      const u = F.rand() * 2 - 1;
      const v = F.rand() * 2 - 1;
      F.a[i] = u * 0.15;
      F.b[i] = v * 0.25;
      F.c[i] = F.rand() * Math.PI * 2;
      F.d[i] = F.rand();
      F.len[i] = 0.45 + F.rand() * 0.5;
    });
    steam.begin();
    for (let i = 0; i < F.n; i++) {
      if (!F.live[i]) continue;
      const p = panSpots[F.src[i]!];
      if (!p) {
        F.kill(i);
        continue;
      }
      F.s[i] += dt * rise;
      const sv = F.s[i]!;
      if (sv >= F.len[i]!) {
        F.kill(i);
        continue;
      }
      let x = p.x + F.a[i]!;
      let z = p.z + F.b[i]!;
      let y = p.y + 0.01;
      let out = 0;
      if (p.gap > 0) {
        // Under a pan: drift out sideways (along x: the rack's open sides) through the tier gap.
        const side = F.a[i]! >= 0 ? 1 : -1;
        const exitX = p.x + side * 0.3;
        const lh = Math.abs(exitX - x);
        if (sv < lh) {
          x += side * sv;
          y += Math.min(p.gap * 0.7, sv * 0.3);
        } else {
          out = sv - lh;
          x = exitX + side * out * 0.25;
          y += p.gap * 0.7 + out;
        }
      } else {
        out = sv;
        y += out;
      }
      const ph = F.c[i]!;
      x += Math.sin(t * 1.3 + ph) * 0.03 * out;
      z += Math.cos(t * 1.1 + ph * 1.7) * 0.03 * out - out * 0.12;
      // The ceiling is there even where the view cuts it away: steam reaching it spreads and is gone.
      const ceil = WI_EXT.h - 0.17;
      if (y > ceil - 0.15) y = ceil - 0.15 + (y - ceil + 0.15) * 0.25;
      const life = sv / F.len[i]!;
      const a = smooth(0, 0.12, life) * (1 - smooth(0.5, 1, life)) * (1 - smooth(ceil - 0.25, ceil, y)) * (0.3 + 0.45 * Math.sqrt(drive));
      steam.push(x, y, z, 0.12 + out * 0.7, a, steamCol[0], steamCol[1], steamCol[2], ph + t * 0.3 * (F.d[i]! - 0.5));
    }
    steam.end();

    // Doorway: the kernel's infiltration as air exchanged. Two streams of equal volume, each through half the door.
    const dT = Math.max(0.1, wiCfg.kitchenAmbientC - s.box_air_c);
    const rhoBox = 353.0 / (s.box_air_c + 273.15);
    const vdot = s.infiltration_w / (rhoBox * 1006.0 * dT);
    const open = wiCfg.doorOpenFrac;
    const speed = open > 0.001 ? clamp(vdot / (open * DOOR_AREA * 0.5), 0, 2) : 0;
    const swing = Math.sin(Math.min(Math.PI / 2, open * 1.92));
    const laneW = (DOOR.x1 - DOOR.x0 - 0.06) * Math.min(1, swing * 1.4);
    const neutral = DOOR.h * 0.5;
    const rate = vdot * PUFFS_PER_M3;
    const contrast = clamp01(dT / 25);
    outFlock.spawn(rate * 0.5, dt, (i) => {
      outFlock.a[i] = DOOR.x1 - 0.03 - outFlock.rand() * laneW;
      outFlock.b[i] = 0.03 + outFlock.rand() * neutral * 0.8;
      outFlock.c[i] = outFlock.rand() * Math.PI * 2;
      outFlock.len[i] = 1.4 + outFlock.rand() * 0.9;
      outFlock.d[i] = 0;
    });
    inFlock.spawn(rate * 0.5, dt, (i) => {
      inFlock.a[i] = DOOR.x1 - 0.03 - inFlock.rand() * laneW;
      inFlock.b[i] = neutral * 1.1 + inFlock.rand() * (DOOR.h - neutral * 1.1 - 0.06);
      inFlock.c[i] = inFlock.rand() * Math.PI * 2;
      inFlock.len[i] = 1.3 + inFlock.rand() * 0.6;
      inFlock.d[i] = 0;
    });
    doorFog.begin();
    const dcx = (DOOR.x0 + DOOR.x1) / 2;
    for (let i = 0; i < outFlock.n; i++) {
      if (!outFlock.live[i]) continue;
      outFlock.s[i] += dt * Math.max(0.15, speed);
      outFlock.d[i] += dt;
      const sv = outFlock.s[i]!;
      if (sv >= outFlock.len[i]! || outFlock.d[i]! > FOG_AGE) {
        outFlock.kill(i);
        continue;
      }
      // From just inside the curtain, out under the neutral plane, slumping to the floor and spreading.
      const z = DOOR.zIn - 0.15 + sv;
      const past = Math.max(0, z - DOOR.zOut + 0.08);
      const size = 0.22 + past * 0.5;
      const floorY = size * 0.32;
      const y = floorY + Math.max(0, outFlock.b[i]! - floorY) * Math.exp(-past * 2.4);
      const x = dcx + (outFlock.a[i]! - dcx) * (1 + past * 0.9) + Math.sin(t * 0.9 + outFlock.c[i]!) * 0.04 * past;
      const life = Math.max(sv / outFlock.len[i]!, outFlock.d[i]! / FOG_AGE);
      const a = smooth(0, 0.3, past) * (1 - smooth(0.55, 1, life)) * 0.46 * contrast;
      doorFog.push(x, y, z, size, a, coldCol[0], coldCol[1], coldCol[2], outFlock.c[i]! + t * 0.2);
    }
    for (let i = 0; i < inFlock.n; i++) {
      if (!inFlock.live[i]) continue;
      inFlock.s[i] += dt * Math.max(0.15, speed);
      inFlock.d[i] += dt;
      const sv = inFlock.s[i]!;
      if (sv >= inFlock.len[i]! || inFlock.d[i]! > FOG_AGE) {
        inFlock.kill(i);
        continue;
      }
      // From the kitchen over the neutral plane and in along the ceiling: it fogs as it meets the box air.
      const z = DOOR.zOut + 0.35 - sv;
      const inside = Math.max(0, DOOR.zOut - z);
      const y = inFlock.b[i]! + inside * 0.25;
      const x = inFlock.a[i]! + Math.sin(t * 0.8 + inFlock.c[i]!) * 0.04 * inside;
      const life = Math.max(sv / inFlock.len[i]!, inFlock.d[i]! / FOG_AGE);
      const a = smooth(0, 0.25, inside) * (1 - smooth(0.6, 1, life)) * 0.24 * contrast;
      doorFog.push(x, Math.min(y, WI_EXT.h - 0.2), z, 0.22 + inside * 0.4, a, warmCol[0], warmCol[1], warmCol[2], inFlock.c[i]! - t * 0.2);
    }
    doorFog.end();

    // A thermal cue: the food's surface takes a warmer cast while the kernel's pans are hot (a tint, not a light).
    const hot = smooth(20, 90, s.pan_core_c);
    warmTint.setRGB(0.55, 0.16, 0.03);
    for (const kind of loadKinds.values())
      for (const m of kind.foodMats) {
        m.emissive.copy(warmTint);
        m.emissiveIntensity = 0.12 * hot;
      }
  }

  const PAN_KG = 2;
  function drawCarrier(t: number, dt: number): void {
    const s = carState;
    const vent = clamp01(carCfg.ventFraction);
    carrier.hooks.set?.("vent", vent);
    carrier.hooks.set?.("pans", clamp(Math.ceil((carCfg.foodMassKg ?? 4.5) / PAN_KG - 1e-6), 1, 6));
    const amb = carCfg.ambientC ?? 22;
    const t0 = carCfg.initialFoodC ?? 84;
    // Vapour leaving the vent: the opening × the headspace's vapour pressure, against a wide-open saturated start.
    const flux = clamp01((vent * s.headspace_rh * psatPa(s.food_temp_c)) / psatPa(t0));
    const buoy = Math.sqrt(Math.max(0, (G * 0.06 * (s.food_temp_c - amb)) / (s.food_temp_c + 273.15)));
    const F = plumeFlock;
    F.spawn(flux * 90, dt, (i) => {
      const a = F.rand() * Math.PI * 2;
      const r = 0.012 + F.rand() * 0.018;
      F.a[i] = Math.cos(a) * r;
      F.b[i] = Math.sin(a) * r;
      F.c[i] = F.rand() * Math.PI * 2;
      F.len[i] = 0.35 + F.rand() * 0.35;
    });
    plume.begin();
    const v = 0.1 + 0.7 * buoy;
    for (let i = 0; i < F.n; i++) {
      if (!F.live[i]) continue;
      F.s[i] += dt * v;
      const sv = F.s[i]!;
      if (sv >= F.len[i]!) {
        F.kill(i);
        continue;
      }
      // Out of the slots, bending up as it rises; spreading and fading.
      const jet = 0.05 * (1 - Math.exp(-sv / 0.03));
      const up = Math.max(0, sv - 0.02);
      const ph = F.c[i]!;
      const x = VENT.x + F.a[i]! * (1 + up * 4) + Math.sin(t * 2.1 + ph + sv * 9) * 0.012 * up * 6;
      const y = VENT.y + F.b[i]! + up * 0.95;
      const z = VENT.z + jet + Math.cos(t * 1.7 + ph) * 0.01 * up * 4;
      const life = sv / F.len[i]!;
      const a = smooth(0, 0.08, life) * (1 - smooth(0.35, 1, life)) * (0.16 + 0.2 * flux);
      plume.push(x, y, z, 0.03 + up * 0.5, a, steamCol[0], steamCol[1], steamCol[2], ph + t * 0.4);
    }
    plume.end();

    // The headspace: a haze that thickens toward saturation, slowly turning over inside the cavity.
    haze.begin();
    const hz = 0.22 * smooth(0.3, 1.0, s.headspace_rh);
    const cx0 = CAR_POS.x + 0.07;
    const cx1 = CAR_POS.x + CAR_CUT_X - 0.04;
    for (let i = 0; i < haze.n; i++) {
      const fx = ((i * 0.618) % 1);
      const fy = ((i * 0.381 + 0.13) % 1);
      const fz = ((i * 0.737 + 0.41) % 1);
      const x = cx0 + (cx1 - cx0) * fx + Math.sin(t * 0.21 + i) * 0.02;
      const y = CAR_POS.y + 0.08 + 0.48 * ((fy + t * 0.012 * (1 + (i % 3))) % 1);
      const z = CAR_POS.z - 0.58 + 0.5 * fz + Math.cos(t * 0.17 + i * 1.3) * 0.02;
      haze.push(x, y, z, 0.16, hz, steamCol[0], steamCol[1], steamCol[2], i + t * 0.05);
    }
    haze.end();
  }

  const dbs = new THREE.Vector2();
  let lastTs = performance.now();
  let tAnim = 0;
  // What the last frame drew, for checking the budget from a console.
  (window as unknown as Record<string, unknown>).__KITCHEN_COLD__ = {
    stats: () => ({ mode, triangles: renderer.info.render.triangles, calls: renderer.info.render.calls, pans: panSpots.length, puffs: steam.count + doorFog.count + plume.count + haze.count }),
  };

  const fitWalkin = clipFitter(camera, wiScene, 0.03);
  const fitCarrier = clipFitter(camera, carScene, 0.03);
  function animate(now: number): void {
    requestAnimationFrame(animate);
    const dt = Math.max(0, Math.min(0.05, (now - lastTs) / 1000));
    lastTs = now;
    if (!inView || document.hidden) return;
    tAnim += dt;

    if (mode === "walkin") {
      if (wiState.time_s < 14400) stepWalkinKernel(wiState, wiCfg, dt * 35.0);
    } else {
      if (carState.time_s < 10800) stepCarrierKernel(carState, carCfg, dt * 30.0);
    }

    chartAccum += dt;
    if (chartAccum >= 0.3) {
      chartAccum = 0;
      if (mode === "walkin") {
        wiHistBox.push(wiState.box_air_c);
        wiHistPan.push(wiState.pan_core_c);
        if (wiHistBox.length > 90) {
          wiHistBox.shift();
          wiHistPan.shift();
        }
      } else {
        carHistFood.push(carState.food_temp_c);
        carHistTg.push(carState.crust_tg_c);
        if (carHistFood.length > 90) {
          carHistFood.shift();
          carHistTg.shift();
        }
      }
      syncColdReadouts();
    }

    if (Math.abs(camera.aspect - framedAspect) > 0.01 && !userMoved) frame();
    renderer.getDrawingBufferSize(dbs);
    puffScale.value = dbs.y * 0.5 * camera.projectionMatrix.elements[5]!;

    if (mode === "walkin") {
      if (loadsDirty) {
        loadsDirty = false;
        layoutLoads();
      }
      drawWalkin(tAnim, dt);
    } else drawCarrier(tAnim, dt);

    controls.update();
    (mode === "walkin" ? fitWalkin : fitCarrier)();
    // The models run on; the 3D view waits for its shaders (both scenes are compiled at mount).
    if (gate.ready()) renderer.render(mode === "walkin" ? wiScene : carScene, camera);
  }
  requestAnimationFrame(animate);
  void gate.refresh(wiScene);
  void gate.refresh(carScene);
}
