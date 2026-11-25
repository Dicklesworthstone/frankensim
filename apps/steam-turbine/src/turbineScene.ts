// Photoreal Chronometer-grade 3D CAD & Multiphysics Visualization of the Siemens SST-600
// Industrial Steam Turbine Digital Twin.
//
// Rendering & Engineering Architecture:
// - HDR HalfFloat EffectComposer (MSAA WebGL2 + Selective UnrealBloomPass + SMAAPass + NeutralToneMapping)
// - Industrial studio photographic lightbox PMREM (overhead softbox + tall left/right specular strip lights)
// - 7 procedural high-DPI machining PBR normal & roughness maps (CNC lathe-turned forged steel,
//   brushed 12Cr/Ti airfoils, sand-cast CrMoV casing, face-milled bright cutaway cross-section steel,
//   Siemens slate-blue satin enamel, hot-rolled matte-black structural baseplate)
// - Split multi-material cutaway caps (exterior shell vs gleaming face-milled cross-section steel!)
// - GPU onBeforeCompile per-pixel 6-mode continuum field shader reading 64×24 WASM PDE DataTextures
//   (Temperature T(r,z), von Mises Stress σ_vM, Steam Pressure P(z), Mach Ma, LP Wetness / HP Creep)
//   with 10-level scientific iso-contours and exact CIE 1931 Planck blackbody incandescence
// - Rich Siemens CAD mechanical details: 68 3D hex split-joint bolts + counterbores, visible 5-pad
//   bronze/Babbitt tilting-pad journal bearings + Kingsbury thrust collar, inter-stage gland carriers,
//   LP Stage 16 stainless mid-span damping wire (rebuilt live on design parameter changes!), bolted
//   raised-face pipe flanges, 32-bolt exhaust exit ring, cutaway 60 MVA synchronous generator revealing
//   laminated M19 stator core + copper end-windings, and movable 3D fs-aeroac studio microphone
// - 6 GPU-instanced tapered-arrow 3D physics vector fields + interactive 3D raycast surface probe!

import * as THREE from "three";
import { OrbitControls } from "three/examples/jsm/controls/OrbitControls.js";
import { EffectComposer } from "three/examples/jsm/postprocessing/EffectComposer.js";
import { RenderPass } from "three/examples/jsm/postprocessing/RenderPass.js";
import { UnrealBloomPass } from "three/examples/jsm/postprocessing/UnrealBloomPass.js";
import { SMAAPass } from "three/examples/jsm/postprocessing/SMAAPass.js";
import { OutputPass } from "three/examples/jsm/postprocessing/OutputPass.js";
import type { SteamTurbineSim } from "./wasm-pkg/fs_steam_turbine";
import {
  CONTINUUM_STRIDE,
  decodeBladeRows,
  type BladeRowSpec,
  type StageTelemetry,
  type TurbineSnapshot,
} from "./snapshotLayout";
import type { CameraPreset, CutawayMode, ShaderMode } from "./experiments";
import {
  buildBlackbodyLut,
  buildIndustrialMachiningTextures,
  buildIndustrialStudioEnvironment,
  FIELD_MODE_SPECS,
  installContinuumFieldShader,
  type ContinuumShaderUniforms,
  type IndustrialTextures,
} from "./turbineShaders";

export type NavMode = "orbit" | "pan";

export interface VectorFieldFlags {
  velocity: boolean;
  heat: boolean;
  forces: boolean;
  stress: boolean;
  whirl: boolean;
  acoustic: boolean;
}

export interface SurfaceProbeReading {
  partName: string;
  zM: number;
  rM: number;
  tempC: number;
  vmMpa: number;
  pBar: number;
  mach: number;
  qMagKwM2: number;
  wetnessPct: number;
  creepLog10: number;
}

interface CameraTargetState {
  pos: THREE.Vector3;
  target: THREE.Vector3;
}

const CAMERA_PRESETS: Record<CameraPreset, { pos: [number, number, number]; target: [number, number, number] }> = {
  // Signature Siemens SST-600 front-left 3/4 cutaway hero angle with full breathing room between HUD drawers
  overview: { pos: [-4.25, 3.25, 9.55], target: [-0.08, 1.16, 0.0] },
  // Close-up of the iconic giant twisted golden LP Stage 16 titanium wheel + damping wire
  lp_stage: { pos: [-2.45, 1.92, 2.68], target: [-1.18, 1.38, 0.0] },
  // HP steam inlet S-bend, ESV hydraulic actuator, top twin-actuator CV bridge & Curtis control wheel
  hp_valves: { pos: [3.05, 2.55, 3.45], target: [1.32, 1.40, 0.14] },
  // Close-up on twin balance pistons, Inconel labyrinth seal strips, tilting-pad front bearing & HP drum
  seals: { pos: [1.62, 1.96, 2.38], target: [0.85, 1.38, 0.0] },
  // Exhaust diffuser, 32-bolt flange, manhole hatch, internal struts & rear tilting-pad pedestal
  exhaust: { pos: [-4.95, 2.25, 3.35], target: [-1.72, 1.35, 0.0] },
  // High overhead CAD inspection view straight down into the open rotor, blades & labyrinth glands
  top_open: { pos: [0.15, 6.25, 2.65], target: [0.0, 1.35, 0.0] },
};

const CONTINUUM_NZ = 64;
const CONTINUUM_NR = 24;

function makeWasmGeometry(sim: SteamTurbineSim, meshId: number): THREE.BufferGeometry {
  const pos = sim.mesh_positions(meshId);
  const norm = sim.mesh_normals(meshId);
  const uvs = sim.mesh_uvs(meshId);
  const idx = sim.mesh_indices(meshId);
  const capStart = sim.mesh_cap_start(meshId);

  const geo = new THREE.BufferGeometry();
  geo.setAttribute("position", new THREE.BufferAttribute(new Float32Array(pos), 3));
  geo.setAttribute("normal", new THREE.BufferAttribute(new Float32Array(norm), 3));
  if (uvs.length > 0) {
    geo.setAttribute("uv", new THREE.BufferAttribute(new Float32Array(uvs), 2));
  }
  geo.setIndex(new THREE.BufferAttribute(new Uint32Array(idx), 1));

  // Split shell exterior vs cutaway cross-section cap faces so section walls gleam in face-milled steel!
  if (capStart > 0 && capStart < idx.length) {
    geo.addGroup(0, capStart, 0);
    geo.addGroup(capStart, idx.length - capStart, 1);
  }
  geo.computeBoundingSphere();
  return geo;
}

function buildTaperedArrowGeometry(): THREE.BufferGeometry {
  // Sleek tapered engineering vector arrow along +Y (unit length 1.0)
  const shaft = new THREE.CylinderGeometry(0.011, 0.015, 0.68, 10, 1);
  shaft.translate(0, 0.34, 0);
  const head = new THREE.ConeGeometry(0.036, 0.32, 12, 1);
  head.translate(0, 0.84, 0);

  const p1 = shaft.attributes.position!;
  const n1 = shaft.attributes.normal!;
  const p2 = head.attributes.position!;
  const n2 = head.attributes.normal!;
  const i1 = shaft.index!;
  const i2 = head.index!;

  const pos = new Float32Array((p1.count + p2.count) * 3);
  const norm = new Float32Array((n1.count + n2.count) * 3);
  pos.set(p1.array as Float32Array, 0);
  pos.set(p2.array as Float32Array, p1.count * 3);
  norm.set(n1.array as Float32Array, 0);
  norm.set(n2.array as Float32Array, n1.count * 3);

  const idx = new Uint16Array(i1.count + i2.count);
  idx.set(i1.array as Uint16Array, 0);
  for (let k = 0; k < i2.count; k++) idx[i1.count + k] = (i2.array[k] ?? 0) + p1.count;

  const out = new THREE.BufferGeometry();
  out.setAttribute("position", new THREE.BufferAttribute(pos, 3));
  out.setAttribute("normal", new THREE.BufferAttribute(norm, 3));
  out.setIndex(new THREE.BufferAttribute(idx, 1));
  shaft.dispose();
  head.dispose();
  return out;
}

function classifyPartAtZR(z: number, r: number): string {
  if (z < 0.0) return "60 MVA Synchronous Generator";
  if (z >= 0.25 && z <= 0.37 && r < 0.30) return "Kingsbury Axial Thrust Bearing";
  if (z >= 0.36 && z <= 0.65 && r < 0.42) return "Front 5-Pad Tilting-Pad Bearing";
  if (z >= 4.18 && z <= 4.72 && r < 0.44) return "Rear 5-Pad Tilting-Pad Bearing";
  if (z >= 0.87 && z <= 1.20 && r >= 0.23 && r <= 0.42) return "HP Balance-Piston Labyrinth Seals";
  if (z >= 3.66 && z <= 3.92 && r >= 0.27 && r <= 1.06) return "LP Stage 16 Ti-6Al-4V Last Blade";
  if (z >= 0.0 && z <= 4.88 && r < 0.28) return "Forged 30CrMoNiV5-11 Rotor Shaft";
  if (z >= 1.20 && z <= 2.30 && r >= 0.28 && r <= 0.50) return "HP Curtis & Reaction Blading (12Cr)";
  if (z >= 2.30 && z <= 3.16 && r >= 0.28 && r <= 0.62) return "IP Reaction Blading (12Cr)";
  if (z >= 3.16 && z <= 3.66 && r >= 0.28 && r <= 0.88) return "LP Expanding Blade Rows";
  if (z >= 3.66 && r >= 0.42) return "S355 Slate-Blue Axial Exhaust Diffuser";
  if (z >= 0.88 && z <= 3.38 && r >= 0.44 && r <= 0.64) return "GX12CrMoV12 Inner Blade Carrier";
  if (z >= 0.60 && z <= 3.68 && r >= 0.56) return "Cast CrMoV Outer Pressure Casing";
  return "SST-600 Structural Assembly";
}

export class TurbineScene {
  private readonly container: HTMLElement;
  private readonly sim: SteamTurbineSim;
  private readonly renderer: THREE.WebGLRenderer;
  private readonly composer: EffectComposer;
  private readonly bloomPass: UnrealBloomPass;
  private readonly smaaPass: SMAAPass;
  private readonly scene: THREE.Scene;
  private readonly camera: THREE.PerspectiveCamera;
  private readonly controls: OrbitControls;
  private readonly envRT: THREE.WebGLRenderTarget;
  private readonly resizeObserver: ResizeObserver;
  private readonly textures: IndustrialTextures;

  // Continuum GPU texture uniform bundle (64 x 24 x 8 f32 fields + Planck LUT)
  private readonly gridZMin: number;
  private readonly gridZSpan: number;
  private readonly gridRMax: number;
  private readonly fieldAData = new Float32Array(CONTINUUM_NZ * CONTINUUM_NR * 4);
  private readonly fieldBData = new Float32Array(CONTINUUM_NZ * CONTINUUM_NR * 4);
  private readonly fieldATex: THREE.DataTexture;
  private readonly fieldBTex: THREE.DataTexture;
  private readonly blackbodyTex: THREE.DataTexture;
  private readonly continuumUniforms: ContinuumShaderUniforms;

  // Root transformed so axisymmetric (x_a, y_a, z_a) -> world (-z_a + 2.72, y_a + 1.38, x_a)
  private readonly turbineRoot = new THREE.Group();
  private readonly rotorAssembly = new THREE.Group();
  private readonly rotorSpinGroup = new THREE.Group();
  private readonly rotorBladesGroup = new THREE.Group();
  private readonly statorVanesLower = new THREE.Group();
  private readonly statorVanesUpperBack = new THREE.Group();
  private readonly statorVanesUpperFront = new THREE.Group();

  private readonly innerCarrierLower = new THREE.Group();
  private readonly innerCarrierUpperBack = new THREE.Group();
  private readonly innerCarrierUpperFront = new THREE.Group();
  private readonly outerCasingLower = new THREE.Group();
  private readonly outerCasingUpperBack = new THREE.Group();
  private readonly outerCasingUpperFront = new THREE.Group();
  private readonly exhaustLower = new THREE.Group();
  private readonly exhaustUpperBack = new THREE.Group();
  private readonly exhaustUpperFront = new THREE.Group();
  private readonly frontBrgUpperBack = new THREE.Group();
  private readonly frontBrgUpperFront = new THREE.Group();
  private readonly rearBrgUpperBack = new THREE.Group();
  private readonly rearBrgUpperFront = new THREE.Group();
  private readonly generatorGroup = new THREE.Group();
  private readonly generatorOuterUpperFront = new THREE.Group();
  private readonly generatorCouplingSpin = new THREE.Group();
  private readonly skidAndPipesGroup = new THREE.Group();
  private readonly micAssembly = new THREE.Group();
  private readonly probeMarker: THREE.Mesh;

  // Cached airfoil BufferGeometries from WASM (mesh IDs 17..23) for live parametric blade rebuilds
  private readonly airfoilGeos: THREE.BufferGeometry[] = [];
  private readonly probeTargets: THREE.Object3D[] = [];
  private bladeProbeTargets: THREE.Object3D[] = [];

  // Animated hydraulic actuator piston rods
  private esvRodMesh!: THREE.Mesh;
  private cvRod1Mesh!: THREE.Mesh;
  private cvRod2Mesh!: THREE.Mesh;
  private extRodMesh!: THREE.Mesh;

  // Photoreal Siemens industrial PBR materials (with GPU continuum field shader injection!)
  private readonly matRotor: THREE.MeshPhysicalMaterial;
  private readonly matBladesHp: THREE.MeshStandardMaterial;
  private readonly matBladesLp: THREE.MeshStandardMaterial;
  private readonly matStatorVanes: THREE.MeshStandardMaterial;
  private readonly matLabyrinth: THREE.MeshStandardMaterial;
  private readonly matInnerCarrier: THREE.MeshStandardMaterial;
  private readonly matOuterCasing: THREE.MeshStandardMaterial;
  private readonly matExhaustBlue: THREE.MeshPhysicalMaterial;
  private readonly matBearingHousing: THREE.MeshStandardMaterial;
  private readonly matMachinedSection: THREE.MeshStandardMaterial;
  private readonly matExhaustSection: THREE.MeshStandardMaterial;
  private readonly matBronzePad: THREE.MeshStandardMaterial;
  private readonly matCopperCoil: THREE.MeshStandardMaterial;
  private readonly matStatorCore: THREE.MeshStandardMaterial;
  private readonly matSkidBlack: THREE.MeshStandardMaterial;
  private readonly matChrome: THREE.MeshStandardMaterial;
  private readonly matPipeSteel: THREE.MeshStandardMaterial;
  private readonly matGenBlue: THREE.MeshPhysicalMaterial;

  // Subtle physical effects (fine rub sparks + wispy LP moisture condensation only when wet)
  private readonly mistPoints: THREE.Points;
  private readonly mistData: { z: Float32Array; r: Float32Array; phi: Float32Array };
  private readonly sparkPoints: THREE.Points;
  private readonly sparkVel: Float32Array;

  // 6 GPU-instanced tapered-arrow 3D Physics Vector Fields
  private readonly vfVelocityGroup = new THREE.Group();
  private readonly vfVelocityInst: THREE.InstancedMesh;
  private readonly vfHeatGroup = new THREE.Group();
  private readonly vfHeatInst: THREE.InstancedMesh;
  private readonly vfForcesGroup = new THREE.Group();
  private readonly vfForcesInst: THREE.InstancedMesh;
  private readonly vfStressGroup = new THREE.Group();
  private readonly vfStressInst: THREE.InstancedMesh;
  private readonly vfWhirlGroup = new THREE.Group();
  private readonly vfWhirlPinFront: THREE.Mesh;
  private readonly vfWhirlPinRear: THREE.Mesh;
  private readonly vfWhirlTrailLine: THREE.Line;
  private readonly vfOilWedgeFront: THREE.Mesh;
  private readonly vfOilWedgeRear: THREE.Mesh;
  private readonly vfAcousticGroup = new THREE.Group();
  private readonly vfAcousticInst: THREE.InstancedMesh;

  // Visual state & interaction
  private cutawayMode: CutawayMode = "half";
  private shaderMode: ShaderMode = "pbr";
  private navMode: NavMode = "orbit";
  private isoContoursEnabled = true;
  private explodeFrac = 0;
  private visualSlowMo = 0.035;
  private exaggerateOrbit = true;
  private showSteamMist = true;
  private visualTheta = 0;
  private micPos = new THREE.Vector3(2.2, 1.4, 2.4);
  private onProbeCallback: ((reading: SurfaceProbeReading | null) => void) | null = null;
  private onFirstInteractCallback: (() => void) | null = null;
  private hasInteracted = false;
  private readonly raycaster = new THREE.Raycaster();
  private readonly pointerNdc = new THREE.Vector2();
  private lastProbeMs = 0;
  private isPointerDown = false;
  private lastTapMs = 0;
  private lastTapX = 0;
  private lastTapY = 0;
  private readonly heldKeys = new Set<string>();
  private shiftHeld = false;
  private zoomHoldDir: -1 | 0 | 1 = 0;
  private boundKeyDown?: (ev: KeyboardEvent) => void;
  private boundKeyUp?: (ev: KeyboardEvent) => void;
  private boundBlur?: () => void;
  private boundPointerUp?: () => void;

  private camLerp: { active: boolean; goal: CameraTargetState } = {
    active: false,
    goal: { pos: new THREE.Vector3(), target: new THREE.Vector3() },
  };
  private readonly tmpMat = new THREE.Matrix4();
  private readonly tmpQuat = new THREE.Quaternion();
  private readonly tmpPos = new THREE.Vector3();
  private readonly tmpDir = new THREE.Vector3();
  private readonly tmpScale = new THREE.Vector3();
  private readonly tmpColor = new THREE.Color();
  private readonly upVec = new THREE.Vector3(0, 1, 0);

  constructor(container: HTMLElement, sim: SteamTurbineSim) {
    this.container = container;
    this.sim = sim;

    const w = container.clientWidth || 1280;
    const h = container.clientHeight || 760;

    this.renderer = new THREE.WebGLRenderer({
      antialias: true,
      powerPreference: "high-performance",
      alpha: false,
    });
    this.renderer.setPixelRatio(Math.min(window.devicePixelRatio || 1, 2));
    this.renderer.setSize(w, h);
    this.renderer.toneMapping = THREE.NeutralToneMapping;
    this.renderer.toneMappingExposure = 1.08;
    this.renderer.shadowMap.enabled = true;
    this.renderer.shadowMap.type = THREE.PCFShadowMap;
    container.appendChild(this.renderer.domElement);

    this.scene = new THREE.Scene();
    this.scene.background = new THREE.Color("#0b0f16");
    this.scene.fog = new THREE.Fog("#0b0f16", 18, 42);

    // Industrial studio lightbox PMREM (watchScene.ts technique!)
    const pmrem = new THREE.PMREMGenerator(this.renderer);
    pmrem.compileEquirectangularShader();
    const studioEnvScene = buildIndustrialStudioEnvironment();
    this.envRT = pmrem.fromScene(studioEnvScene, 0.025);
    this.scene.environment = this.envRT.texture;
    this.scene.environmentIntensity = 1.05;
    pmrem.dispose();

    this.camera = new THREE.PerspectiveCamera(34, w / h, 0.08, 80);
    this.applyPresetImmediate("overview");

    this.controls = new OrbitControls(this.camera, this.renderer.domElement);
    this.controls.enableDamping = true;
    this.controls.dampingFactor = 0.08;
    this.controls.enablePan = true;
    this.controls.enableZoom = true;
    this.controls.enableRotate = true;
    this.controls.screenSpacePanning = true;
    this.controls.zoomToCursor = true;
    this.controls.rotateSpeed = 0.88;
    this.controls.panSpeed = 1.05;
    this.controls.zoomSpeed = 1.15;
    this.controls.minDistance = 0.65;
    this.controls.maxDistance = 20.0;
    this.controls.minPolarAngle = 0.04;
    this.controls.maxPolarAngle = Math.PI * 0.88;
    this.controls.cursor.set(-0.08, 1.24, 0.0);
    this.controls.maxTargetRadius = 6.5;
    this.controls.cursorStyle = "grab";
    this.controls.target.set(...CAMERA_PRESETS.overview.target);
    this.controls.mouseButtons = {
      LEFT: THREE.MOUSE.ROTATE,
      MIDDLE: THREE.MOUSE.PAN,
      RIGHT: THREE.MOUSE.PAN,
    };
    this.controls.touches = {
      ONE: THREE.TOUCH.ROTATE,
      TWO: THREE.TOUCH.DOLLY_PAN,
    };
    this.controls.addEventListener("start", () => {
      this.camLerp.active = false;
      if (!this.hasInteracted) {
        this.hasInteracted = true;
        this.onFirstInteractCallback?.();
      }
    });
    // Ensure native browser gestures (page scroll, Safari pinch-zoom, right-click menu) never fight 3D CAD controls
    this.renderer.domElement.style.touchAction = "none";
    this.renderer.domElement.addEventListener("contextmenu", (e) => e.preventDefault());

    // HDR HalfFloat MSAA EffectComposer + selective UnrealBloom + SMAA + OutputPass
    const rt = new THREE.WebGLRenderTarget(w, h, {
      type: THREE.HalfFloatType,
      samples: this.renderer.capabilities.isWebGL2 ? 4 : 0,
    });
    this.composer = new EffectComposer(this.renderer, rt);
    this.composer.setPixelRatio(Math.min(window.devicePixelRatio || 1, 2));
    this.composer.setSize(w, h);
    this.composer.addPass(new RenderPass(this.scene, this.camera));
    this.bloomPass = new UnrealBloomPass(new THREE.Vector2(w, h), 0.10, 0.30, 1.85);
    this.composer.addPass(this.bloomPass);
    this.smaaPass = new SMAAPass();
    this.composer.addPass(this.smaaPass);
    this.composer.addPass(new OutputPass());

    // Studio key, rim kicker, warm fill & soft bounce
    const hemi = new THREE.HemisphereLight("#d4e7fa", "#141922", 0.55);
    this.scene.add(hemi);

    const key = new THREE.DirectionalLight("#fff7eb", 2.15);
    key.position.set(5.2, 8.6, 6.2);
    key.castShadow = true;
    key.shadow.mapSize.set(2048, 2048);
    key.shadow.camera.near = 1.5;
    key.shadow.camera.far = 24;
    const d = 5.4;
    key.shadow.camera.left = -d;
    key.shadow.camera.right = d;
    key.shadow.camera.top = d;
    key.shadow.camera.bottom = -d;
    key.shadow.bias = -0.0004;
    key.shadow.normalBias = 0.015;
    this.scene.add(key);

    const rim = new THREE.DirectionalLight("#8cc4ee", 1.10);
    rim.position.set(-6.4, 4.8, -5.4);
    this.scene.add(rim);

    const fill = new THREE.DirectionalLight("#ffddc0", 0.48);
    fill.position.set(1.2, 1.6, 6.8);
    this.scene.add(fill);

    // Clean photographic studio floor plinth (NO arcade blue grid!)
    const plinthGeo = new THREE.CylinderGeometry(8.8, 9.2, 0.14, 72);
    const plinthMat = new THREE.MeshStandardMaterial({
      color: "#11161f",
      roughness: 0.78,
      metalness: 0.18,
    });
    const plinth = new THREE.Mesh(plinthGeo, plinthMat);
    plinth.position.set(-0.2, -0.075, 0);
    plinth.receiveShadow = true;
    this.scene.add(plinth);

    const floorGeo = new THREE.PlaneGeometry(56, 56);
    const floorMat = new THREE.MeshStandardMaterial({
      color: "#0b0f16",
      roughness: 0.92,
      metalness: 0.05,
    });
    const floor = new THREE.Mesh(floorGeo, floorMat);
    floor.rotation.x = -Math.PI / 2;
    floor.position.y = -0.145;
    floor.receiveShadow = true;
    this.scene.add(floor);

    // Subtle machined scale ticks embedded along the dark plinth front edge (engineering ruler!)
    const tickGeo = new THREE.BoxGeometry(0.018, 0.003, 0.18);
    const tickMat = new THREE.MeshBasicMaterial({ color: "#2d3f52" });
    for (let zm = 0; zm <= 5; zm++) {
      const tick = new THREE.Mesh(tickGeo, tickMat);
      tick.position.set(2.72 - zm, 0.001, 2.32);
      this.scene.add(tick);
    }

    // Build procedural industrial machining PBR normal & roughness maps + Planck LUT + DataTextures
    this.textures = buildIndustrialMachiningTextures();
    this.blackbodyTex = buildBlackbodyLut();

    this.fieldATex = new THREE.DataTexture(
      this.fieldAData,
      CONTINUUM_NZ,
      CONTINUUM_NR,
      THREE.RGBAFormat,
      THREE.FloatType,
    );
    this.fieldATex.minFilter = THREE.NearestFilter;
    this.fieldATex.magFilter = THREE.NearestFilter;
    this.fieldATex.needsUpdate = true;

    this.fieldBTex = new THREE.DataTexture(
      this.fieldBData,
      CONTINUUM_NZ,
      CONTINUUM_NR,
      THREE.RGBAFormat,
      THREE.FloatType,
    );
    this.fieldBTex.minFilter = THREE.NearestFilter;
    this.fieldBTex.magFilter = THREE.NearestFilter;
    this.fieldBTex.needsUpdate = true;

    const dims = sim.continuum_grid_dims();
    this.gridZMin = dims[2] ?? -0.15;
    const gridZMax = dims[3] ?? 5.55;
    this.gridZSpan = Math.max(1e-3, gridZMax - this.gridZMin);
    this.gridRMax = dims[4] ?? 1.52;

    this.continuumUniforms = {
      uFieldA: { value: this.fieldATex },
      uFieldB: { value: this.fieldBTex },
      uBlackbody: { value: this.blackbodyTex },
      uFieldMode: { value: 0 },
      uGridBounds: { value: new THREE.Vector4(this.gridZMin, this.gridZSpan, this.gridRMax, 1.0) },
      uTurbineInvMat: { value: new THREE.Matrix4() },
    };

    // Photoreal industrial CAD materials with machining micro-structure & GPU continuum shaders!
    this.matRotor = new THREE.MeshPhysicalMaterial({
      color: "#cfd5dc",
      metalness: 0.88,
      roughness: 0.25,
      clearcoat: 0.08,
      clearcoatRoughness: 0.32,
      normalMap: this.textures.turnedSteelNormal,
      normalScale: new THREE.Vector2(0.36, 0.36),
    });
    this.matBladesHp = new THREE.MeshStandardMaterial({
      color: "#c9bfaa",
      metalness: 0.82,
      roughness: 0.31,
      normalMap: this.textures.bladeBrushNormal,
      normalScale: new THREE.Vector2(0.30, 0.30),
      side: THREE.DoubleSide,
    });
    this.matBladesLp = new THREE.MeshStandardMaterial({
      color: "#d4c192", // Iconic Siemens golden champagne Ti-6Al-4V LP last-stage wheel!
      metalness: 0.82,
      roughness: 0.27,
      normalMap: this.textures.bladeBrushNormal,
      normalScale: new THREE.Vector2(0.32, 0.32),
      side: THREE.DoubleSide,
    });
    this.matStatorVanes = new THREE.MeshStandardMaterial({
      color: "#9aa4af",
      metalness: 0.80,
      roughness: 0.30,
      normalMap: this.textures.bladeBrushNormal,
      normalScale: new THREE.Vector2(0.20, 0.20),
      side: THREE.DoubleSide,
    });
    this.matLabyrinth = new THREE.MeshStandardMaterial({
      color: "#87929e",
      metalness: 0.88,
      roughness: 0.25,
      normalMap: this.textures.turnedSteelNormal,
      normalScale: new THREE.Vector2(0.42, 0.42),
    });
    this.matInnerCarrier = new THREE.MeshStandardMaterial({
      color: "#aab2bc",
      metalness: 0.72,
      roughness: 0.42,
      normalMap: this.textures.castSteelNormal,
      normalScale: new THREE.Vector2(0.38, 0.38),
      roughnessMap: this.textures.castSteelRoughness,
      side: THREE.DoubleSide,
    });
    this.matOuterCasing = new THREE.MeshStandardMaterial({
      color: "#c4cacf",
      metalness: 0.62,
      roughness: 0.44,
      normalMap: this.textures.castSteelNormal,
      normalScale: new THREE.Vector2(0.42, 0.42),
      roughnessMap: this.textures.castSteelRoughness,
      side: THREE.DoubleSide,
    });
    this.matExhaustBlue = new THREE.MeshPhysicalMaterial({
      color: "#527893", // Authentic Siemens Energy slate-blue enamel
      metalness: 0.28,
      roughness: 0.38,
      clearcoat: 0.26,
      clearcoatRoughness: 0.34,
      normalMap: this.textures.siemensEnamelNormal,
      normalScale: new THREE.Vector2(0.14, 0.14),
      side: THREE.DoubleSide,
    });
    this.matBearingHousing = new THREE.MeshStandardMaterial({
      color: "#5a6674",
      metalness: 0.72,
      roughness: 0.36,
      normalMap: this.textures.castSteelNormal,
      normalScale: new THREE.Vector2(0.30, 0.30),
      side: THREE.DoubleSide,
    });
    // Crisp face-milled bright cutaway cross-section steel (contrasts starkly with cast exteriors!)
    this.matMachinedSection = new THREE.MeshStandardMaterial({
      color: "#dee4eb",
      metalness: 0.88,
      roughness: 0.20,
      normalMap: this.textures.milledSectionNormal,
      normalScale: new THREE.Vector2(0.35, 0.35),
      side: THREE.DoubleSide,
    });
    this.matExhaustSection = new THREE.MeshStandardMaterial({
      color: "#a7bfd0",
      metalness: 0.76,
      roughness: 0.26,
      normalMap: this.textures.milledSectionNormal,
      normalScale: new THREE.Vector2(0.30, 0.30),
      side: THREE.DoubleSide,
    });
    this.matBronzePad = new THREE.MeshStandardMaterial({
      color: "#cc9652",
      metalness: 0.86,
      roughness: 0.24,
      normalMap: this.textures.turnedSteelNormal,
      normalScale: new THREE.Vector2(0.22, 0.22),
    });
    this.matCopperCoil = new THREE.MeshStandardMaterial({
      color: "#cf7844",
      metalness: 0.88,
      roughness: 0.24,
      normalMap: this.textures.bladeBrushNormal,
      normalScale: new THREE.Vector2(0.32, 0.32),
    });
    this.matStatorCore = new THREE.MeshStandardMaterial({
      color: "#363e49",
      metalness: 0.82,
      roughness: 0.34,
      normalMap: this.textures.turnedSteelNormal,
      normalScale: new THREE.Vector2(0.45, 0.45),
    });
    this.matSkidBlack = new THREE.MeshStandardMaterial({
      color: "#1b2028",
      metalness: 0.36,
      roughness: 0.64,
      normalMap: this.textures.skidPlateNormal,
      normalScale: new THREE.Vector2(0.28, 0.28),
    });
    this.matChrome = new THREE.MeshStandardMaterial({
      color: "#ebf0f5",
      metalness: 0.96,
      roughness: 0.11,
    });
    this.matPipeSteel = new THREE.MeshStandardMaterial({
      color: "#c8cdd2",
      metalness: 0.68,
      roughness: 0.34,
      normalMap: this.textures.castSteelNormal,
      normalScale: new THREE.Vector2(0.22, 0.22),
    });
    this.matGenBlue = new THREE.MeshPhysicalMaterial({
      color: "#43657c",
      metalness: 0.30,
      roughness: 0.40,
      clearcoat: 0.22,
      clearcoatRoughness: 0.38,
      normalMap: this.textures.siemensEnamelNormal,
      normalScale: new THREE.Vector2(0.14, 0.14),
    });

    // Inject GPU 6-mode continuum field & Planck incandescence shader into all turbine metal materials!
    const fieldMaterials: Array<[THREE.MeshStandardMaterial | THREE.MeshPhysicalMaterial, string]> = [
      [this.matRotor, "rotor"],
      [this.matBladesHp, "blhp"],
      [this.matBladesLp, "bllp"],
      [this.matStatorVanes, "stat"],
      [this.matLabyrinth, "lab"],
      [this.matInnerCarrier, "carr"],
      [this.matOuterCasing, "cas"],
      [this.matExhaustBlue, "exh"],
      [this.matBearingHousing, "brg"],
      [this.matMachinedSection, "sec"],
      [this.matExhaustSection, "exsec"],
      [this.matBronzePad, "pad"],
    ];
    for (const [m, key] of fieldMaterials) {
      installContinuumFieldShader(m, this.continuumUniforms, key);
    }

    // Root orientation: Axisymmetric +z (HP=0 -> LP=5.45) maps to World -X (HP RIGHT -> LP LEFT)!
    this.turbineRoot.rotation.y = -Math.PI / 2;
    this.turbineRoot.position.set(2.72, 1.38, 0.0);
    this.scene.add(this.turbineRoot);
    this.turbineRoot.updateMatrixWorld(true);
    this.continuumUniforms.uTurbineInvMat.value.copy(this.turbineRoot.matrixWorld).invert();

    this.rotorSpinGroup.add(this.rotorBladesGroup);
    this.rotorAssembly.add(this.rotorSpinGroup);
    this.innerCarrierLower.add(this.statorVanesLower);
    this.innerCarrierUpperBack.add(this.statorVanesUpperBack);
    this.innerCarrierUpperFront.add(this.statorVanesUpperFront);

    this.turbineRoot.add(
      this.rotorAssembly,
      this.innerCarrierLower,
      this.innerCarrierUpperBack,
      this.innerCarrierUpperFront,
      this.outerCasingLower,
      this.outerCasingUpperBack,
      this.outerCasingUpperFront,
      this.exhaustLower,
      this.exhaustUpperBack,
      this.exhaustUpperFront,
      this.frontBrgUpperBack,
      this.frontBrgUpperFront,
      this.rearBrgUpperBack,
      this.rearBrgUpperFront,
      this.skidAndPipesGroup,
      this.generatorGroup,
    );

    // 1. Load all 24 FrankenSim WASM meshes (with UVs & split machined section caps!) + 34 NURBS blade rows
    this.buildFromWasmGeometry(sim);

    // 2. Build precision Siemens SST-600 mechanical CAD sub-assemblies
    this.buildTiltingPadBearingsAndSeals();
    this.buildBoltedFlangesAndHardware();
    this.buildBaseplateSkidAndActuators();
    this.buildSteamPipingAndValves();
    this.buildExhaustDiffuserDetails();
    this.buildCoupledGenerator();
    this.buildMeasurementMicrophone();

    // 3. Subtle physical LP condensation mist (wet steam only) & labyrinth rub sparks (rub only)
    const mist = this.buildMistParticles();
    this.mistPoints = mist.points;
    this.mistData = mist.data;
    this.turbineRoot.add(this.mistPoints);

    const sparks = this.buildRubSparks();
    this.sparkPoints = sparks.points;
    this.sparkVel = sparks.vel;
    this.turbineRoot.add(this.sparkPoints);

    // 4. Build all 6 GPU-instanced tapered-arrow 3D physics vector fields
    const arrowGeo = buildTaperedArrowGeometry();
    this.vfVelocityInst = this.buildVelocityVectorField(arrowGeo);
    this.vfHeatInst = this.buildHeatFluxVectorField(arrowGeo);
    this.vfForcesInst = this.buildForcesVectorField(arrowGeo);
    this.vfStressInst = this.buildStressVectorField(arrowGeo);
    const whirl = this.buildWhirlOrbitOverlay();
    this.vfWhirlPinFront = whirl.front;
    this.vfWhirlPinRear = whirl.rear;
    this.vfWhirlTrailLine = whirl.trail;
    this.vfOilWedgeFront = whirl.wedgeFront;
    this.vfOilWedgeRear = whirl.wedgeRear;
    this.vfAcousticInst = this.buildAcousticVectorField(arrowGeo);

    this.turbineRoot.add(
      this.vfVelocityGroup,
      this.vfHeatGroup,
      this.vfForcesGroup,
      this.vfStressGroup,
      this.vfWhirlGroup,
      this.vfAcousticGroup,
    );

    // 5. Interactive 3D hover probe reticle
    this.probeMarker = new THREE.Mesh(
      new THREE.TorusGeometry(0.055, 0.008, 12, 28),
      new THREE.MeshBasicMaterial({ color: "#7dd3fc", depthTest: false, transparent: true, opacity: 0.92 }),
    );
    this.probeMarker.visible = false;
    this.probeMarker.renderOrder = 999;
    this.scene.add(this.probeMarker);

    this.bindPointerProbe();

    this.setVectorFields({
      velocity: false,
      heat: false,
      forces: false,
      stress: false,
      whirl: false,
      acoustic: false,
    });
    this.setCutawayMode("half");

    this.resizeObserver = new ResizeObserver(() => this.handleResize());
    this.resizeObserver.observe(container);
  }

  private buildFromWasmGeometry(sim: SteamTurbineSim): void {
    const addMesh = (
      meshId: number,
      mat: THREE.Material | THREE.Material[],
      parent: THREE.Object3D,
      castShadow = true,
      probeable = true,
    ): THREE.Mesh => {
      const geo = makeWasmGeometry(sim, meshId);
      const mesh = new THREE.Mesh(geo, mat);
      mesh.castShadow = castShadow;
      mesh.receiveShadow = true;
      parent.add(mesh);
      if (probeable) this.probeTargets.push(mesh);
      return mesh;
    };

    // 0: Forged CrMoV watertight rotor shaft (AxisymmetricChart, Euler χ = 2)
    addMesh(0, this.matRotor, this.rotorSpinGroup, true);
    // 1: Annular Inconel balance-piston labyrinth seals
    addMesh(1, this.matLabyrinth, this.innerCarrierLower, false);

    // 2..16: Capped cutaway shells using [exteriorMat, brightMilledSectionMat] multi-material groups!
    const innerPair = [this.matInnerCarrier, this.matMachinedSection];
    const outerPair = [this.matOuterCasing, this.matMachinedSection];
    const exhaustPair = [this.matExhaustBlue, this.matExhaustSection];
    const brgPair = [this.matBearingHousing, this.matMachinedSection];

    addMesh(2, innerPair, this.innerCarrierLower);
    addMesh(3, innerPair, this.innerCarrierUpperBack);
    addMesh(4, innerPair, this.innerCarrierUpperFront);

    addMesh(5, outerPair, this.outerCasingLower);
    addMesh(6, outerPair, this.outerCasingUpperBack);
    addMesh(7, outerPair, this.outerCasingUpperFront);

    addMesh(8, exhaustPair, this.exhaustLower);
    addMesh(9, exhaustPair, this.exhaustUpperBack);
    addMesh(10, exhaustPair, this.exhaustUpperFront);

    addMesh(11, brgPair, this.exhaustLower);
    addMesh(12, brgPair, this.rearBrgUpperBack);
    addMesh(13, brgPair, this.rearBrgUpperFront);

    addMesh(14, brgPair, this.outerCasingLower);
    addMesh(15, brgPair, this.frontBrgUpperBack);
    addMesh(16, brgPair, this.frontBrgUpperFront);

    // 17..23: Cache all 7 NURBS 3D twisted airfoil geometries from fs-rep-nurbs
    for (let id = 17; id <= 23; id++) {
      this.airfoilGeos.push(makeWasmGeometry(sim, id));
    }

    this.rebuildBladeRows(sim);
  }

  /**
   * Build or live-rebuild all 34 instanced NURBS rotor/stator blade rows, integral shrouds,
   * and LP Stage 16 stainless mid-span damping wire (called on boot AND whenever the engineer
   * changes LP blade height or applies an fs-dfo optimized design!).
   */
  public rebuildBladeRows(sim: SteamTurbineSim): void {
    const clearGroup = (grp: THREE.Group) => {
      while (grp.children.length > 0) {
        const ch = grp.children[0]!;
        grp.remove(ch);
        if (ch instanceof THREE.Mesh && !(ch instanceof THREE.InstancedMesh)) {
          ch.geometry.dispose();
        }
      }
    };
    clearGroup(this.rotorBladesGroup);
    clearGroup(this.statorVanesLower);
    clearGroup(this.statorVanesUpperBack);
    clearGroup(this.statorVanesUpperFront);
    this.bladeProbeTargets = [];

    const rows: BladeRowSpec[] = decodeBladeRows(sim.blade_rows_flat());
    const dummy = new THREE.Object3D();
    const buckets = new Map<
      string,
      { airfoilId: number; parent: THREE.Group; mat: THREE.Material; matrices: THREE.Matrix4[] }
    >();

    const pushMatrix = (
      key: string,
      airfoilId: number,
      parent: THREE.Group,
      mat: THREE.Material,
      mat4: THREE.Matrix4,
    ) => {
      let b = buckets.get(key);
      if (!b) {
        b = { airfoilId, parent, mat, matrices: [] };
        buckets.set(key, b);
      }
      b.matrices.push(mat4.clone());
    };

    for (const row of rows) {
      const span = Math.max(0.02, row.rTipM - row.rHubM);
      const chord = Math.max(0.015, row.axialChordM);
      const isLpBlade = row.isRotating && row.airfoilId >= 4;

      // Rotating integral shroud ring band with dual machined radial seal lips
      if (row.isRotating && row.hasShroud) {
        const ringGeo = new THREE.CylinderGeometry(
          row.rTipM + 0.005,
          row.rTipM + 0.005,
          chord * 0.92,
          56,
          1,
          true,
        );
        ringGeo.rotateX(Math.PI / 2);
        const ring = new THREE.Mesh(ringGeo, isLpBlade ? this.matBladesLp : this.matBladesHp);
        ring.position.z = row.zCenterM;
        this.rotorBladesGroup.add(ring);
      }

      // Signature Siemens stainless mid-span damping / lacing wire on free-standing LP Stage 16!
      if (row.isRotating && row.stageIndex === 16) {
        const rWire = row.rHubM + span * 0.62;
        const wireGeo = new THREE.TorusGeometry(rWire, 0.0075, 10, 64);
        const wire = new THREE.Mesh(wireGeo, this.matChrome);
        wire.position.z = row.zCenterM;
        this.rotorBladesGroup.add(wire);

        // Second inner snubber tie ring near 36% span
        const wire2Geo = new THREE.TorusGeometry(row.rHubM + span * 0.36, 0.0055, 8, 64);
        const wire2 = new THREE.Mesh(wire2Geo, this.matChrome);
        wire2.position.z = row.zCenterM - 0.012;
        this.rotorBladesGroup.add(wire2);
      }

      for (let k = 0; k < row.bladeCount; k++) {
        let phi = (k / row.bladeCount) * Math.PI * 2;
        if (phi > Math.PI) phi -= Math.PI * 2;

        const c = Math.cos(phi);
        const s = Math.sin(phi);
        dummy.position.set(row.rHubM * c, row.rHubM * s, row.zCenterM);
        dummy.rotation.set(0, 0, phi - Math.PI / 2);
        dummy.scale.set(chord, span, chord);
        dummy.updateMatrix();

        if (row.isRotating) {
          const mat = isLpBlade ? this.matBladesLp : this.matBladesHp;
          pushMatrix(`rot-${row.airfoilId}-${isLpBlade ? "lp" : "hp"}`, row.airfoilId, this.rotorBladesGroup, mat, dummy.matrix);
        } else {
          let parent: THREE.Group;
          let tag: string;
          if (phi < 0) {
            parent = this.statorVanesLower;
            tag = "lo";
          } else if (phi > Math.PI / 2) {
            parent = this.statorVanesUpperBack;
            tag = "ub";
          } else {
            parent = this.statorVanesUpperFront;
            tag = "uf";
          }
          pushMatrix(`stat-${row.airfoilId}-${tag}`, row.airfoilId, parent, this.matStatorVanes, dummy.matrix);
        }
      }
    }

    for (const bucket of buckets.values()) {
      const geo = this.airfoilGeos[bucket.airfoilId];
      if (!geo || bucket.matrices.length === 0) continue;
      const inst = new THREE.InstancedMesh(geo, bucket.mat, bucket.matrices.length);
      for (let i = 0; i < bucket.matrices.length; i++) {
        inst.setMatrixAt(i, bucket.matrices[i]!);
      }
      inst.instanceMatrix.needsUpdate = true;
      inst.castShadow = true;
      inst.receiveShadow = true;
      bucket.parent.add(inst);
      this.bladeProbeTargets.push(inst);
    }
  }

  /**
   * Visible 5-pad phosphor-bronze/Babbitt tilting-pad journal bearing shoes inside both open pedestals,
   * Kingsbury axial thrust collar + sector shoes, and balance-piston stationary gland seal carrier!
   */
  private buildTiltingPadBearingsAndSeals(): void {
    const addTiltingPads = (
      zCenter: number,
      rJournal: number,
      padLength: number,
      lowerParent: THREE.Group,
      upperBackParent: THREE.Group,
      upperFrontParent: THREE.Group,
    ) => {
      const nPads = 5;
      const padArc = (Math.PI * 2) / nPads * 0.72;
      const padThick = 0.036;
      const shape = new THREE.Shape();
      const rIn = rJournal + 0.008;
      const rOut = rIn + padThick;
      shape.absarc(0, 0, rOut, -padArc * 0.5, padArc * 0.5, false);
      shape.absarc(0, 0, rIn, padArc * 0.5, -padArc * 0.5, true);
      shape.closePath();

      const padGeo = new THREE.ExtrudeGeometry(shape, { depth: padLength, bevelEnabled: true, bevelSize: 0.004, bevelThickness: 0.004, bevelSegments: 2 });
      padGeo.translate(0, 0, -padLength * 0.5);

      for (let p = 0; p < nPads; p++) {
        // Place pad 0 at bottom (-PI/2) so lower half carries the vertical rotor gravity load!
        let phi = -Math.PI * 0.5 + (p / nPads) * Math.PI * 2;
        if (phi > Math.PI) phi -= Math.PI * 2;
        const mesh = new THREE.Mesh(padGeo, this.matBronzePad);
        mesh.rotation.z = phi;
        mesh.position.z = zCenter;
        mesh.castShadow = true;
        if (phi < 0.12) lowerParent.add(mesh);
        else if (phi > Math.PI * 0.5) upperBackParent.add(mesh);
        else upperFrontParent.add(mesh);
      }

      // Oil-guard bronze labyrinth wiper rings flanking the bearing shell
      const guardGeo = new THREE.TorusGeometry(rJournal + 0.022, 0.011, 10, 44, Math.PI);
      for (const dz of [-padLength * 0.62, padLength * 0.62]) {
        const gLo = new THREE.Mesh(guardGeo, this.matBronzePad);
        gLo.rotation.z = Math.PI;
        gLo.position.z = zCenter + dz;
        lowerParent.add(gLo);
      }
    };

    // Front journal bearing (aligned to geometry.rs FRONT_BEARING_Z_M = 0.50 m, r_shaft = 0.128 m)
    addTiltingPads(0.49, 0.128, 0.22, this.outerCasingLower, this.frontBrgUpperBack, this.frontBrgUpperFront);
    // Rear journal bearing (aligned to geometry.rs REAR_BEARING_Z_M = 4.44 m, r_shaft = 0.132 m)
    addTiltingPads(4.44, 0.132, 0.22, this.exhaustLower, this.rearBrgUpperBack, this.rearBrgUpperFront);

    // Kingsbury axial thrust bronze sector shoes flanking the forged rotor collar (z = 0.273..0.345 m, r = 0.214 m)
    const shoeGeo = new THREE.BoxGeometry(0.048, 0.052, 0.016);
    for (let s = 0; s < 8; s++) {
      let phi = (s / 8) * Math.PI * 2;
      if (phi > Math.PI) phi -= Math.PI * 2;
      for (const dz of [-0.044, 0.044]) {
        const shoe = new THREE.Mesh(shoeGeo, this.matBronzePad);
        shoe.position.set(0.176 * Math.cos(phi), 0.176 * Math.sin(phi), 0.309 + dz);
        shoe.rotation.z = phi;
        if (phi < 0) this.outerCasingLower.add(shoe);
        else if (phi > Math.PI * 0.5) this.frontBrgUpperBack.add(shoe);
        else this.frontBrgUpperFront.add(shoe);
      }
    }

    // Precision machined T-slot stationary gland-seal carrier teeth flanking twin balance pistons (z = 0.90..1.18)
    for (let i = 0; i < 4; i++) {
      const z1 = 0.905 + i * 0.042;
      const fin1 = new THREE.Mesh(
        new THREE.TorusGeometry(0.264, 0.0055, 8, 42, Math.PI),
        this.matLabyrinth,
      );
      fin1.rotation.z = Math.PI;
      fin1.position.z = z1;
      this.innerCarrierLower.add(fin1);

      const z2 = 1.055 + i * 0.041;
      const fin2 = new THREE.Mesh(
        new THREE.TorusGeometry(0.310, 0.0055, 8, 42, Math.PI),
        this.matLabyrinth,
      );
      fin2.rotation.z = Math.PI;
      fin2.position.z = z2;
      this.innerCarrierLower.add(fin2);
    }
  }

  /** Horizontal split-joint bolted flange decks with machined bolt counterbores + 3-part hex bolts/nuts/washers. */
  private buildBoltedFlangesAndHardware(): void {
    const segments: Array<{ z0: number; z1: number; rInner: number; rOuter: number }> = [
      { z0: 0.14, z1: 0.72, rInner: 0.34, rOuter: 0.56 },
      { z0: 0.72, z1: 2.38, rInner: 0.56, rOuter: 0.90 },
      { z0: 2.38, z1: 3.45, rInner: 0.66, rOuter: 1.00 },
      { z0: 3.45, z1: 3.68, rInner: 0.86, rOuter: 1.12 },
    ];

    const boltPositions: Array<[number, number, number]> = [];
    for (const seg of segments) {
      const len = seg.z1 - seg.z0;
      const zMid = 0.5 * (seg.z0 + seg.z1);
      const width = seg.rOuter - seg.rInner;
      const xMid = seg.rInner + 0.5 * width;

      const plateGeo = new THREE.BoxGeometry(width, 0.088, len);
      // Top machined sealing land on lower flange gleams in bright face-milled steel!
      const landGeo = new THREE.BoxGeometry(width - 0.012, 0.008, len - 0.01);

      for (const side of [1, -1]) {
        const lowerFlange = new THREE.Mesh(plateGeo, this.matOuterCasing);
        lowerFlange.position.set(side * xMid, -0.044, zMid);
        lowerFlange.castShadow = true;
        lowerFlange.receiveShadow = true;
        this.outerCasingLower.add(lowerFlange);

        const land = new THREE.Mesh(landGeo, this.matMachinedSection);
        land.position.set(side * xMid, -0.001, zMid);
        this.outerCasingLower.add(land);

        const upperFlange = new THREE.Mesh(plateGeo, this.matOuterCasing);
        upperFlange.position.set(side * xMid, 0.044, zMid);
        upperFlange.castShadow = true;
        upperFlange.receiveShadow = true;
        if (side > 0) this.outerCasingUpperFront.add(upperFlange);
        else this.outerCasingUpperBack.add(upperFlange);
      }

      const nBolts = Math.max(3, Math.round(len / 0.132));
      for (let i = 0; i < nBolts; i++) {
        const zb = seg.z0 + ((i + 0.5) / nBolts) * len;
        const xb = seg.rOuter - 0.068;
        boltPositions.push([xb, 0.0, zb]);
        boltPositions.push([-xb, 0.0, zb]);
      }
    }

    const frontBolts = boltPositions.filter((p) => p[0] > 0);
    const backBolts = boltPositions.filter((p) => p[0] < 0);

    // 3-part mechanical fastener: threaded stud + wide washer + machined hex nut!
    const studGeo = new THREE.CylinderGeometry(0.015, 0.015, 0.24, 10);
    const nutGeo = new THREE.CylinderGeometry(0.028, 0.028, 0.036, 6);
    const washerGeo = new THREE.CylinderGeometry(0.034, 0.034, 0.008, 16);
    const boreGeo = new THREE.CylinderGeometry(0.022, 0.022, 0.012, 14);
    const boreMat = new THREE.MeshStandardMaterial({ color: "#2b323b", metalness: 0.78, roughness: 0.45 });

    const frontBores = new THREE.InstancedMesh(boreGeo, boreMat, frontBolts.length);
    const frontUpperStuds = new THREE.InstancedMesh(studGeo, this.matChrome, frontBolts.length);
    const frontUpperNuts = new THREE.InstancedMesh(nutGeo, this.matChrome, frontBolts.length);
    const backStuds = new THREE.InstancedMesh(studGeo, this.matChrome, backBolts.length);
    const backNuts = new THREE.InstancedMesh(nutGeo, this.matChrome, backBolts.length);
    const backWashers = new THREE.InstancedMesh(washerGeo, this.matMachinedSection, backBolts.length);

    frontBolts.forEach(([x, y, z], i) => {
      this.tmpMat.makeTranslation(x, y + 0.001, z);
      frontBores.setMatrixAt(i, this.tmpMat);
      this.tmpMat.makeTranslation(x, y + 0.02, z);
      frontUpperStuds.setMatrixAt(i, this.tmpMat);
      this.tmpMat.makeTranslation(x, y + 0.105, z);
      frontUpperNuts.setMatrixAt(i, this.tmpMat);
    });
    backBolts.forEach(([x, y, z], i) => {
      this.tmpMat.makeTranslation(x, y + 0.02, z);
      backStuds.setMatrixAt(i, this.tmpMat);
      this.tmpMat.makeTranslation(x, y + 0.092, z);
      backWashers.setMatrixAt(i, this.tmpMat);
      this.tmpMat.makeTranslation(x, y + 0.112, z);
      backNuts.setMatrixAt(i, this.tmpMat);
    });

    for (const inst of [frontBores, frontUpperStuds, frontUpperNuts, backStuds, backNuts, backWashers]) {
      inst.instanceMatrix.needsUpdate = true;
    }

    this.outerCasingLower.add(frontBores);
    this.outerCasingUpperFront.add(frontUpperStuds, frontUpperNuts);
    this.outerCasingUpperBack.add(backStuds, backNuts, backWashers);
  }

  /** Heavy matte-black welded Siemens structural baseplate skid, anchor bolts, twin hydraulic accumulators & lube lines. */
  private buildBaseplateSkidAndActuators(): void {
    const g = this.skidAndPipesGroup;
    const beamLen = 4.42;
    const beamZ = 1.98;
    const beamY = -1.14;

    const mainBeamGeo = new THREE.BoxGeometry(0.24, 0.46, beamLen);
    const topCapGeo = new THREE.BoxGeometry(0.34, 0.048, beamLen + 0.08);
    const botCapGeo = new THREE.BoxGeometry(0.38, 0.056, beamLen + 0.10);
    const anchorGeo = new THREE.CylinderGeometry(0.018, 0.018, 0.08, 6);

    for (const side of [1, -1]) {
      const bx = side * 0.90;
      const web = new THREE.Mesh(mainBeamGeo, this.matSkidBlack);
      web.position.set(bx, beamY, beamZ);
      web.castShadow = true;
      web.receiveShadow = true;
      g.add(web);

      const topCap = new THREE.Mesh(topCapGeo, this.matSkidBlack);
      topCap.position.set(bx, beamY + 0.24, beamZ);
      g.add(topCap);

      const botCap = new THREE.Mesh(botCapGeo, this.matSkidBlack);
      botCap.position.set(bx, beamY - 0.23, beamZ);
      g.add(botCap);

      const gusGeo = new THREE.BoxGeometry(0.065, 0.42, 0.036);
      for (let gz = 0.0; gz <= 4.05; gz += 0.65) {
        const gus = new THREE.Mesh(gusGeo, this.matSkidBlack);
        gus.position.set(bx + side * 0.135, beamY, gz);
        g.add(gus);

        // Foundation anchor hex bolts along lower skid flange
        const anc = new THREE.Mesh(anchorGeo, this.matChrome);
        anc.position.set(bx + side * 0.15, beamY - 0.19, gz + 0.30);
        g.add(anc);
      }
    }

    const crossGeo = new THREE.BoxGeometry(2.14, 0.42, 0.26);
    for (const cz of [-0.08, 1.38, 2.86, 4.06]) {
      const cross = new THREE.Mesh(crossGeo, this.matSkidBlack);
      cross.position.set(0, beamY, cz);
      cross.castShadow = true;
      cross.receiveShadow = true;
      g.add(cross);
    }

    // Casing support paws & machined sliding-key pedestals rising from skid beams
    const pedGeo = new THREE.BoxGeometry(0.36, 0.72, 0.48);
    const padPlateGeo = new THREE.BoxGeometry(0.40, 0.04, 0.52);
    for (const side of [1, -1]) {
      for (const pz of [0.42, 3.65]) {
        const ped = new THREE.Mesh(pedGeo, this.matSkidBlack);
        ped.position.set(side * 0.84, -0.58, pz);
        ped.castShadow = true;
        g.add(ped);

        const padPlate = new THREE.Mesh(padPlateGeo, this.matMachinedSection);
        padPlate.position.set(side * 0.84, -0.21, pz);
        g.add(padPlate);
      }
    }

    // Twin horizontal hydraulic control-oil accumulators slung beneath the front-right skid beam
    const cylBodyGeo = new THREE.CylinderGeometry(0.095, 0.095, 0.58, 24);
    cylBodyGeo.rotateX(Math.PI / 2);
    const cylCapGeo = new THREE.CylinderGeometry(0.118, 0.118, 0.06, 24);
    cylCapGeo.rotateX(Math.PI / 2);
    for (const dx of [-0.14, 0.14]) {
      const cx = 0.88 + dx;
      const body = new THREE.Mesh(cylBodyGeo, this.matSkidBlack);
      body.position.set(cx, -1.50, 0.78);
      g.add(body);
      const cap1 = new THREE.Mesh(cylCapGeo, this.matChrome);
      cap1.position.set(cx, -1.50, 0.47);
      g.add(cap1);
      const cap2 = new THREE.Mesh(cylCapGeo, this.matSkidBlack);
      cap2.position.set(cx, -1.50, 1.09);
      g.add(cap2);
    }

    // Stainless hydraulic supply line rising from accumulators toward the foreground ESV servo!
    const hydCurve = new THREE.CatmullRomCurve3([
      new THREE.Vector3(1.02, -1.42, 0.62),
      new THREE.Vector3(1.24, -1.18, 0.54),
      new THREE.Vector3(1.38, -0.62, 0.42),
      new THREE.Vector3(1.46, -0.18, 0.34),
    ]);
    const hydLine = new THREE.Mesh(new THREE.TubeGeometry(hydCurve, 24, 0.016, 10, false), this.matChrome);
    g.add(hydLine);
  }

  /** Helper: build a raised-face bolted pipe flange with 8 machined hex bolts. */
  private makeBoltedFlange(radius: number, thickness: number): THREE.Group {
    const fg = new THREE.Group();
    const ring = new THREE.Mesh(new THREE.CylinderGeometry(radius, radius, thickness, 24), this.matPipeSteel);
    fg.add(ring);
    const boltGeo = new THREE.CylinderGeometry(radius * 0.075, radius * 0.075, thickness * 1.9, 6);
    for (let i = 0; i < 8; i++) {
      const a = (i / 8) * Math.PI * 2;
      const b = new THREE.Mesh(boltGeo, this.matChrome);
      b.position.set(radius * 0.78 * Math.cos(a), 0, radius * 0.78 * Math.sin(a));
      fg.add(b);
    }
    return fg;
  }

  /** Foreground S-bend live-steam inlet + ESV diffuser/actuator, top twin CV rack, extraction servo & bolted vertical bleeds. */
  private buildSteamPipingAndValves(): void {
    const g = this.skidAndPipesGroup;

    // 1. Foreground right S-bend live-steam inlet pipe feeding the Emergency Stop Valve (ESV)
    const esvCurve = new THREE.CatmullRomCurve3([
      new THREE.Vector3(1.52, -1.36, 0.30),
      new THREE.Vector3(1.52, -0.78, 0.30),
      new THREE.Vector3(1.52, -0.22, 0.50),
      new THREE.Vector3(1.48, -0.10, 0.86),
      new THREE.Vector3(1.36, -0.06, 1.16),
      new THREE.Vector3(0.95, -0.10, 1.24),
      new THREE.Vector3(0.68, -0.32, 1.24),
    ]);
    const esvPipe = new THREE.Mesh(new THREE.TubeGeometry(esvCurve, 48, 0.090, 20, false), this.matPipeSteel);
    esvPipe.castShadow = true;
    g.add(esvPipe);

    const esvBody = new THREE.Mesh(new THREE.CylinderGeometry(0.150, 0.150, 0.38, 28), this.matPipeSteel);
    esvBody.rotation.x = Math.PI / 2;
    esvBody.position.set(1.50, -0.10, 0.68);
    g.add(esvBody);

    // Bolted flanges on ESV bonnet & inlet collar
    const esvFlangeF = this.makeBoltedFlange(0.182, 0.042);
    esvFlangeF.rotation.x = Math.PI / 2;
    esvFlangeF.position.set(1.50, -0.10, 0.49);
    g.add(esvFlangeF);

    const esvFlangeR = this.makeBoltedFlange(0.182, 0.042);
    esvFlangeR.rotation.x = Math.PI / 2;
    esvFlangeR.position.set(1.50, -0.10, 0.87);
    g.add(esvFlangeR);

    const esvActuator = new THREE.Mesh(new THREE.CylinderGeometry(0.090, 0.090, 0.48, 22), this.matPipeSteel);
    esvActuator.rotation.x = Math.PI / 2;
    esvActuator.position.set(1.50, -0.10, 0.20);
    g.add(esvActuator);

    const esvBlackCap = new THREE.Mesh(new THREE.CylinderGeometry(0.094, 0.094, 0.14, 22), this.matSkidBlack);
    esvBlackCap.rotation.x = Math.PI / 2;
    esvBlackCap.position.set(1.50, -0.10, -0.07);
    g.add(esvBlackCap);

    this.esvRodMesh = new THREE.Mesh(new THREE.CylinderGeometry(0.026, 0.026, 0.24, 14), this.matChrome);
    this.esvRodMesh.rotation.x = Math.PI / 2;
    this.esvRodMesh.position.set(1.50, -0.10, 0.46);
    g.add(this.esvRodMesh);

    const inletCollar = this.makeBoltedFlange(0.165, 0.052);
    inletCollar.position.set(0.95, -0.26, 1.24);
    g.add(inletCollar);

    // 2. Top-rear twin HP Control Valve Actuator Bridge
    const bridgeGroup = new THREE.Group();
    const topPipeCurve = new THREE.CatmullRomCurve3([
      new THREE.Vector3(-0.62, 0.46, 1.38),
      new THREE.Vector3(-0.78, 0.92, 1.38),
      new THREE.Vector3(-0.78, 1.38, 1.44),
      new THREE.Vector3(-0.66, 1.50, 1.72),
    ]);
    const topPipe = new THREE.Mesh(new THREE.TubeGeometry(topPipeCurve, 32, 0.096, 18, false), this.matPipeSteel);
    topPipe.castShadow = true;
    bridgeGroup.add(topPipe);

    const topFlange = this.makeBoltedFlange(0.148, 0.042);
    topFlange.position.set(-0.78, 1.12, 1.40);
    bridgeGroup.add(topFlange);

    const stanchionGeo = new THREE.BoxGeometry(0.072, 1.04, 0.092);
    for (const sz of [1.45, 1.85]) {
      const st = new THREE.Mesh(stanchionGeo, this.matSkidBlack);
      st.position.set(-0.62, 1.02, sz);
      bridgeGroup.add(st);
    }

    const cvCylGeo = new THREE.CylinderGeometry(0.088, 0.088, 0.56, 22);
    cvCylGeo.rotateZ(Math.PI / 2);
    const cvCapGeo = new THREE.CylinderGeometry(0.094, 0.094, 0.12, 22);
    cvCapGeo.rotateZ(Math.PI / 2);
    const cvRodGeo = new THREE.CylinderGeometry(0.025, 0.025, 0.30, 14);
    cvRodGeo.rotateZ(Math.PI / 2);

    const cv1 = new THREE.Mesh(cvCylGeo, this.matPipeSteel);
    cv1.position.set(-0.88, 1.50, 1.72);
    const cv1Cap = new THREE.Mesh(cvCapGeo, this.matSkidBlack);
    cv1Cap.position.set(-1.19, 1.50, 1.72);
    this.cvRod1Mesh = new THREE.Mesh(cvRodGeo, this.matChrome);
    this.cvRod1Mesh.position.set(-0.56, 1.50, 1.72);

    const cv2 = new THREE.Mesh(cvCylGeo, this.matPipeSteel);
    cv2.position.set(-0.82, 1.22, 1.38);
    const cv2Cap = new THREE.Mesh(cvCapGeo, this.matSkidBlack);
    cv2Cap.position.set(-1.12, 1.22, 1.38);
    this.cvRod2Mesh = new THREE.Mesh(cvRodGeo, this.matChrome);
    this.cvRod2Mesh.position.set(-0.50, 1.22, 1.38);

    bridgeGroup.add(cv1, cv1Cap, this.cvRod1Mesh, cv2, cv2Cap, this.cvRod2Mesh);
    this.outerCasingUpperBack.add(bridgeGroup);

    // 3. Mid-casing top controlled 12-bar extraction valve dome + slanted hydraulic servo
    const extDomeGroup = new THREE.Group();
    const dome = new THREE.Mesh(new THREE.CylinderGeometry(0.24, 0.32, 0.50, 28), this.matOuterCasing);
    dome.position.set(-0.26, 0.72, 2.52);
    extDomeGroup.add(dome);

    const domeFlange = this.makeBoltedFlange(0.30, 0.048);
    domeFlange.position.set(-0.26, 0.94, 2.52);
    extDomeGroup.add(domeFlange);

    const extServo = new THREE.Mesh(new THREE.CylinderGeometry(0.11, 0.11, 0.62, 22), this.matPipeSteel);
    extServo.rotation.x = Math.PI * 0.34;
    extServo.position.set(-0.28, 1.02, 2.92);
    extDomeGroup.add(extServo);

    const extCap = new THREE.Mesh(new THREE.CylinderGeometry(0.116, 0.116, 0.14, 22), this.matSkidBlack);
    extCap.rotation.x = Math.PI * 0.34;
    extCap.position.set(-0.28, 1.16, 3.24);
    extDomeGroup.add(extCap);

    this.extRodMesh = new THREE.Mesh(new THREE.CylinderGeometry(0.028, 0.028, 0.32, 14), this.matChrome);
    this.extRodMesh.rotation.x = Math.PI * 0.34;
    this.extRodMesh.position.set(-0.28, 0.88, 2.66);
    extDomeGroup.add(this.extRodMesh);
    this.outerCasingUpperBack.add(extDomeGroup);

    // 4. Five vertical bolted bleed & extraction pipes descending cleanly through the black skid
    const bleedSpecs: Array<[number, number, number]> = [
      [0.36, 1.32, 0.088],
      [0.44, 1.72, 0.072],
      [0.46, 2.08, 0.076],
      [0.42, 2.48, 0.105], // Main 12-bar CHP controlled extraction pipe!
      [0.36, 2.94, 0.082],
    ];
    for (const [bx, bz, br] of bleedSpecs) {
      const pipe = new THREE.Mesh(new THREE.CylinderGeometry(br, br, 1.18, 20), this.matPipeSteel);
      pipe.position.set(bx, -0.92, bz);
      pipe.castShadow = true;
      g.add(pipe);

      const flange = this.makeBoltedFlange(br * 1.54, 0.044);
      flange.position.set(bx, -0.46, bz);
      g.add(flange);
    }
  }

  /** Slate-blue exhaust gusset ribs, 32-bolt condenser exit flange, aerodynamic struts, manhole & chrome vent pipe. */
  private buildExhaustDiffuserDetails(): void {
    // 1. 16 External triangular radial stiffener ribs around the slate-blue conical shoulder
    const ribShape = new THREE.Shape();
    ribShape.moveTo(0.0, 0.0);
    ribShape.lineTo(0.32, 0.0);
    ribShape.lineTo(0.32, 0.88);
    ribShape.lineTo(0.06, 0.88);
    ribShape.closePath();
    const ribGeo = new THREE.ExtrudeGeometry(ribShape, {
      depth: 0.028,
      bevelEnabled: true,
      bevelSize: 0.004,
      bevelThickness: 0.004,
      bevelSegments: 2,
    });
    ribGeo.translate(0, 0, -0.014);

    for (let i = 0; i < 16; i++) {
      let phi = (i / 16) * Math.PI * 2;
      if (phi > Math.PI) phi -= Math.PI * 2;
      const rib = new THREE.Mesh(ribGeo, this.matExhaustBlue);
      const holder = new THREE.Group();
      holder.rotation.z = phi;
      rib.rotation.x = Math.PI / 2;
      rib.position.set(1.00, 0, 4.26);
      holder.add(rib);

      if (phi < 0) this.exhaustLower.add(holder);
      else if (phi > Math.PI / 2) this.exhaustUpperBack.add(holder);
      else this.exhaustUpperFront.add(holder);
    }

    // 2. 32-bolt perimeter condenser flange studs seated on the z = 5.44 m exhaust exit rim (r = 1.270..1.315 m)
    const exitBoltGeo = new THREE.CylinderGeometry(0.018, 0.018, 0.075, 6);
    exitBoltGeo.rotateX(Math.PI / 2);
    for (let i = 0; i < 32; i++) {
      let phi = (i / 32) * Math.PI * 2;
      if (phi > Math.PI) phi -= Math.PI * 2;
      const b = new THREE.Mesh(exitBoltGeo, this.matChrome);
      b.position.set(1.292 * Math.cos(phi), 1.292 * Math.sin(phi), 5.44);
      if (phi < 0) this.exhaustLower.add(b);
      else if (phi > Math.PI / 2) this.exhaustUpperBack.add(b);
      else this.exhaustUpperFront.add(b);
    }

    // 3. Four internal aerodynamic streamlined struts supporting the rear bearing cone inside the exhaust hood
    const strutGeo = new THREE.BoxGeometry(0.068, 0.74, 0.44);
    const strutAngles = [-Math.PI * 0.5, -Math.PI * 0.14, Math.PI * 0.86, Math.PI * 0.52];
    for (const phi of strutAngles) {
      const strut = new THREE.Mesh(strutGeo, this.matExhaustBlue);
      strut.position.set(0.76 * Math.cos(phi), 0.76 * Math.sin(phi), 4.50);
      strut.rotation.z = phi - Math.PI / 2;
      if (phi < 0) this.exhaustLower.add(strut);
      else this.exhaustUpperBack.add(strut);
    }

    // 4. Lower-front circular bolted inspection manhole hatch (10 hex bolts + machined rim)
    const hatchGroup = new THREE.Group();
    const hatchAngle = -0.62;
    const hr = 1.24;
    hatchGroup.position.set(hr * Math.cos(hatchAngle), hr * Math.sin(hatchAngle), 4.78);
    hatchGroup.rotation.z = hatchAngle - Math.PI / 2;

    const neck = new THREE.Mesh(new THREE.CylinderGeometry(0.22, 0.24, 0.16, 28), this.matExhaustBlue);
    const lid = new THREE.Mesh(new THREE.CylinderGeometry(0.265, 0.265, 0.048, 28), this.matExhaustBlue);
    lid.position.y = 0.092;
    hatchGroup.add(neck, lid);
    const hBoltGeo = new THREE.CylinderGeometry(0.014, 0.014, 0.058, 6);
    for (let b = 0; b < 10; b++) {
      const a = (b / 10) * Math.PI * 2;
      const hb = new THREE.Mesh(hBoltGeo, this.matChrome);
      hb.position.set(0.215 * Math.cos(a), 0.112, 0.215 * Math.sin(a));
      hatchGroup.add(hb);
    }
    this.exhaustLower.add(hatchGroup);

    // 5. Top atmospheric relief hatch + tall thin vertical chrome vent pipe
    const ventBase = new THREE.Mesh(new THREE.CylinderGeometry(0.16, 0.18, 0.16, 22), this.matExhaustBlue);
    ventBase.position.set(-0.18, 1.32, 4.92);
    const ventPipe = new THREE.Mesh(new THREE.CylinderGeometry(0.022, 0.022, 0.88, 14), this.matChrome);
    ventPipe.position.set(-0.14, 1.76, 5.10);
    this.exhaustUpperBack.add(ventBase, ventPipe);
  }

  /**
   * Coupled 60 MVA 2-pole synchronous generator with cutaway outer enclosure revealing
   * laminated M19 stator core, gleaming braided copper end-windings, and forged 2-pole field rotor!
   */
  private buildCoupledGenerator(): void {
    const g = this.generatorGroup;

    // Flexible bolted shaft coupling disc + 12 pins + rotating 2-pole field rotor inside generator
    const coupDisc = new THREE.Mesh(new THREE.CylinderGeometry(0.24, 0.24, 0.10, 36), this.matChrome);
    coupDisc.rotation.x = Math.PI / 2;
    coupDisc.position.z = -0.07;
    this.generatorCouplingSpin.add(coupDisc);

    const pinGeo = new THREE.CylinderGeometry(0.016, 0.016, 0.14, 8);
    pinGeo.rotateX(Math.PI / 2);
    for (let i = 0; i < 12; i++) {
      const a = (i / 12) * Math.PI * 2;
      const pin = new THREE.Mesh(pinGeo, this.matSkidBlack);
      pin.position.set(0.185 * Math.cos(a), 0.185 * Math.sin(a), -0.07);
      this.generatorCouplingSpin.add(pin);
    }

    // Internal 2-pole forged field rotor + cooling fan blades turning inside generator cutaway!
    const fieldRotor = new THREE.Mesh(new THREE.CylinderGeometry(0.34, 0.34, 1.48, 32), this.matRotor);
    fieldRotor.rotation.x = Math.PI / 2;
    fieldRotor.position.z = -1.22;
    this.generatorCouplingSpin.add(fieldRotor);

    const retainRingGeo = new THREE.CylinderGeometry(0.365, 0.365, 0.18, 32);
    retainRingGeo.rotateX(Math.PI / 2);
    for (const rz of [-0.54, -1.90]) {
      const rr = new THREE.Mesh(retainRingGeo, this.matChrome);
      rr.position.z = rz;
      this.generatorCouplingSpin.add(rr);
    }
    g.add(this.generatorCouplingSpin);

    // Laminated M19 electrical steel stator core (open front-upper quarter in cutaway!)
    const coreGeo = new THREE.CylinderGeometry(0.68, 0.68, 1.24, 36, 1, false, Math.PI * 0.5, Math.PI * 1.5);
    coreGeo.rotateX(Math.PI / 2);
    const statorCore = new THREE.Mesh(coreGeo, this.matStatorCore);
    statorCore.position.z = -1.22;
    g.add(statorCore);

    // Gleaming braided copper stator end-windings protruding from both ends of M19 core!
    for (const cz of [-0.56, -1.88]) {
      const endRing = new THREE.Mesh(new THREE.TorusGeometry(0.49, 0.088, 14, 36), this.matCopperCoil);
      endRing.position.z = cz;
      g.add(endRing);
    }

    // Main outer stator housing: 270° solid wrap (always shown) + 90° removable front cutaway door!
    const outer270 = new THREE.CylinderGeometry(0.84, 0.84, 1.92, 36, 1, false, Math.PI * 0.5, Math.PI * 1.5);
    outer270.rotateX(Math.PI / 2);
    const barrelBack = new THREE.Mesh(outer270, this.matGenBlue);
    barrelBack.position.z = -1.22;
    barrelBack.castShadow = true;
    barrelBack.receiveShadow = true;
    g.add(barrelBack);

    const outer90 = new THREE.CylinderGeometry(0.84, 0.84, 1.92, 18, 1, false, 0, Math.PI * 0.5);
    outer90.rotateX(Math.PI / 2);
    const barrelFrontDoor = new THREE.Mesh(outer90, this.matGenBlue);
    barrelFrontDoor.position.z = -1.22;
    barrelFrontDoor.castShadow = true;
    this.generatorOuterUpperFront.add(barrelFrontDoor);
    g.add(this.generatorOuterUpperFront);

    // Cooling ribs around back 270°
    for (const fz of [-0.45, -0.78, -1.12, -1.46, -1.80, -2.06]) {
      const fin = new THREE.Mesh(
        new THREE.CylinderGeometry(0.89, 0.89, 0.04, 32, 1, false, Math.PI * 0.5, Math.PI * 1.5),
        this.matGenBlue,
      );
      fin.rotation.x = Math.PI / 2;
      fin.position.z = fz;
      g.add(fin);
    }

    // End shields
    const shieldGeo = new THREE.RingGeometry(0.24, 0.84, 32);
    for (const sz of [-0.26, -2.18]) {
      const sh = new THREE.Mesh(shieldGeo, this.matGenBlue);
      sh.position.z = sz;
      sh.material = this.matGenBlue;
      g.add(sh);
    }

    // Top air-to-water cooler box + 11 kV HV bushing terminal box
    const coolerBox = new THREE.Mesh(new THREE.BoxGeometry(1.18, 0.44, 1.46), this.matGenBlue);
    coolerBox.position.set(0, 0.96, -1.22);
    coolerBox.castShadow = true;
    g.add(coolerBox);

    const hvBox = new THREE.Mesh(new THREE.BoxGeometry(0.46, 0.34, 0.42), this.matOuterCasing);
    hvBox.position.set(0.76, 0.56, -0.82);
    g.add(hvBox);

    // Generator baseplate extension on matte-black skid
    const genSkid = new THREE.Mesh(new THREE.BoxGeometry(1.96, 0.46, 2.24), this.matSkidBlack);
    genSkid.position.set(0, -1.14, -1.24);
    genSkid.castShadow = true;
    genSkid.receiveShadow = true;
    g.add(genSkid);
  }

  /** Movable 3D studio measurement microphone capsule (aligned to fs-aeroac virtual mic coordinates!). */
  private buildMeasurementMicrophone(): void {
    const stand = new THREE.Mesh(new THREE.CylinderGeometry(0.012, 0.016, 1.38, 12), this.matSkidBlack);
    stand.position.y = 0.69;
    const base = new THREE.Mesh(new THREE.CylinderGeometry(0.14, 0.16, 0.024, 20), this.matSkidBlack);
    base.position.y = 0.012;
    const capsule = new THREE.Mesh(new THREE.CylinderGeometry(0.016, 0.022, 0.22, 16), this.matChrome);
    capsule.rotation.x = Math.PI / 2;
    capsule.position.set(0, 1.40, -0.06);
    const tipRing = new THREE.Mesh(
      new THREE.TorusGeometry(0.018, 0.004, 8, 18),
      new THREE.MeshBasicMaterial({ color: "#38bdf8" }),
    );
    tipRing.position.set(0, 1.40, -0.17);
    this.micAssembly.add(base, stand, capsule, tipRing);
    this.micAssembly.visible = false; // Shown cleanly when Acoustics tab, Audio, or Acoustic Vectors are active
    this.scene.add(this.micAssembly);
    this.setMicPosition(this.micPos.x, this.micPos.y, this.micPos.z);
  }

  public setMicVisible(visible: boolean): void {
    this.micAssembly.visible = visible;
  }

  public setMicPosition(standOffX: number, heightY: number, axialZ: number): void {
    this.micPos.set(standOffX, heightY, axialZ);
    // Convert turbine-local (x_a = standOffX, z_a = axialZ) -> world (-z_a + 2.72, 0, x_a)
    this.micAssembly.position.set(-axialZ + 2.72, 0.0, standOffX);
  }

  private buildMistParticles() {
    const count = 220;
    const pos = new Float32Array(count * 3);
    const z = new Float32Array(count);
    const r = new Float32Array(count);
    const phi = new Float32Array(count);
    for (let i = 0; i < count; i++) {
      z[i] = 3.92 + Math.random() * 1.45;
      r[i] = 0.42 + Math.random() * 0.72;
      phi[i] = Math.random() * Math.PI * 2;
    }
    const geo = new THREE.BufferGeometry();
    geo.setAttribute("position", new THREE.BufferAttribute(pos, 3));
    const mat = new THREE.PointsMaterial({
      color: "#b7e4ff",
      size: 0.14,
      map: this.textures.softGlowTex,
      transparent: true,
      opacity: 0.0,
      blending: THREE.AdditiveBlending,
      depthWrite: false,
    });
    return { points: new THREE.Points(geo, mat), data: { z, r, phi } };
  }

  private buildRubSparks() {
    const count = 110;
    const pos = new Float32Array(count * 3);
    const vel = new Float32Array(count * 3);
    for (let i = 0; i < count; i++) pos[i * 3 + 2] = 0.92;
    const geo = new THREE.BufferGeometry();
    geo.setAttribute("position", new THREE.BufferAttribute(pos, 3));
    const mat = new THREE.PointsMaterial({
      color: "#ff9b3d",
      size: 0.095,
      map: this.textures.softGlowTex,
      transparent: true,
      opacity: 0.0,
      blending: THREE.AdditiveBlending,
      depthWrite: false,
    });
    return { points: new THREE.Points(geo, mat), vel };
  }

  // ---------- 6 GPU-Instanced Tapered-Arrow 3D Physics Vector Fields ----------

  private buildVelocityVectorField(arrowGeo: THREE.BufferGeometry): THREE.InstancedMesh {
    const mat = new THREE.MeshStandardMaterial({
      color: "#ffffff",
      emissive: "#15445e",
      emissiveIntensity: 0.55,
      roughness: 0.25,
      metalness: 0.35,
    });
    // 17 stages × 14 front/cutaway azimuthal arrows = 238 crisp Euler velocity-triangle arrows!
    const inst = new THREE.InstancedMesh(arrowGeo, mat, 17 * 14);
    this.vfVelocityGroup.add(inst);
    return inst;
  }

  private buildHeatFluxVectorField(arrowGeo: THREE.BufferGeometry): THREE.InstancedMesh {
    const mat = new THREE.MeshStandardMaterial({
      color: "#ffffff",
      emissive: "#b9380c",
      emissiveIntensity: 0.60,
      roughness: 0.30,
      metalness: 0.25,
    });
    // 16 axial stations × 12 azimuthal Fourier q = -k ∇T arrows = 192 arrows!
    const inst = new THREE.InstancedMesh(arrowGeo, mat, 16 * 12);
    this.vfHeatGroup.add(inst);
    return inst;
  }

  private buildForcesVectorField(arrowGeo: THREE.BufferGeometry): THREE.InstancedMesh {
    const mat = new THREE.MeshStandardMaterial({
      color: "#ffffff",
      emissive: "#0d6b56",
      emissiveIntensity: 0.55,
      roughness: 0.26,
      metalness: 0.35,
    });
    // 17 stages × 10 aerodynamic lift/thrust arrows + 12 balance-piston counter-thrust arrows = 182!
    const inst = new THREE.InstancedMesh(arrowGeo, mat, 17 * 10 + 12);
    this.vfForcesGroup.add(inst);
    return inst;
  }

  private buildStressVectorField(arrowGeo: THREE.BufferGeometry): THREE.InstancedMesh {
    const mat = new THREE.MeshStandardMaterial({
      color: "#ffffff",
      emissive: "#b45309",
      emissiveIntensity: 0.60,
      roughness: 0.26,
      metalness: 0.35,
    });
    // 17 stages × 8 centrifugal radial arrows + 80 casing Lamé hoop arrows = 216 arrows!
    const inst = new THREE.InstancedMesh(arrowGeo, mat, 17 * 8 + 80);
    this.vfStressGroup.add(inst);
    return inst;
  }

  private buildWhirlOrbitOverlay(): {
    front: THREE.Mesh;
    rear: THREE.Mesh;
    trail: THREE.Line;
    wedgeFront: THREE.Mesh;
    wedgeRear: THREE.Mesh;
  } {
    // Journal bearing clearance ring (cyan) + 175 µm balance-piston labyrinth rub ring (crimson)
    const brgRingGeo = new THREE.TorusGeometry(0.28, 0.010, 12, 56);
    const brgRingMat = new THREE.MeshBasicMaterial({ color: "#38bdf8", transparent: true, opacity: 0.85 });
    const frontRing = new THREE.Mesh(brgRingGeo, brgRingMat);
    frontRing.position.z = 0.50;
    const rearRing = new THREE.Mesh(brgRingGeo, brgRingMat);
    rearRing.position.z = 4.44;

    const rubRing = new THREE.Mesh(
      new THREE.TorusGeometry(0.32, 0.009, 10, 56),
      new THREE.MeshBasicMaterial({ color: "#ff5266", transparent: true, opacity: 0.82 }),
    );
    rubRing.position.z = 1.05;

    // Hydrodynamic Ocvirk oil-film pressure crescent wedges inside both bearings
    const wedgeGeo = new THREE.RingGeometry(0.155, 0.255, 28, 1, -Math.PI * 0.85, Math.PI * 0.72);
    const wedgeMat = new THREE.MeshBasicMaterial({
      color: "#f59e0b",
      transparent: true,
      opacity: 0.45,
      side: THREE.DoubleSide,
    });
    const wedgeFront = new THREE.Mesh(wedgeGeo, wedgeMat);
    wedgeFront.position.z = 0.50;
    const wedgeRear = new THREE.Mesh(wedgeGeo, wedgeMat);
    wedgeRear.position.z = 4.44;

    const pinGeo = new THREE.SphereGeometry(0.032, 16, 16);
    const pinMat = new THREE.MeshBasicMaterial({ color: "#facc15" });
    const front = new THREE.Mesh(pinGeo, pinMat);
    front.position.z = 0.50;
    const rear = new THREE.Mesh(pinGeo, pinMat);
    rear.position.z = 4.44;

    // Live 3D 128-point whirl orbit ribbon at front bearing plane
    const trailPos = new Float32Array(128 * 3);
    const trailGeo = new THREE.BufferGeometry();
    trailGeo.setAttribute("position", new THREE.BufferAttribute(trailPos, 3));
    const trail = new THREE.Line(
      trailGeo,
      new THREE.LineBasicMaterial({ color: "#7dd3fc", transparent: true, opacity: 0.92 }),
    );

    this.vfWhirlGroup.add(frontRing, rearRing, rubRing, wedgeFront, wedgeRear, front, rear, trail);
    return { front, rear, trail, wedgeFront, wedgeRear };
  }

  private buildAcousticVectorField(arrowGeo: THREE.BufferGeometry): THREE.InstancedMesh {
    const mat = new THREE.MeshStandardMaterial({
      color: "#ffffff",
      emissive: "#4f46e5",
      emissiveIntensity: 0.65,
      roughness: 0.25,
      metalness: 0.3,
    });
    // 2 dipole/quadrupole acoustic sources × 72 spherical Poynting intensity arrows = 144 arrows!
    const inst = new THREE.InstancedMesh(arrowGeo, mat, 144);
    this.vfAcousticGroup.add(inst);
    return inst;
  }

  // ---------- Interactive 3D CAD Surface Probe ----------

  public setProbeCallback(cb: ((reading: SurfaceProbeReading | null) => void) | null): void {
    this.onProbeCallback = cb;
  }

  public setFirstInteractCallback(cb: (() => void) | null): void {
    this.onFirstInteractCallback = cb;
  }

  public setNavMode(mode: NavMode): void {
    this.navMode = mode;
    if (mode === "pan") {
      this.controls.mouseButtons.LEFT = THREE.MOUSE.PAN;
      this.controls.touches.ONE = THREE.TOUCH.PAN;
    } else {
      this.controls.mouseButtons.LEFT = THREE.MOUSE.ROTATE;
      this.controls.touches.ONE = THREE.TOUCH.ROTATE;
    }
    this.controls.touches.TWO = THREE.TOUCH.DOLLY_PAN;
  }

  public getNavMode(): NavMode {
    return this.navMode;
  }

  /** Hold continuous smooth zoom while pressing +/- nav buttons on mouse or touch. */
  public setZoomHold(dir: -1 | 0 | 1): void {
    this.zoomHoldDir = dir;
    if (dir !== 0) this.camLerp.active = false;
  }

  /** Smoothly dolly the camera closer (< 1) or farther (> 1) around the current target. */
  public zoomByFactor(factor: number): void {
    const offset = this.camera.position.clone().sub(this.controls.target);
    const curDist = offset.length();
    const nextDist = Math.max(this.controls.minDistance, Math.min(this.controls.maxDistance, curDist * factor));
    offset.setLength(nextDist);
    this.camLerp.goal.target.copy(this.controls.target);
    this.camLerp.goal.pos.copy(this.controls.target).add(offset);
    this.camLerp.active = true;
  }

  private isHierarchyVisible(obj: THREE.Object3D): boolean {
    let cur: THREE.Object3D | null = obj;
    while (cur) {
      if (!cur.visible) return false;
      cur = cur.parent;
    }
    return true;
  }

  private allProbeTargets(): THREE.Object3D[] {
    return this.bladeProbeTargets.length > 0
      ? [...this.probeTargets, ...this.bladeProbeTargets]
      : this.probeTargets;
  }

  /** Smoothly focus & frame the exact 3D CAD surface under a double-click / double-tap, or reset to Overview. */
  private focusAtClientPoint(clientX: number, clientY: number): void {
    const rect = this.renderer.domElement.getBoundingClientRect();
    this.pointerNdc.x = ((clientX - rect.left) / rect.width) * 2 - 1;
    this.pointerNdc.y = -((clientY - rect.top) / rect.height) * 2 + 1;
    this.raycaster.setFromCamera(this.pointerNdc, this.camera);

    const hits = this.raycaster.intersectObjects(this.allProbeTargets(), false);
    const hit = hits.find((h) => this.isHierarchyVisible(h.object));
    if (!hit) {
      this.flyToPreset("overview");
      return;
    }

    const dir = this.camera.position.clone().sub(hit.point);
    const curDist = dir.length();
    const targetDist = Math.max(1.35, Math.min(3.4, curDist * 0.56));
    dir.normalize().multiplyScalar(targetDist);
    this.camLerp.goal.target.copy(hit.point);
    this.camLerp.goal.pos.copy(hit.point).add(dir);
    this.camLerp.active = true;
  }

  private bindPointerProbe(): void {
    const el = this.renderer.domElement;

    // Track continuous WASD / Arrow camera panning (Shift+Arrow orbits!) & Q/E zoom
    this.boundKeyDown = (ev: KeyboardEvent) => {
      this.shiftHeld = ev.shiftKey;
      if (ev.target instanceof HTMLInputElement || ev.target instanceof HTMLTextAreaElement) return;
      if (["ArrowUp", "ArrowDown", "ArrowLeft", "ArrowRight", "KeyW", "KeyA", "KeyS", "KeyD", "KeyQ", "KeyE"].includes(ev.code)) {
        this.heldKeys.add(ev.code);
        this.camLerp.active = false;
      }
    };
    this.boundKeyUp = (ev: KeyboardEvent) => {
      this.shiftHeld = ev.shiftKey;
      this.heldKeys.delete(ev.code);
    };
    this.boundBlur = () => {
      this.heldKeys.clear();
      this.shiftHeld = false;
      this.zoomHoldDir = 0;
    };
    this.boundPointerUp = () => {
      this.isPointerDown = false;
      this.zoomHoldDir = 0;
    };
    window.addEventListener("keydown", this.boundKeyDown);
    window.addEventListener("keyup", this.boundKeyUp);
    window.addEventListener("blur", this.boundBlur);
    window.addEventListener("pointerup", this.boundPointerUp);

    el.addEventListener("pointerdown", (ev: PointerEvent) => {
      this.isPointerDown = true;
      // Touch double-tap detection for mobile CAD focus!
      if (ev.pointerType === "touch" && ev.isPrimary) {
        const now = performance.now();
        const dist = Math.hypot(ev.clientX - this.lastTapX, ev.clientY - this.lastTapY);
        if (now - this.lastTapMs < 320 && dist < 28) {
          this.focusAtClientPoint(ev.clientX, ev.clientY);
          this.lastTapMs = 0;
        } else {
          this.lastTapMs = now;
          this.lastTapX = ev.clientX;
          this.lastTapY = ev.clientY;
        }
      }
    });
    el.addEventListener("dblclick", (ev: MouseEvent) => {
      this.focusAtClientPoint(ev.clientX, ev.clientY);
    });
    el.addEventListener("pointerleave", () => {
      this.probeMarker.visible = false;
      this.onProbeCallback?.(null);
    });
    el.addEventListener("pointermove", (ev: PointerEvent) => {
      // Never run hover raycast on touch pointers so smartphone drag/pinch stays 100% zero-overhead!
      if (ev.pointerType === "touch" || this.isPointerDown || !this.onProbeCallback) return;
      const now = performance.now();
      if (now - this.lastProbeMs < 42) return;
      this.lastProbeMs = now;

      const rect = el.getBoundingClientRect();
      this.pointerNdc.x = ((ev.clientX - rect.left) / rect.width) * 2 - 1;
      this.pointerNdc.y = -((ev.clientY - rect.top) / rect.height) * 2 + 1;
      this.raycaster.setFromCamera(this.pointerNdc, this.camera);

      const hits = this.raycaster.intersectObjects(this.allProbeTargets(), false);
      const hit = hits.find((h) => this.isHierarchyVisible(h.object));
      if (!hit) {
        this.probeMarker.visible = false;
        this.onProbeCallback(null);
        return;
      }

      this.probeMarker.visible = true;
      this.probeMarker.position.copy(hit.point);
      if (hit.face) {
        const nWorld = hit.face.normal.clone().transformDirection(hit.object.matrixWorld).normalize();
        this.probeMarker.position.addScaledVector(nWorld, 0.012);
        this.probeMarker.quaternion.setFromUnitVectors(new THREE.Vector3(0, 0, 1), nWorld);
      }

      const localP = hit.point.clone().applyMatrix4(this.continuumUniforms.uTurbineInvMat.value);
      const zM = localP.z;
      const rM = Math.hypot(localP.x, localP.y);
      // probe_continuum returns [temp_c, sigma_vm_mpa, p_bar, mach, q_mag_kw_m2, wetness_pct, creep_log, region_id]
      const raw = this.sim.probe_continuum(zM, rM);
      const tempC = raw[0] ?? 25;
      const vmMpa = raw[1] ?? 0;
      const pBar = raw[2] ?? 0;
      const mach = raw[3] ?? 0;
      const qMagKwM2 = raw[4] ?? 0;
      const wetnessPct = raw[5] ?? 0;
      const creepLog10 = raw[6] ?? -16;

      this.onProbeCallback({
        partName: classifyPartAtZR(zM, rM),
        zM,
        rM,
        tempC,
        vmMpa,
        pBar,
        mach,
        qMagKwM2,
        wetnessPct,
        creepLog10,
      });
    });
  }

  // ---------- Public Visual Controls ----------

  public setCutawayMode(mode: CutawayMode): void {
    this.cutawayMode = mode;
    const showUpperFront = mode === "closed";
    const showUpperBack = mode !== "open";

    this.innerCarrierUpperFront.visible = showUpperFront;
    this.outerCasingUpperFront.visible = showUpperFront;
    this.exhaustUpperFront.visible = showUpperFront;
    this.frontBrgUpperFront.visible = showUpperFront;
    this.rearBrgUpperFront.visible = showUpperFront;
    this.generatorOuterUpperFront.visible = showUpperFront;

    this.innerCarrierUpperBack.visible = showUpperBack;
    this.outerCasingUpperBack.visible = showUpperBack;
    this.exhaustUpperBack.visible = showUpperBack;
    this.frontBrgUpperBack.visible = showUpperBack;
    this.rearBrgUpperBack.visible = showUpperBack;
  }

  public setShaderMode(mode: ShaderMode): void {
    this.shaderMode = mode;
    this.continuumUniforms.uFieldMode.value = FIELD_MODE_SPECS[mode].modeIndex;
  }

  public setIsoContours(enabled: boolean): void {
    this.isoContoursEnabled = enabled;
    this.continuumUniforms.uGridBounds.value.w = enabled ? 1.0 : 0.0;
  }

  public setExplodeFraction(frac: number): void {
    this.explodeFrac = Math.max(0, Math.min(1, frac));
    const e = this.explodeFrac;
    this.outerCasingUpperBack.position.set(-0.35 * e, 1.05 * e, 0);
    this.outerCasingUpperFront.position.set(0.45 * e, 1.05 * e, 0);
    this.innerCarrierUpperBack.position.set(-0.18 * e, 0.52 * e, 0);
    this.innerCarrierUpperFront.position.set(0.24 * e, 0.52 * e, 0);
    this.exhaustUpperBack.position.set(-0.22 * e, 0.85 * e, 0.55 * e);
    this.exhaustUpperFront.position.set(0.28 * e, 0.85 * e, 0.55 * e);
    this.exhaustLower.position.set(0, -0.18 * e, 0.55 * e);
    this.frontBrgUpperBack.position.set(0, 0.62 * e, -0.35 * e);
    this.frontBrgUpperFront.position.set(0, 0.62 * e, -0.35 * e);
    this.generatorGroup.position.set(0, 0, -0.65 * e);
  }

  public setVisualSlowMo(frac: number): void {
    this.visualSlowMo = Math.max(0.01, Math.min(1.0, frac));
  }

  public setShowGenerator(show: boolean): void {
    this.generatorGroup.visible = show;
  }

  public setShowSteam(show: boolean): void {
    this.showSteamMist = show;
    this.mistPoints.visible = show;
  }

  public setExaggerateOrbit(exag: boolean): void {
    this.exaggerateOrbit = exag;
  }

  public setVectorFields(flags: VectorFieldFlags): void {
    this.vfVelocityGroup.visible = flags.velocity;
    this.vfHeatGroup.visible = flags.heat;
    this.vfForcesGroup.visible = flags.forces;
    this.vfStressGroup.visible = flags.stress;
    this.vfWhirlGroup.visible = flags.whirl;
    this.vfAcousticGroup.visible = flags.acoustic;
  }

  private computePresetGoal(preset: CameraPreset): { pos: THREE.Vector3; target: THREE.Vector3 } {
    const spec = CAMERA_PRESETS[preset];
    const target = new THREE.Vector3(...spec.target);
    const pos = new THREE.Vector3(...spec.pos);
    // Portrait smartphones & narrow tablets automatically pull back cleanly so the entire turbine fits!
    const aspect = this.camera.aspect || 1.6;
    if (aspect < 1.05) {
      const pull = Math.min(2.05, 1.08 / Math.max(0.46, aspect));
      pos.sub(target).multiplyScalar(pull).add(target);
    }
    return { pos, target };
  }

  private applyPresetImmediate(preset: CameraPreset): void {
    const goal = this.computePresetGoal(preset);
    this.camera.position.copy(goal.pos);
    if (this.controls) {
      this.controls.target.copy(goal.target);
    }
  }

  public flyToPreset(preset: CameraPreset): void {
    const goal = this.computePresetGoal(preset);
    this.camLerp.goal.pos.copy(goal.pos);
    this.camLerp.goal.target.copy(goal.target);
    this.camLerp.active = true;
  }

  // ---------- Per-Frame Live Multiphysics GPU Sync & HDR Render ----------

  public updateAndRender(
    snap: TurbineSnapshot,
    stages: StageTelemetry[],
    continuumFlat: Float32Array,
    orbitFlat: Float32Array,
    dtWallS: number,
    simRunning: boolean,
    bearingClearanceUm = 220,
  ): void {
    const dt = Math.min(0.05, Math.max(0.001, dtWallS));

    if (this.camLerp.active) {
      const k = 1 - Math.exp(-6.5 * dt);
      this.camera.position.lerp(this.camLerp.goal.pos, k);
      this.controls.target.lerp(this.camLerp.goal.target, k);
      if (
        this.camera.position.distanceTo(this.camLerp.goal.pos) < 0.012 &&
        this.controls.target.distanceTo(this.camLerp.goal.target) < 0.012
      ) {
        this.camLerp.active = false;
      }
    }
    // Smooth continuous press-and-hold +/- zoom on desktop mouse & smartphone touch
    if (this.zoomHoldDir !== 0) {
      const off = this.camera.position.clone().sub(this.controls.target);
      const f = this.zoomHoldDir > 0 ? 1 - 1.15 * dt : 1 + 1.15 * dt;
      const nextD = Math.max(this.controls.minDistance, Math.min(this.controls.maxDistance, off.length() * f));
      off.setLength(nextD);
      this.camera.position.copy(this.controls.target).add(off);
    }
    // Smooth continuous WASD / Arrow panning, Shift+Arrow orbiting & Q/E zoom
    if (this.heldKeys.size > 0) {
      if (
        this.shiftHeld &&
        (this.heldKeys.has("ArrowLeft") ||
          this.heldKeys.has("ArrowRight") ||
          this.heldKeys.has("ArrowUp") ||
          this.heldKeys.has("ArrowDown"))
      ) {
        const off = this.camera.position.clone().sub(this.controls.target);
        const sph = new THREE.Spherical().setFromVector3(off);
        const rotStep = 1.15 * dt;
        if (this.heldKeys.has("ArrowLeft")) sph.theta += rotStep;
        if (this.heldKeys.has("ArrowRight")) sph.theta -= rotStep;
        if (this.heldKeys.has("ArrowUp")) sph.phi = Math.max(this.controls.minPolarAngle, sph.phi - rotStep);
        if (this.heldKeys.has("ArrowDown")) sph.phi = Math.min(this.controls.maxPolarAngle, sph.phi + rotStep);
        off.setFromSpherical(sph);
        this.camera.position.copy(this.controls.target).add(off);
      } else {
        const dist = this.camera.position.distanceTo(this.controls.target);
        const panStep = Math.max(0.9, dist * 0.48) * dt;
        const right = new THREE.Vector3().setFromMatrixColumn(this.camera.matrixWorld, 0).normalize();
        const up = new THREE.Vector3().setFromMatrixColumn(this.camera.matrixWorld, 1).normalize();
        const delta = new THREE.Vector3();
        if (this.heldKeys.has("KeyA") || this.heldKeys.has("ArrowLeft")) delta.addScaledVector(right, -panStep);
        if (this.heldKeys.has("KeyD") || this.heldKeys.has("ArrowRight")) delta.addScaledVector(right, panStep);
        if (this.heldKeys.has("KeyW") || this.heldKeys.has("ArrowUp")) delta.addScaledVector(up, panStep);
        if (this.heldKeys.has("KeyS") || this.heldKeys.has("ArrowDown")) delta.addScaledVector(up, -panStep);
        if (delta.lengthSq() > 0) {
          this.camera.position.add(delta);
          this.controls.target.add(delta);
        }
      }
      if (this.heldKeys.has("KeyE") || this.heldKeys.has("KeyQ")) {
        const off = this.camera.position.clone().sub(this.controls.target);
        const f = this.heldKeys.has("KeyE") ? 1 - 0.9 * dt : 1 + 0.9 * dt;
        const nextD = Math.max(this.controls.minDistance, Math.min(this.controls.maxDistance, off.length() * f));
        off.setLength(nextD);
        this.camera.position.copy(this.controls.target).add(off);
      }
    }

    this.controls.update();
    // Gently keep target within turbine hall bounds so user can never pan off into empty void
    this.controls.target.x = Math.max(-5.5, Math.min(5.5, this.controls.target.x));
    this.controls.target.y = Math.max(-0.4, Math.min(4.2, this.controls.target.y));
    this.controls.target.z = Math.max(-4.5, Math.min(4.5, this.controls.target.z));

    // 1. Rotor rotation & exaggerated 4-DOF lateral whirl orbit displacement
    if (simRunning) {
      this.visualTheta = (this.visualTheta + snap.omegaRadS * this.visualSlowMo * dt) % (Math.PI * 2);
    }
    this.rotorSpinGroup.rotation.z = this.visualTheta;
    this.generatorCouplingSpin.rotation.z = this.visualTheta;

    const orbitScale = this.exaggerateOrbit ? 250e-6 : 1e-6;
    const ox = Math.max(-0.08, Math.min(0.08, snap.rotorCenterXUm * orbitScale));
    const oy = Math.max(-0.08, Math.min(0.08, snap.rotorCenterYUm * orbitScale));
    this.rotorAssembly.position.set(ox, oy, 0);

    // 2. Live electro-hydraulic valve actuator stem strokes
    this.esvRodMesh.position.z = 0.40 + 0.14 * snap.esvOpen;
    this.cvRod1Mesh.position.x = -0.64 + 0.16 * snap.cvLift;
    this.cvRod2Mesh.position.x = -0.58 + 0.16 * snap.cvLift;
    const extStroke = 0.14 * snap.extLift;
    this.extRodMesh.position.y = 0.84 + extStroke * Math.cos(Math.PI * 0.34);
    this.extRodMesh.position.z = 2.60 + extStroke * Math.sin(Math.PI * 0.34);

    // 3. Upload 64×24×8 WASM 2D axisymmetric PDE tensor to GPU FloatType DataTextures!
    this.uploadContinuumTextures(continuumFlat);

    // 4. Subtle physical LP wetness condensation mist & labyrinth rub sparks
    if (this.showSteamMist) {
      this.updateWetMist(snap, dt, simRunning);
    }
    this.updateRubSparks(snap, dt, simRunning);

    // 5. Update active 3D Physics Vector Field overlays directly from 17-stage meanline & 2D PDE!
    this.updateVectorFields(snap, stages, continuumFlat, orbitFlat, bearingClearanceUm);

    this.composer.render();
  }

  private uploadContinuumTextures(flat: Float32Array): void {
    const nCells = CONTINUUM_NZ * CONTINUUM_NR;
    if (flat.length < nCells * CONTINUUM_STRIDE) return;
    for (let i = 0; i < nCells; i++) {
      const s = i * CONTINUUM_STRIDE;
      const d = i * 4;
      this.fieldAData[d] = flat[s] ?? 300; // T_K
      this.fieldAData[d + 1] = flat[s + 1] ?? 0; // σ_vM MPa
      this.fieldAData[d + 2] = flat[s + 2] ?? 0; // P_bar
      this.fieldAData[d + 3] = flat[s + 3] ?? 0; // Mach

      this.fieldBData[d] = flat[s + 4] ?? 0; // q_r kW/m²
      this.fieldBData[d + 1] = flat[s + 5] ?? 0; // q_z kW/m²
      this.fieldBData[d + 2] = flat[s + 6] ?? 0; // wetness (1-x)
      this.fieldBData[d + 3] = flat[s + 7] ?? -16; // log10 creep
    }
    this.fieldATex.needsUpdate = true;
    this.fieldBTex.needsUpdate = true;
  }

  private updateWetMist(snap: TurbineSnapshot, dt: number, running: boolean): void {
    const moisture = Math.max(0, 1 - snap.lpQualityX);
    const mistMat = this.mistPoints.material as THREE.PointsMaterial;
    // Only faintly visible when LP exhaust truly drops below saturation (x < 0.96)
    mistMat.opacity = moisture < 0.04 ? 0.0 : Math.min(0.42, (moisture - 0.03) * 2.6);
    if (mistMat.opacity <= 0.005) return;

    const flowNorm = Math.max(0.05, Math.min(2.0, snap.lpFlowKgS / 35.0));
    const spinNorm = Math.max(0.05, snap.rpm / 3000.0);
    const mPos = this.mistPoints.geometry.attributes.position as THREE.BufferAttribute;
    for (let i = 0; i < this.mistData.z.length; i++) {
      let z = this.mistData.z[i] ?? 4.1;
      let phi = this.mistData.phi[i] ?? 0;
      if (running) {
        z += (1.1 + 2.2 * flowNorm) * dt;
        phi += 0.6 * spinNorm * dt;
        if (z > 5.42) z = 3.92;
        this.mistData.z[i] = z;
        this.mistData.phi[i] = phi;
      }
      const r = this.mistData.r[i] ?? 0.65;
      mPos.setXYZ(i, r * Math.cos(phi), r * Math.sin(phi), z);
    }
    mPos.needsUpdate = true;
  }

  private updateRubSparks(snap: TurbineSnapshot, dt: number, running: boolean): void {
    const active = snap.rubForceKn > 0.05;
    const mat = this.sparkPoints.material as THREE.PointsMaterial;
    mat.opacity = active ? Math.min(0.98, 0.45 + snap.rubForceKn * 0.06) : 0.0;
    if (!active || !running) return;

    const pos = this.sparkPoints.geometry.attributes.position as THREE.BufferAttribute;
    for (let i = 0; i < pos.count; i++) {
      let x = pos.getX(i);
      let y = pos.getY(i);
      let z = pos.getZ(i);
      const r = Math.hypot(x, y);
      if (r < 0.22 || r > 0.78 || Math.random() < 0.08) {
        const a = Math.random() * Math.PI * 2;
        x = 0.30 * Math.cos(a);
        y = 0.30 * Math.sin(a);
        z = 0.84 + Math.random() * 0.36;
        const tSpeed = 2.4 + Math.random() * 3.0;
        this.sparkVel[i * 3] = -Math.sin(a) * tSpeed + (Math.random() - 0.5) * 1.2;
        this.sparkVel[i * 3 + 1] = Math.cos(a) * tSpeed + (Math.random() - 0.5) * 1.2;
        this.sparkVel[i * 3 + 2] = (Math.random() - 0.5) * 1.4;
      } else {
        x += (this.sparkVel[i * 3] ?? 0) * dt;
        y += (this.sparkVel[i * 3 + 1] ?? 0) * dt;
        z += (this.sparkVel[i * 3 + 2] ?? 0) * dt;
      }
      pos.setXYZ(i, x, y, z);
    }
    pos.needsUpdate = true;
  }

  private updateVectorFields(
    snap: TurbineSnapshot,
    stages: StageTelemetry[],
    continuumFlat: Float32Array,
    orbitFlat: Float32Array,
    bearingClearanceUm = 220,
  ): void {
    const dir = this.tmpDir;

    // 1. 17-Stage Euler Velocity Triangle Field (c_axial, c_swirl) directly from stage_meanline.rs!
    if (this.vfVelocityGroup.visible && stages.length >= 17) {
      let idx = 0;
      for (let s = 0; s < 17; s++) {
        const st = stages[s]!;
        const rMid = 0.5 * (st.rHubM + st.rTipM) + 0.05;
        const ca = Math.max(8, st.cAxialMS);
        const cu = 0.5 * (st.cSwirlInMS + st.cSwirlOutMS);
        const vMag = Math.hypot(ca, cu);
        const len = Math.max(0.11, Math.min(0.42, (vMag / 480) * 0.38));
        const uHeat = Math.max(0, Math.min(1, (st.tOutC - 40) / 480));

        for (let k = 0; k < 14; k++) {
          if (idx >= this.vfVelocityInst.count) break;
          // Distribute along front/upper cutaway arc so vectors stay crisply legible
          const phi = -0.55 + (k / 13) * 2.55;
          const cp = Math.cos(phi);
          const sp = Math.sin(phi);
          this.tmpPos.set(rMid * cp, rMid * sp, st.zM);
          dir.set(-cu * sp, cu * cp, ca).normalize();
          this.tmpQuat.setFromUnitVectors(this.upVec, dir);
          this.tmpScale.set(0.76, len, 0.76);
          this.tmpMat.compose(this.tmpPos, this.tmpQuat, this.tmpScale);
          this.vfVelocityInst.setMatrixAt(idx, this.tmpMat);
          // Superheated orange-gold -> expanding turquoise -> wet cyan!
          this.tmpColor.setHSL(0.54 - 0.48 * uHeat, 0.92, 0.58);
          this.vfVelocityInst.setColorAt(idx, this.tmpColor);
          idx++;
        }
      }
      this.vfVelocityInst.instanceMatrix.needsUpdate = true;
      if (this.vfVelocityInst.instanceColor) this.vfVelocityInst.instanceColor.needsUpdate = true;
    }

    // 2. Exact 2D Axisymmetric Fourier Heat Flux Field q = -k(T) ∇T from continuum_2d.rs!
    if (this.vfHeatGroup.visible && continuumFlat.length >= CONTINUUM_NZ * CONTINUUM_NR * CONTINUUM_STRIDE) {
      let idx = 0;
      for (let s = 0; s < 16; s++) {
        const iz = Math.min(CONTINUUM_NZ - 1, Math.round(4 + s * 3.6));
        const z = this.gridZMin + ((iz + 0.5) / CONTINUUM_NZ) * this.gridZSpan;
        const rSample = z > 3.68 ? 1.18 : 0.72;
        const ir = Math.min(CONTINUUM_NR - 1, Math.max(0, Math.floor((rSample / this.gridRMax) * CONTINUUM_NR)));
        const off = (ir * CONTINUUM_NZ + iz) * CONTINUUM_STRIDE;
        const qr = continuumFlat[off + 4] ?? 0;
        const qz = continuumFlat[off + 5] ?? 0;
        const qMag = Math.hypot(qr, qz);
        const len = Math.max(0.10, Math.min(0.42, 0.12 + Math.log1p(qMag * 0.25) * 0.14));
        const hot = Math.max(0, Math.min(1, ((continuumFlat[off] ?? 300) - 310) / 480));

        for (let k = 0; k < 12; k++) {
          if (idx >= this.vfHeatInst.count) break;
          const phi = (k / 12) * Math.PI * 2;
          const cp = Math.cos(phi);
          const sp = Math.sin(phi);
          this.tmpPos.set((rSample + 0.06) * cp, (rSample + 0.06) * sp, z);
          // Map axisymmetric (q_r, q_z) into 3D vector!
          const qrVis = Math.abs(qr) < 0.05 ? 0.5 : qr;
          dir.set(qrVis * cp, qrVis * sp, qz).normalize();
          this.tmpQuat.setFromUnitVectors(this.upVec, dir);
          this.tmpScale.set(0.85, len, 0.85);
          this.tmpMat.compose(this.tmpPos, this.tmpQuat, this.tmpScale);
          this.vfHeatInst.setMatrixAt(idx, this.tmpMat);
          this.tmpColor.setHSL(0.08 - 0.07 * hot, 0.95, 0.56);
          this.vfHeatInst.setColorAt(idx, this.tmpColor);
          idx++;
        }
      }
      this.vfHeatInst.instanceMatrix.needsUpdate = true;
      if (this.vfHeatInst.instanceColor) this.vfHeatInst.instanceColor.needsUpdate = true;
    }

    // 3. Aerodynamic Blade Lift/Torque & Stage Axial Thrust + Balance-Piston Counter-Thrust!
    if (this.vfForcesGroup.visible && stages.length >= 17) {
      let idx = 0;
      for (let s = 0; s < 17; s++) {
        const st = stages[s]!;
        const rMid = 0.5 * (st.rHubM + st.rTipM) + 0.04;
        const fTang = Math.max(0.2, st.powerMw * 1.8);
        const fAx = Math.max(0.1, st.axialThrustKn * 0.08);
        const len = Math.max(0.12, Math.min(0.42, 0.12 + st.powerMw * 0.055));
        for (let k = 0; k < 10; k++) {
          if (idx >= this.vfForcesInst.count - 12) break;
          const phi = -0.3 + (k / 9) * 2.1;
          const cp = Math.cos(phi);
          const sp = Math.sin(phi);
          this.tmpPos.set(rMid * cp, rMid * sp, st.zM);
          // Tangential torque force + downstream axial thrust
          dir.set(-fTang * sp, fTang * cp, fAx).normalize();
          this.tmpQuat.setFromUnitVectors(this.upVec, dir);
          this.tmpScale.set(0.84, len, 0.84);
          this.tmpMat.compose(this.tmpPos, this.tmpQuat, this.tmpScale);
          this.vfForcesInst.setMatrixAt(idx, this.tmpMat);
          this.tmpColor.set("#34d399");
          this.vfForcesInst.setColorAt(idx, this.tmpColor);
          idx++;
        }
      }
      // 12 bright cyan opposing counter-thrust vectors on twin HP balance pistons (z = 0.96 m, pointing -z!)
      const bpLen = Math.max(0.16, Math.min(0.46, 0.16 + Math.abs(snap.axialThrustKn) * 0.0025));
      for (let k = 0; k < 12 && idx < this.vfForcesInst.count; k++) {
        const phi = (k / 12) * Math.PI * 2;
        this.tmpPos.set(0.36 * Math.cos(phi), 0.36 * Math.sin(phi), 0.96);
        dir.set(0, 0, -1); // Upstream counter-balance thrust!
        this.tmpQuat.setFromUnitVectors(this.upVec, dir);
        this.tmpScale.set(0.95, bpLen, 0.95);
        this.tmpMat.compose(this.tmpPos, this.tmpQuat, this.tmpScale);
        this.vfForcesInst.setMatrixAt(idx, this.tmpMat);
        this.tmpColor.set("#38bdf8");
        this.vfForcesInst.setColorAt(idx, this.tmpColor);
        idx++;
      }
      this.vfForcesInst.instanceMatrix.needsUpdate = true;
      if (this.vfForcesInst.instanceColor) this.vfForcesInst.instanceColor.needsUpdate = true;
    }

    // 4. Principal Centrifugal (Radial on Rotating Blades) + Lamé Hoop Stress (Tangential on HP Casing)!
    if (this.vfStressGroup.visible && stages.length >= 17) {
      let idx = 0;
      // Radial centrifugal tension arrows on all 17 stages
      for (let s = 0; s < 17; s++) {
        const st = stages[s]!;
        const rArrow = st.rTipM + 0.02;
        const len = Math.max(0.08, Math.min(0.45, (st.rootVmMpa / 520) * 0.42));
        const sev = Math.max(0, Math.min(1, st.rootVmMpa / 480));
        for (let k = 0; k < 8; k++) {
          if (idx >= 17 * 8) break;
          const phi = -0.3 + (k / 7) * 2.2;
          const cp = Math.cos(phi);
          const sp = Math.sin(phi);
          this.tmpPos.set(rArrow * cp, rArrow * sp, st.zM);
          dir.set(cp, sp, 0).normalize();
          this.tmpQuat.setFromUnitVectors(this.upVec, dir);
          this.tmpScale.set(0.80, len, 0.80);
          this.tmpMat.compose(this.tmpPos, this.tmpQuat, this.tmpScale);
          this.vfStressInst.setMatrixAt(idx, this.tmpMat);
          this.tmpColor.setHSL(0.14 - 0.14 * sev, 0.96, 0.56);
          this.vfStressInst.setColorAt(idx, this.tmpColor);
          idx++;
        }
      }
      // Tangential Lamé pressure + thermal hoop stress around HP/IP outer casing (80 arrows)
      const hoopScale = Math.max(0.14, Math.min(0.42, (snap.hpThermalStressMpa + snap.chestBar * 0.85) / 480));
      for (let ring = 0; ring < 8; ring++) {
        const z = 1.15 + ring * 0.34;
        const r = z < 2.4 ? 0.78 : 0.88;
        for (let k = 0; k < 10 && idx < this.vfStressInst.count; k++) {
          const phi = (k / 10) * Math.PI * 2;
          this.tmpPos.set(r * Math.cos(phi), r * Math.sin(phi), z);
          dir.set(-Math.sin(phi), Math.cos(phi), 0).normalize();
          this.tmpQuat.setFromUnitVectors(this.upVec, dir);
          this.tmpScale.set(0.82, hoopScale, 0.82);
          this.tmpMat.compose(this.tmpPos, this.tmpQuat, this.tmpScale);
          this.vfStressInst.setMatrixAt(idx, this.tmpMat);
          this.tmpColor.set("#f43f5e");
          this.vfStressInst.setColorAt(idx, this.tmpColor);
          idx++;
        }
      }
      this.vfStressInst.instanceMatrix.needsUpdate = true;
      if (this.vfStressInst.instanceColor) this.vfStressInst.instanceColor.needsUpdate = true;
    }

    // 5. 3D Bearing Oil-Film Pressure Crescent + Live 3D Whirl Orbit Ribbon
    if (this.vfWhirlGroup.visible) {
      const scale = 0.28 / Math.max(100, bearingClearanceUm); // Live bearing clearance -> 0.28 m display ring
      const px = Math.max(-0.34, Math.min(0.34, snap.rotorCenterXUm * scale));
      const py = Math.max(-0.34, Math.min(0.34, snap.rotorCenterYUm * scale));
      this.vfWhirlPinFront.position.set(px, py, 0.50);
      this.vfWhirlPinRear.position.set(px * 1.08, py * 1.08, 4.44);

      // Rotate oil-film pressure crescent opposite journal eccentricity attitude angle!
      const att = Math.atan2(py, px);
      this.vfOilWedgeFront.rotation.z = att;
      this.vfOilWedgeRear.rotation.z = att + 0.15;

      const trailAttr = this.vfWhirlTrailLine.geometry.attributes.position as THREE.BufferAttribute;
      const nPts = Math.min(128, Math.floor(orbitFlat.length / 2));
      for (let i = 0; i < 128; i++) {
        const k = Math.min(i, Math.max(0, nPts - 1));
        const ox = Math.max(-0.34, Math.min(0.34, (orbitFlat[k * 2] ?? 0) * scale));
        const oy = Math.max(-0.34, Math.min(0.34, (orbitFlat[k * 2 + 1] ?? 0) * scale));
        trailAttr.setXYZ(i, ox, oy, 0.50);
      }
      trailAttr.needsUpdate = true;
    }

    // 6. Farassat/Curle 3D Acoustic Poynting Intensity Vectors radiating from Valve & Exhaust Dipoles!
    if (this.vfAcousticGroup.visible) {
      const splScale = Math.max(0.14, Math.min(0.45, ((snap.splDb - 52) / 55) * 0.36));
      const pulse = 1.0 + 0.14 * Math.sin(performance.now() * 0.009);
      const sources: Array<[number, number, number, number]> = [
        [0.52, 0.24, 1.28, 0.85], // HP control-valve throttling dipole
        [0.0, 0.12, 4.75, 1.32], // S355 exhaust diffuser open-end dipole
      ];
      let idx = 0;
      for (let s = 0; s < sources.length; s++) {
        const [cx, cy, cz, baseR] = sources[s]!;
        for (let ring = 0; ring < 6; ring++) {
          const theta = ((ring + 0.5) / 6) * Math.PI;
          const st = Math.sin(theta);
          const ct = Math.cos(theta);
          for (let k = 0; k < 12 && idx < this.vfAcousticInst.count; k++) {
            const phi = (k / 12) * Math.PI * 2;
            dir.set(st * Math.cos(phi), ct, st * Math.sin(phi)).normalize();
            const rr = baseR + 0.18 * (s === 0 ? pulse : 1.04 / pulse);
            this.tmpPos.set(cx + rr * dir.x, cy + rr * dir.y, cz + rr * dir.z);
            this.tmpQuat.setFromUnitVectors(this.upVec, dir);
            this.tmpScale.set(0.85, splScale * pulse, 0.85);
            this.tmpMat.compose(this.tmpPos, this.tmpQuat, this.tmpScale);
            this.vfAcousticInst.setMatrixAt(idx, this.tmpMat);
            this.tmpColor.set(s === 0 ? "#38bdf8" : "#a78bfa");
            this.vfAcousticInst.setColorAt(idx, this.tmpColor);
            idx++;
          }
        }
      }
      this.vfAcousticInst.instanceMatrix.needsUpdate = true;
      if (this.vfAcousticInst.instanceColor) this.vfAcousticInst.instanceColor.needsUpdate = true;
    }
  }

  private handleResize(): void {
    const w = this.container.clientWidth || 1;
    const h = this.container.clientHeight || 1;
    this.camera.aspect = w / h;
    this.camera.updateProjectionMatrix();
    this.renderer.setSize(w, h);
    this.composer.setSize(w, h);
    this.bloomPass.setSize(w, h);
    this.smaaPass.setSize(w, h);
  }

  public dispose(): void {
    this.resizeObserver.disconnect();
    if (this.boundKeyDown) window.removeEventListener("keydown", this.boundKeyDown);
    if (this.boundKeyUp) window.removeEventListener("keyup", this.boundKeyUp);
    if (this.boundBlur) window.removeEventListener("blur", this.boundBlur);
    if (this.boundPointerUp) window.removeEventListener("pointerup", this.boundPointerUp);
    this.controls.dispose();
    this.envRT.dispose();
    this.fieldATex.dispose();
    this.fieldBTex.dispose();
    this.blackbodyTex.dispose();
    Object.values(this.textures).forEach((t) => t.dispose());
    this.composer.dispose();
    this.scene.traverse((obj) => {
      if (obj instanceof THREE.Mesh || obj instanceof THREE.Points || obj instanceof THREE.Line) {
        obj.geometry.dispose();
        if (Array.isArray(obj.material)) obj.material.forEach((m) => m.dispose());
        else obj.material.dispose();
      }
    });
    this.renderer.dispose();
    this.renderer.domElement.remove();
  }
}
