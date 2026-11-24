/**
 * Shared Three.js PBR studio kit, procedural PBR material upgrade pass,
 * GPU instanced GLB loader, 29-bone SkeletonUtils rig cloner, scroll-safe
 * OrbitControls, and scientific thermal colormap helpers.
 */

import * as THREE from "three";
import { OrbitControls } from "three/examples/jsm/controls/OrbitControls.js";
import { GLTFLoader } from "three/examples/jsm/loaders/GLTFLoader.js";
import * as SkeletonUtils from "three/examples/jsm/utils/SkeletonUtils.js";
import {
  getBrushedSteelMaps,
  getCastIronMaps,
  getChefTwillBump,
  getQuarryTileMaps,
  getWispyCloudTexture,
} from "./textures";

export const ASSET_BASE = `${import.meta.env.BASE_URL}kitchen-assets`;

const gltfLoader = new GLTFLoader();
const equipmentTemplateCache = new Map<string, Promise<THREE.Group | null>>();
const staffTemplateCache = new Map<string, Promise<{ scene: THREE.Group; bytes: number } | null>>();

/** Exact coordinate transform from Blender Z-up meters (manifest) to glTF2 / Three.js Y-up meters. */
export function blenderToThree(bx: number, by: number, bz: number): THREE.Vector3 {
  return new THREE.Vector3(bx, bz, -by);
}

/** Create an ACESFilmic PBR WebGLRenderer with crisp PCFSoftShadowMap shadows. */
export function createStudioRenderer(canvas: HTMLCanvasElement, alpha = true): THREE.WebGLRenderer | null {
  try {
    const renderer = new THREE.WebGLRenderer({
      canvas,
      antialias: true,
      alpha,
      powerPreference: "high-performance",
    });
    renderer.outputColorSpace = THREE.SRGBColorSpace;
    renderer.toneMapping = THREE.ACESFilmicToneMapping;
    renderer.toneMappingExposure = 1.06;
    renderer.shadowMap.enabled = true;
    renderer.shadowMap.type = THREE.PCFSoftShadowMap;
    renderer.setPixelRatio(Math.min(window.devicePixelRatio || 1, 2));
    return renderer;
  } catch {
    const note = document.createElement("p");
    note.className = "webgl-fallback";
    note.textContent = "Interactive 3D view requires WebGL support in this browser.";
    canvas.replaceWith(note);
    return null;
  }
}

/**
 * Build a balanced warm-culinary studio environment map for PMREM reflections:
 * soft diffused overhead canopy + warm tungsten & cool steel side softboxes + rich
 * dark-umber floor bounce so metallic stainless reads sculpted on EVERY face without
 * blinding #FFFFFF top blowouts or #000000 front-face crush.
 */
function buildCulinaryEnvironmentScene(): THREE.Scene {
  const envScene = new THREE.Scene();
  envScene.background = new THREE.Color(0x231d19);

  // Soft warm studio shell interior so reflections never drop to pure black
  const shell = new THREE.Mesh(
    new THREE.BoxGeometry(24, 14, 24),
    new THREE.MeshBasicMaterial({ color: 0x2a2420, side: THREE.BackSide }),
  );
  envScene.add(shell);

  const addPanel = (
    w: number,
    h: number,
    colorHex: number,
    pos: readonly [number, number, number],
    lookAt: readonly [number, number, number] = [0, 0, 0],
  ) => {
    const m = new THREE.Mesh(
      new THREE.PlaneGeometry(w, h),
      new THREE.MeshBasicMaterial({ color: colorHex, side: THREE.DoubleSide }),
    );
    m.position.set(pos[0], pos[1], pos[2]);
    m.lookAt(lookAt[0], lookAt[1], lookAt[2]);
    envScene.add(m);
  };

  // Controlled diffused overhead softbox (warm satin reflection, never pure white clip)
  addPanel(10, 10, 0xd8cec2, [0, 6.8, 0], [0, 0, 0]);
  // Warm tungsten front-right studio bank illuminating vertical cabinet faces
  addPanel(10, 6, 0xffe4c4, [9.5, 2.2, 8.5], [0, 0, 0]);
  // Cool architectural steel left-front fill bank
  addPanel(9, 6, 0xb8d4ee, [-9.5, 2.4, 7.5], [0, 0, 0]);
  // Warm copper rear rim strip for crisp silhouette speculars
  addPanel(11, 4.5, 0xf2b578, [0, 3.0, -10.5], [0, 0, 0]);
  // Soft warm floor bounce panel so undersides and bevels carry depth
  addPanel(14, 14, 0x3b2f27, [0, -6.5, 0], [0, 0, 0]);

  return envScene;
}

/**
 * Generate deterministic box-projected or top-down planar UV coordinates for Blender
 * bmesh geometries that were exported without a TEXCOORD_0 channel so procedural
 * Canvas PBR maps (brushed stainless, quarry tile, marbling, sourdough crust) render
 * crisply across every triangle without stretching.
 */
export function ensureMeshUVs(
  geom: THREE.BufferGeometry,
  mode: "box" | "planar" = "box",
  worldTileMeters = 0.65,
): void {
  if (!geom || !geom.attributes.position) return;
  if (mode === "box" && geom.attributes.uv) return;

  const pos = geom.attributes.position;
  if (!geom.attributes.normal) geom.computeVertexNormals();
  const norm = geom.attributes.normal;

  const uvs = new Float32Array(pos.count * 2);
  if (mode === "planar") {
    geom.computeBoundingBox();
    const bb = geom.boundingBox;
    const minX = bb ? bb.min.x : -0.1;
    const minZ = bb ? bb.min.z : -0.1;
    const spanX = Math.max(1e-4, bb ? bb.max.x - bb.min.x : 0.2);
    const spanZ = Math.max(1e-4, bb ? bb.max.z - bb.min.z : 0.2);
    for (let i = 0; i < pos.count; i += 1) {
      uvs[i * 2] = (pos.getX(i) - minX) / spanX;
      uvs[i * 2 + 1] = (pos.getZ(i) - minZ) / spanZ;
    }
  } else {
    const inv = 1 / Math.max(0.05, worldTileMeters);
    for (let i = 0; i < pos.count; i += 1) {
      const x = pos.getX(i);
      const y = pos.getY(i);
      const z = pos.getZ(i);
      const nx = Math.abs(norm ? norm.getX(i) : 0);
      const ny = Math.abs(norm ? norm.getY(i) : 1);
      const nz = Math.abs(norm ? norm.getZ(i) : 0);
      if (ny >= nx && ny >= nz) {
        uvs[i * 2] = x * inv;
        uvs[i * 2 + 1] = z * inv;
      } else if (nx >= nz) {
        uvs[i * 2] = z * inv;
        uvs[i * 2 + 1] = y * inv;
      } else {
        uvs[i * 2] = x * inv;
        uvs[i * 2 + 1] = y * inv;
      }
    }
  }
  geom.setAttribute("uv", new THREE.BufferAttribute(uvs, 2));
}

/**
 * Calibrated dark-culinary studio lighting rig (zero blown-out white clipping,
 * crisp sculpted metallic speculars, warm tungsten key + cool steel rim + soft fill).
 */
export function createStudioScene(renderer: THREE.WebGLRenderer, envIntensity = 0.72): THREE.Scene {
  const scene = new THREE.Scene();
  const pmrem = new THREE.PMREMGenerator(renderer);
  const envScene = buildCulinaryEnvironmentScene();
  const envTex = pmrem.fromScene(envScene, 0.04).texture;
  scene.environment = envTex;
  scene.environmentIntensity = envIntensity;
  pmrem.dispose();

  // Balanced warm sky / forged-umber ground ambient (keeps shadows rich, never pure black)
  const hemi = new THREE.HemisphereLight(0xffeed9, 0x45352b, 0.62);
  scene.add(hemi);

  // Sculpted warm tungsten key light casting crisp 2048^2 contact shadows
  const key = new THREE.DirectionalLight(0xffedd8, 1.18);
  key.position.set(8.5, 11.0, 7.5);
  key.castShadow = true;
  key.shadow.mapSize.set(2048, 2048);
  key.shadow.bias = -0.00035;
  key.shadow.normalBias = 0.024;
  key.shadow.camera.near = 0.5;
  key.shadow.camera.far = 48;
  key.shadow.camera.left = -16;
  key.shadow.camera.right = 16;
  key.shadow.camera.top = 16;
  key.shadow.camera.bottom = -16;
  scene.add(key);

  // Cool brushed-steel architectural rim from opposite rear corner for silhouette separation
  const coolRim = new THREE.DirectionalLight(0x93c5fd, 0.74);
  coolRim.position.set(-10.5, 6.2, -8.5);
  scene.add(coolRim);

  // Front-grazing warm studio fill so vertical stainless cabinet doors show crisp brushed detail
  const frontFill = new THREE.DirectionalLight(0xffd6b0, 0.72);
  frontFill.position.set(-6.5, 3.6, 11.5);
  scene.add(frontFill);

  return scene;
}

/**
 * Scroll-safe OrbitControls: wheel zoom is OFF by default so normal vertical page
 * scrolling is NEVER hijacked; clicking/focusing the canvas unlocks zoom until blur/leave.
 */
export function createScrollSafeOrbit(
  camera: THREE.PerspectiveCamera,
  canvas: HTMLCanvasElement,
  target: readonly [number, number, number],
  opts?: { minDistance?: number; maxDistance?: number; maxPolarAngle?: number },
): OrbitControls {
  const controls = new OrbitControls(camera, canvas);
  controls.target.set(target[0], target[1], target[2]);
  controls.enableDamping = true;
  controls.dampingFactor = 0.08;
  controls.enableZoom = false;
  controls.enablePan = false;
  controls.minDistance = opts?.minDistance ?? 0.5;
  controls.maxDistance = opts?.maxDistance ?? 35;
  controls.maxPolarAngle = opts?.maxPolarAngle ?? Math.PI * 0.49;

  const activateZoom = () => {
    controls.enableZoom = true;
  };
  const deactivateZoom = () => {
    controls.enableZoom = false;
  };
  canvas.addEventListener("pointerdown", activateZoom);
  canvas.addEventListener("focus", activateZoom);
  canvas.addEventListener("mouseleave", deactivateZoom);
  canvas.addEventListener("blur", deactivateZoom);

  controls.update();
  return controls;
}

/** Scientific ironbow / blackbody thermal false-color ramp in linear THREE.Color. */
export function thermalColor(tempC: number, minC = 0, maxC = 320, target = new THREE.Color()): THREE.Color {
  const t = Math.max(0, Math.min(1, (tempC - minC) / Math.max(1e-3, maxC - minC)));
  if (t < 0.25) {
    const u = t / 0.25;
    return target.setRGB(0.06 + u * 0.06, 0.14 + u * 0.24, 0.38 + u * 0.34);
  }
  if (t < 0.52) {
    const u = (t - 0.25) / 0.27;
    return target.setRGB(0.12 + u * 0.68, 0.38 - u * 0.22, 0.72 - u * 0.46);
  }
  if (t < 0.8) {
    const u = (t - 0.52) / 0.28;
    return target.setRGB(0.8 + u * 0.18, 0.16 + u * 0.48, 0.26 - u * 0.22);
  }
  const u = (t - 0.8) / 0.2;
  return target.setRGB(0.98 + u * 0.02, 0.64 + u * 0.32, 0.04 + u * 0.72);
}

/** Soft organic wispy cloud sprite texture for steam, smoke, and exhaust plumes. */
export function getSoftParticleTexture(): THREE.CanvasTexture {
  return getWispyCloudTexture();
}

export interface Scene3dFixture {
  readonly archetype_id: string;
  readonly fixture_id: string;
  readonly gltf_node?: string;
  readonly pose: { x_m: number; y_m: number; z_m: number; yaw_rad: number };
  readonly local_extents_m: readonly [number, number, number];
  readonly world_aabb_m?: { min_m: readonly [number, number, number]; max_m: readonly [number, number, number] };
  readonly swept_volumes?: readonly {
    volume_id: string;
    kind: string;
    aabb_world_m: { min_m: readonly [number, number, number]; max_m: readonly [number, number, number] };
    clearance_required_m: number;
  }[];
  readonly triangle_count?: number;
}

export interface Scene3dManifest {
  readonly scene_id: string;
  readonly archetype_id: string;
  readonly room_shell: { width_m: number; depth_m: number; ceiling_height_m: number };
  readonly triangle_count: number;
  readonly fixtures: readonly Scene3dFixture[];
}

export interface SiteMetadata {
  readonly scene_id: string;
  readonly walk_mesh_corridors?: readonly {
    corridor_id: string;
    min_m: readonly [number, number, number];
    max_m: readonly [number, number, number];
  }[];
  readonly aisle_clearances_m?: Record<string, unknown>;
}

export interface LoadedExemplarScene {
  readonly dir: string;
  readonly manifest: Scene3dManifest;
  readonly siteMeta: SiteMetadata | null;
  readonly sceneGroup: THREE.Group;
  readonly allMeshes: THREE.Mesh[];
  readonly hotMeshes: THREE.Mesh[];
  readonly coldMeshes: THREE.Mesh[];
  readonly drawCallCount: number;
  readonly triangleCount: number;
}

function classifyHeatRole(arch: string): "hot" | "cold" | "neutral" {
  const a = arch.toUpperCase();
  if (
    a.includes("RANGE") ||
    a.includes("EURO-ISLAND") ||
    a.includes("GRIDDLE") ||
    a.includes("CHARBROILER") ||
    a.includes("FRYER") ||
    a.includes("INDUCTION") ||
    a.includes("COMBI") ||
    a.includes("HEATLAMP") ||
    a.includes("TOASTER") ||
    a.includes("SALAMANDER") ||
    a.includes("TILT-SKILLET") ||
    a.includes("KETTLE") ||
    a.includes("STEAMTABLE")
  ) {
    return "hot";
  }
  if (
    a.includes("WALKIN") ||
    a.includes("REACHIN") ||
    a.includes("UNDERCOUNTER") ||
    a.includes("BLAST") ||
    a.includes("MAKELINE") ||
    a.includes("ICEMACHINE")
  ) {
    return "cold";
  }
  return "neutral";
}

/**
 * Upgrade exported Blender glTF materials with tactile procedural PBR canvas textures
 * and calibrated metallic/roughness contrast (eliminating pure-white clipped floors and
 * flat untextured stainless faces).
 */
export function upgradeKitchenMaterial(mat: THREE.Material, meshName = ""): void {
  if (!(mat instanceof THREE.MeshStandardMaterial)) return;
  const mName = (mat.name || "").toLowerCase();
  const oName = meshName.toLowerCase();

  if (mName.includes("quarry_tile") || oName.includes("arch_floor")) {
    const q = getQuarryTileMaps();
    mat.color.setHex(0xffffff);
    mat.map = q.map;
    mat.bumpMap = q.bumpMap;
    mat.bumpScale = 0.008;
    mat.roughnessMap = q.roughnessMap ?? null;
    mat.roughness = 0.78;
    mat.metalness = 0.04;
    mat.envMapIntensity = 0.55;
  } else if (mName.includes("cast_iron") || mName.includes("carbon_steel")) {
    const ci = getCastIronMaps();
    mat.color.setHex(mName.includes("bare") ? 0x3b3836 : 0x242120);
    mat.map = ci.map;
    mat.bumpMap = ci.bumpMap;
    mat.bumpScale = 0.004;
    mat.roughness = 0.56;
    mat.metalness = 0.62;
    mat.envMapIntensity = 1.05;
  } else if (mName.includes("copper") || mName.includes("brass")) {
    mat.color.setHex(mName.includes("brass") ? 0xd4a438 : 0xc96c3b);
    mat.roughness = 0.28;
    mat.metalness = 0.86;
    mat.envMapIntensity = 1.35;
  } else if (mName.includes("ss304") || mName.includes("stainless") || mName.includes("aluminum")) {
    const st = getBrushedSteelMaps();
    mat.color.setHex(mName.includes("backsplash") ? 0xb4bcc4 : 0x9fa8b2);
    mat.map = st.map;
    mat.bumpMap = st.bumpMap;
    mat.bumpScale = 0.0028;
    mat.roughnessMap = st.roughnessMap ?? null;
    mat.roughness = 0.4;
    mat.metalness = 0.44;
    mat.envMapIntensity = 0.95;
  } else if (mName.includes("thermal_gas_flame")) {
    mat.color.setHex(0x38bdf8);
    mat.emissive.setHex(0x2563eb);
    mat.emissiveIntensity = 2.6;
    mat.roughness = 0.2;
    mat.metalness = 0.0;
  } else if (mName.includes("thermal_calrod")) {
    mat.color.setHex(0x451a10);
    mat.emissive.setHex(0xf97316);
    mat.emissiveIntensity = 2.2;
    mat.roughness = 0.4;
    mat.metalness = 0.4;
  } else if (mName.includes("thermal_induction")) {
    mat.color.setHex(0x111216);
    mat.roughness = 0.14;
    mat.metalness = 0.15;
    mat.envMapIntensity = 1.3;
  } else if (mName.includes("digital_led")) {
    mat.color.setHex(0x0c2414);
    mat.emissive.setHex(0x22c55e);
    mat.emissiveIntensity = 1.8;
  } else if (mName.includes("liquid_oil")) {
    mat.color.setHex(0xb4690e);
    mat.emissive.setHex(0x451a03);
    mat.emissiveIntensity = 0.35;
    mat.roughness = 0.14;
    mat.metalness = 0.05;
  } else if (mName.includes("tempered_oven_glass") || mName.includes("polycarbonate") || mName.includes("pvc_strip")) {
    mat.color.setHex(0x94b8c8);
    mat.roughness = 0.12;
    mat.metalness = 0.08;
    mat.transparent = true;
    mat.opacity = mName.includes("pvc") ? 0.36 : 0.42;
    mat.depthWrite = false;
  } else if (mName.includes("polymer") || mName.includes("carrier_poly")) {
    mat.color.setHex(mName.includes("carrier") ? 0x2d333b : 0x191817);
    mat.roughness = 0.54;
    mat.metalness = 0.08;
    mat.envMapIntensity = 0.85;
  } else if (mName.includes("hdpe")) {
    mat.color.setHex(0xd6d2c8);
    mat.roughness = 0.56;
    mat.metalness = 0.0;
    mat.envMapIntensity = 0.6;
  } else if (mName.includes("wine_burgundy")) {
    mat.color.setHex(0x3a0812);
    mat.emissive.setHex(0x1e0308);
    mat.emissiveIntensity = 0.25;
    mat.roughness = 0.08;
    mat.metalness = 0.1;
    mat.transparent = true;
    mat.opacity = 0.88;
    mat.envMapIntensity = 1.6;
  } else {
    // Generic safety clamp so untextured bright colors never clip white under studio lighting
    const lum = mat.color.r * 0.299 + mat.color.g * 0.587 + mat.color.b * 0.114;
    if (lum > 0.78) mat.color.multiplyScalar(0.74 / lum);
    mat.roughness = Math.max(mat.roughness, 0.36);
    mat.metalness = Math.min(mat.metalness, 0.72);
    mat.envMapIntensity = 0.95;
  }

  mat.userData.baseColor = mat.color.clone();
  mat.userData.baseMap = mat.map;
  mat.userData.baseEmissive = mat.emissive.clone();
  mat.userData.baseEmissiveIntensity = mat.emissiveIntensity;
  mat.userData.baseRoughness = mat.roughness;
  mat.userData.baseMetalness = mat.metalness;
  mat.needsUpdate = true;
}

/** Load one canonical Blender `.glb` kitchen scene with true GPU InstancedMesh batching (<= 180 draw calls). */
export async function loadExemplarScene(dir: string): Promise<LoadedExemplarScene | null> {
  const manifestUrl = `${ASSET_BASE}/exemplars/${dir}/${dir}.scene3d.v1.json`;
  const siteUrl = `${ASSET_BASE}/exemplars/${dir}/${dir}.site_metadata.json`;
  const glbUrl = `${ASSET_BASE}/exemplars/${dir}/${dir}.glb`;

  const [mResp, sResp, bResp] = await Promise.all([fetch(manifestUrl), fetch(siteUrl), fetch(glbUrl)]);
  if (!mResp.ok || !bResp.ok) return null;

  const manifest = (await mResp.json()) as Scene3dManifest;
  const siteMeta = sResp.ok ? ((await sResp.json()) as SiteMetadata) : null;
  const buf = await bResp.arrayBuffer();
  const gltf = await new Promise<{ scene: THREE.Group }>((resolve, reject) => {
    gltfLoader.parse(buf, "", (g) => resolve(g as unknown as { scene: THREE.Group }), reject);
  });

  const sceneGroup = gltf.scene ?? new THREE.Group();
  sceneGroup.name = `EXEMPLAR_${dir}`;
  sceneGroup.updateMatrixWorld(true);

  const nodeToFx = new Map<string, Scene3dFixture>();
  for (const fx of manifest.fixtures ?? []) {
    if (fx.gltf_node) nodeToFx.set(fx.gltf_node, fx);
  }

  interface RawEntry {
    obj: THREE.Mesh;
    role: "hot" | "cold" | "neutral" | "wall" | "floor";
    fx: Scene3dFixture | null;
  }
  const rawMeshes: RawEntry[] = [];
  const toHide: THREE.Object3D[] = [];
  sceneGroup.traverse((obj) => {
    if (!(obj instanceof THREE.Mesh) || obj instanceof THREE.SkinnedMesh) return;
    let cur: THREE.Object3D | null = obj;
    let fx: Scene3dFixture | null = null;
    let isCeiling = false;
    let isWall = false;
    let isFloor = false;
    while (cur) {
      const nm = (cur.name || "").toLowerCase();
      if (nm.includes("arch_ceiling") || nm.includes("ceiling_slab")) {
        isCeiling = true;
        break;
      }
      if (nm.includes("arch_walls") || nm.includes("walls_frp")) {
        isWall = true;
      }
      if (nm.includes("arch_floor")) {
        isFloor = true;
      }
      if (cur.name && nodeToFx.has(cur.name)) {
        fx = nodeToFx.get(cur.name) ?? null;
        break;
      }
      cur = cur.parent;
    }
    if (isCeiling) {
      toHide.push(obj);
      return;
    }
    rawMeshes.push({
      obj,
      role: isWall ? "wall" : isFloor ? "floor" : fx ? classifyHeatRole(fx.archetype_id) : "neutral",
      fx,
    });
  });
  for (const h of toHide) {
    h.visible = false;
    h.removeFromParent();
  }

  // Group identical (geometry, material, heatRole) primitives into GPU InstancedMesh batches
  const batchGroups = new Map<string, RawEntry[]>();
  const geomIds = new WeakMap<THREE.BufferGeometry, number>();
  const matIds = new WeakMap<THREE.Material, number>();
  let nextG = 1;
  let nextM = 1;
  for (const item of rawMeshes) {
    if (item.role === "wall" || item.role === "floor") continue;
    const m = item.obj.material;
    if (!m || Array.isArray(m) || !item.obj.geometry) continue;
    let gId = geomIds.get(item.obj.geometry);
    if (!gId) {
      gId = nextG++;
      geomIds.set(item.obj.geometry, gId);
    }
    let mId = matIds.get(m);
    if (!mId) {
      mId = nextM++;
      matIds.set(m, mId);
    }
    const key = `${gId}:${mId}:${item.role}`;
    const list = batchGroups.get(key) ?? [];
    list.push(item);
    batchGroups.set(key, list);
  }

  const allMeshes: THREE.Mesh[] = [];
  const hotMeshes: THREE.Mesh[] = [];
  const coldMeshes: THREE.Mesh[] = [];
  const instancedContainer = new THREE.Group();
  instancedContainer.name = "GPU_INSTANCED_STATIC_BATCH";
  const invScene = new THREE.Matrix4().copy(sceneGroup.matrixWorld).invert();
  const localMat = new THREE.Matrix4();
  const handled = new Set<THREE.Mesh>();

  for (const list of batchGroups.values()) {
    if (list.length < 2) continue;
    const first = list[0];
    ensureMeshUVs(first.obj.geometry, "box", 0.65);
    const clonedMat = (first.obj.material as THREE.Material).clone();
    upgradeKitchenMaterial(clonedMat, first.obj.name);
    const inst = new THREE.InstancedMesh(first.obj.geometry, clonedMat, list.length);
    inst.castShadow = true;
    inst.receiveShadow = true;
    inst.frustumCulled = false;
    for (let i = 0; i < list.length; i += 1) {
      localMat.multiplyMatrices(invScene, list[i].obj.matrixWorld);
      inst.setMatrixAt(i, localMat);
      handled.add(list[i].obj);
    }
    inst.instanceMatrix.needsUpdate = true;
    instancedContainer.add(inst);
    allMeshes.push(inst);
    if (first.role === "hot") hotMeshes.push(inst);
    if (first.role === "cold") coldMeshes.push(inst);
  }

  for (const entry of rawMeshes) {
    if (handled.has(entry.obj)) continue;
    if (entry.role === "wall") {
      // Smoked translucent architectural cutaway walls + crisp bronze structural wireframe edge
      entry.obj.castShadow = false;
      entry.obj.receiveShadow = false;
      entry.obj.material = new THREE.MeshStandardMaterial({
        color: 0x2c221b,
        roughness: 0.55,
        metalness: 0.3,
        transparent: true,
        opacity: 0.11,
        depthWrite: false,
        side: THREE.DoubleSide,
      });
      continue;
    }
    if (entry.role === "floor") {
      entry.obj.castShadow = false;
      entry.obj.receiveShadow = true;
      ensureMeshUVs(entry.obj.geometry, "box", 0.62);
      // True 6"x6" commercial terracotta-charcoal quarry tile pattern scaled in world meters
      const q = getQuarryTileMaps();
      const floorColMap = q.map.clone();
      floorColMap.needsUpdate = true;
      floorColMap.repeat.set(1, 1);
      const floorBumpMap = q.bumpMap.clone();
      floorBumpMap.needsUpdate = true;
      floorBumpMap.repeat.set(1, 1);
      const floorRoughMap = q.roughnessMap?.clone() ?? null;
      if (floorRoughMap) {
        floorRoughMap.needsUpdate = true;
        floorRoughMap.repeat.set(1, 1);
      }
      const floorMat = new THREE.MeshStandardMaterial({
        color: 0xffffff,
        map: floorColMap,
        bumpMap: floorBumpMap,
        bumpScale: 0.008,
        roughnessMap: floorRoughMap,
        roughness: 0.82,
        metalness: 0.05,
        envMapIntensity: 0.45,
      });
      floorMat.userData.baseColor = floorMat.color.clone();
      floorMat.userData.baseMap = floorMat.map;
      floorMat.userData.baseEmissive = floorMat.emissive.clone();
      floorMat.userData.baseEmissiveIntensity = 0;
      floorMat.userData.baseRoughness = floorMat.roughness;
      floorMat.userData.baseMetalness = floorMat.metalness;
      entry.obj.material = floorMat;
      allMeshes.push(entry.obj);
      continue;
    }
    ensureMeshUVs(entry.obj.geometry, "box", 0.65);
    entry.obj.castShadow = true;
    entry.obj.receiveShadow = true;
    if (Array.isArray(entry.obj.material)) {
      entry.obj.material = entry.obj.material.map((m) => {
        const c = m.clone();
        upgradeKitchenMaterial(c, entry.obj.name);
        return c;
      });
    } else if (entry.obj.material) {
      const c = entry.obj.material.clone();
      upgradeKitchenMaterial(c, entry.obj.name);
      entry.obj.material = c;
    }
    allMeshes.push(entry.obj);
    if (entry.role === "hot") hotMeshes.push(entry.obj);
    if (entry.role === "cold") coldMeshes.push(entry.obj);
  }

  for (const obj of handled) obj.removeFromParent();
  if (instancedContainer.children.length > 0) sceneGroup.add(instancedContainer);

  let tris = 0;
  for (const m of allMeshes) {
    const g = m.geometry;
    const count = g.index ? g.index.count / 3 : (g.attributes.position?.count ?? 0) / 3;
    const mult = m instanceof THREE.InstancedMesh ? m.count : 1;
    tris += Math.floor(count) * mult;
  }

  return {
    dir,
    manifest,
    siteMeta,
    sceneGroup,
    allMeshes,
    hotMeshes,
    coldMeshes,
    drawCallCount: allMeshes.length,
    triangleCount: tris,
  };
}

/** Load and clone a single Blender equipment / cookware / food `.glb` centered on X/Z with base at Y = 0. */
export async function loadEquipmentModel(glbFile: string): Promise<THREE.Group | null> {
  let prom = equipmentTemplateCache.get(glbFile);
  if (!prom) {
    prom = (async () => {
      const resp = await fetch(`${ASSET_BASE}/equipment/${glbFile}`);
      if (!resp.ok) return null;
      const buf = await resp.arrayBuffer();
      const gltf = await new Promise<{ scene: THREE.Group }>((resolve, reject) => {
        gltfLoader.parse(buf, "", (g) => resolve(g as unknown as { scene: THREE.Group }), reject);
      });
      return gltf.scene ?? null;
    })();
    equipmentTemplateCache.set(glbFile, prom);
  }
  const tpl = await prom;
  if (!tpl) return null;
  const clone = tpl.clone(true);
  clone.traverse((obj) => {
    if (obj instanceof THREE.Mesh) {
      ensureMeshUVs(obj.geometry, "box", 0.55);
      obj.castShadow = true;
      obj.receiveShadow = true;
      if (Array.isArray(obj.material)) {
        obj.material = obj.material.map((m) => {
          const c = m.clone();
          upgradeKitchenMaterial(c, obj.name);
          return c;
        });
      } else if (obj.material) {
        const c = obj.material.clone();
        upgradeKitchenMaterial(c, obj.name);
        obj.material = c;
      }
    }
  });
  clone.updateMatrixWorld(true);
  const box = new THREE.Box3().setFromObject(clone);
  const center = box.getCenter(new THREE.Vector3());
  clone.position.set(-center.x, -box.min.y, -center.z);
  const wrapper = new THREE.Group();
  wrapper.add(clone);
  return wrapper;
}

export interface StaffRigInstance {
  readonly cohort: "p05" | "p50" | "p95";
  readonly root: THREE.Group;
  readonly bones: Map<string, THREE.Bone>;
}

/** Load and deep-clone one of the 3 articulated 29-bone SkinnedMesh chef rigs (`p05`, `p50`, `p95`). */
export async function cloneStaffRig(cohort: "p05" | "p50" | "p95"): Promise<StaffRigInstance | null> {
  let prom = staffTemplateCache.get(cohort);
  if (!prom) {
    prom = (async () => {
      const resp = await fetch(`${ASSET_BASE}/staff/chef_rig_${cohort}.glb`);
      if (!resp.ok) return null;
      const buf = await resp.arrayBuffer();
      const gltf = await new Promise<{ scene: THREE.Group }>((resolve, reject) => {
        gltfLoader.parse(buf, "", (g) => resolve(g as unknown as { scene: THREE.Group }), reject);
      });
      return { scene: gltf.scene, bytes: buf.byteLength };
    })();
    staffTemplateCache.set(cohort, prom);
  }
  const tpl = await prom;
  if (!tpl) return null;
  const root = SkeletonUtils.clone(tpl.scene) as THREE.Group;
  const bones = new Map<string, THREE.Bone>();
  const twillBump = getChefTwillBump();

  root.traverse((obj) => {
    if (obj instanceof THREE.Bone && obj.name) {
      obj.userData.restRot = obj.rotation.clone();
      bones.set(obj.name, obj);
    }
    if (obj instanceof THREE.Mesh) {
      obj.castShadow = true;
      obj.receiveShadow = true;
      const polishChefMat = (m: THREE.Material): THREE.Material => {
        const c = m.clone();
        if (c instanceof THREE.MeshStandardMaterial) {
          const nm = (c.name || "").toLowerCase();
          if (nm.includes("jacket") || nm.includes("apron") || nm.includes("pants")) {
            c.bumpMap = twillBump;
            c.bumpScale = 0.0025;
            c.roughness = Math.max(0.76, c.roughness);
            c.metalness = 0.0;
            if (nm.includes("jacket")) c.color.setHex(0xdad4ca);
            if (nm.includes("apron")) {
              c.color.setHex(0x1f2228);
              c.polygonOffset = true;
              c.polygonOffsetFactor = -2;
              c.polygonOffsetUnits = -2;
            }
            if (nm.includes("pants")) c.color.setHex(0x282c34);
          } else if (nm.includes("piping")) {
            c.color.setHex(0xd97706); // Warm copper brigade piping accent!
            c.roughness = 0.45;
          } else if (nm.includes("mannequin")) {
            c.color.setHex(0xc4b09e);
            c.roughness = 0.58;
            c.metalness = 0.0;
          } else if (nm.includes("clog")) {
            c.color.setHex(0x131316);
            c.roughness = 0.52;
          }
          c.envMapIntensity = 0.65;
          c.needsUpdate = true;
        }
        return c;
      };
      if (Array.isArray(obj.material)) {
        obj.material = obj.material.map(polishChefMat);
      } else if (obj.material) {
        obj.material = polishChefMat(obj.material);
      }
    }
  });
  return { cohort, root, bones };
}

/** Dispose all geometries, materials, and skeletons inside an Object3D subtree. */
export function disposeObjectTree(root: THREE.Object3D | null, preserveSharedGeometry = false): void {
  if (!root) return;
  const geos = new Set<THREE.BufferGeometry>();
  const mats = new Set<THREE.Material>();
  root.traverse((obj) => {
    if (obj instanceof THREE.SkinnedMesh && obj.skeleton) {
      obj.skeleton.dispose();
    }
    if (obj instanceof THREE.Mesh) {
      if (!preserveSharedGeometry && obj.geometry) geos.add(obj.geometry);
      const list = Array.isArray(obj.material) ? obj.material : obj.material ? [obj.material] : [];
      for (const m of list) mats.add(m);
    }
  });
  for (const g of geos) g.dispose();
  for (const m of mats) m.dispose();
}

/** Sync a WebGLRenderer and PerspectiveCamera cleanly to container ResizeObserver measurements. */
export function bindCanvasResize(
  container: HTMLElement,
  renderer: THREE.WebGLRenderer,
  camera: THREE.PerspectiveCamera,
): () => void {
  let lastW = 0;
  let lastH = 0;
  const resize = () => {
    const w = Math.max(1, Math.floor(container.clientWidth));
    const h = Math.max(1, Math.floor(container.clientHeight));
    if (Math.abs(w - lastW) <= 2 && Math.abs(h - lastH) <= 2) return;
    lastW = w;
    lastH = h;
    renderer.setSize(w, h, false);
    camera.aspect = w / h;
    camera.updateProjectionMatrix();
  };
  resize();
  const ro = new ResizeObserver(resize);
  ro.observe(container);
  return () => ro.disconnect();
}

/** Observe element viewport intersection so offscreen WebGL loops pause immediately. */
export function observeInView(el: Element, cb: (visible: boolean) => void, rootMargin = "40px 0px"): () => void {
  const io = new IntersectionObserver(
    (entries) => {
      const entry = entries[0];
      if (entry) cb(entry.isIntersecting);
    },
    { rootMargin },
  );
  io.observe(el);
  return () => io.disconnect();
}

export interface StripSeries {
  readonly label: string;
  readonly color: string;
  readonly values: readonly number[];
  readonly dashed?: boolean;
}

/** Crisp Retina 2D strip-chart renderer with subtle gradient area fill for live physics readouts. */
export function drawStripChart(
  canvas: HTMLCanvasElement,
  series: readonly StripSeries[],
  opts: { minY: number; maxY: number; unit: string; threshold?: { value: number; color: string; label: string } },
): void {
  const dpr = Math.min(window.devicePixelRatio || 1, 2);
  const w = Math.max(220, Math.floor(canvas.clientWidth || 320));
  const h = Math.max(96, Math.floor(canvas.clientHeight || 116));
  if (canvas.width !== Math.floor(w * dpr) || canvas.height !== Math.floor(h * dpr)) {
    canvas.width = Math.floor(w * dpr);
    canvas.height = Math.floor(h * dpr);
  }
  const ctx = canvas.getContext("2d");
  if (!ctx) return;
  ctx.save();
  ctx.scale(dpr, dpr);
  ctx.clearRect(0, 0, w, h);

  const padL = 38;
  const padR = 12;
  const padT = 10;
  const padB = 20;
  const pw = w - padL - padR;
  const ph = h - padT - padB;

  // Gridlines
  ctx.strokeStyle = "rgba(255, 224, 196, 0.09)";
  ctx.lineWidth = 1;
  ctx.font = '500 10px "JetBrains Mono Variable", "JetBrains Mono", monospace';
  ctx.fillStyle = "rgba(214, 198, 184, 0.62)";
  for (let i = 0; i <= 3; i += 1) {
    const f = i / 3;
    const y = padT + ph * (1 - f);
    const val = opts.minY + (opts.maxY - opts.minY) * f;
    ctx.beginPath();
    ctx.moveTo(padL, y);
    ctx.lineTo(padL + pw, y);
    ctx.stroke();
    ctx.fillText(`${Math.round(val)}`, 4, y + 3);
  }

  // Optional threshold line (e.g. smoke point or HACCP 41 F)
  if (opts.threshold) {
    const f = Math.max(0, Math.min(1, (opts.threshold.value - opts.minY) / Math.max(1e-3, opts.maxY - opts.minY)));
    const y = padT + ph * (1 - f);
    ctx.save();
    ctx.setLineDash([4, 4]);
    ctx.strokeStyle = opts.threshold.color;
    ctx.lineWidth = 1.25;
    ctx.beginPath();
    ctx.moveTo(padL, y);
    ctx.lineTo(padL + pw, y);
    ctx.stroke();
    ctx.restore();
    ctx.fillStyle = opts.threshold.color;
    ctx.fillText(opts.threshold.label, padL + 4, Math.max(padT + 9, y - 4));
  }

  for (let si = 0; si < series.length; si += 1) {
    const s = series[si];
    if (s.values.length < 2) continue;
    ctx.save();
    if (s.dashed) ctx.setLineDash([5, 3]);
    ctx.strokeStyle = s.color;
    ctx.lineWidth = 2.1;
    ctx.beginPath();
    const n = s.values.length;
    let lastX = padL;
    let lastY = padT + ph;
    for (let i = 0; i < n; i += 1) {
      const x = padL + (i / (n - 1)) * pw;
      const f = Math.max(0, Math.min(1, (s.values[i] - opts.minY) / Math.max(1e-3, opts.maxY - opts.minY)));
      const y = padT + ph * (1 - f);
      if (i === 0) ctx.moveTo(x, y);
      else ctx.lineTo(x, y);
      lastX = x;
      lastY = y;
    }
    ctx.stroke();

    // Subtle luminous area fill under primary series
    if (si === 0 && !s.dashed) {
      ctx.lineTo(padL + pw, padT + ph);
      ctx.lineTo(padL, padT + ph);
      ctx.closePath();
      const grad = ctx.createLinearGradient(0, padT, 0, padT + ph);
      grad.addColorStop(0, `${s.color}2e`);
      grad.addColorStop(1, `${s.color}00`);
      ctx.fillStyle = grad;
      ctx.fill();
    }

    // Crisp live endpoint dot
    ctx.setLineDash([]);
    ctx.fillStyle = s.color;
    ctx.beginPath();
    ctx.arc(lastX, lastY, 3, 0, Math.PI * 2);
    ctx.fill();
    ctx.restore();
  }

  // Inline legend at bottom
  let lx = padL;
  for (const s of series) {
    ctx.fillStyle = s.color;
    ctx.fillRect(lx, h - 10, 8, 3);
    ctx.fillStyle = "rgba(247, 241, 235, 0.84)";
    ctx.fillText(s.label, lx + 12, h - 5);
    lx += ctx.measureText(s.label).width + 24;
  }
  ctx.fillStyle = "rgba(214, 198, 184, 0.55)";
  ctx.fillText(opts.unit, w - padR - ctx.measureText(opts.unit).width, h - 5);
  ctx.restore();
}
