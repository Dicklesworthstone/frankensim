/**
 * HERO STAGE (#hero) — Full-Bleed Interactive 3D Commercial Kitchen Showcase.
 *
 * Loads all 5 real Blender-authored canonical commercial kitchen `.glb` scenes
 * with GPU InstancedMesh batching (<= 180 draw calls), quarry-tile floors &
 * brushed 304 stainless PBR, perforated anti-fatigue cookline mats, live glowing
 * burner flames & under-hood lighting, articulated 29-bone SkinnedMesh chef rigs
 * performing walk-and-cook station cycles while carrying utensils, 4 live 3D
 * physics overlays (PBR Stainless, Thermal IR, 3D Volumetric Exhaust Plumes,
 * 3D Floor Aisle Clearances & Door Swings), 5 camera presets, and the floating
 * frosted-glass telemetry console.
 */

import * as THREE from "three";
import {
  EXEMPLARS,
  computeKitchenTelemetry,
  type CameraPresetId,
  type ExemplarId,
  type ExemplarSpec,
  type StaffWaypointRoute,
} from "../sim/exemplars";
import {
  bindCanvasResize,
  blenderToThree,
  cloneStaffRig,
  createScrollSafeOrbit,
  createStudioRenderer,
  createStudioScene,
  disposeObjectTree,
  getSoftParticleTexture,
  loadExemplarScene,
  observeInView,
  thermalColor,
  type LoadedExemplarScene,
  type StaffRigInstance,
} from "./kit";
import { getAntiFatigueMatMaps } from "./textures";

type OverlayMode = "pbr" | "thermal" | "plumes" | "aisles";

interface WalkingChef {
  readonly rig: StaffRigInstance;
  readonly route: StaffWaypointRoute;
  readonly points: THREE.Vector3[];
  segmentIdx: number;
  segProgressM: number;
  dwellLeftS: number;
  currentYaw: number;
  readonly stridePhaseOffset: number;
}

interface PlumeParticleSystem {
  readonly points: THREE.Points;
  readonly origins: THREE.Vector3[];
  readonly phases: Float32Array;
  readonly speeds: Float32Array;
  readonly ceilingY: number;
}

interface LiveBurnerBeacon {
  readonly flameGroup: THREE.Group;
  readonly ringMesh: THREE.Mesh;
  readonly baseY: number;
  readonly isInductionCapable: boolean;
}

export interface HeroController {
  setExemplar(id: ExemplarId): Promise<void>;
  getLoadedExemplarId(): ExemplarId;
  getDrawCalls(): number;
}

function shortestAngleDelta(from: number, to: number): number {
  let d = (to - from) % (Math.PI * 2);
  if (d > Math.PI) d -= Math.PI * 2;
  if (d < -Math.PI) d += Math.PI * 2;
  return d;
}

/** Build a lightweight miniature utensil (skillet or spatula) held in a chef's right hand bone. */
function makeHeldChefProp(role: string): THREE.Group {
  const g = new THREE.Group();
  const r = role.toLowerCase();
  const steel = new THREE.MeshStandardMaterial({ color: 0xb8c2cc, roughness: 0.3, metalness: 0.8 });
  const dark = new THREE.MeshStandardMaterial({ color: 0x1c1917, roughness: 0.65, metalness: 0.4 });

  if (r.includes("grill") || r.includes("fry") || r.includes("saucier") || r.includes("hot") || r.includes("rôt")) {
    // Miniature sauté skillet held forward in right hand
    const handle = new THREE.Mesh(new THREE.BoxGeometry(0.02, 0.014, 0.16), steel);
    handle.position.set(0, 0, 0.08);
    const pan = new THREE.Mesh(new THREE.CylinderGeometry(0.11, 0.095, 0.036, 16), dark);
    pan.position.set(0, 0.004, 0.24);
    g.add(handle, pan);
  } else {
    // Stainless half-sheet / gastronorm tray or offset spatula for pass/prep/cold staff
    const tray = new THREE.Mesh(new THREE.BoxGeometry(0.24, 0.018, 0.17), steel);
    tray.position.set(0.02, 0.01, 0.14);
    g.add(tray);
  }
  return g;
}

export function mountHeroKitchen(root: HTMLElement): HeroController {
  const stageEl = root.querySelector<HTMLElement>("[data-hero-stage]")!;
  const canvas = root.querySelector<HTMLCanvasElement>("[data-hero-canvas]")!;
  const statusEl = root.querySelector<HTMLElement>("[data-hero-status]");

  const maybeRenderer = createStudioRenderer(canvas, true);
  if (!maybeRenderer) {
    return {
      setExemplar: async () => {},
      getLoadedExemplarId: () => "mall_qsr_400",
      getDrawCalls: () => 0,
    };
  }
  const renderer = maybeRenderer;

  const scene = createStudioScene(renderer, 0.64);
  const camera = new THREE.PerspectiveCamera(40, 16 / 9, 0.1, 140);
  const initialSpec = EXEMPLARS[0];
  camera.position.set(...initialSpec.cameras.isometric.pos);

  const controls = createScrollSafeOrbit(camera, canvas, initialSpec.cameras.isometric.target, {
    minDistance: 1.6,
    maxDistance: 44,
    maxPolarAngle: Math.PI * 0.485,
  });
  bindCanvasResize(stageEl, renderer, camera);

  // Subtle architectural floor grid around the kitchen footprint
  const floorGrid = new THREE.GridHelper(44, 44, 0x3f2d20, 0x19130f);
  floorGrid.position.set(4, -0.012, -6);
  scene.add(floorGrid);

  const kitchenAnchor = new THREE.Group();
  const detailAnchor = new THREE.Group();
  const staffAnchor = new THREE.Group();
  const overlayAnchor = new THREE.Group();
  scene.add(kitchenAnchor, detailAnchor, staffAnchor, overlayAnchor);

  let currentSpec: ExemplarSpec = initialSpec;
  let loadedScene: LoadedExemplarScene | null = null;
  let activeCameraPreset: CameraPresetId = "isometric";
  const activeOverlays = new Set<OverlayMode>(["pbr", "plumes"]);
  let walkingChefs: WalkingChef[] = [];
  let burnerBeacons: LiveBurnerBeacon[] = [];
  let plumeSystem: PlumeParticleSystem | null = null;
  let aisleGroup: THREE.Group | null = null;
  let paused = false;
  let inView = true;
  let loadSeq = 0;

  let camAnimT = 1;
  const camStartPos = new THREE.Vector3();
  const camGoalPos = new THREE.Vector3(...initialSpec.cameras.isometric.pos);
  const tgtStartPos = new THREE.Vector3();
  const tgtGoalPos = new THREE.Vector3(...initialSpec.cameras.isometric.target);

  let ordersPerHr = initialSpec.nominalOrdersPerHr;
  let electrificationPct = 35;
  let ambientC = 26;

  // DOM bindings
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
    for (const b of camBtns) {
      b.setAttribute("aria-pressed", String(b.dataset.cam === preset));
    }
    const pose = currentSpec.cameras[preset];
    camStartPos.copy(camera.position);
    tgtStartPos.copy(controls.target);
    camGoalPos.set(pose.pos[0], pose.pos[1], pose.pos[2]);
    tgtGoalPos.set(pose.target[0], pose.target[1], pose.target[2]);
    camAnimT = 0;
  }

  function applyThermalOverlay(): void {
    if (!loadedScene) return;
    const isThermal = activeOverlays.has("thermal");
    const loadRatio = ordersPerHr / currentSpec.nominalOrdersPerHr;
    const elecFrac = electrificationPct / 100;
    const hotTempC = 175 + 115 * loadRatio * (1 - 0.34 * elecFrac);
    const hotCol = thermalColor(hotTempC, 0, 320);
    const coldCol = thermalColor(4, 0, 320);
    const neutralCol = thermalColor(ambientC + 10, 0, 320);

    const setMat = (mesh: THREE.Mesh, role: "hot" | "cold" | "neutral") => {
      const mats = Array.isArray(mesh.material) ? mesh.material : [mesh.material];
      for (const m of mats) {
        if (!(m instanceof THREE.MeshStandardMaterial)) continue;
        if (isThermal) {
          const target = role === "hot" ? hotCol : role === "cold" ? coldCol : neutralCol;
          m.map = null;
          m.color.copy(target);
          m.emissive.copy(target);
          m.emissiveIntensity = role === "hot" ? 0.95 : role === "cold" ? 0.32 : 0.10;
          m.roughness = 0.65;
          m.metalness = 0.08;
        } else {
          if (m.userData.baseColor instanceof THREE.Color) m.color.copy(m.userData.baseColor);
          if ("baseMap" in m.userData) m.map = m.userData.baseMap as THREE.Texture | null;
          if (m.userData.baseEmissive instanceof THREE.Color) m.emissive.copy(m.userData.baseEmissive);
          const baseEmI = (m.userData.baseEmissiveIntensity as number | undefined) ?? 0;
          m.emissiveIntensity =
            role === "hot" ? Math.max(baseEmI, 0.08 + 0.14 * loadRatio) : baseEmI;
          if (role === "hot" && m.emissive.getHex() === 0) {
            m.emissive.setHex(elecFrac > 0.65 ? 0x0b3a4a : 0x451a05);
          }
          m.roughness = (m.userData.baseRoughness as number | undefined) ?? 0.38;
          m.metalness = (m.userData.baseMetalness as number | undefined) ?? 0.65;
        }
        m.needsUpdate = true;
      }
    };

    const hotSet = new Set(loadedScene.hotMeshes);
    const coldSet = new Set(loadedScene.coldMeshes);
    for (const mesh of loadedScene.allMeshes) {
      setMat(mesh, hotSet.has(mesh) ? "hot" : coldSet.has(mesh) ? "cold" : "neutral");
    }
  }

  /**
   * Populate realistic architectural & culinary life inside the loaded kitchen:
   * - Architectural copper-bronze room footprint outline frame
   * - Perforated black rubber anti-fatigue mats along primary walk corridors
   * - Live animated flame / induction glow beacons on every cooking deck
   * - Warm recessed under-hood & expo pass point lights + cool walk-in glow
   */
  function buildKitchenDetails(loaded: LoadedExemplarScene): LiveBurnerBeacon[] {
    const beacons: LiveBurnerBeacon[] = [];
    const { width_m: rw, depth_m: rd, ceiling_height_m: rh } = loaded.manifest.room_shell;

    // 1. Subtle architectural bronze envelope frame outlining cutaway room extents
    const roomOutline = new THREE.LineSegments(
      new THREE.EdgesGeometry(new THREE.BoxGeometry(rw, rh, rd)),
      new THREE.LineBasicMaterial({ color: 0x8a6443, transparent: true, opacity: 0.28 }),
    );
    roomOutline.position.copy(blenderToThree(rw * 0.5, rd * 0.5, rh * 0.5));
    detailAnchor.add(roomOutline);

    // 2. Perforated black rubber anti-fatigue mats inside walk corridors
    const matMaps = getAntiFatigueMatMaps();
    const matPbr = new THREE.MeshStandardMaterial({
      map: matMaps.map,
      bumpMap: matMaps.bumpMap,
      bumpScale: 0.006,
      roughness: 0.88,
      metalness: 0.02,
    });
    const corridors = loaded.siteMeta?.walk_mesh_corridors ?? [];
    for (let i = 0; i < Math.min(corridors.length, 3); i += 1) {
      const cor = corridors[i];
      const cw = cor.max_m[0] - cor.min_m[0];
      const cd = cor.max_m[1] - cor.min_m[1];
      const mw = Math.min(cw * 0.72, cw > cd ? Math.min(cw - 0.3, 4.2) : 0.78);
      const md = Math.min(cd * 0.72, cd >= cw ? Math.min(cd - 0.3, 4.2) : 0.78);
      if (mw < 0.45 || md < 0.45) continue;
      const cx = 0.5 * (cor.min_m[0] + cor.max_m[0]);
      const cy = 0.5 * (cor.min_m[1] + cor.max_m[1]);
      const matMesh = new THREE.Mesh(new THREE.BoxGeometry(mw, 0.012, md), matPbr);
      matMesh.position.copy(blenderToThree(cx, cy, 0.008));
      matMesh.receiveShadow = true;
      detailAnchor.add(matMesh);
    }

    // 3. Live burner glow beacons & under-hood warm culinary worklights (capped at 4 lights for fast perf)
    let warmLightsAdded = 0;
    let coldLightsAdded = 0;
    const flameConeGeo = new THREE.ConeGeometry(0.024, 0.068, 7);
    flameConeGeo.translate(0, 0.034, 0);

    for (const fx of loaded.manifest.fixtures) {
      const a = fx.archetype_id.toUpperCase();
      if (!fx.world_aabb_m) continue;
      const mn = fx.world_aabb_m.min_m;
      const mx = fx.world_aabb_m.max_m;
      const cx = 0.5 * (mn[0] + mx[0]);
      const cy = 0.5 * (mn[1] + mx[1]);

      if (a.includes("HOOD-TYPE1") && warmLightsAdded < 2) {
        const hoodLight = new THREE.PointLight(0xffb45c, 2.2, 6.5, 1.5);
        hoodLight.position.copy(blenderToThree(cx, cy, Math.max(1.75, mn[2] - 0.06)));
        detailAnchor.add(hoodLight);
        warmLightsAdded += 1;
      } else if (a.includes("HEATLAMP") && warmLightsAdded < 3) {
        const passLight = new THREE.PointLight(0xff6b2b, 1.8, 3.8, 1.6);
        passLight.position.copy(blenderToThree(cx, cy, Math.max(1.25, mn[2] - 0.04)));
        detailAnchor.add(passLight);
        warmLightsAdded += 1;
      } else if (a.includes("WALKIN") && coldLightsAdded < 1) {
        const coldLight = new THREE.PointLight(0x38bdf8, 1.4, 4.5, 1.6);
        coldLight.position.copy(blenderToThree(cx, cy, Math.min(2.1, mx[2] - 0.3)));
        detailAnchor.add(coldLight);
        coldLightsAdded += 1;
      }

      if (
        a.includes("RANGE") ||
        a.includes("CHARBROILER") ||
        a.includes("INDUCTION") ||
        a.includes("EURO-ISLAND") ||
        a.includes("GRIDDLE")
      ) {
        const deckZ = Math.min(mx[2], 0.94);
        const pos = blenderToThree(cx, cy, deckZ + 0.015);
        const flameGroup = new THREE.Group();
        flameGroup.position.copy(pos);

        // Luminous induction / ember ring on cooking surface
        const ringMesh = new THREE.Mesh(
          new THREE.RingGeometry(0.11, 0.18, 24),
          new THREE.MeshBasicMaterial({
            color: a.includes("INDUCTION") ? 0x22d3ee : 0xf97316,
            transparent: true,
            opacity: 0.65,
            side: THREE.DoubleSide,
            blending: THREE.AdditiveBlending,
            depthWrite: false,
          }),
        );
        ringMesh.rotation.x = -Math.PI / 2;
        flameGroup.add(ringMesh);

        // 4 subtle flickering premixed flame tongues
        for (let k = 0; k < 4; k += 1) {
          const ang = (k / 4) * Math.PI * 2 + 0.4;
          const tongue = new THREE.Mesh(
            flameConeGeo,
            new THREE.MeshBasicMaterial({
              color: k % 2 === 0 ? 0x60a5fa : 0xfb923c,
              transparent: true,
              opacity: 0.72,
              blending: THREE.AdditiveBlending,
              depthWrite: false,
            }),
          );
          tongue.position.set(Math.cos(ang) * 0.11, 0.0, Math.sin(ang) * 0.11);
          flameGroup.add(tongue);
        }
        detailAnchor.add(flameGroup);
        beacons.push({
          flameGroup,
          ringMesh,
          baseY: pos.y,
          isInductionCapable: !a.includes("CHARBROILER"),
        });
      }
    }

    return beacons;
  }

  function buildPlumeSystem(loaded: LoadedExemplarScene): PlumeParticleSystem {
    const origins: THREE.Vector3[] = [];
    for (const fx of loaded.manifest.fixtures) {
      const a = fx.archetype_id.toUpperCase();
      if (
        a.includes("RANGE") ||
        a.includes("GRIDDLE") ||
        a.includes("CHARBROILER") ||
        a.includes("FRYER") ||
        a.includes("EURO-ISLAND") ||
        a.includes("COMBI") ||
        a.includes("KETTLE") ||
        a.includes("TILT-SKILLET")
      ) {
        if (fx.world_aabb_m) {
          const cx = 0.5 * (fx.world_aabb_m.min_m[0] + fx.world_aabb_m.max_m[0]);
          const cy = 0.5 * (fx.world_aabb_m.min_m[1] + fx.world_aabb_m.max_m[1]);
          const cz = fx.world_aabb_m.max_m[2];
          origins.push(blenderToThree(cx, cy, Math.max(0.88, Math.min(1.22, cz))));
        }
      }
    }
    if (origins.length === 0) {
      const [w, d] = currentSpec.roomM;
      origins.push(blenderToThree(w * 0.35, d * 0.45, 0.95));
    }

    const perOrigin = 44;
    const total = origins.length * perOrigin;
    const positions = new Float32Array(total * 3);
    const colors = new Float32Array(total * 3);
    const phases = new Float32Array(total);
    const speeds = new Float32Array(total);

    let seed = 1337;
    const nextRand = () => {
      seed = (seed * 1664525 + 1013904223) >>> 0;
      return seed / 4294967296;
    };

    for (let i = 0; i < total; i += 1) {
      phases[i] = nextRand();
      speeds[i] = 0.52 + nextRand() * 0.58;
      const o = origins[Math.floor(i / perOrigin)];
      positions[i * 3] = o.x;
      positions[i * 3 + 1] = o.y;
      positions[i * 3 + 2] = o.z;
      colors[i * 3] = 0.98;
      colors[i * 3 + 1] = 0.62;
      colors[i * 3 + 2] = 0.24;
    }

    const geo = new THREE.BufferGeometry();
    geo.setAttribute("position", new THREE.BufferAttribute(positions, 3));
    geo.setAttribute("color", new THREE.BufferAttribute(colors, 3));

    const mat = new THREE.PointsMaterial({
      size: 0.46,
      map: getSoftParticleTexture(),
      vertexColors: true,
      transparent: true,
      opacity: 0.48,
      depthWrite: false,
      blending: THREE.AdditiveBlending,
    });
    const points = new THREE.Points(geo, mat);
    points.visible = activeOverlays.has("plumes");
    return {
      points,
      origins,
      phases,
      speeds,
      ceilingY: Math.min(2.25, (loaded.manifest.room_shell?.ceiling_height_m ?? 2.8) - 0.38),
    };
  }

  function buildAisleOverlays(loaded: LoadedExemplarScene): THREE.Group {
    const grp = new THREE.Group();
    grp.visible = activeOverlays.has("aisles");

    const corridors = loaded.siteMeta?.walk_mesh_corridors ?? [];
    for (const cor of corridors) {
      const w = Math.max(0.25, cor.max_m[0] - cor.min_m[0]);
      const d = Math.max(0.25, cor.max_m[1] - cor.min_m[1]);
      const cx = 0.5 * (cor.min_m[0] + cor.max_m[0]);
      const cy = 0.5 * (cor.min_m[1] + cor.max_m[1]);
      const isTight = Math.min(w, d) < 1.05;
      const col = isTight ? 0xf59e0b : 0x34d399;

      const plane = new THREE.Mesh(
        new THREE.PlaneGeometry(w, d),
        new THREE.MeshBasicMaterial({
          color: col,
          transparent: true,
          opacity: 0.22,
          side: THREE.DoubleSide,
          depthWrite: false,
        }),
      );
      plane.rotation.x = -Math.PI / 2;
      const p3 = blenderToThree(cx, cy, 0.025);
      plane.position.copy(p3);
      grp.add(plane);

      const border = new THREE.LineSegments(
        new THREE.EdgesGeometry(new THREE.BoxGeometry(w, 0.03, d)),
        new THREE.LineBasicMaterial({ color: col, transparent: true, opacity: 0.85 }),
      );
      border.position.copy(p3);
      grp.add(border);
    }

    for (const fx of loaded.manifest.fixtures) {
      for (const sw of fx.swept_volumes ?? []) {
        if (sw.kind !== "door_swing") continue;
        const mn = sw.aabb_world_m.min_m;
        const mx = sw.aabb_world_m.max_m;
        const w = Math.max(0.15, mx[0] - mn[0]);
        const d = Math.max(0.15, mx[1] - mn[1]);
        const cx = 0.5 * (mn[0] + mx[0]);
        const cy = 0.5 * (mn[1] + mx[1]);
        const p3 = blenderToThree(cx, cy, 0.045);

        const swingMesh = new THREE.Mesh(
          new THREE.PlaneGeometry(w, d),
          new THREE.MeshBasicMaterial({
            color: 0xfb7185,
            transparent: true,
            opacity: 0.24,
            side: THREE.DoubleSide,
            depthWrite: false,
          }),
        );
        swingMesh.rotation.x = -Math.PI / 2;
        swingMesh.position.copy(p3);
        grp.add(swingMesh);

        const arcRing = new THREE.LineLoop(
          new THREE.EdgesGeometry(new THREE.CircleGeometry(Math.min(w, d) * 0.48, 24)),
          new THREE.LineBasicMaterial({ color: 0xfb7185, transparent: true, opacity: 0.82 }),
        );
        arcRing.rotation.x = -Math.PI / 2;
        arcRing.position.set(p3.x, 0.055, p3.z);
        grp.add(arcRing);
      }
    }
    return grp;
  }

  function syncTelemetryUI(): void {
    const m = computeKitchenTelemetry(currentSpec, ordersPerHr, electrificationPct, ambientC);
    if (rushVal) rushVal.textContent = `${Math.round(ordersPerHr)} ord/hr`;
    if (elecVal) elecVal.textContent = `${Math.round(electrificationPct)}%`;
    if (ambVal) ambVal.textContent = `${Math.round(ambientC)} °C`;

    if (outKw) outKw.textContent = `${m.thermalKw.toFixed(1)} kW`;
    if (outCapture) outCapture.textContent = `${m.hoodCapturePct.toFixed(1)}%`;
    if (outAisle) outAisle.textContent = `${m.minAisleIn.toFixed(0)}" · ${m.minAisleM.toFixed(2)}m`;
    if (outMeta) {
      const draws = loadedScene ? `${loadedScene.drawCallCount} draws` : "loading";
      outMeta.textContent = `${currentSpec.fixtures} fixtures · ${currentSpec.staffCount} staff · ${draws}`;
    }

    if (barGas) barGas.style.flexGrow = String(Math.max(1, m.gasSharePct));
    if (barElec) barElec.style.flexGrow = String(Math.max(1, m.elecSharePct));
    if (barHvac) barHvac.style.flexGrow = String(Math.max(1, m.hvacSharePct));
    if (keyGas) keyGas.textContent = `${(m.thermalKw - m.electricalKw).toFixed(1)} kW`;
    if (keyElec) keyElec.textContent = `${m.electricalKw.toFixed(1)} kW`;
    if (keyHvac) keyHvac.textContent = `${m.hvacCoolingKw.toFixed(1)} kW`;

    applyThermalOverlay();
  }

  async function loadExemplar(id: ExemplarId): Promise<void> {
    const seq = ++loadSeq;
    const spec = EXEMPLARS.find((e) => e.id === id) ?? EXEMPLARS[0];
    currentSpec = spec;

    for (const b of archetypeBtns) {
      b.setAttribute("aria-pressed", String(b.dataset.archetype === id));
    }

    if (rushInput) {
      rushInput.min = String(Math.round(spec.nominalOrdersPerHr * 0.35));
      rushInput.max = String(spec.maxOrdersPerHr);
      ordersPerHr = spec.nominalOrdersPerHr;
      rushInput.value = String(ordersPerHr);
    }
    if (statusEl) {
      statusEl.textContent = `Loading ${spec.shortName} (.glb)…`;
    }

    const loaded = await loadExemplarScene(spec.id);
    if (seq !== loadSeq || !loaded) return;

    for (const anchor of [kitchenAnchor, detailAnchor, staffAnchor, overlayAnchor]) {
      while (anchor.children.length > 0) {
        const child = anchor.children[0];
        anchor.remove(child);
        disposeObjectTree(child, anchor === staffAnchor);
      }
    }
    walkingChefs = [];
    burnerBeacons = [];

    loadedScene = loaded;
    kitchenAnchor.add(loaded.sceneGroup);

    const [rw, rd] = spec.roomM;
    floorGrid.position.set(rw * 0.5, -0.015, -rd * 0.5);

    burnerBeacons = buildKitchenDetails(loaded);
    plumeSystem = buildPlumeSystem(loaded);
    aisleGroup = buildAisleOverlays(loaded);
    overlayAnchor.add(plumeSystem.points, aisleGroup);

    // Clone articulated 29-bone SkinnedMesh staff rigs with held culinary tools & station dwell behavior
    const chefPromises = spec.staffRoutes.map(async (route, idx) => {
      const rig = await cloneStaffRig(route.cohort);
      if (!rig || seq !== loadSeq) return null;
      const pts = route.waypoints.map(([x, y]) => blenderToThree(x, y, 0));

      // Attach a real culinary tool (sauté skillet or stainless tray) to right hand bone!
      const handR = rig.bones.get("hand_R") ?? rig.bones.get("grip_R");
      if (handR) {
        const prop = makeHeldChefProp(route.role);
        prop.position.set(0, 0.04, 0.02);
        prop.rotation.set(Math.PI * 0.42, 0, 0);
        handR.add(prop);
      }

      // Subtle station-role indicator ring under boots
      const halo = new THREE.Mesh(
        new THREE.RingGeometry(0.22, 0.27, 28),
        new THREE.MeshBasicMaterial({
          color: route.accentHex,
          side: THREE.DoubleSide,
          transparent: true,
          opacity: 0.62,
        }),
      );
      halo.rotation.x = -Math.PI / 2;
      halo.position.y = 0.018;
      rig.root.add(halo);

      const p0 = pts[0];
      const p1 = pts[1 % pts.length];
      const initYaw = Math.atan2(p1.x - p0.x, p1.z - p0.z);
      rig.root.position.copy(p0);
      rig.root.rotation.y = initYaw;
      staffAnchor.add(rig.root);

      return {
        rig,
        route,
        points: pts,
        segmentIdx: 0,
        segProgressM: 0,
        dwellLeftS: idx % 2 === 0 ? 1.8 : 0.4,
        currentYaw: initYaw,
        stridePhaseOffset: idx * 2.1,
      };
    });
    const built = await Promise.all(chefPromises);
    if (seq !== loadSeq) return;
    walkingChefs = built.filter((c): c is WalkingChef => c !== null);

    lerpToPreset(activeCameraPreset);
    syncTelemetryUI();

    if (statusEl) {
      statusEl.textContent = `${spec.shortName} · ${loaded.triangleCount.toLocaleString()} tris · ${loaded.drawCallCount} draw calls · ${walkingChefs.length} rigged staff`;
    }
    stageEl.dataset.live = "1";
  }

  for (const btn of archetypeBtns) {
    btn.addEventListener("click", () => {
      void loadExemplar(btn.dataset.archetype as ExemplarId);
    });
  }
  for (const btn of camBtns) {
    btn.addEventListener("click", () => {
      lerpToPreset(btn.dataset.cam as CameraPresetId);
    });
  }
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
      } else if (activeOverlays.has(mode)) {
        activeOverlays.delete(mode);
      } else {
        activeOverlays.add(mode);
      }
      for (const ob of overlayBtns) {
        ob.setAttribute("aria-pressed", String(activeOverlays.has(ob.dataset.overlay as OverlayMode)));
      }
      if (plumeSystem) plumeSystem.points.visible = activeOverlays.has("plumes");
      if (aisleGroup) aisleGroup.visible = activeOverlays.has("aisles");
      applyThermalOverlay();
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

  // Keyboard shortcuts: 1-5 switches canonical kitchen, Space toggles pause (when not typing in an input)
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
  let elapsed = 0;

  function animate(now: number): void {
    requestAnimationFrame(animate);
    const dt = Math.max(0, Math.min(0.05, (now - lastTs) / 1000));
    lastTs = now;
    if (!inView || document.hidden) return;

    if (camAnimT < 1) {
      camAnimT = Math.min(1, camAnimT + dt * 1.9);
      const ease = 1 - Math.pow(1 - camAnimT, 3);
      camera.position.lerpVectors(camStartPos, camGoalPos, ease);
      controls.target.lerpVectors(tgtStartPos, tgtGoalPos, ease);
    }

    if (!paused) {
      elapsed += dt;
      const loadRatio = ordersPerHr / currentSpec.nominalOrdersPerHr;
      const elecFrac = electrificationPct / 100;

      // Animate live cookline burner flames / induction coils
      for (let bi = 0; bi < burnerBeacons.length; bi += 1) {
        const b = burnerBeacons[bi];
        const ringMat = b.ringMesh.material as THREE.MeshBasicMaterial;
        const pulse = 0.72 + 0.28 * Math.sin(elapsed * 5.2 + bi * 1.3);
        const useCyan = b.isInductionCapable && elecFrac > 0.52;
        ringMat.color.setHex(useCyan ? 0x22d3ee : 0xf97316);
        ringMat.opacity = (0.35 + 0.45 * loadRatio) * pulse;
        // Hide open flames as kitchen shifts toward 100% induction!
        const flameStrength = b.isInductionCapable ? Math.max(0.05, 1.0 - elecFrac * 0.92) : 1.0;
        for (let ci = 1; ci < b.flameGroup.children.length; ci += 1) {
          const tongue = b.flameGroup.children[ci];
          tongue.visible = flameStrength > 0.14;
          if (tongue.visible) {
            const flk = 0.7 + 0.4 * Math.sin(elapsed * 18 + bi * 3.1 + ci * 1.7);
            tongue.scale.set(0.85, flk * flameStrength * (0.7 + 0.45 * loadRatio), 0.85);
          }
        }
      }

      // Animate 29-bone SkinnedMesh chefs with smooth yaw, gait biomechanics, and station work dwells!
      for (let idx = 0; idx < walkingChefs.length; idx += 1) {
        const chef = walkingChefs[idx];
        const p0 = chef.points[chef.segmentIdx];
        const p1 = chef.points[(chef.segmentIdx + 1) % chef.points.length];
        const segLen = Math.max(0.05, p0.distanceTo(p1));
        const applyDelta = (name: string, dx: number, dy = 0, dz = 0) => {
          const bone = chef.rig.bones.get(name);
          if (!bone) return;
          const rest = (bone.userData.restRot as THREE.Euler | undefined) ?? new THREE.Euler();
          bone.rotation.set(rest.x + dx, rest.y + dy, rest.z + dz);
        };

        if (chef.dwellLeftS > 0) {
          // STATION WORK POSE: chef pauses beside cookline/pass, turns slightly toward worktop,
          // and performs rhythmic sauté toss / plating / checking motions!
          chef.dwellLeftS -= dt;
          const workPhase = elapsed * 3.8 + chef.stridePhaseOffset;
          const workWave = Math.sin(workPhase);
          const targetWorkYaw = Math.atan2(p1.x - p0.x, p1.z - p0.z) + (idx % 2 === 0 ? -Math.PI * 0.45 : Math.PI * 0.45);
          chef.currentYaw += shortestAngleDelta(chef.currentYaw, targetWorkYaw) * Math.min(1, dt * 4.2);
          chef.rig.root.rotation.y = chef.currentYaw;
          chef.rig.root.position.y = Math.sin(elapsed * 2.2 + idx) * 0.005;

          // Subtle grounded stance + forward work reach & sauté wrist motion
          applyDelta("thigh_L", -0.05);
          applyDelta("thigh_R", 0.04);
          applyDelta("shin_L", 0.08);
          applyDelta("shin_R", 0.06);
          applyDelta("spine_01", 0.08);
          applyDelta("spine_02", 0.09, workWave * 0.05);
          applyDelta("head", 0.14 + workWave * 0.03, Math.sin(workPhase * 0.5) * 0.08);
          applyDelta("upper_arm_R", 0.46 + workWave * 0.16, 0, -0.08);
          applyDelta("forearm_R", 0.52 + Math.cos(workPhase) * 0.18);
          applyDelta("upper_arm_L", 0.34 - workWave * 0.08, 0, 0.08);
          applyDelta("forearm_L", 0.42 + workWave * 0.08);
        } else {
          // SMOOTH WALK CYCLE along aisle corridor
          const speed = chef.route.speedMps * (0.78 + 0.32 * loadRatio);
          chef.segProgressM += speed * dt;
          if (chef.segProgressM >= segLen) {
            chef.segProgressM = 0;
            chef.segmentIdx = (chef.segmentIdx + 1) % chef.points.length;
            chef.dwellLeftS = 2.0 + ((idx + chef.segmentIdx) % 3) * 0.85;
          }
          const a = chef.points[chef.segmentIdx];
          const b = chef.points[(chef.segmentIdx + 1) % chef.points.length];
          const u = Math.min(1, chef.segProgressM / Math.max(0.05, a.distanceTo(b)));
          chef.rig.root.position.lerpVectors(a, b, u);

          const heading = Math.atan2(b.x - a.x, b.z - a.z);
          chef.currentYaw += shortestAngleDelta(chef.currentYaw, heading) * Math.min(1, dt * 6.5);
          chef.rig.root.rotation.y = chef.currentYaw;

          const phase = elapsed * 5.2 * (0.82 + 0.22 * loadRatio) + chef.stridePhaseOffset;
          const stride = Math.sin(phase);
          // Natural double-frequency vertical COM bob
          chef.rig.root.position.y = Math.abs(Math.cos(phase)) * 0.022;

          applyDelta("pelvis", 0.02, -stride * 0.07);
          applyDelta("spine_01", 0.04, stride * 0.05);
          applyDelta("spine_02", 0.03, stride * 0.06);
          applyDelta("head", 0.03, -stride * 0.04);
          applyDelta("thigh_L", stride * 0.34);
          applyDelta("thigh_R", -stride * 0.34);
          applyDelta("shin_L", Math.max(0, -stride * 0.46) + 0.04);
          applyDelta("shin_R", Math.max(0, stride * 0.46) + 0.04);
          // Right hand holds tool steady slightly forward; left arm swings naturally
          applyDelta("upper_arm_R", 0.32 + stride * 0.10, 0, -0.05);
          applyDelta("forearm_R", 0.48 + Math.abs(stride) * 0.08);
          applyDelta("upper_arm_L", -stride * 0.26, 0, 0.05);
          applyDelta("forearm_L", 0.22 + Math.max(0, stride) * 0.18);
        }
      }

      // Animate silky curling 3D volumetric exhaust plumes rising into canopy hoods
      if (plumeSystem && plumeSystem.points.visible) {
        const posAttr = plumeSystem.points.geometry.getAttribute("position") as THREE.BufferAttribute;
        const colAttr = plumeSystem.points.geometry.getAttribute("color") as THREE.BufferAttribute;
        const n = plumeSystem.phases.length;
        const perOrigin = Math.max(1, Math.floor(n / plumeSystem.origins.length));
        const buoyancy = 0.55 + 0.55 * loadRatio * (1 - 0.38 * elecFrac);
        const mat = plumeSystem.points.material as THREE.PointsMaterial;
        mat.opacity = 0.26 + 0.32 * loadRatio * (1 - 0.44 * elecFrac);

        for (let i = 0; i < n; i += 1) {
          plumeSystem.phases[i] = (plumeSystem.phases[i] + dt * 0.34 * plumeSystem.speeds[i] * buoyancy) % 1;
          const p = plumeSystem.phases[i];
          const origin = plumeSystem.origins[Math.floor(i / perOrigin) % plumeSystem.origins.length];
          const rise = Math.max(0.45, plumeSystem.ceilingY - origin.y);
          const swirl = i * 2.39996 + elapsed * 1.25 + p * 2.4;
          // Plume expands gently mid-rise then pulls inward toward hood baffle capture zone
          const capturePinch = 1.0 - 0.32 * Math.max(0, p - 0.65) / 0.35;
          const spread = (0.10 + 0.32 * Math.sin(p * Math.PI * 0.85)) * capturePinch;
          posAttr.setXYZ(
            i,
            origin.x + Math.cos(swirl) * spread,
            origin.y + p * rise,
            origin.z + Math.sin(swirl) * spread,
          );
          // Fade softly at birth & capture ceiling via vertex brightness
          const fade = Math.sin(p * Math.PI);
          const r = ((1 - p) * (0.95 - 0.50 * elecFrac) + p * 0.22) * fade;
          const g = ((1 - p) * (0.56 + 0.26 * elecFrac) + p * 0.84) * fade;
          const b = ((1 - p) * (0.18 + 0.68 * elecFrac) + p * 0.64) * fade;
          colAttr.setXYZ(i, r, g, b);
        }
        posAttr.needsUpdate = true;
        colAttr.needsUpdate = true;
      }
    }

    controls.update();
    renderer.render(scene, camera);
  }

  void loadExemplar(initialSpec.id);
  requestAnimationFrame(animate);

  return {
    setExemplar: loadExemplar,
    getLoadedExemplarId: () => currentSpec.id,
    getDrawCalls: () => loadedScene?.drawCallCount ?? 0,
  };
}
