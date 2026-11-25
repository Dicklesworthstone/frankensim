/**
 * CHAPTER 3 - "THE CAVITY & THE CRUST": a sourdough boule baking in a 6-pan combi-steamer (#oven).
 *
 * The oven is the model library's NSF-COMBI-6PAN with its door open. A window is cut in the air baffle (a section,
 * as a service manual draws one) so the radial fan wheel and its ring element can be seen working behind it, and
 * when the view swings round to the right the right-hand wall and column come away to show them side on. The loaf
 * is ARCH-FOOD-SOURDOUGH-BOULE on a perforated GN 1/1 baking tray whose rim rides the oven's own runner rods
 * (combiCavity gives the rail heights); the oven stands on a 48 in prep table.
 *
 * Everything that moves or changes is set every frame from `stepOvenKernel` (<1e-6 parity with Rust `engine-oven`)
 * and its inputs: the fan wheel turns at the fan setting, the element glows with the cavity set point, the cavity
 * hazes with the steam setting, the loaf springs by the kernel's oven_spring_ratio and browns by its
 * maillard_progress, and the motes riding the fan's air loop are as bright as the kernel's effective h.
 */

import * as THREE from "three";
import { createInitialOvenState, stepOvenKernel, type OvenConfig, type OvenState } from "../sim/physicsKernels";
import { bindCanvasResize, createScrollSafeOrbit, createStudioRenderer, drawStripChart, getSoftParticleTexture, glOf, observeInView, shaderGate } from "./kit";
import { buildFixture, type Ext, type Model } from "./models/all";
import { palette, perforation } from "./models/materials";
import { combiCavity, vessel } from "./models/ovens";
import { Shop } from "./models/shop";

/** Archetype extents, metres (assets/kitchens/specs/equipment_dimensions.tsv). */
const COMBI_ID = "NSF-COMBI-6PAN";
const COMBI_EXT: Ext = { w: 0.847, d: 0.776, h: 0.782 };
const BOULE_ID = "ARCH-FOOD-SOURDOUGH-BOULE";
const BOULE_EXT: Ext = { w: 0.18, d: 0.18, h: 0.088 };
const TABLE_ID = "NSF-PREPTABLE-48";
const TABLE_EXT: Ext = { w: 1.219, d: 0.762, h: 0.914 };

/** The runner level the tray rides on, counted from the bottom. */
const RAIL = 1;
/** A GN 1/1 perforated baking tray: 325 x 530 mm, 20 mm deep, its flat rim 10 mm wide. */
const TRAY = { w: 0.325, d: 0.53, depth: 0.02, lip: 0.01, t: 0.0015 };

/** Simulated seconds per second of page time, and where the bake ends. */
const SIM_RATE = 15;
const BAKE_END_S = 1600;
/** The bake is already this far along when the chapter mounts, so the rise begins within moments of arriving. */
const PREWARM_S = 120;
/** The strip chart takes a point every quarter second of page time. */
const CHART_EVERY_S = 0.25;
const CHART_DT_S = CHART_EVERY_S * SIM_RATE;

/**
 * The tray: drawn aluminium walls, a flat rim that rests on the runner rods, and a perforated floor (3.6 mm holes on
 * 8 mm centres) so the fan's air reaches the loaf's base. Built in its own frame, rim underside at y = 0, centred.
 */
function buildTray(): THREE.Group {
  const s = new Shop(2, "baking_tray");
  const p = palette();
  const { w, d, depth, lip, t } = TRAY;
  const dr = 0.003;
  // Outer bottom edge up the wall, over the rim and its rolled edge, down the inside to the floor.
  const prof: [number, number][] = [
    [lip + dr + 0.006, -depth],
    [lip + dr + 0.003, -depth + 0.0025],
    [lip + 0.001, -0.004],
    [lip - 0.002, -0.0025],
    [0.0015, -0.0025],
    [0, -0.0012],
    [0.0015, 0],
    [lip - 0.002, 0],
    [lip + t, -0.0015],
    [lip + dr + t + 0.003, -depth + t + 0.0025],
    [lip + dr + t + 0.006, -depth + t],
  ];
  s.add(p.alu, vessel(s, w, d, 0.014, prof, [false, false]), [0, 0.0025, 0]);
  // The floor: one perforated sheet across the wall's foot.
  const perf = new THREE.MeshStandardMaterial({ color: 0xc9cdd1, metalness: 0.9, roughness: 0.44, alphaMap: perforation(32, 0.45), alphaTest: 0.5, side: THREE.DoubleSide });
  perf.alphaToCoverage = true;
  perf.userData.uv = "box";
  perf.userData.uvScale = 1 / (0.008 * 32);
  const fw = w - 2 * (lip + dr + 0.006);
  const fd = d - 2 * (lip + dr + 0.006);
  // The wall's foot is square-cornered (the rim's 14 mm radius is used up by the wall's inset); a hair of rounding.
  const sh = new THREE.Shape();
  const r = 0.002;
  sh.moveTo(-fw / 2 + r, -fd / 2);
  sh.lineTo(fw / 2 - r, -fd / 2);
  sh.quadraticCurveTo(fw / 2, -fd / 2, fw / 2, -fd / 2 + r);
  sh.lineTo(fw / 2, fd / 2 - r);
  sh.quadraticCurveTo(fw / 2, fd / 2, fw / 2 - r, fd / 2);
  sh.lineTo(-fw / 2 + r, fd / 2);
  sh.quadraticCurveTo(-fw / 2, fd / 2, -fw / 2, fd / 2 - r);
  sh.lineTo(-fw / 2, -fd / 2 + r);
  sh.quadraticCurveTo(-fw / 2, -fd / 2, -fw / 2 + r, -fd / 2);
  const floor = new THREE.ShapeGeometry(sh, 4);
  floor.rotateX(-Math.PI / 2);
  s.add(perf, floor, [0, 0.0025 - depth + t * 0.5, 0]);
  return s.build(true).group;
}

/**
 * A dark studio for the stainless to reflect: a big warm softbox high in front, a warm strip to the right, a cool
 * strip behind on the left, the rest nearly black, so that steel reads as steel (bright bands on dark) rather than
 * as grey paint.
 */
function darkStudio(r: THREE.WebGLRenderer): THREE.Texture {
  const env = new THREE.Scene();
  env.background = new THREE.Color(0x0c0907);
  const panel = (w: number, h: number, color: number, k: number, at: [number, number, number]) => {
    const m = new THREE.Mesh(new THREE.PlaneGeometry(w, h), new THREE.MeshBasicMaterial({ color: new THREE.Color(color).multiplyScalar(k), side: THREE.DoubleSide }));
    m.position.set(at[0], at[1], at[2]);
    m.lookAt(0, 0, 0);
    env.add(m);
  };
  panel(7, 4, 0xffd6a6, 2.2, [0.5, 5.5, 5.5]);
  panel(1.4, 7, 0xffc890, 2.6, [7, 1.5, 2.5]);
  panel(1.2, 6, 0x9fc0ff, 1.1, [-6.5, 2, -4]);
  panel(9, 1.2, 0xffb070, 0.9, [0, 1.2, -8]);
  panel(14, 14, 0x2a1d14, 1, [0, -6, 0]);
  const pm = new THREE.PMREMGenerator(glOf(r));
  const t = pm.fromScene(env, 0.02).texture;
  pm.dispose();
  return t;
}

export function mountOvenLab(root: HTMLElement): void {
  const stageEl = root.querySelector<HTMLElement>("[data-oven-stage]")!;
  const canvas = root.querySelector<HTMLCanvasElement>("[data-oven-canvas]")!;
  const chartCanvas = root.querySelector<HTMLCanvasElement>("[data-oven-chart]");

  const maybeRenderer = createStudioRenderer(canvas, true);
  if (!maybeRenderer) return;
  const renderer = maybeRenderer;
  renderer.toneMappingExposure = 1.0;
  renderer.localClippingEnabled = true;
  renderer.setClearColor(0x000000, 0);

  const t0 = performance.now();
  const scene = new THREE.Scene();
  scene.environment = darkStudio(renderer);
  scene.environmentIntensity = 1.0;

  /* ---------------------------------------------------------------------------------------------------------- */
  /* The machine                                                                                                 */
  /* ---------------------------------------------------------------------------------------------------------- */

  const cav = combiCavity(COMBI_ID, COMBI_EXT)!;
  const oven: Model = buildFixture(COMBI_ID, COMBI_EXT, 2);
  const tOven = performance.now();
  // World frame: the cavity's centre line on x = 0, the tray's centre (pushed home) on z = 0, the oven on the table.
  const O = new THREE.Vector3(-cav.cx, TABLE_EXT.h, -cav.zc);
  oven.group.position.copy(O);
  scene.add(oven.group);
  oven.hooks.open?.(0.8);
  // The cavity's steam is fog lit dimly by the lamp and the element: a warm grey, not a glow.
  oven.group.traverse((o) => {
    if ((o as THREE.Sprite).isSprite) ((o as THREE.Sprite).material as THREE.SpriteMaterial).color.set(0x9c8670);
  });

  const table = buildFixture(TABLE_ID, TABLE_EXT, 1);
  table.group.position.set(O.x + COMBI_EXT.w / 2 - TABLE_EXT.w / 2, 0, O.z + 0.035);
  scene.add(table.group);

  // The section window in the air baffle: from just left of and below the hub to near the baffle's top right corner,
  // through the baffle's thickness only. Cut on this chapter's own copies of the baffle's materials.
  const fanJ = oven.joints.get("fan");
  const fanY = fanJ ? fanJ.position.y : (cav.box.y0 + cav.box.y1) / 2;
  const win = {
    x0: O.x + cav.cx - 0.045,
    x1: O.x + cav.box.x1 - 0.055,
    y0: O.y + fanY - 0.045,
    y1: O.y + cav.box.y1 - 0.055,
    z0: O.z + cav.box.z0 - 0.0026,
    z1: O.z + cav.box.z0 + 0.0016,
  };
  const winPlanes = [
    new THREE.Plane(new THREE.Vector3(-1, 0, 0), win.x0),
    new THREE.Plane(new THREE.Vector3(1, 0, 0), -win.x1),
    new THREE.Plane(new THREE.Vector3(0, -1, 0), win.y0),
    new THREE.Plane(new THREE.Vector3(0, 1, 0), -win.y1),
    new THREE.Plane(new THREE.Vector3(0, 0, -1), win.z0),
    new THREE.Plane(new THREE.Vector3(0, 0, 1), -win.z1),
  ];
  const pal = palette();
  for (const o of oven.group.children) {
    const m = o as THREE.Mesh;
    if (!m.isMesh) continue;
    const mat = m.material as THREE.MeshStandardMaterial;
    if (mat !== pal.ss && !mat.alphaMap) continue;
    const c = mat.clone();
    c.clippingPlanes = winPlanes;
    c.clipIntersection = true;
    c.clipShadows = true;
    m.material = c;
  }

  /* ---------------------------------------------------------------------------------------------------------- */
  /* The tray and the loaf                                                                                       */
  /* ---------------------------------------------------------------------------------------------------------- */

  const railY = cav.rails[RAIL]!;
  const tray = buildTray();
  tray.position.set(O.x + cav.cx, O.y + railY, O.z + cav.zc);
  scene.add(tray);
  const trayFloorY = O.y + railY + 0.0025 - TRAY.depth + TRAY.t;

  // The loaf is the costliest part to build (its crust is painted texel by texel), so it is made in a task of its own
  // just after the mount, rather than in the same stall as the oven.
  const bouleHolder = new THREE.Group();
  // Forward of the tray's middle, its score toward the door.
  const LOAF_Z = 0.075;
  bouleHolder.position.set(0.005, trayFloorY, LOAF_Z);
  bouleHolder.rotation.y = 0.2;
  scene.add(bouleHolder);
  let boule: Model | null = null;
  let bouleMs = 0;
  function buildBoule(): void {
    const tb = performance.now();
    boule = buildFixture(BOULE_ID, BOULE_EXT, 2);
    // The boule's frame has its joint at (w/2, 0, -d/2); centre it on the holder, which sits on the tray floor.
    boule.group.position.set(-BOULE_EXT.w / 2, 0, BOULE_EXT.d / 2);
    bouleHolder.add(boule.group);
    bouleMs = performance.now() - tb;
  }
  const loafCentre = new THREE.Vector3(0, trayFloorY + 0.035, LOAF_Z);

  /* ---------------------------------------------------------------------------------------------------------- */
  /* Light                                                                                                       */
  /* ---------------------------------------------------------------------------------------------------------- */

  const hemi = new THREE.HemisphereLight(0xffe9d2, 0x2a2018, 0.35);
  scene.add(hemi);
  // Warm key, high in front and a little right: into the cavity through the open door, onto the loaf.
  const key = new THREE.DirectionalLight(0xffe2c0, 2.4);
  key.position.copy(loafCentre).add(new THREE.Vector3(0.55, 1.45, 1.6));
  key.target.position.copy(loafCentre);
  key.castShadow = true;
  key.shadow.mapSize.set(2048, 2048);
  {
    const sc = key.shadow.camera;
    sc.left = -0.75;
    sc.right = 0.75;
    sc.top = 0.75;
    sc.bottom = -0.75;
    sc.near = 0.5;
    sc.far = 4.5;
  }
  key.shadow.bias = -0.0002;
  key.shadow.normalBias = 0.004;
  key.shadow.radius = 4;
  scene.add(key, key.target);
  // Cool rim from behind left, to draw the housing's edges out of the dark.
  const rim = new THREE.DirectionalLight(0xa8c4ff, 0.9);
  rim.position.copy(loafCentre).add(new THREE.Vector3(-1.6, 0.9, -1.4));
  scene.add(rim);
  // The cavity lamp: a halogen behind glass on the left wall, near the top front.
  const lampPos = new THREE.Vector3(O.x + cav.box.x0 + 0.012, O.y + cav.box.y1 - 0.075, O.z + cav.box.z1 - 0.13);
  const lamp = new THREE.SpotLight(0xffb468, 0.9, 1.4, 1.05, 0.7, 2);
  lamp.position.copy(lampPos);
  lamp.target.position.copy(loafCentre);
  scene.add(lamp, lamp.target);
  // The element's glow, behind the baffle round the fan: it lights what faces it (the wheel, the baffle's back, the
  // loaf's back through the eye and window), and nothing that faces away.
  const glowLight = new THREE.PointLight(0xff6a24, 0, 0.24, 2);
  glowLight.position.set(O.x + cav.cx, O.y + fanY, O.z + (fanJ ? fanJ.position.z : cav.box.z0 - 0.034) + 0.012);
  scene.add(glowLight);

  /* ---------------------------------------------------------------------------------------------------------- */
  /* The air                                                                                                     */
  /* ---------------------------------------------------------------------------------------------------------- */

  // The fan draws cavity air in through the baffle's eye and throws it out through the slots down the baffle's sides;
  // it runs forward along each wall, turns in at the door and comes back down the middle, over the loaf, to the eye.
  // Motes ride that loop, each side its own way round, slowed for the eye in proportion to the fan, and as bright as
  // the kernel's effective h (how hard the air is working on the crust).
  const AIR_N = 150;
  const air = {
    x0: O.x + cav.box.x0,
    x1: O.x + cav.box.x1,
    y0: O.y + cav.box.y0 + 0.035,
    y1: O.y + cav.box.y1 - 0.045,
    zb: O.z + cav.box.z0 + 0.03,
    zf: O.z + cav.box.z1 - 0.07,
    cx: O.x + cav.cx,
    eyeY: O.y + fanY,
  };
  const motes: { s: number; l: number; y: number; ph: number; sp: number; tw: number }[] = [];
  {
    let seed = 8821;
    const rnd = () => {
      seed = (seed * 1664525 + 1013904223) >>> 0;
      return seed / 4294967296;
    };
    for (let i = 0; i < AIR_N; i += 1) motes.push({ s: i % 2 ? 1 : -1, l: rnd(), y: air.y0 + rnd() * (air.y1 - air.y0), ph: rnd() * Math.PI * 2, sp: 0.8 + rnd() * 0.4, tw: rnd() * 6.28 });
  }
  const airPos = new Float32Array(AIR_N * 3);
  const airCol = new Float32Array(AIR_N * 3);
  const airGeo = new THREE.BufferGeometry();
  airGeo.setAttribute("position", new THREE.BufferAttribute(airPos, 3).setUsage(THREE.DynamicDrawUsage));
  airGeo.setAttribute("color", new THREE.BufferAttribute(airCol, 3).setUsage(THREE.DynamicDrawUsage));
  const airMat = new THREE.PointsMaterial({ size: 0.014, map: getSoftParticleTexture(), vertexColors: true, transparent: true, depthWrite: false, blending: THREE.AdditiveBlending });
  const airPts = new THREE.Points(airGeo, airMat);
  airPts.frustumCulled = false;
  scene.add(airPts);
  const airTint = new THREE.Color(0xffc68c);
  const BACK = Math.PI * 1.5;
  let airPhase = 0;
  function moveAir(t: number, dt: number, fan: number, hEff: number): void {
    // Slow-motion air speed along the loop, m/s: a still cavity's natural draught plus the fan's.
    const v = 0.05 + 0.32 * fan;
    const hw = (air.x1 - air.x0) / 2;
    const zc = (air.zb + air.zf) / 2;
    const az = (air.zf - air.zb) / 2;
    airPhase += dt * v;
    const k = Math.min(1, hEff / 90);
    for (let i = 0; i < AIR_N; i += 1) {
      const m = motes[i]!;
      const xw = m.s * (hw - 0.035 - 0.04 * m.l);
      const xm = m.s * (0.03 + 0.07 * m.l);
      const xc = (xw + xm) / 2;
      const ax = (xw - xm) / 2;
      const per = Math.PI * (Math.abs(ax) + az);
      const phi = (m.ph + (airPhase * m.sp * Math.PI * 2) / per) % (Math.PI * 2);
      // Gone behind the baffle at the back, between the eye and the slots.
      let d = Math.abs(phi - BACK);
      d = Math.min(d, Math.PI * 2 - d);
      const seen = THREE.MathUtils.smoothstep(d, 0.35, 0.9);
      const toEye = Math.exp(-((d / 0.8) ** 2)) * 0.6;
      airPos[i * 3] = air.cx + xc + ax * Math.cos(phi);
      airPos[i * 3 + 1] = m.y + (air.eyeY - m.y) * toEye;
      airPos[i * 3 + 2] = zc + az * Math.sin(phi);
      const a = 0.42 * k * seen * (0.55 + 0.45 * Math.sin(t * 0.9 + m.tw));
      airCol[i * 3] = airTint.r * a;
      airCol[i * 3 + 1] = airTint.g * a;
      airCol[i * 3 + 2] = airTint.b * a;
    }
    airGeo.attributes.position!.needsUpdate = true;
    airGeo.attributes.color!.needsUpdate = true;
  }

  /* ---------------------------------------------------------------------------------------------------------- */
  /* Camera                                                                                                      */
  /* ---------------------------------------------------------------------------------------------------------- */

  const BASE_FOV = 36;
  const BASE_ASPECT = 1.46;
  const camera = new THREE.PerspectiveCamera(BASE_FOV, 16 / 10, 0.02, 12);
  const target = loafCentre.clone().add(new THREE.Vector3(0.015, 0.06, -0.175));
  camera.position.copy(target).add(new THREE.Vector3(0.075, 0.15, 0.68));
  const controls = createScrollSafeOrbit(camera, canvas, [target.x, target.y, target.z], {
    minDistance: 0.22,
    maxDistance: 2.6,
    maxPolarAngle: Math.PI * 0.6,
  });
  controls.minAzimuthAngle = -Math.PI * 0.6;
  controls.maxAzimuthAngle = Math.PI * 0.6;
  controls.update();
  const gate = shaderGate(renderer, camera);
  bindCanvasResize(stageEl, renderer, camera, gate);
  // Swing round to the right and the right-hand wall and column come away, as a service manual's cutaway does, to
  // show the wheel and element side on behind the baffle; from the front the wall stays, reflecting the lamp.
  let cutOn = false;
  oven.hooks.set?.("cutaway", 0);
  function followCutaway(): void {
    const az = Math.atan2(camera.position.x - controls.target.x, camera.position.z - controls.target.z);
    const want = cutOn ? az > 0.6 : az > 0.75;
    if (want === cutOn) return;
    cutOn = want;
    oven.hooks.set?.("cutaway", cutOn ? 1 : 0);
  }
  let fitAspect = -1;
  /** Keep the horizontal field when the stage is narrower than the composition (a phone held upright). */
  function fitFov(): void {
    if (Math.abs(camera.aspect - fitAspect) < 1e-3) return;
    fitAspect = camera.aspect;
    const h = 2 * Math.atan(Math.tan(THREE.MathUtils.degToRad(BASE_FOV) / 2) * BASE_ASPECT);
    const v = camera.aspect >= BASE_ASPECT ? THREE.MathUtils.degToRad(BASE_FOV) : 2 * Math.atan(Math.tan(h / 2) / camera.aspect);
    camera.fov = THREE.MathUtils.radToDeg(v);
    camera.updateProjectionMatrix();
  }

  /* ---------------------------------------------------------------------------------------------------------- */
  /* The kernel                                                                                                  */
  /* ---------------------------------------------------------------------------------------------------------- */

  const cfg: OvenConfig = {
    cavityTempC: 235,
    steamRh: 0.65,
    fanSpeed: 0.55,
    hydration: 0.76,
    initialMassG: 850,
  };
  let state: OvenState = createInitialOvenState(cfg);
  const histCrust: number[] = [];
  const histCore: number[] = [];
  const HIST_MAX = 90;
  function record(): void {
    histCrust.push(state.crust_temp_c);
    histCore.push(state.core_temp_c);
    if (histCrust.length > HIST_MAX) {
      histCrust.shift();
      histCore.shift();
    }
  }
  function startBake(prewarmS: number): void {
    state = createInitialOvenState(cfg);
    histCrust.length = 0;
    histCore.length = 0;
    record();
    // The same spacing as the live chart's points, so the strip's time axis is uniform.
    const steps = Math.round(prewarmS / CHART_DT_S);
    for (let s = 0; s < steps; s += 1) {
      stepOvenKernel(state, cfg, CHART_DT_S);
      record();
    }
  }
  startBake(PREWARM_S);

  // DOM bindings
  const tempIn = root.querySelector<HTMLInputElement>("[data-oven-temp]");
  const steamIn = root.querySelector<HTMLInputElement>("[data-oven-steam]");
  const fanIn = root.querySelector<HTMLInputElement>("[data-oven-fan]");
  const hydIn = root.querySelector<HTMLInputElement>("[data-oven-hyd]");
  const resetBtn = root.querySelector<HTMLButtonElement>("[data-oven-reset]");

  const tempVal = root.querySelector<HTMLElement>("[data-oven-temp-val]");
  const steamVal = root.querySelector<HTMLElement>("[data-oven-steam-val]");
  const fanVal = root.querySelector<HTMLElement>("[data-oven-fan-val]");
  const hydVal = root.querySelector<HTMLElement>("[data-oven-hyd-val]");

  const outSpring = root.querySelector<HTMLElement>("[data-ov-spring]");
  const outCrust = root.querySelector<HTMLElement>("[data-ov-crust]");
  const outCore = root.querySelector<HTMLElement>("[data-ov-core]");
  const outH = root.querySelector<HTMLElement>("[data-ov-h]");
  const outMaillard = root.querySelector<HTMLElement>("[data-ov-maillard]");
  const outYield = root.querySelector<HTMLElement>("[data-ov-yield]");

  // Take the controls' starting values from the page, so the markup and the kernel agree.
  if (tempIn) cfg.cavityTempC = Number(tempIn.value);
  if (steamIn) cfg.steamRh = Number(steamIn.value) / 100;
  if (fanIn) cfg.fanSpeed = Number(fanIn.value) / 100;
  if (hydIn) cfg.hydration = Number(hydIn.value) / 100;

  function syncOvenUI(): void {
    if (tempVal) tempVal.textContent = `${Math.round(cfg.cavityTempC)} °C`;
    if (steamVal) steamVal.textContent = `${Math.round(cfg.steamRh * 100)}% RH`;
    if (fanVal) fanVal.textContent = `${Math.round(cfg.fanSpeed * 100)}%`;
    if (hydVal) hydVal.textContent = `${Math.round(cfg.hydration * 100)}%`;

    if (outSpring) outSpring.textContent = `${state.oven_spring_ratio.toFixed(2)}×`;
    if (outCrust) outCrust.textContent = `${state.crust_temp_c.toFixed(1)} °C`;
    if (outCore) outCore.textContent = `${state.core_temp_c.toFixed(1)} °C`;
    if (outH) outH.textContent = `${state.effective_h_w_m2_k.toFixed(1)} W/m²K`;
    if (outMaillard) outMaillard.textContent = `${(state.maillard_progress * 100).toFixed(1)}%`;
    if (outYield) outYield.textContent = `${(state.cook_yield_ppm / 10000).toFixed(1)}%`;

    if (chartCanvas) {
      drawStripChart(
        chartCanvas,
        [
          { label: "Crust °C", color: "#eab308", values: histCrust },
          { label: "Core °C", color: "#a5f3fc", values: histCore },
        ],
        {
          minY: 15,
          maxY: 265,
          unit: "°C",
          threshold: { value: 96, color: "#34d399", label: "Crumb set 96°C" },
        },
      );
    }
  }

  tempIn?.addEventListener("input", () => {
    cfg.cavityTempC = Number(tempIn.value);
    syncOvenUI();
  });
  steamIn?.addEventListener("input", () => {
    cfg.steamRh = Number(steamIn.value) / 100;
    syncOvenUI();
  });
  fanIn?.addEventListener("input", () => {
    cfg.fanSpeed = Number(fanIn.value) / 100;
    syncOvenUI();
  });
  hydIn?.addEventListener("input", () => {
    cfg.hydration = Number(hydIn.value) / 100;
    syncOvenUI();
  });
  resetBtn?.addEventListener("click", () => {
    startBake(0);
    syncOvenUI();
  });

  syncOvenUI();
  let inView = true;
  observeInView(stageEl, (visible) => {
    inView = visible;
  });

  /* ---------------------------------------------------------------------------------------------------------- */
  /* The frame: the kernel's state onto the models                                                               */
  /* ---------------------------------------------------------------------------------------------------------- */

  const clamp01 = (v: number) => Math.max(0, Math.min(1, v));
  function showState(): void {
    const oh = oven.hooks;
    oh.set?.("fan", cfg.fanSpeed);
    // The element's glow from the set point over the kernel's working range (95 to 280 °C): a sheathed element holding
    // a cavity at these temperatures runs from a dull red to orange, never to yellow-white.
    const heat = 0.22 + 0.45 * clamp01((cfg.cavityTempC - 95) / (280 - 95));
    oh.set?.("heat", heat);
    // Steam that is well superheated is a clear gas; it shows as haze nearer saturation. So the haze is the steam
    // fraction thinned by the cavity's superheat above 100 °C.
    const superheat = Math.max(0.25, Math.min(1, 1 - (cfg.cavityTempC - 100) / 230));
    oh.set?.("steam", 0.7 * cfg.steamRh * superheat);
    glowLight.intensity = heat <= 0.2 ? 0 : 0.35 * Math.pow((heat - 0.2) / 0.8, 1.6);
    // The loaf reshapes and repaints its crust (a million texels) when told of a change, so it is told in steps finer
    // than the eye can follow: 0.004 of spring ratio, 0.012 of Maillard progress. The readouts carry the exact values.
    const bh = boule?.hooks;
    bh?.set?.("spring", Math.round(state.oven_spring_ratio / 0.004) * 0.004);
    bh?.set?.("crust", Math.round(state.maillard_progress / 0.012) * 0.012);
  }

  let chartAccum = 0;
  let lastTs = performance.now();
  let tAnim = 0;
  showState();

  function animate(now: number): void {
    requestAnimationFrame(animate);
    const dt = Math.max(0, Math.min(0.05, (now - lastTs) / 1000));
    lastTs = now;
    if (!inView || document.hidden) return;
    tAnim += dt;

    if (state.time_s < BAKE_END_S) stepOvenKernel(state, cfg, Math.min(dt * SIM_RATE, BAKE_END_S - state.time_s + 1e-6));
    chartAccum += dt;
    if (chartAccum >= CHART_EVERY_S) {
      chartAccum = 0;
      // Once the bake is done the chart holds its last stretch.
      if (state.time_s < BAKE_END_S) record();
      syncOvenUI();
    }

    showState();
    oven.hooks.tick?.(tAnim, dt, cfg.fanSpeed);
    moveAir(tAnim, dt, cfg.fanSpeed, state.effective_h_w_m2_k);
    fitFov();
    controls.update();
    followCutaway();
    // The bake runs on; the 3D view waits for its shaders.
    if (gate.ready()) renderer.render(scene, camera);
  }
  requestAnimationFrame(animate);
  const tMount = performance.now();
  void gate.refresh(scene);
  // The chapter mounts well before it scrolls into view; the loaf is ready by the time it gets there.
  setTimeout(() => {
    if (!boule) buildBoule();
    showState();
    void gate.refresh(scene);
  }, 0);

  if (import.meta.env.DEV) {
    const w = window as unknown as Record<string, unknown>;
    w.__ovenLab = {
      build: () => ({ ovenMs: tOven - t0, mountMs: tMount - t0, bouleMs, ovenTris: oven.tris, bouleTris: boule?.tris, tableTris: table.tris }),
      render: () => ({ ...renderer.info.render, calls: renderer.info.render.calls }),
      state: () => ({ ...state }),
      /** Run the bake ahead by `s` simulated seconds. */
      jump: (s: number) => {
        const n = Math.max(1, Math.round(s / 10));
        for (let i = 0; i < n; i += 1) {
          stepOvenKernel(state, cfg, Math.min(10, Math.max(0, BAKE_END_S - state.time_s)) || 1e-9);
          record();
        }
        syncOvenUI();
      },
      stage: () => [stageEl.clientWidth, stageEl.clientHeight],
      /** Look from `eye` at `at` (world metres). */
      view: (eye: [number, number, number], at: [number, number, number]) => {
        controls.target.set(at[0], at[1], at[2]);
        camera.position.set(eye[0], eye[1], eye[2]);
        controls.update();
      },
      fov: (f: number) => {
        camera.fov = f;
        camera.updateProjectionMatrix();
      },
      bench: () => {
        const bh = boule!.hooks;
        const t1 = performance.now();
        bh.set?.("crust", Math.min(1, state.maillard_progress + 0.05));
        const t2 = performance.now();
        bh.set?.("spring", Math.min(1.85, state.oven_spring_ratio + 0.05));
        const t3 = performance.now();
        return { crustRepaintMs: t2 - t1, springReshapeAndRepaintMs: t3 - t2 };
      },
      yaw: (r: number) => {
        bouleHolder.rotation.y = r;
      },
      cam: () => ({ eye: camera.position.toArray(), at: controls.target.toArray(), fov: camera.fov }),
    };
  }
}
