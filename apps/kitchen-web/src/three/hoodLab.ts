/**
 * Keep the hood (chapter 06, after the floor plan): three gas broilers under a 10-foot wall canopy hood that stays
 * where it is. The broilers stand flush with the hood's front edge, 6 inches back or 12 inches back, the gap between
 * their backs and the wall open or closed with a steel sheet; the hood's fan pulls what the reader sets. Whether all
 * the smoke goes up the hood is the laboratory's measured capture rate for that setup against the fan (sim/
 * hoodCapture.ts); the smoke and the air it draws in are drawn to show it.
 *
 * The drawing: smoke rises off the grates and leans as the air feeding it pushes it. With the gap open, air comes up
 * the gap behind the broilers and leans the smoke forward, toward the hood's front edge; sealed, the air comes from
 * the front and leans it back. When the fan pulls less than the setup needs, part of the smoke rolls out under the
 * hood's front edge and up into the kitchen, more of it the further short the fan is; else all of it is drawn up
 * into the hood and back to its filters.
 */
import * as THREE from "three";
import { CODE_CFM, COLD_OUT_F, EDGE_IN, SETS, SETUPS, catches, caughtCount, needCfm, perFoot, warmBtuh, type SetIndex, type Setup } from "../sim/hoodCapture";
import { bindCanvasResize, clipFitter, createScrollSafeOrbit, createStudioRenderer, observeInView, shaderGate } from "./kit";
import { buildFixture, type Piece } from "./models/all";
import { palette } from "./models/materials";
import { makeStaff } from "./models/people";
import { Shop, rbox } from "./models/shop";
import { studioEnvironment, studioLights } from "./models/stage";
import { STOVE_STEAK } from "./stoveLab";
import { backdrop, tilePlane } from "./weekLab";

const IN = 0.0254;
const BROILER = { id: "NSF-CHARBROILER-36", w: 0.914, d: 0.813, h: 0.381 } as const;
/** The hood: 10 ft long, as deep as the broilers plus the widest gap behind them, its lower edge 78 in up. */
const HOOD = { id: "NSF-HOOD-TYPE1-CANOPY", w: 10 * 12 * IN, d: BROILER.d + SETS[0].gapIn * IN, h: 0.61 } as const;
/** Each broiler on a stainless equipment stand, its grates at about 36 inches: a work table (models/support.ts
 *  workTable, built 36 inches tall whatever its extents say) brought down to 23 inches. */
const STAND = { id: "NSF-PREPTABLE-48", w: 0.914, d: 0.79, h: 0.584 } as const;
const TABLE_Y = 0.914;
const STEAK = { w: STOVE_STEAK.w, d: STOVE_STEAK.d, h: STOVE_STEAK.h } as const;
const EDGE_Y = EDGE_IN * IN;
/** In the hood's frame (models/support.ts typeOneHood): the make-up air plenum along its front, 0.2 m deep, and the
 *  filter bank sloping up from 0.1 m above its edge, its foot 0.08 m off the wall. Smoke is caught behind the plenum. */
const PLENUM = 0.2;
const CAPTURE_Z = HOOD.d - PLENUM;
const FILTER_Z = 0.36;
const HOOD_TOP_IN = EDGE_Y + HOOD.h - 0.12;
/** The broilers' grates, in their own frame (models/cooking.ts charbroiler): their top at the front, how much higher
 *  at the back, their middle, and front to back. */
const GRATE_YF = 0.318;
const GRATE_RISE = 0.018;
const GRATE_Y = GRATE_YF + GRATE_RISE / 2;
const GRATE_Z0 = -0.072;
const GRATE_Z1 = -0.698;
const COOK_Y = STAND.h + GRATE_Y;
/** The broilers' middles along the hood, and how far their line runs each way. */
const XS = [-BROILER.w, 0, BROILER.w] as const;
const LINE_HALF = 1.5 * BROILER.w;
const HOOD_HALF = HOOD.w / 2;
/** Where the broilers' fronts stand for each setting. */
const frontOf = (at: SetIndex) => HOOD.d - SETS[at].backIn * IN;

/** The colours: the air the smoke draws in, a setup the fan catches and one it does not, the fan. */
const AIR = "#7dd3fc";
const OK = "#4ade80";
const SPILL = "#fb7185";
const FAN = "#fbbf24";

/** The broilers close up on a laptop, a step coarser on a phone: chosen once, so what main.ts builds ahead is used. */
let detailNow: 1 | 2 | null = null;
const hoodDetail = (): 1 | 2 => (detailNow ??= Math.min(window.innerWidth, window.innerHeight) < 600 ? 1 : 2);
export function hoodPieces(): readonly Piece[] {
  const d = hoodDetail();
  return [
    [HOOD.id, HOOD, d],
    [BROILER.id, BROILER, d],
    [STAND.id, STAND, d],
    ["ARCH-FOOD-STEAK-RIBEYE", STEAK, STOVE_STEAK.detail],
  ];
}

const cfm = (v: number) => `${Math.round(v).toLocaleString("en-US")} cfm`;
const btuh = (v: number) => `${(Math.round(v / 1000) * 1000).toLocaleString("en-US")} Btu/h`;
const smooth = (x: number) => {
  const t = Math.min(1, Math.max(0, x));
  return t * t * (3 - 2 * t);
};
/** Where the broilers stand, as a sentence names it. */
const SAID_AT = ["flush with the hood's front edge", "6 inches back from the hood's edge", "12 inches back from the hood's edge"] as const;
const SHORT_AT = ["Flush with the edge", "6 in back", "12 in back"] as const;
/** The share of the smoke drawn escaping when the fan pulls `fan` against a need of `need`: none at or over it. */
const spillShare = (fan: number, need: number) => (fan >= need ? 0 : Math.min(0.62, 0.16 + 1.6 * (1 - fan / need)));
/** The tour: flush and open, back 6, back 12, sealed, then sealed and flush. */
const TOUR: readonly Setup[] = [
  { at: 0, sealed: false },
  { at: 1, sealed: false },
  { at: 2, sealed: false },
  { at: 2, sealed: true },
  { at: 0, sealed: true },
];
const TOUR_S = 5.5;

/** A soft puff: a round fall-off broken up by a little value noise (its red channel kept for shading). */
function puffTexture(size: number, hard: boolean): THREE.CanvasTexture {
  const c = document.createElement("canvas");
  c.width = c.height = size;
  const g = c.getContext("2d")!;
  const img = g.createImageData(size, size);
  const n = 8;
  const grid = Array.from({ length: (n + 1) * (n + 1) }, (_, i) => 0.5 + 0.5 * Math.sin(i * 12.9898 + 78.233 * Math.cos(i * 4.1414)));
  const noise = (u: number, v: number) => {
    const x = u * n;
    const y = v * n;
    const i = Math.min(n - 1, Math.floor(x));
    const j = Math.min(n - 1, Math.floor(y));
    const fx = smooth(x - i);
    const fy = smooth(y - j);
    const a = grid[j * (n + 1) + i]!;
    const b = grid[j * (n + 1) + i + 1]!;
    const c2 = grid[(j + 1) * (n + 1) + i]!;
    const d = grid[(j + 1) * (n + 1) + i + 1]!;
    return a + (b - a) * fx + (c2 - a) * fy + (a - b - c2 + d) * fx * fy;
  };
  for (let y = 0; y < size; y++)
    for (let x = 0; x < size; x++) {
      const u = (x + 0.5) / size;
      const v = (y + 0.5) / size;
      const r = Math.hypot(u - 0.5, v - 0.5) * 2;
      const nz = noise(u, v);
      const base = hard ? 1 - smooth((r - 0.45) / 0.55) : (1 - smooth(r)) * (0.55 + 0.45 * nz);
      const k = (y * size + x) * 4;
      img.data[k] = Math.round(255 * (0.6 + 0.4 * nz));
      img.data[k + 1] = 255;
      img.data[k + 2] = 255;
      img.data[k + 3] = Math.round(255 * Math.max(0, Math.min(1, base)));
    }
  g.putImageData(img, 0, 0);
  const t = new THREE.CanvasTexture(c);
  t.colorSpace = THREE.NoColorSpace;
  return t;
}

/** Charcoal glazed subway tile, 6 by 3 inches in running bond on pale grout, for a wall `w` by `h` (m). */
function tileWall(w: number, h: number): THREE.MeshStandardMaterial {
  const c = document.createElement("canvas");
  c.width = c.height = 512;
  const g = c.getContext("2d")!;
  // 0.608 m a side: four tiles across, eight courses up.
  g.fillStyle = "#4c4f52";
  g.fillRect(0, 0, 512, 512);
  const tw = 128;
  const th = 64;
  let seed = 7;
  const r = () => ((seed = (seed * 16807) % 2147483647) / 2147483647);
  for (let row = 0; row < 8; row++)
    for (let col = -1; col < 5; col++) {
      const x = col * tw + (row % 2 ? tw / 2 : 0);
      const v = 30 + Math.round(r() * 7);
      g.fillStyle = `rgb(${v},${v + 3},${v + 6})`;
      g.fillRect(x + 1.5, row * th + 1.5, tw - 3, th - 3);
      // The glaze pooled a little lighter at the tile's edges.
      g.strokeStyle = `rgba(255,255,255,0.05)`;
      g.lineWidth = 2;
      g.strokeRect(x + 4, row * th + 4, tw - 8, th - 8);
    }
  const t = new THREE.CanvasTexture(c);
  t.colorSpace = THREE.SRGBColorSpace;
  t.wrapS = t.wrapT = THREE.RepeatWrapping;
  t.repeat.set(w / 0.608, h / 0.608);
  t.anisotropy = 4;
  return new THREE.MeshStandardMaterial({ map: t, roughness: 0.22, metalness: 0 });
}

/** Points drawn as soft puffs, each its own size (m) and opacity, cool or warm. */
function puffMaterial(tex: THREE.Texture, cool: string, warm: string): THREE.ShaderMaterial {
  return new THREE.ShaderMaterial({
    uniforms: { uTex: { value: tex }, uScale: { value: 400 }, uCool: { value: new THREE.Color(cool) }, uWarm: { value: new THREE.Color(warm) } },
    vertexShader: `
      attribute float aAlpha;
      attribute float aSize;
      attribute float aWarm;
      uniform float uScale;
      varying float vAlpha;
      varying float vWarm;
      void main() {
        vec4 mv = modelViewMatrix * vec4(position, 1.0);
        gl_Position = projectionMatrix * mv;
        gl_PointSize = aSize * uScale / max(0.05, -mv.z);
        vAlpha = aAlpha;
        vWarm = aWarm;
      }`,
    fragmentShader: `
      uniform sampler2D uTex;
      uniform vec3 uCool;
      uniform vec3 uWarm;
      varying float vAlpha;
      varying float vWarm;
      void main() {
        vec4 t = texture2D(uTex, gl_PointCoord);
        float a = t.a * vAlpha;
        if (a < 0.003) discard;
        vec3 c = mix(uCool, uWarm, vWarm) * (0.78 + 0.22 * t.r);
        gl_FragColor = vec4(c, a);
        #include <colorspace_fragment>
      }`,
    transparent: true,
    depthWrite: false,
    toneMapped: false,
  });
}

export function mountHoodLab(root: HTMLElement): void {
  const mountT0 = performance.now();
  const stageEl = root.querySelector<HTMLElement>("[data-hc-stage]");
  const canvas = root.querySelector<HTMLCanvasElement>("[data-hc-canvas]");
  if (!stageEl || !canvas) return;
  const maybeRenderer = createStudioRenderer(canvas, false);
  if (!maybeRenderer) return;
  const renderer = maybeRenderer;
  renderer.toneMappingExposure = 1.0;
  const phone = window.matchMedia("(max-width: 720px)").matches;
  const detail = hoodDetail();

  const scene = new THREE.Scene();
  scene.background = backdrop();
  scene.fog = new THREE.Fog(0x151510, 8, 18);
  scene.environment = studioEnvironment(renderer);
  scene.environmentIntensity = 0.55;
  const lights = studioLights(scene, new THREE.Vector3(0, 1, 0.8), 4, phone ? 1024 : 2048);
  lights.key.position.set(-3.2, 4.4, 4.6);
  lights.key.intensity = 1.35;
  lights.key.shadow.normalBias = 0.02;
  lights.rim.position.set(3.4, 3.2, 2.6);
  lights.rim.intensity = 0.5;
  lights.hemi.intensity = 0.5;
  const p = palette();
  {
    const floor = tilePlane(16, 6);
    floor.rotation.x = -Math.PI / 2;
    floor.position.set(0, 0, 2.6);
    scene.add(floor);
    // The wall: charcoal glazed tile from the floor to the hood's top, painted above.
    const top = EDGE_Y + HOOD.h;
    const tiled = new THREE.Mesh(new THREE.PlaneGeometry(16, top), tileWall(16, top));
    tiled.position.set(0, top / 2, -0.002);
    tiled.receiveShadow = true;
    scene.add(tiled);
    const paint = new THREE.Mesh(new THREE.PlaneGeometry(16, 1.2), new THREE.MeshStandardMaterial({ color: 0x34312c, roughness: 0.94 }));
    paint.position.set(0, top + 0.6, -0.002);
    scene.add(paint);
    // The room's corner, 2.7 m past the hood's end: the same tile and paint, running out into the kitchen.
    const SIDE = 7;
    const sideTile = new THREE.Mesh(new THREE.PlaneGeometry(SIDE, top), tileWall(SIDE, top));
    sideTile.rotation.y = -Math.PI / 2;
    sideTile.position.set(HOOD_HALF + 2.7, top / 2, SIDE / 2 - 0.002);
    sideTile.receiveShadow = true;
    scene.add(sideTile);
    const sidePaint = new THREE.Mesh(new THREE.PlaneGeometry(SIDE, 1.2), paint.material);
    sidePaint.rotation.y = -Math.PI / 2;
    sidePaint.position.set(HOOD_HALF + 2.7, top + 0.6, SIDE / 2 - 0.002);
    scene.add(sidePaint);
    const s = new Shop(detail, "hood_duct");
    // The exhaust duct from the hood's collar up through the ceiling.
    s.add(p.ssDull, rbox(s, 0.52, 1.2, 0.27, 0.004), [0, top + 0.6, 0.18]);
    const g = s.build(true).group;
    scene.add(g);
  }

  /* The hood, against the wall. */
  const hood = buildFixture(HOOD.id, HOOD, detail);
  hood.group.position.set(-HOOD.w / 2, EDGE_Y, HOOD.d);
  scene.add(hood.group);

  /* The broilers on their stands, and a pair of steaks on each. */
  const broiler = buildFixture(BROILER.id, BROILER, detail);
  broiler.hooks.set?.("heat", 0.85);
  const stand = buildFixture(STAND.id, STAND, detail);
  const steak = buildFixture("ARCH-FOOD-STEAK-RIBEYE", STEAK, STOVE_STEAK.detail);
  steak.hooks.set?.("browning", 1);
  steak.hooks.set?.("core", 60);
  /** Each broiler's stand and steaks ride with it as it slides back. */
  const units = XS.map((x, i) => {
    const g = new THREE.Group();
    const b = i === 0 ? broiler.group : broiler.group.clone();
    b.position.set(x - BROILER.w / 2, STAND.h, 0);
    g.add(b);
    const st = i === 0 ? stand.group : stand.group.clone();
    st.position.set(x - STAND.w / 2, 0, -0.012);
    st.scale.y = STAND.h / TABLE_Y;
    g.add(st);
    for (const [dx, dz, yaw] of [
      [-0.2, -0.33, 0.25],
      [0.17, -0.42, -0.4],
    ] as const) {
      const s = i === 0 && dx < 0 ? steak.group : steak.group.clone();
      const holder = new THREE.Group();
      s.position.set(-STEAK.w / 2, 0, STEAK.d / 2);
      holder.add(s);
      // On the sloped grates, which fall toward the trough at the front, settled a millimetre into the bars.
      const L = GRATE_Z0 - GRATE_Z1;
      holder.position.set(x + dx + i * 0.01, STAND.h + GRATE_YF + ((GRATE_Z0 - dz) / L) * GRATE_RISE - 0.001, dz);
      holder.rotation.set(Math.atan2(GRATE_RISE, L), yaw, 0);
      g.add(holder);
    }
    scene.add(g);
    return g;
  });

  /* The seal: one stainless sheet from the broilers' backs to the wall, along the whole line, at their tops. */
  const sealMat = p.ss;
  const seal = new THREE.Mesh(new THREE.BoxGeometry(2 * LINE_HALF + 0.02, 0.0015, 1), sealMat);
  seal.castShadow = true;
  seal.receiveShadow = true;
  const SEAL_Y = STAND.h + BROILER.h - 0.004;
  scene.add(seal);
  // Its folded front lip down the broilers' backs.
  const lip = new THREE.Mesh(new THREE.BoxGeometry(2 * LINE_HALF + 0.02, 0.05, 0.0015), sealMat);
  lip.castShadow = true;
  scene.add(lip);

  /* The cook at the middle broiler, tongs in hand, turning the steaks. */
  const cook = makeStaff("p50", 2, { seed: 47, apron: 0x2b3a55, sleeves: "cuffed" });
  scene.add(cook.group);
  const COOK_X = 0.05;
  // Locking tongs in the right hand: their hinge end in the fist, the cups 0.185 m on down the grip (as the hero
  // kitchen's grill cook holds them).
  const tongs = new THREE.Group();
  {
    const m = buildFixture("NSF-TONGS-12IN", { w: 0.055, d: 0.305, h: 0.038 }, 1);
    m.group.rotation.y = Math.PI;
    m.group.position.set(0.0275, -0.019, -0.09);
    tongs.add(m.group);
  }
  const grateAt = new THREE.Vector3();

  /* The smoke as soft points and the air it draws in as faint streaks, updated each frame. */
  const NS = phone ? 950 : 1500;
  const NA = phone ? 150 : 240;
  interface Cloud {
    pos: Float32Array;
    alpha: Float32Array;
    size: Float32Array;
    warm: Float32Array;
    geo: THREE.BufferGeometry;
    pts: THREE.Points;
  }
  function cloud(n: number, mat: THREE.ShaderMaterial): Cloud {
    const pos = new Float32Array(n * 3);
    const alpha = new Float32Array(n);
    const size = new Float32Array(n);
    const warm = new Float32Array(n);
    const geo = new THREE.BufferGeometry();
    geo.setAttribute("position", new THREE.BufferAttribute(pos, 3).setUsage(THREE.DynamicDrawUsage));
    geo.setAttribute("aAlpha", new THREE.BufferAttribute(alpha, 1).setUsage(THREE.DynamicDrawUsage));
    geo.setAttribute("aSize", new THREE.BufferAttribute(size, 1).setUsage(THREE.DynamicDrawUsage));
    geo.setAttribute("aWarm", new THREE.BufferAttribute(warm, 1).setUsage(THREE.DynamicDrawUsage));
    geo.boundingSphere = new THREE.Sphere(new THREE.Vector3(0, 1.8, 1), 6);
    const pts = new THREE.Points(geo, mat);
    pts.frustumCulled = false;
    scene.add(pts);
    return { pos, alpha, size, warm, geo, pts };
  }
  const smokeMat = puffMaterial(puffTexture(64, false), "#b4b9c1", "#dccdb9");
  const smoke = cloud(NS, smokeMat);
  smoke.pts.renderOrder = 2;
  /** Each streak: its head where the air is, its tail where it was a moment ago, fading to nothing. */
  const airPos = new Float32Array(NA * 6);
  const airAlpha = new Float32Array(NA * 2);
  const airGeo = new THREE.BufferGeometry();
  airGeo.setAttribute("position", new THREE.BufferAttribute(airPos, 3).setUsage(THREE.DynamicDrawUsage));
  airGeo.setAttribute("aAlpha", new THREE.BufferAttribute(airAlpha, 1).setUsage(THREE.DynamicDrawUsage));
  const airMat = new THREE.ShaderMaterial({
    uniforms: { uColor: { value: new THREE.Color(AIR) } },
    vertexShader: `
      attribute float aAlpha;
      varying float vA;
      void main() {
        vA = aAlpha;
        gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0);
      }`,
    fragmentShader: `
      uniform vec3 uColor;
      varying float vA;
      void main() {
        if (vA < 0.003) discard;
        gl_FragColor = vec4(uColor, vA);
        #include <colorspace_fragment>
      }`,
    transparent: true,
    depthWrite: false,
    toneMapped: false,
  });
  const airLines = new THREE.LineSegments(airGeo, airMat);
  airLines.frustumCulled = false;
  airLines.renderOrder = 1;
  scene.add(airLines);
  // Each smoke puff: its age and life (s), its broiler's middle, whether it escapes, a seed, the time it got out.
  const sAge = new Float32Array(NS).fill(-1);
  const sLife = new Float32Array(NS);
  const sX0 = new Float32Array(NS);
  const sOut = new Uint8Array(NS);
  const sSeed = new Float32Array(NS);
  const sEsc = new Float32Array(NS);
  const sDie = new Float32Array(NS);
  // Each streak of air: where it is and how it moves, its age and life, and whether it comes up the gap (1) or from
  // the front (0).
  const aP = new Float32Array(NA * 3);
  const aV = new Float32Array(NA * 3);
  const aAge = new Float32Array(NA).fill(-1);
  const aLife = new Float32Array(NA);
  const aGap = new Uint8Array(NA);
  const aSeed = new Float32Array(NA);

  /* The camera: from in front and to the left, so the gap behind the broilers shows past the hood's end. */
  const camera = new THREE.PerspectiveCamera(32, 16 / 9, 0.05, 60);
  const TARGET = new THREE.Vector3(-0.4, 1.4, 0.75);
  const VIEW = new THREE.Vector3();
  const keep: THREE.Vector3[] = [];
  const aim = () => {
    VIEW.set(-1, phone ? 0.16 : 0.14, phone ? 0.62 : 0.5).normalize();
    keep.length = 0;
    for (const sx of [-1, 1]) {
      keep.push(new THREE.Vector3(sx * HOOD_HALF, EDGE_Y + HOOD.h, 0), new THREE.Vector3(sx * HOOD_HALF, EDGE_Y, HOOD.d));
      keep.push(new THREE.Vector3(sx * LINE_HALF, 0, HOOD.d));
    }
    // The smoke that escapes, over the hood's front along its length, and the cook's head.
    keep.push(new THREE.Vector3(-0.6, EDGE_Y + HOOD.h + 0.25, HOOD.d + 0.4));
    keep.push(new THREE.Vector3(HOOD_HALF - 0.3, EDGE_Y + HOOD.h + 0.3, HOOD.d + 0.5));
    keep.push(new THREE.Vector3(0, 1.8, HOOD.d + 0.5));
  };
  const controls = createScrollSafeOrbit(camera, canvas, [TARGET.x, TARGET.y, TARGET.z], { minDistance: 1, maxDistance: 14, maxPolarAngle: Math.PI * 0.49 });
  let userMoved = false;
  controls.addEventListener("start", () => {
    userMoved = true;
  });
  const hudEl = stageEl.querySelector<HTMLElement>(".hc-hud");
  const hudRoom = (): number => {
    if (!hudEl) return 0;
    const parts = ["[data-hc-hud-set]", "[data-hc-hud-doing]", "[data-hc-hud-note]"].map((k) => hudEl.querySelector<HTMLElement>(k));
    const was = parts.map((e) => e?.textContent ?? "");
    let most = 0;
    for (const [a, b, c] of longestHud()) {
      if (parts[0]) parts[0].textContent = a;
      if (parts[1]) parts[1].textContent = b;
      if (parts[2]) parts[2].textContent = c;
      most = Math.max(most, hudEl.offsetTop + hudEl.offsetHeight);
    }
    parts.forEach((e, i) => {
      if (e) e.textContent = was[i]!;
    });
    return most;
  };
  const legendEl = stageEl.querySelector<HTMLElement>(".hc-legend");
  const ndc = new THREE.Vector3();
  let fitW = 0;
  let fitH = 0;
  const frame = () => {
    if (userMoved) return;
    aim();
    const w = Math.max(1, stageEl.clientWidth);
    const h = Math.max(1, stageEl.clientHeight);
    fitW = w;
    fitH = h;
    camera.aspect = Math.max(0.2, w / h);
    camera.clearViewOffset();
    camera.updateProjectionMatrix();
    const top = Math.min(0.88, 1 - (2 * (hudRoom() + 10)) / h);
    const bot = -1 + (2 * ((legendEl?.offsetHeight ?? 0) + 18)) / h;
    let dist = 2;
    let lo = 0;
    let hi = 0;
    let left = 0;
    let right = 0;
    for (; dist < 30; dist += 0.05) {
      camera.position.copy(TARGET).addScaledVector(VIEW, dist);
      camera.lookAt(TARGET);
      camera.updateMatrixWorld(true);
      lo = left = Infinity;
      hi = right = -Infinity;
      for (const k of keep) {
        ndc.copy(k).project(camera);
        lo = Math.min(lo, ndc.y);
        hi = Math.max(hi, ndc.y);
        left = Math.min(left, ndc.x);
        right = Math.max(right, ndc.x);
      }
      if (right - left < 1.9 && hi - lo < top - bot) break;
    }
    const shiftY = (top + bot) / 2 - (hi + lo) / 2;
    const shiftX = (left + right) / 2;
    if (Math.abs(shiftY) > 0.005 || Math.abs(shiftX) > 0.005) {
      camera.setViewOffset(w, h, (shiftX * w) / 2, (shiftY * h) / 2, w, h);
      camera.updateProjectionMatrix();
    }
    const fog = scene.fog as THREE.Fog;
    fog.near = dist * 1.5;
    fog.far = dist * 3.4;
    controls.target.copy(TARGET);
    controls.update();
  };
  const gate = shaderGate(renderer, camera);
  bindCanvasResize(stageEl, renderer, camera, gate);
  const sizeScale = () => {
    const hPx = renderer.domElement.height;
    const k = hPx / (2 * Math.tan(THREE.MathUtils.degToRad(camera.fov) / 2));
    smokeMat.uniforms.uScale!.value = k;
  };

  /* ---------------------------------------------------------------------------------------------------------- */
  /* The page.                                                                                                   */
  /* ---------------------------------------------------------------------------------------------------------- */
  const q = <T extends Element = HTMLElement>(s: string) => root.querySelector<T>(s);
  const atBtns = Array.from(root.querySelectorAll<HTMLButtonElement>("[data-hc-at]"));
  const sealBtns = Array.from(root.querySelectorAll<HTMLButtonElement>("[data-hc-seal]"));
  const fanIn = q<HTMLInputElement>("[data-hc-fan]");
  const fanVal = q("[data-hc-fan-val]");
  const statusEl = q("[data-hc-status]");
  const chartEl = q("[data-hc-chart]");
  const hudSet = q("[data-hc-hud-set]");
  const hudDoing = q("[data-hc-hud-doing]");
  const hudNote = q("[data-hc-hud-note]");
  const tiles = {
    need: { v: q("[data-hc-need]"), s: q("[data-hc-need-s]") },
    fan: { box: q("[data-hc-fan-box]"), v: q("[data-hc-fan-t]"), s: q("[data-hc-fan-s]") },
    gap: { v: q("[data-hc-gap]"), s: q("[data-hc-gap-s]") },
    warm: { v: q("[data-hc-warm]"), s: q("[data-hc-warm-s]") },
  };
  const fanKey = q("[data-hc-key-fan]");
  const tagGap = q<HTMLElement>('[data-hc-tag="gap"]');
  const tagLip = q<HTMLElement>('[data-hc-tag="lip"]');

  let setup: Setup = { ...TOUR[0]! };
  let fan = CODE_CFM;
  /** The tour runs until the reader picks where the broilers stand or the seal, moves the fan or takes hold of the
   *  view; with reduced motion asked for, it does not run. */
  let touring = !(window.matchMedia?.("(prefers-reduced-motion: reduce)").matches ?? false);
  let tourT = 0;
  /** What is drawn, eased toward the setup: the broilers' fronts (m), the seal let down (0..1). */
  let frontZ = frontOf(setup.at);
  let sealK = setup.sealed ? 1 : 0;

  const setName = (s: Setup) => `${SHORT_AT[s.at]}, gap ${s.sealed ? "sealed" : "open"}`;
  const hudLines = (s: Setup, f: number): [string, string, string] => [setName(s), `Needs ${cfm(needCfm(s))}; the fan pulls ${cfm(f)}`, catches(f, s) ? "All the smoke caught" : "Smoke escapes"];
  /** The longest the HUD's lines get (the figures are set in even widths). */
  function longestHud(): [string, string, string][] {
    return [hudLines({ at: 0, sealed: true }, 6000), hudLines({ at: 0, sealed: false }, 2000)];
  }

  function renderControls(): void {
    for (const b of atBtns) b.setAttribute("aria-pressed", String(Number(b.dataset.hcAt) === setup.at));
    for (const b of sealBtns) b.setAttribute("aria-pressed", String((b.dataset.hcSeal === "1") === setup.sealed));
    if (fanIn) {
      fanIn.value = String(fan);
      fanIn.style.setProperty("--p", `${(((fan - 2000) / 4000) * 100).toFixed(2)}%`);
    }
    if (fanVal) fanVal.textContent = cfm(fan);
  }
  function renderTiles(): void {
    const need = needCfm(setup);
    const ok = catches(fan, setup);
    if (tiles.need.v) tiles.need.v.textContent = cfm(need);
    if (tiles.need.s) tiles.need.s.textContent = `${perFoot(need)} cfm a foot of hood\nthe lab's capture rate`;
    tiles.fan.box?.setAttribute("data-ok", String(ok));
    if (tiles.fan.v) tiles.fan.v.textContent = ok ? "Caught" : "Escapes";
    if (tiles.fan.s) tiles.fan.s.textContent = ok ? (fan === need ? `at ${cfm(fan)}\nnothing to spare` : `at ${cfm(fan)}\n${cfm(fan - need)} to spare`) : `at ${cfm(fan)}\n${cfm(need - fan)} short`;
    const set = SETS[setup.at];
    if (tiles.gap.v) tiles.gap.v.textContent = setup.sealed ? "Sealed" : `${set.gapIn} in`;
    if (tiles.gap.s) tiles.gap.s.textContent = `${setup.sealed ? `a ${set.gapIn} in gap, closed` : "open to the floor"}\nbroilers ${SHORT_AT[setup.at].toLowerCase()}`;
    if (tiles.warm.v) tiles.warm.v.textContent = btuh(warmBtuh(fan, COLD_OUT_F));
    if (tiles.warm.s) tiles.warm.s.textContent = `the fan's ${cfm(fan)}\nfrom ${COLD_OUT_F} °F to 70 °F`;
  }
  function renderStatus(): void {
    if (!statusEl) return;
    const n = caughtCount(fan);
    const parts: string[] = [];
    const need = needCfm(setup);
    const ok = catches(fan, setup);
    parts.push(
      ok
        ? `With the broilers ${SAID_AT[setup.at]} and the gap behind them ${setup.sealed ? "sealed" : "open"}, a fan pulling ${cfm(fan)} catches all their smoke: the lab measured ${cfm(need)}.`
        : `With the broilers ${SAID_AT[setup.at]} and the gap behind them ${setup.sealed ? "sealed" : "open"}, smoke escapes the hood: the fan pulls ${cfm(fan)}, and the lab measured ${cfm(need)} to catch it all.`,
    );
    // Where the broilers must stand for this fan, open and sealed.
    const firstAt = (sealed: boolean) => SETUPS.find((s) => s.sealed === sealed && catches(fan, s))?.at ?? null;
    const open = firstAt(false);
    const sealed = firstAt(true);
    if (n === 0) parts.push(`At ${cfm(fan)} no setup is caught; the least that catches one is ${cfm(needCfm({ at: 2, sealed: true }))}, 12 inches back and sealed.`);
    else if (n === 6) parts.push(`At ${cfm(fan)} every setup is caught, even flush with the edge and open.`);
    else {
      const where = (at: SetIndex | null, how: string) => (at === null ? `${how}, nowhere under this hood` : at === 0 ? `${how}, anywhere, even flush with the edge` : at === 1 ? `${how}, 6 or 12 inches back` : `${how}, 12 inches back`);
      parts.push(`At ${cfm(fan)} the broilers can stand: ${where(open, "with the gap open")}; ${where(sealed, "with it sealed")}.`);
    }
    // What the worst setup costs over the code's least, and how low the best one goes.
    const worst = needCfm({ at: 0, sealed: false });
    const best = needCfm({ at: 2, sealed: true });
    parts.push(`Flush and open, catching it all takes ${cfm(worst)}, ${(worst - CODE_CFM).toLocaleString("en-US")} more than the code's ${CODE_CFM.toLocaleString("en-US")} and ${btuh(warmBtuh(worst - CODE_CFM, COLD_OUT_F))} more to warm on a ${COLD_OUT_F} °F morning; 12 inches back and sealed, ${cfm(best)} catches it all, for a hood listed to run that low.`);
    statusEl.textContent = parts.join(" ");
    statusEl.dataset.state = ok ? "ok" : "warn";
  }

  /* The chart: the six setups' measured rates, against the fan and the code's least. */
  function renderChart(): void {
    if (!chartEl) return;
    const W = Math.max(280, Math.round(chartEl.clientWidth));
    const narrow = W < 520;
    const H = narrow ? 210 : 230;
    const L = narrow ? 40 : 50;
    const R = W - (narrow ? 8 : 14);
    const T = 18;
    const B = H - 34;
    const top = 6000;
    const y = (v: number) => B - (Math.min(top, Math.max(0, v)) / top) * (B - T);
    let svg = `<svg viewBox="0 0 ${W} ${H}" width="${W}" height="${H}" class="rs-svg bo-svg hc-svg"><defs><pattern id="hc-hatch" width="6" height="6" patternUnits="userSpaceOnUse" patternTransform="rotate(45)"><line x1="0" y1="0" x2="0" y2="6" class="hatch" /></pattern></defs>`;
    for (const v of [0, 2000, 4000, 6000]) {
      svg += `<line class="grid" x1="${L}" x2="${R}" y1="${y(v).toFixed(1)}" y2="${y(v).toFixed(1)}" />`;
      svg += `<text class="ya" x="${L - 6}" y="${(y(v) + 4).toFixed(1)}" text-anchor="end">${v ? `${v / 1000}k` : "0"}</text>`;
    }
    svg += `<text class="ya" x="${L - 6}" y="${T - 8}" text-anchor="end">cfm</text>`;
    const gw = (R - L) / 3;
    const bw = Math.min(46, gw * 0.3);
    for (const s of SETUPS) {
      const need = needCfm(s);
      const ok = catches(fan, s);
      const x0 = L + gw * s.at + gw / 2 + (s.sealed ? 3 : -bw - 3);
      const cur = s.at === setup.at && s.sealed === setup.sealed;
      svg += `<g class="bar${cur ? " sel" : ""}" data-hc-bar="${s.at}${s.sealed ? "s" : "o"}" style="--c:${ok ? OK : SPILL}">`;
      svg += `<rect class="hit" x="${(x0 - 2).toFixed(1)}" y="${T}" width="${(bw + 4).toFixed(1)}" height="${B - T}" />`;
      const box = `x="${x0.toFixed(1)}" y="${y(need).toFixed(1)}" width="${bw.toFixed(1)}" height="${(B - y(need)).toFixed(1)}" rx="2"`;
      svg += `<rect class="b${s.sealed ? " sealed" : ""}" ${box} />`;
      if (s.sealed) svg += `<rect class="stripes" ${box} fill="url(#hc-hatch)" />`;
      svg += `<text class="v" x="${(x0 + bw / 2).toFixed(1)}" y="${(y(need) - 5).toFixed(1)}" text-anchor="middle">${narrow ? `${(need / 1000).toFixed(1)}k` : need.toLocaleString("en-US")}</text>`;
      svg += `</g>`;
    }
    for (let at = 0; at < 3; at++) svg += `<text class="xa" x="${(L + gw * at + gw / 2).toFixed(1)}" y="${H - 18}" text-anchor="middle">${narrow && at === 0 ? "Flush" : SHORT_AT[at]}</text>`;
    svg += `<text class="xa sub" x="${(L + (R - L) / 2).toFixed(1)}" y="${H - 4}" text-anchor="middle">Solid: gap open · striped: sealed</text>`;
    // The code's least, and the fan (named in the keys above, where no bar's figure can sit on them).
    svg += `<line class="fan" x1="${L}" x2="${R}" y1="${y(fan).toFixed(1)}" y2="${y(fan).toFixed(1)}" style="--c:${FAN}" />`;
    svg += `<line class="code" x1="${L}" x2="${R}" y1="${y(CODE_CFM).toFixed(1)}" y2="${y(CODE_CFM).toFixed(1)}" />`;
    svg += `</svg>`;
    if (fanKey) fanKey.textContent = `Fan, ${cfm(fan)}`;
    chartEl.innerHTML = svg;
    for (const g of chartEl.querySelectorAll<SVGGElement>("[data-hc-bar]")) {
      g.addEventListener("click", () => {
        const k = g.dataset.hcBar ?? "0o";
        pick({ at: Number(k[0]) as SetIndex, sealed: k[1] === "s" });
      });
    }
  }
  function renderAll(): void {
    renderControls();
    renderTiles();
    renderStatus();
    renderChart();
    lastHud = "";
  }
  /** When the setup or the fan changes, the smoke still low in the plume goes where the new one sends it: only as
   *  much of it changes course as the share escaping changes, and none already near the hood's edge or out. */
  let fateShare = spillShare(fan, needCfm(setup));
  function refate(): void {
    const sp = spillShare(fan, needCfm(setup));
    const was = fateShare;
    fateShare = sp;
    if (sp === was) return;
    const toCaught = sp < was ? 1 - sp / was : 0;
    const toOut = sp > was ? (sp - was) / (1 - was) : 0;
    for (let i = 0; i < NS; i++) {
      if (sAge[i]! < 0 || sEsc[i]! >= 0 || smoke.pos[i * 3 + 1]! > EDGE_Y - 0.25) continue;
      if (sOut[i]) {
        if (rnd() < toCaught) sOut[i] = 0;
      } else if (rnd() < toOut) sOut[i] = 1;
    }
  }
  function pick(s: Setup): void {
    touring = false;
    show(s);
  }
  function show(s: Setup): void {
    if (s.at === setup.at && s.sealed === setup.sealed) return;
    setup = { ...s };
    refate();
    renderAll();
  }

  /* The readouts on the stage. */
  let lastHud = "";
  function showHud(): void {
    const sig = `${setup.at}${setup.sealed ? "s" : "o"}${fan}`;
    if (sig === lastHud) return;
    lastHud = sig;
    const ok = catches(fan, setup);
    const [a, b, c] = hudLines(setup, fan);
    if (hudSet) hudSet.textContent = a;
    if (hudDoing) hudDoing.textContent = b;
    if (hudNote) {
      hudNote.textContent = c;
      hudNote.dataset.kind = ok ? "caught" : "spill";
    }
  }
  const tagAt = new THREE.Vector3();
  const tagSize = new Map<HTMLElement, { w: number; h: number }>();
  const STEM = 26;
  function placeTag(el: HTMLElement | null, at: THREE.Vector3, text: [string, string], warn: boolean, colour: string, W: number, H: number): void {
    if (!el) return;
    tagAt.copy(at).project(camera);
    const on = tagAt.z < 1 && Math.abs(tagAt.x) < 0.96 && Math.abs(tagAt.y) < 0.96;
    el.hidden = !on;
    if (!on) return;
    const sig = text.join("|") + warn;
    if (el.dataset.last !== sig || !tagSize.has(el)) {
      el.dataset.last = sig;
      const i = el.querySelector("i");
      const b = el.querySelector("b");
      if (i) i.textContent = text[0];
      if (b) b.textContent = text[1];
      el.toggleAttribute("data-boiled", warn);
      el.style.setProperty("--c", colour);
      el.style.setProperty("--stem", `${STEM}px`);
      tagSize.set(el, { w: el.offsetWidth, h: el.offsetHeight });
    }
    const sz = tagSize.get(el)!;
    const px = ((tagAt.x + 1) / 2) * W;
    const py = ((1 - tagAt.y) / 2) * H;
    const x = Math.min(W - sz.w / 2 - 4, Math.max(sz.w / 2 + 4, px));
    el.style.transform = `translate(${x.toFixed(1)}px, ${(py - STEM).toFixed(1)}px)`;
    const sx = `${Math.round(px - x)}px`;
    if (el.dataset.sx !== sx) {
      el.dataset.sx = sx;
      el.style.setProperty("--sx", sx);
    }
  }
  const gapPoint = new THREE.Vector3();
  const lipPoint = new THREE.Vector3();
  function showTags(): void {
    const W = stageEl!.clientWidth;
    const H = stageEl!.clientHeight;
    const back = frontZ - BROILER.d;
    gapPoint.set(-LINE_HALF + 0.12, SEAL_Y + 0.01, Math.max(0.02, back / 2));
    placeTag(tagGap, gapPoint, setup.sealed ? ["Gap behind", "sealed"] : ["Gap behind", `${SETS[setup.at].gapIn} in, open`], false, AIR, W, H);
    lipPoint.set(-0.55, EDGE_Y - 0.02, HOOD.d);
    const ok = catches(fan, setup);
    placeTag(tagLip, lipPoint, ok ? ["Hood's edge", "all caught"] : ["Hood's edge", "smoke escapes"], !ok, ok ? OK : SPILL, W, H);
  }

  /* The stage: the broilers and the seal eased to the setup, the cook, the smoke and the air. */
  const rnd = (() => {
    let s = 0x9e3779b9;
    return () => {
      s ^= s << 13;
      s ^= s >>> 17;
      s ^= s << 5;
      return ((s >>> 0) % 1_000_000) / 1_000_000;
    };
  })();
  let emitAcc = 0;
  let nextUnit = 0;
  let airAcc = 0;
  /** How far forward the smoke leans as it rises (horizontal over vertical): the air up an open gap pushes it
   *  forward, the more the wider the gap; sealed, the air from the front pushes it back. */
  const leanNow = () => (1 - sealK) * (0.1 + 0.55 * Math.max(0, frontZ - BROILER.d)) - sealK * 0.07;
  function emitSmoke(i: number, x0: number, spill: number): void {
    sAge[i] = 0;
    sLife[i] = 4.6 + rnd() * 1.6;
    sX0[i] = x0;
    sOut[i] = rnd() < spill ? 1 : 0;
    sSeed[i] = rnd();
    sEsc[i] = -1;
    sDie[i] = -1;
    const k = i * 3;
    // Most smoke comes off the steaks; the rest off the bars where the fat burns.
    const onSteak = rnd() < 0.55;
    smoke.pos[k] = x0 + (onSteak ? (rnd() < 0.5 ? -0.2 : 0.17) + (rnd() - 0.5) * 0.16 : (rnd() - 0.5) * 0.8);
    smoke.pos[k + 1] = COOK_Y + 0.02 + rnd() * 0.04;
    smoke.pos[k + 2] = frontZ + (onSteak ? -0.33 - rnd() * 0.1 : GRATE_Z0 + (GRATE_Z1 - GRATE_Z0) * rnd());
  }
  function emitAir(i: number): void {
    aAge[i] = 0;
    aSeed[i] = rnd();
    const back = frontZ - BROILER.d;
    const gap = (1 - sealK) * back > 0.05 && rnd() < 0.72 * (1 - sealK);
    aGap[i] = gap ? 1 : 0;
    const k = i * 3;
    aP[k] = (rnd() * 2 - 1) * (LINE_HALF - 0.08);
    if (gap) {
      aLife[i] = 2.8;
      aP[k + 1] = 0.1 + rnd() * 0.45;
      aP[k + 2] = 0.03 + rnd() * Math.max(0.01, back - 0.06);
    } else {
      aLife[i] = 3.2;
      aP[k + 1] = STAND.h + BROILER.h + 0.03 + rnd() * 0.3;
      aP[k + 2] = frontZ + 0.5 + rnd() * 0.4;
      // Not through the cook, who stands in front of the middle broiler.
      if (Math.abs(aP[k]! - COOK_X) < 0.34) aP[k] = COOK_X + Math.sign(aP[k]! - COOK_X || 1) * (0.34 + rnd() * 0.1);
    }
    aV[k] = aV[k + 1] = aV[k + 2] = 0;
  }
  /** The kitchen's stirring: one slow field that every puff moves in, so that neighbours move together into wisps. */
  const swirlX = (x: number, y: number, z: number, t: number, sd: number) => 0.11 * Math.sin(1.3 * t + 3.1 * y + 2.3 * x) + 0.05 * Math.sin(2.9 * t - 5.3 * y + 1.7 * z + sd);
  const swirlZ = (x: number, y: number, t: number, sd: number) => 0.07 * Math.cos(1.1 * t + 2.7 * y - 1.9 * x) + 0.04 * Math.sin(3.3 * t + 4.1 * y + 2 * sd);
  function stepSmoke(dt: number, t: number): void {
    const need = needCfm(setup);
    const spill = spillShare(fan, need);
    const lean = leanNow();
    const draw = 0.35 + 0.55 * (fan / CODE_CFM);
    // Emission: each broiler's share, steady.
    emitAcc += dt * (phone ? 190 : 300);
    for (let i = 0; i < NS && emitAcc >= 1; i++) {
      if (sAge[i]! >= 0) continue;
      emitSmoke(i, XS[nextUnit]!, spill);
      nextUnit = (nextUnit + 1) % 3;
      emitAcc -= 1;
    }
    emitAcc = Math.min(emitAcc, 4);
    for (let i = 0; i < NS; i++) {
      let age = sAge[i]!;
      const k = i * 3;
      if (age < 0) {
        smoke.alpha[i] = 0;
        continue;
      }
      age += dt;
      sAge[i] = age;
      let x = smoke.pos[k]!;
      let y = smoke.pos[k + 1]!;
      let z = smoke.pos[k + 2]!;
      const sd = sSeed[i]!;
      const rise = 0.62 * (0.85 + 0.3 * sd);
      const tx = swirlX(x, y, z, t, sd);
      const tz = swirlZ(x, y, t, sd);
      // Positions are kept in single precision: a puff put at the edge reads back a hair under it, so the test allows it.
      const inHood = y >= EDGE_Y - 0.002 && z < CAPTURE_Z && Math.abs(x) < HOOD_HALF - 0.02;
      let vx: number;
      let vy: number;
      let vz: number;
      if (sEsc[i]! >= 0 || (y >= EDGE_Y - 0.005 && !inHood)) {
        // Out in the kitchen: up the hood's face and away into the room, thinning.
        if (sEsc[i]! < 0) sEsc[i] = age;
        vx = tx * 1.4 + (x - sX0[i]!) * 0.08;
        vy = 0.28 + 0.1 * sd;
        vz = 0.15 + 0.06 * sd + tz;
        if (z < HOOD.d + 0.03) z = HOOD.d + 0.03;
      } else if (inHood) {
        // In the hood's reservoir: up to its top, then back along it to the filters, faster the harder the fan pulls.
        vx = tx * 0.5;
        vy = y < HOOD_TOP_IN ? rise * 0.6 : 0.02;
        vz = -draw * smooth((y - EDGE_Y) / 0.25);
        if (z < FILTER_Z && sDie[i]! < 0) sDie[i] = age;
      } else {
        // Rising off the grates: widening, leaning with the air that feeds it.
        const h = Math.max(0, (y - COOK_Y) / (EDGE_Y - COOK_Y));
        vx = (x - sX0[i]!) * 0.18 + tx;
        vy = rise;
        vz = lean * rise + tz;
        if (sOut[i]) {
          // Escaping: carried forward, under the hood's front plenum, held down until clear of its edge.
          const out = HOOD.d + 0.09;
          if (h > 0.45 && z < out) vz += (out - z) * 4.2 * h;
          if (y > EDGE_Y - 0.12 && z < HOOD.d + 0.04) vy *= 0.12;
        } else {
          // Caught: drawn back as it nears the hood, clear of its plenum, and in from its open ends.
          const lim = CAPTURE_Z - 0.08;
          if (h > 0.35 && z > lim) vz -= (z - lim) * 7 * h;
          if (y > EDGE_Y - 0.1 && z > CAPTURE_Z - 0.04) vy *= 0.15;
          const xl = HOOD_HALF - 0.14;
          if (h > 0.5 && Math.abs(x) > xl) vx -= Math.sign(x) * (Math.abs(x) - xl) * 6 * h;
        }
      }
      x += vx * dt;
      let ny = y + vy * dt;
      z += vz * dt;
      // At the hood's edge, nothing crosses where it should not: what is caught goes up behind the plenum and inside
      // the ends; what escapes stays under the lip until it is clear of it.
      if (!inHood && sEsc[i]! < 0 && ny >= EDGE_Y - 0.005) {
        if (sOut[i]) {
          if (z < HOOD.d + 0.03) ny = Math.min(ny, EDGE_Y - 0.006);
        } else {
          z = Math.min(z, CAPTURE_Z - 0.04);
          x = Math.max(-HOOD_HALF + 0.05, Math.min(HOOD_HALF - 0.05, x));
          ny = Math.max(ny, EDGE_Y + 0.003);
        }
      }
      y = ny;
      // Nothing goes through the wall, or out through the hood's ends while it is inside.
      if (z < 0.03) z = 0.03;
      if (inHood) x = Math.max(-HOOD_HALF + 0.03, Math.min(HOOD_HALF - 0.03, x));
      smoke.pos[k] = x;
      smoke.pos[k + 1] = y;
      smoke.pos[k + 2] = z;
      // How it looks: grown as it rises, thinned as it spreads, gone into the filters or out in the room.
      const up = Math.max(0, y - COOK_Y);
      const esc = sEsc[i]! >= 0 ? age - sEsc[i]! : 0;
      const die = sDie[i]! >= 0 ? age - sDie[i]! : 0;
      let a = 0.2 * smooth(age / 0.3) * (1 - 0.3 * smooth(up / 1.1));
      if (esc > 0) a *= 1 - smooth(esc / 2.8);
      if (die > 0) a *= 1 - smooth(die / 0.5);
      if (age > sLife[i]! - 0.6) a *= Math.max(0, (sLife[i]! - age) / 0.6);
      smoke.alpha[i] = a;
      smoke.size[i] = 0.05 + 0.2 * smooth(up / 1.2) + 0.32 * smooth(esc / 2.4);
      smoke.warm[i] = 1 - smooth(up / 0.5);
      if (age >= sLife[i]! || (die > 0.5 && sDie[i]! >= 0) || (esc > 2.8 && sEsc[i]! >= 0) || y > 3.6) sAge[i] = -1;
    }
    for (const a of ["position", "aAlpha", "aSize", "aWarm"] as const) (smoke.geo.getAttribute(a) as THREE.BufferAttribute).needsUpdate = true;
  }
  function stepAir(dt: number, t: number): void {
    airAcc += dt * (phone ? 50 : 80);
    for (let i = 0; i < NA && airAcc >= 1; i++) {
      if (aAge[i]! >= 0) continue;
      emitAir(i);
      airAcc -= 1;
    }
    airAcc = Math.min(airAcc, 3);
    const back = frontZ - BROILER.d;
    const plumeZ = frontZ - 0.38;
    for (let i = 0; i < NA; i++) {
      let age = aAge[i]!;
      const k = i * 3;
      if (age < 0) {
        airAlpha[2 * i] = airAlpha[2 * i + 1] = 0;
        continue;
      }
      age += dt;
      aAge[i] = age;
      const sd = aSeed[i]!;
      let vy: number;
      let vz: number;
      if (aGap[i]) {
        // Up the gap behind the broilers, then over their backs into the smoke.
        if (aP[k + 1]! < SEAL_Y - 0.03 && aP[k + 2]! < back) {
          vy = 0.6 + 0.15 * sd;
          vz = 0;
        } else {
          vy = 0.28;
          vz = 0.55;
        }
        // Sealed meanwhile: it stops under the sheet and fades.
        if (sealK > 0.5 && aP[k + 1]! > SEAL_Y - 0.06 && aP[k + 2]! < back) {
          vy = 0;
          aAge[i] = Math.max(age, aLife[i]! - 0.4);
        }
      } else {
        // From the kitchen in front, low, toward the foot of the smoke, rising as it nears it.
        vy = 0.06 + 0.25 * smooth((frontZ + 0.1 - aP[k + 2]!) / 0.5);
        vz = -0.42 - 0.1 * sd;
      }
      const vx = 0.04 * Math.sin(t * 1.3 + sd * 20);
      // Eased, so a streak turns rather than kinks.
      const e = Math.min(1, dt * 6);
      aV[k] = aV[k]! + (vx - aV[k]!) * e;
      aV[k + 1] = aV[k + 1]! + (vy - aV[k + 1]!) * e;
      aV[k + 2] = aV[k + 2]! + (vz - aV[k + 2]!) * e;
      for (let c = 0; c < 3; c++) aP[k + c] = aP[k + c]! + aV[k + c]! * dt;
      const done = aGap[i] ? aP[k + 2]! > back + 0.34 : aP[k + 2]! < plumeZ;
      let a = 0.8 * smooth(age / 0.3);
      if (age > aLife[i]! - 0.4) a *= Math.max(0, (aLife[i]! - age) / 0.4);
      // The head where the air is, the tail a fifth of a second back along its way.
      const j = i * 6;
      for (let c = 0; c < 3; c++) {
        airPos[j + c] = aP[k + c]!;
        airPos[j + 3 + c] = aP[k + c]! - aV[k + c]! * 0.2;
      }
      airAlpha[2 * i] = a;
      airAlpha[2 * i + 1] = 0;
      if (age >= aLife[i]! || done) aAge[i] = -1;
    }
    (airGeo.getAttribute("position") as THREE.BufferAttribute).needsUpdate = true;
    (airGeo.getAttribute("aAlpha") as THREE.BufferAttribute).needsUpdate = true;
  }
  function placeLine(t: number): void {
    for (const u of units) u.position.z = frontZ;
    const back = frontZ - BROILER.d;
    // The seal comes down from above the broilers' backs onto them, as long as the gap.
    const drop = (1 - smooth(sealK)) * 0.5;
    seal.visible = lip.visible = sealK > 0.01;
    seal.scale.z = Math.max(0.001, back - 0.002);
    seal.position.set(0, SEAL_Y + drop, back / 2);
    lip.position.set(0, SEAL_Y + drop - 0.025, back - 0.0008);
    // The cook a step in front of the middle broiler, reaching to the steaks.
    cook.place(COOK_X, frontZ + 0.44, Math.PI);
    grateAt.set(0.1 + 0.08 * Math.sin(t * 0.9), COOK_Y + 0.02, frontZ - 0.35);
    cook.work("plate", t, { at: grateAt, surfaceY: COOK_Y });
  }

  renderAll();
  frame();
  sizeScale();
  new ResizeObserver(() => {
    sizeScale();
    if (stageEl.clientWidth === fitW && stageEl.clientHeight === fitH) return;
    tagSize.clear();
    frame();
  }).observe(stageEl);
  void document.fonts?.ready.then(() => frame());
  if (chartEl) {
    let lastW = chartEl.clientWidth;
    new ResizeObserver(() => {
      if (Math.abs(chartEl.clientWidth - lastW) < 4) return;
      lastW = chartEl.clientWidth;
      renderChart();
    }).observe(chartEl);
  }
  for (const b of atBtns) b.addEventListener("click", () => pick({ at: Number(b.dataset.hcAt) as SetIndex, sealed: setup.sealed }));
  for (const b of sealBtns) b.addEventListener("click", () => pick({ at: setup.at, sealed: b.dataset.hcSeal === "1" }));
  canvas.addEventListener("pointerdown", () => {
    touring = false;
  });
  fanIn?.addEventListener("input", () => {
    touring = false;
    fan = Math.max(2000, Math.min(6000, Math.round(Number(fanIn.value) / 100) * 100));
    refate();
    renderAll();
  });
  void cook.ready.then(() => {
    cook.hand("R").add(tongs);
    void gate.refresh(scene);
  });

  let inView = false;
  observeInView(stageEl, (v) => {
    inView = v;
  });
  let last = performance.now();
  const fitClip = clipFitter(camera, scene, 0.05);
  let frozen = false;
  function animate(now: number): void {
    requestAnimationFrame(animate);
    const dt = Math.max(0, Math.min(0.05, (now - last) / 1000));
    last = now;
    if (!inView || document.hidden || !gate.ready()) return;
    const t = now / 1000;
    if (touring && !frozen) {
      tourT += dt;
      const k = Math.floor(tourT / TOUR_S) % TOUR.length;
      show(TOUR[k]!);
    }
    // Ease the broilers and the seal: the seal comes out before the broilers move and goes in once they have stopped;
    // they slide at about a foot a second.
    const zTo = frontOf(setup.at);
    const dz = zTo - frontZ;
    const moving = Math.abs(dz) > 1e-4;
    const sGo = moving ? 0 : setup.sealed ? 1 : 0;
    sealK += Math.sign(sGo - sealK) * Math.min(Math.abs(sGo - sealK), dt * 1.4);
    if (moving && sealK <= 0) frontZ += Math.sign(dz) * Math.min(Math.abs(dz), dt * 0.3);
    placeLine(t);
    stepAir(dt, t);
    stepSmoke(dt, t);
    showHud();
    controls.update();
    fitClip();
    renderer.render(scene, camera);
    showTags();
  }
  requestAnimationFrame(animate);
  void gate.refresh(scene);

  // A hook for checking the page by script.
  Object.assign(stageEl, {
    hood: {
      info: () => {
        let escaping = 0;
        let stray = 0;
        let live = 0;
        for (let i = 0; i < NS; i++) {
          if (sAge[i]! < 0) continue;
          live++;
          if (sEsc[i]! >= 0) escaping++;
          if (sEsc[i]! >= 0 && !sOut[i]) stray++;
        }
        return {
          mountMs: Math.round(mountMs),
          running: inView && !document.hidden && gate.ready(),
          inView,
          gateReady: gate.ready(),
          at: setup.at,
          sealed: setup.sealed,
          fan,
          touring,
          need: needCfm(setup),
          caught: catches(fan, setup),
          caughtCount: caughtCount(fan),
          frontIn: +((HOOD.d - frontZ) / IN).toFixed(2),
          sealK: +sealK.toFixed(2),
          smoke: { live, escaping, stray },
          tags: [tagGap, tagLip].map((el) => (el && !el.hidden ? el.textContent : null)),
          hud: hudEl?.textContent ?? "",
          status: statusEl?.textContent ?? "",
        };
      },
      /** Show a setup (stops the tour), and set the fan. */
      set(s: { at?: SetIndex; sealed?: boolean; fan?: number }) {
        if (s.fan !== undefined) {
          fan = Math.max(2000, Math.min(6000, Math.round(s.fan / 100) * 100));
          refate();
        }
        if (s.at !== undefined || s.sealed !== undefined) pick({ at: s.at ?? setup.at, sealed: s.sealed ?? setup.sealed });
        renderAll();
      },
      /** Jump the broilers and the seal to where the setup has them, at once. */
      settle() {
        frontZ = frontOf(setup.at);
        sealK = setup.sealed ? 1 : 0;
      },
      /** Hold the tour where it is (true), or let it run. */
      hold(on: boolean) {
        frozen = on;
      },
      view(eye: [number, number, number] | null, at?: [number, number, number]) {
        userMoved = !!eye;
        if (eye && at) {
          camera.clearViewOffset();
          camera.position.set(...eye);
          controls.target.set(...at);
          controls.update();
        } else frame();
      },
    },
  });
  const mountMs = performance.now() - mountT0;
}
