// Kitty Hawk scene dressing, three.js layer (bead guzez.13). THIN
// consumer of the tested math in dressing.ts: this file only builds
// meshes and applies poses. Procedural canvas textures keep the app
// self-contained (no asset fetches; the CSP and the size budget both
// stay honest).

import * as THREE from "three";
import {
  DEFAULT_HEADWIND_MPS,
  campLayout,
  emberAt,
  exhaustPuff,
  flagPoint,
  gullAttitude,
  gullFleet,
  gullPose,
  landingDust,
  lcg,
  orvillePose,
  hash01,
  propwashPuff,
  railTies,
  scrubField,
  smokePuff,
  streamerPoint,
  flybyFleet,
  flybyPose,
  type FlybyPath,
  type GullPath,
} from "./dressing.ts";
import { createBrotherFigure } from "./figure3d.ts";
import { strideFreqHz } from "./figure.ts";
import {
  SKY_DOME_FRAG_GLSL,
  SKY_DOME_VERT_GLSL,
  groundHazeColor,
  horizonColor,
  sunDiscColor,
  zenithColor,
} from "./sky/atmosphere.ts";

/* ------------------------------ the sun ---------------------------- */
/* ONE shared sun: the sky-texture disc and every light/glint derive
 * from this single direction, so painted sun and lit shading can
 * never disagree. December 17, ~10:35 a.m. solar time, Kill Devil
 * Hills (~36°N): sun low in the SOUTH-EAST, elevation ≈ 28°. Scene
 * frame: x east, z south, y up — SE means POSITIVE x and z. */
export const SUN_DIRECTION: readonly [number, number, number] = [
  0.438, 0.468, 0.767,
];
/** Warm low-winter-sun tint for lights and glints. */
export const SUN_COLOR = 0xffe3bd;

/* ---------- procedural textures (canvas, deterministic) ---------- */

function canvasTexture(
  size: number,
  draw: (ctx: CanvasRenderingContext2D, size: number) => void,
): THREE.CanvasTexture {
  const canvas = document.createElement("canvas");
  canvas.width = size;
  canvas.height = size;
  const ctx = canvas.getContext("2d");
  if (ctx === null) {
    // Headless/degraded: a 1x1 texture keeps materials valid.
    canvas.width = 1;
    canvas.height = 1;
  } else {
    draw(ctx, size);
  }
  const tex = new THREE.CanvasTexture(canvas);
  tex.colorSpace = THREE.SRGBColorSpace;
  return tex;
}

/** A precomputed sand-surface stroke: ONE LCG pass feeds the albedo,
 * the height field, AND the roughness variation, so the three maps
 * describe the same grains and ripples (they can never disagree). */
interface SandStroke {
  y: number;
  width: number;
  bright: boolean;
  c1y: number;
  c2y: number;
}
interface SandSpeck {
  x: number;
  y: number;
  tone: number;
}

/** Deterministic sand surface recipe (seed 1903, same flock law as
 * everything else in this app). */
function sandRecipe(): { strokes: SandStroke[]; specks: SandSpeck[] } {
  const rand = lcg(1903);
  const strokes: SandStroke[] = [];
  for (let i = 0; i < 46; i += 1) {
    strokes.push({
      y: rand(),
      width: 2 + rand() * 5,
      bright: rand() < 0.5,
      c1y: (rand() - 0.5) * 14,
      c2y: (rand() - 0.5) * 14,
    });
  }
  const specks: SandSpeck[] = [];
  for (let i = 0; i < 9000; i += 1) {
    specks.push({ x: rand(), y: rand(), tone: 150 + Math.floor(rand() * 90) });
  }
  return { strokes, specks };
}

const SAND_SRGB = new THREE.Color("#c2b088");

/** Full PBR sand set (T1.2): albedo with wind-ripple banding, a normal
 * map derived from the SAME ripple height field, and a roughness map
 * that darkens wet-looking troughs. All procedural canvases — no asset
 * fetches, deterministic under the seed. */
function sandMaps(): {
  map: THREE.CanvasTexture;
  normalMap: THREE.CanvasTexture;
  roughnessMap: THREE.CanvasTexture;
} {
  const { strokes, specks } = sandRecipe();
  const size = 512;
  // Height field first (grayscale canvas), then numerically derive the
  // tangent-space normals from it.
  const hCanvas = document.createElement("canvas");
  hCanvas.width = size;
  hCanvas.height = size;
  const hctx = hCanvas.getContext("2d")!;
  if (hctx !== null) {
    hctx.fillStyle = "#808080";
    hctx.fillRect(0, 0, size, size);
    for (const st of strokes) {
      // Crests raised (+), troughs sunk (−) — matches the albedo bands.
      hctx.strokeStyle = st.bright ? "rgba(210,210,210,0.55)" : "rgba(70,70,70,0.45)";
      hctx.lineWidth = st.width * 1.6;
      hctx.beginPath();
      hctx.moveTo(0, st.y * size);
      hctx.bezierCurveTo(size * 0.3, st.y * size + st.c1y, size * 0.7, st.y * size + st.c2y, size, st.y * size);
      hctx.stroke();
    }
    for (const sp of specks) {
      hctx.fillStyle = `rgba(${sp.tone},${sp.tone},${sp.tone},0.25)`;
      hctx.fillRect(sp.x * size, sp.y * size, 1.4, 1.4);
    }
  }
  const data = hctx.getImageData(0, 0, size, size).data;
  const nCanvas = document.createElement("canvas");
  nCanvas.width = size;
  nCanvas.height = size;
  const nctx = nCanvas.getContext("2d")!;
  const rCanvas = document.createElement("canvas");
  rCanvas.width = size;
  rCanvas.height = size;
  const rctx = rCanvas.getContext("2d")!;
  const nImg = nctx.createImageData(size, size);
  const rImg = rctx.createImageData(size, size);
  const at = (x: number, y: number): number =>
    data[(((y + size) % size) * size + ((x + size) % size)) * 4]! / 255;
  const strength = 2.2;
  for (let y = 0; y < size; y++) {
    for (let x = 0; x < size; x++) {
      const dxl = at(x - 1, y);
      const dxr = at(x + 1, y);
      const dyu = at(x, y - 1);
      const dyd = at(x, y + 1);
      let nx = (dxl - dxr) * strength;
      let ny = (dyu - dyd) * strength;
      const nz = 1;
      const len = Math.hypot(nx, ny, nz);
      nx /= len;
      ny /= len;
      const o = (y * size + x) * 4;
      nImg.data[o] = Math.round((nx * 0.5 + 0.5) * 255);
      nImg.data[o + 1] = Math.round((ny * 0.5 + 0.5) * 255);
      nImg.data[o + 2] = Math.round((nz / len * 0.5 + 0.5) * 255);
      nImg.data[o + 3] = 255;
      // Roughness: crests slightly rougher/brighter than troughs.
      const h = at(x, y);
      const rough = Math.round((0.82 + 0.16 * (h - 0.5)) * 255);
      rImg.data[o] = rough;
      rImg.data[o + 1] = rough;
      rImg.data[o + 2] = rough;
      rImg.data[o + 3] = 255;
    }
  }
  nctx.putImageData(nImg, 0, 0);
  rctx.putImageData(rImg, 0, 0);
  // Albedo canvas from the same recipe.
  const aCanvas = document.createElement("canvas");
  aCanvas.width = size;
  aCanvas.height = size;
  const actx = aCanvas.getContext("2d")!;
  if (actx !== null) {
    actx.fillStyle = SAND_SRGB.getStyle();
    actx.fillRect(0, 0, size, size);
    for (const st of strokes) {
      actx.strokeStyle = st.bright ? "rgba(232,218,184,0.20)" : "rgba(150,130,96,0.18)";
      actx.lineWidth = st.width;
      actx.beginPath();
      actx.moveTo(0, st.y * size);
      actx.bezierCurveTo(size * 0.3, st.y * size + st.c1y, size * 0.7, st.y * size + st.c2y, size, st.y * size);
      actx.stroke();
    }
    for (const sp of specks) {
      actx.fillStyle = `rgba(${sp.tone},${sp.tone - 16},${sp.tone - 46},0.28)`;
      actx.fillRect(sp.x * size, sp.y * size, 1.4, 1.4);
    }
  }
  const map = new THREE.CanvasTexture(aCanvas);
  map.colorSpace = THREE.SRGBColorSpace;
  const normalMap = new THREE.CanvasTexture(nCanvas);
  const roughnessMap = new THREE.CanvasTexture(rCanvas);
  return { map, normalMap, roughnessMap };
}

/** Apply tiling parameters consistently to a whole sand map set. */
function tileSand(
  maps: { map: THREE.CanvasTexture; normalMap: THREE.CanvasTexture; roughnessMap: THREE.CanvasTexture },
  repeat: number,
): void {
  for (const tex of [maps.map, maps.normalMap, maps.roughnessMap]) {
    tex.wrapS = THREE.RepeatWrapping;
    tex.wrapT = THREE.RepeatWrapping;
    tex.repeat.set(repeat, repeat);
    tex.anisotropy = 8;
  }
}

/** Weathered plank siding for the camp buildings. */
function plankTexture(): THREE.CanvasTexture {
  const rand = lcg(1911);
  return canvasTexture(256, (ctx, s) => {
    ctx.fillStyle = "#8a7354";
    ctx.fillRect(0, 0, s, s);
    const plank = s / 8;
    for (let i = 0; i < 8; i += 1) {
      const shade = 108 + Math.floor(rand() * 40);
      ctx.fillStyle = `rgb(${shade},${Math.round(shade * 0.82)},${Math.round(shade * 0.6)})`;
      ctx.fillRect(0, i * plank, s, plank - 2);
      // Grain streaks per plank.
      for (let g = 0; g < 5; g += 1) {
        ctx.fillStyle = `rgba(60,44,26,${0.08 + rand() * 0.14})`;
        ctx.fillRect(rand() * s, i * plank + rand() * (plank - 4), 20 + rand() * 90, 1 + rand() * 2);
      }
      // Nail heads at the plank ends.
      ctx.fillStyle = "rgba(40,36,30,0.7)";
      ctx.fillRect(6 + rand() * 8, i * plank + plank / 2, 2, 2);
      ctx.fillRect(s - 10 - rand() * 8, i * plank + plank / 2, 2, 2);
    }
  });
}

/** Cloud card textures (A4). `style` picks the deck:
 * - "cumulus": puffs stacked ABOVE a flat base line (the fair-weather
 *   December look — bright tops, cut bottoms), warm-lit from the SE.
 * - "cirrus": long horizontal fiber streaks, very thin. */
function cloudTexture(seed: number, style: "cumulus" | "cirrus"): THREE.CanvasTexture {
  const rand = lcg(seed);
  return canvasTexture(256, (ctx, s) => {
    ctx.clearRect(0, 0, s, s);
    if (style === "cirrus") {
      for (let i = 0; i < 9; i += 1) {
        const y = s * (0.3 + rand() * 0.4);
        const x0 = rand() * s * 0.5;
        const len = s * (0.3 + rand() * 0.55);
        const g = ctx.createLinearGradient(x0, y, x0 + len, y);
        g.addColorStop(0, "rgba(255,253,248,0)");
        g.addColorStop(0.5, `rgba(255,253,248,${0.16 + rand() * 0.2})`);
        g.addColorStop(1, "rgba(255,253,248,0)");
        ctx.strokeStyle = g;
        ctx.lineWidth = 2 + rand() * 5;
        ctx.beginPath();
        ctx.moveTo(x0, y);
        ctx.bezierCurveTo(x0 + len * 0.3, y - 4, x0 + len * 0.7, y + 4, x0 + len, y);
        ctx.stroke();
      }
      return;
    }
    // Cumulus: flat base at ~62% height, puffs only above it; a faint
    // warm tint on the sun-facing (right) side, cool shadow beneath.
    const base = s * 0.62;
    ctx.fillStyle = "rgba(226,222,214,0.5)";
    ctx.fillRect(s * 0.18, base - 4, s * 0.64, 5);
    for (let i = 0; i < 13; i += 1) {
      const x = s * (0.24 + rand() * 0.52);
      const y = base - rand() * s * 0.3;
      const r = s * (0.07 + rand() * 0.14);
      const warm = x > s * 0.55;
      const puff = ctx.createRadialGradient(x, y, r * 0.15, x, y, r);
      puff.addColorStop(0, warm ? "rgba(255,250,238,0.7)" : "rgba(252,251,247,0.66)");
      puff.addColorStop(0.7, "rgba(238,236,230,0.28)");
      puff.addColorStop(1, "rgba(230,228,222,0)");
      ctx.fillStyle = puff;
      ctx.fillRect(x - r, y - r, r * 2, r * 2);
    }
  });
}

/* ---------------------- shared materials ------------------------- */

const WOOD = new THREE.MeshStandardMaterial({ color: 0x6f5231, roughness: 0.9 });
const WOOD_DARK = new THREE.MeshStandardMaterial({ color: 0x4c3a22, roughness: 0.95 });
const IRON = new THREE.MeshStandardMaterial({ color: 0x3a3f45, roughness: 0.6, metalness: 0.55 });
const SUIT = new THREE.MeshStandardMaterial({ color: 0x2d2c33, roughness: 0.9 });
const SKIN = new THREE.MeshStandardMaterial({ color: 0xc99b78, roughness: 0.8 });
const GULL_BODY = new THREE.MeshStandardMaterial({ color: 0xe8e8e4, roughness: 0.85 });
const GULL_WING = new THREE.MeshStandardMaterial({
  color: 0xd6d6d0,
  roughness: 0.85,
  side: THREE.DoubleSide,
});

function box(
  w: number,
  h: number,
  d: number,
  mat: THREE.Material,
  x = 0,
  y = 0,
  z = 0,
): THREE.Mesh {
  const m = new THREE.Mesh(new THREE.BoxGeometry(w, h, d), mat);
  m.position.set(x, y, z);
  return m;
}

/* ------------------------- environment --------------------------- */

/** Sky dome (BackSide sphere) — replaces the flat clear color. The
 * caller lifts it to the SITE's base elevation: the horizon band sits
 * on the dome's equator, and Huffman Prairie is ~250 m above the
 * KDH sea-level datum (absolute-elevation grids). */
export function buildSky(): THREE.Mesh {
  const geo = new THREE.SphereGeometry(2600, 48, 32);
  // Atmospheric-shader dome (A2): every uniform comes from the tested
  // pure math in sky/atmosphere.ts, and the sun direction IS the
  // scene's SUN_DIRECTION — painted light and cast light cannot
  // disagree. HDR disc feeds the bloom pass; ACES tonemaps once.
  const mat = new THREE.ShaderMaterial({
    vertexShader: SKY_DOME_VERT_GLSL,
    fragmentShader: SKY_DOME_FRAG_GLSL,
    uniforms: {
      uSunDirection: { value: new THREE.Vector3(...SUN_DIRECTION) },
      uZenith: { value: new THREE.Vector3(...zenithColor(SUN_DIRECTION)) },
      uHorizon: { value: new THREE.Vector3(...horizonColor(SUN_DIRECTION)) },
      uGroundHaze: { value: new THREE.Vector3(...groundHazeColor(SUN_DIRECTION)) },
      uSunDisc: { value: new THREE.Vector3(...sunDiscColor(SUN_DIRECTION)) },
      uDiscCos: { value: 0.999989 },
      uGlowStrength: { value: 0.55 },
    },
    side: THREE.BackSide,
    depthWrite: false,
    fog: false,
  });
  const sky = new THREE.Mesh(geo, mat);
  sky.renderOrder = -10;
  return sky;
}

export function buildClouds(baseY: number): THREE.Group {
  const group = new THREE.Group();
  const rand = lcg(1908);
  // Three-layer deck (A4): wind shear means higher cards drift faster.
  // Heights are relative to the site's base elevation — at Huffman's
  // ~250 m datum, absolute cloud bases would sit UNDERGROUND.
  const layers: {
    n: number;
    style: "cumulus" | "cirrus";
    yMin: number;
    ySpan: number;
    wMin: number;
    wSpan: number;
    driftMin: number;
    driftSpan: number;
    opacity: number;
  }[] = [
    { n: 6, style: "cirrus", yMin: 330, ySpan: 110, wMin: 420, wSpan: 380, driftMin: 2.6, driftSpan: 1.4, opacity: 0.5 },
    { n: 7, style: "cumulus", yMin: 210, ySpan: 90, wMin: 240, wSpan: 300, driftMin: 1.8, driftSpan: 1.0, opacity: 0.9 },
    { n: 5, style: "cumulus", yMin: 140, ySpan: 50, wMin: 120, wSpan: 140, driftMin: 1.2, driftSpan: 0.6, opacity: 0.75 },
  ];
  let seed = 200;
  for (const L of layers) {
    for (let i = 0; i < L.n; i += 1) {
      const mat = new THREE.MeshBasicMaterial({
        map: cloudTexture(seed++, L.style),
        transparent: true,
        depthWrite: false,
        opacity: L.opacity,
        fog: false,
      });
      const w = L.wMin + rand() * L.wSpan;
      const h = L.style === "cirrus" ? w * 0.22 : w * (L.n === 5 ? 0.5 : 0.42);
      const cloud = new THREE.Mesh(new THREE.PlaneGeometry(w, h), mat);
      mat.side = THREE.DoubleSide;
      cloud.rotation.x = -Math.PI / 2 + 0.35; // near-horizontal, tipped at the camera
      cloud.position.set(0, baseY + L.yMin + rand() * L.ySpan, -900 + rand() * 1800);
      cloud.userData["baseX"] = -900 + rand() * 1800;
      cloud.userData["driftMps"] = L.driftMin + rand() * L.driftSpan;
      group.add(cloud);
    }
  }
  return group;
}
/** Procedural water-normal map: two octaves of smooth value noise
 * converted to tangent-space normals (same derivation law as the sand
 * set). Seeded — replays render the same sea. */
function waterNormalTexture(): THREE.CanvasTexture {
  const size = 256;
  const rand = lcg(1909);
  const lattice = 17;
  const gridVals: number[] = [];
  for (let i = 0; i < lattice * lattice; i += 1) {
    gridVals.push(rand());
  }
  const sample = (u: number, v: number, freqScale: number, offset: number): number => {
    const fx = ((u % 1) + 1) % 1;
    const fy = ((v % 1) + 1) % 1;
    const gx = fx * (lattice - 1) * freqScale + offset;
    const gy = fy * (lattice - 1) * freqScale + offset * 1.7;
    const ix = Math.floor(gx);
    const iy = Math.floor(gy);
    const tx = gx - ix;
    const ty = gy - iy;
    const sx = tx * tx * (3 - 2 * tx);
    const sy = ty * ty * (3 - 2 * ty);
    const g = (a: number, b: number): number =>
      gridVals[((b % lattice) + lattice) % lattice * lattice + (((a % lattice) + lattice) % lattice)]!;
    const a = g(ix, iy);
    const b = g(ix + 1, iy);
    const c = g(ix, iy + 1);
    const d = g(ix + 1, iy + 1);
    return a * (1 - sx) * (1 - sy) + b * sx * (1 - sy) + c * (1 - sx) * sy + d * sx * sy;
  };
  const hCanvas = document.createElement("canvas");
  hCanvas.width = size;
  hCanvas.height = size;
  const hctx = hCanvas.getContext("2d")!;
  const img = hctx.createImageData(size, size);
  for (let y = 0; y < size; y++) {
    for (let x = 0; x < size; x++) {
      const u = x / size;
      const v = y / size;
      const h =
        sample(u, v, 1, 0) * 0.65 +
        sample(u * 2.7, v * 2.7, 1, 31.7) * 0.35;
      const o = (y * size + x) * 4;
      const b = Math.round(h * 255);
      img.data[o] = b;
      img.data[o + 1] = b;
      img.data[o + 2] = b;
      img.data[o + 3] = 255;
    }
  }
  hctx.putImageData(img, 0, 0);
  // Derive tangent-space normals from the height field.
  const data = img.data;
  const nCanvas = document.createElement("canvas");
  nCanvas.width = size;
  nCanvas.height = size;
  const nctx = nCanvas.getContext("2d")!;
  const nImg = nctx.createImageData(size, size);
  const at = (x: number, y: number): number =>
    data[(((y + size) % size) * size + ((x + size) % size)) * 4]! / 255;
  const strength = 1.6;
  for (let y = 0; y < size; y++) {
    for (let x = 0; x < size; x++) {
      const nx = (at(x - 1, y) - at(x + 1, y)) * strength;
      const ny = (at(x, y - 1) - at(x, y + 1)) * strength;
      const nz = 1;
      const len = Math.hypot(nx, ny, nz);
      const o = (y * size + x) * 4;
      nImg.data[o] = Math.round((nx / len * 0.5 + 0.5) * 255);
      nImg.data[o + 1] = Math.round((ny / len * 0.5 + 0.5) * 255);
      nImg.data[o + 2] = Math.round(255 * (nz / len * 0.5 + 0.5));
      nImg.data[o + 3] = 255;
    }
  }
  nctx.putImageData(nImg, 0, 0);
  const tex = new THREE.CanvasTexture(nCanvas);
  tex.wrapS = THREE.RepeatWrapping;
  tex.wrapT = THREE.RepeatWrapping;
  tex.repeat.set(24, 24);
  return tex;
}

/** Shared ocean time uniform — advanced by animateDressing, pure wall
 * clock t (frame-rate never changes the sea). */
export const OCEAN_TIME: { value: number } = { value: 0 };

/** Sand skirt far beyond the surveyed tile so the horizon is never
 * void, plus (Kill Devil Hills only) the Atlantic to the EAST —
 * Huffman Prairie is landlocked Ohio pasture and gets NO ocean.
 * The sea scrolls its own normal map, swells gently in the vertex
 * stage, and carries an oscillating surf-foam line where it meets the
 * beach (T1.4). */
export function buildOuterGround(
  tileExtentM: number,
  withOcean: boolean,
  baseY: number,
): THREE.Group {
  const group = new THREE.Group();
  const sandMapsKdh = sandMaps();
  tileSand(sandMapsKdh, 160);
  const skirt = new THREE.Mesh(
    new THREE.CircleGeometry(2400, 64),
    new THREE.MeshStandardMaterial({
      map: sandMapsKdh.map,
      normalMap: sandMapsKdh.normalMap,
      roughnessMap: sandMapsKdh.roughnessMap,
      color: 0xcabb95,
      roughness: 1,
      normalScale: new THREE.Vector2(0.55, 0.55),
    }),
  );
  skirt.rotation.x = -Math.PI / 2;
  // Tucked just under the tile's LOWEST elevation (grids store
  // absolute heights: KDH ~0 m, Huffman ~242 m).
  skirt.position.y = baseY - 0.35;
  skirt.receiveShadow = true;
  group.add(skirt);
  if (!withOcean) {
    return group;
  }
  const waterNrm = waterNormalTexture();
  const oceanMat = new THREE.MeshStandardMaterial({
    color: 0x28495c,
    roughness: 0.16,
    metalness: 0.06,
    normalMap: waterNrm,
    normalScale: new THREE.Vector2(1.6, 1.6),
    transparent: true,
    opacity: 0.96,
  });
  oceanMat.onBeforeCompile = (shader) => {
    shader.uniforms["uOceanTime"] = OCEAN_TIME;
    shader.vertexShader =
      "uniform float uOceanTime;\n" +
      shader.vertexShader.replace(
        "#include <begin_vertex>",
        [
          "#include <begin_vertex>",
          // Gentle swell: three crossed sines lift the surface — two
          // long east-west trains plus a short cross chop. The plane
          // is rotated flat, so LOCAL z is WORLD up.
          "float swellA = sin(transformed.x * 0.045 + uOceanTime * 1.05) * 0.30;",
          "float swellB = sin(transformed.x * 0.013 - transformed.y * 0.021 + uOceanTime * 0.62) * 0.42;",
          "float swellC = sin(transformed.y * 0.085 + uOceanTime * 1.7) * 0.10;",
          "transformed.z += swellA + swellB + swellC;",
        ].join("\n"),
      );
  };
  const ocean = new THREE.Mesh(new THREE.PlaneGeometry(2600, 5200, 128, 128), oceanMat);
  ocean.rotation.x = -Math.PI / 2;
  ocean.position.set(tileExtentM / 2 + 1300, -0.15, 0);
  group.add(ocean);
  // Surf line: a long thin foam strip hugging the shoreline, breathing
  // with the swell period (presentation-only).
  const foamTex = canvasTexture(128, (ctx, s) => {
    ctx.clearRect(0, 0, s, s);
    const grad = ctx.createLinearGradient(0, 0, s, 0);
    grad.addColorStop(0, "rgba(255,255,255,0)");
    grad.addColorStop(0.55, "rgba(255,255,255,0.85)");
    grad.addColorStop(1, "rgba(240,246,248,0)");
    ctx.fillStyle = grad;
    ctx.fillRect(0, 0, s, s);
    const rand = lcg(1912);
    for (let i = 0; i < 40; i += 1) {
      ctx.fillStyle = `rgba(255,255,255,${0.1 + rand() * 0.25})`;
      ctx.fillRect(rand() * s, rand() * s, 2 + rand() * 6, 1 + rand() * 2);
    }
  });
  foamTex.wrapS = THREE.RepeatWrapping;
  foamTex.wrapT = THREE.RepeatWrapping;
  foamTex.repeat.set(60, 1);
  const foam = new THREE.Mesh(
    new THREE.PlaneGeometry(46, 2600),
    new THREE.MeshBasicMaterial({
      map: foamTex,
      transparent: true,
      depthWrite: false,
      opacity: 0.55,
    }),
  );
  foam.rotation.x = -Math.PI / 2;
  foam.position.set(tileExtentM / 2 + 26, -0.02, 0);
  foam.userData["foam"] = true;
  group.add(foam);
  group.userData["foamMesh"] = foam;
  group.userData["waterNormal"] = waterNrm;
  return group;
}

/** The surveyed tile gets the same PBR sand set — vertex colors still
 * tint water/sand/dune classes under the albedo. */
export function sandTileMaterial(): THREE.MeshStandardMaterial {
  const maps = sandMaps();
  tileSand(maps, 130);
  return new THREE.MeshStandardMaterial({
    map: maps.map,
    normalMap: maps.normalMap,
    roughnessMap: maps.roughnessMap,
    vertexColors: true,
    roughness: 1,
    normalScale: new THREE.Vector2(0.8, 0.8),
  });
}
/* --------------------------- the rail ---------------------------- */

/** The 60 ft monorail: 2x4 on edge, half-buried ties, and the small
 * starting trestle. Positions are launch-relative; caller places. */
export function buildRail(railLengthM: number): THREE.Group {
  const group = new THREE.Group();
  const rail = box(railLengthM, 0.09, 0.04, WOOD_DARK, railLengthM / 2, 0.16, 0);
  group.add(rail);
  const cap = box(railLengthM, 0.015, 0.09, IRON, railLengthM / 2, 0.21, 0);
  group.add(cap);
  for (const x of railTies(railLengthM)) {
    group.add(box(0.35, 0.06, 1.15, WOOD, x, 0.05, 0));
  }
  // Starting trestle: two A-frames + a cross bench at the tail.
  for (const dz of [-0.45, 0.45]) {
    group.add(box(0.07, 0.5, 0.07, WOOD, -0.6, 0.25, dz));
  }
  group.add(box(0.09, 0.07, 1.15, WOOD, -0.6, 0.52, 0));
  return group;
}

/* --------------------------- the camp ---------------------------- */

function buildBuilding(widthM: number, depthM: number, heightM: number): THREE.Group {
  const group = new THREE.Group();
  const planks = plankTexture();
  const wall = new THREE.MeshStandardMaterial({ map: planks, roughness: 0.95 });
  const body = box(widthM, heightM, depthM, wall, 0, heightM / 2, 0);
  group.add(body);
  // Gable roof: two pitched slabs + dark tar paper & wood battens.
  const half = widthM / 2;
  const rise = heightM * 0.45;
  const slabLen = Math.hypot(half, rise) + 0.2;
  for (const side of [-1, 1]) {
    const slab = box(slabLen, 0.06, depthM + 0.4, WOOD_DARK);
    slab.position.set((side * half) / 2, heightM + rise / 2, 0);
    slab.rotation.z = -side * Math.atan2(rise, half);
    group.add(slab);
    // Roof battens across the slope
    for (let b = -depthM / 2; b <= depthM / 2; b += 1.4) {
      const batten = box(slabLen * 0.98, 0.025, 0.04, WOOD_DARK);
      batten.position.set((side * half) / 2, heightM + rise / 2 + 0.04, b);
      batten.rotation.z = -side * Math.atan2(rise, half);
      group.add(batten);
    }
  }
  // Dark doorway on the +x gable end with wood lintel framing.
  const door = box(0.04, heightM * 0.62, 1.15, WOOD_DARK, widthM / 2 + 0.01, heightM * 0.31, 0);
  door.material = new THREE.MeshStandardMaterial({ color: 0x17120b, roughness: 1 });
  group.add(door);
  const frameTop = box(0.08, 0.08, 1.35, WOOD_DARK, widthM / 2 + 0.02, heightM * 0.63, 0);
  group.add(frameTop);
  group.traverse((o) => {
    if ((o as THREE.Mesh).isMesh) {
      o.castShadow = true;
      o.receiveShadow = true;
    }
  });
  return group;
}

export interface Campfire {
  group: THREE.Group;
  light: THREE.PointLight;
  flame: THREE.Mesh;
}

function buildCampfire(): Campfire {
  const group = new THREE.Group();
  const rand = lcg(1917);
  // Stone ring + charred logs + glowing ash bed.
  const ashBed = new THREE.Mesh(
    new THREE.CircleGeometry(0.58, 16),
    new THREE.MeshStandardMaterial({ color: 0x221a14, roughness: 1 }),
  );
  ashBed.rotation.x = -Math.PI / 2;
  ashBed.position.y = 0.01;
  group.add(ashBed);
  for (let i = 0; i < 9; i += 1) {
    const a = (i / 9) * Math.PI * 2;
    const stone = new THREE.Mesh(
      new THREE.DodecahedronGeometry(0.14 + rand() * 0.07),
      new THREE.MeshStandardMaterial({ color: 0x8f8a80, roughness: 1 }),
    );
    stone.position.set(0.62 * Math.cos(a), 0.08, 0.62 * Math.sin(a));
    stone.castShadow = true;
    group.add(stone);
  }
  for (let i = 0; i < 4; i += 1) {
    const log = new THREE.Mesh(new THREE.CylinderGeometry(0.055, 0.07, 0.8, 6), WOOD_DARK);
    log.rotation.set(Math.PI / 2.3, (i / 4) * Math.PI, 0);
    log.position.set(0, 0.12, 0);
    log.castShadow = true;
    group.add(log);
  }
  const flame = new THREE.Mesh(
    new THREE.ConeGeometry(0.16, 0.5, 7),
    new THREE.MeshBasicMaterial({ color: 0xff9a33, transparent: true, opacity: 0.85 }),
  );
  flame.position.y = 0.38;
  group.add(flame);
  const light = new THREE.PointLight(0xff9440, 14, 22, 2);
  light.position.set(0, 0.9, 0);
  light.castShadow = true;
  group.add(light);
  return { group, light, flame };
}

function buildChair(): THREE.Group {
  const g = new THREE.Group();
  g.add(box(0.42, 0.04, 0.42, WOOD, 0, 0.45, 0));
  g.add(box(0.42, 0.5, 0.04, WOOD, 0, 0.73, -0.19));
  for (const [dx, dz] of [
    [-0.18, -0.18],
    [0.18, -0.18],
    [-0.18, 0.18],
    [0.18, 0.18],
  ] as const) {
    g.add(box(0.04, 0.45, 0.04, WOOD, dx, 0.225, dz));
  }
  g.traverse((o) => {
    if ((o as THREE.Mesh).isMesh) {
      o.castShadow = true;
      o.receiveShadow = true;
    }
  });
  return g;
}

function buildBarrel(): THREE.Mesh {
  const geo = new THREE.CylinderGeometry(0.3, 0.26, 0.85, 12);
  geo.translate(0, 0.425, 0); // base-origin: the camp loop seats y = ground
  const b = new THREE.Mesh(geo, WOOD);
  b.castShadow = true;
  b.receiveShadow = true;
  return b;
}

function buildWorkbench(): THREE.Group {
  const g = new THREE.Group();
  // Heavy timber top & legs
  g.add(box(1.8, 0.07, 0.7, WOOD, 0, 0.85, 0));
  for (const [dx, dz] of [
    [-0.8, -0.28],
    [0.8, -0.28],
    [-0.8, 0.28],
    [0.8, 0.28],
  ] as const) {
    g.add(box(0.07, 0.85, 0.07, WOOD_DARK, dx, 0.42, dz));
  }
  // Bottom tool shelf
  g.add(box(1.6, 0.03, 0.55, WOOD_DARK, 0, 0.25, 0));
  // Tools: hand saw, wood plane, grease pot, hammer, brass oiler
  // 1. Hand saw
  g.add(box(0.55, 0.012, 0.12, IRON, -0.4, 0.89, 0.05));
  g.add(box(0.08, 0.03, 0.14, WOOD_DARK, -0.66, 0.9, 0.05));
  // 2. Wood block plane
  g.add(box(0.24, 0.065, 0.08, WOOD_DARK, 0.15, 0.92, -0.1));
  const blade = box(0.02, 0.07, 0.06, IRON, 0.15, 0.94, -0.1);
  blade.rotation.z = -0.45;
  g.add(blade);
  // 3. Oil / grease can
  const can = new THREE.Mesh(new THREE.CylinderGeometry(0.07, 0.07, 0.16, 10), IRON);
  can.position.set(0.6, 0.97, 0.15);
  g.add(can);
  // 4. Spruce wood shavings curled on the bench
  const SHAVING = new THREE.MeshStandardMaterial({ color: 0xd9c59a, roughness: 0.9 });
  for (let i = 0; i < 4; i += 1) {
    const curl = new THREE.Mesh(new THREE.TorusGeometry(0.025, 0.008, 4, 8, Math.PI * 1.5), SHAVING);
    curl.rotation.set(Math.PI / 2, i * 0.8, 0);
    curl.position.set(-0.05 + i * 0.06, 0.89, 0.08 + (i % 2) * 0.04);
    g.add(curl);
  }
  g.traverse((o) => {
    if ((o as THREE.Mesh).isMesh) {
      o.castShadow = true;
      o.receiveShadow = true;
    }
  });
  return g;
}

function buildToolchest(): THREE.Group {
  const g = new THREE.Group();
  g.add(box(0.9, 0.4, 0.5, WOOD_DARK, 0, 0.2, 0));
  g.add(box(0.94, 0.05, 0.54, WOOD, 0, 0.44, 0));
  g.traverse((o) => {
    if ((o as THREE.Mesh).isMesh) {
      o.castShadow = true;
      o.receiveShadow = true;
    }
  });
  return g;
}

/** Packing crate with diagonal corner planks (T2.4). */
function buildCrate(s = 0.55): THREE.Group {
  const g = new THREE.Group();
  g.add(box(s, s, s, WOOD, 0, s / 2, 0));
  for (const off of [-s / 2 + 0.02, s / 2 - 0.02]) {
    const plankA = box(s * 1.42, 0.07, 0.02, WOOD_DARK);
    plankA.rotation.set(Math.PI / 2, 0, Math.PI / 4);
    plankA.position.set(off * 0.7, s / 2, off);
    g.add(plankA);
    const plankB = box(0.07, s * 1.42, 0.02, WOOD_DARK);
    plankB.rotation.z = Math.PI / 4;
    plankB.position.set(off, s / 2, off * 0.7);
    g.add(plankB);
  }
  return g;
}

/** The magneto/battery box: dark case, glass jar cells, and a copper
 * coil — the engine's ignition rig lived on a crate like this. */
function buildBatteryBox(): THREE.Group {
  const g = new THREE.Group();
  const CASE = new THREE.MeshStandardMaterial({ color: 0x241a10, roughness: 0.9 });
  g.add(box(0.62, 0.4, 0.44, CASE, 0, 0.28, 0));
  const GLASS = new THREE.MeshStandardMaterial({
    color: 0xb8c4b0,
    roughness: 0.15,
    metalness: 0.1,
    transparent: true,
    opacity: 0.55,
  });
  for (const dx of [-0.18, 0, 0.18]) {
    const jar = new THREE.Mesh(new THREE.CylinderGeometry(0.06, 0.06, 0.16, 10), GLASS);
    jar.position.set(dx, 0.54, 0);
    g.add(jar);
  }
  const coil = new THREE.Mesh(new THREE.TorusGeometry(0.09, 0.022, 6, 14), IRON);
  coil.rotation.x = Math.PI / 2;
  coil.position.set(0, 0.66, 0);
  g.add(coil);
  return g;
}

/** The classic thumb-pump oil can. */
function buildOilCan(): THREE.Group {
  const g = new THREE.Group();
  const body = new THREE.Mesh(new THREE.SphereGeometry(0.11, 10, 8), IRON);
  body.scale.y = 0.85;
  body.position.y = 0.11;
  g.add(body);
  const spout = new THREE.Mesh(new THREE.CylinderGeometry(0.012, 0.02, 0.3, 6), IRON);
  spout.rotation.z = -0.9;
  spout.position.set(0.12, 0.2, 0);
  g.add(spout);
  return g;
}

/** Hand-anemometer post: four cups on a cross at the top. */
function buildWindPost(): THREE.Group {
  const g = new THREE.Group();
  const pole = new THREE.Mesh(new THREE.CylinderGeometry(0.03, 0.04, 2.2, 8), WOOD_DARK);
  pole.position.y = 1.1;
  pole.castShadow = true;
  g.add(pole);
  const cross = new THREE.Group();
  cross.position.y = 2.25;
  for (let i = 0; i < 4; i += 1) {
    const a = (i / 4) * Math.PI * 2;
    const arm = new THREE.Mesh(new THREE.CylinderGeometry(0.008, 0.008, 0.26, 5), IRON);
    arm.rotation.z = Math.PI / 2;
    arm.rotation.y = a;
    arm.position.set(Math.cos(a) * 0.13, 0, -Math.sin(a) * 0.13);
    cross.add(arm);
    const cup = new THREE.Mesh(
      new THREE.SphereGeometry(0.045, 8, 6, 0, Math.PI * 2, 0, Math.PI / 2),
      IRON,
    );
    cup.position.set(Math.cos(a) * 0.26, 0, -Math.sin(a) * 0.26);
    cup.rotation.y = a;
    cross.add(cup);
  }
  g.add(cross);
  return g;
}

/* -------------------------- the people --------------------------- */

export interface Figure {
  group: THREE.Group;
  leftLeg: THREE.Mesh;
  rightLeg: THREE.Mesh;
  leftArm: THREE.Group;
  rightArm: THREE.Group;
  glasses: THREE.Mesh;
}

/** A standing 1900s figure (dark suit, flat cap) with poseable limbs
 * — Orville on the ground. ~1.78 m tall. */
export function buildFigure(): Figure {
  const group = new THREE.Group();
  const torso = box(0.34, 0.62, 0.22, SUIT, 0, 1.18, 0);
  group.add(torso);
  const head = new THREE.Mesh(new THREE.SphereGeometry(0.115, 12, 10), SKIN);
  head.position.set(0, 1.66, 0);
  group.add(head);
  const cap = new THREE.Mesh(new THREE.CylinderGeometry(0.13, 0.135, 0.055, 12), SUIT);
  cap.position.set(0, 1.755, 0);
  group.add(cap);
  const brim = box(0.12, 0.015, 0.1, SUIT, 0, 1.73, 0.13);
  group.add(brim);
  const mkLeg = (side: number): THREE.Mesh => {
    const leg = box(0.11, 0.85, 0.13, SUIT, side * 0.09, 0.445, 0);
    leg.geometry.translate(0, -0.36, 0);
    leg.position.y = 0.85;
    return leg;
  };
  const leftLeg = mkLeg(-1);
  const rightLeg = mkLeg(1);
  group.add(leftLeg, rightLeg);
  const mkArm = (side: number): THREE.Group => {
    const arm = new THREE.Group();
    arm.position.set(side * 0.22, 1.44, 0);
    const upper = box(0.09, 0.52, 0.1, SUIT, 0, -0.24, 0);
    arm.add(upper);
    const hand = new THREE.Mesh(new THREE.SphereGeometry(0.05, 8, 6), SKIN);
    hand.position.set(0, -0.52, 0);
    arm.add(hand);
    return arm;
  };
  const leftArm = mkArm(-1);
  const rightArm = mkArm(1);
  group.add(leftArm, rightArm);
  // Field glasses: hidden until the pose raises them.
  const glasses = box(0.14, 0.05, 0.07, IRON, 0, 1.62, 0.16);
  glasses.visible = false;
  group.add(glasses);
  return { group, leftLeg, rightLeg, leftArm, rightArm, glasses };
}

/** Wilbur prone on the lower wing (the pilot the player embodies —
 * visible in chase/wingtip/Daniels views). Lies along +x, head fore. */
export function buildProneWilbur(): THREE.Group {
  const group = new THREE.Group();
  const torso = box(0.62, 0.16, 0.3, SUIT, 0, 0.08, 0);
  group.add(torso);
  const head = new THREE.Mesh(new THREE.SphereGeometry(0.1, 12, 10), SKIN);
  head.position.set(0.4, 0.12, 0);
  group.add(head);
  const cap = new THREE.Mesh(new THREE.CylinderGeometry(0.11, 0.115, 0.045, 12), SUIT);
  cap.position.set(0.4, 0.2, 0);
  group.add(cap);
  for (const side of [-1, 1]) {
    const leg = box(0.5, 0.11, 0.12, SUIT, -0.55, 0.06, side * 0.08);
    group.add(leg);
    const arm = box(0.34, 0.09, 0.09, SUIT, 0.28, 0.06, side * 0.2);
    arm.rotation.y = side * 0.25;
    group.add(arm);
  }
  return group;
}

/* ----------------------------- gulls ----------------------------- */

interface GullRig {
  group: THREE.Group;
  leftWing: THREE.Mesh;
  rightWing: THREE.Mesh;
  path: GullPath;
}

interface GullMesh {
  group: THREE.Group;
  leftWing: THREE.Mesh;
  rightWing: THREE.Mesh;
}

/** Herring-gull silhouette (T3.1): tapered curved wings with dark
 * tips, a forked tail, and a body the wings root into. Wings bake flat
 * in the xz-plane with the ROOT at the body, span along ±z — the
 * mesh's rotation.x stays a clean flap hinge. */
function buildGullMesh(): GullMesh {
  const group = new THREE.Group();
  const body = new THREE.Mesh(new THREE.CapsuleGeometry(0.085, 0.34, 3, 6), GULL_BODY);
  body.rotation.z = Math.PI / 2;
  group.add(body);
  const head = new THREE.Mesh(new THREE.SphereGeometry(0.062, 8, 6), GULL_BODY);
  head.position.set(0.26, 0.015, 0);
  group.add(head);
  const beak = new THREE.Mesh(
    new THREE.ConeGeometry(0.024, 0.11, 6),
    new THREE.MeshStandardMaterial({ color: 0xd9a13a, roughness: 0.7 }),
  );
  beak.rotation.z = -Math.PI / 2;
  beak.position.set(0.35, 0.01, 0);
  group.add(beak);
  const TIP = new THREE.MeshStandardMaterial({
    color: 0x4a4a48,
    roughness: 0.85,
    side: THREE.DoubleSide,
  });
  const mkWing = (side: number): THREE.Mesh => {
    // Tapered wing: inner panel (white) + outer panel (grey, dark tip)
    // built as one shaped BufferGeometry — chord shrinks toward the tip,
    // slight camber bows the panel down at rest.
    const segs = 6;
    const positions: number[] = [];
    const idx: number[] = [];
    for (let s = 0; s <= segs; s += 1) {
      const u = s / segs;
      const spanZ = side * (0.06 + u * 0.52);
      const chord = 0.34 * (1 - 0.55 * u) + 0.03;
      const camberY = -Math.sin(u * Math.PI) * 0.035;
      positions.push(0.1 - u * 0.06, camberY, spanZ - (chord / 2) * 0.9);
      positions.push(0.1 + u * 0.02, camberY, spanZ + (chord / 2) * 1.1);
    }
    for (let s = 0; s < segs; s += 1) {
      const a = s * 2;
      idx.push(a, a + 2, a + 1, a + 1, a + 2, a + 3);
    }
    const geo = new THREE.BufferGeometry();
    geo.setAttribute("position", new THREE.Float32BufferAttribute(positions, 3));
    geo.setIndex(idx);
    geo.computeVertexNormals();
    const wing = new THREE.Mesh(geo, GULL_WING);
    // Dark tip: a small overlay quad at the outer 18% of the span.
    const tip = new THREE.Mesh(new THREE.PlaneGeometry(0.16, 0.13), TIP);
    tip.rotation.x = -Math.PI / 2;
    tip.position.set(0.06, -0.012, side * 0.55);
    wing.add(tip);
    return wing;
  };
  const leftWing = mkWing(-1);
  const rightWing = mkWing(1);
  group.add(leftWing, rightWing);
  // Forked tail.
  const tail = new THREE.Mesh(new THREE.ConeGeometry(0.07, 0.22, 4), GULL_WING);
  tail.rotation.set(Math.PI / 2, 0, Math.PI / 2);
  tail.scale.set(1, 1, 0.4);
  tail.position.set(-0.26, 0.005, 0);
  group.add(tail);
  group.traverse((o) => {
    if ((o as THREE.Mesh).isMesh) {
      o.castShadow = true;
    }
  });
  return { group, leftWing, rightWing };
}

/** The 1903 takeoff dolly (T2.1): two wheeled axles with a plank
 * platform the skids rested on. The machine rides it on the rail and
 * it STAYS BEHIND at liftoff (the scene latches the drop). */
export function buildTakeoffDolly(): THREE.Group {
  const group = new THREE.Group();
  for (const ax of [-0.95, 0.95]) {
    const axle = new THREE.Mesh(new THREE.CylinderGeometry(0.028, 0.028, 1.86, 8), IRON);
    axle.rotation.x = Math.PI / 2;
    axle.position.set(ax, 0.3, 0);
    group.add(axle);
    for (const side of [-0.9, 0.9]) {
      const wheel = new THREE.Mesh(new THREE.TorusGeometry(0.27, 0.03, 8, 20), WOOD_DARK);
      wheel.position.set(ax, 0.3, side);
      group.add(wheel);
      const hub = new THREE.Mesh(new THREE.CylinderGeometry(0.05, 0.05, 0.05, 10), IRON);
      hub.rotation.x = Math.PI / 2;
      hub.position.set(ax, 0.3, side);
      group.add(hub);
      for (let s = 0; s < 4; s += 1) {
        const spoke = new THREE.Mesh(new THREE.CylinderGeometry(0.008, 0.008, 0.5, 4), IRON);
        // Wheel disc is the x-y plane (axle along z): rotating about z
        // fans the spokes radially; the old x+y Euler left all four
        // superimposed along the axle.
        spoke.rotation.z = (s / 4) * Math.PI;
        spoke.position.set(ax, 0.3, side);
        group.add(spoke);
      }
    }
    // Strut rising to the skid plate.
    const strut = new THREE.Mesh(new THREE.BoxGeometry(0.07, 0.34, 1.8), WOOD);
    strut.position.set(ax, 0.56, 0);
    group.add(strut);
  }
  // Cross planks the skid runners rest on.
  for (const px of [-0.95, 0, 0.95]) {
    const plank = new THREE.Mesh(new THREE.BoxGeometry(0.16, 0.04, 1.9), WOOD);
    plank.position.set(px, 0.75, 0);
    group.add(plank);
  }
  group.traverse((o) => {
    const mesh = o as THREE.Mesh;
    if (mesh.isMesh) {
      mesh.castShadow = true;
    }
  });
  return group;
}

function buildGull(path: GullPath): GullRig {
  return { ...buildGullMesh(), path };
}

/* --------------------- atmosphere & life systems -------------------- */

/** Shared sway clock for vegetation (advanced by animateDressing). */
export const SWAY_TIME: { value: number } = { value: 0 };

function softDotTexture(inner: string, outer: string): THREE.CanvasTexture {
  return canvasTexture(64, (ctx, s) => {
    ctx.clearRect(0, 0, s, s);
    const g = ctx.createRadialGradient(s / 2, s / 2, 1, s / 2, s / 2, s / 2);
    g.addColorStop(0, inner);
    g.addColorStop(1, outer);
    ctx.fillStyle = g;
    ctx.fillRect(0, 0, s, s);
  });
}

interface SpritePool {
  group: THREE.Group;
  set(i: number, x: number, y: number, z: number, scale: number, opacity: number): void;
  hideAll(): void;
}

/** A fixed pool of billboards with per-sprite opacity (cloned
 * materials — presentation-only counts, ~dozens). */
function makeSpritePool(n: number, tex: THREE.CanvasTexture, baseScale: number): SpritePool {
  const group = new THREE.Group();
  const items: { spr: THREE.Sprite; mat: THREE.SpriteMaterial }[] = [];
  for (let i = 0; i < n; i += 1) {
    const mat = new THREE.SpriteMaterial({
      map: tex,
      transparent: true,
      opacity: 0,
      depthWrite: false,
    });
    const spr = new THREE.Sprite(mat);
    spr.scale.setScalar(baseScale);
    spr.visible = false;
    group.add(spr);
    items.push({ spr, mat });
  }
  return {
    group,
    set(i, x, y, z, scale, opacity) {
      const it = items[i]!;
      if (opacity <= 0.004) {
        it.spr.visible = false;
        return;
      }
      it.spr.visible = true;
      it.spr.position.set(x, y, z);
      it.spr.scale.setScalar(scale);
      it.mat.opacity = opacity;
    },
    hideAll() {
      for (const it of items) {
        it.spr.visible = false;
      }
    },
  };
}

/** Deterministic scrub field meshes (T1.3): instanced tufts/bushes/
 * pines seated on the sampled terrain, tufts and bushes swaying with
 * SWAY_TIME (presentation-only vertex bend). */
function buildVegetation(
  scrub: readonly { x: number; z: number; rotY: number; scale: number; kind: "tuft" | "bush" | "pine" }[],
  groundY: (xRel: number, zRel: number) => number,
): THREE.Group {
  const group = new THREE.Group();
  const grassMat = new THREE.MeshStandardMaterial({ color: 0x96854f, roughness: 1 });
  const bushMat = new THREE.MeshStandardMaterial({ color: 0x5f6136, roughness: 1, flatShading: true });
  const pineMat = new THREE.MeshStandardMaterial({ color: 0x39543a, roughness: 1, flatShading: true });
  const trunkMat = new THREE.MeshStandardMaterial({ color: 0x4c3a24, roughness: 1 });
  const addSway = (mat: THREE.MeshStandardMaterial): void => {
    mat.onBeforeCompile = (shader) => {
      shader.uniforms["uSway"] = SWAY_TIME;
      shader.vertexShader =
        "uniform float uSway;\n" +
        shader.vertexShader.replace(
          "#include <begin_vertex>",
          [
            "#include <begin_vertex>",
            "#ifdef USE_INSTANCING",
            // Bend proportional to local height; phase varies per instance.
            "float ph = instanceMatrix[3][0] * 0.53 + instanceMatrix[3][2] * 0.41;",
            "transformed.x += sin(uSway * 1.8 + ph) * 0.06 * max(transformed.y, 0.0);",
            "#endif",
          ].join("\n"),
        );
    };
  };
  addSway(grassMat);
  addSway(bushMat);
  const tuftGeo = new THREE.ConeGeometry(0.17, 0.52, 5);
  tuftGeo.translate(0, 0.26, 0);
  const bushGeo = new THREE.IcosahedronGeometry(0.44, 0);
  bushGeo.scale(1, 0.72, 1);
  bushGeo.translate(0, 0.32, 0);
  const pineGeo = new THREE.ConeGeometry(0.62, 1.9, 7);
  pineGeo.translate(0, 1.25, 0);
  const trunkGeo = new THREE.CylinderGeometry(0.07, 0.1, 0.6, 6);
  trunkGeo.translate(0, 0.3, 0);
  const nTuft = scrub.filter((p) => p.kind === "tuft").length;
  const nBush = scrub.filter((p) => p.kind === "bush").length;
  const nPine = scrub.filter((p) => p.kind === "pine").length;
  const tufts = new THREE.InstancedMesh(tuftGeo, grassMat, Math.max(1, nTuft));
  const bushes = new THREE.InstancedMesh(bushGeo, bushMat, Math.max(1, nBush));
  const pines = new THREE.InstancedMesh(pineGeo, pineMat, Math.max(1, nPine));
  const trunks = new THREE.InstancedMesh(trunkGeo, trunkMat, Math.max(1, nPine));
  tufts.castShadow = true;
  bushes.castShadow = true;
  pines.castShadow = true;
  const m = new THREE.Matrix4();
  const q = new THREE.Quaternion();
  const eu = new THREE.Euler();
  const one = new THREE.Vector3(1, 1, 1);
  const pos = new THREE.Vector3();
  let ti = 0;
  let bi = 0;
  let pi = 0;
  for (const p of scrub) {
    pos.set(p.x, groundY(p.x, p.z), p.z);
    q.setFromEuler(eu.set(0, p.rotY, 0));
    one.set(p.scale, p.scale, p.scale);
    m.compose(pos, q, one);
    if (p.kind === "tuft") {
      tufts.setMatrixAt(ti++, m);
    } else if (p.kind === "bush") {
      bushes.setMatrixAt(bi++, m);
    } else {
      pines.setMatrixAt(pi++, m);
      trunks.setMatrixAt(pi - 1, m);
    }
  }
  tufts.count = nTuft;
  bushes.count = nBush;
  pines.count = nPine;
  trunks.count = nPine;
  group.add(tufts, bushes, pines, trunks);
  return group;
}

/** Ground detail scatter (A5): instanced shells + pebbles seeded by
 * hash01 — the close-up cue that separates "textured plane" from
 * "beach". One draw call per kind; static geometry, no animation.
 * Deterministic: same seed -> same beach. */
function buildScatter(
  counts: { shells: number; pebbles: number },
  groundY: (xRel: number, zRel: number) => number,
): THREE.Group {
  const group = new THREE.Group();
  const shellMat = new THREE.MeshStandardMaterial({
    color: 0xe8ddc8,
    roughness: 0.55,
    flatShading: true,
  });
  const pebbleMat = new THREE.MeshStandardMaterial({
    color: 0xb3a284,
    roughness: 0.9,
    flatShading: true,
  });
  const shellGeo = new THREE.ConeGeometry(0.055, 0.09, 6);
  const pebbleGeo = new THREE.IcosahedronGeometry(0.04, 0);
  const shells = new THREE.InstancedMesh(shellGeo, shellMat, Math.max(1, counts.shells));
  const pebbles = new THREE.InstancedMesh(pebbleGeo, pebbleMat, Math.max(1, counts.pebbles));
  shells.receiveShadow = true;
  pebbles.receiveShadow = true;
  const m = new THREE.Matrix4();
  const q = new THREE.Quaternion();
  const eu = new THREE.Euler();
  const pos = new THREE.Vector3();
  const one = new THREE.Vector3();
  let si = 0;
  let pi = 0;
  // Ring-rejection placement outside the launch flat (the corridor
  // floor stays pristine for rail/figures), inside ~420 m.
  for (let n = 0; n < counts.shells + counts.pebbles; n += 1) {
    const a = hash01(n * 2 + 1) * Math.PI * 2;
    const r = 62 + Math.pow(hash01(n * 2 + 2), 0.7) * 358;
    const x = Math.cos(a) * r;
    const z = Math.sin(a) * r;
    pos.set(x, groundY(x, z) - 0.01, z);
    q.setFromEuler(eu.set(hash01(n * 3 + 5) * 0.5 - 0.25, hash01(n * 3 + 6) * Math.PI * 2, hash01(n * 3 + 7) * 0.5 - 0.25));
    const s = 0.6 + hash01(n * 5 + 9) * 0.9;
    one.set(s, s, s);
    m.compose(pos, q, one);
    if (n % 2 === 0 && si < counts.shells) {
      shells.setMatrixAt(si++, m);
    } else if (pi < counts.pebbles) {
      pebbles.setMatrixAt(pi++, m);
    }
  }
  shells.count = si;
  pebbles.count = pi;
  group.add(shells, pebbles);
  return group;
}

/** A downwind ribbon: `segs` spans × 2 columns of vertices whose
 * positions the animate loop rewrites from the pure math each frame. */
function buildRibbon(segs: number, halfWidthM: number, color: number, opacity: number): {
  mesh: THREE.Mesh;
  write: (pts: readonly { x: number; y: number; z: number }[]) => void;
} {
  const geo = new THREE.BufferGeometry();
  const positions = new Float32Array((segs + 1) * 2 * 3);
  geo.setAttribute("position", new THREE.BufferAttribute(positions, 3));
  const idx: number[] = [];
  for (let s = 0; s < segs; s += 1) {
    const a = s * 2;
    idx.push(a, a + 2, a + 1, a + 1, a + 2, a + 3);
  }
  geo.setIndex(idx);
  const mesh = new THREE.Mesh(
    geo,
    new THREE.MeshBasicMaterial({
      color,
      transparent: true,
      opacity,
      depthWrite: false,
      side: THREE.DoubleSide,
    }),
  );
  mesh.frustumCulled = false;
  return {
    mesh,
    write(pts) {
      for (let s = 0; s <= segs; s += 1) {
        const p = pts[Math.min(s, pts.length - 1)]!;
        const o = s * 6;
        positions[o] = p.x;
        positions[o + 1] = p.y;
        positions[o + 2] = p.z - halfWidthM;
        positions[o + 3] = p.x;
        positions[o + 4] = p.y;
        positions[o + 5] = p.z + halfWidthM;
      }
      (geo.getAttribute("position") as THREE.BufferAttribute).needsUpdate = true;
    },
  };
}

/** John T. Daniels and his Korona view camera at the tripod (T2.3):
 * static period figure in shirtsleeves + suspenders and a boater hat,
 * hand on the shutter. Returns the flash sprite + light the animate
 * loop pulses when the machine passes. */
function buildDaniels(): { group: THREE.Group; flash: THREE.Sprite; lamp: THREE.PointLight } {
  const group = new THREE.Group();
  const SHIRT = new THREE.MeshStandardMaterial({ color: 0xd8d2c0, roughness: 0.92 });
  const SUSPENDER = new THREE.MeshStandardMaterial({ color: 0x3a3226, roughness: 0.95 });
  const TROUSERS = new THREE.MeshStandardMaterial({ color: 0x59534a, roughness: 0.95 });
  const HAT = new THREE.MeshStandardMaterial({ color: 0xcfc4a2, roughness: 0.9 });
  // Legs + torso.
  for (const dz of [-0.11, 0.11]) {
    const leg = new THREE.Mesh(new THREE.CylinderGeometry(0.075, 0.065, 0.86, 8), TROUSERS);
    leg.position.set(0, 0.43, dz);
    group.add(leg);
  }
  const torso = new THREE.Mesh(new THREE.CapsuleGeometry(0.17, 0.42, 4, 10), SHIRT);
  torso.position.y = 1.18;
  torso.scale.set(1.15, 1, 0.8);
  group.add(torso);
  for (const dz of [-0.09, 0.09]) {
    const strap = new THREE.Mesh(new THREE.BoxGeometry(0.045, 0.46, 0.02), SUSPENDER);
    strap.position.set(0, 1.26, dz * 1.6);
    strap.rotation.x = 0.08;
    group.add(strap);
  }
  // Head + boater hat.
  const head = new THREE.Mesh(new THREE.SphereGeometry(0.115, 12, 10), SKIN);
  head.position.y = 1.63;
  group.add(head);
  const crown = new THREE.Mesh(new THREE.CylinderGeometry(0.115, 0.125, 0.09, 12), HAT);
  crown.position.y = 1.76;
  group.add(crown);
  const brim = new THREE.Mesh(new THREE.CylinderGeometry(0.19, 0.19, 0.014, 14), HAT);
  brim.position.y = 1.72;
  group.add(brim);
  // Arms: right hand low on the shutter, left cradling the lens bed.
  const armR = new THREE.Mesh(new THREE.CylinderGeometry(0.048, 0.042, 0.5, 8), SHIRT);
  armR.position.set(-0.02, 1.28, 0.3);
  armR.rotation.set(-1.15, 0, 0.15);
  group.add(armR);
  const armL = new THREE.Mesh(new THREE.CylinderGeometry(0.048, 0.042, 0.48, 8), SHIRT);
  armL.position.set(0.05, 1.24, -0.28);
  armL.rotation.set(1.05, 0, -0.12);
  group.add(armL);
  // Korona view camera on a wooden tripod, aimed WEST (-x) at the rail.
  const CAM = new THREE.MeshStandardMaterial({ color: 0x33241a, roughness: 0.7 });
  for (const a of [0.4, 2.5, 4.6]) {
    const leg = new THREE.Mesh(new THREE.CylinderGeometry(0.022, 0.018, 1.3, 6), WOOD_DARK);
    leg.position.set(Math.sin(a) * 0.3, 0.65, Math.cos(a) * 0.3);
    leg.rotation.z = Math.sin(a) * 0.24;
    leg.rotation.x = -Math.cos(a) * 0.24;
    group.add(leg);
  }
  const bodyBox = new THREE.Mesh(new THREE.BoxGeometry(0.34, 0.3, 0.42), CAM);
  bodyBox.position.set(0, 1.42, 0);
  group.add(bodyBox);
  const lensBoard = new THREE.Mesh(new THREE.BoxGeometry(0.03, 0.2, 0.2), IRON);
  lensBoard.position.set(-0.2, 1.42, 0);
  group.add(lensBoard);
  const bellows = new THREE.Mesh(new THREE.BoxGeometry(0.1, 0.22, 0.3), SUSPENDER);
  bellows.position.set(0.1, 1.42, 0);
  group.add(bellows);
  // The flash: a bright billboard + short-lived point light.
  const flashTex = softDotTexture("rgba(255,252,238,0.98)", "rgba(255,240,200,0)");
  const flashMat = new THREE.SpriteMaterial({
    map: flashTex,
    transparent: true,
    opacity: 0,
    depthWrite: false,
  });
  const flash = new THREE.Sprite(flashMat);
  flash.scale.setScalar(1.6);
  flash.position.set(-0.3, 1.5, 0);
  group.add(flash);
  const lamp = new THREE.PointLight(0xfff2cc, 0, 30, 2);
  lamp.position.set(-0.4, 1.6, 0);
  group.add(lamp);
  return { group, flash, lamp };
}

/* ------------------------ assembled diorama ---------------------- */

export interface Dressing {
  group: THREE.Group;
  /** Advance every animated element to scene time t. `machine` is
   * supplied while a real/simulated aircraft state exists (attract
   * mode omits it → no plumes); rpm01/gust01 are clamped [0,1]. */
  animate(
    t: number,
    orville: {
      onRail: boolean;
      aircraftX: number;
      releaseX: number | null;
      releaseT: number | null;
      /** Set once the machine is DOWN (ground contact): he runs to it. */
      landedX?: number | null;
      machine?: { x: number; rpm01: number; gust01: number; dustT?: number };
    },
  ): void;
  /** QoS presentation budget (T0.6): 0 = full, 1 = secondary plumes
   * off, 2 = every particle system hidden. Presentation-only. */
  setParticleLevel(level: 0 | 1 | 2): void;
  /** Orville's world position at scene time t under the same state the
   * animate loop uses (pure recomputation of orvillePose + terrain) —
   * the camera glances AT him without reaching into the rig. Eye-height
   * offset included. */
  orvillePosition(
    t: number,
    orville: {
      onRail: boolean;
      aircraftX: number;
      releaseX: number | null;
      releaseT: number | null;
      landedX?: number | null;
    },
  ): [number, number, number];
}

/** Named contract for the camp flagpole rig. */
export interface Flagpole {
  group: THREE.Group;
  writeFlag: (pts: readonly { x: number; y: number; z: number }[]) => void;
}

/** Camp flagpole: 6 m pole, truck ball, and a cream weathered flag
 * ribbon whose cloth the animate loop rewrites from `flagPoint`. */
function buildFlagpole(): Flagpole {
  const group = new THREE.Group();
  const pole = new THREE.Mesh(new THREE.CylinderGeometry(0.035, 0.05, 6, 8), WOOD_DARK);
  pole.position.y = 3;
  pole.castShadow = true;
  group.add(pole);
  const truck = new THREE.Mesh(new THREE.SphereGeometry(0.06, 8, 6), IRON);
  truck.position.y = 6.03;
  group.add(truck);
  // Cloth hangs from just under the truck, streaming downwind (-x).
  const segs = 8;
  const geo = new THREE.BufferGeometry();
  const positions = new Float32Array((segs + 1) * 2 * 3);
  geo.setAttribute("position", new THREE.BufferAttribute(positions, 3));
  const idx: number[] = [];
  for (let s = 0; s < segs; s += 1) {
    const a = s * 2;
    idx.push(a, a + 2, a + 1, a + 1, a + 2, a + 3);
  }
  geo.setIndex(idx);
  const cloth = new THREE.Mesh(
    geo,
    new THREE.MeshStandardMaterial({
      color: 0xd8cfae,
      side: THREE.DoubleSide,
      roughness: 0.95,
    }),
  );
  cloth.castShadow = true;
  cloth.frustumCulled = false;
  group.add(cloth);
  return {
    group,
    writeFlag(pts) {
      for (let s = 0; s <= segs; s += 1) {
        const p = pts[Math.min(s, pts.length - 1)]!;
        const o = s * 6;
        positions[o] = p.x;
        positions[o + 1] = 5.9 + p.y;
        positions[o + 2] = p.z - 0.28;
        positions[o + 3] = p.x;
        positions[o + 4] = 5.9 + p.y;
        positions[o + 5] = p.z + 0.28;
      }
      (geo.getAttribute("position") as THREE.BufferAttribute).needsUpdate = true;
    },
  };
}
/** Build the full diorama around the launch point. `groundY` samples
 * terrain height for prop placement (launch-relative x/z in metres).
 * `windMps` is the scenario headwind driving every wind-made-visible
 * system; `withLife` gates the vegetation/particle systems for the QoS
 * Critical tier and headless tests. */
export function buildDressing(
  launch: readonly [number, number, number],
  railLengthM: number,
  tileExtentM: number,
  withOcean: boolean,
  baseY: number,
  groundY: (xRel: number, zRel: number) => number,
  windMps: number = DEFAULT_HEADWIND_MPS,
  withLife: boolean = true,
): Dressing {
  const group = new THREE.Group();
  const sky = buildSky();
  sky.position.y = baseY;
  group.add(sky);
  const clouds = buildClouds(baseY);
  group.add(clouds);
  const ground = buildOuterGround(tileExtentM, withOcean, baseY);
  group.add(ground);
  const foamMesh = ground.userData["foamMesh"] as THREE.Mesh | undefined;
  const waterNrm = ground.userData["waterNormal"] as THREE.Texture | undefined;
  const rail = buildRail(railLengthM);
  rail.position.set(launch[0], launch[1], launch[2]);
  group.add(rail);
  const fire = buildCampfire();
  // Camp props, each seated on the sampled terrain.
  for (const p of campLayout()) {
    let obj: THREE.Object3D;
    switch (p.kind) {
      case "hangar":
        obj = buildBuilding(4.8, 12.5, 2.6);
        break;
      case "shack":
        obj = buildBuilding(4.2, 6.8, 2.4);
        break;
      case "campfire":
        obj = fire.group;
        break;
      case "chair":
        obj = buildChair();
        break;
      case "barrel":
        obj = buildBarrel();
        break;
      case "workbench":
        obj = buildWorkbench();
        break;
      case "toolchest":
        obj = buildToolchest();
        break;
      case "crate":
        obj = buildCrate();
        break;
      case "battery":
        obj = buildBatteryBox();
        break;
      case "oilcan":
        obj = buildOilCan();
        break;
      case "windpost":
        obj = buildWindPost();
        break;
    }
    obj.position.set(launch[0] + p.x, launch[1] + groundY(p.x, p.z), launch[2] + p.z);
    obj.rotation.y = p.rotY;
    group.add(obj);
  }
  // guzez.14: the anthropometric Orville (photo face, hinged joints)
  // replaces the v1 block figure; buildFigure stays for reference.
  const orvilleFig = createBrotherFigure("orville");
  group.add(orvilleFig.group);
  const fleet = gullFleet(14, 1903);
  const gulls = fleet.map((path) => {
    const rig = buildGull(path);
    group.add(rig.group);
    return rig;
  });
  // Low flyby birds (T3.2): same gull mesh, straight tracks near the flat.
  const flybys = (withLife ? flybyFleet(4, 1907) : []).map((path: FlybyPath) => {
    const mesh = buildGullMesh();
    mesh.group.visible = false;
    group.add(mesh.group);
    return { mesh, path };
  });
  // --- life systems (T1.3/T1.6/T1.7/T2.3/T2.5) ---
  const dustTex = softDotTexture("rgba(196,178,138,0.55)", "rgba(196,178,138,0)");
  const smokeTex = softDotTexture("rgba(122,116,106,0.5)", "rgba(122,116,106,0)");
  const emberTex = softDotTexture("rgba(255,176,64,0.95)", "rgba(255,84,10,0)");
  let veg: THREE.Group | null = null;
  const streamers = new Map<number, { mesh: THREE.Mesh; write: (pts: readonly { x: number; y: number; z: number }[]) => void }>();
  const STREAMER_SEGS = 10;
  let flagpole: Flagpole | null = null;
  const smoke = makeSpritePool(14, smokeTex, 1.25);
  const embers = makeSpritePool(10, emberTex, 0.16);
  const propwash = makeSpritePool(26, dustTex, 0.9);
  const exhaust = makeSpritePool(10, smokeTex, 0.7);
  const dust = makeSpritePool(14, dustTex, 1.1);
  if (withLife) {
    veg = buildVegetation(scrubField({ tufts: 340, bushes: 56, pines: 12 }, 1903), groundY);
    group.add(veg);
    // Beach detail (A5): shells + pebbles outside the launch flat.
    group.add(buildScatter({ shells: 260, pebbles: 220 }, groundY));
    for (let i = 0; i < 9; i += 1) {
      const r = buildRibbon(STREAMER_SEGS, 0.055, 0xdccda6, 0.42);
      r.mesh.visible = false;
      group.add(r.mesh);
      streamers.set(i, r);
    }
    flagpole = buildFlagpole();
    const fx = -13.5;
    const fz = -7.5;
    flagpole.group.position.set(launch[0] + fx, launch[1] + groundY(fx, fz), launch[2] + fz);
    group.add(flagpole.group);
    // Fire effects anchor at the campfire ring (campLayout fire pos).
    const FIRE_X = -18;
    const FIRE_Z = -11;
    const fy = groundY(FIRE_X, FIRE_Z);
    smoke.group.position.set(launch[0] + FIRE_X, launch[1] + fy, launch[2] + FIRE_Z);
    embers.group.position.copy(smoke.group.position);
    group.add(smoke.group, embers.group);
  }
  group.add(dust.group);
  // Distant flock (T3.3): far gull silhouettes riding the wind east of
  // the flat — cheap billboards, pure-in-t drift.
  const flockTex = canvasTexture(96, (ctx, s) => {
    ctx.clearRect(0, 0, s, s);
    ctx.strokeStyle = "rgba(235,236,238,0.9)";
    ctx.lineWidth = 4;
    ctx.lineCap = "round";
    ctx.beginPath();
    ctx.moveTo(s * 0.1, s * 0.62);
    ctx.quadraticCurveTo(s * 0.34, s * 0.38, s * 0.5, s * 0.55);
    ctx.quadraticCurveTo(s * 0.66, s * 0.38, s * 0.9, s * 0.62);
    ctx.stroke();
  });
  const farFlock = new THREE.Group();
  const flockRand = lcg(1919);
  for (let i = 0; i < 7; i += 1) {
    const mat = new THREE.SpriteMaterial({
      map: flockTex,
      transparent: true,
      opacity: 0.5 + flockRand() * 0.35,
      depthWrite: false,
      fog: true,
    });
    const spr = new THREE.Sprite(mat);
    spr.scale.set(9 + flockRand() * 8, 5 + flockRand() * 3, 1);
    spr.userData["fx"] = 320 + flockRand() * 700;
    spr.userData["fy"] = 45 + flockRand() * 70;
    spr.userData["fz"] = -420 + flockRand() * 840;
    spr.userData["fs"] = 1.2 + flockRand() * 2.2;
    farFlock.add(spr);
  }
  if (withLife) {
    group.add(farFlock);
  }
  group.add(propwash.group, exhaust.group);
  const daniels = buildDaniels();
  const DX = -9;
  const DZ = 13;
  daniels.group.position.set(DX, 0, DZ);
  // He photographs the machine departing +x: the builder aims the Korona
  // west (-x), so flip him to face the flight path (fresh-eyes audit #8).
  daniels.group.rotation.y = Math.PI;
  group.add(daniels.group);
  // Daniels' flash latch: fires once per run release, presentation-only.
  let flashAtT: number | null = null;
  let armedRelease: number | null = null;
  // QoS particle tier witness (set by setParticleLevel, read by animate):
  let plumesHidden = false;
  return {
    group,
    animate(t, orville): void {
      // The sea hangs off ONE wall-clock t (frame rate never changes
      // the ocean): swell phase, shoreward normal scroll, and the
      // breathing surf-foam line below.
      OCEAN_TIME.value = t;
      if (waterNrm !== undefined) {
        // Wind (from the east) drives waves shoreward (-x): scroll the
        // normal field against it; the slow cross term breaks the
        // symmetry so the sparkle is not a conveyor belt.
  // (He was silently parked at the origin before this — inside the
  // rail corridor — because DX/DZ were declared but never applied.)
  daniels.group.position.set(launch[0] + DX, launch[1] + groundY(DX, DZ), launch[2] + DZ);
        const drift = Math.max(2, windMps);
        waterNrm.offset.set(
          (((-t * 0.016 * drift) % 1) + 1) % 1,
          (((t * 0.005 * drift) % 1) + 1) % 1,
        );
      }
      if (foamMesh !== undefined) {
        const foamMat = foamMesh.material as THREE.MeshBasicMaterial;
        // The surf line advances and retreats ~2.4 m with the swell
        // beat and its opacity surges on the breaking phase.
        const surfPhase = t * 0.62;
        foamMesh.position.x = tileExtentM / 2 + 26 + Math.sin(surfPhase + 1.1) * 2.4;
        foamMat.opacity = 0.36 + 0.32 * (0.5 + 0.5 * Math.sin(surfPhase));
        const foamMap = foamMat.map;
        if (foamMap !== null) {
          foamMap.offset.y = (t * 0.045) % 1;
        }
      }
      // Clouds drift east — position is a pure function of t (frame
      // rate never changes the weather).
      for (const cloud of clouds.children) {
        const baseX = cloud.userData["baseX"] as number;
        const drift = cloud.userData["driftMps"] as number;
        const span = 2800;
        cloud.position.x = ((((baseX + drift * t + 1400) % span) + span) % span) - 1400;
      }
      // Campfire flicker (deterministic in t).
      const flick = 0.82 + 0.18 * Math.sin(11 * t) * Math.sin(17.3 * t + 1.2);
      fire.light.intensity = 14 * flick;
      fire.flame.scale.setScalar(0.9 + 0.2 * flick);
      // Orville (guzez.14): the tested gait model drives the hinged
      // figure; while the machine is on the rail within reach, his
      // LEFT arm tracks the actual right wingtip — the Daniels-photo
      const pose = orvillePose(
        t,
        orville.onRail,
        orville.aircraftX,
        orville.releaseX,
        orville.releaseT,
        orville.landedX ?? null,
      );
      orvilleFig.group.position.set(
        launch[0] + pose.x,
        launch[1] + groundY(pose.x, pose.z),
        launch[2] + pose.z,
      );
      orvilleFig.group.rotation.y = pose.headingRad;
      if (pose.glassesUp) {
        // Face the machine he is watching: Ry(phi) sends +x to
        // (cos phi, -sin phi), so phi = atan2(-dz, dx) for target
        // direction (dx, dz) = (aircraftX - x, 0 - z).
        orvilleFig.group.rotation.y = Math.atan2(pose.z, orville.aircraftX - pose.x);
      }
      const chaseGap = orville.aircraftX - pose.x;
      orvilleFig.setGlasses(pose.glassesUp);
      orvilleFig.setGait(pose.gaitRad, pose.gaitRad > 0 ? (orville.landedX ? 4.6 : 3.8) : 0);
      orvilleFig.aimLeftArm(
        orville.onRail && chaseGap < 2.2
          ? [chaseGap, 1.45, 6.15 - pose.z] // wingtip: fore, wing-high, on his left
          : null,
      );
      for (const rig of gulls) {
        const p = gullPose(rig.path, t);
        const att = gullAttitude(rig.path, t);
        rig.group.position.set(launch[0] + p.x, launch[1] + p.y, launch[2] + p.z);
        // Yaw faces travel; roll banks into the orbit; pitch follows
        // the climb rate. Small angles — Euler order is immaterial.
        rig.group.rotation.set(att.rollRad, -p.headingRad, att.pitchRad);
        rig.leftWing.rotation.x = p.flapRad;
        rig.rightWing.rotation.x = -p.flapRad;
      }
      // Flybys: hide while outside their active span (no mesh churn).
      for (const fb of flybys) {
        const fp = flybyPose(fb.path, t);
        if (fp === null) {
          fb.mesh.group.visible = false;
          continue;
        }
        fb.mesh.group.visible = true;
        fb.mesh.group.position.set(launch[0] + fp.x, launch[1] + fp.y, launch[2] + fp.z);
        fb.mesh.group.rotation.set(0, -fp.headingRad, 0);
        fb.mesh.leftWing.rotation.x = fp.flapRad;
        fb.mesh.rightWing.rotation.x = -fp.flapRad;
      }
      // --- wind-made-visible + fire + machine plumes (pure in t) ---
      SWAY_TIME.value = t;
      // Sand streamers: ribbons rewrite every frame from the pure law.
      // QoS Critical hides them — fold the tier in HERE, or this
      // per-frame visible write would resurrect what setParticleLevel hid.
      for (const [i, r] of streamers) {
        const show = withLife && windMps > 3 && !plumesHidden;
        r.mesh.visible = show;
        if (!show) {
          continue;
        }
        const pts: { x: number; y: number; z: number }[] = [];
        for (let s = 0; s <= STREAMER_SEGS; s += 1) {
          const p = streamerPoint(i, s, STREAMER_SEGS, t, windMps);
          // Ride the terrain along the whole ribbon (the anchors run
          // past the launch flat, where dune relief is nonzero).
          p.y += groundY(p.x, p.z);
          pts.push(p);
        }
        r.write(pts);
      }
      // Camp flag cloth.
      if (flagpole !== null && withLife) {
        const pts: { x: number; y: number; z: number }[] = [];
        for (let s = 0; s <= 8; s += 1) {
          pts.push(flagPoint(s, 8, t, windMps));
        }
        flagpole.writeFlag(pts);
      }
      // Fire smoke + embers (strength breathes with the flame flicker).
      const strength = 0.55 + 0.45 * flick;
      for (let i = 0; i < 14; i += 1) {
        if (!withLife) {
          break;
        }
        const p = smokePuff(i, t, windMps, strength);
        smoke.set(i, p.x, p.y, p.z, p.scale * 1.15, p.opacity);
      }
      for (let i = 0; i < 10; i += 1) {
        if (!withLife) {
          break;
        }
        const e = emberAt(i, t);
        embers.set(i, e.x, e.y, e.z, 0.16, e.opacity * 0.85);
      }
      // Machine plumes: propwash sand blast behind the props and oil
      // exhaust from the crank — only while a machine state exists.
      if (orville.machine !== undefined && orville.machine.rpm01 > 0.02 && withLife) {
        const mach = orville.machine;
        for (let i = 0; i < 26; i += 1) {
          const q = propwashPuff(i, t, mach.x, Math.min(1, mach.rpm01 * 1.3));
          propwash.set(
            i,
            launch[0] + q.x,
            launch[1] + 0.18 + q.y,
            launch[2] + q.z,
            q.scale,
            q.opacity,
          );
        }
        for (let i = 0; i < 10; i += 1) {
          const x2 = exhaustPuff(i, t, mach.rpm01);
          exhaust.set(
            i,
            launch[0] + mach.x + x2.dx,
            launch[1] + 1.35 + x2.dy,
            launch[2] + x2.dz,
            x2.scale,
            x2.opacity,
          );
        }
      } else {
        propwash.hideAll();
        exhaust.hideAll();
      }
      // Landing dust burst (T3.5): outward ring from the touchdown spot.
      const dustT = orville.machine?.dustT;
      if (dustT !== undefined && withLife) {
        for (let i = 0; i < 14; i += 1) {
          const d = landingDust(i, dustT, t);
          if (d === null) {
            dust.hideAll(); // past the burst window: no frozen sprites
            continue;
          }
          dust.set(
            i,
            launch[0] + orville.aircraftX + d.dx,
            launch[1] + d.dy,
            launch[2] + d.dz,
            d.scale,
            d.opacity,
          );
        }
      } else {
        dust.hideAll();
      }
      // Distant flock drifts downwind (-x), wrapping across the sky.
      if (withLife) {
        const span = 1400;
        for (const spr of farFlock.children) {
          const fx = spr.userData["fx"] as number;
          const fs = spr.userData["fs"] as number;
          spr.position.set(
            launch[0] +
              ((((fx - windMps * 0.35 * t * fs + span / 2) % span) + span) % span) -
              span / 2,
            launch[1] + (spr.userData["fy"] as number) + 2.2 * Math.sin(t * 0.21 * fs + fs),
            launch[2] + (spr.userData["fz"] as number),
          );
        }
      }
      // Daniels' plate: ONE flash when THIS run's release is latched
      // (the historic first-photo beat); re-arms per new run.
      if (orville.releaseT !== null && armedRelease !== orville.releaseT) {
        armedRelease = orville.releaseT;
        flashAtT = t;
        // Shutter sound seam: fire EXACTLY at the latch — an edge keyed
        // on the fop window could be stepped over by one slow frame.
        window.dispatchEvent(new CustomEvent("wf-flash"));
      }
      const fsince = flashAtT === null ? Number.POSITIVE_INFINITY : t - flashAtT;
      const fop = fsince >= 0 && fsince <= 0.26 ? Math.sin((fsince / 0.26) * Math.PI) : 0;
      (daniels.flash.material as THREE.SpriteMaterial).opacity = fop;
      daniels.lamp.intensity = 30 * fop;
    },
    orvillePosition(t, orville) {
      const pose = orvillePose(
        t,
        orville.onRail,
        orville.aircraftX,
        orville.releaseX,
        orville.releaseT,
        orville.landedX ?? null,
      );
      return [
        launch[0] + pose.x,
        launch[1] + groundY(pose.x, pose.z) + 1.5,
        launch[2] + pose.z,
      ];
    },
    setParticleLevel(level): void {
      // Secondary plumes go first, everything airborne at Critical.
      const reduced = level >= 1;
      const none = level >= 2;
      plumesHidden = none;
      exhaust.group.visible = !reduced && !none;
      dust.group.visible = !reduced && !none;
      embers.group.visible = !none;
      propwash.group.visible = !none;
      smoke.group.visible = !none;
      for (const [, r] of streamers) {
        r.mesh.visible = !none;
      }
      farFlock.visible = !none;
    },
  };
}
