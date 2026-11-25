// Studio photographic HDR lightbox, procedural industrial metal machining PBR normal/roughness
// maps, CIE 1931 Planck blackbody radiance LUT, and GPU onBeforeCompile per-pixel continuum
// false-color field shaders (Temperature, von Mises Stress, Steam Pressure, Mach, Wetness/Creep)
// modeled directly on FrankenSim's Chronometer (apps/pocket-watch/src/watchScene.ts) pipeline.

import * as THREE from "three";

export type FieldShaderMode = "pbr" | "temp" | "stress" | "pressure" | "mach" | "wetcreep";

export interface FieldModeSpec {
  id: FieldShaderMode;
  modeIndex: number;
  label: string;
  shortLabel: string;
  unit: string;
  minLabel: string;
  midLabel: string;
  maxLabel: string;
  gradientCss: string;
}

export const FIELD_MODE_SPECS: Record<FieldShaderMode, FieldModeSpec> = {
  pbr: {
    id: "pbr",
    modeIndex: 0,
    label: "Studio CAD PBR",
    shortLabel: "PBR",
    unit: "",
    minLabel: "",
    midLabel: "",
    maxLabel: "",
    gradientCss: "",
  },
  temp: {
    id: "temp",
    modeIndex: 1,
    label: "Temperature T(r, z)",
    shortLabel: "Temp °C",
    unit: "°C",
    minLabel: "30 °C",
    midLabel: "290 °C",
    maxLabel: "550 °C",
    gradientCss: "linear-gradient(90deg, #070b2b 0%, #3b1782 24%, #b5245a 50%, #f57c16 76%, #fcf79a 100%)",
  },
  stress: {
    id: "stress",
    modeIndex: 2,
    label: "Von Mises Stress σ_vM",
    shortLabel: "Stress MPa",
    unit: "MPa",
    minLabel: "0 MPa",
    midLabel: "280 MPa",
    maxLabel: "560 MPa",
    gradientCss: "linear-gradient(90deg, #1b2852 0%, #1c85b8 25%, #37d484 50%, #f6c424 75%, #ef233c 100%)",
  },
  pressure: {
    id: "pressure",
    modeIndex: 3,
    label: "Steam Pressure P(z)",
    shortLabel: "Press bar",
    unit: "bar (log)",
    minLabel: "0.05 bar",
    midLabel: "2.8 bar",
    maxLabel: "165 bar",
    gradientCss: "linear-gradient(90deg, #1a1040 0%, #2a568f 28%, #20a387 56%, #95d840 82%, #fde725 100%)",
  },
  mach: {
    id: "mach",
    modeIndex: 4,
    label: "Flow Mach Number Ma",
    shortLabel: "Mach",
    unit: "Ma",
    minLabel: "0.0",
    midLabel: "0.70",
    maxLabel: "1.40",
    gradientCss: "linear-gradient(90deg, #0b1e3d 0%, #1482c4 30%, #3dd6b5 54%, #facc15 78%, #ff3b5c 100%)",
  },
  wetcreep: {
    id: "wetcreep",
    modeIndex: 5,
    label: "LP Moisture / HP Creep",
    shortLabel: "Wet/Creep",
    unit: "Risk",
    minLabel: "LP Wet 18%",
    midLabel: "Nominal",
    maxLabel: "HP Creep",
    gradientCss: "linear-gradient(90deg, #00c6ff 0%, #1d5b9e 28%, #1e293b 50%, #d97706 76%, #ff2e4c 100%)",
  },
};

export interface ContinuumShaderUniforms {
  uFieldA: { value: THREE.DataTexture };
  uFieldB: { value: THREE.DataTexture };
  uBlackbody: { value: THREE.DataTexture };
  uFieldMode: { value: number };
  uGridBounds: { value: THREE.Vector4 }; // (zMin, zSpan, rMax, isoContours)
  uTurbineInvMat: { value: THREE.Matrix4 };
}

/**
 * Industrial studio light cyclorama for PMREM image-based lighting (inspired by watchScene.ts):
 * crisp overhead key softbox, tall vertical specular strip lights on both sides to draw clean
 * longitudinal reflection ribbons along the CNC-turned rotor shaft and machined flanges,
 * plus warm concrete bounce so undersides retain depth without turning muddy.
 */
export function buildIndustrialStudioEnvironment(): THREE.Scene {
  const env = new THREE.Scene();
  env.background = new THREE.Color(0x05070b);

  const dome = new THREE.Mesh(
    new THREE.SphereGeometry(65, 48, 24),
    new THREE.MeshBasicMaterial({ color: new THREE.Color(0x0e131b), side: THREE.BackSide }),
  );
  env.add(dome);

  // Subtle industrial hall floor bounce
  const floor = new THREE.Mesh(
    new THREE.CircleGeometry(45, 48),
    new THREE.MeshBasicMaterial({ color: new THREE.Color().setRGB(0.08, 0.085, 0.095) }),
  );
  floor.rotation.x = -Math.PI / 2;
  floor.position.set(0, -24, 0);
  env.add(floor);

  const panel = (
    w: number,
    h: number,
    pos: [number, number, number],
    intensity: number,
    tint: [number, number, number],
  ): void => {
    const m = new THREE.Mesh(
      new THREE.PlaneGeometry(w, h),
      new THREE.MeshBasicMaterial({
        color: new THREE.Color().setRGB(tint[0] * intensity, tint[1] * intensity, tint[2] * intensity),
        side: THREE.DoubleSide,
      }),
    );
    m.position.set(...pos);
    m.lookAt(0, 0, 0);
    env.add(m);
  };

  // Large overhead engineering lightbox softbox (balanced for crisp PBR highlights without glare blowout!)
  panel(44, 18, [2, 38, 12], 2.1, [0.98, 0.99, 1.0]);
  // Front-left warm key strip (catches turned rotor discs & LP twisted blades cleanly)
  panel(10, 42, [-36, 14, 28], 2.35, [1.0, 0.96, 0.90]);
  // Right cool steel kicker strip (separates slate-blue exhaust ribs & machined flanges)
  panel(10, 42, [38, 10, 24], 1.85, [0.84, 0.92, 1.0]);
  // Rear-top rim strip (carves crisp silhouette on twin CV rack & exhaust hood)
  panel(46, 6, [0, 32, -36], 1.75, [0.90, 0.95, 1.0]);
  // Low warm bounce strip illuminating lower split flange & matte black skid
  panel(28, 6, [0, -14, 34], 0.85, [0.96, 0.92, 0.86]);

  return env;
}

/**
 * Exact CIE 1931 2° Planck blackbody radiance LUT (300 K .. 2400 K, 256 texels)
 * identical in physics to apps/pocket-watch/src/watchScene.ts.
 */
export function buildBlackbodyLut(): THREE.DataTexture {
  const g = (x: number, mu: number, s1: number, s2: number): number => {
    const t = (x - mu) / (x < mu ? s1 : s2);
    return Math.exp(-0.5 * t * t);
  };
  const cmf = (l: number): [number, number, number] => [
    1.056 * g(l, 599.8, 37.9, 31.0) + 0.362 * g(l, 442.0, 16.0, 26.7) - 0.065 * g(l, 501.1, 20.4, 26.2),
    0.821 * g(l, 568.8, 46.9, 40.5) + 0.286 * g(l, 530.9, 16.3, 31.1),
    1.217 * g(l, 437.0, 11.8, 36.0) + 0.681 * g(l, 459.0, 26.0, 13.8),
  ];
  const planck = (lNm: number, tK: number): number => {
    const l = lNm * 1e-9;
    return 1.191042972e-16 / (l ** 5 * (Math.exp(1.438776877e-2 / (l * tK)) - 1.0));
  };
  const xyzAt = (tK: number): [number, number, number] => {
    let X = 0, Y = 0, Z = 0;
    for (let l = 380; l <= 780; l += 5) {
      const b = planck(l, tK);
      const [xb, yb, zb] = cmf(l);
      X += b * xb; Y += b * yb; Z += b * zb;
    }
    return [X, Y, Z];
  };
  const yRef = xyzAt(1250.0)[1];
  const n = 256;
  const data = new Float32Array(n * 4);
  for (let i = 0; i < n; i++) {
    const tK = 300.0 + (2100.0 * i) / (n - 1);
    const [X, Y, Z] = xyzAt(tK);
    const inv = 1.0 / Math.max(Y, 1e-300);
    const xn = X * inv, zn = Z * inv;
    let r = Math.max(0, 3.2406 * xn - 1.5372 - 0.4986 * zn);
    let gg = Math.max(0, -0.9689 * xn + 1.8758 + 0.0415 * zn);
    let b = Math.max(0, 0.0557 * xn - 0.204 + 1.057 * zn);
    const lum = Math.max(0.2126 * r + 0.7152 * gg + 0.0722 * b, 1e-9);
    data[i * 4] = r / lum;
    data[i * 4 + 1] = gg / lum;
    data[i * 4 + 2] = b / lum;
    data[i * 4 + 3] = Y / yRef;
  }
  const tex = new THREE.DataTexture(data, n, 1, THREE.RGBAFormat, THREE.FloatType);
  tex.minFilter = THREE.NearestFilter;
  tex.magFilter = THREE.NearestFilter;
  tex.needsUpdate = true;
  return tex;
}

export interface IndustrialTextures {
  turnedSteelNormal: THREE.CanvasTexture;
  bladeBrushNormal: THREE.CanvasTexture;
  castSteelNormal: THREE.CanvasTexture;
  castSteelRoughness: THREE.CanvasTexture;
  milledSectionNormal: THREE.CanvasTexture;
  siemensEnamelNormal: THREE.CanvasTexture;
  skidPlateNormal: THREE.CanvasTexture;
  softGlowTex: THREE.CanvasTexture;
}

function hash2(x: number, y: number): number {
  const s = Math.sin(x * 127.1 + y * 311.7) * 43758.5453;
  return s - Math.floor(s);
}

function makeRepeatTex(canvas: HTMLCanvasElement, rx: number, ry: number): THREE.CanvasTexture {
  const tex = new THREE.CanvasTexture(canvas);
  tex.wrapS = THREE.RepeatWrapping;
  tex.wrapT = THREE.RepeatWrapping;
  tex.repeat.set(rx, ry);
  tex.needsUpdate = true;
  return tex;
}

/**
 * Procedural industrial machining normal & roughness maps (analogous to Chronometer's Geneva
 * stripes & Perlage) so every CrMoV forging, cast steel casing, face-milled cutaway cap,
 * brushed 12Cr/Ti airfoil, and painted Siemens exhaust hood reflects realistic micro-structure!
 */
export function buildIndustrialMachiningTextures(): IndustrialTextures {
  const N = 512;

  // 1. CNC Lathe-Turned Forged Rotor Steel (circumferential feed micro-grooves)
  const cTurn = document.createElement("canvas");
  cTurn.width = N; cTurn.height = N;
  const tCtx = cTurn.getContext("2d")!;
  const tImg = tCtx.createImageData(N, N);
  for (let y = 0; y < N; y++) {
    for (let x = 0; x < N; x++) {
      // u = axial (x), v = circumferential (y) -> grooves run along y (vary sharply with x)
      const fine = Math.sin(x * 1.15) * 0.34 + Math.sin(x * 3.7 + 0.8) * 0.14;
      const chatter = (hash2(x >> 1, y >> 3) - 0.5) * 0.14;
      const nx = fine + chatter;
      const ny = (hash2(x, y >> 2) - 0.5) * 0.06;
      const idx = (y * N + x) * 4;
      tImg.data[idx] = Math.round(Math.max(0, Math.min(255, (nx * 0.5 + 0.5) * 255)));
      tImg.data[idx + 1] = Math.round(Math.max(0, Math.min(255, (ny * 0.5 + 0.5) * 255)));
      tImg.data[idx + 2] = 246;
      tImg.data[idx + 3] = 255;
    }
  }
  tCtx.putImageData(tImg, 0, 0);
  const turnedSteelNormal = makeRepeatTex(cTurn, 3.0, 2.0);

  // 2. Spanwise Brushed 12Cr / Titanium Airfoil Grain
  const cBlade = document.createElement("canvas");
  cBlade.width = 256; cBlade.height = 256;
  const bCtx = cBlade.getContext("2d")!;
  const bImg = bCtx.createImageData(256, 256);
  for (let y = 0; y < 256; y++) {
    for (let x = 0; x < 256; x++) {
      const streak = Math.sin(x * 1.9 + hash2(x, y >> 4) * 1.2) * 0.28;
      const nx = streak;
      const ny = Math.sin(y * 0.45) * 0.05;
      const idx = (y * 256 + x) * 4;
      bImg.data[idx] = Math.round((nx * 0.5 + 0.5) * 255);
      bImg.data[idx + 1] = Math.round((ny * 0.5 + 0.5) * 255);
      bImg.data[idx + 2] = 248;
      bImg.data[idx + 3] = 255;
    }
  }
  bCtx.putImageData(bImg, 0, 0);
  const bladeBrushNormal = makeRepeatTex(cBlade, 2.0, 2.0);

  // 3. Sand-Cast Nodular CrMoV Steel Normal + Roughness Variation (for outer HP/IP casing & carrier)
  const cCast = document.createElement("canvas");
  cCast.width = 256; cCast.height = 256;
  const castCtx = cCast.getContext("2d")!;
  const castImg = castCtx.createImageData(256, 256);
  const cCastR = document.createElement("canvas");
  cCastR.width = 256; cCastR.height = 256;
  const castRCtx = cCastR.getContext("2d")!;
  const castRImg = castRCtx.createImageData(256, 256);
  const field = new Float32Array(256 * 256);
  for (let y = 0; y < 256; y++) {
    for (let x = 0; x < 256; x++) {
      const h1 = hash2(x >> 2, y >> 2);
      const h2 = hash2(x >> 1, y >> 1);
      const h3 = hash2(x, y);
      field[y * 256 + x] = 0.5 * h1 + 0.32 * h2 + 0.18 * h3;
    }
  }
  for (let y = 0; y < 256; y++) {
    for (let x = 0; x < 256; x++) {
      const xm = (x + 255) & 255, xp = (x + 1) & 255;
      const ym = (y + 255) & 255, yp = (y + 1) & 255;
      const gx = (field[y * 256 + xp]! - field[y * 256 + xm]!) * 1.4;
      const gy = (field[yp * 256 + x]! - field[ym * 256 + x]!) * 1.4;
      const idx = (y * 256 + x) * 4;
      castImg.data[idx] = Math.round(Math.max(0, Math.min(255, (gx * 0.5 + 0.5) * 255)));
      castImg.data[idx + 1] = Math.round(Math.max(0, Math.min(255, (gy * 0.5 + 0.5) * 255)));
      castImg.data[idx + 2] = 238;
      castImg.data[idx + 3] = 255;
      const rv = Math.round(135 + field[y * 256 + x]! * 68);
      castRImg.data[idx] = rv;
      castRImg.data[idx + 1] = rv;
      castRImg.data[idx + 2] = rv;
      castRImg.data[idx + 3] = 255;
    }
  }
  castCtx.putImageData(castImg, 0, 0);
  castRCtx.putImageData(castRImg, 0, 0);
  const castSteelNormal = makeRepeatTex(cCast, 4.0, 4.0);
  const castSteelRoughness = makeRepeatTex(cCastR, 4.0, 4.0);

  // 4. Precision Face-Milled Cutaway Section Cross-Hatch (signature CAD cutaway look!)
  const cMill = document.createElement("canvas");
  cMill.width = 256; cMill.height = 256;
  const mCtx = cMill.getContext("2d")!;
  const mImg = mCtx.createImageData(256, 256);
  for (let y = 0; y < 256; y++) {
    for (let x = 0; x < 256; x++) {
      // Crisp 45-deg fly-cutter milling arcs + light cross-feed marks
      const diag = ((x + y) / 256) * 24.0 * Math.PI * 2;
      const micro = Math.sin(diag) * 0.36 + Math.sin(((x - y) / 256) * 12.0 * Math.PI * 2) * 0.08;
      const idx = (y * 256 + x) * 4;
      mImg.data[idx] = Math.round((micro * 0.5 + 0.5) * 255);
      mImg.data[idx + 1] = Math.round((-micro * 0.5 + 0.5) * 255);
      mImg.data[idx + 2] = 244;
      mImg.data[idx + 3] = 255;
    }
  }
  mCtx.putImageData(mImg, 0, 0);
  const milledSectionNormal = makeRepeatTex(cMill, 3.5, 3.5);

  // 5. Siemens Industrial Satin Enamel Paint (exhaust hood micro-grain)
  const cEnam = document.createElement("canvas");
  cEnam.width = 256; cEnam.height = 256;
  const eCtx = cEnam.getContext("2d")!;
  const eImg = eCtx.createImageData(256, 256);
  for (let y = 0; y < 256; y++) {
    for (let x = 0; x < 256; x++) {
      const gx = (Math.sin(x * 0.35 + y * 0.21) * Math.cos(y * 0.41 - x * 0.17)) * 0.16;
      const gy = (Math.cos(x * 0.29 - y * 0.33) * Math.sin(x * 0.19 + y * 0.27)) * 0.16;
      const idx = (y * 256 + x) * 4;
      eImg.data[idx] = Math.round((gx * 0.5 + 0.5) * 255);
      eImg.data[idx + 1] = Math.round((gy * 0.5 + 0.5) * 255);
      eImg.data[idx + 2] = 250;
      eImg.data[idx + 3] = 255;
    }
  }
  eCtx.putImageData(eImg, 0, 0);
  const siemensEnamelNormal = makeRepeatTex(cEnam, 5.0, 5.0);

  // 6. Hot-Rolled Welded Matte Black Skid Plate Grain
  const cSkid = document.createElement("canvas");
  cSkid.width = 256; cSkid.height = 256;
  const kCtx = cSkid.getContext("2d")!;
  const kImg = kCtx.createImageData(256, 256);
  for (let y = 0; y < 256; y++) {
    for (let x = 0; x < 256; x++) {
      const n = (hash2(x >> 1, y >> 1) - 0.5) * 0.32;
      const idx = (y * 256 + x) * 4;
      kImg.data[idx] = Math.round((n * 0.5 + 0.5) * 255);
      kImg.data[idx + 1] = Math.round((n * 0.5 + 0.5) * 255);
      kImg.data[idx + 2] = 244;
      kImg.data[idx + 3] = 255;
    }
  }
  kCtx.putImageData(kImg, 0, 0);
  const skidPlateNormal = makeRepeatTex(cSkid, 6.0, 4.0);

  // 7. Crisp spark / mist radial glow sprite (for genuine rub sparks & subtle LP condensation)
  const cGlow = document.createElement("canvas");
  cGlow.width = 64; cGlow.height = 64;
  const gCtx = cGlow.getContext("2d")!;
  const grad = gCtx.createRadialGradient(32, 32, 0, 32, 32, 32);
  grad.addColorStop(0.0, "rgba(255,255,255,1.0)");
  grad.addColorStop(0.28, "rgba(255,235,200,0.65)");
  grad.addColorStop(1.0, "rgba(255,180,100,0.0)");
  gCtx.fillStyle = grad;
  gCtx.fillRect(0, 0, 64, 64);
  const softGlowTex = new THREE.CanvasTexture(cGlow);

  return {
    turnedSteelNormal,
    bladeBrushNormal,
    castSteelNormal,
    castSteelRoughness,
    milledSectionNormal,
    siemensEnamelNormal,
    skidPlateNormal,
    softGlowTex,
  };
}

/**
 * Inject GPU `onBeforeCompile` per-pixel continuum field colormaps & Planck blackbody radiance
 * into any Three.js MeshStandardMaterial / MeshPhysicalMaterial on the turbine assembly.
 *
 * Transforms every fragment's world position into exact SST-600 axisymmetric `(z_m, r_m)` via
 * `uTurbineInvMat`, bilinearly samples the `64 × 24` WASM PDE textures (`uFieldA`, `uFieldB`),
 * draws calibrated scientific iso-contours, and blends seamlessly with HDR PBR studio lighting!
 */
export function installContinuumFieldShader(
  mat: THREE.MeshStandardMaterial | THREE.MeshPhysicalMaterial,
  uniforms: ContinuumShaderUniforms,
  cacheKeySuffix: string,
): void {
  mat.customProgramCacheKey = () => `fs-sst600-continuum-${cacheKeySuffix}`;
  mat.onBeforeCompile = (shader) => {
    shader.uniforms.uFieldA = uniforms.uFieldA;
    shader.uniforms.uFieldB = uniforms.uFieldB;
    shader.uniforms.uBlackbody = uniforms.uBlackbody;
    shader.uniforms.uFieldMode = uniforms.uFieldMode;
    shader.uniforms.uGridBounds = uniforms.uGridBounds;
    shader.uniforms.uTurbineInvMat = uniforms.uTurbineInvMat;

    shader.vertexShader = shader.vertexShader
      .replace(
        "#include <common>",
        /* glsl */ `
#include <common>
uniform mat4 uTurbineInvMat;
varying vec2 vFsTurbZR;
`,
      )
      .replace(
        "#include <worldpos_vertex>",
        /* glsl */ `
#include <worldpos_vertex>
vec4 fsWorldP = vec4(transformed, 1.0);
#ifdef USE_INSTANCING
  fsWorldP = instanceMatrix * fsWorldP;
#endif
fsWorldP = modelMatrix * fsWorldP;
vec3 fsLocP = (uTurbineInvMat * fsWorldP).xyz;
vFsTurbZR = vec2(fsLocP.z, length(fsLocP.xy));
`,
      );

    shader.fragmentShader = shader.fragmentShader
      .replace(
        "#include <common>",
        /* glsl */ `
#include <common>
uniform sampler2D uFieldA;
uniform sampler2D uFieldB;
uniform sampler2D uBlackbody;
uniform float uFieldMode;
uniform vec4 uGridBounds;
varying vec2 vFsTurbZR;

vec4 fsBilerpTex(sampler2D tex, vec2 uv, vec2 dims) {
  vec2 p = clamp(uv, 0.0, 1.0) * (dims - 1.0);
  ivec2 i0 = ivec2(floor(p));
  ivec2 i1 = min(i0 + ivec2(1), ivec2(dims) - ivec2(1));
  vec2 f = fract(p);
  vec4 a00 = texelFetch(tex, ivec2(i0.x, i0.y), 0);
  vec4 a10 = texelFetch(tex, ivec2(i1.x, i0.y), 0);
  vec4 a01 = texelFetch(tex, ivec2(i0.x, i1.y), 0);
  vec4 a11 = texelFetch(tex, ivec2(i1.x, i1.y), 0);
  return mix(mix(a00, a10, f.x), mix(a01, a11, f.x), f.y);
}

vec3 fsInferno(float u) {
  u = clamp(u, 0.0, 1.0);
  vec3 c0 = vec3(0.02, 0.03, 0.14);
  vec3 c1 = vec3(0.24, 0.08, 0.52);
  vec3 c2 = vec3(0.74, 0.16, 0.36);
  vec3 c3 = vec3(0.97, 0.50, 0.08);
  vec3 c4 = vec3(0.99, 0.96, 0.62);
  if (u < 0.25) return mix(c0, c1, u * 4.0);
  if (u < 0.50) return mix(c1, c2, (u - 0.25) * 4.0);
  if (u < 0.75) return mix(c2, c3, (u - 0.50) * 4.0);
  return mix(c3, c4, (u - 0.75) * 4.0);
}

vec3 fsStressMap(float u) {
  u = clamp(u, 0.0, 1.0);
  vec3 c0 = vec3(0.08, 0.15, 0.34);
  vec3 c1 = vec3(0.10, 0.54, 0.74);
  vec3 c2 = vec3(0.22, 0.84, 0.52);
  vec3 c3 = vec3(0.96, 0.76, 0.14);
  vec3 c4 = vec3(0.96, 0.14, 0.24);
  if (u < 0.25) return mix(c0, c1, u * 4.0);
  if (u < 0.50) return mix(c1, c2, (u - 0.25) * 4.0);
  if (u < 0.75) return mix(c2, c3, (u - 0.50) * 4.0);
  return mix(c3, c4, (u - 0.75) * 4.0);
}

vec3 fsViridis(float u) {
  u = clamp(u, 0.0, 1.0);
  vec3 c0 = vec3(0.11, 0.06, 0.26);
  vec3 c1 = vec3(0.16, 0.35, 0.56);
  vec3 c2 = vec3(0.12, 0.64, 0.53);
  vec3 c3 = vec3(0.58, 0.84, 0.25);
  vec3 c4 = vec3(0.99, 0.90, 0.16);
  if (u < 0.25) return mix(c0, c1, u * 4.0);
  if (u < 0.50) return mix(c1, c2, (u - 0.25) * 4.0);
  if (u < 0.75) return mix(c2, c3, (u - 0.50) * 4.0);
  return mix(c3, c4, (u - 0.75) * 4.0);
}

vec3 fsMachMap(float u) {
  u = clamp(u, 0.0, 1.0);
  vec3 c0 = vec3(0.05, 0.12, 0.26);
  vec3 c1 = vec3(0.08, 0.52, 0.78);
  vec3 c2 = vec3(0.24, 0.85, 0.72);
  vec3 c3 = vec3(0.98, 0.80, 0.10);
  vec3 c4 = vec3(1.00, 0.20, 0.34);
  if (u < 0.25) return mix(c0, c1, u * 4.0);
  if (u < 0.50) return mix(c1, c2, (u - 0.25) * 4.0);
  if (u < 0.75) return mix(c2, c3, (u - 0.50) * 4.0);
  return mix(c3, c4, (u - 0.75) * 4.0);
}

vec3 fsBlackbody(float tK) {
  float x = clamp((tK - 300.0) / 2100.0, 0.0, 1.0) * 255.0;
  int i0 = int(floor(x));
  int i1 = min(i0 + 1, 255);
  vec4 c = mix(texelFetch(uBlackbody, ivec2(i0, 0), 0), texelFetch(uBlackbody, ivec2(i1, 0), 0), fract(x));
  float L = c.a / (1.0 + c.a * 0.4);
  return c.rgb * L;
}
`,
      )
      .replace(
        "#include <color_fragment>",
        /* glsl */ `
#include <color_fragment>
vec2 fsUv = vec2((vFsTurbZR.x - uGridBounds.x) / uGridBounds.y, vFsTurbZR.y / uGridBounds.z);
vec4 fsA = fsBilerpTex(uFieldA, fsUv, vec2(64.0, 24.0));
vec4 fsB = fsBilerpTex(uFieldB, fsUv, vec2(64.0, 24.0));
float fsInDomain = step(-0.02, fsUv.x) * step(fsUv.x, 1.02) * step(fsUv.y, 1.05);
vec3 fsFieldCol = vec3(0.0);
float fsContour = 0.0;
if (uFieldMode > 0.5 && fsInDomain > 0.5) {
  float scalarNorm = 0.0;
  if (uFieldMode < 1.5) {
    // Mode 1: Temperature T(r, z) [30..550 C]
    scalarNorm = clamp((fsA.x - 303.15) / 520.0, 0.0, 1.0);
    fsFieldCol = fsInferno(scalarNorm);
  } else if (uFieldMode < 2.5) {
    // Mode 2: Von Mises Stress [0..560 MPa]
    scalarNorm = clamp(fsA.y / 560.0, 0.0, 1.0);
    fsFieldCol = fsStressMap(scalarNorm);
  } else if (uFieldMode < 3.5) {
    // Mode 3: Steam Pressure P(z) log10 [0.05..165 bar]
    float lp = log(max(fsA.z, 0.04));
    scalarNorm = clamp((lp - log(0.05)) / (log(165.0) - log(0.05)), 0.0, 1.0);
    fsFieldCol = fsViridis(scalarNorm);
  } else if (uFieldMode < 4.5) {
    // Mode 4: Mach Number Ma [0..1.4]
    scalarNorm = clamp(fsA.w / 1.40, 0.0, 1.0);
    fsFieldCol = fsMachMap(scalarNorm);
  } else {
    // Mode 5: Dual LP Wetness (cyan) / HP Norton Creep (amber-crimson)
    float wet = clamp(fsB.z / 0.18, 0.0, 1.0);
    float creep = clamp((fsB.w + 15.0) / 8.0, 0.0, 1.0);
    vec3 baseSlate = vec3(0.14, 0.18, 0.24);
    vec3 wetCol = mix(baseSlate, vec3(0.02, 0.78, 1.0), wet);
    fsFieldCol = mix(wetCol, vec3(1.0, 0.24, 0.18), creep);
    scalarNorm = max(wet, creep);
  }
  // Subtle crisp engineering iso-contour bands (10 levels)
  float band = abs(fract(scalarNorm * 10.0) - 0.5);
  fsContour = (1.0 - smoothstep(0.0, 0.045, band)) * uGridBounds.w;
  vec3 mapped = fsFieldCol * (1.0 - 0.28 * fsContour);
  diffuseColor.rgb = mix(diffuseColor.rgb, mapped, 0.88);
}
`,
      )
      .replace(
        "#include <emissivemap_fragment>",
        /* glsl */ `
#include <emissivemap_fragment>
if (uFieldMode > 0.5 && fsInDomain > 0.5) {
  totalEmissiveRadiance += fsFieldCol * (0.28 - 0.10 * fsContour);
} else if (fsInDomain > 0.5 && fsA.x > 740.0) {
  // Real physical Planck blackbody incandescence when CrMoV steel approaches dull red heat!
  totalEmissiveRadiance += fsBlackbody(fsA.x) * 1.4;
}
`,
      );
  };
}
