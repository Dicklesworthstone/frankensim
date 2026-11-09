import * as THREE from "three";
import { OrbitControls } from "three/examples/jsm/controls/OrbitControls.js";
import { RoomEnvironment } from "three/examples/jsm/environments/RoomEnvironment.js";
import { EffectComposer } from "three/examples/jsm/postprocessing/EffectComposer.js";
import { OutputPass } from "three/examples/jsm/postprocessing/OutputPass.js";
import { RenderPass } from "three/examples/jsm/postprocessing/RenderPass.js";
import { SMAAPass } from "three/examples/jsm/postprocessing/SMAAPass.js";
import { UnrealBloomPass } from "three/examples/jsm/postprocessing/UnrealBloomPass.js";
import type { WatchSnapshot } from "./snapshotLayout";

export type ActiveToolMode = "inspect" | "torch" | "ln2" | "oil" | "magnet" | "air" | "mic";
export type VectorFieldLayer = "heat" | "convection" | "magnetic" | "acoustic" | "stress";
export type CameraPreset =
  | "overview"
  | "tourbillon"
  | "barrel"
  | "side"
  | "escapement"
  | "profile";

export interface SceneCallbacks {
  onTargetMoved: (tool: ActiveToolMode, posMm: [number, number, number]) => void;
}

/** Render quality tier: decides MSAA, shadow resolution, transmission, bloom and the DPR ceiling. */
export type QualityTier = "high" | "medium" | "low";

/** CSS-pixel insets of UI chrome that floats over the full-bleed canvas. */
export interface ViewInsets {
  left: number;
  right: number;
  top: number;
  bottom: number;
}

export interface SceneOptions {
  tier?: QualityTier;
  reducedMotion?: boolean;
}

export interface ToolRenderParams {
  torchPower: number;
  torchPosMm: [number, number, number];
  ln2Rate: number;
  ln2PosMm: [number, number, number];
  oilRate: number;
  oilPosMm: [number, number, number];
  magnetOn: boolean;
  magnetPosMm: [number, number, number];
  airPower?: number;
  airPosMm?: [number, number, number];
  ambientPressureAtm?: number;
  gravityGVec?: [number, number, number];
  micPosMm: [number, number, number];
  caseShellData?: Float32Array;
  spatialGridData?: Float32Array;
  mainspringCoilsMm?: Float32Array;
  clickLeverAngleRad?: number;
  crownWheelAngleRad?: number;
}

const CASE_SHELL_N_THETA = 64;
const CASE_SHELL_N_Z = 8;
const CASE_N_THETA = CASE_SHELL_N_THETA;
const CASE_N_Z = CASE_SHELL_N_Z;
const CASE_SHELL_NODES = CASE_SHELL_N_THETA * CASE_SHELL_N_Z;
const CASE_RENDER_N_THETA = 128;
const CASE_RENDER_N_Z = 22;
const CASE_RENDER_NODES = CASE_RENDER_N_THETA * CASE_RENDER_N_Z;

// The Rust spatial PDE grid size can change between engine builds: the scene infers N from the
// snapshot length (7 floats per cell) and takes the half span from `grid_half_span_mm()` when exposed.
let GRID_N = 48;
let GRID_HALF_SPAN_MM = 28.8;
let GRID_DX_MM = (2.0 * GRID_HALF_SPAN_MM) / GRID_N;
/** Vertices per side of the GPU-displaced molten/solid gold heightfield meshes. */
const HEIGHTFIELD_RES = 160;

/**
 * Blackbody emission lookup (300 K .. 3000 K, 256 texels): linear-sRGB chromaticity of Planck's law
 * integrated against the CIE 1931 2° observer (Wyman-Sloan-Shirley multi-lobe fit), with the
 * photometric luminance Y(T) relative to Y(1337 K) (the melting point of gold) in the alpha channel.
 * Glow therefore scales physically: about 1e-4 of the melt-point brightness at 900 K.
 */
const BLACKBODY_T_MIN = 300.0;
const BLACKBODY_T_MAX = 3000.0;
const BLACKBODY_TEXELS = 256;

function buildBlackbodyLut(): THREE.DataTexture {
  const g = (x: number, mu: number, s1: number, s2: number): number => {
    const t = (x - mu) / (x < mu ? s1 : s2);
    return Math.exp(-0.5 * t * t);
  };
  const cmf = (lambdaNm: number): [number, number, number] => {
    const x =
      1.056 * g(lambdaNm, 599.8, 37.9, 31.0) +
      0.362 * g(lambdaNm, 442.0, 16.0, 26.7) -
      0.065 * g(lambdaNm, 501.1, 20.4, 26.2);
    const y = 0.821 * g(lambdaNm, 568.8, 46.9, 40.5) + 0.286 * g(lambdaNm, 530.9, 16.3, 31.1);
    const z = 1.217 * g(lambdaNm, 437.0, 11.8, 36.0) + 0.681 * g(lambdaNm, 459.0, 26.0, 13.8);
    return [x, y, z];
  };
  const planck = (lambdaNm: number, tK: number): number => {
    const l = lambdaNm * 1e-9;
    const c1 = 1.191042972e-16; // 2 h c^2 [W m^2 / sr]
    const c2 = 1.438776877e-2; // h c / k [m K]
    return c1 / (l ** 5 * (Math.exp(c2 / (l * tK)) - 1.0));
  };
  const xyzAt = (tK: number): [number, number, number] => {
    let X = 0;
    let Y = 0;
    let Z = 0;
    for (let l = 380; l <= 780; l += 4) {
      const b = planck(l, tK);
      const [xb, yb, zb] = cmf(l);
      X += b * xb;
      Y += b * yb;
      Z += b * zb;
    }
    return [X, Y, Z];
  };
  const yRef = xyzAt(1337.33)[1];
  const data = new Float32Array(BLACKBODY_TEXELS * 4);
  for (let i = 0; i < BLACKBODY_TEXELS; i++) {
    const tK = BLACKBODY_T_MIN + ((BLACKBODY_T_MAX - BLACKBODY_T_MIN) * i) / (BLACKBODY_TEXELS - 1);
    const [X, Y, Z] = xyzAt(tK);
    // XYZ -> linear sRGB (D65), chromaticity normalized to unit luminance.
    const inv = 1.0 / Math.max(Y, 1e-300);
    const xn = X * inv;
    const zn = Z * inv;
    let r = 3.2406 * xn - 1.5372 - 0.4986 * zn;
    let gg = -0.9689 * xn + 1.8758 + 0.0415 * zn;
    let b = 0.0557 * xn - 0.204 + 1.057 * zn;
    r = Math.max(0, r);
    gg = Math.max(0, gg);
    b = Math.max(0, b);
    // Normalize to unit Rec.709 luminance so rgb * (Y/Yref) carries the physical brightness.
    const lum = Math.max(0.2126 * r + 0.7152 * gg + 0.0722 * b, 1e-9);
    data[i * 4] = r / lum;
    data[i * 4 + 1] = gg / lum;
    data[i * 4 + 2] = b / lum;
    data[i * 4 + 3] = Y / yRef;
  }
  const tex = new THREE.DataTexture(data, BLACKBODY_TEXELS, 1, THREE.RGBAFormat, THREE.FloatType);
  tex.minFilter = THREE.NearestFilter;
  tex.magFilter = THREE.NearestFilter;
  tex.needsUpdate = true;
  return tex;
}

/** GLSL: blackbody radiance lookup (manual lerp between texels so FloatType never needs filtering). */
const GLSL_BLACKBODY = /* glsl */ `
uniform sampler2D uBlackbody;
uniform float uGlowScale;
uniform float uThermalOverlay;
vec3 blackbodyRadiance(float tK) {
  float x = clamp((tK - ${BLACKBODY_T_MIN.toFixed(1)}) / ${(BLACKBODY_T_MAX - BLACKBODY_T_MIN).toFixed(1)}, 0.0, 1.0) * ${(BLACKBODY_TEXELS - 1).toFixed(1)};
  int i0 = int(floor(x));
  int i1 = min(i0 + 1, ${BLACKBODY_TEXELS - 1});
  vec4 a = texelFetch(uBlackbody, ivec2(i0, 0), 0);
  vec4 b = texelFetch(uBlackbody, ivec2(i1, 0), 0);
  vec4 c = mix(a, b, fract(x));
  return c.rgb * c.a * uGlowScale;
}
vec3 thermalFalseColor(float tK) {
  // Perceptually ordered "inferno"-style ramp from 77 K (LN2) to 1500 K.
  float u = clamp((tK - 77.0) / 1423.0, 0.0, 1.0);
  vec3 c0 = vec3(0.02, 0.02, 0.10);
  vec3 c1 = vec3(0.18, 0.12, 0.55);
  vec3 c2 = vec3(0.72, 0.18, 0.42);
  vec3 c3 = vec3(0.98, 0.52, 0.10);
  vec3 c4 = vec3(0.99, 0.96, 0.62);
  if (u < 0.25) return mix(c0, c1, u / 0.25);
  if (u < 0.5) return mix(c1, c2, (u - 0.25) / 0.25);
  if (u < 0.75) return mix(c2, c3, (u - 0.5) / 0.25);
  return mix(c3, c4, (u - 0.75) / 0.25);
}
`;

/** GLSL: uniform cubic B-spline weights and their derivatives. */
const GLSL_CUBIC = /* glsl */ `
void bsplineWeights(float t, out vec4 w, out vec4 dw) {
  float t2 = t * t;
  float t3 = t2 * t;
  float omt = 1.0 - t;
  w = vec4(omt * omt * omt, 3.0 * t3 - 6.0 * t2 + 4.0, -3.0 * t3 + 3.0 * t2 + 3.0 * t + 1.0, t3) / 6.0;
  dw = vec4(-0.5 * omt * omt, 1.5 * t2 - 2.0 * t, -1.5 * t2 + t + 0.5, 0.5 * t2);
}
void mitchellWeights(float t, out vec4 w) {
  float t2 = t * t;
  float t3 = t2 * t;
  w = vec4(-7.0 * t3 + 15.0 * t2 - 9.0 * t + 1.0,
           21.0 * t3 - 36.0 * t2 + 16.0,
           -21.0 * t3 + 27.0 * t2 + 9.0 * t + 1.0,
           7.0 * t3 - 6.0 * t2) / 18.0;
}
`;

function bSplineCubic(p0: number, p1: number, p2: number, p3: number, t: number): number {
  const oneMinusT = 1.0 - t;
  const t2 = t * t;
  const t3 = t2 * t;
  const b0 = (oneMinusT * oneMinusT * oneMinusT) / 6.0;
  const b1 = (3.0 * t3 - 6.0 * t2 + 4.0) / 6.0;
  const b2 = (-3.0 * t3 + 3.0 * t2 + 3.0 * t + 1.0) / 6.0;
  const b3 = t3 / 6.0;
  return b0 * p0 + b1 * p1 + b2 * p2 + b3 * p3;
}

/**
 * Mitchell-Netravali (B = 1/3, C = 1/3) cubic filter clamped to the local [min, max] convex hull.
 * Combines 1/3 C^2 B-Spline + 2/3 Catmull-Rom so free-surface heightfields are completely free of
 * grid-line creases while preserving full droplet peak height and crisp Young-Laplace meniscus edges.
 */
function mitchellNetravaliClamped(p0: number, p1: number, p2: number, p3: number, t: number): number {
  const t2 = t * t;
  const t3 = t2 * t;
  const w0 = (-7.0 * t3 + 15.0 * t2 - 9.0 * t + 1.0) / 18.0;
  const w1 = (21.0 * t3 - 36.0 * t2 + 16.0) / 18.0;
  const w2 = (-21.0 * t3 + 27.0 * t2 + 9.0 * t + 1.0) / 18.0;
  const w3 = (7.0 * t3 - 6.0 * t2) / 18.0;
  const raw = w0 * p0 + w1 * p1 + w2 * p2 + w3 * p3;
  const lo = Math.min(p0, p1, p2, p3);
  const hi = Math.max(p0, p1, p2, p3);
  return Math.max(lo, Math.min(hi, raw));
}

/**
 * Smoothly sample a 2D scalar field on the `48 x 48` PDE grid at world coordinates `(xMm, yMm)`
 * using convex-hull-clamped Bicubic Mitchell-Netravali (B=1/3, C=1/3) interpolation so free-surface
 * molten & solidified gold pools render with organic Young-Laplace meniscus contours, full peak height,
 * and continuous specular curvature.
 */
function sampleSmoothPdeField(
  gridData: Float32Array,
  xMm: number,
  yMm: number,
  channelOffset: number,
): number {
  const gx = (xMm + GRID_HALF_SPAN_MM) / GRID_DX_MM - 0.5;
  const gy = (yMm + GRID_HALF_SPAN_MM) / GRID_DX_MM - 0.5;
  const ix0 = Math.floor(gx);
  const iy0 = Math.floor(gy);
  const fx = gx - ix0;
  const fy = gy - iy0;
  const ch = channelOffset | 0;

  const readCell = (cx: number, cy: number): number => {
    const xClamp = Math.max(0, Math.min(GRID_N - 1, cx));
    const yClamp = Math.max(0, Math.min(GRID_N - 1, cy));
    return gridData[(yClamp * GRID_N + xClamp) * 7 + ch] ?? 0.0;
  };

  // Evaluate 4x4 Bicubic Mitchell-Netravali patch (convex-hull bounded, zero ringing)
  const rowVals = [0, 0, 0, 0];
  for (let m = -1; m <= 2; m++) {
    const v0 = readCell(ix0 - 1, iy0 + m);
    const v1 = readCell(ix0, iy0 + m);
    const v2 = readCell(ix0 + 1, iy0 + m);
    const v3 = readCell(ix0 + 2, iy0 + m);
    rowVals[m + 1] = mitchellNetravaliClamped(v0, v1, v2, v3, fx);
  }
  const val = mitchellNetravaliClamped(rowVals[0]!, rowVals[1]!, rowVals[2]!, rowVals[3]!, fy);
  return ch === 0 ? Math.max(70.0, val) : Math.max(0.0, val);
}

/**
 * Generate high-DPI procedural horological surface finishing normal & roughness textures:
 * 1. Côtes de Genève (Geneva Stripes) for the 4 contoured movement bridges
 * 2. Perlage (overlapping circular engine-turned graining) for the stepped mainplate recesses
 * 3. Clous de Paris / Barleycorn Guilloché for the silvered exhibition dial
 * 4. Soleil (radial sunburst brushing) for the winding ratchet & crown wheels
 */
function buildHorologicalFinishTextures(): {
  genevaNormal: THREE.CanvasTexture;
  perlageNormal: THREE.CanvasTexture;
  guillocheNormal: THREE.CanvasTexture;
  sunburstNormal: THREE.CanvasTexture;
  softParticleTex: THREE.CanvasTexture;
} {
  // 1. Côtes de Genève (Geneva Stripes)
  const cGeneva = document.createElement("canvas");
  cGeneva.width = 512;
  cGeneva.height = 512;
  const gCtx = cGeneva.getContext("2d")!;
  const gImg = gCtx.createImageData(512, 512);
  for (let y = 0; y < 512; y++) {
    for (let x = 0; x < 512; x++) {
      const u = (x + y * 0.18) / 512;
      const ribPhase = (u * 14.0) % 1.0;
      const microBrush = Math.sin(y * 1.65 + Math.sin(x * 0.25) * 2.0) * 0.12;
      // Sharp rib boundary + gentle cylindrical camber across each Geneva stripe
      const nx = (ribPhase - 0.5) * 0.65 + (ribPhase < 0.04 || ribPhase > 0.96 ? -0.55 : 0.0);
      const ny = microBrush;
      const idx = (y * 512 + x) * 4;
      gImg.data[idx] = Math.round(Math.max(0, Math.min(255, (nx * 0.5 + 0.5) * 255)));
      gImg.data[idx + 1] = Math.round(Math.max(0, Math.min(255, (ny * 0.5 + 0.5) * 255)));
      gImg.data[idx + 2] = 242;
      gImg.data[idx + 3] = 255;
    }
  }
  gCtx.putImageData(gImg, 0, 0);
  const genevaNormal = new THREE.CanvasTexture(cGeneva);
  genevaNormal.wrapS = THREE.RepeatWrapping;
  genevaNormal.wrapT = THREE.RepeatWrapping;
  genevaNormal.repeat.set(2.5, 2.5);

  // 2. Perlage (overlapping circular graining)
  const cPerlage = document.createElement("canvas");
  cPerlage.width = 512;
  cPerlage.height = 512;
  const pCtx = cPerlage.getContext("2d")!;
  const pImg = pCtx.createImageData(512, 512);
  const cellPx = 32;
  for (let y = 0; y < 512; y++) {
    for (let x = 0; x < 512; x++) {
      let bestDist = 999;
      let bestDx = 0;
      let bestDy = 0;
      const gx0 = Math.floor(x / (cellPx * 0.65));
      const gy0 = Math.floor(y / (cellPx * 0.65));
      for (let oy = -1; oy <= 1; oy++) {
        for (let ox = -1; ox <= 1; ox++) {
          const cx = ((gx0 + ox + ((gy0 + oy) & 1 ? 0.5 : 0.0)) * cellPx * 0.65) % 512;
          const cy = ((gy0 + oy) * cellPx * 0.65) % 512;
          let dx = x - cx;
          let dy = y - cy;
          if (dx > 256) dx -= 512;
          if (dx < -256) dx += 512;
          if (dy > 256) dy -= 512;
          if (dy < -256) dy += 512;
          const d = Math.hypot(dx, dy);
          if (d < cellPx * 0.55) {
            bestDist = d;
            bestDx = dx;
            bestDy = dy;
          }
        }
      }
      const rNorm = Math.max(1e-3, bestDist);
      const ringGroove = Math.sin(bestDist * 1.45) * 0.42;
      const tx = -bestDy / rNorm;
      const ty = bestDx / rNorm;
      const nx = tx * ringGroove;
      const ny = ty * ringGroove;
      const idx = (y * 512 + x) * 4;
      pImg.data[idx] = Math.round(Math.max(0, Math.min(255, (nx * 0.5 + 0.5) * 255)));
      pImg.data[idx + 1] = Math.round(Math.max(0, Math.min(255, (ny * 0.5 + 0.5) * 255)));
      pImg.data[idx + 2] = 240;
      pImg.data[idx + 3] = 255;
    }
  }
  pCtx.putImageData(pImg, 0, 0);
  const perlageNormal = new THREE.CanvasTexture(cPerlage);
  perlageNormal.wrapS = THREE.RepeatWrapping;
  perlageNormal.wrapT = THREE.RepeatWrapping;
  perlageNormal.repeat.set(3.5, 3.5);

  // 3. Clous de Paris / Barleycorn Guilloché
  const cGuill = document.createElement("canvas");
  cGuill.width = 512;
  cGuill.height = 512;
  const uCtx = cGuill.getContext("2d")!;
  const uImg = uCtx.createImageData(512, 512);
  for (let y = 0; y < 512; y++) {
    for (let x = 0; x < 512; x++) {
      const k = (2 * Math.PI * 36) / 512;
      const nx = Math.cos((x + y) * k) * 0.45;
      const ny = Math.cos((x - y) * k) * 0.45;
      const idx = (y * 512 + x) * 4;
      uImg.data[idx] = Math.round((nx * 0.5 + 0.5) * 255);
      uImg.data[idx + 1] = Math.round((ny * 0.5 + 0.5) * 255);
      uImg.data[idx + 2] = 235;
      uImg.data[idx + 3] = 255;
    }
  }
  uCtx.putImageData(uImg, 0, 0);
  const guillocheNormal = new THREE.CanvasTexture(cGuill);
  guillocheNormal.wrapS = THREE.RepeatWrapping;
  guillocheNormal.wrapT = THREE.RepeatWrapping;
  guillocheNormal.repeat.set(2, 2);

  // 4. Soleil (Radial Sunburst Brushing)
  const cSun = document.createElement("canvas");
  cSun.width = 256;
  cSun.height = 256;
  const sCtx = cSun.getContext("2d")!;
  const sImg = sCtx.createImageData(256, 256);
  for (let y = 0; y < 256; y++) {
    for (let x = 0; x < 256; x++) {
      const dx = x - 128;
      const dy = y - 128;
      const ang = Math.atan2(dy, dx);
      const ray = Math.sin(ang * 96) * 0.42;
      const nx = -Math.sin(ang) * ray;
      const ny = Math.cos(ang) * ray;
      const idx = (y * 256 + x) * 4;
      sImg.data[idx] = Math.round((nx * 0.5 + 0.5) * 255);
      sImg.data[idx + 1] = Math.round((ny * 0.5 + 0.5) * 255);
      sImg.data[idx + 2] = 240;
      sImg.data[idx + 3] = 255;
    }
  }
  sCtx.putImageData(sImg, 0, 0);
  const sunburstNormal = new THREE.CanvasTexture(cSun);

  // 5. Soft Radial Gaussian Sprite Texture (prevents square WebGL point-sprite artifacts)
  const cPart = document.createElement("canvas");
  cPart.width = 64;
  cPart.height = 64;
  const ptCtx = cPart.getContext("2d")!;
  const radGrad = ptCtx.createRadialGradient(32, 32, 1, 32, 32, 31);
  radGrad.addColorStop(0.0, "rgba(255, 255, 255, 1.0)");
  radGrad.addColorStop(0.35, "rgba(235, 248, 255, 0.65)");
  radGrad.addColorStop(0.72, "rgba(195, 232, 255, 0.18)");
  radGrad.addColorStop(1.0, "rgba(180, 225, 255, 0.0)");
  ptCtx.fillStyle = radGrad;
  ptCtx.fillRect(0, 0, 64, 64);
  const softParticleTex = new THREE.CanvasTexture(cPart);

  return { genevaNormal, perlageNormal, guillocheNormal, sunburstNormal, softParticleTex };
}

/**
 * Create an involute/epicycloidal horological gear wheel with realistic tooth profile,
 * crossed spokes, hub boss, and steel pinion leaf.
 */
function buildToothedGearGeometry(opts: {
  outerRadius: number;
  innerRimRadius: number;
  hubRadius: number;
  teeth: number;
  thickness: number;
  spokes: number;
  clubTeeth?: boolean;
}): THREE.ExtrudeGeometry {
  const { outerRadius, innerRimRadius, hubRadius, teeth, thickness, spokes, clubTeeth } = opts;
  const rootR = outerRadius * (clubTeeth ? 0.84 : 0.905);
  const pitchR = outerRadius * (clubTeeth ? 0.92 : 0.955);
  const shape = new THREE.Shape();
  const step = (Math.PI * 2) / teeth;

  for (let i = 0; i < teeth; i++) {
    const a0 = i * step;
    if (clubTeeth) {
      // Authentic Swiss lever club-tooth escape wheel profile: undercut stalk + angled impulse plane at tip
      const aRoot1 = a0;
      const aStalk = a0 + step * 0.24;
      const aHeel = a0 + step * 0.30;
      const aToe = a0 + step * 0.48;
      const aBack = a0 + step * 0.56;
      const p0x = Math.cos(aRoot1) * rootR;
      const p0y = Math.sin(aRoot1) * rootR;
      if (i === 0) shape.moveTo(p0x, p0y);
      else shape.lineTo(p0x, p0y);
      shape.lineTo(Math.cos(aStalk) * rootR, Math.sin(aStalk) * rootR);
      shape.lineTo(Math.cos(aHeel) * pitchR, Math.sin(aHeel) * pitchR);
      shape.lineTo(Math.cos(aToe) * outerRadius, Math.sin(aToe) * outerRadius);
      shape.lineTo(Math.cos(aBack) * (rootR * 0.98), Math.sin(aBack) * (rootR * 0.98));
    } else {
      const aRoot1 = a0 + step * 0.08;
      const aPitch1 = a0 + step * 0.24;
      const aTip1 = a0 + step * 0.34;
      const aTip2 = a0 + step * 0.50;
      const aPitch2 = a0 + step * 0.60;
      const aRoot2 = a0 + step * 0.76;
      const p0x = Math.cos(a0) * rootR;
      const p0y = Math.sin(a0) * rootR;
      if (i === 0) shape.moveTo(p0x, p0y);
      else shape.lineTo(p0x, p0y);
      shape.lineTo(Math.cos(aRoot1) * rootR, Math.sin(aRoot1) * rootR);
      shape.lineTo(Math.cos(aPitch1) * pitchR, Math.sin(aPitch1) * pitchR);
      shape.lineTo(Math.cos(aTip1) * outerRadius, Math.sin(aTip1) * outerRadius);
      shape.lineTo(Math.cos(aTip2) * outerRadius, Math.sin(aTip2) * outerRadius);
      shape.lineTo(Math.cos(aPitch2) * pitchR, Math.sin(aPitch2) * pitchR);
      shape.lineTo(Math.cos(aRoot2) * rootR, Math.sin(aRoot2) * rootR);
    }
  }
  shape.closePath();

  // Cut out spoked sector windows between hubRadius and innerRimRadius
  if (spokes >= 3 && innerRimRadius > hubRadius + 0.35) {
    const spokeSpan = (Math.PI * 2) / spokes;
    const halfSpokeAngle = Math.max(0.06, 0.14 / Math.max(1.5, innerRimRadius * 0.25));
    for (let s = 0; s < spokes; s++) {
      const wStart = s * spokeSpan + halfSpokeAngle;
      const wEnd = (s + 1) * spokeSpan - halfSpokeAngle;
      if (wEnd <= wStart + 0.08) continue;
      const hole = new THREE.Path();
      const segs = 12;
      hole.moveTo(Math.cos(wStart) * hubRadius, Math.sin(wStart) * hubRadius);
      for (let k = 0; k <= segs; k++) {
        const a = wStart + (wEnd - wStart) * (k / segs);
        hole.lineTo(Math.cos(a) * innerRimRadius, Math.sin(a) * innerRimRadius);
      }
      for (let k = segs; k >= 0; k--) {
        const a = wStart + (wEnd - wStart) * (k / segs);
        hole.lineTo(Math.cos(a) * hubRadius, Math.sin(a) * hubRadius);
      }
      hole.closePath();
      shape.holes.push(hole);
    }
  }

  // Center arbor bore hole
  const bore = new THREE.Path();
  bore.absarc(0, 0, Math.max(0.18, hubRadius * 0.35), 0, Math.PI * 2, true);
  shape.holes.push(bore);

  const geom = new THREE.ExtrudeGeometry(shape, {
    depth: thickness,
    bevelEnabled: true,
    bevelThickness: Math.min(0.045, thickness * 0.16),
    bevelSize: Math.min(0.035, thickness * 0.14),
    bevelSegments: 2,
    curveSegments: 16,
  });
  geom.translate(0, 0, -thickness * 0.5);
  geom.computeVertexNormals();
  return geom;
}

/**
 * Create a slotted horological screw (domed polished head + screwdriver slot) for bridges & chatons.
 */
function createSlottedWatchScrew(
  radiusMm: number,
  heightMm: number,
  headMat: THREE.Material,
  slotMat: THREE.Material,
): THREE.Group {
  const g = new THREE.Group();
  const head = new THREE.Mesh(
    new THREE.CylinderGeometry(radiusMm, radiusMm * 0.94, heightMm, 20),
    headMat,
  );
  head.rotation.x = Math.PI * 0.5;
  g.add(head);

  const slot = new THREE.Mesh(
    new THREE.BoxGeometry(radiusMm * 1.85, radiusMm * 0.24, heightMm * 0.48),
    slotMat,
  );
  slot.position.z = heightMm * 0.32;
  g.add(slot);
  return g;
}

/**
 * Create a 3-screw 18k Gold Chaton + Synthetic Corundum Ruby Jewel Bearing + Steel Arbor Pivot.
 */
function createRubyChatonAssembly(
  outerRadiusMm: number,
  rubyRadiusMm: number,
  goldMat: THREE.Material,
  rubyMat: THREE.Material,
  steelMat: THREE.Material,
  slotMat: THREE.Material,
  withIncablocLyre = false,
): THREE.Group {
  const g = new THREE.Group();

  // Raised beveled 18k gold chaton setting
  const chaton = new THREE.Mesh(
    new THREE.CylinderGeometry(outerRadiusMm * 0.88, outerRadiusMm, 0.34, 28),
    goldMat,
  );
  chaton.rotation.x = Math.PI * 0.5;
  g.add(chaton);

  // Convex chamfered synthetic ruby jewel with olived oil sink
  const ruby = new THREE.Mesh(
    new THREE.CylinderGeometry(rubyRadiusMm, rubyRadiusMm * 0.95, 0.38, 28),
    rubyMat,
  );
  ruby.rotation.x = Math.PI * 0.5;
  ruby.position.z = 0.03;
  g.add(ruby);

  // Polished steel pivot tip gleaming inside the ruby hole
  const pivot = new THREE.Mesh(
    new THREE.CylinderGeometry(rubyRadiusMm * 0.32, rubyRadiusMm * 0.32, 0.44, 16),
    steelMat,
  );
  pivot.rotation.x = Math.PI * 0.5;
  pivot.position.z = 0.05;
  g.add(pivot);

  // 3 micro-retaining screws around the gold chaton rim
  for (let k = 0; k < 3; k++) {
    const a = (k * Math.PI * 2) / 3 + 0.35;
    const sc = createSlottedWatchScrew(outerRadiusMm * 0.22, 0.22, steelMat, slotMat);
    sc.position.set(
      Math.cos(a) * outerRadiusMm * 1.16,
      Math.sin(a) * outerRadiusMm * 1.16,
      0.12,
    );
    sc.rotation.z = a;
    g.add(sc);
  }

  // Optional Incabloc lyre-shaped brass shock absorber spring over the balance endstone
  if (withIncablocLyre) {
    const lyreShape = new THREE.Shape();
    lyreShape.absarc(0, 0, outerRadiusMm * 0.92, -Math.PI * 0.38, Math.PI * 1.38, false);
    const hole = new THREE.Path();
    hole.absarc(0, 0, outerRadiusMm * 0.66, 0, Math.PI * 2, true);
    lyreShape.holes.push(hole);
    const lyreMesh = new THREE.Mesh(
      new THREE.ExtrudeGeometry(lyreShape, { depth: 0.12, bevelEnabled: false }),
      goldMat,
    );
    lyreMesh.position.z = 0.22;
    g.add(lyreMesh);
  }

  return g;
}

/**
 * Construct a 3D rectangular-cross-section spiral ribbon geometry for the 18-coil Mainspring
 * or 14-coil Breguet Overcoil Hairspring.
 */
function buildSpiralRibbonGeometry(
  samples: number,
  ribbonHeightMm: number,
  ribbonThicknessMm: number,
): THREE.BufferGeometry {
  // Each sample station has 4 vertices (top-inner, top-outer, bottom-outer, bottom-inner)
  const vertexCount = samples * 4;
  const positions = new Float32Array(vertexCount * 3);
  const normals = new Float32Array(vertexCount * 3);
  const indices: number[] = [];

  for (let i = 0; i < samples - 1; i++) {
    const b0 = i * 4;
    const b1 = (i + 1) * 4;
    // Top face
    indices.push(b0, b0 + 1, b1 + 1, b0, b1 + 1, b1);
    // Outer face
    indices.push(b0 + 1, b0 + 2, b1 + 2, b0 + 1, b1 + 2, b1 + 1);
    // Bottom face
    indices.push(b0 + 2, b0 + 3, b1 + 3, b0 + 2, b1 + 3, b1 + 2);
    // Inner face
    indices.push(b0 + 3, b0, b1, b0 + 3, b1, b1 + 3);
  }

  const geom = new THREE.BufferGeometry();
  geom.setAttribute("position", new THREE.BufferAttribute(positions, 3));
  geom.setAttribute("normal", new THREE.BufferAttribute(normals, 3));
  geom.setIndex(indices);
  geom.userData = { samples, ribbonHeightMm, ribbonThicknessMm };
  return geom;
}

/**
 * Update a 3D spiral ribbon geometry along a centerline curve `(x(i), y(i), z(i))` with radial normal `(nx(i), ny(i))`.
 */
function updateSpiralRibbonVertices(
  geom: THREE.BufferGeometry,
  evalPoint: (i: number, samples: number) => { x: number; y: number; z: number; nx: number; ny: number },
): void {
  const { samples, ribbonHeightMm, ribbonThicknessMm } = geom.userData as {
    samples: number;
    ribbonHeightMm: number;
    ribbonThicknessMm: number;
  };
  const posAttr = geom.getAttribute("position") as THREE.BufferAttribute;
  const nrmAttr = geom.getAttribute("normal") as THREE.BufferAttribute;
  const arr = posAttr.array as Float32Array;
  const nrm = nrmAttr.array as Float32Array;
  const halfH = ribbonHeightMm * 0.5;
  const halfT = ribbonThicknessMm * 0.5;
  // Analytic corner normals, weighted like an area-weighted average of the strip's side face
  // (height h) and its top/bottom face (thickness t), without a per-frame computeVertexNormals pass.
  const wSide = ribbonHeightMm / Math.hypot(ribbonHeightMm, ribbonThicknessMm);
  const wCap = ribbonThicknessMm / Math.hypot(ribbonHeightMm, ribbonThicknessMm);

  for (let i = 0; i < samples; i++) {
    const { x, y, z, nx, ny } = evalPoint(i, samples);
    const xi = x - nx * halfT;
    const yi = y - ny * halfT;
    const xo = x + nx * halfT;
    const yo = y + ny * halfT;

    const base = i * 12;
    nrm[base] = -nx * wSide;
    nrm[base + 1] = -ny * wSide;
    nrm[base + 2] = wCap;
    nrm[base + 3] = nx * wSide;
    nrm[base + 4] = ny * wSide;
    nrm[base + 5] = wCap;
    nrm[base + 6] = nx * wSide;
    nrm[base + 7] = ny * wSide;
    nrm[base + 8] = -wCap;
    nrm[base + 9] = -nx * wSide;
    nrm[base + 10] = -ny * wSide;
    nrm[base + 11] = -wCap;
    // 0: top-inner
    arr[base] = xi;
    arr[base + 1] = yi;
    arr[base + 2] = z + halfH;
    // 1: top-outer
    arr[base + 3] = xo;
    arr[base + 4] = yo;
    arr[base + 5] = z + halfH;
    // 2: bottom-outer
    arr[base + 6] = xo;
    arr[base + 7] = yo;
    arr[base + 8] = z - halfH;
    // 3: bottom-inner
    arr[base + 9] = xi;
    arr[base + 10] = yi;
    arr[base + 11] = z - halfH;
  }
  posAttr.needsUpdate = true;
  nrmAttr.needsUpdate = true;
  if (!geom.boundingSphere) geom.computeBoundingSphere();
}

/**
 * Build a unit-length 3D arrow glyph oriented along +Y (y in [0, 1]) combining a tapered
 * cylindrical shaft (y in [0, 0.66]) and a sharp conical arrowhead (y in [0.62, 1.0])
 * for GPU-instanced 3D vector field visualization.
 */
function buildUnitArrowGeometry(): THREE.BufferGeometry {
  const radialSegs = 12;
  const rings: Array<{ y: number; r: number }> = [
    { y: 0.0, r: 0.0 },
    { y: 0.0, r: 0.068 },
    { y: 0.64, r: 0.052 },
    { y: 0.64, r: 0.175 },
    { y: 1.0, r: 0.0 },
  ];
  const vertsPerRing = radialSegs + 1;
  const positions = new Float32Array(rings.length * vertsPerRing * 3);
  const indices: number[] = [];

  for (let ir = 0; ir < rings.length; ir++) {
    const { y, r } = rings[ir]!;
    for (let is = 0; is <= radialSegs; is++) {
      const a = (is / radialSegs) * Math.PI * 2;
      const idx = (ir * vertsPerRing + is) * 3;
      positions[idx] = Math.cos(a) * r;
      positions[idx + 1] = y;
      positions[idx + 2] = Math.sin(a) * r;
    }
  }
  for (let ir = 0; ir < rings.length - 1; ir++) {
    for (let is = 0; is < radialSegs; is++) {
      const a = ir * vertsPerRing + is;
      const b = a + 1;
      const c = (ir + 1) * vertsPerRing + is;
      const d = c + 1;
      indices.push(a, c, b, b, c, d);
    }
  }
  const geom = new THREE.BufferGeometry();
  geom.setAttribute("position", new THREE.BufferAttribute(positions, 3));
  geom.setIndex(indices);
  geom.computeVertexNormals();
  return geom;
}

const VEC_FIELD_LATTICE_COUNT = 168;
const VEC_UP = new THREE.Vector3(0, 1, 0);

export class PocketWatchScene {
  private readonly canvas: HTMLCanvasElement;
  private readonly callbacks: SceneCallbacks;
  private readonly renderer: THREE.WebGLRenderer;
  private readonly composer: EffectComposer;
  private readonly bloomPass: UnrealBloomPass;
  private readonly smaaPass: SMAAPass;
  private readonly scene: THREE.Scene;
  private readonly camera: THREE.PerspectiveCamera;
  private readonly controls: OrbitControls;
  private readonly raycaster = new THREE.Raycaster();
  private readonly pointerNdc = new THREE.Vector2();
  private readonly interactionPlane = new THREE.Plane(new THREE.Vector3(0, 0, 1), 0);

  private activeTool: ActiveToolMode = "inspect";
  private isDraggingTool = false;
  private explodedFactor = 0.0;
  private openCaseMode = true;
  private thermalOverlayMode = false;

  // Root layer groups for exploded-view axial separation
  private readonly watchRoot = new THREE.Group();
  private readonly layerCase = new THREE.Group();
  private readonly layerMainplate = new THREE.Group();
  private readonly layerKeyless = new THREE.Group();
  private readonly layerBarrel = new THREE.Group();
  private readonly layerGoingTrain = new THREE.Group();
  private readonly layerBridges = new THREE.Group();
  private readonly layerTourbillon = new THREE.Group();
  private readonly layerDialHands = new THREE.Group();
  private readonly layerCrystalFront = new THREE.Group();

  // Dynamic kinematic groups driven by WASM state
  private readonly crownGroup = new THREE.Group();
  private readonly crownWheelGroup = new THREE.Group();
  private readonly ratchetWheelGroup = new THREE.Group();
  private readonly clickLeverGroup = new THREE.Group();
  private readonly barrelWheelGroup = new THREE.Group();
  private readonly centerWheelGroup = new THREE.Group();
  private readonly thirdWheelGroup = new THREE.Group();
  private readonly tourbillonCageGroup = new THREE.Group();
  private readonly escapeWheelGroup = new THREE.Group();
  private readonly palletForkGroup = new THREE.Group();
  private readonly balanceWheelGroup = new THREE.Group();
  private readonly hourHandGroup = new THREE.Group();
  private readonly minuteHandGroup = new THREE.Group();
  private readonly secondsHandGroup = new THREE.Group();
  private readonly powerReserveHandGroup = new THREE.Group();
  private readonly moonDiscGroup = new THREE.Group();

  // Emergent Continuum Meshes driven directly by WASM arrays
  private caseShellGeometry!: THREE.BufferGeometry;
  private caseShellMesh!: THREE.Mesh;
  private moltenGoldGeometry!: THREE.BufferGeometry;
  private moltenGoldMesh!: THREE.Mesh;
  private solidGoldGeometry!: THREE.BufferGeometry;
  private solidGoldMesh!: THREE.Mesh;
  private mainspringRibbonGeometry!: THREE.BufferGeometry;
  private hairspringRibbonGeometry!: THREE.BufferGeometry;

  // Interactive perturbation 3D apparatus
  private readonly torchGroup = new THREE.Group();
  private readonly torchFlameCore!: THREE.Mesh;
  private readonly torchFlamePlume!: THREE.Mesh;
  private readonly torchPointLight!: THREE.PointLight;
  private readonly ln2Group = new THREE.Group();
  private readonly ln2StreamMesh!: THREE.Mesh;
  private readonly ln2FrostCloud!: THREE.Points;
  private readonly oilGroup = new THREE.Group();
  private readonly oilStreamMesh!: THREE.Mesh;
  private readonly magnetGroup = new THREE.Group();
  private readonly airGroup = new THREE.Group();
  private readonly airStreamMesh!: THREE.Mesh;
  private readonly micGroup = new THREE.Group();

  // Superimposed 3D Physics Vector Fields (Heat, Convection, Magnetic, Acoustic, Stress/Gravity)
  private readonly vectorFieldsGroup = new THREE.Group();
  private vecHeatMesh!: THREE.InstancedMesh;
  private vecConvectionMesh!: THREE.InstancedMesh;
  private vecMagneticMesh!: THREE.InstancedMesh;
  private vecAcousticMesh!: THREE.InstancedMesh;
  private vecStressMesh!: THREE.InstancedMesh;
  private readonly vecFieldVisible: Record<VectorFieldLayer, boolean> = {
    heat: false,
    convection: false,
    magnetic: false,
    acoustic: false,
    stress: false,
  };
  private readonly vecLatticePts: Array<[number, number, number]> = [];
  private readonly tmpMat4 = new THREE.Matrix4();
  private readonly tmpPos = new THREE.Vector3();
  private readonly tmpDir = new THREE.Vector3();
  private readonly tmpQuat = new THREE.Quaternion();
  private readonly tmpScale = new THREE.Vector3();
  private readonly tmpColor = new THREE.Color();

  // Shared PBR Materials
  private readonly matGoldCase: THREE.MeshPhysicalMaterial;
  private readonly matSolidGoldStatic: THREE.MeshPhysicalMaterial;
  private readonly matMoltenGold: THREE.MeshPhysicalMaterial;
  private readonly matSolidifiedGold: THREE.MeshPhysicalMaterial;
  private readonly matRhodiumPlate: THREE.MeshStandardMaterial;
  private readonly matGenevaBridge: THREE.MeshStandardMaterial;
  private readonly matBrassWheel: THREE.MeshStandardMaterial;
  private readonly matPolishedSteel: THREE.MeshStandardMaterial;
  private readonly matSunburstSteel: THREE.MeshStandardMaterial;
  private readonly matBluedSteel: THREE.MeshPhysicalMaterial;
  private readonly matHairspring: THREE.MeshPhysicalMaterial;
  private readonly matTitaniumCage: THREE.MeshStandardMaterial;
  private readonly matGlucydurBalance: THREE.MeshStandardMaterial;
  private readonly matRubyJewel: THREE.MeshPhysicalMaterial;
  private readonly matSapphireCrystal: THREE.MeshPhysicalMaterial;
  private readonly matScrewSlot: THREE.MeshBasicMaterial;
  private readonly matFrostOverlay: THREE.MeshStandardMaterial;
  private readonly frostOverlayMesh!: THREE.Mesh;

  // GPU continuum fields: blackbody LUT, case-shell node textures, spatial PDE grid texture
  private readonly glowUniforms: {
    uBlackbody: { value: THREE.DataTexture };
    uGlowScale: { value: number };
    uThermalOverlay: { value: number };
  };
  private readonly shellBufA = new Float32Array(CASE_SHELL_NODES * 4);
  private readonly shellBufB = new Float32Array(CASE_SHELL_NODES * 4);
  private readonly shellUniforms: {
    uShellA: { value: THREE.DataTexture };
    uShellB: { value: THREE.DataTexture };
  };
  private shellChecksum = Number.NaN;
  private gridTexBuf = new Float32Array(GRID_N * GRID_N * 4);
  private readonly gridUniforms: {
    uGrid: { value: THREE.DataTexture };
    uGridN: { value: number };
    uHalfSpan: { value: number };
  };
  private lastGridData: Float32Array | null = null;
  private lastCaseShellData: Float32Array | null = null;
  private maxShellTempK = 293.15;
  private maxGridTempK = 293.15;

  // Quality, framing, camera animation, and escapement contact highlights
  private tier: QualityTier;
  private readonly reducedMotion: boolean;
  private pixelRatio = 1;
  private readonly keyLight: THREE.DirectionalLight;
  private viewInsets: ViewInsets = { left: 0, right: 0, top: 0, bottom: 0 };
  private cssWidth = 960;
  private cssHeight = 640;
  private currentPreset: CameraPreset = "overview";
  private cameraTween: {
    fromPos: THREE.Vector3;
    toPos: THREE.Vector3;
    fromTarget: THREE.Vector3;
    toTarget: THREE.Vector3;
    t0: number;
    durationMs: number;
  } | null = null;
  private readonly contactMats: THREE.MeshPhysicalMaterial[] = [];
  private readonly contactGlow = new Float32Array(3);
  private aimPointerId: number | null = null;

  constructor(canvas: HTMLCanvasElement, callbacks: SceneCallbacks, opts: SceneOptions = {}) {
    this.canvas = canvas;
    this.callbacks = callbacks;
    this.tier = opts.tier ?? "high";
    this.reducedMotion = opts.reducedMotion ?? false;

    // The composer owns anti-aliasing (MSAA on the HDR target), so the default framebuffer stays
    // single-sampled: no wasted multisample buffer on phones.
    this.renderer = new THREE.WebGLRenderer({
      canvas,
      antialias: false,
      powerPreference: "high-performance",
    });
    this.pixelRatio = Math.min(window.devicePixelRatio || 1, PocketWatchScene.dprCap(this.tier));
    this.renderer.setPixelRatio(this.pixelRatio);
    // Khronos PBR Neutral keeps base colors faithful (gold stays gold instead of drifting orange).
    this.renderer.toneMapping = THREE.NeutralToneMapping;
    this.renderer.toneMappingExposure = 1.0;
    this.renderer.shadowMap.enabled = true;
    this.renderer.shadowMap.type = THREE.PCFSoftShadowMap;

    this.scene = new THREE.Scene();
    this.scene.background = new THREE.Color(0x05070b);

    const pmrem = new THREE.PMREMGenerator(this.renderer);
    const studio = buildStudioEnvironmentScene();
    this.scene.environment = pmrem.fromScene(studio, 0.02).texture;
    this.scene.environmentIntensity = 0.62;
    pmrem.dispose();
    studio.traverse((o) => {
      const m = o as THREE.Mesh;
      if (m.isMesh) {
        m.geometry.dispose();
        (m.material as THREE.Material).dispose();
      }
    });
    void RoomEnvironment;

    this.glowUniforms = {
      uBlackbody: { value: buildBlackbodyLut() },
      uGlowScale: { value: 0.75 },
      uThermalOverlay: { value: 0 },
    };
    const makeFloatTex = (buf: Float32Array, w: number, h: number): THREE.DataTexture => {
      const t = new THREE.DataTexture(buf, w, h, THREE.RGBAFormat, THREE.FloatType);
      t.minFilter = THREE.NearestFilter;
      t.magFilter = THREE.NearestFilter;
      t.generateMipmaps = false;
      t.needsUpdate = true;
      return t;
    };
    this.shellUniforms = {
      uShellA: { value: makeFloatTex(this.shellBufA, CASE_SHELL_N_THETA, CASE_SHELL_N_Z) },
      uShellB: { value: makeFloatTex(this.shellBufB, CASE_SHELL_N_THETA, CASE_SHELL_N_Z) },
    };
    for (let i = 0; i < GRID_N * GRID_N; i++) this.gridTexBuf[i * 4 + 2] = 293.15;
    this.gridUniforms = {
      uGrid: { value: makeFloatTex(this.gridTexBuf, GRID_N, GRID_N) },
      uGridN: { value: GRID_N },
      uHalfSpan: { value: GRID_HALF_SPAN_MM },
    };

    this.camera = new THREE.PerspectiveCamera(34, 1.0, 0.5, 450);
    this.camera.position.set(0, -34, 66);

    // Half-float HDR post chain (RenderPass -> UnrealBloomPass -> OutputPass). MSAA lives on the HDR
    // target; SMAA is only enabled on the low tier where the target is single-sampled. Bloom is
    // gated off whenever nothing is incandescent or burning, so idle frames skip its five blur mips.
    const initW = Math.max(320, canvas.clientWidth || 960);
    const initH = Math.max(260, canvas.clientHeight || 640);
    const hdrTarget = new THREE.WebGLRenderTarget(initW * this.pixelRatio, initH * this.pixelRatio, {
      type: THREE.HalfFloatType,
      format: THREE.RGBAFormat,
      samples: this.tier === "high" ? 4 : this.tier === "medium" ? 2 : 0,
    });
    this.composer = new EffectComposer(this.renderer, hdrTarget);
    this.composer.addPass(new RenderPass(this.scene, this.camera));
    this.bloomPass = new UnrealBloomPass(new THREE.Vector2(initW, initH), 0.08, 0.36, 1.25);
    this.bloomPass.enabled = false;
    this.composer.addPass(this.bloomPass);
    this.smaaPass = new SMAAPass();
    this.smaaPass.enabled = this.tier === "low";
    this.composer.addPass(this.smaaPass);
    this.composer.addPass(new OutputPass());

    this.controls = new OrbitControls(this.camera, this.canvas);
    this.controls.enableDamping = true;
    this.controls.dampingFactor = 0.07;
    this.controls.rotateSpeed = 0.8;
    this.controls.zoomSpeed = 0.9;
    this.controls.minDistance = 4;
    this.controls.maxDistance = 180;
    this.controls.target.set(0, -2, 0);
    this.controls.addEventListener("start", () => {
      this.cameraTween = null;
    });

    // Museum-grade chiaroscuro horological lighting (balanced to prevent metallic blowout)
    const keyLight = new THREE.DirectionalLight(0xfff4e0, 1.15);
    this.keyLight = keyLight;
    keyLight.position.set(32, 42, 58);
    keyLight.castShadow = true;
    const shadowRes = this.tier === "low" ? 1024 : 2048;
    keyLight.shadow.mapSize.set(shadowRes, shadowRes);
    keyLight.shadow.camera.near = 10;
    keyLight.shadow.camera.far = 160;
    const d = 34;
    keyLight.shadow.camera.left = -d;
    keyLight.shadow.camera.right = d;
    keyLight.shadow.camera.top = d;
    keyLight.shadow.camera.bottom = -d;
    keyLight.shadow.bias = -0.0005;
    this.scene.add(keyLight);

    const rimLight = new THREE.DirectionalLight(0x8ab6f9, 0.58);
    rimLight.position.set(-38, -28, 36);
    this.scene.add(rimLight);

    const underFill = new THREE.DirectionalLight(0xffd79e, 0.28);
    underFill.position.set(0, 24, -38);
    this.scene.add(underFill);

    // Procedural Horological Finishes
    const finishes = buildHorologicalFinishTextures();

    // Physically based conductor reflectances (linear-sRGB F0 at normal incidence, from measured
    // n,k data: Au/Ag/Cu 18k yellow alloy, CuZn brass, rhodium plating, carbon steel, Ti-6Al-4V,
    // CuBe2 Glucydur). Dielectrics use their real IOR: corundum (ruby/sapphire) n = 1.768.
    const linear = (r: number, g: number, b: number): THREE.Color =>
      new THREE.Color().setRGB(r, g, b, THREE.LinearSRGBColorSpace);
    const lowTier = this.tier === "low";

    // 18k case: base F0 lives in `color`; the shader multiplies thermal tints / IR false-color and
    // adds per-pixel blackbody emission from the GPU-sampled shell temperature.
    this.matGoldCase = new THREE.MeshPhysicalMaterial({
      color: linear(0.95, 0.76, 0.4),
      metalness: 1.0,
      roughness: 0.2,
      clearcoat: 0.1,
      clearcoatRoughness: 0.12,
      specularIntensity: 1.0,
    });

    this.matSolidGoldStatic = new THREE.MeshPhysicalMaterial({
      color: linear(0.95, 0.76, 0.4),
      metalness: 1.0,
      roughness: 0.17,
      clearcoat: 0.1,
      clearcoatRoughness: 0.12,
    });

    // Liquid gold is a mirror (roughness ~0.05) whose base color whitens slightly when hot.
    this.matMoltenGold = new THREE.MeshPhysicalMaterial({
      color: new THREE.Color(1, 1, 1),
      metalness: 1.0,
      roughness: 0.06,
      clearcoat: 0.5,
      clearcoatRoughness: 0.05,
      transparent: true,
      alphaTest: 0.02,
      depthWrite: true,
    });

    this.matSolidifiedGold = new THREE.MeshPhysicalMaterial({
      color: new THREE.Color(1, 1, 1),
      metalness: 1.0,
      roughness: 0.34,
      clearcoat: 0.08,
      transparent: true,
      alphaTest: 0.02,
      depthWrite: true,
    });

    this.matRhodiumPlate = new THREE.MeshStandardMaterial({
      color: linear(0.66, 0.67, 0.68),
      metalness: 1.0,
      roughness: 0.32,
      normalMap: finishes.perlageNormal,
      normalScale: new THREE.Vector2(0.45, 0.45),
    });

    this.matGenevaBridge = new THREE.MeshStandardMaterial({
      color: linear(0.7, 0.71, 0.72),
      metalness: 1.0,
      roughness: 0.24,
      normalMap: finishes.genevaNormal,
      normalScale: new THREE.Vector2(0.55, 0.55),
    });

    this.matBrassWheel = new THREE.MeshStandardMaterial({
      color: linear(0.91, 0.73, 0.38),
      metalness: 1.0,
      roughness: 0.22,
    });

    this.matPolishedSteel = new THREE.MeshStandardMaterial({
      color: linear(0.58, 0.58, 0.6),
      metalness: 1.0,
      roughness: 0.12,
    });

    this.matSunburstSteel = new THREE.MeshStandardMaterial({
      color: linear(0.56, 0.57, 0.6),
      metalness: 1.0,
      roughness: 0.2,
      normalMap: finishes.sunburstNormal,
      normalScale: new THREE.Vector2(0.52, 0.52),
    });

    // Blued steel is a ~250 nm Fe3O4 temper film on steel: thin-film interference (iridescence)
    // over a steel substrate, tinted by the engine's fs-render interference color.
    this.matBluedSteel = new THREE.MeshPhysicalMaterial({
      color: linear(0.09, 0.2, 0.62),
      metalness: 0.9,
      roughness: 0.16,
      clearcoat: 0.4,
      clearcoatRoughness: 0.08,
      iridescence: lowTier ? 0 : 0.55,
      iridescenceIOR: 2.4,
      iridescenceThicknessRange: [220, 300],
    });

    this.matHairspring = new THREE.MeshPhysicalMaterial({
      color: linear(0.12, 0.25, 0.7),
      metalness: 0.9,
      roughness: 0.16,
      clearcoat: 0.35,
      iridescence: lowTier ? 0 : 0.45,
      iridescenceIOR: 2.4,
      iridescenceThicknessRange: [220, 300],
    });

    this.matTitaniumCage = new THREE.MeshStandardMaterial({
      color: linear(0.54, 0.5, 0.45),
      metalness: 1.0,
      roughness: 0.18,
    });

    this.matGlucydurBalance = new THREE.MeshStandardMaterial({
      color: linear(0.93, 0.66, 0.4),
      metalness: 1.0,
      roughness: 0.16,
    });

    // Synthetic ruby (Cr:Al2O3): transmissive corundum with Cr3+ absorption (green/blue absorbed).
    this.matRubyJewel = new THREE.MeshPhysicalMaterial({
      color: linear(1.0, 0.22, 0.3),
      emissive: linear(0.45, 0.0, 0.04),
      emissiveIntensity: 0.0,
      metalness: 0.0,
      roughness: 0.05,
      transmission: lowTier ? 0 : 0.9,
      ior: 1.768,
      thickness: 0.8,
      attenuationColor: linear(0.95, 0.06, 0.12),
      attenuationDistance: 0.45,
      specularIntensity: 1.0,
      transparent: lowTier,
      opacity: lowTier ? 0.9 : 1.0,
    });

    this.matSapphireCrystal = new THREE.MeshPhysicalMaterial({
      color: new THREE.Color(1, 1, 1),
      metalness: 0.0,
      roughness: 0.03,
      transmission: lowTier ? 0 : 1.0,
      ior: 1.768,
      thickness: 0.65,
      specularIntensity: 1.0,
      transparent: true,
      opacity: lowTier ? 0.14 : 1.0,
      depthWrite: false,
    });

    this.matScrewSlot = new THREE.MeshBasicMaterial({ color: 0x0b101b });

    this.matFrostOverlay = new THREE.MeshStandardMaterial({
      color: new THREE.Color(0xdff4ff),
      roughness: 0.85,
      metalness: 0.05,
      transparent: true,
      opacity: 0.0,
      depthWrite: false,
    });

    // Assemble the complete ~380-part Patek Philippe Grand Complication Pocket Watch
    this.watchRoot.add(
      this.layerCase,
      this.layerMainplate,
      this.layerKeyless,
      this.layerBarrel,
      this.layerGoingTrain,
      this.layerBridges,
      this.layerTourbillon,
      this.layerDialHands,
      this.layerCrystalFront,
    );
    this.scene.add(this.watchRoot);

    this.buildDeformableBassineCaseAndCrown();
    this.buildSteppedMainplateAndJewels();
    this.buildKeylessWorksAndRatchetTrain();
    this.buildMainspringBarrelAndElasticaRibbon();
    this.buildGoingTrainWheels();
    this.buildContouredGenevaBridges();
    this.buildTourbillonEscapementAndBalance();
    this.buildExhibitionGuillocheDialAndHands(finishes.guillocheNormal);
    this.buildContinuumMoltenAndSolidGoldHeightfields();

    // Frost rime overlay disc
    this.frostOverlayMesh = new THREE.Mesh(
      new THREE.CylinderGeometry(23.5, 23.5, 0.2, 48),
      this.matFrostOverlay,
    );
    this.frostOverlayMesh.rotation.x = Math.PI * 0.5;
    this.frostOverlayMesh.position.z = 2.65;
    this.frostOverlayMesh.visible = false;
    this.watchRoot.add(this.frostOverlayMesh);

    // Build Interactive Perturbation Tools (Blowtorch, LN2, Oil, Magnet, Air Jet, Studio Microphone)
    const tools = this.buildPerturbationApparatus(finishes.softParticleTex);
    this.torchFlameCore = tools.torchFlameCore;
    this.torchFlamePlume = tools.torchFlamePlume;
    this.torchPointLight = tools.torchPointLight;
    this.ln2StreamMesh = tools.ln2StreamMesh;
    this.ln2FrostCloud = tools.ln2FrostCloud;
    this.oilStreamMesh = tools.oilStreamMesh;
    this.airStreamMesh = tools.airStreamMesh;

    // Build Superimposed 3D Physics Vector Fields (Heat, Convection, Magnetic, Acoustic, Stress/Gravity)
    this.buildVectorFieldOverlays();

    this.bindPointerEvents();
    this.resize();
    window.addEventListener("resize", () => this.resize());
  }

  /**
   * 1. Emergent 3D Deformable 18k Gold Bassine Case Shell (`128 x 22` Smooth Catmull-Rom Subdivided Surface
   *    evaluated from the `64 x 8` Rust Lagrangian shell continuum nodes)
   *    + 12h Pendant, 24-Flute Onion Crown, Articulated Bow, and Sapphire Crystals.
   */
  private buildDeformableBassineCaseAndCrown(): void {
    const ringVerts = CASE_RENDER_N_THETA * CASE_RENDER_N_Z;
    const totalVerts = ringVerts * 2; // 0..ringVerts-1 = outer surface, ringVerts..2*ringVerts-1 = inner surface
    const positions = new Float32Array(totalVerts * 3);
    const colors = new Float32Array(totalVerts * 3);
    const indices: number[] = [];

    const profileZ = [-3.8, -2.8, -1.6, -0.5, 0.5, 1.6, 2.8, 3.8];
    const profileR = [23.0, 25.0, 26.3, 26.8, 26.8, 26.3, 25.0, 23.0];
    const profileH = [2.1, 2.3, 2.6, 2.8, 2.8, 2.6, 2.3, 2.1];

    // Initialize synthetic default 64x8 buffer and populate the smooth 128x22 surface
    const initData = new Float32Array(CASE_SHELL_NODES * 6);
    for (let iz = 0; iz < CASE_SHELL_N_Z; iz++) {
      for (let it = 0; it < CASE_SHELL_N_THETA; it++) {
        const i = (iz * CASE_SHELL_N_THETA + it) * 6;
        const theta = (it * Math.PI * 2) / CASE_SHELL_N_THETA;
        const dyTop = Math.max(0.0, Math.sin(theta));
        const dAng = Math.abs(theta - Math.PI * 0.5);
        const pendantLobe =
          dAng < 0.22
            ? 2.2 * (1.0 - dAng / 0.22) * (1.0 - dAng / 0.22) * (1.0 - 0.12 * Math.abs(profileZ[iz]!))
            : 0.0;
        const r0 = profileR[iz]! + pendantLobe * dyTop;
        initData[i] = r0 * Math.cos(theta);
        initData[i + 1] = r0 * Math.sin(theta);
        initData[i + 2] = profileZ[iz]!;
        initData[i + 3] = profileH[iz]! + 0.35 * pendantLobe;
        initData[i + 4] = 293.15;
        initData[i + 5] = 0.0;
      }
    }

    // Outer wall quads + Inner wall quads on the subdivided 128 x 22 lattice
    for (let iz = 0; iz < CASE_RENDER_N_Z - 1; iz++) {
      for (let it = 0; it < CASE_RENDER_N_THETA; it++) {
        const itNext = (it + 1) % CASE_RENDER_N_THETA;
        const o00 = iz * CASE_RENDER_N_THETA + it;
        const o10 = iz * CASE_RENDER_N_THETA + itNext;
        const o01 = (iz + 1) * CASE_RENDER_N_THETA + it;
        const o11 = (iz + 1) * CASE_RENDER_N_THETA + itNext;
        indices.push(o00, o10, o11, o00, o11, o01);

        const i00 = ringVerts + o00;
        const i10 = ringVerts + o10;
        const i01 = ringVerts + o01;
        const i11 = ringVerts + o11;
        indices.push(i00, i11, i10, i00, i01, i11);
      }
    }

    // Top bezel rim cap & Bottom bezel rim cap
    for (let it = 0; it < CASE_RENDER_N_THETA; it++) {
      const itNext = (it + 1) % CASE_RENDER_N_THETA;
      const topO0 = (CASE_RENDER_N_Z - 1) * CASE_RENDER_N_THETA + it;
      const topO1 = (CASE_RENDER_N_Z - 1) * CASE_RENDER_N_THETA + itNext;
      const topI0 = ringVerts + topO0;
      const topI1 = ringVerts + topO1;
      indices.push(topO0, topO1, topI1, topO0, topI1, topI0);

      const botO0 = it;
      const botO1 = itNext;
      const botI0 = ringVerts + botO0;
      const botI1 = ringVerts + botO1;
      indices.push(botO0, botI1, botO1, botO0, botI0, botI1);
    }

    this.caseShellGeometry = new THREE.BufferGeometry();
    this.caseShellGeometry.setAttribute("position", new THREE.BufferAttribute(positions, 3));
    this.caseShellGeometry.setAttribute("color", new THREE.BufferAttribute(colors, 3));
    this.caseShellGeometry.setIndex(indices);
    this.updateCaseShellFromContinuum(initData);

    this.caseShellMesh = new THREE.Mesh(this.caseShellGeometry, this.matGoldCase);
    this.caseShellMesh.castShadow = true;
    this.caseShellMesh.receiveShadow = true;
    this.layerCase.add(this.caseShellMesh);

    // 12 o'clock Pendant Post & Collar
    const pendantPost = new THREE.Mesh(
      new THREE.CylinderGeometry(2.4, 3.0, 5.2, 32),
      this.matSolidGoldStatic,
    );
    pendantPost.position.set(0, 28.4, 0);
    this.layerCase.add(pendantPost);

    // 24-Flute Onion Winding Crown
    this.crownGroup.position.set(0, 32.0, 0);
    const crownCore = new THREE.Mesh(
      new THREE.SphereGeometry(3.4, 32, 20),
      this.matSolidGoldStatic,
    );
    crownCore.scale.set(1.0, 0.72, 1.0);
    this.crownGroup.add(crownCore);

    for (let k = 0; k < 24; k++) {
      const ang = (k * Math.PI * 2) / 24;
      const flute = new THREE.Mesh(
        new THREE.CylinderGeometry(0.32, 0.32, 2.2, 8),
        this.matSolidGoldStatic,
      );
      flute.position.set(Math.cos(ang) * 3.25, 0, Math.sin(ang) * 3.25);
      this.crownGroup.add(flute);
    }
    this.layerCase.add(this.crownGroup);

    // Articulated 18k Gold Pocket-Watch Bow (Torus ring arching over the onion crown)
    const bow = new THREE.Mesh(
      new THREE.TorusGeometry(6.2, 0.72, 20, 64, Math.PI * 1.48),
      this.matSolidGoldStatic,
    );
    bow.position.set(0, 30.2, 0);
    bow.rotation.z = -Math.PI * 0.24;
    this.layerCase.add(bow);

    // Rear Exhibition Sapphire Crystal
    const backCrystal = new THREE.Mesh(
      new THREE.CylinderGeometry(22.4, 22.4, 0.65, 64),
      this.matSapphireCrystal,
    );
    backCrystal.rotation.x = Math.PI * 0.5;
    backCrystal.position.z = -3.45;
    this.layerCase.add(backCrystal);

    // Front Exhibition Sapphire Crystal + Thin Bezel Ring (toggled via openCaseMode)
    const frontCrystal = new THREE.Mesh(
      new THREE.CylinderGeometry(22.6, 22.6, 0.65, 64),
      this.matSapphireCrystal,
    );
    frontCrystal.rotation.x = Math.PI * 0.5;
    frontCrystal.position.z = 3.55;
    this.layerCrystalFront.add(frontCrystal);
    this.layerCrystalFront.visible = !this.openCaseMode;
  }

  /**
   * 2. Stepped Multi-Tier C36000 Brass / Rhodium-Gilt Mainplate with Perlage Graining,
   *    Recessed 6h Tourbillon Well, Gold Rim Engraving, and Lower Ruby Chatons.
   */
  private buildSteppedMainplateAndJewels(): void {
    const plateShape = new THREE.Shape();
    plateShape.absarc(0, 0, 22.2, 0, Math.PI * 2, false);
    // Open the 6 o'clock Tourbillon exhibition aperture (centered at x=0, y=-11.0, radius 6.6 mm)
    const tourbHole = new THREE.Path();
    tourbHole.absarc(0, -11.0, 6.55, 0, Math.PI * 2, true);
    plateShape.holes.push(tourbHole);

    const plateGeom = new THREE.ExtrudeGeometry(plateShape, {
      depth: 1.65,
      bevelEnabled: true,
      bevelThickness: 0.14,
      bevelSize: 0.14,
      bevelSegments: 3,
      curveSegments: 72,
    });
    plateGeom.translate(0, 0, -1.95);
    const mainplate = new THREE.Mesh(plateGeom, this.matRhodiumPlate);
    mainplate.receiveShadow = true;
    this.layerMainplate.add(mainplate);

    // Chamfered 18k gold inner bezel ring around the 6 o'clock Tourbillon well
    const tourbBezel = new THREE.Mesh(
      new THREE.TorusGeometry(6.65, 0.24, 16, 64),
      this.matSolidGoldStatic,
    );
    tourbBezel.position.set(0, -11.0, -0.32);
    this.layerMainplate.add(tourbBezel);

    // Fixed 70-tooth fourth wheel stationary ring around the Tourbillon well (for the escape pinion to orbit)
    const fixedFourthGeom = buildToothedGearGeometry({
      outerRadius: 5.45,
      innerRimRadius: 4.85,
      hubRadius: 4.65,
      teeth: 60,
      thickness: 0.32,
      spokes: 0,
    });
    const fixedFourthMesh = new THREE.Mesh(fixedFourthGeom, this.matPolishedSteel);
    fixedFourthMesh.position.set(0, -11.0, -0.65);
    this.layerMainplate.add(fixedFourthMesh);

    // Lower Tourbillon Cock at 6h supporting the lower Incabloc ruby endstone
    const lowerCock = new THREE.Mesh(
      new THREE.BoxGeometry(14.5, 2.2, 0.45),
      this.matPolishedSteel,
    );
    lowerCock.position.set(0, -11.0, -1.75);
    this.layerMainplate.add(lowerCock);

    const lowerTourbChaton = createRubyChatonAssembly(
      1.35,
      0.72,
      this.matSolidGoldStatic,
      this.matRubyJewel,
      this.matPolishedSteel,
      this.matScrewSlot,
      true,
    );
    lowerTourbChaton.position.set(0, -11.0, -1.5);
    this.layerMainplate.add(lowerTourbChaton);

    // Mainplate perimeter locator pillars & jewel chatons
    const plateChatonCoords: Array<[number, number]> = [
      [-5.8, 2.2],
      [5.8, 2.2],
      [-9.5, -3.5],
      [9.5, -3.5],
      [-12.2, 7.5],
      [12.2, 7.5],
    ];
    for (const [cx, cy] of plateChatonCoords) {
      const ch = createRubyChatonAssembly(
        1.1,
        0.58,
        this.matSolidGoldStatic,
        this.matRubyJewel,
        this.matPolishedSteel,
        this.matScrewSlot,
      );
      ch.position.set(cx, cy, -0.25);
      this.layerMainplate.add(ch);
    }
  }

  /**
   * 3. Complete Keyless Works at 12 o'clock:
   *    Winding stem, sliding clutch pinion, winding pinion, contoured setting lever (tirette),
   *    rocking yoke (bascule), curved steel yoke spring, 42-tooth sunburst crown wheel,
   *    78-tooth sunburst ratchet wheel, and dynamic recoiling click (cliquet) + click spring.
   */
  private buildKeylessWorksAndRatchetTrain(): void {
    // Fluted steel winding stem entering at 12h (x=0, y=18..28)
    const stem = new THREE.Mesh(
      new THREE.CylinderGeometry(0.68, 0.68, 12.5, 20),
      this.matPolishedSteel,
    );
    stem.position.set(0, 21.8, 0.35);
    this.layerKeyless.add(stem);

    // Sliding clutch pinion (Breguet ratcheted teeth) & Winding pinion on the stem
    const clutchPinion = new THREE.Mesh(
      new THREE.CylinderGeometry(1.45, 1.45, 1.6, 18),
      this.matPolishedSteel,
    );
    clutchPinion.position.set(0, 18.2, 0.35);
    this.layerKeyless.add(clutchPinion);

    const windingPinion = new THREE.Mesh(
      new THREE.CylinderGeometry(1.75, 1.55, 1.4, 20),
      this.matPolishedSteel,
    );
    windingPinion.position.set(0, 20.6, 0.35);
    this.layerKeyless.add(windingPinion);

    // Contoured steel Setting Lever (Tirette) & Rocking Yoke (Bascule) with long curved spring
    const tiretteCurve = new THREE.CatmullRomCurve3([
      new THREE.Vector3(0.8, 19.2, 0.55),
      new THREE.Vector3(3.6, 18.4, 0.55),
      new THREE.Vector3(6.8, 16.2, 0.55),
      new THREE.Vector3(8.5, 13.8, 0.55),
    ]);
    const tirette = new THREE.Mesh(
      new THREE.TubeGeometry(tiretteCurve, 24, 0.42, 12, false),
      this.matPolishedSteel,
    );
    this.layerKeyless.add(tirette);

    const basculeCurve = new THREE.CatmullRomCurve3([
      new THREE.Vector3(-0.8, 18.0, 0.55),
      new THREE.Vector3(-4.5, 17.2, 0.55),
      new THREE.Vector3(-8.2, 14.8, 0.55),
      new THREE.Vector3(-10.8, 11.5, 0.55),
    ]);
    const bascule = new THREE.Mesh(
      new THREE.TubeGeometry(basculeCurve, 24, 0.46, 12, false),
      this.matPolishedSteel,
    );
    this.layerKeyless.add(bascule);

    // Long S-curved blued-steel Yoke Spring
    const yokeSpringCurve = new THREE.CatmullRomCurve3([
      new THREE.Vector3(-13.5, 7.2, 0.58),
      new THREE.Vector3(-12.8, 10.5, 0.58),
      new THREE.Vector3(-10.4, 12.0, 0.58),
    ]);
    const yokeSpring = new THREE.Mesh(
      new THREE.TubeGeometry(yokeSpringCurve, 20, 0.18, 10, false),
      this.matBluedSteel,
    );
    this.layerKeyless.add(yokeSpring);

    // 42-Tooth Sunburst Crown Wheel (x = -5.8, y = 15.2, z = 1.38)
    this.crownWheelGroup.position.set(-5.8, 15.2, 1.38);
    const crownWheelGeom = buildToothedGearGeometry({
      outerRadius: 3.45,
      innerRimRadius: 1.5,
      hubRadius: 1.35,
      teeth: 42,
      thickness: 0.42,
      spokes: 0,
    });
    const crownWheelMesh = new THREE.Mesh(crownWheelGeom, this.matSunburstSteel);
    this.crownWheelGroup.add(crownWheelMesh);
    const cwCore = createSlottedWatchScrew(1.15, 0.48, this.matBluedSteel, this.matScrewSlot);
    this.crownWheelGroup.add(cwCore);
    this.layerKeyless.add(this.crownWheelGroup);

    // 78-Tooth Sunburst Ratchet Wheel mounted atop the Mainspring Barrel Arbor (x = 0, y = 11.5, z = 1.42)
    this.ratchetWheelGroup.position.set(0, 11.5, 1.42);
    const ratchetGeom = buildToothedGearGeometry({
      outerRadius: 5.45,
      innerRimRadius: 1.8,
      hubRadius: 1.6,
      teeth: 64,
      thickness: 0.44,
      spokes: 0,
    });
    const ratchetMesh = new THREE.Mesh(ratchetGeom, this.matSunburstSteel);
    this.ratchetWheelGroup.add(ratchetMesh);
    const arborScrew = createSlottedWatchScrew(1.35, 0.52, this.matBluedSteel, this.matScrewSlot);
    this.ratchetWheelGroup.add(arborScrew);
    this.layerKeyless.add(this.ratchetWheelGroup);

    // Dynamic Recoil Click Lever (Cliquet) at (x = 4.95, y = 14.8, z = 1.44) + Curved Click Spring
    this.clickLeverGroup.position.set(4.95, 14.8, 1.44);
    const clickPivotScrew = createSlottedWatchScrew(
      0.68,
      0.48,
      this.matBluedSteel,
      this.matScrewSlot,
    );
    this.clickLeverGroup.add(clickPivotScrew);

    const clickBeakCurve = new THREE.CatmullRomCurve3([
      new THREE.Vector3(0, 0, 0),
      new THREE.Vector3(-0.65, -1.1, 0),
      new THREE.Vector3(-1.35, -2.05, 0),
    ]);
    const clickBeak = new THREE.Mesh(
      new THREE.TubeGeometry(clickBeakCurve, 16, 0.34, 10, false),
      this.matPolishedSteel,
    );
    this.clickLeverGroup.add(clickBeak);
    this.layerKeyless.add(this.clickLeverGroup);

    const clickSpringCurve = new THREE.CatmullRomCurve3([
      new THREE.Vector3(7.4, 12.2, 1.42),
      new THREE.Vector3(6.8, 14.2, 1.42),
      new THREE.Vector3(5.3, 14.5, 1.42),
    ]);
    const clickSpring = new THREE.Mesh(
      new THREE.TubeGeometry(clickSpringCurve, 20, 0.16, 10, false),
      this.matBluedSteel,
    );
    this.layerKeyless.add(clickSpring);
  }

  /**
   * 4. Cutaway Mainspring Barrel at 12h (`x = 0, y = 11.5`) & 18-Turn 3D Extruded Elastica Ribbon Mainspring.
   */
  private buildMainspringBarrelAndElasticaRibbon(): void {
    this.barrelWheelGroup.position.set(0, 11.5, 0.25);

    // 78-tooth brass barrel drum with 4 skeletonized inspection windows in the top cover so the 18-turn
    // high-carbon steel mainspring coils are completely visible as they wind and unwind!
    const barrelGearGeom = buildToothedGearGeometry({
      outerRadius: 7.45,
      innerRimRadius: 6.45,
      hubRadius: 1.85,
      teeth: 78,
      thickness: 1.35,
      spokes: 4,
    });
    const barrelDrum = new THREE.Mesh(barrelGearGeom, this.matBrassWheel);
    barrelDrum.castShadow = true;
    this.barrelWheelGroup.add(barrelDrum);

    // Central steel barrel arbor & hook
    const arbor = new THREE.Mesh(
      new THREE.CylinderGeometry(1.65, 1.65, 1.75, 24),
      this.matPolishedSteel,
    );
    arbor.rotation.x = Math.PI * 0.5;
    this.barrelWheelGroup.add(arbor);

    // True 18-Turn 3D Rectangular-Ribbon High-Carbon Steel Mainspring inside the barrel
    const samples = 18 * 24;
    this.mainspringRibbonGeometry = buildSpiralRibbonGeometry(samples, 1.05, 0.085);
    const defaultCoils = new Float32Array(18);
    for (let k = 0; k < 18; k++) {
      defaultCoils[k] = 1.75 + (k / 17) * 4.55;
    }
    this.updateMainspringRibbonFromCoils(defaultCoils);

    const mainspringMesh = new THREE.Mesh(this.mainspringRibbonGeometry, this.matBluedSteel);
    mainspringMesh.position.set(0, 11.5, 0.25);
    this.layerBarrel.add(this.barrelWheelGroup, mainspringMesh);
  }

  private updateMainspringRibbonFromCoils(coilsMm: Float32Array): void {
    updateSpiralRibbonVertices(this.mainspringRibbonGeometry, (i, samples) => {
      const u = i / (samples - 1);
      const coilIdxFloat = u * 17.0;
      const k0 = Math.min(17, Math.floor(coilIdxFloat));
      const k1 = Math.min(17, k0 + 1);
      const frac = coilIdxFloat - k0;
      const r = (coilsMm[k0] ?? 2.0) * (1.0 - frac) + (coilsMm[k1] ?? 6.2) * frac;
      const theta = u * 18.0 * Math.PI * 2.0;
      const nx = Math.cos(theta);
      const ny = Math.sin(theta);
      return { x: r * nx, y: r * ny, z: 0.0, nx, ny };
    });
  }

  /**
   * 5. Complete 4-Stage Going Train (Center Wheel 80T, Third Wheel 75T, Fourth/Intermediate Wheels,
   *    Steel Pinions & Motion Works).
   */
  private buildGoingTrainWheels(): void {
    // Center Wheel (80T, x = 0, y = 0, z = 0.15)
    this.centerWheelGroup.position.set(0, 0, 0.15);
    const centerGeom = buildToothedGearGeometry({
      outerRadius: 6.1,
      innerRimRadius: 5.15,
      hubRadius: 1.35,
      teeth: 64,
      thickness: 0.42,
      spokes: 5,
    });
    const centerMesh = new THREE.Mesh(centerGeom, this.matBrassWheel);
    const centerPinion = new THREE.Mesh(
      new THREE.CylinderGeometry(1.15, 1.15, 2.2, 12),
      this.matPolishedSteel,
    );
    centerPinion.rotation.x = Math.PI * 0.5;
    this.centerWheelGroup.add(centerMesh, centerPinion);

    // Third Wheel (75T, x = -5.8, y = 2.2, z = 0.35)
    this.thirdWheelGroup.position.set(-5.8, 2.2, 0.35);
    const thirdGeom = buildToothedGearGeometry({
      outerRadius: 4.55,
      innerRimRadius: 3.75,
      hubRadius: 1.05,
      teeth: 56,
      thickness: 0.38,
      spokes: 5,
    });
    const thirdMesh = new THREE.Mesh(thirdGeom, this.matBrassWheel);
    const thirdPinion = new THREE.Mesh(
      new THREE.CylinderGeometry(0.88, 0.88, 1.8, 10),
      this.matPolishedSteel,
    );
    thirdPinion.rotation.x = Math.PI * 0.5;
    this.thirdWheelGroup.add(thirdMesh, thirdPinion);

    // Symmetrical right-flank date/moonphase intermediate transmission wheel (x = +5.8, y = 2.2, z = 0.25)
    const auxWheelGroup = new THREE.Group();
    auxWheelGroup.position.set(5.8, 2.2, 0.25);
    const auxGeom = buildToothedGearGeometry({
      outerRadius: 4.2,
      innerRimRadius: 3.45,
      hubRadius: 0.95,
      teeth: 48,
      thickness: 0.36,
      spokes: 5,
    });
    auxWheelGroup.add(new THREE.Mesh(auxGeom, this.matBrassWheel));

    this.layerGoingTrain.add(this.centerWheelGroup, this.thirdWheelGroup, auxWheelGroup);
  }

  /**
   * 6. 4 Contoured Geneva Movement Bridges with Côtes de Genève Stripes, Polished Anglage Bevels,
   *    28 Slotted Blued-Steel Screws, Brass Steady Pins, and 3-Screw Gold Ruby Chatons,
   *    plus the Mirror-Polished Rounded-Steel Patek Upper Tourbillon Bridge (Pont Barreau).
   */
  private buildContouredGenevaBridges(): void {
    const createBridgeMesh = (pts: Array<[number, number]>, zPos: number): THREE.Mesh => {
      const s = new THREE.Shape();
      pts.forEach(([x, y], idx) => {
        if (idx === 0) s.moveTo(x, y);
        else s.lineTo(x, y);
      });
      s.closePath();
      const g = new THREE.ExtrudeGeometry(s, {
        depth: 0.58,
        bevelEnabled: true,
        bevelThickness: 0.08,
        bevelSize: 0.08,
        bevelSegments: 3,
      });
      g.translate(0, 0, zPos);
      const m = new THREE.Mesh(g, this.matGenevaBridge);
      m.castShadow = true;
      m.receiveShadow = true;
      return m;
    };

    // Bridge 1: Left Barrel & Winding Bridge
    const bridgeLeft = createBridgeMesh(
      [
        [-18.8, 6.2],
        [-7.5, 7.8],
        [-6.2, 16.5],
        [-12.5, 18.2],
        [-19.2, 13.0],
      ],
      0.72,
    );
    // Bridge 2: Right Barrel & Complication Bridge
    const bridgeRight = createBridgeMesh(
      [
        [18.8, 6.2],
        [7.5, 7.8],
        [6.2, 16.5],
        [12.5, 18.2],
        [19.2, 13.0],
      ],
      0.72,
    );
    // Bridge 3: S-Curved Going-Train Bridge (Left Flank)
    const bridgeTrainLeft = createBridgeMesh(
      [
        [-19.5, -5.2],
        [-8.8, -3.5],
        [-3.8, 1.2],
        [-7.8, 4.8],
        [-18.5, 3.2],
      ],
      0.78,
    );
    // Bridge 4: Center & Right Flank Bridge
    const bridgeCenterRight = createBridgeMesh(
      [
        [-1.8, -1.6],
        [8.8, -3.5],
        [19.5, -5.2],
        [18.5, 3.2],
        [7.8, 4.8],
        [1.8, 1.6],
      ],
      0.78,
    );

    this.layerBridges.add(bridgeLeft, bridgeRight, bridgeTrainLeft, bridgeCenterRight);

    // 24 Slotted Blued-Steel Bridge Screws across the 4 Geneva bridges
    const screwCoords: Array<[number, number, number]> = [
      [-16.8, 8.2, 1.38],
      [-10.8, 16.2, 1.38],
      [-16.5, 12.2, 1.38],
      [16.8, 8.2, 1.38],
      [10.8, 16.2, 1.38],
      [16.5, 12.2, 1.38],
      [-16.8, -3.2, 1.44],
      [-14.2, 1.8, 1.44],
      [-9.5, -1.8, 1.44],
      [16.8, -3.2, 1.44],
      [14.2, 1.8, 1.44],
      [9.5, -1.8, 1.44],
    ];
    for (const [sx, sy, sz] of screwCoords) {
      const sc = createSlottedWatchScrew(0.76, 0.32, this.matBluedSteel, this.matScrewSlot);
      sc.position.set(sx, sy, sz);
      sc.rotation.z = sx * 0.4 + sy * 0.3;
      this.layerBridges.add(sc);
    }

    // Raised 3-Screw 18k Gold Ruby Chatons on the Bridges (Center wheel, Third wheel, Intermediate wheels)
    const bridgeChatons: Array<[number, number, number]> = [
      [0.0, 0.0, 1.42],
      [-5.8, 2.2, 1.42],
      [5.8, 2.2, 1.42],
      [-11.5, 11.2, 1.36],
      [11.5, 11.2, 1.36],
    ];
    for (const [cx, cy, cz] of bridgeChatons) {
      const ch = createRubyChatonAssembly(
        1.25,
        0.66,
        this.matSolidGoldStatic,
        this.matRubyJewel,
        this.matPolishedSteel,
        this.matScrewSlot,
      );
      ch.position.set(cx, cy, cz);
      this.layerBridges.add(ch);
    }

    // Iconic Patek Mirror-Polished Rounded-Steel Upper Tourbillon Bridge (Pont Barreau) at 6 o'clock
    const tourbBridgeGroup = new THREE.Group();
    tourbBridgeGroup.position.set(0, -11.0, 2.15);

    const barArmLeft = new THREE.Mesh(
      new THREE.CylinderGeometry(0.32, 0.48, 7.6, 20),
      this.matPolishedSteel,
    );
    barArmLeft.rotation.z = Math.PI * 0.5;
    barArmLeft.position.set(-4.6, 0, 0);

    const barArmRight = new THREE.Mesh(
      new THREE.CylinderGeometry(0.32, 0.48, 7.6, 20),
      this.matPolishedSteel,
    );
    barArmRight.rotation.z = -Math.PI * 0.5;
    barArmRight.position.set(4.6, 0, 0);

    const footLeft = new THREE.Mesh(
      new THREE.CylinderGeometry(1.55, 1.65, 0.42, 24),
      this.matPolishedSteel,
    );
    footLeft.rotation.x = Math.PI * 0.5;
    footLeft.position.set(-8.2, 0, -0.1);

    const footRight = new THREE.Mesh(
      new THREE.CylinderGeometry(1.55, 1.65, 0.42, 24),
      this.matPolishedSteel,
    );
    footRight.rotation.x = Math.PI * 0.5;
    footRight.position.set(8.2, 0, -0.1);

    // 4 blued-steel retaining screws + 2 gold locator pins on the Tourbillon bridge feet
    for (const side of [-1, 1]) {
      for (const dy of [-0.75, 0.75]) {
        const sc = createSlottedWatchScrew(0.52, 0.34, this.matBluedSteel, this.matScrewSlot);
        sc.position.set(side * 8.2, dy, 0.15);
        tourbBridgeGroup.add(sc);
      }
    }

    // Center upper Incabloc Ruby Endstone Chaton on the Tourbillon bridge
    const upperTourbChaton = createRubyChatonAssembly(
      1.35,
      0.72,
      this.matSolidGoldStatic,
      this.matRubyJewel,
      this.matPolishedSteel,
      this.matScrewSlot,
      true,
    );
    upperTourbChaton.position.set(0, 0, 0.12);

    tourbBridgeGroup.add(barArmLeft, barArmRight, footLeft, footRight, upperTourbChaton);
    this.layerBridges.add(tourbBridgeGroup);
  }

  /**
   * 7. 48-Part 1-Minute Titanium Tourbillon Carriage at 6h (`x = 0, y = -11.0`),
   *    15-Tooth Club-Tooth Steel Escape Wheel, Counterpoised Swiss Lever Pallet Fork
   *    with Angled Synthetic Ruby Pallet Stones, Double-Roller Glucydur Balance Wheel
   *    with 16 Gold Timing Screws + 4 Microstella Nuts, and 14-Turn 3D Breguet Overcoil Hairspring.
   */
  private buildTourbillonEscapementAndBalance(): void {
    this.tourbillonCageGroup.position.set(0, -11.0, 0.45);

    // Lower & Upper 3-Arm Lyric Titanium Carriage Rings + 3 Waisted Pillars with Gold Nuts
    const lowerCageRing = new THREE.Mesh(
      new THREE.TorusGeometry(5.65, 0.22, 16, 60),
      this.matTitaniumCage,
    );
    lowerCageRing.position.z = -1.1;

    const upperCageRing = new THREE.Mesh(
      new THREE.TorusGeometry(5.65, 0.24, 16, 60),
      this.matTitaniumCage,
    );
    upperCageRing.position.z = 1.15;
    this.tourbillonCageGroup.add(lowerCageRing, upperCageRing);

    for (let k = 0; k < 3; k++) {
      const a = (k * Math.PI * 2) / 3 + Math.PI / 6;
      const px = Math.cos(a) * 5.65;
      const py = Math.sin(a) * 5.65;

      // Turned waisted titanium cage pillar
      const pillar = new THREE.Mesh(
        new THREE.CylinderGeometry(0.28, 0.28, 2.25, 16),
        this.matTitaniumCage,
      );
      pillar.rotation.x = Math.PI * 0.5;
      pillar.position.set(px, py, 0.02);

      // Gold hexagonal nut on top of pillar
      const nut = new THREE.Mesh(
        new THREE.CylinderGeometry(0.45, 0.45, 0.22, 6),
        this.matSolidGoldStatic,
      );
      nut.rotation.x = Math.PI * 0.5;
      nut.position.set(px, py, 1.28);

      // Lower & upper lyric spokes
      const spokeLow = new THREE.Mesh(
        new THREE.BoxGeometry(5.65, 0.32, 0.18),
        this.matTitaniumCage,
      );
      spokeLow.position.set(px * 0.5, py * 0.5, -1.1);
      spokeLow.rotation.z = a;

      const spokeHigh = new THREE.Mesh(
        new THREE.BoxGeometry(5.65, 0.34, 0.20),
        this.matTitaniumCage,
      );
      spokeHigh.position.set(px * 0.5, py * 0.5, 1.15);
      spokeHigh.rotation.z = a;

      this.tourbillonCageGroup.add(pillar, nut, spokeLow, spokeHigh);
    }

    // 15-Tooth Club-Tooth Steel Escape Wheel inside the Tourbillon Carriage (offset at x = -2.35, y = -0.55)
    this.escapeWheelGroup.position.set(-2.35, -0.55, -0.45);
    const escapeGeom = buildToothedGearGeometry({
      outerRadius: 2.15,
      innerRimRadius: 1.55,
      hubRadius: 0.48,
      teeth: 15,
      thickness: 0.22,
      spokes: 4,
      clubTeeth: true,
    });
    const escapeMesh = new THREE.Mesh(escapeGeom, this.matPolishedSteel);
    const escapePinion = new THREE.Mesh(
      new THREE.CylinderGeometry(0.42, 0.42, 1.45, 8),
      this.matPolishedSteel,
    );
    escapePinion.rotation.x = Math.PI * 0.5;
    const escapeJewel = createRubyChatonAssembly(
      0.65,
      0.35,
      this.matSolidGoldStatic,
      this.matRubyJewel,
      this.matPolishedSteel,
      this.matScrewSlot,
    );
    escapeJewel.position.z = 0.55;
    this.escapeWheelGroup.add(escapeMesh, escapePinion, escapeJewel);
    this.tourbillonCageGroup.add(this.escapeWheelGroup);

    // Counterpoised Steel Swiss Lever Pallet Fork (pivot at x = -1.15, y = 1.05, z = -0.38)
    this.palletForkGroup.position.set(-1.15, 1.05, -0.38);
    // Lever arm toward balance roller (0,0) + pallet yoke arm + elliptical counterpoise
    const leverBeam = new THREE.Mesh(
      new THREE.BoxGeometry(2.45, 0.26, 0.18),
      this.matPolishedSteel,
    );
    leverBeam.position.set(0.35, -0.35, 0);
    leverBeam.rotation.z = -0.68;

    // Fork Horn with Crescent Slot & Guard Dart Pin facing the balance roller at (0,0)
    const hornLeft = new THREE.Mesh(
      new THREE.BoxGeometry(0.55, 0.14, 0.20),
      this.matPolishedSteel,
    );
    hornLeft.position.set(1.05, -0.72, 0);
    hornLeft.rotation.z = -0.45;
    const hornRight = new THREE.Mesh(
      new THREE.BoxGeometry(0.55, 0.14, 0.20),
      this.matPolishedSteel,
    );
    hornRight.position.set(0.82, -0.95, 0);
    hornRight.rotation.z = -0.92;
    const guardDart = new THREE.Mesh(
      new THREE.CylinderGeometry(0.05, 0.02, 0.45, 8),
      this.matSolidGoldStatic,
    );
    guardDart.position.set(0.95, -0.85, 0.08);
    guardDart.rotation.z = 0.85;

    // Pallet Yoke Arms + Translucent Angled Synthetic Ruby Entry & Exit Pallet Stones!
    const palletYoke = new THREE.Mesh(
      new THREE.BoxGeometry(1.65, 0.28, 0.18),
      this.matPolishedSteel,
    );
    palletYoke.position.set(-0.68, -0.55, 0);
    palletYoke.rotation.z = 0.52;

    const entryRubyPallet = new THREE.Mesh(
      new THREE.BoxGeometry(0.24, 0.72, 0.26),
      this.matRubyJewel,
    );
    entryRubyPallet.position.set(-0.35, -0.98, 0);
    entryRubyPallet.rotation.z = 0.32;

    const exitRubyPallet = new THREE.Mesh(
      new THREE.BoxGeometry(0.24, 0.72, 0.26),
      this.matRubyJewel,
    );
    exitRubyPallet.position.set(-1.32, -0.45, 0);
    exitRubyPallet.rotation.z = -0.42;

    const forkPivotJewel = new THREE.Mesh(
      new THREE.CylinderGeometry(0.28, 0.28, 0.36, 16),
      this.matRubyJewel,
    );
    forkPivotJewel.rotation.x = Math.PI * 0.5;

    this.palletForkGroup.add(
      leverBeam,
      hornLeft,
      hornRight,
      guardDart,
      palletYoke,
      entryRubyPallet,
      exitRubyPallet,
      forkPivotJewel,
    );
    this.tourbillonCageGroup.add(this.palletForkGroup);

    // Twin Eccentric Brass Banking Pins limiting pallet fork travel
    for (const [bx, by] of [
      [-0.25, 0.55],
      [-0.65, 0.15],
    ] as const) {
      const bankPin = new THREE.Mesh(
        new THREE.CylinderGeometry(0.11, 0.11, 0.65, 12),
        this.matSolidGoldStatic,
      );
      bankPin.rotation.x = Math.PI * 0.5;
      bankPin.position.set(bx, by, -0.38);
      this.tourbillonCageGroup.add(bankPin);
    }

    // Pierced 4-Arm Glucydur Beryllium-Copper Balance Wheel + Double Roller + 16 Gold Screws + 4 Microstella Nuts
    this.balanceWheelGroup.position.set(0, 0, 0.22);

    const balanceRim = new THREE.Mesh(
      new THREE.TorusGeometry(4.55, 0.34, 20, 64),
      this.matGlucydurBalance,
    );
    this.balanceWheelGroup.add(balanceRim);

    // 4 Contoured Balance Arms
    for (let k = 0; k < 4; k++) {
      const a = (k * Math.PI) / 2;
      const arm = new THREE.Mesh(
        new THREE.BoxGeometry(4.5, 0.32, 0.20),
        this.matGlucydurBalance,
      );
      arm.position.set(Math.cos(a) * 2.25, Math.sin(a) * 2.25, 0);
      arm.rotation.z = a;
      this.balanceWheelGroup.add(arm);

      // 4 Slotted Gold Microstella Regulating Nuts facing inward on the rim
      const microstella = new THREE.Mesh(
        new THREE.CylinderGeometry(0.32, 0.32, 0.28, 6),
        this.matSolidGoldStatic,
      );
      microstella.position.set(Math.cos(a + 0.35) * 4.12, Math.sin(a + 0.35) * 4.12, 0);
      microstella.rotation.z = a + 0.35 + Math.PI * 0.5;
      this.balanceWheelGroup.add(microstella);
    }

    // 16 Slotted 18k Gold Poising Screws around the Glucydur Balance Rim
    for (let k = 0; k < 16; k++) {
      const a = (k * Math.PI * 2) / 16;
      const screw = new THREE.Mesh(
        new THREE.CylinderGeometry(0.22, 0.16, 0.55, 12),
        this.matSolidGoldStatic,
      );
      screw.position.set(Math.cos(a) * 4.95, Math.sin(a) * 4.95, 0);
      screw.rotation.z = a - Math.PI * 0.5;
      this.balanceWheelGroup.add(screw);
    }

    // Double-Roller Assembly on the Balance Staff:
    // 1. Large steel impulse roller table + D-shaped vertical synthetic ruby impulse pin (ellipse)
    const largeRoller = new THREE.Mesh(
      new THREE.CylinderGeometry(0.95, 0.95, 0.18, 28),
      this.matPolishedSteel,
    );
    largeRoller.rotation.x = Math.PI * 0.5;
    largeRoller.position.z = -0.42;

    const rubyImpulsePin = new THREE.Mesh(
      new THREE.CylinderGeometry(0.12, 0.12, 0.45, 14),
      this.matRubyJewel,
    );
    rubyImpulsePin.rotation.x = Math.PI * 0.5;
    rubyImpulsePin.position.set(-0.68, 0.42, -0.55);

    // 2. Small safety roller with crescent passing notch
    const safetyRoller = new THREE.Mesh(
      new THREE.CylinderGeometry(0.56, 0.56, 0.22, 24),
      this.matPolishedSteel,
    );
    safetyRoller.rotation.x = Math.PI * 0.5;
    safetyRoller.position.z = -0.62;

    const balanceStaff = new THREE.Mesh(
      new THREE.CylinderGeometry(0.24, 0.24, 2.35, 16),
      this.matPolishedSteel,
    );
    balanceStaff.rotation.x = Math.PI * 0.5;

    this.balanceWheelGroup.add(largeRoller, rubyImpulsePin, safetyRoller, balanceStaff);
    this.tourbillonCageGroup.add(this.balanceWheelGroup);

    // True 14-Turn 3D Rectangular-Ribbon Breguet Overcoil Hairspring (with Terminal Overcoil Riser & Stud!)
    const hsSamples = 14 * 28;
    this.hairspringRibbonGeometry = buildSpiralRibbonGeometry(hsSamples, 0.24, 0.045);
    this.updateHairspringRibbon(0.0, 0.0);
    const hairspringMesh = new THREE.Mesh(this.hairspringRibbonGeometry, this.matHairspring);
    hairspringMesh.position.z = 0.65;
    this.tourbillonCageGroup.add(hairspringMesh);

    // Split Steel Collet at center + Triangular Hairspring Stud & Curb Regulator Pins at outer overcoil
    const collet = new THREE.Mesh(
      new THREE.CylinderGeometry(0.52, 0.52, 0.28, 20),
      this.matPolishedSteel,
    );
    collet.rotation.x = Math.PI * 0.5;
    collet.position.z = 0.65;
    this.tourbillonCageGroup.add(collet);

    this.layerTourbillon.add(this.tourbillonCageGroup);
  }

  private updateHairspringRibbon(balanceAngleRad: number, oilBridgeFraction: number): void {
    const breathAmp = Math.sin(balanceAngleRad * 0.5) * 0.16;
    updateSpiralRibbonVertices(this.hairspringRibbonGeometry, (i, samples) => {
      const u = i / (samples - 1);
      // First 90% of samples: 13 planar Archimedes turns from collet (r=0.58 mm) to outer coil (r=3.45 mm)
      // Final 10% of samples: Breguet terminal overcoil rising +0.34 mm in Z and curving inward to r=2.65 mm!
      if (u <= 0.90) {
        const v = u / 0.90;
        const theta = v * 13.0 * Math.PI * 2.0 + balanceAngleRad * (1.0 - v);
        // Breathing deformation + capillary oil meniscus coil-bridging pinch
        const pinch =
          oilBridgeFraction > 0.01 && v > 0.62 ? Math.sin(v * 38.0) * 0.14 * oilBridgeFraction : 0.0;
        const r = 0.58 + v * 2.85 + breathAmp * Math.sin(v * Math.PI) * Math.cos(theta) + pinch;
        return {
          x: r * Math.cos(theta),
          y: r * Math.sin(theta),
          z: 0.0,
          nx: Math.cos(theta),
          ny: Math.sin(theta),
        };
      } else {
        const w = (u - 0.90) / 0.10;
        const smoothW = w * w * (3.0 - 2.0 * w);
        const thetaStart = 13.0 * Math.PI * 2.0;
        const theta = thetaStart + w * Math.PI * 1.35;
        const r = 3.43 * (1.0 - 0.24 * smoothW);
        const z = 0.34 * Math.sin(Math.min(1.0, w * 1.8) * Math.PI * 0.5);
        return {
          x: r * Math.cos(theta),
          y: r * Math.sin(theta),
          z,
          nx: Math.cos(theta),
          ny: Math.sin(theta),
        };
      }
    });
  }

  /**
   * 8. Skeletonized Exhibition Dial with Clous de Paris Guilloché, Railway Minute Track,
   *    11 Applied 3D Gold Hour Markers, 10h Power-Reserve Sector Subdial, 2h Lapis-Lazuli
   *    & 18k Gold Moonphase Complication, and Sculpted 3D Blued-Steel Breguet Pomme Hands.
   */
  private buildExhibitionGuillocheDialAndHands(guillocheNormal: THREE.CanvasTexture): void {
    const dialMat = new THREE.MeshStandardMaterial({
      color: new THREE.Color(0xf4efe4),
      metalness: 0.25,
      roughness: 0.32,
      normalMap: guillocheNormal,
      normalScale: new THREE.Vector2(0.35, 0.35),
    });

    // Peripheral skeletonized chapter ring (R = 18.6..22.1 mm) leaving the entire movement & bridges visible!
    const ringShape = new THREE.Shape();
    ringShape.absarc(0, 0, 22.1, 0, Math.PI * 2, false);
    const innerCut = new THREE.Path();
    innerCut.absarc(0, 0, 18.6, 0, Math.PI * 2, true);
    ringShape.holes.push(innerCut);

    const chapterRing = new THREE.Mesh(
      new THREE.ExtrudeGeometry(ringShape, { depth: 0.26, bevelEnabled: false, curveSegments: 72 }),
      dialMat,
    );
    chapterRing.position.z = 1.85;
    this.layerDialHands.add(chapterRing);

    // 60 Railway Minute Pearls + 11 Applied 3D 18k Gold Hour Baton/Cabochon Markers (skipping 6h Tourbillon)
    for (let m = 0; m < 60; m++) {
      const a = Math.PI * 0.5 - (m * Math.PI * 2) / 60;
      const isHour = m % 5 === 0;
      const pearl = new THREE.Mesh(
        new THREE.SphereGeometry(isHour ? 0.24 : 0.11, 10, 10),
        this.matSolidGoldStatic,
      );
      pearl.position.set(Math.cos(a) * 21.4, Math.sin(a) * 21.4, 2.14);
      this.layerDialHands.add(pearl);
    }

    for (let h = 1; h <= 12; h++) {
      if (h === 6) continue; // 6 o'clock is the open Tourbillon bridge
      const a = Math.PI * 0.5 - (h * Math.PI * 2) / 12;
      const marker = new THREE.Mesh(
        new THREE.BoxGeometry(0.46, 1.85, 0.24),
        this.matSolidGoldStatic,
      );
      marker.position.set(Math.cos(a) * 20.1, Math.sin(a) * 20.1, 2.18);
      marker.rotation.z = a - Math.PI * 0.5;
      this.layerDialHands.add(marker);
    }

    // 10 o'clock Power-Reserve Sector Subdial (x = -8.6, y = +5.2, z = 1.92)
    const prRingShape = new THREE.Shape();
    prRingShape.absarc(0, 0, 4.1, 0, Math.PI * 2, false);
    const prHole = new THREE.Path();
    prHole.absarc(0, 0, 3.1, 0, Math.PI * 2, true);
    prRingShape.holes.push(prHole);
    const prSubdial = new THREE.Mesh(
      new THREE.ExtrudeGeometry(prRingShape, { depth: 0.18, bevelEnabled: false }),
      dialMat,
    );
    prSubdial.position.set(-8.6, 5.2, 1.92);
    this.layerDialHands.add(prSubdial);

    this.powerReserveHandGroup.position.set(-8.6, 5.2, 2.18);
    const prHand = new THREE.Mesh(
      new THREE.ConeGeometry(0.24, 3.4, 12),
      this.matBluedSteel,
    );
    prHand.position.y = 1.3;
    this.powerReserveHandGroup.add(prHand);
    this.layerDialHands.add(this.powerReserveHandGroup);

    // 2 o'clock Astronomical Moonphase Complication (x = +8.6, y = +5.2, z = 1.90)
    const moonGroup = new THREE.Group();
    moonGroup.position.set(8.6, 5.2, 1.90);
    const nightSkyMat = new THREE.MeshStandardMaterial({
      color: new THREE.Color(0x0c1b48),
      metalness: 0.35,
      roughness: 0.25,
    });
    const skyDisc = new THREE.Mesh(
      new THREE.CylinderGeometry(3.85, 3.85, 0.16, 40),
      nightSkyMat,
    );
    skyDisc.rotation.x = Math.PI * 0.5;
    this.moonDiscGroup.add(skyDisc);

    for (const sign of [-1, 1]) {
      const moonOrb = new THREE.Mesh(
        new THREE.CylinderGeometry(1.25, 1.25, 0.20, 28),
        this.matSolidGoldStatic,
      );
      moonOrb.rotation.x = Math.PI * 0.5;
      moonOrb.position.set(sign * 1.95, 0, 0.03);
      this.moonDiscGroup.add(moonOrb);
    }
    // 16 gold stars on the Lapis-Lazuli disc
    for (let s = 0; s < 16; s++) {
      const a = (s * Math.PI * 2) / 16 + 0.2;
      const r = 1.1 + (s % 3) * 0.85;
      const star = new THREE.Mesh(new THREE.OctahedronGeometry(0.14, 0), this.matSolidGoldStatic);
      star.position.set(Math.cos(a) * r, Math.sin(a) * r, 0.09);
      this.moonDiscGroup.add(star);
    }
    const moonBezel = new THREE.Mesh(
      new THREE.TorusGeometry(3.92, 0.18, 16, 48),
      this.matSolidGoldStatic,
    );
    moonBezel.position.z = 0.08;
    moonGroup.add(this.moonDiscGroup, moonBezel);
    this.layerDialHands.add(moonGroup);

    // Sculpted 3D Blued-Steel Breguet Pomme (Hollow Moon-Tip) Hour & Minute Hands
    const createBreguetPommeHand = (lengthMm: number, pommeRadiusMm: number): THREE.Group => {
      const g = new THREE.Group();
      const shaftLen = lengthMm * 0.68;
      const shaft = new THREE.Mesh(
        new THREE.CylinderGeometry(0.14, 0.34, shaftLen, 16),
        this.matBluedSteel,
      );
      shaft.position.y = shaftLen * 0.5;
      // Hollow eccentric Breguet moon ring near the tip
      const pommeRing = new THREE.Mesh(
        new THREE.TorusGeometry(pommeRadiusMm, 0.16, 14, 32),
        this.matBluedSteel,
      );
      pommeRing.position.y = shaftLen + pommeRadiusMm * 0.85;
      // Needle pointer tip beyond the pomme ring
      const tipLen = lengthMm - (shaftLen + pommeRadiusMm * 1.6);
      const tip = new THREE.Mesh(
        new THREE.ConeGeometry(0.14, Math.max(1.2, tipLen), 14),
        this.matBluedSteel,
      );
      tip.position.y = shaftLen + pommeRadiusMm * 1.65 + Math.max(0.6, tipLen * 0.5);
      // Counterpoise tail
      const tail = new THREE.Mesh(
        new THREE.CylinderGeometry(0.26, 0.10, 2.6, 12),
        this.matBluedSteel,
      );
      tail.position.y = -1.3;
      const boss = new THREE.Mesh(
        new THREE.CylinderGeometry(0.88, 0.88, 0.32, 24),
        this.matSolidGoldStatic,
      );
      boss.rotation.x = Math.PI * 0.5;
      g.add(shaft, pommeRing, tip, tail, boss);
      return g;
    };

    this.hourHandGroup.position.set(0, 0, 2.32);
    this.hourHandGroup.add(createBreguetPommeHand(13.8, 1.05));

    this.minuteHandGroup.position.set(0, 0, 2.58);
    this.minuteHandGroup.add(createBreguetPommeHand(19.2, 0.92));

    // Small-seconds blued-steel hand mounted atop the 6h Tourbillon Carriage axis
    this.secondsHandGroup.position.set(0, -11.0, 2.48);
    const secNeedle = new THREE.Mesh(
      new THREE.CylinderGeometry(0.05, 0.14, 5.8, 10),
      this.matBluedSteel,
    );
    secNeedle.position.y = 2.1;
    this.secondsHandGroup.add(secNeedle);

    this.layerDialHands.add(this.hourHandGroup, this.minuteHandGroup, this.secondsHandGroup);
  }

  /**
   * 9. Continuous 3D Free-Surface Molten & Solidified Gold Heightfield Meshes reconstructed via
   *    C^2 Bicubic B-Spline Interpolation (`84 x 84`) over the `48 x 48` Spatial Finite-Volume
   *    PDE Grid from WASM (`sim.spatial_grid_snapshot()`), with 4-channel RGBA per-pixel iso-contour
   *    meniscus feathering so liquid and frozen gold pools never show Cartesian quad edges.
   */
  private buildContinuumMoltenAndSolidGoldHeightfields(): void {
    const buildGridGeom = (): THREE.BufferGeometry => {
      const vCount = RENDER_GRID_N * RENDER_GRID_N;
      const pos = new Float32Array(vCount * 3);
      const col = new Float32Array(vCount * 4);
      const idx: number[] = [];
      for (let gy = 0; gy < RENDER_GRID_N; gy++) {
        for (let gx = 0; gx < RENDER_GRID_N; gx++) {
          const i = gy * RENDER_GRID_N + gx;
          pos[i * 3] = -GRID_HALF_SPAN_MM + (gx + 0.5) * RENDER_GRID_DX_MM;
          pos[i * 3 + 1] = -GRID_HALF_SPAN_MM + (gy + 0.5) * RENDER_GRID_DX_MM;
          pos[i * 3 + 2] = -10.0;
          col[i * 4] = 0.96;
          col[i * 4 + 1] = 0.78;
          col[i * 4 + 2] = 0.34;
          col[i * 4 + 3] = 0.0;
        }
      }
      for (let gy = 0; gy < RENDER_GRID_N - 1; gy++) {
        for (let gx = 0; gx < RENDER_GRID_N - 1; gx++) {
          const i00 = gy * RENDER_GRID_N + gx;
          const i10 = i00 + 1;
          const i01 = i00 + RENDER_GRID_N;
          const i11 = i01 + 1;
          idx.push(i00, i10, i11, i00, i11, i01);
        }
      }
      const g = new THREE.BufferGeometry();
      g.setAttribute("position", new THREE.BufferAttribute(pos, 3));
      g.setAttribute("color", new THREE.BufferAttribute(col, 4));
      g.setIndex(idx);
      return g;
    };

    this.moltenGoldGeometry = buildGridGeom();
    this.moltenGoldMesh = new THREE.Mesh(this.moltenGoldGeometry, this.matMoltenGold);
    this.moltenGoldMesh.visible = false;

    this.solidGoldGeometry = buildGridGeom();
    this.solidGoldMesh = new THREE.Mesh(this.solidGoldGeometry, this.matSolidifiedGold);
    this.solidGoldMesh.visible = false;

    this.watchRoot.add(this.moltenGoldMesh, this.solidGoldMesh);
  }

  /**
   * Build the 3D interactive physics perturbation tools:
   * - Oxy-Butane Blowtorch nozzle & multi-layered soft additive plasma flame plume
   * - Cryogenic LN2 Dewar nozzle & soft textured Leidenfrost vapor cloud
   * - Horological Oil Precision Oiler & stream
   * - Sintered NdFeB N42 Horseshoe/Bar Super-Magnet
   * - Studio Condenser Measurement Microphone
   */
  private buildPerturbationApparatus(softParticleTex: THREE.CanvasTexture): {
    torchFlameCore: THREE.Mesh;
    torchFlamePlume: THREE.Mesh;
    torchPointLight: THREE.PointLight;
    ln2StreamMesh: THREE.Mesh;
    ln2FrostCloud: THREE.Points;
    oilStreamMesh: THREE.Mesh;
    airStreamMesh: THREE.Mesh;
  } {
    // 1. Blowtorch with multi-layered soft additive plasma cones & surface hotspot halo
    const nozzle = new THREE.Mesh(
      new THREE.CylinderGeometry(1.1, 1.65, 11.5, 28),
      this.matBrassWheel,
    );
    nozzle.position.set(0, -14.2, 8.2);
    nozzle.rotation.x = -0.55;
    const nozzleTip = new THREE.Mesh(
      new THREE.CylinderGeometry(0.88, 1.12, 2.2, 28),
      this.matPolishedSteel,
    );
    nozzleTip.position.set(0, -9.2, 5.1);
    nozzleTip.rotation.x = -0.55;
    this.torchGroup.add(nozzle, nozzleTip);

    const torchFlameCore = new THREE.Mesh(
      new THREE.ConeGeometry(1.15, 9.2, 36, 1, true),
      new THREE.MeshBasicMaterial({
        color: 0xaef4ff,
        transparent: true,
        opacity: 0.86,
        blending: THREE.AdditiveBlending,
        depthWrite: false,
        side: THREE.DoubleSide,
      }),
    );
    torchFlameCore.position.set(0, -5.2, 3.5);
    torchFlameCore.rotation.x = 0.98;

    const innerCoreSheath = new THREE.Mesh(
      new THREE.ConeGeometry(1.75, 11.6, 36, 1, true),
      new THREE.MeshBasicMaterial({
        color: 0x3ec8ff,
        transparent: true,
        opacity: 0.48,
        blending: THREE.AdditiveBlending,
        depthWrite: false,
        side: THREE.DoubleSide,
      }),
    );
    torchFlameCore.add(innerCoreSheath);

    const torchFlamePlume = new THREE.Mesh(
      new THREE.ConeGeometry(2.85, 14.5, 36, 1, true),
      new THREE.MeshBasicMaterial({
        color: 0xff8a24,
        transparent: true,
        opacity: 0.42,
        blending: THREE.AdditiveBlending,
        depthWrite: false,
        side: THREE.DoubleSide,
      }),
    );
    torchFlamePlume.position.set(0, -3.2, 2.2);
    torchFlamePlume.rotation.x = 0.98;

    // Soft 3D spherical plasma impingement tip inside torchFlameCore (never slices watch planes!)
    const plasmaTip = new THREE.Mesh(
      new THREE.SphereGeometry(1.35, 24, 16),
      new THREE.MeshBasicMaterial({
        color: 0xffb454,
        transparent: true,
        opacity: 0.44,
        blending: THREE.AdditiveBlending,
        depthWrite: false,
      }),
    );
    plasmaTip.position.set(0, -4.4, 0);
    torchFlameCore.add(plasmaTip);

    const torchPointLight = new THREE.PointLight(0xff9638, 0.0, 48, 1.25);
    torchPointLight.position.set(0, -2.0, 7.5);
    this.torchGroup.add(torchFlameCore, torchFlamePlume, torchPointLight);
    this.torchGroup.visible = false;
    this.scene.add(this.torchGroup);

    // 2. Liquid Nitrogen (LN2) Cryogenic Nozzle & Soft Billowing Vapor Plume
    const ln2Spout = new THREE.Mesh(
      new THREE.CylinderGeometry(1.3, 1.75, 10, 28),
      this.matPolishedSteel,
    );
    ln2Spout.position.set(0, 0, 14);
    ln2Spout.rotation.x = Math.PI * 0.5;
    const ln2StreamMesh = new THREE.Mesh(
      new THREE.CylinderGeometry(0.55, 2.2, 12, 28, 1, true),
      new THREE.MeshBasicMaterial({
        color: 0xc2f4ff,
        transparent: true,
        opacity: 0.52,
        blending: THREE.AdditiveBlending,
        depthWrite: false,
        side: THREE.DoubleSide,
      }),
    );
    ln2StreamMesh.position.set(0, 0, 5.5);
    ln2StreamMesh.rotation.x = Math.PI * 0.5;

    const cloudPts = new Float32Array(140 * 3);
    for (let i = 0; i < 140; i++) {
      const a = (i * 2.3999632) % (Math.PI * 2);
      const r = Math.sqrt(i / 140) * 11.5;
      cloudPts[i * 3] = Math.cos(a) * r;
      cloudPts[i * 3 + 1] = Math.sin(a) * r;
      cloudPts[i * 3 + 2] = ((i * 7) % 19) * 0.18;
    }
    const cloudGeom = new THREE.BufferGeometry();
    cloudGeom.setAttribute("position", new THREE.BufferAttribute(cloudPts, 3));
    const ln2FrostCloud = new THREE.Points(
      cloudGeom,
      new THREE.PointsMaterial({
        map: softParticleTex,
        color: 0xdff8ff,
        size: 3.4,
        transparent: true,
        opacity: 0.46,
        alphaTest: 0.01,
        depthWrite: false,
        blending: THREE.AdditiveBlending,
      }),
    );
    this.ln2Group.add(ln2Spout, ln2StreamMesh, ln2FrostCloud);
    this.ln2Group.visible = false;
    this.scene.add(this.ln2Group);

    // 3. Watch Oil Dropper
    const oilNeedle = new THREE.Mesh(
      new THREE.CylinderGeometry(0.25, 0.85, 11, 20),
      this.matSolidGoldStatic,
    );
    oilNeedle.position.set(0, 0, 12);
    oilNeedle.rotation.x = Math.PI * 0.5;
    const oilStreamMesh = new THREE.Mesh(
      new THREE.CylinderGeometry(0.22, 0.52, 8.5, 20),
      new THREE.MeshPhysicalMaterial({
        color: 0xffb82e,
        transmission: 0.85,
        roughness: 0.05,
        transparent: true,
        opacity: 0.82,
      }),
    );
    oilStreamMesh.position.set(0, 0, 4.2);
    oilStreamMesh.rotation.x = Math.PI * 0.5;
    this.oilGroup.add(oilNeedle, oilStreamMesh);
    this.oilGroup.visible = false;
    this.scene.add(this.oilGroup);

    // 4. Sintered NdFeB N42 Magnet
    const magNorth = new THREE.Mesh(
      new THREE.BoxGeometry(5.5, 4.2, 4.2),
      new THREE.MeshStandardMaterial({ color: 0xd92638, metalness: 0.5, roughness: 0.3 }),
    );
    magNorth.position.x = -2.75;
    const magSouth = new THREE.Mesh(
      new THREE.BoxGeometry(5.5, 4.2, 4.2),
      new THREE.MeshStandardMaterial({ color: 0x2462d1, metalness: 0.5, roughness: 0.3 }),
    );
    magSouth.position.x = 2.75;
    this.magnetGroup.add(magNorth, magSouth);
    this.magnetGroup.visible = false;
    this.scene.add(this.magnetGroup);

    // 4b. High-Velocity Compressed Air / Argon Shear Nozzle
    const airNozzleBarrel = new THREE.Mesh(
      new THREE.CylinderGeometry(0.95, 1.5, 10.5, 24),
      this.matPolishedSteel,
    );
    airNozzleBarrel.position.set(0, -10.5, 9.5);
    airNozzleBarrel.rotation.x = 0.88;
    const airCollar = new THREE.Mesh(
      new THREE.CylinderGeometry(1.68, 1.68, 2.2, 24),
      new THREE.MeshStandardMaterial({ color: 0x38bdf8, metalness: 0.75, roughness: 0.25 }),
    );
    airCollar.position.set(0, -8.5, 7.8);
    airCollar.rotation.x = 0.88;
    const airStreamMesh = new THREE.Mesh(
      new THREE.ConeGeometry(2.4, 13.5, 32, 1, true),
      new THREE.MeshBasicMaterial({
        color: 0x7dd3fc,
        transparent: true,
        opacity: 0.36,
        blending: THREE.AdditiveBlending,
        depthWrite: false,
        side: THREE.DoubleSide,
      }),
    );
    airStreamMesh.position.set(0, -3.5, 2.5);
    airStreamMesh.rotation.x = 0.88;
    this.airGroup.add(airNozzleBarrel, airCollar, airStreamMesh);
    this.airGroup.visible = false;
    this.scene.add(this.airGroup);

    // 5. Studio Condenser Microphone Capsule
    const micBody = new THREE.Mesh(
      new THREE.CylinderGeometry(1.6, 1.6, 8.5, 24),
      this.matTitaniumCage,
    );
    micBody.position.set(0, 5.5, 5.5);
    micBody.rotation.x = -0.65;
    const micGrille = new THREE.Mesh(
      new THREE.SphereGeometry(1.65, 20, 16),
      new THREE.MeshStandardMaterial({ color: 0xd8b458, metalness: 0.9, roughness: 0.25, wireframe: true }),
    );
    micGrille.position.set(0, 2.2, 2.5);
    this.micGroup.add(micBody, micGrille);
    this.scene.add(this.micGroup);

    return {
      torchFlameCore,
      torchFlamePlume,
      torchPointLight,
      ln2StreamMesh,
      ln2FrostCloud,
      oilStreamMesh,
      airStreamMesh,
    };
  }

  /**
   * Construct 5 GPU-instanced 3D Vector Field layers:
   * - `heat`: Fourier heat conduction flux q = -k * grad(T) & radiative/plume flux
   * - `convection`: Buoyant Rayleigh-Bénard plumes, impinging wall-jets, and Marangoni/free-surface flow
   * - `magnetic`: 3D dipolar magnetic flux density B(r) + ferromagnetic remanence + AC degauss field
   * - `acoustic`: Spherical Helmholtz/Poynting acoustic wave radiation I(r, t) toward the 3D microphone
   * - `stress`: 3D Case Shell viscoplastic/capillary forces + Going-Train torque path + 3D Gravity g
   */
  private buildVectorFieldOverlays(): void {
    const arrowGeom = buildUnitArrowGeometry();
    const makeLayerMesh = (emissiveHex: number): THREE.InstancedMesh => {
      const mat = new THREE.MeshStandardMaterial({
        color: 0xffffff,
        emissive: new THREE.Color(emissiveHex),
        emissiveIntensity: 0.32,
        roughness: 0.24,
        metalness: 0.18,
        transparent: true,
        opacity: 0.88,
      });
      const mesh = new THREE.InstancedMesh(arrowGeom, mat, VEC_FIELD_LATTICE_COUNT);
      mesh.instanceMatrix.setUsage(THREE.DynamicDrawUsage);
      const zeroMat = new THREE.Matrix4().makeScale(0, 0, 0);
      const defaultCol = new THREE.Color(0xffffff);
      for (let i = 0; i < VEC_FIELD_LATTICE_COUNT; i++) {
        mesh.setMatrixAt(i, zeroMat);
        mesh.setColorAt(i, defaultCol);
      }
      mesh.visible = false;
      mesh.frustumCulled = false;
      return mesh;
    };

    this.vecHeatMesh = makeLayerMesh(0xff5518);
    this.vecConvectionMesh = makeLayerMesh(0x18d8b4);
    this.vecMagneticMesh = makeLayerMesh(0xa855f7);
    this.vecAcousticMesh = makeLayerMesh(0x38bdf8);
    this.vecStressMesh = makeLayerMesh(0xf59e0b);

    this.vectorFieldsGroup.add(
      this.vecHeatMesh,
      this.vecConvectionMesh,
      this.vecMagneticMesh,
      this.vecAcousticMesh,
      this.vecStressMesh,
    );
    this.scene.add(this.vectorFieldsGroup);

    // Precompute a golden-spiral + multi-ring 3D sampling lattice (168 points across 2 z-planes + case rim)
    this.vecLatticePts.length = 0;
    const rings: Array<{ r: number; n: number; z: number }> = [
      { r: 4.5, n: 8, z: 2.1 },
      { r: 9.8, n: 16, z: 2.1 },
      { r: 15.2, n: 24, z: 2.1 },
      { r: 20.4, n: 28, z: 2.2 },
      { r: 25.4, n: 36, z: 1.8 },
      { r: 8.5, n: 14, z: 5.4 },
      { r: 16.5, n: 20, z: 5.6 },
      { r: 23.8, n: 22, z: 5.2 },
    ];
    for (const ring of rings) {
      for (let k = 0; k < ring.n; k++) {
        if (this.vecLatticePts.length >= VEC_FIELD_LATTICE_COUNT) break;
        const a = (k / ring.n) * Math.PI * 2 + ring.z * 0.27;
        this.vecLatticePts.push([Math.cos(a) * ring.r, Math.sin(a) * ring.r, ring.z]);
      }
    }
    while (this.vecLatticePts.length < VEC_FIELD_LATTICE_COUNT) {
      this.vecLatticePts.push([0, -11, 2.2]);
    }
  }

  /**
   * Update the high-resolution `128 x 22` Bassine Case Shell vertices by Catmull-Rom spline
   * interpolation of the `64 x 8` Lagrangian continuum nodes computed by the Rust
   * Thermo-Elasto-Viscoplastic + Surface-Tension solver (`sim.case_shell_snapshot()`).
   */
  private updateCaseShellFromContinuum(caseShellData: Float32Array): void {
    if (caseShellData.length < CASE_SHELL_NODES * 6) return;

    const posAttr = this.caseShellGeometry.getAttribute("position") as THREE.BufferAttribute;
    const colAttr = this.caseShellGeometry.getAttribute("color") as THREE.BufferAttribute;
    const pos = posAttr.array as Float32Array;
    const col = colAttr.array as Float32Array;

    const readRawNode = (izClamped: number, tWrapped: number, ch: number): number =>
      caseShellData[(izClamped * CASE_N_THETA + tWrapped) * 6 + ch]!;

    const sampleNodeField = (iz: number, it: number, ch: number): number => {
      const tWrapped = ((it % CASE_N_THETA) + CASE_N_THETA) % CASE_N_THETA;
      // Phantom endpoint reflection for spatial coordinates (ch <= 2) so the C^2 cubic B-spline
      // interpolates the exact top and bottom bezel rims (z = +/-3.8 mm, R = 23.0 mm) with zero gap:
      if (ch <= 2) {
        if (iz < 0) {
          return 2.0 * readRawNode(0, tWrapped, ch) - readRawNode(1, tWrapped, ch);
        }
        if (iz >= CASE_N_Z) {
          return 2.0 * readRawNode(CASE_N_Z - 1, tWrapped, ch) - readRawNode(CASE_N_Z - 2, tWrapped, ch);
        }
      }
      const zClamped = Math.max(0, Math.min(CASE_N_Z - 1, iz));
      return readRawNode(zClamped, tWrapped, ch);
    };

    const sampleShellSmooth = (uTheta: number, vZ: number, ch: number): number => {
      const ft = uTheta * CASE_N_THETA;
      const fz = vZ * (CASE_N_Z - 1);
      const it0 = Math.floor(ft);
      const iz0 = Math.floor(fz);
      const tt = ft - it0;
      const tz = fz - iz0;

      const rowVals = [0, 0, 0, 0];
      for (let dz = -1; dz <= 2; dz++) {
        const p0 = sampleNodeField(iz0 + dz, it0 - 1, ch);
        const p1 = sampleNodeField(iz0 + dz, it0, ch);
        const p2 = sampleNodeField(iz0 + dz, it0 + 1, ch);
        const p3 = sampleNodeField(iz0 + dz, it0 + 2, ch);
        rowVals[dz + 1] = bSplineCubic(p0, p1, p2, p3, tt);
      }
      return bSplineCubic(rowVals[0]!, rowVals[1]!, rowVals[2]!, rowVals[3]!, tz);
    };

    let maxShellTempK = 293.15;

    for (let iz = 0; iz < CASE_RENDER_N_Z; iz++) {
      const vZ = iz / (CASE_RENDER_N_Z - 1);
      const rimPinch = Math.sin(vZ * Math.PI);
      const bezelLip = 0.72 + 0.28 * Math.pow(rimPinch, 0.45);

      for (let it = 0; it < CASE_RENDER_N_THETA; it++) {
        const uTheta = it / CASE_RENDER_N_THETA;
        const idx = iz * CASE_RENDER_N_THETA + it;

        const x = sampleShellSmooth(uTheta, vZ, 0);
        const y = sampleShellSmooth(uTheta, vZ, 1);
        const z = sampleShellSmooth(uTheta, vZ, 2);
        const h = Math.max(0.08, sampleShellSmooth(uTheta, vZ, 3)) * bezelLip;
        const tempK = Math.max(50.0, sampleShellSmooth(uTheta, vZ, 4));
        const meltFrac = Math.max(0.0, Math.min(1.0, sampleShellSmooth(uTheta, vZ, 5)));
        if (tempK > maxShellTempK) maxShellTempK = tempK;

        const rXY = Math.max(1.0, Math.hypot(x, y));
        const nx = x / rXY;
        const ny = y / rXY;

        // Outer shell surface at (x + 0.45*h*nx, y + 0.45*h*ny, z)
        const oIdx = idx * 3;
        pos[oIdx] = x + nx * h * 0.45;
        pos[oIdx + 1] = y + ny * h * 0.45;
        pos[oIdx + 2] = z;

        // Inner shell surface at (x - 0.55*h*nx, y - 0.55*h*ny, z)
        const iIdx = (CASE_RENDER_NODES + idx) * 3;
        pos[iIdx] = x - nx * h * 0.55;
        pos[iIdx + 1] = y - ny * h * 0.55;
        pos[iIdx + 2] = z;

        // Vertex color: 18k Yellow Gold reflectance strictly in [0, 1] (preserves Schlick Fresnel energy conservation)
        let cr = 0.93;
        let cg = 0.75;
        let cb = 0.32;
        if (this.thermalOverlayMode) {
          const u = Math.max(0, Math.min(1, (tempK - 77.36) / 1400.0));
          cr = Math.min(1, u * 2.0);
          cg = Math.sin(u * Math.PI) * 0.85;
          cb = Math.max(0, 1.0 - u * 2.0);
        } else if (tempK > 720.0) {
          const glow = Math.min(1.0, (tempK - 720.0) / 680.0);
          cr = Math.min(1.0, 0.93 + glow * 0.07);
          cg = Math.min(0.95, 0.75 - glow * 0.14 * (1.0 - meltFrac) + glow * 0.16 * meltFrac);
          cb = Math.min(0.72, 0.32 + glow * 0.34 * meltFrac);
        } else if (tempK < 220.0) {
          const cryo = Math.min(1.0, (220.0 - tempK) / 145.0);
          cr = 0.93 - cryo * 0.22;
          cg = 0.75 + cryo * 0.14;
          cb = 0.32 + cryo * 0.52;
        }

        col[oIdx] = cr;
        col[oIdx + 1] = cg;
        col[oIdx + 2] = cb;
        col[iIdx] = cr * 0.90;
        col[iIdx + 1] = cg * 0.90;
        col[iIdx + 2] = cb * 0.90;
      }
    }

    // Subtle uniform blackbody self-emission when case exceeds 820 K
    this.matGoldCase.emissiveIntensity =
      maxShellTempK > 820.0 ? Math.min(0.28, ((maxShellTempK - 820.0) / 650.0) * 0.28) : 0.0;

    posAttr.needsUpdate = true;
    colAttr.needsUpdate = true;
    this.caseShellGeometry.computeVertexNormals();
  }

  /**
   * Update the continuous 3D Molten & Solidified Gold Free-Surface Heightfields using C^2
   * Bicubic B-Spline reconstruction (`84 x 84`) over the `48 x 48 = 2304` Spatial PDE Grid
   * computed in Rust (`sim.spatial_grid_snapshot()`), with per-pixel RGBA iso-contour meniscus
   * feathering so pools never show boxy Cartesian quad boundaries.
   */
  private updateGoldHeightfieldsFromSpatialGrid(gridData: Float32Array): void {
    if (gridData.length < GRID_N * GRID_N * 7) return;

    // Fast coarse scan to avoid 84x84 spline evaluation when no molten or frozen gold exists
    let maxLiqCourse = 0.0;
    let maxFrzCourse = 0.0;
    for (let i = 0; i < GRID_N * GRID_N; i++) {
      const liq = gridData[i * 7 + 3]!;
      const frz = gridData[i * 7 + 4]!;
      if (liq > maxLiqCourse) maxLiqCourse = liq;
      if (frz > maxFrzCourse) maxFrzCourse = frz;
    }

    if (maxLiqCourse <= 0.008 && maxFrzCourse <= 0.008) {
      this.moltenGoldMesh.visible = false;
      this.solidGoldMesh.visible = false;
      return;
    }

    const mPosAttr = this.moltenGoldGeometry.getAttribute("position") as THREE.BufferAttribute;
    const mColAttr = this.moltenGoldGeometry.getAttribute("color") as THREE.BufferAttribute;
    const mPos = mPosAttr.array as Float32Array;
    const mCol = mColAttr.array as Float32Array;

    const sPosAttr = this.solidGoldGeometry.getAttribute("position") as THREE.BufferAttribute;
    const sColAttr = this.solidGoldGeometry.getAttribute("color") as THREE.BufferAttribute;
    const sPos = sPosAttr.array as Float32Array;
    const sCol = sColAttr.array as Float32Array;

    const totalVerts = RENDER_GRID_N * RENDER_GRID_N;
    let anyMolten = false;
    let anySolid = false;

    for (let idx = 0; idx < totalVerts; idx++) {
      const xMm = mPos[idx * 3]!;
      const yMm = mPos[idx * 3 + 1]!;
      const rCenter = Math.hypot(xMm, yMm);

      // Smooth circular radial taper near the outer case rim (r <= 25.0 mm)
      const radialT = Math.max(0.0, Math.min(1.0, (25.0 - rCenter) / 2.6));
      const radialMask = radialT * radialT * (3.0 - 2.0 * radialT);

      if (radialMask <= 0.001) {
        mPos[idx * 3 + 2] = -2.2;
        mCol[idx * 4 + 3] = 0.0;
        sPos[idx * 3 + 2] = -2.2;
        sCol[idx * 4 + 3] = 0.0;
        continue;
      }

      const tempK = sampleSmoothPdeField(gridData, xMm, yMm, 0);
      const rawLiq = sampleSmoothPdeField(gridData, xMm, yMm, 3) * radialMask;
      const rawFrz = sampleSmoothPdeField(gridData, xMm, yMm, 4) * radialMask;

      // Smooth substrate elevation transition between recessed 6h Tourbillon well (-0.65 mm) and mainplate (+0.25 mm)
      const distTourb = Math.hypot(xMm, yMm + 11.0);
      const wellBlend = Math.max(0.0, Math.min(1.0, (distTourb - 5.2) / 2.8));
      const smoothWell = wellBlend * wellBlend * (3.0 - 2.0 * wellBlend);
      const baseZ = -0.65 * (1.0 - smoothWell) + 0.28 * smoothWell;

      // 1. Molten Gold Free-Surface Height & Smooth Per-Pixel Meniscus Alpha
      const liqAlphaT = Math.max(0.0, Math.min(1.0, (rawLiq - 0.006) / 0.042));
      const liqAlpha = liqAlphaT * liqAlphaT * (3.0 - 2.0 * liqAlphaT);
      if (liqAlpha > 0.005) {
        anyMolten = true;
      }
      const liqDome = Math.pow(liqAlpha, 0.65);
      mPos[idx * 3 + 2] =
        baseZ - 0.42 * (1.0 - liqDome) + rawFrz * 1.25 + Math.min(2.4, rawLiq * 1.32) * liqDome;

      const glow = Math.max(0.0, Math.min(1.0, (tempK - 950.0) / 550.0));
      mCol[idx * 4] = Math.min(1.0, 0.96 + glow * 0.04);
      mCol[idx * 4 + 1] = Math.min(0.92, 0.74 + glow * 0.16);
      mCol[idx * 4 + 2] = Math.min(0.68, 0.28 + glow * 0.32);
      mCol[idx * 4 + 3] = liqAlpha;

      // 2. Solidified (Quenched) Gold Free-Surface Height & Smooth Per-Pixel Meniscus Alpha
      const frzAlphaT = Math.max(0.0, Math.min(1.0, (rawFrz - 0.006) / 0.042));
      const frzAlpha = frzAlphaT * frzAlphaT * (3.0 - 2.0 * frzAlphaT);
      if (frzAlpha > 0.005) {
        anySolid = true;
      }
      const frzDome = Math.pow(frzAlpha, 0.65);
      sPos[idx * 3 + 2] =
        baseZ - 0.42 * (1.0 - frzDome) + Math.min(2.4, rawFrz * 1.32) * frzDome;

      const cryo = tempK < 220.0 ? Math.min(1.0, (220.0 - tempK) / 145.0) : 0.0;
      sCol[idx * 4] = 0.88 - cryo * 0.16;
      sCol[idx * 4 + 1] = 0.69 + cryo * 0.12;
      sCol[idx * 4 + 2] = 0.28 + cryo * 0.45;
      sCol[idx * 4 + 3] = frzAlpha;
    }

    this.moltenGoldMesh.visible = anyMolten;
    if (anyMolten) {
      mPosAttr.needsUpdate = true;
      mColAttr.needsUpdate = true;
      this.moltenGoldGeometry.computeVertexNormals();
    }

    this.solidGoldMesh.visible = anySolid;
    if (anySolid) {
      sPosAttr.needsUpdate = true;
      sColAttr.needsUpdate = true;
      this.solidGoldGeometry.computeVertexNormals();
    }
  }

  public renderFrame(snap: WatchSnapshot, tools: ToolRenderParams): void {
    // 1. Update Emergent Continuum Geometries from WASM buffers
    if (tools.caseShellData) {
      this.updateCaseShellFromContinuum(tools.caseShellData);
    }
    if (tools.spatialGridData) {
      this.updateGoldHeightfieldsFromSpatialGrid(tools.spatialGridData);
    }
    if (tools.mainspringCoilsMm) {
      this.updateMainspringRibbonFromCoils(tools.mainspringCoilsMm);
    }

    // 2. Update Horological Multi-Body Rotations
    this.crownGroup.rotation.y = snap.barrelArborAngleRad * 2.5;
    this.crownWheelGroup.rotation.z = -(tools.crownWheelAngleRad ?? snap.barrelArborAngleRad * 1.85);
    this.ratchetWheelGroup.rotation.z = snap.barrelArborAngleRad;
    this.clickLeverGroup.rotation.z = -(tools.clickLeverAngleRad ?? 0.0);
    this.barrelWheelGroup.rotation.z = snap.barrelWheelAngleRad;
    this.centerWheelGroup.rotation.z = snap.centerWheelAngleRad;
    this.thirdWheelGroup.rotation.z = snap.thirdWheelAngleRad;

    this.tourbillonCageGroup.rotation.z = -snap.tourbillonCageAngleRad;
    this.escapeWheelGroup.rotation.z = -snap.escapeWheelVisualRad;
    this.palletForkGroup.rotation.z = snap.palletForkAngleRad;
    this.balanceWheelGroup.rotation.z = snap.balanceAngleRad;

    const oilBridgeFrac = Math.min(0.38, Math.max(0, (snap.nodes[11]!.oilFilmUm - 2.5) / 38.0));
    this.updateHairspringRibbon(snap.balanceAngleRad, oilBridgeFrac);

    // 3. Update Blued-Steel Thin-Film Interference Color on Hairspring & Hands
    this.matHairspring.color.setRGB(
      snap.bluedSteelRgb[0],
      snap.bluedSteelRgb[1],
      snap.bluedSteelRgb[2],
    );

    // 4. Update Dial Hands (12-hour watchDisplayedTimeS + Power Reserve + Moonphase)
    const tSec = snap.watchDisplayedTimeS;
    const secAngle = -((tSec % 60) / 60) * Math.PI * 2;
    const minAngle = -(((tSec / 60) % 60) / 60) * Math.PI * 2;
    const hourAngle = -(((tSec / 3600) % 12) / 12) * Math.PI * 2;
    this.secondsHandGroup.rotation.z = secAngle;
    this.minuteHandGroup.rotation.z = minAngle;
    this.hourHandGroup.rotation.z = hourAngle;
    this.powerReserveHandGroup.rotation.z = (0.5 - snap.mainspringWind) * 2.1;
    this.moonDiscGroup.rotation.z = snap.moonphaseAngleRad;

    // 5. Update Frost Rime Overlay when cooled below 270 K
    const maxFrostUm = Math.max(...snap.nodes.map((n) => n.frostThicknessUm));
    if (maxFrostUm > 1.5) {
      this.frostOverlayMesh.visible = true;
      this.matFrostOverlay.opacity = Math.min(0.62, maxFrostUm / 135.0);
    } else {
      this.frostOverlayMesh.visible = false;
    }

    // 6. Update Perturbation Tool Apparatus Positions, Dynamic Bloom & Intensities
    const maxNodeTempK = Math.max(...snap.nodes.map((n) => n.temperatureK));
    const incandescenceFactor = Math.max(0, Math.min(1, (maxNodeTempK - 820.0) / 550.0));
    this.bloomPass.strength = 0.06 + Math.max(incandescenceFactor * 0.12, tools.torchPower * 0.08);

    this.torchGroup.visible = tools.torchPower > 0.01 || this.activeTool === "torch";
    this.torchGroup.position.set(tools.torchPosMm[0], tools.torchPosMm[1], tools.torchPosMm[2]);
    const flameActive = tools.torchPower > 0.01;
    this.torchFlameCore.visible = flameActive;
    this.torchFlamePlume.visible = flameActive;
    const flameScale = Math.max(0.08, tools.torchPower);
    this.torchFlameCore.scale.set(flameScale, flameScale, flameScale);
    this.torchFlamePlume.scale.set(
      flameScale * (0.96 + Math.sin(performance.now() * 0.045) * 0.05),
      flameScale * (0.94 + Math.cos(performance.now() * 0.062) * 0.07),
      flameScale,
    );
    this.torchPointLight.intensity = tools.torchPower * 2.2 + incandescenceFactor * 0.9;

    this.ln2Group.visible = tools.ln2Rate > 0.01 || this.activeTool === "ln2";
    this.ln2Group.position.set(tools.ln2PosMm[0], tools.ln2PosMm[1], tools.ln2PosMm[2]);
    this.ln2StreamMesh.visible = tools.ln2Rate > 0.01;
    this.ln2FrostCloud.visible = tools.ln2Rate > 0.01 || snap.ln2BoilingPowerW > 1.0;
    this.ln2FrostCloud.rotation.z += 0.03;

    this.oilGroup.visible = tools.oilRate > 0.01 || this.activeTool === "oil";
    this.oilGroup.position.set(tools.oilPosMm[0], tools.oilPosMm[1], tools.oilPosMm[2]);
    this.oilStreamMesh.visible = tools.oilRate > 0.01;

    this.magnetGroup.visible = tools.magnetOn || this.activeTool === "magnet";
    this.magnetGroup.position.set(tools.magnetPosMm[0], tools.magnetPosMm[1], tools.magnetPosMm[2]);

    const airPower = tools.airPower ?? 0.0;
    const airPosMm = tools.airPosMm ?? [0.0, -20.0, 3.0];
    this.airGroup.visible = airPower > 0.01 || this.activeTool === "air";
    this.airGroup.position.set(airPosMm[0], airPosMm[1], airPosMm[2]);
    this.airStreamMesh.visible = airPower > 0.01;
    const airScale = Math.max(0.1, airPower);
    this.airStreamMesh.scale.set(
      airScale * (0.94 + Math.sin(performance.now() * 0.06) * 0.06),
      airScale,
      airScale * (0.94 + Math.cos(performance.now() * 0.05) * 0.06),
    );

    this.micGroup.position.set(tools.micPosMm[0], tools.micPosMm[1], tools.micPosMm[2]);

    // 7. Update Superimposed 3D Physics Vector Fields when any layer is active
    this.updateVectorFieldOverlays(snap, tools);

    this.controls.update();
    this.composer.render();
  }

  private setArrowInstance(
    mesh: THREE.InstancedMesh,
    index: number,
    px: number,
    py: number,
    pz: number,
    vx: number,
    vy: number,
    vz: number,
    lengthMm: number,
    thicknessMm: number,
    rCol: number,
    gCol: number,
    bCol: number,
  ): void {
    const mag = Math.hypot(vx, vy, vz);
    if (mag < 1.0e-6 || lengthMm < 0.15) {
      this.tmpMat4.makeScale(0, 0, 0);
      mesh.setMatrixAt(index, this.tmpMat4);
      return;
    }
    this.tmpPos.set(px, py, pz);
    this.tmpDir.set(vx / mag, vy / mag, vz / mag);
    this.tmpQuat.setFromUnitVectors(VEC_UP, this.tmpDir);
    this.tmpScale.set(thicknessMm, lengthMm, thicknessMm);
    this.tmpMat4.compose(this.tmpPos, this.tmpQuat, this.tmpScale);
    mesh.setMatrixAt(index, this.tmpMat4);
    this.tmpColor.setRGB(
      Math.max(0, Math.min(1, rCol)),
      Math.max(0, Math.min(1, gCol)),
      Math.max(0, Math.min(1, bCol)),
    );
    mesh.setColorAt(index, this.tmpColor);
  }

  /**
   * Compute and update the 5 superimposed 3D physics vector fields directly from the live
   * PDE spatial grid, 3D Lagrangian case shell, magnetic dipole + remanence, acoustic wave
   * radiation, and going-train / gravity forces.
   */
  private updateVectorFieldOverlays(snap: WatchSnapshot, tools: ToolRenderParams): void {
    const tNow = performance.now() * 0.001;
    const grid = tools.spatialGridData;
    const shell = tools.caseShellData;
    const airPower = tools.airPower ?? 0.0;
    const airPos = tools.airPosMm ?? [0.0, -20.0, 3.0];
    const gVec = tools.gravityGVec ?? [0.0, -1.0, 0.0];
    const pAtm = tools.ambientPressureAtm ?? 1.0;

    // 1. Heat Flux Vector Field: q = -k * grad(T) + buoyant/radiative axial heat flux
    if (this.vecFieldVisible.heat && grid && grid.length >= GRID_N * GRID_N * 7) {
      const hStep = 1.6;
      for (let i = 0; i < VEC_FIELD_LATTICE_COUNT; i++) {
        const [x, y, z] = this.vecLatticePts[i]!;
        const tC = sampleSmoothPdeField(grid, x, y, 0);
        const tXp = sampleSmoothPdeField(grid, x + hStep, y, 0);
        const tXm = sampleSmoothPdeField(grid, x - hStep, y, 0);
        const tYp = sampleSmoothPdeField(grid, x, y + hStep, 0);
        const tYm = sampleSmoothPdeField(grid, x, y - hStep, 0);
        // Fourier heat conduction vector q_xy = -grad(T)
        let qx = -(tXp - tXm) / (2.0 * hStep);
        let qy = -(tYp - tYm) / (2.0 * hStep);
        let qz = (tC - 293.15) * 0.22;

        // Add subtle thermoelastic dissipation wave around 6h balance & 12h barrel at room temp
        const dBal = Math.hypot(x, y + 11.0);
        if (Math.abs(tC - 293.15) < 1.5 && Math.hypot(qx, qy) < 0.25) {
          const ang = Math.atan2(y + 11.0, x);
          const pulse = 0.55 + 0.45 * Math.sin(tNow * 15.7 - dBal * 0.45);
          qx += Math.cos(ang) * 0.9 * pulse * Math.exp(-dBal / 16.0);
          qy += Math.sin(ang) * 0.9 * pulse * Math.exp(-dBal / 16.0);
          qz += 0.35 * pulse;
        }

        const qMag = Math.hypot(qx, qy, qz);
        const lenMm = Math.min(5.2, 1.35 + 1.15 * Math.log1p(qMag * 0.45));
        const pulseShift = 0.28 * Math.sin(tNow * 5.5 - i * 0.17);
        const px = x + (qx / Math.max(0.01, qMag)) * pulseShift;
        const py = y + (qy / Math.max(0.01, qMag)) * pulseShift;

        // Color: cryogenic cyan (<285 K) -> warm amber (300-600 K) -> incandescent white-gold (>950 K)
        let rC = 1.0;
        let gC = 0.58;
        let bC = 0.16;
        if (tC < 285.0) {
          const cold = Math.min(1.0, (293.15 - tC) / 140.0);
          rC = 0.18 * (1.0 - cold) + 0.15;
          gC = 0.78 + 0.2 * cold;
          bC = 1.0;
        } else {
          const hot = Math.min(1.0, (tC - 293.15) / 650.0);
          rC = 1.0;
          gC = 0.42 + 0.52 * hot;
          bC = 0.12 + 0.65 * Math.max(0, hot - 0.35);
        }
        this.setArrowInstance(this.vecHeatMesh, i, px, py, z, qx, qy, qz, lenMm, 1.05, rC, gC, bC);
      }
      this.vecHeatMesh.instanceMatrix.needsUpdate = true;
      if (this.vecHeatMesh.instanceColor) this.vecHeatMesh.instanceColor.needsUpdate = true;
    }

    // 2. Convection, Impinging Jet & Free-Surface Flow Field: u_flow(x, y, z)
    if (this.vecFieldVisible.convection && grid && grid.length >= GRID_N * GRID_N * 7) {
      for (let i = 0; i < VEC_FIELD_LATTICE_COUNT; i++) {
        const [x, y, z] = this.vecLatticePts[i]!;
        const tC = sampleSmoothPdeField(grid, x, y, 0);
        const hLiq = sampleSmoothPdeField(grid, x, y, 3);

        // (a) Buoyant Rayleigh-Bénard natural convection plume (scaled by ambient pressure)
        const buoyZ = (tC - 293.15) * 0.035 * Math.sqrt(Math.max(0.02, pAtm));
        let ux = -y * 0.04 * Math.sqrt(Math.max(0.05, pAtm));
        let uy = x * 0.04 * Math.sqrt(Math.max(0.05, pAtm));
        let uz = buoyZ;

        // (b) Oscillating 5 Hz Balance Wheel & 1-Minute Tourbillon aerodynamic swirl near (0, -11)
        const dxTourb = x;
        const dyTourb = y + 11.0;
        const rTourb = Math.hypot(dxTourb, dyTourb);
        if (rTourb < 14.0 && pAtm > 0.02) {
          const swirl =
            (snap.balanceOmegaRadS * 0.035 + 0.45) *
            Math.exp(-(rTourb * rTourb) / (2.0 * 6.5 * 6.5)) *
            Math.sqrt(pAtm);
          ux += (-dyTourb / Math.max(1.0, rTourb)) * swirl;
          uy += (dxTourb / Math.max(1.0, rTourb)) * swirl;
        }

        // (c) Impinging Wall-Jet streamlines from Blowtorch, LN2, and Compressed Air Nozzle
        const addWallJet = (jx: number, jy: number, power: number, tempSign: number): void => {
          if (power <= 0.01) return;
          const dx = x - jx;
          const dy = y - jy;
          const r = Math.hypot(dx, dy);
          const w = Math.exp(-(r * r) / (2.0 * 11.5 * 11.5));
          const radOut = (r / (r + 2.5)) * power * 3.8 * w;
          ux += (dx / Math.max(0.5, r)) * radOut;
          uy += (dy / Math.max(0.5, r)) * radOut;
          uz += (z > 3.8 ? -2.2 : 1.1 * tempSign) * power * w;
        };
        addWallJet(tools.torchPosMm[0], tools.torchPosMm[1], tools.torchPower, 1.0);
        addWallJet(tools.ln2PosMm[0], tools.ln2PosMm[1], tools.ln2Rate, -1.0);
        addWallJet(airPos[0], airPos[1], airPower * 1.35, 0.2);

        // (d) Molten gold Marangoni + Gravity free-surface velocity when liquid gold is present
        if (hLiq > 0.01 && z < 3.5) {
          ux += gVec[0] * 1.8;
          uy += gVec[1] * 1.8;
        }

        const uMag = Math.hypot(ux, uy, uz);
        const lenMm = Math.min(5.0, 1.45 + 1.15 * Math.log1p(uMag));
        const phase = ((tNow * 2.4 + i * 0.13) % 1.0) - 0.5;
        const px = x + (ux / Math.max(0.01, uMag)) * phase * 0.9;
        const py = y + (uy / Math.max(0.01, uMag)) * phase * 0.9;

        const warm = Math.max(0, Math.min(1, (tC - 260.0) / 500.0));
        const rC = 0.12 + 0.85 * warm;
        const gC = 0.88 - 0.22 * Math.abs(warm - 0.4);
        const bC = 0.78 * (1.0 - warm * 0.65);
        this.setArrowInstance(this.vecConvectionMesh, i, px, py, z, ux, uy, uz, lenMm, 1.0, rC, gC, bC);
      }
      this.vecConvectionMesh.instanceMatrix.needsUpdate = true;
      if (this.vecConvectionMesh.instanceColor) this.vecConvectionMesh.instanceColor.needsUpdate = true;
    }

    // 3. 3D Magnetic Flux Density Vector Field: B(r) = B_dipole(r) + B_remanence(r) + B_ferro(r)
    if (this.vecFieldVisible.magnetic) {
      const mx0 = tools.magnetOn ? tools.magnetPosMm[0] : 42.0;
      const my0 = tools.magnetOn ? tools.magnetPosMm[1] : -11.0;
      const mz0 = tools.magnetOn ? tools.magnetPosMm[2] : 3.0;
      const mStrength = tools.magnetOn ? 1.0 : 0.12;
      const remStrength = snap.remanentMagnetizationMt * 0.08;

      for (let i = 0; i < VEC_FIELD_LATTICE_COUNT; i++) {
        const [x, y, z] = this.vecLatticePts[i]!;
        // Exact 3D magnetic dipole field B(r) proportional to (3 (m . r_hat) r_hat - m) / |r|^3
        // with dipole moment m = (-1, 0, 0) along the NdFeB North-South axis
        const rx = x - mx0;
        const ry = y - my0;
        const rz = z - mz0;
        const r2 = Math.max(28.0, rx * rx + ry * ry + rz * rz);
        const r = Math.sqrt(r2);
        const r5 = r2 * r2 * r;
        const mDotR = -rx * mStrength * 4500.0;
        let bx = (3.0 * mDotR * rx - -1.0 * mStrength * 4500.0 * r2) / r5;
        let by = (3.0 * mDotR * ry) / r5;
        let bz = (3.0 * mDotR * rz) / r5;

        // Add residual ferromagnetic remanence dipole centered on the 6h steel pallet fork & balance
        if (remStrength > 0.01 || !tools.magnetOn) {
          const ex = x;
          const ey = y + 11.0;
          const ez = z - 0.5;
          const er2 = Math.max(16.0, ex * ex + ey * ey + ez * ez);
          const er = Math.sqrt(er2);
          const er5 = er2 * er2 * er;
          const effRem = Math.max(0.18, remStrength);
          const emDotR = (ex * 0.85 + ey * 0.52) * effRem * 950.0;
          bx += (3.0 * emDotR * ex - 0.85 * effRem * 950.0 * er2) / er5;
          by += (3.0 * emDotR * ey - 0.52 * effRem * 950.0 * er2) / er5;
          bz += (3.0 * emDotR * ez) / er5;
        }

        const bMag = Math.hypot(bx, by, bz);
        const lenMm = Math.min(4.8, 1.55 + 1.35 * Math.log1p(bMag * 4.5));
        const normB = Math.min(1.0, Math.log1p(bMag * 8.0) / 2.6);
        const rC = 0.52 + 0.48 * normB;
        const gC = 0.26 + 0.58 * (1.0 - Math.abs(normB - 0.5));
        const bC = 1.0;
        this.setArrowInstance(this.vecMagneticMesh, i, x, y, z, bx, by, bz, lenMm, 1.02, rC, gC, bC);
      }
      this.vecMagneticMesh.instanceMatrix.needsUpdate = true;
      if (this.vecMagneticMesh.instanceColor) this.vecMagneticMesh.instanceColor.needsUpdate = true;
    }

    // 4. 3D Acoustic Wave & Poynting Intensity Field: I_ac(r, t) = p(r, t) * v_ac(r, t)
    if (this.vecFieldVisible.acoustic) {
      const airTrans = Math.max(0.12, Math.min(1.3, pAtm / (pAtm + 0.3)));
      for (let i = 0; i < VEC_FIELD_LATTICE_COUNT; i++) {
        const [x, y, z] = this.vecLatticePts[i]!;
        // Primary spherical wave source at 6h Ruby Escapement (0, -11, 0.5)
        const sx = x;
        const sy = y + 11.0;
        const sz = z - 0.5;
        const rEsc = Math.max(2.5, Math.hypot(sx, sy, sz));
        const kWave = 0.38;
        const waveEnvelope =
          (0.65 + 0.35 * Math.cos(kWave * rEsc - tNow * 31.4159)) *
          (12.0 / (rEsc + 6.0)) *
          airTrans;

        let ax = (sx / rEsc) * waveEnvelope;
        let ay = (sy / rEsc) * waveEnvelope;
        let az = ((sz + 2.2) / rEsc) * waveEnvelope;

        // Secondary acoustic dipole from active Blowtorch / LN2 / Air Nozzle
        if (tools.torchPower > 0.01) {
          const tx = x - tools.torchPosMm[0];
          const ty = y - tools.torchPosMm[1];
          const rT = Math.max(3.0, Math.hypot(tx, ty));
          const wT = (tools.torchPower * 4.5) / (rT + 4.0);
          ax += (tx / rT) * wT;
          ay += (ty / rT) * wT;
        }

        const aMag = Math.hypot(ax, ay, az);
        const lenMm = Math.min(4.8, 1.2 + 2.2 * aMag);
        const crest = 0.5 + 0.5 * Math.sin(kWave * rEsc - tNow * 31.4159);
        const rC = 0.18 + 0.72 * crest;
        const gC = 0.78 + 0.2 * crest;
        const bC = 1.0 - 0.35 * crest;
        this.setArrowInstance(this.vecAcousticMesh, i, x, y, z, ax, ay, az, lenMm, 1.0, rC, gC, bC);
      }
      this.vecAcousticMesh.instanceMatrix.needsUpdate = true;
      if (this.vecAcousticMesh.instanceColor) this.vecAcousticMesh.instanceColor.needsUpdate = true;
    }

    // 5. Structural Stress, Viscoplastic Case Force, Going-Train Torque & 3D Gravity Vector Field
    if (this.vecFieldVisible.stress) {
      let idx = 0;
      // (a) 64 perimeter arrows around the 18k Gold Case Shell showing net inward capillary +
      //     viscoplastic slump + thermal expansion + gravity force on the 3D continuum shell
      if (shell && shell.length >= CASE_SHELL_NODES * 6) {
        for (let it = 0; it < CASE_SHELL_N_THETA && idx < 64; it++) {
          const base = (4 * CASE_SHELL_N_THETA + it) * 6;
          const sx = shell[base]!;
          const sy = shell[base + 1]!;
          const sz = shell[base + 2]! + 2.2;
          const tempK = shell[base + 4]!;
          const melt = shell[base + 5]!;
          const rNow = Math.hypot(sx, sy);
          const rRef = 26.8;
          const radialInward = Math.max(0.0, rRef - rNow);
          const nx = -sx / Math.max(1.0, rNow);
          const ny = -sy / Math.max(1.0, rNow);
          const softness = Math.min(1.0, melt + Math.max(0, (tempK - 550.0) / 700.0));
          const fx = nx * (0.65 + 2.5 * softness + radialInward * 0.4) + gVec[0] * 0.45;
          const fy = ny * (0.65 + 2.5 * softness + radialInward * 0.4) + gVec[1] * 0.45;
          const fz = -0.25 * softness + gVec[2] * 0.35;
          const lenMm = Math.min(5.2, 1.65 + 2.4 * softness + radialInward * 0.35);
          const rC = 0.98;
          const gC = Math.max(0.22, 0.78 - 0.55 * softness);
          const bC = 0.15 + 0.25 * softness;
          this.setArrowInstance(this.vecStressMesh, idx++, sx, sy, sz, fx, fy, fz, lenMm, 1.08, rC, gC, bC);
        }
      }
      // (b) Going-train mechanical torque transmission + 3D Gravity load field across interior points
      const torqueScale = Math.max(0.25, snap.mainspringTorqueNm / 0.0082);
      for (; idx < VEC_FIELD_LATTICE_COUNT; idx++) {
        const [x, y, z] = this.vecLatticePts[idx]!;
        // Tangential gear-train force path from 12h Mainspring Barrel (0, 11.5) -> Center (0, 0) -> 6h Tourbillon (0, -11)
        const dBarrel = Math.hypot(x, y - 11.5);
        const dTourb = Math.hypot(x, y + 11.0);
        let fx = gVec[0] * 0.95;
        let fy = gVec[1] * 0.95 - 0.55 * torqueScale * Math.exp(-(x * x) / (2.0 * 9.0 * 9.0));
        let fz = gVec[2] * 0.95;
        if (dBarrel < 9.5) {
          fx += -(y - 11.5) * 0.22 * torqueScale;
          fy += x * 0.22 * torqueScale;
        } else if (dTourb < 9.0) {
          fx += (y + 11.0) * 0.25 * torqueScale;
          fy += -x * 0.25 * torqueScale;
        }
        const fMag = Math.hypot(fx, fy, fz);
        const lenMm = Math.min(4.6, 1.4 + 1.35 * fMag);
        this.setArrowInstance(
          this.vecStressMesh,
          idx,
          x,
          y,
          z,
          fx,
          fy,
          fz,
          lenMm,
          0.98,
          0.96,
          0.74,
          0.22,
        );
      }
      this.vecStressMesh.instanceMatrix.needsUpdate = true;
      if (this.vecStressMesh.instanceColor) this.vecStressMesh.instanceColor.needsUpdate = true;
    }
  }

  public setVectorFieldVisible(layer: VectorFieldLayer, visible: boolean): void {
    this.vecFieldVisible[layer] = visible;
    if (layer === "heat") this.vecHeatMesh.visible = visible;
    else if (layer === "convection") this.vecConvectionMesh.visible = visible;
    else if (layer === "magnetic") this.vecMagneticMesh.visible = visible;
    else if (layer === "acoustic") this.vecAcousticMesh.visible = visible;
    else if (layer === "stress") this.vecStressMesh.visible = visible;
  }

  public toggleVectorField(layer: VectorFieldLayer): boolean {
    const next = !this.vecFieldVisible[layer];
    this.setVectorFieldVisible(layer, next);
    return next;
  }

  public isVectorFieldVisible(layer: VectorFieldLayer): boolean {
    return this.vecFieldVisible[layer];
  }

  public setToolMode(mode: ActiveToolMode): void {
    this.activeTool = mode;
    this.controls.enabled = true;
    this.controls.enableZoom = true;
    this.controls.mouseButtons =
      mode === "inspect"
        ? {
            LEFT: THREE.MOUSE.ROTATE,
            MIDDLE: THREE.MOUSE.DOLLY,
            RIGHT: THREE.MOUSE.PAN,
          }
        : {
            LEFT: null as unknown as THREE.MOUSE,
            MIDDLE: THREE.MOUSE.DOLLY,
            RIGHT: THREE.MOUSE.ROTATE,
          };
  }

  public setExplodedFactor(factor: number): void {
    this.explodedFactor = Math.max(0, Math.min(1, factor));
    const s = this.explodedFactor * 14.0;
    this.layerCase.position.z = -s * 0.45;
    this.layerMainplate.position.z = -s * 0.25;
    this.layerBarrel.position.z = 0;
    this.layerGoingTrain.position.z = s * 0.22;
    this.layerKeyless.position.z = s * 0.36;
    this.layerBridges.position.z = s * 0.52;
    this.layerTourbillon.position.z = s * 0.75;
    this.layerDialHands.position.z = s * 1.05;
    this.layerCrystalFront.position.z = s * 1.35;
  }

  public setOpenCaseMode(open: boolean): void {
    this.openCaseMode = open;
    this.layerCrystalFront.visible = !open;
  }

  public setThermalOverlayMode(enabled: boolean): void {
    this.thermalOverlayMode = enabled;
  }

  public setCameraPreset(preset: CameraPreset): void {
    if (preset === "overview") {
      this.camera.position.set(0, -34, 66);
      this.controls.target.set(0, -2, 0);
    } else if (preset === "tourbillon") {
      this.camera.position.set(0, -22, 22);
      this.controls.target.set(0, -11, 0.5);
    } else if (preset === "barrel") {
      this.camera.position.set(0, -6, 28);
      this.controls.target.set(0, 9.5, 0.5);
    } else if (preset === "escapement") {
      this.camera.position.set(-6, -18, 13);
      this.controls.target.set(-1.5, -10.8, 0.2);
    } else if (preset === "side" || preset === "profile") {
      this.camera.position.set(46, -24, 18);
      this.controls.target.set(0, -2, 0);
    }
    this.controls.update();
  }

  private bindPointerEvents(): void {
    const pickWatchPlane = (ev: PointerEvent): [number, number, number] | null => {
      const rect = this.canvas.getBoundingClientRect();
      this.pointerNdc.x = ((ev.clientX - rect.left) / rect.width) * 2 - 1;
      this.pointerNdc.y = -((ev.clientY - rect.top) / rect.height) * 2 + 1;
      this.raycaster.setFromCamera(this.pointerNdc, this.camera);
      const hit = new THREE.Vector3();
      if (this.raycaster.ray.intersectPlane(this.interactionPlane, hit)) {
        const x = Math.max(-29, Math.min(29, hit.x));
        const y = Math.max(-29, Math.min(29, hit.y));
        const z = this.activeTool === "mic" ? 11.0 : 2.0;
        return [x, y, z];
      }
      return null;
    };

    this.canvas.addEventListener("pointerdown", (ev) => {
      if (this.activeTool === "inspect" || ev.button !== 0) return;
      this.isDraggingTool = true;
      const pt = pickWatchPlane(ev);
      if (pt) this.callbacks.onTargetMoved(this.activeTool, pt);
    });

    window.addEventListener("pointermove", (ev) => {
      if (!this.isDraggingTool || this.activeTool === "inspect") return;
      const pt = pickWatchPlane(ev);
      if (pt) this.callbacks.onTargetMoved(this.activeTool, pt);
    });

    window.addEventListener("pointerup", () => {
      this.isDraggingTool = false;
    });
  }

  public resize(): void {
    const parent = this.canvas.parentElement;
    if (!parent) return;
    const w = Math.max(320, parent.clientWidth);
    const h = Math.max(260, parent.clientHeight);
    this.renderer.setSize(w, h, false);
    this.composer.setSize(w, h);
    this.bloomPass.resolution.set(w, h);
    this.camera.aspect = w / h;
    this.camera.updateProjectionMatrix();
  }
}
