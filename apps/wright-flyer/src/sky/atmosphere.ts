// Atmospheric-sky mathematics for the Wright Flyer visual overhaul
// (sky foundation bead). PURE and deterministic: every color a shader
// or the fog/background consumes comes from here so it can be tested
// headless. No THREE.js types leak into this module — the builders in
// dressing3d/flyerScene are thin consumers.
// Repro: node --test test/atmosphere.test.ts
//
// World frame (terrainMesh.ts): x east, z south, y up. The sun
// direction points FROM the origin TOWARD the sun.
//
// Model: single-scatter-flavored analytic approximation. Direct-beam
// transmittance follows a Kasten-Young airmass with wavelength-dependent
// optical depth (Rayleigh ~ lambda^-4 plus mild maritime aerosol); the
// dome splits into a Rayleigh-dominated zenith against an
// aerosol/transmission-dominated horizon, with the Mie forward glow and
// the photosolar disc layered on top in the fragment stage. It is an
// approximation tuned for ONE afternoon — December 17, 1903, 10:35 a.m.,
// Kill Devil Hills NC (~36 N): LOW sun ~28 degrees over the SE horizon,
// hard winter light, pale haze off the Atlantic. Physically plausible,
// not physically solved.

/** Sun elevation above the horizon for the scenario instant [deg].
 * Dec 17 1903, ~10:35 a.m. local, Kill Devil Hills (~36 N): the sun
 * sits ~28 deg up in the southeast — LOW enough for warm transmission
 * colors, high enough that this is daylight, not sunset. */
export const SUN_ELEVATION_DEG = 28;

/** Linear-light RGB triple, channels 0..1 for ordinary colors, > 1
 * allowed for HDR (the sun disc). LINEAR, never sRGB: consumers feed
 * these straight into three.js color uniforms under the ACES pipeline. */
export type RGB = readonly [number, number, number];

/* ------------------------------ tuning ------------------------------ */

const DEG_RAD = Math.PI / 180;

/** Sea-level optical depth per unit airmass, per channel (linear RGB).
 * Ratio B/R ~= 3 tracks Rayleigh's ~lambda^-4 scaling; the elevated
 * floor over pure Rayleigh (~0.05 at blue) stands in for the mild
 * maritime aerosol the Wrights' accounts and photographs show as a
 * pale coastal haze. Tuned so the 28-degree sun transmits WARM-WHITE,
 * not sunset-red. */
const OPTICAL_DEPTH_RGB: readonly [number, number, number] = [0.1, 0.16, 0.3];

/** Relative air mass along the oblique path (Kasten-Young 1989): the
 * standard low-elevation correction; at 28 deg it gives m ~= 2.1 —
 * twice the vertical column, which is what makes a low sun warm. */
function airmass(elevDeg: number): number {
  const s = Math.sin(elevDeg * DEG_RAD);
  return 1 / (s + 0.15 * Math.pow(elevDeg + 3.885, -1.253));
}

/** Fraction of DIRECT sunlight surviving the path to the observer at
 * the scenario elevation, per channel: T = exp(-tau * m). */
function transmittance(): RGB {
  const m = airmass(SUN_ELEVATION_DEG);
  return [
    Math.exp(-OPTICAL_DEPTH_RGB[0]! * m),
    Math.exp(-OPTICAL_DEPTH_RGB[1]! * m),
    Math.exp(-OPTICAL_DEPTH_RGB[2]! * m),
  ];
}

/** Peak HDR radiance multiplier for the disc core. Deliberately far
 * above 1 so the ACESFilmic tonemap (renderer exposure ~1.12) has
 * headroom to roll the limb off and bloom pass can catch it. */
const DISC_HDR_PEAK = 60;

/** Solar photosphere color at ground level before atmospheric loss
 * (slightly warm white, ~5600 K). Multiplied by transmittance for the
 * observed disc. */
const SOLAR_WHITE: RGB = [1.0, 0.96, 0.9];

/* ---------------------------- validation ---------------------------- */

/** House-style domain refusal: sunDir must be a finite, unit-length
 * RGB (tolerance 1e-3 — shader-fed directions round-trip through
 * normalize()). Every exported function calls this first. */
function checkSunDir(sunDir: RGB): void {
  if (
    !Array.isArray(sunDir) ||
    sunDir.length !== 3 ||
    !sunDir.every((c) => typeof c === "number" && Number.isFinite(c))
  ) {
    throw new RangeError(`sun direction must be a finite [x, y, z] triple, got ${JSON.stringify(sunDir)}`);
  }
  const len = Math.hypot(sunDir[0], sunDir[1], sunDir[2]);
  if (Math.abs(len - 1) > 1e-3) {
    throw new RangeError(`sun direction must be normalized (|v| = 1 +/- 1e-3), got length ${len}`);
  }
}

const clamp01 = (v: number): number => Math.min(Math.max(v, 0), 1);

/** Actual elevation of the supplied direction [deg]; drives the
 * low-sun warmth modulation so callers who nudge the sun still get
 * plausible shifts instead of frozen scenario colors. */
function elevationDegOf(sunDir: RGB): number {
  return Math.asin(clamp01(sunDir[1])) / DEG_RAD;
}

/* ------------------------------ colors ------------------------------ */

/** Zenith dome color (LINEAR rgb): Rayleigh-dominated deep sky.
 * Scatter strength follows the optical-depth ratios (blue ~3x red),
 * then a multiple-scattering white floor lifts the red channel —
 * without it the zenith reads as saturated alpine blue, whereas
 * photographs of a clear December mid-Atlantic sky show a DESATURATED
 * blue-gray. Dimmed slightly as the sun drops below the scenario
 * elevation (less incident flux to scatter). */
export function zenithColor(sunDir: RGB): RGB {
  checkSunDir(sunDir);
  const tau = OPTICAL_DEPTH_RGB;
  const blue = tau[2]!;
  const high = clamp01(
    Math.sin(elevationDegOf(sunDir) * DEG_RAD) / Math.sin(SUN_ELEVATION_DEG * DEG_RAD),
  );
  const luma = 0.42; // zenith radiance anchor; exposure 1.12 lifts this to a readable sky
  const dim = 0.85 + 0.15 * high;
  return [
    luma * (0.75 * (tau[0]! / blue) + 0.25) * dim,
    luma * (0.75 * (tau[1]! / blue) + 0.25) * dim,
    luma * (0.75 * 1.0 + 0.25) * dim,
  ];
}

/** Horizon dome color (LINEAR rgb): warm haze band. Near the horizon
 * the Rayleigh path is optically thick — scattered blue leaves the
 * beam — so what reaches the eye is direct warm sunlight transmitted
 * through the column PLUS neutral-grey aerosol scatter. Warmth grows
 * as the sun descends (longer columns), which at 28 deg is already a
 * clearly warm pale band, strongest toward the sun azimuth (the Mie
 * glow layer in the fragment stage supplies that asymmetry). */
export function horizonColor(sunDir: RGB): RGB {
  checkSunDir(sunDir);
  const t = transmittance();
  const grey = 0.85; // aerosol scatter is nearly spectrally flat
  const aerosol = 0.35; // fraction of horizon radiance from grey aerosol vs transmitted beam
  const warm = clamp01((32 - elevationDegOf(sunDir)) / 32); // 0 at/above 32 deg -> 1 at horizon-level sun
  const luma = 0.95; // horizon band is the brightest sky element after the disc itself
  return [
    luma * (t[0]! * (1 - aerosol) + grey * aerosol) + 0.06 * warm,
    luma * (t[1]! * (1 - aerosol) + grey * aerosol),
    luma * (t[2]! * (1 - aerosol) + grey * aerosol) - 0.1 * warm,
  ];
}

/** Below-horizon dome band (LINEAR rgb): the thick near-ground aerosol
 * layer hanging over the dune sand and the winter Atlantic. Denser
 * than the sky horizon (aerosol density grows downward), so BRIGHTER,
 * and sandy-warm from the terrain it veils. Rendered as the
 * uGroundHaze uniform blended in over the 0..-0.06 dir.y band. */
export function groundHazeColor(sunDir: RGB): RGB {
  checkSunDir(sunDir);
  const t = transmittance();
  const luma = 1.04; // brighter than the sky horizon: looking DOWN the haze column
  const sand = [1.0, 0.94, 0.83] as const; // dune-sand reflectance tint (warm, blue-starved)
  return [
    luma * (t[0]! * 0.55 + sand[0] * 0.45),
    luma * (t[1]! * 0.55 + sand[1] * 0.45),
    luma * (t[2]! * 0.55 + sand[2] * 0.45),
  ];
}

/** Photosolar disc color (LINEAR, HDR > 1 at the core): transmitted
 * solar white at scenario elevation, scaled to DISC_HDR_PEAK so the
 * ACES tonemapper — not this module — decides how blown-out the core
 * reads. */
export function sunDiscColor(sunDir: RGB): RGB {
  checkSunDir(sunDir);
  const t = transmittance();
  return [
    SOLAR_WHITE[0] * t[0]! * DISC_HDR_PEAK,
    SOLAR_WHITE[1] * t[1]! * DISC_HDR_PEAK,
    SOLAR_WHITE[2] * t[2]! * DISC_HDR_PEAK,
  ];
}

/* --------------------------- fog + distance -------------------------- */

/** Exponential aerial-perspective scale [m]. Chosen to complement the
 * scene's linear THREE.Fog(0xd9dee2, 260, 2400): ~40% haze by 760 m
 * (mid-field dunes), ~80% by the 2400 m fog end, so far terrain melts
 * into the sky band instead of cutting off at a plane. */
const AERIAL_SCALE_M = 1500;

/** Linear fog/target color shared by fogColorHex and aerialPerspective:
 * the horizon band mixed toward the (brighter) ground haze, because a
 * horizontal sightline spends most of its length just ABOVE the
 * terrain looking through both layers. Warm-gray by construction —
 * the December coast is haze over sand, not alpine blue. */
function fogLinear(sunDir: RGB): RGB {
  const h = horizonColor(sunDir);
  const g = groundHazeColor(sunDir);
  const gW = 0.35;
  return [
    h[0] * (1 - gW) + g[0] * gW,
    h[1] * (1 - gW) + g[1] * gW,
    h[2] * (1 - gW) + g[2] * gW,
  ];
}

/** sRGB hex for THREE.Fog / scene.background. Encodes fogLinear()
 * linear->sRGB (the sRGB transfer function, IEC 61966-2-1) because the
 * scene sets Color.setHex, which interprets hex as sRGB. */
export function fogColorHex(sunDir: RGB): number {
  checkSunDir(sunDir);
  const enc = (c: number): number =>
    Math.round(
      255 *
        (c <= 0.0031308 ? 12.92 * c : 1.055 * Math.pow(c, 1 / 2.4) - 0.055),
    );
  const f = fogLinear(sunDir);
  return (enc(clamp01(f[0])) << 16) | (enc(clamp01(f[1])) << 8) | enc(clamp01(f[2]));
}

/** Aerial perspective for terrain/material shading: how much of a
 * surface at distM dissolves into the sky. fogFactor01 in [0, 1) is an
 * exponential in distance; colorShift is the LINEAR fog color the
 * surface lerps toward (lerp(surface, colorShift, fogFactor01)).
 * Optional — the global THREE.Fog already covers most cases; this is
 * for per-material custom shaders. Refuses non-finite or negative
 * distances. */
export function aerialPerspective(
  sunDir: RGB,
  distM: number,
): { colorShift: RGB; fogFactor01: number } {
  checkSunDir(sunDir);
  if (!Number.isFinite(distM) || distM < 0) {
    throw new RangeError(`distance must be finite and >= 0, got ${distM}`);
  }
  return { colorShift: fogLinear(sunDir), fogFactor01: 1 - Math.exp(-distM / AERIAL_SCALE_M) };
}

/* ------------------------------- GLSL -------------------------------- */

/** Sky-dome VERTEX stage. The dome mesh is a large inward-facing
 * sphere centred on the camera, so the world-space vertex position IS
 * (an offset along) the view direction; the fragment stage subtracts
 * cameraPosition and normalizes. */
export const SKY_DOME_VERT_GLSL: string = /* glsl */ `
// Wright Flyer sky dome — vertex stage (companion to src/sky/atmosphere.ts).
varying vec3 vWorldPos;

void main() {
  vWorldPos = (modelMatrix * vec4(position, 1.0)).xyz;
  gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0);
}
`;

/** Sky-dome FRAGMENT stage. Uniforms are fed verbatim from the pure
 * functions above. Outputs PREMULTIPLIED-LINEAR color (alpha 1) —
 * tone mapping happens downstream in the renderer (ACESFilmic,
 * exposure ~1.12), never here. */
export const SKY_DOME_FRAG_GLSL: string = /* glsl */ `
// Wright Flyer sky dome — fragment stage (companion to src/sky/atmosphere.ts).
precision highp float;

uniform vec3 uSunDirection; // normalized, world frame (x east, z south, y up)
uniform vec3 uZenith;       // linear rgb, zenithColor()
uniform vec3 uHorizon;      // linear rgb, horizonColor()
uniform vec3 uGroundHaze;   // linear rgb, groundHazeColor()
uniform vec3 uSunDisc;      // linear HDR rgb, sunDiscColor()
uniform float uDiscCos;     // cos of the sun's angular radius (~0.27 deg -> 0.999989)
uniform float uGlowStrength; // Mie forward-glow gain (scene-tuned)

varying vec3 vWorldPos;

// Zenith-gradient exponent: pow(max(dir.y, 0), K_ZENITH). ~1.8 keeps the
// deep colour overhead while letting the bright band hug the horizon —
// flatter exponents smear the haze too high, steeper ones pinch it.
const float K_ZENITH = 1.8;

// Mie phase-function falloff for coastal haze at visible wavelengths.
// Small particles (r ~ wavelength) forward-scatter strongly; n = 6 gives
// a tight glare cone around the sun without a full HG implementation.
const float MIE_N = 6.0;

// Forward-scattered haze tint: the 28-degree sun's transmitted color
// (matches transmittance() ratios) — the glow is warm, never white.
const vec3 GLOW_TINT = vec3(0.81, 0.71, 0.53);

// Disc limb softness: cos(r + e) - cos(r) for e ~ 0.35 deg. Kept as a
// DELTA on uDiscCos rather than a second angle so the pair survives any
// future change of disc radius.
const float DISC_EDGE_DELTA = 0.00006;

void main() {
  vec3 dir = normalize(vWorldPos - cameraPosition);

  // Rayleigh split: horizon band blending to the zenith colour overhead.
  float up = max(dir.y, 0.0);
  vec3 sky = mix(uHorizon, uZenith, pow(up, K_ZENITH));

  // Mie forward glow: brightest looking straight at the sun, decaying
  // with angle — the hazy bright halo a low coastal sun wears.
  float mu = dot(dir, uSunDirection);
  float glow = pow(max(mu, 0.0), MIE_N) * uGlowStrength;
  vec3 col = sky + glow * GLOW_TINT;

  // Photosolar disc: 1 inside the angular radius, rolling off over the
  // soft limb (smoothstep(cos(r + e), cos(r), dot)).
  float disc = smoothstep(uDiscCos - DISC_EDGE_DELTA, uDiscCos, mu);
  col += disc * uSunDisc;

  // Below the horizon the dome becomes the ground-haze band, blended in
  // over dir.y 0 .. -0.06 (about 3.4 deg at dome scale — reads as the
  // thickness of the near-ground aerosol layer, not a hard cut).
  float below = smoothstep(0.0, 0.06, -dir.y);
  col = mix(col, uGroundHaze, below);

  gl_FragColor = vec4(col, 1.0); // premultiplied-linear
}
`;
