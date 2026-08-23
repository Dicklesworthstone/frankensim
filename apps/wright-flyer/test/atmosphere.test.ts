// Atmospheric-sky battery: domain refusals on the sun direction, the
// Rayleigh/Mie color RELATIONSHIPS (zenith cooler+bluer than the warm
// horizon; glow peaked at the sun; disc zero off-axis, saturated
// on-axis), fog hex sanity, aerial-perspective monotonicity, and GLSL
// string integrity. Colors are checked by RELATIONSHIP, not magic
// pixels, so tuning constants can move without rewriting the battery.
// Repro: node --test test/atmosphere.test.ts

import assert from "node:assert/strict";
import { test } from "node:test";

import {
  SKY_DOME_FRAG_GLSL,
  SKY_DOME_VERT_GLSL,
  SUN_ELEVATION_DEG,
  aerialPerspective,
  fogColorHex,
  groundHazeColor,
  horizonColor,
  sunDiscColor,
  zenithColor,
} from "../src/sky/atmosphere.ts";
import type { RGB } from "../src/sky/atmosphere.ts";

function jlog(kase: string, payload: string): void {
  console.log(`{"suite":"wf-app-atmosphere","case":"${kase}",${payload}}`);
}

/** Scenario sun: ~28 deg up over the SE horizon (x east, z south ->
 * SE azimuth is +x,-z; y = sin(28 deg)). */
const SCENARIO_SUN: RGB = [
  Math.cos(SUN_ELEVATION_DEG * (Math.PI / 180)) * Math.SQRT1_2,
  Math.sin(SUN_ELEVATION_DEG * (Math.PI / 180)),
  -Math.cos(SUN_ELEVATION_DEG * (Math.PI / 180)) * Math.SQRT1_2,
];

const luma = (c: RGB): number => 0.2126 * c[0] + 0.7152 * c[1] + 0.0722 * c[2];

test("every exported function refuses a non-finite or non-normalized sunDir", () => {
  const bads: unknown[] = [
    [Number.NaN, 0, 1],
    [0, Infinity, 0],
    [0, 0, undefined],
    "up",
    [0, 0], // wrong arity
    [0, 0, 2], // |v| = 2
    [0.6, 0.6, 0.6], // |v| ~ 1.039, outside 1e-3 tolerance
    [0, 0, 0],
  ];
  const fns = [zenithColor, horizonColor, groundHazeColor, sunDiscColor, fogColorHex];
  for (const bad of bads) {
    for (const f of fns) {
      assert.throws(() => f(bad as RGB), RangeError, `${f.name} must refuse ${JSON.stringify(bad)}`);
    }
    assert.throws(() => aerialPerspective(bad as RGB, 100), RangeError);
  }
  // the scenario direction itself admits everywhere
  assert.ok(Number.isFinite(fogColorHex(SCENARIO_SUN)));
  jlog("domain-refusals", `"refused":${bads.length * (fns.length + 1)}`);
});

test("determinism: same sunDir -> identical colors, every call", () => {
  for (const f of [zenithColor, horizonColor, groundHazeColor, sunDiscColor]) {
    assert.deepEqual(f(SCENARIO_SUN), f(SCENARIO_SUN), f.name);
  }
  assert.equal(fogColorHex(SCENARIO_SUN), fogColorHex(SCENARIO_SUN));
});

test("zenith is cooler, bluer, and darker than the warm horizon", () => {
  const z = zenithColor(SCENARIO_SUN);
  const h = horizonColor(SCENARIO_SUN);
  // bluer: blue-minus-red gap strongly positive overhead...
  assert.ok(z[2] - z[0] > 0.08, `zenith b-r ${z[2] - z[0]} should be clearly positive`);
  // ...and strictly more than at the horizon, which is WARM (r > b)
  assert.ok(h[0] > h[2], `horizon must be warm, got r=${h[0]} b=${h[2]}`);
  assert.ok(z[2] - z[0] > h[2] - h[0], "zenith bluer than horizon");
  // cooler/darker: zenith luminance well under the bright haze band
  assert.ok(luma(z) < luma(h) * 0.8, `zenith luma ${luma(z)} < 0.8x horizon ${luma(h)}`);
  // desaturated blue-gray, not alpine blue: red channel survives overhead
  assert.ok(z[0] > z[2] * 0.3, `zenith too saturated: r=${z[0]} vs b=${z[2]}`);
  jlog("zenith-vs-horizon", `"zenith":[${z}],"horizon":[${h}]`);
});

test("ground haze band is brighter than the sky horizon and stays warm", () => {
  const g = groundHazeColor(SCENARIO_SUN);
  const h = horizonColor(SCENARIO_SUN);
  assert.ok(luma(g) > luma(h), `haze luma ${luma(g)} must exceed horizon ${luma(h)}`);
  assert.ok(g[0] > g[2], "haze sandy-warm (r > b)");
  for (const c of g) {
    assert.ok(Number.isFinite(c) && c >= 0 && c <= 1.5, `channel sane: ${c}`);
  }
});

test("sun disc is HDR (>1) and warm-white", () => {
  const d = sunDiscColor(SCENARIO_SUN);
  for (const c of d) {
    assert.ok(c > 1, `disc channel ${c} must be HDR`);
    assert.ok(Number.isFinite(c));
  }
  assert.ok(d[0] > d[2], "low-sun disc transmits warm");
  // not sunset-red: the core keeps most of its green/blue
  assert.ok(d[1] > d[0] * 0.6 && d[2] > d[0] * 0.4, "disc reads warm-white, not red");
});

test("fog hex is finite, in range, and warm-gray", () => {
  const hex = fogColorHex(SCENARIO_SUN);
  assert.ok(Number.isInteger(hex) && hex >= 0 && hex <= 0xffffff, `hex ${hex}`);
  const r = (hex >> 16) & 0xff;
  const g = (hex >> 8) & 0xff;
  const b = hex & 0xff;
  // warm-gray: light overall, red-leaning over blue, channels close
  assert.ok(r >= b, `warm: r=${r} >= b=${b}`);
  assert.ok(Math.abs(r - b) < 48, "gray-ish, not orange paint");
  assert.ok(Math.min(r, g, b) > 128, "pale December haze, not dusk");
  jlog("fog-hex", `"hex":${hex.toString(16).padStart(6, "0")}`);
});

test("aerial perspective: zero haze at the camera, monotonic, bounded", () => {
  assert.throws(() => aerialPerspective(SCENARIO_SUN, -1), RangeError);
  assert.throws(() => aerialPerspective(SCENARIO_SUN, Number.NaN), RangeError);
  assert.throws(() => aerialPerspective(SCENARIO_SUN, Infinity), RangeError);
  let prev = -1;
  for (let d = 0; d <= 3000; d += 137) {
    const { colorShift, fogFactor01 } = aerialPerspective(SCENARIO_SUN, d);
    assert.ok(fogFactor01 >= prev && fogFactor01 < 1, `monotonic bounded at ${d}: ${fogFactor01}`);
    assert.ok(Number.isFinite(colorShift[0]! + colorShift[1]! + colorShift[2]!));
    prev = fogFactor01;
  }
  assert.ok(aerialPerspective(SCENARIO_SUN, 0).fogFactor01 === 0, "no haze at distance 0");
});

/* ------------------- shader-math property mirrors -------------------- */
/* The glow and disc terms live in SKY_DOME_FRAG_GLSL; the fragment
 * string is static text, so the properties are verified against this
 * JS mirror using the SAME formulas the shader source contains (the
 * substring assertions below fail loudly if shader and mirror drift). */

const MIE_N = 6;
const DISC_EDGE_DELTA = 0.00006;
/** cos of ~5 degrees off-axis: comfortably outside disc + limb. */
const WELL_OFF_AXIS_COS = Math.cos((5 * Math.PI) / 180);

test("Mie glow is strongest along the sun direction and decays off-axis", () => {
  const glowAt = (mu: number): number => Math.pow(Math.max(mu, 0), MIE_N);
  // strictly approaching the sun: mu rises to exactly 1 on-axis
  const mus = [180, 120, 90, 45, 10, 5, 2, 1, 0].map((offDeg) => Math.cos((offDeg * Math.PI) / 180));
  let prev = -1;
  for (const mu of mus) {
    const g = glowAt(mu);
    assert.ok(g >= prev && g <= 1, `glow monotonic toward sun at mu=${mu}`);
    prev = g;
  }
  assert.equal(glowAt(-0.5), 0, "anti-solar glow exactly zero");
  assert.equal(glowAt(1), 1, "glow saturates along sunDir");
  // substring sync check: the shader really applies uGlowStrength this way
  assert.match(SKY_DOME_FRAG_GLSL, /float mu = dot\(dir, uSunDirection\)/);
  assert.match(
    SKY_DOME_FRAG_GLSL,
    /pow\(max\(mu, 0\.0\), MIE_N\) \* uGlowStrength/,
  );
});

test("disc term is zero away from the sun and ~1 at the core", () => {
  const discAt = (mu: number, discCos: number): number => {
    const t = (mu - (discCos - DISC_EDGE_DELTA)) / DISC_EDGE_DELTA;
    const x = Math.min(Math.max(t, 0), 1);
    return x * x * (3 - 2 * x); // GLSL smoothstep
  };
  const discCos = Math.cos((0.27 * Math.PI) / 180); // solar angular radius
  assert.equal(discAt(WELL_OFF_AXIS_COS, discCos), 0, "zero looking away from the sun");
  assert.equal(discAt(-0.5, discCos), 0, "zero anti-solar");
  assert.ok(discAt(1, discCos) > 0.999, "~1 at the very center");
  // monotone through the limb band
  const limbLo = discCos - DISC_EDGE_DELTA / 2;
  const mid = discAt(limbLo, discCos);
  assert.ok(mid > 0 && mid < 1, `limb midpoint inside (0,1): ${mid}`);
  assert.match(
    SKY_DOME_FRAG_GLSL,
    /smoothstep\(uDiscCos - DISC_EDGE_DELTA, uDiscCos, mu\)/,
  );
});

test("GLSL chunks carry the contracted uniforms, varyings, and structure", () => {
  const fragRequired = [
    "uSunDirection",
    "uZenith",
    "uHorizon",
    "uGroundHaze",
    "uSunDisc",
    "uDiscCos",
    "uGlowStrength",
    "vWorldPos", // consumes the vertex stage's varying
  ];
  for (const name of fragRequired) {
    assert.ok(SKY_DOME_FRAG_GLSL.includes(`uniform vec3 ${name}`) ||
      SKY_DOME_FRAG_GLSL.includes(`uniform float ${name}`) ||
      SKY_DOME_FRAG_GLSL.includes(name), `frag missing ${name}`);
  }
  assert.ok(SKY_DOME_VERT_GLSL.includes("varying vec3 vWorldPos"), "vert passes world varying");
  // gradient shape: mix(horizon, zenith, pow(max(dir.y, 0), K))
  // below-horizon blend over the 0..-0.06 band
  assert.match(SKY_DOME_FRAG_GLSL, /up = max\(dir\.y, 0\.0\)/);
  assert.match(SKY_DOME_FRAG_GLSL, /mix\(uHorizon, uZenith, pow\(up, K_ZENITH\)\)/);
  // premultiplied-linear output, no tonemapping in-shader
  assert.match(SKY_DOME_FRAG_GLSL, /vec4\(col, 1\.0\)/);
  assert.ok(!SKY_DOME_FRAG_GLSL.toLowerCase().includes("tonemap"), "tone mapping stays downstream");
  jlog("glsl-integrity", `"fragBytes":${SKY_DOME_FRAG_GLSL.length},"vertBytes":${SKY_DOME_VERT_GLSL.length}`);
});
