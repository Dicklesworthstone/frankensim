// E5.6 QoS-gated post-processing chain: RenderPass -> UnrealBloomPass ->
// vignette -> SMAAPass -> OutputPass, tiered to mirror qos.ts fieldThrottle
// (0 = full chain, 1 = bloom only, 2 = direct renderer.render bypass).
// LAW: this module owns PRESENTATION polish only. It never touches scene
// contents, camera state, or anything physics-bearing — same gate as the
// governor itself.
//
// DOUBLE-TONEMAP CONTRACT (verified against three r172 source, not memory):
// - WebGLPrograms applies renderer.toneMapping to material shaders ONLY when
//   rendering to the default framebuffer (`currentRenderTarget === null`);
//   EffectComposer's internal targets are HalfFloatType, so RenderPass output
//   stays linear-HDR and OutputPass performs the ONE ACESFilmic pass plus the
//   sRGB transfer at the end of the chain.
// - Tier 2 bypasses the composer entirely; renderer.render() hits the default
//   framebuffer directly, so materials tonemap once there too. Both paths are
//   single-tonemap by construction — provided the CALLER KEEPS, on the
//   renderer it hands us:
//     * toneMapping = THREE.ACESFilmicToneMapping, toneMappingExposure = 1.12
//       (OutputPass re-reads both from the renderer EVERY frame — never zero
//       them out "for the composer path");
//     * logarithmicDepthBuffer = true and pixelRatio capped at 2
//       (forwarded here via setPixelRatio, which the composer mirrors into
//       its target sizing);
//     * antialias = true (MSAA covers the tier-2 direct path; the composer
//       path gets edge quality from SMAAPass instead, since MSAA does not
//       survive an offscreen half-float target).
//
// SELF-CHECK (browser verification plan, no unit test by design):
// 1. Boot at tier 0: sun disc / wet-sand glints bloom softly (strength 0.32),
//    corners visibly darker than center but never crushed, edges crisp under
//    SMAA; screenshot should match direct render except at high-contrast
//    edges and corners.
// 2. Flip to tier 1: vignette and SMAA drop out (pass.enabled=false), bloom
//    persists — corners brighten back to flat, aliasing returns on grid
//    lines. No resize glitch on transition.
// 3. Flip to tier 2: identical exposure/color to tier 1 minus bloom — proves
//    no double-tonemap (a second ACES pass would look washed and desaturated
//    instantly).
// 4. Resize window across all tiers: no stretched bloom halos (bloom mip
//    chain resizes via composer.setSize -> pass.setSize) and no SMAA blur.
// 5. Frame-time probe at tier 2 vs tier 0 confirms the Critical tier is
//    materially cheaper (composer disabled, one straight draw).

import type { Camera, Scene } from "three";
import { Vector2, WebGLRenderer } from "three";
import { EffectComposer } from "three/examples/jsm/postprocessing/EffectComposer.js";
import { OutputPass } from "three/examples/jsm/postprocessing/OutputPass.js";
import { RenderPass } from "three/examples/jsm/postprocessing/RenderPass.js";
import { ShaderPass } from "three/examples/jsm/postprocessing/ShaderPass.js";
import { SMAAPass } from "three/examples/jsm/postprocessing/SMAAPass.js";
import { UnrealBloomPass } from "three/examples/jsm/postprocessing/UnrealBloomPass.js";

/**
 * Post tier, aligned with `PresentationProfile.fieldThrottle` (qos.ts):
 * 0 = full chain (bloom + vignette + SMAA), 1 = bloom only, 2 = composer
 * bypassed entirely (Critical tier pays one straight draw, nothing else).
 */
export type PostTier = 0 | 1 | 2;

export interface PostChain {
  /** The SAME canvas the renderer owns; the composer renders to screen. */
  readonly domElement: HTMLCanvasElement;
  /**
   * Apply a presentation tier. Idempotent: re-applying the current tier is a
   * no-op, and transitions only flip `pass.enabled` flags (no rebuild, no
   * reallocation — a mid-flight QoS demotion must not stutter).
   */
  setTier(tier: PostTier): void;
  /** CSS-pixel size; forwarded to the composer, which fans out to every
   * pass (bloom mips, SMAA render targets) with device pixels applied. */
  setSize(width: number, height: number): void;
  /** Forwarded to the composer; it multiplies into all target sizes. Must
   * match the caller's `renderer.setPixelRatio` cap of 2. */
  setPixelRatio(pr: number): void;
  /** Advance one presentation frame; `dtS` is wall-clock seconds. Passing
   * it explicitly keeps the composer off its internal delta clock. */
  render(dtS: number): void;
  /** Release every pass-owned resource AND the composer's ping-pong targets
   * (EffectComposer.dispose only handles its own two targets and copyPass —
   * member passes are OUR responsibility, verified in r172 source). The
   * canvas stays with the renderer/caller. */
  dispose(): void;
}

/** Vignette shader: darkens toward corners, max 0.22 attenuation exactly at
 * the extreme corner, smoothstep onset past 60% of the way out. Static
 * uniforms — nothing is written per frame, so no per-frame allocation. */
const VIGNETTE_SHADER = {
  uniforms: {
    tDiffuse: { value: null },
  },
  vertexShader: /* glsl */ `
    varying vec2 vUv;
    void main() {
      vUv = uv;
      gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0);
    }
  `,
  fragmentShader: /* glsl */ `
    uniform sampler2D tDiffuse;
    varying vec2 vUv;
    void main() {
      vec4 color = texture2D(tDiffuse, vUv);
      // Scale so the extreme corner lands at r = 1.0 exactly (half-diagonal
      // of unit UV space is sqrt(0.5)), making the max shade exact.
      float r = length(vUv - 0.5) * 1.41421356;
      float shade = smoothstep(0.6, 1.0, r);
      color.rgb *= mix(1.0, 1.0 - 0.22, shade);
      gl_FragColor = color;
    }
  `,
} as const;

/** r172's SMAAPass constructor takes (width, height) at RUNTIME (verified in
 * examples/jsm/postprocessing/SMAAPass.js) but @types/three declares zero
 * args. Bridge the drift with one honest local signature. */
type SMAAPassCtor = new (width: number, height: number) => SMAAPass;

const BLOOM_STRENGTH = 0.32;
const BLOOM_RADIUS = 0.55;
const BLOOM_THRESHOLD = 0.85;

export function createPostChain(
  renderer: WebGLRenderer,
  scene: Scene,
  camera: Camera,
): PostChain {
  const composer = new EffectComposer(renderer);

  const renderPass = new RenderPass(scene, camera);

  // Threshold 0.85 keeps bloom to HDR sun disc + specular glints; everything
  // below that line passes through untouched. Resolution is a seed only —
  // composer.setSize drives the mip chain thereafter.
  const bloomPass = new UnrealBloomPass(
    new Vector2(256, 256),
    BLOOM_STRENGTH,
    BLOOM_RADIUS,
    BLOOM_THRESHOLD,
  );

  const vignettePass = new ShaderPass(VIGNETTE_SHADER);

  const SMAAPassWH = SMAAPass as unknown as SMAAPassCtor;
  const size = renderer.getSize(new Vector2());
  const smaaPass = new SMAAPassWH(size.x, size.y);

  // MUST stay last: it owns the single tonemapping + sRGB conversion (see
  // header contract).
  const outputPass = new OutputPass();

  composer.addPass(renderPass);
  composer.addPass(bloomPass);
  composer.addPass(vignettePass);
  composer.addPass(smaaPass);
  composer.addPass(outputPass);

  let tier: PostTier = 0;
  let disposed = false;

  return {
    domElement: renderer.domElement,

    setTier(next: PostTier): void {
      if (disposed || next === tier) return;
      tier = next;
      // Tier 2 leaves every pass enabled=false AND routes render() around
      // the composer; belt and suspenders cost nothing.
      bloomPass.enabled = tier !== 2;
      vignettePass.enabled = tier === 0;
      smaaPass.enabled = tier === 0;
    },

    setSize(width: number, height: number): void {
      if (disposed) return;
      composer.setSize(width, height);
    },

    setPixelRatio(pr: number): void {
      if (disposed) return;
      // Composer.setPixelRatio internally re-runs setSize, so targets and
      // pass resolutions track the caller's cap in one call.
      composer.setPixelRatio(pr);
    },

    render(dtS: number): void {
      if (disposed) return;
      if (tier === 2) {
        // Critical tier: bypass the composer entirely. Materials tonemap ONCE
        // here (default framebuffer), matching the OutputPass path visually.
        renderer.render(scene, camera);
        return;
      }
      composer.render(dtS);
    },

    dispose(): void {
      if (disposed) return;
      disposed = true;
      // Passes first (they own materials/targets the composer knows nothing
      // about), then the composer's own ping-pong buffers and copy pass.
      renderPass.dispose();
      bloomPass.dispose();
      vignettePass.dispose();
      smaaPass.dispose();
      outputPass.dispose();
      composer.dispose();
    },
  };
}
