// Thin renderer interface (plan §8.1): the app talks to THIS surface, never
// to three.js classes directly outside this module, so a future
// WebGPURenderer/TSL migration (E7.6) is a module change, not a rewrite.
// Bead wf-root-guzez.1.2 (E0.2).

import {
  AmbientLight,
  Color,
  DirectionalLight,
  Fog,
  GridHelper,
  Mesh,
  MeshStandardMaterial,
  PerspectiveCamera,
  PlaneGeometry,
  Scene,
  WebGLRenderer,
} from "three";

import type { PresentationProfile } from "./qos.ts";

export interface FlyerRenderer {
  /** Advance one presentation frame; `dtS` is wall-clock presentation time. */
  render(dtS: number): void;
  /** One presentation frame WITHOUT the loop-level try/catch — the
   * scene splits render/renderFrame so a thrown frame surfaces loudly
   * and the rAF loop survives (presentation resilience). */
  renderFrame(dtS: number): void;
  resize(width: number, height: number): void;
  dispose(): void;
  /** E5.6: apply an ATOMIC presentation profile (QoS governor output —
   * presentation only; the physics tier is immutable by construction). */
  applyQuality?(profile: PresentationProfile): void;
}

/**
 * The E0.2 placeholder scene: a dawn-lit ground plane and grid standing in
 * for Kill Devil Hills until E2.3 lands the real site. Deliberately trivial —
 * its only jobs are proving the render loop, the resize path, and the
 * renderer seam.
 */
export function createPlaceholderRenderer(container: HTMLElement): FlyerRenderer {
  const renderer = new WebGLRenderer({ antialias: true });
  renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
  container.appendChild(renderer.domElement);

  const scene = new Scene();
  scene.background = new Color(0x10151d);
  scene.fog = new Fog(0x10151d, 60, 400);

  const camera = new PerspectiveCamera(55, 1, 0.1, 1000);
  camera.position.set(14, 6, 22);
  camera.lookAt(0, 1, 0);

  const sand = new Mesh(
    new PlaneGeometry(600, 600),
    new MeshStandardMaterial({ color: 0x8a7a5c, roughness: 1.0 }),
  );
  sand.rotation.x = -Math.PI / 2;
  scene.add(sand);
  scene.add(new GridHelper(120, 60, 0x39424f, 0x232a35));

  const sun = new DirectionalLight(0xffe8c4, 2.2);
  sun.position.set(-40, 25, 10);
  scene.add(sun, new AmbientLight(0x8fa3c0, 0.5));

  let elapsedS = 0;
  return {
    render(dtS: number): void {
      this.renderFrame(dtS);
    },
    renderFrame(dtS: number): void {
      elapsedS += dtS;
      // A slow dawn drift proves per-frame updates without faking physics.
      camera.position.x = 14 * Math.cos(elapsedS * 0.03);
      camera.position.z = 22 * Math.sin(elapsedS * 0.03) + 10;
      camera.lookAt(0, 1, 0);
      renderer.render(scene, camera);
    },
    resize(width: number, height: number): void {
      camera.aspect = width / height;
      camera.updateProjectionMatrix();
      renderer.setSize(width, height);
    },
    dispose(): void {
      renderer.dispose();
      renderer.domElement.remove();
    },
  };
}
