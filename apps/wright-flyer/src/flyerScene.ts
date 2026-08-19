// E2.2 scene: the parametric airframe in the dawn scene, driven by the
// scripted demo (bead wf-root-guzez.3.2). Implements the FlyerRenderer
// seam so main.ts swaps one factory call.

import * as THREE from "three";
import { buildWrightFlyerAirframe } from "./airframe/parametricAirframe.ts";
import { driveScripted } from "./airframe/applyPose.ts";
import type { FlyerRenderer } from "./renderer.ts";

export function createFlyerSceneRenderer(container: HTMLElement): FlyerRenderer {
  const renderer = new THREE.WebGLRenderer({ antialias: true });
  renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
  container.appendChild(renderer.domElement);
  const scene = new THREE.Scene();
  scene.background = new THREE.Color(0x8fa8c8); // overcast December morning
  scene.fog = new THREE.Fog(0x8fa8c8, 60, 240);
  const camera = new THREE.PerspectiveCamera(50, 1, 0.1, 500);
  camera.position.set(9, 3.2, 11);
  camera.lookAt(0, 1.2, 0);
  const sun = new THREE.DirectionalLight(0xfff1dd, 2.2);
  sun.position.set(-40, 25, 20);
  scene.add(sun, new THREE.HemisphereLight(0xbcd0e8, 0x8a7d64, 0.9));
  const sand = new THREE.Mesh(
    new THREE.PlaneGeometry(600, 600),
    new THREE.MeshStandardMaterial({ color: 0xb8a883, roughness: 1 }),
  );
  sand.rotation.x = -Math.PI / 2;
  scene.add(sand);
  const airframe = buildWrightFlyerAirframe();
  airframe.group.position.y = 1.2; // on the launch dolly height
  scene.add(airframe.group);
  let elapsedS = 0;
  return {
    render(dtS: number): void {
      elapsedS += dtS;
      const pose = driveScripted(airframe, elapsedS);
      airframe.group.rotation.z = 0.02 * Math.sin(elapsedS * 0.8); // idle sway
      if (pose.clamped) {
        console.warn(JSON.stringify({ suite: "wf-scene", event: "control-stop", t: elapsedS }));
      }
      renderer.render(scene, camera);
    },
    resize(width: number, height: number): void {
      renderer.setSize(width, height);
      camera.aspect = width / Math.max(1, height);
      camera.updateProjectionMatrix();
    },
    dispose(): void {
      renderer.dispose();
      container.removeChild(renderer.domElement);
    },
  };
}
