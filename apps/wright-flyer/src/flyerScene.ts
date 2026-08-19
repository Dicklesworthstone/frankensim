// E2.2 scene: the parametric airframe in the dawn scene, driven by the
// scripted demo (bead wf-root-guzez.3.2). Implements the FlyerRenderer
// seam so main.ts swaps one factory call.

import * as THREE from "three";
import { buildWrightFlyerAirframe } from "./airframe/parametricAirframe.ts";
import { driveScripted } from "./airframe/applyPose.ts";
import { arrivalCamera, buildTerrainArrays } from "./terrainMesh.ts";
import kdhGrid from "../../../data/wright-flyer/terrain/kill-devil-hills-17x17-v1.json";
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
  // The REAL Kill Devil Hills tile (E1.3 3DEP grid): heightfield + splat.
  const terrain = buildTerrainArrays(kdhGrid, 96);
  const tGeo = new THREE.BufferGeometry();
  tGeo.setAttribute("position", new THREE.BufferAttribute(terrain.positions, 3));
  tGeo.setAttribute("color", new THREE.BufferAttribute(terrain.colors, 3));
  tGeo.setIndex(new THREE.BufferAttribute(terrain.indices, 1));
  tGeo.computeVertexNormals();
  scene.add(new THREE.Mesh(tGeo, new THREE.MeshStandardMaterial({ vertexColors: true, roughness: 1 })));
  const airframe = buildWrightFlyerAirframe();
  const launch = terrain.launch;
  airframe.group.position.set(launch[0], launch[1] + 1.2, launch[2]);
  scene.add(airframe.group);
  let elapsedS = 0;
  return {
    render(dtS: number): void {
      elapsedS += dtS;
      const pose = driveScripted(airframe, elapsedS);
      airframe.group.rotation.z = 0.02 * Math.sin(elapsedS * 0.8); // idle sway
      // The arrival shot owns the camera for the first 14 s, then holds.
      const cam = arrivalCamera(elapsedS, launch);
      camera.position.set(cam.pos[0], cam.pos[1], cam.pos[2]);
      camera.lookAt(cam.look[0], cam.look[1], cam.look[2]);
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
