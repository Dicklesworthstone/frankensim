// E2.2 scene: the parametric airframe in the dawn scene, driven by the
// scripted demo (bead wf-root-guzez.3.2). Implements the FlyerRenderer
// seam so main.ts swaps one factory call.

import * as THREE from "three";
import { buildWrightFlyerAirframe } from "./airframe/parametricAirframe.ts";
import { driveScripted } from "./airframe/applyPose.ts";
import { buildTerrainArrays } from "./terrainMesh.ts";
import { type CameraPreset, PRESET_KEYS, cameraFor } from "./camera.ts";
import { NEUTRAL, keysFrom, stepCommand } from "./input.ts";
import { hudLines } from "./hud.ts";
import { computePose } from "./airframe/pose.ts";
import { applyPose, scriptedState } from "./airframe/applyPose.ts";
import kdhGrid from "../../../data/wright-flyer/terrain/kill-devil-hills-17x17-v1.json";
import type { FlyerRenderer } from "./renderer.ts";
import type { SimClient } from "./sim/simClient.ts";
import {
  advanceProp,
  controlStateFrom,
  hudInputsFrom,
  phaseBanner,
  worldTransformFrom,
  type SimDriveState,
} from "./sim/snapshotView.ts";

export function createFlyerSceneRenderer(
  container: HTMLElement,
  simClient?: SimClient,
): FlyerRenderer {
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
  // E2.4: input, cameras, HUD. The pilot takes over from the script the
  // first time a control key goes down.
  const down = new Set<string>();
  let preset: CameraPreset = "free";
  let manual = false;
  let command = NEUTRAL;
  const hud = document.createElement("div");
  hud.style.cssText =
    "position:fixed;left:12px;bottom:12px;font:12px/1.5 monospace;color:#f5efe0;" +
    "background:rgba(20,24,30,.55);padding:8px 10px;border-radius:6px;white-space:pre";
  container.appendChild(hud);
  const onKey = (e: KeyboardEvent, isDown: boolean): void => {
    if (isDown && PRESET_KEYS[e.code]) {
      preset = PRESET_KEYS[e.code]!;
      return;
    }
    if (isDown) manual = manual || ["ArrowUp", "ArrowDown", "ArrowLeft", "ArrowRight", "KeyW", "KeyA", "KeyS", "KeyD", "Space"].includes(e.code);
    if (isDown) down.add(e.code);
    else down.delete(e.code);
  };
  const keydown = (e: KeyboardEvent): void => onKey(e, true);
  const keyup = (e: KeyboardEvent): void => onKey(e, false);
  window.addEventListener("keydown", keydown);
  window.addEventListener("keyup", keyup);
  // E5.2b: the sim drive (REAL engine state) supersedes both the
  // script and manual pose play when a SimClient is attached.
  let drive: SimDriveState = { propAngleRad: 0 };
  return {
    render(dtS: number): void {
      elapsedS += dtS;
      const snap = simClient?.sample(performance.now()) ?? null;
      if (snap !== null) {
        drive = advanceProp(drive, snap, dtS);
        const pose = computePose(controlStateFrom(snap, drive));
        applyPose(airframe, pose);
        const world = worldTransformFrom(snap, launch);
        airframe.group.position.set(world.position[0], world.position[1] + 1.2, world.position[2]);
        airframe.group.rotation.set(0, 0, world.pitchRad);
        const cam = cameraFor(preset, elapsedS, launch, [
          world.position[0],
          world.position[1] + 1.2,
          world.position[2],
        ]);
        camera.position.set(cam.pos[0], cam.pos[1], cam.pos[2]);
        camera.lookAt(cam.look[0], cam.look[1], cam.look[2]);
        const hudIn = hudInputsFrom(snap);
        const lines = hudLines({
          airspeedMps: hudIn.airspeedMps,
          elapsedS: hudIn.elapsedS,
          engineRpm: hudIn.engineRpm,
          camera: `${preset} (sim)`,
          pose,
        });
        lines.push(`phase ${hudIn.phase}  h ${snap.hM.toFixed(1)} m  x ${snap.xM.toFixed(1)} m`);
        const banner = phaseBanner(snap, simClient?.envelopeRefusalCode());
        if (banner !== null) {
          lines.push(banner);
        }
        hud.textContent = lines.join("\n");
        renderer.render(scene, camera);
        return;
      }
      let pose;
      if (manual) {
        command = stepCommand(command, keysFrom(down), dtS);
        pose = computePose({
          canardDeg: command.canard * 30,
          warpDeg: command.warp * 8.5,
          rudderDeg: 0,
          coupled: true,
          propAngleRad: scriptedState(elapsedS).propAngleRad,
        });
        applyPose(airframe, pose);
      } else {
        pose = driveScripted(airframe, elapsedS);
      }
      airframe.group.rotation.z = 0.02 * Math.sin(elapsedS * 0.8); // idle sway
      const ac: [number, number, number] = [
        airframe.group.position.x, airframe.group.position.y, airframe.group.position.z,
      ];
      const cam = cameraFor(preset, elapsedS, launch, ac);
      camera.position.set(cam.pos[0], cam.pos[1], cam.pos[2]);
      camera.lookAt(cam.look[0], cam.look[1], cam.look[2]);
      hud.textContent = hudLines({
        airspeedMps: 10.73,
        elapsedS,
        engineRpm: Math.min(1025, elapsedS * 180),
        camera: manual ? `${preset} (manual)` : preset,
        pose,
      }).join("\n");
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
      window.removeEventListener("keydown", keydown);
      window.removeEventListener("keyup", keyup);
      container.removeChild(hud);
      renderer.dispose();
      container.removeChild(renderer.domElement);
    },
  };
}
