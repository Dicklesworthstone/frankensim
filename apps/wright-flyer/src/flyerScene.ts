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
import huffmanGrid from "../../../data/wright-flyer/terrain/huffman-prairie-17x17-v1.json";
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
import { ghostAt, type FlightRecording } from "./sim/replay.ts";
import { CaptionStream, formatCaption } from "./sim/captions.ts";
import { IDLE_INPUTS, phaseDisplay } from "./gauges.ts";
import { createHudDials, createPhaseBanner } from "./hudDials.ts";
import { heightAt } from "./terrainMesh.ts";
import { orvilleReachableX } from "./dressing.ts";
import { buildDressing, buildProneWilbur, sandTileMaterial } from "./dressing3d.ts";

export function createFlyerSceneRenderer(
  container: HTMLElement,
  simClient?: SimClient,
  ghost?: FlightRecording,
): FlyerRenderer {
  const renderer = new THREE.WebGLRenderer({ antialias: true });
  renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
  container.appendChild(renderer.domElement);
  const scene = new THREE.Scene();
  // December-morning haze: the sky dome carries the color now; the fog
  // starts far enough out that the DUNES are visible (the old 60..240 m
  // fog erased the ground the moment the machine climbed).
  scene.background = new THREE.Color(0xc3d0dd);
  scene.fog = new THREE.Fog(0xc3d0dd, 300, 2100);
  // Far plane must clear the sky dome (2600 m) and the Atlantic.
  const camera = new THREE.PerspectiveCamera(50, 1, 0.1, 5000);
  camera.position.set(9, 3.2, 11);
  camera.lookAt(0, 1.2, 0);
  const sun = new THREE.DirectionalLight(0xfff1dd, 2.2);
  sun.position.set(-40, 25, 20);
  scene.add(sun, new THREE.HemisphereLight(0xbcd0e8, 0x8a7d64, 0.9));
  // The REAL site tile (E1.3 3DEP grids): heightfield + splat.
  // E5.4: ?site=huffman loads the Huffman Prairie grid.
  const site = new URLSearchParams(window.location.search).get("site");
  const terrain = buildTerrainArrays(site === "huffman" ? huffmanGrid : kdhGrid, 192);
  const tGeo = new THREE.BufferGeometry();
  tGeo.setAttribute("position", new THREE.BufferAttribute(terrain.positions, 3));
  tGeo.setAttribute("color", new THREE.BufferAttribute(terrain.colors, 3));
  tGeo.setIndex(new THREE.BufferAttribute(terrain.indices, 1));
  tGeo.computeVertexNormals();
  scene.add(new THREE.Mesh(tGeo, sandTileMaterial()));
  const airframe = buildWrightFlyerAirframe();
  const launch = terrain.launch;
  airframe.group.position.set(launch[0], launch[1] + 1.2, launch[2]);
  scene.add(airframe.group);
  // The player's pilot: Wilbur prone on the lower wing at the cradle.
  const wilbur = buildProneWilbur();
  wilbur.position.set(-0.1, 0.12, 0.06);
  airframe.cradleGroup.add(wilbur);
  // The Kitty Hawk diorama (bead guzez.13): sky, clouds, outer sand +
  // Atlantic, the launch rail, the 1903 camp, Orville, and the gulls.
  const grid = site === "huffman" ? huffmanGrid : kdhGrid;
  const tileExtent = (grid.grid_n - 1) * grid.spacing_m;
  const half = tileExtent / 2;
  const dressing = buildDressing(
    launch,
    site === "huffman" ? 30.0 : 18.3,
    tileExtent,
    (xRel, zRel) =>
      heightAt(grid, launch[0] + xRel + half, -(launch[2] + zRel) + half) - launch[1],
  );
  scene.add(dressing.group);
  // Orville's release point is LATCHED the first frame the machine is
  // off the rail (pure pose math needs the release constants).
  let orvilleReleaseX: number | null = null;
  let orvilleReleaseT: number | null = null;
  let elapsedS = 0;
  // E2.4: input, cameras, HUD. The pilot takes over from the script the
  // first time a control key goes down.
  const down = new Set<string>();
  let preset: CameraPreset = "free";
  let manual = false;
  let command = NEUTRAL;
  // Telemetry text panel (the honest raw numbers; T toggles, kept ON
  // by default — the dials never replace the telemetry, they front it).
  const hud = document.createElement("div");
  hud.id = "wf-telemetry";
  hud.style.cssText =
    "position:fixed;right:12px;bottom:12px;font:11px/1.5 monospace;color:#f5efe0;" +
    "background:rgba(20,24,30,.55);padding:8px 10px;border-radius:6px;white-space:pre;z-index:5";
  container.appendChild(hud);
  // Game HUD: analog dial panel + phase banner (tested math in gauges.ts).
  const dials = createHudDials(container);
  const phaseEl = createPhaseBanner(container);
  const onKey = (e: KeyboardEvent, isDown: boolean): void => {
    if (isDown && e.code === "KeyT" && !e.repeat) {
      hud.style.display = hud.style.display === "none" ? "block" : "none";
      return;
    }
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
  // E5.3b-ii: the honest caption stream (every claim labeled).
  const captions = new CaptionStream();
  let lastCaptionTick = -1;
  // E5.2c: the replay ghost — the PREVIOUS run's recorded flight as a
  // translucent twin, driven tick-locked to the live run.
  let ghostFrame: ReturnType<typeof buildWrightFlyerAirframe> | null = null;
  let ghostDrive: SimDriveState = { propAngleRad: 0 };
  if (ghost !== undefined) {
    ghostFrame = buildWrightFlyerAirframe();
    ghostFrame.group.traverse((obj) => {
      const mesh = obj as THREE.Mesh;
      if (mesh.isMesh) {
        const src = mesh.material as THREE.Material;
        const mat = src.clone();
        mat.transparent = true;
        mat.opacity = 0.3;
        mat.depthWrite = false;
        mesh.material = mat;
      }
    });
    ghostFrame.group.position.set(launch[0], launch[1] + 1.2, launch[2]);
    scene.add(ghostFrame.group);
  }
  return {
    render(dtS: number): void {
      elapsedS += dtS;
      const snap = simClient?.sample(performance.now()) ?? null;
      if (snap !== null) {
        if (ghostFrame !== null && ghost !== undefined) {
          const g = ghostAt(ghost, snap.tick);
          if (g !== null) {
            ghostDrive = advanceProp(ghostDrive, g, dtS);
            applyPose(ghostFrame, computePose(controlStateFrom(g, ghostDrive)));
            const gw = worldTransformFrom(g, launch);
            ghostFrame.group.position.set(gw.position[0], gw.position[1] + 1.2, gw.position[2]);
            ghostFrame.group.rotation.set(0, 0, gw.pitchRad);
          }
        }
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
        if (snap.assistActive) {
          lines.push("ASSIST ACTIVE (bounded authority 0.3 of canard stop — model aid, not history)");
        }
        const banner = phaseBanner(snap, simClient?.envelopeRefusalCode());
        if (banner !== null) {
          lines.push(banner);
        }
        // Analog panel + phase card (dials are fronting, not replacing,
        // the telemetry lines above).
        dials.update({
          airspeedMps: hudIn.airspeedMps,
          engineRpm: hudIn.engineRpm,
          elapsedS: hudIn.elapsedS,
          hM: snap.hM,
          thetaRad: snap.thetaRad,
          dcRad: snap.dcRad,
          warpRad: snap.warpRad,
        });
        const disp = phaseDisplay(snap.phase, banner);
        phaseEl.set(disp.text, disp.tone);
        // Diorama: Orville chases while the machine is ON the rail,
        // lets go at the first off-rail frame (latched), then watches.
        const onRail = snap.phase === "on-rail";
        if (!onRail && orvilleReleaseX === null) {
          orvilleReleaseX = orvilleReachableX(hudIn.elapsedS, snap.xM);
          orvilleReleaseT = hudIn.elapsedS;
        }
        dressing.animate(hudIn.elapsedS, {
          onRail,
          aircraftX: snap.xM,
          releaseX: orvilleReleaseX,
          releaseT: orvilleReleaseT,
        });
        // Labeled ride-along captions (latest two).
        if (snap.tick > lastCaptionTick) {
          lastCaptionTick = snap.tick;
          captions.feed(snap);
        }
        for (const c of captions.upTo(snap.tick).slice(-2)) {
          lines.push(formatCaption(c));
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
      // Attract-mode dials: scripted demo numbers where they exist,
      // rest position elsewhere (no fabricated altitude/pitch).
      dials.update({
        ...IDLE_INPUTS,
        airspeedMps: 10.73,
        elapsedS,
        engineRpm: Math.min(1025, elapsedS * 180),
        dcRad: pose.canardRad,
        warpRad: pose.warpTipRad,
      });
      phaseEl.set(null, "info");
      // Attract-mode diorama: the machine idles at the rail head, so
      // Orville stands at the wingtip; gulls and the fire still live.
      dressing.animate(elapsedS, {
        onRail: true,
        aircraftX: 0,
        releaseX: null,
        releaseT: null,
      });
      if (pose.clamped) {
        console.warn(JSON.stringify({ suite: "wf-scene", event: "control-stop", t: elapsedS }));
      }
      renderer.render(scene, camera);
    },
    applyQuality(profile): void {
      // Atomic presentation-only application (E5.6): pixel ratio cap
      // and ghost visibility; terrain re-tessellation is a rebuild
      // concern (deferred, logged by the governor's JSONL).
      renderer.setPixelRatio(Math.min(window.devicePixelRatio, profile.pixelRatioCap));
      if (ghostFrame !== null) {
        ghostFrame.group.visible = profile.ghostVisible;
      }
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
      dials.dispose();
      phaseEl.dispose();
      renderer.dispose();
      container.removeChild(renderer.domElement);
    },
  };
}
