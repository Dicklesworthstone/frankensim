// E2.2 scene: the parametric airframe in the dawn scene, driven by the
// scripted demo (bead wf-root-guzez.3.2). Implements the FlyerRenderer
// seam so main.ts swaps one factory call.

import * as THREE from "three";
import { buildWrightFlyerAirframe } from "./airframe/parametricAirframe.ts";
import { driveScripted } from "./airframe/applyPose.ts";
import { FlightAudio } from "./audio.ts";
import {
  bigHillDetail,
  buildTerrainArrays,
  duneDetail,
  heightAt,
} from "./terrainMesh.ts";
import {
  type CameraPreset,
  type CameraState,
  PRESET_KEYS,
  cameraFor,
  easeCameraToward,
  speedFov,
} from "./camera.ts";
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
import { orvilleReachableX } from "./dressing.ts";
import {
  SUN_COLOR,
  SUN_DIRECTION,
  buildDressing,
  buildTakeoffDolly,
  sandTileMaterial,
} from "./dressing3d.ts";
import { createProneBrother } from "./figure3d.ts";

// ONE soundscape for the whole app — module-level singleton so scene
// rebuilds (R replay) never stack AudioContexts or event listeners.
// The context starts on the first user gesture; M toggles mute (T4.3);
// Daniels' flash publishes `wf-flash` for the shutter click.
const AUDIO = new FlightAudio({ withOcean: true });
const onFirstGesture = (): void => {
  AUDIO.ensureStarted();
  window.removeEventListener("pointerdown", onFirstGesture);
  window.removeEventListener("keydown", onFirstGesture);
};
window.addEventListener("pointerdown", onFirstGesture);
window.addEventListener("keydown", onFirstGesture);
window.addEventListener("wf-flash", () => AUDIO.shutter());


export function createFlyerSceneRenderer(
  container: HTMLElement,
  simClient?: SimClient,
  ghost?: FlightRecording,
): FlyerRenderer {
  // Logarithmic depth: the diorama spans 0.1 m (cockpit) to 2.6 km
  // (sky dome); a linear 24-bit buffer z-fights the ocean/skirt/tile
  // layers at that range.
  const renderer = new THREE.WebGLRenderer({ antialias: true, logarithmicDepthBuffer: true });
  renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
  // Real shadows + filmic dawn grade: the single biggest presentation
  // upgrades available. castShadow flags across the diorama were inert
  // until the shadow map existed; ACES keeps the low-sun highlights
  // from clipping to white.
  renderer.shadowMap.enabled = true;
  renderer.shadowMap.type = THREE.PCFSoftShadowMap;
  renderer.toneMapping = THREE.ACESFilmicToneMapping;
  renderer.toneMappingExposure = 1.12;
  container.appendChild(renderer.domElement);
  const scene = new THREE.Scene();
  // December-morning haze matched to the rebuilt sky's horizon band;
  // starts far enough out that the DUNES stay visible from altitude.
  scene.background = new THREE.Color(0xd7dde2);
  scene.fog = new THREE.Fog(0xd9dee2, 260, 2400);
  // Far plane must clear the sky dome (2600 m) and the Atlantic.
  const camera = new THREE.PerspectiveCamera(50, 1, 0.1, 5000);
  camera.position.set(9, 3.2, 11);
  camera.lookAt(0, 1.2, 0);
  // Presentation-plane smoothing state: every preset target glides
  // toward the lens instead of teleporting (T0.5).
  let camState: CameraState = { pos: [9, 3.2, 11], look: [0, 1.2, 0], fovDeg: 50 };
  // Particle budget from the QoS tier (0 = full, 2 = hidden).
  let particleLevel: 0 | 1 | 2 = 0;
  // ONE sun: the light sits along the SAME direction the sky texture
  // paints its disc from (dressing3d.SUN_DIRECTION), warm like a
  // 10:35 a.m. December sun. Shadow ortho box covers the launch flat,
  // rail corridor, and camp — not the whole 2 km tile.
  const sun = new THREE.DirectionalLight(SUN_COLOR, 2.6);
  sun.position.set(
    SUN_DIRECTION[0] * 160,
    SUN_DIRECTION[1] * 160,
    SUN_DIRECTION[2] * 160,
  );
  sun.castShadow = true;
  sun.shadow.mapSize.set(2048, 2048);
  sun.shadow.camera.left = -80;
  sun.shadow.camera.right = 110;
  sun.shadow.camera.top = 70;
  sun.shadow.camera.bottom = -60;
  sun.shadow.camera.near = 40;
  sun.shadow.camera.far = 340;
  sun.shadow.bias = -0.00035;
  sun.shadow.normalBias = 0.02;
  scene.add(sun, new THREE.HemisphereLight(0xbcd0e8, 0x8a7d64, 0.85));
  // The REAL site tile (E1.3 3DEP grids): heightfield + splat.
  // E5.4: ?site=huffman loads the Huffman Prairie grid.
  const site = new URLSearchParams(window.location.search).get("site");
  // Presentation-only relief layer (terrainMesh.ts): rolling dune sea
  // Big Kill Devil Hill at Kill Devil Hills; gentle swells only for
  // landlocked Huffman pasture. The launch/camp corridor stays at
  // survey height by the mask inside duneDetail.
  const siteRelief = (xRel: number, zRel: number): number =>
    site === "huffman"
      ? duneDetail(xRel, zRel, { amplitudeM: 2.2 })
      : duneDetail(xRel, zRel) + bigHillDetail(xRel, zRel);
  const terrain = buildTerrainArrays(site === "huffman" ? huffmanGrid : kdhGrid, 256, siteRelief);
  const tGeo = new THREE.BufferGeometry();
  tGeo.setAttribute("position", new THREE.BufferAttribute(terrain.positions, 3));
  tGeo.setAttribute("color", new THREE.BufferAttribute(terrain.colors, 3));
  tGeo.setAttribute("uv", new THREE.BufferAttribute(terrain.uvs, 2));
  tGeo.setIndex(new THREE.BufferAttribute(terrain.indices, 1));
  tGeo.computeVertexNormals();
  scene.add(new THREE.Mesh(tGeo, sandTileMaterial()));
  const airframe = buildWrightFlyerAirframe();
  const launch = terrain.launch;
  airframe.group.position.set(launch[0], launch[1] + 1.2, launch[2]);
  scene.add(airframe.group);
  // The player's pilot (guzez.14): the anthropometric Wilbur — PD
  // portrait face, hands coupled to the REAL canard deflection below.
  const wilburFig = createProneBrother("wilbur");
  wilburFig.group.position.set(-0.1, 0.12, 0.06);
  airframe.cradleGroup.add(wilburFig.group);
  // The Kitty Hawk diorama (bead guzez.13): sky, clouds, outer sand +
  // Atlantic, the launch rail, the 1903 camp, Orville, and the gulls.
  const grid = site === "huffman" ? huffmanGrid : kdhGrid;
  const tileExtent = (grid.grid_n - 1) * grid.spacing_m;
  const half = tileExtent / 2;
  // Grids store ABSOLUTE elevations (KDH ~0 m, Huffman ~242 m); the
  // dome, clouds, and sand skirt anchor to the tile's lowest point.
  const baseY = Math.min(...grid.rows_south_to_north.flat());
  // Scenario headwind drives every wind-made-visible system (T1.6-8).
  const windMps = site === "huffman" ? 2.0 : 11.0;
  const dressing = buildDressing(
    launch,
    site === "huffman" ? 30.0 : 18.3,
    tileExtent,
    site !== "huffman", // the Atlantic is a Kill Devil Hills fact
    baseY,
    (xRel, zRel) =>
      heightAt(grid, launch[0] + xRel + half, -(launch[2] + zRel) + half) -
      launch[1] +
      siteRelief(xRel, zRel),
    windMps,
  );
  // Trim prop speed [rad/s] — engine 1025 rpm through the 23:8 chain —
  // the normalization for the plume strength (presentation-only).
  const TRIM_PROP_OMEGA = ((1025 * (8 / 23)) / 60) * 2 * Math.PI;
  scene.add(dressing.group);
  // Orville's release point is LATCHED the first frame the machine is
  // off the rail (pure pose math needs the release constants).
  let orvilleReleaseX: number | null = null;
  let orvilleReleaseT: number | null = null;
  let elapsedS = 0;
  // E2.4: input, cameras, HUD. The pilot takes over from the script the
  // first time a control key goes down. Human mode STARTS in the prone
  // pilot's seat (T3.4) — you are Wilbur.
  const down = new Set<string>();
  const humanMode = new URLSearchParams(window.location.search).get("mode") === "human";
  // Third person BEHIND Wilbur is the default piloting view (booting
  // straight into the prone first-person eye was disorienting: no
  // aircraft in frame, sand rushing). V toggles chase <-> onboard.
  let preset: CameraPreset = humanMode ? "chase" : "free";
  let manual = false;
  let command = NEUTRAL;
  // The takeoff dolly (T2.1): under the skids on the rail; it STAYS on
  // the track at the liftoff spot once the machine is off the rail.
  const dolly = buildTakeoffDolly();
  dolly.position.set(launch[0], launch[1], launch[2]);
  scene.add(dolly);
  let dollyDropped = false;
  // Trim prop speed [rad/s] — engine 1025 rpm through the 23:8 chain —
  // the normalization for the plume strength (presentation-only).
  // Landing latch: once the machine is DOWN, Orville runs to it and
  // the dust burst fires (T2.6/T3.5).
  let landedX: number | null = null;
  let landedT: number | null = null;
  const hud = document.createElement("div");
  hud.id = "wf-telemetry";
  hud.style.cssText =
    "position:fixed;right:12px;bottom:12px;font:11px/1.5 monospace;color:#f5efe0;" +
    "background:rgba(20,24,30,.55);padding:8px 10px;border-radius:6px;white-space:pre;z-index:5";
  hud.style.display = "none"; // raw telemetry opt-in via T (clutter fix)
  container.appendChild(hud);
  // Game HUD: analog dial panel + phase banner (tested math in gauges.ts).
  const dials = createHudDials(container);
  const phaseEl = createPhaseBanner(container);
  // Always-visible controls card (H hides it): the fix for "the
  // controls are confusing" is showing them, in flight, at all times.
  const helpCard = document.createElement("div");
  helpCard.id = "wf-help-card";
  helpCard.style.cssText =
    "position:fixed;left:12px;bottom:12px;font:12px/1.7 monospace;color:#f0e4c8;" +
    "background:rgba(32,22,12,.78);border:1px solid #8a6a38;padding:10px 12px;" +
    "border-radius:6px;white-space:pre;z-index:6";
  helpCard.textContent = simClient
    ? "CONTROLS\n" +
      "S or ↓   pull — nose UP\n" +
      "W or ↑   push — nose DOWN\n" +
      "A/D ←/→  wing warp (bank)\n" +
      "Space    recenter controls\n" +
      "V        camera: behind ↔ pilot's eyes\n" +
      "M sound · H hide this · T telemetry"
    : "1-5 cameras · M sound · H hide";
  container.appendChild(helpCard);
  const onKey = (e: KeyboardEvent, isDown: boolean): void => {
    if (isDown && e.code === "KeyT" && !e.repeat) {
      hud.style.display = hud.style.display === "none" ? "block" : "none";
      return;
    }
    if (isDown && e.code === "KeyM" && !e.repeat) {
      console.info(
        JSON.stringify({ suite: "wf-scene", event: "audio-mute", muted: AUDIO.toggleMute() }),
      );
      return;
    }
    // V: the one camera key a player needs — third person behind
    // Wilbur <-> his own eyes (number keys still reach every preset).
    if (isDown && e.code === "KeyV" && !e.repeat) {
      preset = preset === "onboard" ? "chase" : "onboard";
      return;
    }
    if (isDown && e.code === "KeyH" && !e.repeat) {
      helpCard.style.display = helpCard.style.display === "none" ? "block" : "none";
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
            // Nose-frame mapping (visual-verification fix): the airframe's
            // nose is LOCAL +z; flight runs WORLD +x. Ry(-pi/2) leftmost
            // ('YXZ' order) maps nose->+x and wingspan->z AFTER the pitch.
            ghostFrame.group.rotation.order = "YXZ";
            ghostFrame.group.rotation.set(0, -Math.PI / 2, gw.pitchRad);
          }
        }
        drive = advanceProp(drive, snap, dtS);
        const pose = computePose(controlStateFrom(snap, drive));
        applyPose(airframe, pose);
        const world = worldTransformFrom(snap, launch);
        airframe.group.position.set(world.position[0], world.position[1] + 1.2, world.position[2]);
        // Nose-frame mapping (visual-verification fix): the airframe's
        // nose is LOCAL +z; flight runs WORLD +x. Ry(-pi/2) leftmost
        // ('YXZ' order) maps nose->+x and wingspan->z AFTER the pitch.
        airframe.group.rotation.order = "YXZ";
        airframe.group.rotation.set(0, -Math.PI / 2, world.pitchRad);
        // Cockpit feel (T3.4): the pilot's lever mirrors live canard.
        airframe.elevatorLever.rotation.z = -0.22 - pose.canardRad * 0.9;
        wilburFig.setLever(pose.canardRad);
        const gust01 = Math.min(1, Math.abs(snap.gustWMps) / 1.5);
        // Muslin flutter (T1.9): the wing fabric's weave map trembles
        // with the live gust sample — texture-space only, presentation.
        const muslinTex = airframe.textures[1]!;
        muslinTex.offset.set(
          Math.sin(elapsedS * 11.3) * 0.0045 * gust01,
          Math.cos(elapsedS * 7.1) * 0.003 * gust01,
        );
        const cam = cameraFor(preset, elapsedS, launch, [
          world.position[0],
          world.position[1] + 1.2,
          world.position[2],
        ]);
        camState = easeCameraToward(camState, cam, dtS);
        camera.position.set(camState.pos[0], camState.pos[1], camState.pos[2]);
        camera.lookAt(camState.look[0], camState.look[1], camState.look[2]);
        // Onboard head-bob: the prone pilot's head nods with gusts and
        // pitch rate — cockpit life without touching the sim state.
        if (preset === "onboard") {
          camera.position.y += Math.sin(elapsedS * 8.7) * 0.018 * gust01;
          camera.position.z += Math.sin(elapsedS * 5.3) * 0.012 * gust01;
        }
        const fov = speedFov(50, Math.hypot(snap.uMps, snap.wMps));
        if (Math.abs(fov - camera.fov) > 0.01) {
          camera.fov = fov;
          camera.updateProjectionMatrix();
        }
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
        // Wilbur's hands ride the REAL canard deflection (guzez.14).
        wilburFig.setLever(snap.dcRad);
        // Diorama: Orville chases while the machine is ON the rail,
        // lets go at the first off-rail frame (latched), then watches.
        const onRail = snap.phase === "on-rail";
        if (!onRail && orvilleReleaseX === null) {
          orvilleReleaseX = orvilleReachableX(hudIn.elapsedS, snap.xM);
          orvilleReleaseT = hudIn.elapsedS;
        }
        if (snap.phase === "ended:ground-contact") {
          if (landedX === null) {
            landedX = snap.xM;
            landedT = hudIn.elapsedS;
          }
        }
        dressing.animate(hudIn.elapsedS, {
          onRail,
          aircraftX: snap.xM,
          releaseX: orvilleReleaseX,
          releaseT: orvilleReleaseT,
          landedX,
          // Machine state feeds the plume systems (propwash sand blast
          // + exhaust + landing dust) — presentation-only normalization.
          machine: {
            x: snap.xM,
            rpm01: Math.min(1, Math.max(0, snap.omegaPropRadS / TRIM_PROP_OMEGA)),
            gust01,
            dustT: landedT ?? undefined,
          },
        });
        // The dolly rides the rail under the machine until liftoff,
        // then STAYS VISIBLE at the drop spot for the rest of the run
        // (T2.1 — a real object left on the track, never hidden).
        if (!dollyDropped) {
          if (!onRail) {
            dollyDropped = true;
          } else {
            dolly.position.x = launch[0] + snap.xM;
          }
        }
        // Labeled ride-along captions (latest two).
        if (snap.tick > lastCaptionTick) {
          lastCaptionTick = snap.tick;
          captions.feed(snap);
        }
        for (const c of captions.upTo(snap.tick).slice(-2)) {
          lines.push(formatCaption(c));
        }
        hud.textContent = lines.join("\n");
        // Live mix (T4.3): engine/wind/rail from the real state.
        AUDIO.update({
          propOmegaRadS: snap.omegaPropRadS,
          airspeedMps: Math.hypot(snap.uMps, snap.wMps),
          onRail,
          groundSpeedMps: onRail ? Math.max(0, snap.xM / Math.max(hudIn.elapsedS, 0.5)) : 0,
          nowS: hudIn.elapsedS,
        });
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
      camState = easeCameraToward(camState, cam, dtS);
      camera.position.set(camState.pos[0], camState.pos[1], camState.pos[2]);
      camera.lookAt(camState.look[0], camState.look[1], camState.look[2]);
      const fovIdle = speedFov(50, 10.73);
      if (Math.abs(fovIdle - camera.fov) > 0.01) {
        camera.fov = fovIdle;
        camera.updateProjectionMatrix();
      }
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
        // Attract idle: props are static (no plumes) — no machine state.
      });
      dolly.visible = true;
      dolly.position.x = launch[0];
      // Attract mix: the scripted demo's engine ramp, no rail motion.
      AUDIO.update({
        propOmegaRadS: (Math.min(350, elapsedS * 180) / 60) * 2 * Math.PI,
        airspeedMps: 0,
        onRail: true,
        groundSpeedMps: 0,
        nowS: elapsedS,
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
      // T0.6: shadow + particle budgets ride the existing presentation
      // tiers — Critical (fieldThrottle 2) drops the shadow map and
      // hides the particle systems; Constrained halves their budget.
      const shadows = profile.fieldThrottle < 2;
      if (renderer.shadowMap.enabled !== shadows) {
        renderer.shadowMap.enabled = shadows;
        scene.traverse((obj) => {
          const mesh = obj as THREE.Mesh;
          if (mesh.isMesh && Array.isArray(mesh.material)) {
            mesh.material.forEach((m) => (m.needsUpdate = true));
          } else if (mesh.isMesh) {
            (mesh.material as THREE.Material).needsUpdate = true;
          }
        });
      }
      particleLevel = profile.fieldThrottle;
      dressing.setParticleLevel(particleLevel);
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
      helpCard.remove();
      dials.dispose();
      phaseEl.dispose();
      renderer.dispose();
      container.removeChild(renderer.domElement);
    },
  };
}
