// First Flight — app entry. Bead wf-root-guzez.1.2 (E0.2).
// Scope: capability probe + honest banner, renderer seam, fixed-cadence FPS
// meter. The sim/field workers (E5.*) and real scenes (E2.*) land later.

import { describeCapabilities, probeCapabilities } from "./capability";
import { createFlyerSceneRenderer } from "./flyerScene.ts";
import { SimClient } from "./sim/simClient.ts";
import { MODE_FIXED, MODE_HISTORICAL, dec17Scenario } from "./sim/protocol.ts";
import { recordedToScenario, replayVerdict, type FlightRecording } from "./sim/replay.ts";
import { MODE_HUMAN } from "./sim/protocol.ts";
import { LatencyLedger, toPhysical } from "./sim/humanControls.ts";
import { NEUTRAL, keysFrom, stepCommand, type PilotCommand } from "./input.ts";

function main(): void {
  const app = document.getElementById("app");
  const capabilityText = document.getElementById("capability-text");
  const fpsText = document.getElementById("fps-text");
  if (!app || !capabilityText || !fpsText) {
    throw new Error("index.html scaffold elements missing");
  }

  const caps = probeCapabilities();
  const described = describeCapabilities(caps);
  capabilityText.textContent = described.text;
  capabilityText.className = described.degraded ? "warn" : "";
  // Structured, greppable capability log line (JSONL-style, one object).
  console.info(
    JSON.stringify({
      suite: "wright-flyer-app",
      stage: "capability-probe",
      ...caps,
      degraded: described.degraded,
    }),
  );

  // E5.2b: ?sim=1 drives the scene from the REAL wasm engine (E5.1)
  // via the sim worker; without the flag the E2.2 scripted demo runs.
  // ?mode=historical selects the registered pilot family member 0.
  const params = new URLSearchParams(window.location.search);
  const modeWord = params.get("mode");
  const mode =
    modeWord === "historical" ? MODE_HISTORICAL : modeWord === "human" ? MODE_HUMAN : MODE_FIXED;
  // E5.3a: latency decomposition ledger (device→sent→ack→published→
  // present), JSONL per completed control sample.
  const ledger = new LatencyLedger((line) => console.info(line));
  let simClient: SimClient | undefined;
  let renderer = createFlyerSceneRenderer(app);
  const resize = (): void => renderer.resize(app.clientWidth, app.clientHeight);

  const makeClient = (expected?: FlightRecording): SimClient => {
    const client: SimClient = new SimClient({
      onReady(info): void {
        client.bindAnchor(info.tick0Digest);
        // Clock-sync burst (min-RTT sample wins inside the client).
        for (let i = 0; i < 5; i += 1) {
          client.sendPing();
        }
        // Human mode: the run starts at the first admitted control —
        // send neutral so the rail run begins immediately.
        if (mode === MODE_HUMAN) {
          sendHumanControl(NEUTRAL, performance.now());
        }
        console.info(JSON.stringify({ suite: "wright-flyer-app", stage: "sim-ready", ...info }));
      },
      onRefusal(stage, refusal): void {
        capabilityText.textContent = `sim ${stage} refusal: ${refusal.code} — ${refusal.message}`;
        capabilityText.className = "warn";
        console.warn(
          JSON.stringify({ suite: "wright-flyer-app", stage: "sim-refusal", at: stage, ...refusal }),
        );
      },
      onTerminal(info): void {
        console.info(JSON.stringify({ suite: "wright-flyer-app", stage: "sim-terminal", ...info }));
        // E5.2c: replay identity verdict against the previous run —
        // the engine's chained digest, surfaced on screen and in the log.
        if (expected !== undefined) {
          const verdict = replayVerdict(expected, info.digest);
          capabilityText.textContent =
            verdict.kind === "identical"
              ? `REPLAY IDENTICAL — digest ${info.digest.slice(0, 16)}…`
              : `REPLAY DIVERGED — expected ${verdict.expectedDigest.slice(0, 16)}… observed ${verdict.observedDigest.slice(0, 16)}…`;
          capabilityText.className = verdict.kind === "identical" ? "" : "warn";
          console.info(
            JSON.stringify({ suite: "wright-flyer-app", stage: "replay-verdict", ...verdict }),
          );
        }
        capabilityText.textContent += "  [R replays with ghost]";
      },
      onControlAck(ack): void {
        ledger.acked(ack.sequence, ack.appliedTick, ack.lateByTicks, performance.now());
      },
    });
    return client;
  };

  // E5.3a human-control pump: keyboard-rate transducer → quantized
  // command → physical units → worker (with the latency ledger fed).
  let controlSeq = 0;
  let command: PilotCommand = NEUTRAL;
  let lastSent: PilotCommand | null = null;
  let lastSentMs = 0;
  const heldKeys = new Set<string>();
  const sendHumanControl = (cmd: PilotCommand, deviceMs: number): void => {
    if (simClient === undefined) {
      return;
    }
    controlSeq += 1;
    const phys = toPhysical(cmd);
    ledger.sent(controlSeq, deviceMs, performance.now());
    simClient.sendControl(phys.leverForceN, phys.warpCmdRad, controlSeq, deviceMs);
    lastSent = cmd;
    lastSentMs = performance.now();
  };
  const pumpHuman = (nowMs: number, dtS: number): void => {
    if (mode !== MODE_HUMAN || simClient === undefined) {
      return;
    }
    command = stepCommand(command, keysFrom(heldKeys), dtS);
    const changed =
      lastSent === null || command.canard !== lastSent.canard || command.warp !== lastSent.warp;
    // Send on change, plus a 100 ms heartbeat (the worker holds ZOH).
    if (changed || nowMs - lastSentMs > 100) {
      sendHumanControl(command, nowMs);
    }
  };
  if (mode === MODE_HUMAN) {
    window.addEventListener("keydown", (e: KeyboardEvent) => heldKeys.add(e.code));
    window.addEventListener("keyup", (e: KeyboardEvent) => heldKeys.delete(e.code));
  }

  if (params.get("sim") === "1") {
    simClient = makeClient();
    simClient.start(dec17Scenario(1903n, mode));
    renderer.dispose();
    renderer = createFlyerSceneRenderer(app, simClient);
    // R after a terminal: rebuild the scene with the finished run as a
    // ghost and rerun the SAME scenario for the on-screen overlay test.
    window.addEventListener("keydown", (e: KeyboardEvent) => {
      if (e.code !== "KeyR" || simClient === undefined) {
        return;
      }
      const recording = simClient.takeRecording();
      if (recording === null) {
        return; // run not finished yet
      }
      simClient.dispose();
      simClient = makeClient(recording);
      simClient.start(recordedToScenario(recording.scenario));
      renderer.dispose();
      renderer = createFlyerSceneRenderer(app, simClient, recording);
      resize();
      console.info(
        JSON.stringify({
          suite: "wright-flyer-app",
          stage: "replay-start",
          frames: recording.ticks.length,
          expectedDigest: recording.finalDigest,
        }),
      );
    });
  }
  window.addEventListener("resize", resize);
  resize();

  // FPS meter over a 1-second window (presentation-plane measurement only;
  // sim-tick metrics are E0.8's separate contract).
  let frames = 0;
  let windowStartMs = performance.now();
  let lastMs = windowStartMs;
  let lastPublishedTick = 0;
  const frame = (nowMs: number): void => {
    const dtS = (nowMs - lastMs) / 1000;
    pumpHuman(nowMs, dtS);
    renderer.render(dtS);
    // Latency ledger: publication (new sim tick visible) then present.
    const tick = simClient?.latestTick() ?? 0;
    if (tick > lastPublishedTick) {
      lastPublishedTick = tick;
      ledger.published(tick, nowMs);
    }
    ledger.presented(performance.now());
    lastMs = nowMs;
    frames += 1;
    if (nowMs - windowStartMs >= 1000) {
      const fps = (frames * 1000) / (nowMs - windowStartMs);
      fpsText.textContent = `${fps.toFixed(0)} fps`;
      console.info(
        JSON.stringify({
          suite: "wright-flyer-app",
          stage: "fps-window",
          fps: Number(fps.toFixed(1)),
        }),
      );
      frames = 0;
      windowStartMs = nowMs;
    }
    requestAnimationFrame(frame);
  };
  requestAnimationFrame(frame);
}

main();
