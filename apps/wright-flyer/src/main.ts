// First Flight — app entry. Bead wf-root-guzez.1.2 (E0.2).
// Scope: capability probe + honest banner, renderer seam, fixed-cadence FPS
// meter. The sim/field workers (E5.*) and real scenes (E2.*) land later.

import { describeCapabilities, probeCapabilities } from "./capability";
import { createFlyerSceneRenderer } from "./flyerScene.ts";
import { SimClient } from "./sim/simClient.ts";
import { MODE_FIXED, MODE_HISTORICAL, dec17Scenario } from "./sim/protocol.ts";

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
  let simClient: SimClient | undefined;
  if (params.get("sim") === "1") {
    simClient = new SimClient({
      onReady(info): void {
        simClient?.bindAnchor(info.tick0Digest);
        console.info(
          JSON.stringify({ suite: "wright-flyer-app", stage: "sim-ready", ...info }),
        );
      },
      onRefusal(stage, refusal): void {
        capabilityText.textContent = `sim ${stage} refusal: ${refusal.code} — ${refusal.message}`;
        capabilityText.className = "warn";
        console.warn(
          JSON.stringify({ suite: "wright-flyer-app", stage: "sim-refusal", at: stage, ...refusal }),
        );
      },
      onTerminal(info): void {
        console.info(
          JSON.stringify({ suite: "wright-flyer-app", stage: "sim-terminal", ...info }),
        );
      },
    });
    const mode = params.get("mode") === "historical" ? MODE_HISTORICAL : MODE_FIXED;
    simClient.start(dec17Scenario(1903n, mode));
  }
  const renderer = createFlyerSceneRenderer(app, simClient);
  const resize = (): void =>
    renderer.resize(app.clientWidth, app.clientHeight);
  window.addEventListener("resize", resize);
  resize();

  // FPS meter over a 1-second window (presentation-plane measurement only;
  // sim-tick metrics are E0.8's separate contract).
  let frames = 0;
  let windowStartMs = performance.now();
  let lastMs = windowStartMs;
  const frame = (nowMs: number): void => {
    renderer.render((nowMs - lastMs) / 1000);
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
