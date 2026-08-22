// First Flight — app entry. Bead wf-root-guzez.1.2 (E0.2).
// Scope: capability probe + honest banner, renderer seam, fixed-cadence FPS
// meter. The sim/field workers (E5.*) and real scenes (E2.*) land later.

import { describeCapabilities, probeCapabilities } from "./capability";
import { createFlyerSceneRenderer } from "./flyerScene.ts";
import { SimClient } from "./sim/simClient.ts";
import { MODE_FIXED, MODE_HISTORICAL, dec17Scenario, huffmanScenario } from "./sim/protocol.ts";
import { recordedToScenario, replayVerdict, type FlightRecording } from "./sim/replay.ts";
import { cardLines, computeKpis, kpiRecomputeDivergence } from "./sim/resultsCard.ts";
import { QosGovernor } from "./qos.ts";
import { MODE_HUMAN } from "./sim/protocol.ts";
import { LatencyLedger, toPhysical } from "./sim/humanControls.ts";
import { NEUTRAL, keysFrom, stepCommand, type PilotCommand } from "./input.ts";
import {
  DEFAULT_SELECTION,
  KEY_LINES,
  MODE_CARDS,
  assistAvailable,
  menuQuery,
  type MenuSelection,
} from "./menu.ts";

/** Landing menu overlay (game front door). The scripted demo keeps
 * running behind it as the attract mode; every button just navigates
 * to the URL params the app already honors. */
function buildMenu(container: HTMLElement): void {
  let sel: MenuSelection = DEFAULT_SELECTION;
  const overlay = document.createElement("div");
  overlay.id = "wf-menu";
  const card = document.createElement("div");
  card.id = "wf-menu-card";
  card.innerHTML =
    `<div id="wf-menu-title">FIRST FLIGHT</div>` +
    `<div id="wf-menu-sub">Wright Flyer · December 17, 1903 · a FrankenSim reconstruction</div>`;
  const modeRow = document.createElement("div");
  modeRow.className = "wf-menu-row";
  const modeButtons = new Map<string, HTMLButtonElement>();
  const refresh = (): void => {
    for (const [mode, btn] of modeButtons) {
      btn.classList.toggle("selected", mode === sel.mode);
    }
    const canAssist = assistAvailable(sel);
    assistBtn.classList.toggle("selected", sel.assist && canAssist);
    assistBtn.disabled = !canAssist;
    siteBtn.textContent =
      sel.site === "kdh" ? "SITE: KILL DEVIL HILLS 1903" : "SITE: HUFFMAN PRAIRIE 1904-05 (catapult)";
    go.textContent = `LAUNCH ${sel.site === "huffman" ? "BY CATAPULT" : "INTO THE HEADWIND"} →`;
  };
  for (const c of MODE_CARDS) {
    const btn = document.createElement("button");
    btn.className = "wf-mode-btn";
    btn.innerHTML = `<b>${c.title}</b><span>${c.blurb}</span>`;
    btn.addEventListener("click", () => {
      sel = { ...sel, mode: c.mode };
      refresh();
    });
    modeButtons.set(c.mode, btn);
    modeRow.appendChild(btn);
  }
  card.appendChild(modeRow);
  const optRow = document.createElement("div");
  optRow.className = "wf-menu-row";
  const siteBtn = document.createElement("button");
  siteBtn.className = "wf-opt-btn";
  siteBtn.addEventListener("click", () => {
    sel = { ...sel, site: sel.site === "kdh" ? "huffman" : "kdh" };
    refresh();
  });
  const assistBtn = document.createElement("button");
  assistBtn.className = "wf-opt-btn";
  assistBtn.textContent = "TRAINING ASSIST (bounded aid — not history)";
  assistBtn.addEventListener("click", () => {
    sel = { ...sel, assist: !sel.assist };
    refresh();
  });
  optRow.appendChild(siteBtn);
  optRow.appendChild(assistBtn);
  card.appendChild(optRow);
  const go = document.createElement("button");
  go.id = "wf-go-btn";
  go.addEventListener("click", () => {
    window.location.search = menuQuery(sel);
  });
  card.appendChild(go);
  const keys = document.createElement("pre");
  keys.id = "wf-menu-keys";
  keys.textContent = KEY_LINES.join("\n");
  card.appendChild(keys);
  const demo = document.createElement("a");
  demo.id = "wf-demo-link";
  demo.href = "?demo=1";
  demo.textContent = "just watch the scripted demo";
  card.appendChild(demo);
  overlay.appendChild(card);
  container.appendChild(overlay);
  refresh();
}

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
  // E5.5 results-card overlay (hidden until a run ends).
  const resultsCardEl = document.createElement("pre");
  resultsCardEl.id = "wf-results-card";
  resultsCardEl.style.display = "none";
  document.body.appendChild(resultsCardEl);
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
        // E5.5: the results card — KPIs recomputed from the sim-plane
        // transcript, gated by the recompute check before display.
        const recording = client.takeRecording();
        if (recording !== null) {
          const kpis = computeKpis(recording);
          const divergence = kpiRecomputeDivergence(recording, kpis);
          const site = params.get("site") === "huffman" ? "Huffman Prairie" : "Kill Devil Hills";
          const lines =
            divergence === null
              ? cardLines(kpis, site)
              : [`RESULTS CARD WITHHELD — KPI recompute divergence: ${divergence}`];
          resultsCardEl.textContent = lines.join("\n");
          resultsCardEl.style.display = "block";
          console.info(
            JSON.stringify({ suite: "wright-flyer-app", stage: "results-card", kpis, divergence }),
          );
        }
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
    // Historical mode flies the NONLINEAR-calibrated member 3 (the
    // E5.3b-i registration; members 0-2 are linear-plant tuned and
    // PIO out of the full lifecycle — kept for the H-campaigns).
    // E5.4: ?site=huffman selects the 1904-05 catapult scenario (the
    // near-calm run currently ends with a receipted envelope exit —
    // guzez.5.7.1 — surfaced on the HUD, never silently).
    const member = mode === MODE_HISTORICAL ? 3 : 0;
    simClient.start(
      params.get("site") === "huffman"
        ? huffmanScenario(1903n, mode, member)
        : dec17Scenario(1903n, mode, member, params.get("assist") === "1"),
    );
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
      resultsCardEl.style.display = "none";
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
  // Landing menu: front door when neither a sim run nor the explicit
  // scripted-demo view was requested (the demo keeps playing behind it).
  if (params.get("sim") !== "1" && params.get("demo") !== "1") {
    buildMenu(document.body);
  }
  window.addEventListener("resize", resize);
  resize();

  // FPS meter over a 1-second window (presentation-plane measurement only;
  // sim-tick metrics are E0.8's separate contract).
  let frames = 0;
  let windowStartMs = performance.now();
  let lastMs = windowStartMs;
  // E5.6: the hysteretic QoS governor — presentation only, badge
  // honest, one typed refusal on a persistent budget miss.
  const qos = new QosGovernor();
  const badgeEl = document.createElement("div");
  badgeEl.style.cssText =
    "position:fixed;left:12px;top:12px;display:none;font:12px monospace;" +
    "color:#1a1a1a;background:#e8c76a;padding:4px 8px;border-radius:4px";
  document.body.appendChild(badgeEl);
  let lastPublishedTick = 0;
  const frame = (nowMs: number): void => {
    const dtS = (nowMs - lastMs) / 1000;
    const q = qos.sample(Math.max(0, nowMs - lastMs));
    if (q.changed) {
      renderer.applyQuality?.(q.profile);
      badgeEl.style.display = q.profile.badge !== null ? "block" : "none";
      badgeEl.textContent = q.profile.badge ?? "";
      console.info(
        JSON.stringify({ suite: "wright-flyer-app", stage: "qos", state: q.state, ...q.profile }),
      );
    }
    if (q.refusal !== undefined) {
      console.warn(JSON.stringify({ suite: "wright-flyer-app", stage: "qos-refusal", ...q.refusal }));
      capabilityText.textContent = `${q.refusal.code}: ${q.refusal.message}`;
      capabilityText.className = "warn";
    }
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
