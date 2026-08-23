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
import {
  NEUTRAL,
  cradleFromPointer,
  decayCradle,
  keysFrom,
  sampleGamepad,
  stepCommand,
  type PilotCommand,
} from "./input.ts";
import {
  DEFAULT_SELECTION,
  FLIGHT_CHIPS,
  KEY_LINES,
  MODE_CARDS,
  assistAvailable,
  menuQuery,
  type MenuSelection,
} from "./menu.ts";
import { flightByIndex, flightSeed, missionOutcome } from "./missions/flights.ts";
import { JOURNEY_STAGES, journeyNextUrl, journeyStage } from "./journey.ts";
import { togglePhotoMode } from "./photoMode.ts";
import { scoreTouchdown } from "./landingScore.ts";
import { createInstrumentsPanel, type InstrumentSample } from "./instruments.ts";
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
  const flightChips = new Map<number, HTMLButtonElement>();
  const refresh = (): void => {
    for (const [mode, btn] of modeButtons) {
      btn.classList.toggle("selected", mode === sel.mode);
    }
    const canAssist = assistAvailable(sel);
    assistBtn.classList.toggle("selected", sel.assist && canAssist);
    assistBtn.disabled = !canAssist;
    for (const [id, chip] of flightChips) {
      const usable = sel.site === "kdh";
      chip.classList.toggle("selected", usable && sel.flight === id);
      chip.disabled = !usable;
    }
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
  // Mission chips: the four Dec-17 flights as scenario presets (KDH
  // only — Huffman disables them, mirroring the query-string law).
  const flightRow = document.createElement("div");
  flightRow.className = "wf-menu-row";
  const freeChip = document.createElement("button");
  freeChip.className = "wf-opt-btn";
  freeChip.textContent = "FREE ENSEMBLE";
  freeChip.addEventListener("click", () => {
    sel = { ...sel, flight: undefined };
    refresh();
  });
  flightRow.appendChild(freeChip);
  for (const c of FLIGHT_CHIPS) {
    const chip = document.createElement("button");
    chip.className = "wf-opt-btn";
    chip.textContent = c.label;
    chip.addEventListener("click", () => {
      sel = { ...sel, flight: c.id };
      refresh();
    });
    flightChips.set(c.id, chip);
    flightRow.appendChild(chip);
  }
  card.appendChild(flightRow);
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
  const journeyBtn = document.createElement("button");
  journeyBtn.className = "wf-opt-btn";
  journeyBtn.textContent = "GUIDED JOURNEY — watch a modeled pilot, then fly (assist first)";
  journeyBtn.addEventListener("click", () => {
    window.location.search = JOURNEY_STAGES[0]!.url;
  });
  card.appendChild(journeyBtn);
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
  // Guided journey overlay (?journey=N): caption while the stage runs,
  // continue-link when its run ends (J advances early). Presentation
  // only — the stage chain lives in journey.ts and advances by URL.
  const stage = journeyStage(params.get("journey"));
  let journeyEl: HTMLDivElement | null = null;
  if (stage !== null) {
    journeyEl = document.createElement("div");
    journeyEl.id = "wf-journey";
    const captionEl = document.createElement("div");
    captionEl.className = "wf-journey-caption";
    captionEl.textContent = stage.caption;
    const nextEl = document.createElement("div");
    nextEl.className = "wf-journey-next";
    journeyEl.appendChild(captionEl);
    journeyEl.appendChild(nextEl);
    document.body.appendChild(journeyEl);
    console.info(
      JSON.stringify({ suite: "wright-flyer-app", stage: "journey-stage", index: stage.index }),
    );
  }
  let simClient: SimClient | undefined;
  // Startup-watchdog state (sim block): the frame loop turns a 6 s
  // tickless window into a LOUD typed warning. Human mode is exempt —
  // waiting for the pilot's first control is by design.
  let simStartMs = performance.now();
  let watchdogFired = false;
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
          const base =
            divergence === null
              ? cardLines(kpis, site)
              : [`RESULTS CARD WITHHELD — KPI recompute divergence: ${divergence}`];
          const flightNum = Number(params.get("flight"));
          const mission =
            params.get("site") !== "huffman" && Number.isInteger(flightNum)
              ? flightByIndex(flightNum)
              : null;
          const touchdown = divergence === null ? scoreTouchdown(recording) : null;
          const lines =
            mission !== null && divergence === null
              ? [...base, ...missionOutcome(mission, kpis.downrangeM, kpis.airborneS).lines]
              : base;
          if (touchdown !== null) {
            lines.push(touchdown.line);
          }
          resultsCardEl.textContent = lines.join("\n");
          resultsCardEl.style.display = "block";
          // Telegram slide-in (T3.5): restart the animation each run.
          resultsCardEl.classList.remove("wf-card-in");
          void resultsCardEl.offsetWidth;
          resultsCardEl.classList.add("wf-card-in");
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
        if (stage !== null && journeyEl !== null) {
          const nextUrl = journeyNextUrl(stage.index);
          const slot = journeyEl.children[1] as HTMLElement;
          slot.textContent = "";
          const link = document.createElement("a");
          link.href = nextUrl ?? "?demo=1";
          link.textContent = nextUrl !== null ? `${stage.prompt}  CONTINUE →` : `${stage.prompt}  FINISH`;
          slot.appendChild(link);
        }
      },
      onControlAck(ack): void {
        ledger.acked(ack.sequence, ack.appliedTick, ack.lateByTicks, performance.now());
      },
    });
    return client;
  };

  // E5.3a human-control pump, extended for plan §2.3 devices: three
  // transducer families (keyboard-rate slew, pointer hip-cradle
  // position, gamepad-stick position) feed ONE quantized path — the
  // 1/4096 grid is the trace identity, so engine and replays see the
  // same command shape regardless of device. Priority: cradle (while
  // grabbed or springing back) > gamepad (while deflected, until
  // another device is touched) > keyboard.
  let controlSeq = 0;
  let command: PilotCommand = NEUTRAL;
  let lastSent: PilotCommand | null = null;
  let lastSentMs = 0;
  const heldKeys = new Set<string>();
  const CONTROL_KEY: Record<string, true> = {
    ArrowUp: true, ArrowDown: true, ArrowLeft: true, ArrowRight: true,
    KeyW: true, KeyA: true, KeyS: true, KeyD: true, Space: true,
  };
  let cradleActive = false;
  let cradleCmd: PilotCommand | null = null; // non-null = the cradle owns the pump
  let grabX = 0;
  let grabY = 0;
  let padOwned = false;
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
    let padCmd: PilotCommand | null = null;
    if (typeof navigator.getGamepads === "function") {
      for (const pad of navigator.getGamepads()) {
        const sampled = sampleGamepad(pad);
        if (sampled !== null) {
          padCmd = sampled;
          break;
        }
      }
    }
    if (padCmd !== null && (padCmd.canard !== 0 || padCmd.warp !== 0)) {
      padOwned = true;
    }
    let next: PilotCommand;
    if (cradleCmd !== null) {
      next = cradleActive ? cradleCmd : decayCradle(cradleCmd, dtS);
      cradleCmd = !cradleActive && next.canard === 0 && next.warp === 0 ? null : next;
    } else if (padOwned && padCmd !== null) {
      next = padCmd;
    } else {
      next = stepCommand(command, keysFrom(heldKeys), dtS);
      padOwned = false;
    }
    command = next;
    const changed =
      lastSent === null || command.canard !== lastSent.canard || command.warp !== lastSent.warp;
    // Send on change, plus a 100 ms heartbeat (the worker holds ZOH).
    if (changed || nowMs - lastSentMs > 100) {
      sendHumanControl(command, nowMs);
    }
  };
  if (mode === MODE_HUMAN) {
    window.addEventListener("keydown", (e: KeyboardEvent) => {
      heldKeys.add(e.code);
      if (CONTROL_KEY[e.code] === true) {
        padOwned = false;
        cradleActive = false;
        cradleCmd = null; // last device touched wins
      }
    });
    window.addEventListener("keyup", (e: KeyboardEvent) => heldKeys.delete(e.code));
    // Hip cradle: drag anywhere on the view. Listeners live on the
    // persistent #app container because renderer canvases are swapped
    // by replay/restart, and pointer events bubble from them.
    app.addEventListener("pointerdown", (e: PointerEvent) => {
      if (!e.isPrimary || e.button !== 0) {
        return;
      }
      grabX = e.clientX;
      grabY = e.clientY;
      cradleActive = true;
      cradleCmd = cradleFromPointer(0, 0);
      padOwned = false;
      e.preventDefault();
    });
    window.addEventListener("pointermove", (e: PointerEvent) => {
      if (!cradleActive || !e.isPrimary) {
        return;
      }
      cradleCmd = cradleFromPointer(e.clientX - grabX, e.clientY - grabY);
    });
    const endCradle = (): void => {
      cradleActive = false; // spring-back decay proceeds in the pump
    };
    window.addEventListener("pointerup", endCradle);
    window.addEventListener("pointercancel", endCradle);
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
    // Mission preset (?flight=N): each Dec-17 flight flies its own
    // deterministic wind ensemble; anything else keeps the 1903 base.
    const missionNum = Number(params.get("flight"));
    const missionSeed =
      params.get("site") !== "huffman" && Number.isInteger(missionNum) && missionNum >= 1 && missionNum <= 4
        ? flightSeed(missionNum)
        : 1903n;
    // HONESTY LAW (menu.ts/protocol.ts agreement): the bounded aid flies
    // ONLY when the URL claims it — the menu emits ?assist=1 on its
    // toggle, journey stage 2 passes it explicitly, and stage 3 plus
    // every unparametered run get the authentic machine. A silent
    // default-on here would launder a modern aid into the historical
    // and fixed-controls replays too.
    const assistOn = params.get("assist") === "1";
    simStartMs = performance.now();
    simClient.start(
      params.get("site") === "huffman"
        ? huffmanScenario(missionSeed, mode, member)
        : dec17Scenario(missionSeed, mode, member, assistOn),
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
      simStartMs = performance.now();
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
    // N: fresh relaunch of the SAME scenario — no ghost, no replay
    // verdict, a clean new run (B11 gameplay loop). Works mid-run too;
    // the old client is simply replaced.
    window.addEventListener("keydown", (e: KeyboardEvent) => {
      if (e.code !== "KeyN" || simClient === undefined) {
        return;
      }
      resultsCardEl.style.display = "none";
      simClient.dispose();
      simClient = makeClient();
      simStartMs = performance.now();
      simClient.start(
        params.get("site") === "huffman"
          ? huffmanScenario(missionSeed, mode, member)
          : dec17Scenario(missionSeed, mode, member, assistOn),
      );
      renderer.dispose();
      renderer = createFlyerSceneRenderer(app, simClient);
      resize();
      console.info(
        JSON.stringify({ suite: "wright-flyer-app", stage: "fresh-relaunch" }),
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
  // P toggles photo mode anywhere (HUD hide + plate grade, CSS-only).
  window.addEventListener("keydown", (e: KeyboardEvent) => {
    if (e.code === "KeyP" && !e.repeat) {
      const t = togglePhotoMode(document);
      console.info(
        JSON.stringify({ suite: "wright-flyer-app", stage: "photo-mode", active: t.active }),
      );
    }
  });
  // B10: FIELD INSTRUMENTS panel — the dormant visualization/lesson
  // modules behind one opt-in toggle row, fed by the latest snapshot.
  // The feed slot is declared FIRST: the panel registers its callback
  // synchronously during construction, so assigning it from the hook
  // would otherwise hit the temporal dead zone and kill main().
  let instrumentsFeed: ((s: InstrumentSample | null) => void) | null = null;
  const instruments = createInstrumentsPanel(document.body, {
    onFrame(cb) {
      instrumentsFeed = cb;
    },
  });


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
    // Clamp negatives: the FIRST rAF timestamp can precede the setup
    // capture (frame queued before module eval finished), and a
    // negative dt threw advanceProp's domain refusal, killing the
    // whole loop on frame 1. Presentation-plane only — the sim worker
    // never sees this number.
    const dtS = Math.min(0.25, Math.max(0, (nowMs - lastMs) / 1000));
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
    if (instrumentsFeed !== null && simClient !== undefined) {
      const s = simClient.sample(nowMs);
      if (s !== null) {
        instrumentsFeed({
          tick: s.tick,
          xM: s.xM,
          hM: s.hM,
          uMps: s.uMps,
          wMps: s.wMps,
          qRadS: s.qRadS,
          thetaRad: s.thetaRad,
          dcRad: s.dcRad,
          warpRad: s.warpRad,
          omegaRadS: s.omegaPropRadS,
          gustWMps: s.gustWMps,
          phase: s.phase,
        });
      } else {
        instrumentsFeed(null);
      }
    }
    // Latency ledger: publication (new sim tick visible) then present.
    const tick = simClient?.latestTick() ?? 0;
    if (tick > lastPublishedTick) {
      lastPublishedTick = tick;
      ledger.published(tick, nowMs);
    }
    ledger.presented(performance.now());
    // Startup watchdog (surfaced-failure law): tickless 6 s after a
    // start in an auto-stepping mode means the worker never came up —
    // say so on the honest banner instead of silently demoing.
    if (
      simClient !== undefined &&
      mode !== MODE_HUMAN &&
      !watchdogFired &&
      tick === 0 &&
      nowMs - simStartMs > 6000
    ) {
      watchdogFired = true;
      capabilityText.textContent =
        "SIM WORKER DID NOT START — showing the scripted attract demo. Reload to retry.";
      capabilityText.className = "warn";
      console.error(
        JSON.stringify({
          suite: "wright-flyer-app",
          stage: "sim-worker-startup-timeout",
          waitedMs: nowMs - simStartMs,
        }),
      );
    }
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
