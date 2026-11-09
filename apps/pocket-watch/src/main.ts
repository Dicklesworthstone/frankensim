import "@fontsource-variable/fraunces";
import "@fontsource-variable/inter";
import "@fontsource/jetbrains-mono/400.css";
import "@fontsource/jetbrains-mono/600.css";
import "./styles.css";

import initWasm, { PocketWatchSim } from "./wasm-pkg/fs_pocket_watch";
import { decodeWatchSnapshot, formatHorologicalTime, type WatchSnapshot } from "./snapshotLayout";
import { VirtualMicrophoneAudio } from "./audioEngine";
import { HorologicalTimegrapher, renderNodeLedgerTable } from "./timegrapher";
import { PocketWatchScene, type ActiveToolMode, type CameraPreset } from "./watchScene";

function $<T extends Element = HTMLElement>(sel: string): T {
  const el = document.querySelector<T>(sel);
  if (!el) throw new Error(`Missing DOM element: ${sel}`);
  return el;
}

async function boot(): Promise<void> {
  await initWasm();
  const sim = new PocketWatchSim();
  const micAudio = new VirtualMicrophoneAudio();
  const timegrapher = new HorologicalTimegrapher();

  // Interactive tool state
  let torchPower = 0.0;
  let torchPosMm: [number, number, number] = [0.0, -24.0, 2.0];
  let ln2Rate = 0.0;
  let ln2PosMm: [number, number, number] = [0.0, -16.0, 3.0];
  let oilRate = 0.0;
  let oilPosMm: [number, number, number] = [0.0, -11.0, 2.0];
  let magnetOn = false;
  let magnetDistMm = 28.0;
  let magnetPosMm: [number, number, number] = [28.0, -11.0, 2.5];
  let micPosMm: [number, number, number] = [0.0, -11.0, 11.0];
  let micGain = 1.4;
  let slowMo = false;

  const syncToolsToSim = (): void => {
    sim.set_blowtorch(torchPower, torchPosMm[0], torchPosMm[1], torchPosMm[2]);
    sim.set_liquid_nitrogen(ln2Rate, ln2PosMm[0], ln2PosMm[1], ln2PosMm[2]);
    sim.set_oil_pour(oilRate, oilPosMm[0], oilPosMm[1], oilPosMm[2]);
    sim.set_magnet(magnetOn, 1.0, magnetPosMm[0], magnetPosMm[1], magnetPosMm[2]);
    sim.set_microphone(micPosMm[0], micPosMm[1], micPosMm[2], micGain);
  };
  syncToolsToSim();

  // DOM References
  const canvasViewport = $<HTMLCanvasElement>("#viewport-canvas");
  const canvasScope = $<HTMLCanvasElement>("#canvas-oscilloscope");
  const canvasTimegrapher = $<HTMLCanvasElement>("#canvas-timegrapher");
  const tbodyLedger = $<HTMLElement>("#node-ledger-body");

  const hudWatchTime = $("#hud-watch-time");
  const hudLabTime = $("#hud-lab-time");
  const hudTimeDrift = $("#hud-time-drift");
  const hudRateSDay = $("#hud-rate-s-day");
  const hudStatusPill = $("#hud-status-pill");
  const hudDigest = $("#hud-digest");

  const hudPowerReservePct = $("#hud-power-reserve-pct");
  const hudPowerReserveBar = $("#hud-power-reserve-bar");
  const hudBarrelTorque = $("#hud-barrel-torque");
  const hudBalanceAmp = $("#hud-balance-amp");
  const hudPivotClearance = $("#hud-pivot-clearance");
  const hudTotalBeats = $("#hud-total-beats");

  const sliderExploded = $<HTMLInputElement>("#slider-exploded");
  const chkOpenCase = $<HTMLInputElement>("#chk-open-case");
  const chkInvarSpring = $<HTMLInputElement>("#chk-invar-spring");
  const chkSlowMo = $<HTMLInputElement>("#chk-slow-mo");
  const chkThermalOverlay = $<HTMLInputElement>("#chk-thermal-overlay");

  const sliderTorchPower = $<HTMLInputElement>("#slider-torch-power");
  const lblTorchPower = $("#lbl-torch-power");
  const sliderLn2Rate = $<HTMLInputElement>("#slider-ln2-rate");
  const lblLn2Rate = $("#lbl-ln2-rate");
  const sliderOilRate = $<HTMLInputElement>("#slider-oil-rate");
  const lblOilRate = $("#lbl-oil-rate");
  const chkMagnetOn = $<HTMLInputElement>("#chk-magnet-on");
  const sliderMagnetDist = $<HTMLInputElement>("#slider-magnet-dist");
  const lblBField = $("#lbl-b-field");
  const lblRemanence = $("#lbl-remanence");

  const btnToggleAudio = $<HTMLButtonElement>("#btn-toggle-audio");
  const sliderMicGain = $<HTMLInputElement>("#slider-mic-gain");
  const lblMicSpl = $("#lbl-mic-spl");

  const tgRate = $("#tg-rate");
  const tgAmp = $("#tg-amp");
  const tgBeatErr = $("#tg-beat-err");

  const ovCaseTemp = $("#ov-case-temp");
  const ovGoldMelt = $("#ov-gold-melt");
  const ovFouledGold = $("#ov-fouled-gold");
  const ovBField = $("#ov-b-field");

  const scene = new PocketWatchScene(canvasViewport, {
    onTargetMoved(tool, posMm) {
      if (tool === "torch") {
        torchPosMm = posMm;
        if (torchPower < 0.05) {
          torchPower = 0.75;
          sliderTorchPower.value = "0.75";
        }
      } else if (tool === "ln2") {
        ln2PosMm = posMm;
        if (ln2Rate < 0.05) {
          ln2Rate = 0.85;
          sliderLn2Rate.value = "0.85";
        }
      } else if (tool === "oil") {
        oilPosMm = posMm;
        if (oilRate < 0.05) {
          oilRate = 0.65;
          sliderOilRate.value = "0.65";
        }
      } else if (tool === "magnet") {
        magnetOn = true;
        chkMagnetOn.checked = true;
        magnetPosMm = posMm;
      } else if (tool === "mic") {
        micPosMm = posMm;
      }
      syncToolsToSim();
    },
  });

  // Bind Tool Buttons
  document.querySelectorAll<HTMLButtonElement>(".tool-btn").forEach((btn) => {
    btn.addEventListener("click", () => {
      document.querySelectorAll(".tool-btn").forEach((b) => b.classList.remove("active"));
      btn.classList.add("active");
      const mode = (btn.dataset.tool ?? "inspect") as ActiveToolMode;
      scene.setToolMode(mode);
    });
  });

  // Bind Camera Presets
  document.querySelectorAll<HTMLButtonElement>(".cam-btn").forEach((btn) => {
    btn.addEventListener("click", () => {
      document.querySelectorAll(".cam-btn").forEach((b) => b.classList.remove("active"));
      btn.classList.add("active");
      const preset = (btn.dataset.cam ?? "overview") as CameraPreset;
      scene.setCameraPreset(preset);
    });
  });

  // Crown Winding & Reset
  $("#btn-wind-crown").addEventListener("click", () => {
    sim.wind_crown(1.5);
  });

  const resetAllControls = (): void => {
    torchPower = 0.0;
    sliderTorchPower.value = "0";
    ln2Rate = 0.0;
    sliderLn2Rate.value = "0";
    oilRate = 0.0;
    sliderOilRate.value = "0";
    magnetOn = false;
    chkMagnetOn.checked = false;
    magnetDistMm = 28.0;
    sliderMagnetDist.value = "28";
    magnetPosMm = [28.0, -11.0, 2.5];
    syncToolsToSim();
  };

  $("#btn-reset-watch").addEventListener("click", () => {
    resetAllControls();
    sim.reset();
    sim.set_compensated_hairspring(chkInvarSpring.checked);
    timegrapher.reset();
  });

  // Architecture & View Toggles
  sliderExploded.addEventListener("input", () => {
    scene.setExplodedFactor(Number(sliderExploded.value));
  });
  chkOpenCase.addEventListener("change", () => {
    scene.setOpenCaseMode(chkOpenCase.checked);
  });
  chkInvarSpring.addEventListener("change", () => {
    sim.set_compensated_hairspring(chkInvarSpring.checked);
  });
  chkSlowMo.addEventListener("change", () => {
    slowMo = chkSlowMo.checked;
  });
  chkThermalOverlay.addEventListener("change", () => {
    scene.setThermalOverlayMode(chkThermalOverlay.checked);
  });

  // Perturbation Sliders
  sliderTorchPower.addEventListener("input", () => {
    torchPower = Number(sliderTorchPower.value);
    syncToolsToSim();
  });

  document.querySelectorAll<HTMLButtonElement>("[data-torch-aim]").forEach((chip) => {
    chip.addEventListener("click", () => {
      document.querySelectorAll("[data-torch-aim]").forEach((c) => c.classList.remove("active"));
      chip.classList.add("active");
      const parts = (chip.dataset.torchAim ?? "0,-24,2").split(",").map(Number);
      torchPosMm = [parts[0] ?? 0, parts[1] ?? -24, parts[2] ?? 2];
      syncToolsToSim();
    });
  });

  sliderLn2Rate.addEventListener("input", () => {
    ln2Rate = Number(sliderLn2Rate.value);
    syncToolsToSim();
  });

  sliderOilRate.addEventListener("input", () => {
    oilRate = Number(sliderOilRate.value);
    syncToolsToSim();
  });

  $("#btn-flood-oil").addEventListener("click", () => {
    oilPosMm = [0.0, -11.0, 2.0];
    sim.set_oil_pour(1.0, 0.0, -11.0, 2.0);
    for (let i = 0; i < 16; i++) {
      sim.step(0.05);
    }
    sim.set_oil_pour(oilRate, oilPosMm[0], oilPosMm[1], oilPosMm[2]);
  });

  $("#btn-clean-oil").addEventListener("click", () => {
    oilRate = 0.0;
    sliderOilRate.value = "0";
    sim.set_oil_pour(0.0, oilPosMm[0], oilPosMm[1], oilPosMm[2]);
    sim.clean_excess_oil();
  });

  chkMagnetOn.addEventListener("change", () => {
    magnetOn = chkMagnetOn.checked;
    magnetPosMm = [magnetDistMm, -11.0, 2.5];
    syncToolsToSim();
  });

  sliderMagnetDist.addEventListener("input", () => {
    magnetDistMm = Number(sliderMagnetDist.value);
    magnetPosMm = [magnetDistMm, -11.0, 2.5];
    syncToolsToSim();
  });

  $("#btn-degauss").addEventListener("click", () => {
    magnetOn = false;
    chkMagnetOn.checked = false;
    syncToolsToSim();
    sim.trigger_degauss();
  });

  // 3D Virtual Microphone Controls
  btnToggleAudio.addEventListener("click", async () => {
    const running = await micAudio.toggle(sim);
    btnToggleAudio.textContent = running ? "🔇 Mute 48 kHz Mic" : "🔊 Enable Live 48 kHz Mic";
  });

  sliderMicGain.addEventListener("input", () => {
    micGain = Number(sliderMicGain.value);
    syncToolsToSim();
  });

  document.querySelectorAll<HTMLButtonElement>("[data-mic-pos]").forEach((chip) => {
    chip.addEventListener("click", () => {
      document.querySelectorAll("[data-mic-pos]").forEach((c) => c.classList.remove("active"));
      chip.classList.add("active");
      const parts = (chip.dataset.micPos ?? "0,-11,11").split(",").map(Number);
      micPosMm = [parts[0] ?? 0, parts[1] ?? -11, parts[2] ?? 11];
      syncToolsToSim();
    });
  });

  // Guided Multiphysics Scenarios
  document.querySelectorAll<HTMLButtonElement>("[data-preset]").forEach((btn) => {
    btn.addEventListener("click", () => {
      const preset = btn.dataset.preset;
      switch (preset) {
        case "chronometer":
          resetAllControls();
          sim.reset();
          scene.setCameraPreset("overview");
          timegrapher.reset();
          break;

        case "thermal-drift":
          resetAllControls();
          sim.reset();
          torchPosMm = [0.0, -22.0, 2.0];
          torchPower = 0.42;
          sliderTorchPower.value = "0.42";
          syncToolsToSim();
          // Pre-warm for 1.5 s so thermal expansion drift is immediately visible
          for (let i = 0; i < 30; i++) sim.step(0.05);
          scene.setCameraPreset("tourbillon");
          break;

        case "melt-gold":
          resetAllControls();
          sim.reset();
          torchPosMm = [0.0, -25.0, 1.5];
          torchPower = 1.0;
          sliderTorchPower.value = "1.0";
          syncToolsToSim();
          // Heat the 6h 18k gold case into the mushy/liquid melting plateau so molten gold flows into the Tourbillon
          for (let i = 0; i < 125; i++) sim.step(0.05);
          scene.setCameraPreset("overview");
          break;

        case "melt-then-ln2": {
          // Ensure the gold case has melted into the gears first, then quench with LN2!
          const curSnap = decodeWatchSnapshot(sim.snapshot());
          if (curSnap.nodes[3]!.plasticSlump < 0.08) {
            sim.set_blowtorch(1.0, 0.0, -25.0, 1.5);
            for (let i = 0; i < 140; i++) sim.step(0.05);
          }
          torchPower = 0.0;
          sliderTorchPower.value = "0";
          ln2PosMm = [0.0, -18.0, 2.0];
          ln2Rate = 1.0;
          sliderLn2Rate.value = "1.0";
          syncToolsToSim();
          for (let i = 0; i < 90; i++) sim.step(0.05);
          scene.setCameraPreset("overview");
          break;
        }

        case "cryo-freeze":
          resetAllControls();
          sim.reset();
          ln2PosMm = [0.0, -11.0, 3.0];
          ln2Rate = 1.0;
          sliderLn2Rate.value = "1.0";
          syncToolsToSim();
          for (let i = 0; i < 85; i++) sim.step(0.05);
          scene.setCameraPreset("tourbillon");
          break;

        case "oil-flood":
          resetAllControls();
          sim.reset();
          oilPosMm = [0.0, -11.0, 2.0];
          oilRate = 0.85;
          sliderOilRate.value = "0.85";
          syncToolsToSim();
          for (let i = 0; i < 25; i++) sim.step(0.05);
          scene.setCameraPreset("tourbillon");
          break;

        case "magnet-attack":
          resetAllControls();
          magnetOn = true;
          chkMagnetOn.checked = true;
          magnetDistMm = 13.0;
          sliderMagnetDist.value = "13";
          magnetPosMm = [13.0, -11.0, 2.5];
          syncToolsToSim();
          scene.setCameraPreset("tourbillon");
          break;
      }
    });
  });

  // Main 60 FPS Simulation & Render Loop
  let lastFrameMs = performance.now();
  let frameCounter = 0;

  const updateHud = (snap: WatchSnapshot): void => {
    hudWatchTime.textContent = formatHorologicalTime(snap.watchDisplayedTimeS);
    hudLabTime.textContent = formatHorologicalTime(snap.labTimeS);
    const driftS = snap.watchDisplayedTimeS - snap.labTimeS;
    hudTimeDrift.textContent = `${driftS >= 0 ? "+" : ""}${driftS.toFixed(2)} s`;
    hudTimeDrift.style.color =
      Math.abs(driftS) < 0.25 ? "#67e8c8" : driftS < 0 ? "#ffc168" : "#ff7a8e";

    const effRate = 0.35 * snap.measuredRateSPerDay + 0.65 * snap.theoreticalRateSPerDay;
    const rateStr = `${effRate >= 0 ? "+" : ""}${effRate.toFixed(1)} s/d`;
    hudRateSDay.textContent = rateStr;
    tgRate.textContent = rateStr;
    tgAmp.textContent = `${snap.balanceAmplitudeDeg.toFixed(0)}°`;
    tgBeatErr.textContent = `${snap.measuredBeatErrorMs.toFixed(1)} ms`;

    const statusText = sim.status_label();
    hudStatusPill.textContent = statusText;
    hudStatusPill.className =
      "status-pill " +
      (snap.statusCode === 0
        ? "status-ok"
        : snap.statusCode === 1
          ? "status-warn"
          : "status-danger");

    if (frameCounter % 12 === 0) {
      hudDigest.textContent = sim.state_digest_hex().slice(0, 16) + "…";
    }

    const reservePct = Math.round(snap.mainspringWind * 100);
    const reserveHours = (snap.mainspringWind * 42.0).toFixed(1);
    hudPowerReservePct.textContent = `${reservePct}% · ${reserveHours} h`;
    hudPowerReserveBar.style.width = `${reservePct}%`;
    hudBarrelTorque.textContent = `${(snap.mainspringTorqueNm * 1000).toFixed(2)} mN·m`;
    hudBalanceAmp.textContent = `${snap.balanceAmplitudeDeg.toFixed(1)}°`;
    hudPivotClearance.textContent = `${snap.minPivotClearanceUm.toFixed(2)} µm`;
    hudTotalBeats.textContent = String(Math.floor(snap.totalBeats));

    lblTorchPower.textContent = `${Math.round(torchPower * 100)}% (${snap.blowtorchHeatingPowerW.toFixed(0)} W)`;
    lblLn2Rate.textContent = `${Math.round(ln2Rate * 100)}% (${snap.ln2BoilingPowerW.toFixed(0)} W)`;
    lblOilRate.textContent = `${snap.nodes[9]!.oilFilmUm.toFixed(1)} µm film`;
    lblBField.textContent = `${snap.bFieldEscapementMt.toFixed(1)} mT`;
    lblRemanence.textContent = `${snap.remanentMagnetizationMt.toFixed(1)} mT`;
    lblMicSpl.textContent = `${micAudio.lastRmsSplDb.toFixed(1)} dB`;

    // Stage bottom overlay
    const n6h = snap.nodes[3]!;
    ovCaseTemp.textContent = `${n6h.temperatureK.toFixed(1)} K (${(n6h.temperatureK - 273.15).toFixed(1)} °C)`;
    ovGoldMelt.textContent = `${(n6h.meltFraction * 100).toFixed(1)}% (slump ${(n6h.plasticSlump * 100).toFixed(0)}%)`;

    let totalLiqAu = 0;
    let totalSolidAu = 0;
    for (let i = 4; i < 12; i++) {
      totalLiqAu += snap.nodes[i]!.moltenGoldG;
      totalSolidAu += snap.nodes[i]!.frozenGoldG;
    }
    ovFouledGold.textContent = `${totalLiqAu.toFixed(3)} g / ${totalSolidAu.toFixed(3)} g`;
    ovBField.textContent = `${snap.bFieldEscapementMt.toFixed(1)} mT (rem ${snap.remanentMagnetizationMt.toFixed(1)} mT)`;

    if (frameCounter % 3 === 0) {
      renderNodeLedgerTable(tbodyLedger, snap);
    }
  };

  const tick = (nowMs: number): void => {
    const rawDtS = Math.min(0.05, Math.max(0.001, (nowMs - lastFrameMs) / 1000));
    lastFrameMs = nowMs;
    frameCounter++;

    const simDt = slowMo ? rawDtS * 0.15 : rawDtS;
    sim.step(simDt);
    micAudio.pumpSilentScopeBlock(sim);

    const snap = decodeWatchSnapshot(sim.snapshot());
    updateHud(snap);
    micAudio.renderOscilloscope(canvasScope, micPosMm);
    timegrapher.update(snap, canvasTimegrapher);
    scene.renderFrame(snap, {
      torchPower,
      torchPosMm,
      ln2Rate,
      ln2PosMm,
      oilRate,
      oilPosMm,
      magnetOn,
      magnetPosMm,
      micPosMm,
      caseShellData: sim.case_shell_snapshot(),
      spatialGridData: sim.spatial_grid_snapshot(),
      mainspringCoilsMm: sim.mainspring_coils_mm(),
      clickLeverAngleRad: sim.click_lever_angle_rad(),
      crownWheelAngleRad: sim.crown_wheel_angle_rad(),
    });

    requestAnimationFrame(tick);
  };

  requestAnimationFrame(tick);
}

void boot();
