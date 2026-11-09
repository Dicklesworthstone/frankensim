import "@fontsource-variable/fraunces";
import "@fontsource-variable/inter";
import "@fontsource/jetbrains-mono/400.css";
import "@fontsource/jetbrains-mono/600.css";
import "./styles.css";

import initWasm, { PocketWatchSim } from "./wasm-pkg/fs_pocket_watch";
import { decodeWatchSnapshot, formatHorologicalTime, type WatchSnapshot } from "./snapshotLayout";
import { VirtualMicrophoneAudio } from "./audioEngine";
import { HorologicalTimegrapher, renderNodeLedgerTable } from "./timegrapher";
import {
  PocketWatchScene,
  type ActiveToolMode,
  type CameraPreset,
  type VectorFieldLayer,
} from "./watchScene";

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

  // Interactive tool & simulation state
  let activeToolMode: ActiveToolMode = "inspect";
  let torchPower = 0.0;
  let torchPosMm: [number, number, number] = [0.0, -24.0, 2.0];
  let ln2Rate = 0.0;
  let ln2PosMm: [number, number, number] = [0.0, -16.0, 3.0];
  let oilRate = 0.0;
  let oilPosMm: [number, number, number] = [0.0, -11.0, 2.0];
  let magnetOn = false;
  let magnetDistMm = 28.0;
  let magnetPosMm: [number, number, number] = [28.0, -11.0, 2.5];
  let airPower = 0.0;
  let airPosMm: [number, number, number] = [0.0, -20.0, 3.0];
  let ambientPressureAtm = 1.0;
  let gravityGVec: [number, number, number] = [0.0, -1.0, 0.0];
  let gravityLabel = "1.0G 6h";
  let micPosMm: [number, number, number] = [0.0, -11.0, 11.0];
  let micGain = 1.4;
  let slowMo = false;
  let simPaused = false;

  // DOM References
  const appShell = $("#app-shell");
  const leftDock = $("#left-dock");
  const rightDock = $("#right-dock");
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
  const btnPauseSim = $<HTMLButtonElement>("#btn-pause-sim");
  const btnTheaterMode = $<HTMLButtonElement>("#btn-theater-mode");
  const btnStopPerturbations = $<HTMLButtonElement>("#btn-stop-perturbations");

  const stageToolBanner = $("#stage-tool-banner");
  const stageBannerText = $("#stage-banner-text");
  const btnBannerOrbit = $<HTMLButtonElement>("#btn-banner-orbit");
  const btnBannerStop = $<HTMLButtonElement>("#btn-banner-stop");

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
  const sliderAirPower = $<HTMLInputElement>("#slider-air-power");
  const lblAirPower = $("#lbl-air-power");
  const sliderPressureAtm = $<HTMLInputElement>("#slider-pressure-atm");
  const lblChamberEnv = $("#lbl-chamber-env");
  const btnShock500g = $<HTMLButtonElement>("#btn-shock-500g");

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

  const updatePerturbationBannerAndBadges = (): void => {
    const anyPerturbActive =
      torchPower > 0.01 ||
      ln2Rate > 0.01 ||
      oilRate > 0.01 ||
      magnetOn ||
      airPower > 0.01 ||
      Math.abs(ambientPressureAtm - 1.0) > 0.05;
    btnStopPerturbations.classList.toggle("hidden", !anyPerturbActive);

    // Highlight active subtool parameter box
    document.querySelectorAll(".subtool-box").forEach((b) => b.classList.remove("is-focused"));
    if (activeToolMode === "torch" || torchPower > 0.01) {
      $("#subtool-torch").classList.add("is-focused");
    }
    if (activeToolMode === "ln2" || ln2Rate > 0.01) {
      $("#subtool-ln2").classList.add("is-focused");
    }
    if (activeToolMode === "oil" || oilRate > 0.01) {
      $("#subtool-oil").classList.add("is-focused");
    }
    if (activeToolMode === "magnet" || magnetOn) {
      $("#subtool-magnet").classList.add("is-focused");
    }
    if (activeToolMode === "air" || airPower > 0.01) {
      $("#subtool-air").classList.add("is-focused");
    }

    canvasViewport.classList.toggle("tool-cursor", activeToolMode !== "inspect");

    if (activeToolMode !== "inspect" || anyPerturbActive) {
      stageToolBanner.classList.remove("hidden");
      btnBannerOrbit.classList.toggle("hidden", activeToolMode === "inspect");
      btnBannerStop.classList.toggle("hidden", !anyPerturbActive);

      if (activeToolMode === "torch") {
        stageBannerText.textContent = `🔥 Blowtorch Engaged (${Math.round(torchPower * 100)}%) — Drag across the 3D watch to aim the 2250 K flame jet`;
      } else if (activeToolMode === "ln2") {
        stageBannerText.textContent = `❄️ Cryogenic LN₂ Engaged (${Math.round(ln2Rate * 100)}%) — Drag across the 3D watch to pour 77 K liquid nitrogen`;
      } else if (activeToolMode === "oil") {
        stageBannerText.textContent = `💧 Precision Oiler Engaged (${Math.round(oilRate * 100)}%) — Drag across the 3D watch to deposit Pennzane oil`;
      } else if (activeToolMode === "magnet") {
        stageBannerText.textContent = `🧲 NdFeB Super-Magnet Engaged — Drag across the 3D watch to sweep the 1.32 T magnetic field`;
      } else if (activeToolMode === "air") {
        stageBannerText.textContent = `💨 Air Shear Nozzle Engaged (${Math.round(airPower * 100)}%) — Drag across the 3D watch to blow molten gold & oil`;
      } else if (activeToolMode === "mic") {
        stageBannerText.textContent = `🎙️ 3D Acoustic Probe Mode — Drag across the 3D watch to position the condenser microphone`;
      } else {
        const activeNames: string[] = [];
        if (torchPower > 0.01) activeNames.push(`🔥 Blowtorch ${Math.round(torchPower * 100)}%`);
        if (ln2Rate > 0.01) activeNames.push(`❄️ LN₂ ${Math.round(ln2Rate * 100)}%`);
        if (oilRate > 0.01) activeNames.push(`💧 Oil ${Math.round(oilRate * 100)}%`);
        if (magnetOn) activeNames.push(`🧲 Magnet (${magnetDistMm.toFixed(0)} mm)`);
        if (airPower > 0.01) activeNames.push(`💨 Air Jet ${Math.round(airPower * 100)}%`);
        if (Math.abs(ambientPressureAtm - 1.0) > 0.05) {
          activeNames.push(`🌍 ${ambientPressureAtm.toFixed(2)} atm`);
        }
        stageBannerText.textContent = `Active Perturbation: ${activeNames.join(" + ")} · Orbiting 3D camera`;
      }
    } else {
      stageToolBanner.classList.add("hidden");
    }
  };

  const syncToolsToSim = (): void => {
    sim.set_blowtorch(torchPower, torchPosMm[0], torchPosMm[1], torchPosMm[2]);
    sim.set_liquid_nitrogen(ln2Rate, ln2PosMm[0], ln2PosMm[1], ln2PosMm[2]);
    sim.set_oil_pour(oilRate, oilPosMm[0], oilPosMm[1], oilPosMm[2]);
    sim.set_magnet(magnetOn, 1.0, magnetPosMm[0], magnetPosMm[1], magnetPosMm[2]);
    sim.set_air_jet(airPower, airPosMm[0], airPosMm[1], airPosMm[2]);
    sim.set_ambient_pressure_atm(ambientPressureAtm);
    sim.set_gravity_vector(gravityGVec[0], gravityGVec[1], gravityGVec[2]);
    sim.set_microphone(micPosMm[0], micPosMm[1], micPosMm[2], micGain);
    updatePerturbationBannerAndBadges();
  };

  const scene = new PocketWatchScene(canvasViewport, {
    onTargetMoved(tool, posMm) {
      if (tool === "torch") {
        torchPosMm = posMm;
        if (torchPower < 0.05) {
          torchPower = 0.78;
          sliderTorchPower.value = "0.78";
        }
      } else if (tool === "ln2") {
        ln2PosMm = posMm;
        if (ln2Rate < 0.05) {
          ln2Rate = 0.80;
          sliderLn2Rate.value = "0.80";
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
        magnetDistMm = Math.round(
          Math.max(8, Math.min(80, Math.hypot(posMm[0] + 2.2, posMm[1] + 11.0))),
        );
        sliderMagnetDist.value = String(magnetDistMm);
      } else if (tool === "air") {
        airPosMm = posMm;
        if (airPower < 0.05) {
          airPower = 0.85;
          sliderAirPower.value = "0.85";
        }
      } else if (tool === "mic") {
        micPosMm = posMm;
      }
      syncToolsToSim();
    },
  });
  syncToolsToSim();

  const ALL_VEC_LAYERS: VectorFieldLayer[] = ["heat", "convection", "magnetic", "acoustic", "stress"];
  const btnVecAll = $<HTMLButtonElement>("#btn-vec-all");

  const syncVectorFieldButtons = (): void => {
    let allOn = true;
    document.querySelectorAll<HTMLButtonElement>("[data-vec-field]").forEach((btn) => {
      const layer = btn.dataset.vecField as VectorFieldLayer;
      const isOn = scene.isVectorFieldVisible(layer);
      btn.classList.toggle("active", isOn);
      if (!isOn) allOn = false;
    });
    btnVecAll.classList.toggle("active", allOn);
  };

  const setVectorField = (layer: VectorFieldLayer, visible: boolean): void => {
    scene.setVectorFieldVisible(layer, visible);
    syncVectorFieldButtons();
  };

  document.querySelectorAll<HTMLButtonElement>("[data-vec-field]").forEach((btn) => {
    btn.addEventListener("click", () => {
      const layer = btn.dataset.vecField as VectorFieldLayer;
      scene.toggleVectorField(layer);
      syncVectorFieldButtons();
    });
  });

  btnVecAll.addEventListener("click", () => {
    const anyOff = ALL_VEC_LAYERS.some((l) => !scene.isVectorFieldVisible(l));
    for (const l of ALL_VEC_LAYERS) {
      scene.setVectorFieldVisible(l, anyOff);
    }
    syncVectorFieldButtons();
  });

  const selectToolMode = (mode: ActiveToolMode, autoArm = true): void => {
    activeToolMode = mode;
    document.querySelectorAll<HTMLButtonElement>(".tool-btn").forEach((b) => {
      b.classList.toggle("active", b.dataset.tool === mode);
    });
    scene.setToolMode(mode);

    // Auto-arm a sensible working intensity when selecting a perturbation tool
    // so clicking the tool button immediately activates the physical apparatus!
    if (autoArm) {
      if (mode === "torch" && torchPower < 0.05) {
        torchPower = 0.78;
        sliderTorchPower.value = "0.78";
      } else if (mode === "ln2" && ln2Rate < 0.05) {
        ln2Rate = 0.80;
        sliderLn2Rate.value = "0.80";
      } else if (mode === "oil" && oilRate < 0.05) {
        oilRate = 0.65;
        sliderOilRate.value = "0.65";
      } else if (mode === "magnet" && !magnetOn) {
        magnetOn = true;
        chkMagnetOn.checked = true;
      } else if (mode === "air" && airPower < 0.05) {
        airPower = 0.85;
        sliderAirPower.value = "0.85";
      }
    }
    syncToolsToSim();
  };

  const applyCameraPreset = (preset: CameraPreset): void => {
    document.querySelectorAll<HTMLButtonElement>(".cam-btn").forEach((b) => {
      b.classList.toggle("active", b.dataset.cam === preset);
    });
    scene.setCameraPreset(preset);
  };

  const stopAllPerturbations = (): void => {
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
    airPower = 0.0;
    sliderAirPower.value = "0";
    ambientPressureAtm = 1.0;
    sliderPressureAtm.value = "1";
    gravityGVec = [0.0, -1.0, 0.0];
    gravityLabel = "1.0G 6h";
    document.querySelectorAll("[data-gravity-vec]").forEach((c, idx) => {
      c.classList.toggle("active", idx === 0);
    });
    selectToolMode("inspect", false);
  };

  // Bind Tool Selector Buttons
  document.querySelectorAll<HTMLButtonElement>(".tool-btn").forEach((btn) => {
    btn.addEventListener("click", () => {
      const mode = (btn.dataset.tool ?? "inspect") as ActiveToolMode;
      selectToolMode(mode, true);
    });
  });

  btnStopPerturbations.addEventListener("click", stopAllPerturbations);
  btnBannerStop.addEventListener("click", stopAllPerturbations);
  btnBannerOrbit.addEventListener("click", () => selectToolMode("inspect", false));

  // Bind Camera Presets
  document.querySelectorAll<HTMLButtonElement>(".cam-btn").forEach((btn) => {
    btn.addEventListener("click", () => {
      const preset = (btn.dataset.cam ?? "overview") as CameraPreset;
      applyCameraPreset(preset);
    });
  });

  // Pause / Resume & Theater Studio Mode
  const togglePause = (): void => {
    simPaused = !simPaused;
    btnPauseSim.textContent = simPaused ? "▶ Resume" : "⏸ Pause";
    btnPauseSim.classList.toggle("active", simPaused);
  };
  btnPauseSim.addEventListener("click", togglePause);

  const toggleTheaterMode = (): void => {
    const isTheater = appShell.classList.toggle("theater-mode");
    btnTheaterMode.classList.toggle("active", isTheater);
    btnTheaterMode.textContent = isTheater ? "⤡ Exit Studio" : "⤢ Studio View";
    setTimeout(() => scene.resize(), 60);
    setTimeout(() => scene.resize(), 240);
  };
  btnTheaterMode.addEventListener("click", toggleTheaterMode);

  // Mobile Segmented Tab Navigation
  const setMobileTab = (tabName: string): void => {
    document.querySelectorAll<HTMLButtonElement>(".mob-tab-btn").forEach((b) => {
      b.classList.toggle("active", b.dataset.mobileTab === tabName);
    });
    document.querySelectorAll<HTMLElement>(".panel-card[data-mobile-group]").forEach((card) => {
      card.classList.toggle("mobile-tab-active", card.dataset.mobileGroup === tabName);
    });
    const isMetrology = tabName === "metrology";
    leftDock.classList.toggle("mobile-dock-hidden", isMetrology);
    rightDock.classList.toggle("mobile-dock-hidden", !isMetrology);
    setTimeout(() => scene.resize(), 40);
  };
  setMobileTab("perturb");

  document.querySelectorAll<HTMLButtonElement>(".mob-tab-btn").forEach((btn) => {
    btn.addEventListener("click", () => {
      setMobileTab(btn.dataset.mobileTab ?? "perturb");
    });
  });

  // Crown Winding & Reset
  $("#btn-wind-crown").addEventListener("click", () => {
    sim.wind_crown(1.5);
  });

  const resetWatchAll = (): void => {
    stopAllPerturbations();
    sim.reset();
    sim.set_compensated_hairspring(chkInvarSpring.checked);
    timegrapher.reset();
  };

  $("#btn-reset-watch").addEventListener("click", resetWatchAll);

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
      if (torchPower < 0.05) {
        torchPower = 0.78;
        sliderTorchPower.value = "0.78";
      }
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
    updatePerturbationBannerAndBadges();
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

  // Compressed Air Shear Nozzle, Bell-Jar Pressure, 3D Gravity & 500G Shock Controls
  sliderAirPower.addEventListener("input", () => {
    airPower = Number(sliderAirPower.value);
    syncToolsToSim();
  });

  sliderPressureAtm.addEventListener("input", () => {
    ambientPressureAtm = Number(sliderPressureAtm.value);
    syncToolsToSim();
  });

  btnShock500g.addEventListener("click", () => {
    sim.trigger_mechanical_shock(500.0);
  });

  document.querySelectorAll<HTMLButtonElement>("[data-gravity-vec]").forEach((chip) => {
    chip.addEventListener("click", () => {
      document.querySelectorAll("[data-gravity-vec]").forEach((c) => c.classList.remove("active"));
      chip.classList.add("active");
      const raw = (chip.dataset.gravityVec ?? "0,-1,0,1.0G 6h").split(",");
      gravityGVec = [Number(raw[0] ?? 0), Number(raw[1] ?? -1), Number(raw[2] ?? 0)];
      gravityLabel = raw[3] ?? "1.0G";
      syncToolsToSim();
    });
  });

  // 3D Virtual Microphone Controls
  btnToggleAudio.addEventListener("click", async () => {
    const running = await micAudio.toggle(sim);
    btnToggleAudio.textContent = running ? "🔇 Mute 48 kHz Mic" : "🔊 Enable 48 kHz Mic";
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
      if (simPaused) togglePause();
      document.querySelectorAll(".preset-btn").forEach((b) => b.classList.remove("active"));
      btn.classList.add("active");
      const preset = btn.dataset.preset;
      switch (preset) {
        case "chronometer":
          resetWatchAll();
          applyCameraPreset("overview");
          break;

        case "thermal-drift":
          stopAllPerturbations();
          sim.reset();
          torchPosMm = [0.0, -22.0, 2.0];
          torchPower = 0.42;
          sliderTorchPower.value = "0.42";
          syncToolsToSim();
          for (let i = 0; i < 30; i++) sim.step(0.05);
          applyCameraPreset("tourbillon");
          break;

        case "melt-gold":
          stopAllPerturbations();
          sim.reset();
          torchPosMm = [0.0, -25.0, 1.5];
          torchPower = 1.0;
          sliderTorchPower.value = "1.0";
          syncToolsToSim();
          for (let i = 0; i < 125; i++) sim.step(0.05);
          applyCameraPreset("overview");
          break;

        case "melt-then-ln2": {
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
          applyCameraPreset("overview");
          break;
        }

        case "cryo-freeze":
          stopAllPerturbations();
          sim.reset();
          ln2PosMm = [0.0, -11.0, 3.0];
          ln2Rate = 1.0;
          sliderLn2Rate.value = "1.0";
          syncToolsToSim();
          for (let i = 0; i < 85; i++) sim.step(0.05);
          applyCameraPreset("tourbillon");
          break;

        case "oil-flood":
          stopAllPerturbations();
          sim.reset();
          oilPosMm = [0.0, -11.0, 2.0];
          oilRate = 0.85;
          sliderOilRate.value = "0.85";
          syncToolsToSim();
          for (let i = 0; i < 25; i++) sim.step(0.05);
          applyCameraPreset("tourbillon");
          break;

        case "magnet-attack":
          stopAllPerturbations();
          magnetOn = true;
          chkMagnetOn.checked = true;
          magnetDistMm = 13.0;
          sliderMagnetDist.value = "13";
          magnetPosMm = [13.0, -11.0, 2.5];
          syncToolsToSim();
          setVectorField("magnetic", true);
          applyCameraPreset("tourbillon");
          break;

        case "air-blow-gold":
          stopAllPerturbations();
          sim.reset();
          sim.set_blowtorch(1.0, 0.0, -25.0, 1.5);
          for (let i = 0; i < 120; i++) sim.step(0.05);
          sim.set_blowtorch(0.0, 0.0, -25.0, 1.5);
          airPosMm = [0.0, -24.0, 2.5];
          airPower = 0.92;
          sliderAirPower.value = "0.92";
          syncToolsToSim();
          for (let i = 0; i < 30; i++) sim.step(0.05);
          setVectorField("convection", true);
          applyCameraPreset("overview");
          break;

        case "vacuum-chamber":
          stopAllPerturbations();
          sim.reset();
          ambientPressureAtm = 0.0;
          sliderPressureAtm.value = "0";
          syncToolsToSim();
          for (let i = 0; i < 45; i++) sim.step(0.05);
          setVectorField("convection", true);
          applyCameraPreset("tourbillon");
          break;

        case "shock-centrifuge":
          stopAllPerturbations();
          sim.reset();
          gravityGVec = [0.0, -4.0, 0.0];
          gravityLabel = "4.0G Centr";
          syncToolsToSim();
          sim.trigger_mechanical_shock(500.0);
          for (let i = 0; i < 15; i++) sim.step(0.05);
          setVectorField("stress", true);
          applyCameraPreset("overview");
          break;
      }
    });
  });

  // Global Keyboard Shortcuts (guarded against inputs)
  window.addEventListener("keydown", (ev) => {
    const tag = (ev.target as HTMLElement | null)?.tagName;
    if (tag === "INPUT" || tag === "TEXTAREA" || tag === "SELECT") return;
    if (ev.metaKey || ev.ctrlKey || ev.altKey) return;

    const key = ev.key.toLowerCase();
    if (key === "1" || ev.key === "Escape") {
      selectToolMode("inspect", false);
    } else if (key === "2") {
      selectToolMode("torch", true);
    } else if (key === "3") {
      selectToolMode("ln2", true);
    } else if (key === "4") {
      selectToolMode("oil", true);
    } else if (key === "5") {
      selectToolMode("magnet", true);
    } else if (key === "6") {
      selectToolMode("air", true);
    } else if (key === "7") {
      selectToolMode("mic", false);
    } else if (key === "h") {
      scene.toggleVectorField("heat");
      syncVectorFieldButtons();
    } else if (key === "c") {
      scene.toggleVectorField("convection");
      syncVectorFieldButtons();
    } else if (key === "b") {
      scene.toggleVectorField("magnetic");
      syncVectorFieldButtons();
    } else if (key === "a") {
      scene.toggleVectorField("acoustic");
      syncVectorFieldButtons();
    } else if (key === "g") {
      scene.toggleVectorField("stress");
      syncVectorFieldButtons();
    } else if (key === "v") {
      btnVecAll.click();
    } else if (key === "s") {
      sim.trigger_mechanical_shock(500.0);
    } else if (key === "x") {
      stopAllPerturbations();
    } else if (key === "w") {
      sim.wind_crown(1.5);
    } else if (key === "r") {
      resetWatchAll();
    } else if (key === "f") {
      toggleTheaterMode();
    } else if (ev.code === "Space") {
      ev.preventDefault();
      togglePause();
    }
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
      Math.abs(driftS) < 0.25 ? "#52eba8" : driftS < 0 ? "#ffc56e" : "#ff8295";

    const effRate = 0.35 * snap.measuredRateSPerDay + 0.65 * snap.theoreticalRateSPerDay;
    const rateStr = `${effRate >= 0 ? "+" : ""}${effRate.toFixed(1)} s/d`;
    const dRateDt = sim.rate_sensitivity_s_day_per_k();
    hudRateSDay.textContent = `${rateStr} (${dRateDt.toFixed(1)} s/d/K)`;
    tgRate.textContent = rateStr;
    tgAmp.textContent = `${snap.balanceAmplitudeDeg.toFixed(0)}°`;
    tgBeatErr.textContent = `${snap.measuredBeatErrorMs.toFixed(1)} ms`;

    const statusText = simPaused ? "PAUSED — INSPECTING 3D CALIBRE" : sim.status_label();
    hudStatusPill.textContent = statusText;
    hudStatusPill.className =
      "status-pill " +
      (simPaused
        ? "status-warn"
        : snap.statusCode === 0
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
    lblAirPower.textContent = `${Math.round(airPower * 100)}% (${Math.round(airPower * 120)} m/s)`;
    lblChamberEnv.textContent = `${ambientPressureAtm.toFixed(2)} atm · ${gravityLabel}`;
    const splDb = sim.mic_spl_db();
    const sone = sim.mic_loudness_sone();
    const acum = sim.mic_sharpness_acum();
    lblMicSpl.textContent = `${splDb.toFixed(1)} dB · ${sone.toFixed(2)} sone · ${acum.toFixed(2)} acum`;

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

    if (!simPaused) {
      const simDt = slowMo ? rawDtS * 0.15 : rawDtS;
      sim.step(simDt);
      micAudio.pumpSilentScopeBlock(sim);
    }

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
      airPower,
      airPosMm,
      ambientPressureAtm,
      gravityGVec,
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
