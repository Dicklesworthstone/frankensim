// Main application controller for the FrankenSim Multiphysics Steam Turbine Simulation
// & L4 ASCENT Design Studio.
//
// Boots WebAssembly, binds all plant controls + 6 parametric machine hardware dimensions +
// live fs-dfo (CMA-ES + Nelder-Mead) optimizer, runs the 12-component fs-system convergence
// kernel + 17-stage IF97 meanline + 64×24 axisymmetric thermal/stress PDE, drives the HDR
// Chronometer-grade Three.js CAD viewport with 6 per-pixel GPU false-color shaders, 6 tapered-arrow
// 3D physics vector fields, interactive 3D surface probe, 48 kHz AudioWorklet, and Canvas2D metrology.

import "./styles.css";
import initWasm, { SteamTurbineSim } from "./wasm-pkg/fs_steam_turbine";
import {
  decodeDesignParams,
  decodeFemCertificate,
  decodeGeometryReceipt,
  decodeOptReport,
  decodeSnapshot,
  decodeStageProfile,
  STATUS_LABELS,
  type DesignOptimizationReport,
  type DesignParams,
  type FemShellCertificate,
  type GeometryReceipt,
  type StageTelemetry,
  type TurbineSnapshot,
} from "./snapshotLayout";
import { TurbineScene, type SurfaceProbeReading, type VectorFieldFlags } from "./turbineScene";
import { TurbineAudioEngine } from "./audioEngine";
import {
  drawCampbellDiagram,
  drawMollierDiagram,
  drawOptConvergenceChart,
  drawOrbitPlot,
  drawStageProfileChart,
  type StageChartMode,
} from "./metrology";
import { EXPERIMENTS, type CameraPreset, type CutawayMode, type ShaderMode } from "./experiments";
import { FIELD_MODE_SPECS } from "./turbineShaders";

function byId<T extends HTMLElement>(id: string): T {
  const el = document.getElementById(id);
  if (!el) throw new Error(`Missing DOM element #${id}`);
  return el as T;
}

function renderMetricGrid(
  container: HTMLElement,
  items: Array<{ k: string; v: string; tone?: "ok" | "warn" | "crit" }>,
): void {
  container.innerHTML = items
    .map(
      (it) =>
        `<div class="metric-tile ${it.tone ?? ""}"><span class="m-k">${it.k}</span><span class="m-v">${it.v}</span></div>`,
    )
    .join("");
}

async function boot(): Promise<void> {
  const bootOverlay = byId<HTMLDivElement>("boot-overlay");
  const bootStatus = byId<HTMLParagraphElement>("boot-status");

  bootStatus.textContent = "Loading FrankenSim WebAssembly module (26 crates · IF97 · FEM · fs-dfo)…";
  await initWasm();

  bootStatus.textContent = "Building AxisymmetricChart rotor, NurbsSurface blades, fs-plate FEM & 12-component graph…";
  const sim = new SteamTurbineSim();
  const geomReceipt: GeometryReceipt = decodeGeometryReceipt(sim.geometry_receipt_flat());
  const femCert: FemShellCertificate = decodeFemCertificate(sim.fem_shell_certificate_flat());
  let mollierDome = sim.mollier_dome();
  let campbellTable = sim.campbell_data();

  bootStatus.textContent = "Compiling HDR studio shaders & uploading 24 FrankenSim meshes + 34 NURBS rows…";
  const viewport = byId<HTMLDivElement>("viewport");
  const scene = new TurbineScene(viewport, sim);
  const audio = new TurbineAudioEngine();

  // State
  let running = true;
  let simRate = 1.0;
  let cutaway: CutawayMode = "half";
  let shader: ShaderMode = "pbr";
  let stageChartMode: StageChartMode = "pt";
  let activeRightTab: "thermo" | "rotor" | "audio" | "ledger" = "thermo";
  let activeExpId: string | null = "nominal";
  let selectedOptObj = 0;
  let lastOptReport: DesignOptimizationReport | null = null;
  let fastForwardRemainingS = 0;
  let frameCounter = 0;
  let lastWallMs = performance.now();

  // DOM References
  const kpiRpm = byId<HTMLSpanElement>("kpi-rpm");
  const kpiMw = byId<HTMLSpanElement>("kpi-mw");
  const kpiEff = byId<HTMLSpanElement>("kpi-eff");
  const kpiCond = byId<HTMLSpanElement>("kpi-cond");
  const kpiVib = byId<HTMLSpanElement>("kpi-vib");
  const kpiTemp = byId<HTMLSpanElement>("kpi-temp");
  const statusPill = byId<HTMLDivElement>("status-pill");
  const mkpiText = byId<HTMLSpanElement>("mkpi-text");
  const mkpiStatus = byId<HTMLSpanElement>("mkpi-status");

  const btnPlay = byId<HTMLButtonElement>("btn-play");
  const btnAudio = byId<HTMLButtonElement>("btn-audio");
  const btnReset = byId<HTMLButtonElement>("btn-reset");
  const btnBreaker = byId<HTMLButtonElement>("btn-breaker");
  const btnTrip = byId<HTMLButtonElement>("btn-trip");
  const btnResetTrip = byId<HTMLButtonElement>("btn-reset-trip");
  const btnColdStart = byId<HTMLButtonElement>("btn-cold-start");
  const btnQuickCutaway = byId<HTMLButtonElement>("btn-quick-cutaway");
  const btnQuickField = byId<HTMLButtonElement>("btn-quick-field");

  const expToast = byId<HTMLDivElement>("exp-toast");
  const expToastTag = byId<HTMLSpanElement>("exp-toast-tag");
  const expToastTitle = byId<HTMLElement>("exp-toast-title");
  const expToastDesc = byId<HTMLParagraphElement>("exp-toast-desc");

  const colorbarCard = byId<HTMLDivElement>("colorbar-card");
  const cbTitle = byId<HTMLSpanElement>("cb-title");
  const cbUnit = byId<HTMLSpanElement>("cb-unit");
  const cbBar = byId<HTMLDivElement>("cb-bar");
  const cbMin = byId<HTMLSpanElement>("cb-min");
  const cbMid = byId<HTMLSpanElement>("cb-mid");
  const cbMax = byId<HTMLSpanElement>("cb-max");

  const probeCard = byId<HTMLDivElement>("probe-card");
  const prPart = byId<HTMLElement>("pr-part");
  const prZR = byId<HTMLSpanElement>("pr-zr");
  const prGrid = byId<HTMLDivElement>("pr-grid");

  const canvasMollier = byId<HTMLCanvasElement>("canvas-mollier");
  const canvasStages = byId<HTMLCanvasElement>("canvas-stages");
  const canvasCampbell = byId<HTMLCanvasElement>("canvas-campbell");
  const canvasOrbit = byId<HTMLCanvasElement>("canvas-orbit");
  const canvasScope = byId<HTMLCanvasElement>("canvas-scope");
  const canvasFft = byId<HTMLCanvasElement>("canvas-fft");
  const canvasOpt = byId<HTMLCanvasElement>("canvas-opt");

  const badgeMoisture = byId<HTMLSpanElement>("badge-moisture");
  const badgeAbscissa = byId<HTMLSpanElement>("badge-abscissa");
  const badgeEproc = byId<HTMLSpanElement>("badge-eproc");
  const badgeSpl = byId<HTMLSpanElement>("badge-spl");
  const valDigest = byId<HTMLElement>("val-digest");

  const gridThermo = byId<HTMLDivElement>("grid-thermo");
  const gridRotor = byId<HTMLDivElement>("grid-rotor");
  const gridAudio = byId<HTMLDivElement>("grid-audio");
  const gridFemAudio = byId<HTMLDivElement>("grid-fem-audio");
  const gridLedger = byId<HTMLDivElement>("grid-ledger");
  const gridFem = byId<HTMLDivElement>("grid-fem");
  const gridGeom = byId<HTMLDivElement>("grid-geom");
  const gridThermal = byId<HTMLDivElement>("grid-thermal");

  // Populate static L2 MORPH Geometry Certificate & fs-plate::shell FEM Certificate once!
  renderMetricGrid(gridGeom, [
    { k: "Rotor Forging Mass", v: `${geomReceipt.shaftMassKg.toFixed(0)} kg`, tone: "ok" },
    { k: "17 Blade Rows Mass", v: `${geomReceipt.bladesMassKg.toFixed(0)} kg` },
    { k: "Total Rotor Mass", v: `${geomReceipt.totalMassKg.toFixed(0)} kg` },
    { k: "Exact Polar J_p", v: `${geomReceipt.polarInertiaKgM2.toFixed(1)} kg·m²` },
    { k: "Transverse J_d", v: `${geomReceipt.transverseInertiaKgM2.toFixed(0)} kg·m²` },
    { k: "Hausdorff Bound", v: `${(geomReceipt.rotorHausdorffM * 1000).toFixed(2)} mm`, tone: "ok" },
    {
      k: "Watertight / χ",
      v: `${geomReceipt.rotorWatertight ? "YES" : "NO"} · χ=${geomReceipt.rotorEulerChar}`,
      tone: geomReceipt.rotorWatertight ? "ok" : "crit",
    },
    { k: "Instanced Scene Tris", v: `${geomReceipt.totalSceneTriangles.toLocaleString()}` },
  ]);

  const femTiles: Array<{ k: string; v: string; tone?: "ok" | "warn" | "crit" }> = [
    {
      k: "Sylvester FEM Cert",
      v: femCert.certifiedFromFem ? "CERTIFIED (fs-plate)" : "FALLBACK",
      tone: femCert.certifiedFromFem ? "ok" : "warn",
    },
    { k: "Hood Mode f₀ / f₁", v: `${femCert.f0Hz.toFixed(0)} / ${femCert.f1Hz.toFixed(0)} Hz` },
    { k: "Hood Mode f₂", v: `${femCert.f2Hz.toFixed(0)} Hz` },
    { k: "Stiffness tr(K)", v: `${femCert.stiffnessTraceNM.toExponential(2)} N/m` },
    { k: "Levine-Schwinger R", v: `${femCert.diffuserRadResistance.toFixed(3)} (fs-duct)`, tone: "ok" },
    { k: "Shell Wall / Steel", v: "28 mm · S355J2" },
  ];
  renderMetricGrid(gridFem, femTiles);
  renderMetricGrid(gridFemAudio, femTiles.slice(0, 4));

  // ---------- Interactive 3D CAD Hover Probe Readout ----------
  scene.setProbeCallback((reading: SurfaceProbeReading | null) => {
    if (!reading) {
      probeCard.classList.add("hidden");
      return;
    }
    probeCard.classList.remove("hidden");
    prPart.textContent = reading.partName;
    prZR.textContent = `z = ${reading.zM.toFixed(2)} m · r = ${reading.rM.toFixed(2)} m`;
    prGrid.innerHTML = [
      ["Temp", `${reading.tempC.toFixed(1)} °C`],
      ["Stress σ", `${reading.vmMpa.toFixed(0)} MPa`],
      ["Press P", `${reading.pBar >= 1 ? reading.pBar.toFixed(1) : reading.pBar.toFixed(2)} bar`],
      ["Mach", `${reading.mach.toFixed(2)}`],
      ["Flux |q|", `${reading.qMagKwM2.toFixed(1)} kW/m²`],
      ["Wetness", `${reading.wetnessPct.toFixed(1)}%`],
    ]
      .map(([k, v]) => `<div class="probe-cell"><span>${k}</span><strong>${v}</strong></div>`)
      .join("");
  });

  // ---------- Tab & Mobile Drawer Wiring ----------
  const switchLeftTab = (tab: string) => {
    document.querySelectorAll<HTMLButtonElement>("[data-left-tab]").forEach((b) => {
      const active = b.dataset.leftTab === tab;
      b.classList.toggle("active", active);
      b.setAttribute("aria-selected", String(active));
    });
    document.querySelectorAll<HTMLElement>("#panel-left .tab-pane").forEach((p) => {
      p.classList.toggle("active", p.id === `left-${tab}`);
    });
    if (tab === "design") {
      drawOptConvergenceChart(canvasOpt, lastOptReport?.historyObj ?? []);
    }
  };

  let hudReady = false;
  const switchRightTab = (tab: "thermo" | "rotor" | "audio" | "ledger") => {
    activeRightTab = tab;
    document.querySelectorAll<HTMLButtonElement>("[data-right-tab]").forEach((b) => {
      const active = b.dataset.rightTab === tab;
      b.classList.toggle("active", active);
      b.setAttribute("aria-selected", String(active));
    });
    document.querySelectorAll<HTMLElement>("#panel-right .tab-pane").forEach((p) => {
      p.classList.toggle("active", p.id === `right-${tab}`);
    });
    if (tab === "audio") {
      audio.prefillIfNeeded(sim);
    }
    if (hudReady) {
      updateHud(decodeSnapshot(sim.snapshot()), decodeStageProfile(sim.stage_profile_flat()));
    }
  };

  document.querySelectorAll<HTMLButtonElement>("[data-left-tab]").forEach((b) => {
    b.addEventListener("click", () => switchLeftTab(b.dataset.leftTab ?? "controls"));
  });
  document.querySelectorAll<HTMLButtonElement>("[data-right-tab]").forEach((b) => {
    b.addEventListener("click", () => {
      switchRightTab((b.dataset.rightTab ?? "thermo") as "thermo" | "rotor" | "audio" | "ledger");
    });
  });

  const panelLeft = byId<HTMLElement>("panel-left");
  const panelRight = byId<HTMLElement>("panel-right");
  const appRoot = byId<HTMLDivElement>("app");
  byId<HTMLButtonElement>("btn-drawer-left").addEventListener("click", () => {
    panelLeft.classList.toggle("mobile-open");
    panelRight.classList.remove("mobile-open");
  });
  byId<HTMLButtonElement>("btn-drawer-right").addEventListener("click", () => {
    panelRight.classList.toggle("mobile-open");
    panelLeft.classList.remove("mobile-open");
  });
  // Tapping the 3D viewport on smartphones cleanly dismisses open drawers for full-screen touch orbiting!
  viewport.addEventListener("pointerdown", () => {
    if (panelLeft.classList.contains("mobile-open") || panelRight.classList.contains("mobile-open")) {
      panelLeft.classList.remove("mobile-open");
      panelRight.classList.remove("mobile-open");
    }
  });

  // ---------- Intuitive CAD Viewport Navigation Pill (Mouse + Trackpad + Smartphone Touch) ----------
  const btnNavOrbit = byId<HTMLButtonElement>("btn-nav-orbit");
  const btnNavPan = byId<HTMLButtonElement>("btn-nav-pan");
  const btnNavZoomIn = byId<HTMLButtonElement>("btn-nav-zoomin");
  const btnNavZoomOut = byId<HTMLButtonElement>("btn-nav-zoomout");
  const btnNavFit = byId<HTMLButtonElement>("btn-nav-fit");
  const btnZen = byId<HTMLButtonElement>("btn-zen");
  const navHint = byId<HTMLDivElement>("nav-hint");

  const syncNavHint = () => {
    const touchUi =
      window.innerWidth <= 1020 ||
      window.matchMedia("(pointer: coarse)").matches ||
      "ontouchstart" in window;
    navHint.textContent = touchUi
      ? "1 Finger: Rotate/Pan · 2 Fingers: Pinch Zoom & Pan · Dbl-tap: Focus"
      : "Drag: Orbit · Right/Shift+Drag: Pan · Scroll: Zoom to cursor · Dbl-click: Focus";
  };
  syncNavHint();
  window.addEventListener("resize", syncNavHint);
  scene.setFirstInteractCallback(() => {
    navHint.classList.add("faded");
  });

  const applyNavMode = (mode: "orbit" | "pan") => {
    scene.setNavMode(mode);
    btnNavOrbit.classList.toggle("active", mode === "orbit");
    btnNavPan.classList.toggle("active", mode === "pan");
  };
  btnNavOrbit.addEventListener("click", () => applyNavMode("orbit"));
  btnNavPan.addEventListener("click", () => applyNavMode("pan"));
  const bindZoomButton = (btn: HTMLButtonElement, dir: 1 | -1, tapFactor: number) => {
    let holdTimer: number | null = null;
    let holding = false;
    const stop = () => {
      if (holdTimer !== null) {
        window.clearTimeout(holdTimer);
        holdTimer = null;
      }
      scene.setZoomHold(0);
    };
    btn.addEventListener("pointerdown", (e) => {
      e.preventDefault();
      holding = false;
      holdTimer = window.setTimeout(() => {
        holding = true;
        scene.setZoomHold(dir);
      }, 170);
    });
    btn.addEventListener("pointerup", () => {
      const wasHolding = holding;
      stop();
      if (!wasHolding) scene.zoomByFactor(tapFactor);
    });
    btn.addEventListener("pointerleave", stop);
    btn.addEventListener("pointercancel", stop);
  };
  bindZoomButton(btnNavZoomIn, 1, 0.76);
  bindZoomButton(btnNavZoomOut, -1, 1.32);
  btnNavFit.addEventListener("click", () => selectCamera("overview"));
  btnZen.addEventListener("click", () => {
    const zen = appRoot.classList.toggle("zen-mode");
    btnZen.classList.toggle("active", zen);
    btnZen.textContent = zen ? "Show UI" : "Full 3D";
  });

  // ---------- Control Room Sliders & Buttons ----------
  const bindSlider = (
    sliderId: string,
    labelId: string,
    format: (v: number) => string,
    onChange: (v: number) => void,
  ) => {
    const sl = byId<HTMLInputElement>(sliderId);
    const lb = byId<HTMLElement>(labelId);
    const sync = () => {
      const v = Number(sl.value);
      lb.textContent = format(v);
      onChange(v);
    };
    sl.addEventListener("input", sync);
    return {
      get: () => Number(sl.value),
      set: (v: number) => {
        sl.value = String(v);
        lb.textContent = format(v);
      },
    };
  };

  const ctlLoad = bindSlider("sl-load", "val-load", (v) => `${v.toFixed(0)}%`, (v) => sim.set_load_setpoint(v / 100));
  const ctlExt = bindSlider("sl-ext", "val-ext", (v) => `${v.toFixed(0)}%`, (v) => sim.set_extraction_demand(v / 100));
  const ctlPBoiler = bindSlider("sl-pboiler", "val-pboiler", (v) => `${v.toFixed(0)} bar`, (v) => sim.set_live_steam_bar(v));
  const ctlTBoiler = bindSlider("sl-tboiler", "val-tboiler", (v) => `${v.toFixed(0)} °C`, (v) => sim.set_live_steam_temp_c(v));
  const ctlRpmRef = bindSlider("sl-rpmref", "val-rpmref", (v) => `${v.toFixed(0)} rpm`, (v) => sim.set_speed_ref_rpm(v));
  const ctlExc = bindSlider("sl-exc", "val-exc", (v) => `${(v / 100).toFixed(2)} pu`, (v) => sim.set_excitation_pu(v / 100));
  const ctlAir = bindSlider("sl-air", "val-air", (v) => `${v.toFixed(2)}%`, (v) => sim.set_condenser_air_fraction(v / 100));
  const ctlCwFlow = bindSlider("sl-cwflow", "val-cwflow", (v) => `${v.toFixed(0)}%`, (v) => sim.set_cw_flow_frac(v / 100));
  const ctlCwTemp = bindSlider("sl-cwtemp", "val-cwtemp", (v) => `${v.toFixed(1)} °C`, (v) => sim.set_cw_inlet_temp_c(v));
  const ctlUnbal = bindSlider("sl-unbal", "val-unbal", (v) => `${Math.round(v).toLocaleString()} g·mm`, (v) => sim.set_rotor_unbalance_g_mm(v));
  const ctlWhip = bindSlider("sl-whip", "val-whip", (v) => v.toFixed(2), (v) => {
    sim.set_oil_whip_destabilize(v);
    sim.refresh_campbell();
    campbellTable = sim.campbell_data();
  });
  const ctlOilP = bindSlider("sl-oilp", "val-oilp", (v) => `${v.toFixed(0)}%`, (v) => {
    sim.set_lube_oil_pressure_frac(v / 100);
    sim.refresh_campbell();
    campbellTable = sim.campbell_data();
  });

  // ---------- L4 ASCENT Design Studio Sliders & fs-dfo Optimizer ----------
  const btnDesApply = byId<HTMLButtonElement>("btn-des-apply");
  let activeDesign: DesignParams = decodeDesignParams(sim.design_params_flat());
  const updateDesignDirtyUi = () => {
    const dirty =
      Math.abs(desTip.get() - activeDesign.tipClearanceMm) > 0.012 ||
      Math.abs(desReact.get() - activeDesign.hpReaction) > 0.012 ||
      Math.abs(desLpBl.get() - activeDesign.lpLastBladeHM) > 0.002 ||
      Math.abs(desBrgC.get() - activeDesign.bearingClearanceUm) > 1.2 ||
      Math.abs(desPrel.get() - activeDesign.bearingPreload) > 0.012 ||
      Math.abs(desWall.get() - activeDesign.casingWallMm) > 1.2;
    btnDesApply.classList.toggle("state-off", dirty);
    btnDesApply.classList.toggle("state-on", !dirty);
    btnDesApply.textContent = dirty ? "Apply Geometry *" : "Apply Geometry";
  };

  const desTip = bindSlider("sl-des-tip", "val-des-tip", (v) => `${v.toFixed(2)} mm`, updateDesignDirtyUi);
  const desReact = bindSlider("sl-des-react", "val-des-react", (v) => v.toFixed(2), updateDesignDirtyUi);
  const desLpBl = bindSlider("sl-des-lpbl", "val-des-lpbl", (v) => `${v.toFixed(3)} m`, updateDesignDirtyUi);
  const desBrgC = bindSlider("sl-des-brgc", "val-des-brgc", (v) => `${v.toFixed(0)} µm`, updateDesignDirtyUi);
  const desPrel = bindSlider("sl-des-prel", "val-des-prel", (v) => v.toFixed(2), updateDesignDirtyUi);
  const desWall = bindSlider("sl-des-wall", "val-des-wall", (v) => `${v.toFixed(0)} mm`, updateDesignDirtyUi);

  const syncDesignSliders = (p: DesignParams) => {
    activeDesign = p;
    desTip.set(p.tipClearanceMm);
    desReact.set(p.hpReaction);
    desLpBl.set(p.lpLastBladeHM);
    desBrgC.set(p.bearingClearanceUm);
    desPrel.set(p.bearingPreload);
    desWall.set(p.casingWallMm);
    updateDesignDirtyUi();
  };
  syncDesignSliders(activeDesign);

  const applyDesignToEngine = (p: DesignParams) => {
    sim.apply_design_params(
      p.tipClearanceMm,
      p.hpReaction,
      p.lpLastBladeHM,
      p.bearingClearanceUm,
      p.bearingPreload,
      p.casingWallMm,
    );
    sim.refresh_campbell();
    campbellTable = sim.campbell_data();
    scene.rebuildBladeRows(sim);
    syncDesignSliders(decodeDesignParams(sim.design_params_flat()));
    updateHud(decodeSnapshot(sim.snapshot()), decodeStageProfile(sim.stage_profile_flat()));
  };

  btnDesApply.addEventListener("click", () => {
    applyDesignToEngine({
      tipClearanceMm: desTip.get(),
      hpReaction: desReact.get(),
      lpLastBladeHM: desLpBl.get(),
      bearingClearanceUm: desBrgC.get(),
      bearingPreload: desPrel.get(),
      casingWallMm: desWall.get(),
    });
  });

  byId<HTMLButtonElement>("btn-des-reset").addEventListener("click", () => {
    applyDesignToEngine({
      tipClearanceMm: 0.75,
      hpReaction: 0.48,
      lpLastBladeHM: 0.527,
      bearingClearanceUm: 220,
      bearingPreload: 0.32,
      casingWallMm: 95,
    });
  });

  document.querySelectorAll<HTMLButtonElement>("#seg-opt-obj button").forEach((b) => {
    b.addEventListener("click", () => {
      selectedOptObj = Number(b.dataset.obj ?? 0);
      document.querySelectorAll<HTMLButtonElement>("#seg-opt-obj button").forEach((x) => x.classList.remove("active"));
      b.classList.add("active");
    });
  });

  const btnRunOpt = byId<HTMLButtonElement>("btn-run-opt");
  const btnApplyOpt = byId<HTMLButtonElement>("btn-apply-opt");
  const optSummary = byId<HTMLDivElement>("opt-summary");

  const runDesignOpt = () => {
    btnRunOpt.disabled = true;
    btnRunOpt.textContent = "Optimizing 6D Space (fs-dfo)…";
    const t0 = performance.now();
    const raw = sim.run_design_optimization(selectedOptObj, 120);
    const ms = performance.now() - t0;
    const rep = decodeOptReport(raw);
    lastOptReport = rep;
    btnRunOpt.disabled = false;
    btnRunOpt.textContent = "Run fs-dfo CMA-ES + Nelder-Mead";
    btnApplyOpt.disabled = false;

    drawOptConvergenceChart(canvasOpt, rep.historyObj);

    const b = rep.baseline;
    const o = rep.optimized;
    const dElec = o.elecEffPct - b.elecEffPct;
    const dChp = o.chpEffPct - b.chpEffPct;
    const dMw = o.elecPowerMw - b.elecPowerMw;
    const dRe = o.spectralAbscissaPerS - b.spectralAbscissaPerS;
    const fmtD = (v: number, unit: string) => `${v >= 0 ? "+" : ""}${v.toFixed(2)}${unit}`;

    optSummary.innerHTML = `
      <div><strong>fs-dfo converged</strong> in ${rep.totalEvals} evals (${rep.cmaGenerations} CMA-ES gens + NM polish, ${ms.toFixed(0)} ms)</div>
      <table class="opt-delta-table">
        <thead><tr><th>Metric</th><th>Base</th><th>Opt</th><th>Δ</th></tr></thead>
        <tbody>
          <tr><td>Elec η</td><td>${b.elecEffPct.toFixed(2)}%</td><td>${o.elecEffPct.toFixed(2)}%</td><td class="${dElec >= 0 ? "good" : "warn"}">${fmtD(dElec, "%")}</td></tr>
          <tr><td>CHP η</td><td>${b.chpEffPct.toFixed(2)}%</td><td>${o.chpEffPct.toFixed(2)}%</td><td class="${dChp >= 0 ? "good" : "warn"}">${fmtD(dChp, "%")}</td></tr>
          <tr><td>Output</td><td>${b.elecPowerMw.toFixed(1)}</td><td>${o.elecPowerMw.toFixed(1)} MW</td><td class="${dMw >= 0 ? "good" : "warn"}">${fmtD(dMw, "")}</td></tr>
          <tr><td>Blade σ</td><td>${b.maxBladeStressMpa.toFixed(0)}</td><td>${o.maxBladeStressMpa.toFixed(0)} MPa</td><td class="${o.maxBladeStressMpa <= 480 ? "good" : "warn"}">SF ${o.minGoodmanSf.toFixed(2)}</td></tr>
          <tr><td>max Re(λ)</td><td>${b.spectralAbscissaPerS.toFixed(1)}</td><td>${o.spectralAbscissaPerS.toFixed(1)}/s</td><td class="${dRe <= 0 ? "good" : "warn"}">${fmtD(dRe, "")}</td></tr>
        </tbody>
      </table>
      <div style="margin-top:5px; font-family:var(--font-mono); font-size:9.5px; color:#86efac;">
        Tip ${rep.optimizedParams.tipClearanceMm.toFixed(2)}mm · R ${rep.optimizedParams.hpReaction.toFixed(2)} · LP ${rep.optimizedParams.lpLastBladeHM.toFixed(3)}m · Brg ${rep.optimizedParams.bearingClearanceUm.toFixed(0)}µm · Pre ${rep.optimizedParams.bearingPreload.toFixed(2)} · Wall ${rep.optimizedParams.casingWallMm.toFixed(0)}mm
      </div>
    `;
    return rep;
  };

  btnRunOpt.addEventListener("click", () => {
    runDesignOpt();
  });

  btnApplyOpt.addEventListener("click", () => {
    if (!lastOptReport) return;
    applyDesignToEngine(lastOptReport.optimizedParams);
  });

  // ---------- Visual & Mic Sliders ----------
  bindSlider("sl-explode", "val-explode", (v) => `${v.toFixed(0)}%`, (v) => scene.setExplodeFraction(v / 100));
  bindSlider("sl-slowmo", "val-slowmo", (v) => (v >= 99 ? "1× real" : `1/${Math.max(1, Math.round(100 / v))}×`), (v) =>
    scene.setVisualSlowMo(v / 100),
  );
  let micZ = 2.4;
  let micX = 2.2;
  bindSlider("sl-micz", "val-micz", (v) => `${v.toFixed(1)} m`, (v) => {
    micZ = v;
    sim.set_mic_pos(micX, 1.4, micZ);
    scene.setMicPosition(micX, 1.4, micZ);
  });
  bindSlider("sl-micx", "val-micx", (v) => `${v.toFixed(1)} m`, (v) => {
    micX = v;
    sim.set_mic_pos(micX, 1.4, micZ);
    scene.setMicPosition(micX, 1.4, micZ);
  });
  bindSlider("sl-gain", "val-gain", (v) => `${v.toFixed(1)}×`, (v) => sim.set_mic_gain(v));

  let breakerClosed = true;
  const syncBreakerBtn = (closed: boolean) => {
    breakerClosed = closed;
    btnBreaker.textContent = closed ? "11 kV Breaker: CLOSED" : "11 kV Breaker: OPEN";
    btnBreaker.classList.toggle("state-on", closed);
    btnBreaker.classList.toggle("state-off", !closed);
  };

  btnBreaker.addEventListener("click", () => {
    syncBreakerBtn(!breakerClosed);
    sim.set_grid_breaker(breakerClosed);
  });
  btnTrip.addEventListener("click", () => sim.set_manual_trip(true));
  btnResetTrip.addEventListener("click", () => {
    sim.set_manual_trip(false);
    sim.reset_trips();
  });

  const resetControlSlidersToNominal = () => {
    ctlLoad.set(78);
    ctlExt.set(40);
    ctlPBoiler.set(110);
    ctlTBoiler.set(525);
    ctlRpmRef.set(3000);
    ctlExc.set(105);
    ctlAir.set(0.05);
    ctlCwFlow.set(100);
    ctlCwTemp.set(20);
    ctlUnbal.set(14000);
    ctlWhip.set(0);
    ctlOilP.set(100);
    syncBreakerBtn(true);
  };

  btnColdStart.addEventListener("click", () => {
    fastForwardRemainingS = 0;
    sim.reset(true);
    resetControlSlidersToNominal();
    ctlLoad.set(20);
    ctlExt.set(0);
    syncBreakerBtn(false);
    sim.set_load_setpoint(0.20);
    sim.refresh_campbell();
    campbellTable = sim.campbell_data();
    updateHud(decodeSnapshot(sim.snapshot()), decodeStageProfile(sim.stage_profile_flat()));
  });

  btnReset.addEventListener("click", () => {
    fastForwardRemainingS = 0;
    activeExpId = "nominal";
    expToast.classList.add("hidden");
    sim.reset(false);
    resetControlSlidersToNominal();
    syncDesignSliders(decodeDesignParams(sim.design_params_flat()));
    scene.rebuildBladeRows(sim);
    sim.refresh_campbell();
    campbellTable = sim.campbell_data();
    mollierDome = sim.mollier_dome();
    renderExperimentCards();
    updateHud(decodeSnapshot(sim.snapshot()), decodeStageProfile(sim.stage_profile_flat()));
  });

  btnPlay.addEventListener("click", () => {
    running = !running;
    btnPlay.textContent = running ? "Pause" : "Resume";
    btnPlay.classList.toggle("primary", running);
  });

  btnAudio.addEventListener("click", async () => {
    const on = await audio.toggle(sim);
    btnAudio.textContent = on ? "Audio Live" : "Audio Off";
    btnAudio.classList.toggle("active-audio", on);
  });

  // ---------- Visual Modes, 6 GPU Shaders, Colorbar, 6 Vector Fields & Cameras ----------
  const applyCutaway = (m: CutawayMode) => {
    cutaway = m;
    scene.setCutawayMode(m);
    document.querySelectorAll<HTMLButtonElement>("#seg-cutaway button").forEach((b) => {
      b.classList.toggle("active", b.dataset.cutaway === m);
    });
    const label = m === "half" ? "Half" : m === "open" ? "Open" : "Closed";
    btnQuickCutaway.textContent = `Cutaway: ${label}`;
  };

  const SHADER_ORDER: ShaderMode[] = ["pbr", "temp", "stress", "pressure", "mach", "wetcreep"];
  const applyShader = (s: ShaderMode) => {
    shader = s;
    scene.setShaderMode(s);
    const spec = FIELD_MODE_SPECS[s];
    document.querySelectorAll<HTMLButtonElement>("#seg-shader button").forEach((b) => {
      b.classList.toggle("active", b.dataset.shader === s);
    });
    btnQuickField.textContent = `Field: ${spec.shortLabel}`;
    btnQuickField.classList.toggle("active", s !== "pbr");

    if (s === "pbr") {
      colorbarCard.classList.add("hidden");
    } else {
      colorbarCard.classList.remove("hidden");
      cbTitle.textContent = spec.label;
      cbUnit.textContent = spec.unit;
      cbBar.style.background = spec.gradientCss;
      cbMin.textContent = spec.minLabel;
      cbMid.textContent = spec.midLabel;
      cbMax.textContent = spec.maxLabel;
    }
  };

  document.querySelectorAll<HTMLButtonElement>("#seg-cutaway button").forEach((b) => {
    b.addEventListener("click", () => applyCutaway((b.dataset.cutaway ?? "half") as CutawayMode));
  });
  document.querySelectorAll<HTMLButtonElement>("#seg-shader button").forEach((b) => {
    b.addEventListener("click", () => applyShader((b.dataset.shader ?? "pbr") as ShaderMode));
  });
  document.querySelectorAll<HTMLButtonElement>("#seg-timescale button").forEach((b) => {
    b.addEventListener("click", () => {
      simRate = Number(b.dataset.rate ?? 1);
      document.querySelectorAll<HTMLButtonElement>("#seg-timescale button").forEach((x) => x.classList.remove("active"));
      b.classList.add("active");
    });
  });
  document.querySelectorAll<HTMLButtonElement>("#seg-stage-mode button").forEach((b) => {
    b.addEventListener("click", () => {
      stageChartMode = (b.dataset.stmode ?? "pt") as StageChartMode;
      document.querySelectorAll<HTMLButtonElement>("#seg-stage-mode button").forEach((x) => x.classList.remove("active"));
      b.classList.add("active");
      drawStageProfileChart(canvasStages, decodeStageProfile(sim.stage_profile_flat()), stageChartMode);
    });
  });

  btnQuickCutaway.addEventListener("click", () => {
    const order: CutawayMode[] = ["half", "open", "closed"];
    applyCutaway(order[(order.indexOf(cutaway) + 1) % order.length]!);
  });
  btnQuickField.addEventListener("click", () => {
    applyShader(SHADER_ORDER[(SHADER_ORDER.indexOf(shader) + 1) % SHADER_ORDER.length]!);
  });

  byId<HTMLInputElement>("chk-iso").addEventListener("change", (e) => {
    scene.setIsoContours((e.target as HTMLInputElement).checked);
  });
  byId<HTMLInputElement>("chk-generator").addEventListener("change", (e) => {
    scene.setShowGenerator((e.target as HTMLInputElement).checked);
  });
  byId<HTMLInputElement>("chk-steam").addEventListener("change", (e) => {
    scene.setShowSteam((e.target as HTMLInputElement).checked);
  });
  byId<HTMLInputElement>("chk-orbit-exag").addEventListener("change", (e) => {
    scene.setExaggerateOrbit((e.target as HTMLInputElement).checked);
  });

  const vfIds: Array<[keyof VectorFieldFlags, string]> = [
    ["velocity", "vf-velocity"],
    ["heat", "vf-heat"],
    ["forces", "vf-forces"],
    ["stress", "vf-stress"],
    ["whirl", "vf-whirl"],
    ["acoustic", "vf-acoustic"],
  ];
  const syncVectorFields = () => {
    const flags: VectorFieldFlags = {
      velocity: byId<HTMLInputElement>("vf-velocity").checked,
      heat: byId<HTMLInputElement>("vf-heat").checked,
      forces: byId<HTMLInputElement>("vf-forces").checked,
      stress: byId<HTMLInputElement>("vf-stress").checked,
      whirl: byId<HTMLInputElement>("vf-whirl").checked,
      acoustic: byId<HTMLInputElement>("vf-acoustic").checked,
    };
    scene.setVectorFields(flags);
  };
  const setVectorFlags = (partial: Partial<VectorFieldFlags>) => {
    for (const [key, id] of vfIds) {
      if (partial[key] !== undefined) {
        byId<HTMLInputElement>(id).checked = Boolean(partial[key]);
      }
    }
    syncVectorFields();
  };
  for (const [, id] of vfIds) {
    byId<HTMLInputElement>(id).addEventListener("change", syncVectorFields);
  }

  const selectCamera = (preset: CameraPreset) => {
    scene.flyToPreset(preset);
    document.querySelectorAll<HTMLButtonElement>("[data-cam]").forEach((b) => {
      b.classList.toggle("active", b.dataset.cam === preset);
    });
  };
  document.querySelectorAll<HTMLButtonElement>("[data-cam]").forEach((b) => {
    b.addEventListener("click", () => selectCamera((b.dataset.cam ?? "overview") as CameraPreset));
  });

  // ---------- 9 Guided Experiments & Fast-Forward ----------
  const expList = byId<HTMLDivElement>("experiment-list");
  byId<HTMLButtonElement>("btn-exp-dismiss").addEventListener("click", () => expToast.classList.add("hidden"));

  const renderExperimentCards = () => {
    expList.innerHTML = "";
    for (const exp of EXPERIMENTS) {
      const btn = document.createElement("button");
      btn.className = `exp-card${activeExpId === exp.id ? " active" : ""}`;
      btn.innerHTML = `
        <div class="exp-card-head">
          <strong>${exp.title}</strong>
          <span class="exp-tag">${exp.tag}</span>
        </div>
        <p>${exp.summary}</p>
      `;
      btn.addEventListener("click", () => {
        activeExpId = exp.id;
        fastForwardRemainingS = 0;
        resetControlSlidersToNominal();
        exp.apply(sim);
        if (exp.id === "load_rejection") syncBreakerBtn(false);
        if (exp.id === "air_ingress") ctlAir.set(18);
        if (exp.id === "oil_whip") {
          ctlWhip.set(1.95);
          ctlOilP.set(42);
          ctlUnbal.set(75000);
        }
        if (exp.id === "cold_start") {
          syncBreakerBtn(false);
          ctlLoad.set(28);
          ctlExt.set(0);
          ctlUnbal.set(18000);
        }
        if (exp.id === "wet_steam") {
          ctlPBoiler.set(125);
          ctlTBoiler.set(365);
          ctlExt.set(5);
        }
        if (exp.id === "labyrinth_rub") {
          ctlUnbal.set(280000);
          ctlWhip.set(1.8);
          ctlOilP.set(50);
        }
        if (exp.id === "chp_swing") {
          ctlLoad.set(102);
          ctlExt.set(90);
        }
        if (exp.id === "lube_failure") {
          ctlOilP.set(24);
          ctlCwTemp.set(35);
        }

        campbellTable = sim.campbell_data();
        applyCutaway(exp.cutaway);
        applyShader(exp.shader);
        selectCamera(exp.camera);
        switchRightTab(exp.rightTab);

        expToastTag.textContent = exp.tag;
        expToastTitle.textContent = exp.title;
        expToastDesc.textContent = exp.whatToWatch;
        expToast.classList.remove("hidden");
        renderExperimentCards();
      });
      expList.appendChild(btn);
    }
  };
  renderExperimentCards();

  document.querySelectorAll<HTMLButtonElement>("[data-ff]").forEach((b) => {
    b.addEventListener("click", () => {
      fastForwardRemainingS += Number(b.dataset.ff ?? 2);
    });
  });

  // ---------- Keyboard Shortcuts ----------
  window.addEventListener("keydown", (ev) => {
    if (ev.target instanceof HTMLInputElement || ev.target instanceof HTMLTextAreaElement) return;
    if (ev.code === "Space") {
      ev.preventDefault();
      btnPlay.click();
    } else if (ev.key === "m" || ev.key === "M") {
      btnAudio.click();
    } else if (ev.key === "r" || ev.key === "R") {
      btnReset.click();
    } else if (ev.key === "c" || ev.key === "C") {
      btnQuickCutaway.click();
    } else if (ev.key === "t" || ev.key === "T") {
      btnQuickField.click();
    } else if (ev.key === "f" || ev.key === "F") {
      selectCamera("overview");
    } else if (ev.key === "h" || ev.key === "H") {
      btnZen.click();
    } else if (ev.key === "+" || ev.key === "=") {
      scene.zoomByFactor(0.80);
    } else if (ev.key === "-" || ev.key === "_") {
      scene.zoomByFactor(1.25);
    } else if (ev.key >= "1" && ev.key <= "6") {
      const presets: CameraPreset[] = ["overview", "lp_stage", "hp_valves", "seals", "exhaust", "top_open"];
      selectCamera(presets[Number(ev.key) - 1]!);
    }
  });

  document.addEventListener("visibilitychange", () => {
    if (document.hidden) {
      void audio.suspend();
    } else {
      void audio.resume();
      lastWallMs = performance.now();
    }
  });

  // Expose debug handle on window for automated Playwright / live QA verification
  (window as unknown as Record<string, unknown>).__TURBINE__ = {
    sim,
    scene,
    setCutaway: applyCutaway,
    setShader: applyShader,
    setNavMode: applyNavMode,
    setVectorFlags,
    selectCamera,
    switchRightTab,
    switchLeftTab,
    runDesignOpt,
    applyOptimized: () => lastOptReport && applyDesignToEngine(lastOptReport.optimizedParams),
    getSnapshot: () => decodeSnapshot(sim.snapshot()),
    getStageProfile: () => decodeStageProfile(sim.stage_profile_flat()),
    getFemCert: () => femCert,
  };

  // ---------- Main Animation & Telemetry Loop ----------
  const updateHud = (snap: TurbineSnapshot, stages: StageTelemetry[]) => {
    kpiRpm.innerHTML = `${Math.round(snap.rpm).toLocaleString()}<small>rpm</small>`;
    kpiMw.innerHTML = `${snap.elecMw.toFixed(1)}<small>MW</small>`;
    kpiEff.innerHTML = `${snap.chpEffPct.toFixed(0)} / ${snap.elecEffPct.toFixed(0)}<small>%</small>`;
    kpiCond.innerHTML = `${snap.condKpa.toFixed(1)}<small>kPa</small>`;
    const maxVib = Math.max(snap.frontVibUmPp, snap.rearVibUmPp);
    kpiVib.innerHTML = `${maxVib.toFixed(0)}<small>µm pp</small>`;
    kpiTemp.innerHTML = `${snap.tHpRotorC.toFixed(0)}<small>°C</small>`;

    const st = STATUS_LABELS[snap.statusCode] ?? STATUS_LABELS[0]!;
    statusPill.textContent = st.text;
    statusPill.className = `status-pill status-${st.tone}`;
    mkpiText.textContent = `${Math.round(snap.rpm).toLocaleString()} rpm · ${snap.elecMw.toFixed(1)} MW · ${snap.elecEffPct.toFixed(0)}% η · ${maxVib.toFixed(0)} µm`;
    mkpiStatus.textContent = st.text;
    mkpiStatus.className = `status-pill status-${st.tone}`;

    if (activeRightTab === "thermo") {
      badgeMoisture.textContent = `x = ${snap.lpQualityX.toFixed(3)}`;
      drawMollierDiagram(canvasMollier, mollierDome, sim.mollier_path(), snap);
      drawStageProfileChart(canvasStages, stages, stageChartMode);
      const maxVm = Math.max(...stages.map((s) => s.rootVmMpa), 0);
      const maxMach = Math.max(...stages.map((s) => s.machMax), 0);
      renderMetricGrid(gridThermo, [
        { k: "Live / Chest P", v: `${snap.liveSteamBar.toFixed(0)} / ${snap.chestBar.toFixed(1)} bar` },
        { k: "12-bar Ext P / T", v: `${snap.extHeaderBar.toFixed(1)} bar · ${snap.extHeaderC.toFixed(0)}°C` },
        { k: "Live / Ext / LP Flow", v: `${snap.liveFlowKgS.toFixed(1)} / ${snap.extFlowKgS.toFixed(1)} / ${snap.lpFlowKgS.toFixed(1)} kg/s` },
        { k: "HP / IP-LP η_is", v: `${snap.hpEtaPct.toFixed(1)}% / ${snap.iplpEtaPct.toFixed(1)}%` },
        { k: "Peak Blade σ_vM / Ma", v: `${maxVm.toFixed(0)} MPa · Ma ${maxMach.toFixed(2)}`, tone: maxVm > 450 ? "warn" : "ok" },
        {
          k: "Condenser U / h_cond",
          v: `${snap.uCondWm2K.toFixed(0)} / ${snap.hCondWm2K.toFixed(0)} W/m²K`,
          tone: snap.uCondWm2K < 1400 ? "warn" : "ok",
        },
        { k: "Condenser Duty", v: `${snap.condDutyMw.toFixed(1)} MW (${snap.cwOutC.toFixed(1)}°C)` },
        {
          k: "CV / ESV / Ext Valve",
          v: `${(snap.cvLift * 100).toFixed(0)}% / ${(snap.esvOpen * 100).toFixed(0)}% / ${(snap.extLift * 100).toFixed(0)}%`,
        },
        { k: "Gen δ / Current", v: `${snap.loadAngleDeg.toFixed(1)}° · ${(snap.statorCurrentA / 1000).toFixed(2)} kA` },
        { k: "Axial Thrust / Windage", v: `${snap.axialThrustKn.toFixed(1)} kN · ${snap.windageKw.toFixed(0)} kW` },
      ]);
    } else if (activeRightTab === "rotor") {
      badgeAbscissa.textContent = `Re(λ) = ${snap.spectralAbscissaPerS.toFixed(1)}/s`;
      badgeEproc.textContent = `log E = ${snap.isoLogE.toFixed(2)}`;
      drawCampbellDiagram(canvasCampbell, campbellTable, snap);
      drawOrbitPlot(canvasOrbit, sim.orbit_history(), snap, activeDesign.bearingClearanceUm);
      renderMetricGrid(gridRotor, [
        {
          k: "Stability max Re(λ)",
          v: `${snap.spectralAbscissaPerS.toFixed(2)} 1/s`,
          tone: snap.spectralAbscissaPerS > 0 ? "crit" : snap.spectralAbscissaPerS > -5 ? "warn" : "ok",
        },
        { k: "Whirl Freq / 1st Crit", v: `${snap.dominantWhirlHz.toFixed(1)} Hz · ${snap.crit1Rpm.toFixed(0)} rpm` },
        {
          k: "Front / Rear Vib",
          v: `${snap.frontVibUmPp.toFixed(1)} / ${snap.rearVibUmPp.toFixed(1)} µm`,
          tone: maxVib > 140 ? "crit" : maxVib > 80 ? "warn" : "ok",
        },
        {
          k: "Min Oil Film h_min",
          v: `${snap.frontFilmUm.toFixed(1)} / ${snap.rearFilmUm.toFixed(1)} µm`,
          tone: snap.frontFilmUm < 20 ? "crit" : snap.frontFilmUm < 35 ? "warn" : "ok",
        },
        { k: "Alford K_xy", v: `${snap.alfordKxyMnM.toFixed(2)} MN/m` },
        {
          k: "Labyrinth Rub Force",
          v: `${snap.rubForceKn.toFixed(2)} kN`,
          tone: snap.rubForceKn > 0.1 ? "crit" : "ok",
        },
        { k: "Oil Visc / Ecc ε", v: `${snap.oilViscMpaS.toFixed(1)} mPa·s · ${snap.frontEccentricity.toFixed(2)}` },
        { k: "Bearing + Rub Heat", v: `${snap.bearingHeatKw.toFixed(1)} kW` },
      ]);
    } else if (activeRightTab === "audio") {
      badgeSpl.textContent = `${snap.splDb.toFixed(1)} dB`;
      audio.drawScopeAndFft(canvasScope, canvasFft, snap);
      renderMetricGrid(gridAudio, [
        { k: "SPL (20 µPa ref)", v: `${snap.splDb.toFixed(1)} dB` },
        { k: "ISO 532-1 Loudness", v: `${snap.loudnessSone.toFixed(1)} sone` },
        { k: "DIN 45692 Sharpness", v: `${snap.sharpnessAcum.toFixed(2)} acum` },
        { k: "LP Last-Stage BPF", v: `${snap.bpfLpHz.toFixed(0)} Hz (44×)` },
      ]);
    } else if (activeRightTab === "ledger") {
      const relE = Math.abs(snap.ledgerEnergyResJ) / (Math.abs(snap.ledgerBoundaryMj) * 1e6 + 1);
      renderMetricGrid(gridLedger, [
        {
          k: "First-Law Residual",
          v: `${relE.toExponential(2)} rel`,
          tone: relE < 1e-5 ? "ok" : "warn",
        },
        {
          k: "Stream Mass Residual",
          v: `${snap.ledgerMassResKg.toExponential(2)} kg`,
          tone: Math.abs(snap.ledgerMassResKg) < 1e-5 ? "ok" : "warn",
        },
        {
          k: "2nd-Law Entropy Gen",
          v: `+${(snap.ledgerEntropyGenJK * 1e-3).toFixed(1)} kJ/K`,
          tone: snap.ledgerEntropyGenJK >= 0 ? "ok" : "crit",
        },
        { k: "Sim Time / Steps", v: `${snap.timeS.toFixed(2)} s · ${snap.steps.toFixed(0)} st` },
        { k: "Cumulative Energy", v: `${snap.ledgerBoundaryMj.toFixed(2)} MJ` },
        { k: "Last Newton Iters", v: `${snap.lastNewtonIters.toFixed(0)} iters` },
      ]);
      valDigest.textContent = sim.state_digest_hex();

      renderMetricGrid(gridThermal, [
        { k: "HP Rotor / Casing", v: `${snap.tHpRotorC.toFixed(1)} / ${snap.tHpCasingC.toFixed(1)} °C` },
        { k: "IP Rotor / Casing", v: `${snap.tIpRotorC.toFixed(1)} / ${snap.tIpCasingC.toFixed(1)} °C` },
        { k: "LP Hood / Oil Sump", v: `${snap.tLpHoodC.toFixed(1)} / ${snap.tLubeOilC.toFixed(1)} °C` },
        {
          k: "Front / Rear Babbitt",
          v: `${snap.tFrontBrgC.toFixed(1)} / ${snap.tRearBrgC.toFixed(1)} °C`,
          tone: Math.max(snap.tFrontBrgC, snap.tRearBrgC) > 105 ? "crit" : "ok",
        },
        { k: "Diff Axial Expansion", v: `${snap.diffExpansionMm.toFixed(2)} mm` },
        { k: "HP Thermal Stress", v: `${snap.hpThermalStressMpa.toFixed(1)} MPa` },
        { k: "Norton Creep Rate", v: `${snap.hpCreepRatePerS.toExponential(1)} /s` },
        { k: "Thermal Derating", v: `${(snap.deratingFrac * 100).toFixed(0)}%` },
      ]);
    }
  };

  bootOverlay.classList.add("fade-out");
  setTimeout(() => bootOverlay.remove(), 340);

  const tick = (nowMs: number) => {
    const dtWallS = Math.min(0.05, Math.max(0.001, (nowMs - lastWallMs) * 0.001));
    lastWallMs = nowMs;

    if (fastForwardRemainingS > 0) {
      const tStart = performance.now();
      while (fastForwardRemainingS > 0 && performance.now() - tStart < 10.0) {
        const chunk = Math.min(0.05, fastForwardRemainingS);
        if (!sim.step(chunk)) {
          fastForwardRemainingS = 0;
          break;
        }
        fastForwardRemainingS -= chunk;
      }
    } else if (running) {
      sim.step(dtWallS * simRate);
    }

    const acousticOn = byId<HTMLInputElement>("vf-acoustic").checked;
    const micActive = activeRightTab === "audio" || acousticOn || btnAudio.classList.contains("active-audio");
    scene.setMicVisible(micActive);
    audio.pump(sim, running || fastForwardRemainingS > 0, activeRightTab === "audio" || acousticOn);
    const snap = decodeSnapshot(sim.snapshot());
    const stages = decodeStageProfile(sim.stage_profile_flat());
    const continuumFlat = sim.continuum_field_flat();
    const orbitFlat = sim.orbit_history();
    scene.updateAndRender(snap, stages, continuumFlat, orbitFlat, dtWallS, running, activeDesign.bearingClearanceUm);

    if (frameCounter % 3 === 0) {
      updateHud(snap, stages);
    }
    frameCounter++;
    requestAnimationFrame(tick);
  };

  hudReady = true;
  updateHud(decodeSnapshot(sim.snapshot()), decodeStageProfile(sim.stage_profile_flat()));
  requestAnimationFrame(tick);
}

boot().catch((err) => {
  const bootStatus = document.getElementById("boot-status");
  if (bootStatus) {
    bootStatus.textContent = `Boot failed: ${err instanceof Error ? err.message : String(err)}`;
  }
  console.error(err);
});
