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
  type QualityTier,
  type VectorFieldLayer,
  type ViewInsets,
} from "./watchScene";

/** Must match the html.is-compact rules in styles.css (short landscape screens stay "wide"). */
const COMPACT_QUERY = "(max-width: 999px) and (min-height: 561px)";
const SHORT_QUERY = "(max-height: 560px)";
const STORE_COACH = "fs-pocket-watch.coach.v1";
const STORE_QUALITY = "fs-pocket-watch.quality";
const STORE_INSPECTOR = "fs-pocket-watch.inspector";

type Vec3 = [number, number, number];
type PanelName = "bench" | "experiments" | "movement" | "metrology";
type SheetState = "peek" | "half" | "full";

/**
 * Engine methods that newer builds add. Each one is feature-detected so the front end runs on any
 * engine build; absent methods fall back to what the snapshot already carries.
 */
interface OptionalEngineApi {
  escapement_state?: () => Float64Array;
  timegrapher_drain?: () => Float64Array;
  grid_n?: () => number;
  grid_half_span_mm?: () => number;
  ambient_temperature_k?: () => number;
}

function $<T extends Element = HTMLElement>(sel: string): T {
  const el = document.querySelector<T>(sel);
  if (!el) throw new Error(`Missing DOM element: ${sel}`);
  return el;
}

function $$<T extends Element = HTMLElement>(sel: string): T[] {
  return Array.from(document.querySelectorAll<T>(sel));
}

function storeGet(key: string): string | null {
  try {
    return window.localStorage.getItem(key);
  } catch {
    return null;
  }
}

function storeSet(key: string, value: string): void {
  try {
    window.localStorage.setItem(key, value);
  } catch {
    /* private mode or blocked storage: the preference just isn't remembered */
  }
}

function setText(el: Element, text: string): void {
  if (el.textContent !== text) el.textContent = text;
}

function setAttr(el: Element, name: string, value: string): void {
  if (el.getAttribute(name) !== value) el.setAttribute(name, value);
}

function signed(v: number, digits: number): string {
  const s = Math.abs(v).toFixed(digits);
  return `${v < 0 && Number(s) !== 0 ? "−" : "+"}${s}`;
}

/** Keep a range input's filled-track percentage in sync with its value. */
function syncRangeFill(input: HTMLInputElement): void {
  const min = Number(input.min || 0);
  const max = Number(input.max || 1);
  const pct = ((Number(input.value) - min) / Math.max(1e-9, max - min)) * 100;
  input.style.setProperty("--fill", `${Math.max(0, Math.min(100, pct)).toFixed(1)}%`);
}

function setRange(input: HTMLInputElement, value: number): void {
  input.value = String(value);
  syncRangeFill(input);
}

function detectQualityTier(): QualityTier {
  const q = new URLSearchParams(window.location.search).get("quality");
  if (q === "high" || q === "medium" || q === "low") return q;
  const coarse = window.matchMedia("(pointer: coarse)").matches;
  const minSide = Math.min(window.screen.width, window.screen.height);
  const cores = navigator.hardwareConcurrency || 4;
  const mem = (navigator as unknown as { deviceMemory?: number }).deviceMemory ?? 8;
  if (coarse && minSide < 600) return "low"; // phones
  if (coarse || cores <= 4 || mem <= 4) return "medium"; // tablets, small laptops
  return "high";
}

async function boot(): Promise<void> {
  const splash = $("#boot-splash");
  const bootStatus = $("#boot-status");
  try {
    await initWasm();
  } catch (err) {
    setText(bootStatus, "This browser could not start the WebAssembly physics engine.");
    console.error(err);
    return;
  }
  setText(bootStatus, "Lighting the bench…");

  const root = document.documentElement;
  const shell = $("#app-shell");
  const sim = new PocketWatchSim();
  const ext = sim as unknown as OptionalEngineApi;
  const hasEscState = typeof ext.escapement_state === "function";
  const hasDrain = typeof ext.timegrapher_drain === "function";
  const micAudio = new VirtualMicrophoneAudio();
  const timegrapher = new HorologicalTimegrapher();

  const reducedMotionMq = window.matchMedia("(prefers-reduced-motion: reduce)");
  const compactMq = window.matchMedia(COMPACT_QUERY);
  let isCompact = compactMq.matches;
  root.classList.toggle("is-compact", isCompact);
  // Short landscape (phones on their side) keeps the wide layout but folds the view dock into a menu.
  const shortMq = window.matchMedia(SHORT_QUERY);
  let isShort = shortMq.matches;
  root.classList.toggle("is-short", isShort);

  const detectedTier = detectQualityTier();
  const storedQuality = storeGet(STORE_QUALITY);
  const qualityChoice = storedQuality === "high" || storedQuality === "medium" || storedQuality === "low" ? storedQuality : "auto";
  let tier: QualityTier = qualityChoice === "auto" ? detectedTier : qualityChoice;
  const applyTierClass = (): void => {
    root.classList.remove("tier-high", "tier-medium", "tier-low");
    root.classList.add(`tier-${tier}`);
  };
  applyTierClass();

  // ---------------------------------------------------------------------------------------------
  // Tool & environment state
  let activeToolMode: ActiveToolMode = "inspect";
  let torchPower = 0.0;
  let torchPosMm: Vec3 = [0.0, -24.0, 2.0];
  let ln2Rate = 0.0;
  let ln2PosMm: Vec3 = [0.0, -16.0, 3.0];
  let oilRate = 0.0;
  let oilPosMm: Vec3 = [0.0, -11.0, 2.0];
  let magnetOn = false;
  let magnetDistMm = 28.0;
  let magnetPosMm: Vec3 = [28.0, -11.0, 2.5];
  let airPower = 0.0;
  let airPosMm: Vec3 = [0.0, -20.0, 3.0];
  let ambientPressureAtm = 1.0;
  let gravityGVec: Vec3 = [0.0, -1.0, 0.0];
  let gravityLabel = "1.0 g · 6h down";
  let micPosMm: Vec3 = [0.0, -11.0, 11.0];
  let micGain = 1.4;
  let slowMo = false;
  let simPaused = false;

  // ---------------------------------------------------------------------------------------------
  // DOM
  const canvasViewport = $<HTMLCanvasElement>("#viewport-canvas");
  const canvasScope = $<HTMLCanvasElement>("#canvas-oscilloscope");
  const canvasTimegrapher = $<HTMLCanvasElement>("#canvas-timegrapher");
  const tbodyLedger = $<HTMLElement>("#node-ledger-body");
  const topBar = $("#top-bar");
  const toolRail = $("#tool-rail");
  const inspector = $("#inspector");
  const inspectorBody = $("#inspector-body");
  const viewDock = $("#view-dock");
  const btnViewMenu = $<HTMLButtonElement>("#btn-view-menu");

  const hudWatchTime = $("#hud-watch-time");
  const hudLabTime = $("#hud-lab-time");
  const hudTimeDrift = $("#hud-time-drift");
  const hudRateSDay = $("#hud-rate-s-day");
  const hudAmp = $("#hud-amp");
  const hudBeatError = $("#hud-beat-error");
  const hudStatusPill = $("#hud-status-pill");
  const statusHead = hudStatusPill.querySelector(".status-head")!;
  const statusDetail = hudStatusPill.querySelector(".status-detail")!;
  const statusShort = hudStatusPill.querySelector(".status-short")!;
  const hudDigest = $("#hud-digest");
  const hudPerf = $("#hud-perf");
  const btnPauseSim = $<HTMLButtonElement>("#btn-pause-sim");
  const btnTheaterMode = $<HTMLButtonElement>("#btn-theater-mode");
  const btnInspectorToggle = $<HTMLButtonElement>("#btn-inspector-toggle");
  const btnStopPerturbations = $<HTMLButtonElement>("#btn-stop-perturbations");

  const stageToolBanner = $("#stage-tool-banner");
  const stageBannerTitle = $("#stage-banner-title");
  const stageBannerText = $("#stage-banner-text");
  const bannerIconUse = $<SVGUseElement>("#banner-icon-use");
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
  const selQuality = $<HTMLSelectElement>("#sel-quality");

  const sliderTorchPower = $<HTMLInputElement>("#slider-torch-power");
  const lblTorchPower = $("#lbl-torch-power");
  const sliderLn2Rate = $<HTMLInputElement>("#slider-ln2-rate");
  const lblLn2Rate = $("#lbl-ln2-rate");
  const sliderOilRate = $<HTMLInputElement>("#slider-oil-rate");
  const lblOilRate = $("#lbl-oil-rate");
  const chkMagnetOn = $<HTMLInputElement>("#chk-magnet-on");
  const sliderMagnetDist = $<HTMLInputElement>("#slider-magnet-dist");
  const lblMagnetDist = $("#lbl-magnet-dist");
  const lblBField = $("#lbl-b-field");
  const lblRemanence = $("#lbl-remanence");
  const sliderAirPower = $<HTMLInputElement>("#slider-air-power");
  const lblAirPower = $("#lbl-air-power");
  const sliderPressureAtm = $<HTMLInputElement>("#slider-pressure-atm");
  const lblPressure = $("#lbl-pressure");
  const lblChamberEnv = $("#lbl-chamber-env");
  const btnShock500g = $<HTMLButtonElement>("#btn-shock-500g");

  const sliderMicGain = $<HTMLInputElement>("#slider-mic-gain");
  const lblMicGain = $("#lbl-mic-gain");
  const lblMicSpl = $("#lbl-mic-spl");

  const tgRate = $("#tg-rate");
  const tgAmp = $("#tg-amp");
  const tgBeatErr = $("#tg-beat-err");
  const tgSource = $("#tg-source");

  const ovCaseTemp = $("#ov-case-temp");
  const ovGoldMelt = $("#ov-gold-melt");
  const ovFouledGold = $("#ov-fouled-gold");
  const ovBField = $("#ov-b-field");

  const escSource = $("#esc-source");
  const escNote = $("#esc-note");
  const escFields = {
    impulse: $("#esc-impulse"),
    eff: $("#esc-eff"),
    lift: $("#esc-lift"),
    q: $("#esc-q"),
    torque: $("#esc-torque"),
    friction: $("#esc-friction"),
    air: $("#esc-air"),
    reserve: $("#esc-reserve"),
  };
  const contactEls = $$<HTMLElement>(".contact[data-bit]");
  const chrSens = $("#chr-sens");
  const chrInterval = $("#chr-interval");
  const chrCosc = $("#chr-cosc");
  const chrDrag = $("#chr-drag");

  const ffOverlay = $("#ff-overlay");
  const ffTitle = $("#ff-title");
  const ffDetail = $("#ff-detail");
  const ffBarFill = $("#ff-bar-fill");
  const coach = $("#coach");
  const helpDialog = $<HTMLDialogElement>("#help-dialog");
  const announcer = $("#sr-announcer");

  const boundEls = new Map<string, HTMLElement[]>();
  for (const el of $$<HTMLElement>("[data-bind]")) {
    const key = el.dataset.bind ?? "";
    const list = boundEls.get(key) ?? [];
    list.push(el);
    boundEls.set(key, list);
  }
  const bindText = (key: string, text: string): void => {
    for (const el of boundEls.get(key) ?? []) setText(el, text);
  };
  const announce = (msg: string): void => {
    announcer.textContent = "";
    window.setTimeout(() => (announcer.textContent = msg), 30);
  };

  const allRanges = $$<HTMLInputElement>('input[type="range"]');
  for (const r of allRanges) {
    syncRangeFill(r);
    r.addEventListener("input", () => syncRangeFill(r));
  }

  // ---------------------------------------------------------------------------------------------
  // Scene
  const scene = new PocketWatchScene(
    canvasViewport,
    {
      onTargetMoved(tool, posMm) {
        dismissCoach();
        if (tool === "torch") {
          torchPosMm = posMm;
          if (torchPower < 0.05) {
            torchPower = 0.78;
            setRange(sliderTorchPower, torchPower);
          }
        } else if (tool === "ln2") {
          ln2PosMm = posMm;
          if (ln2Rate < 0.05) {
            ln2Rate = 0.8;
            setRange(sliderLn2Rate, ln2Rate);
          }
        } else if (tool === "oil") {
          oilPosMm = posMm;
          if (oilRate < 0.05) {
            oilRate = 0.65;
            setRange(sliderOilRate, oilRate);
          }
        } else if (tool === "magnet") {
          magnetOn = true;
          chkMagnetOn.checked = true;
          magnetPosMm = posMm;
          magnetDistMm = Math.round(Math.max(8, Math.min(80, Math.hypot(posMm[0] + 2.2, posMm[1] + 11.0))));
          setRange(sliderMagnetDist, magnetDistMm);
        } else if (tool === "air") {
          airPosMm = posMm;
          if (airPower < 0.05) {
            airPower = 0.85;
            setRange(sliderAirPower, airPower);
          }
        } else if (tool === "mic") {
          micPosMm = posMm;
        }
        syncToolsToSim();
      },
    },
    { tier, reducedMotion: reducedMotionMq.matches },
  );

  if (typeof ext.grid_half_span_mm === "function") {
    const h = ext.grid_half_span_mm();
    if (Number.isFinite(h)) scene.setGridHalfSpanMm(h);
  }

  // ---------------------------------------------------------------------------------------------
  // Perturbation state -> engine, banner and badges
  const anyPerturbationActive = (): boolean =>
    torchPower > 0.01 ||
    ln2Rate > 0.01 ||
    oilRate > 0.01 ||
    magnetOn ||
    airPower > 0.01 ||
    Math.abs(ambientPressureAtm - 1.0) > 0.05 ||
    Math.abs(gravityGVec[0]) + Math.abs(gravityGVec[1] + 1) + Math.abs(gravityGVec[2]) > 0.01;

  const TOOL_META: Record<Exclude<ActiveToolMode, "inspect">, { icon: string; tone: string; name: string; wide: string; compact: string }> = {
    torch: { icon: "#i-flame", tone: "torch", name: "Blowtorch", wide: "Drag on the watch to aim the 2250 K flame · right-drag to orbit", compact: "Drag to aim the flame · two fingers to orbit" },
    ln2: { icon: "#i-snow", tone: "cryo", name: "Liquid nitrogen", wide: "Drag on the watch to pour 77 K nitrogen · right-drag to orbit", compact: "Drag to pour · two fingers to orbit" },
    oil: { icon: "#i-drop", tone: "oil", name: "Watch oil", wide: "Drag on the watch to place the oiler · right-drag to orbit", compact: "Drag to place the oiler · two fingers to orbit" },
    magnet: { icon: "#i-magnet", tone: "magnet", name: "NdFeB magnet", wide: "Drag on the watch to sweep the magnet · right-drag to orbit", compact: "Drag to move the magnet · two fingers to orbit" },
    air: { icon: "#i-wind", tone: "air", name: "Air jet", wide: "Drag on the watch to aim the jet · right-drag to orbit", compact: "Drag to aim the jet · two fingers to orbit" },
    mic: { icon: "#i-mic", tone: "mic", name: "Microphone", wide: "Drag on the watch to place the microphone · right-drag to orbit", compact: "Drag to place the mic · two fingers to orbit" },
  };

  const updatePerturbationBannerAndBadges = (): void => {
    const anyActive = anyPerturbationActive();
    btnStopPerturbations.classList.toggle("hidden", !anyActive);

    const focus: Record<string, boolean> = {
      torch: activeToolMode === "torch" || torchPower > 0.01,
      ln2: activeToolMode === "ln2" || ln2Rate > 0.01,
      oil: activeToolMode === "oil" || oilRate > 0.01,
      magnet: activeToolMode === "magnet" || magnetOn,
      air: activeToolMode === "air" || airPower > 0.01,
      chamber: Math.abs(ambientPressureAtm - 1.0) > 0.05 || Math.abs(gravityGVec[1] + 1) + Math.abs(gravityGVec[0]) + Math.abs(gravityGVec[2]) > 0.01,
    };
    for (const [k, on] of Object.entries(focus)) {
      $(`#subtool-${k}`).classList.toggle("is-focused", on);
    }
    const live: Record<string, boolean> = {
      torch: torchPower > 0.01,
      ln2: ln2Rate > 0.01,
      oil: oilRate > 0.01,
      magnet: magnetOn,
      air: airPower > 0.01,
    };
    for (const b of $$<HTMLButtonElement>(".tool-btn[data-tool]")) {
      b.classList.toggle("is-live", !!live[b.dataset.tool ?? ""] && b.dataset.tool !== activeToolMode);
    }

    canvasViewport.classList.toggle("tool-cursor", activeToolMode !== "inspect");

    if (activeToolMode !== "inspect") {
      const meta = TOOL_META[activeToolMode];
      const pct =
        activeToolMode === "torch" ? torchPower : activeToolMode === "ln2" ? ln2Rate : activeToolMode === "oil" ? oilRate : activeToolMode === "air" ? airPower : -1;
      setText(stageBannerTitle, pct >= 0 ? `${meta.name} · ${Math.round(pct * 100)}%` : meta.name);
      setText(stageBannerText, isCompact ? meta.compact : meta.wide);
      bannerIconUse.setAttribute("href", meta.icon);
      stageToolBanner.dataset.tone = meta.tone;
      btnBannerOrbit.classList.remove("hidden");
      btnBannerStop.classList.toggle("hidden", !anyActive);
      stageToolBanner.classList.remove("hidden");
    } else if (anyActive) {
      const names: string[] = [];
      if (torchPower > 0.01) names.push(`torch ${Math.round(torchPower * 100)}%`);
      if (ln2Rate > 0.01) names.push(`LN₂ ${Math.round(ln2Rate * 100)}%`);
      if (oilRate > 0.01) names.push(`oil ${Math.round(oilRate * 100)}%`);
      if (magnetOn) names.push(`magnet at ${magnetDistMm.toFixed(0)} mm`);
      if (airPower > 0.01) names.push(`air ${Math.round(airPower * 100)}%`);
      if (Math.abs(ambientPressureAtm - 1.0) > 0.05) names.push(`${ambientPressureAtm.toFixed(2)} atm`);
      if (focus.chamber && Math.abs(ambientPressureAtm - 1.0) <= 0.05) names.push(gravityLabel);
      setText(stageBannerTitle, "Perturbations running");
      setText(stageBannerText, names.join(" · "));
      bannerIconUse.setAttribute("href", "#i-sparkle");
      stageToolBanner.dataset.tone = "gold";
      btnBannerOrbit.classList.add("hidden");
      btnBannerStop.classList.remove("hidden");
      stageToolBanner.classList.remove("hidden");
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
    setText(lblMagnetDist, `${magnetDistMm.toFixed(0)} mm`);
    setText(lblPressure, `${ambientPressureAtm.toFixed(2)} atm`);
    setText(lblMicGain, `${micGain.toFixed(2)}×`);
    updatePerturbationBannerAndBadges();
  };
  syncToolsToSim();
  // The page's switch, not the engine's constructor default, decides the hairspring alloy.
  sim.set_compensated_hairspring(chkInvarSpring.checked);

  // ---------------------------------------------------------------------------------------------
  // Layout: inspector panel / bottom sheet, focus mode, view insets for camera framing
  let layoutDirty = true;
  let reframeOnLayout = true;
  // Chrome changed size (sheet snap, inspector toggle, rotation): refit the camera preset, unless
  // the user has orbited or zoomed since, in which case only the projection window moves.
  let refitPending = false;
  let sheetState: SheetState = "peek";
  let inspectorOpen = (storeGet(STORE_INSPECTOR) ?? (window.innerWidth >= 1100 ? "open" : "closed")) === "open";
  let theater = false;
  const insetsCur: ViewInsets = { left: 0, right: 0, top: 0, bottom: 0 };
  const insetsTarget: ViewInsets = { left: 0, right: 0, top: 0, bottom: 0 };

  const setInspectorOpen = (open: boolean, persist = true): void => {
    inspectorOpen = open;
    shell.dataset.inspector = open ? "open" : "closed";
    setAttr(btnInspectorToggle, "aria-expanded", String(open));
    setAttr(btnInspectorToggle, "aria-pressed", String(open));
    if (persist) storeSet(STORE_INSPECTOR, open ? "open" : "closed");
    layoutDirty = true;
    refitPending = true;
  };
  setInspectorOpen(inspectorOpen, false);

  const setSheet = (state: SheetState): void => {
    sheetState = state;
    shell.dataset.sheet = state;
    layoutDirty = true;
    refitPending = true;
  };

  const sheetVisibleHeight = (state: SheetState): number => {
    const h = inspector.offsetHeight;
    if (state === "full") return h;
    if (state === "half") return Math.min(h, window.innerHeight * 0.5);
    return peekHeight;
  };
  // Peek = drag handle + tab row, measured so CSS and the camera framing agree.
  let peekHeight = 60;
  const measurePeek = (): void => {
    const handle = $("#sheet-handle");
    const head = inspector.querySelector<HTMLElement>(".inspector-head")!;
    const h = Math.round(handle.offsetHeight + head.offsetHeight);
    if (h > 20 && h !== peekHeight) {
      peekHeight = h;
      root.style.setProperty("--peek-h", `${h}px`);
    }
  };

  const computeLayout = (): void => {
    const W = shell.clientWidth;
    const H = shell.clientHeight;
    const tb = topBar.getBoundingClientRect();
    root.style.setProperty("--topbar-bottom", `${Math.round(tb.bottom)}px`);
    const top = tb.bottom + 6;
    if (isCompact) {
      measurePeek();
      const toolbarH = toolRail.offsetHeight;
      const visible = sheetVisibleHeight(sheetState === "full" ? "half" : sheetState);
      insetsTarget.top = top;
      insetsTarget.bottom = toolbarH + visible + 8;
      insetsTarget.left = 0;
      insetsTarget.right = 0;
      root.style.removeProperty("--stage-cx");
    } else {
      const rail = toolRail.getBoundingClientRect();
      insetsTarget.top = theater ? Math.min(top, 24) : top;
      insetsTarget.left = theater ? 0 : rail.right + 6;
      const edge = Math.max(0, W - inspector.getBoundingClientRect().right);
      insetsTarget.right = inspectorOpen && !theater ? inspector.offsetWidth + edge + 8 : 0;
      const dockH = viewDock.offsetHeight;
      root.style.setProperty("--dock-h", `${dockH}px`);
      // On short screens the dock is a popover menu, so it takes no room at the bottom.
      insetsTarget.bottom = theater || isShort ? 0 : dockH + 20;
      const cx = insetsTarget.left + (W - insetsTarget.left - insetsTarget.right) / 2;
      root.style.setProperty("--stage-cx", `${Math.round(cx)}px`);
      root.style.setProperty("--stage-r", `${Math.round(Math.max(insetsTarget.right, edge + 8))}px`);
    }
    // While the first-run coach is up, frame the whole watch beside or above it, not behind it.
    if (!coach.classList.contains("hidden")) {
      const cr = coach.getBoundingClientRect();
      if (cr.height > 0) {
        if (cr.left > W * 0.45) insetsTarget.right = Math.max(insetsTarget.right, W - cr.left + 10);
        else if (cr.top > 0) insetsTarget.bottom = Math.max(insetsTarget.bottom, H - cr.top + 10);
      }
    }
    // Never let the chrome squeeze the stage below a usable window.
    const minH = H * 0.32;
    if (H - insetsTarget.top - insetsTarget.bottom < minH) {
      insetsTarget.bottom = Math.max(0, H - insetsTarget.top - minH);
    }
  };

  const stepInsets = (dt: number): void => {
    if (layoutDirty) {
      computeLayout();
      layoutDirty = false;
      if (reframeOnLayout) {
        Object.assign(insetsCur, insetsTarget);
        scene.setViewInsets({ ...insetsCur });
        scene.reframe(false);
        reframeOnLayout = false;
        refitPending = false;
        return;
      }
      if (refitPending) {
        refitPending = false;
        if (scene.isAtPreset()) scene.reframe(true, { ...insetsTarget });
      }
    }
    const k = reducedMotionMq.matches ? 1 : 1 - Math.exp(-dt * 10);
    let moved = false;
    for (const key of ["left", "right", "top", "bottom"] as const) {
      const d = insetsTarget[key] - insetsCur[key];
      if (Math.abs(d) > 0.25) {
        insetsCur[key] += Math.abs(d) < 1 ? d : d * k;
        moved = true;
      }
    }
    if (moved) scene.setViewInsets({ ...insetsCur });
  };

  if (typeof ResizeObserver !== "undefined") {
    const ro = new ResizeObserver(() => (layoutDirty = true));
    ro.observe(shell);
    ro.observe(topBar);
    ro.observe(viewDock);
  }
  window.addEventListener("resize", () => {
    layoutDirty = true;
    refitPending = true;
  });
  inspector.addEventListener("transitionend", () => (layoutDirty = true));
  compactMq.addEventListener("change", () => {
    isCompact = compactMq.matches;
    root.classList.toggle("is-compact", isCompact);
    shell.classList.remove("view-open");
    setAttr(btnViewMenu, "aria-expanded", "false");
    layoutDirty = true;
    reframeOnLayout = true;
    updatePerturbationBannerAndBadges();
  });
  shortMq.addEventListener("change", () => {
    isShort = shortMq.matches;
    root.classList.toggle("is-short", isShort);
    shell.classList.remove("view-open");
    setAttr(btnViewMenu, "aria-expanded", "false");
    layoutDirty = true;
    reframeOnLayout = true;
  });

  const toggleTheaterMode = (): void => {
    theater = !shell.classList.contains("theater-mode");
    shell.classList.toggle("theater-mode", theater);
    setAttr(btnTheaterMode, "aria-pressed", String(theater));
    layoutDirty = true;
    refitPending = true;
  };
  btnTheaterMode.addEventListener("click", toggleTheaterMode);
  btnInspectorToggle.addEventListener("click", () => {
    if (theater) toggleTheaterMode();
    setInspectorOpen(!inspectorOpen);
  });
  $("#btn-inspector-close").addEventListener("click", () => setInspectorOpen(false));

  // Tabs (roving tabindex, arrow keys)
  const tabs = $$<HTMLButtonElement>('[role="tab"][data-panel]');
  let currentPanel: PanelName = "bench";
  const panelVisible = (name: PanelName): boolean =>
    currentPanel === name && !theater && (isCompact ? sheetState !== "peek" : inspectorOpen);

  const setPanel = (name: PanelName, reveal = true): void => {
    currentPanel = name;
    for (const t of tabs) {
      const on = t.dataset.panel === name;
      setAttr(t, "aria-selected", String(on));
      t.tabIndex = on ? 0 : -1;
      const panel = document.getElementById(t.getAttribute("aria-controls") ?? "");
      if (panel) panel.hidden = !on;
    }
    if (reveal) {
      if (isCompact) {
        if (sheetState === "peek") setSheet("half");
      } else if (!inspectorOpen) {
        setInspectorOpen(true);
      }
    }
    inspectorBody.scrollTop = 0;
  };

  for (const t of tabs) {
    t.addEventListener("click", () => {
      const name = (t.dataset.panel ?? "bench") as PanelName;
      if (isCompact && name === currentPanel && sheetState !== "peek") {
        setSheet("peek");
        return;
      }
      setPanel(name, true);
    });
    t.addEventListener("keydown", (ev) => {
      const i = tabs.indexOf(t);
      let j = -1;
      if (ev.key === "ArrowRight" || ev.key === "ArrowDown") j = (i + 1) % tabs.length;
      else if (ev.key === "ArrowLeft" || ev.key === "ArrowUp") j = (i - 1 + tabs.length) % tabs.length;
      else if (ev.key === "Home") j = 0;
      else if (ev.key === "End") j = tabs.length - 1;
      if (j >= 0) {
        ev.preventDefault();
        const next = tabs[j]!;
        setPanel((next.dataset.panel ?? "bench") as PanelName, true);
        next.focus();
      }
    });
  }

  const revealSubtool = (id: string): void => {
    const box = document.getElementById(id);
    if (!box || !panelVisible("bench")) return;
    const top =
      box.getBoundingClientRect().top - inspectorBody.getBoundingClientRect().top + inspectorBody.scrollTop - 12;
    inspectorBody.scrollTo({ top, behavior: reducedMotionMq.matches ? "auto" : "smooth" });
  };

  // Bottom sheet drag (handle + tab header)
  {
    let startY = 0;
    let startOffset = 0;
    let dragging = false;
    let pid = -1;
    let suppressClick = false;
    const samples: { t: number; y: number }[] = [];
    const currentOffset = (): number => {
      const m = new DOMMatrixReadOnly(getComputedStyle(inspector).transform);
      return m.m42;
    };
    const offsets = (): Record<SheetState, number> => {
      const h = inspector.offsetHeight;
      return { full: 0, half: h - sheetVisibleHeight("half"), peek: h - sheetVisibleHeight("peek") };
    };
    const onDown = (ev: PointerEvent): void => {
      if (!isCompact || (ev.pointerType === "mouse" && ev.button !== 0)) return;
      startY = ev.clientY;
      startOffset = currentOffset();
      pid = ev.pointerId;
      dragging = false;
      samples.length = 0;
      samples.push({ t: performance.now(), y: ev.clientY });
    };
    const onMove = (ev: PointerEvent): void => {
      if (ev.pointerId !== pid) return;
      const dy = ev.clientY - startY;
      if (!dragging && Math.abs(dy) > 7) {
        dragging = true;
        inspector.classList.add("is-dragging");
        (ev.currentTarget as HTMLElement).setPointerCapture(pid);
      }
      if (!dragging) return;
      const o = offsets();
      const y = Math.max(o.full - 24, Math.min(o.peek + 12, startOffset + dy));
      inspector.style.transform = `translateY(${y}px)`;
      samples.push({ t: performance.now(), y: ev.clientY });
      if (samples.length > 6) samples.shift();
    };
    const onUp = (ev: PointerEvent): void => {
      if (ev.pointerId !== pid) return;
      pid = -1;
      if (!dragging) return;
      dragging = false;
      suppressClick = true;
      window.setTimeout(() => (suppressClick = false), 60);
      const first = samples[0]!;
      const last = samples[samples.length - 1]!;
      const v = (last.y - first.y) / Math.max(1, last.t - first.t); // px/ms
      const projected = currentOffset() + v * 160;
      const o = offsets();
      let best: SheetState = "peek";
      let bestD = Infinity;
      for (const s of ["peek", "half", "full"] as const) {
        const d = Math.abs(o[s] - projected);
        if (d < bestD) {
          bestD = d;
          best = s;
        }
      }
      inspector.classList.remove("is-dragging");
      inspector.style.transform = "";
      setSheet(best);
    };
    for (const el of [$("#sheet-handle"), inspector.querySelector<HTMLElement>(".inspector-head")!]) {
      el.addEventListener("pointerdown", onDown);
      el.addEventListener("pointermove", onMove);
      el.addEventListener("pointerup", onUp);
      el.addEventListener("pointercancel", onUp);
    }
    inspector.addEventListener(
      "click",
      (ev) => {
        if (suppressClick) {
          ev.stopPropagation();
          ev.preventDefault();
        }
      },
      true,
    );
    $("#sheet-handle").addEventListener("click", () => {
      if (!isCompact) return;
      setSheet(sheetState === "peek" ? "half" : sheetState === "half" ? "full" : "peek");
    });
  }

  // View menu (compact)
  const setViewMenu = (open: boolean): void => {
    shell.classList.toggle("view-open", open);
    setAttr(btnViewMenu, "aria-expanded", String(open));
  };
  btnViewMenu.addEventListener("click", () => {
    const open = !shell.classList.contains("view-open");
    if (open && isCompact && sheetState !== "peek") setSheet("peek");
    setViewMenu(open);
  });
  canvasViewport.addEventListener("pointerdown", () => {
    if (shell.classList.contains("view-open")) setViewMenu(false);
  });

  // ---------------------------------------------------------------------------------------------
  // Vector fields & camera
  const ALL_VEC_LAYERS: VectorFieldLayer[] = ["heat", "convection", "magnetic", "acoustic", "stress"];
  const btnVecAll = $<HTMLButtonElement>("#btn-vec-all");

  const syncVectorFieldButtons = (): void => {
    let allOn = true;
    for (const btn of $$<HTMLButtonElement>("[data-vec-field]")) {
      const isOn = scene.isVectorFieldVisible(btn.dataset.vecField as VectorFieldLayer);
      btn.classList.toggle("active", isOn);
      setAttr(btn, "aria-pressed", String(isOn));
      if (!isOn) allOn = false;
    }
    btnVecAll.classList.toggle("active", allOn);
    setAttr(btnVecAll, "aria-pressed", String(allOn));
  };

  const setVectorField = (layer: VectorFieldLayer, visible: boolean): void => {
    scene.setVectorFieldVisible(layer, visible);
    syncVectorFieldButtons();
  };

  for (const btn of $$<HTMLButtonElement>("[data-vec-field]")) {
    btn.addEventListener("click", () => {
      scene.toggleVectorField(btn.dataset.vecField as VectorFieldLayer);
      syncVectorFieldButtons();
    });
  }

  const toggleAllFields = (): void => {
    const anyOff = ALL_VEC_LAYERS.some((l) => !scene.isVectorFieldVisible(l));
    for (const l of ALL_VEC_LAYERS) scene.setVectorFieldVisible(l, anyOff);
    syncVectorFieldButtons();
  };
  btnVecAll.addEventListener("click", toggleAllFields);

  const applyCameraPreset = (preset: CameraPreset): void => {
    // Fit against the layout the chrome is heading to (e.g. a sheet that is about to collapse),
    // not the one it is leaving.
    computeLayout();
    layoutDirty = false;
    Object.assign(insetsCur, insetsTarget);
    scene.setViewInsets({ ...insetsCur });
    for (const b of $$<HTMLButtonElement>(".cam-btn")) {
      const on = b.dataset.cam === preset;
      b.classList.toggle("active", on);
      setAttr(b, "aria-pressed", String(on));
    }
    scene.setCameraPreset(preset);
  };

  for (const btn of $$<HTMLButtonElement>(".cam-btn")) {
    btn.addEventListener("click", () => {
      applyCameraPreset((btn.dataset.cam ?? "overview") as CameraPreset);
      if (isCompact) setViewMenu(false);
    });
  }

  // ---------------------------------------------------------------------------------------------
  // Tools
  const toolButtons = $$<HTMLButtonElement>(".tool-btn[data-tool]");

  const selectToolMode = (mode: ActiveToolMode, autoArm = true): void => {
    activeToolMode = mode;
    for (const b of toolButtons) {
      const on = b.dataset.tool === mode;
      b.classList.toggle("active", on);
      setAttr(b, "aria-checked", String(on));
      b.tabIndex = on ? 0 : -1;
    }
    scene.setToolMode(mode);
    if (mode !== "inspect") dismissCoach();

    // Selecting a perturbation tool arms a sensible working intensity right away.
    if (autoArm) {
      if (mode === "torch" && torchPower < 0.05) {
        torchPower = 0.78;
        setRange(sliderTorchPower, torchPower);
      } else if (mode === "ln2" && ln2Rate < 0.05) {
        ln2Rate = 0.8;
        setRange(sliderLn2Rate, ln2Rate);
      } else if (mode === "oil" && oilRate < 0.05) {
        oilRate = 0.65;
        setRange(sliderOilRate, oilRate);
      } else if (mode === "magnet" && !magnetOn) {
        magnetOn = true;
        chkMagnetOn.checked = true;
      } else if (mode === "air" && airPower < 0.05) {
        airPower = 0.85;
        setRange(sliderAirPower, airPower);
      }
    }
    syncToolsToSim();

    if (mode === "mic") {
      if (panelVisible(currentPanel) && currentPanel !== "metrology") setPanel("metrology", false);
    } else if (mode !== "inspect") {
      if (panelVisible(currentPanel) && currentPanel !== "bench") setPanel("bench", false);
      window.requestAnimationFrame(() => revealSubtool(`subtool-${mode}`));
    }
  };

  for (const btn of toolButtons) {
    btn.tabIndex = btn.classList.contains("active") ? 0 : -1;
    btn.addEventListener("click", () => selectToolMode((btn.dataset.tool ?? "inspect") as ActiveToolMode, true));
    btn.addEventListener("keydown", (ev) => {
      const i = toolButtons.indexOf(btn);
      let j = -1;
      if (ev.key === "ArrowDown" || ev.key === "ArrowRight") j = (i + 1) % toolButtons.length;
      else if (ev.key === "ArrowUp" || ev.key === "ArrowLeft") j = (i - 1 + toolButtons.length) % toolButtons.length;
      if (j >= 0) {
        ev.preventDefault();
        const next = toolButtons[j]!;
        selectToolMode((next.dataset.tool ?? "inspect") as ActiveToolMode, true);
        next.focus();
      }
    });
  }

  const setChipGroupActive = (selector: string, chip: Element | null): void => {
    for (const c of $$(selector)) {
      const on = c === chip;
      c.classList.toggle("active", on);
      setAttr(c, "aria-pressed", String(on));
    }
  };

  const stopAllPerturbations = (): void => {
    torchPower = 0.0;
    setRange(sliderTorchPower, 0);
    ln2Rate = 0.0;
    setRange(sliderLn2Rate, 0);
    oilRate = 0.0;
    setRange(sliderOilRate, 0);
    magnetOn = false;
    chkMagnetOn.checked = false;
    magnetDistMm = 28.0;
    setRange(sliderMagnetDist, 28);
    magnetPosMm = [28.0, -11.0, 2.5];
    airPower = 0.0;
    setRange(sliderAirPower, 0);
    ambientPressureAtm = 1.0;
    setRange(sliderPressureAtm, 1);
    gravityGVec = [0.0, -1.0, 0.0];
    gravityLabel = "1.0 g · 6h down";
    setChipGroupActive("[data-gravity-vec]", document.querySelector("[data-gravity-vec]"));
    selectToolMode("inspect", false);
  };

  btnStopPerturbations.addEventListener("click", stopAllPerturbations);
  btnBannerStop.addEventListener("click", stopAllPerturbations);
  btnBannerOrbit.addEventListener("click", () => selectToolMode("inspect", false));

  // ---------------------------------------------------------------------------------------------
  // Pause, audio, wind, reset, shock
  const pauseIcon = btnPauseSim.querySelector("use")!;
  const setPaused = (paused: boolean): void => {
    simPaused = paused;
    setAttr(btnPauseSim, "aria-pressed", String(paused));
    setAttr(btnPauseSim, "aria-label", paused ? "Resume simulation" : "Pause simulation");
    btnPauseSim.title = paused ? "Resume simulation (Space)" : "Pause simulation (Space)";
    pauseIcon.setAttribute("href", paused ? "#i-play" : "#i-pause");
  };
  const togglePause = (): void => setPaused(!simPaused);
  btnPauseSim.addEventListener("click", togglePause);

  const audioButtons = $$<HTMLButtonElement>('[data-action="toggle-audio"]');
  const syncAudioButtons = (): void => {
    const on = micAudio.isRunning();
    for (const b of audioButtons) {
      setAttr(b, "aria-pressed", String(on));
      const use = b.querySelector("use");
      if (b.id === "btn-toggle-audio") {
        use?.setAttribute("href", on ? "#i-sound" : "#i-mute");
        setAttr(b, "aria-label", on ? "Mute the virtual microphone" : "Listen to the virtual microphone");
      } else {
        use?.setAttribute("href", on ? "#i-mute" : "#i-sound");
        const lbl = b.querySelector(".audio-label");
        if (lbl) setText(lbl, on ? "Mute" : "Listen");
      }
    }
  };
  for (const b of audioButtons) {
    b.addEventListener("click", async () => {
      await micAudio.toggle(sim);
      syncAudioButtons();
      announce(micAudio.isRunning() ? "Microphone audio on" : "Microphone audio muted");
    });
  }

  const windCrown = (): void => {
    sim.wind_crown(1.5);
    announce("Crown wound one and a half turns");
  };
  const shock = (): void => {
    sim.trigger_mechanical_shock(500.0);
    announce("500 g drop shock delivered");
  };
  $("#btn-wind-crown").addEventListener("click", windCrown);
  btnShock500g.addEventListener("click", shock);

  const resetWatchAll = (): void => {
    cancelFastForward();
    stopAllPerturbations();
    sim.reset();
    sim.set_compensated_hairspring(chkInvarSpring.checked);
    syncToolsToSim();
    timegrapher.reset();
  };
  $("#btn-reset-watch").addEventListener("click", resetWatchAll);

  for (const b of $$<HTMLButtonElement>("[data-action]")) {
    const a = b.dataset.action;
    if (a === "shock") b.addEventListener("click", shock);
    else if (a === "wind") b.addEventListener("click", windCrown);
    else if (a === "reset") b.addEventListener("click", () => {
      resetWatchAll();
      announce("Watch reset to 20 °C");
    });
  }

  // Display toggles
  sliderExploded.addEventListener("input", () => scene.setExplodedFactor(Number(sliderExploded.value)));
  chkOpenCase.addEventListener("change", () => scene.setOpenCaseMode(chkOpenCase.checked));
  chkInvarSpring.addEventListener("change", () => sim.set_compensated_hairspring(chkInvarSpring.checked));
  chkSlowMo.addEventListener("change", () => (slowMo = chkSlowMo.checked));
  chkThermalOverlay.addEventListener("change", () => scene.setThermalOverlayMode(chkThermalOverlay.checked));
  selQuality.value = qualityChoice;
  selQuality.addEventListener("change", () => {
    const v = selQuality.value;
    storeSet(STORE_QUALITY, v);
    tier = v === "high" || v === "medium" || v === "low" ? v : detectedTier;
    scene.setQualityTier(tier);
    applyTierClass();
  });

  // Perturbation controls
  sliderTorchPower.addEventListener("input", () => {
    torchPower = Number(sliderTorchPower.value);
    syncToolsToSim();
  });
  for (const chip of $$<HTMLButtonElement>("[data-torch-aim]")) {
    chip.addEventListener("click", () => {
      setChipGroupActive("[data-torch-aim]", chip);
      const parts = (chip.dataset.torchAim ?? "0,-24,2").split(",").map(Number);
      torchPosMm = [parts[0] ?? 0, parts[1] ?? -24, parts[2] ?? 2];
      if (torchPower < 0.05) {
        torchPower = 0.78;
        setRange(sliderTorchPower, torchPower);
      }
      syncToolsToSim();
    });
  }
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
    startFastForward(
      "Flooding the hairspring",
      [{ steps: 16, before: () => sim.set_oil_pour(1.0, 0.0, -11.0, 2.0) }],
      () => sim.set_oil_pour(oilRate, oilPosMm[0], oilPosMm[1], oilPosMm[2]),
    );
  });
  $("#btn-clean-oil").addEventListener("click", () => {
    oilRate = 0.0;
    setRange(sliderOilRate, 0);
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
  sliderAirPower.addEventListener("input", () => {
    airPower = Number(sliderAirPower.value);
    syncToolsToSim();
  });
  sliderPressureAtm.addEventListener("input", () => {
    ambientPressureAtm = Number(sliderPressureAtm.value);
    syncToolsToSim();
  });
  for (const chip of $$<HTMLButtonElement>("[data-gravity-vec]")) {
    chip.addEventListener("click", () => {
      setChipGroupActive("[data-gravity-vec]", chip);
      const raw = (chip.dataset.gravityVec ?? "0,-1,0,1.0 g").split(",");
      gravityGVec = [Number(raw[0] ?? 0), Number(raw[1] ?? -1), Number(raw[2] ?? 0)];
      gravityLabel = raw[3] ?? "1.0 g";
      syncToolsToSim();
    });
  }

  // Microphone
  sliderMicGain.addEventListener("input", () => {
    micGain = Number(sliderMicGain.value);
    syncToolsToSim();
  });
  for (const chip of $$<HTMLButtonElement>("[data-mic-pos]")) {
    chip.addEventListener("click", () => {
      setChipGroupActive("[data-mic-pos]", chip);
      const parts = (chip.dataset.micPos ?? "0,-11,11").split(",").map(Number);
      micPosMm = [parts[0] ?? 0, parts[1] ?? -11, parts[2] ?? 11];
      syncToolsToSim();
    });
  }

  // ---------------------------------------------------------------------------------------------
  // Fast-forward queue: long experiments advance lab time in time-budgeted chunks per frame so the
  // page stays responsive and you watch the physics unfold instead of waiting on a frozen tab.
  interface FfPhase {
    steps: number;
    dt?: number;
    before?: () => void;
  }
  interface FfJob {
    title: string;
    phases: FfPhase[];
    phase: number;
    stepInPhase: number;
    done: number;
    total: number;
    labS: number;
    totalLabS: number;
    after?: () => void;
  }
  let ff: FfJob | null = null;
  let heavyPollNow = false;
  /** Engine time the fast-forward spent in the last tick, and an EMA of the rest of a frame. */
  let ffSpentMs = 0;
  let frameRestMs = 8;

  function startFastForward(title: string, phases: FfPhase[], after?: () => void): void {
    cancelFastForward();
    const total = phases.reduce((n, p) => n + p.steps, 0);
    const totalLabS = phases.reduce((n, p) => n + p.steps * (p.dt ?? 0.05), 0);
    ff = { title, phases, phase: 0, stepInPhase: 0, done: 0, total, labS: 0, totalLabS, after };
    if (simPaused) setPaused(false);
    setText(ffTitle, title);
    ffOverlay.classList.remove("hidden");
    announce(`${title}: fast-forwarding ${totalLabS.toFixed(1)} seconds of lab time`);
  }

  function cancelFastForward(): void {
    if (!ff) return;
    ff = null;
    ffOverlay.classList.add("hidden");
  }

  const runFastForward = (): void => {
    const job = ff;
    if (!job) return;
    // Give the engine ~60 % of each frame: 1.5x whatever rendering and the page take, so a slow GPU
    // (or a software rasterizer) does not starve the fast-forward, while frames stay under ~0.2 s.
    const budgetMs = Math.min(120, Math.max(isCompact ? 12 : 20, 1.5 * frameRestMs));
    const t0 = performance.now();
    do {
      const ph = job.phases[job.phase];
      if (!ph) break;
      if (job.stepInPhase === 0 && ph.before) ph.before();
      if (job.stepInPhase < ph.steps) {
        const dt = ph.dt ?? 0.05;
        sim.step(dt);
        job.stepInPhase++;
        job.done++;
        job.labS += dt;
      }
      if (job.stepInPhase >= ph.steps) {
        job.phase++;
        job.stepInPhase = 0;
      }
    } while (performance.now() - t0 < budgetMs);
    ffSpentMs = performance.now() - t0;
    heavyPollNow = true;
    const frac = job.total > 0 ? job.done / job.total : 1;
    ffBarFill.style.width = `${(frac * 100).toFixed(1)}%`;
    setText(ffDetail, `${job.labS.toFixed(1)} / ${job.totalLabS.toFixed(1)} s lab time`);
    if (job.phase >= job.phases.length) {
      ff = null;
      ffOverlay.classList.add("hidden");
      job.after?.();
      announce(`${job.title} complete`);
    }
  };
  $("#btn-ff-cancel").addEventListener("click", () => {
    cancelFastForward();
    announce("Stopped; the watch continues in real time from here");
  });

  // ---------------------------------------------------------------------------------------------
  // Lay the watch in a chosen pose (gravity in the watch frame) and light the matching chip.
  const setPose = (vec: Vec3, label: string): void => {
    gravityGVec = vec;
    gravityLabel = label;
    setChipGroupActive("[data-gravity-vec]", document.querySelector(`[data-gravity-vec^="${vec.join(",")},"]`));
  };

  // Guided experiments. Durations come from the engine's own time scales (native calibration runs
  // in engine-watch/src/lib.rs, `presets_*`): the whole case has to soak up heat before gold melts,
  // and LN₂ boils in a Leidenfrost film, so these are tens of lab seconds, not a few.
  const FF_DT = 0.1;
  const ffSteps = (labS: number): number => Math.round(labS / FF_DT);
  const runPreset = (preset: string): void => {
    dismissCoach();
    cancelFastForward();
    // Each experiment starts from a clean view; the ones that need an overlay turn it on below.
    for (const layer of ["heat", "convection", "magnetic", "acoustic", "stress"] as const) {
      scene.setVectorFieldVisible(layer, false);
    }
    syncVectorFieldButtons();
    switch (preset) {
      case "chronometer":
        resetWatchAll();
        applyCameraPreset("overview");
        break;

      case "thermal-drift":
        // ~0.1 kW absorbed at the 6h rim: the balance warms ~120 K in 10 s and the watch loses.
        resetWatchAll();
        torchPosMm = [0.0, -22.0, 6.0];
        torchPower = 0.25;
        setRange(sliderTorchPower, torchPower);
        syncToolsToSim();
        applyCameraPreset("tourbillon");
        startFastForward("Thermal rate drift", [{ steps: ffSteps(10), dt: FF_DT }]);
        break;

      case "melt-gold":
        // Full flame on the 6h rim: the case soaks for about a minute before the band slumps.
        resetWatchAll();
        torchPosMm = [0.0, -25.0, 1.5];
        torchPower = 1.0;
        setRange(sliderTorchPower, torchPower);
        syncToolsToSim();
        applyCameraPreset("overview");
        startFastForward("Melting the 18k case", [{ steps: ffSteps(75), dt: FF_DT }]);
        break;

      case "melt-then-ln2": {
        const curSnap = decodeWatchSnapshot(sim.snapshot());
        const phases: FfPhase[] = [];
        if ((curSnap.nodes[3]?.temperatureK ?? 0) < 950) {
          phases.push({
            steps: ffSteps(70),
            dt: FF_DT,
            before: () => {
              torchPosMm = [0.0, -25.0, 1.5];
              torchPower = 1.0;
              setRange(sliderTorchPower, torchPower);
              syncToolsToSim();
            },
          });
        }
        phases.push({
          steps: ffSteps(30),
          dt: FF_DT,
          before: () => {
            // Lay it dial-up so the nitrogen pools in the tourbillon aperture instead of running off.
            torchPower = 0.0;
            setRange(sliderTorchPower, 0);
            setPose([0, 0, -1], "1.0 g · dial up");
            ln2PosMm = [0.0, -11.0, 3.0];
            ln2Rate = 1.0;
            setRange(sliderLn2Rate, ln2Rate);
            syncToolsToSim();
          },
        });
        applyCameraPreset("overview");
        startFastForward("Melt, then quench with LN₂", phases);
        break;
      }

      case "cryo-freeze":
        // Dial-up so the LN₂ pools over the balance; the oil congeals near 190 K and seizes it.
        resetWatchAll();
        setPose([0, 0, -1], "1.0 g · dial up");
        ln2PosMm = [0.0, -11.0, 3.0];
        ln2Rate = 1.0;
        setRange(sliderLn2Rate, ln2Rate);
        syncToolsToSim();
        applyCameraPreset("tourbillon");
        startFastForward("Cryogenic seizure", [{ steps: ffSteps(50), dt: FF_DT }]);
        break;

      case "oil-flood":
        resetWatchAll();
        oilPosMm = [0.0, -11.0, 2.0];
        oilRate = 0.85;
        setRange(sliderOilRate, oilRate);
        syncToolsToSim();
        applyCameraPreset("tourbillon");
        startFastForward("Flooding the hairspring", [{ steps: ffSteps(7), dt: FF_DT }]);
        break;

      case "magnet-attack":
        // From a fresh running watch, so the stall it shows is the magnet's and nothing earlier.
        resetWatchAll();
        magnetOn = true;
        chkMagnetOn.checked = true;
        magnetDistMm = 13.0;
        setRange(sliderMagnetDist, 13);
        magnetPosMm = [13.0, -11.0, 2.5];
        syncToolsToSim();
        setVectorField("magnetic", true);
        applyCameraPreset("tourbillon");
        break;

      case "air-blow-gold":
        // Melt the rim, kill the flame, then a compressed-air jet: forced convection pulls heat out
        // and its crossflow drags on the balance rim.
        resetWatchAll();
        applyCameraPreset("overview");
        startFastForward(
          "Air-jet quench",
          [
            {
              steps: ffSteps(70),
              dt: FF_DT,
              before: () => {
                torchPosMm = [0.0, -25.0, 1.5];
                torchPower = 1.0;
                setRange(sliderTorchPower, torchPower);
                syncToolsToSim();
              },
            },
            {
              steps: ffSteps(6),
              dt: FF_DT,
              before: () => {
                torchPower = 0.0;
                setRange(sliderTorchPower, 0);
                airPosMm = [0.0, -24.0, 2.5];
                airPower = 0.92;
                setRange(sliderAirPower, airPower);
                syncToolsToSim();
                setVectorField("convection", true);
              },
            },
          ],
        );
        break;

      case "vacuum-chamber":
        resetWatchAll();
        ambientPressureAtm = 0.0;
        setRange(sliderPressureAtm, 0);
        syncToolsToSim();
        setVectorField("convection", true);
        applyCameraPreset("tourbillon");
        startFastForward("Evacuating the bell jar", [{ steps: ffSteps(10), dt: FF_DT }]);
        break;

      case "shock-centrifuge": {
        resetWatchAll();
        setPose([0, -4, 0], "4.0 g · centrifuge");
        syncToolsToSim();
        sim.trigger_mechanical_shock(500.0);
        setVectorField("stress", true);
        applyCameraPreset("overview");
        startFastForward("Drop shock under 4 g", [{ steps: ffSteps(2), dt: FF_DT }]);
        break;
      }
    }
  };

  for (const btn of $$<HTMLButtonElement>("[data-preset]")) {
    btn.addEventListener("click", () => {
      if (simPaused) setPaused(false);
      for (const b of $$(".preset-btn")) {
        const on = b === btn;
        b.classList.toggle("active", on);
        setAttr(b, "aria-pressed", String(on));
      }
      // Let the watch be seen: drop the sheet on small screens (before the preset frames the camera).
      if (isCompact) setSheet("peek");
      runPreset(btn.dataset.preset ?? "chronometer");
    });
  }

  // ---------------------------------------------------------------------------------------------
  // Coach & help
  function dismissCoach(): void {
    if (coach.classList.contains("hidden")) return;
    coach.classList.add("hidden");
    storeSet(STORE_COACH, "done");
    layoutDirty = true;
    refitPending = true;
  }
  // Shown under the boot splash so the first layout already frames the watch above it.
  if (storeGet(STORE_COACH) !== "done") coach.classList.remove("hidden");
  $("#btn-coach-dismiss").addEventListener("click", dismissCoach);
  $("#btn-coach-experiments").addEventListener("click", () => {
    dismissCoach();
    setPanel("experiments", true);
  });

  const openHelp = (): void => {
    if (!helpDialog.open) helpDialog.showModal();
  };
  $("#btn-help").addEventListener("click", openHelp);
  $("#btn-help-close").addEventListener("click", () => helpDialog.close());
  helpDialog.addEventListener("click", (ev) => {
    if (ev.target === helpDialog) helpDialog.close(); // backdrop click
  });

  // ---------------------------------------------------------------------------------------------
  // Keyboard shortcuts (ignored while typing or with modifiers)
  window.addEventListener("keydown", (ev) => {
    if (helpDialog.open) return;
    const target = ev.target as HTMLElement | null;
    const tag = target?.tagName;
    if (tag === "INPUT" || tag === "TEXTAREA" || tag === "SELECT") return;
    if (ev.metaKey || ev.ctrlKey || ev.altKey) return;
    if (target?.getAttribute("role") === "tab" && ev.key.startsWith("Arrow")) return;

    const key = ev.key.toLowerCase();
    if (ev.key === "Escape") {
      if (shell.classList.contains("view-open")) setViewMenu(false);
      else if (activeToolMode !== "inspect") selectToolMode("inspect", false);
      else if (isCompact && sheetState !== "peek") setSheet("peek");
      else if (theater) toggleTheaterMode();
      return;
    }
    if (ev.key === "?") {
      openHelp();
      return;
    }
    // Camera: arrows pan (Shift orbits), + and − zoom, 0 or Home reframes. Only when no control
    // has focus, so arrows still move sliders, tabs and the tool rail.
    const free = !ev.defaultPrevented && (target === document.body || target === canvasViewport || target?.id === "stage");
    if (free && ev.key.startsWith("Arrow")) {
      ev.preventDefault();
      const dir = ev.key === "ArrowLeft" ? [1, 0] : ev.key === "ArrowRight" ? [-1, 0] : ev.key === "ArrowUp" ? [0, 1] : [0, -1];
      if (ev.shiftKey) scene.nudgeCamera("orbit", dir[0]! * 0.12, dir[1]! * 0.09);
      else scene.nudgeCamera("pan", dir[0]! * 36, dir[1]! * 36);
      return;
    }
    if (ev.key === "+" || ev.key === "=") {
      scene.nudgeCamera("zoom", 1.18);
      return;
    }
    if (ev.key === "-" || ev.key === "_") {
      scene.nudgeCamera("zoom", 1 / 1.18);
      return;
    }
    if (ev.key === "0" || ev.key === "Home") {
      scene.resetView();
      return;
    }
    const toolKeys: Record<string, ActiveToolMode> = { "1": "inspect", "2": "torch", "3": "ln2", "4": "oil", "5": "magnet", "6": "air", "7": "mic" };
    const tool = toolKeys[key];
    if (tool) {
      selectToolMode(tool, tool !== "inspect" && tool !== "mic");
      return;
    }
    const fieldKeys: Record<string, VectorFieldLayer> = { h: "heat", c: "convection", b: "magnetic", a: "acoustic", g: "stress" };
    const field = fieldKeys[key];
    if (field) {
      scene.toggleVectorField(field);
      syncVectorFieldButtons();
      return;
    }
    if (key === "v") toggleAllFields();
    else if (key === "s") shock();
    else if (key === "x") stopAllPerturbations();
    else if (key === "w") windCrown();
    else if (key === "r") resetWatchAll();
    else if (key === "f") toggleTheaterMode();
    else if (key === "i") {
      if (isCompact) setSheet(sheetState === "peek" ? "half" : "peek");
      else setInspectorOpen(!inspectorOpen);
    } else if (ev.code === "Space") {
      ev.preventDefault();
      togglePause();
    }
  });

  // ---------------------------------------------------------------------------------------------
  // HUD
  let contactAccum = 0;
  let lastEsc: Float64Array | null = null;
  let frameCounter = 0;

  // One word for the compact pill, so certified and perturbed differ in text, not colour alone.
  const STATUS_SHORT: Readonly<Record<string, string>> = {
    "CHRONOMETER CERTIFIED": "Certified",
    "MULTIPHYSICS PERTURBED": "Perturbed",
    "MAINSPRING UNWOUND": "Unwound",
    "FOULED BY LIQUID MOLTEN GOLD": "Molten gold",
    "ENCASED IN SOLIDIFIED GOLD": "In gold",
    "THERMAL / CRYOGENIC CLEARANCE LOSS": "No clearance",
    "PALLET FORK MAGNETICALLY PINNED": "Magnetized",
    "VISCOUS OIL FLOOD OVERDAMPING": "Oil flood",
  };
  const parseStatus = (label: string): [string, string] => {
    const [head, detail] = label.split(" — ");
    const cap = (s: string): string => (s ? s.charAt(0) + s.slice(1).toLowerCase() : "");
    return [cap(head ?? label), cap(detail ?? "")];
  };

  const updateHud = (snap: WatchSnapshot): void => {
    setText(hudWatchTime, formatHorologicalTime(snap.watchDisplayedTimeS));
    setText(hudLabTime, formatHorologicalTime(snap.labTimeS));
    const driftS = snap.watchDisplayedTimeS - snap.labTimeS;
    setText(hudTimeDrift, `${signed(driftS, 2)} s`);
    setAttr(hudTimeDrift, "data-tone", Math.abs(driftS) < 0.5 ? "ok" : Math.abs(driftS) < 5 ? "warn" : "bad");

    const tg = timegrapher.readings(snap.labTimeS);
    // A stopped balance has no rate or beat error; say so instead of extrapolating a fit.
    const stopped = tg.stopped || (snap.statusCode >= 2 && tg.rateSDay === null);
    // Until a dozen beats are in, there is no measured rate: the free-oscillator figure would
    // leave out the escapement's own error (which regulation cancels), so show a dash instead.
    const engineRate = snap.totalBeats >= 12 ? snap.measuredRateSPerDay : null;
    const rate = tg.rateSDay ?? engineRate;
    const beatErr = tg.beatErrorMs ?? snap.measuredBeatErrorMs;
    const rateStr = stopped
      ? "stopped"
      : rate === null
        ? "— s/d"
        : `${signed(rate, Math.abs(rate) >= 100 ? 0 : 1)}${Math.abs(rate) >= 10_000 ? "" : " s/d"}`;
    setText(hudRateSDay, rateStr);
    setAttr(
      hudRateSDay,
      "data-tone",
      stopped ? "bad" : rate === null ? "" : rate >= -4 && rate <= 6 ? "ok" : Math.abs(rate) < 60 ? "warn" : "bad",
    );
    const amp = snap.balanceAmplitudeDeg;
    setText(hudAmp, `${amp.toFixed(0)}°`);
    setAttr(hudAmp, "data-tone", amp >= 220 && amp <= 320 ? "ok" : amp >= 150 ? "warn" : "bad");
    const beatStr = stopped ? "—" : `${beatErr.toFixed(1)} ms`;
    setText(hudBeatError, beatStr);
    setAttr(hudBeatError, "data-tone", stopped ? "bad" : beatErr < 0.6 ? "ok" : beatErr < 1.5 ? "warn" : "bad");
    setText(tgRate, rateStr);
    setText(tgAmp, `${amp.toFixed(0)}°`);
    setText(tgBeatErr, beatStr);
    setText(
      tgSource,
      tg.source === "impacts" ? "lock impacts" : tg.source === "crossings" ? "balance timing" : "engine readout",
    );

    const label = sim.status_label();
    let [head, detail] = simPaused ? ["Paused", "Simulation frozen"] : parseStatus(label);
    let title = simPaused ? "Paused" : label;
    let short = simPaused ? "Time frozen" : (STATUS_SHORT[label.split(" — ")[1] ?? ""] ?? detail);
    // With the magnet gone the lever stays where the field threw it: across its banking, the guard
    // pin riding the safety roller. Say so, since degaussing or a knock will not bring it back.
    if (!simPaused && label.startsWith("STALLED — PALLET FORK MAGNETICALLY") && snap.bFieldEscapementMt < 1) {
      detail = "Lever overbanked · reset to recover";
      short = "Overbanked";
      title = `${label}. The field is gone, but it dragged the lever across its banking; the guard pin now rides the safety roller and holds it there. A watchmaker would have to open the case: press Reset.`;
    }
    setText(statusHead, head);
    setText(statusDetail, detail);
    setText(statusShort, short);
    hudStatusPill.title = title;
    const cls =
      "status-pill " + (simPaused ? "status-warn" : snap.statusCode === 0 ? "status-ok" : snap.statusCode === 1 ? "status-warn" : "status-danger");
    if (hudStatusPill.className !== cls) hudStatusPill.className = cls;

    const reserveH = lastEsc ? (lastEsc[13] ?? snap.mainspringWind * 42) : snap.mainspringWind * 42;
    const reservePct = Math.round(snap.mainspringWind * 100);
    setText(hudPowerReservePct, `${reservePct}% · ${reserveH.toFixed(1)} h`);
    setAttr(
      hudPowerReservePct,
      "title",
      reserveH < 0.05 && reservePct > 5
        ? `${reservePct}% of the slip turns are still wound, but the spring can no longer hold 150° of amplitude: no power reserve (a mainspring held above ~700 K loses its set)`
        : `${reservePct}% of the slip turns wound · ${reserveH.toFixed(1)} h until the amplitude falls below ~150°`,
    );
    hudPowerReserveBar.style.width = `${reservePct}%`;
    setText(hudBarrelTorque, `${(snap.mainspringTorqueNm * 1000).toFixed(2)} mN·m`);
    setText(hudBalanceAmp, `${amp.toFixed(1)}°`);
    setText(hudPivotClearance, `${snap.minPivotClearanceUm.toFixed(2)} µm`);
    setText(hudTotalBeats, Math.floor(snap.totalBeats).toLocaleString("en-US"));

    setText(lblTorchPower, `${Math.round(torchPower * 100)}% · ${snap.blowtorchHeatingPowerW.toFixed(0)} W`);
    setText(lblLn2Rate, `${Math.round(ln2Rate * 100)}% · ${snap.ln2BoilingPowerW.toFixed(0)} W`);
    setText(lblOilRate, `${(snap.nodes[9]?.oilFilmUm ?? 0).toFixed(1)} µm film`);
    setText(lblBField, `${snap.bFieldEscapementMt.toFixed(1)} mT`);
    setText(lblRemanence, `${snap.remanentMagnetizationMt.toFixed(1)} mT`);
    setText(lblAirPower, `${Math.round(airPower * 100)}% · ${Math.round(airPower * 150)} m/s`);
    setText(lblChamberEnv, `${ambientPressureAtm.toFixed(2)} atm · ${gravityLabel.split(" · ")[0] ?? gravityLabel}`);

    const n6h = snap.nodes[3];
    if (n6h) {
      const caseT = `${n6h.temperatureK.toFixed(1)} K · ${(n6h.temperatureK - 273.15).toFixed(1)} °C`;
      const melt = `${(n6h.meltFraction * 100).toFixed(1)}% · slump ${(n6h.plasticSlump * 100).toFixed(0)}%`;
      setText(ovCaseTemp, caseT);
      setText(ovGoldMelt, melt);
      bindText("caseTemp", caseT);
      bindText("goldMelt", melt);
    }
    let liq = 0;
    let sol = 0;
    for (let i = 4; i < 12; i++) {
      liq += snap.nodes[i]?.moltenGoldG ?? 0;
      sol += snap.nodes[i]?.frozenGoldG ?? 0;
    }
    const gold = `${liq.toFixed(3)} / ${sol.toFixed(3)} g`;
    const bField = `${snap.bFieldEscapementMt.toFixed(1)} mT · J ${snap.remanentMagnetizationMt.toFixed(snap.remanentMagnetizationMt >= 100 ? 0 : 1)} mT`;
    setText(ovFouledGold, gold);
    setText(ovBField, bField);
    bindText("fouledGold", gold);
    bindText("bField", bField);

    if (panelVisible("metrology")) {
      setText(
        lblMicSpl,
        `${sim.mic_spl_db().toFixed(1)} dB SPL · ${sim.mic_loudness_sone().toFixed(2)} sone · ${sim.mic_sharpness_acum().toFixed(2)} acum`,
      );
    }

    if (panelVisible("movement")) {
      if (frameCounter % 4 === 0) setText(hudDigest, `${sim.state_digest_hex().slice(0, 16)}…`);
      setText(chrSens, `${signed(sim.rate_sensitivity_s_day_per_k(), 2)} s/d/K`);
      const lo = sim.rate_interval_lo_s_day();
      const hi = sim.rate_interval_hi_s_day();
      const dg = Math.max(Math.abs(lo), Math.abs(hi)) >= 100 ? 0 : 1;
      setText(chrInterval, `${signed(lo, dg)} … ${signed(hi, dg)} s/d`);
      setText(chrCosc, sim.cosc_e_value().toFixed(2));
      setText(chrDrag, `${sim.balance_aero_drag_s_inv().toExponential(2)} s⁻¹`);
      if (lastEsc) {
        const v = lastEsc;
        setText(escFields.impulse, `${(v[7] ?? 0).toFixed(2)} µJ`);
        setText(escFields.eff, `${((v[8] ?? 0) * 100).toFixed(0)}%`);
        setText(escFields.lift, `${(v[9] ?? 0).toFixed(1)}°`);
        setText(escFields.q, (v[10] ?? 0).toFixed(0));
        setText(escFields.torque, `${(v[11] ?? 0).toFixed(2)} µN·m`);
        setText(escFields.friction, `${(v[14] ?? 0).toFixed(1)} nN·m`);
        setText(escFields.air, `${((v[15] ?? 0) * 100).toFixed(0)}%`);
        setText(escFields.reserve, `${(v[13] ?? 0).toFixed(1)} h`);
      }
      for (const el of contactEls) {
        el.classList.toggle("on", (contactAccum & (1 << Number(el.dataset.bit))) !== 0);
      }
    }
    contactAccum = 0;
    frameCounter++;
  };

  if (!hasEscState) {
    escNote.classList.remove("hidden");
    setText(escSource, "not reported");
    $("#esc-kv").classList.add("hidden");
    $(".contacts").classList.add("hidden");
  } else {
    setText(escSource, "engine live");
  }

  // ---------------------------------------------------------------------------------------------
  // Adaptive resolution: hold the frame rate by trading pixels, never physics.
  const gov = { ema: 16.7, slowSince: 0, fastSince: 0, lastChange: 0, fps: 60 };
  // ?debug exposes the live objects for profiling from the console.
  if (new URLSearchParams(window.location.search).has("debug")) {
    (window as unknown as { __pw: unknown }).__pw = { scene, sim, gov };
  }
  const governDpr = (frameMs: number, now: number): void => {
    if (frameMs > 250) return; // tab switch or breakpoint
    gov.ema += (frameMs - gov.ema) * 0.06;
    gov.fps = 1000 / gov.ema;
    if (now - gov.lastChange < 1500) return;
    if (gov.ema > 23.5) {
      gov.fastSince = 0;
      if (!gov.slowSince) gov.slowSince = now;
      if (now - gov.slowSince > 1200) {
        const before = scene.getPixelRatio();
        const after = scene.setPixelRatio(before - 0.15);
        if (after !== before) layoutDirty = true;
        gov.lastChange = now;
        gov.slowSince = 0;
      }
    } else if (gov.ema < 17.8) {
      gov.slowSince = 0;
      if (!gov.fastSince) gov.fastSince = now;
      if (now - gov.fastSince > 5000) {
        scene.setPixelRatio(scene.getPixelRatio() + 0.1);
        gov.lastChange = now;
        gov.fastSince = 0;
      }
    } else {
      gov.slowSince = 0;
      gov.fastSince = 0;
    }
  };

  // ---------------------------------------------------------------------------------------------
  // Main loop
  let rafId = 0;
  let lastFrameMs = performance.now();
  let lastHudMs = 0;
  let lastHeavyMs = -1e9;
  let lastPlotMs = 0;
  let lastLedgerMs = 0;
  let firstFrame = true;

  const tick = (nowMs: number): void => {
    rafId = window.requestAnimationFrame(tick);
    const rawMs = nowMs - lastFrameMs;
    lastFrameMs = nowMs;
    const dtReal = Math.min(0.05, Math.max(0.001, rawMs / 1000));
    // The resolution governor and the fast-forward budget both look at the frame minus the
    // engine work the previous tick spent fast-forwarding.
    const restMs = Math.max(0, rawMs - ffSpentMs);
    ffSpentMs = 0;
    if (restMs < 500) frameRestMs += (restMs - frameRestMs) * 0.2;
    governDpr(restMs, nowMs);

    if (ff) {
      runFastForward();
    } else if (!simPaused) {
      sim.step(slowMo ? dtReal * 0.15 : dtReal);
    }
    if (!simPaused || ff) micAudio.pump(sim, dtReal);

    if (hasDrain) {
      const drained = ext.timegrapher_drain!();
      if (drained.length >= 3) timegrapher.ingestImpacts(drained);
    }
    const snap = decodeWatchSnapshot(sim.snapshot());
    if (!hasDrain) timegrapher.ingestBalance(snap.labTimeS, snap.balanceAngleRad, snap.balanceOmegaRadS);
    if (hasEscState) {
      lastEsc = ext.escapement_state!();
      const mask = lastEsc[6] ?? 0;
      scene.setEscapementContacts(mask);
      contactAccum |= mask;
    }

    // The continuum fields (case shell, spatial grid, mainspring coils) are polled at 25 Hz; the
    // scene keeps the last upload on the GPU and only re-uploads when the data actually changed.
    const heavy = heavyPollNow || nowMs - lastHeavyMs >= 40;
    if (heavy) {
      lastHeavyMs = nowMs;
      heavyPollNow = false;
    }
    stepInsets(dtReal);
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
      micVisible: activeToolMode === "mic" || micAudio.isRunning() || panelVisible("metrology"),
      caseShellData: heavy ? sim.case_shell_snapshot() : undefined,
      spatialGridData: heavy ? sim.spatial_grid_snapshot() : undefined,
      mainspringCoilsMm: heavy ? sim.mainspring_coils_mm() : undefined,
      clickLeverAngleRad: sim.click_lever_angle_rad(),
      crownWheelAngleRad: sim.crown_wheel_angle_rad(),
    });

    if (nowMs - lastHudMs >= 100) {
      lastHudMs = nowMs;
      updateHud(snap);
      setText(hudPerf, `${gov.fps.toFixed(0)} fps · ${scene.getPixelRatio().toFixed(2)}× · ${tier}`);
    }
    if (panelVisible("metrology")) {
      if (nowMs - lastPlotMs >= 33) {
        lastPlotMs = nowMs;
        micAudio.renderOscilloscope(canvasScope, micPosMm);
        timegrapher.draw(canvasTimegrapher, snap);
      }
      if (nowMs - lastLedgerMs >= 250) {
        lastLedgerMs = nowMs;
        renderNodeLedgerTable(tbodyLedger, snap);
      }
    }

    if (firstFrame) {
      firstFrame = false;
      splash.classList.add("is-done");
      window.setTimeout(() => splash.remove(), 700);
    }
  };

  // Pause everything while the page is hidden: no GPU work, no audio, no physics catch-up burst.
  document.addEventListener("visibilitychange", () => {
    if (document.hidden) {
      if (rafId) window.cancelAnimationFrame(rafId);
      rafId = 0;
      void micAudio.suspend();
    } else if (!rafId) {
      lastFrameMs = performance.now();
      void micAudio.resume();
      rafId = window.requestAnimationFrame(tick);
    }
  });

  if (import.meta.env.DEV) {
    // Dev-only handle for console / automated checks.
    (window as unknown as { __pocketWatch?: unknown }).__pocketWatch = { sim, scene, timegrapher, micAudio };
  }

  setPanel("bench", false);
  syncVectorFieldButtons();
  syncAudioButtons();
  updatePerturbationBannerAndBadges();
  rafId = window.requestAnimationFrame(tick);
}

void boot();
