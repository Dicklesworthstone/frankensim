// 9 Guided Physics Experiments for the Steam Turbine Simulation.
// Every preset injects genuine boundary conditions or physical faults into the
// 12-component fs-system convergence graph; all responses are emergent.

import type { SteamTurbineSim } from "./wasm-pkg/fs_steam_turbine";

export type CutawayMode = "half" | "open" | "closed";
export type ShaderMode = "pbr" | "temp" | "stress" | "pressure" | "mach" | "wetcreep";
export type CameraPreset = "overview" | "lp_stage" | "hp_valves" | "seals" | "exhaust" | "top_open";

export interface ExperimentSpec {
  id: string;
  title: string;
  tag: string;
  summary: string;
  whatToWatch: string;
  camera: CameraPreset;
  cutaway: CutawayMode;
  shader: ShaderMode;
  rightTab: "thermo" | "rotor" | "audio" | "ledger";
  apply: (sim: SteamTurbineSim) => void;
}

export const EXPERIMENTS: ExperimentSpec[] = [
  {
    id: "nominal",
    title: "Nominal 50 MW CHP Design Point",
    tag: "DESIGN",
    summary:
      "110 bar / 525 °C live steam synchronized at 3,000 rpm with 12-bar CHP process extraction and 7.5 kPa vacuum condensing.",
    whatToWatch:
      "Observe ~49 MW output, ~41% CHP utilization vs ~28% electric efficiency (rising >53% in CHP Swing), ~93% LP exhaust dryness, and exact fs-system::Ledger energy closure.",
    camera: "overview",
    cutaway: "half",
    shader: "pbr",
    rightTab: "thermo",
    apply: (sim) => {
      sim.reset(false);
    },
  },
  {
    id: "load_rejection",
    title: "Full Load Rejection (11 kV Breaker Trip)",
    tag: "TRANSIENT",
    summary:
      "Opens the 11 kV synchronous generator tie breaker at full 50 MW load, unloading 160 kN·m of electromagnetic braking torque instantaneously.",
    whatToWatch:
      "Watch shaft speed surge past 3,060 rpm as the 4% droop + fs-control SpeedController slams the HP control valve shut to catch the rotor.",
    camera: "hp_valves",
    cutaway: "half",
    shader: "pbr",
    rightTab: "thermo",
    apply: (sim) => {
      sim.reset(false);
      sim.set_grid_breaker(false);
    },
  },
  {
    id: "air_ingress",
    title: "Condenser Air In-Leakage (Colburn-Hougen)",
    tag: "MULTIPHASE",
    summary:
      "Injects 18% non-condensable air into the vacuum condenser shell, blanketing the 3,600 tube bundle via fs-multiphase CondensationModelV1.",
    whatToWatch:
      "Shell-side condensation HTC collapses, condenser pressure climbs from ~7.5 kPa toward >25 kPa, Mollier end-point lifts, and MW output falls.",
    camera: "exhaust",
    cutaway: "half",
    shader: "pbr",
    rightTab: "thermo",
    apply: (sim) => {
      sim.reset(false);
      sim.set_condenser_air_fraction(0.18);
    },
  },
  {
    id: "oil_whip",
    title: "Sub-Synchronous Oil Whip & Unbalance",
    tag: "ROTORDYNAMICS",
    summary:
      "Degrades tilting-pad bearing preload & oil supply pressure while adding 75,000 g·mm blade fouling unbalance on the 4-DOF gyroscopic rotor.",
    whatToWatch:
      "4-DOF fs-mbd spectral abscissa max Re(λ) climbs past zero into Hopf oil-whip instability, the whirl orbit blooms with 0.46× sub-synchronous motion, and audio develops a ~23 Hz growl.",
    camera: "seals",
    cutaway: "half",
    shader: "pbr",
    rightTab: "rotor",
    apply: (sim) => {
      sim.reset(false);
      sim.set_oil_whip_destabilize(1.95);
      sim.set_lube_oil_pressure_frac(0.42);
      sim.set_rotor_unbalance_g_mm(75_000);
      sim.refresh_campbell();
    },
  },
  {
    id: "cold_start",
    title: "Cold Start Roll-Up Through 1st Critical Speed",
    tag: "THERMOMECH",
    summary:
      "Starts cold (45 °C metal, 143 rpm turning roll, breaker open) and accelerates smoothly across the ~1,800 rpm first lateral critical speed.",
    whatToWatch:
      "Switch to IR Thermal Map and Campbell tab: watch the thinner CrMoV rotor heat faster than the thick cast HP casing, building positive differential axial expansion.",
    camera: "overview",
    cutaway: "half",
    shader: "temp",
    rightTab: "rotor",
    apply: (sim) => {
      sim.reset(true);
      sim.set_load_setpoint(0.28);
      sim.set_rotor_unbalance_g_mm(18_000);
      sim.refresh_campbell();
    },
  },
  {
    id: "wet_steam",
    title: "Wet Steam & Baumann LP Erosion Regime",
    tag: "IF97 WET",
    summary:
      "Drops boiler superheat temperature from 525 °C down to 365 °C at 125 bar so expansion dives deep into the IAPWS-IF97 Region 4 wet dome.",
    whatToWatch:
      "LP last-stage dryness falls below 84%, Baumann moisture penalty cuts IP/LP efficiency, dense wet-mist particles fill the blue exhaust diffuser, and crackle hiss rises.",
    camera: "lp_stage",
    cutaway: "half",
    shader: "wetcreep",
    rightTab: "thermo",
    apply: (sim) => {
      sim.reset(false);
      sim.set_live_steam_bar(125);
      sim.set_live_steam_temp_c(365);
      sim.set_extraction_demand(0.05);
    },
  },
  {
    id: "labyrinth_rub",
    title: "Severe Blade-Loss Unbalance & Labyrinth Rub",
    tag: "FAULT",
    summary:
      "Injects 280,000 g·mm single-blade-loss unbalance + bearing degradation so shaft excursion strikes the 175 µm Inconel balance-piston labyrinth teeth.",
    whatToWatch:
      "Watch rub sparks at the balance piston, Stribeck rub normal force & braking torque spike, fs-eproc ISO 7919 e-value surges, and vibration trip latches the ESV.",
    camera: "seals",
    cutaway: "open",
    shader: "pbr",
    rightTab: "rotor",
    apply: (sim) => {
      sim.reset(false);
      sim.set_rotor_unbalance_g_mm(280_000);
      sim.set_oil_whip_destabilize(1.8);
      sim.set_lube_oil_pressure_frac(0.50);
      sim.refresh_campbell();
    },
  },
  {
    id: "chp_swing",
    title: "High-Extraction Industrial CHP Swing",
    tag: "CHP",
    summary:
      "Opens the 12-bar controlled extraction valve to 90% at 102% admission to supply maximum district/process steam while maintaining 50 Hz sync.",
    whatToWatch:
      "Extraction flow rises sharply, LP condensing duty drops, and overall CHP combined thermal+electric utilization climbs sharply.",
    camera: "overview",
    cutaway: "half",
    shader: "pbr",
    rightTab: "thermo",
    apply: (sim) => {
      sim.reset(false);
      sim.set_load_setpoint(1.02);
      sim.set_extraction_demand(0.90);
    },
  },
  {
    id: "lube_failure",
    title: "Lube Oil Cooler Loss & Thermal Derating",
    tag: "THERMAL",
    summary:
      "Starves lube oil pressure/cooler flow to 24% with warm cooling water so thin-film viscous shear heats the Babbitt journal bearings past 110 °C.",
    whatToWatch:
      "Watch oil film thickness h_min collapse, Babbitt temperature climb in IR mode, and fs-control::ThermalDeratingPolicy progressively derate and trip the unit.",
    camera: "exhaust",
    cutaway: "half",
    shader: "temp",
    rightTab: "ledger",
    apply: (sim) => {
      sim.reset(false);
      sim.set_lube_oil_pressure_frac(0.24);
      sim.set_cw_inlet_temp_c(35.0);
      sim.refresh_campbell();
    },
  },
];
