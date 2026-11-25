//! Composition of the Siemens SST-600 steam turbine digital twin on FrankenSim's
//! convergence kernel (`fs_system::SystemBuilder` -> `fs_system::System`).
//!
//! All 12 physics component instances are connected purely through typed `Power`,
//! `Stream` (`water-if97-v1`), `Heat`, and `Signal` connections, stepped atomically via
//! `System::step(dt, &cx)`, and audited every step by `System::ledger()`.

use core::f64::consts::PI;
use std::cell::RefCell;
use std::rc::Rc;

use fs_exec::Cx;
use fs_system::{NodeId, RouteTemperature, SolveConfig, System, SystemBuilder, heat_sink};

use crate::components::{
    BearingSystemState, CondenserState, GeneratorState, GovernorControlComponent, GovernorState,
    JournalBearingComponent, PlenumState, RotorShaftComponent, RotorShaftState, StageGroupState,
    SteamControlValveComponent, SteamPlenumComponent, SurfaceCondenserComponent,
    SynchronousGeneratorComponent, TH_FRONT_BEARING, TH_LP_HOOD, TH_REAR_BEARING,
    ThermalNetworkState, TurbineCasingThermalComponent, TurbineStageGroupComponent, ValveState,
};
use crate::continuum_2d::AxisymmetricContinuum2D;
use crate::geometry::TurbineGeometryCatalog;
use crate::rotordynamics::LateralRotorDynamics;
use crate::stage_meanline::{StageMeanlineReport, TurbineDesignParams, solve_17_stage_meanline};
use crate::steam::{eval_ph, eval_pt};

/// Live operator setpoints, boundary conditions, and injectable physical faults shared
/// with the sealed `System` boundary closures.
#[derive(Debug, Clone, Copy)]
pub struct OperatorControls {
    /// Live boiler steam pressure [Pa] (nominal 110 bar = 11.0 MPa; SST-600 max 165 bar).
    pub live_steam_p_pa: f64,
    /// Live boiler steam temperature [K] (nominal 525 °C = 798.15 K; SST-600 max 565 °C).
    pub live_steam_t_k: f64,
    /// Governor speed reference [rad/s] (nominal 3,000 rpm = 314.159 rad/s).
    pub speed_ref_rad_s: f64,
    /// Load / admission setpoint fraction `[0.0, 1.10]` (nominal 0.78 -> ~50 MW).
    pub load_setpoint: f64,
    /// Controlled 12-bar CHP process extraction valve opening `[0.0, 1.0]` (nominal 0.40).
    pub extraction_demand: f64,
    /// Extraction process header downstream backpressure [Pa] (nominal 10.5 bar = 1.05 MPa).
    pub extraction_header_p_pa: f64,
    /// Condenser cooling-water inlet temperature [K] (nominal 20 °C = 293.15 K).
    pub cw_inlet_temp_k: f64,
    /// Condenser cooling-water pump flow fraction `[0.15, 1.5]` (nominal 1.0).
    pub cw_flow_frac: f64,
    /// Condenser non-condensable air in-leakage mass fraction `[0.0, 0.40]` (nominal 0.0005).
    pub condenser_air_fraction: f64,
    /// Grid tie breaker state (`true` = synchronized to 50 Hz bus; `false` = Load Rejection!).
    pub grid_breaker_closed: bool,
    /// Generator field excitation [pu] (nominal 1.05).
    pub excitation_pu: f64,
    /// Grid electrical frequency [Hz] (nominal 50.0 Hz).
    pub grid_freq_hz: f64,
    /// Lube-oil supply pressure / cooler flow fraction `[0.15, 1.3]` (nominal 1.0).
    pub lube_oil_pressure_frac: f64,
    /// Bearing oil-whip / preload degradation factor `[0.0, 2.5]` (0 = healthy tilting-pad/elliptical).
    pub oil_whip_destabilize: f64,
    /// Rotor mass unbalance [g·mm] (nominal 14,000 g·mm ISO G0.52 precision balance; blade loss up to 350,000 g·mm).
    pub rotor_unbalance_g_mm: f64,
    /// Manual emergency trip button latch.
    pub manual_trip: bool,
    /// Momentary trip reset pulse.
    pub reset_trip_pulse: bool,
    // --- Cached 1-step sampled feedback for actuator/sensor boundaries ---
    pub sampled_omega_rad_s: f64,
    pub sampled_hp_valve_lift: f64,
    pub sampled_ext_valve_lift: f64,
    pub sampled_bearing_max_t_k: f64,
    pub sampled_lp_hood_t_k: f64,
    pub sampled_max_vib_um_pp: f64,
    pub sampled_rub_torque_nm: f64,
    pub sampled_axial_thrust_n: f64,
}

impl Default for OperatorControls {
    fn default() -> Self {
        Self {
            live_steam_p_pa: 11.0e6,        // 110 bar
            live_steam_t_k: 798.15,         // 525 °C
            speed_ref_rad_s: 100.0 * PI,    // 3,000 rpm
            load_setpoint: 0.78,            // ~50 MW condensing/extraction
            extraction_demand: 0.40,        // ~12 bar CHP bleed open 40%
            extraction_header_p_pa: 1.05e6, // 10.5 bar process header
            cw_inlet_temp_k: 293.15,        // 20 °C river/tower water
            cw_flow_frac: 1.0,
            condenser_air_fraction: 0.0005,
            grid_breaker_closed: true,
            excitation_pu: 1.05,
            grid_freq_hz: 50.0,
            lube_oil_pressure_frac: 1.0,
            oil_whip_destabilize: 0.0,
            rotor_unbalance_g_mm: 14_000.0,
            manual_trip: false,
            reset_trip_pulse: false,
            sampled_omega_rad_s: 100.0 * PI,
            sampled_hp_valve_lift: 0.78,
            sampled_ext_valve_lift: 0.40,
            sampled_bearing_max_t_k: 352.0,
            sampled_lp_hood_t_k: 332.0,
            sampled_max_vib_um_pp: 15.0,
            sampled_rub_torque_nm: 0.0,
            sampled_axial_thrust_n: 38_000.0,
        }
    }
}

/// Node handles into the sealed `fs_system::System` graph for zero-copy state inspection.
#[derive(Debug, Clone, Copy)]
pub struct TurbineNodeIds {
    pub hp_valve: NodeId,
    pub hp_chest: NodeId,
    pub hp_stages: NodeId,
    pub ext_header: NodeId,
    pub ext_valve: NodeId,
    pub ip_lp_stages: NodeId,
    pub condenser: NodeId,
    pub rotor_shaft: NodeId,
    pub bearings: NodeId,
    pub generator: NodeId,
    pub thermal: NodeId,
    pub governor: NodeId,
}

/// Sealed convergence-kernel turbine plant + coupled 4-DOF lateral rotordynamics +
/// 17-stage row-by-row IAPWS-IF97 meanline + 2D axisymmetric `(r, z)` thermal/stress PDE.
pub struct AssembledTurbinePlant {
    pub system: System,
    pub nodes: TurbineNodeIds,
    pub controls: Rc<RefCell<OperatorControls>>,
    pub lateral: LateralRotorDynamics,
    pub design: TurbineDesignParams,
    pub meanline: StageMeanlineReport,
    pub continuum: AxisymmetricContinuum2D,
    pub last_newton_iters: u32,
    pub last_shortenings: u32,
}

impl AssembledTurbinePlant {
    /// Build and seal the complete 12-component Siemens SST-600 `fs_system::System`
    /// using default nominal engineering design parameters.
    pub fn build(geom: &TurbineGeometryCatalog, cold_start: bool) -> Self {
        Self::build_with_design(geom, cold_start, TurbineDesignParams::default())
    }

    /// Build and seal the complete 12-component Siemens SST-600 `fs_system::System`
    /// parameterized by live engineer-editable [`TurbineDesignParams`].
    pub fn build_with_design(
        geom: &TurbineGeometryCatalog,
        cold_start: bool,
        design: TurbineDesignParams,
    ) -> Self {
        let design = design.clamped();
        let mut ctrl = OperatorControls::default();
        if cold_start {
            ctrl.speed_ref_rad_s = 100.0 * PI;
            ctrl.load_setpoint = 0.08;
            ctrl.extraction_demand = 0.0;
            ctrl.grid_breaker_closed = false;
            ctrl.sampled_omega_rad_s = 15.0; // Turning-gear / initial steam roll (~143 rpm)
            ctrl.sampled_hp_valve_lift = 0.08;
            ctrl.sampled_ext_valve_lift = 0.0;
            ctrl.sampled_bearing_max_t_k = 315.0;
            ctrl.sampled_lp_hood_t_k = 302.0;
        }
        let controls = Rc::new(RefCell::new(ctrl));
        let lateral = LateralRotorDynamics::new(geom.rotor_props);
        let (f_front, f_rear) = lateral.static_bearing_loads_n();

        // Initial IF97 states for steam plenums
        let live_st = eval_pt(ctrl.live_steam_p_pa, ctrl.live_steam_t_k);
        let chest_st = if cold_start {
            eval_ph(0.6e6, live_st.h_j_kg - 250_000.0)
        } else {
            eval_ph(9.6e6, live_st.h_j_kg)
        };
        let ext_st = if cold_start {
            eval_ph(0.15e6, chest_st.h_j_kg - 100_000.0)
        } else {
            eval_ph(1.35e6, 3_010_000.0)
        };
        let p_cond_init = if cold_start { 12_000.0 } else { 7_500.0 };

        let omega_init = ctrl.sampled_omega_rad_s;
        // Total shaft polar inertia = exact forged SST-600 rotor + 17 NURBS blade rows (~515 kg·m²)
        // + 60 MVA 2-pole synchronous generator field rotor & exciter (~2,100 kg·m² -> H ≈ 2.2 s, textbook!).
        let j_total = geom.rotor_props.polar_inertia_kg_m2 + 2_100.0;

        let hp_valve_c = SteamControlValveComponent::new(0.0078, 0.546);
        let hp_chest_c = SteamPlenumComponent::new(0.30);
        let hp_stages_c = TurbineStageGroupComponent::hp_section_with_design(&design);
        let ext_header_c = SteamPlenumComponent::new(0.95);
        let ext_valve_c = SteamControlValveComponent::new(0.024, 0.546);
        let ip_lp_stages_c = TurbineStageGroupComponent::ip_lp_section_with_design(&design);
        let condenser_c = SurfaceCondenserComponent::new();
        let rotor_shaft_c = RotorShaftComponent::new(j_total);
        let bearings_c = JournalBearingComponent::with_design(f_front, f_rear, &design);
        let generator_c = SynchronousGeneratorComponent::new();
        let casing_mass_scale = (design.casing_wall_mm / 95.0).clamp(0.65, 1.65);
        let thermal_c = TurbineCasingThermalComponent::new(
            geom.rotor_props.total_mass_kg,
            (geom.outer_hp_ip_casing.mass_props.mass + geom.inner_carrier.mass_props.mass)
                * casing_mass_scale,
            geom.exhaust_diffuser.mass_props.mass,
        );
        let governor_c = GovernorControlComponent::new();

        let brg_init = bearings_c.initial_state(omega_init);
        let gen_init = generator_c.initial_state(omega_init, if cold_start { 0.0 } else { 48.0e6 });
        let th_init = if cold_start {
            TurbineCasingThermalComponent::cold_start_state()
        } else {
            TurbineCasingThermalComponent::warm_nominal_state()
        };
        let gov_init = GovernorControlComponent::initial_state(
            ctrl.sampled_hp_valve_lift,
            ctrl.sampled_ext_valve_lift,
        );

        let mut b = SystemBuilder::new();
        let hp_valve = b
            .add("hp_valve", hp_valve_c, ValveState::default())
            .expect("add hp_valve");
        let hp_chest = b
            .add(
                "hp_chest",
                hp_chest_c,
                PlenumState::from_steam_point(chest_st, 0.30),
            )
            .expect("add hp_chest");
        let hp_stages = b
            .add("hp_stages", hp_stages_c, StageGroupState::default())
            .expect("add hp_stages");
        let ext_header = b
            .add(
                "ext_header",
                ext_header_c,
                PlenumState::from_steam_point(ext_st, 0.95),
            )
            .expect("add ext_header");
        let ext_valve = b
            .add("ext_valve", ext_valve_c, ValveState::default())
            .expect("add ext_valve");
        let ip_lp_stages = b
            .add("ip_lp_stages", ip_lp_stages_c, StageGroupState::default())
            .expect("add ip_lp_stages");
        let condenser = b
            .add(
                "condenser",
                condenser_c,
                CondenserState::at_vacuum_pressure(p_cond_init, 85.0),
            )
            .expect("add condenser");
        let rotor_shaft = b
            .add_with_tangent(
                "rotor_shaft",
                rotor_shaft_c,
                RotorShaftState {
                    omega_rad_s: omega_init,
                    theta_rad: 0.0,
                    net_torque_nm: 0.0,
                },
            )
            .expect("add rotor_shaft");
        let bearings = b
            .add("bearings", bearings_c, brg_init)
            .expect("add bearings");
        let generator = b
            .add("generator", generator_c, gen_init)
            .expect("add generator");
        let thermal = b
            .add_with_tangent("thermal", thermal_c, th_init)
            .expect("add thermal");
        let governor = b
            .add("governor", governor_c, gov_init)
            .expect("add governor");

        // --- Port references ---
        let p = |b: &SystemBuilder, n: NodeId, id: &str| b.port(n, id).expect("valid port");
        let hp_v_in = p(&b, hp_valve, "inlet");
        let hp_v_out = p(&b, hp_valve, "outlet");
        let hp_v_lift = p(&b, hp_valve, "lift");
        let chest_st_p = p(&b, hp_chest, "steam");
        let chest_q_p = p(&b, hp_chest, "wall_heat");
        let hp_st_in = p(&b, hp_stages, "inlet");
        let hp_st_out = p(&b, hp_stages, "outlet");
        let hp_st_sh = p(&b, hp_stages, "shaft");
        let hp_st_q = p(&b, hp_stages, "casing_heat");
        let ext_st_p = p(&b, ext_header, "steam");
        let ext_q_p = p(&b, ext_header, "wall_heat");
        let ext_v_in = p(&b, ext_valve, "inlet");
        let ext_v_out = p(&b, ext_valve, "outlet");
        let ext_v_lift = p(&b, ext_valve, "lift");
        let iplp_in = p(&b, ip_lp_stages, "inlet");
        let iplp_out = p(&b, ip_lp_stages, "outlet");
        let iplp_sh = p(&b, ip_lp_stages, "shaft");
        let iplp_q = p(&b, ip_lp_stages, "casing_heat");
        let cond_in = p(&b, condenser, "steam_in");
        let cond_dr = p(&b, condenser, "condensate_out");
        let cond_q = p(&b, condenser, "cooling_heat");
        let cond_ctrl = p(&b, condenser, "controls");
        let shaft_p = p(&b, rotor_shaft, "shaft");
        let brg_sh = p(&b, bearings, "shaft");
        let brg_q = p(&b, bearings, "oil_heat");
        let brg_flt = p(&b, bearings, "faults");
        let gen_sh = p(&b, generator, "shaft");
        let gen_el = p(&b, generator, "grid_elec");
        let gen_q = p(&b, generator, "cooler_heat");
        let gen_ctrl = p(&b, generator, "controls");
        let th_hp_q = p(&b, thermal, "hp_heat");
        let th_iplp_q = p(&b, thermal, "iplp_heat");
        let th_brg_q = p(&b, thermal, "bearing_heat");
        let th_amb_q = p(&b, thermal, "ambient_out");
        let th_ctrl = p(&b, thermal, "controls");
        let gov_sense = p(&b, governor, "sense");

        // --- 1. Steam Flowpath Connections (water-if97-v1) ---
        b.common_effort(chest_st_p, &[hp_v_out, hp_st_in])
            .expect("connect hp_chest");
        b.common_effort(ext_st_p, &[hp_st_out, ext_v_in, iplp_in])
            .expect("connect ext_header");
        b.bond(cond_in, iplp_out).expect("connect condenser inlet");

        // --- 2. Common-Flow Mechanical Shaft Junction ---
        b.common_flow(shaft_p, &[hp_st_sh, iplp_sh, brg_sh, gen_sh])
            .expect("connect shaft common_flow");

        // --- 3. Thermal Bonds & External Heat Sink Routes ---
        b.bond(th_hp_q, hp_st_q).expect("bond hp casing heat");
        b.bond(th_iplp_q, iplp_q).expect("bond ip/lp casing heat");
        b.bond(th_brg_q, brg_q).expect("bond bearing oil heat");
        b.heat_route(
            th_amb_q,
            &[],
            Some(heat_sink("ambient_and_oil_cooler", 298.15).expect("sink")),
            RouteTemperature::Fixed(298.15),
        )
        .expect("route ambient heat");
        b.heat_route(
            cond_q,
            &[],
            Some(heat_sink("condenser_cooling_water", 293.15).expect("sink")),
            RouteTemperature::Fixed(293.15),
        )
        .expect("route condenser heat");
        b.heat_route(
            gen_q,
            &[],
            Some(heat_sink("generator_cooler", 308.15).expect("sink")),
            RouteTemperature::Fixed(308.15),
        )
        .expect("route generator heat");

        // --- 4. System Boundaries reading shared OperatorControls ---
        {
            let c = Rc::clone(&controls);
            b.boundary(
                hp_v_in,
                Box::new(move |_t0, _dt, vals| {
                    let st_c = c.borrow();
                    let live = eval_pt(st_c.live_steam_p_pa, st_c.live_steam_t_k);
                    vals[0] = live.p_pa;
                    vals[1] = live.h_j_kg;
                    vals[2] = 1.0;
                }),
            )
            .expect("boiler live steam boundary");
        }
        {
            let c = Rc::clone(&controls);
            b.boundary(
                hp_v_lift,
                Box::new(move |_t0, _dt, vals| {
                    vals[0] = c.borrow().sampled_hp_valve_lift;
                }),
            )
            .expect("hp valve lift boundary");
        }
        {
            let c = Rc::clone(&controls);
            b.boundary(
                ext_v_lift,
                Box::new(move |_t0, _dt, vals| {
                    vals[0] = c.borrow().sampled_ext_valve_lift;
                }),
            )
            .expect("ext valve lift boundary");
        }
        {
            let c = Rc::clone(&controls);
            b.boundary(
                ext_v_out,
                Box::new(move |_t0, _dt, vals| {
                    let p_ext = c.borrow().extraction_header_p_pa.clamp(1.0e5, 6.0e6);
                    vals[0] = p_ext;
                    vals[1] = 2_950_000.0;
                    vals[2] = 1.0;
                }),
            )
            .expect("extraction process header boundary");
        }
        b.boundary(
            cond_dr,
            Box::new(|_t0, _dt, vals| {
                // Hotwell condensate extraction pump sink
                vals[0] = 5_000.0;
                vals[1] = 165_000.0;
                vals[2] = 1.0;
            }),
        )
        .expect("condenser hotwell boundary");
        b.boundary(chest_q_p, Box::new(|_t0, _dt, vals| vals[0] = 0.0))
            .expect("chest adiabatic boundary");
        b.boundary(ext_q_p, Box::new(|_t0, _dt, vals| vals[0] = 0.0))
            .expect("ext header adiabatic boundary");
        b.boundary(gen_el, Box::new(|_t0, _dt, vals| vals[0] = 11_000.0))
            .expect("11 kV grid bus boundary");
        {
            let c = Rc::clone(&controls);
            b.boundary(
                cond_ctrl,
                Box::new(move |_t0, _dt, vals| {
                    let st_c = c.borrow();
                    vals[0] = st_c.cw_inlet_temp_k;
                    vals[1] = st_c.cw_flow_frac;
                    vals[2] = st_c.condenser_air_fraction;
                }),
            )
            .expect("condenser controls boundary");
        }
        {
            let c = Rc::clone(&controls);
            b.boundary(
                brg_flt,
                Box::new(move |_t0, _dt, vals| {
                    let st_c = c.borrow();
                    vals[0] = st_c.lube_oil_pressure_frac;
                    vals[1] = st_c.oil_whip_destabilize;
                    vals[2] = st_c.sampled_rub_torque_nm;
                    vals[3] = st_c.sampled_axial_thrust_n;
                }),
            )
            .expect("bearing faults boundary");
        }
        {
            let c = Rc::clone(&controls);
            b.boundary(
                gen_ctrl,
                Box::new(move |_t0, _dt, vals| {
                    let st_c = c.borrow();
                    vals[0] = if st_c.grid_breaker_closed { 1.0 } else { 0.0 };
                    vals[1] = st_c.excitation_pu;
                    vals[2] = st_c.grid_freq_hz;
                }),
            )
            .expect("generator controls boundary");
        }
        {
            let c = Rc::clone(&controls);
            b.boundary(
                th_ctrl,
                Box::new(move |_t0, _dt, vals| {
                    let st_c = c.borrow();
                    vals[0] = st_c.lube_oil_pressure_frac;
                    vals[1] = st_c.cw_inlet_temp_k;
                }),
            )
            .expect("thermal controls boundary");
        }
        {
            let c = Rc::clone(&controls);
            b.boundary(
                gov_sense,
                Box::new(move |_t0, _dt, vals| {
                    let st_c = c.borrow();
                    vals[0] = st_c.sampled_omega_rad_s;
                    vals[1] = st_c.speed_ref_rad_s;
                    vals[2] = st_c.load_setpoint;
                    vals[3] = st_c.extraction_demand;
                    vals[4] = st_c.sampled_bearing_max_t_k;
                    vals[5] = st_c.sampled_lp_hood_t_k;
                    vals[6] = st_c.sampled_max_vib_um_pp;
                    vals[7] = if st_c.manual_trip { 1.0 } else { 0.0 };
                    vals[8] = if st_c.reset_trip_pulse { 1.0 } else { 0.0 };
                }),
            )
            .expect("governor sensor boundary");
        }

        let mut system = b.seal().expect("seal SST-600 convergence System");
        system.set_config(SolveConfig {
            relative_tolerance: 1.0e-7,
            absolute_tolerance: 1.0e-5,
            max_iterations: 30,
            max_shortenings: 16,
            interface_energy_tolerance: 1.0e-4,
        });

        // Seed positive physical initial guesses for every looped temperature / steam / shaft wire!
        let _ = system.guess(hp_v_out, &[chest_st.p_pa, chest_st.h_j_kg, 1.0]);
        let _ = system.guess(hp_st_in, &[chest_st.p_pa, chest_st.h_j_kg, 1.0]);
        let _ = system.guess(hp_st_out, &[ext_st.p_pa, ext_st.h_j_kg, 1.0]);
        let _ = system.guess(ext_v_in, &[ext_st.p_pa, ext_st.h_j_kg, 1.0]);
        let _ = system.guess(iplp_in, &[ext_st.p_pa, ext_st.h_j_kg, 1.0]);
        let _ = system.guess(iplp_out, &[p_cond_init, 2_570_000.0, 1.0]);
        let _ = system.guess(hp_st_q, &[th_init.temps_k[1]]);
        let _ = system.guess(iplp_q, &[th_init.temps_k[3]]);
        let _ = system.guess(brg_q, &[0.5 * (th_init.temps_k[5] + th_init.temps_k[6])]);
        let _ = system.guess(hp_st_sh, &[-omega_init]);
        let _ = system.guess(iplp_sh, &[-omega_init]);
        let _ = system.guess(brg_sh, &[-omega_init]);
        let _ = system.guess(gen_sh, &[-omega_init]);

        let nodes = TurbineNodeIds {
            hp_valve,
            hp_chest,
            hp_stages,
            ext_header,
            ext_valve,
            ip_lp_stages,
            condenser,
            rotor_shaft,
            bearings,
            generator,
            thermal,
            governor,
        };

        let meanline = solve_17_stage_meanline(
            &geom.blade_rows,
            &design,
            chest_st.p_pa,
            chest_st.h_j_kg,
            ext_st.p_pa,
            p_cond_init,
            if cold_start { 4.0 } else { 52.0 },
            if cold_start { 4.0 } else { 44.0 },
            omega_init,
        );
        let mut continuum = AxisymmetricContinuum2D::new(cold_start, &design);
        continuum.step(
            0.02,
            &meanline,
            &design,
            omega_init,
            chest_st.p_pa * 1.0e-5,
            th_init.temps_k[TH_FRONT_BEARING],
            th_init.temps_k[TH_REAR_BEARING],
            0.0,
            if cold_start { 4.0 } else { 52.0 },
        );

        Self {
            system,
            nodes,
            controls,
            lateral,
            design,
            meanline,
            continuum,
            last_newton_iters: 0,
            last_shortenings: 0,
        }
    }

    /// Advance the entire coupled `fs_system::System` + 4-DOF gyroscopic `LateralRotorDynamics`
    /// + 17-stage meanline + 2D axisymmetric `(r, z)` thermal/stress PDE over `dt_s`.
    pub fn step(&mut self, dt_s: f64, geom: &TurbineGeometryCatalog, cx: &Cx<'_>) -> Result<(), String> {
        let total = dt_s.clamp(1.0e-5, 2.0);
        let n_sub = ((total / 0.020).ceil() as usize).clamp(1, 100);
        let h = total / n_sub as f64;

        for _ in 0..n_sub {
            // 1. Atomic convergence-kernel multiphysics step!
            let rep = self
                .system
                .step(h, cx)
                .map_err(|e| format!("fs_system::System::step failed: {e}"))?;
            self.last_newton_iters = rep.iterations;
            self.last_shortenings = rep.shortenings;

            // 2. Read accepted states from the sealed System
            let shaft = *self.system.state::<RotorShaftState>(self.nodes.rotor_shaft).expect("shaft");
            let hp = *self.system.state::<StageGroupState>(self.nodes.hp_stages).expect("hp");
            let iplp = *self.system.state::<StageGroupState>(self.nodes.ip_lp_stages).expect("iplp");
            let brg = *self.system.state::<BearingSystemState>(self.nodes.bearings).expect("brg");
            let th = *self.system.state::<ThermalNetworkState>(self.nodes.thermal).expect("thermal");
            let gov = *self.system.state::<GovernorState>(self.nodes.governor).expect("gov");

            // 3. Advance 4-DOF gyroscopic lateral whirl & Stribeck rub using live bearing (K, C) + Alford Kxy!
            let alford_total = hp.alford_kxy_n_per_m + iplp.alford_kxy_n_per_m;
            let unbal = self.controls.borrow().rotor_unbalance_g_mm;
            self.lateral.step(
                rep.dt_s,
                shaft.theta_rad,
                shaft.omega_rad_s,
                &brg.front,
                &brg.rear,
                alford_total,
                unbal,
                th.diff_expansion_mm,
            );

            // 4. Push accepted feedback into the sampled sensor/actuator boundary registers
            let mut c = self.controls.borrow_mut();
            c.sampled_omega_rad_s = shaft.omega_rad_s;
            c.sampled_hp_valve_lift = gov.cv_lift * gov.esv_open;
            c.sampled_ext_valve_lift = gov.ext_lift;
            c.sampled_bearing_max_t_k = th.temps_k[TH_FRONT_BEARING].max(th.temps_k[TH_REAR_BEARING]);
            c.sampled_lp_hood_t_k = th.temps_k[TH_LP_HOOD];
            c.sampled_max_vib_um_pp = self.lateral.front_vib_um_pp.max(self.lateral.rear_vib_um_pp);
            c.sampled_rub_torque_nm = self.lateral.rub_torque_nm;
            c.sampled_axial_thrust_n = hp.axial_thrust_n + iplp.axial_thrust_n;
            c.reset_trip_pulse = false;
        }

        // 5. Solve full 17-stage row-by-row meanline & step 64x24 2D axisymmetric thermal/stress PDE once per frame slice
        let chest = self.hp_chest_state();
        let ext = self.ext_header_state();
        let cond = self.condenser_state();
        let hp = self.hp_stages_state();
        let iplp = self.ip_lp_stages_state();
        let shaft = self.rotor_shaft_state();
        let th = self.thermal_state();
        let hp_v = self.hp_valve_state();

        self.meanline = solve_17_stage_meanline(
            &geom.blade_rows,
            &self.design,
            chest.p_pa,
            chest.h_j_kg,
            ext.p_pa,
            cond.p_cond_pa,
            hp.m_dot_kg_s,
            iplp.m_dot_kg_s,
            shaft.omega_rad_s,
        );
        let rub_power_w = self.lateral.rub_torque_nm * shaft.omega_rad_s.abs();
        self.continuum.step(
            total,
            &self.meanline,
            &self.design,
            shaft.omega_rad_s,
            chest.p_pa * 1.0e-5,
            th.temps_k[TH_FRONT_BEARING],
            th.temps_k[TH_REAR_BEARING],
            rub_power_w,
            hp_v.m_dot_kg_s.max(0.0),
        );

        Ok(())
    }

    // Convenience typed accepted-state accessors
    #[must_use]
    pub fn hp_valve_state(&self) -> ValveState {
        *self.system.state::<ValveState>(self.nodes.hp_valve).expect("hp_valve")
    }
    #[must_use]
    pub fn hp_chest_state(&self) -> PlenumState {
        *self.system.state::<PlenumState>(self.nodes.hp_chest).expect("hp_chest")
    }
    #[must_use]
    pub fn hp_stages_state(&self) -> StageGroupState {
        *self.system.state::<StageGroupState>(self.nodes.hp_stages).expect("hp_stages")
    }
    #[must_use]
    pub fn ext_header_state(&self) -> PlenumState {
        *self.system.state::<PlenumState>(self.nodes.ext_header).expect("ext_header")
    }
    #[must_use]
    pub fn ext_valve_state(&self) -> ValveState {
        *self.system.state::<ValveState>(self.nodes.ext_valve).expect("ext_valve")
    }
    #[must_use]
    pub fn ip_lp_stages_state(&self) -> StageGroupState {
        *self.system.state::<StageGroupState>(self.nodes.ip_lp_stages).expect("ip_lp_stages")
    }
    #[must_use]
    pub fn condenser_state(&self) -> CondenserState {
        *self.system.state::<CondenserState>(self.nodes.condenser).expect("condenser")
    }
    #[must_use]
    pub fn rotor_shaft_state(&self) -> RotorShaftState {
        *self.system.state::<RotorShaftState>(self.nodes.rotor_shaft).expect("rotor_shaft")
    }
    #[must_use]
    pub fn bearing_state(&self) -> BearingSystemState {
        *self.system.state::<BearingSystemState>(self.nodes.bearings).expect("bearings")
    }
    #[must_use]
    pub fn generator_state(&self) -> GeneratorState {
        *self.system.state::<GeneratorState>(self.nodes.generator).expect("generator")
    }
    #[must_use]
    pub fn thermal_state(&self) -> ThermalNetworkState {
        *self.system.state::<ThermalNetworkState>(self.nodes.thermal).expect("thermal")
    }
    #[must_use]
    pub fn governor_state(&self) -> GovernorState {
        *self.system.state::<GovernorState>(self.nodes.governor).expect("governor")
    }
}
