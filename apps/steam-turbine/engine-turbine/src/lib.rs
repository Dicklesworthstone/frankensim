//! FrankenSim Industrial Steam Turbine Simulation Engine
//! (native Rust + `wasm32-unknown-unknown` WebAssembly).
//!
//! Architecture:
//! - [`geometry`]: Exact 3D CAD solids built from FrankenSim L2 primitives
//!   (`fs_rep_frep::AxisymmetricChart` + `tessellate_axisymmetric` + `mass_properties` +
//!   `surface_area`, and `fs_rep_nurbs::NurbsSurface` with exact `partials` normals across
//!   7 twisted airfoil families and 34 blade/vane rows).
//! - [`materials`] & [`matdb_src`]: 21 compile-time embedded `data/matdb/seed-v1` packs,
//!   `fs_render::thermal_surface` blackbody chromaticity, and `fs_material::creep::NortonCreep`.
//! - [`steam`]: Exact IAPWS-IF97 (`fs_thermochem::water::if97`) thermodynamic inversions,
//!   canonical `water-if97-v1` `StreamBasis`, Stodola choking, and Mollier dome sampling.
//! - [`components`] & [`system_assembly`]: 12 physics component instances implementing
//!   `fs_ports::component::Component` composed strictly through `fs_system::SystemBuilder`
//!   and audited every step by `fs_system::System::ledger()`.
//! - [`rotordynamics`]: 4-DOF gyroscopic lateral whirl, Campbell diagram, oil/steam whip
//!   eigenvalues (`fs_mbd::rotordynamics`), Stribeck labyrinth rub, and `fs_eproc::BettingEProcess`.
//! - [`acoustics`]: 48 kHz physical audio synthesizer (`fs_modal`, `fs_aeroac`, `fs_psycho`).

pub mod acoustics;
pub mod components;
pub mod continuum_2d;
pub mod design_opt;
pub mod geometry;
pub mod matdb_src;
pub mod materials;
pub mod rotordynamics;
pub mod stage_meanline;
pub mod steam;
pub mod system_assembly;

use core::f64::consts::PI;
use std::sync::OnceLock;

#[cfg(target_arch = "wasm32")]
use wasm_bindgen::prelude::*;

use fs_blake3::DomainHasher;

use crate::acoustics::{TurbineAcousticParams, TurbineAcousticSynth};
use crate::components::{
    TH_FRONT_BEARING, TH_HP_CASING, TH_HP_ROTOR, TH_IP_CASING, TH_IP_ROTOR, TH_LP_HOOD,
    TH_LUBE_OIL, TH_REAR_BEARING,
};
use crate::geometry::{ExportMesh, TurbineGeometryCatalog, build_turbine_geometry};
use crate::steam::{eval_pt, sample_mollier_dome};
use crate::system_assembly::AssembledTurbinePlant;

/// Total number of `f64` fields in the flat telemetry snapshot returned by [`SteamTurbineSim::snapshot`].
pub const SNAPSHOT_LEN: usize = 81;

/// Run a closure with a deterministic, clock-free execution context safe on
/// both native targets and `wasm32-unknown-unknown` (never touches `Instant::now`).
pub fn with_clock_free_cx<R>(f: impl FnOnce(&fs_exec::Cx<'_>) -> R) -> R {
    static GATE: OnceLock<fs_exec::CancelGate> = OnceLock::new();
    static POOL: OnceLock<fs_alloc::ArenaPool> = OnceLock::new();
    let gate = GATE.get_or_init(fs_exec::CancelGate::new_clock_free);
    let pool = POOL.get_or_init(|| fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()));
    let arena = pool.arena();
    let cx = fs_exec::Cx::new(
        gate,
        &arena,
        fs_exec::StreamKey {
            seed: 0x5354_3630_305F_0001,
            kernel_id: 0x5353_5436_3030_0001,
            tile: 0,
            iteration: 0,
        },
        fs_exec::Budget::INFINITE,
        fs_exec::ExecMode::Deterministic,
    );
    f(&cx)
}

static GEOM_CATALOG: OnceLock<TurbineGeometryCatalog> = OnceLock::new();

/// Shared geometry catalog built once on first use from `AxisymmetricChart` and `NurbsSurface`.
pub fn geometry_catalog() -> &'static TurbineGeometryCatalog {
    GEOM_CATALOG.get_or_init(|| with_clock_free_cx(build_turbine_geometry))
}

/// Ordered flat list of all 24 exported FrankenSim meshes for WebGL/Three.js boot loading.
///
/// Mesh IDs:
/// - `0`: `rotor_forging` (watertight `AxisymmetricChart`, Euler χ = 2)
/// - `1`: `labyrinth_seals` (`AxisymmetricChart`)
/// - `2, 3, 4`: `inner_carrier` (`lower_half`, `upper_back`, `upper_front`)
/// - `5, 6, 7`: `outer_hp_ip_casing` (`lower_half`, `upper_back`, `upper_front`)
/// - `8, 9, 10`: `exhaust_diffuser` (`lower_half`, `upper_back`, `upper_front`)
/// - `11, 12, 13`: `rear_bearing_cone` (`lower_half`, `upper_back`, `upper_front`)
/// - `14, 15, 16`: `front_bearing_housing` (`lower_half`, `upper_back`, `upper_front`)
/// - `17..=23`: 7 `fs_rep_nurbs::NurbsSurface` 3D twisted airfoil families (Curtis, HP, IP1, IP2, LP mid, LP giant last stage, Stator vane)
fn mesh_by_id(cat: &TurbineGeometryCatalog, id: usize) -> Option<&ExportMesh> {
    match id {
        0 => Some(&cat.rotor_mesh),
        1 => Some(&cat.labyrinth_seals),
        2 => Some(&cat.inner_carrier.lower_half),
        3 => Some(&cat.inner_carrier.upper_back),
        4 => Some(&cat.inner_carrier.upper_front),
        5 => Some(&cat.outer_hp_ip_casing.lower_half),
        6 => Some(&cat.outer_hp_ip_casing.upper_back),
        7 => Some(&cat.outer_hp_ip_casing.upper_front),
        8 => Some(&cat.exhaust_diffuser.lower_half),
        9 => Some(&cat.exhaust_diffuser.upper_back),
        10 => Some(&cat.exhaust_diffuser.upper_front),
        11 => Some(&cat.rear_bearing_cone.lower_half),
        12 => Some(&cat.rear_bearing_cone.upper_back),
        13 => Some(&cat.rear_bearing_cone.upper_front),
        14 => Some(&cat.front_bearing_housing.lower_half),
        15 => Some(&cat.front_bearing_housing.upper_back),
        16 => Some(&cat.front_bearing_housing.upper_front),
        17..=23 => cat.nurbs_airfoils.get(id - 17),
        _ => None,
    }
}

/// Top-level interactive steam turbine model exported to native Rust & WebAssembly.
#[cfg_attr(target_arch = "wasm32", wasm_bindgen)]
pub struct SteamTurbineSim {
    plant: AssembledTurbinePlant,
    acoustics: TurbineAcousticSynth,
    snapshot_buf: Vec<f64>,
    audio_buf: Vec<f32>,
    mollier_dome: Vec<f64>,
    last_error: String,
}

impl Default for SteamTurbineSim {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg_attr(target_arch = "wasm32", wasm_bindgen)]
impl SteamTurbineSim {
    /// Create a warm, synchronized 50 MW steam turbine running at 3,000 rpm.
    #[cfg_attr(target_arch = "wasm32", wasm_bindgen(constructor))]
    #[must_use]
    pub fn new() -> Self {
        let cat = geometry_catalog();
        let mut sim = Self {
            plant: AssembledTurbinePlant::build(cat, false),
            acoustics: TurbineAcousticSynth::new(),
            snapshot_buf: vec![0.0; SNAPSHOT_LEN],
            audio_buf: vec![0.0; 2048],
            mollier_dome: sample_mollier_dome(44),
            last_error: String::new(),
        };
        // Take one tiny settling step + initial acoustic probe so all telemetry & SPL are live at t=0.02 s
        let _ = sim.step(0.02);
        let _ = sim.render_audio_block(512);
        sim
    }

    /// Restart the operating state: nominal warm 50 MW synchronized operation (`cold_start = false`)
    /// or a cold roll-up from turning gear (`cold_start = true`). The machine's hardware (the
    /// active design) and the microphone placement are kept; only the operating state restarts.
    pub fn reset(&mut self, cold_start: bool) {
        let cat = geometry_catalog();
        let design = self.plant.design;
        let (mic_pos, mic_gain) = (self.acoustics.mic_pos_m, self.acoustics.mic_gain);
        self.plant = AssembledTurbinePlant::build_with_design(cat, cold_start, design);
        self.acoustics = TurbineAcousticSynth::new();
        self.acoustics.set_mic_pos(mic_pos);
        self.acoustics.mic_gain = mic_gain;
        self.last_error.clear();
        let _ = self.step(0.02);
        let _ = self.render_audio_block(512);
    }

    /// Restore everything the page can change: the default design, the default microphone and
    /// warm nominal operation.
    pub fn reset_to_defaults(&mut self) {
        let cat = geometry_catalog();
        self.plant = AssembledTurbinePlant::build(cat, false);
        self.acoustics = TurbineAcousticSynth::new();
        self.last_error.clear();
        let _ = self.step(0.02);
        let _ = self.render_audio_block(512);
    }

    /// The protection limits the engine applies, for the page to draw and colour against:
    /// `[vib_zone_ab_um, vib_alarm_um, vib_trip_um, bearing_derate_c, bearing_alarm_c,
    ///   bearing_trip_c, overspeed_trip_rpm]`.
    #[must_use]
    pub fn limits(&self) -> Vec<f64> {
        vec![
            crate::rotordynamics::ISO_7919_ZONE_AB_UM_PP,
            crate::rotordynamics::ISO_7919_ALARM_UM_PP,
            crate::rotordynamics::ISO_7919_TRIP_UM_PP,
            crate::components::BEARING_DERATE_START_K - 273.15,
            crate::components::BEARING_ALARM_K - 273.15,
            crate::components::BEARING_TRIP_K - 273.15,
            crate::components::OVERSPEED_TRIP_RAD_S * 30.0 / PI,
        ]
    }

    /// Microphone position [m] and gain, as the engine holds them: `[x, y, z, gain]`.
    #[must_use]
    pub fn mic_settings(&self) -> Vec<f64> {
        let p = self.acoustics.mic_pos_m;
        vec![p[0], p[1], p[2], self.acoustics.mic_gain]
    }

    /// Advance the coupled `fs_system::System` + 4-DOF rotordynamics + 17-stage meanline + 2D continuum by `dt_s`.
    pub fn step(&mut self, dt_s: f64) -> bool {
        let dt = if dt_s.is_finite() && dt_s > 0.0 {
            dt_s.clamp(1.0e-5, 2.0)
        } else {
            return false;
        };
        let cat = geometry_catalog();
        let res = with_clock_free_cx(|cx| self.plant.step(dt, cat, cx));
        match res {
            Ok(()) => {
                self.last_error.clear();
                self.pack_snapshot();
                true
            }
            Err(err) => {
                self.last_error = err;
                false
            }
        }
    }

    /// Last solver error string (empty when healthy).
    #[must_use]
    pub fn last_error(&self) -> String {
        self.last_error.clone()
    }

    // -----------------------------------------------------------------------
    // Interactive Control Room & Fault Setters
    // -----------------------------------------------------------------------

    pub fn set_live_steam_bar(&mut self, p_bar: f64) {
        if p_bar.is_finite() {
            self.plant.controls.borrow_mut().live_steam_p_pa = (p_bar * 1.0e5).clamp(15.0e5, 165.0e5);
        }
    }

    pub fn set_live_steam_temp_c(&mut self, t_c: f64) {
        if t_c.is_finite() {
            self.plant.controls.borrow_mut().live_steam_t_k = (t_c + 273.15).clamp(480.0, 845.0);
        }
    }

    pub fn set_speed_ref_rpm(&mut self, rpm: f64) {
        if rpm.is_finite() {
            self.plant.controls.borrow_mut().speed_ref_rad_s = (rpm.clamp(100.0, 3_600.0)) * PI / 30.0;
        }
    }

    pub fn set_load_setpoint(&mut self, frac: f64) {
        if frac.is_finite() {
            self.plant.controls.borrow_mut().load_setpoint = frac.clamp(0.0, 1.15);
        }
    }

    pub fn set_extraction_demand(&mut self, frac: f64) {
        if frac.is_finite() {
            self.plant.controls.borrow_mut().extraction_demand = frac.clamp(0.0, 1.0);
        }
    }

    pub fn set_cw_inlet_temp_c(&mut self, t_c: f64) {
        if t_c.is_finite() {
            self.plant.controls.borrow_mut().cw_inlet_temp_k = (t_c + 273.15).clamp(278.15, 323.15);
        }
    }

    pub fn set_cw_flow_frac(&mut self, frac: f64) {
        if frac.is_finite() {
            self.plant.controls.borrow_mut().cw_flow_frac = frac.clamp(0.15, 1.5);
        }
    }

    pub fn set_condenser_air_fraction(&mut self, w_air: f64) {
        if w_air.is_finite() {
            self.plant.controls.borrow_mut().condenser_air_fraction = w_air.clamp(0.0, 0.40);
        }
    }

    pub fn set_grid_breaker(&mut self, closed: bool) {
        self.plant.controls.borrow_mut().grid_breaker_closed = closed;
    }

    pub fn set_excitation_pu(&mut self, pu: f64) {
        if pu.is_finite() {
            self.plant.controls.borrow_mut().excitation_pu = pu.clamp(0.4, 1.5);
        }
    }

    pub fn set_lube_oil_pressure_frac(&mut self, frac: f64) {
        if frac.is_finite() {
            self.plant.controls.borrow_mut().lube_oil_pressure_frac = frac.clamp(0.15, 1.35);
        }
    }

    pub fn set_oil_whip_destabilize(&mut self, level: f64) {
        if level.is_finite() {
            self.plant.controls.borrow_mut().oil_whip_destabilize = level.clamp(0.0, 2.5);
        }
    }

    pub fn set_rotor_unbalance_g_mm(&mut self, g_mm: f64) {
        if g_mm.is_finite() {
            self.plant.controls.borrow_mut().rotor_unbalance_g_mm = g_mm.clamp(0.0, 400_000.0);
        }
    }

    pub fn set_manual_trip(&mut self, tripped: bool) {
        self.plant.controls.borrow_mut().manual_trip = tripped;
    }

    pub fn reset_trips(&mut self) {
        let mut c = self.plant.controls.borrow_mut();
        c.manual_trip = false;
        c.reset_trip_pulse = true;
    }

    pub fn set_mic_pos(&mut self, x: f64, y: f64, z: f64) {
        if x.is_finite() && y.is_finite() && z.is_finite() {
            self.acoustics.set_mic_pos([x, y, z]);
        }
    }

    pub fn set_mic_gain(&mut self, gain: f64) {
        if gain.is_finite() {
            self.acoustics.mic_gain = gain.clamp(0.0, 3.0);
        }
    }

    /// Force an immediate Campbell diagram recomputation using live bearing & Alford matrices.
    pub fn refresh_campbell(&mut self) {
        let shaft = self.plant.rotor_shaft_state();
        let brg = self.plant.bearing_state();
        let hp = self.plant.hp_stages_state();
        let iplp = self.plant.ip_lp_stages_state();
        self.plant.lateral.refresh_eigenvalues_and_campbell(
            &brg.front,
            &brg.rear,
            shaft.omega_rad_s,
            hp.alford_kxy_n_per_m + iplp.alford_kxy_n_per_m,
            true,
        );
    }

    // -----------------------------------------------------------------------
    // Telemetry, Metrology, Audio & Geometry Exports
    // -----------------------------------------------------------------------

    /// Flat 81-element telemetry snapshot array (see `snapshotLayout.ts` for field mapping).
    #[must_use]
    pub fn snapshot(&self) -> Vec<f64> {
        self.snapshot_buf.clone()
    }

    /// Exact IAPWS-IF97 Region 4 saturation dome `[s_f, h_f, s_g, h_g, ...]` (`kJ/kg-K`, `kJ/kg`).
    #[must_use]
    pub fn mollier_dome(&self) -> Vec<f64> {
        self.mollier_dome.clone()
    }

    /// Live 5-station steam expansion path on the Mollier (h-s) diagram:
    /// `[s_live, h_live, s_chest, h_chest, s_ext, h_ext, s_exhaust, h_exhaust, s_exhaust_is, h_exhaust_is]`
    /// in `[kJ/(kg·K), kJ/kg]`.
    #[must_use]
    pub fn mollier_path(&self) -> Vec<f64> {
        let ctrl = self.plant.controls.borrow();
        let live = eval_pt(ctrl.live_steam_p_pa, ctrl.live_steam_t_k);
        let chest = self.plant.hp_chest_state();
        let ext = self.plant.ext_header_state();
        let iplp = self.plant.ip_lp_stages_state();
        let cond = self.plant.condenser_state();
        // Isentropic expansion reference starts strictly vertically from HP Chest (s_chest) down to P_cond
        let exh_is = crate::steam::eval_ps(cond.p_cond_pa, chest.s_j_kg_k);
        vec![
            live.s_j_kg_k * 1.0e-3,
            live.h_j_kg * 1.0e-3,
            chest.s_j_kg_k * 1.0e-3,
            chest.h_j_kg * 1.0e-3,
            ext.s_j_kg_k * 1.0e-3,
            ext.h_j_kg * 1.0e-3,
            iplp.s_out_j_kg_k * 1.0e-3,
            iplp.h_out_j_kg * 1.0e-3,
            exh_is.s_j_kg_k * 1.0e-3,
            exh_is.h_j_kg * 1.0e-3,
        ]
    }

    /// Campbell diagram data `[rpm, f1_hz, f2_hz, f3_hz, f4_hz, max_re_per_s, ...]` (24 speeds × 6).
    #[must_use]
    pub fn campbell_data(&self) -> Vec<f64> {
        self.plant.lateral.campbell_table.clone()
    }

    /// Recent 128-point shaft whirl orbit `[x0_um, y0_um, ...]` at the front journal bearing.
    #[must_use]
    pub fn orbit_history(&self) -> Vec<f32> {
        self.plant.lateral.ordered_orbit_flat()
    }

    /// Render `frames` mono 48 kHz physical audio samples in `[-1.0, 1.0]`.
    pub fn render_audio_block(&mut self, frames: usize) -> Vec<f32> {
        let n = frames.clamp(1, 8192);
        if self.audio_buf.len() != n {
            self.audio_buf.resize(n, 0.0);
        }
        let shaft = self.plant.rotor_shaft_state();
        let hp_v = self.plant.hp_valve_state();
        let iplp = self.plant.ip_lp_stages_state();
        let gen_st = self.plant.generator_state();
        let th = self.plant.thermal_state();
        let lat = &self.plant.lateral;

        let whirl_intensity =
            ((lat.spectral_abscissa_per_s + 2.5) / 8.0).clamp(0.0, 1.0)
                .max(((lat.seal_excursion_um - 65.0) / 220.0).clamp(0.0, 1.0));
        let params = TurbineAcousticParams {
            omega_rad_s: shaft.omega_rad_s,
            steam_flow_kg_s: hp_v.m_dot_kg_s.max(0.0),
            valve_dp_pa: hp_v.dp_pa,
            valve_lift: hp_v.lift,
            exhaust_flow_intensity: (iplp.m_dot_kg_s.abs() / 42.0).clamp(0.0, 1.6),
            exhaust_moisture: (1.0 - iplp.quality_out).clamp(0.0, 0.25),
            max_vib_um_pp: lat.front_vib_um_pp.max(lat.rear_vib_um_pp),
            whirl_hz: lat.dominant_whirl_hz,
            whirl_intensity,
            rub_force_n: lat.rub_normal_force_n,
            gen_power_w: if gen_st.breaker_closed { gen_st.elec_power_w } else { 0.0 },
            elec_freq_hz: self.plant.controls.borrow().grid_freq_hz,
            casing_temp_k: th.temps_k[TH_HP_CASING],
        };
        self.acoustics.render_block(&params, &mut self.audio_buf);
        self.pack_snapshot();
        self.audio_buf.clone()
    }

    /// Deterministic BLAKE3 content digest of the current accepted multiphysics state.
    #[must_use]
    pub fn state_digest_hex(&self) -> String {
        let mut h = DomainHasher::new("fs-steam-turbine.state-digest.v2");
        // The acoustic slots 76..=79 change only when audio is rendered, so they are left out:
        // the digest is of the plant state alone.
        for (i, &v) in self.snapshot_buf.iter().enumerate() {
            if (76..=79).contains(&i) {
                continue;
            }
            h.update(&v.to_bits().to_le_bytes());
        }
        for &q in &self.plant.lateral.q {
            h.update(&q.to_bits().to_le_bytes());
        }
        let digest = h.finalize();
        let mut s = String::with_capacity(64);
        for b in digest.as_bytes() {
            use core::fmt::Write;
            let _ = write!(&mut s, "{b:02x}");
        }
        s
    }

    // -----------------------------------------------------------------------
    // FrankenSim L2 Geometry Exports (AxisymmetricChart + NurbsSurface)
    // -----------------------------------------------------------------------

    /// Total number of exported base meshes (`24`: 17 revolved solids/cutaway shells + 7 NURBS airfoils).
    #[must_use]
    pub fn mesh_count(&self) -> usize {
        24
    }

    #[must_use]
    pub fn mesh_positions(&self, mesh_id: usize) -> Vec<f32> {
        mesh_by_id(geometry_catalog(), mesh_id)
            .map(|m| m.positions.clone())
            .unwrap_or_default()
    }

    #[must_use]
    pub fn mesh_normals(&self, mesh_id: usize) -> Vec<f32> {
        mesh_by_id(geometry_catalog(), mesh_id)
            .map(|m| m.normals.clone())
            .unwrap_or_default()
    }

    #[must_use]
    pub fn mesh_indices(&self, mesh_id: usize) -> Vec<u32> {
        mesh_by_id(geometry_catalog(), mesh_id)
            .map(|m| m.indices.clone())
            .unwrap_or_default()
    }

    #[must_use]
    pub fn mesh_uvs(&self, mesh_id: usize) -> Vec<f32> {
        mesh_by_id(geometry_catalog(), mesh_id)
            .map(|m| m.uvs.clone())
            .unwrap_or_default()
    }

    #[must_use]
    pub fn mesh_cap_start(&self, mesh_id: usize) -> usize {
        mesh_by_id(geometry_catalog(), mesh_id)
            .map(|m| m.cap_index_start)
            .unwrap_or(0)
    }

    /// Export all 34 blade/vane rows as flat `[stage_idx, family_id, is_rotor, z_m, r_hub_m, r_tip_m, blade_count, axial_chord_m, has_shroud, ...]`,
    /// dynamically reflecting the live `lp_last_blade_h_m` engineering design parameter!
    #[must_use]
    pub fn blade_rows_flat(&self) -> Vec<f64> {
        let cat = geometry_catalog();
        let dh_lp = self.plant.design.lp_last_blade_h_m - 0.527;
        let mut out = Vec::with_capacity(cat.blade_rows.len() * 9);
        for r in &cat.blade_rows {
            let r_tip = if r.stage_index == 16 {
                (r.r_hub_m + self.plant.design.lp_last_blade_h_m).clamp(0.80, 1.08)
            } else if r.stage_index == 15 {
                (r.r_tip_m + 0.45 * dh_lp).clamp(0.62, 0.84)
            } else {
                r.r_tip_m
            };
            out.extend_from_slice(&[
                r.stage_index as f64,
                r.family_id as f64,
                if r.is_rotor { 1.0 } else { 0.0 },
                r.z_center_m,
                r.r_hub_m,
                r_tip,
                f64::from(r.blade_count),
                r.axial_chord_m,
                if r.has_shroud { 1.0 } else { 0.0 },
            ]);
        }
        out
    }

    // -----------------------------------------------------------------------
    // 17-Stage Meanline, 2D Axisymmetric Continuum PDE, FEM Shell & L4 fs-dfo Optimizer
    // -----------------------------------------------------------------------

    /// Row-by-row 17-stage meanline profile (`17 × 16 = 272` floats):
    /// `[z, r_hub, r_tip, P_bar, T_C, x, Mach, eta_pct, MW, c_a, c_u1, c_u2, R, sigma_cf, sigma_vm, thrust_kN]`.
    #[must_use]
    pub fn stage_profile_flat(&self) -> Vec<f64> {
        self.plant.meanline.to_flat_vec()
    }

    /// Packed 2D axisymmetric `(64 × 24 × 8)` continuum field tensor for GPU false-color DataTextures:
    /// `[T_K, sigma_vM_MPa, P_bar, Mach, q_r_kWm2, q_z_kWm2, wetness, creep_log]`.
    #[must_use]
    pub fn continuum_field_flat(&self) -> Vec<f32> {
        self.plant.continuum.export_flat().to_vec()
    }

    /// Axisymmetric 2D continuum grid bounds `[nz, nr, z_min_m, z_max_m, r_max_m]`.
    #[must_use]
    pub fn continuum_grid_dims(&self) -> Vec<f64> {
        vec![
            crate::continuum_2d::GRID_NZ as f64,
            crate::continuum_2d::GRID_NR as f64,
            crate::continuum_2d::GRID_Z_MIN_M,
            crate::continuum_2d::GRID_Z_MAX_M,
            crate::continuum_2d::GRID_R_MAX_M,
        ]
    }

    /// Live engineering probe at any axisymmetric `(z_m, r_m)` coordinate on the 3D CAD mesh:
    /// returns `[temp_c, sigma_vm_mpa, p_bar, mach, q_mag_kw_m2, wetness_pct, creep_log, region_id]`.
    #[must_use]
    pub fn probe_continuum(&self, z_m: f64, r_m: f64) -> Vec<f64> {
        self.plant.continuum.probe_at(z_m, r_m).to_vec()
    }

    /// Certified FEM exhaust-hood structural shell eigenmodes from `fs-plate::shell::modes_shell`
    /// and Levine-Schwinger diffuser radiation resistance from `fs-duct`:
    /// `[f0_hz, f1_hz, f2_hz, stiffness_trace_n_m, rad_resistance, certified_from_fem]`.
    #[must_use]
    pub fn fem_shell_certificate_flat(&self) -> Vec<f64> {
        let fem = crate::acoustics::exhaust_hood_fem_spectrum();
        vec![
            fem.freqs_hz[0],
            fem.freqs_hz[1],
            fem.freqs_hz[2],
            fem.stiffness_trace_n_m,
            fem.diffuser_rad_resistance,
            if fem.certified_from_fem { 1.0 } else { 0.0 },
        ]
    }

    /// Current engineer-editable design parameters `[tip_mm, hp_react, lp_h_m, brg_c_um, brg_preload, wall_mm]`.
    #[must_use]
    pub fn design_params_flat(&self) -> Vec<f64> {
        let d = self.plant.design;
        vec![
            d.tip_clearance_mm,
            d.hp_reaction,
            d.lp_last_blade_h_m,
            d.bearing_clearance_um,
            d.bearing_preload,
            d.casing_wall_mm,
        ]
    }

    /// Apply new physical design parameters, rebuilding the coupled 12-component `fs-system` plant,
    /// 17-stage meanline, 2D axisymmetric continuum, and 4-DOF Campbell diagram while keeping live operator controls!
    pub fn apply_design_params(
        &mut self,
        tip_clearance_mm: f64,
        hp_reaction: f64,
        lp_last_blade_h_m: f64,
        bearing_clearance_um: f64,
        bearing_preload: f64,
        casing_wall_mm: f64,
    ) {
        let next_design = crate::stage_meanline::TurbineDesignParams {
            tip_clearance_mm,
            hp_reaction,
            lp_last_blade_h_m,
            bearing_clearance_um,
            bearing_preload,
            casing_wall_mm,
        }
        .clamped();
        let saved_controls = *self.plant.controls.borrow();
        let cat = geometry_catalog();
        self.plant = AssembledTurbinePlant::build_with_design(cat, false, next_design);
        {
            let mut c = self.plant.controls.borrow_mut();
            c.live_steam_p_pa = saved_controls.live_steam_p_pa;
            c.live_steam_t_k = saved_controls.live_steam_t_k;
            c.speed_ref_rad_s = saved_controls.speed_ref_rad_s;
            c.load_setpoint = saved_controls.load_setpoint;
            c.extraction_demand = saved_controls.extraction_demand;
            c.cw_inlet_temp_k = saved_controls.cw_inlet_temp_k;
            c.cw_flow_frac = saved_controls.cw_flow_frac;
            c.condenser_air_fraction = saved_controls.condenser_air_fraction;
            c.grid_breaker_closed = saved_controls.grid_breaker_closed;
            c.excitation_pu = saved_controls.excitation_pu;
            c.lube_oil_pressure_frac = saved_controls.lube_oil_pressure_frac;
            c.oil_whip_destabilize = saved_controls.oil_whip_destabilize;
            c.rotor_unbalance_g_mm = saved_controls.rotor_unbalance_g_mm;
            c.manual_trip = saved_controls.manual_trip;
        }
        // Settle cleanly over 0.4 s so live MW / efficiency KPIs reflect the new machine design immediately!
        for _ in 0..8 {
            let _ = self.step(0.05);
        }
        self.refresh_campbell();
        let _ = self.render_audio_block(512);
    }

    /// Run FrankenSim's L4 ASCENT `fs-dfo::cma::cmaes` + `fs-dfo::nelder_mead` multi-variable
    /// design optimizer over the live operating envelope!
    ///
    /// `objective_id`: `0` = Max Net Electric Efficiency, `1` = Max CHP Utilization, `2` = Max Stability & Creep Life.
    #[must_use]
    pub fn run_design_optimization(&self, objective_id: u32, max_evals: u32) -> Vec<f64> {
        let cat = geometry_catalog();
        let ctrl = self.plant.controls.borrow();
        let chest = self.plant.hp_chest_state();
        let ext = self.plant.ext_header_state();
        let cond = self.plant.condenser_state();
        let hp = self.plant.hp_stages_state();
        let iplp = self.plant.ip_lp_stages_state();
        let shaft = self.plant.rotor_shaft_state();

        let rep = crate::design_opt::optimize_turbine_design(
            cat,
            self.plant.design,
            crate::design_opt::DesignObjective::from_u32(objective_id),
            max_evals as usize,
            ctrl.live_steam_p_pa,
            ctrl.live_steam_t_k,
            chest.p_pa.max(1.5e6),
            ext.p_pa.max(0.3e6),
            cond.p_cond_pa.max(2_000.0),
            hp.m_dot_kg_s.max(20.0),
            iplp.m_dot_kg_s.max(15.0),
            shaft.omega_rad_s.max(100.0),
        );
        crate::design_opt::pack_opt_report(&rep)
    }

    /// Export FrankenSim geometry provenance & exact Green's-theorem rotor inertia metrics:
    /// `[shaft_mass_kg, blades_mass_kg, total_mass_kg, com_z_m, polar_inertia_kg_m2,
    ///   transverse_inertia_kg_m2, shaft_surface_area_m2, rotor_vertices, rotor_triangles,
    ///   rotor_hausdorff_m, rotor_euler_chi, rotor_watertight, total_base_tris, total_scene_tris]`.
    #[must_use]
    pub fn geometry_receipt_flat(&self) -> Vec<f64> {
        let cat = geometry_catalog();
        let p = cat.rotor_props;
        vec![
            p.shaft_mass_kg,
            p.blades_mass_kg,
            p.total_mass_kg,
            p.com_z_m,
            p.polar_inertia_kg_m2,
            p.transverse_inertia_kg_m2,
            p.shaft_surface_area_m2,
            p.rotor_vertex_count as f64,
            p.rotor_triangle_count as f64,
            p.rotor_hausdorff_bound_m,
            f64::from(p.rotor_euler_characteristic),
            if p.rotor_watertight { 1.0 } else { 0.0 },
            cat.total_base_triangles as f64,
            cat.total_scene_triangles as f64,
        ]
    }

    fn pack_snapshot(&mut self) {
        let ctrl = self.plant.controls.borrow();
        let shaft = self.plant.rotor_shaft_state();
        let hp_v = self.plant.hp_valve_state();
        let chest = self.plant.hp_chest_state();
        let hp = self.plant.hp_stages_state();
        let ext_h = self.plant.ext_header_state();
        let ext_v = self.plant.ext_valve_state();
        let iplp = self.plant.ip_lp_stages_state();
        let cond = self.plant.condenser_state();
        let brg = self.plant.bearing_state();
        let gen_st = self.plant.generator_state();
        let th = self.plant.thermal_state();
        let gov = self.plant.governor_state();
        let lat = &self.plant.lateral;
        let ledger = self.plant.system.ledger();

        let live = eval_pt(ctrl.live_steam_p_pa, ctrl.live_steam_t_k);
        let rpm = shaft.omega_rad_s * 30.0 / PI;
        let total_shaft_w = hp.shaft_power_w + iplp.shaft_power_w;
        // Net steam thermal input (live steam minus exact IAPWS-IF97 saturated hotwell condensate h_f(P_cond))
        let h_f_cond = crate::steam::sat_pair_at_p(cond.p_cond_pa.clamp(1_500.0, 4.0e5))
            .liquid()
            .specific_enthalpy()
            .value();
        // Efficiencies are this instant's outputs over the live steam's heat input. With
        // (almost) no live steam, after a trip or at standstill, the ratio is undefined: it is
        // reported as NaN (the page shows "—") rather than clamped into a plausible number.
        let q_boiler_in_w = hp_v.m_dot_kg_s.max(0.0) * (live.h_j_kg - h_f_cond);
        let live_steam_flowing = hp_v.m_dot_kg_s > 2.5;
        let ext_heat_w = (ext_v.m_dot_kg_s.max(0.0) * (ext_h.h_j_kg - h_f_cond)).max(0.0);
        let (elec_eff_pct, chp_eff_pct) = if live_steam_flowing {
            (
                gen_st.elec_power_w / q_boiler_in_w * 100.0,
                (gen_st.elec_power_w + ext_heat_w) / q_boiler_in_w * 100.0,
            )
        } else {
            (f64::NAN, f64::NAN)
        };

        // Status code: 0 = SYNC_NOMINAL, 1 = ISLAND_OR_RUNUP, 2 = ALARM_VIB_OR_TEMP,
        //              3 = TRIP_OVERSPEED, 4 = TRIP_THERMAL, 5 = TRIP_VIBRATION, 6 = TRIP_MANUAL,
        //              7 = TRIP_POLE_SLIP (generator out of step, off the bus until the breaker is opened)
        let max_vib = lat.front_vib_um_pp.max(lat.rear_vib_um_pp);
        let max_brg_t = th.temps_k[TH_FRONT_BEARING].max(th.temps_k[TH_REAR_BEARING]);
        let vib_alarm = max_vib > crate::rotordynamics::ISO_7919_ALARM_UM_PP;
        let brg_alarm = max_brg_t > crate::components::BEARING_ALARM_K;
        let whirl_unstable = lat.spectral_abscissa_per_s > 0.0;
        // Every alarm and trip as its own bit, so the page can say which one fired.
        let flag_bits = [
            vib_alarm,
            brg_alarm,
            whirl_unstable,
            gov.overspeed_tripped,
            gov.thermal_state.tripped,
            gov.vibration_tripped,
            gov.manual_trip_flag(),
            gen_st.out_of_step_tripped,
            gov.derating_fraction < 1.0,
        ];
        let alarm_flags = flag_bits
            .iter()
            .enumerate()
            .fold(0_u32, |acc, (i, &on)| if on { acc | (1 << i) } else { acc });
        let status_code = if gov.overspeed_tripped {
            3.0
        } else if gov.thermal_state.tripped {
            4.0
        } else if gov.vibration_tripped {
            5.0
        } else if gov.manual_trip_flag() {
            6.0
        } else if gen_st.out_of_step_tripped {
            7.0
        } else if vib_alarm || brg_alarm || whirl_unstable {
            2.0
        } else if !gen_st.breaker_closed || rpm < 2_850.0 {
            1.0
        } else {
            0.0
        };

        let s = &mut self.snapshot_buf;
        s[0] = self.plant.system.time_s();
        s[1] = self.plant.system.steps() as f64;
        s[2] = rpm;
        s[3] = shaft.theta_rad;
        s[4] = shaft.omega_rad_s;
        s[5] = gen_st.elec_power_w * 1.0e-6;         // Electrical power [MW]
        s[6] = total_shaft_w * 1.0e-6;               // Shaft mechanical power [MW]
        s[7] = elec_eff_pct;                         // Net electrical Rankine efficiency [%]
        s[8] = chp_eff_pct;                          // CHP combined heat + power utilization [%]
        s[9] = status_code;

        // Steam path pressures [bar], temperatures [°C], flows [kg/s], qualities
        s[10] = ctrl.live_steam_p_pa * 1.0e-5;
        s[11] = ctrl.live_steam_t_k - 273.15;
        s[12] = hp_v.m_dot_kg_s;
        s[13] = chest.p_pa * 1.0e-5;
        s[14] = chest.t_k - 273.15;
        s[15] = ext_h.p_pa * 1.0e-5;
        s[16] = ext_h.t_k - 273.15;
        s[17] = ext_v.m_dot_kg_s;                    // Controlled process extraction flow [kg/s]
        s[18] = iplp.m_dot_kg_s;                     // LP exhaust steam flow [kg/s]
        s[19] = cond.p_cond_pa * 1.0e-3;             // Condenser pressure [kPa]
        s[20] = cond.t_sat_k - 273.15;               // Condenser saturation temp [°C]
        s[21] = iplp.quality_out;                    // LP exhaust dryness fraction x [0..1]
        s[22] = hp.eta_is * 100.0;                   // HP section isentropic efficiency [%]
        s[23] = iplp.eta_is * 100.0;                 // IP/LP section isentropic efficiency [%]
        s[24] = hp.shaft_power_w * 1.0e-6;           // HP power [MW]
        s[25] = iplp.shaft_power_w * 1.0e-6;         // IP/LP power [MW]
        s[26] = (hp.axial_thrust_n + iplp.axial_thrust_n) * 1.0e-3; // Net axial thrust [kN]
        s[27] = (hp.windage_w + iplp.windage_w) * 1.0e-3;           // Windage heating [kW]

        // Condenser heat transfer (fs_multiphase + fs_convection)
        s[28] = cond.q_cond_w * 1.0e-6;              // Condenser duty [MW]
        s[29] = cond.h_cond_w_m2_k;                  // Shell-side Nusselt/Colburn-Hougen HTC [W/(m²·K)]
        s[30] = cond.h_cw_w_m2_k;                    // Tube-side Gnielinski CW HTC [W/(m²·K)]
        s[31] = cond.u_overall_w_m2_k;               // Overall condenser U [W/(m²·K)]
        s[32] = cond.t_cw_out_k - 273.15;            // CW outlet temp [°C]
        s[33] = cond.air_mass_fraction * 100.0;      // Non-condensable air [%]

        // Governor & Valve Positions
        s[34] = gov.cv_lift;
        s[35] = gov.esv_open;
        s[36] = gov.ext_lift;
        s[37] = gov.derating_fraction;

        // Synchronous Generator (fs-phs + M19 + IACS Cu)
        s[38] = gen_st.load_angle_rad.to_degrees();  // Power angle δ [deg]
        s[39] = gen_st.line_current_a;               // Stator RMS current [A]
        s[40] = gen_st.reactive_power_var * 1.0e-6;  // Reactive power [MVAr]
        s[41] = gen_st.total_loss_w * 1.0e-3;        // Generator Cu + M19 loss [kW]
        s[42] = if gen_st.breaker_closed { 1.0 } else { 0.0 };

        // Rotordynamics (fs-mbd + fs-tribo + fs-eproc)
        s[43] = lat.front_vib_um_pp;
        s[44] = lat.rear_vib_um_pp;
        s[45] = lat.seal_excursion_um;
        s[46] = lat.rub_normal_force_n * 1.0e-3;     // Labyrinth rub force [kN]
        s[47] = lat.spectral_abscissa_per_s;         // Stability margin max Re(λ) [1/s]
        s[48] = lat.dominant_whirl_hz;               // Least-damped whirl mode [Hz]
        s[49] = lat.critical_speeds_rpm[0];          // 1st critical speed [rpm]
        s[50] = lat.critical_speeds_rpm[1];          // 2nd critical speed [rpm]
        s[51] = (hp.alford_kxy_n_per_m + iplp.alford_kxy_n_per_m) * 1.0e-6; // Alford K_xy [MN/m]
        s[52] = lat.iso_log_e_value;                 // fs-eproc ISO 7919 log e-value
        s[53] = lat.q[0] * 1.0e6;                    // Instantaneous rotor center x [µm]
        s[54] = lat.q[1] * 1.0e6;                    // Instantaneous rotor center y [µm]

        // Bearings (fs-tribo + Shell Oil S2)
        s[55] = brg.front.min_film_thickness_m * 1.0e6; // Front min oil film h_min [µm]
        s[56] = brg.rear.min_film_thickness_m * 1.0e6;  // Rear min oil film h_min [µm]
        s[57] = brg.front.eccentricity;
        s[58] = brg.front.viscosity_pa_s * 1.0e3;       // Oil dynamic viscosity [mPa·s]
        s[59] = brg.total_heat_w * 1.0e-3;              // Total bearing + rub heat [kW]

        // 8-Node Thermal Network & Norton Creep (fs-convection + fs-material)
        s[60] = th.temps_k[TH_HP_ROTOR] - 273.15;
        s[61] = th.temps_k[TH_HP_CASING] - 273.15;
        s[62] = th.temps_k[TH_IP_ROTOR] - 273.15;
        s[63] = th.temps_k[TH_IP_CASING] - 273.15;
        s[64] = th.temps_k[TH_LP_HOOD] - 273.15;
        s[65] = th.temps_k[TH_FRONT_BEARING] - 273.15;
        s[66] = th.temps_k[TH_REAR_BEARING] - 273.15;
        s[67] = th.temps_k[TH_LUBE_OIL] - 273.15;
        s[68] = th.diff_expansion_mm;
        s[69] = th.hp_thermal_stress_mpa;
        s[70] = th.hp_creep_rate_per_s;

        // Convergence Kernel fs-system::Ledger Conservation Audit
        s[71] = ledger.energy_residual_j();
        s[72] = ledger.mass_residual_kg().unwrap_or(0.0);
        s[73] = ledger.entropy_generated_j_per_k();
        s[74] = ledger.boundary_energy_in_j() * 1.0e-6; // Cumulative boundary energy [MJ]
        s[75] = self.plant.last_newton_iters as f64;

        // 48 kHz Acoustics & Psychoacoustics (fs-modal + fs-aeroac + fs-psycho)
        s[76] = self.acoustics.last_spl_db;
        s[77] = self.acoustics.last_loudness_sone;
        s[78] = self.acoustics.last_sharpness_acum;
        s[79] = self.acoustics.last_bpf_lp_hz;
        s[80] = f64::from(alarm_flags);
    }
}

trait ManualTripHelper {
    fn manual_trip_flag(&self) -> bool;
}
impl ManualTripHelper for crate::components::GovernorState {
    fn manual_trip_flag(&self) -> bool {
        self.manual_tripped
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn nominal_50mw_steady_state_closes_first_and_second_laws_in_fs_system_ledger() {
        let mut sim = SteamTurbineSim::new();
        for step_i in 0..40 {
            let s0 = sim.snapshot();
            if !sim.step(0.05) {
                panic!(
                    "step {step_i} failed: {} | prev: rpm={:.1} MW={:.2} p_ch={:.2} p_ex={:.2} p_c={:.2} t_fb={:.2} t_rb={:.2} t_oil={:.2} q_brg={:.1}kW min_h={:.2}um",
                    sim.last_error(),
                    s0[2], s0[5], s0[13], s0[15], s0[19], s0[65], s0[66], s0[67], s0[59], s0[55]
                );
            }
        }
        let s = sim.snapshot();
        let rpm = s[2];
        let elec_mw = s[5];
        let eff_pct = s[7];
        let chp_pct = s[8];
        let cond_kpa = s[19];
        let exh_x = s[21];
        let energy_res_j = s[71];
        let mass_res_kg = s[72];
        let entropy_gen = s[73];
        let boundary_mj = s[74];

        // Emergent physical steady state (NOT scripted inputs!):
        assert!((rpm - 3_000.0).abs() < 35.0, "rpm={rpm}");
        assert!(elec_mw > 25.0 && elec_mw < 80.0, "elec_mw={elec_mw}");
        assert!(eff_pct > 25.0 && eff_pct < 48.0, "eff_pct={eff_pct}");
        assert!(chp_pct > eff_pct, "CHP efficiency must exceed pure electric");
        assert!(cond_kpa > 4.0 && cond_kpa < 20.0, "cond_kpa={cond_kpa}");
        assert!(exh_x > 0.80 && exh_x < 0.99, "LP exhaust quality x={exh_x}");

        // Exact fs-system::Ledger First Law, Mass Conservation, and Second Law!
        let rel_energy_err = energy_res_j.abs() / (boundary_mj.abs() * 1.0e6 + 1.0);
        assert!(
            rel_energy_err < 1.0e-5,
            "First-Law ledger residual too large: {energy_res_j} J vs {boundary_mj} MJ (rel={rel_energy_err})"
        );
        assert!(
            mass_res_kg.abs() < 1.0e-5,
            "Stream mass ledger residual too large: {mass_res_kg} kg"
        );
        assert!(
            entropy_gen > 0.0,
            "Second-Law entropy generation must be strictly positive, got {entropy_gen}"
        );
    }

    /// Wet exhaust at a fixed inflow: the shell must settle where its duty equals
    /// `m_in (h_in - h_f)` and the condensate equals the inflow, and its internal energy must
    /// stay exactly `M u_g(p)` (the old model condensed `q / h_fg` and lost the wet fraction).
    #[test]
    fn condenser_wet_exhaust_settles_on_the_enthalpy_balance() {
        use crate::components::{CondenserState, SurfaceCondenserComponent};
        use fs_ports::component::{Component, Proposal};
        let cond = SurfaceCondenserComponent::new();
        let mut st = CondenserState::at_vacuum_pressure(7_500.0, cond.vapor_volume_m3);
        let m_in = 38.0;
        let mut h_in = 0.0;
        for _ in 0..3_000 {
            // Inlet steam of quality 0.90 at the current shell pressure.
            let pair = crate::steam::sat_pair_at_p(st.p_cond_pa);
            let h_f = pair.liquid().specific_enthalpy().value();
            let h_g = pair.vapor().specific_enthalpy().value();
            h_in = h_f + 0.90 * (h_g - h_f);
            let inputs = [
                m_in,
                m_in * h_in,
                m_in,
                5_000.0,
                165_000.0,
                1.0,
                293.15,
                293.15,
                1.0,
                0.0005,
            ];
            let mut outputs = [0.0; 7];
            let mut next = st;
            let prop = with_clock_free_cx(|cx| {
                cond.propose(&st, 0.02, &inputs, &mut next, &mut outputs, cx)
            })
            .expect("condenser proposal");
            assert!(
                matches!(prop, Proposal::Step(_)),
                "condenser refused a 20 ms step at {} Pa",
                st.p_cond_pa
            );
            st = next;
        }
        let pair = crate::steam::sat_pair_at_p(st.p_cond_pa);
        let h_f = pair.liquid().specific_enthalpy().value();
        let u_g = pair.vapor().specific_internal_energy().value();
        let duty_expected = m_in * (h_in - h_f);
        assert!(
            (st.q_cond_w - duty_expected).abs() < 2.0e-3 * duty_expected,
            "steady duty {:.4} MW vs m_in (h_in - h_f) = {:.4} MW",
            st.q_cond_w * 1.0e-6,
            duty_expected * 1.0e-6
        );
        assert!(
            (st.m_condensate_kg_s - m_in).abs() < 2.0e-3 * m_in,
            "condensate {} vs inflow {m_in}",
            st.m_condensate_kg_s
        );
        assert!(
            (st.total_u_j - st.vapor_mass_kg * u_g).abs() < 1.0e-6 * st.total_u_j.abs(),
            "U = {} must equal M u_g = {}",
            st.total_u_j,
            st.vapor_mass_kg * u_g
        );
        assert!(
            st.s_gen_w_k > 0.0,
            "condensing to colder water must produce entropy: {}",
            st.s_gen_w_k
        );
    }

    /// Entropy production of every valve, stage section and the condenser is computed, not
    /// clamped, and must not go negative through a load rejection, an ESV trip or air ingress.
    #[test]
    fn component_entropy_production_never_negative_through_transients() {
        fn check(sim: &SteamTurbineSim, label: &str) {
            let p = &sim.plant;
            let parts = [
                ("hp_valve", p.hp_valve_state().s_gen_w_k),
                ("ext_valve", p.ext_valve_state().s_gen_w_k),
                ("hp_stages", p.hp_stages_state().s_gen_w_k),
                ("ip_lp_stages", p.ip_lp_stages_state().s_gen_w_k),
                ("condenser", p.condenser_state().s_gen_w_k),
            ];
            for (name, s) in parts {
                assert!(
                    s.is_finite() && s > -1.0e-3,
                    "{label}: {name} entropy production {s} W/K"
                );
            }
        }
        let mut sim = SteamTurbineSim::new();
        for _ in 0..20 {
            assert!(sim.step(0.05));
            check(&sim, "nominal");
        }
        sim.set_grid_breaker(false);
        for _ in 0..40 {
            assert!(sim.step(0.05), "{}", sim.last_error());
            check(&sim, "load rejection");
        }
        let mut sim = SteamTurbineSim::new();
        sim.set_manual_trip(true);
        for _ in 0..40 {
            assert!(sim.step(0.05), "{}", sim.last_error());
            check(&sim, "ESV trip");
        }
        let mut sim = SteamTurbineSim::new();
        sim.set_condenser_air_fraction(0.18);
        for _ in 0..40 {
            assert!(sim.step(0.10), "{}", sim.last_error());
            check(&sim, "air ingress");
        }
    }

    /// The plant clock must advance by exactly the requested time, however the convergence
    /// kernel splits it.
    #[test]
    fn plant_clock_advances_by_the_requested_time() {
        let mut sim = SteamTurbineSim::new();
        let t0 = sim.plant.system.time_s();
        for _ in 0..30 {
            assert!(sim.step(0.07));
        }
        let advanced = sim.plant.system.time_s() - t0;
        assert!(
            (advanced - 30.0 * 0.07).abs() < 1.0e-9,
            "advanced {advanced} s for 2.1 s requested"
        );
    }

    /// The 2-D metal field must advance by the whole plant step: ten 1 s steps and a hundred
    /// 0.1 s steps of a cold start must warm the metal by the same amount (the field once
    /// clamped each call to 0.5 s, so 1 s steps lost half the heating).
    #[test]
    fn continuum_field_covers_the_whole_plant_step() {
        fn mean_metal_k(sim: &SteamTurbineSim) -> f64 {
            let c = &sim.plant.continuum;
            let (mut sum, mut n) = (0.0, 0.0);
            for (t, r) in c.temp_k.iter().zip(c.region.iter()) {
                if *r != crate::continuum_2d::CellRegion::SteamOrAir {
                    sum += t;
                    n += 1.0;
                }
            }
            sum / n
        }
        let mut coarse = SteamTurbineSim::new();
        coarse.reset(true);
        let mut fine = SteamTurbineSim::new();
        fine.reset(true);
        let t0 = mean_metal_k(&coarse);
        for _ in 0..10 {
            assert!(coarse.step(1.0), "{}", coarse.last_error());
        }
        for _ in 0..100 {
            assert!(fine.step(0.1), "{}", fine.last_error());
        }
        let (dc, df) = (mean_metal_k(&coarse) - t0, mean_metal_k(&fine) - t0);
        assert!(
            df.abs() > 1.0e-3,
            "the cold-start field must warm measurably in 10 s: {df} K"
        );
        assert!(
            (dc - df).abs() < 0.05 * df.abs(),
            "1 s steps warmed {dc} K, 0.1 s steps {df} K"
        );
    }

    /// The two stage sections take their efficiency from the row-by-row meanline, so their
    /// shaft power must match the meanline's stage powers less labyrinth leakage (and, for the
    /// IP/LP section, the unrecovered leaving loss). The sections once used their own typed
    /// efficiencies, and the HP section ran 9% below its own stages.
    #[test]
    fn stage_section_power_agrees_with_the_row_by_row_meanline() {
        let mut sim = SteamTurbineSim::new();
        for _ in 0..20 {
            assert!(sim.step(0.05), "{}", sim.last_error());
        }
        let p = &sim.plant;
        let (d, ml) = (p.design, &p.meanline);
        let (hp, ip) = (p.hp_stages_state(), p.ip_lp_stages_state());
        let ml_hp: f64 = ml.stages[0..6].iter().map(|s| s.power_w).sum();
        let ml_ip: f64 =
            ml.stages[6..17].iter().map(|s| s.power_w).sum::<f64>() - ml.lp_leaving_loss_kw * 1.0e3;
        let hp_expect = ml_hp * (1.0 - d.hp_leakage_fraction());
        let ip_expect = ml_ip * (1.0 - d.iplp_leakage_fraction());
        eprintln!(
            "HP {:.3} vs {:.3} MW, IP/LP {:.3} vs {:.3} MW, eta {:.4} {:.4}",
            hp.shaft_power_w * 1e-6,
            hp_expect * 1e-6,
            ip.shaft_power_w * 1e-6,
            ip_expect * 1e-6,
            hp.eta_is,
            ip.eta_is
        );
        assert!(
            (hp.shaft_power_w - hp_expect).abs() < 0.01 * hp_expect,
            "HP section {} W vs meanline {} W",
            hp.shaft_power_w,
            hp_expect
        );
        assert!(
            (ip.shaft_power_w - ip_expect).abs() < 0.015 * ip_expect,
            "IP/LP section {} W vs meanline {} W",
            ip.shaft_power_w,
            ip_expect
        );
    }

    /// A longer last-stage blade enlarges the exit annulus, slows the leaving flow and cuts the
    /// leaving loss; the plant must turn that into IP/LP shaft power at the same valve positions.
    /// The blade length reaches the plant only through the meanline's section efficiency.
    #[test]
    fn longer_last_stage_blades_raise_plant_power_through_the_meanline() {
        let run = |blade_h_m: f64| {
            let mut sim = SteamTurbineSim::new();
            let d = sim.plant.design;
            sim.apply_design_params(
                d.tip_clearance_mm,
                d.hp_reaction,
                blade_h_m,
                d.bearing_clearance_um,
                d.bearing_preload,
                d.casing_wall_mm,
            );
            for _ in 0..40 {
                assert!(sim.step(0.05), "{}", sim.last_error());
            }
            (
                sim.plant.ip_lp_stages_state().shaft_power_w,
                sim.plant.meanline.lp_leaving_loss_kw,
            )
        };
        let (p_short, leave_short) = run(0.46);
        let (p_long, leave_long) = run(0.60);
        eprintln!(
            "0.46 m: {:.3} MW, leaving {:.0} kW; 0.60 m: {:.3} MW, leaving {:.0} kW",
            p_short * 1e-6,
            leave_short,
            p_long * 1e-6,
            leave_long
        );
        assert!(
            leave_long < leave_short,
            "leaving loss {leave_long} kW vs {leave_short} kW"
        );
        assert!(
            p_long > p_short + 0.5 * (leave_short - leave_long) * 1.0e3,
            "IP/LP power {p_long} W vs {p_short} W"
        );
    }

    /// Excitation moves reactive power, not real power: from the machine equations the generator
    /// must absorb VArs under-excited and export them over-excited, its reactive power must match
    /// the classical `(E V cos(delta) - V^2) / X`, and its line current must carry `|S|`.
    #[test]
    fn generator_reactive_power_and_current_follow_the_machine_equations() {
        let gen_model = crate::components::SynchronousGeneratorComponent::new();
        let x = gen_model.reactance_ohm();
        let v = gen_model.nominal_voltage_v;
        let mut rows = Vec::new();
        for exc in [0.8, 1.3] {
            let mut sim = SteamTurbineSim::new();
            sim.set_excitation_pu(exc);
            for _ in 0..60 {
                assert!(sim.step(0.05), "{}", sim.last_error());
            }
            let g = sim.plant.generator_state();
            let q_classical = (exc * v * v * g.load_angle_rad.cos() - v * v) / x;
            let s_mag = g.elec_power_w.hypot(g.reactive_power_var);
            let s_from_i = 3.0_f64.sqrt() * v * g.line_current_a;
            eprintln!(
                "exc {exc}: P {:.2} MW, Q {:.2} MVAr (classical {:.2}), delta {:.2} deg, |S| {:.2} vs sqrt3 V I {:.2} MVA",
                g.elec_power_w * 1e-6,
                g.reactive_power_var * 1e-6,
                q_classical * 1e-6,
                g.load_angle_rad.to_degrees(),
                s_mag * 1e-6,
                s_from_i * 1e-6
            );
            assert!(
                (g.reactive_power_var - q_classical).abs() < 0.03 * q_classical.abs() + 0.5e6,
                "Q {} vs classical {}",
                g.reactive_power_var,
                q_classical
            );
            assert!(
                (s_mag - s_from_i).abs() < 0.015 * s_mag,
                "|S| {s_mag} vs sqrt3 V I {s_from_i}"
            );
            rows.push((g.elec_power_w, g.reactive_power_var, g.load_angle_rad));
        }
        let ((p_lo, q_lo, d_lo), (p_hi, q_hi, d_hi)) = (rows[0], rows[1]);
        assert!(
            q_lo < -20.0e6 && q_hi > 20.0e6,
            "Q under/over-excited: {q_lo} / {q_hi} VAr"
        );
        assert!(
            (p_lo - p_hi).abs() < 0.01 * p_hi,
            "real power must not follow excitation: {p_lo} vs {p_hi}"
        );
        assert!(
            d_lo > d_hi,
            "weaker field needs a larger load angle: {d_lo} vs {d_hi}"
        );
    }

    /// After a step in excitation the whole plant must swing against the bus near the
    /// electromechanical frequency sqrt(K_s / J), K_s = E V cos(delta) / (X omega_s), and the
    /// swing must decay. The machine once had a 230 pu damper and could not swing at all.
    /// In the plant the governor's 4% droop acts through the steam chest's lag (about 0.2 s),
    /// which above its corner frequency works as an extra spring: the plant swings about 7%
    /// faster than the bare machine, so the band is [0.97, 1.15] of sqrt(K_s / J). The bare
    /// machine is held to 1% by `generator_alone_swings_and_decays_as_its_constants_predict`.
    #[test]
    fn generator_swings_at_the_electromechanical_frequency_and_decays() {
        let mut sim = SteamTurbineSim::new();
        for _ in 0..40 {
            assert!(sim.step(0.05), "{}", sim.last_error());
        }
        let exc = 1.15;
        sim.set_excitation_pu(exc);
        let dt = 0.005;
        let mut delta = Vec::new();
        for _ in 0..600 {
            assert!(sim.step(dt), "{}", sim.last_error());
            delta.push(sim.plant.generator_state().load_angle_rad);
        }
        let peaks: Vec<usize> = (1..delta.len() - 1)
            .filter(|&k| delta[k] > delta[k - 1] && delta[k] >= delta[k + 1])
            .collect();
        assert!(peaks.len() >= 4, "expected a swing, got peaks at {peaks:?}");
        let period_s = (peaks[peaks.len() - 1] - peaks[0]) as f64 * dt / (peaks.len() - 1) as f64;
        let gen_model = crate::components::SynchronousGeneratorComponent::new();
        let v = gen_model.nominal_voltage_v;
        let delta_eq = *delta.last().unwrap();
        let k_s =
            exc * v * v * delta_eq.cos() / (gen_model.reactance_ohm() * gen_model.sync_omega_rad_s);
        let geom = crate::geometry_catalog();
        let j = geom.rotor_props.polar_inertia_kg_m2
            + crate::system_assembly::GENERATOR_ROTOR_INERTIA_KG_M2;
        let f_expected = (k_s / j).sqrt() / (2.0 * std::f64::consts::PI);
        let f_seen = 1.0 / period_s;
        eprintln!(
            "swing {f_seen:.3} Hz vs sqrt(K_s/J) {f_expected:.3} Hz (K_s {k_s:.3e} N m/rad, J {j:.0} kg m^2)"
        );
        assert!(
            f_seen > 0.97 * f_expected && f_seen < 1.15 * f_expected,
            "swing {f_seen} Hz vs {f_expected} Hz"
        );
        let amp = |k: usize| delta[k] - delta_eq;
        assert!(
            amp(peaks[3]) < 0.7 * amp(peaks[0]),
            "swing must decay: {} -> {}",
            amp(peaks[0]),
            amp(peaks[3])
        );
    }

    /// The generator alone on a rigid rotor of the plant's inertia, driven by a constant torque
    /// that balances it at 50 MW: kicked 3 degrees off its operating point, the load angle must
    /// ring at the damped frequency sqrt(K_s / J - sigma^2), K_s = E V cos(delta_0) / (X omega_s),
    /// and decay at sigma = D / (2 J), D being the amortisseur coefficient.
    #[test]
    fn generator_alone_swings_and_decays_as_its_constants_predict() {
        use crate::components::{GeneratorState, SynchronousGeneratorComponent};
        use fs_ports::component::Component;
        let g = SynchronousGeneratorComponent::new();
        let (ws, v, exc, dt) = (g.sync_omega_rad_s, g.nominal_voltage_v, 1.05, 1.0e-3);
        let j = geometry_catalog().rotor_props.polar_inertia_kg_m2
            + crate::system_assembly::GENERATOR_ROTOR_INERTIA_KG_M2;
        let step = |st: &GeneratorState, omega: f64| {
            let inputs = [-omega, v, 308.15, 1.0, exc, 50.0];
            let mut next = *st;
            let mut out = [0.0; 3];
            with_clock_free_cx(|cx| g.propose(st, dt, &inputs, &mut next, &mut out, cx))
                .expect("generator step");
            (next, out[0])
        };
        let st0 = g.initial_state(ws, 50.0e6, exc);
        let delta0 = st0.load_angle_rad;
        let tau_mech = -step(&st0, ws).1;
        let mut st = st0;
        st.load_angle_rad += 3.0_f64.to_radians();
        let mut omega = ws;
        let mut dev = Vec::new();
        for _ in 0..4_000 {
            let (next, tau_gen) = step(&st, omega);
            omega += dt * (tau_mech + tau_gen) / j;
            st = next;
            dev.push(st.load_angle_rad - delta0);
        }
        let peaks: Vec<usize> = (1..dev.len() - 1)
            .filter(|&k| dev[k] > dev[k - 1] && dev[k] >= dev[k + 1])
            .collect();
        assert!(peaks.len() >= 5, "expected a swing, got peaks at {peaks:?}");
        let n = peaks.len() - 1;
        let period_s = (peaks[n] - peaks[0]) as f64 * dt / n as f64;
        let sigma_seen = (dev[peaks[0]] / dev[peaks[n]]).ln() / (n as f64 * period_s);

        let k_s = exc * v * v * delta0.cos() / (g.reactance_ohm() * ws);
        let sigma = g.amortisseur_damping_nm_s_per_rad / (2.0 * j);
        let f_expected = (k_s / j - sigma * sigma).sqrt() / (2.0 * std::f64::consts::PI);
        let f_seen = 1.0 / period_s;
        eprintln!(
            "bare machine: {f_seen:.4} Hz vs {f_expected:.4} Hz, decay {sigma_seen:.4} vs {sigma:.4} 1/s"
        );
        assert!(
            (f_seen - f_expected).abs() < 0.01 * f_expected,
            "swing {f_seen} Hz vs {f_expected} Hz"
        );
        assert!(
            (sigma_seen - sigma).abs() < 0.10 * sigma,
            "decay {sigma_seen} 1/s vs {sigma} 1/s"
        );
    }

    /// A pole slip (|delta| > pi) trips the machine off the bus, it stays off while the breaker
    /// stays closed, and opening the breaker clears the latch.
    #[test]
    fn generator_pole_slip_trips_until_the_breaker_is_opened() {
        use crate::components::SynchronousGeneratorComponent;
        use fs_ports::component::Component;
        let g = SynchronousGeneratorComponent::new();
        let ws = g.sync_omega_rad_s;
        let mut st = g.initial_state(ws, 50.0e6, 1.05);
        st.load_angle_rad = 3.0;
        let run = |st: &crate::components::GeneratorState, omega: f64, breaker: f64| {
            let inputs = [-omega, 11_000.0, 308.15, breaker, 1.05, 50.0];
            let mut next = *st;
            let mut out = [0.0; 3];
            with_clock_free_cx(|cx| g.propose(st, 0.02, &inputs, &mut next, &mut out, cx))
                .expect("generator step");
            next
        };
        let st = run(&st, ws + 20.0, 1.0); // delta 3.0 + 0.4 rad: past pi
        // Current flows during the slipping step itself; the trip takes effect at its end.
        assert!(
            st.out_of_step_tripped && !st.breaker_closed,
            "slip must trip: {st:?}"
        );
        let st = run(&st, ws, 1.0);
        assert!(
            st.out_of_step_tripped && st.line_current_a == 0.0,
            "trip must latch: {st:?}"
        );
        let st = run(&st, ws, 0.0);
        assert!(
            !st.out_of_step_tripped,
            "opening the breaker must clear the trip: {st:?}"
        );
        let st = run(&st, ws, 1.0);
        assert!(
            st.breaker_closed && !st.out_of_step_tripped,
            "re-closing must reconnect: {st:?}"
        );
    }

    #[test]
    fn load_rejection_triggers_speed_surge_and_governor_catch() {
        let mut sim = SteamTurbineSim::new();
        for _ in 0..20 {
            assert!(sim.step(0.05));
        }
        let rpm_before = sim.snapshot()[2];
        // Open the 11 kV grid tie breaker at full 50 MW load!
        sim.set_grid_breaker(false);
        let mut peak_rpm = rpm_before;
        for _ in 0..80 {
            assert!(sim.step(0.05), "step failed during load rejection: {}", sim.last_error());
            peak_rpm = peak_rpm.max(sim.snapshot()[2]);
        }
        let s_after = sim.snapshot();
        assert!(
            peak_rpm > rpm_before + 80.0,
            "Load rejection must accelerate rotor from steam stored in chest/stages: peak={peak_rpm}, before={rpm_before}"
        );
        assert!(
            s_after[5].abs() < 0.2,
            "Open breaker must drop exported electrical MW to zero, got {}",
            s_after[5]
        );
        assert!(
            s_after[34] < 0.35 || s_after[35] < 0.5,
            "Governor CV or ESV must slam closed to catch overspeed: cv={}, esv={}",
            s_after[34],
            s_after[35]
        );
    }

    #[test]
    fn condenser_air_ingress_degrades_htc_raises_vacuum_pressure_and_cuts_power() {
        let mut sim = SteamTurbineSim::new();
        // MW is averaged over the last 0.5 s / 1 s: the generator swings against the bus (2.5 Hz,
        // damping ratio ~0.12) as the warm start settles, so one instant is not the level.
        let mut base_mw = 0.0;
        for k in 0..30 {
            assert!(sim.step(0.05));
            if k >= 20 {
                base_mw += sim.snapshot()[5] / 10.0;
            }
        }
        let mut base = sim.snapshot();
        base[5] = base_mw;
        // Inject 18% non-condensable air into the condenser shell (Colburn-Hougen blanketing!)
        sim.set_condenser_air_fraction(0.18);
        let mut degraded_mw = 0.0;
        for k in 0..120 {
            assert!(
                sim.step(0.10),
                "air ingress step failed: {}",
                sim.last_error()
            );
            if k >= 110 {
                degraded_mw += sim.snapshot()[5] / 10.0;
            }
        }
        let mut degraded = sim.snapshot();
        degraded[5] = degraded_mw;
        // Thresholds restated 2026-11-25 with the condenser energy fix (frankensim-89y9b.1): the
        // shell now rejects m (h_in - h_f), about 7% less than the old q = m h_fg, so the same
        // air raises the pressure less (8.8 -> 13.7 kPa) and costs 1.46 MW, not more.
        assert!(
            degraded[29] < base[29] * 0.6,
            "Colburn-Hougen condensation HTC must plunge under air ingress: before={}, after={}",
            base[29],
            degraded[29]
        );
        assert!(
            degraded[19] > base[19] + 3.0,
            "Condenser backpressure must rise under air ingress: before={} kPa, after={} kPa",
            base[19],
            degraded[19]
        );
        // Restated again with frankensim-89y9b.6/.9: the meanline's leaving loss falls as the exhaust
        // densifies and offsets part of the lost isentropic drop (see the node test), so 18% air
        // costs about 1.0 MW (51.34 -> 50.34 MW, window means).
        assert!(
            degraded[5] < base[5] - 0.8,
            "Lost condenser vacuum must reduce emergent electrical MW: before={}, after={}",
            base[5],
            degraded[5]
        );
    }

    #[test]
    fn oil_whip_and_unbalance_destabilize_4dof_rotor_eigenvalues_and_grow_vibration() {
        let mut sim = SteamTurbineSim::new();
        for _ in 0..20 {
            assert!(sim.step(0.05));
        }
        let healthy = sim.snapshot();
        // Inject oil-pressure loss + bearing preload loss + blade deposit unbalance
        sim.set_lube_oil_pressure_frac(0.35);
        sim.set_oil_whip_destabilize(2.2);
        sim.set_rotor_unbalance_g_mm(120_000.0);
        for _ in 0..60 {
            assert!(sim.step(0.05));
        }
        let faulty = sim.snapshot();
        assert!(
            faulty[47] > healthy[47],
            "4-DOF spectral abscissa max Re(lambda) must move toward instability: healthy={}, faulty={}",
            healthy[47],
            faulty[47]
        );
        assert!(
            faulty[43].max(faulty[44]) > healthy[43].max(healthy[44]) + 30.0,
            "Shaft orbit vibration must grow strongly under oil-whip/unbalance fault: healthy={}, faulty={}",
            healthy[43].max(healthy[44]),
            faulty[43].max(faulty[44])
        );
    }

    #[test]
    fn manual_esv_trip_near_zero_flow_coasts_cleanly_and_closes_ledger() {
        let mut sim = SteamTurbineSim::new();
        for _ in 0..15 {
            assert!(sim.step(0.05));
        }
        let s0 = sim.snapshot();
        assert!(
            s0[7] > 20.0 && s0[7] < 45.0 && s0[8] > s0[7] && s0[8] < 100.0,
            "nominal electric / CHP efficiency out of band: {} / {}",
            s0[7],
            s0[8]
        );
        sim.set_manual_trip(true);
        sim.set_grid_breaker(false);
        let mut rpm_at_1s = 0.0_f64;
        for step_i in 0..80 {
            assert!(
                sim.step(0.05),
                "ESV trip coastdown step {step_i} failed: {}",
                sim.last_error()
            );
            if step_i == 20 {
                rpm_at_1s = sim.snapshot()[2];
            }
        }
        let s = sim.snapshot();
        // Speed must stay below 110% overspeed (3,300 rpm) and coast smoothly downward
        // (neither runaway from oversized steam chests nor instantaneous crash from a windage sign bug!)
        assert!(
            s[2] > 2_400.0 && s[2] < 3_250.0,
            "4s ESV coastdown RPM out of physical band: rpm_4s={}, rpm_1s={}",
            s[2],
            rpm_at_1s
        );
        assert!(
            s[2] < rpm_at_1s,
            "Rotor must coast downward once trapped steam clears: rpm_1s={}, rpm_4s={}",
            rpm_at_1s,
            s[2]
        );
        assert!(s[35] < 0.05, "ESV must be closed after trip, got {}", s[35]);
        assert!(s[12] < 0.5, "Live steam flow must collapse after ESV trip, got {}", s[12]);
        // With no live steam the efficiency ratio is undefined; it used to show the 65 % / 96 %
        // display clamps.
        assert!(s[7].is_nan() && s[8].is_nan(), "efficiency without live steam: {} / {}", s[7], s[8]);
        let rel_energy_err = s[71].abs() / (s[74].abs() * 1.0e6 + 1.0);
        assert!(
            rel_energy_err < 1.0e-5,
            "First-Law ledger must stay closed through near-zero-flow ESV trip: rel={rel_energy_err}"
        );
    }

    #[test]
    fn blake3_determinism_and_48khz_audio_no_clipping() {
        let mut a = SteamTurbineSim::new();
        let mut b = SteamTurbineSim::new();
        for _ in 0..25 {
            assert!(a.step(0.04));
            assert!(b.step(0.04));
        }
        assert_eq!(a.state_digest_hex(), b.state_digest_hex());
        assert_eq!(a.state_digest_hex().len(), 64);

        let pcm = a.render_audio_block(2048);
        assert_eq!(pcm.len(), 2048);
        let peak = pcm.iter().map(|v| v.abs()).fold(0.0_f32, f32::max);
        assert!(peak > 0.01 && peak <= 0.99, "unexpected audio peak: {peak}");
        assert!(pcm.iter().all(|v| v.is_finite()));
        let s = a.snapshot();
        assert!(s[76] > 40.0 && s[76] < 120.0, "SPL dB out of range: {}", s[76]);
    }

    #[test]
    fn fem_shell_modes_2d_continuum_17_stage_and_design_opt_all_work_end_to_end() {
        let mut sim = SteamTurbineSim::new();
        for _ in 0..10 {
            assert!(sim.step(0.05));
        }

        // 1. FEM shell certification via fs-plate::shell::modes_shell + fs-duct
        let fem = sim.fem_shell_certificate_flat();
        assert_eq!(fem.len(), 6);
        assert!(fem[0] > 40.0 && fem[1] > fem[0] && fem[2] > fem[1], "FEM shell modes: {fem:?}");
        assert_eq!(fem[5], 1.0, "exhaust hood shell modes must be certified by fs-plate::shell::modes_shell");

        // 2. 17-stage meanline export & 2D axisymmetric continuum PDE export + probe
        let stages = sim.stage_profile_flat();
        assert_eq!(stages.len(), 17 * 16);
        assert!(stages[3] > stages[16 * 16 + 3], "pressure must drop monotonically from stage 0 to 16");

        let field = sim.continuum_field_flat();
        assert_eq!(field.len(), 64 * 24 * 8);
        assert!(field.iter().all(|v| v.is_finite()));
        let hp_probe = sim.probe_continuum(1.45, 0.22);
        let lp_probe = sim.probe_continuum(4.40, 0.90);
        assert!(hp_probe[0] > lp_probe[0] + 150.0, "HP rotor must be much hotter than LP hood: {hp_probe:?} vs {lp_probe:?}");

        // 3. Parametric design rebuild & L4 fs-dfo optimization — AND verify applying the optimum
        // genuinely raises the LIVE 12-component fs-system electrical MW & net electric efficiency!
        sim.apply_design_params(1.45, 0.25, 0.44, 270.0, 0.10, 130.0);
        for _ in 0..20 {
            assert!(sim.step(0.05));
        }
        let live_before = sim.snapshot();
        let opt = sim.run_design_optimization(0, 96);
        assert!(opt.len() >= 28);
        assert!(opt[16] > opt[6], "fs-dfo optimizer must improve electric efficiency: {} -> {}", opt[6], opt[16]);
        // Apply optimized geometry directly into the live fs-system convergence plant!
        sim.apply_design_params(opt[0], opt[1], opt[2], opt[3], opt[4], opt[5]);
        for _ in 0..20 {
            assert!(sim.step(0.05));
        }
        let live_after = sim.snapshot();
        assert!(
            live_after[5] > live_before[5] + 0.5 && live_after[7] > live_before[7] + 0.3,
            "Applying fs-dfo optimum must boost LIVE fs-system plant output & efficiency: {:.2} MW ({:.2}%) -> {:.2} MW ({:.2}%)",
            live_before[5],
            live_before[7],
            live_after[5],
            live_after[7]
        );
    }

    /// Orbit samples are taken inside the 0.5 ms rotor sub-steps, so one revolution at
    /// 3,000 rpm is drawn by about 40 points that sweep the full circle, and the vibration level
    /// is the peak-to-peak of a revolution. One sample per 20 ms plant step (one revolution)
    /// used to draw a stroboscopic dot or a slow false precession, and the level was 2|r| at an
    /// instant.
    #[test]
    fn orbit_resolves_each_revolution_and_vibration_is_its_peak_to_peak() {
        let mut sim = SteamTurbineSim::new();
        for _ in 0..10 {
            assert!(sim.step(0.02), "{}", sim.last_error());
        }
        let orbit = sim.orbit_history();
        let n = orbit.len() / 2;
        assert_eq!(n, crate::rotordynamics::ORBIT_RING_LEN);
        let omega = sim.plant.rotor_shaft_state().omega_rad_s;
        let per_rev = ((2.0 * PI / omega) / 0.0005).round() as usize;
        let pts: Vec<(f64, f64)> = (n - per_rev..n)
            .map(|i| (f64::from(orbit[2 * i]), f64::from(orbit[2 * i + 1])))
            .collect();
        let cx = pts.iter().map(|p| p.0).sum::<f64>() / pts.len() as f64;
        let cy = pts.iter().map(|p| p.1).sum::<f64>() / pts.len() as f64;
        let mut swept = 0.0;
        for w in pts.windows(2) {
            let a0 = (w[0].1 - cy).atan2(w[0].0 - cx);
            let a1 = (w[1].1 - cy).atan2(w[1].0 - cx);
            swept += (a1 - a0 + PI).rem_euclid(2.0 * PI) - PI;
        }
        let expect = 2.0 * PI * (per_rev - 1) as f64 / per_rev as f64;
        let (lo_x, hi_x) = pts
            .iter()
            .fold((f64::MAX, f64::MIN), |(l, h), p| (l.min(p.0), h.max(p.0)));
        let (lo_y, hi_y) = pts
            .iter()
            .fold((f64::MAX, f64::MIN), |(l, h), p| (l.min(p.1), h.max(p.1)));
        let pp = (hi_x - lo_x).max(hi_y - lo_y);
        let front = sim.snapshot()[43];
        eprintln!(
            "{per_rev} samples/rev, swept {:.3} rad (expect {:.3}), pk-pk from ring {pp:.2} um, reported {front:.2} um",
            swept, expect
        );
        assert!(
            (swept.abs() - expect).abs() < 0.25 * 2.0 * PI,
            "one revolution of orbit samples must sweep the circle: {swept} rad vs {expect} rad"
        );
        assert!(
            (front - pp).abs() <= 0.05 * pp + 0.5,
            "reported front S(p-p) {front} um must be the revolution's peak-to-peak {pp} um"
        );
    }

    /// Every alarm and trip has its own bit in snapshot slot 80, and the bits agree with the
    /// values they are judged on. The status code alone used one value for three alarms.
    #[test]
    fn alarm_flags_name_each_alarm_and_trip() {
        let check = |sim: &SteamTurbineSim, seen: &mut u32| {
            let s = sim.snapshot();
            let flags = s[80] as u32;
            *seen |= flags;
            let max_vib = s[43].max(s[44]);
            let status = s[9] as u32;
            assert_eq!(
                flags & 1 != 0,
                max_vib > crate::rotordynamics::ISO_7919_ALARM_UM_PP,
                "vibration bit vs {max_vib} um"
            );
            assert_eq!(
                flags & 4 != 0,
                s[47] > 0.0,
                "whirl bit vs Re(lambda) {}",
                s[47]
            );
            match status {
                2 => assert!(
                    flags & 0b111 != 0 && flags & 0b1111_1000 == 0,
                    "alarm status, flags {flags:#b}"
                ),
                5 => assert!(
                    flags & (1 << 5) != 0,
                    "vibration trip status, flags {flags:#b}"
                ),
                0 | 1 => assert!(flags & 0b1111_1111 == 0, "healthy status, flags {flags:#b}"),
                _ => {}
            }
        };
        let mut seen = 0_u32;
        let mut sim = SteamTurbineSim::new();
        for _ in 0..10 {
            assert!(sim.step(0.05), "{}", sim.last_error());
            check(&sim, &mut seen);
        }
        assert_eq!(
            seen & 0b1111_1111,
            0,
            "nominal running raises no alarm: {seen:#b}"
        );
        sim.set_rotor_unbalance_g_mm(280_000.0);
        for _ in 0..60 {
            assert!(sim.step(0.05), "{}", sim.last_error());
            check(&sim, &mut seen);
        }
        let s = sim.snapshot();
        eprintln!(
            "after heavy unbalance: vib {:.0}/{:.0} um, status {}, flags seen {seen:#b}",
            s[43], s[44], s[9]
        );
        assert!(
            seen & 1 != 0,
            "heavy unbalance must raise the vibration alarm bit: {seen:#b}"
        );
    }

    /// Restarting the operating state keeps the machine's hardware and the microphone; only
    /// "defaults" restores them. Cold Start and the experiment cards used to rebuild the
    /// default machine behind a page still showing the custom one.
    #[test]
    fn reset_keeps_the_design_and_microphone_and_defaults_restore_them() {
        let defaults = SteamTurbineSim::new();
        let mut sim = SteamTurbineSim::new();
        sim.apply_design_params(1.2, 0.40, 0.58, 260.0, 0.2, 110.0);
        sim.set_mic_pos(4.0, 1.4, 6.0);
        sim.set_mic_gain(2.0);
        let custom = sim.design_params_flat();
        assert_ne!(custom, defaults.design_params_flat());
        for cold in [true, false] {
            sim.reset(cold);
            assert_eq!(
                sim.design_params_flat(),
                custom,
                "reset(cold = {cold}) changed the design"
            );
            assert_eq!(
                sim.mic_settings(),
                vec![4.0, 1.4, 6.0, 2.0],
                "reset(cold = {cold}) moved the microphone"
            );
        }
        sim.reset_to_defaults();
        assert_eq!(sim.design_params_flat(), defaults.design_params_flat());
        assert_eq!(sim.mic_settings(), defaults.mic_settings());
        assert_eq!(
            sim.state_digest_hex(),
            defaults.state_digest_hex(),
            "defaults must be the boot state"
        );
    }

    /// Blade tip clearance costs efficiency through each row's tip loss in the meanline. The
    /// section's labyrinth bypass used to scale with the same clearance, charging it twice;
    /// the plant-to-meanline power ratio must not move with the clearance.
    #[test]
    fn blade_tip_clearance_is_charged_once() {
        let run = |tip_mm: f64| {
            let mut sim = SteamTurbineSim::new();
            let d = sim.plant.design;
            sim.apply_design_params(
                tip_mm,
                d.hp_reaction,
                d.lp_last_blade_h_m,
                d.bearing_clearance_um,
                d.bearing_preload,
                d.casing_wall_mm,
            );
            for _ in 0..20 {
                assert!(sim.step(0.05), "{}", sim.last_error());
            }
            let p = &sim.plant;
            let ml_hp: f64 = p.meanline.stages[0..6].iter().map(|s| s.power_w).sum();
            (p.hp_stages_state().shaft_power_w, ml_hp)
        };
        let (plant_a, ml_a) = run(0.75);
        let (plant_b, ml_b) = run(1.50);
        let drift = (plant_b / ml_b) / (plant_a / ml_a) - 1.0;
        eprintln!(
            "0.75 mm: plant {:.3} / meanline {:.3} MW; 1.50 mm: {:.3} / {:.3} MW; ratio drift {drift:.5}",
            plant_a * 1e-6,
            ml_a * 1e-6,
            plant_b * 1e-6,
            ml_b * 1e-6
        );
        assert!(
            ml_b < ml_a,
            "a wider tip clearance must cost meanline power: {ml_a} -> {ml_b}"
        );
        assert!(
            drift.abs() < 0.003,
            "plant/meanline power ratio moved {drift} with tip clearance"
        );
    }
}
