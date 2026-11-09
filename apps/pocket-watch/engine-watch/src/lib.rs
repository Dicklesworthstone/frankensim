//! FrankenSim Grand Complication Pocket Watch (`Calibre FS-1903 Tourbillon`)
//! coupled multiphysics & horological engine (Native Rust + WebAssembly).

pub mod acoustics;
pub mod horology;
pub mod materials;
pub mod multiphysics;

#[cfg(target_arch = "wasm32")]
use wasm_bindgen::prelude::*;

use crate::acoustics::WatchAcousticSynth;
use crate::horology::{HorologyState, MechanicalAcousticEvent};
use crate::multiphysics::{MultiphysicsState, NODE_COUNT};

/// Number of `f64` header fields before the 12-node telemetry block.
pub const SNAPSHOT_HEADER_LEN: usize = 32;

/// Number of `f64` fields per node in the snapshot buffer.
pub const SNAPSHOT_NODE_STRIDE: usize = 12;

/// Total length of the flat `f64` snapshot array (`32 + 12 * 12 = 176`).
pub const SNAPSHOT_LEN: usize = SNAPSHOT_HEADER_LEN + NODE_COUNT * SNAPSHOT_NODE_STRIDE;

/// Top-level coupled simulation engine for the FrankenSim Grand Complication Pocket Watch.
#[cfg_attr(target_arch = "wasm32", wasm_bindgen)]
pub struct PocketWatchSim {
    horology: HorologyState,
    multiphysics: MultiphysicsState,
    acoustics: WatchAcousticSynth,
    pending_events: Vec<MechanicalAcousticEvent>,
    snapshot_buf: Vec<f64>,
    audio_buf: Vec<f32>,
}

impl Default for PocketWatchSim {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg_attr(target_arch = "wasm32", wasm_bindgen)]
impl PocketWatchSim {
    /// Construct a new, chronometer-calibrated Patek-style Tourbillon pocket watch at $293.15\text{ K}$.
    #[cfg_attr(target_arch = "wasm32", wasm_bindgen(constructor))]
    #[must_use]
    pub fn new() -> Self {
        let mut sim = Self {
            horology: HorologyState::new(),
            multiphysics: MultiphysicsState::new(),
            acoustics: WatchAcousticSynth::new(),
            pending_events: Vec::with_capacity(64),
            snapshot_buf: vec![0.0; SNAPSHOT_LEN],
            audio_buf: vec![0.0; 2048],
        };
        sim.pack_snapshot();
        sim
    }

    /// Reset the entire watch and all materials back to pristine $20.0^\circ\text{C}$ chronometer state.
    pub fn reset(&mut self) {
        self.horology = HorologyState::new();
        self.multiphysics = MultiphysicsState::new();
        self.acoustics = WatchAcousticSynth::new();
        self.pending_events.clear();
        self.pack_snapshot();
    }

    /// Wind the mainspring crown by `delta_turns` (where `8.5` turns = full 42-hour wind).
    pub fn wind_crown(&mut self, delta_turns: f64) {
        self.horology
            .wind_crown(delta_turns, &mut self.pending_events);
        self.flush_events_to_acoustics();
        self.pack_snapshot();
    }

    /// Configure the virtual oxy-butane blowtorch ($T_{\text{flame}} = 2250\text{ K}$) in 3D watch coordinates ($\text{mm}$).
    pub fn set_blowtorch(&mut self, power: f64, x_mm: f64, y_mm: f64, z_mm: f64) {
        self.multiphysics.inputs.blowtorch_power = power.clamp(0.0, 1.0);
        self.multiphysics.inputs.blowtorch_pos_mm = [x_mm, y_mm, z_mm];
    }

    /// Configure the cryogenic liquid nitrogen ($\text{LN}_2$ at $77.36\text{ K}$) pour stream in 3D watch coordinates ($\text{mm}$).
    pub fn set_liquid_nitrogen(&mut self, pour_rate: f64, x_mm: f64, y_mm: f64, z_mm: f64) {
        self.multiphysics.inputs.ln2_pour_rate = pour_rate.clamp(0.0, 1.0);
        self.multiphysics.inputs.ln2_pos_mm = [x_mm, y_mm, z_mm];
    }

    /// Configure the watch oil dropper / flood stream (`pennzane-shf-x-2000-bearing-oil`) in 3D watch coordinates ($\text{mm}$).
    pub fn set_oil_pour(&mut self, pour_rate: f64, x_mm: f64, y_mm: f64, z_mm: f64) {
        self.multiphysics.inputs.oil_pour_rate = pour_rate.clamp(0.0, 1.0);
        self.multiphysics.inputs.oil_pos_mm = [x_mm, y_mm, z_mm];
    }

    /// Ultrasonic clean excess flooded oil and restore factory $1.6\text{ }\mu\text{m}$ jewel lubrication.
    pub fn clean_excess_oil(&mut self) {
        self.multiphysics.ultrasonic_clean_and_relubricate();
        self.pack_snapshot();
    }

    /// Configure the external sintered $\text{NdFeB N42}$ super-magnet ($B_r = 1.32\text{ T}$) in 3D watch coordinates ($\text{mm}$).
    pub fn set_magnet(
        &mut self,
        enabled: bool,
        strength: f64,
        x_mm: f64,
        y_mm: f64,
        z_mm: f64,
    ) {
        self.multiphysics.inputs.magnet_enabled = enabled;
        self.multiphysics.inputs.magnet_strength = strength.clamp(0.0, 2.0);
        self.multiphysics.inputs.magnet_pos_mm = [x_mm, y_mm, z_mm];
    }

    /// Fire the alternating-field degaussing coil to wipe out residual ferromagnetic remanence $M_{\text{rem}}$.
    pub fn trigger_degauss(&mut self) {
        self.multiphysics.trigger_degauss();
        self.pack_snapshot();
    }

    /// Toggle between traditional High-Carbon Blued Steel (`false`) and Compensated Invar/Fe-36Ni (`true`) hairspring.
    pub fn set_compensated_hairspring(&mut self, compensated: bool) {
        self.multiphysics.inputs.compensated_hairspring = compensated;
    }

    /// Position the 3D studio measurement microphone in watch coordinates ($\text{mm}$) and set gain.
    pub fn set_microphone(&mut self, x_mm: f64, y_mm: f64, z_mm: f64, gain: f64) {
        self.acoustics.mic_pos_mm = [x_mm, y_mm, z_mm];
        self.acoustics.mic_gain = gain.clamp(0.0, 4.0);
    }

    /// Advance the coupled multiphysics and horological simulation by `dt_s` seconds.
    pub fn step(&mut self, dt_s: f64) {
        let clamped_dt = dt_s.clamp(1.0e-4, 0.10);
        // Step in <= 10 ms slices so thermal/phase and escapement coupling stay tightly locked
        let slices = ((clamped_dt / 0.010).ceil() as usize).clamp(1, 12);
        let sub_dt = clamped_dt / (slices as f64);

        for _ in 0..slices {
            self.multiphysics.step(sub_dt);
            self.horology.step(
                sub_dt,
                &self.multiphysics.modifiers,
                &mut self.pending_events,
            );
            self.flush_events_to_acoustics();
        }
        self.pack_snapshot();
    }

    /// Render `frames` mono audio samples at $48\text{ kHz}$ from the 3D virtual microphone
    /// and return them as a `Vec<f32>` (converted to `Float32Array` across `wasm-bindgen`).
    #[must_use]
    pub fn render_audio_block(&mut self, frames: usize) -> Vec<f32> {
        let n = frames.clamp(1, 8192);
        if self.audio_buf.len() != n {
            self.audio_buf.resize(n, 0.0);
        }
        self.acoustics
            .render_samples(&self.multiphysics, &mut self.audio_buf);
        self.audio_buf.clone()
    }

    /// Return the current 176-element `f64` state snapshot for 3D rendering and telemetry.
    #[must_use]
    pub fn snapshot(&self) -> Vec<f64> {
        self.snapshot_buf.clone()
    }

    /// Return the 64x8 (512-node) 3D Thermo-Elasto-Viscoplastic + Surface-Tension Case Shell snapshot
    /// (`512 * 6 = 3072` `f32`s): `[x_mm, y_mm, z_mm, thickness_mm, temperature_k, melt_fraction]`
    /// for each 3D Lagrangian node `(it, iz)` on the Bassine 18k gold case.
    #[must_use]
    pub fn case_shell_snapshot(&self) -> Vec<f32> {
        let mut out = Vec::with_capacity(self.multiphysics.case_shell.len() * 6);
        for sn in &self.multiphysics.case_shell {
            out.push(sn.pos_mm[0] as f32);
            out.push(sn.pos_mm[1] as f32);
            out.push(sn.pos_mm[2] as f32);
            out.push(sn.thickness_mm as f32);
            out.push(sn.temperature_k as f32);
            out.push(sn.melt_fraction as f32);
        }
        out
    }

    /// Return the 24x24 (576-cell) 2D Finite-Volume Spatial PDE grid snapshot (`576 * 7 = 4032` `f32`s):
    /// `[temperature_k, melt_fraction, case_slump_frac, liquid_gold_mm, frozen_gold_mm, frost_um, oil_um]`
    /// for each cell `(ix, iy)` in row-major order.
    #[must_use]
    pub fn spatial_grid_snapshot(&self) -> Vec<f32> {
        let mut out = Vec::with_capacity(self.multiphysics.grid.len() * 7);
        for cell in &self.multiphysics.grid {
            out.push(cell.temperature_k as f32);
            out.push(cell.melt_fraction as f32);
            out.push(cell.case_slump_frac as f32);
            out.push(cell.liquid_gold_mm as f32);
            out.push(cell.frozen_gold_mm as f32);
            out.push(cell.frost_um as f32);
            out.push(cell.oil_um as f32);
        }
        out
    }

    /// Return the 18 discretized elastica mainspring coil radii ($\text{mm}$) inside the barrel.
    #[must_use]
    pub fn mainspring_coils_mm(&self) -> Vec<f32> {
        self.horology
            .mainspring_coil_radii_mm
            .iter()
            .map(|&r| r as f32)
            .collect()
    }

    /// Return the dynamic recoil click lever (*cliquet*) deflection angle ($\text{rad}$).
    #[must_use]
    pub fn click_lever_angle_rad(&self) -> f64 {
        self.horology.click_lever_angle_rad
    }

    /// Return the keyless works crown wheel rotation angle ($\text{rad}$).
    #[must_use]
    pub fn crown_wheel_angle_rad(&self) -> f64 {
        self.horology.crown_wheel_angle_rad
    }

    /// Human-readable status banner for the current horological & multiphysics state.
    #[must_use]
    pub fn status_label(&self) -> String {
        self.horology.status.label().to_string()
    }

    /// Deterministic BLAKE3 digest of the current 176-float state snapshot (`fs-blake3`).
    #[must_use]
    pub fn state_digest_hex(&self) -> String {
        let mut bytes = Vec::with_capacity(self.snapshot_buf.len() * 8);
        for &val in &self.snapshot_buf {
            let canonical = if val == 0.0 { 0.0 } else { val };
            bytes.extend_from_slice(&canonical.to_le_bytes());
        }
        fs_blake3::hash_bytes(&bytes).to_hex()
    }

    fn flush_events_to_acoustics(&mut self) {
        for ev in self.pending_events.drain(..) {
            self.acoustics.trigger_event(ev, &self.multiphysics);
        }
    }

    fn pack_snapshot(&mut self) {
        let h = &self.horology;
        let m = &self.multiphysics.modifiers;
        let b = &mut self.snapshot_buf;

        b[0] = h.lab_time_s;
        b[1] = h.watch_displayed_time_s;
        b[2] = h.mainspring_wind;
        b[3] = h.mainspring_torque_nm;
        b[4] = h.barrel_arbor_angle_rad;
        b[5] = h.barrel_wheel_angle_rad;
        b[6] = h.center_wheel_angle_rad;
        b[7] = h.third_wheel_angle_rad;
        b[8] = h.tourbillon_cage_angle_rad;
        b[9] = h.escape_wheel_visual_rad;
        b[10] = h.pallet_fork_angle_rad;
        b[11] = h.balance_angle_rad;
        b[12] = h.balance_omega_rad_s;
        b[13] = h.balance_amplitude_deg;
        b[14] = h.total_beats as f64;
        b[15] = h.measured_rate_s_per_day;
        b[16] = m.theoretical_rate_s_per_day;
        b[17] = h.measured_beat_error_ms;
        b[18] = h.moonphase_angle_rad;
        b[19] = (h.status as u8) as f64;
        b[20] = m.balance_inertia_ratio;
        b[21] = m.hairspring_stiffness_ratio;
        b[22] = m.pivot_friction_multiplier;
        b[23] = m.min_pivot_clearance_um;
        b[24] = m.b_field_escapement_mt;
        b[25] = m.remanent_magnetization_mt;
        b[26] = m.gold_fouling_jam_factor;
        b[27] = m.ln2_boiling_power_w;
        b[28] = m.blowtorch_heating_power_w;
        b[29] = m.blued_steel_rgb[0];
        b[30] = m.blued_steel_rgb[1];
        b[31] = m.blued_steel_rgb[2];

        for (idx, node) in self.multiphysics.nodes.iter().enumerate() {
            let base = SNAPSHOT_HEADER_LEN + idx * SNAPSHOT_NODE_STRIDE;
            b[base] = node.temperature_k;
            b[base + 1] = node.enthalpy_j_kg;
            b[base + 2] = node.melt_fraction;
            b[base + 3] = node.plastic_slump;
            b[base + 4] = node.molten_gold_g;
            b[base + 5] = node.frozen_gold_g;
            b[base + 6] = node.ln2_film_g;
            b[base + 7] = node.frost_thickness_um;
            b[base + 8] = node.oil_film_um;
            b[base + 9] = node.incandescence_rgb[0];
            b[base + 10] = node.incandescence_rgb[1];
            b[base + 11] = node.incandescence_rgb[2];
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn full_watch_sim_deterministic_replay_and_audio_synthesis() {
        let mut sim_a = PocketWatchSim::new();
        let mut sim_b = PocketWatchSim::new();

        for _ in 0..50 {
            sim_a.step(0.02);
            sim_b.step(0.02);
        }
        assert_eq!(sim_a.state_digest_hex(), sim_b.state_digest_hex());

        let audio = sim_a.render_audio_block(1024);
        assert_eq!(audio.len(), 1024);
        let max_sample = audio.iter().fold(0.0_f32, |acc, &x| acc.max(x.abs()));
        assert!(max_sample > 1.0e-4, "expected audible tick-tock modal response");
    }

    #[test]
    fn full_scenario_blowtorch_melt_then_ln2_freeze_stops_watch_permanently() {
        let mut sim = PocketWatchSim::new();
        // Verify watch is running initially
        for _ in 0..25 {
            sim.step(0.02);
        }
        let beats_before = sim.snapshot()[14];
        assert!(beats_before >= 2.0);

        // Apply full blowtorch to 6h gold case until gold melts and fouls the tourbillon
        sim.set_blowtorch(1.0, 0.0, -25.0, 0.0);
        for _ in 0..160 {
            sim.step(0.05);
        }
        let snap_melted = sim.snapshot();
        assert!(
            snap_melted[26] > 0.3,
            "expected gold fouling jam factor > 0.3, got {}",
            snap_melted[26]
        );

        // Now turn off blowtorch and pour Liquid Nitrogen (LN2) to quench & solidify in deformed shape
        sim.set_blowtorch(0.0, 0.0, -25.0, 0.0);
        sim.set_liquid_nitrogen(1.0, 0.0, -18.0, 2.0);
        for _ in 0..120 {
            sim.step(0.05);
        }
        assert_eq!(
            sim.status_label(),
            "LOCKED — ENCASED IN SOLIDIFIED GOLD"
        );
        let beats_locked = sim.snapshot()[14];
        for _ in 0..25 {
            sim.step(0.02);
        }
        assert_eq!(
            sim.snapshot()[14],
            beats_locked,
            "watch must not tick once locked in solidified gold"
        );
    }
}
