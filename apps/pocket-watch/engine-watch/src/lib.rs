//! FrankenSim Grand Complication Pocket Watch (`Calibre FS-1903 Tourbillon`):
//! a coupled multiphysics and horological engine (native Rust + WebAssembly).
//!
//! The pieces, in the order a step runs them:
//! - [`multiphysics`]: the heat, phase-change, liquid-film and magnetic fields
//!   on a 64 × 64 plan grid of the watch, and the case-band shell. It samples
//!   what each moving part sits in (temperature, gold, oil, frost, field) into a
//!   [`horology::MovementEnvironment`].
//! - [`horology`]: the mainspring, going train, tourbillon cage and the
//!   contact-dynamics Swiss lever [`escapement`] driving the balance. Every
//!   rate, amplitude and stoppage is what that mechanism does in that
//!   environment; contact impacts are reported as [`horology::AcousticEvent`]s.
//! - [`acoustics`]: the case and movement modes, rung by those impacts at
//!   their own lab times, heard at a virtual microphone.

pub mod acoustics;
pub mod escapement;
pub mod horology;
pub mod matdb_src;
pub mod materials;
pub mod multiphysics;

#[cfg(target_arch = "wasm32")]
use wasm_bindgen::prelude::*;

use crate::acoustics::WatchAcousticSynth;
use crate::horology::{AcousticEvent, Horology};
use crate::multiphysics::{GRID_HALF_SPAN_MM, GRID_N, MultiphysicsState, NODE_COUNT};

/// Number of `f64` header fields before the 12-node telemetry block.
pub const SNAPSHOT_HEADER_LEN: usize = 32;

/// Number of `f64` fields per node in the snapshot buffer.
pub const SNAPSHOT_NODE_STRIDE: usize = 12;

/// Total length of the flat `f64` snapshot array (`32 + 12 * 12 = 176`).
pub const SNAPSHOT_LEN: usize = SNAPSHOT_HEADER_LEN + NODE_COUNT * SNAPSHOT_NODE_STRIDE;

/// Longest coupling slice between the field solve and the movement [s].
const COUPLING_SLICE_S: f64 = 0.010;

/// Top-level coupled simulation of the pocket watch.
#[cfg_attr(target_arch = "wasm32", wasm_bindgen)]
pub struct PocketWatchSim {
    horology: Horology,
    multiphysics: MultiphysicsState,
    acoustics: WatchAcousticSynth,
    pending_events: Vec<AcousticEvent>,
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
    /// A regulated, wound tourbillon pocket watch at 20 °C.
    #[cfg_attr(target_arch = "wasm32", wasm_bindgen(constructor))]
    #[must_use]
    pub fn new() -> Self {
        let mut sim = Self {
            horology: Horology::new(),
            multiphysics: MultiphysicsState::new(),
            acoustics: WatchAcousticSynth::new(),
            pending_events: Vec::with_capacity(64),
            snapshot_buf: vec![0.0; SNAPSHOT_LEN],
            audio_buf: vec![0.0; 2048],
        };
        sim.sync_lab_time();
        sim.refresh_derived();
        sim.pack_snapshot();
        sim
    }

    /// Back to the pristine watch at 20 °C (the hairspring choice is kept).
    pub fn reset(&mut self) {
        let compensated = self.horology.compensated;
        let mut mp = MultiphysicsState::new();
        mp.inputs.compensated_hairspring = compensated;
        let mut h = Horology::new();
        h.set_compensated(compensated);
        self.horology = h;
        self.multiphysics = mp;
        self.acoustics = WatchAcousticSynth::new();
        self.pending_events.clear();
        self.sync_lab_time();
        self.refresh_derived();
        self.pack_snapshot();
    }

    /// Turn the crown by `delta_turns` (the crown-to-barrel gearing and the
    /// click are simulated; the bridle slips past full wind).
    pub fn wind_crown(&mut self, delta_turns: f64) {
        let turns = if delta_turns.is_finite() { delta_turns.clamp(-20.0, 20.0) } else { 0.0 };
        self.horology.wind_crown(turns, &self.multiphysics.environment, &mut self.pending_events);
        self.flush_events_to_acoustics();
        self.pack_snapshot();
    }

    /// Air–propane torch: power 0–1 and nozzle position in watch coordinates [mm].
    pub fn set_blowtorch(&mut self, power: f64, x_mm: f64, y_mm: f64, z_mm: f64) {
        self.multiphysics.inputs.blowtorch_power = unit(power);
        self.multiphysics.inputs.blowtorch_pos_mm = pos(x_mm, y_mm, z_mm);
    }

    /// Liquid-nitrogen pour: rate 0–1 and spout position [mm].
    pub fn set_liquid_nitrogen(&mut self, pour_rate: f64, x_mm: f64, y_mm: f64, z_mm: f64) {
        self.multiphysics.inputs.ln2_pour_rate = unit(pour_rate);
        self.multiphysics.inputs.ln2_pos_mm = pos(x_mm, y_mm, z_mm);
    }

    /// Watch-oil dropper: rate 0–1 and nozzle position [mm].
    pub fn set_oil_pour(&mut self, pour_rate: f64, x_mm: f64, y_mm: f64, z_mm: f64) {
        self.multiphysics.inputs.oil_pour_rate = unit(pour_rate);
        self.multiphysics.inputs.oil_pos_mm = pos(x_mm, y_mm, z_mm);
    }

    /// Ultrasonic clean: strip flooded oil, frost and LN₂ from the movement
    /// and re-oil the jewels.
    pub fn clean_excess_oil(&mut self) {
        self.multiphysics.ultrasonic_clean_and_relubricate();
        self.horology.relubricate();
        self.pack_snapshot();
    }

    /// N42 NdFeB magnet: enabled, strength (0–2, scales the remanence) and
    /// position [mm].
    pub fn set_magnet(&mut self, enabled: bool, strength: f64, x_mm: f64, y_mm: f64, z_mm: f64) {
        let i = &mut self.multiphysics.inputs;
        i.magnet_enabled = enabled;
        i.magnet_strength = if strength.is_finite() { strength.clamp(0.0, 2.0) } else { 0.0 };
        i.magnet_pos_mm = pos(x_mm, y_mm, z_mm);
        self.pack_snapshot();
    }

    /// Fire the 50 Hz degaussing coil (a 3 s geometrically decaying alternating field).
    pub fn trigger_degauss(&mut self) {
        self.multiphysics.trigger_degauss();
        self.pack_snapshot();
    }

    /// Hairspring alloy: carbon steel (`false`) or a Nivarox-type
    /// thermocompensating alloy (`true`). Each is regulated at the factory.
    pub fn set_compensated_hairspring(&mut self, compensated: bool) {
        self.multiphysics.inputs.compensated_hairspring = compensated;
        self.horology.set_compensated(compensated);
        self.refresh_derived();
        self.pack_snapshot();
    }

    /// Measurement microphone position [mm] and gain (0–4).
    pub fn set_microphone(&mut self, x_mm: f64, y_mm: f64, z_mm: f64, gain: f64) {
        self.acoustics.mic_pos_mm = pos(x_mm, y_mm, z_mm);
        self.acoustics.mic_gain = if gain.is_finite() { gain.clamp(0.0, 4.0) } else { 1.0 };
    }

    /// Compressed-air nozzle: power 0–1 and position [mm].
    pub fn set_air_jet(&mut self, power: f64, x_mm: f64, y_mm: f64, z_mm: f64) {
        self.multiphysics.inputs.air_jet_power = unit(power);
        self.multiphysics.inputs.air_jet_pos_mm = pos(x_mm, y_mm, z_mm);
        self.refresh_derived();
    }

    /// Bell-jar pressure [atm], 0 (vacuum) to 5.
    pub fn set_ambient_pressure_atm(&mut self, pressure_atm: f64) {
        let p = if pressure_atm.is_finite() { pressure_atm.clamp(0.0, 5.0) } else { 1.0 };
        self.multiphysics.inputs.ambient_pressure_atm = p;
        self.multiphysics.wake();
        self.refresh_derived();
        self.pack_snapshot();
    }

    /// Gravity / centrifuge acceleration in g, watch frame (x to 3 o'clock,
    /// y to 12 o'clock, z out of the dial).
    pub fn set_gravity_vector(&mut self, gx: f64, gy: f64, gz: f64) {
        let c = |v: f64| if v.is_finite() { v.clamp(-6.0, 6.0) } else { 0.0 };
        self.multiphysics.inputs.gravity_g_vec = [c(gx), c(gy), c(gz)];
        self.multiphysics.wake();
        self.refresh_derived();
        self.pack_snapshot();
    }

    /// A drop or knock: a half-sine acceleration pulse of `g_peak` (10–5000 g)
    /// applied to the movement and the liquids on it.
    pub fn trigger_mechanical_shock(&mut self, g_peak: f64) {
        let g = if g_peak.is_finite() { g_peak.clamp(10.0, 5000.0) } else { 10.0 };
        self.multiphysics.apply_mechanical_shock(g);
        self.horology.apply_mechanical_shock(g, &mut self.pending_events);
        self.flush_events_to_acoustics();
        self.pack_snapshot();
    }

    /// Advance the coupled simulation by `dt_s` seconds (clamped to 0.1 ms –
    /// 100 ms), in slices of at most 10 ms: the fields first (placing the
    /// lever and cage for the magnet), then the movement in the environment
    /// the fields leave, then the sound sources it produced.
    pub fn step(&mut self, dt_s: f64) {
        let dt = if dt_s.is_finite() { dt_s.clamp(1.0e-4, 0.10) } else { 1.0e-4 };
        let slices = (dt / COUPLING_SLICE_S).ceil().max(1.0) as usize;
        let sub = dt / slices as f64;
        for _ in 0..slices {
            let cage = self.horology.tourbillon_cage_angle_rad;
            let lever = self.horology.esc.q[1];
            self.multiphysics.step(sub, cage, lever);
            self.horology.step(sub, &self.multiphysics.environment, &mut self.pending_events);
            self.flush_events_to_acoustics();
        }
        self.sync_lab_time();
        self.pack_snapshot();
    }

    /// Render `frames` (1–8192) mono samples at 48 kHz from the microphone.
    #[must_use]
    pub fn render_audio_block(&mut self, frames: usize) -> Vec<f32> {
        let n = frames.clamp(1, 8192);
        if self.audio_buf.len() != n {
            self.audio_buf.resize(n, 0.0);
        }
        self.acoustics.render_samples(&self.multiphysics, &mut self.audio_buf);
        self.audio_buf.clone()
    }

    /// The 176-value state snapshot (32 header values, then 12 per node).
    #[must_use]
    pub fn snapshot(&self) -> Vec<f64> {
        self.snapshot_buf.clone()
    }

    /// The 16-value escapement state: θ, θ̇, φ, φ̇, ψ, ψ̇ (balance, lever and
    /// escape wheel in the cage frame), the active-contact mask (bit 0 entry
    /// stone, 1 exit stone, 2 pin/fork, 3 banking, 4 guard), the last impulse
    /// energy [µJ], escapement efficiency, measured lift [deg], balance Q,
    /// escape-wheel drive [µN m], barrel torque [mN m], power reserve [h],
    /// balance pivot friction [nN m] and the air share of the losses.
    #[must_use]
    pub fn escapement_state(&self) -> Vec<f64> {
        self.horology.escapement_state().to_vec()
    }

    /// Lock drops since the last call, as `[lab_time_s, side, energy_nJ]`
    /// triples (side 0 = entry, 1 = exit), for the timegrapher.
    pub fn timegrapher_drain(&mut self) -> Vec<f64> {
        self.horology.drain_timegrapher()
    }

    /// Cells per side of the plan grid behind `spatial_grid_snapshot`.
    #[must_use]
    pub fn grid_n(&self) -> usize {
        GRID_N
    }

    /// Half the side of the plan grid [mm].
    #[must_use]
    pub fn grid_half_span_mm(&self) -> f64 {
        GRID_HALF_SPAN_MM
    }

    /// The case-band shell, 64 × 8 nodes × `[x, y, z, thickness, T, melt]`
    /// (mm, K, 0–1).
    #[must_use]
    pub fn case_shell_snapshot(&self) -> Vec<f32> {
        let mut out = Vec::with_capacity(self.multiphysics.case_shell.len() * 6);
        for sn in &self.multiphysics.case_shell {
            out.extend_from_slice(&[
                sn.pos_mm[0] as f32,
                sn.pos_mm[1] as f32,
                sn.pos_mm[2] as f32,
                sn.thickness_mm as f32,
                sn.temperature_k as f32,
                sn.melt_fraction as f32,
            ]);
        }
        out
    }

    /// The plan grid, `grid_n()²` cells row-major from (−span, −span), each
    /// `[T K, melt 0–1, case slump 0–1, liquid gold mm, frozen gold mm,
    /// frost µm, oil µm]`.
    #[must_use]
    pub fn spatial_grid_snapshot(&self) -> Vec<f32> {
        let mut out = Vec::with_capacity(self.multiphysics.grid.len() * 7);
        for c in &self.multiphysics.grid {
            out.extend_from_slice(&[
                c.temperature_k as f32,
                c.melt_fraction as f32,
                c.case_slump_frac() as f32,
                c.liquid_gold_mm() as f32,
                c.frozen_gold_mm() as f32,
                c.frost_um() as f32,
                (c.oil_depth_m() * 1.0e6) as f32,
            ]);
        }
        out
    }

    /// Mainspring coil radii inside the barrel [mm], innermost first.
    #[must_use]
    pub fn mainspring_coils_mm(&self) -> Vec<f32> {
        self.horology.mainspring_coil_radii_mm.iter().map(|&r| r as f32).collect()
    }

    /// Click (ratchet pawl) angle [rad].
    #[must_use]
    pub fn click_lever_angle_rad(&self) -> f64 {
        self.horology.click_lever_angle_rad
    }

    /// Crown-wheel angle [rad].
    #[must_use]
    pub fn crown_wheel_angle_rad(&self) -> f64 {
        self.horology.crown_wheel_angle_rad
    }

    /// Status banner: which mechanism dominates the movement's behaviour.
    #[must_use]
    pub fn status_label(&self) -> String {
        self.horology.status.label().to_string()
    }

    /// Sound pressure level at the microphone [dB re 20 µPa].
    #[must_use]
    pub fn mic_spl_db(&self) -> f64 {
        self.acoustics.last_spl_db
    }

    /// ISO 532-1 loudness at the microphone [sone] (`fs-psycho`).
    #[must_use]
    pub fn mic_loudness_sone(&self) -> f64 {
        self.acoustics.last_loudness_sone
    }

    /// ISO 532-1 loudness level at the microphone [phon] (`fs-psycho`).
    #[must_use]
    pub fn mic_loudness_phon(&self) -> f64 {
        self.acoustics.last_loudness_phon
    }

    /// DIN 45692 sharpness at the microphone [acum] (`fs-psycho`).
    #[must_use]
    pub fn mic_sharpness_acum(&self) -> f64 {
        self.acoustics.last_sharpness_acum
    }

    /// Rate change per kelvin of the balance and hairspring [s/day/K], from
    /// the alloy's thermoelastic coefficient and the balance's expansion.
    #[must_use]
    pub fn rate_sensitivity_s_day_per_k(&self) -> f64 {
        self.horology.rate_sensitivity_s_day_per_k
    }

    /// Lower end of the ±2σ interval on the measured rate [s/day].
    #[must_use]
    pub fn rate_interval_lo_s_day(&self) -> f64 {
        self.horology.rate_interval().0
    }

    /// Upper end of the ±2σ interval on the measured rate [s/day].
    #[must_use]
    pub fn rate_interval_hi_s_day(&self) -> f64 {
        self.horology.rate_interval().1
    }

    /// Air damping of the balance as a decay rate [1/s].
    #[must_use]
    pub fn balance_aero_drag_s_inv(&self) -> f64 {
        self.horology.balance_aero_drag_s_inv
    }

    /// Anytime-valid e-value that the watch is outside the COSC −4/+6 s/day
    /// band (`fs-eproc` betting e-process).
    #[must_use]
    pub fn cosc_e_value(&self) -> f64 {
        self.horology.cosc_e_value
    }

    /// BLAKE3 digest of the snapshot (`fs-blake3`), for replay checks.
    #[must_use]
    pub fn state_digest_hex(&self) -> String {
        let mut bytes = Vec::with_capacity(self.snapshot_buf.len() * 8);
        for &val in &self.snapshot_buf {
            let canonical = if val == 0.0 { 0.0 } else { val };
            bytes.extend_from_slice(&canonical.to_le_bytes());
        }
        fs_blake3::hash_bytes(&bytes).to_hex()
    }
}

impl PocketWatchSim {
    /// The field state (for native tools and tests).
    #[must_use]
    pub fn multiphysics(&self) -> &MultiphysicsState {
        &self.multiphysics
    }

    /// The movement state (for native tools and tests).
    #[must_use]
    pub fn horology(&self) -> &Horology {
        &self.horology
    }

    /// Re-derive the environment and the horology readouts after an input
    /// changed, so a paused watch reports the new state at once.
    fn refresh_derived(&mut self) {
        self.multiphysics.refresh_environment();
        self.horology.refresh_derived(&self.multiphysics.environment);
    }

    fn sync_lab_time(&mut self) {
        self.acoustics.advance_sim_clock(self.horology.lab_time_s);
    }

    fn flush_events_to_acoustics(&mut self) {
        for ev in self.pending_events.drain(..) {
            self.acoustics.push_event(ev);
        }
    }

    fn pack_snapshot(&mut self) {
        let h = &self.horology;
        let m = &self.multiphysics;
        let b = &mut self.snapshot_buf;

        b[0] = h.lab_time_s;
        b[1] = h.watch_displayed_time_s;
        b[2] = h.wind_fraction();
        b[3] = h.mainspring_torque_nm;
        b[4] = h.barrel_arbor_angle_rad;
        b[5] = h.barrel_wheel_angle_rad;
        b[6] = h.center_wheel_angle_rad;
        b[7] = h.third_wheel_angle_rad;
        b[8] = h.tourbillon_cage_angle_rad;
        b[9] = h.esc.q[2];
        b[10] = h.esc.q[1];
        b[11] = h.esc.q[0];
        b[12] = h.esc.v[0];
        b[13] = h.balance_amplitude_deg;
        b[14] = h.total_beats as f64;
        b[15] = h.measured_rate_s_per_day;
        b[16] = h.theoretical_rate_s_per_day;
        b[17] = h.measured_beat_error_ms;
        b[18] = h.moonphase_angle_rad;
        b[19] = f64::from(h.status as u8);
        b[20] = h.inertia_ratio;
        b[21] = h.stiffness_ratio;
        b[22] = h.pivot_friction_multiplier;
        b[23] = h.min_clearance_um;
        b[24] = m.b_field_escapement_mt;
        b[25] = m.remanent_magnetization_mt;
        b[26] = h.gold_jam_factor;
        b[27] = m.ln2_boiling_power_w;
        b[28] = m.blowtorch_heating_power_w;
        b[29] = m.blued_steel_rgb[0];
        b[30] = m.blued_steel_rgb[1];
        b[31] = m.blued_steel_rgb[2];

        for (idx, node) in m.nodes.iter().enumerate() {
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

fn unit(v: f64) -> f64 {
    if v.is_finite() { v.clamp(0.0, 1.0) } else { 0.0 }
}

fn pos(x: f64, y: f64, z: f64) -> [f64; 3] {
    let c = |v: f64| if v.is_finite() { v.clamp(-500.0, 500.0) } else { 0.0 };
    [c(x), c(y), c(z)]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn deterministic_replay_and_audible_ticks() {
        let mut a = PocketWatchSim::new();
        let mut b = PocketWatchSim::new();
        let mut peak = 0.0_f32;
        for _ in 0..60 {
            a.step(0.02);
            b.step(0.02);
            let audio = a.render_audio_block(960);
            assert_eq!(audio.len(), 960);
            peak = audio.iter().fold(peak, |acc, &x| acc.max(x.abs()));
        }
        assert_eq!(a.state_digest_hex(), b.state_digest_hex());
        assert!(peak > 1.0e-4, "the escapement impacts ring the case: peak {peak}");
        assert!(a.snapshot()[14] >= 4.0, "beats counted");
        assert!(a.mic_spl_db().is_finite());
        assert!(a.mic_loudness_sone() >= 0.0);
        assert!(a.mic_sharpness_acum() >= 0.0);
        assert!(a.rate_interval_lo_s_day() <= a.rate_interval_hi_s_day());
        assert!(a.balance_aero_drag_s_inv() > 0.0);
        assert!(a.cosc_e_value() > 0.0);
        let esc = a.escapement_state();
        assert_eq!(esc.len(), 16);
        assert!(esc.iter().all(|v| v.is_finite()));
        let drained = a.timegrapher_drain();
        assert_eq!(drained.len() % 3, 0);
        assert_eq!(a.spatial_grid_snapshot().len(), a.grid_n() * a.grid_n() * 7);
        assert_eq!(a.case_shell_snapshot().len(), 512 * 6);
    }

    /// Print a time series of the coupled state (scenario calibration).
    fn trace(sim: &mut PocketWatchSim, tag: &str, secs: f64, every: f64) {
        let t_end = sim.horology().lab_time_s + secs;
        let mut next = sim.horology().lab_time_s;
        let clock = std::time::Instant::now();
        let t0 = sim.horology().lab_time_s;
        while sim.horology().lab_time_s < t_end - 1e-9 {
            sim.step(0.05);
            let h = sim.horology();
            if h.lab_time_s >= next - 1e-9 {
                next += every;
                let m = sim.multiphysics();
                let env = &m.environment;
                let case_t = m.nodes[..4].iter().map(|n| n.temperature_k).fold(0.0, f64::max);
                let melt = m.nodes[..4].iter().map(|n| n.melt_fraction).fold(0.0, f64::max);
                println!(
                    "{tag} t={:6.2} {:<44} amp={:6.1} rate={:8.2} beats={:5} jam={:.2} \
                     Tcase={:6.0} melt={:.2} Tbal={:6.1} Tesc={:6.1} goldL={:.3}mm goldS={:.3}mm oil={:.1}um frost={:.1}um \
                     drip={:.2}g B={:.1}mT P_ln2={:.1}W P_torch={:.1}W clr={:.2}um wall={:.2}x",
                    h.lab_time_s,
                    h.status.label(),
                    h.balance_amplitude_deg,
                    h.measured_rate_s_per_day,
                    h.total_beats,
                    h.gold_jam_factor,
                    case_t,
                    melt,
                    env.balance.temperature_k,
                    env.escapement.temperature_k,
                    env.escapement.liquid_gold_m * 1e3,
                    env.escapement.solid_gold_m * 1e3,
                    env.escapement.oil_m * 1e6,
                    env.escapement.frost_m * 1e6,
                    m.dripped_gold_g,
                    m.b_field_escapement_mt,
                    m.ln2_boiling_power_w,
                    m.blowtorch_heating_power_w,
                    h.min_clearance_um,
                    clock.elapsed().as_secs_f64() / (h.lab_time_s - t0).max(1e-9),
                );
            }
        }
    }

    #[test]
    #[ignore = "scenario calibration dump"]
    fn scenario_torch_hanging_6h_then_ln2() {
        let mut sim = PocketWatchSim::new();
        trace(&mut sim, "idle", 3.0, 1.0);
        sim.set_blowtorch(1.0, 0.0, -25.0, 1.5);
        trace(&mut sim, "torch6h", 60.0, 3.0);
        sim.set_blowtorch(0.0, 0.0, -25.0, 1.5);
        sim.set_liquid_nitrogen(1.0, 0.0, -18.0, 2.0);
        trace(&mut sim, "ln2", 30.0, 2.0);
    }

    #[test]
    #[ignore = "scenario calibration dump"]
    fn scenario_torch_flat_6h_then_ln2() {
        let mut sim = PocketWatchSim::new();
        sim.set_gravity_vector(0.0, 0.0, -1.0);
        sim.set_blowtorch(1.0, 0.0, -25.0, 1.5);
        trace(&mut sim, "flat6h", 60.0, 3.0);
        sim.set_blowtorch(0.0, 0.0, -25.0, 1.5);
        sim.set_liquid_nitrogen(1.0, 0.0, -11.0, 3.0);
        trace(&mut sim, "flatln2", 30.0, 2.0);
    }

    #[test]
    #[ignore = "scenario calibration dump"]
    fn scenario_torch_hanging_12h() {
        let mut sim = PocketWatchSim::new();
        sim.set_blowtorch(1.0, 0.0, 25.0, 1.5);
        trace(&mut sim, "torch12h", 60.0, 3.0);
        sim.set_blowtorch(0.0, 0.0, 25.0, 1.5);
        sim.set_air_jet(0.92, 0.0, 20.0, 2.5);
        trace(&mut sim, "air12h", 15.0, 1.5);
    }

    #[test]
    #[ignore = "scenario calibration dump"]
    fn scenario_torch_on_tourbillon_then_ln2() {
        let mut sim = PocketWatchSim::new();
        sim.set_blowtorch(0.42, 0.0, -22.0, 2.0);
        trace(&mut sim, "drift", 20.0, 2.0);
        sim.set_blowtorch(1.0, 0.0, -11.0, 3.0);
        trace(&mut sim, "torchT", 20.0, 2.0);
        sim.set_blowtorch(0.0, 0.0, -11.0, 3.0);
        sim.set_liquid_nitrogen(1.0, 0.0, -11.0, 3.0);
        trace(&mut sim, "cryoT", 20.0, 2.0);
    }

    #[test]
    #[ignore = "scenario calibration dump"]
    fn scenario_oil_magnet_vacuum_shock() {
        let mut sim = PocketWatchSim::new();
        sim.set_oil_pour(0.85, 0.0, -11.0, 2.0);
        trace(&mut sim, "oil", 20.0, 2.0);
        sim.set_oil_pour(0.0, 0.0, -11.0, 2.0);
        sim.clean_excess_oil();
        trace(&mut sim, "clean", 4.0, 2.0);
        sim.set_compensated_hairspring(false);
        sim.set_magnet(true, 1.0, 13.0, -11.0, 2.5);
        trace(&mut sim, "magnet", 8.0, 1.0);
        sim.set_magnet(false, 1.0, 13.0, -11.0, 2.5);
        trace(&mut sim, "remanent", 6.0, 1.0);
        sim.trigger_degauss();
        trace(&mut sim, "degauss", 6.0, 1.0);
        sim.set_ambient_pressure_atm(0.0);
        trace(&mut sim, "vacuum", 12.0, 2.0);
        sim.set_ambient_pressure_atm(1.0);
        sim.set_gravity_vector(0.0, -4.0, 0.0);
        sim.trigger_mechanical_shock(500.0);
        trace(&mut sim, "shock", 6.0, 0.5);
    }

    fn esc_line(sim: &PocketWatchSim, tag: &str) {
        let h = sim.horology();
        let e = h.escapement_state();
        let l = &h.loads;
        println!(
            "{tag} t={:7.3} th={:7.2} th_eq={:6.2} w={:7.2} phi={:6.2} psi={:8.3} mask={:05b} amp={:6.1} k={:.3e} \
             mag={:+.3} tau_lev={:+.2e} drive={:.2e} frozen={} beats={} {}",
            h.lab_time_s,
            e[0] / DEG_TEST,
            l.theta_eq / DEG_TEST,
            e[1],
            e[2] / DEG_TEST,
            e[4],
            e[6] as u32,
            h.balance_amplitude_deg,
            l.k_hairspring,
            sim.multiphysics().hairspring_magnetisation,
            l.lever_torque,
            l.escape_drive,
            h.esc.frozen,
            h.total_beats,
            h.status.label(),
        );
    }
    const DEG_TEST: f64 = core::f64::consts::PI / 180.0;

    #[test]
    #[ignore = "escapement diagnostic dump"]
    fn diag_magnet_steel_spring() {
        let mut sim = PocketWatchSim::new();
        sim.set_compensated_hairspring(false);
        for _ in 0..20 {
            sim.step(0.05);
        }
        esc_line(&sim, "pre");
        sim.set_magnet(true, 1.0, 13.0, -11.0, 2.5);
        for i in 0..80 {
            sim.step(0.05);
            if i % 4 == 0 {
                esc_line(&sim, "mag");
            }
        }
        sim.set_magnet(false, 1.0, 13.0, -11.0, 2.5);
        for i in 0..60 {
            sim.step(0.05);
            if i % 4 == 0 {
                esc_line(&sim, "off");
            }
        }
        sim.trigger_degauss();
        for i in 0..80 {
            sim.step(0.05);
            if i % 8 == 0 {
                esc_line(&sim, "dg");
            }
        }
        sim.trigger_mechanical_shock(500.0);
        for i in 0..40 {
            sim.step(0.05);
            if i % 4 == 0 {
                esc_line(&sim, "shock");
            }
        }
    }

    #[test]
    #[ignore = "escapement diagnostic dump"]
    fn diag_self_start_from_rest() {
        // Stop the balance dead at every angle across its swing and see whether
        // the escapement restarts it (a Swiss lever in beat is self-starting).
        for deg in [-40.0, -20.0, -10.0, -5.0, 0.0, 5.0, 10.0, 20.0, 40.0] {
            let mut sim = PocketWatchSim::new();
            for _ in 0..10 {
                sim.step(0.05);
            }
            sim.horology.esc.q[0] = sim.horology.loads.theta_eq + deg * DEG_TEST;
            sim.horology.esc.v = [0.0; 3];
            for _ in 0..60 {
                sim.step(0.05);
            }
            esc_line(&sim, &format!("start{deg:+}"));
        }
    }

    #[test]
    #[ignore = "front-end preset calibration"]
    fn presets_thermal_vacuum_shock_oil() {
        let mut sim = PocketWatchSim::new();
        sim.set_blowtorch(0.25, 0.0, -22.0, 6.0);
        trace(&mut sim, "drift25z6", 30.0, 5.0);
        let mut sim = PocketWatchSim::new();
        sim.set_ambient_pressure_atm(0.0);
        trace(&mut sim, "vacuum", 12.0, 3.0);
        let mut sim = PocketWatchSim::new();
        sim.set_gravity_vector(0.0, -4.0, 0.0);
        sim.trigger_mechanical_shock(500.0);
        trace(&mut sim, "shock4g", 3.0, 0.5);
        let mut sim = PocketWatchSim::new();
        sim.set_oil_pour(0.85, 0.0, -11.0, 2.0);
        trace(&mut sim, "oilflood", 8.0, 1.0);
    }

    #[test]
    #[ignore = "front-end preset calibration"]
    fn presets_melt_flat_then_ln2() {
        let mut sim = PocketWatchSim::new();
        sim.set_gravity_vector(0.0, 0.0, -1.0);
        sim.set_blowtorch(1.0, 0.0, -25.0, 1.5);
        trace(&mut sim, "meltflat", 150.0, 10.0);
        sim.set_blowtorch(0.0, 0.0, -25.0, 1.5);
        sim.set_liquid_nitrogen(1.0, 0.0, -14.0, 3.0);
        trace(&mut sim, "ln2flat", 60.0, 10.0);
    }

    #[test]
    #[ignore = "front-end preset calibration"]
    fn presets_melt_tilted_into_well() {
        for (name, g, ty) in [("tilt60", [0.0, 0.87, -0.5], -24.0), ("tilt30", [0.0, 0.5, -0.87], -24.0)] {
            let mut sim = PocketWatchSim::new();
            sim.set_gravity_vector(g[0], g[1], g[2]);
            sim.set_blowtorch(1.0, 0.0, ty, 1.5);
            trace(&mut sim, name, 120.0, 10.0);
            sim.set_blowtorch(0.0, 0.0, ty, 1.5);
            sim.set_liquid_nitrogen(1.0, 0.0, -11.0, 3.0);
            trace(&mut sim, &format!("{name}ln2"), 60.0, 10.0);
        }
    }

    #[test]
    #[ignore = "front-end preset calibration"]
    fn presets_air_blow_gold() {
        for (name, torch_on) in [("blow_torchon", true), ("blow_torchoff", false)] {
            let mut sim = PocketWatchSim::new();
            sim.set_blowtorch(1.0, 0.0, -25.0, 1.5);
            trace(&mut sim, name, 70.0, 35.0);
            if !torch_on {
                sim.set_blowtorch(0.0, 0.0, -25.0, 1.5);
            }
            sim.set_air_jet(0.92, 0.0, -24.0, 2.5);
            trace(&mut sim, &format!("{name}_jet"), 8.0, 1.0);
        }
    }

    #[test]
    #[ignore = "front-end preset calibration"]
    fn presets_magnet_compensated() {
        let mut sim = PocketWatchSim::new();
        sim.set_magnet(true, 1.0, 13.0, -11.0, 2.5);
        trace(&mut sim, "magcomp", 8.0, 1.0);
    }

    #[test]
    #[ignore = "front-end preset calibration"]
    fn presets_melt_hanging_12h() {
        let mut sim = PocketWatchSim::new();
        sim.set_blowtorch(1.0, 0.0, 25.0, 1.5);
        trace(&mut sim, "melt12", 150.0, 10.0);
    }

    #[test]
    #[ignore = "front-end preset calibration"]
    fn presets_cryo_flat() {
        let mut sim = PocketWatchSim::new();
        sim.set_gravity_vector(0.0, 0.0, -1.0);
        sim.set_liquid_nitrogen(1.0, 0.0, -11.0, 3.0);
        trace(&mut sim, "cryoflat", 90.0, 10.0);
    }

    #[test]
    fn steel_spring_swap_keeps_the_watch_running() {
        let mut sim = PocketWatchSim::new();
        sim.set_compensated_hairspring(false);
        for _ in 0..50 {
            sim.step(0.02);
        }
        let s = sim.snapshot();
        assert!(s[13] > 150.0, "amplitude {}", s[13]);
        assert!(sim.status_label().starts_with("RUNNING"), "{}", sim.status_label());
    }
}
