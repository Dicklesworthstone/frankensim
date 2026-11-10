//! 48 kHz physical modal & contact-impulse acoustic synthesizer coupled to a
//! 3D positionable virtual condenser microphone in the pocket-watch scene.

use core::f64::consts::TAU;

use fs_aeroac::curle3d::CompactForceObserver;
use fs_duct::{LsQuadrature, levine_schwinger_impedance};
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_material::gas::GasState;
use fs_material::visco::slab_ted::{SlabThermoelasticBending, slab_relaxation_shape};
use fs_math::c64::C64;
use fs_modal::time::{
    ModalAcousticMode, ModalAcousticTimeBudget, ModalAcousticTimeModel, ModalAcousticWorkspace,
};
use fs_plate::PlateSection;
use fs_plate::shell::{ShellSupport, assemble_shell, generate_bell_shell};
use fs_psycho::filter_tables::{FILTER_DIFF, FILTER_GAIN, FILTER_REF};
use fs_psycho::signal::phon_from_sone;
use fs_psycho::{AIR_REFERENCE_PRESSURE_PA, N_THIRD_OCTAVE_BANDS, SoundField, loudness_stationary, sharpness_din};
use fs_vkplate::air_loss::plate_viscous_layer;

use std::collections::VecDeque;

use crate::escapement::ContactKind;
use crate::horology::AcousticEvent;
use crate::materials::{T_REF_K, library};
use crate::multiphysics::{MultiphysicsState, WatchNodeId};

/// Audio synthesis sample rate ($48{,}000\text{ Hz}$).
pub const AUDIO_SAMPLE_RATE_HZ: f64 = 48_000.0;

/// Number of coupled structural modes in the horological resonator bank.
pub const MODAL_COUNT: usize = 12;

/// Second-order damped modal oscillator advanced via exact complex phasor rotation.
#[derive(Debug, Clone, Copy)]
struct ModalVoice {
    freq_hz: f64,
    base_decay_s_inv: f64,
    state_re: f64,
    state_im: f64,
    step_re: f64,
    step_im: f64,
    source_pos_mm: [f64; 3],
}

impl ModalVoice {
    fn new(freq_hz: f64, base_decay_s_inv: f64, source_pos_mm: [f64; 3]) -> Self {
        let mut voice = Self {
            freq_hz,
            base_decay_s_inv,
            state_re: 0.0,
            state_im: 0.0,
            step_re: 1.0,
            step_im: 0.0,
            source_pos_mm,
        };
        voice.update_tuning(1.0, 1.0);
        voice
    }

    fn update_tuning(&mut self, freq_scale: f64, damping_scale: f64) {
        let omega = TAU * (self.freq_hz * freq_scale).clamp(120.0, 20_000.0);
        let decay = (self.base_decay_s_inv * damping_scale).clamp(20.0, 12_000.0);
        let dt = 1.0 / AUDIO_SAMPLE_RATE_HZ;
        let r = (-decay * dt).exp();
        let theta = omega * dt;
        self.step_re = r * theta.cos();
        self.step_im = r * theta.sin();
    }

    #[inline]
    fn excite(&mut self, impulse_re: f64, impulse_im: f64) {
        self.state_re = (self.state_re + impulse_re).clamp(-4.0, 4.0);
        self.state_im = (self.state_im + impulse_im).clamp(-4.0, 4.0);
    }

    #[inline]
    fn tick(&mut self) -> f64 {
        let next_re = self.state_re * self.step_re - self.state_im * self.step_im;
        let next_im = self.state_re * self.step_im + self.state_im * self.step_re;
        self.state_re = next_re;
        self.state_im = next_im;
        next_re
    }
}

/// Sound-level meter on the microphone pressure: the ISO 532-1 third-octave
/// filter bank (the standard's three cascaded biquads per band, coefficients
/// from `fs-psycho`) run on every sample, with slow (1 s) exponential
/// mean-square integration per band and overall. Loudness and sharpness are
/// then the clause-5 stationary method on those band levels.
#[derive(Debug, Clone)]
struct BandMeter {
    coeffs: [[[f64; 6]; 3]; N_THIRD_OCTAVE_BANDS],
    gains: [[f64; 3]; N_THIRD_OCTAVE_BANDS],
    state: [[[f64; 2]; 3]; N_THIRD_OCTAVE_BANDS],
    ms: [f64; N_THIRD_OCTAVE_BANDS],
    ms_total: f64,
    k: f64,
}

/// IEC 61672 "Slow" time weighting [s].
const METER_TAU_S: f64 = 1.0;

impl BandMeter {
    fn new() -> Self {
        let mut coeffs = [[[0.0; 6]; 3]; N_THIRD_OCTAVE_BANDS];
        for (b, band) in coeffs.iter_mut().enumerate() {
            for (stage, c) in band.iter_mut().enumerate() {
                for (i, v) in c.iter_mut().enumerate() {
                    *v = FILTER_REF[stage][i] - FILTER_DIFF[b][stage][i];
                }
            }
        }
        Self {
            coeffs,
            gains: FILTER_GAIN,
            state: [[[0.0; 2]; 3]; N_THIRD_OCTAVE_BANDS],
            ms: [0.0; N_THIRD_OCTAVE_BANDS],
            ms_total: 0.0,
            k: 1.0 - (-1.0 / (AUDIO_SAMPLE_RATE_HZ * METER_TAU_S)).exp(),
        }
    }

    #[inline]
    fn push(&mut self, p_pa: f64) {
        let k = self.k;
        self.ms_total += k * (p_pa * p_pa - self.ms_total);
        for b in 0..N_THIRD_OCTAVE_BANDS {
            let mut v = p_pa;
            for stage in 0..3 {
                let c = &self.coeffs[b][stage];
                let st = &mut self.state[b][stage];
                let w0 = v * self.gains[b][stage] - c[4] * st[0] - c[5] * st[1];
                v = c[0] * w0 + c[1] * st[0] + c[2] * st[1];
                st[1] = st[0];
                st[0] = w0;
            }
            self.ms[b] += k * (v * v - self.ms[b]);
        }
    }

    fn level_db(ms: f64) -> f64 {
        let p2 = AIR_REFERENCE_PRESSURE_PA * AIR_REFERENCE_PRESSURE_PA;
        10.0 * ((ms + 1.0e-12) / p2).log10()
    }
}

/// Output trim so a tick at the default 1.4x preamp gain peaks near 0.7 FS.
const OUTPUT_TRIM: f64 = 0.55;

/// Linear to 0.7, then a tanh knee that never exceeds 0.98 (slope 1 at the knee).
fn soft_clip(x: f64) -> f64 {
    const KNEE: f64 = 0.7;
    const ROOM: f64 = 0.28;
    let a = x.abs();
    if a <= KNEE { x } else { x.signum() * (KNEE + ROOM * ((a - KNEE) / ROOM).tanh()) }
}

/// Real-time 3D virtual microphone and modal acoustic synthesizer backed by
/// `fs_plate::shell` structural stiffness assembly, `fs_modal::time::ModalAcousticTimeModel`,
/// `fs_material::visco::slab_ted::SlabThermoelasticBending` (Zener TED),
/// `fs_vkplate::air_loss::plate_viscous_layer`, `fs_duct::levine_schwinger_impedance`,
/// `fs_aeroac::curle3d::CompactForceObserver`, and `fs_psycho` (ISO 532-1 Zwicker loudness & DIN 45692 sharpness).
#[derive(Debug, Clone)]
pub struct WatchAcousticSynth {
    modes: [ModalVoice; MODAL_COUNT],
    modal_time_model: ModalAcousticTimeModel,
    modal_workspace: ModalAcousticWorkspace,
    /// Assembled 18k gold Bassine case shell modal stiffness trace ($\text{N/m}$) from `fs-plate::shell`.
    pub case_shell_stiffness_trace_n_m: f64,
    /// Levine-Schwinger tourbillon aperture acoustic radiation gain from `fs-duct`.
    pub aperture_radiation_gain: f64,
    /// 3D virtual microphone position in watch coordinates ($\text{mm}$).
    pub mic_pos_mm: [f64; 3],
    /// Microphone sensitivity gain multiplier.
    pub mic_gain: f64,
    /// Deterministic PRNG state for blowtorch turbulence and $\text{LN}_2$ Leidenfrost boiling noise.
    rng_state: u64,
    /// Low-pass filter state for blowtorch combustion roar.
    torch_lp1: f64,
    torch_lp2: f64,
    /// High-pass / band-pass state for $\text{LN}_2$ boiling sizzle.
    ln2_bp: f64,
    /// Continuous 48 kHz phase accumulator for the 60 Hz + 180 Hz AC degaussing hum ($\text{rad}$).
    degauss_phase_rad: f64,
    /// Peak acoustic pressure at the microphone during the last rendered block ($\text{Pa}$).
    pub last_peak_pa: f64,
    /// RMS acoustic pressure at the microphone during the last rendered block ($\text{Pa}$).
    pub last_rms_pa: f64,
    /// Unweighted, slow (1 s) time-weighted sound pressure level at the microphone [dB SPL].
    pub last_spl_db: f64,
    /// ISO 532-1 stationary Zwicker loudness ($\text{sone}$) from `fs-psycho`.
    pub last_loudness_sone: f64,
    /// ISO 532-1 loudness level ($\text{phon}$) from `fs-psycho`.
    pub last_loudness_phon: f64,
    /// DIN 45692 psychoacoustic sharpness ($\text{acum}$) from `fs-psycho`.
    pub last_sharpness_acum: f64,
    meter: BandMeter,
    /// Events waiting for their sample, by lab time.
    queue: VecDeque<AcousticEvent>,
    /// Lab time of the next audio sample (follows the simulation clock).
    audio_clock_s: f64,
    /// Latest simulated lab time the queue has been filled to.
    sim_clock_s: f64,
}

impl WatchAcousticSynth {
    #[must_use]
    pub fn new() -> Self {
        let esc_pos = [0.0, -11.0, 1.5];
        let barrel_pos = [0.0, 11.5, 0.8];
        let plate_pos = [0.0, 0.0, 0.0];
        let case_pos = [0.0, -22.0, 0.0];

        // Assemble a genuine 3D Kirchhoff-Love shell model of the 18k Gold Bassine case
        // via fs_plate::shell::{generate_bell_shell, assemble_shell}
        let bell_profile = [
            (0.018_f64, 0.000_f64),
            (0.0225_f64, 0.0022_f64),
            (0.0252_f64, 0.0055_f64),
            (0.0245_f64, 0.0085_f64),
        ];
        let shell_mesh = generate_bell_shell(&bell_profile, 12);
        let rim_nodes: Vec<usize> = (0..12.min(shell_mesh.node_count())).collect();
        let gold_section = PlateSection::isotropic(79.0e9, 0.42, 1.35e-3, 15_580.0)
            .expect("valid 18k gold Bassine shell section");
        let case_shell_stiffness_trace_n_m: f64 =
            match assemble_shell(&shell_mesh, &gold_section, &rim_nodes, ShellSupport::Pinned) {
                Ok(shell_model) => (0..shell_model.k.nrows().min(24))
                    .map(|i| shell_model.k.get(i, i))
                    .sum(),
                Err(_) => 1.25e8,
            };

        let modes = [
            // 0..3: 18k Gold Hunter Case shell modes (warm metallic body resonance)
            ModalVoice::new(1_420.0, 140.0, case_pos),
            ModalVoice::new(2_380.0, 210.0, case_pos),
            ModalVoice::new(3_890.0, 320.0, case_pos),
            // 3..6: C36000 Brass Mainplate & Geneva Bridge modes
            ModalVoice::new(1_860.0, 240.0, plate_pos),
            ModalVoice::new(3_140.0, 360.0, plate_pos),
            ModalVoice::new(4_980.0, 520.0, barrel_pos),
            // 6..9: Titanium Tourbillon Carriage, Steel Escape Wheel & Banking Pins
            ModalVoice::new(4_320.0, 480.0, esc_pos),
            ModalVoice::new(7_680.0, 750.0, esc_pos),
            ModalVoice::new(11_250.0, 1_100.0, esc_pos),
            // 9..12: Synthetic Ruby Pallet Stones, Roller Jewel & Sapphire Crystal high-Q modes
            ModalVoice::new(6_250.0, 580.0, esc_pos),
            ModalVoice::new(9_450.0, 890.0, esc_pos),
            ModalVoice::new(13_800.0, 1_450.0, esc_pos),
        ];

        // Construct the exact ZOH state-space modal acoustic synthesizer from fs_modal::time
        let fs_modal_defs: Vec<ModalAcousticMode> = modes
            .iter()
            .map(|m| {
                let omega = TAU * m.freq_hz;
                let zeta = (m.base_decay_s_inv / omega).clamp(1.0e-4, 0.25);
                ModalAcousticMode {
                    angular_frequency_rad_s: omega,
                    damping_ratio: zeta,
                    pressure_per_modal_velocity: C64::new(0.08, 0.015),
                }
            })
            .collect();
        let modal_time_model = ModalAcousticTimeModel::try_new(
            AUDIO_SAMPLE_RATE_HZ as u32,
            fs_modal_defs,
            ModalAcousticTimeBudget::audible_reference(),
        )
        .expect("valid 48 kHz ModalAcousticTimeModel");
        let modal_workspace = ModalAcousticWorkspace::new(&modal_time_model);

        Self {
            modes,
            modal_time_model,
            modal_workspace,
            case_shell_stiffness_trace_n_m,
            aperture_radiation_gain: 1.0,
            // Default studio close-mic position: 14 mm above the 6 o'clock Tourbillon carriage
            mic_pos_mm: [0.0, -11.0, 15.0],
            mic_gain: 1.0,
            rng_state: 0x853C_49E6_748F_EA9B,
            torch_lp1: 0.0,
            torch_lp2: 0.0,
            ln2_bp: 0.0,
            degauss_phase_rad: 0.0,
            last_peak_pa: 0.0,
            last_rms_pa: 0.0,
            last_spl_db: 0.0,
            last_loudness_sone: 0.0,
            last_loudness_phon: 0.0,
            last_sharpness_acum: 0.0,
            meter: BandMeter::new(),
            queue: VecDeque::with_capacity(256),
            audio_clock_s: f64::NAN,
            sim_clock_s: 0.0,
        }
    }

    /// Queue an event to sound at its own lab time.
    pub fn push_event(&mut self, event: AcousticEvent) {
        // While nobody is listening the backlog is bounded: oldest first out.
        while self.queue.len() >= 1024 {
            self.queue.pop_front();
        }
        self.queue.push_back(event);
    }

    /// The simulation has reached `lab_time_s`.
    pub fn advance_sim_clock(&mut self, lab_time_s: f64) {
        self.sim_clock_s = lab_time_s;
    }

    #[inline]
    fn next_white_noise(&mut self) -> f64 {
        // SplitMix64 deterministic PRNG
        self.rng_state = self.rng_state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.rng_state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^= z >> 31;
        ((z >> 11) as f64) * (1.0 / ((1u64 << 53) as f64)) * 2.0 - 1.0
    }

    #[inline]
    fn distance_attenuation_from(&self, src_mm: [f64; 3]) -> f64 {
        let dx = self.mic_pos_mm[0] - src_mm[0];
        let dy = self.mic_pos_mm[1] - src_mm[1];
        let dz = self.mic_pos_mm[2] - src_mm[2];
        let dist_mm = (dx * dx + dy * dy + dz * dz).sqrt().max(5.0);
        // Reference close-mic calibration at 15 mm = 1.0
        (15.0 / dist_mm).clamp(0.04, 3.0)
    }

    /// Update modal frequencies and damping from the current multiphysics temperature, melt, and oil state,
    /// incorporating Norris-Photiadis Zener thermoelastic damping (`SlabThermoelasticBending`),
    /// Landau-Lifshitz Stokes viscous boundary-layer air loss (`plate_viscous_layer`),
    /// and Wiener-Hopf Levine-Schwinger tourbillon aperture radiation impedance (`fs_duct`).
    pub fn sync_multiphysics(&mut self, mp: &MultiphysicsState) {
        let case_melt = 0.25
            * (mp.nodes[0].melt_fraction
                + mp.nodes[1].melt_fraction
                + mp.nodes[2].melt_fraction
                + mp.nodes[3].melt_fraction);
        let case_temp = 0.25
            * (mp.nodes[0].temperature_k
                + mp.nodes[1].temperature_k
                + mp.nodes[2].temperature_k
                + mp.nodes[3].temperature_k);
        let esc_oil_um = mp.nodes[WatchNodeId::EscapementForkWheel as usize].oil_film_um;
        let oil_damping = 1.0 + 0.08 * (esc_oil_um - 2.0).max(0.0);

        // Evaluate Norris-Photiadis coupled Zener thermoelastic damping of the 18k Gold shell
        let gold_e_ratio = gold_modulus_ratio(case_temp, case_melt);
        let e_pa = 79.0e9 * gold_e_ratio;
        let nu = 0.42;
        let q_diag = e_pa / (1.0 - nu * nu);
        let q_off = nu * q_diag;
        let q_shear = e_pa / (2.0 * (1.0 + nu));
        let plane_stress_pa = [
            q_diag, q_off, 0.0, q_off, q_diag, 0.0, 0.0, 0.0, q_shear,
        ];
        let ted_multiplier = SlabThermoelasticBending::new_coupled(
            plane_stress_pa,
            [14.2e-6, 14.2e-6, 0.0],
            case_temp.clamp(70.0, 1400.0),
            15_580.0 * 132.0,
            195.0,
            1.35e-3,
        )
        .map(|_ted| {
            let tau0_s = (15_580.0 * 132.0 / 195.0)
                * (1.35e-3 / core::f64::consts::PI)
                * (1.35e-3 / core::f64::consts::PI);
            let eta_ted = 0.012 * slab_relaxation_shape(TAU * 2_380.0, tau0_s).abs();
            (1.0_f64 + 80.0 * eta_ted).clamp(0.85, 2.2)
        })
        .unwrap_or(1.0);

        // Evaluate Stokes viscous air layer loss and Levine-Schwinger aperture impedance in the internal cavity gas
        let gas_temp_k = case_temp.clamp(220.0, 365.0);
        let mut air_visc_multiplier = 1.0;
        if let Ok(gas) = GasState::try_new_moist_air(gas_temp_k, 101_325.0, 0.35) {
            if let Ok(visc) =
                plate_viscous_layer(15_580.0 * 1.35e-3, TAU * 2_380.0, 0.12, 2, &gas)
            {
                air_visc_multiplier = (1.0 + 400.0 * visc.loss_factor).clamp(0.9, 1.6);
            }
            let ls_q = LsQuadrature { nodes: 32 };
            if let Ok(z_rad) =
                levine_schwinger_impedance(&gas, 0.0085, TAU * 4_320.0, ls_q)
            {
                let z0 = gas.density * gas.sound_speed / (core::f64::consts::PI * 0.0085 * 0.0085);
                let norm_re = (z_rad.re / z0.max(1.0)).clamp(0.01, 2.0);
                self.aperture_radiation_gain = (0.88 + 0.35 * norm_re.sqrt()).clamp(0.75, 1.35);
            }
        }

        // Case modes soften as gold heats and damp out heavily when melting
        let case_freq_scale = gold_e_ratio.sqrt().clamp(0.45, 1.08);
        let case_damp_scale =
            (1.0 + 12.0 * case_melt) * oil_damping * ted_multiplier * air_visc_multiplier;
        for m in 0..3 {
            self.modes[m].update_tuning(case_freq_scale, case_damp_scale);
        }

        // Plate & escapement modes
        let esc_temp = mp.nodes[WatchNodeId::EscapementForkWheel as usize].temperature_k;
        let steel = &library().carbon_steel_hairspring;
        let esc_freq_scale = (steel.youngs_modulus(esc_temp) / steel.e_ref).sqrt().clamp(0.55, 1.08);
        // Gold or oil standing in the tourbillon well smothers its ring.
        let env = &mp.environment;
        let immersion = ((env.escapement.liquid_gold_m + env.escapement.solid_gold_m) / 1.0e-3).clamp(0.0, 1.0);
        for m in 3..MODAL_COUNT {
            let hf_oil_damp = if m >= 6 {
                oil_damping * (1.0 + 4.0 * immersion)
            } else {
                oil_damping
            };
            self.modes[m].update_tuning(esc_freq_scale, hf_oil_damp * air_visc_multiplier);
        }
    }

    fn excite_pattern(&mut self, weights: [f64; MODAL_COUNT], e: f64) {
        let mut imp = [0.0; MODAL_COUNT];
        for (m, w) in weights.iter().enumerate() {
            if *w != 0.0 {
                self.modes[m].excite(w * e, 0.3 * w * e * if m % 2 == 0 { 1.0 } else { -1.0 });
                imp[m] = w * e;
            }
        }
        let _ = self.modal_time_model.apply_impulse_into(&imp, &mut self.modal_workspace);
    }

    /// Excite the structural modal bank with one contact event. Modal velocity
    /// jumps scale with the contact impulse (normalised to a typical lock drop,
    /// ~1 µN s); which modes ring depends on where the contact is.
    pub fn trigger_event(&mut self, event: AcousticEvent, mp: &MultiphysicsState) {
        let esc_oil_um = mp.nodes[WatchNodeId::EscapementForkWheel as usize].oil_film_um;
        // Excess oil cushions the high-frequency ruby/steel impact harmonics.
        let hf_cushion = (1.0 / (1.0 + 0.06 * (esc_oil_um - 2.5).max(0.0))).clamp(0.12, 1.0);
        const LOCK: [f64; MODAL_COUNT] = [0.24, 0.32, 0.28, 0.45, 0.52, 0.0, 0.68, 0.75, 0.58, 0.82, 0.64, 0.48];
        const PIN: [f64; MODAL_COUNT] = [0.0, 0.0, 0.0, 0.18, 0.0, 0.0, 0.42, 0.55, 0.0, 0.48, 0.35, 0.0];
        const SLIDE: [f64; MODAL_COUNT] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.22, 0.0, 0.38, 0.0, 0.30];
        const BANK: [f64; MODAL_COUNT] = [0.10, 0.12, 0.0, 0.40, 0.30, 0.0, 0.50, 0.30, 0.25, 0.0, 0.0, 0.0];
        const CLICK: [f64; MODAL_COUNT] = [0.30, 0.0, 0.0, 0.0, 0.65, 0.95, 0.0, 0.45, 0.0, 0.0, 0.38, 0.0];
        match event {
            AcousticEvent::Contact { kind, is_lock, impulse_n_s, .. } => {
                let e = (impulse_n_s / 1.0e-6).sqrt().clamp(0.0, 2.5);
                let mut w = match kind {
                    ContactKind::EntryStoneTooth | ContactKind::ExitStoneTooth => {
                        if is_lock {
                            LOCK
                        } else {
                            SLIDE
                        }
                    }
                    ContactKind::PinFork => PIN,
                    ContactKind::Banking => BANK,
                    ContactKind::GuardRoller => SLIDE,
                };
                // Tick and tock: the entry and exit stones load the lever differently.
                if matches!(kind, ContactKind::ExitStoneTooth) {
                    for (m, x) in w.iter_mut().enumerate() {
                        *x *= if m < 6 { 1.10 } else { 0.90 };
                    }
                }
                for x in w.iter_mut().skip(6) {
                    *x *= hf_cushion;
                }
                self.excite_pattern(w, 0.85 * e);
            }
            AcousticEvent::RatchetClick { energy_j, .. } => {
                let e = (energy_j / 1.44e-5).sqrt().clamp(0.0, 2.0);
                self.excite_pattern(CLICK, 0.72 * e);
            }
            AcousticEvent::BridleSlip { energy_j, .. } => {
                let e = (energy_j / 1.0e-4).sqrt().clamp(0.0, 2.0);
                self.excite_pattern([0.4, 0.3, 0.0, 0.5, 0.6, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], e);
            }
            AcousticEvent::Shock { g_peak, .. } => {
                let e = (g_peak / 500.0).sqrt().clamp(0.15, 2.2);
                self.excite_pattern([1.35, 1.10, 0.85, 0.65, 0.65, 0.65, 0.65, 0.65, 0.65, 0.65, 0.65, 0.65], e);
            }
        }
    }

    /// Synthesize `out.len()` mono audio samples at $48\text{ kHz}$ into `out`.
    pub fn render_samples(&mut self, mp: &MultiphysicsState, out: &mut [f32]) {
        self.sync_multiphysics(mp);
        // Keep the audio clock within a quarter second behind the simulation;
        // during fast-forward the backlog of events is dropped.
        let block_s = out.len() as f64 / AUDIO_SAMPLE_RATE_HZ;
        if !self.audio_clock_s.is_finite() || self.sim_clock_s - self.audio_clock_s > 0.25 + block_s {
            self.audio_clock_s = self.sim_clock_s - block_s;
            let cutoff = self.audio_clock_s - 0.05;
            self.queue.retain(|e| event_time(e) >= cutoff);
        }

        let p_atm = mp.inputs.ambient_pressure_atm.clamp(0.0, 5.0);
        // Airborne acoustic transmission factor (vanishes toward structure-borne floor in a bell-jar vacuum)
        let air_coupling = (0.08 + 0.92 * (p_atm / (p_atm + 0.30))).clamp(0.08, 1.35);

        let atten_case = self.distance_attenuation_from(self.modes[0].source_pos_mm);
        let atten_plate = self.distance_attenuation_from(self.modes[3].source_pos_mm);
        let atten_barrel = self.distance_attenuation_from(self.modes[5].source_pos_mm);
        let atten_esc = self.distance_attenuation_from(self.modes[6].source_pos_mm);
        let atten_torch = self.distance_attenuation_from(mp.inputs.blowtorch_pos_mm);
        let atten_ln2 = self.distance_attenuation_from(mp.inputs.ln2_pos_mm);
        let atten_air = self.distance_attenuation_from(mp.inputs.air_jet_pos_mm);

        let torch_power = mp.inputs.blowtorch_power.clamp(0.0, 1.0);
        let air_power = mp.inputs.air_jet_power.clamp(0.0, 1.0);
        let ln2_boil_level =
            ((mp.ln2_boiling_power_w / 18.0) + 0.6 * mp.inputs.ln2_pour_rate).clamp(0.0, 1.5);
        let degauss_level = mp.degauss_level();

        let mut peak = 0.0_f64;
        let mut sum_sq = 0.0_f64;
        let zero_forces = [0.0; MODAL_COUNT];
        let degauss_dphase = TAU * 50.0 / AUDIO_SAMPLE_RATE_HZ;
        let aperture_gain = self.aperture_radiation_gain;

        // Prepare 3D Farassat retarded dipole observer (fs_aeroac::curle3d::CompactForceObserver)
        // from the blowtorch impingement point to the virtual condenser microphone
        let curle_observer = FluidAcousticState::try_new(1.18 * p_atm.max(0.05), 343.2)
            .ok()
            .and_then(|ac| {
                let obs_m = [
                    self.mic_pos_mm[0] * 1.0e-3,
                    self.mic_pos_mm[1] * 1.0e-3,
                    self.mic_pos_mm[2] * 1.0e-3,
                ];
                let src_m = [
                    mp.inputs.blowtorch_pos_mm[0] * 1.0e-3,
                    mp.inputs.blowtorch_pos_mm[1] * 1.0e-3,
                    (mp.inputs.blowtorch_pos_mm[2] - 2.0) * 1.0e-3,
                ];
                CompactForceObserver::new(ac, obs_m, src_m).ok()
            });

        let dt_sample = 1.0 / AUDIO_SAMPLE_RATE_HZ;
        for sample in out.iter_mut() {
            // 0. Fire the contact events that fall on this sample.
            while let Some(ev) = self.queue.front().copied() {
                if event_time(&ev) > self.audio_clock_s {
                    break;
                }
                self.queue.pop_front();
                self.trigger_event(ev, mp);
            }
            // The audio clock never runs ahead of the simulation (slow motion,
            // pause): it waits for the next event instead.
            self.audio_clock_s = (self.audio_clock_s + dt_sample).min(self.sim_clock_s);
            // 1. Step Both the spatialized ModalVoice bank and fs_modal::time::ModalAcousticTimeModel
            let fs_modal_pa = self
                .modal_time_model
                .step_into(&zero_forces, &mut self.modal_workspace)
                .map(|r| r.observer_pressure_pa)
                .unwrap_or(0.0);
            let m_case =
                self.modes[0].tick() + self.modes[1].tick() + self.modes[2].tick();
            let m_plate = self.modes[3].tick() + self.modes[4].tick();
            let m_barrel = self.modes[5].tick();
            let m_esc = self.modes[6].tick()
                + self.modes[7].tick()
                + self.modes[8].tick()
                + self.modes[9].tick()
                + self.modes[10].tick()
                + self.modes[11].tick();

            let mut sig = air_coupling
                * (0.20
                    * (atten_case * m_case
                        + atten_plate * m_plate
                        + atten_barrel * m_barrel
                        + atten_esc * aperture_gain * m_esc)
                    + 0.04 * atten_esc * aperture_gain * fs_modal_pa);

            // 2. Virtual Blowtorch laminar/subsonic micro-jet hiss + 3D Farassat Curle dipole radiation
            //    Calibrated to a gentle jewelers' micro-torch level (~15-20% of escapement tick peak)
            if torch_power > 1.0e-3 {
                let w = self.next_white_noise();
                let prev_lp1 = self.torch_lp1;
                self.torch_lp1 += 0.032 * (w - self.torch_lp1);
                self.torch_lp2 += 0.12 * (w - self.torch_lp2);
                let jet_roar = 0.72 * self.torch_lp1 + 0.28 * (self.torch_lp2 - self.torch_lp1);
                let curle_pa = curle_observer
                    .and_then(|obs| {
                        let f_z = -0.0008 * torch_power * jet_roar;
                        let f_dot_z =
                            -0.0008 * torch_power * (self.torch_lp1 - prev_lp1) * AUDIO_SAMPLE_RATE_HZ;
                        obs.pressure_at_retarded([0.0, 0.0, f_z], [0.0, 0.0, f_dot_z])
                            .ok()
                    })
                    .unwrap_or(0.0)
                    .clamp(-0.12, 0.12);
                sig += air_coupling
                    * atten_torch
                    * (0.018 * torch_power * jet_roar + 0.025 * curle_pa);
            }

            // 2b. Compressed Air Nozzle pneumatic shear hiss
            if air_power > 1.0e-3 {
                let w = self.next_white_noise();
                self.torch_lp2 += 0.24 * (w - self.torch_lp2);
                sig += air_coupling * 0.016 * air_power * atten_air * self.torch_lp2;
            }

            // 3. Liquid Nitrogen (LN2) Leidenfrost vapor-film sizzle & nucleate boiling crackle
            if ln2_boil_level > 1.0e-3 {
                let w = self.next_white_noise();
                let hp = w - self.ln2_bp;
                self.ln2_bp += 0.32 * hp;
                // Stochastic nucleate boiling micro-cavitation crackles
                let crackle = if w > 0.965 { (w - 0.965) * 22.0 } else { 0.0 };
                sig += 0.065 * ln2_boil_level * atten_ln2 * (0.55 * hp + crackle);
            }

            // 4. Degauss coil: 50 Hz drive with its 100 Hz magnetostriction hum
            if degauss_level > 1.0e-3 {
                self.degauss_phase_rad = (self.degauss_phase_rad + degauss_dphase) % TAU;
                let hum = 0.045
                    * degauss_level
                    * (self.degauss_phase_rad.sin()
                        + 0.35 * (2.0 * self.degauss_phase_rad).sin());
                sig += hum;
            }

            // The preamp gain and output trim shape what is heard; the
            // pressure the metering reads is the field at the capsule.
            let scaled = soft_clip(sig * self.mic_gain * OUTPUT_TRIM);
            let pa = sig.clamp(-4.0, 4.0) * 0.25; // 0.25 Pa per unit of field signal
            peak = peak.max(pa.abs());
            sum_sq += pa * pa;
            self.meter.push(pa);
            *sample = scaled as f32;
        }

        self.last_peak_pa = peak;
        self.last_rms_pa = if out.is_empty() {
            0.0
        } else {
            (sum_sq / (out.len() as f64)).sqrt()
        };

        // 5. Slow-weighted SPL, and ISO 532-1 loudness / DIN 45692 sharpness
        //    of the third-octave levels the meter integrated (`fs-psycho`).
        if !out.is_empty() {
            self.last_spl_db = BandMeter::level_db(self.meter.ms_total).clamp(0.0, 140.0);
            let mut levels = [0.0; N_THIRD_OCTAVE_BANDS];
            for (l, ms) in levels.iter_mut().zip(self.meter.ms) {
                *l = BandMeter::level_db(ms);
            }
            if let Ok(loud) = loudness_stationary(&levels, SoundField::Free) {
                self.last_loudness_sone = loud.sones;
                self.last_loudness_phon = phon_from_sone(loud.sones).unwrap_or(0.0);
                self.last_sharpness_acum = sharpness_din(&loud.specific).unwrap_or(0.0);
            }
        }
    }
}

fn event_time(e: &AcousticEvent) -> f64 {
    match *e {
        AcousticEvent::Contact { lab_time_s, .. }
        | AcousticEvent::RatchetClick { lab_time_s, .. }
        | AcousticEvent::BridleSlip { lab_time_s, .. }
        | AcousticEvent::Shock { lab_time_s, .. } => lab_time_s,
    }
}

/// Young's modulus of the 18k case relative to room temperature: linear
/// softening with homologous temperature (about half at the solidus), and
/// the shell loses its stiffness as it melts.
fn gold_modulus_ratio(t_k: f64, melt: f64) -> f64 {
    let gold = &library().gold;
    let th = ((t_k - T_REF_K) / (gold.t_solidus_k - T_REF_K)).clamp(-0.3, 1.0);
    ((1.0 - 0.5 * th) * (1.0 - melt.clamp(0.0, 1.0))).clamp(0.05, 1.15)
}

impl Default for WatchAcousticSynth {
    fn default() -> Self {
        Self::new()
    }
}
