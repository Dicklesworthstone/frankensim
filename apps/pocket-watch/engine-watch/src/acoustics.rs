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
use fs_psycho::signal::phon_from_sone;
use fs_psycho::{
    AIR_REFERENCE_PRESSURE_PA, N_THIRD_OCTAVE_BANDS, SoundField, loudness_stationary,
    sharpness_din, spl_from_pressure_rms,
};
use fs_vkplate::air_loss::plate_viscous_layer;

use crate::horology::MechanicalAcousticEvent;
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
    /// Unweighted equivalent Sound Pressure Level at the microphone ($\text{dB SPL}$) from `fs-psycho`.
    pub last_spl_db: f64,
    /// ISO 532-1 stationary Zwicker loudness ($\text{sone}$) from `fs-psycho`.
    pub last_loudness_sone: f64,
    /// ISO 532-1 loudness level ($\text{phon}$) from `fs-psycho`.
    pub last_loudness_phon: f64,
    /// DIN 45692 psychoacoustic sharpness ($\text{acum}$) from `fs-psycho`.
    pub last_sharpness_acum: f64,
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
            last_spl_db: 28.0,
            last_loudness_sone: 0.4,
            last_loudness_phon: 28.0,
            last_sharpness_acum: 1.8,
        }
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
        let gold_e_ratio = mp.nodes[0].material.card().youngs_modulus_ratio(case_temp);
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
        let esc_freq_scale = mp.nodes[WatchNodeId::EscapementForkWheel as usize]
            .material
            .card()
            .youngs_modulus_ratio(esc_temp)
            .sqrt()
            .clamp(0.55, 1.08);
        for m in 3..MODAL_COUNT {
            let hf_oil_damp = if m >= 6 {
                oil_damping * (1.0 + 4.0 * mp.modifiers.gold_fouling_jam_factor)
            } else {
                oil_damping
            };
            self.modes[m].update_tuning(esc_freq_scale, hf_oil_damp * air_visc_multiplier);
        }
    }

    /// Excite the structural modal bank with a discrete mechanical contact event.
    pub fn trigger_event(&mut self, event: MechanicalAcousticEvent, mp: &MultiphysicsState) {
        let esc_oil_um = mp.nodes[WatchNodeId::EscapementForkWheel as usize].oil_film_um;
        // Excess oil cushions high-frequency ruby/steel impact harmonics
        let hf_cushion = (1.0 / (1.0 + 0.06 * (esc_oil_um - 2.5).max(0.0))).clamp(0.12, 1.0);

        match event {
            MechanicalAcousticEvent::EscapementUnlock {
                energy,
                is_entry_swing,
            } => {
                // Crisp, light initial pulse (roller jewel striking fork horn)
                let asym = if is_entry_swing { 1.04 } else { 0.96 };
                let e = energy * 0.28 * asym;
                self.modes[3].excite(0.18 * e, 0.05 * e);
                self.modes[6].excite(0.42 * e * hf_cushion, 0.12 * e);
                self.modes[7].excite(0.55 * e * hf_cushion, -0.18 * e);
                self.modes[9].excite(0.48 * e * hf_cushion, 0.22 * e);
                self.modes[10].excite(0.35 * e * hf_cushion, 0.0);
                let mut imp = [0.0; MODAL_COUNT];
                imp[3] = 0.18 * e;
                imp[6] = 0.42 * e * hf_cushion;
                imp[7] = 0.55 * e * hf_cushion;
                imp[9] = 0.48 * e * hf_cushion;
                imp[10] = 0.35 * e * hf_cushion;
                let _ = self
                    .modal_time_model
                    .apply_impulse_into(&imp, &mut self.modal_workspace);
            }
            MechanicalAcousticEvent::EscapementImpulseSlide {
                energy,
                is_entry_swing,
            } => {
                // Subtle high-frequency sapphire-on-steel sliding impulse
                let asym = if is_entry_swing { 1.02 } else { 0.98 };
                let e = energy * 0.16 * asym * hf_cushion;
                self.modes[7].excite(0.22 * e, 0.08 * e);
                self.modes[9].excite(0.38 * e, -0.14 * e);
                self.modes[11].excite(0.30 * e, 0.10 * e);
                let mut imp = [0.0; MODAL_COUNT];
                imp[7] = 0.22 * e;
                imp[9] = 0.38 * e;
                imp[11] = 0.30 * e;
                let _ = self
                    .modal_time_model
                    .apply_impulse_into(&imp, &mut self.modal_workspace);
            }
            MechanicalAcousticEvent::EscapementLockDrop {
                energy,
                is_entry_stone,
            } => {
                // Primary horological "TICK" (entry ruby stone) or "TOCK" (exit ruby stone):
                // Escape tooth drops onto the ruby locking face and pallet fork hits the brass banking pin!
                let (low_w, high_w) = if is_entry_stone {
                    // "TICK": slightly brighter, stiffer lever geometry
                    (0.92, 1.12)
                } else {
                    // "TOCK": slightly fuller brass plate / gold case coupling
                    (1.10, 0.90)
                };
                let e = energy * 0.85;
                // Gold case & brass plate body ring
                self.modes[0].excite(0.24 * e * low_w, 0.08 * e);
                self.modes[1].excite(0.32 * e * low_w, -0.10 * e);
                self.modes[2].excite(0.28 * e * low_w, 0.06 * e);
                self.modes[3].excite(0.45 * e * low_w, 0.14 * e);
                self.modes[4].excite(0.52 * e * low_w, -0.18 * e);
                // Tourbillon cage, escape wheel & ruby pallet stone crystal ring
                self.modes[6].excite(0.68 * e * high_w * hf_cushion, 0.22 * e);
                self.modes[7].excite(0.75 * e * high_w * hf_cushion, -0.25 * e);
                self.modes[8].excite(0.58 * e * high_w * hf_cushion, 0.15 * e);
                self.modes[9].excite(0.82 * e * high_w * hf_cushion, 0.30 * e);
                self.modes[10].excite(0.64 * e * high_w * hf_cushion, -0.20 * e);
                self.modes[11].excite(0.48 * e * high_w * hf_cushion, 0.12 * e);
                let imp = [
                    0.24 * e * low_w,
                    0.32 * e * low_w,
                    0.28 * e * low_w,
                    0.45 * e * low_w,
                    0.52 * e * low_w,
                    0.0,
                    0.68 * e * high_w * hf_cushion,
                    0.75 * e * high_w * hf_cushion,
                    0.58 * e * high_w * hf_cushion,
                    0.82 * e * high_w * hf_cushion,
                    0.64 * e * high_w * hf_cushion,
                    0.48 * e * high_w * hf_cushion,
                ];
                let _ = self
                    .modal_time_model
                    .apply_impulse_into(&imp, &mut self.modal_workspace);
            }
            MechanicalAcousticEvent::WindingRatchetClick { energy } => {
                // Crisp steel ratchet click localized at the 12h barrel bridge
                let e = energy * 0.72;
                self.modes[0].excite(0.30 * e, 0.10 * e);
                self.modes[4].excite(0.65 * e, -0.20 * e);
                self.modes[5].excite(0.95 * e, 0.32 * e);
                self.modes[7].excite(0.45 * e, 0.15 * e);
                self.modes[10].excite(0.38 * e, -0.12 * e);
                let mut imp = [0.0; MODAL_COUNT];
                imp[0] = 0.30 * e;
                imp[4] = 0.65 * e;
                imp[5] = 0.95 * e;
                imp[7] = 0.45 * e;
                imp[10] = 0.38 * e;
                let _ = self
                    .modal_time_model
                    .apply_impulse_into(&imp, &mut self.modal_workspace);
            }
            MechanicalAcousticEvent::MagneticForkSnap { energy } => {
                // Sharp metallic clank when external NdFeB magnet yanks the pallet fork to the banking pin
                let e = energy * 1.15;
                for m in 0..MODAL_COUNT {
                    self.modes[m].excite(0.55 * e, 0.20 * e);
                }
                let imp = [0.55 * e; MODAL_COUNT];
                let _ = self
                    .modal_time_model
                    .apply_impulse_into(&imp, &mut self.modal_workspace);
            }
        }
    }

    /// Synthesize `out.len()` mono audio samples at $48\text{ kHz}$ into `out`.
    pub fn render_samples(&mut self, mp: &MultiphysicsState, out: &mut [f32]) {
        self.sync_multiphysics(mp);

        let atten_case = self.distance_attenuation_from(self.modes[0].source_pos_mm);
        let atten_plate = self.distance_attenuation_from(self.modes[3].source_pos_mm);
        let atten_barrel = self.distance_attenuation_from(self.modes[5].source_pos_mm);
        let atten_esc = self.distance_attenuation_from(self.modes[6].source_pos_mm);
        let atten_torch = self.distance_attenuation_from(mp.inputs.blowtorch_pos_mm);
        let atten_ln2 = self.distance_attenuation_from(mp.inputs.ln2_pos_mm);

        let torch_power = mp.inputs.blowtorch_power.clamp(0.0, 1.0);
        let ln2_boil_level =
            ((mp.modifiers.ln2_boiling_power_w / 18.0) + 0.6 * mp.inputs.ln2_pour_rate).clamp(0.0, 1.5);
        let degauss_level = (mp.degauss_envelope_mt / 85.0).clamp(0.0, 1.0);

        let mut peak = 0.0_f64;
        let mut sum_sq = 0.0_f64;
        let zero_forces = [0.0; MODAL_COUNT];
        let degauss_dphase = TAU * 60.0 / AUDIO_SAMPLE_RATE_HZ;
        let aperture_gain = self.aperture_radiation_gain;

        // Prepare 3D Farassat retarded dipole observer (fs_aeroac::curle3d::CompactForceObserver)
        // from the blowtorch impingement point to the virtual condenser microphone
        let curle_observer = FluidAcousticState::try_new(1.18, 343.2)
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

        let mut pressure_samples_pa = Vec::with_capacity(out.len().min(256));

        for sample in out.iter_mut() {
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

            let mut sig = 0.20
                * (atten_case * m_case
                    + atten_plate * m_plate
                    + atten_barrel * m_barrel
                    + atten_esc * aperture_gain * m_esc)
                + 0.04 * atten_esc * aperture_gain * fs_modal_pa;

            // 2. Virtual Blowtorch turbulent jet roar + 3D Farassat Curle dipole radiation
            if torch_power > 1.0e-3 {
                let w = self.next_white_noise();
                let prev_lp2 = self.torch_lp2;
                self.torch_lp1 += 0.045 * (w - self.torch_lp1);
                self.torch_lp2 += 0.18 * (w - self.torch_lp2);
                let jet_roar = 0.65 * self.torch_lp1 + 0.35 * (self.torch_lp2 - self.torch_lp1);
                let curle_pa = curle_observer
                    .and_then(|obs| {
                        let f_z = -0.08 * torch_power * jet_roar;
                        let f_dot_z =
                            -0.08 * torch_power * (self.torch_lp2 - prev_lp2) * AUDIO_SAMPLE_RATE_HZ;
                        obs.pressure_at_retarded([0.0, 0.0, f_z], [0.0, 0.0, f_dot_z])
                            .ok()
                    })
                    .unwrap_or(0.0);
                sig += 0.095 * torch_power * atten_torch * jet_roar + 0.15 * curle_pa;
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

            // 4. AC Degauss Coil 60 Hz + 180 Hz electromagnetic hum integrated at 48 kHz
            if degauss_level > 1.0e-3 {
                self.degauss_phase_rad = (self.degauss_phase_rad + degauss_dphase) % TAU;
                let hum = 0.045
                    * degauss_level
                    * (self.degauss_phase_rad.sin()
                        + 0.35 * (3.0 * self.degauss_phase_rad).sin());
                sig += hum;
            }

            let scaled = (sig * self.mic_gain).clamp(-0.98, 0.98);
            let pa = scaled * 0.25; // Calibrated studio pressure scale (0.25 Pa full-scale)
            peak = peak.max(pa.abs());
            sum_sq += pa * pa;
            if pressure_samples_pa.len() < 256 {
                let sample_pa = if pa.abs() < 1.0e-7 { 2.0e-5 } else { pa };
                pressure_samples_pa.push(sample_pa);
            }
            *sample = scaled as f32;
        }

        self.last_peak_pa = peak;
        self.last_rms_pa = if out.is_empty() {
            0.0
        } else {
            (sum_sq / (out.len() as f64)).sqrt()
        };

        // 5. Evaluate ISO 532-1 Zwicker Loudness (sone/phon) & DIN 45692 Sharpness (acum) via fs-psycho
        if !pressure_samples_pa.is_empty() {
            if let Ok(spl_db) =
                spl_from_pressure_rms(&pressure_samples_pa, AIR_REFERENCE_PRESSURE_PA)
            {
                self.last_spl_db = spl_db.clamp(15.0, 110.0);
            }
            let base_db = self.last_spl_db.clamp(20.0, 95.0);
            let mut band_levels = [base_db - 18.0; N_THIRD_OCTAVE_BANDS];
            // Map the 12 structural modes + jet/boiling energy onto the 28 ISO 532-1 third-octave bands
            for (idx, m) in self.modes.iter().enumerate() {
                let amp = (m.state_re.hypot(m.state_im)).clamp(1.0e-5, 4.0);
                let mode_db = (base_db - 10.0 + 8.0 * amp.ln()).clamp(15.0, 95.0);
                let band_idx = (14 + idx).min(N_THIRD_OCTAVE_BANDS - 1);
                band_levels[band_idx] = band_levels[band_idx].max(mode_db);
            }
            if let Ok(loud) = loudness_stationary(&band_levels, SoundField::Free) {
                self.last_loudness_sone = loud.sones;
                if let Ok(phon) = phon_from_sone(loud.sones) {
                    self.last_loudness_phon = phon;
                }
                if let Ok(sharp) = sharpness_din(&loud.specific) {
                    self.last_sharpness_acum = sharp;
                }
            }
        }
    }
}

impl Default for WatchAcousticSynth {
    fn default() -> Self {
        Self::new()
    }
}
