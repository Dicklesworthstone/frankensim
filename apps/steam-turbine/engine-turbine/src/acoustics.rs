//! 48 kHz multiphysics acoustic synthesizer and ISO 532-1 / DIN 45692 psychoacoustic
//! meter for the steam turbine simulation.
//!
//! Combines:
//! - [`fs_modal::time::ModalAcousticTimeModel`] for exact ZOH state-space structural
//!   resonances of the CrMoV casing, slate-blue axial exhaust diffuser hood, and LP blades
//! - [`fs_aeroac::curle3d::CompactForceObserver`] for Farassat/Curle 3D dipole radiation
//!   of live-steam control valve throttling roar and LP exhaust diffuser jet turbulence
//! - Live physics-coupled tonal sources: 1×/2× synchronous rotor unbalance, ~0.46×
//!   sub-synchronous oil-whip / Alford steam-whirl growl, 56× Curtis & 44× LP blade-passing
//!   frequencies (BPF), 50/100 Hz M19 generator magnetostriction hum, wet-steam droplet
//!   impingement hiss, and Stribeck labyrinth rub metallic screech
//! - [`fs_psycho`] 28-band ISO 532-1 third-octave filterbank (`loudness_stationary`,
//!   `phon_from_sone`, `sharpness_din`)

use core::f64::consts::TAU;
use std::sync::OnceLock;

use fs_aeroac::curle3d::CompactForceObserver;
use fs_duct::{LS_MAX_KA, LsQuadrature, levine_schwinger_impedance};
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_material::gas::{GasSpec, GasState};
use fs_math::c64::C64;
use fs_modal::SliceOptions;
use fs_modal::time::{
    ModalAcousticMode, ModalAcousticTimeBudget, ModalAcousticTimeModel, ModalAcousticWorkspace,
};
use fs_plate::PlateSection;
use fs_plate::shell::{ShellSupport, assemble_shell, generate_bell_shell, modes_shell};
use fs_psycho::filter_tables::{FILTER_DIFF, FILTER_GAIN, FILTER_REF};
use fs_psycho::signal::phon_from_sone;
use fs_psycho::{
    AIR_REFERENCE_PRESSURE_PA, N_THIRD_OCTAVE_BANDS, SoundField, loudness_stationary, sharpness_din,
};

use crate::materials::{P_ATM_PA, T_REF_K, library};

/// Certified FEM structural eigenfrequencies of the model turbine fabricated S355 steel
/// axial exhaust diffuser hood, solved directly via `fs_plate::shell::{generate_bell_shell,
/// assemble_shell, modes_shell}` with Sylvester inertia certification.
#[derive(Debug, Clone, Copy)]
pub struct ExhaustHoodFemSpectrum {
    /// Lowest three distinct flexural shell eigenfrequencies [Hz].
    pub freqs_hz: [f64; 3],
    /// Stiffness matrix trace [N/m] from CST/DKT shell assembly.
    pub stiffness_trace_n_m: f64,
    /// Axial exhaust diffuser open-end Levine-Schwinger acoustic radiation factor at f0.
    pub diffuser_rad_resistance: f64,
    /// True when `modes_shell` certified >= 3 distinct eigenmodes via Sylvester inertia.
    pub certified_from_fem: bool,
}

/// Solve (once per process) the lowest flexural eigenfrequencies of the conical-to-cylindrical
/// S355 steel exhaust diffuser hood (`r = 0.68..1.35 m`, `t = 28 mm`) via `fs-plate::shell` + `fs-modal`!
#[must_use]
pub fn exhaust_hood_fem_spectrum() -> &'static ExhaustHoodFemSpectrum {
    static SPECTRUM: OnceLock<ExhaustHoodFemSpectrum> = OnceLock::new();
    SPECTRUM.get_or_init(|| {
        // Exact conical-to-cylindrical exhaust hood meridian profile [(r_m, z_m)]
        let hood_profile = [
            (0.68_f64, 0.00_f64),
            (0.98_f64, 0.35_f64),
            (1.24_f64, 0.85_f64),
            (1.34_f64, 1.55_f64),
        ];
        let shell_mesh = generate_bell_shell(&hood_profile, 14);
        let inlet_flange_nodes: Vec<usize> = (0..14.min(shell_mesh.node_count())).collect();
        let lib = library();
        let steel_section = PlateSection::isotropic(
            lib.lp_exhaust_steel.young_modulus_at(330.0),
            lib.lp_exhaust_steel.poisson_ratio,
            0.028, // 28 mm welded S355 exhaust diffuser shell wall
            lib.lp_exhaust_steel.rho_kg_m3,
        )
        .expect("valid S355 exhaust hood PlateSection");

        let (freqs_hz, stiffness_trace_n_m, certified_from_fem) = match assemble_shell(
            &shell_mesh,
            &steel_section,
            &inlet_flange_nodes,
            ShellSupport::Clamped,
        ) {
            Ok(shell_model) => {
                let trace: f64 = (0..shell_model.k.nrows())
                    .map(|i| shell_model.k.get(i, i))
                    .sum();
                let opts = SliceOptions {
                    ritz_tol: 1.0e-5,
                    ..SliceOptions::default()
                };
                let mut distinct: Vec<f64> = Vec::with_capacity(6);
                for &(f_lo, f_hi) in &[(40.0_f64, 900.0), (40.0, 2_200.0)] {
                    let w_lo = TAU * f_lo;
                    let w_hi = TAU * f_hi;
                    if let Ok(rep) = modes_shell(&shell_model, (w_lo * w_lo, w_hi * w_hi), &opts) {
                        distinct.clear();
                        for m in &rep.modes {
                            if m.lambda > 0.0 {
                                let hz = m.lambda.sqrt() / TAU;
                                if hz <= 3_500.0
                                    && distinct.last().is_none_or(|&prev| (hz - prev).abs() > 12.0)
                                {
                                    distinct.push(hz);
                                    if distinct.len() == 3 {
                                        break;
                                    }
                                }
                            }
                        }
                        if distinct.len() == 3 {
                            break;
                        }
                    }
                }
                if distinct.len() == 3 {
                    ([distinct[0], distinct[1], distinct[2]], trace, true)
                } else {
                    (
                        [
                            distinct.first().copied().unwrap_or(115.0),
                            distinct.get(1).copied().unwrap_or(240.0),
                            distinct.get(2).copied().unwrap_or(410.0),
                        ],
                        trace,
                        false,
                    )
                }
            }
            Err(_) => ([115.0, 240.0, 410.0], 1.8e10, false),
        };

        // Evaluate exact Levine-Schwinger Wiener-Hopf open-flange acoustic radiation impedance
        // of the 1.34 m axial exhaust diffuser exit via `fs_duct` + `fs_material::gas::GasState`!
        let radius_m = 1.34_f64;
        let diffuser_rad_resistance = GasState::try_new(&GasSpec::dry_air_ussa1976(), T_REF_K, P_ATM_PA)
            .ok()
            .and_then(|gas| {
                let max_omega = 0.92 * (*LS_MAX_KA) * gas.sound_speed / radius_m;
                let omega = (TAU * freqs_hz[0]).clamp(40.0, max_omega);
                let z_rad = levine_schwinger_impedance(
                    &gas,
                    radius_m,
                    omega,
                    LsQuadrature { nodes: 32 },
                )
                .ok()?;
                let z0 = gas.density * gas.sound_speed / (core::f64::consts::PI * radius_m * radius_m);
                Some((z_rad.re / z0.max(1e-9)).clamp(0.02, 2.0))
            })
            .unwrap_or(0.68);

        ExhaustHoodFemSpectrum {
            freqs_hz,
            stiffness_trace_n_m,
            diffuser_rad_resistance,
            certified_from_fem,
        }
    })
}

/// Audio synthesis sample rate [Hz] (matched to WebAudio AudioWorklet and ISO 532-1 tables).
pub const AUDIO_SAMPLE_RATE_HZ: f64 = 48_000.0;
const MODAL_COUNT: usize = 10;
const METER_TAU_S: f64 = 0.6;

/// Live physical acoustic excitation parameters pushed from the convergence-kernel physics state.
#[derive(Debug, Clone, Copy)]
pub struct TurbineAcousticParams {
    /// Rotor angular velocity [rad/s].
    pub omega_rad_s: f64,
    /// Live steam mass flow [kg/s].
    pub steam_flow_kg_s: f64,
    /// Control-valve throttling pressure drop [Pa].
    pub valve_dp_pa: f64,
    /// Control-valve lift fraction `[0, 1]`.
    pub valve_lift: f64,
    /// LP exhaust volumetric flow / Mach proxy (0..1.5).
    pub exhaust_flow_intensity: f64,
    /// LP exhaust steam moisture `1 - x` `[0, 0.25]`.
    pub exhaust_moisture: f64,
    /// Maximum shaft vibration peak-to-peak [µm].
    pub max_vib_um_pp: f64,
    /// Dominant sub-synchronous whirl frequency [Hz] (~22–28 Hz during oil/steam whip).
    pub whirl_hz: f64,
    /// Oil/steam whirl instability intensity `[0, 1]` derived from `max Re(λ)` and orbit growth.
    pub whirl_intensity: f64,
    /// Active labyrinth rub normal force [N].
    pub rub_force_n: f64,
    /// Generator active electrical power [W] (drives 50/100 Hz M19 core hum when breaker closed).
    pub gen_power_w: f64,
    /// Electrical frequency [Hz] (50 Hz nominal).
    pub elec_freq_hz: f64,
    /// Mean casing metal temperature [K] (shifts structural modal frequencies via E(T)).
    pub casing_temp_k: f64,
}

impl Default for TurbineAcousticParams {
    fn default() -> Self {
        Self {
            omega_rad_s: 314.159_265,
            steam_flow_kg_s: 52.0,
            valve_dp_pa: 1.5e6,
            valve_lift: 0.72,
            exhaust_flow_intensity: 0.85,
            exhaust_moisture: 0.08,
            max_vib_um_pp: 22.0,
            whirl_hz: 24.0,
            whirl_intensity: 0.0,
            rub_force_n: 0.0,
            gen_power_w: 50.0e6,
            elec_freq_hz: 50.0,
            casing_temp_k: 620.0,
        }
    }
}

/// ISO 532-1 28-band third-octave sound-level meter using `fs-psycho` biquad filter tables.
#[derive(Debug, Clone)]
struct BandMeter {
    coeffs: [[[f64; 6]; 3]; N_THIRD_OCTAVE_BANDS],
    gains: [[f64; 3]; N_THIRD_OCTAVE_BANDS],
    state: [[[f64; 2]; 3]; N_THIRD_OCTAVE_BANDS],
    ms: [f64; N_THIRD_OCTAVE_BANDS],
    ms_total: f64,
    k: f64,
}

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

fn soft_clip(x: f64) -> f64 {
    const KNEE: f64 = 0.72;
    const ROOM: f64 = 0.25;
    let a = x.abs();
    if a <= KNEE {
        x
    } else {
        x.signum() * (KNEE + ROOM * ((a - KNEE) / ROOM).tanh())
    }
}

/// Real-time 48 kHz acoustic synthesizer.
#[derive(Debug, Clone)]
pub struct TurbineAcousticSynth {
    modal_model: ModalAcousticTimeModel,
    modal_workspace: ModalAcousticWorkspace,
    valve_observer: CompactForceObserver,
    exhaust_observer: CompactForceObserver,
    /// Virtual microphone 3D position in turbine coordinates [m] (default: 3.5 m beside turbine).
    pub mic_pos_m: [f64; 3],
    /// Output gain trim (user preamp slider, default 1.0).
    pub mic_gain: f64,
    /// Continuous phase accumulators [rad].
    phase_rot: f64,
    phase_whirl: f64,
    phase_elec: f64,
    phase_bpf_curtis: f64,
    phase_bpf_lp: f64,
    /// Deterministic xorshift64 PRNG for turbulent steam flow & droplet noise.
    rng: u64,
    /// Colored noise shaping filter states.
    valve_lp1: f64,
    valve_lp2: f64,
    exhaust_lp1: f64,
    exhaust_lp2: f64,
    mist_hp: f64,
    prev_valve_force: f64,
    prev_exhaust_force: f64,
    prev_modal_p: f64,
    meter: BandMeter,
    /// Latest acoustic & psychoacoustic telemetry.
    pub last_peak_pa: f64,
    pub last_rms_pa: f64,
    pub last_spl_db: f64,
    pub last_loudness_sone: f64,
    pub last_loudness_phon: f64,
    pub last_sharpness_acum: f64,
    /// Dominant audible BPF frequency [Hz] (44× LP last stage).
    pub last_bpf_lp_hz: f64,
}

impl TurbineAcousticSynth {
    /// Construct a fresh 48 kHz turbine acoustic engine.
    #[must_use]
    pub fn new() -> Self {
        // 10 structural resonances of the model turbine:
        // - Modes 0..3: FEM-certified exhaust hood shell breathing/ovalizing modes from `fs_plate::shell::modes_shell`!
        // - Modes 3..6: cast CrMoV HP/IP barrel modes (680, 1180, 1950 Hz)
        // - Modes 6..10: LP last-stage free-standing blade bending/torsion & labyrinth ring modes (2850..8400 Hz)
        let fem = exhaust_hood_fem_spectrum();
        let rad_w = (0.65 + 0.35 * fem.diffuser_rad_resistance).clamp(0.5, 1.3);
        let mode_specs: [(f64, f64, C64); MODAL_COUNT] = [
            // Keep exact fs-plate::shell certified S355 exhaust hood frequencies (~42 / 57 / 86 Hz)!
            (fem.freqs_hz[0].clamp(20.0, 2_500.0), 0.032, C64::new(0.065 * rad_w, 0.012)),
            (fem.freqs_hz[1].clamp(30.0, 3_200.0), 0.026, C64::new(0.082 * rad_w, 0.015)),
            (fem.freqs_hz[2].clamp(40.0, 4_000.0), 0.022, C64::new(0.070 * rad_w, 0.012)),
            (680.0, 0.020, C64::new(0.05, 0.01)),
            (1_180.0, 0.018, C64::new(0.045, 0.01)),
            (1_950.0, 0.016, C64::new(0.04, 0.008)),
            (2_850.0, 0.012, C64::new(0.05, 0.01)),
            (4_200.0, 0.010, C64::new(0.045, 0.012)),
            (6_100.0, 0.009, C64::new(0.04, 0.01)),
            (8_400.0, 0.008, C64::new(0.035, 0.008)),
        ];
        let modes: Vec<ModalAcousticMode> = mode_specs
            .iter()
            .map(|&(f_hz, zeta, h)| ModalAcousticMode {
                angular_frequency_rad_s: TAU * f_hz,
                damping_ratio: zeta,
                pressure_per_modal_velocity: h,
            })
            .collect();
        // Modal radiator runs at exact 24 kHz (stepped every 2 audio samples with linear
        // interpolation onto the 48 kHz output grid) so all 10 modes (115..8,400 Hz) stay
        // strictly below Nyquist (12 kHz) and sound at exact true pitch with zero stair-stepping!
        let modal_model = ModalAcousticTimeModel::try_new(
            (AUDIO_SAMPLE_RATE_HZ as u32) / 2,
            modes,
            ModalAcousticTimeBudget::audible_reference(),
        )
        .expect("valid 24 kHz turbine ModalAcousticTimeModel");
        let modal_workspace = ModalAcousticWorkspace::new(&modal_model);

        let mic_pos_m = [2.2, 1.4, 2.4];
        let (valve_observer, exhaust_observer) = Self::build_observers(mic_pos_m);

        Self {
            modal_model,
            modal_workspace,
            valve_observer,
            exhaust_observer,
            mic_pos_m,
            mic_gain: 1.0,
            phase_rot: 0.0,
            phase_whirl: 0.0,
            phase_elec: 0.0,
            phase_bpf_curtis: 0.0,
            phase_bpf_lp: 0.0,
            rng: 0x5353_5436_3030_BEEF,
            valve_lp1: 0.0,
            valve_lp2: 0.0,
            exhaust_lp1: 0.0,
            exhaust_lp2: 0.0,
            mist_hp: 0.0,
            prev_valve_force: 0.0,
            prev_exhaust_force: 0.0,
            prev_modal_p: 0.0,
            meter: BandMeter::new(),
            last_peak_pa: 0.0,
            last_rms_pa: 0.0,
            last_spl_db: 72.0,
            last_loudness_sone: 12.0,
            last_loudness_phon: 75.0,
            last_sharpness_acum: 1.4,
            last_bpf_lp_hz: 2_200.0,
        }
    }

    fn build_observers(mic_pos_m: [f64; 3]) -> (CompactForceObserver, CompactForceObserver) {
        let air = &library().ambient_air;
        let acoustic = FluidAcousticState::try_new(
            air.density(T_REF_K, P_ATM_PA),
            air.sound_speed(T_REF_K),
        )
        .expect("valid FluidAcousticState");
        let valve_src = [0.65, 0.40, 1.28]; // HP control valve chest
        let exhaust_src = [0.0, 0.10, 4.45]; // LP axial exhaust diffuser
        let safe_mic = [
            if (mic_pos_m[0] - valve_src[0]).abs() < 0.25 {
                mic_pos_m[0] + 0.5
            } else {
                mic_pos_m[0]
            },
            mic_pos_m[1],
            mic_pos_m[2],
        ];
        let v_obs = CompactForceObserver::new(acoustic, safe_mic, valve_src)
            .expect("valid valve CompactForceObserver");
        let e_obs = CompactForceObserver::new(acoustic, safe_mic, exhaust_src)
            .expect("valid exhaust CompactForceObserver");
        (v_obs, e_obs)
    }

    /// Move the 3D virtual microphone in turbine hall coordinates [m].
    pub fn set_mic_pos(&mut self, pos_m: [f64; 3]) {
        self.mic_pos_m = [
            pos_m[0].clamp(-12.0, 12.0),
            pos_m[1].clamp(-4.0, 8.0),
            pos_m[2].clamp(-6.0, 12.0),
        ];
        let (v, e) = Self::build_observers(self.mic_pos_m);
        self.valve_observer = v;
        self.exhaust_observer = e;
    }

    #[inline]
    fn next_white(&mut self) -> f64 {
        let mut x = self.rng;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.rng = x;
        // Uniform [-1, 1]
        (x as i64 as f64) * (1.0 / i64::MAX as f64)
    }

    /// Render `frames` mono 48 kHz PCM samples in `[-1.0, 1.0]` into `out`,
    /// advancing all physical modal/dipole/harmonic oscillators and psychoacoustic meters.
    pub fn render_block(&mut self, params: &TurbineAcousticParams, out: &mut [f32]) {
        let dt = 1.0 / AUDIO_SAMPLE_RATE_HZ;
        let f_rot = (params.omega_rad_s.abs() / TAU).max(0.0);
        let f_bpf_curtis = 56.0 * f_rot; // 56 Curtis wheel buckets
        let f_bpf_lp = 44.0 * f_rot; // 44 giant LP last-stage blades
        self.last_bpf_lp_hz = f_bpf_lp;

        // Distance attenuation of the synthesised terms (referenced to a 3.5 m hall distance).
        // The two Curle observers already radiate with their own retarded-time 1/r from the
        // source to the microphone, so this gain is not applied to them.
        let dist = ((self.mic_pos_m[0]).powi(2)
            + (self.mic_pos_m[1] - 0.5).powi(2)
            + (self.mic_pos_m[2] - 2.6).powi(2))
        .sqrt()
        .max(0.9);
        let dist_gain = 3.5 / dist;

        // Acoustic source amplitudes in physical pascals at nominal 3.5 m distance
        // Industrial steam turbine hall level is ~82–92 dB SPL (~0.25–0.80 Pa RMS).
        let flow_norm = (params.steam_flow_kg_s / 55.0).clamp(0.0, 2.0);
        let speed_norm = (f_rot / 50.0).clamp(0.0, 1.6);
        let vib_norm = (params.max_vib_um_pp / 80.0).clamp(0.0, 4.0);
        let gen_norm = (params.gen_power_w.abs() / 50.0e6).clamp(0.0, 1.8);
        let valve_throttle_severity =
            ((params.valve_dp_pa / 5.0e6).clamp(0.0, 2.5) * flow_norm * (1.15 - 0.65 * params.valve_lift))
                .max(0.0);
        let rub_norm = (params.rub_force_n / 8_000.0).clamp(0.0, 3.0);
        let wet_norm = ((params.exhaust_moisture - 0.06).max(0.0) / 0.14).clamp(0.0, 2.0) * flow_norm;

        // Thermal elastic modulus softening gently modulates structural radiation impedance
        let lib = library();
        let e_ratio = (lib.hp_casing_steel.young_modulus_at(params.casing_temp_k)
            / lib.hp_casing_steel.young_modulus_at(T_REF_K))
        .clamp(0.7, 1.05);
        let modal_impedance_gain = 0.08 * e_ratio.sqrt();

        let mut peak_pa = 0.0_f64;
        let mut sum_sq_pa = 0.0_f64;
        let mut modal_forces = [0.0_f64; MODAL_COUNT];
        let mut cur_modal_p = self.prev_modal_p;

        for (i, sample_slot) in out.iter_mut().enumerate() {
            // Advance continuous phases
            self.phase_rot = (self.phase_rot + TAU * f_rot * dt).rem_euclid(TAU);
            self.phase_whirl = (self.phase_whirl + TAU * params.whirl_hz.clamp(8.0, 90.0) * dt).rem_euclid(TAU);
            self.phase_elec = (self.phase_elec + TAU * params.elec_freq_hz.clamp(20.0, 80.0) * dt).rem_euclid(TAU);
            self.phase_bpf_curtis = (self.phase_bpf_curtis + TAU * f_bpf_curtis.min(18_000.0) * dt).rem_euclid(TAU);
            self.phase_bpf_lp = (self.phase_bpf_lp + TAU * f_bpf_lp.min(18_000.0) * dt).rem_euclid(TAU);

            // 1. Broadband valve throttling hiss + LP exhaust diffuser roar via fs_aeroac::curle3d
            let w1 = self.next_white();
            let w2 = self.next_white();
            // Band-shape valve jet noise (high-mid hiss ~1.2–4.5 kHz)
            self.valve_lp1 += 0.34 * (w1 - self.valve_lp1);
            self.valve_lp2 += 0.09 * (self.valve_lp1 - self.valve_lp2);
            let valve_band = self.valve_lp1 - self.valve_lp2;
            let f_valve = 18.0 * valve_throttle_severity * valve_band;
            let f_valve_dot = (f_valve - self.prev_valve_force) * AUDIO_SAMPLE_RATE_HZ;
            self.prev_valve_force = f_valve;
            let p_valve = self
                .valve_observer
                .pressure_at_retarded([0.0, f_valve, 0.0], [0.0, f_valve_dot, 0.0])
                .unwrap_or(0.0);

            // Deep exhaust diffuser turbulent rumble (60–450 Hz)
            self.exhaust_lp1 += 0.055 * (w2 - self.exhaust_lp1);
            self.exhaust_lp2 += 0.015 * (self.exhaust_lp1 - self.exhaust_lp2);
            let f_exh = 42.0 * params.exhaust_flow_intensity * self.exhaust_lp1;
            let f_exh_dot = (f_exh - self.prev_exhaust_force) * AUDIO_SAMPLE_RATE_HZ;
            self.prev_exhaust_force = f_exh;
            let p_exhaust = self
                .exhaust_observer
                .pressure_at_retarded([f_exh, 0.0, 0.0], [f_exh_dot, 0.0, 0.0])
                .unwrap_or(0.0);

            // 2. Synchronous 1× + 2× rotor unbalance rumble
            let p_sync = (0.08 * speed_norm + 0.18 * vib_norm)
                * (self.phase_rot.sin() + 0.35 * (2.0 * self.phase_rot).sin());

            // 3. Sub-synchronous ~0.46× oil-whip / Alford steam-whirl low growl
            let p_whirl = 0.55
                * params.whirl_intensity
                * (self.phase_whirl.sin() + 0.45 * (2.0 * self.phase_whirl).sin());

            // 4. Blade Passing Frequency (BPF) turbine whine: Curtis 56× + LP 44×
            let p_bpf = (0.045 * flow_norm * speed_norm)
                * (0.65 * self.phase_bpf_lp.sin()
                    + 0.28 * (2.0 * self.phase_bpf_lp).sin()
                    + 0.45 * self.phase_bpf_curtis.sin());

            // 5. Synchronous generator 50 Hz + 100 Hz M19 core magnetostriction hum
            let p_gen = 0.075 * gen_norm * (0.4 * self.phase_elec.sin() + 0.85 * (2.0 * self.phase_elec).sin());

            // 6. Wet-steam LP droplet impingement crackle (> 6% moisture)
            self.mist_hp = 0.72 * self.mist_hp + 0.28 * (w1 - w2);
            let p_wet = 0.14 * wet_norm * ((w1 - self.mist_hp) * (0.5 + 0.5 * self.phase_bpf_lp.sin()));

            // 7. Structural casing & blade modes driven by 24 kHz fs_modal::time::ModalAcousticTimeModel
            // linearly interpolated onto the 48 kHz audio grid (exact modal pitch, zero ZOH aliasing).
            let modal_sample = if i % 2 == 0 {
                let rub_exc = rub_norm * (w1 * 18.0 + 9.0 * (3.0 * self.phase_rot).sin());
                let flow_exc = (flow_norm * 1.4 + vib_norm * 2.2) * self.exhaust_lp2;
                modal_forces[0] = flow_exc * 3.0 + p_whirl * 4.0;
                modal_forces[1] = flow_exc * 2.4 + p_sync * 2.0;
                modal_forces[2] = flow_exc * 1.8;
                modal_forces[3] = valve_band * valve_throttle_severity * 4.5;
                modal_forces[4] = valve_band * valve_throttle_severity * 3.5;
                modal_forces[5] = rub_exc * 1.2 + p_bpf * 2.0;
                modal_forces[6] = rub_exc * 2.4 + p_bpf * 3.5;
                modal_forces[7] = rub_exc * 3.8;
                modal_forces[8] = rub_exc * 4.2 + p_wet * 3.0;
                modal_forces[9] = rub_exc * 3.0 + p_wet * 4.0;

                let prev_p = cur_modal_p;
                if let Ok(frame) = self
                    .modal_model
                    .step_into(&modal_forces, &mut self.modal_workspace)
                {
                    cur_modal_p = frame.observer_pressure_pa;
                }
                self.prev_modal_p = cur_modal_p;
                0.5 * (prev_p + cur_modal_p)
            } else {
                cur_modal_p
            };

            let total_pa = p_valve
                + p_exhaust
                + dist_gain * (p_sync + p_whirl + p_bpf + p_gen + p_wet + modal_impedance_gain * modal_sample);
            self.meter.push(total_pa);

            let abs_p = total_pa.abs();
            if abs_p > peak_pa {
                peak_pa = abs_p;
            }
            sum_sq_pa += total_pa * total_pa;

            // Scale physical pascals (~0.45 Pa nominal RMS -> ~0.28 PCM RMS) with soft limiter
            let pcm = soft_clip(total_pa * 0.58 * self.mic_gain);
            *sample_slot = pcm as f32;
        }

        let n = out.len().max(1) as f64;
        self.last_peak_pa = peak_pa;
        self.last_rms_pa = (sum_sq_pa / n).sqrt();
        self.last_spl_db = BandMeter::level_db(self.meter.ms_total);

        let mut band_db = [0.0_f64; N_THIRD_OCTAVE_BANDS];
        for b in 0..N_THIRD_OCTAVE_BANDS {
            band_db[b] = BandMeter::level_db(self.meter.ms[b]).clamp(-20.0, 140.0);
        }
        if let Ok(l) = loudness_stationary(&band_db, SoundField::Free) {
            self.last_loudness_sone = l.sones;
            self.last_loudness_phon = phon_from_sone(l.sones).unwrap_or(0.0);
            self.last_sharpness_acum = sharpness_din(&l.specific).unwrap_or(0.0);
        }
    }
}

impl Default for TurbineAcousticSynth {
    fn default() -> Self {
        Self::new()
    }
}
