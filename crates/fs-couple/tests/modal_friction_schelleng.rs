//! Known-answer test of the single-port friction route (`FrictionModalSystem`,
//! the route every one-port bowed card runs) inside Schelleng's window
//! (IRP pexfp.299).
//!
//! **Known answer.** J. C. Schelleng, "The bowed string and the player",
//! J. Acoust. Soc. Am. 53, 26-41 (1973), as restated by Woodhouse on Euphonics
//! section 9.3: "Helmholtz motion is only possible within the wedge-shaped
//! region between the two lines". Below the minimum force the motion "is the
//! double-slipping motion". The bounds, from Euphonics section 9.3.1:
//! - Eq. (4): `F_max = 2 Z0 v_b / (beta (mu_s - mu_d))`.
//! - Eq. (9): `F_min = Z0^2 v_b / (2 R beta^2 (mu_s - mu_d))`.
//! - Eq. (10): `F_max / F_min = 4 beta R / Z0`.
//!
//! There the bridge is a dashpot of impedance `R`, "assumed to be much greater
//! than the impedance of the string" `Z0 = sqrt(T m)`.
//!
//! **String.** An ideal string (no bending), fixed at both ends. Its modes are
//! mass-normalized, `sqrt(2/(m L)) sin(k pi x/L)`, with the bridge at `x = 0`
//! and the bow at `x = beta L`. The dashpot enters as modal damping. It reflects
//! velocity waves with `r = (Z0 - R)/(Z0 + R)`, so every mode decays at
//! `sigma = (c/2L) ln((R + Z0)/(R - Z0))`, and `zeta_k = sigma/omega_k`. These
//! are the dashpot string's exact decay rates, but the bridge point itself does
//! not move. This is the construction of `bowed_stick_slip_fixture.rs`
//! (pexfp.152), where a genuine dashpot connection could not be built on this
//! route. That fixture bows through the `bowed_string` runtime's own load line.
//! This file bows the same string through `FrictionModalSystem::step`.
//!
//! **Start.** Schelleng's bounds are bounds on sustaining Helmholtz motion, so
//! each run starts from ideal Helmholtz motion at mid-stick: `q_k = 0` and
//! `v_k = sqrt(m L/2) (Y omega_1/k) (-1)^(k+1)`, with
//! `Y = 2 v_b/(beta pi omega_1)`. This is the Fourier series of the duty-beta
//! bow-point velocity.
//!
//! One printed case starts instead from rest (`q_k = v_k = 0`), with the bow
//! at full speed and force from the first sample: a step start. Schelleng's
//! wedge is a steady-state result, and no attack-transient answer for this
//! schedule is cited here. So the case is printed, not gated, together with
//! `settled_by_period`: the period from which every later onset interval is
//! within 2 % of `1/f0`.
//!
//! **Pass.** In the last 20 of 60 periods:
//! - exactly one stuck-to-sliding onset per period;
//! - every onset interval within 2 % of the string's period `1/f0`;
//! - the fundamental is the strongest of the first eight harmonics of the
//!   bridge force `T u_x(0, t)`.
//!
//! **Forces.**
//! - The pass is gated at the window's geometric mean, `sqrt(F_min F_max)`.
//!   Eq. (10) gives a ratio of 16 here, so that force is `4 F_min`.
//! - On this equal-decay string the minimum measured in pexfp.152 lies near
//!   twice Eq. (9). `2 F_min` is therefore printed, not gated.
//! - `sqrt(F_min F_max)` from rest is printed, not gated (see **Start**).
//! - The planted negative is `F_min/4`. It must not give Helmholtz motion.
//! - Every case also prints a diagnostic line, not gated: each sliding run that
//!   starts in the tail, with its length in samples and its peak
//!   `|slip speed|`, and the number per period of runs shorter than 5 samples.
//!
//! **No-Claim.** A pass shows that this route sustains Helmholtz motion inside
//! Schelleng's wedge for a resistive termination. It does not show that the
//! route reproduces `F_min`. Nor does it show anything about a card whose
//! termination is spring-dominated, such as the violin E5 card.

use fs_couple::modal_acoustic_time::{
    ModalAcousticMode, ModalAcousticState, ModalAcousticTimeBudget, ModalAcousticTimeModel,
};
use fs_couple::render::schedule::force::coupled::contact::friction::{
    FrictionDrive, FrictionModalSystem, ModalFrictionConfig, ModalFrictionPort, ModalFrictionRegime,
};
use fs_couple::render::schedule::force::coupled::{
    CoupledModalSystem, ModalAttachment, ModalCouplingConfig,
};
use fs_exec::CancelGate;
use fs_math::c64::C64;
use fs_math::det;
use fs_tribo::FrictionLaw;

// The ideal string of `bowed_stick_slip_fixture.rs`. These are the cello D
// string of the Galluzzo bowing rig (Euphonics 9.6.3 Table 2), used here as
// fixture values.
const F0_HZ: f64 = 146.8;
const TENSION_N: f64 = 111.0;
const DENSITY_KG_M: f64 = 2.7e-3;
/// Bridge resistance over string impedance: `R >> Z0`, as Eq. (9) assumes.
const R_OVER_Z0: f64 = 20.0;
const BOW_SPEED_M_S: f64 = 0.1;
/// Static and asymptotic kinetic limits of Smith & Woodhouse (2000) Eq. (2).
const MU_S: f64 = 1.2;
const MU_D: f64 = 0.35;
/// Stribeck width, 8x below the Helmholtz slip speed `v_b (1/beta - 1)`
/// (0.4 m/s), so the sliding coefficient is `mu_d` to within
/// `(mu_s - mu_d) e^-8`.
const DECAY_M_S: f64 = 0.05;
const BETA: f64 = 0.2;
const RATE_HZ: u32 = 48_000;
/// 64 modes at 48 kHz: the fixture's first rung.
const MODES: usize = 64;
const RUN_PERIODS: f64 = 60.0;
const TAIL_PERIODS: usize = 20;
const INTERVAL_TOLERANCE: f64 = 0.02;
const HARMONICS: usize = 8;
/// Sliding runs shorter than this are counted as toggles in the printed
/// diagnostic. A Helmholtz slip lasts `beta / f0`, about 65 samples here.
const SHORT_RUN_SAMPLES: usize = 5;

fn z0() -> f64 {
    det::sqrt(TENSION_N * DENSITY_KG_M)
}
fn wave_speed() -> f64 {
    det::sqrt(TENSION_N / DENSITY_KG_M)
}
fn length_m() -> f64 {
    wave_speed() / (2.0 * F0_HZ)
}
fn resistance() -> f64 {
    R_OVER_Z0 * z0()
}
/// Euphonics 9.3.1 Eq. (4).
fn schelleng_fmax() -> f64 {
    2.0 * z0() * BOW_SPEED_M_S / (BETA * (MU_S - MU_D))
}
/// Euphonics 9.3.1 Eq. (9).
fn schelleng_fmin() -> f64 {
    z0() * z0() * BOW_SPEED_M_S / (2.0 * resistance() * BETA * BETA * (MU_S - MU_D))
}
#[allow(clippy::cast_precision_loss)]
fn omega(k: usize) -> f64 {
    k as f64 * core::f64::consts::PI * wave_speed() / length_m()
}
/// Mass-normalized mode-shape amplitude `sqrt(2/(m L))` [1/sqrt(kg)].
fn shape_scale() -> f64 {
    det::sqrt(2.0 / (DENSITY_KG_M * length_m()))
}

/// Where a run starts.
#[derive(Clone, Copy)]
enum Start {
    /// Ideal Helmholtz motion at mid-stick.
    Helmholtz,
    /// The string at rest.
    Rest,
}

/// The dashpot string as one modal component, at `start`.
fn dashpot_string(start: Start) -> ModalAcousticTimeModel {
    let sigma = F0_HZ * det::ln((R_OVER_Z0 + 1.0) / (R_OVER_Z0 - 1.0));
    let modes = (1..=MODES)
        .map(|k| ModalAcousticMode {
            angular_frequency_rad_s: omega(k),
            damping_ratio: sigma / omega(k),
            pressure_per_modal_velocity: C64::new(0.0, 0.0),
        })
        .collect();
    let mut model = ModalAcousticTimeModel::try_new(
        RATE_HZ,
        modes,
        ModalAcousticTimeBudget::audible_reference(),
    )
    .unwrap();
    let amplitude = 2.0 * BOW_SPEED_M_S / (BETA * core::f64::consts::PI * omega(1));
    let scale = det::sqrt(DENSITY_KG_M * length_m() / 2.0);
    let initial: Vec<ModalAcousticState> = match start {
        Start::Helmholtz => (1..=MODES)
            .map(|k| {
                #[allow(clippy::cast_precision_loss)]
                let kk = k as f64;
                let sign = if k % 2 == 1 { 1.0 } else { -1.0 };
                ModalAcousticState {
                    displacement_m_sqrt_kg: 0.0,
                    velocity_m_sqrt_kg_per_s: scale * amplitude * omega(1) / kk * sign,
                }
            })
            .collect(),
        Start::Rest => vec![ModalAcousticState::default(); MODES],
    };
    model.restore_states(&initial).unwrap();
    model
}

fn bowed_string(start: Start) -> FrictionModalSystem {
    let gate = CancelGate::new();
    let budgets = ModalCouplingConfig {
        max_modes: MODES,
        max_connections: 4,
        max_setup_terms: 1 << 16,
        nyquist_guard_fraction: 0.9,
        maximum_total_energy_j: 100.0,
        maximum_abs_pressure_pa: 1e6,
        maximum_abs_connection_force_n: 1e6,
        solve_relative_tolerance: 1e-11,
        energy_absolute_tolerance_j: 1e-12,
        energy_relative_tolerance: 1e-8,
    };
    let network =
        CoupledModalSystem::new(vec![dashpot_string(start)], vec![], budgets, &gate).unwrap();
    let bow = ModalAttachment {
        component: 0,
        shapes: (1..=MODES)
            .map(|k| {
                #[allow(clippy::cast_precision_loss)]
                let kk = k as f64;
                shape_scale() * det::sin(kk * core::f64::consts::PI * BETA)
            })
            .collect(),
    };
    let port = ModalFrictionPort {
        left: bow,
        right: None,
        law: FrictionLaw::Stribeck {
            static_mu: MU_S,
            kinetic_mu: MU_D,
            characteristic_speed: DECAY_M_S,
            viscous_per_speed: 0.0,
        },
    };
    let config = ModalFrictionConfig {
        max_iterations: 80,
        maximum_force_n: 100.0,
        maximum_slip_speed_m_s: 10.0,
        force_absolute_tolerance_n: 1e-11,
        force_relative_tolerance: 1e-10,
        sticking_tolerance_m: 1e-12,
    };
    FrictionModalSystem::new(network, port, config, &gate).unwrap()
}

/// Bridge force `T u_x(0, t)` from the mass-normalized displacements.
fn bridge_force(states: &[ModalAcousticState]) -> f64 {
    let slope = shape_scale() * core::f64::consts::PI / length_m();
    states
        .iter()
        .enumerate()
        .map(|(i, s)| {
            #[allow(clippy::cast_precision_loss)]
            let k = (i + 1) as f64;
            TENSION_N * slope * k * s.displacement_m_sqrt_kg
        })
        .sum()
}

struct Verdict {
    onsets: usize,
    worst_interval_error: f64,
    mean_period_s: f64,
    /// Harmonic magnitudes of the bridge force, fundamental first.
    harmonics: Vec<f64>,
    sliding_fraction: f64,
    /// Over the whole run, the period from which every later onset interval
    /// is within `INTERVAL_TOLERANCE` of `1/f0` (`None` if the last is not).
    settled_by_period: Option<usize>,
    /// Each sliding run that starts in the tail: its length in samples and its
    /// peak `|slip speed|` [m/s].
    runs: Vec<(usize, f64)>,
}

impl Verdict {
    fn helmholtz(&self) -> bool {
        let strongest = self.harmonics[1..].iter().all(|&h| self.harmonics[0] > h);
        (self.onsets == TAIL_PERIODS || self.onsets + 1 == TAIL_PERIODS)
            && self.worst_interval_error <= INTERVAL_TOLERANCE
            && strongest
    }
    fn json(&self, label: &str, force: f64) -> String {
        let levels: Vec<String> = self
            .harmonics
            .iter()
            .map(|h| {
                format!(
                    "{:.2}",
                    20.0 * det::ln(h / self.harmonics[0]) / det::ln(10.0)
                )
            })
            .collect();
        let settled = self
            .settled_by_period
            .map_or_else(|| "null".to_string(), |p| p.to_string());
        format!(
            "{{\"test\":\"modal_friction_schelleng\",\"case\":\"{label}\",\"force_n\":{force},\"fmin_eq9_n\":{},\"fmax_eq4_n\":{},\"onsets\":{},\"tail_periods\":{TAIL_PERIODS},\"worst_interval_error\":{},\"mean_period_s\":{},\"string_period_s\":{},\"harmonics_db_re_f1\":[{}],\"sliding_fraction\":{},\"settled_by_period\":{settled},\"helmholtz\":{}}}",
            schelleng_fmin(),
            schelleng_fmax(),
            self.onsets,
            self.worst_interval_error,
            self.mean_period_s,
            1.0 / F0_HZ,
            levels.join(","),
            self.sliding_fraction,
            self.helmholtz()
        )
    }
    /// Printed diagnostic only: the tail's sliding runs and the short ones per period.
    #[allow(clippy::cast_precision_loss)]
    fn runs_json(&self, label: &str) -> String {
        let short = self.runs.iter().filter(|r| r.0 < SHORT_RUN_SAMPLES).count();
        let runs: Vec<String> = self
            .runs
            .iter()
            .map(|(len, peak)| format!("[{len},{peak:.4}]"))
            .collect();
        format!(
            "{{\"test\":\"modal_friction_schelleng\",\"diagnostic\":\"sliding_runs\",\"case\":\"{label}\",\"short_run_samples\":{SHORT_RUN_SAMPLES},\"short_runs_per_period\":{},\"runs_len_peak_slip_m_s\":[{}]}}",
            short as f64 / TAIL_PERIODS as f64,
            runs.join(",")
        )
    }
}

/// Bow the string for `RUN_PERIODS` from `from` at `force` and classify the
/// last `TAIL_PERIODS`.
#[allow(
    clippy::cast_precision_loss,
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss
)]
fn bow(force: f64, from: Start) -> Verdict {
    let mut system = bowed_string(from);
    let external = vec![0.0; MODES];
    let drive = FrictionDrive {
        speed_m_s: BOW_SPEED_M_S,
        normal_force_n: force,
    };
    let samples = (f64::from(RATE_HZ) * RUN_PERIODS / F0_HZ).round() as usize;
    let tail = (f64::from(RATE_HZ) * TAIL_PERIODS as f64 / F0_HZ).round() as usize;
    let mut sliding = Vec::with_capacity(samples);
    let mut force_at_bridge = Vec::with_capacity(samples);
    let mut slip_speed = Vec::with_capacity(samples);
    for _ in 0..samples {
        let frame = system.step(&external, drive).unwrap();
        let (regime, slip_m) = (frame.regime, frame.slip_distance_m);
        sliding.push(matches!(regime, ModalFrictionRegime::Sliding));
        slip_speed.push(slip_m.abs() * f64::from(RATE_HZ));
        force_at_bridge.push(bridge_force(system.components()[0].states()));
    }
    let start = samples - tail;
    let onsets: Vec<usize> = (start.max(1)..samples)
        .filter(|&n| sliding[n] && !sliding[n - 1])
        .collect();
    let step_s = 1.0 / f64::from(RATE_HZ);
    let worst_interval_error = onsets
        .windows(2)
        .map(|w| ((w[1] - w[0]) as f64 * step_s * F0_HZ - 1.0).abs())
        .fold(if onsets.len() < 2 { f64::INFINITY } else { 0.0 }, f64::max);
    let mean_period_s = if onsets.len() < 2 {
        f64::NAN
    } else {
        (onsets[onsets.len() - 1] - onsets[0]) as f64 * step_s / (onsets.len() - 1) as f64
    };
    // Harmonics over a whole number of measured periods (first to last onset),
    // so each is evaluated at an exact multiple of the motion's own frequency.
    let harmonics = if onsets.len() < 2 {
        vec![0.0; HARMONICS]
    } else {
        let (a, b) = (onsets[0], onsets[onsets.len() - 1]);
        let period_samples = (b - a) as f64 / (onsets.len() - 1) as f64;
        (1..=HARMONICS)
            .map(|h| {
                let (mut re, mut im) = (0.0, 0.0);
                for (j, f) in force_at_bridge[a..b].iter().enumerate() {
                    let phase = 2.0 * core::f64::consts::PI * h as f64 * j as f64 / period_samples;
                    let (s, c) = det::sin_cos(phase);
                    re += f * c;
                    im -= f * s;
                }
                det::hypot(re, im)
            })
            .collect()
    };
    let sliding_fraction = sliding[start..].iter().filter(|&&s| s).count() as f64 / tail as f64;
    // Walk back from the run's last onset while each interval is one period.
    let all: Vec<usize> = (1..samples)
        .filter(|&n| sliding[n] && !sliding[n - 1])
        .collect();
    let mut first_settled = None;
    for w in all.windows(2).rev() {
        if ((w[1] - w[0]) as f64 * step_s * F0_HZ - 1.0).abs() > INTERVAL_TOLERANCE {
            break;
        }
        first_settled = Some(w[0]);
    }
    let settled_by_period = first_settled.map(|n| (n as f64 * step_s * F0_HZ) as usize);
    let runs = onsets
        .iter()
        .map(|&n| {
            let end = (n..samples).find(|&m| !sliding[m]).unwrap_or(samples);
            (
                end - n,
                slip_speed[n..end].iter().copied().fold(0.0, f64::max),
            )
        })
        .collect();
    Verdict {
        onsets: onsets.len(),
        worst_interval_error,
        mean_period_s,
        harmonics,
        sliding_fraction,
        settled_by_period,
        runs,
    }
}

#[test]
fn schelleng_window_holds_the_fixture_ratio() {
    let ratio = schelleng_fmax() / schelleng_fmin();
    assert!(
        (ratio - 4.0 * BETA * R_OVER_Z0).abs() <= 1e-12 * ratio,
        "Eq. (10): {ratio}"
    );
}

#[test]
fn single_port_route_sustains_helmholtz_inside_schellengs_window() {
    let force = det::sqrt(schelleng_fmin() * schelleng_fmax());
    let verdict = bow(force, Start::Helmholtz);
    println!("{}", verdict.json("geometric_mean", force));
    println!("{}", verdict.runs_json("geometric_mean"));
    let printed = 2.0 * schelleng_fmin();
    let two = bow(printed, Start::Helmholtz);
    println!("{}", two.json("two_fmin_not_gated", printed));
    println!("{}", two.runs_json("two_fmin_not_gated"));
    let rest = bow(force, Start::Rest);
    println!("{}", rest.json("geometric_mean_from_rest_not_gated", force));
    println!("{}", rest.runs_json("geometric_mean_from_rest_not_gated"));
    assert!(
        verdict.helmholtz(),
        "no Helmholtz motion at sqrt(F_min F_max) = {force} N: {}",
        verdict.json("geometric_mean", force)
    );
}

#[test]
fn planted_negative_below_minimum_force_is_not_helmholtz() {
    let force = schelleng_fmin() / 4.0;
    let verdict = bow(force, Start::Helmholtz);
    println!("{}", verdict.json("quarter_fmin", force));
    println!("{}", verdict.runs_json("quarter_fmin"));
    assert!(
        !verdict.helmholtz(),
        "Helmholtz motion survived at F_min/4 = {force} N: {}",
        verdict.json("quarter_fmin", force)
    );
}
