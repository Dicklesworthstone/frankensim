//! Bowed-string e2e fixture: a modal exact-ZOH string driven by a Stribeck
//! friction island at a bow station (bead frankensim-music-v8-root-3ez8g.7.5).
//!
//! The Helmholtz corner must EMERGE from the friction island acting on the
//! modal string; nothing here injects a corner, a sawtooth, or a forced
//! oscillator (the program NEVER list forbids that neighbor explicitly).
//!
//! Physics summary:
//! - String: mass-normalized retained modes `omega_k = k pi c / L`,
//!   `c = sqrt(T/mu)`, unit-mass shapes `phi_k(x) = sin(k pi x / L) /
//!   sqrt(mu L / 2)`, stepped by the exact-ZOH [`ModalAcousticTimeModel`].
//! - Bow: constant-velocity driver with EVENT-BASED COULOMB STICTION.
//!   While stuck the contact is pinned to `v_str = v_bow` by exactly the
//!   force the modal string demands (acceleration cancellation plus a
//!   velocity-error correction), breaking away when that force exceeds
//!   `mu_static * F_n`; while slipping the traction follows the kinetic
//!   Stribeck curve `mu_k + (mu_s - mu_k) exp(-(v/v0)^2)`. Capture fires
//!   when the relative velocity changes sign between sub-steps or is
//!   already inside `stiction_m_s`: a regularized ramp alone can never
//!   capture this contact (point mass ~2e-5 kg crosses any velocity window
//!   within one sub-step).
//! - Bridge: rigid v1 logs only the string; the second configuration drives
//!   a one-port [`CompactBody`] with the transmitted bridge force (rigid,
//!   massless bridge approximation; the reaction is one-way and disclosed).
//!
//! Determinism class: ONE-HOST. The kinetic curve uses platform `exp`
//! through the existing [`StribeckFriction`] coupling layer; bit-replay is
//! guaranteed only on the same host/library, which every gate states.

use crate::modal_acoustic_time::{
    ModalAcousticMode, ModalAcousticTimeBudget, ModalAcousticTimeError, ModalAcousticTimeModel,
};
use crate::stribeck_friction::StribeckFriction;
use crate::thin_plate::CompactBody;

/// Typed refusal for a bow gesture outside its physical domain.
#[derive(Clone, Debug, PartialEq)]
pub enum BowGestureError {
    /// A field was non-finite.
    NonFinite {
        /// Offending field name.
        what: &'static str,
    },
    /// Normal force must be strictly positive.
    NonPositiveNormalForce,
    /// Station fraction must lie strictly inside the string.
    StationOutOfRange,
}

/// Gesture inputs of a bowed run, all caller-supplied and admitted.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BowGesture {
    /// Bow speed [m/s], positive in this fixture.
    pub v_bow_m_s: f64,
    /// Normal bow force [N], strictly positive.
    pub normal_force_n: f64,
    /// Bow station as a fraction of string length measured from the bridge,
    /// strictly inside `(0, 1)`; violin-like values are `0.05..=0.15`.
    pub station_fraction: f64,
}

impl BowGesture {
    /// Admit a gesture, refusing non-finite fields, non-positive normal
    /// force, and out-of-range stations with typed errors.
    ///
    /// # Errors
    /// Every physical-domain violation documented by [`BowGestureError`].
    pub fn admit(
        v_bow_m_s: f64,
        normal_force_n: f64,
        station_fraction: f64,
    ) -> Result<Self, BowGestureError> {
        if !v_bow_m_s.is_finite() {
            return Err(BowGestureError::NonFinite { what: "v_bow_m_s" });
        }
        if !normal_force_n.is_finite() {
            return Err(BowGestureError::NonFinite {
                what: "normal_force_n",
            });
        }
        if !station_fraction.is_finite() {
            return Err(BowGestureError::NonFinite {
                what: "station_fraction",
            });
        }
        if !(normal_force_n > 0.0) {
            return Err(BowGestureError::NonPositiveNormalForce);
        }
        if !(station_fraction > 0.0 && station_fraction < 1.0) {
            return Err(BowGestureError::StationOutOfRange);
        }
        Ok(Self {
            v_bow_m_s,
            normal_force_n,
            station_fraction,
        })
    }
}

/// Shared string card (same family as the bakeoff string, more modes).
#[derive(Clone, Debug)]
pub struct BowedStringCard {
    /// Speaking length [m].
    pub length_m: f64,
    /// Tension [N].
    pub tension_n: f64,
    /// Linear mass density [kg/m].
    pub linear_density_kg_m: f64,
    /// Retained transverse modes (all below the Nyquist guard).
    pub mode_count: usize,
    /// Per-mode viscous damping ratios, length `mode_count`.
    pub zetas: Vec<f64>,
    /// Sample rate [Hz].
    pub sample_rate_hz: u32,
}

impl BowedStringCard {
    /// Wave speed `sqrt(T/mu)` [m/s].
    #[must_use]
    pub fn wave_speed_m_s(&self) -> f64 {
        (self.tension_n / self.linear_density_kg_m).sqrt()
    }

    /// Transverse fundamental [Hz].
    #[must_use]
    pub fn fundamental_hz(&self) -> f64 {
        self.wave_speed_m_s() / (2.0 * self.length_m)
    }

    /// Validate the card before admission.
    ///
    /// # Errors
    /// Refuses a zeta count disagreeing with `mode_count` or non-physical
    /// scalars.
    pub fn validate(&self) -> Result<(), BowedRunError> {
        if self.zetas.len() != self.mode_count {
            return Err(BowedRunError::InvalidCard {
                what: format!(
                    "zeta count {} != mode_count {}",
                    self.zetas.len(),
                    self.mode_count
                ),
            });
        }
        if !(self.length_m > 0.0
            && self.tension_n > 0.0
            && self.linear_density_kg_m > 0.0
            && self.sample_rate_hz > 0
            && self.mode_count > 0)
        {
            return Err(BowedRunError::InvalidCard {
                what: "card scalars must be positive".to_string(),
            });
        }
        Ok(())
    }
}

/// Friction island variants. Under the event-stiction coupling the
/// [`StribeckFriction`] fields mean: `mu_static` caps the pinning force,
/// `mu_dynamic` is the sliding floor, and `stiction_m_s` is the KINETIC
/// decay width of the drop — it must span cm/s so a sliding operating point
/// samples the negative-slope region and pumps the oscillation.
///
/// `ViscousOnly` is the FALSIFIER law: purely viscous opposition has no
/// stiction window and no Stribeck drop, so stick-slip gates must fail.
#[derive(Clone, Copy, Debug)]
pub enum FrictionIsland {
    /// The regularized Stribeck rung used across the music stack.
    Stribeck(StribeckFriction),
    /// Pure viscous opposition `c * v_rel`; falsifier only.
    ViscousOnly {
        /// Viscous coefficient [N s/m].
        viscous_n_s_per_m: f64,
    },
}

impl FrictionIsland {
    fn capture_tol_m_s(self) -> f64 {
        match self {
            Self::Stribeck(law) => law.stiction_m_s.max(0.02),
            Self::ViscousOnly { .. } => 0.0,
        }
    }

    fn kinetic_traction(self, v_rel_m_s: f64, normal_force_n: f64) -> f64 {
        match self {
            Self::Stribeck(law) => law.traction(v_rel_m_s, normal_force_n),
            Self::ViscousOnly {
                viscous_n_s_per_m: cst,
            } => cst * v_rel_m_s,
        }
    }
}

/// Bridge termination configuration of a run.
#[derive(Clone)]
pub enum Termination {
    /// Rigid bridge; nothing but the string is logged.
    Rigid,
    /// The bridge force drives this one-port body (one-way rigid-bridge
    /// approximation, disclosed in the run receipt).
    PlateOnePort(Box<CompactBody>),
}

/// Full run configuration after admission.
#[derive(Clone)]
pub struct BowedRunConfig {
    /// String card.
    pub card: BowedStringCard,
    /// Friction island.
    pub island: FrictionIsland,
    /// Admitted gesture.
    pub gesture: BowGesture,
    /// Number of samples to run.
    pub steps: usize,
    /// Exact-ZOH sub-steps per sample refreshing the friction coupling.
    pub subsamples: usize,
    /// Bridge termination.
    pub termination: Termination,
    /// Listener distance for radiation logging [m] (plate configuration).
    pub listener_m: f64,
}

/// One completed bowed run: per-sample histories at the bow point and the
/// bridge, ready for gate analysis. Everything is plain `f64` in fixed
/// order, so equal configurations replay bitwise on one host.
#[derive(Clone, Debug)]
pub struct BowedRunLog {
    /// String transverse velocity at the bow station [m/s] per sample.
    pub bow_point_velocity_m_s: Vec<f64>,
    /// Relative velocity `v_bow - v_string` at the station [m/s] per sample.
    pub relative_velocity_m_s: Vec<f64>,
    /// Transmitted bridge force [N] per sample.
    pub bridge_force_n: Vec<f64>,
    /// Body velocity [m/s] per sample; empty for [`Termination::Rigid`].
    pub body_velocity_m_s: Vec<f64>,
    /// Radiated pressure at the listener [Pa]; empty for rigid.
    pub radiated_pressure_pa: Vec<f64>,
    /// Final total modal energy [J].
    pub final_total_energy_j: f64,
    /// Peak total modal energy observed during the run [J].
    pub peak_total_energy_j: f64,
}

/// Typed refusal from run admission or an underlying model error.
#[derive(Debug)]
pub enum BowedRunError {
    /// The gesture failed admission.
    Gesture(BowGestureError),
    /// The card is inconsistent (zeta count, non-physical values).
    InvalidCard {
        /// Human-readable cause.
        what: String,
    },
    /// The modal runtime refused admission or a step.
    Model(ModalAcousticTimeError),
    /// A state limit was exceeded mid-run; earlier samples remain valid.
    LimitExceeded(ModalAcousticTimeError),
}

impl std::fmt::Display for BowedRunError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Gesture(e) => write!(f, "gesture refused: {e:?}"),
            Self::InvalidCard { what } => write!(f, "card refused: {what}"),
            Self::Model(e) => write!(f, "model refused: {e:?}"),
            Self::LimitExceeded(e) => write!(f, "state limit exceeded: {e:?}"),
        }
    }
}

impl std::error::Error for BowedRunError {}

fn unit_shape(k: usize, x_fraction: f64, length_m: f64, mu: f64) -> f64 {
    let norm = (mu * length_m * 0.5).sqrt();
    ((k + 1) as f64 * core::f64::consts::PI * x_fraction).sin() / norm
}

fn unit_shape_slope_at_bridge(k: usize, length_m: f64, mu: f64) -> f64 {
    let norm = (mu * length_m * 0.5).sqrt();
    let kk = (k + 1) as f64;
    let slope = kk * core::f64::consts::PI / length_m * (kk * core::f64::consts::PI).cos();
    slope / norm
}

/// Run the admitted configuration, returning per-sample histories.
///
/// Coupling discipline: the island interaction refreshes on
/// `config.subsamples` exact-ZOH sub-steps per audio sample
/// (`step_duration`). Same host, same order, bitwise replayable; no
/// wall-clock enters anywhere.
///
/// A state-limit refusal aborts the run and is returned as
/// [`BowedRunError::LimitExceeded`]; callers treat blown limits as data
/// (the raucous regime can genuinely blow a displacement ceiling).
///
/// # Errors
/// See [`BowedRunError`].
pub fn run_bowed(config: &BowedRunConfig) -> Result<BowedRunLog, BowedRunError> {
    config.card.validate()?;
    BowGesture::admit(
        config.gesture.v_bow_m_s,
        config.gesture.normal_force_n,
        config.gesture.station_fraction,
    )
    .map_err(BowedRunError::Gesture)?;
    let card = &config.card;
    let mu = card.linear_density_kg_m;
    let c = card.wave_speed_m_s();

    let modes: Vec<ModalAcousticMode> = (0..card.mode_count)
        .map(|k| ModalAcousticMode {
            angular_frequency_rad_s: (k + 1) as f64 * core::f64::consts::PI * c / card.length_m,
            damping_ratio: card.zetas[k],
            pressure_per_modal_velocity: fs_math::c64::C64::new(1.0, 0.0),
        })
        .collect();
    let mut model = ModalAcousticTimeModel::try_new(
        card.sample_rate_hz,
        modes,
        ModalAcousticTimeBudget::audible_reference(),
    )
    .map_err(BowedRunError::Model)?;

    let x_bow = config.gesture.station_fraction;
    let shapes_at_bow: Vec<f64> = (0..card.mode_count)
        .map(|k| unit_shape(k, x_bow, card.length_m, mu))
        .collect();
    let slopes_at_bridge: Vec<f64> = (0..card.mode_count)
        .map(|k| unit_shape_slope_at_bridge(k, card.length_m, mu))
        .collect();
    let w_point: f64 = shapes_at_bow.iter().map(|phi| phi * phi).sum();

    let mut log = BowedRunLog {
        bow_point_velocity_m_s: Vec::with_capacity(config.steps),
        relative_velocity_m_s: Vec::with_capacity(config.steps),
        bridge_force_n: Vec::with_capacity(config.steps),
        body_velocity_m_s: Vec::new(),
        radiated_pressure_pa: Vec::new(),
        final_total_energy_j: 0.0,
        peak_total_energy_j: f64::MIN,
    };
    let mut body: Option<CompactBody> = match &config.termination {
        Termination::PlateOnePort(boxed) => Some(boxed.as_ref().clone()),
        Termination::Rigid => None,
    };
    if body.is_some() {
        log.body_velocity_m_s.reserve(config.steps);
        log.radiated_pressure_pa.reserve(config.steps);
    }
    // Start from rest; the bow spins the string up from silence so every
    // periodicity in the logs is EMERGENT, never an initial-condition echo.
    let mut forces = vec![0.0_f64; card.mode_count];
    let dt = model.sample_period_s();
    let subsamples = config.subsamples.max(1);
    let sub_dt = dt / subsamples as f64;
    let capture_tol = config.island.capture_tol_m_s();
    // Capture basin: both sides of a sign flip must be slower than this,
    // or the crossing is a chatter spike rather than a capturable corner.
    const CAPTURE_BASIN_M_S: f64 = 0.15;
    let mut stuck = false;
    let mut prev_v_rel = f64::NAN;
    let mut peak_energy_j = f64::MIN;

    for _ in 0..config.steps {
        let mut v_str = 0.0_f64;
        let mut v_rel = 0.0_f64;
        for _ in 0..subsamples {
            v_str = model
                .states()
                .iter()
                .zip(&shapes_at_bow)
                .map(|(s, phi)| phi * s.velocity_m_sqrt_kg_per_s)
                .sum();
            v_rel = config.gesture.v_bow_m_s - v_str;
            let flip_speed = prev_v_rel.abs().max(v_rel.abs());
            let flipped = prev_v_rel.is_finite()
                && !stuck
                && prev_v_rel.signum() != v_rel.signum()
                && flip_speed <= CAPTURE_BASIN_M_S;
            prev_v_rel = v_rel;
            let traction = match config.island {
                FrictionIsland::Stribeck(law) => {
                    let hold_cap = law.mu_static * config.gesture.normal_force_n;
                    // Total pinning force: cancel contact acceleration AND
                    // drive the velocity error to zero this sub-step.
                    // Minimum-norm projection through the mode shapes.
                    let pin = || -> f64 {
                        let accel_hold: f64 = model
                            .modes()
                            .iter()
                            .zip(model.states())
                            .zip(&shapes_at_bow)
                            .map(|((m, s), phi)| {
                                phi
                                    * (2.0 * m.damping_ratio * m.angular_frequency_rad_s
                                        * s.velocity_m_sqrt_kg_per_s
                                        + m.angular_frequency_rad_s
                                            * m.angular_frequency_rad_s
                                            * s.displacement_m_sqrt_kg)
                            })
                            .sum::<f64>()
                            / w_point;
                        accel_hold + v_rel / (w_point * sub_dt)
                    };
                    if stuck {
                        let p = pin();
                        if p.abs() <= hold_cap {
                            p
                        } else {
                            stuck = false;
                            law.traction(v_rel, config.gesture.normal_force_n)
                        }
                    } else if flipped || v_rel.abs() <= capture_tol {
                        let p = pin();
                        if p.abs() <= hold_cap {
                            stuck = true;
                            p
                        } else {
                            law.traction(v_rel, config.gesture.normal_force_n)
                        }
                    } else {
                        law.traction(v_rel, config.gesture.normal_force_n)
                    }
                }
                FrictionIsland::ViscousOnly {
                    viscous_n_s_per_m: cst,
                } => cst * v_rel,
            };
            for (q, phi) in forces.iter_mut().zip(&shapes_at_bow) {
                *q = traction * phi;
            }
            let stepped = model
                .step_duration(&forces, sub_dt)
                .map_err(BowedRunError::LimitExceeded)?;
            peak_energy_j = peak_energy_j.max(stepped.total_modal_energy_j);
        }

        let bridge_force: f64 = model
            .states()
            .iter()
            .zip(&slopes_at_bridge)
            .map(|(s, slope)| card.tension_n * slope * s.displacement_m_sqrt_kg)
            .sum();

        log.bow_point_velocity_m_s.push(v_str);
        log.relative_velocity_m_s.push(v_rel);
        log.bridge_force_n.push(bridge_force);

        if let Some(body) = body.as_mut() {
            let acc = body.drive(bridge_force, dt);
            log.body_velocity_m_s.push(body.volume_velocity());
            log.radiated_pressure_pa
                .push(body.radiate(acc, AIR_DENSITY_KG_M3, config.listener_m));
        }
    }
    log.final_total_energy_j = frame_total_energy(&model);
    log.peak_total_energy_j = peak_energy_j;
    Ok(log)
}

const AIR_DENSITY_KG_M3: f64 = 1.204;

fn frame_total_energy(model: &ModalAcousticTimeModel) -> f64 {
    model
        .modes()
        .iter()
        .zip(model.states())
        .map(|(m, s)| {
            0.5 * s.velocity_m_sqrt_kg_per_s.powi(2)
                + 0.5
                    * m.angular_frequency_rad_s
                    * m.angular_frequency_rad_s
                    * s.displacement_m_sqrt_kg.powi(2)
        })
        .sum()
}

/// Bitwise content hash of a run log (f64 bits, fixed order): the replay
/// identity used by the determinism gate.
#[must_use]
pub fn run_log_hash(log: &BowedRunLog) -> fs_blake3::ContentHash {
    let mut bytes = Vec::with_capacity(
        (log.bow_point_velocity_m_s.len()
            + log.relative_velocity_m_s.len()
            + log.bridge_force_n.len()
            + log.body_velocity_m_s.len()
            + log.radiated_pressure_pa.len())
            * 8
            + 16,
    );
    for series in [
        &log.bow_point_velocity_m_s,
        &log.relative_velocity_m_s,
        &log.bridge_force_n,
        &log.body_velocity_m_s,
        &log.radiated_pressure_pa,
    ] {
        for value in series.iter() {
            bytes.extend_from_slice(&value.to_bits().to_le_bytes());
        }
    }
    bytes.extend_from_slice(&log.final_total_energy_j.to_bits().to_le_bytes());
    bytes.extend_from_slice(&log.peak_total_energy_j.to_bits().to_le_bytes());
    fs_blake3::hash_bytes(&bytes)
}

/// Goertzel magnitude of `samples` at `frequency_hz` for `sample_rate_hz`.
#[must_use]
pub fn goertzel_magnitude(samples: &[f64], frequency_hz: f64, sample_rate_hz: f64) -> f64 {
    if samples.len() < 2 || !frequency_hz.is_finite() || !(sample_rate_hz > 0.0) {
        return 0.0;
    }
    let n = samples.len() as f64;
    let omega = core::f64::consts::TAU * frequency_hz / sample_rate_hz;
    let coeff = 2.0 * omega.cos();
    let (mut s1, mut s2) = (0.0_f64, 0.0_f64);
    for x in samples {
        let s0 = x + coeff * s1 - s2;
        s2 = s1;
        s1 = s0;
    }
    (s1 * s1 + s2 * s2 - coeff * s1 * s2).sqrt() / (0.5 * n)
}

/// Quadratic-refined peak frequency [Hz] of a Goertzel scan over
/// `f_min..=f_max` at `probe_count` bins, refined by three-point parabolic
/// interpolation.
#[must_use]
pub fn refined_peak_hz(
    samples: &[f64],
    f_min: f64,
    f_max: f64,
    probe_count: usize,
    sample_rate_hz: f64,
) -> f64 {
    if probe_count < 3 || !(f_max > f_min) || samples.len() < 2 {
        return f_min;
    }
    let step = (f_max - f_min) / (probe_count as f64 - 1.0);
    let mut best = 0usize;
    let mut best_mag = -1.0_f64;
    let mut mags = Vec::with_capacity(probe_count);
    for i in 0..probe_count {
        let f = f_min + step * i as f64;
        let m = goertzel_magnitude(samples, f, sample_rate_hz);
        if m > best_mag {
            best_mag = m;
            best = i;
        }
        mags.push(m);
    }
    if best == 0 || best + 1 == probe_count {
        return f_min + step * best as f64;
    }
    let alpha = mags[best - 1];
    let beta = mags[best];
    let gamma = mags[best + 1];
    let denom = alpha - 2.0 * beta + gamma;
    let delta = if denom.abs() <= f64::MIN_POSITIVE {
        0.0
    } else {
        ((alpha - gamma) / (2.0 * denom)).clamp(-1.0, 1.0)
    };
    f_min + step * (best as f64 + delta)
}

/// Count contiguous slipping intervals per nominal period, averaged over
/// the analysis window. Helmholtz motion gives ONE flyback per period;
/// doubled slipping gives two; raucous motion gives many or none.
///
/// Slip means `|v_rel| > slip_threshold_m_s`.
#[must_use]
pub fn mean_slip_intervals_per_period(
    relative_velocity_m_s: &[f64],
    fundamental_hz: f64,
    sample_rate_hz: f64,
    slip_threshold_m_s: f64,
) -> f64 {
    if relative_velocity_m_s.len() < 16 || !(fundamental_hz > 0.0) || !(sample_rate_hz > 0.0) {
        return f64::NAN;
    }
    let period_samples = (sample_rate_hz / fundamental_hz).round();
    if !(period_samples >= 4.0) {
        return f64::NAN;
    }
    let period = period_samples as usize;
    let windows = relative_velocity_m_s.len() / period;
    if windows == 0 {
        return f64::NAN;
    }
    let mut total_intervals = 0.0_f64;
    for w in 0..windows {
        let start = w * period;
        let mut intervals = 0u32;
        let mut in_slip = false;
        for i in 0..period {
            let slipping = relative_velocity_m_s[start + i].abs() > slip_threshold_m_s;
            if slipping && !in_slip {
                intervals += 1;
            }
            in_slip = slipping;
        }
        total_intervals += f64::from(intervals);
    }
    total_intervals / windows as f64
}

/// Fraction of samples spent slipping (`|v_rel| > threshold`) in the window.
#[must_use]
pub fn slip_fraction(relative_velocity_m_s: &[f64], slip_threshold_m_s: f64) -> f64 {
    if relative_velocity_m_s.is_empty() {
        return f64::NAN;
    }
    let slipped = relative_velocity_m_s
        .iter()
        .filter(|v| v.abs() > slip_threshold_m_s)
        .count();
    slipped as f64 / relative_velocity_m_s.len() as f64
}

/// Pitch of the strongest spectral component in `[f_min, f_max]`, expressed
/// in cents relative to `reference_hz` (positive = sharper).
#[must_use]
pub fn cents_deviation(
    samples: &[f64],
    f_min: f64,
    f_max: f64,
    reference_hz: f64,
    sample_rate_hz: f64,
) -> f64 {
    let peak = refined_peak_hz(samples, f_min, f_max, 241, sample_rate_hz);
    1200.0 * (peak / reference_hz).log2()
}

/// Measured emergent-gate observables for one completed run window.
#[derive(Clone, Copy, Debug)]
pub struct GateMetrics {
    /// Refined dominant-peak frequency [Hz] near the fundamental.
    pub peak_hz: f64,
    /// Nominal string fundamental [Hz].
    pub fundamental_hz: f64,
    /// Fundamental-bin magnitude over the semitone-neighbor magnitude.
    pub peak_to_semitone_ratio: f64,
    /// Fraction of window samples slipping.
    pub slip_frac: f64,
    /// Mean contiguous slip intervals per nominal period.
    pub intervals_per_period: f64,
}

fn tail_window(series: &[f64], fraction: f64) -> &[f64] {
    let start = series.len() - ((series.len() as f64 * fraction) as usize);
    &series[start..]
}

/// Compute [`GateMetrics`] over the steady tail (last 40%) of a run log.
#[must_use]
pub fn gate_metrics(log: &BowedRunLog, card: &BowedStringCard) -> GateMetrics {
    let sr = f64::from(card.sample_rate_hz);
    let f1 = card.fundamental_hz();
    let win = tail_window(&log.bow_point_velocity_m_s, 0.4);
    let semitone = f1 * 2.0_f64.powf(1.0 / 12.0);
    let at_f1 = goertzel_magnitude(win, f1, sr);
    let at_neighbor = goertzel_magnitude(win, semitone, sr);
    let lo = f1 * 0.9;
    let hi = f1 * 1.12;
    let probes = 97_usize;
    let mut peak = lo;
    let mut best_mag = -1.0_f64;
    for i in 0..probes {
        let f = lo + (hi - lo) * i as f64 / (probes as f64 - 1.0);
        let m = goertzel_magnitude(win, f, sr);
        if m > best_mag {
            best_mag = m;
            peak = f;
        }
    }
    let rel = tail_window(&log.relative_velocity_m_s, 0.4);
    let mean_abs_rel = rel.iter().map(|v| v.abs()).sum::<f64>() / rel.len() as f64;
    let slip_threshold = 0.15 * mean_abs_rel + 1.0e-4;
    GateMetrics {
        peak_hz: peak,
        fundamental_hz: f1,
        peak_to_semitone_ratio: if at_neighbor > 0.0 {
            at_f1 / at_neighbor
        } else {
            f64::INFINITY
        },
        slip_frac: slip_fraction(rel, slip_threshold),
        intervals_per_period: mean_slip_intervals_per_period(rel, f1, sr, slip_threshold),
    }
}

/// Regime classification of measured metrics. The PLAYABLE band sits
/// between a minimum and a maximum bow force: above it motion goes
/// raucous/aperiodic, below it decays toward surface sound.
///
/// Declared rig tolerances (not literature constants): pitch within 6% of
/// the string fundamental, fundamental dominance over its semitone
/// neighbor above 2x, bounded slip fraction, single flyback interval per
/// period within measurement slack.
#[must_use]
pub fn classify(m: &GateMetrics) -> &'static str {
    let near_fundamental = (m.peak_hz / m.fundamental_hz - 1.0).abs() < 0.06;
    let selective = m.peak_to_semitone_ratio > 2.0;
    let bounded_slip = m.slip_frac > 0.02 && m.slip_frac < 0.85;
    let single_flyback = m.intervals_per_period > 0.5 && m.intervals_per_period < 1.6;
    if near_fundamental && selective && bounded_slip && single_flyback {
        "playable"
    } else if m.slip_frac >= 0.85 || m.intervals_per_period >= 1.6 {
        "raucous"
    } else {
        "surface-or-dead"
    }
}
