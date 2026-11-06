//! Spontaneous ground-state radiative recombination reciprocal to a photon source.
//!
//! The Milne relation is sigma_rr=(g_A/g_I)*(E_gamma/(m_e*c*v))^2*sigma_pi,
//! with E_gamma=chi+epsilon and epsilon=m_e*v^2/2. There is no additional
//! electron-spin factor. Averaging v*sigma_rr over Maxwellian electrons gives
//! alpha=C/(kT)^(3/2)*integral E_gamma^2*sigma_pi*exp(-epsilon/kT) d epsilon,
//! C=(g_A/g_I)*sqrt(2/pi)/(m_e^(3/2)*c^2). Inserting epsilon gives the captured
//! kinetic-energy moment L; emitted photon energy is chi*alpha+L.
//! Reference: Nahar & Pradhan, ApJ 447 (1995), Eqs.(1)-(2):
//! <https://www.astronomy.ohio-state.edu/nahar.1/papers/apj-unified-sis.pdf>.
//!
//! Only this original ground-state channel is reversed. Nonresonant input
//! describes direct radiative capture; any supplied resonances remain in the
//! reciprocal result. No additional dielectronic channel, excited-level cascade,
//! stimulated recombination, three-body capture, recoil or radiation transport
//! is inferred. Ions are stationary massive targets. This is neither a total
//! Case A nor Case B recombination coefficient. An optically
//! thin consumer may let photons escape; a trapping consumer must transport them.
//! The source's declared zero below its onset is retained, including the
//! rounded Verner onset; it is not replaced with the thermochemical gap.
//!
//! Positive Gauss rules split at every source knot and thermal-scale panels.
//! Quadrature differences are estimates, not certificates. A separately sourced
//! GLOBAL cross-section envelope bounds all omitted Maxwellian moments, including
//! energies above source coverage. Insufficient coverage or work refuses rather
//! than silently treating unintegrated photons/electrons as absent.

use super::AtomicIonizationErrorV1;
use super::photoionization::{
    AtomicPhotoionizationV1, PhotoionizationCrossSectionPointV1, PhotoionizationErrorV1,
};
use crate::electron_collision::{
    AVOGADRO_PER_MOL, BOLTZMANN_J_PER_K, ELECTRON_MASS_KG, MaxwellianRateControlsV1,
    integrate::{W8, W16, X8, X16},
};
use fs_matdb::Provenance;
use fs_math::det::{exp, ln};
use fs_qty::{Area, Energy, Temperature};

const LIGHT_SPEED_M_S: f64 = 299_792_458.0;

/// V1 Maxwellian admission ceiling kT/(m_e*c²), not a physical-error bound.
pub const MAX_RADIATIVE_RECOMBINATION_THERMAL_RATIO_V1: f64 = 1e-3;

/// Refusal before any state change or complete numerical result.
#[derive(Clone, Debug, PartialEq)]
pub enum RadiativeRecombinationErrorV1 {
    /// Invalid envelope, provenance, temperature or controls.
    Invalid(&'static str),
    /// Original photon source refused an evaluation.
    Photoionization(PhotoionizationErrorV1),
    /// Original electron caloric source refused the temperature.
    Atomic(AtomicIonizationErrorV1),
    /// Total Gauss-panel budget exhausted, including retries.
    QuadratureBudget,
    /// The declared envelope cannot resolve both moments inside source support.
    UnresolvedTail {
        /// Maximum omitted-tail bound relative to its retained moment.
        relative_bound: f64,
    },
    /// Positive moments or source scales cannot be represented in binary64.
    Arithmetic(&'static str),
    /// Cancellation before publication.
    Cancelled,
}
impl core::fmt::Display for RadiativeRecombinationErrorV1 {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "radiative recombination: {self:?}")
    }
}
impl std::error::Error for RadiativeRecombinationErrorV1 {}

type Error = RadiativeRecombinationErrorV1;

/// Immutable reciprocal source; every successful evaluation allocates nothing.
#[derive(Clone, Debug, PartialEq)]
pub struct AtomicRadiativeRecombinationV1 {
    source: AtomicPhotoionizationV1,
    global_cross_section_upper_bound: Area,
    envelope_provenance: Provenance,
    scale: f64,
}
impl AtomicRadiativeRecombinationV1 {
    /// Assert 0<=sigma_pi(E)<=bound for ALL photon energies, including beyond
    /// the original source. The known interval is checked; the continuation
    /// remains an explicit caller/source assertion, retained with its provenance.
    /// A zero bound is admitted only for an identically zero original source.
    pub fn new(
        source: AtomicPhotoionizationV1,
        global_cross_section_upper_bound: Area,
        envelope_provenance: Provenance,
    ) -> Result<Self, Error> {
        let bound = global_cross_section_upper_bound.value();
        if !bound.is_finite()
            || bound < 0.0
            || envelope_provenance.source.trim().is_empty()
            || envelope_provenance.license.trim().is_empty()
        {
            return Err(Error::Invalid("global photon envelope and its provenance"));
        }
        let scale = source
            .maximum_cross_section_in_domain()
            .map_err(Error::Photoionization)?
            .value();
        if bound < scale {
            return Err(Error::Invalid(
                "global envelope below original cross section",
            ));
        }
        Ok(Self {
            source,
            global_cross_section_upper_bound,
            envelope_provenance,
            scale,
        })
    }

    /// Original photon source, including its complete atomic species bank.
    pub const fn source(&self) -> &AtomicPhotoionizationV1 {
        &self.source
    }
    /// Caller/source assertion for all energies, never used as an extrapolation.
    pub const fn global_cross_section_upper_bound(&self) -> Area {
        self.global_cross_section_upper_bound
    }
    /// Original evidence/usage declaration for that global envelope.
    pub const fn envelope_provenance(&self) -> &Provenance {
        &self.envelope_provenance
    }

    /// Maxwellian moments at the actual electron temperature, per ion-electron
    /// number-density pair. For molar concentrations q=NA*alpha*c_I*c_e is the
    /// recombination progress in mol/(m³ s); L*(NA*c_I)*(NA*c_e) is the electron
    /// kinetic-energy loss in W/m³. It is not 3*kT/2 times the event rate.
    ///
    /// Fixed-field derivatives are ordinary binary mass action; this evaluation
    /// neither consumes populations nor owns photon escape/thermal evolution.
    /// Original electron caloric temperature support is enforced here, together
    /// with the explicit V1 nonrelativistic Maxwellian thermal-energy ceiling.
    /// Both alpha and L must meet the tolerance. Positive underflow refuses;
    /// only an explicitly zero global envelope produces exactly zero moments.
    pub fn evaluate(
        &self,
        electron_temperature: Temperature,
        controls: MaxwellianRateControlsV1,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<RadiativeRecombinationEvaluationV1, Error> {
        if cancelled() {
            return Err(Error::Cancelled);
        }
        controls
            .validate()
            .map_err(|_| Error::Invalid("Maxwellian integration controls"))?;
        if !positive(electron_temperature.value()) {
            return Err(Error::Invalid("electron temperature"));
        }
        let channel = self.source.channel();
        channel
            .admit_source_temperature(channel.electron_index(), electron_temperature)
            .map_err(Error::Atomic)?;
        let thermal = BOLTZMANN_J_PER_K * electron_temperature.value();
        if !positive(thermal)
            || thermal / (ELECTRON_MASS_KG * LIGHT_SPEED_M_S * LIGHT_SPEED_M_S)
                > MAX_RADIATIVE_RECOMBINATION_THERMAL_RATIO_V1
        {
            return Err(Error::Invalid(
                "nonrelativistic Maxwellian thermal-energy limit",
            ));
        }
        if self.scale == 0.0 {
            if self.global_cross_section_upper_bound.value() != 0.0 {
                return Err(Error::UnresolvedTail {
                    relative_bound: f64::INFINITY,
                });
            }
            if cancelled() {
                return Err(Error::Cancelled);
            }
            return Ok(RadiativeRecombinationEvaluationV1 {
                coefficient_m3_s: 0.0,
                captured_kinetic_energy_j_m3_s: 0.0,
                emitted_photon_energy_j_m3_s: 0.0,
                relative_quadrature_error_estimate: 0.0,
                relative_tail_bound: 0.0,
            });
        }
        let gap = channel.ionization_energy_j_mol() / AVOGADRO_PER_MOL;
        let onset = self.source.threshold().value();
        let delta = onset - gap;
        let photon = onset + thermal;
        let kinetic = delta + thermal;
        let onset_ratio = delta / thermal;
        if !positive(thermal)
            || !positive(photon)
            || !positive(kinetic)
            || !onset_ratio.is_finite()
            || onset_ratio < 0.0
        {
            return Err(Error::Arithmetic("Maxwellian source scales"));
        }
        let scales = Scales {
            thermal,
            onset,
            p: [onset / photon, thermal / photon],
            k: [delta / kinetic, thermal / kinetic],
        };
        let (moments, error, tail) = self.integrate(scales, controls, &mut cancelled)?;
        let weights = channel.statistical_weights();
        let log_factor = ln(f64::from(weights.neutral) / f64::from(weights.ion))
            + 0.5 * ln(2.0 / core::f64::consts::PI)
            - 1.5 * ln(ELECTRON_MASS_KG)
            - 2.0 * ln(LIGHT_SPEED_M_S)
            + ln(self.scale)
            + 2.0 * ln(photon)
            - 0.5 * ln(thermal)
            - onset_ratio;
        let alpha = exp(log_factor + ln(moments[0]));
        let captured = exp(log_factor + ln(kinetic) + ln(moments[1]));
        let emitted = gap * alpha + captured;
        if !positive(alpha)
            || !positive(captured)
            || !positive(emitted)
            || !positive(captured / alpha)
        {
            return Err(Error::Arithmetic("positive recombination moments"));
        }
        if cancelled() {
            return Err(Error::Cancelled);
        }
        Ok(RadiativeRecombinationEvaluationV1 {
            coefficient_m3_s: alpha,
            captured_kinetic_energy_j_m3_s: captured,
            emitted_photon_energy_j_m3_s: emitted,
            relative_quadrature_error_estimate: error,
            relative_tail_bound: tail,
        })
    }

    fn integrate(
        &self,
        scales: Scales,
        controls: MaxwellianRateControlsV1,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<([f64; 2], f64, f64), Error> {
        let source_end = (self.source.maximum_energy().value() - scales.onset) / scales.thermal;
        if source_end <= 0.0 || source_end.is_nan() {
            return Err(Error::UnresolvedTail {
                relative_bound: f64::INFINITY,
            });
        }
        let mut end = source_end.min(32.0);
        let mut panels_used = 0_usize;
        loop {
            let mut refinement = 1_usize;
            let (values, error) = loop {
                let mut sums = [Sum::default(); 2];
                let mut errors = [Sum::default(); 2];
                let mut segment = |lower: f64,
                                   upper: f64,
                                   points: Option<[PhotoionizationCrossSectionPointV1; 2]>|
                 -> Result<(), Error> {
                    let lower = lower.max(0.0).min(end);
                    let upper = upper.max(0.0).min(end);
                    if lower >= upper {
                        return Ok(());
                    }
                    let count = (((upper - lower) / 4.0).ceil() as usize)
                        .max(1)
                        .checked_mul(refinement)
                        .ok_or(Error::QuadratureBudget)?;
                    if count > controls.maximum_subdivisions.saturating_sub(panels_used) {
                        return Err(Error::QuadratureBudget);
                    }
                    panels_used += count;
                    let width = (upper - lower) / count as f64;
                    for panel in 0..count {
                        if cancelled() {
                            return Err(Error::Cancelled);
                        }
                        let a = lower + panel as f64 * width;
                        let b = if panel + 1 == count {
                            upper
                        } else {
                            lower + (panel + 1) as f64 * width
                        };
                        let fine = self.rule(scales, points, a, b, &X16, &W16)?;
                        let coarse = self.rule(scales, points, a, b, &X8, &W8)?;
                        for j in 0..2 {
                            sums[j].add(fine[j]);
                            errors[j].add((fine[j] - coarse[j]).abs());
                        }
                    }
                    Ok(())
                };
                if let Some(points) = self.source.points() {
                    for pair in points.windows(2) {
                        let lower = (pair[0].photon_energy.value() - scales.onset) / scales.thermal;
                        if lower >= end {
                            break;
                        }
                        let upper = (pair[1].photon_energy.value() - scales.onset) / scales.thermal;
                        segment(lower, upper, Some([pair[0], pair[1]]))?;
                    }
                } else {
                    segment(0.0, end, None)?;
                }
                let values = sums.map(|sum| sum.value);
                if values.iter().any(|v| !v.is_finite() || *v < 0.0) {
                    return Err(Error::Arithmetic("recombination integral"));
                }
                let error = (0..2)
                    .map(|j| {
                        if values[j] == 0.0 {
                            0.0
                        } else {
                            errors[j].value / values[j]
                        }
                    })
                    .fold(0.0, f64::max);
                if error <= controls.relative_tolerance * 0.25 {
                    break (values, error);
                }
                refinement = refinement.checked_mul(2).ok_or(Error::QuadratureBudget)?;
            };
            let tail = self.relative_tail(end, scales, values);
            if tail <= controls.relative_tolerance * 0.25 {
                return Ok((values, error, tail));
            }
            if end == source_end || end >= 1024.0 {
                return Err(Error::UnresolvedTail {
                    relative_bound: tail,
                });
            }
            end = (2.0 * end).min(source_end);
        }
    }

    fn rule<const N: usize>(
        &self,
        scales: Scales,
        points: Option<[PhotoionizationCrossSectionPointV1; 2]>,
        lower: f64,
        upper: f64,
        nodes: &[f64; N],
        weights: &[f64; N],
    ) -> Result<[f64; 2], Error> {
        let half = 0.5 * (upper - lower);
        let midpoint = lower + half;
        let mut sums = [Sum::default(); 2];
        for (&node, &weight) in nodes.iter().zip(weights) {
            for sign in [-1.0, 1.0] {
                let y = midpoint + sign * half * node;
                let sigma = if let Some([left, right]) = points {
                    // Evaluate within this source segment without another binary
                    // search or losing thermal offsets when onset >> kT.
                    let f = ((scales.onset - left.photon_energy.value()) + scales.thermal * y)
                        / (right.photon_energy.value() - left.photon_energy.value());
                    let f = f.clamp(0.0, 1.0);
                    let a = left.cross_section.value();
                    let b = right.cross_section.value();
                    if a <= b {
                        a + f * (b - a)
                    } else {
                        b + (1.0 - f) * (a - b)
                    }
                } else {
                    let energy = (scales.onset + scales.thermal * y)
                        .min(self.source.maximum_energy().value());
                    self.source
                        .cross_section(Energy::new(energy))
                        .map_err(Error::Photoionization)?
                        .value()
                };
                let photon_ratio = scales.p[0] + scales.p[1] * y;
                let kinetic_ratio = scales.k[0] + scales.k[1] * y;
                let value = weight * (sigma / self.scale) * photon_ratio * photon_ratio * exp(-y);
                sums[0].add(value);
                sums[1].add(value * kinetic_ratio);
            }
        }
        Ok(sums.map(|sum| sum.value * half))
    }

    fn relative_tail(&self, end: f64, scales: Scales, values: [f64; 2]) -> f64 {
        // Exact constant-envelope incomplete-gamma tails, with positive
        // coefficients. Work in log space so exp(-end) cannot hide a tail.
        let f = scales.p[0] + scales.p[1] * end;
        let q = scales.p[1];
        let g = scales.k[0] + scales.k[1] * end;
        let r = scales.k[1];
        let rate = f * f + 2.0 * f * q + 2.0 * q * q;
        let energy = g * rate + r * (f * f + 4.0 * f * q + 6.0 * q * q);
        let log_scale = ln(self.global_cross_section_upper_bound.value()) - ln(self.scale) - end;
        [rate, energy]
            .into_iter()
            .zip(values)
            .map(|(polynomial, value)| {
                if value == 0.0 {
                    f64::INFINITY
                } else {
                    // Numerical outward allowance, not an interval certificate.
                    let bound = exp(log_scale + ln(polynomial) - ln(value));
                    (bound * (1.0 + 128.0 * f64::EPSILON)).max(f64::from_bits(1))
                }
            })
            .fold(0.0, f64::max)
    }
}

#[derive(Clone, Copy)]
struct Scales {
    thermal: f64,
    onset: f64,
    p: [f64; 2],
    k: [f64; 2],
}

#[derive(Clone, Copy, Default)]
struct Sum {
    value: f64,
    correction: f64,
}
impl Sum {
    fn add(&mut self, value: f64) {
        let corrected = value - self.correction;
        let next = self.value + corrected;
        self.correction = (next - self.value) - corrected;
        self.value = next;
    }
}
fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}

/// Retained Maxwellian moments. A common relative tail bound covers all three
/// positive moments, conditional on the explicitly asserted global envelope.
/// Numerical quadrature estimates do not include physical source uncertainty.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RadiativeRecombinationEvaluationV1 {
    /// Event coefficient alpha [m³/s], for number-density mass action.
    pub coefficient_m3_s: f64,
    /// Captured electron kinetic-energy coefficient L [J m³/s].
    pub captured_kinetic_energy_j_m3_s: f64,
    /// Emitted photon-energy coefficient chi*alpha+L [J m³/s].
    pub emitted_photon_energy_j_m3_s: f64,
    /// Largest numerical Gauss 8/16 difference relative to either moment.
    pub relative_quadrature_error_estimate: f64,
    /// Largest conditional omitted-tail bound relative to either moment.
    pub relative_tail_bound: f64,
}
impl RadiativeRecombinationEvaluationV1 {
    /// Mean electron kinetic energy captured per event [J], undefined if dark.
    pub fn mean_captured_energy_j(&self) -> Option<f64> {
        (self.coefficient_m3_s > 0.0)
            .then(|| self.captured_kinetic_energy_j_m3_s / self.coefficient_m3_s)
    }
}
