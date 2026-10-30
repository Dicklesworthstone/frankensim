//! Exact constant-coefficient exponential rolling paths and their physical work.
//!
//! The scalar path solves `q' = drive - decay*q` on a held interval. Signed
//! rate, displacement, squared-rate integral and integrated displacement are
//! kept separately. In particular, a mean rate cannot supply viscous heat and
//! an RMS rate cannot supply the signed moment impulse.
//!
//! Weak decay uses bounded analytic series, including cancellations in the
//! endpoint and first moments. The squared-rate integral is a positive Gram
//! form. Strong decay uses endpoint interpolation and its positive variance;
//! transient integrals are retained even when `decay*duration` overflows.
//! Arithmetic overflow or unresolved positivity refuses; no negative heat is
//! clamped. These are floating-point evaluations of the analytic model, not
//! interval-certified bounds on roundoff or the material approximation.

use core::fmt;

use fs_math::det;

use super::super::{
    RollingGeneralizedWork, RollingKinematics, RollingLossCheckpoint, RollingLossDissipation,
    RollingLossError, RollingLossLaw, RollingLossState, RollingLossStateKind, RollingLossStep,
    RollingLossUncertainty, RollingLossWrench, RollingPatchReceipt, RollingWorkOwnership,
    checked_add, checked_mul, validate_interface,
};
use super::ViscoelasticSphereResponse;
use crate::InterfaceSystemRef;

/// Invalid path arithmetic, original constitutive refusal, or cancellation.
#[derive(Debug, Clone, PartialEq)]
pub enum ExponentialPathError {
    Invalid(&'static str),
    Constitutive(RollingLossError),
    Cancelled,
}
impl From<RollingLossError> for ExponentialPathError {
    fn from(error: RollingLossError) -> Self {
        Self::Constitutive(error)
    }
}
impl fmt::Display for ExponentialPathError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "exponential rolling path: {self:?}")
    }
}
impl std::error::Error for ExponentialPathError {}

/// Immutable analytic path `q' = drive - decay*q` with its original moments.
/// Units of q are caller-owned; material integration requires rad/s.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ExponentialRatePath {
    initial_rate: f64,
    drive_rate_per_s: f64,
    decay_rate_per_s: f64,
    duration_s: f64,
    endpoint_rate: f64,
    rate_integral: f64,
    rate_squared_integral: f64,
    integrated_displacement: f64,
}

impl ExponentialRatePath {
    /// All supplied values are finite; decay and duration are nonnegative.
    /// Work is fixed at 32 analytic terms for `decay*duration <= 1`, with
    /// cancellation polled every eight terms, or constant work above that.
    pub fn new(
        initial_rate: f64,
        drive_rate_per_s: f64,
        decay_rate_per_s: f64,
        duration_s: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, ExponentialPathError> {
        poll(&mut cancelled)?;
        if [initial_rate, drive_rate_per_s, decay_rate_per_s, duration_s]
            .iter()
            .any(|x| !x.is_finite())
            || decay_rate_per_s < 0.0
            || duration_s < 0.0
        {
            return Err(invalid(
                "finite exponential path with nonnegative decay/duration",
            ));
        }
        let mut result = Self {
            initial_rate,
            drive_rate_per_s,
            decay_rate_per_s,
            duration_s,
            endpoint_rate: initial_rate,
            rate_integral: 0.0,
            rate_squared_integral: 0.0,
            integrated_displacement: 0.0,
        };
        if duration_s == 0.0 {
            return Ok(result);
        }
        let z = decay_rate_per_s * duration_s;
        let moments = if z <= 1.0 {
            weak_moments(
                initial_rate,
                drive_rate_per_s,
                duration_s,
                z,
                &mut cancelled,
            )?
        } else {
            strong_moments(
                initial_rate,
                drive_rate_per_s,
                decay_rate_per_s,
                duration_s,
                z,
            )?
        };
        result.endpoint_rate = moments[0];
        result.rate_integral = moments[1];
        result.rate_squared_integral = moments[2];
        result.integrated_displacement = moments[3];
        poll(&mut cancelled)?;
        Ok(result)
    }

    #[must_use]
    pub const fn initial_rate(&self) -> f64 {
        self.initial_rate
    }
    #[must_use]
    pub const fn drive_rate_per_s(&self) -> f64 {
        self.drive_rate_per_s
    }
    #[must_use]
    pub const fn decay_rate_per_s(&self) -> f64 {
        self.decay_rate_per_s
    }
    #[must_use]
    pub const fn duration_s(&self) -> f64 {
        self.duration_s
    }
    #[must_use]
    pub const fn endpoint_rate(&self) -> f64 {
        self.endpoint_rate
    }
    /// Signed `integral_0^duration q(t) dt`.
    #[must_use]
    pub const fn rate_integral(&self) -> f64 {
        self.rate_integral
    }
    /// Nonnegative `integral_0^duration q(t)^2 dt`.
    #[must_use]
    pub const fn rate_squared_integral(&self) -> f64 {
        self.rate_squared_integral
    }
    /// `integral_0^duration (duration-t)*q(t) dt`, useful for native positions.
    #[must_use]
    pub const fn integrated_displacement(&self) -> f64 {
        self.integrated_displacement
    }
}

/// Original material-owned path result. The wrench in `rolling` is the actual
/// ENDPOINT wrench, while heat, work and state contain the complete path integral.
#[derive(Debug, Clone, PartialEq)]
pub struct ViscoelasticSpherePathStep {
    pub path: ExponentialRatePath,
    pub rolling: RollingLossStep,
    /// Signed `-A*N*R * integral omega dt` [N m s].
    pub rolling_moment_impulse_nms: f64,
    pub response_before: ViscoelasticSphereResponse,
    pub response_after: ViscoelasticSphereResponse,
    pub temperature_k: f64,
    pub excitation_frequency_hz: f64,
    pub spin_rate_rad_s: f64,
}

impl RollingLossLaw {
    /// Integrate the original viscoelastic sphere law along an angular path.
    ///
    /// The patch, source card, temperature, frequency, axis and work owner are
    /// held. The exponential scalar rate is monotone, so its two endpoints
    /// bound every absolute speed used by the original quasistatic gates.
    /// Spin is an explicit held coordinate and produces no loss in this law.
    /// Other rolling/contour laws refuse this integration rather than inherit
    /// an unsupported exponential or quadratic-power assumption.
    ///
    /// The native owner restores accepted checkpoint lineage as for `advance`;
    /// this constitutive owner alone increments cumulative heat and produces the
    /// new original lineage-bound checkpoint. No caller-supplied heat enters.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_exponential_sphere_path(
        &self,
        patch: &RollingPatchReceipt,
        interface: &InterfaceSystemRef,
        path: &ExponentialRatePath,
        spin_rate_rad_s: f64,
        temperature_k: f64,
        excitation_frequency_hz: f64,
        ownership: &RollingWorkOwnership,
        state: &RollingLossState,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ViscoelasticSpherePathStep, ExponentialPathError> {
        poll(&mut cancelled)?;
        let Self::ViscoelasticSphere(card) = self else {
            return Err(invalid(
                "exponential path requires original viscoelastic sphere law",
            ));
        };
        patch.validate()?;
        validate_interface(interface)?;
        self.validate()?;
        ownership.validate_for(patch.patch_id(), self.required_channel())?;
        state.validate()?;
        let radius = card.parameters().radius_m;
        let kinematics = |omega| {
            RollingKinematics::new(
                checked_mul(radius, omega, "path_contour_speed_mps")?,
                omega,
                spin_rate_rad_s,
                temperature_k,
                excitation_frequency_hz,
                path.duration_s,
            )
        };
        let response_before = card.response(patch, interface, kinematics(path.initial_rate)?)?;
        poll(&mut cancelled)?;
        let response_after = card.response(patch, interface, kinematics(path.endpoint_rate)?)?;
        let coefficient = checked_mul(
            checked_mul(
                card.rolling_time_s(),
                radius,
                "path_material_rolling_damping",
            )?,
            patch.normal_load_n(),
            "path_material_rolling_damping",
        )?;
        if coefficient == 0.0 && card.rolling_time_s() > 0.0 && patch.normal_load_n() > 0.0 {
            return Err(invalid("unrepresentable original path rolling damping"));
        }
        let rolling_moment_impulse_nms = product(&[-coefficient, path.rate_integral])?;
        let heat = product(&[coefficient, path.rate_squared_integral])?;
        if heat < 0.0 {
            return Err(invalid("nonpassive original rolling path heat"));
        }
        let next_state = RollingLossState {
            cumulative_heat_j: checked_add(state.cumulative_heat_j, heat, "cumulative_heat_j")?,
            accepted_steps: state
                .accepted_steps
                .checked_add(1)
                .ok_or_else(|| invalid("path accepted-step counter"))?,
        };
        next_state.validate()?;
        let checkpoint = RollingLossCheckpoint {
            patch: patch.clone(),
            interface: interface.clone(),
            law: self.clone(),
            ownership: ownership.clone(),
            state: next_state.clone(),
        };
        let average_power = if path.duration_s == 0.0 {
            0.0
        } else {
            number(-heat / path.duration_s)?
        };
        let rolling = RollingLossStep {
            state: if response_before.rolling_moment_nm == 0.0
                && response_after.rolling_moment_nm == 0.0
            {
                RollingLossStateKind::Quiescent
            } else {
                RollingLossStateKind::RollingHysteresis
            },
            wrench: RollingLossWrench {
                contour_force_n: 0.0,
                rolling_moment_nm: response_after.rolling_moment_nm,
                spin_moment_nm: 0.0,
            },
            stored_energy_j: 0.0,
            storage_change_j: 0.0,
            dissipation: RollingLossDissipation {
                contour_deformation_heat_j: 0.0,
                rolling_hysteresis_heat_j: heat,
                total_heat_j: heat,
            },
            generalized_work: RollingGeneralizedWork {
                ownership: ownership.clone(),
                work_into_interface_j: heat,
                endpoint_body_power_w: -response_after.dissipated_power_w,
                reconstructed_body_power_w: average_power,
            },
            next_state,
            checkpoint,
            applicability: self.applicability_kind(),
            uncertainty: RollingLossUncertainty::ModelFormNoNumericalCertificate {
                curvature_approximation_authority: patch.curvature.approximation_authority(),
            },
            patch_authority: patch.authority,
            interface_authority: interface.provenance().authority(),
            card_authority: self.authority(),
        };
        poll(&mut cancelled)?;
        Ok(ViscoelasticSpherePathStep {
            path: *path,
            rolling,
            rolling_moment_impulse_nms,
            response_before,
            response_after,
            temperature_k,
            excitation_frequency_hz,
            spin_rate_rad_s,
        })
    }
}

// Taylor coefficients are integrated analytically, not by time sampling. The
// truncation after 32 terms is below 1e-26 on z<=1, including phi(2z); ordinary
// floating arithmetic dominates. Compensated summation retains weak-decay
// contributions when the zero-decay endpoint or first moment cancels exactly.
fn weak_moments(
    initial: f64,
    drive: f64,
    duration: f64,
    z: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<[f64; 4], ExponentialPathError> {
    let released = product(&[drive, duration])?;
    let scale = initial.abs().max(released.abs());
    if scale == 0.0 {
        return Ok([0.0; 4]);
    }
    let q = initial / scale;
    let b = released / scale;
    let mut endpoint = Sum::default();
    let mut first = Sum::default();
    let mut second = Sum::default();
    let mut aa = Sum::default();
    let mut ab = Sum::default();
    let mut bb = Sum::default();
    let mut term = 1.0;
    let mut two = 1.0;
    for k in 0..32 {
        if k % 8 == 0 {
            poll(cancelled)?;
        }
        let a = (k + 1) as f64;
        let c = (k + 2) as f64;
        let d = (k + 3) as f64;
        endpoint.add(term * (q + b / a));
        first.add(term * (q / a + b / (a * c)));
        second.add(term * (q / (a * c) + b / (a * c * d)));
        aa.add(term * two / a);
        ab.add(term * (2.0 * two - 1.0) / (a * c));
        bb.add(term * (4.0 * two - 2.0) / (a * c * d));
        term *= -z / a;
        two *= 2.0;
    }
    let aa = aa.value();
    let ab = ab.value();
    let remainder = if z == 0.0 {
        1.0 / 12.0
    } else {
        bb.value() - (ab / aa) * ab
    };
    if aa <= 0.0 || remainder <= 0.0 {
        return Err(invalid("unresolved exponential path Gram positivity"));
    }
    let projected = q + (ab / aa) * b;
    let norm = number(aa * projected * projected + remainder * b * b)?;
    Ok([
        product(&[scale, endpoint.value()])?,
        product(&[scale, duration, first.value()])?,
        product(&[scale, scale, duration, norm])?,
        product(&[scale, duration, duration, second.value()])?,
    ])
}

fn strong_moments(
    initial: f64,
    drive: f64,
    decay: f64,
    duration: f64,
    z: f64,
) -> Result<[f64; 4], ExponentialPathError> {
    let exponential = if z.is_infinite() { 0.0 } else { det::exp(-z) };
    let difference = 1.0 - exponential;
    let integral_exp = difference / decay;
    let endpoint = number(initial.mul_add(exponential, product(&[drive, integral_exp])?))?;
    // q(t) is a convex interpolation of its two endpoints. Coefficients below
    // are integrated in seconds, preserving a short initial transient even
    // when the dimensionless decay*time product is not representable.
    let first_initial = number((integral_exp - exponential * duration) / difference)?;
    let first_final = number(duration - first_initial)?;
    let first = number(product(&[initial, first_initial])? + product(&[endpoint, first_final])?)?;
    let variance = number(
        (0.5 * (1.0 + exponential) * integral_exp - (integral_exp / duration) * integral_exp)
            / (difference * difference),
    )?;
    if first_initial < 0.0 || first_final < 0.0 || variance <= 0.0 {
        return Err(invalid("unresolved exponential path variance"));
    }
    let scale = initial.abs().max(endpoint.abs());
    let rate_difference = if scale == 0.0 {
        0.0
    } else {
        endpoint / scale - initial / scale
    };
    let squared = number(
        product(&[first, first, 1.0 / duration])?
            + product(&[scale, scale, variance, rate_difference * rate_difference])?,
    )?;
    let second_initial = number(
        ((1.0 - integral_exp / duration) / decay - 0.5 * exponential * duration) / difference,
    )?;
    let second_final = number(0.5 * duration - second_initial)?;
    if second_initial < 0.0 || second_final < 0.0 {
        return Err(invalid("unresolved integrated exponential displacement"));
    }
    let second = number(
        product(&[initial, second_initial, duration])?
            + product(&[endpoint, second_final, duration])?,
    )?;
    Ok([endpoint, first, squared, second])
}

#[derive(Default)]
struct Sum {
    value: f64,
    correction: f64,
}
impl Sum {
    fn add(&mut self, term: f64) {
        let next = self.value + term;
        self.correction += if self.value.abs() >= term.abs() {
            (self.value - next) + term
        } else {
            (term - next) + self.value
        };
        self.value = next;
    }
    fn value(&self) -> f64 {
        self.value + self.correction
    }
}

// At most four factors. Pairing the smallest and largest magnitudes avoids
// overflowing q^2 before a small duration/decay coefficient can scale it, and
// avoids underflowing a small rate*time before a later large factor restores it.
fn product(factors: &[f64]) -> Result<f64, ExponentialPathError> {
    let mut work = [1.0; 4];
    let mut n = factors.len();
    let mut sign = 1.0;
    for (index, &factor) in factors.iter().enumerate() {
        number(factor)?;
        if factor == 0.0 {
            return Ok(0.0);
        }
        sign *= factor.signum();
        work[index] = factor.abs();
    }
    while n > 1 {
        work[..n].sort_by(f64::total_cmp);
        let value = number(work[0] * work[n - 1])?;
        work[0] = value;
        n -= 1;
    }
    number(sign * work[0])
}
fn number(value: f64) -> Result<f64, ExponentialPathError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid("unrepresentable exponential path arithmetic"))
    }
}
fn invalid(field: &'static str) -> ExponentialPathError {
    ExponentialPathError::Invalid(field)
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), ExponentialPathError> {
    if cancelled() {
        Err(ExponentialPathError::Cancelled)
    } else {
        Ok(())
    }
}
