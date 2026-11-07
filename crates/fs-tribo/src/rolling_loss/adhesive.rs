//! Fixed-axis adhesive rolling with elastic memory and irreversible yield.
//!
//! Dominik and Tielens (1995), DOI 10.1080/01418619508243800, equation (25),
//! give M = -4 Fc (a/a0)^(3/2) xi, with xi = R theta and Fc = 3 pi W R / 2.
//! Here W is WORK of adhesion, twice their single-surface energy gamma.
//! The retained spring has k = 4 Fc (a/a0)^(3/2)/R and U = k xi^2/2.
//! This owner explicitly selects a constant-moment, elastic/perfectly plastic
//! continuation at an independently supplied critical displacement. It is not
//! a universal microscopic law for contact migration or a calibrated material.
//!
//! The normal owner supplies one attached JKR radius and its applicability.
//! Normal area, source conditions and the rolling axis remain fixed. Rolling
//! storage must not be discarded on detachment or changed-area rebinding; those
//! operations, arbitrary vector rolling, spin, sliding and heat allocation are
//! outside this owner. One signed scalar increment is a monotone physical path.
//! Its exact elastic/yield split is allocation-free and independent of step size.

use core::f64::consts::PI;
use fs_math::{dd::Dd, det};

/// Explicit SI parameters; the normal owner establishes the JKR approximation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AdhesiveRollingParameters {
    /// Reduced sphere radius; the sphere radius for a rigid plane [m].
    pub effective_radius_m: f64,
    /// Work of adhesion W, not a single-surface energy [J/m2].
    pub work_of_adhesion_j_m2: f64,
    pub reduced_modulus_pa: f64,
    /// Actual attached radius, including at zero or tensile resultant load [m].
    pub contact_radius_m: f64,
    /// Independent yield shift; there is no hidden atomic-length default [m].
    pub critical_rolling_displacement_m: f64,
    /// Caller small-shift ceiling xi_c/a, strictly between zero and one.
    pub maximum_critical_displacement_to_contact_radius: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub enum AdhesiveRollingError {
    InvalidInput { field: &'static str },
    OutsideSmallShift { ratio: f64, limit: f64 },
    StateMismatch,
    Unrepresentable { field: &'static str },
}
impl core::fmt::Display for AdhesiveRollingError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "adhesive rolling refused: {self:?}")
    }
}
impl std::error::Error for AdhesiveRollingError {}

/// Immutable history bound to all original physical parameters and limits.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AdhesiveRollingState {
    parameters: AdhesiveRollingParameters,
    elastic_fraction: Dd,
    angle: Dd,
}
impl AdhesiveRollingState {
    pub fn elastic_shift_m(self) -> f64 {
        self.parameters.critical_rolling_displacement_m * self.elastic_fraction.to_f64()
    }
    /// Signed total angle, rounded once for observation; low bits remain retained.
    pub const fn accumulated_angle_rad(self) -> f64 {
        self.angle.to_f64()
    }
}

/// Work on the original monotone scalar path; resisting moments oppose loading.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AdhesiveRollingStep {
    pub state: AdhesiveRollingState,
    pub endpoint_moment_n_m: f64,
    /// Integral of resisting moment divided by angle; endpoint moment at zero step.
    pub mean_moment_n_m: f64,
    pub storage_before_j: f64,
    pub storage_after_j: f64,
    /// Factored change, retaining work even when rounded endpoint energies agree.
    pub storage_change_j: f64,
    pub irreversible_loss_j: f64,
    /// Applied moment work, positive into the contact: change of storage + loss.
    pub applied_work_j: f64,
    pub plastic_rolling_increment_m: f64,
}

/// Frozen normal patch with exact scalar return and separate storage/loss.
#[derive(Clone, Debug, PartialEq)]
pub struct AdhesiveRollingContact {
    parameters: AdhesiveRollingParameters,
    stiffness: f64,
    critical_moment: f64,
    critical_angle: f64,
    yield_energy: f64,
    equilibrium_radius: f64,
}
impl AdhesiveRollingContact {
    pub fn new(p: AdhesiveRollingParameters) -> Result<Self, AdhesiveRollingError> {
        for (field, value) in [
            ("effective radius", p.effective_radius_m),
            ("work of adhesion", p.work_of_adhesion_j_m2),
            ("reduced modulus", p.reduced_modulus_pa),
            ("attached radius", p.contact_radius_m),
            ("critical displacement", p.critical_rolling_displacement_m),
            (
                "small-shift limit",
                p.maximum_critical_displacement_to_contact_radius,
            ),
        ] {
            if !(value.is_finite() && value > 0.0) {
                return Err(AdhesiveRollingError::InvalidInput { field });
            }
        }
        let limit = p.maximum_critical_displacement_to_contact_radius;
        if limit >= 1.0 {
            return Err(AdhesiveRollingError::InvalidInput {
                field: "small-shift limit",
            });
        }
        let ratio = p.critical_rolling_displacement_m / p.contact_radius_m;
        if !ratio.is_finite() || ratio > limit {
            return Err(AdhesiveRollingError::OutsideSmallShift { ratio, limit });
        }
        // Eliminate a0^3 = 9 pi W R^2/(2 E*) before evaluating stiffness.
        // Separate square roots avoid squaring a radius or multiplying W*E*a.
        let stiffness = positive(
            (p.contact_radius_m / p.effective_radius_m)
                * det::sqrt(8.0 * PI)
                * det::sqrt(p.work_of_adhesion_j_m2)
                * det::sqrt(p.reduced_modulus_pa)
                * det::sqrt(p.contact_radius_m),
            "stiffness",
        )?;
        let critical_force =
            positive(stiffness * p.critical_rolling_displacement_m, "yield force")?;
        let critical_moment = positive(critical_force * p.effective_radius_m, "yield moment")?;
        let critical_angle = positive(
            p.critical_rolling_displacement_m / p.effective_radius_m,
            "yield angle",
        )?;
        let yield_energy = positive(
            (0.5 * critical_force) * p.critical_rolling_displacement_m,
            "yield storage",
        )?;
        let root = positive(
            (p.work_of_adhesion_j_m2 / p.reduced_modulus_pa) / p.effective_radius_m,
            "equilibrium radius ratio cubed",
        )?;
        let equilibrium_radius = positive(
            p.effective_radius_m * det::pow(4.5 * PI * root, 1.0 / 3.0),
            "equilibrium radius",
        )?;
        Ok(Self {
            parameters: p,
            stiffness,
            critical_moment,
            critical_angle,
            yield_energy,
            equilibrium_radius,
        })
    }
    pub const fn parameters(&self) -> AdhesiveRollingParameters {
        self.parameters
    }
    pub const fn stiffness_n_per_m(&self) -> f64 {
        self.stiffness
    }
    pub const fn critical_moment_n_m(&self) -> f64 {
        self.critical_moment
    }
    pub const fn critical_rolling_angle_rad(&self) -> f64 {
        self.critical_angle
    }
    pub const fn equilibrium_contact_radius_m(&self) -> f64 {
        self.equilibrium_radius
    }
    pub const fn initial_state(&self) -> AdhesiveRollingState {
        AdhesiveRollingState {
            parameters: self.parameters,
            elastic_fraction: Dd::ZERO,
            angle: Dd::ZERO,
        }
    }
    pub fn stored_energy_j(
        &self,
        state: &AdhesiveRollingState,
    ) -> Result<f64, AdhesiveRollingError> {
        self.admit(state)?;
        finite(
            (Dd::from_f64(self.yield_energy) * state.elastic_fraction * state.elastic_fraction)
                .to_f64(),
            "stored energy",
        )
    }
    fn admit(&self, state: &AdhesiveRollingState) -> Result<(), AdhesiveRollingError> {
        if state.parameters != self.parameters {
            return Err(AdhesiveRollingError::StateMismatch);
        }
        Ok(())
    }

    /// Advance an explicitly prescribed fixed-axis rolling increment [rad].
    /// Sub-ULP history increments are retained with the existing double-double
    /// arithmetic. Failure returns no successor and changes no accepted state.
    pub fn advance(
        &self,
        old: &AdhesiveRollingState,
        angle_increment_rad: f64,
    ) -> Result<AdhesiveRollingStep, AdhesiveRollingError> {
        self.admit(old)?;
        if !angle_increment_rad.is_finite() {
            return Err(AdhesiveRollingError::InvalidInput {
                field: "angle increment",
            });
        }
        let angle = old.angle + Dd::from_f64(angle_increment_rad);
        finite(angle.hi, "accumulated angle")?;
        finite(angle.lo, "accumulated angle residual")?;
        if angle_increment_rad != 0.0 && angle == old.angle {
            return Err(AdhesiveRollingError::Unrepresentable {
                field: "accumulated angle increment",
            });
        }
        let increment = Dd::from_f64(angle_increment_rad) / Dd::from_f64(self.critical_angle);
        finite(increment.hi, "normalized rolling increment")?;
        finite(increment.lo, "normalized rolling increment residual")?;
        if angle_increment_rad != 0.0 && increment == Dd::ZERO {
            return Err(AdhesiveRollingError::Unrepresentable {
                field: "rolling increment",
            });
        }
        let trial = old.elastic_fraction + increment;
        if increment != Dd::ZERO && trial == old.elastic_fraction {
            return Err(AdhesiveRollingError::Unrepresentable {
                field: "elastic history increment",
            });
        }
        let mut fraction = trial;
        if Dd::ONE.lt(trial) {
            fraction = Dd::ONE;
        }
        if trial.lt(-Dd::ONE) {
            fraction = -Dd::ONE;
        }
        let plastic = trial - fraction;
        let change = fraction - old.elastic_fraction;
        let energy = Dd::from_f64(self.yield_energy);
        let storage_change = energy * change * (fraction + old.elastic_fraction);
        let loss = energy * Dd::from_f64(2.0) * plastic.abs();
        let work = storage_change + loss;
        // Pure elasticity uses its analytic midpoint, including zero/tiny steps.
        // On a yielding path, the nonnegative factored numerator avoids loss of
        // the small plastic work after a nearly complete elastic reversal.
        let mean_fraction = if plastic == Dd::ZERO {
            (old.elastic_fraction + fraction) * Dd::from_f64(0.5)
        } else {
            (change * (fraction + old.elastic_fraction) * Dd::from_f64(0.5) + plastic.abs())
                / increment
        };
        let state = AdhesiveRollingState {
            parameters: self.parameters,
            elastic_fraction: fraction,
            angle,
        };
        let step = AdhesiveRollingStep {
            state,
            endpoint_moment_n_m: -self.critical_moment * fraction.to_f64(),
            mean_moment_n_m: -self.critical_moment * mean_fraction.to_f64(),
            storage_before_j: self.stored_energy_j(old)?,
            storage_after_j: self.stored_energy_j(&state)?,
            storage_change_j: storage_change.to_f64(),
            irreversible_loss_j: loss.to_f64(),
            applied_work_j: work.to_f64(),
            plastic_rolling_increment_m: self.parameters.critical_rolling_displacement_m
                * plastic.to_f64(),
        };
        for (field, value) in [
            ("endpoint moment", step.endpoint_moment_n_m),
            ("mean moment", step.mean_moment_n_m),
            ("storage change", step.storage_change_j),
            ("rolling loss", step.irreversible_loss_j),
            ("applied work", step.applied_work_j),
            ("plastic increment", step.plastic_rolling_increment_m),
        ] {
            finite(value, field)?;
        }
        Ok(step)
    }
}
fn finite(value: f64, field: &'static str) -> Result<f64, AdhesiveRollingError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(AdhesiveRollingError::Unrepresentable { field })
    }
}
fn positive(value: f64, field: &'static str) -> Result<f64, AdhesiveRollingError> {
    if value.is_finite() && value > 0.0 {
        Ok(value)
    } else {
        Err(AdhesiveRollingError::Unrepresentable { field })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn parameters() -> AdhesiveRollingParameters {
        let r = 2e-6;
        let w = 0.08;
        let e = 2e9;
        AdhesiveRollingParameters {
            effective_radius_m: r,
            work_of_adhesion_j_m2: w,
            reduced_modulus_pa: e,
            contact_radius_m: (4.5 * PI * w * r * r / e).cbrt(),
            critical_rolling_displacement_m: 2e-10,
            maximum_critical_displacement_to_contact_radius: 0.01,
        }
    }
    fn close(a: f64, b: f64) {
        assert!(
            (a - b).abs() <= 3e-13 * a.abs().max(b.abs()).max(1e-100),
            "{a} != {b}"
        );
    }
    fn balance(step: &AdhesiveRollingStep, angle: f64) {
        close(
            step.applied_work_j,
            step.storage_change_j + step.irreversible_loss_j,
        );
        close(step.applied_work_j, -step.mean_moment_n_m * angle);
        assert!(step.irreversible_loss_j >= 0.0);
    }
    #[test]
    fn g1_elastic_energy_and_primary_zero_load_stiffness() {
        let p = parameters();
        let model = AdhesiveRollingContact::new(p).unwrap();
        close(model.equilibrium_contact_radius_m(), p.contact_radius_m);
        close(
            model.stiffness_n_per_m(),
            6.0 * PI * p.work_of_adhesion_j_m2,
        );
        for sign in [-1.0, 1.0] {
            let angle = sign * 0.4 * model.critical_rolling_angle_rad();
            let step = model.advance(&model.initial_state(), angle).unwrap();
            close(
                step.endpoint_moment_n_m,
                -sign * 0.4 * model.critical_moment_n_m(),
            );
            close(step.mean_moment_n_m, 0.5 * step.endpoint_moment_n_m);
            close(
                step.storage_after_j,
                0.5 * model.stiffness_n_per_m() * (0.4 * p.critical_rolling_displacement_m).powi(2),
            );
            assert_eq!(step.irreversible_loss_j, 0.0);
            balance(&step, angle);
            let reverse = model.advance(&step.state, -angle).unwrap();
            assert_eq!(reverse.state.elastic_shift_m(), 0.0);
            close(reverse.applied_work_j, -step.applied_work_j);
            assert_eq!(reverse.irreversible_loss_j, 0.0);
        }
    }
    #[test]
    fn g1_yield_plateau_reversal_and_exact_work() {
        let model = AdhesiveRollingContact::new(parameters()).unwrap();
        let critical = model.critical_rolling_angle_rad();
        let forward = model
            .advance(&model.initial_state(), 3.0 * critical)
            .unwrap();
        close(forward.endpoint_moment_n_m, -model.critical_moment_n_m());
        close(
            forward.irreversible_loss_j,
            2.0 * critical * model.critical_moment_n_m(),
        );
        balance(&forward, 3.0 * critical);
        let unload = model.advance(&forward.state, -critical).unwrap();
        assert_eq!(unload.irreversible_loss_j, 0.0);
        close(unload.applied_work_j, -forward.storage_after_j);
        let reverse = model.advance(&forward.state, -2.5 * critical).unwrap();
        close(reverse.endpoint_moment_n_m, model.critical_moment_n_m());
        close(
            reverse.irreversible_loss_j,
            0.5 * critical * model.critical_moment_n_m(),
        );
        close(reverse.storage_change_j, 0.0);
        balance(&reverse, -2.5 * critical);
    }
    #[test]
    fn g3_monotone_subdivision_preserves_history_work_and_loss() {
        let model = AdhesiveRollingContact::new(parameters()).unwrap();
        let d = 0.125 * model.critical_rolling_angle_rad();
        let one = model.advance(&model.initial_state(), 32.0 * d).unwrap();
        let mut state = model.initial_state();
        let (mut work, mut loss) = (0.0, 0.0);
        for _ in 0..32 {
            let step = model.advance(&state, d).unwrap();
            state = step.state;
            work += step.applied_work_j;
            loss += step.irreversible_loss_j;
        }
        close(state.elastic_shift_m(), one.state.elastic_shift_m());
        close(
            state.accumulated_angle_rad(),
            one.state.accumulated_angle_rad(),
        );
        close(work, one.applied_work_j);
        close(loss, one.irreversible_loss_j);
    }
    #[test]
    fn g1_tensile_pull_off_retains_half_zero_load_rolling_stiffness() {
        let mut p = parameters();
        let original = AdhesiveRollingContact::new(p).unwrap();
        p.contact_radius_m /= 4.0_f64.cbrt();
        let pull_off = AdhesiveRollingContact::new(p).unwrap();
        close(
            pull_off.stiffness_n_per_m(),
            0.5 * original.stiffness_n_per_m(),
        );
        assert!(pull_off.critical_moment_n_m() > 0.0);
    }
    #[test]
    fn g3_zero_step_and_sub_ulp_increments_preserve_elastic_work() {
        let model = AdhesiveRollingContact::new(parameters()).unwrap();
        let critical = model.critical_rolling_angle_rad();
        let start = model
            .advance(&model.initial_state(), 0.5 * critical)
            .unwrap()
            .state;
        let zero = model.advance(&start, 0.0).unwrap();
        assert_eq!(zero.state, start);
        assert_eq!(zero.mean_moment_n_m, zero.endpoint_moment_n_m);
        assert_eq!(zero.applied_work_j, 0.0);
        let tiny = critical * 1e-20;
        let step = model.advance(&start, tiny).unwrap();
        assert_ne!(step.state, start);
        assert_eq!(step.storage_after_j, step.storage_before_j);
        assert!(step.storage_change_j > 0.0);
        balance(&step, tiny);
        assert!(model.advance(&step.state, critical * 1e-60).is_err());
        let undone = model.advance(&step.state, -tiny).unwrap();
        assert_eq!(undone.state.elastic_shift_m(), start.elastic_shift_m());
        assert_eq!(
            undone.state.accumulated_angle_rad(),
            start.accumulated_angle_rad()
        );
        assert!(
            (undone.state.elastic_fraction - start.elastic_fraction)
                .abs()
                .to_f64()
                < 1e-30
        );
        close(undone.applied_work_j, -step.applied_work_j);
    }
    #[test]
    fn g0_original_state_small_shift_and_nonfinite_refusals() {
        let mut p = parameters();
        let model = AdhesiveRollingContact::new(p).unwrap();
        let state = model.initial_state();
        p.work_of_adhesion_j_m2 *= 2.0;
        let other = AdhesiveRollingContact::new(p).unwrap();
        assert_eq!(
            other.advance(&state, 0.0),
            Err(AdhesiveRollingError::StateMismatch)
        );
        p.critical_rolling_displacement_m = p.contact_radius_m * 0.1;
        assert!(matches!(
            AdhesiveRollingContact::new(p),
            Err(AdhesiveRollingError::OutsideSmallShift { .. })
        ));
        for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert!(model.advance(&state, bad).is_err());
            let mut p = parameters();
            p.reduced_modulus_pa = bad;
            assert!(AdhesiveRollingContact::new(p).is_err());
        }
        assert_eq!(state, model.initial_state());
        assert!(model.advance(&state, f64::MAX).is_err());
        let mut p = parameters();
        p.contact_radius_m = 0.0;
        assert!(AdhesiveRollingContact::new(p).is_err());
    }
}
