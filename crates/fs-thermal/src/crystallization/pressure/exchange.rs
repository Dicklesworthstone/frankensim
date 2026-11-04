//! Implicit two-parcel contact on the original enthalpy and phase-history owners.
//!
//! Positive Q leaves the left parcel and enters the right parcel. Each candidate
//! advances BOTH accepted inputs from their original state with heats -Q,+Q,
//! using their existing adaptive crystallization kernels. The contact equation is
//! backward Euler: Q = G*dt*(T_left_end-T_right_end). G is the caller's admitted
//! total conductance [W/K] for this interval; no conductivity mixture rule,
//! contact geometry, material identity, or changing-interface law is inferred.
//!
//! The explicit finite heat bracket must cover a root and admit every evaluated
//! phase/thermal path. Neither heat direction nor uniqueness is inferred from
//! initial temperatures: adiabatic crystallization can change their ordering.
//! The inner phase adaptivity resolves its prescribed linear enthalpy path; it
//! does not estimate this backward-Euler contact step's time-discretization error.
//! A contact residual, a caloric error and an energy balance are distinct tests.

use super::super::{CrystallizationAdaptivity, CrystallizationError, check_cancelled};
use super::parcel::IsobaricCrystallizationParcel;
use super::{
    IsobaricCrystallizationControls, IsobaricCrystallizationStep, finite, invalid, product, sum,
};

/// Original per-body accuracy and bounded phase work. The two sides may differ.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct IsobaricExchangeSideControls {
    pub caloric: IsobaricCrystallizationControls,
    pub adaptivity: CrystallizationAdaptivity,
}

impl IsobaricExchangeSideControls {
    pub fn validate(self) -> Result<(), CrystallizationError> {
        self.caloric.validate()?;
        self.adaptivity.validate()
    }
}

/// Explicit contact solve bounds, independent of each body's phase solver.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct IsobaricExchangeControls {
    /// Signed total heat interval [lower,upper], J; positive means left to right.
    /// A singleton is permitted only if that candidate satisfies the residual.
    pub heat_bracket_j: [f64; 2],
    /// All candidate evaluations, including both bracket endpoints. Each calls
    /// at most two original parcel advances with their own bounded trial caps.
    pub maximum_evaluations: usize,
    /// Absolute backward-Euler contact residual, J. Not a time-error estimate.
    pub residual_tolerance_j: f64,
    /// Combined noncancelling caloric/work/arithmetic defects, J. The internal
    /// contact-law residual is reported separately and is never added as heat.
    pub energy_tolerance_j: f64,
}

impl IsobaricExchangeControls {
    pub fn validate(self) -> Result<(), CrystallizationError> {
        let [lower, upper] = self.heat_bracket_j;
        if !lower.is_finite()
            || !upper.is_finite()
            || lower > upper
            || self.maximum_evaluations == 0
            || !self.residual_tolerance_j.is_finite()
            || self.residual_tolerance_j < 0.0
            || !self.energy_tolerance_j.is_finite()
            || self.energy_tolerance_j < 0.0
        {
            return Err(invalid("invalid bounded isobaric exchange controls"));
        }
        Ok(())
    }
}

/// One accepted simultaneous proposal; internal heat cancels across the pair.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct IsobaricExchangeStep {
    /// Actual represented increment of both initially equal parcel clocks.
    pub duration_s: f64,
    pub conductance_w_k: f64,
    pub heat_left_to_right_j: f64,
    /// Signed Q-G*dt*(T_left_end-T_right_end), J.
    pub residual_j: f64,
    pub evaluated_candidates: usize,
    pub combined_enthalpy_change_j: f64,
    pub combined_internal_energy_change_j: f64,
    pub combined_boundary_work_j: f64,
    /// Combined delta_U+work, with zero externally supplied heat.
    pub energy_residual_j: f64,
    /// Noncancelling original per-body defects plus combination roundoff.
    pub absolute_energy_residual_j: f64,
    pub left: IsobaricCrystallizationStep,
    pub right: IsobaricCrystallizationStep,
}

/// Solve one backward-Euler contact step without mutating either accepted input.
///
/// A safeguarded scalar bracket only solves for shared Q; all phase evolution,
/// remelting, density, pressure work and enthalpy inversion retain their original
/// owners. Distinct source lifetimes and per-body controls permit unlike media.
/// Endpoint/source failures and cancellation propagate unchanged. A bracket
/// failure, exhausted evaluations, unrepresentable progress or final energy
/// refusal publishes neither body. No partial interval is returned.
#[allow(clippy::too_many_arguments, clippy::type_complexity)]
pub fn exchange_heat<'ll, 'ls, 'rl, 'rs>(
    left: &IsobaricCrystallizationParcel<'ll, 'ls>,
    right: &IsobaricCrystallizationParcel<'rl, 'rs>,
    conductance_w_k: f64,
    duration_s: f64,
    left_controls: IsobaricExchangeSideControls,
    right_controls: IsobaricExchangeSideControls,
    controls: IsobaricExchangeControls,
    mut cancelled: impl FnMut() -> bool,
) -> Result<
    (
        IsobaricCrystallizationParcel<'ll, 'ls>,
        IsobaricCrystallizationParcel<'rl, 'rs>,
        IsobaricExchangeStep,
    ),
    CrystallizationError,
> {
    check_cancelled(&mut cancelled)?;
    controls.validate()?;
    left_controls.validate()?;
    right_controls.validate()?;
    if !conductance_w_k.is_finite() || conductance_w_k < 0.0 {
        return Err(invalid(
            "isobaric exchange conductance must be finite and nonnegative",
        ));
    }
    if !duration_s.is_finite() || duration_s <= 0.0 || left.time_s() != right.time_s() {
        return Err(invalid(
            "isobaric exchange requires positive duration and equal parcel clocks",
        ));
    }
    let final_time = finite(left.time_s() + duration_s)?;
    let represented_duration = final_time - left.time_s();
    if represented_duration <= 0.0 {
        return Err(invalid(
            "isobaric exchange elapsed-time increment is unrepresentable",
        ));
    }
    let conductance_time = product(conductance_w_k, represented_duration)?;
    let mut evaluate = |heat: f64| -> Result<Candidate<'ll, 'ls, 'rl, 'rs>, CrystallizationError> {
        check_cancelled(&mut cancelled)?;
        let (next_left, left_step) = left.advance_heat(
            -heat,
            duration_s,
            left_controls.caloric,
            left_controls.adaptivity,
            &mut cancelled,
        )?;
        let (next_right, right_step) = right.advance_heat(
            heat,
            duration_s,
            right_controls.caloric,
            right_controls.adaptivity,
            &mut cancelled,
        )?;
        if next_left.time_s() != final_time
            || next_right.time_s() != final_time
            || left_step.duration_s != represented_duration
            || right_step.duration_s != represented_duration
        {
            return Err(invalid("isobaric exchange candidate clocks disagree"));
        }
        let difference =
            finite(next_left.state().temperature_k() - next_right.state().temperature_k())?;
        let residual = sum(&[heat, -product(conductance_time, difference)?])?;
        check_cancelled(&mut cancelled)?;
        Ok(Candidate {
            left: next_left,
            right: next_right,
            left_step,
            right_step,
            heat,
            residual,
        })
    };
    let [mut lower, mut upper] = controls.heat_bracket_j;
    // No-conductance still advances both phase histories adiabatically.
    if conductance_w_k == 0.0 {
        if lower > 0.0 || upper < 0.0 {
            return Err(invalid(
                "zero-conductance exchange bracket must contain zero heat",
            ));
        }
        let candidate = evaluate(0.0)?;
        return candidate.finish(conductance_w_k, 1, controls, &mut cancelled);
    }
    let lower_candidate = evaluate(lower)?;
    if lower_candidate.residual.abs() <= controls.residual_tolerance_j {
        return lower_candidate.finish(conductance_w_k, 1, controls, &mut cancelled);
    }
    if lower == upper {
        return Err(invalid(
            "isobaric exchange heat bracket does not enclose a root",
        ));
    }
    if controls.maximum_evaluations == 1 {
        return Err(CrystallizationError::NonConverged);
    }
    let upper_candidate = evaluate(upper)?;
    if upper_candidate.residual.abs() <= controls.residual_tolerance_j {
        return upper_candidate.finish(conductance_w_k, 2, controls, &mut cancelled);
    }
    let mut lower_residual = lower_candidate.residual;
    let mut upper_residual = upper_candidate.residual;
    if lower_residual.is_sign_positive() == upper_residual.is_sign_positive() {
        return Err(invalid(
            "isobaric exchange heat bracket does not enclose a root",
        ));
    }
    for evaluations in 3..=controls.maximum_evaluations {
        // Alternate interpolation and bisection so a curved residual cannot
        // indefinitely pin one endpoint. Scaling avoids overflowing |f0|+|f1|.
        let midpoint = midpoint(lower, upper);
        let proposal = if evaluations % 2 == 1 {
            let scale = lower_residual.abs().max(upper_residual.abs());
            let a = lower_residual.abs() / scale;
            let b = upper_residual.abs() / scale;
            let weight = a / (a + b);
            let width = upper - lower;
            let interpolated = if width.is_finite() {
                lower + weight * width
            } else {
                (1.0 - weight) * lower + weight * upper
            };
            if interpolated.is_finite() && interpolated > lower && interpolated < upper {
                interpolated
            } else {
                midpoint
            }
        } else {
            midpoint
        };
        if proposal <= lower || proposal >= upper || !proposal.is_finite() {
            return Err(CrystallizationError::NonConverged);
        }
        let candidate = evaluate(proposal)?;
        if candidate.residual.abs() <= controls.residual_tolerance_j {
            return candidate.finish(conductance_w_k, evaluations, controls, &mut cancelled);
        }
        if candidate.residual.is_sign_positive() == lower_residual.is_sign_positive() {
            lower = proposal;
            lower_residual = candidate.residual;
        } else {
            upper = proposal;
            upper_residual = candidate.residual;
        }
    }
    Err(CrystallizationError::NonConverged)
}

struct Candidate<'ll, 'ls, 'rl, 'rs> {
    left: IsobaricCrystallizationParcel<'ll, 'ls>,
    right: IsobaricCrystallizationParcel<'rl, 'rs>,
    left_step: IsobaricCrystallizationStep,
    right_step: IsobaricCrystallizationStep,
    heat: f64,
    residual: f64,
}

impl<'ll, 'ls, 'rl, 'rs> Candidate<'ll, 'ls, 'rl, 'rs> {
    #[allow(clippy::type_complexity)]
    fn finish(
        self,
        conductance_w_k: f64,
        evaluated_candidates: usize,
        controls: IsobaricExchangeControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<
        (
            IsobaricCrystallizationParcel<'ll, 'ls>,
            IsobaricCrystallizationParcel<'rl, 'rs>,
            IsobaricExchangeStep,
        ),
        CrystallizationError,
    > {
        let left = self.left_step;
        let right = self.right_step;
        let combined_enthalpy_change_j = sum(&[left.enthalpy_change_j, right.enthalpy_change_j])?;
        let combined_internal_energy_change_j = sum(&[
            left.internal_energy_change_j,
            right.internal_energy_change_j,
        ])?;
        let combined_boundary_work_j = sum(&[left.boundary_work_j, right.boundary_work_j])?;
        let enthalpy_roundoff = sum(&[
            combined_enthalpy_change_j,
            -left.enthalpy_change_j,
            -right.enthalpy_change_j,
        ])?;
        let internal_roundoff = sum(&[
            combined_internal_energy_change_j,
            -left.internal_energy_change_j,
            -right.internal_energy_change_j,
        ])?;
        let work_roundoff = sum(&[
            combined_boundary_work_j,
            -left.boundary_work_j,
            -right.boundary_work_j,
        ])?;
        let energy_residual_j =
            sum(&[combined_internal_energy_change_j, combined_boundary_work_j])?;
        let absolute_energy_residual_j = sum(&[
            left.absolute_energy_residual_j,
            right.absolute_energy_residual_j,
            enthalpy_roundoff.abs(),
            internal_roundoff.abs(),
            work_roundoff.abs(),
        ])?
        .max(energy_residual_j.abs())
        .max(combined_enthalpy_change_j.abs());
        if absolute_energy_residual_j > controls.energy_tolerance_j {
            return Err(invalid(
                "combined isobaric exchange energy residual exceeds its tolerance",
            ));
        }
        check_cancelled(cancelled)?;
        Ok((
            self.left,
            self.right,
            IsobaricExchangeStep {
                duration_s: left.duration_s,
                conductance_w_k,
                heat_left_to_right_j: self.heat,
                residual_j: self.residual,
                evaluated_candidates,
                combined_enthalpy_change_j,
                combined_internal_energy_change_j,
                combined_boundary_work_j,
                energy_residual_j,
                absolute_energy_residual_j,
                left,
                right,
            },
        ))
    }
}

fn midpoint(lower: f64, upper: f64) -> f64 {
    let width = upper - lower;
    if width.is_finite() {
        lower + 0.5 * width
    } else {
        0.5 * lower + 0.5 * upper
    }
}

#[cfg(test)]
#[path = "exchange/tests.rs"]
mod tests;
