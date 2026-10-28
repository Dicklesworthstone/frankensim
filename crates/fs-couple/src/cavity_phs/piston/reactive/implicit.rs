//! Implicit homogeneous reaction extent on the original conservative gas chart.
//!
//! Backward Euler solves `xi - dt V r(n0 + nu xi, U0, V) = 0` in moles.
//! Every trial reuses the original exact stoichiometry, source-bound rate and
//! detailed-balance reverse closure. The original gas owner reconstructs its
//! composition, T, p, transport and acoustic response at the SAME mass and
//! total internal energy; formation enthalpy is never added again as heat.
//!
//! This scalar rung admits one reaction, including inerts and initially absent
//! products. It is first order and can use the existing dyadic subcycling error
//! estimate. It is not a coupled implicit mechanism, equilibrium projection,
//! ignition-time accuracy certificate or unconditional nonlinear stability
//! theorem. Source support is required at all queried bracket states. The
//! accepted state must independently pass the existing entropy/balance gates.
//! Use [`coupled`] for the simultaneous multi-reaction endpoint equations.

use super::*;
use crate::cavity_phs::piston::thermal::feedback::{self, ScalarError};

/// Simultaneous endpoint extents on the native Newton--Krylov owner.
pub mod coupled;

/// Opt-in endpoint integration, retained across all existing gas consumers.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ImplicitReactionSolve {
    /// Signed extent range [mol], containing zero. It is intersected with the
    /// CURRENT donor-consumption allowance on each step; species are never
    /// clipped. Use this range to restrict the caloric/rate source domain.
    pub extent_bracket_mol: [f64; 2],
    /// Maximum absolute `xi - dt V r_end` [mol], independently of species and
    /// energy reconstruction tolerances. A small bracket alone never accepts.
    pub extent_tolerance_mol: f64,
    /// Interior safeguarded secant/bisection trials in 1..=128, plus at most
    /// two endpoint evaluations. Each trial retains the gas solver's own cap.
    pub maximum_iterations: usize,
}
impl ImplicitReactionSolve {
    fn validate(self) -> Result<(), GasReactionError> {
        let [lo, hi] = self.extent_bracket_mol;
        if !lo.is_finite() || !hi.is_finite() || lo >= hi || lo > 0.0 || hi < 0.0
            || !positive(self.extent_tolerance_mol)
            || !(1..=128).contains(&self.maximum_iterations)
        {
            return Err(invalid("implicit reaction needs a signed molar bracket, positive extent tolerance and 1..=128 iterations"));
        }
        Ok(())
    }
}

/// Actual nonlinear endpoint equation. The surrounding original reaction
/// report retains the endpoint rate, material source, species and energy.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ImplicitReactionReport {
    pub duration_s: f64,
    pub admitted_bracket_mol: [f64; 2],
    pub residual_mol: f64,
    pub trials: usize,
}

impl GasReaction {
    /// Select a conservative backward-Euler extent instead of explicit Euler.
    /// This choice reaches parcel, piston, nozzle, diffusion and finite thermal
    /// network consumers through their existing `parcel_candidate` call.
    /// Existing optional subcycling compares this method's actual trajectories.
    /// Binding a mechanism with multiple channels is refused at construction;
    /// no channel is silently omitted or advanced by a different time policy.
    pub fn with_implicit(mut self, policy: ImplicitReactionSolve) -> Result<Self, GasReactionError> {
        policy.validate()?;
        if self.implicit.is_some() || self.implicit_mechanism.is_some() || !self.additional_channels.is_empty() {
            return Err(invalid("implicit extent requires one reaction without an existing implicit policy"));
        }
        self.implicit = Some(policy);
        Ok(self)
    }

    #[must_use]
    pub const fn implicit(&self) -> Option<ImplicitReactionSolve> { self.implicit }
}

struct Trial {
    gas: ClosedGasParcel,
    actual: Vec<f64>,
    changes: Vec<f64>,
    expected: Vec<f64>,
    consumed: f64,
    extent_mol: f64,
    report: ImplicitReactionReport,
}

fn equation_residual(
    extent: f64, value: f64, volume_m3: f64, duration_s: f64,
) -> Result<f64, GasReactionError> {
    let predicted = (value * volume_m3) * duration_s;
    let residual = extent - predicted;
    if !predicted.is_finite() || (value != 0.0 && predicted == 0.0) || !residual.is_finite() {
        return Err(invalid("implicit reaction endpoint equation is unrepresentable"));
    }
    Ok(residual)
}

pub(super) fn candidate(
    law: &GasReaction, initial_gas: &ClosedGasParcel, duration_s: f64,
    controls: GasReactionSolve, policy: ImplicitReactionSolve,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<GasReactionCandidate, GasReactionError> {
    poll(cancelled)?;
    policy.validate()?;
    if !law.additional_channels.is_empty() {
        return Err(invalid("coupled implicit reaction mechanisms are not admitted"));
    }
    let initial = law.inventory(initial_gas)?;
    let components = law.catalog.thermodynamics().components();
    let equilibrium = law.kinetics.reverse_rate_closure().equilibrium();
    let mut coefficients = vec![0.0; initial.len()];
    let [mut lower, mut upper] = policy.extent_bracket_mol;
    for term in equilibrium.terms() {
        poll(cancelled)?;
        let index = components.binary_search_by(|s| s.species().cmp(term.species()))
            .map_err(|_| invalid("bound reaction species disappeared"))?;
        let coefficient = term.coefficient() as f64;
        coefficients[index] = coefficient;
        // Stay a few rounding units inside the physical donor allowance.
        // This only restricts the nonlinear search; it never alters species.
        let amount = initial[index] * (controls.maximum_consumed_fraction * (1.0 - 8.0 * f64::EPSILON));
        let bound = amount / coefficient.abs();
        if !bound.is_finite() {
            return Err(invalid("implicit donor-extent bound is unrepresentable"));
        }
        if coefficient < 0.0 { upper = upper.min(bound); }
        else { lower = lower.max(-bound); }
    }
    let bracket = [lower, upper];
    let trial = feedback::solve_scalar(bracket, policy.extent_tolerance_mol,
        policy.maximum_iterations, |extent, trials, cancelled| {
        poll(cancelled)?;
        let mut changes = Vec::with_capacity(initial.len());
        let mut expected = Vec::with_capacity(initial.len());
        let mut consumed = 0.0_f64;
        for (i, (&amount, &coefficient)) in initial.iter().zip(&coefficients).enumerate() {
            poll(cancelled)?;
            let change = coefficient * extent;
            let next = amount + change;
            if !change.is_finite() || !next.is_finite() || next < 0.0
                || (change != 0.0 && next == amount)
            { return Err(invalid("implicit extent lost a species increment or made an invalid inventory")); }
            if change < 0.0 {
                let fraction = -change / amount;
                if !fraction.is_finite() || fraction > controls.maximum_consumed_fraction {
                    return Err(GasReactionError::ConsumptionLimit {
                        species: components[i].species().clone(), fraction,
                        maximum: controls.maximum_consumed_fraction,
                    });
                }
                consumed = consumed.max(fraction);
            }
            changes.push(change);
            expected.push(next);
        }
        let gas = law.reconstruct_inventory(initial_gas, &changes, &expected, controls, cancelled)?;
        let actual = law.inventory(&gas)?;
        let value = law.evaluate_channel_trial(&law.kinetics,
            law.rate_material.as_ref(), &gas, &actual, cancelled)?
            .ok_or_else(||
                invalid("implicit reaction rate is log-only; no finite endpoint equation"))?
            .value();
        let residual = equation_residual(extent, value, gas.volume_m3(), duration_s)?;
        let report = ImplicitReactionReport {
            duration_s, admitted_bracket_mol: bracket, residual_mol: residual, trials,
        };
        Ok((Trial {
            gas, actual, changes, expected, consumed, extent_mol: extent, report,
        }, residual))
    }, cancelled).map_err(|error| match error {
        ScalarError::Evaluation(error) => error,
        ScalarError::Invalid(message) => invalid(message),
        ScalarError::Cancelled => GasReactionError::Cancelled,
        ScalarError::OutsideBracket { lower, upper } => GasReactionError::ImplicitOutsideBracket {
            lower_residual_mol: lower, upper_residual_mol: upper,
        },
        ScalarError::NotConverged { iterations, residual } => GasReactionError::ImplicitNotConverged {
            iterations, residual_mol: residual,
        },
    })?;
    // Build complete source evidence only at the accepted endpoint, using the
    // same inventory and live source evaluation as the temporary trial.
    let (rate, rate_material) = law.evaluate_channel(&law.kinetics,
        law.rate_material.as_ref(), &trial.gas, &trial.actual, cancelled)?;
    let value = rate.net_progress_rate().ok_or_else(||
        invalid("implicit reaction rate is log-only; no finite endpoint equation"))?.value();
    let residual = equation_residual(
        trial.extent_mol, value, trial.gas.volume_m3(), duration_s,
    )?;
    if residual.abs() > policy.extent_tolerance_mol {
        return Err(GasReactionError::ImplicitNotConverged {
            iterations: trial.report.trials.saturating_sub(2), residual_mol: residual,
        });
    }
    let progress = mechanism::Progress {
        primary: mechanism::ReactionChannelStep {
            reaction_id: equilibrium.reaction().clone(), rate, rate_material,
            extent_mol: trial.extent_mol,
        },
        additional: Vec::new(), changes: trial.changes,
        expected: trial.expected, consumed: trial.consumed,
    };
    // Reuse the exact explicit-path admission, once for the accepted nonlinear
    // endpoint. Rejected/root-search trials never consume accepted inventory.
    let mut result = law.admit_progress(initial_gas, &initial, trial.gas,
        progress, controls, cancelled)?;
    result.step.implicit = Some(ImplicitReactionReport {
        residual_mol: residual, ..trial.report
    });
    poll(cancelled)?;
    Ok(result)
}

