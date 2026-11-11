//! Bounded adaptive demixing with independent composition and interface accuracy.
//!
//! Full and two-half explicit candidates use the same original Maxwell–Stefan
//! kernel. Initial rates are reused across coarse/first-half trials and retries.
//! Only the accepted fine trajectory contributes physical transfer receipts.
//! The complete returned interval rechecks all original inventories and total
//! bulk-plus-interface free energy; an instantaneous dissipation rate is not a
//! finite-time stability certificate.

use super::*;
use crate::miscibility::adaptive::trial_endpoint;
use super::super::adaptive::{
    MulticomponentMiscibilityTimeControl, admit_time, count_evaluations, count_step,
    error_ratio as composition_error_ratio, factor, retryable,
};
use super::super::solve::Rates;

pub use super::super::adaptive::MulticomponentMiscibilityAdvanceError
    as MulticomponentPhaseFieldAdvanceError;

/// Independent local composition and whole-graph interface-energy allowances.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MulticomponentPhaseFieldTimeControl {
    pub composition: MulticomponentMiscibilityTimeControl,
    /// Allowance on the sum of absolute coarse/fine per-face energy differences.
    pub absolute_interfacial_energy_tolerance_j: f64,
    /// Scale is the sum of max(coarse, fine) energies of each physical face.
    pub relative_interfacial_energy_tolerance: f64,
}

/// Accepted physical path, actual completed work and complete-window balances.
#[derive(Clone, Debug, PartialEq)]
pub struct MulticomponentPhaseFieldAdvance {
    pub duration_s: f64,
    pub attempted_intervals: usize,
    pub accepted_intervals: usize,
    pub rejected_intervals: usize,
    pub completed_steps: usize,
    /// Shared initial rates count once, including across rejected retries.
    /// An incomplete refused evaluation may do additional bounded work.
    pub completed_maxwell_stefan_evaluations: usize,
    pub maximum_accepted_error_ratio: f64,
    pub minimum_accepted_step_s: f64,
    pub maximum_accepted_step_s: f64,
    pub suggested_next_step_s: f64,
    /// Canonical face then original species. Coarse/rejected transfers excluded.
    pub transferred_mol: Vec<Vec<f64>>,
    pub component_residual_mol: Vec<f64>,
    pub maximum_capacity_residual_mol: f64,
    pub bulk_free_energy_change_j: f64,
    pub interfacial_energy_change_j: f64,
    pub free_energy_change_j: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct MulticomponentPhaseFieldAdvanceCandidate {
    field: MulticomponentPhaseField,
    advance: MulticomponentPhaseFieldAdvance,
}
impl MulticomponentPhaseFieldAdvanceCandidate {
    pub const fn field(&self) -> &MulticomponentPhaseField {
        &self.field
    }
    pub const fn advance(&self) -> &MulticomponentPhaseFieldAdvance {
        &self.advance
    }
    pub fn into_parts(
        self,
    ) -> (MulticomponentPhaseField, MulticomponentPhaseFieldAdvance) {
        (self.field, self.advance)
    }
}

impl MulticomponentPhaseField {
    /// Advance one closed demixing interval atomically. Refine positivity,
    /// total-free-energy and temporal-error refusals within the offered limits.
    /// No clipping, equilibrium replacement or hidden interface smoothing occurs.
    /// These are local time-error estimates, not a spatial/global error bound or
    /// an unconditional stability claim for a stiff fourth-order equation.
    pub fn advance_candidate(
        &self,
        duration_s: f64,
        time: MulticomponentPhaseFieldTimeControl,
        controls: MulticomponentMiscibilitySolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> core::result::Result<
        MulticomponentPhaseFieldAdvanceCandidate,
        MulticomponentPhaseFieldAdvanceError,
    > {
        poll(&mut cancelled)?;
        admit_controls(controls)?;
        admit_time(duration_s, time.composition)?;
        if !nonnegative(time.absolute_interfacial_energy_tolerance_j)
            || !nonnegative(time.relative_interfacial_energy_tolerance)
            || !positive(
                time.absolute_interfacial_energy_tolerance_j
                    + time.relative_interfacial_energy_tolerance,
            )
        {
            return Err(MulticomponentMiscibilityError::Invalid(
                "finite nonnegative and nonzero interface time-error tolerances required",
            )
            .into());
        }
        let target = finite(self.elapsed_s() + duration_s, "phase-field target time")?;
        if target <= self.elapsed_s() {
            return Err(MulticomponentMiscibilityError::Invalid(
                "requested phase-field interval does not advance time",
            )
            .into());
        }
        let t = time.composition;
        let n = self.network().species_count();
        let faces = self.network().faces().len();
        let mut field = self.clone();
        let mut transferred = vec![Sum::default(); faces * n];
        let mut dt = t.initial_step_s.min(duration_s);
        let mut attempts = 0usize;
        let mut accepted = 0usize;
        let mut rejected = 0usize;
        let mut completed_steps = 0usize;
        let mut completed_evaluations = 0usize;
        let mut cached_rates: Option<Rates> = None;
        let mut maximum_ratio = 0.0_f64;
        let mut minimum_step = f64::INFINITY;
        let mut maximum_step = 0.0_f64;
        while field.elapsed_s() < target {
            poll(&mut cancelled)?;
            if attempts == t.maximum_attempts {
                return Err(MulticomponentPhaseFieldAdvanceError::WorkLimit {
                    attempted_intervals: attempts,
                    accepted_intervals: accepted,
                    elapsed_s: field.elapsed_s() - self.elapsed_s(),
                });
            }
            let remaining = target - field.elapsed_s();
            let endpoint = trial_endpoint(field.elapsed_s(), dt, target);
            dt = endpoint - field.elapsed_s();
            let floor = t.minimum_step_s.min(remaining);
            let half = 0.5 * dt;
            if !positive(half) || field.elapsed_s() + half <= field.elapsed_s() {
                return Err(MulticomponentMiscibilityError::Invalid(
                    "adaptive phase-field interval cannot be represented",
                )
                .into());
            }
            attempts += 1;
            let trial = (|| -> Result<_> {
                if cached_rates.is_none() {
                    let rates = solve::rates(&field, controls, &mut cancelled)?;
                    count_evaluations(&rates, &mut completed_evaluations)?;
                    cached_rates = Some(rates);
                }
                let rates = cached_rates
                    .as_ref()
                    .expect("rates prepared for unchanged phase field");
                let coarse =
                    solve::candidate_from_rates(&field, dt, rates, controls, &mut cancelled)?;
                count_step(&mut completed_steps)?;
                let first =
                    solve::candidate_from_rates(&field, half, rates, controls, &mut cancelled)?;
                count_step(&mut completed_steps)?;
                let midpoint_rates = solve::rates(first.field(), controls, &mut cancelled)?;
                count_evaluations(&midpoint_rates, &mut completed_evaluations)?;
                let second = solve::candidate_from_rates(
                    first.field(),
                    endpoint - first.field().elapsed_s(),
                    &midpoint_rates,
                    controls,
                    &mut cancelled,
                )?;
                count_step(&mut completed_steps)?;
                if coarse.field().elapsed_s() != endpoint || second.field().elapsed_s() != endpoint {
                    return Err(MulticomponentMiscibilityError::Invalid(
                        "adaptive phase-field endpoint is unrepresentable",
                    ));
                }
                let ratio = composition_error_ratio(
                    &coarse.field().bulk,
                    &second.field().bulk,
                    t,
                    &mut cancelled,
                )?
                .max(interface_error_ratio(
                    coarse.field(),
                    second.field(),
                    time,
                    &mut cancelled,
                )?);
                Ok((first, second, ratio))
            })();
            let (first, second, ratio) = match trial {
                Ok(result) => result,
                Err(error) if retryable(&error) && dt > floor => {
                    rejected += 1;
                    dt = (0.5 * dt).max(floor);
                    continue;
                }
                Err(error) => return Err(error.into()),
            };
            if ratio > 1.0 {
                rejected += 1;
                if dt <= floor {
                    return Err(MulticomponentPhaseFieldAdvanceError::TemporalAccuracy {
                        elapsed_s: field.elapsed_s() - self.elapsed_s(),
                        duration_s: dt,
                        error_ratio: ratio,
                    });
                }
                dt = (dt * factor(ratio)).max(floor);
                continue;
            }
            for report in [first.step(), second.step()] {
                for (face, amounts) in report.transferred_mol.iter().enumerate() {
                    poll(&mut cancelled)?;
                    for (species, &amount) in amounts.iter().enumerate() {
                        transferred[face * n + species].add(amount)?;
                    }
                }
            }
            field = second.into_parts().0;
            cached_rates = None;
            accepted += 1;
            maximum_ratio = maximum_ratio.max(ratio);
            minimum_step = minimum_step.min(dt);
            maximum_step = maximum_step.max(dt);
            dt = (dt * factor(ratio)).clamp(t.minimum_step_s, t.maximum_step_s);
        }
        let (
            component_residual_mol,
            maximum_capacity_residual_mol,
            bulk_free_energy_change_j,
            interfacial_energy_change_j,
            free_energy_change_j,
        ) = self.balances(&field, controls, &mut cancelled)?;
        let mut transferred_mol = Vec::with_capacity(faces);
        for amounts in transferred.chunks(n) {
            poll(&mut cancelled)?;
            transferred_mol.push(
                amounts.iter().map(Sum::value).collect::<Result<Vec<_>>>()?,
            );
        }
        poll(&mut cancelled)?;
        Ok(MulticomponentPhaseFieldAdvanceCandidate {
            field,
            advance: MulticomponentPhaseFieldAdvance {
                duration_s,
                attempted_intervals: attempts,
                accepted_intervals: accepted,
                rejected_intervals: rejected,
                completed_steps,
                completed_maxwell_stefan_evaluations: completed_evaluations,
                maximum_accepted_error_ratio: maximum_ratio,
                minimum_accepted_step_s: minimum_step,
                maximum_accepted_step_s: maximum_step,
                suggested_next_step_s: dt,
                transferred_mol,
                component_residual_mol,
                maximum_capacity_residual_mol,
                bulk_free_energy_change_j,
                interfacial_energy_change_j,
                free_energy_change_j,
            },
        })
    }
}

/// Whole-graph L1 edge-energy error. Comparing total interface energies alone
/// would allow over-resolved and under-resolved interfaces to cancel.
fn interface_error_ratio(
    a: &MulticomponentPhaseField,
    b: &MulticomponentPhaseField,
    time: MulticomponentPhaseFieldTimeControl,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64> {
    let mut discrepancy = Sum::default();
    let mut scale = Sum::default();
    for face in 0..a.network().faces().len() {
        poll(cancelled)?;
        discrepancy.add(solve::edge_energy_change_j(a, b, face)?.abs())?;
        scale.add(solve::edge_energy_j(a, face)?.max(solve::edge_energy_j(b, face)?))?;
    }
    let tolerance = finite(
        time.absolute_interfacial_energy_tolerance_j
            + time.relative_interfacial_energy_tolerance * scale.value()?,
        "phase-field interface temporal tolerance",
    )?;
    let discrepancy = discrepancy.value()?;
    Ok(if discrepancy == 0.0 {
        0.0
    } else if tolerance == 0.0 {
        f64::INFINITY
    } else {
        discrepancy / tolerance
    })
}

#[cfg(test)]
mod tests;
