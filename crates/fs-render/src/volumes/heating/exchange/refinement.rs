//! Bounded time refinement of native radiative cooling/heating trajectories.
//!
//! Every trial retains the SAME sampled optical operator. Differences assess
//! temporal integration conditional on that operator, never Monte Carlo,
//! spectral, spatial, or constitutive accuracy. No extrapolated thermal state
//! or independently evolved temperature is introduced.

use super::{CheckedSum, RadiativeExchangeError, VolumeRadiativeExchange};
use fs_math::det;
use fs_thermal::enthalpy::{DepositionError, DepositionSimulation, EnergyAudit, StepConfig};

/// Time/error allowances for a coarse/fine/finer sequence. Discarded trials
/// consume the same total work allowances as the accepted trial.
#[derive(Clone, Copy, Debug)]
pub struct RadiativeExchangeRefinement {
    /// Starting maximum native timestep, before native stability/event limits.
    pub initial_max_dt_s: f64,
    /// Smallest permitted maximum timestep during refinement.
    pub minimum_max_dt_s: f64,
    /// At least two and at most 64 successive halvings.
    pub maximum_halvings: usize,
    /// Total actual native substeps, including all discarded trajectories.
    pub maximum_total_substeps: usize,
    /// Total radiation source preparations over every trajectory.
    pub maximum_total_source_evaluations: usize,
    /// Absolute final specific-enthalpy difference allowance.
    pub enthalpy_absolute_tolerance_j_kg: f64,
    /// Absolute final native-temperature difference allowance.
    pub temperature_absolute_tolerance_k: f64,
    /// Absolute interval bath/internal/unresolved-energy difference allowance.
    pub radiative_energy_absolute_tolerance_j: f64,
    /// Relative allowance uses changes from the original accepted state.
    pub relative_tolerance: f64,
    /// Original native joule-conservation gate, separate from time accuracy.
    pub conservation_tolerance_j: f64,
}

/// Measured time discrepancies on the finest accepted native trajectory.
/// Two successive comparisons must meet the declared tolerances AND change
/// actual native step counts. A smaller requested dt that leaves a CFL-limited
/// trajectory unchanged cannot be evidence of time convergence.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RadiativeExchangeConvergence {
    /// The accepted trajectory's original thermal conservation audit.
    pub thermal: EnergyAudit,
    /// Number of maximum-timestep halvings actually evaluated.
    pub halvings: usize,
    /// Includes discarded trial work.
    pub total_substeps: usize,
    /// Includes discarded trial source preparations.
    pub total_source_evaluations: usize,
    /// Maximum timestep used by the accepted original trajectory.
    pub accepted_max_dt_s: f64,
    /// Actual substeps in the preceding, coarser trajectory.
    pub penultimate_substeps: usize,
    /// Largest final specific-enthalpy difference across native cells.
    pub maximum_enthalpy_difference_j_kg: f64,
    /// Largest final native-temperature difference across cells.
    pub maximum_temperature_difference_k: f64,
    /// Largest interval radiation-account difference, including each lane.
    pub maximum_radiative_energy_difference_j: f64,
    /// Largest discrepancy divided by its requested allowance, at most one.
    pub normalized_difference: f64,
    /// The preceding comparison must also be within its allowance.
    pub previous_normalized_difference: f64,
    /// Native h/T coordinate resolution divided by the requested allowance.
    /// A value above one refuses convergence; no tolerance is widened.
    pub maximum_coordinate_resolution_ratio: f64,
    /// Both final comparisons lie within a documented binary64 roundoff scale.
    pub roundoff_limited: bool,
    /// log2 of successive normalized discrepancies, when resolved and positive.
    /// This is observed behavior, not a global error or smoothness certificate.
    pub observed_order: Option<f64>,
}

/// Refinement refusal leaves the original thermal state and observer untouched.
#[derive(Clone, Debug, PartialEq)]
pub enum RadiativeExchangeRefinementError {
    /// Invalid time or error controls.
    Invalid(&'static str),
    /// The shared refinement work allowance is exhausted.
    Budget,
    /// Allowed refinements did not demonstrate temporal agreement.
    Convergence {
        /// Last measured normalized discrepancy, or infinity before comparison.
        normalized_difference: f64,
        /// Last evaluated maximum timestep.
        timestep_s: f64,
        /// Number of halvings evaluated.
        halvings: usize,
    },
    /// Original material, transport, thermal, arithmetic, or cancellation error.
    Exchange(RadiativeExchangeError),
}
impl From<RadiativeExchangeError> for RadiativeExchangeRefinementError {
    fn from(error: RadiativeExchangeError) -> Self {
        Self::Exchange(error)
    }
}
impl From<DepositionError> for RadiativeExchangeRefinementError {
    fn from(error: DepositionError) -> Self {
        Self::Exchange(error.into())
    }
}
impl From<super::super::OpticalHeatingError> for RadiativeExchangeRefinementError {
    fn from(error: super::super::OpticalHeatingError) -> Self {
        Self::Exchange(error.into())
    }
}
impl core::fmt::Display for RadiativeExchangeRefinementError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "radiative time refinement refused: {self:?}")
    }
}
impl core::error::Error for RadiativeExchangeRefinementError {}

type Result<T> = core::result::Result<T, RadiativeExchangeRefinementError>;

#[derive(Clone, Copy)]
struct Difference {
    enthalpy: f64,
    temperature: f64,
    radiative: f64,
    normalized: f64,
    coordinate_resolution: f64,
    roundoff_limited: bool,
    actual_refinement: bool,
}
impl Default for Difference {
    fn default() -> Self {
        Self {
            enthalpy: 0.0,
            temperature: 0.0,
            radiative: 0.0,
            normalized: 0.0,
            coordinate_resolution: 0.0,
            roundoff_limited: true,
            actual_refinement: false,
        }
    }
}
#[derive(Default)]
struct Work {
    substeps: usize,
    evaluations: usize,
}

impl<const N: usize> VolumeRadiativeExchange<N> {
    /// Refine complete original thermal/observer advances from the same input
    /// state. Accept only two successive tolerance-qualified comparisons with
    /// increasing actual substep counts and decreasing resolved discrepancy.
    /// Roundoff may relax decreasing-order evidence, never the caller's error
    /// allowances or the requirement to refine actual native timesteps.
    ///
    /// Spectral support, material domains, source/error allowances and native
    /// event ordering remain those of the original owners. A domain, budget,
    /// cancellation or conservation failure refuses the entire operation;
    /// retries are reserved for temporal disagreement. Discontinuities must be
    /// native events or call endpoints. The fixed sampled transport operator is
    /// shared in Arc, so this operation neither reruns photon histories nor
    /// interprets their correlated coefficient errors as thermal-state errors.
    pub fn advance_refined<'m>(
        &self,
        thermal: &DepositionSimulation<'m>,
        target_s: f64,
        controls: RadiativeExchangeRefinement,
        mut checkpoint: impl FnMut() -> core::result::Result<(), DepositionError>,
    ) -> Result<(DepositionSimulation<'m>, Self, RadiativeExchangeConvergence)> {
        checkpoint()?;
        self.check_clock(thermal)?;
        if !target_s.is_finite()
            || target_s <= thermal.time_s()
            || !positive(controls.initial_max_dt_s)
            || !positive(controls.minimum_max_dt_s)
            || controls.minimum_max_dt_s > controls.initial_max_dt_s
            || !(2..=64).contains(&controls.maximum_halvings)
            || !positive(controls.enthalpy_absolute_tolerance_j_kg)
            || !positive(controls.temperature_absolute_tolerance_k)
            || !positive(controls.radiative_energy_absolute_tolerance_j)
            || !controls.relative_tolerance.is_finite()
            || !(0.0..1.0).contains(&controls.relative_tolerance)
            || !super::nonnegative(controls.conservation_tolerance_j)
        {
            return Err(RadiativeExchangeRefinementError::Invalid(
                "positive time interval and error bounds, 2..=64 halvings and finite relative tolerance in [0,1) required",
            ));
        }
        let mut work = Work::default();
        let mut dt = controls.initial_max_dt_s;
        let mut previous = run(
            thermal,
            self,
            target_s,
            dt,
            controls,
            &mut work,
            &mut checkpoint,
        )?;
        let mut prior_difference: Option<Difference> = None;
        for halving in 1..=controls.maximum_halvings {
            let next_dt = dt * 0.5;
            if !positive(next_dt) || next_dt < controls.minimum_max_dt_s {
                return Err(convergence(prior_difference, dt, halving - 1));
            }
            dt = next_dt;
            let mut current = run(
                thermal,
                self,
                target_s,
                dt,
                controls,
                &mut work,
                &mut checkpoint,
            )?;
            let mut difference = compare(
                thermal,
                self,
                &previous.0,
                &previous.1,
                &current.0,
                &current.1,
                controls,
                &mut checkpoint,
            )?;
            difference.actual_refinement = current.2.substeps > previous.2.substeps;
            if let Some(prior) = prior_difference {
                let roundoff = prior.roundoff_limited && difference.roundoff_limited;
                let decreasing =
                    difference.normalized == 0.0 || difference.normalized < prior.normalized;
                if prior.actual_refinement
                    && difference.actual_refinement
                    && prior.normalized <= 1.0
                    && difference.normalized <= 1.0
                    && prior.coordinate_resolution <= 1.0
                    && difference.coordinate_resolution <= 1.0
                    && (decreasing || roundoff)
                {
                    let observed_order = if !prior.roundoff_limited
                        && !difference.roundoff_limited
                        && prior.normalized > 0.0
                        && difference.normalized > 0.0
                    {
                        Some(
                            (det::ln(prior.normalized) - det::ln(difference.normalized))
                                / det::ln(2.0),
                        )
                    } else {
                        None
                    };
                    let report = RadiativeExchangeConvergence {
                        thermal: current.2,
                        halvings: halving,
                        total_substeps: work.substeps,
                        total_source_evaluations: work.evaluations,
                        accepted_max_dt_s: dt,
                        penultimate_substeps: previous.2.substeps,
                        maximum_enthalpy_difference_j_kg: difference.enthalpy,
                        maximum_temperature_difference_k: difference.temperature,
                        maximum_radiative_energy_difference_j: difference.radiative,
                        normalized_difference: difference.normalized,
                        previous_normalized_difference: prior.normalized,
                        maximum_coordinate_resolution_ratio: difference.coordinate_resolution,
                        roundoff_limited: roundoff,
                        observed_order,
                    };
                    // The study's remaining allowance does not permanently
                    // tighten the accepted observer's subsequent advances.
                    current.1.evaluation_limit = self.evaluation_limit;
                    checkpoint()?;
                    return Ok((current.0, current.1, report));
                }
            }
            prior_difference = Some(difference);
            previous = current;
        }
        Err(convergence(prior_difference, dt, controls.maximum_halvings))
    }
}

fn convergence(
    difference: Option<Difference>,
    timestep_s: f64,
    halvings: usize,
) -> RadiativeExchangeRefinementError {
    RadiativeExchangeRefinementError::Convergence {
        normalized_difference: difference
            .map_or(f64::INFINITY, |d| d.normalized.max(d.coordinate_resolution)),
        timestep_s,
        halvings,
    }
}

#[allow(clippy::too_many_arguments)]
fn run<'m, const N: usize>(
    thermal: &DepositionSimulation<'m>,
    optics: &VolumeRadiativeExchange<N>,
    target: f64,
    dt: f64,
    controls: RadiativeExchangeRefinement,
    work: &mut Work,
    checkpoint: &mut dyn FnMut() -> core::result::Result<(), DepositionError>,
) -> Result<(
    DepositionSimulation<'m>,
    VolumeRadiativeExchange<N>,
    EnergyAudit,
)> {
    let substeps = controls
        .maximum_total_substeps
        .checked_sub(work.substeps)
        .ok_or(RadiativeExchangeRefinementError::Budget)?;
    let evaluations = controls
        .maximum_total_source_evaluations
        .checked_sub(work.evaluations)
        .ok_or(RadiativeExchangeRefinementError::Budget)?;
    if substeps == 0 || evaluations == 0 {
        return Err(RadiativeExchangeRefinementError::Budget);
    }
    let limited = optics.with_work_limit(evaluations);
    let result = thermal.advance_observed_limited(
        &limited,
        target,
        StepConfig {
            max_dt_s: dt,
            energy_tolerance_j: controls.conservation_tolerance_j,
        },
        substeps,
        checkpoint,
    )?;
    work.substeps = work
        .substeps
        .checked_add(result.2.substeps)
        .ok_or(RadiativeExchangeRefinementError::Budget)?;
    work.evaluations = work
        .evaluations
        .checked_add(result.1.evaluations)
        .ok_or(RadiativeExchangeRefinementError::Budget)?;
    Ok(result)
}

#[allow(clippy::too_many_arguments)]
fn compare<const N: usize>(
    initial: &DepositionSimulation<'_>,
    initial_optics: &VolumeRadiativeExchange<N>,
    coarse: &DepositionSimulation<'_>,
    coarse_optics: &VolumeRadiativeExchange<N>,
    fine: &DepositionSimulation<'_>,
    fine_optics: &VolumeRadiativeExchange<N>,
    controls: RadiativeExchangeRefinement,
    checkpoint: &mut dyn FnMut() -> core::result::Result<(), DepositionError>,
) -> Result<Difference> {
    checkpoint()?;
    if coarse.cells().len() != fine.cells().len() || coarse.time_s() != fine.time_s() {
        return Err(RadiativeExchangeRefinementError::Invalid(
            "native candidate topology/timeline mismatch",
        ));
    }
    let mut result = Difference::default();
    for (&id, f) in fine.cells() {
        checkpoint()?;
        let c = coarse
            .cells()
            .get(&id)
            .ok_or(RadiativeExchangeRefinementError::Invalid(
                "native cell mismatch",
            ))?;
        let dh = finite((f.specific_h_j_kg - c.specific_h_j_kg).abs())?;
        let h_scale = match initial.cells().get(&id) {
            Some(base) => finite(
                (f.specific_h_j_kg - base.specific_h_j_kg)
                    .abs()
                    .max((c.specific_h_j_kg - base.specific_h_j_kg).abs()),
            )?,
            // A new event-born cell has no accepted input datum for a relative
            // scale; its arbitrary formation energy cannot become an allowance.
            None => 0.0,
        };
        let tf = fine.temperature_k(id)?;
        let tc = coarse.temperature_k(id)?;
        let dt = finite((tf - tc).abs())?;
        let t_scale = if initial.cells().contains_key(&id) {
            let t0 = initial.temperature_k(id)?;
            finite((tf - t0).abs().max((tc - t0).abs()))?
        } else {
            0.0
        };
        result.enthalpy = result.enthalpy.max(dh);
        result.temperature = result.temperature.max(dt);
        for (difference, scale, magnitude, absolute) in [
            (
                dh,
                h_scale,
                f.specific_h_j_kg.abs().max(c.specific_h_j_kg.abs()),
                controls.enthalpy_absolute_tolerance_j_kg,
            ),
            (
                dt,
                t_scale,
                tf.abs().max(tc.abs()),
                controls.temperature_absolute_tolerance_k,
            ),
        ] {
            result.normalized = result.normalized.max(normalize(
                difference,
                scale,
                absolute,
                controls.relative_tolerance,
            )?);
            result.roundoff_limited &= within_roundoff(difference, magnitude);
            // Equal rounded coordinates cannot establish accuracy finer than
            // their binary64 resolution. This is a floor, not an error bound.
            let floor = (f64::EPSILON * magnitude).max(f64::from_bits(1));
            result.coordinate_resolution = result.coordinate_resolution.max(normalize(
                floor,
                scale,
                absolute,
                controls.relative_tolerance,
            )?);
        }
    }
    for lane in 0..N {
        checkpoint()?;
        for (c, f, base) in [
            (
                coarse_optics.sums.bath[lane],
                fine_optics.sums.bath[lane],
                initial_optics.sums.bath[lane],
            ),
            (
                coarse_optics.sums.internal[lane],
                fine_optics.sums.internal[lane],
                initial_optics.sums.internal[lane],
            ),
        ] {
            compare_energy(&mut result, c, f, base, controls)?;
        }
    }
    for (c, f, base) in [
        (
            coarse_optics.sums.external,
            fine_optics.sums.external,
            initial_optics.sums.external,
        ),
        (
            coarse_optics.sums.unresolved,
            fine_optics.sums.unresolved,
            initial_optics.sums.unresolved,
        ),
    ] {
        compare_energy(&mut result, c, f, base, controls)?;
    }
    Ok(result)
}

fn compare_energy(
    result: &mut Difference,
    coarse: CheckedSum,
    fine: CheckedSum,
    base: CheckedSum,
    controls: RadiativeExchangeRefinement,
) -> Result<()> {
    let difference = sum_difference(fine, coarse)?.abs();
    let scale = sum_difference(fine, base)?
        .abs()
        .max(sum_difference(coarse, base)?.abs());
    result.radiative = result.radiative.max(difference);
    result.roundoff_limited &= within_roundoff(difference, scale);
    result.normalized = result.normalized.max(normalize(
        difference,
        scale,
        controls.radiative_energy_absolute_tolerance_j,
        controls.relative_tolerance,
    )?);
    Ok(())
}

// Compare compensated interval increments before collapsing the accounts to a
// single f64, so a tiny new heat increment survives a huge historic fluence.
fn sum_difference(a: CheckedSum, b: CheckedSum) -> Result<f64> {
    let mut result = CheckedSum::default();
    result.add(finite(a.sum - b.sum)?)?;
    result.add(a.correction)?;
    result.add(-b.correction)?;
    Ok(result.value()?)
}
fn normalize(difference: f64, scale: f64, absolute: f64, relative: f64) -> Result<f64> {
    let allowance = finite(absolute + relative * scale)?;
    Ok(difference / allowance)
}
fn finite(value: f64) -> Result<f64> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(RadiativeExchangeRefinementError::Invalid(
            "unrepresentable refinement coordinate",
        ))
    }
}
fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn within_roundoff(difference: f64, scale: f64) -> bool {
    difference <= 64.0 * f64::EPSILON * scale
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn interval_accounts_preserve_heat_below_historic_fluence_ulp() {
        let mut base = CheckedSum::default();
        base.add(1e25).unwrap();
        let mut coarse = base;
        let mut fine = base;
        coarse.add(-0.125).unwrap();
        fine.add(-0.25).unwrap();
        assert_eq!(base.value().unwrap(), fine.value().unwrap());
        assert_eq!(sum_difference(coarse, base).unwrap(), -0.125);
        assert_eq!(sum_difference(fine, base).unwrap(), -0.25);
        assert_eq!(sum_difference(fine, coarse).unwrap(), -0.125);
    }
}
