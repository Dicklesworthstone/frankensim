//! Time refinement of complete native thermal/optical candidates.
use super::*;
use fs_thermal::enthalpy::{EnergyAudit, StepConfig};

/// All tolerances refer to temporal discretization, separately from the native
/// joule conservation gate. Relative scales use changes from the input state,
/// so a large arbitrary caloric reference cannot hide an inaccurate trajectory.
#[derive(Clone, Copy, Debug)]
pub struct FeedbackRefinement {
    pub initial_max_dt_s: f64,
    pub minimum_max_dt_s: f64,
    /// At least two: a coarse/fine/finer sequence is always evaluated.
    pub maximum_halvings: usize,
    /// Includes discarded coarse/fine work; each run also retains native and
    /// optical-source budgets. Work is never renewed invisibly during retries.
    pub maximum_total_substeps: usize,
    pub enthalpy_absolute_tolerance_j_kg: f64,
    pub temperature_absolute_tolerance_k: f64,
    pub optical_energy_absolute_tolerance_j: f64,
    pub relative_tolerance: f64,
    pub conservation_tolerance_j: f64,
}

/// Measured differences, not a rigorous global error bound. First-order explicit
/// source integration is claimed only on smooth supported laws and trajectories;
/// constitutive corners and event clipping can change the observed order.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FeedbackConvergence {
    pub thermal: EnergyAudit,
    pub halvings: usize,
    pub total_substeps: usize,
    pub total_source_evaluations: usize,
    pub accepted_max_dt_s: f64,
    pub maximum_enthalpy_difference_j_kg: f64,
    pub maximum_temperature_difference_k: f64,
    pub maximum_optical_energy_difference_j: f64,
    pub normalized_difference: f64,
    /// Every compared coordinate differs by at most 64 binary64 epsilons of
    /// its represented magnitude. Optical accounts use interval increments,
    /// never historic cumulative energy. This diagnostic only relaxes the
    /// monotonic-refinement check; the caller's error tolerances still apply.
    pub roundoff_limited: bool,
    /// log2(previous discrepancy / final discrepancy); absent for exact-zero
    /// or roundoff-limited differences. Refinement must decrease a resolved
    /// discrepancy to accept.
    pub observed_order: Option<f64>,
}

struct Difference {
    enthalpy: f64,
    temperature: f64,
    optical: f64,
    normalized: f64,
    roundoff_limited: bool,
}
impl Default for Difference {
    fn default() -> Self {
        Self {
            enthalpy: 0.0,
            temperature: 0.0,
            optical: 0.0,
            normalized: 0.0,
            roundoff_limited: true,
        }
    }
}

#[derive(Default)]
struct Work {
    substeps: usize,
    evaluations: usize,
}

impl<const N: usize> TemperatureOpticalHeating<N> {
    /// Integrate the actual native conduction/enthalpy/event owner repeatedly
    /// with halved maximum steps. Every run starts at the same input state.
    /// Publish only the finest original trajectory after decreasing successive
    /// discrepancies meet the declared tolerances. If both comparisons reach
    /// the documented roundoff floor, tolerance-qualified nonmonotone roundoff
    /// is also admitted without an observed-order claim. No extrapolated state or
    /// independently integrated temperature is ever constructed.
    ///
    /// Any cancellation, material-domain, conservation or budget error refuses
    /// the whole call. Smaller steps are tried only for temporal disagreement.
    /// Boundary/source discontinuities must be native events or call endpoints.
    /// This estimates time error for this prescribed geometry and quadrature;
    /// it cannot assess spatial, spectral, source-data or constitutive errors.
    pub fn advance_refined<'m>(
        &self,
        thermal: &DepositionSimulation<'m>,
        target_s: f64,
        controls: FeedbackRefinement,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(DepositionSimulation<'m>, Self, FeedbackConvergence), OpticalFeedbackError> {
        checkpoint()?;
        self.check_clock(thermal)?;
        if !positive(controls.initial_max_dt_s)
            || !positive(controls.minimum_max_dt_s)
            || controls.minimum_max_dt_s > controls.initial_max_dt_s
            || controls.maximum_halvings < 2
            || controls.maximum_halvings > 64
            || !positive(controls.enthalpy_absolute_tolerance_j_kg)
            || !positive(controls.temperature_absolute_tolerance_k)
            || !positive(controls.optical_energy_absolute_tolerance_j)
            || !controls.relative_tolerance.is_finite()
            || !(0.0..1.0).contains(&controls.relative_tolerance)
            || !nonnegative(controls.conservation_tolerance_j)
        {
            return Err(invalid(
                "positive time/error bounds, at least two halvings and finite relative tolerance in [0,1) required",
            ));
        }
        let mut work = Work::default();
        let mut dt_s = controls.initial_max_dt_s;
        let mut previous = run(
            thermal,
            self,
            target_s,
            dt_s,
            controls,
            &mut work,
            &mut checkpoint,
        )?;
        let mut previous_difference: Option<(f64, bool)> = None;
        for halving in 1..=controls.maximum_halvings {
            let next_dt = dt_s * 0.5;
            if !positive(next_dt) || next_dt < controls.minimum_max_dt_s {
                return Err(OpticalFeedbackError::Convergence {
                    normalized_error: previous_difference.map_or(f64::INFINITY, |p| p.0),
                    timestep_s: dt_s,
                });
            }
            dt_s = next_dt;
            let current = run(
                thermal,
                self,
                target_s,
                dt_s,
                controls,
                &mut work,
                &mut checkpoint,
            )?;
            let difference = compare(
                thermal,
                self,
                &previous.0,
                &previous.1,
                &current.0,
                &current.1,
                controls,
                &mut checkpoint,
            )?;
            if let Some((prior, prior_roundoff)) = previous_difference {
                let decreasing = difference.normalized == 0.0 || difference.normalized < prior;
                let unresolved_order = prior_roundoff && difference.roundoff_limited;
                if difference.normalized <= 1.0 && (decreasing || unresolved_order) {
                    let observed_order = if !prior_roundoff
                        && !difference.roundoff_limited
                        && prior > 0.0
                        && difference.normalized > 0.0
                    {
                        Some(
                            (fs_math::det::ln(prior) - fs_math::det::ln(difference.normalized))
                                / fs_math::det::ln(2.0),
                        )
                    } else {
                        None
                    };
                    let report = FeedbackConvergence {
                        thermal: current.2,
                        halvings: halving,
                        total_substeps: work.substeps,
                        total_source_evaluations: work.evaluations,
                        accepted_max_dt_s: dt_s,
                        maximum_enthalpy_difference_j_kg: difference.enthalpy,
                        maximum_temperature_difference_k: difference.temperature,
                        maximum_optical_energy_difference_j: difference.optical,
                        normalized_difference: difference.normalized,
                        roundoff_limited: difference.roundoff_limited,
                        observed_order,
                    };
                    return Ok((current.0, current.1, report));
                }
            }
            previous_difference = Some((difference.normalized, difference.roundoff_limited));
            previous = current;
        }
        Err(OpticalFeedbackError::Convergence {
            normalized_error: previous_difference.map_or(f64::INFINITY, |p| p.0),
            timestep_s: dt_s,
        })
    }
}

#[allow(clippy::too_many_arguments)]
fn run<'m, const N: usize>(
    thermal: &DepositionSimulation<'m>,
    optics: &TemperatureOpticalHeating<N>,
    target_s: f64,
    dt_s: f64,
    controls: FeedbackRefinement,
    work: &mut Work,
    checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
) -> Result<
    (
        DepositionSimulation<'m>,
        TemperatureOpticalHeating<N>,
        EnergyAudit,
    ),
    OpticalFeedbackError,
> {
    let remaining = controls
        .maximum_total_substeps
        .checked_sub(work.substeps)
        .ok_or(OpticalFeedbackError::Budget)?;
    let result = thermal.advance_observed_limited(
        optics,
        target_s,
        StepConfig {
            max_dt_s: dt_s,
            energy_tolerance_j: controls.conservation_tolerance_j,
        },
        remaining,
        checkpoint,
    )?;
    work.substeps = work
        .substeps
        .checked_add(result.2.substeps)
        .ok_or(OpticalFeedbackError::Budget)?;
    work.evaluations = work
        .evaluations
        .checked_add(result.1.evaluations)
        .ok_or(OpticalFeedbackError::Budget)?;
    Ok(result)
}

#[allow(clippy::too_many_arguments)]
fn compare<const N: usize>(
    initial: &DepositionSimulation<'_>,
    initial_optics: &TemperatureOpticalHeating<N>,
    coarse: &DepositionSimulation<'_>,
    coarse_optics: &TemperatureOpticalHeating<N>,
    fine: &DepositionSimulation<'_>,
    fine_optics: &TemperatureOpticalHeating<N>,
    controls: FeedbackRefinement,
    checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
) -> Result<Difference, OpticalFeedbackError> {
    checkpoint()?;
    if coarse.cells().len() != fine.cells().len() || coarse.time_s() != fine.time_s() {
        return Err(invalid(
            "refinement candidates disagree on native timeline/topology",
        ));
    }
    let mut difference = Difference::default();
    for (&id, f) in fine.cells() {
        checkpoint()?;
        let c = coarse
            .cells()
            .get(&id)
            .ok_or(invalid("refinement cell mismatch"))?;
        let dh = finite((f.specific_h_j_kg - c.specific_h_j_kg).abs())?;
        // New native births use the absolute allowance until an accepted input
        // state exists; their arbitrary formation-energy datum is never a scale.
        let h_scale = initial.cells().get(&id).map_or(0.0, |base| {
            (f.specific_h_j_kg - base.specific_h_j_kg)
                .abs()
                .max((c.specific_h_j_kg - base.specific_h_j_kg).abs())
        });
        let tf = fine.temperature_k(id)?;
        let tc = coarse.temperature_k(id)?;
        let dt = finite((tf - tc).abs())?;
        let t_scale = if initial.cells().contains_key(&id) {
            let t0 = initial.temperature_k(id)?;
            (tf - t0).abs().max((tc - t0).abs())
        } else {
            0.0
        };
        difference.enthalpy = difference.enthalpy.max(dh);
        difference.temperature = difference.temperature.max(dt);
        difference.roundoff_limited &=
            within_roundoff(dh, f.specific_h_j_kg.abs().max(c.specific_h_j_kg.abs()))
                && within_roundoff(dt, tf.abs().max(tc.abs()));
        difference.normalized = difference.normalized.max(normalize(
            dh,
            h_scale,
            controls.enthalpy_absolute_tolerance_j_kg,
            controls.relative_tolerance,
        )?);
        difference.normalized = difference.normalized.max(normalize(
            dt,
            t_scale,
            controls.temperature_absolute_tolerance_k,
            controls.relative_tolerance,
        )?);
    }
    for lane in 0..N {
        checkpoint()?;
        for (a, b, base) in [
            (
                coarse_optics.sums.absorbed[lane],
                fine_optics.sums.absorbed[lane],
                initial_optics.sums.absorbed[lane],
            ),
            (
                coarse_optics.sums.transmitted[lane],
                fine_optics.sums.transmitted[lane],
                initial_optics.sums.transmitted[lane],
            ),
        ] {
            compare_energy(&mut difference, a, b, base, controls)?;
        }
    }
    compare_energy(
        &mut difference,
        coarse_optics.sums.thermal,
        fine_optics.sums.thermal,
        initial_optics.sums.thermal,
        controls,
    )?;
    Ok(difference)
}

fn compare_energy(
    difference: &mut Difference,
    coarse: CheckedSum,
    fine: CheckedSum,
    base: CheckedSum,
    controls: FeedbackRefinement,
) -> Result<(), OpticalFeedbackError> {
    let change = sum_difference(fine, coarse)?.abs();
    let scale = sum_difference(fine, base)?
        .abs()
        .max(sum_difference(coarse, base)?.abs());
    difference.optical = difference.optical.max(change);
    difference.roundoff_limited &= within_roundoff(change, scale);
    difference.normalized = difference.normalized.max(normalize(
        change,
        scale,
        controls.optical_energy_absolute_tolerance_j,
        controls.relative_tolerance,
    )?);
    Ok(())
}
fn sum_difference(a: CheckedSum, b: CheckedSum) -> Result<f64, OpticalFeedbackError> {
    let mut difference = CheckedSum::default();
    difference.add(a.sum - b.sum)?;
    difference.add(a.correction)?;
    difference.add(-b.correction)?;
    Ok(difference.value()?)
}
fn normalize(
    difference: f64,
    scale: f64,
    absolute: f64,
    relative: f64,
) -> Result<f64, OpticalFeedbackError> {
    let allowance = finite(absolute + relative * scale)?;
    Ok(difference / allowance)
}
fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}

fn within_roundoff(difference: f64, scale: f64) -> bool {
    difference <= 64.0 * f64::EPSILON * scale
}
