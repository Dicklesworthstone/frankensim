//! One physical exponential-midpoint schedule for finite inventories and pumps.
//!
//! Numerical propagation stays in the original positive reaction and optical
//! owners. This loop only freezes the authored pair law, retains corrector
//! histories, and controls their joint physical-time refinement.

use super::*;

pub(crate) trait PairHistory {
    fn populations(&self) -> Vec<[f64; 3]>;

    /// Integrated physical channels, and per-stage allowances that add over
    /// stages as those channels do, only: never endpoint stores or residuals.
    fn visit_integrated(
        &mut self,
        visitor: &mut impl FnMut(&mut f64) -> Result<(), SingletTripletError>,
    ) -> Result<(), SingletTripletError>;
}

pub(crate) struct Evolution<S> {
    pub step: S,
    pub substeps: usize,
    pub discrepancy_j: f64,
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn midpoint<S, C, F, A, E>(
    initial: &[[f64; 3]],
    laws: &[Option<TripletAnnihilation>],
    duration: f64,
    fastest_original: f64,
    limits: AnnihilationLimits,
    optical_tolerance: Option<(f64, f64)>,
    mut frozen_step: F,
    mut finalize: A,
    mut components: E,
    cancelled: &mut C,
) -> Result<Evolution<S>, SingletTripletError>
where
    S: PairHistory,
    C: FnMut() -> bool,
    F: FnMut(&[[f64; 3]], &[f64], f64, &mut C) -> Result<Option<S>, SingletTripletError>,
    A: FnMut(&mut S, &mut C) -> Result<(), SingletTripletError>,
    E: FnMut(&S) -> Result<Vec<f64>, SingletTripletError>,
{
    cap("finite annihilation substeps", 2, limits.maximum_substeps)?;
    let mut base_substeps = 1usize;
    let mut previous: Option<(usize, Vec<f64>)> = None;
    loop {
        poll(cancelled)?;
        cap(
            "finite annihilation substeps",
            base_substeps,
            limits.maximum_substeps,
        )?;
        let durations = segment_durations(
            duration,
            fastest_original,
            base_substeps,
            limits.maximum_substeps,
            cancelled,
        )?;
        let actual_substeps = durations.len();
        let mut state = initial.to_vec();
        let mut integrated: Vec<Sum> = Vec::new();
        let mut last = None;
        let mut valid = true;
        for h in durations {
            poll(cancelled)?;
            let half = h * 0.5;
            if half == 0.0 || half * 2.0 != h {
                return Err(SingletTripletError::Unrepresentable(
                    "finite pair predictor duration",
                ));
            }
            let initial_rates = frozen_rates(&state, laws)?;
            if initial_rates.iter().any(|rate| rate * h > 1.0) {
                valid = false;
                break;
            }
            let Some(predictor) = frozen_step(&state, &initial_rates, half, cancelled)? else {
                valid = false;
                break;
            };
            let midpoint_rates = frozen_rates(&predictor.populations(), laws)?;
            if midpoint_rates.iter().any(|rate| rate * h > 1.0) {
                valid = false;
                break;
            }
            let Some(mut corrector) = frozen_step(&state, &midpoint_rates, h, cancelled)? else {
                valid = false;
                break;
            };
            let mut index = 0;
            corrector.visit_integrated(&mut |value| {
                if index == integrated.len() {
                    integrated.push(Sum::default());
                }
                integrated[index].add(*value)?;
                index += 1;
                Ok(())
            })?;
            if index != integrated.len() {
                return Err(SingletTripletError::Invalid(
                    "changing finite pair history topology",
                ));
            }
            state = corrector.populations();
            last = Some(corrector);
        }
        if valid {
            let mut step = last.ok_or(SingletTripletError::Invalid("empty finite pair grid"))?;
            let mut index = 0;
            step.visit_integrated(&mut |value| {
                *value = integrated[index].value()?;
                index += 1;
                Ok(())
            })?;
            finalize(&mut step, cancelled)?;
            let current = components(&step)?;
            if let Some((old_count, old)) = &previous
                && actual_substeps > *old_count
                && !(base_substeps == 2 && *old_count > 1)
            {
                let mut converged = true;
                let mut discrepancy = 0.0_f64;
                for (&before, &after) in old.iter().zip(&current) {
                    let difference = (after - before).abs();
                    let scale = before.abs().max(after.abs());
                    discrepancy = discrepancy.max(difference);
                    converged &= difference
                        <= limits.absolute_energy_tolerance_j + limits.relative_tolerance * scale;
                    if let Some((relative, absolute)) = optical_tolerance {
                        converged &= difference <= absolute + relative * scale;
                    }
                }
                if converged {
                    poll(cancelled)?;
                    return Ok(Evolution {
                        step,
                        substeps: actual_substeps,
                        discrepancy_j: discrepancy,
                    });
                }
            }
            previous = Some((actual_substeps, current));
        } else {
            // A discarded exposure or photon-coupling trial is not a coarse
            // physical history. Its work remains charged by the original owner.
            previous = None;
        }
        base_substeps =
            base_substeps
                .checked_mul(2)
                .ok_or(SingletTripletError::Unrepresentable(
                    "finite pair grid counter",
                ))?;
    }
}

fn segment_durations(
    total: f64,
    fastest: f64,
    base: usize,
    maximum: usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<f64>, SingletTripletError> {
    let duration = total / base as f64;
    if duration == 0.0 || duration * base as f64 != total {
        return Err(SingletTripletError::Unrepresentable(
            "finite pair split duration",
        ));
    }
    let mut tail = Vec::new();
    let mut remaining = duration;
    // Match the existing physical endpoint-tail pattern: refinement also
    // halves the smallest true-time interval, rather than reusing a fixed tail.
    while (remaining * base as f64) * fastest > 1.0 {
        poll(cancelled)?;
        let half = remaining * 0.5;
        if half == 0.0 || half * 2.0 != remaining {
            return Err(SingletTripletError::Unrepresentable(
                "finite pair endpoint duration",
            ));
        }
        remaining = half;
        tail.push(remaining);
        let count = base
            .checked_add(tail.len())
            .ok_or(SingletTripletError::Unrepresentable(
                "finite pair endpoint work",
            ))?;
        cap(
            "finite annihilation substeps including endpoint tail",
            count,
            maximum,
        )?;
    }
    let mut durations = vec![duration; base - 1];
    durations.extend(tail);
    durations.push(remaining);
    Ok(durations)
}
