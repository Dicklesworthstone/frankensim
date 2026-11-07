//! Threshold observations with independent Gaussian sensor noise.
//!
//! For Y = h X + epsilon, let s = h P h' + R, g = P h'/sqrt(s),
//! a = (threshold - h m)/sqrt(s), and lambda = phi(a)/Pr(Z >= a).
//! Gaussian conditioning and total covariance give
//! m_event = m + g lambda, P_event = P_Y + g g' Var(Z | Z >= a).
//! P_Y is the existing Joseph-form covariance conditioned on a point reading.
//! Reversing h and the threshold handles the lower tail. These are exact
//! moment formulas for ONE Gaussian-prior event, not a Gaussian posterior law.

use super::{CensorDirection, ObservationBatch, ObservationState};
use crate::{
    AssimError, Belief, Observation, WorkPlan, WorkProgress, assimilate_checked,
    assimilation_work_plan, checked_dot, clone_belief_for_update, joseph_covariance,
    validate_assimilation_work, validated_canonical_observations,
};
use fs_exec::Cx;

const PHASE: &str = "censored-observation-moments";

pub(super) fn assimilate_independent(
    prior: &Belief,
    batch: &ObservationBatch,
    cx: &Cx<'_>,
) -> Result<Belief, AssimError> {
    let active = batch
        .records
        .iter()
        .filter(|record| !matches!(record.state, ObservationState::Missing { .. }))
        .count();
    validate_assimilation_work(prior.dim(), active)?;
    let mut observations = Vec::with_capacity(active);
    let mut directions = Vec::with_capacity(active);
    for record in &batch.records {
        let (observation, threshold, direction) = match &record.state {
            ObservationState::Available(observation) => (observation, observation.value(), None),
            ObservationState::Missing { .. } => continue,
            ObservationState::Censored {
                observation,
                threshold,
                direction,
            } => (observation, *threshold, Some(*direction)),
            ObservationState::Saturated {
                observation,
                lower,
                upper,
            } if lower < upper => {
                let direction = if observation.value() == *lower {
                    CensorDirection::AtOrBelow
                } else {
                    CensorDirection::AtOrAbove
                };
                (observation, observation.value(), Some(direction))
            }
            _ => {
                return Err(AssimError::InvalidPathologyParameter {
                    parameter: "independent threshold event",
                });
            }
        };
        if observation.operator().len() != prior.dim() {
            return Err(AssimError::DimMismatch {
                state: prior.dim(),
                operator: observation.operator().len(),
            });
        }
        // Never sort or update an event by the placeholder, unobserved value.
        observations.push(Observation::new(
            observation.operator().to_vec(),
            threshold,
            observation.noise_var(),
            observation.instrument(),
        )?);
        directions.push(direction);
    }
    let base = assimilation_work_plan(
        prior,
        &observations,
        0,
        true,
        false,
        None,
        cx.mode().name().len(),
    )?;
    let n = prior.dim() as u128;
    let event_count = directions
        .iter()
        .filter(|direction| direction.is_some())
        .count() as u128;
    // Fixed 256-level tail recurrence, prediction/gain work and the final
    // rank-one covariance term augment the ordinary Joseph/PSD work envelope.
    let extra = event_count * (4_096 + 8 * n * n + 32 * n);
    let plan = WorkPlan::checked(
        base.validation_psd,
        base.record_materialization,
        base.canonical_ordering,
        0,
        base.joseph_update + extra,
        base.psd_revalidation,
        0,
    )?;
    let mut progress = WorkProgress::new(cx, plan)?;
    progress.checkpoint(PHASE)?;
    let canonical = validated_canonical_observations(&observations, prior.dim(), &mut progress)?;
    let mut belief = clone_belief_for_update(prior, &mut progress)?;
    // All point observations preserve Gaussianity and therefore precede every
    // non-Gaussian event. Relative order reuses the established core sorter.
    for index in canonical.ordered_indices() {
        if directions[index].is_none() {
            belief = assimilate_checked(&belief, &observations[index], &mut progress)?;
        }
    }
    for index in canonical.ordered_indices() {
        if let Some(direction) = directions[index] {
            belief = assimilate_event(&belief, &observations[index], direction, &mut progress)?;
        }
    }
    progress.checkpoint("censored-observation-finalize")?;
    Ok(belief)
}

fn assimilate_event(
    prior: &Belief,
    observation: &Observation,
    direction: CensorDirection,
    progress: &mut WorkProgress<'_, '_>,
) -> Result<Belief, AssimError> {
    let h = observation.operator();
    let mut ph = Vec::with_capacity(prior.dim());
    for row in prior.covariance() {
        ph.push(checked_dot(
            row,
            h,
            "censored covariance-times-operator",
            PHASE,
            progress,
        )?);
    }
    let variance = checked_dot(h, &ph, "censored predictive variance", PHASE, progress)?
        + observation.noise_var();
    progress.scalar(PHASE, 1)?;
    finite(variance)?;
    if variance <= 0.0 {
        return Err(AssimError::SingularInnovation);
    }
    let scale = variance.sqrt();
    let prediction = checked_dot(h, prior.mean(), "censored prediction", PHASE, progress)?;
    let sign = match direction {
        CensorDirection::AtOrAbove => 1.0,
        CensorDirection::AtOrBelow => -1.0,
    };
    let residual = observation.value() - prediction;
    // Divide before subtracting only if the otherwise finite physical values
    // overflow their difference. No nonfinite standardized event is admitted.
    let standardized = sign
        * if residual.is_finite() {
            residual / scale
        } else {
            observation.value() / scale - prediction / scale
        };
    let (lambda, tail_scale) = normal_tail_moments(finite(standardized)?, progress)?;
    let mut gain = Vec::with_capacity(prior.dim());
    let mut mean = Vec::with_capacity(prior.dim());
    let mut residual_factor = Vec::with_capacity(prior.dim());
    for (&cross, &prior_mean) in ph.iter().zip(prior.mean()) {
        let standardized_cross = finite(cross / scale)?;
        gain.push(finite(cross / variance)?);
        mean.push(finite(
            (sign * standardized_cross).mul_add(lambda, prior_mean),
        )?);
        residual_factor.push(finite(standardized_cross * tail_scale)?);
        progress.scalar(PHASE, 6)?;
    }
    // Law of total covariance: both terms are PSD. In particular this avoids
    // subtracting nearly equal covariances when an event is deep in the tail.
    let mut covariance = joseph_covariance(prior, h, observation.noise_var(), &gain, progress)?;
    for row in 0..prior.dim() {
        for column in row..prior.dim() {
            let entry = finite(
                residual_factor[row].mul_add(residual_factor[column], covariance[row][column]),
            )?;
            covariance[row][column] = entry;
            covariance[column][row] = entry;
            progress.scalar(PHASE, 1)?;
        }
    }
    Belief::from_covariance_preserving_update(mean, covariance, progress)
}

fn finite(value: f64) -> Result<f64, AssimError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(AssimError::NonFiniteComputation { stage: PHASE })
    }
}

/// Mean and standard deviation of a standard normal conditioned on Z >= a.
/// The dimensionless variance may underflow while the physical covariance is
/// representable, so retain its square root until physical scaling.
/// Fixed-depth numerical evaluation, not an interval-certified enclosure.
fn normal_tail_moments(
    a: f64,
    progress: &mut WorkProgress<'_, '_>,
) -> Result<(f64, f64), AssimError> {
    if a >= 1.5 {
        return positive_tail_moments(a, progress);
    }
    let density = (-0.5 * a * a).exp() * 0.398_942_280_401_432_7;
    let probability = if a <= -1.5 {
        // Symmetry avoids subtracting an almost-one CDF in either rare tail.
        1.0 - density / positive_tail_moments(-a, progress)?.0
    } else {
        // Integral_0^a phi(t)dt = phi(a) * a * Sum a^(2k)/(2k+1)!!.
        // |a| < 1.5 keeps this positive-term series rapidly convergent.
        let mut term = a;
        let mut sum = term;
        for k in 1..=32 {
            term *= a * a / f64::from(2 * k + 1);
            sum += term;
            progress.scalar(PHASE, 5)?;
        }
        (-density).mul_add(sum, 0.5)
    };
    let lambda = density / probability;
    let variance = (-lambda).mul_add(lambda - a, 1.0);
    progress.scalar(PHASE, 12)?;
    checked_moments(lambda, variance.sqrt())
}

fn positive_tail_moments(
    a: f64,
    progress: &mut WorkProgress<'_, '_>,
) -> Result<(f64, f64), AssimError> {
    // Laplace's continued fraction, equivalent to NIST DLMF 7.9.1:
    // lambda(a) = a + 1/(a + 2/(a + 3/(a + ...))).
    // https://dlmf.nist.gov/7.9.E1
    // Differentiate the fraction alongside its value. With b_k=k/(a+b_(k+1)),
    // Var=-b_1'=b_1^2*(1+b_2'). Form sqrt(Var)=b_1*sqrt(1+b_2') directly:
    // b_1^2 can underflow before a large physical scale restores the variance.
    // No exp(-a*a/2) or tail probability is formed on this rare-event path.
    let mut tail = 0.0;
    let mut derivative = 0.0;
    for k in (2..=256).rev() {
        let denominator = a + tail;
        let next = f64::from(k) / denominator;
        progress.scalar(PHASE, 8)?;
        derivative = -(next / denominator) * (1.0 + derivative);
        tail = next;
    }
    let excess = 1.0 / (a + tail);
    progress.scalar(PHASE, 8)?;
    checked_moments(a + excess, excess * (1.0 + derivative).sqrt())
}

fn checked_moments(mean: f64, standard_deviation: f64) -> Result<(f64, f64), AssimError> {
    if mean.is_finite()
        && mean >= 0.0
        && standard_deviation.is_finite()
        && (0.0..=1.0).contains(&standard_deviation)
    {
        Ok((mean, standard_deviation))
    } else {
        Err(AssimError::NonFiniteComputation {
            stage: "normal conditional moments",
        })
    }
}
