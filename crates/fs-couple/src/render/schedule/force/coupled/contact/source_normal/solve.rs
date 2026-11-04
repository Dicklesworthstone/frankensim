//! Scalar original-network condensation with the source's analytic work slope.
use super::{
    CancelGate, NormalLawSource, SourceNormalModalConfig, SourceNormalModalError, finite,
    force_tolerance, mean_force, poll,
};

pub(super) struct Solution {
    pub force_n: f64,
    pub approach_m: f64,
    pub iterations: usize,
}
#[derive(Clone, Copy)]
struct Trial {
    approach: f64,
    force: f64,
    residual: f64,
    slope: f64,
}

pub(super) fn positive_slope(
    source: &NormalLawSource,
    from: f64,
    to: f64,
    compliance: f64,
) -> Result<f64, SourceNormalModalError> {
    let value = slope(source, from, to, compliance)?;
    if value <= 0.0 {
        return Err(SourceNormalModalError::Unresolved(
            "source-normal step needs positive connected mechanical slope",
        ));
    }
    Ok(value)
}
fn slope(
    source: &NormalLawSource,
    from: f64,
    to: f64,
    compliance: f64,
) -> Result<f64, SourceNormalModalError> {
    let derivative =
        source
            .mean_force_end_derivative(from, to)?
            .ok_or(SourceNormalModalError::Unresolved(
                "original source path derivative unavailable",
            ))?;
    Ok(finite(1.0 + compliance * derivative)?)
}

pub(super) fn solve(
    source: &NormalLawSource,
    from: f64,
    free_increment: f64,
    compliance: f64,
    config: SourceNormalModalConfig,
    gate: &CancelGate,
) -> Result<Solution, SourceNormalModalError> {
    // DMT's only negative tangent is its original linear attraction tail.
    // A positive lower bound makes the complete finite-range load line
    // monotone, including both original corners. A coarse ambiguous step
    // refuses; the solve never substitutes a monotone attractive surrogate.
    if let NormalLawSource::DmtSphere(law) = source {
        if config.minimum_approach_m < 0.0 && config.maximum_approach_m > -law.interaction_range_m()
        {
            let minimum = law
                .evaluate(0.0)
                .map_err(fs_contact::normal_patch::NormalLawSourceError::Dmt)?
                .one_sided_tangents_n_per_m[0];
            if finite(1.0 + compliance * (0.5 * minimum))? <= 0.0 {
                return Err(SourceNormalModalError::Unresolved(
                    "DMT finite-range step needs positive original load-line slope",
                ));
            }
        }
    }
    let evaluate = |to: f64| -> Result<Trial, SourceNormalModalError> {
        poll(Some(gate))?;
        let force = mean_force(source, from, to)?;
        Ok(Trial {
            approach: to,
            force,
            residual: finite((to - from) - free_increment + compliance * force)?,
            slope: slope(source, from, to, compliance)?,
        })
    };
    let initial = evaluate(from)?;
    if initial.slope <= 0.0 {
        return Err(SourceNormalModalError::Unresolved(
            "initial source-normal mechanical slope needs a smaller step",
        ));
    }
    let mut iterations = 0;
    if converged(initial, from, free_increment, compliance, config)? {
        return finish(initial, iterations, config);
    }
    let (mut lo, mut hi) = if initial.residual < 0.0 {
        let hi = evaluate(config.maximum_approach_m)?;
        if hi.residual < 0.0 {
            return Err(SourceNormalModalError::Unresolved(
                "source-normal root exceeds the admitted upper approach",
            ));
        }
        (initial, hi)
    } else {
        let mut lo = evaluate(config.minimum_approach_m)?;
        if lo.slope <= 0.0 {
            // JKR F'(delta) increases on its original attached branch, so
            // delta + S*Fbar is convex. Locate the positive-slope side of its
            // unique minimum before bracketing the root connected to `from`.
            // This excludes the lower unstable root without inventing a snap.
            if !matches!(source, NormalLawSource::JkrAttached(_)) {
                return Err(SourceNormalModalError::Unresolved(
                    "source-normal lower bracket has unresolved slope",
                ));
            }
            let mut left = lo;
            let mut right = initial;
            loop {
                tick(&mut iterations, config, gate)?;
                let middle = left.approach.midpoint(right.approach);
                if middle == left.approach || middle == right.approach {
                    return Err(SourceNormalModalError::Unresolved(
                        "JKR stable source-normal root cannot be resolved",
                    ));
                }
                let trial = evaluate(middle)?;
                if trial.slope > 0.0 {
                    if trial.residual <= 0.0 {
                        lo = trial;
                        break;
                    }
                    right = trial;
                } else {
                    left = trial;
                }
            }
        }
        if lo.residual > 0.0 {
            return Err(SourceNormalModalError::Unresolved(
                "source-normal root crosses the admitted lower approach or fold",
            ));
        }
        (lo, initial)
    };
    if lo.slope <= 0.0 || hi.slope <= 0.0 {
        return Err(SourceNormalModalError::Unresolved(
            "source-normal bracket is not on the connected stable branch",
        ));
    }
    for endpoint in [lo, hi] {
        if converged(endpoint, from, free_increment, compliance, config)? {
            return finish(endpoint, iterations, config);
        }
    }
    let mut trial = initial;
    loop {
        tick(&mut iterations, config, gate)?;
        let newton = trial.approach - trial.residual / trial.slope;
        let mut to = if newton.is_finite() && newton > lo.approach && newton < hi.approach {
            newton
        } else {
            lo.approach.midpoint(hi.approach)
        };
        if to == trial.approach {
            to = lo.approach.midpoint(hi.approach);
        }
        if to == lo.approach || to == hi.approach {
            return Err(SourceNormalModalError::Unresolved(
                "source-normal scalar root precision",
            ));
        }
        trial = evaluate(to)?;
        if trial.slope <= 0.0 {
            return Err(SourceNormalModalError::Unresolved(
                "source-normal trial leaves the connected stable branch",
            ));
        }
        if converged(trial, from, free_increment, compliance, config)? {
            return finish(trial, iterations, config);
        }
        if trial.residual < 0.0 {
            lo = trial;
        } else {
            hi = trial;
        }
    }
}

fn tick(
    iterations: &mut usize,
    config: SourceNormalModalConfig,
    gate: &CancelGate,
) -> Result<(), SourceNormalModalError> {
    poll(Some(gate))?;
    if *iterations == config.maximum_iterations {
        return Err(SourceNormalModalError::Unresolved(
            "source-normal scalar refinement budget",
        ));
    }
    *iterations += 1;
    Ok(())
}
fn converged(
    trial: Trial,
    from: f64,
    free_increment: f64,
    compliance: f64,
    config: SourceNormalModalConfig,
) -> Result<bool, SourceNormalModalError> {
    let applied = finite((free_increment - (trial.approach - from)) / compliance)?;
    let force = force_tolerance(applied, trial.force, config)?;
    let coordinate_scale = from
        .abs()
        .max(trial.approach.abs())
        .max(free_increment.abs());
    let coordinate = finite(
        config.coordinate_absolute_tolerance_m
            + config.coordinate_relative_tolerance * coordinate_scale,
    )?;
    // Resolve the scalar root below the outer native compatibility allowance.
    // Its cancellation scale includes the physical contact-force displacement.
    let scale = coordinate_scale.max((compliance * trial.force).abs());
    let roundoff = 8.0 * f64::EPSILON * scale;
    let tolerance = (compliance * force).min(coordinate).min(roundoff);
    Ok(trial.residual.abs() <= tolerance)
}
fn finish(
    trial: Trial,
    iterations: usize,
    config: SourceNormalModalConfig,
) -> Result<Solution, SourceNormalModalError> {
    super::limit(
        "sourced normal held force",
        trial.force.abs(),
        config.maximum_force_n,
    )?;
    Ok(Solution {
        force_n: trial.force,
        approach_m: trial.approach,
        iterations,
    })
}
