//! Scalar bath clock with exact exposure-event boundaries. Within each segment,
//! H=H_left+(H_right-H_left)w² removes a virgin pure-diffusion square-root cusp.
use super::*;

pub(super) fn solve(
    evaluator: &mut Evaluator<'_>,
    target: f64,
    tolerance: f64,
    numerics: PopulationNumerics,
    work: &mut Work<'_, impl FnMut() -> bool>,
) -> Result<(f64, f64)> {
    let mut events = Vec::new();
    for (i, p) in evaluator.batch.layout.populations.iter().enumerate() {
        if p.original_amount.value() > 0.0 && evaluator.batch.penetrations[i] < 1.0 {
            let remaining =
                evaluator.batch.layout.curves[i].remaining(evaluator.batch.penetrations[i])?;
            positive(
                remaining,
                "unfinished core has no representable remaining exposure",
            )?;
            events.push(remaining);
        }
    }
    events.sort_by(f64::total_cmp);
    events.dedup();
    let maximum = *events.last().ok_or(PopulationBatchError::Invalid(
        "no active core for positive progress",
    ))?;
    // Only evaluate the finite concentration at full solid completion. Do NOT
    // integrate its remote endpoint before finding an early finite-time state.
    let completion_concentration = evaluator.evaluate(maximum, work)?;
    let initial_concentration = evaluator.evaluate(0.0, work)?;
    let exposure_bound = initial_concentration * target;
    let mut low = 0.0;
    let mut high = if exposure_bound.is_finite() {
        maximum.min(exposure_bound)
    } else {
        maximum
    };
    positive(high, "finite-time exposure bound loses range")?;
    for _ in 0..numerics.max_iterations {
        work.tick()?;
        let middle = low + 0.5 * (high - low);
        if middle == low || middle == high {
            return Err(PopulationBatchError::Arithmetic(
                "time inversion reached exposure resolution",
            ));
        }
        let dose = middle;
        let concentration = evaluator.evaluate(dose, work)?;
        if concentration <= 0.0 {
            high = middle;
            continue;
        }
        let (time, error) = integrate(evaluator, dose, &events, tolerance * 0.125, work)?;
        let residual = (time - target).abs() + error;
        if residual <= tolerance {
            return Ok((dose, residual));
        }
        if time < target {
            if completion_concentration > 0.0 {
                // C(H) is nonincreasing: remaining time is at most remaining
                // exposure divided by its positive completion concentration.
                let upper = (maximum - dose) / completion_concentration;
                if upper.is_finite() && time + error + upper <= target {
                    let (completed, estimate) =
                        integrate(evaluator, maximum, &events, tolerance * 0.125, work)?;
                    let hold = target - completed;
                    if hold < 0.0 {
                        return Err(PopulationBatchError::Arithmetic(
                            "inconsistent completion-time estimate",
                        ));
                    }
                    let extra = completion_concentration * hold;
                    let dose = maximum + extra;
                    if !dose.is_finite() || (hold > 0.0 && (extra == 0.0 || dose == maximum)) {
                        return Err(PopulationBatchError::Arithmetic(
                            "post-completion exposure loses range",
                        ));
                    }
                    let rounding = ((dose - maximum) / completion_concentration - hold).abs();
                    if estimate + rounding > tolerance {
                        return Err(ShrinkingCoreError::Budget.into());
                    }
                    return Ok((dose, estimate + rounding));
                }
            }
            low = middle;
        } else {
            high = middle;
        }
    }
    Err(ShrinkingCoreError::Budget.into())
}

fn integrand(
    evaluator: &mut Evaluator<'_>,
    z: f64,
    doses: [f64; 2],
    work: &mut Work<'_, impl FnMut() -> bool>,
) -> Result<f64> {
    work.tick()?;
    work.integrands += 1;
    // The square-root coordinate removes the pure-diffusion virgin-core cusp.
    // All admitted initial positive-time batches start with positive reagent.
    if z == 0.0 {
        return Ok(0.0);
    }
    let width = doses[1] - doses[0];
    // Preserve the exact physical event dose at each endpoint; sqrt followed
    // by squaring could spuriously put a completed core one ulp before its end.
    let dose = if z == 1.0 {
        doses[1]
    } else {
        width.mul_add(z * z, doses[0])
    };
    let concentration = evaluator.evaluate(dose, work)?;
    if concentration <= 0.0 {
        return Err(PopulationBatchError::Arithmetic(
            "quadrature reached reagent exhaustion",
        ));
    }
    let value = 2.0 * z * (width / concentration);
    Ok(positive(value, "time integrand loses range")?)
}

fn integrate(
    evaluator: &mut Evaluator<'_>,
    end: f64,
    events: &[f64],
    tolerance: f64,
    work: &mut Work<'_, impl FnMut() -> bool>,
) -> Result<(f64, f64)> {
    let mut boundaries = vec![0.0];
    for &event in events {
        if event > 0.0 && event < end {
            boundaries.push(event);
        }
    }
    boundaries.push(end);
    let per_segment = tolerance / (boundaries.len() - 1) as f64;
    positive(per_segment, "quadrature tolerance loses range")?;
    let mut time = Sum::default();
    let mut error = Sum::default();
    for interval in boundaries.windows(2) {
        let doses = [interval[0], interval[1]];
        let values = [
            0.0,
            integrand(evaluator, 0.5, doses, work)?,
            integrand(evaluator, 1.0, doses, work)?,
        ];
        let whole = simpson(1.0, values)?;
        let (value, estimate) = refine(
            evaluator,
            0.0,
            1.0,
            doses,
            values,
            whole,
            per_segment,
            0,
            work,
        )?;
        time.add(value)?;
        error.add(estimate)?;
    }
    Ok((time.value()?, error.value()?))
}

fn simpson(width: f64, values: [f64; 3]) -> Result<f64> {
    // Scale before adding to avoid overflowing a representable weighted mean.
    let average = values[0] / 6.0 + (2.0 / 3.0) * values[1] + values[2] / 6.0;
    let value = width * average;
    if !value.is_finite() || value < 0.0 || (width > 0.0 && average > 0.0 && value == 0.0) {
        Err(PopulationBatchError::Arithmetic(
            "quadrature area loses range",
        ))
    } else {
        Ok(value)
    }
}

#[allow(clippy::too_many_arguments)]
fn refine(
    evaluator: &mut Evaluator<'_>,
    a: f64,
    b: f64,
    doses: [f64; 2],
    values: [f64; 3],
    whole: f64,
    tolerance: f64,
    depth: usize,
    work: &mut Work<'_, impl FnMut() -> bool>,
) -> Result<(f64, f64)> {
    let middle = a + 0.5 * (b - a);
    let left = a + 0.5 * (middle - a);
    let right = middle + 0.5 * (b - middle);
    if depth >= 48 || left == a || right == b {
        return Err(ShrinkingCoreError::Budget.into());
    }
    let fl = integrand(evaluator, left, doses, work)?;
    let fr = integrand(evaluator, right, doses, work)?;
    let lv = [values[0], fl, values[1]];
    let rv = [values[1], fr, values[2]];
    let lhs = simpson(middle - a, lv)?;
    let rhs = simpson(b - middle, rv)?;
    let sum = lhs + rhs;
    let delta = sum - whole;
    let error = delta.abs() / 15.0;
    if !sum.is_finite() || !error.is_finite() {
        return Err(PopulationBatchError::Arithmetic(
            "quadrature estimate loses range",
        ));
    }
    if error <= tolerance {
        let corrected = sum + delta / 15.0;
        nonnegative(corrected, "corrected quadrature area loses range")?;
        return Ok((corrected, error));
    }
    let (l, le) = refine(
        evaluator,
        a,
        middle,
        doses,
        lv,
        lhs,
        tolerance * 0.5,
        depth + 1,
        work,
    )?;
    let (r, re) = refine(
        evaluator,
        middle,
        b,
        doses,
        rv,
        rhs,
        tolerance * 0.5,
        depth + 1,
        work,
    )?;
    let result = (l + r, le + re);
    nonnegative(result.0, "time integral loses range")?;
    nonnegative(result.1, "quadrature error sum loses range")?;
    Ok(result)
}
