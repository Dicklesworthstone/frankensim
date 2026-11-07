//! Full-coherency angular refinement over the original Mie coefficient owner.
//! The error is sampled componentwise against the direct series, scaled by the
//! local unpolarized density or 1/(4 pi). It is not a continuous error bound.

use super::{
    Cx, MAX_PHASE_KNOTS, MieAngularSolution, MieError, MiePhaseTable, SphericalPhaseTable,
    SphericalPhaseValue, det, poll,
};

pub(super) fn prepare<S: MieAngularSolution + ?Sized>(
    solution: &S,
    cx: &Cx<'_>,
    tolerance: f64,
    max_knots: usize,
) -> Result<MiePhaseTable, MieError> {
    poll(cx)?;
    if !tolerance.is_finite() || tolerance <= 0.0 || tolerance > 0.1 {
        return Err(MieError::Parameters);
    }
    if solution.cross_sections().scattering_m2 == 0.0 {
        return Err(MieError::NoScattering);
    }
    let intervals = (4 * solution.multipoles()).max(16);
    if max_knots > MAX_PHASE_KNOTS || max_knots <= intervals {
        return Err(MieError::AngularResolution);
    }
    let mut knots = Vec::new();
    knots
        .try_reserve_exact(max_knots)
        .map_err(|_| MieError::Budget)?;
    for index in 0..=intervals {
        poll(cx)?;
        let cosine = if index == 0 {
            -1.0
        } else if index == intervals {
            1.0
        } else {
            -det::cos(core::f64::consts::PI * index as f64 / intervals as f64)
        };
        knots.push((cosine, solution.polarized_value(cosine)?));
    }
    loop {
        let mut next = Vec::new();
        next.try_reserve_exact(max_knots)
            .map_err(|_| MieError::Budget)?;
        next.push(knots[0]);
        let (mut maximum_error, mut normalization) = (0.0_f64, 0.0);
        for pair in knots.windows(2) {
            poll(cx)?;
            let [(lo, left), (hi, right)] = [pair[0], pair[1]];
            let mut error = 0.0_f64;
            for fraction in [0.25, 0.5, 0.75] {
                let actual = solution.polarized_value((hi - lo).mul_add(fraction, lo))?;
                let interpolated = left
                    .interpolate(right, fraction)
                    .map_err(|_| MieError::Numerical)?;
                error = error.max(discrepancy(actual, interpolated));
            }
            maximum_error = maximum_error.max(error);
            normalization += core::f64::consts::PI
                * (hi - lo)
                * (left.unpolarized_per_sr() + right.unpolarized_per_sr());
            if error > tolerance * 0.5 {
                if next.len() + 1 >= max_knots {
                    return Err(MieError::AngularResolution);
                }
                let middle = (hi - lo).mul_add(0.5, lo);
                next.push((middle, solution.polarized_value(middle)?));
            }
            if next.len() == max_knots {
                return Err(MieError::AngularResolution);
            }
            next.push((hi, right));
        }
        if next.len() != knots.len() {
            knots = next;
            continue;
        }
        let polarization = SphericalPhaseTable::try_new(&knots).map_err(|_| MieError::Numerical)?;
        let phase = polarization.scalar_phase().clone();
        let normalization_error = (normalization - 1.0).abs();
        let mean_cosine_error = (phase.mean_cosine() - solution.cross_sections().mean_cosine).abs();
        // Each retained power is at most twice M00. The normalization term
        // therefore carries the same factor two in the componentwise bound.
        let normalized_error = maximum_error + 2.0 * normalization_error / normalization;
        if normalized_error <= tolerance && mean_cosine_error <= tolerance {
            poll(cx)?;
            return Ok(MiePhaseTable {
                phase,
                polarization: Some(polarization),
                normalization_error,
                mean_cosine_error,
                max_sampled_relative_error: normalized_error,
            });
        }
        if 2 * knots.len() - 1 > max_knots {
            return Err(MieError::AngularResolution);
        }
        let mut refined = Vec::new();
        refined
            .try_reserve_exact(max_knots)
            .map_err(|_| MieError::Budget)?;
        refined.push(knots[0]);
        for pair in knots.windows(2) {
            poll(cx)?;
            let middle = (pair[1].0 - pair[0].0).mul_add(0.5, pair[0].0);
            refined.push((middle, solution.polarized_value(middle)?));
            refined.push(pair[1]);
        }
        knots = refined;
    }
}

pub(super) fn discrepancy(a: SphericalPhaseValue, b: SphericalPhaseValue) -> f64 {
    let scale = a
        .unpolarized_per_sr()
        .max(b.unpolarized_per_sr())
        .max(1.0 / (4.0 * core::f64::consts::PI));
    let x = a.powers_per_sr();
    let y = b.powers_per_sr();
    let z = a.cross_per_sr() - b.cross_per_sr();
    [
        (x[0] - y[0]).abs(),
        (x[1] - y[1]).abs(),
        z.re.abs(),
        z.im.abs(),
    ]
    .into_iter()
    .fold(0.0_f64, f64::max)
        / scale
}

#[cfg(test)]
mod tests;
