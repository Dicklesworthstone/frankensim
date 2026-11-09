//! Retained low parts for stiff three-state cyclic reaction histories.
//!
//! This is the same positive shifted series and semigroup as the general
//! owner, evaluated with the existing double-double primitive. Rounding a
//! nearly stationary base mode to binary64 before many squarings can turn a
//! nanosecond lifetime into spurious population drift on a millisecond clock.
//! Generator construction, scalar damping, counts and source/residence maps
//! retain their low parts until the final conversion and original audits.
//! No column is normalized, no population is clipped, and no fast lifetime
//! becomes a physical integration step. Only a bounded three-node cyclic
//! block whose scaling amplifies binary64 rounding selects this path.

use super::{DecayError, DecayLimits, ReactionPlan, charge, poll, positive_or_zero, source};
use fs_math::dd::Dd;

const N: usize = 3;
const CELLS: usize = N * N;
type Matrix = [Dd; CELLS];

pub(super) fn needs_extended_precision(plan: &ReactionPlan, cyclic: &[bool]) -> bool {
    if plan.nodes != N || !cyclic.iter().any(|&value| value) {
        return false;
    }
    let mut amplified_rounding = f64::EPSILON;
    for _ in 0..plan.squarings {
        amplified_rounding *= 2.0;
    }
    amplified_rounding > plan.tolerance * 0.125
}

#[derive(Clone, Copy)]
struct Blocks {
    p: Matrix,
    c: Matrix,
    f: Matrix,
    g: Matrix,
    h: Matrix,
}

impl Blocks {
    fn initial() -> Self {
        let mut result = Self {
            p: [Dd::ZERO; CELLS],
            c: [Dd::ZERO; CELLS],
            f: [Dd::ZERO; CELLS],
            g: [Dd::ZERO; CELLS],
            h: [Dd::ZERO; CELLS],
        };
        for i in 0..N {
            result.p[i * N + i] = Dd::ONE;
        }
        result
    }
}

pub(super) fn prepare(
    mut plan: ReactionPlan,
    shift: f64,
    time: f64,
    limits: DecayLimits,
    with_source: bool,
    with_residence: bool,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(ReactionPlan, Option<source::SourceMaps>), DecayError> {
    charge(&mut plan.work, 64 * (2 * CELLS + N), limits.maximum_work)?;
    poll(cancelled)?;
    let time = Dd::from_f64(time);
    let mu = Dd::from_f64(shift) * time;
    let removal: [Dd; N] = core::array::from_fn(|i| Dd::from_f64(plan.rates[i]) * time);
    let mut scaled = [Dd::ZERO; CELLS];
    let mut block_norm = mu.to_f64();
    for j in 0..N {
        if plan.rates[j] > 0.0 && removal[j].hi == 0.0 {
            return Err(DecayError::Unrepresentable("scaled direct removal"));
        }
        let mut norm = removal[j];
        for i in 0..N {
            let at = i * N + j;
            scaled[at] = removal[j] * Dd::from_f64(plan.yields[at]);
            if plan.rates[j] > 0.0 && plan.yields[at] > 0.0 && scaled[at].hi == 0.0 {
                return Err(DecayError::Unrepresentable(
                    "scaled direct reaction transfer",
                ));
            }
            if i == j {
                scaled[at] = scaled[at] + (mu - removal[i]);
            }
            if !nonnegative(scaled[at]) {
                return Err(DecayError::Unrepresentable(
                    "scaled positive reaction operator",
                ));
            }
            norm = norm + scaled[at];
        }
        block_norm = block_norm.max(norm.to_f64());
    }
    if !positive_or_zero(block_norm) || block_norm > 0.500_000_000_000_001 {
        return Err(DecayError::Unrepresentable(
            "scaled positive reaction operator",
        ));
    }
    let mut term = Blocks::initial();
    let mut sum = term;
    let mut exponential_term = Dd::ONE;
    let mut exponential = Dd::ONE;
    let mut norm_term = 1.0;
    let mut base_tolerance = limits.relative_tolerance * 0.125;
    for _ in 0..plan.squarings {
        base_tolerance *= 0.5;
    }
    let source_blocks = usize::from(with_source);
    let residence_blocks = usize::from(with_residence);
    let mut converged = false;
    for k in 1..=limits.maximum_series_terms {
        // Include the arithmetic of each retained low part in the original
        // scalar-work budget. Fixed-size temporaries require no extra nodes.
        charge(
            &mut plan.work,
            64 * (N * N * N + (2 + 4 * source_blocks + 2 * residence_blocks) * CELLS + 4),
            limits.maximum_work,
        )?;
        poll(cancelled)?;
        let inverse = Dd::ONE / Dd::from_f64(k as f64);
        let previous = term;
        if with_source {
            for at in 0..CELLS {
                term.f[at] = (mu * previous.f[at] + previous.p[at]) * inverse;
                term.g[at] = (mu * previous.g[at] + previous.c[at]) * inverse;
                sum.f[at] = sum.f[at] + term.f[at];
                sum.g[at] = sum.g[at] + term.g[at];
                if with_residence {
                    term.h[at] = (mu * previous.h[at] + previous.f[at]) * inverse;
                    sum.h[at] = sum.h[at] + term.h[at];
                }
            }
        }
        term.p = multiply(&scaled, &previous.p);
        for i in 0..N {
            for j in 0..N {
                let at = i * N + j;
                term.p[at] = term.p[at] * inverse;
                term.c[at] = (removal[i] * previous.p[at] + mu * previous.c[at]) * inverse;
                sum.p[at] = sum.p[at] + term.p[at];
                sum.c[at] = sum.c[at] + term.c[at];
            }
        }
        // A binary64 exp(-mu) would reintroduce the very base-mode error
        // retained low parts prevent. The reciprocal of the positive scalar
        // exp(mu) series shares the same bounded tail criterion.
        exponential_term = (exponential_term * mu) * inverse;
        exponential = exponential + exponential_term;
        plan.series_terms = k;
        let previous_norm_term = norm_term;
        norm_term *= block_norm / k as f64;
        let tail = norm_term * (block_norm / (k + 1) as f64) / (1.0 - block_norm / (k + 2) as f64);
        let source_tail = norm_term / (1.0 - block_norm / (k + 1) as f64);
        let residence_tail = 0.5 * previous_norm_term / (1.0 - block_norm / k as f64);
        if k >= N + 1 + source_blocks + residence_blocks
            && covered(&sum.p, tail, base_tolerance)
            && covered(&sum.c, tail, base_tolerance)
            && (!with_source
                || (covered(&sum.f, source_tail, base_tolerance)
                    && covered(&sum.g, source_tail, base_tolerance)))
            && (!with_residence || covered(&sum.h, residence_tail, base_tolerance))
        {
            converged = true;
            break;
        }
    }
    if !converged {
        return Err(DecayError::SeriesNotConverged {
            terms: plan.series_terms,
        });
    }
    let damping = Dd::ONE / exponential;
    for block in [&mut sum.p, &mut sum.c, &mut sum.f, &mut sum.g, &mut sum.h] {
        for entry in block {
            *entry = *entry * damping;
        }
    }
    let half = Dd::from_f64(0.5);
    let quarter = Dd::from_f64(0.25);
    for _ in 0..plan.squarings {
        charge(
            &mut plan.work,
            64 * ((2 + 2 * source_blocks + residence_blocks) * (N * N * N + CELLS) + CELLS),
            limits.maximum_work,
        )?;
        poll(cancelled)?;
        let previous = sum;
        sum.p = multiply(&previous.p, &previous.p);
        sum.c = multiply(&previous.c, &previous.p);
        for at in 0..CELLS {
            sum.c[at] = previous.c[at] + sum.c[at];
        }
        if with_source {
            sum.f = multiply(&previous.p, &previous.f);
            sum.g = multiply(&previous.c, &previous.f);
            for at in 0..CELLS {
                sum.f[at] = half * (previous.f[at] + sum.f[at]);
                sum.g[at] = previous.g[at] + half * sum.g[at];
            }
            if with_residence {
                sum.h = multiply(&previous.p, &previous.h);
                for at in 0..CELLS {
                    sum.h[at] = quarter * (previous.h[at] + previous.f[at] + sum.h[at]);
                }
            }
        }
    }
    charge(
        &mut plan.work,
        64 * (5 * CELLS + 2 * N * N * (N + 3)),
        limits.maximum_work,
    )?;
    plan.transition = rounded(&sum.p)?;
    plan.integrated = rounded(&sum.c)?;
    let maps = if with_source {
        let mut maps = source::SourceMaps::identity(N, with_residence);
        maps.population = rounded(&sum.f)?;
        maps.removals = rounded(&sum.g)?;
        if let Some(residence) = &mut maps.residence {
            residence.population = rounded(&sum.h)?;
        }
        Some(maps)
    } else {
        None
    };
    for column in 0..N {
        poll(cancelled)?;
        let initial: [f64; N] = core::array::from_fn(|i| f64::from(i == column));
        let endpoint: [f64; N] = core::array::from_fn(|i| plan.transition[i * N + column]);
        let removals: [f64; N] = core::array::from_fn(|i| plan.integrated[i * N + column]);
        plan.check_balance(&initial, None, &endpoint, &removals, column)?;
        if let Some(maps) = &maps {
            let endpoint: [f64; N] = core::array::from_fn(|i| maps.population[i * N + column]);
            let removals: [f64; N] = core::array::from_fn(|i| maps.removals[i * N + column]);
            plan.check_balance(&initial, None, &endpoint, &removals, column)?;
        }
    }
    poll(cancelled)?;
    Ok((plan, maps))
}

fn multiply(a: &Matrix, b: &Matrix) -> Matrix {
    let mut result = [Dd::ZERO; CELLS];
    for i in 0..N {
        for j in 0..N {
            for k in 0..N {
                result[i * N + j] = result[i * N + j] + a[i * N + k] * b[k * N + j];
            }
        }
    }
    result
}

fn covered(matrix: &Matrix, tail: f64, tolerance: f64) -> bool {
    matrix
        .iter()
        .all(|entry| entry.hi == 0.0 || tail <= tolerance * entry.hi)
}

fn nonnegative(value: Dd) -> bool {
    value.hi.is_finite() && value.lo.is_finite() && !value.lt(Dd::ZERO)
}

fn rounded(matrix: &Matrix) -> Result<Vec<f64>, DecayError> {
    if !matrix.iter().copied().all(nonnegative) {
        return Err(DecayError::Unrepresentable(
            "nonnegative extended reaction map",
        ));
    }
    Ok(matrix.iter().map(|entry| entry.to_f64()).collect())
}

#[cfg(test)]
#[path = "small/tests.rs"]
mod tests;
