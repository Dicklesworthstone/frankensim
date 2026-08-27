//! Baseline NSGA-III (epic bead `frankensim-epic-ascent-7tv.24.13`): a
//! deterministic gradient-free population optimizer with its complete
//! DEFECT-MATRIX battery frozen BEFORE any adaptive-quality claim rides
//! on it.
//!
//! Design contract obeyed here (see crate CONTRACT.md house rules):
//! determinism first — every ordering decision goes through
//! `f64::total_cmp` plus an explicit index tie-break so the result is a
//! pure function of the input SET; refusals are typed and loud;
//! duplicate objective vectors collapse stably; budget accounting is
//! integral (evaluations never exceed the authored ceiling; a partial
//! generation is discarded rather than half-reported); the returned
//! population always carries its non-dominated FRONT indices so
//! callers can audit what was claimed.
//!
//! AUTHORITY: Estimate-class. This module provides no convergence
//! certificate and no quality claim about optima; those live with the
//! KKT-bearing tracing engines in [`crate::pareto`].

use crate::stop::StopRule;

/// Typed refusals from NSGA-III setup and execution.
#[derive(Debug, Clone, PartialEq)]
pub enum NsgaError {
    /// The initial population is empty.
    PopulationEmpty,
    /// Individuals disagree on objective dimensionality.
    ObjectiveCountMismatch {
        /// Dimension of the first individual.
        expected: usize,
        /// Offending individual position (0-based).
        at: usize,
        /// Its observed dimension.
        found: usize,
    },
    /// A non-finite objective or decision value entered the pipeline.
    NonFinite {
        /// Individual position.
        individual: usize,
        /// Component position inside that vector.
        component: usize,
        /// Which vector tripped: objectives or decisions.
        kind: NonFiniteKind,
    },
    /// The reference-direction request is malformed.
    ReferenceInvalid {
        /// Human-readable refusal reason (authored, receipt-visible).
        what: String,
    },
    /// The authored evaluation budget cannot host even one complete
    /// generation, so nothing could ever converge honestly.
    BudgetInfeasible {
        /// Minimum individuals per generation.
        minimum_population: usize,
        /// Authored budget.
        budget: usize,
    },
}

/// Which vector carried the non-finite sample.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NonFiniteKind {
    /// Objective vector component.
    Objective,
    /// Decision vector component.
    Decision,
}

impl core::fmt::Display for NsgaError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::PopulationEmpty => write!(f, "initial population is empty"),
            Self::ObjectiveCountMismatch {
                expected,
                at,
                found,
            } => write!(
                f,
                "individual {at} has {found} objectives; expected {expected}"
            ),
            Self::NonFinite {
                individual,
                component,
                kind,
            } => match kind {
                NonFiniteKind::Objective => write!(
                    f,
                    "individual {individual} objective {component} is non-finite"
                ),
                NonFiniteKind::Decision => write!(
                    f,
                    "individual {individual} decision {component} is non-finite"
                ),
            },
            Self::ReferenceInvalid { what } => write!(f, "reference invalid: {what}"),
            Self::BudgetInfeasible {
                minimum_population,
                budget,
            } => write!(
                f,
                "budget {budget} cannot host one generation of \
                 {minimum_population} individuals"
            ),
        }
    }
}

impl core::error::Error for NsgaError {}

/// One candidate: an uninterpreted decision vector and its objective
/// vector. `x` exists so tie-breaks and duplicate policy can key on
/// decision identity independently of objective bits.
#[derive(Debug, Clone, PartialEq)]
pub struct NsgaIndividual {
    /// Decision variables (free-form; box-checked only when bounds are
    /// supplied through [`NsgaConfig`]).
    pub x: Vec<f64>,
    /// Minimization objectives.
    pub f: Vec<f64>,
}

impl NsgaIndividual {
    /// Total-order bit key used wherever two components must compare
    /// deterministically without float semantics leaking (-0.0/+0.0,
    /// NaN exclusion happens at admission).
    fn lt_bit(a: &[f64], b: &[f64]) -> bool {
        for (x, y) in a.iter().zip(b.iter()) {
            match x.total_cmp(y) {
                core::cmp::Ordering::Less => return true,
                core::cmp::Ordering::Greater => return false,
                core::cmp::Ordering::Equal => {}
            }
        }
        a.len() < b.len()
    }

    /// Lexicographic total-order comparison by objectives, then
    /// decisions, then length. Canonical across enumeration order.
    #[must_use]
    pub fn canonical_key_ordering(&self, other: &Self) -> core::cmp::Ordering {
        if !core::ptr::eq(self, other) {
            let fo = Self::vec_ordering(&self.f, &other.f);
            if fo != core::cmp::Ordering::Equal {
                return fo;
            }
            let xo = Self::vec_ordering(&self.x, &other.x);
            if xo != core::cmp::Ordering::Equal {
                return xo;
            }
        }
        core::cmp::Ordering::Equal
    }

    fn vec_ordering(a: &[f64], b: &[f64]) -> core::cmp::Ordering {
        let n = a.len().min(b.len());
        for i in 0..n {
            match a[i].total_cmp(&b[i]) {
                core::cmp::Ordering::Equal => {}
                other => return other,
            }
        }
        a.len().cmp(&b.len())
    }

    /// Strict Pareto dominance under minimization: no worse everywhere,
    /// strictly better somewhere.
    #[must_use]
    pub fn dominates(&self, other: &Self) -> bool {
        let mut strictly_better = false;
        for (fi, fj) in self.f.iter().zip(other.f.iter()) {
            if fi > fj {
                return false;
            }
            if fi < fj {
                strictly_better = true;
            }
        }
        strictly_better
    }
}

/// Das–Dennis reference directions on the unit simplex, generated from
/// integer partitions (exact, order-canonical: lex ascending) so the
/// reference SET never depends on floating accumulation order.
///
/// # Errors
/// [`NsgaError::ReferenceInvalid`] when `divisions == 0`, `m < 2`, or
/// the partition family would be empty.
#[must_use]
pub fn das_dennis(divisions: usize, m: usize) -> Vec<Vec<f64>> {
    if divisions == 0 || m < 2 {
        return Vec::new();
    }
    // Enumerate all m-tuples of non-negative ints summing to
    // `divisions`, lex ascending — recursive depth-first.
    let mut dirs_int: Vec<Vec<usize>> = Vec::new();
    let mut prefix = vec![0usize; m];
    fn walk(m: usize, i: usize, left: usize, prefix: &mut [usize], out: &mut Vec<Vec<usize>>) {
        if i == m - 1 {
            prefix[i] = left;
            out.push(prefix.to_vec());
            return;
        }
        for v in 0..=left {
            prefix[i] = v;
            walk(m, i + 1, left - v, prefix, out);
        }
        prefix[i] = 0;
    }
    walk(m, 0, divisions, &mut prefix, &mut dirs_int);
    let denom = divisions as f64;
    dirs_int
        .into_iter()
        .map(|t| t.into_iter().map(|v| v as f64 / denom).collect())
        .collect()
}

/// Cardinality C(p+m−1, m−1): the exhaustive permutation oracle's
/// independent closed form.
#[must_use]
pub fn das_dennis_cardinality(divisions: usize, m: usize) -> usize {
    let n = divisions + m.saturating_sub(1);
    let k = m.saturating_sub(1);
    (0..k).fold(1usize, |acc, i| acc * (n - i) / (i + 1))
}

/// Outcome of the ideal/extreme/intercept normalization stage.
#[derive(Debug, Clone)]
struct Normalization {
    /// Per-objective translation subtracted before scaling.
    ideal: Vec<f64>,
    /// Per-objective scale divisors; never zero or negative.
    scales: Vec<f64>,
    /// Whether the extreme-axis system was singular and theASF-free
    /// nadir fallback engaged (recorded for receipts).
    singular_fallback: bool,
}

/// Compute normalization following NSGA-III §IV.B: translate by ideal
/// point; find axis extremes by minimizing the Achievement Scalarizing
/// Function over the union population; solve the m×m intercept system.
/// Any degenerate/singular detour falls back to the plain nadir scale
/// (max minus ideal), keeping the outcome deterministic either way.
///
/// Returns `None` only when the population itself is empty (caller
/// refuses earlier, so this stays internal).
fn normalize_population(
    pop: &[NsgaIndividual],
    m: usize,
) -> Normalization {
    let mut ideal = vec![f64::INFINITY; m];
    for ind in pop {
        for j in 0..m {
            if ind.f[j] < ideal[j] {
                ideal[j] = ind.f[j];
            }
        }
    }
    // Translated objectives T = f - ideal (>= 0 componentwise except
    // exact zeros on the ideal carrier).
    let translated: Vec<Vec<f64>> = pop
        .iter()
        .map(|ind| (0..m).map(|j| ind.f[j] - ideal[j]).collect())
        .collect();
    // ASF extremes: for each axis j pick argmin_k max_l (T[k][l]/eps_l)
    // with eps_l = 1e-6 on l != j, eps_j = 1.0. Ties resolve to the
    // lowest population index (canonical).
    let eps_non_axis = 1.0e-6;
    let mut extremes: Vec<Option<usize>> = vec![None; m];
    for j in 0..m {
        let mut best: Option<(f64, usize)> = None;
        for (idx, t) in translated.iter().enumerate() {
            let asf = (0..m)
                .map(|l| {
                    if l == j {
                        t[l]
                    } else {
                        t[l] / eps_non_axis
                    }
                })
                .fold(f64::NEG_INFINITY, f64::max);
            best = match best {
                None => Some((asf, idx)),
                Some((bv, _)) if asf < bv || (asf == bv && idx < best.unwrap_or(idx)) => {
                    Some((asf, idx))
                }
                other => other,
            };
        }
        extremes[j] = best.map(|(_, i)| i);
    }
    // Intercept solve: given rows E_j = T(extreme_j) (translated),
    // find alpha with E^T alpha = ones; intercepts = 1/alpha.
    // Gaussian elimination with partial pivoting over translated rows;
    // any singularity/degeneracy routes to fallback.
    let mut singular_fallback = false;
    let mut scales = vec![1.0f64; m];
    'solve: {
        // Build matrix A[j][l] = T(extreme_j)[l].
        let mut a: Vec<Vec<f64>> = (0..m)
            .map(|j| match extremes[j] {
                Some(i) => translated[i].clone(),
                None => {
                    singular_fallback = true;
                    break 'solve;
                }
            })
            .collect();
        // Duplicate-row check: identical extremes guarantee singularity.
        for j in 0..m {
            for k2 in j + 1..m {
                if a[j] == a[k2] {
                    singular_fallback = true;
                    break 'solve;
                }
            }
        }
        // Solve A^T z = 1 via elimination on A^T (m x m), pivoting on
        // magnitude, absolute pivot floor relative to column maxima.
        let mut at = vec![vec![0.0f64; m]; m];
        for j in 0..m {
            for l in 0..m {
                at[l][j] = a[j][l];
            }
        }
        for col in 0..m {
            let mut piv = col;
            for r in col..m {
                if at[r][col].abs() > at[piv][col].abs() {
                    piv = r;
                }
            }
            if at[piv][col].abs() <= 1.0e-12 * at[col]
                .iter()
                .fold(0.0f64, |mx, &v| mx.max(v.abs()))
                .max(1.0e-300)
            {
                singular_fallback = true;
                break 'solve;
            }
            at.swap(col, piv);
            for r in col + 1..m {
                let factor = at[r][col] / at[col][col];
                for cc in col..m {
                    at[r][cc] -= factor * at[col][cc];
                }
            }
        }
        let mut rhs = vec![1.0f64; m];
        for r in (0..m).rev() {
            let mut s = rhs[r];
            for cc in r + 1..m {
                s -= at[r][cc] * rhs[cc];
            }
            rhs[r] = s / at[r][r];
        }
        // intercepts_l = 1/z_l must be positive.
        for (l, z) in rhs.iter().enumerate() {
            if !z.is_finite() || *z <= 0.0 {
                singular_fallback = true;
                break 'solve;
            }
            scales[l] = 1.0 / z;
        }
    }
    if singular_fallback {
        // Nadir fallback: scale by max translated extent per axis; a
        // degenerate axis (all-zero span) keeps scale 1.0 to stay finite.
        for l in 0..m {
            let span = translated.iter().map(|t| t[l]).fold(0.0f64, f64::max);
            scales[l] = if span > 0.0 && span.is_finite() {
                span
            } else {
                1.0
            };
        }
    }
    Normalization {
        ideal,
        scales,
        singular_fallback,
    }
}

/// Association distance of a translated point to reference direction r:
/// perpendicular component d_perp = |T' − (T'·r̂)r̂| in the NORMALIZED
/// space. Implemented exactly as the paper defines; both association
/// and niching consume this.
fn perpendicular_distance(t_norm: &[f64], dir: &[f64]) -> f64 {
    let dot: f64 = t_norm.iter().zip(dir.iter()).map(|(t, d)| t * d).sum();
    let sq: f64 = t_norm
        .iter()
        .zip(dir.iter())
        .map(|(t, d)| {
            let v = t - dot * d;
            v * v
        })
        .sum();
    sq.sqrt()
}

/// Pair up an individual index with its best (distance, reference)
/// association under canonical ties: smallest distance wins; equal-to-
/// bit distance prefers smaller reference index.
fn associate_one(
    t_norm: &[f64],
    dirs: &[Vec<f64>],
) -> (usize, f64) {
    let mut best: Option<(usize, f64)> = None;
    for (j, d) in dirs.iter().enumerate() {
        let dist = perpendicular_distance(t_norm, d);
        best = match best {
            None => Some((j, dist)),
            Some((bj, bd)) => {
                if dist < bd || (dist == bd && j < bj) {
                    Some((j, dist))
                } else {
                    Some((bj, bd))
                }
            }
        };
    }
    best.unwrap_or((0, 0.0))
}

/// Selection within the last partial front: fill remaining slots from
/// niches with the FEWEST already-selected members; within a niche,
/// smallest perpendicular distance; full ties prefer lower original
/// index. Mirrors Deb et al.'s niching with canonical resolution.
fn niche_fill(
    last_front: &[usize],
    slots: usize,
    assoc: &[(usize, f64)],
    niche_counts: &mut Vec<usize>,
) -> Vec<usize> {
    let mut chosen = Vec::with_capacity(slots);
    let mut pending: Vec<usize> = last_front.to_vec();
    // Order candidate scan canonically ONCE: by (niche_count, dist,
    // index) re-evaluated after every pick — implemented as repeated
    // argmin, O(slots * front), fine for baseline budgets.
    while chosen.len() < slots && !pending.is_empty() {
        let mut best_pick: Option<(usize, usize, f64)> = None; // (pending_pos, niche, dist)
        for (pos, &idx) in pending.iter().enumerate() {
            let (niche, dist) = assoc[idx];
            let key_cnt = niche_counts[niche];
            let cur = (key_cnt, niche, dist);
            let better = match best_pick {
                None => true,
                Some((_, bn, bd)) => {
                    let (_, bcnt,) = (
                        niche_counts[bn],
                        (),
                    );
                    let _ = bcnt;
                    // Compare against CURRENT stored best state below.
                    let best_niche_count_now = current_best_count(niche_counts, best_pick);
                    (key_cnt, niche, dist) < (best_niche_count_now, bn, bd)
                }
            };
            if better {
                let _ = cur;
                best_pick = Some((pos, niche, dist));
            }
        }
        if let Some((pos, niche, _)) = best_pick {
            let idx = pending.remove(pos);
            niche_counts[niche] += 1;
            chosen.push(idx);
        } else {
            break;
        }
    }
    chosen
}

// Helper kept separate so the nested comparator above stays readable.
fn current_best_count(
    niche_counts: &[usize],
    best_pick: Option<(usize, usize, f64)>,
) -> usize {
    best_pick.map_or(usize::MAX, |(_, niche, _)| niche_counts[niche])
}
