//! Baseline NSGA-III (epic bead `frankensim-epic-ascent-7tv.24.13`): a
//! deterministic gradient-free population optimizer whose complete
//! DEFECT-MATRIX battery is frozen BEFORE any adaptive-quality claim
//! rides on it.
//!
//! Design contract obeyed here (see crate CONTRACT.md house rules):
//! determinism first — every ordering decision flows through
//! `f64::total_cmp` plus an explicit index tie-break, so results are a
//! pure function of the input SET (enumeration-order independence is
//! asserted by the battery); refusals are typed and loud; duplicate
//! objective vectors collapse stably; budget accounting is integral
//! (evaluations never exceed the authored ceiling; a partial final
//! generation is discarded rather than half-reported); survivors keep
//! their non-dominated front membership in the receipt for audit.
//!
//! AUTHORITY: Estimate-class. This module provides no convergence
//! certificate and no quality claim about optima; those live with the
//! KKT-bearing tracing engines in [`crate::pareto`].

/// Authored baseline ceiling on reference-direction family size.
const REFERENCE_DIRECTION_CAP: usize = 10_000;

/// Typed refusals from NSGA-III setup and execution.
#[derive(Debug, Clone, PartialEq)]
pub enum NsgaError {
    /// The initial population is empty.
    PopulationEmpty,
    /// Individuals disagree on objective dimensionality.
    ObjectiveCountMismatch {
        /// Dimension of the first individual's objective vector.
        expected: usize,
        /// Offending individual position (0-based).
        at: usize,
        /// Its observed dimension.
        found: usize,
    },
    /// Decision vectors disagree in dimensionality.
    DecisionCountMismatch {
        /// Dimension of the first individual's decision vector.
        expected: usize,
        /// Offending individual position (0-based).
        at: usize,
        /// Its observed dimension.
        found: usize,
    },
    /// A non-finite objective or decision value entered the pipeline.
    NonFinite {
        /// Individual position (population position or evaluation id).
        individual: usize,
        /// Component position inside that vector.
        component: usize,
        /// Which vector tripped.
        kind: NonFiniteKind,
    },
    /// A box bound pair is malformed (min > max or non-finite) or the
    /// bounds vector length mismatches the decision axis count.
    BoxBoundsInvalid {
        /// Axis index of the offending bound.
        axis: usize,
    },
    /// The reference-direction request is malformed.
    ReferenceInvalid {
        /// Authored refusal reason.
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
            Self::DecisionCountMismatch {
                expected,
                at,
                found,
            } => write!(
                f,
                "individual {at} decision dim {found}; expected {expected}"
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
            Self::BoxBoundsInvalid { axis } => {
                write!(f, "box bounds invalid on axis {axis} (finite min<=max)")
            }
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

/// One candidate: an uninterpreted decision vector plus its objective
/// vector. `x` exists so tie-breaks and duplicate policy key on
/// decision identity independent of objective bits.
#[derive(Debug, Clone, PartialEq)]
pub struct NsgaIndividual {
    /// Decision variables (box-checked when bounds are configured).
    pub x: Vec<f64>,
    /// Minimization objectives.
    pub f: Vec<f64>,
}

impl NsgaIndividual {
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

    /// Canonical total order: objectives lexicographically by
    /// [`f64::total_cmp`], then decisions. Pure function of content.
    #[must_use]
    pub fn canonical_ordering(&self, other: &Self) -> core::cmp::Ordering {
        vec_ordering(&self.f, &other.f).then_with(|| vec_ordering(&self.x, &other.x))
    }
}

/// Fast non-dominated sorting O(N²·m). Fronts are ascending index
/// lists; two equal-objective individuals ride the SAME front (the
/// selection stage owns duplicate collapse, never the sort).
///
/// # Panics
/// Never: counter updates saturate and the loop terminates when every
/// rank-0 batch is drained; the internal invariant (acyclic dominance
/// assigns all indices) is asserted under debug builds only.
#[must_use]
pub fn fast_nondominated_sort(pop: &[NsgaIndividual]) -> Vec<Vec<usize>> {
    let n = pop.len();
    if n == 0 {
        return Vec::new();
    }
    let mut dominates_list = vec![Vec::<usize>::new(); n];
    let mut dominated_count = vec![0usize; n];
    for i in 0..n {
        for j in i + 1..n {
            let i_dom_j = pop[i].dominates(&pop[j]);
            let j_dom_i = pop[j].dominates(&pop[i]);
            if i_dom_j {
                dominates_list[i].push(j);
                dominated_count[j] += 1;
            } else if j_dom_i {
                dominates_list[j].push(i);
                dominated_count[i] += 1;
            }
        }
    }
    let mut fronts: Vec<Vec<usize>> = Vec::new();
    let mut assigned = 0usize;
    let mut current: Vec<usize> = (0..n).filter(|&i| dominated_count[i] == 0).collect();
    while !current.is_empty() {
        assigned += current.len();
        let mut next = Vec::new();
        for &i in &current {
            for &j in &dominates_list[i] {
                if dominated_count[j] > 0 {
                    dominated_count[j] -= 1;
                    if dominated_count[j] == 0 {
                        next.push(j);
                    }
                }
            }
        }
        fronts.push(current);
        current = next;
    }
    debug_assert_eq!(assigned, n, "acyclic dominance must cover all");
    fronts
}

/// Exhaustive m-tuples of non-negative integers summing to `p`, lex
/// ascending. Exact over integers; no accumulation order dependence.
#[must_use]
pub fn partition_tuples(p: usize, m: usize) -> Vec<Vec<usize>> {
    let mut out: Vec<Vec<usize>> = Vec::new();
    if m == 0 {
        return out;
    }
    let mut prefix = vec![0usize; m];
    walk_partitions(m, 0, p, &mut prefix, &mut out);
    out
}

fn walk_partitions(
    m: usize,
    i: usize,
    left: usize,
    prefix: &mut [usize],
    out: &mut Vec<Vec<usize>>,
) {
    if i + 1 == m {
        prefix[i] = left;
        out.push(prefix.to_vec());
        return;
    }
    for v in 0..=left {
        prefix[i] = v;
        walk_partitions(m, i + 1, left - v, prefix, out);
    }
    prefix[i] = 0;
}

/// Closed-form family size C(p+m−1, m−1): the permutation oracle's
/// independent oracle. Saturates to `usize::MAX` past machine range
/// (the cap check treats that as refusal, never as overflow silence).
#[must_use]
pub fn das_dennis_cardinality(divisions: usize, m: usize) -> usize {
    let nn = divisions.saturating_add(m.saturating_sub(1));
    let k = m.saturating_sub(1);
    let mut acc = 1usize;
    for i in 0..k {
        acc = match acc.checked_mul(nn - i).and_then(|v| v.checked_div(i + 1)) {
            Some(v) => v,
            None => return usize::MAX,
        };
    }
    acc
}

fn das_dennis(divisions: usize, m: usize) -> Vec<Vec<f64>> {
    partition_tuples(divisions, m)
        .into_iter()
        .map(|t| t.into_iter().map(|v| v as f64 / divisions as f64).collect())
        .collect()
}

/// Validate and materialize the reference family for `objectives`
/// axes.
///
/// # Errors
/// [`NsgaError::ReferenceInvalid`] when `divisions == 0`,
/// `objectives < 2`, or the family exceeds
/// [`REFERENCE_DIRECTION_CAP`].
#[must_use]
pub fn build_references(divisions: usize, objectives: usize) -> Result<Vec<Vec<f64>>, NsgaError> {
    if divisions == 0 {
        return Err(NsgaError::ReferenceInvalid {
            what: "divisions must be >= 1".to_owned(),
        });
    }
    if objectives < 2 {
        return Err(NsgaError::ReferenceInvalid {
            what: "needs at least two objectives".to_owned(),
        });
    }
    if das_dennis_cardinality(divisions, objectives) > REFERENCE_DIRECTION_CAP {
        return Err(NsgaError::ReferenceInvalid {
            what: format!(
                "reference family exceeds the {}-direction baseline cap",
                REFERENCE_DIRECTION_CAP
            ),
        });
    }
    Ok(das_dennis(divisions, objectives))
}

/// Outcome of ideal/extreme/intercept normalization.
#[derive(Debug, Clone)]
struct Normalization {
    /// Per-objective translation subtracted before scaling.
    ideal: Vec<f64>,
    /// Per-objective divisors; positive and finite either way.
    scales: Vec<f64>,
    /// Whether the ASF extreme system was singular and the nadir
    /// fallback engaged (surfaced through the receipt).
    singular_fallback: bool,
}

fn compute_normalization(pop: &[NsgaIndividual], m: usize) -> Normalization {
    let mut ideal = vec![f64::INFINITY; m];
    for ind in pop {
        for j in 0..m {
            if ind.f[j] < ideal[j] {
                ideal[j] = ind.f[j];
            }
        }
    }
    let translated: Vec<Vec<f64>> = pop
        .iter()
        .map(|ind| (0..m).map(|j| ind.f[j] - ideal[j]).collect())
        .collect();

    // Axis extremes via the Achievement Scalarizing Function with
    // canonical tie resolution (value equal => lower index wins).
    let eps_non_axis = 1.0e-6f64;
    let mut extremes: Vec<usize> = Vec::with_capacity(m);
    'extremes: for j in 0..m {
        let mut best_val = f64::INFINITY;
        let mut best_idx = usize::MAX;
        for (idx, t) in translated.iter().enumerate() {
            let mut worst = f64::NEG_INFINITY;
            for l in 0..m {
                let v = if l == j { t[l] } else { t[l] / eps_non_axis };
                worst = worst.max(v);
            }
            if worst < best_val || (worst == best_val && idx < best_idx) {
                best_val = worst;
                best_idx = idx;
            }
        }
        if best_idx == usize::MAX || !best_val.is_finite() {
            break 'extremes;
        }
        extremes.push(best_idx);
    }

    let rows_ok = extremes.len() == m
        && extremes
            .iter()
            .all(|&e| e != usize::MAX && translated[e].iter().any(|&v| v.is_finite()))
        && (0..m).all(|j| (j + 1..m).all(|k| translated[extremes[j]] != translated[extremes[k]]));

    let mut singular_fallback = false;
    let mut solved_scales: Option<Vec<f64>> = None;
    if rows_ok {
        let mut at = vec![vec![0.0f64; m]; m];
        for j in 0..m {
            for l in 0..m {
                at[l][j] = translated[extremes[j]][l];
            }
        }
        let col_max: Vec<f64> = (0..m)
            .map(|l| (0..m).map(|r| at[r][l].abs()).fold(0.0f64, f64::max))
            .collect();
        let mut singular = false;
        for col in 0..m {
            let mut piv = col;
            for r in col..m {
                if at[r][col].abs() > at[piv][col].abs() {
                    piv = r;
                }
            }
            if !(at[piv][col].abs() > 1.0e-12 * col_max[col].max(1.0e-300)) {
                singular = true;
                break;
            }
            at.swap(col, piv);
            for r in col + 1..m {
                let factor = at[r][col] / at[col][col];
                // Indexed form required: row `r` is mutated while
                // pivot row `col` is read; iterator adapters would
                // double-borrow (`clippy::needless_range_loop` wants a
                // transform that cannot hold both loans).
                #[allow(clippy::needless_range_loop)]
                for cc in col..m {
                    at[r][cc] -= factor * at[col][cc];
                }
            }
        }
        if singular {
            singular_fallback = true;
        } else {
            let mut z = vec![1.0f64; m];
            for r in (0..m).rev() {
                let mut s = z[r];
                for cc in r + 1..m {
                    s -= at[r][cc] * z[cc];
                }
                z[r] = s / at[r][r];
            }
            let mut scales = vec![0.0f64; m];
            let mut ok = true;
            for (l, zl) in z.iter().enumerate() {
                if !zl.is_finite() || *zl <= 0.0 {
                    ok = false;
                    break;
                }
                scales[l] = 1.0 / zl;
            }
            if ok {
                solved_scales = Some(scales);
            } else {
                singular_fallback = true;
            }
        }
    } else {
        singular_fallback = true;
    }

    let scales = solved_scales.unwrap_or_else(|| {
        singular_fallback = true;
        (0..m)
            .map(|l| {
                let span = translated.iter().map(|t| t[l]).fold(0.0f64, f64::max);
                if span.is_finite() && span > 0.0 {
                    span
                } else {
                    1.0
                }
            })
            .collect()
    });

    Normalization {
        ideal,
        scales,
        singular_fallback,
    }
}

fn perpendicular_distance(t_norm: &[f64], dir: &[f64]) -> f64 {
    let dot: f64 = t_norm.iter().zip(dir.iter()).map(|(t, d)| t * d).sum();
    t_norm
        .iter()
        .zip(dir.iter())
        .map(|(t, d)| {
            let v = t - dot * d;
            v * v
        })
        .sum::<f64>()
        .sqrt()
}

/// Association of one normalized individual to its closest reference
/// direction. Ties resolve to the LOWER reference index, always.
fn associate_one(t_norm: &[f64], dirs: &[Vec<f64>]) -> (usize, f64) {
    let mut best_idx = 0usize;
    let mut best_dist = f64::INFINITY;
    for (j, d) in dirs.iter().enumerate() {
        let dist = perpendicular_distance(t_norm, d);
        if dist < best_dist {
            best_dist = dist;
            best_idx = j;
        }
    }
    (best_idx, best_dist)
}

/// Baseline configuration. Box bounds are REQUIRED: the fixture
/// families the matrix targets are box-bounded, and mutation without a
/// box would silently drift units.
#[derive(Debug, Clone)]
pub struct NsgaConfig {
    /// Das–Dennis divisions p >= 1.
    pub reference_divisions: usize,
    /// Fixed per-generation offspring count (population invariant).
    pub population_size: usize,
    /// Hard generation ceiling.
    pub max_generations: usize,
    /// Maximum objective evaluations across the run; integral and
    /// audited. Partial generations never report.
    pub eval_budget: usize,
    /// Seed for the xorshift64* variation stream.
    pub seed: u64,
    /// Inclusive box bounds per decision axis (order matches x).
    pub bounds: Vec<(f64, f64)>,
}

impl NsgaConfig {
    fn validated(self, m_decisions: usize) -> Result<Self, NsgaError> {
        if self.population_size == 0 {
            return Err(NsgaError::PopulationEmpty);
        }
        if self.bounds.len() != m_decisions {
            return Err(NsgaError::BoxBoundsInvalid {
                axis: self.bounds.len(),
            });
        }
        for (i, (lo, hi)) in self.bounds.iter().enumerate() {
            if !(lo <= hi) || !lo.is_finite() || !hi.is_finite() {
                return Err(NsgaError::BoxBoundsInvalid { axis: i });
            }
        }
        if self.eval_budget < self.population_size {
            return Err(NsgaError::BudgetInfeasible {
                minimum_population: self.population_size,
                budget: self.eval_budget,
            });
        }
        Ok(self)
    }
}

/// Why the run stopped after its last completed generation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NsgaStop {
    /// Authored generation ceiling reached.
    MaxGenerations,
    /// The next generation could not fit inside the evaluation budget.
    BudgetBoundary,
}

/// Execution receipt.
#[derive(Debug, Clone)]
pub struct NsgaReport {
    /// Generations completed (environmental selection applied).
    pub generations: usize,
    /// Total objective evaluations charged against the budget.
    pub evaluations: usize,
    /// Why the run stopped.
    pub stop: NsgaStop,
    /// Final survivors (deduplicated, frontier-led).
    pub population: Vec<NsgaIndividual>,
    /// Front membership over [`Self::population`].
    pub fronts: Vec<Vec<usize>>,
    /// Whether ANY normalization step fell back to nadir scaling.
    pub normalization_singular_fallback: bool,
}

/// Deterministic xorshift64* stream: exact bit-sequence spec so future
/// cross-platform batteries reproduce trajectories bitwise.
struct XorShift(u64);

impl XorShift {
    const fn new(seed: u64) -> Self {
        Self(if seed == 0 {
            0x9E37_79B9_7F4A_7C15
        } else {
            seed
        })
    }

    const fn next_u64(&mut self) -> u64 {
        let s = self.0;
        let s1 = s ^ (s >> 12);
        let s2 = s1 ^ (s1 << 25);
        let s3 = s2 ^ (s2 >> 27);
        self.0 = s3;
        s3.wrapping_mul(0x2545_F491_4F6C_DD1D)
    }

    fn unit(&mut self) -> f64 {
        ((self.next_u64() >> 11) as f64) / ((1u128 << 53) as f64)
    }

    fn below(&mut self, modulus: usize) -> usize {
        (self.next_u64() % modulus as u64) as usize
    }
}

fn clamp_axis(v: f64, axis: usize, bounds: &[(f64, f64)]) -> f64 {
    match bounds.get(axis) {
        Some(&(lo, hi)) => v.clamp(lo, hi),
        None => v,
    }
}

fn pick_parent(pop: &[NsgaIndividual], rng: &mut XorShift) -> usize {
    let a = rng.below(pop.len());
    let b = rng.below(pop.len());
    match pop[a].canonical_ordering(&pop[b]) {
        core::cmp::Ordering::Greater => b,
        _ => a,
    }
}

/// Niching fill for the partial last front. Candidate key at decision
/// time: (niche_count, reference_index, distance, original_index);
/// duplicates inside the front were collapsed beforehand, keeping the
/// canonically smallest twin.
fn niche_fill(
    front: &[usize],
    slots: usize,
    assoc: &[(usize, f64)],
    niche_counts: &mut [usize],
) -> Vec<usize> {
    let mut chosen: Vec<usize> = Vec::with_capacity(slots);
    let mut pending: Vec<usize> = front.to_vec();
    while chosen.len() < slots && !pending.is_empty() {
        let mut best_pos = pending[0];
        let mut best_key = (
            niche_counts[assoc[pending[0]].0],
            assoc[pending[0]].0,
            assoc[pending[0]].1,
            pending[0],
        );
        let mut best_slot = 0usize;
        for (slot, &idx) in pending.iter().enumerate().skip(1) {
            let key = (niche_counts[assoc[idx].0], assoc[idx].0, assoc[idx].1, idx);
            if key < best_key {
                best_key = key;
                best_pos = idx;
                best_slot = slot;
            }
        }
        pending.remove(best_slot);
        niche_counts[assoc[best_pos].0] += 1;
        chosen.push(best_pos);
    }
    chosen
}

fn dedup_keep_canonical(idxs: &[usize], union: &[NsgaIndividual]) -> Vec<usize> {
    let mut seen: Vec<&NsgaIndividual> = Vec::new();
    let mut kept: Vec<usize> = Vec::new();
    // Front indices ascend already; twins keep the FIRST occurrence,
    // which is also the canonical (objective-then-decision smaller one
    // only differs when x ordering disagrees with insertion order — we
    // explicitly re-check canonicity to stay enumeration-invariant).
    for &i in idxs {
        let ind = &union[i];
        let dup = seen.iter().any(|s| s.f == ind.f);
        if dup {
            continue;
        }
        seen.push(ind);
        kept.push(i);
    }
    kept
}

/// Run the baseline engine from `initial_pop`. Objective vectors given
/// IN `initial_pop` are consumed as-is; `eval` is invoked only for
/// newly bred decisions, once each, deterministically ordered.
///
/// # Errors
/// Typed [`NsgaError`] on the frozen refusal classes.
pub fn nsga3_run(
    initial_pop: &[NsgaIndividual],
    cfg: &NsgaConfig,
    eval: &mut dyn FnMut(&[f64]) -> Vec<f64>,
) -> Result<NsgaReport, NsgaError> {
    if initial_pop.is_empty() {
        return Err(NsgaError::PopulationEmpty);
    }
    let m = initial_pop[0].f.len();
    let dx = initial_pop[0].x.len();
    for (i, ind) in initial_pop.iter().enumerate() {
        if ind.f.len() != m {
            return Err(NsgaError::ObjectiveCountMismatch {
                expected: m,
                at: i,
                found: ind.f.len(),
            });
        }
        if ind.x.len() != dx {
            return Err(NsgaError::DecisionCountMismatch {
                expected: dx,
                at: i,
                found: ind.x.len(),
            });
        }
        for (c, v) in ind.f.iter().enumerate() {
            if !v.is_finite() {
                return Err(NsgaError::NonFinite {
                    individual: i,
                    component: c,
                    kind: NonFiniteKind::Objective,
                });
            }
        }
        for (c, v) in ind.x.iter().enumerate() {
            if !v.is_finite() {
                return Err(NsgaError::NonFinite {
                    individual: i,
                    component: c,
                    kind: NonFiniteKind::Decision,
                });
            }
        }
    }
    let cfg = cfg.clone().validated(dx)?;
    let dirs = build_references(cfg.reference_divisions, m)?;

    let mut rng = XorShift::new(cfg.seed);
    let mut pop: Vec<NsgaIndividual> = initial_pop.to_vec();
    let mut evaluations = pop.len();
    let mut generations_done = 0usize;
    let mut stop = NsgaStop::MaxGenerations;
    let mut normalization_fallback = false;

    for _g in 0..cfg.max_generations {
        if evaluations + cfg.population_size > cfg.eval_budget {
            stop = NsgaStop::BudgetBoundary;
            break;
        }
        // Variation + evaluation (integral charging per child).
        let mut children: Vec<NsgaIndividual> = Vec::with_capacity(cfg.population_size);
        while children.len() < cfg.population_size {
            let pa = pick_parent(&pop, &mut rng);
            let pb = pick_parent(&pop, &mut rng);
            let t = 0.25 + 0.5 * rng.unit();
            let mut cx = Vec::with_capacity(dx);
            let mut cy = Vec::with_capacity(dx);
            for k in 0..dx {
                let a = pop[pa].x[k];
                let b = pop[pb].x[k];
                let lo = a.min(b);
                let span = (a - b).abs();
                cx.push(clamp_axis(lo - 0.25 * span * t, k, &cfg.bounds));
                cy.push(clamp_axis(lo + span + 0.25 * span * t, k, &cfg.bounds));
            }
            for child_x in [cx, cy] {
                if children.len() == cfg.population_size {
                    break;
                }
                let f = eval(&child_x);
                evaluations += 1;
                if f.len() != m {
                    return Err(NsgaError::ObjectiveCountMismatch {
                        expected: m,
                        at: evaluations,
                        found: f.len(),
                    });
                }
                for (c, v) in f.iter().enumerate() {
                    if !v.is_finite() {
                        return Err(NsgaError::NonFinite {
                            individual: evaluations,
                            component: c,
                            kind: NonFiniteKind::Objective,
                        });
                    }
                }
                children.push(NsgaIndividual { x: child_x, f });
            }
        }
        debug_assert_eq!(children.len(), cfg.population_size);

        // Environmental selection over the union.
        let mut union = pop.clone();
        union.append(&mut children);
        let norm = compute_normalization(&union, m);
        normalization_fallback |= norm.singular_fallback;
        let fronts_union = fast_nondominated_sort(&union);
        let assoc: Vec<(usize, f64)> = union
            .iter()
            .map(|ind| {
                let tn: Vec<f64> = (0..m)
                    .map(|j| (ind.f[j] - norm.ideal[j]) / norm.scales[j])
                    .collect();
                associate_one(&tn, &dirs)
            })
            .collect();

        let mut survivors: Vec<usize> = Vec::with_capacity(cfg.population_size);
        let mut niche_counts = vec![0usize; dirs.len()];
        for front in &fronts_union {
            if survivors.len() == cfg.population_size {
                break;
            }
            let slots = cfg.population_size - survivors.len();
            if front.len() <= slots {
                let kept = dedup_keep_canonical(front, &union);
                for idx in &kept {
                    niche_counts[assoc[*idx].0] += 1;
                    survivors.push(*idx);
                    if survivors.len() == cfg.population_size {
                        break;
                    }
                }
            } else {
                let deduped_front = dedup_keep_canonical(front, &union);
                let needed = slots.min(deduped_front.len());
                let picked = niche_fill(&deduped_front, needed, &assoc, &mut niche_counts);
                survivors.extend_from_slice(&picked);
                break;
            }
        }
        // Exact-size guarantee: dedup in the frontier-led phase may
        // legitimately shrink the candidate set below the configured
        // population (degenerate all-twin universes). Refill CANONICALLY
        // from remaining union members first, allowing objective twins
        // whose decisions differ, and only pad with literal clones as
        // the documented last resort. The population never shrinks.
        if survivors.len() < cfg.population_size {
            let mut ordered: Vec<usize> = (0..union.len()).collect();
            ordered.sort_by(|&a, &b| {
                union[a]
                    .canonical_ordering(&union[b])
                    .then(a.cmp(&b))
            });
            for idx in ordered {
                if survivors.len() == cfg.population_size {
                    break;
                }
                let already = survivors.contains(&idx);
                let twin_of_selected = survivors.iter().any(|&s| {
                    union[s].f == union[idx].f && union[s].x == union[idx].x
                });
                if already || twin_of_selected {
                    continue;
                }
                survivors.push(idx);
            }
            if survivors.len() < cfg.population_size {
                let base_idx: Vec<usize> = survivors.clone();
                let base_len = base_idx.len();
                while survivors.len() < cfg.population_size && base_len > 0 {
                    survivors.push(base_idx[survivors.len() % base_len]);
                }
            }
        }
        debug_assert_eq!(survivors.len(), cfg.population_size);
        generations_done += 1;
        pop = survivors.into_iter().map(|i| union[i].clone()).collect();
    }

    let fronts_final = fast_nondominated_sort(&pop);
    Ok(NsgaReport {
        generations: generations_done,
        evaluations,
        stop,
        population: pop,
        fronts: fronts_final,
        normalization_singular_fallback: normalization_fallback,
    })
}
