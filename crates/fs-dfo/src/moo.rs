//! Multi-objective machinery (plan §9.9): NSGA-II with DETERMINISTIC
//! tie-breaking everywhere (index order — bitwise-replayable runs),
//! exact hypervolume for m <= 4 (2D sweep + recursive exclusive
//! contributions above), Monte Carlo hypervolume beyond that,
//! least-contributor archives, NSGA-III/MOEA/D many-objective lanes,
//! and knee-point detection. Canonical empirical CVaR lives in
//! fs-robust and is re-exported by the crate root. All randomness flows
//! through fs-rand streams.

use fs_rand::StreamKey;

/// One evaluated individual.
#[derive(Debug, Clone)]
pub struct Individual {
    /// Decision vector.
    pub x: Vec<f64>,
    /// Objective vector (minimization).
    pub f: Vec<f64>,
}

fn check_variation_inputs(dim: usize, bounds: (f64, f64), eta_c: f64, eta_m: f64, p_mut: f64) {
    assert!(
        dim > 0,
        "box-search optimizers need at least one decision variable"
    );
    assert!(
        bounds.0.is_finite() && bounds.1.is_finite() && bounds.0 <= bounds.1,
        "box bounds must be finite and ordered"
    );
    assert!(
        eta_c.is_finite() && eta_c >= 0.0 && eta_m.is_finite() && eta_m >= 0.0,
        "SBX and mutation distribution indices must be finite and non-negative"
    );
    assert!(
        p_mut.is_finite() && (0.0..=1.0).contains(&p_mut),
        "mutation probability must be in [0, 1]"
    );
}

fn objective_dimension(pop: &[Individual], label: &str) -> usize {
    assert!(!pop.is_empty(), "{label} needs a non-empty population");
    let first = &pop[0];
    let m = first.f.len();
    assert!(m > 0, "objective vectors must not be empty");
    for ind in pop {
        assert_eq!(
            ind.f.len(),
            m,
            "all objective vectors in a population must have the same dimension"
        );
    }
    m
}

fn direction_dimension(directions: &[Vec<f64>], label: &str) -> usize {
    assert!(!directions.is_empty(), "{label} must not be empty");
    let m = directions[0].len();
    assert!(m > 0, "{label} vectors must not be empty");
    for dir in directions {
        assert_eq!(
            dir.len(),
            m,
            "{label} vectors must all have the same dimension"
        );
        assert!(
            dir.iter().all(|v| v.is_finite() && *v >= 0.0) && dir.iter().any(|v| *v > 0.0),
            "{label} vectors must be finite, non-negative, and not all zero"
        );
    }
    m
}

/// `a` Pareto-dominates `b` (minimization).
#[must_use]
pub fn dominates(a: &[f64], b: &[f64]) -> bool {
    assert_eq!(
        a.len(),
        b.len(),
        "dominance comparison requires matching objective dimensions"
    );
    assert!(
        !a.is_empty(),
        "dominance comparison needs at least one objective"
    );
    let mut strictly = false;
    for (ai, bi) in a.iter().zip(b) {
        if ai > bi {
            return false;
        }
        if ai < bi {
            strictly = true;
        }
    }
    strictly
}

/// Fast non-dominated sort: returns front index per individual
/// (0 = non-dominated). Deterministic (index-ordered scans).
#[must_use]
pub fn non_dominated_sort(pop: &[Individual]) -> Vec<usize> {
    let n = pop.len();
    if n == 0 {
        return Vec::new();
    }
    objective_dimension(pop, "non-dominated sort");
    let mut dominated_by = vec![0usize; n]; // count of dominators
    let mut dominates_list: Vec<Vec<usize>> = vec![Vec::new(); n];
    for i in 0..n {
        for j in 0..n {
            if i != j && dominates(&pop[i].f, &pop[j].f) {
                dominates_list[i].push(j);
            }
        }
    }
    for (j, count) in dominated_by.iter_mut().enumerate() {
        *count = (0..n)
            .filter(|&i| i != j && dominates(&pop[i].f, &pop[j].f))
            .count();
    }
    let mut front = vec![usize::MAX; n];
    let mut current: Vec<usize> = (0..n).filter(|&i| dominated_by[i] == 0).collect();
    let mut level = 0usize;
    while !current.is_empty() {
        let mut next = Vec::new();
        for &i in &current {
            front[i] = level;
        }
        for &i in &current {
            for &j in &dominates_list[i] {
                dominated_by[j] -= 1;
                if dominated_by[j] == 0 {
                    next.push(j);
                }
            }
        }
        next.sort_unstable();
        next.dedup();
        current = next;
        level += 1;
    }
    front
}

/// Crowding distance within one front (boundary points get ∞).
#[must_use]
pub fn crowding_distance(front: &[&Individual]) -> Vec<f64> {
    let n = front.len();
    if n == 0 {
        return Vec::new();
    }
    let m = front[0].f.len();
    let mut dist = vec![0.0f64; n];
    for obj in 0..m {
        let mut order: Vec<usize> = (0..n).collect();
        order.sort_by(|&a, &b| front[a].f[obj].total_cmp(&front[b].f[obj]).then(a.cmp(&b)));
        let lo = front[order[0]].f[obj];
        let hi = front[order[n - 1]].f[obj];
        let span = (hi - lo).max(1e-30);
        dist[order[0]] = f64::INFINITY;
        dist[order[n - 1]] = f64::INFINITY;
        for w in 1..n - 1 {
            dist[order[w]] += (front[order[w + 1]].f[obj] - front[order[w - 1]].f[obj]) / span;
        }
    }
    dist
}

/// NSGA-II configuration.
#[derive(Debug, Clone, Copy)]
pub struct NsgaParams {
    /// Population size (even).
    pub pop: usize,
    /// Generations.
    pub generations: usize,
    /// SBX distribution index.
    pub eta_c: f64,
    /// Polynomial-mutation distribution index.
    pub eta_m: f64,
    /// Per-variable mutation probability.
    pub p_mut: f64,
    /// Root seed for all derived study streams.
    pub seed: u64,
}

/// Run NSGA-II on a box-bounded problem. Returns the final
/// population's first front. Deterministic per seed.
pub fn nsga2(
    objectives: &mut dyn FnMut(&[f64]) -> Vec<f64>,
    dim: usize,
    bounds: (f64, f64),
    params: &NsgaParams,
) -> Vec<Individual> {
    assert!(params.pop > 0, "NSGA-II population must be non-empty");
    check_variation_inputs(dim, bounds, params.eta_c, params.eta_m, params.p_mut);
    let (lo, hi) = bounds;
    let mut stream = StreamKey {
        seed: params.seed,
        kernel: 0x05A2,
        tile: 0,
    }
    .stream();
    let mut pop: Vec<Individual> = (0..params.pop)
        .map(|_| {
            let x: Vec<f64> = (0..dim)
                .map(|_| (hi - lo).mul_add(stream.next_f64(), lo))
                .collect();
            let f = objectives(&x);
            Individual { x, f }
        })
        .collect();
    for _ in 0..params.generations {
        // Selection pool: rank + crowding over the current population.
        let fronts = non_dominated_sort(&pop);
        let crowd = crowding_for_population(&pop, &fronts);
        // Binary tournaments → offspring.
        let mut offspring = Vec::with_capacity(params.pop);
        while offspring.len() < params.pop {
            let pick = |s: &mut fs_rand::Stream| -> usize {
                let a = s.next_below(pop.len() as u64) as usize;
                let b = s.next_below(pop.len() as u64) as usize;
                if fronts[a] < fronts[b]
                    || (fronts[a] == fronts[b] && crowd[a] > crowd[b])
                    || (fronts[a] == fronts[b]
                        && crowd[a].to_bits() == crowd[b].to_bits()
                        && a <= b)
                {
                    a
                } else {
                    b
                }
            };
            let p1 = pick(&mut stream);
            let p2 = pick(&mut stream);
            let (mut c1, mut c2) = sbx(&pop[p1].x, &pop[p2].x, params.eta_c, lo, hi, &mut stream);
            mutate(&mut c1, params.eta_m, params.p_mut, lo, hi, &mut stream);
            mutate(&mut c2, params.eta_m, params.p_mut, lo, hi, &mut stream);
            let f1 = objectives(&c1);
            offspring.push(Individual { x: c1, f: f1 });
            if offspring.len() < params.pop {
                let f2 = objectives(&c2);
                offspring.push(Individual { x: c2, f: f2 });
            }
        }
        // Environmental selection over parents ∪ offspring.
        pop.extend(offspring);
        let fronts = non_dominated_sort(&pop);
        let crowd = crowding_for_population(&pop, &fronts);
        let mut order: Vec<usize> = (0..pop.len()).collect();
        order.sort_by(|&a, &b| {
            fronts[a]
                .cmp(&fronts[b])
                .then(crowd[b].total_cmp(&crowd[a]))
                .then(a.cmp(&b))
        });
        order.truncate(params.pop);
        order.sort_unstable();
        pop = order.into_iter().map(|i| pop[i].clone()).collect();
    }
    let fronts = non_dominated_sort(&pop);
    pop.into_iter()
        .zip(fronts)
        .filter(|(_, r)| *r == 0)
        .map(|(ind, _)| ind)
        .collect()
}

fn crowding_for_population(pop: &[Individual], fronts: &[usize]) -> Vec<f64> {
    let mut crowd = vec![0.0f64; pop.len()];
    let max_front = fronts.iter().copied().max().unwrap_or(0);
    for level in 0..=max_front {
        let idx: Vec<usize> = (0..pop.len()).filter(|&i| fronts[i] == level).collect();
        if idx.is_empty() {
            continue;
        }
        if idx.len() <= 2 {
            for &i in &idx {
                crowd[i] = f64::INFINITY;
            }
            continue;
        }
        let members: Vec<&Individual> = idx.iter().map(|&i| &pop[i]).collect();
        let d = crowding_distance(&members);
        for (k, &i) in idx.iter().enumerate() {
            crowd[i] = d[k];
        }
    }
    crowd
}

/// Simulated binary crossover (SBX).
fn sbx(
    p1: &[f64],
    p2: &[f64],
    eta: f64,
    lo: f64,
    hi: f64,
    stream: &mut fs_rand::Stream,
) -> (Vec<f64>, Vec<f64>) {
    let mut c1 = p1.to_vec();
    let mut c2 = p2.to_vec();
    for i in 0..p1.len() {
        if stream.next_f64() < 0.9 {
            let u = stream.next_f64();
            let beta = if u <= 0.5 {
                fs_math::det::pow(2.0 * u, 1.0 / (eta + 1.0))
            } else {
                fs_math::det::pow(1.0 / (2.0 * (1.0 - u)), 1.0 / (eta + 1.0))
            };
            let a = f64::midpoint((1.0 + beta) * p1[i], (1.0 - beta) * p2[i]);
            let b = f64::midpoint((1.0 - beta) * p1[i], (1.0 + beta) * p2[i]);
            c1[i] = a.clamp(lo, hi);
            c2[i] = b.clamp(lo, hi);
        }
    }
    (c1, c2)
}

/// Polynomial mutation.
fn mutate(x: &mut [f64], eta: f64, p_mut: f64, lo: f64, hi: f64, stream: &mut fs_rand::Stream) {
    for xi in x.iter_mut() {
        if stream.next_f64() < p_mut {
            let u = stream.next_f64();
            let delta = if u < 0.5 {
                fs_math::det::pow(2.0 * u, 1.0 / (eta + 1.0)) - 1.0
            } else {
                1.0 - fs_math::det::pow(2.0 * (1.0 - u), 1.0 / (eta + 1.0))
            };
            *xi = (delta * (hi - lo) + *xi).clamp(lo, hi);
        }
    }
}

/// Exact hypervolume of a minimization front w.r.t. a reference point
/// (all points must dominate `reference`). 2D: dominated-strip sweep.
/// Higher m (≤ 4 intended): recursive exclusive contributions.
#[must_use]
pub fn hypervolume(front: &[Vec<f64>], reference: &[f64]) -> f64 {
    if reference.is_empty() {
        return 0.0;
    }
    let pts: Vec<Vec<f64>> = front
        .iter()
        .filter(|p| p.len() == reference.len() && p.iter().zip(reference).all(|(a, r)| a < r))
        .cloned()
        .collect();
    if pts.is_empty() {
        return 0.0;
    }
    hv_recursive(&pts, reference)
}

fn hv_recursive(pts: &[Vec<f64>], reference: &[f64]) -> f64 {
    let m = reference.len();
    if m == 1 {
        let best = pts
            .iter()
            .map(|p| p[0])
            .min_by(f64::total_cmp)
            .unwrap_or(reference[0]);
        return (reference[0] - best).max(0.0);
    }
    if m == 2 {
        // Sort by f1 ascending, deterministic tie-break on f2.
        let mut sorted = pts.to_vec();
        sorted.sort_by(|a, b| a[0].total_cmp(&b[0]).then(a[1].total_cmp(&b[1])));
        let mut hv = 0.0f64;
        let mut best_f2 = reference[1];
        for p in &sorted {
            if p[1] < best_f2 {
                hv += (reference[0] - p[0]) * (best_f2 - p[1]);
                best_f2 = p[1];
            }
        }
        return hv;
    }
    // Inclusion-exclusion by slicing on the last objective (WFG-style
    // exclusive contributions; exponential worst case — m ≤ 4 scope).
    let mut order: Vec<usize> = (0..pts.len()).collect();
    order.sort_by(|&a, &b| pts[a][m - 1].total_cmp(&pts[b][m - 1]).then(a.cmp(&b)));
    let mut hv = 0.0f64;
    for (k, &i) in order.iter().enumerate() {
        let z_hi = if k + 1 < order.len() {
            pts[order[k + 1]][m - 1]
        } else {
            reference[m - 1]
        };
        let z_lo = pts[i][m - 1];
        if z_hi <= z_lo {
            continue;
        }
        // Slab [z_lo, z_hi): points active are those with f_m ≤ z_lo.
        let active: Vec<Vec<f64>> = order[..=k]
            .iter()
            .map(|&j| pts[j][..m - 1].to_vec())
            .collect();
        hv += (z_hi - z_lo) * hv_recursive(&active, &reference[..m - 1]);
    }
    hv
}

/// Knee point of a 2-objective front: the member with maximum
/// perpendicular distance to the chord between the two extremes
/// (objectives normalized to [0,1] first). Returns the index into
/// `front`.
#[must_use]
pub fn knee_point(front: &[Vec<f64>]) -> usize {
    let n = front.len();
    assert!(n >= 3, "knee needs at least 3 points");
    // The perpendicular-distance formula reads p[0] and p[1]; a front of
    // <2-objective points would index out of bounds. Fail closed with the
    // crate's structured-panic error model rather than panic on `p[1]`.
    assert!(
        front.iter().all(|p| p.len() >= 2),
        "knee needs 2-objective points"
    );
    let (mut lo0, mut hi0) = (f64::INFINITY, f64::NEG_INFINITY);
    let (mut lo1, mut hi1) = (f64::INFINITY, f64::NEG_INFINITY);
    for p in front {
        lo0 = lo0.min(p[0]);
        hi0 = hi0.max(p[0]);
        lo1 = lo1.min(p[1]);
        hi1 = hi1.max(p[1]);
    }
    let norm = |p: &Vec<f64>| -> (f64, f64) {
        (
            (p[0] - lo0) / (hi0 - lo0).max(1e-30),
            (p[1] - lo1) / (hi1 - lo1).max(1e-30),
        )
    };
    // Extremes: best f0 and best f1 (normalized (0, y) and (x, 0)
    // corners of the front).
    let a = (0.0f64, 1.0f64);
    let b = (1.0f64, 0.0f64);
    let mut best = 0usize;
    let mut best_d = f64::NEG_INFINITY;
    for (i, p) in front.iter().enumerate() {
        let (x, y) = norm(p);
        // Distance to the line a—b: |(b0−a0)(a1−y) − (a0−x)(b1−a1)|/‖b−a‖.
        let d = ((b.0 - a.0) * (a.1 - y) - (a.0 - x) * (b.1 - a.1)).abs()
            / fs_math::det::sqrt((b.0 - a.0).powi(2) + (b.1 - a.1).powi(2));
        if d > best_d {
            best_d = d;
            best = i;
        }
    }
    best
}

/// Monte Carlo hypervolume estimate for m > 4 (where exact recursion
/// is exponential): scrambled-Sobol points in the reference box,
/// counting the dominated fraction. Deterministic per seed; the
/// standard error scales as √(p(1−p)/n) and the SAMPLE COUNT is the
/// caller's accuracy knob (reported honestly by returning the hit
/// count alongside).
#[must_use]
pub fn mc_hypervolume(
    front: &[Vec<f64>],
    reference: &[f64],
    samples: usize,
    seed: u64,
) -> (f64, usize) {
    let m = reference.len();
    if m == 0 {
        return (0.0, 0);
    }
    assert!(samples > 0, "MC hypervolume needs at least one sample");
    assert!(
        u32::try_from(samples).is_ok(),
        "MC hypervolume sample count must fit the Sobol index range"
    );
    let pts: Vec<&Vec<f64>> = front
        .iter()
        .filter(|p| p.len() == m && p.iter().zip(reference).all(|(a, r)| a < r))
        .collect();
    if pts.is_empty() {
        return (0.0, 0);
    }
    // Box lower corner: componentwise min of the front (nothing below
    // it is dominated-relevant beyond the points themselves).
    let mut lo = vec![f64::INFINITY; m];
    for p in &pts {
        for (l, v) in lo.iter_mut().zip(p.iter()) {
            *l = l.min(*v);
        }
    }
    let vol_box: f64 = lo
        .iter()
        .zip(reference)
        .map(|(l, r)| (r - l).max(0.0))
        .product();
    if vol_box <= 0.0 {
        return (0.0, 0);
    }
    let kq = m.min(fs_rand::qmc::MAX_SOBOL_DIM);
    let sobol = fs_rand::qmc::Sobol::scrambled(kq, seed);
    let mut tail = StreamKey {
        seed,
        kernel: 0x0871,
        tile: 0,
    }
    .stream();
    let mut pt = vec![0.0f64; kq];
    let mut hits = 0usize;
    let mut y = vec![0.0f64; m];
    for s in 0..samples {
        sobol.point((s + 1) as u32, &mut pt);
        for d in 0..m {
            let u = if d < kq { pt[d] } else { tail.next_f64() };
            y[d] = (reference[d] - lo[d]).mul_add(u, lo[d]);
        }
        if pts.iter().any(|p| p.iter().zip(&y).all(|(a, b)| a <= b)) {
            hits += 1;
        }
    }
    (vol_box * hits as f64 / samples as f64, hits)
}

/// Bounded Pareto archive with least-hypervolume-contributor
/// eviction (exact contributions via [`hypervolume`] — the m ≤ 4
/// regime; MC-contribution eviction joins with its consumer).
pub struct HvArchive {
    /// Archive members (mutually non-dominated).
    pub members: Vec<Individual>,
    /// Capacity.
    pub capacity: usize,
    /// Reference point.
    pub reference: Vec<f64>,
}

impl HvArchive {
    /// Empty archive.
    #[must_use]
    pub fn new(capacity: usize, reference: Vec<f64>) -> HvArchive {
        assert!(capacity >= 2, "an archive below 2 keeps no front");
        HvArchive {
            members: Vec::new(),
            capacity,
            reference,
        }
    }

    /// Current archive hypervolume.
    #[must_use]
    pub fn hv(&self) -> f64 {
        let pts: Vec<Vec<f64>> = self.members.iter().map(|i| i.f.clone()).collect();
        hypervolume(&pts, &self.reference)
    }

    /// Insert: dominated candidates are a NO-OP (returns false);
    /// otherwise dominated members are evicted, the candidate joins,
    /// and if over capacity the LEAST exclusive contributor leaves
    /// (deterministic index tie-break).
    pub fn insert(&mut self, cand: Individual) -> bool {
        if self
            .members
            .iter()
            .any(|mem| dominates(&mem.f, &cand.f) || mem.f == cand.f)
        {
            return false;
        }
        self.members.retain(|mem| !dominates(&cand.f, &mem.f));
        self.members.push(cand);
        if self.members.len() > self.capacity {
            let pts: Vec<Vec<f64>> = self.members.iter().map(|i| i.f.clone()).collect();
            let full = hypervolume(&pts, &self.reference);
            let mut worst = (f64::INFINITY, 0usize);
            for k in 0..pts.len() {
                let rest: Vec<Vec<f64>> = pts
                    .iter()
                    .enumerate()
                    .filter(|&(j, _)| j != k)
                    .map(|(_, p)| p.clone())
                    .collect();
                let contribution = full - hypervolume(&rest, &self.reference);
                if contribution < worst.0 {
                    worst = (contribution, k);
                }
            }
            self.members.remove(worst.1);
        }
        true
    }
}

/// Das–Dennis systematic reference directions on the (m−1)-simplex
/// with `p` divisions: all non-negative integer m-tuples summing to
/// p, scaled by 1/p. Deterministic lexicographic enumeration;
/// C(p+m−1, m−1) directions.
#[must_use]
pub fn das_dennis(m: usize, p: usize) -> Vec<Vec<f64>> {
    assert!(m > 0, "Das-Dennis directions need at least one objective");
    assert!(p > 0, "Das-Dennis directions need at least one division");
    let mut out = Vec::new();
    let mut stack = vec![(Vec::<usize>::new(), p)];
    while let Some((prefix, rest)) = stack.pop() {
        if prefix.len() == m - 1 {
            let mut dir: Vec<f64> = prefix.iter().map(|&k| k as f64 / p as f64).collect();
            dir.push(rest as f64 / p as f64);
            out.push(dir);
            continue;
        }
        // Push in reverse so the popped order is lexicographic.
        for k in (0..=rest).rev() {
            let mut next = prefix.clone();
            next.push(k);
            stack.push((next, rest - k));
        }
    }
    out
}

/// Perpendicular distance from a (normalized) objective vector to the
/// ray through a reference direction.
fn perp_distance(f: &[f64], dir: &[f64]) -> f64 {
    let dd: f64 = dir.iter().map(|d| d * d).sum();
    if dd < 1e-30 {
        return f.iter().map(|v| v * v).sum::<f64>().sqrt();
    }
    let t: f64 = f.iter().zip(dir).map(|(a, b)| a * b).sum::<f64>() / dd;
    f.iter()
        .zip(dir)
        .map(|(a, b)| {
            let r = t.mul_add(-b, *a);
            r * r
        })
        .sum::<f64>()
        .sqrt()
}

/// NSGA-III for many-objective minimization: reference-direction
/// niching replaces crowding distance. Normalization uses
/// achievement-scalarizing-function extreme points and the resulting
/// hyperplane intercepts, with the legacy floored first-front maxima as
/// a deterministic fail-closed fallback. Deterministic throughout
/// (index tie-breaks). The cubic intercept solve is bounded to 64
/// objectives; larger objective vectors retain the legacy normalization.
pub fn nsga3(
    objectives: &mut dyn FnMut(&[f64]) -> Vec<f64>,
    dim: usize,
    bounds: (f64, f64),
    directions: &[Vec<f64>],
    params: &NsgaParams,
) -> Vec<Individual> {
    assert!(params.pop > 0, "NSGA-III population must be non-empty");
    let direction_m = direction_dimension(directions, "NSGA-III reference directions");
    check_variation_inputs(dim, bounds, params.eta_c, params.eta_m, params.p_mut);
    let (lo, hi) = bounds;
    let mut stream = StreamKey {
        seed: params.seed,
        kernel: 0x05A3,
        tile: 0,
    }
    .stream();
    let mut pop: Vec<Individual> = (0..params.pop)
        .map(|_| {
            let x: Vec<f64> = (0..dim)
                .map(|_| (hi - lo).mul_add(stream.next_f64(), lo))
                .collect();
            let f = objectives(&x);
            Individual { x, f }
        })
        .collect();
    let objective_m = objective_dimension(&pop, "NSGA-III");
    assert_eq!(
        direction_m, objective_m,
        "NSGA-III reference-direction dimension must match objective dimension"
    );
    for _ in 0..params.generations {
        let fronts = non_dominated_sort(&pop);
        // Rank-only binary tournament (NSGA-III drops crowding).
        let mut offspring = Vec::with_capacity(params.pop);
        while offspring.len() < params.pop {
            let pick = |s: &mut fs_rand::Stream| -> usize {
                let a = s.next_below(pop.len() as u64) as usize;
                let b = s.next_below(pop.len() as u64) as usize;
                if fronts[a] < fronts[b] || (fronts[a] == fronts[b] && a <= b) {
                    a
                } else {
                    b
                }
            };
            let p1 = pick(&mut stream);
            let p2 = pick(&mut stream);
            let (mut c1, mut c2) = sbx(&pop[p1].x, &pop[p2].x, params.eta_c, lo, hi, &mut stream);
            mutate(&mut c1, params.eta_m, params.p_mut, lo, hi, &mut stream);
            mutate(&mut c2, params.eta_m, params.p_mut, lo, hi, &mut stream);
            let f1 = objectives(&c1);
            offspring.push(Individual { x: c1, f: f1 });
            if offspring.len() < params.pop {
                let f2 = objectives(&c2);
                offspring.push(Individual { x: c2, f: f2 });
            }
        }
        pop.extend(offspring);
        pop = nsga3_select(&pop, directions, params.pop);
    }
    let fronts = non_dominated_sort(&pop);
    pop.into_iter()
        .zip(fronts)
        .filter(|(_, r)| *r == 0)
        .map(|(ind, _)| ind)
        .collect()
}

/// NSGA-III environmental selection to `target` members.
fn nsga3_select(pop: &[Individual], directions: &[Vec<f64>], target: usize) -> Vec<Individual> {
    let m = objective_dimension(pop, "NSGA-III selection");
    assert_eq!(
        direction_dimension(directions, "NSGA-III reference directions"),
        m,
        "NSGA-III reference-direction dimension must match objective dimension"
    );
    let fronts = non_dominated_sort(pop);
    let max_front = fronts.iter().copied().max().unwrap_or(0);
    let mut accepted: Vec<usize> = Vec::with_capacity(target);
    let mut partial: Vec<usize> = Vec::new();
    for level in 0..=max_front {
        let members: Vec<usize> = (0..pop.len()).filter(|&i| fronts[i] == level).collect();
        if accepted.len() + members.len() <= target {
            accepted.extend(members);
        } else {
            partial = members;
            break;
        }
        if accepted.len() == target {
            break;
        }
    }
    if accepted.len() == target || partial.is_empty() {
        return accepted.into_iter().map(|i| pop[i].clone()).collect();
    }
    // Normalize over accepted ∪ partial: ideal point plus ASF extreme
    // points and hyperplane intercepts, with the legacy floored-maxima
    // behavior as the fail-closed fallback.
    let considered: Vec<usize> = accepted.iter().chain(&partial).copied().collect();
    let (ideal, span) = nsga3_normalization(pop, &fronts, &considered, m);
    let associate = |i: usize| nsga3_association(pop, directions, &ideal, &span, i);
    let mut niche = vec![0usize; directions.len()];
    for &i in &accepted {
        niche[associate(i).0] += 1;
    }
    let mut pool: Vec<(usize, usize, f64)> = partial
        .iter()
        .map(|&i| {
            let (k, d) = associate(i);
            (i, k, d)
        })
        .collect();
    while accepted.len() < target {
        // Direction with members in the pool and MINIMUM niche count
        // (index tie-break).
        let mut best_dir = usize::MAX;
        let mut best_count = usize::MAX;
        for &(_, k, _) in &pool {
            if niche[k] < best_count {
                best_count = niche[k];
                best_dir = k;
            }
        }
        // Among pool members of that direction: min perpendicular
        // distance if the niche is empty, else first by index.
        let mut chosen = usize::MAX;
        let mut chosen_pos = usize::MAX;
        let mut chosen_d = f64::INFINITY;
        for (pos, &(i, k, d)) in pool.iter().enumerate() {
            if k != best_dir {
                continue;
            }
            let better = if best_count == 0 {
                d < chosen_d || (d.to_bits() == chosen_d.to_bits() && i < chosen)
            } else {
                i < chosen
            };
            if better {
                chosen = i;
                chosen_pos = pos;
                chosen_d = d;
            }
        }
        accepted.push(chosen);
        pool.swap_remove(chosen_pos);
        niche[best_dir] += 1;
        if pool.is_empty() {
            break;
        }
    }
    accepted.into_iter().map(|i| pop[i].clone()).collect()
}

fn nsga3_normalization(
    pop: &[Individual],
    fronts: &[usize],
    considered: &[usize],
    m: usize,
) -> (Vec<f64>, Vec<f64>) {
    const MIN_SPAN: f64 = 1e-30;
    const MAX_HYPERPLANE_ENTRIES: usize = 4_096;

    let mut ideal = vec![f64::INFINITY; m];
    for &i in considered {
        for (d, v) in ideal.iter_mut().zip(&pop[i].f) {
            *d = d.min(*v);
        }
    }
    let mut fallback_span = vec![MIN_SPAN; m];
    for &i in considered {
        if fronts[i] == 0 {
            for (d, (v, id)) in fallback_span.iter_mut().zip(pop[i].f.iter().zip(&ideal)) {
                *d = d.max(v - id);
            }
        }
    }

    let matrix_entries = match m.checked_mul(m) {
        Some(entries) if entries <= MAX_HYPERPLANE_ENTRIES => entries,
        _ => return (ideal, fallback_span),
    };
    if matrix_entries == 0 {
        return (ideal, fallback_span);
    }
    let Some(extreme_indices) = nsga3_extreme_indices(pop, considered, &ideal, &fallback_span, m)
    else {
        return (ideal, fallback_span);
    };
    let Some(span) = nsga3_hyperplane_span(pop, &extreme_indices, &ideal, &fallback_span) else {
        return (ideal, fallback_span);
    };
    (ideal, span)
}

/// Select one considered extreme per objective with a scale-equilibrated
/// Deb-Jain-style ASF.
/// Translated objectives are first divided by the legacy floored
/// first-front-max fallback
/// spans. Algebraically this is equivariant under positive objective
/// scaling above the documented span floor; rounding and ASF tie
/// boundaries remain explicit no-claim cases. Off-axis objectives
/// receive a `1e-6` weight. Iteration and explicit index tie-breaking
/// make equal ASF values replayable.
fn nsga3_extreme_indices(
    pop: &[Individual],
    considered: &[usize],
    ideal: &[f64],
    fallback_span: &[f64],
    m: usize,
) -> Option<Vec<usize>> {
    const OFF_AXIS_WEIGHT: f64 = 1e-6;

    if considered.len() < m
        || fallback_span.len() != m
        || fallback_span
            .iter()
            .any(|span| !span.is_finite() || *span <= 0.0)
    {
        return None;
    }

    let mut extremes = Vec::with_capacity(m);
    for axis in 0..m {
        extremes.push(nsga3_asf_extreme_index(
            pop,
            considered,
            ideal,
            fallback_span,
            axis,
            OFF_AXIS_WEIGHT,
        )?);
    }

    for i in 0..extremes.len() {
        if extremes[..i].contains(&extremes[i]) {
            return None;
        }
    }
    Some(extremes)
}

fn nsga3_asf_extreme_index(
    pop: &[Individual],
    considered: &[usize],
    ideal: &[f64],
    fallback_span: &[f64],
    axis: usize,
    off_axis_weight: f64,
) -> Option<usize> {
    let mut best_index = usize::MAX;
    let mut best_asf = f64::INFINITY;
    for &i in considered {
        let mut asf = f64::NEG_INFINITY;
        for objective in 0..ideal.len() {
            let translated = pop[i].f[objective] - ideal[objective];
            if !translated.is_finite() || translated < 0.0 {
                return None;
            }
            let normalized = if translated == 0.0 {
                0.0
            } else {
                translated / fallback_span[objective]
            };
            if !normalized.is_finite() {
                return None;
            }
            let weight = if objective == axis {
                1.0
            } else {
                off_axis_weight
            };
            let component = normalized / weight;
            if !component.is_finite() {
                return None;
            }
            asf = asf.max(component);
        }
        let order = asf.total_cmp(&best_asf);
        if order.is_lt() || (order.is_eq() && i < best_index) {
            best_index = i;
            best_asf = asf;
        }
    }
    (best_index != usize::MAX).then_some(best_index)
}

/// Derive objective spans from the hyperplane through the ASF extreme
/// points. Columns are scaled by the maxima fallback before elimination;
/// this is algebraically equivalent to the unscaled system and improves
/// scale equivariance above the floor without claiming bitwise invariance
/// at rounding or tie boundaries.
fn nsga3_hyperplane_span(
    pop: &[Individual],
    extreme_indices: &[usize],
    ideal: &[f64],
    fallback_span: &[f64],
) -> Option<Vec<f64>> {
    const MAX_HYPERPLANE_ENTRIES: usize = 4_096;

    let m = ideal.len();
    if extreme_indices.len() != m
        || fallback_span.len() != m
        || fallback_span
            .iter()
            .any(|span| !span.is_finite() || *span <= 0.0)
    {
        return None;
    }
    match m.checked_mul(m) {
        Some(1..=MAX_HYPERPLANE_ENTRIES) => {}
        _ => return None,
    }

    let mut matrix = vec![vec![0.0; m]; m];
    for (row, &i) in extreme_indices.iter().enumerate() {
        for objective in 0..m {
            let translated = pop[i].f[objective] - ideal[objective];
            let scaled = translated / fallback_span[objective];
            if !translated.is_finite() || translated < 0.0 || !scaled.is_finite() {
                return None;
            }
            matrix[row][objective] = scaled;
        }
    }

    let coefficients = nsga3_solve_hyperplane(&matrix)?;
    let mut span = Vec::with_capacity(m);
    for (fallback, coefficient) in fallback_span.iter().zip(coefficients) {
        if !coefficient.is_finite() || coefficient <= 0.0 {
            return None;
        }
        let intercept = fallback / coefficient;
        if !intercept.is_finite() || intercept <= 0.0 {
            return None;
        }
        span.push(intercept);
    }
    Some(span)
}

/// Solve `A x = 1` with deterministic partial-pivot LU factorization.
/// The matrix is already column-scaled. Small relative pivots, a large
/// one-norm condition estimate, non-finite arithmetic, or a poor
/// scale-relative residual refuse the intercept path so the caller can
/// use its exact legacy fallback.
fn nsga3_solve_hyperplane(matrix: &[Vec<f64>]) -> Option<Vec<f64>> {
    const PIVOT_RATIO_TOLERANCE: f64 = 1e-12;
    const CONDITION_ERROR_LIMIT: f64 = 1e-8;

    let n = matrix.len();
    if n == 0 || matrix.iter().any(|row| row.len() != n) {
        return None;
    }
    let matrix_scale = matrix
        .iter()
        .flatten()
        .map(|value| value.abs())
        .fold(0.0f64, f64::max);
    if !matrix_scale.is_finite() || matrix_scale == 0.0 {
        return None;
    }

    let original = matrix.to_vec();
    let mut lu = original.clone();
    let mut swaps = vec![0usize; n];
    let mut smallest_pivot = f64::INFINITY;
    let mut largest_pivot = 0.0f64;

    for column in 0..n {
        let mut pivot_row = column;
        let mut pivot_abs = lu[column][column].abs();
        for row in column + 1..n {
            let candidate = lu[row][column].abs();
            if candidate > pivot_abs {
                pivot_row = row;
                pivot_abs = candidate;
            }
        }
        if !pivot_abs.is_finite() || pivot_abs == 0.0 {
            return None;
        }
        lu.swap(column, pivot_row);
        swaps[column] = pivot_row;
        smallest_pivot = smallest_pivot.min(pivot_abs);
        largest_pivot = largest_pivot.max(pivot_abs);

        for row in column + 1..n {
            let factor = lu[row][column] / lu[column][column];
            if !factor.is_finite() {
                return None;
            }
            lu[row][column] = factor;
            for entry in column + 1..n {
                lu[row][entry] = factor.mul_add(-lu[column][entry], lu[row][entry]);
                if !lu[row][entry].is_finite() {
                    return None;
                }
            }
        }
    }

    if smallest_pivot <= largest_pivot * PIVOT_RATIO_TOLERANCE {
        return None;
    }

    let matrix_norm_1 = (0..n)
        .map(|column| original.iter().map(|row| row[column].abs()).sum::<f64>())
        .fold(0.0f64, f64::max);
    let mut inverse_norm_1 = 0.0f64;
    for column in 0..n {
        let mut basis = vec![0.0; n];
        basis[column] = 1.0;
        let inverse_column = nsga3_lu_solve(&lu, &swaps, &basis)?;
        inverse_norm_1 = inverse_norm_1.max(inverse_column.iter().map(|v| v.abs()).sum());
    }
    let condition_1 = matrix_norm_1 * inverse_norm_1;
    if !condition_1.is_finite() || condition_1 * f64::EPSILON * n as f64 > CONDITION_ERROR_LIMIT {
        return None;
    }

    let solution = nsga3_lu_solve(&lu, &swaps, &vec![1.0; n])?;
    let residual_tolerance = 512.0 * f64::EPSILON * n as f64;
    for row in &original {
        let projected: f64 = row
            .iter()
            .zip(&solution)
            .map(|(coefficient, value)| coefficient * value)
            .sum();
        let magnitude: f64 = row
            .iter()
            .zip(&solution)
            .map(|(coefficient, value)| (coefficient * value).abs())
            .sum();
        if !projected.is_finite()
            || !magnitude.is_finite()
            || (projected - 1.0).abs() > residual_tolerance * (1.0 + magnitude)
        {
            return None;
        }
    }
    Some(solution)
}

fn nsga3_lu_solve(lu: &[Vec<f64>], swaps: &[usize], rhs: &[f64]) -> Option<Vec<f64>> {
    let n = lu.len();
    if rhs.len() != n || swaps.len() != n {
        return None;
    }
    let mut solution = rhs.to_vec();
    for (column, &pivot_row) in swaps.iter().enumerate() {
        solution.swap(column, pivot_row);
    }
    for row in 0..n {
        let known: f64 = lu[row][..row]
            .iter()
            .zip(&solution[..row])
            .map(|(coefficient, value)| coefficient * value)
            .sum();
        solution[row] -= known;
        if !solution[row].is_finite() {
            return None;
        }
    }
    for row in (0..n).rev() {
        let known: f64 = lu[row][row + 1..]
            .iter()
            .zip(&solution[row + 1..])
            .map(|(coefficient, value)| coefficient * value)
            .sum();
        solution[row] = (solution[row] - known) / lu[row][row];
        if !solution[row].is_finite() {
            return None;
        }
    }
    Some(solution)
}

fn nsga3_association(
    pop: &[Individual],
    directions: &[Vec<f64>],
    ideal: &[f64],
    span: &[f64],
    i: usize,
) -> (usize, f64) {
    let fnorm: Vec<f64> = pop[i]
        .f
        .iter()
        .zip(ideal.iter().zip(span))
        .map(|(v, (id, sp))| (v - id) / sp.max(1e-30))
        .collect();
    let mut best = (0usize, f64::INFINITY);
    for (k, dir) in directions.iter().enumerate() {
        let d = perp_distance(&fnorm, dir);
        if d < best.1 {
            best = (k, d);
        }
    }
    best
}

/// MOEA/D configuration.
#[derive(Debug, Clone, Copy)]
pub struct MoeadParams {
    /// Neighborhood size T.
    pub neighbors: usize,
    /// Max replacements per child (bounded-replacement variant —
    /// preserves diversity).
    pub max_replace: usize,
    /// Generations (each visits every subproblem once).
    pub generations: usize,
    /// SBX distribution index.
    pub eta_c: f64,
    /// Polynomial-mutation distribution index.
    pub eta_m: f64,
    /// Per-variable mutation probability.
    pub p_mut: f64,
    /// Root seed for all derived study streams.
    pub seed: u64,
}

/// Tchebycheff scalarization g(f | w, z*) = max_i w_i·|f_i − z*_i|
/// (weights floored at 1e−6 so no objective is invisible).
fn tchebycheff(f: &[f64], w: &[f64], ideal: &[f64]) -> f64 {
    f.iter()
        .zip(w.iter().zip(ideal))
        .map(|(fi, (wi, zi))| wi.max(1e-6) * (fi - zi).abs())
        .fold(f64::NEG_INFINITY, f64::max)
}

/// MOEA/D with Tchebycheff decomposition: one subproblem per weight
/// vector (use [`das_dennis`]), T-nearest-weight neighborhoods,
/// neighborhood mating, ideal-point tracking, and bounded
/// neighborhood replacement. Deterministic per seed. Returns the
/// final population's non-dominated front.
pub fn moead(
    objectives: &mut dyn FnMut(&[f64]) -> Vec<f64>,
    dim: usize,
    bounds: (f64, f64),
    weights: &[Vec<f64>],
    params: &MoeadParams,
) -> Vec<Individual> {
    let n = weights.len();
    let weight_m = direction_dimension(weights, "MOEA/D weight vectors");
    assert!(
        params.neighbors > 0,
        "MOEA/D neighborhood size must be positive"
    );
    check_variation_inputs(dim, bounds, params.eta_c, params.eta_m, params.p_mut);
    let t = params.neighbors.min(n);
    let (lo, hi) = bounds;
    // T nearest weights per subproblem (Euclidean, index tie-break).
    let mut hood: Vec<Vec<usize>> = Vec::with_capacity(n);
    for i in 0..n {
        let mut order: Vec<usize> = (0..n).collect();
        order.sort_by(|&a, &b| {
            let da: f64 = weights[i]
                .iter()
                .zip(&weights[a])
                .map(|(p, q)| (p - q) * (p - q))
                .sum();
            let db: f64 = weights[i]
                .iter()
                .zip(&weights[b])
                .map(|(p, q)| (p - q) * (p - q))
                .sum();
            da.total_cmp(&db).then(a.cmp(&b))
        });
        order.truncate(t);
        hood.push(order);
    }
    let mut stream = StreamKey {
        seed: params.seed,
        kernel: 0x0D0E,
        tile: 0,
    }
    .stream();
    let mut pop: Vec<Individual> = (0..n)
        .map(|_| {
            let x: Vec<f64> = (0..dim)
                .map(|_| (hi - lo).mul_add(stream.next_f64(), lo))
                .collect();
            let f = objectives(&x);
            Individual { x, f }
        })
        .collect();
    let m = objective_dimension(&pop, "MOEA/D");
    assert_eq!(
        weight_m, m,
        "MOEA/D weight-vector dimension must match objective dimension"
    );
    let mut ideal = vec![f64::INFINITY; m];
    for ind in &pop {
        for (z, v) in ideal.iter_mut().zip(&ind.f) {
            *z = z.min(*v);
        }
    }
    for _ in 0..params.generations {
        for neighbors in &hood {
            // Neighborhood mating.
            let a = neighbors[stream.next_below(t as u64) as usize];
            let b = neighbors[stream.next_below(t as u64) as usize];
            let (mut c1, _) = sbx(&pop[a].x, &pop[b].x, params.eta_c, lo, hi, &mut stream);
            mutate(&mut c1, params.eta_m, params.p_mut, lo, hi, &mut stream);
            let fc = objectives(&c1);
            for (z, v) in ideal.iter_mut().zip(&fc) {
                *z = z.min(*v);
            }
            // Bounded neighborhood replacement.
            let mut replaced = 0usize;
            for &j in neighbors {
                if replaced >= params.max_replace {
                    break;
                }
                if tchebycheff(&fc, &weights[j], &ideal)
                    < tchebycheff(&pop[j].f, &weights[j], &ideal)
                {
                    pop[j] = Individual {
                        x: c1.clone(),
                        f: fc.clone(),
                    };
                    replaced += 1;
                }
            }
        }
    }
    let fronts = non_dominated_sort(&pop);
    pop.into_iter()
        .zip(fronts)
        .filter(|(_, r)| *r == 0)
        .map(|(ind, _)| ind)
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn individual(objectives: [f64; 3]) -> Individual {
        Individual {
            x: Vec::new(),
            f: objectives.to_vec(),
        }
    }

    fn intercept_fixture() -> (Vec<Individual>, Vec<usize>, Vec<usize>) {
        let pop = vec![
            individual([4.0, 1.0, 2.0]),
            individual([2.0, 5.0, 1.0]),
            individual([1.0, 2.0, 6.0]),
            individual([0.0, 3.0, 3.0]),
            individual([3.0, 0.0, 3.0]),
            individual([3.0, 3.0, 0.0]),
        ];
        let fronts = vec![0; pop.len()];
        let considered = (0..pop.len()).collect();
        (pop, fronts, considered)
    }

    fn legacy_normalization(
        pop: &[Individual],
        fronts: &[usize],
        considered: &[usize],
    ) -> (Vec<f64>, Vec<f64>) {
        let m = pop[0].f.len();
        let mut ideal = vec![f64::INFINITY; m];
        for &i in considered {
            for (d, v) in ideal.iter_mut().zip(&pop[i].f) {
                *d = d.min(*v);
            }
        }
        let mut span = vec![1e-30f64; m];
        for &i in considered {
            if fronts[i] == 0 {
                for (d, (v, id)) in span.iter_mut().zip(pop[i].f.iter().zip(&ideal)) {
                    *d = d.max(v - id);
                }
            }
        }
        (ideal, span)
    }

    #[test]
    fn nsga3_normalization_uses_independent_hyperplane_intercepts() {
        let (pop, fronts, considered) = intercept_fixture();
        let (ideal, span) = nsga3_normalization(&pop, &fronts, &considered, 3);

        assert_eq!(ideal, vec![0.0, 0.0, 0.0]);
        let (_, fallback_span) = legacy_normalization(&pop, &fronts, &considered);
        assert_eq!(
            nsga3_extreme_indices(&pop, &considered, &ideal, &fallback_span, 3),
            Some(vec![0, 1, 2])
        );
        let expected = [99.0 / 17.0, 9.0, 99.0 / 10.0];
        for (value, expected) in span.iter().zip(expected) {
            assert!((value - expected).abs() <= 32.0 * f64::EPSILON * expected);
        }
        assert_ne!(
            span.iter().map(|value| value.to_bits()).collect::<Vec<_>>(),
            fallback_span
                .iter()
                .map(|value| value.to_bits())
                .collect::<Vec<_>>()
        );
    }

    #[test]
    fn nsga3_normalization_collinear_extremes_use_exact_maxima_fallback() {
        let pop = vec![
            individual([2.0, 0.0, 0.0]),
            individual([0.0, 2.0, 0.0]),
            individual([1.0, 1.0, 0.0]),
        ];
        let fronts = vec![0; pop.len()];
        let considered: Vec<usize> = (0..pop.len()).collect();
        let (ideal, span) = nsga3_normalization(&pop, &fronts, &considered, 3);
        let (legacy_ideal, legacy_span) = legacy_normalization(&pop, &fronts, &considered);

        assert_eq!(ideal, vec![0.0, 0.0, 0.0]);
        assert_eq!(
            ideal
                .iter()
                .map(|value| value.to_bits())
                .collect::<Vec<_>>(),
            legacy_ideal
                .iter()
                .map(|value| value.to_bits())
                .collect::<Vec<_>>()
        );
        assert_eq!(
            nsga3_extreme_indices(&pop, &considered, &ideal, &legacy_span, 3),
            Some(vec![0, 1, 2])
        );
        assert_eq!(
            span.iter().map(|value| value.to_bits()).collect::<Vec<_>>(),
            legacy_span
                .iter()
                .map(|value| value.to_bits())
                .collect::<Vec<_>>()
        );
    }

    #[test]
    fn nsga3_normalization_preserves_associations_under_dyadic_affine_scaling() {
        let (pop, fronts, considered) = intercept_fixture();
        let (ideal, span) = nsga3_normalization(&pop, &fronts, &considered, 3);
        let scales = [2.0, 0.5, 4.0];
        let shifts = [8.0, -3.0, 16.0];
        let transformed: Vec<Individual> = pop
            .iter()
            .map(|individual| Individual {
                x: Vec::new(),
                f: individual
                    .f
                    .iter()
                    .enumerate()
                    .map(|(objective, value)| value.mul_add(scales[objective], shifts[objective]))
                    .collect(),
            })
            .collect();
        let (transformed_ideal, transformed_span) =
            nsga3_normalization(&transformed, &fronts, &considered, 3);
        let (_, fallback_span) = legacy_normalization(&pop, &fronts, &considered);
        let (_, transformed_fallback_span) =
            legacy_normalization(&transformed, &fronts, &considered);
        let directions = das_dennis(3, 3);

        assert_ne!(
            span.iter().map(|value| value.to_bits()).collect::<Vec<_>>(),
            fallback_span
                .iter()
                .map(|value| value.to_bits())
                .collect::<Vec<_>>()
        );
        assert_ne!(
            transformed_span
                .iter()
                .map(|value| value.to_bits())
                .collect::<Vec<_>>(),
            transformed_fallback_span
                .iter()
                .map(|value| value.to_bits())
                .collect::<Vec<_>>()
        );
        assert_eq!(
            nsga3_extreme_indices(&pop, &considered, &ideal, &fallback_span, 3),
            nsga3_extreme_indices(
                &transformed,
                &considered,
                &transformed_ideal,
                &transformed_fallback_span,
                3,
            )
        );
        for objective in 0..3 {
            let expected_ideal = ideal[objective].mul_add(scales[objective], shifts[objective]);
            assert_eq!(
                transformed_ideal[objective].to_bits(),
                expected_ideal.to_bits()
            );
            let expected_span = span[objective] * scales[objective];
            assert!(
                (transformed_span[objective] - expected_span).abs()
                    <= 128.0 * f64::EPSILON * expected_span.abs().max(1.0)
            );
        }
        for i in &considered {
            for objective in 0..3 {
                let normalized = (pop[*i].f[objective] - ideal[objective]) / span[objective];
                let transformed_normalized = (transformed[*i].f[objective]
                    - transformed_ideal[objective])
                    / transformed_span[objective];
                assert!(
                    (normalized - transformed_normalized).abs()
                        <= 128.0 * f64::EPSILON * normalized.abs().max(1.0)
                );
            }
            let base_association = nsga3_association(&pop, &directions, &ideal, &span, *i);
            let transformed_association = nsga3_association(
                &transformed,
                &directions,
                &transformed_ideal,
                &transformed_span,
                *i,
            );
            assert_eq!(base_association.0, transformed_association.0);
            assert!(
                (base_association.1 - transformed_association.1).abs()
                    <= 256.0 * f64::EPSILON * base_association.1.abs().max(1.0)
            );
        }
    }

    #[test]
    fn nsga3_normalization_asf_ties_choose_lowest_population_index() {
        let (base, _, _) = intercept_fixture();
        let mut pop = vec![base[0].clone(), base[0].clone()];
        pop.extend(base.into_iter().skip(1));
        let fronts = vec![0; pop.len()];
        let mut considered: Vec<usize> = (1..pop.len()).collect();
        considered.push(0);
        let (ideal, _) = nsga3_normalization(&pop, &fronts, &considered, 3);
        let (_, fallback_span) = legacy_normalization(&pop, &fronts, &considered);

        assert_eq!(
            nsga3_extreme_indices(&pop, &considered, &ideal, &fallback_span, 3),
            Some(vec![0, 2, 3])
        );
    }

    #[test]
    fn nsga3_normalization_asf_scans_all_considered_points() {
        let pop = vec![
            individual([1.0, 1.0, 0.0]),
            individual([0.0, 1.0, 0.0]),
            individual([10.0, 0.0, 10.0]),
        ];
        let fronts = vec![1, 0, 0];
        let considered = vec![1, 2, 0];
        let (ideal, fallback_span) = legacy_normalization(&pop, &fronts, &considered);

        assert_eq!(
            nsga3_asf_extreme_index(&pop, &considered, &ideal, &fallback_span, 0, 1e-6,),
            Some(0)
        );
    }

    #[test]
    fn nsga3_normalization_duplicate_extremes_use_exact_maxima_fallback() {
        let pop = vec![
            individual([0.0, 0.0, 1.0]),
            individual([1.0, 1.0, 0.0]),
            individual([2.0, 2.0, 2.0]),
        ];
        let fronts = vec![0, 0, 1];
        let considered: Vec<usize> = (0..pop.len()).collect();
        let (legacy_ideal, legacy_span) = legacy_normalization(&pop, &fronts, &considered);
        let (ideal, span) = nsga3_normalization(&pop, &fronts, &considered, 3);

        assert_eq!(
            nsga3_extreme_indices(&pop, &considered, &ideal, &legacy_span, 3),
            None
        );
        assert_eq!(
            ideal
                .iter()
                .map(|value| value.to_bits())
                .collect::<Vec<_>>(),
            legacy_ideal
                .iter()
                .map(|value| value.to_bits())
                .collect::<Vec<_>>()
        );
        assert_eq!(
            span.iter().map(|value| value.to_bits()).collect::<Vec<_>>(),
            legacy_span
                .iter()
                .map(|value| value.to_bits())
                .collect::<Vec<_>>()
        );
    }

    #[test]
    fn nsga3_normalization_tiny_pivot_hyperplane_is_refused() {
        let delta = 2.0f64.powi(-48);
        let matrix = vec![
            vec![1.0, 0.0, 0.0],
            vec![0.0, 1.0, 0.0],
            vec![1.0, 1.0, delta],
        ];
        assert_eq!(nsga3_solve_hyperplane(&matrix), None);
    }

    #[test]
    fn nsga3_normalization_condition_gate_refuses_nonzero_pivots() {
        let delta = 2.0f64.powi(-26);
        let matrix = vec![
            vec![1.0, 0.0, 0.0],
            vec![0.0, 1.0, 0.0],
            vec![1.0, 1.0, delta],
        ];
        assert_eq!(nsga3_solve_hyperplane(&matrix), None);
    }

    #[test]
    fn nsga3_normalization_lu_applies_row_permutation_to_nonuniform_rhs() {
        let lu = vec![vec![1.0, 3.0], vec![0.0, 2.0]];
        let solution = nsga3_lu_solve(&lu, &[1, 1], &[5.0, 7.0]).expect("LU is valid");
        assert_eq!(
            solution
                .iter()
                .map(|value| value.to_bits())
                .collect::<Vec<_>>(),
            [-0.5f64, 2.5]
                .iter()
                .map(|value| value.to_bits())
                .collect::<Vec<_>>()
        );
    }

    #[test]
    fn nsga3_normalization_large_objective_count_uses_exact_bounded_fallback() {
        let m = 65;
        let pop: Vec<Individual> = (0..m)
            .map(|axis| {
                let mut f = vec![1.0; m];
                f[axis] = 0.0;
                Individual { x: Vec::new(), f }
            })
            .collect();
        let fronts = vec![0; pop.len()];
        let considered: Vec<usize> = (0..pop.len()).collect();
        let (legacy_ideal, legacy_span) = legacy_normalization(&pop, &fronts, &considered);
        let (ideal, span) = nsga3_normalization(&pop, &fronts, &considered, m);

        assert_eq!(
            nsga3_hyperplane_span(&pop, &considered, &legacy_ideal, &legacy_span),
            None
        );
        assert_eq!(
            ideal
                .iter()
                .map(|value| value.to_bits())
                .collect::<Vec<_>>(),
            legacy_ideal
                .iter()
                .map(|value| value.to_bits())
                .collect::<Vec<_>>()
        );
        assert_eq!(
            span.iter().map(|value| value.to_bits()).collect::<Vec<_>>(),
            legacy_span
                .iter()
                .map(|value| value.to_bits())
                .collect::<Vec<_>>()
        );
    }

    #[test]
    fn nsga3_normalization_zero_and_near_zero_spans_stay_finite() {
        let tiny = f64::from_bits(1);
        let pop = vec![
            individual([2.0, 0.0, tiny]),
            individual([0.0, 2.0, tiny]),
            individual([1.0, 1.0, f64::from_bits(2)]),
        ];
        let fronts = vec![0; pop.len()];
        let considered: Vec<usize> = (0..pop.len()).collect();
        let (ideal, span) = nsga3_normalization(&pop, &fronts, &considered, 3);
        let directions = das_dennis(3, 2);

        assert!(ideal.iter().all(|value| value.is_finite()));
        assert!(span.iter().all(|value| value.is_finite() && *value > 0.0));
        for i in considered {
            let (_, distance) = nsga3_association(&pop, &directions, &ideal, &span, i);
            assert!(distance.is_finite());
        }
    }
}
