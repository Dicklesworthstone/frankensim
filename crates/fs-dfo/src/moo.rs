//! Multi-objective machinery (plan §9.9): NSGA-II with DETERMINISTIC
//! tie-breaking everywhere (index order — bitwise-replayable runs),
//! exact hypervolume for m ≤ 4 (2D sweep + recursive exclusive
//! contributions above), knee-point detection, and CVaR through the
//! Rockafellar–Uryasev reformulation. All randomness flows through
//! fs-rand streams.

use fs_rand::StreamKey;

/// One evaluated individual.
#[derive(Debug, Clone)]
pub struct Individual {
    /// Decision vector.
    pub x: Vec<f64>,
    /// Objective vector (minimization).
    pub f: Vec<f64>,
}

/// `a` Pareto-dominates `b` (minimization).
#[must_use]
pub fn dominates(a: &[f64], b: &[f64]) -> bool {
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
    /// Master seed.
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

/// CVaR_β of a loss SAMPLE SET via the Rockafellar–Uryasev
/// reformulation: CVaR = min_α α + E[max(0, L − α)]/(1 − β). The
/// minimizer α* is the β-quantile; on samples the minimum is attained
/// at an order statistic, so we evaluate all candidates exactly
/// (deterministic; no inner optimizer needed).
#[must_use]
pub fn cvar_rockafellar_uryasev(losses: &[f64], beta: f64) -> (f64, f64) {
    assert!(beta.is_finite() && beta > 0.0 && beta < 1.0);
    assert!(!losses.is_empty(), "CVaR needs at least one loss sample");
    assert!(
        losses.iter().all(|l| l.is_finite()),
        "CVaR loss samples must be finite"
    );
    let n = losses.len() as f64;
    let mut sorted = losses.to_vec();
    sorted.sort_by(f64::total_cmp);
    let mut best = (f64::INFINITY, 0.0f64);
    let mut suffix_sum = vec![0.0f64; sorted.len() + 1];
    for i in (0..sorted.len()).rev() {
        suffix_sum[i] = suffix_sum[i + 1] + sorted[i];
    }
    for (i, &alpha) in sorted.iter().enumerate() {
        let tail_count = sorted.len() - i - 1;
        let excess = (suffix_sum[i + 1] - alpha * tail_count as f64).max(0.0);
        let val = alpha + excess / (n * (1.0 - beta));
        if val < best.0 {
            best = (val, alpha);
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
    let pts: Vec<&Vec<f64>> = front
        .iter()
        .filter(|p| p.iter().zip(reference).all(|(a, r)| a < r))
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
        sobol.point(
            u32::try_from(s + 1).expect("sample count fits u32"),
            &mut pt,
        );
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
/// niching replaces crowding distance. Normalization uses the ideal
/// point and the FIRST-front per-objective maxima as the nadir
/// estimate (the full ASF extreme-point construction is the recorded
/// refinement). Deterministic throughout (index tie-breaks).
pub fn nsga3(
    objectives: &mut dyn FnMut(&[f64]) -> Vec<f64>,
    dim: usize,
    bounds: (f64, f64),
    directions: &[Vec<f64>],
    params: &NsgaParams,
) -> Vec<Individual> {
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
            let (mut c1, mut c2) =
                sbx(&pop[p1].x, &pop[p2].x, params.eta_c, lo, hi, &mut stream);
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
        pop = nsga3_select(pop, directions, params.pop);
    }
    let fronts = non_dominated_sort(&pop);
    pop.into_iter()
        .zip(fronts)
        .filter(|(_, r)| *r == 0)
        .map(|(ind, _)| ind)
        .collect()
}

/// NSGA-III environmental selection to `target` members.
fn nsga3_select(pop: Vec<Individual>, directions: &[Vec<f64>], target: usize) -> Vec<Individual> {
    let m = pop[0].f.len();
    let fronts = non_dominated_sort(&pop);
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
    // Normalize over accepted ∪ partial: ideal point + first-front max.
    let considered: Vec<usize> = accepted.iter().chain(&partial).copied().collect();
    let mut ideal = vec![f64::INFINITY; m];
    for &i in &considered {
        for (d, v) in ideal.iter_mut().zip(&pop[i].f) {
            *d = d.min(*v);
        }
    }
    let mut span = vec![1e-30f64; m];
    for &i in &considered {
        if fronts[i] == 0 {
            for (d, (v, id)) in span.iter_mut().zip(pop[i].f.iter().zip(&ideal)) {
                *d = d.max(v - id);
            }
        }
    }
    let normalize = |i: usize| -> Vec<f64> {
        pop[i]
            .f
            .iter()
            .zip(ideal.iter().zip(&span))
            .map(|(v, (id, sp))| (v - id) / sp.max(1e-30))
            .collect()
    };
    // Associate everyone to the nearest direction.
    let associate = |i: usize| -> (usize, f64) {
        let fnorm = normalize(i);
        let mut best = (0usize, f64::INFINITY);
        for (k, dir) in directions.iter().enumerate() {
            let d = perp_distance(&fnorm, dir);
            if d < best.1 {
                best = (k, d);
            }
        }
        best
    };
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
                d < chosen_d || (d == chosen_d && i < chosen)
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
