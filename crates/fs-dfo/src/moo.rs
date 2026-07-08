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
    let pts: Vec<Vec<f64>> = front
        .iter()
        .filter(|p| p.iter().zip(reference).all(|(a, r)| a < r))
        .cloned()
        .collect();
    if pts.is_empty() {
        return 0.0;
    }
    hv_recursive(&pts, reference)
}

fn hv_recursive(pts: &[Vec<f64>], reference: &[f64]) -> f64 {
    let m = reference.len();
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
    assert!(beta > 0.0 && beta < 1.0);
    let n = losses.len() as f64;
    let mut sorted = losses.to_vec();
    sorted.sort_by(f64::total_cmp);
    let mut best = (f64::INFINITY, 0.0f64);
    for &alpha in &sorted {
        let excess: f64 = sorted
            .iter()
            .rev()
            .take_while(|&&l| l > alpha)
            .map(|l| l - alpha)
            .sum();
        let val = alpha + excess / (n * (1.0 - beta));
        if val < best.0 {
            best = (val, alpha);
        }
    }
    best
}
