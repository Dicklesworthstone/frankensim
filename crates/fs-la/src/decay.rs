//! Prepared positive evolution of a finite acyclic first-order transfer network.
//!
//! A parent with rate lambda leaves its population once and enters exactly one
//! daughter, selected by the declared branch fractions. Stable nodes have zero
//! rate and no outgoing edges. This is a numerical owner, not a nuclide catalog:
//! identities, physical rate data, emitted products and energy belong to callers.
//!
//! For column populations, N' = A N, where A_ii = -lambda_i and
//! A_daughter,parent += lambda_parent * fraction. Preparation retains BOTH
//! P = exp(A t) and C = integral diag(lambda) exp(A s) ds. Thus C N0
//! counts each parent's transformations, including transformations of daughters
//! created during the interval. Counts are never recovered by subtracting two
//! nearly equal inventories. P and C use the caller's original node order.
//!
//! Uniformization at max(lambda)*h <= 1/2 uses only nonnegative terms, followed
//! by triangular semigroup composition P(2h)=P(h)P(h) and
//! C(2h)=C(h)+C(h)P(h). Exact scalar exp/expm1 diagonals are reinserted at
//! EVERY scale: a slow parent's tiny loss must not disappear because another
//! parent required a much smaller initial interval. Equal rates need no special
//! case or Bateman divided difference. Work grows with the logarithm of the
//! fastest dimensionless rate, not with the number of fastest lifetimes.
//!
//! The positive series has a componentwise tail test: after term k >= n, a
//! path has at most n non-self transitions in the augmented population/count
//! system, so subsequent term ratios are bounded by mu/(k+1-n). This bounds
//! truncation in exact arithmetic, not floating-point error. Final stochastic,
//! per-node balance and visit-count checks are numerical admission tests; no
//! interval enclosure, adjoint, branching-data accuracy or physical validation
//! is claimed. Unrepresentable direct base transfers refuse instead of silently
//! dropping a declared positive branch. Unrepresentable late exponential tails
//! and smaller multi-hop series terms may underflow to zero. Application refuses
//! inventories that could magnify an unrepresented coefficient into a resolved
//! contribution; it never treats a large-population tail as necessarily zero.

use fs_math::det;

/// Hard limit for the bounded dense prepared maps.
pub const MAX_DECAY_NODES: usize = 64;
/// Hard limit for explicit branches, including repeated parent/daughter pairs.
pub const MAX_DECAY_EDGES: usize = MAX_DECAY_NODES * MAX_DECAY_NODES;
const MAX_SQUARINGS: usize = 2048;
const MAX_SERIES_TERMS: usize = 256;

/// One explicitly supplied branch. Several branches may share a daughter.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DecayEdge {
    pub parent: usize,
    pub daughter: usize,
    /// Finite positive fraction; all fractions of a decaying parent sum to one.
    pub branch_fraction: f64,
}

/// Numerical resources and admission accuracy for one prepared interval.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DecayLimits {
    pub maximum_nodes: usize,
    pub maximum_edges: usize,
    pub maximum_squarings: usize,
    pub maximum_series_terms: usize,
    /// Conservative scalar-term work, charged before each bounded work block.
    pub maximum_work: u64,
    /// Componentwise tail target and final relative balance admission. This is
    /// not a certified forward error for every entry of the exponential.
    pub relative_tolerance: f64,
}

impl Default for DecayLimits {
    fn default() -> Self {
        Self {
            maximum_nodes: MAX_DECAY_NODES,
            maximum_edges: MAX_DECAY_EDGES,
            maximum_squarings: MAX_SQUARINGS,
            maximum_series_terms: 128,
            maximum_work: 1_000_000_000,
            relative_tolerance: 2.0e-12,
        }
    }
}

/// Preparation/application refusal; no caller-owned population is mutated.
#[derive(Clone, Debug, PartialEq)]
pub enum DecayError {
    Invalid(&'static str),
    Budget { quantity: &'static str, required: u64, maximum: u64 },
    Cancelled,
    Unrepresentable(&'static str),
    SeriesNotConverged { terms: usize },
    Balance { row: usize, column: usize, residual: f64, allowance: f64 },
}
impl core::fmt::Display for DecayError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "first-order decay refused: {self:?}")
    }
}
impl std::error::Error for DecayError {}

/// Original-order endpoint populations and independently integrated counts.
#[derive(Clone, Debug, PartialEq)]
pub struct DecayAdvance {
    pub populations: Vec<f64>,
    /// decays[i] = integral lambda_i N_i(s) ds; stable entries are zero.
    pub decays: Vec<f64>,
    /// Signed final minus initial total population, before any correction.
    pub population_balance_residual: f64,
}

/// Immutable interval propagator, reusable for many initial inventories.
#[derive(Clone, Debug)]
pub struct DecayPlan {
    nodes: usize,
    duration_s: f64,
    transition: Vec<f64>,
    integrated: Vec<f64>,
    squarings: usize,
    series_terms: usize,
    work: u64,
    tolerance: f64,
    // Original-order graph bounds distinguish true structural zeros from
    // coefficients lost to underflow. Empty only for an identity interval.
    tails: Vec<TailBound>,
    rates: Vec<f64>,
}

#[derive(Clone, Copy)]
struct Incoming {
    parent: usize,
    fraction: f64,
    probability: f64,
}

#[derive(Clone, Copy, Debug)]
struct TailBound {
    reachable: bool,
    minimum_rate: f64,
    maximum_steps: usize,
}
impl TailBound {
    const NONE: Self = Self { reachable: false, minimum_rate: f64::INFINITY, maximum_steps: 0 };

    // On every path from the source through this transient target, at most
    // L+1 waits have rates >= minimum_rate. Remaining in the target is bounded
    // by the survival of Erlang(L+1,minimum_rate), independently of branching.
    // Only a large margin below the smallest representable positive value
    // permits a vanished coefficient to be used with an inventory above one.
    fn scaled_tail_vanishes(self, population: f64, duration: f64) -> bool {
        if self.minimum_rate == 0.0 { return false; }
        let exposure = self.minimum_rate * duration;
        if exposure.is_infinite() { return true; }
        if exposure <= 0.0 || exposure < self.maximum_steps as f64 { return false; }
        let mut backward_sum = 1.0;
        let mut term = 1.0;
        let mut log_factorial = 0.0;
        for k in (1..=self.maximum_steps).rev() {
            term *= k as f64 / exposure;
            backward_sum += term;
            log_factorial += det::ln(k as f64);
        }
        let log_bound = det::ln(population) - exposure
            + self.maximum_steps as f64 * det::ln(exposure)
            - log_factorial + det::ln(backward_sum);
        // ln(minimum_subnormal/2) is approximately -745.133. The additional
        // margin dominates rounding in these <=64-term logarithmic bounds.
        log_bound < -746.0
    }
}

#[derive(Clone, Copy, Default)]
struct Sum {
    value: f64,
    correction: f64,
}
impl Sum {
    fn add(&mut self, value: f64) {
        let next = self.value + value;
        self.correction += if self.value.abs() >= value.abs() {
            (self.value - next) + value
        } else {
            (value - next) + self.value
        };
        self.value = next;
    }
    fn product(&mut self, a: f64, b: f64) {
        let product = a * b;
        self.add(product);
        self.add(a.mul_add(b, -product));
    }
    fn total(self) -> f64 { self.value + self.correction }
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), DecayError> {
    if cancelled() { Err(DecayError::Cancelled) } else { Ok(()) }
}
fn charge(used: &mut u64, amount: usize, maximum: u64) -> Result<(), DecayError> {
    let required = used.checked_add(amount as u64)
        .ok_or(DecayError::Unrepresentable("preparation work counter"))?;
    if required > maximum {
        return Err(DecayError::Budget { quantity: "preparation work", required, maximum });
    }
    *used = required;
    Ok(())
}
fn cap(quantity: &'static str, required: usize, maximum: usize) -> Result<(), DecayError> {
    if required > maximum {
        Err(DecayError::Budget { quantity, required: required as u64, maximum: maximum as u64 })
    } else { Ok(()) }
}
fn positive_or_zero(value: f64) -> bool { value.is_finite() && value >= 0.0 }

fn diagonal(rate: f64, time: f64) -> (f64, f64) {
    if rate == 0.0 || time == 0.0 { return (1.0, 0.0); }
    let exposure = rate * time;
    (det::exp(-exposure), -det::expm1(-exposure))
}

impl DecayPlan {
    /// Prepare an unchanged-rate, closed, acyclic branching network. Fractions
    /// must sum to one within 32 machine epsilons; they are never renormalized.
    /// Every positive-rate node needs outgoing branches to represented nodes.
    /// Node order is unrestricted; the stable topological order is internal.
    /// Zero duration and completely stable networks produce the identity.
    /// cancelled is polled between rows and scaling/series stages.
    pub fn prepare(
        rates_per_s: &[f64],
        edges: &[DecayEdge],
        duration_s: f64,
        limits: DecayLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, DecayError> {
        poll(&mut cancelled)?;
        let n = rates_per_s.len();
        if limits.maximum_nodes == 0 || limits.maximum_nodes > MAX_DECAY_NODES
            || limits.maximum_edges > MAX_DECAY_EDGES
            || limits.maximum_squarings > MAX_SQUARINGS
            || limits.maximum_series_terms == 0 || limits.maximum_series_terms > MAX_SERIES_TERMS
            || limits.maximum_work == 0
            || !limits.relative_tolerance.is_finite()
            || limits.relative_tolerance <= 0.0 || limits.relative_tolerance >= 1.0
        {
            return Err(DecayError::Invalid("finite bounded decay limits required"));
        }
        cap("nodes", n, limits.maximum_nodes)?;
        cap("edges", edges.len(), limits.maximum_edges)?;
        if n == 0 || !positive_or_zero(duration_s)
            || rates_per_s.iter().any(|&rate| !positive_or_zero(rate))
        {
            return Err(DecayError::Invalid("nonempty finite nonnegative rates and duration required"));
        }
        let mut work = 0;
        charge(&mut work, n + edges.len(), limits.maximum_work)?;
        let mut indegree = vec![0usize; n];
        let mut fractions = vec![Sum::default(); n];
        for edge in edges {
            if edge.parent >= n || edge.daughter >= n || edge.parent == edge.daughter
                || !edge.branch_fraction.is_finite() || edge.branch_fraction <= 0.0
                || edge.branch_fraction > 1.0 || rates_per_s[edge.parent] == 0.0
            {
                return Err(DecayError::Invalid("each branch needs distinct existing nodes, a decaying parent and a positive fraction"));
            }
            indegree[edge.daughter] += 1;
            fractions[edge.parent].add(edge.branch_fraction);
        }
        for (i, &rate) in rates_per_s.iter().enumerate() {
            let sum = fractions[i].total();
            if rate > 0.0 && (sum - 1.0).abs() > 32.0 * f64::EPSILON {
                return Err(DecayError::Invalid("all decaying-parent branch fractions must sum to one"));
            }
        }
        let mut order = Vec::with_capacity(n);
        let mut selected = vec![false; n];
        for _ in 0..n {
            poll(&mut cancelled)?;
            charge(&mut work, n + edges.len(), limits.maximum_work)?;
            let next = (0..n).find(|&i| !selected[i] && indegree[i] == 0)
                .ok_or(DecayError::Invalid("decay network must be acyclic"))?;
            selected[next] = true;
            order.push(next);
            for edge in edges.iter().filter(|edge| edge.parent == next) {
                indegree[edge.daughter] -= 1;
            }
        }
        let mut position = vec![0usize; n];
        for (i, &original) in order.iter().enumerate() { position[original] = i; }
        let rates: Vec<_> = order.iter().map(|&i| rates_per_s[i]).collect();
        let maximum_rate = rates.iter().copied().fold(0.0_f64, f64::max);
        let mut transition = vec![0.0; n * n];
        let mut integrated = vec![0.0; n * n];
        for i in 0..n { transition[i * n + i] = 1.0; }
        if duration_s == 0.0 || maximum_rate == 0.0 {
            poll(&mut cancelled)?;
            return Ok(Self { nodes: n, duration_s, transition, integrated, squarings: 0,
                series_terms: 0, work, tolerance: limits.relative_tolerance,
                tails: Vec::new(), rates: rates_per_s.to_vec() });
        }
        let mut h = duration_s;
        let mut squarings = 0;
        while maximum_rate * h > 0.5 {
            poll(&mut cancelled)?;
            squarings += 1;
            cap("squarings", squarings, limits.maximum_squarings)?;
            let next = 0.5 * h;
            if next == 0.0 || next * 2.0 != h {
                return Err(DecayError::Unrepresentable("exact scaled duration"));
            }
            h = next;
        }
        let mu = maximum_rate * h;
        if !mu.is_finite() || mu <= 0.0 {
            return Err(DecayError::Unrepresentable("positive scaled rate interval"));
        }
        let ratios: Vec<_> = rates.iter().map(|&rate| rate / maximum_rate).collect();
        for (i, &rate) in rates.iter().enumerate() {
            if rate > 0.0 && (ratios[i] == 0.0 || rate * h == 0.0) {
                return Err(DecayError::Unrepresentable("rate hierarchy at the scaled interval"));
            }
        }
        let mut incoming = vec![Vec::new(); n];
        for edge in edges {
            let parent = position[edge.parent];
            let daughter = position[edge.daughter];
            let probability = ratios[parent] * edge.branch_fraction;
            if probability == 0.0 || (rates[parent] * h) * edge.branch_fraction == 0.0 {
                return Err(DecayError::Unrepresentable("positive direct branch at the scaled interval"));
            }
            incoming[daughter].push(Incoming { parent, fraction: edge.branch_fraction, probability });
        }
        // Original input order remains the tie-break among identical endpoints.
        for row in &mut incoming { row.sort_by_key(|edge| edge.parent); }
        let mut tails = vec![TailBound::NONE; n * n];
        for i in 0..n {
            poll(&mut cancelled)?;
            charge(&mut work, 4 * n * (incoming[i].len() + 1), limits.maximum_work)?;
            tails[i * n + i] = TailBound { reachable: true, minimum_rate: rates[i], maximum_steps: 0 };
            for edge in &incoming[i] {
                for j in 0..=edge.parent {
                    let previous = tails[edge.parent * n + j];
                    if previous.reachable {
                        let bound = &mut tails[i * n + j];
                        bound.reachable = true;
                        bound.minimum_rate = bound.minimum_rate.min(previous.minimum_rate.min(rates[i]));
                        bound.maximum_steps = bound.maximum_steps.max(previous.maximum_steps + 1);
                    }
                }
            }
        }
        let mut pterm = vec![0.0; n * n];
        let mut cterm = vec![0.0; n * n];
        let mut psum = vec![Sum::default(); n * n];
        let mut csum = vec![Sum::default(); n * n];
        for i in 0..n {
            pterm[i * n + i] = det::exp(-mu);
            psum[i * n + i].add(pterm[i * n + i]);
        }
        let mut pnext = vec![0.0; n * n];
        let mut cnext = vec![0.0; n * n];
        // At most n genuine transitions multiply any individual path. Reserve
        // most of the tolerance for floating arithmetic and semigroup staging.
        let series_tolerance = limits.relative_tolerance.min(1.0e-15)
            / (8.0 * n as f64 * (squarings + 1) as f64);
        let mut series_terms = 0;
        for k in 1..=limits.maximum_series_terms {
            poll(&mut cancelled)?;
            let weight = mu / k as f64;
            for i in 0..n {
                poll(&mut cancelled)?;
                charge(&mut work, 8 * (i + 1) * (incoming[i].len() + 3), limits.maximum_work)?;
                for j in 0..=i {
                    let at = i * n + j;
                    let mut p = Sum::default();
                    p.product(1.0 - ratios[i], pterm[at]);
                    for edge in &incoming[i] {
                        p.product(edge.probability, pterm[edge.parent * n + j]);
                    }
                    pnext[at] = weight * p.total();
                    let mut c = Sum::default();
                    c.add(cterm[at]);
                    c.product(ratios[i], pterm[at]);
                    cnext[at] = weight * c.total();
                    psum[at].add(pnext[at]);
                    csum[at].add(cnext[at]);
                }
            }
            core::mem::swap(&mut pterm, &mut pnext);
            core::mem::swap(&mut cterm, &mut cnext);
            if k >= n {
                let ratio = mu / (k + 1 - n) as f64;
                let tail = ratio / (1.0 - ratio);
                let converged = (0..n * n).all(|at| {
                    tail * pterm[at] <= series_tolerance * psum[at].total()
                        && tail * cterm[at] <= series_tolerance * csum[at].total()
                });
                if converged { series_terms = k; break; }
            }
        }
        if series_terms == 0 {
            return Err(DecayError::SeriesNotConverged { terms: limits.maximum_series_terms });
        }
        for at in 0..n * n {
            transition[at] = psum[at].total();
            integrated[at] = csum[at].total();
        }
        for i in 0..n {
            (transition[i * n + i], integrated[i * n + i]) = diagonal(rates[i], h);
        }
        pnext.fill(0.0);
        cnext.fill(0.0);
        for _ in 0..squarings {
            poll(&mut cancelled)?;
            h *= 2.0;
            for i in 0..n {
                poll(&mut cancelled)?;
                charge(&mut work, 16 * (i + 1) * (i + 2), limits.maximum_work)?;
                for j in 0..i {
                    let mut p = Sum::default();
                    let mut c = Sum::default();
                    c.add(integrated[i * n + j]);
                    for k in j..=i {
                        p.product(transition[i * n + k], transition[k * n + j]);
                        c.product(integrated[i * n + k], transition[k * n + j]);
                    }
                    pnext[i * n + j] = p.total();
                    cnext[i * n + j] = c.total();
                }
                (pnext[i * n + i], cnext[i * n + i]) = diagonal(rates[i], h);
            }
            core::mem::swap(&mut transition, &mut pnext);
            core::mem::swap(&mut integrated, &mut cnext);
        }
        validate_maps(&transition, &integrated, &incoming, n, limits.relative_tolerance,
            &mut work, limits.maximum_work, &mut cancelled)?;
        let mut original_tails = vec![TailBound::NONE; n * n];
        for i in 0..n { for j in 0..n {
            pnext[order[i] * n + order[j]] = transition[i * n + j];
            cnext[order[i] * n + order[j]] = integrated[i * n + j];
            original_tails[order[i] * n + order[j]] = tails[i * n + j];
        }}
        poll(&mut cancelled)?;
        Ok(Self { nodes: n, duration_s, transition: pnext, integrated: cnext, squarings,
            series_terms, work, tolerance: limits.relative_tolerance,
            tails: original_tails, rates: rates_per_s.to_vec() })
    }

    #[must_use] pub const fn node_count(&self) -> usize { self.nodes }
    #[must_use] pub const fn duration_s(&self) -> f64 { self.duration_s }
    /// Row-major original-order population propagator.
    #[must_use] pub fn transition_matrix(&self) -> &[f64] { &self.transition }
    /// Row-major original-order map from initial populations to parent counts.
    #[must_use] pub fn integrated_decay_matrix(&self) -> &[f64] { &self.integrated }
    #[must_use] pub const fn squarings(&self) -> usize { self.squarings }
    #[must_use] pub const fn series_terms(&self) -> usize { self.series_terms }
    #[must_use] pub const fn preparation_work(&self) -> u64 { self.work }

    /// Apply both prepared maps without changing the supplied inventory.
    /// Same-node counts remain available even when its survivor rounds to its
    /// initial population. Individual populations and their total must fit f64.
    ///
    /// An inventory above one must not magnify a subnormal or vanished map
    /// coefficient into a resolved contribution. Such an application refuses
    /// unless a zero transient coefficient has an independent Erlang upper
    /// bound below representability even after multiplication by that inventory.
    /// This is deliberately conservative; no logarithmic inventory is inferred.
    pub fn apply(&self, populations: &[f64]) -> Result<DecayAdvance, DecayError> {
        self.apply_with_cancel(populations, || false)
    }

    /// As apply, with bounded cooperative row checkpoints.
    pub fn apply_with_cancel(&self, populations: &[f64], mut cancelled: impl FnMut() -> bool)
        -> Result<DecayAdvance, DecayError>
    {
        poll(&mut cancelled)?;
        let n = self.nodes;
        if populations.len() != n || populations.iter().any(|&value| !positive_or_zero(value)) {
            return Err(DecayError::Invalid("one finite nonnegative population per original node required"));
        }
        let mut initial = Sum::default();
        for &population in populations { initial.add(population); }
        let initial = initial.total();
        if !initial.is_finite() { return Err(DecayError::Unrepresentable("total initial population")); }
        let mut result = vec![0.0; n];
        let mut decays = vec![0.0; n];
        for i in 0..n {
            poll(&mut cancelled)?;
            let mut p = Sum::default();
            let mut c = Sum::default();
            for (j, &population) in populations.iter().enumerate() {
                let at = i * n + j;
                if population > 1.0 && !self.tails.is_empty() && self.tails[at].reachable {
                    if self.transition[at].is_subnormal() || self.integrated[at].is_subnormal() {
                        return Err(DecayError::Unrepresentable("inventory amplifies a subnormal prepared coefficient"));
                    }
                    if self.transition[at] == 0.0
                        && !self.tails[at].scaled_tail_vanishes(population, self.duration_s)
                    {
                        return Err(DecayError::Unrepresentable("inventory may recover a vanished exponential tail"));
                    }
                    if self.integrated[at] == 0.0 && self.rates[i] > 0.0 {
                        return Err(DecayError::Unrepresentable("inventory may recover a vanished integrated count"));
                    }
                }
                p.product(self.transition[at], population);
                c.product(self.integrated[at], population);
            }
            result[i] = p.total();
            decays[i] = c.total();
            if !positive_or_zero(result[i]) || !positive_or_zero(decays[i]) {
                return Err(DecayError::Unrepresentable("endpoint population or integrated count"));
            }
        }
        let mut final_total = Sum::default();
        for &population in &result { final_total.add(population); }
        let residual = final_total.total() - initial;
        let allowance = self.tolerance * initial;
        if !residual.is_finite() || residual.abs() > allowance {
            return Err(DecayError::Balance { row: n, column: n, residual, allowance });
        }
        poll(&mut cancelled)?;
        Ok(DecayAdvance { populations: result, decays, population_balance_residual: residual })
    }
}

#[allow(clippy::too_many_arguments)]
fn validate_maps(
    p: &[f64], c: &[f64], incoming: &[Vec<Incoming>], n: usize, tolerance: f64,
    work: &mut u64, maximum_work: u64, cancelled: &mut impl FnMut() -> bool,
) -> Result<(), DecayError> {
    for j in 0..n {
        poll(cancelled)?;
        charge(work, 8 * n * (incoming.iter().map(Vec::len).sum::<usize>() + n), maximum_work)?;
        let mut column = Sum::default();
        for i in 0..n {
            let at = i * n + j;
            if !positive_or_zero(p[at]) || !positive_or_zero(c[at]) {
                return Err(DecayError::Unrepresentable("positive prepared population/count map"));
            }
            column.add(p[at]);
            if c[at] > 1.0 + tolerance {
                return Err(DecayError::Balance { row: i, column: j,
                    residual: c[at] - 1.0, allowance: tolerance });
            }
            let mut generated = Sum::default();
            for edge in &incoming[i] { generated.product(edge.fraction, c[edge.parent * n + j]); }
            let generated = generated.total();
            let initial = if i == j { 1.0 } else { 0.0 };
            let mut balance = Sum::default();
            balance.add(p[at]); balance.add(c[at]); balance.add(-initial); balance.add(-generated);
            let residual = balance.total();
            let allowance = tolerance * (p[at] + c[at] + initial + generated);
            if residual.abs() > allowance {
                return Err(DecayError::Balance { row: i, column: j, residual, allowance });
            }
        }
        let residual = column.total() - 1.0;
        if residual.abs() > tolerance {
            return Err(DecayError::Balance { row: n, column: j, residual, allowance: tolerance });
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests;
