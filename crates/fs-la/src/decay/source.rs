//! Constant nonnegative production integrated by the original decay owner.
//!
//! Alongside P(t) and C(t), preparation retains Fbar = integral(P)/t and
//! Gbar = integral(C)/t. Applied to J = source_rate*t, these return surviving
//! injected populations and their transformations, including daughters born
//! and transformed during the interval. Injection is uniform in time.
//!
//! These normalized maps have coefficients in [0,1] for the admitted DAG.
//! They avoid storing time-dimensional maps or forming an inverse tiny rate.
//! The original nonnegative uniformization also integrates the source blocks;
//! no quadrature, invented source half-life or instantaneous injection is used.
//! Work, topology, cancellation and numerical admission share the original
//! owner. Zero/subnormal coefficients cannot be amplified by large injections;
//! unlike exhausted initial inventories, a continuing feed has no vanished-tail
//! exemption. Numerical range refusal is conservative, not a forward-error bound.

use super::{
    DecayAdvance, DecayEdge, DecayError, DecayLimits, DecayPlan, Sum, charge, det, poll,
    positive_or_zero,
};

/// Reusable constant-source interval and its original homogeneous propagator.
#[derive(Clone, Debug)]
pub struct DecaySourcePlan {
    homogeneous: DecayPlan,
    source_population: Vec<f64>,
    source_decays: Vec<f64>,
}

impl DecaySourcePlan {
    /// Prepare both initial-inventory and uniform-source responses. Source
    /// rates are not part of preparation; callers provide total injected counts
    /// per application and retain responsibility for a representable rate*time.
    pub fn prepare(
        rates_per_s: &[f64],
        edges: &[DecayEdge],
        duration_s: f64,
        limits: DecayLimits,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self, DecayError> {
        let (homogeneous, maps) =
            DecayPlan::prepare_impl(rates_per_s, edges, duration_s, limits, true, cancelled)?;
        let maps = maps.expect("source maps requested by this preparation");
        Ok(Self {
            homogeneous,
            source_population: maps.population,
            source_decays: maps.decays,
        })
    }

    /// Same initial-population/count owner, including complete preparation work.
    #[must_use]
    pub fn homogeneous_plan(&self) -> &DecayPlan {
        &self.homogeneous
    }

    /// Row-major Fbar: surviving population per uniformly injected constituent.
    #[must_use]
    pub fn source_population_matrix(&self) -> &[f64] {
        &self.source_population
    }

    /// Row-major Gbar: parent transformations per uniformly injected constituent.
    #[must_use]
    pub fn source_decay_matrix(&self) -> &[f64] {
        &self.source_decays
    }

    /// Apply N1=P*N0+Fbar*J and events=C*N0+Gbar*J. Every injected count is
    /// distributed uniformly over this plan's interval. A zero-length interval
    /// requires zero injection; use ordinary inventory construction for a dose.
    pub fn apply(&self, populations: &[f64], injected: &[f64]) -> Result<DecayAdvance, DecayError> {
        self.apply_with_cancel(populations, injected, || false)
    }

    /// Transactional application with bounded row checkpoints. The residual
    /// includes both initial and injected counts; inputs are never mutated.
    pub fn apply_with_cancel(
        &self,
        populations: &[f64],
        injected: &[f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<DecayAdvance, DecayError> {
        poll(&mut cancelled)?;
        let plan = &self.homogeneous;
        let n = plan.nodes;
        if populations.len() != n
            || injected.len() != n
            || !populations
                .iter()
                .chain(injected)
                .copied()
                .all(positive_or_zero)
        {
            return Err(DecayError::Invalid(
                "one finite nonnegative initial population and injection per node required",
            ));
        }
        if plan.duration_s == 0.0 && injected.iter().any(|&value| value != 0.0) {
            return Err(DecayError::Invalid(
                "zero duration cannot contain uniform injection",
            ));
        }
        let mut total = Sum::default();
        for &population in populations.iter().chain(injected) {
            total.add(population);
        }
        let total = total.total();
        if !total.is_finite() {
            return Err(DecayError::Unrepresentable(
                "initial plus injected population",
            ));
        }
        let mut answer = plan.apply_with_cancel(populations, &mut cancelled)?;
        for i in 0..n {
            poll(&mut cancelled)?;
            let mut population = Sum::default();
            let mut count = Sum::default();
            population.add(answer.populations[i]);
            count.add(answer.decays[i]);
            for (j, &injection) in injected.iter().enumerate() {
                let at = i * n + j;
                if injection > 1.0 && !plan.tails.is_empty() && plan.tails[at].reachable {
                    if self.source_population[at] == 0.0
                        || self.source_population[at].is_subnormal()
                        || self.source_decays[at].is_subnormal()
                        || (self.source_decays[at] == 0.0 && plan.rates[i] > 0.0)
                    {
                        return Err(DecayError::Unrepresentable(
                            "injection may amplify an unrepresented source coefficient",
                        ));
                    }
                }
                population.product(self.source_population[at], injection);
                count.product(self.source_decays[at], injection);
            }
            answer.populations[i] = population.total();
            answer.decays[i] = count.total();
            if !positive_or_zero(answer.populations[i]) || !positive_or_zero(answer.decays[i]) {
                return Err(DecayError::Unrepresentable(
                    "source-fed population or transformation count",
                ));
            }
        }
        let mut residual = Sum::default();
        for &population in &answer.populations {
            residual.add(population);
        }
        for &population in populations.iter().chain(injected) {
            residual.add(-population);
        }
        let residual = residual.total();
        let allowance = plan.tolerance * total;
        if !residual.is_finite() || residual.abs() > allowance {
            return Err(DecayError::Balance {
                row: n,
                column: n,
                residual,
                allowance,
            });
        }
        answer.population_balance_residual = residual;
        poll(&mut cancelled)?;
        Ok(answer)
    }
}

pub(super) struct SourceTerms {
    nodes: usize,
    population_term: Vec<f64>,
    decay_term: Vec<f64>,
    population_sum: Vec<Sum>,
    decay_sum: Vec<Sum>,
}

impl SourceTerms {
    pub(super) fn new(nodes: usize) -> Self {
        Self {
            nodes,
            population_term: vec![0.0; nodes * nodes],
            decay_term: vec![0.0; nodes * nodes],
            population_sum: vec![Sum::default(); nodes * nodes],
            decay_sum: vec![Sum::default(); nodes * nodes],
        }
    }

    // Right-multiplying the augmented uniformization powers cancels 1/mu:
    // Z_k=(mu*Z_{k-1}+pterm_{k-1})/k, W_k=(mu*W_{k-1}+cterm_{k-1})/k.
    // All four terms already include their Poisson weights. No extra GEMM.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn step(
        &mut self,
        pterm: &[f64],
        cterm: &[f64],
        mu: f64,
        k: usize,
        work: &mut u64,
        maximum_work: u64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), DecayError> {
        let weight = mu / k as f64;
        for i in 0..self.nodes {
            poll(cancelled)?;
            charge(work, 16 * (i + 1), maximum_work)?;
            for j in 0..=i {
                let at = i * self.nodes + j;
                let mut p = Sum::default();
                p.product(weight, self.population_term[at]);
                p.add(pterm[at] / k as f64);
                self.population_term[at] = p.total();
                let mut c = Sum::default();
                c.product(weight, self.decay_term[at]);
                c.add(cterm[at] / k as f64);
                self.decay_term[at] = c.total();
                self.population_sum[at].add(self.population_term[at]);
                self.decay_sum[at].add(self.decay_term[at]);
            }
        }
        Ok(())
    }

    pub(super) fn converged(&self, tail: f64, tolerance: f64) -> bool {
        (0..self.nodes * self.nodes).all(|at| {
            tail * self.population_term[at] <= tolerance * self.population_sum[at].total()
                && tail * self.decay_term[at] <= tolerance * self.decay_sum[at].total()
        })
    }

    pub(super) fn finish(self, rates: &[f64], time: f64) -> SourceMaps {
        let mut maps = SourceMaps {
            population: self.population_sum.into_iter().map(Sum::total).collect(),
            decays: self.decay_sum.into_iter().map(Sum::total).collect(),
            next_population: self.population_term,
            next_decays: self.decay_term,
        };
        maps.refresh_diagonals(rates, time);
        maps
    }
}

pub(super) struct SourceMaps {
    pub(super) population: Vec<f64>,
    pub(super) decays: Vec<f64>,
    next_population: Vec<f64>,
    next_decays: Vec<f64>,
}

impl SourceMaps {
    pub(super) fn identity(n: usize) -> Self {
        let mut population = vec![0.0; n * n];
        for i in 0..n {
            population[i * n + i] = 1.0;
        }
        Self {
            population,
            decays: vec![0.0; n * n],
            next_population: Vec::new(),
            next_decays: Vec::new(),
        }
    }

    fn refresh_diagonals(&mut self, rates: &[f64], time: f64) {
        for (i, &rate) in rates.iter().enumerate() {
            (
                self.population[i * rates.len() + i],
                self.decays[i * rates.len() + i],
            ) = source_diagonal(rate, time);
        }
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn double(
        &mut self,
        p: &[f64],
        c: &[f64],
        rates: &[f64],
        time: f64,
        work: &mut u64,
        maximum_work: u64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), DecayError> {
        let n = rates.len();
        for i in 0..n {
            poll(cancelled)?;
            charge(work, 16 * (i + 1) * (i + 2), maximum_work)?;
            for j in 0..i {
                let at = i * n + j;
                let mut f = Sum::default();
                let mut cf = Sum::default();
                f.add(self.population[at]);
                for k in j..=i {
                    f.product(p[i * n + k], self.population[k * n + j]);
                    cf.product(c[i * n + k], self.population[k * n + j]);
                }
                self.next_population[at] = 0.5 * f.total();
                let mut g = Sum::default();
                g.add(self.decays[at]);
                g.add(0.5 * cf.total());
                self.next_decays[at] = g.total();
            }
        }
        core::mem::swap(&mut self.population, &mut self.next_population);
        core::mem::swap(&mut self.decays, &mut self.next_decays);
        self.refresh_diagonals(rates, time);
        Ok(())
    }

    pub(super) fn reorder(&mut self, order: &[usize]) {
        let n = order.len();
        for i in 0..n {
            for j in 0..n {
                self.next_population[order[i] * n + order[j]] = self.population[i * n + j];
                self.next_decays[order[i] * n + order[j]] = self.decays[i * n + j];
            }
        }
        core::mem::swap(&mut self.population, &mut self.next_population);
        core::mem::swap(&mut self.decays, &mut self.next_decays);
    }
}

// phi=(1-exp(-x))/x and its independently evaluated complement. Reinsert
// these at every scale, including when a much faster unrelated node sets h.
fn source_diagonal(rate: f64, time: f64) -> (f64, f64) {
    if rate == 0.0 || time == 0.0 {
        return (1.0, 0.0);
    }
    let x = rate * time;
    if x < 0.01 {
        let lost = x
            * (0.5
                + x * (-1.0 / 6.0
                    + x * (1.0 / 24.0
                        + x * (-1.0 / 120.0
                            + x * (1.0 / 720.0 + x * (-1.0 / 5040.0 + x / 40320.0))))));
        (1.0 - lost, lost)
    } else {
        let retained = if x.is_infinite() {
            (1.0 / rate) / time
        } else {
            -det::expm1(-x) / x
        };
        (retained, 1.0 - retained)
    }
}

#[cfg(test)]
#[path = "source/tests.rs"]
mod tests;
