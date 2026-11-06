//! Constant production through the original positive reaction/count owner.
//!
//! For an interval of length t, Fbar = integral(P)/t and Gbar = integral(C)/t
//! map uniformly injected totals J into surviving populations and removals.
//! Thus N1=P*N0+Fbar*J and removals=C*N0+Gbar*J retain actual production ages,
//! repeated visits in cycles, finite sinks and declared product multiplicity.
//! These maps need not have coefficients below one or conserve population.
//!
//! Source terms reuse the same nonnegative shifted powers and semigroup
//! composition as `ReactionPlan`; no inverse generator, quadrature or added
//! numerical nuclides are used. Noncyclic scalar source diagonals are restored
//! at every scale. Source rates are not part of preparation, which is reusable
//! across independent initial inventories and nonnegative injection vectors.
//! A continuing source has fresh material even after old inventories exhaust:
//! unrepresented source coefficients receive no vanished-tail exemption.

use super::{
    DecayError, DecayLimits, ReactionAdvance, ReactionEdge, ReactionPlan, Sum, charge, det, poll,
    positive_or_zero,
};

/// Reusable unchanged-rate reaction history with uniform external production.
#[derive(Clone, Debug)]
pub struct ReactionSourcePlan {
    homogeneous: ReactionPlan,
    source_population: Vec<f64>,
    source_removals: Vec<f64>,
}

impl ReactionSourcePlan {
    /// Prepare both responses under the original node/edge/work limits. Callers
    /// supply representable integrated injection counts per application, so no
    /// source vector or fabricated source species becomes part of this map.
    pub fn prepare(
        rates_per_s: &[f64],
        edges: &[ReactionEdge],
        duration_s: f64,
        limits: DecayLimits,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self, DecayError> {
        let (homogeneous, maps) =
            ReactionPlan::prepare_impl(rates_per_s, edges, duration_s, limits, true, cancelled)?;
        let maps = maps.expect("source maps requested by this preparation");
        Ok(Self {
            homogeneous,
            source_population: maps.population,
            source_removals: maps.removals,
        })
    }

    /// Original population/removal maps, including the full preparation work.
    #[must_use]
    pub fn homogeneous_plan(&self) -> &ReactionPlan {
        &self.homogeneous
    }

    /// Row-major Fbar: endpoint populations per uniformly injected constituent.
    #[must_use]
    pub fn source_population_matrix(&self) -> &[f64] {
        &self.source_population
    }

    /// Row-major Gbar: integrated removals per uniformly injected constituent.
    #[must_use]
    pub fn source_removal_matrix(&self) -> &[f64] {
        &self.source_removals
    }

    /// Apply to initial populations and totals injected uniformly over the
    /// interval. Zero duration requires zero injection; an instantaneous dose
    /// belongs in the initial inventory rather than this constant source.
    pub fn apply(
        &self,
        populations: &[f64],
        injected: &[f64],
    ) -> Result<ReactionAdvance, DecayError> {
        self.apply_with_cancel(populations, injected, || false)
    }

    /// Evaluate an independent candidate with bounded row checkpoints and
    /// original-order, source-aware reaction balances. Inputs are unchanged.
    pub fn apply_with_cancel(
        &self,
        populations: &[f64],
        injected: &[f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ReactionAdvance, DecayError> {
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
                "one finite nonnegative initial population and injection per reaction node required",
            ));
        }
        if injected.iter().all(|&value| value == 0.0) {
            return plan.apply_with_cancel(populations, cancelled);
        }
        if plan.duration_s == 0.0 {
            return Err(DecayError::Invalid(
                "zero duration cannot contain uniform reaction injection",
            ));
        }
        let mut answer = plan.apply_with_cancel(populations, &mut cancelled)?;
        for i in 0..n {
            poll(&mut cancelled)?;
            let mut population = Sum::default();
            let mut removals = Sum::default();
            population.add(answer.populations[i]);
            removals.add(answer.removals[i]);
            for (j, &injection) in injected.iter().enumerate() {
                let at = i * n + j;
                if injection > 1.0 && plan.reachable[at] {
                    let p = self.source_population[at];
                    let c = self.source_removals[at];
                    if p == 0.0
                        || p.is_subnormal()
                        || c.is_subnormal()
                        || (c == 0.0 && plan.rates[i] > 0.0)
                    {
                        return Err(DecayError::Unrepresentable(
                            "injection may amplify an unrepresented reaction source coefficient",
                        ));
                    }
                }
                population.product(self.source_population[at], injection);
                removals.product(self.source_removals[at], injection);
            }
            answer.populations[i] = population.total();
            answer.removals[i] = removals.total();
        }
        // Keep initial and supplied counts separate: adding a tiny feed to a
        // large initial inventory first can erase its entire balance debit.
        (answer.balance_residuals, answer.population_balance_residual) = plan.check_balance(
            populations,
            Some(injected),
            &answer.populations,
            &answer.removals,
            n,
        )?;
        poll(&mut cancelled)?;
        Ok(answer)
    }
}

pub(super) struct SourceTerms {
    nodes: usize,
    population_term: Vec<f64>,
    removal_term: Vec<f64>,
    population_sum: Vec<Sum>,
    removal_sum: Vec<Sum>,
}

impl SourceTerms {
    pub(super) fn new(nodes: usize) -> Self {
        Self {
            nodes,
            population_term: vec![0.0; nodes * nodes],
            removal_term: vec![0.0; nodes * nodes],
            population_sum: vec![Sum::default(); nodes * nodes],
            removal_sum: vec![Sum::default(); nodes * nodes],
        }
    }

    // Right-multiplication of the shifted block powers needs no new GEMM:
    // Z_k=(mu*Z_{k-1}+Pterm_{k-1})/k, W_k=(mu*W_{k-1}+Cterm_{k-1})/k.
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
        for i in 0..self.nodes {
            poll(cancelled)?;
            charge(work, 16 * self.nodes, maximum_work)?;
            for j in 0..self.nodes {
                let at = i * self.nodes + j;
                let mut p = Sum::default();
                p.product(mu, self.population_term[at]);
                p.add(pterm[at]);
                self.population_term[at] = p.total() / k as f64;
                let mut c = Sum::default();
                c.product(mu, self.removal_term[at]);
                c.add(cterm[at]);
                self.removal_term[at] = c.total() / k as f64;
                self.population_sum[at].add(self.population_term[at]);
                self.removal_sum[at].add(self.removal_term[at]);
            }
        }
        Ok(())
    }

    pub(super) fn converged(&self, tail: f64, tolerance: f64) -> bool {
        self.population_sum
            .iter()
            .chain(&self.removal_sum)
            .all(|entry| {
                let value = entry.total();
                value == 0.0 || tail <= tolerance * value
            })
    }

    pub(super) fn finish(
        self,
        rates: &[f64],
        net_rates: &[f64],
        cyclic: &[bool],
        time: f64,
        damping: f64,
    ) -> Result<SourceMaps, DecayError> {
        let mut maps = SourceMaps {
            population: self
                .population_sum
                .into_iter()
                .map(|entry| damping * entry.total())
                .collect(),
            removals: self
                .removal_sum
                .into_iter()
                .map(|entry| damping * entry.total())
                .collect(),
            next_population: self.population_term,
            next_removals: self.removal_term,
        };
        maps.refresh_diagonals(rates, net_rates, cyclic, time)?;
        Ok(maps)
    }
}

pub(super) struct SourceMaps {
    pub(super) population: Vec<f64>,
    pub(super) removals: Vec<f64>,
    next_population: Vec<f64>,
    next_removals: Vec<f64>,
}

impl SourceMaps {
    pub(super) fn identity(n: usize) -> Self {
        let mut population = vec![0.0; n * n];
        for i in 0..n {
            population[i * n + i] = 1.0;
        }
        Self {
            population,
            removals: vec![0.0; n * n],
            next_population: Vec::new(),
            next_removals: Vec::new(),
        }
    }

    fn refresh_diagonals(
        &mut self,
        rates: &[f64],
        net_rates: &[f64],
        cyclic: &[bool],
        time: f64,
    ) -> Result<(), DecayError> {
        for i in 0..rates.len() {
            if !cyclic[i] {
                (
                    self.population[i * rates.len() + i],
                    self.removals[i * rates.len() + i],
                ) = source_diagonal(rates[i], net_rates[i], time)?;
            }
        }
        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn double(
        &mut self,
        p: &[f64],
        c: &[f64],
        rates: &[f64],
        net_rates: &[f64],
        cyclic: &[bool],
        time: f64,
        work: &mut u64,
        maximum_work: u64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), DecayError> {
        let n = rates.len();
        for i in 0..n {
            poll(cancelled)?;
            charge(work, n * (2 * n + 4), maximum_work)?;
            for j in 0..n {
                let at = i * n + j;
                let mut f = Sum::default();
                let mut cf = Sum::default();
                f.add(self.population[at]);
                for k in 0..n {
                    f.product(p[i * n + k], self.population[k * n + j]);
                    cf.product(c[i * n + k], self.population[k * n + j]);
                }
                self.next_population[at] = 0.5 * f.total();
                let mut g = Sum::default();
                g.add(self.removals[at]);
                g.add(0.5 * cf.total());
                self.next_removals[at] = g.total();
            }
        }
        core::mem::swap(&mut self.population, &mut self.next_population);
        core::mem::swap(&mut self.removals, &mut self.next_removals);
        self.refresh_diagonals(rates, net_rates, cyclic, time)
    }
}

// Fbar=phi1(-r*t), Gbar=lambda*t*phi2(-r*t), with r=lambda-production_self.
// Independent phi2 evaluation retains tiny removal counts and net-zero self
// returns. Large positive exposures avoid forming lambda*t or an infinite r*t.
fn source_diagonal(rate: f64, net_rate: f64, time: f64) -> Result<(f64, f64), DecayError> {
    let x = net_rate * time;
    let (population, removals) = if x.abs() < 0.01 {
        let phi2 = 0.5
            + x * (-1.0 / 6.0
                + x * (1.0 / 24.0
                    + x * (-1.0 / 120.0 + x * (1.0 / 720.0 + x * (-1.0 / 5040.0 + x / 40320.0)))));
        (1.0 - x * phi2, (rate * time) * phi2)
    } else {
        let phi1 = if x.is_infinite() && x > 0.0 {
            (1.0 / net_rate) / time
        } else {
            -det::expm1(-x) / x
        };
        (phi1, (rate / net_rate) * (1.0 - phi1))
    };
    if !positive_or_zero(population) || !positive_or_zero(removals) {
        return Err(DecayError::Unrepresentable(
            "scalar reaction-source population or removal count",
        ));
    }
    Ok((population, removals))
}

#[cfg(test)]
#[path = "source/tests.rs"]
mod tests;
