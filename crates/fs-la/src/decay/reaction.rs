//! Prepared first-order reactions with explicit removal and product multiplicity.
//!
//! For column populations, A = (Y - I) diag(lambda). An edge supplies one
//! entry of Y; duplicate entries add. Columns need not sum to one: finite
//! sinks, cycles and several products per removal are represented without
//! renormalizing their yields. The caller owns reaction identities and units.
//!
//! Preparation retains P = exp(A t) and C = integral diag(lambda) exp(A s) ds.
//! Applying C counts removals of daughters born DURING the interval. A caller
//! with competing channels partitions C_i by the channel rate / lambda_i;
//! it must not infer events by subtracting endpoint populations.
//! Optional population integrals J = integral exp(A s) ds [s] reuse the
//! normalized population map from the uniform-production owner below. Stable
//! receiving nodes retain their residence without an artificial removal clock
//! or division of their zero removal count. Ordinary preparation is unchanged.
//!
//! Nonnegative uniformization of the population/count block is followed by
//! P(2h)=P(h)P(h), C(2h)=C(h)+C(h)P(h). Scalar survival and count diagonals
//! are refreshed at every scale for nodes with no nontrivial return cycle.
//! This preserves slow removal in the presence of unrelated stiff components.
//! Cyclic blocks retain their full return paths and never receive that scalar
//! replacement. The positive-series tail test and final node/total balances
//! are numerical admission, not a certified floating-point forward-error bound.
//! Extreme cyclic stiffness, multiplication or unrepresented transfers may
//! refuse; no clipping, yield normalization or invented conservation is used.

use super::{
    DecayError, DecayLimits, MAX_DECAY_EDGES, MAX_DECAY_NODES, MAX_SERIES_TERMS, MAX_SQUARINGS,
    Sum, cap, charge, det, poll, positive_or_zero,
};

/// Uniform production through the same positive reaction/count propagation.
#[path = "reaction/source.rs"]
pub mod source;

/// An expected number of products per removal of a parent. Repeated endpoints
/// add, self-products are allowed, and a parent with no edges is a finite sink.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ReactionEdge {
    pub parent: usize,
    pub daughter: usize,
    pub yield_per_removal: f64,
}

/// Endpoint populations and counts from one unchanged-rate interval.
#[derive(Clone, Debug, PartialEq)]
pub struct ReactionAdvance {
    pub populations: Vec<f64>,
    /// Integral lambda_i N_i(s) ds, including repeated visits in a cycle.
    pub removals: Vec<f64>,
    /// N1_i - N0_i - injected_i + removals_i - sum_j Y_ij removals_j,
    /// uncorrected. The homogeneous plan has zero injected counts.
    pub balance_residuals: Vec<f64>,
    /// Total population change less injected counts and
    /// sum_i (sum_j Y_ji - 1) removals_i.
    /// This does not assume population conservation for multiplying reactions.
    pub population_balance_residual: f64,
}

/// Immutable bounded dense maps, reusable across independent inventories.
#[derive(Clone, Debug)]
pub struct ReactionPlan {
    nodes: usize,
    duration_s: f64,
    rates: Vec<f64>,
    yields: Vec<f64>,
    net_yields: Vec<f64>,
    transition: Vec<f64>,
    integrated: Vec<f64>,
    population_integral: Option<Vec<f64>>,
    reachable: Vec<bool>,
    tails: Vec<Tail>,
    squarings: usize,
    series_terms: usize,
    work: u64,
    tolerance: f64,
}

#[derive(Clone, Copy, Debug)]
struct Tail {
    minimum_rate: f64,
    steps: usize,
    log_weight: f64,
}

impl Tail {
    const UNKNOWN: Self = Self {
        minimum_rate: 0.0,
        steps: 0,
        log_weight: 0.0,
    };

    // For an acyclic path, the remaining transient population is bounded by
    // the path's multiplicity times an Erlang survival tail. Summing path
    // weights also bounds a network with yields greater than one. Cyclic paths
    // receive no exhaustion exemption. The margin is deliberately conservative.
    fn vanishes(self, population: f64, duration: f64) -> bool {
        if self.minimum_rate <= 0.0 || !self.log_weight.is_finite() {
            return false;
        }
        let exposure = self.minimum_rate * duration;
        if exposure.is_infinite() {
            return true;
        }
        if exposure <= 0.0 || exposure < self.steps as f64 {
            return false;
        }
        let mut sum = 1.0;
        let mut term = 1.0;
        let mut log_factorial = 0.0;
        for k in (1..=self.steps).rev() {
            term *= k as f64 / exposure;
            sum += term;
            log_factorial += det::ln(k as f64);
        }
        det::ln(population) + self.log_weight - exposure + self.steps as f64 * det::ln(exposure)
            - log_factorial
            + det::ln(sum)
            < -746.0
    }
}

impl ReactionPlan {
    /// Prepare nonnegative, finite, unchanged removal rates and product yields.
    /// Zero-rate nodes have no active production even if an inactive edge was
    /// supplied. Zero duration gives identity populations and zero counts.
    /// All resources, numerical admission and errors use the original decay
    /// owner's limits; no existing closed/acyclic contract is relaxed.
    pub fn prepare(
        removal_rates_per_s: &[f64],
        edges: &[ReactionEdge],
        duration_s: f64,
        limits: DecayLimits,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self, DecayError> {
        Self::prepare_impl(
            removal_rates_per_s,
            edges,
            duration_s,
            limits,
            false,
            false,
            cancelled,
        )
        .map(|(plan, _)| plan)
    }

    /// Retain J = integral exp(A s) ds [s], including zero-removal receiving
    /// nodes. The existing uniform-source series, scalar diagonals, semigroup
    /// and numerical admission supply J/t; this wrapper scales that same map.
    /// Source-removal preparation and the additional scaling/application work
    /// are explicit costs of this opt-in path. Ordinary/source preparation
    /// retain their original behavior and work limits.
    pub fn prepare_with_population_integrals(
        removal_rates_per_s: &[f64],
        edges: &[ReactionEdge],
        duration_s: f64,
        limits: DecayLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, DecayError> {
        let (mut plan, maps) = Self::prepare_impl(
            removal_rates_per_s,
            edges,
            duration_s,
            limits,
            true,
            false,
            &mut cancelled,
        )?;
        let integral = maps
            .expect("source map requested for population integrals")
            .population;
        let n = plan.nodes;
        charge(&mut plan.work, 6 * n * n, limits.maximum_work)?;
        plan.population_integral =
            Some(plan.scale_population_integral(integral, &plan.integrated, &mut cancelled)?);
        poll(&mut cancelled)?;
        Ok(plan)
    }

    // Both J/t and K/t^2 become population-time maps after multiplying by
    // duration. Their independent removal maps admit the same row identity.
    fn scale_population_integral(
        &self,
        mut integral: Vec<f64>,
        removals: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<f64>, DecayError> {
        let n = self.nodes;
        let duration_s = self.duration_s;
        for i in 0..n {
            poll(cancelled)?;
            for j in 0..n {
                let at = i * n + j;
                let coefficient = integral[at];
                // Duration is an extensive multiplier just as injected count
                // is for ReactionSourcePlan. Do not amplify vanished or
                // subnormal normalized entries into an apparently resolved J.
                if duration_s > 1.0
                    && self.reachable[at]
                    && (coefficient == 0.0 || coefficient.is_subnormal())
                {
                    return Err(DecayError::Unrepresentable(
                        "duration may amplify an unrepresented population integral coefficient",
                    ));
                }
                let value = duration_s * coefficient;
                if !positive_or_zero(value)
                    || (duration_s > 0.0 && coefficient > 0.0 && value == 0.0)
                {
                    return Err(DecayError::Unrepresentable(
                        "dimensional population integral",
                    ));
                }
                integral[at] = value;
                // Independently integrated removals must remain consistent
                // with residence where a real removal rate exists. Zero-rate
                // rows keep the source owner's direct positive population map.
                let mut residual = Sum::default();
                residual.add(removals[at]);
                residual.product(-self.rates[i], value);
                let mut scale = Sum::default();
                scale.add(removals[at]);
                scale.product(self.rates[i], value);
                admit_balance(i, j, residual.total(), scale.total(), self.tolerance)?;
            }
        }
        Ok(integral)
    }

    fn prepare_impl(
        removal_rates_per_s: &[f64],
        edges: &[ReactionEdge],
        duration_s: f64,
        limits: DecayLimits,
        with_source: bool,
        with_source_integrals: bool,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(Self, Option<source::SourceMaps>), DecayError> {
        poll(&mut cancelled)?;
        let n = removal_rates_per_s.len();
        if limits.maximum_nodes == 0
            || limits.maximum_nodes > MAX_DECAY_NODES
            || limits.maximum_edges > MAX_DECAY_EDGES
            || limits.maximum_squarings > MAX_SQUARINGS
            || limits.maximum_series_terms == 0
            || limits.maximum_series_terms > MAX_SERIES_TERMS
            || limits.maximum_work == 0
            || !limits.relative_tolerance.is_finite()
            || limits.relative_tolerance <= 0.0
            || limits.relative_tolerance >= 1.0
        {
            return Err(DecayError::Invalid(
                "finite bounded reaction limits required",
            ));
        }
        cap("reaction nodes", n, limits.maximum_nodes)?;
        cap("reaction edges", edges.len(), limits.maximum_edges)?;
        if n == 0
            || !positive_or_zero(duration_s)
            || !removal_rates_per_s.iter().copied().all(positive_or_zero)
        {
            return Err(DecayError::Invalid(
                "nonempty finite nonnegative rates and duration required",
            ));
        }
        let mut work = 0;
        charge(&mut work, n * n + edges.len(), limits.maximum_work)?;
        let mut sums = vec![Sum::default(); n * n];
        for edge in edges {
            poll(&mut cancelled)?;
            if edge.parent >= n
                || edge.daughter >= n
                || !edge.yield_per_removal.is_finite()
                || edge.yield_per_removal <= 0.0
            {
                return Err(DecayError::Invalid(
                    "valid reaction endpoints and positive yields required",
                ));
            }
            sums[edge.daughter * n + edge.parent].add(edge.yield_per_removal);
        }
        let yields: Vec<_> = sums.into_iter().map(Sum::total).collect();
        if !yields.iter().copied().all(positive_or_zero) {
            return Err(DecayError::Unrepresentable("summed reaction yield"));
        }
        let mut net_yields = vec![0.0; n];
        for j in 0..n {
            let mut sum = Sum::default();
            sum.add(-1.0);
            for i in 0..n {
                sum.add(yields[i * n + j]);
            }
            net_yields[j] = sum.total();
            if !net_yields[j].is_finite() {
                return Err(DecayError::Unrepresentable("net reaction multiplicity"));
            }
        }
        let mut plan = Self {
            nodes: n,
            duration_s,
            rates: removal_rates_per_s.to_vec(),
            yields,
            net_yields,
            transition: vec![0.0; n * n],
            integrated: vec![0.0; n * n],
            population_integral: None,
            reachable: vec![false; n * n],
            tails: vec![Tail::UNKNOWN; n * n],
            squarings: 0,
            series_terms: 0,
            work,
            tolerance: limits.relative_tolerance,
        };
        for i in 0..n {
            plan.transition[i * n + i] = 1.0;
            plan.reachable[i * n + i] = true;
        }
        let shift = removal_rates_per_s.iter().copied().fold(0.0, f64::max);
        if duration_s == 0.0 || shift == 0.0 {
            if with_source {
                charge(
                    &mut plan.work,
                    (2 + usize::from(with_source_integrals)) * n * n,
                    limits.maximum_work,
                )?;
            }
            poll(&mut cancelled)?;
            return Ok((
                plan,
                with_source.then(|| source::SourceMaps::identity(n, with_source_integrals)),
            ));
        }
        let mut production = vec![0.0; n * n];
        let mut maximum_production = 0.0_f64;
        for j in 0..n {
            poll(&mut cancelled)?;
            let mut sum = Sum::default();
            for i in 0..n {
                let at = i * n + j;
                production[at] = plan.yields[at] * plan.rates[j];
                if !positive_or_zero(production[at])
                    || (plan.yields[at] > 0.0 && plan.rates[j] > 0.0 && production[at] == 0.0)
                {
                    return Err(DecayError::Unrepresentable("reaction production rate"));
                }
                if production[at] > 0.0 {
                    plan.reachable[at] = true;
                }
                sum.add(production[at]);
            }
            let total = sum.total();
            if !positive_or_zero(total) {
                return Err(DecayError::Unrepresentable(
                    "total reaction production rate",
                ));
            }
            maximum_production = maximum_production.max(total);
        }
        // The 1-norm of the shifted population/count block is no larger than
        // h*(shift + maximum column production). Halving also handles a rate
        // times duration that overflows before scaling; no inverse tiny rate.
        let mut time = duration_s;
        while shift * time + maximum_production * time > 0.5 {
            poll(&mut cancelled)?;
            cap(
                "reaction squarings",
                plan.squarings + 1,
                limits.maximum_squarings,
            )?;
            let half = time * 0.5;
            if half == 0.0 || half + half != time {
                return Err(DecayError::Unrepresentable("reaction base interval"));
            }
            time = half;
            plan.squarings += 1;
        }
        charge(&mut plan.work, n * n * n, limits.maximum_work)?;
        for k in 0..n {
            poll(&mut cancelled)?;
            for i in 0..n {
                for j in 0..n {
                    plan.reachable[i * n + j] |=
                        plan.reachable[i * n + k] && plan.reachable[k * n + j];
                }
            }
        }
        let mut cyclic = vec![false; n];
        for i in 0..n {
            cyclic[i] =
                (0..n).any(|j| j != i && plan.reachable[i * n + j] && plan.reachable[j * n + i]);
        }
        let net_rates: Vec<_> = (0..n)
            .map(|i| plan.rates[i] - production[i * n + i])
            .collect();
        prepare_tails(
            &mut plan,
            &production,
            &net_rates,
            &cyclic,
            limits,
            &mut cancelled,
        )?;
        let mu = shift * time;
        let mut scaled = vec![0.0; n * n];
        let mut scaled_removal = vec![0.0; n];
        for i in 0..n {
            scaled_removal[i] = plan.rates[i] * time;
            if plan.rates[i] > 0.0 && scaled_removal[i] == 0.0 {
                return Err(DecayError::Unrepresentable("scaled direct removal"));
            }
            for j in 0..n {
                let at = i * n + j;
                scaled[at] = production[at] * time;
                if production[at] > 0.0 && scaled[at] == 0.0 {
                    return Err(DecayError::Unrepresentable(
                        "scaled direct reaction transfer",
                    ));
                }
                if i == j {
                    scaled[at] += (shift - plan.rates[i]) * time;
                }
            }
        }
        let mut block_norm = mu;
        for j in 0..n {
            let mut sum = Sum::default();
            sum.add(scaled_removal[j]);
            for i in 0..n {
                sum.add(scaled[i * n + j]);
            }
            block_norm = block_norm.max(sum.total());
        }
        if !positive_or_zero(block_norm) || block_norm > 0.500_000_000_000_001 {
            return Err(DecayError::Unrepresentable(
                "scaled positive reaction operator",
            ));
        }
        let incoming: Vec<Vec<_>> = (0..n)
            .map(|i| {
                (0..n)
                    .filter_map(|j| (scaled[i * n + j] > 0.0).then_some((j, scaled[i * n + j])))
                    .collect()
            })
            .collect();
        let nonzero_terms: usize = incoming.iter().map(Vec::len).sum();
        let mut term_p = plan.transition.clone();
        let mut term_c = vec![0.0; n * n];
        let mut next_p = vec![0.0; n * n];
        let mut next_c = vec![0.0; n * n];
        let mut sum_p = vec![Sum::default(); n * n];
        let mut sum_c = vec![Sum::default(); n * n];
        for i in 0..n {
            sum_p[i * n + i].add(1.0);
        }
        let mut source_terms = if with_source {
            charge(
                &mut plan.work,
                (6 + 3 * usize::from(with_source_integrals)) * n * n,
                limits.maximum_work,
            )?;
            Some(source::SourceTerms::new(n, with_source_integrals))
        } else {
            None
        };
        let mut norm_term = 1.0;
        // A tighter base-series allowance accounts for the number of repeated
        // compositions. Roundoff is still assessed by the final balances.
        let mut base_tolerance = limits.relative_tolerance * 0.125;
        for _ in 0..plan.squarings {
            base_tolerance *= 0.5;
        }
        let mut converged = false;
        for k in 1..=limits.maximum_series_terms {
            if let Some(source) = &mut source_terms {
                source.step(
                    &term_p,
                    &term_c,
                    mu,
                    k,
                    &mut plan.work,
                    limits.maximum_work,
                    &mut cancelled,
                )?;
            }
            charge(
                &mut plan.work,
                n * (nonzero_terms + 2 * n),
                limits.maximum_work,
            )?;
            for i in 0..n {
                poll(&mut cancelled)?;
                for j in 0..n {
                    let at = i * n + j;
                    let mut product = Sum::default();
                    for &(p, coefficient) in &incoming[i] {
                        product.product(coefficient, term_p[p * n + j]);
                    }
                    next_p[at] = product.total() / k as f64;
                    let mut count = Sum::default();
                    count.product(scaled_removal[i], term_p[at]);
                    count.product(mu, term_c[at]);
                    next_c[at] = count.total() / k as f64;
                    sum_p[at].add(next_p[at]);
                    sum_c[at].add(next_c[at]);
                }
            }
            core::mem::swap(&mut term_p, &mut next_p);
            core::mem::swap(&mut term_c, &mut next_c);
            plan.series_terms = k;
            let previous_norm_term = norm_term;
            norm_term *= block_norm / k as f64;
            let tail =
                norm_term * (block_norm / (k + 1) as f64) / (1.0 - block_norm / (k + 2) as f64);
            // Requiring at least n+1 terms exposes every shortest path in both
            // blocks. A common absolute tail compared to every represented
            // positive entry avoids dropping small but resolved branches.
            if k >= n + 1 + usize::from(with_source) + usize::from(with_source_integrals)
                && sum_p.iter().chain(&sum_c).all(|entry| {
                    let value = entry.total();
                    value == 0.0 || tail <= base_tolerance * value
                })
                && source_terms.as_ref().is_none_or(|source| {
                    // An integrated block's kth term is bounded by
                    // block_norm^(k-1)/(k-1)!, one order below P/C.
                    let source_tail = norm_term / (1.0 - block_norm / (k + 1) as f64);
                    // The second integral's first omitted term is bounded
                    // by block_norm^(k-1)/(2*(k-1)!). Retaining the previous
                    // norm term avoids dividing by a tiny or zero norm.
                    let residence_tail = 0.5 * previous_norm_term / (1.0 - block_norm / k as f64);
                    source.converged(source_tail, residence_tail, base_tolerance)
                })
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
        let damping = det::exp(-mu);
        plan.transition = sum_p
            .into_iter()
            .map(|entry| damping * entry.total())
            .collect();
        plan.integrated = sum_c
            .into_iter()
            .map(|entry| damping * entry.total())
            .collect();
        restore_diagonals(&mut plan, &net_rates, &cyclic, time)?;
        if with_source {
            charge(
                &mut plan.work,
                (2 + usize::from(with_source_integrals)) * n * n
                    + (1 + usize::from(with_source_integrals)) * n,
                limits.maximum_work,
            )?;
        }
        let mut source_maps = source_terms
            .map(|source| source.finish(&plan.rates, &net_rates, &cyclic, time, damping))
            .transpose()?;
        for _ in 0..plan.squarings {
            if let Some(source) = &mut source_maps {
                source.double(
                    &plan.transition,
                    &plan.integrated,
                    &plan.rates,
                    &net_rates,
                    &cyclic,
                    time + time,
                    &mut plan.work,
                    limits.maximum_work,
                    &mut cancelled,
                )?;
            }
            charge(&mut plan.work, 2 * n * n * n, limits.maximum_work)?;
            for i in 0..n {
                poll(&mut cancelled)?;
                for j in 0..n {
                    let mut population = Sum::default();
                    let mut count = Sum::default();
                    count.add(plan.integrated[i * n + j]);
                    for k in 0..n {
                        population.product(plan.transition[i * n + k], plan.transition[k * n + j]);
                        count.product(plan.integrated[i * n + k], plan.transition[k * n + j]);
                    }
                    next_p[i * n + j] = population.total();
                    next_c[i * n + j] = count.total();
                }
            }
            core::mem::swap(&mut plan.transition, &mut next_p);
            core::mem::swap(&mut plan.integrated, &mut next_c);
            time += time;
            restore_diagonals(&mut plan, &net_rates, &cyclic, time)?;
        }
        charge(&mut plan.work, n * n * (n + 3), limits.maximum_work)?;
        for j in 0..n {
            poll(&mut cancelled)?;
            let initial: Vec<_> = (0..n).map(|i| if i == j { 1.0 } else { 0.0 }).collect();
            let populations: Vec<_> = (0..n).map(|i| plan.transition[i * n + j]).collect();
            let removals: Vec<_> = (0..n).map(|i| plan.integrated[i * n + j]).collect();
            plan.check_balance(&initial, None, &populations, &removals, j)?;
        }
        if let Some(source) = &source_maps {
            charge(&mut plan.work, n * n * (n + 3), limits.maximum_work)?;
            for j in 0..n {
                poll(&mut cancelled)?;
                let injected: Vec<_> = (0..n).map(|i| if i == j { 1.0 } else { 0.0 }).collect();
                let populations: Vec<_> = (0..n).map(|i| source.population[i * n + j]).collect();
                let removals: Vec<_> = (0..n).map(|i| source.removals[i * n + j]).collect();
                // A source column satisfies Fbar-I+(I-Y)Gbar=0, the same
                // unit-input balance as one homogeneous map column.
                plan.check_balance(&injected, None, &populations, &removals, j)?;
            }
        }
        poll(&mut cancelled)?;
        Ok((plan, source_maps))
    }

    #[must_use]
    pub const fn node_count(&self) -> usize {
        self.nodes
    }
    #[must_use]
    pub const fn duration_s(&self) -> f64 {
        self.duration_s
    }
    /// Row-major P: endpoint populations per initial constituent.
    #[must_use]
    pub fn transition_matrix(&self) -> &[f64] {
        &self.transition
    }
    /// Row-major C: independently integrated removals per initial constituent.
    #[must_use]
    pub fn integrated_removal_matrix(&self) -> &[f64] {
        &self.integrated
    }
    /// Optional row-major J [s]: integrated population per initial constituent.
    /// Includes residence in stable receiving nodes. Use the checked application
    /// method when a large inventory could amplify an unrepresented coefficient.
    #[must_use]
    pub fn population_integral_matrix(&self) -> Option<&[f64]> {
        self.population_integral.as_deref()
    }
    #[must_use]
    pub const fn squarings(&self) -> usize {
        self.squarings
    }
    #[must_use]
    pub const fn series_terms(&self) -> usize {
        self.series_terms
    }
    #[must_use]
    pub const fn preparation_work(&self) -> u64 {
        self.work
    }

    pub fn apply(&self, populations: &[f64]) -> Result<ReactionAdvance, DecayError> {
        self.apply_with_cancel(populations, || false)
    }

    /// Apply without mutating input. Large populations may not silently amplify
    /// an underflowed map coefficient; only proved exhausted DAG tails are exempt.
    pub fn apply_with_cancel(
        &self,
        populations: &[f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ReactionAdvance, DecayError> {
        poll(&mut cancelled)?;
        let n = self.nodes;
        if populations.len() != n || !populations.iter().copied().all(positive_or_zero) {
            return Err(DecayError::Invalid(
                "one finite nonnegative population per reaction node required",
            ));
        }
        let mut endpoint = vec![0.0; n];
        let mut removals = vec![0.0; n];
        for i in 0..n {
            poll(&mut cancelled)?;
            let mut amount = Sum::default();
            let mut count = Sum::default();
            for (j, &initial) in populations.iter().enumerate() {
                let at = i * n + j;
                if self.duration_s > 0.0 && initial > 1.0 && self.reachable[at] {
                    let p = self.transition[at];
                    if (p == 0.0 || p.is_subnormal())
                        && !(p == 0.0 && self.tails[at].vanishes(initial, self.duration_s))
                    {
                        return Err(DecayError::Unrepresentable(
                            "inventory may amplify an unrepresented reaction population",
                        ));
                    }
                    let c = self.integrated[at];
                    if self.rates[i] > 0.0 && (c == 0.0 || c.is_subnormal()) {
                        return Err(DecayError::Unrepresentable(
                            "inventory may amplify an unrepresented removal count",
                        ));
                    }
                }
                amount.product(self.transition[at], initial);
                count.product(self.integrated[at], initial);
            }
            endpoint[i] = amount.total();
            removals[i] = count.total();
        }
        let (balance_residuals, population_balance_residual) =
            self.check_balance(populations, None, &endpoint, &removals, n)?;
        poll(&mut cancelled)?;
        Ok(ReactionAdvance {
            populations: endpoint,
            removals,
            balance_residuals,
            population_balance_residual,
        })
    }

    /// Apply the optional J map to obtain each node's population-time integral.
    pub fn integrate_populations(&self, populations: &[f64]) -> Result<Vec<f64>, DecayError> {
        self.integrate_populations_with_cancel(populations, || false)
    }

    /// Return J*N0 without changing the plan or inventory. A zero-removal node
    /// may have a nonzero integral even when its integrated removal count is zero.
    /// Cumulative residence has no exhausted-tail exemption: missing positive
    /// coefficients or products refuse instead of silently erasing fluence.
    pub fn integrate_populations_with_cancel(
        &self,
        populations: &[f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Vec<f64>, DecayError> {
        poll(&mut cancelled)?;
        let integral = self
            .population_integral
            .as_ref()
            .ok_or(DecayError::Invalid(
                "population integral maps were not prepared",
            ))?;
        let n = self.nodes;
        if populations.len() != n || !populations.iter().copied().all(positive_or_zero) {
            return Err(DecayError::Invalid(
                "one finite nonnegative population per reaction node required",
            ));
        }
        let mut answer = Vec::with_capacity(n);
        for i in 0..n {
            poll(&mut cancelled)?;
            let mut residence = Sum::default();
            self.add_population_integral(&mut residence, integral, i, populations)?;
            let value = residence.total();
            if !positive_or_zero(value) {
                return Err(DecayError::Unrepresentable("integrated population"));
            }
            answer.push(value);
        }
        poll(&mut cancelled)?;
        Ok(answer)
    }

    // Shared checked accumulation for initial and uniformly injected residence.
    // The two vectors stay separate, including when a small feed would vanish
    // if it were first added to a much larger initial population.
    fn add_population_integral(
        &self,
        sum: &mut Sum,
        integral: &[f64],
        row: usize,
        populations: &[f64],
    ) -> Result<(), DecayError> {
        for (j, &population) in populations.iter().enumerate() {
            let at = row * self.nodes + j;
            let coefficient = integral[at];
            if self.duration_s > 0.0
                && population > 0.0
                && self.reachable[at]
                && (coefficient == 0.0 || (population > 1.0 && coefficient.is_subnormal()))
            {
                return Err(DecayError::Unrepresentable(
                    "inventory requires an unrepresented population integral coefficient",
                ));
            }
            let contribution = coefficient * population;
            if !positive_or_zero(contribution)
                || (coefficient > 0.0 && population > 0.0 && contribution == 0.0)
            {
                return Err(DecayError::Unrepresentable(
                    "integrated population contribution",
                ));
            }
            sum.product(coefficient, population);
        }
        Ok(())
    }

    fn check_balance(
        &self,
        initial: &[f64],
        injected: Option<&[f64]>,
        populations: &[f64],
        removals: &[f64],
        column: usize,
    ) -> Result<(Vec<f64>, f64), DecayError> {
        let n = self.nodes;
        if !populations
            .iter()
            .chain(removals)
            .copied()
            .all(positive_or_zero)
        {
            return Err(DecayError::Unrepresentable(
                "nonnegative reaction maps or inventory",
            ));
        }
        let mut residuals = Vec::with_capacity(n);
        for i in 0..n {
            let mut residual = Sum::default();
            let mut scale = Sum::default();
            for value in [populations[i], -initial[i], removals[i]] {
                residual.add(value);
                scale.add(value.abs());
            }
            if let Some(injected) = injected {
                residual.add(-injected[i]);
                scale.add(injected[i]);
            }
            for j in 0..n {
                residual.product(-self.yields[i * n + j], removals[j]);
                scale.product(self.yields[i * n + j], removals[j]);
            }
            let value = residual.total();
            admit_balance(i, column, value, scale.total(), self.tolerance)?;
            residuals.push(value);
        }
        // This form retains conservation for zero-net-yield cycles without
        // subtracting two enormous, nearly equal cyclic event totals.
        let mut residual = Sum::default();
        let mut scale = Sum::default();
        for i in 0..n {
            residual.add(populations[i]);
            residual.add(-initial[i]);
            if let Some(injected) = injected {
                residual.add(-injected[i]);
                scale.add(injected[i]);
            }
            residual.product(-self.net_yields[i], removals[i]);
            scale.add(populations[i]);
            scale.add(initial[i]);
            scale.product(self.net_yields[i].abs(), removals[i]);
        }
        let value = residual.total();
        admit_balance(n, column, value, scale.total(), self.tolerance)?;
        Ok((residuals, value))
    }
}

fn admit_balance(
    row: usize,
    column: usize,
    residual: f64,
    scale: f64,
    tolerance: f64,
) -> Result<(), DecayError> {
    let allowance = tolerance * scale;
    if !residual.is_finite() || !scale.is_finite() || residual.abs() > allowance {
        Err(DecayError::Balance {
            row,
            column,
            residual,
            allowance,
        })
    } else {
        Ok(())
    }
}

fn restore_diagonals(
    plan: &mut ReactionPlan,
    net_rates: &[f64],
    cyclic: &[bool],
    time: f64,
) -> Result<(), DecayError> {
    for i in 0..plan.nodes {
        if cyclic[i] {
            continue;
        }
        let rate = net_rates[i];
        let exposure = rate * time;
        let population = det::exp(-exposure);
        let count = if plan.rates[i] == 0.0 {
            0.0
        } else if rate == 0.0 {
            plan.rates[i] * time
        } else if exposure.abs() < 0.5 {
            // exprel avoids the inverse of a tiny rate, including a net
            // rate almost cancelled by self-return products.
            let mean = if exposure == 0.0 {
                1.0
            } else {
                -det::expm1(-exposure) / exposure
            };
            (plan.rates[i] * time) * mean
        } else {
            (plan.rates[i] / rate) * -det::expm1(-exposure)
        };
        if !positive_or_zero(population) || !positive_or_zero(count) {
            return Err(DecayError::Unrepresentable(
                "scalar reaction survival or count",
            ));
        }
        plan.transition[i * plan.nodes + i] = population;
        plan.integrated[i * plan.nodes + i] = count;
    }
    Ok(())
}

fn prepare_tails(
    plan: &mut ReactionPlan,
    production: &[f64],
    net_rates: &[f64],
    cyclic: &[bool],
    limits: DecayLimits,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), DecayError> {
    let n = plan.nodes;
    for i in 0..n {
        if !cyclic[i] {
            plan.tails[i * n + i] = Tail {
                minimum_rate: net_rates[i],
                steps: 0,
                log_weight: 0.0,
            };
        }
    }
    if cyclic.iter().any(|&value| value) {
        return Ok(());
    }
    charge(&mut plan.work, n * n * (n + 1), limits.maximum_work)?;
    let mut order = Vec::with_capacity(n);
    let mut taken = vec![false; n];
    for _ in 0..n {
        poll(cancelled)?;
        let node = (0..n).find(|&i| {
            !taken[i] && (0..n).all(|j| j == i || production[i * n + j] == 0.0 || taken[j])
        });
        let Some(node) = node else {
            return Err(DecayError::Invalid("reaction DAG ordering"));
        };
        taken[node] = true;
        order.push(node);
    }
    for source in 0..n {
        poll(cancelled)?;
        let mut known = vec![false; n];
        known[source] = true;
        for &parent in &order {
            if !known[parent] {
                continue;
            }
            let tail = plan.tails[parent * n + source];
            for daughter in 0..n {
                let at = daughter * n + parent;
                if daughter == parent || production[at] == 0.0 {
                    continue;
                }
                let candidate = Tail {
                    minimum_rate: tail.minimum_rate.min(net_rates[daughter]),
                    steps: tail.steps + 1,
                    log_weight: if net_rates[parent] > 0.0 {
                        tail.log_weight + det::ln(production[at]) - det::ln(net_rates[parent])
                    } else {
                        f64::INFINITY
                    },
                };
                let target = &mut plan.tails[daughter * n + source];
                if known[daughter] {
                    target.minimum_rate = target.minimum_rate.min(candidate.minimum_rate);
                    target.steps = target.steps.max(candidate.steps);
                    // max(a,b)+ln(2) bounds log(exp(a)+exp(b)); an upper
                    // bound is sufficient and does not discard a tiny path.
                    target.log_weight =
                        target.log_weight.max(candidate.log_weight) + core::f64::consts::LN_2;
                } else {
                    *target = candidate;
                }
                known[daughter] = true;
            }
        }
    }
    Ok(())
}

#[cfg(test)]
#[path = "reaction/tests.rs"]
mod tests;
