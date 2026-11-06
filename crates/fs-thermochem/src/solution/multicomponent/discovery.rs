//! Bounded deterministic phase-count and seed discovery for multicomponent liquids.
//!
//! The explicit equilibrium owner remains the authority for every accepted phase
//! candidate. This layer only proposes a bounded sequence of physically admissible
//! seeds, runs that existing solver, and keeps the lowest accepted extensive Gibbs
//! energy. It never promotes finite multistart search into a global-minimum proof.
//!
//! Seeds preserve the feed support exactly. Species-rich starts are convex blends
//! of the actual feed with a canonical species vertex, so every present trace stays
//! positive in every phase without an absolute concentration floor. A spinodal
//! feed also gets a two-phase seed along its complete smallest-curvature mode.
//! Exact feed absences remain absent in all generated seeds.

use super::equilibrium::{
    MulticomponentEquilibriumCandidateV1, MulticomponentEquilibriumControlV1,
    MulticomponentEquilibriumErrorV1, MulticomponentEquilibriumSeedV1,
};
use super::stability::{
    MulticomponentLocalStabilityKindV1, MulticomponentLocalStabilityV1,
    MulticomponentStabilityErrorV1,
};
use super::{
    MulticomponentRegularSolutionV1, MulticomponentSolutionErrorV1,
    MulticomponentSolutionStateV1, sum,
};
use crate::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;
use fs_qty::Temperature;

/// Search policy. Candidate solves are sequential so peak scratch is still owned
/// by the existing equilibrium control rather than multiplied by seed count.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MulticomponentPhaseDiscoveryControlV1 {
    /// Maximum phase count considered, including the homogeneous candidate.
    pub maximum_phases: usize,
    /// Whole-call bound on candidate solves, including the homogeneous attempt.
    pub maximum_seed_trials: usize,
    /// Three deterministic convex enrichments used for canonical species anchors.
    pub seed_enrichments: [f64; 3],
    /// Fraction of the exact positivity-limited step used for a spinodal mode.
    pub spinodal_step_fraction: f64,
    /// Reduced-energy tie allowance. Near ties prefer fewer phases.
    pub reduced_selection_tolerance: f64,
    pub equilibrium: MulticomponentEquilibriumControlV1,
}
impl Default for MulticomponentPhaseDiscoveryControlV1 {
    fn default() -> Self {
        Self {
            maximum_phases: 4,
            maximum_seed_trials: 64,
            seed_enrichments: [0.85, 0.60, 0.35],
            spinodal_step_fraction: 0.75,
            reduced_selection_tolerance: 2e-13,
            equilibrium: MulticomponentEquilibriumControlV1::default(),
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub enum MulticomponentPhaseDiscoveryErrorV1 {
    Invalid(&'static str),
    Thermodynamics(MulticomponentSolutionErrorV1),
    Stability(MulticomponentStabilityErrorV1),
    NoAdmissibleCandidate { attempted_candidates: usize },
    Unrepresentable(&'static str),
    Cancelled,
}
impl core::fmt::Display for MulticomponentPhaseDiscoveryErrorV1 {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "multicomponent phase discovery: {self:?}")
    }
}
impl std::error::Error for MulticomponentPhaseDiscoveryErrorV1 {}
impl From<MulticomponentSolutionErrorV1> for MulticomponentPhaseDiscoveryErrorV1 {
    fn from(error: MulticomponentSolutionErrorV1) -> Self {
        match error {
            MulticomponentSolutionErrorV1::Cancelled => Self::Cancelled,
            other => Self::Thermodynamics(other),
        }
    }
}
impl From<MulticomponentStabilityErrorV1> for MulticomponentPhaseDiscoveryErrorV1 {
    fn from(error: MulticomponentStabilityErrorV1) -> Self {
        match error {
            MulticomponentStabilityErrorV1::Cancelled => Self::Cancelled,
            other => Self::Stability(other),
        }
    }
}
type Result<T> = core::result::Result<T, MulticomponentPhaseDiscoveryErrorV1>;

/// Bounded-search evidence. Counts are actual work, not projected maxima.
#[derive(Clone, Debug, PartialEq)]
pub struct MulticomponentPhaseDiscoveryReportV1 {
    feed_local_stability: MulticomponentLocalStabilityV1,
    attempted_candidates: usize,
    converged_candidates: usize,
    rejected_candidates: usize,
    selected_phase_count: usize,
    maximum_phase_count_attempted: usize,
    search_exhausted: bool,
    most_negative_seed_tangent_plane_distance: Option<f64>,
    selected_reduced_gibbs_improvement: f64,
}
impl MulticomponentPhaseDiscoveryReportV1 {
    pub const fn feed_local_stability(&self) -> &MulticomponentLocalStabilityV1 {
        &self.feed_local_stability
    }
    pub const fn attempted_candidates(&self) -> usize {
        self.attempted_candidates
    }
    pub const fn converged_candidates(&self) -> usize {
        self.converged_candidates
    }
    pub const fn rejected_candidates(&self) -> usize {
        self.rejected_candidates
    }
    pub const fn selected_phase_count(&self) -> usize {
        self.selected_phase_count
    }
    pub const fn maximum_phase_count_attempted(&self) -> usize {
        self.maximum_phase_count_attempted
    }
    pub const fn search_exhausted(&self) -> bool {
        self.search_exhausted
    }
    /// Minimum finite tangent-plane distance among generated nonhomogeneous seeds.
    /// A negative value is an instability witness; a nonnegative value is not a
    /// global stability certificate because only finitely many trials were made.
    pub const fn most_negative_seed_tangent_plane_distance(&self) -> Option<f64> {
        self.most_negative_seed_tangent_plane_distance
    }
    /// (G_feed-G_selected)/(n_feed R T). Positive means the selected candidate
    /// lowered the conserved feed's mixing Gibbs energy.
    pub const fn selected_reduced_gibbs_improvement(&self) -> f64 {
        self.selected_reduced_gibbs_improvement
    }
}

/// Best fully admitted candidate found by the declared bounded deterministic search.
#[derive(Clone, Debug, PartialEq)]
pub struct MulticomponentDiscoveredEquilibriumV1 {
    candidate: MulticomponentEquilibriumCandidateV1,
    controls: MulticomponentPhaseDiscoveryControlV1,
    report: MulticomponentPhaseDiscoveryReportV1,
}
impl MulticomponentDiscoveredEquilibriumV1 {
    pub const fn candidate(&self) -> &MulticomponentEquilibriumCandidateV1 {
        &self.candidate
    }
    pub const fn controls(&self) -> MulticomponentPhaseDiscoveryControlV1 {
        self.controls
    }
    pub const fn report(&self) -> &MulticomponentPhaseDiscoveryReportV1 {
        &self.report
    }
    pub fn into_candidate(self) -> MulticomponentEquilibriumCandidateV1 {
        self.candidate
    }
}

impl MulticomponentRegularSolutionV1 {
    /// Discover the best accepted local phase candidate among a bounded,
    /// deterministic family of phase counts and seeds.
    ///
    /// The search always attempts the homogeneous candidate first. If the feed
    /// is spinodally unstable, one two-phase seed follows the complete smallest
    /// curvature direction. Remaining starts use lexicographic combinations of
    /// present species as rich-phase anchors at each declared enrichment.
    ///
    /// Every returned candidate has already passed the explicit equilibrium
    /// solver's species conservation, chemical-potential equality, local phase
    /// stability, phase separation, and free-energy checks. Search exhaustion,
    /// however, is not evidence that an untried lower-Gibbs candidate does not
    /// exist.
    pub fn discover_equilibrium_candidate(
        &self,
        temperature: Temperature,
        component_moles: &[f64],
        controls: MulticomponentPhaseDiscoveryControlV1,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MulticomponentDiscoveredEquilibriumV1> {
        poll(&mut cancelled)?;
        validate_controls(controls)?;
        if component_moles.len() != self.species().len()
            || component_moles
                .iter()
                .any(|value| !value.is_finite() || *value < 0.0)
        {
            return Err(MulticomponentPhaseDiscoveryErrorV1::Invalid(
                "feed moles",
            ));
        }
        let total_moles = positive(sum(component_moles.iter().copied()), "total feed moles")?;
        let fractions = normalized(component_moles, total_moles)?;
        let feed = self.evaluate(temperature, &fractions, &mut cancelled)?;
        let active: Vec<_> = component_moles
            .iter()
            .enumerate()
            .filter_map(|(i, &value)| (value > 0.0).then_some(i))
            .collect();
        if active.is_empty() {
            return Err(MulticomponentPhaseDiscoveryErrorV1::Invalid(
                "positive feed inventory",
            ));
        }
        let feed_stability = feed.local_stability(controls.equilibrium.stability, &mut cancelled)?;
        let maximum_phases = controls
            .maximum_phases
            .min(active.len())
            .min(controls.equilibrium.max_phases);
        let rt = positive(R * temperature.value(), "RT")?;
        let energy_tolerance = finite(
            controls.reduced_selection_tolerance * total_moles * rt,
            "selection energy tolerance",
        )?;

        let mut attempted = 0usize;
        let mut converged = 0usize;
        let mut rejected = 0usize;
        let mut maximum_attempted = 0usize;
        let mut best: Option<MulticomponentEquilibriumCandidateV1> = None;
        let mut best_energy = f64::INFINITY;
        let mut minimum_tpd: Option<f64> = None;
        let mut search_exhausted = false;

        let homogeneous = MulticomponentEquilibriumSeedV1 {
            phase_mole_fractions: vec![feed.mole_fractions().to_vec()],
            phase_weights: vec![1.0],
        };
        attempt(
            self,
            temperature,
            component_moles,
            &homogeneous,
            controls,
            energy_tolerance,
            &mut attempted,
            &mut converged,
            &mut rejected,
            &mut maximum_attempted,
            &mut best,
            &mut best_energy,
            &mut cancelled,
        )?;

        if attempted < controls.maximum_seed_trials
            && maximum_phases >= 2
            && feed_stability.kind() == MulticomponentLocalStabilityKindV1::Spinodal
        {
            if let Some(seed) = spinodal_seed(
                &feed,
                feed_stability.direction(),
                controls.spinodal_step_fraction,
            )? {
                update_tpd(&feed, &seed, &mut minimum_tpd, &mut cancelled)?;
                attempt(
                    self,
                    temperature,
                    component_moles,
                    &seed,
                    controls,
                    energy_tolerance,
                    &mut attempted,
                    &mut converged,
                    &mut rejected,
                    &mut maximum_attempted,
                    &mut best,
                    &mut best_energy,
                    &mut cancelled,
                )?;
            }
        }

        'phase_counts: for phase_count in 2..=maximum_phases {
            let mut combination: Vec<usize> = (0..phase_count).collect();
            loop {
                for enrichment in controls.seed_enrichments {
                    if attempted >= controls.maximum_seed_trials {
                        search_exhausted = true;
                        break 'phase_counts;
                    }
                    let anchors: Vec<_> = combination.iter().map(|&i| active[i]).collect();
                    let seed = anchored_seed(&feed, &anchors, enrichment)?;
                    update_tpd(&feed, &seed, &mut minimum_tpd, &mut cancelled)?;
                    attempt(
                        self,
                        temperature,
                        component_moles,
                        &seed,
                        controls,
                        energy_tolerance,
                        &mut attempted,
                        &mut converged,
                        &mut rejected,
                        &mut maximum_attempted,
                        &mut best,
                        &mut best_energy,
                        &mut cancelled,
                    )?;
                }
                if !next_combination(&mut combination, active.len()) {
                    break;
                }
            }
        }

        poll(&mut cancelled)?;
        let candidate = best.ok_or(
            MulticomponentPhaseDiscoveryErrorV1::NoAdmissibleCandidate {
                attempted_candidates: attempted,
            },
        )?;
        let selected_phase_count = candidate.phases().len();
        let selected_reduced_gibbs_improvement = finite(
            (candidate.report().feed_mixing_gibbs_j
                - candidate.report().equilibrium_mixing_gibbs_j)
                / (total_moles * rt),
            "selected reduced Gibbs improvement",
        )?;
        Ok(MulticomponentDiscoveredEquilibriumV1 {
            candidate,
            controls,
            report: MulticomponentPhaseDiscoveryReportV1 {
                feed_local_stability: feed_stability,
                attempted_candidates: attempted,
                converged_candidates: converged,
                rejected_candidates: rejected,
                selected_phase_count,
                maximum_phase_count_attempted: maximum_attempted,
                search_exhausted,
                most_negative_seed_tangent_plane_distance: minimum_tpd,
                selected_reduced_gibbs_improvement,
            },
        })
    }
}

#[allow(clippy::too_many_arguments)]
fn attempt(
    model: &MulticomponentRegularSolutionV1,
    temperature: Temperature,
    component_moles: &[f64],
    seed: &MulticomponentEquilibriumSeedV1,
    controls: MulticomponentPhaseDiscoveryControlV1,
    energy_tolerance: f64,
    attempted: &mut usize,
    converged: &mut usize,
    rejected: &mut usize,
    maximum_attempted: &mut usize,
    best: &mut Option<MulticomponentEquilibriumCandidateV1>,
    best_energy: &mut f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<()> {
    if *attempted >= controls.maximum_seed_trials {
        return Ok(());
    }
    *attempted += 1;
    *maximum_attempted = (*maximum_attempted).max(seed.phase_mole_fractions.len());
    let candidate = match model.equilibrium_candidate(
        temperature,
        component_moles,
        seed,
        controls.equilibrium,
        &mut *cancelled,
    ) {
        Ok(candidate) => candidate,
        Err(MulticomponentEquilibriumErrorV1::Cancelled) => {
            return Err(MulticomponentPhaseDiscoveryErrorV1::Cancelled);
        }
        Err(_) => {
            *rejected += 1;
            return Ok(());
        }
    };
    *converged += 1;
    let energy = candidate.report().equilibrium_mixing_gibbs_j;
    if !energy.is_finite() {
        return Err(MulticomponentPhaseDiscoveryErrorV1::Unrepresentable(
            "candidate Gibbs energy",
        ));
    }
    let better = if best.is_none() {
        true
    } else if energy < *best_energy - energy_tolerance {
        true
    } else {
        (energy - *best_energy).abs() <= energy_tolerance
            && best.as_ref().is_some_and(|current| candidate.phases().len() < current.phases().len())
    };
    if better {
        *best_energy = energy;
        *best = Some(candidate);
    }
    Ok(())
}

fn anchored_seed(
    feed: &MulticomponentSolutionStateV1,
    anchors: &[usize],
    enrichment: f64,
) -> Result<MulticomponentEquilibriumSeedV1> {
    let residual = 1.0 - enrichment;
    let mut phases = Vec::with_capacity(anchors.len());
    for &anchor in anchors {
        let mut x: Vec<_> = feed
            .mole_fractions()
            .iter()
            .map(|&value| value * residual)
            .collect();
        x[anchor] = finite(x[anchor] + enrichment, "enriched seed composition")?;
        if x.iter().any(|value| !value.is_finite() || *value < 0.0)
            || feed
                .mole_fractions()
                .iter()
                .zip(&x)
                .any(|(&original, &trial)| (original > 0.0) != (trial > 0.0))
        {
            return Err(MulticomponentPhaseDiscoveryErrorV1::Unrepresentable(
                "seed support",
            ));
        }
        phases.push(x);
    }
    Ok(MulticomponentEquilibriumSeedV1 {
        phase_mole_fractions: phases,
        phase_weights: vec![1.0; anchors.len()],
    })
}

fn spinodal_seed(
    feed: &MulticomponentSolutionStateV1,
    direction: &[f64],
    fraction: f64,
) -> Result<Option<MulticomponentEquilibriumSeedV1>> {
    let x = feed.mole_fractions();
    if direction.len() != x.len() {
        return Err(MulticomponentPhaseDiscoveryErrorV1::Invalid(
            "stability direction dimensions",
        ));
    }
    let mut plus_limit = f64::INFINITY;
    let mut minus_limit = f64::INFINITY;
    for (&value, &delta) in x.iter().zip(direction) {
        if delta < 0.0 {
            plus_limit = plus_limit.min(value / -delta);
        } else if delta > 0.0 {
            minus_limit = minus_limit.min(value / delta);
        }
    }
    let limit = plus_limit.min(minus_limit);
    if !limit.is_finite() || limit <= 0.0 {
        return Ok(None);
    }
    let step = positive(limit * fraction, "spinodal seed step")?;
    let mut plus = Vec::with_capacity(x.len());
    let mut minus = Vec::with_capacity(x.len());
    for (&value, &delta) in x.iter().zip(direction) {
        let upper = finite(value + step * delta, "spinodal positive seed")?;
        let lower = finite(value - step * delta, "spinodal negative seed")?;
        if (value > 0.0 && (upper <= 0.0 || lower <= 0.0))
            || (value == 0.0 && (upper != 0.0 || lower != 0.0))
        {
            return Ok(None);
        }
        plus.push(upper);
        minus.push(lower);
    }
    Ok(Some(MulticomponentEquilibriumSeedV1 {
        phase_mole_fractions: vec![plus, minus],
        phase_weights: vec![1.0, 1.0],
    }))
}

fn update_tpd(
    feed: &MulticomponentSolutionStateV1,
    seed: &MulticomponentEquilibriumSeedV1,
    minimum: &mut Option<f64>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<()> {
    for trial in &seed.phase_mole_fractions {
        poll(cancelled)?;
        let value = feed.tangent_plane_distance(trial, &mut *cancelled)?;
        if !value.is_finite() {
            return Err(MulticomponentPhaseDiscoveryErrorV1::Unrepresentable(
                "seed tangent-plane distance",
            ));
        }
        *minimum = Some((*minimum).map_or(value, |current| current.min(value)));
    }
    Ok(())
}

fn next_combination(combination: &mut [usize], n: usize) -> bool {
    let k = combination.len();
    for i in (0..k).rev() {
        let maximum = n - k + i;
        if combination[i] < maximum {
            combination[i] += 1;
            for j in i + 1..k {
                combination[j] = combination[j - 1] + 1;
            }
            return true;
        }
    }
    false
}

fn normalized(amounts: &[f64], total: f64) -> Result<Vec<f64>> {
    amounts
        .iter()
        .map(|&amount| {
            let x = finite(amount / total, "feed mole fraction")?;
            if amount > 0.0 && x == 0.0 {
                return Err(MulticomponentPhaseDiscoveryErrorV1::Unrepresentable(
                    "positive feed trace",
                ));
            }
            Ok(x)
        })
        .collect()
}

fn validate_controls(control: MulticomponentPhaseDiscoveryControlV1) -> Result<()> {
    if control.maximum_phases == 0
        || control.maximum_phases > control.equilibrium.max_phases
        || control.maximum_seed_trials == 0
        || control.maximum_seed_trials > 4096
        || control
            .seed_enrichments
            .iter()
            .any(|value| !value.is_finite() || *value <= 0.0 || *value >= 1.0)
        || !control.spinodal_step_fraction.is_finite()
        || control.spinodal_step_fraction <= 0.0
        || control.spinodal_step_fraction >= 1.0
        || !control.reduced_selection_tolerance.is_finite()
        || control.reduced_selection_tolerance < 0.0
        || control.reduced_selection_tolerance > 1e-6
    {
        return Err(MulticomponentPhaseDiscoveryErrorV1::Invalid(
            "discovery controls",
        ));
    }
    Ok(())
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<()> {
    if cancelled() {
        Err(MulticomponentPhaseDiscoveryErrorV1::Cancelled)
    } else {
        Ok(())
    }
}
fn finite(value: f64, label: &'static str) -> Result<f64> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(MulticomponentPhaseDiscoveryErrorV1::Unrepresentable(label))
    }
}
fn positive(value: f64, label: &'static str) -> Result<f64> {
    if value.is_finite() && value > 0.0 {
        Ok(value)
    } else {
        Err(MulticomponentPhaseDiscoveryErrorV1::Unrepresentable(label))
    }
}

#[cfg(test)]
mod tests;
