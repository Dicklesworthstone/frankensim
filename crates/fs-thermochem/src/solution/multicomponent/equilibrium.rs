//! Finite-inventory, isothermal liquid phase-equilibrium candidates.
//!
//! The caller supplies a phase count and strictly positive seeds on the feed's
//! present-species support. A bounded damped Newton search transfers actual
//! species between phases, minimizing the extensive Gibbs energy of the SAME
//! regular-solution owner. Species-column allocations preserve the inventory;
//! the largest allocation is reconstructed, so a real trace is never obtained
//! by subtracting its rounded dominant component from one.
//!
//! Returned candidates pass chemical-potential equality, species accounting,
//! energy and complete local-curvature checks. This is a local, seeded search:
//! it neither discovers the correct phase count nor proves a global minimum.
//! In particular a locally stable one-phase candidate can be metastable. There
//! is no concentration floor, automatic phase insertion, phase disappearance,
//! interfacial energy, nucleation law or spatial transport in this API.

use super::stability::{
    MulticomponentLocalStabilityKindV1, MulticomponentLocalStabilityV1,
    MulticomponentStabilityControlV1, MulticomponentStabilityErrorV1,
};
use super::{
    MAX_REGULAR_SOLUTION_SPECIES_V1, MulticomponentRegularSolutionV1,
    MulticomponentSolutionErrorV1, MulticomponentSolutionStateV1, Sum, sum,
};
use crate::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;
use fs_math::det::sqrt;
use fs_qty::Temperature;

/// Compositions use the model's canonical species order. Weights are positive
/// relative phase amounts. They need not reproduce the feed: each species
/// column is normalized to its offered inventory before the search begins.
/// Every present feed species must be positive in every seed; feed absences
/// must be absent from every seed. No floor or missing pair fit is invented.
#[derive(Clone, Debug, PartialEq)]
pub struct MulticomponentEquilibriumSeedV1 {
    pub phase_mole_fractions: Vec<Vec<f64>>,
    pub phase_weights: Vec<f64>,
}

/// Bounded work and acceptance policy for one explicit phase count and seed.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MulticomponentEquilibriumControlV1 {
    pub max_phases: usize,
    pub max_iterations: usize,
    pub max_line_search_steps: usize,
    pub max_regularization_steps: usize,
    /// Peak owned numerical scratch, including both dense Newton matrices and
    /// the phase states. The local-stability scratch has a separate inner bound.
    pub max_workspace_bytes: usize,
    /// A smaller phase is refused, never deleted or redistributed silently.
    pub minimum_phase_fraction: f64,
    /// Required infinity-norm composition distance between distinct phases.
    pub minimum_phase_separation: f64,
    /// Maximum span of mu_i/(RT) over ALL phases for every present species.
    pub chemical_potential_tolerance: f64,
    /// Species-wise relative allowance, without an absolute trace allowance.
    pub relative_species_balance_tolerance: f64,
    /// One final allowance in G/(n_feed RT), against BOTH feed and seed energy.
    pub reduced_free_energy_tolerance: f64,
    pub stability: MulticomponentStabilityControlV1,
}
impl Default for MulticomponentEquilibriumControlV1 {
    fn default() -> Self {
        Self {
            max_phases: 8,
            max_iterations: 100,
            max_line_search_steps: 40,
            max_regularization_steps: 20,
            max_workspace_bytes: 16 * 1024 * 1024,
            minimum_phase_fraction: 1e-12,
            minimum_phase_separation: 1e-8,
            chemical_potential_tolerance: 1e-10,
            relative_species_balance_tolerance: 64.0 * f64::EPSILON,
            reduced_free_energy_tolerance: 2e-13,
            stability: MulticomponentStabilityControlV1::default(),
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub enum MulticomponentEquilibriumErrorV1 {
    Invalid(&'static str),
    WorkspaceBudget { required_bytes: usize },
    IterationLimit { completed_iterations: usize },
    LineSearchLimit { completed_iterations: usize },
    RegularizationLimit { completed_iterations: usize },
    PhaseDisappeared { phase: usize },
    CoincidentPhases { first: usize, second: usize },
    LocallyUnstablePhase { phase: usize },
    UnresolvedPhaseStability { phase: usize },
    Acceptance(&'static str),
    Unrepresentable(&'static str),
    Thermodynamics(MulticomponentSolutionErrorV1),
    Stability(MulticomponentStabilityErrorV1),
    Cancelled,
}
impl core::fmt::Display for MulticomponentEquilibriumErrorV1 {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "multicomponent liquid equilibrium candidate: {self:?}")
    }
}
impl std::error::Error for MulticomponentEquilibriumErrorV1 {}
impl From<MulticomponentSolutionErrorV1> for MulticomponentEquilibriumErrorV1 {
    fn from(error: MulticomponentSolutionErrorV1) -> Self {
        match error {
            MulticomponentSolutionErrorV1::Cancelled => Self::Cancelled,
            other => Self::Thermodynamics(other),
        }
    }
}
impl From<MulticomponentStabilityErrorV1> for MulticomponentEquilibriumErrorV1 {
    fn from(error: MulticomponentStabilityErrorV1) -> Self {
        match error {
            MulticomponentStabilityErrorV1::Cancelled => Self::Cancelled,
            other => Self::Stability(other),
        }
    }
}
type Result<T> = core::result::Result<T, MulticomponentEquilibriumErrorV1>;

/// Residuals and actually performed work; no global equilibrium certificate.
#[derive(Clone, Debug, PartialEq)]
pub struct MulticomponentEquilibriumReportV1 {
    pub maximum_reduced_chemical_potential_residual: f64,
    pub maximum_relative_species_balance_residual: f64,
    pub feed_mixing_gibbs_j: f64,
    pub seed_mixing_gibbs_j: f64,
    pub equilibrium_mixing_gibbs_j: f64,
    pub completed_iterations: usize,
    pub evaluated_states: usize,
    pub rejected_steps: usize,
    pub factorizations: usize,
    pub stability_sweeps: usize,
}

/// A phase with actual finite moles, the original thermodynamic owner, and its
/// full local-curvature assessment. Species order is unchanged from the feed.
#[derive(Clone, Debug, PartialEq)]
pub struct MulticomponentEquilibriumPhaseV1 {
    thermodynamics: MulticomponentSolutionStateV1,
    component_moles: Box<[f64]>,
    total_moles: f64,
    mole_fraction_of_feed: f64,
    local_stability: MulticomponentLocalStabilityV1,
}
impl MulticomponentEquilibriumPhaseV1 {
    pub const fn thermodynamics(&self) -> &MulticomponentSolutionStateV1 {
        &self.thermodynamics
    }
    pub fn component_moles(&self) -> &[f64] {
        &self.component_moles
    }
    pub const fn total_moles(&self) -> f64 {
        self.total_moles
    }
    pub const fn mole_fraction_of_feed(&self) -> f64 {
        self.mole_fraction_of_feed
    }
    pub const fn local_stability(&self) -> &MulticomponentLocalStabilityV1 {
        &self.local_stability
    }
}

/// Complete atomic candidate, retaining the original offered mole inventory.
/// Pure-liquid standard terms cancel because every species is conserved at the
/// same temperature and common caller-qualified pressure. Density, volume and
/// transport properties must be resolved separately for EACH resulting phase.
#[derive(Clone, Debug, PartialEq)]
pub struct MulticomponentEquilibriumCandidateV1 {
    feed: MulticomponentSolutionStateV1,
    component_moles: Box<[f64]>,
    total_moles: f64,
    phases: Box<[MulticomponentEquilibriumPhaseV1]>,
    controls: MulticomponentEquilibriumControlV1,
    report: MulticomponentEquilibriumReportV1,
}
impl MulticomponentEquilibriumCandidateV1 {
    pub const fn feed(&self) -> &MulticomponentSolutionStateV1 {
        &self.feed
    }
    pub fn component_moles(&self) -> &[f64] {
        &self.component_moles
    }
    pub const fn total_moles(&self) -> f64 {
        self.total_moles
    }
    pub fn phases(&self) -> &[MulticomponentEquilibriumPhaseV1] {
        &self.phases
    }
    pub const fn controls(&self) -> MulticomponentEquilibriumControlV1 {
        self.controls
    }
    pub const fn report(&self) -> &MulticomponentEquilibriumReportV1 {
        &self.report
    }
}

impl MulticomponentRegularSolutionV1 {
    /// Search one explicit, non-disappearing phase count. Phase count must not
    /// exceed the number of present feed species or the caller's bound. Dense
    /// damped Newton steps cost O(d^3), d=(phases-1)*present_species; scratch is
    /// preflighted before allocation. Trial states use this original law.
    ///
    /// Acceptance requires all present-species chemical potentials to match,
    /// each phase to be locally stable, distinct phase compositions, and no
    /// increase of extensive Gibbs energy versus either the homogeneous feed
    /// or the feasible initial seed, within ONE stated final energy allowance.
    /// Errors and cancellation expose no partially optimized phase inventory.
    pub fn equilibrium_candidate(
        &self,
        temperature: Temperature,
        component_moles: &[f64],
        seed: &MulticomponentEquilibriumSeedV1,
        controls: MulticomponentEquilibriumControlV1,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MulticomponentEquilibriumCandidateV1> {
        poll(&mut cancelled)?;
        validate_controls(controls)?;
        let n = self.species().len();
        let p = seed.phase_mole_fractions.len();
        if component_moles.len() != n || component_moles.iter().any(|v| !v.is_finite() || *v < 0.0)
        {
            return Err(MulticomponentEquilibriumErrorV1::Invalid("feed moles"));
        }
        let k = component_moles.iter().filter(|&&v| v > 0.0).count();
        if p == 0 || p > k || p > controls.max_phases {
            return Err(MulticomponentEquilibriumErrorV1::Invalid("phase count"));
        }
        if seed.phase_weights.len() != p
            || seed
                .phase_weights
                .iter()
                .any(|w| !w.is_finite() || *w <= 0.0)
            || seed.phase_mole_fractions.iter().any(|x| x.len() != n)
            || controls.minimum_phase_fraction * p as f64 >= 1.0
        {
            return Err(MulticomponentEquilibriumErrorV1::Invalid("phase seeds"));
        }
        let d = (p - 1) * k;
        // Include both current/trial states and the final phase objects; the
        // dense search scratch is released before the stability assessments.
        let base_bytes = (4 * p * k + 5 * d + 4 * n + 3 * k) * 8
            + 3 * p * core::mem::size_of::<MulticomponentSolutionStateV1>();
        let stability_bytes = (k * k + 2 * (k - 1) * (k - 1) + 10 * k + 2 * n) * 8;
        let required_bytes = base_bytes + (2 * d * d * 8).max(stability_bytes);
        if required_bytes > controls.max_workspace_bytes {
            return Err(MulticomponentEquilibriumErrorV1::WorkspaceBudget { required_bytes });
        }
        let total_moles = positive(sum(component_moles.iter().copied()), "total feed moles")?;
        let fractions = normalized_amounts(component_moles, total_moles)?;
        let feed = self.evaluate(temperature, &fractions, &mut cancelled)?;
        let rt = positive(R * temperature.value(), "RT")?;
        let active: Vec<_> = component_moles
            .iter()
            .enumerate()
            .filter_map(|(i, &m)| (m > 0.0).then_some(i))
            .collect();
        let z: Vec<_> = active.iter().map(|&i| feed.mole_fractions()[i]).collect();
        let roots: Vec<_> = z.iter().map(|&x| sqrt(x)).collect();
        let mut allocations = seed_allocations(self, temperature, seed, &active, &mut cancelled)?;
        let mut evaluated_states = 1 + p;
        let mut current = evaluate_allocations(
            self,
            temperature,
            &active,
            &z,
            &allocations,
            controls,
            rt,
            &mut cancelled,
        )?;
        evaluated_states += p;
        let seed_reduced_g = objective(&current);
        let mut completed = 0;
        let mut rejected = 0;
        let mut factorizations = 0;

        // Work arrays are reused for all iterations and backtracking trials.
        {
            let mut matrix = vec![0.0; d * d];
            let mut factor = vec![0.0; d * d];
            let mut gradient = vec![0.0; d];
            let mut direction = vec![0.0; d];
            let mut trial_allocations = vec![0.0; p * k];
            let mut variables = Vec::with_capacity(d);
            let mut anchors = vec![0; k];
            while chemical_residual(&current, &active) > controls.chemical_potential_tolerance {
                poll(&mut cancelled)?;
                if completed == controls.max_iterations {
                    return Err(MulticomponentEquilibriumErrorV1::IterationLimit {
                        completed_iterations: completed,
                    });
                }
                variables.clear();
                for i in 0..k {
                    anchors[i] = largest_allocation(&allocations, p, k, i);
                    for phase in 0..p {
                        if phase != anchors[i] {
                            variables.push(Variable {
                                component: i,
                                phase,
                            });
                        }
                    }
                }
                assemble_newton(
                    &current,
                    &active,
                    &roots,
                    &allocations,
                    &anchors,
                    &variables,
                    &mut matrix,
                    &mut gradient,
                    &mut cancelled,
                )?;
                newton_direction(
                    &matrix,
                    &gradient,
                    &mut factor,
                    &mut direction,
                    controls,
                    completed,
                    &mut factorizations,
                    &mut cancelled,
                )?;
                let slope = finite(
                    sum(gradient.iter().zip(&direction).map(|(g, v)| g * v)),
                    "Newton slope",
                )?;
                if slope >= 0.0 {
                    return Err(MulticomponentEquilibriumErrorV1::Unrepresentable(
                        "descent direction",
                    ));
                }
                // Use trial storage for allocation derivatives first. Every
                // species transfer is paired with its largest phase, protecting
                // tiny phase/species allocations from subtractive cancellation.
                trial_allocations.fill(0.0);
                for (v, &step) in variables.iter().zip(&direction) {
                    let delta = finite(step / roots[v.component], "allocation direction")?;
                    trial_allocations[v.phase * k + v.component] = delta;
                    trial_allocations[anchors[v.component] * k + v.component] -= delta;
                }
                let mut alpha = 1.0_f64;
                for (&a, &da) in allocations.iter().zip(&trial_allocations) {
                    if da < 0.0 {
                        alpha = alpha.min(0.99 * (-a / da));
                    }
                }
                positive(alpha, "positive line-search step")?;
                let old_g = objective(&current);
                let old_residual = chemical_residual(&current, &active);
                let mut accepted = None;
                let mut disappearing = None;
                for _ in 0..controls.max_line_search_steps {
                    poll(&mut cancelled)?;
                    trial_allocations.copy_from_slice(&allocations);
                    for (v, &step) in variables.iter().zip(&direction) {
                        let index = v.phase * k + v.component;
                        trial_allocations[index] =
                            allocations[index] + alpha * (step / roots[v.component]);
                    }
                    restore_anchors(&mut trial_allocations, p, k, &anchors)?;
                    let trial = match evaluate_allocations(
                        self,
                        temperature,
                        &active,
                        &z,
                        &trial_allocations,
                        controls,
                        rt,
                        &mut cancelled,
                    ) {
                        Ok(value) => value,
                        Err(MulticomponentEquilibriumErrorV1::PhaseDisappeared { phase }) => {
                            disappearing = Some(phase);
                            rejected += 1;
                            alpha *= 0.5;
                            continue;
                        }
                        Err(error) => return Err(error),
                    };
                    evaluated_states += p;
                    let trial_g = objective(&trial);
                    let trial_residual = chemical_residual(&trial, &active);
                    let roundoff = 32.0 * f64::EPSILON * old_g.abs().max(trial_g.abs()).max(1.0);
                    if trial_g <= old_g + 1e-4 * alpha * slope
                        || (trial_g <= old_g + roundoff && trial_residual < old_residual)
                    {
                        accepted = Some(trial);
                        break;
                    }
                    rejected += 1;
                    alpha *= 0.5;
                }
                current = match accepted {
                    Some(value) => value,
                    None => {
                        return Err(match disappearing {
                            Some(phase) => {
                                MulticomponentEquilibriumErrorV1::PhaseDisappeared { phase }
                            }
                            None => MulticomponentEquilibriumErrorV1::LineSearchLimit {
                                completed_iterations: completed,
                            },
                        });
                    }
                };
                allocations.copy_from_slice(&trial_allocations);
                completed += 1;
            }
        }

        // Rebuild each final state from the actual offered mole inventory,
        // rather than scaling an already rounded normalized inventory back up.
        let amounts = actual_phase_amounts(component_moles, &active, &allocations, p)?;
        let mut phases = Vec::with_capacity(p);
        let mut stability_sweeps = 0;
        for (phase, amounts) in amounts.into_iter().enumerate() {
            poll(&mut cancelled)?;
            let phase_moles = positive(sum(amounts.iter().copied()), "phase moles")?;
            let fraction = positive(phase_moles / total_moles, "phase fraction")?;
            if fraction < controls.minimum_phase_fraction {
                return Err(MulticomponentEquilibriumErrorV1::PhaseDisappeared { phase });
            }
            let x = normalized_amounts(&amounts, phase_moles)?;
            let thermodynamics = self.evaluate(temperature, &x, &mut cancelled)?;
            evaluated_states += 1;
            let local_stability =
                thermodynamics.local_stability(controls.stability, &mut cancelled)?;
            stability_sweeps += local_stability.completed_sweeps();
            match local_stability.kind() {
                MulticomponentLocalStabilityKindV1::Spinodal => {
                    return Err(MulticomponentEquilibriumErrorV1::LocallyUnstablePhase { phase });
                }
                MulticomponentLocalStabilityKindV1::Unresolved => {
                    return Err(MulticomponentEquilibriumErrorV1::UnresolvedPhaseStability {
                        phase,
                    });
                }
                _ => {}
            }
            phases.push(MulticomponentEquilibriumPhaseV1 {
                thermodynamics,
                component_moles: amounts.into_boxed_slice(),
                total_moles: phase_moles,
                mole_fraction_of_feed: fraction,
                local_stability,
            });
        }
        for first in 0..p {
            for second in first + 1..p {
                let distance = phases[first]
                    .thermodynamics
                    .mole_fractions()
                    .iter()
                    .zip(phases[second].thermodynamics.mole_fractions())
                    .map(|(a, b)| (a - b).abs())
                    .fold(0.0, f64::max);
                if distance <= controls.minimum_phase_separation {
                    return Err(MulticomponentEquilibriumErrorV1::CoincidentPhases {
                        first,
                        second,
                    });
                }
            }
        }
        let mut balance = 0.0_f64;
        for (i, &original) in component_moles.iter().enumerate() {
            let final_amount = sum(phases.iter().map(|phase| phase.component_moles[i]));
            if original == 0.0 {
                if final_amount != 0.0 {
                    return Err(MulticomponentEquilibriumErrorV1::Acceptance(
                        "absent species",
                    ));
                }
            } else {
                balance = balance.max((final_amount - original).abs() / original);
            }
        }
        let residual = phase_chemical_residual(&phases, &active);
        if balance > controls.relative_species_balance_tolerance {
            return Err(MulticomponentEquilibriumErrorV1::Acceptance(
                "species inventory",
            ));
        }
        if residual > controls.chemical_potential_tolerance {
            return Err(MulticomponentEquilibriumErrorV1::Acceptance(
                "chemical-potential equality",
            ));
        }
        let final_reduced_g = finite(
            sum(phases.iter().map(|phase| {
                phase.mole_fraction_of_feed * (phase.thermodynamics.mixing_gibbs_j_per_mol() / rt)
            })),
            "final reduced Gibbs energy",
        )?;
        let feed_reduced_g = feed.mixing_gibbs_j_per_mol() / rt;
        if final_reduced_g > feed_reduced_g + controls.reduced_free_energy_tolerance
            || final_reduced_g > seed_reduced_g + controls.reduced_free_energy_tolerance
        {
            return Err(MulticomponentEquilibriumErrorV1::Acceptance(
                "free energy increased",
            ));
        }
        let report = MulticomponentEquilibriumReportV1 {
            maximum_reduced_chemical_potential_residual: residual,
            maximum_relative_species_balance_residual: balance,
            feed_mixing_gibbs_j: finite(
                feed.mixing_gibbs_j_per_mol() * total_moles,
                "feed Gibbs energy",
            )?,
            seed_mixing_gibbs_j: finite((seed_reduced_g * rt) * total_moles, "seed Gibbs energy")?,
            equilibrium_mixing_gibbs_j: finite(
                sum(phases.iter().map(|phase| {
                    phase.thermodynamics.mixing_gibbs_j_per_mol() * phase.total_moles
                })),
                "equilibrium Gibbs energy",
            )?,
            completed_iterations: completed,
            evaluated_states,
            rejected_steps: rejected,
            factorizations,
            stability_sweeps,
        };
        poll(&mut cancelled)?;
        Ok(MulticomponentEquilibriumCandidateV1 {
            feed,
            component_moles: component_moles.into(),
            total_moles,
            phases: phases.into_boxed_slice(),
            controls,
            report,
        })
    }
}

fn validate_controls(c: MulticomponentEquilibriumControlV1) -> Result<()> {
    if !(1..=MAX_REGULAR_SOLUTION_SPECIES_V1).contains(&c.max_phases)
        || c.max_iterations > 10_000
        || !(1..=256).contains(&c.max_line_search_steps)
        || !(1..=64).contains(&c.max_regularization_steps)
        || !c.minimum_phase_fraction.is_finite()
        || c.minimum_phase_fraction <= 0.0
        || c.minimum_phase_fraction >= 1.0
        || !c.minimum_phase_separation.is_finite()
        || !(0.0..1.0).contains(&c.minimum_phase_separation)
        || !c.chemical_potential_tolerance.is_finite()
        || !(4.0 * f64::EPSILON..=1e-3).contains(&c.chemical_potential_tolerance)
        || !c.relative_species_balance_tolerance.is_finite()
        || !(0.0..=1e-8).contains(&c.relative_species_balance_tolerance)
        || !c.reduced_free_energy_tolerance.is_finite()
        || !(0.0..=1e-6).contains(&c.reduced_free_energy_tolerance)
        || c.stability.max_sweeps == 0
        || c.stability.max_sweeps > 256
        || !c.stability.relative_eigensolve_tolerance.is_finite()
        || c.stability.relative_eigensolve_tolerance <= 0.0
        || c.stability.relative_eigensolve_tolerance >= 1.0
        || !c.stability.stability_tolerance.is_finite()
        || c.stability.stability_tolerance < 0.0
    {
        return Err(MulticomponentEquilibriumErrorV1::Invalid(
            "equilibrium controls",
        ));
    }
    Ok(())
}

fn normalized_amounts(amounts: &[f64], total: f64) -> Result<Vec<f64>> {
    amounts
        .iter()
        .map(|&amount| {
            let x = finite(amount / total, "mole fraction")?;
            if amount > 0.0 && x == 0.0 {
                return Err(MulticomponentEquilibriumErrorV1::Unrepresentable(
                    "positive mole fraction",
                ));
            }
            Ok(x)
        })
        .collect()
}

fn seed_allocations(
    model: &MulticomponentRegularSolutionV1,
    temperature: Temperature,
    seed: &MulticomponentEquilibriumSeedV1,
    active: &[usize],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<f64>> {
    let p = seed.phase_weights.len();
    let k = active.len();
    let mut allocations = vec![0.0; p * k];
    let weight_scale = seed.phase_weights.iter().copied().fold(0.0, f64::max);
    for phase in 0..p {
        let state = model.evaluate(
            temperature,
            &seed.phase_mole_fractions[phase],
            &mut *cancelled,
        )?;
        for (i, &x) in state.mole_fractions().iter().enumerate() {
            if active.binary_search(&i).is_ok() != (x > 0.0) {
                return Err(MulticomponentEquilibriumErrorV1::Invalid(
                    "seed must have exactly the feed support",
                ));
            }
        }
        for (i, &species) in active.iter().enumerate() {
            allocations[phase * k + i] = positive(
                (seed.phase_weights[phase] / weight_scale) * state.mole_fractions()[species],
                "positive seed allocation",
            )?;
        }
    }
    let mut anchors = vec![0; k];
    for i in 0..k {
        let column = positive(
            sum((0..p).map(|phase| allocations[phase * k + i])),
            "seed species column",
        )?;
        for phase in 0..p {
            allocations[phase * k + i] =
                positive(allocations[phase * k + i] / column, "normalized allocation")?;
        }
        anchors[i] = largest_allocation(&allocations, p, k, i);
    }
    restore_anchors(&mut allocations, p, k, &anchors)?;
    Ok(allocations)
}

fn largest_allocation(allocations: &[f64], p: usize, k: usize, i: usize) -> usize {
    let mut largest = 0;
    for phase in 1..p {
        if allocations[phase * k + i] > allocations[largest * k + i] {
            largest = phase;
        }
    }
    largest
}

fn restore_anchors(allocations: &mut [f64], p: usize, k: usize, anchors: &[usize]) -> Result<()> {
    for i in 0..k {
        let others = sum((0..p)
            .filter(|&phase| phase != anchors[i])
            .map(|phase| allocations[phase * k + i]));
        allocations[anchors[i] * k + i] = positive(1.0 - others, "conserved largest allocation")?;
    }
    if allocations
        .iter()
        .any(|a| !a.is_finite() || *a <= 0.0 || *a > 1.0)
    {
        return Err(MulticomponentEquilibriumErrorV1::Unrepresentable(
            "positive phase allocations",
        ));
    }
    Ok(())
}

struct PhaseEvaluation {
    state: MulticomponentSolutionStateV1,
    fraction: f64,
    reduced_g: f64,
}

#[allow(clippy::too_many_arguments)]
fn evaluate_allocations(
    model: &MulticomponentRegularSolutionV1,
    temperature: Temperature,
    active: &[usize],
    z: &[f64],
    allocations: &[f64],
    controls: MulticomponentEquilibriumControlV1,
    rt: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<PhaseEvaluation>> {
    let k = active.len();
    let mut phases = Vec::with_capacity(allocations.len() / k);
    let mut x = [0.0; MAX_REGULAR_SOLUTION_SPECIES_V1];
    for (phase, column) in allocations.chunks_exact(k).enumerate() {
        poll(cancelled)?;
        let fraction = positive(
            sum(z.iter().zip(column).map(|(z, a)| z * a)),
            "normalized phase amount",
        )?;
        if fraction < controls.minimum_phase_fraction {
            return Err(MulticomponentEquilibriumErrorV1::PhaseDisappeared { phase });
        }
        for (i, &species) in active.iter().enumerate() {
            x[species] = positive((z[i] * column[i]) / fraction, "positive phase composition")?;
        }
        let state = model.evaluate(temperature, &x[..model.species().len()], &mut *cancelled)?;
        let reduced_g = finite(state.mixing_gibbs_j_per_mol() / rt, "reduced Gibbs energy")?;
        phases.push(PhaseEvaluation {
            state,
            fraction,
            reduced_g,
        });
    }
    Ok(phases)
}

fn objective(phases: &[PhaseEvaluation]) -> f64 {
    sum(phases.iter().map(|p| p.fraction * p.reduced_g))
}

fn chemical_residual(phases: &[PhaseEvaluation], active: &[usize]) -> f64 {
    chemical_span(active, phases.iter().map(|p| &p.state))
}

fn phase_chemical_residual(phases: &[MulticomponentEquilibriumPhaseV1], active: &[usize]) -> f64 {
    chemical_span(active, phases.iter().map(|p| &p.thermodynamics))
}

fn chemical_span<'a>(
    active: &[usize],
    states: impl Iterator<Item = &'a MulticomponentSolutionStateV1> + Clone,
) -> f64 {
    let mut maximum = 0.0_f64;
    for &i in active {
        let mut low = f64::INFINITY;
        let mut high = f64::NEG_INFINITY;
        for state in states.clone() {
            let mu = state.reduced_chemical_potentials()[i].expect("present species");
            low = low.min(mu);
            high = high.max(mu);
        }
        maximum = maximum.max(high - low);
    }
    maximum
}

#[derive(Clone, Copy)]
struct Variable {
    component: usize,
    phase: usize,
}

#[allow(clippy::too_many_arguments)]
fn assemble_newton(
    phases: &[PhaseEvaluation],
    active: &[usize],
    roots: &[f64],
    allocations: &[f64],
    anchors: &[usize],
    variables: &[Variable],
    matrix: &mut [f64],
    gradient: &mut [f64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<()> {
    let d = variables.len();
    let k = active.len();
    for (row, v) in variables.iter().enumerate() {
        poll(cancelled)?;
        let i = active[v.component];
        let a = anchors[v.component];
        gradient[row] = finite(
            roots[v.component]
                * (phases[v.phase].state.reduced_chemical_potentials()[i]
                    .expect("present species")
                    - phases[a].state.reduced_chemical_potentials()[i].expect("present species")),
            "scaled chemical-potential gradient",
        )?;
        for (column, w) in variables.iter().enumerate().take(row + 1) {
            let j = active[w.component];
            let b = anchors[w.component];
            let mut value = Sum::default();
            for (pv, sv) in [(v.phase, 1.0), (a, -1.0)] {
                for (pw, sw) in [(w.phase, 1.0), (b, -1.0)] {
                    if pv != pw {
                        continue;
                    }
                    let phase = &phases[pv];
                    // Entropy-scaled extensive Hessian: the diagonal ideal
                    // term is 1/allocation, not 1/x or 1/(trace mole amount).
                    // gamma_i = (chi*x)_i - g_ex eliminates a cancellation.
                    let ideal = if v.component == w.component {
                        1.0 / allocations[pv * k + v.component]
                    } else {
                        0.0
                    };
                    let response = phase.state.chi(i, j)
                        - phase.state.log_activity_coefficients()[i]
                        - phase.state.log_activity_coefficients()[j]
                        - 1.0;
                    let h = ideal
                        + ((roots[v.component] * roots[w.component]) / phase.fraction) * response;
                    value.add((sv * sw) * h);
                }
            }
            let h = finite(value.value(), "scaled extensive Gibbs Hessian")?;
            matrix[row * d + column] = h;
            matrix[column * d + row] = h;
        }
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn newton_direction(
    matrix: &[f64],
    gradient: &[f64],
    factor: &mut [f64],
    direction: &mut [f64],
    controls: MulticomponentEquilibriumControlV1,
    completed: usize,
    factorizations: &mut usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<()> {
    let d = gradient.len();
    let scale = matrix.iter().map(|v| v.abs()).fold(1.0, f64::max);
    let mut shift = 0.0;
    for attempt in 0..controls.max_regularization_steps {
        poll(cancelled)?;
        *factorizations += 1;
        factor.copy_from_slice(matrix);
        let mut success = true;
        for row in 0..d {
            poll(cancelled)?;
            for column in 0..=row {
                let dot = sum((0..column).map(|j| factor[row * d + j] * factor[column * d + j]));
                let mut value = factor[row * d + column] - dot;
                if row == column {
                    value += shift;
                    if !value.is_finite() || value <= 16.0 * f64::EPSILON * scale {
                        success = false;
                        break;
                    }
                    factor[row * d + column] = sqrt(value);
                } else {
                    factor[row * d + column] = value / factor[column * d + column];
                }
            }
            if !success {
                break;
            }
        }
        if success {
            for row in 0..d {
                let dot = sum((0..row).map(|j| factor[row * d + j] * direction[j]));
                direction[row] = finite(
                    (-gradient[row] - dot) / factor[row * d + row],
                    "forward Newton solve",
                )?;
            }
            for row in (0..d).rev() {
                let dot = sum((row + 1..d).map(|j| factor[j * d + row] * direction[j]));
                direction[row] = finite(
                    (direction[row] - dot) / factor[row * d + row],
                    "backward Newton solve",
                )?;
            }
            return Ok(());
        }
        shift = if attempt == 0 {
            1e-8 * scale
        } else {
            shift * 10.0
        };
    }
    Err(MulticomponentEquilibriumErrorV1::RegularizationLimit {
        completed_iterations: completed,
    })
}

fn actual_phase_amounts(
    feed: &[f64],
    active: &[usize],
    allocations: &[f64],
    p: usize,
) -> Result<Vec<Vec<f64>>> {
    let k = active.len();
    let mut amounts = vec![vec![0.0; feed.len()]; p];
    for (i, &species) in active.iter().enumerate() {
        let anchor = largest_allocation(allocations, p, k, i);
        for phase in 0..p {
            if phase != anchor {
                amounts[phase][species] = positive(
                    feed[species] * allocations[phase * k + i],
                    "actual phase species moles",
                )?;
            }
        }
        let others = sum((0..p)
            .filter(|&phase| phase != anchor)
            .map(|phase| amounts[phase][species]));
        amounts[anchor][species] =
            positive(feed[species] - others, "conserved phase species moles")?;
    }
    Ok(amounts)
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<()> {
    if cancelled() {
        Err(MulticomponentEquilibriumErrorV1::Cancelled)
    } else {
        Ok(())
    }
}
fn finite(value: f64, label: &'static str) -> Result<f64> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(MulticomponentEquilibriumErrorV1::Unrepresentable(label))
    }
}
fn positive(value: f64, label: &'static str) -> Result<f64> {
    if value.is_finite() && value > 0.0 {
        Ok(value)
    } else {
        Err(MulticomponentEquilibriumErrorV1::Unrepresentable(label))
    }
}

#[cfg(test)]
mod tests;
