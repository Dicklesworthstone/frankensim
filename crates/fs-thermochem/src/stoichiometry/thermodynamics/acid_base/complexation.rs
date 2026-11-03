//! Coupled finite ligand, metal-family and aqueous proton balances.
//!
//! One neutral acid `HL -> H+ + L-` supplies the free monovalent ligand.
//! Each original metal family has consecutive association steps
//! `M_(j-1) + q_j L- -> M_j`; q can differ between steps. Original selected
//! Gibbs values supply the step constants, with no inferred complexes.
//!
//! At a trial free ligand activity l, member weights are beta_j*l^n_j.
//! A reduced ORIGINAL acid/base system treats free L- and charged metal
//! members as spectators and solves charge balance for h. Neutral HL then
//! has amount n_L*h/Ka. The required total n_L+n_HL+sum(n_j*n_Mj) increases
//! with n_L: added ligand lowers metal charge and raises h, and mean bound
//! ligand increases with l. One bounded logarithmic root closes the finite
//! total. The original full acid/base owner supplies the final aqueous report.
//! Near floating logarithmic resolution, a well-resolved component difference
//! total-minus-bound can instead set that report's uncomplexed pool; all actual
//! member, charge, free-ligand bracket and mass-action checks still apply.
//!
//! Mass action, analytical totals and charge balance follow the PHREEQC
//! formulation: <https://pubs.usgs.gov/tm/06/a43/>. This explicit ideal model
//! does not supply metal redox, precipitation, kinetic rates, nonideal activity
//! corrections or stability against omitted complexes. The common n_ref is
//! an aqueous activity-reference AMOUNT, not a solution-volume conversion.

use super::*;
use fs_math::det;

/// At most eight consecutive association steps in one metal component.
pub const MAX_ASSOCIATIONS: usize = 8;
/// At most 32 explicitly declared, mutually disjoint metal components.
pub const MAX_METAL_FAMILIES: usize = 32;
const LN_10: f64 = core::f64::consts::LN_10;

/// Invalid source/topology, unresolved finite balances or bounded-work refusal.
#[derive(Clone, Debug, PartialEq)]
pub enum LigandSpeciationError {
    /// Original aqueous phase/source/charge/speciation owner refused.
    Aqueous(AcidBaseError),
    /// Original exact reaction or selected Gibbs source refused.
    Source(ReactionBalanceError),
    /// Shared original species definitions or reaction identities disagree.
    Network(ReactionNetworkError),
    /// Invalid input, overlapping components or incompatible association topology.
    Invalid(&'static str),
    /// Positive stocks, residuals or requested numerical resolution lose range.
    Arithmetic(&'static str),
    /// The declared free-ligand floor excludes the finite-balance root.
    OutsideBracket,
    /// Outer work or complete aqueous-solve budget was exhausted.
    Budget,
    /// Cancellation occurred before a complete result could be published.
    Cancelled,
}
impl core::fmt::Display for LigandSpeciationError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "ideal metal/ligand speciation: {self:?}")
    }
}
impl std::error::Error for LigandSpeciationError {}
impl From<AcidBaseError> for LigandSpeciationError {
    fn from(e: AcidBaseError) -> Self {
        if e == AcidBaseError::Cancelled {
            Self::Cancelled
        } else {
            Self::Aqueous(e)
        }
    }
}
impl From<ReactionBalanceError> for LigandSpeciationError {
    fn from(e: ReactionBalanceError) -> Self {
        Self::Source(e)
    }
}
impl From<ReactionNetworkError> for LigandSpeciationError {
    fn from(e: ReactionNetworkError) -> Self {
        Self::Network(e)
    }
}
type Result<T> = core::result::Result<T, LigandSpeciationError>;

/// One original exactly conserved aqueous parent + q ligand -> complex step.
#[derive(Clone, Debug)]
pub struct LigandAssociation {
    reaction: ConservedReaction,
    standard: StandardReactionGibbs,
    parent: SpeciesId,
    ligand: SpeciesId,
    complex: SpeciesId,
    count: usize,
    parent_charge: i32,
    complex_charge: i32,
}
impl LigandAssociation {
    /// Bind exactly three distinct aqueous species with parent -1, ligand -q
    /// and complex +1, for q=1..=8. The ligand must have charge -1; original
    /// exact elemental and charge definitions stay in the conserved equation.
    pub fn from_sources(
        reaction: ConservedReaction,
        parent: SpeciesId,
        ligand: SpeciesId,
        complex: SpeciesId,
        sources: &[FormationGibbsSource<'_>],
    ) -> Result<Self> {
        let coefficient = |id: &SpeciesId| {
            reaction
                .terms()
                .iter()
                .find(|t| &t.species == id)
                .map(|t| t.coefficient)
        };
        if parent == ligand
            || parent == complex
            || ligand == complex
            || reaction.terms().len() != 3
            || coefficient(&parent) != Some(-1)
            || coefficient(&complex) != Some(1)
            || !matches!(coefficient(&ligand), Some(n) if (-8..=-1).contains(&n))
            || charge(&reaction, &ligand)? != -1
        {
            return Err(LigandSpeciationError::Invalid(
                "require aqueous parent + 1..=8 monovalent ligands -> complex",
            ));
        }
        let count = (-coefficient(&ligand).expect("admitted ligand")) as usize;
        let parent_charge = charge(&reaction, &parent)?;
        let complex_charge = charge(&reaction, &complex)?;
        if complex_charge != parent_charge - count as i32 {
            return Err(LigandSpeciationError::Invalid(
                "association charge definitions disagree",
            ));
        }
        let standard = reaction.standard_gibbs_change(sources)?;
        for source in sources {
            require_phase(source, &standard, "source-phase-aqueous", false)?;
        }
        Ok(Self {
            reaction,
            standard,
            parent,
            ligand,
            complex,
            count,
            parent_charge,
            complex_charge,
        })
    }
    /// Original conserved association equation.
    pub fn reaction(&self) -> &ConservedReaction {
        &self.reaction
    }
    /// Original stepwise Gibbs sources and natural logarithm of K.
    pub fn standard(&self) -> &StandardReactionGibbs {
        &self.standard
    }
    /// Parent member in this association step.
    pub fn parent(&self) -> &SpeciesId {
        &self.parent
    }
    /// Shared monovalent ligand identity.
    pub fn ligand(&self) -> &SpeciesId {
        &self.ligand
    }
    /// Product member of this association step.
    pub fn complex(&self) -> &SpeciesId {
        &self.complex
    }
    /// Number of ligand molecules added by this step.
    pub const fn ligand_count(&self) -> usize {
        self.count
    }
}

/// One explicit metal component and its complete ordered association chain.
#[derive(Clone, Debug)]
pub struct MetalLigandFamily {
    id: String,
    steps: Vec<LigandAssociation>,
    species: Vec<SpeciesId>,
    charges: Vec<i32>,
    ligand_counts: Vec<usize>,
    log_beta: Vec<f64>,
}
impl MetalLigandFamily {
    /// The root species defines one unit of this component. Ligand counts are
    /// molecules added relative to that root; omitted intermediates are never
    /// inferred. Every step must continue the preceding product exactly.
    pub fn new(id: impl Into<String>, steps: Vec<LigandAssociation>) -> Result<Self> {
        let id = id.into();
        if id.trim().is_empty()
            || id.len() > 256
            || id.chars().any(char::is_control)
            || steps.is_empty()
            || steps.len() > MAX_ASSOCIATIONS
        {
            return Err(LigandSpeciationError::Invalid(
                "bounded identity and 1..=8 association steps required",
            ));
        }
        let first = &steps[0];
        let mut species = vec![first.parent.clone()];
        let mut charges = vec![first.parent_charge];
        let mut ligand_counts = vec![0];
        let mut log_beta = vec![0.0];
        for step in &steps {
            if step.ligand != first.ligand
                || species.last() != Some(&step.parent)
                || species.contains(&step.complex)
                || step.complex == first.ligand
                || step.parent == first.ligand
                || charges.last() != Some(&step.parent_charge)
            {
                return Err(LigandSpeciationError::Invalid(
                    "metal family must be a distinct ordered association chain",
                ));
            }
            species.push(step.complex.clone());
            charges.push(step.complex_charge);
            ligand_counts.push(ligand_counts.last().copied().unwrap() + step.count);
            log_beta.push(add([
                *log_beta.last().unwrap(),
                step.standard.log_equilibrium_constant().value(),
            ])?);
        }
        ReactionNetwork::new(steps.iter().map(|s| s.reaction.clone()).collect())?;
        consistent_sources(steps.iter().map(|s| &s.standard))?;
        Ok(Self {
            id,
            steps,
            species,
            charges,
            ligand_counts,
            log_beta,
        })
    }
    /// Caller-declared analytical component identity.
    pub fn id(&self) -> &str {
        &self.id
    }
    /// Original stepwise association equations and sources.
    pub fn steps(&self) -> &[LigandAssociation] {
        &self.steps
    }
    /// Complete ordered root/complex species axis.
    pub fn species(&self) -> &[SpeciesId] {
        &self.species
    }
    /// Original exact charge of every ordered member.
    pub fn charges(&self) -> &[i32] {
        &self.charges
    }
    /// Added ligand count for every ordered member, starting with zero.
    pub fn ligand_counts(&self) -> &[usize] {
        &self.ligand_counts
    }
    /// Shared original ligand identity.
    pub fn ligand(&self) -> &SpeciesId {
        &self.steps[0].ligand
    }
}

#[derive(Clone, Debug)]
struct LigandLayout {
    aqueous: IdealAcidBaseSystem,
    reduced: IdealAcidBaseSystem,
    ligand_family: String,
    ligand: SpeciesId,
    log_ka: f64,
    families: Vec<MetalLigandFamily>,
    network: ReactionNetwork,
}
/// Immutable coupled ideal aqueous and finite metal/ligand equilibrium model.
#[derive(Clone, Debug)]
pub struct IdealLigandSystem {
    inner: Arc<LigandLayout>,
}
impl IdealLigandSystem {
    /// Join 0..=32 disjoint metal components to one neutral HL/H+/L- acid pair.
    /// Every shared species retains identical original selected Gibbs answers;
    /// all original equations also join the exact species-definition owner.
    pub fn new(
        aqueous: IdealAcidBaseSystem,
        ligand_acid_family_id: impl Into<String>,
        mut families: Vec<MetalLigandFamily>,
    ) -> Result<Self> {
        let ligand_family = ligand_acid_family_id.into();
        let acid = aqueous
            .families()
            .iter()
            .find(|f| f.id() == ligand_family)
            .ok_or(LigandSpeciationError::Invalid(
                "missing ligand acid/base family",
            ))?;
        if acid.steps().len() != 1
            || acid.steps()[0].acid_charge != 0
            || charge(acid.steps()[0].reaction(), acid.steps()[0].base())? != -1
            || families.len() > MAX_METAL_FAMILIES
        {
            return Err(LigandSpeciationError::Invalid(
                "require one neutral HL -> H+ + L- family and at most32 metal families",
            ));
        }
        let ligand = acid.steps()[0].base.clone();
        let log_ka = acid.steps()[0].standard.log_equilibrium_constant().value();
        families.sort_by(|a, b| a.id.cmp(&b.id));
        if families.windows(2).any(|p| p[0].id == p[1].id) {
            return Err(LigandSpeciationError::Invalid(
                "duplicate metal component identity",
            ));
        }
        let mut occupied: BTreeSet<_> = aqueous
            .reaction_network()
            .species()
            .iter()
            .cloned()
            .collect();
        let mut reactions = aqueous.reaction_network().reactions().to_vec();
        let mut standards = vec![aqueous.water().standard()];
        for family in aqueous.families() {
            for step in family.steps() {
                standards.push(step.standard());
            }
        }
        for family in &families {
            if family.ligand() != &ligand || aqueous.families().iter().any(|f| f.id() == family.id)
            {
                return Err(LigandSpeciationError::Invalid(
                    "metal families require the chosen ligand and distinct component IDs",
                ));
            }
            for species in &family.species {
                if !occupied.insert(species.clone()) {
                    return Err(LigandSpeciationError::Invalid(
                        "metal and aqueous components overlap",
                    ));
                }
            }
            for step in &family.steps {
                reactions.push(step.reaction.clone());
                standards.push(&step.standard);
            }
        }
        consistent_sources(standards)?;
        let network = ReactionNetwork::new(reactions)?;
        let reduced = IdealAcidBaseSystem::new(
            aqueous.water().clone(),
            aqueous
                .families()
                .iter()
                .filter(|f| f.id() != ligand_family)
                .cloned()
                .collect(),
        )?;
        Ok(Self {
            inner: Arc::new(LigandLayout {
                aqueous,
                reduced,
                ligand_family,
                ligand,
                log_ka,
                families,
                network,
            }),
        })
    }
    /// Original complete aqueous acid/base model.
    pub fn aqueous_model(&self) -> &IdealAcidBaseSystem {
        &self.inner.aqueous
    }
    /// Explicit neutral-acid/monovalent-ligand analytical family ID.
    pub fn ligand_family_id(&self) -> &str {
        &self.inner.ligand_family
    }
    /// Original metal families in canonical component-ID order.
    pub fn families(&self) -> &[MetalLigandFamily] {
        &self.inner.families
    }
    /// Original admitted complete reaction network and shared definitions.
    pub fn reaction_network(&self) -> &ReactionNetwork {
        &self.inner.network
    }
}

/// Numerical resolution and work bounds; not experimental uncertainty.
#[derive(Clone, Copy, Debug)]
pub struct LigandSpeciationOptions {
    /// Inner original aqueous pH/charge settings.
    pub aqueous: SpeciationOptions,
    /// Positive lower allowed free-ligand amount when an outer root is needed.
    pub minimum_free_ligand: Amount,
    /// Maximum final natural-log free-ligand bracket width.
    pub log_ligand_tolerance: Dimensionless,
    /// Positive absolute ligand-balance allowance [mol].
    pub absolute_ligand_tolerance: Amount,
    /// Additional fraction of the supplied finite ligand total, at most0.01.
    pub relative_ligand_tolerance: Dimensionless,
    /// Maximum recomputed association/acid-base ln(Q/K), including inner pH spread.
    pub log_mass_action_tolerance: Dimensionless,
    /// Bracketing/bisection rounds, at most256; an ambiguous round may probe quarters.
    pub max_iterations: usize,
    /// Complete inner aqueous solves, including final restoration, at most4096.
    pub max_aqueous_solves: usize,
}
impl Default for LigandSpeciationOptions {
    fn default() -> Self {
        Self {
            aqueous: SpeciationOptions {
                ph_tolerance: Dimensionless::new(1e-12),
                absolute_charge_tolerance: Amount::new(1e-15),
                relative_charge_tolerance: Dimensionless::new(1e-12),
                ..SpeciationOptions::default()
            },
            minimum_free_ligand: Amount::new(1e-100),
            log_ligand_tolerance: Dimensionless::new(1e-9),
            absolute_ligand_tolerance: Amount::new(1e-15),
            relative_ligand_tolerance: Dimensionless::new(1e-10),
            log_mass_action_tolerance: Dimensionless::new(1e-8),
            max_iterations: 256,
            max_aqueous_solves: 1024,
        }
    }
}
/// One original metal member and its resolved amount/activity.
#[derive(Clone, Debug)]
pub struct MetalLigandSpecies {
    /// Original phase-specific identity.
    pub species: SpeciesId,
    /// Original exact charge number, including zero for neutral complexes.
    pub charge: i32,
    /// Added ligand molecules relative to the family root.
    pub ligand_count: usize,
    /// Independently represented member amount [mol].
    pub amount: Amount,
    /// Rounded family fraction; may underflow while amount/log_fraction survive.
    pub fraction: Dimensionless,
    /// Log fraction when positive mathematically; None for exact zero-ligand exclusion.
    pub log_fraction: Option<Dimensionless>,
    /// Ideal relative activity on the common n_ref basis.
    pub activity: Dimensionless,
}
/// One finite analytical metal component across its original declared members.
#[derive(Clone, Debug)]
pub struct MetalSpeciation {
    /// Original analytical component identity.
    pub id: String,
    /// Supplied component total [mol].
    pub total: Amount,
    /// Original complete member axis in association order.
    pub members: Vec<MetalLigandSpecies>,
    /// Ligand molecules bound by this component [mol].
    pub bound_ligand: Amount,
    /// Independently reduced sum(member amounts)-total [mol].
    pub component_balance_residual: Amount,
}
/// Complete simultaneous ideal aqueous and finite metal/ligand state.
#[derive(Clone, Debug)]
pub struct LigandSpeciation {
    model: IdealLigandSystem,
    /// Supplied free+neutral-acid+bound ligand total [mol].
    pub ligand_total: Amount,
    /// Independently reduced final neutral HL plus free L member amounts [mol].
    pub uncomplexed_ligand: Amount,
    /// Actual free L- from the final original aqueous report [mol].
    pub free_ligand: Amount,
    /// Sum of ligand molecules bound by all metal families [mol].
    pub bound_ligand: Amount,
    /// Complete original aqueous report; charged metal members appear as spectators.
    pub aqueous: AcidBaseSpeciation,
    /// Complete metal-family inventories, including neutral complexes.
    pub metals: Vec<MetalSpeciation>,
    /// Independently reduced uncomplexed+bound-supplied ligand [mol].
    pub ligand_balance_residual: Amount,
    /// Required-ligand variation across the retained inner pH bracket [mol].
    /// Zero when the final uncomplexed pool follows the well-resolved analytical
    /// component difference instead of an inferred neutral-acid amount.
    pub ligand_balance_numerical_spread: Amount,
    /// Largest recomputed ln(Q/K) residual, excluding equations with empty components.
    pub max_log_mass_action_residual: Dimensionless,
    /// Final numerical free-ligand bracket [mol]; not an interval certificate.
    pub free_ligand_bracket: [Amount; 2],
    /// Outer bracket/bisection rounds.
    pub iterations: usize,
    /// Complete original aqueous solves, including final report restoration.
    pub aqueous_solves: usize,
    /// Sum of original inner pH iteration counts.
    pub aqueous_iterations: usize,
}
impl LigandSpeciation {
    /// Original source/equation owners and component layout.
    pub fn model(&self) -> &IdealLigandSystem {
        &self.model
    }
}

impl IdealLigandSystem {
    /// Solve shared free ligand, pH, every metal component and other acid/base
    /// families together. Every metal/other-acid total must occur exactly once;
    /// fixed spectators cannot duplicate any owned species. No partial result
    /// survives failure or cancellation. The ligand total counts ligands added
    /// relative to each explicitly declared metal-family root.
    #[allow(clippy::too_many_arguments)]
    pub fn speciate(
        &self,
        ligand_total: Amount,
        metal_totals: &[(String, Amount)],
        other_acid_totals: &[(String, Amount)],
        spectators: &[SpectatorIon],
        amount_at_unit_activity: Amount,
        options: LigandSpeciationOptions,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<LigandSpeciation> {
        let basis = require_positive(
            amount_at_unit_activity.value(),
            "positive common aqueous activity amount required",
        )?;
        let total = require_nonnegative(ligand_total.value(), "invalid finite ligand total")?;
        let tolerance = settings(options, total)?;
        if metal_totals.len() != self.families().len() || spectators.len() > 256 {
            return Err(LigandSpeciationError::Invalid(
                "supply every metal total exactly once and at most256 spectators",
            ));
        }
        let mut indexed = BTreeMap::new();
        for (id, n) in metal_totals {
            require_nonnegative(n.value(), "invalid metal component total")?;
            if indexed.insert(id.as_str(), n.value()).is_some()
                || !self.families().iter().any(|f| f.id() == id)
            {
                return Err(LigandSpeciationError::Invalid(
                    "duplicate or unknown metal total",
                ));
            }
        }
        let ordered: Vec<_> = self.families().iter().map(|f| indexed[f.id()]).collect();
        if spectators.iter().any(|s| {
            self.reaction_network()
                .species()
                .binary_search(&s.species)
                .is_ok()
        }) {
            return Err(LigandSpeciationError::Invalid(
                "spectator duplicates an owned aqueous/metal species",
            ));
        }
        let mut work = Work {
            solves: 0,
            iterations: 0,
            maximum: options.max_aqueous_solves,
            cancelled: &mut cancelled,
        };
        work.poll()?;
        let (aqueous, metals, mut bracket, spread, iterations) =
            if total == 0.0 || ordered.iter().all(|n| *n == 0.0) {
                let initial_metals = self.distribute(0.0, &ordered, basis, &mut work)?;
                let ions = metal_ions(spectators, &initial_metals);
                let mut pools = other_acid_totals.to_vec();
                pools.push((self.inner.ligand_family.clone(), ligand_total));
                let aqueous = work.solve(
                    &self.inner.aqueous,
                    &pools,
                    &ions,
                    amount_at_unit_activity,
                    options.aqueous,
                )?;
                let free = ligand_pool(&aqueous, &self.inner.ligand_family, &self.inner.ligand)?.1;
                let metals = self.distribute(free, &ordered, basis, &mut work)?;
                (aqueous, metals, [free; 2], 0.0, 0)
            } else {
                let floor = options.minimum_free_ligand.value();
                if floor >= total {
                    return Err(LigandSpeciationError::OutsideBracket);
                }
                let mut hi = det::ln(total);
                let log_floor = det::ln(floor);
                let mut probe = hi;
                let mut step = core::f64::consts::LN_2;
                let mut iterations = 0;
                let mut lo;
                // The exact inventory bound n_L<=total supplies the upper endpoint.
                // Expand down only as far as needed, avoiding irrelevant underflowing
                // high-order complex tails at an unnecessarily remote search floor.
                loop {
                    advance_iteration(&mut iterations, options.max_iterations)?;
                    probe = (probe - step).max(log_floor);
                    let n = if probe == log_floor {
                        floor
                    } else {
                        exp_positive(probe, "free ligand bracket loses range")?
                    };
                    let trial = self.trial(
                        n,
                        total,
                        &ordered,
                        other_acid_totals,
                        spectators,
                        amount_at_unit_activity,
                        options,
                        &mut work,
                    )?;
                    if trial.residual + trial.spread < 0.0 {
                        lo = det::ln(n);
                        break;
                    }
                    if trial.residual - trial.spread > 0.0 {
                        hi = det::ln(n);
                    }
                    if probe == log_floor {
                        return Err(if trial.residual - trial.spread > 0.0 {
                            LigandSpeciationError::OutsideBracket
                        } else {
                            LigandSpeciationError::Arithmetic(
                                "inner pH bracket cannot resolve free-ligand floor",
                            )
                        });
                    }
                    step *= 2.0;
                }
                let (accepted, balanced_pool);
                loop {
                    advance_iteration(&mut iterations, options.max_iterations)?;
                    let mid = lo.midpoint(hi);
                    if mid <= lo || mid >= hi {
                        return Err(LigandSpeciationError::Arithmetic(
                            "free ligand root reached logarithmic resolution",
                        ));
                    }
                    let n = exp_positive(mid, "positive free ligand loses range")?;
                    let trial = self.trial(
                        n,
                        total,
                        &ordered,
                        other_acid_totals,
                        spectators,
                        amount_at_unit_activity,
                        options,
                        &mut work,
                    )?;
                    if hi - lo <= options.log_ligand_tolerance.value()
                        && trial.residual.abs() + trial.spread <= tolerance
                    {
                        (accepted, balanced_pool) = (trial, None);
                        break;
                    }
                    // Near charge cancellation, one free-ligand ULP can imply
                    // a large change in inferred HL even when the free-ligand
                    // root is already sharply bracketed. Eliminate the total
                    // analytically BEFORE restoring the original full aqueous
                    // solve, provided the difference itself is well resolved.
                    // A trace uncomplexed pool must keep its independent sum;
                    // final member, bracket and mass-action checks are mandatory.
                    if hi - lo <= options.log_ligand_tolerance.value()
                        && hi - lo <= 64.0 * f64::EPSILON * lo.abs().max(hi.abs()).max(1.0)
                    {
                        let bound = add_iter(trial.metals.iter().map(|m| m.bound_ligand.value()))?;
                        if let Some(pool) = resolved_component_difference(
                            total,
                            bound,
                            options.log_mass_action_tolerance.value(),
                        ) {
                            (accepted, balanced_pool) = (trial, Some(pool));
                            break;
                        }
                    }
                    if trial.residual + trial.spread < 0.0 {
                        lo = mid;
                        continue;
                    }
                    if trial.residual - trial.spread > 0.0 {
                        hi = mid;
                        continue;
                    }
                    if hi - lo <= options.log_ligand_tolerance.value() {
                        return Err(LigandSpeciationError::Arithmetic(
                            "inner pH resolution cannot meet ligand tolerance",
                        ));
                    }
                    // Do not discard a possible sign when the inner pH bracket
                    // straddles zero. Quarter probes can shrink both sides safely.
                    let before = [lo, hi];
                    for quarter in [lo.midpoint(mid), mid.midpoint(hi)] {
                        if quarter <= lo || quarter >= hi {
                            continue;
                        }
                        let t = self.trial(
                            exp_positive(quarter, "free ligand quarter loses range")?,
                            total,
                            &ordered,
                            other_acid_totals,
                            spectators,
                            amount_at_unit_activity,
                            options,
                            &mut work,
                        )?;
                        if t.residual + t.spread < 0.0 {
                            lo = quarter;
                        } else if t.residual - t.spread > 0.0 {
                            hi = quarter;
                        }
                    }
                    if before == [lo, hi] {
                        return Err(LigandSpeciationError::Arithmetic(
                            "inner pH bracket leaves ligand sign unresolved",
                        ));
                    }
                }
                let uncomplexed = match balanced_pool {
                    Some(pool) => pool,
                    None => add([accepted.free, accepted.neutral])?,
                };
                let mut pools = other_acid_totals.to_vec();
                pools.push((self.inner.ligand_family.clone(), Amount::new(uncomplexed)));
                let ions = metal_ions(spectators, &accepted.metals);
                let aqueous = work.solve(
                    &self.inner.aqueous,
                    &pools,
                    &ions,
                    amount_at_unit_activity,
                    options.aqueous,
                )?;
                (
                    aqueous,
                    accepted.metals,
                    [
                        exp_positive(lo, "lower ligand bracket loses range")?,
                        exp_positive(hi, "upper ligand bracket loses range")?,
                    ],
                    if balanced_pool.is_some() {
                        0.0
                    } else {
                        accepted.spread
                    },
                    iterations,
                )
            };
        let (uncomplexed, free) =
            ligand_pool(&aqueous, &self.inner.ligand_family, &self.inner.ligand)?;
        let bound = add_iter(metals.iter().map(|m| m.bound_ligand.value()))?;
        let residual = add([uncomplexed, bound, -total])?;
        if residual.abs() + spread > tolerance {
            return Err(LigandSpeciationError::Arithmetic(
                "final ligand inventory does not meet requested tolerance",
            ));
        }
        bracket[0] = bracket[0].min(free);
        bracket[1] = bracket[1].max(free);
        if free > 0.0
            && det::ln(bracket[1]) - det::ln(bracket[0]) > options.log_ligand_tolerance.value()
        {
            return Err(LigandSpeciationError::Arithmetic(
                "restored aqueous report exceeds free-ligand bracket tolerance",
            ));
        }
        let mut maximum = aqueous.max_log_mass_action_residual.value();
        let mut largest_count = 0usize;
        if free > 0.0 {
            let ln_l = det::ln(free) - det::ln(basis);
            for (family, metal) in self.families().iter().zip(&metals) {
                work.poll()?;
                if metal.total.value() == 0.0 {
                    continue;
                }
                for (j, step) in family.steps.iter().enumerate() {
                    let r = add([
                        det::ln(metal.members[j + 1].activity.value()),
                        -det::ln(metal.members[j].activity.value()),
                        -(step.count as f64) * ln_l,
                        -step.standard.log_equilibrium_constant().value(),
                    ])?;
                    maximum = maximum.max(r.abs());
                    largest_count = largest_count.max(step.count);
                }
            }
        }
        let ph_spread = aqueous
            .ph_bracket
            .iter()
            .map(|p| (p.value() - aqueous.ph.value()).abs())
            .fold(0.0, f64::max);
        if maximum + (largest_count as f64) * LN_10 * ph_spread
            > options.log_mass_action_tolerance.value()
        {
            return Err(LigandSpeciationError::Arithmetic(
                "final aqueous/metal mass action is insufficiently resolved",
            ));
        }
        work.poll()?;
        Ok(LigandSpeciation {
            model: self.clone(),
            ligand_total,
            uncomplexed_ligand: Amount::new(uncomplexed),
            free_ligand: Amount::new(free),
            bound_ligand: Amount::new(bound),
            aqueous,
            metals,
            ligand_balance_residual: Amount::new(residual),
            ligand_balance_numerical_spread: Amount::new(spread),
            max_log_mass_action_residual: Dimensionless::new(maximum),
            free_ligand_bracket: bracket.map(Amount::new),
            iterations,
            aqueous_solves: work.solves,
            aqueous_iterations: work.iterations,
        })
    }

    #[allow(clippy::too_many_arguments)]
    fn trial(
        &self,
        free: f64,
        total: f64,
        metal_totals: &[f64],
        other_totals: &[(String, Amount)],
        spectators: &[SpectatorIon],
        basis: Amount,
        options: LigandSpeciationOptions,
        work: &mut Work<'_, impl FnMut() -> bool>,
    ) -> Result<Trial> {
        let metals = self.distribute(free, metal_totals, basis.value(), work)?;
        let mut ions = metal_ions(spectators, &metals);
        ions.push(SpectatorIon {
            species: self.inner.ligand.clone(),
            charge: -1,
            amount: Amount::new(free),
        });
        let aq = work.solve(
            &self.inner.reduced,
            other_totals,
            &ions,
            basis,
            options.aqueous,
        )?;
        let neutral = exp_positive(
            add([
                det::ln(free),
                det::ln(aq.proton_activity.value()),
                -self.inner.log_ka,
            ])?,
            "neutral ligand amount loses range",
        )?;
        let bound = add_iter(metals.iter().map(|m| m.bound_ligand.value()))?;
        let residual = add([free, neutral, bound, -total])?;
        let spread = aq
            .ph_bracket
            .iter()
            .map(|p| neutral * det::expm1(LN_10 * (aq.ph.value() - p.value())).abs())
            .fold(0.0, f64::max);
        if !spread.is_finite() {
            return Err(LigandSpeciationError::Arithmetic(
                "inner pH ligand spread loses range",
            ));
        }
        Ok(Trial {
            free,
            neutral,
            metals,
            residual,
            spread,
        })
    }
    fn distribute(
        &self,
        free: f64,
        totals: &[f64],
        basis: f64,
        work: &mut Work<'_, impl FnMut() -> bool>,
    ) -> Result<Vec<MetalSpeciation>> {
        let mut result = Vec::with_capacity(self.families().len());
        for (family, &total) in self.families().iter().zip(totals) {
            work.poll()?;
            let logs = if free == 0.0 {
                vec![0.0; family.species.len()]
            } else {
                let l = det::ln(free) - det::ln(basis);
                family
                    .log_beta
                    .iter()
                    .zip(&family.ligand_counts)
                    .map(|(&b, &n)| add([b, (n as f64) * l]))
                    .collect::<Result<Vec<_>>>()?
            };
            let max = logs.iter().copied().fold(f64::NEG_INFINITY, f64::max);
            let denominator = if free == 0.0 {
                0.0
            } else {
                det::ln(add_iter(logs.iter().map(|w| det::exp(*w - max)))?)
            };
            let mut members = Vec::with_capacity(logs.len());
            for (j, &log) in logs.iter().enumerate() {
                let excluded = free == 0.0 && j > 0;
                let lf = if free == 0.0 {
                    0.0
                } else {
                    add([log, -max, -denominator])?
                };
                let fraction = if excluded { 0.0 } else { det::exp(lf) };
                let amount = if excluded || total == 0.0 {
                    0.0
                } else {
                    let ordinary = total * fraction;
                    if ordinary > 0.0 {
                        ordinary
                    } else {
                        exp_positive(
                            add([det::ln(total), lf])?,
                            "positive metal member loses range",
                        )?
                    }
                };
                let activity = if amount == 0.0 {
                    0.0
                } else {
                    computed_positive(amount / basis, "positive metal activity loses range")?
                };
                members.push(MetalLigandSpecies {
                    species: family.species[j].clone(),
                    charge: family.charges[j],
                    ligand_count: family.ligand_counts[j],
                    amount: Amount::new(amount),
                    fraction: Dimensionless::new(fraction),
                    log_fraction: if excluded {
                        None
                    } else {
                        Some(Dimensionless::new(lf))
                    },
                    activity: Dimensionless::new(activity),
                });
            }
            let balance = add_iter(
                members
                    .iter()
                    .map(|m| m.amount.value())
                    .chain(std::iter::once(-total)),
            )?;
            if balance.abs() > 128.0 * f64::EPSILON * total {
                return Err(LigandSpeciationError::Arithmetic(
                    "metal component balance loses closure",
                ));
            }
            let bound = add_iter(
                members
                    .iter()
                    .map(|m| m.amount.value() * (m.ligand_count as f64)),
            )?;
            result.push(MetalSpeciation {
                id: family.id.clone(),
                total: Amount::new(total),
                members,
                bound_ligand: Amount::new(bound),
                component_balance_residual: Amount::new(balance),
            });
        }
        Ok(result)
    }
}

struct Trial {
    free: f64,
    neutral: f64,
    metals: Vec<MetalSpeciation>,
    residual: f64,
    spread: f64,
}

// This guard only selects a numerical formulation; acceptance still requires
// the final independently recomputed component, charge and equilibrium checks.
fn resolved_component_difference(total: f64, bound: f64, log_accuracy: f64) -> Option<f64> {
    let pool = total - bound;
    (pool.is_finite() && pool > 0.0 && 128.0 * f64::EPSILON * (total / pool) <= log_accuracy / 16.0)
        .then_some(pool)
}
struct Work<'a, F: FnMut() -> bool> {
    solves: usize,
    iterations: usize,
    maximum: usize,
    cancelled: &'a mut F,
}
impl<F: FnMut() -> bool> Work<'_, F> {
    fn poll(&mut self) -> Result<()> {
        if (self.cancelled)() {
            Err(LigandSpeciationError::Cancelled)
        } else {
            Ok(())
        }
    }
    fn solve(
        &mut self,
        model: &IdealAcidBaseSystem,
        totals: &[(String, Amount)],
        ions: &[SpectatorIon],
        basis: Amount,
        options: SpeciationOptions,
    ) -> Result<AcidBaseSpeciation> {
        self.poll()?;
        if self.solves >= self.maximum {
            return Err(LigandSpeciationError::Budget);
        }
        self.solves += 1;
        let answer = model.speciate(totals, ions, basis, options, &mut self.cancelled)?;
        self.iterations += answer.iterations;
        self.poll()?;
        Ok(answer)
    }
}
fn metal_ions(spectators: &[SpectatorIon], metals: &[MetalSpeciation]) -> Vec<SpectatorIon> {
    spectators
        .iter()
        .cloned()
        .chain(
            metals
                .iter()
                .flat_map(|m| m.members.iter())
                .filter(|m| m.charge != 0)
                .map(|m| SpectatorIon {
                    species: m.species.clone(),
                    charge: m.charge,
                    amount: m.amount,
                }),
        )
        .collect()
}
fn ligand_pool(aqueous: &AcidBaseSpeciation, id: &str, ligand: &SpeciesId) -> Result<(f64, f64)> {
    let family =
        aqueous
            .families
            .iter()
            .find(|f| f.id == id)
            .ok_or(LigandSpeciationError::Invalid(
                "missing restored ligand family",
            ))?;
    let free = family.forms.iter().find(|f| &f.species == ligand).ok_or(
        LigandSpeciationError::Invalid("missing restored free ligand"),
    )?;
    Ok((
        add_iter(family.forms.iter().map(|f| f.amount.value()))?,
        free.amount.value(),
    ))
}
fn settings(o: LigandSpeciationOptions, total: f64) -> Result<f64> {
    require_positive(
        o.minimum_free_ligand.value(),
        "positive free ligand floor required",
    )?;
    for x in [o.log_ligand_tolerance, o.log_mass_action_tolerance] {
        require_positive(x.value(), "positive logarithmic tolerance required")?;
    }
    require_positive(
        o.absolute_ligand_tolerance.value(),
        "positive absolute ligand tolerance required",
    )?;
    if !o.relative_ligand_tolerance.value().is_finite()
        || !(0.0..=0.01).contains(&o.relative_ligand_tolerance.value())
        || !(1..=256).contains(&o.max_iterations)
        || !(1..=4096).contains(&o.max_aqueous_solves)
    {
        return Err(LigandSpeciationError::Invalid(
            "invalid relative ligand tolerance or work bounds",
        ));
    }
    add([
        o.absolute_ligand_tolerance.value(),
        o.relative_ligand_tolerance.value() * total,
    ])
}
fn advance_iteration(iterations: &mut usize, maximum: usize) -> Result<()> {
    if *iterations >= maximum {
        Err(LigandSpeciationError::Budget)
    } else {
        *iterations += 1;
        Ok(())
    }
}
fn require_nonnegative(x: f64, reason: &'static str) -> Result<f64> {
    if x.is_finite() && x >= 0.0 {
        Ok(x)
    } else {
        Err(LigandSpeciationError::Invalid(reason))
    }
}
fn require_positive(x: f64, reason: &'static str) -> Result<f64> {
    if x.is_finite() && x > 0.0 {
        Ok(x)
    } else {
        Err(LigandSpeciationError::Invalid(reason))
    }
}
fn computed_positive(x: f64, reason: &'static str) -> Result<f64> {
    if x.is_finite() && x > 0.0 {
        Ok(x)
    } else {
        Err(LigandSpeciationError::Arithmetic(reason))
    }
}
fn exp_positive(x: f64, reason: &'static str) -> Result<f64> {
    computed_positive(det::exp(x), reason)
}
fn add<const N: usize>(values: [f64; N]) -> Result<f64> {
    add_iter(values.into_iter())
}
fn add_iter(values: impl IntoIterator<Item = f64>) -> Result<f64> {
    let mut s = Sum::default();
    for v in values {
        s.add(v)?;
    }
    Ok(s.finish()?)
}

#[cfg(test)]
#[path = "complexation/tests.rs"]
mod tests;
