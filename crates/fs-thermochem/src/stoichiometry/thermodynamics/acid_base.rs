//! Coupled ideal aqueous acid/base families and water autoprotolysis.
//!
//! Each conserved one-proton step is S_j -> H+ + S_(j+1). Its source-bound
//! thermodynamic K supplies the relative weights w_j = product(K_1..K_j)/aH^j.
//! All families share the SAME proton activity, and water supplies aOH=Kw/aH.
//! One monotone electroneutrality equation closes the analytical family totals
//! and caller-declared spectator ions. Sequential single-reaction equilibria
//! cannot provide this coupled closure.
//!
//! The caller supplies a common amount at unit activity: a_i=n_i/n_ref. For
//! unit-molality references, n_ref is solvent mass times 1 mol/kg, not volume
//! times 1 mol/m³. The solvent basis stays fixed and water activity is one;
//! water is not a finite reactant inventory here. Family totals conserve their
//! non-proton components, but no closed H/O solvent balance is asserted.
//! pH=-log10(aH), following IUPAC <https://goldbook.iupac.org/terms/view/P04524>.
//! Mass-action and charge-balance conventions follow USGS PHREEQC:
//! <https://wwwbrr.cr.usgs.gov/projects/GWC_coupled/phreeqc.v1/html/phqc_11.htm>.
//!
//! This is an explicit ideal activity approximation, not measured pH, a
//! concentrated-electrolyte model, complexation, redox, precipitation, kinetics,
//! or plant calibration. In particular HCN+CN- is the uncomplexed cyanide
//! family; it is not total cyanide including gold/copper complexes or WAD-CN.

use super::{FormationGibbsSource, StandardReactionGibbs, Sum};
use crate::{
    SpeciesId,
    stoichiometry::{
        ConservedReaction, ReactionBalanceError,
        network::{ReactionNetwork, ReactionNetworkError},
    },
};
use fs_matdb::MaterialAnswer;
use fs_qty::{Amount, Dimensionless, Dims, QuantitySpec};
use std::{
    collections::{BTreeMap, BTreeSet},
    sync::Arc,
};

/// Finite shared ligand competition between source-bound aqueous metal families.
pub mod complexation;
/// Finite gas inventories coupled to the same aqueous charge and proton balance.
pub mod gas_liquid;
/// Coupled charge balance, acid/base families and one finite hydrolyzing solid.
pub mod precipitation;
#[path = "acid_base/solve.rs"]
mod solve;

/// At most eight successive one-proton dissociations in one analytical family.
pub const MAX_DISSOCIATIONS: usize = 8;
/// Resource limit; the existing conserved network also bounds total species.
pub const MAX_FAMILIES: usize = 32;

/// No unsupported source or failed solve returns a partially speciated state.
#[derive(Clone, Debug, PartialEq)]
pub enum AcidBaseError {
    /// Original exact chemistry or material-source admission.
    Balance(ReactionBalanceError),
    /// Original shared-species identity admission.
    Network(ReactionNetworkError),
    /// Original source-bound neutral gas/aqueous partition refused.
    GasLiquid(super::gas_liquid::GasLiquidError),
    /// Invalid topology, activity basis, inventory or numerical limit.
    Invalid(&'static str),
    /// A required positive species or balance cannot be represented.
    Arithmetic(&'static str),
    /// The caller's pH interval does not bracket electroneutrality.
    OutsideBracket,
    /// The requested charge/pH accuracy exceeded the iteration limit.
    Budget,
    /// Cancellation was requested before publishing the result.
    Cancelled,
}
impl core::fmt::Display for AcidBaseError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "ideal aqueous acid/base speciation: {self:?}")
    }
}
impl std::error::Error for AcidBaseError {}
impl From<ReactionBalanceError> for AcidBaseError {
    fn from(e: ReactionBalanceError) -> Self {
        Self::Balance(e)
    }
}
impl From<ReactionNetworkError> for AcidBaseError {
    fn from(e: ReactionNetworkError) -> Self {
        Self::Network(e)
    }
}
impl From<super::gas_liquid::GasLiquidError> for AcidBaseError {
    fn from(e: super::gas_liquid::GasLiquidError) -> Self {
        Self::GasLiquid(e)
    }
}
type Result<T> = core::result::Result<T, AcidBaseError>;

/// One original conserved one-proton dissociation and its selected Gibbs sources.
#[derive(Clone, Debug)]
pub struct AcidDissociation {
    reaction: ConservedReaction,
    standard: StandardReactionGibbs,
    acid: SpeciesId,
    base: SpeciesId,
    proton: SpeciesId,
    acid_charge: i32,
}
impl AcidDissociation {
    /// Require exactly acid -> base + H+, all aqueous, with the original exact
    /// H atom and +1 charge definition. Reversed or rescaled equations refuse.
    pub fn from_sources(
        reaction: ConservedReaction,
        acid: SpeciesId,
        base: SpeciesId,
        proton: SpeciesId,
        sources: &[FormationGibbsSource<'_>],
    ) -> Result<Self> {
        require_terms(&reaction, &[(&acid, -1), (&base, 1), (&proton, 1)])?;
        require_formula(&reaction, &proton, &[("H", 1)], 1)?;
        let acid_charge = charge(&reaction, &acid)?;
        let base_charge = charge(&reaction, &base)?;
        if base_charge != acid_charge - 1 {
            return Err(AcidBaseError::Invalid(
                "acid and base must differ by one proton charge",
            ));
        }
        let standard = reaction.standard_gibbs_change(sources)?;
        for source in sources {
            require_phase(source, &standard, "source-phase-aqueous", false)?;
        }
        Ok(Self {
            reaction,
            standard,
            acid,
            base,
            proton,
            acid_charge,
        })
    }
    /// Original equation and its exact atom/charge certificate.
    pub fn reaction(&self) -> &ConservedReaction {
        &self.reaction
    }
    /// Original selected sources, reference state and natural logarithm of Ka.
    pub fn standard(&self) -> &StandardReactionGibbs {
        &self.standard
    }
    /// Protonated member of this step.
    pub fn acid(&self) -> &SpeciesId {
        &self.acid
    }
    /// Member after losing one proton.
    pub fn base(&self) -> &SpeciesId {
        &self.base
    }
    /// Shared aqueous hydrogen ion.
    pub fn proton(&self) -> &SpeciesId {
        &self.proton
    }
}

/// One ordered, non-overlapping analytical family, including every intermediate.
#[derive(Clone, Debug)]
pub struct AcidFamily {
    id: String,
    steps: Vec<AcidDissociation>,
    species: Vec<SpeciesId>,
}
impl AcidFamily {
    /// The caller names the analytical pool; names do not infer constituents.
    /// One to eight original steps must form an unbranched dissociation chain.
    pub fn new(id: impl Into<String>, steps: Vec<AcidDissociation>) -> Result<Self> {
        let id = id.into();
        if id.trim().is_empty()
            || id.len() > 256
            || steps.is_empty()
            || steps.len() > MAX_DISSOCIATIONS
        {
            return Err(AcidBaseError::Invalid(
                "bounded family identity and 1..=8 proton steps required",
            ));
        }
        let first = &steps[0];
        let mut species = vec![first.acid.clone()];
        for step in &steps {
            if step.proton != first.proton
                || species.last() != Some(&step.acid)
                || species.contains(&step.base)
                || step.base == first.proton
                || step.acid == first.proton
            {
                return Err(AcidBaseError::Invalid(
                    "family must be a distinct ordered one-proton chain",
                ));
            }
            species.push(step.base.clone());
        }
        ReactionNetwork::new(steps.iter().map(|s| s.reaction.clone()).collect())?;
        consistent_sources(steps.iter().map(|s| &s.standard))?;
        Ok(Self { id, steps, species })
    }
    /// Explicit analytical-pool identity.
    pub fn id(&self) -> &str {
        &self.id
    }
    /// Original consecutive dissociation steps, most protonated first.
    pub fn steps(&self) -> &[AcidDissociation] {
        &self.steps
    }
    /// Complete ordered member axis for this analytical total.
    pub fn species(&self) -> &[SpeciesId] {
        &self.species
    }
}

/// Source-bound H2O(l) -> H+(aq) + OH-(aq), with unit water activity.
#[derive(Clone, Debug)]
pub struct WaterAutoprotolysis {
    reaction: ConservedReaction,
    standard: StandardReactionGibbs,
    water: SpeciesId,
    proton: SpeciesId,
    hydroxide: SpeciesId,
}
impl WaterAutoprotolysis {
    /// Require exact water/proton/hydroxide definitions and the original pure
    /// liquid/aqueous source phases. This does not admit a finite solvent stock.
    pub fn from_sources(
        reaction: ConservedReaction,
        water: SpeciesId,
        proton: SpeciesId,
        hydroxide: SpeciesId,
        sources: &[FormationGibbsSource<'_>],
    ) -> Result<Self> {
        require_terms(&reaction, &[(&water, -1), (&proton, 1), (&hydroxide, 1)])?;
        require_formula(&reaction, &water, &[("H", 2), ("O", 1)], 0)?;
        require_formula(&reaction, &proton, &[("H", 1)], 1)?;
        require_formula(&reaction, &hydroxide, &[("H", 1), ("O", 1)], -1)?;
        let standard = reaction.standard_gibbs_change(sources)?;
        for source in sources {
            let solvent = *source.species == water;
            require_phase(
                source,
                &standard,
                if solvent {
                    "source-phase-liquid"
                } else {
                    "source-phase-aqueous"
                },
                solvent,
            )?;
        }
        Ok(Self {
            reaction,
            standard,
            water,
            proton,
            hydroxide,
        })
    }
    /// Original solvent dissociation and its exact conservation certificate.
    pub fn reaction(&self) -> &ConservedReaction {
        &self.reaction
    }
    /// Original reference state, property receipts and natural logarithm of Kw.
    pub fn standard(&self) -> &StandardReactionGibbs {
        &self.standard
    }
    /// Fixed-activity solvent identity.
    pub fn water(&self) -> &SpeciesId {
        &self.water
    }
    /// Shared aqueous proton identity.
    pub fn proton(&self) -> &SpeciesId {
        &self.proton
    }
    /// Shared aqueous hydroxide identity.
    pub fn hydroxide(&self) -> &SpeciesId {
        &self.hydroxide
    }
}

#[derive(Clone, Debug)]
struct Layout {
    water: WaterAutoprotolysis,
    families: Vec<AcidFamily>,
    network: ReactionNetwork,
}

/// Immutable source-backed ideal acid/base model; no mutable trial inventories.
#[derive(Clone, Debug)]
pub struct IdealAcidBaseSystem {
    inner: Arc<Layout>,
}
impl IdealAcidBaseSystem {
    /// Join zero to 32 independent analytical families sharing the original
    /// water/proton reference. Shared source receipts must be identical; every
    /// reaction also passes original network species-definition admission.
    pub fn new(water: WaterAutoprotolysis, mut families: Vec<AcidFamily>) -> Result<Self> {
        if families.len() > MAX_FAMILIES {
            return Err(AcidBaseError::Invalid("too many acid/base families"));
        }
        families.sort_by(|a, b| a.id.cmp(&b.id));
        if families.windows(2).any(|p| p[0].id == p[1].id) {
            return Err(AcidBaseError::Invalid(
                "duplicate analytical family identity",
            ));
        }
        let mut species = BTreeSet::from([
            water.water.clone(),
            water.proton.clone(),
            water.hydroxide.clone(),
        ]);
        let mut reactions = vec![water.reaction.clone()];
        let mut standards = vec![&water.standard];
        for family in &families {
            if family.steps[0].proton != water.proton {
                return Err(AcidBaseError::Invalid(
                    "families must share the water proton identity",
                ));
            }
            for id in &family.species {
                if !species.insert(id.clone()) {
                    return Err(AcidBaseError::Invalid(
                        "analytical families overlap or contain water/proton/hydroxide",
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
        Ok(Self {
            inner: Arc::new(Layout {
                water,
                families,
                network,
            }),
        })
    }
    /// Original water source and fixed-solvent convention.
    pub fn water(&self) -> &WaterAutoprotolysis {
        &self.inner.water
    }
    /// Original analytical families in canonical identity order.
    pub fn families(&self) -> &[AcidFamily] {
        &self.inner.families
    }
    /// Original equations; retained for identity and exact conservation checks.
    /// Its simultaneous finite-stock advance is NOT used for speciation.
    pub fn reaction_network(&self) -> &ReactionNetwork {
        &self.inner.network
    }
}

/// An explicit nonreacting aqueous ion, kept outside every acid/base family.
/// Charge is a caller declaration; no charge or counterion is parsed from a name.
#[derive(Clone, Debug, PartialEq)]
pub struct SpectatorIon {
    /// Distinct identity retained in the result.
    pub species: SpeciesId,
    /// Signed ionic charge number, bounded to -32..=32 and nonzero.
    pub charge: i32,
    /// Fixed nonnegative amount [mol].
    pub amount: Amount,
}

/// Numerical limits for the monotone charge balance, separate from uncertainty.
#[derive(Clone, Copy, Debug)]
pub struct SpeciationOptions {
    /// Inclusive increasing pH bracket. Unsupported brackets refuse, not widen.
    pub ph_bracket: [Dimensionless; 2],
    /// Positive maximum numerical pH bracket width.
    pub ph_tolerance: Dimensionless,
    /// Positive absolute charge-equivalent residual allowance [mol].
    pub absolute_charge_tolerance: Amount,
    /// Relative residual allowance against total unsigned ionic charge.
    pub relative_charge_tolerance: Dimensionless,
    /// At most 1024 bisections, with bounded per-family cancellation checks.
    pub max_iterations: usize,
}
impl Default for SpeciationOptions {
    fn default() -> Self {
        Self {
            ph_bracket: [Dimensionless::new(-2.0), Dimensionless::new(16.0)],
            ph_tolerance: Dimensionless::new(1e-10),
            absolute_charge_tolerance: Amount::new(1e-12),
            relative_charge_tolerance: Dimensionless::new(1e-10),
            max_iterations: 128,
        }
    }
}

/// One member's amount, fraction and ideal activity on the declared common basis.
#[derive(Clone, Debug, PartialEq)]
pub struct AcidBaseSpecies {
    /// Original phase-specific member identity.
    pub species: SpeciesId,
    /// Nonnegative member amount [mol].
    pub amount: Amount,
    /// Equilibrium family fraction; also meaningful for an empty family.
    pub fraction: Dimensionless,
    /// Ideal relative activity n/n_ref; zero for an empty family.
    pub activity: Dimensionless,
}
/// Explicit analytical total and all its source-defined equilibrium members.
#[derive(Clone, Debug, PartialEq)]
pub struct FamilySpeciation {
    /// Original analytical-pool identity.
    pub id: String,
    /// Original conserved analytical amount, not total metal/cyanide inferred.
    pub total: Amount,
    /// Members in original protonation order.
    pub forms: Vec<AcidBaseSpecies>,
}
/// One simultaneous ideal-model solution, retaining every original source owner.
#[derive(Clone, Debug)]
pub struct AcidBaseSpeciation {
    model: IdealAcidBaseSystem,
    /// Same explicit amount scale used by every aqueous participant [mol].
    pub amount_at_unit_activity: Amount,
    /// -log10(aH), under the admitted ideal-activity approximation.
    pub ph: Dimensionless,
    /// Relative proton activity, not mol/m³.
    pub proton_activity: Dimensionless,
    /// Relative hydroxide activity, not mol/m³.
    pub hydroxide_activity: Dimensionless,
    /// Free proton amount [mol] on the declared solvent basis.
    pub proton_amount: Amount,
    /// Free hydroxide amount [mol] on the declared solvent basis.
    pub hydroxide_amount: Amount,
    /// Canonically ordered analytical-family results.
    pub families: Vec<FamilySpeciation>,
    /// Unchanged caller-declared nonreacting ions, in canonical order.
    pub spectators: Vec<SpectatorIon>,
    /// Sum(z_i n_i) [mol charge equivalents], not a charge certificate.
    pub charge_balance_residual: Amount,
    /// Largest recomputed absolute ln(Q/K); empty-family equations are omitted.
    pub max_log_mass_action_residual: Dimensionless,
    /// Final numerical pH bracket; not an interval-arithmetic certificate.
    pub ph_bracket: [Dimensionless; 2],
    /// Number of bisections used after the two endpoint evaluations.
    pub iterations: usize,
}
impl AcidBaseSpeciation {
    /// Original water, family equations and source receipts.
    pub fn model(&self) -> &IdealAcidBaseSystem {
        &self.model
    }
}

fn require_terms(reaction: &ConservedReaction, terms: &[(&SpeciesId, i128)]) -> Result<()> {
    let ids: BTreeSet<_> = terms.iter().map(|(id, _)| *id).collect();
    if ids.len() != terms.len()
        || reaction.terms().len() != terms.len()
        || terms.iter().any(|(id, nu)| {
            !reaction
                .terms()
                .iter()
                .any(|t| t.species == **id && t.coefficient == *nu)
        })
    {
        return Err(AcidBaseError::Invalid(
            "requires original unscaled one-proton dissociation equation",
        ));
    }
    Ok(())
}
fn require_formula(
    reaction: &ConservedReaction,
    id: &SpeciesId,
    atoms: &[(&str, i128)],
    charge: i128,
) -> Result<()> {
    let index =
        reaction
            .terms
            .iter()
            .position(|t| &t.species == id)
            .ok_or(AcidBaseError::Invalid(
                "missing declared acid/base participant",
            ))?;
    let definition = &reaction.definitions[index];
    if definition.charge != charge
        || definition.elements.len() != atoms.len()
        || atoms.iter().any(|&(symbol, count)| {
            !definition
                .elements
                .iter()
                .any(|(e, n)| e.as_str() == symbol && *n == count)
        })
    {
        return Err(AcidBaseError::Invalid(
            "water/proton/hydroxide exact atom or charge identity differs",
        ));
    }
    Ok(())
}
fn charge(reaction: &ConservedReaction, id: &SpeciesId) -> Result<i32> {
    let index =
        reaction
            .terms
            .iter()
            .position(|t| &t.species == id)
            .ok_or(AcidBaseError::Invalid(
                "missing acid/base charge definition",
            ))?;
    let value = reaction.definitions[index].charge;
    if !(-32..=32).contains(&value) {
        return Err(AcidBaseError::Invalid(
            "acid/base charge exceeds bounded model",
        ));
    }
    Ok(value as i32)
}
fn require_phase(
    source: &FormationGibbsSource<'_>,
    standard: &StandardReactionGibbs,
    phase: &str,
    pure: bool,
) -> Result<()> {
    let evidence = standard
        .terms()
        .iter()
        .find(|t| t.term.species == *source.species)
        .ok_or(AcidBaseError::Invalid("missing selected species source"))?;
    let claim = source
        .pack
        .claims()
        .claim(evidence.answer.receipt.selected)
        .ok_or(AcidBaseError::Invalid(
            "missing selected acid/base source claim",
        ))?;
    let owns = |axis: &str| {
        claim.validity.bound(axis) == Some((1.0, 1.0))
            && claim.validity.axis_quantities().get(axis)
                == Some(&QuantitySpec::dimensional(Dims::NONE))
            && source.at.axes().get(axis) == Some(&1.0)
            && source.at.axis_quantities().get(axis) == Some(&QuantitySpec::dimensional(Dims::NONE))
    };
    if !owns(phase) || (pure && !owns("source-pure-species")) {
        return Err(AcidBaseError::Invalid(
            "acid/base source must own its aqueous or pure-liquid phase",
        ));
    }
    Ok(())
}
fn consistent_sources<'a>(
    standards: impl IntoIterator<Item = &'a StandardReactionGibbs>,
) -> Result<()> {
    let mut first: Option<&StandardReactionGibbs> = None;
    let mut selected: BTreeMap<&SpeciesId, (&str, &MaterialAnswer)> = BTreeMap::new();
    for standard in standards {
        if let Some(reference) = first {
            if standard.temperature() != reference.temperature()
                || standard.reference_pressure() != reference.reference_pressure()
                || standard.reference_axis() != reference.reference_axis()
            {
                return Err(AcidBaseError::Invalid(
                    "all acid/base references must share temperature pressure and elemental convention",
                ));
            }
        } else {
            first = Some(standard);
        }
        for term in standard.terms() {
            if let Some((pack, answer)) = selected.get(&term.term.species) {
                if *pack != term.pack_id || **answer != term.answer {
                    return Err(AcidBaseError::Invalid(
                        "shared species must retain the identical original selected source",
                    ));
                }
            } else {
                selected.insert(&term.term.species, (&term.pack_id, &term.answer));
            }
        }
    }
    Ok(())
}
fn positive(value: f64, reason: &'static str) -> Result<f64> {
    if value.is_finite() && value > 0.0 {
        Ok(value)
    } else {
        Err(AcidBaseError::Arithmetic(reason))
    }
}
fn nonnegative(value: f64, reason: &'static str) -> Result<f64> {
    if value.is_finite() && value >= 0.0 {
        Ok(value)
    } else {
        Err(AcidBaseError::Invalid(reason))
    }
}

#[cfg(test)]
#[path = "acid_base/tests.rs"]
mod tests;
