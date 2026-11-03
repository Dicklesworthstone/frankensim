//! Phase-independent reaction amounts and source-bound standard reaction heat.
//!
//! Exact element/charge balance stays in `fs_qty::chemistry`. These closures
//! do not select a reaction pathway, predict conversion, integrate kinetics,
//! or supply aqueous activities. In particular, a mineral is never admitted
//! through the ideal-gas species/EOS lane to perform this bookkeeping.

pub mod minerals;
/// Shared-inventory stages of independently admitted, exactly conserved reactions.
pub mod network;

use std::collections::BTreeMap;

use fs_matdb::{
    MatDbError, MaterialAnswer, NormalizedPack, PropertyKey, QueryPoint, SelectionPolicy,
};
use fs_qty::{Amount, Dims, ElectricCharge, Energy, Mass, MolarMass, QuantitySpec};

use crate::{
    ChargeVector, ChemistryError, ConservationCertificate, ElementalMatrix, MassAmountBasis,
    MolarEnergyQuantityV1, ReactionId, SpeciesId, StoichiometricMatrix, verify_conservation,
};

/// Faraday constant [C/mol], from the exact SI elementary charge and Avogadro
/// constant. Reference: https://physics.nist.gov/cuu/Constants/Table/allascii.txt
pub const FARADAY_C_PER_MOL: f64 = 1.602_176_634e-19 * 6.022_140_76e23;
const MAX_SPECIES: usize = 256;
const MAX_REACTIONS: usize = 256;
const MAX_EXACT_COEFFICIENT: i128 = 1_i128 << 53;

/// A malformed equation, missing source, or unrepresentable calculation.
#[derive(Debug, Clone, PartialEq)]
pub enum ReactionBalanceError {
    /// Exact chemistry admission failed in its original owner.
    Chemistry(ChemistryError),
    /// Property selection or source-receipt verification failed.
    Material(MatDbError),
    /// Invalid physical, reference-state or resource input.
    Invalid(&'static str),
    /// A required species is absent, duplicated, or not a participant.
    Species {
        /// Species whose input cannot be associated with this equation.
        species: SpeciesId,
        /// Required repair at the caller boundary.
        reason: &'static str,
    },
    /// A declared forward extent would consume more material than is present.
    ExtentExceedsInventory {
        /// Requested amount of the integer reaction as written [mol].
        requested: Amount,
        /// Conservatively admitted material-availability bound [mol].
        maximum: Amount,
    },
    /// A finite nonzero result cannot be represented without range loss.
    Arithmetic(&'static str),
}

impl core::fmt::Display for ReactionBalanceError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "reaction balance: {self:?}")
    }
}
impl std::error::Error for ReactionBalanceError {}
impl From<ChemistryError> for ReactionBalanceError {
    fn from(error: ChemistryError) -> Self {
        Self::Chemistry(error)
    }
}
impl From<MatDbError> for ReactionBalanceError {
    fn from(error: MatDbError) -> Self {
        Self::Material(error)
    }
}
type Result<T> = core::result::Result<T, ReactionBalanceError>;

/// One participant; negative coefficients consume, positive coefficients form.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ReactionTerm {
    /// Canonical species identity, including phase when the source declares it.
    pub species: SpeciesId,
    /// Signed coefficient in the integer equation as written.
    pub coefficient: i128,
}

/// Immutable, exactly conserved reaction column with canonical species order.
#[derive(Debug, Clone)]
pub struct ConservedReaction {
    id: ReactionId,
    terms: Vec<ReactionTerm>,
    definitions: Vec<SpeciesDefinition>,
    conservation: ConservationCertificate,
}

// Aligned with active terms; absent/zero element rows have the same identity.
// Retained from fs-qty's admitted artifacts, not reparsed from species labels.
#[derive(Debug, Clone, PartialEq, Eq)]
struct SpeciesDefinition {
    elements: Vec<(crate::ElementId, i128)>,
    charge: i128,
}

impl ConservedReaction {
    /// Select one reaction from supplied chemical bookkeeping artifacts.
    ///
    /// All columns must conserve elements and charge. The selected column must
    /// have both reactants and products and coefficients exactly convertible
    /// to binary64. Zero-coefficient spectators need no property sources.
    pub fn new(
        elemental: &ElementalMatrix,
        stoichiometric: &StoichiometricMatrix,
        charge: &ChargeVector,
        id: ReactionId,
    ) -> Result<Self> {
        if elemental.row_count() > MAX_SPECIES
            || stoichiometric.row_count() > MAX_SPECIES
            || stoichiometric.column_count() > MAX_REACTIONS
        {
            return Err(ReactionBalanceError::Invalid(
                "reaction layout exceeds 256-axis bound",
            ));
        }
        let conservation = verify_conservation(elemental, stoichiometric, charge)?;
        let column = stoichiometric
            .reactions()
            .binary_search(&id)
            .map_err(|_| ReactionBalanceError::Invalid("reaction id is absent from matrix"))?;
        let mut terms = Vec::new();
        let mut definitions = Vec::new();
        for (row, species) in stoichiometric.species().iter().enumerate() {
            let coefficient =
                stoichiometric
                    .get(row, column)
                    .ok_or(ReactionBalanceError::Invalid(
                        "invalid stoichiometric matrix shape",
                    ))?;
            if !(-MAX_EXACT_COEFFICIENT..=MAX_EXACT_COEFFICIENT).contains(&coefficient) {
                return Err(ReactionBalanceError::Invalid(
                    "coefficient exceeds exact binary64 integer range",
                ));
            }
            if coefficient != 0 {
                let mut elements = Vec::new();
                for (element_row, element) in elemental.elements().iter().enumerate() {
                    let count =
                        elemental
                            .get(element_row, row)
                            .ok_or(ReactionBalanceError::Invalid(
                                "invalid elemental matrix shape",
                            ))?;
                    if count != 0 {
                        elements.push((element.clone(), count));
                    }
                }
                definitions.push(SpeciesDefinition {
                    elements,
                    charge: charge
                        .get(row)
                        .ok_or(ReactionBalanceError::Invalid("invalid charge vector shape"))?,
                });
                terms.push(ReactionTerm {
                    species: species.clone(),
                    coefficient,
                });
            }
        }
        if !terms.iter().any(|term| term.coefficient < 0)
            || !terms.iter().any(|term| term.coefficient > 0)
        {
            return Err(ReactionBalanceError::Invalid(
                "reaction must consume and form species",
            ));
        }
        Ok(Self {
            id,
            terms,
            definitions,
            conservation,
        })
    }

    /// Original selected reaction identity.
    #[must_use]
    pub fn id(&self) -> &ReactionId {
        &self.id
    }
    /// Active nonzero participants in canonical species order.
    #[must_use]
    pub fn terms(&self) -> &[ReactionTerm] {
        &self.terms
    }
    /// Exact element/charge certificate for the supplied matrix artifacts.
    #[must_use]
    pub fn conservation(&self) -> &ConservationCertificate {
        &self.conservation
    }

    /// Nonnegative amounts consumed/formed at a caller-supplied extent [mol].
    /// The coefficient sign retains direction; this does not predict extent.
    pub fn amounts(&self, extent: Amount) -> Result<Vec<ReactionAmount>> {
        nonnegative(extent.value(), "extent must be finite and nonnegative")?;
        self.terms
            .iter()
            .map(|term| {
                let amount = product(
                    extent.value(),
                    (term.coefficient as f64).abs(),
                    "species amount",
                )?;
                Ok(ReactionAmount {
                    term: term.clone(),
                    amount: Amount::new(amount),
                })
            })
            .collect()
    }

    /// Extent corresponding to the mass of ONE reactant actually reacted.
    /// Feed purity and conversion must be applied explicitly by the caller.
    pub fn extent_from_reacted_mass(
        &self,
        species: &SpeciesId,
        mass: Mass,
        molar_mass: MolarMass,
    ) -> Result<Amount> {
        let term = self
            .terms
            .iter()
            .find(|term| &term.species == species && term.coefficient < 0)
            .ok_or_else(|| ReactionBalanceError::Species {
                species: species.clone(),
                reason: "expected a consumed species",
            })?;
        let basis = MassAmountBasis::from_mass(mass, molar_mass)?;
        let amount = basis.amount().value();
        if mass.value() > 0.0 && amount == 0.0 {
            return Err(ReactionBalanceError::Arithmetic("reacted amount underflow"));
        }
        let extent = amount / -(term.coefficient as f64);
        representable(extent, amount > 0.0, "reaction extent")?;
        Ok(Amount::new(extent))
    }

    /// Stoichiometric upper bound from all and only the reactant inventories.
    /// This is a material-availability bound, not an equilibrium or recovery.
    /// A nearest-rounded ratio is moved down one binary64 step when its exact
    /// coefficient product exceeds stock, so the returned extent is feasible.
    pub fn limiting_extent(&self, inventories: &[(SpeciesId, Amount)]) -> Result<Amount> {
        let reactants: Vec<_> = self
            .terms
            .iter()
            .filter(|term| term.coefficient < 0)
            .collect();
        if inventories.len() != reactants.len() {
            return Err(ReactionBalanceError::Invalid(
                "supply exactly one inventory per reactant",
            ));
        }
        let mut by_species = BTreeMap::new();
        for (species, amount) in inventories {
            nonnegative(amount.value(), "inventory must be finite and nonnegative")?;
            if by_species.insert(species, amount.value()).is_some() {
                return Err(ReactionBalanceError::Species {
                    species: species.clone(),
                    reason: "duplicate inventory",
                });
            }
        }
        let mut limit = f64::INFINITY;
        for term in reactants {
            let amount =
                by_species
                    .get(&term.species)
                    .ok_or_else(|| ReactionBalanceError::Species {
                        species: term.species.clone(),
                        reason: "missing reactant inventory",
                    })?;
            let scale = -(term.coefficient as f64);
            let mut extent = amount / scale;
            representable(extent, *amount > 0.0, "inventory extent")?;
            // Ordinary multiplication can round an overspend back onto stock.
            // FMA preserves its sign. The coefficient is an exact integer, so
            // a nonzero residual cannot be smaller than the least subnormal.
            if extent.mul_add(scale, -*amount) > 0.0 {
                extent = extent.next_down();
            }
            representable(extent, *amount > 0.0, "inventory extent")?;
            limit = limit.min(extent);
        }
        Ok(Amount::new(limit))
    }

    /// Advance a complete inventory by a declared nonnegative forward extent.
    ///
    /// Supply exactly one amount for every active species, including explicit
    /// zero for absent products. Spectators are managed outside this reaction.
    /// Input order is arbitrary; the returned inventory is canonical. Inputs
    /// are never mutated, and an excessive extent is refused rather than
    /// clipped. Element/charge conservation is exact for the equation; final
    /// amounts are binary64 evaluations of `initial + coefficient * extent`.
    /// Nonzero changes below inventory resolution refuse instead of vanishing.
    pub fn advance_inventory(
        &self,
        inventories: &[(SpeciesId, Amount)],
        extent: Amount,
    ) -> Result<Vec<(SpeciesId, Amount)>> {
        nonnegative(extent.value(), "extent must be finite and nonnegative")?;
        if inventories.len() != self.terms.len() {
            return Err(ReactionBalanceError::Invalid(
                "supply exactly one inventory per active species",
            ));
        }
        let mut by_species = BTreeMap::new();
        for (species, amount) in inventories {
            nonnegative(amount.value(), "inventory must be finite and nonnegative")?;
            if by_species.insert(species, amount.value()).is_some() {
                return Err(ReactionBalanceError::Species {
                    species: species.clone(),
                    reason: "duplicate inventory",
                });
            }
        }
        let mut initial = Vec::with_capacity(self.terms.len());
        let mut reactants = Vec::new();
        for term in &self.terms {
            let amount =
                *by_species
                    .get(&term.species)
                    .ok_or_else(|| ReactionBalanceError::Species {
                        species: term.species.clone(),
                        reason: "missing active-species inventory",
                    })?;
            initial.push((term.species.clone(), Amount::new(amount)));
            if term.coefficient < 0 {
                reactants.push((term.species.clone(), Amount::new(amount)));
            }
        }
        if extent.value() == 0.0 {
            return Ok(initial);
        }
        let maximum = self.limiting_extent(&reactants)?;
        if extent.value() > maximum.value() {
            return Err(ReactionBalanceError::ExtentExceedsInventory {
                requested: extent,
                maximum,
            });
        }
        let mut result = Vec::with_capacity(initial.len());
        for ((species, before), flow) in initial.into_iter().zip(self.amounts(extent)?) {
            let after = (flow.term.coefficient as f64).mul_add(extent.value(), before.value());
            representable(after, false, "final species inventory")?;
            if after < 0.0 {
                return Err(ReactionBalanceError::Arithmetic("negative final inventory"));
            }
            if after == before.value() {
                return Err(ReactionBalanceError::Arithmetic(
                    "reaction advance is below inventory resolution",
                ));
            }
            result.push((species, Amount::new(after)));
        }
        Ok(result)
    }

    /// Sum standard formation enthalpies in canonical species order.
    ///
    /// Every active species must have exactly one source at the same typed
    /// temperature, pressure and elemental-reference flag. Phase is retained
    /// in each query. The species-to-pack association is explicitly supplied
    /// by the caller; ordinary property packs do not authenticate formulas.
    /// Returns products minus reactants; negative means exothermic at constant
    /// pressure. No heat of mixing/dilution, kinetics or uncertainty propagation
    /// is inferred. Original uncertainty/evidence and receipts remain available.
    pub fn standard_enthalpy_change(
        &self,
        sources: &[FormationEnthalpySource<'_>],
    ) -> Result<ReactionEnthalpy> {
        if sources.len() != self.terms.len() {
            return Err(ReactionBalanceError::Invalid(
                "supply one formation source per active species",
            ));
        }
        let mut by_species = BTreeMap::new();
        for source in sources {
            source.validate_reference()?;
            if by_species.insert(source.species, source).is_some() {
                return Err(ReactionBalanceError::Species {
                    species: source.species.clone(),
                    reason: "duplicate formation source",
                });
            }
        }
        let first = sources
            .first()
            .ok_or(ReactionBalanceError::Invalid("no formation sources"))?;
        let mut terms = Vec::with_capacity(self.terms.len());
        let mut sum = 0.0;
        let mut correction = 0.0;
        for term in &self.terms {
            let source =
                by_species
                    .get(&term.species)
                    .ok_or_else(|| ReactionBalanceError::Species {
                        species: term.species.clone(),
                        reason: "missing formation source",
                    })?;
            if source.reference_axis != first.reference_axis
                || source.at.axes().get("temperature") != first.at.axes().get("temperature")
                || source.at.axes().get("pressure") != first.at.axes().get("pressure")
            {
                return Err(ReactionBalanceError::Invalid(
                    "formation reference states differ",
                ));
            }
            let answer = source.pack.claims().query_typed(
                &PropertyKey::with_quantity(
                    "standard-molar-formation-enthalpy",
                    QuantitySpec::dimensional(Dims([2, 1, -2, 0, 0, -1])),
                ),
                source.at,
                SelectionPolicy::SingleClaimOnly,
            )?;
            source.pack.claims().verify_receipt(&answer.receipt)?;
            source.validate_selected(&answer)?;
            let value = product(
                term.coefficient as f64,
                answer.evidence.value.value,
                "reaction enthalpy term",
            )?;
            // Neumaier reduction keeps near-cancelling formation references
            // without making an experimental/numerical error-bound claim.
            let next = sum + value;
            representable(next, false, "reaction enthalpy sum")?;
            correction += if sum.abs() >= value.abs() {
                (sum - next) + value
            } else {
                (value - next) + sum
            };
            sum = next;
            terms.push(FormationEnthalpyTerm {
                term: term.clone(),
                pack_id: source.pack.pack_id().to_owned(),
                answer,
            });
        }
        let value = sum + correction;
        representable(value, false, "reaction enthalpy")?;
        Ok(ReactionEnthalpy {
            molar: MolarEnergyQuantityV1::new(value),
            terms,
        })
    }
}

/// Absolute molar amount; its term's sign identifies consumption or production.
#[derive(Debug, Clone)]
pub struct ReactionAmount {
    /// Species and signed coefficient specifying the direction of transfer.
    pub term: ReactionTerm,
    /// Nonnegative absolute amount consumed or formed [mol].
    pub amount: Amount,
}

/// Source query for one explicitly associated species/phase. The reference
/// axis is a common dimensionless flag such as `source-elemental-reference-janaf`.
pub struct FormationEnthalpySource<'a> {
    /// Caller-declared chemical identity associated with the offered pack.
    pub species: &'a SpeciesId,
    /// Original normalized property source; no gas-species association is used.
    pub pack: &'a NormalizedPack,
    /// Complete typed reference-state and phase coordinates.
    pub at: &'a QueryPoint,
    /// Declared `source-elemental-reference-*` flag shared by every source.
    pub reference_axis: &'a str,
}

impl FormationEnthalpySource<'_> {
    fn validate_reference(&self) -> Result<()> {
        use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
        for (axis, kind) in [
            ("temperature", QuantityKind::AbsoluteTemperature),
            ("pressure", QuantityKind::Pressure),
        ] {
            let quantity = QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static));
            if self.at.axis_quantities().get(axis) != Some(&quantity)
                || !self.at.axes().get(axis).is_some_and(|value| *value > 0.0)
            {
                return Err(ReactionBalanceError::Invalid(
                    "formation source needs positive typed temperature and pressure",
                ));
            }
        }
        if !self
            .reference_axis
            .starts_with("source-elemental-reference-")
            || self.reference_axis.len() == "source-elemental-reference-".len()
            || self.at.axes().get(self.reference_axis) != Some(&1.0)
            || self.at.axis_quantities().get(self.reference_axis)
                != Some(&QuantitySpec::dimensional(Dims::NONE))
        {
            return Err(ReactionBalanceError::Invalid(
                "formation source needs explicit common elemental reference",
            ));
        }
        Ok(())
    }

    fn validate_selected(&self, answer: &MaterialAnswer) -> Result<()> {
        let claim = self.pack.claims().claim(answer.receipt.selected).ok_or(
            ReactionBalanceError::Invalid("selected formation source is absent"),
        )?;
        // Extra caller coordinates are not source authority. Require the
        // selected claim itself to declare each common reference coordinate.
        for axis in [
            "temperature",
            "pressure",
            self.reference_axis,
            "source-standard-state",
        ] {
            if claim.validity.bound(axis).is_none()
                || claim.validity.axis_quantities().get(axis) != self.at.axis_quantities().get(axis)
            {
                return Err(ReactionBalanceError::Invalid(
                    "formation reference is not declared by selected source",
                ));
            }
        }
        for axis in [self.reference_axis, "source-standard-state"] {
            if claim.validity.bound(axis) != Some((1.0, 1.0))
                || self.at.axes().get(axis) != Some(&1.0)
                || claim.validity.axis_quantities().get(axis)
                    != Some(&QuantitySpec::dimensional(Dims::NONE))
            {
                return Err(ReactionBalanceError::Invalid(
                    "source must explicitly fix its standard elemental reference",
                ));
            }
        }
        let phases = [
            "source-phase-solid",
            "source-phase-liquid",
            "source-phase-aqueous",
            "source-phase-gas",
        ];
        let declared = phases
            .iter()
            .filter(|axis| {
                claim.validity.bound(axis) == Some((1.0, 1.0))
                    && claim.validity.axis_quantities().get(**axis)
                        == Some(&QuantitySpec::dimensional(Dims::NONE))
                    && self.at.axes().get(**axis) == Some(&1.0)
            })
            .count();
        if declared != 1 {
            return Err(ReactionBalanceError::Invalid(
                "source must declare exactly one standard-state phase",
            ));
        }
        Ok(())
    }
}

/// Source evidence for one load-bearing term. No new authority is minted.
#[derive(Debug, Clone)]
pub struct FormationEnthalpyTerm {
    /// Active participant and signed integer coefficient.
    pub term: ReactionTerm,
    /// Source pack identity retained alongside its original property receipt.
    pub pack_id: String,
    /// Original evidence sample and verified usage receipt, including uncertainty.
    pub answer: MaterialAnswer,
}

/// Standard reaction enthalpy and its complete active source dependency set.
#[derive(Debug, Clone)]
pub struct ReactionEnthalpy {
    molar: MolarEnergyQuantityV1,
    terms: Vec<FormationEnthalpyTerm>,
}
impl ReactionEnthalpy {
    /// Standard enthalpy change [J/mol of reaction extent], products minus reactants.
    #[must_use]
    pub fn molar(&self) -> MolarEnergyQuantityV1 {
        self.molar
    }
    /// Complete, canonical set of load-bearing source samples and receipts.
    #[must_use]
    pub fn terms(&self) -> &[FormationEnthalpyTerm] {
        &self.terms
    }
    /// Enthalpy change [J] for nonnegative extent; do not add this twice to a
    /// total-energy formulation already containing chemical formation energy.
    pub fn at_extent(&self, extent: Amount) -> Result<Energy> {
        nonnegative(extent.value(), "extent must be finite and nonnegative")?;
        Ok(Energy::new(product(
            self.molar.value(),
            extent.value(),
            "reaction heat",
        )?))
    }
}

/// Faraday product yield for integrated charge. `electrons_per_formula` is the
/// positive electron count per deposited formula unit; efficiency is explicit
/// in [0,1]. For Cu(II) -> Cu use 2. No voltage, energy, kinetics, current
/// distribution or current-efficiency prediction follows from this identity.
pub fn faraday_yield(
    charge: ElectricCharge,
    molar_mass: MolarMass,
    electrons_per_formula: u32,
    current_efficiency: f64,
) -> Result<MassAmountBasis> {
    nonnegative(charge.value(), "charge must be finite and nonnegative")?;
    if electrons_per_formula == 0
        || !current_efficiency.is_finite()
        || !(0.0..=1.0).contains(&current_efficiency)
    {
        return Err(ReactionBalanceError::Invalid(
            "electron count must be positive and efficiency in [0,1]",
        ));
    }
    let effective = product(charge.value(), current_efficiency, "effective charge")?;
    let amount = effective / (FARADAY_C_PER_MOL * f64::from(electrons_per_formula));
    representable(amount, effective != 0.0, "Faraday amount")?;
    let basis = MassAmountBasis::from_amount(Amount::new(amount), molar_mass)?;
    representable(basis.mass().value(), amount != 0.0, "Faraday mass")?;
    Ok(basis)
}

fn nonnegative(value: f64, reason: &'static str) -> Result<()> {
    if !value.is_finite() || value < 0.0 {
        return Err(ReactionBalanceError::Invalid(reason));
    }
    Ok(())
}
fn representable(value: f64, nonzero: bool, operation: &'static str) -> Result<()> {
    if !value.is_finite() || (nonzero && value == 0.0) {
        return Err(ReactionBalanceError::Arithmetic(operation));
    }
    Ok(())
}
fn product(a: f64, b: f64, operation: &'static str) -> Result<f64> {
    let value = a * b;
    representable(value, a != 0.0 && b != 0.0, operation)?;
    Ok(value)
}

#[cfg(test)]
mod tests;
