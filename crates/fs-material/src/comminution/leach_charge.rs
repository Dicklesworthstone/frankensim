//! Collect a size-resolved stream into explicit, source-retaining mineral amounts.
//!
//! Collection duration converts kg/s to kg; it is not milling residence time or
//! leach contact time. Each mineral's mass fraction is supplied independently of
//! the stream's elemental assays. No mineralogy, liberation, or reactive grain
//! radius is inferred from an assay, mineral name, or bulk particle diameter.
//! A caller may explicitly declare liberated spheres to use class diameter / 2.
//!
//! The original stream, whole-stream elemental accounts, carrier water, source
//! queries, and property receipts remain inspectable. Unallocated dry mass has
//! unknown mineralogy: it is not certified inert, and the original elemental
//! accounts are not assigned to it. Without elemental mineralogy data this
//! adapter cannot certify that mineral allocations reproduce elemental assays.
//! It does not initialize a solution, select kinetics, or perform reactions.

use crate::separation::ParticleStream;
use fs_blake3::ContentHash;
use fs_matdb::{
    MatDbError, MaterialAnswer, NormalizedPack, PropertyKey, Provenance, QueryPoint,
    SelectionPolicy,
};
use fs_qty::chemistry::{ChemistryError, ElementId, MassAmountBasis, SpeciesId};
use fs_qty::{Amount, Dimensionless, Length, Mass, MolarMass, QuantitySpec, Time};
use std::collections::{BTreeMap, BTreeSet};

/// Maximum number of original particle classes in one collection operation.
pub const MAX_COLLECTED_CLASSES: usize = 128;
/// Maximum explicitly declared mineral populations in one collection operation.
pub const MAX_MINERAL_POPULATIONS: usize = 1024;

/// A refused source, allocation, or arithmetic operation; inputs are never repaired.
#[derive(Debug, Clone, PartialEq)]
pub enum LeachChargeError {
    /// Original typed material query or receipt failure.
    Material(MatDbError),
    /// Original typed mass-to-amount conversion failure.
    Chemistry(ChemistryError),
    /// Invalid or unsupported physical/source input.
    Invalid(&'static str),
    /// Population identity used more than once.
    DuplicatePopulation(String),
    /// Allocation refers to a class absent from the original stream.
    UnknownClass(String),
    /// Exact sum of declared binary64 fractions exceeds unity.
    OverallocatedClass(String),
    /// Overflow or erasure of a positive quantity by underflow.
    Arithmetic(&'static str),
}
impl core::fmt::Display for LeachChargeError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Material(e) => write!(f, "leach-charge material source: {e}"),
            Self::Chemistry(e) => write!(f, "leach-charge amount basis: {e}"),
            Self::Invalid(e) | Self::Arithmetic(e) => write!(f, "leach charge: {e}"),
            Self::DuplicatePopulation(id) => write!(f, "duplicate leach population {id:?}"),
            Self::UnknownClass(id) => write!(f, "unknown leach feed class {id:?}"),
            Self::OverallocatedClass(id) => write!(f, "mineral fractions exceed unity in {id:?}"),
        }
    }
}
impl std::error::Error for LeachChargeError {}
impl From<MatDbError> for LeachChargeError {
    fn from(value: MatDbError) -> Self {
        Self::Material(value)
    }
}
impl From<ChemistryError> for LeachChargeError {
    fn from(value: ChemistryError) -> Self {
        Self::Chemistry(value)
    }
}
type Result<T> = core::result::Result<T, LeachChargeError>;

/// Immutable positive molar mass selected through the real material database.
#[derive(Debug, Clone, PartialEq)]
pub struct MolarMassSource {
    pack_id: String,
    pack_hash: ContentHash,
    query: QueryPoint,
    answer: MaterialAnswer,
}
impl MolarMassSource {
    /// Resolve the typed `molar-mass` property, retaining and verifying its receipt.
    ///
    /// The caller remains responsible for the declared species-to-pack mapping;
    /// neither a pack name nor a species label authenticates chemical identity.
    pub fn resolve(
        pack: &NormalizedPack,
        query: &QueryPoint,
        policy: SelectionPolicy,
    ) -> Result<Self> {
        let answer = pack.claims().query_typed(
            &PropertyKey::with_quantity("molar-mass", QuantitySpec::dimensional(MolarMass::DIMS)),
            query,
            policy,
        )?;
        pack.claims().verify_receipt(&answer.receipt)?;
        require(
            answer.evidence.value.value.is_finite() && answer.evidence.value.value > 0.0,
            "source molar mass must be finite and positive",
        )?;
        Ok(Self {
            pack_id: pack.pack_id().to_owned(),
            pack_hash: pack.content_hash(),
            query: query.clone(),
            answer,
        })
    }
    /// Selected molar mass, in kg/mol.
    #[must_use]
    pub fn molar_mass(&self) -> MolarMass {
        MolarMass::new(self.answer.evidence.value.value)
    }
    /// Original normalized pack identifier.
    #[must_use]
    pub fn pack_id(&self) -> &str {
        &self.pack_id
    }
    /// Content hash of the complete source pack.
    #[must_use]
    pub const fn pack_hash(&self) -> ContentHash {
        self.pack_hash
    }
    /// Original typed query conditions.
    #[must_use]
    pub const fn query(&self) -> &QueryPoint {
        &self.query
    }
    /// Original evidence, uncertainty, and verified usage receipt.
    #[must_use]
    pub const fn answer(&self) -> &MaterialAnswer {
        &self.answer
    }
}

/// Explicit interpretation of the mineral radius used by a later reaction model.
#[derive(Debug, Clone, PartialEq)]
pub enum ReactiveRadius {
    /// Independently sourced reactive grain radius; no size-class relation is inferred.
    Sourced {
        radius: Length,
        provenance: Provenance,
    },
    /// Caller explicitly declares liberated spherical grains with class diameter.
    LiberatedSphere { provenance: Provenance },
}
impl ReactiveRadius {
    fn resolve(&self, diameter: Length) -> Result<Length> {
        let (radius, provenance) = match self {
            Self::Sourced { radius, provenance } => (*radius, provenance),
            Self::LiberatedSphere { provenance } => {
                (Length::new(product(diameter.value(), 0.5)?), provenance)
            }
        };
        source(provenance)?;
        require(
            radius.value().is_finite() && radius.value() > 0.0,
            "reactive radius must be finite and positive",
        )?;
        Ok(radius)
    }
}

/// Caller-declared mineral mass fraction of one original dry size class.
///
/// Different population IDs may identify distinct grain-radius populations of
/// the same species in one class; all their fractions spend the same class budget.
#[derive(Debug, Clone, PartialEq)]
pub struct MineralClassAllocation {
    /// Caller identity used by subsequent reaction populations.
    pub population_id: String,
    /// Exact original `ParticleClass` identity.
    pub class_id: String,
    /// Caller-declared solid species, not inferred from a pack or elemental assay.
    pub solid: SpeciesId,
    /// Fraction of this class's total dry mass, in [0, 1].
    pub mass_fraction: Dimensionless,
    /// Artifact-backed mineralogy or explicitly authored scenario evidence.
    pub fraction_source: Provenance,
    /// Actual typed source selection used to convert mineral mass to amount.
    pub molar_mass: MolarMassSource,
    /// Explicit radius evidence or declared liberated-spherical interpretation.
    pub radius: ReactiveRadius,
}

/// Immutable mineral amount and its complete original allocation/source context.
#[derive(Debug, Clone, PartialEq)]
pub struct CollectedMineralPopulation {
    allocation: MineralClassAllocation,
    mass_basis: MassAmountBasis,
    radius: Length,
}
impl CollectedMineralPopulation {
    /// Caller population identity.
    #[must_use]
    pub fn id(&self) -> &str {
        &self.allocation.population_id
    }
    /// Original class identity.
    #[must_use]
    pub fn class_id(&self) -> &str {
        &self.allocation.class_id
    }
    /// Caller-declared solid species.
    #[must_use]
    pub const fn solid(&self) -> &SpeciesId {
        &self.allocation.solid
    }
    /// Collected mineral amount in mol, independent of any later leach duration.
    #[must_use]
    pub const fn amount(&self) -> Amount {
        self.mass_basis.amount()
    }
    /// Collected mineral mass in kg.
    #[must_use]
    pub const fn mass(&self) -> Mass {
        self.mass_basis.mass()
    }
    /// Explicitly resolved reactive grain radius in metres.
    #[must_use]
    pub const fn radius(&self) -> Length {
        self.radius
    }
    /// Original declarations, including fraction and radius provenance.
    #[must_use]
    pub const fn allocation(&self) -> &MineralClassAllocation {
        &self.allocation
    }
    /// Typed mass / molar-mass / amount relationship owned by `fs-qty`.
    #[must_use]
    pub const fn mass_basis(&self) -> MassAmountBasis {
        self.mass_basis
    }
}

/// Whole collected class accounts; assays are independent of mineral allocations.
#[derive(Debug, Clone, PartialEq)]
pub struct CollectedClassAccount {
    class_id: String,
    dry_mass: Mass,
    allocated_mass: Mass,
    unallocated_mass: Mass,
    elements: BTreeMap<ElementId, Mass>,
}
impl CollectedClassAccount {
    /// Original class identity.
    #[must_use]
    pub fn class_id(&self) -> &str {
        &self.class_id
    }
    /// Original dry mass collected in this class.
    #[must_use]
    pub const fn dry_mass(&self) -> Mass {
        self.dry_mass
    }
    /// Sum of explicitly allocated mineral masses (subject to binary64 rounding).
    #[must_use]
    pub const fn allocated_mass(&self) -> Mass {
        self.allocated_mass
    }
    /// Dry mass whose mineralogy remains unallocated; not certified inert mass.
    #[must_use]
    pub const fn unallocated_mass(&self) -> Mass {
        self.unallocated_mass
    }
    /// Original whole-class assay masses, without assignment to any mineral.
    #[must_use]
    pub const fn elements(&self) -> &BTreeMap<ElementId, Mass> {
        &self.elements
    }
}

/// Collected charge with immutable original stream and independent mass accounts.
#[derive(Debug, Clone, PartialEq)]
pub struct CollectedLeachCharge {
    stream: ParticleStream,
    collection_duration: Time,
    populations: Vec<CollectedMineralPopulation>,
    classes: Vec<CollectedClassAccount>,
    dry_mass: Mass,
    unallocated_mass: Mass,
    water_mass: Mass,
    elements: BTreeMap<ElementId, Mass>,
}
impl CollectedLeachCharge {
    /// Collect constant original stream rates over finite nonnegative time.
    ///
    /// Zero time, zero-flow classes, zero fractions, and an empty allocation list
    /// are explicit zero/unallocated cases. Every declaration is still validated.
    /// Fractions are summed with an error-free expansion: an over-allocation is
    /// refused even if ordinary addition rounds its sum to 1. No normalization
    /// is performed. Derived masses round independently to binary64; this is not
    /// a mineral-to-element assay reconciliation or a certified error interval.
    pub fn collect(
        stream: &ParticleStream,
        collection_duration: Time,
        mut allocations: Vec<MineralClassAllocation>,
    ) -> Result<Self> {
        require(
            stream.classes().len() <= MAX_COLLECTED_CLASSES,
            "too many collected size classes",
        )?;
        require(
            allocations.len() <= MAX_MINERAL_POPULATIONS,
            "too many mineral populations",
        )?;
        let duration = collection_duration.value();
        require(
            duration.is_finite() && duration >= 0.0,
            "collection duration must be finite and nonnegative",
        )?;
        allocations.sort_by(|a, b| a.population_id.cmp(&b.population_id));
        let by_class: BTreeMap<_, _> = stream.classes().iter().map(|c| (c.id(), c)).collect();
        let mut ids = BTreeSet::new();
        let mut remaining: BTreeMap<&str, Vec<f64>> =
            by_class.keys().map(|&id| (id, vec![1.0])).collect();
        let mut populations = Vec::with_capacity(allocations.len());
        for allocation in allocations {
            require(
                !allocation.population_id.trim().is_empty()
                    && allocation.population_id.len() <= 256
                    && !allocation.population_id.chars().any(char::is_control),
                "population identity must contain 1 to 256 printable bytes",
            )?;
            if !ids.insert(allocation.population_id.clone()) {
                return Err(LeachChargeError::DuplicatePopulation(
                    allocation.population_id,
                ));
            }
            let class = by_class
                .get(allocation.class_id.as_str())
                .ok_or_else(|| LeachChargeError::UnknownClass(allocation.class_id.clone()))?;
            let fraction = allocation.mass_fraction.value();
            require(
                fraction.is_finite() && (0.0..=1.0).contains(&fraction),
                "mineral mass fraction must lie in [0, 1]",
            )?;
            source(&allocation.fraction_source)?;
            let remainder = remaining
                .get_mut(class.id())
                .expect("validated original class");
            expansion_add(remainder, -fraction)?;
            if remainder.last().is_some_and(|&largest| largest < 0.0) {
                return Err(LeachChargeError::OverallocatedClass(class.id().to_owned()));
            }
            let radius = allocation.radius.resolve(class.diameter())?;
            let mass = Mass::new(product(
                product(class.dry_solids().value(), duration)?,
                fraction,
            )?);
            let mass_basis = MassAmountBasis::from_mass(mass, allocation.molar_mass.molar_mass())?;
            if mass.value() > 0.0 && mass_basis.amount().value() == 0.0 {
                return Err(LeachChargeError::Arithmetic(
                    "positive mineral amount underflows",
                ));
            }
            populations.push(CollectedMineralPopulation {
                allocation,
                mass_basis,
                radius,
            });
        }
        let mut classes = Vec::with_capacity(by_class.len());
        for class in stream.classes() {
            let dry_mass = Mass::new(product(class.dry_solids().value(), duration)?);
            let unallocated_fraction = expansion_value(&remaining[class.id()])?;
            let elements = class
                .elements()
                .iter()
                .map(|(id, flow)| {
                    product(flow.value(), duration).map(|mass| (id.clone(), Mass::new(mass)))
                })
                .collect::<Result<BTreeMap<_, _>>>()?;
            classes.push(CollectedClassAccount {
                class_id: class.id().to_owned(),
                dry_mass,
                allocated_mass: Mass::new(sum(populations
                    .iter()
                    .filter(|p| p.class_id() == class.id())
                    .map(|p| p.mass().value()))?),
                unallocated_mass: Mass::new(product(dry_mass.value(), unallocated_fraction)?),
                elements,
            });
        }
        let elements = classes[0]
            .elements
            .keys()
            .map(|id| {
                sum(classes.iter().map(|c| c.elements[id].value()))
                    .map(|v| (id.clone(), Mass::new(v)))
            })
            .collect::<Result<BTreeMap<_, _>>>()?;
        Ok(Self {
            stream: stream.clone(),
            collection_duration,
            populations,
            dry_mass: Mass::new(sum(classes.iter().map(|c| c.dry_mass.value()))?),
            unallocated_mass: Mass::new(sum(classes.iter().map(|c| c.unallocated_mass.value()))?),
            water_mass: Mass::new(product(stream.water().value(), duration)?),
            classes,
            elements,
        })
    }
    /// Exact immutable original stream, including size coordinates and assay axes.
    #[must_use]
    pub const fn original_stream(&self) -> &ParticleStream {
        &self.stream
    }
    /// Collection interval; no reaction or mill residence interval is implied.
    #[must_use]
    pub const fn collection_duration(&self) -> Time {
        self.collection_duration
    }
    /// Populations in canonical caller-ID order, including explicitly zero amounts.
    #[must_use]
    pub fn populations(&self) -> &[CollectedMineralPopulation] {
        &self.populations
    }
    /// Original class order with separate allocation and whole-class assay accounts.
    #[must_use]
    pub fn class_accounts(&self) -> &[CollectedClassAccount] {
        &self.classes
    }
    /// Total original dry mass in kg.
    #[must_use]
    pub const fn dry_mass(&self) -> Mass {
        self.dry_mass
    }
    /// Total dry mass without a declared mineral population, in kg.
    #[must_use]
    pub const fn unallocated_mass(&self) -> Mass {
        self.unallocated_mass
    }
    /// Original carrier water in kg; no aqueous speciation is inferred.
    #[must_use]
    pub const fn water_mass(&self) -> Mass {
        self.water_mass
    }
    /// Original whole-stream element masses, independent of mineral allocations.
    #[must_use]
    pub const fn elements(&self) -> &BTreeMap<ElementId, Mass> {
        &self.elements
    }
}

fn require(condition: bool, message: &'static str) -> Result<()> {
    if condition {
        Ok(())
    } else {
        Err(LeachChargeError::Invalid(message))
    }
}
fn source(provenance: &Provenance) -> Result<()> {
    require(
        !provenance.source.trim().is_empty()
            && !provenance.license.trim().is_empty()
            && provenance.artifact.is_some(),
        "fraction and radius sources require nonblank source/license and an artifact",
    )
}
fn product(a: f64, b: f64) -> Result<f64> {
    let value = a * b;
    if !value.is_finite() || (a > 0.0 && b > 0.0 && value == 0.0) {
        Err(LeachChargeError::Arithmetic(
            "collected quantity overflows or positive quantity underflows",
        ))
    } else {
        Ok(if value == 0.0 { 0.0 } else { value })
    }
}

// Error-free nonoverlapping expansion of binary64 inputs. FastTwoSum retains
// the rounding residual after ordering by magnitude; the largest part gives
// the exact sign. This keeps a sub-ulp excess over unity visible and preserves
// a small positive remainder even when the allocated sum rounds to unity.
fn expansion_add(parts: &mut Vec<f64>, mut x: f64) -> Result<()> {
    let mut retained = 0;
    for index in 0..parts.len() {
        let mut y = parts[index];
        if x.abs() < y.abs() {
            core::mem::swap(&mut x, &mut y);
        }
        let high = x + y;
        if !high.is_finite() {
            return Err(LeachChargeError::Arithmetic(
                "collected mass total overflows",
            ));
        }
        let low = y - (high - x);
        if low != 0.0 {
            parts[retained] = low;
            retained += 1;
        }
        x = high;
    }
    parts.truncate(retained);
    if x != 0.0 {
        parts.push(x);
    }
    Ok(())
}
fn expansion_value(parts: &[f64]) -> Result<f64> {
    let value: f64 = parts.iter().sum();
    if !value.is_finite() || (!parts.is_empty() && value == 0.0) {
        Err(LeachChargeError::Arithmetic(
            "collected mass sum is not representable",
        ))
    } else {
        Ok(value)
    }
}
fn sum(values: impl Iterator<Item = f64>) -> Result<f64> {
    let mut parts = Vec::new();
    for value in values {
        expansion_add(&mut parts, value)?;
    }
    expansion_value(&parts)
}

#[cfg(test)]
#[path = "leach_charge/tests.rs"]
mod tests;
