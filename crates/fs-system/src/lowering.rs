//! Owner-extensible lowering registry: from canonical assembly descriptions to
//! sealed kernel [`System`] graphs, static [`ExecutionPlan`]s, and auditable
//! [`LoweringReceipt`]s (bead `frankensim-rc-1109-root-3gc12.11.2`).
//!
//! # Architecture
//!
//! Rather than a monolithic `if/else` dispatcher over device slots, each
//! physics owner registers constructors keyed by canonical structural tuples:
//!
//! - **Parts**: `(DraftPartKind, BTreeSet<PhysicsDomain>, FidelityRungKind)` ->
//!   [`PartConstructor`]
//! - **Interfaces**: `(InterfaceKindTag, card_family)` -> [`InterfaceConstructor`]
//! - **Joints**: [`JointLoweringKey`] -> [`JointConstructor`]
//! - **Gear meshes**: [`GearMeshLoweringKey`] -> [`GearMeshConstructor`]
//! - **Acoustic families**: [`AcousticFamilyKey`] -> [`AcousticFamilyConstructor`]
//!
//! Lowering resolves geometry fingerprints and `fs-matdb` material/interface
//! cards (collecting real [`PropertyUsageReceipt`]s), populates one
//! [`SystemBuilder`], finalizes shared port junctions, runs kernel admission
//! (`SystemBuilder::seal`), hands the sealed graph to the structural planner
//! (`System::execution_plan_for_dt` + `System::causalization_audit`), and
//! verifies receipt completeness (no parameter without an explicit source).

#![deny(unsafe_code)]

use std::collections::{BTreeMap, BTreeSet};
use std::f64::consts::PI;
use std::fmt;

use fs_blake3::{ContentHash, hash_domain};
use fs_matdb::{
    ClaimSet, InterfaceSystemCard, MatDbError, MaterialCard, PropertyKey, PropertySample,
    PropertyUsageReceipt, Provenance, QueryPoint, SelectionPolicy,
};
use fs_ports::PortKind;
use fs_ports::component::{Carrier, PortRole};
use fs_qty::semantic::QuantityKind::{AbsoluteTemperature as AbsTemp, ThermalConductivity};
use fs_qty::semantic::{SemanticType, ValueForm};
use fs_qty::{Dims, QtyAny, QuantitySpec};
use fs_scenario::assembly::draft::{
    DraftAssembly, DraftGeometrySource, DraftInterfaceDecl, DraftInterfaceKind, DraftPartDecl,
    DraftPartKind, FidelityRung, GearMeshDecl, JointDecl, JointKind, PhysicsDomain, PrimitiveShape,
    SourceDecl, SourceTarget, SpiralSpringGenerator, StepPolicy, ThirdBodyDecl, TimeHistory,
};
use fs_scenario::assembly::{GeometrySource, MaterialRef};
use fs_scenario::entity::{Entity, EntityCatalog, EntityId, EntityKind};

use crate::elements::{
    Conductance, Damper, ElementError, FlowSource, Force, HeatCapacity, Heater, Mass, Motion,
    SourceAccount, Spring, Storage, Waveform,
};
use crate::{
    AdmissionError, CausalizationAudit, ExecutionPlan, HeatSink, NodeId, PortRef, RouteTemperature,
    System, SystemBuilder, heat_sink,
};

/// Owner constructors: thin adapters over the physics owners' public
/// Components, registered by [`owners::register_owner_defaults`].
#[cfg(feature = "slices")]
pub mod owners;

/// Schema identifier for serialized [`LoweringReceipt`] payloads stored in
/// `fs-ledger`.
pub const LOWERING_RECEIPT_SCHEMA: &str = "frankensim.lowering_receipt.v1";

/// Dispatch family of a card-less [`DraftInterfaceKind::HeatRoute`].
pub const HEAT_ROUTE_FAMILY: &str = "heat-route";
/// Dispatch family of a card-less [`DraftInterfaceKind::ShearFilm`].
pub const SHEAR_FILM_FAMILY: &str = "shear-film";
/// Dispatch family of a card-less [`DraftInterfaceKind::TemperatureSource`].
pub const TEMPERATURE_SOURCE_FAMILY: &str = "temperature-source";
/// Dispatch family of a card-less [`DraftInterfaceKind::PistonFace`].
pub const PISTON_FACE_FAMILY: &str = "piston-face";

/// BLAKE3 domain for parametric geometry fingerprints minted during lowering.
pub const LOWERING_GEOMETRY_HASH_DOMAIN: &str = "org.frankensim.fs-system.lowering.geometry.v1";

/// BLAKE3 domain for canonical [`LoweringReceipt`] content hashes.
pub const LOWERING_RECEIPT_HASH_DOMAIN: &str = "org.frankensim.fs-system.lowering.receipt.v1";

/// BLAKE3 domain for [`PowerPortSchemaCard`] content identities.
pub const POWER_PORT_SCHEMA_HASH_DOMAIN: &str =
    "org.frankensim.fs-system.lowering.power-port-schema.v1";

const MASS_DIMS: Dims = Dims([0, 1, 0, 0, 0, 0]);
const INERTIA_DIMS: Dims = Dims([2, 1, 0, 0, 0, 0]);
const LENGTH_DIMS: Dims = Dims([1, 0, 0, 0, 0, 0]);
const VOLUME_DIMS: Dims = Dims([3, 0, 0, 0, 0, 0]);
const STIFFNESS_DIMS: Dims = Dims([0, 1, -2, 0, 0, 0]);
const TORSIONAL_STIFFNESS_DIMS: Dims = Dims([2, 1, -2, 0, 0, 0]);
const DAMPING_DIMS: Dims = Dims([0, 1, -1, 0, 0, 0]);
const TORSIONAL_DAMPING_DIMS: Dims = Dims([2, 1, -1, 0, 0, 0]);
const HEAT_CAPACITY_DIMS: Dims = Dims([2, 1, -2, -1, 0, 0]);
const CONDUCTANCE_DIMS: Dims = Dims([2, 1, -3, -1, 0, 0]);
const THERMAL_CONDUCTIVITY_DIMS: Dims = Dims([1, 1, -3, -1, 0, 0]);
const AREA_DIMS: Dims = Dims([2, 0, 0, 0, 0, 0]);
const CONTACT_RESISTANCE_DIMS: Dims = Dims([0, -1, 3, 1, 0, 0]);
const TEMPERATURE_DIMS: Dims = Dims([0, 0, 0, 1, 0, 0]);
const POWER_DIMS: Dims = Dims([2, 1, -3, 0, 0, 0]);
const TIME_DIMS: Dims = Dims([0, 0, 1, 0, 0, 0]);
const FREQUENCY_DIMS: Dims = Dims([0, 0, -1, 0, 0, 0]);
const VELOCITY_DIMS: Dims = Dims([1, 0, -1, 0, 0, 0]);

// ============================================================================
// Registry Keys
// ============================================================================

/// Discriminant of a [`FidelityRung`] used as the third element of a
/// [`PartLoweringKey`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum FidelityRungKind {
    /// Lumped-parameter ODE/DAE model.
    Lumped,
    /// Truncated modal basis.
    Modal,
    /// Reduced-order state-space basis.
    ReducedOrder,
    /// 1D finite-volume / finite-cell discretisation.
    OneDField,
    /// Finite-element discretisation.
    FiniteElement,
    /// Lattice-Boltzmann discretisation.
    LatticeBoltzmann,
    /// Boundary-element discretisation.
    BoundaryElement,
}

impl FidelityRungKind {
    /// Stable canonical label for diagnostics and receipts.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Lumped => "lumped",
            Self::Modal => "modal",
            Self::ReducedOrder => "reduced_order",
            Self::OneDField => "one_d_field",
            Self::FiniteElement => "finite_element",
            Self::LatticeBoltzmann => "lattice_boltzmann",
            Self::BoundaryElement => "boundary_element",
        }
    }

    /// Relative fidelity rank from cheapest (`0`) to highest-resolution (`6`).
    #[must_use]
    pub const fn rank(self) -> u8 {
        match self {
            Self::Lumped => 0,
            Self::Modal => 1,
            Self::ReducedOrder => 2,
            Self::OneDField => 3,
            Self::FiniteElement => 4,
            Self::LatticeBoltzmann => 5,
            Self::BoundaryElement => 6,
        }
    }
}

impl From<&FidelityRung> for FidelityRungKind {
    fn from(rung: &FidelityRung) -> Self {
        match rung {
            FidelityRung::Lumped => Self::Lumped,
            FidelityRung::Modal { .. } => Self::Modal,
            FidelityRung::ReducedOrder { .. } => Self::ReducedOrder,
            FidelityRung::OneDField { .. } => Self::OneDField,
            FidelityRung::FiniteElement { .. } => Self::FiniteElement,
            FidelityRung::LatticeBoltzmann { .. } => Self::LatticeBoltzmann,
            FidelityRung::BoundaryElement { .. } => Self::BoundaryElement,
        }
    }
}

impl fmt::Display for FidelityRungKind {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_str())
    }
}

/// Strategy used to pick one rung from a part's ordered `FidelityPolicy::rungs`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum FidelitySelectionPolicy {
    /// Pick the first rung in the draft's ordered preference ladder that has a
    /// registered constructor.
    #[default]
    FirstRegistered,
    /// Pick the highest-rank registered rung in the ladder.
    HighestRegistered,
    /// Pick the lowest-rank (cheapest) registered rung in the ladder.
    LowestRegistered,
}

impl FidelitySelectionPolicy {
    /// Stable canonical tag for receipts.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::FirstRegistered => "first_registered",
            Self::HighestRegistered => "highest_registered",
            Self::LowestRegistered => "lowest_registered",
        }
    }
}

/// Canonical dispatch key for part constructors:
/// `(part kind, physics set, fidelity rung kind)`.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct PartLoweringKey {
    /// Physical part class.
    pub kind: DraftPartKind,
    /// Active physics domain set.
    pub physics: BTreeSet<PhysicsDomain>,
    /// Fidelity rung class.
    pub rung: FidelityRungKind,
}

impl PartLoweringKey {
    /// Construct a part lowering key from a kind, iterator of physics domains,
    /// and fidelity rung class.
    #[must_use]
    pub fn new(
        kind: DraftPartKind,
        physics: impl IntoIterator<Item = PhysicsDomain>,
        rung: FidelityRungKind,
    ) -> Self {
        Self {
            kind,
            physics: physics.into_iter().collect(),
            rung,
        }
    }
}

impl fmt::Display for PartLoweringKey {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "({}, {{", self.kind.as_str())?;
        for (i, domain) in self.physics.iter().enumerate() {
            if i > 0 {
                f.write_str(", ")?;
            }
            f.write_str(domain.as_str())?;
        }
        write!(f, "}}, {})", self.rung)
    }
}

/// Discriminant of a [`DraftInterfaceKind`] used in [`InterfaceLoweringKey`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum InterfaceKindTag {
    /// Typed power-conjugate port schema.
    PowerPort,
    /// Mechanical / tribological contact interface.
    MechanicalContact,
    /// Bonded / adhesive interface.
    Glue,
    /// Thermal interface conductance / contact resistance.
    ThermalContact,
    /// Electrical contact interface.
    ElectricalContact,
    /// Card-less heat route.
    HeatRoute,
    /// Card-less temperature source.
    TemperatureSource,
    /// Card-less mechanical-to-fluid piston face.
    PistonFace,
}

impl InterfaceKindTag {
    /// Stable canonical label.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::PowerPort => "power_port",
            Self::MechanicalContact => "mechanical_contact",
            Self::Glue => "glue",
            Self::ThermalContact => "thermal_contact",
            Self::ElectricalContact => "electrical_contact",
            Self::HeatRoute => "heat_route",
            Self::TemperatureSource => "temperature_source",
            Self::PistonFace => "piston_face",
        }
    }
}

impl From<&DraftInterfaceKind> for InterfaceKindTag {
    fn from(kind: &DraftInterfaceKind) -> Self {
        match kind {
            DraftInterfaceKind::PowerPort { .. } => Self::PowerPort,
            DraftInterfaceKind::MechanicalContact { .. } | DraftInterfaceKind::ShearFilm { .. } => {
                Self::MechanicalContact
            }
            DraftInterfaceKind::Glue { .. } => Self::Glue,
            DraftInterfaceKind::ThermalContact { .. } => Self::ThermalContact,
            DraftInterfaceKind::ElectricalContact { .. } => Self::ElectricalContact,
            DraftInterfaceKind::HeatRoute { .. } => Self::HeatRoute,
            DraftInterfaceKind::TemperatureSource => Self::TemperatureSource,
            DraftInterfaceKind::PistonFace { .. } => Self::PistonFace,
        }
    }
}

impl fmt::Display for InterfaceKindTag {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_str())
    }
}

/// Canonical dispatch key for interface constructors:
/// `(interface kind, card family)`.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct InterfaceLoweringKey {
    /// Interface class tag.
    pub kind: InterfaceKindTag,
    /// Card family name declared by the resolved card entry.
    pub card_family: String,
}

impl InterfaceLoweringKey {
    /// Construct an interface lowering key.
    #[must_use]
    pub fn new(kind: InterfaceKindTag, card_family: impl Into<String>) -> Self {
        Self {
            kind,
            card_family: card_family.into(),
        }
    }
}

impl fmt::Display for InterfaceLoweringKey {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "({}, family=\"{}\")", self.kind, self.card_family)
    }
}

/// Canonical dispatch key for kinematic joint constructors.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct JointLoweringKey(pub JointKind);

impl JointLoweringKey {
    /// Stable canonical label.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self.0 {
            JointKind::Fixed => "fixed",
            JointKind::Revolute => "revolute",
            JointKind::Prismatic => "prismatic",
            JointKind::Helical => "helical",
        }
    }

    const fn discriminant(self) -> u8 {
        match self.0 {
            JointKind::Fixed => 0,
            JointKind::Revolute => 1,
            JointKind::Prismatic => 2,
            JointKind::Helical => 3,
        }
    }
}

impl PartialOrd for JointLoweringKey {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for JointLoweringKey {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.discriminant().cmp(&other.discriminant())
    }
}

impl fmt::Display for JointLoweringKey {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "joint({})", self.as_str())
    }
}

/// Canonical dispatch key for gear-mesh constructors, keyed by mesh card family.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct GearMeshLoweringKey {
    /// Card family name of the resolved gear-mesh card.
    pub card_family: String,
}

impl GearMeshLoweringKey {
    /// Construct a gear-mesh lowering key.
    #[must_use]
    pub fn new(card_family: impl Into<String>) -> Self {
        Self {
            card_family: card_family.into(),
        }
    }
}

impl fmt::Display for GearMeshLoweringKey {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "gear_mesh(family=\"{}\")", self.card_family)
    }
}

/// Canonical dispatch key for acoustic assembly family entries (re-expressing
/// `AcousticAssembly` realization family by family).
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct AcousticFamilyKey(pub String);

impl AcousticFamilyKey {
    /// Construct an acoustic family key (e.g. `"string"`, `"plate"`, `"duct"`,
    /// `"cavity"`, `"reed"`, `"bow"`, `"obstacle"`).
    #[must_use]
    pub fn new(family: impl Into<String>) -> Self {
        Self(family.into())
    }

    /// Borrow the family label.
    #[must_use]
    pub fn as_str(&self) -> &str {
        &self.0
    }
}

impl fmt::Display for AcousticFamilyKey {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "acoustic_family(\"{}\")", self.0)
    }
}

// ============================================================================
// Card Store & Provenance
// ============================================================================

/// Schema card for a [`DraftInterfaceKind::PowerPort`] interface, binding a
/// canonical family tag, power domain, and provenance to a content address.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PowerPortSchemaCard {
    /// Card family used to dispatch into [`LoweringRegistry`].
    pub family: String,
    /// Human-readable schema name (e.g. `"mechanical-1d-bond"`).
    pub schema_name: String,
    /// Power domain carried across the bond.
    pub domain: PortKind,
    /// Description of sign convention and pairing.
    pub description: String,
    /// Provenance record (source + license required).
    pub provenance: Provenance,
}

impl PowerPortSchemaCard {
    /// Validate required non-empty fields and provenance.
    ///
    /// # Errors
    /// Returns [`LoweringError::Invalid`] when required metadata is blank.
    pub fn validate(&self) -> Result<(), LoweringError> {
        if self.family.trim().is_empty() {
            return Err(LoweringError::Invalid {
                item: "power_port_schema".to_owned(),
                detail: "card family must be non-empty".to_owned(),
            });
        }
        if self.schema_name.trim().is_empty() {
            return Err(LoweringError::Invalid {
                item: format!("power_port_schema:{}", self.family),
                detail: "schema_name must be non-empty".to_owned(),
            });
        }
        if self.provenance.source.trim().is_empty() || self.provenance.license.trim().is_empty() {
            return Err(LoweringError::Invalid {
                item: format!("power_port_schema:{}", self.schema_name),
                detail: "provenance source and license must be non-empty".to_owned(),
            });
        }
        Ok(())
    }

    /// Deterministic BLAKE3 content identity of this power-port schema card.
    #[must_use]
    pub fn content_hash(&self) -> ContentHash {
        let mut payload = Vec::new();
        push_len_prefixed(&mut payload, self.family.as_bytes());
        push_len_prefixed(&mut payload, self.schema_name.as_bytes());
        push_len_prefixed(&mut payload, format!("{:?}", self.domain).as_bytes());
        push_len_prefixed(&mut payload, self.description.as_bytes());
        push_len_prefixed(&mut payload, self.provenance.source.as_bytes());
        push_len_prefixed(&mut payload, self.provenance.license.as_bytes());
        if let Some(artifact) = self.provenance.artifact {
            push_len_prefixed(&mut payload, &artifact.0);
        }
        hash_domain(POWER_PORT_SCHEMA_HASH_DOMAIN, &payload)
    }
}

/// Content-addressed interface card binding an [`InterfaceSystemCard`] to its
/// lowering family tag (e.g. `"thermal-conductance"`, `"viscous-damper"`,
/// `"involute-gear-mesh"`).
#[derive(Debug, Clone, PartialEq)]
pub struct InterfaceCardEntry {
    /// Interface law family name used by [`InterfaceLoweringKey`].
    pub family: String,
    /// Immutable `fs-matdb` interface system card holding claims and models.
    pub card: InterfaceSystemCard,
}

impl InterfaceCardEntry {
    /// Canonical content hash of the underlying [`InterfaceSystemCard`].
    #[must_use]
    pub fn content_hash(&self) -> ContentHash {
        self.card.content_hash()
    }
}

/// Immutable-by-hash card catalog consulted during lowering to resolve
/// [`MaterialCard`], [`InterfaceCardEntry`], and [`PowerPortSchemaCard`]
/// references into typed parameters and [`PropertyUsageReceipt`]s.
#[derive(Debug, Clone, Default)]
pub struct CardStore {
    materials: BTreeMap<ContentHash, MaterialCard>,
    interfaces: BTreeMap<ContentHash, InterfaceCardEntry>,
    power_port_schemas: BTreeMap<ContentHash, PowerPortSchemaCard>,
}

impl CardStore {
    /// Create an empty card store.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Insert a validated [`MaterialCard`] indexed by its canonical
    /// `content_hash()`, returning that hash.
    pub fn insert_material(&mut self, card: MaterialCard) -> ContentHash {
        let hash = card.content_hash();
        self.materials.insert(hash, card);
        hash
    }

    /// Insert an [`InterfaceSystemCard`] tagged with its dispatch `family`,
    /// indexed by `card.content_hash()`.
    ///
    /// # Errors
    /// Refuses a blank family name.
    pub fn insert_interface(
        &mut self,
        family: impl Into<String>,
        card: InterfaceSystemCard,
    ) -> Result<ContentHash, LoweringError> {
        let family = family.into();
        if family.trim().is_empty() {
            return Err(LoweringError::Invalid {
                item: "card_store.interface".to_owned(),
                detail: "interface card family must be non-empty".to_owned(),
            });
        }
        let hash = card.content_hash();
        self.interfaces
            .insert(hash, InterfaceCardEntry { family, card });
        Ok(hash)
    }

    /// Insert a [`PowerPortSchemaCard`] indexed by its content hash.
    ///
    /// # Errors
    /// Propagates [`PowerPortSchemaCard::validate`].
    pub fn insert_power_port_schema(
        &mut self,
        card: PowerPortSchemaCard,
    ) -> Result<ContentHash, LoweringError> {
        card.validate()?;
        let hash = card.content_hash();
        self.power_port_schemas.insert(hash, card);
        Ok(hash)
    }

    /// Look up a [`MaterialCard`] by content hash.
    #[must_use]
    pub fn material(&self, hash: &ContentHash) -> Option<&MaterialCard> {
        self.materials.get(hash)
    }

    /// Look up an [`InterfaceCardEntry`] by content hash.
    #[must_use]
    pub fn interface(&self, hash: &ContentHash) -> Option<&InterfaceCardEntry> {
        self.interfaces.get(hash)
    }

    /// Look up a [`PowerPortSchemaCard`] by content hash.
    #[must_use]
    pub fn power_port_schema(&self, hash: &ContentHash) -> Option<&PowerPortSchemaCard> {
        self.power_port_schemas.get(hash)
    }
}

/// Auditable origin of a single numeric parameter consumed while lowering a
/// component, interface, joint, source, or acoustic member.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ParameterSource {
    /// Value evaluated from an `fs-matdb` claim with a persisted
    /// [`PropertyUsageReceipt`].
    MatDbClaim {
        /// Content hash of the material or interface card queried.
        card_hash: ContentHash,
        /// Property name queried.
        property: String,
        /// Canonical content hash of the matching [`PropertyUsageReceipt`].
        receipt_hash: ContentHash,
        /// Selected claim content hash.
        selected_claim: ContentHash,
        /// Whether the selected claim is backed by observations.
        observation_backed: bool,
    },
    /// Value read from a validated `ConstitutiveModelCard` parameter block.
    ConstitutiveModelParam {
        /// Content hash of the owning material or interface card.
        card_hash: ContentHash,
        /// Content hash of the constitutive model card.
        model_hash: ContentHash,
        /// Law identifier.
        law: String,
        /// Law version.
        law_version: u32,
        /// Parameter name in the model's parameter block.
        parameter: String,
    },
    /// Value derived deterministically from resolved geometry (and optionally
    /// from already-receipted material properties named in `formula`).
    GeometryDerived {
        /// Deterministic geometry fingerprint of the owning part.
        geometry_fingerprint: ContentHash,
        /// Generator identifier or `"imported"`.
        generator: String,
        /// Explicit algebraic derivation (e.g. `"mass = density * box_volume"`).
        formula: String,
    },
    /// Value taken directly from an explicitly declared field in the admitted
    /// assembly or environment (e.g. `"joint:slide.stiffness"` or
    /// `"environment.ambient_temperature"`).
    DeclaredAssemblyField {
        /// Item scope (e.g. `"joint:slide"`).
        item: String,
        /// Field name on that item.
        field: String,
    },
}

/// One parameter recorded in a [`LoweringReceipt`], pairing value and SI
/// dimensions with its mandatory [`ParameterSource`].
#[derive(Debug, Clone, PartialEq)]
pub struct ParameterProvenance {
    /// Parameter name used by the constructed component/connection.
    pub parameter: String,
    /// Numeric value in SI base units (must be finite).
    pub value: f64,
    /// Physical dimensions of `value`.
    pub dims: Dims,
    /// Non-optional provenance source.
    pub source: ParameterSource,
}

/// Per-part record inside a [`LoweringReceipt`].
#[derive(Debug, Clone, PartialEq)]
pub struct ComponentLoweringReceipt {
    /// Part name in the assembly draft.
    pub part_name: String,
    /// Persistent catalog entity identity.
    pub entity: EntityId,
    /// Physical part kind.
    pub part_kind: DraftPartKind,
    /// Active physics domain set.
    pub physics: BTreeSet<PhysicsDomain>,
    /// Resolved geometry fingerprint.
    pub geometry_fingerprint: ContentHash,
    /// Fidelity rung actually chosen by the registry.
    pub chosen_fidelity: FidelityRung,
    /// Human- and machine-readable rationale for why `chosen_fidelity` was selected.
    pub fidelity_rationale: String,
    /// Names of kernel nodes created for this part.
    pub created_nodes: Vec<String>,
    /// All `fs-matdb` [`PropertyUsageReceipt`]s consumed when lowering this part.
    pub property_receipts: Vec<PropertyUsageReceipt>,
    /// Every numeric parameter used to construct this part, each with a source.
    pub parameter_sources: Vec<ParameterProvenance>,
}

/// Per-interface record inside a [`LoweringReceipt`].
#[derive(Debug, Clone, PartialEq)]
pub struct InterfaceLoweringReceipt {
    /// Interface name in the assembly draft.
    pub interface_name: String,
    /// Catalog interface entity identity.
    pub entity: EntityId,
    /// Interface class tag.
    pub kind: InterfaceKindTag,
    /// Resolved law/schema card hash (`None` for a card-less heat route).
    pub card_hash: Option<ContentHash>,
    /// Card family dispatched to the constructor.
    pub card_family: String,
    /// Optional third-body material card hash.
    pub third_body_card: Option<ContentHash>,
    /// Names of kernel nodes created for this interface.
    pub created_nodes: Vec<String>,
    /// All `fs-matdb` [`PropertyUsageReceipt`]s consumed.
    pub property_receipts: Vec<PropertyUsageReceipt>,
    /// Every numeric parameter used, each with a source.
    pub parameter_sources: Vec<ParameterProvenance>,
}

/// Per-joint record inside a [`LoweringReceipt`].
#[derive(Debug, Clone, PartialEq)]
pub struct JointLoweringReceipt {
    /// Joint name in the assembly draft.
    pub joint_name: String,
    /// Kinematic joint class.
    pub kind: JointKind,
    /// Parent entity (`root` for ground).
    pub parent: EntityId,
    /// Child entity.
    pub child: EntityId,
    /// Names of kernel nodes created for this joint.
    pub created_nodes: Vec<String>,
    /// Every numeric parameter used, each with a source.
    pub parameter_sources: Vec<ParameterProvenance>,
}

/// Per-gear-mesh record inside a [`LoweringReceipt`].
#[derive(Debug, Clone, PartialEq)]
pub struct GearMeshLoweringReceipt {
    /// Gear-mesh name in the assembly draft.
    pub mesh_name: String,
    /// Driving gear entity.
    pub driver: EntityId,
    /// Driven gear entity.
    pub driven: EntityId,
    /// Mesh card content hash.
    pub mesh_card: ContentHash,
    /// Mesh card family dispatched.
    pub card_family: String,
    /// Kinematic transmission ratio applied.
    pub ratio: f64,
    /// Names of kernel nodes created.
    pub created_nodes: Vec<String>,
    /// All `fs-matdb` [`PropertyUsageReceipt`]s consumed.
    pub property_receipts: Vec<PropertyUsageReceipt>,
    /// Every numeric parameter used, each with a source.
    pub parameter_sources: Vec<ParameterProvenance>,
}

/// Per-source record inside a [`LoweringReceipt`].
#[derive(Debug, Clone, PartialEq)]
pub struct SourceLoweringReceipt {
    /// Source name in the assembly draft.
    pub source_name: String,
    /// Target port or parameter description.
    pub target_description: String,
    /// Names of kernel nodes created (empty for direct boundary bindings).
    pub created_nodes: Vec<String>,
    /// Every numeric parameter used, each with a source.
    pub parameter_sources: Vec<ParameterProvenance>,
}

/// Per-acoustic-member record inside a [`LoweringReceipt`].
#[derive(Debug, Clone, PartialEq)]
pub struct AcousticMemberLoweringReceipt {
    /// Acoustic member name.
    pub member_name: String,
    /// Dispatched acoustic family key.
    pub family: String,
    /// Names of kernel nodes created.
    pub created_nodes: Vec<String>,
    /// All `fs-matdb` [`PropertyUsageReceipt`]s consumed.
    pub property_receipts: Vec<PropertyUsageReceipt>,
    /// Every numeric parameter used, each with a source.
    pub parameter_sources: Vec<ParameterProvenance>,
}

/// Complete provenance receipt emitted for one lowered assembly, suitable for
/// canonical JSON serialization and persistence in `fs-ledger`.
#[derive(Debug, Clone, PartialEq)]
pub struct LoweringReceipt {
    /// Stable schema identifier ([`LOWERING_RECEIPT_SCHEMA`]).
    pub schema: &'static str,
    /// Assembly name.
    pub assembly_name: String,
    /// Deterministic content hash of the input assembly description.
    pub assembly_hash: ContentHash,
    /// Canonical scenario IR hash bound to the assembly.
    pub scenario_ir_hash: ContentHash,
    /// Static execution plan content hash produced after sealing.
    pub plan_hash: ContentHash,
    /// Deterministic BLAKE3 digest of this receipt's canonical JSON payload.
    pub receipt_hash: ContentHash,
    /// Component lowering records in declaration order.
    pub components: Vec<ComponentLoweringReceipt>,
    /// Interface lowering records in declaration order.
    pub interfaces: Vec<InterfaceLoweringReceipt>,
    /// Joint lowering records in declaration order.
    pub joints: Vec<JointLoweringReceipt>,
    /// Gear-mesh lowering records in declaration order.
    pub gear_meshes: Vec<GearMeshLoweringReceipt>,
    /// External source lowering records in declaration order.
    pub sources: Vec<SourceLoweringReceipt>,
    /// Acoustic family member records in declaration order.
    pub acoustic_members: Vec<AcousticMemberLoweringReceipt>,
}

impl LoweringReceipt {
    /// Total count of `fs-matdb` [`PropertyUsageReceipt`]s across all lowered
    /// parts, interfaces, gear meshes, and acoustic members.
    #[must_use]
    pub fn total_property_receipts(&self) -> usize {
        self.components
            .iter()
            .map(|c| c.property_receipts.len())
            .sum::<usize>()
            + self
                .interfaces
                .iter()
                .map(|i| i.property_receipts.len())
                .sum::<usize>()
            + self
                .gear_meshes
                .iter()
                .map(|g| g.property_receipts.len())
                .sum::<usize>()
            + self
                .acoustic_members
                .iter()
                .map(|a| a.property_receipts.len())
                .sum::<usize>()
    }

    /// Audit receipt completeness: every parameter must have a valid non-empty
    /// source and every `MatDbClaim` must match a [`PropertyUsageReceipt`]
    /// stored on the same item.
    ///
    /// # Errors
    /// [`LoweringError::MissingParameterSource`] on any missing, non-finite, or
    /// unlinked parameter provenance entry.
    pub fn verify_completeness(&self) -> Result<(), LoweringError> {
        for c in &self.components {
            let scope = format!("part:{}", c.part_name);
            if c.geometry_fingerprint.0 == [0_u8; 32] {
                return Err(LoweringError::MissingParameterSource {
                    item: scope,
                    parameter: "geometry_fingerprint".to_owned(),
                    detail: "geometry fingerprint must be nonzero".to_owned(),
                });
            }
            if c.fidelity_rationale.trim().is_empty() {
                return Err(LoweringError::MissingParameterSource {
                    item: scope,
                    parameter: "fidelity_rationale".to_owned(),
                    detail: "fidelity selection rationale must be non-empty".to_owned(),
                });
            }
            verify_param_list(&scope, &c.parameter_sources, &c.property_receipts)?;
        }
        for i in &self.interfaces {
            let scope = format!("interface:{}", i.interface_name);
            verify_param_list(&scope, &i.parameter_sources, &i.property_receipts)?;
        }
        for j in &self.joints {
            let scope = format!("joint:{}", j.joint_name);
            verify_param_list(&scope, &j.parameter_sources, &[])?;
        }
        for g in &self.gear_meshes {
            let scope = format!("gear_mesh:{}", g.mesh_name);
            verify_param_list(&scope, &g.parameter_sources, &g.property_receipts)?;
        }
        for s in &self.sources {
            let scope = format!("source:{}", s.source_name);
            verify_param_list(&scope, &s.parameter_sources, &[])?;
        }
        for a in &self.acoustic_members {
            let scope = format!("acoustic:{}", a.member_name);
            verify_param_list(&scope, &a.parameter_sources, &a.property_receipts)?;
        }
        Ok(())
    }

    /// Render the canonical deterministic JSON document suitable for
    /// `fs_ledger::Ledger::record_lowering_receipt`.
    #[must_use]
    pub fn canonical_json(&self) -> String {
        let mut out = String::with_capacity(1024);
        out.push_str("{\"schema\":");
        push_json_str(&mut out, self.schema);
        out.push_str(",\"receipt_hash\":");
        push_json_str(&mut out, &self.receipt_hash.to_hex());
        out.push_str(",\"body\":");
        out.push_str(&self.canonical_body_json());
        out.push('}');
        out
    }

    fn canonical_body_json(&self) -> String {
        let mut out = String::with_capacity(1024);
        out.push_str("{\"assembly_name\":");
        push_json_str(&mut out, &self.assembly_name);
        out.push_str(",\"assembly_hash\":");
        push_json_str(&mut out, &self.assembly_hash.to_hex());
        out.push_str(",\"scenario_ir_hash\":");
        push_json_str(&mut out, &self.scenario_ir_hash.to_hex());
        out.push_str(",\"plan_hash\":");
        push_json_str(&mut out, &self.plan_hash.to_hex());

        out.push_str(",\"components\":");
        push_json_components(&mut out, &self.components);
        out.push_str(",\"interfaces\":");
        push_json_interfaces(&mut out, &self.interfaces);
        out.push_str(",\"joints\":");
        push_json_joints(&mut out, &self.joints);
        out.push_str(",\"gear_meshes\":");
        push_json_gear_meshes(&mut out, &self.gear_meshes);
        out.push_str(",\"sources\":");
        push_json_sources(&mut out, &self.sources);
        out.push_str(",\"acoustic_members\":");
        push_json_acoustic_members(&mut out, &self.acoustic_members);
        out.push('}');
        out
    }

    fn finalize_hash(&mut self) {
        let body = self.canonical_body_json();
        self.receipt_hash = hash_domain(LOWERING_RECEIPT_HASH_DOMAIN, body.as_bytes());
    }
}

fn push_json_components(out: &mut String, components: &[ComponentLoweringReceipt]) {
    out.push('[');
    for (idx, c) in components.iter().enumerate() {
        if idx > 0 {
            out.push(',');
        }
        out.push_str("{\"part_name\":");
        push_json_str(out, &c.part_name);
        out.push_str(",\"entity\":");
        push_json_str(out, &c.entity.token());
        out.push_str(",\"part_kind\":");
        push_json_str(out, c.part_kind.as_str());
        out.push_str(",\"physics\":[");
        for (p_idx, p) in c.physics.iter().enumerate() {
            if p_idx > 0 {
                out.push(',');
            }
            push_json_str(out, p.as_str());
        }
        out.push_str("],\"geometry_fingerprint\":");
        push_json_str(out, &c.geometry_fingerprint.to_hex());
        out.push_str(",\"chosen_fidelity\":");
        push_json_fidelity(out, &c.chosen_fidelity);
        out.push_str(",\"fidelity_rationale\":");
        push_json_str(out, &c.fidelity_rationale);
        out.push_str(",\"created_nodes\":");
        push_json_str_slice(out, &c.created_nodes);
        out.push_str(",\"property_receipts\":");
        push_json_property_receipts(out, &c.property_receipts);
        out.push_str(",\"parameter_sources\":");
        push_json_param_sources(out, &c.parameter_sources);
        out.push('}');
    }
    out.push(']');
}

fn push_json_interfaces(out: &mut String, interfaces: &[InterfaceLoweringReceipt]) {
    out.push('[');
    for (idx, iface) in interfaces.iter().enumerate() {
        if idx > 0 {
            out.push(',');
        }
        out.push_str("{\"interface_name\":");
        push_json_str(out, &iface.interface_name);
        out.push_str(",\"entity\":");
        push_json_str(out, &iface.entity.token());
        out.push_str(",\"kind\":");
        push_json_str(out, iface.kind.as_str());
        out.push_str(",\"card_hash\":");
        match iface.card_hash {
            Some(h) => push_json_str(out, &h.to_hex()),
            None => out.push_str("null"),
        }
        out.push_str(",\"card_family\":");
        push_json_str(out, &iface.card_family);
        out.push_str(",\"third_body_card\":");
        match iface.third_body_card {
            Some(h) => push_json_str(out, &h.to_hex()),
            None => out.push_str("null"),
        }
        out.push_str(",\"created_nodes\":");
        push_json_str_slice(out, &iface.created_nodes);
        out.push_str(",\"property_receipts\":");
        push_json_property_receipts(out, &iface.property_receipts);
        out.push_str(",\"parameter_sources\":");
        push_json_param_sources(out, &iface.parameter_sources);
        out.push('}');
    }
    out.push(']');
}

fn push_json_joints(out: &mut String, joints: &[JointLoweringReceipt]) {
    out.push('[');
    for (idx, j) in joints.iter().enumerate() {
        if idx > 0 {
            out.push(',');
        }
        out.push_str("{\"joint_name\":");
        push_json_str(out, &j.joint_name);
        out.push_str(",\"kind\":");
        push_json_str(out, JointLoweringKey(j.kind).as_str());
        out.push_str(",\"parent\":");
        push_json_str(out, &j.parent.token());
        out.push_str(",\"child\":");
        push_json_str(out, &j.child.token());
        out.push_str(",\"created_nodes\":");
        push_json_str_slice(out, &j.created_nodes);
        out.push_str(",\"parameter_sources\":");
        push_json_param_sources(out, &j.parameter_sources);
        out.push('}');
    }
    out.push(']');
}

fn push_json_gear_meshes(out: &mut String, gear_meshes: &[GearMeshLoweringReceipt]) {
    out.push('[');
    for (idx, g) in gear_meshes.iter().enumerate() {
        if idx > 0 {
            out.push(',');
        }
        out.push_str("{\"mesh_name\":");
        push_json_str(out, &g.mesh_name);
        out.push_str(",\"driver\":");
        push_json_str(out, &g.driver.token());
        out.push_str(",\"driven\":");
        push_json_str(out, &g.driven.token());
        out.push_str(",\"mesh_card\":");
        push_json_str(out, &g.mesh_card.to_hex());
        out.push_str(",\"card_family\":");
        push_json_str(out, &g.card_family);
        out.push_str(",\"ratio\":");
        push_json_f64(out, g.ratio);
        out.push_str(",\"created_nodes\":");
        push_json_str_slice(out, &g.created_nodes);
        out.push_str(",\"property_receipts\":");
        push_json_property_receipts(out, &g.property_receipts);
        out.push_str(",\"parameter_sources\":");
        push_json_param_sources(out, &g.parameter_sources);
        out.push('}');
    }
    out.push(']');
}

fn push_json_sources(out: &mut String, sources: &[SourceLoweringReceipt]) {
    out.push('[');
    for (idx, s) in sources.iter().enumerate() {
        if idx > 0 {
            out.push(',');
        }
        out.push_str("{\"source_name\":");
        push_json_str(out, &s.source_name);
        out.push_str(",\"target\":");
        push_json_str(out, &s.target_description);
        out.push_str(",\"created_nodes\":");
        push_json_str_slice(out, &s.created_nodes);
        out.push_str(",\"parameter_sources\":");
        push_json_param_sources(out, &s.parameter_sources);
        out.push('}');
    }
    out.push(']');
}

fn push_json_acoustic_members(out: &mut String, members: &[AcousticMemberLoweringReceipt]) {
    out.push('[');
    for (idx, a) in members.iter().enumerate() {
        if idx > 0 {
            out.push(',');
        }
        out.push_str("{\"member_name\":");
        push_json_str(out, &a.member_name);
        out.push_str(",\"family\":");
        push_json_str(out, &a.family);
        out.push_str(",\"created_nodes\":");
        push_json_str_slice(out, &a.created_nodes);
        out.push_str(",\"property_receipts\":");
        push_json_property_receipts(out, &a.property_receipts);
        out.push_str(",\"parameter_sources\":");
        push_json_param_sources(out, &a.parameter_sources);
        out.push('}');
    }
    out.push(']');
}

fn verify_param_list(
    scope: &str,
    params: &[ParameterProvenance],
    receipts: &[PropertyUsageReceipt],
) -> Result<(), LoweringError> {
    if params.is_empty() {
        return Err(LoweringError::MissingParameterSource {
            item: scope.to_owned(),
            parameter: "<all>".to_owned(),
            detail: "lowered item recorded zero parameter provenance entries".to_owned(),
        });
    }
    for p in params {
        if p.parameter.trim().is_empty() {
            return Err(LoweringError::MissingParameterSource {
                item: scope.to_owned(),
                parameter: "<empty>".to_owned(),
                detail: "parameter name must be non-empty".to_owned(),
            });
        }
        if !p.value.is_finite() {
            return Err(LoweringError::MissingParameterSource {
                item: scope.to_owned(),
                parameter: p.parameter.clone(),
                detail: format!("parameter value {} is not finite", p.value),
            });
        }
        verify_param_source(scope, p, receipts)?;
    }
    Ok(())
}

fn verify_param_source(
    scope: &str,
    p: &ParameterProvenance,
    receipts: &[PropertyUsageReceipt],
) -> Result<(), LoweringError> {
    match &p.source {
        ParameterSource::MatDbClaim {
            card_hash,
            property,
            receipt_hash,
            selected_claim,
            ..
        } => {
            if card_hash.0 == [0_u8; 32]
                || receipt_hash.0 == [0_u8; 32]
                || selected_claim.0 == [0_u8; 32]
                || property.trim().is_empty()
            {
                return Err(LoweringError::MissingParameterSource {
                    item: scope.to_owned(),
                    parameter: p.parameter.clone(),
                    detail: "MatDbClaim source requires nonzero hashes and non-empty property"
                        .to_owned(),
                });
            }
            let matched = receipts.iter().any(|r| {
                r.content_hash() == *receipt_hash
                    && r.property == *property
                    && r.selected.0 == *selected_claim
            });
            if !matched {
                return Err(LoweringError::MissingParameterSource {
                    item: scope.to_owned(),
                    parameter: p.parameter.clone(),
                    detail: format!(
                        "MatDbClaim cites receipt {} for property '{}' which is absent from item property_receipts",
                        receipt_hash.to_hex(),
                        property
                    ),
                });
            }
        }
        ParameterSource::ConstitutiveModelParam {
            card_hash,
            model_hash,
            law,
            parameter,
            ..
        } => {
            if card_hash.0 == [0_u8; 32]
                || model_hash.0 == [0_u8; 32]
                || law.trim().is_empty()
                || parameter.trim().is_empty()
            {
                return Err(LoweringError::MissingParameterSource {
                    item: scope.to_owned(),
                    parameter: p.parameter.clone(),
                    detail: "ConstitutiveModelParam requires nonzero hashes, law, and parameter"
                        .to_owned(),
                });
            }
        }
        ParameterSource::GeometryDerived {
            geometry_fingerprint,
            generator,
            formula,
        } => {
            if geometry_fingerprint.0 == [0_u8; 32]
                || generator.trim().is_empty()
                || formula.trim().is_empty()
            {
                return Err(LoweringError::MissingParameterSource {
                    item: scope.to_owned(),
                    parameter: p.parameter.clone(),
                    detail: "GeometryDerived requires nonzero fingerprint, generator, and formula"
                        .to_owned(),
                });
            }
        }
        ParameterSource::DeclaredAssemblyField { item, field } => {
            if item.trim().is_empty() || field.trim().is_empty() {
                return Err(LoweringError::MissingParameterSource {
                    item: scope.to_owned(),
                    parameter: p.parameter.clone(),
                    detail: "DeclaredAssemblyField requires non-empty item and field".to_owned(),
                });
            }
        }
    }
    Ok(())
}

// ============================================================================
// Errors (Do 3: Named Refusals Stating Missing Constructor)
// ============================================================================

/// Structured, teaching refusal from [`LoweringRegistry`] lowering.
#[derive(Debug)]
pub enum LoweringError {
    /// No registered part constructor matches any fidelity rung in the part's
    /// ladder.
    MissingPartConstructor {
        /// Offending part name.
        part: String,
        /// Part kind.
        kind: DraftPartKind,
        /// Active physics domains.
        physics: BTreeSet<PhysicsDomain>,
        /// Rungs attempted in preference order.
        tried_rungs: Vec<FidelityRungKind>,
    },
    /// No registered interface constructor matches `(kind, card_family)`.
    MissingInterfaceConstructor {
        /// Offending interface name.
        interface: String,
        /// Interface class tag.
        kind: InterfaceKindTag,
        /// Card family resolved from [`CardStore`].
        card_family: String,
    },
    /// No registered joint constructor matches `kind`.
    MissingJointConstructor {
        /// Offending joint name.
        joint: String,
        /// Missing joint class.
        kind: JointKind,
    },
    /// No registered gear-mesh constructor matches `card_family`.
    MissingGearMeshConstructor {
        /// Offending gear-mesh name.
        mesh: String,
        /// Card family resolved from [`CardStore`].
        card_family: String,
    },
    /// No registered acoustic family constructor matches `family`.
    MissingAcousticFamilyConstructor {
        /// Offending acoustic member name.
        member: String,
        /// Missing acoustic family key.
        family: String,
    },
    /// A content-addressed material, interface, or schema card is missing from
    /// [`CardStore`].
    MissingCard {
        /// Offending assembly item name.
        item: String,
        /// Card table expected (`"material"`, `"interface"`, `"power_port_schema"`).
        card_kind: &'static str,
        /// Missing content hash.
        card_hash: ContentHash,
    },
    /// A constructor attempted to query or default a property marked
    /// `UnknownProperty` on a `MaterialRef`.
    UnknownPropertyDefaulted {
        /// Offending item name.
        item: String,
        /// Property marked unknown.
        property: String,
    },
    /// Receipt completeness gate refused an un-sourced or invalid parameter.
    MissingParameterSource {
        /// Offending item scope.
        item: String,
        /// Offending parameter name.
        parameter: String,
        /// Teaching reason.
        detail: String,
    },
    /// A joint, interface, or source referenced a port not exposed by the
    /// lowered target entity.
    MissingPort {
        /// Offending item scope.
        item: String,
        /// Target entity.
        entity: EntityId,
        /// Missing port name.
        port: String,
    },
    /// Failure querying an `fs-matdb` card (extrapolation, missing property,
    /// ambiguity, etc.).
    MatDb {
        /// Offending item scope.
        item: String,
        /// Underlying `fs-matdb` refusal.
        source: Box<MatDbError>,
    },
    /// Failure constructing a kernel element.
    Element {
        /// Offending item scope.
        item: String,
        /// Underlying element refusal.
        source: ElementError,
    },
    /// Failure during kernel graph admission (`SystemBuilder::seal` or
    /// connection wiring).
    Admission(Box<AdmissionError>),
    /// General structural or lowering validation refusal.
    Invalid {
        /// Offending item scope.
        item: String,
        /// Teaching detail.
        detail: String,
    },
}

impl fmt::Display for LoweringError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MissingPartConstructor {
                part,
                kind,
                physics,
                tried_rungs,
            } => {
                write!(
                    f,
                    "lowering refusal for part '{part}': missing part constructor for ({}, {{",
                    kind.as_str()
                )?;
                for (i, d) in physics.iter().enumerate() {
                    if i > 0 {
                        f.write_str(", ")?;
                    }
                    f.write_str(d.as_str())?;
                }
                f.write_str("}) across fidelity rungs [")?;
                for (i, r) in tried_rungs.iter().enumerate() {
                    if i > 0 {
                        f.write_str(", ")?;
                    }
                    write!(f, "{r}")?;
                }
                f.write_str("]")
            }
            Self::MissingInterfaceConstructor {
                interface,
                kind,
                card_family,
            } => write!(
                f,
                "lowering refusal for interface '{interface}': missing interface constructor for ({kind}, family=\"{card_family}\")"
            ),
            Self::MissingJointConstructor { joint, kind } => write!(
                f,
                "lowering refusal for joint '{joint}': missing joint constructor for {}",
                JointLoweringKey(*kind)
            ),
            Self::MissingGearMeshConstructor { mesh, card_family } => write!(
                f,
                "lowering refusal for gear_mesh '{mesh}': missing gear-mesh constructor for family=\"{card_family}\""
            ),
            Self::MissingAcousticFamilyConstructor { member, family } => write!(
                f,
                "lowering refusal for acoustic member '{member}': missing acoustic family constructor for \"{family}\""
            ),
            Self::MissingCard {
                item,
                card_kind,
                card_hash,
            } => write!(
                f,
                "lowering refusal for {item}: missing {card_kind} card {}",
                card_hash.to_hex()
            ),
            Self::UnknownPropertyDefaulted { item, property } => write!(
                f,
                "lowering refusal for {item}: property '{property}' is explicitly declared unknown and cannot be defaulted"
            ),
            Self::MissingParameterSource {
                item,
                parameter,
                detail,
            } => write!(
                f,
                "lowering receipt refusal for {item} parameter '{parameter}': {detail}"
            ),
            Self::MissingPort { item, entity, port } => write!(
                f,
                "lowering refusal for {item}: entity {} has no exposed port '{port}'",
                entity.token()
            ),
            Self::MatDb { item, source } => {
                write!(f, "lowering matdb refusal for {item}: {source}")
            }
            Self::Element { item, source } => {
                write!(f, "lowering element refusal for {item}: {source}")
            }
            Self::Admission(err) => write!(f, "lowering graph admission refusal: {err}"),
            Self::Invalid { item, detail } => {
                write!(f, "lowering refusal for {item}: {detail}")
            }
        }
    }
}

impl std::error::Error for LoweringError {}

impl From<AdmissionError> for LoweringError {
    fn from(value: AdmissionError) -> Self {
        Self::Admission(Box::new(value))
    }
}

// ============================================================================
// Geometry Fingerprinting & Parametric Mensuration Helpers
// ============================================================================

/// Compute the deterministic BLAKE3 fingerprint of a [`DraftGeometrySource`].
///
/// Imported charts preserve their exact catalog fingerprint; parametric
/// generators hash their canonical generator name and ordered dimensioned
/// parameter table under [`LOWERING_GEOMETRY_HASH_DOMAIN`].
#[must_use]
pub fn resolve_geometry_fingerprint(geom: &DraftGeometrySource) -> ContentHash {
    match geom {
        DraftGeometrySource::Standard(GeometrySource::Imported { fingerprint, .. }) => {
            fingerprint.hash()
        }
        DraftGeometrySource::Standard(src @ GeometrySource::Parametric { .. }) => {
            hash_parametric_geometry(src)
        }
        DraftGeometrySource::Gear(g) => hash_parametric_geometry(&g.to_geometry_source()),
        DraftGeometrySource::SpiralSpring(s) => hash_parametric_geometry(&s.to_geometry_source()),
        DraftGeometrySource::Primitive(p) => hash_parametric_geometry(&p.to_geometry_source()),
    }
}

/// Human-readable generator label for receipts.
#[must_use]
pub fn geometry_generator_label(geom: &DraftGeometrySource) -> String {
    match geom {
        DraftGeometrySource::Standard(GeometrySource::Parametric { generator, .. }) => {
            generator.clone()
        }
        DraftGeometrySource::Standard(GeometrySource::Imported {
            representation,
            ledger_row,
            ..
        }) => format!("imported:{representation:?}:{ledger_row}"),
        DraftGeometrySource::Gear(_) => "Gear.draft.v1".to_owned(),
        DraftGeometrySource::SpiralSpring(_) => "SpiralSpring.draft.v1".to_owned(),
        DraftGeometrySource::Primitive(p) => match p {
            PrimitiveShape::Box { .. } => "Primitive.Box.draft.v1".to_owned(),
            PrimitiveShape::Cylinder { .. } => "Primitive.Cylinder.draft.v1".to_owned(),
            PrimitiveShape::Sphere { .. } => "Primitive.Sphere.draft.v1".to_owned(),
            PrimitiveShape::Tube { .. } => "Primitive.Tube.draft.v1".to_owned(),
        },
    }
}

/// Exact analytical volume [m^3] for closed-form parametric primitives and gears.
#[must_use]
pub fn geometry_volume_m3(geom: &DraftGeometrySource) -> Option<(f64, &'static str)> {
    match geom {
        DraftGeometrySource::Primitive(PrimitiveShape::Box { lx, ly, lz }) => {
            Some((lx.value * ly.value * lz.value, "volume = lx * ly * lz"))
        }
        DraftGeometrySource::Primitive(PrimitiveShape::Cylinder { radius, height }) => Some((
            PI * radius.value * radius.value * height.value,
            "volume = pi * radius^2 * height",
        )),
        DraftGeometrySource::Primitive(PrimitiveShape::Sphere { radius }) => Some((
            (4.0 / 3.0) * PI * radius.value.powi(3),
            "volume = (4/3) * pi * radius^3",
        )),
        DraftGeometrySource::Primitive(PrimitiveShape::Tube {
            inner_radius,
            outer_radius,
            length,
        }) => Some((
            PI * (outer_radius.value * outer_radius.value
                - inner_radius.value * inner_radius.value)
                * length.value,
            "volume = pi * (r_outer^2 - r_inner^2) * length",
        )),
        DraftGeometrySource::Gear(g) => {
            let r_pitch = g.pitch_radius_m();
            let r_bore = g.bore_radius.value;
            Some((
                PI * (r_pitch * r_pitch - r_bore * r_bore) * g.face_width.value,
                "volume = pi * (r_pitch^2 - r_bore^2) * face_width",
            ))
        }
        DraftGeometrySource::SpiralSpring(s) => {
            let mean_r = 0.5 * (s.inner_radius.value + s.outer_radius.value);
            let active_length = 2.0 * PI * mean_r * s.turns;
            Some((
                active_length * s.blade_thickness.value * s.blade_height.value,
                "volume = (2 * pi * r_mean * turns) * blade_thickness * blade_height",
            ))
        }
        DraftGeometrySource::Standard(_) => None,
    }
}

/// Analytical axial polar moment of inertia `I_zz / rho` [m^5] for rotational
/// primitives (`Gear`, `Cylinder`, `Tube`) so rotational inertia `J = rho * I_zz_over_rho`.
#[must_use]
pub fn geometry_polar_inertia_over_rho(geom: &DraftGeometrySource) -> Option<(f64, &'static str)> {
    match geom {
        DraftGeometrySource::Gear(g) => {
            let r_pitch = g.pitch_radius_m();
            let r_bore = g.bore_radius.value;
            let vol = PI * (r_pitch * r_pitch - r_bore * r_bore) * g.face_width.value;
            Some((
                0.5 * vol * (r_pitch * r_pitch + r_bore * r_bore),
                "J = 0.5 * rho * V * (r_pitch^2 + r_bore^2)",
            ))
        }
        DraftGeometrySource::Primitive(PrimitiveShape::Cylinder { radius, height }) => {
            let vol = PI * radius.value * radius.value * height.value;
            Some((
                0.5 * vol * radius.value * radius.value,
                "J = 0.5 * rho * V * radius^2",
            ))
        }
        DraftGeometrySource::Primitive(PrimitiveShape::Tube {
            inner_radius,
            outer_radius,
            length,
        }) => {
            let r_o = outer_radius.value;
            let r_i = inner_radius.value;
            let vol = PI * (r_o * r_o - r_i * r_i) * length.value;
            Some((
                0.5 * vol * (r_o * r_o + r_i * r_i),
                "J = 0.5 * rho * V * (r_outer^2 + r_inner^2)",
            ))
        }
        _ => None,
    }
}

/// Compute axial bar/rod geometric factor `A / L` [m] for a `SpringBeam` part.
#[must_use]
pub fn geometry_axial_area_over_length_m(
    geom: &DraftGeometrySource,
) -> Option<(f64, &'static str)> {
    match geom {
        DraftGeometrySource::Primitive(PrimitiveShape::Box { lx, ly, lz }) => {
            Some(((ly.value * lz.value) / lx.value, "k = E * (ly * lz) / lx"))
        }
        DraftGeometrySource::Primitive(PrimitiveShape::Cylinder { radius, height }) => Some((
            (PI * radius.value * radius.value) / height.value,
            "k = E * (pi * radius^2) / height",
        )),
        DraftGeometrySource::Primitive(PrimitiveShape::Tube {
            inner_radius,
            outer_radius,
            length,
        }) => Some((
            (PI * (outer_radius.value.powi(2) - inner_radius.value.powi(2))) / length.value,
            "k = E * pi * (r_outer^2 - r_inner^2) / length",
        )),
        _ => None,
    }
}

/// Compute torsional spring geometric factor `b * h^3 / (12 * L_coil)` [m^3]
/// for a [`SpiralSpringGenerator`], so torsional stiffness is `k_theta = E * factor`.
#[must_use]
pub fn spiral_spring_stiffness_factor_m3(s: &SpiralSpringGenerator) -> f64 {
    let mean_r = 0.5 * (s.inner_radius.value + s.outer_radius.value);
    let active_length = 2.0 * PI * mean_r * s.turns;
    (s.blade_height.value * s.blade_thickness.value.powi(3)) / (12.0 * active_length)
}

fn hash_parametric_geometry(src: &GeometrySource) -> ContentHash {
    let mut payload = Vec::new();
    match src {
        GeometrySource::Parametric { generator, params } => {
            push_len_prefixed(&mut payload, b"parametric");
            push_len_prefixed(&mut payload, generator.as_bytes());
            payload.extend_from_slice(&(params.len() as u64).to_le_bytes());
            for (key, qty) in params {
                push_len_prefixed(&mut payload, key.name.as_bytes());
                for d in key.dims.0 {
                    payload.push(d as u8);
                }
                payload.extend_from_slice(&qty.value.to_bits().to_le_bytes());
            }
        }
        GeometrySource::Imported { fingerprint, .. } => {
            return fingerprint.hash();
        }
    }
    hash_domain(LOWERING_GEOMETRY_HASH_DOMAIN, &payload)
}

// ============================================================================
// Graph Assembly & Deferred Junctions
// ============================================================================

struct HeatRouteSpec {
    source: PortRef,
    sinks: Vec<(PortRef, f64)>,
    remainder: Option<HeatSink>,
    temperature: RouteTemperature,
}

struct SealedLoweringGraph {
    system: System,
    node_map: BTreeMap<String, NodeId>,
    entity_nodes: BTreeMap<EntityId, Vec<NodeId>>,
    entity_ports: BTreeMap<(EntityId, String), PortRef>,
    observers: BTreeMap<EntityId, Observer>,
}

/// Mutable graph construction surface passed to registered constructors.
///
/// Wraps [`SystemBuilder`] while tracking exposed entity ports and aggregating
/// multi-member junctions (`common_flow`, `common_effort`, `heat_route`) so
/// independent joints, interfaces, and sources targeting the same part compose
/// automatically without knowing how many peers attach to that port.
pub struct LoweringGraph {
    builder: SystemBuilder,
    node_map: BTreeMap<String, NodeId>,
    entity_nodes: BTreeMap<EntityId, Vec<NodeId>>,
    entity_ports: BTreeMap<(EntityId, String), PortRef>,
    port_meta: BTreeMap<PortRef, (Carrier, PortRole, usize)>,
    common_flow_hubs: BTreeMap<PortRef, Vec<(PortRef, f64)>>,
    common_effort_hubs: BTreeMap<PortRef, Vec<(PortRef, f64)>>,
    heat_routes: Vec<HeatRouteSpec>,
    boundary_zero_ports: BTreeSet<PortRef>,
    ground_entities: BTreeSet<EntityId>,
    initial_guesses: BTreeMap<PortRef, Vec<f64>>,
    observers: BTreeMap<EntityId, Observer>,
    consumed_sources: BTreeMap<usize, SourceLoweringReceipt>,
}

impl LoweringGraph {
    fn new() -> Self {
        Self {
            builder: SystemBuilder::new(),
            node_map: BTreeMap::new(),
            entity_nodes: BTreeMap::new(),
            entity_ports: BTreeMap::new(),
            port_meta: BTreeMap::new(),
            common_flow_hubs: BTreeMap::new(),
            common_effort_hubs: BTreeMap::new(),
            heat_routes: Vec::new(),
            boundary_zero_ports: BTreeSet::new(),
            ground_entities: BTreeSet::new(),
            initial_guesses: BTreeMap::new(),
            observers: BTreeMap::new(),
            consumed_sources: BTreeMap::new(),
        }
    }

    /// Borrow the underlying [`SystemBuilder`] for direct component addition or
    /// custom connections.
    pub fn builder_mut(&mut self) -> &mut SystemBuilder {
        &mut self.builder
    }

    /// Borrow the underlying [`SystemBuilder`] immutably.
    #[must_use]
    pub fn builder(&self) -> &SystemBuilder {
        &self.builder
    }

    /// Add a [`fs_ports::component::Component`] with [`fs_ports::component::PortTangent`],
    /// recording its node name in the lowering index.
    ///
    /// # Errors
    /// Propagates [`AdmissionError`].
    pub fn add_with_tangent<C>(
        &mut self,
        name: &str,
        owner: Option<EntityId>,
        component: C,
        state: C::State,
    ) -> Result<NodeId, LoweringError>
    where
        C: fs_ports::component::PortTangent + 'static,
        C::State: 'static,
    {
        let specs: Vec<(Carrier, PortRole, usize)> = component
            .ports()
            .iter()
            .map(|p| (p.carrier(), p.role(), p.width()))
            .collect();
        let node = self.builder.add_with_tangent(name, component, state)?;
        self.register_node_meta(name, owner, node, &specs);
        Ok(node)
    }

    /// Add a [`fs_ports::component::Component`] with both
    /// [`fs_ports::component::PortTangent`] and
    /// [`fs_ports::component::StageResidual`], recording its node name in the
    /// lowering index.
    ///
    /// # Errors
    /// Propagates [`AdmissionError`].
    pub fn add_with_tangent_and_residual<C>(
        &mut self,
        name: &str,
        owner: Option<EntityId>,
        component: C,
        state: C::State,
    ) -> Result<NodeId, LoweringError>
    where
        C: fs_ports::component::PortTangent + fs_ports::component::StageResidual + 'static,
        C::State: 'static,
    {
        let specs: Vec<(Carrier, PortRole, usize)> = component
            .ports()
            .iter()
            .map(|p| (p.carrier(), p.role(), p.width()))
            .collect();
        let node = self
            .builder
            .add_with_tangent_and_stage_residual(name, component, state)?;
        self.register_node_meta(name, owner, node, &specs);
        Ok(node)
    }

    /// Add a plain [`fs_ports::component::Component`], recording its node name
    /// in the lowering index.
    ///
    /// # Errors
    /// Propagates [`AdmissionError`].
    pub fn add<C>(
        &mut self,
        name: &str,
        owner: Option<EntityId>,
        component: C,
        state: C::State,
    ) -> Result<NodeId, LoweringError>
    where
        C: fs_ports::component::Component + 'static,
        C::State: 'static,
    {
        let specs: Vec<(Carrier, PortRole, usize)> = component
            .ports()
            .iter()
            .map(|p| (p.carrier(), p.role(), p.width()))
            .collect();
        let node = self.builder.add(name, component, state)?;
        self.register_node_meta(name, owner, node, &specs);
        Ok(node)
    }

    /// Add a non-smooth [`fs_ports::component::LoadLine`] component to the
    /// lowering graph, recording its node name and port metadata.
    ///
    /// # Errors
    /// Propagates [`AdmissionError`].
    pub fn add_on_line<C>(
        &mut self,
        name: &str,
        owner: Option<EntityId>,
        component: C,
        state: C::State,
    ) -> Result<NodeId, LoweringError>
    where
        C: fs_ports::component::LoadLine + 'static,
        C::State: 'static,
    {
        let specs: Vec<(Carrier, PortRole, usize)> = component
            .ports()
            .iter()
            .map(|p| (p.carrier(), p.role(), p.width()))
            .collect();
        let node = self.builder.add_on_line(name, component, state)?;
        self.register_node_meta(name, owner, node, &specs);
        Ok(node)
    }

    /// Connect a signal output port to a signal input port of equal width.
    ///
    /// # Errors
    /// Propagates [`AdmissionError`].
    pub fn connect_signal(&mut self, from: PortRef, to: PortRef) -> Result<(), LoweringError> {
        Ok(self.builder.signal(from, to)?)
    }

    fn register_node_meta(
        &mut self,
        name: &str,
        owner: Option<EntityId>,
        node: NodeId,
        specs: &[(Carrier, PortRole, usize)],
    ) {
        self.node_map.insert(name.to_owned(), node);
        if let Some(entity) = owner {
            self.entity_nodes.entry(entity).or_default().push(node);
        }
        for (port, &(carrier, role, width)) in specs.iter().enumerate() {
            self.port_meta
                .insert(PortRef { node: node.0, port }, (carrier, role, width));
        }
    }

    /// Expose a named port on an assembly entity so joints, interfaces, gear
    /// meshes, and sources can look it up by `(entity, logical_port_name)`.
    pub fn expose_port(&mut self, entity: EntityId, port_name: impl Into<String>, port: PortRef) {
        self.entity_ports.insert((entity, port_name.into()), port);
    }

    /// Look up an exposed port on `entity`.
    #[must_use]
    pub fn entity_port(&self, entity: EntityId, port_name: &str) -> Option<PortRef> {
        self.entity_ports
            .get(&(entity, port_name.to_owned()))
            .copied()
    }

    /// Look up an exposed port on `entity`, walking up its parent chain in
    /// `catalog` when `entity` is a child surface/region feature of a part.
    #[must_use]
    pub fn resolve_entity_port(
        &self,
        catalog: &EntityCatalog,
        entity: EntityId,
        port_name: &str,
    ) -> Option<PortRef> {
        let mut cur = Some(entity);
        for _ in 0..=catalog.entities().len() {
            let id = cur?;
            if let Some(port) = self.entity_port(id, port_name) {
                return Some(port);
            }
            cur = catalog.get(id).and_then(Entity::parent);
        }
        None
    }

    /// Mark `entity` (for example a stationary `DrivenBody` frame/chassis part)
    /// as a zero-velocity reference ground so surfaces owned by it terminate
    /// mechanical/electrical contact endpoints at zero boundary input.
    pub fn mark_ground_entity(&mut self, entity: EntityId) {
        self.ground_entities.insert(entity);
    }

    /// Return whether `entity` represents the ground assembly root, a marked
    /// ground entity, or a surface/region child whose ancestor is grounded.
    #[must_use]
    pub fn is_ground_feature(
        &self,
        catalog: &EntityCatalog,
        entity: EntityId,
        root: EntityId,
    ) -> bool {
        if entity == root || self.ground_entities.contains(&entity) {
            return true;
        }
        let mut cur = catalog.get(entity).and_then(Entity::parent);
        if cur == Some(root) {
            return entity.kind() == EntityKind::Surface;
        }
        for _ in 0..=catalog.entities().len() {
            let Some(id) = cur else { break };
            if self.ground_entities.contains(&id) {
                return true;
            }
            if id == root {
                return false;
            }
            cur = catalog.get(id).and_then(Entity::parent);
        }
        false
    }

    /// Register an explicit starting guess for a looped port input (applied
    /// right after [`SystemBuilder::seal`]).
    pub fn set_initial_guess(&mut self, port: PortRef, values: Vec<f64>) {
        self.initial_guesses.insert(port, values);
    }

    /// Look up the `(Carrier, PortRole)` recorded for a port added through this
    /// graph.
    #[must_use]
    pub fn port_meta(&self, port: PortRef) -> Option<(Carrier, PortRole)> {
        self.port_meta
            .get(&port)
            .map(|&(carrier, role, _)| (carrier, role))
    }

    /// Attach a `FlowIn` power port to a shared `common_flow` junction driven by
    /// `setter` (an `EffortIn` port such as a rigid body's mechanical port),
    /// scaled by `ratio` (`1.0` for unit-ratio direct connection).
    pub fn attach_common_flow(&mut self, setter: PortRef, member: PortRef, ratio: f64) {
        self.common_flow_hubs
            .entry(setter)
            .or_default()
            .push((member, ratio));
    }

    /// Attach an `EffortIn` port to a shared `common_effort` junction driven by
    /// `setter` (a `FlowIn` port such as a `HeatCapacity`'s thermal port or an
    /// electrical/acoustic `Storage` port), scaled by `ratio`.
    pub fn attach_common_effort(&mut self, setter: PortRef, member: PortRef, ratio: f64) {
        self.common_effort_hubs
            .entry(setter)
            .or_default()
            .push((member, ratio));
    }

    /// Schedule a `heat_route` from a dissipative heat source port to one or
    /// more thermal sinks and/or an external remainder [`HeatSink`].
    pub fn add_heat_route(
        &mut self,
        source: PortRef,
        sinks: Vec<(PortRef, f64)>,
        remainder: Option<HeatSink>,
        temperature: RouteTemperature,
    ) {
        self.heat_routes.push(HeatRouteSpec {
            source,
            sinks,
            remainder,
            temperature,
        });
    }

    /// Bind a port to prescribed zero input (e.g. ground velocity = 0 m/s for a
    /// spring or damper anchored to the assembly root).
    pub fn bind_zero_boundary(&mut self, port: PortRef) {
        self.boundary_zero_ports.insert(port);
    }

    /// Look up a `PortRef` on a previously added node.
    ///
    /// # Errors
    /// Propagates [`AdmissionError::UnknownPort`].
    pub fn port(&self, node: NodeId, id: &str) -> Result<PortRef, LoweringError> {
        Ok(self.builder.port(node, id)?)
    }

    fn finalize_and_seal(
        mut self,
        ambient_temp_k: f64,
    ) -> Result<SealedLoweringGraph, LoweringError> {
        for (setter, members) in self.common_flow_hubs {
            let unit_ratios = members.iter().all(|&(_, r)| (r - 1.0).abs() == 0.0);
            if unit_ratios {
                let ports: Vec<PortRef> = members.iter().map(|&(p, _)| p).collect();
                self.builder.common_flow(setter, &ports)?;
            } else {
                self.builder.common_flow_ratios(setter, &members)?;
            }
        }
        for (setter, members) in self.common_effort_hubs {
            let unit_ratios = members.iter().all(|&(_, r)| (r - 1.0).abs() == 0.0);
            if unit_ratios {
                let ports: Vec<PortRef> = members.iter().map(|&(p, _)| p).collect();
                self.builder.common_effort(setter, &ports)?;
            } else {
                self.builder.common_effort_ratios(setter, &members)?;
            }
        }
        for route in self.heat_routes {
            self.builder.heat_route(
                route.source,
                &route.sinks,
                route.remainder,
                route.temperature,
            )?;
        }
        for port in self.boundary_zero_ports {
            self.builder.boundary(
                port,
                Box::new(|_t_s: f64, _dt_s: f64, out: &mut [f64]| out.fill(0.0)),
            )?;
        }
        let mut system = self.builder.seal()?;
        // Seed positive ambient temperature starting guesses for looped heat
        // effort-in ports so thermal loops start inside the valid T > 0 K
        // domain without altering any conservation law.
        if ambient_temp_k.is_finite() && ambient_temp_k > 0.0 {
            for (&port, &(carrier, role, _)) in &self.port_meta {
                if carrier == Carrier::Heat && role == PortRole::EffortIn {
                    let _ = system.guess(port, &[ambient_temp_k]);
                }
            }
        }
        for (port, values) in self.initial_guesses {
            let _ = system.guess(port, &values);
        }
        Ok(SealedLoweringGraph {
            system,
            node_map: self.node_map,
            entity_nodes: self.entity_nodes,
            entity_ports: self.entity_ports,
            observers: self.observers,
        })
    }
}

// ============================================================================
// Constructor Contexts & Output Records
// ============================================================================

/// Ambient query state and claim-selection policy passed to `fs-matdb` lookups
/// through [`ProvenanceRecorder`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MatDbQueryEnv {
    /// Ambient temperature [K] supplied when the item's query point omits `T`.
    pub ambient_temp_k: f64,
    /// Ambient pressure [Pa] supplied when the item's query point omits `P`.
    pub ambient_pressure_pa: f64,
    /// Material claim selection policy.
    pub policy: SelectionPolicy,
}

/// Identifying reference to a parameter sourced from a `ConstitutiveModelCard`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ConstitutiveModelRef<'a> {
    /// Content hash of the enclosing material or interface card.
    pub card_hash: ContentHash,
    /// Content hash of the constitutive model card.
    pub model_hash: ContentHash,
    /// Constitutive law name.
    pub law: &'a str,
    /// Constitutive law schema version.
    pub law_version: u32,
    /// Parameter name inside the constitutive model.
    pub parameter: &'a str,
}

/// Provenance-tracking query helper shared across constructor contexts so every
/// `fs-matdb` lookup automatically records its [`PropertyUsageReceipt`] and
/// [`ParameterProvenance`].
#[derive(Debug, Default)]
pub struct ProvenanceRecorder {
    property_receipts: Vec<PropertyUsageReceipt>,
    parameter_sources: Vec<ParameterProvenance>,
}

impl ProvenanceRecorder {
    /// Query a property from a [`MaterialCard`] referenced by `mat_ref`,
    /// checking [`MaterialRef::require_known`] first and recording both the
    /// [`PropertyUsageReceipt`] and a [`ParameterProvenance`] entry.
    ///
    /// # Errors
    /// Refuses missing cards, explicitly unknown properties, out-of-domain
    /// queries, missing properties, or ambiguous selections.
    pub fn query_material_property(
        &mut self,
        item_scope: &str,
        cards: &CardStore,
        mat_ref: &MaterialRef,
        property: &str,
        parameter_name: &str,
        env: MatDbQueryEnv,
    ) -> Result<(PropertySample, PropertyUsageReceipt), LoweringError> {
        if mat_ref.require_known(property).is_err() {
            return Err(LoweringError::UnknownPropertyDefaulted {
                item: item_scope.to_owned(),
                property: property.to_owned(),
            });
        }
        let card = cards
            .material(&mat_ref.card)
            .ok_or(LoweringError::MissingCard {
                item: item_scope.to_owned(),
                card_kind: "material",
                card_hash: mat_ref.card,
            })?;
        let point = build_query_point(
            item_scope,
            &mat_ref.query_point,
            env.ambient_temp_k,
            env.ambient_pressure_pa,
            (card.claims(), property),
        )?;
        let k_spec =
            QuantitySpec::semantic(SemanticType::new(ThermalConductivity, ValueForm::Static));
        let k_list = card.claims().claims_for(property);
        let is_k = matches!(property, "thermal-conductivity" | "thermal_conductivity")
            && !k_list.is_empty()
            && k_list.iter().all(|(_, c)| c.key.quantity() == k_spec);
        let answer = if is_k {
            let key = PropertyKey::with_quantity(property, k_spec);
            card.claims().query_typed(&key, &point, env.policy)
        } else {
            card.claims().query(property, &point, env.policy)
        }
        .map_err(|source| LoweringError::MatDb {
            item: item_scope.to_owned(),
            source: Box::new(source),
        })?;
        let sample = answer.evidence.value.clone();
        let receipt = answer.receipt;
        self.record_matdb_claim(parameter_name, &sample, mat_ref.card, &receipt);
        Ok((sample, receipt))
    }

    /// Query a property from an [`InterfaceSystemCard`], automatically
    /// recording both its [`PropertyUsageReceipt`] and [`ParameterProvenance`].
    ///
    /// # Errors
    /// Propagates [`LoweringError::MatDb`] on query refusals.
    pub fn query_interface_property(
        &mut self,
        item_scope: &str,
        entry: &InterfaceCardEntry,
        query_state: &[(fs_scenario::assembly::ParameterKey, fs_qty::QtyAny)],
        property: &str,
        parameter_name: &str,
        env: MatDbQueryEnv,
    ) -> Result<(PropertySample, PropertyUsageReceipt), LoweringError> {
        let point = build_query_point(
            item_scope,
            query_state,
            env.ambient_temp_k,
            env.ambient_pressure_pa,
            (entry.card.claims(), property),
        )?;
        let answer = entry
            .card
            .claims()
            .query(property, &point, env.policy)
            .map_err(|source| LoweringError::MatDb {
                item: item_scope.to_owned(),
                source: Box::new(source),
            })?;
        let sample = answer.evidence.value.clone();
        let receipt = answer.receipt;
        self.record_matdb_claim(parameter_name, &sample, entry.content_hash(), &receipt);
        Ok((sample, receipt))
    }

    /// Query a property from a [`ThirdBodyDecl`]'s material card.
    ///
    /// # Errors
    /// Refuses missing cards or matdb query failures.
    pub fn query_third_body_property(
        &mut self,
        item_scope: &str,
        cards: &CardStore,
        third_body: &ThirdBodyDecl,
        property: &str,
        parameter_name: &str,
        env: MatDbQueryEnv,
    ) -> Result<(PropertySample, PropertyUsageReceipt), LoweringError> {
        let card = cards
            .material(&third_body.material_card)
            .ok_or(LoweringError::MissingCard {
                item: item_scope.to_owned(),
                card_kind: "material",
                card_hash: third_body.material_card,
            })?;
        let point = build_query_point(
            item_scope,
            &third_body.query_point,
            env.ambient_temp_k,
            env.ambient_pressure_pa,
            (card.claims(), property),
        )?;
        let answer = card
            .claims()
            .query(property, &point, env.policy)
            .map_err(|source| LoweringError::MatDb {
                item: item_scope.to_owned(),
                source: Box::new(source),
            })?;
        let sample = answer.evidence.value.clone();
        let receipt = answer.receipt;
        self.record_matdb_claim(parameter_name, &sample, third_body.material_card, &receipt);
        Ok((sample, receipt))
    }

    /// Record a parameter backed by an already-produced [`PropertyUsageReceipt`].
    pub fn record_matdb_claim(
        &mut self,
        parameter: impl Into<String>,
        sample: &PropertySample,
        card_hash: ContentHash,
        receipt: &PropertyUsageReceipt,
    ) {
        let receipt_hash = receipt.content_hash();
        self.property_receipts.push(receipt.clone());
        self.parameter_sources.push(ParameterProvenance {
            parameter: parameter.into(),
            value: sample.value,
            dims: sample.dims,
            source: ParameterSource::MatDbClaim {
                card_hash,
                property: receipt.property.clone(),
                receipt_hash,
                selected_claim: receipt.selected.0,
                observation_backed: receipt.observation_backed,
            },
        });
    }

    /// Record a parameter derived from geometry (and/or receipted properties).
    pub fn record_geometry_derived(
        &mut self,
        parameter: impl Into<String>,
        value: f64,
        dims: Dims,
        geometry_fingerprint: ContentHash,
        generator: impl Into<String>,
        formula: impl Into<String>,
    ) {
        self.parameter_sources.push(ParameterProvenance {
            parameter: parameter.into(),
            value,
            dims,
            source: ParameterSource::GeometryDerived {
                geometry_fingerprint,
                generator: generator.into(),
                formula: formula.into(),
            },
        });
    }

    /// Record a parameter taken directly from an explicitly declared assembly
    /// or environment field.
    pub fn record_declared_field(
        &mut self,
        parameter: impl Into<String>,
        value: f64,
        dims: Dims,
        item: impl Into<String>,
        field: impl Into<String>,
    ) {
        self.parameter_sources.push(ParameterProvenance {
            parameter: parameter.into(),
            value,
            dims,
            source: ParameterSource::DeclaredAssemblyField {
                item: item.into(),
                field: field.into(),
            },
        });
    }

    /// Record a parameter from a `ConstitutiveModelCard`.
    pub fn record_model_parameter(
        &mut self,
        parameter: impl Into<String>,
        value: f64,
        dims: Dims,
        model: ConstitutiveModelRef<'_>,
    ) {
        self.parameter_sources.push(ParameterProvenance {
            parameter: parameter.into(),
            value,
            dims,
            source: ParameterSource::ConstitutiveModelParam {
                card_hash: model.card_hash,
                model_hash: model.model_hash,
                law: model.law.to_owned(),
                law_version: model.law_version,
                parameter: model.parameter.to_owned(),
            },
        });
    }

    /// Consume recorded property receipts and parameter provenances.
    #[must_use]
    pub fn finish(self) -> (Vec<PropertyUsageReceipt>, Vec<ParameterProvenance>) {
        (self.property_receipts, self.parameter_sources)
    }
}

fn build_query_point(
    item_scope: &str,
    declared: &[(fs_scenario::assembly::ParameterKey, fs_qty::QtyAny)],
    ambient_temp_k: f64,
    ambient_pressure_pa: f64,
    (claims, property): (&ClaimSet, &str),
) -> Result<QueryPoint, LoweringError> {
    let abs_t = QuantitySpec::semantic(SemanticType::new(AbsTemp, ValueForm::Static));
    let list = claims.claims_for(property);
    let typed_t = !list.is_empty()
        && list
            .iter()
            .all(|(_, c)| c.validity.axis_quantities().get("T") == Some(&abs_t));
    // A declared key takes the one spec its bounding claims agree on; never a guess. A bare
    // quantity has dims alone, so the one semantic kind it takes is absolute temperature
    // (kelvin); under another it stays untyped and the query refuses the axis.
    let spec = |name: &str, dims: Dims| {
        let spec_of = |c: &fs_matdb::PropertyClaim| c.validity.axis_quantities().get(name).copied();
        let mut bounding = list
            .iter()
            .filter(|(_, c)| c.validity.bound(name).is_some());
        let first = bounding.next().and_then(|(_, c)| spec_of(c));
        let agreed = bounding.all(|(_, c)| spec_of(c) == first);
        match first {
            _ if !agreed => Err("its claims disagree on its quantity"),
            Some(s) if s.dims() != dims => Err("its dims are not its claims' axis dims"),
            Some(s) if s.semantic_type().is_none() || s == abs_t => Ok(Some(s)),
            _ => Ok(None),
        }
        .map_err(|why| LoweringError::Invalid {
            item: item_scope.to_owned(),
            detail: format!("query key `{name}`: {why}"),
        })
    };
    let mut point = QueryPoint::new();
    let mut has_temp = false;
    let mut has_pressure = false;
    for (k, q) in declared {
        if k.name == "T" {
            has_temp = true;
        }
        if k.name == "P" {
            has_pressure = true;
        }
        point = match spec(&k.name, q.dims)? {
            Some(s) => point.with_quantity(&k.name, s, q.value),
            None => point.with(&k.name, q.value),
        }
        .map_err(|source| LoweringError::MatDb {
            item: item_scope.to_owned(),
            source: Box::new(source),
        })?;
    }
    if !has_temp && ambient_temp_k.is_finite() {
        point = if typed_t {
            point.with_quantity("T", abs_t, ambient_temp_k)
        } else {
            point.with("T", ambient_temp_k)
        }
        .map_err(|source| LoweringError::MatDb {
            item: item_scope.to_owned(),
            source: Box::new(source),
        })?;
    }
    if !has_pressure && ambient_pressure_pa.is_finite() {
        point = point
            .with("P", ambient_pressure_pa)
            .map_err(|source| LoweringError::MatDb {
                item: item_scope.to_owned(),
                source: Box::new(source),
            })?;
    }
    Ok(point)
}

/// Context passed to a registered [`PartConstructor`].
pub struct PartLoweringContext<'a> {
    /// Assembly being lowered.
    pub draft: &'a DraftAssembly,
    /// Part declaration being lowered.
    pub part: &'a DraftPartDecl,
    /// Chosen fidelity rung from `part.fidelity.rungs`.
    pub chosen_rung: &'a FidelityRung,
    /// Deterministic geometry fingerprint of `part.geometry`.
    pub geometry_fingerprint: ContentHash,
    /// Card store for `fs-matdb` lookups.
    pub cards: &'a CardStore,
    /// Material claim selection policy.
    pub selection_policy: SelectionPolicy,
}

impl PartLoweringContext<'_> {
    /// Ambient temperature [K] from the assembly environment.
    #[must_use]
    pub fn ambient_temperature_k(&self) -> f64 {
        self.draft
            .input()
            .environment
            .base
            .ambient_temperature
            .value
    }

    /// Ambient pressure [Pa] from the assembly environment.
    #[must_use]
    pub fn ambient_pressure_pa(&self) -> f64 {
        self.draft.input().environment.base.ambient_pressure.value
    }

    /// Ambient state and claim selection policy for `fs-matdb` lookups.
    #[must_use]
    pub fn query_env(&self) -> MatDbQueryEnv {
        MatDbQueryEnv {
            ambient_temp_k: self.ambient_temperature_k(),
            ambient_pressure_pa: self.ambient_pressure_pa(),
            policy: self.selection_policy,
        }
    }

    /// Convenience helper to query a property from `part.materials[material_idx]`.
    ///
    /// # Errors
    /// Refuses missing material index, missing card, unknown property, or
    /// `fs-matdb` query errors.
    pub fn query_material(
        &self,
        recorder: &mut ProvenanceRecorder,
        material_idx: usize,
        property: &str,
        parameter_name: &str,
    ) -> Result<(PropertySample, PropertyUsageReceipt), LoweringError> {
        let scope = format!("part:{}", self.part.name);
        let mat_ref =
            self.part
                .materials
                .get(material_idx)
                .ok_or_else(|| LoweringError::Invalid {
                    item: scope.clone(),
                    detail: format!("material reference index {material_idx} is missing"),
                })?;
        recorder.query_material_property(
            &scope,
            self.cards,
            mat_ref,
            property,
            parameter_name,
            self.query_env(),
        )
    }
}

/// Output returned by a [`PartConstructor`].
#[derive(Debug, Default)]
pub struct LoweredPartOutput {
    /// Names of kernel nodes created for this part.
    pub created_nodes: Vec<String>,
    /// All `fs-matdb` property receipts consumed.
    pub property_receipts: Vec<PropertyUsageReceipt>,
    /// All parameter provenances consumed.
    pub parameter_sources: Vec<ParameterProvenance>,
}

/// Resolved card passed to an [`InterfaceConstructor`].
#[derive(Debug, Clone, Copy)]
pub enum ResolvedInterfaceCard<'a> {
    /// `fs-matdb` interface system card entry (contact, thermal, electrical, glue).
    System(&'a InterfaceCardEntry),
    /// Power-port schema card entry.
    PowerPort(&'a PowerPortSchemaCard),
    /// No card: a [`DraftInterfaceKind::HeatRoute`], [`DraftInterfaceKind::ShearFilm`] or
    /// [`DraftInterfaceKind::TemperatureSource`].
    None,
}

/// Context passed to a registered [`InterfaceConstructor`].
pub struct InterfaceLoweringContext<'a> {
    /// Assembly being lowered.
    pub draft: &'a DraftAssembly,
    /// Interface declaration being lowered.
    pub interface: &'a DraftInterfaceDecl,
    /// Resolved card hash (`None` for a card-less heat route).
    pub card_hash: Option<ContentHash>,
    /// Resolved card entry.
    pub card: ResolvedInterfaceCard<'a>,
    /// Card store (for third-body lookups).
    pub cards: &'a CardStore,
    /// Material claim selection policy.
    pub selection_policy: SelectionPolicy,
}

impl InterfaceLoweringContext<'_> {
    /// Ambient temperature [K] from the assembly environment.
    #[must_use]
    pub fn ambient_temperature_k(&self) -> f64 {
        self.draft
            .input()
            .environment
            .base
            .ambient_temperature
            .value
    }

    /// Ambient pressure [Pa] from the assembly environment.
    #[must_use]
    pub fn ambient_pressure_pa(&self) -> f64 {
        self.draft.input().environment.base.ambient_pressure.value
    }

    /// Ambient state and claim selection policy for `fs-matdb` lookups.
    #[must_use]
    pub fn query_env(&self) -> MatDbQueryEnv {
        MatDbQueryEnv {
            ambient_temp_k: self.ambient_temperature_k(),
            ambient_pressure_pa: self.ambient_pressure_pa(),
            policy: self.selection_policy,
        }
    }
}

/// Output returned by an [`InterfaceConstructor`].
#[derive(Debug, Default)]
pub struct LoweredInterfaceOutput {
    /// Names of kernel nodes created for this interface.
    pub created_nodes: Vec<String>,
    /// All `fs-matdb` property receipts consumed.
    pub property_receipts: Vec<PropertyUsageReceipt>,
    /// All parameter provenances consumed.
    pub parameter_sources: Vec<ParameterProvenance>,
}

/// Context passed to a registered [`JointConstructor`].
pub struct JointLoweringContext<'a> {
    /// Assembly being lowered.
    pub draft: &'a DraftAssembly,
    /// Joint declaration being lowered.
    pub joint: &'a JointDecl,
    /// Card store.
    pub cards: &'a CardStore,
}

/// Output returned by a [`JointConstructor`].
#[derive(Debug, Default)]
pub struct LoweredJointOutput {
    /// Names of kernel nodes created for this joint.
    pub created_nodes: Vec<String>,
    /// All parameter provenances consumed.
    pub parameter_sources: Vec<ParameterProvenance>,
}

/// Context passed to a registered [`GearMeshConstructor`].
pub struct GearMeshLoweringContext<'a> {
    /// Assembly being lowered.
    pub draft: &'a DraftAssembly,
    /// Gear-mesh declaration being lowered.
    pub mesh: &'a GearMeshDecl,
    /// Resolved mesh card entry.
    pub card: &'a InterfaceCardEntry,
    /// Card store.
    pub cards: &'a CardStore,
    /// Material claim selection policy.
    pub selection_policy: SelectionPolicy,
}

/// Output returned by a [`GearMeshConstructor`].
#[derive(Debug, Default)]
pub struct LoweredGearMeshOutput {
    /// Kinematic transmission ratio applied.
    pub ratio: f64,
    /// Names of kernel nodes created.
    pub created_nodes: Vec<String>,
    /// All `fs-matdb` property receipts consumed.
    pub property_receipts: Vec<PropertyUsageReceipt>,
    /// All parameter provenances consumed.
    pub parameter_sources: Vec<ParameterProvenance>,
}

// ============================================================================
// AcousticAssembly Family Registry Types (Do 5)
// ============================================================================

/// One acoustic member described by family key + explicitly sourced parameters
/// (re-expressing `AcousticAssembly` members family by family without a
/// monolithic slot table).
#[derive(Debug, Clone, PartialEq)]
pub struct AcousticMemberSpec {
    /// Unique member name (e.g. `"e4_string"`, `"spruce_top"`, `"bore_duct"`).
    pub name: String,
    /// Family dispatch key (e.g. `"string"`, `"plate"`, `"duct"`, `"cavity"`,
    /// `"reed"`, `"bow"`, `"obstacle"`).
    pub family: String,
    /// Optional `fs-matdb` material card reference for sourced properties.
    pub material_card: Option<ContentHash>,
    /// Dimensioned geometric / tuning parameters with SI units.
    pub parameters: BTreeMap<String, (f64, Dims)>,
    /// Target member name to couple mechanically/acoustically to, if any.
    pub coupled_to: Option<String>,
}

/// Acoustic assembly description lowered via registered family entries in
/// [`LoweringRegistry`].
#[derive(Debug, Clone, PartialEq)]
pub struct AcousticAssemblySpec {
    /// Instrument / device name.
    pub name: String,
    /// Sample rate [Hz].
    pub sample_rate_hz: u32,
    /// Ambient temperature [K].
    pub ambient_temperature_k: f64,
    /// Ambient pressure [Pa].
    pub ambient_pressure_pa: f64,
    /// Ordered list of acoustic members.
    pub members: Vec<AcousticMemberSpec>,
}

/// Context passed to a registered [`AcousticFamilyConstructor`].
pub struct AcousticFamilyContext<'a> {
    /// Whole acoustic spec being lowered.
    pub spec: &'a AcousticAssemblySpec,
    /// Member being lowered.
    pub member: &'a AcousticMemberSpec,
    /// Card store.
    pub cards: &'a CardStore,
    /// Material claim selection policy.
    pub selection_policy: SelectionPolicy,
    /// Ports exposed by previously lowered acoustic members (`member_name -> PortRef`).
    pub member_ports: &'a BTreeMap<String, PortRef>,
}

/// Output returned by an [`AcousticFamilyConstructor`].
#[derive(Debug, Default)]
pub struct LoweredAcousticMemberOutput {
    /// Names of kernel nodes created.
    pub created_nodes: Vec<String>,
    /// Optional primary coupling port exposed by this member to subsequent members.
    pub exposed_port: Option<PortRef>,
    /// All `fs-matdb` property receipts consumed.
    pub property_receipts: Vec<PropertyUsageReceipt>,
    /// All parameter provenances consumed.
    pub parameter_sources: Vec<ParameterProvenance>,
}

/// Constructor callback registered for a [`PartLoweringKey`].
pub type PartConstructor = Box<
    dyn Fn(
        &PartLoweringContext<'_>,
        &mut LoweringGraph,
    ) -> Result<LoweredPartOutput, LoweringError>,
>;

/// Constructor callback registered for an [`InterfaceLoweringKey`].
pub type InterfaceConstructor = Box<
    dyn Fn(
        &InterfaceLoweringContext<'_>,
        &mut LoweringGraph,
    ) -> Result<LoweredInterfaceOutput, LoweringError>,
>;

/// Constructor callback registered for a [`JointLoweringKey`].
pub type JointConstructor = Box<
    dyn Fn(
        &JointLoweringContext<'_>,
        &mut LoweringGraph,
    ) -> Result<LoweredJointOutput, LoweringError>,
>;

/// Constructor callback registered for a [`GearMeshLoweringKey`].
pub type GearMeshConstructor = Box<
    dyn Fn(
        &GearMeshLoweringContext<'_>,
        &mut LoweringGraph,
    ) -> Result<LoweredGearMeshOutput, LoweringError>,
>;

/// Constructor callback registered for an [`AcousticFamilyKey`].
pub type AcousticFamilyConstructor = Box<
    dyn Fn(
        &AcousticFamilyContext<'_>,
        &mut LoweringGraph,
    ) -> Result<LoweredAcousticMemberOutput, LoweringError>,
>;

/// Which parts, interfaces and sources of a draft one owner Component takes
/// whole, as indices into the draft's `parts`, `interfaces` and `sources`.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct GroupSelection {
    /// Indices into `DraftAssemblyInput::parts`, in assembly order.
    pub parts: Vec<usize>,
    /// Indices into `DraftAssemblyInput::interfaces`, in assembly order.
    pub interfaces: Vec<usize>,
    /// Indices into `DraftAssemblyInput::sources`, in assembly order.
    pub sources: Vec<usize>,
}

/// Context passed to a registered [`GroupConstructor`].
pub struct GroupLoweringContext<'a> {
    /// Assembly being lowered.
    pub draft: &'a DraftAssembly,
    /// The items this constructor lowers.
    pub selection: &'a GroupSelection,
    /// Card store for `fs-matdb` lookups.
    pub cards: &'a CardStore,
    /// Material claim selection policy.
    pub selection_policy: SelectionPolicy,
}

/// Output of a [`GroupConstructor`]: one record per selected item, in the
/// selection's order.
#[derive(Debug, Default)]
pub struct LoweredGroupOutput {
    /// One record per selected part.
    pub parts: Vec<LoweredPartOutput>,
    /// One record per selected interface.
    pub interfaces: Vec<LoweredInterfaceOutput>,
    /// One receipt per selected source.
    pub sources: Vec<SourceLoweringReceipt>,
}

/// Selector of an owner group: the selections it lowers in one draft.
pub type GroupSelector =
    Box<dyn Fn(&DraftAssembly, &CardStore) -> Result<Vec<GroupSelection>, LoweringError>>;

/// Constructor of an owner group: lowers one [`GroupSelection`].
pub type GroupConstructor = Box<
    dyn Fn(
        &GroupLoweringContext<'_>,
        &mut LoweringGraph,
    ) -> Result<LoweredGroupOutput, LoweringError>,
>;

struct GroupEntry {
    name: String,
    selector: GroupSelector,
    constructor: GroupConstructor,
}

/// Record that `group` lowered item `index`; refuse an item two groups take.
fn claim<T>(
    claimed: &mut BTreeMap<usize, (String, T)>,
    index: usize,
    group: &str,
    item: impl FnOnce() -> String,
    out: T,
) -> Result<(), LoweringError> {
    if let Some((first, _)) = claimed.get(&index) {
        return Err(LoweringError::Invalid {
            item: item(),
            detail: format!("selected by owner groups '{first}' and '{group}'"),
        });
    }
    claimed.insert(index, (group.to_owned(), out));
    Ok(())
}

// ============================================================================
// Lowering Output
// ============================================================================

/// Complete output of lowering an assembly: the sealed kernel [`System`], its
/// static [`ExecutionPlan`], its [`CausalizationAudit`], its complete
/// [`LoweringReceipt`], and lookup tables from assembly entities/names to
/// kernel [`NodeId`]s and [`PortRef`]s.
pub struct LoweredSystem {
    /// Admitted, sealed composition kernel system ready to step or analyze.
    pub system: System,
    /// Static structural execution plan computed after sealing.
    pub plan: ExecutionPlan,
    /// Static causalization audit computed after sealing.
    pub causalization_audit: CausalizationAudit,
    /// Complete provenance receipt ready to store in `fs-ledger`.
    pub receipt: LoweringReceipt,
    /// Mapping from created kernel node name to [`NodeId`].
    pub node_map: BTreeMap<String, NodeId>,
    /// Mapping from catalog [`EntityId`] to all kernel [`NodeId`]s owned by it.
    pub entity_nodes: BTreeMap<EntityId, Vec<NodeId>>,
    /// Mapping from `(EntityId, logical_port_name)` to [`PortRef`].
    pub part_ports: BTreeMap<(EntityId, String), PortRef>,
    /// Readers of each part's slot in an owner node's state (a body in a
    /// network, a cell in a timeline), registered by the owner constructors:
    /// `observers.get(&part).and_then(|read| read(system))`. A part with a node
    /// of its own has none.
    pub observers: BTreeMap<EntityId, Observer>,
}

/// One part's motion or temperature, read from its slot in its owner node.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Observed {
    /// A rigid body's centre-of-mass velocity (m/s) and angular velocity
    /// (rad/s), both in the assembly frame.
    Body {
        /// Centre-of-mass velocity, m/s.
        velocity_m_s: [f64; 3],
        /// Angular velocity, rad/s.
        angular_velocity_rad_s: [f64; 3],
    },
    /// A cell's temperature, K, computed from its stored enthalpy.
    Cell(f64),
    /// A 1-D thermal conductance's steady heat rate [W] entering port `"a"`.
    HeatRate(f64),
}

/// Reads one part's [`Observed`] value from a state of its lowered system.
pub type Observer = Box<dyn Fn(&System) -> Option<Observed>>;

// ============================================================================
// Lowering Registry
// ============================================================================

/// Owner-extensible registry translating canonical assembly drafts and acoustic
/// specs into sealed [`System`] graphs, static [`ExecutionPlan`]s, and
/// auditable [`LoweringReceipt`]s.
pub struct LoweringRegistry {
    parts: BTreeMap<PartLoweringKey, PartConstructor>,
    interfaces: BTreeMap<InterfaceLoweringKey, InterfaceConstructor>,
    joints: BTreeMap<JointLoweringKey, JointConstructor>,
    gear_meshes: BTreeMap<GearMeshLoweringKey, GearMeshConstructor>,
    acoustic_families: BTreeMap<AcousticFamilyKey, AcousticFamilyConstructor>,
    groups: Vec<GroupEntry>,
    fidelity_policy: FidelitySelectionPolicy,
    selection_policy: SelectionPolicy,
}

impl Default for LoweringRegistry {
    fn default() -> Self {
        Self::new()
    }
}

impl LoweringRegistry {
    /// Create an empty registry with no pre-registered constructors.
    #[must_use]
    pub fn new() -> Self {
        Self {
            parts: BTreeMap::new(),
            interfaces: BTreeMap::new(),
            joints: BTreeMap::new(),
            gear_meshes: BTreeMap::new(),
            acoustic_families: BTreeMap::new(),
            groups: Vec::new(),
            fidelity_policy: FidelitySelectionPolicy::FirstRegistered,
            selection_policy: SelectionPolicy::SingleClaimOnly,
        }
    }

    /// Create a registry populated with the kernel's default lumped
    /// multi-domain constructors (mechanical, rotational, thermal, contacts,
    /// joints, gear meshes, and acoustic family adapters).
    ///
    /// Physics owners can extend or override any entry with the `register_*`
    /// methods.
    #[must_use]
    pub fn with_kernel_defaults() -> Self {
        let mut reg = Self::new();
        register_kernel_defaults(&mut reg);
        reg
    }

    /// Configure how a rung is chosen from each part's `FidelityPolicy::rungs`.
    #[must_use]
    pub const fn with_fidelity_policy(mut self, policy: FidelitySelectionPolicy) -> Self {
        self.fidelity_policy = policy;
        self
    }

    /// Configure the `fs-matdb` [`SelectionPolicy`] used when querying material
    /// and interface cards.
    #[must_use]
    pub const fn with_selection_policy(mut self, policy: SelectionPolicy) -> Self {
        self.selection_policy = policy;
        self
    }

    /// Register an owner group: every selection `selector` makes in a draft
    /// lowers whole through `constructor`, before any per-item constructor
    /// runs. Groups lower in registration order; registering a name again
    /// replaces that group in its place.
    pub fn register_group<S, F>(&mut self, name: impl Into<String>, selector: S, constructor: F)
    where
        S: Fn(&DraftAssembly, &CardStore) -> Result<Vec<GroupSelection>, LoweringError> + 'static,
        F: Fn(
                &GroupLoweringContext<'_>,
                &mut LoweringGraph,
            ) -> Result<LoweredGroupOutput, LoweringError>
            + 'static,
    {
        let entry = GroupEntry {
            name: name.into(),
            selector: Box::new(selector),
            constructor: Box::new(constructor),
        };
        match self
            .groups
            .iter_mut()
            .find(|group| group.name == entry.name)
        {
            Some(slot) => *slot = entry,
            None => self.groups.push(entry),
        }
    }

    /// Register a part constructor for `(part kind, physics set, fidelity rung)`.
    pub fn register_part<F>(&mut self, key: PartLoweringKey, constructor: F)
    where
        F: Fn(
                &PartLoweringContext<'_>,
                &mut LoweringGraph,
            ) -> Result<LoweredPartOutput, LoweringError>
            + 'static,
    {
        self.parts.insert(key, Box::new(constructor));
    }

    /// Register an interface constructor for `(interface kind, card family)`.
    pub fn register_interface<F>(&mut self, key: InterfaceLoweringKey, constructor: F)
    where
        F: Fn(
                &InterfaceLoweringContext<'_>,
                &mut LoweringGraph,
            ) -> Result<LoweredInterfaceOutput, LoweringError>
            + 'static,
    {
        self.interfaces.insert(key, Box::new(constructor));
    }

    /// Register a kinematic joint constructor for `kind`.
    pub fn register_joint<F>(&mut self, kind: JointKind, constructor: F)
    where
        F: Fn(
                &JointLoweringContext<'_>,
                &mut LoweringGraph,
            ) -> Result<LoweredJointOutput, LoweringError>
            + 'static,
    {
        self.joints
            .insert(JointLoweringKey(kind), Box::new(constructor));
    }

    /// Register a gear-mesh constructor for `card_family`.
    pub fn register_gear_mesh<F>(&mut self, card_family: impl Into<String>, constructor: F)
    where
        F: Fn(
                &GearMeshLoweringContext<'_>,
                &mut LoweringGraph,
            ) -> Result<LoweredGearMeshOutput, LoweringError>
            + 'static,
    {
        self.gear_meshes
            .insert(GearMeshLoweringKey::new(card_family), Box::new(constructor));
    }

    /// Register an acoustic assembly family constructor (Do 5).
    pub fn register_acoustic_family<F>(&mut self, family: impl Into<String>, constructor: F)
    where
        F: Fn(
                &AcousticFamilyContext<'_>,
                &mut LoweringGraph,
            ) -> Result<LoweredAcousticMemberOutput, LoweringError>
            + 'static,
    {
        self.acoustic_families
            .insert(AcousticFamilyKey::new(family), Box::new(constructor));
    }

    /// List all registered part keys in canonical order.
    #[must_use]
    pub fn registered_part_keys(&self) -> Vec<&PartLoweringKey> {
        self.parts.keys().collect()
    }

    /// List all registered interface keys in canonical order.
    #[must_use]
    pub fn registered_interface_keys(&self) -> Vec<&InterfaceLoweringKey> {
        self.interfaces.keys().collect()
    }

    /// List all registered joint keys in canonical order.
    #[must_use]
    pub fn registered_joint_keys(&self) -> Vec<JointLoweringKey> {
        self.joints.keys().copied().collect()
    }

    /// List all registered acoustic family keys in canonical order.
    #[must_use]
    pub fn registered_acoustic_families(&self) -> Vec<&AcousticFamilyKey> {
        self.acoustic_families.keys().collect()
    }

    /// Lower an admitted [`DraftAssembly`] against `cards`:
    ///
    /// 1. Resolve geometry fingerprints and pick each part's fidelity rung.
    /// 2. Dispatch parts, joints, gear meshes, interfaces, and sources onto one
    ///    [`LoweringGraph`].
    /// 3. Finalize shared junctions and run kernel admission (`seal()`).
    /// 4. Hand the sealed [`System`] to the planner (`execution_plan_for_dt` +
    ///    `causalization_audit`).
    /// 5. Assemble and completeness-audit the [`LoweringReceipt`].
    ///
    /// # Errors
    /// Returns a named [`LoweringError`] stating the offending item and missing
    /// constructor, card, source, or admission failure.
    #[allow(clippy::too_many_lines)]
    pub fn lower(
        &self,
        draft: &DraftAssembly,
        cards: &CardStore,
    ) -> Result<LoweredSystem, LoweringError> {
        let input = draft.input();
        if let Some(controller) = input.controllers.first() {
            return Err(LoweringError::Invalid {
                item: format!("controller:{}", controller.name),
                detail: "controller lowering is not registered; construct and connect a signal-port controller explicitly instead of omitting its feedback law".to_owned(),
            });
        }
        let mut graph = LoweringGraph::new();
        let mut component_receipts = Vec::with_capacity(input.parts.len());
        let mut joint_receipts = Vec::with_capacity(input.joints.len());
        let mut gear_mesh_receipts = Vec::with_capacity(input.gear_meshes.len());
        let mut interface_receipts = Vec::with_capacity(input.interfaces.len());
        let mut source_receipts = Vec::with_capacity(input.sources.len());

        // 0. Owner groups, in registration order: each selection lowers
        //    whole into one owner Component, and its items skip the
        //    per-item constructors below.
        let mut group_parts = BTreeMap::new();
        let mut group_interfaces = BTreeMap::new();
        let mut group_sources = BTreeMap::new();
        let beyond = |items: &[usize], count: usize| items.iter().any(|&i| i >= count);
        for group in &self.groups {
            for selection in (group.selector)(draft, cards)? {
                let scope = format!("group:{}", group.name);
                if beyond(&selection.parts, input.parts.len())
                    || beyond(&selection.interfaces, input.interfaces.len())
                    || beyond(&selection.sources, input.sources.len())
                {
                    return Err(LoweringError::Invalid {
                        item: scope,
                        detail: "selection names an item the draft does not have".to_owned(),
                    });
                }
                let ctx = GroupLoweringContext {
                    draft,
                    selection: &selection,
                    cards,
                    selection_policy: self.selection_policy,
                };
                let out = (group.constructor)(&ctx, &mut graph)?;
                if out.parts.len() != selection.parts.len()
                    || out.interfaces.len() != selection.interfaces.len()
                    || out.sources.len() != selection.sources.len()
                {
                    return Err(LoweringError::Invalid {
                        item: scope,
                        detail: "one record per selected item is required".to_owned(),
                    });
                }
                for (&i, part_out) in selection.parts.iter().zip(out.parts) {
                    let item = || format!("part:{}", input.parts[i].name);
                    claim(&mut group_parts, i, &group.name, item, part_out)?;
                }
                for (&i, iface_out) in selection.interfaces.iter().zip(out.interfaces) {
                    let item = || format!("interface:{}", input.interfaces[i].name);
                    claim(&mut group_interfaces, i, &group.name, item, iface_out)?;
                }
                for (&i, receipt) in selection.sources.iter().zip(out.sources) {
                    let item = || format!("source:{}", input.sources[i].name);
                    claim(&mut group_sources, i, &group.name, item, receipt)?;
                }
            }
        }

        // 1. Lower parts.
        for (index, part) in input.parts.iter().enumerate() {
            if let Some((group, out)) = group_parts.remove(&index) {
                let chosen_fidelity =
                    part.fidelity
                        .rungs
                        .first()
                        .cloned()
                        .ok_or_else(|| LoweringError::Invalid {
                            item: format!("part:{}", part.name),
                            detail: "part declares no fidelity rung".to_owned(),
                        })?;
                component_receipts.push(ComponentLoweringReceipt {
                    part_name: part.name.clone(),
                    entity: part.entity,
                    part_kind: part.kind,
                    physics: part.physics.clone(),
                    geometry_fingerprint: resolve_geometry_fingerprint(&part.geometry),
                    chosen_fidelity,
                    fidelity_rationale: format!("lowered whole by owner group '{group}'"),
                    created_nodes: out.created_nodes,
                    property_receipts: out.property_receipts,
                    parameter_sources: out.parameter_sources,
                });
                continue;
            }
            let (chosen_rung, rationale, ctor) = self.select_part_constructor(part)?;
            let geometry_fingerprint = resolve_geometry_fingerprint(&part.geometry);
            let ctx = PartLoweringContext {
                draft,
                part,
                chosen_rung,
                geometry_fingerprint,
                cards,
                selection_policy: self.selection_policy,
            };
            let out = ctor(&ctx, &mut graph)?;
            component_receipts.push(ComponentLoweringReceipt {
                part_name: part.name.clone(),
                entity: part.entity,
                part_kind: part.kind,
                physics: part.physics.clone(),
                geometry_fingerprint,
                chosen_fidelity: chosen_rung.clone(),
                fidelity_rationale: rationale,
                created_nodes: out.created_nodes,
                property_receipts: out.property_receipts,
                parameter_sources: out.parameter_sources,
            });
        }

        // 2. Lower kinematic joints.
        for joint in &input.joints {
            let key = JointLoweringKey(joint.kind);
            let ctor =
                self.joints
                    .get(&key)
                    .ok_or_else(|| LoweringError::MissingJointConstructor {
                        joint: joint.name.clone(),
                        kind: joint.kind,
                    })?;
            let ctx = JointLoweringContext {
                draft,
                joint,
                cards,
            };
            let out = ctor(&ctx, &mut graph)?;
            joint_receipts.push(JointLoweringReceipt {
                joint_name: joint.name.clone(),
                kind: joint.kind,
                parent: joint.parent,
                child: joint.child,
                created_nodes: out.created_nodes,
                parameter_sources: out.parameter_sources,
            });
        }

        // 3. Lower gear meshes.
        for mesh in &input.gear_meshes {
            let entry =
                cards
                    .interface(&mesh.mesh_card)
                    .ok_or_else(|| LoweringError::MissingCard {
                        item: format!("gear_mesh:{}", mesh.name),
                        card_kind: "interface",
                        card_hash: mesh.mesh_card,
                    })?;
            let key = GearMeshLoweringKey::new(entry.family.clone());
            let ctor = self.gear_meshes.get(&key).ok_or_else(|| {
                LoweringError::MissingGearMeshConstructor {
                    mesh: mesh.name.clone(),
                    card_family: entry.family.clone(),
                }
            })?;
            let ctx = GearMeshLoweringContext {
                draft,
                mesh,
                card: entry,
                cards,
                selection_policy: self.selection_policy,
            };
            let out = ctor(&ctx, &mut graph)?;
            gear_mesh_receipts.push(GearMeshLoweringReceipt {
                mesh_name: mesh.name.clone(),
                driver: mesh.driver,
                driven: mesh.driven,
                mesh_card: mesh.mesh_card,
                card_family: entry.family.clone(),
                ratio: out.ratio,
                created_nodes: out.created_nodes,
                property_receipts: out.property_receipts,
                parameter_sources: out.parameter_sources,
            });
        }

        // 4. Lower card-bound interfaces.
        for (index, iface) in input.interfaces.iter().enumerate() {
            let kind_tag = InterfaceKindTag::from(&iface.kind);
            let scope = format!("interface:{}", iface.name);
            let (card_hash, card_family, resolved_card, third_body_card) =
                match &iface.kind {
                    DraftInterfaceKind::PowerPort { schema } => {
                        let card = cards.power_port_schema(schema).ok_or_else(|| {
                            LoweringError::MissingCard {
                                item: scope.clone(),
                                card_kind: "power_port_schema",
                                card_hash: *schema,
                            }
                        })?;
                        (
                            Some(*schema),
                            card.family.clone(),
                            ResolvedInterfaceCard::PowerPort(card),
                            None,
                        )
                    }
                    DraftInterfaceKind::MechanicalContact {
                        law_card,
                        third_body,
                        ..
                    }
                    | DraftInterfaceKind::ThermalContact {
                        law_card,
                        third_body,
                        ..
                    }
                    | DraftInterfaceKind::ElectricalContact {
                        law_card,
                        third_body,
                        ..
                    } => {
                        let entry = cards.interface(law_card).ok_or_else(|| {
                            LoweringError::MissingCard {
                                item: scope.clone(),
                                card_kind: "interface",
                                card_hash: *law_card,
                            }
                        })?;
                        (
                            Some(*law_card),
                            entry.family.clone(),
                            ResolvedInterfaceCard::System(entry),
                            third_body.as_ref().map(|tb| tb.material_card),
                        )
                    }
                    DraftInterfaceKind::Glue { law_card } => {
                        let entry = cards.interface(law_card).ok_or_else(|| {
                            LoweringError::MissingCard {
                                item: scope.clone(),
                                card_kind: "interface",
                                card_hash: *law_card,
                            }
                        })?;
                        (
                            Some(*law_card),
                            entry.family.clone(),
                            ResolvedInterfaceCard::System(entry),
                            None,
                        )
                    }
                    DraftInterfaceKind::HeatRoute { .. } => (
                        None,
                        HEAT_ROUTE_FAMILY.to_owned(),
                        ResolvedInterfaceCard::None,
                        None,
                    ),
                    DraftInterfaceKind::ShearFilm { third_body, .. } => (
                        None,
                        SHEAR_FILM_FAMILY.to_owned(),
                        ResolvedInterfaceCard::None,
                        Some(third_body.material_card),
                    ),
                    DraftInterfaceKind::TemperatureSource => (
                        None,
                        TEMPERATURE_SOURCE_FAMILY.to_owned(),
                        ResolvedInterfaceCard::None,
                        None,
                    ),
                    DraftInterfaceKind::PistonFace { .. } => (
                        None,
                        PISTON_FACE_FAMILY.to_owned(),
                        ResolvedInterfaceCard::None,
                        None,
                    ),
                };
            let key = InterfaceLoweringKey::new(kind_tag, card_family.clone());
            let ctor = self.interfaces.get(&key).ok_or_else(|| {
                LoweringError::MissingInterfaceConstructor {
                    interface: iface.name.clone(),
                    kind: kind_tag,
                    card_family: card_family.clone(),
                }
            });
            let ctx = InterfaceLoweringContext {
                draft,
                interface: iface,
                card_hash,
                card: resolved_card,
                cards,
                selection_policy: self.selection_policy,
            };
            let out = match group_interfaces.remove(&index) {
                Some((_, out)) => out,
                None => (ctor?)(&ctx, &mut graph)?,
            };
            interface_receipts.push(InterfaceLoweringReceipt {
                interface_name: iface.name.clone(),
                entity: iface.entity,
                kind: kind_tag,
                card_hash,
                card_family,
                third_body_card,
                created_nodes: out.created_nodes,
                property_receipts: out.property_receipts,
                parameter_sources: out.parameter_sources,
            });
        }

        // 5. Lower declared external sources.
        for (index, src) in input.sources.iter().enumerate() {
            if let Some((_, rec)) = group_sources.remove(&index) {
                source_receipts.push(rec);
                continue;
            }
            if let Some(rec) = graph.consumed_sources.remove(&index) {
                source_receipts.push(rec);
                continue;
            }
            let rec = lower_source(src, &input.catalog, &mut graph)?;
            source_receipts.push(rec);
        }

        // 6. Seal graph (runs unit, port-pairing, and power-closure admission)
        //    and hand to the structural planner.
        let ambient_temp_k = input.environment.base.ambient_temperature.value;
        let SealedLoweringGraph {
            system,
            node_map,
            entity_nodes,
            entity_ports: part_ports,
            observers,
        } = graph.finalize_and_seal(ambient_temp_k)?;
        let probe_dt = match input.time_plan.step_policy {
            StepPolicy::Fixed { dt } => dt.value,
            StepPolicy::Adaptive { dt_initial, .. } => dt_initial.value,
        };
        let plan = system.execution_plan_for_dt(probe_dt);
        let causalization_audit = system.causalization_audit();

        // 7. Assemble and audit receipt completeness.
        let mut receipt = LoweringReceipt {
            schema: LOWERING_RECEIPT_SCHEMA,
            assembly_name: input.name.clone(),
            assembly_hash: draft.content_hash(),
            scenario_ir_hash: draft.scenario_ir_hash(),
            plan_hash: ContentHash(plan.plan_hash),
            receipt_hash: ContentHash([0_u8; 32]),
            components: component_receipts,
            interfaces: interface_receipts,
            joints: joint_receipts,
            gear_meshes: gear_mesh_receipts,
            sources: source_receipts,
            acoustic_members: Vec::new(),
        };
        receipt.verify_completeness()?;
        receipt.finalize_hash();

        Ok(LoweredSystem {
            system,
            plan,
            causalization_audit,
            receipt,
            node_map,
            entity_nodes,
            part_ports,
            observers,
        })
    }

    /// Lower an [`AcousticAssemblySpec`] family by family through the
    /// registered [`AcousticFamilyConstructor`] entries (Do 5).
    ///
    /// # Errors
    /// Refuses missing acoustic family constructors, invalid sample rate,
    /// incomplete parameter provenance, or kernel admission errors.
    pub fn lower_acoustic_spec(
        &self,
        spec: &AcousticAssemblySpec,
        cards: &CardStore,
    ) -> Result<LoweredSystem, LoweringError> {
        if spec.name.trim().is_empty() {
            return Err(LoweringError::Invalid {
                item: "acoustic_spec".to_owned(),
                detail: "assembly name must be non-empty".to_owned(),
            });
        }
        if spec.sample_rate_hz == 0 {
            return Err(LoweringError::Invalid {
                item: format!("acoustic_spec:{}", spec.name),
                detail: "sample_rate_hz must be positive".to_owned(),
            });
        }
        if spec.members.is_empty() {
            return Err(LoweringError::Invalid {
                item: format!("acoustic_spec:{}", spec.name),
                detail: "acoustic spec requires at least one member".to_owned(),
            });
        }
        let mut graph = LoweringGraph::new();
        let mut member_ports = BTreeMap::new();
        let mut acoustic_receipts = Vec::with_capacity(spec.members.len());

        for member in &spec.members {
            let key = AcousticFamilyKey::new(member.family.clone());
            let ctor = self.acoustic_families.get(&key).ok_or_else(|| {
                LoweringError::MissingAcousticFamilyConstructor {
                    member: member.name.clone(),
                    family: member.family.clone(),
                }
            })?;
            let ctx = AcousticFamilyContext {
                spec,
                member,
                cards,
                selection_policy: self.selection_policy,
                member_ports: &member_ports,
            };
            let out = ctor(&ctx, &mut graph)?;
            if let Some(port) = out.exposed_port {
                member_ports.insert(member.name.clone(), port);
            }
            acoustic_receipts.push(AcousticMemberLoweringReceipt {
                member_name: member.name.clone(),
                family: member.family.clone(),
                created_nodes: out.created_nodes,
                property_receipts: out.property_receipts,
                parameter_sources: out.parameter_sources,
            });
        }

        let SealedLoweringGraph {
            system,
            node_map,
            entity_nodes,
            entity_ports: part_ports,
            observers,
        } = graph.finalize_and_seal(spec.ambient_temperature_k)?;
        let dt_s = 1.0 / f64::from(spec.sample_rate_hz);
        let plan = system.execution_plan_for_dt(dt_s);
        let causalization_audit = system.causalization_audit();

        let assembly_hash = hash_domain(
            "org.frankensim.fs-system.lowering.acoustic-spec.v1",
            format!("{spec:?}").as_bytes(),
        );
        let mut receipt = LoweringReceipt {
            schema: LOWERING_RECEIPT_SCHEMA,
            assembly_name: spec.name.clone(),
            assembly_hash,
            scenario_ir_hash: assembly_hash,
            plan_hash: ContentHash(plan.plan_hash),
            receipt_hash: ContentHash([0_u8; 32]),
            components: Vec::new(),
            interfaces: Vec::new(),
            joints: Vec::new(),
            gear_meshes: Vec::new(),
            sources: Vec::new(),
            acoustic_members: acoustic_receipts,
        };
        receipt.verify_completeness()?;
        receipt.finalize_hash();

        Ok(LoweredSystem {
            system,
            plan,
            causalization_audit,
            receipt,
            node_map,
            entity_nodes,
            part_ports,
            observers,
        })
    }

    fn select_part_constructor<'a>(
        &'a self,
        part: &'a DraftPartDecl,
    ) -> Result<(&'a FidelityRung, String, &'a PartConstructor), LoweringError> {
        let mut tried = Vec::with_capacity(part.fidelity.rungs.len());
        let mut matches: Vec<(
            usize,
            &'a FidelityRung,
            FidelityRungKind,
            &'a PartConstructor,
        )> = Vec::new();
        for (idx, rung) in part.fidelity.rungs.iter().enumerate() {
            let rung_kind = FidelityRungKind::from(rung);
            tried.push(rung_kind);
            let key = PartLoweringKey {
                kind: part.kind,
                physics: part.physics.clone(),
                rung: rung_kind,
            };
            if let Some(ctor) = self.parts.get(&key) {
                matches.push((idx, rung, rung_kind, ctor));
            }
        }
        let chosen = match self.fidelity_policy {
            FidelitySelectionPolicy::FirstRegistered => matches.first().copied(),
            FidelitySelectionPolicy::HighestRegistered => matches
                .iter()
                .max_by_key(|(_, _, rk, _)| rk.rank())
                .copied(),
            FidelitySelectionPolicy::LowestRegistered => matches
                .iter()
                .min_by_key(|(_, _, rk, _)| rk.rank())
                .copied(),
        };
        let Some((chosen_idx, chosen_rung, chosen_kind, ctor)) = chosen else {
            return Err(LoweringError::MissingPartConstructor {
                part: part.name.clone(),
                kind: part.kind,
                physics: part.physics.clone(),
                tried_rungs: tried,
            });
        };
        let ladder: Vec<&'static str> = tried.iter().map(|r| r.as_str()).collect();
        let rationale = format!(
            "selected ladder[{chosen_idx}] ({chosen_kind}) via {} from [{}]; relative_tol={:.3e}",
            self.fidelity_policy.as_str(),
            ladder.join(", "),
            part.fidelity.error_budget.relative_tol
        );
        Ok((chosen_rung, rationale, ctor))
    }
}

// ============================================================================
// Source Lowering
// ============================================================================

fn lower_source(
    src: &SourceDecl,
    catalog: &EntityCatalog,
    graph: &mut LoweringGraph,
) -> Result<SourceLoweringReceipt, LoweringError> {
    let scope = format!("source:{}", src.name);
    let mut recorder = ProvenanceRecorder::default();
    let (entity, port_name) = match &src.target {
        SourceTarget::Port { entity, port } => (*entity, port.as_str()),
        SourceTarget::Parameter(path) => {
            return Err(LoweringError::Invalid {
                item: scope,
                detail: format!(
                    "parameter-target source on {}.{} requires a parameter-modulating component; target a port instead",
                    path.entity.token(),
                    path.parameter.name
                ),
            });
        }
    };
    let target_port = graph
        .resolve_entity_port(catalog, entity, port_name)
        .ok_or_else(|| LoweringError::MissingPort {
            item: scope.clone(),
            entity,
            port: port_name.to_owned(),
        })?;
    let (carrier, role, width) =
        graph
            .port_meta
            .get(&target_port)
            .copied()
            .ok_or_else(|| LoweringError::Invalid {
                item: scope.clone(),
                detail: format!("target port '{port_name}' has no registered carrier metadata"),
            })?;
    let target_description = format!("{}.{}", entity.token(), port_name);

    // Record every numeric source parameter explicitly.
    let history_dims = record_time_history_provenance(&scope, &src.history, &mut recorder);
    if carrier == Carrier::Stream {
        return Err(LoweringError::Invalid {
            item: scope,
            detail: format!(
                "scalar source cannot drive stream port '{port_name}': mass, enthalpy and carried quantities require a stream component with explicit coordinate units"
            ),
        });
    }
    if width != 1 {
        return Err(LoweringError::Invalid {
            item: scope,
            detail: format!(
                "scalar source cannot drive {width}-coefficient port '{port_name}'; declare a component with the required field or signal shape"
            ),
        });
    }
    let expected_dims = match (carrier, role) {
        (Carrier::Power(kind), PortRole::EffortIn) => kind.canonical_effort_dimensions(),
        (Carrier::Power(kind), PortRole::FlowIn) => kind.canonical_flow_dimensions(),
        (Carrier::Heat, PortRole::EffortIn) => TEMPERATURE_DIMS,
        (Carrier::Heat, PortRole::FlowIn) => POWER_DIMS,
        // Signal ports currently carry no physical-unit declaration. Only a
        // dimensionless source can enter without discarding authored units.
        (Carrier::Signal, PortRole::SignalIn) => Dims::NONE,
        (_, PortRole::SignalOut) => {
            return Err(LoweringError::Invalid {
                item: scope,
                detail: format!("cannot drive output-only signal port '{port_name}' with a source"),
            });
        }
        _ => {
            return Err(LoweringError::Invalid {
                item: scope,
                detail: format!("port '{port_name}' has no supported scalar source input"),
            });
        }
    };
    if history_dims != expected_dims {
        return Err(LoweringError::Invalid {
            item: scope,
            detail: format!(
                "source history dimensions {:?} do not match port '{port_name}' input dimensions {:?}",
                history_dims.0, expected_dims.0
            ),
        });
    }

    match (carrier, role) {
        (Carrier::Power(kind), PortRole::EffortIn) => {
            // Attach a Force (effort source) to the target's common_flow hub.
            let waveform = time_history_to_waveform(&scope, &src.history)?;
            let node_name = format!("source.{}", src.name);
            let node = graph.add_with_tangent_and_residual(
                &node_name,
                None,
                Force::in_domain(kind, waveform),
                SourceAccount::default(),
            )?;
            let src_port = graph.port(node, "p")?;
            graph.attach_common_flow(target_port, src_port, 1.0);
            let (_, parameter_sources) = recorder.finish();
            Ok(SourceLoweringReceipt {
                source_name: src.name.clone(),
                target_description,
                created_nodes: vec![node_name],
                parameter_sources,
            })
        }
        (Carrier::Heat, PortRole::FlowIn) => {
            // Attach a Joule Heater to the thermal body's common_effort hub.
            let waveform = time_history_to_waveform(&scope, &src.history)?;
            let heater = Heater::new(waveform).map_err(|source| LoweringError::Element {
                item: scope.clone(),
                source,
            })?;
            let node_name = format!("source.{}", src.name);
            let node = graph.add_with_tangent_and_residual(
                &node_name,
                None,
                heater,
                SourceAccount::default(),
            )?;
            let src_port = graph.port(node, "heat")?;
            graph.attach_common_effort(target_port, src_port, 1.0);
            let (_, parameter_sources) = recorder.finish();
            Ok(SourceLoweringReceipt {
                source_name: src.name.clone(),
                target_description,
                created_nodes: vec![node_name],
                parameter_sources,
            })
        }
        (_, PortRole::EffortIn | PortRole::FlowIn | PortRole::SignalIn) => {
            // Prescribe direct scheduled boundary on the target port.
            let (boundary_fn, events) = time_history_to_boundary(&src.history);
            graph
                .builder_mut()
                .scheduled_boundary(target_port, boundary_fn, events)?;
            let (_, parameter_sources) = recorder.finish();
            Ok(SourceLoweringReceipt {
                source_name: src.name.clone(),
                target_description,
                created_nodes: Vec::new(),
                parameter_sources,
            })
        }
        (_, PortRole::SignalOut) => Err(LoweringError::Invalid {
            item: scope,
            detail: format!("cannot drive output-only signal port '{port_name}' with a source"),
        }),
    }
}

fn record_time_history_provenance(
    scope: &str,
    history: &TimeHistory,
    recorder: &mut ProvenanceRecorder,
) -> Dims {
    match history {
        TimeHistory::Constant(v) => {
            recorder.record_declared_field("value", v.value, v.dims, scope, "history.constant");
            v.dims
        }
        TimeHistory::Harmonic {
            amplitude,
            frequency,
            phase,
            offset,
        } => {
            recorder.record_declared_field(
                "amplitude",
                amplitude.value,
                amplitude.dims,
                scope,
                "history.harmonic.amplitude",
            );
            recorder.record_declared_field(
                "frequency",
                frequency.value,
                FREQUENCY_DIMS,
                scope,
                "history.harmonic.frequency",
            );
            recorder.record_declared_field(
                "phase",
                phase.value,
                Dims::NONE,
                scope,
                "history.harmonic.phase",
            );
            recorder.record_declared_field(
                "offset",
                offset.value,
                offset.dims,
                scope,
                "history.harmonic.offset",
            );
            amplitude.dims
        }
        TimeHistory::PiecewiseLinear { knots } => {
            let dims = knots.first().map_or(Dims::NONE, |(_, v)| v.dims);
            for (idx, (t, v)) in knots.iter().enumerate() {
                recorder.record_declared_field(
                    format!("knot[{idx}].t"),
                    t.value,
                    TIME_DIMS,
                    scope,
                    format!("history.piecewise_linear.knots[{idx}].t"),
                );
                recorder.record_declared_field(
                    format!("knot[{idx}].v"),
                    v.value,
                    v.dims,
                    scope,
                    format!("history.piecewise_linear.knots[{idx}].v"),
                );
            }
            dims
        }
        TimeHistory::Step {
            t_switch,
            before,
            after,
        } => {
            recorder.record_declared_field(
                "t_switch",
                t_switch.value,
                TIME_DIMS,
                scope,
                "history.step.t_switch",
            );
            recorder.record_declared_field(
                "before",
                before.value,
                before.dims,
                scope,
                "history.step.before",
            );
            recorder.record_declared_field(
                "after",
                after.value,
                after.dims,
                scope,
                "history.step.after",
            );
            before.dims
        }
    }
}

fn time_history_to_waveform(scope: &str, history: &TimeHistory) -> Result<Waveform, LoweringError> {
    let err = |source| LoweringError::Element {
        item: scope.to_owned(),
        source,
    };
    match history {
        TimeHistory::Constant(v) => Waveform::constant(v.value).map_err(err),
        TimeHistory::Harmonic {
            amplitude,
            frequency,
            phase,
            offset,
        } => {
            Waveform::sine(amplitude.value, frequency.value, phase.value, offset.value).map_err(err)
        }
        TimeHistory::PiecewiseLinear { knots } => {
            let table: Vec<(f64, f64)> = knots.iter().map(|(t, v)| (t.value, v.value)).collect();
            Waveform::table(table).map_err(err)
        }
        TimeHistory::Step {
            t_switch,
            before,
            after,
        } => {
            // The run starts at zero. A switch at zero has already happened;
            // otherwise the held waveform preserves the declared jump and lets
            // the source's ordinary event negotiation end a step at its time.
            if t_switch.value == 0.0 {
                Waveform::constant(after.value).map_err(err)
            } else {
                Waveform::held_table(vec![(0.0, before.value), (t_switch.value, after.value)])
                    .map_err(err)
            }
        }
    }
}

fn time_history_to_boundary(history: &TimeHistory) -> (crate::BoundaryFn, Vec<f64>) {
    let cloned = history.clone();
    let events = match &cloned {
        TimeHistory::Constant(_) | TimeHistory::Harmonic { .. } => Vec::new(),
        TimeHistory::PiecewiseLinear { knots } => knots.iter().map(|(t, _)| t.value).collect(),
        TimeHistory::Step { t_switch, .. } => vec![t_switch.value],
    };
    let func: crate::BoundaryFn = Box::new(move |t_s: f64, dt_s: f64, out: &mut [f64]| {
        let val = evaluate_time_history(&cloned, t_s + 0.5 * dt_s);
        out.fill(val);
    });
    (func, events)
}

fn evaluate_time_history(history: &TimeHistory, t_s: f64) -> f64 {
    match history {
        TimeHistory::Constant(v) => v.value,
        TimeHistory::Harmonic {
            amplitude,
            frequency,
            phase,
            offset,
        } => {
            offset.value + amplitude.value * (2.0 * PI * frequency.value * t_s + phase.value).sin()
        }
        TimeHistory::PiecewiseLinear { knots } => {
            if t_s <= knots[0].0.value {
                return knots[0].1.value;
            }
            let last = knots.len() - 1;
            if t_s >= knots[last].0.value {
                return knots[last].1.value;
            }
            for pair in knots.windows(2) {
                let (t0, v0) = (pair[0].0.value, pair[0].1.value);
                let (t1, v1) = (pair[1].0.value, pair[1].1.value);
                if t_s <= t1 {
                    let frac = (t_s - t0) / (t1 - t0);
                    return v0 + frac * (v1 - v0);
                }
            }
            knots[last].1.value
        }
        TimeHistory::Step {
            t_switch,
            before,
            after,
        } => {
            if t_s < t_switch.value {
                before.value
            } else {
                after.value
            }
        }
    }
}

// ============================================================================
// Built-In Kernel Defaults (Lumped Mechanics, Thermal, Contacts, Joints, Acoustic)
// ============================================================================

/// Populate `registry` with the composition kernel's default lumped constructors
/// and acoustic family adapters.
pub fn register_kernel_defaults(registry: &mut LoweringRegistry) {
    // 1. RigidBody + {Mechanics} + Lumped (translational mass or gear inertia).
    registry.register_part(
        PartLoweringKey::new(
            DraftPartKind::RigidBody,
            [PhysicsDomain::Mechanics],
            FidelityRungKind::Lumped,
        ),
        lower_lumped_rigid_body,
    );

    // 2. RigidBody + {Mechanics, Thermal} + Lumped (coupled mass + heat capacity).
    registry.register_part(
        PartLoweringKey::new(
            DraftPartKind::RigidBody,
            [PhysicsDomain::Mechanics, PhysicsDomain::Thermal],
            FidelityRungKind::Lumped,
        ),
        lower_lumped_rigid_thermal_body,
    );

    // 3. Solid + {Thermal} + Lumped (lumped thermal body).
    registry.register_part(
        PartLoweringKey::new(
            DraftPartKind::Solid,
            [PhysicsDomain::Thermal],
            FidelityRungKind::Lumped,
        ),
        lower_lumped_thermal_solid,
    );

    // 3b. Solid + {Thermal} + OneDField (1-cell 1-D steady conductance).
    registry.register_part(
        PartLoweringKey::new(
            DraftPartKind::Solid,
            [PhysicsDomain::Thermal],
            FidelityRungKind::OneDField,
        ),
        lower_one_d_thermal_conductance,
    );

    // 4. SpringBeam + {Mechanics} + Lumped (axial or spiral torsional spring).
    registry.register_part(
        PartLoweringKey::new(
            DraftPartKind::SpringBeam,
            [PhysicsDomain::Mechanics],
            FidelityRungKind::Lumped,
        ),
        lower_lumped_spring_beam,
    );

    // 4b. DrivenBody + {Mechanics} + Lumped (stationary chassis / ground anchor).
    registry.register_part(
        PartLoweringKey::new(
            DraftPartKind::DrivenBody,
            [PhysicsDomain::Mechanics],
            FidelityRungKind::Lumped,
        ),
        lower_lumped_driven_frame,
    );

    // 5. Joints: Prismatic, Revolute, Fixed.
    registry.register_joint(JointKind::Prismatic, lower_prismatic_joint);
    registry.register_joint(JointKind::Revolute, lower_revolute_joint);
    registry.register_joint(JointKind::Fixed, lower_fixed_joint);

    // 6. Gear mesh: involute-gear-mesh.
    registry.register_gear_mesh("involute-gear-mesh", lower_involute_gear_mesh);

    // 7. Card-bound interfaces.
    let tc = InterfaceKindTag::ThermalContact;
    registry.register_interface(
        InterfaceLoweringKey::new(tc, "thermal-conductance"),
        |c, g| lower_thermal_interface(c, g, false),
    );
    registry.register_interface(
        InterfaceLoweringKey::new(tc, "thermal-contact-resistance"),
        |c, g| lower_thermal_interface(c, g, true),
    );
    registry.register_interface(
        InterfaceLoweringKey::new(InterfaceKindTag::MechanicalContact, "viscous-damper"),
        lower_viscous_damper_interface,
    );
    registry.register_interface(
        InterfaceLoweringKey::new(InterfaceKindTag::MechanicalContact, "kelvin-voigt"),
        lower_kelvin_voigt_interface,
    );
    registry.register_interface(
        InterfaceLoweringKey::new(InterfaceKindTag::Glue, "elastic-glue"),
        lower_elastic_glue_interface,
    );
    registry.register_interface(
        InterfaceLoweringKey::new(InterfaceKindTag::PowerPort, "power-bond"),
        lower_power_bond_interface,
    );

    // 8. AcousticAssembly family entries (Do 5).
    registry.register_acoustic_family("string", lower_acoustic_string_family);
    registry.register_acoustic_family("plate", lower_acoustic_plate_family);
    registry.register_acoustic_family("duct", lower_acoustic_duct_family);
    registry.register_acoustic_family("cavity", lower_acoustic_cavity_family);
}

fn lumped_initial_motion(
    part: &DraftPartDecl,
    scope: &str,
    rot: bool,
    rec: &mut ProvenanceRecorder,
) -> Result<Motion, LoweringError> {
    let Some(init) = &part.initial else {
        return Ok(Motion {
            position_m: 0.0,
            velocity_m_s: 0.0,
        });
    };
    let refuse = |field: &str, val: String| LoweringError::Invalid {
        item: scope.to_owned(),
        detail: format!("lumped 1-DOF rigid body refuses non-zero {field} = {val}"),
    };
    let q = init.orientation;
    if !q[0].is_finite() || (q[0].abs() - 1.0).abs() > 1.0e-12 || q[1..].iter().any(|c| *c != 0.0) {
        return Err(refuse("initial.orientation", format!("{q:?}")));
    }
    let zero = |name: &str, vals: &[QtyAny], skip: Option<usize>| -> Result<(), LoweringError> {
        for (ax, qty) in vals.iter().enumerate() {
            if skip != Some(ax) && qty.value != 0.0 {
                return Err(refuse(
                    &format!("initial.{name}[{ax}]"),
                    qty.value.to_string(),
                ));
            }
        }
        Ok(())
    };
    if rot {
        zero("position", &init.position, None)?;
        zero("linear_velocity", &init.linear_velocity, None)?;
        zero("angular_velocity", &init.angular_velocity, Some(2))?;
        let wz = init.angular_velocity[2].value;
        rec.record_declared_field(
            "initial_angular_velocity_rad_s",
            wz,
            FREQUENCY_DIMS,
            scope,
            "initial.angular_velocity[2]",
        );
        return Ok(Motion {
            position_m: 0.0,
            velocity_m_s: wz,
        });
    }
    zero("position", &init.position, Some(0))?;
    zero("linear_velocity", &init.linear_velocity, Some(0))?;
    zero("angular_velocity", &init.angular_velocity, None)?;
    let (x, vx) = (init.position[0].value, init.linear_velocity[0].value);
    rec.record_declared_field(
        "initial_position_m",
        x,
        LENGTH_DIMS,
        scope,
        "initial.position[0]",
    );
    rec.record_declared_field(
        "initial_velocity_m_s",
        vx,
        VELOCITY_DIMS,
        scope,
        "initial.linear_velocity[0]",
    );
    Ok(Motion {
        position_m: x,
        velocity_m_s: vx,
    })
}

fn lumped_initial_temp(
    ctx: &PartLoweringContext<'_>,
    scope: &str,
    rec: &mut ProvenanceRecorder,
) -> f64 {
    let (t, src, field) = match ctx.part.initial.as_ref().and_then(|i| i.temperature) {
        Some(t) => (t.value, scope, "initial.temperature"),
        None => (
            ctx.ambient_temperature_k(),
            "environment",
            "ambient_temperature",
        ),
    };
    rec.record_declared_field("initial_temperature_k", t, TEMPERATURE_DIMS, src, field);
    t
}

fn lower_lumped_rigid_body(
    ctx: &PartLoweringContext<'_>,
    graph: &mut LoweringGraph,
) -> Result<LoweredPartOutput, LoweringError> {
    let scope = format!("part:{}", ctx.part.name);
    let mut recorder = ProvenanceRecorder::default();
    let (rho, _) = ctx.query_material(&mut recorder, 0, "density", "density")?;
    let gen_label = geometry_generator_label(&ctx.part.geometry);

    let (component, node_name, port_name) =
        if let DraftGeometrySource::Gear(_) = &ctx.part.geometry {
            let (i_over_rho, formula) = geometry_polar_inertia_over_rho(&ctx.part.geometry)
                .ok_or_else(|| LoweringError::Invalid {
                    item: scope.clone(),
                    detail: "cannot compute rotational inertia for gear geometry".to_owned(),
                })?;
            let inertia = rho.value * i_over_rho;
            recorder.record_geometry_derived(
                "moment_of_inertia",
                inertia,
                INERTIA_DIMS,
                ctx.geometry_fingerprint,
                &gen_label,
                formula,
            );
            let comp = Mass::in_domain(PortKind::RotationalTorqueAngularVelocity, inertia)
                .map_err(|source| LoweringError::Element {
                    item: scope.clone(),
                    source,
                })?;
            (comp, ctx.part.name.clone(), "rotation")
        } else {
            let (vol, vol_formula) =
                geometry_volume_m3(&ctx.part.geometry).ok_or_else(|| LoweringError::Invalid {
                    item: scope.clone(),
                    detail: "lumped rigid body requires closed-form volume geometry".to_owned(),
                })?;
            let mass_kg = rho.value * vol;
            recorder.record_geometry_derived(
                "volume",
                vol,
                VOLUME_DIMS,
                ctx.geometry_fingerprint,
                &gen_label,
                vol_formula,
            );
            recorder.record_geometry_derived(
                "mass",
                mass_kg,
                MASS_DIMS,
                ctx.geometry_fingerprint,
                &gen_label,
                "mass = density * volume",
            );
            let comp = Mass::new(mass_kg).map_err(|source| LoweringError::Element {
                item: scope.clone(),
                source,
            })?;
            (comp, ctx.part.name.clone(), "mechanics")
        };

    let init = lumped_initial_motion(ctx.part, &scope, port_name == "rotation", &mut recorder)?;
    let node =
        graph.add_with_tangent_and_residual(&node_name, Some(ctx.part.entity), component, init)?;
    let p = graph.port(node, "p")?;
    graph.expose_port(ctx.part.entity, port_name, p);
    if port_name == "rotation" {
        graph.expose_port(ctx.part.entity, "mechanics", p);
    }

    let (property_receipts, parameter_sources) = recorder.finish();
    Ok(LoweredPartOutput {
        created_nodes: vec![node_name],
        property_receipts,
        parameter_sources,
    })
}

fn lower_lumped_rigid_thermal_body(
    ctx: &PartLoweringContext<'_>,
    graph: &mut LoweringGraph,
) -> Result<LoweredPartOutput, LoweringError> {
    let scope = format!("part:{}", ctx.part.name);
    let mut recorder = ProvenanceRecorder::default();
    let (rho, _) = ctx.query_material(&mut recorder, 0, "density", "density")?;
    let (cp, _) = ctx.query_material(
        &mut recorder,
        0,
        "specific_heat_capacity",
        "specific_heat_capacity",
    )?;
    let gen_label = geometry_generator_label(&ctx.part.geometry);
    let (vol, vol_formula) =
        geometry_volume_m3(&ctx.part.geometry).ok_or_else(|| LoweringError::Invalid {
            item: scope.clone(),
            detail: "lumped thermo-mechanical rigid body requires closed-form volume geometry"
                .to_owned(),
        })?;
    let mass_kg = rho.value * vol;
    let capacity_j_k = mass_kg * cp.value;

    recorder.record_geometry_derived(
        "volume",
        vol,
        VOLUME_DIMS,
        ctx.geometry_fingerprint,
        &gen_label,
        vol_formula,
    );
    recorder.record_geometry_derived(
        "mass",
        mass_kg,
        MASS_DIMS,
        ctx.geometry_fingerprint,
        &gen_label,
        "mass = density * volume",
    );
    recorder.record_geometry_derived(
        "heat_capacity",
        capacity_j_k,
        HEAT_CAPACITY_DIMS,
        ctx.geometry_fingerprint,
        &gen_label,
        "heat_capacity = mass * specific_heat_capacity",
    );
    let init = lumped_initial_motion(ctx.part, &scope, false, &mut recorder)?;
    let t_init = lumped_initial_temp(ctx, &scope, &mut recorder);

    let mech_name = format!("{}.mechanics", ctx.part.name);
    let therm_name = format!("{}.thermal", ctx.part.name);
    let mass = Mass::new(mass_kg).map_err(|source| LoweringError::Element {
        item: scope.clone(),
        source,
    })?;
    let cap = HeatCapacity::new(capacity_j_k).map_err(|source| LoweringError::Element {
        item: scope,
        source,
    })?;

    let mech_node =
        graph.add_with_tangent_and_residual(&mech_name, Some(ctx.part.entity), mass, init)?;
    let therm_node =
        graph.add_with_tangent_and_residual(&therm_name, Some(ctx.part.entity), cap, t_init)?;

    let mech_port = graph.port(mech_node, "p")?;
    let therm_port = graph.port(therm_node, "heat")?;
    graph.expose_port(ctx.part.entity, "mechanics", mech_port);
    graph.expose_port(ctx.part.entity, "thermal", therm_port);

    let (property_receipts, parameter_sources) = recorder.finish();
    Ok(LoweredPartOutput {
        created_nodes: vec![mech_name, therm_name],
        property_receipts,
        parameter_sources,
    })
}

fn lower_lumped_thermal_solid(
    ctx: &PartLoweringContext<'_>,
    graph: &mut LoweringGraph,
) -> Result<LoweredPartOutput, LoweringError> {
    let scope = format!("part:{}", ctx.part.name);
    let mut recorder = ProvenanceRecorder::default();
    let (rho, _) = ctx.query_material(&mut recorder, 0, "density", "density")?;
    let (cp, _) = ctx.query_material(
        &mut recorder,
        0,
        "specific_heat_capacity",
        "specific_heat_capacity",
    )?;
    let gen_label = geometry_generator_label(&ctx.part.geometry);
    let (vol, vol_formula) =
        geometry_volume_m3(&ctx.part.geometry).ok_or_else(|| LoweringError::Invalid {
            item: scope.clone(),
            detail: "lumped thermal solid requires closed-form volume geometry".to_owned(),
        })?;
    let capacity_j_k = rho.value * vol * cp.value;

    recorder.record_geometry_derived(
        "volume",
        vol,
        VOLUME_DIMS,
        ctx.geometry_fingerprint,
        &gen_label,
        vol_formula,
    );
    recorder.record_geometry_derived(
        "heat_capacity",
        capacity_j_k,
        HEAT_CAPACITY_DIMS,
        ctx.geometry_fingerprint,
        &gen_label,
        "heat_capacity = density * volume * specific_heat_capacity",
    );
    let t_init = lumped_initial_temp(ctx, &scope, &mut recorder);

    let cap = HeatCapacity::new(capacity_j_k).map_err(|source| LoweringError::Element {
        item: scope,
        source,
    })?;
    let node =
        graph.add_with_tangent_and_residual(&ctx.part.name, Some(ctx.part.entity), cap, t_init)?;
    let port = graph.port(node, "heat")?;
    graph.expose_port(ctx.part.entity, "thermal", port);

    let (property_receipts, parameter_sources) = recorder.finish();
    Ok(LoweredPartOutput {
        created_nodes: vec![ctx.part.name.clone()],
        property_receipts,
        parameter_sources,
    })
}

fn lower_one_d_thermal_conductance(
    ctx: &PartLoweringContext<'_>,
    graph: &mut LoweringGraph,
) -> Result<LoweredPartOutput, LoweringError> {
    let scope = format!("part:{}", ctx.part.name);
    let inv = |detail| LoweringError::Invalid {
        item: scope.clone(),
        detail,
    };
    if *ctx.chosen_rung != (FidelityRung::OneDField { cells: 1 }) {
        return Err(inv(format!(
            "requires (one-d-field 1), got {:?}",
            ctx.chosen_rung
        )));
    }
    let DraftGeometrySource::Primitive(PrimitiveShape::Box { lx, ly, lz }) = &ctx.part.geometry
    else {
        let g = geometry_generator_label(&ctx.part.geometry);
        return Err(inv(format!("requires box geometry, got `{g}`")));
    };
    let prop = ctx
        .part
        .materials
        .first()
        .and_then(|m| ctx.cards.material(&m.card))
        .filter(|c| !c.claims().claims_for("thermal-conductivity").is_empty())
        .map_or("thermal_conductivity", |_| "thermal-conductivity");
    let mut rec = ProvenanceRecorder::default();
    let (k, _) = ctx.query_material(&mut rec, 0, prop, "thermal_conductivity")?;
    if k.dims != THERMAL_CONDUCTIVITY_DIMS {
        return Err(inv(format!(
            "conductivity needs W/(m*K), got {:?}",
            k.dims.0
        )));
    }
    let (aol, gf, gl) = (
        (ly.value * lz.value) / lx.value,
        ctx.geometry_fingerprint,
        "Primitive.Box.draft.v1",
    );
    let g_w_k = k.value * aol;
    rec.record_geometry_derived("area_over_length", aol, LENGTH_DIMS, gf, gl, "(ly*lz)/lx");
    rec.record_geometry_derived(
        "thermal_conductance",
        g_w_k,
        CONDUCTANCE_DIMS,
        gf,
        gl,
        "k*A/L",
    );
    let cond = Conductance::new(g_w_k).map_err(|source| LoweringError::Element {
        item: scope,
        source,
    })?;
    let node =
        graph.add_with_tangent_and_residual(&ctx.part.name, Some(ctx.part.entity), cond, ())?;
    let (pa, pb) = (graph.port(node, "a")?, graph.port(node, "b")?);
    graph.expose_port(ctx.part.entity, "a", pa);
    graph.expose_port(ctx.part.entity, "b", pb);
    graph.observers.insert(
        ctx.part.entity,
        Box::new(move |sys| sys.last_output(pa).first().copied().map(Observed::HeatRate)),
    );
    let (property_receipts, parameter_sources) = rec.finish();
    Ok(LoweredPartOutput {
        created_nodes: vec![ctx.part.name.clone()],
        property_receipts,
        parameter_sources,
    })
}

fn lower_lumped_spring_beam(
    ctx: &PartLoweringContext<'_>,
    graph: &mut LoweringGraph,
) -> Result<LoweredPartOutput, LoweringError> {
    let scope = format!("part:{}", ctx.part.name);
    let mut recorder = ProvenanceRecorder::default();
    let (modulus, _) = ctx.query_material(&mut recorder, 0, "youngs_modulus", "youngs_modulus")?;
    let gen_label = geometry_generator_label(&ctx.part.geometry);

    let spring = if let DraftGeometrySource::SpiralSpring(s) = &ctx.part.geometry {
        let factor = spiral_spring_stiffness_factor_m3(s);
        let k_theta = modulus.value * factor;
        recorder.record_geometry_derived(
            "torsional_stiffness",
            k_theta,
            TORSIONAL_STIFFNESS_DIMS,
            ctx.geometry_fingerprint,
            &gen_label,
            "k_theta = E * blade_height * blade_thickness^3 / (12 * 2 * pi * r_mean * turns)",
        );
        Spring::in_domain(PortKind::RotationalTorqueAngularVelocity, k_theta, 0.0).map_err(
            |source| LoweringError::Element {
                item: scope.clone(),
                source,
            },
        )?
    } else {
        let (area_over_len, formula) = geometry_axial_area_over_length_m(&ctx.part.geometry)
            .ok_or_else(|| LoweringError::Invalid {
                item: scope.clone(),
                detail: "lumped spring_beam requires Box, Cylinder, Tube, or SpiralSpring geometry"
                    .to_owned(),
            })?;
        let k = modulus.value * area_over_len;
        recorder.record_geometry_derived(
            "stiffness",
            k,
            STIFFNESS_DIMS,
            ctx.geometry_fingerprint,
            &gen_label,
            formula,
        );
        Spring::new(k, 0.0).map_err(|source| LoweringError::Element {
            item: scope.clone(),
            source,
        })?
    };

    let node =
        graph.add_with_tangent_and_residual(&ctx.part.name, Some(ctx.part.entity), spring, 0.0)?;
    let port_a = graph.port(node, "a")?;
    let port_b = graph.port(node, "b")?;
    graph.expose_port(ctx.part.entity, "a", port_a);
    graph.expose_port(ctx.part.entity, "b", port_b);

    let (property_receipts, parameter_sources) = recorder.finish();
    Ok(LoweredPartOutput {
        created_nodes: vec![ctx.part.name.clone()],
        property_receipts,
        parameter_sources,
    })
}

fn lower_lumped_driven_frame(
    ctx: &PartLoweringContext<'_>,
    graph: &mut LoweringGraph,
) -> Result<LoweredPartOutput, LoweringError> {
    let (scope, mut rec) = (
        format!("part:{}", ctx.part.name),
        ProvenanceRecorder::default(),
    );
    let drives: Vec<_> = ctx.draft.input().sources.iter().enumerate().filter(|(_, s)| {
        matches!(&s.target, SourceTarget::Port { entity, port } if *entity == ctx.part.entity && port == "mechanics")
    }).collect();
    if drives.is_empty() {
        rec.record_declared_field(
            "prescribed_velocity_m_s",
            0.0,
            VELOCITY_DIMS,
            &scope,
            "driven_body.stationary_frame",
        );
        graph.mark_ground_entity(ctx.part.entity);
        let (property_receipts, parameter_sources) = rec.finish();
        return Ok(LoweredPartOutput {
            created_nodes: vec![],
            property_receipts,
            parameter_sources,
        });
    }
    if drives.len() > 1 {
        return Err(LoweringError::Invalid {
            item: scope,
            detail: "driven body accepts at most one mechanics velocity source".into(),
        });
    }
    let (src_idx, src) = drives[0];
    let (src_scope, mut src_rec) = (
        format!("source:{}", src.name),
        ProvenanceRecorder::default(),
    );
    let dims = record_time_history_provenance(&src_scope, &src.history, &mut src_rec);
    if dims != VELOCITY_DIMS {
        return Err(LoweringError::Invalid {
            item: src_scope,
            detail: format!(
                "driven body velocity source needs dims [1,0,-1,0,0,0], got {:?}",
                dims.0
            ),
        });
    }
    let wf = time_history_to_waveform(&src_scope, &src.history)?;
    let v0 = evaluate_time_history(&src.history, 0.0);
    rec.record_declared_field(
        "initial_velocity_m_s",
        v0,
        VELOCITY_DIMS,
        &src_scope,
        "driven_body.velocity_t0",
    );
    let node = graph.add_with_tangent_and_residual(
        &ctx.part.name,
        Some(ctx.part.entity),
        FlowSource::new(wf).publishing(),
        SourceAccount::default(),
    )?;
    graph.expose_port(ctx.part.entity, "mechanics", graph.port(node, "p")?);
    graph.expose_port(ctx.part.entity, "span", graph.port(node, "span")?);
    let (_, src_params) = src_rec.finish();
    graph.consumed_sources.insert(
        src_idx,
        SourceLoweringReceipt {
            source_name: src.name.clone(),
            target_description: format!("{}.mechanics", ctx.part.entity.token()),
            created_nodes: vec![ctx.part.name.clone()],
            parameter_sources: src_params,
        },
    );
    let (property_receipts, parameter_sources) = rec.finish();
    Ok(LoweredPartOutput {
        created_nodes: vec![ctx.part.name.clone()],
        property_receipts,
        parameter_sources,
    })
}

fn lower_prismatic_joint(
    ctx: &JointLoweringContext<'_>,
    graph: &mut LoweringGraph,
) -> Result<LoweredJointOutput, LoweringError> {
    lower_1dof_joint(
        ctx,
        graph,
        PortKind::MechanicalForceVelocity,
        "mechanics",
        STIFFNESS_DIMS,
        DAMPING_DIMS,
    )
}

fn lower_revolute_joint(
    ctx: &JointLoweringContext<'_>,
    graph: &mut LoweringGraph,
) -> Result<LoweredJointOutput, LoweringError> {
    lower_1dof_joint(
        ctx,
        graph,
        PortKind::RotationalTorqueAngularVelocity,
        "rotation",
        TORSIONAL_STIFFNESS_DIMS,
        TORSIONAL_DAMPING_DIMS,
    )
}

fn lower_fixed_joint(
    ctx: &JointLoweringContext<'_>,
    graph: &mut LoweringGraph,
) -> Result<LoweredJointOutput, LoweringError> {
    let scope = format!("joint:{}", ctx.joint.name);
    let mut recorder = ProvenanceRecorder::default();
    let stiffness = ctx.joint.stiffness.ok_or_else(|| LoweringError::Invalid {
        item: scope.clone(),
        detail:
            "lumped fixed joint requires explicit finite stiffness to realize as a compliant tie"
                .to_owned(),
    })?;
    recorder.record_declared_field(
        "stiffness",
        stiffness.value,
        stiffness.dims,
        &scope,
        "stiffness",
    );
    let child_port = graph
        .entity_port(ctx.joint.child, "mechanics")
        .ok_or_else(|| LoweringError::MissingPort {
            item: scope.clone(),
            entity: ctx.joint.child,
            port: "mechanics".to_owned(),
        })?;
    let spring = Spring::new(stiffness.value, 0.0).map_err(|source| LoweringError::Element {
        item: scope.clone(),
        source,
    })?;
    let node_name = format!("joint.{}.tie", ctx.joint.name);
    let node = graph.add_with_tangent_and_residual(&node_name, None, spring, 0.0)?;
    let a = graph.port(node, "a")?;
    let b = graph.port(node, "b")?;
    graph.attach_common_flow(child_port, a, 1.0);
    if ctx.joint.parent == ctx.draft.input().root {
        graph.bind_zero_boundary(b);
    } else {
        let parent_port = graph
            .entity_port(ctx.joint.parent, "mechanics")
            .ok_or_else(|| LoweringError::MissingPort {
                item: scope,
                entity: ctx.joint.parent,
                port: "mechanics".to_owned(),
            })?;
        graph.attach_common_flow(parent_port, b, 1.0);
    }
    let (_, parameter_sources) = recorder.finish();
    Ok(LoweredJointOutput {
        created_nodes: vec![node_name],
        parameter_sources,
    })
}

fn lower_1dof_joint(
    ctx: &JointLoweringContext<'_>,
    graph: &mut LoweringGraph,
    domain: PortKind,
    preferred_port: &str,
    stiffness_dims: Dims,
    damping_dims: Dims,
) -> Result<LoweredJointOutput, LoweringError> {
    let scope = format!("joint:{}", ctx.joint.name);
    let mut recorder = ProvenanceRecorder::default();
    let mut created_nodes = Vec::new();
    let root = ctx.draft.input().root;

    let child_port = graph
        .entity_port(ctx.joint.child, preferred_port)
        .or_else(|| graph.entity_port(ctx.joint.child, "mechanics"))
        .ok_or_else(|| LoweringError::MissingPort {
            item: scope.clone(),
            entity: ctx.joint.child,
            port: preferred_port.to_owned(),
        })?;
    let parent_port = if ctx.joint.parent == root {
        None
    } else {
        Some(
            graph
                .entity_port(ctx.joint.parent, preferred_port)
                .or_else(|| graph.entity_port(ctx.joint.parent, "mechanics"))
                .ok_or_else(|| LoweringError::MissingPort {
                    item: scope.clone(),
                    entity: ctx.joint.parent,
                    port: preferred_port.to_owned(),
                })?,
        )
    };

    if let Some(k) = ctx.joint.stiffness
        && k.value > 0.0
    {
        recorder.record_declared_field("stiffness", k.value, stiffness_dims, &scope, "stiffness");
        let spring =
            Spring::in_domain(domain, k.value, 0.0).map_err(|source| LoweringError::Element {
                item: scope.clone(),
                source,
            })?;
        let name = format!("joint.{}.spring", ctx.joint.name);
        let node = graph.add_with_tangent_and_residual(&name, None, spring, 0.0)?;
        let a = graph.port(node, "a")?;
        let b = graph.port(node, "b")?;
        graph.attach_common_flow(child_port, a, 1.0);
        match parent_port {
            Some(p) => graph.attach_common_flow(p, b, 1.0),
            None => graph.bind_zero_boundary(b),
        }
        created_nodes.push(name);
    }

    if let Some(c) = ctx.joint.damping
        && c.value > 0.0
    {
        recorder.record_declared_field("damping", c.value, damping_dims, &scope, "damping");
        let damper =
            Damper::in_domain(domain, c.value).map_err(|source| LoweringError::Element {
                item: scope.clone(),
                source,
            })?;
        let name = format!("joint.{}.damper", ctx.joint.name);
        let node = graph.add_with_tangent_and_residual(&name, None, damper, ())?;
        let a = graph.port(node, "a")?;
        let b = graph.port(node, "b")?;
        let heat = graph.port(node, "heat")?;
        graph.attach_common_flow(child_port, a, 1.0);
        match parent_port {
            Some(p) => graph.attach_common_flow(p, b, 1.0),
            None => graph.bind_zero_boundary(b),
        }
        route_dissipated_heat(
            graph,
            &ctx.draft.input().catalog,
            root,
            heat,
            ctx.joint.child,
            (ctx.joint.parent != root).then_some(ctx.joint.parent),
            ctx.draft.input().environment.base.ambient_temperature.value,
            &scope,
            &mut recorder,
        )?;
        created_nodes.push(name);
    }

    if created_nodes.is_empty() {
        return Err(LoweringError::Invalid {
            item: scope,
            detail: "lumped 1-DOF joint requires positive stiffness and/or damping to lower into kernel elements".to_owned(),
        });
    }

    let (_, parameter_sources) = recorder.finish();
    Ok(LoweredJointOutput {
        created_nodes,
        parameter_sources,
    })
}

fn lower_involute_gear_mesh(
    ctx: &GearMeshLoweringContext<'_>,
    graph: &mut LoweringGraph,
) -> Result<LoweredGearMeshOutput, LoweringError> {
    let scope = format!("gear_mesh:{}", ctx.mesh.name);
    let mut recorder = ProvenanceRecorder::default();
    let ratio = -(f64::from(ctx.mesh.teeth_driver) / f64::from(ctx.mesh.teeth_driven));
    recorder.record_declared_field(
        "teeth_driver",
        f64::from(ctx.mesh.teeth_driver),
        Dims::NONE,
        &scope,
        "teeth_driver",
    );
    recorder.record_declared_field(
        "teeth_driven",
        f64::from(ctx.mesh.teeth_driven),
        Dims::NONE,
        &scope,
        "teeth_driven",
    );
    recorder.record_declared_field(
        "center_distance",
        ctx.mesh.center_distance.value,
        LENGTH_DIMS,
        &scope,
        "center_distance",
    );

    let driver_port = graph
        .entity_port(ctx.mesh.driver, "rotation")
        .or_else(|| graph.entity_port(ctx.mesh.driver, "mechanics"))
        .ok_or_else(|| LoweringError::MissingPort {
            item: scope.clone(),
            entity: ctx.mesh.driver,
            port: "rotation".to_owned(),
        })?;
    let driven_port = graph
        .entity_port(ctx.mesh.driven, "rotation")
        .or_else(|| graph.entity_port(ctx.mesh.driven, "mechanics"))
        .ok_or_else(|| LoweringError::MissingPort {
            item: scope.clone(),
            entity: ctx.mesh.driven,
            port: "rotation".to_owned(),
        })?;

    let (k_mesh, _) = recorder.query_interface_property(
        &scope,
        ctx.card,
        &[],
        "mesh_stiffness",
        "mesh_stiffness",
        MatDbQueryEnv {
            ambient_temp_k: ctx.draft.input().environment.base.ambient_temperature.value,
            ambient_pressure_pa: ctx.draft.input().environment.base.ambient_pressure.value,
            policy: ctx.selection_policy,
        },
    )?;

    // Compliant gear mesh: torsional spring with ratio on the driven side so
    // power cancels identically across the teeth.
    let spring = Spring::in_domain(PortKind::RotationalTorqueAngularVelocity, k_mesh.value, 0.0)
        .map_err(|source| LoweringError::Element {
            item: scope.clone(),
            source,
        })?;
    let name = format!("gear_mesh.{}.compliance", ctx.mesh.name);
    let node = graph.add_with_tangent_and_residual(&name, None, spring, 0.0)?;
    let a = graph.port(node, "a")?;
    let b = graph.port(node, "b")?;
    graph.attach_common_flow(driver_port, a, 1.0);
    graph.attach_common_flow(driven_port, b, ratio);

    let (property_receipts, parameter_sources) = recorder.finish();
    Ok(LoweredGearMeshOutput {
        ratio,
        created_nodes: vec![name],
        property_receipts,
        parameter_sources,
    })
}

/// A `thermal-conductance` card's G [W/K], or `A/R''` of a resistance card at
/// `:area`, recorded in `rec`: the kernel's card-only lumped contact.
pub(crate) fn query_thermal_conductance(
    scope: &str,
    rec: &mut ProvenanceRecorder,
    entry: &InterfaceCardEntry,
    state: &[(fs_scenario::assembly::ParameterKey, fs_qty::QtyAny)],
    resistance_law: bool,
    env: MatDbQueryEnv,
) -> Result<f64, LoweringError> {
    let inv = |d: &str| LoweringError::Invalid {
        item: scope.to_owned(),
        detail: d.to_owned(),
    };
    let (area, prop) = if resistance_law {
        let a = state
            .iter()
            .find(|(k, _)| k.name == "area")
            .filter(|(_, q)| q.dims == AREA_DIMS && q.value.is_finite() && q.value > 0.0)
            .ok_or_else(|| inv("thermal-contact-resistance requires positive :area [m^2]"))?;
        (Some(a.1.value), "area-specific-thermal-contact-resistance")
    } else {
        (None, "thermal_conductance")
    };
    let s = rec
        .query_interface_property(scope, entry, state, prop, prop, env)?
        .0;
    let Some(a) = area else {
        return Ok(s.value);
    };
    let g = a / s.value;
    if s.dims != CONTACT_RESISTANCE_DIMS || s.value <= 0.0 || !g.is_finite() {
        return Err(inv("contact resistance must be positive [m^2*K/W]"));
    }
    rec.record_declared_field("area", a, AREA_DIMS, scope, "state.area");
    rec.record_declared_field("thermal_conductance", g, CONDUCTANCE_DIMS, scope, "A/R''");
    Ok(g)
}

fn wire_thermal_ends(
    ctx: &InterfaceLoweringContext<'_>,
    graph: &mut LoweringGraph,
    recorder: &mut ProvenanceRecorder,
    ports: [PortRef; 2],
    scope: &str,
) -> Result<(), LoweringError> {
    let (cat, root, t_amb) = (
        &ctx.draft.input().catalog,
        ctx.draft.input().root,
        ctx.ambient_temperature_k(),
    );
    let resolve = |feat| {
        if let Some(p) = graph.resolve_entity_port(cat, feat, "thermal") {
            Ok(Some(p))
        } else if graph.is_ground_feature(cat, feat, root) {
            Ok(None)
        } else {
            Err(LoweringError::MissingPort {
                item: scope.to_owned(),
                entity: feat,
                port: "thermal".to_owned(),
            })
        }
    };
    let (pa, pb) = (
        resolve(ctx.interface.feature_a)?,
        resolve(ctx.interface.feature_b)?,
    );
    if pa.is_none() && pb.is_none() {
        return Err(LoweringError::Invalid {
            item: scope.to_owned(),
            detail: "thermal interface cannot connect two ground features".to_owned(),
        });
    }
    for (target, port) in [(pa, ports[0]), (pb, ports[1])] {
        if let Some(t) = target {
            graph.attach_common_effort(t, port, 1.0);
        } else {
            recorder.record_declared_field(
                "ambient_bath_temperature_k",
                t_amb,
                TEMPERATURE_DIMS,
                "environment",
                "ambient_temperature",
            );
            graph.add_heat_route(
                port,
                Vec::new(),
                Some(ambient_heat_sink(scope, t_amb)?),
                RouteTemperature::Fixed(t_amb),
            );
        }
    }
    Ok(())
}

fn lower_thermal_interface(
    ctx: &InterfaceLoweringContext<'_>,
    graph: &mut LoweringGraph,
    resistance_law: bool,
) -> Result<LoweredInterfaceOutput, LoweringError> {
    let scope = format!("interface:{}", ctx.interface.name);
    let ResolvedInterfaceCard::System(entry) = ctx.card else {
        return Err(LoweringError::Invalid {
            item: scope,
            detail: "thermal-conductance requires an InterfaceSystemCard".to_owned(),
        });
    };
    let (state, third_body) = match &ctx.interface.kind {
        DraftInterfaceKind::ThermalContact {
            state, third_body, ..
        } => (state.as_slice(), third_body.as_ref()),
        _ => (&[][..], None),
    };
    let (mut recorder, env) = (ProvenanceRecorder::default(), ctx.query_env());
    let mut g_eff =
        query_thermal_conductance(&scope, &mut recorder, entry, state, resistance_law, env)?;
    if let Some(tb) = third_body {
        // Optional third-body TIM contributes series thermal resistance when
        // thermal_conductivity is present on the third-body card.
        let (k_tim, _) = recorder.query_third_body_property(
            &scope,
            ctx.cards,
            tb,
            "thermal_conductivity",
            "third_body_thermal_conductivity",
            ctx.query_env(),
        )?;
        recorder.record_declared_field(
            "third_body_film_thickness",
            tb.film_thickness.value,
            LENGTH_DIMS,
            &scope,
            "third_body.film_thickness",
        );
        // Combine per-unit-area film conductance k_tim / thickness in series
        // using contact_area from interface state if provided, else 1 m^2 nominal.
        let area = state
            .iter()
            .find(|(k, _)| k.name == "contact_area")
            .map_or(1.0, |(_, q)| {
                recorder.record_declared_field(
                    "contact_area",
                    q.value,
                    q.dims,
                    &scope,
                    "state.contact_area",
                );
                q.value
            });
        let g_tim = k_tim.value * area / tb.film_thickness.value;
        g_eff = (g_eff * g_tim) / (g_eff + g_tim);
        recorder.record_declared_field(
            "effective_thermal_conductance",
            g_eff,
            CONDUCTANCE_DIMS,
            &scope,
            "series_combination(thermal_conductance, k_tim * A / t)",
        );
    }

    let cond = Conductance::new(g_eff).map_err(|source| LoweringError::Element {
        item: scope.clone(),
        source,
    })?;
    let node_name = format!("interface.{}", ctx.interface.name);
    let node =
        graph.add_with_tangent_and_residual(&node_name, Some(ctx.interface.entity), cond, ())?;
    let ports = [graph.port(node, "a")?, graph.port(node, "b")?];
    wire_thermal_ends(ctx, graph, &mut recorder, ports, &scope)?;

    let (property_receipts, parameter_sources) = recorder.finish();
    Ok(LoweredInterfaceOutput {
        created_nodes: vec![node_name],
        property_receipts,
        parameter_sources,
    })
}

fn lower_viscous_damper_interface(
    ctx: &InterfaceLoweringContext<'_>,
    graph: &mut LoweringGraph,
) -> Result<LoweredInterfaceOutput, LoweringError> {
    let scope = format!("interface:{}", ctx.interface.name);
    let ResolvedInterfaceCard::System(entry) = ctx.card else {
        return Err(LoweringError::Invalid {
            item: scope,
            detail: "viscous-damper requires an InterfaceSystemCard".to_owned(),
        });
    };
    let state = match &ctx.interface.kind {
        DraftInterfaceKind::MechanicalContact { state, .. } => state.as_slice(),
        _ => &[],
    };
    let mut recorder = ProvenanceRecorder::default();
    let (c_sample, _) = recorder.query_interface_property(
        &scope,
        entry,
        state,
        "damping_coefficient",
        "damping_coefficient",
        ctx.query_env(),
    )?;

    let catalog = &ctx.draft.input().catalog;
    let root = ctx.draft.input().root;
    let damper = Damper::new(c_sample.value).map_err(|source| LoweringError::Element {
        item: scope.clone(),
        source,
    })?;
    let node_name = format!("interface.{}", ctx.interface.name);
    let node =
        graph.add_with_tangent_and_residual(&node_name, Some(ctx.interface.entity), damper, ())?;
    let a = graph.port(node, "a")?;
    let b = graph.port(node, "b")?;
    let heat = graph.port(node, "heat")?;

    attach_mech_endpoint(graph, catalog, ctx.interface.feature_a, root, a, &scope)?;
    attach_mech_endpoint(graph, catalog, ctx.interface.feature_b, root, b, &scope)?;
    route_dissipated_heat(
        graph,
        catalog,
        root,
        heat,
        ctx.interface.feature_a,
        Some(ctx.interface.feature_b),
        ctx.ambient_temperature_k(),
        &scope,
        &mut recorder,
    )?;

    let (property_receipts, parameter_sources) = recorder.finish();
    Ok(LoweredInterfaceOutput {
        created_nodes: vec![node_name],
        property_receipts,
        parameter_sources,
    })
}

fn lower_kelvin_voigt_interface(
    ctx: &InterfaceLoweringContext<'_>,
    graph: &mut LoweringGraph,
) -> Result<LoweredInterfaceOutput, LoweringError> {
    let scope = format!("interface:{}", ctx.interface.name);
    let ResolvedInterfaceCard::System(entry) = ctx.card else {
        return Err(LoweringError::Invalid {
            item: scope,
            detail: "kelvin-voigt requires an InterfaceSystemCard".to_owned(),
        });
    };
    let state = match &ctx.interface.kind {
        DraftInterfaceKind::MechanicalContact { state, .. } => state.as_slice(),
        _ => &[],
    };
    let mut recorder = ProvenanceRecorder::default();
    let (k_sample, _) = recorder.query_interface_property(
        &scope,
        entry,
        state,
        "contact_stiffness",
        "contact_stiffness",
        ctx.query_env(),
    )?;
    let (c_sample, _) = recorder.query_interface_property(
        &scope,
        entry,
        state,
        "damping_coefficient",
        "damping_coefficient",
        ctx.query_env(),
    )?;

    let catalog = &ctx.draft.input().catalog;
    let root = ctx.draft.input().root;
    let spring = Spring::new(k_sample.value, 0.0).map_err(|source| LoweringError::Element {
        item: scope.clone(),
        source,
    })?;
    let damper = Damper::new(c_sample.value).map_err(|source| LoweringError::Element {
        item: scope.clone(),
        source,
    })?;
    let s_name = format!("interface.{}.spring", ctx.interface.name);
    let d_name = format!("interface.{}.damper", ctx.interface.name);
    let s_node =
        graph.add_with_tangent_and_residual(&s_name, Some(ctx.interface.entity), spring, 0.0)?;
    let d_node =
        graph.add_with_tangent_and_residual(&d_name, Some(ctx.interface.entity), damper, ())?;

    let (sa, sb) = (graph.port(s_node, "a")?, graph.port(s_node, "b")?);
    let (da, db, d_heat) = (
        graph.port(d_node, "a")?,
        graph.port(d_node, "b")?,
        graph.port(d_node, "heat")?,
    );
    attach_mech_endpoint(graph, catalog, ctx.interface.feature_a, root, sa, &scope)?;
    attach_mech_endpoint(graph, catalog, ctx.interface.feature_b, root, sb, &scope)?;
    attach_mech_endpoint(graph, catalog, ctx.interface.feature_a, root, da, &scope)?;
    attach_mech_endpoint(graph, catalog, ctx.interface.feature_b, root, db, &scope)?;
    route_dissipated_heat(
        graph,
        catalog,
        root,
        d_heat,
        ctx.interface.feature_a,
        Some(ctx.interface.feature_b),
        ctx.ambient_temperature_k(),
        &scope,
        &mut recorder,
    )?;

    let (property_receipts, parameter_sources) = recorder.finish();
    Ok(LoweredInterfaceOutput {
        created_nodes: vec![s_name, d_name],
        property_receipts,
        parameter_sources,
    })
}

fn lower_elastic_glue_interface(
    ctx: &InterfaceLoweringContext<'_>,
    graph: &mut LoweringGraph,
) -> Result<LoweredInterfaceOutput, LoweringError> {
    let scope = format!("interface:{}", ctx.interface.name);
    let ResolvedInterfaceCard::System(entry) = ctx.card else {
        return Err(LoweringError::Invalid {
            item: scope,
            detail: "elastic-glue requires an InterfaceSystemCard".to_owned(),
        });
    };
    let mut recorder = ProvenanceRecorder::default();
    let (k_sample, _) = recorder.query_interface_property(
        &scope,
        entry,
        &[],
        "bond_stiffness",
        "bond_stiffness",
        ctx.query_env(),
    )?;

    let catalog = &ctx.draft.input().catalog;
    let root = ctx.draft.input().root;
    let spring = Spring::new(k_sample.value, 0.0).map_err(|source| LoweringError::Element {
        item: scope.clone(),
        source,
    })?;
    let node_name = format!("interface.{}", ctx.interface.name);
    let node =
        graph.add_with_tangent_and_residual(&node_name, Some(ctx.interface.entity), spring, 0.0)?;
    let a = graph.port(node, "a")?;
    let b = graph.port(node, "b")?;
    attach_mech_endpoint(graph, catalog, ctx.interface.feature_a, root, a, &scope)?;
    attach_mech_endpoint(graph, catalog, ctx.interface.feature_b, root, b, &scope)?;

    let (property_receipts, parameter_sources) = recorder.finish();
    Ok(LoweredInterfaceOutput {
        created_nodes: vec![node_name],
        property_receipts,
        parameter_sources,
    })
}

fn lower_power_bond_interface(
    ctx: &InterfaceLoweringContext<'_>,
    graph: &mut LoweringGraph,
) -> Result<LoweredInterfaceOutput, LoweringError> {
    let scope = format!("interface:{}", ctx.interface.name);
    let ResolvedInterfaceCard::PowerPort(schema) = ctx.card else {
        return Err(LoweringError::Invalid {
            item: scope,
            detail: "power-bond requires a PowerPortSchemaCard".to_owned(),
        });
    };
    let mut recorder = ProvenanceRecorder::default();
    recorder.record_declared_field(
        "bond_ratio",
        1.0,
        Dims::NONE,
        format!("power_port_schema:{}", schema.schema_name),
        "unit_power_bond",
    );

    let catalog = &ctx.draft.input().catalog;
    let root = ctx.draft.input().root;
    let port_a = first_power_port(graph, catalog, ctx.interface.feature_a).ok_or_else(|| {
        LoweringError::MissingPort {
            item: scope.clone(),
            entity: ctx.interface.feature_a,
            port: "mechanics|a".to_owned(),
        }
    })?;
    if graph.is_ground_feature(catalog, ctx.interface.feature_b, root) {
        graph.bind_zero_boundary(port_a);
    } else {
        let port_b =
            first_power_port(graph, catalog, ctx.interface.feature_b).ok_or_else(|| {
                LoweringError::MissingPort {
                    item: scope.clone(),
                    entity: ctx.interface.feature_b,
                    port: "mechanics|b".to_owned(),
                }
            })?;
        match (graph.port_meta(port_a), graph.port_meta(port_b)) {
            (Some((_, PortRole::EffortIn)), Some((_, PortRole::FlowIn))) => {
                graph.attach_common_flow(port_a, port_b, 1.0);
            }
            (Some((_, PortRole::FlowIn)), Some((_, PortRole::EffortIn))) => {
                graph.attach_common_flow(port_b, port_a, 1.0);
            }
            _ => {
                graph.builder_mut().bond(port_a, port_b)?;
            }
        }
    }

    let (property_receipts, parameter_sources) = recorder.finish();
    Ok(LoweredInterfaceOutput {
        created_nodes: Vec::new(),
        property_receipts,
        parameter_sources,
    })
}

fn first_power_port(
    graph: &LoweringGraph,
    catalog: &EntityCatalog,
    entity: EntityId,
) -> Option<PortRef> {
    graph
        .resolve_entity_port(catalog, entity, "mechanics")
        .or_else(|| graph.resolve_entity_port(catalog, entity, "rotation"))
        .or_else(|| graph.resolve_entity_port(catalog, entity, "a"))
        .or_else(|| graph.resolve_entity_port(catalog, entity, "b"))
}

fn attach_mech_endpoint(
    graph: &mut LoweringGraph,
    catalog: &EntityCatalog,
    entity: EntityId,
    root: EntityId,
    member_port: PortRef,
    scope: &str,
) -> Result<(), LoweringError> {
    if graph.is_ground_feature(catalog, entity, root) {
        graph.bind_zero_boundary(member_port);
        return Ok(());
    }
    let target = graph
        .resolve_entity_port(catalog, entity, "mechanics")
        .or_else(|| graph.resolve_entity_port(catalog, entity, "rotation"))
        .ok_or_else(|| LoweringError::MissingPort {
            item: scope.to_owned(),
            entity,
            port: "mechanics".to_owned(),
        })?;
    graph.attach_common_flow(target, member_port, 1.0);
    Ok(())
}

fn ambient_heat_sink(scope: &str, kelvin: f64) -> Result<HeatSink, LoweringError> {
    heat_sink("ambient", kelvin).map_err(|err| LoweringError::Invalid {
        item: scope.to_owned(),
        detail: format!("invalid ambient heat sink temperature {kelvin} K: {err:?}"),
    })
}

#[allow(clippy::too_many_arguments)]
fn route_dissipated_heat(
    graph: &mut LoweringGraph,
    catalog: &EntityCatalog,
    root: EntityId,
    heat_source: PortRef,
    entity_a: EntityId,
    entity_b: Option<EntityId>,
    ambient_temp_k: f64,
    scope: &str,
    recorder: &mut ProvenanceRecorder,
) -> Result<(), LoweringError> {
    let therm_a = if graph.is_ground_feature(catalog, entity_a, root) {
        None
    } else {
        graph.resolve_entity_port(catalog, entity_a, "thermal")
    };
    let therm_b = entity_b
        .filter(|&e| !graph.is_ground_feature(catalog, e, root))
        .and_then(|e| graph.resolve_entity_port(catalog, e, "thermal"));
    match (therm_a, therm_b) {
        (Some(ta), None) => {
            // Single thermal body: attach via common_effort so it composes
            // cleanly with any Conductances/Heaters on the same body.
            graph.attach_common_effort(ta, heat_source, 1.0);
        }
        (None, Some(tb)) => {
            graph.attach_common_effort(tb, heat_source, 1.0);
        }
        (Some(ta), Some(tb)) => {
            // Two thermal bodies: split 50/50 via explicit heat_route.
            recorder.record_declared_field(
                "heat_partition_a",
                0.5,
                Dims::NONE,
                "lowering.heat_route",
                "symmetric_two_body_split",
            );
            graph.add_heat_route(
                heat_source,
                vec![(ta, 0.5), (tb, 0.5)],
                None,
                RouteTemperature::SinkEquivalent,
            );
        }
        (None, None) => {
            // No thermal body on either side: route explicitly to named ambient
            // HeatSink so entropy and energy close without silent drops.
            recorder.record_declared_field(
                "ambient_sink_temperature_k",
                ambient_temp_k,
                TEMPERATURE_DIMS,
                "environment",
                "ambient_temperature",
            );
            graph.add_heat_route(
                heat_source,
                Vec::new(),
                Some(ambient_heat_sink(scope, ambient_temp_k)?),
                RouteTemperature::Fixed(ambient_temp_k),
            );
        }
    }
    Ok(())
}

// ============================================================================
// AcousticAssembly Family Constructors (Do 5)
// ============================================================================

fn acoustic_param(
    member: &AcousticMemberSpec,
    name: &str,
    recorder: &mut ProvenanceRecorder,
) -> Result<f64, LoweringError> {
    let scope = format!("acoustic:{}", member.name);
    let &(value, dims) = member
        .parameters
        .get(name)
        .ok_or_else(|| LoweringError::Invalid {
            item: scope.clone(),
            detail: format!(
                "acoustic member '{}' (family '{}') requires parameter '{name}'",
                member.name, member.family
            ),
        })?;
    recorder.record_declared_field(name, value, dims, &scope, name);
    Ok(value)
}

fn lower_acoustic_string_family(
    ctx: &AcousticFamilyContext<'_>,
    graph: &mut LoweringGraph,
) -> Result<LoweredAcousticMemberOutput, LoweringError> {
    let scope = format!("acoustic:{}", ctx.member.name);
    let mut recorder = ProvenanceRecorder::default();
    let length_m = acoustic_param(ctx.member, "length_m", &mut recorder)?;
    let tension_n = acoustic_param(ctx.member, "tension_n", &mut recorder)?;
    let linear_density = acoustic_param(ctx.member, "linear_density_kg_m", &mut recorder)?;
    let damping = acoustic_param(ctx.member, "modal_damping_n_s_m", &mut recorder)?;

    // Fundamental pinned-string modal mass (0.5 * mu * L) and stiffness (pi^2 * T / (2 * L)).
    let m_eff = 0.5 * linear_density * length_m;
    let k_eff = 0.5 * PI * PI * tension_n / length_m;
    recorder.record_declared_field(
        "modal_mass",
        m_eff,
        MASS_DIMS,
        &scope,
        "0.5 * linear_density_kg_m * length_m",
    );
    recorder.record_declared_field(
        "modal_stiffness",
        k_eff,
        STIFFNESS_DIMS,
        &scope,
        "0.5 * pi^2 * tension_n / length_m",
    );

    let m_name = format!("{}.mass", ctx.member.name);
    let s_name = format!("{}.spring", ctx.member.name);
    let d_name = format!("{}.damper", ctx.member.name);
    let m_node = graph.add_with_tangent_and_residual(
        &m_name,
        None,
        Mass::new(m_eff).map_err(|source| LoweringError::Element {
            item: scope.clone(),
            source,
        })?,
        Motion {
            position_m: 0.0,
            velocity_m_s: 0.0,
        },
    )?;
    let s_node = graph.add_with_tangent_and_residual(
        &s_name,
        None,
        Spring::new(k_eff, 0.0).map_err(|source| LoweringError::Element {
            item: scope.clone(),
            source,
        })?,
        0.0,
    )?;
    let d_node = graph.add_with_tangent_and_residual(
        &d_name,
        None,
        Damper::new(damping).map_err(|source| LoweringError::Element {
            item: scope.clone(),
            source,
        })?,
        (),
    )?;

    let m_p = graph.port(m_node, "p")?;
    let (sa, sb) = (graph.port(s_node, "a")?, graph.port(s_node, "b")?);
    let (da, db, d_heat) = (
        graph.port(d_node, "a")?,
        graph.port(d_node, "b")?,
        graph.port(d_node, "heat")?,
    );
    graph.attach_common_flow(m_p, sa, 1.0);
    graph.attach_common_flow(m_p, da, 1.0);
    if let Some(target_name) = &ctx.member.coupled_to {
        let target_port =
            ctx.member_ports
                .get(target_name)
                .copied()
                .ok_or_else(|| LoweringError::Invalid {
                    item: scope.clone(),
                    detail: format!(
                        "coupled acoustic member '{target_name}' has not been lowered yet"
                    ),
                })?;
        graph.attach_common_flow(target_port, sb, 1.0);
        graph.attach_common_flow(target_port, db, 1.0);
    } else {
        graph.bind_zero_boundary(sb);
        graph.bind_zero_boundary(db);
    }
    graph.add_heat_route(
        d_heat,
        Vec::new(),
        Some(ambient_heat_sink(&scope, ctx.spec.ambient_temperature_k)?),
        RouteTemperature::Fixed(ctx.spec.ambient_temperature_k),
    );

    let (property_receipts, parameter_sources) = recorder.finish();
    Ok(LoweredAcousticMemberOutput {
        created_nodes: vec![m_name, s_name, d_name],
        exposed_port: Some(m_p),
        property_receipts,
        parameter_sources,
    })
}

fn lower_acoustic_plate_family(
    ctx: &AcousticFamilyContext<'_>,
    graph: &mut LoweringGraph,
) -> Result<LoweredAcousticMemberOutput, LoweringError> {
    let scope = format!("acoustic:{}", ctx.member.name);
    let mut recorder = ProvenanceRecorder::default();
    let modal_mass = acoustic_param(ctx.member, "modal_mass_kg", &mut recorder)?;
    let modal_freq = acoustic_param(ctx.member, "modal_frequency_hz", &mut recorder)?;
    let modal_damping = acoustic_param(ctx.member, "modal_damping_n_s_m", &mut recorder)?;
    let omega = 2.0 * PI * modal_freq;
    let k_modal = modal_mass * omega * omega;
    recorder.record_declared_field(
        "modal_stiffness",
        k_modal,
        STIFFNESS_DIMS,
        &scope,
        "modal_mass_kg * (2 * pi * modal_frequency_hz)^2",
    );

    let m_name = format!("{}.modal_mass", ctx.member.name);
    let s_name = format!("{}.modal_spring", ctx.member.name);
    let d_name = format!("{}.modal_damper", ctx.member.name);
    let m_node = graph.add_with_tangent_and_residual(
        &m_name,
        None,
        Mass::new(modal_mass).map_err(|source| LoweringError::Element {
            item: scope.clone(),
            source,
        })?,
        Motion {
            position_m: 0.0,
            velocity_m_s: 0.0,
        },
    )?;
    let s_node = graph.add_with_tangent_and_residual(
        &s_name,
        None,
        Spring::new(k_modal, 0.0).map_err(|source| LoweringError::Element {
            item: scope.clone(),
            source,
        })?,
        0.0,
    )?;
    let d_node = graph.add_with_tangent_and_residual(
        &d_name,
        None,
        Damper::new(modal_damping).map_err(|source| LoweringError::Element {
            item: scope.clone(),
            source,
        })?,
        (),
    )?;

    let m_p = graph.port(m_node, "p")?;
    let (sa, sb) = (graph.port(s_node, "a")?, graph.port(s_node, "b")?);
    let (da, db, d_heat) = (
        graph.port(d_node, "a")?,
        graph.port(d_node, "b")?,
        graph.port(d_node, "heat")?,
    );
    graph.attach_common_flow(m_p, sa, 1.0);
    graph.attach_common_flow(m_p, da, 1.0);
    graph.bind_zero_boundary(sb);
    graph.bind_zero_boundary(db);
    graph.add_heat_route(
        d_heat,
        Vec::new(),
        Some(ambient_heat_sink(&scope, ctx.spec.ambient_temperature_k)?),
        RouteTemperature::Fixed(ctx.spec.ambient_temperature_k),
    );

    let (property_receipts, parameter_sources) = recorder.finish();
    Ok(LoweredAcousticMemberOutput {
        created_nodes: vec![m_name, s_name, d_name],
        exposed_port: Some(m_p),
        property_receipts,
        parameter_sources,
    })
}

fn lower_acoustic_duct_family(
    ctx: &AcousticFamilyContext<'_>,
    graph: &mut LoweringGraph,
) -> Result<LoweredAcousticMemberOutput, LoweringError> {
    let scope = format!("acoustic:{}", ctx.member.name);
    let mut recorder = ProvenanceRecorder::default();
    let inertance = acoustic_param(ctx.member, "acoustic_inertance_kg_m4", &mut recorder)?;
    let resistance = acoustic_param(ctx.member, "viscothermal_resistance", &mut recorder)?;

    let l_name = format!("{}.inertance", ctx.member.name);
    let r_name = format!("{}.resistance", ctx.member.name);
    let l_node = graph.add_with_tangent_and_residual(
        &l_name,
        None,
        Mass::in_domain(PortKind::FluidPressureFlux, inertance).map_err(|source| {
            LoweringError::Element {
                item: scope.clone(),
                source,
            }
        })?,
        Motion {
            position_m: 0.0,
            velocity_m_s: 0.0,
        },
    )?;
    let r_node = graph.add_with_tangent_and_residual(
        &r_name,
        None,
        Damper::in_domain(PortKind::FluidPressureFlux, resistance).map_err(|source| {
            LoweringError::Element {
                item: scope.clone(),
                source,
            }
        })?,
        (),
    )?;
    let l_p = graph.port(l_node, "p")?;
    let (ra, rb, r_heat) = (
        graph.port(r_node, "a")?,
        graph.port(r_node, "b")?,
        graph.port(r_node, "heat")?,
    );
    graph.attach_common_flow(l_p, ra, 1.0);
    graph.bind_zero_boundary(rb);
    graph.add_heat_route(
        r_heat,
        Vec::new(),
        Some(ambient_heat_sink(&scope, ctx.spec.ambient_temperature_k)?),
        RouteTemperature::Fixed(ctx.spec.ambient_temperature_k),
    );
    if let Some(target_name) = &ctx.member.coupled_to {
        let target_port =
            ctx.member_ports
                .get(target_name)
                .copied()
                .ok_or_else(|| LoweringError::Invalid {
                    item: scope,
                    detail: format!(
                        "coupled acoustic member '{target_name}' has not been lowered yet"
                    ),
                })?;
        // Duct inertance `l_p` (EffortIn) joins series FlowIn compliance/cavity
        // and resistance ports on one common_flow junction.
        graph.attach_common_flow(l_p, target_port, 1.0);
    }

    let (property_receipts, parameter_sources) = recorder.finish();
    Ok(LoweredAcousticMemberOutput {
        created_nodes: vec![l_name, r_name],
        exposed_port: Some(l_p),
        property_receipts,
        parameter_sources,
    })
}

fn lower_acoustic_cavity_family(
    ctx: &AcousticFamilyContext<'_>,
    graph: &mut LoweringGraph,
) -> Result<LoweredAcousticMemberOutput, LoweringError> {
    let scope = format!("acoustic:{}", ctx.member.name);
    let mut recorder = ProvenanceRecorder::default();
    let compliance = acoustic_param(ctx.member, "acoustic_compliance_m5_n", &mut recorder)?;

    let c_name = format!("{}.compliance", ctx.member.name);
    let c_node = graph.add_with_tangent_and_residual(
        &c_name,
        None,
        Storage::in_domain(PortKind::FluidPressureFlux, compliance).map_err(|source| {
            LoweringError::Element {
                item: scope,
                source,
            }
        })?,
        0.0,
    )?;
    let c_p = graph.port(c_node, "p")?;

    let (property_receipts, parameter_sources) = recorder.finish();
    Ok(LoweredAcousticMemberOutput {
        created_nodes: vec![c_name],
        exposed_port: Some(c_p),
        property_receipts,
        parameter_sources,
    })
}

// ============================================================================
// Canonical JSON & Hashing Helpers (Zero External Dependencies)
// ============================================================================

fn push_len_prefixed(buf: &mut Vec<u8>, bytes: &[u8]) {
    buf.extend_from_slice(&(bytes.len() as u64).to_le_bytes());
    buf.extend_from_slice(bytes);
}

fn push_json_str(out: &mut String, s: &str) {
    out.push('"');
    for ch in s.chars() {
        match ch {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            c if (c as u32) < 0x20 => {
                use std::fmt::Write;
                let _ = write!(out, "\\u{:04x}", c as u32);
            }
            c => out.push(c),
        }
    }
    out.push('"');
}

fn push_json_f64(out: &mut String, v: f64) {
    use std::fmt::Write;
    if v.is_finite() {
        let _ = write!(out, "{v:.17e}");
    } else {
        out.push_str("null");
    }
}

fn push_json_str_slice(out: &mut String, items: &[String]) {
    out.push('[');
    for (idx, item) in items.iter().enumerate() {
        if idx > 0 {
            out.push(',');
        }
        push_json_str(out, item);
    }
    out.push(']');
}

fn push_json_dims(out: &mut String, dims: Dims) {
    use std::fmt::Write;
    out.push('[');
    for (idx, d) in dims.0.iter().enumerate() {
        if idx > 0 {
            out.push(',');
        }
        let _ = write!(out, "{d}");
    }
    out.push(']');
}

fn push_json_fidelity(out: &mut String, rung: &FidelityRung) {
    use std::fmt::Write;
    out.push_str("{\"kind\":");
    push_json_str(out, FidelityRungKind::from(rung).as_str());
    match rung {
        FidelityRung::Lumped => {}
        FidelityRung::Modal { mode_count } => {
            let _ = write!(out, ",\"mode_count\":{mode_count}");
        }
        FidelityRung::ReducedOrder { basis_dim } => {
            let _ = write!(out, ",\"basis_dim\":{basis_dim}");
        }
        FidelityRung::OneDField { cells } => {
            let _ = write!(out, ",\"cells\":{cells}");
        }
        FidelityRung::FiniteElement { order, h_target } => {
            let _ = write!(out, ",\"order\":{order},\"h_target_m\":");
            push_json_f64(out, h_target.value);
        }
        FidelityRung::LatticeBoltzmann { dx } => {
            out.push_str(",\"dx_m\":");
            push_json_f64(out, dx.value);
        }
        FidelityRung::BoundaryElement { panels } => {
            let _ = write!(out, ",\"panels\":{panels}");
        }
    }
    out.push('}');
}

fn push_json_property_receipts(out: &mut String, receipts: &[PropertyUsageReceipt]) {
    use std::fmt::Write;
    out.push('[');
    for (idx, r) in receipts.iter().enumerate() {
        if idx > 0 {
            out.push(',');
        }
        out.push_str("{\"receipt_hash\":");
        push_json_str(out, &r.content_hash().to_hex());
        out.push_str(",\"property\":");
        push_json_str(out, &r.property);
        out.push_str(",\"selected_claim\":");
        push_json_str(out, &r.selected.0.to_hex());
        out.push_str(",\"policy\":");
        push_json_str(out, r.policy);
        let _ = write!(
            out,
            ",\"observation_backed\":{},\"evaluator_version\":{},\"considered_count\":{},\"in_domain_count\":{},\"query_point\":{{",
            r.observation_backed,
            r.evaluator_version,
            r.considered.len(),
            r.in_domain.len()
        );
        for (q_idx, (axis, val)) in r.query_point.iter().enumerate() {
            if q_idx > 0 {
                out.push(',');
            }
            push_json_str(out, axis);
            out.push(':');
            push_json_f64(out, *val);
        }
        out.push_str("}}");
    }
    out.push(']');
}

fn push_json_param_sources(out: &mut String, params: &[ParameterProvenance]) {
    use std::fmt::Write;
    out.push('[');
    for (idx, p) in params.iter().enumerate() {
        if idx > 0 {
            out.push(',');
        }
        out.push_str("{\"parameter\":");
        push_json_str(out, &p.parameter);
        out.push_str(",\"value\":");
        push_json_f64(out, p.value);
        out.push_str(",\"dims\":");
        push_json_dims(out, p.dims);
        out.push_str(",\"source\":");
        match &p.source {
            ParameterSource::MatDbClaim {
                card_hash,
                property,
                receipt_hash,
                selected_claim,
                observation_backed,
            } => {
                out.push_str("{\"kind\":\"matdb_claim\",\"card_hash\":");
                push_json_str(out, &card_hash.to_hex());
                out.push_str(",\"property\":");
                push_json_str(out, property);
                out.push_str(",\"receipt_hash\":");
                push_json_str(out, &receipt_hash.to_hex());
                out.push_str(",\"selected_claim\":");
                push_json_str(out, &selected_claim.to_hex());
                let _ = write!(out, ",\"observation_backed\":{observation_backed}}}");
            }
            ParameterSource::ConstitutiveModelParam {
                card_hash,
                model_hash,
                law,
                law_version,
                parameter,
            } => {
                out.push_str("{\"kind\":\"constitutive_model_param\",\"card_hash\":");
                push_json_str(out, &card_hash.to_hex());
                out.push_str(",\"model_hash\":");
                push_json_str(out, &model_hash.to_hex());
                out.push_str(",\"law\":");
                push_json_str(out, law);
                let _ = write!(out, ",\"law_version\":{law_version},\"parameter\":");
                push_json_str(out, parameter);
                out.push('}');
            }
            ParameterSource::GeometryDerived {
                geometry_fingerprint,
                generator,
                formula,
            } => {
                out.push_str("{\"kind\":\"geometry_derived\",\"geometry_fingerprint\":");
                push_json_str(out, &geometry_fingerprint.to_hex());
                out.push_str(",\"generator\":");
                push_json_str(out, generator);
                out.push_str(",\"formula\":");
                push_json_str(out, formula);
                out.push('}');
            }
            ParameterSource::DeclaredAssemblyField { item, field } => {
                out.push_str("{\"kind\":\"declared_assembly_field\",\"item\":");
                push_json_str(out, item);
                out.push_str(",\"field\":");
                push_json_str(out, field);
                out.push('}');
            }
        }
        out.push('}');
    }
    out.push(']');
}

#[cfg(test)]
mod source_admission_tests {
    use super::*;
    use crate::elements::{GasVolume, SignalLimiter};
    use fs_qty::QtyAny;
    use fs_scenario::EntityDeclaration;

    #[test]
    fn scalar_sources_preserve_signal_units_and_refuse_vector_or_stream_broadcast() {
        let mut catalog = EntityCatalog::new();
        let root = catalog
            .declare(EntityDeclaration::assembly("sources"))
            .unwrap();
        let entity = catalog
            .declare(EntityDeclaration::part(root, "receiver"))
            .unwrap();
        let mut graph = LoweringGraph::new();
        for (name, width) in [("scalar", 1), ("vector", 2)] {
            let node = graph
                .add_with_tangent(
                    name,
                    Some(entity),
                    SignalLimiter::new(vec![0.0; width], vec![1.0; width]).unwrap(),
                    (),
                )
                .unwrap();
            let port = graph.port(node, "input").unwrap();
            graph.expose_port(entity, name, port);
        }
        let gas = GasVolume::new(287.0, 1.4, 0).unwrap();
        let contents = gas.contents(101_325.0, 300.0, 1.0, &[]).unwrap();
        let node = graph
            .add_with_tangent("gas", Some(entity), gas, contents)
            .unwrap();
        let port = graph.port(node, "stream").unwrap();
        graph.expose_port(entity, "stream", port);

        for (port, dims, reason) in [
            ("scalar", TEMPERATURE_DIMS, "dimensions"),
            ("vector", Dims::NONE, "2-coefficient"),
            ("stream", Dims([0, 1, -1, 0, 0, 0]), "stream component"),
        ] {
            let source = SourceDecl {
                name: "drive".to_owned(),
                target: SourceTarget::Port {
                    entity,
                    port: port.to_owned(),
                },
                history: TimeHistory::Constant(QtyAny::new(0.5, dims)),
            };
            let error = lower_source(&source, &catalog, &mut graph)
                .expect_err("scalar source cannot discard dimensions or invent a vector shape");
            assert!(
                matches!(error, LoweringError::Invalid { ref item, ref detail }
                if item == "source:drive" && detail.contains(port) && detail.contains(reason))
            );
        }
        let source = SourceDecl {
            name: "dimensionless".to_owned(),
            target: SourceTarget::Port {
                entity,
                port: "scalar".to_owned(),
            },
            history: TimeHistory::Constant(QtyAny::dimensionless(0.5)),
        };
        lower_source(&source, &catalog, &mut graph)
            .expect("dimensionless scalar signal remains supported");
    }
}
