//! Feature-gated I01.3 lowering from neutral [`fs_couple::PortSchema`] and
//! [`fs_couple::StreamPort`] declarations into the multi-field
//! [`crate::system`] IR.
//!
//! This module is intentionally a type-and-accounting compiler, not a numeric
//! interface solver. It re-derives the power dimensions, binds the complete
//! port schema into the resulting [`crate::SystemId`], makes orientation
//! reversal an explicit algebraic sign, and refuses source/loss terms without
//! exactly one ownership disposition. Numeric contraction, quadrature, port
//! adapter truth, and closed-window conservation remain external obligations.

use std::collections::BTreeMap;
use std::num::NonZeroUsize;

use fs_couple::{
    BoundaryTreatment, ConservationRole, ConservativeJunction, DissipationEvidence, DissipationLaw,
    DissipativeRelation, FieldMeasureSide, PORT_SCHEMA_VERSION, PortKind, PortOrientation,
    PortSchema, PortValueShape, PowerPairing, STREAM_PORT_VERSION, SourceClass, SourceOrReservoir,
    StableId, StorageElement, StoragePotential, StreamConstituentId, StreamEnergyChart, StreamPort,
    StreamStressWorkConvention,
};
use fs_iface::SpaceType;
use fs_qty::Dims;

use crate::expr::Space;
use crate::system::{
    AdmittedSystem, AtomSignature, BlockEquation, ConventionRef, CoordinateConvention, FieldDecl,
    FieldQuantity, MAX_SYSTEM_EXTENSION_BYTES, SpatialSupport, StateOwnership, SystemDef,
    SystemExpr, SystemId, SystemTypeError,
};

/// Canonical schema label for [`PortEquationReceipt::to_json`].
pub const PORT_EQUATION_RECEIPT_SCHEMA_V1: &str = "fs-opdsl-port-equation-receipt-v1";

/// Canonical schema label for [`StreamEquationReceipt::to_json`].
pub const STREAM_EQUATION_RECEIPT_SCHEMA_V1: &str = "fs-opdsl-stream-equation-receipt-v1";

/// Maximum port equations admitted in one deterministic batch.
pub const MAX_PORT_EQUATIONS: usize = 4_096;

const RAW_VECTOR_DEGREE: u8 = 255;
const POWER_DIMS: Dims = Dims([2, 1, -3, 0, 0, 0]);
const MASS_FLOW_DIMS: Dims = Dims([0, 1, -1, 0, 0, 0]);
const AMOUNT_FLOW_DIMS: Dims = Dims([0, 0, -1, 0, 0, 1]);
const MOMENTUM_FLOW_DIMS: Dims = Dims([1, 1, -2, 0, 0, 0]);
const ENTROPY_FLOW_DIMS: Dims = Dims([2, 1, -3, -1, 0, 0]);
const PORT_EXTENSION_VERSION: u32 = 1;
const PRIMITIVE_BINDING_VERSION: u32 = 1;
const STREAM_EXTENSION_VERSION: u32 = 1;
const STORAGE_ACTION_BINDING_VERSION: u32 = 1;
const LOSS_OWNERSHIP_DOMAIN_V1: &str = "org.frankensim.fs-opdsl.loss-ownership.v1";

/// Nominal content identity for one concretely owned dissipative term.
///
/// The compiler derives it from the complete source schema, discretization,
/// dissipative role, and declared concrete owner after normalizing only the
/// algebraic equation sense. Reversal therefore does not invent a second
/// physical loss owner, while any physical/schema change moves the identity.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct LossOwnershipId(fs_blake3::ContentHash);

impl LossOwnershipId {
    /// Exact digest bytes.
    #[must_use]
    pub fn as_bytes(&self) -> &[u8; 32] {
        self.0.as_bytes()
    }

    /// Lowercase hexadecimal rendering.
    #[must_use]
    pub fn to_hex(self) -> String {
        self.0.to_hex()
    }

    /// Parse an exact 64-digit hexadecimal transport. Parsing adds no
    /// authority; callers compare it with a freshly compiled receipt.
    #[must_use]
    pub fn parse_hex(value: &str) -> Option<Self> {
        fs_blake3::ContentHash::from_hex(value).map(Self)
    }
}

impl core::fmt::Display for LossOwnershipId {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        core::fmt::Display::fmt(&self.0, f)
    }
}

/// Whether the compiled equation follows or reverses the schema's declared
/// positive orientation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PortEquationSense {
    /// Preserve the declared positive sense.
    AsDeclared,
    /// Apply the explicit negative of the declared power contribution.
    Reversed,
}

impl PortEquationSense {
    /// Exact multiplier inserted into the generated expression.
    #[must_use]
    pub const fn sign(self) -> i8 {
        match self {
            Self::AsDeclared => 1,
            Self::Reversed => -1,
        }
    }

    const fn as_str(self) -> &'static str {
        match self {
            Self::AsDeclared => "as-declared",
            Self::Reversed => "reversed",
        }
    }
}

/// Accounting role of one generated power term.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AccountingTermKind {
    /// Lossless/reversible interface exchange. No source or loss owner exists.
    Reversible,
    /// Stored-energy contribution. A concrete owner is mandatory.
    Storage,
    /// Source or reservoir contribution.
    Source,
    /// Irreversible production/loss contribution.
    Dissipation,
}

impl AccountingTermKind {
    const fn as_str(self) -> &'static str {
        match self {
            Self::Reversible => "reversible",
            Self::Storage => "storage",
            Self::Source => "source",
            Self::Dissipation => "dissipation",
        }
    }
}

/// Explicit ownership status for an accounting term.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum OwnershipDisposition {
    /// Only reversible terms may carry no ownership concept.
    NotApplicable,
    /// Exactly one stable component/operator owns the term.
    Owned(StableId),
    /// A source/loss is intentionally unowned under a retained rationale ID.
    ExplicitlyUnowned {
        /// Durable reason, policy, or scope-exclusion identifier.
        rationale: StableId,
    },
}

impl OwnershipDisposition {
    /// Borrow the unique owner when one exists.
    #[must_use]
    pub fn owner(&self) -> Option<&StableId> {
        match self {
            Self::Owned(owner) => Some(owner),
            Self::NotApplicable | Self::ExplicitlyUnowned { .. } => None,
        }
    }

    fn diagnostic(&self) -> String {
        match self {
            Self::NotApplicable => "not-applicable".to_string(),
            Self::Owned(owner) => format!("owned:{}", owner.as_str()),
            Self::ExplicitlyUnowned { rationale } => {
                format!("explicitly-unowned:{}", rationale.as_str())
            }
        }
    }
}

/// Discretization selected for a neutral port shape.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PortDiscretization {
    /// Scalar/vector/tensor port coordinates represented as one raw block.
    Lumped,
    /// Field-duality port with explicit nonzero effort/flow dof counts.
    Field {
        /// Total effort-coordinate dofs, including all components.
        effort_dofs: NonZeroUsize,
        /// Total flow-coordinate dofs, including all components.
        flow_dofs: NonZeroUsize,
    },
}

impl PortDiscretization {
    /// Canonical lumped port discretization.
    #[must_use]
    pub const fn lumped() -> Self {
        Self::Lumped
    }

    /// Construct a field discretization without admitting empty vectors.
    ///
    /// # Errors
    /// [`PortEquationError::ZeroFieldDofs`] names the empty side.
    pub fn field(effort_dofs: usize, flow_dofs: usize) -> Result<Self, PortEquationError> {
        let effort_dofs = NonZeroUsize::new(effort_dofs)
            .ok_or(PortEquationError::ZeroFieldDofs { variable: "effort" })?;
        let flow_dofs = NonZeroUsize::new(flow_dofs)
            .ok_or(PortEquationError::ZeroFieldDofs { variable: "flow" })?;
        Ok(Self::Field {
            effort_dofs,
            flow_dofs,
        })
    }
}

/// Which power-conjugate coordinate is produced by a storage gradient.
///
/// The neutral [`StorageElement`] identifies the state schema, coordinate
/// count, and constitutive-gradient operator, but it deliberately does not
/// guess whether a domain's co-energy variable is the port effort or flow.
/// The compiler caller must make that binding explicit.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StorageGradientTarget {
    /// The constitutive gradient produces the port effort coordinate.
    Effort,
    /// The constitutive gradient produces the port flow coordinate.
    Flow,
}

impl StorageGradientTarget {
    const fn as_str(self) -> &'static str {
        match self {
            Self::Effort => "effort",
            Self::Flow => "flow",
        }
    }
}

/// Compiler-side type declaration for one storage constitutive-gradient
/// action.
///
/// `state_coordinate_dims` supplies the physical dimensions omitted by the
/// neutral storage descriptor. The state coordinate count still comes from
/// [`StorageElement::state_dimension`], and the output space comes from the
/// selected admitted port side, so neither cardinality is caller-repeated.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct StorageStateAction {
    state_coordinate_dims: Dims,
    gradient_target: StorageGradientTarget,
}

impl StorageStateAction {
    /// Declare the state-coordinate dimensions and the port side produced by
    /// the constitutive gradient.
    #[must_use]
    pub const fn new(state_coordinate_dims: Dims, gradient_target: StorageGradientTarget) -> Self {
        Self {
            state_coordinate_dims,
            gradient_target,
        }
    }

    /// Physical dimensions shared by the opaque state coordinates.
    #[must_use]
    pub const fn state_coordinate_dims(self) -> Dims {
        self.state_coordinate_dims
    }

    /// Port coordinate produced by the external gradient action.
    #[must_use]
    pub const fn gradient_target(self) -> StorageGradientTarget {
        self.gradient_target
    }
}

/// Closed port primitive family retained by a compiled equation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PortPrimitiveKind {
    /// One side of an admitted lossless two-port junction.
    ConservativeJunction,
    /// One side of an evidence-bound reversible skew coupling.
    ReversibleSkewCoupling,
    /// Stored-energy descriptor with a named constitutive gradient.
    Storage,
    /// Evidence-bound irreversible constitutive relation.
    Dissipation,
    /// Included source or external reservoir boundary.
    SourceOrReservoir,
}

impl PortPrimitiveKind {
    const fn as_str(self) -> &'static str {
        match self {
            Self::ConservativeJunction => "conservative-junction",
            Self::ReversibleSkewCoupling => "reversible-skew-coupling",
            Self::Storage => "storage",
            Self::Dissipation => "dissipation",
            Self::SourceOrReservoir => "source-or-reservoir",
        }
    }
}

/// Exact side of a [`ReversibleSkewCoupling`].
///
/// The ordered roles are semantic: the forward action maps side A effort to
/// side B flow, while the adjoint action maps side B effort to side A flow and
/// receives the explicit negative skew sign.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReversibleSkewSide {
    /// Adjoint-action side; its generated power contribution is negated.
    A,
    /// Forward-action side; its generated power contribution keeps its sign.
    B,
}

impl ReversibleSkewSide {
    const fn as_str(self) -> &'static str {
        match self {
            Self::A => "a-adjoint-negative",
            Self::B => "b-forward-positive",
        }
    }
}

/// Compiler-facing descriptor for a reversible cross coupling.
///
/// The structural convention is
/// `flow_b = forward_operator(effort_a)` and
/// `flow_a = -adjoint_operator(effort_b)`. `skew_adjoint_evidence` is the
/// durable external evidence reference asserting that the two actions are
/// adjoints under the declared port pairings. This compiler identity-binds
/// those references; it does not execute or authenticate them.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ReversibleSkewCoupling {
    id: StableId,
    port_a: PortSchema,
    port_b: PortSchema,
    forward_operator: StableId,
    adjoint_operator: StableId,
    skew_adjoint_evidence: StableId,
}

impl ReversibleSkewCoupling {
    /// Construct an ordered reversible coupling descriptor.
    ///
    /// # Errors
    /// Refuses reused port identities and relation/port identity aliasing.
    pub fn new(
        id: StableId,
        port_a: PortSchema,
        port_b: PortSchema,
        forward_operator: StableId,
        adjoint_operator: StableId,
        skew_adjoint_evidence: StableId,
    ) -> Result<Self, PortEquationError> {
        if port_a.id() == port_b.id() {
            return Err(PortEquationError::DuplicatePortId {
                port: port_a.id().as_str().to_string(),
            });
        }
        if &id == port_a.id() || &id == port_b.id() {
            let port = if &id == port_a.id() {
                port_a.id()
            } else {
                port_b.id()
            };
            return Err(PortEquationError::PrimitivePortIdentityAlias {
                primitive: id.as_str().to_string(),
                port: port.as_str().to_string(),
            });
        }
        Ok(Self {
            id,
            port_a,
            port_b,
            forward_operator,
            adjoint_operator,
            skew_adjoint_evidence,
        })
    }

    /// Stable coupling identity.
    #[must_use]
    pub const fn id(&self) -> &StableId {
        &self.id
    }

    /// Ordered A/B power ports.
    #[must_use]
    pub const fn ports(&self) -> (&PortSchema, &PortSchema) {
        (&self.port_a, &self.port_b)
    }

    /// Action mapping side A effort to side B flow.
    #[must_use]
    pub const fn forward_operator(&self) -> &StableId {
        &self.forward_operator
    }

    /// Adjoint action mapping side B effort to side A flow before negation.
    #[must_use]
    pub const fn adjoint_operator(&self) -> &StableId {
        &self.adjoint_operator
    }

    /// External evidence for the declared adjoint relation.
    #[must_use]
    pub const fn skew_adjoint_evidence(&self) -> &StableId {
        &self.skew_adjoint_evidence
    }
}

/// Exact side of an admitted [`ConservativeJunction`].
///
/// The side is identity-bearing because the neutral junction's scalar seed
/// declares side A's flow as supplied and side B's flow as its negative.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ConservativeJunctionSide {
    /// The junction's declared A side; its contribution keeps its sign.
    A,
    /// The junction's declared B side; its contribution is explicitly negated.
    B,
}

impl ConservativeJunctionSide {
    const fn as_str(self) -> &'static str {
        match self {
            Self::A => "a",
            Self::B => "b",
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
enum PortPrimitiveBinding {
    ConservativeJunction {
        junction: ConservativeJunction,
        side: ConservativeJunctionSide,
    },
    ReversibleSkewCoupling {
        coupling: ReversibleSkewCoupling,
        side: ReversibleSkewSide,
    },
    Storage {
        primitive: StorageElement,
        action: StorageStateAction,
    },
    Dissipation(DissipativeRelation),
    SourceOrReservoir(SourceOrReservoir),
}

impl PortPrimitiveBinding {
    fn kind(&self) -> PortPrimitiveKind {
        match self {
            Self::ConservativeJunction { .. } => PortPrimitiveKind::ConservativeJunction,
            Self::ReversibleSkewCoupling { .. } => PortPrimitiveKind::ReversibleSkewCoupling,
            Self::Storage { .. } => PortPrimitiveKind::Storage,
            Self::Dissipation(_) => PortPrimitiveKind::Dissipation,
            Self::SourceOrReservoir(_) => PortPrimitiveKind::SourceOrReservoir,
        }
    }

    fn id(&self) -> &StableId {
        match self {
            Self::ConservativeJunction { junction, .. } => junction.id(),
            Self::ReversibleSkewCoupling { coupling, .. } => coupling.id(),
            Self::Storage { primitive, .. } => primitive.id(),
            Self::Dissipation(primitive) => primitive.id(),
            Self::SourceOrReservoir(primitive) => primitive.id(),
        }
    }

    fn junction_side(&self) -> Option<ConservativeJunctionSide> {
        match self {
            Self::ConservativeJunction { side, .. } => Some(*side),
            Self::ReversibleSkewCoupling { .. }
            | Self::Storage { .. }
            | Self::Dissipation(_)
            | Self::SourceOrReservoir(_) => None,
        }
    }

    fn reversible_skew_side(&self) -> Option<ReversibleSkewSide> {
        match self {
            Self::ReversibleSkewCoupling { side, .. } => Some(*side),
            Self::ConservativeJunction { .. }
            | Self::Storage { .. }
            | Self::Dissipation(_)
            | Self::SourceOrReservoir(_) => None,
        }
    }

    fn reversible_skew_coupling(&self) -> Option<&ReversibleSkewCoupling> {
        match self {
            Self::ReversibleSkewCoupling { coupling, .. } => Some(coupling),
            Self::ConservativeJunction { .. }
            | Self::Storage { .. }
            | Self::Dissipation(_)
            | Self::SourceOrReservoir(_) => None,
        }
    }

    fn storage(&self) -> Option<(&StorageElement, StorageStateAction)> {
        match self {
            Self::Storage { primitive, action } => Some((primitive, *action)),
            Self::ConservativeJunction { .. }
            | Self::ReversibleSkewCoupling { .. }
            | Self::Dissipation(_)
            | Self::SourceOrReservoir(_) => None,
        }
    }
}

/// Complete request to compile one admitted port schema.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PortEquationSpec {
    schema: PortSchema,
    discretization: PortDiscretization,
    sense: PortEquationSense,
    term_kind: AccountingTermKind,
    ownership: OwnershipDisposition,
    primitive: Option<PortPrimitiveBinding>,
}

impl PortEquationSpec {
    /// Bind one port schema to its discretization, algebraic sense, and exact
    /// accounting ownership declaration. Validation occurs transactionally in
    /// [`compile_port_equation`] or [`compile_port_equations`].
    #[must_use]
    pub fn new(
        schema: PortSchema,
        discretization: PortDiscretization,
        sense: PortEquationSense,
        term_kind: AccountingTermKind,
        ownership: OwnershipDisposition,
    ) -> Self {
        Self {
            schema,
            discretization,
            sense,
            term_kind,
            ownership,
            primitive: None,
        }
    }

    /// Construct the two signed power contributions of an admitted lossless
    /// junction.
    ///
    /// The returned requests retain the neutral descriptor's A/B order: side A
    /// is as declared and side B is explicitly reversed. They intentionally
    /// compile as two independent system fragments because the v1 system IR
    /// has no semantic interface-side field identity. The upstream junction's
    /// shared-effort/opposite-flow constraint is identity-bound but is not
    /// re-emitted or executed by this structural lowering.
    #[must_use]
    pub fn from_conservative_junction(
        primitive: ConservativeJunction,
        discretization: PortDiscretization,
    ) -> [Self; 2] {
        let (port_a, port_b) = primitive.ports();
        let schema_a = port_a.clone();
        let schema_b = port_b.clone();
        [
            Self {
                schema: schema_a,
                discretization,
                sense: PortEquationSense::AsDeclared,
                term_kind: AccountingTermKind::Reversible,
                ownership: OwnershipDisposition::NotApplicable,
                primitive: Some(PortPrimitiveBinding::ConservativeJunction {
                    junction: primitive.clone(),
                    side: ConservativeJunctionSide::A,
                }),
            },
            Self {
                schema: schema_b,
                discretization,
                sense: PortEquationSense::Reversed,
                term_kind: AccountingTermKind::Reversible,
                ownership: OwnershipDisposition::NotApplicable,
                primitive: Some(PortPrimitiveBinding::ConservativeJunction {
                    junction: primitive,
                    side: ConservativeJunctionSide::B,
                }),
            },
        ]
    }

    /// Construct the two signed power contributions of an evidence-bound
    /// reversible skew coupling.
    ///
    /// Side A carries the explicit negative adjoint-action contribution and
    /// side B carries the positive forward-action contribution. The operator
    /// references and skew-adjoint evidence are identity-bearing, but their
    /// numeric action is not executed by this structural lowering.
    #[must_use]
    pub fn from_reversible_skew_coupling(
        coupling: ReversibleSkewCoupling,
        discretization: PortDiscretization,
    ) -> [Self; 2] {
        let (port_a, port_b) = coupling.ports();
        let schema_a = port_a.clone();
        let schema_b = port_b.clone();
        [
            Self {
                schema: schema_a,
                discretization,
                sense: PortEquationSense::Reversed,
                term_kind: AccountingTermKind::Reversible,
                ownership: OwnershipDisposition::NotApplicable,
                primitive: Some(PortPrimitiveBinding::ReversibleSkewCoupling {
                    coupling: coupling.clone(),
                    side: ReversibleSkewSide::A,
                }),
            },
            Self {
                schema: schema_b,
                discretization,
                sense: PortEquationSense::AsDeclared,
                term_kind: AccountingTermKind::Reversible,
                ownership: OwnershipDisposition::NotApplicable,
                primitive: Some(PortPrimitiveBinding::ReversibleSkewCoupling {
                    coupling,
                    side: ReversibleSkewSide::B,
                }),
            },
        ]
    }

    /// Construct a stored-energy port equation whose owner, state schema,
    /// coordinate count, and constitutive-gradient reference come from one
    /// admitted primitive.
    ///
    /// The compiler caller supplies the physical dimensions of the opaque
    /// state coordinates and explicitly selects whether the gradient produces
    /// effort or flow. Compilation turns that declaration into an owned state
    /// field and a typed external atom application; it does not execute or
    /// authenticate the referenced operator.
    #[must_use]
    pub fn from_storage(
        primitive: StorageElement,
        action: StorageStateAction,
        discretization: PortDiscretization,
        sense: PortEquationSense,
    ) -> Self {
        Self {
            schema: primitive.port().clone(),
            discretization,
            sense,
            term_kind: AccountingTermKind::Storage,
            ownership: OwnershipDisposition::Owned(primitive.id().clone()),
            primitive: Some(PortPrimitiveBinding::Storage { primitive, action }),
        }
    }

    /// Construct an irreversible port equation whose constitutive operator and
    /// sign evidence are bound into semantic identity.
    #[must_use]
    pub fn from_dissipation(
        primitive: DissipativeRelation,
        discretization: PortDiscretization,
        sense: PortEquationSense,
    ) -> Self {
        Self {
            schema: primitive.port().clone(),
            discretization,
            sense,
            term_kind: AccountingTermKind::Dissipation,
            ownership: OwnershipDisposition::Owned(primitive.id().clone()),
            primitive: Some(PortPrimitiveBinding::Dissipation(primitive)),
        }
    }

    /// Construct a source/reservoir port equation whose class and explicit
    /// signed accounting boundary are bound into semantic identity.
    #[must_use]
    pub fn from_source_or_reservoir(
        primitive: SourceOrReservoir,
        discretization: PortDiscretization,
        sense: PortEquationSense,
    ) -> Self {
        Self {
            schema: primitive.port().clone(),
            discretization,
            sense,
            term_kind: AccountingTermKind::Source,
            ownership: OwnershipDisposition::Owned(primitive.id().clone()),
            primitive: Some(PortPrimitiveBinding::SourceOrReservoir(primitive)),
        }
    }

    /// Borrow the neutral source schema.
    #[must_use]
    pub const fn schema(&self) -> &PortSchema {
        &self.schema
    }

    /// Retained admitted primitive family, when constructed from a closed
    /// neutral primitive or compiler-facing reversible-skew descriptor.
    #[must_use]
    pub fn primitive_kind(&self) -> Option<PortPrimitiveKind> {
        self.primitive.as_ref().map(PortPrimitiveBinding::kind)
    }

    /// Stable relation identity of the retained primitive, if present.
    #[must_use]
    pub fn primitive_id(&self) -> Option<&StableId> {
        self.primitive.as_ref().map(PortPrimitiveBinding::id)
    }

    /// Retained A/B role when this request comes from a conservative junction.
    #[must_use]
    pub fn conservative_junction_side(&self) -> Option<ConservativeJunctionSide> {
        self.primitive
            .as_ref()
            .and_then(PortPrimitiveBinding::junction_side)
    }

    /// Retained A/B role when this request comes from a reversible skew
    /// coupling.
    #[must_use]
    pub fn reversible_skew_side(&self) -> Option<ReversibleSkewSide> {
        self.primitive
            .as_ref()
            .and_then(PortPrimitiveBinding::reversible_skew_side)
    }

    fn reversible_skew_coupling(&self) -> Option<&ReversibleSkewCoupling> {
        self.primitive
            .as_ref()
            .and_then(PortPrimitiveBinding::reversible_skew_coupling)
    }

    fn storage(&self) -> Option<(&StorageElement, StorageStateAction)> {
        self.primitive
            .as_ref()
            .and_then(PortPrimitiveBinding::storage)
    }
}

/// Exactly-once ownership of the energy coordinate carried by a stream.
///
/// Unlike [`OwnershipDisposition`], this type has no `NotApplicable` state:
/// every admitted [`StreamPort`] carries energy, so its compiler request must
/// either name one owner or retain an explicit durable reason why ownership is
/// outside the compiled model.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum StreamEnergyOwnership {
    /// Exactly one stable component/operator owns the stream-energy term.
    Owned(StableId),
    /// Ownership is outside the compiled model under a retained rationale ID.
    ExplicitlyUnowned {
        /// Durable reason, policy, or scope-exclusion identifier.
        rationale: StableId,
    },
}

impl StreamEnergyOwnership {
    /// Borrow the unique owner when one exists.
    #[must_use]
    pub fn owner(&self) -> Option<&StableId> {
        match self {
            Self::Owned(owner) => Some(owner),
            Self::ExplicitlyUnowned { .. } => None,
        }
    }

    fn diagnostic(&self) -> String {
        match self {
            Self::Owned(owner) => format!("owned:{}", owner.as_str()),
            Self::ExplicitlyUnowned { rationale } => {
                format!("explicitly-unowned:{}", rationale.as_str())
            }
        }
    }
}

/// Publicly observable selected chart family of one admitted stream.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StreamEnergyChartKind {
    /// Canonical moving-stream enthalpy chart.
    MovingStreamEnthalpy,
    /// Internal energy with exact pressure/deviatoric Cauchy-work crosswalk.
    InternalEnergyCauchyWork,
    /// Exact Euler/Legendre conjugate-potential chart.
    ConjugatePotential,
}

impl StreamEnergyChartKind {
    fn from_chart(chart: &StreamEnergyChart) -> Self {
        match chart {
            StreamEnergyChart::MovingStreamEnthalpy(_) => Self::MovingStreamEnthalpy,
            StreamEnergyChart::InternalEnergyCauchyWork(_) => Self::InternalEnergyCauchyWork,
            StreamEnergyChart::ConjugatePotential(_) => Self::ConjugatePotential,
        }
    }

    const fn as_str(self) -> &'static str {
        match self {
            Self::MovingStreamEnthalpy => "moving-stream-enthalpy",
            Self::InternalEnergyCauchyWork => "internal-energy-cauchy-work",
            Self::ConjugatePotential => "conjugate-potential",
        }
    }
}

/// Complete request to lower one already-admitted neutral stream bundle.
///
/// `source_receipt` is the caller's durable identity for the upstream chart
/// admission. The compiler additionally binds every public context/rate field
/// into its own [`SystemId`], but does not claim to re-run or authenticate the
/// upstream thermodynamic evidence.
#[derive(Debug, Clone, PartialEq)]
pub struct StreamEquationSpec {
    stream: StreamPort,
    source_receipt: StableId,
    energy_ownership: StreamEnergyOwnership,
}

impl StreamEquationSpec {
    /// Bind an admitted stream to its upstream receipt and exactly-once energy
    /// ownership declaration.
    #[must_use]
    pub fn new(
        stream: StreamPort,
        source_receipt: StableId,
        energy_ownership: StreamEnergyOwnership,
    ) -> Self {
        Self {
            stream,
            source_receipt,
            energy_ownership,
        }
    }

    /// Borrow the admitted neutral stream bundle.
    #[must_use]
    pub const fn stream(&self) -> &StreamPort {
        &self.stream
    }

    /// Borrow the upstream stream-admission receipt identity.
    #[must_use]
    pub const fn source_receipt(&self) -> &StableId {
        &self.source_receipt
    }

    /// Borrow the exactly-once stream-energy ownership declaration.
    #[must_use]
    pub const fn energy_ownership(&self) -> &StreamEnergyOwnership {
        &self.energy_ownership
    }
}

/// One declaration admitted by the cross-kind interface compiler.
#[derive(Debug, Clone, PartialEq)]
pub enum InterfaceEquationSpec {
    /// Power-conjugate effort/flow declaration.
    Port(PortEquationSpec),
    /// Multi-conserved-flux stream declaration.
    Stream(StreamEquationSpec),
}

impl InterfaceEquationSpec {
    fn port_id(&self) -> &StableId {
        match self {
            Self::Port(spec) => spec.schema.id(),
            Self::Stream(spec) => spec.stream.binding().port_id(),
        }
    }

    fn owner(&self) -> Option<&StableId> {
        match self {
            Self::Port(spec) => spec.ownership.owner(),
            Self::Stream(spec) => spec.energy_ownership.owner(),
        }
    }
}

/// Structured refusal from port-equation lowering.
#[derive(Debug, PartialEq)]
pub enum PortEquationError {
    /// A deterministic batch was empty.
    EmptyBatch,
    /// The request exceeded the static equation ceiling.
    TooManyEquations {
        /// Requested equations.
        count: usize,
        /// Static ceiling.
        cap: usize,
    },
    /// Two source declarations used the same port identity.
    DuplicatePortId {
        /// Duplicated stable port identity.
        port: String,
    },
    /// A closed primitive reused one of its port identities.
    PrimitivePortIdentityAlias {
        /// Stable primitive identity.
        primitive: String,
        /// Aliased stable port identity.
        port: String,
    },
    /// One stable interface identity requested both an effort/flow energy term
    /// and a stream-carried energy term.
    DuplicateEnergyCarrier {
        /// Stable boundary identity carrying energy twice.
        port: String,
    },
    /// One owner was assigned to two independently generated terms.
    DuplicateOwnership {
        /// Duplicated owner identity.
        owner: String,
        /// First port using the owner.
        first_port: String,
        /// Second port using the owner.
        second_port: String,
    },
    /// A field discretization had zero dofs on one side.
    ZeroFieldDofs {
        /// `effort` or `flow`.
        variable: &'static str,
    },
    /// Lumped and field port shapes used incompatible discretizations.
    DiscretizationMismatch {
        /// Shape expected by the schema.
        expected: &'static str,
        /// Discretization supplied by the caller.
        actual: &'static str,
    },
    /// Field dofs did not contain a whole number of component tuples.
    FieldComponentMismatch {
        /// `effort` or `flow`.
        variable: &'static str,
        /// Total dofs.
        dofs: usize,
        /// Components per field point.
        components: usize,
    },
    /// Tensor component count overflowed the platform index type.
    ShapeExtentOverflow,
    /// The upstream schema version is not the version this compiler binds.
    PortSchemaVersionMismatch {
        /// Supported schema version.
        expected: u16,
        /// Received schema version.
        actual: u16,
    },
    /// The upstream stream version is not the version this compiler binds.
    StreamSchemaVersionMismatch {
        /// Supported stream version.
        expected: u16,
        /// Received stream version.
        actual: u16,
    },
    /// A future/upstream schema no longer re-derived exact power dimensions.
    SchemaPowerDrift {
        /// Effort dimensions.
        effort: Dims,
        /// Flow dimensions.
        flow: Dims,
        /// Integration-measure dimensions.
        measure: Dims,
    },
    /// The ownership declaration contradicted the accounting term kind.
    OwnershipMismatch {
        /// Declared accounting role.
        term_kind: AccountingTermKind,
        /// Stable diagnostic for the supplied disposition.
        disposition: String,
    },
    /// Identity-bearing compiler metadata exceeded the system extension cap.
    CompilerMetadataTooLarge {
        /// Estimated bytes.
        bytes: usize,
        /// Static cap.
        cap: usize,
    },
    /// A bounded metadata allocation was refused.
    Resource,
    /// The underlying system IR refused the generated structure.
    System(Box<SystemTypeError>),
}

impl core::fmt::Display for PortEquationError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::EmptyBatch => f.write_str("port-equation batch is empty"),
            Self::TooManyEquations { count, cap } => {
                write!(f, "port-equation count {count} exceeds cap {cap}")
            }
            Self::DuplicatePortId { port } => {
                write!(f, "port identity {port} appears more than once")
            }
            Self::PrimitivePortIdentityAlias { primitive, port } => write!(
                f,
                "primitive identity {primitive} aliases its port identity {port}"
            ),
            Self::DuplicateEnergyCarrier { port } => write!(
                f,
                "interface {port} carries energy through both an effort/flow term and a stream bundle"
            ),
            Self::DuplicateOwnership {
                owner,
                first_port,
                second_port,
            } => write!(
                f,
                "accounting owner {owner} is assigned to both ports {first_port} and {second_port}"
            ),
            Self::ZeroFieldDofs { variable } => {
                write!(f, "field port has zero {variable} dofs")
            }
            Self::DiscretizationMismatch { expected, actual } => write!(
                f,
                "port shape requires {expected} discretization, received {actual}"
            ),
            Self::FieldComponentMismatch {
                variable,
                dofs,
                components,
            } => write!(
                f,
                "field port {variable} dof count {dofs} is not divisible by component count {components}"
            ),
            Self::ShapeExtentOverflow => {
                f.write_str("tensor port component count overflowed usize")
            }
            Self::PortSchemaVersionMismatch { expected, actual } => write!(
                f,
                "port schema version {actual} is unsupported; expected {expected}"
            ),
            Self::StreamSchemaVersionMismatch { expected, actual } => write!(
                f,
                "stream schema version {actual} is unsupported; expected {expected}"
            ),
            Self::SchemaPowerDrift {
                effort,
                flow,
                measure,
            } => write!(
                f,
                "port schema no longer re-derives power dimensions: effort {effort:?}, flow {flow:?}, measure {measure:?}"
            ),
            Self::OwnershipMismatch {
                term_kind,
                disposition,
            } => write!(
                f,
                "{} term has inadmissible ownership disposition {disposition}",
                term_kind.as_str()
            ),
            Self::CompilerMetadataTooLarge { bytes, cap } => write!(
                f,
                "port compiler metadata estimate {bytes} bytes exceeds cap {cap}"
            ),
            Self::Resource => f.write_str("port compiler metadata allocation was refused"),
            Self::System(error) => write!(f, "generated system IR refused: {error}"),
        }
    }
}

impl std::error::Error for PortEquationError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::System(error) => Some(error.as_ref()),
            _ => None,
        }
    }
}

impl From<SystemTypeError> for PortEquationError {
    fn from(value: SystemTypeError) -> Self {
        Self::System(Box::new(value))
    }
}

/// Structural proof receipt for one generated port power equation.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PortEquationReceipt {
    port_id: String,
    port_schema_version: u16,
    system_identity: SystemId,
    kind: PortKind,
    pairing: PowerPairing,
    effort_dims: Dims,
    flow_dims: Dims,
    measure_dims: Dims,
    product_dims: Dims,
    sense: PortEquationSense,
    term_kind: AccountingTermKind,
    ownership: OwnershipDisposition,
    loss_ownership_id: Option<LossOwnershipId>,
    primitive_kind: Option<PortPrimitiveKind>,
    primitive_id: Option<String>,
    conservative_junction_side: Option<ConservativeJunctionSide>,
    reversible_skew_side: Option<ReversibleSkewSide>,
    reversible_skew_forward_operator: Option<String>,
    reversible_skew_adjoint_operator: Option<String>,
    reversible_skew_evidence: Option<String>,
    storage_state_schema: Option<String>,
    storage_state_dimension: Option<usize>,
    storage_state_action: Option<StorageStateAction>,
    storage_gradient_operator: Option<String>,
}

impl PortEquationReceipt {
    /// Stable source port identity.
    #[must_use]
    pub fn port_id(&self) -> &str {
        &self.port_id
    }

    /// Neutral source schema version bound by this compiler.
    #[must_use]
    pub const fn port_schema_version(&self) -> u16 {
        self.port_schema_version
    }

    /// Identity of the generated, fully admitted system fragment.
    #[must_use]
    pub const fn system_identity(&self) -> SystemId {
        self.system_identity
    }

    /// Physical effort/flow port kind.
    #[must_use]
    pub const fn kind(&self) -> PortKind {
        self.kind
    }

    /// Exact contraction declared by the neutral schema.
    #[must_use]
    pub const fn pairing(&self) -> PowerPairing {
        self.pairing
    }

    /// Effort-coordinate dimensions before measure application.
    #[must_use]
    pub const fn effort_dims(&self) -> Dims {
        self.effort_dims
    }

    /// Flow-coordinate dimensions before measure application.
    #[must_use]
    pub const fn flow_dims(&self) -> Dims {
        self.flow_dims
    }

    /// Integration-measure dimensions used by the contraction.
    #[must_use]
    pub const fn measure_dims(&self) -> Dims {
        self.measure_dims
    }

    /// Explicit declared/reversed algebraic sense.
    #[must_use]
    pub const fn sense(&self) -> PortEquationSense {
        self.sense
    }

    /// Exact multiplier inserted for orientation handling.
    #[must_use]
    pub const fn sign(&self) -> i8 {
        self.sense.sign()
    }

    /// Declared accounting role.
    #[must_use]
    pub const fn term_kind(&self) -> AccountingTermKind {
        self.term_kind
    }

    /// Exact ownership disposition retained by the compiler.
    #[must_use]
    pub const fn ownership(&self) -> &OwnershipDisposition {
        &self.ownership
    }

    /// Compiler-derived nominal identity for a concretely owned dissipative
    /// term. Explicitly unowned losses and non-loss roles return `None`.
    #[must_use]
    pub const fn loss_ownership_id(&self) -> Option<LossOwnershipId> {
        self.loss_ownership_id
    }

    /// Retained closed primitive family, if this request was constructed from
    /// an admitted neutral primitive or reversible-skew descriptor.
    #[must_use]
    pub const fn primitive_kind(&self) -> Option<PortPrimitiveKind> {
        self.primitive_kind
    }

    /// Stable relation identity of the retained primitive.
    #[must_use]
    pub fn primitive_id(&self) -> Option<&str> {
        self.primitive_id.as_deref()
    }

    /// Exact A/B role when the equation is one side of a conservative
    /// junction.
    #[must_use]
    pub const fn conservative_junction_side(&self) -> Option<ConservativeJunctionSide> {
        self.conservative_junction_side
    }

    /// Exact signed A/B action role for a reversible skew coupling.
    #[must_use]
    pub const fn reversible_skew_side(&self) -> Option<ReversibleSkewSide> {
        self.reversible_skew_side
    }

    /// Forward action mapping side A effort to side B flow.
    #[must_use]
    pub fn reversible_skew_forward_operator(&self) -> Option<&str> {
        self.reversible_skew_forward_operator.as_deref()
    }

    /// Adjoint action mapping side B effort to side A flow before negation.
    #[must_use]
    pub fn reversible_skew_adjoint_operator(&self) -> Option<&str> {
        self.reversible_skew_adjoint_operator.as_deref()
    }

    /// External evidence reference for the declared skew-adjoint relation.
    #[must_use]
    pub fn reversible_skew_evidence(&self) -> Option<&str> {
        self.reversible_skew_evidence.as_deref()
    }

    /// State-schema identity retained by a closed storage action.
    #[must_use]
    pub fn storage_state_schema(&self) -> Option<&str> {
        self.storage_state_schema.as_deref()
    }

    /// Exact upstream state-coordinate count consumed by the storage
    /// gradient.
    #[must_use]
    pub const fn storage_state_dimension(&self) -> Option<usize> {
        self.storage_state_dimension
    }

    /// Compiler-side physical state dimensions and explicit effort/flow
    /// gradient target.
    #[must_use]
    pub const fn storage_state_action(&self) -> Option<StorageStateAction> {
        self.storage_state_action
    }

    /// External constitutive-gradient operator bound into the generated atom.
    #[must_use]
    pub fn storage_gradient_operator(&self) -> Option<&str> {
        self.storage_gradient_operator.as_deref()
    }

    /// Re-derived power dimensions after the pairing measure is applied.
    #[must_use]
    pub const fn product_dims(&self) -> Dims {
        self.product_dims
    }

    /// Deterministic diagnostic transport. This is a receipt view, not a
    /// substitute for the typed [`SystemId`].
    #[must_use]
    pub fn to_json(&self) -> String {
        let loss_ownership_id = self
            .loss_ownership_id
            .map_or_else(|| "null".to_string(), |id| format!("\"{}\"", id.to_hex()));
        let reversible_skew_fields = match (
            self.reversible_skew_side,
            self.reversible_skew_forward_operator.as_ref(),
            self.reversible_skew_adjoint_operator.as_ref(),
            self.reversible_skew_evidence.as_ref(),
        ) {
            (Some(side), Some(forward), Some(adjoint), Some(evidence)) => format!(
                ",\"reversible_skew_side\":\"{}\",\
                 \"reversible_skew_forward_operator\":\"{forward}\",\
                 \"reversible_skew_adjoint_operator\":\"{adjoint}\",\
                 \"reversible_skew_evidence\":\"{evidence}\"",
                side.as_str()
            ),
            _ => String::new(),
        };
        let storage_fields = match (
            self.storage_state_schema.as_ref(),
            self.storage_state_dimension,
            self.storage_state_action,
            self.storage_gradient_operator.as_ref(),
        ) {
            (Some(schema), Some(dimension), Some(action), Some(operator)) => format!(
                ",\"storage_state_schema\":\"{schema}\",\
                 \"storage_state_dimension\":{dimension},\
                 \"storage_state_dims\":{},\
                 \"storage_gradient_target\":\"{}\",\
                 \"storage_gradient_operator\":\"{operator}\"",
                dims_json(action.state_coordinate_dims()),
                action.gradient_target().as_str(),
            ),
            _ => String::new(),
        };
        let primitive_fields = self
            .primitive_kind
            .zip(self.primitive_id.as_ref())
            .map_or_else(String::new, |(kind, id)| {
                let side = self
                    .conservative_junction_side
                    .map_or_else(String::new, |side| {
                        format!(",\"conservative_junction_side\":\"{}\"", side.as_str())
                    });
                format!(
                    ",\"primitive_kind\":\"{}\",\"primitive_id\":\"{id}\"{side}{reversible_skew_fields}{storage_fields}",
                    kind.as_str()
                )
            });
        format!(
            "{{\"schema\":\"{}\",\"port_id\":\"{}\",\
             \"compiler_version\":\"{}\",\"feature\":\"port-equations\",\
             \"port_schema_version\":{},\"system_id\":\"{}\",\
             \"port_kind\":\"{}\",\"pairing\":\"{}\",\
             \"effort_dims\":{},\"flow_dims\":{},\"measure_dims\":{},\
             \"product_dims\":{},\"sense\":\"{}\",\"sign\":{},\
             \"term_kind\":\"{}\",\"ownership\":\"{}\",\
             \"loss_ownership_id\":{}{},\
             \"authority\":\"structural-generated\",\
             \"no_claim\":\"numeric contraction, quadrature, adapter truth, referenced operator/evidence truth, and closed-window conservation remain external\"}}",
            PORT_EQUATION_RECEIPT_SCHEMA_V1,
            self.port_id,
            crate::VERSION,
            self.port_schema_version,
            self.system_identity,
            port_kind_name(self.kind),
            pairing_name(self.pairing),
            dims_json(self.effort_dims),
            dims_json(self.flow_dims),
            dims_json(self.measure_dims),
            dims_json(self.product_dims),
            self.sense.as_str(),
            self.sign(),
            self.term_kind.as_str(),
            self.ownership.diagnostic(),
            loss_ownership_id,
            primitive_fields,
        )
    }
}

/// Structural receipt for one admitted stream-bundle lowering.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StreamEquationReceipt {
    port_id: String,
    stream_schema_version: u16,
    source_receipt: String,
    system_identity: SystemId,
    chart_kind: StreamEnergyChartKind,
    constituent_count: usize,
    energy_flow_bits: u64,
    energy_ownership: StreamEnergyOwnership,
}

impl StreamEquationReceipt {
    /// Stable source stream identity.
    #[must_use]
    pub fn port_id(&self) -> &str {
        &self.port_id
    }

    /// Neutral source stream version bound by this compiler.
    #[must_use]
    pub const fn stream_schema_version(&self) -> u16 {
        self.stream_schema_version
    }

    /// Upstream chart-admission receipt supplied by the caller.
    #[must_use]
    pub fn source_receipt(&self) -> &str {
        &self.source_receipt
    }

    /// Identity of the generated, fully admitted multi-field fragment.
    #[must_use]
    pub const fn system_identity(&self) -> SystemId {
        self.system_identity
    }

    /// Selected upstream stream-energy chart family.
    #[must_use]
    pub const fn chart_kind(&self) -> StreamEnergyChartKind {
        self.chart_kind
    }

    /// Number of canonical constituent amount-flow coordinates.
    #[must_use]
    pub const fn constituent_count(&self) -> usize {
        self.constituent_count
    }

    /// Exact IEEE-754 bits of the admitted signed stream energy rate.
    #[must_use]
    pub const fn energy_flow_bits(&self) -> u64 {
        self.energy_flow_bits
    }

    /// Exactly-once stream-energy ownership declaration.
    #[must_use]
    pub const fn energy_ownership(&self) -> &StreamEnergyOwnership {
        &self.energy_ownership
    }

    /// Deterministic diagnostic transport. The upstream receipt remains an
    /// asserted provenance reference, not authenticated evidence.
    #[must_use]
    pub fn to_json(&self) -> String {
        format!(
            "{{\"schema\":\"{}\",\"port_id\":\"{}\",\
             \"compiler_version\":\"{}\",\"feature\":\"port-equations\",\
             \"stream_schema_version\":{},\"source_receipt\":\"{}\",\
             \"system_id\":\"{}\",\"chart_kind\":\"{}\",\
             \"constituent_count\":{},\"energy_flow_bits\":\"{:016x}\",\
             \"energy_ownership\":\"{}\",\
             \"authority\":\"structural-generated\",\
             \"no_claim\":\"upstream chart evidence, numeric boundary application, and closed-window conservation remain external\"}}",
            STREAM_EQUATION_RECEIPT_SCHEMA_V1,
            self.port_id,
            crate::VERSION,
            self.stream_schema_version,
            self.source_receipt,
            self.system_identity,
            self.chart_kind.as_str(),
            self.constituent_count,
            self.energy_flow_bits,
            self.energy_ownership.diagnostic(),
        )
    }
}

/// One generated and fully admitted system fragment plus its receipt.
#[derive(Debug)]
pub struct CompiledPortEquation {
    system: AdmittedSystem,
    receipt: PortEquationReceipt,
}

impl CompiledPortEquation {
    /// Borrow the admitted system fragment.
    #[must_use]
    pub const fn system(&self) -> &AdmittedSystem {
        &self.system
    }

    /// Borrow the structural lowering receipt.
    #[must_use]
    pub const fn receipt(&self) -> &PortEquationReceipt {
        &self.receipt
    }

    /// Consume the wrapper and return the admitted system fragment.
    #[must_use]
    pub fn into_system(self) -> AdmittedSystem {
        self.system
    }
}

/// One generated stream field bundle plus its structural receipt.
#[derive(Debug)]
pub struct CompiledStreamEquation {
    system: AdmittedSystem,
    receipt: StreamEquationReceipt,
}

impl CompiledStreamEquation {
    /// Borrow the admitted multi-field fragment.
    #[must_use]
    pub const fn system(&self) -> &AdmittedSystem {
        &self.system
    }

    /// Borrow the structural lowering receipt.
    #[must_use]
    pub const fn receipt(&self) -> &StreamEquationReceipt {
        &self.receipt
    }

    /// Consume the wrapper and return the admitted system fragment.
    #[must_use]
    pub fn into_system(self) -> AdmittedSystem {
        self.system
    }
}

/// One cross-kind interface compilation result.
#[derive(Debug)]
pub enum CompiledInterfaceEquation {
    /// Generated effort/flow power equation.
    Port(CompiledPortEquation),
    /// Generated external stream field bundle.
    Stream(CompiledStreamEquation),
}

impl CompiledInterfaceEquation {
    /// Stable source interface identity.
    #[must_use]
    pub fn port_id(&self) -> &str {
        match self {
            Self::Port(compiled) => compiled.receipt.port_id(),
            Self::Stream(compiled) => compiled.receipt.port_id(),
        }
    }
}

/// Canonically port-ID-ordered result of a transactional compile batch.
#[derive(Debug)]
pub struct PortEquationBatch {
    equations: Vec<CompiledPortEquation>,
}

/// Canonically port-ID-ordered result of a transactional cross-kind batch.
#[derive(Debug)]
pub struct InterfaceEquationBatch {
    equations: Vec<CompiledInterfaceEquation>,
}

impl InterfaceEquationBatch {
    /// Generated declarations in canonical stable-port order.
    #[must_use]
    pub fn equations(&self) -> &[CompiledInterfaceEquation] {
        &self.equations
    }

    /// Consume the batch.
    #[must_use]
    pub fn into_equations(self) -> Vec<CompiledInterfaceEquation> {
        self.equations
    }
}

impl PortEquationBatch {
    /// Generated equations in canonical source-port order.
    #[must_use]
    pub fn equations(&self) -> &[CompiledPortEquation] {
        &self.equations
    }

    /// Consume the batch.
    #[must_use]
    pub fn into_equations(self) -> Vec<CompiledPortEquation> {
        self.equations
    }
}

/// Compile one neutral port declaration transactionally.
///
/// # Errors
/// Any schema, discretization, ownership, metadata, or system-IR refusal.
pub fn compile_port_equation(
    spec: PortEquationSpec,
) -> Result<CompiledPortEquation, PortEquationError> {
    validate_ownership(spec.term_kind, &spec.ownership)?;
    compile_one(spec)
}

/// Lower one already-admitted neutral stream into five typed external field
/// blocks: mass, constituent amount, momentum, energy, and entropy flow.
///
/// Use [`compile_interface_equations`] when effort/flow declarations may share
/// the same accounting boundary; that cross-kind entry point additionally
/// enforces stream-energy exclusivity.
///
/// # Errors
/// Any stream version, metadata, resource, or system-IR refusal.
pub fn compile_stream_equation(
    spec: StreamEquationSpec,
) -> Result<CompiledStreamEquation, PortEquationError> {
    compile_stream_one(spec)
}

/// Compile a deterministic mixed batch of effort/flow and stream declarations.
///
/// Input order is irrelevant. Reusing one stable port identity for both kinds
/// refuses before generation because the stream's admitted energy coordinate
/// would duplicate the effort/flow power term. Duplicate same-kind identities
/// and duplicate concrete accounting owners also refuse transactionally.
///
/// # Errors
/// Any batch, exclusivity, ownership, schema, metadata, resource, or system
/// refusal.
pub fn compile_interface_equations(
    mut specs: Vec<InterfaceEquationSpec>,
) -> Result<InterfaceEquationBatch, PortEquationError> {
    if specs.is_empty() {
        return Err(PortEquationError::EmptyBatch);
    }
    if specs.len() > MAX_PORT_EQUATIONS {
        return Err(PortEquationError::TooManyEquations {
            count: specs.len(),
            cap: MAX_PORT_EQUATIONS,
        });
    }
    specs.sort_by(|left, right| left.port_id().cmp(right.port_id()));
    for pair in specs.windows(2) {
        if pair[0].port_id() == pair[1].port_id() {
            let port = pair[0].port_id().as_str().to_string();
            return match (&pair[0], &pair[1]) {
                (InterfaceEquationSpec::Port(_), InterfaceEquationSpec::Stream(_))
                | (InterfaceEquationSpec::Stream(_), InterfaceEquationSpec::Port(_)) => {
                    Err(PortEquationError::DuplicateEnergyCarrier { port })
                }
                _ => Err(PortEquationError::DuplicatePortId { port }),
            };
        }
    }

    let mut owner_ports: BTreeMap<String, String> = BTreeMap::new();
    for spec in &specs {
        match spec {
            InterfaceEquationSpec::Port(spec) => {
                validate_ownership(spec.term_kind, &spec.ownership)?;
                extension_estimate(spec)?;
            }
            InterfaceEquationSpec::Stream(spec) => {
                validate_stream_version(spec)?;
                stream_extension_estimate(spec)?;
            }
        }
        if let Some(owner) = spec.owner()
            && let Some(first_port) = owner_ports.insert(
                owner.as_str().to_string(),
                spec.port_id().as_str().to_string(),
            )
        {
            return Err(PortEquationError::DuplicateOwnership {
                owner: owner.as_str().to_string(),
                first_port,
                second_port: spec.port_id().as_str().to_string(),
            });
        }
    }

    let mut equations = Vec::new();
    equations
        .try_reserve_exact(specs.len())
        .map_err(|_| PortEquationError::Resource)?;
    for spec in specs {
        equations.push(match spec {
            InterfaceEquationSpec::Port(spec) => {
                CompiledInterfaceEquation::Port(compile_one(spec)?)
            }
            InterfaceEquationSpec::Stream(spec) => {
                CompiledInterfaceEquation::Stream(compile_stream_one(spec)?)
            }
        });
    }
    Ok(InterfaceEquationBatch { equations })
}

/// Compile a deterministic batch. Input order is irrelevant; duplicate port
/// identities and duplicate concrete owners refuse before any output is
/// returned.
///
/// # Errors
/// Any batch, ownership, schema, discretization, metadata, or system refusal.
pub fn compile_port_equations(
    mut specs: Vec<PortEquationSpec>,
) -> Result<PortEquationBatch, PortEquationError> {
    if specs.is_empty() {
        return Err(PortEquationError::EmptyBatch);
    }
    if specs.len() > MAX_PORT_EQUATIONS {
        return Err(PortEquationError::TooManyEquations {
            count: specs.len(),
            cap: MAX_PORT_EQUATIONS,
        });
    }
    specs.sort_by(|left, right| left.schema.id().cmp(right.schema.id()));
    for pair in specs.windows(2) {
        if pair[0].schema.id() == pair[1].schema.id() {
            return Err(PortEquationError::DuplicatePortId {
                port: pair[0].schema.id().as_str().to_string(),
            });
        }
    }

    let mut owner_ports: BTreeMap<String, String> = BTreeMap::new();
    for spec in &specs {
        validate_ownership(spec.term_kind, &spec.ownership)?;
        extension_estimate(spec)?;
        if let Some(owner) = spec.ownership.owner()
            && let Some(first_port) = owner_ports.insert(
                owner.as_str().to_string(),
                spec.schema.id().as_str().to_string(),
            )
        {
            return Err(PortEquationError::DuplicateOwnership {
                owner: owner.as_str().to_string(),
                first_port,
                second_port: spec.schema.id().as_str().to_string(),
            });
        }
    }

    let mut equations = Vec::new();
    equations
        .try_reserve_exact(specs.len())
        .map_err(|_| PortEquationError::Resource)?;
    for spec in specs {
        equations.push(compile_one(spec)?);
    }
    Ok(PortEquationBatch { equations })
}

fn validate_ownership(
    term_kind: AccountingTermKind,
    ownership: &OwnershipDisposition,
) -> Result<(), PortEquationError> {
    let admitted = match term_kind {
        AccountingTermKind::Reversible => {
            matches!(ownership, OwnershipDisposition::NotApplicable)
        }
        AccountingTermKind::Storage => matches!(ownership, OwnershipDisposition::Owned(_)),
        AccountingTermKind::Source | AccountingTermKind::Dissipation => matches!(
            ownership,
            OwnershipDisposition::Owned(_) | OwnershipDisposition::ExplicitlyUnowned { .. }
        ),
    };
    if admitted {
        Ok(())
    } else {
        Err(PortEquationError::OwnershipMismatch {
            term_kind,
            disposition: ownership.diagnostic(),
        })
    }
}

fn compile_one(spec: PortEquationSpec) -> Result<CompiledPortEquation, PortEquationError> {
    if spec.schema.version() != PORT_SCHEMA_VERSION {
        return Err(PortEquationError::PortSchemaVersionMismatch {
            expected: PORT_SCHEMA_VERSION,
            actual: spec.schema.version(),
        });
    }
    let (effort_space, flow_space) = resolve_spaces(&spec.schema, spec.discretization)?;
    let measure_dims = pairing_measure(spec.schema.power_pairing());
    let product_dims = spec
        .schema
        .effort_dimensions()
        .checked_plus(spec.schema.flow_dimensions())
        .and_then(|sum| sum.checked_plus(measure_dims))
        .ok_or(PortEquationError::SchemaPowerDrift {
            effort: spec.schema.effort_dimensions(),
            flow: spec.schema.flow_dimensions(),
            measure: measure_dims,
        })?;
    if product_dims != POWER_DIMS {
        return Err(PortEquationError::SchemaPowerDrift {
            effort: spec.schema.effort_dimensions(),
            flow: spec.schema.flow_dimensions(),
            measure: measure_dims,
        });
    }

    let extension = encode_extension(&spec, effort_space, flow_space)?;
    let loss_ownership_id = derive_loss_ownership_id(&spec, effort_space, flow_space)?;
    let coordinates = coordinate_convention(&spec.schema)?;
    let clock = ConventionRef::new(spec.schema.timestamp().clock().as_str().to_string())?;
    let mut system = SystemDef::new().with_extension(extension)?;
    let effort = system.declare_field(FieldDecl {
        name: "port-effort".to_string(),
        space: effort_space,
        quantity: FieldQuantity::Dimensional(effort_space.dims),
        coordinates: coordinates.clone(),
        clock: clock.clone(),
        support: SpatialSupport::BoundaryTrace,
        state: StateOwnership::External,
    })?;
    let flow = system.declare_field(FieldDecl {
        name: "port-flow".to_string(),
        space: flow_space,
        quantity: FieldQuantity::Dimensional(flow_space.dims),
        coordinates: coordinates.clone(),
        clock: clock.clone(),
        support: SpatialSupport::BoundaryTrace,
        state: StateOwnership::External,
    })?;
    let power = system.declare_field(FieldDecl {
        name: "port-power".to_string(),
        space: Space {
            degree: 0,
            n: 1,
            dims: POWER_DIMS,
        },
        quantity: FieldQuantity::Dimensional(POWER_DIMS),
        coordinates: coordinates.clone(),
        clock: clock.clone(),
        support: SpatialSupport::BoundaryTrace,
        state: StateOwnership::External,
    })?;
    let pairing = SystemExpr::PortPair {
        kind: spec.schema.kind(),
        effort: Box::new(SystemExpr::FieldRef(effort)),
        flow: Box::new(SystemExpr::FieldRef(flow)),
        measure_dims,
    };
    let rhs = match spec.sense {
        PortEquationSense::AsDeclared => pairing,
        PortEquationSense::Reversed => SystemExpr::Scale(-1.0, Box::new(pairing)),
    };
    system.add_equation(BlockEquation {
        name: "port-power-balance".to_string(),
        target: power,
        rhs,
    })?;
    if let Some((storage, action)) = spec.storage() {
        let state_space = Space {
            degree: RAW_VECTOR_DEGREE,
            n: storage.state_dimension().get(),
            dims: action.state_coordinate_dims(),
        };
        let state = system.declare_field(FieldDecl {
            name: "storage-state".to_string(),
            space: state_space,
            quantity: FieldQuantity::Dimensional(state_space.dims),
            coordinates: coordinates.clone(),
            clock: clock.clone(),
            support: SpatialSupport::Interior,
            state: StateOwnership::Owned { slot: 0 },
        })?;
        let (gradient_target, gradient_space) = match action.gradient_target() {
            StorageGradientTarget::Effort => (effort, effort_space),
            StorageGradientTarget::Flow => (flow, flow_space),
        };
        let gradient_atom = system.register_atom(AtomSignature {
            name: "storage-constitutive-gradient".to_string(),
            content: ConventionRef::new(storage.constitutive_gradient().as_str().to_string())?,
            in_space: state_space,
            out_space: gradient_space,
        });
        system.add_equation(BlockEquation {
            name: "storage-constitutive-gradient".to_string(),
            target: gradient_target,
            rhs: SystemExpr::Apply {
                atom: gradient_atom,
                arg: Box::new(SystemExpr::FieldRef(state)),
            },
        })?;
    }
    let system = system.admit()?;
    let primitive_kind = spec.primitive_kind();
    let primitive_id = spec.primitive_id().map(|id| id.as_str().to_string());
    let conservative_junction_side = spec.conservative_junction_side();
    let reversible_skew_side = spec.reversible_skew_side();
    let (
        reversible_skew_forward_operator,
        reversible_skew_adjoint_operator,
        reversible_skew_evidence,
    ) = spec
        .reversible_skew_coupling()
        .map_or((None, None, None), |coupling| {
            (
                Some(coupling.forward_operator().as_str().to_string()),
                Some(coupling.adjoint_operator().as_str().to_string()),
                Some(coupling.skew_adjoint_evidence().as_str().to_string()),
            )
        });
    let (
        storage_state_schema,
        storage_state_dimension,
        storage_state_action,
        storage_gradient_operator,
    ) = spec
        .storage()
        .map_or((None, None, None, None), |(storage, action)| {
            (
                Some(storage.state_schema().as_str().to_string()),
                Some(storage.state_dimension().get()),
                Some(action),
                Some(storage.constitutive_gradient().as_str().to_string()),
            )
        });
    let receipt = PortEquationReceipt {
        port_id: spec.schema.id().as_str().to_string(),
        port_schema_version: spec.schema.version(),
        system_identity: system.identity(),
        kind: spec.schema.kind(),
        pairing: spec.schema.power_pairing(),
        effort_dims: spec.schema.effort_dimensions(),
        flow_dims: spec.schema.flow_dimensions(),
        measure_dims,
        product_dims,
        sense: spec.sense,
        term_kind: spec.term_kind,
        ownership: spec.ownership,
        loss_ownership_id,
        primitive_kind,
        primitive_id,
        conservative_junction_side,
        reversible_skew_side,
        reversible_skew_forward_operator,
        reversible_skew_adjoint_operator,
        reversible_skew_evidence,
        storage_state_schema,
        storage_state_dimension,
        storage_state_action,
        storage_gradient_operator,
    };
    Ok(CompiledPortEquation { system, receipt })
}

fn validate_stream_version(spec: &StreamEquationSpec) -> Result<(), PortEquationError> {
    if spec.stream.version() == STREAM_PORT_VERSION {
        Ok(())
    } else {
        Err(PortEquationError::StreamSchemaVersionMismatch {
            expected: STREAM_PORT_VERSION,
            actual: spec.stream.version(),
        })
    }
}

fn compile_stream_one(
    spec: StreamEquationSpec,
) -> Result<CompiledStreamEquation, PortEquationError> {
    validate_stream_version(&spec)?;
    let extension = encode_stream_extension(&spec)?;
    let binding = spec.stream.binding();
    let coordinates = CoordinateConvention {
        basis: ConventionRef::new(binding.coordinates().basis().as_str().to_string())?,
        frame: ConventionRef::new(binding.coordinates().frame().as_str().to_string())?,
        orientation: binding.coordinates().orientation(),
    };
    let clock = ConventionRef::new(binding.timestamp().clock().as_str().to_string())?;
    let port_id = binding.port_id().as_str().to_string();
    let constituent_count = spec.stream.constituent_flows().len();
    let mut system = SystemDef::new().with_extension(extension)?;
    for (name, n, dims) in [
        ("stream-mass-flow", 1, MASS_FLOW_DIMS),
        (
            "stream-constituent-flow",
            constituent_count,
            AMOUNT_FLOW_DIMS,
        ),
        ("stream-momentum-flow", 3, MOMENTUM_FLOW_DIMS),
        ("stream-energy-flow", 1, POWER_DIMS),
        ("stream-entropy-flow", 1, ENTROPY_FLOW_DIMS),
    ] {
        let space = Space {
            degree: RAW_VECTOR_DEGREE,
            n,
            dims,
        };
        system.declare_field(FieldDecl {
            name: name.to_string(),
            space,
            quantity: FieldQuantity::Dimensional(dims),
            coordinates: coordinates.clone(),
            clock: clock.clone(),
            support: SpatialSupport::BoundaryTrace,
            state: StateOwnership::External,
        })?;
    }
    let system = system.admit()?;
    let stream_schema_version = spec.stream.version();
    let source_receipt = spec.source_receipt.as_str().to_string();
    let chart_kind = StreamEnergyChartKind::from_chart(spec.stream.energy_chart());
    let energy_flow_bits = spec.stream.energy_flow().value().to_bits();
    let receipt = StreamEquationReceipt {
        port_id,
        stream_schema_version,
        source_receipt,
        system_identity: system.identity(),
        chart_kind,
        constituent_count,
        energy_flow_bits,
        energy_ownership: spec.energy_ownership,
    };
    Ok(CompiledStreamEquation { system, receipt })
}

fn derive_loss_ownership_id(
    spec: &PortEquationSpec,
    effort_space: Space,
    flow_space: Space,
) -> Result<Option<LossOwnershipId>, PortEquationError> {
    let AccountingTermKind::Dissipation = spec.term_kind else {
        return Ok(None);
    };
    let OwnershipDisposition::Owned(_) = &spec.ownership else {
        return Ok(None);
    };
    let mut canonical = spec.clone();
    canonical.sense = PortEquationSense::AsDeclared;
    let payload = encode_extension(&canonical, effort_space, flow_space)?;
    Ok(Some(LossOwnershipId(fs_blake3::hash_domain(
        LOSS_OWNERSHIP_DOMAIN_V1,
        &payload,
    ))))
}

fn resolve_spaces(
    schema: &PortSchema,
    discretization: PortDiscretization,
) -> Result<(Space, Space), PortEquationError> {
    let (effort_degree, effort_dofs, flow_degree, flow_dofs) =
        match (schema.shape(), discretization) {
            (PortValueShape::Scalar, PortDiscretization::Lumped) => {
                (RAW_VECTOR_DEGREE, 1, RAW_VECTOR_DEGREE, 1)
            }
            (PortValueShape::Vector(components), PortDiscretization::Lumped) => (
                RAW_VECTOR_DEGREE,
                components.get(),
                RAW_VECTOR_DEGREE,
                components.get(),
            ),
            (PortValueShape::Tensor { rows, columns }, PortDiscretization::Lumped) => {
                let components = rows
                    .get()
                    .checked_mul(columns.get())
                    .ok_or(PortEquationError::ShapeExtentOverflow)?;
                (RAW_VECTOR_DEGREE, components, RAW_VECTOR_DEGREE, components)
            }
            (
                PortValueShape::Field {
                    components,
                    effort_space,
                    flow_space,
                },
                PortDiscretization::Field {
                    effort_dofs,
                    flow_dofs,
                },
            ) => {
                for (variable, dofs) in [("effort", effort_dofs.get()), ("flow", flow_dofs.get())] {
                    if dofs % components.get() != 0 {
                        return Err(PortEquationError::FieldComponentMismatch {
                            variable,
                            dofs,
                            components: components.get(),
                        });
                    }
                }
                (
                    effort_space.form_degree(),
                    effort_dofs.get(),
                    flow_space.form_degree(),
                    flow_dofs.get(),
                )
            }
            (PortValueShape::Field { .. }, PortDiscretization::Lumped) => {
                return Err(PortEquationError::DiscretizationMismatch {
                    expected: "field",
                    actual: "lumped",
                });
            }
            (_, PortDiscretization::Field { .. }) => {
                return Err(PortEquationError::DiscretizationMismatch {
                    expected: "lumped",
                    actual: "field",
                });
            }
        };
    Ok((
        Space {
            degree: effort_degree,
            n: effort_dofs,
            dims: schema.effort_dimensions(),
        },
        Space {
            degree: flow_degree,
            n: flow_dofs,
            dims: schema.flow_dimensions(),
        },
    ))
}

fn coordinate_convention(schema: &PortSchema) -> Result<CoordinateConvention, PortEquationError> {
    Ok(CoordinateConvention {
        basis: ConventionRef::new(schema.coordinates().basis().as_str().to_string())?,
        frame: ConventionRef::new(schema.coordinates().frame().as_str().to_string())?,
        orientation: schema.coordinates().orientation(),
    })
}

fn pairing_measure(pairing: PowerPairing) -> Dims {
    match pairing {
        PowerPairing::ScalarProduct | PowerPairing::EuclideanDot => Dims::NONE,
        PowerPairing::FieldDuality {
            measure_dimensions, ..
        } => measure_dimensions,
    }
}

fn encode_stream_extension(spec: &StreamEquationSpec) -> Result<Vec<u8>, PortEquationError> {
    let estimated = stream_extension_estimate(spec)?;
    let mut bytes = Vec::new();
    bytes
        .try_reserve_exact(estimated)
        .map_err(|_| PortEquationError::Resource)?;
    bytes.extend_from_slice(&STREAM_EXTENSION_VERSION.to_le_bytes());
    push_text(&mut bytes, "stream-bundle");
    push_text(&mut bytes, crate::VERSION);
    bytes.extend_from_slice(&spec.stream.version().to_le_bytes());
    push_text(&mut bytes, spec.source_receipt.as_str());

    let binding = spec.stream.binding();
    push_text(&mut bytes, binding.port_id().as_str());
    push_text(&mut bytes, binding.state_schema().as_str());
    push_text(&mut bytes, binding.constituent_basis().as_str());
    bytes.extend_from_slice(&(binding.constituent_axis().len() as u64).to_le_bytes());
    for constituent in binding.constituent_axis() {
        encode_stream_constituent_id(&mut bytes, constituent);
    }
    push_text(&mut bytes, binding.chemical_reference_state().as_str());
    push_text(&mut bytes, binding.coordinates().basis().as_str());
    push_text(&mut bytes, binding.coordinates().frame().as_str());
    bytes.push(orientation_tag(binding.coordinates().orientation()));
    push_text(&mut bytes, binding.timestamp().clock().as_str());
    bytes.extend_from_slice(&binding.timestamp().tick().to_le_bytes());
    push_text(&mut bytes, binding.gravity_datum().as_str());
    bytes.push(match binding.stress_convention() {
        StreamStressWorkConvention::CauchyTensionPositiveOutwardPower => 0,
    });

    push_f64_bits(&mut bytes, spec.stream.mass_flow().value());
    bytes.extend_from_slice(&(spec.stream.constituent_flows().len() as u64).to_le_bytes());
    for constituent in spec.stream.constituent_flows() {
        encode_stream_constituent_id(&mut bytes, constituent.id());
        push_f64_bits(&mut bytes, constituent.amount_flow().value());
    }
    for component in spec.stream.momentum_flow() {
        push_f64_bits(&mut bytes, component.value());
    }
    push_f64_bits(&mut bytes, spec.stream.energy_flow().value());
    push_f64_bits(&mut bytes, spec.stream.entropy_flow().value());
    bytes.push(stream_chart_kind_tag(StreamEnergyChartKind::from_chart(
        spec.stream.energy_chart(),
    )));
    bytes.extend_from_slice(&(spec.stream.conservation_roles().len() as u64).to_le_bytes());
    for role in spec.stream.conservation_roles() {
        bytes.push(conservation_role_tag(*role));
    }
    encode_stream_ownership(&mut bytes, &spec.energy_ownership);
    if bytes.len() > MAX_SYSTEM_EXTENSION_BYTES {
        return Err(PortEquationError::CompilerMetadataTooLarge {
            bytes: bytes.len(),
            cap: MAX_SYSTEM_EXTENSION_BYTES,
        });
    }
    Ok(bytes)
}

fn stream_extension_estimate(spec: &StreamEquationSpec) -> Result<usize, PortEquationError> {
    let binding = spec.stream.binding();
    let ownership_bytes = match &spec.energy_ownership {
        StreamEnergyOwnership::Owned(owner) => owner.as_str().len(),
        StreamEnergyOwnership::ExplicitlyUnowned { rationale } => rationale.as_str().len(),
    };
    let mut variable_bytes = 0usize;
    for len in [
        "stream-bundle".len(),
        crate::VERSION.len(),
        spec.source_receipt.as_str().len(),
        binding.port_id().as_str().len(),
        binding.state_schema().as_str().len(),
        binding.constituent_basis().as_str().len(),
        binding.chemical_reference_state().as_str().len(),
        binding.coordinates().basis().as_str().len(),
        binding.coordinates().frame().as_str().len(),
        binding.timestamp().clock().as_str().len(),
        binding.gravity_datum().as_str().len(),
        ownership_bytes,
    ] {
        variable_bytes = checked_metadata_add(variable_bytes, len)?;
    }
    for constituent in binding.constituent_axis() {
        variable_bytes = checked_metadata_add(variable_bytes, stream_constituent_len(constituent))?;
        variable_bytes = checked_metadata_add(variable_bytes, 9)?;
    }
    for constituent in spec.stream.constituent_flows() {
        variable_bytes =
            checked_metadata_add(variable_bytes, stream_constituent_len(constituent.id()))?;
        variable_bytes = checked_metadata_add(variable_bytes, 17)?;
    }
    let estimated = checked_metadata_add(variable_bytes, 512)?;
    if estimated > MAX_SYSTEM_EXTENSION_BYTES {
        return Err(PortEquationError::CompilerMetadataTooLarge {
            bytes: estimated,
            cap: MAX_SYSTEM_EXTENSION_BYTES,
        });
    }
    Ok(estimated)
}

fn checked_metadata_add(left: usize, right: usize) -> Result<usize, PortEquationError> {
    left.checked_add(right)
        .ok_or(PortEquationError::CompilerMetadataTooLarge {
            bytes: usize::MAX,
            cap: MAX_SYSTEM_EXTENSION_BYTES,
        })
}

fn stream_constituent_len(constituent: &StreamConstituentId) -> usize {
    match constituent {
        StreamConstituentId::Species(id) => id.as_str().len(),
        StreamConstituentId::Element(id) => id.as_str().len(),
    }
}

fn encode_stream_constituent_id(bytes: &mut Vec<u8>, constituent: &StreamConstituentId) {
    match constituent {
        StreamConstituentId::Species(id) => {
            bytes.push(0);
            push_text(bytes, id.as_str());
        }
        StreamConstituentId::Element(id) => {
            bytes.push(1);
            push_text(bytes, id.as_str());
        }
    }
}

fn encode_stream_ownership(bytes: &mut Vec<u8>, ownership: &StreamEnergyOwnership) {
    match ownership {
        StreamEnergyOwnership::Owned(owner) => {
            bytes.push(0);
            push_text(bytes, owner.as_str());
        }
        StreamEnergyOwnership::ExplicitlyUnowned { rationale } => {
            bytes.push(1);
            push_text(bytes, rationale.as_str());
        }
    }
}

fn push_f64_bits(bytes: &mut Vec<u8>, value: f64) {
    bytes.extend_from_slice(&value.to_bits().to_le_bytes());
}

fn encode_extension(
    spec: &PortEquationSpec,
    effort_space: Space,
    flow_space: Space,
) -> Result<Vec<u8>, PortEquationError> {
    let estimated = extension_estimate(spec)?;
    let mut bytes = Vec::new();
    bytes
        .try_reserve_exact(estimated)
        .map_err(|_| PortEquationError::Resource)?;
    bytes.extend_from_slice(&PORT_EXTENSION_VERSION.to_le_bytes());
    push_text(&mut bytes, crate::VERSION);
    bytes.extend_from_slice(&spec.schema.version().to_le_bytes());
    push_text(&mut bytes, spec.schema.id().as_str());
    bytes.push(port_kind_tag(spec.schema.kind()));
    push_dims(&mut bytes, spec.schema.effort_dimensions());
    push_dims(&mut bytes, spec.schema.flow_dimensions());
    encode_shape(&mut bytes, spec.schema.shape());
    push_text(&mut bytes, spec.schema.coordinates().basis().as_str());
    push_text(&mut bytes, spec.schema.coordinates().frame().as_str());
    bytes.push(orientation_tag(spec.schema.coordinates().orientation()));
    encode_pairing(&mut bytes, spec.schema.power_pairing());
    push_text(&mut bytes, spec.schema.timestamp().clock().as_str());
    bytes.extend_from_slice(&spec.schema.timestamp().tick().to_le_bytes());
    bytes.extend_from_slice(&(spec.schema.conservation_roles().len() as u64).to_le_bytes());
    for role in spec.schema.conservation_roles() {
        bytes.push(conservation_role_tag(*role));
    }
    bytes.push(spec.sense.sign() as u8);
    bytes.push(accounting_kind_tag(spec.term_kind));
    encode_ownership(&mut bytes, &spec.ownership);
    bytes.push(effort_space.degree);
    bytes.extend_from_slice(&(effort_space.n as u64).to_le_bytes());
    bytes.push(flow_space.degree);
    bytes.extend_from_slice(&(flow_space.n as u64).to_le_bytes());
    if let Some(primitive) = &spec.primitive {
        encode_primitive_binding(&mut bytes, primitive);
    }
    if bytes.len() > MAX_SYSTEM_EXTENSION_BYTES {
        return Err(PortEquationError::CompilerMetadataTooLarge {
            bytes: bytes.len(),
            cap: MAX_SYSTEM_EXTENSION_BYTES,
        });
    }
    Ok(bytes)
}

fn extension_estimate(spec: &PortEquationSpec) -> Result<usize, PortEquationError> {
    let ownership_bytes = match &spec.ownership {
        OwnershipDisposition::NotApplicable => 0,
        OwnershipDisposition::Owned(owner) => owner.as_str().len(),
        OwnershipDisposition::ExplicitlyUnowned { rationale } => rationale.as_str().len(),
    };
    let mut variable_bytes = 0usize;
    for len in [
        crate::VERSION.len(),
        spec.schema.id().as_str().len(),
        spec.schema.coordinates().basis().as_str().len(),
        spec.schema.coordinates().frame().as_str().len(),
        spec.schema.timestamp().clock().as_str().len(),
        ownership_bytes,
    ] {
        variable_bytes =
            variable_bytes
                .checked_add(len)
                .ok_or(PortEquationError::CompilerMetadataTooLarge {
                    bytes: usize::MAX,
                    cap: MAX_SYSTEM_EXTENSION_BYTES,
                })?;
    }
    if let Some(primitive) = &spec.primitive {
        variable_bytes =
            checked_metadata_add(variable_bytes, primitive_binding_variable_bytes(primitive)?)?;
    }
    let estimated =
        variable_bytes
            .checked_add(256)
            .ok_or(PortEquationError::CompilerMetadataTooLarge {
                bytes: usize::MAX,
                cap: MAX_SYSTEM_EXTENSION_BYTES,
            })?;
    if estimated > MAX_SYSTEM_EXTENSION_BYTES {
        return Err(PortEquationError::CompilerMetadataTooLarge {
            bytes: estimated,
            cap: MAX_SYSTEM_EXTENSION_BYTES,
        });
    }
    Ok(estimated)
}

fn encode_primitive_binding(bytes: &mut Vec<u8>, primitive: &PortPrimitiveBinding) {
    bytes.push(0x80);
    bytes.extend_from_slice(&PRIMITIVE_BINDING_VERSION.to_le_bytes());
    match primitive {
        PortPrimitiveBinding::ConservativeJunction { junction, side } => {
            bytes.push(3);
            push_text(bytes, junction.id().as_str());
            bytes.push(match side {
                ConservativeJunctionSide::A => 0,
                ConservativeJunctionSide::B => 1,
            });
            let (port_a, port_b) = junction.ports();
            push_text(bytes, port_a.id().as_str());
            push_text(bytes, port_b.id().as_str());
        }
        PortPrimitiveBinding::ReversibleSkewCoupling { coupling, side } => {
            bytes.push(4);
            push_text(bytes, coupling.id().as_str());
            bytes.push(match side {
                ReversibleSkewSide::A => 0,
                ReversibleSkewSide::B => 1,
            });
            let (port_a, port_b) = coupling.ports();
            push_text(bytes, port_a.id().as_str());
            push_text(bytes, port_b.id().as_str());
            push_text(bytes, coupling.forward_operator().as_str());
            push_text(bytes, coupling.adjoint_operator().as_str());
            push_text(bytes, coupling.skew_adjoint_evidence().as_str());
        }
        PortPrimitiveBinding::Storage { primitive, action } => {
            // Tag 5 is the first storage binding that carries a typed state
            // action. Legacy tag 0 bound only opaque metadata and is never
            // emitted by this compiler revision.
            bytes.push(5);
            bytes.extend_from_slice(&STORAGE_ACTION_BINDING_VERSION.to_le_bytes());
            push_text(bytes, primitive.id().as_str());
            bytes.push(match primitive.potential() {
                StoragePotential::Hamiltonian => 0,
                StoragePotential::FreeEnergy => 1,
            });
            push_text(bytes, primitive.state_schema().as_str());
            bytes.extend_from_slice(&(primitive.state_dimension().get() as u64).to_le_bytes());
            push_dims(bytes, action.state_coordinate_dims());
            bytes.push(match action.gradient_target() {
                StorageGradientTarget::Effort => 0,
                StorageGradientTarget::Flow => 1,
            });
            push_text(bytes, primitive.constitutive_gradient().as_str());
        }
        PortPrimitiveBinding::Dissipation(dissipation) => {
            bytes.push(1);
            push_text(bytes, dissipation.id().as_str());
            bytes.push(dissipation_law_tag(dissipation.law()));
            push_text(bytes, dissipation.constitutive_operator().as_str());
            match dissipation.evidence() {
                DissipationEvidence::Monotonicity(evidence) => {
                    bytes.push(0);
                    push_text(bytes, evidence.as_str());
                }
                DissipationEvidence::NonnegativeProduction(evidence) => {
                    bytes.push(1);
                    push_text(bytes, evidence.as_str());
                }
            }
        }
        PortPrimitiveBinding::SourceOrReservoir(source) => {
            bytes.push(2);
            push_text(bytes, source.id().as_str());
            bytes.push(source_class_tag(source.class()));
            push_text(bytes, source.boundary().id().as_str());
            bytes.push(boundary_treatment_tag(source.boundary().treatment()));
            push_text(bytes, source.boundary().coordinates().basis().as_str());
            push_text(bytes, source.boundary().coordinates().frame().as_str());
            bytes.push(orientation_tag(
                source.boundary().coordinates().orientation(),
            ));
        }
    }
}

fn primitive_binding_variable_bytes(
    primitive: &PortPrimitiveBinding,
) -> Result<usize, PortEquationError> {
    match primitive {
        PortPrimitiveBinding::ConservativeJunction { junction, .. } => {
            let (port_a, port_b) = junction.ports();
            [
                junction.id().as_str().len(),
                port_a.id().as_str().len(),
                port_b.id().as_str().len(),
            ]
            .into_iter()
            .try_fold(64usize, checked_metadata_add)
        }
        PortPrimitiveBinding::ReversibleSkewCoupling { coupling, .. } => {
            let (port_a, port_b) = coupling.ports();
            [
                coupling.id().as_str().len(),
                port_a.id().as_str().len(),
                port_b.id().as_str().len(),
                coupling.forward_operator().as_str().len(),
                coupling.adjoint_operator().as_str().len(),
                coupling.skew_adjoint_evidence().as_str().len(),
            ]
            .into_iter()
            .try_fold(64usize, checked_metadata_add)
        }
        PortPrimitiveBinding::Storage { primitive, .. } => [
            primitive.id().as_str().len(),
            primitive.state_schema().as_str().len(),
            primitive.constitutive_gradient().as_str().len(),
        ]
        .into_iter()
        .try_fold(64usize, checked_metadata_add),
        PortPrimitiveBinding::Dissipation(dissipation) => [
            dissipation.id().as_str().len(),
            dissipation.constitutive_operator().as_str().len(),
            match dissipation.evidence() {
                DissipationEvidence::Monotonicity(evidence)
                | DissipationEvidence::NonnegativeProduction(evidence) => evidence.as_str().len(),
            },
        ]
        .into_iter()
        .try_fold(64usize, checked_metadata_add),
        PortPrimitiveBinding::SourceOrReservoir(source) => [
            source.id().as_str().len(),
            source.boundary().id().as_str().len(),
            source.boundary().coordinates().basis().as_str().len(),
            source.boundary().coordinates().frame().as_str().len(),
        ]
        .into_iter()
        .try_fold(64usize, checked_metadata_add),
    }
}

fn push_text(bytes: &mut Vec<u8>, value: &str) {
    bytes.extend_from_slice(&(value.len() as u64).to_le_bytes());
    bytes.extend_from_slice(value.as_bytes());
}

fn push_dims(bytes: &mut Vec<u8>, dims: Dims) {
    bytes.extend(dims.0.map(|exponent| exponent as u8));
}

fn encode_shape(bytes: &mut Vec<u8>, shape: PortValueShape) {
    match shape {
        PortValueShape::Scalar => bytes.push(0),
        PortValueShape::Vector(components) => {
            bytes.push(1);
            bytes.extend_from_slice(&(components.get() as u64).to_le_bytes());
        }
        PortValueShape::Tensor { rows, columns } => {
            bytes.push(2);
            bytes.extend_from_slice(&(rows.get() as u64).to_le_bytes());
            bytes.extend_from_slice(&(columns.get() as u64).to_le_bytes());
        }
        PortValueShape::Field {
            components,
            effort_space,
            flow_space,
        } => {
            bytes.push(3);
            bytes.extend_from_slice(&(components.get() as u64).to_le_bytes());
            bytes.push(space_type_tag(effort_space));
            bytes.push(space_type_tag(flow_space));
        }
    }
}

fn encode_pairing(bytes: &mut Vec<u8>, pairing: PowerPairing) {
    match pairing {
        PowerPairing::ScalarProduct => bytes.push(0),
        PowerPairing::EuclideanDot => bytes.push(1),
        PowerPairing::FieldDuality {
            measure_dimensions,
            measure_side,
        } => {
            bytes.push(2);
            push_dims(bytes, measure_dimensions);
            bytes.push(match measure_side {
                FieldMeasureSide::Effort => 0,
                FieldMeasureSide::Flow => 1,
            });
        }
    }
}

fn encode_ownership(bytes: &mut Vec<u8>, ownership: &OwnershipDisposition) {
    match ownership {
        OwnershipDisposition::NotApplicable => bytes.push(0),
        OwnershipDisposition::Owned(owner) => {
            bytes.push(1);
            push_text(bytes, owner.as_str());
        }
        OwnershipDisposition::ExplicitlyUnowned { rationale } => {
            bytes.push(2);
            push_text(bytes, rationale.as_str());
        }
    }
}

const fn port_kind_tag(kind: PortKind) -> u8 {
    match kind {
        PortKind::MechanicalForceVelocity => 0,
        PortKind::FluidPressureFlux => 1,
        PortKind::ThermalTemperatureEntropy => 2,
        PortKind::RotationalTorqueAngularVelocity => 3,
        PortKind::ElectricalVoltageCurrent => 4,
        PortKind::MagneticMmfFluxRate => 5,
        PortKind::ChemicalPotentialAmountFlow => 6,
    }
}

const fn port_kind_name(kind: PortKind) -> &'static str {
    match kind {
        PortKind::MechanicalForceVelocity => "mechanical-force-velocity",
        PortKind::FluidPressureFlux => "fluid-pressure-flux",
        PortKind::ThermalTemperatureEntropy => "thermal-temperature-entropy",
        PortKind::RotationalTorqueAngularVelocity => "rotational-torque-angular-velocity",
        PortKind::ElectricalVoltageCurrent => "electrical-voltage-current",
        PortKind::MagneticMmfFluxRate => "magnetic-mmf-flux-rate",
        PortKind::ChemicalPotentialAmountFlow => "chemical-potential-amount-flow",
    }
}

const fn pairing_name(pairing: PowerPairing) -> &'static str {
    match pairing {
        PowerPairing::ScalarProduct => "scalar-product",
        PowerPairing::EuclideanDot => "euclidean-dot",
        PowerPairing::FieldDuality { .. } => "field-duality",
    }
}

const fn orientation_tag(orientation: PortOrientation) -> u8 {
    match orientation {
        PortOrientation::OutwardFromOwner => 0,
        PortOrientation::AlongFrame => 1,
        PortOrientation::AgainstFrame => 2,
    }
}

const fn accounting_kind_tag(kind: AccountingTermKind) -> u8 {
    match kind {
        AccountingTermKind::Reversible => 0,
        AccountingTermKind::Storage => 1,
        AccountingTermKind::Source => 2,
        AccountingTermKind::Dissipation => 3,
    }
}

const fn dissipation_law_tag(law: DissipationLaw) -> u8 {
    match law {
        DissipationLaw::Resistive => 0,
        DissipationLaw::Frictional => 1,
        DissipationLaw::Viscous => 2,
        DissipationLaw::Conductive => 3,
        DissipationLaw::Plastic => 4,
    }
}

const fn source_class_tag(class: SourceClass) -> u8 {
    match class {
        SourceClass::PrescribedEffort => 0,
        SourceClass::PrescribedFlow => 1,
        SourceClass::Reservoir => 2,
    }
}

const fn boundary_treatment_tag(treatment: BoundaryTreatment) -> u8 {
    match treatment {
        BoundaryTreatment::IncludedSourceTerm => 0,
        BoundaryTreatment::ExternalReservoirExchange => 1,
    }
}

const fn stream_chart_kind_tag(kind: StreamEnergyChartKind) -> u8 {
    match kind {
        StreamEnergyChartKind::MovingStreamEnthalpy => 0,
        StreamEnergyChartKind::InternalEnergyCauchyWork => 1,
        StreamEnergyChartKind::ConjugatePotential => 2,
    }
}

const fn conservation_role_tag(role: ConservationRole) -> u8 {
    match role {
        ConservationRole::Energy => 0,
        ConservationRole::Mass => 1,
        ConservationRole::Amount => 2,
        ConservationRole::LinearMomentum => 3,
        ConservationRole::AngularMomentum => 4,
        ConservationRole::Entropy => 5,
        ConservationRole::ElectricCharge => 6,
    }
}

const fn space_type_tag(space: SpaceType) -> u8 {
    match space {
        SpaceType::HGrad => 0,
        SpaceType::HCurl => 1,
        SpaceType::HDiv => 2,
        SpaceType::L2 => 3,
    }
}

fn dims_json(dims: Dims) -> String {
    format!(
        "[{},{},{},{},{},{}]",
        dims.0[0], dims.0[1], dims.0[2], dims.0[3], dims.0[4], dims.0[5]
    )
}
