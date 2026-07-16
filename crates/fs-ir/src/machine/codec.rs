//! Canonical FrankenScript AST transport for admitted Machine graph/behavior.
//!
//! This module owns syntax only. It decodes explicit, versioned ASTs into an
//! authority-free [`MachineGraphDraft`] or [`MachineBehaviorDraft`], then
//! delegates every semantic decision to the existing admission boundaries.
//! Behavior syntax binds one exact graph identity and cannot be replayed
//! against another graph. The codec never treats collection position as
//! identity and never publishes admitted identity on a syntax or semantic
//! refusal.

use core::fmt;
use core::num::NonZeroU64;

use fs_blake3::identity::StrongIdentity;
use fs_qty::Dims;
use fs_qty::semantic::{
    AngleDomain, CompositionBasis, QuantityKind, SemanticType, StrainBasis, StrainComponent,
    ValueForm,
};

use crate::VersionedProgram;
use crate::ast::{Node, NodeKind, Span};

use super::semantics::{
    AdmittedMachineBehavior, BodyMotion, ConditionBinding, ConditionHistoryRef, ConditionSource,
    ConditionTarget, ConditionValueRef, CorrelationModelRef, CrossingSemantics, DependenceMember,
    DependenceModel, DependenceSpec, DistributionRef, EventDependency, EventId, EventOrder,
    EventSpec, EventWitnessRef, FiniteNonNegative, GuardOrientation, GuardRef, HistoryContinuity,
    MAX_MACHINE_BEHAVIOR_CONDITIONS, MAX_MACHINE_BEHAVIOR_DEPENDENCES, MAX_MACHINE_BEHAVIOR_EVENTS,
    MAX_MACHINE_BEHAVIOR_MOTIONS, MAX_MACHINE_BEHAVIOR_NESTED_REFERENCES,
    MAX_MACHINE_BEHAVIOR_STATE_CONTRACTS, MAX_MACHINE_BEHAVIOR_TOLERANCES, MachineBehaviorDraft,
    MachineBehaviorRefusal, MotionBinding, MotionPathRef, NoClaimRef, OutcomeSetRef, ParameterRef,
    ResetMapRef, ResetSemantics, SimultaneityGroupRef, StateSlotContract, ToleranceId,
    ToleranceLawRef, ToleranceSemantics, ToleranceSpec, ToleranceTarget,
};
use super::{
    AdmittedMachineGraph, BodyId, ClockId, ClockSpec, ContactFeatureId, FrameBinding,
    InterfaceBinding, InterfaceCardRef, InterfaceId, InterfaceOrientation,
    MAX_MACHINE_ENTITY_KEY_BYTES, MAX_MACHINE_GRAPH_CLOCKS, MAX_MACHINE_GRAPH_INTERFACES,
    MAX_MACHINE_GRAPH_MATERIALS, MAX_MACHINE_GRAPH_OWNED_ELEMENTS, MAX_MACHINE_GRAPH_PORTS,
    MAX_MACHINE_GRAPH_RELATIONS, MAX_MACHINE_GRAPH_SUBSYSTEMS, MAX_MACHINE_GRAPH_TERMINALS,
    MachineClock, MachineElementId, MachineGraphDraft, MachineGraphIdV1, MachineGraphRefusal,
    MachineIdError, MachineReferenceError, MaterialBinding, MaterialCardRef, MaterialTarget,
    ModelRef, OrientationParity, PortEnergyRole, PortId, PortSpec, RelationId, RelationMode,
    RelationSpec, SolvePolicyRef, StateSlotId, SubsystemId, SubsystemSpec, SurfacePatchId,
    TerminalCausality, TerminalId, TerminalQuantitySpec, TerminalShape, TerminalSpec,
};

/// Version of the canonical Machine-graph FrankenScript form.
pub const MACHINE_GRAPH_AST_SCHEMA_VERSION_V1: u32 = 1;
/// Root symbol for the canonical Machine-graph FrankenScript form.
pub const MACHINE_GRAPH_AST_HEAD_V1: &str = "machine-graph-v1";
/// Maximum total generic AST nodes inspected by one Machine-graph decode.
pub const MAX_MACHINE_GRAPH_AST_NODES: usize = 262_144;
/// Maximum aggregate string/symbol/keyword/quantity-text bytes in one decode.
pub const MAX_MACHINE_GRAPH_AST_TEXT_BYTES: usize = 16 * 1_024 * 1_024;

/// Closed syntax/resource refusal vocabulary for the Machine-graph codec.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u8)]
pub enum MachineGraphCodecRule {
    /// The caller supplied a forged or otherwise invalid generic AST.
    InvalidAst = 1,
    /// A list had the wrong head, arity, or atom kind.
    UnexpectedForm = 2,
    /// A closed enum/tag token was unknown.
    UnknownTag = 3,
    /// A numeric atom was noncanonical, out of range, or zero when forbidden.
    InvalidNumber = 4,
    /// A durable Machine identifier was not a canonical role-specific key.
    InvalidIdentifier = 5,
    /// An opaque external reference was malformed or all-zero.
    InvalidReference = 6,
    /// A public collection or aggregate owned-element bound was exceeded.
    ResourceLimit = 7,
}

impl MachineGraphCodecRule {
    /// Stable structured diagnostic code.
    #[must_use]
    pub const fn code(self) -> &'static str {
        match self {
            Self::InvalidAst => "MachineGraphCodecInvalidAst",
            Self::UnexpectedForm => "MachineGraphCodecUnexpectedForm",
            Self::UnknownTag => "MachineGraphCodecUnknownTag",
            Self::InvalidNumber => "MachineGraphCodecInvalidNumber",
            Self::InvalidIdentifier => "MachineGraphCodecInvalidIdentifier",
            Self::InvalidReference => "MachineGraphCodecInvalidReference",
            Self::ResourceLimit => "MachineGraphCodecResourceLimit",
        }
    }
}

/// One bounded, path-addressed Machine-graph syntax refusal.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MachineGraphCodecError {
    rule: MachineGraphCodecRule,
    span: Span,
    path: Box<str>,
    detail: Box<str>,
    hint: Box<str>,
}

impl MachineGraphCodecError {
    /// Closed refusal rule.
    #[must_use]
    pub const fn rule(&self) -> MachineGraphCodecRule {
        self.rule
    }

    /// Stable structured diagnostic code.
    #[must_use]
    pub const fn code(&self) -> &'static str {
        self.rule.code()
    }

    /// Exact source span of the offending AST node.
    #[must_use]
    pub const fn span(&self) -> Span {
        self.span
    }

    /// Deterministic structural path within the Machine form.
    #[must_use]
    pub fn path(&self) -> &str {
        &self.path
    }

    /// Human-readable refusal detail.
    #[must_use]
    pub fn detail(&self) -> &str {
        &self.detail
    }

    /// Actionable correction hint.
    #[must_use]
    pub fn hint(&self) -> &str {
        &self.hint
    }
}

impl fmt::Display for MachineGraphCodecError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{} at {} (bytes {}..{}): {}; hint: {}",
            self.code(),
            self.path,
            self.span.start,
            self.span.end,
            self.detail,
            self.hint
        )
    }
}

impl std::error::Error for MachineGraphCodecError {}

/// Distinguishes syntax refusal from semantic Machine-graph refusal.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MachineGraphAstAdmissionError {
    /// The generic AST could not be decoded as the v1 Machine graph grammar.
    Codec(MachineGraphCodecError),
    /// The decoded authority-free draft failed the existing graph admission.
    Graph(MachineGraphRefusal),
}

impl MachineGraphAstAdmissionError {
    /// Stable top-level refusal code without collapsing syntax and semantics.
    #[must_use]
    pub const fn code(&self) -> &'static str {
        match self {
            Self::Codec(error) => error.code(),
            Self::Graph(_) => "MachineGraphRefused",
        }
    }
}

impl fmt::Display for MachineGraphAstAdmissionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Codec(error) => error.fmt(f),
            Self::Graph(error) => error.fmt(f),
        }
    }
}

impl std::error::Error for MachineGraphAstAdmissionError {}

impl From<MachineGraphCodecError> for MachineGraphAstAdmissionError {
    fn from(value: MachineGraphCodecError) -> Self {
        Self::Codec(value)
    }
}

/// Decode one canonical v1 Machine-graph AST into an authority-free draft.
///
/// The root grammar is positional and closed:
///
/// ```text
/// (machine-graph-v1
///   (clocks ...)
///   (subsystems ...)
///   (terminals ...)
///   (ports ...)
///   (relations ...)
///   (materials ...)
///   (interfaces ...))
/// ```
///
/// Integer-sized semantic fields are encoded as canonical decimal strings so
/// the syntax can retain the complete `u64` domain instead of inheriting the
/// generic AST's signed-integer limit.
///
/// # Errors
/// Refuses malformed ASTs, unknown closed tags, noncanonical numbers or
/// digests, invalid role-specific IDs/references, and oversized collections.
pub fn parse_machine_graph_v1(node: &Node) -> Result<MachineGraphDraft, MachineGraphCodecError> {
    preflight_ast(node)?;
    preflight_graph_caps(node)?;
    node.validate().map_err(|error| MachineGraphCodecError {
        rule: MachineGraphCodecRule::InvalidAst,
        span: error.span,
        path: validation_path(&error.detail),
        detail: error.to_string().into_boxed_str(),
        hint: "repair the generic FrankenScript AST before Machine decoding"
            .to_string()
            .into_boxed_str(),
    })?;

    let root = exact_form(node, MACHINE_GRAPH_AST_HEAD_V1, 7, "$")?;
    let clocks = section_items(&root[0], "clocks", MAX_MACHINE_GRAPH_CLOCKS, "$[1]")?;
    let subsystems = section_items(&root[1], "subsystems", MAX_MACHINE_GRAPH_SUBSYSTEMS, "$[2]")?;
    let terminals = section_items(&root[2], "terminals", MAX_MACHINE_GRAPH_TERMINALS, "$[3]")?;
    let ports = section_items(&root[3], "ports", MAX_MACHINE_GRAPH_PORTS, "$[4]")?;
    let relations = section_items(&root[4], "relations", MAX_MACHINE_GRAPH_RELATIONS, "$[5]")?;
    let materials = section_items(&root[5], "materials", MAX_MACHINE_GRAPH_MATERIALS, "$[6]")?;
    let interfaces = section_items(&root[6], "interfaces", MAX_MACHINE_GRAPH_INTERFACES, "$[7]")?;

    let mut decoded_clocks = reserved_vec(clocks.len(), &root[0], "$[1]")?;
    for (index, item) in clocks.iter().enumerate() {
        decoded_clocks.push(parse_clock(item, &format!("$[1][{}]", index + 1))?);
    }

    let mut owned_elements = 0_usize;
    let mut decoded_subsystems = reserved_vec(subsystems.len(), &root[1], "$[2]")?;
    for (index, item) in subsystems.iter().enumerate() {
        decoded_subsystems.push(parse_subsystem(
            item,
            &format!("$[2][{}]", index + 1),
            &mut owned_elements,
        )?);
    }

    let mut decoded_terminals = reserved_vec(terminals.len(), &root[2], "$[3]")?;
    for (index, item) in terminals.iter().enumerate() {
        decoded_terminals.push(parse_terminal(item, &format!("$[3][{}]", index + 1))?);
    }

    let mut decoded_ports = reserved_vec(ports.len(), &root[3], "$[4]")?;
    for (index, item) in ports.iter().enumerate() {
        decoded_ports.push(parse_port(item, &format!("$[4][{}]", index + 1))?);
    }

    let mut decoded_relations = reserved_vec(relations.len(), &root[4], "$[5]")?;
    for (index, item) in relations.iter().enumerate() {
        decoded_relations.push(parse_relation(item, &format!("$[5][{}]", index + 1))?);
    }

    let mut decoded_materials = reserved_vec(materials.len(), &root[5], "$[6]")?;
    for (index, item) in materials.iter().enumerate() {
        decoded_materials.push(parse_material(item, &format!("$[6][{}]", index + 1))?);
    }

    let mut decoded_interfaces = reserved_vec(interfaces.len(), &root[6], "$[7]")?;
    for (index, item) in interfaces.iter().enumerate() {
        decoded_interfaces.push(parse_interface(item, &format!("$[7][{}]", index + 1))?);
    }

    Ok(MachineGraphDraft {
        clocks: decoded_clocks,
        subsystems: decoded_subsystems,
        terminals: decoded_terminals,
        ports: decoded_ports,
        relations: decoded_relations,
        materials: decoded_materials,
        interfaces: decoded_interfaces,
    })
}

/// Decode the program body of a version-enforced FrankenScript artifact.
///
/// The outer [`VersionedProgram`] has already enforced the global IR version;
/// this boundary separately enforces the local Machine graph v1 grammar.
///
/// # Errors
/// Returns the same bounded syntax/resource refusals as
/// [`parse_machine_graph_v1`].
pub fn parse_machine_graph_program_v1(
    program: &VersionedProgram,
) -> Result<MachineGraphDraft, MachineGraphCodecError> {
    parse_machine_graph_v1(program.program())
}

/// Decode and semantically admit one v1 Machine-graph AST.
///
/// # Errors
/// [`MachineGraphAstAdmissionError::Codec`] names syntax/resource refusal;
/// [`MachineGraphAstAdmissionError::Graph`] retains the complete deterministic
/// finding set from the existing semantic admission boundary.
pub fn admit_machine_graph_ast_v1(
    node: &Node,
) -> Result<AdmittedMachineGraph, MachineGraphAstAdmissionError> {
    parse_machine_graph_v1(node)?
        .admit()
        .map_err(MachineGraphAstAdmissionError::Graph)
}

/// Encode one admitted graph as the canonical v1 Machine-graph AST.
///
/// # Errors
/// Returns a structured internal-boundary error if an admitted value cannot be
/// represented by the declared grammar or the synthesized AST fails its
/// generic invariants.
pub fn write_machine_graph_v1(
    graph: &AdmittedMachineGraph,
) -> Result<Node, MachineGraphCodecError> {
    let root = form(
        MACHINE_GRAPH_AST_HEAD_V1,
        vec![
            section("clocks", graph.clocks().iter().map(write_clock).collect()),
            section(
                "subsystems",
                graph.subsystems().iter().map(write_subsystem).collect(),
            ),
            section(
                "terminals",
                graph.terminals().iter().map(write_terminal).collect(),
            ),
            section("ports", graph.ports().iter().map(write_port).collect()),
            section(
                "relations",
                graph.relations().iter().map(write_relation).collect(),
            ),
            section(
                "materials",
                graph.materials().iter().map(write_material).collect(),
            ),
            section(
                "interfaces",
                graph.interfaces().iter().map(write_interface).collect(),
            ),
        ],
    );
    root.validate().map_err(|error| MachineGraphCodecError {
        rule: MachineGraphCodecRule::InvalidAst,
        span: error.span,
        path: validation_path(&error.detail),
        detail: format!("canonical Machine graph writer produced an invalid AST: {error}")
            .into_boxed_str(),
        hint: "treat this as a codec implementation defect"
            .to_string()
            .into_boxed_str(),
    })?;
    Ok(root)
}

/// Encode one admitted graph into a current-version FrankenScript artifact.
///
/// # Errors
/// Refuses any internal representation defect rather than emitting an
/// unversioned or structurally invalid artifact.
pub fn write_machine_graph_program_v1(
    graph: &AdmittedMachineGraph,
) -> Result<VersionedProgram, MachineGraphCodecError> {
    let node = write_machine_graph_v1(graph)?;
    VersionedProgram::try_current(node).map_err(|error| MachineGraphCodecError {
        rule: MachineGraphCodecRule::InvalidAst,
        span: error.span,
        path: validation_path(&error.detail),
        detail: format!("canonical Machine graph could not enter the IR envelope: {error}")
            .into_boxed_str(),
        hint: "treat this as a codec implementation defect"
            .to_string()
            .into_boxed_str(),
    })
}

/// Version of the canonical Machine-behavior FrankenScript form.
pub const MACHINE_BEHAVIOR_AST_SCHEMA_VERSION_V1: u32 = 1;
/// Root symbol for the canonical Machine-behavior FrankenScript form.
pub const MACHINE_BEHAVIOR_AST_HEAD_V1: &str = "machine-behavior-v1";
const MAX_CANONICAL_BEHAVIOR_AST_NODES_PER_ENTRY: usize = 64;
const MAX_CANONICAL_BEHAVIOR_AST_NODES_PER_NESTED_REFERENCE: usize = 8;
const MACHINE_BEHAVIOR_AST_ENVELOPE_NODES: usize = 64;
const MACHINE_BEHAVIOR_IDENTITY_BYTES_V1: usize = 24 * 1_024 * 1_024;
/// Maximum total generic AST nodes inspected by one Machine-behavior decode.
///
/// The bound covers the fixed root/section envelope, 64 nodes for every
/// top-level entry admitted by the six public semantic caps, and eight nodes
/// for every nested history-event, dependency, reset-write, or dependence
/// member. Those per-row maxima exceed the largest canonical v1 forms, so an
/// admitted behavior writer cannot emit an AST that its decoder rejects.
pub const MAX_MACHINE_BEHAVIOR_AST_NODES: usize = MACHINE_BEHAVIOR_AST_ENVELOPE_NODES
    + MAX_CANONICAL_BEHAVIOR_AST_NODES_PER_ENTRY
        * (MAX_MACHINE_BEHAVIOR_STATE_CONTRACTS
            + MAX_MACHINE_BEHAVIOR_CONDITIONS
            + MAX_MACHINE_BEHAVIOR_MOTIONS
            + MAX_MACHINE_BEHAVIOR_EVENTS
            + MAX_MACHINE_BEHAVIOR_TOLERANCES
            + MAX_MACHINE_BEHAVIOR_DEPENDENCES)
    + MAX_CANONICAL_BEHAVIOR_AST_NODES_PER_NESTED_REFERENCE
        * MAX_MACHINE_BEHAVIOR_NESTED_REFERENCES;
/// Maximum aggregate atom text inspected by one Machine-behavior decode.
///
/// Eight times the canonical identity-byte envelope covers lowercase-hex and
/// decimal expansion of binary fields plus every closed grammar head/tag in
/// the derived node envelope.
pub const MAX_MACHINE_BEHAVIOR_AST_TEXT_BYTES: usize = 8 * MACHINE_BEHAVIOR_IDENTITY_BYTES_V1;

/// Closed syntax/resource refusal vocabulary for the Machine-behavior codec.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u8)]
pub enum MachineBehaviorCodecRule {
    /// The caller supplied a forged or otherwise invalid generic AST.
    InvalidAst = 1,
    /// A list had the wrong head, arity, or atom kind.
    UnexpectedForm = 2,
    /// A closed enum/tag token was unknown.
    UnknownTag = 3,
    /// A numeric atom was noncanonical, out of range, or semantically invalid.
    InvalidNumber = 4,
    /// A durable Machine identifier was not a canonical role-specific key.
    InvalidIdentifier = 5,
    /// An identity or opaque external reference was malformed or all-zero.
    InvalidReference = 6,
    /// A public collection or aggregate nested-reference bound was exceeded.
    ResourceLimit = 7,
}

impl MachineBehaviorCodecRule {
    /// Stable structured diagnostic code.
    #[must_use]
    pub const fn code(self) -> &'static str {
        match self {
            Self::InvalidAst => "MachineBehaviorCodecInvalidAst",
            Self::UnexpectedForm => "MachineBehaviorCodecUnexpectedForm",
            Self::UnknownTag => "MachineBehaviorCodecUnknownTag",
            Self::InvalidNumber => "MachineBehaviorCodecInvalidNumber",
            Self::InvalidIdentifier => "MachineBehaviorCodecInvalidIdentifier",
            Self::InvalidReference => "MachineBehaviorCodecInvalidReference",
            Self::ResourceLimit => "MachineBehaviorCodecResourceLimit",
        }
    }
}

/// One bounded, path-addressed Machine-behavior syntax refusal.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MachineBehaviorCodecError {
    rule: MachineBehaviorCodecRule,
    span: Span,
    path: Box<str>,
    detail: Box<str>,
    hint: Box<str>,
}

impl MachineBehaviorCodecError {
    /// Closed refusal rule.
    #[must_use]
    pub const fn rule(&self) -> MachineBehaviorCodecRule {
        self.rule
    }

    /// Stable structured diagnostic code.
    #[must_use]
    pub const fn code(&self) -> &'static str {
        self.rule.code()
    }

    /// Exact source span of the offending AST node.
    #[must_use]
    pub const fn span(&self) -> Span {
        self.span
    }

    /// Deterministic structural path within the Machine form.
    #[must_use]
    pub fn path(&self) -> &str {
        &self.path
    }

    /// Human-readable refusal detail.
    #[must_use]
    pub fn detail(&self) -> &str {
        &self.detail
    }

    /// Actionable correction hint.
    #[must_use]
    pub fn hint(&self) -> &str {
        &self.hint
    }

    fn from_graph(error: MachineGraphCodecError) -> Self {
        let rule = match error.rule {
            MachineGraphCodecRule::InvalidAst => MachineBehaviorCodecRule::InvalidAst,
            MachineGraphCodecRule::UnexpectedForm => MachineBehaviorCodecRule::UnexpectedForm,
            MachineGraphCodecRule::UnknownTag => MachineBehaviorCodecRule::UnknownTag,
            MachineGraphCodecRule::InvalidNumber => MachineBehaviorCodecRule::InvalidNumber,
            MachineGraphCodecRule::InvalidIdentifier => MachineBehaviorCodecRule::InvalidIdentifier,
            MachineGraphCodecRule::InvalidReference => MachineBehaviorCodecRule::InvalidReference,
            MachineGraphCodecRule::ResourceLimit => MachineBehaviorCodecRule::ResourceLimit,
        };
        Self {
            rule,
            span: error.span,
            path: error.path,
            detail: error.detail,
            hint: error.hint,
        }
    }
}

impl fmt::Display for MachineBehaviorCodecError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{} at {} (bytes {}..{}): {}; hint: {}",
            self.code(),
            self.path,
            self.span.start,
            self.span.end,
            self.detail,
            self.hint
        )
    }
}

impl std::error::Error for MachineBehaviorCodecError {}

/// Decoded behavior syntax plus its explicit base-graph binding.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DecodedMachineBehaviorV1 {
    base_graph: MachineGraphIdV1,
    draft: MachineBehaviorDraft,
}

impl DecodedMachineBehaviorV1 {
    /// Exact graph identity declared by the syntax.
    #[must_use]
    pub const fn base_graph(&self) -> MachineGraphIdV1 {
        self.base_graph
    }

    /// Borrow the authority-free decoded overlay.
    #[must_use]
    pub const fn draft(&self) -> &MachineBehaviorDraft {
        &self.draft
    }

    /// Recover the exact graph binding and authority-free overlay.
    #[must_use]
    pub fn into_parts(self) -> (MachineGraphIdV1, MachineBehaviorDraft) {
        (self.base_graph, self.draft)
    }

    /// Check the explicit graph binding and delegate semantic publication to
    /// the existing behavior-admission boundary.
    ///
    /// # Errors
    /// Refuses a graph-identity mismatch or the complete semantic finding set.
    pub fn admit_against(
        self,
        graph: &AdmittedMachineGraph,
    ) -> Result<AdmittedMachineBehavior, MachineBehaviorAstAdmissionError> {
        if self.base_graph != graph.identity() {
            return Err(MachineBehaviorAstAdmissionError::BaseGraphMismatch {
                declared: self.base_graph,
                provided: graph.identity(),
            });
        }
        self.draft
            .admit_against(graph)
            .map_err(MachineBehaviorAstAdmissionError::Behavior)
    }
}

/// Distinguishes syntax, graph-binding, and semantic behavior refusal.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MachineBehaviorAstAdmissionError {
    /// The generic AST could not be decoded as the v1 behavior grammar.
    Codec(MachineBehaviorCodecError),
    /// The syntax names a different admitted graph than the caller supplied.
    BaseGraphMismatch {
        /// Identity embedded in the behavior syntax.
        declared: MachineGraphIdV1,
        /// Identity of the graph offered for admission.
        provided: MachineGraphIdV1,
    },
    /// The decoded authority-free overlay failed existing semantic admission.
    Behavior(MachineBehaviorRefusal),
}

impl MachineBehaviorAstAdmissionError {
    /// Stable top-level refusal code without collapsing authority boundaries.
    #[must_use]
    pub const fn code(&self) -> &'static str {
        match self {
            Self::Codec(error) => error.code(),
            Self::BaseGraphMismatch { .. } => "MachineBehaviorBaseGraphMismatch",
            Self::Behavior(_) => "MachineBehaviorRefused",
        }
    }
}

impl fmt::Display for MachineBehaviorAstAdmissionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Codec(error) => error.fmt(f),
            Self::BaseGraphMismatch { declared, provided } => write!(
                f,
                "MachineBehaviorBaseGraphMismatch: syntax binds {}, provided graph is {}",
                digest_hex(*declared.as_bytes()),
                digest_hex(*provided.as_bytes())
            ),
            Self::Behavior(error) => error.fmt(f),
        }
    }
}

impl std::error::Error for MachineBehaviorAstAdmissionError {}

impl From<MachineBehaviorCodecError> for MachineBehaviorAstAdmissionError {
    fn from(value: MachineBehaviorCodecError) -> Self {
        Self::Codec(value)
    }
}

/// Decode one canonical v1 Machine-behavior AST into an authority-free draft
/// paired with its exact base-graph identity.
///
/// The root grammar is positional and closed:
///
/// ```text
/// (machine-behavior-v1
///   (base-graph "<lowercase-hex-id>")
///   (state-contracts ...)
///   (conditions ...)
///   (motions ...)
///   (events ...)
///   (tolerances ...)
///   (dependences ...))
/// ```
///
/// # Errors
/// Refuses malformed ASTs, unknown closed tags, invalid identifiers or
/// references, noncanonical numbers, and oversized collections.
pub fn parse_machine_behavior_v1(
    node: &Node,
) -> Result<DecodedMachineBehaviorV1, MachineBehaviorCodecError> {
    preflight_behavior_ast(node)?;
    preflight_behavior_caps(node)?;
    node.validate().map_err(|error| MachineBehaviorCodecError {
        rule: MachineBehaviorCodecRule::InvalidAst,
        span: error.span,
        path: validation_path(&error.detail),
        detail: error.to_string().into_boxed_str(),
        hint: "repair the generic FrankenScript AST before Machine behavior decoding"
            .to_string()
            .into_boxed_str(),
    })?;

    let root = behavior_exact_form(node, MACHINE_BEHAVIOR_AST_HEAD_V1, 7, "$")?;
    let base = behavior_exact_form(&root[0], "base-graph", 1, "$[1]")?;
    let base_digest = behavior_parse_digest(&base[0], "$[1][1]")?;
    if base_digest == [0; 32] {
        return Err(behavior_reference_error(
            &base[0],
            "$[1][1]",
            "base graph identity must not be all zero",
        ));
    }
    let base_graph = MachineGraphIdV1::parse_slice(&base_digest).ok_or_else(|| {
        behavior_reference_error(
            &base[0],
            "$[1][1]",
            "base graph identity must contain exactly 32 bytes",
        )
    })?;

    let state_contracts = behavior_section_items(
        &root[1],
        "state-contracts",
        MAX_MACHINE_BEHAVIOR_STATE_CONTRACTS,
        "$[2]",
    )?;
    let conditions = behavior_section_items(
        &root[2],
        "conditions",
        MAX_MACHINE_BEHAVIOR_CONDITIONS,
        "$[3]",
    )?;
    let motions =
        behavior_section_items(&root[3], "motions", MAX_MACHINE_BEHAVIOR_MOTIONS, "$[4]")?;
    let events = behavior_section_items(&root[4], "events", MAX_MACHINE_BEHAVIOR_EVENTS, "$[5]")?;
    let tolerances = behavior_section_items(
        &root[5],
        "tolerances",
        MAX_MACHINE_BEHAVIOR_TOLERANCES,
        "$[6]",
    )?;
    let dependences = behavior_section_items(
        &root[6],
        "dependences",
        MAX_MACHINE_BEHAVIOR_DEPENDENCES,
        "$[7]",
    )?;

    let mut decoded_state_contracts =
        behavior_reserved_vec(state_contracts.len(), &root[1], "$[2]")?;
    for (index, item) in state_contracts.iter().enumerate() {
        decoded_state_contracts.push(parse_state_contract(item, &format!("$[2][{}]", index + 1))?);
    }
    let mut decoded_conditions = behavior_reserved_vec(conditions.len(), &root[2], "$[3]")?;
    for (index, item) in conditions.iter().enumerate() {
        decoded_conditions.push(parse_condition(item, &format!("$[3][{}]", index + 1))?);
    }
    let mut decoded_motions = behavior_reserved_vec(motions.len(), &root[3], "$[4]")?;
    for (index, item) in motions.iter().enumerate() {
        decoded_motions.push(parse_motion(item, &format!("$[4][{}]", index + 1))?);
    }
    let mut decoded_events = behavior_reserved_vec(events.len(), &root[4], "$[5]")?;
    for (index, item) in events.iter().enumerate() {
        decoded_events.push(parse_event(item, &format!("$[5][{}]", index + 1))?);
    }
    let mut decoded_tolerances = behavior_reserved_vec(tolerances.len(), &root[5], "$[6]")?;
    for (index, item) in tolerances.iter().enumerate() {
        decoded_tolerances.push(parse_tolerance(item, &format!("$[6][{}]", index + 1))?);
    }
    let mut decoded_dependences = behavior_reserved_vec(dependences.len(), &root[6], "$[7]")?;
    for (index, item) in dependences.iter().enumerate() {
        decoded_dependences.push(parse_dependence(item, &format!("$[7][{}]", index + 1))?);
    }

    Ok(DecodedMachineBehaviorV1 {
        base_graph,
        draft: MachineBehaviorDraft {
            state_contracts: decoded_state_contracts,
            conditions: decoded_conditions,
            motions: decoded_motions,
            events: decoded_events,
            tolerances: decoded_tolerances,
            dependences: decoded_dependences,
        },
    })
}

/// Decode a version-enforced FrankenScript program body as Machine behavior.
///
/// # Errors
/// Returns the same bounded syntax/resource refusals as
/// [`parse_machine_behavior_v1`].
pub fn parse_machine_behavior_program_v1(
    program: &VersionedProgram,
) -> Result<DecodedMachineBehaviorV1, MachineBehaviorCodecError> {
    parse_machine_behavior_v1(program.program())
}

/// Decode and semantically admit one Machine-behavior AST against an exact
/// admitted graph.
///
/// # Errors
/// Preserves syntax, explicit graph-binding, and semantic refusal boundaries.
pub fn admit_machine_behavior_ast_v1(
    node: &Node,
    graph: &AdmittedMachineGraph,
) -> Result<AdmittedMachineBehavior, MachineBehaviorAstAdmissionError> {
    parse_machine_behavior_v1(node)?.admit_against(graph)
}

/// Encode one admitted behavior as the canonical v1 Machine-behavior AST.
///
/// # Errors
/// Returns a structured internal-boundary error if the admitted value cannot
/// be represented by the declared grammar.
pub fn write_machine_behavior_v1(
    behavior: &AdmittedMachineBehavior,
) -> Result<Node, MachineBehaviorCodecError> {
    let root = form(
        MACHINE_BEHAVIOR_AST_HEAD_V1,
        vec![
            form(
                "base-graph",
                vec![string_node(&digest_hex(*behavior.base_graph().as_bytes()))],
            ),
            section(
                "state-contracts",
                behavior
                    .state_contracts()
                    .iter()
                    .map(write_state_contract)
                    .collect(),
            ),
            section(
                "conditions",
                behavior.conditions().iter().map(write_condition).collect(),
            ),
            section(
                "motions",
                behavior.motions().iter().map(write_motion).collect(),
            ),
            section(
                "events",
                behavior.events().iter().map(write_event).collect(),
            ),
            section(
                "tolerances",
                behavior.tolerances().iter().map(write_tolerance).collect(),
            ),
            section(
                "dependences",
                behavior
                    .dependences()
                    .iter()
                    .map(write_dependence)
                    .collect(),
            ),
        ],
    );
    preflight_behavior_ast(&root)?;
    preflight_behavior_caps(&root)?;
    root.validate().map_err(|error| MachineBehaviorCodecError {
        rule: MachineBehaviorCodecRule::InvalidAst,
        span: error.span,
        path: validation_path(&error.detail),
        detail: format!("canonical Machine behavior writer produced an invalid AST: {error}")
            .into_boxed_str(),
        hint: "treat this as a codec implementation defect"
            .to_string()
            .into_boxed_str(),
    })?;
    Ok(root)
}

/// Encode one admitted behavior into a current-version FrankenScript artifact.
///
/// # Errors
/// Refuses an internal representation defect rather than emitting invalid or
/// unversioned behavior syntax.
pub fn write_machine_behavior_program_v1(
    behavior: &AdmittedMachineBehavior,
) -> Result<VersionedProgram, MachineBehaviorCodecError> {
    let node = write_machine_behavior_v1(behavior)?;
    VersionedProgram::try_current(node).map_err(|error| MachineBehaviorCodecError {
        rule: MachineBehaviorCodecRule::InvalidAst,
        span: error.span,
        path: validation_path(&error.detail),
        detail: format!("canonical Machine behavior could not enter the IR envelope: {error}")
            .into_boxed_str(),
        hint: "treat this as a codec implementation defect"
            .to_string()
            .into_boxed_str(),
    })
}

fn preflight_ast(node: &Node) -> Result<(), MachineGraphCodecError> {
    let mut stack = Vec::new();
    stack
        .try_reserve_exact(256)
        .map_err(|_| resource_error(node, "$", "AST preflight stack allocation was refused"))?;
    stack.push(node);
    let mut visited = 0_usize;
    let mut text_bytes = 0_usize;
    while let Some(current) = stack.pop() {
        visited = visited
            .checked_add(1)
            .ok_or_else(|| resource_error(current, "$", "AST node count overflowed usize"))?;
        if visited > MAX_MACHINE_GRAPH_AST_NODES {
            return Err(resource_error(
                current,
                "$",
                &format!(
                    "AST node count exceeds cap {MAX_MACHINE_GRAPH_AST_NODES} before semantic decoding"
                ),
            ));
        }
        let added_text = match &current.kind {
            NodeKind::Qty { text, .. }
            | NodeKind::Str(text)
            | NodeKind::Symbol(text)
            | NodeKind::Keyword(text) => text.len(),
            NodeKind::Int(_)
            | NodeKind::Float(_)
            | NodeKind::Count { .. }
            | NodeKind::Seed(_)
            | NodeKind::List(_) => 0,
        };
        text_bytes = text_bytes
            .checked_add(added_text)
            .ok_or_else(|| resource_error(current, "$", "AST text byte count overflowed usize"))?;
        if text_bytes > MAX_MACHINE_GRAPH_AST_TEXT_BYTES {
            return Err(resource_error(
                current,
                "$",
                &format!(
                    "AST text bytes exceed cap {MAX_MACHINE_GRAPH_AST_TEXT_BYTES} before semantic decoding"
                ),
            ));
        }
        if let NodeKind::List(items) = &current.kind {
            let projected = visited
                .checked_add(stack.len())
                .and_then(|count| count.checked_add(items.len()))
                .ok_or_else(|| resource_error(current, "$", "AST node count overflowed usize"))?;
            if projected > MAX_MACHINE_GRAPH_AST_NODES {
                return Err(resource_error(
                    current,
                    "$",
                    &format!(
                        "AST node count exceeds cap {MAX_MACHINE_GRAPH_AST_NODES} before semantic decoding"
                    ),
                ));
            }
            stack.try_reserve(items.len()).map_err(|_| {
                resource_error(current, "$", "AST preflight stack growth was refused")
            })?;
            stack.extend(items);
        }
    }
    Ok(())
}

fn preflight_graph_caps(node: &Node) -> Result<(), MachineGraphCodecError> {
    let NodeKind::List(root) = &node.kind else {
        return Ok(());
    };
    if !matches!(root.first().map(|item| &item.kind), Some(NodeKind::Symbol(head)) if head == MACHINE_GRAPH_AST_HEAD_V1)
    {
        return Ok(());
    }

    let sections = [
        (1_usize, "clocks", MAX_MACHINE_GRAPH_CLOCKS),
        (2, "subsystems", MAX_MACHINE_GRAPH_SUBSYSTEMS),
        (3, "terminals", MAX_MACHINE_GRAPH_TERMINALS),
        (4, "ports", MAX_MACHINE_GRAPH_PORTS),
        (5, "relations", MAX_MACHINE_GRAPH_RELATIONS),
        (6, "materials", MAX_MACHINE_GRAPH_MATERIALS),
        (7, "interfaces", MAX_MACHINE_GRAPH_INTERFACES),
    ];
    for (index, head, cap) in sections {
        let Some(section) = root.get(index) else {
            continue;
        };
        let Some(items) = recognized_form_items(section, head) else {
            continue;
        };
        if items.len() > cap {
            return Err(resource_error(
                section,
                &format!("$[{index}]"),
                &format!(
                    "section {head} contains {} entries, above cap {cap}",
                    items.len()
                ),
            ));
        }
    }

    let Some(subsystems) = root
        .get(2)
        .and_then(|section| recognized_form_items(section, "subsystems"))
    else {
        return Ok(());
    };
    let mut owned = 0_usize;
    for (subsystem_index, subsystem) in subsystems.iter().enumerate() {
        let Some(fields) = recognized_form_items(subsystem, "subsystem") else {
            continue;
        };
        for (field_index, head) in [
            (2_usize, "bodies"),
            (3, "surface-patches"),
            (4, "contact-features"),
            (5, "state-slots"),
        ] {
            let Some(field) = fields.get(field_index) else {
                continue;
            };
            let Some(items) = recognized_form_items(field, head) else {
                continue;
            };
            owned = owned.checked_add(items.len()).ok_or_else(|| {
                resource_error(
                    field,
                    "$[2]",
                    "aggregate owned-element count overflowed usize",
                )
            })?;
            if owned > MAX_MACHINE_GRAPH_OWNED_ELEMENTS {
                return Err(resource_error(
                    field,
                    &format!("$[2][{}][{}]", subsystem_index + 1, field_index + 1),
                    &format!(
                        "aggregate owned-element count {owned} exceeds cap {MAX_MACHINE_GRAPH_OWNED_ELEMENTS}"
                    ),
                ));
            }
        }
    }
    Ok(())
}

pub(super) fn recognized_form_items<'a>(node: &'a Node, expected_head: &str) -> Option<&'a [Node]> {
    let NodeKind::List(items) = &node.kind else {
        return None;
    };
    if matches!(items.first().map(|item| &item.kind), Some(NodeKind::Symbol(head)) if head == expected_head)
    {
        Some(&items[1..])
    } else {
        None
    }
}

fn preflight_behavior_ast(node: &Node) -> Result<(), MachineBehaviorCodecError> {
    let mut stack = Vec::new();
    stack.try_reserve_exact(256).map_err(|_| {
        behavior_resource_error(node, "$", "AST preflight stack allocation was refused")
    })?;
    stack.push(node);
    let mut visited = 0_usize;
    let mut text_bytes = 0_usize;
    while let Some(current) = stack.pop() {
        visited = visited.checked_add(1).ok_or_else(|| {
            behavior_resource_error(current, "$", "AST node count overflowed usize")
        })?;
        if visited > MAX_MACHINE_BEHAVIOR_AST_NODES {
            return Err(behavior_resource_error(
                current,
                "$",
                &format!(
                    "AST node count exceeds cap {MAX_MACHINE_BEHAVIOR_AST_NODES} before semantic decoding"
                ),
            ));
        }
        let added_text = match &current.kind {
            NodeKind::Qty { text, .. }
            | NodeKind::Str(text)
            | NodeKind::Symbol(text)
            | NodeKind::Keyword(text) => text.len(),
            NodeKind::Int(_)
            | NodeKind::Float(_)
            | NodeKind::Count { .. }
            | NodeKind::Seed(_)
            | NodeKind::List(_) => 0,
        };
        text_bytes = text_bytes.checked_add(added_text).ok_or_else(|| {
            behavior_resource_error(current, "$", "AST text byte count overflowed usize")
        })?;
        if text_bytes > MAX_MACHINE_BEHAVIOR_AST_TEXT_BYTES {
            return Err(behavior_resource_error(
                current,
                "$",
                &format!(
                    "AST text bytes exceed cap {MAX_MACHINE_BEHAVIOR_AST_TEXT_BYTES} before semantic decoding"
                ),
            ));
        }
        if let NodeKind::List(items) = &current.kind {
            let projected = visited
                .checked_add(stack.len())
                .and_then(|count| count.checked_add(items.len()))
                .ok_or_else(|| {
                    behavior_resource_error(current, "$", "AST node count overflowed usize")
                })?;
            if projected > MAX_MACHINE_BEHAVIOR_AST_NODES {
                return Err(behavior_resource_error(
                    current,
                    "$",
                    &format!(
                        "AST node count exceeds cap {MAX_MACHINE_BEHAVIOR_AST_NODES} before semantic decoding"
                    ),
                ));
            }
            stack.try_reserve(items.len()).map_err(|_| {
                behavior_resource_error(current, "$", "AST preflight stack growth was refused")
            })?;
            stack.extend(items);
        }
    }
    Ok(())
}

fn preflight_behavior_caps(node: &Node) -> Result<(), MachineBehaviorCodecError> {
    let NodeKind::List(root) = &node.kind else {
        return Ok(());
    };
    if !matches!(root.first().map(|item| &item.kind), Some(NodeKind::Symbol(head)) if head == MACHINE_BEHAVIOR_AST_HEAD_V1)
    {
        return Ok(());
    }

    let sections = [
        (
            2_usize,
            "state-contracts",
            MAX_MACHINE_BEHAVIOR_STATE_CONTRACTS,
        ),
        (3, "conditions", MAX_MACHINE_BEHAVIOR_CONDITIONS),
        (4, "motions", MAX_MACHINE_BEHAVIOR_MOTIONS),
        (5, "events", MAX_MACHINE_BEHAVIOR_EVENTS),
        (6, "tolerances", MAX_MACHINE_BEHAVIOR_TOLERANCES),
        (7, "dependences", MAX_MACHINE_BEHAVIOR_DEPENDENCES),
    ];
    for (index, head, cap) in sections {
        let Some(section) = root.get(index) else {
            continue;
        };
        let Some(items) = recognized_form_items(section, head) else {
            continue;
        };
        if items.len() > cap {
            return Err(behavior_resource_error(
                section,
                &format!("$[{index}]"),
                &format!(
                    "section {head} contains {} entries, above cap {cap}",
                    items.len()
                ),
            ));
        }
    }

    let mut nested = 0_usize;
    if let Some(conditions) = root
        .get(3)
        .and_then(|section| recognized_form_items(section, "conditions"))
    {
        for (condition_index, condition) in conditions.iter().enumerate() {
            let Some(fields) = recognized_form_items(condition, "condition") else {
                continue;
            };
            let Some(source) = fields.get(5) else {
                continue;
            };
            let Some(history_fields) = recognized_form_items(source, "history") else {
                continue;
            };
            let Some(continuity) = history_fields.get(1) else {
                continue;
            };
            let Some(events) = recognized_form_items(continuity, "reset-at-events") else {
                continue;
            };
            behavior_add_nested(
                &mut nested,
                events.len(),
                continuity,
                &format!("$[3][{}][6][2]", condition_index + 1),
            )?;
        }
    }
    if let Some(events) = root
        .get(5)
        .and_then(|section| recognized_form_items(section, "events"))
    {
        for (event_index, event) in events.iter().enumerate() {
            let Some(fields) = recognized_form_items(event, "event") else {
                continue;
            };
            if let Some(dependencies) = fields
                .get(5)
                .and_then(|field| recognized_form_items(field, "dependencies"))
            {
                behavior_add_nested(
                    &mut nested,
                    dependencies.len(),
                    &fields[5],
                    &format!("$[5][{}][6]", event_index + 1),
                )?;
            }
            let Some(reset) = fields.get(6) else {
                continue;
            };
            let writes = recognized_form_items(reset, "deterministic")
                .and_then(|items| items.get(1).map(|writes| (writes, 2_usize)))
                .or_else(|| {
                    recognized_form_items(reset, "set-valued")
                        .and_then(|items| items.get(2).map(|writes| (writes, 3_usize)))
                });
            let Some((writes, writes_index)) = writes else {
                continue;
            };
            let Some(items) = recognized_form_items(writes, "writes") else {
                continue;
            };
            behavior_add_nested(
                &mut nested,
                items.len(),
                writes,
                &format!("$[5][{}][7][{writes_index}]", event_index + 1),
            )?;
        }
    }
    if let Some(dependences) = root
        .get(7)
        .and_then(|section| recognized_form_items(section, "dependences"))
    {
        for (dependence_index, dependence) in dependences.iter().enumerate() {
            let Some(fields) = recognized_form_items(dependence, "dependence") else {
                continue;
            };
            let Some(members) = fields
                .first()
                .and_then(|field| recognized_form_items(field, "members"))
            else {
                continue;
            };
            behavior_add_nested(
                &mut nested,
                members.len(),
                &fields[0],
                &format!("$[7][{}][1]", dependence_index + 1),
            )?;
        }
    }
    Ok(())
}

fn behavior_add_nested(
    total: &mut usize,
    added: usize,
    node: &Node,
    path: &str,
) -> Result<(), MachineBehaviorCodecError> {
    *total = total.checked_add(added).ok_or_else(|| {
        behavior_resource_error(
            node,
            path,
            "aggregate nested-reference count overflowed usize",
        )
    })?;
    if *total > MAX_MACHINE_BEHAVIOR_NESTED_REFERENCES {
        return Err(behavior_resource_error(
            node,
            path,
            &format!(
                "aggregate nested-reference count {} exceeds cap {MAX_MACHINE_BEHAVIOR_NESTED_REFERENCES}",
                *total
            ),
        ));
    }
    Ok(())
}

fn parse_clock(node: &Node, path: &str) -> Result<ClockSpec, MachineGraphCodecError> {
    let args = exact_form(node, "clock", 2, path)?;
    Ok(ClockSpec {
        id: parse_id(&args[0], "clock-id", &format!("{path}[1]"), |key| {
            ClockId::new(key)
        })?,
        clock: parse_clock_mode(&args[1], &format!("{path}[2]"))?,
    })
}

fn parse_clock_mode(node: &Node, path: &str) -> Result<MachineClock, MachineGraphCodecError> {
    let (head, args) = raw_form(node, path)?;
    match head {
        "continuous" => {
            require_arity(node, args, 0, path)?;
            Ok(MachineClock::Continuous)
        }
        "event-driven" => {
            require_arity(node, args, 0, path)?;
            Ok(MachineClock::EventDriven)
        }
        "periodic" => {
            require_arity(node, args, 2, path)?;
            Ok(MachineClock::Periodic {
                period_ns: parse_nonzero_u64(&args[0], &format!("{path}[1]"))?,
                phase_ns: parse_u64(&args[1], &format!("{path}[2]"))?,
            })
        }
        _ => Err(unknown_tag(node, path, "clock mode", head)),
    }
}

fn parse_subsystem(
    node: &Node,
    path: &str,
    owned_elements: &mut usize,
) -> Result<SubsystemSpec, MachineGraphCodecError> {
    let args = exact_form(node, "subsystem", 6, path)?;
    Ok(SubsystemSpec {
        id: parse_id(&args[0], "subsystem-id", &format!("{path}[1]"), |key| {
            SubsystemId::new(key)
        })?,
        model: parse_reference(
            &args[1],
            "model-ref",
            &format!("{path}[2]"),
            |namespace, version, digest| ModelRef::new(namespace, version, digest),
        )?,
        bodies: parse_owned_ids(
            &args[2],
            "bodies",
            "body-id",
            &format!("{path}[3]"),
            owned_elements,
            |key| BodyId::new(key),
        )?,
        surface_patches: parse_owned_ids(
            &args[3],
            "surface-patches",
            "surface-patch-id",
            &format!("{path}[4]"),
            owned_elements,
            |key| SurfacePatchId::new(key),
        )?,
        contact_features: parse_owned_ids(
            &args[4],
            "contact-features",
            "contact-feature-id",
            &format!("{path}[5]"),
            owned_elements,
            |key| ContactFeatureId::new(key),
        )?,
        state_slots: parse_owned_ids(
            &args[5],
            "state-slots",
            "state-slot-id",
            &format!("{path}[6]"),
            owned_elements,
            |key| StateSlotId::new(key),
        )?,
    })
}

fn parse_terminal(node: &Node, path: &str) -> Result<TerminalSpec, MachineGraphCodecError> {
    let args = exact_form(node, "terminal", 7, path)?;
    Ok(TerminalSpec {
        id: parse_id(&args[0], "terminal-id", &format!("{path}[1]"), |key| {
            TerminalId::new(key)
        })?,
        owner: parse_id(&args[1], "subsystem-id", &format!("{path}[2]"), |key| {
            SubsystemId::new(key)
        })?,
        quantity: parse_quantity(&args[2], &format!("{path}[3]"))?,
        shape: parse_shape(&args[3], &format!("{path}[4]"))?,
        causality: parse_causality(&args[4], &format!("{path}[5]"))?,
        clock: parse_id(&args[5], "clock-id", &format!("{path}[6]"), |key| {
            ClockId::new(key)
        })?,
        frame: parse_frame(&args[6], &format!("{path}[7]"))?,
    })
}

pub(super) fn parse_quantity(
    node: &Node,
    path: &str,
) -> Result<TerminalQuantitySpec, MachineGraphCodecError> {
    let (head, args) = raw_form(node, path)?;
    match head {
        "dims" => {
            require_arity(node, args, 6, path)?;
            let mut dims = [0_i8; 6];
            for (index, value) in args.iter().enumerate() {
                let NodeKind::Int(exponent) = &value.kind else {
                    return Err(unexpected(
                        value,
                        &format!("{path}[{}]", index + 1),
                        "dimension exponent must be an integer atom",
                        "use six signed i8 exponents in [m, kg, s, K, A, mol] order",
                    ));
                };
                dims[index] = i8::try_from(*exponent).map_err(|_| {
                    number_error(
                        &args[index],
                        &format!("{path}[{}]", index + 1),
                        "dimension exponent is outside i8",
                    )
                })?;
            }
            Ok(TerminalQuantitySpec::Dimensional(Dims(dims)))
        }
        "semantic" => {
            require_arity(node, args, 2, path)?;
            Ok(TerminalQuantitySpec::Semantic(SemanticType::new(
                parse_quantity_kind(&args[0], &format!("{path}[1]"))?,
                parse_value_form(&args[1], &format!("{path}[2]"))?,
            )))
        }
        _ => Err(unknown_tag(node, path, "terminal quantity", head)),
    }
}

fn parse_quantity_kind(node: &Node, path: &str) -> Result<QuantityKind, MachineGraphCodecError> {
    let (head, args) = raw_form(node, path)?;
    let simple = match head {
        "absolute-temperature" => Some(QuantityKind::AbsoluteTemperature),
        "temperature-difference" => Some(QuantityKind::TemperatureDifference),
        "torque" => Some(QuantityKind::Torque),
        "energy" => Some(QuantityKind::Energy),
        "pressure" => Some(QuantityKind::Pressure),
        "stress" => Some(QuantityKind::Stress),
        "mass" => Some(QuantityKind::Mass),
        "amount" => Some(QuantityKind::Amount),
        "molar-mass" => Some(QuantityKind::MolarMass),
        "mass-concentration" => Some(QuantityKind::MassConcentration),
        "amount-concentration" => Some(QuantityKind::AmountConcentration),
        "entropy" => Some(QuantityKind::Entropy),
        "heat-capacity" => Some(QuantityKind::HeatCapacity),
        "acoustic-pressure" => Some(QuantityKind::AcousticPressure),
        "acoustic-power" => Some(QuantityKind::AcousticPower),
        _ => None,
    };
    if let Some(kind) = simple {
        require_arity(node, args, 0, path)?;
        return Ok(kind);
    }
    match head {
        "angle" => {
            require_arity(node, args, 1, path)?;
            Ok(QuantityKind::Angle(parse_angle_domain(
                &args[0],
                &format!("{path}[1]"),
            )?))
        }
        "angular-velocity" => {
            require_arity(node, args, 1, path)?;
            Ok(QuantityKind::AngularVelocity(parse_angle_domain(
                &args[0],
                &format!("{path}[1]"),
            )?))
        }
        "strain" => {
            require_arity(node, args, 2, path)?;
            Ok(QuantityKind::Strain {
                basis: match symbol(&args[0], &format!("{path}[1]"))? {
                    "tensor" => StrainBasis::Tensor,
                    "engineering" => StrainBasis::Engineering,
                    value => {
                        return Err(unknown_tag(
                            &args[0],
                            &format!("{path}[1]"),
                            "strain basis",
                            value,
                        ));
                    }
                },
                component: match symbol(&args[1], &format!("{path}[2]"))? {
                    "normal" => StrainComponent::Normal,
                    "shear" => StrainComponent::Shear,
                    value => {
                        return Err(unknown_tag(
                            &args[1],
                            &format!("{path}[2]"),
                            "strain component",
                            value,
                        ));
                    }
                },
            })
        }
        "composition" => {
            require_arity(node, args, 1, path)?;
            let basis = match symbol(&args[0], &format!("{path}[1]"))? {
                "mass-fraction" => CompositionBasis::MassFraction,
                "mole-fraction" => CompositionBasis::MoleFraction,
                "volume-fraction" => CompositionBasis::VolumeFraction,
                value => {
                    return Err(unknown_tag(
                        &args[0],
                        &format!("{path}[1]"),
                        "composition basis",
                        value,
                    ));
                }
            };
            Ok(QuantityKind::Composition(basis))
        }
        _ => Err(unknown_tag(node, path, "semantic quantity kind", head)),
    }
}

fn parse_angle_domain(node: &Node, path: &str) -> Result<AngleDomain, MachineGraphCodecError> {
    match symbol(node, path)? {
        "mechanical" => Ok(AngleDomain::Mechanical),
        "electrical" => Ok(AngleDomain::Electrical),
        value => Err(unknown_tag(node, path, "angle domain", value)),
    }
}

fn parse_value_form(node: &Node, path: &str) -> Result<ValueForm, MachineGraphCodecError> {
    match symbol(node, path)? {
        "static" => Ok(ValueForm::Static),
        "instantaneous" => Ok(ValueForm::Instantaneous),
        "peak" => Ok(ValueForm::Peak),
        "rms" => Ok(ValueForm::Rms),
        value => Err(unknown_tag(node, path, "semantic value form", value)),
    }
}

pub(super) fn parse_shape(
    node: &Node,
    path: &str,
) -> Result<TerminalShape, MachineGraphCodecError> {
    let (head, args) = raw_form(node, path)?;
    match head {
        "scalar" => {
            require_arity(node, args, 0, path)?;
            Ok(TerminalShape::Scalar)
        }
        "vector" => {
            require_arity(node, args, 1, path)?;
            Ok(TerminalShape::Vector {
                components: parse_nonzero_u64(&args[0], &format!("{path}[1]"))?,
            })
        }
        "tensor" => {
            require_arity(node, args, 2, path)?;
            Ok(TerminalShape::Tensor {
                rows: parse_nonzero_u64(&args[0], &format!("{path}[1]"))?,
                columns: parse_nonzero_u64(&args[1], &format!("{path}[2]"))?,
            })
        }
        "field-trace" => {
            require_arity(node, args, 1, path)?;
            Ok(TerminalShape::FieldTrace {
                components: parse_nonzero_u64(&args[0], &format!("{path}[1]"))?,
            })
        }
        _ => Err(unknown_tag(node, path, "terminal shape", head)),
    }
}

fn parse_causality(node: &Node, path: &str) -> Result<TerminalCausality, MachineGraphCodecError> {
    match symbol(node, path)? {
        "input" => Ok(TerminalCausality::Input),
        "output" => Ok(TerminalCausality::Output),
        "external-input" => Ok(TerminalCausality::ExternalInput),
        value => Err(unknown_tag(node, path, "terminal causality", value)),
    }
}

pub(super) fn parse_frame(node: &Node, path: &str) -> Result<FrameBinding, MachineGraphCodecError> {
    let args = exact_form(node, "frame", 2, path)?;
    let orientation = match symbol(&args[1], &format!("{path}[2]"))? {
        "preserving" => OrientationParity::Preserving,
        "reversing" => OrientationParity::Reversing,
        value => {
            return Err(unknown_tag(
                &args[1],
                &format!("{path}[2]"),
                "frame orientation",
                value,
            ));
        }
    };
    parse_id(&args[0], "frame-binding", &format!("{path}[1]"), |key| {
        FrameBinding::new(key, orientation)
    })
}

fn parse_port(node: &Node, path: &str) -> Result<PortSpec, MachineGraphCodecError> {
    let args = exact_form(node, "port", 5, path)?;
    let energy_role = match symbol(&args[4], &format!("{path}[5]"))? {
        "into-subsystem" => PortEnergyRole::IntoSubsystem,
        "out-of-subsystem" => PortEnergyRole::OutOfSubsystem,
        value => {
            return Err(unknown_tag(
                &args[4],
                &format!("{path}[5]"),
                "port energy role",
                value,
            ));
        }
    };
    Ok(PortSpec {
        id: parse_id(&args[0], "port-id", &format!("{path}[1]"), |key| {
            PortId::new(key)
        })?,
        owner: parse_id(&args[1], "subsystem-id", &format!("{path}[2]"), |key| {
            SubsystemId::new(key)
        })?,
        effort: parse_id(&args[2], "terminal-id", &format!("{path}[3]"), |key| {
            TerminalId::new(key)
        })?,
        flow: parse_id(&args[3], "terminal-id", &format!("{path}[4]"), |key| {
            TerminalId::new(key)
        })?,
        energy_role,
    })
}

fn parse_relation(node: &Node, path: &str) -> Result<RelationSpec, MachineGraphCodecError> {
    let args = exact_form(node, "relation", 4, path)?;
    Ok(RelationSpec {
        id: parse_id(&args[0], "relation-id", &format!("{path}[1]"), |key| {
            RelationId::new(key)
        })?,
        source: parse_id(&args[1], "terminal-id", &format!("{path}[2]"), |key| {
            TerminalId::new(key)
        })?,
        target: parse_id(&args[2], "terminal-id", &format!("{path}[3]"), |key| {
            TerminalId::new(key)
        })?,
        mode: parse_relation_mode(&args[3], &format!("{path}[4]"))?,
    })
}

fn parse_relation_mode(node: &Node, path: &str) -> Result<RelationMode, MachineGraphCodecError> {
    let (head, args) = raw_form(node, path)?;
    match head {
        "algebraic" if args.is_empty() => Ok(RelationMode::Algebraic { solve_policy: None }),
        "algebraic" if args.len() == 1 => Ok(RelationMode::Algebraic {
            solve_policy: Some(parse_reference(
                &args[0],
                "solve-policy-ref",
                &format!("{path}[1]"),
                |namespace, version, digest| SolvePolicyRef::new(namespace, version, digest),
            )?),
        }),
        "algebraic" => Err(arity_error(node, path, "algebraic", "0 or 1", args.len())),
        "stateful" => {
            require_arity(node, args, 1, path)?;
            Ok(RelationMode::Stateful {
                state_slot: parse_id(&args[0], "state-slot-id", &format!("{path}[1]"), |key| {
                    StateSlotId::new(key)
                })?,
            })
        }
        _ => Err(unknown_tag(node, path, "relation mode", head)),
    }
}

fn parse_material(node: &Node, path: &str) -> Result<MaterialBinding, MachineGraphCodecError> {
    let args = exact_form(node, "material", 2, path)?;
    let (head, target_args) = raw_form(&args[0], &format!("{path}[1]"))?;
    require_arity(&args[0], target_args, 1, &format!("{path}[1]"))?;
    let target = match head {
        "body" => MaterialTarget::Body(parse_id(
            &target_args[0],
            "body-id",
            &format!("{path}[1][1]"),
            |key| BodyId::new(key),
        )?),
        "surface-patch" => MaterialTarget::SurfacePatch(parse_id(
            &target_args[0],
            "surface-patch-id",
            &format!("{path}[1][1]"),
            |key| SurfacePatchId::new(key),
        )?),
        _ => {
            return Err(unknown_tag(
                &args[0],
                &format!("{path}[1]"),
                "material target",
                head,
            ));
        }
    };
    Ok(MaterialBinding {
        target,
        material: parse_reference(
            &args[1],
            "material-card-ref",
            &format!("{path}[2]"),
            |namespace, version, digest| MaterialCardRef::new(namespace, version, digest),
        )?,
    })
}

fn parse_interface(node: &Node, path: &str) -> Result<InterfaceBinding, MachineGraphCodecError> {
    let args = exact_form(node, "interface", 5, path)?;
    let orientation = match symbol(&args[4], &format!("{path}[5]"))? {
        "aligned" => InterfaceOrientation::Aligned,
        "opposed" => InterfaceOrientation::Opposed,
        value => {
            return Err(unknown_tag(
                &args[4],
                &format!("{path}[5]"),
                "interface orientation",
                value,
            ));
        }
    };
    Ok(InterfaceBinding {
        id: parse_id(&args[0], "interface-id", &format!("{path}[1]"), |key| {
            InterfaceId::new(key)
        })?,
        negative: parse_id(&args[1], "port-id", &format!("{path}[2]"), |key| {
            PortId::new(key)
        })?,
        positive: parse_id(&args[2], "port-id", &format!("{path}[3]"), |key| {
            PortId::new(key)
        })?,
        interface: parse_reference(
            &args[3],
            "interface-card-ref",
            &format!("{path}[4]"),
            |namespace, version, digest| InterfaceCardRef::new(namespace, version, digest),
        )?,
        orientation,
    })
}

fn parse_state_contract(
    node: &Node,
    path: &str,
) -> Result<StateSlotContract, MachineBehaviorCodecError> {
    let args = behavior_exact_form(node, "state-contract", 6, path)?;
    Ok(StateSlotContract {
        id: behavior_parse_id(&args[0], "state-slot-id", &format!("{path}[1]"), |key| {
            StateSlotId::new(key)
        })?,
        owner: behavior_parse_id(&args[1], "subsystem-id", &format!("{path}[2]"), |key| {
            SubsystemId::new(key)
        })?,
        quantity: behavior_parse_quantity(&args[2], &format!("{path}[3]"))?,
        shape: behavior_parse_shape(&args[3], &format!("{path}[4]"))?,
        clock: behavior_parse_id(&args[4], "clock-id", &format!("{path}[5]"), |key| {
            ClockId::new(key)
        })?,
        frame: behavior_parse_frame(&args[5], &format!("{path}[6]"))?,
    })
}

fn parse_condition(node: &Node, path: &str) -> Result<ConditionBinding, MachineBehaviorCodecError> {
    let args = behavior_exact_form(node, "condition", 6, path)?;
    Ok(ConditionBinding {
        target: parse_condition_target(&args[0], &format!("{path}[1]"))?,
        quantity: behavior_parse_quantity(&args[1], &format!("{path}[2]"))?,
        shape: behavior_parse_shape(&args[2], &format!("{path}[3]"))?,
        clock: behavior_parse_id(&args[3], "clock-id", &format!("{path}[4]"), |key| {
            ClockId::new(key)
        })?,
        frame: behavior_parse_frame(&args[4], &format!("{path}[5]"))?,
        source: parse_condition_source(&args[5], &format!("{path}[6]"))?,
    })
}

fn parse_condition_target(
    node: &Node,
    path: &str,
) -> Result<ConditionTarget, MachineBehaviorCodecError> {
    let (head, args) = behavior_raw_form(node, path)?;
    match head {
        "initial" => {
            behavior_require_arity(node, args, 1, path)?;
            Ok(ConditionTarget::Initial(behavior_parse_id(
                &args[0],
                "state-slot-id",
                &format!("{path}[1]"),
                |key| StateSlotId::new(key),
            )?))
        }
        "boundary" => {
            behavior_require_arity(node, args, 1, path)?;
            Ok(ConditionTarget::Boundary(behavior_parse_id(
                &args[0],
                "terminal-id",
                &format!("{path}[1]"),
                |key| TerminalId::new(key),
            )?))
        }
        _ => Err(behavior_unknown_tag(node, path, "condition target", head)),
    }
}

fn parse_condition_source(
    node: &Node,
    path: &str,
) -> Result<ConditionSource, MachineBehaviorCodecError> {
    let (head, args) = behavior_raw_form(node, path)?;
    match head {
        "fixed" => {
            behavior_require_arity(node, args, 1, path)?;
            Ok(ConditionSource::Fixed(behavior_parse_reference(
                &args[0],
                "condition-value-ref",
                &format!("{path}[1]"),
                |namespace, version, digest| ConditionValueRef::new(namespace, version, digest),
            )?))
        }
        "history" => {
            behavior_require_arity(node, args, 2, path)?;
            Ok(ConditionSource::History {
                history: behavior_parse_reference(
                    &args[0],
                    "condition-history-ref",
                    &format!("{path}[1]"),
                    |namespace, version, digest| {
                        ConditionHistoryRef::new(namespace, version, digest)
                    },
                )?,
                continuity: parse_history_continuity(&args[1], &format!("{path}[2]"))?,
            })
        }
        "distribution" => {
            behavior_require_arity(node, args, 1, path)?;
            Ok(ConditionSource::Distribution(behavior_parse_reference(
                &args[0],
                "distribution-ref",
                &format!("{path}[1]"),
                |namespace, version, digest| DistributionRef::new(namespace, version, digest),
            )?))
        }
        _ => Err(behavior_unknown_tag(node, path, "condition source", head)),
    }
}

fn parse_history_continuity(
    node: &Node,
    path: &str,
) -> Result<HistoryContinuity, MachineBehaviorCodecError> {
    let (head, args) = behavior_raw_form(node, path)?;
    match head {
        "continuous" => {
            behavior_require_arity(node, args, 0, path)?;
            Ok(HistoryContinuity::Continuous)
        }
        "reset-at-events" => {
            let mut events = behavior_reserved_vec(args.len(), node, path)?;
            for (index, event) in args.iter().enumerate() {
                events.push(behavior_parse_id(
                    event,
                    "event-id",
                    &format!("{path}[{}]", index + 1),
                    |key| EventId::new(key),
                )?);
            }
            Ok(HistoryContinuity::ResetAtEvents { events })
        }
        _ => Err(behavior_unknown_tag(node, path, "history continuity", head)),
    }
}

fn parse_motion(node: &Node, path: &str) -> Result<MotionBinding, MachineBehaviorCodecError> {
    let args = behavior_exact_form(node, "motion", 4, path)?;
    Ok(MotionBinding {
        body: behavior_parse_id(&args[0], "body-id", &format!("{path}[1]"), |key| {
            BodyId::new(key)
        })?,
        clock: behavior_parse_id(&args[1], "clock-id", &format!("{path}[2]"), |key| {
            ClockId::new(key)
        })?,
        reference_frame: behavior_parse_frame(&args[2], &format!("{path}[3]"))?,
        motion: parse_body_motion(&args[3], &format!("{path}[4]"))?,
    })
}

fn parse_body_motion(node: &Node, path: &str) -> Result<BodyMotion, MachineBehaviorCodecError> {
    let (head, args) = behavior_raw_form(node, path)?;
    match head {
        "static" => {
            behavior_require_arity(node, args, 0, path)?;
            Ok(BodyMotion::Static)
        }
        "prescribed" => {
            behavior_require_arity(node, args, 1, path)?;
            Ok(BodyMotion::Prescribed {
                path: behavior_parse_reference(
                    &args[0],
                    "motion-path-ref",
                    &format!("{path}[1]"),
                    |namespace, version, digest| MotionPathRef::new(namespace, version, digest),
                )?,
            })
        }
        _ => Err(behavior_unknown_tag(node, path, "body motion", head)),
    }
}

fn parse_event(node: &Node, path: &str) -> Result<EventSpec, MachineBehaviorCodecError> {
    let args = behavior_exact_form(node, "event", 8, path)?;
    Ok(EventSpec {
        id: behavior_parse_id(&args[0], "event-id", &format!("{path}[1]"), |key| {
            EventId::new(key)
        })?,
        clock: behavior_parse_id(&args[1], "clock-id", &format!("{path}[2]"), |key| {
            ClockId::new(key)
        })?,
        guard: behavior_parse_reference(
            &args[2],
            "guard-ref",
            &format!("{path}[3]"),
            |namespace, version, digest| GuardRef::new(namespace, version, digest),
        )?,
        orientation: parse_guard_orientation(&args[3], &format!("{path}[4]"))?,
        crossing: parse_crossing(&args[4], &format!("{path}[5]"))?,
        dependencies: parse_event_dependencies(&args[5], &format!("{path}[6]"))?,
        reset: parse_reset(&args[6], &format!("{path}[7]"))?,
        order: parse_event_order(&args[7], &format!("{path}[8]"))?,
    })
}

fn parse_guard_orientation(
    node: &Node,
    path: &str,
) -> Result<GuardOrientation, MachineBehaviorCodecError> {
    match behavior_symbol(node, path)? {
        "negative-to-positive" => Ok(GuardOrientation::NegativeToPositive),
        "positive-to-negative" => Ok(GuardOrientation::PositiveToNegative),
        "bidirectional" => Ok(GuardOrientation::Bidirectional),
        value => Err(behavior_unknown_tag(node, path, "guard orientation", value)),
    }
}

fn parse_crossing(node: &Node, path: &str) -> Result<CrossingSemantics, MachineBehaviorCodecError> {
    let (head, args) = behavior_raw_form(node, path)?;
    behavior_require_arity(node, args, 1, path)?;
    match head {
        "transverse" => Ok(CrossingSemantics::Transverse(behavior_parse_reference(
            &args[0],
            "event-witness-ref",
            &format!("{path}[1]"),
            |namespace, version, digest| EventWitnessRef::new(namespace, version, digest),
        )?)),
        "grazing" => Ok(CrossingSemantics::Grazing(behavior_parse_reference(
            &args[0],
            "event-witness-ref",
            &format!("{path}[1]"),
            |namespace, version, digest| EventWitnessRef::new(namespace, version, digest),
        )?)),
        "unknown" => Ok(CrossingSemantics::Unknown(behavior_parse_reference(
            &args[0],
            "no-claim-ref",
            &format!("{path}[1]"),
            |namespace, version, digest| NoClaimRef::new(namespace, version, digest),
        )?)),
        _ => Err(behavior_unknown_tag(node, path, "crossing semantics", head)),
    }
}

fn parse_event_dependencies(
    node: &Node,
    path: &str,
) -> Result<Vec<EventDependency>, MachineBehaviorCodecError> {
    let items = behavior_section_items(
        node,
        "dependencies",
        MAX_MACHINE_BEHAVIOR_NESTED_REFERENCES,
        path,
    )?;
    let mut dependencies = behavior_reserved_vec(items.len(), node, path)?;
    for (index, item) in items.iter().enumerate() {
        let item_path = format!("{path}[{}]", index + 1);
        let (head, args) = behavior_raw_form(item, &item_path)?;
        behavior_require_arity(item, args, 1, &item_path)?;
        dependencies.push(match head {
            "state" => EventDependency::State(behavior_parse_id(
                &args[0],
                "state-slot-id",
                &format!("{item_path}[1]"),
                |key| StateSlotId::new(key),
            )?),
            "terminal" => EventDependency::Terminal(behavior_parse_id(
                &args[0],
                "terminal-id",
                &format!("{item_path}[1]"),
                |key| TerminalId::new(key),
            )?),
            _ => {
                return Err(behavior_unknown_tag(
                    item,
                    &item_path,
                    "event dependency",
                    head,
                ));
            }
        });
    }
    Ok(dependencies)
}

fn parse_reset(node: &Node, path: &str) -> Result<ResetSemantics, MachineBehaviorCodecError> {
    let (head, args) = behavior_raw_form(node, path)?;
    match head {
        "deterministic" => {
            behavior_require_arity(node, args, 2, path)?;
            Ok(ResetSemantics::Deterministic {
                map: behavior_parse_reference(
                    &args[0],
                    "reset-map-ref",
                    &format!("{path}[1]"),
                    |namespace, version, digest| ResetMapRef::new(namespace, version, digest),
                )?,
                writes: parse_reset_writes(&args[1], &format!("{path}[2]"))?,
            })
        }
        "set-valued" => {
            behavior_require_arity(node, args, 3, path)?;
            Ok(ResetSemantics::SetValued {
                relation: behavior_parse_reference(
                    &args[0],
                    "reset-map-ref",
                    &format!("{path}[1]"),
                    |namespace, version, digest| ResetMapRef::new(namespace, version, digest),
                )?,
                outcomes: behavior_parse_reference(
                    &args[1],
                    "outcome-set-ref",
                    &format!("{path}[2]"),
                    |namespace, version, digest| OutcomeSetRef::new(namespace, version, digest),
                )?,
                writes: parse_reset_writes(&args[2], &format!("{path}[3]"))?,
            })
        }
        "terminal" => {
            behavior_require_arity(node, args, 1, path)?;
            Ok(ResetSemantics::Terminal {
                relation: behavior_parse_reference(
                    &args[0],
                    "reset-map-ref",
                    &format!("{path}[1]"),
                    |namespace, version, digest| ResetMapRef::new(namespace, version, digest),
                )?,
            })
        }
        "unknown" => {
            behavior_require_arity(node, args, 1, path)?;
            Ok(ResetSemantics::Unknown {
                no_claim: behavior_parse_reference(
                    &args[0],
                    "no-claim-ref",
                    &format!("{path}[1]"),
                    |namespace, version, digest| NoClaimRef::new(namespace, version, digest),
                )?,
            })
        }
        _ => Err(behavior_unknown_tag(node, path, "reset semantics", head)),
    }
}

fn parse_reset_writes(
    node: &Node,
    path: &str,
) -> Result<Vec<StateSlotId>, MachineBehaviorCodecError> {
    let items =
        behavior_section_items(node, "writes", MAX_MACHINE_BEHAVIOR_NESTED_REFERENCES, path)?;
    let mut writes = behavior_reserved_vec(items.len(), node, path)?;
    for (index, item) in items.iter().enumerate() {
        writes.push(behavior_parse_id(
            item,
            "state-slot-id",
            &format!("{path}[{}]", index + 1),
            |key| StateSlotId::new(key),
        )?);
    }
    Ok(writes)
}

fn parse_event_order(node: &Node, path: &str) -> Result<EventOrder, MachineBehaviorCodecError> {
    let (head, args) = behavior_raw_form(node, path)?;
    match head {
        "total-priority" => {
            behavior_require_arity(node, args, 1, path)?;
            let value = behavior_parse_u64(&args[0], &format!("{path}[1]"))?;
            let microstep = u32::try_from(value).map_err(|_| {
                behavior_number_error(&args[0], &format!("{path}[1]"), "microstep is outside u32")
            })?;
            Ok(EventOrder::TotalPriority { microstep })
        }
        "set-valued" => {
            behavior_require_arity(node, args, 1, path)?;
            Ok(EventOrder::SetValued {
                group: behavior_parse_reference(
                    &args[0],
                    "simultaneity-group-ref",
                    &format!("{path}[1]"),
                    |namespace, version, digest| {
                        SimultaneityGroupRef::new(namespace, version, digest)
                    },
                )?,
            })
        }
        _ => Err(behavior_unknown_tag(node, path, "event order", head)),
    }
}

fn parse_tolerance(node: &Node, path: &str) -> Result<ToleranceSpec, MachineBehaviorCodecError> {
    let args = behavior_exact_form(node, "tolerance", 6, path)?;
    Ok(ToleranceSpec {
        id: behavior_parse_id(&args[0], "tolerance-id", &format!("{path}[1]"), |key| {
            ToleranceId::new(key)
        })?,
        target: parse_tolerance_target(&args[1], &format!("{path}[2]"))?,
        parameter: behavior_parse_reference(
            &args[2],
            "parameter-ref",
            &format!("{path}[3]"),
            |namespace, version, digest| ParameterRef::new(namespace, version, digest),
        )?,
        quantity: behavior_parse_quantity(&args[3], &format!("{path}[4]"))?,
        shape: behavior_parse_shape(&args[4], &format!("{path}[5]"))?,
        semantics: parse_tolerance_semantics(&args[5], &format!("{path}[6]"))?,
    })
}

fn parse_tolerance_target(
    node: &Node,
    path: &str,
) -> Result<ToleranceTarget, MachineBehaviorCodecError> {
    let (head, args) = behavior_raw_form(node, path)?;
    match head {
        "subsystem" => {
            behavior_require_arity(node, args, 1, path)?;
            Ok(ToleranceTarget::Subsystem(behavior_parse_id(
                &args[0],
                "subsystem-id",
                &format!("{path}[1]"),
                |key| SubsystemId::new(key),
            )?))
        }
        "element" => {
            behavior_require_arity(node, args, 1, path)?;
            Ok(ToleranceTarget::Element(parse_machine_element(
                &args[0],
                &format!("{path}[1]"),
            )?))
        }
        _ => Err(behavior_unknown_tag(node, path, "tolerance target", head)),
    }
}

fn parse_machine_element(
    node: &Node,
    path: &str,
) -> Result<MachineElementId, MachineBehaviorCodecError> {
    let (head, args) = behavior_raw_form(node, path)?;
    behavior_require_arity(node, args, 1, path)?;
    match head {
        "body" => Ok(MachineElementId::Body(behavior_parse_id(
            &args[0],
            "body-id",
            &format!("{path}[1]"),
            |key| BodyId::new(key),
        )?)),
        "surface-patch" => Ok(MachineElementId::SurfacePatch(behavior_parse_id(
            &args[0],
            "surface-patch-id",
            &format!("{path}[1]"),
            |key| SurfacePatchId::new(key),
        )?)),
        "contact-feature" => Ok(MachineElementId::ContactFeature(behavior_parse_id(
            &args[0],
            "contact-feature-id",
            &format!("{path}[1]"),
            |key| ContactFeatureId::new(key),
        )?)),
        "terminal" => Ok(MachineElementId::Terminal(behavior_parse_id(
            &args[0],
            "terminal-id",
            &format!("{path}[1]"),
            |key| TerminalId::new(key),
        )?)),
        "port" => Ok(MachineElementId::Port(behavior_parse_id(
            &args[0],
            "port-id",
            &format!("{path}[1]"),
            |key| PortId::new(key),
        )?)),
        "state-slot" => Ok(MachineElementId::StateSlot(behavior_parse_id(
            &args[0],
            "state-slot-id",
            &format!("{path}[1]"),
            |key| StateSlotId::new(key),
        )?)),
        _ => Err(behavior_unknown_tag(node, path, "machine element", head)),
    }
}

fn parse_tolerance_semantics(
    node: &Node,
    path: &str,
) -> Result<ToleranceSemantics, MachineBehaviorCodecError> {
    let (head, args) = behavior_raw_form(node, path)?;
    match head {
        "bounded" => {
            behavior_require_arity(node, args, 3, path)?;
            Ok(ToleranceSemantics::Bounded {
                minus: parse_finite_nonnegative(&args[0], &format!("{path}[1]"))?,
                plus: parse_finite_nonnegative(&args[1], &format!("{path}[2]"))?,
                law: behavior_parse_reference(
                    &args[2],
                    "tolerance-law-ref",
                    &format!("{path}[3]"),
                    |namespace, version, digest| ToleranceLawRef::new(namespace, version, digest),
                )?,
            })
        }
        "random" => {
            behavior_require_arity(node, args, 3, path)?;
            Ok(ToleranceSemantics::Random {
                scale: parse_finite_nonnegative(&args[0], &format!("{path}[1]"))?,
                law: behavior_parse_reference(
                    &args[1],
                    "tolerance-law-ref",
                    &format!("{path}[2]"),
                    |namespace, version, digest| ToleranceLawRef::new(namespace, version, digest),
                )?,
                marginal: behavior_parse_reference(
                    &args[2],
                    "distribution-ref",
                    &format!("{path}[3]"),
                    |namespace, version, digest| DistributionRef::new(namespace, version, digest),
                )?,
            })
        }
        _ => Err(behavior_unknown_tag(
            node,
            path,
            "tolerance semantics",
            head,
        )),
    }
}

fn parse_finite_nonnegative(
    node: &Node,
    path: &str,
) -> Result<FiniteNonNegative, MachineBehaviorCodecError> {
    let NodeKind::Float(value) = &node.kind else {
        return Err(behavior_unexpected(
            node,
            path,
            "expected a finite nonnegative float atom",
            "use an explicit decimal float such as 0.0 or 1.25",
        ));
    };
    FiniteNonNegative::new(*value)
        .map_err(|error| behavior_number_error(node, path, &error.to_string()))
}

fn parse_dependence(node: &Node, path: &str) -> Result<DependenceSpec, MachineBehaviorCodecError> {
    let args = behavior_exact_form(node, "dependence", 2, path)?;
    let members = behavior_section_items(
        &args[0],
        "members",
        MAX_MACHINE_BEHAVIOR_NESTED_REFERENCES,
        &format!("{path}[1]"),
    )?;
    let mut decoded_members =
        behavior_reserved_vec(members.len(), &args[0], &format!("{path}[1]"))?;
    for (index, member) in members.iter().enumerate() {
        decoded_members.push(parse_dependence_member(
            member,
            &format!("{path}[1][{}]", index + 1),
        )?);
    }
    Ok(DependenceSpec {
        members: decoded_members,
        model: parse_dependence_model(&args[1], &format!("{path}[2]"))?,
    })
}

fn parse_dependence_member(
    node: &Node,
    path: &str,
) -> Result<DependenceMember, MachineBehaviorCodecError> {
    let (head, args) = behavior_raw_form(node, path)?;
    behavior_require_arity(node, args, 1, path)?;
    match head {
        "condition" => Ok(DependenceMember::Condition(parse_condition_target(
            &args[0],
            &format!("{path}[1]"),
        )?)),
        "tolerance" => Ok(DependenceMember::Tolerance(behavior_parse_id(
            &args[0],
            "tolerance-id",
            &format!("{path}[1]"),
            |key| ToleranceId::new(key),
        )?)),
        _ => Err(behavior_unknown_tag(node, path, "dependence member", head)),
    }
}

fn parse_dependence_model(
    node: &Node,
    path: &str,
) -> Result<DependenceModel, MachineBehaviorCodecError> {
    let (head, args) = behavior_raw_form(node, path)?;
    match head {
        "independent" => {
            behavior_require_arity(node, args, 0, path)?;
            Ok(DependenceModel::Independent)
        }
        "correlated" => {
            behavior_require_arity(node, args, 1, path)?;
            Ok(DependenceModel::Correlated(behavior_parse_reference(
                &args[0],
                "correlation-model-ref",
                &format!("{path}[1]"),
                |namespace, version, digest| CorrelationModelRef::new(namespace, version, digest),
            )?))
        }
        _ => Err(behavior_unknown_tag(node, path, "dependence model", head)),
    }
}

fn parse_owned_ids<T, F>(
    node: &Node,
    head: &str,
    role: &str,
    path: &str,
    owned_elements: &mut usize,
    mut constructor: F,
) -> Result<Vec<T>, MachineGraphCodecError>
where
    F: FnMut(&str) -> Result<T, MachineIdError>,
{
    let items = section_items(node, head, MAX_MACHINE_GRAPH_OWNED_ELEMENTS, path)?;
    let next = owned_elements.checked_add(items.len()).ok_or_else(|| {
        resource_error(node, path, "aggregate owned-element count overflowed usize")
    })?;
    if next > MAX_MACHINE_GRAPH_OWNED_ELEMENTS {
        return Err(resource_error(
            node,
            path,
            &format!(
                "aggregate owned-element count {next} exceeds cap {MAX_MACHINE_GRAPH_OWNED_ELEMENTS}"
            ),
        ));
    }
    let mut decoded = reserved_vec(items.len(), node, path)?;
    for (index, item) in items.iter().enumerate() {
        decoded.push(parse_id(
            item,
            role,
            &format!("{path}[{}]", index + 1),
            &mut constructor,
        )?);
    }
    *owned_elements = next;
    Ok(decoded)
}

pub(super) fn parse_id<T, F>(
    node: &Node,
    role: &str,
    path: &str,
    constructor: F,
) -> Result<T, MachineGraphCodecError>
where
    F: FnOnce(&str) -> Result<T, MachineIdError>,
{
    let key = string(node, path)?;
    if key.len() > MAX_MACHINE_ENTITY_KEY_BYTES {
        return Err(MachineGraphCodecError {
            rule: MachineGraphCodecRule::InvalidIdentifier,
            span: node.span,
            path: path.to_string().into_boxed_str(),
            detail: format!(
                "{role} exceeds the {MAX_MACHINE_ENTITY_KEY_BYTES}-byte Machine key limit"
            )
            .into_boxed_str(),
            hint: "use a bounded canonical Machine identifier key"
                .to_string()
                .into_boxed_str(),
        });
    }
    constructor(key).map_err(|error| id_error(node, path, role, &error))
}

pub(super) fn parse_reference<T, F>(
    node: &Node,
    role: &str,
    path: &str,
    constructor: F,
) -> Result<T, MachineGraphCodecError>
where
    F: FnOnce(&str, NonZeroU64, [u8; 32]) -> Result<T, MachineReferenceError>,
{
    let args = exact_form(node, "ref", 3, path)?;
    let namespace = string(&args[0], &format!("{path}[1]"))?;
    if namespace.len() > MAX_MACHINE_ENTITY_KEY_BYTES {
        return Err(MachineGraphCodecError {
            rule: MachineGraphCodecRule::InvalidReference,
            span: args[0].span,
            path: format!("{path}[1]").into_boxed_str(),
            detail: format!(
                "{role} namespace exceeds the {MAX_MACHINE_ENTITY_KEY_BYTES}-byte Machine key limit"
            )
            .into_boxed_str(),
            hint: "use a bounded canonical external-reference namespace"
                .to_string()
                .into_boxed_str(),
        });
    }
    let version = parse_nonzero_u64(&args[1], &format!("{path}[2]"))?;
    let digest = parse_digest(&args[2], &format!("{path}[3]"))?;
    constructor(namespace, version, digest).map_err(|error| match error {
        MachineReferenceError::Namespace(source) => MachineGraphCodecError {
            rule: MachineGraphCodecRule::InvalidReference,
            span: args[0].span,
            path: format!("{path}[1]").into_boxed_str(),
            detail: format!("invalid {role} namespace: {source}").into_boxed_str(),
            hint: "use a bounded canonical external-reference namespace"
                .to_string()
                .into_boxed_str(),
        },
        MachineReferenceError::ZeroDigest {
            role: reference_role,
        } => MachineGraphCodecError {
            rule: MachineGraphCodecRule::InvalidReference,
            span: args[2].span,
            path: format!("{path}[3]").into_boxed_str(),
            detail: format!("{reference_role} semantic digest must not be all zero")
                .into_boxed_str(),
            hint: "supply the exact nonzero 32-byte semantic digest as lowercase hex"
                .to_string()
                .into_boxed_str(),
        },
    })
}

pub(super) fn parse_u64(node: &Node, path: &str) -> Result<u64, MachineGraphCodecError> {
    let text = string(node, path)?;
    if text.is_empty()
        || (text.len() > 1 && text.starts_with('0'))
        || !text.bytes().all(|byte| byte.is_ascii_digit())
    {
        return Err(number_error(
            node,
            path,
            "expected canonical unsigned decimal string without sign or leading zero",
        ));
    }
    text.parse::<u64>()
        .map_err(|_| number_error(node, path, "unsigned decimal is outside u64"))
}

fn parse_nonzero_u64(node: &Node, path: &str) -> Result<NonZeroU64, MachineGraphCodecError> {
    let value = parse_u64(node, path)?;
    NonZeroU64::new(value).ok_or_else(|| number_error(node, path, "value must be nonzero"))
}

pub(super) fn parse_digest(node: &Node, path: &str) -> Result<[u8; 32], MachineGraphCodecError> {
    let text = string(node, path)?;
    if text.len() != 64
        || !text
            .bytes()
            .all(|byte| byte.is_ascii_digit() || (b'a'..=b'f').contains(&byte))
    {
        return Err(MachineGraphCodecError {
            rule: MachineGraphCodecRule::InvalidReference,
            span: node.span,
            path: path.to_string().into_boxed_str(),
            detail: "semantic digest must be exactly 64 lowercase hexadecimal characters"
                .to_string()
                .into_boxed_str(),
            hint: "encode the exact 32-byte digest as lowercase hex"
                .to_string()
                .into_boxed_str(),
        });
    }
    let mut digest = [0_u8; 32];
    let bytes = text.as_bytes();
    for (index, slot) in digest.iter_mut().enumerate() {
        let high = hex_nibble(bytes[index * 2]).ok_or_else(|| {
            invalid_digest(
                node,
                path,
                "semantic digest contains a non-hexadecimal byte",
            )
        })?;
        let low = hex_nibble(bytes[index * 2 + 1]).ok_or_else(|| {
            invalid_digest(
                node,
                path,
                "semantic digest contains a non-hexadecimal byte",
            )
        })?;
        *slot = (high << 4) | low;
    }
    Ok(digest)
}

fn hex_nibble(byte: u8) -> Option<u8> {
    match byte {
        b'0'..=b'9' => Some(byte - b'0'),
        b'a'..=b'f' => Some(byte - b'a' + 10),
        _ => None,
    }
}

fn invalid_digest(node: &Node, path: &str, detail: &str) -> MachineGraphCodecError {
    MachineGraphCodecError {
        rule: MachineGraphCodecRule::InvalidReference,
        span: node.span,
        path: path.to_string().into_boxed_str(),
        detail: detail.to_string().into_boxed_str(),
        hint: "encode the exact 32-byte digest as lowercase hex"
            .to_string()
            .into_boxed_str(),
    }
}

fn exact_form<'a>(
    node: &'a Node,
    expected_head: &str,
    expected_args: usize,
    path: &str,
) -> Result<&'a [Node], MachineGraphCodecError> {
    let (head, args) = raw_form(node, path)?;
    if head != expected_head {
        return Err(unexpected(
            node,
            path,
            &format!("expected ({expected_head} ...), found ({head} ...)"),
            &format!("use the exact {expected_head} form in canonical field order"),
        ));
    }
    require_arity(node, args, expected_args, path)?;
    Ok(args)
}

fn raw_form<'a>(
    node: &'a Node,
    path: &str,
) -> Result<(&'a str, &'a [Node]), MachineGraphCodecError> {
    let NodeKind::List(items) = &node.kind else {
        return Err(unexpected(
            node,
            path,
            "expected a list form",
            "wrap the value in the documented headed form",
        ));
    };
    let Some(first) = items.first() else {
        return Err(unexpected(
            node,
            path,
            "empty lists are not valid Machine forms",
            "start the list with its documented symbol",
        ));
    };
    let NodeKind::Symbol(head) = &first.kind else {
        return Err(unexpected(
            first,
            &format!("{path}[0]"),
            "Machine form head must be a symbol",
            "use the documented literal head symbol",
        ));
    };
    Ok((head, &items[1..]))
}

fn section_items<'a>(
    node: &'a Node,
    expected_head: &str,
    max_items: usize,
    path: &str,
) -> Result<&'a [Node], MachineGraphCodecError> {
    let (head, items) = raw_form(node, path)?;
    if head != expected_head {
        return Err(unexpected(
            node,
            path,
            &format!("expected section {expected_head}, found {head}"),
            "retain all seven sections in canonical graph order, including empty sections",
        ));
    }
    if items.len() > max_items {
        return Err(resource_error(
            node,
            path,
            &format!(
                "section {expected_head} contains {} entries, above cap {max_items}",
                items.len()
            ),
        ));
    }
    Ok(items)
}

fn require_arity(
    node: &Node,
    args: &[Node],
    expected: usize,
    path: &str,
) -> Result<(), MachineGraphCodecError> {
    if args.len() == expected {
        Ok(())
    } else {
        let head = node.head().unwrap_or("<non-form>");
        Err(arity_error(
            node,
            path,
            head,
            &expected.to_string(),
            args.len(),
        ))
    }
}

fn arity_error(
    node: &Node,
    path: &str,
    head: &str,
    expected: &str,
    actual: usize,
) -> MachineGraphCodecError {
    unexpected(
        node,
        path,
        &format!("form {head} expects {expected} arguments, found {actual}"),
        "use the exact positional v1 grammar; omitted values require explicit empty sections/forms",
    )
}

pub(super) fn string<'a>(node: &'a Node, path: &str) -> Result<&'a str, MachineGraphCodecError> {
    match &node.kind {
        NodeKind::Str(value) => Ok(value),
        _ => Err(unexpected(
            node,
            path,
            "expected a string atom",
            "quote identifiers, canonical keys, decimal u64 values, and digests",
        )),
    }
}

pub(super) fn symbol<'a>(node: &'a Node, path: &str) -> Result<&'a str, MachineGraphCodecError> {
    match &node.kind {
        NodeKind::Symbol(value) => Ok(value),
        _ => Err(unexpected(
            node,
            path,
            "expected a closed-tag symbol",
            "use one documented lowercase tag without quotes",
        )),
    }
}

pub(super) fn reserved_vec<T>(
    count: usize,
    node: &Node,
    path: &str,
) -> Result<Vec<T>, MachineGraphCodecError> {
    let mut values = Vec::new();
    values
        .try_reserve_exact(count)
        .map_err(|_| resource_error(node, path, "allocation for decoded entries was refused"))?;
    Ok(values)
}

fn unexpected(node: &Node, path: &str, detail: &str, hint: &str) -> MachineGraphCodecError {
    MachineGraphCodecError {
        rule: MachineGraphCodecRule::UnexpectedForm,
        span: node.span,
        path: path.to_string().into_boxed_str(),
        detail: detail.to_string().into_boxed_str(),
        hint: hint.to_string().into_boxed_str(),
    }
}

fn unknown_tag(node: &Node, path: &str, family: &str, value: &str) -> MachineGraphCodecError {
    MachineGraphCodecError {
        rule: MachineGraphCodecRule::UnknownTag,
        span: node.span,
        path: path.to_string().into_boxed_str(),
        detail: format!("unknown {family} tag {value:?}").into_boxed_str(),
        hint: "use one tag from the closed v1 Machine grammar"
            .to_string()
            .into_boxed_str(),
    }
}

fn number_error(node: &Node, path: &str, detail: &str) -> MachineGraphCodecError {
    MachineGraphCodecError {
        rule: MachineGraphCodecRule::InvalidNumber,
        span: node.span,
        path: path.to_string().into_boxed_str(),
        detail: detail.to_string().into_boxed_str(),
        hint: "use the exact bounded numeric spelling required by the v1 grammar"
            .to_string()
            .into_boxed_str(),
    }
}

fn id_error(node: &Node, path: &str, role: &str, error: &MachineIdError) -> MachineGraphCodecError {
    MachineGraphCodecError {
        rule: MachineGraphCodecRule::InvalidIdentifier,
        span: node.span,
        path: path.to_string().into_boxed_str(),
        detail: format!("invalid {role}: {error}").into_boxed_str(),
        hint: "use a bounded canonical role-specific Machine key"
            .to_string()
            .into_boxed_str(),
    }
}

fn resource_error(node: &Node, path: &str, detail: &str) -> MachineGraphCodecError {
    MachineGraphCodecError {
        rule: MachineGraphCodecRule::ResourceLimit,
        span: node.span,
        path: path.to_string().into_boxed_str(),
        detail: detail.to_string().into_boxed_str(),
        hint: "split the machine or reduce the declared collection before admission"
            .to_string()
            .into_boxed_str(),
    }
}

fn behavior_exact_form<'a>(
    node: &'a Node,
    expected_head: &str,
    expected_args: usize,
    path: &str,
) -> Result<&'a [Node], MachineBehaviorCodecError> {
    let (head, args) = behavior_raw_form(node, path)?;
    if head != expected_head {
        return Err(behavior_unexpected(
            node,
            path,
            &format!("expected ({expected_head} ...), found ({head} ...)"),
            &format!("use the exact {expected_head} form in canonical field order"),
        ));
    }
    behavior_require_arity(node, args, expected_args, path)?;
    Ok(args)
}

fn behavior_raw_form<'a>(
    node: &'a Node,
    path: &str,
) -> Result<(&'a str, &'a [Node]), MachineBehaviorCodecError> {
    let NodeKind::List(items) = &node.kind else {
        return Err(behavior_unexpected(
            node,
            path,
            "expected a list form",
            "wrap the value in the documented headed form",
        ));
    };
    let Some(first) = items.first() else {
        return Err(behavior_unexpected(
            node,
            path,
            "empty lists are not valid Machine forms",
            "start the list with its documented symbol",
        ));
    };
    let NodeKind::Symbol(head) = &first.kind else {
        return Err(behavior_unexpected(
            first,
            &format!("{path}[0]"),
            "Machine form head must be a symbol",
            "use the documented literal head symbol",
        ));
    };
    Ok((head, &items[1..]))
}

fn behavior_section_items<'a>(
    node: &'a Node,
    expected_head: &str,
    max_items: usize,
    path: &str,
) -> Result<&'a [Node], MachineBehaviorCodecError> {
    let (head, items) = behavior_raw_form(node, path)?;
    if head != expected_head {
        return Err(behavior_unexpected(
            node,
            path,
            &format!("expected section {expected_head}, found {head}"),
            "retain every behavior section in canonical order, including empty sections",
        ));
    }
    if items.len() > max_items {
        return Err(behavior_resource_error(
            node,
            path,
            &format!(
                "section {expected_head} contains {} entries, above cap {max_items}",
                items.len()
            ),
        ));
    }
    Ok(items)
}

fn behavior_require_arity(
    node: &Node,
    args: &[Node],
    expected: usize,
    path: &str,
) -> Result<(), MachineBehaviorCodecError> {
    if args.len() == expected {
        Ok(())
    } else {
        Err(behavior_unexpected(
            node,
            path,
            &format!(
                "form {} expects {expected} arguments, found {}",
                node.head().unwrap_or("<non-form>"),
                args.len()
            ),
            "use the exact positional v1 grammar; omitted values require explicit empty forms",
        ))
    }
}

fn behavior_string<'a>(node: &'a Node, path: &str) -> Result<&'a str, MachineBehaviorCodecError> {
    string(node, path).map_err(MachineBehaviorCodecError::from_graph)
}

fn behavior_symbol<'a>(node: &'a Node, path: &str) -> Result<&'a str, MachineBehaviorCodecError> {
    symbol(node, path).map_err(MachineBehaviorCodecError::from_graph)
}

fn behavior_reserved_vec<T>(
    count: usize,
    node: &Node,
    path: &str,
) -> Result<Vec<T>, MachineBehaviorCodecError> {
    let mut values = Vec::new();
    values.try_reserve_exact(count).map_err(|_| {
        behavior_resource_error(node, path, "allocation for decoded entries was refused")
    })?;
    Ok(values)
}

fn behavior_parse_id<T, F>(
    node: &Node,
    role: &str,
    path: &str,
    constructor: F,
) -> Result<T, MachineBehaviorCodecError>
where
    F: FnOnce(&str) -> Result<T, MachineIdError>,
{
    let key = behavior_string(node, path)?;
    constructor(key).map_err(|error| MachineBehaviorCodecError {
        rule: MachineBehaviorCodecRule::InvalidIdentifier,
        span: node.span,
        path: path.to_string().into_boxed_str(),
        detail: format!("invalid {role}: {error}").into_boxed_str(),
        hint: "use a bounded canonical role-specific Machine key"
            .to_string()
            .into_boxed_str(),
    })
}

fn behavior_parse_reference<T, F>(
    node: &Node,
    role: &str,
    path: &str,
    constructor: F,
) -> Result<T, MachineBehaviorCodecError>
where
    F: FnOnce(&str, NonZeroU64, [u8; 32]) -> Result<T, MachineReferenceError>,
{
    parse_reference(node, role, path, constructor).map_err(MachineBehaviorCodecError::from_graph)
}

fn behavior_parse_quantity(
    node: &Node,
    path: &str,
) -> Result<TerminalQuantitySpec, MachineBehaviorCodecError> {
    parse_quantity(node, path).map_err(MachineBehaviorCodecError::from_graph)
}

fn behavior_parse_shape(
    node: &Node,
    path: &str,
) -> Result<TerminalShape, MachineBehaviorCodecError> {
    parse_shape(node, path).map_err(MachineBehaviorCodecError::from_graph)
}

fn behavior_parse_frame(
    node: &Node,
    path: &str,
) -> Result<FrameBinding, MachineBehaviorCodecError> {
    parse_frame(node, path).map_err(MachineBehaviorCodecError::from_graph)
}

fn behavior_parse_u64(node: &Node, path: &str) -> Result<u64, MachineBehaviorCodecError> {
    parse_u64(node, path).map_err(MachineBehaviorCodecError::from_graph)
}

fn behavior_parse_digest(node: &Node, path: &str) -> Result<[u8; 32], MachineBehaviorCodecError> {
    parse_digest(node, path).map_err(MachineBehaviorCodecError::from_graph)
}

fn behavior_unexpected(
    node: &Node,
    path: &str,
    detail: &str,
    hint: &str,
) -> MachineBehaviorCodecError {
    MachineBehaviorCodecError {
        rule: MachineBehaviorCodecRule::UnexpectedForm,
        span: node.span,
        path: path.to_string().into_boxed_str(),
        detail: detail.to_string().into_boxed_str(),
        hint: hint.to_string().into_boxed_str(),
    }
}

fn behavior_unknown_tag(
    node: &Node,
    path: &str,
    family: &str,
    value: &str,
) -> MachineBehaviorCodecError {
    MachineBehaviorCodecError {
        rule: MachineBehaviorCodecRule::UnknownTag,
        span: node.span,
        path: path.to_string().into_boxed_str(),
        detail: format!("unknown {family} tag {value:?}").into_boxed_str(),
        hint: "use one tag from the closed v1 Machine behavior grammar"
            .to_string()
            .into_boxed_str(),
    }
}

fn behavior_number_error(node: &Node, path: &str, detail: &str) -> MachineBehaviorCodecError {
    MachineBehaviorCodecError {
        rule: MachineBehaviorCodecRule::InvalidNumber,
        span: node.span,
        path: path.to_string().into_boxed_str(),
        detail: detail.to_string().into_boxed_str(),
        hint: "use the exact bounded numeric spelling required by the v1 behavior grammar"
            .to_string()
            .into_boxed_str(),
    }
}

fn behavior_reference_error(node: &Node, path: &str, detail: &str) -> MachineBehaviorCodecError {
    MachineBehaviorCodecError {
        rule: MachineBehaviorCodecRule::InvalidReference,
        span: node.span,
        path: path.to_string().into_boxed_str(),
        detail: detail.to_string().into_boxed_str(),
        hint: "supply the exact nonzero identity as 64 lowercase hexadecimal characters"
            .to_string()
            .into_boxed_str(),
    }
}

fn behavior_resource_error(node: &Node, path: &str, detail: &str) -> MachineBehaviorCodecError {
    MachineBehaviorCodecError {
        rule: MachineBehaviorCodecRule::ResourceLimit,
        span: node.span,
        path: path.to_string().into_boxed_str(),
        detail: detail.to_string().into_boxed_str(),
        hint: "split the behavior or reduce the declared collection before admission"
            .to_string()
            .into_boxed_str(),
    }
}

pub(super) fn sym(value: &str) -> Node {
    Node::synthetic(NodeKind::Symbol(value.to_string()))
}

pub(super) fn string_node(value: &str) -> Node {
    Node::synthetic(NodeKind::Str(value.to_string()))
}

fn int_node(value: i64) -> Node {
    Node::synthetic(NodeKind::Int(value))
}

fn float_node(value: f64) -> Node {
    Node::synthetic(NodeKind::Float(value))
}

pub(super) fn form(head: &str, args: Vec<Node>) -> Node {
    let mut items = Vec::with_capacity(args.len() + 1);
    items.push(sym(head));
    items.extend(args);
    Node::synthetic(NodeKind::List(items))
}

pub(super) fn section(head: &str, entries: Vec<Node>) -> Node {
    form(head, entries)
}

pub(super) fn u64_node(value: u64) -> Node {
    string_node(&value.to_string())
}

pub(super) fn digest_hex(digest: [u8; 32]) -> String {
    const HEX: &[u8; 16] = b"0123456789abcdef";
    let mut text = String::with_capacity(64);
    for byte in digest {
        text.push(char::from(HEX[usize::from(byte >> 4)]));
        text.push(char::from(HEX[usize::from(byte & 0x0f)]));
    }
    text
}

pub(super) fn validation_path(detail: &str) -> Box<str> {
    detail
        .strip_prefix("invalid AST at ")
        .and_then(|rest| rest.split_once(':').map(|(path, _)| path))
        .filter(|path| path.starts_with('$'))
        .unwrap_or("$")
        .to_string()
        .into_boxed_str()
}

pub(super) fn reference_node(namespace: &str, version: NonZeroU64, digest: [u8; 32]) -> Node {
    form(
        "ref",
        vec![
            string_node(namespace),
            u64_node(version.get()),
            string_node(&digest_hex(digest)),
        ],
    )
}

fn write_clock(clock: &ClockSpec) -> Node {
    let mode = match clock.clock {
        MachineClock::Continuous => form("continuous", Vec::new()),
        MachineClock::Periodic {
            period_ns,
            phase_ns,
        } => form(
            "periodic",
            vec![u64_node(period_ns.get()), u64_node(phase_ns)],
        ),
        MachineClock::EventDriven => form("event-driven", Vec::new()),
    };
    form("clock", vec![string_node(clock.id.canonical_key()), mode])
}

fn write_subsystem(subsystem: &SubsystemSpec) -> Node {
    form(
        "subsystem",
        vec![
            string_node(subsystem.id.canonical_key()),
            reference_node(
                subsystem.model.namespace(),
                subsystem.model.schema_version(),
                subsystem.model.semantic_digest(),
            ),
            section(
                "bodies",
                subsystem
                    .bodies
                    .iter()
                    .map(|id| string_node(id.canonical_key()))
                    .collect(),
            ),
            section(
                "surface-patches",
                subsystem
                    .surface_patches
                    .iter()
                    .map(|id| string_node(id.canonical_key()))
                    .collect(),
            ),
            section(
                "contact-features",
                subsystem
                    .contact_features
                    .iter()
                    .map(|id| string_node(id.canonical_key()))
                    .collect(),
            ),
            section(
                "state-slots",
                subsystem
                    .state_slots
                    .iter()
                    .map(|id| string_node(id.canonical_key()))
                    .collect(),
            ),
        ],
    )
}

fn write_terminal(terminal: &TerminalSpec) -> Node {
    let causality = match terminal.causality {
        TerminalCausality::Input => "input",
        TerminalCausality::Output => "output",
        TerminalCausality::ExternalInput => "external-input",
    };
    form(
        "terminal",
        vec![
            string_node(terminal.id.canonical_key()),
            string_node(terminal.owner.canonical_key()),
            write_quantity(terminal.quantity),
            write_shape(terminal.shape),
            sym(causality),
            string_node(terminal.clock.canonical_key()),
            form(
                "frame",
                vec![
                    string_node(terminal.frame.canonical_key()),
                    sym(match terminal.frame.orientation() {
                        OrientationParity::Preserving => "preserving",
                        OrientationParity::Reversing => "reversing",
                    }),
                ],
            ),
        ],
    )
}

pub(super) fn write_quantity(quantity: TerminalQuantitySpec) -> Node {
    match quantity {
        TerminalQuantitySpec::Dimensional(Dims(dims)) => form(
            "dims",
            dims.into_iter()
                .map(|value| int_node(i64::from(value)))
                .collect(),
        ),
        TerminalQuantitySpec::Semantic(semantic) => form(
            "semantic",
            vec![
                write_quantity_kind(semantic.kind()),
                sym(match semantic.form() {
                    ValueForm::Static => "static",
                    ValueForm::Instantaneous => "instantaneous",
                    ValueForm::Peak => "peak",
                    ValueForm::Rms => "rms",
                }),
            ],
        ),
    }
}

fn write_quantity_kind(kind: QuantityKind) -> Node {
    match kind {
        QuantityKind::AbsoluteTemperature => form("absolute-temperature", Vec::new()),
        QuantityKind::TemperatureDifference => form("temperature-difference", Vec::new()),
        QuantityKind::Angle(domain) => form("angle", vec![write_angle_domain(domain)]),
        QuantityKind::AngularVelocity(domain) => {
            form("angular-velocity", vec![write_angle_domain(domain)])
        }
        QuantityKind::Torque => form("torque", Vec::new()),
        QuantityKind::Energy => form("energy", Vec::new()),
        QuantityKind::Pressure => form("pressure", Vec::new()),
        QuantityKind::Stress => form("stress", Vec::new()),
        QuantityKind::Strain { basis, component } => form(
            "strain",
            vec![
                sym(match basis {
                    StrainBasis::Tensor => "tensor",
                    StrainBasis::Engineering => "engineering",
                }),
                sym(match component {
                    StrainComponent::Normal => "normal",
                    StrainComponent::Shear => "shear",
                }),
            ],
        ),
        QuantityKind::Composition(basis) => form(
            "composition",
            vec![sym(match basis {
                CompositionBasis::MassFraction => "mass-fraction",
                CompositionBasis::MoleFraction => "mole-fraction",
                CompositionBasis::VolumeFraction => "volume-fraction",
            })],
        ),
        QuantityKind::Mass => form("mass", Vec::new()),
        QuantityKind::Amount => form("amount", Vec::new()),
        QuantityKind::MolarMass => form("molar-mass", Vec::new()),
        QuantityKind::MassConcentration => form("mass-concentration", Vec::new()),
        QuantityKind::AmountConcentration => form("amount-concentration", Vec::new()),
        QuantityKind::Entropy => form("entropy", Vec::new()),
        QuantityKind::HeatCapacity => form("heat-capacity", Vec::new()),
        QuantityKind::AcousticPressure => form("acoustic-pressure", Vec::new()),
        QuantityKind::AcousticPower => form("acoustic-power", Vec::new()),
    }
}

fn write_angle_domain(domain: AngleDomain) -> Node {
    sym(match domain {
        AngleDomain::Mechanical => "mechanical",
        AngleDomain::Electrical => "electrical",
    })
}

pub(super) fn write_shape(shape: TerminalShape) -> Node {
    match shape {
        TerminalShape::Scalar => form("scalar", Vec::new()),
        TerminalShape::Vector { components } => form("vector", vec![u64_node(components.get())]),
        TerminalShape::Tensor { rows, columns } => form(
            "tensor",
            vec![u64_node(rows.get()), u64_node(columns.get())],
        ),
        TerminalShape::FieldTrace { components } => {
            form("field-trace", vec![u64_node(components.get())])
        }
    }
}

fn write_port(port: &PortSpec) -> Node {
    form(
        "port",
        vec![
            string_node(port.id.canonical_key()),
            string_node(port.owner.canonical_key()),
            string_node(port.effort.canonical_key()),
            string_node(port.flow.canonical_key()),
            sym(match port.energy_role {
                PortEnergyRole::IntoSubsystem => "into-subsystem",
                PortEnergyRole::OutOfSubsystem => "out-of-subsystem",
            }),
        ],
    )
}

fn write_relation(relation: &RelationSpec) -> Node {
    let mode = match &relation.mode {
        RelationMode::Algebraic { solve_policy } => form(
            "algebraic",
            solve_policy
                .iter()
                .map(|policy| {
                    reference_node(
                        policy.namespace(),
                        policy.schema_version(),
                        policy.semantic_digest(),
                    )
                })
                .collect(),
        ),
        RelationMode::Stateful { state_slot } => {
            form("stateful", vec![string_node(state_slot.canonical_key())])
        }
    };
    form(
        "relation",
        vec![
            string_node(relation.id.canonical_key()),
            string_node(relation.source.canonical_key()),
            string_node(relation.target.canonical_key()),
            mode,
        ],
    )
}

fn write_material(material: &MaterialBinding) -> Node {
    let target = match &material.target {
        MaterialTarget::Body(id) => form("body", vec![string_node(id.canonical_key())]),
        MaterialTarget::SurfacePatch(id) => {
            form("surface-patch", vec![string_node(id.canonical_key())])
        }
    };
    form(
        "material",
        vec![
            target,
            reference_node(
                material.material.namespace(),
                material.material.schema_version(),
                material.material.semantic_digest(),
            ),
        ],
    )
}

fn write_interface(interface: &InterfaceBinding) -> Node {
    form(
        "interface",
        vec![
            string_node(interface.id.canonical_key()),
            string_node(interface.negative.canonical_key()),
            string_node(interface.positive.canonical_key()),
            reference_node(
                interface.interface.namespace(),
                interface.interface.schema_version(),
                interface.interface.semantic_digest(),
            ),
            sym(match interface.orientation {
                InterfaceOrientation::Aligned => "aligned",
                InterfaceOrientation::Opposed => "opposed",
            }),
        ],
    )
}

pub(super) fn write_frame(frame: &FrameBinding) -> Node {
    form(
        "frame",
        vec![
            string_node(frame.canonical_key()),
            sym(match frame.orientation() {
                OrientationParity::Preserving => "preserving",
                OrientationParity::Reversing => "reversing",
            }),
        ],
    )
}

fn write_state_contract(contract: &StateSlotContract) -> Node {
    form(
        "state-contract",
        vec![
            string_node(contract.id.canonical_key()),
            string_node(contract.owner.canonical_key()),
            write_quantity(contract.quantity),
            write_shape(contract.shape),
            string_node(contract.clock.canonical_key()),
            write_frame(&contract.frame),
        ],
    )
}

fn write_condition(condition: &ConditionBinding) -> Node {
    form(
        "condition",
        vec![
            write_condition_target(&condition.target),
            write_quantity(condition.quantity),
            write_shape(condition.shape),
            string_node(condition.clock.canonical_key()),
            write_frame(&condition.frame),
            write_condition_source(&condition.source),
        ],
    )
}

fn write_condition_target(target: &ConditionTarget) -> Node {
    match target {
        ConditionTarget::Initial(id) => form("initial", vec![string_node(id.canonical_key())]),
        ConditionTarget::Boundary(id) => form("boundary", vec![string_node(id.canonical_key())]),
    }
}

fn write_condition_source(source: &ConditionSource) -> Node {
    match source {
        ConditionSource::Fixed(value) => form(
            "fixed",
            vec![reference_node(
                value.namespace(),
                value.schema_version(),
                value.semantic_digest(),
            )],
        ),
        ConditionSource::History {
            history,
            continuity,
        } => form(
            "history",
            vec![
                reference_node(
                    history.namespace(),
                    history.schema_version(),
                    history.semantic_digest(),
                ),
                match continuity {
                    HistoryContinuity::Continuous => form("continuous", Vec::new()),
                    HistoryContinuity::ResetAtEvents { events } => form(
                        "reset-at-events",
                        events
                            .iter()
                            .map(|event| string_node(event.canonical_key()))
                            .collect(),
                    ),
                },
            ],
        ),
        ConditionSource::Distribution(distribution) => form(
            "distribution",
            vec![reference_node(
                distribution.namespace(),
                distribution.schema_version(),
                distribution.semantic_digest(),
            )],
        ),
    }
}

fn write_motion(motion: &MotionBinding) -> Node {
    let body_motion = match &motion.motion {
        BodyMotion::Static => form("static", Vec::new()),
        BodyMotion::Prescribed { path } => form(
            "prescribed",
            vec![reference_node(
                path.namespace(),
                path.schema_version(),
                path.semantic_digest(),
            )],
        ),
    };
    form(
        "motion",
        vec![
            string_node(motion.body.canonical_key()),
            string_node(motion.clock.canonical_key()),
            write_frame(&motion.reference_frame),
            body_motion,
        ],
    )
}

fn write_event(event: &EventSpec) -> Node {
    form(
        "event",
        vec![
            string_node(event.id.canonical_key()),
            string_node(event.clock.canonical_key()),
            reference_node(
                event.guard.namespace(),
                event.guard.schema_version(),
                event.guard.semantic_digest(),
            ),
            sym(match event.orientation {
                GuardOrientation::NegativeToPositive => "negative-to-positive",
                GuardOrientation::PositiveToNegative => "positive-to-negative",
                GuardOrientation::Bidirectional => "bidirectional",
            }),
            write_crossing(&event.crossing),
            section(
                "dependencies",
                event
                    .dependencies
                    .iter()
                    .map(|dependency| match dependency {
                        EventDependency::State(id) => {
                            form("state", vec![string_node(id.canonical_key())])
                        }
                        EventDependency::Terminal(id) => {
                            form("terminal", vec![string_node(id.canonical_key())])
                        }
                    })
                    .collect(),
            ),
            write_reset(&event.reset),
            write_event_order(&event.order),
        ],
    )
}

fn write_crossing(crossing: &CrossingSemantics) -> Node {
    match crossing {
        CrossingSemantics::Transverse(witness) => form(
            "transverse",
            vec![reference_node(
                witness.namespace(),
                witness.schema_version(),
                witness.semantic_digest(),
            )],
        ),
        CrossingSemantics::Grazing(witness) => form(
            "grazing",
            vec![reference_node(
                witness.namespace(),
                witness.schema_version(),
                witness.semantic_digest(),
            )],
        ),
        CrossingSemantics::Unknown(no_claim) => form(
            "unknown",
            vec![reference_node(
                no_claim.namespace(),
                no_claim.schema_version(),
                no_claim.semantic_digest(),
            )],
        ),
    }
}

fn write_reset(reset: &ResetSemantics) -> Node {
    match reset {
        ResetSemantics::Deterministic { map, writes } => form(
            "deterministic",
            vec![
                reference_node(map.namespace(), map.schema_version(), map.semantic_digest()),
                write_reset_writes(writes),
            ],
        ),
        ResetSemantics::SetValued {
            relation,
            outcomes,
            writes,
        } => form(
            "set-valued",
            vec![
                reference_node(
                    relation.namespace(),
                    relation.schema_version(),
                    relation.semantic_digest(),
                ),
                reference_node(
                    outcomes.namespace(),
                    outcomes.schema_version(),
                    outcomes.semantic_digest(),
                ),
                write_reset_writes(writes),
            ],
        ),
        ResetSemantics::Terminal { relation } => form(
            "terminal",
            vec![reference_node(
                relation.namespace(),
                relation.schema_version(),
                relation.semantic_digest(),
            )],
        ),
        ResetSemantics::Unknown { no_claim } => form(
            "unknown",
            vec![reference_node(
                no_claim.namespace(),
                no_claim.schema_version(),
                no_claim.semantic_digest(),
            )],
        ),
    }
}

fn write_reset_writes(writes: &[StateSlotId]) -> Node {
    section(
        "writes",
        writes
            .iter()
            .map(|id| string_node(id.canonical_key()))
            .collect(),
    )
}

fn write_event_order(order: &EventOrder) -> Node {
    match order {
        EventOrder::TotalPriority { microstep } => {
            form("total-priority", vec![u64_node(u64::from(*microstep))])
        }
        EventOrder::SetValued { group } => form(
            "set-valued",
            vec![reference_node(
                group.namespace(),
                group.schema_version(),
                group.semantic_digest(),
            )],
        ),
    }
}

fn write_tolerance(tolerance: &ToleranceSpec) -> Node {
    form(
        "tolerance",
        vec![
            string_node(tolerance.id.canonical_key()),
            write_tolerance_target(&tolerance.target),
            reference_node(
                tolerance.parameter.namespace(),
                tolerance.parameter.schema_version(),
                tolerance.parameter.semantic_digest(),
            ),
            write_quantity(tolerance.quantity),
            write_shape(tolerance.shape),
            write_tolerance_semantics(&tolerance.semantics),
        ],
    )
}

fn write_tolerance_target(target: &ToleranceTarget) -> Node {
    match target {
        ToleranceTarget::Subsystem(id) => form("subsystem", vec![string_node(id.canonical_key())]),
        ToleranceTarget::Element(id) => form("element", vec![write_machine_element(id)]),
    }
}

pub(super) fn write_machine_element(element: &MachineElementId) -> Node {
    let (head, key) = match element {
        MachineElementId::Body(id) => ("body", id.canonical_key()),
        MachineElementId::SurfacePatch(id) => ("surface-patch", id.canonical_key()),
        MachineElementId::ContactFeature(id) => ("contact-feature", id.canonical_key()),
        MachineElementId::Terminal(id) => ("terminal", id.canonical_key()),
        MachineElementId::Port(id) => ("port", id.canonical_key()),
        MachineElementId::StateSlot(id) => ("state-slot", id.canonical_key()),
    };
    form(head, vec![string_node(key)])
}

fn write_tolerance_semantics(semantics: &ToleranceSemantics) -> Node {
    match semantics {
        ToleranceSemantics::Bounded { minus, plus, law } => form(
            "bounded",
            vec![
                float_node(minus.get()),
                float_node(plus.get()),
                reference_node(law.namespace(), law.schema_version(), law.semantic_digest()),
            ],
        ),
        ToleranceSemantics::Random {
            scale,
            law,
            marginal,
        } => form(
            "random",
            vec![
                float_node(scale.get()),
                reference_node(law.namespace(), law.schema_version(), law.semantic_digest()),
                reference_node(
                    marginal.namespace(),
                    marginal.schema_version(),
                    marginal.semantic_digest(),
                ),
            ],
        ),
    }
}

fn write_dependence(dependence: &DependenceSpec) -> Node {
    let members = dependence
        .members
        .iter()
        .map(|member| match member {
            DependenceMember::Condition(target) => {
                form("condition", vec![write_condition_target(target)])
            }
            DependenceMember::Tolerance(id) => {
                form("tolerance", vec![string_node(id.canonical_key())])
            }
        })
        .collect();
    let model = match &dependence.model {
        DependenceModel::Independent => form("independent", Vec::new()),
        DependenceModel::Correlated(correlation) => form(
            "correlated",
            vec![reference_node(
                correlation.namespace(),
                correlation.schema_version(),
                correlation.semantic_digest(),
            )],
        ),
    };
    form("dependence", vec![section("members", members), model])
}
