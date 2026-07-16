//! Canonical FrankenScript AST transport for admitted Machine assurance.
//!
//! This module owns syntax only. It decodes an explicit, versioned AST into an
//! authority-free [`MachineAssuranceDraft`] bound to exact graph and behavior
//! identities. Serialized V&V receipt rows are commitments, never authority:
//! admission still requires caller-supplied [`AdmittedVvCase`] values, verifies
//! them through the existing assurance boundary, and only then compares the
//! derived bindings with the declared rows.

use core::cmp::Ordering;
use core::fmt;

use fs_blake3::ContentHash;
use fs_blake3::identity::StrongIdentity;
use fs_evidence::vv::{
    AdmittedVvCase, ArtifactId, ArtifactKind, ArtifactRef, AssumptionId, MAX_VV_ID_BYTES, QoiId,
    UnitId,
};

use crate::VersionedProgram;
use crate::ast::{Node, NodeKind, Span};

use super::assurance::{
    AccountingBoundaryRef, AccountingEntry, AccountingIntervalRef, AccountingOrientation,
    AccountingPolicyRef, AccountingRole, AccountingTarget, AccountingWindow, AccountingWindowId,
    AdmittedMachineAssurance, BalanceKind, BalanceLawRef, CalibrationRef, ContextBinding,
    ContextQoiKey, CostErrorModelRef, DecisionBudgetRef, EscalationAction, EscalationSpec,
    EscalationTriggerRef, ExperimentId, ExperimentSpec, FalsifierRef, FaultContainmentRef,
    FaultCoverage, FaultId, FaultInjectionRef, FaultModelRef, FaultSpec, FidelityPolicy,
    FidelityRung, FidelityRungId, FixedReplayRef, HazardId, HazardSpec, LossOwnershipRef,
    MAX_MACHINE_ASSURANCE_ACCOUNTING_WINDOWS, MAX_MACHINE_ASSURANCE_CONTEXTS,
    MAX_MACHINE_ASSURANCE_EXPERIMENTS, MAX_MACHINE_ASSURANCE_FAULTS,
    MAX_MACHINE_ASSURANCE_FIDELITY_RUNGS, MAX_MACHINE_ASSURANCE_HAZARDS,
    MAX_MACHINE_ASSURANCE_NESTED_REFERENCES, MAX_MACHINE_ASSURANCE_SENSORS, MachineAssuranceDraft,
    MachineAssuranceRefusal, MachineScope, ModelCrosswalkRef, NoClaimRef, ObservationTarget,
    ObservationTiming, OperatingEnvelopeRef, QoiBinding, QoiDefinitionRef, QoiInput, QoiTarget,
    SafetyCaseRef, SafetyRequirementRef, SamplingBridgeRef, SensorExposure, SensorId,
    SensorInstrumentBinding, SensorModelRef, SensorSpec, StateTransferRef, UnitQuantityBridgeRef,
    ValidityDomainRef, VvCaseBinding,
};
use super::codec::{
    MachineGraphCodecError, MachineGraphCodecRule, digest_hex, form, parse_digest, parse_frame,
    parse_id, parse_quantity, parse_reference, parse_shape, parse_u64, recognized_form_items,
    reference_node, reserved_vec, section, string, string_node, sym, symbol, u64_node,
    validation_path, write_frame, write_machine_element, write_quantity, write_shape,
};
use super::semantics::{AdmittedMachineBehavior, MachineBehaviorIdV1};
use super::{
    AdmittedMachineGraph, BodyId, ClockId, ContactFeatureId, InterfaceId, MachineElementId,
    MachineGraphIdV1, ModelRef, PortId, RelationId, StateSlotId, SubsystemId, SurfacePatchId,
    TerminalId,
};

/// Version of the canonical Machine-assurance FrankenScript form.
pub const MACHINE_ASSURANCE_AST_SCHEMA_VERSION_V1: u32 = 1;
/// Root symbol for the canonical Machine-assurance FrankenScript form.
pub const MACHINE_ASSURANCE_AST_HEAD_V1: &str = "machine-assurance-v1";
/// Maximum generic AST nodes inspected by one assurance decode.
///
/// The bound is derived from every public top-level and aggregate nested cap,
/// with a conservative per-row envelope. It intentionally exceeds the graph
/// codec's smaller historical generic-AST limit.
pub const MAX_MACHINE_ASSURANCE_AST_NODES: usize = 5_832_832;
/// Maximum aggregate text bytes inspected by one assurance decode.
pub const MAX_MACHINE_ASSURANCE_AST_TEXT_BYTES: usize = 256 * 1_024 * 1_024;
/// Maximum bytes accepted in any one assurance AST text atom.
pub const MAX_MACHINE_ASSURANCE_AST_ATOM_BYTES: usize = MAX_VV_ID_BYTES;
/// Maximum serialized V&V binding commitments in one assurance artifact.
pub const MAX_MACHINE_ASSURANCE_VV_CASE_RECEIPTS: usize = MAX_MACHINE_ASSURANCE_CONTEXTS;

/// Closed syntax/resource refusal vocabulary for the assurance codec.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u8)]
pub enum MachineAssuranceCodecRule {
    /// The caller supplied a forged or otherwise invalid generic AST.
    InvalidAst = 1,
    /// A list had the wrong head, arity, or atom kind.
    UnexpectedForm = 2,
    /// A closed enum/tag token was unknown.
    UnknownTag = 3,
    /// A numeric atom was noncanonical or outside its declared domain.
    InvalidNumber = 4,
    /// A Machine or fs-evidence identifier was invalid.
    InvalidIdentifier = 5,
    /// An identity, artifact reference, or opaque external reference was malformed.
    InvalidReference = 6,
    /// A public collection or aggregate nested-reference bound was exceeded.
    ResourceLimit = 7,
}

impl MachineAssuranceCodecRule {
    /// Stable structured diagnostic code.
    #[must_use]
    pub const fn code(self) -> &'static str {
        match self {
            Self::InvalidAst => "MachineAssuranceCodecInvalidAst",
            Self::UnexpectedForm => "MachineAssuranceCodecUnexpectedForm",
            Self::UnknownTag => "MachineAssuranceCodecUnknownTag",
            Self::InvalidNumber => "MachineAssuranceCodecInvalidNumber",
            Self::InvalidIdentifier => "MachineAssuranceCodecInvalidIdentifier",
            Self::InvalidReference => "MachineAssuranceCodecInvalidReference",
            Self::ResourceLimit => "MachineAssuranceCodecResourceLimit",
        }
    }
}

/// One bounded, path-addressed Machine-assurance syntax refusal.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MachineAssuranceCodecError {
    rule: MachineAssuranceCodecRule,
    span: Span,
    path: Box<str>,
    detail: Box<str>,
    hint: Box<str>,
}

impl MachineAssuranceCodecError {
    /// Closed refusal rule.
    #[must_use]
    pub const fn rule(&self) -> MachineAssuranceCodecRule {
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

    /// Deterministic structural path within the assurance form.
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
        let rule = match error.rule() {
            MachineGraphCodecRule::InvalidAst => MachineAssuranceCodecRule::InvalidAst,
            MachineGraphCodecRule::UnexpectedForm => MachineAssuranceCodecRule::UnexpectedForm,
            MachineGraphCodecRule::UnknownTag => MachineAssuranceCodecRule::UnknownTag,
            MachineGraphCodecRule::InvalidNumber => MachineAssuranceCodecRule::InvalidNumber,
            MachineGraphCodecRule::InvalidIdentifier => {
                MachineAssuranceCodecRule::InvalidIdentifier
            }
            MachineGraphCodecRule::InvalidReference => MachineAssuranceCodecRule::InvalidReference,
            MachineGraphCodecRule::ResourceLimit => MachineAssuranceCodecRule::ResourceLimit,
        };
        Self {
            rule,
            span: error.span(),
            path: error.path().to_string().into_boxed_str(),
            detail: error.detail().to_string().into_boxed_str(),
            hint: error.hint().to_string().into_boxed_str(),
        }
    }
}

impl fmt::Display for MachineAssuranceCodecError {
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

impl std::error::Error for MachineAssuranceCodecError {}

impl From<MachineGraphCodecError> for MachineAssuranceCodecError {
    fn from(value: MachineGraphCodecError) -> Self {
        Self::from_graph(value)
    }
}

/// Decoded assurance syntax plus explicit base identities and non-authoritative
/// receipt commitments.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DecodedMachineAssuranceV1 {
    base_graph: MachineGraphIdV1,
    base_behavior: MachineBehaviorIdV1,
    declared_vv_cases: Vec<VvCaseBinding>,
    draft: MachineAssuranceDraft,
}

impl DecodedMachineAssuranceV1 {
    /// Exact graph identity declared by the syntax.
    #[must_use]
    pub const fn base_graph(&self) -> MachineGraphIdV1 {
        self.base_graph
    }

    /// Exact behavior identity declared by the syntax.
    #[must_use]
    pub const fn base_behavior(&self) -> MachineBehaviorIdV1 {
        self.base_behavior
    }

    /// Non-authoritative receipt commitments declared by the syntax.
    #[must_use]
    pub fn declared_vv_cases(&self) -> &[VvCaseBinding] {
        &self.declared_vv_cases
    }

    /// Borrow the authority-free decoded overlay.
    #[must_use]
    pub const fn draft(&self) -> &MachineAssuranceDraft {
        &self.draft
    }

    /// Recover every decoded part without publishing semantic authority.
    #[must_use]
    pub fn into_parts(
        self,
    ) -> (
        MachineGraphIdV1,
        MachineBehaviorIdV1,
        Vec<VvCaseBinding>,
        MachineAssuranceDraft,
    ) {
        (
            self.base_graph,
            self.base_behavior,
            self.declared_vv_cases,
            self.draft,
        )
    }

    /// Check exact base identities, perform real semantic admission with
    /// caller-supplied V&V authority, then compare derived receipt bindings.
    ///
    /// # Errors
    /// Refuses either base mismatch before semantic inspection, preserves the
    /// complete assurance refusal, and rejects missing, duplicate, changed, or
    /// extra serialized receipt commitments after successful admission.
    pub fn admit_against(
        self,
        graph: &AdmittedMachineGraph,
        behavior: &AdmittedMachineBehavior,
        vv_cases: &[AdmittedVvCase],
    ) -> Result<AdmittedMachineAssurance, MachineAssuranceAstAdmissionError> {
        if self.base_graph != graph.identity() {
            return Err(MachineAssuranceAstAdmissionError::BaseGraphMismatch {
                declared: self.base_graph,
                provided: graph.identity(),
            });
        }
        if self.base_behavior != behavior.identity() {
            return Err(MachineAssuranceAstAdmissionError::BaseBehaviorMismatch {
                declared: self.base_behavior,
                provided: behavior.identity(),
            });
        }
        let admitted = self
            .draft
            .admit_against(graph, behavior, vv_cases)
            .map_err(MachineAssuranceAstAdmissionError::Assurance)?;
        compare_vv_case_bindings(&self.declared_vv_cases, admitted.vv_cases())?;
        Ok(admitted)
    }
}

/// Distinguishes syntax, identity binding, evidence commitment, and semantic
/// assurance refusals.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MachineAssuranceAstAdmissionError {
    /// The generic AST could not be decoded as the v1 assurance grammar.
    Codec(MachineAssuranceCodecError),
    /// The syntax names a different admitted graph than the caller supplied.
    BaseGraphMismatch {
        /// Identity embedded in the assurance syntax.
        declared: MachineGraphIdV1,
        /// Identity of the graph offered for admission.
        provided: MachineGraphIdV1,
    },
    /// The syntax names a different admitted behavior than the caller supplied.
    BaseBehaviorMismatch {
        /// Identity embedded in the assurance syntax.
        declared: MachineBehaviorIdV1,
        /// Identity of the behavior offered for admission.
        provided: MachineBehaviorIdV1,
    },
    /// Serialized receipt commitments disagree with verified derived bindings.
    VvCaseBindingMismatch {
        /// Deterministic mismatch explanation; the syntax grants no authority.
        detail: Box<str>,
    },
    /// The authority-free draft failed the existing assurance admission.
    Assurance(MachineAssuranceRefusal),
}

impl MachineAssuranceAstAdmissionError {
    /// Stable top-level refusal code without collapsing authority boundaries.
    #[must_use]
    pub const fn code(&self) -> &'static str {
        match self {
            Self::Codec(error) => error.code(),
            Self::BaseGraphMismatch { .. } => "MachineAssuranceBaseGraphMismatch",
            Self::BaseBehaviorMismatch { .. } => "MachineAssuranceBaseBehaviorMismatch",
            Self::VvCaseBindingMismatch { .. } => "MachineAssuranceVvCaseBindingMismatch",
            Self::Assurance(_) => "MachineAssuranceRefused",
        }
    }
}

impl fmt::Display for MachineAssuranceAstAdmissionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Codec(error) => error.fmt(f),
            Self::BaseGraphMismatch { declared, provided } => write!(
                f,
                "MachineAssuranceBaseGraphMismatch: syntax binds {}, provided graph is {}",
                digest_hex(*declared.as_bytes()),
                digest_hex(*provided.as_bytes())
            ),
            Self::BaseBehaviorMismatch { declared, provided } => write!(
                f,
                "MachineAssuranceBaseBehaviorMismatch: syntax binds {}, provided behavior is {}",
                digest_hex(*declared.as_bytes()),
                digest_hex(*provided.as_bytes())
            ),
            Self::VvCaseBindingMismatch { detail } => {
                write!(f, "MachineAssuranceVvCaseBindingMismatch: {detail}")
            }
            Self::Assurance(error) => error.fmt(f),
        }
    }
}

impl std::error::Error for MachineAssuranceAstAdmissionError {}

impl From<MachineAssuranceCodecError> for MachineAssuranceAstAdmissionError {
    fn from(value: MachineAssuranceCodecError) -> Self {
        Self::Codec(value)
    }
}

/// Decode one canonical v1 Machine-assurance AST into an authority-free draft.
///
/// # Errors
/// Refuses malformed ASTs, unknown closed tags, invalid identifiers or
/// references, noncanonical numbers, and oversized collections.
pub fn parse_machine_assurance_v1(
    node: &Node,
) -> Result<DecodedMachineAssuranceV1, MachineAssuranceCodecError> {
    preflight_assurance_ast(node)?;
    preflight_assurance_caps(node)?;
    node.validate()
        .map_err(|error| MachineAssuranceCodecError {
            rule: MachineAssuranceCodecRule::InvalidAst,
            span: error.span,
            path: validation_path(&error.detail),
            detail: error.to_string().into_boxed_str(),
            hint: "repair the generic FrankenScript AST before Machine assurance decoding"
                .to_string()
                .into_boxed_str(),
        })?;

    let root = exact_form(node, MACHINE_ASSURANCE_AST_HEAD_V1, 10, "$")?;
    let base_graph = parse_graph_identity(&root[0], "$[1]")?;
    let base_behavior = parse_behavior_identity(&root[1], "$[2]")?;
    let vv_rows = section_items(
        &root[2],
        "vv-case-receipts",
        MAX_MACHINE_ASSURANCE_VV_CASE_RECEIPTS,
        "$[3]",
    )?;
    let sensor_rows = section_items(&root[3], "sensors", MAX_MACHINE_ASSURANCE_SENSORS, "$[4]")?;
    let experiment_rows = section_items(
        &root[4],
        "experiments",
        MAX_MACHINE_ASSURANCE_EXPERIMENTS,
        "$[5]",
    )?;
    let context_rows = section_items(&root[5], "contexts", MAX_MACHINE_ASSURANCE_CONTEXTS, "$[6]")?;
    let hazard_rows = section_items(&root[6], "hazards", MAX_MACHINE_ASSURANCE_HAZARDS, "$[7]")?;
    let fault_rows = section_items(&root[7], "faults", MAX_MACHINE_ASSURANCE_FAULTS, "$[8]")?;
    let accounting_rows = section_items(
        &root[8],
        "accounting-windows",
        MAX_MACHINE_ASSURANCE_ACCOUNTING_WINDOWS,
        "$[9]",
    )?;

    Ok(DecodedMachineAssuranceV1 {
        base_graph,
        base_behavior,
        declared_vv_cases: parse_rows(vv_rows, &root[2], "$[3]", parse_vv_case_binding)?,
        draft: MachineAssuranceDraft {
            sensors: parse_rows(sensor_rows, &root[3], "$[4]", parse_sensor)?,
            experiments: parse_rows(experiment_rows, &root[4], "$[5]", parse_experiment)?,
            contexts: parse_rows(context_rows, &root[5], "$[6]", parse_context)?,
            hazards: parse_rows(hazard_rows, &root[6], "$[7]", parse_hazard)?,
            faults: parse_rows(fault_rows, &root[7], "$[8]", parse_fault)?,
            accounting_windows: parse_rows(
                accounting_rows,
                &root[8],
                "$[9]",
                parse_accounting_window,
            )?,
            fidelity: parse_fidelity(&root[9], "$[10]")?,
        },
    })
}

/// Decode a version-enforced FrankenScript program body as Machine assurance.
///
/// # Errors
/// Returns the same bounded syntax/resource refusals as
/// [`parse_machine_assurance_v1`].
pub fn parse_machine_assurance_program_v1(
    program: &VersionedProgram,
) -> Result<DecodedMachineAssuranceV1, MachineAssuranceCodecError> {
    parse_machine_assurance_v1(program.program())
}

/// Decode and semantically admit one Machine-assurance AST.
///
/// # Errors
/// Preserves syntax, both identity bindings, V&V commitment comparison, and
/// semantic refusal boundaries.
pub fn admit_machine_assurance_ast_v1(
    node: &Node,
    graph: &AdmittedMachineGraph,
    behavior: &AdmittedMachineBehavior,
    vv_cases: &[AdmittedVvCase],
) -> Result<AdmittedMachineAssurance, MachineAssuranceAstAdmissionError> {
    parse_machine_assurance_v1(node)?.admit_against(graph, behavior, vv_cases)
}

/// Encode one admitted assurance overlay as canonical v1 Machine-assurance AST.
///
/// # Errors
/// Returns a structured internal-boundary error if an admitted value cannot be
/// represented by the declared grammar.
pub fn write_machine_assurance_v1(
    assurance: &AdmittedMachineAssurance,
) -> Result<Node, MachineAssuranceCodecError> {
    let root = form(
        MACHINE_ASSURANCE_AST_HEAD_V1,
        vec![
            form(
                "base-graph",
                vec![string_node(&digest_hex(*assurance.base_graph().as_bytes()))],
            ),
            form(
                "base-behavior",
                vec![string_node(&digest_hex(
                    *assurance.base_behavior().as_bytes(),
                ))],
            ),
            section(
                "vv-case-receipts",
                assurance
                    .vv_cases()
                    .iter()
                    .map(write_vv_case_binding)
                    .collect(),
            ),
            section(
                "sensors",
                assurance.sensors().iter().map(write_sensor).collect(),
            ),
            section(
                "experiments",
                assurance
                    .experiments()
                    .iter()
                    .map(write_experiment)
                    .collect(),
            ),
            section(
                "contexts",
                assurance.contexts().iter().map(write_context).collect(),
            ),
            section(
                "hazards",
                assurance.hazards().iter().map(write_hazard).collect(),
            ),
            section(
                "faults",
                assurance.faults().iter().map(write_fault).collect(),
            ),
            section(
                "accounting-windows",
                assurance
                    .accounting_windows()
                    .iter()
                    .map(write_accounting_window)
                    .collect(),
            ),
            write_fidelity(assurance.fidelity()),
        ],
    );
    preflight_assurance_ast(&root)?;
    preflight_assurance_caps(&root)?;
    root.validate()
        .map_err(|error| MachineAssuranceCodecError {
            rule: MachineAssuranceCodecRule::InvalidAst,
            span: error.span,
            path: validation_path(&error.detail),
            detail: format!("canonical Machine assurance writer produced an invalid AST: {error}")
                .into_boxed_str(),
            hint: "treat this as a codec implementation defect"
                .to_string()
                .into_boxed_str(),
        })?;
    Ok(root)
}

/// Encode one admitted assurance overlay into a current-version artifact.
///
/// # Errors
/// Refuses an internal representation defect rather than emitting invalid or
/// unversioned assurance syntax.
pub fn write_machine_assurance_program_v1(
    assurance: &AdmittedMachineAssurance,
) -> Result<VersionedProgram, MachineAssuranceCodecError> {
    let node = write_machine_assurance_v1(assurance)?;
    VersionedProgram::try_current(node).map_err(|error| MachineAssuranceCodecError {
        rule: MachineAssuranceCodecRule::InvalidAst,
        span: error.span,
        path: validation_path(&error.detail),
        detail: format!("canonical Machine assurance could not enter the IR envelope: {error}")
            .into_boxed_str(),
        hint: "treat this as a codec implementation defect"
            .to_string()
            .into_boxed_str(),
    })
}

fn preflight_assurance_ast(node: &Node) -> Result<(), MachineAssuranceCodecError> {
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
        if visited > MAX_MACHINE_ASSURANCE_AST_NODES {
            return Err(resource_error(
                current,
                "$",
                &format!(
                    "AST node count exceeds cap {MAX_MACHINE_ASSURANCE_AST_NODES} before semantic decoding"
                ),
            ));
        }
        let added_text = match &current.kind {
            NodeKind::Qty { text, .. }
            | NodeKind::Str(text)
            | NodeKind::Symbol(text)
            | NodeKind::Keyword(text) => {
                if text.len() > MAX_MACHINE_ASSURANCE_AST_ATOM_BYTES {
                    return Err(resource_error(
                        current,
                        "$",
                        &format!(
                            "AST text atom bytes {} exceed cap {MAX_MACHINE_ASSURANCE_AST_ATOM_BYTES}",
                            text.len()
                        ),
                    ));
                }
                text.len()
            }
            NodeKind::Int(_)
            | NodeKind::Float(_)
            | NodeKind::Count { .. }
            | NodeKind::Seed(_)
            | NodeKind::List(_) => 0,
        };
        text_bytes = text_bytes
            .checked_add(added_text)
            .ok_or_else(|| resource_error(current, "$", "AST text byte count overflowed usize"))?;
        if text_bytes > MAX_MACHINE_ASSURANCE_AST_TEXT_BYTES {
            return Err(resource_error(
                current,
                "$",
                &format!(
                    "AST text bytes exceed cap {MAX_MACHINE_ASSURANCE_AST_TEXT_BYTES} before semantic decoding"
                ),
            ));
        }
        if let NodeKind::List(items) = &current.kind {
            let projected = visited
                .checked_add(stack.len())
                .and_then(|count| count.checked_add(items.len()))
                .ok_or_else(|| resource_error(current, "$", "AST node count overflowed usize"))?;
            if projected > MAX_MACHINE_ASSURANCE_AST_NODES {
                return Err(resource_error(
                    current,
                    "$",
                    &format!(
                        "AST node count exceeds cap {MAX_MACHINE_ASSURANCE_AST_NODES} before semantic decoding"
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

fn preflight_assurance_caps(node: &Node) -> Result<(), MachineAssuranceCodecError> {
    let Some(root) = recognized_form_items(node, MACHINE_ASSURANCE_AST_HEAD_V1) else {
        return Ok(());
    };
    if root.len() < 10 {
        return Ok(());
    }
    cap_recognized_section(
        &root[2],
        "vv-case-receipts",
        MAX_MACHINE_ASSURANCE_VV_CASE_RECEIPTS,
        "$[3]",
    )?;
    cap_recognized_section(&root[3], "sensors", MAX_MACHINE_ASSURANCE_SENSORS, "$[4]")?;
    let experiments = cap_recognized_section(
        &root[4],
        "experiments",
        MAX_MACHINE_ASSURANCE_EXPERIMENTS,
        "$[5]",
    )?;
    let contexts =
        cap_recognized_section(&root[5], "contexts", MAX_MACHINE_ASSURANCE_CONTEXTS, "$[6]")?;
    let hazards =
        cap_recognized_section(&root[6], "hazards", MAX_MACHINE_ASSURANCE_HAZARDS, "$[7]")?;
    let faults = cap_recognized_section(&root[7], "faults", MAX_MACHINE_ASSURANCE_FAULTS, "$[8]")?;
    let accounting = cap_recognized_section(
        &root[8],
        "accounting-windows",
        MAX_MACHINE_ASSURANCE_ACCOUNTING_WINDOWS,
        "$[9]",
    )?;

    let mut nested = 0_usize;
    if let Some(rows) = experiments {
        for (index, row) in rows.iter().enumerate() {
            if let Some(args) = recognized_form_items(row, "experiment") {
                if args.len() >= 5 {
                    add_nested_section(
                        &mut nested,
                        &args[3],
                        "instruments",
                        &format!("$[5][{}][4]", index + 1),
                    )?;
                    add_nested_section(
                        &mut nested,
                        &args[4],
                        "qois",
                        &format!("$[5][{}][5]", index + 1),
                    )?;
                }
            }
        }
    }
    if let Some(rows) = contexts {
        for (index, row) in rows.iter().enumerate() {
            if let Some(args) = recognized_form_items(row, "context") {
                if args.len() >= 3 {
                    let path = format!("$[6][{}][3]", index + 1);
                    if let Some(qois) = add_nested_section(&mut nested, &args[2], "qois", &path)? {
                        for (qoi_index, qoi) in qois.iter().enumerate() {
                            if let Some(qoi_args) = recognized_form_items(qoi, "qoi") {
                                if qoi_args.len() >= 2 {
                                    add_nested_section(
                                        &mut nested,
                                        &qoi_args[1],
                                        "inputs",
                                        &format!("{path}[{}][2]", qoi_index + 1),
                                    )?;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    if let Some(rows) = hazards {
        for (index, row) in rows.iter().enumerate() {
            if let Some(args) = recognized_form_items(row, "hazard") {
                if args.len() >= 7 {
                    add_nested_section(
                        &mut nested,
                        &args[2],
                        "scope",
                        &format!("$[7][{}][3]", index + 1),
                    )?;
                    add_nested_section(
                        &mut nested,
                        &args[6],
                        "assumptions",
                        &format!("$[7][{}][7]", index + 1),
                    )?;
                }
            }
        }
    }
    if let Some(rows) = faults {
        for (index, row) in rows.iter().enumerate() {
            if let Some(args) = recognized_form_items(row, "fault") {
                if args.len() >= 3 {
                    add_nested_section(
                        &mut nested,
                        &args[1],
                        "affected",
                        &format!("$[8][{}][2]", index + 1),
                    )?;
                    add_nested_section(
                        &mut nested,
                        &args[2],
                        "hazards",
                        &format!("$[8][{}][3]", index + 1),
                    )?;
                }
            }
        }
    }
    if let Some(rows) = accounting {
        for (index, row) in rows.iter().enumerate() {
            if let Some(args) = recognized_form_items(row, "accounting-window") {
                if args.len() >= 8 {
                    add_nested_section(
                        &mut nested,
                        &args[7],
                        "entries",
                        &format!("$[9][{}][8]", index + 1),
                    )?;
                }
            }
        }
    }
    if let Some(fidelity) = recognized_form_items(&root[9], "fidelity") {
        if fidelity.len() >= 3 {
            add_nested_section(&mut nested, &fidelity[0], "baselines", "$[10][1]")?;
            let rungs = cap_recognized_section(
                &fidelity[1],
                "rungs",
                MAX_MACHINE_ASSURANCE_FIDELITY_RUNGS,
                "$[10][2]",
            )?;
            add_nested_section(&mut nested, &fidelity[2], "escalations", "$[10][3]")?;
            if let Some(rows) = rungs {
                for (index, row) in rows.iter().enumerate() {
                    if let Some(args) = recognized_form_items(row, "rung") {
                        if args.len() >= 8 {
                            add_nested_section(
                                &mut nested,
                                &args[6],
                                "falsifiers",
                                &format!("$[10][2][{}][7]", index + 1),
                            )?;
                            add_nested_section(
                                &mut nested,
                                &args[7],
                                "qois",
                                &format!("$[10][2][{}][8]", index + 1),
                            )?;
                        }
                    }
                }
            }
        }
    }
    Ok(())
}

fn cap_recognized_section<'a>(
    node: &'a Node,
    head: &str,
    cap: usize,
    path: &str,
) -> Result<Option<&'a [Node]>, MachineAssuranceCodecError> {
    let Some(items) = recognized_form_items(node, head) else {
        return Ok(None);
    };
    if items.len() > cap {
        return Err(resource_error(
            node,
            path,
            &format!("section {head} count {} exceeds cap {cap}", items.len()),
        ));
    }
    Ok(Some(items))
}

fn add_nested_section<'a>(
    aggregate: &mut usize,
    node: &'a Node,
    head: &str,
    path: &str,
) -> Result<Option<&'a [Node]>, MachineAssuranceCodecError> {
    let Some(items) = recognized_form_items(node, head) else {
        return Ok(None);
    };
    if items.len() > MAX_MACHINE_ASSURANCE_NESTED_REFERENCES {
        return Err(resource_error(
            node,
            path,
            &format!(
                "section {head} count {} exceeds aggregate nested-reference cap {MAX_MACHINE_ASSURANCE_NESTED_REFERENCES}",
                items.len()
            ),
        ));
    }
    *aggregate = aggregate.checked_add(items.len()).ok_or_else(|| {
        resource_error(
            node,
            path,
            "aggregate nested-reference count overflowed usize",
        )
    })?;
    if *aggregate > MAX_MACHINE_ASSURANCE_NESTED_REFERENCES {
        return Err(resource_error(
            node,
            path,
            &format!(
                "aggregate nested-reference count {aggregate} exceeds cap {MAX_MACHINE_ASSURANCE_NESTED_REFERENCES}"
            ),
        ));
    }
    Ok(Some(items))
}

fn exact_form<'a>(
    node: &'a Node,
    expected_head: &str,
    arity: usize,
    path: &str,
) -> Result<&'a [Node], MachineAssuranceCodecError> {
    let (head, args) = raw_form(node, path)?;
    if head != expected_head {
        return Err(unexpected(
            node,
            path,
            &format!("expected form {expected_head}, found {head}"),
            "use the exact positional v1 assurance grammar",
        ));
    }
    if args.len() != arity {
        return Err(unexpected(
            node,
            path,
            &format!(
                "form {expected_head} expects {arity} arguments, found {}",
                args.len()
            ),
            "use the exact positional v1 grammar; omitted values require explicit empty forms or sections",
        ));
    }
    Ok(args)
}

fn raw_form<'a>(
    node: &'a Node,
    path: &str,
) -> Result<(&'a str, &'a [Node]), MachineAssuranceCodecError> {
    let NodeKind::List(items) = &node.kind else {
        return Err(unexpected(
            node,
            path,
            "expected a list form",
            "use one documented assurance form",
        ));
    };
    let Some(head) = items.first() else {
        return Err(unexpected(
            node,
            path,
            "empty lists are not valid assurance forms",
            "start the form with its documented lowercase symbol",
        ));
    };
    let head = symbol(head, &format!("{path}[0]")).map_err(MachineAssuranceCodecError::from)?;
    Ok((head, &items[1..]))
}

fn section_items<'a>(
    node: &'a Node,
    expected_head: &str,
    cap: usize,
    path: &str,
) -> Result<&'a [Node], MachineAssuranceCodecError> {
    let (head, items) = raw_form(node, path)?;
    if head != expected_head {
        return Err(unexpected(
            node,
            path,
            &format!("expected section {expected_head}, found {head}"),
            "use every required section exactly once and in canonical order",
        ));
    }
    if items.len() > cap {
        return Err(resource_error(
            node,
            path,
            &format!(
                "section {expected_head} count {} exceeds cap {cap}",
                items.len()
            ),
        ));
    }
    Ok(items)
}

fn parse_rows<T>(
    rows: &[Node],
    section_node: &Node,
    path: &str,
    parser: fn(&Node, &str) -> Result<T, MachineAssuranceCodecError>,
) -> Result<Vec<T>, MachineAssuranceCodecError> {
    let mut decoded =
        reserved_vec(rows.len(), section_node, path).map_err(MachineAssuranceCodecError::from)?;
    for (index, row) in rows.iter().enumerate() {
        decoded.push(parser(row, &format!("{path}[{}]", index + 1))?);
    }
    Ok(decoded)
}

fn parse_graph_identity(
    node: &Node,
    path: &str,
) -> Result<MachineGraphIdV1, MachineAssuranceCodecError> {
    let args = exact_form(node, "base-graph", 1, path)?;
    let digest = parse_identity_digest(&args[0], &format!("{path}[1]"), "base graph")?;
    MachineGraphIdV1::parse_slice(&digest).ok_or_else(|| {
        reference_error(
            &args[0],
            &format!("{path}[1]"),
            "base graph identity must contain exactly 32 bytes",
        )
    })
}

fn parse_behavior_identity(
    node: &Node,
    path: &str,
) -> Result<MachineBehaviorIdV1, MachineAssuranceCodecError> {
    let args = exact_form(node, "base-behavior", 1, path)?;
    let digest = parse_identity_digest(&args[0], &format!("{path}[1]"), "base behavior")?;
    MachineBehaviorIdV1::parse_slice(&digest).ok_or_else(|| {
        reference_error(
            &args[0],
            &format!("{path}[1]"),
            "base behavior identity must contain exactly 32 bytes",
        )
    })
}

fn parse_identity_digest(
    node: &Node,
    path: &str,
    role: &str,
) -> Result<[u8; 32], MachineAssuranceCodecError> {
    let digest = parse_digest(node, path).map_err(MachineAssuranceCodecError::from)?;
    if digest == [0; 32] {
        return Err(reference_error(
            node,
            path,
            &format!("{role} identity must not be all zero"),
        ));
    }
    Ok(digest)
}

fn parse_u32(node: &Node, path: &str) -> Result<u32, MachineAssuranceCodecError> {
    let value = parse_u64(node, path).map_err(MachineAssuranceCodecError::from)?;
    u32::try_from(value).map_err(|_| {
        number_error(
            node,
            path,
            "unsigned decimal is outside the complete u32 domain",
        )
    })
}

fn parse_evidence_id<T>(
    node: &Node,
    path: &str,
    role: &str,
    constructor: impl FnOnce(String) -> Result<T, fs_evidence::vv::VvErrors>,
) -> Result<T, MachineAssuranceCodecError> {
    let value = string(node, path).map_err(MachineAssuranceCodecError::from)?;
    if value.len() > MAX_VV_ID_BYTES {
        return Err(MachineAssuranceCodecError {
            rule: MachineAssuranceCodecRule::InvalidIdentifier,
            span: node.span,
            path: path.to_string().into_boxed_str(),
            detail: format!("{role} exceeds the {MAX_VV_ID_BYTES}-byte fs-evidence limit")
                .into_boxed_str(),
            hint: "use the bounded canonical fs-evidence identifier grammar"
                .to_string()
                .into_boxed_str(),
        });
    }
    constructor(value.to_string()).map_err(|error| MachineAssuranceCodecError {
        rule: MachineAssuranceCodecRule::InvalidIdentifier,
        span: node.span,
        path: path.to_string().into_boxed_str(),
        detail: format!("invalid {role}: {error}").into_boxed_str(),
        hint: "use the bounded canonical fs-evidence identifier grammar"
            .to_string()
            .into_boxed_str(),
    })
}

fn parse_content_hash(node: &Node, path: &str) -> Result<ContentHash, MachineAssuranceCodecError> {
    let digest = parse_digest(node, path).map_err(MachineAssuranceCodecError::from)?;
    Ok(ContentHash(digest))
}

fn compare_vv_case_bindings(
    declared: &[VvCaseBinding],
    admitted: &[VvCaseBinding],
) -> Result<(), MachineAssuranceAstAdmissionError> {
    let mut declared = declared.to_vec();
    let mut admitted = admitted.to_vec();
    declared.sort_by(compare_vv_binding);
    admitted.sort_by(compare_vv_binding);
    if declared
        .windows(2)
        .any(|pair| pair[0].context == pair[1].context)
    {
        return Err(MachineAssuranceAstAdmissionError::VvCaseBindingMismatch {
            detail: "serialized receipt commitments contain a duplicate Context-of-Use reference; commitments cannot mint evidence authority"
                .to_string()
                .into_boxed_str(),
        });
    }
    if declared != admitted {
        return Err(MachineAssuranceAstAdmissionError::VvCaseBindingMismatch {
            detail: format!(
                "{} serialized receipt commitments do not exactly match {} bindings derived from caller-supplied admitted V&V cases",
                declared.len(),
                admitted.len()
            )
            .into_boxed_str(),
        });
    }
    Ok(())
}

fn compare_vv_binding(left: &VvCaseBinding, right: &VvCaseBinding) -> Ordering {
    left.context
        .cmp(&right.context)
        .then_with(|| left.schema_version.cmp(&right.schema_version))
        .then_with(|| left.ruleset_version.cmp(&right.ruleset_version))
        .then_with(|| left.case_hash.cmp(&right.case_hash))
        .then_with(|| left.receipt_hash.cmp(&right.receipt_hash))
}

fn unexpected(node: &Node, path: &str, detail: &str, hint: &str) -> MachineAssuranceCodecError {
    MachineAssuranceCodecError {
        rule: MachineAssuranceCodecRule::UnexpectedForm,
        span: node.span,
        path: path.to_string().into_boxed_str(),
        detail: detail.to_string().into_boxed_str(),
        hint: hint.to_string().into_boxed_str(),
    }
}

fn unknown_tag(node: &Node, path: &str, family: &str, value: &str) -> MachineAssuranceCodecError {
    MachineAssuranceCodecError {
        rule: MachineAssuranceCodecRule::UnknownTag,
        span: node.span,
        path: path.to_string().into_boxed_str(),
        detail: format!("unknown {family} tag {value:?}").into_boxed_str(),
        hint: format!("use one documented closed {family} tag").into_boxed_str(),
    }
}

fn number_error(node: &Node, path: &str, detail: &str) -> MachineAssuranceCodecError {
    MachineAssuranceCodecError {
        rule: MachineAssuranceCodecRule::InvalidNumber,
        span: node.span,
        path: path.to_string().into_boxed_str(),
        detail: detail.to_string().into_boxed_str(),
        hint: "use the exact canonical unsigned-decimal domain"
            .to_string()
            .into_boxed_str(),
    }
}

fn reference_error(node: &Node, path: &str, detail: &str) -> MachineAssuranceCodecError {
    MachineAssuranceCodecError {
        rule: MachineAssuranceCodecRule::InvalidReference,
        span: node.span,
        path: path.to_string().into_boxed_str(),
        detail: detail.to_string().into_boxed_str(),
        hint: "supply the exact canonical identity or content-bound reference"
            .to_string()
            .into_boxed_str(),
    }
}

fn resource_error(node: &Node, path: &str, detail: &str) -> MachineAssuranceCodecError {
    MachineAssuranceCodecError {
        rule: MachineAssuranceCodecRule::ResourceLimit,
        span: node.span,
        path: path.to_string().into_boxed_str(),
        detail: detail.to_string().into_boxed_str(),
        hint: "split the assurance artifact or reduce the declared collection before admission"
            .to_string()
            .into_boxed_str(),
    }
}

macro_rules! parse_assurance_ref {
    ($node:expr, $role:literal, $path:expr, $ty:ty) => {
        parse_reference($node, $role, $path, |namespace, version, digest| {
            <$ty>::new(namespace, version, digest)
        })
        .map_err(MachineAssuranceCodecError::from)
    };
}

macro_rules! write_ref {
    ($value:expr) => {
        reference_node(
            $value.namespace(),
            $value.schema_version(),
            $value.semantic_digest(),
        )
    };
}

fn parse_artifact_ref(node: &Node, path: &str) -> Result<ArtifactRef, MachineAssuranceCodecError> {
    let args = exact_form(node, "artifact-ref", 3, path)?;
    let kind_path = format!("{path}[1]");
    let kind = match symbol(&args[0], &kind_path).map_err(MachineAssuranceCodecError::from)? {
        "context-of-use" => ArtifactKind::ContextOfUse,
        "validation-plan" => ArtifactKind::ValidationPlan,
        "experiment-artifact" => ArtifactKind::ExperimentArtifact,
        "calibration-split" => ArtifactKind::CalibrationSplit,
        "solution-verification-receipt" => ArtifactKind::SolutionVerificationReceipt,
        "prediction-assessment" => ArtifactKind::PredictionAssessment,
        "assumptions-ledger" => ArtifactKind::AssumptionsLedger,
        value => return Err(unknown_tag(&args[0], &kind_path, "artifact kind", value)),
    };
    Ok(ArtifactRef::new(
        kind,
        parse_evidence_id(
            &args[1],
            &format!("{path}[2]"),
            "artifact-id",
            ArtifactId::try_new,
        )?,
        parse_content_hash(&args[2], &format!("{path}[3]"))?,
    ))
}

fn write_artifact_ref(reference: &ArtifactRef) -> Node {
    form(
        "artifact-ref",
        vec![
            sym(reference.kind().slug()),
            string_node(reference.id().as_str()),
            string_node(&reference.hash().to_hex()),
        ],
    )
}

fn parse_vv_case_binding(
    node: &Node,
    path: &str,
) -> Result<VvCaseBinding, MachineAssuranceCodecError> {
    let args = exact_form(node, "vv-case-receipt", 5, path)?;
    Ok(VvCaseBinding {
        context: parse_artifact_ref(&args[0], &format!("{path}[1]"))?,
        schema_version: parse_u32(&args[1], &format!("{path}[2]"))?,
        ruleset_version: parse_u32(&args[2], &format!("{path}[3]"))?,
        case_hash: parse_content_hash(&args[3], &format!("{path}[4]"))?,
        receipt_hash: parse_content_hash(&args[4], &format!("{path}[5]"))?,
    })
}

fn write_vv_case_binding(binding: &VvCaseBinding) -> Node {
    form(
        "vv-case-receipt",
        vec![
            write_artifact_ref(&binding.context),
            u64_node(u64::from(binding.schema_version)),
            u64_node(u64::from(binding.ruleset_version)),
            string_node(&binding.case_hash.to_hex()),
            string_node(&binding.receipt_hash.to_hex()),
        ],
    )
}

fn parse_sensor(node: &Node, path: &str) -> Result<SensorSpec, MachineAssuranceCodecError> {
    let args = exact_form(node, "sensor", 11, path)?;
    Ok(SensorSpec {
        id: parse_id(&args[0], "sensor-id", &format!("{path}[1]"), SensorId::new)
            .map_err(MachineAssuranceCodecError::from)?,
        owner: parse_id(
            &args[1],
            "subsystem-id",
            &format!("{path}[2]"),
            SubsystemId::new,
        )
        .map_err(MachineAssuranceCodecError::from)?,
        target: parse_observation_target(&args[2], &format!("{path}[3]"))?,
        quantity: parse_quantity(&args[3], &format!("{path}[4]"))
            .map_err(MachineAssuranceCodecError::from)?,
        shape: parse_shape(&args[4], &format!("{path}[5]"))
            .map_err(MachineAssuranceCodecError::from)?,
        clock: parse_id(&args[5], "clock-id", &format!("{path}[6]"), ClockId::new)
            .map_err(MachineAssuranceCodecError::from)?,
        frame: parse_frame(&args[6], &format!("{path}[7]"))
            .map_err(MachineAssuranceCodecError::from)?,
        timing: parse_observation_timing(&args[7], &format!("{path}[8]"))?,
        model: parse_assurance_ref!(
            &args[8],
            "sensor-model-ref",
            &format!("{path}[9]"),
            SensorModelRef
        )?,
        calibration: parse_assurance_ref!(
            &args[9],
            "sensor-calibration-ref",
            &format!("{path}[10]"),
            CalibrationRef
        )?,
        exposure: parse_sensor_exposure(&args[10], &format!("{path}[11]"))?,
    })
}

fn parse_observation_target(
    node: &Node,
    path: &str,
) -> Result<ObservationTarget, MachineAssuranceCodecError> {
    let (head, args) = raw_form(node, path)?;
    if args.len() != 1 {
        return Err(unexpected(
            node,
            path,
            &format!("observation target {head} expects one identifier"),
            "use (terminal \"id\") or (state \"id\")",
        ));
    }
    match head {
        "terminal" => Ok(ObservationTarget::Terminal(
            parse_id(
                &args[0],
                "terminal-id",
                &format!("{path}[1]"),
                TerminalId::new,
            )
            .map_err(MachineAssuranceCodecError::from)?,
        )),
        "state" => Ok(ObservationTarget::State(
            parse_id(
                &args[0],
                "state-slot-id",
                &format!("{path}[1]"),
                StateSlotId::new,
            )
            .map_err(MachineAssuranceCodecError::from)?,
        )),
        _ => Err(unknown_tag(node, path, "observation target", head)),
    }
}

fn parse_observation_timing(
    node: &Node,
    path: &str,
) -> Result<ObservationTiming, MachineAssuranceCodecError> {
    let (head, args) = raw_form(node, path)?;
    match head {
        "direct" if args.is_empty() => Ok(ObservationTiming::Direct),
        "modeled-resampling" if args.len() == 1 => Ok(ObservationTiming::ModeledResampling {
            bridge: parse_assurance_ref!(
                &args[0],
                "sampling-bridge-ref",
                &format!("{path}[1]"),
                SamplingBridgeRef
            )?,
        }),
        "direct" | "modeled-resampling" => Err(unexpected(
            node,
            path,
            &format!("observation timing {head} has invalid arity {}", args.len()),
            "use (direct) or (modeled-resampling (ref ...))",
        )),
        _ => Err(unknown_tag(node, path, "observation timing", head)),
    }
}

fn parse_sensor_exposure(
    node: &Node,
    path: &str,
) -> Result<SensorExposure, MachineAssuranceCodecError> {
    let (head, args) = raw_form(node, path)?;
    match head {
        "plant-signal" if args.len() == 1 => Ok(SensorExposure::PlantSignal {
            output: parse_id(
                &args[0],
                "terminal-id",
                &format!("{path}[1]"),
                TerminalId::new,
            )
            .map_err(MachineAssuranceCodecError::from)?,
        }),
        "experiment-only" if args.is_empty() => Ok(SensorExposure::ExperimentOnly),
        "plant-signal" | "experiment-only" => Err(unexpected(
            node,
            path,
            &format!("sensor exposure {head} has invalid arity {}", args.len()),
            "use (plant-signal \"terminal-id\") or (experiment-only)",
        )),
        _ => Err(unknown_tag(node, path, "sensor exposure", head)),
    }
}

fn write_sensor(sensor: &SensorSpec) -> Node {
    form(
        "sensor",
        vec![
            string_node(sensor.id.canonical_key()),
            string_node(sensor.owner.canonical_key()),
            match &sensor.target {
                ObservationTarget::Terminal(id) => {
                    form("terminal", vec![string_node(id.canonical_key())])
                }
                ObservationTarget::State(id) => {
                    form("state", vec![string_node(id.canonical_key())])
                }
            },
            write_quantity(sensor.quantity),
            write_shape(sensor.shape),
            string_node(sensor.clock.canonical_key()),
            write_frame(&sensor.frame),
            match &sensor.timing {
                ObservationTiming::Direct => form("direct", Vec::new()),
                ObservationTiming::ModeledResampling { bridge } => {
                    form("modeled-resampling", vec![write_ref!(bridge)])
                }
            },
            write_ref!(&sensor.model),
            write_ref!(&sensor.calibration),
            match &sensor.exposure {
                SensorExposure::PlantSignal { output } => {
                    form("plant-signal", vec![string_node(output.canonical_key())])
                }
                SensorExposure::ExperimentOnly => form("experiment-only", Vec::new()),
            },
        ],
    )
}

fn parse_experiment(node: &Node, path: &str) -> Result<ExperimentSpec, MachineAssuranceCodecError> {
    let args = exact_form(node, "experiment", 5, path)?;
    let instruments = section_items(
        &args[3],
        "instruments",
        MAX_MACHINE_ASSURANCE_NESTED_REFERENCES,
        &format!("{path}[4]"),
    )?;
    let qois = section_items(
        &args[4],
        "qois",
        MAX_MACHINE_ASSURANCE_NESTED_REFERENCES,
        &format!("{path}[5]"),
    )?;
    Ok(ExperimentSpec {
        id: parse_id(
            &args[0],
            "experiment-id",
            &format!("{path}[1]"),
            ExperimentId::new,
        )
        .map_err(MachineAssuranceCodecError::from)?,
        artifact: parse_artifact_ref(&args[1], &format!("{path}[2]"))?,
        context: parse_artifact_ref(&args[2], &format!("{path}[3]"))?,
        instruments: parse_rows(
            instruments,
            &args[3],
            &format!("{path}[4]"),
            parse_sensor_instrument,
        )?,
        qois: parse_evidence_ids(
            qois,
            &args[4],
            &format!("{path}[5]"),
            "qoi-id",
            QoiId::try_new,
        )?,
    })
}

fn parse_sensor_instrument(
    node: &Node,
    path: &str,
) -> Result<SensorInstrumentBinding, MachineAssuranceCodecError> {
    let args = exact_form(node, "instrument", 2, path)?;
    Ok(SensorInstrumentBinding {
        sensor: parse_id(&args[0], "sensor-id", &format!("{path}[1]"), SensorId::new)
            .map_err(MachineAssuranceCodecError::from)?,
        instrument: parse_evidence_id(
            &args[1],
            &format!("{path}[2]"),
            "artifact-id",
            ArtifactId::try_new,
        )?,
    })
}

fn parse_evidence_ids<T>(
    rows: &[Node],
    section_node: &Node,
    path: &str,
    role: &str,
    constructor: fn(String) -> Result<T, fs_evidence::vv::VvErrors>,
) -> Result<Vec<T>, MachineAssuranceCodecError> {
    let mut decoded =
        reserved_vec(rows.len(), section_node, path).map_err(MachineAssuranceCodecError::from)?;
    for (index, row) in rows.iter().enumerate() {
        decoded.push(parse_evidence_id(
            row,
            &format!("{path}[{}]", index + 1),
            role,
            constructor,
        )?);
    }
    Ok(decoded)
}

fn write_experiment(experiment: &ExperimentSpec) -> Node {
    form(
        "experiment",
        vec![
            string_node(experiment.id.canonical_key()),
            write_artifact_ref(&experiment.artifact),
            write_artifact_ref(&experiment.context),
            section(
                "instruments",
                experiment
                    .instruments
                    .iter()
                    .map(|binding| {
                        form(
                            "instrument",
                            vec![
                                string_node(binding.sensor.canonical_key()),
                                string_node(binding.instrument.as_str()),
                            ],
                        )
                    })
                    .collect(),
            ),
            section(
                "qois",
                experiment
                    .qois
                    .iter()
                    .map(|qoi| string_node(qoi.as_str()))
                    .collect(),
            ),
        ],
    )
}

fn parse_context(node: &Node, path: &str) -> Result<ContextBinding, MachineAssuranceCodecError> {
    let args = exact_form(node, "context", 4, path)?;
    let qois = section_items(
        &args[2],
        "qois",
        MAX_MACHINE_ASSURANCE_NESTED_REFERENCES,
        &format!("{path}[3]"),
    )?;
    Ok(ContextBinding {
        context: parse_artifact_ref(&args[0], &format!("{path}[1]"))?,
        validation_plan: parse_artifact_ref(&args[1], &format!("{path}[2]"))?,
        qois: parse_rows(qois, &args[2], &format!("{path}[3]"), parse_qoi_binding)?,
        budget: parse_assurance_ref!(
            &args[3],
            "decision-budget-ref",
            &format!("{path}[4]"),
            DecisionBudgetRef
        )?,
    })
}

fn parse_qoi_binding(node: &Node, path: &str) -> Result<QoiBinding, MachineAssuranceCodecError> {
    let args = exact_form(node, "qoi", 5, path)?;
    let inputs = section_items(
        &args[1],
        "inputs",
        MAX_MACHINE_ASSURANCE_NESTED_REFERENCES,
        &format!("{path}[2]"),
    )?;
    Ok(QoiBinding {
        id: parse_evidence_id(&args[0], &format!("{path}[1]"), "qoi-id", QoiId::try_new)?,
        inputs: parse_rows(inputs, &args[1], &format!("{path}[2]"), parse_qoi_input)?,
        unit: parse_evidence_id(&args[2], &format!("{path}[3]"), "unit-id", UnitId::try_new)?,
        definition: parse_assurance_ref!(
            &args[3],
            "qoi-definition-ref",
            &format!("{path}[4]"),
            QoiDefinitionRef
        )?,
        unit_bridge: parse_assurance_ref!(
            &args[4],
            "unit-quantity-bridge-ref",
            &format!("{path}[5]"),
            UnitQuantityBridgeRef
        )?,
    })
}

fn parse_qoi_input(node: &Node, path: &str) -> Result<QoiInput, MachineAssuranceCodecError> {
    let args = exact_form(node, "input", 3, path)?;
    Ok(QoiInput {
        target: parse_qoi_target(&args[0], &format!("{path}[1]"))?,
        quantity: parse_quantity(&args[1], &format!("{path}[2]"))
            .map_err(MachineAssuranceCodecError::from)?,
        shape: parse_shape(&args[2], &format!("{path}[3]"))
            .map_err(MachineAssuranceCodecError::from)?,
    })
}

fn parse_qoi_target(node: &Node, path: &str) -> Result<QoiTarget, MachineAssuranceCodecError> {
    let (head, args) = raw_form(node, path)?;
    if args.len() != 1 {
        return Err(unexpected(
            node,
            path,
            &format!("QoI target {head} expects one identifier"),
            "use (sensor \"id\"), (terminal \"id\"), or (state \"id\")",
        ));
    }
    match head {
        "sensor" => Ok(QoiTarget::Sensor(
            parse_id(&args[0], "sensor-id", &format!("{path}[1]"), SensorId::new)
                .map_err(MachineAssuranceCodecError::from)?,
        )),
        "terminal" => Ok(QoiTarget::Terminal(
            parse_id(
                &args[0],
                "terminal-id",
                &format!("{path}[1]"),
                TerminalId::new,
            )
            .map_err(MachineAssuranceCodecError::from)?,
        )),
        "state" => Ok(QoiTarget::State(
            parse_id(
                &args[0],
                "state-slot-id",
                &format!("{path}[1]"),
                StateSlotId::new,
            )
            .map_err(MachineAssuranceCodecError::from)?,
        )),
        _ => Err(unknown_tag(node, path, "QoI target", head)),
    }
}

fn write_context(context: &ContextBinding) -> Node {
    form(
        "context",
        vec![
            write_artifact_ref(&context.context),
            write_artifact_ref(&context.validation_plan),
            section("qois", context.qois.iter().map(write_qoi_binding).collect()),
            write_ref!(&context.budget),
        ],
    )
}

fn write_qoi_binding(binding: &QoiBinding) -> Node {
    form(
        "qoi",
        vec![
            string_node(binding.id.as_str()),
            section(
                "inputs",
                binding.inputs.iter().map(write_qoi_input).collect(),
            ),
            string_node(binding.unit.as_str()),
            write_ref!(&binding.definition),
            write_ref!(&binding.unit_bridge),
        ],
    )
}

fn write_qoi_input(input: &QoiInput) -> Node {
    form(
        "input",
        vec![
            match &input.target {
                QoiTarget::Sensor(id) => form("sensor", vec![string_node(id.canonical_key())]),
                QoiTarget::Terminal(id) => form("terminal", vec![string_node(id.canonical_key())]),
                QoiTarget::State(id) => form("state", vec![string_node(id.canonical_key())]),
            },
            write_quantity(input.quantity),
            write_shape(input.shape),
        ],
    )
}

fn parse_hazard(node: &Node, path: &str) -> Result<HazardSpec, MachineAssuranceCodecError> {
    let args = exact_form(node, "hazard", 8, path)?;
    let scopes = section_items(
        &args[2],
        "scope",
        MAX_MACHINE_ASSURANCE_NESTED_REFERENCES,
        &format!("{path}[3]"),
    )?;
    let assumptions = section_items(
        &args[6],
        "assumptions",
        MAX_MACHINE_ASSURANCE_NESTED_REFERENCES,
        &format!("{path}[7]"),
    )?;
    Ok(HazardSpec {
        id: parse_id(&args[0], "hazard-id", &format!("{path}[1]"), HazardId::new)
            .map_err(MachineAssuranceCodecError::from)?,
        context: parse_artifact_ref(&args[1], &format!("{path}[2]"))?,
        scope: parse_rows(scopes, &args[2], &format!("{path}[3]"), parse_machine_scope)?,
        requirement: parse_assurance_ref!(
            &args[3],
            "safety-requirement-ref",
            &format!("{path}[4]"),
            SafetyRequirementRef
        )?,
        operating_envelope: parse_assurance_ref!(
            &args[4],
            "operating-envelope-ref",
            &format!("{path}[5]"),
            OperatingEnvelopeRef
        )?,
        safety_case: parse_assurance_ref!(
            &args[5],
            "safety-case-ref",
            &format!("{path}[6]"),
            SafetyCaseRef
        )?,
        assumptions: parse_evidence_ids(
            assumptions,
            &args[6],
            &format!("{path}[7]"),
            "assumption-id",
            AssumptionId::try_new,
        )?,
        fault_coverage: parse_fault_coverage(&args[7], &format!("{path}[8]"))?,
    })
}

fn parse_fault_coverage(
    node: &Node,
    path: &str,
) -> Result<FaultCoverage, MachineAssuranceCodecError> {
    let (head, args) = raw_form(node, path)?;
    match head {
        "modeled" if args.is_empty() => Ok(FaultCoverage::Modeled),
        "unmodeled" if args.len() == 1 => Ok(FaultCoverage::Unmodeled(parse_assurance_ref!(
            &args[0],
            "assurance-no-claim-ref",
            &format!("{path}[1]"),
            NoClaimRef
        )?)),
        "modeled" | "unmodeled" => Err(unexpected(
            node,
            path,
            &format!("fault coverage {head} has invalid arity {}", args.len()),
            "use (modeled) or (unmodeled (ref ...))",
        )),
        _ => Err(unknown_tag(node, path, "fault coverage", head)),
    }
}

fn parse_machine_scope(
    node: &Node,
    path: &str,
) -> Result<MachineScope, MachineAssuranceCodecError> {
    let (head, args) = raw_form(node, path)?;
    match head {
        "whole-machine" if args.is_empty() => Ok(MachineScope::WholeMachine),
        "subsystem" if args.len() == 1 => Ok(MachineScope::Subsystem(
            parse_id(
                &args[0],
                "subsystem-id",
                &format!("{path}[1]"),
                SubsystemId::new,
            )
            .map_err(MachineAssuranceCodecError::from)?,
        )),
        "element" if args.len() == 1 => Ok(MachineScope::Element(parse_machine_element_local(
            &args[0],
            &format!("{path}[1]"),
        )?)),
        "relation" if args.len() == 1 => Ok(MachineScope::Relation(
            parse_id(
                &args[0],
                "relation-id",
                &format!("{path}[1]"),
                RelationId::new,
            )
            .map_err(MachineAssuranceCodecError::from)?,
        )),
        "interface" if args.len() == 1 => Ok(MachineScope::Interface(
            parse_id(
                &args[0],
                "interface-id",
                &format!("{path}[1]"),
                InterfaceId::new,
            )
            .map_err(MachineAssuranceCodecError::from)?,
        )),
        "whole-machine" | "subsystem" | "element" | "relation" | "interface" => Err(unexpected(
            node,
            path,
            &format!("machine scope {head} has invalid arity {}", args.len()),
            "use one exact documented Machine scope form",
        )),
        _ => Err(unknown_tag(node, path, "machine scope", head)),
    }
}

fn parse_machine_element_local(
    node: &Node,
    path: &str,
) -> Result<MachineElementId, MachineAssuranceCodecError> {
    let (head, args) = raw_form(node, path)?;
    if args.len() != 1 {
        return Err(unexpected(
            node,
            path,
            &format!("machine element {head} expects one identifier"),
            "use one exact documented Machine element form",
        ));
    }
    let id_path = format!("{path}[1]");
    match head {
        "body" => Ok(MachineElementId::Body(
            parse_id(&args[0], "body-id", &id_path, BodyId::new)
                .map_err(MachineAssuranceCodecError::from)?,
        )),
        "surface-patch" => Ok(MachineElementId::SurfacePatch(
            parse_id(&args[0], "surface-patch-id", &id_path, SurfacePatchId::new)
                .map_err(MachineAssuranceCodecError::from)?,
        )),
        "contact-feature" => Ok(MachineElementId::ContactFeature(
            parse_id(
                &args[0],
                "contact-feature-id",
                &id_path,
                ContactFeatureId::new,
            )
            .map_err(MachineAssuranceCodecError::from)?,
        )),
        "terminal" => Ok(MachineElementId::Terminal(
            parse_id(&args[0], "terminal-id", &id_path, TerminalId::new)
                .map_err(MachineAssuranceCodecError::from)?,
        )),
        "port" => Ok(MachineElementId::Port(
            parse_id(&args[0], "port-id", &id_path, PortId::new)
                .map_err(MachineAssuranceCodecError::from)?,
        )),
        "state-slot" => Ok(MachineElementId::StateSlot(
            parse_id(&args[0], "state-slot-id", &id_path, StateSlotId::new)
                .map_err(MachineAssuranceCodecError::from)?,
        )),
        _ => Err(unknown_tag(node, path, "machine element", head)),
    }
}

fn write_hazard(hazard: &HazardSpec) -> Node {
    form(
        "hazard",
        vec![
            string_node(hazard.id.canonical_key()),
            write_artifact_ref(&hazard.context),
            section(
                "scope",
                hazard.scope.iter().map(write_machine_scope).collect(),
            ),
            write_ref!(&hazard.requirement),
            write_ref!(&hazard.operating_envelope),
            write_ref!(&hazard.safety_case),
            section(
                "assumptions",
                hazard
                    .assumptions
                    .iter()
                    .map(|id| string_node(id.as_str()))
                    .collect(),
            ),
            match &hazard.fault_coverage {
                FaultCoverage::Modeled => form("modeled", Vec::new()),
                FaultCoverage::Unmodeled(no_claim) => form("unmodeled", vec![write_ref!(no_claim)]),
            },
        ],
    )
}

fn write_machine_scope(scope: &MachineScope) -> Node {
    match scope {
        MachineScope::WholeMachine => form("whole-machine", Vec::new()),
        MachineScope::Subsystem(id) => form("subsystem", vec![string_node(id.canonical_key())]),
        MachineScope::Element(id) => form("element", vec![write_machine_element(id)]),
        MachineScope::Relation(id) => form("relation", vec![string_node(id.canonical_key())]),
        MachineScope::Interface(id) => form("interface", vec![string_node(id.canonical_key())]),
    }
}

fn parse_fault(node: &Node, path: &str) -> Result<FaultSpec, MachineAssuranceCodecError> {
    let args = exact_form(node, "fault", 6, path)?;
    let affected = section_items(
        &args[1],
        "affected",
        MAX_MACHINE_ASSURANCE_NESTED_REFERENCES,
        &format!("{path}[2]"),
    )?;
    let hazards = section_items(
        &args[2],
        "hazards",
        MAX_MACHINE_ASSURANCE_NESTED_REFERENCES,
        &format!("{path}[3]"),
    )?;
    let mut parsed_hazards = reserved_vec(hazards.len(), &args[2], &format!("{path}[3]"))
        .map_err(MachineAssuranceCodecError::from)?;
    for (index, hazard) in hazards.iter().enumerate() {
        parsed_hazards.push(
            parse_id(
                hazard,
                "hazard-id",
                &format!("{path}[3][{}]", index + 1),
                HazardId::new,
            )
            .map_err(MachineAssuranceCodecError::from)?,
        );
    }
    Ok(FaultSpec {
        id: parse_id(&args[0], "fault-id", &format!("{path}[1]"), FaultId::new)
            .map_err(MachineAssuranceCodecError::from)?,
        affected: parse_rows(
            affected,
            &args[1],
            &format!("{path}[2]"),
            parse_machine_scope,
        )?,
        hazards: parsed_hazards,
        model: parse_assurance_ref!(
            &args[3],
            "fault-model-ref",
            &format!("{path}[4]"),
            FaultModelRef
        )?,
        containment: parse_assurance_ref!(
            &args[4],
            "fault-containment-ref",
            &format!("{path}[5]"),
            FaultContainmentRef
        )?,
        injection: parse_assurance_ref!(
            &args[5],
            "fault-injection-ref",
            &format!("{path}[6]"),
            FaultInjectionRef
        )?,
    })
}

fn write_fault(fault: &FaultSpec) -> Node {
    form(
        "fault",
        vec![
            string_node(fault.id.canonical_key()),
            section(
                "affected",
                fault.affected.iter().map(write_machine_scope).collect(),
            ),
            section(
                "hazards",
                fault
                    .hazards
                    .iter()
                    .map(|id| string_node(id.canonical_key()))
                    .collect(),
            ),
            write_ref!(&fault.model),
            write_ref!(&fault.containment),
            write_ref!(&fault.injection),
        ],
    )
}

fn parse_accounting_window(
    node: &Node,
    path: &str,
) -> Result<AccountingWindow, MachineAssuranceCodecError> {
    let args = exact_form(node, "accounting-window", 9, path)?;
    let entries = section_items(
        &args[7],
        "entries",
        MAX_MACHINE_ASSURANCE_NESTED_REFERENCES,
        &format!("{path}[8]"),
    )?;
    Ok(AccountingWindow {
        id: parse_id(
            &args[0],
            "accounting-window-id",
            &format!("{path}[1]"),
            AccountingWindowId::new,
        )
        .map_err(MachineAssuranceCodecError::from)?,
        context: parse_artifact_ref(&args[1], &format!("{path}[2]"))?,
        clock: parse_id(&args[2], "clock-id", &format!("{path}[3]"), ClockId::new)
            .map_err(MachineAssuranceCodecError::from)?,
        balance: parse_balance_kind(&args[3], &format!("{path}[4]"))?,
        quantity: parse_quantity(&args[4], &format!("{path}[5]"))
            .map_err(MachineAssuranceCodecError::from)?,
        boundary: parse_assurance_ref!(
            &args[5],
            "accounting-boundary-ref",
            &format!("{path}[6]"),
            AccountingBoundaryRef
        )?,
        interval: parse_assurance_ref!(
            &args[6],
            "accounting-interval-ref",
            &format!("{path}[7]"),
            AccountingIntervalRef
        )?,
        entries: parse_rows(
            entries,
            &args[7],
            &format!("{path}[8]"),
            parse_accounting_entry,
        )?,
        audit_policy: parse_assurance_ref!(
            &args[8],
            "accounting-policy-ref",
            &format!("{path}[9]"),
            AccountingPolicyRef
        )?,
    })
}

fn parse_balance_kind(node: &Node, path: &str) -> Result<BalanceKind, MachineAssuranceCodecError> {
    let (head, args) = raw_form(node, path)?;
    let no_args = |value| {
        if args.is_empty() {
            Ok(value)
        } else {
            Err(unexpected(
                node,
                path,
                &format!("balance kind {head} expects no arguments"),
                "use the exact documented balance-kind form",
            ))
        }
    };
    match head {
        "energy" => no_args(BalanceKind::Energy),
        "enthalpy" => no_args(BalanceKind::Enthalpy),
        "linear-momentum" => no_args(BalanceKind::LinearMomentum),
        "angular-momentum" => no_args(BalanceKind::AngularMomentum),
        "mass" => no_args(BalanceKind::Mass),
        "electric-charge" => no_args(BalanceKind::ElectricCharge),
        "amount-of-substance" => no_args(BalanceKind::AmountOfSubstance),
        "entropy" => no_args(BalanceKind::Entropy),
        "exergy" => no_args(BalanceKind::Exergy),
        "species" if args.len() == 1 => Ok(BalanceKind::Species(parse_assurance_ref!(
            &args[0],
            "balance-law-ref",
            &format!("{path}[1]"),
            BalanceLawRef
        )?)),
        "elements" if args.len() == 1 => Ok(BalanceKind::Elements(parse_assurance_ref!(
            &args[0],
            "balance-law-ref",
            &format!("{path}[1]"),
            BalanceLawRef
        )?)),
        "custom" if args.len() == 1 => Ok(BalanceKind::Custom(parse_assurance_ref!(
            &args[0],
            "balance-law-ref",
            &format!("{path}[1]"),
            BalanceLawRef
        )?)),
        "species" | "elements" | "custom" => Err(unexpected(
            node,
            path,
            &format!("balance kind {head} expects one balance-law reference"),
            "use the exact documented balance-kind form",
        )),
        _ => Err(unknown_tag(node, path, "balance kind", head)),
    }
}

fn parse_accounting_entry(
    node: &Node,
    path: &str,
) -> Result<AccountingEntry, MachineAssuranceCodecError> {
    let args = exact_form(node, "entry", 5, path)?;
    Ok(AccountingEntry {
        target: parse_accounting_target(&args[0], &format!("{path}[1]"))?,
        role: parse_accounting_role(&args[1], &format!("{path}[2]"))?,
        orientation: parse_accounting_orientation(&args[2], &format!("{path}[3]"))?,
        policy: parse_assurance_ref!(
            &args[3],
            "accounting-policy-ref",
            &format!("{path}[4]"),
            AccountingPolicyRef
        )?,
        loss_ownership: parse_loss_ownership(&args[4], &format!("{path}[5]"))?,
    })
}

fn parse_accounting_target(
    node: &Node,
    path: &str,
) -> Result<AccountingTarget, MachineAssuranceCodecError> {
    let (head, args) = raw_form(node, path)?;
    if args.len() != 1 {
        return Err(unexpected(
            node,
            path,
            &format!("accounting target {head} expects one identifier"),
            "use one exact documented accounting-target form",
        ));
    }
    match head {
        "relation" => Ok(AccountingTarget::Relation(
            parse_id(
                &args[0],
                "relation-id",
                &format!("{path}[1]"),
                RelationId::new,
            )
            .map_err(MachineAssuranceCodecError::from)?,
        )),
        "port" => Ok(AccountingTarget::Port(
            parse_id(&args[0], "port-id", &format!("{path}[1]"), PortId::new)
                .map_err(MachineAssuranceCodecError::from)?,
        )),
        "interface" => Ok(AccountingTarget::Interface(
            parse_id(
                &args[0],
                "interface-id",
                &format!("{path}[1]"),
                InterfaceId::new,
            )
            .map_err(MachineAssuranceCodecError::from)?,
        )),
        "state" => Ok(AccountingTarget::State(
            parse_id(
                &args[0],
                "state-slot-id",
                &format!("{path}[1]"),
                StateSlotId::new,
            )
            .map_err(MachineAssuranceCodecError::from)?,
        )),
        _ => Err(unknown_tag(node, path, "accounting target", head)),
    }
}

fn parse_accounting_role(
    node: &Node,
    path: &str,
) -> Result<AccountingRole, MachineAssuranceCodecError> {
    match symbol(node, path).map_err(MachineAssuranceCodecError::from)? {
        "storage" => Ok(AccountingRole::Storage),
        "dissipation" => Ok(AccountingRole::Dissipation),
        "included-source" => Ok(AccountingRole::IncludedSource),
        "external-exchange" => Ok(AccountingRole::ExternalExchange),
        "stream" => Ok(AccountingRole::Stream),
        value => Err(unknown_tag(node, path, "accounting role", value)),
    }
}

fn parse_accounting_orientation(
    node: &Node,
    path: &str,
) -> Result<AccountingOrientation, MachineAssuranceCodecError> {
    match symbol(node, path).map_err(MachineAssuranceCodecError::from)? {
        "stored-increase-positive" => Ok(AccountingOrientation::StoredIncreasePositive),
        "nonnegative-loss" => Ok(AccountingOrientation::NonnegativeLoss),
        "into-boundary-positive" => Ok(AccountingOrientation::IntoBoundaryPositive),
        "out-of-boundary-positive" => Ok(AccountingOrientation::OutOfBoundaryPositive),
        value => Err(unknown_tag(node, path, "accounting orientation", value)),
    }
}

fn parse_loss_ownership(
    node: &Node,
    path: &str,
) -> Result<Option<LossOwnershipRef>, MachineAssuranceCodecError> {
    let (head, args) = raw_form(node, path)?;
    match head {
        "no-loss-owner" if args.is_empty() => Ok(None),
        "loss-owner" if args.len() == 1 => Ok(Some(parse_assurance_ref!(
            &args[0],
            "loss-ownership-ref",
            &format!("{path}[1]"),
            LossOwnershipRef
        )?)),
        "no-loss-owner" | "loss-owner" => Err(unexpected(
            node,
            path,
            &format!("loss ownership {head} has invalid arity {}", args.len()),
            "use (no-loss-owner) or (loss-owner (ref ...))",
        )),
        _ => Err(unknown_tag(node, path, "loss ownership", head)),
    }
}

fn write_accounting_window(window: &AccountingWindow) -> Node {
    form(
        "accounting-window",
        vec![
            string_node(window.id.canonical_key()),
            write_artifact_ref(&window.context),
            string_node(window.clock.canonical_key()),
            write_balance_kind(&window.balance),
            write_quantity(window.quantity),
            write_ref!(&window.boundary),
            write_ref!(&window.interval),
            section(
                "entries",
                window.entries.iter().map(write_accounting_entry).collect(),
            ),
            write_ref!(&window.audit_policy),
        ],
    )
}

fn write_balance_kind(balance: &BalanceKind) -> Node {
    match balance {
        BalanceKind::Energy => form("energy", Vec::new()),
        BalanceKind::Enthalpy => form("enthalpy", Vec::new()),
        BalanceKind::LinearMomentum => form("linear-momentum", Vec::new()),
        BalanceKind::AngularMomentum => form("angular-momentum", Vec::new()),
        BalanceKind::Mass => form("mass", Vec::new()),
        BalanceKind::ElectricCharge => form("electric-charge", Vec::new()),
        BalanceKind::AmountOfSubstance => form("amount-of-substance", Vec::new()),
        BalanceKind::Species(law) => form("species", vec![write_ref!(law)]),
        BalanceKind::Elements(law) => form("elements", vec![write_ref!(law)]),
        BalanceKind::Entropy => form("entropy", Vec::new()),
        BalanceKind::Exergy => form("exergy", Vec::new()),
        BalanceKind::Custom(law) => form("custom", vec![write_ref!(law)]),
    }
}

fn write_accounting_entry(entry: &AccountingEntry) -> Node {
    form(
        "entry",
        vec![
            match &entry.target {
                AccountingTarget::Relation(id) => {
                    form("relation", vec![string_node(id.canonical_key())])
                }
                AccountingTarget::Port(id) => form("port", vec![string_node(id.canonical_key())]),
                AccountingTarget::Interface(id) => {
                    form("interface", vec![string_node(id.canonical_key())])
                }
                AccountingTarget::State(id) => form("state", vec![string_node(id.canonical_key())]),
            },
            sym(match entry.role {
                AccountingRole::Storage => "storage",
                AccountingRole::Dissipation => "dissipation",
                AccountingRole::IncludedSource => "included-source",
                AccountingRole::ExternalExchange => "external-exchange",
                AccountingRole::Stream => "stream",
            }),
            sym(match entry.orientation {
                AccountingOrientation::StoredIncreasePositive => "stored-increase-positive",
                AccountingOrientation::NonnegativeLoss => "nonnegative-loss",
                AccountingOrientation::IntoBoundaryPositive => "into-boundary-positive",
                AccountingOrientation::OutOfBoundaryPositive => "out-of-boundary-positive",
            }),
            write_ref!(&entry.policy),
            match &entry.loss_ownership {
                Some(owner) => form("loss-owner", vec![write_ref!(owner)]),
                None => form("no-loss-owner", Vec::new()),
            },
        ],
    )
}

fn parse_fidelity(node: &Node, path: &str) -> Result<FidelityPolicy, MachineAssuranceCodecError> {
    let args = exact_form(node, "fidelity", 4, path)?;
    let baselines = section_items(
        &args[0],
        "baselines",
        MAX_MACHINE_ASSURANCE_NESTED_REFERENCES,
        &format!("{path}[1]"),
    )?;
    let rungs = section_items(
        &args[1],
        "rungs",
        MAX_MACHINE_ASSURANCE_FIDELITY_RUNGS,
        &format!("{path}[2]"),
    )?;
    let escalations = section_items(
        &args[2],
        "escalations",
        MAX_MACHINE_ASSURANCE_NESTED_REFERENCES,
        &format!("{path}[3]"),
    )?;
    let fixed = exact_form(&args[3], "fixed-replay", 1, &format!("{path}[4]"))?;

    let mut parsed_baselines = reserved_vec(baselines.len(), &args[0], &format!("{path}[1]"))
        .map_err(MachineAssuranceCodecError::from)?;
    for (index, baseline) in baselines.iter().enumerate() {
        parsed_baselines.push(
            parse_id(
                baseline,
                "fidelity-rung-id",
                &format!("{path}[1][{}]", index + 1),
                FidelityRungId::new,
            )
            .map_err(MachineAssuranceCodecError::from)?,
        );
    }
    Ok(FidelityPolicy {
        baselines: parsed_baselines,
        rungs: parse_rows(rungs, &args[1], &format!("{path}[2]"), parse_fidelity_rung)?,
        escalations: parse_rows(
            escalations,
            &args[2],
            &format!("{path}[3]"),
            parse_escalation,
        )?,
        fixed_replay: parse_assurance_ref!(
            &fixed[0],
            "fixed-replay-ref",
            &format!("{path}[4][1]"),
            FixedReplayRef
        )?,
    })
}

fn parse_fidelity_rung(
    node: &Node,
    path: &str,
) -> Result<FidelityRung, MachineAssuranceCodecError> {
    let args = exact_form(node, "rung", 8, path)?;
    let falsifiers = section_items(
        &args[6],
        "falsifiers",
        MAX_MACHINE_ASSURANCE_NESTED_REFERENCES,
        &format!("{path}[7]"),
    )?;
    let qois = section_items(
        &args[7],
        "qois",
        MAX_MACHINE_ASSURANCE_NESTED_REFERENCES,
        &format!("{path}[8]"),
    )?;
    let mut parsed_falsifiers = reserved_vec(falsifiers.len(), &args[6], &format!("{path}[7]"))
        .map_err(MachineAssuranceCodecError::from)?;
    for (index, falsifier) in falsifiers.iter().enumerate() {
        parsed_falsifiers.push(parse_assurance_ref!(
            falsifier,
            "falsifier-ref",
            &format!("{path}[7][{}]", index + 1),
            FalsifierRef
        )?);
    }
    Ok(FidelityRung {
        id: parse_id(
            &args[0],
            "fidelity-rung-id",
            &format!("{path}[1]"),
            FidelityRungId::new,
        )
        .map_err(MachineAssuranceCodecError::from)?,
        subsystem: parse_id(
            &args[1],
            "subsystem-id",
            &format!("{path}[2]"),
            SubsystemId::new,
        )
        .map_err(MachineAssuranceCodecError::from)?,
        model: parse_assurance_ref!(&args[2], "model-ref", &format!("{path}[3]"), ModelRef)?,
        model_crosswalk: parse_assurance_ref!(
            &args[3],
            "model-crosswalk-ref",
            &format!("{path}[4]"),
            ModelCrosswalkRef
        )?,
        validity_domain: parse_assurance_ref!(
            &args[4],
            "validity-domain-ref",
            &format!("{path}[5]"),
            ValidityDomainRef
        )?,
        cost_error_model: parse_assurance_ref!(
            &args[5],
            "cost-error-model-ref",
            &format!("{path}[6]"),
            CostErrorModelRef
        )?,
        falsifiers: parsed_falsifiers,
        qois: parse_rows(qois, &args[7], &format!("{path}[8]"), parse_context_qoi)?,
    })
}

fn parse_context_qoi(node: &Node, path: &str) -> Result<ContextQoiKey, MachineAssuranceCodecError> {
    let args = exact_form(node, "context-qoi", 2, path)?;
    Ok(ContextQoiKey {
        context: parse_evidence_id(
            &args[0],
            &format!("{path}[1]"),
            "artifact-id",
            ArtifactId::try_new,
        )?,
        qoi: parse_evidence_id(&args[1], &format!("{path}[2]"), "qoi-id", QoiId::try_new)?,
    })
}

fn parse_escalation(node: &Node, path: &str) -> Result<EscalationSpec, MachineAssuranceCodecError> {
    let args = exact_form(node, "escalation", 3, path)?;
    Ok(EscalationSpec {
        from: parse_id(
            &args[0],
            "fidelity-rung-id",
            &format!("{path}[1]"),
            FidelityRungId::new,
        )
        .map_err(MachineAssuranceCodecError::from)?,
        trigger: parse_assurance_ref!(
            &args[1],
            "escalation-trigger-ref",
            &format!("{path}[2]"),
            EscalationTriggerRef
        )?,
        action: parse_escalation_action(&args[2], &format!("{path}[3]"))?,
    })
}

fn parse_escalation_action(
    node: &Node,
    path: &str,
) -> Result<EscalationAction, MachineAssuranceCodecError> {
    let (head, args) = raw_form(node, path)?;
    match head {
        "escalate" if args.len() == 3 => Ok(EscalationAction::Escalate {
            target: parse_id(
                &args[0],
                "fidelity-rung-id",
                &format!("{path}[1]"),
                FidelityRungId::new,
            )
            .map_err(MachineAssuranceCodecError::from)?,
            transfer: parse_assurance_ref!(
                &args[1],
                "state-transfer-ref",
                &format!("{path}[2]"),
                StateTransferRef
            )?,
            crosswalk: parse_assurance_ref!(
                &args[2],
                "model-crosswalk-ref",
                &format!("{path}[3]"),
                ModelCrosswalkRef
            )?,
        }),
        "refuse" if args.len() == 1 => Ok(EscalationAction::Refuse(parse_assurance_ref!(
            &args[0],
            "assurance-no-claim-ref",
            &format!("{path}[1]"),
            NoClaimRef
        )?)),
        "escalate" | "refuse" => Err(unexpected(
            node,
            path,
            &format!("escalation action {head} has invalid arity {}", args.len()),
            "use (escalate \"target\" (ref ...) (ref ...)) or (refuse (ref ...))",
        )),
        _ => Err(unknown_tag(node, path, "escalation action", head)),
    }
}

fn write_fidelity(fidelity: &FidelityPolicy) -> Node {
    form(
        "fidelity",
        vec![
            section(
                "baselines",
                fidelity
                    .baselines
                    .iter()
                    .map(|id| string_node(id.canonical_key()))
                    .collect(),
            ),
            section(
                "rungs",
                fidelity.rungs.iter().map(write_fidelity_rung).collect(),
            ),
            section(
                "escalations",
                fidelity.escalations.iter().map(write_escalation).collect(),
            ),
            form("fixed-replay", vec![write_ref!(&fidelity.fixed_replay)]),
        ],
    )
}

fn write_fidelity_rung(rung: &FidelityRung) -> Node {
    form(
        "rung",
        vec![
            string_node(rung.id.canonical_key()),
            string_node(rung.subsystem.canonical_key()),
            write_ref!(&rung.model),
            write_ref!(&rung.model_crosswalk),
            write_ref!(&rung.validity_domain),
            write_ref!(&rung.cost_error_model),
            section(
                "falsifiers",
                rung.falsifiers
                    .iter()
                    .map(|value| write_ref!(value))
                    .collect(),
            ),
            section("qois", rung.qois.iter().map(write_context_qoi).collect()),
        ],
    )
}

fn write_context_qoi(key: &ContextQoiKey) -> Node {
    form(
        "context-qoi",
        vec![
            string_node(key.context.as_str()),
            string_node(key.qoi.as_str()),
        ],
    )
}

fn write_escalation(escalation: &EscalationSpec) -> Node {
    form(
        "escalation",
        vec![
            string_node(escalation.from.canonical_key()),
            write_ref!(&escalation.trigger),
            match &escalation.action {
                EscalationAction::Escalate {
                    target,
                    transfer,
                    crosswalk,
                } => form(
                    "escalate",
                    vec![
                        string_node(target.canonical_key()),
                        write_ref!(transfer),
                        write_ref!(crosswalk),
                    ],
                ),
                EscalationAction::Refuse(no_claim) => form("refuse", vec![write_ref!(no_claim)]),
            },
        ],
    )
}
