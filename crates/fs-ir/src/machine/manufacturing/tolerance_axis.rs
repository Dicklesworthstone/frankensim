//! Machine-IR binding for correlated manufacturing-tolerance axes.
//!
//! `fs-toleralloc` deliberately retains positional axes without claiming that
//! their labels match a Machine-IR behavior declaration. This L6 adapter closes
//! that gap: it binds the canonical correlated dependence members of one
//! admitted [`AdmittedMachineBehavior`] to the exact positional terms and
//! factor carried by one non-forgeable [`CorrelatedStackReceipt`]. It also
//! closes the behavior and manufacturing receipts over the same Machine graph.
//! Because the manufacturing state uses a content-addressed coordinate while
//! `fs-toleralloc` retains a caller-supplied semantic digest, the adapter never
//! treats those byte domains as interchangeable: an additional versioned
//! artifact explicitly names the caller's cross-domain coordinate-link policy.
//!
//! The crosswalk is structural. It does not establish that a tolerance's
//! quantity scale is the stack term's standard deviation, prove unit closure,
//! authenticate the external model digest, validate a manufacturing
//! population, or promote first-order moments into nonlinear, tail, quantile,
//! reliability, or gear-backlash authority.

use core::fmt;

use fs_blake3::identity::{
    CanonicalEncoder, CanonicalError, CanonicalLimits, CanonicalSchema, Field, FieldSpec,
    IdentityReceipt, NeverCancel, ProblemSemanticId, StrongIdentity, WireType,
};
use fs_evidence::ColorRank;
use fs_toleralloc::{AdmittedCorrelationModel, CorrelatedStackReceipt, CorrelatedStackTerm};

use crate::IR_VERSION;

use super::{
    AdmittedMachineManufacturingStateV1, MachineManufacturingStateIdV1, ManufacturingArtifactRefV1,
    append_bytes,
};
use crate::machine::semantics::{
    AdmittedMachineBehavior, CorrelationModelRef, DependenceMember, DependenceModel,
    MachineBehaviorIdV1, ToleranceId, ToleranceSemantics, ToleranceTarget,
};
use crate::machine::{
    BodyId, DependentBinding, DependentKind, MachineElementId, MachineGraphIdV1, TerminalShape,
};

/// Identity-schema version for the Machine-IR correlated-tolerance crosswalk.
pub const MACHINE_TOLERANCE_AXIS_CROSSWALK_SCHEMA_VERSION_V1: u32 = 1;

const TOLERANCE_AXIS_CROSSWALK_IDENTITY_LIMITS: CanonicalLimits =
    CanonicalLimits::new(2 * 1_024 * 1_024, 1 * 1_024 * 1_024, 10, 512, 256);

/// Closed refusal vocabulary for correlated-tolerance axis binding.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MachineToleranceAxisCrosswalkErrorV1 {
    /// The behavior and manufacturing state extend different Machine graphs.
    GraphMismatch {
        /// Graph named by the admitted behavior.
        behavior: MachineGraphIdV1,
        /// Graph named by the admitted manufacturing state.
        manufacturing: MachineGraphIdV1,
    },
    /// The behavior has no random-source dependence declaration to bind.
    MissingDependence,
    /// An independent declaration cannot be reinterpreted as correlated.
    IndependentDependence,
    /// Machine-IR and `fs-toleralloc` expose different positional dimensions.
    AxisDimensionMismatch {
        /// Canonical Machine-IR dependence-member count.
        behavior: usize,
        /// Retained correlated-stack term count.
        stack: usize,
    },
    /// A random condition occupies a model position; v1 has no admitted
    /// projection from mixed uncertainty into a manufacturing-only stack.
    ConditionAxisUnsupported {
        /// Canonical zero-based model position.
        index: usize,
    },
    /// The correlated group contains no Machine-IR tolerance to crosswalk.
    NoToleranceAxes,
    /// The behavior correlation reference differs from the stack model.
    BehaviorCorrelationModelMismatch,
    /// A dependence tolerance did not resolve to its admitted specification.
    MissingToleranceSpec {
        /// Canonical zero-based model position.
        index: usize,
        /// Missing durable tolerance.
        tolerance: ToleranceId,
    },
    /// A future or malformed behavior exposed a non-random dependence member.
    NonRandomTolerance {
        /// Canonical zero-based model position.
        index: usize,
        /// Offending durable tolerance.
        tolerance: ToleranceId,
    },
    /// V1 cannot bind vector/tensor components without basis/order semantics.
    NonScalarTolerance {
        /// Canonical zero-based model position.
        index: usize,
        /// Offending durable tolerance.
        tolerance: ToleranceId,
    },
    /// V1 manufacturing history is body-level, so other target kinds refuse.
    UnsupportedToleranceTarget {
        /// Canonical zero-based model position.
        index: usize,
        /// Offending durable tolerance.
        tolerance: ToleranceId,
    },
    /// The body-bound tolerance has no retained manufacturing process history.
    MissingManufacturingHistory {
        /// Canonical zero-based model position.
        index: usize,
        /// Offending durable tolerance.
        tolerance: ToleranceId,
        /// Body absent from the manufacturing-state process set.
        body: BodyId,
    },
    /// A positional stack label does not exactly name its bound tolerance ID.
    AxisNameMismatch {
        /// Canonical zero-based model position.
        index: usize,
        /// Machine-IR tolerance assigned to that position.
        tolerance: ToleranceId,
        /// Exact bounded label retained by `fs-toleralloc`.
        stack_name: String,
    },
    /// Canonical identity publication failed.
    Identity(CanonicalError),
}

impl MachineToleranceAxisCrosswalkErrorV1 {
    /// Stable diagnostic code.
    #[must_use]
    pub const fn code(&self) -> &'static str {
        match self {
            Self::GraphMismatch { .. } => "MachineToleranceAxisGraphMismatch",
            Self::MissingDependence => "MachineToleranceAxisMissingDependence",
            Self::IndependentDependence => "MachineToleranceAxisIndependentDependence",
            Self::AxisDimensionMismatch { .. } => "MachineToleranceAxisDimensionMismatch",
            Self::ConditionAxisUnsupported { .. } => "MachineToleranceAxisConditionAxisUnsupported",
            Self::NoToleranceAxes => "MachineToleranceAxisNoToleranceAxes",
            Self::BehaviorCorrelationModelMismatch => {
                "MachineToleranceAxisBehaviorCorrelationModelMismatch"
            }
            Self::MissingToleranceSpec { .. } => "MachineToleranceAxisMissingToleranceSpec",
            Self::NonRandomTolerance { .. } => "MachineToleranceAxisNonRandomTolerance",
            Self::NonScalarTolerance { .. } => "MachineToleranceAxisNonScalarTolerance",
            Self::UnsupportedToleranceTarget { .. } => {
                "MachineToleranceAxisUnsupportedToleranceTarget"
            }
            Self::MissingManufacturingHistory { .. } => {
                "MachineToleranceAxisMissingManufacturingHistory"
            }
            Self::AxisNameMismatch { .. } => "MachineToleranceAxisNameMismatch",
            Self::Identity(_) => "MachineToleranceAxisIdentity",
        }
    }
}

impl fmt::Display for MachineToleranceAxisCrosswalkErrorV1 {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::GraphMismatch {
                behavior,
                manufacturing,
            } => write!(
                formatter,
                "behavior graph {behavior} differs from manufacturing graph {manufacturing}"
            ),
            Self::MissingDependence => {
                formatter.write_str("behavior has no random dependence declaration")
            }
            Self::IndependentDependence => formatter
                .write_str("an independent Machine-IR dependence cannot bind a correlated stack"),
            Self::AxisDimensionMismatch { behavior, stack } => write!(
                formatter,
                "Machine-IR has {behavior} correlated axes but the stack retains {stack} terms"
            ),
            Self::ConditionAxisUnsupported { index } => write!(
                formatter,
                "correlated position {index} is a random condition; v1 requires a tolerance-only model"
            ),
            Self::NoToleranceAxes => {
                formatter.write_str("correlated behavior contains no tolerance axes to bind")
            }
            Self::BehaviorCorrelationModelMismatch => formatter.write_str(
                "Machine-IR correlation reference differs from the correlated-stack model",
            ),
            Self::MissingToleranceSpec { index, tolerance } => write!(
                formatter,
                "correlated position {index} names missing tolerance {tolerance}"
            ),
            Self::NonRandomTolerance { index, tolerance } => write!(
                formatter,
                "correlated position {index} names non-random tolerance {tolerance}"
            ),
            Self::NonScalarTolerance { index, tolerance } => write!(
                formatter,
                "correlated position {index} names non-scalar tolerance {tolerance}"
            ),
            Self::UnsupportedToleranceTarget { index, tolerance } => write!(
                formatter,
                "correlated position {index} names tolerance {tolerance} without a body target"
            ),
            Self::MissingManufacturingHistory {
                index,
                tolerance,
                body,
            } => write!(
                formatter,
                "correlated position {index} names tolerance {tolerance} on body {body} without manufacturing history"
            ),
            Self::AxisNameMismatch {
                index,
                tolerance,
                stack_name,
            } => write!(
                formatter,
                "correlated position {index} binds tolerance {tolerance} but stack label is {stack_name}"
            ),
            Self::Identity(error) => write!(
                formatter,
                "Machine tolerance-axis crosswalk identity refused: {error}"
            ),
        }
    }
}

impl std::error::Error for MachineToleranceAxisCrosswalkErrorV1 {}

impl From<CanonicalError> for MachineToleranceAxisCrosswalkErrorV1 {
    fn from(error: CanonicalError) -> Self {
        Self::Identity(error)
    }
}

/// Canonical identity schema for one exact Machine/tolerance-stack crosswalk.
pub enum MachineToleranceAxisCrosswalkIdentitySchemaV1 {}

impl CanonicalSchema for MachineToleranceAxisCrosswalkIdentitySchemaV1 {
    const DOMAIN: &'static str = "org.frankensim.fs-ir.machine.tolerance-axis-crosswalk.v1";
    const NAME: &'static str = "admitted-machine-tolerance-axis-crosswalk";
    const VERSION: u32 = MACHINE_TOLERANCE_AXIS_CROSSWALK_SCHEMA_VERSION_V1;
    const CONTEXT: &'static str = "one Machine graph, behavior, manufacturing state, exact fs-toleralloc factor and positional term receipt, and canonical tolerance-axis association";
    const FIELDS: &'static [FieldSpec] = &[
        FieldSpec::required("crosswalk-schema-version", WireType::U64),
        FieldSpec::required("frankenscript-ir-version", WireType::U64),
        FieldSpec::required("machine-graph", WireType::Bytes),
        FieldSpec::required("machine-behavior", WireType::Bytes),
        FieldSpec::required("manufacturing-state", WireType::Bytes),
        FieldSpec::required("correlation-model", WireType::Bytes),
        FieldSpec::required("correlation-coordinate-link", WireType::Bytes),
        FieldSpec::required("stack-terms", WireType::OrderedBytes),
        FieldSpec::required("tolerance-axis-bindings", WireType::OrderedBytes),
        FieldSpec::required("first-order-moments", WireType::Bytes),
    ];
}

/// Strong semantic identity of an admitted correlated-tolerance crosswalk.
pub type MachineToleranceAxisCrosswalkIdV1 =
    ProblemSemanticId<MachineToleranceAxisCrosswalkIdentitySchemaV1>;

/// One Machine tolerance's exact position in the correlated factor/term order.
///
/// Positions are explicit even though v1 admits only tolerance-only correlated
/// groups, preserving the external factor's axis convention in the receipt.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MachineToleranceAxisBindingV1 {
    position: usize,
    tolerance: ToleranceId,
    body: BodyId,
}

impl MachineToleranceAxisBindingV1 {
    /// Zero-based position in both the Machine dependence and stack receipt.
    #[must_use]
    pub const fn position(&self) -> usize {
        self.position
    }

    /// Exact durable Machine-IR tolerance bound at that position.
    #[must_use]
    pub const fn tolerance(&self) -> &ToleranceId {
        &self.tolerance
    }

    /// Durable body whose process history supplies the as-built attachment.
    #[must_use]
    pub const fn body(&self) -> &BodyId {
        &self.body
    }
}

/// Non-forgeable positional binding between Machine-IR and `fs-toleralloc`.
#[derive(Debug, Clone, PartialEq)]
pub struct AdmittedMachineToleranceAxisCrosswalkV1 {
    graph: MachineGraphIdV1,
    behavior: MachineBehaviorIdV1,
    manufacturing_state: MachineManufacturingStateIdV1,
    correlation_coordinate_link: ManufacturingArtifactRefV1,
    tolerance_axes: Vec<MachineToleranceAxisBindingV1>,
    stack: CorrelatedStackReceipt,
    receipt: IdentityReceipt<MachineToleranceAxisCrosswalkIdV1>,
}

impl AdmittedMachineToleranceAxisCrosswalkV1 {
    /// Exact Machine graph shared by both source receipts.
    #[must_use]
    pub const fn graph(&self) -> MachineGraphIdV1 {
        self.graph
    }

    /// Exact admitted behavior whose canonical dependence order is bound.
    #[must_use]
    pub const fn behavior(&self) -> MachineBehaviorIdV1 {
        self.behavior
    }

    /// Exact admitted manufacturing-state receipt participating in the binding.
    #[must_use]
    pub const fn manufacturing_state(&self) -> MachineManufacturingStateIdV1 {
        self.manufacturing_state
    }

    /// Versioned caller assertion linking content-addressed manufacturing and
    /// semantic `fs-toleralloc` correlation coordinates.
    #[must_use]
    pub const fn correlation_coordinate_link(&self) -> &ManufacturingArtifactRefV1 {
        &self.correlation_coordinate_link
    }

    /// Tolerance IDs in the exact positional order of the retained factor.
    #[must_use]
    pub fn tolerance_axes(&self) -> &[MachineToleranceAxisBindingV1] {
        &self.tolerance_axes
    }

    /// Exact non-forgeable stack receipt, including model, terms, and moments.
    #[must_use]
    pub const fn stack_receipt(&self) -> &CorrelatedStackReceipt {
        &self.stack
    }

    /// Domain-separated semantic identity of the complete crosswalk.
    #[must_use]
    pub const fn identity(&self) -> MachineToleranceAxisCrosswalkIdV1 {
        self.receipt.id()
    }

    /// Complete canonical preimage for collision adjudication.
    #[must_use]
    pub const fn identity_receipt(&self) -> IdentityReceipt<MachineToleranceAxisCrosswalkIdV1> {
        self.receipt
    }

    /// Body-attached lineage dependents for every bound tolerance axis.
    ///
    /// Generic Machine lineage records move or invalidate these axis
    /// attachments only. They do not mutate this receipt or mint a successor
    /// crosswalk: any graph/body/behavior/manufacturing transition requires
    /// explicit readmission against new admitted endpoints. If only one of
    /// several bodies is ambiguous, the lineage invalidation lists only the
    /// affected axis attachments, while the old aggregate receipt remains
    /// bound to its original endpoints and is not reusable as a successor.
    #[must_use]
    pub fn lineage_dependents(&self) -> Vec<DependentBinding> {
        self.tolerance_axes
            .iter()
            .map(|axis| {
                DependentBinding::new(
                    DependentKind::ManufacturingToleranceAxis,
                    axis.tolerance.canonical_key(),
                    MachineElementId::Body(axis.body.clone()),
                )
                .expect("admitted tolerance keys remain canonical")
            })
            .collect()
    }
}

impl AdmittedMachineManufacturingStateV1 {
    /// Bind one admitted behavior's canonical random-tolerance order to an
    /// exact `fs-toleralloc` correlated-stack receipt.
    ///
    /// # Errors
    /// Refuses graph/model gaps, absent or independent dependence semantics,
    /// absent tolerance axes, mixed random conditions, non-scalar/non-body
    /// tolerances, missing body process history, dimension gaps, positional
    /// label gaps, or canonical identity publication failure.
    pub fn bind_correlated_tolerance_axes(
        &self,
        behavior: &AdmittedMachineBehavior,
        stack: &CorrelatedStackReceipt,
        correlation_coordinate_link: ManufacturingArtifactRefV1,
    ) -> Result<AdmittedMachineToleranceAxisCrosswalkV1, MachineToleranceAxisCrosswalkErrorV1> {
        bind_correlated_tolerance_axes(behavior, self, stack, correlation_coordinate_link)
    }
}

/// Bind one admitted Machine behavior and manufacturing state to an exact
/// correlated-stack receipt.
///
/// # Errors
/// Returns a structured refusal without publishing an identity whenever graph,
/// model, dimension, axis-kind, or positional-name closure is absent.
#[allow(clippy::too_many_lines)] // One ordered fail-closed admission and identity pass.
pub fn bind_correlated_tolerance_axes(
    behavior: &AdmittedMachineBehavior,
    manufacturing: &AdmittedMachineManufacturingStateV1,
    stack: &CorrelatedStackReceipt,
    correlation_coordinate_link: ManufacturingArtifactRefV1,
) -> Result<AdmittedMachineToleranceAxisCrosswalkV1, MachineToleranceAxisCrosswalkErrorV1> {
    if behavior.base_graph() != manufacturing.graph() {
        return Err(MachineToleranceAxisCrosswalkErrorV1::GraphMismatch {
            behavior: behavior.base_graph(),
            manufacturing: manufacturing.graph(),
        });
    }

    let dependence = behavior
        .dependences()
        .first()
        .ok_or(MachineToleranceAxisCrosswalkErrorV1::MissingDependence)?;
    let correlation = match &dependence.model {
        DependenceModel::Independent => {
            return Err(MachineToleranceAxisCrosswalkErrorV1::IndependentDependence);
        }
        DependenceModel::Correlated(model) => model,
    };

    if dependence.members.len() != stack.terms().len() {
        return Err(
            MachineToleranceAxisCrosswalkErrorV1::AxisDimensionMismatch {
                behavior: dependence.members.len(),
                stack: stack.terms().len(),
            },
        );
    }

    let mut tolerance_axes = Vec::with_capacity(dependence.members.len());
    for (index, member) in dependence.members.iter().enumerate() {
        match member {
            DependenceMember::Condition(_) => {
                return Err(
                    MachineToleranceAxisCrosswalkErrorV1::ConditionAxisUnsupported { index },
                );
            }
            DependenceMember::Tolerance(tolerance) => {
                let Some(specification) = behavior
                    .tolerances()
                    .iter()
                    .find(|candidate| candidate.id == *tolerance)
                else {
                    return Err(MachineToleranceAxisCrosswalkErrorV1::MissingToleranceSpec {
                        index,
                        tolerance: tolerance.clone(),
                    });
                };
                if !matches!(&specification.semantics, ToleranceSemantics::Random { .. }) {
                    return Err(MachineToleranceAxisCrosswalkErrorV1::NonRandomTolerance {
                        index,
                        tolerance: tolerance.clone(),
                    });
                }
                if specification.shape != TerminalShape::Scalar {
                    return Err(MachineToleranceAxisCrosswalkErrorV1::NonScalarTolerance {
                        index,
                        tolerance: tolerance.clone(),
                    });
                }
                let ToleranceTarget::Element(MachineElementId::Body(body)) = &specification.target
                else {
                    return Err(
                        MachineToleranceAxisCrosswalkErrorV1::UnsupportedToleranceTarget {
                            index,
                            tolerance: tolerance.clone(),
                        },
                    );
                };
                if !manufacturing
                    .process_steps()
                    .iter()
                    .any(|step| step.body() == body)
                {
                    return Err(
                        MachineToleranceAxisCrosswalkErrorV1::MissingManufacturingHistory {
                            index,
                            tolerance: tolerance.clone(),
                            body: body.clone(),
                        },
                    );
                }
                tolerance_axes.push(MachineToleranceAxisBindingV1 {
                    position: index,
                    tolerance: tolerance.clone(),
                    body: body.clone(),
                });
            }
        }
    }
    if tolerance_axes.is_empty() {
        return Err(MachineToleranceAxisCrosswalkErrorV1::NoToleranceAxes);
    }
    if !behavior_model_matches(correlation, stack.model()) {
        return Err(MachineToleranceAxisCrosswalkErrorV1::BehaviorCorrelationModelMismatch);
    }
    for binding in &tolerance_axes {
        let term = &stack.terms()[binding.position];
        if binding.tolerance.canonical_key() != term.name {
            return Err(MachineToleranceAxisCrosswalkErrorV1::AxisNameMismatch {
                index: binding.position,
                tolerance: binding.tolerance.clone(),
                stack_name: term.name.clone(),
            });
        }
    }

    let model_row = correlation_model_row(stack.model());
    let coordinate_link_row = correlation_coordinate_link.canonical_row();
    let term_rows: Vec<Vec<u8>> = stack.terms().iter().map(stack_term_row).collect();
    let binding_rows: Vec<Vec<u8>> = tolerance_axes
        .iter()
        .map(tolerance_axis_binding_row)
        .collect();
    let moments_row = stack_moments_row(stack);
    let graph = behavior.base_graph();
    let behavior_id = behavior.identity();
    let manufacturing_state = manufacturing.identity();
    let receipt = CanonicalEncoder::<MachineToleranceAxisCrosswalkIdV1, _>::new(
        TOLERANCE_AXIS_CROSSWALK_IDENTITY_LIMITS,
        NeverCancel,
    )?
    .u64(
        Field::new(0, "crosswalk-schema-version"),
        u64::from(MACHINE_TOLERANCE_AXIS_CROSSWALK_SCHEMA_VERSION_V1),
    )?
    .u64(
        Field::new(1, "frankenscript-ir-version"),
        u64::from(IR_VERSION),
    )?
    .bytes(Field::new(2, "machine-graph"), graph.as_bytes())?
    .bytes(Field::new(3, "machine-behavior"), behavior_id.as_bytes())?
    .bytes(
        Field::new(4, "manufacturing-state"),
        manufacturing_state.as_bytes(),
    )?
    .bytes(Field::new(5, "correlation-model"), &model_row)?
    .bytes(
        Field::new(6, "correlation-coordinate-link"),
        &coordinate_link_row,
    )?
    .ordered_bytes(
        Field::new(7, "stack-terms"),
        term_rows.len() as u64,
        term_rows.iter().map(Vec::as_slice),
    )?
    .ordered_bytes(
        Field::new(8, "tolerance-axis-bindings"),
        binding_rows.len() as u64,
        binding_rows.iter().map(Vec::as_slice),
    )?
    .bytes(Field::new(9, "first-order-moments"), &moments_row)?
    .finish()?;

    Ok(AdmittedMachineToleranceAxisCrosswalkV1 {
        graph,
        behavior: behavior_id,
        manufacturing_state,
        correlation_coordinate_link,
        tolerance_axes,
        stack: stack.clone(),
        receipt,
    })
}

fn behavior_model_matches(
    reference: &CorrelationModelRef,
    model: &AdmittedCorrelationModel,
) -> bool {
    reference.namespace() == model.namespace()
        && reference.schema_version() == model.schema_version()
        && reference.semantic_digest() == model.semantic_digest()
}

fn correlation_model_row(model: &AdmittedCorrelationModel) -> Vec<u8> {
    let mut row = Vec::with_capacity(96 + model.lower_factor().len() * 8);
    append_bytes(&mut row, model.namespace().as_bytes());
    row.extend_from_slice(&model.schema_version().get().to_le_bytes());
    row.extend_from_slice(&model.semantic_digest());
    row.extend_from_slice(&(model.dimension() as u64).to_le_bytes());
    row.extend_from_slice(&model.max_row_norm_defect().to_bits().to_le_bytes());
    row.extend_from_slice(&(model.lower_factor().len() as u64).to_le_bytes());
    for value in model.lower_factor() {
        row.extend_from_slice(&value.to_bits().to_le_bytes());
    }
    row
}

fn stack_term_row(term: &CorrelatedStackTerm) -> Vec<u8> {
    let mut row = Vec::with_capacity(term.name.len() + 32);
    append_bytes(&mut row, term.name.as_bytes());
    row.extend_from_slice(&term.signed_sensitivity.to_bits().to_le_bytes());
    row.push(color_rank_tag(term.sensitivity_color));
    row.extend_from_slice(&term.standard_deviation.to_bits().to_le_bytes());
    row
}

fn tolerance_axis_binding_row(binding: &MachineToleranceAxisBindingV1) -> Vec<u8> {
    let mut row = Vec::with_capacity(
        binding.tolerance.canonical_key().len() + binding.body.canonical_key().len() + 104,
    );
    row.extend_from_slice(&(binding.position as u64).to_le_bytes());
    row.extend_from_slice(binding.tolerance.identity().as_bytes());
    append_bytes(&mut row, binding.tolerance.canonical_key().as_bytes());
    row.extend_from_slice(binding.body.identity().as_bytes());
    append_bytes(&mut row, binding.body.canonical_key().as_bytes());
    row
}

fn color_rank_tag(rank: ColorRank) -> u8 {
    match rank {
        ColorRank::Estimated => 1,
        ColorRank::Validated => 2,
        ColorRank::Verified => 3,
    }
}

fn stack_moments_row(stack: &CorrelatedStackReceipt) -> Vec<u8> {
    let values = [
        stack.independent_standard_deviation(),
        stack.independent_variance(),
        stack.correlated_standard_deviation(),
        stack.correlated_variance(),
        stack.correlation_variance_delta(),
    ];
    let mut row = Vec::with_capacity(values.len() * 8);
    for value in values {
        row.extend_from_slice(&value.to_bits().to_le_bytes());
    }
    row
}
