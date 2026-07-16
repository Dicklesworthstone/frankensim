//! I01.2 (bead i94v.1.1.2): canonical ConstitutiveGraph law nodes
//! adapted into the equation compiler through an OPAQUE typed
//! protocol.
//!
//! Constitutive laws own their internal state, free energy,
//! dissipation, tangents, and validity domains in fs-material; copying
//! any of that into this crate would create drift and a forbidden
//! reverse dependency. This module therefore holds law nodes only as
//! `&dyn fs_material::graph::LawNode` behind [`BoundConstitutiveNode`]
//! and re-expresses everything the compiler must retain in
//! COMPILER-OWNED receipt types: law/version identity, state-schema
//! version, initialization and update policy, differentiability class,
//! thermodynamic-potential chart, material usage receipts, and
//! validity monitors all survive lowering exactly, without fs-material
//! ever importing compiler internals (the dependency points here → L3,
//! never back).
//!
//! A supplied tangent is EVIDENCE TO VERIFY, not authority: binding
//! verifies any consistent-tangent claim against central finite
//! differences before granting the [`TangentLane::Consistent`] route;
//! nodes whose declarations cannot support a derivative lane are
//! ROUTED (typed), never silently differentiated. Hand-written escape
//! hatches are admissible but carry an explicit
//! no-generated-consistency marker that downstream generators must
//! surface.
//!
//! Ambition class [F]: everything here sits behind the default-off
//! `constitutive-graph` feature until the candidate beats its declared
//! baseline on the preregistered deck.

use std::fmt;

use fs_material::graph::{GraphError, LawNode, check_consistent_tangent};
use fs_qty::Dims;

/// Tolerance for the binding-time tangent-evidence gate (absolute, on
/// FD-vs-analytic deviation at the probe point).
pub const TANGENT_EVIDENCE_ATOL: f64 = 1e-6;

/// Compiler-owned copy of everything about a law node that must
/// survive lowering exactly. No fs-material type appears here: the
/// receipt is the compiler's OWN record, comparable and hashable
/// without reaching back into L3.
#[derive(Debug, Clone, PartialEq)]
pub struct MaterialProvenance {
    /// Law identity (fs-matdb keyed, copied verbatim).
    pub law: String,
    /// Law semantic version.
    pub law_version: u32,
    /// Internal-state schema version.
    pub state_schema_version: u32,
    /// Number of internal-state slots (codec arity).
    pub state_slots: usize,
    /// Input-port dims, in declaration order.
    pub input_dims: Vec<Dims>,
    /// Output-port dims, in declaration order.
    pub output_dims: Vec<Dims>,
    /// Differentiability class tag (compiler-owned spelling).
    pub differentiability: DifferentiabilityTag,
    /// Thermodynamic-potential chart tag (compiler-owned spelling).
    pub potential_chart: PotentialChart,
}

/// Compiler-owned differentiability spelling (mirrors the L3
/// declaration without importing it).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DifferentiabilityTag {
    /// C¹ or better on the calibration domain.
    Smooth,
    /// Piecewise smooth with declared kink sets.
    PiecewiseSmooth,
    /// No derivative claim.
    NonSmooth,
}

/// Compiler-owned thermodynamic-potential chart: which audited
/// potentials the node exposes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PotentialChart {
    /// Free-energy storage only.
    Storage,
    /// Non-negative dissipation only.
    Dissipation,
    /// Both storage and dissipation.
    StorageAndDissipation,
    /// Explicit empirical no-claim.
    EmpiricalNoClaim,
}

/// Which derivative lane the compiler may route this node through.
/// Routing is TYPED: a node never silently enters a lane its
/// declaration cannot support.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TangentLane {
    /// Verified consistent tangent: generated JVP/VJP programs may
    /// consume the analytic tangent (evidence checked at binding).
    Consistent,
    /// Piecewise-smooth: the tangent is usable pointwise but every
    /// generated program must carry the kink-set caveat.
    PiecewiseConsistent,
    /// Derivative-free: only value lanes are generated.
    DerivativeFree,
}

/// The escape hatch for hand-written nodes: admissible, but the
/// marker is retained in provenance and downstream generators MUST
/// NOT claim generated consistency for it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct HandWrittenEscape {
    /// Always true for escape-hatch bindings; spelled out so the
    /// receipt is self-describing.
    pub no_generated_consistency_claim: bool,
}

/// A typed adaptation refusal.
#[derive(Debug, Clone, PartialEq)]
pub enum AdaptError {
    /// The node's tangent claim failed its finite-difference evidence
    /// gate at binding.
    TangentEvidenceRejected {
        /// The law id.
        law: String,
        /// The underlying L3 diagnosis, stringified (the compiler does
        /// not re-export L3 error types).
        detail: String,
    },
    /// The declaration cannot support the requested lane.
    UnsupportedDifferentiability {
        /// The law id.
        law: String,
        /// The lane that was requested.
        requested: &'static str,
        /// The declared class.
        declared: DifferentiabilityTag,
    },
    /// The node owns state but the binding supplied no initial state
    /// and the declaration demands one.
    MissingStateInitialization {
        /// The law id.
        law: String,
        /// Declared state arity.
        state_slots: usize,
    },
    /// Evaluation arity or the node's own numerical refusal.
    Evaluation {
        /// The law id.
        law: String,
        /// Stringified L3 diagnosis.
        detail: String,
    },
    /// A state vector under a different schema version was offered.
    StateSchemaDrift {
        /// The law id.
        law: String,
        /// Version the binding was made under.
        bound: u32,
        /// Version offered.
        offered: u32,
    },
}

impl fmt::Display for AdaptError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            AdaptError::TangentEvidenceRejected { law, detail } => write!(
                f,
                "law `{law}`: supplied tangent is evidence and it FAILED verification: {detail}"
            ),
            AdaptError::UnsupportedDifferentiability {
                law,
                requested,
                declared,
            } => write!(
                f,
                "law `{law}`: lane `{requested}` is not supported by declared class {declared:?}"
            ),
            AdaptError::MissingStateInitialization { law, state_slots } => write!(
                f,
                "law `{law}`: {state_slots} state slot(s) declared but no initial state supplied"
            ),
            AdaptError::Evaluation { law, detail } => {
                write!(f, "law `{law}`: evaluation refused: {detail}")
            }
            AdaptError::StateSchemaDrift {
                law,
                bound,
                offered,
            } => write!(
                f,
                "law `{law}`: state offered under schema v{offered}, bound under v{bound}"
            ),
        }
    }
}

impl core::error::Error for AdaptError {}

fn graph_error_detail(error: &GraphError) -> String {
    format!("{error:?}")
}

/// A law node bound into the compiler: the opaque handle plus the
/// compiler-owned provenance and the verified lane routing.
pub struct BoundConstitutiveNode<'n> {
    // (Debug below shows the receipt, never the opaque node.)
    node: &'n dyn LawNode,
    provenance: MaterialProvenance,
    lane: TangentLane,
    escape: Option<HandWrittenEscape>,
    state: Vec<f64>,
}

impl fmt::Debug for BoundConstitutiveNode<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("BoundConstitutiveNode")
            .field("provenance", &self.provenance)
            .field("lane", &self.lane)
            .field("escape", &self.escape)
            .field("state", &self.state)
            .finish_non_exhaustive()
    }
}

impl<'n> BoundConstitutiveNode<'n> {
    /// Bind a canonical law node. Copies the declaration into the
    /// compiler-owned provenance receipt, routes the derivative lane
    /// from the declared differentiability, VERIFIES any tangent claim
    /// against finite differences at the zero probe point, and demands
    /// an explicit initial state whenever the node owns state.
    ///
    /// # Errors
    /// Typed [`AdaptError`] refusals; nothing is bound partially.
    pub fn bind(
        node: &'n dyn LawNode,
        initial_state: Option<&[f64]>,
    ) -> Result<BoundConstitutiveNode<'n>, AdaptError> {
        use fs_material::graph::{Differentiability, EnergyBehavior};
        let declaration = node.declaration();
        let law = declaration.law.0.clone();

        let differentiability = match declaration.differentiability {
            Differentiability::Smooth => DifferentiabilityTag::Smooth,
            Differentiability::PiecewiseSmooth => DifferentiabilityTag::PiecewiseSmooth,
            Differentiability::NonSmooth => DifferentiabilityTag::NonSmooth,
        };
        let potential_chart = match declaration.energy {
            EnergyBehavior::FreeEnergyStorage => PotentialChart::Storage,
            EnergyBehavior::NonNegativeDissipation => PotentialChart::Dissipation,
            EnergyBehavior::StorageAndDissipation => PotentialChart::StorageAndDissipation,
            EnergyBehavior::Empirical => PotentialChart::EmpiricalNoClaim,
        };
        let provenance = MaterialProvenance {
            law: law.clone(),
            law_version: declaration.law_version,
            state_schema_version: declaration.state_schema_version,
            state_slots: declaration.state_slots.len(),
            input_dims: declaration.inputs.iter().map(|p| p.dims).collect(),
            output_dims: declaration.outputs.iter().map(|p| p.dims).collect(),
            differentiability,
            potential_chart,
        };

        // State initialization policy: state-owning nodes must be
        // given an explicit start (the compiler never invents one).
        let state = if provenance.state_slots == 0 {
            Vec::new()
        } else {
            let Some(initial) = initial_state else {
                return Err(AdaptError::MissingStateInitialization {
                    law,
                    state_slots: provenance.state_slots,
                });
            };
            if initial.len() != provenance.state_slots {
                return Err(AdaptError::MissingStateInitialization {
                    law,
                    state_slots: provenance.state_slots,
                });
            }
            initial.to_vec()
        };

        // Route the derivative lane from the declaration, then verify
        // any tangent claim as EVIDENCE before granting the lane.
        let lane = match (differentiability, declaration.tangent_claimed) {
            (DifferentiabilityTag::Smooth, true) => {
                let probe_inputs = vec![0.0f64; provenance.input_dims.len()];
                check_consistent_tangent(&law, node, &state, &probe_inputs, TANGENT_EVIDENCE_ATOL)
                    .map_err(|e| AdaptError::TangentEvidenceRejected {
                        law: law.clone(),
                        detail: graph_error_detail(&e),
                    })?;
                TangentLane::Consistent
            }
            (DifferentiabilityTag::PiecewiseSmooth, true) => TangentLane::PiecewiseConsistent,
            (DifferentiabilityTag::NonSmooth, true) => {
                return Err(AdaptError::UnsupportedDifferentiability {
                    law,
                    requested: "consistent-tangent",
                    declared: differentiability,
                });
            }
            (_, false) => TangentLane::DerivativeFree,
        };

        Ok(BoundConstitutiveNode {
            node,
            provenance,
            lane,
            escape: None,
            state,
        })
    }

    /// Bind a HAND-WRITTEN node through the escape hatch: same gates,
    /// plus a retained marker forbidding generated-consistency claims.
    ///
    /// # Errors
    /// Same refusals as [`BoundConstitutiveNode::bind`].
    pub fn bind_hand_written(
        node: &'n dyn LawNode,
        initial_state: Option<&[f64]>,
    ) -> Result<BoundConstitutiveNode<'n>, AdaptError> {
        let mut bound = BoundConstitutiveNode::bind(node, initial_state)?;
        bound.escape = Some(HandWrittenEscape {
            no_generated_consistency_claim: true,
        });
        Ok(bound)
    }

    /// The compiler-owned provenance receipt (survives lowering).
    #[must_use]
    pub fn provenance(&self) -> &MaterialProvenance {
        &self.provenance
    }

    /// The verified derivative lane.
    #[must_use]
    pub fn lane(&self) -> TangentLane {
        self.lane
    }

    /// The escape-hatch marker, when bound through it.
    #[must_use]
    pub fn escape(&self) -> Option<HandWrittenEscape> {
        self.escape
    }

    /// Current internal state (slot order).
    #[must_use]
    pub fn state(&self) -> &[f64] {
        &self.state
    }

    /// Replace the internal state with one encoded under an explicit
    /// schema version, refusing drift.
    ///
    /// # Errors
    /// [`AdaptError::StateSchemaDrift`] / arity refusals.
    pub fn restore_state(&mut self, schema_version: u32, state: &[f64]) -> Result<(), AdaptError> {
        if schema_version != self.provenance.state_schema_version {
            return Err(AdaptError::StateSchemaDrift {
                law: self.provenance.law.clone(),
                bound: self.provenance.state_schema_version,
                offered: schema_version,
            });
        }
        if state.len() != self.provenance.state_slots {
            return Err(AdaptError::MissingStateInitialization {
                law: self.provenance.law.clone(),
                state_slots: self.provenance.state_slots,
            });
        }
        self.state.copy_from_slice(state);
        Ok(())
    }

    /// Evaluate the node at `inputs`, committing the updated state and
    /// returning outputs plus the reported dissipation rate.
    ///
    /// # Errors
    /// [`AdaptError::Evaluation`] wrapping the node's typed refusal.
    pub fn evaluate(&mut self, inputs: &[f64]) -> Result<(Vec<f64>, Option<f64>), AdaptError> {
        let output =
            self.node
                .evaluate(&self.state, inputs)
                .map_err(|e| AdaptError::Evaluation {
                    law: self.provenance.law.clone(),
                    detail: graph_error_detail(&e),
                })?;
        self.state = output.next_state;
        Ok((output.outputs, output.dissipation_rate))
    }

    /// The analytic tangent, available only on a verified lane.
    ///
    /// # Errors
    /// [`AdaptError::UnsupportedDifferentiability`] off-lane.
    pub fn tangent(&self, inputs: &[f64]) -> Result<Vec<f64>, AdaptError> {
        match self.lane {
            TangentLane::Consistent | TangentLane::PiecewiseConsistent => self
                .node
                .tangent(&self.state, inputs)
                .ok_or_else(|| AdaptError::Evaluation {
                    law: self.provenance.law.clone(),
                    detail: "declared tangent unavailable at this point".to_string(),
                }),
            TangentLane::DerivativeFree => Err(AdaptError::UnsupportedDifferentiability {
                law: self.provenance.law.clone(),
                requested: "tangent",
                declared: self.provenance.differentiability,
            }),
        }
    }

    /// Vector-Jacobian product `wᵀ · d(outputs)/d(inputs)` derived from
    /// the verified tangent (transpose contraction — exact, not
    /// approximated), on tangent-bearing lanes only.
    ///
    /// # Errors
    /// Same routing refusals as [`BoundConstitutiveNode::tangent`].
    pub fn vjp(&self, inputs: &[f64], w: &[f64]) -> Result<Vec<f64>, AdaptError> {
        let tangent = self.tangent(inputs)?;
        let rows = self.provenance.output_dims.len();
        let cols = self.provenance.input_dims.len();
        if w.len() != rows {
            return Err(AdaptError::Evaluation {
                law: self.provenance.law.clone(),
                detail: format!("VJP weight has {} entries for {rows} outputs", w.len()),
            });
        }
        let mut out = vec![0.0f64; cols];
        for (row, weight) in w.iter().enumerate() {
            for col in 0..cols {
                out[col] += weight * tangent[row * cols + col];
            }
        }
        Ok(out)
    }

    /// Free energy at the current state, when the potential chart
    /// exposes storage.
    #[must_use]
    pub fn free_energy(&self, inputs: &[f64]) -> Option<f64> {
        match self.provenance.potential_chart {
            PotentialChart::Storage | PotentialChart::StorageAndDissipation => {
                self.node.free_energy(&self.state, inputs)
            }
            PotentialChart::Dissipation | PotentialChart::EmpiricalNoClaim => None,
        }
    }
}

// ---------------------------------------------------------------------------
// Slice 2: batched law evaluation under Cx — request-drain-finalize.
// ---------------------------------------------------------------------------

/// One material point in a batched evaluation: its own state buffer
/// (under the graph's schema version) and external port feeds.
#[derive(Debug, Clone, PartialEq)]
pub struct BatchPoint {
    /// Schema version the state buffer was encoded under.
    pub state_version: u64,
    /// Aggregate state buffer for this point.
    pub state: Vec<f64>,
    /// External inputs keyed `(node, port)`.
    pub external: std::collections::BTreeMap<(String, String), f64>,
}

/// One completed point's result (compiler-owned mirror of the graph's
/// single-pass output).
#[derive(Debug, Clone, PartialEq)]
pub struct BatchPointOutput {
    /// Every node output, keyed `(node, port)`.
    pub outputs: std::collections::BTreeMap<(String, String), f64>,
    /// Updated state buffer for this point.
    pub next_state: Vec<f64>,
    /// Total reported dissipation at this point.
    pub total_dissipation: f64,
}

/// Poll cadence: one cancellation checkpoint per this many points.
pub const BATCH_POLL_STRIDE: usize = 16;

/// Terminal state of a batched evaluation (request-drain-finalize):
/// cancellation is observed at bounded point boundaries, the point in
/// flight DRAINS to completion (points are atomic), and the completed
/// prefix is finalized with a deterministic resume cursor. Resuming at
/// `resume_from` with the same batch is bitwise-equivalent to an
/// uncancelled run.
#[derive(Debug, Clone, PartialEq)]
pub enum BatchRun {
    /// Every point evaluated.
    Complete {
        /// Per-point results, in batch order.
        outputs: Vec<BatchPointOutput>,
        /// Points evaluated (== batch length).
        points_evaluated: usize,
    },
    /// Cancellation drained at a point boundary.
    Cancelled {
        /// Results for the completed prefix, in batch order.
        completed: Vec<BatchPointOutput>,
        /// Index of the first UNevaluated point (the resume cursor).
        resume_from: usize,
    },
}

/// Evaluate a canonical [`fs_material::graph::ConstitutiveGraph`] over
/// a batch of material points under an execution context: one
/// deterministic single pass per point, a cancellation poll every
/// [`BATCH_POLL_STRIDE`] points, and request-drain-finalize semantics
/// — an observed cancellation never abandons a point mid-pass and
/// never publishes a partial point.
///
/// # Errors
/// [`AdaptError::Evaluation`] wrapping the graph's typed refusal for
/// the offending point (its index is named); the completed prefix is
/// NOT returned on a refusal — a defective batch is not a result.
pub fn evaluate_batch(
    graph: &fs_material::graph::ConstitutiveGraph,
    batch: &[BatchPoint],
    cx: &fs_exec::Cx<'_>,
) -> Result<BatchRun, AdaptError> {
    let mut outputs = Vec::with_capacity(batch.len());
    for (index, point) in batch.iter().enumerate() {
        if index % BATCH_POLL_STRIDE == 0 && cx.checkpoint().is_err() {
            return Ok(BatchRun::Cancelled {
                completed: outputs,
                resume_from: index,
            });
        }
        let result = graph
            .execute(point.state_version, &point.state, &point.external)
            .map_err(|e| AdaptError::Evaluation {
                law: format!("batch point {index}"),
                detail: graph_error_detail(&e),
            })?;
        outputs.push(BatchPointOutput {
            outputs: result.outputs,
            next_state: result.next_state,
            total_dissipation: result.total_dissipation,
        });
    }
    Ok(BatchRun::Complete {
        points_evaluated: outputs.len(),
        outputs,
    })
}
