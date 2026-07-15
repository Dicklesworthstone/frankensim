//! DWR GOAL-ORIENTED ACCEPT TEST (addendum Proposal 9, bead lmp4.4;
//! [F] — behind the `dwr-accept` feature): dual-weighted-residual
//! estimates target the QUERY's actual quantity of interest, which
//! sharpens the accept test enormously — but DWR constants are NOT
//! guaranteed, so a DWR-only accept carries ESTIMATED color. Promotion
//! to VERIFIED additionally requires a typed proof that the supplied dual is
//! the exact dual of THIS query. The current v0 API cannot express that proof,
//! so even an independently reverified Cauchy–Schwarz energy-product bracket
//! remains an ESTIMATED diagnostic and never promotes or vetoes acceptance.
//! False certification is worse than temporarily losing the promotion path.

use fs_evidence::Color;
use fs_exec::Cx;
use fs_ivl::expansion::two_diff;
use fs_verify::estimator::{
    EstimatorFamily, VERIFIER_POLL_POLICY_VERSION, VERIFIER_POLL_STRIDE_WORK_UNITS,
    VERIFIER_WORK_PLAN_VERSION, VerifierCheckpointKind, VerifierPhase, VerifierProgress,
    VerifierWorkPlan, verify_with_checkpoint,
};
use fs_verify::fem1d::{
    MAX_FEM1D_MESH_NODES, MAX_FEM1D_POLY_COEFFICIENTS, MAX_FEM1D_PROBLEM_CANONICAL_IDENTITY_BYTES,
    MmsProblem, gauss5,
};
use fs_verify::interval::up;

/// Maximum mesh nodes admitted by the in-process rigorous bracket verifier.
pub const MAX_BRACKET_MESH_NODES: usize = MAX_FEM1D_MESH_NODES;
/// Maximum coarse mesh/candidate nodes admitted by the DWR execution path.
pub const MAX_DWR_MESH_NODES: usize = MAX_FEM1D_MESH_NODES;
/// Maximum manufactured-solution polynomial coefficients admitted by DWR.
pub const MAX_DWR_POLY_COEFFICIENTS: usize = MAX_FEM1D_POLY_COEFFICIENTS;
/// Maximum bounded logical work units admitted by one DWR workflow.
pub const MAX_DWR_WORK_UNITS: usize = 100_000_000;
/// Maximum UTF-8 bytes admitted in a QoI provenance label.
pub const MAX_DWR_QOI_BYTES: usize = 4_096;
/// Version of the complete checked work-plan encoding bound into evidence.
pub const DWR_WORK_PLAN_VERSION: u32 = 2;
/// Version of the cancellation-poll policy bound into evidence.
pub const DWR_POLL_POLICY_VERSION: u32 = 3;
/// Maximum bounded items processed between cancellation checkpoints.
pub const DWR_POLL_STRIDE_ITEMS: usize = 256;
const DWR_POLL_STRIDE_IDENTITY: u64 = DWR_POLL_STRIDE_ITEMS as u64;
/// Version of the retained DWR execution/evidence identity.
pub const DWR_EVIDENCE_IDENTITY_VERSION: u32 = 5;
const MAX_DWR_REFINED_NODES: usize = MAX_DWR_MESH_NODES * 2 - 1;

const DWR_OUTPUT_IDENTITY_SCHEMA: &[u8] = b"fs-adjoint-dwr-output-identity-v5";
const DWR_BRACKET_IDENTITY_SCHEMA: &[u8] = b"fs-adjoint-dwr-bracket-identity-v5";
const DWR_ACCEPT_IDENTITY_SCHEMA: &[u8] = b"fs-adjoint-dwr-accept-identity-v5";

const DWR_INITIAL_PHASE: &str = "dwr.initial";
const DWR_VALIDATE_POLYNOMIAL_PHASE: &str = "dwr.validate-polynomial";
const DWR_VALIDATE_CANDIDATE_PHASE: &str = "dwr.validate-candidate";
const DWR_VALIDATE_MESH_PHASE: &str = "dwr.validate-mesh";
const DWR_VALIDATE_CELLS_PHASE: &str = "dwr.validate-cells";
const DWR_VALIDATE_FORCING_PHASE: &str = "dwr.validate-forcing";
const DWR_PRIMAL_PHASE: &str = "dwr.primal-integral";
const DWR_REFINE_PHASE: &str = "dwr.refine";
const DWR_DUAL_ASSEMBLY_PHASE: &str = "dwr.dual-assembly";
const DWR_THOMAS_VALIDATE_PHASE: &str = "dwr.thomas-validate";
const DWR_THOMAS_FORWARD_PHASE: &str = "dwr.thomas-forward";
const DWR_THOMAS_BACK_PHASE: &str = "dwr.thomas-back";
const DWR_DUAL_PUBLICATION_PHASE: &str = "dwr.dual-publication";
const DWR_RESIDUAL_PHASE: &str = "dwr.residual";
const DWR_OUTPUT_VALIDATE_PHASE: &str = "dwr.output-validation";
const DWR_IDENTITY_PHASE: &str = "dwr.identity";
const DWR_PUBLICATION_PHASE: &str = "dwr.publish";
const BRACKET_INITIAL_PHASE: &str = "dwr-bracket.initial";
const BRACKET_PRIMAL_VALIDATE_PHASE: &str = "dwr-bracket.primal-validation";
const BRACKET_PRIMAL_VERIFY_PHASE: &str = "dwr-bracket.primal-verifier";
const BRACKET_DUAL_VALIDATE_PHASE: &str = "dwr-bracket.dual-validation";
const BRACKET_DUAL_VERIFY_PHASE: &str = "dwr-bracket.dual-verifier";
const BRACKET_IDENTITY_PHASE: &str = "dwr-bracket.identity";
const BRACKET_PUBLICATION_PHASE: &str = "dwr-bracket.publish";
const ACCEPT_INITIAL_PHASE: &str = "dwr-accept.initial";
const ACCEPT_IDENTITY_PHASE: &str = "dwr-accept.identity";
const ACCEPT_PUBLICATION_PHASE: &str = "dwr-accept.publish";

/// A retained, domain-separated identity for a DWR execution or decision.
///
/// The root binds the execution mode, all budget fields, the complete
/// [`fs_exec::StreamKey`], the workflow's complete checked work plan, and the
/// versioned fixed-stride polling policy in addition to its scientific inputs
/// and result.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct DwrEvidenceIdentity(fs_blake3::ContentHash);

impl DwrEvidenceIdentity {
    /// Raw 32-byte retained root.
    #[must_use]
    pub fn as_bytes(&self) -> &[u8; 32] {
        self.0.as_bytes()
    }

    /// Lowercase hexadecimal retained root.
    #[must_use]
    pub fn to_hex(self) -> String {
        self.0.to_hex()
    }
}

impl core::fmt::Display for DwrEvidenceIdentity {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        self.0.fmt(formatter)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct WorkProgress {
    completed_work_units: u128,
    planned_work_units: u128,
    polls_remaining: u32,
}

impl WorkProgress {
    fn new(planned_work_units: u128, cx: &Cx<'_>) -> Self {
        Self {
            completed_work_units: 0,
            planned_work_units,
            polls_remaining: cx.budget().poll_quota,
        }
    }

    fn advance(&mut self, units: u128) -> Result<(), DwrError> {
        self.completed_work_units = self
            .completed_work_units
            .checked_add(units)
            .filter(|completed| *completed <= self.planned_work_units)
            .ok_or(DwrError::WorkPlanOverflow)?;
        Ok(())
    }

    fn advance_bracket(&mut self, units: u128) -> Result<(), BracketError> {
        self.completed_work_units = self
            .completed_work_units
            .checked_add(units)
            .filter(|completed| *completed <= self.planned_work_units)
            .ok_or(BracketError::WorkPlanOverflow)?;
        Ok(())
    }

    fn finish_dwr(&self, phase: &'static str) -> Result<(), DwrError> {
        if self.completed_work_units == self.planned_work_units {
            return Ok(());
        }
        Err(DwrError::WorkPlanMismatch {
            phase,
            completed_work_units: self.completed_work_units,
            planned_work_units: self.planned_work_units,
        })
    }

    fn finish_bracket(&self, phase: &'static str) -> Result<(), BracketError> {
        if self.completed_work_units == self.planned_work_units {
            return Ok(());
        }
        Err(BracketError::WorkPlanMismatch {
            phase,
            completed_work_units: self.completed_work_units,
            planned_work_units: self.planned_work_units,
        })
    }
}

fn dwr_checkpoint(
    phase: &'static str,
    progress: &mut WorkProgress,
    cx: &Cx<'_>,
) -> Result<(), DwrError> {
    if progress.polls_remaining == 0 {
        return Err(DwrError::Cancelled {
            phase,
            completed_work_units: progress.completed_work_units,
            planned_work_units: progress.planned_work_units,
        });
    }
    if progress.polls_remaining != u32::MAX {
        progress.polls_remaining -= 1;
    }
    cx.checkpoint().map_err(|_| DwrError::Cancelled {
        phase,
        completed_work_units: progress.completed_work_units,
        planned_work_units: progress.planned_work_units,
    })
}

fn bracket_checkpoint(
    phase: &'static str,
    progress: &mut WorkProgress,
    cx: &Cx<'_>,
) -> Result<(), BracketError> {
    if progress.polls_remaining == 0 {
        return Err(BracketError::Cancelled {
            phase,
            completed_work_units: progress.completed_work_units,
            planned_work_units: progress.planned_work_units,
        });
    }
    if progress.polls_remaining != u32::MAX {
        progress.polls_remaining -= 1;
    }
    cx.checkpoint().map_err(|_| BracketError::Cancelled {
        phase,
        completed_work_units: progress.completed_work_units,
        planned_work_units: progress.planned_work_units,
    })
}

fn poll_dwr_scan(
    index: usize,
    phase: &'static str,
    progress: &mut WorkProgress,
    cx: &Cx<'_>,
) -> Result<(), DwrError> {
    if index != 0 && index.is_multiple_of(DWR_POLL_STRIDE_ITEMS) {
        dwr_checkpoint(phase, progress, cx)?;
    }
    Ok(())
}

fn poll_bracket_scan(
    index: usize,
    phase: &'static str,
    progress: &mut WorkProgress,
    cx: &Cx<'_>,
) -> Result<(), BracketError> {
    if index != 0 && index.is_multiple_of(DWR_POLL_STRIDE_ITEMS) {
        bracket_checkpoint(phase, progress, cx)?;
    }
    Ok(())
}

fn hash_execution_header(
    hasher: &mut fs_blake3::Blake3,
    schema: &[u8],
    work_plan_fields: &[u128],
    cx: &Cx<'_>,
) {
    hasher.update(schema);
    hasher.update(&DWR_EVIDENCE_IDENTITY_VERSION.to_le_bytes());
    hasher.update(&[match cx.mode() {
        fs_exec::ExecMode::Deterministic => 0,
        fs_exec::ExecMode::Fast => 1,
    }]);
    let stream = cx.stream_key();
    for value in [stream.seed, stream.kernel_id, stream.tile, stream.iteration] {
        hasher.update(&value.to_le_bytes());
    }
    let budget = cx.budget();
    match budget.deadline {
        Some(deadline) => {
            hasher.update(&[1]);
            hasher.update(&deadline.as_nanos().to_le_bytes());
        }
        None => hasher.update(&[0]),
    }
    hasher.update(&budget.poll_quota.to_le_bytes());
    match budget.cost_quota {
        Some(cost_quota) => {
            hasher.update(&[1]);
            hasher.update(&cost_quota.to_le_bytes());
        }
        None => hasher.update(&[0]),
    }
    hasher.update(&[budget.priority]);
    hasher.update(&DWR_WORK_PLAN_VERSION.to_le_bytes());
    for field in work_plan_fields {
        hasher.update(&field.to_le_bytes());
    }
    hasher.update(&DWR_POLL_POLICY_VERSION.to_le_bytes());
    hasher.update(&DWR_POLL_STRIDE_IDENTITY.to_le_bytes());
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct VerifierPolicyIdentity {
    work_plan_version: u32,
    poll_policy_version: u32,
    poll_stride_work_units: u128,
}

const CURRENT_VERIFIER_POLICY_IDENTITY: VerifierPolicyIdentity = VerifierPolicyIdentity {
    work_plan_version: VERIFIER_WORK_PLAN_VERSION,
    poll_policy_version: VERIFIER_POLL_POLICY_VERSION,
    poll_stride_work_units: VERIFIER_POLL_STRIDE_WORK_UNITS,
};

fn hash_bracket_execution_header(
    hasher: &mut fs_blake3::Blake3,
    work_plan_fields: &[u128; 23],
    cx: &Cx<'_>,
    verifier_policy: VerifierPolicyIdentity,
) {
    hash_execution_header(hasher, DWR_BRACKET_IDENTITY_SCHEMA, work_plan_fields, cx);
    hasher.update(&verifier_policy.work_plan_version.to_le_bytes());
    hasher.update(&verifier_policy.poll_policy_version.to_le_bytes());
    hasher.update(&verifier_policy.poll_stride_work_units.to_le_bytes());
}

/// A QoI query: what the caller actually asked.
#[derive(Debug, Clone, PartialEq)]
pub struct DwrQuery {
    /// The quantity of interest (provenance label).
    pub qoi: String,
    /// The tolerance the answer must meet.
    pub tolerance: f64,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct DwrWorkPlan {
    mesh_nodes: u128,
    candidate_values: u128,
    polynomial_coefficients: u128,
    forcing_coefficients: u128,
    problem_identity_bytes: u128,
    coarse_cells: u128,
    refined_nodes: u128,
    refined_cells: u128,
    free_dual_nodes: u128,
    validation_work_units: u128,
    primal_work_units: u128,
    refinement_work_units: u128,
    initialization_work_units: u128,
    dual_assembly_work_units: u128,
    thomas_work_units: u128,
    dual_publication_work_units: u128,
    residual_work_units: u128,
    output_validation_work_units: u128,
    identity_work_units: u128,
    finalization_work_units: u128,
    planned_work_units: u128,
}

impl DwrWorkPlan {
    fn preflight(
        mesh_nodes: usize,
        candidate_values: usize,
        polynomial_coefficients: usize,
        forcing_coefficients: usize,
        problem_identity_bytes: usize,
    ) -> Result<Self, DwrError> {
        let mesh_nodes_input = mesh_nodes;
        let polynomial_coefficients_input = polynomial_coefficients;
        if !(2..=MAX_DWR_MESH_NODES).contains(&mesh_nodes) {
            return Err(DwrError::MeshNodeCount {
                count: mesh_nodes,
                minimum: 2,
                maximum: MAX_DWR_MESH_NODES,
            });
        }
        if candidate_values > MAX_DWR_MESH_NODES {
            return Err(DwrError::CandidateNodeCount {
                count: candidate_values,
                maximum: MAX_DWR_MESH_NODES,
            });
        }
        if candidate_values != mesh_nodes {
            return Err(DwrError::CandidateLengthMismatch {
                mesh_nodes,
                candidate_values,
            });
        }
        if !(1..=MAX_DWR_POLY_COEFFICIENTS).contains(&polynomial_coefficients) {
            return Err(DwrError::PolynomialCoefficientCount {
                count: polynomial_coefficients,
                minimum: 1,
                maximum: MAX_DWR_POLY_COEFFICIENTS,
            });
        }
        if !(1..=MAX_DWR_POLY_COEFFICIENTS).contains(&forcing_coefficients) {
            return Err(DwrError::PolynomialCoefficientCount {
                count: forcing_coefficients,
                minimum: 1,
                maximum: MAX_DWR_POLY_COEFFICIENTS,
            });
        }
        if problem_identity_bytes > MAX_FEM1D_PROBLEM_CANONICAL_IDENTITY_BYTES {
            return Err(DwrError::WorkPlanOverflow);
        }

        let refined_nodes_usize = refined_node_count(mesh_nodes)?;
        let mesh_nodes = u128::try_from(mesh_nodes).map_err(|_| DwrError::WorkPlanOverflow)?;
        let candidate_values =
            u128::try_from(candidate_values).map_err(|_| DwrError::WorkPlanOverflow)?;
        let polynomial_coefficients =
            u128::try_from(polynomial_coefficients).map_err(|_| DwrError::WorkPlanOverflow)?;
        let forcing_coefficients =
            u128::try_from(forcing_coefficients).map_err(|_| DwrError::WorkPlanOverflow)?;
        let problem_identity_bytes =
            u128::try_from(problem_identity_bytes).map_err(|_| DwrError::WorkPlanOverflow)?;
        let refined_nodes =
            u128::try_from(refined_nodes_usize).map_err(|_| DwrError::WorkPlanOverflow)?;
        let coarse_cells = mesh_nodes
            .checked_sub(1)
            .ok_or(DwrError::WorkPlanOverflow)?;
        let refined_cells = refined_nodes
            .checked_sub(1)
            .ok_or(DwrError::WorkPlanOverflow)?;
        let free_dual_nodes = refined_nodes
            .checked_sub(2)
            .ok_or(DwrError::WorkPlanOverflow)?;
        let validation_work_units = polynomial_coefficients
            .checked_add(candidate_values)
            .and_then(|total| total.checked_add(mesh_nodes))
            .and_then(|total| total.checked_add(coarse_cells))
            .and_then(|total| total.checked_add(forcing_coefficients))
            .ok_or(DwrError::WorkPlanOverflow)?;
        let primal_work_units = coarse_cells;
        let refinement_work_units = coarse_cells;
        let initialization_work_units = free_dual_nodes
            .checked_mul(5)
            .and_then(|units| units.checked_add(refined_nodes))
            .and_then(|units| units.checked_add(coarse_cells))
            .ok_or(DwrError::WorkPlanOverflow)?;
        let dual_assembly_work_units = refined_cells;
        let thomas_work_units = if free_dual_nodes == 0 {
            0
        } else {
            free_dual_nodes
                .checked_mul(6)
                .and_then(|units| units.checked_sub(1))
                .ok_or(DwrError::WorkPlanOverflow)?
        };
        let dual_publication_work_units = free_dual_nodes;
        let residual_work_units = coarse_cells;
        let output_validation_work_units = coarse_cells;
        let identity_work_units = problem_identity_bytes
            .checked_add(coarse_cells)
            .ok_or(DwrError::WorkPlanOverflow)?;
        let finalization_work_units = 1;
        let planned_work_units = [
            validation_work_units,
            primal_work_units,
            refinement_work_units,
            initialization_work_units,
            dual_assembly_work_units,
            thomas_work_units,
            dual_publication_work_units,
            residual_work_units,
            output_validation_work_units,
            identity_work_units,
            finalization_work_units,
        ]
        .into_iter()
        .try_fold(0_u128, |total, units| total.checked_add(units))
        .ok_or(DwrError::WorkPlanOverflow)?;
        if planned_work_units
            > u128::try_from(MAX_DWR_WORK_UNITS).map_err(|_| DwrError::WorkPlanOverflow)?
        {
            return Err(DwrError::WorkBudgetExceeded {
                mesh_nodes: mesh_nodes_input,
                polynomial_coefficients: polynomial_coefficients_input,
                estimated_work: usize::try_from(planned_work_units).ok(),
                maximum: MAX_DWR_WORK_UNITS,
            });
        }
        Ok(Self {
            mesh_nodes,
            candidate_values,
            polynomial_coefficients,
            forcing_coefficients,
            problem_identity_bytes,
            coarse_cells,
            refined_nodes,
            refined_cells,
            free_dual_nodes,
            validation_work_units,
            primal_work_units,
            refinement_work_units,
            initialization_work_units,
            dual_assembly_work_units,
            thomas_work_units,
            dual_publication_work_units,
            residual_work_units,
            output_validation_work_units,
            identity_work_units,
            finalization_work_units,
            planned_work_units,
        })
    }

    fn identity_fields(self) -> [u128; 21] {
        [
            self.mesh_nodes,
            self.candidate_values,
            self.polynomial_coefficients,
            self.forcing_coefficients,
            self.problem_identity_bytes,
            self.coarse_cells,
            self.refined_nodes,
            self.refined_cells,
            self.free_dual_nodes,
            self.validation_work_units,
            self.primal_work_units,
            self.refinement_work_units,
            self.initialization_work_units,
            self.dual_assembly_work_units,
            self.thomas_work_units,
            self.dual_publication_work_units,
            self.residual_work_units,
            self.output_validation_work_units,
            self.identity_work_units,
            self.finalization_work_units,
            self.planned_work_units,
        ]
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct BracketWorkPlan {
    primal_mesh_nodes: u128,
    primal_polynomial_coefficients: u128,
    primal_problem_identity_bytes: u128,
    primal_validation_work_units: u128,
    primal_verifier_plan: VerifierWorkPlan,
    dual_mesh_nodes: u128,
    dual_polynomial_coefficients: u128,
    dual_problem_identity_bytes: u128,
    dual_validation_work_units: u128,
    dual_verifier_plan: VerifierWorkPlan,
    identity_work_units: u128,
    finalization_work_units: u128,
    planned_work_units: u128,
}

impl BracketWorkPlan {
    #[allow(clippy::too_many_lines)] // One auditable checked sum for both nested verifier factors.
    fn preflight(
        primal_problem: &MmsProblem,
        primal_candidate: &[f64],
        dual_problem: &MmsProblem,
        dual_candidate: &[f64],
    ) -> Result<Self, BracketError> {
        fn factor_shape(
            factor: &'static str,
            problem: &MmsProblem,
            candidate: &[f64],
        ) -> Result<(u128, u128, u128, u128, VerifierWorkPlan), BracketError> {
            let verifier_plan =
                VerifierWorkPlan::for_inputs(problem, candidate).map_err(|reason| {
                    BracketError::InvalidInput {
                        factor,
                        reason: reason.id(),
                    }
                })?;
            let mesh_nodes = problem.mesh().len();
            let polynomial_coefficients = problem.exact_solution().coefficients().len();
            let mesh_nodes =
                u128::try_from(mesh_nodes).map_err(|_| BracketError::InvalidInput {
                    factor,
                    reason: "work plan overflowed",
                })?;
            let polynomial_coefficients =
                u128::try_from(polynomial_coefficients).map_err(|_| {
                    BracketError::InvalidInput {
                        factor,
                        reason: "work plan overflowed",
                    }
                })?;
            let problem_identity_bytes =
                u128::try_from(problem.canonical_bytes().len()).map_err(|_| {
                    BracketError::InvalidInput {
                        factor,
                        reason: "work plan overflowed",
                    }
                })?;
            let cells = mesh_nodes
                .checked_sub(1)
                .ok_or(BracketError::InvalidInput {
                    factor,
                    reason: "work plan overflowed",
                })?;
            let validation = mesh_nodes
                .checked_mul(2)
                .and_then(|units| units.checked_add(cells))
                .ok_or(BracketError::InvalidInput {
                    factor,
                    reason: "work plan overflowed",
                })?;
            Ok((
                mesh_nodes,
                polynomial_coefficients,
                problem_identity_bytes,
                validation,
                verifier_plan,
            ))
        }

        let (
            primal_mesh_nodes,
            primal_polynomial_coefficients,
            primal_problem_identity_bytes,
            primal_validation_work_units,
            primal_verifier_plan,
        ) = factor_shape("primal", primal_problem, primal_candidate)?;
        let (
            dual_mesh_nodes,
            dual_polynomial_coefficients,
            dual_problem_identity_bytes,
            dual_validation_work_units,
            dual_verifier_plan,
        ) = factor_shape("dual", dual_problem, dual_candidate)?;
        let identity_work_units = primal_problem_identity_bytes
            .checked_add(dual_problem_identity_bytes)
            .ok_or(BracketError::WorkPlanOverflow)?;
        let finalization_work_units = 1;
        let planned_work_units = primal_validation_work_units
            .checked_add(primal_verifier_plan.planned_work_units())
            .and_then(|work| work.checked_add(dual_validation_work_units))
            .and_then(|work| work.checked_add(dual_verifier_plan.planned_work_units()))
            .and_then(|work| work.checked_add(identity_work_units))
            .and_then(|work| work.checked_add(finalization_work_units))
            .ok_or(BracketError::WorkBudgetExceeded {
                planned_work_units: u128::MAX,
                maximum: MAX_DWR_WORK_UNITS,
            })?;
        let maximum_work_units =
            u128::try_from(MAX_DWR_WORK_UNITS).map_err(|_| BracketError::WorkBudgetExceeded {
                planned_work_units,
                maximum: MAX_DWR_WORK_UNITS,
            })?;
        if planned_work_units > maximum_work_units {
            return Err(BracketError::WorkBudgetExceeded {
                planned_work_units,
                maximum: MAX_DWR_WORK_UNITS,
            });
        }
        Ok(Self {
            primal_mesh_nodes,
            primal_polynomial_coefficients,
            primal_problem_identity_bytes,
            primal_validation_work_units,
            primal_verifier_plan,
            dual_mesh_nodes,
            dual_polynomial_coefficients,
            dual_problem_identity_bytes,
            dual_validation_work_units,
            dual_verifier_plan,
            identity_work_units,
            finalization_work_units,
            planned_work_units,
        })
    }

    fn identity_fields(self) -> [u128; 23] {
        let primal_verifier = self.primal_verifier_plan.identity_fields();
        let dual_verifier = self.dual_verifier_plan.identity_fields();
        [
            self.primal_mesh_nodes,
            self.primal_polynomial_coefficients,
            self.primal_problem_identity_bytes,
            self.primal_validation_work_units,
            primal_verifier[0],
            primal_verifier[1],
            primal_verifier[2],
            primal_verifier[3],
            primal_verifier[4],
            primal_verifier[5],
            self.dual_mesh_nodes,
            self.dual_polynomial_coefficients,
            self.dual_problem_identity_bytes,
            self.dual_validation_work_units,
            dual_verifier[0],
            dual_verifier[1],
            dual_verifier[2],
            dual_verifier[3],
            dual_verifier[4],
            dual_verifier[5],
            self.identity_work_units,
            self.finalization_work_units,
            self.planned_work_units,
        ]
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct AcceptWorkPlan {
    qoi_bytes: u128,
    decision_work_units: u128,
    finalization_work_units: u128,
    planned_work_units: u128,
}

impl AcceptWorkPlan {
    fn preflight(query: &DwrQuery) -> Result<Self, DwrError> {
        if query.qoi.len() > MAX_DWR_QOI_BYTES {
            return Err(DwrError::QoiLabelTooLong {
                bytes: query.qoi.len(),
                maximum: MAX_DWR_QOI_BYTES,
            });
        }
        let qoi_bytes = u128::try_from(query.qoi.len()).map_err(|_| DwrError::WorkPlanOverflow)?;
        let decision_work_units = 1;
        let finalization_work_units = 1;
        let planned_work_units = qoi_bytes
            .checked_add(decision_work_units)
            .and_then(|work| work.checked_add(finalization_work_units))
            .ok_or(DwrError::WorkPlanOverflow)?;
        Ok(Self {
            qoi_bytes,
            decision_work_units,
            finalization_work_units,
            planned_work_units,
        })
    }

    fn identity_fields(self) -> [u128; 4] {
        [
            self.qoi_bytes,
            self.decision_work_units,
            self.finalization_work_units,
            self.planned_work_units,
        ]
    }
}

/// An independently reverified primal/dual energy-product diagnostic.
///
/// The fields are sealed: safe downstream code cannot assert that an arbitrary
/// number is certified. A bracket can only be created by
/// [`Bracket::cauchy_schwarz`], which reruns the equilibrated-flux verifier on
/// the exact problem/candidate pairs. It is not yet a QoI-error certificate:
/// the v0 query type does not bind the dual problem to the requested functional.
///
/// ```compile_fail
/// use fs_adjoint::dwr_accept::Bracket;
///
/// let forged = Bracket {
///     bound: 0.0,
///     source: "caller assertion".to_string(),
/// };
/// ```
#[derive(Debug, Clone, PartialEq)]
pub struct Bracket {
    /// The outward-rounded product of the two energy-error upper bounds.
    bound: f64,
    /// Where the bound came from (audit trail).
    source: String,
    /// Retained identity of the verifier inputs, outputs, and execution policy.
    evidence_identity: DwrEvidenceIdentity,
}

/// Why a rigorous bracket could not be issued.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum BracketError {
    /// Cancellation was observed at a deterministic work boundary. No partial
    /// bracket is published.
    Cancelled {
        /// Stable phase at the observing checkpoint.
        phase: &'static str,
        /// Completed bounded logical work units.
        completed_work_units: u128,
        /// Complete preflighted logical work units.
        planned_work_units: u128,
    },
    /// A problem/candidate pair is malformed or exceeds the verifier envelope.
    InvalidInput {
        /// `primal` or `dual`.
        factor: &'static str,
        /// Stable refusal reason.
        reason: &'static str,
    },
    /// The independently rerun verifier panicked; no authority escaped.
    VerifierPanicked {
        /// `primal` or `dual`.
        factor: &'static str,
    },
    /// The verifier did not return a complete finite equilibrated certificate.
    VerifierRefused {
        /// `primal` or `dual`.
        factor: &'static str,
    },
    /// The outward-rounded product is not a finite usable QoI bound.
    ProductOverflow,
    /// The complete bracket work plan exceeds the public DWR work cap.
    WorkBudgetExceeded {
        /// Complete planned work units.
        planned_work_units: u128,
        /// Admitted maximum.
        maximum: usize,
    },
    /// Checked bracket progress or work-shape arithmetic overflowed.
    WorkPlanOverflow,
    /// Completed bracket work did not exactly match its preflighted plan.
    WorkPlanMismatch {
        /// Stable phase at which the mismatch was detected.
        phase: &'static str,
        /// Completed bounded logical work units.
        completed_work_units: u128,
        /// Expected bounded logical work units.
        planned_work_units: u128,
    },
}

impl core::fmt::Display for BracketError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Cancelled {
                phase,
                completed_work_units,
                planned_work_units,
            } => write!(
                f,
                "DWR bracket cancelled during {phase} after {completed_work_units}/{planned_work_units} work units"
            ),
            Self::InvalidInput { factor, reason } => {
                write!(f, "{factor} bracket input refused: {reason}")
            }
            Self::VerifierPanicked { factor } => {
                write!(f, "{factor} equilibrated verifier panicked")
            }
            Self::VerifierRefused { factor } => {
                write!(
                    f,
                    "{factor} equilibrated verifier produced no finite certificate"
                )
            }
            Self::ProductOverflow => f.write_str("Cauchy-Schwarz bracket product is not finite"),
            Self::WorkBudgetExceeded {
                planned_work_units,
                maximum,
            } => write!(
                f,
                "DWR bracket work plan {planned_work_units} exceeds maximum {maximum}"
            ),
            Self::WorkPlanOverflow => f.write_str("DWR bracket work plan overflowed"),
            Self::WorkPlanMismatch {
                phase,
                completed_work_units,
                planned_work_units,
            } => write!(
                f,
                "DWR bracket work plan mismatch during {phase}: completed {completed_work_units}/{planned_work_units} work units"
            ),
        }
    }
}

impl std::error::Error for BracketError {}

struct VerifiedFactor {
    report: fs_verify::estimator::VerifierReport,
    candidate_digest: fs_blake3::ContentHash,
}

fn bracket_verifier_checkpoint_phase(
    factor: &'static str,
    fallback: &'static str,
    phase: VerifierPhase,
) -> &'static str {
    match (factor, phase) {
        ("primal", VerifierPhase::Validation) => "dwr-bracket.primal-verifier.validation",
        ("primal", VerifierPhase::Tightness) => "dwr-bracket.primal-verifier.tightness",
        ("primal", VerifierPhase::Equilibrated) => "dwr-bracket.primal-verifier.equilibrated",
        ("primal", VerifierPhase::Hash) => "dwr-bracket.primal-verifier.hash",
        ("primal", VerifierPhase::Finalization) => "dwr-bracket.primal-verifier.finalization",
        ("dual", VerifierPhase::Validation) => "dwr-bracket.dual-verifier.validation",
        ("dual", VerifierPhase::Tightness) => "dwr-bracket.dual-verifier.tightness",
        ("dual", VerifierPhase::Equilibrated) => "dwr-bracket.dual-verifier.equilibrated",
        ("dual", VerifierPhase::Hash) => "dwr-bracket.dual-verifier.hash",
        ("dual", VerifierPhase::Finalization) => "dwr-bracket.dual-verifier.finalization",
        _ => fallback,
    }
}

#[derive(Debug, Default)]
struct VerifierObservationState {
    last_completed_work_units: u128,
    next_phase_index: usize,
    active_phase: Option<VerifierPhase>,
    saw_refusal_flush: bool,
    saw_publication: bool,
}

fn verifier_phase_position(plan: VerifierWorkPlan, phase: VerifierPhase) -> Option<(usize, u128)> {
    let after_validation = plan.validation_work_units();
    let after_tightness = after_validation.checked_add(plan.tightness_work_units())?;
    let after_equilibrated = after_tightness.checked_add(plan.equilibrated_work_units())?;
    let after_hash = after_equilibrated.checked_add(plan.hash_work_units())?;
    match phase {
        VerifierPhase::Validation => Some((0, 0)),
        VerifierPhase::Tightness => Some((1, after_validation)),
        VerifierPhase::Equilibrated => Some((2, after_tightness)),
        VerifierPhase::Hash => Some((3, after_equilibrated)),
        VerifierPhase::Finalization => Some((4, after_hash)),
        _ => None,
    }
}

fn verifier_trace_mismatch(
    phase: &'static str,
    snapshot: VerifierProgress,
    plan: VerifierWorkPlan,
) -> BracketError {
    BracketError::WorkPlanMismatch {
        phase,
        completed_work_units: snapshot.completed_work_units,
        planned_work_units: plan.planned_work_units(),
    }
}

fn verifier_observation_complete(
    state: &VerifierObservationState,
    plan: VerifierWorkPlan,
    refused: bool,
) -> bool {
    if refused {
        state.saw_refusal_flush && !state.saw_publication
    } else {
        state.last_completed_work_units == plan.planned_work_units()
            && state.saw_publication
            && !state.saw_refusal_flush
    }
}

#[allow(clippy::too_many_arguments)]
fn observe_verifier_progress(
    factor: &'static str,
    fallback_phase: &'static str,
    verifier_plan: VerifierWorkPlan,
    snapshot: VerifierProgress,
    state: &mut VerifierObservationState,
    progress: &mut WorkProgress,
    cx: &Cx<'_>,
) -> Result<(), BracketError> {
    let phase = bracket_verifier_checkpoint_phase(factor, fallback_phase, snapshot.phase);
    let Some((phase_index, phase_prefix)) = verifier_phase_position(verifier_plan, snapshot.phase)
    else {
        return Err(verifier_trace_mismatch(phase, snapshot, verifier_plan));
    };
    if snapshot.planned_work_units != verifier_plan.planned_work_units()
        || snapshot.completed_work_units < state.last_completed_work_units
        || snapshot.completed_work_units > verifier_plan.planned_work_units()
        || state.saw_refusal_flush
        || state.saw_publication
    {
        return Err(verifier_trace_mismatch(phase, snapshot, verifier_plan));
    }
    let delta = snapshot
        .completed_work_units
        .checked_sub(state.last_completed_work_units)
        .ok_or_else(|| verifier_trace_mismatch(phase, snapshot, verifier_plan))?;
    let next_work_boundary = state
        .last_completed_work_units
        .checked_div(VERIFIER_POLL_STRIDE_WORK_UNITS)
        .and_then(|completed_boundaries| completed_boundaries.checked_add(1))
        .and_then(|next_boundary| next_boundary.checked_mul(VERIFIER_POLL_STRIDE_WORK_UNITS))
        .ok_or_else(|| verifier_trace_mismatch(phase, snapshot, verifier_plan))?;
    let reaches_next_work_boundary = snapshot.completed_work_units == next_work_boundary;
    if delta > VERIFIER_POLL_STRIDE_WORK_UNITS
        || snapshot.completed_work_units > next_work_boundary
        || (reaches_next_work_boundary && snapshot.kind != VerifierCheckpointKind::WorkBoundary)
        || (snapshot.kind == VerifierCheckpointKind::WorkBoundary && !reaches_next_work_boundary)
    {
        return Err(verifier_trace_mismatch(phase, snapshot, verifier_plan));
    }
    match snapshot.kind {
        VerifierCheckpointKind::PhaseEntry => {
            if phase_index != state.next_phase_index
                || snapshot.completed_work_units != phase_prefix
            {
                return Err(verifier_trace_mismatch(phase, snapshot, verifier_plan));
            }
            state.next_phase_index += 1;
            state.active_phase = Some(snapshot.phase);
        }
        VerifierCheckpointKind::WorkBoundary => {
            if state.active_phase != Some(snapshot.phase)
                || delta == 0
                || snapshot.completed_work_units == 0
                || !snapshot
                    .completed_work_units
                    .is_multiple_of(VERIFIER_POLL_STRIDE_WORK_UNITS)
            {
                return Err(verifier_trace_mismatch(phase, snapshot, verifier_plan));
            }
        }
        VerifierCheckpointKind::RefusalFlush => {
            if state.active_phase != Some(snapshot.phase) {
                return Err(verifier_trace_mismatch(phase, snapshot, verifier_plan));
            }
            state.saw_refusal_flush = true;
        }
        VerifierCheckpointKind::Publication => {
            if snapshot.phase != VerifierPhase::Finalization
                || state.active_phase != Some(VerifierPhase::Finalization)
                || state.next_phase_index != 5
                || snapshot.completed_work_units != verifier_plan.planned_work_units()
            {
                return Err(verifier_trace_mismatch(phase, snapshot, verifier_plan));
            }
            state.saw_publication = true;
        }
        _ => return Err(verifier_trace_mismatch(phase, snapshot, verifier_plan)),
    }
    progress.advance_bracket(delta)?;
    state.last_completed_work_units = snapshot.completed_work_units;
    bracket_checkpoint(phase, progress, cx)
}

#[allow(clippy::too_many_arguments)]
#[allow(clippy::too_many_lines)] // One fail-closed nested verifier transaction and audit gate.
fn verify_factor(
    factor: &'static str,
    validation_phase: &'static str,
    verifier_phase: &'static str,
    validation_work_units: u128,
    verifier_plan: VerifierWorkPlan,
    problem: &MmsProblem,
    candidate: &[f64],
    progress: &mut WorkProgress,
    cx: &Cx<'_>,
) -> Result<VerifiedFactor, BracketError> {
    bracket_checkpoint(validation_phase, progress, cx)?;
    let validation_start = progress.completed_work_units;
    let mut candidate_hasher = fs_blake3::Blake3::new();
    candidate_hasher.update(b"fs-adjoint-dwr-bracket-candidate-v1");
    candidate_hasher.update(factor.as_bytes());
    for (index, value) in candidate.iter().enumerate() {
        poll_bracket_scan(index, validation_phase, progress, cx)?;
        if !value.is_finite() {
            return Err(BracketError::InvalidInput {
                factor,
                reason: "candidate contains a non-finite value",
            });
        }
        candidate_hasher.update(&value.to_bits().to_le_bytes());
        progress.advance_bracket(1)?;
    }
    if candidate.first().map(|value| value.to_bits()) != Some(0.0_f64.to_bits())
        || candidate.last().map(|value| value.to_bits()) != Some(0.0_f64.to_bits())
    {
        return Err(BracketError::InvalidInput {
            factor,
            reason: "candidate endpoints must be canonical homogeneous +0.0",
        });
    }
    bracket_checkpoint(validation_phase, progress, cx)?;
    for (index, value) in problem.mesh().iter().enumerate() {
        poll_bracket_scan(index, validation_phase, progress, cx)?;
        if !value.is_finite() {
            return Err(BracketError::InvalidInput {
                factor,
                reason: "mesh coordinates must be finite and strictly increasing",
            });
        }
        progress.advance_bracket(1)?;
    }
    bracket_checkpoint(validation_phase, progress, cx)?;
    for (index, pair) in problem.mesh().windows(2).enumerate() {
        poll_bracket_scan(index, validation_phase, progress, cx)?;
        if pair[0] >= pair[1] {
            return Err(BracketError::InvalidInput {
                factor,
                reason: "mesh coordinates must be finite and strictly increasing",
            });
        }
        progress.advance_bracket(1)?;
    }
    let completed_validation_work_units = progress
        .completed_work_units
        .checked_sub(validation_start)
        .ok_or(BracketError::WorkPlanOverflow)?;
    if completed_validation_work_units != validation_work_units {
        return Err(BracketError::WorkPlanMismatch {
            phase: validation_phase,
            completed_work_units: completed_validation_work_units,
            planned_work_units: validation_work_units,
        });
    }

    let mut verifier_observation = VerifierObservationState::default();
    let run = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        verify_with_checkpoint(problem, candidate, f64::MAX, |snapshot| {
            observe_verifier_progress(
                factor,
                verifier_phase,
                verifier_plan,
                snapshot,
                &mut verifier_observation,
                progress,
                cx,
            )
        })
    }));
    let report = match run {
        Err(_) => return Err(BracketError::VerifierPanicked { factor }),
        Ok(Err(error)) => return Err(error),
        Ok(Ok(report)) => report,
    };
    if !verifier_observation_complete(
        &verifier_observation,
        verifier_plan,
        report.refusal.is_some(),
    ) {
        return Err(BracketError::WorkPlanMismatch {
            phase: verifier_phase,
            completed_work_units: verifier_observation.last_completed_work_units,
            planned_work_units: verifier_plan.planned_work_units(),
        });
    }
    let color_matches = matches!(
        &report.color,
        Some(Color::Verified { lo, hi })
            if lo.to_bits() == 0.0_f64.to_bits()
                && hi.to_bits() == report.bound.hi.to_bits()
    );
    if report.refusal.is_some()
        || report.family != EstimatorFamily::EquilibratedFlux.id()
        || report.tolerance.to_bits() != f64::MAX.to_bits()
        || !report.accept
        || !color_matches
        || report.bound.lo.is_nan()
        || !report.bound.hi.is_finite()
        || report.bound.lo < 0.0
        || report.bound.hi < report.bound.lo
    {
        return Err(BracketError::VerifierRefused { factor });
    }
    Ok(VerifiedFactor {
        report,
        candidate_digest: candidate_hasher.finalize(),
    })
}

impl Bracket {
    /// The Cauchy–Schwarz bracket from two equilibrated energy-norm
    /// enclosures: `|a(e_u, e_z)| ≤ ‖e_u‖_E · ‖e_z‖_E`, outward-rounded.
    /// Both reports are independently recomputed here; public report fields are
    /// never accepted as authority.
    ///
    /// # Errors
    /// [`BracketError`] when either input is malformed, verification panics or
    /// refuses, or the product overflows.
    pub fn cauchy_schwarz(
        primal_problem: &MmsProblem,
        primal_candidate: &[f64],
        dual_problem: &MmsProblem,
        dual_candidate: &[f64],
        cx: &Cx<'_>,
    ) -> Result<Bracket, BracketError> {
        let plan = BracketWorkPlan::preflight(
            primal_problem,
            primal_candidate,
            dual_problem,
            dual_candidate,
        )?;
        let mut progress = WorkProgress::new(plan.planned_work_units, cx);
        bracket_checkpoint(BRACKET_INITIAL_PHASE, &mut progress, cx)?;
        let primal = verify_factor(
            "primal",
            BRACKET_PRIMAL_VALIDATE_PHASE,
            BRACKET_PRIMAL_VERIFY_PHASE,
            plan.primal_validation_work_units,
            plan.primal_verifier_plan,
            primal_problem,
            primal_candidate,
            &mut progress,
            cx,
        )?;
        let dual = verify_factor(
            "dual",
            BRACKET_DUAL_VALIDATE_PHASE,
            BRACKET_DUAL_VERIFY_PHASE,
            plan.dual_validation_work_units,
            plan.dual_verifier_plan,
            dual_problem,
            dual_candidate,
            &mut progress,
            cx,
        )?;
        let raw_bound = primal.report.bound.hi * dual.report.bound.hi;
        let bound = up(raw_bound);
        if !bound.is_finite() || bound < 0.0 {
            return Err(BracketError::ProductOverflow);
        }
        let source = format!(
            "cauchy-schwarz(equilibrated primal {:.3e} flux {:016x} x equilibrated dual {:.3e} flux {:016x})",
            primal.report.bound.hi,
            primal.report.flux_hash,
            dual.report.bound.hi,
            dual.report.flux_hash
        );
        let mut hasher = fs_blake3::Blake3::new();
        hash_bracket_execution_header(
            &mut hasher,
            &plan.identity_fields(),
            cx,
            CURRENT_VERIFIER_POLICY_IDENTITY,
        );
        for problem in [primal_problem, dual_problem] {
            hasher.update(&problem.identity().version().to_le_bytes());
            hasher.update(&problem.identity().root().to_le_bytes());
            let canonical = problem.canonical_bytes();
            let canonical_len =
                u64::try_from(canonical.len()).map_err(|_| BracketError::WorkPlanOverflow)?;
            hasher.update(&canonical_len.to_le_bytes());
            for chunk in canonical.chunks(DWR_POLL_STRIDE_ITEMS) {
                bracket_checkpoint(BRACKET_IDENTITY_PHASE, &mut progress, cx)?;
                hasher.update(chunk);
                progress.advance_bracket(
                    u128::try_from(chunk.len()).map_err(|_| BracketError::WorkPlanOverflow)?,
                )?;
            }
        }
        hasher.update(primal.candidate_digest.as_bytes());
        hasher.update(dual.candidate_digest.as_bytes());
        for report in [&primal.report, &dual.report] {
            hasher.update(&report.bound.lo.to_bits().to_le_bytes());
            hasher.update(&report.bound.hi.to_bits().to_le_bytes());
            hasher.update(&report.flux_hash.to_le_bytes());
        }
        hasher.update(&bound.to_bits().to_le_bytes());
        let source_bytes = source.as_bytes();
        let source_len =
            u64::try_from(source_bytes.len()).map_err(|_| BracketError::WorkPlanOverflow)?;
        hasher.update(&source_len.to_le_bytes());
        hasher.update(source_bytes);
        let bracket = Bracket {
            bound,
            source,
            evidence_identity: DwrEvidenceIdentity(hasher.finalize()),
        };
        progress.advance_bracket(plan.finalization_work_units)?;
        progress.finish_bracket(BRACKET_PUBLICATION_PHASE)?;
        bracket_checkpoint(BRACKET_PUBLICATION_PHASE, &mut progress, cx)?;
        Ok(bracket)
    }

    /// Outward-rounded energy-error product. This is diagnostic until a typed
    /// QoI-dual relation is verified.
    #[must_use]
    pub fn bound(&self) -> f64 {
        self.bound
    }

    /// Deterministic verifier audit label.
    #[must_use]
    pub fn source(&self) -> &str {
        &self.source
    }

    /// Retained execution/evidence identity for this exact bracket.
    #[must_use]
    pub const fn evidence_identity(&self) -> DwrEvidenceIdentity {
        self.evidence_identity
    }
}

/// The colored accept outcome.
#[derive(Debug, Clone, PartialEq)]
pub struct AcceptOutcome {
    /// Was the query discharged?
    accepted: bool,
    /// The color the answer carries.
    color: Color,
    /// True when malformed public inputs prevented an accept/reject decision.
    /// False includes both a valid accept and a valid over-tolerance rejection.
    refused: bool,
    /// The audit trail.
    audit: String,
    /// Retained identity of the decision, inputs, and execution policy.
    evidence_identity: DwrEvidenceIdentity,
}

struct AcceptDraft {
    accepted: bool,
    color: Color,
    refused: bool,
    audit: String,
}

fn malformed_refusal(estimator: &str, audit: String) -> AcceptDraft {
    AcceptDraft {
        accepted: false,
        color: Color::Estimated {
            estimator: estimator.to_string(),
            dispersion: f64::INFINITY,
        },
        refused: true,
        audit,
    }
}

/// The accept test. Color logic (mechanical, auditable):
/// - no acceptance path → rejected (estimated color on the estimate);
/// - DWR-only accept (`|η| ≤ tol`, no valid bracket) → ESTIMATED;
/// - a sealed energy-product bracket is retained in the audit, but cannot
///   promote or veto because the v0 query does not bind its dual relation.
pub fn accept(
    query: &DwrQuery,
    dwr_abs: f64,
    bracket: Option<&Bracket>,
    cx: &Cx<'_>,
) -> Result<AcceptOutcome, DwrError> {
    let plan = AcceptWorkPlan::preflight(query)?;
    let mut progress = WorkProgress::new(plan.planned_work_units, cx);
    dwr_checkpoint(ACCEPT_INITIAL_PHASE, &mut progress, cx)?;
    let draft = if !query.tolerance.is_finite() || query.tolerance <= 0.0 {
        malformed_refusal(
            "dwr-invalid-tolerance",
            format!(
                "REFUSED: tolerance must be finite and positive, got {:.3e}",
                query.tolerance
            ),
        )
    } else if let Some(b) = bracket
        && (!b.bound.is_finite() || b.bound < 0.0)
    {
        malformed_refusal(
            "dwr-invalid-guaranteed-bracket",
            format!(
                "REFUSED: guaranteed bracket from {} has invalid bound {:.3e}",
                b.source, b.bound
            ),
        )
    } else if !dwr_abs.is_finite() || dwr_abs < 0.0 {
        malformed_refusal(
            "dwr-invalid-estimate",
            format!(
                "REFUSED: DWR absolute estimate must be finite and non-negative, got {dwr_abs:.3e}"
            ),
        )
    } else if dwr_abs <= query.tolerance {
        let bracket_note = bracket.map_or_else(
            || "no energy-product diagnostic".to_string(),
            |value| {
                format!(
                    "energy-product diagnostic {:.3e} from {}; QoI-dual relation unverified",
                    value.bound, value.source
                )
            },
        );
        AcceptDraft {
            accepted: true,
            color: Color::Estimated {
                estimator: if bracket.is_some() {
                    "dwr-with-unbound-energy-diagnostic".to_string()
                } else {
                    "dwr-unbracketed".to_string()
                },
                dispersion: dwr_abs,
            },
            refused: false,
            audit: format!(
                "estimated-only accept: dwr {:.3e} <= tol {:.3e}; {bracket_note}",
                dwr_abs, query.tolerance,
            ),
        }
    } else {
        let bracket_note = bracket.map_or_else(
            || "no energy-product diagnostic".to_string(),
            |value| {
                format!(
                    "energy-product diagnostic {:.3e} from {}; QoI-dual relation unverified",
                    value.bound, value.source
                )
            },
        );
        AcceptDraft {
            accepted: false,
            color: Color::Estimated {
                estimator: "dwr-rejected".to_string(),
                dispersion: dwr_abs,
            },
            refused: false,
            audit: format!(
                "rejected: dwr {:.3e} > tol {:.3e}; {bracket_note}",
                dwr_abs, query.tolerance,
            ),
        }
    };
    progress.advance(plan.decision_work_units)?;

    dwr_checkpoint(ACCEPT_IDENTITY_PHASE, &mut progress, cx)?;
    let mut hasher = fs_blake3::Blake3::new();
    hash_execution_header(
        &mut hasher,
        DWR_ACCEPT_IDENTITY_SCHEMA,
        &plan.identity_fields(),
        cx,
    );
    for chunk in query.qoi.as_bytes().chunks(DWR_POLL_STRIDE_ITEMS) {
        dwr_checkpoint(ACCEPT_IDENTITY_PHASE, &mut progress, cx)?;
        hasher.update(chunk);
        progress.advance(u128::try_from(chunk.len()).map_err(|_| DwrError::WorkPlanOverflow)?)?;
    }
    hasher.update(&query.tolerance.to_bits().to_le_bytes());
    hasher.update(&dwr_abs.to_bits().to_le_bytes());
    match bracket {
        Some(value) => {
            hasher.update(&[1]);
            hasher.update(value.evidence_identity.as_bytes());
        }
        None => hasher.update(&[0]),
    }
    hasher.update(&[u8::from(draft.accepted), u8::from(draft.refused)]);
    let color_bytes = draft.color.canonical_bytes();
    let color_len = u64::try_from(color_bytes.len()).map_err(|_| DwrError::WorkPlanOverflow)?;
    hasher.update(&color_len.to_le_bytes());
    hasher.update(&color_bytes);
    let audit_bytes = draft.audit.as_bytes();
    let audit_len = u64::try_from(audit_bytes.len()).map_err(|_| DwrError::WorkPlanOverflow)?;
    hasher.update(&audit_len.to_le_bytes());
    hasher.update(audit_bytes);
    let outcome = AcceptOutcome {
        accepted: draft.accepted,
        color: draft.color,
        refused: draft.refused,
        audit: draft.audit,
        evidence_identity: DwrEvidenceIdentity(hasher.finalize()),
    };
    progress.advance(plan.finalization_work_units)?;
    progress.finish_dwr(ACCEPT_PUBLICATION_PHASE)?;
    dwr_checkpoint(ACCEPT_PUBLICATION_PHASE, &mut progress, cx)?;
    Ok(outcome)
}

impl AcceptOutcome {
    /// Whether the query was discharged.
    #[must_use]
    pub const fn accepted(&self) -> bool {
        self.accepted
    }

    /// Evidence color carried by this exact decision.
    #[must_use]
    pub const fn color(&self) -> &Color {
        &self.color
    }

    /// Whether malformed inputs prevented an accept/reject decision.
    #[must_use]
    pub const fn refused(&self) -> bool {
        self.refused
    }

    /// Deterministic audit trail bound into the retained identity.
    #[must_use]
    pub fn audit(&self) -> &str {
        &self.audit
    }

    /// Retained execution/evidence identity for this exact decision.
    #[must_use]
    pub const fn evidence_identity(&self) -> DwrEvidenceIdentity {
        self.evidence_identity
    }
}

/// The 1-D reference DWR estimator for integral QoIs
/// `J(u) = ∫_{w_lo}^{w_hi} u dx` over an fs-verify problem: the dual
/// `−z″ = 1_{[w_lo, w_hi]}` solves by P1 FEM on the ONCE-REFINED mesh
/// (the enriched dual), and the estimate is the dual-weighted residual
/// `η = r(z_f − I_h z_f)` with per-COARSE-element indicators.
#[derive(Debug, Clone, PartialEq)]
pub struct DwrOutput {
    /// `J(u_h)`.
    j_primal: f64,
    /// The signed estimate `η ≈ J(u) − J(u_h)`.
    eta: f64,
    /// Per-coarse-element |indicator| (refinement guidance).
    indicators: Vec<f64>,
    /// Retained identity of the estimator inputs, outputs, and execution policy.
    evidence_identity: DwrEvidenceIdentity,
}

impl DwrOutput {
    /// `J(u_h)` for the exact retained execution.
    #[must_use]
    pub const fn j_primal(&self) -> f64 {
        self.j_primal
    }

    /// Signed DWR estimate for the exact retained execution.
    #[must_use]
    pub const fn eta(&self) -> f64 {
        self.eta
    }

    /// Per-coarse-element absolute indicators.
    #[must_use]
    pub fn indicators(&self) -> &[f64] {
        &self.indicators
    }

    /// Retained execution/evidence identity for this exact DWR estimate.
    #[must_use]
    pub const fn evidence_identity(&self) -> DwrEvidenceIdentity {
        self.evidence_identity
    }
}

/// Why the public DWR execution path refused an input or derived state.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum DwrError {
    /// Cancellation was observed at a deterministic work boundary. No partial
    /// estimate or accept outcome is published.
    Cancelled {
        /// Stable phase at the observing checkpoint.
        phase: &'static str,
        /// Completed bounded logical work units.
        completed_work_units: u128,
        /// Complete preflighted logical work units.
        planned_work_units: u128,
    },
    /// The QoI provenance label exceeds the bounded acceptance envelope.
    QoiLabelTooLong {
        /// Supplied UTF-8 bytes.
        bytes: usize,
        /// Admitted maximum.
        maximum: usize,
    },
    /// Coarse mesh node count is outside the bounded execution envelope.
    MeshNodeCount {
        /// Supplied node count.
        count: usize,
        /// Required minimum.
        minimum: usize,
        /// Admitted maximum.
        maximum: usize,
    },
    /// Candidate node count exceeds the bounded execution envelope.
    CandidateNodeCount {
        /// Supplied value count.
        count: usize,
        /// Admitted maximum.
        maximum: usize,
    },
    /// Candidate and mesh shapes differ.
    CandidateLengthMismatch {
        /// Mesh node count.
        mesh_nodes: usize,
        /// Candidate value count.
        candidate_values: usize,
    },
    /// Manufactured-solution polynomial size is outside the bounded envelope.
    PolynomialCoefficientCount {
        /// Supplied coefficient count.
        count: usize,
        /// Required minimum.
        minimum: usize,
        /// Admitted maximum.
        maximum: usize,
    },
    /// A manufactured-solution coefficient is NaN or infinite.
    NonFinitePolynomialCoefficient {
        /// Coefficient index.
        index: usize,
    },
    /// A candidate nodal value is NaN or infinite.
    NonFiniteCandidate {
        /// Candidate index.
        index: usize,
    },
    /// Candidate endpoints are not canonical homogeneous `+0.0` values.
    CandidateBoundary,
    /// A mesh coordinate is NaN or infinite.
    NonFiniteMeshNode {
        /// Mesh-node index.
        index: usize,
    },
    /// A mesh cell is not strictly increasing.
    NonIncreasingMeshCell {
        /// Left node / cell index.
        cell: usize,
    },
    /// Subtracting a cell's finite endpoints produced a non-finite width.
    NonFiniteCellWidth {
        /// Left node / cell index.
        cell: usize,
    },
    /// The once-refined midpoint is not strictly inside its coarse cell.
    NonInteriorMidpoint {
        /// Coarse cell index.
        cell: usize,
    },
    /// A coarse or refined cell has a non-finite reciprocal width.
    NonFiniteReciprocal {
        /// Coarse cell index.
        cell: usize,
        /// `None` for the coarse cell, `Some(0|1)` for a refined half.
        refined_half: Option<u8>,
    },
    /// The QoI integration window is non-finite or inverted.
    InvalidQoiWindow {
        /// Stable refusal reason.
        reason: &'static str,
    },
    /// Computing `2 * nodes - 1` overflowed or exceeded the refined cap.
    RefinedMeshSizeOverflow {
        /// Supplied coarse node count.
        mesh_nodes: usize,
    },
    /// The mesh/polynomial cross-product exceeds the bounded execution budget.
    WorkBudgetExceeded {
        /// Supplied coarse node count.
        mesh_nodes: usize,
        /// Supplied manufactured-solution coefficient count.
        polynomial_coefficients: usize,
        /// Estimated work, or `None` when the estimate itself overflowed.
        estimated_work: Option<usize>,
        /// Admitted maximum.
        maximum: usize,
    },
    /// Checked work-shape or progress arithmetic overflowed.
    WorkPlanOverflow,
    /// Completed DWR work did not exactly match its preflighted plan.
    WorkPlanMismatch {
        /// Stable phase at which the mismatch was detected.
        phase: &'static str,
        /// Completed bounded logical work units.
        completed_work_units: u128,
        /// Expected bounded logical work units.
        planned_work_units: u128,
    },
    /// A bounded scientific vector could not reserve its complete capacity.
    AllocationFailed {
        /// Stable allocation phase.
        phase: &'static str,
        /// Requested element count.
        elements: usize,
    },
    /// Tridiagonal storage lengths are inconsistent.
    LinearSystemShape,
    /// Assembly or elimination produced a non-finite linear-system value.
    NonFiniteLinearSystem {
        /// Stable component/stage name.
        component: &'static str,
        /// Row index.
        index: usize,
    },
    /// Thomas elimination encountered a non-positive or non-finite pivot in a
    /// system whose exact stiffness matrix is positive definite.
    InvalidLinearPivot {
        /// Pivot row.
        row: usize,
    },
    /// A derived quadrature, slope, residual, or output value is non-finite.
    NonFiniteDerived {
        /// Stable derived quantity name.
        quantity: &'static str,
        /// Cell/coefficient index when applicable.
        index: Option<usize>,
    },
    /// A clipped P1 integral rounded to zero without an exact-zero proof and
    /// therefore cannot safely drive an acceptance decision in binary64.
    UnresolvedZeroIntegral {
        /// Stable integral role.
        quantity: &'static str,
        /// Mesh-cell index.
        cell: usize,
    },
}

impl core::fmt::Display for DwrError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Cancelled {
                phase,
                completed_work_units,
                planned_work_units,
            } => write!(
                f,
                "DWR workflow cancelled during {phase} after {completed_work_units}/{planned_work_units} work units"
            ),
            Self::QoiLabelTooLong { bytes, maximum } => {
                write!(f, "DWR QoI label has {bytes} bytes; maximum is {maximum}")
            }
            Self::MeshNodeCount {
                count,
                minimum,
                maximum,
            } => write!(
                f,
                "DWR mesh has {count} nodes; expected {minimum}..={maximum}"
            ),
            Self::CandidateNodeCount { count, maximum } => {
                write!(f, "DWR candidate has {count} values; maximum is {maximum}")
            }
            Self::CandidateLengthMismatch {
                mesh_nodes,
                candidate_values,
            } => write!(
                f,
                "DWR candidate length {candidate_values} differs from mesh length {mesh_nodes}"
            ),
            Self::PolynomialCoefficientCount {
                count,
                minimum,
                maximum,
            } => write!(
                f,
                "DWR polynomial has {count} coefficients; expected {minimum}..={maximum}"
            ),
            Self::NonFinitePolynomialCoefficient { index } => {
                write!(f, "DWR polynomial coefficient {index} is non-finite")
            }
            Self::NonFiniteCandidate { index } => {
                write!(f, "DWR candidate value {index} is non-finite")
            }
            Self::CandidateBoundary => {
                f.write_str("DWR candidate endpoints must be canonical homogeneous +0.0")
            }
            Self::NonFiniteMeshNode { index } => {
                write!(f, "DWR mesh node {index} is non-finite")
            }
            Self::NonIncreasingMeshCell { cell } => {
                write!(f, "DWR mesh cell {cell} is not strictly increasing")
            }
            Self::NonFiniteCellWidth { cell } => {
                write!(f, "DWR mesh cell {cell} has a non-finite width")
            }
            Self::NonInteriorMidpoint { cell } => write!(
                f,
                "DWR mesh cell {cell} has no representable strictly interior midpoint"
            ),
            Self::NonFiniteReciprocal { cell, refined_half } => match refined_half {
                Some(half) => write!(
                    f,
                    "DWR refined half {half} of coarse cell {cell} has a non-finite reciprocal width"
                ),
                None => write!(
                    f,
                    "DWR coarse cell {cell} has a non-finite reciprocal width"
                ),
            },
            Self::InvalidQoiWindow { reason } => {
                write!(f, "DWR QoI window refused: {reason}")
            }
            Self::RefinedMeshSizeOverflow { mesh_nodes } => write!(
                f,
                "DWR refined mesh size overflowed its bound for {mesh_nodes} coarse nodes"
            ),
            Self::WorkBudgetExceeded {
                mesh_nodes,
                polynomial_coefficients,
                estimated_work,
                maximum,
            } => match estimated_work {
                Some(work) => write!(
                    f,
                    "DWR work estimate {work} for {mesh_nodes} mesh nodes x {polynomial_coefficients} coefficients exceeds {maximum}"
                ),
                None => write!(
                    f,
                    "DWR work estimate overflowed for {mesh_nodes} mesh nodes x {polynomial_coefficients} coefficients (maximum {maximum})"
                ),
            },
            Self::WorkPlanOverflow => f.write_str("DWR complete work plan overflowed"),
            Self::WorkPlanMismatch {
                phase,
                completed_work_units,
                planned_work_units,
            } => write!(
                f,
                "DWR work plan mismatch during {phase}: completed {completed_work_units}/{planned_work_units} work units"
            ),
            Self::AllocationFailed { phase, elements } => write!(
                f,
                "DWR allocation during {phase} could not reserve {elements} elements"
            ),
            Self::LinearSystemShape => {
                f.write_str("DWR tridiagonal linear-system shapes are inconsistent")
            }
            Self::NonFiniteLinearSystem { component, index } => write!(
                f,
                "DWR linear-system {component} is non-finite at row {index}"
            ),
            Self::InvalidLinearPivot { row } => {
                write!(f, "DWR Thomas pivot {row} is non-positive or non-finite")
            }
            Self::NonFiniteDerived { quantity, index } => match index {
                Some(index) => write!(f, "DWR derived {quantity} is non-finite at index {index}"),
                None => write!(f, "DWR derived {quantity} is non-finite"),
            },
            Self::UnresolvedZeroIntegral { quantity, cell } => write!(
                f,
                "DWR clipped {quantity} on cell {cell} is zero in binary64 without an exact-zero proof"
            ),
        }
    }
}

impl std::error::Error for DwrError {}

fn refined_node_count(mesh_nodes: usize) -> Result<usize, DwrError> {
    mesh_nodes
        .checked_mul(2)
        .and_then(|count| count.checked_sub(1))
        .filter(|&count| count <= MAX_DWR_REFINED_NODES)
        .ok_or(DwrError::RefinedMeshSizeOverflow { mesh_nodes })
}

struct ValidatedDwrInputs {
    refined_nodes: usize,
    candidate_digest: fs_blake3::ContentHash,
}

fn validate_dwr_inputs(
    problem: &MmsProblem,
    candidate: &[f64],
    w_lo: f64,
    w_hi: f64,
    progress: &mut WorkProgress,
    cx: &Cx<'_>,
) -> Result<ValidatedDwrInputs, DwrError> {
    let mesh = problem.mesh();
    let coefficients = problem.exact_solution().coefficients();
    if !w_lo.is_finite() || !w_hi.is_finite() {
        return Err(DwrError::InvalidQoiWindow {
            reason: "endpoints must be finite",
        });
    }
    if w_lo >= w_hi {
        return Err(DwrError::InvalidQoiWindow {
            reason: "lower endpoint must be strictly below upper endpoint",
        });
    }

    dwr_checkpoint(DWR_VALIDATE_POLYNOMIAL_PHASE, progress, cx)?;
    for (index, value) in coefficients.iter().enumerate() {
        poll_dwr_scan(index, DWR_VALIDATE_POLYNOMIAL_PHASE, progress, cx)?;
        if !value.is_finite() {
            return Err(DwrError::NonFinitePolynomialCoefficient { index });
        }
        progress.advance(1)?;
    }

    dwr_checkpoint(DWR_VALIDATE_CANDIDATE_PHASE, progress, cx)?;
    let mut candidate_hasher = fs_blake3::Blake3::new();
    candidate_hasher.update(b"fs-adjoint-dwr-candidate-v1");
    for (index, value) in candidate.iter().enumerate() {
        poll_dwr_scan(index, DWR_VALIDATE_CANDIDATE_PHASE, progress, cx)?;
        if !value.is_finite() {
            return Err(DwrError::NonFiniteCandidate { index });
        }
        candidate_hasher.update(&value.to_bits().to_le_bytes());
        progress.advance(1)?;
    }
    if candidate[0].to_bits() != 0.0_f64.to_bits()
        || candidate[candidate.len() - 1].to_bits() != 0.0_f64.to_bits()
    {
        return Err(DwrError::CandidateBoundary);
    }

    dwr_checkpoint(DWR_VALIDATE_MESH_PHASE, progress, cx)?;
    for (index, value) in mesh.iter().enumerate() {
        poll_dwr_scan(index, DWR_VALIDATE_MESH_PHASE, progress, cx)?;
        if !value.is_finite() {
            return Err(DwrError::NonFiniteMeshNode { index });
        }
        progress.advance(1)?;
    }

    dwr_checkpoint(DWR_VALIDATE_CELLS_PHASE, progress, cx)?;
    for (cell, nodes) in mesh.windows(2).enumerate() {
        poll_dwr_scan(cell, DWR_VALIDATE_CELLS_PHASE, progress, cx)?;
        let (x0, x1) = (nodes[0], nodes[1]);
        if x0 >= x1 {
            return Err(DwrError::NonIncreasingMeshCell { cell });
        }
        let width = x1 - x0;
        if !width.is_finite() {
            return Err(DwrError::NonFiniteCellWidth { cell });
        }
        if !(1.0 / width).is_finite() {
            return Err(DwrError::NonFiniteReciprocal {
                cell,
                refined_half: None,
            });
        }
        let midpoint = f64::midpoint(x0, x1);
        if !(x0 < midpoint && midpoint < x1) {
            return Err(DwrError::NonInteriorMidpoint { cell });
        }
        for (half, half_width) in [(0_u8, midpoint - x0), (1_u8, x1 - midpoint)] {
            if !half_width.is_finite() || !(1.0 / half_width).is_finite() {
                return Err(DwrError::NonFiniteReciprocal {
                    cell,
                    refined_half: Some(half),
                });
            }
        }
        progress.advance(1)?;
    }
    Ok(ValidatedDwrInputs {
        refined_nodes: refined_node_count(mesh.len())?,
        candidate_digest: candidate_hasher.finalize(),
    })
}

fn non_finite_derived(quantity: &'static str, index: Option<usize>) -> DwrError {
    DwrError::NonFiniteDerived { quantity, index }
}

fn zeroed_vec(
    elements: usize,
    phase: &'static str,
    progress: &mut WorkProgress,
    cx: &Cx<'_>,
) -> Result<Vec<f64>, DwrError> {
    dwr_checkpoint(phase, progress, cx)?;
    let mut values = Vec::new();
    values
        .try_reserve_exact(elements)
        .map_err(|_| DwrError::AllocationFailed { phase, elements })?;
    dwr_checkpoint(phase, progress, cx)?;
    for index in 0..elements {
        poll_dwr_scan(index, phase, progress, cx)?;
        values.push(0.0);
        progress.advance(1)?;
    }
    Ok(values)
}

fn thomas_solve(
    sub: &[f64],
    diag: &[f64],
    sup: &[f64],
    rhs: &mut [f64],
    progress: &mut WorkProgress,
    cx: &Cx<'_>,
) -> Result<(), DwrError> {
    let n = rhs.len();
    if sub.len() != n || diag.len() != n || sup.len() != n {
        return Err(DwrError::LinearSystemShape);
    }
    dwr_checkpoint(DWR_THOMAS_VALIDATE_PHASE, progress, cx)?;
    if n == 0 {
        return Ok(());
    }
    for (component, values) in [
        ("subdiagonal", sub),
        ("diagonal", diag),
        ("superdiagonal", sup),
    ] {
        dwr_checkpoint(DWR_THOMAS_VALIDATE_PHASE, progress, cx)?;
        for (index, value) in values.iter().enumerate() {
            poll_dwr_scan(index, DWR_THOMAS_VALIDATE_PHASE, progress, cx)?;
            if !value.is_finite() {
                return Err(DwrError::NonFiniteLinearSystem { component, index });
            }
            progress.advance(1)?;
        }
    }
    dwr_checkpoint(DWR_THOMAS_VALIDATE_PHASE, progress, cx)?;
    for (index, value) in rhs.iter().enumerate() {
        poll_dwr_scan(index, DWR_THOMAS_VALIDATE_PHASE, progress, cx)?;
        if !value.is_finite() {
            return Err(DwrError::NonFiniteLinearSystem {
                component: "right-hand side",
                index,
            });
        }
        progress.advance(1)?;
    }
    let mut c = zeroed_vec(n, DWR_THOMAS_FORWARD_PHASE, progress, cx)?;
    dwr_checkpoint(DWR_THOMAS_FORWARD_PHASE, progress, cx)?;
    let mut d = diag[0];
    if !d.is_finite() || d <= 0.0 {
        return Err(DwrError::InvalidLinearPivot { row: 0 });
    }
    if n > 1 {
        c[0] = sup[0] / d;
        if !c[0].is_finite() {
            return Err(DwrError::NonFiniteLinearSystem {
                component: "forward coefficient",
                index: 0,
            });
        }
    }
    rhs[0] /= d;
    if !rhs[0].is_finite() {
        return Err(DwrError::NonFiniteLinearSystem {
            component: "forward right-hand side",
            index: 0,
        });
    }
    progress.advance(1)?;
    for i in 1..n {
        poll_dwr_scan(i, DWR_THOMAS_FORWARD_PHASE, progress, cx)?;
        d = diag[i] - sub[i] * c[i - 1];
        if !d.is_finite() || d <= 0.0 {
            return Err(DwrError::InvalidLinearPivot { row: i });
        }
        if i < n - 1 {
            c[i] = sup[i] / d;
            if !c[i].is_finite() {
                return Err(DwrError::NonFiniteLinearSystem {
                    component: "forward coefficient",
                    index: i,
                });
            }
        }
        rhs[i] = (rhs[i] - sub[i] * rhs[i - 1]) / d;
        if !rhs[i].is_finite() {
            return Err(DwrError::NonFiniteLinearSystem {
                component: "forward right-hand side",
                index: i,
            });
        }
        progress.advance(1)?;
    }
    dwr_checkpoint(DWR_THOMAS_BACK_PHASE, progress, cx)?;
    for i in (0..n - 1).rev() {
        let completed = n - 2 - i;
        if completed != 0 && completed.is_multiple_of(DWR_POLL_STRIDE_ITEMS) {
            dwr_checkpoint(DWR_THOMAS_BACK_PHASE, progress, cx)?;
        }
        rhs[i] -= c[i] * rhs[i + 1];
        if !rhs[i].is_finite() {
            return Err(DwrError::NonFiniteLinearSystem {
                component: "back-substitution",
                index: i,
            });
        }
        progress.advance(1)?;
    }
    Ok(())
}

/// P1 FEM solve of `−z″ = w` (zero Dirichlet BC) on `mesh`, with
/// `w = 1` on `[w_lo, w_hi]` — deterministic Thomas solve.
fn dual_solve(
    mesh: &[f64],
    w_lo: f64,
    w_hi: f64,
    progress: &mut WorkProgress,
    cx: &Cx<'_>,
) -> Result<Vec<f64>, DwrError> {
    let n = mesh.len();
    let free = n.saturating_sub(2);
    dwr_checkpoint(DWR_DUAL_ASSEMBLY_PHASE, progress, cx)?;
    if free == 0 {
        return zeroed_vec(n, DWR_DUAL_PUBLICATION_PHASE, progress, cx);
    }
    let mut sub = zeroed_vec(free, DWR_DUAL_ASSEMBLY_PHASE, progress, cx)?;
    let mut diag = zeroed_vec(free, DWR_DUAL_ASSEMBLY_PHASE, progress, cx)?;
    let mut sup = zeroed_vec(free, DWR_DUAL_ASSEMBLY_PHASE, progress, cx)?;
    let mut rhs = zeroed_vec(free, DWR_DUAL_ASSEMBLY_PHASE, progress, cx)?;
    let interior_nodes = 1..=free;
    for e in 0..n - 1 {
        poll_dwr_scan(e, DWR_DUAL_ASSEMBLY_PHASE, progress, cx)?;
        let h = mesh[e + 1] - mesh[e];
        let k = 1.0 / h;
        if !h.is_finite() || h <= 0.0 || !k.is_finite() {
            return Err(DwrError::NonFiniteReciprocal {
                cell: e / 2,
                refined_half: Some((e % 2) as u8),
            });
        }
        for (a, b, v) in [(e, e, k), (e + 1, e + 1, k), (e, e + 1, -k), (e + 1, e, -k)] {
            if interior_nodes.contains(&a) && interior_nodes.contains(&b) {
                let (i, j) = (a - 1, b - 1);
                if i == j {
                    diag[i] += v;
                } else if j == i + 1 {
                    sup[i] += v;
                } else {
                    sub[i] += v;
                }
            }
        }
        // Load: use the closed-form P1 integral on the element/window
        // intersection, evaluated in binary64. Sampling the discontinuous
        // indicator at whole-cell Gauss nodes can alias a narrow
        // non-mesh-aligned window to zero.
        for (node, endpoint_values) in [(e, [1.0, 0.0]), (e + 1, [0.0, 1.0])] {
            if interior_nodes.contains(&node) {
                let contribution = integrate_clipped_p1(
                    mesh[e],
                    mesh[e + 1],
                    endpoint_values,
                    w_lo,
                    w_hi,
                    "dual load",
                    e,
                )?;
                let updated = rhs[node - 1] + contribution;
                if !updated.is_finite() {
                    return Err(DwrError::NonFiniteLinearSystem {
                        component: "assembled right-hand side",
                        index: node - 1,
                    });
                }
                rhs[node - 1] = updated;
            }
        }
        progress.advance(1)?;
    }
    thomas_solve(&sub, &diag, &sup, &mut rhs, progress, cx)?;
    let mut z = zeroed_vec(n, DWR_DUAL_PUBLICATION_PHASE, progress, cx)?;
    dwr_checkpoint(DWR_DUAL_PUBLICATION_PHASE, progress, cx)?;
    for (index, value) in rhs.iter().copied().enumerate() {
        poll_dwr_scan(index, DWR_DUAL_PUBLICATION_PHASE, progress, cx)?;
        z[index + 1] = value;
        progress.advance(1)?;
    }
    Ok(z)
}

fn integrate_clipped_p1(
    x0: f64,
    x1: f64,
    endpoint_values: [f64; 2],
    w_lo: f64,
    w_hi: f64,
    quantity: &'static str,
    cell: usize,
) -> Result<f64, DwrError> {
    let a = x0.max(w_lo);
    let b = x1.min(w_hi);
    if a >= b {
        return Ok(0.0);
    }
    let h = x1 - x0;
    let interpolate = |x: f64| {
        let xi = (x - x0) / h;
        (1.0 - xi) * endpoint_values[0] + xi * endpoint_values[1]
    };
    let value_a = interpolate(a);
    let value_b = interpolate(b);
    let average = f64::midpoint(value_a, value_b);
    let width = b - a;
    let integral = width * average;
    if !value_a.is_finite()
        || !value_b.is_finite()
        || !average.is_finite()
        || !width.is_finite()
        || width <= 0.0
        || !integral.is_finite()
    {
        return Err(non_finite_derived(quantity, Some(cell)));
    }
    // `average != 0.0` is not a sufficient underflow guard: midpoint(0,
    // MIN_POSITIVE_SUBNORMAL) itself rounds to zero. Nor may rounded clipped
    // endpoint values authenticate their own cancellation: interpolation can
    // round a nonzero exact P1 value to zero or manufacture additive inverses.
    // Admit a zero only when it follows directly from the public P1 data:
    // either the field is identically zero, or its original endpoint values
    // are exact additive inverses and the clipped subinterval is exactly
    // symmetric about the cell midpoint. `two_diff` compares the complete
    // error-free expansions rather than rounded offsets. Other zero results
    // remain fail-closed.
    let zero_is_proven = (endpoint_values[0] == 0.0 && endpoint_values[1] == 0.0)
        || (endpoint_values[0] == -endpoint_values[1] && two_diff(a, x0) == two_diff(x1, b));
    if integral == 0.0 && !zero_is_proven {
        return Err(DwrError::UnresolvedZeroIntegral { quantity, cell });
    }
    Ok(integral)
}

fn refine(
    mesh: &[f64],
    refined_nodes: usize,
    progress: &mut WorkProgress,
    cx: &Cx<'_>,
) -> Result<Vec<f64>, DwrError> {
    if refined_node_count(mesh.len())? != refined_nodes {
        return Err(DwrError::RefinedMeshSizeOverflow {
            mesh_nodes: mesh.len(),
        });
    }
    let mut out = Vec::new();
    dwr_checkpoint(DWR_REFINE_PHASE, progress, cx)?;
    out.try_reserve_exact(refined_nodes)
        .map_err(|_| DwrError::AllocationFailed {
            phase: DWR_REFINE_PHASE,
            elements: refined_nodes,
        })?;
    dwr_checkpoint(DWR_REFINE_PHASE, progress, cx)?;
    for e in 0..mesh.len() - 1 {
        poll_dwr_scan(e, DWR_REFINE_PHASE, progress, cx)?;
        let midpoint = f64::midpoint(mesh[e], mesh[e + 1]);
        if !(mesh[e] < midpoint && midpoint < mesh[e + 1]) {
            return Err(DwrError::NonInteriorMidpoint { cell: e });
        }
        out.push(mesh[e]);
        out.push(midpoint);
        progress.advance(1)?;
    }
    let Some(&last) = mesh.last() else {
        return Err(DwrError::MeshNodeCount {
            count: 0,
            minimum: 2,
            maximum: MAX_DWR_MESH_NODES,
        });
    };
    out.push(last);
    Ok(out)
}

/// Run the 1-D goal-oriented estimate (see [`DwrOutput`]).
///
/// # Errors
/// [`DwrError`] when public inputs exceed the bounded execution envelope or any
/// input/derived numerical value cannot support a finite Estimated diagnostic.
pub fn dwr_integral_qoi(
    problem: &MmsProblem,
    candidate: &[f64],
    w_lo: f64,
    w_hi: f64,
    cx: &Cx<'_>,
) -> Result<DwrOutput, DwrError> {
    let mesh = problem.mesh();
    let f = problem.forcing();
    let plan = DwrWorkPlan::preflight(
        mesh.len(),
        candidate.len(),
        problem.exact_solution().coefficients().len(),
        f.coefficients().len(),
        problem.canonical_bytes().len(),
    )?;
    let mut progress = WorkProgress::new(plan.planned_work_units, cx);
    dwr_checkpoint(DWR_INITIAL_PHASE, &mut progress, cx)?;
    let validated = validate_dwr_inputs(problem, candidate, w_lo, w_hi, &mut progress, cx)?;

    dwr_checkpoint(DWR_VALIDATE_FORCING_PHASE, &mut progress, cx)?;
    for (index, value) in f.coefficients().iter().enumerate() {
        poll_dwr_scan(index, DWR_VALIDATE_FORCING_PHASE, &mut progress, cx)?;
        if !value.is_finite() {
            return Err(non_finite_derived("forcing coefficient", Some(index)));
        }
        progress.advance(1)?;
    }
    if progress.completed_work_units != plan.validation_work_units {
        return Err(DwrError::WorkPlanMismatch {
            phase: DWR_VALIDATE_FORCING_PHASE,
            completed_work_units: progress.completed_work_units,
            planned_work_units: plan.validation_work_units,
        });
    }

    // J(u_h): evaluate the closed-form P1 integral in binary64 on every
    // cell/window intersection. Whole-cell indicator sampling can miss narrow
    // windows.
    dwr_checkpoint(DWR_PRIMAL_PHASE, &mut progress, cx)?;
    let mut j_primal = 0.0f64;
    for e in 0..mesh.len() - 1 {
        poll_dwr_scan(e, DWR_PRIMAL_PHASE, &mut progress, cx)?;
        let contribution = integrate_clipped_p1(
            mesh[e],
            mesh[e + 1],
            [candidate[e], candidate[e + 1]],
            w_lo,
            w_hi,
            "primal QoI",
            e,
        )?;
        let updated = j_primal + contribution;
        if !updated.is_finite() {
            return Err(non_finite_derived("primal QoI", Some(e)));
        }
        j_primal = updated;
        progress.advance(1)?;
    }
    // Enriched dual on the refined mesh.
    let fine = refine(mesh, validated.refined_nodes, &mut progress, cx)?;
    let z = dual_solve(&fine, w_lo, w_hi, &mut progress, cx)?;
    // Coarse-node interpolant of z, subtracted (Galerkin orthogonality
    // makes the coarse part vanish; the fine remainder drives η).
    let mut eta = 0.0f64;
    let mut indicators = zeroed_vec(mesh.len() - 1, DWR_RESIDUAL_PHASE, &mut progress, cx)?;
    dwr_checkpoint(DWR_RESIDUAL_PHASE, &mut progress, cx)?;
    for e in 0..mesh.len() - 1 {
        poll_dwr_scan(e, DWR_RESIDUAL_PHASE, &mut progress, cx)?;
        let (x0, x1) = (mesh[e], mesh[e + 1]);
        let slope = (candidate[e + 1] - candidate[e]) / (x1 - x0);
        if !slope.is_finite() {
            return Err(non_finite_derived("primal slope", Some(e)));
        }
        let (z0, z1) = (z[2 * e], z[2 * e + 2]);
        let mut local = 0.0f64;
        // Two fine halves of the coarse element.
        for half in 0..2usize {
            let (fa, fb) = (fine[2 * e + half], fine[2 * e + half + 1]);
            let (za, zb) = (z[2 * e + half], z[2 * e + half + 1]);
            let zslope = (zb - za) / (fb - fa);
            // Coarse interpolant of z on this fine piece.
            let islope = (z1 - z0) / (x1 - x0);
            if !zslope.is_finite() || !islope.is_finite() {
                return Err(non_finite_derived("dual slope", Some(e)));
            }
            for (gx, gw) in gauss5(fa, fb) {
                let xi_f = (gx - fa) / (fb - fa);
                let zf = (1.0 - xi_f) * za + xi_f * zb;
                let zi = z0 + (gx - x0) * islope;
                // r(v) = ∫ f v − ∫ u_h′ v′ with v = z_f − I_h z_f.
                let forcing = f.eval(gx);
                let contribution = gw * (forcing * (zf - zi) - slope * (zslope - islope));
                let updated = local + contribution;
                if !gx.is_finite()
                    || !gw.is_finite()
                    || !xi_f.is_finite()
                    || !zf.is_finite()
                    || !zi.is_finite()
                    || !forcing.is_finite()
                    || !contribution.is_finite()
                    || !updated.is_finite()
                {
                    return Err(non_finite_derived("dual-weighted residual", Some(e)));
                }
                local = updated;
            }
        }
        let updated_eta = eta + local;
        if !updated_eta.is_finite() {
            return Err(non_finite_derived("global DWR estimate", Some(e)));
        }
        eta = updated_eta;
        indicators[e] = local.abs();
        progress.advance(1)?;
    }
    if !j_primal.is_finite() || !eta.is_finite() {
        return Err(non_finite_derived("DWR output", None));
    }

    dwr_checkpoint(DWR_OUTPUT_VALIDATE_PHASE, &mut progress, cx)?;
    for (index, indicator) in indicators.iter().enumerate() {
        poll_dwr_scan(index, DWR_OUTPUT_VALIDATE_PHASE, &mut progress, cx)?;
        if !indicator.is_finite() {
            return Err(non_finite_derived("DWR output", Some(index)));
        }
        progress.advance(1)?;
    }

    dwr_checkpoint(DWR_IDENTITY_PHASE, &mut progress, cx)?;
    let mut hasher = fs_blake3::Blake3::new();
    hash_execution_header(
        &mut hasher,
        DWR_OUTPUT_IDENTITY_SCHEMA,
        &plan.identity_fields(),
        cx,
    );
    hasher.update(&problem.identity().version().to_le_bytes());
    hasher.update(&problem.identity().root().to_le_bytes());
    let problem_identity = problem.canonical_bytes();
    let problem_identity_len =
        u64::try_from(problem_identity.len()).map_err(|_| DwrError::WorkPlanOverflow)?;
    hasher.update(&problem_identity_len.to_le_bytes());
    for chunk in problem_identity.chunks(DWR_POLL_STRIDE_ITEMS) {
        dwr_checkpoint(DWR_IDENTITY_PHASE, &mut progress, cx)?;
        hasher.update(chunk);
        progress.advance(u128::try_from(chunk.len()).map_err(|_| DwrError::WorkPlanOverflow)?)?;
    }
    hasher.update(validated.candidate_digest.as_bytes());
    hasher.update(&w_lo.to_bits().to_le_bytes());
    hasher.update(&w_hi.to_bits().to_le_bytes());
    hasher.update(&j_primal.to_bits().to_le_bytes());
    hasher.update(&eta.to_bits().to_le_bytes());
    // Reset the documented 256-item cancellation bound between the
    // variable-length canonical problem bytes and the independent indicator
    // stream. Without this boundary, a final problem chunk and the first 256
    // indicators could consume almost 512 bounded items under one poll.
    dwr_checkpoint(DWR_IDENTITY_PHASE, &mut progress, cx)?;
    for (index, indicator) in indicators.iter().enumerate() {
        poll_dwr_scan(index, DWR_IDENTITY_PHASE, &mut progress, cx)?;
        hasher.update(&indicator.to_bits().to_le_bytes());
        progress.advance(1)?;
    }
    let output = DwrOutput {
        j_primal,
        eta,
        indicators,
        evidence_identity: DwrEvidenceIdentity(hasher.finalize()),
    };
    progress.advance(plan.finalization_work_units)?;
    progress.finish_dwr(DWR_PUBLICATION_PHASE)?;
    dwr_checkpoint(DWR_PUBLICATION_PHASE, &mut progress, cx)?;
    Ok(output)
}

#[cfg(test)]
mod execution_tests {
    use super::*;

    fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
        let gate = fs_exec::CancelGate::new_clock_free();
        let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
        pool.scope(|arena| {
            let cx = Cx::new(
                &gate,
                arena,
                fs_exec::StreamKey {
                    seed: 7,
                    kernel_id: 11,
                    tile: 13,
                    iteration: 17,
                },
                fs_exec::Budget::INFINITE,
                fs_exec::ExecMode::Deterministic,
            );
            f(&cx)
        })
    }

    #[test]
    fn zero_interior_linear_and_dual_systems_are_total() {
        with_cx(|cx| {
            let mut progress = WorkProgress::new(2, cx);
            let mut rhs = Vec::new();
            thomas_solve(&[], &[], &[], &mut rhs, &mut progress, cx)
                .expect("empty system is solved");
            let dual =
                dual_solve(&[0.0, 1.0], 0.0, 1.0, &mut progress, cx).expect("boundary-only dual");
            assert_eq!(dual, vec![0.0, 0.0]);
        });
    }

    #[test]
    fn hostile_nonuniform_mesh_refuses_a_rounded_negative_spd_pivot() {
        with_cx(|cx| {
            // The exact 1-D Dirichlet stiffness matrix is SPD, hence every
            // exact Thomas pivot is positive. On this admitted finite mesh the
            // binary64 recurrence suffers enough cancellation to produce a
            // negative later pivot. Continuing would silently solve a
            // numerically indefinite system and corrupt the DWR estimate.
            let coarse = [
                -1.0e175, -1.0e80, -1.0e-225, 1.0e-280, 1.0e-255, 1.0e-250, 1.0e-235, 1.0e-65,
            ];
            let mut progress = WorkProgress::new(100_000, cx);
            let fine = refine(
                &coarse,
                refined_node_count(coarse.len()).expect("small hostile mesh"),
                &mut progress,
                cx,
            )
            .expect("hostile mesh still has finite interior midpoints");
            // Use an admitted finite window disjoint from the mesh. Its zero
            // load isolates stiffness elimination; a full-domain load on this
            // deliberately extreme mesh overflows the forward right-hand side
            // before the hostile pivot is reached.
            let error = dual_solve(&fine, 1.0, 2.0, &mut progress, cx)
                .expect_err("a non-positive SPD elimination pivot must fail closed");
            assert!(matches!(error, DwrError::InvalidLinearPivot { .. }));
        });
    }

    #[test]
    fn incomplete_work_plans_fail_closed_before_publication() {
        with_cx(|cx| {
            let progress = WorkProgress::new(2, cx);
            assert_eq!(
                progress.finish_dwr("test.dwr-publish"),
                Err(DwrError::WorkPlanMismatch {
                    phase: "test.dwr-publish",
                    completed_work_units: 0,
                    planned_work_units: 2,
                })
            );
            assert_eq!(
                progress.finish_bracket("test.bracket-publish"),
                Err(BracketError::WorkPlanMismatch {
                    phase: "test.bracket-publish",
                    completed_work_units: 0,
                    planned_work_units: 2,
                })
            );
        });
    }

    #[test]
    fn bracket_header_binds_each_nested_plan_field_and_verifier_policy_field() {
        with_cx(|cx| {
            let mut fields = [0_u128; 23];
            for (index, field) in fields.iter_mut().enumerate() {
                *field = u128::try_from(index).expect("small retained field index") + 1;
            }
            let root = |fields: &[u128; 23], policy: VerifierPolicyIdentity| {
                let mut hasher = fs_blake3::Blake3::new();
                hash_bracket_execution_header(&mut hasher, fields, cx, policy);
                hasher.finalize()
            };
            let baseline = root(&fields, CURRENT_VERIFIER_POLICY_IDENTITY);
            assert_eq!(
                baseline.to_hex(),
                "e53dcf11159b73d5516753630b8f0fd6572d08a8f0e60a6b6fba9672e6ceedac",
                "fixed bracket work-policy preimage changed"
            );
            for index in (4..=9).chain(14..=19) {
                let mut changed = fields;
                changed[index] += 1;
                assert_ne!(
                    baseline,
                    root(&changed, CURRENT_VERIFIER_POLICY_IDENTITY),
                    "nested verifier plan field {index} was omitted"
                );
            }
            for changed_policy in [
                VerifierPolicyIdentity {
                    work_plan_version: CURRENT_VERIFIER_POLICY_IDENTITY.work_plan_version + 1,
                    ..CURRENT_VERIFIER_POLICY_IDENTITY
                },
                VerifierPolicyIdentity {
                    poll_policy_version: CURRENT_VERIFIER_POLICY_IDENTITY.poll_policy_version + 1,
                    ..CURRENT_VERIFIER_POLICY_IDENTITY
                },
                VerifierPolicyIdentity {
                    poll_stride_work_units: CURRENT_VERIFIER_POLICY_IDENTITY.poll_stride_work_units
                        + 1,
                    ..CURRENT_VERIFIER_POLICY_IDENTITY
                },
            ] {
                assert_ne!(baseline, root(&fields, changed_policy));
            }
        });
    }

    fn large_verifier_plan() -> VerifierWorkPlan {
        let mesh: Vec<f64> = (0_u32..=512)
            .map(|index| f64::from(index) / 512.0)
            .collect();
        let candidate = vec![0.0; mesh.len()];
        let problem = MmsProblem::new(
            "verifier-observation-policy",
            fs_verify::fem1d::Poly::new(vec![0.0]).expect("zero polynomial"),
            mesh,
        )
        .expect("admitted verifier observation fixture");
        VerifierWorkPlan::for_inputs(&problem, &candidate).expect("admitted verifier work plan")
    }

    #[test]
    #[allow(clippy::too_many_lines)] // One retained trace covers every fail-closed policy edge.
    fn nested_verifier_observation_policy_fails_closed() {
        let plan = large_verifier_plan();
        with_cx(|cx| {
            let entry = VerifierProgress {
                kind: VerifierCheckpointKind::PhaseEntry,
                phase: VerifierPhase::Validation,
                completed_work_units: 0,
                planned_work_units: plan.planned_work_units(),
            };

            let mut state = VerifierObservationState::default();
            let mut outer = WorkProgress::new(plan.planned_work_units(), cx);
            observe_verifier_progress(
                "primal",
                BRACKET_PRIMAL_VERIFY_PHASE,
                plan,
                entry,
                &mut state,
                &mut outer,
                cx,
            )
            .expect("canonical validation entry");
            let premature_tightness = VerifierProgress {
                kind: VerifierCheckpointKind::PhaseEntry,
                phase: VerifierPhase::Tightness,
                completed_work_units: 1,
                planned_work_units: plan.planned_work_units(),
            };
            assert!(matches!(
                observe_verifier_progress(
                    "primal",
                    BRACKET_PRIMAL_VERIFY_PHASE,
                    plan,
                    premature_tightness,
                    &mut state,
                    &mut outer,
                    cx,
                ),
                Err(BracketError::WorkPlanMismatch { .. })
            ));
            let omitted_boundary = VerifierProgress {
                kind: VerifierCheckpointKind::RefusalFlush,
                phase: VerifierPhase::Validation,
                completed_work_units: VERIFIER_POLL_STRIDE_WORK_UNITS,
                planned_work_units: plan.planned_work_units(),
            };
            assert!(matches!(
                observe_verifier_progress(
                    "primal",
                    BRACKET_PRIMAL_VERIFY_PHASE,
                    plan,
                    omitted_boundary,
                    &mut state,
                    &mut outer,
                    cx,
                ),
                Err(BracketError::WorkPlanMismatch { .. })
            ));
            let missing_boundary = VerifierProgress {
                kind: VerifierCheckpointKind::WorkBoundary,
                phase: VerifierPhase::Validation,
                completed_work_units: 512,
                planned_work_units: plan.planned_work_units(),
            };
            assert!(matches!(
                observe_verifier_progress(
                    "primal",
                    BRACKET_PRIMAL_VERIFY_PHASE,
                    plan,
                    missing_boundary,
                    &mut state,
                    &mut outer,
                    cx,
                ),
                Err(BracketError::WorkPlanMismatch { .. })
            ));
            let mut skipped_boundary_state = VerifierObservationState {
                last_completed_work_units: 200,
                next_phase_index: 1,
                active_phase: Some(VerifierPhase::Validation),
                saw_refusal_flush: false,
                saw_publication: false,
            };
            let mut outer = WorkProgress::new(plan.planned_work_units(), cx);
            let callback_after_skipped_boundary = VerifierProgress {
                kind: VerifierCheckpointKind::RefusalFlush,
                phase: VerifierPhase::Validation,
                completed_work_units: 300,
                planned_work_units: plan.planned_work_units(),
            };
            assert!(matches!(
                observe_verifier_progress(
                    "primal",
                    BRACKET_PRIMAL_VERIFY_PHASE,
                    plan,
                    callback_after_skipped_boundary,
                    &mut skipped_boundary_state,
                    &mut outer,
                    cx,
                ),
                Err(BracketError::WorkPlanMismatch { .. })
            ));

            let mut state = VerifierObservationState::default();
            let mut outer = WorkProgress::new(plan.planned_work_units(), cx);
            observe_verifier_progress(
                "primal",
                BRACKET_PRIMAL_VERIFY_PHASE,
                plan,
                entry,
                &mut state,
                &mut outer,
                cx,
            )
            .expect("canonical validation entry");
            let malformed_boundary = VerifierProgress {
                completed_work_units: 255,
                ..missing_boundary
            };
            assert!(matches!(
                observe_verifier_progress(
                    "primal",
                    BRACKET_PRIMAL_VERIFY_PHASE,
                    plan,
                    malformed_boundary,
                    &mut state,
                    &mut outer,
                    cx,
                ),
                Err(BracketError::WorkPlanMismatch { .. })
            ));
            let canonical_boundary = VerifierProgress {
                completed_work_units: VERIFIER_POLL_STRIDE_WORK_UNITS,
                ..missing_boundary
            };
            observe_verifier_progress(
                "primal",
                BRACKET_PRIMAL_VERIFY_PHASE,
                plan,
                canonical_boundary,
                &mut state,
                &mut outer,
                cx,
            )
            .expect("canonical global work boundary");

            let missing_publication = VerifierObservationState {
                last_completed_work_units: plan.planned_work_units(),
                next_phase_index: 5,
                active_phase: Some(VerifierPhase::Finalization),
                saw_refusal_flush: false,
                saw_publication: false,
            };
            assert!(!verifier_observation_complete(
                &missing_publication,
                plan,
                false
            ));
            let mut published = VerifierObservationState {
                last_completed_work_units: plan.planned_work_units() - 1,
                ..missing_publication
            };
            let mut outer = WorkProgress::new(plan.planned_work_units(), cx);
            observe_verifier_progress(
                "primal",
                BRACKET_PRIMAL_VERIFY_PHASE,
                plan,
                VerifierProgress {
                    kind: VerifierCheckpointKind::Publication,
                    phase: VerifierPhase::Finalization,
                    completed_work_units: plan.planned_work_units(),
                    planned_work_units: plan.planned_work_units(),
                },
                &mut published,
                &mut outer,
                cx,
            )
            .expect("canonical publication callback");
            assert!(verifier_observation_complete(&published, plan, false));
        });
    }

    #[test]
    fn hostile_maximum_and_plus_one_work_shapes_are_checked_without_allocation() {
        assert!(matches!(
            refined_node_count(usize::MAX),
            Err(DwrError::RefinedMeshSizeOverflow {
                mesh_nodes: usize::MAX
            })
        ));
        let maximum = DwrWorkPlan::preflight(
            MAX_DWR_MESH_NODES,
            MAX_DWR_MESH_NODES,
            MAX_DWR_POLY_COEFFICIENTS,
            MAX_DWR_POLY_COEFFICIENTS,
            MAX_FEM1D_PROBLEM_CANONICAL_IDENTITY_BYTES,
        )
        .expect("maximum public shape is admitted without allocation");
        assert_eq!(maximum.planned_work_units, 45_004_590);
        assert!(maximum.planned_work_units <= MAX_DWR_WORK_UNITS as u128);
        assert!(matches!(
            DwrWorkPlan::preflight(
                MAX_DWR_MESH_NODES + 1,
                MAX_DWR_MESH_NODES + 1,
                MAX_DWR_POLY_COEFFICIENTS,
                MAX_DWR_POLY_COEFFICIENTS,
                MAX_FEM1D_PROBLEM_CANONICAL_IDENTITY_BYTES,
            ),
            Err(DwrError::MeshNodeCount { .. })
        ));
        assert!(matches!(
            DwrWorkPlan::preflight(2, 2, MAX_DWR_POLY_COEFFICIENTS + 1, 1, 0),
            Err(DwrError::PolynomialCoefficientCount { .. })
        ));
        assert_eq!(DWR_POLL_STRIDE_ITEMS, 256);
        assert_eq!(DWR_POLL_POLICY_VERSION, 3);
    }
}
