//! G5 study-scale replay and seeded-failure self-test for trust-region Newton.
//!
//! The production driver solves the existing six-dimensional Rosenbrock
//! fixture with an exact matrix-free Hessian-vector callback. The receipt binds
//! every objective and Hessian-vector callback input/output under one global
//! ordinal, reconstructs each outer iteration's Steihaug and final-model roles,
//! and binds every public `TrustRegionReport` field, configuration value, and
//! independently reconstructed nonpositive Steihaug-curvature witness. A
//! separate algebraic oracle recomputes every Rosenbrock value, gradient, and
//! Hessian-vector product. A same-input repeat must reproduce the receipt byte
//! for byte. Deterministic red mutations cover the returned decision,
//! objective/gradient evidence, trust configuration/accounting, callback
//! ordinal/role/segment metadata, and negative-curvature claim; stale and
//! correctly resealed forms fail through distinct typed refusal paths.
//!
//! This is one objective/Hessian pair. It does not claim all objectives,
//! approximate-Hessian parity, cancellation, checkpointing, cross-ISA equality,
//! ledger persistence, or performance.

use core::cell::RefCell;

use fs_ascent::{TrustRegionReport, trust_region_newton};
use fs_obs::ident::{IdentityBuilder, ReplayIdentity, check_version};
use fs_obs::{Emitter, EventKind, Severity};

const SUITE: &str = "fs-ascent/trust-region-study-replay";
// `ConformanceCase` requires a seed field even for non-random fixtures. This is
// a schema sentinel, not a claim that randomness influences the study.
const EVENT_SEED_SENTINEL: u64 = 0;
const MUTATION_SEED: u64 = 0x5452_5553_545F_5244;
const DIMENSION: usize = 6;
const GRADIENT_TOLERANCE: f64 = 1e-7;
const MAX_ITERATIONS: usize = 300;
const OBJECTIVE_ORACLE_VERSION: &str = "rosenbrock-chain-independent-v1";
const HESSIAN_ORACLE_VERSION: &str = "rosenbrock-tridiagonal-independent-v1";
const TRUST_CONFIG_VERSION: &str = "trust-region-newton-steihaug-v1";
const INITIAL_TRUST_RADIUS: f64 = 1.0;
const STEIHAUG_RELATIVE_TOLERANCE: f64 = 1e-8;
const STEIHAUG_GRADIENT_NORM_FLOOR: f64 = 1e-30;
const STEIHAUG_MAX_STEPS_PER_DIMENSION: usize = 2;
const SHRINK_RATIO_THRESHOLD: f64 = 0.25;
const GROW_RATIO_THRESHOLD: f64 = 0.75;
const ACCEPT_RATIO_THRESHOLD: f64 = 1e-4;
const RADIUS_SHRINK_FACTOR: f64 = 0.25;
const RADIUS_GROW_FACTOR: f64 = 2.0;
const BOUNDARY_RELATIVE_TOLERANCE: f64 = 1e-10;
const MODEL_DECREASE_ZERO_THRESHOLD: f64 = 1e-300;
const MAX_TRUST_RADIUS: f64 = 1e8;
const MIN_TRUST_RADIUS: f64 = 1e-14;
const NEGATIVE_CURVATURE_THRESHOLD: f64 = 0.0;
const SEMANTIC_COMPARISON_EPSILON_MULTIPLIER: f64 = 4096.0;
const FINAL_OBJECTIVE_LIMIT: f64 = 1e-10;
const FINAL_COORDINATE_ERROR_LIMIT: f64 = 1e-4;
const START: [f64; DIMENSION] = [-1.2; DIMENSION];

#[derive(Clone, Debug, PartialEq, Eq)]
struct ObjectiveCall {
    ordinal: usize,
    point_bits: Vec<u64>,
    objective_bits: u64,
    gradient_bits: Vec<u64>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum HessianVectorRole {
    Steihaug,
    ModelDecrease,
}

impl HessianVectorRole {
    fn identity_tag(self) -> &'static str {
        match self {
            Self::Steihaug => "steihaug-search-direction",
            Self::ModelDecrease => "final-model-decrease-step",
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct HessianVectorCall {
    ordinal: usize,
    outer_iteration: usize,
    role: HessianVectorRole,
    point_bits: Vec<u64>,
    direction_bits: Vec<u64>,
    product_bits: Vec<u64>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct CapturedHessianVectorCall {
    ordinal: usize,
    point_bits: Vec<u64>,
    direction_bits: Vec<u64>,
    product_bits: Vec<u64>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
enum CapturedCall {
    Objective(ObjectiveCall),
    HessianVector(CapturedHessianVectorCall),
}

#[derive(Debug)]
struct RunRecord {
    report: TrustRegionReport,
    objective_calls: Vec<ObjectiveCall>,
    hessian_vector_calls: Vec<HessianVectorCall>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct TrustRegionConfigPayload {
    dimension: usize,
    gradient_tolerance_bits: u64,
    maximum_iterations: usize,
    initial_radius_bits: u64,
    steihaug_relative_tolerance_bits: u64,
    steihaug_gradient_norm_floor_bits: u64,
    steihaug_max_steps_per_dimension: usize,
    shrink_ratio_threshold_bits: u64,
    grow_ratio_threshold_bits: u64,
    accept_ratio_threshold_bits: u64,
    radius_shrink_factor_bits: u64,
    radius_grow_factor_bits: u64,
    boundary_relative_tolerance_bits: u64,
    model_decrease_zero_threshold_bits: u64,
    maximum_radius_bits: u64,
    minimum_radius_bits: u64,
    negative_curvature_threshold_bits: u64,
    semantic_comparison_epsilon_multiplier_bits: u64,
    final_objective_limit_bits: u64,
    final_coordinate_error_limit_bits: u64,
}

impl TrustRegionConfigPayload {
    fn canonical() -> Self {
        Self {
            dimension: DIMENSION,
            gradient_tolerance_bits: GRADIENT_TOLERANCE.to_bits(),
            maximum_iterations: MAX_ITERATIONS,
            initial_radius_bits: INITIAL_TRUST_RADIUS.to_bits(),
            steihaug_relative_tolerance_bits: STEIHAUG_RELATIVE_TOLERANCE.to_bits(),
            steihaug_gradient_norm_floor_bits: STEIHAUG_GRADIENT_NORM_FLOOR.to_bits(),
            steihaug_max_steps_per_dimension: STEIHAUG_MAX_STEPS_PER_DIMENSION,
            shrink_ratio_threshold_bits: SHRINK_RATIO_THRESHOLD.to_bits(),
            grow_ratio_threshold_bits: GROW_RATIO_THRESHOLD.to_bits(),
            accept_ratio_threshold_bits: ACCEPT_RATIO_THRESHOLD.to_bits(),
            radius_shrink_factor_bits: RADIUS_SHRINK_FACTOR.to_bits(),
            radius_grow_factor_bits: RADIUS_GROW_FACTOR.to_bits(),
            boundary_relative_tolerance_bits: BOUNDARY_RELATIVE_TOLERANCE.to_bits(),
            model_decrease_zero_threshold_bits: MODEL_DECREASE_ZERO_THRESHOLD.to_bits(),
            maximum_radius_bits: MAX_TRUST_RADIUS.to_bits(),
            minimum_radius_bits: MIN_TRUST_RADIUS.to_bits(),
            negative_curvature_threshold_bits: NEGATIVE_CURVATURE_THRESHOLD.to_bits(),
            semantic_comparison_epsilon_multiplier_bits: SEMANTIC_COMPARISON_EPSILON_MULTIPLIER
                .to_bits(),
            final_objective_limit_bits: FINAL_OBJECTIVE_LIMIT.to_bits(),
            final_coordinate_error_limit_bits: FINAL_COORDINATE_ERROR_LIMIT.to_bits(),
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct NegativeCurvatureWitness {
    hessian_call_index: usize,
    call_ordinal: usize,
    outer_iteration: usize,
    role: HessianVectorRole,
    point_bits: Vec<u64>,
    direction_bits: Vec<u64>,
    independently_recomputed_product_bits: Vec<u64>,
    recorded_quadratic_form_bits: u64,
    quadratic_form_bits: u64,
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct ReceiptPayload {
    suite: &'static str,
    mutation_seed: u64,
    event_seed_sentinel: u64,
    fs_math_version: &'static str,
    fs_obs_version: &'static str,
    objective_oracle_version: &'static str,
    hessian_oracle_version: &'static str,
    trust_config_version: &'static str,
    config: TrustRegionConfigPayload,
    start_bits: Vec<u64>,
    objective_calls: Vec<ObjectiveCall>,
    hessian_vector_calls: Vec<HessianVectorCall>,
    report_x_bits: Vec<u64>,
    report_f_bits: u64,
    report_gradient_norm_bits: u64,
    report_iterations: usize,
    report_evaluations: usize,
    report_hessian_vector_evaluations: usize,
    report_negative_curvature_hits: usize,
    negative_curvature_witness: NegativeCurvatureWitness,
}

impl ReceiptPayload {
    #[allow(clippy::too_many_lines)] // One canonical field-order audit for the complete receipt.
    fn identity(&self) -> ReplayIdentity {
        let mut builder = IdentityBuilder::new("fs-ascent-trust-region-study-receipt-v3")
            .str("suite", self.suite)
            .u64("mutation-seed", self.mutation_seed)
            .u64("event-seed-sentinel", self.event_seed_sentinel)
            .str("fs-ascent-version", fs_ascent::VERSION)
            .str("fs-math-version", self.fs_math_version)
            .str("fs-obs-version", self.fs_obs_version)
            .str("engine", "trust_region_newton/Steihaug-CG")
            .str("objective", "Rosenbrock-chain")
            .str("hessian-vector", "exact-analytic")
            .str("randomness", "none")
            .str("objective-oracle-version", self.objective_oracle_version)
            .str("hessian-oracle-version", self.hessian_oracle_version)
            .str("trust-config-version", self.trust_config_version)
            .u64("dimension", self.config.dimension as u64)
            .u64(
                "gradient-tolerance-bits",
                self.config.gradient_tolerance_bits,
            )
            .u64("maximum-iterations", self.config.maximum_iterations as u64)
            .u64("initial-radius-bits", self.config.initial_radius_bits)
            .u64(
                "steihaug-relative-tolerance-bits",
                self.config.steihaug_relative_tolerance_bits,
            )
            .u64(
                "steihaug-gradient-norm-floor-bits",
                self.config.steihaug_gradient_norm_floor_bits,
            )
            .u64(
                "steihaug-max-steps-per-dimension",
                self.config.steihaug_max_steps_per_dimension as u64,
            )
            .u64(
                "shrink-ratio-threshold-bits",
                self.config.shrink_ratio_threshold_bits,
            )
            .u64(
                "grow-ratio-threshold-bits",
                self.config.grow_ratio_threshold_bits,
            )
            .u64(
                "accept-ratio-threshold-bits",
                self.config.accept_ratio_threshold_bits,
            )
            .u64(
                "radius-shrink-factor-bits",
                self.config.radius_shrink_factor_bits,
            )
            .u64(
                "radius-grow-factor-bits",
                self.config.radius_grow_factor_bits,
            )
            .u64(
                "boundary-relative-tolerance-bits",
                self.config.boundary_relative_tolerance_bits,
            )
            .u64(
                "model-decrease-zero-threshold-bits",
                self.config.model_decrease_zero_threshold_bits,
            )
            .u64("maximum-radius-bits", self.config.maximum_radius_bits)
            .u64("minimum-radius-bits", self.config.minimum_radius_bits)
            .u64(
                "negative-curvature-threshold-bits",
                self.config.negative_curvature_threshold_bits,
            )
            .u64(
                "semantic-comparison-epsilon-multiplier-bits",
                self.config.semantic_comparison_epsilon_multiplier_bits,
            )
            .u64(
                "final-objective-limit-bits",
                self.config.final_objective_limit_bits,
            )
            .u64(
                "final-coordinate-error-limit-bits",
                self.config.final_coordinate_error_limit_bits,
            )
            .u64("start-values", self.start_bits.len() as u64);
        for &value_bits in &self.start_bits {
            builder = builder.u64("start-value-bits", value_bits);
        }

        builder = builder.u64("objective-calls", self.objective_calls.len() as u64);
        for call in &self.objective_calls {
            builder = builder
                .u64("objective-call-ordinal", call.ordinal as u64)
                .u64("objective-point-values", call.point_bits.len() as u64);
            for &value_bits in &call.point_bits {
                builder = builder.u64("objective-point-bits", value_bits);
            }
            builder = builder
                .u64("objective-value-bits", call.objective_bits)
                .u64("gradient-values", call.gradient_bits.len() as u64);
            for &value_bits in &call.gradient_bits {
                builder = builder.u64("gradient-value-bits", value_bits);
            }
        }

        builder = builder.u64(
            "hessian-vector-calls",
            self.hessian_vector_calls.len() as u64,
        );
        for call in &self.hessian_vector_calls {
            builder = builder
                .u64("hessian-call-ordinal", call.ordinal as u64)
                .u64("hessian-outer-iteration", call.outer_iteration as u64)
                .str("hessian-call-role", call.role.identity_tag())
                .u64("hessian-point-values", call.point_bits.len() as u64);
            for &value_bits in &call.point_bits {
                builder = builder.u64("hessian-point-bits", value_bits);
            }
            builder = builder.u64("hessian-direction-values", call.direction_bits.len() as u64);
            for &value_bits in &call.direction_bits {
                builder = builder.u64("hessian-direction-bits", value_bits);
            }
            builder = builder.u64("hessian-product-values", call.product_bits.len() as u64);
            for &value_bits in &call.product_bits {
                builder = builder.u64("hessian-product-bits", value_bits);
            }
        }

        builder = builder.u64("report-x-values", self.report_x_bits.len() as u64);
        for &value_bits in &self.report_x_bits {
            builder = builder.u64("report-x-bits", value_bits);
        }
        builder = builder
            .u64("report-objective-bits", self.report_f_bits)
            .u64("report-gradient-norm-bits", self.report_gradient_norm_bits)
            .u64("report-iterations", self.report_iterations as u64)
            .u64("report-evaluations", self.report_evaluations as u64)
            .u64(
                "report-hessian-vector-evaluations",
                self.report_hessian_vector_evaluations as u64,
            )
            .u64(
                "report-negative-curvature-hits",
                self.report_negative_curvature_hits as u64,
            );
        let witness = &self.negative_curvature_witness;
        builder = builder
            .u64(
                "negative-curvature-hessian-call-index",
                witness.hessian_call_index as u64,
            )
            .u64(
                "negative-curvature-call-ordinal",
                witness.call_ordinal as u64,
            )
            .u64(
                "negative-curvature-outer-iteration",
                witness.outer_iteration as u64,
            )
            .str("negative-curvature-call-role", witness.role.identity_tag())
            .u64(
                "negative-curvature-point-values",
                witness.point_bits.len() as u64,
            );
        for &value_bits in &witness.point_bits {
            builder = builder.u64("negative-curvature-point-bits", value_bits);
        }
        builder = builder.u64(
            "negative-curvature-direction-values",
            witness.direction_bits.len() as u64,
        );
        for &value_bits in &witness.direction_bits {
            builder = builder.u64("negative-curvature-direction-bits", value_bits);
        }
        builder = builder.u64(
            "negative-curvature-product-values",
            witness.independently_recomputed_product_bits.len() as u64,
        );
        for &value_bits in &witness.independently_recomputed_product_bits {
            builder = builder.u64("negative-curvature-product-bits", value_bits);
        }
        builder
            .u64(
                "negative-curvature-recorded-quadratic-form-bits",
                witness.recorded_quadratic_form_bits,
            )
            .u64(
                "negative-curvature-quadratic-form-bits",
                witness.quadratic_form_bits,
            )
            .finish()
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct RetainedReceipt {
    payload: ReceiptPayload,
    declared_identity: ReplayIdentity,
}

impl RetainedReceipt {
    fn new(payload: ReceiptPayload) -> Self {
        let declared_identity = payload.identity();
        Self {
            payload,
            declared_identity,
        }
    }

    fn reseal(&mut self) {
        self.declared_identity = self.payload.identity();
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum SemanticRefusal {
    FixtureMetadataMismatch,
    TrustConfigurationMismatch,
    DimensionMismatch,
    NonFiniteEvidence,
    ObjectiveValueMismatch,
    ObjectiveGradientMismatch,
    HessianVectorMismatch,
    ReportObjectiveMismatch,
    ReportGradientMismatch,
    FinalOptimalityFailure,
    AccountingMismatch,
    TraceOrdinalMismatch,
    TraceRoleMismatch,
    TraceSegmentMismatch,
    NegativeCurvatureMismatch,
    NegativeCurvatureWitnessMismatch,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum MergeRefusal {
    UnsupportedIdentityVersion,
    PayloadIdentityMismatch,
    PayloadSemanticsMismatch(SemanticRefusal),
    ReferenceIdentityMismatch,
}

fn admit_receipt(
    reference: &ReplayIdentity,
    candidate: &RetainedReceipt,
) -> Result<(), MergeRefusal> {
    check_version(candidate.declared_identity.version())
        .map_err(|_| MergeRefusal::UnsupportedIdentityVersion)?;
    if candidate.payload.identity() != candidate.declared_identity {
        return Err(MergeRefusal::PayloadIdentityMismatch);
    }
    validate_semantics(&candidate.payload).map_err(MergeRefusal::PayloadSemanticsMismatch)?;
    if &candidate.declared_identity != reference {
        return Err(MergeRefusal::ReferenceIdentityMismatch);
    }
    Ok(())
}

fn bits(values: &[f64]) -> Vec<u64> {
    values.iter().map(|value| value.to_bits()).collect()
}

fn fixture_rosenbrock(x: &[f64]) -> (f64, Vec<f64>) {
    assert_eq!(x.len(), DIMENSION);
    let mut objective = 0.0;
    let mut gradient = vec![0.0; DIMENSION];
    for index in 0..DIMENSION - 1 {
        let center = 1.0 - x[index];
        let valley = x[index + 1] - x[index] * x[index];
        objective += center.mul_add(center, 100.0 * valley * valley);
        gradient[index] += (-2.0f64).mul_add(center, -400.0 * x[index] * valley);
        gradient[index + 1] += 200.0 * valley;
    }
    (objective, gradient)
}

fn fixture_rosenbrock_hessian_vector(x: &[f64], direction: &[f64]) -> Vec<f64> {
    assert_eq!(x.len(), DIMENSION);
    assert_eq!(direction.len(), DIMENSION);
    let mut product = vec![0.0; DIMENSION];
    for index in 0..DIMENSION - 1 {
        let xi = x[index];
        let valley = x[index + 1] - xi * xi;
        let diagonal = 2.0 - 400.0 * valley + 800.0 * xi * xi;
        let off_diagonal = -400.0 * xi;
        product[index] += diagonal.mul_add(direction[index], off_diagonal * direction[index + 1]);
        product[index + 1] += off_diagonal.mul_add(direction[index], 200.0 * direction[index + 1]);
    }
    product
}

fn decode_finite(bits: &[u64]) -> Result<Vec<f64>, SemanticRefusal> {
    if bits.len() != DIMENSION {
        return Err(SemanticRefusal::DimensionMismatch);
    }
    let values: Vec<f64> = bits.iter().map(|&value| f64::from_bits(value)).collect();
    if !values.iter().all(|value| value.is_finite()) {
        return Err(SemanticRefusal::NonFiniteEvidence);
    }
    Ok(values)
}

fn oracle_rosenbrock(x: &[f64]) -> Result<(f64, Vec<f64>), SemanticRefusal> {
    if x.len() != DIMENSION {
        return Err(SemanticRefusal::DimensionMismatch);
    }
    let mut objective = 0.0;
    let mut gradient = vec![0.0; DIMENSION];
    for index in 0..DIMENSION - 1 {
        let linear_residual = x[index] - 1.0;
        let valley_residual = x[index] * x[index] - x[index + 1];
        objective += linear_residual * linear_residual + 100.0 * valley_residual * valley_residual;
        gradient[index] += 2.0 * linear_residual + 400.0 * x[index] * valley_residual;
        gradient[index + 1] -= 200.0 * valley_residual;
    }
    if !objective.is_finite() || !gradient.iter().all(|value| value.is_finite()) {
        return Err(SemanticRefusal::NonFiniteEvidence);
    }
    Ok((objective, gradient))
}

fn oracle_rosenbrock_hessian_vector(
    x: &[f64],
    direction: &[f64],
) -> Result<Vec<f64>, SemanticRefusal> {
    if x.len() != DIMENSION || direction.len() != DIMENSION {
        return Err(SemanticRefusal::DimensionMismatch);
    }
    let mut diagonal = vec![0.0; DIMENSION];
    let mut upper = vec![0.0; DIMENSION - 1];
    for index in 0..DIMENSION - 1 {
        diagonal[index] += 1200.0 * x[index] * x[index] - 400.0 * x[index + 1] + 2.0;
        diagonal[index + 1] += 200.0;
        upper[index] = -400.0 * x[index];
    }
    let mut product: Vec<f64> = diagonal
        .iter()
        .zip(direction)
        .map(|(entry, component)| entry * component)
        .collect();
    for index in 0..DIMENSION - 1 {
        product[index] += upper[index] * direction[index + 1];
        product[index + 1] += upper[index] * direction[index];
    }
    if !product.iter().all(|value| value.is_finite()) {
        return Err(SemanticRefusal::NonFiniteEvidence);
    }
    Ok(product)
}

fn inf_norm(values: &[f64]) -> f64 {
    values
        .iter()
        .map(|value| value.abs())
        .fold(0.0f64, f64::max)
}

fn dot(left: &[f64], right: &[f64]) -> f64 {
    left.iter().zip(right).map(|(a, b)| a * b).sum()
}

fn approximately_equal(actual: f64, expected: f64) -> bool {
    if !actual.is_finite() || !expected.is_finite() {
        return false;
    }
    let scale = actual.abs().max(expected.abs()).max(1.0);
    (actual - expected).abs() <= SEMANTIC_COMPARISON_EPSILON_MULTIPLIER * f64::EPSILON * scale
}

fn vectors_approximately_equal(actual: &[f64], expected: &[f64]) -> bool {
    actual.len() == expected.len()
        && actual
            .iter()
            .zip(expected)
            .all(|(&actual, &expected)| approximately_equal(actual, expected))
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum TraceEventLocator {
    Objective(usize),
    HessianVector(usize),
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct IterationSegment {
    steihaug_call_indices: Vec<usize>,
    model_decrease_call_index: usize,
}

fn classify_captured_hessian_calls(
    objective_calls: &[ObjectiveCall],
    calls: Vec<CapturedHessianVectorCall>,
    report_iterations: usize,
) -> Result<Vec<HessianVectorCall>, SemanticRefusal> {
    if objective_calls.len() != report_iterations.saturating_add(1) {
        return Err(SemanticRefusal::AccountingMismatch);
    }

    let mut cursor = 0_usize;
    let mut classified = Vec::with_capacity(calls.len());
    for outer_iteration in 0..report_iterations {
        let opening_ordinal = objective_calls[outer_iteration].ordinal;
        let trial_ordinal = objective_calls[outer_iteration + 1].ordinal;
        let segment_start = cursor;
        while cursor < calls.len() && calls[cursor].ordinal < trial_ordinal {
            if calls[cursor].ordinal <= opening_ordinal {
                return Err(SemanticRefusal::TraceOrdinalMismatch);
            }
            cursor += 1;
        }
        if cursor == segment_start {
            return Err(SemanticRefusal::TraceSegmentMismatch);
        }

        let model_index = cursor - 1;
        for (call_index, call) in calls[segment_start..cursor].iter().enumerate() {
            classified.push(HessianVectorCall {
                ordinal: call.ordinal,
                outer_iteration,
                role: if segment_start + call_index == model_index {
                    HessianVectorRole::ModelDecrease
                } else {
                    HessianVectorRole::Steihaug
                },
                point_bits: call.point_bits.clone(),
                direction_bits: call.direction_bits.clone(),
                product_bits: call.product_bits.clone(),
            });
        }
    }
    if cursor != calls.len() {
        return Err(SemanticRefusal::TraceOrdinalMismatch);
    }
    Ok(classified)
}

#[allow(clippy::too_many_lines)] // The parser keeps the whole callback protocol auditable in one state machine.
fn reconstruct_iteration_segments(
    start_bits: &[u64],
    objective_calls: &[ObjectiveCall],
    hessian_vector_calls: &[HessianVectorCall],
    report_iterations: usize,
) -> Result<Vec<IterationSegment>, SemanticRefusal> {
    if objective_calls.len() != report_iterations.saturating_add(1) {
        return Err(SemanticRefusal::AccountingMismatch);
    }
    if objective_calls
        .first()
        .is_none_or(|call| call.point_bits != start_bits)
    {
        return Err(SemanticRefusal::TraceSegmentMismatch);
    }

    let event_count = objective_calls
        .len()
        .checked_add(hessian_vector_calls.len())
        .ok_or(SemanticRefusal::AccountingMismatch)?;
    let mut events = vec![None; event_count];
    for (index, call) in objective_calls.iter().enumerate() {
        let Some(slot) = events.get_mut(call.ordinal) else {
            return Err(SemanticRefusal::TraceOrdinalMismatch);
        };
        if slot.replace(TraceEventLocator::Objective(index)).is_some() {
            return Err(SemanticRefusal::TraceOrdinalMismatch);
        }
    }
    for (index, call) in hessian_vector_calls.iter().enumerate() {
        let Some(slot) = events.get_mut(call.ordinal) else {
            return Err(SemanticRefusal::TraceOrdinalMismatch);
        };
        if slot
            .replace(TraceEventLocator::HessianVector(index))
            .is_some()
        {
            return Err(SemanticRefusal::TraceOrdinalMismatch);
        }
    }
    if events.iter().any(Option::is_none) {
        return Err(SemanticRefusal::TraceOrdinalMismatch);
    }

    let mut cursor = 0_usize;
    let mut previous_current_point: Option<Vec<u64>> = None;
    let mut segments = Vec::with_capacity(report_iterations);
    for outer_iteration in 0..report_iterations {
        if events.get(cursor).copied().flatten()
            != Some(TraceEventLocator::Objective(outer_iteration))
        {
            return Err(SemanticRefusal::TraceOrdinalMismatch);
        }
        cursor += 1;

        let mut segment_calls = Vec::new();
        while cursor < events.len() {
            match events[cursor] {
                Some(TraceEventLocator::HessianVector(call_index)) => {
                    segment_calls.push(call_index);
                    cursor += 1;
                }
                Some(TraceEventLocator::Objective(objective_index)) => {
                    if objective_index != outer_iteration + 1 {
                        return Err(SemanticRefusal::TraceOrdinalMismatch);
                    }
                    break;
                }
                None => return Err(SemanticRefusal::TraceOrdinalMismatch),
            }
        }
        if cursor == events.len() {
            return Err(SemanticRefusal::TraceOrdinalMismatch);
        }
        let Some((&model_decrease_call_index, steihaug_call_indices)) = segment_calls.split_last()
        else {
            return Err(SemanticRefusal::TraceSegmentMismatch);
        };

        for &call_index in steihaug_call_indices {
            let call = &hessian_vector_calls[call_index];
            if call.outer_iteration != outer_iteration {
                return Err(SemanticRefusal::TraceSegmentMismatch);
            }
            if call.role != HessianVectorRole::Steihaug {
                return Err(SemanticRefusal::TraceRoleMismatch);
            }
        }
        let model_call = &hessian_vector_calls[model_decrease_call_index];
        if model_call.outer_iteration != outer_iteration {
            return Err(SemanticRefusal::TraceSegmentMismatch);
        }
        if model_call.role != HessianVectorRole::ModelDecrease {
            return Err(SemanticRefusal::TraceRoleMismatch);
        }
        if segment_calls
            .iter()
            .any(|&call_index| hessian_vector_calls[call_index].point_bits != model_call.point_bits)
        {
            return Err(SemanticRefusal::TraceSegmentMismatch);
        }

        // The opening objective event is the previous iteration's trial. The
        // next current point is therefore either that trial (accepted) or the
        // previous current point (rejected); no third state is admissible.
        if outer_iteration == 0 {
            if model_call.point_bits != objective_calls[0].point_bits {
                return Err(SemanticRefusal::TraceSegmentMismatch);
            }
        } else if previous_current_point.as_deref() != Some(model_call.point_bits.as_slice())
            && objective_calls[outer_iteration].point_bits != model_call.point_bits
        {
            return Err(SemanticRefusal::TraceSegmentMismatch);
        }

        // Production calls Hv(x, p) exactly once after Steihaug, then evaluates
        // fg(x + p). This exact addition is the independent role discriminator
        // that prevents a search-direction Hv from masquerading as the final
        // model-decrease call.
        let current_point = decode_finite(&model_call.point_bits)?;
        let step = decode_finite(&model_call.direction_bits)?;
        let trial_point = decode_finite(&objective_calls[outer_iteration + 1].point_bits)?;
        if current_point
            .iter()
            .zip(&step)
            .zip(&trial_point)
            .any(|((&current, &step), &trial)| (current + step).to_bits() != trial.to_bits())
        {
            return Err(SemanticRefusal::TraceSegmentMismatch);
        }

        previous_current_point = Some(model_call.point_bits.clone());
        segments.push(IterationSegment {
            steihaug_call_indices: steihaug_call_indices.to_vec(),
            model_decrease_call_index,
        });
    }
    if cursor.checked_add(1) != Some(events.len())
        || events.get(cursor).copied().flatten()
            != Some(TraceEventLocator::Objective(report_iterations))
    {
        return Err(SemanticRefusal::TraceOrdinalMismatch);
    }
    Ok(segments)
}

fn independently_count_negative_curvature_triggers(
    calls: &[HessianVectorCall],
    segments: &[IterationSegment],
) -> Result<Vec<usize>, SemanticRefusal> {
    let mut triggers = Vec::new();
    for segment in segments {
        if calls[segment.model_decrease_call_index].role != HessianVectorRole::ModelDecrease {
            return Err(SemanticRefusal::TraceRoleMismatch);
        }
        for (position, &call_index) in segment.steihaug_call_indices.iter().enumerate() {
            let call = &calls[call_index];
            let point = decode_finite(&call.point_bits)?;
            let direction = decode_finite(&call.direction_bits)?;
            let recorded_product = decode_finite(&call.product_bits)?;
            if inf_norm(&direction) == 0.0 {
                return Err(SemanticRefusal::NegativeCurvatureMismatch);
            }
            let product = oracle_rosenbrock_hessian_vector(&point, &direction)?;
            if !vectors_approximately_equal(&recorded_product, &product) {
                return Err(SemanticRefusal::HessianVectorMismatch);
            }
            let recorded_quadratic_form = dot(&direction, &recorded_product);
            let independent_quadratic_form = dot(&direction, &product);
            if !recorded_quadratic_form.is_finite() || !independent_quadratic_form.is_finite() {
                return Err(SemanticRefusal::NonFiniteEvidence);
            }
            let recorded_nonpositive = recorded_quadratic_form <= NEGATIVE_CURVATURE_THRESHOLD;
            let independent_nonpositive =
                independent_quadratic_form <= NEGATIVE_CURVATURE_THRESHOLD;
            if recorded_nonpositive != independent_nonpositive {
                return Err(SemanticRefusal::NegativeCurvatureMismatch);
            }
            if recorded_nonpositive {
                if position + 1 != segment.steihaug_call_indices.len() {
                    return Err(SemanticRefusal::NegativeCurvatureMismatch);
                }
                triggers.push(call_index);
            }
        }
    }
    Ok(triggers)
}

fn derive_negative_curvature_witness(
    calls: &[HessianVectorCall],
    trigger_indices: &[usize],
) -> Result<NegativeCurvatureWitness, SemanticRefusal> {
    let &hessian_call_index = trigger_indices
        .first()
        .ok_or(SemanticRefusal::NegativeCurvatureWitnessMismatch)?;
    let call = &calls[hessian_call_index];
    if call.role != HessianVectorRole::Steihaug {
        return Err(SemanticRefusal::NegativeCurvatureWitnessMismatch);
    }
    let point = decode_finite(&call.point_bits)?;
    let direction = decode_finite(&call.direction_bits)?;
    let recorded_product = decode_finite(&call.product_bits)?;
    let product = oracle_rosenbrock_hessian_vector(&point, &direction)?;
    if !vectors_approximately_equal(&recorded_product, &product) {
        return Err(SemanticRefusal::HessianVectorMismatch);
    }
    let quadratic_form = dot(&direction, &product);
    let recorded_quadratic_form = dot(&direction, &recorded_product);
    if !quadratic_form.is_finite()
        || !recorded_quadratic_form.is_finite()
        || quadratic_form > NEGATIVE_CURVATURE_THRESHOLD
        || recorded_quadratic_form > NEGATIVE_CURVATURE_THRESHOLD
    {
        return Err(SemanticRefusal::NegativeCurvatureWitnessMismatch);
    }
    Ok(NegativeCurvatureWitness {
        hessian_call_index,
        call_ordinal: call.ordinal,
        outer_iteration: call.outer_iteration,
        role: call.role,
        point_bits: call.point_bits.clone(),
        direction_bits: call.direction_bits.clone(),
        independently_recomputed_product_bits: bits(&product),
        recorded_quadratic_form_bits: recorded_quadratic_form.to_bits(),
        quadratic_form_bits: quadratic_form.to_bits(),
    })
}

#[allow(clippy::too_many_lines)] // Admission intentionally walks every evidence layer in gate order.
fn validate_semantics(payload: &ReceiptPayload) -> Result<(), SemanticRefusal> {
    if payload.suite != SUITE
        || payload.mutation_seed != MUTATION_SEED
        || payload.event_seed_sentinel != EVENT_SEED_SENTINEL
        || payload.fs_math_version != fs_math::VERSION
        || payload.fs_obs_version != fs_obs::VERSION
        || payload.objective_oracle_version != OBJECTIVE_ORACLE_VERSION
        || payload.hessian_oracle_version != HESSIAN_ORACLE_VERSION
        || payload.trust_config_version != TRUST_CONFIG_VERSION
    {
        return Err(SemanticRefusal::FixtureMetadataMismatch);
    }
    if payload.config != TrustRegionConfigPayload::canonical() {
        return Err(SemanticRefusal::TrustConfigurationMismatch);
    }
    if payload.start_bits != bits(&START) {
        return Err(SemanticRefusal::FixtureMetadataMismatch);
    }
    if payload.objective_calls.is_empty() || payload.hessian_vector_calls.is_empty() {
        return Err(SemanticRefusal::AccountingMismatch);
    }
    if payload.report_iterations == 0
        || payload.report_iterations > MAX_ITERATIONS
        || payload.report_iterations.checked_add(1) != Some(payload.report_evaluations)
        || payload.report_evaluations != payload.objective_calls.len()
        || payload.report_hessian_vector_evaluations != payload.hessian_vector_calls.len()
    {
        return Err(SemanticRefusal::AccountingMismatch);
    }

    for call in &payload.objective_calls {
        let point = decode_finite(&call.point_bits)?;
        let reported_gradient = decode_finite(&call.gradient_bits)?;
        let reported_objective = f64::from_bits(call.objective_bits);
        if !reported_objective.is_finite() {
            return Err(SemanticRefusal::NonFiniteEvidence);
        }
        let (objective, gradient) = oracle_rosenbrock(&point)?;
        if !approximately_equal(reported_objective, objective) {
            return Err(SemanticRefusal::ObjectiveValueMismatch);
        }
        if !vectors_approximately_equal(&reported_gradient, &gradient) {
            return Err(SemanticRefusal::ObjectiveGradientMismatch);
        }
    }

    for call in &payload.hessian_vector_calls {
        let point = decode_finite(&call.point_bits)?;
        let direction = decode_finite(&call.direction_bits)?;
        let recorded_product = decode_finite(&call.product_bits)?;
        let product = oracle_rosenbrock_hessian_vector(&point, &direction)?;
        if !vectors_approximately_equal(&recorded_product, &product) {
            return Err(SemanticRefusal::HessianVectorMismatch);
        }
    }

    let segments = reconstruct_iteration_segments(
        &payload.start_bits,
        &payload.objective_calls,
        &payload.hessian_vector_calls,
        payload.report_iterations,
    )?;
    let trigger_indices =
        independently_count_negative_curvature_triggers(&payload.hessian_vector_calls, &segments)?;
    if trigger_indices.is_empty() || payload.report_negative_curvature_hits != trigger_indices.len()
    {
        return Err(SemanticRefusal::NegativeCurvatureMismatch);
    }
    let expected_witness =
        derive_negative_curvature_witness(&payload.hessian_vector_calls, &trigger_indices)?;
    if payload.negative_curvature_witness != expected_witness {
        return Err(SemanticRefusal::NegativeCurvatureWitnessMismatch);
    }
    let witness_form = f64::from_bits(payload.negative_curvature_witness.quadratic_form_bits);
    let recorded_witness_form = f64::from_bits(
        payload
            .negative_curvature_witness
            .recorded_quadratic_form_bits,
    );
    if !witness_form.is_finite()
        || !recorded_witness_form.is_finite()
        || witness_form > NEGATIVE_CURVATURE_THRESHOLD
        || recorded_witness_form > NEGATIVE_CURVATURE_THRESHOLD
    {
        return Err(SemanticRefusal::NegativeCurvatureWitnessMismatch);
    }

    let report_x = decode_finite(&payload.report_x_bits)?;
    if !payload
        .objective_calls
        .iter()
        .any(|call| call.point_bits == payload.report_x_bits)
    {
        return Err(SemanticRefusal::ReportObjectiveMismatch);
    }
    let report_objective = f64::from_bits(payload.report_f_bits);
    let report_gradient_norm = f64::from_bits(payload.report_gradient_norm_bits);
    if !report_objective.is_finite()
        || !report_gradient_norm.is_finite()
        || report_objective < 0.0
        || report_gradient_norm < 0.0
    {
        return Err(SemanticRefusal::NonFiniteEvidence);
    }
    let (objective, gradient) = oracle_rosenbrock(&report_x)?;
    let gradient_norm = inf_norm(&gradient);
    if !approximately_equal(report_objective, objective) {
        return Err(SemanticRefusal::ReportObjectiveMismatch);
    }
    if !approximately_equal(report_gradient_norm, gradient_norm) {
        return Err(SemanticRefusal::ReportGradientMismatch);
    }
    if report_objective >= FINAL_OBJECTIVE_LIMIT
        || report_gradient_norm >= GRADIENT_TOLERANCE
        || report_x
            .iter()
            .any(|value| (value - 1.0).abs() >= FINAL_COORDINATE_ERROR_LIMIT)
    {
        return Err(SemanticRefusal::FinalOptimalityFailure);
    }

    Ok(())
}

fn run_once(start: &[f64]) -> RunRecord {
    // Both callbacks append to one stream. Its length is the single global
    // call ordinal, so objective/Hv ordering cannot be reconstructed from two
    // independently ordered vectors after the fact.
    let captured_calls = RefCell::new(Vec::new());
    let report = {
        let mut objective = |x: &[f64]| {
            let (value, gradient) = fixture_rosenbrock(x);
            let mut calls = captured_calls.borrow_mut();
            let ordinal = calls.len();
            calls.push(CapturedCall::Objective(ObjectiveCall {
                ordinal,
                point_bits: bits(x),
                objective_bits: value.to_bits(),
                gradient_bits: bits(&gradient),
            }));
            (value, gradient)
        };
        let mut hessian_vector = |x: &[f64], direction: &[f64]| {
            let product = fixture_rosenbrock_hessian_vector(x, direction);
            let mut calls = captured_calls.borrow_mut();
            let ordinal = calls.len();
            calls.push(CapturedCall::HessianVector(CapturedHessianVectorCall {
                ordinal,
                point_bits: bits(x),
                direction_bits: bits(direction),
                product_bits: bits(&product),
            }));
            product
        };
        trust_region_newton(
            start,
            &mut objective,
            &mut hessian_vector,
            GRADIENT_TOLERANCE,
            MAX_ITERATIONS,
        )
    };
    let mut objective_calls = Vec::new();
    let mut captured_hessian_vector_calls = Vec::new();
    for call in captured_calls.into_inner() {
        match call {
            CapturedCall::Objective(call) => objective_calls.push(call),
            CapturedCall::HessianVector(call) => captured_hessian_vector_calls.push(call),
        }
    }
    let hessian_vector_calls = classify_captured_hessian_calls(
        &objective_calls,
        captured_hessian_vector_calls,
        report.iters,
    )
    .expect("production trust-region callbacks must follow the declared iteration protocol");
    RunRecord {
        report,
        objective_calls,
        hessian_vector_calls,
    }
}

fn receipt(start: &[f64], run: &RunRecord) -> RetainedReceipt {
    let start_bits = bits(start);
    let segments = reconstruct_iteration_segments(
        &start_bits,
        &run.objective_calls,
        &run.hessian_vector_calls,
        run.report.iters,
    )
    .expect("retained study callbacks must reconstruct into exact outer-iteration segments");
    let trigger_indices =
        independently_count_negative_curvature_triggers(&run.hessian_vector_calls, &segments)
            .expect("retained study curvature roles must admit independent classification");
    let negative_curvature_witness =
        derive_negative_curvature_witness(&run.hessian_vector_calls, &trigger_indices)
            .expect("retained study must contain an independent negative-curvature witness");
    RetainedReceipt::new(ReceiptPayload {
        suite: SUITE,
        mutation_seed: MUTATION_SEED,
        event_seed_sentinel: EVENT_SEED_SENTINEL,
        fs_math_version: fs_math::VERSION,
        fs_obs_version: fs_obs::VERSION,
        objective_oracle_version: OBJECTIVE_ORACLE_VERSION,
        hessian_oracle_version: HESSIAN_ORACLE_VERSION,
        trust_config_version: TRUST_CONFIG_VERSION,
        config: TrustRegionConfigPayload::canonical(),
        start_bits,
        objective_calls: run.objective_calls.clone(),
        hessian_vector_calls: run.hessian_vector_calls.clone(),
        report_x_bits: bits(&run.report.x),
        report_f_bits: run.report.f.to_bits(),
        report_gradient_norm_bits: run.report.grad_norm.to_bits(),
        report_iterations: run.report.iters,
        report_evaluations: run.report.evals,
        report_hessian_vector_evaluations: run.report.hv_evals,
        report_negative_curvature_hits: run.report.negative_curvature_hits,
        negative_curvature_witness,
    })
}

fn mutate_returned_decision(receipt: &RetainedReceipt) -> (RetainedReceipt, usize, u64) {
    let mut mutant = receipt.clone();
    let coordinate = (MUTATION_SEED as usize) % mutant.payload.report_x_bits.len();
    // Keep the mutation inside the mantissa while selecting a high enough bit
    // that the independently recomputed objective must observably move.
    let mask = 1_u64 << (40 + ((MUTATION_SEED >> 8) % 12));
    mutant.payload.report_x_bits[coordinate] ^= mask;
    assert!(
        f64::from_bits(mutant.payload.report_x_bits[coordinate]).is_finite(),
        "mantissa-only mutation must remain a finite wire-valid decision"
    );
    mutant.reseal();
    (mutant, coordinate, mask)
}

fn assert_stale_and_resealed_refusal(
    reference: &RetainedReceipt,
    mutant: &RetainedReceipt,
    expected: SemanticRefusal,
) {
    assert_ne!(
        mutant.declared_identity, reference.declared_identity,
        "red mutation must move the canonical receipt identity"
    );
    let mut stale = mutant.clone();
    stale.declared_identity = reference.declared_identity.clone();
    assert_eq!(
        admit_receipt(&reference.declared_identity, &stale),
        Err(MergeRefusal::PayloadIdentityMismatch),
        "stale mutation must fail the payload-versus-declared-identity gate"
    );
    assert_eq!(
        admit_receipt(&reference.declared_identity, mutant),
        Err(MergeRefusal::PayloadSemanticsMismatch(expected)),
        "correctly resealed mutation must fail the typed semantic gate"
    );
}

fn emit_receipt(
    reference: &RetainedReceipt,
    mutant: &RetainedReceipt,
    coordinate: usize,
    mask: u64,
) {
    let witness = &reference.payload.negative_curvature_witness;
    let json = format!(
        "{{\"randomness\":\"none\",\"event_seed_sentinel\":{EVENT_SEED_SENTINEL},\"mutation_seed\":{MUTATION_SEED},\
         \"reference_identity\":\"{}\",\"mutant_identity\":\"{}\",\
         \"mutated_coordinate\":{coordinate},\"mantissa_mask\":\"{mask:#018x}\",\
         \"negative_curvature_hessian_call_index\":{},\"negative_curvature_call_ordinal\":{},\
         \"negative_curvature_outer_iteration\":{},\"negative_curvature_role\":\"{}\",\
         \"negative_curvature_recorded_quadratic_form_bits\":\"{:#018x}\",\
         \"negative_curvature_independent_quadratic_form_bits\":\"{:#018x}\",\
         \"reported_negative_curvature_hits\":{},\
         \"merge_refusal\":\"payload-semantics-mismatch\"}}",
        reference.declared_identity.hex(),
        mutant.declared_identity.hex(),
        witness.hessian_call_index,
        witness.call_ordinal,
        witness.outer_iteration,
        witness.role.identity_tag(),
        witness.recorded_quadratic_form_bits,
        witness.quadratic_form_bits,
        reference.payload.report_negative_curvature_hits,
    );
    let mut emitter = Emitter::new(SUITE, "rosenbrock-exact-hessian-vector");
    let receipt_event = emitter.emit(
        Severity::Info,
        EventKind::Custom {
            name: "trust-region-study-replay-receipt".to_string(),
            json,
        },
        None,
    );
    let receipt_line = receipt_event.to_jsonl();
    fs_obs::validate_line(&receipt_line)
        .expect("trust-region receipt must use the fs-obs wire schema");
    println!("{receipt_line}");

    let verdict = emitter.emit(
        Severity::Info,
        EventKind::ConformanceCase {
            suite: SUITE.to_string(),
            case: "rosenbrock-exact-hessian-vector".to_string(),
            pass: true,
            detail: format!(
                "the deterministic non-random Rosenbrock fixture replayed the global callback order, every outer-iteration Hv role, and every report bit; independent objective, gradient, Hessian-vector, final-optimality, and Steihaug-only negative-curvature oracles passed; witness ordinal {} in outer iteration {} had role {}, recorded/independent dTHd bits {:#018x}/{:#018x}, and production exactly matched {} independently counted negative-curvature triggers; mutation seed {MUTATION_SEED:#018x} flipped coordinate {coordinate} mask {mask:#018x}, produced stable identity {}, and typed stale/resealed gates refused every red family",
                witness.call_ordinal,
                witness.outer_iteration,
                witness.role.identity_tag(),
                witness.recorded_quadratic_form_bits,
                witness.quadratic_form_bits,
                reference.payload.report_negative_curvature_hits,
                mutant.declared_identity.hex(),
            ),
            seed: EVENT_SEED_SENTINEL,
        },
        None,
    );
    fs_obs::lint_failure_record(&verdict)
        .expect("trust-region seeded-failure verdict must be replayable");
    let verdict_line = verdict.to_jsonl();
    fs_obs::validate_line(&verdict_line)
        .expect("trust-region verdict must use the fs-obs wire schema");
    println!("{verdict_line}");
}

#[test]
#[allow(clippy::too_many_lines)] // One test exercises the reference plus every resealed red family.
fn trust_region_study_replays_and_rejects_seeded_red_mutation() {
    let start = START.to_vec();
    let reference_run = run_once(&start);
    let reference = receipt(&start, &reference_run);
    assert_eq!(
        validate_semantics(&reference.payload),
        Ok(()),
        "reference receipt failed its independent semantic oracle"
    );
    admit_receipt(&reference.declared_identity, &reference)
        .expect("the internally consistent reference receipt must admit");
    let replay = receipt(&start, &run_once(&start));
    assert_eq!(
        replay, reference,
        "the complete callback trace and public report must replay exactly"
    );

    let (mutant, coordinate, mask) = mutate_returned_decision(&reference);
    let (mutant_repeat, repeat_coordinate, repeat_mask) = mutate_returned_decision(&reference);
    assert_eq!((coordinate, mask), (repeat_coordinate, repeat_mask));
    assert_eq!(
        mutant, mutant_repeat,
        "the seeded red mutation and evidence identity must be stable"
    );
    assert_ne!(mutant.declared_identity, reference.declared_identity);
    assert_stale_and_resealed_refusal(
        &reference,
        &mutant,
        SemanticRefusal::ReportObjectiveMismatch,
    );
    assert_eq!(
        admit_receipt(&mutant.declared_identity, &reference),
        Err(MergeRefusal::ReferenceIdentityMismatch),
        "a valid payload must still fail against a different reference identity"
    );

    let mut objective_payload = reference.payload.clone();
    let objective = f64::from_bits(objective_payload.objective_calls[0].objective_bits);
    objective_payload.objective_calls[0].objective_bits = (objective + 1.0).to_bits();
    let objective_mutant = RetainedReceipt::new(objective_payload);
    assert_stale_and_resealed_refusal(
        &reference,
        &objective_mutant,
        SemanticRefusal::ObjectiveValueMismatch,
    );

    let mut gradient_payload = reference.payload.clone();
    let gradient = f64::from_bits(gradient_payload.objective_calls[0].gradient_bits[0]);
    gradient_payload.objective_calls[0].gradient_bits[0] = (gradient + 1.0).to_bits();
    let gradient_mutant = RetainedReceipt::new(gradient_payload);
    assert_stale_and_resealed_refusal(
        &reference,
        &gradient_mutant,
        SemanticRefusal::ObjectiveGradientMismatch,
    );

    let mut radius_payload = reference.payload.clone();
    radius_payload.config.initial_radius_bits = 0.5f64.to_bits();
    let radius_mutant = RetainedReceipt::new(radius_payload);
    assert_stale_and_resealed_refusal(
        &reference,
        &radius_mutant,
        SemanticRefusal::TrustConfigurationMismatch,
    );

    let mut accounting_payload = reference.payload.clone();
    accounting_payload.report_hessian_vector_evaluations += 1;
    let accounting_mutant = RetainedReceipt::new(accounting_payload);
    assert_stale_and_resealed_refusal(
        &reference,
        &accounting_mutant,
        SemanticRefusal::AccountingMismatch,
    );

    let mut ordinal_payload = reference.payload.clone();
    ordinal_payload.hessian_vector_calls[0].ordinal = ordinal_payload.objective_calls[0].ordinal;
    let ordinal_mutant = RetainedReceipt::new(ordinal_payload);
    assert_stale_and_resealed_refusal(
        &reference,
        &ordinal_mutant,
        SemanticRefusal::TraceOrdinalMismatch,
    );

    let mut role_payload = reference.payload.clone();
    let model_call = role_payload
        .hessian_vector_calls
        .iter_mut()
        .find(|call| call.role == HessianVectorRole::ModelDecrease)
        .expect("reference trace contains a final model-decrease call");
    model_call.role = HessianVectorRole::Steihaug;
    let role_mutant = RetainedReceipt::new(role_payload);
    assert_stale_and_resealed_refusal(&reference, &role_mutant, SemanticRefusal::TraceRoleMismatch);

    let mut segment_payload = reference.payload.clone();
    let corrupted_iteration = segment_payload.hessian_vector_calls[0]
        .outer_iteration
        .saturating_add(1);
    segment_payload.hessian_vector_calls[0].outer_iteration = corrupted_iteration;
    let segment_mutant = RetainedReceipt::new(segment_payload);
    assert_stale_and_resealed_refusal(
        &reference,
        &segment_mutant,
        SemanticRefusal::TraceSegmentMismatch,
    );

    let mut negative_claim_payload = reference.payload.clone();
    negative_claim_payload.report_negative_curvature_hits = 0;
    let negative_claim_mutant = RetainedReceipt::new(negative_claim_payload);
    assert_stale_and_resealed_refusal(
        &reference,
        &negative_claim_mutant,
        SemanticRefusal::NegativeCurvatureMismatch,
    );

    let mut witness_payload = reference.payload.clone();
    let witness_form = f64::from_bits(
        witness_payload
            .negative_curvature_witness
            .quadratic_form_bits,
    );
    witness_payload
        .negative_curvature_witness
        .quadratic_form_bits = (witness_form.abs() + 1.0).to_bits();
    let witness_mutant = RetainedReceipt::new(witness_payload);
    assert_stale_and_resealed_refusal(
        &reference,
        &witness_mutant,
        SemanticRefusal::NegativeCurvatureWitnessMismatch,
    );

    emit_receipt(&reference, &mutant, coordinate, mask);
}
