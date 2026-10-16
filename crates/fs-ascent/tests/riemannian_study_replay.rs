//! G5 study-scale replay and seeded-failure self-test for Riemannian L-BFGS.
//!
//! The production engine minimizes a seeded Rayleigh-quotient fixture on the
//! sphere. The retained receipt binds the canonical typed study configuration,
//! regenerated fixture, every objective callback, complete public optimizer
//! state and report, a dense-Jacobi witness, and algebraically independent
//! objective, gradient, sphere-norm, and residual evidence. A separately
//! regenerated same-seed run must reproduce the receipt byte for byte. A
//! deterministic red mutation flips one mantissa bit in the returned decision,
//! remains finite and wire-decodable after resealing, and is refused by the
//! semantic gate even if its own identity is presented as the reference.
//!
//! The emitted same-ISA identity is intentionally left for the central runtime
//! proof pass to pin after `--nocapture`; this source does not guess it. This is
//! one dimensionless dense row-major sphere fixture. It does not claim all
//! manifolds, persisted checkpoints, cancellation recovery, cross-ISA bitwise
//! equality, cryptographic authenticity, ledger persistence, or performance.

use core::cell::RefCell;

use fs_ascent::{RiemannianLbfgs, RiemannianReport, StopReason, StopRule};
use fs_obs::ident::{IdentityBuilder, ReplayIdentity, check_version};
use fs_obs::{Emitter, EventKind, Severity};
use fs_opt::Manifold;
use fs_rand::StreamKey;

const SUITE: &str = "fs-ascent/riemannian-study-replay";
// The same known-convergent logical fixture coordinates used by the original
// sphere battery, retained here instead of silently inventing another corpus.
const INPUT_SEED: u64 = 41;
const MUTATION_SEED: u64 = 0x5245_445F_4249_5401;
const RNG_KERNEL: u32 = 0xA5C3;
const MATRIX_TILE: u32 = 20;
const START_TILE: u32 = 21;
const DIMENSION: usize = 12;
const MEMORY: usize = 10;
// Must sit ABOVE the double-precision gradient floor for this fixture, or the
// study can never stop on the gradient criterion.
//
// `StopRule::GradNorm` compares `inf_norm(&g)` (riemann.rs:489). Measured on the
// canonical seed-41 fixture, the run stalls at `inf_norm = 4.29e-8`
// (`l2 = 7.13e-8`, `f = -3.0844`) after 19 iterations and 62 of 3000
// evaluations. That is not budget exhaustion and not a defect: for a smooth
// objective the attainable gradient floor is about `sqrt(eps * |f|)`, here
// `sqrt(2.22e-16 * 3.084) = 2.6e-8`, the same order as the observed stall. Below
// that the line search cannot certify a Wolfe step, and `riemann.rs:517`
// deliberately reports `Stall` rather than "mint an unrequested GradNorm reason".
//
// The original `1e-9` was a code-first guess that had never been executed; it
// demands a relative gradient near `3e-10`, which f64 cannot deliver here, so the
// study stalled every time. `1e-6` clears the measured floor by ~23x while still
// pinning a genuinely converged eigenvector (relative gradient ~3e-7).
const GRADIENT_TOLERANCE: f64 = 1e-6;
// Bounded BELOW by the gradient gate, not independent of it. For a Rayleigh
// quotient the Riemannian gradient is `2(Ax - lambda x)`, so the intrinsic
// residual is about half the gradient norm and can never be driven far below it.
// The original `1e-10` demanded a residual three orders under the achievable
// gradient floor (~4e-8 inf-norm), so it was unsatisfiable even at the stall
// point, let alone at any GradNorm stop. Measured here: `1.30e-7` against a
// `1e-6` gradient gate.
//
// This gate is deliberately the weakest of the four quality checks. The one with
// real teeth is the INDEPENDENT Jacobi eigensolver comparison below: the returned
// objective matches the true minimum eigenvalue to `2.75e-14` (measured), which
// is the expected quadratic eigenvalue accuracy `~residual^2/gap` for a residual
// of `1.3e-7`. A wrong eigenvector would fail that gate by many orders regardless
// of what this residual bound is set to.
const SCALED_INTRINSIC_RESIDUAL_TOLERANCE: f64 = 1e-6;
const JACOBI_SCALED_RESIDUAL_TOLERANCE: f64 = 1e-12;
const OBJECTIVE_ORACLE_ABS_TOLERANCE: f64 = 1e-7;
const MANIFOLD_VIOLATION_TOLERANCE: f64 = 1e-14;
const TANGENT_REPLAY_ULP_FACTOR: f64 = 8.0;
const GRADIENT_L2_REPLAY_ULP_FACTOR: f64 = 64.0;
const CALLBACK_ORACLE_ULP_FACTOR: f64 = 256.0;
const SPHERE_NORM_REPLAY_ULP_FACTOR: f64 = 16.0;
const EVALUATION_BUDGET: usize = 3_000;
const MAX_ITERATIONS: usize = 1_000;
const CONFIG_SCHEMA_VERSION: u64 = 1;
const FIXTURE_ORACLE_VERSION: &str = "streamkey-symmetric-rayleigh-v1";
const CALLBACK_ORACLE_VERSION: &str = "compensated-upper-triangle-quadratic-v1";
const NORM_ORACLE_VERSION: &str = "serial-deterministic-hypot-v1";
const RESIDUAL_ORACLE_VERSION: &str = "quotient-with-explicit-xTx-v2";
const JACOBI_WITNESS_VERSION: &str = "jacobi-minimum-plus-independent-residual-v1";

#[derive(Clone, Debug, PartialEq, Eq)]
struct StudyConfig {
    schema_version: u64,
    ambient_dimension: usize,
    memory: usize,
    gradient_tolerance_bits: u64,
    evaluation_budget: usize,
    maximum_iterations: usize,
    scaled_intrinsic_residual_tolerance_bits: u64,
    jacobi_scaled_residual_tolerance_bits: u64,
    objective_oracle_abs_tolerance_bits: u64,
    manifold_violation_tolerance_bits: u64,
    tangent_replay_ulp_factor_bits: u64,
    gradient_l2_replay_ulp_factor_bits: u64,
    callback_oracle_ulp_factor_bits: u64,
    sphere_norm_replay_ulp_factor_bits: u64,
    rng_kernel: u32,
    matrix_tile: u32,
    start_tile: u32,
}

impl StudyConfig {
    fn canonical() -> Self {
        Self {
            schema_version: CONFIG_SCHEMA_VERSION,
            ambient_dimension: DIMENSION,
            memory: MEMORY,
            gradient_tolerance_bits: GRADIENT_TOLERANCE.to_bits(),
            evaluation_budget: EVALUATION_BUDGET,
            maximum_iterations: MAX_ITERATIONS,
            scaled_intrinsic_residual_tolerance_bits: SCALED_INTRINSIC_RESIDUAL_TOLERANCE.to_bits(),
            jacobi_scaled_residual_tolerance_bits: JACOBI_SCALED_RESIDUAL_TOLERANCE.to_bits(),
            objective_oracle_abs_tolerance_bits: OBJECTIVE_ORACLE_ABS_TOLERANCE.to_bits(),
            manifold_violation_tolerance_bits: MANIFOLD_VIOLATION_TOLERANCE.to_bits(),
            tangent_replay_ulp_factor_bits: TANGENT_REPLAY_ULP_FACTOR.to_bits(),
            gradient_l2_replay_ulp_factor_bits: GRADIENT_L2_REPLAY_ULP_FACTOR.to_bits(),
            callback_oracle_ulp_factor_bits: CALLBACK_ORACLE_ULP_FACTOR.to_bits(),
            sphere_norm_replay_ulp_factor_bits: SPHERE_NORM_REPLAY_ULP_FACTOR.to_bits(),
            rng_kernel: RNG_KERNEL,
            matrix_tile: MATRIX_TILE,
            start_tile: START_TILE,
        }
    }

    fn gradient_tolerance(&self) -> f64 {
        f64::from_bits(self.gradient_tolerance_bits)
    }

    fn scaled_intrinsic_residual_tolerance(&self) -> f64 {
        f64::from_bits(self.scaled_intrinsic_residual_tolerance_bits)
    }

    fn jacobi_scaled_residual_tolerance(&self) -> f64 {
        f64::from_bits(self.jacobi_scaled_residual_tolerance_bits)
    }

    fn objective_oracle_abs_tolerance(&self) -> f64 {
        f64::from_bits(self.objective_oracle_abs_tolerance_bits)
    }

    fn manifold_violation_tolerance(&self) -> f64 {
        f64::from_bits(self.manifold_violation_tolerance_bits)
    }

    fn tangent_replay_ulp_factor(&self) -> f64 {
        f64::from_bits(self.tangent_replay_ulp_factor_bits)
    }

    fn gradient_l2_replay_ulp_factor(&self) -> f64 {
        f64::from_bits(self.gradient_l2_replay_ulp_factor_bits)
    }

    fn callback_oracle_ulp_factor(&self) -> f64 {
        f64::from_bits(self.callback_oracle_ulp_factor_bits)
    }

    fn sphere_norm_replay_ulp_factor(&self) -> f64 {
        f64::from_bits(self.sphere_norm_replay_ulp_factor_bits)
    }

    fn stop_rule(&self) -> StopRule {
        StopRule::Any(vec![
            StopRule::GradNorm(self.gradient_tolerance()),
            StopRule::Budget(self.evaluation_budget),
        ])
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct Fixture {
    matrix_bits: Vec<u64>,
    start_bits: Vec<u64>,
}

impl Fixture {
    fn matrix(&self) -> Vec<f64> {
        self.matrix_bits
            .iter()
            .map(|&value| f64::from_bits(value))
            .collect()
    }

    fn start(&self) -> Vec<f64> {
        self.start_bits
            .iter()
            .map(|&value| f64::from_bits(value))
            .collect()
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct ObjectiveCall {
    ordinal: usize,
    point_bits: Vec<u64>,
    objective_bits: u64,
    ambient_gradient_bits: Vec<u64>,
    independent_sphere_norm_bits: u64,
    independent_sphere_violation_bits: u64,
}

#[derive(Debug)]
struct RunRecord {
    state: RiemannianLbfgs,
    report: RiemannianReport,
    objective_calls: Vec<ObjectiveCall>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct JacobiWitness {
    minimum_eigenvalue_bits: u64,
    minimum_eigenvector_bits: Vec<u64>,
    independent_rayleigh_quotient_bits: u64,
    independent_residual_inf_bits: u64,
    independent_residual_l2_bits: u64,
    independent_scaled_residual_bits: u64,
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct OracleSnapshot {
    callback_count: usize,
    callback_max_sphere_violation_bits: u64,
    final_objective_bits: u64,
    final_ambient_gradient_bits: Vec<u64>,
    final_tangent_gradient_bits: Vec<u64>,
    final_gradient_inf_bits: u64,
    final_gradient_l2_bits: u64,
    final_rayleigh_quotient_bits: u64,
    final_residual_inf_bits: u64,
    final_residual_l2_bits: u64,
    final_scaled_residual_bits: u64,
    jacobi: JacobiWitness,
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct ReceiptPayload {
    input_seed: u64,
    config: StudyConfig,
    matrix_bits: Vec<u64>,
    start_bits: Vec<u64>,
    objective_calls: Vec<ObjectiveCall>,
    final_x_bits: Vec<u64>,
    final_f_bits: u64,
    final_gradient_bits: Vec<u64>,
    history_bits: Vec<u64>,
    state_iterations: usize,
    state_evaluations: usize,
    stop_reason: &'static str,
    report_gradient_norm_bits: u64,
    report_gradient_l2_norm_bits: u64,
    report_f_bits: u64,
    report_iterations: usize,
    report_evaluations: usize,
    report_worst_violation_bits: u64,
    oracle: OracleSnapshot,
}

impl ReceiptPayload {
    #[allow(clippy::too_many_lines)] // One canonical field-order audit for the full receipt.
    fn identity(&self) -> ReplayIdentity {
        let mut builder = IdentityBuilder::new("fs-ascent-riemannian-study-receipt-v3")
            .str("suite", SUITE)
            .str("fs-ascent-version", fs_ascent::VERSION)
            .str("fs-opt-version", fs_opt::VERSION)
            .str("fs-rand-version", fs_rand::VERSION)
            .str("fs-la-version", fs_la::VERSION)
            .str("fs-math-version", fs_math::VERSION)
            .str("fs-obs-version", fs_obs::VERSION)
            .u64(
                "fs-obs-identity-schema-version",
                u64::from(fs_obs::ident::IDENT_SCHEMA_VERSION),
            )
            .str("engine", "RiemannianLbfgs")
            .str("manifold", "Sphere")
            .str("units", "all-fixture-and-solver-quantities-are-dimensionless")
            .str(
                "matrix-layout",
                "dense-row-major-symmetric;index=row*ambient-dimension+column",
            )
            .str(
                "determinism-scope",
                "fixed-StreamKey-coordinates;serial-callback-order;same-ISA-bit-stable",
            )
            .str(
                "no-claim-boundary",
                "single-sphere-fixture-only;no-all-manifold-coverage;no-cross-ISA-bitwise-equality;no-checkpoint-or-cancellation-recovery;no-cryptographic-authenticity;no-ledger-persistence;no-performance-claim",
            )
            .str("stop-rule", "Any(GradNorm,Budget);first-satisfied-child-attribution")
            .str("fixture-oracle-version", FIXTURE_ORACLE_VERSION)
            .str("callback-oracle-version", CALLBACK_ORACLE_VERSION)
            .str("norm-oracle-version", NORM_ORACLE_VERSION)
            .str("residual-oracle-version", RESIDUAL_ORACLE_VERSION)
            .str("jacobi-witness-version", JACOBI_WITNESS_VERSION)
            .u64("config-schema-version", self.config.schema_version)
            .u64("ambient-dimension", self.config.ambient_dimension as u64)
            .u64("memory", self.config.memory as u64)
            .u64(
                "gradient-tolerance-bits",
                self.config.gradient_tolerance_bits,
            )
            .u64(
                "evaluation-budget",
                self.config.evaluation_budget as u64,
            )
            .u64(
                "maximum-iterations",
                self.config.maximum_iterations as u64,
            )
            .u64(
                "scaled-intrinsic-residual-tolerance-bits",
                self.config.scaled_intrinsic_residual_tolerance_bits,
            )
            .u64(
                "jacobi-scaled-residual-tolerance-bits",
                self.config.jacobi_scaled_residual_tolerance_bits,
            )
            .u64(
                "objective-oracle-absolute-tolerance-bits",
                self.config.objective_oracle_abs_tolerance_bits,
            )
            .u64(
                "manifold-violation-tolerance-bits",
                self.config.manifold_violation_tolerance_bits,
            )
            .u64(
                "tangent-replay-ulp-factor-bits",
                self.config.tangent_replay_ulp_factor_bits,
            )
            .u64(
                "gradient-l2-replay-ulp-factor-bits",
                self.config.gradient_l2_replay_ulp_factor_bits,
            )
            .u64(
                "callback-oracle-ulp-factor-bits",
                self.config.callback_oracle_ulp_factor_bits,
            )
            .u64(
                "sphere-norm-replay-ulp-factor-bits",
                self.config.sphere_norm_replay_ulp_factor_bits,
            )
            .u64("input-seed", self.input_seed)
            .u64("mutation-seed", MUTATION_SEED)
            .u64("rng-kernel", u64::from(self.config.rng_kernel))
            .u64("matrix-tile", u64::from(self.config.matrix_tile))
            .u64("start-tile", u64::from(self.config.start_tile))
            .u64("matrix-values", self.matrix_bits.len() as u64);
        for &value_bits in &self.matrix_bits {
            builder = builder.u64("matrix-value-bits", value_bits);
        }
        builder = builder.u64("start-values", self.start_bits.len() as u64);
        for &value_bits in &self.start_bits {
            builder = builder.u64("start-value-bits", value_bits);
        }
        builder = builder.u64("objective-calls", self.objective_calls.len() as u64);
        for call in &self.objective_calls {
            builder = builder
                .u64("objective-call-ordinal", call.ordinal as u64)
                .u64("objective-call-point-values", call.point_bits.len() as u64);
            for &value_bits in &call.point_bits {
                builder = builder.u64("objective-call-point-bits", value_bits);
            }
            builder = builder
                .u64("objective-call-value-bits", call.objective_bits)
                .u64(
                    "objective-call-gradient-values",
                    call.ambient_gradient_bits.len() as u64,
                );
            for &value_bits in &call.ambient_gradient_bits {
                builder = builder.u64("objective-call-gradient-bits", value_bits);
            }
            builder = builder
                .u64(
                    "objective-call-independent-sphere-norm-bits",
                    call.independent_sphere_norm_bits,
                )
                .u64(
                    "objective-call-independent-sphere-violation-bits",
                    call.independent_sphere_violation_bits,
                );
        }
        builder = builder.u64("final-values", self.final_x_bits.len() as u64);
        for &value_bits in &self.final_x_bits {
            builder = builder.u64("final-value-bits", value_bits);
        }
        builder = builder.u64("final-objective-bits", self.final_f_bits).u64(
            "final-gradient-values",
            self.final_gradient_bits.len() as u64,
        );
        for &value_bits in &self.final_gradient_bits {
            builder = builder.u64("final-gradient-value-bits", value_bits);
        }
        builder = builder.u64("history-values", self.history_bits.len() as u64);
        for &value_bits in &self.history_bits {
            builder = builder.u64("history-value-bits", value_bits);
        }
        builder = builder
            .u64("state-iterations", self.state_iterations as u64)
            .u64("state-evaluations", self.state_evaluations as u64)
            .str("stop-reason", self.stop_reason)
            .u64("report-gradient-norm-bits", self.report_gradient_norm_bits)
            .u64(
                "report-gradient-l2-norm-bits",
                self.report_gradient_l2_norm_bits,
            )
            .u64("report-objective-bits", self.report_f_bits)
            .u64("report-iterations", self.report_iterations as u64)
            .u64("report-evaluations", self.report_evaluations as u64)
            .u64(
                "report-worst-violation-bits",
                self.report_worst_violation_bits,
            )
            .u64("oracle-callback-count", self.oracle.callback_count as u64)
            .u64(
                "oracle-callback-max-sphere-violation-bits",
                self.oracle.callback_max_sphere_violation_bits,
            )
            .u64(
                "oracle-final-objective-bits",
                self.oracle.final_objective_bits,
            )
            .u64(
                "oracle-final-ambient-gradient-values",
                self.oracle.final_ambient_gradient_bits.len() as u64,
            );
        for &value_bits in &self.oracle.final_ambient_gradient_bits {
            builder = builder.u64("oracle-final-ambient-gradient-bits", value_bits);
        }
        builder = builder.u64(
            "oracle-final-tangent-gradient-values",
            self.oracle.final_tangent_gradient_bits.len() as u64,
        );
        for &value_bits in &self.oracle.final_tangent_gradient_bits {
            builder = builder.u64("oracle-final-tangent-gradient-bits", value_bits);
        }
        builder = builder
            .u64(
                "oracle-final-gradient-inf-bits",
                self.oracle.final_gradient_inf_bits,
            )
            .u64(
                "oracle-final-gradient-l2-bits",
                self.oracle.final_gradient_l2_bits,
            )
            .u64(
                "oracle-final-rayleigh-quotient-bits",
                self.oracle.final_rayleigh_quotient_bits,
            )
            .u64(
                "oracle-final-residual-inf-bits",
                self.oracle.final_residual_inf_bits,
            )
            .u64(
                "oracle-final-residual-l2-bits",
                self.oracle.final_residual_l2_bits,
            )
            .u64(
                "oracle-final-scaled-residual-bits",
                self.oracle.final_scaled_residual_bits,
            )
            .u64(
                "jacobi-minimum-eigenvalue-bits",
                self.oracle.jacobi.minimum_eigenvalue_bits,
            )
            .u64(
                "jacobi-minimum-eigenvector-values",
                self.oracle.jacobi.minimum_eigenvector_bits.len() as u64,
            );
        for &value_bits in &self.oracle.jacobi.minimum_eigenvector_bits {
            builder = builder.u64("jacobi-minimum-eigenvector-bits", value_bits);
        }
        builder
            .u64(
                "jacobi-independent-rayleigh-quotient-bits",
                self.oracle.jacobi.independent_rayleigh_quotient_bits,
            )
            .u64(
                "jacobi-independent-residual-inf-bits",
                self.oracle.jacobi.independent_residual_inf_bits,
            )
            .u64(
                "jacobi-independent-residual-l2-bits",
                self.oracle.jacobi.independent_residual_l2_bits,
            )
            .u64(
                "jacobi-independent-scaled-residual-bits",
                self.oracle.jacobi.independent_scaled_residual_bits,
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
enum MergeRefusal {
    UnsupportedIdentityVersion,
    PayloadIdentityMismatch,
    SemanticMismatch(SemanticRefusal),
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
    validate_semantics(&candidate.payload).map_err(MergeRefusal::SemanticMismatch)?;
    if &candidate.declared_identity != reference {
        return Err(MergeRefusal::ReferenceIdentityMismatch);
    }
    Ok(())
}

fn bits(values: &[f64]) -> Vec<u64> {
    values.iter().map(|value| value.to_bits()).collect()
}

fn stop_reason_name(reason: &StopReason) -> &'static str {
    match reason {
        StopReason::GradNorm => "gradient-norm",
        StopReason::ObjectiveBelow => "objective-below",
        StopReason::Budget => "budget",
        StopReason::Stall => "stall",
        StopReason::Composite => "composite",
        StopReason::IterationCap => "iteration-cap",
    }
}

fn stream_vector(seed: u64, kernel: u32, tile: u32, length: usize) -> Vec<f64> {
    let mut stream = StreamKey { seed, kernel, tile }.stream();
    (0..length)
        .map(|_| 2.0f64.mul_add(stream.next_f64(), -1.0))
        .collect()
}

fn seeded_fixture(seed: u64, config: &StudyConfig) -> Fixture {
    let matrix_len = config
        .ambient_dimension
        .checked_mul(config.ambient_dimension)
        .expect("fixture matrix dimensions must not overflow");
    let raw = stream_vector(seed, config.rng_kernel, config.matrix_tile, matrix_len);
    let mut matrix = vec![0.0; matrix_len];
    for row in 0..config.ambient_dimension {
        for column in 0..config.ambient_dimension {
            matrix[row * config.ambient_dimension + column] = raw
                [row * config.ambient_dimension + column]
                + raw[column * config.ambient_dimension + row];
        }
        matrix[row * config.ambient_dimension + row] += 2.0;
    }

    let mut start = stream_vector(
        seed,
        config.rng_kernel,
        config.start_tile,
        config.ambient_dimension,
    );
    let norm = fs_math::det::sqrt(start.iter().map(|value| value * value).sum());
    assert!(
        norm.is_finite() && norm > 0.0,
        "seeded start must normalize"
    );
    for value in &mut start {
        *value /= norm;
    }
    Fixture {
        matrix_bits: bits(&matrix),
        start_bits: bits(&start),
    }
}

fn engine_rayleigh(matrix: &[f64], x: &[f64], dimension: usize) -> (f64, Vec<f64>) {
    let matrix_len = dimension
        .checked_mul(dimension)
        .expect("objective matrix dimensions must not overflow");
    assert_eq!(matrix.len(), matrix_len);
    assert_eq!(x.len(), dimension);
    let mut ax = vec![0.0; dimension];
    for row in 0..dimension {
        for column in 0..dimension {
            ax[row] = matrix[row * dimension + column].mul_add(x[column], ax[row]);
        }
    }
    let objective = x
        .iter()
        .zip(&ax)
        .map(|(coordinate, product)| coordinate * product)
        .sum();
    let gradient = ax.into_iter().map(|value| 2.0 * value).collect();
    (objective, gradient)
}

#[derive(Clone, Copy, Debug, Default)]
struct CompensatedSum {
    sum: f64,
    correction: f64,
}

impl CompensatedSum {
    fn add(&mut self, value: f64) {
        let next = self.sum + value;
        if self.sum.abs() >= value.abs() {
            self.correction += (self.sum - next) + value;
        } else {
            self.correction += (value - next) + self.sum;
        }
        self.sum = next;
    }

    fn total(self) -> f64 {
        self.sum + self.correction
    }
}

fn checked_stable_l2(values: &[f64]) -> Option<f64> {
    if values.is_empty() || !values.iter().all(|value| value.is_finite()) {
        return None;
    }
    let norm = values.iter().fold(0.0, |accumulator, value| {
        fs_math::det::hypot(accumulator, *value)
    });
    norm.is_finite().then_some(norm)
}

fn checked_inf_norm(values: &[f64]) -> Option<f64> {
    if values.is_empty() || !values.iter().all(|value| value.is_finite()) {
        return None;
    }
    Some(values.iter().map(|value| value.abs()).fold(0.0, f64::max))
}

#[derive(Clone, Debug)]
struct IndependentQuadratic {
    objective: f64,
    action: Vec<f64>,
    ambient_gradient: Vec<f64>,
}

fn independent_quadratic(matrix: &[f64], x: &[f64], dimension: usize) -> IndependentQuadratic {
    let mut objective = CompensatedSum::default();
    let mut action = vec![CompensatedSum::default(); dimension];
    for row in 0..dimension {
        let diagonal = matrix[row * dimension + row];
        objective.add(diagonal * x[row] * x[row]);
        action[row].add(diagonal * x[row]);
        for column in row + 1..dimension {
            let coefficient = matrix[row * dimension + column];
            objective.add(2.0 * coefficient * x[row] * x[column]);
            action[row].add(coefficient * x[column]);
            action[column].add(coefficient * x[row]);
        }
    }
    let action: Vec<f64> = action.into_iter().map(CompensatedSum::total).collect();
    let ambient_gradient = action.iter().map(|value| 2.0 * value).collect();
    IndependentQuadratic {
        objective: objective.total(),
        action,
        ambient_gradient,
    }
}

#[derive(Clone, Copy, Debug)]
struct ResidualEvidence {
    quotient: f64,
    inf: f64,
    l2: f64,
    scaled: f64,
}

fn independent_residual(quadratic: &IndependentQuadratic, x: &[f64]) -> ResidualEvidence {
    let mut denominator = CompensatedSum::default();
    for coordinate in x {
        denominator.add(coordinate * coordinate);
    }
    let denominator = denominator.total();
    assert!(denominator.is_finite() && denominator > 0.0);
    // The quotient is x^T A x / x^T x. Omitting the denominator silently
    // assumes an exactly unit-norm floating-point point and is not acceptable
    // evidence even though the manifold keeps it close to one.
    let quotient = quadratic.objective / denominator;
    let residual: Vec<f64> = quadratic
        .action
        .iter()
        .zip(x)
        .map(|(product, coordinate)| (-quotient).mul_add(*coordinate, *product))
        .collect();
    let inf = checked_inf_norm(&residual).expect("oracle residual must be finite and nonempty");
    let l2 = checked_stable_l2(&residual).expect("oracle residual must be finite and nonempty");
    let action_l2 =
        checked_stable_l2(&quadratic.action).expect("oracle action must be finite and nonempty");
    let x_l2 = checked_stable_l2(x).expect("oracle point must be finite and nonempty");
    let scale = action_l2.max(quotient.abs() * x_l2).max(f64::MIN_POSITIVE);
    ResidualEvidence {
        quotient,
        inf,
        l2,
        scaled: l2 / scale,
    }
}

fn project_sphere_gradient(ambient_gradient: &[f64], x: &[f64]) -> Vec<f64> {
    let mut normal_component = CompensatedSum::default();
    for (gradient, coordinate) in ambient_gradient.iter().zip(x) {
        normal_component.add(gradient * coordinate);
    }
    let normal_component = normal_component.total();
    ambient_gradient
        .iter()
        .zip(x)
        .map(|(gradient, coordinate)| (-normal_component).mul_add(*coordinate, *gradient))
        .collect()
}

fn jacobi_witness(matrix: &[f64], config: &StudyConfig) -> JacobiWitness {
    let (eigenvalues, eigenvectors) = fs_la::eigen::jacobi_eigh(matrix, config.ambient_dimension);
    let minimum_eigenvalue = *eigenvalues
        .first()
        .expect("canonical Jacobi fixture must have an eigenvalue");
    let minimum_eigenvector: Vec<f64> = (0..config.ambient_dimension)
        .map(|row| eigenvectors[row * config.ambient_dimension])
        .collect();
    let quadratic = independent_quadratic(matrix, &minimum_eigenvector, config.ambient_dimension);
    let residual = independent_residual(&quadratic, &minimum_eigenvector);
    JacobiWitness {
        minimum_eigenvalue_bits: minimum_eigenvalue.to_bits(),
        minimum_eigenvector_bits: bits(&minimum_eigenvector),
        independent_rayleigh_quotient_bits: residual.quotient.to_bits(),
        independent_residual_inf_bits: residual.inf.to_bits(),
        independent_residual_l2_bits: residual.l2.to_bits(),
        independent_scaled_residual_bits: residual.scaled.to_bits(),
    }
}

fn oracle_snapshot(
    matrix: &[f64],
    final_x: &[f64],
    calls: &[ObjectiveCall],
    config: &StudyConfig,
) -> OracleSnapshot {
    let quadratic = independent_quadratic(matrix, final_x, config.ambient_dimension);
    let tangent_gradient = project_sphere_gradient(&quadratic.ambient_gradient, final_x);
    let gradient_inf =
        checked_inf_norm(&tangent_gradient).expect("oracle gradient must be finite and nonempty");
    let gradient_l2 =
        checked_stable_l2(&tangent_gradient).expect("oracle gradient must be finite and nonempty");
    let residual = independent_residual(&quadratic, final_x);
    let callback_max_sphere_violation = calls
        .iter()
        .map(|call| f64::from_bits(call.independent_sphere_violation_bits))
        .fold(0.0, f64::max);
    OracleSnapshot {
        callback_count: calls.len(),
        callback_max_sphere_violation_bits: callback_max_sphere_violation.to_bits(),
        final_objective_bits: quadratic.objective.to_bits(),
        final_ambient_gradient_bits: bits(&quadratic.ambient_gradient),
        final_tangent_gradient_bits: bits(&tangent_gradient),
        final_gradient_inf_bits: gradient_inf.to_bits(),
        final_gradient_l2_bits: gradient_l2.to_bits(),
        final_rayleigh_quotient_bits: residual.quotient.to_bits(),
        final_residual_inf_bits: residual.inf.to_bits(),
        final_residual_l2_bits: residual.l2.to_bits(),
        final_scaled_residual_bits: residual.scaled.to_bits(),
        jacobi: jacobi_witness(matrix, config),
    }
}

fn run_once(config: &StudyConfig, fixture: &Fixture) -> RunRecord {
    let matrix = fixture.matrix();
    let start = fixture.start();
    let calls = RefCell::new(Vec::new());
    let mut objective = |x: &[f64]| {
        let (value, ambient_gradient) = engine_rayleigh(&matrix, x, config.ambient_dimension);
        let sphere_norm = checked_stable_l2(x)
            .expect("every objective point must be a finite nonempty sphere point");
        let mut calls = calls.borrow_mut();
        let ordinal = calls.len();
        calls.push(ObjectiveCall {
            ordinal,
            point_bits: bits(x),
            objective_bits: value.to_bits(),
            ambient_gradient_bits: bits(&ambient_gradient),
            independent_sphere_norm_bits: sphere_norm.to_bits(),
            independent_sphere_violation_bits: (sphere_norm - 1.0).abs().to_bits(),
        });
        (value, ambient_gradient)
    };
    let ambient = u32::try_from(config.ambient_dimension)
        .expect("fixture dimension must fit the manifold API");
    let mut state = RiemannianLbfgs::new(
        Manifold::Sphere { ambient },
        &start,
        config.memory,
        &mut objective,
    );
    let stop = config.stop_rule();
    let report = state.run(&mut objective, &stop, config.maximum_iterations);
    drop(objective);
    RunRecord {
        state,
        report,
        objective_calls: calls.into_inner(),
    }
}

fn receipt(
    input_seed: u64,
    config: &StudyConfig,
    fixture: &Fixture,
    run: &RunRecord,
) -> RetainedReceipt {
    let ambient = u32::try_from(config.ambient_dimension)
        .expect("fixture dimension must fit the manifold API");
    assert!(matches!(
        &run.state.manifold,
        Manifold::Sphere { ambient: actual } if *actual == ambient
    ));
    let matrix = fixture.matrix();
    RetainedReceipt::new(ReceiptPayload {
        input_seed,
        config: config.clone(),
        matrix_bits: fixture.matrix_bits.clone(),
        start_bits: fixture.start_bits.clone(),
        objective_calls: run.objective_calls.clone(),
        final_x_bits: bits(&run.state.x),
        final_f_bits: run.state.f.to_bits(),
        final_gradient_bits: bits(&run.state.g),
        history_bits: bits(&run.state.history),
        state_iterations: run.state.iters,
        state_evaluations: run.state.evals,
        stop_reason: stop_reason_name(&run.report.reason),
        report_gradient_norm_bits: run.report.grad_norm.to_bits(),
        report_gradient_l2_norm_bits: run.report.grad_l2_norm.to_bits(),
        report_f_bits: run.report.f.to_bits(),
        report_iterations: run.report.iters,
        report_evaluations: run.report.evals,
        report_worst_violation_bits: run.report.worst_violation.to_bits(),
        oracle: oracle_snapshot(&matrix, &run.state.x, &run.objective_calls, config),
    })
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum SemanticRefusal {
    CanonicalConfigMismatch,
    InputSeedMismatch,
    FixtureDimensionMismatch,
    FixtureRegenerationMismatch,
    MatrixSymmetryMismatch,
    NonFiniteEvidence,
    CallbackOrdinalMismatch,
    CallbackAccountingMismatch,
    CallbackOracleMismatch,
    CallbackNormMismatch,
    StateDimensionMismatch,
    StateAccountingMismatch,
    StateHistoryMismatch,
    ReportMismatch,
    ReportedNormMismatch,
    ReturnedPointNotEvaluated,
    TerminalCallbackMismatch,
    ReturnedGradientMismatch,
    OracleSnapshotMismatch,
    QualityGateMismatch,
}

fn decode_finite(bits: &[u64], expected: usize) -> Result<Vec<f64>, SemanticRefusal> {
    if bits.len() != expected {
        return Err(SemanticRefusal::StateDimensionMismatch);
    }
    let values: Vec<f64> = bits.iter().map(|&value| f64::from_bits(value)).collect();
    if !values.iter().all(|value| value.is_finite()) {
        return Err(SemanticRefusal::NonFiniteEvidence);
    }
    Ok(values)
}

fn oracle_close(left: f64, right: f64, ulp_factor: f64) -> bool {
    left.is_finite()
        && right.is_finite()
        && (left - right).abs() <= ulp_factor * f64::EPSILON * left.abs().max(right.abs()).max(1.0)
}

fn relative_norm_close(left: f64, right: f64, ulp_factor: f64) -> bool {
    left.is_finite()
        && right.is_finite()
        && !left.is_sign_negative()
        && !right.is_sign_negative()
        && (left - right).abs()
            <= ulp_factor * f64::EPSILON * left.abs().max(right.abs()).max(f64::MIN_POSITIVE)
}

#[allow(clippy::too_many_lines)] // Each receipt field has an explicit semantic obligation.
fn validate_semantics(payload: &ReceiptPayload) -> Result<(), SemanticRefusal> {
    if payload.config != StudyConfig::canonical() {
        return Err(SemanticRefusal::CanonicalConfigMismatch);
    }
    if payload.input_seed != INPUT_SEED {
        return Err(SemanticRefusal::InputSeedMismatch);
    }
    let config = &payload.config;
    let dimension = config.ambient_dimension;
    let matrix_len = dimension
        .checked_mul(dimension)
        .ok_or(SemanticRefusal::FixtureDimensionMismatch)?;
    if payload.matrix_bits.len() != matrix_len || payload.start_bits.len() != dimension {
        return Err(SemanticRefusal::FixtureDimensionMismatch);
    }
    let regenerated = seeded_fixture(payload.input_seed, config);
    if payload.matrix_bits != regenerated.matrix_bits
        || payload.start_bits != regenerated.start_bits
    {
        return Err(SemanticRefusal::FixtureRegenerationMismatch);
    }
    let matrix = decode_finite(&payload.matrix_bits, matrix_len)?;
    let start = decode_finite(&payload.start_bits, dimension)?;
    for row in 0..dimension {
        for column in row + 1..dimension {
            if payload.matrix_bits[row * dimension + column]
                != payload.matrix_bits[column * dimension + row]
            {
                return Err(SemanticRefusal::MatrixSymmetryMismatch);
            }
        }
    }
    let start_norm = checked_stable_l2(&start).ok_or(SemanticRefusal::NonFiniteEvidence)?;
    if (start_norm - 1.0).abs() > config.manifold_violation_tolerance() {
        return Err(SemanticRefusal::QualityGateMismatch);
    }

    if payload.objective_calls.is_empty()
        || payload.objective_calls.len() != payload.state_evaluations
        || payload.objective_calls.len() != payload.report_evaluations
        || payload.oracle.callback_count != payload.objective_calls.len()
    {
        return Err(SemanticRefusal::CallbackAccountingMismatch);
    }
    let mut callback_max_sphere_violation = 0.0f64;
    for (expected_ordinal, call) in payload.objective_calls.iter().enumerate() {
        if call.ordinal != expected_ordinal {
            return Err(SemanticRefusal::CallbackOrdinalMismatch);
        }
        let point = decode_finite(&call.point_bits, dimension)?;
        let recorded_gradient = decode_finite(&call.ambient_gradient_bits, dimension)?;
        let recorded_objective = f64::from_bits(call.objective_bits);
        if !recorded_objective.is_finite() {
            return Err(SemanticRefusal::NonFiniteEvidence);
        }
        let independent_norm =
            checked_stable_l2(&point).ok_or(SemanticRefusal::NonFiniteEvidence)?;
        let independent_violation = (independent_norm - 1.0).abs();
        if call.independent_sphere_norm_bits != independent_norm.to_bits()
            || call.independent_sphere_violation_bits != independent_violation.to_bits()
        {
            return Err(SemanticRefusal::CallbackNormMismatch);
        }
        callback_max_sphere_violation = callback_max_sphere_violation.max(independent_violation);
        let oracle = independent_quadratic(&matrix, &point, dimension);
        if !oracle_close(
            recorded_objective,
            oracle.objective,
            config.callback_oracle_ulp_factor(),
        ) || recorded_gradient.iter().zip(&oracle.ambient_gradient).any(
            |(recorded, expected)| {
                !oracle_close(*recorded, *expected, config.callback_oracle_ulp_factor())
            },
        ) {
            return Err(SemanticRefusal::CallbackOracleMismatch);
        }
    }
    // Reject off-manifold callbacks before constructing residual evidence. In
    // particular, `independent_residual` deliberately requires a nonzero point;
    // an untrusted all-zero terminal callback must become a typed refusal here,
    // not an assertion panic during admission.
    if callback_max_sphere_violation >= config.manifold_violation_tolerance() {
        return Err(SemanticRefusal::QualityGateMismatch);
    }

    let final_x = decode_finite(&payload.final_x_bits, dimension)?;
    let final_gradient = decode_finite(&payload.final_gradient_bits, dimension)?;
    let history = decode_finite(
        &payload.history_bits,
        payload
            .state_iterations
            .checked_add(1)
            .ok_or(SemanticRefusal::StateAccountingMismatch)?,
    )?;
    let final_f = f64::from_bits(payload.final_f_bits);
    let report_f = f64::from_bits(payload.report_f_bits);
    let report_gradient_inf = f64::from_bits(payload.report_gradient_norm_bits);
    let report_gradient_l2 = f64::from_bits(payload.report_gradient_l2_norm_bits);
    let report_worst_violation = f64::from_bits(payload.report_worst_violation_bits);
    if ![
        final_f,
        report_f,
        report_gradient_inf,
        report_gradient_l2,
        report_worst_violation,
    ]
    .iter()
    .all(|value| value.is_finite())
    {
        return Err(SemanticRefusal::NonFiniteEvidence);
    }
    if report_gradient_inf.is_sign_negative()
        || report_gradient_l2.is_sign_negative()
        || report_worst_violation.is_sign_negative()
    {
        return Err(SemanticRefusal::ReportedNormMismatch);
    }
    if payload.state_iterations != payload.report_iterations
        || payload.state_evaluations != payload.report_evaluations
        || payload.state_iterations > config.maximum_iterations
        || payload.state_evaluations > config.evaluation_budget
        || payload.state_evaluations < payload.state_iterations.saturating_add(1)
        || history.last().map(|value| value.to_bits()) != Some(payload.final_f_bits)
    {
        return Err(SemanticRefusal::StateAccountingMismatch);
    }
    if payload.stop_reason != "gradient-norm"
        || payload.report_f_bits != payload.final_f_bits
        || payload.report_evaluations >= config.evaluation_budget
    {
        return Err(SemanticRefusal::ReportMismatch);
    }
    let first_call = payload
        .objective_calls
        .first()
        .ok_or(SemanticRefusal::StateHistoryMismatch)?;
    if first_call.point_bits != payload.start_bits
        || first_call.objective_bits != payload.history_bits[0]
    {
        return Err(SemanticRefusal::StateHistoryMismatch);
    }
    let mut callback_cursor = 1usize;
    let terminal_callback_index = payload.objective_calls.len() - 1;
    for &retained_objective_bits in payload
        .history_bits
        .iter()
        .skip(1)
        .take(payload.history_bits.len().saturating_sub(2))
    {
        if callback_cursor >= terminal_callback_index {
            return Err(SemanticRefusal::StateHistoryMismatch);
        }
        let relative_index = payload.objective_calls[callback_cursor..terminal_callback_index]
            .iter()
            .position(|call| call.objective_bits == retained_objective_bits)
            .ok_or(SemanticRefusal::StateHistoryMismatch)?;
        callback_cursor = callback_cursor
            .checked_add(relative_index)
            .and_then(|index| index.checked_add(1))
            .ok_or(SemanticRefusal::StateHistoryMismatch)?;
    }
    let recomputed_gradient_inf =
        checked_inf_norm(&final_gradient).ok_or(SemanticRefusal::NonFiniteEvidence)?;
    let recomputed_gradient_l2 =
        checked_stable_l2(&final_gradient).ok_or(SemanticRefusal::NonFiniteEvidence)?;
    if payload.report_gradient_norm_bits != recomputed_gradient_inf.to_bits()
        || !relative_norm_close(
            report_gradient_l2,
            recomputed_gradient_l2,
            config.gradient_l2_replay_ulp_factor(),
        )
    {
        return Err(SemanticRefusal::ReportedNormMismatch);
    }

    let final_call = payload
        .objective_calls
        .last()
        .ok_or(SemanticRefusal::ReturnedPointNotEvaluated)?;
    if final_call.point_bits != payload.final_x_bits
        || final_call.objective_bits != payload.final_f_bits
    {
        let appeared_earlier = payload.objective_calls[..terminal_callback_index]
            .iter()
            .any(|call| {
                call.point_bits == payload.final_x_bits
                    && call.objective_bits == payload.final_f_bits
            });
        return Err(if appeared_earlier {
            // A GradNorm terminal state is checked immediately after accepting
            // the final line-search trial, so no later objective callback is
            // possible. Finding the returned state only in the prefix proves a
            // fabricated post-terminal trace suffix.
            SemanticRefusal::TerminalCallbackMismatch
        } else {
            SemanticRefusal::ReturnedPointNotEvaluated
        });
    }
    let final_ambient_gradient = decode_finite(&final_call.ambient_gradient_bits, dimension)?;
    let normal_component: f64 = final_ambient_gradient
        .iter()
        .zip(&final_x)
        .map(|(gradient, coordinate)| gradient * coordinate)
        .sum();
    let callback_tangent_gradient: Vec<f64> = final_ambient_gradient
        .iter()
        .zip(&final_x)
        .map(|(gradient, coordinate)| (-normal_component).mul_add(*coordinate, *gradient))
        .collect();
    // Compare the recomputation to the returned gradient within the same ULP
    // idiom the tangency residual below already uses, NOT bitwise.
    //
    // Production projects through the authoritative fs-opt manifold runtime
    // (fs-opt/src/eval.rs:1563): it rescales by `max_abs`, accumulates the normal
    // component with a COMPENSATED sum, and divides by a compensated `||x||^2`.
    // The recomputation above deliberately uses the plain textbook form
    // `g - <g,x> x` with an uncompensated sum and no norm division. Those are two
    // different arithmetic routes to the same mathematical projection, so
    // demanding BIT equality asserted formula identity rather than tangency, and
    // refused every correct run once the optimizer moved onto the fs-opt runtime.
    //
    // The oracle keeps its teeth: an untangent, misscaled, or wrong-point
    // gradient still fails here and at the residual gate below. The exactness
    // this bead actually requires is run-to-run REPLAY, which is asserted
    // separately by comparing two independent runs' complete receipts.
    let returned_gradient_tolerance =
        config.gradient_l2_replay_ulp_factor() * f64::EPSILON * recomputed_gradient_l2.max(1.0);
    if callback_tangent_gradient.len() != final_gradient.len()
        || callback_tangent_gradient
            .iter()
            .zip(&final_gradient)
            .any(|(recomputed, returned)| {
                (recomputed - returned).abs() > returned_gradient_tolerance
            })
    {
        return Err(SemanticRefusal::ReturnedGradientMismatch);
    }
    let mut tangent_residual = CompensatedSum::default();
    for (coordinate, gradient) in final_x.iter().zip(&final_gradient) {
        tangent_residual.add(coordinate * gradient);
    }
    if tangent_residual.total().abs()
        > config.tangent_replay_ulp_factor() * f64::EPSILON * recomputed_gradient_l2.max(1.0)
    {
        return Err(SemanticRefusal::ReturnedGradientMismatch);
    }

    let expected_oracle = oracle_snapshot(&matrix, &final_x, &payload.objective_calls, config);
    if payload.oracle != expected_oracle {
        return Err(SemanticRefusal::OracleSnapshotMismatch);
    }
    let oracle_callback_max = f64::from_bits(payload.oracle.callback_max_sphere_violation_bits);
    if oracle_callback_max.to_bits() != callback_max_sphere_violation.to_bits()
        || (report_worst_violation - callback_max_sphere_violation).abs()
            > config.sphere_norm_replay_ulp_factor() * f64::EPSILON
    {
        return Err(SemanticRefusal::CallbackNormMismatch);
    }
    let final_scaled_residual = f64::from_bits(payload.oracle.final_scaled_residual_bits);
    let jacobi_scaled_residual =
        f64::from_bits(payload.oracle.jacobi.independent_scaled_residual_bits);
    let final_rayleigh = f64::from_bits(payload.oracle.final_rayleigh_quotient_bits);
    let jacobi_minimum = f64::from_bits(payload.oracle.jacobi.minimum_eigenvalue_bits);
    let jacobi_rayleigh = f64::from_bits(payload.oracle.jacobi.independent_rayleigh_quotient_bits);
    let oracle_gradient_inf = f64::from_bits(payload.oracle.final_gradient_inf_bits);
    if ![
        final_scaled_residual,
        jacobi_scaled_residual,
        final_rayleigh,
        jacobi_minimum,
        jacobi_rayleigh,
        oracle_gradient_inf,
    ]
    .iter()
    .all(|value| value.is_finite())
        || final_scaled_residual < 0.0
        || jacobi_scaled_residual < 0.0
        || final_scaled_residual > config.scaled_intrinsic_residual_tolerance()
        || jacobi_scaled_residual > config.jacobi_scaled_residual_tolerance()
        || !oracle_close(
            jacobi_rayleigh,
            jacobi_minimum,
            config.callback_oracle_ulp_factor(),
        )
        || (final_rayleigh - jacobi_minimum).abs() >= config.objective_oracle_abs_tolerance()
        || (final_f - jacobi_minimum).abs() >= config.objective_oracle_abs_tolerance()
        || report_gradient_inf > config.gradient_tolerance()
        || oracle_gradient_inf
            > config.gradient_tolerance() + config.callback_oracle_ulp_factor() * f64::EPSILON
        || report_worst_violation >= config.manifold_violation_tolerance()
    {
        return Err(SemanticRefusal::QualityGateMismatch);
    }
    Ok(())
}

fn mutate_returned_decision(receipt: &RetainedReceipt) -> (RetainedReceipt, usize, u64) {
    let mut mutant = receipt.clone();
    let coordinate_count = u64::try_from(mutant.payload.final_x_bits.len())
        .expect("decision dimension must fit the mutation seed coordinate space");
    assert!(
        coordinate_count > 0,
        "mutation requires a nonempty decision"
    );
    let coordinate_start = MUTATION_SEED % coordinate_count;
    let mantissa_start = (MUTATION_SEED >> 8) % 52;
    let candidate_count = coordinate_count
        .checked_mul(52)
        .expect("mutation candidate count must fit u64");
    for offset in 0..candidate_count {
        let coordinate_delta = offset / 52;
        let coordinate_u64 = coordinate_start
            .checked_add(coordinate_delta)
            .expect("mutation coordinate scan must not overflow")
            % coordinate_count;
        let coordinate = usize::try_from(coordinate_u64)
            .expect("reduced mutation coordinate must fit usize on every target");
        let mantissa_bit = (mantissa_start + offset % 52) % 52;
        let mask = 1_u64 << mantissa_bit;
        let mut candidate_bits = mutant.payload.final_x_bits.clone();
        candidate_bits[coordinate] ^= mask;
        if !f64::from_bits(candidate_bits[coordinate]).is_finite()
            || mutant
                .payload
                .objective_calls
                .iter()
                .any(|call| call.point_bits == candidate_bits)
        {
            continue;
        }
        mutant.payload.final_x_bits = candidate_bits;
        mutant.reseal();
        return (mutant, coordinate, mask);
    }
    panic!("seeded mantissa scan must find a finite point absent from the callback trace")
}

fn emit_receipt(
    reference: &RetainedReceipt,
    mutant: &RetainedReceipt,
    coordinate: usize,
    mask: u64,
) {
    let json = format!(
        "{{\"input_seed\":{INPUT_SEED},\"mutation_seed\":{MUTATION_SEED},\
         \"reference_identity\":\"{}\",\"mutant_identity\":\"{}\",\
         \"mutated_coordinate\":{coordinate},\"mantissa_mask\":\"{mask:#018x}\",\
         \"merge_refusal\":\"semantic-returned-point-not-evaluated\"}}",
        reference.declared_identity.hex(),
        mutant.declared_identity.hex(),
    );
    let mut emitter = Emitter::new(SUITE, "sphere-rayleigh-full-study");
    let receipt_event = emitter.emit(
        Severity::Info,
        EventKind::Custom {
            name: "riemannian-study-replay-receipt".to_string(),
            json,
        },
        None,
    );
    let receipt_line = receipt_event.to_jsonl();
    fs_obs::validate_line(&receipt_line)
        .expect("Riemannian study receipt must use the fs-obs wire schema");
    println!("{receipt_line}");

    let verdict = emitter.emit(
        Severity::Info,
        EventKind::ConformanceCase {
            suite: SUITE.to_string(),
            case: "sphere-rayleigh-full-study".to_string(),
            pass: true,
            detail: format!(
                "input seed {INPUT_SEED:#018x} regenerated the fixture and replayed every callback/state/report/oracle bit; mutation seed {MUTATION_SEED:#018x} flipped coordinate {coordinate} mask {mask:#018x}, produced stable identity {}, and semantic admission refused the resealed result",
                mutant.declared_identity.hex(),
            ),
            seed: INPUT_SEED,
        },
        None,
    );
    fs_obs::lint_failure_record(&verdict)
        .expect("Riemannian seeded-failure verdict must be replayable");
    let verdict_line = verdict.to_jsonl();
    fs_obs::validate_line(&verdict_line)
        .expect("Riemannian seeded-failure verdict must use the fs-obs wire schema");
    println!("{verdict_line}");
}

#[test]
fn riemannian_sphere_study_replays_and_rejects_seeded_red_mutation() {
    let config = StudyConfig::canonical();
    let reference_fixture = seeded_fixture(INPUT_SEED, &config);
    let reference_run = run_once(&config, &reference_fixture);
    assert_eq!(
        reference_run.report.reason,
        StopReason::GradNorm,
        "the retained study must converge before exhausting its declared budget: report={:?}",
        reference_run.report,
    );
    assert!(
        reference_run.report.grad_norm <= config.gradient_tolerance(),
        "reported componentwise gradient norm exceeds its stop gate: {:?}",
        reference_run.report,
    );
    assert!(
        reference_run.report.evals < config.evaluation_budget,
        "gradient convergence must occur strictly before the hard evaluation cap: {:?}",
        reference_run.report,
    );
    assert_eq!(
        reference_run.objective_calls.len(),
        reference_run.state.evals,
        "the callback trace must account for every production evaluation"
    );

    let reference = receipt(INPUT_SEED, &config, &reference_fixture, &reference_run);
    admit_receipt(&reference.declared_identity, &reference)
        .expect("the independently validated reference receipt must admit");

    // Do not feed the first run's materialized matrix or start into the replay:
    // regenerate both from the logical seed and declared StreamKey coordinates.
    let replay_fixture = seeded_fixture(INPUT_SEED, &config);
    assert_eq!(
        replay_fixture, reference_fixture,
        "same logical fixture seed must regenerate every matrix/start bit"
    );
    let replay_run = run_once(&config, &replay_fixture);
    let replay = receipt(INPUT_SEED, &config, &replay_fixture, &replay_run);
    assert_eq!(
        replay, reference,
        "same logical seed must replay every callback, state, report, and oracle bit"
    );

    let foreign_reference = IdentityBuilder::new("fs-ascent-riemannian-foreign-reference-v1")
        .str("suite", SUITE)
        .finish();
    assert_eq!(
        admit_receipt(&foreign_reference, &reference),
        Err(MergeRefusal::ReferenceIdentityMismatch),
        "a semantically valid receipt must still match the selected reference identity"
    );

    let mut alternate_seed_mutant = reference.clone();
    alternate_seed_mutant.payload.input_seed ^= 1;
    alternate_seed_mutant.reseal();
    assert_eq!(
        admit_receipt(
            &alternate_seed_mutant.declared_identity,
            &alternate_seed_mutant
        ),
        Err(MergeRefusal::SemanticMismatch(
            SemanticRefusal::InputSeedMismatch
        )),
        "this fixed study must not let a resealed alternate seed redefine its own reference corpus"
    );

    // A terminal GradNorm report performs no callback after the accepted final
    // trial. Prove that a fully resealed suffix forgery cannot authorize itself:
    // all callback values remain oracle-valid and all counters/snapshots are
    // updated, so only the causal terminal-callback invariant can reject it.
    let mut post_terminal_mutant = reference.clone();
    let mut trailing_call = post_terminal_mutant
        .payload
        .objective_calls
        .first()
        .expect("reference receipt must retain its initial callback")
        .clone();
    trailing_call.ordinal = post_terminal_mutant.payload.objective_calls.len();
    assert_ne!(
        trailing_call.point_bits, post_terminal_mutant.payload.final_x_bits,
        "the seeded start must differ from the converged terminal point"
    );
    post_terminal_mutant
        .payload
        .objective_calls
        .push(trailing_call);
    post_terminal_mutant.payload.state_evaluations = post_terminal_mutant
        .payload
        .state_evaluations
        .checked_add(1)
        .expect("red callback accounting must not overflow");
    post_terminal_mutant.payload.report_evaluations = post_terminal_mutant
        .payload
        .report_evaluations
        .checked_add(1)
        .expect("red report accounting must not overflow");
    let post_terminal_oracle = oracle_snapshot(
        &reference_fixture.matrix(),
        &reference_run.state.x,
        &post_terminal_mutant.payload.objective_calls,
        &config,
    );
    post_terminal_mutant.payload.oracle = post_terminal_oracle;
    post_terminal_mutant.reseal();
    assert_eq!(
        admit_receipt(
            &post_terminal_mutant.declared_identity,
            &post_terminal_mutant
        ),
        Err(MergeRefusal::SemanticMismatch(
            SemanticRefusal::TerminalCallbackMismatch
        )),
        "an oracle-valid callback appended after the returned GradNorm state must fail closed"
    );

    let mut signed_zero_mutant = reference.clone();
    signed_zero_mutant.payload.report_worst_violation_bits = (-0.0f64).to_bits();
    signed_zero_mutant.reseal();
    assert_eq!(
        admit_receipt(&signed_zero_mutant.declared_identity, &signed_zero_mutant),
        Err(MergeRefusal::SemanticMismatch(
            SemanticRefusal::ReportedNormMismatch
        )),
        "a norm-like certificate must reject noncanonical negative zero even when resealed"
    );

    // Exercise the admission ordering with a finite, algebraically valid but
    // off-manifold zero callback. This must be a typed quality refusal before
    // the nonzero-denominator residual oracle is reached, never a panic.
    let mut zero_callback_mutant = reference.clone();
    let terminal_call = zero_callback_mutant
        .payload
        .objective_calls
        .last_mut()
        .expect("reference receipt must retain a terminal callback");
    terminal_call.point_bits.fill(0.0f64.to_bits());
    terminal_call.objective_bits = 0.0f64.to_bits();
    terminal_call.ambient_gradient_bits.fill(0.0f64.to_bits());
    terminal_call.independent_sphere_norm_bits = 0.0f64.to_bits();
    terminal_call.independent_sphere_violation_bits = 1.0f64.to_bits();
    zero_callback_mutant.reseal();
    assert_eq!(
        admit_receipt(
            &zero_callback_mutant.declared_identity,
            &zero_callback_mutant
        ),
        Err(MergeRefusal::SemanticMismatch(
            SemanticRefusal::QualityGateMismatch
        )),
        "an off-manifold terminal callback must refuse without reaching an assertion-based residual oracle"
    );

    let (mutant, coordinate, mask) = mutate_returned_decision(&reference);
    let (mutant_repeat, repeat_coordinate, repeat_mask) = mutate_returned_decision(&reference);
    assert_eq!((coordinate, mask), (repeat_coordinate, repeat_mask));
    assert_eq!(
        mutant, mutant_repeat,
        "the seeded red mutation and its evidence identity must be stable"
    );
    assert_ne!(
        mutant.declared_identity, reference.declared_identity,
        "one returned-decision bit must move the retained evidence identity"
    );
    let mut stale_identity_mutant = mutant.clone();
    stale_identity_mutant.declared_identity = reference.declared_identity.clone();
    assert_eq!(
        admit_receipt(&reference.declared_identity, &stale_identity_mutant),
        Err(MergeRefusal::PayloadIdentityMismatch),
        "changing result bits without resealing must fail the payload-integrity check"
    );
    assert_eq!(
        admit_receipt(&reference.declared_identity, &mutant),
        Err(MergeRefusal::SemanticMismatch(
            SemanticRefusal::ReturnedPointNotEvaluated
        )),
        "a resealed returned point absent from the callback trace must fail semantics before reference comparison"
    );
    assert_eq!(
        admit_receipt(&mutant.declared_identity, &mutant),
        Err(MergeRefusal::SemanticMismatch(
            SemanticRefusal::ReturnedPointNotEvaluated
        )),
        "a resealed semantic mutant must not authorize itself as its own reference"
    );

    emit_receipt(&reference, &mutant, coordinate, mask);
}
