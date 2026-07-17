//! G5 replay and seeded-failure self-tests for all constrained ASCENT engines.
//!
//! Augmented Lagrangian, log-barrier interior point, and active-set SQP solve
//! the same analytic equality-plus-active-inequality fixture. Per engine, a
//! retained receipt binds the exact ordered objective, constraint, and
//! Jacobian-transpose callback trace together with every public report and KKT
//! field. An independent analytic oracle reconstructs the objective,
//! constraint values and gradients, multipliers, raw KKT components, and a
//! scaled KKT residual at the returned point. Same-input repeats must reproduce
//! the receipt byte for byte. Engine/lane-namespaced deterministic red
//! mutations target the returned decision, each constraint callback value/JT
//! family, multiplier and KKT claims, and the declared solver configuration;
//! their seeds causally choose only the within-lane coordinate and direction.
//! Stale and correctly resealed forms are refused by distinct typed admission
//! errors, and a resealed semantic mutant cannot authorize itself as its own
//! reference.
//!
//! This is one small dense fixture. It does not claim all constrained problems,
//! large-scale sparse behavior, cancellation, checkpointing, cross-ISA
//! equality, external-solver parity, cryptographic authenticity, ledger
//! persistence, or performance.

use core::cell::RefCell;

use fs_ascent::auglag::ConstrainedProblem;
use fs_ascent::{KktResidual, augmented_lagrangian, interior_point, sqp};
use fs_obs::ident::{IdentityBuilder, ReplayIdentity, check_version};
use fs_obs::{Emitter, EventKind, Severity};

const SUITE: &str = "fs-ascent/constrained-study-replay";
// `ConformanceCase` requires a seed field. These optimizers are deterministic
// and accept no seed, so zero is an explicitly null event-schema sentinel, not
// an optimizer input or a replay-causality claim.
const EVENT_SEED_NONE: u64 = 0;
const MUTATION_SEED: u64 = 0x434F_4E53_5452_5244;
const DIMENSION: usize = 2;
const EQUALITY_COUNT: usize = 1;
const INEQUALITY_COUNT: usize = 1;
const START: [f64; DIMENSION] = [0.0, 0.0];
const ANALYTIC_OPTIMUM: [f64; DIMENSION] = [1.2, 0.8];
const ANALYTIC_MULTIPLIERS: [f64; 2] = [0.4, 1.2];
const ANALYTIC_X_TOLERANCE: f64 = 1e-4;
const ANALYTIC_MULTIPLIER_TOLERANCE: f64 = 1e-3;
// The production certificate and this oracle evaluate the same finite affine
// fixture through independent code. The absolute and relative coefficients are
// each sixty-four epsilons (at unit scale their sum is 128 epsilons), admitting
// only a small reassociation envelope. The much looser solver tolerance remains
// a separate convergence gate and cannot mask a materially wrong certificate.
const ORACLE_ABS_TOLERANCE: f64 = 64.0 * f64::EPSILON;
const ORACLE_REL_TOLERANCE: f64 = 64.0 * f64::EPSILON;

const CONFIG_SCHEMA_VERSION: u64 = 1;
const OBJECTIVE_ORACLE_VERSION: u64 = 1;
const CONSTRAINT_ORACLE_VERSION: u64 = 1;
const CALLBACK_ORACLE_VERSION: u64 = 2;
const KKT_ORACLE_VERSION: u64 = 1;
const SCALED_KKT_ORACLE_VERSION: u64 = 1;
const ACCEPTANCE_GATE_VERSION: u64 = 2;
const MUTATION_PROTOCOL_VERSION: u64 = 3;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Engine {
    AugmentedLagrangian,
    InteriorPoint,
    Sqp,
}

impl Engine {
    const ALL: [Self; 3] = [Self::AugmentedLagrangian, Self::InteriorPoint, Self::Sqp];

    const fn name(self) -> &'static str {
        match self {
            Self::AugmentedLagrangian => "augmented-lagrangian",
            Self::InteriorPoint => "interior-point",
            Self::Sqp => "active-set-sqp",
        }
    }

    const fn tag(self) -> u64 {
        match self {
            Self::AugmentedLagrangian => 0x414C,
            Self::InteriorPoint => 0x4950,
            Self::Sqp => 0x5351,
        }
    }

    const fn tolerance(self) -> f64 {
        match self {
            Self::InteriorPoint => 1e-6,
            Self::AugmentedLagrangian | Self::Sqp => 1e-7,
        }
    }

    const fn iteration_cap(self) -> usize {
        match self {
            Self::AugmentedLagrangian => 40,
            Self::InteriorPoint | Self::Sqp => 60,
        }
    }

    const fn schedule(self) -> &'static str {
        match self {
            Self::AugmentedLagrangian => {
                "cold-zero-multipliers;phr-inequality;mu0=10;grow-if-feasibility>0.25*previous;mu-growth=10;mu-cap=1e10;inner-memory=10;inner-grad-gate=0.1*tolerance;inner-lbfgs-cap=300;update-lambda=lambda+mu*ce;update-nu=max(0,nu+mu*ci);certificate-after-update"
            }
            Self::InteriorPoint => {
                "phase1-if-any-ci>=-margin;phase1-margin=1e-9;phase1-beta=30;phase1-memory=8;phase1-gates=objective<-2margin-or-grad<1e-12;phase1-cap=300;mu0=1;mu-factor=0.2;rho0=10;rho-factor=2;inner-memory=10;inner-grad=max(0.1mu,0.1tol);inner-cap=400;update-lambda=lambda+rho*ce;nu=mu/max(-ci,1e-300);certificate-before-mu-rho-update"
            }
            Self::Sqp => {
                "initial-hessian=identity;dense-kkt-lu;powell-damped-bfgs;initial-active-ci>-1e-8;activate-ci>-1e-10;clamp-active-nu=max(multiplier,0);drop-multiplier<-1e-10;converge-if-not-dropped-and-step-inf<tol-and-kkt<tol;merit-weight=10;backtrack-cap=40;alpha-factor=0.5;accept-merit<m0-1e-12"
            }
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct StudyConfig {
    kkt_tolerance_bits: u64,
    iteration_cap: usize,
    analytic_x_tolerance_bits: u64,
    analytic_multiplier_tolerance_bits: u64,
    oracle_abs_tolerance_bits: u64,
    oracle_rel_tolerance_bits: u64,
    config_schema_version: u64,
    objective_oracle_version: u64,
    constraint_oracle_version: u64,
    callback_oracle_version: u64,
    kkt_oracle_version: u64,
    scaled_kkt_oracle_version: u64,
    acceptance_gate_version: u64,
    mutation_protocol_version: u64,
}

impl StudyConfig {
    fn for_engine(engine: Engine) -> Self {
        Self {
            kkt_tolerance_bits: engine.tolerance().to_bits(),
            iteration_cap: engine.iteration_cap(),
            analytic_x_tolerance_bits: ANALYTIC_X_TOLERANCE.to_bits(),
            analytic_multiplier_tolerance_bits: ANALYTIC_MULTIPLIER_TOLERANCE.to_bits(),
            oracle_abs_tolerance_bits: ORACLE_ABS_TOLERANCE.to_bits(),
            oracle_rel_tolerance_bits: ORACLE_REL_TOLERANCE.to_bits(),
            config_schema_version: CONFIG_SCHEMA_VERSION,
            objective_oracle_version: OBJECTIVE_ORACLE_VERSION,
            constraint_oracle_version: CONSTRAINT_ORACLE_VERSION,
            callback_oracle_version: CALLBACK_ORACLE_VERSION,
            kkt_oracle_version: KKT_ORACLE_VERSION,
            scaled_kkt_oracle_version: SCALED_KKT_ORACLE_VERSION,
            acceptance_gate_version: ACCEPTANCE_GATE_VERSION,
            mutation_protocol_version: MUTATION_PROTOCOL_VERSION,
        }
    }

    fn kkt_tolerance(&self) -> f64 {
        f64::from_bits(self.kkt_tolerance_bits)
    }

    fn analytic_x_tolerance(&self) -> f64 {
        f64::from_bits(self.analytic_x_tolerance_bits)
    }

    fn analytic_multiplier_tolerance(&self) -> f64 {
        f64::from_bits(self.analytic_multiplier_tolerance_bits)
    }

    fn oracle_abs_tolerance(&self) -> f64 {
        f64::from_bits(self.oracle_abs_tolerance_bits)
    }

    fn oracle_rel_tolerance(&self) -> f64 {
        f64::from_bits(self.oracle_rel_tolerance_bits)
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct CallbackCall {
    kind: &'static str,
    point_bits: Vec<u64>,
    weight_bits: Vec<u64>,
    scalar_bits: Option<u64>,
    output_bits: Vec<u64>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct ReportPayload {
    x_bits: Vec<u64>,
    f_bits: u64,
    kkt_bits: [u64; 4],
    lambda_bits: Vec<u64>,
    nu_bits: Vec<u64>,
    iterations: usize,
    evaluations: usize,
    converged: bool,
}

impl ReportPayload {
    #[allow(clippy::too_many_arguments)] // Mirrors the closed public report surface exactly.
    fn from_parts(
        x: &[f64],
        f: f64,
        kkt: &KktResidual,
        lambda: &[f64],
        nu: &[f64],
        iterations: usize,
        evaluations: usize,
        converged: bool,
    ) -> Self {
        Self {
            x_bits: bits(x),
            f_bits: f.to_bits(),
            kkt_bits: [
                kkt.stationarity.to_bits(),
                kkt.feasibility.to_bits(),
                kkt.dual_feasibility.to_bits(),
                kkt.complementarity.to_bits(),
            ],
            lambda_bits: bits(lambda),
            nu_bits: bits(nu),
            iterations,
            evaluations,
            converged,
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct OracleSnapshot {
    objective_bits: u64,
    objective_gradient_bits: Vec<u64>,
    equality_value_bits: Vec<u64>,
    equality_gradient_bits: Vec<Vec<u64>>,
    inequality_value_bits: Vec<u64>,
    inequality_gradient_bits: Vec<Vec<u64>>,
    stationarity_vector_bits: Vec<u64>,
    kkt_bits: [u64; 4],
    scaled_kkt_bits: u64,
}

#[derive(Debug)]
struct RunRecord {
    callbacks: Vec<CallbackCall>,
    report: ReportPayload,
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct ReceiptPayload {
    engine: Engine,
    config: StudyConfig,
    start_bits: Vec<u64>,
    callbacks: Vec<CallbackCall>,
    report: ReportPayload,
    oracle: OracleSnapshot,
}

impl ReceiptPayload {
    #[allow(clippy::too_many_lines)] // Every acceptance and oracle field is intentionally bound.
    fn identity(&self) -> ReplayIdentity {
        let mut builder = IdentityBuilder::new("fs-ascent-constrained-study-receipt-v3")
            .str("suite", SUITE)
            .str("fs-ascent-version", fs_ascent::VERSION)
            .str("fs-la-version", fs_la::VERSION)
            .str("fs-math-version", fs_math::VERSION)
            .str("fs-obs-version", fs_obs::VERSION)
            .u64(
                "fs-obs-identity-schema-version",
                u64::from(fs_obs::ident::IDENT_SCHEMA_VERSION),
            )
            .str("engine", self.engine.name())
            .u64("engine-tag", self.engine.tag())
            .str("engine-schedule", self.engine.schedule())
            .str(
                "units",
                "coordinates-objective-constraints-and-multipliers-are-dimensionless",
            )
            .str("constraint-convention", "equalities-zero;inequalities-less-than-or-equal-zero")
            .str(
                "kkt-component-order",
                "stationarity;primal-feasibility;dual-feasibility;complementarity",
            )
            .str("optimizer-randomness", "none-fixed-input-no-algorithm-seed")
            .str(
                "event-seed-semantics",
                "zero-is-null-ConformanceCase-schema-sentinel-not-optimizer-input",
            )
            .u64("event-seed-none", EVENT_SEED_NONE)
            .u64("mutation-base-seed", MUTATION_SEED)
            .str(
                "mutation-seed-semantics",
                "engine-and-lane-namespaced;causally-selects-only-within-lane-coordinate-and-direction",
            )
            .str(
                "mutation-lanes",
                "returned-decision;equality-callback-value;equality-callback-jt;inequality-callback-value;inequality-callback-jt;equality-multiplier;inequality-multiplier;kkt-claim;scaled-kkt-claim;schedule-config",
            )
            .u64("dimension", DIMENSION as u64)
            .u64("equality-count", EQUALITY_COUNT as u64)
            .u64("inequality-count", INEQUALITY_COUNT as u64)
            .f64_bits("objective-center-x", 2.0)
            .f64_bits("objective-center-y", 1.0)
            .f64_bits("objective-x-quadratic-weight", 1.0)
            .f64_bits("objective-y-quadratic-weight", 1.0)
            .f64_bits("equality-x-coefficient", 1.0)
            .f64_bits("equality-y-coefficient", 1.0)
            .f64_bits("equality-rhs", 2.0)
            .f64_bits("inequality-x-coefficient", 1.0)
            .f64_bits("inequality-y-coefficient", 0.0)
            .f64_bits("inequality-upper-bound", 1.2)
            .f64_bits("analytic-optimum-x", ANALYTIC_OPTIMUM[0])
            .f64_bits("analytic-optimum-y", ANALYTIC_OPTIMUM[1])
            .f64_bits("analytic-lambda", ANALYTIC_MULTIPLIERS[0])
            .f64_bits("analytic-nu", ANALYTIC_MULTIPLIERS[1])
            .u64(
                "config-schema-version",
                self.config.config_schema_version,
            )
            .f64_bits("kkt-tolerance", self.config.kkt_tolerance())
            .u64("iteration-cap", self.config.iteration_cap as u64)
            .f64_bits(
                "analytic-x-tolerance",
                self.config.analytic_x_tolerance(),
            )
            .f64_bits(
                "analytic-multiplier-tolerance",
                self.config.analytic_multiplier_tolerance(),
            )
            .f64_bits(
                "oracle-absolute-tolerance",
                self.config.oracle_abs_tolerance(),
            )
            .f64_bits(
                "oracle-relative-tolerance",
                self.config.oracle_rel_tolerance(),
            )
            .u64(
                "objective-oracle-version",
                self.config.objective_oracle_version,
            )
            .u64(
                "constraint-oracle-version",
                self.config.constraint_oracle_version,
            )
            .u64(
                "callback-oracle-version",
                self.config.callback_oracle_version,
            )
            .u64("kkt-oracle-version", self.config.kkt_oracle_version)
            .u64(
                "scaled-kkt-oracle-version",
                self.config.scaled_kkt_oracle_version,
            )
            .u64(
                "acceptance-gate-version",
                self.config.acceptance_gate_version,
            )
            .u64(
                "mutation-protocol-version",
                self.config.mutation_protocol_version,
            )
            .str(
                "oracle-comparison-gate",
                "abs+rel*max(abs(reported),abs(oracle));64-epsilon-coefficients",
            )
            .str(
                "scaled-kkt-definition",
                "max(stationarity/(1+gradf_inf+abs(lambda)*gradce_inf+abs(nu)*gradci_inf),abs(ce)/(1+max(x_inf,2)),max(ci,0)/(1+max(x_inf,1.2)),dual/(1+abs(nu)),complementarity/(1+abs(nu)*max(abs(ci),1)))-v1",
            )
            .str(
                "convergence-gate",
                "reported-converged-must-equal-both-reported-raw-kkt-and-independent-raw-kkt-strict-below-tolerance-gates",
            )
            .str(
                "admission-gate-order",
                "identity-version;payload-identity;retained-reference-identity;independent-semantic-oracle",
            )
            .str(
                "analytic-gate",
                "componentwise-open-absolute-bands-for-x-lambda-nu-plus-strictly-positive-active-nu",
            )
            .str(
                "callback-gate",
                "every-callback-input-output-recomputed;value-kinds-required-at-returned-point;JT-kinds-required-at-returned-point-with-exact-reported-multipliers",
            )
            .str(
                "accounting-gate",
                "iterations-in-1-through-cap;positive-reported-evals;objective-callback-minus-report-offset=AL3,IP2-before-cap-or-2..3-at-cap,SQP1-in-loop-or-3-with-terminal-adjacent-final-f-and-kkt-objectives",
            )
            .u64("start-values", self.start_bits.len() as u64);
        for &value_bits in &self.start_bits {
            builder = builder.u64("start-value-bits", value_bits);
        }

        builder = builder.u64("callback-calls", self.callbacks.len() as u64);
        for call in &self.callbacks {
            builder = builder
                .str("callback-kind", call.kind)
                .u64("callback-point-values", call.point_bits.len() as u64);
            for &value_bits in &call.point_bits {
                builder = builder.u64("callback-point-bits", value_bits);
            }
            builder = builder.u64("callback-weight-values", call.weight_bits.len() as u64);
            for &value_bits in &call.weight_bits {
                builder = builder.u64("callback-weight-bits", value_bits);
            }
            builder = builder.flag("callback-has-scalar", call.scalar_bits.is_some());
            if let Some(scalar_bits) = call.scalar_bits {
                builder = builder.u64("callback-scalar-bits", scalar_bits);
            }
            builder = builder.u64("callback-output-values", call.output_bits.len() as u64);
            for &value_bits in &call.output_bits {
                builder = builder.u64("callback-output-bits", value_bits);
            }
        }

        builder = builder.u64("report-x-values", self.report.x_bits.len() as u64);
        for &value_bits in &self.report.x_bits {
            builder = builder.u64("report-x-bits", value_bits);
        }
        builder = builder
            .u64("report-objective-bits", self.report.f_bits)
            .u64("report-kkt-values", self.report.kkt_bits.len() as u64);
        for &value_bits in &self.report.kkt_bits {
            builder = builder.u64("report-kkt-bits", value_bits);
        }
        builder = builder.u64("report-lambda-values", self.report.lambda_bits.len() as u64);
        for &value_bits in &self.report.lambda_bits {
            builder = builder.u64("report-lambda-bits", value_bits);
        }
        builder = builder.u64("report-nu-values", self.report.nu_bits.len() as u64);
        for &value_bits in &self.report.nu_bits {
            builder = builder.u64("report-nu-bits", value_bits);
        }
        builder = builder
            .u64("report-iterations", self.report.iterations as u64)
            .u64("report-evaluations", self.report.evaluations as u64)
            .flag("report-converged", self.report.converged)
            .u64("oracle-objective-bits", self.oracle.objective_bits)
            .u64(
                "oracle-objective-gradient-values",
                self.oracle.objective_gradient_bits.len() as u64,
            );
        for &value_bits in &self.oracle.objective_gradient_bits {
            builder = builder.u64("oracle-objective-gradient-bits", value_bits);
        }
        builder = builder.u64(
            "oracle-equality-values",
            self.oracle.equality_value_bits.len() as u64,
        );
        for &value_bits in &self.oracle.equality_value_bits {
            builder = builder.u64("oracle-equality-value-bits", value_bits);
        }
        builder = builder.u64(
            "oracle-equality-gradient-rows",
            self.oracle.equality_gradient_bits.len() as u64,
        );
        for row in &self.oracle.equality_gradient_bits {
            builder = builder.u64("oracle-equality-gradient-values", row.len() as u64);
            for &value_bits in row {
                builder = builder.u64("oracle-equality-gradient-bits", value_bits);
            }
        }
        builder = builder.u64(
            "oracle-inequality-values",
            self.oracle.inequality_value_bits.len() as u64,
        );
        for &value_bits in &self.oracle.inequality_value_bits {
            builder = builder.u64("oracle-inequality-value-bits", value_bits);
        }
        builder = builder.u64(
            "oracle-inequality-gradient-rows",
            self.oracle.inequality_gradient_bits.len() as u64,
        );
        for row in &self.oracle.inequality_gradient_bits {
            builder = builder.u64("oracle-inequality-gradient-values", row.len() as u64);
            for &value_bits in row {
                builder = builder.u64("oracle-inequality-gradient-bits", value_bits);
            }
        }
        builder = builder.u64(
            "oracle-stationarity-vector-values",
            self.oracle.stationarity_vector_bits.len() as u64,
        );
        for &value_bits in &self.oracle.stationarity_vector_bits {
            builder = builder.u64("oracle-stationarity-vector-bits", value_bits);
        }
        builder = builder.u64("oracle-kkt-values", self.oracle.kkt_bits.len() as u64);
        for &value_bits in &self.oracle.kkt_bits {
            builder = builder.u64("oracle-kkt-bits", value_bits);
        }
        builder
            .u64("oracle-scaled-kkt-bits", self.oracle.scaled_kkt_bits)
            .str(
                "no-claims",
                "all-constrained-problems;large-scale-sparse;cancellation;checkpointing;cross-ISA;external-solver-parity;cryptographic-authenticity;authenticated-ledger;persistence;performance",
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
    UnsupportedIdentityVersion { declared: u32 },
    PayloadIdentityMismatch { declared: u64, computed: u64 },
    ReferenceIdentityMismatch { expected: u64, found: u64 },
    SemanticOracleMismatch,
}

fn admit_receipt(
    reference: &ReplayIdentity,
    candidate: &RetainedReceipt,
) -> Result<(), MergeRefusal> {
    check_version(candidate.declared_identity.version()).map_err(|_| {
        MergeRefusal::UnsupportedIdentityVersion {
            declared: candidate.declared_identity.version(),
        }
    })?;
    let computed = candidate.payload.identity();
    if computed != candidate.declared_identity {
        return Err(MergeRefusal::PayloadIdentityMismatch {
            declared: candidate.declared_identity.root(),
            computed: computed.root(),
        });
    }
    if &candidate.declared_identity != reference {
        return Err(MergeRefusal::ReferenceIdentityMismatch {
            expected: reference.root(),
            found: candidate.declared_identity.root(),
        });
    }
    if semantic_mismatch(&candidate.payload).is_some() {
        return Err(MergeRefusal::SemanticOracleMismatch);
    }
    Ok(())
}

fn bits(values: &[f64]) -> Vec<u64> {
    values.iter().map(|value| value.to_bits()).collect()
}

fn bit_rows(rows: &[Vec<f64>]) -> Vec<Vec<u64>> {
    rows.iter().map(|row| bits(row)).collect()
}

fn decode_finite(label: &str, value_bits: &[u64], expected_len: usize) -> Result<Vec<f64>, String> {
    if value_bits.len() != expected_len {
        return Err(format!(
            "{label}-length:{}!=expected-{expected_len}",
            value_bits.len()
        ));
    }
    let values: Vec<f64> = value_bits.iter().copied().map(f64::from_bits).collect();
    if values.iter().any(|value| !value.is_finite()) {
        return Err(format!("{label}-non-finite:{value_bits:016x?}"));
    }
    Ok(values)
}

fn objective_oracle(x: &[f64]) -> (f64, Vec<f64>) {
    let dx = x[0] - 2.0;
    let dy = x[1] - 1.0;
    (dx * dx + dy * dy, vec![2.0 * dx, 2.0 * dy])
}

#[allow(clippy::type_complexity)] // Mirrors values plus Jacobian rows for both constraint kinds.
fn constraint_oracle(x: &[f64]) -> (Vec<f64>, Vec<Vec<f64>>, Vec<f64>, Vec<Vec<f64>>) {
    (
        vec![x[0] + x[1] - 2.0],
        vec![vec![1.0, 1.0]],
        vec![x[0] - 1.2],
        vec![vec![1.0, 0.0]],
    )
}

fn inf_norm(values: &[f64]) -> f64 {
    values.iter().map(|value| value.abs()).fold(0.0, f64::max)
}

#[allow(clippy::too_many_lines)] // The independent KKT derivation is intentionally explicit.
fn independent_oracle(report: &ReportPayload) -> Result<OracleSnapshot, String> {
    let x = decode_finite("report-x", &report.x_bits, DIMENSION)?;
    let lambda = decode_finite("report-lambda", &report.lambda_bits, EQUALITY_COUNT)?;
    let nu = decode_finite("report-nu", &report.nu_bits, INEQUALITY_COUNT)?;
    if nu.iter().any(|multiplier| *multiplier < 0.0) {
        return Err(format!(
            "negative-inequality-multiplier:{:016x?}",
            report.nu_bits
        ));
    }

    let (objective, objective_gradient) = objective_oracle(&x);
    let (equality_values, equality_gradients, inequality_values, inequality_gradients) =
        constraint_oracle(&x);
    let stationarity = vec![
        (objective_gradient[0] + lambda[0]) + nu[0],
        objective_gradient[1] + lambda[0],
    ];
    let stationarity_residual = inf_norm(&stationarity);
    let primal_feasibility = inf_norm(&equality_values).max(
        inequality_values
            .iter()
            .map(|value| value.max(0.0))
            .fold(0.0, f64::max),
    );
    let dual_feasibility = nu
        .iter()
        .map(|multiplier| (-multiplier).max(0.0))
        .fold(0.0, f64::max);
    let complementarity = inequality_values
        .iter()
        .zip(&nu)
        .map(|(constraint, multiplier)| (constraint * multiplier).abs())
        .fold(0.0, f64::max);

    // Scale each KKT block by the magnitude of the terms that form it. This
    // remains dimensionless for this fixture and prevents a large objective or
    // multiplier from making a small absolute residual look artificially bad.
    let stationarity_scale = 1.0
        + inf_norm(&objective_gradient)
        + lambda[0].abs() * inf_norm(&equality_gradients[0])
        + nu[0].abs() * inf_norm(&inequality_gradients[0]);
    let decision_scale = 1.0 + inf_norm(&x).max(2.0);
    let inequality_scale = 1.0 + inf_norm(&x).max(1.2);
    let scaled_primal = (equality_values[0].abs() / decision_scale)
        .max(inequality_values[0].max(0.0) / inequality_scale);
    let scaled_dual = dual_feasibility / (1.0 + nu[0].abs());
    let scaled_complementarity =
        complementarity / (1.0 + nu[0].abs() * inequality_values[0].abs().max(1.0));
    let scaled_kkt = (stationarity_residual / stationarity_scale)
        .max(scaled_primal)
        .max(scaled_dual)
        .max(scaled_complementarity);
    let all_oracle_values = [
        objective,
        stationarity_residual,
        primal_feasibility,
        dual_feasibility,
        complementarity,
        scaled_kkt,
    ];
    if all_oracle_values.iter().any(|value| !value.is_finite()) {
        return Err("independent-oracle-produced-non-finite-value".to_string());
    }

    Ok(OracleSnapshot {
        objective_bits: objective.to_bits(),
        objective_gradient_bits: bits(&objective_gradient),
        equality_value_bits: bits(&equality_values),
        equality_gradient_bits: bit_rows(&equality_gradients),
        inequality_value_bits: bits(&inequality_values),
        inequality_gradient_bits: bit_rows(&inequality_gradients),
        stationarity_vector_bits: bits(&stationarity),
        kkt_bits: [
            stationarity_residual.to_bits(),
            primal_feasibility.to_bits(),
            dual_feasibility.to_bits(),
            complementarity.to_bits(),
        ],
        scaled_kkt_bits: scaled_kkt.to_bits(),
    })
}

fn oracle_close(reported: f64, oracle: f64, config: &StudyConfig) -> bool {
    if !reported.is_finite() || !oracle.is_finite() {
        return false;
    }
    let scale = reported.abs().max(oracle.abs());
    let tolerance = config
        .oracle_rel_tolerance()
        .mul_add(scale, config.oracle_abs_tolerance());
    (reported - oracle).abs() <= tolerance
}

fn callback_mismatch(index: usize, call: &CallbackCall) -> Option<String> {
    let point = match decode_finite("callback-point", &call.point_bits, DIMENSION) {
        Ok(point) => point,
        Err(error) => return Some(format!("callback[{index}]-{error}")),
    };
    let (objective, objective_gradient) = objective_oracle(&point);
    let (equality_values, equality_gradients, inequality_values, inequality_gradients) =
        constraint_oracle(&point);
    let mismatch = match call.kind {
        "objective-gradient" => {
            if !call.weight_bits.is_empty()
                || call.scalar_bits != Some(objective.to_bits())
                || call.output_bits != bits(&objective_gradient)
            {
                Some("objective-value-gradient-or-shape")
            } else {
                None
            }
        }
        "equality" => {
            if !call.weight_bits.is_empty()
                || call.scalar_bits.is_some()
                || call.output_bits != bits(&equality_values)
            {
                Some("equality-value-or-shape")
            } else {
                None
            }
        }
        "equality-jt" => match decode_finite(
            "callback-equality-weight",
            &call.weight_bits,
            EQUALITY_COUNT,
        ) {
            Ok(weights)
                if call.scalar_bits.is_none()
                    && call.output_bits
                        == bits(&[
                            weights[0] * equality_gradients[0][0],
                            weights[0] * equality_gradients[0][1],
                        ]) =>
            {
                None
            }
            _ => Some("equality-gradient-action-or-shape"),
        },
        "inequality" => {
            if !call.weight_bits.is_empty()
                || call.scalar_bits.is_some()
                || call.output_bits != bits(&inequality_values)
            {
                Some("inequality-value-or-shape")
            } else {
                None
            }
        }
        "inequality-jt" => match decode_finite(
            "callback-inequality-weight",
            &call.weight_bits,
            INEQUALITY_COUNT,
        ) {
            Ok(weights)
                if call.scalar_bits.is_none()
                    && call.output_bits
                        == bits(&[
                            weights[0] * inequality_gradients[0][0],
                            // The production fixture canonicalizes the structural
                            // zero to +0.0 even for a negative probe weight.
                            0.0,
                        ]) =>
            {
                None
            }
            _ => Some("inequality-gradient-action-or-shape"),
        },
        _ => Some("unknown-callback-kind"),
    };
    mismatch.map(|reason| format!("callback[{index}]-{}-{reason}", call.kind))
}

#[allow(clippy::too_many_lines)] // Every acceptance claim is independently audited here.
fn semantic_mismatch(payload: &ReceiptPayload) -> Option<String> {
    let expected_config = StudyConfig::for_engine(payload.engine);
    if payload.config != expected_config {
        return Some(format!(
            "declared-config-does-not-match-{}-fixture",
            payload.engine.name()
        ));
    }
    if payload.start_bits != bits(&START) {
        return Some("declared-start-does-not-match-fixed-fixture".to_string());
    }
    if payload.callbacks.is_empty() {
        return Some("empty-callback-trace".to_string());
    }
    for (index, call) in payload.callbacks.iter().enumerate() {
        if let Some(mismatch) = callback_mismatch(index, call) {
            return Some(mismatch);
        }
    }

    let oracle = match independent_oracle(&payload.report) {
        Ok(oracle) => oracle,
        Err(error) => return Some(error),
    };
    if payload.oracle != oracle {
        return Some("retained-independent-oracle-snapshot-mismatch".to_string());
    }
    let report_f = f64::from_bits(payload.report.f_bits);
    if !report_f.is_finite() || payload.report.f_bits != oracle.objective_bits {
        return Some(format!(
            "reported-objective=0x{:016x};oracle=0x{:016x}",
            payload.report.f_bits, oracle.objective_bits
        ));
    }
    for (component, (&reported_bits, &oracle_bits)) in payload
        .report
        .kkt_bits
        .iter()
        .zip(&oracle.kkt_bits)
        .enumerate()
    {
        let reported = f64::from_bits(reported_bits);
        let oracle_value = f64::from_bits(oracle_bits);
        if reported < 0.0 || !oracle_close(reported, oracle_value, &payload.config) {
            return Some(format!(
                "reported-kkt[{component}]=0x{reported_bits:016x};oracle=0x{oracle_bits:016x}"
            ));
        }
    }

    let tolerance = payload.config.kkt_tolerance();
    let reported_kkt_converged = payload
        .report
        .kkt_bits
        .iter()
        .map(|bits| f64::from_bits(*bits))
        .all(|residual| residual.is_finite() && residual < tolerance);
    let oracle_converged = oracle
        .kkt_bits
        .iter()
        .map(|bits| f64::from_bits(*bits))
        .all(|residual| residual.is_finite() && residual < tolerance);
    if payload.report.converged != reported_kkt_converged
        || payload.report.converged != oracle_converged
        || !oracle_converged
    {
        return Some(format!(
            "convergence-claim:{};reported-kkt-gate:{reported_kkt_converged};independent-kkt-gate:{oracle_converged}",
            payload.report.converged,
        ));
    }
    let scaled_kkt = f64::from_bits(oracle.scaled_kkt_bits);
    if !scaled_kkt.is_finite() || scaled_kkt >= tolerance {
        return Some(format!(
            "scaled-kkt={scaled_kkt:e}>=tolerance={tolerance:e}"
        ));
    }

    let x = match decode_finite("report-x", &payload.report.x_bits, DIMENSION) {
        Ok(x) => x,
        Err(error) => return Some(error),
    };
    for (coordinate, (&actual, &expected)) in x.iter().zip(&ANALYTIC_OPTIMUM).enumerate() {
        if (actual - expected).abs() >= payload.config.analytic_x_tolerance() {
            return Some(format!(
                "analytic-optimum[{coordinate}]:{actual:e}!~{expected:e}"
            ));
        }
    }
    let lambda = match decode_finite("report-lambda", &payload.report.lambda_bits, EQUALITY_COUNT) {
        Ok(lambda) => lambda,
        Err(error) => return Some(error),
    };
    let nu = match decode_finite("report-nu", &payload.report.nu_bits, INEQUALITY_COUNT) {
        Ok(nu) => nu,
        Err(error) => return Some(error),
    };
    for (label, actual, expected) in [
        ("lambda", lambda[0], ANALYTIC_MULTIPLIERS[0]),
        ("nu", nu[0], ANALYTIC_MULTIPLIERS[1]),
    ] {
        if (actual - expected).abs() >= payload.config.analytic_multiplier_tolerance() {
            return Some(format!("analytic-{label}:{actual:e}!~{expected:e}"));
        }
    }
    if nu[0] <= 0.0 {
        return Some("active-inequality-multiplier-is-not-positive".to_string());
    }
    if !(1..=payload.config.iteration_cap).contains(&payload.report.iterations) {
        return Some(format!(
            "reported-iterations:{} not in 1..={}",
            payload.report.iterations, payload.config.iteration_cap
        ));
    }
    if payload.report.evaluations == 0 {
        return Some("reported-evaluations-is-zero".to_string());
    }
    let objective_callbacks = payload
        .callbacks
        .iter()
        .filter(|call| call.kind == "objective-gradient")
        .count();
    let objective_callback_indices: Vec<usize> = payload
        .callbacks
        .iter()
        .enumerate()
        .filter_map(|(index, call)| (call.kind == "objective-gradient").then_some(index))
        .collect();
    let terminal_adjacent_objective_pair = match objective_callback_indices.as_slice() {
        [.., earlier, later] => {
            *earlier + 1 == *later
                && payload.callbacks[*earlier].point_bits == payload.report.x_bits
                && payload.callbacks[*later].point_bits == payload.report.x_bits
        }
        _ => false,
    };
    let Some(uncounted_objective_callbacks) =
        objective_callbacks.checked_sub(payload.report.evaluations)
    else {
        return Some(format!(
            "objective-callbacks:{objective_callbacks}<reported-evaluations:{}",
            payload.report.evaluations
        ));
    };
    let accounting_matches = match payload.engine {
        // AL validates twice before the warm driver, then performs the final
        // objective read after its counted KKT check.
        Engine::AugmentedLagrangian => uncounted_objective_callbacks == 3,
        // IP has one initial validation; convergence may occur in-loop (+final
        // objective) or at the post-cap final KKT/objective pair.
        Engine::InteriorPoint if payload.report.iterations < payload.config.iteration_cap => {
            uncounted_objective_callbacks == 2
        }
        Engine::InteriorPoint => (2..=3).contains(&uncounted_objective_callbacks),
        // SQP has one initial validation. Its post-loop path performs a final
        // objective read immediately before the KKT objective read, so the
        // terminal adjacent pair distinguishes the otherwise uncounted +2.
        Engine::Sqp if terminal_adjacent_objective_pair => uncounted_objective_callbacks == 3,
        Engine::Sqp => uncounted_objective_callbacks == 1,
    };
    if !accounting_matches {
        return Some(format!(
            "{}-objective-accounting:callbacks-{objective_callbacks};reported-{};uncounted-{uncounted_objective_callbacks}",
            payload.engine.name(),
            payload.report.evaluations
        ));
    }
    for required_kind in ["objective-gradient", "equality", "inequality"] {
        if !payload
            .callbacks
            .iter()
            .any(|call| call.kind == required_kind && call.point_bits == payload.report.x_bits)
        {
            return Some(format!(
                "returned-point-missing-callback-kind:{required_kind}"
            ));
        }
    }
    for (required_kind, expected_weights) in [
        ("equality-jt", payload.report.lambda_bits.as_slice()),
        ("inequality-jt", payload.report.nu_bits.as_slice()),
    ] {
        if !payload.callbacks.iter().any(|call| {
            call.kind == required_kind
                && call.point_bits == payload.report.x_bits
                && call.weight_bits.as_slice() == expected_weights
                && call
                    .weight_bits
                    .iter()
                    .any(|bits| *bits & 0x7fff_ffff_ffff_ffff != 0)
        }) {
            return Some(format!(
                "returned-point-missing-nonzero-gradient-action:{required_kind}"
            ));
        }
    }
    None
}

fn push_call(
    calls: &RefCell<Vec<CallbackCall>>,
    kind: &'static str,
    point: &[f64],
    weights: &[f64],
    scalar: Option<f64>,
    output: &[f64],
) {
    calls.borrow_mut().push(CallbackCall {
        kind,
        point_bits: bits(point),
        weight_bits: bits(weights),
        scalar_bits: scalar.map(f64::to_bits),
        output_bits: bits(output),
    });
}

#[allow(clippy::too_many_lines)] // One shared instrumented fixture across three report types.
fn run_once(engine: Engine, config: &StudyConfig) -> RunRecord {
    let calls = RefCell::new(Vec::new());
    let report = {
        let mut objective = |x: &[f64]| {
            let dx = x[0] - 2.0;
            let dy = x[1] - 1.0;
            let value = dx * dx + dy * dy;
            let gradient = vec![2.0 * dx, 2.0 * dy];
            push_call(&calls, "objective-gradient", x, &[], Some(value), &gradient);
            (value, gradient)
        };
        let equality = |x: &[f64]| {
            let value = vec![x[0] + x[1] - 2.0];
            push_call(&calls, "equality", x, &[], None, &value);
            value
        };
        let equality_jt = |x: &[f64], weights: &[f64]| {
            let value = vec![weights[0], weights[0]];
            push_call(&calls, "equality-jt", x, weights, None, &value);
            value
        };
        let inequality = |x: &[f64]| {
            let value = vec![x[0] - 1.2];
            push_call(&calls, "inequality", x, &[], None, &value);
            value
        };
        let inequality_jt = |x: &[f64], weights: &[f64]| {
            let value = vec![weights[0], 0.0];
            push_call(&calls, "inequality-jt", x, weights, None, &value);
            value
        };
        let mut problem = ConstrainedProblem {
            fg: &mut objective,
            ce: &equality,
            ce_jt: &equality_jt,
            ci: &inequality,
            ci_jt: &inequality_jt,
        };
        match engine {
            Engine::AugmentedLagrangian => {
                let report = augmented_lagrangian(
                    &mut problem,
                    &START,
                    config.kkt_tolerance(),
                    config.iteration_cap,
                );
                ReportPayload::from_parts(
                    &report.x,
                    report.f,
                    &report.kkt,
                    &report.lambda,
                    &report.nu,
                    report.outer_iters,
                    report.evals,
                    report.converged,
                )
            }
            Engine::InteriorPoint => {
                let report = interior_point(
                    &mut problem,
                    &START,
                    config.kkt_tolerance(),
                    config.iteration_cap,
                );
                ReportPayload::from_parts(
                    &report.x,
                    report.f,
                    &report.kkt,
                    &report.lambda,
                    &report.nu,
                    report.outer_iters,
                    report.evals,
                    report.converged,
                )
            }
            Engine::Sqp => {
                let report = sqp(
                    &mut problem,
                    &START,
                    config.kkt_tolerance(),
                    config.iteration_cap,
                );
                ReportPayload::from_parts(
                    &report.x,
                    report.f,
                    &report.kkt,
                    &report.lambda,
                    &report.nu,
                    report.iters,
                    report.evals,
                    report.converged,
                )
            }
        }
    };
    RunRecord {
        callbacks: calls.into_inner(),
        report,
    }
}

fn receipt(engine: Engine, config: &StudyConfig, run: &RunRecord) -> RetainedReceipt {
    let oracle = independent_oracle(&run.report)
        .expect("production report must be well-formed enough for the independent oracle");
    RetainedReceipt::new(ReceiptPayload {
        engine,
        config: config.clone(),
        start_bits: bits(&START),
        callbacks: run.callbacks.clone(),
        report: run.report.clone(),
        oracle,
    })
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum MutationKind {
    ReturnedDecision,
    EqualityValue,
    EqualityGradient,
    InequalityValue,
    InequalityGradient,
    EqualityMultiplier,
    InequalityMultiplier,
    KktClaim,
    ScaledKktClaim,
    ScheduleConfig,
}

impl MutationKind {
    const ALL: [Self; 10] = [
        Self::ReturnedDecision,
        Self::EqualityValue,
        Self::EqualityGradient,
        Self::InequalityValue,
        Self::InequalityGradient,
        Self::EqualityMultiplier,
        Self::InequalityMultiplier,
        Self::KktClaim,
        Self::ScaledKktClaim,
        Self::ScheduleConfig,
    ];

    const fn name(self) -> &'static str {
        match self {
            Self::ReturnedDecision => "returned-decision",
            Self::EqualityValue => "equality-callback-value",
            Self::EqualityGradient => "equality-callback-jt",
            Self::InequalityValue => "inequality-callback-value",
            Self::InequalityGradient => "inequality-callback-jt",
            Self::EqualityMultiplier => "equality-multiplier",
            Self::InequalityMultiplier => "inequality-multiplier",
            Self::KktClaim => "kkt-claim",
            Self::ScaledKktClaim => "scaled-kkt-claim",
            Self::ScheduleConfig => "schedule-config",
        }
    }

    const fn tag(self) -> u64 {
        match self {
            Self::ReturnedDecision => 0x01,
            Self::EqualityValue => 0x02,
            Self::EqualityGradient => 0x03,
            Self::InequalityValue => 0x04,
            Self::InequalityGradient => 0x05,
            Self::EqualityMultiplier => 0x06,
            Self::InequalityMultiplier => 0x07,
            Self::KktClaim => 0x08,
            Self::ScaledKktClaim => 0x09,
            Self::ScheduleConfig => 0x0A,
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct MutationProof {
    kind: MutationKind,
    mutation_seed: u64,
    target: String,
    before: u64,
    after: u64,
    resealed_identity: ReplayIdentity,
    stale_refusal: MergeRefusal,
    resealed_refusal: MergeRefusal,
    self_reference_refusal: MergeRefusal,
    semantic_mismatch: String,
}

fn mutate_f64(slot: &mut u64, mutation_seed: u64, magnitude: f64) -> (u64, u64) {
    let before = *slot;
    let direction: f64 = if mutation_seed & 1 == 0 { 1.0 } else { -1.0 };
    let after_value = direction.mul_add(magnitude, f64::from_bits(before));
    assert!(
        after_value.is_finite(),
        "deterministic semantic mutation must remain wire-finite"
    );
    let after = after_value.to_bits();
    assert_ne!(before, after, "semantic mutation must change its target");
    *slot = after;
    (before, after)
}

fn returned_callback_index(payload: &ReceiptPayload, kind: &str) -> usize {
    let expected_weights = match kind {
        "equality-jt" => Some(payload.report.lambda_bits.as_slice()),
        "inequality-jt" => Some(payload.report.nu_bits.as_slice()),
        _ => None,
    };
    payload
        .callbacks
        .iter()
        .position(|call| {
            call.kind == kind
                && call.point_bits == payload.report.x_bits
                && expected_weights.is_none_or(|weights| {
                    call.weight_bits.as_slice() == weights
                        && call
                            .weight_bits
                            .iter()
                            .any(|bits| *bits & 0x7fff_ffff_ffff_ffff != 0)
                })
        })
        .unwrap_or_else(|| panic!("reference receipt lacks returned-point callback {kind}"))
}

#[allow(clippy::too_many_lines)] // One explicit independent lane per semantic claim family.
fn apply_mutation(
    payload: &mut ReceiptPayload,
    kind: MutationKind,
    mutation_seed: u64,
) -> (String, u64, u64) {
    match kind {
        MutationKind::ReturnedDecision => {
            let coordinate = (mutation_seed as usize) % payload.report.x_bits.len();
            let (before, after) = mutate_f64(
                &mut payload.report.x_bits[coordinate],
                mutation_seed,
                4.0 * payload.config.analytic_x_tolerance(),
            );
            (format!("report.x[{coordinate}]"), before, after)
        }
        MutationKind::EqualityValue => {
            let callback = returned_callback_index(payload, "equality");
            let (before, after) = mutate_f64(
                &mut payload.callbacks[callback].output_bits[0],
                mutation_seed,
                0.125,
            );
            (
                format!("callbacks[{callback}].equality.output[0]"),
                before,
                after,
            )
        }
        MutationKind::EqualityGradient => {
            let callback = returned_callback_index(payload, "equality-jt");
            let coordinate = (mutation_seed as usize) % DIMENSION;
            let (before, after) = mutate_f64(
                &mut payload.callbacks[callback].output_bits[coordinate],
                mutation_seed,
                0.125,
            );
            (
                format!("callbacks[{callback}].equality_jt.output[{coordinate}]"),
                before,
                after,
            )
        }
        MutationKind::InequalityValue => {
            let callback = returned_callback_index(payload, "inequality");
            let (before, after) = mutate_f64(
                &mut payload.callbacks[callback].output_bits[0],
                mutation_seed,
                0.125,
            );
            (
                format!("callbacks[{callback}].inequality.output[0]"),
                before,
                after,
            )
        }
        MutationKind::InequalityGradient => {
            let callback = returned_callback_index(payload, "inequality-jt");
            let coordinate = (mutation_seed as usize) % DIMENSION;
            let (before, after) = mutate_f64(
                &mut payload.callbacks[callback].output_bits[coordinate],
                mutation_seed,
                0.125,
            );
            (
                format!("callbacks[{callback}].inequality_jt.output[{coordinate}]"),
                before,
                after,
            )
        }
        MutationKind::EqualityMultiplier => {
            let (before, after) =
                mutate_f64(&mut payload.report.lambda_bits[0], mutation_seed, 0.05);
            ("report.lambda[0]".to_string(), before, after)
        }
        MutationKind::InequalityMultiplier => {
            let (before, after) = mutate_f64(&mut payload.report.nu_bits[0], mutation_seed, 0.05);
            ("report.nu[0]".to_string(), before, after)
        }
        MutationKind::KktClaim => {
            let component = (mutation_seed as usize) % payload.report.kkt_bits.len();
            let (before, after) = mutate_f64(
                &mut payload.report.kkt_bits[component],
                mutation_seed,
                4.0 * payload.config.kkt_tolerance(),
            );
            (format!("report.kkt[{component}]"), before, after)
        }
        MutationKind::ScaledKktClaim => {
            let (before, after) = mutate_f64(
                &mut payload.oracle.scaled_kkt_bits,
                mutation_seed,
                4.0 * payload.config.kkt_tolerance(),
            );
            ("oracle.scaled_kkt".to_string(), before, after)
        }
        MutationKind::ScheduleConfig => {
            let before = payload.config.iteration_cap as u64;
            payload.config.iteration_cap = if mutation_seed & 1 == 0 {
                payload.config.iteration_cap + 1
            } else {
                payload.config.iteration_cap - 1
            };
            let after = payload.config.iteration_cap as u64;
            ("config.iteration_cap".to_string(), before, after)
        }
    }
}

fn mutation_proof(reference: &RetainedReceipt, kind: MutationKind) -> MutationProof {
    let mutation_seed = MUTATION_SEED ^ reference.payload.engine.tag() ^ kind.tag();
    let mut stale = reference.clone();
    let (target, before, after) = apply_mutation(&mut stale.payload, kind, mutation_seed);
    let stale_refusal = admit_receipt(&reference.declared_identity, &stale)
        .expect_err("unsealed semantic corruption must fail payload identity admission");
    let mut resealed = stale;
    resealed.reseal();
    let resealed_identity = resealed.declared_identity.clone();
    let resealed_refusal = admit_receipt(&reference.declared_identity, &resealed)
        .expect_err("resealed semantic corruption must fail retained-reference admission");
    let semantic_mismatch = semantic_mismatch(&resealed.payload)
        .expect("each red mutation must independently violate a semantic oracle or gate");
    let self_reference_refusal = admit_receipt(&resealed.declared_identity, &resealed)
        .expect_err("a resealed semantic mutant must not authorize itself as its own reference");
    MutationProof {
        kind,
        mutation_seed,
        target,
        before,
        after,
        resealed_identity,
        stale_refusal,
        resealed_refusal,
        self_reference_refusal,
        semantic_mismatch,
    }
}

fn emit_receipt(engine: Engine, reference: &RetainedReceipt, mutations: &[MutationProof]) {
    let mutation_json = mutations
        .iter()
        .map(|mutation| {
            format!(
                concat!(
                    "{{\"kind\":\"{}\",\"seed\":\"0x{:016x}\",\"target\":\"{}\",",
                    "\"before\":\"0x{:016x}\",\"after\":\"0x{:016x}\",",
                    "\"resealed_payload_identity\":\"{}\",",
                    "\"stale_refusal\":\"{:?}\",\"resealed_refusal\":\"{:?}\",",
                    "\"self_reference_refusal\":\"{:?}\",",
                    "\"semantic_mismatch\":\"{}\"}}"
                ),
                mutation.kind.name(),
                mutation.mutation_seed,
                mutation.target,
                mutation.before,
                mutation.after,
                mutation.resealed_identity.hex(),
                mutation.stale_refusal,
                mutation.resealed_refusal,
                mutation.self_reference_refusal,
                mutation.semantic_mismatch,
            )
        })
        .collect::<Vec<_>>()
        .join(",");
    let json = format!(
        concat!(
            "{{\"engine\":\"{}\",\"optimizer_seed\":null,",
            "\"event_seed_semantics\":\"fixed-input-null-sentinel\",",
            "\"reference_identity\":\"{}\",\"scaled_kkt_bits\":\"0x{:016x}\",",
            "\"mutations\":[{}]}}"
        ),
        engine.name(),
        reference.declared_identity.hex(),
        reference.payload.oracle.scaled_kkt_bits,
        mutation_json,
    );
    let mut emitter = Emitter::new(SUITE, engine.name());
    let receipt_event = emitter.emit(
        Severity::Info,
        EventKind::Custom {
            name: "constrained-study-replay-receipt".to_string(),
            json,
        },
        None,
    );
    let receipt_line = receipt_event.to_jsonl();
    fs_obs::validate_line(&receipt_line)
        .expect("constrained study receipt must use the fs-obs wire schema");
    let retained = receipt_event.content_identity_receipt();
    receipt_event
        .admit_content_identity(&retained)
        .expect("fresh constrained receipt identity must admit exactly");
    println!("{receipt_line}");

    let verdict = emitter.emit(
        Severity::Info,
        EventKind::ConformanceCase {
            suite: SUITE.to_string(),
            case: engine.name().to_string(),
            pass: true,
            detail: format!(
                "fixed deterministic input (no optimizer seed; event seed zero is a null schema sentinel) replayed the complete {} callback/report/oracle receipt; {} engine/lane-namespaced mutation seeds causally selected only within-lane coordinates/directions, and every lane reproduced stable stale PayloadIdentityMismatch, resealed ReferenceIdentityMismatch, and mutant-as-own-reference SemanticOracleMismatch refusals",
                engine.name(),
                mutations.len(),
            ),
            seed: EVENT_SEED_NONE,
        },
        None,
    );
    fs_obs::lint_failure_record(&verdict)
        .expect("constrained seeded-failure verdict must be replayable");
    let verdict_line = verdict.to_jsonl();
    fs_obs::validate_line(&verdict_line)
        .expect("constrained verdict must use the fs-obs wire schema");
    println!("{verdict_line}");
}

fn assert_quality(receipt: &RetainedReceipt) {
    assert_eq!(
        semantic_mismatch(&receipt.payload),
        None,
        "{} failed independent objective/constraint/KKT acceptance",
        receipt.payload.engine.name(),
    );
}

#[test]
fn constrained_families_replay_and_reject_seeded_red_mutations() {
    for engine in Engine::ALL {
        let config = StudyConfig::for_engine(engine);
        let reference_run = run_once(engine, &config);
        let reference = receipt(engine, &config, &reference_run);
        assert_quality(&reference);
        admit_receipt(&reference.declared_identity, &reference)
            .expect("the internally consistent reference receipt must admit");

        let replay = receipt(engine, &config, &run_once(engine, &config));
        assert_quality(&replay);
        assert_eq!(
            replay,
            reference,
            "{} callback trace and report did not replay",
            engine.name(),
        );
        assert_eq!(
            admit_receipt(&reference.declared_identity, &replay),
            Ok(()),
            "{} replay must admit against the retained reference",
            engine.name(),
        );

        let mut mutations = Vec::with_capacity(MutationKind::ALL.len());
        for kind in MutationKind::ALL {
            let mutation = mutation_proof(&reference, kind);
            let repeated = mutation_proof(&reference, kind);
            assert_eq!(
                mutation,
                repeated,
                "{} {} stale/resealed mutation proof was not deterministic",
                engine.name(),
                kind.name(),
            );
            assert_eq!(
                mutation.mutation_seed,
                MUTATION_SEED ^ engine.tag() ^ kind.tag(),
                "mutation seed namespace must be derived from the explicit engine and lane"
            );
            assert!(
                matches!(
                    mutation.stale_refusal,
                    MergeRefusal::PayloadIdentityMismatch { declared, computed }
                        if declared == reference.declared_identity.root()
                            && computed == mutation.resealed_identity.root()
                ),
                "{} {} stale payload reached reference admission",
                engine.name(),
                kind.name(),
            );
            assert!(
                matches!(
                    mutation.resealed_refusal,
                    MergeRefusal::ReferenceIdentityMismatch { expected, found }
                        if expected == reference.declared_identity.root()
                            && found == mutation.resealed_identity.root()
                ),
                "{} {} resealed mutant matched the retained reference",
                engine.name(),
                kind.name(),
            );
            assert_eq!(
                mutation.self_reference_refusal,
                MergeRefusal::SemanticOracleMismatch,
                "{} {} resealed semantic mutant self-authorized",
                engine.name(),
                kind.name(),
            );
            assert_ne!(
                mutation.resealed_identity, reference.declared_identity,
                "semantic mutation must move the canonical receipt identity"
            );
            assert_ne!(mutation.before, mutation.after);
            assert!(
                !mutation.semantic_mismatch.is_empty(),
                "red mutation must retain its independent semantic diagnosis"
            );
            mutations.push(mutation);
        }

        emit_receipt(engine, &reference, &mutations);
    }
}
