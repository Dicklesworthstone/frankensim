//! G5 full-study replay and seeded-failure self-test for Pareto tracing.
//!
//! The production weighted-sum continuation traces the known convex quadratic
//! front, then production epsilon-constraint continuation traces the known
//! concave Fonseca-Fleming front. The retained receipt binds both schedules,
//! starts, every objective callback input/value/gradient, and every public
//! `ParetoPoint` decision/objective/gradient/KKT field. A separately executed
//! repeat must reproduce the receipt byte for byte. A test-local semantic
//! oracle also recomputes both objective/gradient pairs from every returned
//! decision, checks each point against its declared schedule entry, and
//! establishes the existence of a nonnegative one-constraint KKT witness
//! independently of the production report.
//! Deterministic red mutations cover a finite returned-decision bit flip,
//! point/schedule permutations, callback phase/order and returned-point pair
//! coverage, and both invalid and finite-positive corrupted public report
//! fields; even self-consistently resealed forms must fail closed.
//!
//! This is the two-objective tracing family only. It does not claim
//! tri-objective behavior, the full WFG transformation stack, cancellation,
//! checkpointing, cross-ISA equality, persistence, or performance.

use core::cell::RefCell;

use fs_ascent::{ParetoPoint, epsilon_constraint_sweep, weighted_sum_sweep};
use fs_obs::ident::{IdentityBuilder, ReplayIdentity, check_version};
use fs_obs::{Emitter, EventKind, Severity};

const SUITE: &str = "fs-ascent/pareto-study-replay";
const MUTATION_SEED: u64 = 0x5041_5245_544F_5244;
const EPSILON_TOLERANCE: f64 = 1e-7;
const WEIGHTED_GRADIENT_LIMIT: f64 = 1e-9;
const CLOSED_FORM_LIMIT: f64 = 1e-7;
const KKT_RESIDUAL_LIMIT: f64 = 1e-5;
const PARETO_SET_ALIGNMENT_LIMIT: f64 = 1e-4;
const EPSILON_COVERAGE_MINIMUM_SPREAD: f64 = 0.6;
const KKT_ROUNDOFF_SCALE: f64 = 128.0;
const WEIGHTED_DIMENSION: usize = 3;
const EPSILON_DIMENSION: usize = 2;
const SEMANTIC_ORACLE_VERSION: &str =
    "pareto-independent-objective-gradient-paired-callback-schedule-kkt-witness-v2";
const OBJECTIVE_FIXTURE_VERSION: &str =
    "quadratic-three-dimensional-v1+fonseca-fleming-two-dimensional-v1";
const UNIT_CONVENTION: &str = "dimensionless-decisions-objectives-gradients-and-kkt-residuals";
const WEIGHTED_START: [f64; 3] = [0.5, 0.5, 0.5];
const EPSILON_START: [f64; 2] = [0.0, 0.0];

#[derive(Clone, Debug, PartialEq, Eq)]
struct ObjectiveCall {
    phase: &'static str,
    objective: &'static str,
    point_bits: Vec<u64>,
    value_bits: u64,
    gradient_bits: Vec<u64>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct PointPayload {
    x_bits: Vec<u64>,
    objective_bits: [u64; 2],
    kkt_bits: Option<[u64; 4]>,
    gradient_norm_bits: u64,
}

impl From<&ParetoPoint> for PointPayload {
    fn from(point: &ParetoPoint) -> Self {
        Self {
            x_bits: bits(&point.x),
            objective_bits: [point.f[0].to_bits(), point.f[1].to_bits()],
            kkt_bits: point.kkt.as_ref().map(|kkt| {
                [
                    kkt.stationarity.to_bits(),
                    kkt.feasibility.to_bits(),
                    kkt.dual_feasibility.to_bits(),
                    kkt.complementarity.to_bits(),
                ]
            }),
            gradient_norm_bits: point.grad_norm.to_bits(),
        }
    }
}

#[derive(Debug)]
struct RunRecord {
    objective_calls: Vec<ObjectiveCall>,
    weighted_points: Vec<PointPayload>,
    epsilon_points: Vec<PointPayload>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct ReceiptPayload {
    weights_bits: Vec<u64>,
    weighted_start_bits: Vec<u64>,
    epsilon_bits: Vec<u64>,
    epsilon_start_bits: Vec<u64>,
    objective_calls: Vec<ObjectiveCall>,
    weighted_points: Vec<PointPayload>,
    epsilon_points: Vec<PointPayload>,
}

impl ReceiptPayload {
    fn identity(&self) -> ReplayIdentity {
        let mut builder = IdentityBuilder::new("fs-ascent-pareto-study-receipt-v2")
            .str("suite", SUITE)
            .str("fs-ascent-version", fs_ascent::VERSION)
            .str("fs-math-version", fs_math::VERSION)
            .str("fs-obs-version", fs_obs::VERSION)
            .str("semantic-oracle-version", SEMANTIC_ORACLE_VERSION)
            .str("objective-fixture-version", OBJECTIVE_FIXTURE_VERSION)
            .str("unit-convention", UNIT_CONVENTION)
            .str("family", "two-objective-pareto-tracing")
            .str("weighted-engine", "weighted_sum_sweep/L-BFGS")
            .str(
                "epsilon-engine",
                "epsilon_constraint_sweep/augmented-lagrangian",
            )
            .u64("weighted-dimension", WEIGHTED_DIMENSION as u64)
            .u64("epsilon-dimension", EPSILON_DIMENSION as u64)
            .u64("mutation-seed", MUTATION_SEED)
            .f64_bits("epsilon-tolerance", EPSILON_TOLERANCE)
            .f64_bits("weighted-gradient-limit", WEIGHTED_GRADIENT_LIMIT)
            .f64_bits("closed-form-limit", CLOSED_FORM_LIMIT)
            .f64_bits("kkt-residual-limit", KKT_RESIDUAL_LIMIT)
            .f64_bits("pareto-set-alignment-limit", PARETO_SET_ALIGNMENT_LIMIT)
            .f64_bits(
                "epsilon-coverage-minimum-spread",
                EPSILON_COVERAGE_MINIMUM_SPREAD,
            )
            .f64_bits("kkt-roundoff-scale", KKT_ROUNDOFF_SCALE)
            .u64("weights", self.weights_bits.len() as u64);
        for &value_bits in &self.weights_bits {
            builder = builder.u64("weight-bits", value_bits);
        }
        builder = builder.u64(
            "weighted-start-values",
            self.weighted_start_bits.len() as u64,
        );
        for &value_bits in &self.weighted_start_bits {
            builder = builder.u64("weighted-start-bits", value_bits);
        }
        builder = builder.u64("epsilons", self.epsilon_bits.len() as u64);
        for &value_bits in &self.epsilon_bits {
            builder = builder.u64("epsilon-bits", value_bits);
        }
        builder = builder.u64("epsilon-start-values", self.epsilon_start_bits.len() as u64);
        for &value_bits in &self.epsilon_start_bits {
            builder = builder.u64("epsilon-start-bits", value_bits);
        }

        builder = builder.u64("objective-calls", self.objective_calls.len() as u64);
        for call in &self.objective_calls {
            builder = builder
                .str("call-phase", call.phase)
                .str("call-objective", call.objective)
                .u64("call-point-values", call.point_bits.len() as u64);
            for &value_bits in &call.point_bits {
                builder = builder.u64("call-point-bits", value_bits);
            }
            builder = builder
                .u64("call-value-bits", call.value_bits)
                .u64("call-gradient-values", call.gradient_bits.len() as u64);
            for &value_bits in &call.gradient_bits {
                builder = builder.u64("call-gradient-bits", value_bits);
            }
        }

        builder = append_points(builder, "weighted", &self.weighted_points);
        append_points(builder, "epsilon", &self.epsilon_points).finish()
    }
}

fn append_points(
    mut builder: IdentityBuilder,
    path: &'static str,
    points: &[PointPayload],
) -> IdentityBuilder {
    builder = builder
        .str("point-path", path)
        .u64("points", points.len() as u64);
    for point in points {
        builder = builder.u64("point-x-values", point.x_bits.len() as u64);
        for &value_bits in &point.x_bits {
            builder = builder.u64("point-x-bits", value_bits);
        }
        for &value_bits in &point.objective_bits {
            builder = builder.u64("point-objective-bits", value_bits);
        }
        builder = builder.flag("point-has-kkt", point.kkt_bits.is_some());
        if let Some(kkt_bits) = point.kkt_bits {
            for value_bits in kkt_bits {
                builder = builder.u64("point-kkt-bits", value_bits);
            }
        }
        builder = builder.u64("point-gradient-norm-bits", point.gradient_norm_bits);
    }
    builder
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
    if &candidate.declared_identity != reference {
        return Err(MergeRefusal::ReferenceIdentityMismatch);
    }
    Ok(())
}

fn bits(values: &[f64]) -> Vec<u64> {
    values.iter().map(|value| value.to_bits()).collect()
}

#[derive(Clone, Debug)]
struct ObjectivePair {
    values: [f64; 2],
    gradients: [Vec<f64>; 2],
}

impl ObjectivePair {
    fn is_finite(&self) -> bool {
        self.values.iter().all(|value| value.is_finite())
            && self
                .gradients
                .iter()
                .flatten()
                .all(|value| value.is_finite())
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum SemanticRefusal {
    FixtureMetadataMismatch,
    ScheduleMismatch,
    DecisionDimensionMismatch,
    NonFiniteEvidence,
    ObjectiveMismatch,
    CallbackMismatch,
    CallbackOrderMismatch,
    CallbackCoverageMissing,
    ReturnedPointNotEvaluated,
    UnexpectedWeightedKkt,
    MissingEpsilonKkt,
    InvalidCertificateResidual,
    GradientResidualMismatch,
    EpsilonInfeasible,
    KktResidualMismatch,
    QualityRegression,
}

fn decode_finite(bits: &[u64], dimension: usize) -> Result<Vec<f64>, SemanticRefusal> {
    if bits.len() != dimension {
        return Err(SemanticRefusal::DecisionDimensionMismatch);
    }
    let values: Vec<f64> = bits.iter().map(|&value| f64::from_bits(value)).collect();
    if !values.iter().all(|value| value.is_finite()) {
        return Err(SemanticRefusal::NonFiniteEvidence);
    }
    Ok(values)
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

fn comparison_slack(values: &[f64]) -> f64 {
    let scale = values
        .iter()
        .map(|value| value.abs())
        .fold(1.0f64, f64::max);
    KKT_ROUNDOFF_SCALE * f64::EPSILON * scale
}

fn relative_product_slack(values: &[f64]) -> f64 {
    let scale = values
        .iter()
        .map(|value| value.abs())
        .fold(0.0f64, f64::max);
    (KKT_ROUNDOFF_SCALE * f64::EPSILON * scale).max(KKT_ROUNDOFF_SCALE * f64::from_bits(1))
}

fn quadratic_pair(x: &[f64]) -> ObjectivePair {
    let f1: f64 = x.iter().map(|coordinate| coordinate * coordinate).sum();
    let f2: f64 = x
        .iter()
        .map(|coordinate| (coordinate - 1.0) * (coordinate - 1.0))
        .sum();
    ObjectivePair {
        values: [f1, f2],
        gradients: [
            x.iter().map(|coordinate| 2.0 * coordinate).collect(),
            x.iter()
                .map(|coordinate| 2.0 * (coordinate - 1.0))
                .collect(),
        ],
    }
}

fn fonseca_fleming_pair(x: &[f64]) -> ObjectivePair {
    let center = 1.0 / fs_math::det::sqrt(2.0);
    let squared1: f64 = x
        .iter()
        .map(|coordinate| (coordinate - center) * (coordinate - center))
        .sum();
    let squared2: f64 = x
        .iter()
        .map(|coordinate| (coordinate + center) * (coordinate + center))
        .sum();
    let exponential1 = fs_math::det::exp(-squared1);
    let exponential2 = fs_math::det::exp(-squared2);
    ObjectivePair {
        values: [1.0 - exponential1, 1.0 - exponential2],
        gradients: [
            x.iter()
                .map(|coordinate| 2.0 * (coordinate - center) * exponential1)
                .collect(),
            x.iter()
                .map(|coordinate| 2.0 * (coordinate + center) * exponential2)
                .collect(),
        ],
    }
}

fn validate_callback_receipt(calls: &[ObjectiveCall]) -> Result<(), SemanticRefusal> {
    let mut seen = [false; 4];
    let mut epsilon_phase_started = false;
    let mut weighted_expects_f1 = true;
    for call in calls {
        match call.phase {
            "weighted" => {
                if epsilon_phase_started
                    || (weighted_expects_f1 && call.objective != "f1")
                    || (!weighted_expects_f1 && call.objective != "f2")
                {
                    return Err(SemanticRefusal::CallbackOrderMismatch);
                }
                weighted_expects_f1 = !weighted_expects_f1;
            }
            "epsilon" => epsilon_phase_started = true,
            _ => return Err(SemanticRefusal::CallbackMismatch),
        }
        let (slot, dimension, objective, expected) = match (call.phase, call.objective) {
            ("weighted", "f1") => (
                0,
                WEIGHTED_DIMENSION,
                0,
                quadratic_pair as fn(&[f64]) -> ObjectivePair,
            ),
            ("weighted", "f2") => (
                1,
                WEIGHTED_DIMENSION,
                1,
                quadratic_pair as fn(&[f64]) -> ObjectivePair,
            ),
            ("epsilon", "f1") => (
                2,
                EPSILON_DIMENSION,
                0,
                fonseca_fleming_pair as fn(&[f64]) -> ObjectivePair,
            ),
            ("epsilon", "f2") => (
                3,
                EPSILON_DIMENSION,
                1,
                fonseca_fleming_pair as fn(&[f64]) -> ObjectivePair,
            ),
            _ => return Err(SemanticRefusal::CallbackMismatch),
        };
        let point = decode_finite(&call.point_bits, dimension)?;
        let oracle = expected(&point);
        if !oracle.is_finite() {
            return Err(SemanticRefusal::NonFiniteEvidence);
        }
        if call.value_bits != oracle.values[objective].to_bits()
            || call.gradient_bits != bits(&oracle.gradients[objective])
        {
            return Err(SemanticRefusal::CallbackMismatch);
        }
        seen[slot] = true;
    }
    if !weighted_expects_f1 {
        return Err(SemanticRefusal::CallbackOrderMismatch);
    }
    if !seen.into_iter().all(|was_seen| was_seen) {
        return Err(SemanticRefusal::CallbackCoverageMissing);
    }
    Ok(())
}

fn returned_point_evaluation_pairs(
    calls: &[ObjectiveCall],
    phase: &'static str,
    point_bits: &[u64],
) -> usize {
    calls
        .windows(2)
        .filter(|pair| {
            pair[0].phase == phase
                && pair[1].phase == phase
                && pair[0].objective == "f1"
                && pair[1].objective == "f2"
                && pair[0].point_bits.as_slice() == point_bits
                && pair[1].point_bits.as_slice() == point_bits
        })
        .count()
}

fn record_call(
    calls: &RefCell<Vec<ObjectiveCall>>,
    phase: &'static str,
    objective: &'static str,
    point: &[f64],
    value: f64,
    gradient: &[f64],
) {
    calls.borrow_mut().push(ObjectiveCall {
        phase,
        objective,
        point_bits: bits(point),
        value_bits: value.to_bits(),
        gradient_bits: bits(gradient),
    });
}

fn weights() -> Vec<f64> {
    (1..10).map(|index| f64::from(index) / 10.0).collect()
}

fn epsilons() -> Vec<f64> {
    (0..8)
        .map(|index| 0.1f64.mul_add(f64::from(index), 0.15))
        .collect()
}

fn validate_weighted_points(
    points: &[PointPayload],
    schedule: &[f64],
    calls: &[ObjectiveCall],
) -> Result<(), SemanticRefusal> {
    if points.len() != schedule.len() {
        return Err(SemanticRefusal::ScheduleMismatch);
    }
    let mut worst_closed_form = 0.0f64;
    for (index, (point, &weight)) in points.iter().zip(schedule).enumerate() {
        let x = decode_finite(&point.x_bits, WEIGHTED_DIMENSION)?;
        let required_occurrences = points[..=index]
            .iter()
            .filter(|candidate| candidate.x_bits.as_slice() == point.x_bits.as_slice())
            .count();
        if returned_point_evaluation_pairs(calls, "weighted", &point.x_bits) < required_occurrences
        {
            return Err(SemanticRefusal::ReturnedPointNotEvaluated);
        }
        let oracle = quadratic_pair(&x);
        if !oracle.is_finite() {
            return Err(SemanticRefusal::NonFiniteEvidence);
        }
        if point.objective_bits != [oracle.values[0].to_bits(), oracle.values[1].to_bits()] {
            return Err(SemanticRefusal::ObjectiveMismatch);
        }
        if point.kkt_bits.is_some() {
            return Err(SemanticRefusal::UnexpectedWeightedKkt);
        }
        let reported_gradient = f64::from_bits(point.gradient_norm_bits);
        if !reported_gradient.is_finite() || reported_gradient < 0.0 {
            return Err(SemanticRefusal::InvalidCertificateResidual);
        }
        let scalarized_gradient: Vec<f64> = oracle.gradients[0]
            .iter()
            .zip(&oracle.gradients[1])
            .map(|(g1, g2)| weight.mul_add(*g1, (1.0 - weight) * g2))
            .collect();
        if point.gradient_norm_bits != inf_norm(&scalarized_gradient).to_bits() {
            return Err(SemanticRefusal::GradientResidualMismatch);
        }
        if reported_gradient >= WEIGHTED_GRADIENT_LIMIT {
            return Err(SemanticRefusal::QualityRegression);
        }

        let expected_f1 = 3.0 * (1.0 - weight) * (1.0 - weight);
        let expected_f2 = 3.0 * weight * weight;
        worst_closed_form = worst_closed_form
            .max((oracle.values[0] - expected_f1).abs())
            .max((oracle.values[1] - expected_f2).abs());
    }
    if worst_closed_form >= CLOSED_FORM_LIMIT {
        return Err(SemanticRefusal::QualityRegression);
    }
    Ok(())
}

fn reconstruct_epsilon_kkt(
    oracle: &ObjectivePair,
    epsilon: f64,
) -> Result<([f64; 4], f64), SemanticRefusal> {
    let constraint_gradient = &oracle.gradients[0];
    let objective_gradient = &oracle.gradients[1];
    let denominator = dot(constraint_gradient, constraint_gradient);
    if !denominator.is_finite() || denominator <= 0.0 {
        return Err(SemanticRefusal::KktResidualMismatch);
    }

    // This nonnegative least-squares multiplier is independent of the
    // augmented-Lagrangian multiplier hidden inside the production report. It
    // is the unique one-row KKT witness minimizing the stationarity 2-norm.
    let multiplier = (-dot(objective_gradient, constraint_gradient) / denominator).max(0.0);
    let stationarity_vector: Vec<f64> = objective_gradient
        .iter()
        .zip(constraint_gradient)
        .map(|(objective, constraint)| multiplier.mul_add(*constraint, *objective))
        .collect();
    let constraint = oracle.values[0] - epsilon;
    let residuals = [
        inf_norm(&stationarity_vector),
        constraint.max(0.0),
        (-multiplier).max(0.0),
        (multiplier * constraint).abs(),
    ];
    if !multiplier.is_finite()
        || residuals
            .iter()
            .any(|residual| !residual.is_finite() || *residual < 0.0)
    {
        return Err(SemanticRefusal::InvalidCertificateResidual);
    }
    Ok((residuals, multiplier))
}

fn public_kkt_witness_exists(oracle: &ObjectivePair, constraint: f64, reported: [f64; 4]) -> bool {
    let objective_gradient = &oracle.gradients[1];
    let constraint_gradient = &oracle.gradients[0];
    let mut candidates = Vec::with_capacity(2 * objective_gradient.len() + 2);
    candidates.push(0.0);

    // When the constraint is not exactly active, complementarity identifies
    // the only possible nonnegative multiplier up to floating-point roundoff.
    if constraint != 0.0 {
        candidates.push(reported[3] / constraint.abs());
    }

    // If the constraint is exactly active, complementarity is identically
    // zero and does not identify the multiplier.  At any witness for which
    // ||g2 + nu*g1||_inf == reported stationarity, at least one component
    // attains either +reported or -reported.  Solving those two affine
    // equations for every nonzero row gives a complete finite candidate set.
    for (&objective, &constraint_row) in objective_gradient.iter().zip(constraint_gradient) {
        if constraint_row != 0.0 {
            candidates.push((reported[0] - objective) / constraint_row);
            candidates.push((-reported[0] - objective) / constraint_row);
        }
    }

    candidates.into_iter().any(|multiplier| {
        if !multiplier.is_finite() || multiplier < 0.0 {
            return false;
        }
        let multiplier_gradient: Vec<f64> = constraint_gradient
            .iter()
            .map(|gradient| multiplier * gradient)
            .collect();
        if multiplier_gradient.iter().any(|value| !value.is_finite()) {
            return false;
        }
        let stationarity_vector: Vec<f64> = objective_gradient
            .iter()
            .zip(&multiplier_gradient)
            .map(|(objective, constraint)| objective + constraint)
            .collect();
        let stationarity = inf_norm(&stationarity_vector);
        let complementarity = (multiplier * constraint).abs();
        if !stationarity.is_finite() || !complementarity.is_finite() {
            return false;
        }
        let stationarity_slack = comparison_slack(&[
            reported[0],
            stationarity,
            inf_norm(objective_gradient),
            inf_norm(&multiplier_gradient),
        ]);
        let complementarity_slack = relative_product_slack(&[reported[3], complementarity]);
        (reported[0] - stationarity).abs() <= stationarity_slack
            && (reported[3] - complementarity).abs() <= complementarity_slack
    })
}

fn validate_epsilon_points(
    points: &[PointPayload],
    schedule: &[f64],
    calls: &[ObjectiveCall],
) -> Result<(), SemanticRefusal> {
    if points.len() != schedule.len() {
        return Err(SemanticRefusal::ScheduleMismatch);
    }
    let mut lowest_f1 = f64::INFINITY;
    let mut highest_f1 = f64::NEG_INFINITY;
    for (index, (point, &epsilon)) in points.iter().zip(schedule).enumerate() {
        let x = decode_finite(&point.x_bits, EPSILON_DIMENSION)?;
        let required_occurrences = points[..=index]
            .iter()
            .filter(|candidate| candidate.x_bits.as_slice() == point.x_bits.as_slice())
            .count();
        if returned_point_evaluation_pairs(calls, "epsilon", &point.x_bits) < required_occurrences {
            return Err(SemanticRefusal::ReturnedPointNotEvaluated);
        }
        let oracle = fonseca_fleming_pair(&x);
        if !oracle.is_finite() {
            return Err(SemanticRefusal::NonFiniteEvidence);
        }
        if point.objective_bits != [oracle.values[0].to_bits(), oracle.values[1].to_bits()] {
            return Err(SemanticRefusal::ObjectiveMismatch);
        }
        if oracle.values[0] > epsilon + EPSILON_TOLERANCE {
            return Err(SemanticRefusal::EpsilonInfeasible);
        }
        if (x[0] - x[1]).abs() >= PARETO_SET_ALIGNMENT_LIMIT {
            return Err(SemanticRefusal::QualityRegression);
        }

        let kkt_bits = point.kkt_bits.ok_or(SemanticRefusal::MissingEpsilonKkt)?;
        let reported = kkt_bits.map(f64::from_bits);
        if reported
            .iter()
            .any(|residual| !residual.is_finite() || *residual < 0.0)
        {
            return Err(SemanticRefusal::InvalidCertificateResidual);
        }
        let reported_gradient = f64::from_bits(point.gradient_norm_bits);
        if !reported_gradient.is_finite() || reported_gradient < 0.0 {
            return Err(SemanticRefusal::InvalidCertificateResidual);
        }
        if point.gradient_norm_bits != kkt_bits[0] {
            return Err(SemanticRefusal::GradientResidualMismatch);
        }
        if reported
            .iter()
            .any(|residual| *residual >= KKT_RESIDUAL_LIMIT)
        {
            return Err(SemanticRefusal::QualityRegression);
        }

        let (independent, _independent_multiplier) = reconstruct_epsilon_kkt(&oracle, epsilon)?;
        let constraint = oracle.values[0] - epsilon;

        // Feasibility and dual feasibility have no hidden-multiplier
        // ambiguity in this one-inequality production path. They are computed
        // from the same retained objective bits and the production multiplier
        // is projected nonnegative, so require exact canonical bits rather
        // than allowing an arbitrary small positive lie through a tolerance.
        if kkt_bits[1] != independent[1].to_bits() || kkt_bits[2] != independent[2].to_bits() {
            return Err(SemanticRefusal::KktResidualMismatch);
        }

        // The public API does not expose the augmented-Lagrangian multiplier,
        // but stationarity and complementarity must still be jointly
        // realizable by at least one finite nonnegative multiplier. Checking
        // witness existence closes the old loophole where each public field
        // could be perturbed independently while remaining below the quality
        // cap and satisfying only loose necessary bounds.
        if !public_kkt_witness_exists(&oracle, constraint, reported) {
            return Err(SemanticRefusal::KktResidualMismatch);
        }

        // The reconstructed multiplier minimizes the stationarity 2-norm over
        // all ν >= 0. Any valid production multiplier with reported infinity
        // norm r therefore implies ||r_independent||∞ <= sqrt(n) * r.
        let dimension_factor = fs_math::det::sqrt(2.0);
        let lower_bound_slack = comparison_slack(&[
            independent[0],
            reported[0],
            inf_norm(&oracle.gradients[0]),
            inf_norm(&oracle.gradients[1]),
        ]);
        if independent[0] > dimension_factor.mul_add(reported[0], lower_bound_slack) {
            return Err(SemanticRefusal::KktResidualMismatch);
        }

        lowest_f1 = lowest_f1.min(oracle.values[0]);
        highest_f1 = highest_f1.max(oracle.values[0]);
    }
    if highest_f1 - lowest_f1 <= EPSILON_COVERAGE_MINIMUM_SPREAD {
        return Err(SemanticRefusal::QualityRegression);
    }
    Ok(())
}

fn validate_semantics(payload: &ReceiptPayload) -> Result<(), SemanticRefusal> {
    if payload.weighted_start_bits != bits(&WEIGHTED_START)
        || payload.epsilon_start_bits != bits(&EPSILON_START)
    {
        return Err(SemanticRefusal::FixtureMetadataMismatch);
    }
    let weighted_schedule = weights();
    let epsilon_schedule = epsilons();
    if payload.weights_bits != bits(&weighted_schedule)
        || payload.epsilon_bits != bits(&epsilon_schedule)
    {
        return Err(SemanticRefusal::ScheduleMismatch);
    }
    validate_callback_receipt(&payload.objective_calls)?;
    validate_weighted_points(
        &payload.weighted_points,
        &weighted_schedule,
        &payload.objective_calls,
    )?;
    validate_epsilon_points(
        &payload.epsilon_points,
        &epsilon_schedule,
        &payload.objective_calls,
    )
}

fn run_once() -> RunRecord {
    let calls = RefCell::new(Vec::new());
    let weighted_points = {
        let f1 = |x: &[f64]| {
            let value: f64 = x.iter().map(|coordinate| coordinate * coordinate).sum();
            let gradient: Vec<f64> = x.iter().map(|coordinate| 2.0 * coordinate).collect();
            record_call(&calls, "weighted", "f1", x, value, &gradient);
            (value, gradient)
        };
        let f2 = |x: &[f64]| {
            let value: f64 = x
                .iter()
                .map(|coordinate| (coordinate - 1.0) * (coordinate - 1.0))
                .sum();
            let gradient: Vec<f64> = x
                .iter()
                .map(|coordinate| 2.0 * (coordinate - 1.0))
                .collect();
            record_call(&calls, "weighted", "f2", x, value, &gradient);
            (value, gradient)
        };
        weighted_sum_sweep(&f1, &f2, &weights(), &WEIGHTED_START)
            .iter()
            .map(PointPayload::from)
            .collect()
    };

    let epsilon_points = {
        let center = 1.0 / fs_math::det::sqrt(2.0);
        let f1 = |x: &[f64]| {
            let squared: f64 = x
                .iter()
                .map(|coordinate| (coordinate - center) * (coordinate - center))
                .sum();
            let exponential = fs_math::det::exp(-squared);
            let value = 1.0 - exponential;
            let gradient: Vec<f64> = x
                .iter()
                .map(|coordinate| 2.0 * (coordinate - center) * exponential)
                .collect();
            record_call(&calls, "epsilon", "f1", x, value, &gradient);
            (value, gradient)
        };
        let f2 = |x: &[f64]| {
            let squared: f64 = x
                .iter()
                .map(|coordinate| (coordinate + center) * (coordinate + center))
                .sum();
            let exponential = fs_math::det::exp(-squared);
            let value = 1.0 - exponential;
            let gradient: Vec<f64> = x
                .iter()
                .map(|coordinate| 2.0 * (coordinate + center) * exponential)
                .collect();
            record_call(&calls, "epsilon", "f2", x, value, &gradient);
            (value, gradient)
        };
        epsilon_constraint_sweep(&f1, &f2, &epsilons(), &EPSILON_START, EPSILON_TOLERANCE)
            .iter()
            .map(PointPayload::from)
            .collect()
    };

    RunRecord {
        objective_calls: calls.into_inner(),
        weighted_points,
        epsilon_points,
    }
}

fn receipt(run: &RunRecord) -> RetainedReceipt {
    RetainedReceipt::new(ReceiptPayload {
        weights_bits: bits(&weights()),
        weighted_start_bits: bits(&WEIGHTED_START),
        epsilon_bits: bits(&epsilons()),
        epsilon_start_bits: bits(&EPSILON_START),
        objective_calls: run.objective_calls.clone(),
        weighted_points: run.weighted_points.clone(),
        epsilon_points: run.epsilon_points.clone(),
    })
}

fn mutate_returned_decision(receipt: &RetainedReceipt) -> (RetainedReceipt, usize, usize, u64) {
    let mut mutant = receipt.clone();
    let point = (MUTATION_SEED as usize) % mutant.payload.epsilon_points.len();
    let coordinate =
        ((MUTATION_SEED >> 8) as usize) % mutant.payload.epsilon_points[point].x_bits.len();
    let mask = 1_u64 << ((MUTATION_SEED >> 16) % 52);
    mutant.payload.epsilon_points[point].x_bits[coordinate] ^= mask;
    assert!(
        f64::from_bits(mutant.payload.epsilon_points[point].x_bits[coordinate]).is_finite(),
        "mantissa-only mutation must remain a finite wire-valid decision"
    );
    mutant.reseal();
    (mutant, point, coordinate, mask)
}

fn emit_receipt(
    reference: &RetainedReceipt,
    mutant: &RetainedReceipt,
    point: usize,
    coordinate: usize,
    mask: u64,
) {
    let json = format!(
        "{{\"fixture\":\"deterministic-two-objective-tracing\",\"mutation_seed\":{MUTATION_SEED},\
         \"reference_identity\":\"{}\",\"mutant_identity\":\"{}\",\
         \"mutated_path\":\"epsilon\",\"mutated_point\":{point},\
         \"mutated_coordinate\":{coordinate},\"mantissa_mask\":\"{mask:#018x}\",\
         \"semantic_oracle\":\"{SEMANTIC_ORACLE_VERSION}\",\
         \"semantic_red_cases\":[\"point-permutation\",\"paired-schedule-permutation\",\"callback-phase-order\",\"callback-pair-coverage\",\"callback-value\",\"callback-gradient\",\"objective-corruption\",\"negative-kkt\",\"finite-positive-kkt-components\"],\
         \"merge_refusal\":\"reference-identity-mismatch\"}}",
        reference.declared_identity.hex(),
        mutant.declared_identity.hex(),
    );
    let mut emitter = Emitter::new(SUITE, "two-objective-tracing");
    let receipt_event = emitter.emit(
        Severity::Info,
        EventKind::Custom {
            name: "pareto-study-replay-receipt".to_string(),
            json,
        },
        None,
    );
    let receipt_line = receipt_event.to_jsonl();
    fs_obs::validate_line(&receipt_line)
        .expect("Pareto study receipt must use the fs-obs wire schema");
    println!("{receipt_line}");

    let verdict = emitter.emit(
        Severity::Info,
        EventKind::ConformanceCase {
            suite: SUITE.to_string(),
            case: "two-objective-tracing".to_string(),
            pass: true,
            detail: format!(
                "the deterministic weighted and epsilon fixtures replayed every callback/result receipt; the independent objective/gradient, phase-ordered returned-point callback-pair, schedule-feasibility, and nonnegative one-row KKT-witness oracle admitted every reference point; point/schedule, callback order/coverage/value/gradient, objective-report, negative-KKT, and finite-positive KKT-component semantic red cases were refused after resealing; mutation seed {MUTATION_SEED:#018x} flipped epsilon point {point} coordinate {coordinate} mask {mask:#018x}, produced stable identity {}, and both merge gates refused it",
                mutant.declared_identity.hex(),
            ),
            seed: MUTATION_SEED,
        },
        None,
    );
    fs_obs::lint_failure_record(&verdict)
        .expect("Pareto seeded-failure verdict must be replayable");
    let verdict_line = verdict.to_jsonl();
    fs_obs::validate_line(&verdict_line).expect("Pareto verdict must use the fs-obs wire schema");
    println!("{verdict_line}");
}

fn assert_quality(payload: &ReceiptPayload) {
    assert_eq!(
        validate_semantics(payload),
        Ok(()),
        "Pareto study receipt failed its independent semantic oracle"
    );
}

fn assert_resealed_semantic_refusal(
    reference: &RetainedReceipt,
    payload: ReceiptPayload,
    expected: SemanticRefusal,
) {
    let mutant = RetainedReceipt::new(payload);
    assert_ne!(
        mutant.declared_identity, reference.declared_identity,
        "semantic red mutation must move the canonical receipt identity"
    );
    assert_eq!(
        validate_semantics(&mutant.payload),
        Err(expected),
        "self-consistently resealed semantic mutant was not refused for the expected reason"
    );
    assert_eq!(
        admit_receipt(&reference.declared_identity, &mutant),
        Err(MergeRefusal::ReferenceIdentityMismatch),
        "self-consistently resealed semantic mutant bypassed the reference identity gate"
    );
}

fn alternate_positive_residual(original: f64) -> f64 {
    let high = 0.75 * KKT_RESIDUAL_LIMIT;
    let low = 0.25 * KKT_RESIDUAL_LIMIT;
    if original.to_bits() == high.to_bits() {
        low
    } else {
        high
    }
}

#[test]
fn pareto_tracing_replays_and_rejects_seeded_red_mutation() {
    let reference_run = run_once();
    let reference = receipt(&reference_run);
    assert_quality(&reference.payload);
    admit_receipt(&reference.declared_identity, &reference)
        .expect("the internally consistent reference receipt must admit");

    let replay = receipt(&run_once());
    assert_eq!(
        replay, reference,
        "complete Pareto callback and result receipts did not replay"
    );

    let (mutant, point, coordinate, mask) = mutate_returned_decision(&reference);
    let (mutant_repeat, repeat_point, repeat_coordinate, repeat_mask) =
        mutate_returned_decision(&reference);
    assert_eq!(
        (point, coordinate, mask),
        (repeat_point, repeat_coordinate, repeat_mask)
    );
    assert_eq!(mutant, mutant_repeat, "seeded mutation was not stable");
    assert_ne!(mutant.declared_identity, reference.declared_identity);
    assert!(
        validate_semantics(&mutant.payload).is_err(),
        "returned-coordinate mutation must fail the independent semantic oracle"
    );
    let mut stale_identity_mutant = mutant.clone();
    stale_identity_mutant.declared_identity = reference.declared_identity.clone();
    assert_eq!(
        admit_receipt(&reference.declared_identity, &stale_identity_mutant),
        Err(MergeRefusal::PayloadIdentityMismatch)
    );
    assert_eq!(
        admit_receipt(&reference.declared_identity, &mutant),
        Err(MergeRefusal::ReferenceIdentityMismatch)
    );

    let mut point_permutation = reference.payload.clone();
    let last_epsilon_point = point_permutation.epsilon_points.len() - 1;
    point_permutation.epsilon_points.swap(0, last_epsilon_point);
    let point_permutation = RetainedReceipt::new(point_permutation);
    assert_ne!(
        point_permutation.declared_identity, reference.declared_identity,
        "point permutation must move the canonical receipt identity"
    );
    assert!(
        matches!(
            validate_semantics(&point_permutation.payload),
            Err(SemanticRefusal::EpsilonInfeasible) | Err(SemanticRefusal::KktResidualMismatch)
        ),
        "point permutation must fail the schedule-aware epsilon/KKT oracle"
    );

    let mut paired_permutation = reference.payload.clone();
    let last_epsilon = paired_permutation.epsilon_bits.len() - 1;
    paired_permutation.epsilon_bits.swap(0, last_epsilon);
    paired_permutation.epsilon_points.swap(0, last_epsilon);
    assert_resealed_semantic_refusal(
        &reference,
        paired_permutation,
        SemanticRefusal::ScheduleMismatch,
    );

    let mut callback_reordering = reference.payload.clone();
    let relocated_weighted_pair: Vec<ObjectiveCall> =
        callback_reordering.objective_calls.drain(0..2).collect();
    callback_reordering
        .objective_calls
        .extend(relocated_weighted_pair);
    assert_resealed_semantic_refusal(
        &reference,
        callback_reordering,
        SemanticRefusal::CallbackOrderMismatch,
    );

    let mut callback_gap = reference.payload.clone();
    let uncovered_point = callback_gap.weighted_points[0].x_bits.clone();
    callback_gap.objective_calls.retain(|call| {
        call.phase != "weighted" || call.point_bits.as_slice() != uncovered_point.as_slice()
    });
    assert_resealed_semantic_refusal(
        &reference,
        callback_gap,
        SemanticRefusal::ReturnedPointNotEvaluated,
    );

    let mut callback_value_corruption = reference.payload.clone();
    let callback = callback_value_corruption
        .objective_calls
        .iter_mut()
        .find(|call| call.phase == "epsilon" && call.objective == "f1")
        .expect("reference receipt must contain an epsilon/f1 callback");
    callback.value_bits = (f64::from_bits(callback.value_bits) + 0.25).to_bits();
    assert_resealed_semantic_refusal(
        &reference,
        callback_value_corruption,
        SemanticRefusal::CallbackMismatch,
    );

    let mut callback_gradient_corruption = reference.payload.clone();
    let callback = callback_gradient_corruption
        .objective_calls
        .iter_mut()
        .find(|call| call.phase == "epsilon" && call.objective == "f2")
        .expect("reference receipt must contain an epsilon/f2 callback");
    callback.gradient_bits[0] = (f64::from_bits(callback.gradient_bits[0]) + 0.25).to_bits();
    assert_resealed_semantic_refusal(
        &reference,
        callback_gradient_corruption,
        SemanticRefusal::CallbackMismatch,
    );

    let mut objective_corruption = reference.payload.clone();
    let original_objective =
        f64::from_bits(objective_corruption.epsilon_points[0].objective_bits[0]);
    objective_corruption.epsilon_points[0].objective_bits[0] =
        (original_objective + 0.25).to_bits();
    assert_resealed_semantic_refusal(
        &reference,
        objective_corruption,
        SemanticRefusal::ObjectiveMismatch,
    );

    let mut certificate_corruption = reference.payload.clone();
    certificate_corruption.epsilon_points[0]
        .kkt_bits
        .as_mut()
        .expect("reference epsilon point must retain KKT evidence")[0] = (-1e-12f64).to_bits();
    assert_resealed_semantic_refusal(
        &reference,
        certificate_corruption,
        SemanticRefusal::InvalidCertificateResidual,
    );

    let mut stationarity_corruption = reference.payload.clone();
    let epsilon_point = &mut stationarity_corruption.epsilon_points[0];
    let kkt_bits = epsilon_point
        .kkt_bits
        .as_mut()
        .expect("reference epsilon point must retain KKT evidence");
    let corrupted_stationarity = alternate_positive_residual(f64::from_bits(kkt_bits[0]));
    kkt_bits[0] = corrupted_stationarity.to_bits();
    epsilon_point.gradient_norm_bits = corrupted_stationarity.to_bits();
    assert_resealed_semantic_refusal(
        &reference,
        stationarity_corruption,
        SemanticRefusal::KktResidualMismatch,
    );

    let mut feasibility_corruption = reference.payload.clone();
    let kkt_bits = feasibility_corruption.epsilon_points[0]
        .kkt_bits
        .as_mut()
        .expect("reference epsilon point must retain KKT evidence");
    kkt_bits[1] = alternate_positive_residual(f64::from_bits(kkt_bits[1])).to_bits();
    assert_resealed_semantic_refusal(
        &reference,
        feasibility_corruption,
        SemanticRefusal::KktResidualMismatch,
    );

    let mut dual_feasibility_corruption = reference.payload.clone();
    let kkt_bits = dual_feasibility_corruption.epsilon_points[0]
        .kkt_bits
        .as_mut()
        .expect("reference epsilon point must retain KKT evidence");
    kkt_bits[2] = alternate_positive_residual(f64::from_bits(kkt_bits[2])).to_bits();
    assert_resealed_semantic_refusal(
        &reference,
        dual_feasibility_corruption,
        SemanticRefusal::KktResidualMismatch,
    );

    let mut complementarity_corruption = reference.payload.clone();
    let kkt_bits = complementarity_corruption.epsilon_points[0]
        .kkt_bits
        .as_mut()
        .expect("reference epsilon point must retain KKT evidence");
    kkt_bits[3] = alternate_positive_residual(f64::from_bits(kkt_bits[3])).to_bits();
    assert_resealed_semantic_refusal(
        &reference,
        complementarity_corruption,
        SemanticRefusal::KktResidualMismatch,
    );

    emit_receipt(&reference, &mutant, point, coordinate, mask);
}
