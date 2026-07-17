//! WFG4 tri-objective conformance battery (7tv.21.11).
//!
//! The fixture implements the normalized `M=3, k=4, l=20` WFG4
//! definition from the corrected WFG toolkit and composes it with the
//! production [`fs_dfo::nsga3`] engine.  The retained receipt binds the
//! complete benchmark definition, optimization budget, result front,
//! quality gates, and bitwise replay result.  This is a dimensionless
//! in-repository conformance fixture, not an external COCO/WFG campaign.

use fs_dfo::{Individual, NsgaParams, das_dennis, hypervolume, nsga3};
use fs_obs::ident::{IdentityBuilder, ReplayIdentity};
use fs_obs::{Emitter, EventKind, Severity};

const SUITE: &str = "fs-dfo-wfg4";
const CASE: &str = "wfg4-m3-k4-l20";
const INPUT_SEED: u64 = 0x5746_4734_0003_0018;
const NSGA3_STREAM_KERNEL: u32 = 0x05A3;
const NSGA3_STREAM_TILE: u32 = 0;
const JMETAL_REVISION: &str = "ea7e882f6b8f94b99535921674e62cda7986f20e";

const OBJECTIVES: usize = 3;
const POSITION_PARAMETERS: usize = 4;
const DISTANCE_PARAMETERS: usize = 20;
const DIMENSION: usize = POSITION_PARAMETERS + DISTANCE_PARAMETERS;
const REFERENCE_DIVISIONS: usize = 12;
const EXPECTED_REFERENCE_DIRECTIONS: usize = 91;
const PROBE_COUNT: usize = 8;

const S_MULTI_A: f64 = 30.0;
const S_MULTI_B: f64 = 10.0;
const S_MULTI_CENTER: f64 = 0.35;
const MUTANT_CENTER: f64 = 0.50;
const DISTANCE_SCALE: f64 = 1.0;
const OBJECTIVE_SCALES: [f64; OBJECTIVES] = [2.0, 4.0, 6.0];
const HYPERVOLUME_REFERENCE: [f64; OBJECTIVES] = [3.0, 5.0, 7.0];

const POPULATION: usize = 92;
const GENERATIONS: usize = 400;
const ETA_C: f64 = 30.0;
const ETA_M: f64 = 20.0;
const MUTATION_PROBABILITY: f64 = 1.0 / DIMENSION as f64;
const EXPECTED_EVALUATIONS: usize = POPULATION * (GENERATIONS + 1);

// Coarse conformance gates, deliberately separated from performance claims.
// The quarter-scale mean-distance ceiling was calibrated from the first
// retained central run (0.23171459361058497 at the fixed 36,892-evaluation
// budget), leaving 0.01828540638941503 absolute headroom.
const MAX_MEAN_DISTANCE: f64 = 0.25;
const MAX_WORST_DISTANCE: f64 = 0.50;
const MIN_DIRECTION_COVERAGE: f64 = 0.25;
const MIN_HYPERVOLUME: f64 = 20.0;
const MIN_FRONT_SIZE: usize = 24;
const MAX_SHAPE_RESIDUAL: f64 = 5.0e-12;
const CORRECTION_EPSILON: f64 = 1.0e-10;
const TRANSFORM_TOLERANCE: f64 = 1.0e-12;
const CANONICAL_OBJECTIVE_TOLERANCE: f64 = 1.0e-11;
const MIN_MUTANT_DISTANCE: f64 = 0.05;

#[derive(Debug, Clone, Copy)]
struct Wfg4Point {
    reduced: [f64; OBJECTIVES],
    objectives: [f64; OBJECTIVES],
}

impl Wfg4Point {
    fn distance(self) -> f64 {
        self.reduced[OBJECTIVES - 1]
    }
}

#[derive(Debug)]
struct Metrics {
    front_size: usize,
    covered_directions: usize,
    direction_coverage: f64,
    mean_distance: f64,
    worst_distance: f64,
    hypervolume: f64,
    maximum_shape_residual: f64,
    objectives_recompute_bitwise: bool,
}

#[allow(clippy::struct_excessive_bools)]
#[derive(Debug, Clone, Copy)]
struct Verdicts {
    evaluator: bool,
    mutant: bool,
    convergence: bool,
    replay: bool,
}

impl Verdicts {
    fn pass(self) -> bool {
        self.evaluator && self.mutant && self.convergence && self.replay
    }
}

fn usize_u64(value: usize) -> u64 {
    u64::try_from(value).expect("WFG4 fixture cardinality fits u64")
}

/// WFG's reference implementation snaps only roundoff-sized excursions.
fn correct_to_01(value: f64) -> f64 {
    assert!(value.is_finite(), "WFG4 transformation must stay finite");
    if (-CORRECTION_EPSILON..0.0).contains(&value) {
        0.0
    } else if (1.0..=1.0 + CORRECTION_EPSILON).contains(&value) {
        1.0
    } else {
        assert!(
            (0.0..=1.0).contains(&value),
            "WFG4 transformation escaped [0,1]: {value:.17e}"
        );
        value
    }
}

/// Canonical WFG `s_multi(y, A=30, B=10, C)` transformation.
fn s_multi_with_center(y: f64, center: f64) -> f64 {
    assert!((0.0..=1.0).contains(&y), "normalized WFG4 input");
    assert!((0.0..1.0).contains(&center), "WFG4 center");

    // `floor(center - y)` is exactly 0 on the left branch and -1 on
    // the right branch.  Spelling out that branch avoids platform libm.
    let denominator = if y <= center {
        2.0 * center
    } else {
        2.0 * (center - 1.0)
    };
    let ratio = (y - center).abs() / denominator;
    let tmp1 = (4.0 * S_MULTI_A + 2.0) * core::f64::consts::PI * (0.5 - ratio);
    let tmp2 = 4.0 * S_MULTI_B * ratio * ratio;
    correct_to_01((1.0 + fs_math::det::cos(tmp1) + tmp2) / (S_MULTI_B + 2.0))
}

fn mean(values: &[f64]) -> f64 {
    assert!(!values.is_empty(), "WFG4 reduction slice");
    values.iter().sum::<f64>() / values.len() as f64
}

fn approximately_equal(left: f64, right: f64, tolerance: f64) -> bool {
    (left - right).abs() <= tolerance
}

/// Concave WFG shape for `M=3`, indexed in objective order.
fn concave_shape(x: &[f64; OBJECTIVES]) -> [f64; OBJECTIVES] {
    let theta0 = x[0] * core::f64::consts::FRAC_PI_2;
    let theta1 = x[1] * core::f64::consts::FRAC_PI_2;
    let sin0 = fs_math::det::sin(theta0);
    let cos0 = fs_math::det::cos(theta0);
    let sin1 = fs_math::det::sin(theta1);
    let cos1 = fs_math::det::cos(theta1);
    [sin0 * sin1, sin0 * cos1, cos0]
}

fn evaluate_with_center(y: &[f64], center: f64) -> Wfg4Point {
    assert_eq!(y.len(), DIMENSION, "normalized WFG4 dimension");
    assert!(
        y.iter()
            .all(|value| value.is_finite() && (0.0..=1.0).contains(value)),
        "normalized WFG4 decisions must be finite and in [0,1]"
    );

    let transformed: Vec<f64> = y
        .iter()
        .map(|&value| s_multi_with_center(value, center))
        .collect();
    let reduced = [
        mean(&transformed[..POSITION_PARAMETERS / 2]),
        mean(&transformed[POSITION_PARAMETERS / 2..POSITION_PARAMETERS]),
        mean(&transformed[POSITION_PARAMETERS..]),
    ];

    // WFG4 uses A=[1,1], so the standard x-vector reconstruction is
    // exactly t2 for every component.
    let shape = concave_shape(&reduced);
    let distance = DISTANCE_SCALE * reduced[OBJECTIVES - 1];
    let objectives =
        core::array::from_fn(|index| OBJECTIVE_SCALES[index].mul_add(shape[index], distance));
    Wfg4Point {
        reduced,
        objectives,
    }
}

fn wfg4(y: &[f64]) -> Vec<f64> {
    evaluate_with_center(y, S_MULTI_CENTER).objectives.to_vec()
}

fn reconstructed_shape(point: Wfg4Point) -> [f64; OBJECTIVES] {
    let distance = point.distance();
    core::array::from_fn(|index| (point.objectives[index] - distance) / OBJECTIVE_SCALES[index])
}

fn maximum_shape_residual(points: &[Wfg4Point]) -> f64 {
    points.iter().fold(0.0, |maximum, point| {
        let shape = reconstructed_shape(*point);
        let squared_norm = shape.iter().map(|value| value * value).sum::<f64>();
        maximum.max((squared_norm - 1.0).abs())
    })
}

fn direction_coverage(points: &[Wfg4Point], directions: &[Vec<f64>]) -> (usize, f64) {
    let mut hit = vec![false; directions.len()];
    for point in points {
        let shape = reconstructed_shape(*point);
        let mut best = (usize::MAX, f64::INFINITY);
        for (index, direction) in directions.iter().enumerate() {
            let direction_norm = direction.iter().map(|value| value * value).sum::<f64>();
            let projection = shape
                .iter()
                .zip(direction)
                .map(|(value, axis)| value * axis)
                .sum::<f64>()
                / direction_norm;
            let squared_distance = shape
                .iter()
                .zip(direction)
                .map(|(value, axis)| {
                    let residual = projection.mul_add(-axis, *value);
                    residual * residual
                })
                .sum::<f64>();
            if squared_distance < best.1 {
                best = (index, squared_distance);
            }
        }
        hit[best.0] = true;
    }
    let covered = hit.into_iter().filter(|value| *value).count();
    (covered, covered as f64 / directions.len() as f64)
}

fn metrics(front: &[Individual], directions: &[Vec<f64>]) -> Metrics {
    let mut objectives_recompute_bitwise = true;
    let evaluated: Vec<Wfg4Point> = front
        .iter()
        .map(|individual| {
            let point = evaluate_with_center(&individual.x, S_MULTI_CENTER);
            objectives_recompute_bitwise &= individual
                .f
                .iter()
                .zip(point.objectives)
                .all(|(retained, recomputed)| retained.to_bits() == recomputed.to_bits());
            point
        })
        .collect();
    let distances: Vec<f64> = evaluated.iter().map(|point| point.distance()).collect();
    let mean_distance = distances.iter().sum::<f64>() / distances.len() as f64;
    let worst_distance = distances.iter().copied().fold(0.0, f64::max);
    let (covered_directions, direction_coverage) = direction_coverage(&evaluated, directions);
    let objectives: Vec<Vec<f64>> = front
        .iter()
        .map(|individual| individual.f.clone())
        .collect();
    Metrics {
        front_size: front.len(),
        covered_directions,
        direction_coverage,
        mean_distance,
        worst_distance,
        hypervolume: hypervolume(&objectives, &HYPERVOLUME_REFERENCE),
        maximum_shape_residual: maximum_shape_residual(&evaluated),
        objectives_recompute_bitwise,
    }
}

fn front_identity(campaign: &ReplayIdentity, front: &[Individual]) -> ReplayIdentity {
    let mut builder = IdentityBuilder::new("fs-dfo-wfg4-nsga3-front-v1")
        .child("campaign", campaign)
        .u64("front-size", usize_u64(front.len()));
    for (front_index, individual) in front.iter().enumerate() {
        builder = builder
            .u64("front-index", usize_u64(front_index))
            .u64("decision-count", usize_u64(individual.x.len()))
            .u64("objective-count", usize_u64(individual.f.len()));
        for &value in &individual.x {
            builder = builder.f64_bits("decision", value);
        }
        for &value in &individual.f {
            builder = builder.f64_bits("objective", value);
        }
    }
    builder.finish()
}

fn first_front_mismatch(left: &[Individual], right: &[Individual]) -> Option<String> {
    if left.len() != right.len() {
        return Some(format!("front-length:{}!={}", left.len(), right.len()));
    }
    for (front_index, (a, b)) in left.iter().zip(right).enumerate() {
        if a.x.len() != b.x.len() {
            return Some(format!(
                "decision-length[{front_index}]:{}!={}",
                a.x.len(),
                b.x.len()
            ));
        }
        if a.f.len() != b.f.len() {
            return Some(format!(
                "objective-length[{front_index}]:{}!={}",
                a.f.len(),
                b.f.len()
            ));
        }
        for (component, (x, y)) in a.x.iter().zip(&b.x).enumerate() {
            if x.to_bits() != y.to_bits() {
                return Some(format!(
                    "decision[{front_index},{component}]:0x{:016x}!=0x{:016x}",
                    x.to_bits(),
                    y.to_bits()
                ));
            }
        }
        for (component, (x, y)) in a.f.iter().zip(&b.f).enumerate() {
            if x.to_bits() != y.to_bits() {
                return Some(format!(
                    "objective[{front_index},{component}]:0x{:016x}!=0x{:016x}",
                    x.to_bits(),
                    y.to_bits()
                ));
            }
        }
    }
    None
}

fn campaign_identity(direction_count: usize) -> ReplayIdentity {
    IdentityBuilder::new("fs-dfo-wfg4-nsga3-config-v1")
        .str("problem", "WFG4-normalized")
        .str("source", "Huband-et-al-EMO-2005-corrected-WFG-toolkit")
        .str(
            "reference-implementation",
            "jMetal-WFG4-sMulti-rSum-concave",
        )
        .str("reference-revision", JMETAL_REVISION)
        .str("units", "dimensionless-normalized-decisions-and-objectives")
        .u64("objectives", usize_u64(OBJECTIVES))
        .u64("position-parameters-k", usize_u64(POSITION_PARAMETERS))
        .u64("distance-parameters-l", usize_u64(DISTANCE_PARAMETERS))
        .u64("dimension", usize_u64(DIMENSION))
        .f64_bits("decision-lower-bound", 0.0)
        .f64_bits("decision-upper-bound", 1.0)
        .f64_bits("s-multi-a", S_MULTI_A)
        .f64_bits("s-multi-b", S_MULTI_B)
        .f64_bits("s-multi-center", S_MULTI_CENTER)
        .f64_bits("correct-to-01-epsilon", CORRECTION_EPSILON)
        .f64_bits("transform-tolerance", TRANSFORM_TOLERANCE)
        .f64_bits(
            "canonical-objective-tolerance",
            CANONICAL_OBJECTIVE_TOLERANCE,
        )
        .f64_bits("distance-scale-d", DISTANCE_SCALE)
        .f64_bits("objective-scale-1", OBJECTIVE_SCALES[0])
        .f64_bits("objective-scale-2", OBJECTIVE_SCALES[1])
        .f64_bits("objective-scale-3", OBJECTIVE_SCALES[2])
        .str("position-reduction", "equal-rSum-over-0..2-and-2..4")
        .str("distance-reduction", "equal-rSum-over-4..24")
        .str("shape", "WFG-concave-M3")
        .str("optimizer", "fs-dfo-nsga3")
        .u64("input-seed", INPUT_SEED)
        .u64("optimizer-stream-kernel", u64::from(NSGA3_STREAM_KERNEL))
        .u64("optimizer-stream-tile", u64::from(NSGA3_STREAM_TILE))
        .u64(
            "stream-semantics-version",
            u64::from(fs_rand::STREAM_SEMANTICS_VERSION),
        )
        .u64("population", usize_u64(POPULATION))
        .u64("generations", usize_u64(GENERATIONS))
        .f64_bits("eta-c", ETA_C)
        .f64_bits("eta-m", ETA_M)
        .f64_bits("mutation-probability", MUTATION_PROBABILITY)
        .u64(
            "expected-evaluations-per-run",
            usize_u64(EXPECTED_EVALUATIONS),
        )
        .u64("replay-runs", 2)
        .u64("reference-divisions", usize_u64(REFERENCE_DIVISIONS))
        .u64("reference-direction-count", usize_u64(direction_count))
        .u64(
            "expected-reference-direction-count",
            usize_u64(EXPECTED_REFERENCE_DIRECTIONS),
        )
        .u64("deterministic-probe-count", usize_u64(PROBE_COUNT))
        .str(
            "deterministic-probe-rule",
            "y[index]=((37*index+19*probe)%101)/100",
        )
        .f64_bits("hypervolume-reference-1", HYPERVOLUME_REFERENCE[0])
        .f64_bits("hypervolume-reference-2", HYPERVOLUME_REFERENCE[1])
        .f64_bits("hypervolume-reference-3", HYPERVOLUME_REFERENCE[2])
        .f64_bits("maximum-mean-distance", MAX_MEAN_DISTANCE)
        .f64_bits("maximum-worst-distance", MAX_WORST_DISTANCE)
        .f64_bits("minimum-direction-coverage", MIN_DIRECTION_COVERAGE)
        .f64_bits("minimum-hypervolume", MIN_HYPERVOLUME)
        .u64("minimum-front-size", usize_u64(MIN_FRONT_SIZE))
        .f64_bits("maximum-shape-residual", MAX_SHAPE_RESIDUAL)
        .f64_bits("wrong-center-mutant", MUTANT_CENTER)
        .f64_bits("minimum-mutant-distance", MIN_MUTANT_DISTANCE)
        .str("capabilities", "safe-rust;strict-fs-math;keyed-fs-rand")
        .str("execution-context", "single-threaded-direct-test-no-Cx")
        .str("fs-dfo-version", fs_dfo::VERSION)
        .str("fs-math-version", fs_math::VERSION)
        .str("fs-rand-version", fs_rand::VERSION)
        .str("fs-obs-version", fs_obs::VERSION)
        .finish()
}

#[allow(clippy::too_many_arguments)]
fn result_identity(
    campaign: &ReplayIdentity,
    front_identity: &ReplayIdentity,
    replay_front_identity: &ReplayIdentity,
    canonical: Wfg4Point,
    mutant: Wfg4Point,
    metrics: &Metrics,
    evaluations: usize,
    replay_evaluations: usize,
    transform_samples: [f64; 3],
    probe_shape_residual: f64,
    replay_mismatch: Option<&str>,
    verdicts: Verdicts,
) -> ReplayIdentity {
    let mut builder = IdentityBuilder::new("fs-dfo-wfg4-nsga3-result-v1")
        .child("campaign", campaign)
        .child("front", front_identity)
        .child("replay-front", replay_front_identity)
        .flag("evaluator-pass", verdicts.evaluator)
        .flag("wrong-center-mutant-caught", verdicts.mutant)
        .flag("convergence-pass", verdicts.convergence)
        .flag("bitwise-replay-pass", verdicts.replay)
        .u64("evaluations", usize_u64(evaluations))
        .u64("replay-evaluations", usize_u64(replay_evaluations))
        .u64("front-size", usize_u64(metrics.front_size))
        .u64(
            "covered-reference-directions",
            usize_u64(metrics.covered_directions),
        )
        .f64_bits("direction-coverage", metrics.direction_coverage)
        .f64_bits("mean-distance", metrics.mean_distance)
        .f64_bits("worst-distance", metrics.worst_distance)
        .f64_bits("hypervolume", metrics.hypervolume)
        .f64_bits("maximum-shape-residual", metrics.maximum_shape_residual)
        .flag(
            "objectives-recompute-bitwise",
            metrics.objectives_recompute_bitwise,
        )
        .f64_bits("s-multi-at-zero", transform_samples[0])
        .f64_bits("s-multi-at-one", transform_samples[1])
        .f64_bits("s-multi-at-center", transform_samples[2])
        .f64_bits("probe-output-shape-residual", probe_shape_residual)
        .f64_bits("canonical-distance", canonical.distance())
        .f64_bits("mutant-distance", mutant.distance())
        .flag("replay-mismatch-present", replay_mismatch.is_some())
        .str("replay-first-mismatch", replay_mismatch.unwrap_or("none"));
    for value in canonical.objectives {
        builder = builder.f64_bits("canonical-objective", value);
    }
    for value in mutant.objectives {
        builder = builder.f64_bits("mutant-objective", value);
    }
    builder.finish()
}

fn emit_case(emitter: &mut Emitter, case: &str, pass: bool, detail: String, seed: u64) {
    let event = emitter.emit(
        if pass {
            Severity::Info
        } else {
            Severity::Error
        },
        EventKind::ConformanceCase {
            suite: SUITE.to_string(),
            case: case.to_string(),
            pass,
            detail,
            seed,
        },
        None,
    );
    fs_obs::lint_failure_record(&event).expect("WFG4 verdict must be replayable");
    let line = event.to_jsonl();
    fs_obs::validate_line(&line).expect("WFG4 verdict must use the fs-obs wire schema");
    println!("{line}");
}

fn emit_benchmark(emitter: &mut Emitter, metric: &str, value: f64) {
    let event = emitter.emit(
        Severity::Info,
        EventKind::BenchmarkResult {
            kernel: CASE.to_string(),
            metric: metric.to_string(),
            value,
            machine: 0,
        },
        None,
    );
    let line = event.to_jsonl();
    fs_obs::validate_line(&line).expect("WFG4 benchmark row must use the fs-obs wire schema");
    println!("{line}");
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn emit_receipt(
    emitter: &mut Emitter,
    campaign: &ReplayIdentity,
    result: &ReplayIdentity,
    front_identity: &ReplayIdentity,
    replay_front_identity: &ReplayIdentity,
    canonical: Wfg4Point,
    mutant: Wfg4Point,
    metrics: &Metrics,
    evaluations: usize,
    replay_evaluations: usize,
    direction_count: usize,
    transform_samples: [f64; 3],
    probe_shape_residual: f64,
    replay_mismatch: Option<&str>,
    verdicts: Verdicts,
) {
    let replay_mismatch_json =
        replay_mismatch.map_or_else(|| "null".to_string(), |value| format!("\"{value}\""));
    let event = emitter.emit(
        Severity::Info,
        EventKind::Custom {
            name: "wfg4-nsga3-conformance-receipt".to_string(),
            json: format!(
                "{{\"campaign_identity\":\"{}\",\"result_identity\":\"{}\",\
                 \"front_identity\":\"{}\",\"replay_front_identity\":\"{}\",\
                 \"problem\":\"WFG4\",\"normalized\":true,\
                 \"source_revision\":\"{JMETAL_REVISION}\",\
                 \"units\":\"dimensionless-normalized-decisions-and-objectives\",\
                 \"objectives\":{OBJECTIVES},\
                 \"k\":{POSITION_PARAMETERS},\"l\":{DISTANCE_PARAMETERS},\
                 \"dimension\":{DIMENSION},\"input_seed\":{INPUT_SEED},\
                 \"optimizer_stream_kernel\":{NSGA3_STREAM_KERNEL},\
                 \"optimizer_stream_tile\":{NSGA3_STREAM_TILE},\
                 \"stream_semantics_version\":{},\"population\":{POPULATION},\
                 \"generations\":{GENERATIONS},\"eta_c\":{ETA_C},\
                 \"eta_m\":{ETA_M},\"mutation_probability\":{MUTATION_PROBABILITY},\
                 \"evaluations\":{evaluations},\
                 \"replay_evaluations\":{replay_evaluations},\
                 \"expected_evaluations_per_run\":{EXPECTED_EVALUATIONS},\
                 \"reference_divisions\":{REFERENCE_DIVISIONS},\
                 \"reference_directions\":{direction_count},\
                 \"s_multi\":{{\"a\":{S_MULTI_A},\"b\":{S_MULTI_B},\
                 \"center\":{S_MULTI_CENTER},\
                 \"correction_epsilon\":{CORRECTION_EPSILON},\
                 \"actual_at_zero_one_center\":[{},{},{}]}},\
                 \"distance_scale\":{DISTANCE_SCALE},\
                 \"objective_scales\":[{},{},{}],\
                 \"hypervolume_reference\":[{},{},{}],\
                 \"canonical_distance\":{},\
                 \"probe_output_shape_residual\":{probe_shape_residual},\
                 \"mutant_center\":{MUTANT_CENTER},\"mutant_distance\":{},\
                 \"front_size\":{},\"covered_directions\":{},\
                 \"direction_coverage\":{},\"mean_distance\":{},\
                 \"worst_distance\":{},\"hypervolume\":{},\
                 \"maximum_shape_residual\":{},\
                 \"objectives_recompute_bitwise\":{},\
                 \"first_replay_mismatch\":{replay_mismatch_json},\"gates\":{{\
                 \"minimum_front_size\":{MIN_FRONT_SIZE},\
                 \"minimum_direction_coverage\":{MIN_DIRECTION_COVERAGE},\
                 \"maximum_mean_distance\":{MAX_MEAN_DISTANCE},\
                 \"maximum_worst_distance\":{MAX_WORST_DISTANCE},\
                 \"minimum_hypervolume\":{MIN_HYPERVOLUME},\
                 \"maximum_shape_residual\":{MAX_SHAPE_RESIDUAL},\
                 \"minimum_mutant_distance\":{MIN_MUTANT_DISTANCE}}},\
                 \"versions\":{{\"fs_dfo\":\"{}\",\"fs_math\":\"{}\",\
                 \"fs_rand\":\"{}\",\"fs_obs\":\"{}\"}},\
                 \"no_claims\":[\"all-WFG-suite\",\"external-COCO-parity\",\
                 \"executable-external-toolkit-oracle\",\"performance\",\
                 \"cross-ISA-execution\",\"Cx-cancellation\"],\
                 \"pass\":{}}}",
                campaign.hex(),
                result.hex(),
                front_identity.hex(),
                replay_front_identity.hex(),
                fs_rand::STREAM_SEMANTICS_VERSION,
                transform_samples[0],
                transform_samples[1],
                transform_samples[2],
                OBJECTIVE_SCALES[0],
                OBJECTIVE_SCALES[1],
                OBJECTIVE_SCALES[2],
                HYPERVOLUME_REFERENCE[0],
                HYPERVOLUME_REFERENCE[1],
                HYPERVOLUME_REFERENCE[2],
                canonical.distance(),
                mutant.distance(),
                metrics.front_size,
                metrics.covered_directions,
                metrics.direction_coverage,
                metrics.mean_distance,
                metrics.worst_distance,
                metrics.hypervolume,
                metrics.maximum_shape_residual,
                metrics.objectives_recompute_bitwise,
                fs_dfo::VERSION,
                fs_math::VERSION,
                fs_rand::VERSION,
                fs_obs::VERSION,
                verdicts.pass(),
            ),
        },
        None,
    );
    let line = event.to_jsonl();
    fs_obs::validate_line(&line).expect("WFG4 receipt must use the fs-obs wire schema");
    println!("{line}");
}

#[test]
#[allow(clippy::too_many_lines)]
fn wfg4_m3_nsga3_conformance_and_replay() {
    let canonical_input = vec![S_MULTI_CENTER; DIMENSION];
    let canonical = evaluate_with_center(&canonical_input, S_MULTI_CENTER);
    let mutant = evaluate_with_center(&canonical_input, MUTANT_CENTER);

    let probe_points: Vec<Wfg4Point> = (0..PROBE_COUNT)
        .map(|probe| {
            let input: Vec<f64> = (0..DIMENSION)
                .map(|index| ((index * 37 + probe * 19) % 101) as f64 / 100.0)
                .collect();
            evaluate_with_center(&input, S_MULTI_CENTER)
        })
        .collect();
    let transform_samples = [
        s_multi_with_center(0.0, S_MULTI_CENTER),
        s_multi_with_center(1.0, S_MULTI_CENTER),
        s_multi_with_center(S_MULTI_CENTER, S_MULTI_CENTER),
    ];
    let probe_shape_residual = maximum_shape_residual(&probe_points);
    let directions = das_dennis(OBJECTIVES, REFERENCE_DIVISIONS);
    let evaluator_pass = approximately_equal(transform_samples[0], 1.0, TRANSFORM_TOLERANCE)
        && approximately_equal(transform_samples[1], 1.0, TRANSFORM_TOLERANCE)
        && approximately_equal(transform_samples[2], 0.0, TRANSFORM_TOLERANCE)
        && approximately_equal(canonical.distance(), 0.0, TRANSFORM_TOLERANCE)
        && canonical
            .objectives
            .iter()
            .zip([0.0, 0.0, 6.0])
            .all(|(actual, expected)| {
                approximately_equal(*actual, expected, CANONICAL_OBJECTIVE_TOLERANCE)
            })
        && probe_shape_residual <= MAX_SHAPE_RESIDUAL
        && directions.len() == EXPECTED_REFERENCE_DIRECTIONS;
    let mutant_pass = mutant.distance() > MIN_MUTANT_DISTANCE
        && mutant
            .objectives
            .iter()
            .zip(canonical.objectives)
            .any(|(wrong, correct)| wrong.to_bits() != correct.to_bits());

    let params = NsgaParams {
        pop: POPULATION,
        generations: GENERATIONS,
        eta_c: ETA_C,
        eta_m: ETA_M,
        p_mut: MUTATION_PROBABILITY,
        seed: INPUT_SEED,
    };

    let mut evaluations = 0usize;
    let front = {
        let mut objective = |x: &[f64]| {
            evaluations += 1;
            wfg4(x)
        };
        nsga3(&mut objective, DIMENSION, (0.0, 1.0), &directions, &params)
    };
    assert!(!front.is_empty(), "WFG4 NSGA-III front must not be empty");
    let metrics = metrics(&front, &directions);

    let mut replay_evaluations = 0usize;
    let replay = {
        let mut replay_objective = |x: &[f64]| {
            replay_evaluations += 1;
            wfg4(x)
        };
        nsga3(
            &mut replay_objective,
            DIMENSION,
            (0.0, 1.0),
            &directions,
            &params,
        )
    };

    let convergence_pass = evaluations == EXPECTED_EVALUATIONS
        && metrics.front_size >= MIN_FRONT_SIZE
        && metrics.mean_distance <= MAX_MEAN_DISTANCE
        && metrics.worst_distance <= MAX_WORST_DISTANCE
        && metrics.direction_coverage >= MIN_DIRECTION_COVERAGE
        && metrics.hypervolume >= MIN_HYPERVOLUME
        && metrics.maximum_shape_residual <= MAX_SHAPE_RESIDUAL
        && metrics.objectives_recompute_bitwise;
    let replay_mismatch = first_front_mismatch(&front, &replay);
    let replay_pass = replay_evaluations == EXPECTED_EVALUATIONS && replay_mismatch.is_none();
    let verdicts = Verdicts {
        evaluator: evaluator_pass,
        mutant: mutant_pass,
        convergence: convergence_pass,
        replay: replay_pass,
    };

    let campaign = campaign_identity(directions.len());
    let retained_front = front_identity(&campaign, &front);
    let replay_front = front_identity(&campaign, &replay);
    let result = result_identity(
        &campaign,
        &retained_front,
        &replay_front,
        canonical,
        mutant,
        &metrics,
        evaluations,
        replay_evaluations,
        transform_samples,
        probe_shape_residual,
        replay_mismatch.as_deref(),
        verdicts,
    );
    let mut emitter = Emitter::new(SUITE, CASE);
    emit_receipt(
        &mut emitter,
        &campaign,
        &result,
        &retained_front,
        &replay_front,
        canonical,
        mutant,
        &metrics,
        evaluations,
        replay_evaluations,
        directions.len(),
        transform_samples,
        probe_shape_residual,
        replay_mismatch.as_deref(),
        verdicts,
    );
    emit_benchmark(&mut emitter, "hypervolume", metrics.hypervolume);
    emit_benchmark(&mut emitter, "mean_front_distance", metrics.mean_distance);
    emit_case(
        &mut emitter,
        "wfg4-evaluator",
        evaluator_pass,
        format!(
            "campaign={}; result={}; actual s_multi(0,1,C)={transform_samples:?}; \
             canonical objectives={:?}; maximum output-reconstructed concave-sphere \
             residual={probe_shape_residual:.3e}; \
             directions={}/{EXPECTED_REFERENCE_DIRECTIONS}",
            campaign.hex(),
            result.hex(),
            canonical.objectives,
            directions.len(),
        ),
        0,
    );
    emit_case(
        &mut emitter,
        "wrong-center-mutant",
        mutant_pass,
        format!(
            "campaign={}; result={}; correct center={S_MULTI_CENTER} gives distance={:.6}; \
             wrong center={MUTANT_CENTER} gives distance={:.6}",
            campaign.hex(),
            result.hex(),
            canonical.distance(),
            mutant.distance(),
        ),
        0,
    );
    emit_case(
        &mut emitter,
        "nsga3-convergence-budget",
        convergence_pass,
        format!(
            "campaign={}; result={}; evaluations={evaluations}/{EXPECTED_EVALUATIONS}; \
             front={}; distance mean/worst={:.6}/{:.6}; coverage={}/{}={:.4}; \
             hypervolume={:.6}; output shape residual={:.3e}; objective recompute={}",
            campaign.hex(),
            result.hex(),
            metrics.front_size,
            metrics.mean_distance,
            metrics.worst_distance,
            metrics.covered_directions,
            directions.len(),
            metrics.direction_coverage,
            metrics.hypervolume,
            metrics.maximum_shape_residual,
            metrics.objectives_recompute_bitwise,
        ),
        INPUT_SEED,
    );
    emit_case(
        &mut emitter,
        "nsga3-bitwise-replay",
        replay_pass,
        format!(
            "campaign={}; result={}; front={}; replay-front={}; seed={INPUT_SEED}; \
             stream=(kernel=0x{NSGA3_STREAM_KERNEL:04x},tile={NSGA3_STREAM_TILE}); \
             original/replay evaluations={evaluations}/{replay_evaluations}; \
             first mismatch={replay_mismatch:?}",
            campaign.hex(),
            result.hex(),
            retained_front.hex(),
            replay_front.hex(),
        ),
        INPUT_SEED,
    );

    assert!(
        verdicts.pass(),
        "WFG4 battery failed: evaluator={evaluator_pass}; mutant={mutant_pass}; \
         convergence={convergence_pass} ({metrics:?}); replay={replay_pass}"
    );
}
