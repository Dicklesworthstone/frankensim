//! Gauntlet G3 relations for the live optimization-step surface.
//!
//! These checks supplement the fixed manifold and budget pins in the
//! conformance battery; they do not replace those cases.

use asupersync::types::Budget;
use fs_exec::{CancelGate, Cx, ExecMode, StreamKey};
use fs_opt::{DescentOptions, DescentStop, EvalLimit, Manifold, descend_fn};
use fs_propcheck::metamorphic::{
    RelationCase, RelationObservation, Tolerance, check_relation, unit_rescaling,
};

type QuadraticCase = (f64, f64, f64);

const STEPS: u32 = 3;
const EXPECTED_EVALS: u64 = 1 + 2 * STEPS as u64 + 1;

#[derive(Debug)]
struct DescentScaleReceipt {
    x: f64,
    f0: f64,
    f_final: f64,
    evals: u64,
    steps_taken: u32,
    stop: DescentStop,
    budget_stopped: bool,
    work_upper_bound: u64,
    workspace_upper_bound_bytes: u64,
}

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 0x2ACE_0401,
                kernel_id: 1,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx)
    })
}

fn scale_for_exponent(exponent: i64) -> f64 {
    fs_math::det::powi(
        2.0,
        i32::try_from(exponent).expect("generated exponent fits i32"),
    )
}

#[test]
fn g3_descend_fn_is_equivariant_under_power_of_two_unit_rescaling() {
    with_cx(|cx| {
        let operator = |&(x0, target, fd_h): &QuadraticCase| {
            let objective = |x: &[f64]| {
                let residual = x[0] - target;
                residual * residual
            };
            let report = descend_fn(
                Manifold::Rn { dim: 1 },
                &objective,
                &[x0],
                DescentOptions {
                    steps: STEPS,
                    lr: 0.125,
                    fd_h,
                    ..DescentOptions::default()
                },
                EvalLimit::Unlimited,
                cx,
            )
            .expect("generated quadratic descent is admitted");
            DescentScaleReceipt {
                x: report.x[0],
                f0: report.f0,
                f_final: report.f_final,
                evals: report.evals,
                steps_taken: report.steps_taken,
                stop: report.stop,
                budget_stopped: report.budget_stopped,
                work_upper_bound: report.work_upper_bound,
                workspace_upper_bound_bytes: report.workspace_upper_bound_bytes,
            }
        };
        let relation = unit_rescaling(
            "quadratic-descent-power-of-two-units",
            Tolerance::AbsoluteRelative {
                max_abs: 2.0e-12,
                max_relative: 2.0e-12,
            },
            |&(x0, target, fd_h): &QuadraticCase, &exponent: &i64| {
                let scale = scale_for_exponent(exponent);
                (x0 * scale, target * scale, fd_h * scale)
            },
            |base: &DescentScaleReceipt,
             transformed: &DescentScaleReceipt,
             &exponent: &i64,
             tolerance: Tolerance| {
                let scale = scale_for_exponent(exponent);
                let scale_sq = scale * scale;
                let x = tolerance.evaluate_scalar(base.x * scale, transformed.x);
                let f0 = tolerance.evaluate_scalar(base.f0 * scale_sq, transformed.f0);
                let f_final =
                    tolerance.evaluate_scalar(base.f_final * scale_sq, transformed.f_final);
                let discrete_receipts_match = base.evals == EXPECTED_EVALS
                    && transformed.evals == EXPECTED_EVALS
                    && base.steps_taken == STEPS
                    && transformed.steps_taken == STEPS
                    && base.stop == DescentStop::StepLimit
                    && transformed.stop == DescentStop::StepLimit
                    && !base.budget_stopped
                    && !transformed.budget_stopped
                    && base.work_upper_bound == transformed.work_upper_bound
                    && base.workspace_upper_bound_bytes == transformed.workspace_upper_bound_bytes;
                let discrete_margin = if discrete_receipts_match { 0.0 } else { -1.0 };
                RelationObservation::new(
                    x.margin()
                        .min(f0.margin())
                        .min(f_final.margin())
                        .min(discrete_margin),
                    "scaled descent preserves x, quadratic objectives, and discrete receipts",
                )
            },
        );
        let exponents = [-4_i64, -3, -2, -1, 1, 2, 3, 4];

        check_relation(
            "fs-opt::descend_fn",
            0x2ACE_0402,
            256,
            |stream| {
                let target = stream.f64_in(-4.0, 4.0);
                let magnitude = stream.f64_in(0.25, 2.0);
                let sign = if stream.next_u64().is_multiple_of(2) {
                    -1.0
                } else {
                    1.0
                };
                let exponent = exponents[(stream.next_u64() % exponents.len() as u64) as usize];
                RelationCase::new(
                    (
                        target + sign * magnitude,
                        target,
                        fs_math::det::powi(2.0, -10),
                    ),
                    exponent,
                )
            },
            &operator,
            &relation,
        );
    });
}

#[test]
fn g3_relative_landed_step_closure_is_power_of_two_rescaling_invariant() {
    with_cx(|cx| {
        let run = |scale: f64| {
            let objective = |x: &[f64]| x[0] * x[0];
            descend_fn(
                Manifold::Rn { dim: 1 },
                &objective,
                &[scale],
                DescentOptions {
                    steps: 5,
                    lr: 0.125,
                    fd_h: 1e-6 * scale,
                    closure_threshold: 0.3,
                    ..DescentOptions::default()
                },
                EvalLimit::Unlimited,
                cx,
            )
            .expect("rescaled closure case")
        };

        let base = run(1.0);
        let scaled = run(8.0);
        for report in [&base, &scaled] {
            assert_eq!(report.stop, DescentStop::ClosureThreshold);
            assert_eq!(report.steps_taken, 0);
            assert_eq!(report.evals, 3);
            assert!(!report.budget_stopped);
            assert_eq!(report.f_final.to_bits(), report.f0.to_bits());
        }
        assert_eq!(scaled.x[0].to_bits(), (base.x[0] * 8.0).to_bits());
        assert_eq!(scaled.f0.to_bits(), (base.f0 * 64.0).to_bits());
        assert_eq!(scaled.work_upper_bound, base.work_upper_bound);
        assert_eq!(
            scaled.workspace_upper_bound_bytes,
            base.workspace_upper_bound_bytes
        );
    });
}
