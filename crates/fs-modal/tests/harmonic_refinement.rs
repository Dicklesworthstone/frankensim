//! Coupled original-response residual correction at a physical antiresonance.
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_math::c64::C64;
use fs_modal::harmonic::{
    HarmonicFrequencyControls, HarmonicResponseControls, PreparedHarmonicResponse,
};

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 37,
                kernel_id: 25,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}

fn controls() -> HarmonicResponseControls {
    HarmonicResponseControls {
        maximum_coordinates: 2,
        // Existing dense-owner bound: correction scratch fits this admission.
        maximum_payload_bytes: 144 * 4 + 128 * 2,
        maximum_relative_residual: 2e-11,
    }
}

fn prepare(w: f64, cx: &Cx<'_>) -> PreparedHarmonicResponse {
    PreparedHarmonicResponse::try_new(
        &[1.0, 0.0, 0.0, 1.0],
        &[0.2, 0.0, 0.0, 0.0],
        &[2.0, -1.0, -1.0, 4.0],
        2,
        w,
        controls(),
        cx,
    )
    .unwrap()
}

fn analytic(w: f64) -> ([C64; 2], [C64; 2]) {
    // Invert the physical two-coordinate pencil analytically. Fused detuning
    // retains the authored frequency, as in the original operator assembly.
    let a = C64::new((-w).mul_add(w, 2.0), -0.2 * w);
    let b = (-w).mul_add(w, 4.0);
    let aw = C64::new(-2.0 * w, -0.2);
    let bw = -2.0 * w;
    let determinant = a.scale(b) - C64::ONE;
    let determinant_w = aw.scale(b) + a.scale(bw);
    let q = [C64::from_re(b) / determinant, C64::ONE / determinant];
    let qw = [
        (C64::from_re(bw) - determinant_w * q[0]) / determinant,
        -(determinant_w * q[1]) / determinant,
    ];
    (q, qw)
}

fn relative(actual: C64, expected: C64) {
    relative_with_bound(actual, expected, 2e-13);
}

fn weak_relative(actual: C64, expected: C64) {
    // The unchanged gate may accept a first solve without refinement. In the
    // unforced row q1=b*q2, eta=|q1-b*q2|/(|q1|+|b*q2|), so its relative
    // weak-coordinate error is bounded by 2*eta/(1-eta), plus the separately
    // checked strong-coordinate oracle allowance. Do not demand accuracy
    // beyond the caller's actual residual contract from an accepted solve.
    let limit = controls().maximum_relative_residual;
    relative_with_bound(actual, expected, 2.0 * limit / (1.0 - limit) + 2e-13);
}

fn relative_with_bound(actual: C64, expected: C64, bound: f64) {
    if expected == C64::ZERO {
        assert_eq!(actual, C64::ZERO);
    } else {
        assert!(
            (actual - expected).abs() <= bound * expected.abs(),
            "{actual:?} != {expected:?}"
        );
    }
}

#[test]
fn weak_unforced_row_meets_the_original_gate_on_both_sides_of_the_notch() {
    with_cx(|cx| {
        // The first frequency previously refused with residual 2.2151e-11
        // despite the well-conditioned physical pencil and authored 2e-11 gate.
        for w in [1.99999961853027353, 2.00000038146972647] {
            let prepared = prepare(w, cx);
            let (q, _) = analytic(w);
            for phase in [C64::ONE, C64::new(0.6, 0.8), C64::new(0.0, 1.0)] {
                for exponent in [-400, 0, 400] {
                    let drive = phase.scale(2.0_f64.powi(exponent));
                    let response = prepared.solve(&[drive, C64::ZERO], cx).unwrap();
                    assert!(response.relative_residual <= controls().maximum_relative_residual);
                    weak_relative(response.displacement[0], drive * q[0]);
                    relative(response.displacement[1], drive * q[1]);
                }
            }
        }
    });
}

#[test]
fn exact_physical_notch_and_complex_drive_phase_remain_exact() {
    with_cx(|cx| {
        let prepared = prepare(2.0, cx);
        for phase in [C64::ONE, C64::new(0.6, 0.8), C64::new(0.0, 1.0)] {
            for exponent in [-400, 0, 400] {
                let drive = phase.scale(2.0_f64.powi(exponent));
                let response = prepared.solve(&[drive, C64::ZERO], cx).unwrap();
                assert_eq!(response.displacement[0], C64::ZERO);
                relative(response.displacement[1], -drive);
                assert!(response.relative_residual <= controls().maximum_relative_residual);
            }
        }
    });
}

#[test]
fn same_factor_frequency_tangent_crosses_the_weak_and_exact_notch() {
    with_cx(|cx| {
        for w in [1.99999961853027353, 2.0, 2.00000038146972647] {
            let prepared = prepare(w, cx);
            let (q, qw) = analytic(w);
            for exponent in [-400, 0, 400] {
                let drive = C64::new(0.6, 0.8).scale(2.0_f64.powi(exponent));
                let response = prepared
                    .solve_with_frequency_derivative(
                        &[drive, C64::ZERO],
                        None,
                        None,
                        HarmonicFrequencyControls {
                            maximum_additional_payload_bytes: 1024,
                            maximum_work: 1_000_000,
                        },
                        cx,
                    )
                    .unwrap();
                assert!(
                    response.response.relative_residual <= controls().maximum_relative_residual
                );
                assert!(
                    response.relative_derivative_residual <= controls().maximum_relative_residual
                );
                weak_relative(response.response.displacement[0], drive * q[0]);
                relative(response.response.displacement[1], drive * q[1]);
                for i in 0..2 {
                    relative(response.displacement_derivative[i], drive * qw[i]);
                    relative(
                        response.velocity_derivative[i],
                        C64::new(0.0, -1.0) * drive * (q[i] + qw[i].scale(w)),
                    );
                }
            }
        }
    });
}
