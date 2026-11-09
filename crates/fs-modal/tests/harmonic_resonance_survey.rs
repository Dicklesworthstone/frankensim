//! G1/G3: observable peaks and antiresonances of the original forced system.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_math::c64::C64;
use fs_modal::harmonic::resonance::{
    HarmonicTransferSample, ResonanceControls, ResonanceError, TransferExtremumKind,
    survey_transfer_extrema, survey_transfer_extrema_with_cancel,
};
use fs_modal::harmonic::{
    HarmonicFrequencyControls, HarmonicResponseControls, PreparedHarmonicResponse,
};
use std::cell::{Cell, RefCell};
use std::collections::BTreeMap;

fn controls() -> ResonanceControls {
    ResonanceControls {
        absolute_angular_tolerance_rad_s: 1e-11,
        relative_angular_tolerance: 1e-12,
        maximum_evaluations: 512,
    }
}

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 91,
                kernel_id: 34,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}

#[test]
fn coupled_passive_masses_keep_both_resonances_and_the_exact_transfer_zero() {
    with_cx(|cx| {
        let mass = [1.0, 0.0, 0.0, 1.0];
        let stiffness = [5.0, -2.0, -2.0, 4.0];
        let damping = [0.25, 0.0, 0.0, 0.0];
        let force = [C64::new(0.6, 0.8), C64::ZERO];
        let grid = [0.0, 1.0, 1.4, 1.7, 2.0, 2.3, 2.8, 3.2];
        let result = survey_transfer_extrema(&grid, controls(), 3, cx, |omega| {
            PreparedHarmonicResponse::try_new(
                &mass,
                &damping,
                &stiffness,
                2,
                omega,
                HarmonicResponseControls {
                    maximum_coordinates: 2,
                    maximum_payload_bytes: 1 << 20,
                    maximum_relative_residual: 1e-11,
                },
                cx,
            )?
            .solve_with_frequency_derivative(
                &force,
                None,
                None,
                HarmonicFrequencyControls {
                    maximum_additional_payload_bytes: 1 << 20,
                    maximum_work: 1_000_000,
                },
                cx,
            )?
            .displacement_transfer(&[1.0, 0.0])
        })
        .unwrap();

        assert_eq!(result.extrema.len(), 3);
        assert_eq!(result.grid_points, grid.len());
        assert_eq!(result.interval_rad_s, (0.0, 3.2));
        assert!(result.evaluations <= controls().maximum_evaluations);
        // Independent original-equation oracle: H=a/D, a=4-w^2,
        // D=(5-w^2-i*c*w)*a-4. Away from a=0 its magnitude peaks satisfy
        // a^4+(1-c^2/2)*a^3+4*a-16=0. These roots were evaluated with
        // 70-digit decimal arithmetic, not a sampled response or modal fit.
        let expected = [1.560_109_172_196_526_3, 2.0, 2.559_209_928_171_815];
        let kinds = [
            TransferExtremumKind::Maximum,
            TransferExtremumKind::Minimum,
            TransferExtremumKind::Maximum,
        ];
        for ((extremum, &frequency), &kind) in result.extrema.iter().zip(&expected).zip(&kinds) {
            assert_eq!(extremum.kind, kind);
            assert!((extremum.angular_frequency_rad_s - frequency).abs() < 2e-10);
            let w = extremum.angular_frequency_rad_s;
            let a = 4.0 - w * w;
            let denominator = C64::new((5.0 - w * w) * a - 4.0, -0.25 * w * a);
            let reference = force[0] * C64::new(a, 0.0) / denominator;
            assert!((extremum.sample.amplitude - reference).abs() < 2e-11);
        }
        assert!(result.extrema[1].sample.amplitude.abs() < 2e-10);
    });
}

// A positive transfer with exactly represented stationary frequencies 1,2,3.
// The middle extremum is a nonzero notch, the other two are equal peaks.
fn polynomial(omega: f64) -> HarmonicTransferSample {
    let x = omega - 2.0;
    let a = x * x - 1.0;
    HarmonicTransferSample {
        amplitude: C64::new(5.0 - a * a, 0.0),
        angular_derivative: C64::new(-4.0 * x * a, 0.0),
    }
}

const GRID: [f64; 7] = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5];

#[test]
fn exactly_stationary_grid_samples_are_shared_and_counted_once() {
    let calls = RefCell::new(BTreeMap::<u64, usize>::new());
    let result = survey_transfer_extrema_with_cancel(
        &GRID,
        controls(),
        3,
        || false,
        |omega| {
            *calls.borrow_mut().entry(omega.to_bits()).or_default() += 1;
            Ok::<_, &'static str>(polynomial(omega))
        },
    )
    .unwrap();
    assert_eq!(result.extrema.len(), 3);
    for (extremum, expected) in result.extrema.iter().zip([1.0, 2.0, 3.0]) {
        assert_eq!(extremum.angular_frequency_rad_s, expected);
    }
    for frequency in GRID {
        assert_eq!(calls.borrow()[&frequency.to_bits()], 1);
    }
    assert_eq!(result.evaluations, calls.borrow().values().sum::<usize>());
    assert_eq!(
        result.evaluations,
        GRID.len() + result.extrema.iter().map(|e| e.evaluations).sum::<usize>()
    );
}

#[test]
fn silent_or_flat_observations_do_not_invent_modes() {
    for amplitude in [C64::ZERO, C64::new(0.6, 0.8)] {
        let result = survey_transfer_extrema_with_cancel(
            &GRID,
            controls(),
            1,
            || false,
            |_| {
                Ok::<_, &'static str>(HarmonicTransferSample {
                    amplitude,
                    angular_derivative: C64::ZERO,
                })
            },
        )
        .unwrap();
        assert!(result.extrema.is_empty());
        assert_eq!(result.evaluations, GRID.len());
    }
}

#[test]
fn global_evaluation_and_extremum_limits_refuse_instead_of_truncating() {
    let calls = Cell::new(0);
    let limited = ResonanceControls {
        maximum_evaluations: GRID.len() + 1,
        ..controls()
    };
    let error = survey_transfer_extrema_with_cancel(
        &GRID,
        limited,
        3,
        || false,
        |omega| {
            calls.set(calls.get() + 1);
            Ok::<_, &'static str>(polynomial(omega))
        },
    )
    .unwrap_err();
    assert_eq!(
        error,
        ResonanceError::EvaluationBudget {
            maximum_evaluations: limited.maximum_evaluations
        }
    );
    assert_eq!(calls.get(), limited.maximum_evaluations);

    let error = survey_transfer_extrema_with_cancel(
        &GRID,
        controls(),
        2,
        || false,
        |omega| Ok::<_, &'static str>(polynomial(omega)),
    )
    .unwrap_err();
    assert_eq!(error, ResonanceError::ExtremumBudget { maximum_extrema: 2 });
}

#[test]
fn native_source_errors_and_cancellation_publish_no_partial_survey() {
    let error = survey_transfer_extrema_with_cancel(
        &GRID,
        controls(),
        3,
        || false,
        |omega| {
            if omega == 1.5 {
                Err("material response unavailable")
            } else {
                Ok(polynomial(omega))
            }
        },
    )
    .unwrap_err();
    assert_eq!(
        error,
        ResonanceError::Evaluation("material response unavailable")
    );

    let calls = Cell::new(0);
    let cancelled = survey_transfer_extrema_with_cancel(
        &GRID,
        controls(),
        3,
        || calls.get() == 3,
        |omega| {
            calls.set(calls.get() + 1);
            Ok::<_, &'static str>(polynomial(omega))
        },
    )
    .unwrap_err();
    assert_eq!(cancelled, ResonanceError::Cancelled);
    assert_eq!(calls.get(), 3);
}

#[test]
fn malformed_frequency_grids_refuse_before_the_first_physical_evaluation() {
    for grid in [
        &[][..],
        &[1.0][..],
        &[1.0, 1.0][..],
        &[2.0, 1.0][..],
        &[-1.0, 1.0][..],
        &[0.0, f64::NAN][..],
        &[0.0, f64::INFINITY][..],
    ] {
        let calls = Cell::new(0);
        let error = survey_transfer_extrema_with_cancel(
            grid,
            controls(),
            3,
            || false,
            |omega| {
                calls.set(calls.get() + 1);
                Ok::<_, &'static str>(polynomial(omega))
            },
        )
        .unwrap_err();
        assert!(matches!(error, ResonanceError::Invalid(_)));
        assert_eq!(calls.get(), 0);
    }
}
