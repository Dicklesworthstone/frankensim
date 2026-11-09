//! G1/G3: complete physical design derivatives of every resolved transfer feature.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_math::c64::C64;
use fs_modal::harmonic::resonance::{
    HarmonicTransferParameterJet, HarmonicTransferSample, ResonanceControls, ResonanceError,
    TransferExtremumKind, survey_transfer_extrema_with_parameter_derivative,
    survey_transfer_extrema_with_parameter_derivative_with_cancel,
};
use fs_modal::harmonic::{
    HarmonicFrequencyControls, HarmonicFrequencyParameterControls, HarmonicResponseControls,
    HarmonicResponseError, PreparedHarmonicResponse,
};
use fs_modal::quadratic_modes::QuadraticParameterDirection;
use std::cell::{Cell, RefCell};
use std::collections::BTreeMap;

fn controls() -> ResonanceControls {
    ResonanceControls {
        absolute_angular_tolerance_rad_s: 1e-12,
        relative_angular_tolerance: 1e-12,
        maximum_evaluations: 512,
    }
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1.0),
        "{actual:?} != {expected:?}",
    );
}

fn near_complex(actual: C64, expected: C64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1.0),
        "{actual:?} != {expected:?}",
    );
}

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 113,
                kernel_id: 293,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}

#[derive(Clone, Copy)]
enum Direction {
    AttachmentSpring,
    PhysicalDilation,
    SourceGainPhase,
}

impl Direction {
    fn data(self) -> ([f64; 4], [f64; 4], C64) {
        match self {
            Self::AttachmentSpring => ([1.0, -1.0, -1.0, 1.0], [0.0; 4], C64::ZERO),
            Self::PhysicalDilation => (
                [10.0, -4.0, -4.0, 8.0],
                [0.25, 0.0, 0.0, 0.0],
                C64::from_re(2.0),
            ),
            Self::SourceGainPhase => ([0.0; 4], [0.0; 4], C64::new(0.2, 0.3)),
        }
    }
}

fn owner(w: f64, cx: &Cx<'_>) -> Result<PreparedHarmonicResponse, HarmonicResponseError> {
    PreparedHarmonicResponse::try_new(
        &[1.0, 0.0, 0.0, 1.0],
        &[0.25, 0.0, 0.0, 0.0],
        &[5.0, -2.0, -2.0, 4.0],
        2,
        w,
        HarmonicResponseControls {
            maximum_coordinates: 2,
            maximum_payload_bytes: 1 << 20,
            maximum_relative_residual: 2e-11,
        },
        cx,
    )
}

fn first(w: f64, cx: &Cx<'_>) -> Result<HarmonicTransferSample, HarmonicResponseError> {
    owner(w, cx)?
        .solve_with_frequency_derivative(
            &[C64::ONE, C64::ZERO],
            None,
            None,
            HarmonicFrequencyControls {
                maximum_additional_payload_bytes: 1 << 20,
                maximum_work: 1_000_000,
            },
            cx,
        )?
        .displacement_transfer(&[1.0, 0.0])
}

fn complete(
    w: f64,
    direction: Direction,
    cx: &Cx<'_>,
) -> Result<HarmonicTransferParameterJet, HarmonicResponseError> {
    let (dk, dc, ft) = direction.data();
    owner(w, cx)?
        .solve_with_frequency_parameter_derivative(
            &[C64::ONE, C64::ZERO],
            QuadraticParameterDirection {
                mass: &[0.0; 4],
                damping: &dc,
                stiffness: &dk,
            },
            None,
            None,
            Some(&[ft, C64::ZERO]),
            None,
            HarmonicFrequencyParameterControls {
                maximum_additional_payload_bytes: 1 << 20,
                maximum_work: 1_000_000,
            },
            cx,
        )?
        .displacement_transfer_parameter(&[1.0, 0.0])
}

/// Independent original two-coordinate quotient H=b/((a-i*c*w)*b-g^2).
fn physical_reference(w: f64, direction: Direction) -> (C64, C64, C64) {
    let (dk, dc, ft) = direction.data();
    let b = 4.0 - w * w;
    let a = C64::new(5.0 - w * w, -0.25 * w);
    let aw = C64::new(-2.0 * w, -0.25);
    let at = C64::new(dk[0], -dc[0] * w);
    let denominator = a.scale(b) - C64::from_re(4.0);
    let dw = aw.scale(b) + a.scale(-2.0 * w);
    let dt = at.scale(b) + a.scale(dk[3]) + C64::from_re(4.0 * dk[1]);
    let h = C64::from_re(b) / denominator;
    let hw =
        C64::from_re(-2.0 * w) / denominator - C64::from_re(b) * dw / (denominator * denominator);
    let ht = (ft.scale(b) + C64::from_re(dk[3])) / denominator
        - C64::from_re(b) * dt / (denominator * denominator);
    (h, hw, ht)
}

/// Differentiate the exact peak polynomial of the physical two-mass response.
fn attachment_peak_shift(w: f64) -> f64 {
    let y = w * w;
    let a = 5.0 - y;
    let b = 4.0 - y;
    let r = a * b - 4.0;
    let s = b * b + 4.0;
    // Q=2*(a*b-g^2)*(b^2+g^2)-c^2*b^3.
    // The spring direction changes both diagonals and g: a_t=b_t=g_t=1.
    let q_t = 2.0 * (a + b - 4.0) * s + 2.0 * r * (2.0 * b + 4.0) - 3.0 * 0.25 * 0.25 * b * b;
    let q_y = -2.0 * (a + b) * s - 4.0 * b * r + 3.0 * 0.25 * 0.25 * b * b;
    -q_t / (2.0 * w * q_y)
}

#[test]
fn native_two_mass_survey_differentiates_both_peaks_and_the_exact_cancellation_notch() {
    with_cx(|cx| {
        let grid = [0.0, 1.0, 1.75, 2.0, 2.25, 2.75, 3.2];
        // Roots of the independent original peak polynomial, evaluated at
        // 70-digit precision in the existing physical survey oracle.
        let frequencies = [1.560_109_172_196_526_3, 2.0, 2.559_209_928_171_815];
        for direction in [
            Direction::AttachmentSpring,
            Direction::PhysicalDilation,
            Direction::SourceGainPhase,
        ] {
            let first_calls = Cell::new(0);
            let notch_calls = Cell::new(0);
            let full_calls = Cell::new(0);
            let survey = survey_transfer_extrema_with_parameter_derivative(
                &grid,
                controls(),
                3,
                cx,
                |w| {
                    first_calls.set(first_calls.get() + 1);
                    if w == 2.0 {
                        notch_calls.set(notch_calls.get() + 1);
                    }
                    first(w, cx)
                },
                |w| {
                    full_calls.set(full_calls.get() + 1);
                    complete(w, direction, cx)
                },
            )
            .unwrap();
            assert_eq!(survey.extrema.len(), 3);
            assert_eq!(survey.grid_points, grid.len());
            assert_eq!(survey.interval_rad_s, (0.0, 3.2));
            assert_eq!(full_calls.get(), 3);
            assert_eq!(notch_calls.get(), 1);
            assert_eq!(survey.evaluations, first_calls.get() + full_calls.get());
            assert_eq!(
                survey.evaluations,
                grid.len() + survey.extrema.iter().map(|r| r.evaluations).sum::<usize>()
            );

            for (i, response) in survey.extrema.iter().enumerate() {
                let w = response.extremum.angular_frequency_rad_s;
                let kind = if i == 1 {
                    TransferExtremumKind::Minimum
                } else {
                    TransferExtremumKind::Maximum
                };
                assert_eq!(response.extremum.kind, kind);
                near(w, frequencies[i], 3e-11);
                let expected_shift = match direction {
                    Direction::AttachmentSpring if i == 1 => 0.25,
                    Direction::AttachmentSpring => attachment_peak_shift(w),
                    // K(theta)=(1+theta)^2*K, C(theta)=(1+theta)*C,
                    // F(theta)=(1+theta)^2*F => H(w,theta)=H(w/(1+theta),0).
                    Direction::PhysicalDilation => w,
                    Direction::SourceGainPhase => 0.0,
                };
                near(
                    response.angular_frequency_derivative_rad_s,
                    expected_shift,
                    3e-9,
                );
                let (h, hw, ht) = physical_reference(w, direction);
                near_complex(response.extremum.sample.amplitude, h, 2e-10);
                near_complex(
                    response.amplitude_derivative,
                    ht + hw.scale(expected_shift),
                    2e-8,
                );
                if i == 1 {
                    assert_eq!(w, 2.0);
                    assert_eq!(response.extremum.sample.amplitude, C64::ZERO);
                    assert!(response.normalized_curvature > 0.0);
                    assert!(response.log_amplitude_derivative.is_none());
                    assert!(response.phase_derivative.is_none());
                } else {
                    assert!(response.normalized_curvature < 0.0);
                    let expected_log_phase = (ht + hw.scale(expected_shift)) / h;
                    near(
                        response.log_amplitude_derivative.unwrap(),
                        expected_log_phase.re,
                        2e-8,
                    );
                    near(
                        response.phase_derivative.unwrap(),
                        expected_log_phase.im,
                        2e-8,
                    );
                }
            }
        }
    });
}

const GRID: [f64; 7] = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5];

/// Two exact magnitude maxima and one simple exact zero, with a moving
/// coherent source/receiver law H(w,t)=exp(alpha*t)*S*h(w-2-t/4).
fn manufactured(w: f64, units: f64) -> HarmonicTransferParameterJet {
    let x = w - 2.0;
    let scale = C64::new(0.75, 0.5).scale(units);
    let alpha = C64::new(0.5, 0.25);
    let h = scale.scale(x * (3.0 - x * x));
    let hw = scale.scale(3.0 - 3.0 * x * x);
    let hww = scale.scale(-6.0 * x);
    HarmonicTransferParameterJet {
        sample: HarmonicTransferSample {
            amplitude: h,
            angular_derivative: hw,
        },
        angular_second_derivative: hww,
        parameter_derivative: alpha * h - hw.scale(0.25),
        angular_parameter_derivative: alpha * hw - hww.scale(0.25),
    }
}

#[test]
fn exact_notches_and_complex_gain_phase_survive_extreme_units_with_one_grid_cache() {
    for exponent in [-700, 0, 700] {
        let units = 2.0_f64.powi(exponent);
        let first_calls = RefCell::new(BTreeMap::<u64, usize>::new());
        let full_calls = RefCell::new(BTreeMap::<u64, usize>::new());
        let survey = survey_transfer_extrema_with_parameter_derivative_with_cancel(
            &GRID,
            controls(),
            3,
            || false,
            |w| {
                *first_calls.borrow_mut().entry(w.to_bits()).or_default() += 1;
                Ok::<_, ()>(manufactured(w, units).sample)
            },
            |w| {
                *full_calls.borrow_mut().entry(w.to_bits()).or_default() += 1;
                Ok::<_, ()>(manufactured(w, units))
            },
        )
        .unwrap();
        assert_eq!(survey.extrema.len(), 3);
        assert_eq!(survey.grid_points, GRID.len());
        assert_eq!(survey.interval_rad_s, (0.5, 3.5));
        assert_eq!(full_calls.borrow().len(), 3);
        assert!(first_calls.borrow().values().all(|&count| count == 1));
        assert!(full_calls.borrow().values().all(|&count| count == 1));
        for w in GRID {
            assert_eq!(first_calls.borrow()[&w.to_bits()], 1);
        }
        assert_eq!(
            survey.evaluations,
            first_calls.borrow().values().sum::<usize>() + 3
        );
        assert_eq!(
            survey.evaluations,
            GRID.len() + survey.extrema.iter().map(|r| r.evaluations).sum::<usize>()
        );
        for (i, response) in survey.extrema.iter().enumerate() {
            let w = (i + 1) as f64;
            assert_eq!(response.extremum.angular_frequency_rad_s, w);
            assert_eq!(response.evaluations, response.extremum.evaluations + 1);
            near(response.angular_frequency_derivative_rad_s, 0.25, 2e-13);
            let h = manufactured(w, 1.0).sample.amplitude;
            near_complex(
                response.amplitude_derivative.scale(1.0 / units),
                C64::new(0.5, 0.25) * h,
                2e-13,
            );
            if i == 1 {
                assert_eq!(response.extremum.kind, TransferExtremumKind::Minimum);
                assert_eq!(response.extremum.sample.amplitude, C64::ZERO);
                assert_eq!(response.amplitude_derivative, C64::ZERO);
                assert_eq!(response.log_amplitude_derivative, None);
                assert_eq!(response.phase_derivative, None);
            } else {
                assert_eq!(response.extremum.kind, TransferExtremumKind::Maximum);
                near(response.log_amplitude_derivative.unwrap(), 0.5, 2e-13);
                near(response.phase_derivative.unwrap(), 0.25, 2e-13);
            }
        }
    }
}

#[test]
fn flat_silent_and_monotone_transfers_use_no_final_parameter_jets() {
    for case in 0..3 {
        let first_calls = Cell::new(0);
        let survey = survey_transfer_extrema_with_parameter_derivative_with_cancel(
            &GRID,
            ResonanceControls {
                maximum_evaluations: GRID.len(),
                ..controls()
            },
            1,
            || false,
            |w| {
                first_calls.set(first_calls.get() + 1);
                Ok::<_, ()>(match case {
                    0 => HarmonicTransferSample {
                        amplitude: C64::ZERO,
                        angular_derivative: C64::ZERO,
                    },
                    1 => HarmonicTransferSample {
                        amplitude: C64::new(0.75, 0.5),
                        angular_derivative: C64::ZERO,
                    },
                    _ => HarmonicTransferSample {
                        amplitude: C64::from_re(1.0 + w),
                        angular_derivative: C64::ONE,
                    },
                })
            },
            |_| -> Result<HarmonicTransferParameterJet, ()> { unreachable!() },
        )
        .unwrap();
        assert!(survey.extrema.is_empty());
        assert_eq!(first_calls.get(), GRID.len());
        assert_eq!(survey.evaluations, GRID.len());
    }
}

#[test]
fn all_candidate_jets_are_reserved_before_refinement_with_one_total_cap() {
    let first_calls = Cell::new(0);
    let full_calls = Cell::new(0);
    let cap = GRID.len() + 2;
    let error = survey_transfer_extrema_with_parameter_derivative_with_cancel(
        &GRID,
        ResonanceControls {
            maximum_evaluations: cap,
            ..controls()
        },
        3,
        || false,
        |w| {
            first_calls.set(first_calls.get() + 1);
            Ok::<_, ()>(manufactured(w, 1.0).sample)
        },
        |w| {
            full_calls.set(full_calls.get() + 1);
            Ok::<_, ()>(manufactured(w, 1.0))
        },
    )
    .unwrap_err();
    assert_eq!(
        error,
        ResonanceError::EvaluationBudget {
            maximum_evaluations: cap
        }
    );
    assert_eq!(first_calls.get(), GRID.len());
    assert_eq!(full_calls.get(), 0);

    let full_calls = Cell::new(0);
    let error = survey_transfer_extrema_with_parameter_derivative_with_cancel(
        &GRID,
        controls(),
        2,
        || false,
        |w| Ok::<_, ()>(manufactured(w, 1.0).sample),
        |w| {
            full_calls.set(full_calls.get() + 1);
            Ok::<_, ()>(manufactured(w, 1.0))
        },
    )
    .unwrap_err();
    assert_eq!(error, ResonanceError::ExtremumBudget { maximum_extrema: 2 });
    assert_eq!(full_calls.get(), 0);

    let successful = survey_transfer_extrema_with_parameter_derivative_with_cancel(
        &GRID,
        controls(),
        3,
        || false,
        |w| Ok::<_, ()>(manufactured(w, 1.0).sample),
        |w| Ok::<_, ()>(manufactured(w, 1.0)),
    )
    .unwrap();
    for cap in [successful.evaluations, successful.evaluations - 1] {
        let first_calls = Cell::new(0);
        let full_calls = Cell::new(0);
        let result = survey_transfer_extrema_with_parameter_derivative_with_cancel(
            &GRID,
            ResonanceControls {
                maximum_evaluations: cap,
                ..controls()
            },
            3,
            || false,
            |w| {
                first_calls.set(first_calls.get() + 1);
                Ok::<_, ()>(manufactured(w, 1.0).sample)
            },
            |w| {
                full_calls.set(full_calls.get() + 1);
                Ok::<_, ()>(manufactured(w, 1.0))
            },
        );
        if cap == successful.evaluations {
            assert_eq!(result.unwrap().evaluations, cap);
            assert_eq!(full_calls.get(), 3);
        } else {
            assert_eq!(
                result.unwrap_err(),
                ResonanceError::EvaluationBudget {
                    maximum_evaluations: cap
                }
            );
            assert!(full_calls.get() < 3);
        }
        assert!(first_calls.get() <= cap - 3);
        assert!(first_calls.get() + full_calls.get() <= cap);
    }
}

#[test]
fn inconsistent_exact_zero_or_curvature_and_late_failures_publish_no_partial_survey() {
    for case in 0..4 {
        let result = survey_transfer_extrema_with_parameter_derivative_with_cancel(
            &GRID,
            controls(),
            3,
            || false,
            |w| Ok::<_, ()>(manufactured(w, 1.0).sample),
            |w| {
                let mut jet = manufactured(w, 1.0);
                match case {
                    0 if w == 2.0 => jet.sample.amplitude = C64::from_re(f64::from_bits(1)),
                    1 if w == 1.0 => {
                        jet.sample.angular_derivative = C64::from_re(f64::from_bits(1))
                    }
                    2 => jet.angular_second_derivative = -jet.angular_second_derivative,
                    3 => jet.angular_parameter_derivative = C64::new(f64::NAN, 0.0),
                    _ => {}
                }
                Ok::<_, ()>(jet)
            },
        );
        assert!(matches!(result, Err(ResonanceError::Invalid(_))));
    }

    let calls = Cell::new(0);
    assert_eq!(
        survey_transfer_extrema_with_parameter_derivative_with_cancel(
            &GRID,
            controls(),
            3,
            || false,
            |w| Ok::<_, &'static str>(manufactured(w, 1.0).sample),
            |w| {
                calls.set(calls.get() + 1);
                if calls.get() == 2 {
                    Err("original notch constitutive refusal")
                } else {
                    Ok(manufactured(w, 1.0))
                }
            },
        )
        .unwrap_err(),
        ResonanceError::Evaluation("original notch constitutive refusal"),
    );
    assert_eq!(calls.get(), 2);

    let calls = Cell::new(0);
    let cancelled = Cell::new(false);
    assert_eq!(
        survey_transfer_extrema_with_parameter_derivative_with_cancel(
            &GRID,
            controls(),
            3,
            || cancelled.get(),
            |w| Ok::<_, ()>(manufactured(w, 1.0).sample),
            |w| {
                calls.set(calls.get() + 1);
                if calls.get() == 2 {
                    cancelled.set(true);
                }
                Ok::<_, ()>(manufactured(w, 1.0))
            },
        )
        .unwrap_err(),
        ResonanceError::Cancelled,
    );
    assert_eq!(calls.get(), 2);

    assert_eq!(
        survey_transfer_extrema_with_parameter_derivative_with_cancel(
            &GRID,
            controls(),
            3,
            || false,
            |w| if w == 1.5 {
                Err("original sampled pole")
            } else {
                Ok(manufactured(w, 1.0).sample)
            },
            |_| -> Result<HarmonicTransferParameterJet, &'static str> { unreachable!() },
        )
        .unwrap_err(),
        ResonanceError::Evaluation("original sampled pole"),
    );
    assert_eq!(
        survey_transfer_extrema_with_parameter_derivative_with_cancel(
            &GRID,
            controls(),
            3,
            || true,
            |_| -> Result<HarmonicTransferSample, ()> { unreachable!() },
            |_| -> Result<HarmonicTransferParameterJet, ()> { unreachable!() },
        )
        .unwrap_err(),
        ResonanceError::Cancelled,
    );
}
