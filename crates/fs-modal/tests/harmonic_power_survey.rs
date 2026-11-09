//! G1/G3/G4: physical signed-power surveys, sensitivities and bounded refusal.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_math::c64::C64;
use fs_modal::harmonic::resonance::{
    HarmonicPowerCurvatureSample, HarmonicPowerParameterJet, HarmonicPowerSample,
    PowerExtremumKind, ResonanceControls, ResonanceError, survey_power_extrema,
    survey_power_extrema_with_cancel, survey_power_extrema_with_parameter_derivative,
    survey_power_extrema_with_parameter_derivative_with_cancel,
};
use fs_modal::harmonic::{
    HarmonicFrequencyControls, HarmonicFrequencyCurvatureControls,
    HarmonicFrequencyParameterControls, HarmonicPortControls,
    HarmonicPortFrequencyParameterDirection, HarmonicPortParameterDirection,
    HarmonicResponseControls, HarmonicResponseError, PreparedHarmonicResponse,
};
use fs_modal::quadratic_modes::QuadraticParameterDirection;
use std::cell::{Cell, RefCell};
use std::collections::BTreeSet;

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 103,
                kernel_id: 313,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}

fn controls() -> ResonanceControls {
    ResonanceControls {
        absolute_angular_tolerance_rad_s: 1e-11,
        relative_angular_tolerance: 1e-11,
        maximum_evaluations: 384,
    }
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1.0),
        "{actual:?} != {expected:?}"
    );
}

const MASS: [f64; 4] = [1.0, 0.0, 0.0, 1.0];
const DAMPING: [f64; 4] = [0.04, 0.0, 0.0, 0.04];
const STIFFNESS: [f64; 4] = [2.0, -1.0, -1.0, 2.0];
const ZERO: [f64; 4] = [0.0; 4];
const ZERO_Z: [C64; 4] = [C64::ZERO; 4];
const LOAD: [C64; 4] = [
    C64::new(0.08, 0.0),
    C64::ZERO,
    C64::ZERO,
    C64::new(0.08, 0.0),
];
const FORCE: [C64; 2] = [C64::ONE, C64::ZERO];

fn owner(w: f64, cx: &Cx<'_>) -> Result<PreparedHarmonicResponse, HarmonicResponseError> {
    PreparedHarmonicResponse::try_with_port_impedance(
        &MASS,
        &DAMPING,
        &STIFFNESS,
        2,
        &MASS,
        2,
        &LOAD,
        w,
        HarmonicResponseControls {
            maximum_coordinates: 2,
            maximum_payload_bytes: 1 << 20,
            maximum_relative_residual: 2e-11,
        },
        HarmonicPortControls {
            maximum_ports: 2,
            maximum_projection_terms: 100_000,
            maximum_payload_bytes: 1 << 20,
        },
        cx,
    )
}

fn physical_sample(w: f64, cx: &Cx<'_>) -> Result<HarmonicPowerSample, HarmonicResponseError> {
    owner(w, cx)?
        .solve_with_frequency_derivative(
            &FORCE,
            None,
            Some(&ZERO_Z),
            HarmonicFrequencyControls {
                maximum_additional_payload_bytes: 1 << 20,
                maximum_work: 1_000_000,
            },
            cx,
        )?
        .port_power_sample()
}

fn physical_curvature(
    w: f64,
    cx: &Cx<'_>,
) -> Result<HarmonicPowerCurvatureSample, HarmonicResponseError> {
    owner(w, cx)?
        .solve_with_frequency_curvature(
            &FORCE,
            None,
            None,
            Some(&ZERO_Z),
            Some(&ZERO_Z),
            HarmonicFrequencyCurvatureControls {
                maximum_additional_payload_bytes: 1 << 20,
                maximum_work: 1_000_000,
            },
            cx,
        )?
        .port_power_curvature_sample()
}

fn physical_jet(w: f64, cx: &Cx<'_>) -> Result<HarmonicPowerParameterJet, HarmonicResponseError> {
    // K(theta)=exp(2theta)K, C(theta)=exp(theta)C, Z(theta)=exp(theta)Z,
    // F(theta)=exp((1/2+0.7i)theta)F. The actual absorbed-power curve is
    // exactly P(w,theta)=P(w*exp(-theta),0), retaining source phase.
    owner(w, cx)?
        .solve_with_frequency_parameter_and_port_derivative(
            &FORCE,
            QuadraticParameterDirection {
                mass: &ZERO,
                damping: &DAMPING,
                stiffness: &[4.0, -2.0, -2.0, 4.0],
            },
            None,
            None,
            Some(&[C64::new(0.5, 0.7), C64::ZERO]),
            None,
            &ZERO_Z,
            &ZERO_Z,
            HarmonicPortFrequencyParameterDirection {
                parameter: HarmonicPortParameterDirection {
                    force_map: &ZERO,
                    impedance: &LOAD,
                },
                impedance_angular_parameter_derivative: &ZERO_Z,
            },
            HarmonicFrequencyParameterControls {
                maximum_additional_payload_bytes: 1 << 20,
                maximum_work: 1_000_000,
            },
            cx,
        )?
        .port_power_parameter_jet()
}

/// Independent original normal-coordinate expression for the two coupled masses.
/// K has eigenvalues 1 and 3; each normalized mode receives squared force 1/2.
fn physical_reference(w: f64) -> [f64; 3] {
    let s = w * w;
    let c = 0.12;
    let mut value = 0.0;
    let mut slope = 0.0;
    let mut curvature = 0.0;
    for lambda in [1.0, 3.0] {
        let denominator = (lambda - s) * (lambda - s) + c * c * s;
        let numerator = lambda * lambda - s * s;
        value += s / denominator;
        slope += numerator / (denominator * denominator);
        curvature += -2.0 * s / (denominator * denominator)
            - 2.0 * numerator * (2.0 * s + c * c - 2.0 * lambda)
                / (denominator * denominator * denominator);
    }
    [
        0.02 * value,
        0.04 * w * slope,
        0.04 * (slope + 2.0 * s * curvature),
    ]
}

fn reference_stationary(mut a: f64, mut b: f64) -> f64 {
    let positive = physical_reference(a)[1] > 0.0;
    assert_ne!(positive, physical_reference(b)[1] > 0.0);
    for _ in 0..100 {
        let middle = a + 0.5 * (b - a);
        if middle == a || middle == b {
            break;
        }
        if (physical_reference(middle)[1] > 0.0) == positive {
            a = middle;
        } else {
            b = middle;
        }
    }
    a + 0.5 * (b - a)
}

#[test]
fn coupled_absorbed_power_survey_resolves_two_resonances_and_the_intervening_minimum() {
    with_cx(|cx| {
        let grid = [0.3, 0.7, 1.05, 1.25, 1.45, 1.6, 1.8, 2.2];
        let expected = [
            reference_stationary(0.9, 1.1),
            reference_stationary(1.25, 1.6),
            reference_stationary(1.65, 1.85),
        ];
        let first_calls = RefCell::new(BTreeSet::new());
        let final_calls = Cell::new(0);
        let result = survey_power_extrema(
            &grid,
            controls(),
            3,
            cx,
            |w| {
                assert!(
                    first_calls.borrow_mut().insert(w.to_bits()),
                    "repeated original solve at {w}"
                );
                physical_sample(w, cx)
            },
            |w| {
                final_calls.set(final_calls.get() + 1);
                physical_curvature(w, cx)
            },
        )
        .unwrap();
        assert_eq!(result.extrema.len(), 3);
        assert_eq!(result.grid_points, grid.len());
        assert_eq!(result.interval_rad_s, (grid[0], grid[grid.len() - 1]));
        assert_eq!(final_calls.get(), 3);
        assert_eq!(
            result.evaluations,
            first_calls.borrow().len() + final_calls.get()
        );
        assert_eq!(
            result.evaluations,
            grid.len() + result.extrema.iter().map(|e| e.evaluations).sum::<usize>()
        );
        for (index, (extremum, &frequency)) in result.extrema.iter().zip(&expected).enumerate() {
            assert_eq!(
                extremum.kind,
                if index == 1 {
                    PowerExtremumKind::Minimum
                } else {
                    PowerExtremumKind::Maximum
                }
            );
            near(extremum.angular_frequency_rad_s, frequency, 3e-11);
            let reference = physical_reference(frequency);
            near(extremum.sample.average_power_w, reference[0], 3e-11);
            near(extremum.angular_second_derivative, reference[2], 2e-8);
        }

        let first_calls = Cell::new(0);
        let final_calls = Cell::new(0);
        let sensitivity = survey_power_extrema_with_parameter_derivative(
            &grid,
            controls(),
            3,
            cx,
            |w| {
                first_calls.set(first_calls.get() + 1);
                physical_sample(w, cx)
            },
            |w| {
                final_calls.set(final_calls.get() + 1);
                physical_jet(w, cx)
            },
        )
        .unwrap();
        assert_eq!(sensitivity.extrema.len(), 3);
        assert_eq!(final_calls.get(), 3);
        assert_eq!(
            sensitivity.evaluations,
            first_calls.get() + final_calls.get()
        );
        for (response, &frequency) in sensitivity.extrema.iter().zip(&expected) {
            near(response.extremum.angular_frequency_rad_s, frequency, 3e-11);
            near(response.angular_frequency_derivative_rad_s, frequency, 2e-8);
            near(response.average_power_derivative, 0.0, 2e-8);
        }
    });
}

const GRID: [f64; 9] = [1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0];

fn polynomial(w: f64, units: f64, offset: f64, orientation: f64) -> HarmonicPowerParameterJet {
    let x = w - 2.0;
    let y = x * x - 0.25;
    let power = units * (offset + orientation * y * y);
    let slope = units * orientation * 4.0 * x * y;
    let curvature = units * orientation * (12.0 * x * x - 1.0);
    HarmonicPowerParameterJet {
        sample: HarmonicPowerSample {
            average_power_w: power,
            angular_derivative: slope,
        },
        angular_second_derivative: curvature,
        // Complete direction combining power gain, additive signed work and
        // a translation of the physical frequency response.
        parameter_derivative: 0.4 * power + 0.7 * units - 0.3 * slope,
        angular_parameter_derivative: 0.4 * slope - 0.3 * curvature,
    }
}

fn curvature(jet: HarmonicPowerParameterJet) -> HarmonicPowerCurvatureSample {
    HarmonicPowerCurvatureSample {
        sample: jet.sample,
        angular_second_derivative: jet.angular_second_derivative,
    }
}

#[test]
fn signed_zero_extrema_and_exact_grid_samples_keep_complete_parameter_laws_across_power_units() {
    for exponent in [-700, 0, 700] {
        let units = 2.0_f64.powi(exponent);
        for offset in [-3.0, 0.0] {
            for orientation in [-1.0, 1.0] {
                let first_calls = RefCell::new(BTreeSet::new());
                let final_calls = Cell::new(0);
                let result = survey_power_extrema_with_parameter_derivative_with_cancel(
                    &GRID,
                    controls(),
                    3,
                    || false,
                    |w| {
                        assert!(first_calls.borrow_mut().insert(w.to_bits()));
                        Ok::<_, ()>(polynomial(w, units, offset, orientation).sample)
                    },
                    |w| {
                        final_calls.set(final_calls.get() + 1);
                        Ok::<_, ()>(polynomial(w, units, offset, orientation))
                    },
                )
                .unwrap();
                assert_eq!(result.extrema.len(), 3);
                assert_eq!(result.evaluations, first_calls.borrow().len() + 3);
                assert_eq!(final_calls.get(), 3);
                for (response, frequency) in result.extrema.iter().zip([1.5, 2.0, 2.5]) {
                    let expected = polynomial(frequency, units, offset, orientation);
                    assert_eq!(response.extremum.angular_frequency_rad_s, frequency);
                    assert_eq!(
                        response.extremum.kind,
                        if expected.angular_second_derivative > 0.0 {
                            PowerExtremumKind::Minimum
                        } else {
                            PowerExtremumKind::Maximum
                        }
                    );
                    assert_eq!(
                        response.extremum.sample.average_power_w / units,
                        expected.sample.average_power_w / units
                    );
                    assert_eq!(response.extremum.sample.angular_derivative, 0.0);
                    near(response.angular_frequency_derivative_rad_s, 0.3, 1e-14);
                    near(
                        response.average_power_derivative / units,
                        0.4 * expected.sample.average_power_w / units + 0.7,
                        1e-14,
                    );
                    assert!(first_calls.borrow().contains(&frequency.to_bits()));
                }
            }
        }
    }
}

fn coarse(maximum_evaluations: usize) -> ResonanceControls {
    ResonanceControls {
        absolute_angular_tolerance_rad_s: 1.0,
        relative_angular_tolerance: 0.0,
        maximum_evaluations,
    }
}

#[test]
fn one_global_budget_reserves_every_final_jet_before_refinement_and_refuses_truncation() {
    // All three roots are cached grid midpoints within the authored frequency
    // resolution. Nine first samples plus three final curvatures are exact.
    for cap in [8, 11, 12] {
        let first_calls = Cell::new(0);
        let final_calls = Cell::new(0);
        let result = survey_power_extrema_with_cancel(
            &GRID,
            coarse(cap),
            3,
            || false,
            |w| {
                first_calls.set(first_calls.get() + 1);
                Ok::<_, ()>(polynomial(w, 1.0, 0.0, 1.0).sample)
            },
            |w| {
                final_calls.set(final_calls.get() + 1);
                Ok::<_, ()>(curvature(polynomial(w, 1.0, 0.0, 1.0)))
            },
        );
        if cap == 12 {
            let result = result.unwrap();
            assert_eq!(first_calls.get(), GRID.len());
            assert_eq!(final_calls.get(), 3);
            assert_eq!(result.evaluations, cap);
            assert!(result.extrema.iter().all(|value| value.evaluations == 1));
        } else {
            assert_eq!(
                result.unwrap_err(),
                ResonanceError::EvaluationBudget {
                    maximum_evaluations: cap
                }
            );
            assert_eq!(
                first_calls.get(),
                if cap < GRID.len() { 0 } else { GRID.len() }
            );
            assert_eq!(final_calls.get(), 0);
        }
    }

    let first_calls = Cell::new(0);
    let final_calls = Cell::new(0);
    let result = survey_power_extrema_with_parameter_derivative_with_cancel(
        &GRID,
        coarse(12),
        2,
        || false,
        |w| {
            first_calls.set(first_calls.get() + 1);
            Ok::<_, ()>(polynomial(w, 1.0, 0.0, 1.0).sample)
        },
        |w| {
            final_calls.set(final_calls.get() + 1);
            Ok::<_, ()>(polynomial(w, 1.0, 0.0, 1.0))
        },
    );
    assert_eq!(
        result.unwrap_err(),
        ResonanceError::ExtremumBudget { maximum_extrema: 2 }
    );
    assert_eq!(first_calls.get(), GRID.len());
    assert_eq!(final_calls.get(), 0);
}

#[test]
fn empty_and_grid_unresolved_power_bands_are_not_fabricated_extrema_and_invalid_grids_refuse() {
    for power in [-3.0, 0.0, 2.0] {
        let result = survey_power_extrema_with_cancel(
            &[0.0, 1.0, 2.0],
            coarse(3),
            1,
            || false,
            |_| {
                Ok::<_, ()>(HarmonicPowerSample {
                    average_power_w: power,
                    angular_derivative: 0.0,
                })
            },
            |_| unreachable!(),
        )
        .unwrap();
        assert!(result.extrema.is_empty());
        assert_eq!(result.evaluations, 3);
    }
    // A stationary band endpoint has only one-sided evidence.
    let endpoint = survey_power_extrema_with_cancel(
        &[0.0, 1.0, 2.0],
        coarse(3),
        1,
        || false,
        |w| {
            Ok::<_, ()>(HarmonicPowerSample {
                average_power_w: w * w,
                angular_derivative: 2.0 * w,
            })
        },
        |_| unreachable!(),
    )
    .unwrap();
    assert!(endpoint.extrema.is_empty());
    // Equal-sign slopes can hide a minimum/maximum pair between grid points.
    // The survey reports exactly the authored resolution, without pretending
    // that an empty list proves absence of physical extrema inside the band.
    let unresolved = survey_power_extrema_with_cancel(
        &[1.25, 2.25],
        coarse(3),
        3,
        || false,
        |w| Ok::<_, ()>(polynomial(w, 1.0, 0.0, 1.0).sample),
        |_| unreachable!(),
    )
    .unwrap();
    assert!(unresolved.extrema.is_empty());
    assert_eq!(unresolved.evaluations, 2);

    let invalid: [&[f64]; 7] = [
        &[],
        &[1.0],
        &[1.0, 1.0],
        &[2.0, 1.0],
        &[-1.0, 1.0],
        &[0.0, f64::NAN],
        &[0.0, f64::INFINITY],
    ];
    for grid in invalid {
        let calls = Cell::new(0);
        let result = survey_power_extrema_with_cancel(
            grid,
            controls(),
            3,
            || false,
            |w| {
                calls.set(calls.get() + 1);
                Ok::<_, ()>(polynomial(w, 1.0, 0.0, 1.0).sample)
            },
            |_| unreachable!(),
        );
        assert!(matches!(result, Err(ResonanceError::Invalid(_))));
        assert_eq!(calls.get(), 0);
    }
    assert!(matches!(
        survey_power_extrema_with_cancel(
            &GRID,
            controls(),
            0,
            || false,
            |_| Ok::<_, ()>(polynomial(1.0, 1.0, 0.0, 1.0).sample),
            |_| unreachable!(),
        ),
        Err(ResonanceError::Invalid(_))
    ));
    assert!(matches!(
        survey_power_extrema_with_cancel(
            &GRID,
            controls(),
            3,
            || false,
            |_| Ok::<_, ()>(HarmonicPowerSample {
                average_power_w: f64::NAN,
                angular_derivative: 1.0
            }),
            |_| unreachable!(),
        ),
        Err(ResonanceError::Invalid(_))
    ));
}

#[test]
fn cancellation_and_original_errors_abort_the_complete_power_survey_at_each_phase() {
    let calls = Cell::new(0);
    let result = survey_power_extrema_with_cancel(
        &GRID,
        controls(),
        3,
        || true,
        |w| {
            calls.set(calls.get() + 1);
            Ok::<_, ()>(polynomial(w, 1.0, 0.0, 1.0).sample)
        },
        |_| unreachable!(),
    );
    assert_eq!(result.unwrap_err(), ResonanceError::Cancelled);
    assert_eq!(calls.get(), 0);

    let cancelled = Cell::new(false);
    let final_calls = Cell::new(0);
    let result = survey_power_extrema_with_cancel(
        &GRID,
        controls(),
        3,
        || cancelled.get(),
        |w| {
            if w == 1.5 {
                cancelled.set(true);
            }
            Ok::<_, ()>(polynomial(w, 1.0, 0.0, 1.0).sample)
        },
        |w| {
            final_calls.set(final_calls.get() + 1);
            Ok::<_, ()>(curvature(polynomial(w, 1.0, 0.0, 1.0)))
        },
    );
    assert_eq!(result.unwrap_err(), ResonanceError::Cancelled);
    assert_eq!(final_calls.get(), 0);

    let cancelled = Cell::new(false);
    let final_calls = Cell::new(0);
    let result = survey_power_extrema_with_parameter_derivative_with_cancel(
        &GRID,
        coarse(12),
        3,
        || cancelled.get(),
        |w| Ok::<_, ()>(polynomial(w, 1.0, 0.0, 1.0).sample),
        |w| {
            final_calls.set(final_calls.get() + 1);
            cancelled.set(true);
            let mut jet = polynomial(w, 1.0, 0.0, 1.0);
            jet.angular_second_derivative = f64::NAN;
            Ok::<_, ()>(jet)
        },
    );
    // Poll immediately after the physical callback, before interpreting its jet.
    assert_eq!(result.unwrap_err(), ResonanceError::Cancelled);
    assert_eq!(final_calls.get(), 1);

    let final_calls = Cell::new(0);
    let result = survey_power_extrema_with_cancel(
        &GRID,
        coarse(12),
        3,
        || false,
        |w| Ok::<_, &'static str>(polynomial(w, 1.0, 0.0, 1.0).sample),
        |w| {
            final_calls.set(final_calls.get() + 1);
            if final_calls.get() == 2 {
                Err("second original curvature refused")
            } else {
                Ok(curvature(polynomial(w, 1.0, 0.0, 1.0)))
            }
        },
    );
    assert_eq!(
        result.unwrap_err(),
        ResonanceError::Evaluation("second original curvature refused")
    );
    assert_eq!(final_calls.get(), 2);

    let final_calls = Cell::new(0);
    let result = survey_power_extrema_with_cancel(
        &GRID,
        controls(),
        3,
        || false,
        |w| {
            if GRID.contains(&w) {
                Ok(polynomial(w, 1.0, 0.0, 1.0).sample)
            } else {
                Err("original refinement solve refused")
            }
        },
        |w| {
            final_calls.set(final_calls.get() + 1);
            Ok(curvature(polynomial(w, 1.0, 0.0, 1.0)))
        },
    );
    assert_eq!(
        result.unwrap_err(),
        ResonanceError::Evaluation("original refinement solve refused")
    );
    assert_eq!(final_calls.get(), 0);
}

#[test]
fn every_power_survey_result_requires_consistent_complete_oriented_original_jets() {
    for case in 0..6 {
        let result = survey_power_extrema_with_parameter_derivative_with_cancel(
            &GRID,
            coarse(12),
            3,
            || false,
            |w| Ok::<_, ()>(polynomial(w, 1.0, 0.0, 1.0).sample),
            |w| {
                let mut jet = polynomial(w, 1.0, 0.0, 1.0);
                match case {
                    0 => jet.angular_second_derivative = 0.0,
                    1 => jet.angular_second_derivative = -jet.angular_second_derivative,
                    2 => jet.sample.average_power_w = f64::from_bits(1),
                    3 => jet.sample.angular_derivative = f64::from_bits(1),
                    4 => jet.parameter_derivative = f64::NAN,
                    _ => jet.angular_parameter_derivative = f64::INFINITY,
                }
                Ok::<_, ()>(jet)
            },
        );
        assert!(matches!(result, Err(ResonanceError::Invalid(_))));
    }

    let flat = |w: f64| {
        let d = w - 2.0;
        HarmonicPowerCurvatureSample {
            sample: HarmonicPowerSample {
                average_power_w: d * d * d * d,
                angular_derivative: 4.0 * d * d * d,
            },
            angular_second_derivative: 12.0 * d * d,
        }
    };
    assert!(matches!(
        survey_power_extrema_with_cancel(
            &GRID,
            coarse(12),
            3,
            || false,
            |w| Ok::<_, ()>(flat(w).sample),
            |w| Ok::<_, ()>(flat(w)),
        ),
        Err(ResonanceError::Invalid(_))
    ));
}

#[test]
fn sampled_and_unsampled_signed_power_poles_refuse_a_false_finite_survey() {
    let pole = 2.0_f64.sqrt();
    for sign in [-1.0, 1.0] {
        let evaluate = |w: f64| {
            if w == pole {
                return Err("physical pole");
            }
            let d = w - pole;
            Ok(HarmonicPowerCurvatureSample {
                sample: HarmonicPowerSample {
                    average_power_w: sign / (d * d),
                    angular_derivative: -2.0 * sign / (d * d * d),
                },
                angular_second_derivative: 6.0 * sign / (d * d * d * d),
            })
        };
        let result = survey_power_extrema_with_cancel(
            &[0.5, 1.0, 1.8, 2.3],
            ResonanceControls {
                absolute_angular_tolerance_rad_s: 1e-7,
                ..controls()
            },
            1,
            || false,
            |w| evaluate(w).map(|value| value.sample),
            evaluate,
        );
        assert!(matches!(
            result,
            Err(ResonanceError::Invalid(
                "nonzero physical power curvature with the original extremum orientation required"
            ))
        ));
        let result = survey_power_extrema_with_cancel(
            &[0.5, 1.0, pole, 1.8, 2.3],
            controls(),
            1,
            || false,
            |w| evaluate(w).map(|value| value.sample),
            |_| unreachable!(),
        );
        assert_eq!(
            result.unwrap_err(),
            ResonanceError::Evaluation("physical pole")
        );
    }
}
