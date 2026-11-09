//! Actual signed physical-power extrema and material/attachment sensitivity.
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_math::c64::C64;
use fs_modal::harmonic::resonance::{
    HarmonicPowerCurvatureSample, HarmonicPowerParameterJet, HarmonicPowerSample,
    PowerExtremumKind, ResonanceControls, ResonanceError, refine_power_extremum_with_cancel,
    refine_power_extremum_with_parameter_derivative,
    refine_power_extremum_with_parameter_derivative_with_cancel,
};
use fs_modal::harmonic::{
    HarmonicFrequencyControls, HarmonicFrequencyParameterControls, HarmonicPortControls,
    HarmonicPortFrequencyParameterDirection, HarmonicPortParameterDirection,
    HarmonicResponseControls, HarmonicResponseError, PreparedHarmonicResponse,
};
use fs_modal::quadratic_modes::QuadraticParameterDirection;
use std::cell::Cell;

fn controls() -> ResonanceControls {
    ResonanceControls {
        absolute_angular_tolerance_rad_s: 1e-11,
        relative_angular_tolerance: 1e-11,
        maximum_evaluations: 192,
    }
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1.0),
        "{actual:?} != {expected:?}"
    );
}

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 89,
                kernel_id: 277,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}

const M: f64 = 1.7;
const C: f64 = 0.07;
const K: f64 = 4.3;
const MT: f64 = 0.12;
const CT: f64 = 0.03;
const KT: f64 = 0.4;
const B: f64 = 1.3;
const BT: f64 = -0.08;
const LM: f64 = 0.7;
const LC: f64 = 0.4;
const LK: f64 = 0.2;
const LMT: f64 = 0.11;
const LCT: f64 = -0.07;
const LKT: f64 = 0.05;
const F: f64 = 2.0;
const FT: C64 = C64::new(0.3, 0.4);

fn impedance(w: f64) -> [C64; 5] {
    [
        C64::new(LC, LK / w - w * LM),
        C64::new(0.0, -LK / (w * w) - LM),
        C64::new(0.0, 2.0 * LK / (w * w * w)),
        C64::new(LCT, LKT / w - w * LMT),
        C64::new(0.0, -LKT / (w * w) - LMT),
    ]
}

fn owner(w: f64, cx: &Cx<'_>) -> Result<PreparedHarmonicResponse, HarmonicResponseError> {
    PreparedHarmonicResponse::try_with_port_impedance(
        &[M],
        &[C],
        &[K],
        1,
        &[B],
        1,
        &[impedance(w)[0]],
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

#[derive(Clone, Copy)]
enum Channel {
    Input,
    Viscous,
    Port,
}

fn sample(
    w: f64,
    channel: Channel,
    cx: &Cx<'_>,
) -> Result<HarmonicPowerSample, HarmonicResponseError> {
    let response = owner(w, cx)?.solve_with_frequency_derivative(
        &[C64::from_re(F)],
        None,
        Some(&[impedance(w)[1]]),
        HarmonicFrequencyControls {
            maximum_additional_payload_bytes: 1 << 20,
            maximum_work: 1_000_000,
        },
        cx,
    )?;
    match channel {
        Channel::Input => Ok(response.input_power_sample()),
        Channel::Viscous => Ok(response.viscous_power_sample()),
        Channel::Port => response.port_power_sample(),
    }
}

fn jet(
    w: f64,
    channel: Channel,
    cx: &Cx<'_>,
) -> Result<HarmonicPowerParameterJet, HarmonicResponseError> {
    let z = impedance(w);
    let response = owner(w, cx)?.solve_with_frequency_parameter_and_port_derivative(
        &[C64::from_re(F)],
        QuadraticParameterDirection {
            mass: &[MT],
            damping: &[CT],
            stiffness: &[KT],
        },
        None,
        None,
        Some(&[FT]),
        None,
        &[z[1]],
        &[z[2]],
        HarmonicPortFrequencyParameterDirection {
            parameter: HarmonicPortParameterDirection {
                force_map: &[BT],
                impedance: &[z[3]],
            },
            impedance_angular_parameter_derivative: &[z[4]],
        },
        HarmonicFrequencyParameterControls {
            maximum_additional_payload_bytes: 1 << 20,
            maximum_work: 1_000_000,
        },
        cx,
    )?;
    match channel {
        Channel::Input => Ok(response.input_power_parameter_jet()),
        Channel::Viscous => Ok(response.viscous_power_parameter_jet()),
        Channel::Port => response.port_power_parameter_jet(),
    }
}

#[test]
fn native_absorbed_source_and_viscous_peaks_match_original_loaded_oscillator() {
    with_cx(|cx| {
        // Independently combine the original physical inertia, damping and
        // spring, including BOTH occurrences of the changing attachment.
        let mass = M + B * B * LM;
        let damping = C + B * B * LC;
        let stiffness = K + B * B * LK;
        let mass_t = MT + 2.0 * B * BT * LM + B * B * LMT;
        let damping_t = CT + 2.0 * B * BT * LC + B * B * LCT;
        let stiffness_t = KT + 2.0 * B * BT * LK + B * B * LKT;
        let frequency = (stiffness / mass).sqrt();
        let frequency_t = 0.5 * frequency * (stiffness_t / stiffness - mass_t / mass);
        let displacement_peak = (stiffness / mass - damping * damping / (2.0 * mass * mass)).sqrt();
        assert!(frequency - displacement_peak > 0.005);

        for channel in [Channel::Input, Channel::Viscous, Channel::Port] {
            let first_calls = Cell::new(0);
            let jet_calls = Cell::new(0);
            let result = refine_power_extremum_with_parameter_derivative(
                (0.8, 1.8),
                controls(),
                cx,
                |w| {
                    first_calls.set(first_calls.get() + 1);
                    sample(w, channel, cx)
                },
                |w| {
                    jet_calls.set(jet_calls.get() + 1);
                    jet(w, channel, cx)
                },
            )
            .unwrap();
            let (coefficient, coefficient_t) = match channel {
                Channel::Input => (damping, damping_t),
                Channel::Viscous => (C, CT),
                Channel::Port => (B * B * LC, 2.0 * B * BT * LC + B * B * LCT),
            };
            let peak_power = 0.5 * F * F * coefficient / (damping * damping);
            let peak_power_t = peak_power
                * (2.0 * FT.re / F + coefficient_t / coefficient - 2.0 * damping_t / damping);
            assert_eq!(result.extremum.kind, PowerExtremumKind::Maximum);
            near(result.extremum.angular_frequency_rad_s, frequency, 2e-11);
            near(result.angular_frequency_derivative_rad_s, frequency_t, 2e-9);
            near(result.extremum.sample.average_power_w, peak_power, 2e-11);
            near(result.average_power_derivative, peak_power_t, 2e-9);
            assert_eq!(jet_calls.get(), 1);
            assert_eq!(result.extremum.evaluations, first_calls.get() + 1);
        }
    });
}

fn parabola(
    w: f64,
    center: f64,
    offset: f64,
    orientation: f64,
    units: f64,
) -> HarmonicPowerParameterJet {
    let distance = w - center;
    HarmonicPowerParameterJet {
        sample: HarmonicPowerSample {
            average_power_w: units * (offset + orientation * distance * distance),
            angular_derivative: units * 2.0 * orientation * distance,
        },
        angular_second_derivative: units * 2.0 * orientation,
        parameter_derivative: units * (0.7 - 2.0 * orientation * distance * 0.3),
        angular_parameter_derivative: units * -2.0 * orientation * 0.3,
    }
}

#[test]
fn signed_and_zero_power_extrema_keep_orientation_units_and_exact_stationary_sample() {
    for exponent in [-700, 0, 700] {
        let units = f64::from_bits(((1023 + exponent) as u64) << 52);
        for offset in [-3.0, 0.0, 3.0] {
            for orientation in [-1.0, 1.0] {
                let exact_calls = Cell::new(0);
                let final_calls = Cell::new(0);
                let result = refine_power_extremum_with_parameter_derivative_with_cancel(
                    (1.0, 3.0),
                    controls(),
                    || false,
                    |w| {
                        if w == 2.0 {
                            exact_calls.set(exact_calls.get() + 1);
                        }
                        Ok::<_, ()>(parabola(w, 2.0, offset, orientation, units).sample)
                    },
                    |w| {
                        final_calls.set(final_calls.get() + 1);
                        Ok::<_, ()>(parabola(w, 2.0, offset, orientation, units))
                    },
                )
                .unwrap();
                assert_eq!(
                    result.extremum.kind,
                    if orientation < 0.0 {
                        PowerExtremumKind::Maximum
                    } else {
                        PowerExtremumKind::Minimum
                    }
                );
                assert_eq!(result.extremum.angular_frequency_rad_s, 2.0);
                assert_eq!(result.extremum.sample.average_power_w / units, offset);
                assert_eq!(result.extremum.sample.angular_derivative, 0.0);
                near(result.angular_frequency_derivative_rad_s, 0.3, 1e-14);
                near(result.average_power_derivative / units, 0.7, 1e-14);
                assert_eq!(exact_calls.get(), 1);
                assert_eq!(final_calls.get(), 1);
            }
        }
    }
}

#[test]
fn total_power_derivative_retains_motion_of_the_returned_finite_resolution_sample() {
    let result = refine_power_extremum_with_parameter_derivative_with_cancel(
        (1.0, 2.0),
        ResonanceControls {
            absolute_angular_tolerance_rad_s: 0.3,
            ..controls()
        },
        || false,
        |w| Ok::<_, ()>(parabola(w, 1.3, -2.0, -1.0, 1.0).sample),
        |w| Ok::<_, ()>(parabola(w, 1.3, -2.0, -1.0, 1.0)),
    )
    .unwrap();
    let at_sample = parabola(
        result.extremum.angular_frequency_rad_s,
        1.3,
        -2.0,
        -1.0,
        1.0,
    );
    assert!(result.extremum.sample.angular_derivative.abs() > 0.1);
    assert!((at_sample.parameter_derivative - 0.7).abs() > 0.01);
    near(result.average_power_derivative, 0.7, 1e-14);
}

#[test]
fn divergent_passive_and_signed_power_poles_refuse_false_stationary_extrema() {
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
        assert!(matches!(
            refine_power_extremum_with_cancel(
                (0.7, 2.1),
                ResonanceControls {
                    absolute_angular_tolerance_rad_s: 1e-7,
                    ..controls()
                },
                || false,
                |w| evaluate(w).map(|jet| jet.sample),
                evaluate,
            ),
            Err(ResonanceError::Invalid(
                "nonzero physical power curvature with the original extremum orientation required"
            ))
        ));
    }
    // A sampled pole retains the original evaluator's exact error.
    assert_eq!(
        refine_power_extremum_with_cancel(
            (1.0, 3.0),
            controls(),
            || false,
            |w| if w == 2.0 {
                Err("sampled pole")
            } else {
                Ok(parabola(w, 2.0, 1.0, -1.0, 1.0).sample)
            },
            |_| unreachable!(),
        )
        .unwrap_err(),
        ResonanceError::Evaluation("sampled pole"),
    );
}

#[test]
fn final_curvature_and_original_signed_samples_are_required_without_a_zero_floor() {
    for case in 0..5 {
        let result = refine_power_extremum_with_parameter_derivative_with_cancel(
            (1.0, 3.0),
            controls(),
            || false,
            |w| Ok::<_, ()>(parabola(w, 2.0, 0.0, -1.0, 1.0).sample),
            |w| {
                let mut jet = parabola(w, 2.0, 0.0, -1.0, 1.0);
                match case {
                    0 => jet.angular_second_derivative = 0.0,
                    1 => jet.angular_second_derivative = 2.0,
                    2 => jet.sample.average_power_w = f64::from_bits(1),
                    3 => jet.sample.angular_derivative = f64::from_bits(1),
                    _ => jet.parameter_derivative = f64::INFINITY,
                }
                Ok::<_, ()>(jet)
            },
        );
        assert!(matches!(result, Err(ResonanceError::Invalid(_))));
    }
    // A genuinely flat quartic extremum has oriented first slopes but zero
    // physical curvature, so its implicit parameter sensitivity is undefined.
    assert!(matches!(
        refine_power_extremum_with_cancel(
            (1.0, 3.0),
            controls(),
            || false,
            |w| {
                let d = w - 2.0;
                Ok::<_, ()>(HarmonicPowerSample {
                    average_power_w: 1.0 - d * d * d * d,
                    angular_derivative: -4.0 * d * d * d,
                })
            },
            |w| {
                let d = w - 2.0;
                Ok::<_, ()>(HarmonicPowerCurvatureSample {
                    sample: HarmonicPowerSample {
                        average_power_w: 1.0 - d * d * d * d,
                        angular_derivative: -4.0 * d * d * d,
                    },
                    angular_second_derivative: -12.0 * d * d,
                })
            },
        ),
        Err(ResonanceError::Invalid(
            "nonzero physical power curvature with the original extremum orientation required"
        )),
    ));
    // A stationary boundary alone supplies no two-sided extremum evidence.
    assert!(matches!(
        refine_power_extremum_with_cancel(
            (2.0, 3.0),
            controls(),
            || false,
            |w| Ok::<_, ()>(parabola(w, 2.0, 0.0, -1.0, 1.0).sample),
            |_| unreachable!(),
        ),
        Err(ResonanceError::InvalidBracket(_)),
    ));
}

#[test]
fn signed_power_query_reserves_complete_jet_budget_and_cancellation_is_atomic() {
    let first_calls = Cell::new(0);
    let final_calls = Cell::new(0);
    let result = refine_power_extremum_with_parameter_derivative_with_cancel(
        (1.0, 3.0),
        ResonanceControls {
            maximum_evaluations: 4,
            ..controls()
        },
        || false,
        |w| {
            first_calls.set(first_calls.get() + 1);
            Ok::<_, ()>(parabola(w, 2.0, 1.0, -1.0, 1.0).sample)
        },
        |w| {
            final_calls.set(final_calls.get() + 1);
            Ok::<_, ()>(parabola(w, 2.0, 1.0, -1.0, 1.0))
        },
    );
    assert!(matches!(
        result,
        Err(ResonanceError::EvaluationBudget {
            maximum_evaluations: 4
        })
    ));
    assert_eq!(first_calls.get(), 3);
    assert_eq!(final_calls.get(), 0);

    let cancelled = Cell::new(false);
    let result = refine_power_extremum_with_parameter_derivative_with_cancel(
        (1.0, 3.0),
        controls(),
        || cancelled.get(),
        |w| Ok::<_, ()>(parabola(w, 2.0, 1.0, -1.0, 1.0).sample),
        |w| {
            cancelled.set(true);
            Ok::<_, ()>(parabola(w, 2.0, 1.0, -1.0, 1.0))
        },
    );
    assert!(matches!(result, Err(ResonanceError::Cancelled)));

    assert_eq!(
        refine_power_extremum_with_parameter_derivative_with_cancel(
            (1.0, 3.0),
            controls(),
            || false,
            |w| Ok(parabola(w, 2.0, 1.0, -1.0, 1.0).sample),
            |_| Err::<HarmonicPowerParameterJet, _>("native mixed solve refused"),
        )
        .unwrap_err(),
        ResonanceError::Evaluation("native mixed solve refused"),
    );
    let simple = refine_power_extremum_with_cancel(
        (1.0, 3.0),
        controls(),
        || false,
        |w| Ok::<_, ()>(parabola(w, 2.0, -3.0, 1.0, 1.0).sample),
        |w| {
            let jet = parabola(w, 2.0, -3.0, 1.0, 1.0);
            Ok::<_, ()>(HarmonicPowerCurvatureSample {
                sample: jet.sample,
                angular_second_derivative: jet.angular_second_derivative,
            })
        },
    )
    .unwrap();
    assert_eq!(simple.kind, PowerExtremumKind::Minimum);
    assert_eq!(simple.sample.average_power_w, -3.0);
    assert_eq!(simple.angular_frequency_rad_s, 2.0);
}

#[test]
fn native_port_power_adapters_refuse_absent_loads() {
    with_cx(|cx| {
        let owner = PreparedHarmonicResponse::try_new(
            &[1.0],
            &[0.1],
            &[1.0],
            1,
            0.8,
            HarmonicResponseControls {
                maximum_coordinates: 2,
                maximum_payload_bytes: 1 << 20,
                maximum_relative_residual: 2e-11,
            },
            cx,
        )
        .unwrap();
        let zero = [0.0];
        let result = owner
            .solve_with_frequency_parameter_derivative(
                &[C64::ONE],
                QuadraticParameterDirection {
                    mass: &zero,
                    damping: &zero,
                    stiffness: &zero,
                },
                None,
                None,
                None,
                None,
                HarmonicFrequencyParameterControls {
                    maximum_additional_payload_bytes: 1 << 20,
                    maximum_work: 1_000_000,
                },
                cx,
            )
            .unwrap();
        assert!(
            result
                .frequency_response
                .frequency_response
                .port_power_sample()
                .is_err()
        );
        assert!(
            result
                .frequency_response
                .port_power_curvature_sample()
                .is_err()
        );
        assert!(result.port_power_parameter_jet().is_err());
        let source = result.frequency_response.input_power_curvature_sample();
        let viscous = result.frequency_response.viscous_power_curvature_sample();
        near(
            source.sample.average_power_w,
            viscous.sample.average_power_w,
            1e-12,
        );
        near(
            source.angular_second_derivative,
            viscous.angular_second_derivative,
            1e-11,
        );
    });
}
