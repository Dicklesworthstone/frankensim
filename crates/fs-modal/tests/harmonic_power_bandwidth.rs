//! G1/G3: actual physical-power half-height widths and parameter sensitivities.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_math::c64::C64;
use fs_modal::harmonic::resonance::{
    HarmonicPowerParameterJet, HarmonicPowerSample, PowerExtremum, PowerExtremumKind,
    PowerExtremumParameterResponse, ResonanceControls, ResonanceError, refine_power_bandwidth,
    refine_power_bandwidth_with_cancel, refine_power_bandwidth_with_parameter_derivative,
    refine_power_bandwidth_with_parameter_derivative_with_cancel,
    refine_power_extremum_with_parameter_derivative,
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
        absolute_angular_tolerance_rad_s: 1e-12,
        relative_angular_tolerance: 1e-12,
        maximum_evaluations: 256,
    }
}

fn near(actual: f64, expected: f64, tolerance: f64) {
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
                seed: 107,
                kernel_id: 283,
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
const B: f64 = 1.3;
const LM: f64 = 0.7;
const LC: f64 = 0.4;
const LK: f64 = 0.2;
const F: C64 = C64::new(2.0, -1.0);

#[derive(Clone, Copy)]
enum Channel {
    Input,
    Viscous,
    Port,
}

#[derive(Clone, Copy)]
enum Direction {
    Complete,
    Gain,
    Phase,
}

impl Direction {
    fn values(self) -> [f64; 7] {
        match self {
            Self::Complete => [0.12, 0.03, 0.4, -0.08, 0.11, -0.07, 0.05],
            Self::Gain | Self::Phase => [0.0; 7],
        }
    }

    fn force(self) -> C64 {
        match self {
            Self::Complete => C64::new(0.3, 0.4),
            Self::Gain => F.scale(0.7),
            Self::Phase => C64::new(-F.im, F.re).scale(0.9),
        }
    }
}

fn impedance(w: f64) -> [C64; 3] {
    [
        C64::new(LC, LK / w - w * LM),
        C64::new(0.0, -LK / (w * w) - LM),
        C64::new(0.0, 2.0 * LK / (w * w * w)),
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

fn sample(
    w: f64,
    channel: Channel,
    cx: &Cx<'_>,
) -> Result<HarmonicPowerSample, HarmonicResponseError> {
    let response = owner(w, cx)?.solve_with_frequency_derivative(
        &[F],
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
    direction: Direction,
    cx: &Cx<'_>,
) -> Result<HarmonicPowerParameterJet, HarmonicResponseError> {
    let [mt, ct, kt, bt, lmt, lct, lkt] = direction.values();
    let z = impedance(w);
    let response = owner(w, cx)?.solve_with_frequency_parameter_and_port_derivative(
        &[F],
        QuadraticParameterDirection {
            mass: &[mt],
            damping: &[ct],
            stiffness: &[kt],
        },
        None,
        None,
        Some(&[direction.force()]),
        None,
        &[z[1]],
        &[z[2]],
        HarmonicPortFrequencyParameterDirection {
            parameter: HarmonicPortParameterDirection {
                force_map: &[bt],
                impedance: &[C64::new(lct, lkt / w - w * lmt)],
            },
            impedance_angular_parameter_derivative: &[C64::new(0.0, -lkt / (w * w) - lmt)],
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
fn native_loaded_power_width_and_q_retain_all_physical_directions() {
    with_cx(|cx| {
        let mass = M + B * B * LM;
        let damping = C + B * B * LC;
        let stiffness = K + B * B * LK;
        let [mt, ct, kt, bt, lmt, lct, lkt] = Direction::Complete.values();
        let mass_t = mt + 2.0 * B * bt * LM + B * B * lmt;
        let damping_t = ct + 2.0 * B * bt * LC + B * B * lct;
        let stiffness_t = kt + 2.0 * B * bt * LK + B * B * lkt;
        let peak_w = (stiffness / mass).sqrt();
        let peak_w_t = 0.5 * peak_w * (stiffness_t / stiffness - mass_t / mass);
        // Exact half-height roots of w^2/((K-Mw^2)^2+C^2*w^2).
        // They are asymmetric about the velocity/physical-power maximum.
        let root = (damping * damping + 4.0 * mass * stiffness).sqrt();
        let root_t = (damping * damping_t + 2.0 * (mass_t * stiffness + mass * stiffness_t)) / root;
        let lower = (root - damping) / (2.0 * mass);
        let upper = (root + damping) / (2.0 * mass);
        let lower_t = (root_t - damping_t) / (2.0 * mass) - lower * mass_t / mass;
        let upper_t = (root_t + damping_t) / (2.0 * mass) - upper * mass_t / mass;
        let width = damping / mass;
        let width_t = damping_t / mass - damping * mass_t / (mass * mass);
        let q_t = peak_w_t / width - peak_w * width_t / (width * width);
        assert!(((lower + upper) * 0.5 - peak_w).abs() > 0.005);

        for channel in [Channel::Input, Channel::Viscous, Channel::Port] {
            let peak = refine_power_extremum_with_parameter_derivative(
                (0.8, 1.8),
                controls(),
                cx,
                |w| sample(w, channel, cx),
                |w| jet(w, channel, Direction::Complete, cx),
            )
            .unwrap();
            let first_calls = Cell::new(0);
            let edge_calls = Cell::new(0);
            let result = refine_power_bandwidth_with_parameter_derivative(
                &peak,
                (0.3, peak.extremum.angular_frequency_rad_s),
                (peak.extremum.angular_frequency_rad_s, 2.2),
                controls(),
                cx,
                |w| {
                    first_calls.set(first_calls.get() + 1);
                    sample(w, channel, cx)
                },
                |w| {
                    edge_calls.set(edge_calls.get() + 1);
                    jet(w, channel, Direction::Complete, cx)
                },
            )
            .unwrap();
            near(result.bandwidth.lower_angular_frequency_rad_s, lower, 3e-12);
            near(result.bandwidth.upper_angular_frequency_rad_s, upper, 3e-12);
            near(result.bandwidth.angular_bandwidth_rad_s, width, 4e-12);
            near(
                result.bandwidth.response_quality_factor,
                peak_w / width,
                3e-11,
            );
            near(
                result.lower_angular_frequency_derivative_rad_s,
                lower_t,
                2e-9,
            );
            near(
                result.upper_angular_frequency_derivative_rad_s,
                upper_t,
                2e-9,
            );
            near(result.angular_bandwidth_derivative_rad_s, width_t, 3e-9);
            near(result.response_quality_factor_derivative, q_t, 2e-8);
            assert!(
                result.bandwidth.lower_bracket_rad_s.1 < result.bandwidth.upper_bracket_rad_s.0
            );
            assert_eq!(edge_calls.get(), 2);
            assert_eq!(result.evaluations, first_calls.get() + 2);
            assert_eq!(result.bandwidth.evaluations, first_calls.get());
        }

        // The context-aware value-only entry observes the same original power.
        let peak = refine_power_extremum_with_parameter_derivative(
            (0.8, 1.8),
            controls(),
            cx,
            |w| sample(w, Channel::Port, cx),
            |w| jet(w, Channel::Port, Direction::Complete, cx),
        )
        .unwrap();
        let measured = refine_power_bandwidth(
            &peak.extremum,
            (0.3, peak.extremum.angular_frequency_rad_s),
            (peak.extremum.angular_frequency_rad_s, 2.2),
            controls(),
            cx,
            |w| sample(w, Channel::Port, cx),
        )
        .unwrap();
        near(measured.angular_bandwidth_rad_s, width, 4e-12);
    });
}

#[test]
fn common_source_gain_and_phase_do_not_change_a_native_absorbed_power_width() {
    with_cx(|cx| {
        for direction in [Direction::Gain, Direction::Phase] {
            let peak = refine_power_extremum_with_parameter_derivative(
                (0.8, 1.8),
                controls(),
                cx,
                |w| sample(w, Channel::Port, cx),
                |w| jet(w, Channel::Port, direction, cx),
            )
            .unwrap();
            let result = refine_power_bandwidth_with_parameter_derivative(
                &peak,
                (0.3, peak.extremum.angular_frequency_rad_s),
                (peak.extremum.angular_frequency_rad_s, 2.2),
                controls(),
                cx,
                |w| sample(w, Channel::Port, cx),
                |w| jet(w, Channel::Port, direction, cx),
            )
            .unwrap();
            near(result.lower_angular_frequency_derivative_rad_s, 0.0, 2e-9);
            near(result.upper_angular_frequency_derivative_rad_s, 0.0, 2e-9);
            near(result.angular_bandwidth_derivative_rad_s, 0.0, 3e-9);
            near(result.response_quality_factor_derivative, 0.0, 2e-8);
        }
    });
}

/// Analytic power line with independently changing gain, centre and half-width.
fn line(w: f64, power_units: f64, frequency_units: f64) -> HarmonicPowerParameterJet {
    let (centre, half_width, centre_t, half_width_t, gain_t) = (4.0, 0.25, 0.3, -0.04, 0.7);
    let x = (w / frequency_units - centre) / half_width;
    let x_t = -(centre_t + x * half_width_t) / half_width;
    let d = 1.0 + x * x;
    let p_x = -2.0 * x / (d * d);
    let p_xx = (6.0 * x * x - 2.0) / (d * d * d);
    HarmonicPowerParameterJet {
        sample: HarmonicPowerSample {
            average_power_w: power_units / d,
            angular_derivative: (power_units / frequency_units) * p_x / half_width,
        },
        angular_second_derivative: ((power_units / frequency_units) / frequency_units) * p_xx
            / (half_width * half_width),
        parameter_derivative: power_units * (gain_t / d + p_x * x_t),
        angular_parameter_derivative: (power_units / frequency_units)
            * ((gain_t * p_x + p_xx * x_t) / half_width
                - p_x * half_width_t / (half_width * half_width)),
    }
}

fn line_peak(power_units: f64, frequency_units: f64) -> PowerExtremumParameterResponse {
    let w = 4.0 * frequency_units;
    let jet = line(w, power_units, frequency_units);
    PowerExtremumParameterResponse {
        extremum: PowerExtremum {
            kind: PowerExtremumKind::Maximum,
            angular_frequency_rad_s: w,
            sample: jet.sample,
            angular_second_derivative: jet.angular_second_derivative,
            bracket_rad_s: (w, w),
            evaluations: 0,
        },
        angular_frequency_derivative_rad_s: 0.3 * frequency_units,
        average_power_derivative: 0.7 * power_units,
    }
}

#[test]
fn power_width_and_derivatives_survive_extreme_power_and_frequency_units() {
    for (power_exponent, frequency_exponent) in
        [(-700, 0), (700, 0), (-600, -600), (0, 0), (600, 600)]
    {
        let p = 2.0_f64.powi(power_exponent);
        let t = 2.0_f64.powi(frequency_exponent);
        let peak = line_peak(p, t);
        let result = refine_power_bandwidth_with_parameter_derivative_with_cancel(
            &peak,
            (3.5 * t, 4.0 * t),
            (4.0 * t, 4.5 * t),
            ResonanceControls {
                absolute_angular_tolerance_rad_s: 1e-12 * t,
                ..controls()
            },
            || false,
            |w| Ok::<_, ()>(line(w, p, t).sample),
            |w| Ok::<_, ()>(line(w, p, t)),
        )
        .unwrap();
        near(
            result.bandwidth.lower_angular_frequency_rad_s / t,
            3.75,
            3e-12,
        );
        near(
            result.bandwidth.upper_angular_frequency_rad_s / t,
            4.25,
            3e-12,
        );
        near(result.bandwidth.angular_bandwidth_rad_s / t, 0.5, 5e-12);
        near(result.bandwidth.response_quality_factor, 8.0, 1e-10);
        near(
            result.lower_angular_frequency_derivative_rad_s / t,
            0.34,
            2e-10,
        );
        near(
            result.upper_angular_frequency_derivative_rad_s / t,
            0.26,
            2e-10,
        );
        near(result.angular_bandwidth_derivative_rad_s / t, -0.08, 3e-10);
        near(result.response_quality_factor_derivative, 1.88, 2e-9);
    }
}

#[test]
fn one_budget_reserves_two_edge_jets_and_reuses_the_shared_peak_endpoint() {
    let peak = line_peak(1.0, 1.0);
    let first_calls = Cell::new(0);
    let endpoint_calls = Cell::new(0);
    let edge_calls = Cell::new(0);
    let result = refine_power_bandwidth_with_parameter_derivative_with_cancel(
        &peak,
        (3.5, 4.0),
        (4.0, 4.5),
        controls(),
        || false,
        |w| {
            first_calls.set(first_calls.get() + 1);
            if w == 4.0 {
                endpoint_calls.set(endpoint_calls.get() + 1);
            }
            Ok::<_, ()>(line(w, 1.0, 1.0).sample)
        },
        |w| {
            edge_calls.set(edge_calls.get() + 1);
            Ok::<_, ()>(line(w, 1.0, 1.0))
        },
    )
    .unwrap();
    assert_eq!(endpoint_calls.get(), 1);
    assert_eq!(edge_calls.get(), 2);
    assert_eq!(result.evaluations, first_calls.get() + 2);
    assert_eq!(result.bandwidth.evaluations, first_calls.get());

    for budget in [result.evaluations, result.evaluations - 1] {
        let first_calls = Cell::new(0);
        let edge_calls = Cell::new(0);
        let limited = refine_power_bandwidth_with_parameter_derivative_with_cancel(
            &peak,
            (3.5, 4.0),
            (4.0, 4.5),
            ResonanceControls {
                maximum_evaluations: budget,
                ..controls()
            },
            || false,
            |w| {
                first_calls.set(first_calls.get() + 1);
                Ok::<_, ()>(line(w, 1.0, 1.0).sample)
            },
            |w| {
                edge_calls.set(edge_calls.get() + 1);
                Ok::<_, ()>(line(w, 1.0, 1.0))
            },
        );
        if budget == result.evaluations {
            assert_eq!(limited.unwrap().evaluations, budget);
            assert_eq!(edge_calls.get(), 2);
        } else {
            assert_eq!(
                limited.unwrap_err(),
                ResonanceError::EvaluationBudget {
                    maximum_evaluations: budget
                }
            );
            assert_eq!(edge_calls.get(), 0);
        }
        assert!(first_calls.get() + edge_calls.get() <= budget);
    }
}

#[test]
fn missing_crossings_false_peaks_and_unresolved_narrow_power_lines_refuse() {
    let peak = line_peak(1.0, 1.0).extremum;
    for case in 0..7 {
        let mut invalid = peak;
        match case {
            0 => invalid.kind = PowerExtremumKind::Minimum,
            1 => invalid.sample.average_power_w = 0.0,
            2 => invalid.sample.average_power_w = -1.0,
            3 => invalid.angular_frequency_rad_s = 0.0,
            4 => invalid.angular_second_derivative = 0.0,
            5 => invalid.angular_second_derivative = 1.0,
            _ => invalid.sample.angular_derivative = f64::NAN,
        }
        assert!(matches!(
            refine_power_bandwidth_with_cancel(
                &invalid,
                (3.5, 4.0),
                (4.0, 4.5),
                controls(),
                || false,
                |_| -> Result<HarmonicPowerSample, ()> { unreachable!() },
            ),
            Err(ResonanceError::Invalid(_))
        ));
    }
    assert!(matches!(
        refine_power_bandwidth_with_cancel(
            &peak,
            (3.9, 4.0),
            (4.0, 4.5),
            controls(),
            || false,
            |w| Ok::<_, ()>(line(w, 1.0, 1.0).sample),
        ),
        Err(ResonanceError::InvalidBracket(_))
    ));
    // A zero/incorrect tangent cannot certify a simple oriented crossing.
    for slope_scale in [0.0, -1.0] {
        assert!(matches!(
            refine_power_bandwidth_with_cancel(
                &peak,
                (3.5, 4.0),
                (4.0, 4.5),
                controls(),
                || false,
                |w| {
                    let mut sample = line(w, 1.0, 1.0).sample;
                    sample.angular_derivative *= slope_scale;
                    Ok::<_, ()>(sample)
                },
            ),
            Err(ResonanceError::InvalidBracket(_))
        ));
    }

    let half_width = 1e-8;
    let narrow = |w: f64| {
        let x = (w - 4.0) / half_width;
        let d = 1.0 + x * x;
        Ok::<_, ()>(HarmonicPowerSample {
            average_power_w: 1.0 / d,
            angular_derivative: -2.0 * x / (half_width * d * d),
        })
    };
    let narrow_peak = PowerExtremum {
        angular_second_derivative: -2.0 / (half_width * half_width),
        ..peak
    };
    assert!(matches!(
        refine_power_bandwidth_with_cancel(
            &narrow_peak,
            (3.9, 4.0),
            (4.0, 4.1),
            ResonanceControls {
                absolute_angular_tolerance_rad_s: 1e-3,
                relative_angular_tolerance: 0.0,
                ..controls()
            },
            || false,
            narrow,
        ),
        Err(ResonanceError::InvalidBracket(
            "power crossing brackets do not resolve a positive linewidth"
        ))
    ));
    let resolved = refine_power_bandwidth_with_cancel(
        &narrow_peak,
        (3.9, 4.0),
        (4.0, 4.1),
        ResonanceControls {
            absolute_angular_tolerance_rad_s: 1e-14,
            relative_angular_tolerance: 0.0,
            ..controls()
        },
        || false,
        narrow,
    )
    .unwrap();
    assert!((resolved.angular_bandwidth_rad_s - 2.0 * half_width).abs() < 2e-14);
    assert!(resolved.lower_bracket_rad_s.1 < resolved.upper_bracket_rad_s.0);
}

#[test]
fn inconsistent_edge_data_source_failure_and_final_call_cancellation_are_atomic() {
    let peak = line_peak(1.0, 1.0);
    for case in 0..3 {
        let result = refine_power_bandwidth_with_parameter_derivative_with_cancel(
            &peak,
            (3.5, 4.0),
            (4.0, 4.5),
            controls(),
            || false,
            |w| Ok::<_, ()>(line(w, 1.0, 1.0).sample),
            |w| {
                let mut jet = line(w, 1.0, 1.0);
                match case {
                    0 => jet.sample.average_power_w *= 1.001,
                    1 => jet.sample.angular_derivative *= 1.001,
                    _ => jet.parameter_derivative = f64::INFINITY,
                }
                Ok::<_, ()>(jet)
            },
        );
        assert!(matches!(result, Err(ResonanceError::Invalid(_))));
    }

    let cancelled = Cell::new(false);
    let edges = Cell::new(0);
    assert_eq!(
        refine_power_bandwidth_with_parameter_derivative_with_cancel(
            &peak,
            (3.5, 4.0),
            (4.0, 4.5),
            controls(),
            || cancelled.get(),
            |w| Ok::<_, ()>(line(w, 1.0, 1.0).sample),
            |w| {
                edges.set(edges.get() + 1);
                cancelled.set(true);
                Ok::<_, ()>(line(w, 1.0, 1.0))
            },
        )
        .unwrap_err(),
        ResonanceError::Cancelled,
    );
    assert_eq!(edges.get(), 1);

    let edges = Cell::new(0);
    assert_eq!(
        refine_power_bandwidth_with_parameter_derivative_with_cancel(
            &peak,
            (3.5, 4.0),
            (4.0, 4.5),
            controls(),
            || false,
            |w| Ok::<_, &'static str>(line(w, 1.0, 1.0).sample),
            |w| {
                edges.set(edges.get() + 1);
                if edges.get() == 2 {
                    Err("original edge constitutive refusal")
                } else {
                    Ok(line(w, 1.0, 1.0))
                }
            },
        )
        .unwrap_err(),
        ResonanceError::Evaluation("original edge constitutive refusal"),
    );
    assert_eq!(edges.get(), 2);
    assert_eq!(
        refine_power_bandwidth_with_cancel(
            &peak.extremum,
            (3.5, 4.0),
            (4.0, 4.5),
            controls(),
            || false,
            |_| Err::<HarmonicPowerSample, _>("original sampled pole"),
        )
        .unwrap_err(),
        ResonanceError::Evaluation("original sampled pole"),
    );
    assert_eq!(
        refine_power_bandwidth_with_cancel(
            &peak.extremum,
            (3.5, 4.0),
            (4.0, 4.5),
            controls(),
            || true,
            |_| -> Result<HarmonicPowerSample, ()> { unreachable!() },
        )
        .unwrap_err(),
        ResonanceError::Cancelled,
    );
}
