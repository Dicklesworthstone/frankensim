//! Source-owned observed peaks: no pHS conversion or second search algorithm.
use core::convert::Infallible;
use std::cell::Cell;

use fs_math::c64::C64;
use fs_phs::harmonic::resonance::{
    HarmonicObservation, ResonanceControls, ResonanceSearchError, search_resonance,
};
use fs_phs::harmonic::{HarmonicBudget, HarmonicError, HarmonicModel};

fn controls() -> ResonanceControls {
    ResonanceControls {
        absolute_frequency_tolerance_rad_s: 1e-10,
        relative_frequency_tolerance: 0.0,
        maximum_evaluations: 256,
        measure_half_power_bandwidth: true,
    }
}

// Complete scalar mechanical source: (k-m*omega^2-i*c*omega)q=F.
// Its observed velocity, not its free pole, is supplied to the search.
fn mobility(omega: f64, force: C64, damping: f64, tangent: bool) -> (C64, C64, f64) {
    let dynamic = C64::new(200.0 - 2.0 * omega * omega, -damping * omega);
    let dynamic_derivative = C64::new(-4.0 * omega, -damping);
    let displacement = force / dynamic;
    let velocity = C64::new(0.0, -omega) * displacement;
    let mut residual =
        (dynamic * displacement - force).abs() / (dynamic.abs() * displacement.abs() + force.abs());
    let derivative = if tangent {
        let dq = -(dynamic_derivative / dynamic) * displacement;
        residual = residual.max(
            (dynamic * dq + dynamic_derivative * displacement).abs()
                / (dynamic.abs() * dq.abs() + dynamic_derivative.abs() * displacement.abs()),
        );
        C64::new(0.0, -omega) * dq + C64::new(0.0, -1.0) * displacement
    } else {
        // A source can skip its tangent entirely at half-height queries.
        C64::new(f64::NAN, f64::NAN)
    };
    (velocity, derivative, residual)
}

#[test]
fn external_mechanical_source_resolves_high_q_and_retains_physical_query_accounts() {
    let damping = 0.002;
    for scale in [1e-120, 1.0, 1e120] {
        let force = C64::new(scale, -0.3 * scale);
        let queries = Cell::new(0);
        let tangents = Cell::new(0);
        let maximum_residual = Cell::new(0.0_f64);
        let result = search_resonance(
            (8.0, 12.0),
            controls(),
            || false,
            |omega, need_derivative, cancelled| {
                assert!(!cancelled());
                queries.set(queries.get() + 1);
                tangents.set(tangents.get() + usize::from(need_derivative));
                let sample = mobility(omega, force, damping, need_derivative);
                maximum_residual.set(maximum_residual.get().max(sample.2));
                Ok::<_, Infallible>(sample)
            },
        )
        .unwrap();
        let band = result.half_power.unwrap();
        assert!((result.peak.angular_frequency_rad_s - 10.0).abs() <= 1e-10);
        assert!((band.bandwidth_rad_s - damping / 2.0).abs() <= 2e-10);
        assert!((band.response_quality_factor / 10_000.0 - 1.0).abs() <= 2e-7);
        let exact = force.scale(1.0 / damping);
        assert!((result.observation - exact).abs() <= 1e-6 * exact.abs());
        assert_eq!(result.evaluations, queries.get());
        assert!(tangents.get() > 2 && tangents.get() < queries.get());
        assert_eq!(result.maximum_relative_residual, maximum_residual.get());
        assert!(result.maximum_relative_residual <= 1e-14);
    }
}

#[test]
fn source_feedthrough_interference_moves_the_observed_peak_and_both_half_heights() {
    // y=1+i*h/(gamma-i*(omega-center)). The direct term changes the
    // amplitude maximum away from the pole center and makes asymmetric edges.
    let (center, gamma, h) = (100.0_f64, 2.0_f64, 2.0_f64);
    let result = search_resonance(
        (center - 10.0 * gamma, center + 0.8 * gamma),
        controls(),
        || false,
        |omega, need_derivative, _| {
            let denominator = C64::new(gamma, center - omega);
            let value = C64::ONE + C64::new(0.0, h) / denominator;
            let derivative = if need_derivative {
                -C64::from_re(h) / (denominator * denominator)
            } else {
                C64::ZERO
            };
            Ok::<_, Infallible>((value, derivative, 0.0))
        },
    )
    .unwrap();
    let offset = 0.5 * (h - (h * h + 4.0 * gamma * gamma).sqrt());
    assert!((result.peak.angular_frequency_rad_s - center - offset).abs() <= 1e-10);
    assert!(result.peak.angular_frequency_rad_s < center - 0.5 * gamma);
    let peak_squared = (gamma * gamma + (offset - h).powi(2)) / (gamma * gamma + offset * offset);
    let coefficient = 1.0 - 0.5 * peak_squared;
    let root = (h * h - coefficient * (h * h + coefficient * gamma * gamma)).sqrt();
    let band = result.half_power.unwrap();
    assert!(
        (band.lower.angular_frequency_rad_s - center - (h + root) / coefficient).abs() <= 1e-10
    );
    assert!(
        (band.upper.angular_frequency_rad_s - center - (h - root) / coefficient).abs() <= 1e-10
    );
}

#[derive(Debug)]
struct NativeRefusal {
    message: String,
    port: usize,
}

#[test]
fn native_source_error_survives_without_clone_equality_or_error_trait_bounds() {
    let mut queries = 0;
    let result = search_resonance(
        (8.0, 12.0),
        controls(),
        || false,
        |omega, need_derivative, _| {
            queries += 1;
            if !need_derivative {
                return Err(NativeRefusal {
                    message: "original load refused its half-height query".to_owned(),
                    port: 19,
                });
            }
            Ok(mobility(omega, C64::ONE, 0.2, true))
        },
    );
    assert!(queries > 3);
    match result {
        Err(ResonanceSearchError::Source(error)) => {
            assert_eq!(error.port, 19);
            assert_eq!(error.message, "original load refused its half-height query");
        }
        other => panic!("native refusal was replaced or a partial peak escaped: {other:?}"),
    }
}

#[test]
fn source_evaluation_cap_and_cancellation_stop_before_another_physical_query() {
    let mut queries = 0;
    let result = search_resonance(
        (8.0, 12.0),
        ResonanceControls {
            maximum_evaluations: 3,
            ..controls()
        },
        || false,
        |omega, need_derivative, _| {
            queries += 1;
            Ok::<_, Infallible>(mobility(omega, C64::ONE, 0.2, need_derivative))
        },
    );
    assert!(matches!(
        result,
        Err(ResonanceSearchError::EvaluationLimit { maximum: 3 })
    ));
    assert_eq!(queries, 3);

    let stop = Cell::new(false);
    let mut queries = 0;
    let result = search_resonance(
        (8.0, 12.0),
        controls(),
        || stop.get(),
        |omega, need_derivative, _| {
            queries += 1;
            if queries == 2 {
                stop.set(true);
            }
            Ok::<_, Infallible>(mobility(omega, C64::ONE, 0.2, need_derivative))
        },
    );
    assert!(matches!(result, Err(ResonanceSearchError::Cancelled)));
    assert_eq!(queries, 2);
}

#[test]
fn malformed_source_observations_and_controls_refuse_without_a_partial_result() {
    for sample in [
        (C64::new(f64::NAN, 0.0), C64::ONE, 0.0),
        (C64::ONE, C64::new(0.0, f64::INFINITY), 0.0),
        (C64::ONE, C64::ONE, f64::NAN),
        (C64::ONE, C64::ONE, -1e-12),
    ] {
        let mut queries = 0;
        let result = search_resonance(
            (8.0, 12.0),
            controls(),
            || false,
            |_, _, _| {
                queries += 1;
                Ok::<_, Infallible>(sample)
            },
        );
        assert!(matches!(result, Err(ResonanceSearchError::Invalid(_))));
        assert_eq!(queries, 1);
    }
    let result = search_resonance::<Infallible, _, _>(
        (12.0, 8.0),
        controls(),
        || false,
        |_, _, _| panic!("invalid controls evaluated the source"),
    );
    assert!(matches!(result, Err(ResonanceSearchError::Invalid(_))));
}

#[test]
fn existing_phs_method_and_external_evaluator_share_identical_peak_and_query_semantics() {
    let system = fs_phs::mass_spring_damper(2.0, 200.0, 0.8).unwrap();
    let model = HarmonicModel::from_ode(
        &system,
        HarmonicBudget {
            max_dimension: 4,
            max_dense_bytes: 1_000_000,
            max_dense_work: 1_000_000,
            relative_residual_tolerance: 1e-12,
        },
        || false,
    )
    .unwrap();
    let input = [C64::new(0.8, -0.3)];
    let expected = model
        .driven_resonance(
            &input,
            HarmonicObservation::Output(&[1.0]),
            (8.0, 12.0),
            controls(),
            || false,
        )
        .unwrap();
    let actual = search_resonance(
        (8.0, 12.0),
        controls(),
        || false,
        |omega, need_derivative, cancelled| -> Result<_, HarmonicError> {
            let factor = model.factor(omega, &mut *cancelled)?;
            if need_derivative {
                let result =
                    factor.response_with_frequency_derivative(&input, None, &mut *cancelled)?;
                Ok((
                    result.response.output[0],
                    result.output_derivative[0],
                    result
                        .relative_derivative_residual
                        .max(result.response.relative_residual),
                ))
            } else {
                let result = factor.response(&input, &mut *cancelled)?;
                Ok((result.output[0], C64::ZERO, result.relative_residual))
            }
        },
    )
    .unwrap();
    assert_eq!(actual, expected);
}
