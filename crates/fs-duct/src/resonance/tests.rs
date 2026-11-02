use super::*;
use std::cell::Cell;
use std::collections::BTreeSet;

fn receptance(omega: f64, natural: f64, damping_ratio: f64) -> C64 {
    C64::ONE
        / C64::new(
            natural * natural - omega * omega,
            -2.0 * damping_ratio * natural * omega,
        )
}

#[test]
fn oscillator_peak_and_half_power_are_observable_not_natural_frequency() {
    let natural = 731.3;
    let zeta: f64 = 0.2;
    let x = 1.0 - 2.0 * zeta * zeta;
    let split = 2.0 * zeta * (1.0 - zeta * zeta).sqrt();
    let expected = natural * x.sqrt();
    let half = [natural * (x - split).sqrt(), natural * (x + split).sqrt()];
    let result = refine_response_resonance(
        [0.4 * natural, 1.5 * natural],
        ResonanceOptions::default(),
        |w| Ok(receptance(w, natural, zeta)),
        || false,
    )
    .unwrap();
    assert!((result.angular_frequency_rad_s / expected - 1.0).abs() < 2e-8);
    assert!((result.angular_frequency_rad_s / natural - 1.0).abs() > 0.03);
    for (&actual, &want) in result.half_power_frequencies_rad_s.iter().zip(&half) {
        assert!((actual / want - 1.0).abs() < 1e-9);
    }
    let q = expected / (half[1] - half[0]);
    assert!((result.bandwidth_quality_factor / q - 1.0).abs() < 3e-8);
    assert!((result.bandwidth_quality_factor * 2.0 * zeta - 1.0).abs() > 0.05);
    assert_eq!(
        result.response,
        receptance(result.angular_frequency_rad_s, natural, zeta)
    );
}

#[test]
fn very_high_q_off_grid_peak_and_gain_changes_preserve_width_without_repeated_solves() {
    let natural = 1234.56789123;
    let quality = 1.0e7;
    let damping = natural / quality;
    // Velocity/force has an exact peak at natural and bandwidth C/M, unlike
    // displacement/force. A coarse grid undershoots this peak by five decades.
    let mobility =
        |w: f64| C64::new(0.0, -w) / C64::new((natural - w) * (natural + w), -damping * w);
    let coarse = (0..9)
        .map(|i| mobility(1000.0 + f64::from(i) * 62.5).abs())
        .fold(0.0_f64, f64::max);
    assert!(coarse * damping < 1.0e-4);
    for gain in [1e-200, 1.0, 1e200] {
        let mut frequencies = BTreeSet::new();
        let result = refine_response_resonance(
            // The frequency tolerance follows the current peak region, not
            // this extremely small outer lower endpoint.
            [1e-150, 1500.0],
            ResonanceOptions::default(),
            |w| {
                assert!(frequencies.insert(w.to_bits()), "repeated expensive solve");
                Ok(mobility(w).scale(gain))
            },
            || false,
        )
        .unwrap();
        assert!((result.angular_frequency_rad_s - natural).abs() < 1e-6 * damping);
        assert!(
            (result.bandwidth_quality_factor / quality - 1.0).abs() < 3e-6,
            "gain={gain:e}: Q={}",
            result.bandwidth_quality_factor
        );
        assert!(
            result.peak_bracket_rad_s[1] - result.peak_bracket_rad_s[0]
                <= 1e-6 * result.bandwidth_rad_s
        );
        assert_eq!(result.evaluations, frequencies.len());
        assert!(
            result.evaluations < 200,
            "{} evaluations",
            result.evaluations
        );
    }
}

#[test]
fn actual_viscothermal_cylinder_resolves_physical_width_on_a_tolerance_ladder() {
    use crate::{Duct, LossModel, Segment, Termination, modal::PreparedModalDuct};
    use fs_material::gas::{GasSpec, GasState};
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap();
    let (radius, length) = (0.0075, 0.5);
    let duct = Duct {
        segments: vec![Segment::Cylinder { radius, length }],
    };
    let prepared = PreparedModalDuct::new(&duct, 1, 1).unwrap();
    let natural = core::f64::consts::PI * gas.sound_speed / length;
    let mut previous: Option<ResponseResonance> = None;
    for tolerance in [1e-4, 1e-6, 1e-8] {
        let response = |w| {
            prepared
                .input_impedance(&gas, w, LossModel::WideTube, Termination::Closed)
                .map(|r| r.plane_impedance)
        };
        let result = refine_response_resonance(
            [0.9 * natural, 1.1 * natural],
            ResonanceOptions {
                relative_frequency_tolerance: tolerance,
                relative_bandwidth_tolerance: tolerance,
                max_evaluations: 256,
            },
            response,
            || false,
        )
        .unwrap();
        let w = result.angular_frequency_rad_s;
        let nu = gas.dynamic_viscosity / gas.density;
        let alpha = (nu * w / 2.0).sqrt() / (radius * gas.sound_speed)
            * (1.0 + (gas.gamma - 1.0) / gas.prandtl.sqrt());
        let kirchhoff_q = (w / gas.sound_speed) / (2.0 * alpha);
        assert!((result.bandwidth_quality_factor / kirchhoff_q - 1.0).abs() < 0.03);
        assert!(w < natural, "physical viscothermal dispersion lowers pitch");
        for &crossing in &result.half_power_frequencies_rad_s {
            let ratio = response(crossing).unwrap().abs() / result.response.abs();
            assert!((ratio * core::f64::consts::SQRT_2 - 1.0).abs() < 2.0 * tolerance);
        }
        if let Some(before) = previous {
            assert!(
                result.peak_bracket_rad_s[1] - result.peak_bracket_rad_s[0]
                    < before.peak_bracket_rad_s[1] - before.peak_bracket_rad_s[0]
            );
            assert!(
                (result.bandwidth_quality_factor / before.bandwidth_quality_factor - 1.0).abs()
                    < 2e-4
            );
        }
        assert!(result.evaluations < 200);
        previous = Some(result);
    }
}

#[test]
fn missing_peak_crossings_and_observed_overlapping_peaks_refuse() {
    let options = ResonanceOptions::default();
    assert_eq!(
        refine_response_resonance([1.0, 2.0], options, |w| Ok(C64::from_re(w)), || false),
        Err(ResonanceError::PeakNotBracketed)
    );
    assert_eq!(
        refine_response_resonance([1.0, 2.0], options, |_| Ok(C64::ZERO), || false),
        Err(ResonanceError::PeakNotBracketed)
    );
    let clipped = refine_response_resonance(
        [99.9, 100.1],
        options,
        |w| Ok(receptance(w, 100.0, 0.01)),
        || false,
    );
    assert!(matches!(
        clipped,
        Err(ResonanceError::HalfPowerNotBracketed { .. })
    ));
    let overlap = refine_response_resonance(
        [1.0, 2.0],
        options,
        |w| {
            let a = (w - 1.1) / 0.12;
            let b = (w - 1.8) / 0.15;
            Ok(C64::from_re(
                0.2 + 2.0 * det::exp(-a * a) + 10.0 * det::exp(-b * b),
            ))
        },
        || false,
    );
    assert_eq!(overlap, Err(ResonanceError::NotIsolated));
}

#[test]
fn lossless_pole_never_acquires_a_finite_bandwidth_quality_factor() {
    let result = refine_response_resonance(
        [1000.0, 1500.0],
        ResonanceOptions::default(),
        |w| Ok(receptance(w, 1234.56789123, 0.0)),
        || false,
    );
    assert!(
        matches!(
            result,
            Err(ResonanceError::ResolutionLimit
                | ResonanceError::NonFiniteResponse { .. }
                | ResonanceError::EvaluationLimit { .. })
        ),
        "lossless singularity must remain unresolved: {result:?}"
    );
}

#[test]
fn adjacent_cached_equal_samples_cannot_spin_past_the_work_budget() {
    let mut samples = Samples {
        evaluate: |_| Ok(C64::ONE),
        cancelled: || false,
        limit: 8,
        values: BTreeMap::new(),
    };
    let left = samples.at(1.0).unwrap();
    let right = samples.at(1.0_f64.next_up()).unwrap();
    let mut maximum = Maximum {
        bracket: [left.omega, right.omega],
        left,
        right,
    };
    assert!(matches!(
        maximum.refine(&mut samples, 1e-8, 1e-20),
        Err(ResonanceError::ResolutionLimit)
    ));
    assert_eq!(samples.values.len(), 2);
}

#[test]
fn limits_cancellation_and_original_model_failures_never_publish_partial_results() {
    let calls = Cell::new(0);
    let response = |w| {
        calls.set(calls.get() + 1);
        Ok(receptance(w, 100.0, 0.02))
    };
    let limited = ResonanceOptions {
        max_evaluations: 7,
        ..ResonanceOptions::default()
    };
    assert_eq!(
        refine_response_resonance([80.0, 120.0], limited, response, || false),
        Err(ResonanceError::EvaluationLimit { evaluations: 7 })
    );
    assert_eq!(calls.get(), 7);
    calls.set(0);
    assert_eq!(
        refine_response_resonance(
            [80.0, 120.0],
            ResonanceOptions::default(),
            response,
            || calls.get() >= 7
        ),
        Err(ResonanceError::Cancelled { evaluations: 7 })
    );
    assert_eq!(calls.get(), 7);
    let original = DuctError::OutOfTable {
        omega: 80.0,
        omega_lo: 90.0,
        omega_hi: 110.0,
    };
    assert_eq!(
        refine_response_resonance(
            [80.0, 120.0],
            ResonanceOptions::default(),
            |_| Err(original.clone()),
            || false
        ),
        Err(ResonanceError::Response {
            omega: 80.0,
            source: original
        })
    );
    for interval in [
        [0.0, 1.0],
        [2.0, 1.0],
        [1.0, f64::INFINITY],
        [f64::NAN, 2.0],
    ] {
        assert!(matches!(
            refine_response_resonance(
                interval,
                ResonanceOptions::default(),
                |_| panic!("invalid bounds must not evaluate the model"),
                || false
            ),
            Err(ResonanceError::InvalidInput { .. })
        ));
    }
}
