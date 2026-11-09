//! G1: callback queries retain the actual damped mechanical line shape.

use fs_math::c64::C64;
use fs_modal::harmonic::resonance::{
    HarmonicTransferSample, ResonanceControls, ResonanceError, TransferExtremumKind,
    refine_half_power_bandwidth_with_cancel, refine_transfer_extremum_with_cancel,
};
use std::cell::{Cell, RefCell};
use std::collections::BTreeSet;

const MASS: f64 = 2.0;
const STIFFNESS: f64 = 50.0;
const DAMPING: f64 = 5.0;

fn controls() -> ResonanceControls {
    ResonanceControls {
        absolute_angular_tolerance_rad_s: 1e-10,
        relative_angular_tolerance: 1e-12,
        maximum_evaluations: 256,
    }
}

/// The complete scalar equation under exp(-i*w*t), with its exact derivative.
/// No frequency differences, pole fit, modal Q or search result enters it.
fn displacement(
    angular_frequency: f64,
    force: C64,
) -> Result<HarmonicTransferSample, &'static str> {
    let denominator = C64::new(
        STIFFNESS - MASS * angular_frequency * angular_frequency,
        -DAMPING * angular_frequency,
    );
    if denominator == C64::ZERO {
        return Err("singular original mechanical equation");
    }
    Ok(HarmonicTransferSample {
        amplitude: force / denominator,
        angular_derivative: force * C64::new(2.0 * MASS * angular_frequency, DAMPING)
            / (denominator * denominator),
    })
}

fn exact_peak_and_crossings() -> (f64, f64, f64) {
    // Completing the square in s=w^2 gives the displacement maximum and the
    // two roots of |K-M*w^2-i*C*w|^2 = 2*minimum. These crossings are asymmetric
    // around the actual displacement peak and do not have linewidth C/M.
    let rate = DAMPING / MASS;
    let natural_square = STIFFNESS / MASS;
    let peak_square = natural_square - 0.5 * rate * rate;
    let spread = rate * (natural_square - 0.25 * rate * rate).sqrt();
    (
        peak_square.sqrt(),
        (peak_square - spread).sqrt(),
        (peak_square + spread).sqrt(),
    )
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 2e-9 * (1.0 + expected.abs()),
        "actual {actual:.16e}, expected {expected:.16e}"
    );
}

#[test]
fn displacement_peak_retains_signed_complex_phase_and_the_original_frequency_tangent() {
    let (expected_peak, _, _) = exact_peak_and_crossings();
    for force in [C64::new(0.6, 0.8), C64::new(-0.6, -0.8)] {
        let peak = refine_transfer_extremum_with_cancel(
            (4.0, 5.1),
            controls(),
            || false,
            |w| displacement(w, force),
        )
        .unwrap();
        assert_eq!(peak.kind, TransferExtremumKind::Maximum);
        close(peak.angular_frequency_rad_s, expected_peak);
        assert!(peak.bracket_rad_s.0 <= expected_peak);
        assert!(expected_peak <= peak.bracket_rad_s.1);
        assert!(peak.evaluations <= controls().maximum_evaluations);
        let expected = displacement(peak.angular_frequency_rad_s, force).unwrap();
        close(peak.sample.amplitude.re, expected.amplitude.re);
        close(peak.sample.amplitude.im, expected.amplitude.im);
        close(
            peak.sample.angular_derivative.re,
            expected.angular_derivative.re,
        );
        close(
            peak.sample.angular_derivative.im,
            expected.angular_derivative.im,
        );
        // The actual displacement maximum is below the undamped natural
        // frequency; the returned complex transfer retains the force sign.
        assert!(peak.angular_frequency_rad_s < (STIFFNESS / MASS).sqrt() - 0.2);
        assert_eq!(
            peak.sample.amplitude.re.is_sign_positive(),
            force.re.is_sign_negative()
        );
    }
}

#[test]
fn true_asymmetric_half_power_crossings_determine_the_response_quality_factor() {
    let force = C64::new(-0.8, 0.6);
    let peak = refine_transfer_extremum_with_cancel(
        (4.0, 5.1),
        controls(),
        || false,
        |w| displacement(w, force),
    )
    .unwrap();
    let bandwidth = refine_half_power_bandwidth_with_cancel(
        &peak,
        (2.0, 4.0),
        (5.0, 7.0),
        controls(),
        || false,
        |w| displacement(w, force),
    )
    .unwrap();
    let (expected_peak, expected_lower, expected_upper) = exact_peak_and_crossings();
    let expected_width = expected_upper - expected_lower;
    close(bandwidth.lower_angular_frequency_rad_s, expected_lower);
    close(bandwidth.upper_angular_frequency_rad_s, expected_upper);
    close(bandwidth.angular_bandwidth_rad_s, expected_width);
    close(
        bandwidth.response_quality_factor,
        expected_peak / expected_width,
    );
    assert!(bandwidth.lower_bracket_rad_s.0 <= expected_lower);
    assert!(expected_lower <= bandwidth.lower_bracket_rad_s.1);
    assert!(bandwidth.upper_bracket_rad_s.0 <= expected_upper);
    assert!(expected_upper <= bandwidth.upper_bracket_rad_s.1);
    for crossing in [
        bandwidth.lower_angular_frequency_rad_s,
        bandwidth.upper_angular_frequency_rad_s,
    ] {
        let measured = displacement(crossing, force).unwrap().amplitude.norm_sq();
        close(measured / peak.sample.amplitude.norm_sq(), 0.5);
    }
    let lower_distance = expected_peak - expected_lower;
    let upper_distance = expected_upper - expected_peak;
    assert!(lower_distance > upper_distance + 0.2);
    assert!((bandwidth.angular_bandwidth_rad_s - DAMPING / MASS).abs() > 0.1);
    let modal_quality = (MASS * STIFFNESS).sqrt() / DAMPING;
    assert!((bandwidth.response_quality_factor - modal_quality).abs() > 0.1);
}

#[test]
fn stationary_resonance_reuses_its_physical_solve_within_a_tight_query_budget() {
    // The exact response of a damped first-order mode in a rotating frame:
    // (gamma - i*(omega-center))*H = gamma. Its peak and all bisection
    // midpoints are exactly representable. Each callback may be a full coupled
    // factorization, so resolving its two-sided bracket must not repeatedly
    // solve the unchanged peak.
    let calls = RefCell::new(BTreeSet::new());
    let tolerance = 2.0_f64.powi(-34);
    let result = refine_transfer_extremum_with_cancel(
        (6.0, 10.0),
        ResonanceControls {
            absolute_angular_tolerance_rad_s: tolerance,
            relative_angular_tolerance: 0.0,
            maximum_evaluations: 76,
        },
        || false,
        |omega| {
            assert!(
                calls.borrow_mut().insert(omega.to_bits()),
                "repeated original physical solve at {omega}"
            );
            let denominator = C64::new(0.5, 8.0 - omega);
            let amplitude = C64::from_re(0.5) / denominator;
            Ok::<_, &'static str>(HarmonicTransferSample {
                amplitude,
                angular_derivative: C64::new(0.0, 1.0) * amplitude / denominator,
            })
        },
    )
    .unwrap();
    assert_eq!(result.angular_frequency_rad_s, 8.0);
    assert_eq!(result.sample.amplitude, C64::ONE);
    assert!(result.bracket_rad_s.0 < 8.0 && result.bracket_rad_s.1 > 8.0);
    assert!(result.bracket_rad_s.1 - result.bracket_rad_s.0 <= tolerance);
    assert_eq!(result.evaluations, calls.borrow().len());
    assert!(result.evaluations <= 76);
}

#[test]
fn callback_queries_preserve_cancellation_evaluation_caps_and_original_refusals() {
    let force = C64::new(0.6, 0.8);
    let calls = Cell::new(0);
    let immediate = refine_transfer_extremum_with_cancel(
        (4.0, 5.1),
        controls(),
        || true,
        |w| {
            calls.set(calls.get() + 1);
            displacement(w, force)
        },
    );
    assert!(matches!(immediate, Err(ResonanceError::Cancelled)));
    assert_eq!(calls.get(), 0);

    let requested = Cell::new(false);
    let during_evaluation = refine_transfer_extremum_with_cancel(
        (4.0, 5.1),
        controls(),
        || requested.get(),
        |w| {
            calls.set(calls.get() + 1);
            requested.set(true);
            displacement(w, force)
        },
    );
    assert!(matches!(during_evaluation, Err(ResonanceError::Cancelled)));
    assert_eq!(calls.get(), 1);

    calls.set(0);
    let mut limited = controls();
    limited.maximum_evaluations = 3;
    let exhausted = refine_transfer_extremum_with_cancel(
        (4.0, 5.1),
        limited,
        || false,
        |w| {
            calls.set(calls.get() + 1);
            displacement(w, force)
        },
    );
    assert!(matches!(
        exhausted,
        Err(ResonanceError::EvaluationBudget {
            maximum_evaluations: 3,
        })
    ));
    assert_eq!(calls.get(), 3);

    let refused = refine_transfer_extremum_with_cancel(
        (4.0, 5.1),
        controls(),
        || false,
        |_| Err::<HarmonicTransferSample, _>("original physical source refused"),
    );
    assert!(matches!(
        refused,
        Err(ResonanceError::Evaluation(
            "original physical source refused"
        ))
    ));

    let peak = refine_transfer_extremum_with_cancel(
        (4.0, 5.1),
        controls(),
        || false,
        |w| displacement(w, force),
    )
    .unwrap();
    assert!(matches!(
        refine_half_power_bandwidth_with_cancel(
            &peak,
            (2.0, 4.0),
            (5.0, 7.0),
            controls(),
            || true,
            |w| displacement(w, force),
        ),
        Err(ResonanceError::Cancelled)
    ));
}
