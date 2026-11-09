//! Coherent receiver sensitivity, physical delays and cancellation before rounding.

use fs_math::c64::C64;
use fs_modal::harmonic::HarmonicResponseError;
use fs_modal::harmonic::resonance::{
    HarmonicReceiverParameterJet, HarmonicTransferFieldParameterJet, HarmonicTransferParameterJet,
    ResonanceControls, refine_transfer_extremum_with_parameter_derivative_with_cancel,
};

fn controls() -> ResonanceControls {
    ResonanceControls {
        absolute_angular_tolerance_rad_s: 1e-12,
        relative_angular_tolerance: 1e-12,
        maximum_evaluations: 200,
    }
}

fn close(actual: C64, expected: C64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1.0),
        "{actual:?} != {expected:?}",
    );
}

/// A damped oscillator observed after a physical propagation delay tau.
/// Theta is tau, so the mechanical field itself is independent of theta.
fn delayed_oscillator(
    w: f64,
    tau: f64,
) -> Result<HarmonicTransferParameterJet, HarmonicResponseError> {
    let d = C64::new(4.0 - w * w, -0.2 * w);
    let dw = C64::new(-2.0 * w, -0.2);
    let q = C64::ONE / d;
    let qw = -dw / (d * d);
    let qww = C64::from_re(2.0) / (d * d) + (dw * dw).scale(2.0) / (d * d * d);
    let (s, c) = (w * tau).sin_cos();
    let r = C64::new(c, s);
    let rw = C64::new(0.0, tau) * r;
    let rww = r.scale(-tau * tau);
    let rt = C64::new(0.0, w) * r;
    let rwt = C64::new(-w * tau, 1.0) * r;
    HarmonicTransferParameterJet::project_complex_with_derivatives(
        HarmonicReceiverParameterJet {
            weights: &[r],
            angular_derivative: Some(&[rw]),
            angular_second_derivative: Some(&[rww]),
            parameter_derivative: Some(&[rt]),
            angular_parameter_derivative: Some(&[rwt]),
        },
        HarmonicTransferFieldParameterJet {
            values: &[q],
            angular_derivative: &[qw],
            angular_second_derivative: &[qww],
            parameter_derivative: &[C64::ZERO],
            angular_parameter_derivative: &[C64::ZERO],
        },
    )
}

#[test]
fn physical_receiver_delay_changes_phase_without_moving_the_resonance() {
    let tau = 0.37;
    let response = refine_transfer_extremum_with_parameter_derivative_with_cancel(
        (1.7, 2.3),
        controls(),
        || false,
        |w| delayed_oscillator(w, tau).map(|jet| jet.sample),
        |w| delayed_oscillator(w, tau),
    )
    .unwrap();
    let w = response.extremum.angular_frequency_rad_s;
    assert!((w - (4.0_f64 - 0.2 * 0.2 / 2.0).sqrt()).abs() < 3e-12);
    assert!(response.angular_frequency_derivative_rad_s.abs() < 2e-13);
    assert!(response.log_amplitude_derivative.unwrap().abs() < 2e-13);
    assert!((response.phase_derivative.unwrap() - w).abs() < 2e-13);
    close(
        response.amplitude_derivative,
        C64::new(0.0, w) * response.extremum.sample.amplitude,
        3e-13,
    );
}

#[test]
fn mixed_receiver_terms_match_a_closed_delay_quotient() {
    for w in [0.0, 0.7, 1.99, 2.0, 5.0] {
        let tau = 0.37;
        let jet = delayed_oscillator(w, tau).unwrap();
        let d = C64::new(4.0 - w * w, -0.2 * w);
        let dw = C64::new(-2.0 * w, -0.2);
        let (s, c) = (w * tau).sin_cos();
        let r = C64::new(c, s);
        // Differentiate exp(i*w*tau)/(4-w*w-0.2*i*w) symbolically.
        let expected_mixed = r * (C64::new(-tau * w, 1.0) * d - C64::new(0.0, w) * dw) / (d * d);
        close(jet.angular_parameter_derivative, expected_mixed, 5e-13);
        close(jet.parameter_derivative, C64::new(0.0, w) * r / d, 3e-13);
    }
}

#[test]
fn finite_coherent_jets_survive_individually_overflowing_products() {
    let big = 1e200;
    let values = [C64::from_re(big), C64::from_re(big), C64::from_re(3.0)];
    let angular = [C64::from_re(big), C64::from_re(big), C64::from_re(5.0)];
    let second = [C64::from_re(big), C64::from_re(big), C64::from_re(7.0)];
    let parameter = [C64::from_re(big), C64::from_re(big), C64::from_re(11.0)];
    let mixed = [C64::from_re(big), C64::from_re(big), C64::from_re(13.0)];
    let field = HarmonicTransferFieldParameterJet {
        values: &values,
        angular_derivative: &angular,
        angular_second_derivative: &second,
        parameter_derivative: &parameter,
        angular_parameter_derivative: &mixed,
    };
    let real = HarmonicTransferParameterJet::project_real(&[big, -big, 2.0], field).unwrap();
    assert_eq!(real.sample.amplitude, C64::from_re(6.0));
    assert_eq!(real.sample.angular_derivative, C64::from_re(10.0));
    assert_eq!(real.angular_second_derivative, C64::from_re(14.0));
    assert_eq!(real.parameter_derivative, C64::from_re(22.0));
    assert_eq!(real.angular_parameter_derivative, C64::from_re(26.0));
    let complex = HarmonicTransferParameterJet::project_complex(
        &[C64::new(big, big), C64::new(-big, -big), C64::new(2.0, 2.0)],
        field,
    )
    .unwrap();
    assert_eq!(complex.sample.amplitude, C64::new(6.0, 6.0));
    assert_eq!(complex.angular_parameter_derivative, C64::new(26.0, 26.0));
}

#[test]
fn invalid_receiver_and_field_jets_refuse_before_publication() {
    let zero = [C64::ZERO];
    let field = HarmonicTransferFieldParameterJet {
        values: &zero,
        angular_derivative: &zero,
        angular_second_derivative: &zero,
        parameter_derivative: &zero,
        angular_parameter_derivative: &zero,
    };
    assert!(HarmonicTransferParameterJet::project_real(&[], field).is_err());
    assert!(HarmonicTransferParameterJet::project_real(&[f64::NAN], field).is_err());
    assert!(
        HarmonicTransferParameterJet::project_complex_with_derivatives(
            HarmonicReceiverParameterJet {
                weights: &[C64::ONE],
                angular_derivative: None,
                angular_second_derivative: None,
                parameter_derivative: None,
                angular_parameter_derivative: Some(&[]),
            },
            field,
        )
        .is_err()
    );
    assert!(
        HarmonicTransferParameterJet::project_real(
            &[1.0],
            HarmonicTransferFieldParameterJet {
                angular_parameter_derivative: &[C64::new(f64::INFINITY, 0.0)],
                ..field
            },
        )
        .is_err()
    );
}
