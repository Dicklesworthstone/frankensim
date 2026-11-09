//! First-transfer observations retain finite coherent sums before raw overflow.
use fs_math::c64::C64;
use fs_modal::harmonic::HarmonicResponseError;
use fs_modal::harmonic::resonance::{
    HarmonicReceiverParameterJet, HarmonicTransferFieldParameterJet, HarmonicTransferParameterJet,
    HarmonicTransferSample,
};

#[test]
fn real_receiver_retains_weak_phasors_between_individually_overflowing_products() {
    let weight = 2.0_f64.powi(600);
    let field = 2.0_f64.powi(450);
    let weak = 2.0_f64.powi(950);
    assert!((weight * field).is_infinite());
    let sample = HarmonicTransferSample::project_real(
        &[weight, 1.0, weight],
        &[
            C64::new(field, -field),
            C64::new(weak, -0.5 * weak),
            C64::new(-field, field),
        ],
        &[
            C64::new(field, field),
            C64::new(0.25 * weak, 0.125 * weak),
            C64::new(-field, -field),
        ],
    )
    .unwrap();
    // The weak coherent term occurs BETWEEN the cancelling large terms.
    // Its phase and derivative survive without first forming weight*field.
    assert_eq!(sample.amplitude, C64::new(weak, -0.5 * weak));
    assert_eq!(
        sample.angular_derivative,
        C64::new(0.25 * weak, 0.125 * weak),
    );
}

#[test]
fn complex_receiver_and_receiver_slope_match_the_complete_jet_through_overflowing_cancellation() {
    let weight = 2.0_f64.powi(600);
    let field = 2.0_f64.powi(450);
    let weak = 2.0_f64.powi(950);
    let receiver = [
        C64::new(weight, weight),
        C64::new(1.0, 2.0),
        C64::new(-weight, -weight),
    ];
    let receiver_slope = [
        C64::new(-weight, -weight),
        C64::new(3.0, -4.0),
        C64::new(weight, weight),
    ];
    let values = [
        C64::new(field, 0.5 * field),
        C64::new(weak, -0.5 * weak),
        C64::new(field, 0.5 * field),
    ];
    let derivatives = [
        C64::new(field, 0.5 * field),
        C64::new(-0.25 * weak, 0.125 * weak),
        C64::new(field, 0.5 * field),
    ];
    let fixed = HarmonicTransferSample::project_complex(&receiver, &values, &derivatives).unwrap();
    // (1+2i)*(1-i/2) = 2+3i/2: the receiver is NOT conjugated.
    assert_eq!(fixed.amplitude, C64::new(2.0 * weak, 1.5 * weak));
    assert_eq!(
        fixed.angular_derivative,
        C64::new(-0.5 * weak, -0.375 * weak),
    );

    let sample = HarmonicTransferSample::project_complex_with_derivative(
        &receiver,
        &receiver_slope,
        &values,
        &derivatives,
    )
    .unwrap();
    // Both r*x_omega and r_omega*x overflow separately in the first port,
    // while their complete derivative is finite. The middle port contributes
    // (-1/2-3i/8)+(1-11i/2) = 1/2-47i/8.
    assert_eq!(sample.amplitude, fixed.amplitude);
    assert_eq!(
        sample.angular_derivative,
        C64::new(0.5 * weak, -5.875 * weak),
    );

    let zero = [C64::ZERO; 3];
    let complete = HarmonicTransferParameterJet::project_complex_with_derivatives(
        HarmonicReceiverParameterJet {
            weights: &receiver,
            angular_derivative: Some(&receiver_slope),
            angular_second_derivative: None,
            parameter_derivative: None,
            angular_parameter_derivative: None,
        },
        HarmonicTransferFieldParameterJet {
            values: &values,
            angular_derivative: &derivatives,
            angular_second_derivative: &zero,
            parameter_derivative: &zero,
            angular_parameter_derivative: &zero,
        },
    )
    .unwrap();
    assert_eq!(sample, complete.sample);
}

#[test]
fn subnormal_receiver_products_round_after_the_complete_coherent_sum() {
    let tiny = f64::from_bits(1);
    assert_eq!(0.5 * tiny, 0.0);
    let sample = HarmonicTransferSample::project_real(
        &[tiny, tiny],
        &[C64::new(0.5, -0.5); 2],
        &[C64::new(1.5, -1.5); 2],
    )
    .unwrap();
    assert_eq!(sample.amplitude, C64::new(tiny, -tiny));
    // Two exact 1.5*tiny contributions sum to 3*tiny; rounding each first
    // would instead give 4*tiny under round-to-nearest ties-to-even.
    assert_eq!(
        sample.angular_derivative,
        C64::new(f64::from_bits(3), -f64::from_bits(3)),
    );

    let receiver = [C64::new(tiny, tiny); 2];
    let field = [C64::new(0.5, 0.5); 2];
    let complex = HarmonicTransferSample::project_complex_with_derivative(
        &receiver, &receiver, &field, &field,
    )
    .unwrap();
    assert_eq!(complex.amplitude, C64::new(0.0, f64::from_bits(2)));
    assert_eq!(complex.angular_derivative, C64::new(0.0, f64::from_bits(4)),);
}

#[test]
fn receiver_projection_refuses_incomplete_nonfinite_and_unrepresentable_observations() {
    let real = |r: &[f64], x: &[C64], dx: &[C64]| {
        matches!(
            HarmonicTransferSample::project_real(r, x, dx),
            Err(HarmonicResponseError::Invalid(_))
        )
    };
    let complex = |r: &[C64], dr: &[C64], x: &[C64], dx: &[C64]| {
        matches!(
            HarmonicTransferSample::project_complex_with_derivative(r, dr, x, dx),
            Err(HarmonicResponseError::Invalid(_))
        )
    };
    assert!(real(&[], &[], &[]));
    assert!(real(&[], &[C64::ONE], &[C64::ONE]));
    assert!(real(&[1.0], &[C64::ONE], &[]));
    assert!(real(&[f64::NAN], &[C64::ONE], &[C64::ONE]));
    assert!(real(&[0.0], &[C64::new(f64::INFINITY, 0.0)], &[C64::ONE]));
    assert!(real(&[0.0], &[C64::ONE], &[C64::new(0.0, f64::NAN)]));
    assert!(complex(&[], &[], &[], &[]));
    assert!(complex(&[], &[C64::ZERO], &[C64::ONE], &[C64::ONE]));
    assert!(complex(&[C64::ONE], &[], &[C64::ONE], &[C64::ONE]));
    assert!(complex(&[C64::ONE], &[C64::ZERO], &[C64::ONE], &[]));
    assert!(complex(
        &[C64::new(f64::NAN, 0.0)],
        &[C64::ZERO],
        &[C64::ONE],
        &[C64::ONE],
    ));
    assert!(complex(
        &[C64::ONE],
        &[C64::new(0.0, f64::INFINITY)],
        &[C64::ONE],
        &[C64::ONE],
    ));
    assert!(complex(
        &[C64::ZERO],
        &[C64::ZERO],
        &[C64::new(f64::NAN, 0.0)],
        &[C64::ONE],
    ));
    assert!(complex(
        &[C64::ZERO],
        &[C64::ZERO],
        &[C64::ONE],
        &[C64::new(0.0, f64::INFINITY)],
    ));

    let weight = 2.0_f64.powi(600);
    let field = C64::from_re(2.0_f64.powi(450));
    assert!(real(&[weight], &[field], &[C64::ZERO]));
    assert!(complex(
        &[C64::ZERO],
        &[C64::from_re(weight)],
        &[field],
        &[C64::ZERO],
    ));
}
