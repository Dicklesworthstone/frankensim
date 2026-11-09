//! G1/G3 complete work chain rules and scaled signed physical contraction.
use fs_math::c64::C64;
use fs_modal::harmonic::resonance::HarmonicTransferFieldParameterJet;
use fs_modal::harmonic::{HarmonicResponseError, harmonic_mean_power_parameter_jet_with_cancel};

fn field(values: &[[C64; 1]; 5]) -> HarmonicTransferFieldParameterJet<'_> {
    HarmonicTransferFieldParameterJet {
        values: &values[0],
        angular_derivative: &values[1],
        angular_second_derivative: &values[2],
        parameter_derivative: &values[3],
        angular_parameter_derivative: &values[4],
    }
}

#[test]
fn g1_resistive_heat_retains_source_and_material_directions_without_reactive_heat() {
    let (w, g, tau, resistance, inertance) = (1.9, 0.7, 0.8, 1.3, 2.0);
    let (dg, dtau, dr, dl) = (-0.1, 0.3, 0.2, -0.4);
    // u=g*(1+i*w*tau), e=(R+i*w*L)*u. Independent physical heat is
    // P=R*g²*(1+w²*tau²)/2, with no contribution from the inertance.
    let u = C64::new(g, g * w * tau);
    let uw = C64::new(0.0, g * tau);
    let ut = C64::new(dg, w * (dg * tau + g * dtau));
    let uwt = C64::new(0.0, dg * tau + g * dtau);
    let z = C64::new(resistance, w * inertance);
    let zw = C64::new(0.0, inertance);
    let zt = C64::new(dr, w * dl);
    let zwt = C64::new(0.0, dl);
    let flow = [[u], [uw], [C64::ZERO], [ut], [uwt]];
    let effort = [
        [z * u],
        [zw * u + z * uw],
        [(zw * uw).scale(2.0)],
        [zt * u + z * ut],
        [zwt * u + zw * ut + zt * uw + z * uwt],
    ];
    let actual =
        harmonic_mean_power_parameter_jet_with_cancel(field(&flow), field(&effort), || false)
            .unwrap();
    let amplitude = 0.5 * resistance * g * g;
    let amplitude_t = 0.5 * dr * g * g + resistance * g * dg;
    let expected = [
        amplitude * (1.0 + w * w * tau * tau),
        2.0 * amplitude * w * tau * tau,
        2.0 * amplitude * tau * tau,
        amplitude_t * (1.0 + w * w * tau * tau) + 2.0 * amplitude * w * w * tau * dtau,
        2.0 * amplitude_t * w * tau * tau + 4.0 * amplitude * w * tau * dtau,
    ];
    for (a, b) in [
        actual.sample.average_power_w,
        actual.sample.angular_derivative,
        actual.angular_second_derivative,
        actual.parameter_derivative,
        actual.angular_parameter_derivative,
    ]
    .into_iter()
    .zip(expected)
    {
        assert!((a - b).abs() < 3e-14 * b.abs().max(1.0), "{a} != {b}");
    }
}

#[test]
fn g3_equal_and_opposite_large_work_and_reactive_derivatives_cancel_before_range_recovery() {
    let huge = 2.0_f64.powi(700);
    let values = [C64::from_re(huge); 3];
    let effort = [
        C64::from_re(huge),
        C64::from_re(-huge),
        C64::from_re(2.0 / huge),
    ];
    let zero = [C64::ZERO; 3];
    let flow = HarmonicTransferFieldParameterJet {
        values: &values,
        angular_derivative: &values,
        angular_second_derivative: &values,
        parameter_derivative: &values,
        angular_parameter_derivative: &values,
    };
    let effort = HarmonicTransferFieldParameterJet {
        values: &effort,
        angular_derivative: &zero,
        angular_second_derivative: &zero,
        parameter_derivative: &zero,
        angular_parameter_derivative: &zero,
    };
    let actual = harmonic_mean_power_parameter_jet_with_cancel(flow, effort, || false).unwrap();
    assert_eq!(
        [
            actual.sample.average_power_w,
            actual.sample.angular_derivative,
            actual.angular_second_derivative,
            actual.parameter_derivative,
            actual.angular_parameter_derivative,
        ],
        [1.0; 5]
    );
    assert!(matches!(
        harmonic_mean_power_parameter_jet_with_cancel(flow, effort, || true),
        Err(HarmonicResponseError::Cancelled)
    ));
    let incomplete = HarmonicTransferFieldParameterJet {
        parameter_derivative: &zero[..2],
        ..flow
    };
    assert!(harmonic_mean_power_parameter_jet_with_cancel(incomplete, effort, || false).is_err());
}
