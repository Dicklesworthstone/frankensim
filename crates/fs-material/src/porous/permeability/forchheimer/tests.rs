use super::*;

fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 8e-13 * actual.abs().max(expected.abs()),
        "{actual} != {expected}"
    );
}

#[test]
fn manufactured_quadratic_balance_across_viscous_and_inertial_regimes() {
    for (linear, quadratic, expected) in [
        (3.0, 2.0, 7.0),
        (3.0, 0.0, 7.0),
        (1e200, 1e200, 0.5),
        (1e-200, 1e-200, 1e100),
        (1.0, 1e-200, 1.0),
        (1e-300, 1.0, 1e150),
    ] {
        let driving = linear * expected + (quadratic * expected) * expected;
        let solved = positive_quadratic_flow(linear, quadratic, driving).unwrap();
        near(solved.flow, expected);
        let (viscous, inertial) = loss_fractions(linear, quadratic, solved.flow).unwrap();
        near(viscous + inertial, 1.0);
        let reconstructed =
            (linear * solved.flow) / driving + ((quadratic * solved.flow) * solved.flow) / driving;
        near(reconstructed, 1.0);
    }
}

#[test]
fn actual_inertial_conductance_survives_unopposed_viscous_overflow() {
    let flow = positive_quadratic_flow(1e-310, 1e4, 1.0).unwrap().flow;
    near(
        secant_volume_conductance(1e-310, 1e4, flow, 1.0).unwrap(),
        0.01,
    );
    near(
        secant_volume_conductance(f64::MAX, f64::MAX, 1.0, 1.0).unwrap(),
        0.5 / f64::MAX,
    );
    near(
        secant_volume_conductance(0.75, 0.75, 1.0, f64::MAX).unwrap(),
        f64::MAX / 1.5,
    );
}

#[test]
fn quadratic_root_avoids_overflow_and_cancellation() {
    // The unopposed Darcy solution overflows; inertial resistance makes the
    // actual flow representable. a² and b*c also need not be representable.
    let solved = positive_quadratic_flow(1e-300, 1.0, 1e300).unwrap();
    near(solved.flow, 1e150);
    let solved = positive_quadratic_flow(1e200, 1e200, 1e200).unwrap();
    near(solved.flow * solved.flow + solved.flow, 1.0);
    let tiny = positive_quadratic_flow(1.0, 1.0, 1e-300).unwrap();
    near(tiny.flow, 1e-300);
    let (viscous, inertial) = loss_fractions(f64::MAX, f64::MAX, 1.0).unwrap();
    assert_eq!((viscous, inertial), (0.5, 0.5));
    for (a, b, c) in [
        (0.0, 1.0, 1.0),
        (1.0, -1.0, 1.0),
        (1.0, 1.0, 0.0),
        (f64::INFINITY, 1.0, 1.0),
        (1.0, f64::NAN, 1.0),
    ] {
        assert!(positive_quadratic_flow(a, b, c).is_err());
    }
}

#[test]
fn tensor_admission_does_not_hide_transverse_flow_behind_unrelated_contrast() {
    let s = core::f64::consts::FRAC_1_SQRT_2;
    let axes = [[s, 0.0, s], [0.0, 1.0, 0.0], [-s, 0.0, s]];
    let tensor = IntrinsicPermeability::orthotropic([1.0, 1e-300, 2.0], axes).unwrap();
    assert!(!is_eigendirection(&tensor, [1.0, 0.0, 0.0]).unwrap());
    assert!(is_eigendirection(&tensor, [s, 0.0, s]).unwrap());
    let isotropic = IntrinsicPermeability::orthotropic([1e-300; 3], axes).unwrap();
    assert!(is_eigendirection(&isotropic, [1.0, 0.0, 0.0]).unwrap());
    let high_contrast = IntrinsicPermeability::orthotropic(
        [1.0, 1e300, 1.0],
        [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
    )
    .unwrap();
    // A small geometric projection onto an extremely permeable direction can
    // produce dominant transverse flow; normalizing by n alone would miss it.
    assert!(!is_eigendirection(&high_contrast, [1.0, 1e-150, 0.0]).unwrap());
}
