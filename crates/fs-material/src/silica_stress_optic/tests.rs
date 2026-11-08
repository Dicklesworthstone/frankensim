use super::*;

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}"
    );
}

#[test]
fn original_three_lines_retain_si_coefficient_and_specimen_sigma() {
    for (wavelength, mean, sigma) in [
        (193.0, 5.1e-12, 3.0e-13),
        (248.0, 4.2e-12, 1.0e-13),
        (633.0, 3.5e-12, 1.0e-13),
    ] {
        let law = Corning7980StressOptics::at_source_wavelength(wavelength).unwrap();
        assert_eq!(law.wavelength_nm(), wavelength);
        near(law.coefficient_per_pa(), mean, 1.0e-27);
        near(law.specimen_standard_deviation_per_pa(), sigma, 1.0e-28);
        // Original nm/cm/MPa units: a 1 MPa difference over 1 cm produces
        // 35, 42 or 51 nm of differential optical path, not metres.
        let sample = law
            .project(
                &[1.0e6, 0.0, 0.0, 0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
            )
            .unwrap();
        near(
            sample.horizontal_minus_vertical() * 0.01,
            mean * 1.0e4,
            3.0e-23,
        );
    }
    let identity = Corning7980StressOptics::source_identity();
    assert_ne!(identity, ContentHash([0; 32]));
    assert_eq!(identity, Corning7980StressOptics::source_identity());
    for wavelength in [
        0.0,
        -1.0,
        192.9,
        200.0,
        550.0,
        632.8,
        633.001,
        f64::NAN,
        f64::INFINITY,
    ] {
        assert_eq!(
            Corning7980StressOptics::at_source_wavelength(wavelength),
            Err(SilicaStressOpticError::UnsupportedSourceWavelength)
        );
    }
}

#[test]
fn principal_stress_shear_and_oblique_projection_follow_original_voigt_convention() {
    let law = Corning7980StressOptics::at_source_wavelength(633.0).unwrap();
    let x = [1.0, 0.0, 0.0];
    let y = [0.0, 1.0, 0.0];
    let length = 633.0e-9 / (4.0 * 3.5e-12 * 400_000.0);
    for sign in [-1.0, 1.0] {
        let sample = law
            .project(&[sign * 400_000.0, 0.0, 0.0, 0.0, 0.0, 0.0], x, y)
            .unwrap();
        near(sample.horizontal_minus_vertical(), sign * 1.4e-6, 3.0e-22);
        assert_eq!(sample.twice_cross(), 0.0);
        let phase = sample.retardance(length).unwrap();
        near(phase.phase_rad, core::f64::consts::FRAC_PI_2, 5.0e-16);
        assert_eq!(phase.cos_double_slow_axis, sign);
        assert_eq!(phase.sin_double_slow_axis, 0.0);
    }
    // xy is stored ONCE as tensor shear. Its principal difference is 2*tau.
    let shear = law
        .project(&[0.0, 0.0, 0.0, 200_000.0, 0.0, 0.0], x, y)
        .unwrap();
    assert_eq!(shear.horizontal_minus_vertical(), 0.0);
    near(shear.twice_cross(), 1.4e-6, 3.0e-22);
    let phase = shear.retardance(length).unwrap();
    near(phase.phase_rad, core::f64::consts::FRAC_PI_2, 5.0e-16);
    assert_eq!(phase.cos_double_slow_axis, 0.0);
    assert_eq!(phase.sin_double_slow_axis, 1.0);

    // Independent tensor: sigma = pressure*I + tau*a*a^T, a=(1,2,2)/3.
    let (pressure, tau) = (80.0e6, 400_000.0);
    let stress = [
        pressure + tau / 9.0,
        pressure + 4.0 * tau / 9.0,
        pressure + 4.0 * tau / 9.0,
        2.0 * tau / 9.0,
        4.0 * tau / 9.0,
        2.0 * tau / 9.0,
    ];
    let h = [1.0 / 2.0_f64.sqrt(), 1.0 / 2.0_f64.sqrt(), 0.0];
    let v = [
        -1.0 / 6.0_f64.sqrt(),
        1.0 / 6.0_f64.sqrt(),
        2.0 / 6.0_f64.sqrt(),
    ];
    let original = law.project(&stress, h, v).unwrap();
    near(
        original.horizontal_minus_vertical(),
        3.5e-12 * tau / 27.0,
        8.0e-20,
    );
    near(
        original.twice_cross(),
        3.5e-12 * tau * 5.0 / (3.0 * 3.0_f64.sqrt()),
        8.0e-20,
    );
    let c = 3.0_f64.sqrt() / 2.0;
    let turned_h = core::array::from_fn(|i| c * h[i] + 0.5 * v[i]);
    let turned_v = core::array::from_fn(|i| -0.5 * h[i] + c * v[i]);
    let turned = law.project(&stress, turned_h, turned_v).unwrap();
    near(
        turned.horizontal_minus_vertical(),
        0.5 * original.horizontal_minus_vertical() + c * original.twice_cross(),
        1.0e-21,
    );
    near(
        turned.twice_cross(),
        -c * original.horizontal_minus_vertical() + 0.5 * original.twice_cross(),
        1.0e-21,
    );
    let hydro = law
        .project(&[pressure, pressure, pressure, 0.0, 0.0, 0.0], h, v)
        .unwrap();
    assert_eq!(hydro.horizontal_minus_vertical(), 0.0);
    assert_eq!(hydro.twice_cross(), 0.0);
    assert_eq!(hydro.retardance(0.4).unwrap().phase_rad, 0.0);
}

#[test]
fn invalid_stress_frame_and_path_refuse_before_creating_an_optical_state() {
    let law = Corning7980StressOptics::at_source_wavelength(633.0).unwrap();
    let x = [1.0, 0.0, 0.0];
    let y = [0.0, 1.0, 0.0];
    for value in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        let mut stress = [0.0; 6];
        stress[4] = value;
        assert_eq!(
            law.project(&stress, x, y),
            Err(SilicaStressOpticError::InvalidStress)
        );
    }
    for (h, v) in [
        (x, x),
        ([2.0, 0.0, 0.0], y),
        ([0.0; 3], y),
        (x, [0.0, f64::NAN, 0.0]),
    ] {
        assert_eq!(
            law.project(&[0.0; 6], h, v),
            Err(SilicaStressOpticError::InvalidTransverseFrame)
        );
    }
    let zero = law.project(&[0.0; 6], x, y).unwrap();
    for distance in [-1.0, f64::NAN, f64::INFINITY] {
        assert_eq!(
            zero.retardance(distance),
            Err(SilicaStressOpticError::InvalidDistance)
        );
    }
    assert_eq!(
        law.project(&[f64::MAX, 0.0, -f64::MAX, 0.0, 0.0, 0.0], x, y),
        Err(SilicaStressOpticError::NumericalRange)
    );
}
