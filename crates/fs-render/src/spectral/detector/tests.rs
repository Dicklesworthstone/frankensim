use super::*;

fn relative(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1.0e-300),
        "{actual} != {expected}"
    );
}

#[test]
fn g1_infrared_flat_channels_keep_radiometric_and_photon_units() {
    let sensor = QuantumEfficiencySensor::try_new(vec![
        (1200.0, [1.0, 0.25, 0.0]),
        (1800.0, [1.0, 0.25, 0.0]),
    ])
    .unwrap();
    assert_eq!(sensor.wavelength_bounds_nm(), (1200.0, 1800.0));
    let mut power = [0.0; 3];
    let mut electrons = [0.0; 3];
    for u in [0.125, 0.375, 0.625, 0.875] {
        let draw = sensor.sample_wavelength(u).unwrap();
        relative(draw.wavelength_nm(), 1200.0 + 600.0 * u, 2e-15);
        for c in 0..3 {
            power[c] += 3.0 * draw.radiance_weights_nm()[c] / 4.0;
            electrons[c] += 3.0 * draw.electron_weights_nm_per_joule()[c] / 4.0;
        }
    }
    let h = 6.626_070_15e-34;
    let c = 299_792_458.0;
    relative(power[0], 1800.0, 2e-15);
    relative(power[1], 450.0, 2e-15);
    relative(electrons[0], 3.0 * 600.0 * 1500.0e-9 / (h * c), 3e-15);
    relative(electrons[1], electrons[0] * 0.25, 2e-15);
    assert_eq!(power[2], 0.0);
    assert_eq!(electrons[2], 0.0);
}

#[test]
fn g1_linear_efficiency_inverts_its_actual_probability_and_counts_longer_photons() {
    let rising = QuantumEfficiencySensor::try_new(vec![(8000.0, [0.0]), (14000.0, [0.8])]).unwrap();
    let falling =
        QuantumEfficiencySensor::try_new(vec![(8000.0, [0.8]), (14000.0, [0.0])]).unwrap();
    for u in [0.0, 0.01, 0.25, 0.5, 0.81, 0.999] {
        let up = rising.sample_wavelength(u).unwrap();
        let down = falling.sample_wavelength(u).unwrap();
        let x = (up.wavelength_nm() - 8000.0) / 6000.0;
        let y = (down.wavelength_nm() - 8000.0) / 6000.0;
        assert!((x * x - u).abs() < 5e-15);
        assert!((2.0 * y - y * y - u).abs() < 5e-15);
        relative(up.radiance_weights_nm()[0], 2400.0, 2e-15);
        relative(down.radiance_weights_nm()[0], 2400.0, 2e-15);
    }
    // Independent first moment of a linear QE ramp: mean wavelength is
    // a+2*(b-a)/3, not its midpoint or its vacuum frequency counterpart.
    let mut mean = Dd::ZERO;
    let count = 65536;
    for i in 0..count {
        let draw = rising
            .sample_wavelength((i as f64 + 0.5) / count as f64)
            .unwrap();
        mean = mean + Dd::from_f64(draw.electron_weights_nm_per_joule()[0] / count as f64);
    }
    let expected = 2400.0 * 12000.0e-9 / (6.626_070_15e-34 * 299_792_458.0);
    relative(mean.to_f64(), expected, 1e-8);
}

#[test]
fn g1_disjoint_channels_do_not_sample_dark_gaps_or_drop_a_channel() {
    let sensor = QuantumEfficiencySensor::try_new(vec![
        (900.0, [0.0, 0.0]),
        (1000.0, [1.0, 0.0]),
        (1100.0, [0.0, 0.0]),
        (1500.0, [0.0, 0.0]),
        (1600.0, [0.0, 0.5]),
        (1700.0, [0.0, 0.0]),
    ])
    .unwrap();
    let mut integral = [Dd::ZERO; 2];
    let count = 6000;
    for i in 0..count {
        let draw = sensor
            .sample_wavelength((i as f64 + 0.5) / count as f64)
            .unwrap();
        assert!(!(1100.0..=1500.0).contains(&draw.wavelength_nm()));
        for c in 0..2 {
            integral[c] = integral[c] + Dd::from_f64(draw.radiance_weights_nm()[c] / count as f64);
        }
    }
    relative(integral[0].to_f64(), 100.0, 3e-15);
    relative(integral[1].to_f64(), 50.0, 3e-15);
    assert_eq!(sensor.quantum_efficiency(1300.0).unwrap(), [0.0; 2]);
    for wavelength in [550.0, 899.0, 1701.0, f64::NAN] {
        assert_eq!(
            sensor.quantum_efficiency(wavelength),
            Err(SpectralSensorError::InvalidWavelength)
        );
    }
}

#[test]
fn g3_weak_efficiency_changes_signal_without_changing_wavelengths() {
    let a =
        QuantumEfficiencySensor::try_new(vec![(1000.0, [0.2, 0.8]), (2000.0, [0.8, 0.2])]).unwrap();
    let b = QuantumEfficiencySensor::try_new(vec![
        (1000.0, [2e-201, 8e-201]),
        (2000.0, [8e-201, 2e-201]),
    ])
    .unwrap();
    for u in [0.0, 0.125, 0.5, 0.999] {
        let x = a.sample_wavelength(u).unwrap();
        let y = b.sample_wavelength(u).unwrap();
        relative(x.wavelength_nm(), y.wavelength_nm(), 3e-15);
        for c in 0..2 {
            relative(
                y.radiance_weights_nm()[c],
                x.radiance_weights_nm()[c] * 1e-200,
                3e-15,
            );
            relative(
                y.electron_weights_nm_per_joule()[c],
                x.electron_weights_nm_per_joule()[c] * 1e-200,
                3e-15,
            );
        }
    }
}

#[test]
fn g0_invalid_detector_sources_and_variates_refuse() {
    for knots in [
        vec![],
        vec![(1000.0, [0.5])],
        vec![(0.0, [0.5]), (2000.0, [0.5])],
        vec![(2000.0, [0.5]), (1000.0, [0.5])],
        vec![(1000.0, [0.5]), (1000.0, [0.5])],
        vec![(1000.0, [0.5]), (f64::INFINITY, [0.5])],
        vec![(1000.0, [-0.1]), (2000.0, [0.5])],
        vec![(1000.0, [1.1]), (2000.0, [0.5])],
        vec![(1000.0, [f64::NAN]), (2000.0, [0.5])],
        vec![(1000.0, [0.0]), (2000.0, [0.0])],
        vec![(1000.0, [0.5]); MAX_SENSOR_KNOTS + 1],
    ] {
        assert!(QuantumEfficiencySensor::try_new(knots).is_err());
    }
    assert!(QuantumEfficiencySensor::<0>::try_new(vec![(1000.0, []), (2000.0, [])]).is_err());
    assert!(
        QuantumEfficiencySensor::<17>::try_new(vec![(1000.0, [0.5; 17]), (2000.0, [0.5; 17])])
            .is_err()
    );
    let sensor = QuantumEfficiencySensor::try_new(vec![(1000.0, [0.5]), (2000.0, [0.5])]).unwrap();
    for u in [-1.0, 1.0, f64::NAN, f64::INFINITY] {
        assert_eq!(
            sensor.sample_wavelength(u),
            Err(SpectralSensorError::InvalidVariate)
        );
    }
}

#[test]
fn g0_positive_subnormal_channels_refuse_underflow_instead_of_becoming_dark() {
    let tiny = f64::from_bits(1);
    let weak = QuantumEfficiencySensor::try_new(vec![
        (1000.0, [tiny, 1.0, 1.0]),
        (2000.0, [tiny, 1.0, 1.0]),
    ])
    .unwrap();
    // QE/p = QE * 1000 nm is representable, but QE/(sum QE) is not.
    assert!(tiny * 1000.0 > 0.0);
    assert_eq!(
        weak.sample_wavelength(0.5),
        Err(SpectralSensorError::NumericalRange)
    );
    assert_eq!(
        weak.quantum_efficiency(1500.0),
        Err(SpectralSensorError::NumericalRange)
    );
    assert_eq!(weak.quantum_efficiency(1000.0).unwrap()[0], tiny);

    let dark = QuantumEfficiencySensor::try_new(vec![
        (1000.0, [0.0, 1.0, 1.0]),
        (2000.0, [0.0, 1.0, 1.0]),
    ])
    .unwrap();
    let draw = dark.sample_wavelength(0.5).unwrap();
    assert_eq!(draw.radiance_weights_nm()[0], 0.0);
    assert_eq!(draw.electron_weights_nm_per_joule()[0], 0.0);
    assert_eq!(dark.quantum_efficiency(1500.0).unwrap()[0], 0.0);

    // The same refusal applies to a positive one-sided channel at an exactly
    // zero-density endpoint; ordinary ramps retain their finite analytic limit.
    let ramp =
        QuantumEfficiencySensor::try_new(vec![(1000.0, [0.0; 3]), (2000.0, [tiny, 1.0, 1.0])])
            .unwrap();
    assert_eq!(
        ramp.sample_wavelength(0.0),
        Err(SpectralSensorError::NumericalRange)
    );
    assert_eq!(ramp.quantum_efficiency(1000.0).unwrap(), [0.0; 3]);
    assert_eq!(
        ramp.quantum_efficiency(1500.0),
        Err(SpectralSensorError::NumericalRange)
    );

    // A weak supported endpoint must not vanish when normalizing the proposal,
    // even if a different interval gives the sensor an ordinary total area.
    assert!(matches!(
        QuantumEfficiencySensor::try_new(vec![
            (1000.0, [tiny, 0.0, 0.0]),
            (2000.0, [0.0; 3]),
            (3000.0, [0.0, 1.0, 1.0]),
        ]),
        Err(SpectralSensorError::NumericalRange)
    ));
}
