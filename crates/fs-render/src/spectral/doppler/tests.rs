use super::*;

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1.0),
        "actual={actual:e}, expected={expected:e}"
    );
}

fn velocity(beta: [f64; 3]) -> [f64; 3] {
    beta.map(|x| x * VACUUM_LIGHT_SPEED_M_S)
}

#[test]
fn g1_longitudinal_source_and_receiver_have_distinct_signs() {
    let ray = [1.0, 0.0, 0.0];
    for (source, receiver, expected) in [
        ([0.6, 0.0, 0.0], [0.0; 3], 2.0),
        ([-0.6, 0.0, 0.0], [0.0; 3], 0.5),
        ([0.0; 3], [0.6, 0.0, 0.0], 0.5),
        ([0.0; 3], [-0.6, 0.0, 0.0], 2.0),
    ] {
        let view = VacuumDoppler::try_new(ray, velocity(source), velocity(receiver)).unwrap();
        close(view.frequency_ratio(), expected, 2.0e-15);
        close(
            view.source_wavelength_nm(600.0).unwrap(),
            600.0 * expected,
            2.0e-15,
        );
        close(
            view.receiver_wavelength_nm(600.0).unwrap(),
            600.0 / expected,
            2.0e-15,
        );
    }
}

#[test]
fn g1_transverse_time_dilation_and_aberration() {
    let transverse = velocity([0.6, 0.0, 0.0]);
    let source = VacuumDoppler::try_new([0.0, 1.0, 0.0], transverse, [0.0; 3]).unwrap();
    close(source.frequency_ratio(), 0.8, 2.0e-15);
    let receiver = VacuumDoppler::try_new([0.0, 1.0, 0.0], [0.0; 3], transverse).unwrap();
    close(receiver.frequency_ratio(), 1.25, 2.0e-15);
    for (actual, expected) in receiver
        .receiver_direction()
        .into_iter()
        .zip([-0.6, 0.8, 0.0])
    {
        close(actual, expected, 2.0e-15);
    }
}

#[test]
fn g1_wavelength_radiance_has_fifth_power_jacobian() {
    let view =
        VacuumDoppler::try_new([1.0, 0.0, 0.0], velocity([0.6, 0.0, 0.0]), [0.0; 3]).unwrap();
    // A constant source spectral-radiance band [800,1600] becomes [400,800].
    // Per-nm radiance gains D^5; the corresponding band integral gains D^4.
    let radiance = view.receiver_radiance_w_m2_sr_nm(3.0).unwrap();
    close(radiance, 96.0, 5.0e-15);
    close(radiance * 400.0 / (3.0 * 800.0), 16.0, 5.0e-15);
    close(view.receiver_solid_angle_ratio(), 0.25, 2.0e-15);
}

#[test]
fn g3_comoving_endpoints_preserve_frequency_and_radiance_exactly() {
    let v = velocity([0.4, -0.3, 0.2]);
    let view = VacuumDoppler::try_new([2.0, 3.0, -1.0], v, v).unwrap();
    assert_eq!(view.frequency_ratio(), 1.0);
    assert_eq!(view.source_direction(), view.receiver_direction());
    assert_eq!(
        view.receiver_radiance_w_m2_sr_nm(1.234_567)
            .unwrap()
            .to_bits(),
        1.234_567_f64.to_bits()
    );
    assert_eq!(view.receiver_radiance_w_m2_sr_nm(0.0).unwrap(), 0.0);
}

#[test]
fn g3_boost_and_inverse_restore_null_ray_and_frequency() {
    let ray = [0.2, -0.4, 0.7];
    let v = velocity([0.4, 0.2, -0.1]);
    let first = VacuumDoppler::try_new(ray, [0.0; 3], v).unwrap();
    let inverse =
        VacuumDoppler::try_new(first.receiver_direction(), [0.0; 3], v.map(|x| -x)).unwrap();
    close(
        first.frequency_ratio() * inverse.frequency_ratio(),
        1.0,
        3.0e-15,
    );
    for (actual, expected) in inverse
        .receiver_direction()
        .into_iter()
        .zip(first.source_direction())
    {
        close(actual, expected, 3.0e-15);
    }
    let swapped = VacuumDoppler::try_new(ray, v, [0.0; 3]).unwrap();
    close(
        first.frequency_ratio() * swapped.frequency_ratio(),
        1.0,
        2.0e-15,
    );
}

#[test]
fn g3_near_light_speed_and_large_direction_scale_remain_finite() {
    let v = velocity([1.0 - f64::EPSILON, 0.0, 0.0]);
    let forward = VacuumDoppler::try_new([1.0e300, 0.0, 0.0], v, [0.0; 3]).unwrap();
    assert!(forward.frequency_ratio().is_finite() && forward.frequency_ratio() > 1.0e7);
    assert_eq!(forward.source_direction(), [1.0, 0.0, 0.0]);
    let normal = VacuumDoppler::try_new([1.0, 0.0, 0.0], v, [0.0; 3]).unwrap();
    assert_eq!(forward, normal);
}

#[test]
fn g3_invalid_inputs_and_nonrepresentable_outputs_are_explicit() {
    assert_eq!(
        VacuumDoppler::try_new([0.0; 3], [0.0; 3], [0.0; 3]),
        Err(DopplerError::InvalidDirection)
    );
    assert_eq!(
        VacuumDoppler::try_new([f64::NAN, 1.0, 0.0], [0.0; 3], [0.0; 3]),
        Err(DopplerError::InvalidDirection)
    );
    for v in [
        [VACUUM_LIGHT_SPEED_M_S, 0.0, 0.0],
        [f64::INFINITY, 0.0, 0.0],
        velocity([0.8, 0.8, 0.0]),
    ] {
        assert_eq!(
            VacuumDoppler::try_new([1.0, 0.0, 0.0], v, [0.0; 3]),
            Err(DopplerError::InvalidVelocity)
        );
    }
    let view =
        VacuumDoppler::try_new([1.0, 0.0, 0.0], velocity([0.6, 0.0, 0.0]), [0.0; 3]).unwrap();
    assert_eq!(
        view.source_wavelength_nm(0.0),
        Err(DopplerError::InvalidWavelength)
    );
    assert_eq!(
        view.source_wavelength_nm(f64::MAX),
        Err(DopplerError::NumericalRange)
    );
    assert_eq!(
        view.receiver_radiance_w_m2_sr_nm(-1.0),
        Err(DopplerError::InvalidRadiance)
    );
    assert_eq!(
        view.receiver_radiance_w_m2_sr_nm(f64::MAX),
        Err(DopplerError::NumericalRange)
    );
}
