use super::*;

const C: f64 = VACUUM_LIGHT_SPEED_M_S;

fn dd(value: f64) -> Dd {
    Dd::from_f64(value)
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1e-300),
        "actual {actual:.17e}, expected {expected:.17e}"
    );
}

fn source(velocity: [f64; 3]) -> UniformVacuumWorldline {
    UniformVacuumWorldline::try_new(
        Dd::ZERO,
        Dd::ZERO,
        [0.0; 3],
        velocity,
        [dd(-1e30), dd(1e30)],
    )
    .unwrap()
}

fn receiver(time: Dd, position: [f64; 3], velocity: [f64; 3]) -> VacuumReception {
    VacuumReception::try_new(time, position, velocity).unwrap()
}

#[test]
fn g1_stationary_null_ray_retains_reference_proper_clock() {
    let motion = UniformVacuumWorldline::try_new(
        Dd::ZERO,
        dd(10.0),
        [0.0; 3],
        [0.0; 3],
        [dd(0.0), dd(20.0)],
    )
    .unwrap();
    let ray = motion
        .observe(receiver(dd(4.0), [3.0 * C, 0.0, 0.0], [0.0; 3]), 0.0)
        .unwrap();
    assert_eq!(ray.emission_coordinate_time_s(), dd(1.0));
    assert_eq!(ray.emission_proper_time_s(), dd(11.0));
    assert_eq!(ray.flight_time_s(), dd(3.0));
    assert_eq!(ray.range_m(), 3.0 * C);
    assert_eq!(ray.world_propagation_direction(), [1.0, 0.0, 0.0]);
    assert_eq!(ray.emission_position_m(), [Dd::ZERO; 3]);
    assert_eq!(ray.doppler().frequency_ratio(), 1.0);
    assert_eq!(ray.coordinate_time_rate(), 1.0);
    assert_eq!(ray.light_cone_residual_m(), 0.0);
}

#[test]
fn g1_longitudinal_approach_and_recession_select_different_emission_events() {
    for (beta, delay, proper, ratio) in [(0.6, 2.5, -2.0, 2.0), (-0.6, 0.625, -0.5, 0.5)] {
        let ray = source([beta * C, 0.0, 0.0])
            .observe(receiver(Dd::ZERO, [C, 0.0, 0.0], [0.0; 3]), 1e-6)
            .unwrap();
        close(ray.flight_time_s().to_f64(), delay, 2e-15);
        close(ray.emission_proper_time_s().to_f64(), proper, 2e-15);
        close(ray.range_m(), C * delay, 2e-15);
        close(ray.doppler().frequency_ratio(), ratio, 2e-15);
        close(ray.coordinate_time_rate(), 1.0 / (1.0 - beta), 3e-15);
        close(
            ray.emission_position_m()[0].to_f64(),
            -beta * C * delay,
            2e-15,
        );
    }
}

#[test]
fn g1_transverse_present_position_requires_an_oblique_retarded_ray() {
    let ray = source([0.6 * C, 0.0, 0.0])
        .observe(receiver(Dd::ZERO, [0.0, C, 0.0], [0.0; 3]), 1e-6)
        .unwrap();
    close(ray.flight_time_s().to_f64(), 1.25, 2e-15);
    close(ray.emission_proper_time_s().to_f64(), -1.0, 2e-15);
    let direction = ray.world_propagation_direction();
    close(direction[0], 0.6, 2e-15);
    close(direction[1], 0.8, 2e-15);
    // Using the instantaneous transverse position would give D=0.8, which
    // is not the Doppler shift of the photon that actually reaches this event.
    close(ray.doppler().frequency_ratio(), 1.25, 3e-15);
    assert!(ray.doppler().source_direction()[0].abs() < 2e-15);
    close(ray.doppler().source_direction()[1], 1.0, 2e-15);
}

#[test]
fn g1_receiver_motion_derivative_matches_both_coordinate_and_proper_clocks() {
    let velocity_source = [0.3 * C, -0.2 * C, 0.1 * C];
    let velocity_receiver = [-0.1 * C, 0.15 * C, 0.05 * C];
    let center = [2.0 * C, C, -0.3 * C];
    let motion = UniformVacuumWorldline::try_new(
        dd(3.0),
        dd(2.0),
        [0.4 * C, -0.2 * C, 0.1 * C],
        velocity_source,
        [dd(-20.0), dd(20.0)],
    )
    .unwrap();
    let at = |offset: f64| {
        let position = core::array::from_fn(|i| center[i] + velocity_receiver[i] * offset);
        motion
            .observe(
                receiver(dd(4.0) + dd(offset), position, velocity_receiver),
                1e-6,
            )
            .unwrap()
    };
    let h = 1e-5;
    let before = at(-h);
    let after = at(h);
    let middle = at(0.0);
    let coordinate = (after.emission_coordinate_time_s() - before.emission_coordinate_time_s())
        .to_f64()
        / (2.0 * h);
    close(coordinate, middle.coordinate_time_rate(), 4e-10);
    let receiver_clock_rate = (1.0_f64 - 0.1 * 0.1 - 0.15 * 0.15 - 0.05 * 0.05).sqrt();
    let proper = (after.emission_proper_time_s() - before.emission_proper_time_s()).to_f64()
        / (2.0 * h * receiver_clock_rate);
    close(proper, middle.proper_time_rate(), 4e-10);
    assert_eq!(
        middle.proper_time_rate(),
        middle.doppler().frequency_ratio()
    );
}

#[test]
fn g3_dd_epochs_keep_resolvable_flights_and_refuse_lost_delays() {
    let epoch = dd(1e20);
    let motion =
        UniformVacuumWorldline::try_new(epoch, Dd::ZERO, [0.0; 3], [0.0; 3], [dd(-1e30), dd(1e30)])
            .unwrap();
    let reception = epoch + dd(0.25);
    let ray = motion
        .observe(receiver(reception, [C * 1e-6, 0.0, 0.0], [0.0; 3]), 1e-6)
        .unwrap();
    assert_eq!(ray.emission_coordinate_time_s().hi, epoch.hi);
    assert!(ray.emission_coordinate_time_s().lo < reception.lo);
    close(ray.flight_time_s().to_f64(), 1e-6, 1e-10);
    close(ray.emission_proper_time_s().to_f64(), 0.25 - 1e-6, 2e-15);

    // Two terms are not arbitrary precision: the low term itself cannot hold
    // a nanosecond difference at this epoch/offset combination.
    let epoch = dd(2.0_f64.powi(100));
    let motion =
        UniformVacuumWorldline::try_new(epoch, Dd::ZERO, [0.0; 3], [0.0; 3], [dd(-1e50), dd(1e50)])
            .unwrap();
    let unresolved = receiver(epoch + dd(2.0_f64.powi(45)), [1.0, 0.0, 0.0], [0.0; 3]);
    assert!(matches!(
        motion.observe(unresolved, 1e-6),
        Err(RetardedVacuumError::NumericalRange | RetardedVacuumError::Residual { .. })
    ));
}

#[test]
fn g3_near_c_geometry_and_proper_rates_use_the_same_physical_si_velocity() {
    let speed = C.next_down();
    for sign in [1.0, -1.0] {
        let v = sign * speed;
        let ray = source([v, 0.0, 0.0])
            .observe(receiver(Dd::ZERO, [C, 0.0, 0.0], [0.0; 3]), 1e-5)
            .unwrap();
        let expected_delay = C / (C - v);
        let expected_frequency = ((C + v) / (C - v)).sqrt();
        close(ray.flight_time_s().to_f64(), expected_delay, 3e-15);
        close(ray.proper_time_rate(), expected_frequency, 4e-15);
        let direct = VacuumDoppler::try_new([1.0, 0.0, 0.0], [v, 0.0, 0.0], [0.0; 3]).unwrap();
        assert_eq!(ray.doppler(), direct);
    }
}

#[test]
fn g3_unresolved_material_clock_is_refused_even_when_null_geometry_is_resolved() {
    let make = |proper_epoch| {
        UniformVacuumWorldline::try_new(
            Dd::ZERO,
            proper_epoch,
            [0.0; 3],
            [0.0; 3],
            [dd(-1e50), dd(1e50)],
        )
        .unwrap()
    };
    let observation = receiver(Dd::ZERO, [1.0, 0.0, 0.0], [0.0; 3]);
    let resolved = make(dd(1e20)).observe(observation, 1e-6).unwrap();
    close(resolved.flight_time_s().to_f64(), 1.0 / C, 2e-15);
    assert!(resolved.emission_proper_time_s().lo < 0.0);
    assert!(resolved.proper_clock_residual_s() <= 1e-6 / C);
    let unresolved = make(dd(2.0_f64.powi(100)) + dd(2.0_f64.powi(45)));
    assert!(matches!(
        unresolved.observe(observation, 1e-6),
        Err(RetardedVacuumError::ClockResolution { .. })
    ));
}

#[test]
fn g3_bounded_histories_coincidence_and_invalid_events_refuse_explicitly() {
    let bounded = UniformVacuumWorldline::try_new(
        Dd::ZERO,
        Dd::ZERO,
        [0.0; 3],
        [0.0; 3],
        [dd(-1.0), Dd::ZERO],
    )
    .unwrap();
    assert!(
        bounded
            .observe(receiver(Dd::ZERO, [C, 0.0, 0.0], [0.0; 3]), 0.0)
            .is_ok()
    );
    assert_eq!(
        bounded.observe(receiver(Dd::ZERO, [2.0 * C, 0.0, 0.0], [0.0; 3]), 1e-6),
        Err(RetardedVacuumError::OutsideHistory)
    );
    assert_eq!(
        bounded.observe(receiver(Dd::ZERO, [0.0; 3], [0.0; 3]), 1e-6),
        Err(RetardedVacuumError::Coincident)
    );
    assert!(matches!(
        VacuumReception::try_new(dd(f64::NAN), [0.0; 3], [0.0; 3]),
        Err(RetardedVacuumError::InvalidClock)
    ));
    assert!(matches!(
        VacuumReception::try_new(Dd::ZERO, [f64::INFINITY, 0.0, 0.0], [0.0; 3]),
        Err(RetardedVacuumError::InvalidPosition)
    ));
    assert!(matches!(
        VacuumReception::try_new(Dd::ZERO, [0.0; 3], [C, 0.0, 0.0]),
        Err(RetardedVacuumError::Doppler(DopplerError::InvalidVelocity))
    ));
    assert!(matches!(
        UniformVacuumWorldline::try_new(
            Dd::ZERO,
            Dd::ZERO,
            [0.0; 3],
            [0.0; 3],
            [dd(1.0), Dd::ZERO]
        ),
        Err(RetardedVacuumError::InvalidHistory)
    ));
    for tolerance in [-1.0, f64::NAN, f64::INFINITY] {
        assert_eq!(
            bounded.observe(receiver(Dd::ZERO, [C, 0.0, 0.0], [0.0; 3]), tolerance),
            Err(RetardedVacuumError::InvalidTolerance)
        );
    }
}

#[test]
fn g1_uniform_detector_clock_resolves_longitudinal_and_transverse_doppler() {
    let detector = |velocity| {
        UniformVacuumWorldline::try_new(
            Dd::ZERO,
            dd(10.0),
            [3.0 * C, 0.0, 0.0],
            velocity,
            [dd(8.0), dd(12.0)],
        )
        .unwrap()
    };
    for (beta, ratio) in [(-0.6, 2.0), (0.6, 0.5)] {
        let proper = dd(10.0) + dd(0.8);
        let event = detector([beta * C, 0.0, 0.0])
            .reception_at_proper_time(proper, 1e-6)
            .unwrap();
        close(event.coordinate_time_s().to_f64(), 1.0, 3e-15);
        close(event.position_m()[0].to_f64(), (3.0 + beta) * C, 3e-15);
        let ray = source([0.0; 3]).observe(event, 1e-6).unwrap();
        close(ray.emission_proper_time_s().to_f64(), -2.0 - beta, 3e-15);
        close(ray.proper_time_rate(), ratio, 3e-15);
        assert_eq!(ray.reception_proper_time_s(), Some(proper));
        assert_eq!(
            ray.reception_clock_residual_s(),
            event.proper_clock_residual_s()
        );
    }
    let transverse = detector([0.0, 0.6 * C, 0.0])
        .reception_at_proper_time(dd(10.0), 1e-6)
        .unwrap();
    let ray = source([0.0; 3]).observe(transverse, 1e-6).unwrap();
    // A transverse moving receiver sees gamma, whereas the transverse moving
    // source law uses inverse gamma. The two endpoint clocks must not be swapped.
    close(ray.proper_time_rate(), 1.25, 3e-15);
    close(ray.doppler().receiver_direction()[1], -0.6, 3e-15);

    let common = [0.6 * C, 0.0, 0.0];
    let event = detector(common)
        .reception_at_proper_time(dd(10.8), 1e-6)
        .unwrap();
    assert_eq!(
        source(common)
            .observe(event, 1e-6)
            .unwrap()
            .proper_time_rate(),
        1.0
    );
}

#[test]
fn g1_detector_proper_clock_derivative_is_the_lorentz_frequency_ratio() {
    let emitter = source([0.3 * C, -0.2 * C, 0.1 * C]);
    let detector = UniformVacuumWorldline::try_new(
        Dd::ZERO,
        dd(12.0),
        [3.0 * C, C, -0.2 * C],
        [-0.1 * C, 0.25 * C, 0.1 * C],
        [dd(11.0), dd(14.0)],
    )
    .unwrap();
    let at = |tau| {
        emitter
            .observe(detector.reception_at_proper_time(tau, 1e-6).unwrap(), 1e-6)
            .unwrap()
    };
    let center = dd(12.6);
    let h = dd(1e-5);
    let derivative = (at(center + h).emission_proper_time_s()
        - at(center - h).emission_proper_time_s())
        / (dd(2.0) * h);
    close(derivative.to_f64(), at(center).proper_time_rate(), 1e-9);
}

#[test]
fn g3_detector_geometry_keeps_dd_displacements_and_rechecks_clock_loss() {
    let make = |coordinate_epoch, velocity| {
        UniformVacuumWorldline::try_new(
            coordinate_epoch,
            Dd::ZERO,
            [1e20, 0.0, 0.0],
            velocity,
            [Dd::ZERO, dd(1.0)],
        )
        .unwrap()
    };
    let detector = make(Dd::ZERO, [1.0, 0.0, 0.0]);
    let reception = detector.reception_at_proper_time(dd(0.25), 1e-6).unwrap();
    close((reception.position_m()[0] - dd(1e20)).to_f64(), 0.25, 3e-15);
    let ray = make(Dd::ZERO, [0.0; 3]).observe(reception, 1e-6).unwrap();
    close(ray.range_m(), 0.25, 3e-15);

    let unresolved_epoch = dd(2.0_f64.powi(100)) + dd(2.0_f64.powi(45));
    let unresolved = make(unresolved_epoch, [0.0; 3]);
    assert!(matches!(
        unresolved.reception_at_proper_time(dd(0.001), 1e-6),
        Err(RetardedVacuumError::ClockResolution { .. })
    ));
    let loose = unresolved.reception_at_proper_time(dd(0.001), 1e9).unwrap();
    assert!(loose.proper_clock_residual_s() > 0.0);
    assert!(matches!(
        source([0.0; 3]).observe(loose, 1e-6),
        Err(RetardedVacuumError::ClockResolution { .. })
    ));
}

#[test]
fn g3_uniform_receiver_support_and_clock_admission_are_explicit() {
    let detector = UniformVacuumWorldline::try_new(
        dd(5.0),
        dd(10.0),
        [C, 0.0, 0.0],
        [0.0; 3],
        [dd(9.0), dd(11.0)],
    )
    .unwrap();
    for (proper, coordinate) in [(9.0, 4.0), (11.0, 6.0)] {
        assert_eq!(
            detector
                .reception_at_proper_time(dd(proper), 0.0)
                .unwrap()
                .coordinate_time_s(),
            dd(coordinate)
        );
    }
    for proper in [8.99, 11.01] {
        assert_eq!(
            detector.reception_at_proper_time(dd(proper), 1e-6),
            Err(RetardedVacuumError::OutsideHistory)
        );
    }
    for proper in [f64::NAN, f64::INFINITY] {
        assert_eq!(
            detector.reception_at_proper_time(dd(proper), 1e-6),
            Err(RetardedVacuumError::InvalidClock)
        );
    }
    for tolerance in [-1.0, f64::NAN, f64::INFINITY] {
        assert_eq!(
            detector.reception_at_proper_time(dd(10.0), tolerance),
            Err(RetardedVacuumError::InvalidTolerance)
        );
    }
    let authored = receiver(dd(5.0), [C, 0.0, 0.0], [0.0; 3]);
    assert_eq!(authored.proper_time_s(), None);
    assert_eq!(
        source([0.0; 3])
            .observe(authored, 0.0)
            .unwrap()
            .reception_proper_time_s(),
        None
    );
}
