use super::*;
use crate::spectral::doppler::DopplerError;

const C: f64 = VACUUM_LIGHT_SPEED_M_S;

fn dd(value: f64) -> Dd {
    Dd::from_f64(value)
}

fn knot(time: f64, x: f64) -> VacuumTrajectoryKnot {
    VacuumTrajectoryKnot {
        coordinate_time_s: dd(time),
        position_m: [x * C, 0.0, 0.0],
    }
}

fn reception(time: f64, x: f64) -> VacuumReception {
    VacuumReception::try_new(dd(time), [x * C, 0.0, 0.0], [0.0; 3]).unwrap()
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1.0),
        "actual {actual:.17e}, expected {expected:.17e}"
    );
}

#[test]
fn g1_retarded_segment_retains_all_previous_proper_clock_rates() {
    let motion = PiecewiseVacuumWorldline::try_new(
        &[knot(0.0, 0.0), knot(2.0, 1.2), knot(4.0, 1.2)],
        dd(10.0),
        2,
    )
    .unwrap();
    assert_eq!(motion.segment_count(), 2);
    assert_eq!(motion.coordinate_history_s(), [Dd::ZERO, dd(4.0)]);
    close(motion.proper_history_s()[1].to_f64(), 13.6, 2e-15);
    for (received, emitted, proper, ratio) in [(3.4, 1.0, 10.8, 2.0), (4.8, 3.0, 12.6, 1.0)] {
        let ray = motion.observe(reception(received, 3.0), 1e-6).unwrap();
        close(ray.emission_coordinate_time_s().to_f64(), emitted, 3e-15);
        close(ray.emission_proper_time_s().to_f64(), proper, 3e-15);
        close(ray.doppler().frequency_ratio(), ratio, 3e-15);
        assert!(ray.light_cone_residual_m() <= 1e-6);
    }
}

#[test]
fn g1_velocity_reversal_changes_the_received_shift_at_the_retarded_turn() {
    let motion = PiecewiseVacuumWorldline::try_new(
        &[knot(0.0, 0.0), knot(2.0, 1.2), knot(4.0, 0.0)],
        Dd::ZERO,
        2,
    )
    .unwrap();
    let approaching = motion.observe(reception(3.4, 3.0), 1e-6).unwrap();
    let receding = motion.observe(reception(4.6, 3.0), 1e-6).unwrap();
    close(approaching.doppler().frequency_ratio(), 2.0, 3e-15);
    close(receding.doppler().frequency_ratio(), 0.5, 3e-15);
    close(receding.emission_coordinate_time_s().to_f64(), 2.5, 3e-15);
    close(receding.emission_position_m()[0].to_f64(), 0.9 * C, 3e-15);
    close(receding.emission_proper_time_s().to_f64(), 2.0, 3e-15);
}

#[test]
fn g3_knot_right_limit_and_closed_history_do_not_extrapolate() {
    let motion = PiecewiseVacuumWorldline::try_new(
        &[knot(0.0, 0.0), knot(2.0, 1.0), knot(4.0, 1.0)],
        Dd::ZERO,
        2,
    )
    .unwrap();
    let first = motion.observe(reception(3.0, 3.0), 1e-6).unwrap();
    let turn = motion.observe(reception(4.0, 3.0), 1e-6).unwrap();
    let last = motion.observe(reception(6.0, 3.0), 1e-6).unwrap();
    assert_eq!(first.emission_coordinate_time_s(), Dd::ZERO);
    assert_eq!(turn.emission_coordinate_time_s(), dd(2.0));
    assert_eq!(last.emission_coordinate_time_s(), dd(4.0));
    close(first.proper_time_rate(), 3.0_f64.sqrt(), 3e-15);
    assert_eq!(turn.proper_time_rate(), 1.0);
    assert_eq!(last.proper_time_rate(), 1.0);
    close(
        turn.emission_proper_time_s().to_f64(),
        3.0_f64.sqrt(),
        3e-15,
    );
    for time in [2.99, 6.01] {
        assert_eq!(
            motion.observe(reception(time, 3.0), 1e-6),
            Err(RetardedVacuumError::OutsideHistory)
        );
    }
    assert_eq!(
        motion.observe(reception(0.0, 0.0), 1e-6),
        Err(RetardedVacuumError::Coincident)
    );
}

#[test]
fn g3_subdividing_uniform_motion_preserves_geometry_and_proper_clock() {
    let knots: Vec<_> = (0..=64)
        .map(|i| knot(f64::from(i) / 16.0, f64::from(i) / 32.0))
        .collect();
    let motion = PiecewiseVacuumWorldline::try_new(&knots, dd(7.0), 64).unwrap();
    let uniform = UniformVacuumWorldline::try_new(
        Dd::ZERO,
        dd(7.0),
        [0.0; 3],
        [0.5 * C, 0.0, 0.0],
        motion.proper_history_s(),
    )
    .unwrap();
    for time in [3.1, 3.9, 4.4, 4.99] {
        let receiver =
            VacuumReception::try_new(dd(time), [3.0 * C, 0.1 * C, 0.0], [-0.1 * C, 0.2 * C, 0.0])
                .unwrap();
        let actual = motion.observe(receiver, 1e-6).unwrap();
        let expected = uniform.observe(receiver, 1e-6).unwrap();
        close(
            actual.emission_coordinate_time_s().to_f64(),
            expected.emission_coordinate_time_s().to_f64(),
            3e-15,
        );
        close(
            actual.emission_proper_time_s().to_f64(),
            expected.emission_proper_time_s().to_f64(),
            3e-15,
        );
        close(
            actual.proper_time_rate(),
            expected.proper_time_rate(),
            3e-15,
        );
        for (a, b) in actual
            .world_propagation_direction()
            .into_iter()
            .zip(expected.world_propagation_direction())
        {
            close(a, b, 3e-15);
        }
    }
}

#[test]
fn g3_dd_epochs_retain_offsets_and_refuse_unresolved_material_clock() {
    let epoch = dd(1e20);
    let knots: Vec<_> = (0..=2)
        .map(|i| VacuumTrajectoryKnot {
            coordinate_time_s: epoch + dd(f64::from(i)),
            position_m: [0.0; 3],
        })
        .collect();
    let receiver = VacuumReception::try_new(epoch + dd(0.5), [1.0, 0.0, 0.0], [0.0; 3]).unwrap();
    let resolved = PiecewiseVacuumWorldline::try_new(&knots, epoch, 2)
        .unwrap()
        .observe(receiver, 1e-6)
        .unwrap();
    close(resolved.flight_time_s().to_f64() * C, 1.0, 1e-6);
    close(
        (resolved.emission_proper_time_s() - epoch).to_f64(),
        0.5 - 1.0 / C,
        1e-14,
    );
    let unresolved = dd(2.0_f64.powi(100)) + dd(2.0_f64.powi(45));
    let motion = PiecewiseVacuumWorldline::try_new(&knots, unresolved, 2).unwrap();
    assert!(matches!(
        motion.observe(receiver, 1e-6),
        Err(RetardedVacuumError::ClockResolution { .. })
    ));
}

#[test]
fn g3_invalid_and_superluminal_segments_refuse_before_observation() {
    assert!(matches!(
        PiecewiseVacuumWorldline::try_new(&[knot(0.0, 0.0)], Dd::ZERO, 1),
        Err(RetardedVacuumError::InvalidTrajectory)
    ));
    assert!(matches!(
        PiecewiseVacuumWorldline::try_new(&[knot(0.0, 0.0), knot(1.0, 0.0)], Dd::ZERO, 0),
        Err(RetardedVacuumError::TrajectoryBudget)
    ));
    for end in [0.0, -1.0] {
        assert!(matches!(
            PiecewiseVacuumWorldline::try_new(&[knot(0.0, 0.0), knot(end, 0.0)], Dd::ZERO, 1),
            Err(RetardedVacuumError::InvalidTrajectory)
        ));
    }
    for displacement in [1.0, 1.1] {
        assert!(matches!(
            PiecewiseVacuumWorldline::try_new(
                &[knot(0.0, 0.0), knot(1.0, displacement)],
                Dd::ZERO,
                1
            ),
            Err(RetardedVacuumError::Doppler(DopplerError::InvalidVelocity))
        ));
    }
    let valid =
        PiecewiseVacuumWorldline::try_new(&[knot(0.0, 0.0), knot(1.0, 0.0)], Dd::ZERO, 1).unwrap();
    for tolerance in [-1.0, f64::NAN, f64::INFINITY] {
        assert_eq!(
            valid.observe(reception(1.0, 0.5), tolerance),
            Err(RetardedVacuumError::InvalidTolerance)
        );
    }
}

#[test]
fn g1_detector_clock_inversion_preserves_preceding_segment_time_dilation() {
    let detector = PiecewiseVacuumWorldline::try_new(
        &[knot(0.0, 3.0), knot(2.0, 1.8), knot(4.0, 1.8)],
        dd(10.0),
        2,
    )
    .unwrap();
    let source = UniformVacuumWorldline::try_new(
        Dd::ZERO,
        Dd::ZERO,
        [0.0; 3],
        [0.0; 3],
        [dd(-10.0), dd(10.0)],
    )
    .unwrap();
    for (proper, coordinate, emitted, ratio) in [(10.8, 1.0, -1.4, 2.0), (12.6, 3.0, 1.2, 1.0)] {
        let event = detector.reception_at_proper_time(dd(proper), 1e-6).unwrap();
        close(event.coordinate_time_s().to_f64(), coordinate, 3e-15);
        assert_eq!(event.proper_time_s(), Some(dd(proper)));
        let ray = source.observe(event, 1e-6).unwrap();
        close(ray.emission_proper_time_s().to_f64(), emitted, 3e-15);
        close(ray.proper_time_rate(), ratio, 3e-15);
        let h = dd(1e-5);
        let emitted_at = |tau| {
            source
                .observe(detector.reception_at_proper_time(tau, 1e-6).unwrap(), 1e-6)
                .unwrap()
                .emission_proper_time_s()
        };
        let derivative = (emitted_at(dd(proper) + h) - emitted_at(dd(proper) - h)) / (dd(2.0) * h);
        close(derivative.to_f64(), ratio, 1e-9);
    }
}

#[test]
fn g3_receiver_knots_use_outgoing_velocity_and_finite_clock_support() {
    let detector = PiecewiseVacuumWorldline::try_new(
        &[knot(0.0, 3.0), knot(2.0, 2.0), knot(4.0, 3.0)],
        dd(7.0),
        2,
    )
    .unwrap();
    let source = UniformVacuumWorldline::try_new(
        Dd::ZERO,
        Dd::ZERO,
        [0.0; 3],
        [0.0; 3],
        [dd(-10.0), dd(10.0)],
    )
    .unwrap();
    for (index, time, ratio) in [
        (0, 0.0, 3.0_f64.sqrt()),
        (1, 2.0, 1.0 / 3.0_f64.sqrt()),
        (2, 4.0, 1.0 / 3.0_f64.sqrt()),
    ] {
        let proper = detector.proper_knot_time(index).unwrap();
        let event = detector.reception_at_proper_time(proper, 1e-6).unwrap();
        assert_eq!(event.coordinate_time_s(), dd(time));
        close(
            source.observe(event, 1e-6).unwrap().proper_time_rate(),
            ratio,
            3e-15,
        );
    }
    let support = detector.proper_history_s();
    for proper in [support[0] - dd(0.01), support[1] + dd(0.01)] {
        assert_eq!(
            detector.reception_at_proper_time(proper, 1e-6),
            Err(RetardedVacuumError::OutsideHistory)
        );
    }
    for tolerance in [-1.0, f64::NAN, f64::INFINITY] {
        assert_eq!(
            detector.reception_at_proper_time(support[0], tolerance),
            Err(RetardedVacuumError::InvalidTolerance)
        );
    }
    assert_eq!(
        detector.reception_at_proper_time(dd(f64::NAN), 1e-6),
        Err(RetardedVacuumError::InvalidClock)
    );
    let epoch = dd(2.0_f64.powi(100)) + dd(2.0_f64.powi(45));
    let unresolved = PiecewiseVacuumWorldline::try_new(
        &[knot(0.0, 3.0), knot(2.0, 2.0), knot(4.0, 3.0)],
        epoch,
        2,
    )
    .unwrap();
    assert!(matches!(
        unresolved.reception_at_proper_time(unresolved.proper_history_s()[1], 1e-6),
        Err(RetardedVacuumError::ClockResolution { .. })
    ));
}

#[test]
fn g3_receiver_subdivision_preserves_uniform_events_and_common_motion_shift() {
    let knots: Vec<_> = (0..=64)
        .map(|i| knot(f64::from(i) / 16.0, 3.0 + f64::from(i) / 32.0))
        .collect();
    let detector = PiecewiseVacuumWorldline::try_new(&knots, dd(7.0), 64).unwrap();
    let velocity = [0.5 * C, 0.0, 0.0];
    let uniform = UniformVacuumWorldline::try_new(
        Dd::ZERO,
        dd(7.0),
        [3.0 * C, 0.0, 0.0],
        velocity,
        detector.proper_history_s(),
    )
    .unwrap();
    let source = UniformVacuumWorldline::try_new(
        Dd::ZERO,
        Dd::ZERO,
        [0.0; 3],
        velocity,
        [dd(-20.0), dd(20.0)],
    )
    .unwrap();
    for proper in [7.2, 8.3, 10.2] {
        let actual = detector.reception_at_proper_time(dd(proper), 1e-6).unwrap();
        let expected = uniform.reception_at_proper_time(dd(proper), 1e-6).unwrap();
        close(
            actual.coordinate_time_s().to_f64(),
            expected.coordinate_time_s().to_f64(),
            3e-15,
        );
        close(
            actual.position_m()[0].to_f64(),
            expected.position_m()[0].to_f64(),
            3e-15,
        );
        let ray = source.observe(actual, 1e-6).unwrap();
        let uniform_ray = source.observe(expected, 1e-6).unwrap();
        close(
            ray.emission_proper_time_s().to_f64(),
            uniform_ray.emission_proper_time_s().to_f64(),
            3e-15,
        );
        assert_eq!(ray.proper_time_rate(), 1.0);
        assert_eq!(ray.reception_proper_time_s(), Some(dd(proper)));
    }
}
