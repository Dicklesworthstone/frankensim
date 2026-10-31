//! Independent retarded Green-function and density laws for prescribed SI force.
use fs_aeroac::curle3d::CompactForceObserver;
use fs_material::fluid_acoustics::FluidAcousticState;
use std::f64::consts::PI;

fn medium(rho: f64, c: f64) -> FluidAcousticState {
    FluidAcousticState::try_new(rho, c).unwrap()
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 2e-13 * expected.abs().max(1.0),
        "actual={actual:e}, expected={expected:e}"
    );
}

#[test]
fn force_dipole_has_separate_radial_laws_and_signed_angular_response() {
    let source = [0.0; 3];
    let observe =
        |position| CompactForceObserver::new(medium(1.2, 300.0), position, source).unwrap();
    let force = [8.0, 0.0, 0.0];
    let rate = [24.0, 0.0, 0.0];
    let at_two = observe([2.0, 0.0, 0.0]);
    let at_four = observe([4.0, 0.0, 0.0]);
    let near = at_two.pressure_at_retarded(force, [0.0; 3]).unwrap();
    let far = at_two.far_field_at_retarded(rate).unwrap();
    close(near, 8.0 / (4.0 * PI * 4.0));
    close(far, 24.0 / (4.0 * PI * 300.0 * 2.0));
    close(
        at_four.pressure_at_retarded(force, [0.0; 3]).unwrap(),
        near / 4.0,
    );
    close(at_four.far_field_at_retarded(rate).unwrap(), far / 2.0);
    close(
        at_two.pressure_at_retarded(force, rate).unwrap(),
        near + far,
    );
    close(
        observe([-2.0, 0.0, 0.0])
            .pressure_at_retarded(force, rate)
            .unwrap(),
        -(near + far),
    );
    close(
        observe([0.0, 2.0, 0.0])
            .pressure_at_retarded(force, rate)
            .unwrap(),
        0.0,
    );
    close(
        observe([1.2, 1.6, 0.0])
            .pressure_at_retarded(force, rate)
            .unwrap(),
        0.6 * (near + far),
    );
    let moved =
        CompactForceObserver::new(medium(1.2, 300.0), [8.0, 5.0, -2.0], [8.0, 3.0, -2.0]).unwrap();
    close(
        moved
            .pressure_at_retarded([0.0, 8.0, 0.0], [0.0, 24.0, 0.0])
            .unwrap(),
        near + far,
    );
}

#[test]
fn full_field_matches_spatial_derivative_of_retarded_green_function() {
    let source = [0.4, -0.2, 0.3];
    let position = [3.0, 4.0, 1.0];
    let observer = CompactForceObserver::new(medium(1.2, 320.0), position, source).unwrap();
    let time = 0.024;
    let omega = 123.0;
    let history = |t: f64| {
        let (s, c) = (omega * t).sin_cos();
        (
            [2.0 * s, -3.0 * c, 0.5 * s],
            [2.0 * omega * c, 3.0 * omega * s, 0.5 * omega * c],
        )
    };
    let (force, rate) = history(time - observer.delay_s());
    let pressure = observer.pressure_at_retarded(force, rate).unwrap();
    // Independently differentiate -F_i(t-r/c)/(4 pi r) in observer space.
    // This couples the retarded phase, induction sign and all three directions.
    let potential = |x: [f64; 3], component: usize| {
        let r =
            ((x[0] - source[0]).powi(2) + (x[1] - source[1]).powi(2) + (x[2] - source[2]).powi(2))
                .sqrt();
        history(time - r / 320.0).0[component] / (4.0 * PI * r)
    };
    let h = 2e-5;
    let mut derivative = 0.0;
    for axis in 0..3 {
        let mut plus = position;
        let mut minus = position;
        plus[axis] += h;
        minus[axis] -= h;
        derivative -= (potential(plus, axis) - potential(minus, axis)) / (2.0 * h);
    }
    assert!((pressure - derivative).abs() < 3e-9 * pressure.abs());
    let (unretarded_force, unretarded_rate) = history(time);
    assert!(
        (pressure
            - observer
                .pressure_at_retarded(unretarded_force, unretarded_rate)
                .unwrap())
        .abs()
            > 0.1 * pressure.abs(),
        "omitting travel time must be observable"
    );
}

#[test]
fn prescribed_force_is_density_independent_but_propagates_at_material_sound_speed() {
    let prepare =
        |rho, c| CompactForceObserver::new(medium(rho, c), [2.0, 0.0, 0.0], [0.0; 3]).unwrap();
    let force = [3.0, 4.0, 5.0];
    let rate = [400.0, 500.0, 600.0];
    let baseline = prepare(1.0, 300.0);
    let expected = baseline.pressure_at_retarded(force, rate).unwrap();
    for rho in [0.01, 0.1, 10.0, 1000.0] {
        let observer = prepare(rho, 300.0);
        assert_eq!(observer.acoustic_state().density_kg_m3(), rho);
        assert_eq!(
            observer
                .pressure_at_retarded(force, rate)
                .unwrap()
                .to_bits(),
            expected.to_bits()
        );
        assert_eq!(observer.delay_s(), baseline.delay_s());
    }
    let faster = prepare(1.0, 600.0);
    close(faster.delay_s(), baseline.delay_s() / 2.0);
    close(
        faster.far_field_at_retarded(rate).unwrap(),
        baseline.far_field_at_retarded(rate).unwrap() / 2.0,
    );
    close(
        faster.pressure_at_retarded(force, [0.0; 3]).unwrap(),
        baseline.pressure_at_retarded(force, [0.0; 3]).unwrap(),
    );
}

#[test]
fn singular_geometry_and_nonfinite_sources_refuse() {
    for position in [[0.0; 3], [f64::NAN, 0.0, 0.0], [1e-200, 0.0, 0.0]] {
        assert!(CompactForceObserver::new(medium(1.2, 300.0), position, [0.0; 3]).is_err());
    }
    let observer =
        CompactForceObserver::new(medium(1.2, 300.0), [1.0, 0.0, 0.0], [0.0; 3]).unwrap();
    assert!(
        observer
            .pressure_at_retarded([f64::INFINITY, 0.0, 0.0], [0.0; 3])
            .is_err()
    );
    assert!(
        observer
            .far_field_at_retarded([0.0, f64::NAN, 0.0])
            .is_err()
    );
}
