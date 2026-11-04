//! Independent retarded Green-function and density laws for prescribed SI force.
use fs_aeroac::curle3d::{CompactForceObserver, CompactSourcePhasor};
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_math::c64::C64;
use std::f64::consts::{PI, TAU};

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

/// G1: transient fields recover the original harmonic owner from the near
/// field to the radiation zone, including the physical receiver-time phase.
#[test]
fn transient_force_field_matches_harmonic_velocity_and_cycle_mean_flux() {
    let frequency = 173.0;
    let omega = TAU * frequency;
    let force = [C64::new(2.0, -0.7), C64::new(-0.4, 0.9), C64::new(1.1, 0.3)];
    let real_at = |phasor: C64, time: f64| {
        let (sin, cos) = (omega * time).sin_cos();
        phasor.re * cos + phasor.im * sin
    };
    for rho in [1.3, 997.0] {
        for distance in [0.003, 0.2, 7.0] {
            let observer = CompactForceObserver::new(
                medium(rho, 317.0),
                [0.6 * distance, 0.8 * distance, 0.0],
                [0.0; 3],
            )
            .unwrap();
            let harmonic = observer
                .prepare_harmonic(frequency)
                .unwrap()
                .evaluate(&[CompactSourcePhasor::force_dipole(force)])
                .unwrap();
            let mut average_flux = [0.0; 3];
            let mut average_magnitude = [0.0; 3];
            for sample in 0..64 {
                let time = (0.137 + sample as f64 / 64.0) / frequency;
                let source_time = time - observer.delay_s();
                let impulse = force.map(|value| {
                    real_at(C64::new(-value.im / omega, value.re / omega), source_time)
                });
                let current_force = force.map(|value| real_at(value, source_time));
                let rate = force.map(|value| {
                    real_at(C64::new(omega * value.im, -omega * value.re), source_time)
                });
                let field = observer
                    .field_at_retarded(impulse, current_force, rate)
                    .unwrap();
                close(field.pressure_pa, real_at(harmonic.pressure_pa, time));
                for axis in 0..3 {
                    let expected = real_at(harmonic.particle_velocity_m_s[axis], time);
                    assert!(
                        (field.particle_velocity_m_s[axis] - expected).abs()
                            < 3e-12 * harmonic.particle_velocity_m_s[axis].abs().max(1e-12)
                    );
                    close(
                        field.instantaneous_intensity_w_m2[axis],
                        field.pressure_pa * field.particle_velocity_m_s[axis],
                    );
                    average_flux[axis] += field.instantaneous_intensity_w_m2[axis] / 64.0;
                    average_magnitude[axis] +=
                        field.instantaneous_intensity_w_m2[axis].abs() / 64.0;
                }
            }
            for axis in 0..3 {
                assert!(
                    (average_flux[axis] - harmonic.active_intensity_w_m2[axis]).abs()
                        < 2e-12 * average_magnitude[axis].max(1e-16),
                    "cycle flux {average_flux:?}, harmonic {:?}",
                    harmonic.active_intensity_w_m2
                );
            }
        }
    }
}

/// G1: a dipole's pressure node is not a zero-velocity point. The impulse and
/// force contribute independently with 1/r³ and 1/r² transverse spreading.
#[test]
fn pressure_null_retains_transverse_impulse_and_force_velocity() {
    let rho = 1.2;
    let c = 300.0;
    for distance in [1.0, 2.0] {
        let observer =
            CompactForceObserver::new(medium(rho, c), [distance, 0.0, 0.0], [0.0; 3]).unwrap();
        let field = observer
            .field_at_retarded([0.0, 2.0, 0.0], [0.0, 3.0, 0.0], [0.0, 5.0, 0.0])
            .unwrap();
        assert_eq!(field.pressure_pa, 0.0);
        assert_eq!(field.particle_velocity_m_s[0], 0.0);
        assert_eq!(field.particle_velocity_m_s[2], 0.0);
        close(
            field.particle_velocity_m_s[1],
            -(2.0 / distance.powi(3) + 3.0 / (c * distance * distance)) / (4.0 * PI * rho),
        );
        assert_eq!(field.instantaneous_intensity_w_m2, [0.0; 3]);
    }
}

/// G1: independently differentiated spatial pressure and time velocity obey
/// the source-free momentum equation, using one coherent polynomial history.
#[test]
fn transient_vector_field_obeys_source_free_momentum() {
    let rho = 1.3;
    let c = 317.0;
    let source = [0.1, -0.2, 0.05];
    let position = [1.2, 0.7, -0.4];
    let history = |t: f64| {
        (
            [t * t * t, -0.5 * t * t, 2.0 * t],
            [3.0 * t * t, -t, 2.0],
            [6.0 * t, -1.0, 0.0],
        )
    };
    let observe = |x, t| {
        let observer = CompactForceObserver::new(medium(rho, c), x, source).unwrap();
        let (impulse, force, rate) = history(t - observer.delay_s());
        observer.field_at_retarded(impulse, force, rate).unwrap()
    };
    let time = 0.31;
    let dt = 1e-5;
    let forward = observe(position, time + dt);
    let backward = observe(position, time - dt);
    for axis in 0..3 {
        let h = 2e-5;
        let mut plus = position;
        let mut minus = position;
        plus[axis] += h;
        minus[axis] -= h;
        let acceleration = (forward.particle_velocity_m_s[axis]
            - backward.particle_velocity_m_s[axis])
            / (2.0 * dt);
        let gradient =
            (observe(plus, time).pressure_pa - observe(minus, time).pressure_pa) / (2.0 * h);
        assert!((rho * acceleration + gradient).abs() < 5e-10 * gradient.abs().max(1.0));
    }
}

/// G3: source amplitude is linear, energy flux quadratic, and material density
/// affects velocity/flux without being spuriously inserted into force pressure.
#[test]
fn transient_amplitude_and_density_preserve_signed_instantaneous_flux() {
    let observe = |rho: f64, amplitude: f64| {
        CompactForceObserver::new(medium(rho, 10.0), [1.0, 0.0, 0.0], [0.0; 3])
            .unwrap()
            .field_at_retarded(
                [-amplitude, 0.2 * amplitude, 0.0],
                [amplitude, 0.3 * amplitude, 0.0],
                [0.0; 3],
            )
            .unwrap()
    };
    let baseline = observe(1.0, 1.0);
    assert!(baseline.pressure_pa > 0.0);
    assert!(baseline.instantaneous_intensity_w_m2[0] < 0.0);
    let amplitude = -3.5;
    let scaled = observe(1.0, amplitude);
    let denser = observe(4.0, 1.0);
    close(scaled.pressure_pa, amplitude * baseline.pressure_pa);
    assert_eq!(denser.pressure_pa, baseline.pressure_pa);
    for axis in 0..3 {
        close(
            scaled.particle_velocity_m_s[axis],
            amplitude * baseline.particle_velocity_m_s[axis],
        );
        close(
            scaled.instantaneous_intensity_w_m2[axis],
            amplitude * amplitude * baseline.instantaneous_intensity_w_m2[axis],
        );
        close(
            denser.particle_velocity_m_s[axis],
            baseline.particle_velocity_m_s[axis] / 4.0,
        );
        close(
            denser.instantaneous_intensity_w_m2[axis],
            baseline.instantaneous_intensity_w_m2[axis] / 4.0,
        );
    }
}

#[test]
fn transient_field_refuses_nonfinite_history_and_unrepresentable_coefficients() {
    let observer =
        CompactForceObserver::new(medium(1.2, 300.0), [1.0, 0.0, 0.0], [0.0; 3]).unwrap();
    for history in 0..3 {
        let mut inputs = [[0.0; 3]; 3];
        inputs[history][1] = f64::NAN;
        assert!(
            observer
                .field_at_retarded(inputs[0], inputs[1], inputs[2])
                .is_err()
        );
    }
    // Pressure remains representable even when its velocity coefficients are
    // not: asking for a field must refuse, without narrowing the old API.
    let tiny_distance =
        CompactForceObserver::new(medium(1.0, 1.0), [1e-110, 0.0, 0.0], [0.0; 3]).unwrap();
    assert!(
        tiny_distance
            .pressure_at_retarded([0.0; 3], [0.0; 3])
            .is_ok()
    );
    assert!(
        tiny_distance
            .field_at_retarded([0.0; 3], [0.0; 3], [0.0; 3])
            .is_err()
    );
}

/// G3: large Cartesian components cancel coherently before pressure is formed;
/// a small surviving source component must not be lost to an ordinary sum.
#[test]
fn radial_force_reduction_preserves_small_amplitude_after_cancellation() {
    let observer = CompactForceObserver::new(medium(1.2, 300.0), [1.0; 3], [0.0; 3]).unwrap();
    let expected = observer
        .pressure_at_retarded([0.0, 1.0, 0.0], [0.0; 3])
        .unwrap();
    let actual = observer
        .pressure_at_retarded([1e16, 1.0, -1e16], [0.0; 3])
        .unwrap();
    close(actual, expected);
    assert!(actual > 0.0);
}
