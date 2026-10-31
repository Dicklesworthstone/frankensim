//! Independent Jeffcott, gyroscopic, full-coordinate and work references.
use fs_math::c64::C64;
use fs_mbd::rotordynamics::{MAX_CAMPBELL_SAMPLES, MAX_ROTOR_DOFS, RotorBearingModel, RotorError};
use std::cell::Cell;

fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * (1.0 + b.abs()),
        "{a:e} != {b:e}"
    );
}
fn complex(a: C64, b: C64, tolerance: f64) {
    near(a.re, b.re, tolerance);
    near(a.im, b.im, tolerance);
}
fn isotropic(mass: f64, damping: f64, stiffness: f64, gyro: f64) -> RotorBearingModel {
    RotorBearingModel::new(
        2,
        &[mass, 0.0, 0.0, mass],
        &[damping, 0.0, 0.0, damping],
        &[stiffness, 0.0, 0.0, stiffness],
        &[0.0, -gyro, gyro, 0.0],
        || false,
    )
    .unwrap()
}

#[test]
fn jeffcott_unbalance_matches_complex_amplitude_phase_and_bearing_work() {
    let (mass, damping, stiffness, eccentricity) = (2.0, 0.6, 8.0, 0.03);
    let model = isotropic(mass, damping, stiffness, 0.0);
    for speed in [0.0, 0.5, 1.75, 2.0, 3.0] {
        // exp(+i omega t): [1,-i] produces [cos(omega t),sin(omega t)].
        let magnitude = mass * eccentricity * speed * speed;
        let force = [C64::from_re(magnitude), C64::new(0.0, -magnitude)];
        let response = model
            .at_speed(speed, || false)
            .unwrap()
            .harmonic_response(speed, &force, 1e-12, || false)
            .unwrap();
        let reference =
            C64::from_re(magnitude) / C64::new(stiffness - mass * speed * speed, damping * speed);
        complex(response.displacement[0], reference, 2e-13);
        complex(
            response.displacement[1],
            C64::new(0.0, -1.0) * reference,
            2e-13,
        );
        let viscous = 0.5
            * damping
            * speed
            * speed
            * response
                .displacement
                .iter()
                .map(|q| q.norm_sq())
                .sum::<f64>();
        near(response.mean_input_power_w, viscous, 2e-13);
        near(response.mean_velocity_force_power_w, viscous, 2e-13);
        near(response.mean_gyroscopic_power_w, 0.0, 2e-13);
        near(response.mean_stiffness_force_power_w, 0.0, 2e-13);
        near(response.mean_inertial_power_w, 0.0, 2e-13);
        near(response.mean_power_residual_w, 0.0, 2e-13);
    }
    let poles = model
        .at_speed(7.0, || false)
        .unwrap()
        .poles(|| false)
        .unwrap();
    let decay = -damping / (2.0 * mass);
    let frequency = (stiffness / mass - decay * decay).sqrt();
    assert_eq!(poles.poles_per_s.len(), 4);
    for pole in poles.poles_per_s {
        near(pole.re, decay, 2e-13);
        near(pole.im.abs(), frequency, 2e-13);
    }
    // The undamped critical speed is sqrt(k/m)=2. At the exact pole the
    // original dynamic stiffness is singular; no invented damping hides it.
    let conservative = isotropic(mass, 0.0, stiffness, 0.0);
    assert!(matches!(
        conservative
            .at_speed(2.0, || false)
            .unwrap()
            .harmonic_response(2.0, &[C64::ONE, C64::ZERO], 1e-12, || false),
        Err(RotorError::Linear(_))
    ));
}

#[test]
fn gyroscopic_splitting_and_signed_whirl_match_circular_mode_equations() {
    let (mass, stiffness, gyro) = (1.0, 4.0, 0.5);
    let model = isotropic(mass, 0.0, stiffness, gyro);
    for speed in [-3.0, 0.0, 3.0] {
        let point = model.at_speed(speed, || false).unwrap();
        let report = point.poles(|| false).unwrap();
        let mut positive: Vec<_> = report
            .poles_per_s
            .iter()
            .filter(|p| p.im > 0.0)
            .map(|p| {
                near(p.re, 0.0, 2e-11);
                p.im
            })
            .collect();
        positive.sort_by(f64::total_cmp);
        let radicand = ((speed * gyro).powi(2) + 4.0 * mass * stiffness).sqrt();
        let split = (speed * gyro).abs();
        near(positive[0], (radicand - split) / (2.0 * mass), 2e-11);
        near(positive[1], (radicand + split) / (2.0 * mass), 2e-11);
        let frequency = 0.9;
        let response = point
            .harmonic_response(frequency, &[C64::ONE, C64::new(0.0, -1.0)], 1e-12, || false)
            .unwrap();
        // G=[[0,-g],[g,0]], so G[1,-i]=+i*g[1,-i].
        let amplitude = 1.0 / (stiffness - mass * frequency * frequency - speed * gyro * frequency);
        complex(response.displacement[0], C64::from_re(amplitude), 2e-12);
        complex(response.displacement[1], C64::new(0.0, -amplitude), 2e-12);
        near(response.mean_input_power_w, 0.0, 2e-12);
        near(response.mean_gyroscopic_power_w, 0.0, 2e-12);
        let reverse = point
            .harmonic_response(-frequency, &[C64::ONE, C64::new(0.0, 1.0)], 1e-12, || false)
            .unwrap();
        for (a, b) in response.displacement.iter().zip(&reverse.displacement) {
            complex(a.conj(), *b, 2e-12);
        }
    }
}

#[test]
fn cross_coupled_bearings_retain_full_response_and_signed_work() {
    let mass = [2.0, 0.25, 0.25, 1.5];
    let damping = [-0.2, 0.3, -0.1, 0.4];
    let stiffness = [5.0, 1.0, -0.5, 7.0];
    let gyro = [0.0, -0.6, 0.6, 0.0];
    let model = RotorBearingModel::new(2, &mass, &damping, &stiffness, &gyro, || false).unwrap();
    assert_eq!(model.damping(), damping);
    assert_eq!(model.stiffness(), stiffness);
    let (spin, frequency) = (2.0, 1.3);
    let force = [C64::new(1.0, 0.2), C64::new(-0.3, 0.7)];
    let response = model
        .at_speed(spin, || false)
        .unwrap()
        .harmonic_response(frequency, &force, 1e-12, || false)
        .unwrap();
    let d: [C64; 4] = core::array::from_fn(|i| {
        C64::new(
            stiffness[i] - frequency * frequency * mass[i],
            frequency * (damping[i] + spin * gyro[i]),
        )
    });
    let determinant = d[0] * d[3] - d[1] * d[2];
    let reference = [
        (d[3] * force[0] - d[1] * force[1]) / determinant,
        (d[0] * force[1] - d[2] * force[0]) / determinant,
    ];
    for i in 0..2 {
        complex(response.displacement[i], reference[i], 2e-12);
    }
    assert!(response.mean_stiffness_force_power_w.abs() > 1e-8);
    near(response.mean_gyroscopic_power_w, 0.0, 2e-12);
    near(response.mean_power_residual_w, 0.0, 2e-12);
    let active = isotropic(1.0, -0.5, 4.0, 0.0);
    let point = active.at_speed(0.0, || false).unwrap();
    assert!(
        point
            .poles(|| false)
            .unwrap()
            .poles_per_s
            .iter()
            .all(|p| p.re > 0.0)
    );
    let response = point
        .harmonic_response(1.0, &[C64::ONE, C64::ZERO], 1e-12, || false)
        .unwrap();
    assert!(response.mean_velocity_force_power_w < 0.0);
    assert!(response.mean_input_power_w < 0.0);
}

#[test]
fn full_mass_coordinate_congruence_preserves_poles_physical_response_and_power() {
    let original = RotorBearingModel::new(
        2,
        &[2.0, 0.0, 0.0, 3.0],
        &[1.0, 0.0, 0.0, 2.0],
        &[5.0, 1.0, 1.0, 9.0],
        &[0.0, -0.75, 0.75, 0.0],
        || false,
    )
    .unwrap();
    // q=T z, T=[[2,1],[0,3]]. Every operator transforms by T^T A T;
    // transformed generalized effort is T^T f. Mass is genuinely coupled.
    let transformed = RotorBearingModel::new(
        2,
        &[8.0, 4.0, 4.0, 29.0],
        &[4.0, 2.0, 2.0, 19.0],
        &[20.0, 16.0, 16.0, 92.0],
        &[0.0, -4.5, 4.5, 0.0],
        || false,
    )
    .unwrap();
    let a = original.at_speed(1.7, || false).unwrap();
    let b = transformed.at_speed(1.7, || false).unwrap();
    let mut original_poles = a.poles(|| false).unwrap().poles_per_s;
    let mut transformed_poles = b.poles(|| false).unwrap().poles_per_s;
    // Conjugate pairs can differ in their last real-part bit; compare by their
    // separated imaginary frequencies instead of assigning identities by QR order.
    original_poles.sort_by(|x, y| x.im.total_cmp(&y.im));
    transformed_poles.sort_by(|x, y| x.im.total_cmp(&y.im));
    for (x, y) in original_poles.iter().zip(transformed_poles) {
        complex(*x, y, 2e-10);
    }
    let force = [C64::new(1.0, 0.2), C64::new(-0.4, 0.7)];
    let transformed_force = [force[0].scale(2.0), force[0] + force[1].scale(3.0)];
    let ra = a.harmonic_response(1.2, &force, 1e-12, || false).unwrap();
    let rb = b
        .harmonic_response(1.2, &transformed_force, 1e-12, || false)
        .unwrap();
    complex(
        ra.displacement[0],
        rb.displacement[0].scale(2.0) + rb.displacement[1],
        2e-12,
    );
    complex(ra.displacement[1], rb.displacement[1].scale(3.0), 2e-12);
    near(ra.mean_input_power_w, rb.mean_input_power_w, 2e-12);
    near(
        ra.mean_velocity_force_power_w,
        rb.mean_velocity_force_power_w,
        2e-12,
    );
}

#[test]
fn fixed_coefficient_campbell_samples_preserve_order_and_cancel_without_partial_results() {
    let model = isotropic(1.0, 0.05, 4.0, 0.4);
    let speeds = [-4.0, 0.0, 2.0, 4.0];
    let samples = model.campbell_samples(&speeds, || false).unwrap();
    for (sample, &speed) in samples.iter().zip(&speeds) {
        assert_eq!(sample.spin_rad_s, speed);
        assert_eq!(
            sample,
            &model
                .at_speed(speed, || false)
                .unwrap()
                .poles(|| false)
                .unwrap()
        );
    }
    assert_eq!(samples, model.campbell_samples(&speeds, || false).unwrap());
    assert!(model.campbell_samples(&[], || false).is_err());
    assert!(
        model
            .campbell_samples(&vec![1.0; MAX_CAMPBELL_SAMPLES + 1], || false)
            .is_err()
    );
    let point = model.at_speed(1.0, || false).unwrap();
    let polls = Cell::new(0);
    assert!(matches!(
        point.poles(|| {
            polls.set(polls.get() + 1);
            polls.get() == 12
        }),
        Err(RotorError::Cancelled)
    ));
    assert_eq!(polls.get(), 12);
    let force = [C64::ONE, C64::new(0.0, -1.0)];
    let before = model.clone();
    polls.set(0);
    assert!(matches!(
        point.harmonic_response(1.3, &force, 1e-12, || {
            polls.set(polls.get() + 1);
            polls.get() == 6
        }),
        Err(RotorError::Cancelled)
    ));
    assert_eq!(model, before);
    assert_eq!(force, [C64::ONE, C64::new(0.0, -1.0)]);
    assert_eq!(
        point
            .harmonic_response(1.3, &force, 1e-12, || false)
            .unwrap(),
        point
            .harmonic_response(1.3, &force, 1e-12, || false)
            .unwrap()
    );
}

#[test]
fn malformed_mass_gyro_loads_and_unrepresentable_response_refuse() {
    let i = [1.0, 0.0, 0.0, 1.0];
    let z = [0.0; 4];
    for n in [0, MAX_ROTOR_DOFS + 1] {
        assert!(RotorBearingModel::new(n, &i, &z, &i, &z, || false).is_err());
    }
    assert!(matches!(
        RotorBearingModel::new(2, &[1.0, 2.0, 2.0, 1.0], &z, &i, &z, || false),
        Err(RotorError::MassNotSpd)
    ));
    assert!(RotorBearingModel::new(2, &[1.0, 0.2, 0.0, 1.0], &z, &i, &z, || false).is_err());
    assert!(RotorBearingModel::new(2, &i, &z, &i, &[0.0, 1.0, 1.0, 0.0], || false).is_err());
    assert!(RotorBearingModel::new(2, &i, &z, &i, &[1.0, 0.0, 0.0, 0.0], || false).is_err());
    assert!(RotorBearingModel::new(2, &i, &[f64::NAN; 4], &i, &z, || false).is_err());
    assert!(matches!(
        RotorBearingModel::new(2, &i, &z, &i, &z, || true),
        Err(RotorError::Cancelled)
    ));
    let model = isotropic(1.0, 0.1, 4.0, 0.5);
    assert!(model.at_speed(f64::INFINITY, || false).is_err());
    let point = model.at_speed(1.0, || false).unwrap();
    assert!(point.harmonic_response(1.0, &[], 1e-12, || false).is_err());
    assert!(
        point
            .harmonic_response(1.0, &[C64::ONE; 2], 0.0, || false)
            .is_err()
    );
    assert!(
        point
            .harmonic_response(f64::NAN, &[C64::ONE; 2], 1e-12, || false)
            .is_err()
    );
    assert!(
        point
            .harmonic_response(1e308, &[C64::ONE; 2], 1e-12, || false)
            .is_err()
    );
}
