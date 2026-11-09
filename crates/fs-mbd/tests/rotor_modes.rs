//! Independent physical checks of original-coordinate rotor modes and residues.

use fs_math::c64::C64;
use fs_mbd::rotordynamics::{
    RotorBearingModel,
    resonance::poles::{
        PoleOptions, PoleResponseOptions, QuadraticObservation, QuadraticResponseOptions,
        RotorModeError,
    },
};

fn pole_options() -> PoleOptions {
    PoleOptions {
        spectral_scale: 3.0,
        step_tolerance: 2e-13,
        residual_tolerance: 2e-12,
        max_iterations: 40,
        max_dimension: 8,
        max_dense_bytes: 16 * 1024 * 1024,
        max_work: 100_000_000,
    }
}

fn response_options() -> QuadraticResponseOptions {
    QuadraticResponseOptions {
        numerical: PoleResponseOptions {
            spectral_scale: 3.0,
            pole_residual_tolerance: 3e-12,
            pole_step_tolerance: 3e-12,
            linear_residual_tolerance: 3e-12,
            max_dimension: 8,
            max_vectors: 4,
            max_dense_bytes: 16 * 1024 * 1024,
            max_work: 100_000_000,
        },
        max_observations: 4,
        max_physical_bytes: 16 * 1024 * 1024,
        max_physical_work: 100_000_000,
    }
}

fn close(actual: C64, expected: C64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * (1.0 + expected.abs()),
        "actual={actual:?}, expected={expected:?}"
    );
}

#[test]
fn signed_gyro_whirl_shapes_residues_and_spin_slopes_match_scalar_polynomials() {
    let model = RotorBearingModel::new(
        2,
        &[1.0, 0.0, 0.0, 1.0],
        &[0.0; 4],
        &[4.0, 0.0, 0.0, 4.0],
        &[0.0, -1.0, 1.0, 0.0],
        || false,
    )
    .unwrap();
    for spin in [-3.0_f64, 3.0] {
        let point = model.at_speed(spin, || false).unwrap();
        let prepared = point.prepare_resonance_modes(pole_options(), || false).unwrap();
        for handedness in [-1.0_f64, 1.0] {
            // J*[1,-i*h] = i*h*[1,-i*h]. The original polynomial
            // reduces exactly to s^2 + i*h*Omega*s + 4.
            let frequency = (-handedness * spin + (spin * spin + 16.0).sqrt()) / 2.0;
            let pole = C64::new(0.0, frequency);
            let shape = [C64::ONE, C64::new(0.0, -handedness)];
            let mode = prepared
                .refine(pole + C64::new(0.013, 0.021), &shape, || false)
                .unwrap();
            close(mode.physical.laplace_frequency_per_s, pole, 2e-11);
            assert!(mode.physical.relative_residual < 2e-12);
            // Equal-magnitude circular coordinates may select either numerical
            // phase pivot. Unit mass and the physical relative phase are invariant.
            let mass: f64 = mode.physical.displacement.iter().map(|q| q.norm_sq()).sum();
            assert!((mass - 1.0).abs() < 2e-12);
            close(mode.physical.displacement[1] / mode.physical.displacement[0], shape[1], 2e-11);
            for (i, &coordinate) in mode.physical.displacement.iter().enumerate() {
                close(mode.physical.velocity[i], pole * coordinate, 2e-11);
            }
            let factor = prepared.prepare_response(&mode, response_options(), || false).unwrap();
            let derivative = factor.spin_derivative(|| false).unwrap();
            let denominator = pole.scale(2.0) + C64::new(0.0, handedness * spin);
            let expected_slope = -(C64::new(0.0, handedness) * pole) / denominator;
            close(derivative.laplace_frequency_per_s, expected_slope, 2e-11);
            assert!(derivative.polynomial_relative_residual < 3e-12);
            assert!(derivative.normalization_relative_residual < 3e-12);
            for (i, &coordinate) in factor.response().right_displacement().iter().enumerate() {
                close(derivative.displacement[i], C64::ZERO, 2e-11);
                close(derivative.velocity[i], expected_slope * coordinate, 2e-11);
            }
            // A circular source reduces to the same scalar oscillator, with
            // residue F/(2*s+i*h*Omega). Reuse the sensitivity's SAME factor.
            let residues = factor.forcing_residues(&[&shape], &[], || false).unwrap();
            for (i, &coordinate) in shape.iter().enumerate() {
                close(residues[0].displacement[i], coordinate / denominator, 2e-11);
                close(residues[0].velocity[i], pole * coordinate / denominator, 2e-11);
            }
        }
    }
}

fn triangular_model(scales: [f64; 2]) -> RotorBearingModel {
    let transform = |matrix: [f64; 4]| -> Vec<f64> {
        (0..4)
            .map(|ij| scales[ij / 2] * matrix[ij] * scales[ij % 2])
            .collect()
    };
    RotorBearingModel::new(
        2,
        &transform([1.0, 0.0, 0.0, 1.0]),
        &transform([0.2, 0.0, 0.0, 0.2]),
        &transform([4.0, 3.0, 0.0, 9.0]),
        &[0.0; 4],
        || false,
    )
    .unwrap()
}

#[test]
fn nonnormal_left_participation_excites_a_right_orthogonal_force_and_dark_receivers_stay_dark() {
    let model = triangular_model([1.0, 1.0]);
    let point = model.at_speed(-7.0, || false).unwrap();
    let prepared = point.prepare_resonance_modes(pole_options(), || false).unwrap();
    let root_imaginary = 3.99_f64.sqrt();
    let pole = C64::new(-0.1, root_imaginary);
    let mode = prepared
        .refine(
            pole + C64::new(0.02, -0.01),
            &[C64::ONE, C64::new(0.03, 0.01)],
            || false,
        )
        .unwrap();
    close(mode.physical.displacement[0], C64::ONE, 2e-11);
    close(mode.physical.displacement[1], C64::ZERO, 2e-11);
    let factor = prepared.prepare_response(&mode, response_options(), || false).unwrap();
    let e1 = [C64::ONE, C64::ZERO];
    let e2 = [C64::ZERO, C64::ONE];
    let zero = [C64::ZERO; 2];
    let displacement = [C64::ONE, C64::new(0.0, -0.4)];
    let velocity = [C64::new(0.2, 0.1), C64::from_re(0.3)];
    let observations = [
        QuadraticObservation { displacement: &displacement, velocity: &velocity },
        QuadraticObservation { displacement: &e2, velocity: &zero },
    ];
    let forces: [&[C64]; 3] = [&e1, &e2, &zero];
    let residues = factor.forcing_residues(&forces, &observations, || false).unwrap();
    let denominator = pole.scale(2.0) + C64::from_re(0.2);
    // Inverting [[a,3],[0,a+5]] gives q1=-3/[a*(a+5)] for F=e2.
    // Its first-mode residue is -3/[5*a'(s0)], although phi^H F=0.
    for (result, numerator) in residues.iter().zip([1.0, -0.6, 0.0]) {
        let expected = C64::from_re(numerator) / denominator;
        close(result.displacement[0], expected, 2e-11);
        close(result.displacement[1], C64::ZERO, 2e-11);
        close(result.velocity[0], pole * expected, 2e-11);
        close(result.participation, expected, 2e-11);
        close(result.observations[0], (displacement[0] + velocity[0] * pole) * expected, 2e-11);
        close(result.observations[1], C64::ZERO, 2e-11);
        assert!(result.linear_residual < 3e-12);
    }
    assert!(residues[1].displacement[0].abs() > 0.1);
    let omega_residue = C64::new(residues[1].displacement[0].im, -residues[1].displacement[0].re);
    close(omega_residue, C64::from_re(0.3 / root_imaginary), 2e-11);
    close(mode.angular_frequency_rad_s(), C64::new(root_imaginary, 0.1), 2e-11);
    let slope = factor.spin_derivative(|| false).unwrap();
    close(slope.laplace_frequency_per_s, C64::ZERO, 2e-11);
    for value in slope.displacement.iter().chain(&slope.velocity) {
        close(*value, C64::ZERO, 2e-11);
    }
    let replay = factor.forcing_residues(&forces, &observations, || false).unwrap();
    for (first, second) in residues.iter().zip(replay) {
        assert_eq!(first.displacement, second.displacement);
        assert_eq!(first.observations, second.observations);
    }
}

#[test]
fn original_force_and_sensor_residues_are_invariant_under_signed_coordinate_congruence() {
    let pole = C64::new(-0.1, 3.99_f64.sqrt());
    let physical_residue = C64::from_re(-0.6) / (pole.scale(2.0) + C64::from_re(0.2));
    for scales in [[-2.0, 8.0], [-1_048_576.0, 1.0 / 1_048_576.0]] {
        // q_physical=S*q_new: coefficients S^T A S, force S^T F,
        // and receiver l^T S preserve mechanical work and measured motion.
        let model = triangular_model(scales);
        let point = model.at_speed(2.0, || false).unwrap();
        let prepared = point.prepare_resonance_modes(pole_options(), || false).unwrap();
        let mode = prepared.refine(pole, &[C64::from_re(1.0 / scales[0]), C64::ZERO], || false).unwrap();
        let modal_mass = mode.physical.displacement[0].norm_sq() * scales[0] * scales[0]
            + mode.physical.displacement[1].norm_sq() * scales[1] * scales[1];
        assert!((modal_mass - 1.0).abs() < 2e-12);
        let factor = prepared.prepare_response(&mode, response_options(), || false).unwrap();
        let force = [C64::ZERO, C64::from_re(scales[1])];
        let receiver = [C64::from_re(scales[0]), C64::ZERO];
        let zero = [C64::ZERO; 2];
        let residues = factor.forcing_residues(
            &[&force],
            &[QuadraticObservation { displacement: &receiver, velocity: &zero }],
            || false,
        ).unwrap();
        close(residues[0].displacement[0].scale(scales[0]), physical_residue, 2e-11);
        close(residues[0].displacement[1].scale(scales[1]), C64::ZERO, 2e-11);
        close(residues[0].observations[0], physical_residue, 2e-11);
    }
}

#[test]
fn complete_scalar_residues_reconstruct_plus_i_harmonics_and_retain_active_growth() {
    let force = [C64::new(1.2, -0.7)];
    for damping in [-2.0_f64, 2.0] {
        let model = RotorBearingModel::new(1, &[4.0], &[damping], &[36.0], &[0.0], || false).unwrap();
        let point = model.at_speed(7.0, || false).unwrap();
        let prepared = point.prepare_resonance_modes(pole_options(), || false).unwrap();
        let real = -damping / 8.0;
        let imaginary = (9.0 - real * real).sqrt();
        let mut residues = Vec::new();
        for sign in [-1.0, 1.0] {
            let pole = C64::new(real, sign * imaginary);
            let mode = prepared.refine(
                pole + C64::new(0.01, -0.02),
                &[C64::new(1.0, 0.3)],
                || false,
            ).unwrap();
            close(mode.physical.laplace_frequency_per_s, pole, 2e-11);
            close(mode.physical.displacement[0], C64::from_re(0.5), 2e-11);
            assert_eq!(mode.physical.laplace_frequency_per_s.re > 0.0, damping < 0.0);
            assert_eq!(mode.angular_frequency_rad_s().im > 0.0, damping > 0.0);
            let factor = prepared.prepare_response(&mode, response_options(), || false).unwrap();
            let mut result = factor.forcing_residues(&[&force], &[], || false).unwrap();
            let result = result.pop().unwrap();
            close(result.displacement[0], force[0] / (pole.scale(8.0) + C64::from_re(damping)), 2e-11);
            residues.push(result);
        }
        for frequency in [0.0, 1.0, 2.95, 5.0] {
            let harmonic = point.harmonic_response(frequency, &force, 2e-12, || false).unwrap();
            let s = C64::new(0.0, frequency);
            let mut displacement = C64::ZERO;
            let mut velocity = C64::ZERO;
            for residue in &residues {
                displacement = displacement + residue.displacement[0] / (s - residue.laplace_pole_per_s);
                velocity = velocity + residue.velocity[0] / (s - residue.laplace_pole_per_s);
            }
            close(displacement, harmonic.displacement[0], 2e-10);
            close(velocity, harmonic.velocity[0], 2e-10);
        }
    }
}

#[test]
fn source_rechecks_resource_caps_multiplicity_and_cancellation_refuse_without_poisoning_factor() {
    let model = RotorBearingModel::new(1, &[1.0], &[0.2], &[4.0], &[0.0], || false).unwrap();
    let point = model.at_speed(1.0, || false).unwrap();
    let mut tight = pole_options();
    tight.max_dense_bytes = 1;
    assert!(matches!(point.prepare_resonance_modes(tight, || false), Err(RotorModeError::Mode(_))));
    assert!(matches!(point.prepare_resonance_modes(pole_options(), || true), Err(RotorModeError::Cancelled)));
    let prepared = point.prepare_resonance_modes(pole_options(), || false).unwrap();
    let pole = C64::new(-0.1, 3.99_f64.sqrt());
    assert!(prepared.refine(pole, &[C64::ONE, C64::ONE], || false).is_err());
    let mut polls = 0;
    assert!(matches!(
        prepared.refine(pole, &[C64::ONE], || { polls += 1; polls >= 4 }),
        Err(RotorModeError::Cancelled)
    ));
    let mode = prepared.refine(pole, &[C64::ONE], || false).unwrap();
    let mut wrong_speed = mode.clone();
    wrong_speed.spin_rad_s += 1.0;
    assert!(matches!(
        prepared.prepare_response(&wrong_speed, response_options(), || false),
        Err(RotorModeError::Invalid(_))
    ));
    let mut wrong_pole = mode.clone();
    wrong_pole.physical.laplace_frequency_per_s.re += 0.2;
    assert!(prepared.prepare_response(&wrong_pole, response_options(), || false).is_err());
    let mut options = response_options();
    options.numerical.max_vectors = 1;
    let factor = prepared.prepare_response(&mode, options, || false).unwrap();
    let force = [C64::ONE];
    let before = factor.forcing_residues(&[&force], &[], || false).unwrap();
    assert!(factor.forcing_residues(&[&force, &force], &[], || false).is_err());
    assert!(matches!(factor.forcing_residues(&[&force], &[], || true), Err(RotorModeError::Cancelled)));
    assert!(matches!(factor.spin_derivative(|| true), Err(RotorModeError::Cancelled)));
    let after = factor.forcing_residues(&[&force], &[], || false).unwrap();
    assert_eq!(before[0].displacement, after[0].displacement);
    assert_eq!(before[0].velocity, after[0].velocity);
    // Exactly repeated undamped roots give a singular normalized border.
    let repeated_model = RotorBearingModel::new(
        2, &[1.0, 0.0, 0.0, 1.0], &[0.0; 4], &[4.0, 0.0, 0.0, 4.0], &[0.0; 4], || false,
    ).unwrap();
    let repeated_point = repeated_model.at_speed(0.0, || false).unwrap();
    let repeated = repeated_point.prepare_resonance_modes(pole_options(), || false).unwrap();
    assert!(repeated.refine(C64::new(0.0, 2.0), &[C64::ONE, C64::ZERO], || false).is_err());
}
