//! Original coupled motion at exact quiet-receiver frequencies.
use fs_math::c64::C64;
use fs_mbd::rotordynamics::{
    RotorBearingModel, RotorError,
    relaxation::MaxwellDampingBranch,
    resonance::{AntiresonanceControls, ResonanceSearchError, RotorObservation},
};

fn model() -> RotorBearingModel {
    RotorBearingModel::new(
        2,
        &[1.0, 0.0, 0.0, 1.0],
        &[0.4, 0.0, 0.0, 0.0],
        &[7.0, -4.0, -4.0, 4.0],
        &[0.0; 4],
        || false,
    )
    .unwrap()
}

fn controls(recovery_power: Option<f64>) -> AntiresonanceControls {
    AntiresonanceControls {
        absolute_frequency_tolerance_rad_s: 1e-10,
        relative_frequency_tolerance: 0.0,
        maximum_evaluations: 200,
        recovery_power,
    }
}

fn branch() -> MaxwellDampingBranch {
    MaxwellDampingBranch {
        strain: vec![0.0, 1.0],
        stiffness: 2.0,
        relaxation_time_s: 0.5,
    }
}

// Independent inverse of the physical 2x2 dynamic stiffness, never a production
// response, tangent, spectrum, peak or notch solve.
fn oracle(omega: f64, force: [C64; 2], memory: bool) -> [C64; 2] {
    let d11 = C64::new(7.0 - omega * omega, 0.4 * omega);
    let h = if memory {
        C64::new(0.0, omega) / C64::new(1.0, 0.5 * omega)
    } else {
        C64::ZERO
    };
    let d22 = C64::from_re(4.0 - omega * omega) + h;
    let determinant = d11 * d22 - C64::from_re(16.0);
    [
        (d22 * force[0] + force[1].scale(4.0)) / determinant,
        (force[0].scale(4.0) + d11 * force[1]) / determinant,
    ]
}

fn crossing(mut lo: f64, mut hi: f64, force: [C64; 2], level: f64, memory: bool) -> f64 {
    let negative = oracle(lo, force, memory)[0].norm_sq() < level;
    assert_ne!(negative, oracle(hi, force, memory)[0].norm_sq() < level);
    for _ in 0..80 {
        let middle = 0.5 * lo + 0.5 * hi;
        if (oracle(middle, force, memory)[0].norm_sq() < level) == negative {
            lo = middle;
        } else {
            hi = middle;
        }
    }
    0.5 * lo + 0.5 * hi
}

#[test]
fn original_tuned_absorber_zero_retains_finite_internal_motion_and_recovery_band() {
    let model = model();
    let point = model.at_speed(7.0, || false).unwrap();
    let force = [C64::new(1.0, 0.5), C64::ZERO];
    let level = 0.0005;
    let notch = point
        .driven_antiresonance(
            (1.8, 2.2),
            &force,
            RotorObservation::Displacement(&[1.0, 0.0]),
            2e-12,
            controls(Some(level)),
            || false,
        )
        .unwrap();
    assert!((notch.minimum.angular_frequency_rad_s - 2.0).abs() < 2e-10);
    assert!(notch.observation.abs() < 1e-10);
    assert!((notch.observation_derivative - force[0].scale(0.25)).abs() < 2e-9);
    let state = point
        .harmonic_response(2.0, &force, 2e-12, || false)
        .unwrap();
    assert!(state.displacement[0].abs() < 1e-14);
    assert!((state.displacement[1] + force[0].scale(0.25)).abs() < 1e-13);
    let recovery = notch.recovery.unwrap();
    let lower = crossing(1.8, 2.0, force, level, false);
    let upper = crossing(2.0, 2.2, force, level, false);
    assert!((recovery.lower.angular_frequency_rad_s - lower).abs() < 3e-10);
    assert!((recovery.upper.angular_frequency_rad_s - upper).abs() < 3e-10);
    assert!((recovery.bandwidth_rad_s - (upper - lower)).abs() < 5e-10);
    assert_eq!(recovery.power_level, level);
    assert!(notch.maximum_relative_residual < 2e-12);
}

#[test]
fn hereditary_coherent_notch_keeps_original_branch_motion_stress_and_heat() {
    let model = model();
    let branches = [branch()];
    let point = model
        .at_speed_with_maxwell_memory(-4.0, &branches, || false)
        .unwrap();
    // At omega=2, q=[0,1] exactly; the second fixed force supplies finite
    // material heat while its phase cancels the first receiver's motion.
    let force = [C64::from_re(-4.0), C64::new(1.0, 1.0)];
    let level = 0.0025;
    let notch = point
        .driven_antiresonance(
            (1.9, 2.1),
            &force,
            RotorObservation::Displacement(&[1.0, 0.0]),
            3e-12,
            controls(Some(level)),
            || false,
        )
        .unwrap();
    assert!((notch.minimum.angular_frequency_rad_s - 2.0).abs() < 2e-10);
    assert!(notch.observation.abs() < 2e-10);
    let response = point
        .harmonic_response(2.0, &force, 3e-12, || false)
        .unwrap();
    assert!(response.motion.displacement[0].abs() < 1e-13);
    assert!((response.motion.displacement[1] - C64::ONE).abs() < 2e-13);
    assert!((response.branch_stress[0] - C64::new(1.0, 1.0)).abs() < 3e-13);
    assert!((response.branch_mean_heat_w[0] - 1.0).abs() < 5e-13);
    assert!(response.mean_power_residual_w.abs() < 1e-12);
    let recovery = notch.recovery.unwrap();
    let lower = crossing(1.9, 2.0, force, level, true);
    let upper = crossing(2.0, 2.1, force, level, true);
    assert!((recovery.lower.angular_frequency_rad_s - lower).abs() < 3e-10);
    assert!((recovery.upper.angular_frequency_rad_s - upper).abs() < 3e-10);
    assert!((recovery.bandwidth_rad_s - (upper - lower)).abs() < 5e-10);
    let velocity = point
        .driven_antiresonance(
            (1.9, 2.1),
            &force,
            RotorObservation::Velocity(&[1.0, 0.0]),
            3e-12,
            controls(None),
            || false,
        )
        .unwrap();
    assert!((velocity.minimum.angular_frequency_rad_s - 2.0).abs() < 2e-10);
    assert!(
        (velocity.observation_derivative - C64::new(0.0, 2.0) * notch.observation_derivative).abs()
            < 2e-8
    );
}

#[test]
fn physical_notch_queries_retain_empty_memory_equivalence_and_refusals() {
    let model = model();
    let point = model.at_speed(3.0, || false).unwrap();
    let memory = model
        .at_speed_with_maxwell_memory(3.0, &[], || false)
        .unwrap();
    let force = [C64::ONE, C64::ZERO];
    let receiver = RotorObservation::Displacement(&[1.0, 0.0]);
    let a = point
        .driven_antiresonance((1.8, 2.2), &force, receiver, 2e-12, controls(None), || {
            false
        })
        .unwrap();
    let b = memory
        .driven_antiresonance((1.8, 2.2), &force, receiver, 2e-12, controls(None), || {
            false
        })
        .unwrap();
    assert!((a.minimum.angular_frequency_rad_s - b.minimum.angular_frequency_rad_s).abs() < 2e-10);
    assert!((a.observation - b.observation).abs() < 1e-10);
    assert!(
        point
            .driven_antiresonance(
                (1.8, 2.2),
                &[C64::ZERO; 2],
                receiver,
                2e-12,
                controls(None),
                || false
            )
            .is_err()
    );
    assert!(
        memory
            .driven_antiresonance(
                (1.8, 2.2),
                &force,
                RotorObservation::Displacement(&[1.0]),
                2e-12,
                controls(None),
                || false
            )
            .is_err()
    );
    assert!(
        point
            .driven_antiresonance(
                (1.8, 2.2),
                &force,
                receiver,
                2e-12,
                controls(Some(1e4)),
                || false
            )
            .is_err()
    );
    assert!(matches!(
        memory.driven_antiresonance((1.8, 2.2), &force, receiver, 2e-12, controls(None), || true),
        Err(ResonanceSearchError::Source(RotorError::Cancelled))
    ));
    let mut small = controls(None);
    small.maximum_evaluations = 2;
    assert!(
        point
            .driven_antiresonance((1.8, 2.2), &force, receiver, 2e-12, small, || false)
            .is_err()
    );
}
