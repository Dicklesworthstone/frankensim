//! G1/G3 nonlinear motion: an independent Duffing period, work and momentum.
use super::*;
use fs_solid::membrane::dynamics::MembraneMotion;

fn initial(stretch: f64, amplitude: f64) -> Vec<[f64; 3]> {
    let mut u = vec![[0.0; 3]; 5];
    for (dof, value) in prescribed(stretch) {
        u[dof / 3][dof % 3] = value;
    }
    u[0][2] = amplitude;
    u
}

// Integrate the exact energy first integral, not another time-stepping scheme.
// m q'' + k q + beta q^3 = 0 for the existing four-triangle pyramid.
fn exact_period(stretch: f64, amplitude: f64) -> f64 {
    let m = 2.0 * 0.2 * 0.1_f64.powi(2) / 3.0;
    let k = 2.0 * 8000.0 * (stretch * stretch - 1.0);
    let beta = 4.0 * 8000.0 / 0.1_f64.powi(2);
    let intervals = 4096;
    let h = std::f64::consts::FRAC_PI_2 / intervals as f64;
    let integrand = |theta: f64| {
        1.0 / (k + 0.5 * beta * amplitude * amplitude * (1.0 + theta.sin().powi(2))).sqrt()
    };
    let sum: f64 = (0..=intervals)
        .map(|i| {
            let weight = if i == 0 || i == intervals {
                1.0
            } else if i % 2 == 0 {
                2.0
            } else {
                4.0
            };
            weight * integrand(i as f64 * h)
        })
        .sum();
    4.0 * m.sqrt() * h * sum / 3.0
}

fn crossing_period(stretch: f64, amplitude: f64, subdivisions: usize) -> (f64, f64) {
    let p = patch(8000.0, 8000.0, 4000.0, 0.2, 0.0);
    let zero = vec![[0.0; 3]; 5];
    let mut motion = MembraneMotion::new(
        &p,
        &prescribed(stretch),
        &initial(stretch, amplitude),
        &zero,
        0.0,
        || false,
    )
    .unwrap();
    let period = exact_period(stretch, amplitude);
    let dt = period / subdivisions as f64;
    let energy = motion.strain_energy_j();
    let baseline = p
        .evaluate(&initial(stretch, 0.0), || false)
        .unwrap()
        .strain_energy_j;
    let mut maximum_defect = 0.0_f64;
    for _ in 0..subdivisions {
        let old_q = motion.displacement_m()[0][2];
        let old_time = motion.time_s();
        let step = motion.step(dt, &zero, || false).unwrap();
        maximum_defect =
            maximum_defect.max((step.strain_energy_j + step.kinetic_energy_j - energy).abs());
        let new_q = motion.displacement_m()[0][2];
        if new_q <= 0.0 {
            return (
                4.0 * (old_time + dt * old_q / (old_q - new_q)),
                maximum_defect / (energy - baseline),
            );
        }
    }
    panic!("elastic membrane failed to complete its first quarter-cycle");
}

#[test]
fn nonlinear_membrane_period_converges_to_energy_integral_at_second_order() {
    for stretch in [1.0, 1.02] {
        let exact = exact_period(stretch, 0.008);
        let mut errors = Vec::new();
        for subdivisions in [80, 160, 320] {
            let (period, energy_error) = crossing_period(stretch, 0.008, subdivisions);
            errors.push((period / exact - 1.0).abs());
            assert!(
                energy_error < 0.004,
                "unresolved elastic energy: {energy_error}"
            );
        }
        assert!(
            errors[0] / errors[1] > 3.5 && errors[1] / errors[2] > 3.5,
            "{errors:?}"
        );
        assert!(errors[2] < 1e-4, "{errors:?}");
    }
    // With no prestretch the tangent is zero, but finite-amplitude cubic
    // elasticity still restores motion. Doubling amplitude halves its period.
    let (small, _) = crossing_period(1.0, 0.004, 320);
    let (large, _) = crossing_period(1.0, 0.008, 320);
    close(small, 2.0 * large, 1e-10, 1e-10);
}

#[test]
fn source_stiffness_and_reference_mass_set_the_nonlinear_time_scale() {
    let base = patch(8000.0, 8000.0, 4000.0, 0.2, 0.0);
    let stiff = patch(32000.0, 32000.0, 16000.0, 0.2, 0.0);
    let heavy = patch(8000.0, 8000.0, 4000.0, 0.8, 0.0);
    let u = initial(1.02, 0.006);
    let zero = vec![[0.0; 3]; 5];
    let fixed = prescribed(1.02);
    let mut a = MembraneMotion::new(&base, &fixed, &u, &zero, 0.0, || false).unwrap();
    let mut b = MembraneMotion::new(&stiff, &fixed, &u, &zero, 0.0, || false).unwrap();
    let mut c = MembraneMotion::new(&heavy, &fixed, &u, &zero, 0.0, || false).unwrap();
    for _ in 0..200 {
        a.step(1e-5, &zero, || false).unwrap();
        b.step(0.5e-5, &zero, || false).unwrap();
        c.step(2e-5, &zero, || false).unwrap();
        close(
            a.displacement_m()[0][2],
            b.displacement_m()[0][2],
            1e-14,
            1e-13,
        );
        close(
            a.displacement_m()[0][2],
            c.displacement_m()[0][2],
            1e-14,
            1e-13,
        );
        close(
            2.0 * a.velocity_m_s()[0][2],
            b.velocity_m_s()[0][2],
            1e-12,
            1e-13,
        );
        close(
            0.5 * a.velocity_m_s()[0][2],
            c.velocity_m_s()[0][2],
            1e-12,
            1e-13,
        );
    }
    assert_eq!(base.reference_mass_kg(), stiff.reference_mass_kg());
    assert_eq!(4.0 * base.reference_mass_kg(), heavy.reference_mass_kg());
}

fn momenta(motion: &MembraneMotion<'_>) -> ([f64; 3], [f64; 3]) {
    let mut linear = [0.0; 3];
    let mut angular = [0.0; 3];
    for (i, &mass) in motion.patch().nodal_mass_kg().iter().enumerate() {
        let p = motion.velocity_m_s()[i].map(|v| mass * v);
        let u = motion.displacement_m()[i];
        let reference = motion.patch().mesh().nodes[i];
        let x = [reference[0] + u[0], reference[1] + u[1], u[2]];
        for axis in 0..3 {
            linear[axis] += p[axis];
        }
        angular[0] += x[1] * p[2] - x[2] * p[1];
        angular[1] += x[2] * p[0] - x[0] * p[2];
        angular[2] += x[0] * p[1] - x[1] * p[0];
    }
    (linear, angular)
}

#[test]
fn free_motion_preserves_linear_and_angular_momentum_and_reverses() {
    let p = patch(8000.0, 2000.0, 500.0, 0.2, 0.37);
    let u = displaced();
    let v = vec![
        [0.1, 0.2, -0.1],
        [0.0, 0.1, 0.0],
        [0.2, 0.0, 0.1],
        [-0.1, 0.2, 0.0],
        [0.0, -0.1, 0.2],
    ];
    let zero = vec![[0.0; 3]; 5];
    let mut motion = MembraneMotion::new(&p, &[], &u, &v, 0.0, || false).unwrap();
    let before = momenta(&motion);
    for _ in 0..100 {
        motion.step(1e-6, &zero, || false).unwrap();
    }
    let after = momenta(&motion);
    for axis in 0..3 {
        close(after.0[axis], before.0[axis], 2e-15, 1e-12);
        close(after.1[axis], before.1[axis], 2e-16, 1e-12);
    }
    let mut resumed = motion.clone();
    assert_eq!(
        motion.step(1e-6, &zero, || false).unwrap(),
        resumed.step(1e-6, &zero, || false).unwrap()
    );
    assert_eq!(motion.displacement_m(), resumed.displacement_m());
    for _ in 0..101 {
        motion.step(-1e-6, &zero, || false).unwrap();
    }
    for (a, b) in motion
        .displacement_m()
        .iter()
        .flatten()
        .zip(u.iter().flatten())
    {
        close(*a, *b, 1e-14, 1e-12);
    }
    for (a, b) in motion
        .velocity_m_s()
        .iter()
        .flatten()
        .zip(v.iter().flatten())
    {
        close(*a, *b, 1e-11, 1e-11);
    }
}

#[test]
fn held_load_work_and_support_impulse_close_the_actual_step() {
    let p = patch(8000.0, 8000.0, 4000.0, 0.2, 0.0);
    let u = initial(1.02, 0.006);
    let zero = vec![[0.0; 3]; 5];
    let fixed = prescribed(1.02);
    let mut motion = MembraneMotion::new(&p, &fixed, &u, &zero, 0.0, || false).unwrap();
    let mut loads = zero;
    loads[0][2] = -2.0;
    loads[1][0] = 5.0; // a loaded stationary support does no work
    let before = momenta(&motion).0;
    let old_energy = motion.strain_energy_j() + motion.kinetic_energy_j();
    let step = motion.step(1e-5, &loads, || false).unwrap();
    let after = momenta(&motion).0;
    close(
        step.dead_load_work_j,
        -2.0 * (motion.displacement_m()[0][2] - u[0][2]),
        0.0,
        1e-14,
    );
    close(
        step.energy_defect_j,
        motion.strain_energy_j() + motion.kinetic_energy_j() - old_energy - step.dead_load_work_j,
        1e-15,
        1e-10,
    );
    for axis in 0..3 {
        let applied_impulse = 1e-5 * loads.iter().map(|f| f[axis]).sum::<f64>();
        close(
            after[axis] - before[axis],
            applied_impulse + step.support_impulse_n_s[axis],
            1e-15,
            1e-10,
        );
    }
    for (dof, value) in fixed {
        assert_eq!(motion.displacement_m()[dof / 3][dof % 3], value);
        assert_eq!(motion.velocity_m_s()[dof / 3][dof % 3], 0.0);
    }
}

#[test]
fn refused_trials_leave_motion_energy_and_checkpoint_unchanged() {
    let p = patch(8000.0, 8000.0, 4000.0, 0.2, 0.0);
    let u = initial(1.02, 0.006);
    let zero = vec![[0.0; 3]; 5];
    let fixed = prescribed(1.02);
    let initial = MembraneMotion::new(&p, &fixed, &u, &zero, 0.0, || false).unwrap();
    for cancel_at in [1, 10, 20, 30] {
        let mut motion = initial.clone();
        let mut calls = 0;
        let result = motion.step(1e-5, &zero, || {
            calls += 1;
            calls == cancel_at
        });
        assert!(
            matches!(result, Err(MembraneError::Cancelled)),
            "poll {cancel_at}: {result:?}"
        );
        assert_eq!(motion.displacement_m(), initial.displacement_m());
        assert_eq!(motion.velocity_m_s(), initial.velocity_m_s());
        assert_eq!(motion.time_s(), initial.time_s());
        assert_eq!(motion.strain_energy_j(), initial.strain_energy_j());
        assert_eq!(motion.kinetic_energy_j(), initial.kinetic_energy_j());
        let mut clean = initial.clone();
        assert_eq!(
            motion.step(1e-5, &zero, || false).unwrap(),
            clean.step(1e-5, &zero, || false).unwrap()
        );
    }
    let mut motion = initial.clone();
    assert!(motion.step(1.0, &zero, || false).is_err());
    assert_eq!(motion.displacement_m(), initial.displacement_m());
    assert_eq!(motion.velocity_m_s(), initial.velocity_m_s());
    assert!(motion.step(f64::NAN, &zero, || false).is_err());
    let mut invalid_velocity = zero;
    invalid_velocity[1][2] = 1.0;
    assert!(MembraneMotion::new(&p, &fixed, &u, &invalid_velocity, 0.0, || false).is_err());
}
