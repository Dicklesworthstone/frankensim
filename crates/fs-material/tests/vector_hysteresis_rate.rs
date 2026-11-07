//! Intrinsic finite-rate switching on the original energetic vector owner.
//! Scaled analytic fixtures below are not measured magnetic material claims.

use fs_material::field_storage::hysteresis::{
    HysteresisBranch,
    vector::{
        VectorHysteresis, VectorHysteresisControls, VectorHysteresisError, VectorHysteresisState,
        VectorHysteresisStep,
    },
};

type Vector = [f64; 3];

fn branch(chi: f64, pinning: f64, saturation: f64) -> HysteresisBranch {
    HysteresisBranch {
        susceptibility_h_m: chi,
        pinning_field_a_m: pinning,
        saturation_polarization_t: saturation,
    }
}
fn controls() -> VectorHysteresisControls {
    VectorHysteresisControls {
        relative_tolerance: 1e-12,
        absolute_field_tolerance_a_m: 1e-12,
        absolute_flux_tolerance_t: 1e-12,
        work_tolerance_j_m3: 2e-11,
        maximum_branch_iterations: 100,
        maximum_inverse_iterations: 60,
        maximum_backtracks: 40,
    }
}
fn law(pinning: f64, viscosity: f64) -> VectorHysteresis {
    VectorHysteresis::new(0.01, 30.0, vec![branch(0.2, pinning, 1.0)])
        .unwrap()
        .with_viscosities_m_s_h(vec![viscosity])
        .unwrap()
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.16e} != {expected:.16e}; tolerance {tolerance:e}"
    );
}
fn vector_near(actual: Vector, expected: Vector, tolerance: f64) {
    for i in 0..3 {
        near(actual[i], expected[i], tolerance);
    }
}
fn dot(a: Vector, b: Vector) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}
fn norm(a: Vector) -> f64 {
    a[0].hypot(a[1]).hypot(a[2])
}
fn sub(a: Vector, b: Vector) -> Vector {
    std::array::from_fn(|i| a[i] - b[i])
}
fn step(
    material: &VectorHysteresis,
    previous: &VectorHysteresisState,
    field: Vector,
    duration: f64,
) -> VectorHysteresisStep {
    material
        .candidate_field_timed(previous, field, duration, controls(), &mut || false)
        .unwrap()
}
fn audit(step: &VectorHysteresisStep) {
    assert!(step.pinning_heat_j_m3 >= 0.0);
    assert!(step.viscous_heat_j_m3 >= 0.0);
    assert!(step.numerical_dissipation_j_m3 >= 0.0);
    assert_eq!(
        step.irreversible_heat_j_m3,
        step.pinning_heat_j_m3 + step.viscous_heat_j_m3
    );
    near(
        step.endpoint_work_j_m3,
        step.energy_change_j_m3 + step.irreversible_heat_j_m3 + step.numerical_dissipation_j_m3,
        3e-11,
    );
    near(step.energy_residual_j_m3, 0.0, 3e-11);
}

#[test]
fn zero_viscosity_preserves_original_history_and_inverse_bit_for_bit() {
    let original = VectorHysteresis::new(
        0.01,
        30.0,
        vec![branch(0.2, 2.0, 1.0), branch(0.1, 0.0, 0.5)],
    )
    .unwrap();
    let zero = original
        .clone()
        .with_viscosities_m_s_h(vec![0.0, 0.0])
        .unwrap();
    assert_eq!(zero, original);
    assert!(!zero.is_rate_dependent());
    assert_eq!(zero.viscosities_m_s_h(), None);
    let mut previous = original.initial_state();
    for (index, h) in [
        [0.5, 0.0, 0.0],
        [5.0, 0.0, 0.0],
        [3.0, 7.0, 0.5],
        [-8.0, 3.0, 1.0],
        [0.0, 0.0, 0.0],
        [2.0, -4.0, 0.5],
    ]
    .into_iter()
    .enumerate()
    {
        let duration = 0.1 + index as f64;
        let expected = original
            .candidate_field(&previous, h, controls(), &mut || false)
            .unwrap();
        let actual = step(&zero, &previous, h, duration);
        assert_eq!(actual, expected);
        let inverse = original
            .candidate_flux_density(
                &previous,
                actual.state.flux_density_t(),
                controls(),
                &mut || false,
            )
            .unwrap();
        let timed_inverse = zero
            .candidate_flux_density_timed(
                &previous,
                actual.state.flux_density_t(),
                duration,
                controls(),
                &mut || false,
            )
            .unwrap();
        assert_eq!(inverse, timed_inverse);
        assert_eq!(actual.viscous_heat_j_m3, 0.0);
        previous = actual.state;
    }
}

#[test]
fn collinear_holds_reversals_and_rate_changes_match_backward_euler_kinetics() {
    let material = law(2.0, 2.0);
    let mut previous = material.initial_state();
    let mut expected_z: f64 = 0.0;
    for (h, duration) in [(5.0, 0.1), (5.0, 0.1), (5.0, 0.4), (2.0, 0.7), (-5.0, 0.2)] {
        let old_z = expected_z;
        let drive = h - old_z / 0.2;
        if drive.abs() > 2.0 {
            // Independent scalar backward-Euler recurrence on the active
            // dry-friction side: (eta + dt/chi) z = eta z_old + dt (H +/- hc).
            expected_z =
                (2.0 * old_z + duration * (h - 2.0 * drive.signum())) / (2.0 + duration / 0.2);
        }
        let candidate = step(&material, &previous, [h, 0.0, 0.0], duration);
        let dz = expected_z - old_z;
        vector_near(
            candidate.state.polarizations_t()[0],
            [expected_z, 0.0, 0.0],
            2e-14,
        );
        near(candidate.pinning_heat_j_m3, 2.0 * dz.abs(), 2e-14);
        near(candidate.viscous_heat_j_m3, 2.0 * dz * dz / duration, 2e-14);
        near(
            candidate.state.stored_energy_j_m3(),
            0.005 * h * h + expected_z * expected_z / 0.4,
            2e-14,
        );
        let dh = h - previous.field_a_m()[0];
        near(
            candidate.numerical_dissipation_j_m3,
            0.005 * dh * dh + dz * dz / 0.4,
            2e-14,
        );
        near(
            candidate.incremental_potential_change_j_m3().unwrap(),
            candidate.energy_change_j_m3
                + candidate.pinning_heat_j_m3
                + 0.5 * candidate.viscous_heat_j_m3,
            2e-14,
        );
        audit(&candidate);
        previous = candidate.state;
    }
}

#[test]
fn fixed_field_relaxation_and_physical_heat_converge_to_continuous_solution() {
    let material = law(2.0, 2.0);
    let final_time: f64 = 0.8;
    let relaxation_time: f64 = 2.0 * 0.2;
    let equilibrium: f64 = 0.2 * (5.0 - 2.0);
    let exact_z = equilibrium * (1.0 - (-final_time / relaxation_time).exp());
    let exact_viscous =
        equilibrium * equilibrium / 0.4 * (1.0 - (-2.0 * final_time / relaxation_time).exp());
    let mut errors = vec![];
    for intervals in [8, 16, 32] {
        let duration = final_time / f64::from(intervals);
        let mut previous = material.initial_state();
        let (mut pinning, mut viscous, mut numerical) = (0.0, 0.0, 0.0);
        for _ in 0..intervals {
            let candidate = step(&material, &previous, [5.0, 0.0, 0.0], duration);
            assert!(candidate.state.polarizations_t()[0][0] > previous.polarizations_t()[0][0]);
            audit(&candidate);
            pinning += candidate.pinning_heat_j_m3;
            viscous += candidate.viscous_heat_j_m3;
            numerical += candidate.numerical_dissipation_j_m3;
            previous = candidate.state;
        }
        let z = previous.polarizations_t()[0][0];
        near(pinning, 2.0 * z, 3e-14);
        // The initial discontinuous background-field step contributes 0.125
        // independently of the hold integration. Only branch remainder tends
        // to zero as the physical relaxation is temporally resolved.
        errors.push((
            (z - exact_z).abs(),
            (viscous - exact_viscous).abs(),
            numerical - 0.125,
        ));
    }
    for pair in errors.windows(2) {
        assert!(pair[1].0 < 0.56 * pair[0].0);
        assert!(pair[1].1 < 0.56 * pair[0].1);
        assert!(pair[1].2 < 0.56 * pair[0].2);
    }
}

#[test]
fn smooth_zero_drive_branch_retains_finite_rate_compliance() {
    let material = law(0.0, 2.0);
    let result = step(&material, &material.initial_state(), [0.0; 3], 0.5);
    assert_eq!(result.state.polarizations_t()[0], [0.0; 3]);
    assert_eq!(result.irreversible_heat_j_m3, 0.0);
    for i in 0..3 {
        for j in 0..3 {
            near(
                result.incremental_permeability_h_m[i][j],
                if i == j {
                    0.01 + 1.0 / (5.0 + 2.0 / 0.5)
                } else {
                    0.0
                },
                2e-16,
            );
        }
    }
}

#[test]
fn rotating_saturated_kinetics_has_consistent_tangent_inverse_and_half_viscous_potential() {
    for pinning in [0.0, 2.0] {
        let material = law(pinning, 0.4);
        let previous = step(&material, &material.initial_state(), [9.0, 0.0, 0.0], 1.0).state;
        near(norm(previous.polarizations_t()[0]), 1.0, 2e-14);
        let saved = previous.clone();
        let duration = 0.2;
        for h in [[6.0, 0.0, 0.0], [3.0, 3.0, 0.5], [3.0, 8.0, 0.5]] {
            let candidate = step(&material, &previous, h, duration);
            audit(&candidate);
            let z = candidate.state.polarizations_t()[0];
            let dz = sub(z, previous.polarizations_t()[0]);
            let distance = norm(dz);
            if distance > 1e-12 {
                // Independent KKT check in original physical units. Viscosity
                // does not change chi, hc, the saturation radius, or storage.
                let force = std::array::from_fn(|i| {
                    h[i] - z[i] / 0.2 - (0.4 / duration) * dz[i] - pinning * dz[i] / distance
                });
                let reaction = if norm(z) > 1.0 - 1e-12 {
                    dot(force, z)
                } else {
                    0.0
                };
                assert!(reaction >= -3e-11);
                assert!(norm(std::array::from_fn(|i| force[i] - reaction * z[i])) < 3e-11);
            }
            assert!(norm(z) <= 1.0 + 2e-14);
            let inverse = candidate.incremental_reluctivity_m_h().unwrap();
            let b = candidate.state.flux_density_t();
            for column in 0..3 {
                let delta = 1e-5;
                let mut plus_h = h;
                let mut minus_h = h;
                plus_h[column] += delta;
                minus_h[column] -= delta;
                let high_h = step(&material, &previous, plus_h, duration);
                let low_h = step(&material, &previous, minus_h, duration);
                for row in 0..3 {
                    near(
                        (high_h.state.flux_density_t()[row] - low_h.state.flux_density_t()[row])
                            / (2.0 * delta),
                        candidate.incremental_permeability_h_m[row][column],
                        2e-8,
                    );
                    assert_eq!(
                        candidate.incremental_permeability_h_m[row][column],
                        candidate.incremental_permeability_h_m[column][row]
                    );
                }
                let mut plus_b = b;
                let mut minus_b = b;
                plus_b[column] += delta;
                minus_b[column] -= delta;
                let high_b = material
                    .candidate_flux_density_timed(
                        &previous,
                        plus_b,
                        duration,
                        controls(),
                        &mut || false,
                    )
                    .unwrap();
                let low_b = material
                    .candidate_flux_density_timed(
                        &previous,
                        minus_b,
                        duration,
                        controls(),
                        &mut || false,
                    )
                    .unwrap();
                near(
                    (high_b.incremental_potential_j_m3().unwrap()
                        - low_b.incremental_potential_j_m3().unwrap())
                        / (2.0 * delta),
                    h[column],
                    3e-7,
                );
                for row in 0..3 {
                    near(
                        (high_b.state.field_a_m()[row] - low_b.state.field_a_m()[row])
                            / (2.0 * delta),
                        inverse[row][column],
                        3e-6,
                    );
                }
            }
        }
        assert_eq!(previous, saved);
    }
}

#[test]
fn mixed_rate_branches_preserve_vector_rotation_and_separate_heat_owners() {
    let material = VectorHysteresis::new(
        0.01,
        30.0,
        vec![branch(0.2, 2.0, 1.0), branch(0.1, 0.0, 0.5)],
    )
    .unwrap()
    .with_viscosities_m_s_h(vec![0.4, 0.0])
    .unwrap();
    let rotate = |h: Vector| [-h[1], h[0], h[2]];
    let mut original = material.initial_state();
    let mut rotated = material.initial_state();
    for (h, duration) in [
        ([9.0, 0.0, 0.0], 0.1),
        ([3.0, 8.0, 0.0], 0.4),
        ([-5.0, -9.0, 3.0], 0.7),
        ([0.0; 3], 0.2),
        ([0.0; 3], 0.5),
    ] {
        let a = step(&material, &original, h, duration);
        let b = step(&material, &rotated, rotate(h), duration);
        vector_near(
            b.state.flux_density_t(),
            rotate(a.state.flux_density_t()),
            3e-12,
        );
        near(a.pinning_heat_j_m3, b.pinning_heat_j_m3, 3e-12);
        near(a.viscous_heat_j_m3, b.viscous_heat_j_m3, 3e-12);
        let dz = sub(a.state.polarizations_t()[0], original.polarizations_t()[0]);
        near(a.viscous_heat_j_m3, 0.4 / duration * dot(dz, dz), 3e-12);
        audit(&a);
        audit(&b);
        original = a.state;
        rotated = b.state;
    }
}

#[test]
fn extreme_stiffness_keeps_physical_work_and_never_inflates_pinning_tolerance() {
    let extreme = law(2.0, 1e300);
    let candidate = step(&extreme, &extreme.initial_state(), [5.0, 0.0, 0.0], 1.0);
    // dz² underflows here, but dz.(eta*dz/dt) and the physical force do not.
    near(candidate.viscous_heat_j_m3 / 9e-300, 1.0, 3e-15);
    assert!(candidate.maximum_constitutive_residual_a_m < 1e-12);
    audit(&candidate);

    let material = law(2.0, 1.0);
    let previous = step(&material, &material.initial_state(), [5.0, 0.0, 0.0], 1e20).state;
    let saved = previous.clone();
    // The required slip is below representable resolution about this remanence.
    // A shifted enormous effective field must not relax the actual force gate.
    assert!(matches!(
        material
            .candidate_field_timed(&previous, [8.0, 0.0, 0.0], 1e-20, controls(), &mut || false),
        Err(VectorHysteresisError::Invalid(_))
    ));
    assert_eq!(previous, saved);
}

#[test]
fn duration_profile_and_cancellation_admission_preserve_accepted_history() {
    let material = law(2.0, 0.4);
    let previous = step(&material, &material.initial_state(), [5.0, 0.0, 0.0], 0.3).state;
    let saved = previous.clone();
    assert_eq!(material.viscosities_m_s_h(), Some([0.4].as_slice()));
    assert!(material.is_rate_dependent());
    for invalid in [vec![], vec![-1.0], vec![f64::NAN], vec![f64::INFINITY]] {
        assert!(material.clone().with_viscosities_m_s_h(invalid).is_err());
    }
    for duration in [0.0, -1.0, f64::NAN, f64::INFINITY, f64::from_bits(1)] {
        assert!(
            material
                .candidate_field_timed(
                    &previous,
                    [5.0, 0.0, 0.0],
                    duration,
                    controls(),
                    &mut || false
                )
                .is_err()
        );
    }
    assert!(
        material
            .candidate_field(&previous, [5.0, 0.0, 0.0], controls(), &mut || false)
            .is_err()
    );
    assert!(
        material
            .candidate_flux_density(
                &previous,
                previous.flux_density_t(),
                controls(),
                &mut || false
            )
            .is_err()
    );
    for foreign in [law(2.0, 0.0), law(2.0, 0.5)] {
        assert!(
            foreign
                .candidate_field_timed(&previous, [5.0, 0.0, 0.0], 0.3, controls(), &mut || false)
                .is_err()
        );
    }
    assert!(matches!(
        material.candidate_field_timed(&previous, [5.0, 0.0, 0.0], 0.3, controls(), &mut || true),
        Err(VectorHysteresisError::Cancelled)
    ));
    assert_eq!(previous, saved);
}
