//! Manufactured midpoint memory trajectories with independent ODE/work references.
use super::*;

fn dot<const N: usize>(a: &[f64; N], b: &[f64; N]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}

#[test]
fn g1_fixed_strain_maxwell_decay_retains_signed_tensor_memory_and_releases_heat() {
    let direction = [1.0, -1.0, 0.0, 2.0, 0.0, 0.0];
    let branch =
        core::array::from_fn(|i| core::array::from_fn(|j| 1e10 * direction[i] * direction[j]));
    let (card, mechanical, dielectric) = fixture(branch, [0.0, 0.0, 1.0], |_| {});
    let material = bound(&card, &mechanical, &dielectric);
    let zero = material.zero_memory().unwrap();
    let strain = [1e-5, -2e-5, 0.0, 4e-6, 0.0, 0.0];
    let stress = apply(&branch, &strain);
    let amplitude_squared = dot(&strain, &stress);
    let ramp_dt = 0.25 * MECHANICAL_TIME_S;
    let ramp_factor = 1.0 / (1.0 + ramp_dt / (2.0 * MECHANICAL_TIME_S));
    // The API consumes HALF the endpoint strain increment. This finite ramp
    // initializes the branch from physical strain, without guessing Gram axes.
    let (mut memory, ramp) = material
        .midpoint_candidate(&zero, ramp_dt, strain.map(|s| 0.5 * s), [0.0; 3])
        .unwrap()
        .into_parts();
    let initial_energy = 0.5 * amplitude_squared * ramp_factor * ramp_factor;
    near(ramp.memory_energy_after_j_m3, initial_energy, 1e-11, 1e-12);
    near(
        ramp.mechanical_work_j_m3,
        0.5 * amplitude_squared * ramp_factor,
        1e-11,
        1e-12,
    );
    near(
        ramp.mechanical_heat_j_m3,
        ramp_dt * amplitude_squared * ramp_factor * ramp_factor / (4.0 * MECHANICAL_TIME_S),
        1e-11,
        1e-12,
    );
    let dt = 0.4 * MECHANICAL_TIME_S;
    let decay = (2.0 * MECHANICAL_TIME_S - dt) / (2.0 * MECHANICAL_TIME_S + dt);
    for index in 1..=8 {
        let (next, step) = material
            .midpoint_candidate(&memory, dt, [0.0; 6], [0.0; 3])
            .unwrap()
            .into_parts();
        let before_factor = ramp_factor * decay.powi(index - 1);
        let after_factor = ramp_factor * decay.powi(index);
        let midpoint_factor = 0.5 * (before_factor + after_factor);
        let before_energy = 0.5 * amplitude_squared * before_factor * before_factor;
        let after_energy = 0.5 * amplitude_squared * after_factor * after_factor;
        let midpoint_energy = 0.5 * amplitude_squared * midpoint_factor * midpoint_factor;
        let expected_heat = 2.0 * dt * midpoint_energy / MECHANICAL_TIME_S;
        let recovered = material.memory_response(&next).unwrap();
        for i in 0..6 {
            near(
                step.midpoint.additional_stress_mandel_pa[i],
                stress[i] * midpoint_factor,
                1e-8,
                1e-12,
            );
            near(
                step.endpoint.additional_stress_mandel_pa[i],
                stress[i] * after_factor,
                1e-8,
                1e-12,
            );
            near(
                recovered.additional_stress_mandel_pa[i],
                stress[i] * after_factor,
                1e-8,
                1e-12,
            );
        }
        near(step.memory_energy_before_j_m3, before_energy, 1e-11, 1e-12);
        near(step.memory_energy_after_j_m3, after_energy, 1e-11, 1e-12);
        near(
            step.midpoint.mechanical_memory_energy_j_m3,
            midpoint_energy,
            1e-11,
            1e-12,
        );
        near(step.mechanical_heat_j_m3, expected_heat, 1e-11, 1e-12);
        assert!(step.mechanical_heat_j_m3 > 0.0);
        assert!(step.memory_energy_after_j_m3 < step.memory_energy_before_j_m3);
        assert_eq!(step.mechanical_work_j_m3, 0.0);
        assert_eq!(step.dielectric_work_j_m3, 0.0);
        assert_eq!(step.dielectric_heat_j_m3, 0.0);
        near(
            after_energy + step.mechanical_heat_j_m3,
            before_energy,
            1e-11,
            1e-12,
        );
        near(step.energy_residual_j_m3, 0.0, 1e-11, 0.0);
        memory = next;
    }
}

#[test]
fn g1_debye_field_charging_accounts_for_polarization_storage_heat_and_source_work() {
    let direction = [0.6, 0.8, 0.0];
    let (card, mechanical, dielectric) = fixture(diagonal_branch(), direction, |_| {});
    let material = bound(&card, &mechanical, &dielectric);
    let mut memory = material.zero_memory().unwrap();
    let field = 2000.0;
    let equilibrium = PERMITTIVITY_INCREMENT_F_M * field;
    let dt = 0.5 * DIELECTRIC_TIME_S;
    let decay = (2.0 * DIELECTRIC_TIME_S - dt) / (2.0 * DIELECTRIC_TIME_S + dt);
    let first_polarization = 0.5 * (1.0 - decay) * equilibrium;
    let mut previous_polarization = 0.0;
    let mut accumulated_work = 0.0;
    let mut accumulated_heat = 0.0;
    for index in 1..=8 {
        // First ramp E from zero; subsequent intervals hold the same field.
        // Only polarization memory belongs to this material subaccount.
        let midpoint_field = if index == 1 { 0.5 * field } else { field };
        let (next, step) = material
            .midpoint_candidate(&memory, dt, [0.0; 6], direction.map(|n| n * midpoint_field))
            .unwrap()
            .into_parts();
        let polarization = equilibrium + (first_polarization - equilibrium) * decay.powi(index - 1);
        let increment = polarization - previous_polarization;
        let energy = polarization * polarization / (2.0 * PERMITTIVITY_INCREMENT_F_M);
        let expected_heat =
            DIELECTRIC_TIME_S * increment * increment / (PERMITTIVITY_INCREMENT_F_M * dt);
        let expected_work = midpoint_field * increment;
        near(next.polarization_c_m2()[0], polarization, 1e-18, 1e-12);
        for i in 0..3 {
            near(
                step.endpoint.polarization_c_m2[i],
                direction[i] * polarization,
                1e-18,
                1e-12,
            );
            near(
                step.midpoint.polarization_c_m2[i],
                direction[i] * 0.5 * (polarization + previous_polarization),
                1e-18,
                1e-12,
            );
        }
        near(step.memory_energy_after_j_m3, energy, 1e-14, 1e-12);
        near(step.dielectric_heat_j_m3, expected_heat, 1e-14, 1e-12);
        near(step.dielectric_work_j_m3, expected_work, 1e-14, 1e-12);
        assert!(step.dielectric_heat_j_m3 > 0.0);
        assert!(polarization > previous_polarization && polarization < equilibrium);
        assert_eq!(step.mechanical_heat_j_m3, 0.0);
        assert_eq!(step.mechanical_work_j_m3, 0.0);
        assert_eq!(step.endpoint.mechanical_memory_energy_j_m3, 0.0);
        near(step.energy_residual_j_m3, 0.0, 1e-13, 0.0);
        accumulated_work += expected_work;
        accumulated_heat += step.dielectric_heat_j_m3;
        near(accumulated_work, energy + accumulated_heat, 1e-13, 1e-12);
        previous_polarization = polarization;
        memory = next;
    }
}

#[test]
fn g3_midpoint_memory_converges_quadratically_to_both_exact_relaxation_trajectories() {
    let branch = diagonal_branch();
    let (card, mechanical, dielectric) = fixture(branch, [0.0, 0.0, 1.0], |_| {});
    let material = bound(&card, &mechanical, &dielectric);
    let strain = [1e-5, 0.0, 0.0, 0.0, 0.0, 0.0];
    let field = [0.0, 0.0, 2000.0];
    let ramp_dt = 0.1 * MECHANICAL_TIME_S;
    let (initial, _) = material
        .midpoint_candidate(
            &material.zero_memory().unwrap(),
            ramp_dt,
            strain.map(|s| 0.5 * s),
            field.map(|e| 0.5 * e),
        )
        .unwrap()
        .into_parts();
    let initial_stress = branch[0][0] * strain[0] / (1.0 + ramp_dt / (2.0 * MECHANICAL_TIME_S));
    let equilibrium_polarization = PERMITTIVITY_INCREMENT_F_M * field[2];
    let initial_polarization =
        equilibrium_polarization * ramp_dt / (2.0 * DIELECTRIC_TIME_S + ramp_dt);
    let duration = 2.0 * MECHANICAL_TIME_S;
    let exact_stress = initial_stress * fs_math::det::exp(-duration / MECHANICAL_TIME_S);
    let exact_polarization = equilibrium_polarization
        + (initial_polarization - equilibrium_polarization)
            * fs_math::det::exp(-duration / DIELECTRIC_TIME_S);
    let mut previous_errors: Option<[f64; 2]> = None;
    for count in [8, 16, 32] {
        let mut memory = initial.clone();
        for _ in 0..count {
            let (next, step) = material
                .midpoint_candidate(&memory, duration / f64::from(count), [0.0; 6], field)
                .unwrap()
                .into_parts();
            assert!(step.mechanical_heat_j_m3 >= 0.0);
            assert!(step.dielectric_heat_j_m3 >= 0.0);
            near(step.energy_residual_j_m3, 0.0, 1e-11, 0.0);
            memory = next;
        }
        let response = material.memory_response(&memory).unwrap();
        let errors = [
            (response.additional_stress_mandel_pa[0] - exact_stress).abs() / initial_stress,
            (response.polarization_c_m2[2] - exact_polarization).abs()
                / (equilibrium_polarization - initial_polarization),
        ];
        if let Some(previous) = previous_errors {
            for (error, prior) in errors.into_iter().zip(previous) {
                assert!(
                    error > 0.23 * prior && error < 0.27 * prior,
                    "midpoint refinement at {count} steps: {error:e} after {prior:e}"
                );
            }
        }
        previous_errors = Some(errors);
    }
    let errors = previous_errors.unwrap();
    assert!(errors[0] < 1e-4 && errors[1] < 5e-6, "{errors:?}");
}

#[test]
fn g0_refused_midpoint_candidates_preserve_memory_and_require_the_original_broadband_model() {
    let (card, mechanical, dielectric) = fixture(diagonal_branch(), [0.0, 0.0, 1.0], |_| {});
    let material = bound(&card, &mechanical, &dielectric);
    let dt = 0.25 * MECHANICAL_TIME_S;
    let (memory, _) = material
        .midpoint_candidate(
            &material.zero_memory().unwrap(),
            dt,
            [5e-6, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 1000.0],
        )
        .unwrap()
        .into_parts();
    let snapshot = memory.clone();
    let increment = [-2e-6, 1e-6, 0.0, 0.0, 0.0, 0.0];
    let field = [300.0, 200.0, 1500.0];
    let expected = material
        .midpoint_candidate(&memory, dt, increment, field)
        .unwrap();
    for (bad_dt, bad_increment, bad_field) in [
        (0.0, increment, field),
        (-1.0, increment, field),
        (f64::NAN, increment, field),
        (f64::INFINITY, increment, field),
        (dt, [f64::NAN; 6], field),
        (dt, increment, [f64::INFINITY; 3]),
        (dt, [f64::MAX; 6], field),
    ] {
        assert!(
            material
                .midpoint_candidate(&memory, bad_dt, bad_increment, bad_field)
                .is_err()
        );
        assert_eq!(memory, snapshot);
        let retry = material
            .midpoint_candidate(&memory, dt, increment, field)
            .unwrap();
        assert_eq!(retry.memory(), expected.memory());
    }
    let (other_card, m, d) = fixture(diagonal_branch(), [0.0, 0.0, 1.0], |claim| {
        if claim.key.name() == "piezo_memory_tau" {
            claim.value = PropertyValue::Scalar {
                value: 0.003,
                dims: Time::DIMS,
            };
        }
    });
    let other = bound(&other_card, &m, &d);
    let rotated = material
        .rotated(ContentHash([31; 32]), rotation([0.0, 0.0, 1.0], 0.4))
        .unwrap();
    let dielectric_only = RelaxingTensorPiezoelectric::resolve(
        base_law(&card),
        &card,
        vec![],
        vec![dielectric],
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    for incompatible in [&other, &rotated, &dielectric_only] {
        assert!(incompatible.memory_response(&memory).is_err());
        assert!(
            incompatible
                .midpoint_candidate(&memory, dt, increment, field)
                .is_err()
        );
        assert_eq!(memory, snapshot);
    }
    let (limited_card, m, d) = fixture(diagonal_branch(), [0.0, 0.0, 1.0], |claim| {
        if claim.key.name() == "piezo_memory_tau" {
            claim.validity = claim.validity.clone().with("omega", 0.0, 1000.0);
        }
    });
    let limited = bound(&limited_card, &m, &d);
    assert_eq!(limited.frequency_band_rad_s(), [0.0, 1000.0]);
    assert!(limited.frequency_coefficients(100.0).is_ok());
    assert!(limited.zero_memory().is_err());
    assert!(
        limited
            .midpoint_candidate(&memory, dt, increment, field)
            .is_err()
    );
    assert_eq!(memory, snapshot);
}

#[test]
fn g1_g3_equilibrium_memory_and_static_response_share_the_actual_native_inventory() {
    let direction = [0.6, 0.8, 0.0];
    let field = [1200.0, -700.0, 300.0];
    let (card, mechanical, dielectric) = fixture(diagonal_branch(), direction, |_| {});
    for count in [0, 1, 2] {
        let material = RelaxingTensorPiezoelectric::resolve(
            base_law(&card),
            &card,
            vec![mechanical.clone()],
            vec![dielectric.clone(); count],
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
        let memory = material.equilibrium_memory(field).unwrap();
        assert_eq!(memory.model_identity(), material.identity());
        assert_eq!(memory.internal_variables(), material.internal_variables());
        assert_eq!(memory.internal_variables(), 6 + count);
        assert_eq!(memory.mechanical_memory_sqrt_pa(), &[[0.0; 6]]);
        for &polarization in memory.polarization_c_m2() {
            near(
                polarization,
                PERMITTIVITY_INCREMENT_F_M * 160.0,
                1e-20,
                3e-13,
            );
        }
        // This is an exact inventory transfer, including the represented pole
        // energy, not agreement of independently rounded analytic formulas.
        let response = material.memory_response(&memory).unwrap();
        assert_eq!(response, material.equilibrium_response(field).unwrap());
        assert_eq!(response.additional_stress_mandel_pa, [0.0; 6]);
        assert_eq!(response.mechanical_memory_energy_j_m3, 0.0);
        near(
            response.memory_energy_j_m3,
            0.5 * count as f64 * PERMITTIVITY_INCREMENT_F_M * 160.0 * 160.0,
            1e-20,
            3e-13,
        );
        assert_eq!(
            material.equilibrium_memory([0.0; 3]).unwrap(),
            material.zero_memory().unwrap()
        );

        let q = rotation([1.0 / 3.0_f64.sqrt(); 3], 0.73);
        let moved = material.rotated(ContentHash([84; 32]), q).unwrap();
        let moved_memory = moved.equilibrium_memory(apply(&q, &field)).unwrap();
        let moved_response = moved.memory_response(&moved_memory).unwrap();
        for (actual, expected) in moved_response
            .polarization_c_m2
            .iter()
            .zip(apply(&q, &response.polarization_c_m2))
        {
            near(*actual, expected, 1e-20, 3e-12);
        }
        near(
            moved_response.memory_energy_j_m3,
            response.memory_energy_j_m3,
            1e-20,
            3e-12,
        );
        assert!(material.memory_response(&moved_memory).is_err());
        assert!(moved.memory_response(&memory).is_err());
    }
}

#[test]
fn g1_relaxed_polarization_survives_a_hold_then_discharges_with_physical_heat() {
    let (card, mechanical, dielectric) = fixture(diagonal_branch(), [0.0, 0.0, 1.0], |_| {});
    let material = bound(&card, &mechanical, &dielectric);
    let field = [0.0, 0.0, 2000.0];
    let initial = material.equilibrium_memory(field).unwrap();
    let polarization = PERMITTIVITY_INCREMENT_F_M * field[2];
    let energy = 0.5 * polarization * field[2];
    for ratio in [1e-6, 0.4, 4.0] {
        let dt = ratio * DIELECTRIC_TIME_S;
        let mut memory = initial.clone();
        for _ in 0..8 {
            let history = material.history_correction(&memory, dt, field).unwrap();
            assert_eq!(history.mechanical_rhs_stress_mandel_pa, [0.0; 6]);
            for increment in history.polarization_history_increment_c_m2 {
                near(increment, 0.0, 3e-13 * polarization, 0.0);
            }
            let (next, step) = material
                .midpoint_candidate(&memory, dt, [0.0; 6], field)
                .unwrap()
                .into_parts();
            near(next.polarization_c_m2()[0], polarization, 0.0, 3e-13);
            near(step.memory_energy_after_j_m3, energy, 0.0, 3e-13);
            assert_eq!(step.mechanical_heat_j_m3, 0.0);
            // The original kernel reports its actual represented-state
            // friction, including roundoff; equilibrium does not zero heat.
            assert!(step.dielectric_heat_j_m3 >= 0.0);
            assert!(step.dielectric_heat_j_m3 <= 1e-24 * energy);
            near(step.energy_residual_j_m3, 0.0, 3e-13 * energy, 0.0);
            memory = next;
        }
    }
    let dt = 0.4 * DIELECTRIC_TIME_S;
    let decay = (2.0 * DIELECTRIC_TIME_S - dt) / (2.0 * DIELECTRIC_TIME_S + dt);
    let (next, step) = material
        .midpoint_candidate(&initial, dt, [0.0; 6], [0.0; 3])
        .unwrap()
        .into_parts();
    near(
        next.polarization_c_m2()[0],
        decay * polarization,
        0.0,
        3e-13,
    );
    near(
        step.memory_energy_after_j_m3,
        decay * decay * energy,
        0.0,
        3e-13,
    );
    near(
        step.dielectric_heat_j_m3,
        (1.0 - decay * decay) * energy,
        0.0,
        3e-13,
    );
    assert_eq!(step.dielectric_work_j_m3, 0.0);
    near(step.energy_residual_j_m3, 0.0, 3e-13 * energy, 0.0);
}

#[test]
fn g0_g1_equilibrium_memory_requires_broadband_sources_and_representable_storage() {
    let (card, mechanical, dielectric) = fixture(diagonal_branch(), [0.0, 0.0, 1.0], |_| {});
    let material = bound(&card, &mechanical, &dielectric);
    // E² overflows here, but the physical pole state and its storage fit.
    let large_field = [0.0, 0.0, 1e155];
    let memory = material.equilibrium_memory(large_field).unwrap();
    let response = material.memory_response(&memory).unwrap();
    assert_eq!(
        response,
        material.equilibrium_response(large_field).unwrap()
    );
    near(
        response.memory_energy_j_m3,
        (0.5 * PERMITTIVITY_INCREMENT_F_M * large_field[2]) * large_field[2],
        0.0,
        3e-13,
    );
    for field in [
        [f64::NAN, 0.0, 0.0],
        [0.0, f64::INFINITY, 0.0],
        [0.0, 0.0, 1e308],
        [0.0, 0.0, 1e-200],
    ] {
        assert!(material.equilibrium_memory(field).is_err());
    }
    for lower in [0.0, 10.0] {
        let (card, mechanical, dielectric) = fixture(diagonal_branch(), [0.0, 0.0, 1.0], |claim| {
            if claim.key.name() == "piezo_polarization_tau" {
                claim.validity = claim.validity.clone().with("omega", lower, 1000.0);
            }
        });
        let limited = bound(&card, &mechanical, &dielectric);
        assert_eq!(limited.frequency_band_rad_s(), [lower, 1000.0]);
        assert_eq!(
            limited.equilibrium_response([0.0, 0.0, 2000.0]).is_ok(),
            lower == 0.0
        );
        assert!(limited.equilibrium_memory([0.0, 0.0, 2000.0]).is_err());
    }
}
