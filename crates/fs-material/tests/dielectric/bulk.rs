//! G0/G1/G3 independent bulk polarization, energy and time-domain oracles.
use super::*;
use fs_material::field_storage::dielectric::bulk::ResolvedDebyeMaterial;

fn source() -> MaterialCard {
    card(
        &[
            ("epsilon_inf", ELECTRIC_PERMITTIVITY_DIMS, 2e-6),
            ("delta1", ELECTRIC_PERMITTIVITY_DIMS, 3e-6),
            ("tau1", Time::DIMS, 0.1),
            ("delta2", ELECTRIC_PERMITTIVITY_DIMS, 5e-6),
            ("tau2", Time::DIMS, 0.7),
        ],
        ValidityDomain::unconstrained(),
    )
}

fn law(poles: &[DebyePoleProperties]) -> ResolvedDebyeMaterial {
    ResolvedDebyeMaterial::resolve(
        &source(),
        &QueryPoint::new(),
        MaterialPropertySelection::SingleClaimOnly,
        "epsilon_inf",
        poles,
    )
    .unwrap()
}

fn close(a: f64, b: f64, scale: f64) {
    assert!((a - b).abs() <= 5e-13 * scale, "{a} != {b}, scale={scale}");
}

#[test]
fn g1_bulk_poles_preserve_sources_and_uniform_circuit_storage() {
    let source = source();
    let names = properties();
    let capacitor = resolve(&source, &names);
    let material = ResolvedDebyeMaterial::try_from_material(
        capacitor.material().clone(),
        "epsilon_inf",
        &names,
    )
    .unwrap();
    assert_eq!(material.material(), capacitor.material());
    assert_eq!(material.high_frequency_property(), "epsilon_inf");
    assert_eq!(material.pole_properties(), names);
    near(material.high_frequency_permittivity_f_m(), 2e-6);
    near(material.equilibrium_permittivity_f_m(), 1e-5);
    for (bulk, lumped) in material.poles().iter().zip(capacitor.poles()) {
        assert_eq!(
            bulk.permittivity_increment_f_m,
            lumped.permittivity_increment_f_m
        );
        assert_eq!(bulk.relaxation_time_s, lumped.relaxation_time_s);
    }
    let field = 500.0;
    let polarization = [[2e-3], [-1e-3]];
    let volume = capacitor.area_m2() * capacitor.gap_m();
    let bulk_energy = material
        .stored_energy_density_j_m3([field], &polarization)
        .unwrap();
    let charges: Vec<_> = polarization
        .iter()
        .map(|p| p[0] * capacitor.area_m2())
        .collect();
    let uniform_energy = capacitor
        .stored_energy_j(
            material.high_frequency_permittivity_f_m() * field * capacitor.area_m2(),
            &charges,
        )
        .unwrap();
    near(bulk_energy * volume, uniform_energy);
    near(
        material.displacement_c_m2([field], &polarization).unwrap()[0],
        2e-3,
    );
}

#[test]
fn g1_midpoint_elimination_matches_local_ode_and_closes_vector_work() {
    let material = law(&properties());
    let old_field = [200.0, -100.0, 40.0];
    let new_field = [-80.0, 300.0, 90.0];
    let old_polarization = [[1e-4, -2e-4, 3e-4], [-4e-4, 5e-4, -1e-4]];
    let dt = 0.017;
    let effective = material.midpoint_permittivity_f_m(dt).unwrap();
    let history = material
        .midpoint_history_c_m2(dt, old_field, &old_polarization)
        .unwrap();
    let step = material
        .step_midpoint(dt, old_field, new_field, &old_polarization)
        .unwrap();
    for (k, pole) in material.poles().iter().enumerate() {
        for i in 0..3 {
            let p0 = old_polarization[k][i];
            let p1 = step.polarization_after_c_m2[k][i];
            let lhs = pole.relaxation_time_s * (p1 - p0) / dt + 0.5 * (p0 + p1);
            let rhs = pole.permittivity_increment_f_m * 0.5 * (old_field[i] + new_field[i]);
            close(lhs, rhs, 1e-3);
        }
    }
    for i in 0..3 {
        close(
            step.displacement_after_c_m2[i],
            effective * new_field[i] + history[i],
            1e-3,
        );
    }
    let expected_heat: f64 = material
        .poles()
        .iter()
        .enumerate()
        .map(|(k, pole)| {
            (0..3)
                .map(|i| {
                    let delta = step.polarization_after_c_m2[k][i] - old_polarization[k][i];
                    pole.relaxation_time_s * delta * delta / (dt * pole.permittivity_increment_f_m)
                })
                .sum::<f64>()
        })
        .sum();
    near(step.heat_j_m3, expected_heat);
    assert!(step.heat_j_m3 > 0.0);
    close(
        step.stored_change_j_m3,
        step.stored_after_j_m3 - step.stored_before_j_m3,
        1.0,
    );
    close(
        step.stored_change_j_m3 + step.heat_j_m3,
        step.supplied_work_j_m3,
        1.0,
    );
    close(step.energy_residual_j_m3, 0.0, 1.0);
    assert_eq!(
        step,
        material
            .step_midpoint(dt, old_field, new_field, &old_polarization)
            .unwrap()
    );
}

#[test]
fn g1_polarization_decay_is_second_order_and_converts_storage_to_heat() {
    let material = law(&properties()[..1]);
    let total_time = 0.1;
    let initial = 3e-4;
    let mut errors = Vec::new();
    for count in [8, 16, 32] {
        let dt = total_time / f64::from(count);
        let mut p = vec![[initial]];
        let stored = material.stored_energy_density_j_m3([0.0], &p).unwrap();
        let mut heat = 0.0;
        for _ in 0..count {
            let step = material.step_midpoint(dt, [0.0], [0.0], &p).unwrap();
            heat += step.heat_j_m3;
            assert_eq!(step.supplied_work_j_m3, 0.0);
            assert!(step.heat_j_m3 > 0.0);
            p = step.polarization_after_c_m2;
        }
        let final_storage = material.stored_energy_density_j_m3([0.0], &p).unwrap();
        near(final_storage + heat, stored);
        errors.push((p[0][0] - initial * (-1.0_f64).exp()).abs());
    }
    assert!(
        errors[0] / errors[1] > 3.9 && errors[1] / errors[2] > 3.9,
        "{errors:?}"
    );
}

#[test]
fn g1_time_domain_poles_match_discrete_harmonic_debye_response() {
    let material = law(&properties());
    let dt = 0.004;
    let omega = 17.0;
    // Midpoint maps the sampled harmonic derivative to i*(2/dt)*tan(omega*dt/2).
    let discrete_omega = (2.0 / dt) * (0.5_f64 * omega * dt).tan();
    let field = |phase: f64| [phase.cos(), phase.sin()];
    let exact = |phase: f64| -> Vec<[f64; 2]> {
        material
            .poles()
            .iter()
            .map(|pole| {
                let x = discrete_omega * pole.relaxation_time_s;
                let scale = pole.permittivity_increment_f_m / (1.0 + x * x);
                [
                    scale * (phase.cos() + x * phase.sin()),
                    scale * (phase.sin() - x * phase.cos()),
                ]
            })
            .collect()
    };
    let mut p = exact(0.0);
    for index in 0..80 {
        let start = f64::from(index) * dt * omega;
        let end = f64::from(index + 1) * dt * omega;
        let step = material
            .step_midpoint(dt, field(start), field(end), &p)
            .unwrap();
        for (actual, expected) in step.polarization_after_c_m2.iter().zip(exact(end)) {
            for i in 0..2 {
                close(actual[i], expected[i], 1e-5);
            }
        }
        close(step.energy_residual_j_m3, 0.0, 1e-5);
        assert!(step.heat_j_m3 > 0.0);
        p = step.polarization_after_c_m2;
    }
}

#[test]
fn g3_isotropic_rotation_and_zero_pole_limits_preserve_storage_and_work() {
    let material = law(&properties());
    let rotate = |v: [f64; 3]| [-v[1], v[0], v[2]];
    let e0 = [2.0, -3.0, 5.0];
    let e1 = [-7.0, 11.0, 13.0];
    let p = [[2e-6, -3e-6, 1e-6], [-5e-6, 7e-6, 3e-6]];
    let a = material.step_midpoint(0.03, e0, e1, &p).unwrap();
    let b = material
        .step_midpoint(0.03, rotate(e0), rotate(e1), &p.map(rotate))
        .unwrap();
    near(a.heat_j_m3, b.heat_j_m3);
    near(a.stored_after_j_m3, b.stored_after_j_m3);
    near(a.supplied_work_j_m3, b.supplied_work_j_m3);
    for (left, right) in a
        .polarization_after_c_m2
        .iter()
        .zip(&b.polarization_after_c_m2)
    {
        assert_eq!(rotate(*left), *right);
    }
    let empty = law(&[]);
    assert_eq!(empty.midpoint_permittivity_f_m(0.3).unwrap(), 2e-6);
    assert_eq!(empty.midpoint_history_c_m2(0.3, e0, &[]).unwrap(), [0.0; 3]);
    let result = empty.step_midpoint(0.3, e0, e1, &[]).unwrap();
    assert!(result.polarization_after_c_m2.is_empty());
    assert_eq!(result.heat_j_m3, 0.0);
    assert_eq!(result.displacement_after_c_m2, e1.map(|e| 2e-6 * e));
    assert_eq!(result.stored_change_j_m3, result.supplied_work_j_m3);
    assert_eq!(result.energy_residual_j_m3, 0.0);

    // An unchanged large bias must not erase the independently resolved work.
    let biased = empty.step_midpoint(0.3, [1e16], [1e16 + 2.0], &[]).unwrap();
    near(biased.stored_change_j_m3, 4e10);
    assert_eq!(biased.energy_residual_j_m3, 0.0);
    assert!(
        (biased.stored_after_j_m3 - biased.stored_before_j_m3 - biased.stored_change_j_m3).abs()
            > 1e6
    );
}

#[test]
fn g0_invalid_sources_histories_and_intervals_refuse_without_mutation() {
    let material = law(&properties());
    let history = [[1e-4, 2e-4], [-3e-4, 4e-4]];
    for dt in [-1.0, 0.0, f64::NAN, f64::INFINITY] {
        assert!(material.midpoint_permittivity_f_m(dt).is_err());
        assert!(
            material
                .midpoint_history_c_m2(dt, [1.0, 2.0], &history)
                .is_err()
        );
        assert!(
            material
                .step_midpoint(dt, [1.0, 2.0], [3.0, 4.0], &history)
                .is_err()
        );
    }
    assert!(
        material
            .step_midpoint(0.1, [1.0, 2.0], [3.0, 4.0], &history[..1])
            .is_err()
    );
    assert!(
        material
            .step_midpoint(0.1, [1.0, 2.0], [f64::NAN, 4.0], &history)
            .is_err()
    );
    assert!(
        material
            .stored_energy_density_j_m3([f64::INFINITY, 0.0], &history)
            .is_err()
    );
    assert_eq!(history, [[1e-4, 2e-4], [-3e-4, 4e-4]]);
    let bad = card(
        &[
            ("epsilon_inf", ELECTRIC_PERMITTIVITY_DIMS, 2e-6),
            ("delta1", ELECTRIC_PERMITTIVITY_DIMS, -3e-6),
            ("tau1", Time::DIMS, 0.1),
        ],
        ValidityDomain::unconstrained(),
    );
    assert!(
        ResolvedDebyeMaterial::resolve(
            &bad,
            &QueryPoint::new(),
            MaterialPropertySelection::SingleClaimOnly,
            "epsilon_inf",
            &properties()[..1]
        )
        .is_err()
    );
    let duplicate = [properties()[0].clone(), properties()[0].clone()];
    assert!(
        ResolvedDebyeMaterial::resolve(
            &source(),
            &QueryPoint::new(),
            MaterialPropertySelection::SingleClaimOnly,
            "epsilon_inf",
            &duplicate
        )
        .is_err()
    );
}
