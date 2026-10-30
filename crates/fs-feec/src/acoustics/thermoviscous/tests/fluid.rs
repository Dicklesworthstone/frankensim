use super::*;

// Manufactured single-phase fluid; these numbers are not a measured water card.
fn liquid() -> (ThermoviscousMedium, ThermoviscousThermodynamics) {
    let medium = ThermoviscousMedium {
        density_kg_m3: 1000.0,
        pressure_pa: 101325.0,
        temperature_k: 300.0,
        cp_j_kg_k: 4000.0,
        dynamic_viscosity_pa_s: 0.001,
        thermal_conductivity_w_m_k: 0.6,
        bulk_viscosity_pa_s: 0.002,
    };
    let alpha = 2.0e-4;
    let beta_s = 1.0 / (medium.density_kg_m3 * 1400.0_f64.powi(2));
    let beta_t =
        beta_s + alpha * alpha * medium.temperature_k / (medium.density_kg_m3 * medium.cp_j_kg_k);
    (
        medium,
        ThermoviscousThermodynamics {
            isothermal_compressibility_pa_inv: beta_t,
            isobaric_expansion_per_k: alpha,
        },
    )
}

#[test]
fn g1_general_fluid_recovers_original_ideal_gas_operator_and_fields() {
    let size = [0.004, 0.0008, 0.001];
    let (mesh, pos) = box_mesh([1, 2, 2], size);
    let boundary = boundaries(&mesh, &pos, size, true);
    let gas = air();
    let old = ThermoviscousFe::new(&mesh, &pos, gas, &boundary, limits(), &mut || false).unwrap();
    let general = ThermoviscousFe::new_with_thermodynamics(
        &mesh,
        &pos,
        gas,
        ThermoviscousThermodynamics::ideal_gas(gas).unwrap(),
        &boundary,
        limits(),
        &mut || false,
    )
    .unwrap();
    assert_eq!(old.constant, general.constant);
    for (&a, &b) in old.mass.iter().zip(&general.mass) {
        assert!((a - b).abs() <= 8.0 * f64::EPSILON * a.abs().max(b.abs()));
    }
    for frequency in [80.0, 1700.0] {
        let cp = C64::new(970.0, -45.0);
        let drive = [C64::new(1.0, 0.2), C64::new(-0.3, 0.1)];
        let a = old
            .prepare_with_heat_capacity(frequency, cp, &mut || false)
            .unwrap()
            .response(&drive, &mut || false)
            .unwrap();
        let b = general
            .prepare_with_heat_capacity(frequency, cp, &mut || false)
            .unwrap()
            .response(&drive, &mut || false)
            .unwrap();
        assert_work(&b);
        for (x, y) in a.inward_flow_m3_s.iter().zip(&b.inward_flow_m3_s) {
            assert!((*x - *y).abs() < 2e-10 * x.abs());
        }
    }
}

#[test]
fn g1_liquid_storage_has_the_actual_adiabatic_compressibility_and_wave_speed() {
    let size = [0.004, 0.001, 0.0012];
    let volume = size.iter().product::<f64>();
    let (mesh, pos) = box_mesh([1, 2, 2], size);
    let boundary = boundaries(&mesh, &pos, size, false);
    let (medium, thermo) = liquid();
    let model = ThermoviscousFe::new_with_thermodynamics(
        &mesh,
        &pos,
        medium,
        thermo,
        &boundary,
        limits(),
        &mut || false,
    )
    .unwrap();
    let expected_beta = 1.0 / (medium.density_kg_m3 * 1400.0_f64.powi(2));
    assert!((model.fe.medium().sound_speed_m_s / 1400.0 - 1.0).abs() < 2e-15);
    // Independent adiabatic uniform-volume solution: T'/T0 = alpha*p/(rho Cp).
    // A constant hierarchical P2 field has unit vertex coefficients, zero edges.
    let theta = thermo.isobaric_expansion_per_k / (medium.density_kg_m3 * medium.cp_j_kg_k);
    let mut state = vec![0.0; model.n];
    state[model.pressure_start..model.temperature_start].fill(1.0);
    for vertex in 0..mesh.vertex_count {
        state[model.temperature_start + model.temperature_map[vertex].unwrap()] = theta;
    }
    let mass_action = (0..model.n)
        .map(|i| {
            (0..model.n)
                .map(|j| model.mass[i * model.n + j] * state[j])
                .sum::<f64>()
        })
        .collect::<Vec<_>>();
    let compliance: f64 = mass_action[model.pressure_start..model.temperature_start]
        .iter()
        .sum();
    assert!((compliance / (volume * expected_beta) - 1.0).abs() < 2e-13);
    for i in model.temperature_start..model.n {
        let scale = (0..model.n)
            .map(|j| (model.mass[i * model.n + j] * state[j]).abs())
            .sum::<f64>();
        assert!(mass_action[i].abs() <= 64.0 * f64::EPSILON * scale);
    }
    for i in 0..model.n {
        for j in 0..model.n {
            let (a, b) = (model.mass[i * model.n + j], model.mass[j * model.n + i]);
            assert!((a - b).abs() <= 16.0 * f64::EPSILON * a.abs().max(b.abs()));
        }
    }
}

#[test]
fn g3_liquid_traction_work_and_temperature_sign_follow_the_material_eos() {
    let size = [0.004, 0.001, 0.0012];
    let (mesh, pos) = box_mesh([1, 2, 2], size);
    let boundary = boundaries(&mesh, &pos, size, true);
    let (medium, thermo) = liquid();
    let build = |m, t| {
        ThermoviscousFe::new_with_thermodynamics(
            &mesh,
            &pos,
            m,
            t,
            &boundary,
            limits(),
            &mut || false,
        )
        .unwrap()
    };
    let model = build(medium, thermo);
    let changed_reference = build(
        ThermoviscousMedium {
            pressure_pa: 8e7,
            ..medium
        },
        thermo,
    );
    let opposite = build(
        medium,
        ThermoviscousThermodynamics {
            isobaric_expansion_per_k: -thermo.isobaric_expansion_per_k,
            ..thermo
        },
    );
    let drive = [C64::new(1.0, 0.2), C64::new(-0.5, 0.1)];
    let solve = |model: &ThermoviscousFe<'_>| {
        model
            .prepare(400.0, &mut || false)
            .unwrap()
            .response(&drive, &mut || false)
            .unwrap()
    };
    let a = solve(&model);
    let b = solve(&changed_reference);
    let c = solve(&opposite);
    assert_work(&a);
    assert_work(&c);
    // The declared pressure is not converted to a fictitious liquid gas constant.
    assert_eq!(a.inward_flow_m3_s, b.inward_flow_m3_s);
    assert_eq!(a.temperature_k, b.temperature_k);
    for (x, y) in a.inward_flow_m3_s.iter().zip(&c.inward_flow_m3_s) {
        assert!((*x - *y).abs() <= 1e-11 * x.abs().max(1e-30));
    }
    let maximum = a
        .temperature_k
        .iter()
        .map(|z| z.abs())
        .fold(0.0_f64, f64::max);
    assert!(maximum > 0.0);
    for (x, y) in a.temperature_k.iter().zip(&c.temperature_k) {
        assert!((*x + *y).abs() <= 1e-11 * maximum);
    }
    let y = model
        .admittance_with_heat_capacity(400.0, C64::from_re(medium.cp_j_kg_k), &mut || false)
        .unwrap();
    assert!((y[1] - y[2]).abs() < 1e-9 * y[1].abs());
}

#[test]
fn g1_liquid_viscous_loss_converges_to_independent_rectangular_poiseuille() {
    let size = [0.05, 0.0005, 0.001];
    let [length, height, width] = size;
    let (medium, thermo) = liquid();
    let sum = (1..=501)
        .step_by(2)
        .map(|m| {
            let m = m as f64;
            (m * std::f64::consts::PI * width / (2.0 * height)).tanh() / m.powi(5)
        })
        .sum::<f64>();
    let expected = width * height.powi(3) / (12.0 * medium.dynamic_viscosity_pa_s * length)
        * (1.0 - 192.0 * height / (std::f64::consts::PI.powi(5) * width) * sum);
    let mut errors = Vec::new();
    for across in [2, 3, 4] {
        let (mesh, pos) = box_mesh([1, across, across], size);
        let boundary = boundaries(&mesh, &pos, size, true);
        let model = ThermoviscousFe::new_with_thermodynamics(
            &mesh,
            &pos,
            medium,
            thermo,
            &boundary,
            limits(),
            &mut || false,
        )
        .unwrap();
        let field = model
            .prepare(1e-5, &mut || false)
            .unwrap()
            .response(&[C64::from_re(0.5), C64::from_re(-0.5)], &mut || false)
            .unwrap();
        assert_work(&field);
        errors.push((field.inward_flow_m3_s[0].re / expected - 1.0).abs());
    }
    assert!(errors[1] < errors[0] && errors[2] < errors[1], "{errors:?}");
    assert!(errors[2] < 0.015, "{errors:?}");
}

#[test]
fn general_fluid_refuses_unstable_storage_and_active_caloric_response() {
    let size = [0.002, 0.001, 0.001];
    let (mesh, pos) = box_mesh([1, 2, 2], size);
    let boundary = boundaries(&mesh, &pos, size, true);
    let (medium, thermo) = liquid();
    for alpha in [f64::NAN, 1.0] {
        assert!(
            ThermoviscousFe::new_with_thermodynamics(
                &mesh,
                &pos,
                medium,
                ThermoviscousThermodynamics {
                    isobaric_expansion_per_k: alpha,
                    ..thermo
                },
                &boundary,
                limits(),
                &mut || false,
            )
            .is_err()
        );
    }
    let model = ThermoviscousFe::new_with_thermodynamics(
        &mesh,
        &pos,
        medium,
        thermo,
        &boundary,
        limits(),
        &mut || false,
    )
    .unwrap();
    let threshold = thermo.isobaric_expansion_per_k.powi(2) * medium.temperature_k
        / (medium.density_kg_m3 * thermo.isothermal_compressibility_pa_inv);
    assert!(
        model
            .prepare_with_heat_capacity(100.0, C64::from_re(0.5 * threshold), &mut || false)
            .is_err()
    );
    assert!(
        model
            .prepare_with_heat_capacity(100.0, C64::new(medium.cp_j_kg_k, 1.0), &mut || false)
            .is_err()
    );
}
