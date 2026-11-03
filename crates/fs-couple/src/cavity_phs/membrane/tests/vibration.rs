//! Coupled modal regressions use the original gas, membrane and eigen owners.
use super::*;
use fs_solid::membrane::vibration::{VibrationError, VibrationSettings};

fn modal_controls() -> VibrationSettings {
    VibrationSettings {
        max_free_dofs: 16,
        relative_residual_tolerance: 1e-8,
    }
}

#[test]
fn zero_jump_modes_include_native_gas_stiffness_and_reciprocal_pressure() {
    let (patch, reference, fixed) = pyramid();
    let prescribed: Vec<_> = fixed
        .iter()
        .enumerate()
        .filter_map(|(i, &value)| value.map(|value| (i, value)))
        .collect();
    let loads = [[0.0; 3]; 5];
    let held = patch
        .solve_vibration_with_pressure(
            &loads,
            0.0,
            &prescribed,
            &reference,
            controls().membrane,
            modal_controls(),
            || false,
        )
        .unwrap();
    let structural = 2.0 * STIFFNESS * (STRETCH.powi(2) - 1.0);
    let mass = patch.nodal_mass_kg()[0];
    near(held.modes()[0].lambda, structural / mass, 1e-8, 2e-12);

    // Analytic Cp(T), not a fitted mode or a fixed-gamma production substitute.
    for slope in [0.0, 0.002] {
        let source = source_gas(0.75, slope);
        let model = GasMembrane::new(
            &patch,
            source.clone(),
            &reference,
            &fixed,
            P0,
            [0.0007, 0.0015],
            || false,
        )
        .unwrap();
        let vibration = model
            .solve_vibration(&loads, &reference, controls(), modal_controls(), || false)
            .unwrap();
        let basis = vibration.basis();
        let gamma_at_400 = 1.0 + 1.0 / (2.25 + 0.75 * slope * 400.0);
        let gas_stiffness = gamma_at_400 * P0 / V0 * volume_gradient().powi(2);
        assert_eq!(basis.modes().len(), 1);
        near(
            basis.modes()[0].lambda,
            (structural + gas_stiffness) / mass,
            1e-7,
            2e-12,
        );
        assert!(basis.modes()[0].lambda > 2.0 * held.modes()[0].lambda);
        let b = basis.volume_coupling()[0];
        near(b * b, volume_gradient().powi(2) / mass, 1e-14, 2e-12);
        let q = [b.signum() * 1e-6];
        let volume = vibration.volume_perturbation_m3(&q).unwrap();
        assert!(volume > 0.0);
        near(volume, b * q[0], 1e-20, 1e-14);
        near(
            vibration.pressure_perturbation_pa(&q).unwrap(),
            -gamma_at_400 * P0 / V0 * volume,
            1e-10,
            2e-12,
        );
        let displacement = basis.reconstruct(&q).unwrap();
        let volume_from_nodes: f64 = basis
            .pressure_load_per_pa()
            .iter()
            .zip(&displacement)
            .flat_map(|(force, displacement)| force.iter().zip(displacement))
            .map(|(force, displacement)| force * displacement)
            .sum();
        near(volume, volume_from_nodes, 1e-20, 2e-12);
        let external_pressure = 123.0;
        let force: Vec<_> = basis
            .pressure_load_per_pa()
            .iter()
            .map(|value| value.map(|x| external_pressure * x))
            .collect();
        let modal_force = basis.project_nodal_force(&force).unwrap();
        near(
            modal_force[0] * q[0],
            external_pressure * volume,
            1e-14,
            2e-12,
        );
        assert_eq!(
            vibration.solution().equilibrium.gas.mass_kg(),
            source.mass_kg()
        );
        assert_eq!(model.reference_gas(), &MembraneGas::Ideal(source));
    }
}

#[test]
fn deformed_six_coordinate_modes_use_full_geometric_and_gas_tangent() {
    let (patch, reference, fixed) =
        chart(Mesh2::triangles(0.3, 0.2, 3, 2), &[15, 16, 17, 18, 19, 20]);
    let model = GasMembrane::new(
        &patch,
        source_gas(0.75, 0.002),
        &reference,
        &fixed,
        80000.0,
        [0.0007, 0.0015],
        || false,
    )
    .unwrap();
    let mut u = reference.clone();
    u[5][0] += 0.0004;
    u[5][1] -= 0.0003;
    u[5][2] = 0.006;
    u[6][0] -= 0.0006;
    u[6][1] += 0.0002;
    u[6][2] = 0.009;
    // Exact balancing dead forces give a controlled, genuinely curved loaded
    // equilibrium. Every in-plane and normal interior coordinate remains free.
    let initial = model
        .evaluate(&u, &vec![[0.0; 3]; u.len()], gas_controls(), || false)
        .unwrap();
    let vibration = model
        .solve_vibration(
            &initial.residual_n,
            &u,
            controls(),
            modal_controls(),
            || false,
        )
        .unwrap();
    let basis = vibration.basis();
    let state = &vibration.solution().equilibrium;
    assert_eq!(basis.modes().len(), 6);
    let mut geometric_norm = 0.0_f64;
    for mode in basis.modes() {
        let action = model
            .tangent_apply(&state.displacement_m, &mode.phi, gas_controls(), || false)
            .unwrap();
        let mut residual = 0.0_f64;
        let mut scale = 0.0_f64;
        for ((&dof, &value), &phi) in basis.free_dofs().iter().zip(&action).zip(&mode.phi) {
            let inertia = mode.lambda * patch.nodal_mass_kg()[dof / 3] * phi;
            residual = residual.hypot(value - inertia);
            scale = scale.hypot(inertia);
        }
        assert!(
            residual <= 2e-8 * scale,
            "eigen residual {residual:e} / {scale:e}"
        );
        let mut direction = vec![[0.0; 3]; u.len()];
        for (&dof, &phi) in basis.free_dofs().iter().zip(&mode.phi) {
            direction[dof / 3][dof % 3] = phi;
        }
        let geometric = patch
            .pressure_tangent_apply(
                &state.displacement_m,
                &direction,
                state.pressure.pressure_pa - 80000.0,
                || false,
            )
            .unwrap();
        for &dof in basis.free_dofs() {
            geometric_norm = geometric_norm.max(geometric[dof / 3][dof % 3].abs());
        }
    }
    assert!(
        geometric_norm > 1.0,
        "geometric pressure stiffness must be exercised"
    );
    let coordinates = [1e-7, -2e-7, 3e-7, 2e-7, -1e-7, 4e-7];
    near(
        vibration.pressure_perturbation_pa(&coordinates).unwrap(),
        state.pressure.pressure_volume_derivative_pa_m3
            * vibration.volume_perturbation_m3(&coordinates).unwrap(),
        1e-12,
        2e-12,
    );
}

#[test]
fn real_gas_modes_keep_the_original_departure_bulk_modulus() {
    let ideal = source_gas(0.75, 0.001);
    let source = ThermochemicalPengRobinsonModel::new(
        ideal.model().clone(),
        vec![
            PengRobinsonSpeciesBinding {
                species: SpeciesId::new("A").unwrap(),
                critical_temperature: Temperature::new(190.564),
                critical_pressure: Pressure::new(4.5992e6),
                acentric_factor: 0.01142,
            },
            PengRobinsonSpeciesBinding {
                species: SpeciesId::new("B").unwrap(),
                critical_temperature: Temperature::new(150.687),
                critical_pressure: Pressure::new(4.863e6),
                acentric_factor: -0.00219,
            },
        ],
        vec![0.0, 0.03, 0.03, 0.0],
    )
    .unwrap();
    let point = QueryPoint::new()
        .with_quantity(
            "temperature",
            quantity(QuantityKind::AbsoluteTemperature),
            300.0,
        )
        .unwrap()
        .with_quantity("pressure", quantity(QuantityKind::Pressure), 5e6)
        .unwrap()
        .with("specimen-lot", 42.0)
        .unwrap();
    let source = ClosedRealGasParcel::new(source, point, V0, PengRobinsonPhaseV1::Vapor).unwrap();
    let pressure = source.pressure_response().unwrap();
    let ideal =
        ClosedGasParcel::new(source.model().source().clone(), source.point().clone(), V0).unwrap();
    assert!(
        (pressure.bulk_modulus_pa - ideal.pressure_response().unwrap().bulk_modulus_pa).abs() > 1.0
    );
    let (patch, reference, fixed) = pyramid();
    let model = GasMembrane::new(
        &patch,
        source.clone(),
        &reference,
        &fixed,
        pressure.pressure_pa,
        [0.0007, 0.0015],
        || false,
    )
    .unwrap();
    let settings = GasMembraneSolve {
        gas: GasEntropySolve {
            minimum_temperature_k: 280.0,
            maximum_temperature_k: 450.0,
            ..gas_controls()
        },
        ..controls()
    };
    let vibration = model
        .solve_vibration(
            &[[0.0; 3]; 5],
            &reference,
            settings,
            modal_controls(),
            || false,
        )
        .unwrap();
    let expected = (2.0 * STIFFNESS * (STRETCH.powi(2) - 1.0)
        + pressure.bulk_modulus_pa / V0 * volume_gradient().powi(2))
        / patch.nodal_mass_kg()[0];
    near(vibration.basis().modes()[0].lambda, expected, 1e-5, 2e-12);
    assert_eq!(
        vibration.solution().equilibrium.gas,
        MembraneGas::PengRobinson(source.clone())
    );
    assert_eq!(vibration.solution().equilibrium.pressure, pressure);
    assert_eq!(model.reference_gas(), &MembraneGas::PengRobinson(source));
}

#[test]
fn local_modal_admission_refuses_boundaries_budget_and_cancellation_atomically() {
    let (patch, reference, fixed) = pyramid();
    let source = source_gas(0.75, 0.0);
    let model = GasMembrane::new(
        &patch,
        source.clone(),
        &reference,
        &fixed,
        P0,
        [0.0007, 0.0015],
        || false,
    )
    .unwrap();
    let loads = [[0.0; 3]; 5];
    let endpoint = GasMembraneSolve {
        gas: GasEntropySolve {
            minimum_temperature_k: 400.0,
            ..gas_controls()
        },
        ..controls()
    };
    assert!(model.solve(&loads, &reference, endpoint, || false).is_ok());
    assert!(matches!(
        model.solve_vibration(&loads, &reference, endpoint, modal_controls(), || false),
        Err(GasMembraneError::Invalid(_))
    ));
    let volume_endpoint = GasMembrane::new(
        &patch,
        source.clone(),
        &reference,
        &fixed,
        P0,
        [V0, 0.0015],
        || false,
    )
    .unwrap();
    assert!(
        volume_endpoint
            .solve(&loads, &reference, controls(), || false)
            .is_ok()
    );
    assert!(matches!(
        volume_endpoint.solve_vibration(&loads, &reference, controls(), modal_controls(), || false),
        Err(GasMembraneError::Invalid(_))
    ));
    assert!(matches!(
        model.solve_vibration(
            &loads,
            &reference,
            controls(),
            VibrationSettings {
                max_free_dofs: 0,
                ..modal_controls()
            },
            || false,
        ),
        Err(GasMembraneError::Membrane(MembraneError::Budget(_)))
    ));
    assert!(matches!(
        model.solve_vibration(&loads, &reference, controls(), modal_controls(), || true),
        Err(GasMembraneError::Cancelled)
    ));
    let admitted = model
        .solve_vibration(&loads, &reference, controls(), modal_controls(), || false)
        .unwrap();
    assert!(admitted.volume_perturbation_m3(&[]).is_err());
    assert!(admitted.pressure_perturbation_pa(&[f64::NAN]).is_err());
    assert_eq!(model.reference_gas(), &MembraneGas::Ideal(source));
}

#[test]
fn shared_modal_owner_rejects_unstable_and_nonconservative_pencils() {
    let (patch, reference, _) = pyramid();
    assert!(matches!(
        patch.modes_from_tangent(&reference, &[2], &[-1.0], modal_controls(), || false),
        Err(VibrationError::UnresolvedMode { .. })
    ));
    let mut outside = reference.clone();
    outside[0][2] = 1.0;
    assert!(matches!(
        patch.modes_from_tangent(&outside, &[2], &[1.0], modal_controls(), || false),
        Err(VibrationError::Membrane(MembraneError::State { .. }))
    ));
    let free = [15, 16, 17, 18, 19, 20];
    let (patch, reference, _) = chart(Mesh2::triangles(0.3, 0.2, 3, 2), &free);
    let mut matrix = vec![0.0; 36];
    for i in 0..6 {
        matrix[6 * i + i] = 1.0;
    }
    matrix[1] = 1.0;
    assert!(matches!(
        patch.modes_from_tangent(&reference, &free, &matrix, modal_controls(), || false),
        Err(VibrationError::Membrane(MembraneError::Invalid(_)))
    ));
}
