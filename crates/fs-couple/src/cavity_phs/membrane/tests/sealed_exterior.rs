//! Two finite chambers exercise the original equilibrium, modal and time owners.
use super::*;
use fs_solid::membrane::vibration::VibrationSettings;
use fs_time::galpha::ImplicitSolveConfig;

fn exterior(gas: impl Into<MembraneGas>) -> SealedMembraneExterior {
    SealedMembraneExterior {
        gas: gas.into(),
        volume_interval_m3: [0.0007, 0.0015],
        entropy: gas_controls(),
    }
}
fn mode_controls() -> VibrationSettings {
    VibrationSettings {
        max_free_dofs: 16,
        relative_residual_tolerance: 1e-8,
    }
}
fn time_controls() -> GasMembraneMotionControls {
    GasMembraneMotionControls {
        gas: gas_controls(),
        time: ImplicitSolveConfig {
            newton: controls().membrane.newton,
            max_newton_iterations: 40,
        },
        maximum_evaluations: 1000,
        maximum_steps: 32,
        maximum_krylov_restart: 8,
        energy_tolerance_j: 1e-4,
    }
}

#[test]
fn equal_pressure_chambers_add_both_native_stiffnesses_and_opposed_pressure_modes() {
    let (patch, reference, fixed) = pyramid();
    let source = source_gas(0.75, 0.002);
    let model = GasMembrane::with_sealed_exterior(
        &patch,
        source.clone(),
        &reference,
        &fixed,
        [0.0007, 0.0015],
        exterior(source.clone()),
        || false,
    )
    .unwrap();
    let modes = model
        .solve_vibration(
            &[[0.0; 3]; 5],
            &reference,
            controls(),
            mode_controls(),
            || false,
        )
        .unwrap();
    let state = &modes.solution().equilibrium;
    near(
        state.pressure.pressure_pa,
        state.exterior_pressure_pa,
        1e-9,
        1e-14,
    );
    assert_eq!(
        state.exterior_gas.as_ref().unwrap().gas.mass_kg(),
        source.mass_kg()
    );
    // Independent analytic Cp(T) of the synthetic mixture, at T=400 K.
    let gamma = 1.0 + 1.0 / (2.25 + 0.75 * 0.002 * 400.0);
    let gas_stiffness = gamma * P0 / V0 * volume_gradient().powi(2);
    let structural = 2.0 * STIFFNESS * (STRETCH.powi(2) - 1.0);
    near(
        modes.basis().modes()[0].lambda,
        (structural + 2.0 * gas_stiffness) / patch.nodal_mass_kg()[0],
        1e-7,
        3e-12,
    );
    let q = [modes.basis().volume_coupling()[0].signum() * 1e-6];
    let dv = modes.volume_perturbation_m3(&q).unwrap();
    let dp = gamma * P0 / V0 * dv;
    near(
        modes.pressure_perturbation_pa(&q).unwrap(),
        -dp,
        1e-10,
        2e-12,
    );
    near(
        modes.exterior_pressure_perturbation_pa(&q).unwrap(),
        dp,
        1e-10,
        2e-12,
    );
    near(
        modes.pressure_jump_perturbation_pa(&q).unwrap(),
        -2.0 * dp,
        1e-10,
        2e-12,
    );
    let held = GasMembrane::new(
        &patch,
        source,
        &reference,
        &fixed,
        P0,
        [0.0007, 0.0015],
        || false,
    )
    .unwrap();
    let held_modes = held
        .solve_vibration(
            &[[0.0; 3]; 5],
            &reference,
            controls(),
            mode_controls(),
            || false,
        )
        .unwrap();
    assert!(held.sealed_exterior().is_none());
    assert!(held_modes.solution().equilibrium.exterior_gas.is_none());
    assert_eq!(
        held_modes.exterior_pressure_perturbation_pa(&q).unwrap(),
        0.0
    );
    assert_eq!(
        held_modes.pressure_jump_perturbation_pa(&q).unwrap(),
        held_modes.pressure_perturbation_pa(&q).unwrap()
    );
}

#[test]
fn curved_two_chamber_tangent_is_the_gradient_of_both_gas_energies() {
    let (patch, reference, fixed) =
        chart(Mesh2::triangles(0.3, 0.2, 3, 2), &[15, 16, 17, 18, 19, 20]);
    let model = GasMembrane::with_sealed_exterior(
        &patch,
        source_gas(0.75, 0.002),
        &reference,
        &fixed,
        [0.0007, 0.0015],
        exterior(source_gas(0.75, 0.002)),
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
    let loads = vec![[0.0; 3]; u.len()];
    let value = model
        .evaluate(&u, &loads, gas_controls(), || false)
        .unwrap();
    let opposing = value.exterior_gas.as_ref().unwrap();
    near(
        value.gas.volume_m3() + opposing.gas.volume_m3(),
        2.0 * V0,
        1e-18,
        0.0,
    );
    assert!(value.pressure.pressure_pa < P0 && opposing.pressure.pressure_pa > P0);
    assert_eq!(value.exterior_pressure_work_j, 0.0);
    let h = 2e-7;
    let mut geometry = 0.0_f64;
    for (column, &dof) in model.free_dofs().iter().enumerate() {
        let mut lo = u.clone();
        let mut hi = u.clone();
        lo[dof / 3][dof % 3] -= h;
        hi[dof / 3][dof % 3] += h;
        let a = model
            .evaluate(&lo, &loads, gas_controls(), || false)
            .unwrap();
        let b = model
            .evaluate(&hi, &loads, gas_controls(), || false)
            .unwrap();
        near(
            (b.potential_energy_change_j - a.potential_energy_change_j) / (2.0 * h),
            value.residual_n[dof / 3][dof % 3],
            3e-4,
            2e-5,
        );
        let mut direction = vec![0.0; 6];
        direction[column] = 1.0;
        let tangent = model
            .tangent_apply(&u, &direction, gas_controls(), || false)
            .unwrap();
        for (&row, &actual) in model.free_dofs().iter().zip(&tangent) {
            near(
                (b.residual_n[row / 3][row % 3] - a.residual_n[row / 3][row % 3]) / (2.0 * h),
                actual,
                1e-3,
                1e-5,
            );
        }
        let mut full = vec![[0.0; 3]; u.len()];
        full[dof / 3][dof % 3] = 1.0;
        let geometric = patch
            .pressure_tangent_apply(
                &u,
                &full,
                value.pressure.pressure_pa - value.exterior_pressure_pa,
                || false,
            )
            .unwrap();
        for &row in model.free_dofs() {
            geometry = geometry.max(geometric[row / 3][row % 3].abs());
        }
    }
    assert!(
        geometry > 1.0,
        "nonzero follower pressure Hessian is exercised"
    );
}

fn real_gas() -> ClosedRealGasParcel {
    let ideal = source_gas(0.75, 0.001);
    let model = ThermochemicalPengRobinsonModel::new(
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
    ClosedRealGasParcel::new(model, point, V0, PengRobinsonPhaseV1::Vapor).unwrap()
}

#[test]
fn distinct_ideal_and_real_sources_keep_disjoint_temperature_controls_and_mode_admission() {
    let (patch, reference, fixed) = pyramid();
    let first = source_gas(0.75, 0.002);
    let second = real_gas();
    let first_controls = GasEntropySolve {
        minimum_temperature_k: 390.0,
        maximum_temperature_k: 410.0,
        ..gas_controls()
    };
    let second_controls = GasEntropySolve {
        minimum_temperature_k: 290.0,
        maximum_temperature_k: 310.0,
        ..gas_controls()
    };
    let model = GasMembrane::with_sealed_exterior(
        &patch,
        first.clone(),
        &reference,
        &fixed,
        [0.0007, 0.0015],
        SealedMembraneExterior {
            entropy: second_controls,
            ..exterior(second.clone())
        },
        || false,
    )
    .unwrap();
    let mut u = reference.clone();
    u[0][2] = 1e-5;
    let value = model
        .evaluate(&u, &[[0.0; 3]; 5], first_controls, || false)
        .unwrap();
    let mut native_first = first;
    let mut native_second = second.clone();
    let first_step = native_first
        .set_volume_isentropic(value.gas.volume_m3(), first_controls, || false)
        .unwrap();
    let opposing = value.exterior_gas.as_ref().unwrap();
    let second_step = native_second
        .set_volume_isentropic(opposing.gas.volume_m3(), second_controls, || false)
        .unwrap();
    assert_eq!(value.gas, MembraneGas::Ideal(native_first));
    assert_eq!(value.gas_volume, first_step);
    assert_eq!(opposing.gas, MembraneGas::PengRobinson(native_second));
    assert_eq!(opposing.gas_volume, second_step);
    assert_eq!(value.exterior_pressure_pa, opposing.pressure.pressure_pa);
    let loads = model
        .evaluate(&reference, &[[0.0; 3]; 5], first_controls, || false)
        .unwrap()
        .residual_n;
    let settings = GasMembraneSolve {
        gas: first_controls,
        ..controls()
    };
    let modes = model
        .solve_vibration(&loads, &reference, settings, mode_controls(), || false)
        .unwrap();
    assert_eq!(
        modes
            .solution()
            .equilibrium
            .exterior_gas
            .as_ref()
            .unwrap()
            .gas
            .mass_kg(),
        second.mass_kg()
    );
    let endpoint = GasMembrane::with_sealed_exterior(
        &patch,
        source_gas(0.75, 0.002),
        &reference,
        &fixed,
        [0.0007, 0.0015],
        SealedMembraneExterior {
            entropy: GasEntropySolve {
                minimum_temperature_k: 300.0,
                ..second_controls
            },
            ..exterior(second)
        },
        || false,
    )
    .unwrap();
    assert!(matches!(
        endpoint.solve_vibration(&loads, &reference, settings, mode_controls(), || false),
        Err(GasMembraneError::Invalid(_))
    ));
}

#[test]
fn motion_closes_both_gas_accounts_and_opposing_source_refusal_rolls_back_together() {
    let (patch, reference, fixed) = pyramid();
    let source = source_gas(0.75, 0.002);
    let model = GasMembrane::with_sealed_exterior(
        &patch,
        source.clone(),
        &reference,
        &fixed,
        [0.0007, 0.0015],
        exterior(source.clone()),
        || false,
    )
    .unwrap();
    let mut initial = reference.clone();
    initial[0][2] = 0.004;
    let mut velocity = [[0.0; 3]; 5];
    velocity[0][2] = 0.15;
    let loads = [[0.0; 3]; 5];
    let mut motion =
        GasMembraneMotion::new(&model, &initial, &velocity, &loads, time_controls(), || {
            false
        })
        .unwrap();
    for _ in 0..4 {
        let before = motion.accepted().clone();
        let before_kinetic = motion.kinetic_energy_j().unwrap();
        let frame = motion
            .step(1e-5, &loads, time_controls(), || false)
            .unwrap();
        let after = motion.accepted();
        let old_back = before.exterior_gas.as_ref().unwrap();
        let back = after.exterior_gas.as_ref().unwrap();
        near(
            after.gas.volume_m3() + back.gas.volume_m3(),
            2.0 * V0,
            1e-18,
            0.0,
        );
        assert_eq!(after.gas.mass_kg(), source.mass_kg());
        assert_eq!(back.gas.mass_kg(), source.mass_kg());
        near(
            frame.exterior_gas_energy_change_j,
            back.gas_volume.stored_energy_change_j - old_back.gas_volume.stored_energy_change_j,
            1e-14,
            1e-12,
        );
        assert_eq!(
            frame.work_by_exterior_gas_j,
            -frame.exterior_gas_energy_change_j
        );
        assert_eq!(frame.work_by_gas_j, -frame.gas_energy_change_j);
        assert_eq!(frame.exterior_pressure_work_j, 0.0);
        let independently_measured = motion.kinetic_energy_j().unwrap() - before_kinetic
            + after.strain_energy_j
            - before.strain_energy_j
            + after.gas_volume.stored_energy_change_j
            - before.gas_volume.stored_energy_change_j
            + back.gas_volume.stored_energy_change_j
            - old_back.gas_volume.stored_energy_change_j;
        near(
            frame.integration_defect_j,
            independently_measured,
            3e-13,
            1e-8,
        );
        assert!(frame.gas_trials < frame.evaluations);
    }
    // Only the opposing source has the narrow bracket. The first gas admits
    // this shape; failure must not leak its already computed trial state.
    let limited = GasMembrane::with_sealed_exterior(
        &patch,
        source.clone(),
        &reference,
        &fixed,
        [0.0007, 0.0015],
        SealedMembraneExterior {
            entropy: GasEntropySolve {
                minimum_temperature_k: 400.0 - 1e-6,
                maximum_temperature_k: 400.0 + 1e-6,
                ..gas_controls()
            },
            ..exterior(source)
        },
        || false,
    )
    .unwrap();
    let mut refused = reference.clone();
    refused[0][2] = 1e-6;
    assert!(matches!(
        limited.evaluate(&refused, &loads, gas_controls(), || false),
        Err(GasMembraneError::Gas(
            GasEntropyError::OutsideBracket { .. }
        ))
    ));
    velocity[0][2] = 5.0;
    let mut motion = GasMembraneMotion::new(
        &limited,
        &reference,
        &velocity,
        &loads,
        time_controls(),
        || false,
    )
    .unwrap();
    let state = motion.state().clone();
    let accepted = motion.accepted().clone();
    assert!(
        motion
            .step(1e-4, &loads, time_controls(), || false)
            .is_err()
    );
    assert_eq!(motion.state(), &state);
    assert_eq!(motion.accepted().gas, accepted.gas);
    assert_eq!(motion.accepted().gas_volume, accepted.gas_volume);
    assert_eq!(motion.accepted().exterior_gas, accepted.exterior_gas);
    assert_eq!(motion.accepted().displacement_m, accepted.displacement_m);
    assert_eq!(motion.accepted().residual_n, accepted.residual_n);
    assert_eq!(
        motion.accepted().potential_energy_change_j,
        accepted.potential_energy_change_j
    );
}
