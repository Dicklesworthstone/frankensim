//! Physical CHT temperatures are Kelvin; mathematical solver increments are signed.

use super::*;
use fs_lbm::conjugate::CellSink;

fn assert_invalid_field(error: ChtError, expected: &str) {
    match error {
        ChtError::InvalidInput { field, .. } => assert_eq!(field, expected),
        other => panic!("expected invalid {expected}, got {other:?}"),
    }
}

#[test]
fn absolute_temperature_inputs_require_positive_kelvin() {
    let domain = VoxelDomain::new(1, 1, 1, 1.0).unwrap();
    let flow = FlowField::quiescent(&domain);
    let mut faces = [ThermalFace::Adiabatic; 6];
    faces[0] = ThermalFace::Temperature(300.0);
    let valid = ThermalSetup::new(faces);
    let solve = |setup: &ThermalSetup| {
        solve_energy(
            &domain,
            &unit_fluid(),
            &[],
            &flow,
            setup,
            &EnergyConfig::default(),
            &CancelGate::new(),
        )
    };
    for bad in [0.0, -1.0, f64::INFINITY, f64::NAN] {
        for (rule, field) in [
            (ThermalFace::Temperature(bad), "thermal.face.temperature"),
            (
                ThermalFace::Convective {
                    h: 1.0,
                    ambient: bad,
                },
                "thermal.face.ambient",
            ),
            (
                ThermalFace::Inflow { temperature: bad },
                "thermal.face.inflow",
            ),
            (
                ThermalFace::Outflow {
                    backflow_temperature: bad,
                },
                "thermal.face.backflow",
            ),
        ] {
            let mut setup = valid.clone();
            setup.faces[1] = rule;
            assert_invalid_field(solve(&setup).unwrap_err(), field);
        }
        let mut held = valid.clone();
        held.fixed_temperature.push((0, bad));
        assert_invalid_field(solve(&held).unwrap_err(), "thermal.fixed_temperature");
        let mut sunk = valid.clone();
        sunk.cell_sinks.push(CellSink {
            cell: 0,
            conductance_w_k: 1.0,
            temperature_k: bad,
        });
        assert_invalid_field(solve(&sunk).unwrap_err(), "thermal.cell_sinks.temperature_k");
    }
}

#[test]
fn steady_cooling_refuses_nonpositive_endpoints_without_clamping() {
    // A unit cube attached through one half cell has conductance 2 W/K.
    // Cooling is supported, but the physical result must remain above 0 K.
    let domain = VoxelDomain::from_fn(1, 1, 1, 1.0, |_| Voxel::Solid(0)).unwrap();
    let solids = [SolidMaterial::new("unit conductor", 1.0)];
    let flow = FlowField::quiescent(&domain);
    let mut faces = [ThermalFace::Adiabatic; 6];
    faces[0] = ThermalFace::Temperature(300.0);
    let mut setup = ThermalSetup::new(faces);
    let config = EnergyConfig::default();
    let gate = CancelGate::new();
    for power in [-600.0, -1000.0] {
        setup.power_w = vec![power];
        assert_invalid_field(
            solve_energy(
                &domain, &unit_fluid(), &solids, &flow, &setup, &config, &gate,
            )
            .unwrap_err(),
            "thermal.temperature_k",
        );
        assert_invalid_field(
            config
                .solve_with_conductivity(
                    &domain,
                    &unit_fluid(),
                    &solids,
                    &flow,
                    &setup,
                    &[[1.0; 3]],
                    &gate,
                )
                .unwrap_err(),
            "thermal.temperature_k",
        );
        assert_eq!(setup.power_w, [power]);
    }
    setup.power_w = vec![-100.0];
    let solution = solve_energy(
        &domain, &unit_fluid(), &solids, &flow, &setup, &config, &gate,
    )
    .unwrap();
    assert_eq!(solution.temperature, [250.0]);
    assert_eq!(solution.report.balance.boundary_outflow_w, -100.0);
    assert_eq!(solution.report.balance.residual_w, 0.0);
}

#[test]
fn transient_cooling_refuses_before_publishing_a_nonpositive_step() {
    // Isolated unit heat capacity: each one-second step removes 0.4 J.
    let domain = VoxelDomain::from_fn(1, 1, 1, 1.0, |_| Voxel::Solid(0)).unwrap();
    let solids = [SolidMaterial::new("unit thermal mass", 1.0).with_heat_capacity(1.0)];
    let flow = FlowField::quiescent(&domain);
    let mut setup = ThermalSetup::new([ThermalFace::Adiabatic; 6]);
    setup.power_w = vec![-0.4];
    let initial = [1.0];
    let mut config = TransientConfig {
        time_step_s: 1.0,
        steps: 3,
        energy: EnergyConfig::default(),
    };
    let gate = CancelGate::new();
    assert_invalid_field(
        march_energy(
            &domain,
            &unit_fluid(),
            &solids,
            &flow,
            &setup,
            &initial,
            |_| 1.0,
            &config,
            None,
            &gate,
        )
        .unwrap_err(),
        "thermal.temperature_k",
    );
    assert_eq!(initial, [1.0]);
    assert_eq!(setup.power_w, [-0.4]);

    config.steps = 2;
    let accepted = march_energy(
        &domain,
        &unit_fluid(),
        &solids,
        &flow,
        &setup,
        &initial,
        |_| 1.0,
        &config,
        None,
        &gate,
    )
    .unwrap();
    assert!((accepted.temperature[0] - 0.2).abs() < 1e-14);
    assert_eq!(accepted.records.len(), 2);
    assert_eq!(
        accepted.final_energy.as_ref().unwrap().temperature,
        accepted.temperature,
    );
    for record in accepted.records {
        assert!((record.stored_energy_change_j + 0.4).abs() < 1e-14);
        assert!(record.closure_j.abs() < 1e-14);
    }

    config.steps = 0;
    for bad in [0.0, -1.0, f64::INFINITY, f64::NAN] {
        assert_invalid_field(
            march_energy(
                &domain,
                &unit_fluid(),
                &solids,
                &flow,
                &setup,
                &[bad],
                |_| 1.0,
                &config,
                None,
                &gate,
            )
            .unwrap_err(),
            "transient.initial_temperature",
        );
    }
}

#[test]
fn nonlinear_conduction_can_recover_from_a_signed_linear_proposal() {
    // The initial frozen k=1 proposal is 300 - 700/2 = -50 K.
    // The actual k(T) grows on cooling. Its positive physical endpoint obeys
    // 2 * integral_T^300 k(theta) dtheta = 700 W, giving
    // delta + delta^2/200 = 350.
    // Rejecting intermediate linear proposals would lose this supported solution.
    let domain = VoxelDomain::from_fn(1, 1, 1, 1.0, |_| Voxel::Solid(0)).unwrap();
    let solids = [SolidMaterial::new("cold-conducting solid", 1.0)
        .with_conductivity_table(&[(100.0, 3.0), (300.0, 1.0)])];
    let mut faces = [ThermalFace::Adiabatic; 6];
    faces[0] = ThermalFace::Temperature(300.0);
    let mut setup = ThermalSetup::new(faces);
    setup.power_w = vec![-700.0];
    let solution = solve_energy(
        &domain,
        &unit_fluid(),
        &solids,
        &FlowField::quiescent(&domain),
        &setup,
        &EnergyConfig::default(),
        &CancelGate::new(),
    )
    .unwrap();
    let exact = 400.0 - 100.0 * 8.0_f64.sqrt();
    assert!(
        (solution.temperature[0] - exact).abs() < 1e-8,
        "{} vs {exact}",
        solution.temperature[0],
    );
    assert!(solution.report.balance.relative_residual < 1e-9);
}
