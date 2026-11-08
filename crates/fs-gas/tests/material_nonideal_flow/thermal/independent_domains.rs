//! G1/G4: independent wall/gas domains retain the original coupled heat balance.

use super::material_conduction::temperature_source;
use super::*;
use fs_gas::multicomponent::network::wall_exchange::{
    MaterialWallConduction, WallConductionEvaluation,
};
use fs_gas::multicomponent::wall_exchange::WallExchangeError;
use fs_qty::Area;

fn attach(
    flow: &RealGasEuler1d,
    input: RealGasWallContact<'static>,
    films: &[f64],
) -> RealGasWallOperator<'static> {
    RealGasWallOperator::new(flow, vec![contact(WallThermal::Adiabatic); 3], || false)
        .unwrap()
        .with_shared_walls(
            vec![shared_wall(input, films)],
            SharedWallControl::default(),
            || false,
        )
        .unwrap()
}

fn step(
    flow: &mut RealGasEuler1d,
    walls: &mut RealGasWallOperator<'_>,
    dt: f64,
) -> Result<StepReport, RealGasFlowError> {
    flow.step_with_walls(
        Time::new(dt),
        StepControl {
            maximum_halvings: 0,
            ..controls()
        },
        walls,
        WallControl {
            relative_energy_tolerance: 1e-12,
            ..wall_controls()
        },
        || false,
    )
}

fn wall_temperature(walls: &RealGasWallOperator<'_>) -> f64 {
    let storage = &walls.shared_wall(0).unwrap().storage;
    if let Some(source) = &storage.heat_storage {
        return source.temperature_k(storage.thermal).unwrap();
    }
    let WallThermal::Lumped { temperature_k, .. } = storage.thermal else {
        panic!("finite wall required")
    };
    temperature_k
}

fn refused(
    flow: &mut RealGasEuler1d,
    walls: &mut RealGasWallOperator<'_>,
    dt: f64,
) -> RealGasFlowError {
    let saved = (
        flow.states().to_vec(),
        flow.primitives().to_vec(),
        flow.time(),
        walls.contacts().to_vec(),
        walls.shared_wall(0).unwrap().clone(),
        walls.last_report(),
        walls.last_shared_wall_report(),
    );
    let error = step(flow, walls, dt).unwrap_err();
    assert!(matches!(error, RealGasFlowError::StepRetries { .. }));
    assert_eq!(
        (
            flow.states().to_vec(),
            flow.primitives().to_vec(),
            flow.time(),
            walls.contacts().to_vec(),
            walls.shared_wall(0).unwrap().clone(),
            walls.last_report(),
            walls.last_shared_wall_report()
        ),
        saved
    );
    error
}

#[test]
fn small_cold_and_hot_walls_outside_gas_support_have_valid_coupled_endpoints() {
    let (dt, film) = (1e-6, 1e10);
    for (wall0, gas0) in [(250.0, 350.0), (1200.0, 900.0)] {
        let input = |capacity| {
            contact(WallThermal::Lumped {
                temperature_k: wall0,
                heat_capacity_j_per_k: capacity,
                conductance_w_per_k: 0.0,
            })
        };
        let mut flow = uniform(gas(), gas0, 0.0);
        let mass = 20.0 * flow.mesh().volumes()[0];
        let held = root(wall0.min(gas0), wall0.max(gas0), |t| {
            mass * (internal_energy(20.0, t) - internal_energy(20.0, gas0))
                - dt * film * (wall0 - t)
        });
        assert!(
            !(300.0..=1000.0).contains(&held),
            "initial held-wall trial must be unsupported"
        );
        let wall_at =
            |t: f64, capacity: f64| t + capacity / (capacity + 3.0 * dt * film) * (wall0 - t);
        let endpoint = |capacity: f64| {
            root(wall0.min(gas0), wall0.max(gas0), |t| {
                3.0 * mass * (internal_energy(20.0, t) - internal_energy(20.0, gas0))
                    + capacity * (wall_at(t, capacity) - wall0)
            })
        };
        let expected = endpoint(1.0);
        assert!((300.0..=1000.0).contains(&expected));
        let mut walls = attach(&flow, input(1.0), &[film; 3]);
        let report = step(&mut flow, &mut walls, dt).unwrap();
        assert_eq!(report.elapsed, Time::new(dt));
        for primitive in flow.primitives() {
            near(primitive.temperature, expected, 5e-8);
        }
        near(wall_temperature(&walls), wall_at(expected, 1.0), 5e-8);
        near(report.wall_heat, wall0 - wall_at(expected, 1.0), 3e-7);
        check_wall_balance(&flow, &walls, report);

        // A genuine unsupported endpoint must not be clipped into the gas
        // chart or accepted because the wall has its own wider domain.
        assert!(!(300.0..=1000.0).contains(&endpoint(1e8)));
        let mut outside = uniform(gas(), gas0, 0.0);
        let mut large = attach(&outside, input(1e8), &[film; 3]);
        refused(&mut outside, &mut large, dt);

        let mut insulated = uniform(gas(), gas0, 0.0);
        let initial = insulated.states().to_vec();
        let mut disconnected = attach(&insulated, input(1.0), &[0.0; 3]);
        let report = step(&mut insulated, &mut disconnected, dt).unwrap();
        assert_eq!(insulated.states(), initial);
        assert_eq!(wall_temperature(&disconnected), wall0);
        assert_eq!(report.wall_heat, 0.0);
    }
}

#[test]
fn lagged_original_conductivity_uses_its_prepared_series_value_and_nonlinear_storage() {
    let (dt, film, wall0, gas0) = (1e-4, 1e10, 250.0, 350.0);
    let (transport, caloric) = temperature_source([0.2, 0.1], [0.0, 1.0], [200.0, 1000.0]);
    let layer = MaterialWallConduction::new(transport, Length::new(1e-5), Area::new(4.0)).unwrap();
    let input = RealGasWallContact {
        contact: caloric.contact(wall0, 0.0, WallDrag::None).unwrap(),
        heat_storage: Some(caloric),
    };
    let mut flow = uniform(gas(), gas0, 0.0);
    let mass = 20.0 * flow.mesh().volumes()[0];
    let conductance = 1.0 / (1.0 / film + 1e-5 / (4.0 * (0.2 + 0.1 * wall0)));
    let gas_at = |wall: f64| {
        root(wall.min(gas0), wall.max(gas0), |t| {
            mass * (internal_energy(20.0, t) - internal_energy(20.0, gas0))
                - dt * conductance * (wall - t)
        })
    };
    assert!(gas_at(wall0) < 300.0);
    let expected_wall = root(wall0, gas0, |wall| {
        0.5 * (wall * wall - wall0 * wall0)
            + 3.0 * mass * (internal_energy(20.0, gas_at(wall)) - internal_energy(20.0, gas0))
    });
    let mut walls = attach(&flow, input.clone(), &[film; 3])
        .with_shared_wall_conduction(vec![vec![Some(layer.clone()); 3]], || false)
        .unwrap();
    let report = step(&mut flow, &mut walls, dt).unwrap();
    assert_eq!(report.elapsed, Time::new(dt));
    for primitive in flow.primitives() {
        near(primitive.temperature, gas_at(expected_wall), 5e-8);
    }
    near(wall_temperature(&walls), expected_wall, 5e-8);
    near(
        walls.last_report().unwrap().wall_storage_change_j,
        0.5 * (expected_wall * expected_wall - wall0 * wall0),
        3e-5,
    );
    assert!(
        walls
            .shared_wall(0)
            .unwrap()
            .ports
            .iter()
            .all(|port| port.conductance_w_per_k == film)
    );
    check_wall_balance(&flow, &walls, report);

    // The supported interval is derived for constant or prepared lagged G.
    // Outside-domain anchors with an implicit/profile law must refuse explicitly.
    for changing in [
        layer
            .clone()
            .with_evaluation(WallConductionEvaluation::ImplicitEndpoint),
        layer.with_integrated_profile().unwrap(),
    ] {
        let mut original = uniform(gas(), gas0, 0.0);
        let mut unsupported = attach(&original, input.clone(), &[film; 3])
            .with_shared_wall_conduction(vec![vec![Some(changing); 3]], || false)
            .unwrap();
        let RealGasFlowError::StepRetries { last, .. } =
            refused(&mut original, &mut unsupported, dt)
        else {
            unreachable!()
        };
        assert!(matches!(*last, RealGasFlowError::Thermal(ref error)
            if matches!(**error, WallExchangeError::Invalid(_))));
    }
}

#[test]
fn nonlinear_wall_bounds_retain_actual_gas_densities_volumes_and_kinetic_energy() {
    let mut gas_domain = domain();
    gas_domain.energy.minimum_temperature_k = 400.0;
    let gas = bind(
        common::model([0.0; 2], [0.0; 2], [200.0, 2000.0]),
        gas_domain,
    );
    let initial = [
        (12.0, 30.0, 450.0),
        (20.0, -10.0, 650.0),
        (30.0, 15.0, 800.0),
    ]
    .into_iter()
    .map(|(rho, v, t)| state(&gas, rho, v, t))
    .collect::<Vec<_>>();
    let mesh = Mesh1d::new(
        vec![0.0, 0.2, 0.5, 1.0],
        vec![0.07; 4],
        vec![0.014, 0.021, 0.035],
    )
    .unwrap();
    let mut wave = RealGasEuler1d::new_on_mesh(
        gas.clone(),
        mesh.clone(),
        initial.clone(),
        RealGasBoundaries::Periodic,
        || false,
    )
    .unwrap();
    let mut flow =
        RealGasEuler1d::new_on_mesh(gas, mesh, initial, RealGasBoundaries::Periodic, || false)
            .unwrap();
    let (dt, wall0) = (1e-6, 350.0);
    let films = [1e9, 2e8, 7e7];
    let mut walls = attach(&flow, finite_wall(wall0, 0.0), &films);
    let transport = wave.step(Time::new(dt), controls(), || false).unwrap();
    assert_eq!(transport.elapsed, Time::new(dt));
    let gas_at = |wall: f64, i: usize| {
        let (rho, t0) = (wave.states()[i].density, wave.primitives()[i].temperature);
        let mass = rho * wave.mesh().volumes()[i];
        root(wall.min(t0), wall.max(t0), |t| {
            mass * (internal_energy(rho, t) - internal_energy(rho, t0)) - dt * films[i] * (wall - t)
        })
    };
    assert!(gas_at(wall0, 0) < 400.0);
    let heat = |wall: f64| {
        (0..3)
            .map(|i| dt * films[i] * (wall - gas_at(wall, i)))
            .sum::<f64>()
    };
    let expected = root(wall0, 850.0, |wall| {
        0.375 * (wall * wall - wall0 * wall0) + heat(wall)
    });
    let frozen = root(wall0, 850.0, |wall| {
        0.75 * wall0 * (wall - wall0) + heat(wall)
    });
    assert!(
        (expected - frozen).abs() > 1.0,
        "original nonlinear storage must matter"
    );
    let report = step(&mut flow, &mut walls, dt).unwrap();
    assert_eq!(report.elapsed, Time::new(dt));
    for i in 0..3 {
        assert_eq!(flow.states()[i].density, wave.states()[i].density);
        assert_eq!(flow.states()[i].momentum, wave.states()[i].momentum);
        near(flow.primitives()[i].temperature, gas_at(expected, i), 8e-8);
    }
    near(wall_temperature(&walls), expected, 8e-8);
    near(report.wall_heat, heat(expected), 5e-5);
    near(
        walls.last_report().unwrap().wall_storage_change_j,
        0.375 * (expected * expected - wall0 * wall0),
        5e-5,
    );
    check_wall_balance(&flow, &walls, report);
}
