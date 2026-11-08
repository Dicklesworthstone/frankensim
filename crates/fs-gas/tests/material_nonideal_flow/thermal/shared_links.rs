//! G1/G4: ordered heat links retain one original inventory per wall region.

use super::*;
use fs_gas::dynamics::nonideal::SharedWallLink;

const DT: f64 = 1e-5;
const MASS: f64 = 0.75;

fn links() -> Vec<SharedWallLink> {
    vec![
        SharedWallLink {
            walls: [0, 1],
            conductance_w_per_k: 1e8,
        },
        SharedWallLink {
            walls: [1, 2],
            conductance_w_per_k: 5e7,
        },
    ]
}

fn make(linked: bool) -> (RealGasEuler1d, RealGasWallOperator<'static>) {
    let source = MaterialWallCaloric::from_thermal_law(
        LatentEnthalpy,
        wall_cell(100_000.0),
        [0.0, 200_000.0],
    )
    .unwrap();
    let latent = RealGasWallContact {
        contact: source.initial_contact(0.0, WallDrag::None).unwrap(),
        heat_storage: Some(source),
    };
    let regions = [finite_wall(600.0, 0.0), latent, finite_wall(400.0, 0.0)]
        .into_iter()
        .enumerate()
        .map(|(cell, input)| {
            let mut region = shared_wall(input, &[2e6]);
            region.ports[0].gas = SharedWallGas::Branch { branch: 0, cell };
            region
        })
        .collect();
    let flow = uniform(gas(), 500.0, 0.0);
    let walls = RealGasWallOperator::new(&flow, vec![contact(WallThermal::Adiabatic); 3], || false)
        .unwrap()
        .with_shared_walls(regions, SharedWallControl::default(), || false)
        .unwrap();
    let walls = if linked {
        walls.with_shared_wall_links(links(), || false).unwrap()
    } else {
        walls
    };
    (flow, walls)
}

fn states(walls: &RealGasWallOperator<'_>) -> [WallThermal; 3] {
    [0, 1, 2].map(|i| walls.shared_wall(i).unwrap().storage.thermal)
}

fn enthalpy(thermal: WallThermal) -> f64 {
    let WallThermal::Caloric {
        mass_kg,
        specific_energy_j_kg,
        ..
    } = thermal
    else {
        panic!("original finite caloric inventory required");
    };
    assert_eq!(mass_kg, MASS);
    specific_energy_j_kg
}

// Independent closed charts, not a production material inverse or capacity cache.
fn temperature(region: usize, h: f64) -> f64 {
    if region != 1 {
        return (2.0 * h).sqrt();
    }
    if h < 50_000.0 {
        450.0 + h / 500.0
    } else if h <= 150_000.0 {
        550.0
    } else {
        550.0 + (h - 150_000.0) / 500.0
    }
}

#[test]
fn nonlinear_and_latent_regions_exchange_original_heat_then_feed_back_into_gas() {
    let (mut unlinked, mut independent) = make(false);
    let (mut flow, mut walls) = make(true);
    let heat_controls = WallControl {
        relative_energy_tolerance: 1e-12,
        ..wall_controls()
    };
    let wave = unlinked
        .step_with_walls(
            Time::new(DT),
            controls(),
            &mut independent,
            heat_controls,
            || false,
        )
        .unwrap();
    let report = flow
        .step_with_walls(Time::new(DT), controls(), &mut walls, heat_controls, || {
            false
        })
        .unwrap();
    assert_eq!(report.elapsed, Time::new(DT));
    assert_eq!(flow.states(), unlinked.states());
    assert_eq!(flow.primitives(), unlinked.primitives());
    assert_eq!(report.wall_heat, wave.wall_heat);
    assert_eq!(walls.shared_wall_links().unwrap(), links());
    let entering = states(&independent);
    let mut stage = entering;
    let mut expected = entering.map(enthalpy);
    let receipts = walls.shared_wall_link_reports().unwrap();
    assert_eq!(receipts.len(), 2);
    for (link, receipt) in links().iter().zip(receipts) {
        let [a, b] = link.walls;
        assert_eq!(receipt.elapsed, report.elapsed);
        assert_eq!(receipt.initial_state, [stage[a], stage[b]]);
        let integrated = report.elapsed.value() * link.conductance_w_per_k;
        let initial_flux = integrated * (temperature(a, expected[a]) - temperature(b, expected[b]));
        assert!(initial_flux > 0.0);
        let heat = root(0.0, initial_flux, |q| {
            q - integrated
                * (temperature(a, expected[a] - q / MASS) - temperature(b, expected[b] + q / MASS))
        });
        expected[a] -= heat / MASS;
        expected[b] += heat / MASS;
        near(receipt.heat_a_to_b_j, heat, 1e-6);
        for (side, region) in [a, b].into_iter().enumerate() {
            near(
                enthalpy(receipt.endpoint_state[side]),
                expected[region],
                2e-6,
            );
            near(
                receipt.endpoint_temperature_k[side],
                temperature(region, expected[region]),
                1e-8,
            );
            near(
                receipt.wall_energy_change_j[side],
                if side == 0 { -heat } else { heat },
                1e-6,
            );
            stage[region] = receipt.endpoint_state[side];
        }
        let flux = integrated
            * (receipt.endpoint_temperature_k[0] - receipt.endpoint_temperature_k[1]);
        let scale = receipt.wall_energy_change_j.into_iter()
            .fold(flux.abs(), |scale, change| scale.max(change.abs()));
        let allowance = heat_controls.absolute_energy_tolerance_j
            + heat_controls.relative_energy_tolerance * scale;
        near(receipt.energy_residual_j, 0.0, allowance);
        near(receipt.flux_residual_j, 0.0, allowance);
    }
    assert_eq!(receipts[1].initial_state[0], receipts[0].endpoint_state[1]);
    assert_ne!(receipts[1].initial_state[0], entering[1]);
    assert_eq!(states(&walls), stage);
    assert!((50_000.0..150_000.0).contains(&expected[1]));
    assert!((expected[1] - enthalpy(entering[1])).abs() > 1000.0);
    near(receipts[0].endpoint_temperature_k[1], 550.0, 0.0);
    near(receipts[1].endpoint_temperature_k[0], 550.0, 0.0);
    let actual_storage: f64 = states(&walls)
        .map(enthalpy)
        .into_iter()
        .zip([180_000.0, 100_000.0, 80_000.0])
        .map(|(h, initial)| MASS * (h - initial))
        .sum();
    near(
        walls.last_report().unwrap().wall_storage_change_j,
        actual_storage,
        3e-7,
    );
    assert_eq!(walls.last_report().unwrap().external_heat_to_gas_j, 0.0);
    check_wall_balance(&flow, &walls, report);

    // The next gas stage sees the changed finite regions; the first gas stage
    // did not silently become a simultaneous solve or receive link heat twice.
    let before = states(&walls).map(enthalpy);
    unlinked
        .step_with_walls(
            Time::new(DT),
            controls(),
            &mut independent,
            heat_controls,
            || false,
        )
        .unwrap();
    let next = flow
        .step_with_walls(Time::new(DT), controls(), &mut walls, heat_controls, || {
            false
        })
        .unwrap();
    assert!(unlinked.primitives()[0].temperature - flow.primitives()[0].temperature > 0.1);
    assert!(flow.primitives()[2].temperature - unlinked.primitives()[2].temperature > 0.1);
    let storage: f64 = states(&walls)
        .map(enthalpy)
        .into_iter()
        .zip(before)
        .map(|(h, old)| MASS * (h - old))
        .sum();
    near(
        walls.last_report().unwrap().wall_storage_change_j,
        storage,
        3e-7,
    );
    near(
        next.after.energy - next.before.energy + storage - next.boundary_change.energy,
        0.0,
        5e-7,
    );
    near(
        walls.last_report().unwrap().gas_wall_storage_residual_j,
        0.0,
        3e-7,
    );
    assert_eq!(flow.time(), Time::new(2.0 * DT));
}

#[test]
fn late_link_cancellation_and_invalid_links_preserve_accepted_owners() {
    let (mut good, mut good_walls) = make(true);
    let (mut flow, mut walls) = make(true);
    let mut no_retry = controls();
    no_retry.maximum_halvings = 0;
    let heat_controls = WallControl {
        relative_energy_tolerance: 1e-12,
        ..wall_controls()
    };
    for (flow, walls) in [(&mut good, &mut good_walls), (&mut flow, &mut walls)] {
        flow.step_with_walls(Time::new(DT), no_retry, walls, heat_controls, || false)
            .unwrap();
    }
    let snapshot = (
        flow.states().to_vec(),
        flow.primitives().to_vec(),
        flow.time(),
        states(&walls),
        walls.last_report(),
        walls.last_shared_wall_report(),
        walls.shared_wall_link_reports().unwrap().to_vec(),
    );
    let mut polls = 0;
    good.step_with_walls(
        Time::new(DT),
        no_retry,
        &mut good_walls,
        heat_controls,
        || {
            polls += 1;
            false
        },
    )
    .unwrap();
    let mut seen = 0;
    assert_eq!(
        flow.step_with_walls(Time::new(DT), no_retry, &mut walls, heat_controls, || {
            seen += 1;
            seen == polls
        }),
        Err(RealGasFlowError::Cancelled)
    );
    assert_eq!(
        (
            flow.states().to_vec(),
            flow.primitives().to_vec(),
            flow.time(),
            states(&walls),
            walls.last_report(),
            walls.last_shared_wall_report(),
            walls.shared_wall_link_reports().unwrap().to_vec()
        ),
        snapshot
    );

    for link in [
        SharedWallLink {
            walls: [0, 0],
            conductance_w_per_k: 1.0,
        },
        SharedWallLink {
            walls: [0, 3],
            conductance_w_per_k: 1.0,
        },
        SharedWallLink {
            walls: [0, 1],
            conductance_w_per_k: -1.0,
        },
        SharedWallLink {
            walls: [0, 1],
            conductance_w_per_k: f64::NAN,
        },
    ] {
        let (_, fresh) = make(false);
        assert!(fresh.with_shared_wall_links(vec![link], || false).is_err());
    }
    let (_, fresh) = make(false);
    assert!(fresh.with_shared_wall_links(Vec::new(), || false).is_err());
    assert!(walls.with_shared_wall_links(links(), || false).is_err());
}
