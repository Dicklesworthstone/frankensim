//! Independent calorimetric and original-inventory shared-wall regressions.

use super::*;
use crate::dynamics::hydrostatic::IsothermalAtmosphere;
use crate::dynamics::network::{NetworkGravity, Plenum};
use crate::dynamics::{Boundaries, Boundary, Reconstruction, StepControl};
use fs_qty::{Acceleration, Length};
use fs_thermal::enthalpy::{
    Boundary as ThermalBoundary, CellGeometry, CellId, DepositedCell, DepositionError, MaterialId,
    ThermalLaw, ThermalStorageBasis,
};

pub(super) fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "actual={a:.16e}, expected={b:.16e}, tolerance={tolerance:.3e}"
    );
}
pub(super) fn gas(linear: f64, formation: f64) -> FrozenGas {
    let original = crate::dynamics::tests::gas(linear, formation);
    let mut domain = original.domain();
    domain.absolute_energy_tolerance = 1e-13;
    domain.relative_energy_tolerance = 1e-14;
    FrozenGas::new(original.source().clone(), domain, || false).unwrap()
}
pub(super) fn controls() -> WallControl {
    WallControl {
        absolute_energy_tolerance_j: 1e-11,
        relative_energy_tolerance: 1e-12,
        ..WallControl::default()
    }
}
pub(super) fn periodic(gas: &FrozenGas, temperature: f64, velocity: f64) -> Euler1d {
    Euler1d::new(
        gas.clone(),
        Length::new(0.25),
        0.7,
        vec![gas.conserved(1.0, velocity, temperature).unwrap(); 4],
        Boundaries::Periodic,
        || false,
    )
    .unwrap()
}
pub(super) fn base(network: &GasNetwork1d, drag: WallDrag) -> NetworkWallOperator<'static> {
    NetworkWallOperator::new(
        network,
        vec![
            PlenumThermalContact {
                thermal: WallThermal::Adiabatic,
                heat_storage: None,
            };
            network.plenums().len()
        ],
        || false,
    )
    .unwrap()
    .with_branch_contacts(
        network,
        network
            .branches()
            .iter()
            .map(|b| {
                Some(vec![
                    BranchWallContact {
                        contact: WallContact {
                            drag,
                            thermal: WallThermal::Adiabatic
                        },
                        heat_storage: None
                    };
                    b.states().len()
                ])
            })
            .collect(),
        || false,
    )
    .unwrap()
}
pub(super) fn port(branch: usize, cell: usize, conductance: f64) -> SharedWallPort {
    SharedWallPort {
        gas: SharedWallGas::Branch { branch, cell },
        conductance_w_per_k: conductance,
    }
}
pub(super) fn lumped(
    temperature: f64,
    capacity: f64,
    ports: Vec<SharedWallPort>,
) -> SharedWall<'static> {
    SharedWall {
        storage: PlenumThermalContact {
            thermal: WallThermal::Lumped {
                temperature_k: temperature,
                heat_capacity_j_per_k: capacity,
                conductance_w_per_k: 0.0,
            },
            heat_storage: None,
        },
        ports,
    }
}
pub(super) fn wall_temperature(walls: &NetworkWallOperator<'_>, index: usize) -> f64 {
    let wall = &walls.shared_wall(index).unwrap().storage;
    match &wall.heat_storage {
        Some(source) => source.temperature_k(wall.thermal).unwrap(),
        None => match wall.thermal {
            WallThermal::Lumped { temperature_k, .. } => temperature_k,
            _ => panic!("finite wall"),
        },
    }
}
pub(super) fn attach(
    network: &GasNetwork1d,
    drag: WallDrag,
    walls: Vec<SharedWall<'static>>,
) -> NetworkWallOperator<'static> {
    base(network, drag)
        .with_shared_walls(network, walls, SharedWallControl::default(), || false)
        .unwrap()
}
pub(super) fn step(
    network: &mut GasNetwork1d,
    walls: &mut NetworkWallOperator<'_>,
    dt: f64,
) -> crate::dynamics::network::NetworkStepReport {
    network
        .step_with_walls(
            Time::new(dt),
            StepControl::default(),
            walls,
            controls(),
            || false,
        )
        .unwrap()
}
pub(super) fn fields(network: &GasNetwork1d) -> Vec<Vec<Conserved>> {
    network
        .branches()
        .iter()
        .map(|b| b.states().to_vec())
        .collect()
}
fn both_streams(formation: f64) -> (GasNetwork1d, NetworkWallOperator<'static>) {
    let gas = gas(0.4, formation);
    let network = GasNetwork1d::new(
        vec![periodic(&gas, 1.0, 0.0), periodic(&gas, 5.0, 0.0)],
        vec![],
        vec![],
        || false,
    )
    .unwrap();
    let ports = (0..2)
        .flat_map(|b| (0..4).map(move |c| port(b, c, 1e4)))
        .collect();
    let walls = attach(&network, WallDrag::None, vec![lumped(3.0, 1.4, ports)]);
    (network, walls)
}

#[test]
fn one_inventory_matches_independent_three_body_backward_euler_with_drag() {
    let gas = gas(0.0, 0.0);
    let mut network = GasNetwork1d::new(
        vec![periodic(&gas, 2.0, 0.7)],
        vec![Plenum::new(0.3, gas.conserved(1.4, 0.0, 7.0).unwrap()).unwrap()],
        vec![],
        || false,
    )
    .unwrap();
    let mut walls = attach(
        &network,
        WallDrag::Linear { rate_per_s: 2.0 },
        vec![lumped(
            3.5,
            0.9,
            vec![
                port(0, 0, 2.3),
                SharedWallPort {
                    gas: SharedWallGas::Plenum(0),
                    conductance_w_per_k: 1.1,
                },
            ],
        )],
    );
    let before = network.totals().energy + 0.9 * 3.5;
    let report = step(&mut network, &mut walls, 1.0);
    let dt = report.elapsed.value();
    assert!(
        dt < 1.0,
        "wall must use actual CFL-accepted transport duration"
    );
    let v = 0.7 * (-2.0 * dt).exp();
    let t1 = 2.0 + 0.5 * (0.7 * 0.7 - v * v) / 2.5;
    let (c1, c2, cw) = (0.175 * 2.5, 0.3 * 1.4 * 2.5, 0.9);
    let (a1, a2) = (dt * 2.3, dt * 1.1);
    let (k1, k2) = (a1 * c1 / (a1 + c1), a2 * c2 / (a2 + c2));
    let tw = (cw * 3.5 + k1 * t1 + k2 * 5.0) / (cw + k1 + k2);
    let (tg1, tg2) = (
        (c1 * t1 + a1 * tw) / (c1 + a1),
        (c2 * 5.0 + a2 * tw) / (c2 + a2),
    );
    let actual1 = gas
        .primitive(network.branches()[0].states()[0], || false)
        .unwrap();
    let actual2 = gas
        .primitive(network.plenums()[0].state(), || false)
        .unwrap();
    close(actual1.temperature, tg1, 2e-10);
    close(actual1.velocity, v, 2e-14);
    close(actual2.temperature, tg2, 2e-10);
    close(wall_temperature(&walls, 0), tw, 2e-10);
    let heat = report.thermal.unwrap().source;
    close(
        heat.heat_to_gas_j,
        c1 * (tg1 - t1) + c2 * (tg2 - 5.0),
        2e-10,
    );
    close(heat.wall_energy_change_j, cw * (tw - 3.5), 2e-10);
    assert_eq!(heat.external_heat_to_gas_j, 0.0);
    close(heat.drag_heating_j, 0.5 * 0.7 * (0.7 * 0.7 - v * v), 3e-14);
    let branch = network.accepted_branch_reports().unwrap()[0];
    close(branch.wall_heat, c1 * (tg1 - t1), 2e-10);
    close(branch.wall_momentum_impulse, 0.7 * (v - 0.7), 3e-14);
    close(branch.conservation_residual.energy, 0.0, 1e-13);
    close(
        network.totals().energy + cw * wall_temperature(&walls, 0),
        before,
        3e-10,
    );
    let coupled = walls.last_shared_wall_report().unwrap();
    assert!(coupled.contact_solves > 2 && coupled.source_evaluations > coupled.contact_solves);
    assert!(coupled.source_evaluations <= SharedWallControl::default().maximum_source_evaluations);
    close(coupled.energy_residual_j, 0.0, 2e-11);
    close(coupled.maximum_contact_residual_j, 0.0, 2e-11);
}

#[test]
fn stiff_shared_wall_reaches_original_nonlinear_equilibrium_without_crossing_or_gauge_loss() {
    let (mut network, mut walls) = both_streams(0.0);
    let (mut shifted, mut shifted_walls) = both_streams(1e18);
    let initial = network.totals().energy + 1.4 * 3.0;
    // Independently integrate Cv=2.5+0.4T for BOTH distinct gas amounts.
    let a: f64 = 1.4 * 0.2;
    let b: f64 = 1.4 * 2.5 + 1.4;
    let total = 0.7 * (2.5 * 1.0 + 0.2 * 1.0) + 0.7 * (2.5 * 5.0 + 0.2 * 25.0) + 1.4 * 3.0;
    let equilibrium = 2.0 * total / (b + (b * b + 4.0 * a * total).sqrt());
    let entropy = |n: &GasNetwork1d, w: &NetworkWallOperator<'_>| {
        let gas_s: f64 = n
            .branches()
            .iter()
            .flat_map(|b| b.states().iter().zip(b.mesh().volumes()))
            .map(|(q, v)| {
                let t = n.gas().primitive(*q, || false).unwrap().temperature;
                q.density * v * (2.5 * t.ln() + 0.4 * t)
            })
            .sum();
        gas_s + 1.4 * wall_temperature(w, 0).ln()
    };
    let mut previous_entropy = entropy(&network, &walls);
    for _ in 0..8 {
        let report = step(&mut network, &mut walls, 0.01);
        let shifted_report = step(&mut shifted, &mut shifted_walls, 0.01);
        assert_eq!(fields(&network), fields(&shifted));
        assert_eq!(walls.shared_wall(0), shifted_walls.shared_wall(0));
        assert_eq!(report.thermal, shifted_report.thermal);
        let tw = wall_temperature(&walls, 0);
        assert!((1.0..=5.0).contains(&tw));
        for branch in network.branches() {
            for q in branch.states() {
                let t = network.gas().primitive(*q, || false).unwrap().temperature;
                assert!((1.0 - 1e-12..=5.0 + 1e-12).contains(&t));
                assert_eq!(q.momentum, 0.0);
            }
        }
        let next_entropy = entropy(&network, &walls);
        assert!(next_entropy >= previous_entropy - 2e-11);
        previous_entropy = next_entropy;
        close(network.totals().energy + 1.4 * tw, initial, 3e-10);
        close(
            report.thermal.unwrap().gas_wall_conservation_residual_j,
            0.0,
            3e-11,
        );
    }
    close(wall_temperature(&walls, 0), equilibrium, 3e-8);
    for branch in network.branches() {
        for q in branch.states() {
            close(
                network.gas().primitive(*q, || false).unwrap().temperature,
                equilibrium,
                3e-8,
            );
        }
    }
}

#[test]
fn shared_wall_retains_exact_hydrostatic_rest_at_large_common_potential_datum() {
    for reconstruction in [Reconstruction::FirstOrder, Reconstruction::Muscl] {
        let gas = gas(0.4, 1e18);
        let mesh = Mesh1d::new(
            vec![0.0, 0.1, 0.3, 0.6, 1.0],
            vec![0.7; 5],
            vec![0.07, 0.14, 0.21, 0.28],
        )
        .unwrap();
        let branch = Euler1d::new_isothermal_atmosphere(
            gas.clone(),
            mesh,
            IsothermalAtmosphere {
                acceleration: Acceleration::new(-0.5),
                temperature: Temperature::new(2.0),
                reference_pressure: Pressure::new(3.0),
                reference_position: Length::new(0.0),
            },
            Boundaries::Ends {
                left: Boundary::Reflecting,
                right: Boundary::Reflecting,
            },
            || false,
        )
        .unwrap();
        let node = branch.hydrostatic_reference_face(0).unwrap().0;
        let mut network = GasNetwork1d::new_with_gravity(
            vec![branch],
            vec![Plenum::new(0.3, node).unwrap()],
            vec![],
            NetworkGravity {
                branch_offsets_j_kg: vec![2.0_f64.powi(40)],
                plenum_potentials_j_kg: vec![2.0_f64.powi(40)],
            },
            || false,
        )
        .unwrap();
        let mut ports: Vec<_> = (0..4).map(|c| port(0, c, 1e9)).collect();
        ports.push(SharedWallPort {
            gas: SharedWallGas::Plenum(0),
            conductance_w_per_k: 1e9,
        });
        let mut walls = attach(
            &network,
            WallDrag::Linear { rate_per_s: 1e4 },
            vec![lumped(2.0, 0.7, ports)],
        );
        let initial = fields(&network);
        let initial_wall = walls.shared_wall(0).unwrap().clone();
        for _ in 0..20 {
            let report = network
                .step_with_walls(
                    Time::new(0.003),
                    StepControl {
                        reconstruction,
                        ..StepControl::default()
                    },
                    &mut walls,
                    controls(),
                    || false,
                )
                .unwrap();
            assert_eq!(fields(&network), initial);
            assert_eq!(network.plenums()[0].state(), node);
            assert_eq!(walls.shared_wall(0), Some(&initial_wall));
            let source = report.thermal.unwrap();
            assert_eq!(source.source.heat_to_gas_j, 0.0);
            assert_eq!(source.source.wall_energy_change_j, 0.0);
            assert_eq!(source.source.drag_heating_j, 0.0);
            assert_eq!(source.gas_wall_potential_conservation_residual_j, Some(0.0));
            assert_eq!(
                walls
                    .last_shared_wall_report()
                    .unwrap()
                    .maximum_contact_residual_j,
                0.0
            );
            assert_eq!(network.accepted_branch_reports().unwrap()[0].wall_heat, 0.0);
        }
    }
}

struct Latent;
impl ThermalLaw for Latent {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        if !(0.0..=10.0).contains(&h) {
            return Err(DepositionError::MaterialDomain);
        }
        Ok(if h < 2.0 {
            1.0 + h
        } else if h <= 5.0 {
            3.0
        } else {
            3.0 + 0.5 * (h - 5.0)
        })
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0
    }
}
pub(super) fn latent(range: [f64; 2], ports: Vec<SharedWallPort>) -> SharedWall<'static> {
    let source = MaterialWallCaloric::from_thermal_law(
        Latent,
        DepositedCell {
            id: CellId(77),
            segment_id: 1,
            body_id: 4,
            material_id: MaterialId(29),
            mass_kg: 0.2,
            specific_h_j_kg: 4.0,
            geometry: CellGeometry {
                start_m: [0.0; 3],
                end_m: [1.0, 0.0, 0.0],
                volume_m3: 0.001,
                surface_area_m2: 1.0,
            },
            boundary: ThermalBoundary {
                ambient_k: 3.0,
                htc_w_m2_k: 0.0,
                reservoir: None,
                source_w: 0.0,
            },
        },
        range,
    )
    .unwrap();
    SharedWall {
        storage: PlenumThermalContact {
            thermal: source.initial_contact(0.0, WallDrag::None).unwrap().thermal,
            heat_storage: Some(source),
        },
        ports,
    }
}
fn nonlinear_reservoir_temperature(initial: f64, mass: f64, integrated: f64) -> f64 {
    let b = 2.5 + integrated / mass;
    let rhs = 2.5 * initial + 0.2 * initial * initial + 3.0 * integrated / mass;
    2.0 * rhs / (b + (b * b + 0.8 * rhs).sqrt())
}

#[test]
fn original_latent_inventory_is_debited_once_for_all_gas_interfaces() {
    let gas = gas(0.4, 1e18);
    let mut network = GasNetwork1d::new(
        vec![periodic(&gas, 1.5, 0.0), periodic(&gas, 2.5, 0.0)],
        vec![],
        vec![],
        || false,
    )
    .unwrap();
    let mut walls = attach(
        &network,
        WallDrag::None,
        vec![latent([0.0, 10.0], vec![port(0, 0, 2.0), port(1, 0, 1.7)])],
    );
    let before = network.totals().energy + 0.2 * 4.0;
    let report = step(&mut network, &mut walls, 0.01);
    let dt = report.elapsed.value();
    let mut total_heat = 0.0;
    for (b, initial, g) in [(0, 1.5, 2.0), (1, 2.5, 1.7)] {
        let expected = nonlinear_reservoir_temperature(initial, 0.175, dt * g);
        let actual = gas
            .primitive(network.branches()[b].states()[0], || false)
            .unwrap()
            .temperature;
        close(actual, expected, 2e-10);
        let q =
            0.175 * (2.5 * (expected - initial) + 0.2 * (expected * expected - initial * initial));
        total_heat += q;
        close(
            network.accepted_branch_reports().unwrap()[b].wall_heat,
            q,
            2e-10,
        );
    }
    let wall = walls.shared_wall(0).unwrap();
    let WallThermal::Caloric {
        mass_kg,
        specific_energy_j_kg,
        ..
    } = wall.storage.thermal
    else {
        panic!("original latent energy coordinate");
    };
    assert_eq!(mass_kg, 0.2);
    close(specific_energy_j_kg, 4.0 - total_heat / 0.2, 3e-10);
    assert!((2.0..=5.0).contains(&specific_energy_j_kg));
    assert_eq!(wall_temperature(&walls, 0), 3.0);
    let source = wall.storage.heat_storage.as_ref().unwrap();
    assert_eq!(source.initial_cell().unwrap().id, CellId(77));
    assert_eq!(source.storage_basis(), ThermalStorageBasis::Enthalpy);
    close(
        network.totals().energy + 0.2 * specific_energy_j_kg,
        before,
        3e-10,
    );
    close(
        report.thermal.unwrap().source.wall_energy_change_j,
        -total_heat,
        2e-10,
    );
}

#[test]
fn common_wall_heat_drives_real_field_motion_and_retains_closed_energy() {
    let gas = gas(0.4, 0.0);
    let mut network = GasNetwork1d::new(
        vec![periodic(&gas, 1.0, 0.0), periodic(&gas, 5.0, 0.0)],
        vec![],
        vec![],
        || false,
    )
    .unwrap();
    let mut walls = attach(
        &network,
        WallDrag::Linear { rate_per_s: 0.3 },
        vec![lumped(3.0, 0.6, vec![port(0, 0, 2.0), port(1, 0, 2.0)])],
    );
    let before = network.totals().energy + 0.6 * 3.0;
    let mut drag_heat = 0.0;
    for _ in 0..30 {
        let report = step(&mut network, &mut walls, 0.002);
        drag_heat += report.thermal.unwrap().source.drag_heating_j;
        close(
            report.thermal.unwrap().gas_wall_conservation_residual_j,
            0.0,
            2e-11,
        );
        for branch in network.accepted_branch_reports().unwrap() {
            close(branch.conservation_residual.energy, 0.0, 2e-12);
            close(branch.conservation_residual.momentum, 0.0, 2e-12);
        }
    }
    assert!(
        fields(&network)
            .iter()
            .flatten()
            .any(|q| q.momentum.abs() > 1e-6)
    );
    assert!(drag_heat > 0.0);
    close(
        network.totals().energy + 0.6 * wall_temperature(&walls, 0),
        before,
        5e-10,
    );
}

#[test]
fn shared_work_domain_and_duplicate_refusals_preserve_all_accepted_storage() {
    let gas = gas(0.4, 0.0);
    let mut network = GasNetwork1d::new(
        vec![periodic(&gas, 1.5, 0.0), periodic(&gas, 2.0, 0.0)],
        vec![],
        vec![],
        || false,
    )
    .unwrap();
    assert!(
        base(&network, WallDrag::None)
            .with_shared_walls(
                &network,
                vec![lumped(3.0, 0.7, vec![port(0, 0, 1.0), port(0, 0, 1.0)])],
                SharedWallControl::default(),
                || false
            )
            .is_err()
    );
    let mut overlaps = base(&network, WallDrag::None);
    overlaps.contacts[0].thermal = WallThermal::Lumped {
        temperature_k: 3.0,
        heat_capacity_j_per_k: 0.7,
        conductance_w_per_k: 1.0,
    };
    assert!(
        overlaps
            .with_shared_walls(
                &network,
                vec![lumped(3.0, 0.7, vec![port(0, 0, 1.0)])],
                SharedWallControl::default(),
                || false
            )
            .is_err()
    );
    let original = fields(&network);
    let mut exhausted = base(&network, WallDrag::None)
        .with_shared_walls(
            &network,
            vec![lumped(3.0, 0.7, vec![port(0, 0, 1.0), port(1, 0, 1.0)])],
            SharedWallControl {
                maximum_source_evaluations: 1,
            },
            || false,
        )
        .unwrap();
    let step_control = StepControl {
        maximum_halvings: 0,
        ..StepControl::default()
    };
    assert!(
        network
            .step_with_walls(
                Time::new(0.01),
                step_control,
                &mut exhausted,
                controls(),
                || false
            )
            .is_err()
    );
    assert_eq!(fields(&network), original);
    assert_eq!(network.time(), Time::new(0.0));
    assert_eq!(exhausted.last_shared_wall_report(), None);
    let mut walls = attach(
        &network,
        WallDrag::None,
        vec![
            lumped(4.0, 0.7, vec![port(0, 0, 1.0), port(0, 1, 1.0)]),
            latent([3.999, 4.0], vec![port(1, 0, 10.0), port(1, 1, 10.0)]),
        ],
    );
    let stored: Vec<_> = (0..2)
        .map(|i| walls.shared_wall(i).unwrap().clone())
        .collect();
    assert!(
        network
            .step_with_walls(
                Time::new(0.01),
                step_control,
                &mut walls,
                controls(),
                || false
            )
            .is_err()
    );
    assert_eq!(fields(&network), original);
    assert_eq!(network.time(), Time::new(0.0));
    for (i, wall) in stored.iter().enumerate() {
        assert_eq!(walls.shared_wall(i), Some(wall));
    }
    assert_eq!(walls.last_shared_wall_report(), None);
    assert_eq!(network.last_report(), None);
}

#[test]
fn final_shared_source_cancellation_preserves_history_and_replays() {
    let make = || {
        let gas = gas(0.4, 0.0);
        let network = GasNetwork1d::new(
            vec![periodic(&gas, 1.5, 0.1), periodic(&gas, 4.0, -0.1)],
            vec![],
            vec![],
            || false,
        )
        .unwrap();
        let walls = attach(
            &network,
            WallDrag::Linear { rate_per_s: 0.2 },
            vec![lumped(3.0, 0.8, vec![port(0, 0, 2.0), port(1, 0, 1.7)])],
        );
        (network, walls)
    };
    let (mut good, mut good_walls) = make();
    let (mut cancelled, mut cancelled_walls) = make();
    step(&mut good, &mut good_walls, 0.003);
    step(&mut cancelled, &mut cancelled_walls, 0.003);
    let old_fields = fields(&cancelled);
    let old_wall = cancelled_walls.shared_wall(0).unwrap().clone();
    let old_report = cancelled.last_report();
    let old_source = cancelled_walls.last_shared_wall_report();
    let old_branches = cancelled.accepted_branch_reports().unwrap().to_vec();
    let old_time = cancelled.time();
    let mut calls = 0;
    good.step_with_walls(
        Time::new(0.003),
        StepControl::default(),
        &mut good_walls,
        controls(),
        || {
            calls += 1;
            false
        },
    )
    .unwrap();
    assert!(calls > 20);
    let mut count = 0;
    let refused = cancelled.step_with_walls(
        Time::new(0.003),
        StepControl::default(),
        &mut cancelled_walls,
        controls(),
        || {
            count += 1;
            count == calls
        },
    );
    assert_eq!(refused, Err(Error::Cancelled));
    assert_eq!(fields(&cancelled), old_fields);
    assert_eq!(cancelled_walls.shared_wall(0), Some(&old_wall));
    assert_eq!(cancelled.last_report(), old_report);
    assert_eq!(cancelled_walls.last_shared_wall_report(), old_source);
    assert_eq!(cancelled.accepted_branch_reports().unwrap(), old_branches);
    assert_eq!(cancelled.time(), old_time);
    step(&mut cancelled, &mut cancelled_walls, 0.003);
    assert_eq!(fields(&cancelled), fields(&good));
    assert_eq!(cancelled_walls.shared_wall(0), good_walls.shared_wall(0));
    assert_eq!(cancelled.last_report(), good.last_report());
    assert_eq!(
        cancelled_walls.last_shared_wall_report(),
        good_walls.last_shared_wall_report()
    );
}
