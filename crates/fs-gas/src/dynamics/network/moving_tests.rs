//! Compression, heat and gravity checks using independent finite-inventory laws.

use super::*;
use crate::dynamics::{Boundary, StepControl};
use crate::multicomponent::wall_exchange::{WallControl, WallThermal};
use fs_qty::Length;

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:.16e}, expected={expected:.16e}, tolerance={tolerance:.3e}"
    );
}

fn gas(linear_cp: f64, formation: f64) -> FrozenGas {
    let original = crate::dynamics::tests::gas(linear_cp, formation);
    let mut domain = original.domain();
    domain.absolute_energy_tolerance = 1e-13;
    domain.relative_energy_tolerance = 1e-14;
    FrozenGas::new(original.source().clone(), domain, || false).unwrap()
}

fn isolated(linear_cp: f64, formation: f64, temperatures: &[f64]) -> GasNetwork1d {
    let gas = gas(linear_cp, formation);
    let branch = Euler1d::new(
        gas.clone(),
        Length::new(100.0),
        0.001,
        vec![gas.conserved(1.0, 0.0, 2.0).unwrap(); 2],
        Boundaries::Periodic,
        || false,
    )
    .unwrap();
    let plenums = temperatures
        .iter()
        .map(|&t| Plenum::new(0.4, gas.conserved(1.0, 0.0, t).unwrap()).unwrap())
        .collect();
    GasNetwork1d::new(vec![branch], plenums, vec![], || false).unwrap()
}

fn temperature(network: &GasNetwork1d, p: usize) -> f64 {
    network
        .gas()
        .primitive(network.plenums()[p].state(), || false)
        .unwrap()
        .temperature
}

fn connected() -> GasNetwork1d {
    let gas = gas(0.4, -1e18);
    let initial = gas.conserved(1.0, 0.0, 2.0).unwrap();
    let branch = Euler1d::new(
        gas,
        Length::new(0.025),
        0.05,
        vec![initial; 8],
        Boundaries::Ends { left: Boundary::Reflecting, right: Boundary::Reflecting },
        || false,
    )
    .unwrap();
    GasNetwork1d::new(
        vec![branch],
        vec![Plenum::new(0.04, initial).unwrap(), Plenum::new(0.06, initial).unwrap()],
        vec![
            PlenumConnection { branch: 0, end: BranchEnd::Left, plenum: 0 },
            PlenumConnection { branch: 0, end: BranchEnd::Right, plenum: 1 },
        ],
        || false,
    )
    .unwrap()
}

// Manufactured R=1 and Cv=2.5+b*T, integrated independently of NASA routines.
fn isentropic_temperature(linear_cp: f64, volume_ratio: f64) -> f64 {
    let (mut lo, mut hi) = (0.1_f64, 10.0_f64);
    for _ in 0..100 {
        let t = 0.5 * (lo + hi);
        let entropy = 2.5 * (t / 2.0).ln() + linear_cp * (t - 2.0) + volume_ratio.ln();
        if entropy < 0.0 {
            lo = t;
        } else {
            hi = t;
        }
    }
    0.5 * (lo + hi)
}

#[test]
fn g1_compression_and_expansion_converge_to_original_variable_cp_isentropes() {
    for linear_cp in [0.0, 0.4] {
        for ratio in [0.5_f64, 2.0] {
            let mut errors = Vec::new();
            for steps in [20, 40, 80] {
                let mut network = isolated(linear_cp, 1e18, &[2.0]);
                let initial = network.totals();
                let dt = 0.2 / steps as f64;
                let rate = 0.4 * (ratio - 1.0) / 0.2;
                let mut work = 0.0;
                for _ in 0..steps {
                    let report = network.step_with_volume_rates(
                        Time::new(dt), StepControl::default(), &[rate], || false,
                    ).unwrap();
                    assert_eq!(report.elapsed.value(), dt);
                    assert_eq!(report.boundary_change, NetworkTotals::default());
                    close(report.conservation_residual.energy, 0.0, 2e-13);
                    let motion = network.accepted_plenum_boundary_work().unwrap()[0];
                    assert_eq!(motion.final_volume_m3, network.plenums()[0].volume_m3());
                    assert_eq!(motion.work_on_gas_j, report.boundary_work_on_gas_j);
                    work += report.boundary_work_on_gas_j;
                }
                let t = temperature(&network, 0);
                let volume = network.plenums()[0].volume_m3();
                close(volume, 0.4 * ratio, 2e-14);
                close(volume * network.plenums()[0].state().density, 0.4, 2e-14);
                close(network.totals().mass, initial.mass, 3e-14);
                close(network.totals().energy - initial.energy, work, 3e-12);
                let physical_change = 0.4
                    * (2.5 * (t - 2.0) + 0.5 * linear_cp * (t * t - 4.0));
                close(work, physical_change, 3e-12);
                errors.push((t - isentropic_temperature(linear_cp, ratio)).abs());
            }
            assert!(errors[0] / errors[1] > 3.5, "{linear_cp} {ratio}: {errors:?}");
            assert!(errors[1] / errors[2] > 3.5, "{linear_cp} {ratio}: {errors:?}");
        }
    }
}

#[test]
fn g0_compression_drives_port_waves_and_closes_global_mass_and_work() {
    let mut network = connected();
    let before = network.totals();
    let mut work = 0.0;
    let mut delivered = [0.0; 2];
    for _ in 0..40 {
        let report = network.step_with_volume_rates(
            Time::new(0.002), StepControl::default(), &[-0.1, 0.0], || false,
        ).unwrap();
        close(report.conservation_residual.mass, 0.0, 1e-14);
        close(report.conservation_residual.energy, 0.0, 2e-13);
        assert_eq!(report.boundary_change, NetworkTotals::default());
        work += report.boundary_work_on_gas_j;
        for (i, exchange) in network.accepted_port_exchanges().unwrap().iter().enumerate() {
            delivered[i] += exchange.mass_to_plenum;
            let branch = network.accepted_branch_reports().unwrap()[0];
            let boundary = if i == 0 {
                branch.left_boundary_change
            } else {
                branch.right_boundary_change
            };
            assert_eq!(exchange.energy_to_plenum, -boundary.energy);
        }
    }
    assert!(delivered[0] < 0.0, "compressed chamber must discharge");
    assert!(network.branches()[0].states().iter().any(|q| q.momentum > 0.0));
    assert!(work > 0.0);
    close(network.totals().mass, before.mass, 2e-14);
    close(network.totals().energy - before.energy, work, 2e-12);
    for (i, p) in network.plenums().iter().enumerate() {
        let initial_mass = if i == 0 { 0.04 } else { 0.06 };
        close(p.volume_m3() * p.state().density - initial_mass, delivered[i], 2e-14);
    }
}

#[derive(Debug, PartialEq)]
struct Snapshot {
    branches: Vec<Vec<Conserved>>,
    plenums: Vec<Plenum>,
    time: Time,
    report: Option<NetworkStepReport>,
    branches_report: Option<Vec<StepReport>>,
    ports: Option<Vec<PortExchange>>,
    work: Option<Vec<PlenumBoundaryWork>>,
}

fn snapshot(n: &GasNetwork1d) -> Snapshot {
    Snapshot {
        branches: n.branches().iter().map(|b| b.states().to_vec()).collect(),
        plenums: n.plenums().to_vec(),
        time: n.time(),
        report: n.last_report(),
        branches_report: n.accepted_branch_reports().map(<[_]>::to_vec),
        ports: n.accepted_port_exchanges().map(<[_]>::to_vec),
        work: n.accepted_plenum_boundary_work().map(<[_]>::to_vec),
    }
}

#[test]
fn g3_zero_rates_preserve_rigid_arithmetic_and_formation_gauges_do_not_change_compression() {
    let (mut rigid, mut zero) = (connected(), connected());
    for _ in 0..3 {
        rigid.step(Time::new(0.001), StepControl::default(), || false).unwrap();
        zero.step_with_volume_rates(
            Time::new(0.001), StepControl::default(), &[0.0, -0.0], || false,
        ).unwrap();
        assert_eq!(snapshot(&rigid), snapshot(&zero));
    }
    let mut ordinary = isolated(0.4, 0.0, &[2.0]);
    let mut shifted = isolated(0.4, -1e18, &[2.0]);
    for _ in 0..4 {
        ordinary.step_with_volume_rates(
            Time::new(0.01), StepControl::default(), &[-1.0], || false,
        ).unwrap();
        shifted.step_with_volume_rates(
            Time::new(0.01), StepControl::default(), &[-1.0], || false,
        ).unwrap();
        assert_eq!(snapshot(&ordinary), snapshot(&shifted));
    }
}

// Independent closed-cell Heun pressure quadrature with constant R=1, Cv=2.5.
fn compressed_temperature(t: f64, mass: f64, volume: f64, rate: f64, dt: f64) -> f64 {
    let p0 = mass * t / volume;
    let stage_t = t - dt * p0 * rate / (mass * 2.5);
    let p1 = mass * stage_t / (volume + dt * rate);
    t - 0.5 * dt * (p0 + p1) * rate / (mass * 2.5)
}

fn thermal_system(shared: bool) -> (GasNetwork1d, NetworkWallOperator<'static>) {
    let network = isolated(0.0, 1e18, &[2.0, 3.0]);
    let mut walls = NetworkWallOperator::new(
        &network,
        [0.8, 0.5].map(|g| PlenumThermalContact {
            thermal: if shared {
                WallThermal::Adiabatic
            } else {
                WallThermal::Lumped {
                    temperature_k: 1.5,
                    heat_capacity_j_per_k: 1.2,
                    conductance_w_per_k: g,
                }
            },
            heat_storage: None,
        }).to_vec(),
        || false,
    ).unwrap();
    if shared {
        walls = walls.with_shared_walls(
            &network,
            vec![SharedWall {
                storage: PlenumThermalContact {
                    thermal: WallThermal::Lumped {
                        temperature_k: 1.5,
                        heat_capacity_j_per_k: 1.2,
                        conductance_w_per_k: 0.0,
                    },
                    heat_storage: None,
                },
                ports: [0.8, 0.5].into_iter().enumerate().map(|(p, conductance_w_per_k)| {
                    SharedWallPort { gas: SharedWallGas::Plenum(p), conductance_w_per_k }
                }).collect(),
            }],
            SharedWallControl::default(),
            || false,
        ).unwrap();
    }
    (network, walls)
}

fn wall_temperature(thermal: WallThermal) -> f64 {
    match thermal {
        WallThermal::Lumped { temperature_k, .. } => temperature_k,
        _ => panic!("expected finite constant-capacity wall"),
    }
}

fn thermal_control() -> WallControl {
    WallControl {
        absolute_energy_tolerance_j: 1e-11,
        relative_energy_tolerance: 1e-12,
        ..WallControl::default()
    }
}

#[test]
fn g0_moving_chambers_exchange_heat_using_their_actual_endpoint_mass_and_shared_wall() {
    let dt = 0.02;
    let rates = [-1.0, 0.5];
    let conductances = [0.8, 0.5];
    for shared in [false, true] {
        let (mut network, mut walls) = thermal_system(shared);
        for _ in 0..3 {
            let compressed: [f64; 2] = core::array::from_fn(|i| compressed_temperature(
                temperature(&network, i), 0.4, network.plenums()[i].volume_m3(), rates[i], dt,
            ));
            let old_wall: [f64; 2] = core::array::from_fn(|i| wall_temperature(if shared {
                walls.shared_wall(0).unwrap().storage.thermal
            } else {
                walls.contacts()[i].thermal
            }));
            // Cgas=mass*Cv=1, independent of each endpoint volume.
            let effective = conductances.map(|g| dt * g / (1.0 + dt * g));
            let expected_wall: [f64; 2] = core::array::from_fn(|i| if shared {
                (1.2 * old_wall[0] + effective[0] * compressed[0] + effective[1] * compressed[1])
                    / (1.2 + effective[0] + effective[1])
            } else {
                (1.2 * old_wall[i] + effective[i] * compressed[i]) / (1.2 + effective[i])
            });
            let report = network.step_with_walls_and_volume_rates(
                Time::new(dt), StepControl::default(), &rates, &mut walls,
                thermal_control(), || false,
            ).unwrap();
            for i in 0..2 {
                let expected_t = (compressed[i] + dt * conductances[i] * expected_wall[i])
                    / (1.0 + dt * conductances[i]);
                close(temperature(&network, i), expected_t, 2e-10);
                let accepted_wall = wall_temperature(if shared {
                    walls.shared_wall(0).unwrap().storage.thermal
                } else {
                    walls.contacts()[i].thermal
                });
                close(accepted_wall, expected_wall[i], 2e-10);
            }
            close(report.conservation_residual.energy, 0.0, 2e-11);
            close(report.thermal.unwrap().gas_wall_conservation_residual_j, 0.0, 2e-11);
        }
        // The accepted volume binding must also support subsequent rigid steps.
        network.step_with_walls(
            Time::new(dt), StepControl::default(), &mut walls, thermal_control(), || false,
        ).unwrap();
    }
}

#[test]
fn g0_volume_change_preserves_compact_gravitational_mass_and_accounts_for_work_once() {
    for datum in [0.0, 1e12] {
        let original = isolated(0.4, -1e18, &[2.0]);
        let mut network = GasNetwork1d::new_with_gravity(
            original.branches, original.plenums, vec![],
            NetworkGravity {
                branch_offsets_j_kg: vec![datum],
                plenum_potentials_j_kg: vec![datum + 3.0],
            },
            || false,
        ).unwrap();
        for _ in 0..4 {
            let report = network.step_with_volume_rates(
                Time::new(0.02), StepControl::default(), &[-1.0], || false,
            ).unwrap();
            let gravity = report.gravity.unwrap();
            close(gravity.before_centered_potential_energy, 1.2, 1e-14);
            close(gravity.after_centered_potential_energy, 1.2, 1e-14);
            close(gravity.conservation_residual, 0.0, 3e-13);
            assert!(report.boundary_work_on_gas_j > 0.0);
            close(report.conservation_residual.energy, 0.0, 3e-13);
        }
    }
}

#[test]
fn g4_refusal_and_late_cancellation_keep_geometry_heat_and_all_accepted_receipts() {
    let run = |network: &mut GasNetwork1d, walls: &mut NetworkWallOperator<'_>, stop: usize| {
        let mut calls = 0;
        let result = network.step_with_walls_and_volume_rates(
            Time::new(0.02), StepControl::default(), &[-1.0, 0.5], walls, thermal_control(), || {
                calls += 1;
                calls == stop
            },
        );
        (result, calls)
    };
    let (mut probe, mut probe_walls) = thermal_system(true);
    run(&mut probe, &mut probe_walls, usize::MAX).0.unwrap();
    let (result, calls) = run(&mut probe, &mut probe_walls, usize::MAX);
    result.unwrap();
    for stop in [1, calls / 2, calls] {
        let (mut network, mut walls) = thermal_system(true);
        run(&mut network, &mut walls, usize::MAX).0.unwrap();
        let before = snapshot(&network);
        let wall_before = (
            walls.contacts().to_vec(), walls.shared_wall(0).cloned(),
            walls.last_report(), walls.last_shared_wall_report(),
        );
        assert!(matches!(run(&mut network, &mut walls, stop).0, Err(Error::Cancelled)));
        assert_eq!(snapshot(&network), before);
        assert_eq!(
            (walls.contacts().to_vec(), walls.shared_wall(0).cloned(),
                walls.last_report(), walls.last_shared_wall_report()),
            wall_before,
        );
        run(&mut network, &mut walls, usize::MAX).0.unwrap();
        assert_eq!(network.time(), probe.time());
        for i in 0..2 {
            assert_eq!(network.plenums()[i].volume_m3(), probe.plenums()[i].volume_m3());
            close(temperature(&network, i), temperature(&probe, i), 2e-10);
            close(network.plenums()[i].state().density, probe.plenums()[i].state().density, 1e-14);
        }
        close(
            wall_temperature(walls.shared_wall(0).unwrap().storage.thermal),
            wall_temperature(probe_walls.shared_wall(0).unwrap().storage.thermal),
            2e-10,
        );
    }
    let mut network = isolated(0.0, 0.0, &[49.99]);
    let before = snapshot(&network);
    for rates in [&[][..], &[f64::NAN][..], &[f64::MIN_POSITIVE][..], &[-1.0][..]] {
        assert!(network.step_with_volume_rates(
            Time::new(0.02), StepControl { maximum_halvings: 0, ..StepControl::default() },
            rates, || false,
        ).is_err());
        assert_eq!(snapshot(&network), before);
    }
}

#[test]
fn g0_volume_change_cfl_prevents_collapse_and_uses_the_actual_shortened_duration() {
    let mut network = isolated(0.0, 0.0, &[2.0]);
    let report = network.step_with_volume_rates(
        Time::new(1.0), StepControl::default(), &[-20.0], || false,
    ).unwrap();
    assert!(report.elapsed.value() < 0.02);
    assert!(report.maximum_cfl <= StepControl::default().cfl * (1.0 + 1e-13));
    close(network.plenums()[0].volume_m3(), 0.4 - 20.0 * report.elapsed.value(), 1e-15);
    close(network.totals().mass, report.before.mass, 1e-14);
    close(report.conservation_residual.energy, 0.0, 2e-13);
}
