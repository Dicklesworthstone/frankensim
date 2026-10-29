//! Independent species, physical-energy and characteristic-network checks.
use super::{
    BranchEnd, MixtureNetwork1d, NetworkStepReport, Plenum, PlenumConnection, PortExchange,
};
use crate::multicomponent::characteristic_tests::caloric_gas;
use crate::multicomponent::{
    Boundaries, Boundary, Cell, Conserved, Error, Mesh1d, MixtureEuler1d, MixtureGas, StepControl,
    StepReport,
};
use fs_qty::{Length, Time};

fn walls() -> Boundaries {
    Boundaries::Ends {
        left: Boundary::Reflecting,
        right: Boundary::Reflecting,
    }
}
fn branch(
    gas: &MixtureGas,
    n: usize,
    area: f64,
    initial: impl Fn(f64) -> ([f64; 2], f64, f64),
) -> MixtureEuler1d {
    let states = (0..n)
        .map(|i| {
            let (d, u, p) = initial((i as f64 + 0.5) / n as f64);
            gas.conserved(&d, u, p).unwrap()
        })
        .collect();
    MixtureEuler1d::new(
        gas.clone(),
        Length::new(1.0 / n as f64),
        area,
        states,
        walls(),
        || false,
    )
    .unwrap()
}
fn owned(cell: Cell<'_>) -> Conserved {
    Conserved {
        partial_densities: cell.partial_densities.to_vec(),
        momentum: cell.momentum,
        energy: cell.energy,
    }
}
fn close(actual: f64, expected: f64, scale: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * scale,
        "actual={actual:.15e}, expected={expected:.15e}, scale={scale:.9e}, tolerance={tolerance:.3e}"
    );
}
// Integrate the actual fields, including the source's formation references;
// this deliberately does not call network.totals or reuse report residuals.
fn integrated(network: &MixtureNetwork1d) -> ([f64; 2], f64) {
    let mut mass = [0.0; 2];
    let mut physical_energy = 0.0;
    for branch in network.branches() {
        for (i, volume) in branch.mesh().volumes().iter().enumerate() {
            let q = branch.cell(i).unwrap();
            for k in 0..2 {
                mass[k] += volume * q.partial_densities[k];
            }
            physical_energy += volume * branch.gas().absolute_energy_density(q).unwrap();
        }
    }
    for plenum in network.plenums() {
        let q = plenum.state();
        for k in 0..2 {
            mass[k] += plenum.volume_m3() * q.partial_densities[k];
        }
        physical_energy += plenum.volume_m3() * network.gas().absolute_energy_density(q).unwrap();
    }
    (mass, physical_energy)
}
fn connected(gas: &MixtureGas, n: usize, reversed: bool) -> MixtureNetwork1d {
    let inlet = branch(gas, n, 0.1, |x| {
        let x = if reversed { 1.0 - x } else { x };
        ([1.2, 0.0], 0.0, if x > 0.75 { 1.6 } else { 1.0 })
    });
    let daughter = |area| branch(gas, n, area, |_| ([0.0, 2.0], 0.0, 1.0));
    MixtureNetwork1d::new(
        vec![inlet, daughter(0.04), daughter(0.06)],
        vec![Plenum::new(0.03, gas.conserved(&[0.0, 2.0], 0.0, 1.0).unwrap()).unwrap()],
        (0..3)
            .map(|index| PlenumConnection {
                branch: index,
                end: if (index == 0) != reversed {
                    BranchEnd::Right
                } else {
                    BranchEnd::Left
                },
                plenum: 0,
            })
            .collect(),
        || false,
    )
    .unwrap()
}

#[test]
fn different_compositions_and_tapered_branches_preserve_static_pressure_equilibrium() {
    let gas = caloric_gas(true, 0.0);
    let a = MixtureEuler1d::new_on_mesh(
        gas.clone(),
        Mesh1d::new(
            vec![0.0, 0.2, 0.4, 1.0],
            vec![0.1, 0.2, 0.15, 0.1],
            vec![0.03, 0.035, 0.075],
        )
        .unwrap(),
        [[1.2, 0.0], [0.5, 0.7], [0.0, 1.5]]
            .map(|d| gas.conserved(&d, 0.0, 1.0).unwrap())
            .to_vec(),
        walls(),
        || false,
    )
    .unwrap();
    let b = branch(&gas, 8, 0.1, |_| ([0.8, 0.2], 0.0, 1.0));
    let mut network = MixtureNetwork1d::new(
        vec![a, b],
        vec![Plenum::new(0.03, gas.conserved(&[0.1, 0.9], 0.0, 1.0).unwrap()).unwrap()],
        vec![
            PlenumConnection {
                branch: 0,
                end: BranchEnd::Right,
                plenum: 0,
            },
            PlenumConnection {
                branch: 1,
                end: BranchEnd::Left,
                plenum: 0,
            },
        ],
        || false,
    )
    .unwrap();
    let before = integrated(&network);
    let original_plenum = owned(network.plenums()[0].state());
    for _ in 0..20 {
        let report = network
            .step(Time::new(0.003), StepControl::default(), || false)
            .unwrap();
        assert_eq!(&report.boundary_change.species_masses[..2], &[0.0; 2]);
        assert_eq!(report.boundary_change.energy, 0.0);
        for branch in network.branches() {
            for i in 0..branch.cell_count() {
                let primitive = gas.primitive(branch.cell(i).unwrap(), || false).unwrap();
                close(primitive.pressure, 1.0, 1.0, 2e-12);
                close(primitive.velocity, 0.0, 1.0, 2e-12);
            }
        }
    }
    let after = integrated(&network);
    for k in 0..2 {
        close(after.0[k], before.0[k], 1.0, 2e-13);
        close(
            network.plenums()[0].state().partial_densities[k],
            original_plenum.partial_densities[k],
            1.0,
            2e-12,
        );
    }
    close(after.1, before.1, before.1.abs().max(1.0), 2e-12);
}

#[test]
fn finite_pulse_mixes_species_feeds_both_daughters_and_closes_physical_storage() {
    let gas = caloric_gas(true, 5.0);
    let mut network = connected(&gas, 32, false);
    let initial = integrated(&network);
    let mut peaks = [1.0_f64; 2];
    let end = 0.35;
    for _ in 0..10_000 {
        if network.time().value() >= end {
            break;
        }
        let old = owned(network.plenums()[0].state());
        let report = network
            .step(
                Time::new(end - network.time().value()),
                StepControl::default(),
                || false,
            )
            .unwrap();
        let ports = network.accepted_port_exchanges().unwrap();
        let branch_reports = network.accepted_branch_reports().unwrap();
        let q = network.plenums()[0].state();
        let volume = network.plenums()[0].volume_m3();
        for k in 0..2 {
            let received: f64 = ports.iter().map(|p| p.species_masses_to_plenum[k]).sum();
            close(
                volume * (q.partial_densities[k] - old.partial_densities[k]),
                received,
                1.0,
                2e-13,
            );
            assert_eq!(report.boundary_change.species_masses[k], 0.0);
            close(
                report.conservation_residual.species_masses[k],
                0.0,
                1.0,
                2e-13,
            );
        }
        let energy_from_ports: f64 = ports
            .iter()
            .map(|p| p.absolute_energy_to_plenum(&gas).unwrap())
            .sum();
        let actual_energy_change = volume
            * (gas.absolute_energy_density(q).unwrap()
                - gas.absolute_energy_density(old.as_cell()).unwrap());
        close(actual_energy_change, energy_from_ports, 1.0, 4e-12);
        assert_eq!(q.momentum, 0.0);
        assert_eq!(report.boundary_change.energy, 0.0);
        assert!(report.maximum_cfl <= 0.400_000_000_000_01);
        for (i, connection) in network.connections().iter().enumerate() {
            let owner = &branch_reports[connection.branch];
            let end = match connection.end {
                BranchEnd::Left => owner.left_boundary_change,
                BranchEnd::Right => owner.right_boundary_change,
            };
            for k in 0..2 {
                assert_eq!(ports[i].species_masses_to_plenum[k], -end.species_masses[k]);
            }
            assert_eq!(ports[i].energy_to_plenum, -end.energy);
            assert_eq!(ports[i].momentum_to_plenum, -end.momentum);
        }
        for (j, branch) in network.branches()[1..].iter().enumerate() {
            for i in 0..branch.cell_count() {
                let q = branch.cell(i).unwrap();
                assert!(q.partial_densities.iter().all(|v| *v >= 0.0));
                let primitive = gas.primitive(q, || false).unwrap();
                peaks[j] = peaks[j].max(primitive.pressure);
                assert!(primitive.temperature > 0.001 && primitive.temperature < 50.0);
            }
        }
    }
    assert_eq!(network.time().value(), end);
    let after = integrated(&network);
    for k in 0..2 {
        close(after.0[k], initial.0[k], 1.0, 3e-12);
    }
    close(after.1, initial.1, initial.1.abs().max(1.0), 3e-11);
    let reported = network.totals();
    for k in 0..2 {
        close(reported.species_masses[k], after.0[k], 1.0, 2e-13);
    }
    close(
        reported.absolute_energy(&gas).unwrap(),
        after.1,
        after.1.abs().max(1.0),
        2e-12,
    );
    assert!(
        network.plenums()[0].state().partial_densities[0] > 1e-4,
        "feeder species did not enter the plenum"
    );
    for (j, branch) in network.branches()[1..].iter().enumerate() {
        assert!(
            peaks[j] > 1.02,
            "daughter {j} received no finite pressure wave: {peaks:?}"
        );
        assert!(
            branch.totals().species_masses[0] > 1e-6,
            "daughter {j} received no feeder species"
        );
    }
    eprintln!(
        "mixture splitter peaks={peaks:?}, plenum species={:?}",
        network.plenums()[0].state().partial_densities
    );
}

#[test]
fn compact_species_volume_controls_the_shared_step_and_incompatible_sources_refuse() {
    let gas = caloric_gas(false, 0.0);
    let make = || branch(&gas, 4, 0.1, |_| ([1.0, 0.0], 0.0, 1.0));
    let connections = vec![
        PlenumConnection {
            branch: 0,
            end: BranchEnd::Right,
            plenum: 0,
        },
        PlenumConnection {
            branch: 1,
            end: BranchEnd::Left,
            plenum: 0,
        },
    ];
    let mut network = MixtureNetwork1d::new(
        vec![make(), make()],
        vec![Plenum::new(1e-5, gas.conserved(&[0.0, 1.01], 0.0, 1.01).unwrap()).unwrap()],
        connections.clone(),
        || false,
    )
    .unwrap();
    let report = network
        .step(Time::new(0.1), StepControl::default(), || false)
        .unwrap();
    assert!(report.elapsed.value() > 0.0 && report.elapsed.value() < 1e-4);
    assert!(report.maximum_cfl <= 0.400_000_000_000_01);
    for branch in network.branches() {
        assert_eq!(branch.time(), network.time());
    }
    for receipt in network.accepted_port_exchanges().unwrap() {
        assert!(receipt.mass_to_plenum() < 0.0);
        assert_eq!(receipt.species_masses_to_plenum[0], 0.0);
    }
    let other = caloric_gas(true, 0.0);
    let incompatible = branch(&other, 4, 0.1, |_| ([1.0, 0.0], 0.0, 1.0));
    assert!(matches!(
        MixtureNetwork1d::new(vec![make(), incompatible], Vec::new(), Vec::new(), || false),
        Err(Error::Invalid(_))
    ));
    assert!(matches!(
        MixtureNetwork1d::new(
            vec![make()],
            vec![Plenum::new(0.1, gas.conserved(&[0.0, 1.0], 0.0, 1.0).unwrap()).unwrap()],
            vec![connections[0], connections[0]],
            || false
        ),
        Err(Error::Invalid(_))
    ));
}

type Snapshot = (
    Vec<Vec<Conserved>>,
    Vec<Conserved>,
    f64,
    Option<NetworkStepReport>,
    Option<Vec<PortExchange>>,
    Option<Vec<StepReport>>,
);
fn snapshot(network: &MixtureNetwork1d) -> Snapshot {
    (
        network
            .branches()
            .iter()
            .map(|b| {
                (0..b.cell_count())
                    .map(|i| owned(b.cell(i).unwrap()))
                    .collect()
            })
            .collect(),
        network.plenums().iter().map(|p| owned(p.state())).collect(),
        network.time().value(),
        network.last_report(),
        network.accepted_port_exchanges().map(<[_]>::to_vec),
        network.accepted_branch_reports().map(<[_]>::to_vec),
    )
}

#[test]
fn mixture_network_cancellation_is_atomic_and_orientation_and_formation_gauges_preserve_physics() {
    let gas = caloric_gas(true, 0.0);
    let mut probe = connected(&gas, 8, false);
    probe
        .step(Time::new(0.002), StepControl::default(), || false)
        .unwrap();
    let mut count = 0;
    probe
        .step(Time::new(0.002), StepControl::default(), || {
            count += 1;
            false
        })
        .unwrap();
    assert!(count > 100);
    for stop in [1, count / 2, count] {
        let mut network = connected(&gas, 8, false);
        network
            .step(Time::new(0.002), StepControl::default(), || false)
            .unwrap();
        let before = snapshot(&network);
        let mut visits = 0;
        let result = network.step(Time::new(0.002), StepControl::default(), || {
            visits += 1;
            visits >= stop
        });
        assert_eq!(
            result,
            Err(Error::Cancelled),
            "stop={stop}, visits={visits}"
        );
        assert_eq!(snapshot(&network), before);
        network
            .step(Time::new(0.002), StepControl::default(), || false)
            .unwrap();
        assert_eq!(network.time(), probe.time());
        // Reused caloric guesses may end on different final rounding bits;
        // accepted physical state must still agree within the source budget.
        for (a, b) in network.branches().iter().zip(probe.branches()) {
            for i in 0..a.cell_count() {
                let qa = a.cell(i).unwrap();
                let qb = b.cell(i).unwrap();
                for k in 0..2 {
                    close(qa.partial_densities[k], qb.partial_densities[k], 1.0, 2e-12);
                }
                close(qa.momentum, qb.momentum, 1.0, 2e-12);
                close(qa.energy, qb.energy, qa.energy.abs().max(1.0), 2e-12);
            }
        }
        let a = network.plenums()[0].state();
        let b = probe.plenums()[0].state();
        for k in 0..2 {
            close(a.partial_densities[k], b.partial_densities[k], 1.0, 2e-12);
        }
        close(a.energy, b.energy, a.energy.abs().max(1.0), 2e-12);
    }
    let mut ordinary = connected(&gas, 8, false);
    for _ in 0..10 {
        ordinary
            .step(Time::new(0.002), StepControl::default(), || false)
            .unwrap();
    }
    for (reference, reversed) in [(-1e18, false), (1e18, false), (0.0, true)] {
        let shifted = caloric_gas(true, reference);
        let mut other = connected(&shifted, 8, reversed);
        for _ in 0..10 {
            other
                .step(Time::new(0.002), StepControl::default(), || false)
                .unwrap();
        }
        if !reversed {
            assert_eq!(snapshot(&other), snapshot(&ordinary));
        } else {
            assert_eq!(other.time(), ordinary.time());
            for (a, b) in ordinary.branches().iter().zip(other.branches()) {
                for i in 0..a.cell_count() {
                    let qa = a.cell(i).unwrap();
                    let qb = b.cell(a.cell_count() - 1 - i).unwrap();
                    for k in 0..2 {
                        close(qa.partial_densities[k], qb.partial_densities[k], 1.0, 2e-11);
                    }
                    close(qa.momentum, -qb.momentum, 1.0, 2e-11);
                    close(qa.energy, qb.energy, qa.energy.abs().max(1.0), 2e-11);
                }
            }
            for (a, b) in ordinary
                .accepted_port_exchanges()
                .unwrap()
                .iter()
                .zip(other.accepted_port_exchanges().unwrap())
            {
                for k in 0..2 {
                    close(
                        a.species_masses_to_plenum[k],
                        b.species_masses_to_plenum[k],
                        1.0,
                        2e-12,
                    );
                }
                close(a.energy_to_plenum, b.energy_to_plenum, 1.0, 2e-12);
                close(a.momentum_to_plenum, -b.momentum_to_plenum, 1.0, 2e-12);
            }
        }
    }
}
