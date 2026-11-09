//! Independent acceleration histories, finite-port work and atomic retry checks.

use super::*;
use crate::body_force::{BodyForce, ForceControl};
use crate::multicomponent::network::NetworkBodyForces;

const FORCE: ForceControl = ForceControl { maximum_mach: 100.0 };

fn periodic(gas: &MixtureGas, density: [f64; 2], velocity: f64) -> MixtureEuler1d {
    MixtureEuler1d::new(
        gas.clone(), Length::new(0.125), 0.7,
        vec![gas.conserved(&density, velocity, 1.0).unwrap(); 8],
        Boundaries::Periodic, || false,
    ).unwrap()
}

fn uniform_pair(gas: &MixtureGas, velocity: f64) -> MixtureNetwork1d {
    MixtureNetwork1d::new(
        vec![periodic(gas, [0.8, 0.2], velocity), periodic(gas, [0.3, 1.7], -velocity)],
        vec![], vec![], || false,
    ).unwrap()
}

#[test]
fn uniform_branch_histories_match_independent_impulse_work_and_temperature() {
    let gas = caloric_gas(true, 0.0);
    let mut network = uniform_pair(&gas, 0.3);
    let original = network.branches().iter().map(|b| {
        let p = gas.primitive(b.cell(0).unwrap(), || false).unwrap();
        (p.density, p.velocity, p.temperature, b.totals())
    }).collect::<Vec<_>>();
    let forces = [
        BodyForce::UniformAccelerationRamp { reference_time: 0.0, acceleration: 0.8, jerk: -0.4 },
        BodyForce::UniformForceDensityRamp { reference_time: 0.0, force_density: 0.6, rate: 0.2 },
    ];
    let mut work = 0.0;
    for _ in 0..20 {
        let report = network.step_with_body_forces(
            Time::new(0.002), StepControl::default(),
            NetworkBodyForces { branches: &forces, controls: FORCE }, || false,
        ).unwrap();
        work += report.body_force_work_j;
        close(report.conservation_residual.energy, 0.0, 1.0, 4e-12);
        assert_eq!(report.boundary_change.energy, 0.0);
        assert_eq!(&report.conservation_residual.species_masses[..2], &[0.0; 2]);
        for (index, branch) in network.branches().iter().enumerate() {
            let t = branch.time().value();
            let (rho, u0, temperature, initial) = original[index];
            let du = if index == 0 { 0.8 * t - 0.2 * t * t }
                else { (0.6 * t + 0.1 * t * t) / rho };
            let expected = u0 + du;
            for i in 0..branch.cell_count() {
                let p = gas.primitive(branch.cell(i).unwrap(), || false).unwrap();
                close(p.velocity, expected, 1.0, 5e-13);
                close(p.temperature, temperature, 1.0, 2e-11);
            }
            close(branch.totals().momentum - initial.momentum, initial.mass() * du, 1.0, 3e-13);
        }
    }
    let exact_work: f64 = network.branches().iter().zip(&original).map(|(b, (_, u0, _, initial))| {
        let u = gas.primitive(b.cell(0).unwrap(), || false).unwrap().velocity;
        0.5 * initial.mass() * (u * u - u0 * u0)
    }).sum();
    close(work, exact_work, 1.0, 3e-12);
}

#[test]
fn accelerating_connected_branches_use_actual_stage_mass_and_reciprocal_ports() {
    let gas = caloric_gas(true, 0.0);
    let mut network = connected(&gas, 8, false);
    let forces = [
        BodyForce::UniformAcceleration(3.0),
        BodyForce::UniformForceDensity(-2.0),
        BodyForce::UniformAcceleration(0.7),
    ];
    let initial = integrated(&network);
    let mut source_work = 0.0;
    for _ in 0..12 {
        let report = network.step_with_body_forces(
            Time::new(0.001), StepControl::default(),
            NetworkBodyForces { branches: &forces, controls: FORCE }, || false,
        ).unwrap();
        let h = report.elapsed.value();
        source_work += report.body_force_work_j;
        let branches = network.accepted_branch_reports().unwrap();
        let expected_impulses = [
            1.5 * h * (branches[0].before.mass() + branches[0].after.mass()),
            -2.0 * h * 0.04,
            0.35 * h * (branches[2].before.mass() + branches[2].after.mass()),
        ];
        for (branch, &expected) in branches.iter().zip(&expected_impulses) {
            close(branch.source_change.impulse, expected, 1.0, 3e-14);
            close(branch.conservation_residual.momentum, 0.0, 1.0, 3e-12);
            close(branch.conservation_residual.energy, 0.0, 1.0, 3e-12);
            assert!(branch.maximum_cfl <= StepControl::default().cfl * (1.0 + 1e-12));
        }
        for (port, connection) in network.accepted_port_exchanges().unwrap().iter()
            .zip(network.connections())
        {
            let branch = &branches[connection.branch];
            let end = match connection.end {
                BranchEnd::Left => branch.left_boundary_change,
                BranchEnd::Right => branch.right_boundary_change,
            };
            for k in 0..2 { assert_eq!(port.species_masses_to_plenum[k], -end.species_masses[k]); }
            assert_eq!(port.energy_to_plenum, -end.energy);
            assert_eq!(port.momentum_to_plenum, -end.momentum);
        }
        let summed: f64 = branches.iter().map(|r| r.source_change.work).sum();
        close(report.body_force_work_j, summed, 1.0, 3e-15);
    }
    let after = integrated(&network);
    for k in 0..2 { close(after.0[k], initial.0[k], 1.0, 4e-12); }
    close(after.1 - initial.1, source_work, 1.0, 2e-11);
    assert!(network.plenums()[0].state().partial_densities[0] > 0.0);
}

#[test]
fn accelerated_network_accounts_reservoir_enthalpy_and_moving_plenum_work_once() {
    let gas = caloric_gas(true, 0.0);
    let mut network = connected(&gas, 8, false);
    let bath = gas.conserved(&[1.1, 0.2], 0.0, 1.5).unwrap();
    let reservoirs = [[Some(bath.as_cell()), None], [None; 2], [None; 2]];
    let rates = [-0.004];
    let forces = [BodyForce::UniformAcceleration(1.0); 3];
    let initial = integrated(&network);
    let initial_volume = network.plenums()[0].volume_m3();
    let report = network.step_with_body_forces_and_boundaries(
        Time::new(0.001), StepControl::default(),
        NetworkBodyForces { branches: &forces, controls: FORCE },
        Some(&reservoirs), Some(&rates), || false,
    ).unwrap();
    let after = integrated(&network);
    let boundary_energy = report.boundary_change.absolute_energy(network.gas()).unwrap();
    close(after.1 - initial.1, boundary_energy + report.boundary_work_on_gas_j
        + report.body_force_work_j, 1.0, 4e-12);
    assert!(report.boundary_work_on_gas_j > 0.0);
    close(network.plenums()[0].volume_m3(),
        initial_volume + rates[0] * report.elapsed.value(), 1.0, 2e-16);
    for k in 0..2 {
        close(after.0[k] - initial.0[k], report.boundary_change.species_masses[k], 1.0, 3e-13);
    }
}

#[test]
fn interior_history_speed_peak_shortens_the_whole_network_and_uses_accepted_clock() {
    let gas = caloric_gas(false, 0.0);
    let mut network = uniform_pair(&gas, 20.0);
    let forces = [
        BodyForce::UniformAccelerationRamp {
            reference_time: 0.0, acceleration: 20_000.0, jerk: -40_000_000.0,
        },
        BodyForce::UniformAcceleration(0.0),
    ];
    // First source interval at h=.002 has zero net impulse but reaches +25 m/s.
    // The initial and final 20 m/s pass CFL; that interior peak must not be hidden.
    let report = network.step_with_body_forces(
        Time::new(0.002), StepControl::default(),
        NetworkBodyForces { branches: &forces, controls: FORCE }, || false,
    ).unwrap();
    assert!(report.halvings >= 1);
    assert_eq!(report.elapsed.value(), 0.001);
    for b in network.branches() { assert_eq!(b.time(), report.elapsed); }
    let t = report.elapsed.value();
    let expected = 20.0 + 20_000.0 * t - 20_000_000.0 * t * t;
    let p = gas.primitive(network.branches()[0].cell(0).unwrap(), || false).unwrap();
    close(p.velocity, expected, 1.0, 2e-12);
    let reports = network.accepted_branch_reports().unwrap();
    assert_eq!(reports[0].halvings, reports[1].halvings);
    assert!(report.maximum_cfl <= 0.4 * (1.0 + 1e-12));
}

#[test]
fn zero_force_reuses_original_arithmetic_after_active_buffers_and_refusals_publish_nothing() {
    let gas = caloric_gas(true, 0.0);
    let mut ordinary = connected(&gas, 8, false);
    let mut forced = connected(&gas, 8, false);
    let active = [BodyForce::UniformAcceleration(1.3); 3];
    for network in [&mut ordinary, &mut forced] {
        network.step_with_body_forces(
            Time::new(0.001), StepControl::default(),
            NetworkBodyForces { branches: &active, controls: FORCE }, || false,
        ).unwrap();
    }
    let zero = [BodyForce::UniformAcceleration(0.0); 3];
    let a = ordinary.step(Time::new(0.001), StepControl::default(), || false).unwrap();
    let b = forced.step_with_body_forces(
        Time::new(0.001), StepControl::default(),
        NetworkBodyForces { branches: &zero, controls: FORCE }, || false,
    ).unwrap();
    assert_eq!(a, b);
    assert_eq!(snapshot(&ordinary), snapshot(&forced));
    assert_eq!(b.body_force_work_j, 0.0);

    let before = snapshot(&forced);
    let malformed = [BodyForce::Acceleration(&[1.0]), active[1], active[2]];
    assert!(forced.step_with_body_forces(
        Time::new(0.001), StepControl::default(),
        NetworkBodyForces { branches: &malformed, controls: FORCE }, || false,
    ).is_err());
    assert_eq!(snapshot(&forced), before);
    let mut visits = 0;
    ordinary.step_with_body_forces(
        Time::new(0.001), StepControl::default(),
        NetworkBodyForces { branches: &active, controls: FORCE },
        || { visits += 1; false },
    ).unwrap();
    let mut polls = 0;
    assert_eq!(forced.step_with_body_forces(
        Time::new(0.001), StepControl::default(),
        NetworkBodyForces { branches: &active, controls: FORCE },
        || { polls += 1; polls == visits },
    ), Err(Error::Cancelled));
    assert_eq!(snapshot(&forced), before);
    assert_eq!(polls, visits, "late cancellation reaches the final atomic publication poll");
}
