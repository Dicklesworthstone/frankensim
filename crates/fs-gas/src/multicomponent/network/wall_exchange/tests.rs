//! G1 thermal relaxation, G0 source balances and G4 atomic wall publication.
use super::{NetworkWallOperator, NetworkWallReport};
use crate::multicomponent::network::dissipation::tests::{close, fields, gas};
use crate::multicomponent::network::{MixtureNetwork1d, Plenum};
use crate::multicomponent::wall_exchange::{
    WallContact, WallControl, WallDrag, WallExchangeError, WallThermal,
};
use crate::multicomponent::{Boundaries, MixtureEuler1d, MixtureGas, StepControl};
use fs_qty::{Length, Time};

// The original manufactured source has R_specific=1 and Cv=2.5 J/(kg K).
// Each branch cell has V=0.2 and rho=1; the plenum has V=0.1 and rho=2.
// Each physical gas volume therefore has heat capacity 0.5 J/K.
fn network(g: &MixtureGas, velocities: &[f64], with_plenum: bool) -> MixtureNetwork1d {
    let branches = velocities
        .iter()
        .map(|&velocity| {
            MixtureEuler1d::new(
                g.clone(),
                Length::new(0.5),
                0.4,
                vec![g.conserved(&[1.0], velocity, 2.0).unwrap(); 2],
                Boundaries::Periodic,
                || false,
            )
            .unwrap()
        })
        .collect();
    let plenums = if with_plenum {
        vec![Plenum::new(0.1, g.conserved(&[2.0], 0.0, 4.0).unwrap()).unwrap()]
    } else {
        Vec::new()
    };
    MixtureNetwork1d::new(branches, plenums, Vec::new(), || false).unwrap()
}
fn controls() -> WallControl {
    WallControl {
        absolute_energy_tolerance_j: 1e-12,
        relative_energy_tolerance: 1e-12,
        maximum_iterations: 80,
    }
}
fn lumped() -> WallContact {
    WallContact {
        drag: WallDrag::None,
        thermal: WallThermal::Lumped {
            temperature_k: 5.0,
            heat_capacity_j_per_k: 1.0,
            conductance_w_per_k: 0.3,
        },
    }
}
fn uniform_operator(n: &MixtureNetwork1d, contact: WallContact) -> NetworkWallOperator {
    NetworkWallOperator::new(
        n,
        n.branches()
            .iter()
            .map(|b| vec![contact; b.cell_count()])
            .collect(),
        vec![contact; n.plenums().len()],
        || false,
    )
    .unwrap()
}
fn wall_temperature(contact: WallContact) -> f64 {
    match contact.thermal {
        WallThermal::Lumped { temperature_k, .. }
        | WallThermal::Reservoir { temperature_k, .. } => temperature_k,
        WallThermal::Adiabatic => panic!("test expected thermal wall storage"),
    }
}
fn relaxation_error(steps: usize) -> f64 {
    let g = gas(1, false, 0.0);
    let mut n = network(&g, &[0.0], true);
    let mut wall = uniform_operator(&n, lumped());
    let duration = 0.4;
    for _ in 0..steps {
        wall.advance_for(
            &mut n,
            Time::new(duration / steps as f64),
            controls(),
            || false,
        )
        .unwrap();
    }
    // Independent continuum two-body solution: equilibrium=4 K and
    // relaxation rate=G*(1/Cgas+1/Cwall)=0.9/s.
    let expected = 4.0 - 2.0 * (-0.9_f64 * duration).exp();
    let actual = g
        .primitive(n.branches()[0].cell(0).unwrap(), || false)
        .unwrap()
        .temperature;
    (actual - expected).abs()
}

#[test]
fn finite_walls_match_two_body_heat_exchange_and_refine_toward_exponential() {
    let g = gas(1, false, 0.0);
    let mut n = network(&g, &[0.0], true);
    let before = fields(&n);
    let mut wall = uniform_operator(&n, lumped());
    let dt = 0.2;
    let report = wall
        .advance_for(&mut n, Time::new(dt), controls(), || false)
        .unwrap();
    // Independent backward-Euler two-capacity system:
    // Q=dt*G*(Tw0-Tg0)/(1+dt*G*(1/Cgas+1/Cwall)).
    let q = dt * 0.3 * (5.0 - 2.0) / (1.0 + dt * 0.3 * (2.0 + 1.0));
    let expected_gas = 2.0 + q / 0.5;
    let expected_wall = 5.0 - q;
    let mut independent_gas_change = 0.0;
    let mut independent_wall_change = 0.0;
    for i in 0..2 {
        let cell = n.branches()[0].cell(i).unwrap();
        close(g.primitive(cell, || false).unwrap().temperature, expected_gas, 2e-10);
        close(
            wall_temperature(wall.branch_contacts(0).unwrap()[i]),
            expected_wall,
            2e-10,
        );
        assert_eq!(cell.partial_densities, before.0[0][i].partial_densities);
        assert_eq!(cell.momentum, 0.0);
        independent_gas_change +=
            n.branches()[0].mesh().volumes()[i] * (cell.energy - before.0[0][i].energy);
        independent_wall_change += wall_temperature(wall.branch_contacts(0).unwrap()[i]) - 5.0;
    }
    let cell = n.plenums()[0].state();
    close(g.primitive(cell, || false).unwrap().temperature, expected_gas, 2e-10);
    assert_eq!(cell.partial_densities, before.1[0].partial_densities);
    assert_eq!(cell.momentum, 0.0);
    close(wall_temperature(wall.plenum_contacts()[0]), expected_wall, 2e-10);
    independent_gas_change += n.plenums()[0].volume_m3() * (cell.energy - before.1[0].energy);
    independent_wall_change += wall_temperature(wall.plenum_contacts()[0]) - 5.0;
    close(independent_gas_change + independent_wall_change, 0.0, 3e-11);
    close(independent_gas_change, 3.0 * q, 3e-11);
    close(report.heat_to_gas_j, independent_gas_change, 3e-11);
    close(report.wall_energy_change_j, independent_wall_change, 3e-11);
    assert_eq!(report.external_heat_to_gas_j, 0.0);
    assert_eq!(report.drag_heating_j, 0.0);
    close(report.energy_residual_j, 0.0, 3e-11);
    assert_eq!(report.elapsed, Time::new(dt));
    assert_eq!(wall.last_report(), Some(report));
    assert_eq!(n.time(), Time::new(0.0));
    let coarse = relaxation_error(10);
    let fine = relaxation_error(20);
    assert!(fine < 0.005, "fine thermal error={fine:e}");
    assert!(fine < 0.56 * coarse, "coarse={coarse:e}, fine={fine:e}");
}

#[test]
fn linear_and_darcy_wall_drag_convert_kinetic_energy_once_and_keep_branch_impulses() {
    let dt = 0.2_f64;
    let initial_speed = 2.0_f64;
    for (drag, expected_speed) in [
        (
            WallDrag::Linear { rate_per_s: 3.0 },
            initial_speed * (-3.0 * dt).exp(),
        ),
        (
            WallDrag::Darcy {
                friction_factor: 0.1,
                hydraulic_diameter_m: 0.2,
            },
            initial_speed / (1.0 + 0.1 * initial_speed * dt / (2.0 * 0.2)),
        ),
    ] {
        for formation in [0.0, 1e12] {
            let g = gas(1, false, formation);
            let mut n = network(&g, &[initial_speed, -initial_speed], false);
            let before = fields(&n);
            let mut wall = uniform_operator(&n, WallContact {
                drag,
                thermal: WallThermal::Adiabatic,
            });
            let report = wall
                .advance_for(&mut n, Time::new(dt), controls(), || false)
                .unwrap();
            let mut kinetic_loss = 0.0;
            for branch in 0..2 {
                let sign = if branch == 0 { 1.0 } else { -1.0 };
                let mut impulse = 0.0;
                for i in 0..2 {
                    let cell = n.branches()[branch].cell(i).unwrap();
                    let old = before.0[branch][i].as_cell();
                    let primitive = g.primitive(cell, || false).unwrap();
                    let volume = n.branches()[branch].mesh().volumes()[i];
                    close(primitive.velocity, sign * expected_speed, 2e-12);
                    close(
                        primitive.temperature,
                        2.0 + (initial_speed * initial_speed - expected_speed * expected_speed) / 5.0,
                        2e-10,
                    );
                    assert_eq!(cell.partial_densities, old.partial_densities);
                    assert_eq!(cell.energy, old.energy, "drag added total energy twice");
                    impulse += volume * (old.momentum - cell.momentum);
                    kinetic_loss += volume * 0.5
                        * (old.momentum * old.momentum - cell.momentum * cell.momentum)
                        / cell.density();
                }
                close(wall.accepted_branch_impulses().unwrap()[branch], impulse, 3e-12);
                close(impulse, sign * 0.4 * (initial_speed - expected_speed), 3e-12);
            }
            assert!(wall.accepted_branch_impulses().unwrap()[0] > 0.0);
            assert!(wall.accepted_branch_impulses().unwrap()[1] < 0.0);
            close(report.drag_heating_j, kinetic_loss, 3e-12);
            assert!(report.drag_heating_j > 0.0);
            assert_eq!(report.heat_to_gas_j, 0.0);
            assert_eq!(report.wall_energy_change_j, 0.0);
            assert_eq!(report.external_heat_to_gas_j, 0.0);
            close(report.energy_residual_j, 0.0, 3e-12);
            assert_eq!(n.time(), Time::new(0.0));
        }
    }
}

#[test]
fn reservoir_heat_closes_external_balance_and_inadmissible_heat_never_publishes() {
    let g = gas(1, false, 0.0);
    let mut n = network(&g, &[0.0], true);
    let before = fields(&n);
    let contact = WallContact {
        drag: WallDrag::None,
        thermal: WallThermal::Reservoir {
            temperature_k: 8.0,
            conductance_w_per_k: 0.4,
        },
    };
    let mut wall = uniform_operator(&n, contact);
    let report = wall
        .advance_for(&mut n, Time::new(0.25), controls(), || false)
        .unwrap();
    let after = fields(&n);
    let mut energy = 0.0;
    for i in 0..2 {
        let cell = n.branches()[0].cell(i).unwrap();
        close(g.primitive(cell, || false).unwrap().temperature, 3.0, 2e-10);
        energy += n.branches()[0].mesh().volumes()[i]
            * (after.0[0][i].energy - before.0[0][i].energy);
        assert_eq!(wall.branch_contacts(0).unwrap()[i], contact);
    }
    close(g.primitive(n.plenums()[0].state(), || false).unwrap().temperature, 3.0, 2e-10);
    energy += n.plenums()[0].volume_m3() * (after.1[0].energy - before.1[0].energy);
    close(energy, 1.5, 3e-11);
    close(report.heat_to_gas_j, energy, 3e-11);
    close(report.external_heat_to_gas_j, energy, 3e-11);
    assert_eq!(report.wall_energy_change_j, 0.0);
    close(report.energy_residual_j, 0.0, 3e-11);
    let invalid = WallContact {
        thermal: WallThermal::Reservoir {
            temperature_k: 8.0,
            conductance_w_per_k: -0.4,
        },
        ..contact
    };
    assert!(matches!(
        NetworkWallOperator::new(&n, vec![vec![invalid; 2]], vec![contact], || false),
        Err(WallExchangeError::Invalid(_))
    ));
    let hot = WallContact {
        thermal: WallThermal::Reservoir {
            temperature_k: 1e6,
            conductance_w_per_k: 1.0,
        },
        ..contact
    };
    let mut hot_wall = uniform_operator(&n, hot);
    let unchanged = fields(&n);
    let clock = n.time();
    assert!(hot_wall.advance_for(&mut n, Time::new(1.0), controls(), || false).is_err());
    assert_eq!(fields(&n), unchanged);
    assert_eq!(n.time(), clock);
    assert!(hot_wall.last_report().is_none());
    assert_eq!(hot_wall.branch_contacts(0).unwrap(), [hot, hot]);
    assert_eq!(hot_wall.plenum_contacts(), [hot]);
}

#[test]
fn late_wall_cancellation_preserves_storage_temperatures_and_reports_then_reuses_work() {
    let g = gas(1, false, 0.0);
    let make = || {
        let mut n = network(&g, &[0.0], true);
        let mut wall = uniform_operator(&n, lumped());
        wall.advance_for(&mut n, Time::new(0.02), controls(), || false).unwrap();
        n.step(Time::new(0.001), StepControl::default(), || false).unwrap();
        (n, wall)
    };
    let (mut reference, mut reference_wall) = make();
    let mut polls = 0;
    reference_wall.advance_for(&mut reference, Time::new(0.03), controls(), || {
        polls += 1;
        false
    }).unwrap();
    assert!(polls > 1);
    let (mut n, mut wall) = make();
    let before = fields(&n);
    let clock = n.time();
    let network_report = n.last_report();
    let ports = n.accepted_port_exchanges().unwrap().to_vec();
    let branches = n.accepted_branch_reports().unwrap().to_vec();
    let wall_report: Option<NetworkWallReport> = wall.last_report();
    let branch_contacts = wall.branch_contacts(0).unwrap().to_vec();
    let plenum_contacts = wall.plenum_contacts().to_vec();
    let impulses = wall.accepted_branch_impulses().unwrap().to_vec();
    let mut visits = 0;
    let result = wall.advance_for(&mut n, Time::new(0.03), controls(), || {
        visits += 1;
        visits >= polls
    });
    assert_eq!(result, Err(WallExchangeError::Cancelled));
    assert_eq!(fields(&n), before);
    assert_eq!(n.time(), clock);
    assert_eq!(n.last_report(), network_report);
    assert_eq!(n.accepted_port_exchanges().unwrap(), ports);
    assert_eq!(n.accepted_branch_reports().unwrap(), branches);
    assert_eq!(wall.last_report(), wall_report);
    assert_eq!(wall.branch_contacts(0).unwrap(), branch_contacts);
    assert_eq!(wall.plenum_contacts(), plenum_contacts);
    assert_eq!(wall.accepted_branch_impulses().unwrap(), impulses);
    wall.advance_for(&mut n, Time::new(0.03), controls(), || false).unwrap();
    let current = fields(&n);
    let expected = fields(&reference);
    for (a, b) in current.0.iter().flatten().chain(&current.1)
        .zip(expected.0.iter().flatten().chain(&expected.1))
    {
        assert_eq!(a.partial_densities, b.partial_densities);
        close(a.momentum, b.momentum, 3e-12);
        close(a.energy, b.energy, 3e-11);
    }
    for (a, b) in wall.branch_contacts(0).unwrap().iter().chain(wall.plenum_contacts())
        .zip(reference_wall.branch_contacts(0).unwrap().iter().chain(reference_wall.plenum_contacts()))
    {
        close(wall_temperature(*a), wall_temperature(*b), 3e-11);
    }
    assert_eq!(n.time(), reference.time());
    assert!(n.last_report().is_none(), "stale advection report survived source mutation");
}
