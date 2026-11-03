//! Finite support momentum, heat composition and original transaction checks.
use super::*;
use crate::multicomponent::network::dissipation::tests::law;
use crate::multicomponent::network::wall_exchange::CircularLaminarWall;
use crate::multicomponent::wall_exchange::SlidingSupport;

fn support(mass: f64, velocity: f64) -> SlidingSupport {
    SlidingSupport {
        mass_kg: mass,
        velocity_m_s: velocity,
    }
}

fn finite(
    n: &MixtureNetwork1d,
    contact: WallContact,
    mass: f64,
    velocity: f64,
) -> NetworkWallOperator<'static> {
    uniform_operator(n, contact)
        .with_finite_supports(
            n.branches()
                .iter()
                .map(|b| vec![Some(support(mass, velocity)); b.cell_count()])
                .collect(),
            1e-11,
            || false,
        )
        .unwrap()
}

fn contact() -> WallContact {
    WallContact {
        drag: WallDrag::Linear { rate_per_s: 3.0 },
        ..lumped()
    }
}

#[test]
fn finite_support_recoil_and_original_heat_close_actual_endpoint_inventories() {
    let g = gas(1, false, 1e12);
    for drag in [
        WallDrag::Linear { rate_per_s: 3.0 },
        WallDrag::Darcy {
            friction_factor: 0.4,
            hydraulic_diameter_m: 0.2,
        },
        WallDrag::LinearQuadratic {
            linear_rate_per_s: 3.0,
            quadratic_per_m: 1.0,
            wall_velocity_m_s: 999.0,
        },
    ] {
        let mut n = network(&g, &[2.0], false);
        let before = fields(&n);
        let mut wall = finite(&n, WallContact { drag, ..lumped() }, 0.1, 0.5);
        let dt = 0.2;
        let report = wall
            .advance_for(&mut n, Time::new(dt), controls(), || false)
            .unwrap();
        let (a, b) = match drag {
            WallDrag::Linear { .. } => (3.0_f64, 0.0),
            WallDrag::Darcy { .. } => (0.0, 1.0),
            _ => (3.0, 1.0),
        };
        // Independent exact relative law: mu=m*M/(m+M), ds/dt=-(m/mu)(a*s+b*s²).
        let (m, mass, slip0) = (0.2_f64, 0.1_f64, 1.5_f64);
        let mu = m * mass / (m + mass);
        let tau = dt * m / mu;
        let slip = if a == 0.0 {
            slip0 / (1.0 + b * slip0 * tau)
        } else {
            let decay = (-a * tau).exp();
            slip0 * decay / (1.0 + b * slip0 * (1.0 - decay) / a)
        };
        let u = 2.0 - mu / m * (slip0 - slip);
        let wall_u = 0.5 + mu / mass * (slip0 - slip);
        let heat = 0.5 * mu * (slip0 * slip0 - slip * slip);
        let hot_gas = 2.0 + heat / (m * 2.5);
        // Original two-capacity backward Euler contact AFTER mechanical heating.
        let q = dt * 0.3 * (5.0 - hot_gas) / (1.0 + dt * 0.3 * (1.0 / (m * 2.5) + 1.0));
        let mut gas_change = 0.0;
        let mut support_change = 0.0;
        let mut thermal_change = 0.0;
        for i in 0..2 {
            let cell = n.branches()[0].cell(i).unwrap();
            let p = g.primitive(cell, || false).unwrap();
            let actual_support = wall.branch_supports(0).unwrap()[i].unwrap();
            let exchange = wall.accepted_branch_support_exchanges(0).unwrap()[i].unwrap();
            close(p.velocity, u, 2e-12);
            close(actual_support.velocity_m_s, wall_u, 2e-12);
            close(p.temperature, hot_gas + q / (m * 2.5), 3e-10);
            close(
                wall_temperature(wall.branch_contacts(0).unwrap()[i]),
                5.0 - q,
                3e-10,
            );
            assert_eq!(cell.partial_densities, before.0[0][i].partial_densities);
            let impulse = mass * (actual_support.velocity_m_s - 0.5);
            let kinetic = impulse * (0.5 * actual_support.velocity_m_s + 0.25);
            close(exchange.support_impulse_ns, impulse, 1e-13);
            close(exchange.kinetic_energy_change_j, kinetic, 1e-13);
            close(exchange.gas_impulse_ns + impulse, 0.0, 2e-13);
            close(exchange.gas_energy_change_j + kinetic, 0.0, 2e-12);
            close(exchange.drag_heating_j, heat, 2e-12);
            gas_change += 0.2 * (cell.energy - before.0[0][i].energy);
            support_change += kinetic;
            thermal_change += wall_temperature(wall.branch_contacts(0).unwrap()[i]) - 5.0;
        }
        close(gas_change + support_change + thermal_change, 0.0, 4e-11);
        close(report.mechanical_work_to_gas_j + support_change, 0.0, 4e-12);
        close(report.drag_heating_j, 2.0 * heat, 4e-12);
        close(report.heat_to_gas_j, 2.0 * q, 4e-11);
        close(report.energy_residual_j, 0.0, 4e-11);
        assert_eq!(n.time(), Time::new(0.0));
    }
}

#[test]
fn both_source_halves_retain_one_support_inventory_and_publish_together() {
    let g = gas(1, false, 0.0);
    let n = network(&g, &[2.0], false);
    let before = fields(&n);
    let adiabatic = WallContact {
        thermal: WallThermal::Adiabatic,
        ..contact()
    };
    let mut reference = network(&g, &[2.0], false);
    let mut once = finite(&reference, adiabatic, 0.1, 0.5);
    once.advance_for(&mut reference, Time::new(0.4), controls(), || false)
        .unwrap();

    let mut wall = finite(&n, adiabatic, 0.1, 0.5);
    let mut trial = network(&g, &[2.0], false);
    wall.begin_trial(&mut || false).unwrap();
    wall.advance_for_trial(&mut trial, Time::new(0.2), controls(), &mut || false)
        .unwrap();
    assert!(wall.last_report().is_none());
    assert!(wall.accepted_branch_support_exchanges(0).is_none());
    assert_eq!(
        wall.branch_supports(0).unwrap(),
        &[Some(support(0.1, 0.5)); 2]
    );
    assert_eq!(fields(&n), before);
    wall.advance_for_trial(&mut trial, Time::new(0.2), controls(), &mut || false)
        .unwrap();
    wall.publish_trial_report();
    for i in 0..2 {
        close(
            trial.branches()[0].cell(i).unwrap().momentum,
            reference.branches()[0].cell(i).unwrap().momentum,
            3e-13,
        );
        close(
            trial.branches()[0].cell(i).unwrap().energy,
            reference.branches()[0].cell(i).unwrap().energy,
            3e-12,
        );
        close(
            wall.branch_supports(0).unwrap()[i].unwrap().velocity_m_s,
            once.branch_supports(0).unwrap()[i].unwrap().velocity_m_s,
            3e-13,
        );
        let s = wall.branch_supports(0).unwrap()[i].unwrap();
        let receipt = wall.accepted_branch_support_exchanges(0).unwrap()[i].unwrap();
        close(
            receipt.kinetic_energy_change_j,
            0.1 * (s.velocity_m_s - 0.5) * (0.5 * s.velocity_m_s + 0.25),
            1e-14,
        );
        close(receipt.energy_residual_j, 0.0, 3e-12);
    }
    assert_eq!(wall.last_report().unwrap().elapsed, Time::new(0.4));
    assert_eq!(trial.time(), Time::new(0.0));
}

#[test]
fn every_source_cancellation_preserves_supports_gas_heat_and_replay() {
    let g = gas(1, false, 0.0);
    let mut expected = network(&g, &[2.0], false);
    let mut reference = finite(&expected, contact(), 0.1, 0.5);
    let mut polls = 0;
    reference
        .advance_for(&mut expected, Time::new(0.2), controls(), || {
            polls += 1;
            false
        })
        .unwrap();
    assert!(polls > 10);
    for stop in 1..=polls {
        let mut n = network(&g, &[2.0], false);
        let before = fields(&n);
        let mut wall = finite(&n, contact(), 0.1, 0.5);
        let original_contacts = wall.branch_contacts(0).unwrap().to_vec();
        let mut count = 0;
        assert_eq!(
            wall.advance_for(&mut n, Time::new(0.2), controls(), || {
                count += 1;
                count == stop
            }),
            Err(WallExchangeError::Cancelled)
        );
        assert_eq!(fields(&n), before);
        assert_eq!(wall.branch_contacts(0).unwrap(), original_contacts);
        assert_eq!(
            wall.branch_supports(0).unwrap(),
            &[Some(support(0.1, 0.5)); 2]
        );
        assert!(wall.last_report().is_none());
        wall.advance_for(&mut n, Time::new(0.2), controls(), || false)
            .unwrap();
        assert_eq!(fields(&n), fields(&expected));
        assert_eq!(wall.branch_supports(0), reference.branch_supports(0));
        assert_eq!(wall.branch_contacts(0), reference.branch_contacts(0));
        assert_eq!(
            wall.accepted_branch_support_exchanges(0),
            reference.accepted_branch_support_exchanges(0)
        );
        assert_eq!(wall.last_report(), reference.last_report());
    }
}

fn geometry() -> CircularLaminarWall {
    CircularLaminarWall {
        diameter: Length::new((0.4 / core::f64::consts::FRAC_PI_4).sqrt()),
        nusselt_number: 0.0,
        reynolds_number: [0.0, 2300.0],
        maximum_mach_number: 0.3,
    }
}

#[test]
fn original_transport_uses_recoil_relative_slip_at_both_endpoints() {
    let g = gas(1, false, 0.0);
    let mu = 0.004;
    let rate = 32.0 * mu / geometry().diameter.value().powi(2);
    for boost in [0.0, 20.0] {
        let mut n = network(&g, &[0.4 + boost], false);
        let mut wall = finite(
            &n,
            WallContact {
                drag: WallDrag::Linear { rate_per_s: 999.0 },
                thermal: WallThermal::Adiabatic,
            },
            0.1,
            0.2 + boost,
        )
        .with_circular_laminar_transport(law(&g, mu, 0.0, 0.0, 0.0), vec![Some(geometry())], || {
            false
        })
        .unwrap();
        wall.advance_for(&mut n, Time::new(0.2), controls(), || false)
            .unwrap();
        let removed = 0.2 * (1.0 - (-3.0 * rate * 0.2).exp());
        let wall_u = 0.2 + boost + (2.0 / 3.0) * removed;
        for i in 0..2 {
            close(
                g.primitive(n.branches()[0].cell(i).unwrap(), || false)
                    .unwrap()
                    .velocity,
                0.4 + boost - removed / 3.0,
                3e-12,
            );
            close(
                wall.branch_supports(0).unwrap()[i].unwrap().velocity_m_s,
                wall_u,
                3e-12,
            );
            match wall.branch_contacts(0).unwrap()[i].drag {
                WallDrag::MovingLinear {
                    rate_per_s,
                    wall_velocity_m_s,
                } => {
                    close(rate_per_s, rate, 2e-13);
                    close(wall_velocity_m_s, wall_u, 2e-13);
                }
                _ => panic!("lost original material rate or actual support velocity"),
            }
        }
    }

    let mut n = network(&g, &[0.4], false);
    let before = fields(&n);
    let mut domain = geometry();
    domain.reynolds_number[0] = 20.0;
    let mut wall = finite(
        &n,
        WallContact {
            drag: WallDrag::None,
            thermal: WallThermal::Adiabatic,
        },
        0.01,
        0.2,
    )
    .with_circular_laminar_transport(law(&g, mu, 0.0, 0.0, 0.0), vec![Some(domain)], || false)
    .unwrap();
    // Initial Re is ~35.7. Recoil reduces the actual endpoint slip to Re~12.4,
    // outside the authored domain. Keeping the OLD support velocity would
    // incorrectly admit Re~34.6 at this endpoint.
    assert!(matches!(
        wall.advance_for(&mut n, Time::new(0.2), controls(), || false),
        Err(WallExchangeError::TransportDomain {
            quantity: "Reynolds number",
            ..
        })
    ));
    assert_eq!(fields(&n), before);
    assert_eq!(
        wall.branch_supports(0).unwrap(),
        &[Some(support(0.01, 0.2)); 2]
    );
    wall.advance_for(&mut n, Time::new(0.01), controls(), || false)
        .unwrap();
}

#[test]
fn unselected_cells_keep_original_bits_and_invalid_supports_refuse() {
    let g = gas(1, false, 0.0);
    let mut reference = network(&g, &[2.0], false);
    let mut original = uniform_operator(&reference, contact());
    let expected = original
        .advance_for(&mut reference, Time::new(0.2), controls(), || false)
        .unwrap();
    let mut n = network(&g, &[2.0], false);
    let mut wall = uniform_operator(&n, contact())
        .with_finite_supports(vec![vec![None; 2]], 1e-11, || false)
        .unwrap();
    let report = wall
        .advance_for(&mut n, Time::new(0.2), controls(), || false)
        .unwrap();
    assert_eq!(fields(&n), fields(&reference));
    assert_eq!(report, expected);
    assert_eq!(wall.branch_contacts(0), original.branch_contacts(0));
    assert_eq!(
        wall.accepted_branch_support_exchanges(0).unwrap(),
        &[None; 2]
    );
    for input in [
        vec![],
        vec![vec![Some(support(0.1, 0.0))]],
        vec![vec![Some(support(0.0, 0.0)); 2]],
        vec![vec![Some(support(0.1, f64::NAN)); 2]],
    ] {
        assert!(matches!(
            uniform_operator(&n, contact()).with_finite_supports(input, 1e-11, || false),
            Err(WallExchangeError::Invalid(_))
        ));
    }
}

// Full original transport/shared-storage integration; excluded only from the
// standalone native source harness, never from the crate's normal test suite.
#[test]
fn complete_molecular_wave_step_carries_recoil_on_the_same_accepted_clock() {
    use crate::multicomponent::network::dissipating::{
        DissipatingNetwork1d, DissipatingStepControl,
    };
    let g = gas(1, false, 0.0);
    let n = network(&g, &[0.4], false);
    let mut wall = finite(
        &n,
        WallContact {
            drag: WallDrag::Linear { rate_per_s: 3.0 },
            thermal: WallThermal::Adiabatic,
        },
        0.1,
        0.2,
    );
    let mut flow =
        DissipatingNetwork1d::new(n, law(&g, 0.004, 0.0, 0.02, 0.0), vec![], || false).unwrap();
    let dt = 0.01;
    let report = flow
        .step_with_wall(
            &mut wall,
            Time::new(dt),
            DissipatingStepControl::default(),
            controls(),
            || false,
        )
        .unwrap();
    let halves = report.wall.unwrap();
    assert_eq!(flow.time(), report.elapsed);
    assert_eq!(wall.last_report().unwrap().elapsed, report.elapsed);
    assert_eq!(
        halves.first.elapsed.value() + halves.second.elapsed.value(),
        report.elapsed.value()
    );
    let s = 0.2 * (-9.0 * report.elapsed.value()).exp();
    let u = 0.4 - (0.2 - s) / 3.0;
    let wall_u = 0.2 + (2.0 / 3.0) * (0.2 - s);
    for i in 0..2 {
        close(
            flow.gas()
                .primitive(flow.branches()[0].cell(i).unwrap(), || false)
                .unwrap()
                .velocity,
            u,
            3e-12,
        );
        close(
            wall.branch_supports(0).unwrap()[i].unwrap().velocity_m_s,
            wall_u,
            3e-12,
        );
    }
    let support_change = wall
        .accepted_branch_support_exchanges(0)
        .unwrap()
        .iter()
        .map(|r| r.unwrap().kinetic_energy_change_j)
        .sum::<f64>();
    close(
        report.after.energy - report.before.energy + support_change,
        0.0,
        3e-11,
    );
    close(report.absolute_energy_residual_j, 0.0, 3e-11);
}

#[test]
fn independent_supports_retain_one_shared_thermal_inventory() {
    use crate::multicomponent::network::wall_exchange::{
        SharedWall, SharedWallControl, SharedWallGas, SharedWallPort, SharedWallStorage,
    };
    let g = gas(1, false, 0.0);
    let mut n = network(&g, &[2.0], false);
    let before = fields(&n);
    let mut wall = finite(
        &n,
        WallContact {
            drag: WallDrag::Linear { rate_per_s: 3.0 },
            thermal: WallThermal::Adiabatic,
        },
        0.1,
        0.5,
    )
    .with_shared_walls(
        &n,
        vec![SharedWall {
            storage: SharedWallStorage {
                thermal: WallThermal::Lumped {
                    temperature_k: 5.0,
                    heat_capacity_j_per_k: 1.0,
                    conductance_w_per_k: 0.0,
                },
                heat_storage: None,
            },
            ports: (0..2)
                .map(|cell| SharedWallPort {
                    gas: SharedWallGas::Branch { branch: 0, cell },
                    conductance_w_per_k: 0.3,
                })
                .collect(),
        }],
        SharedWallControl::default(),
        || false,
    )
    .unwrap();
    let report = wall
        .advance_for(&mut n, Time::new(0.2), controls(), || false)
        .unwrap();
    let gas_change = (0..2)
        .map(|i| 0.2 * (n.branches()[0].cell(i).unwrap().energy - before.0[0][i].energy))
        .sum::<f64>();
    let support_change = wall
        .accepted_branch_support_exchanges(0)
        .unwrap()
        .iter()
        .map(|r| r.unwrap().kinetic_energy_change_j)
        .sum::<f64>();
    let WallThermal::Lumped { temperature_k, .. } = wall.shared_wall(0).unwrap().storage.thermal
    else {
        panic!("lost shared thermal inventory")
    };
    close(
        gas_change + support_change + temperature_k - 5.0,
        0.0,
        3e-11,
    );
    close(report.wall_energy_change_j, temperature_k - 5.0, 3e-11);
    close(report.energy_residual_j, 0.0, 3e-11);
}
