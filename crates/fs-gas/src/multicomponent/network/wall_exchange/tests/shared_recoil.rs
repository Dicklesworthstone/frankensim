//! G1/G3 shared inertia and G4 original network transaction regressions.
use super::*;
use crate::multicomponent::Mesh1d;
use crate::multicomponent::network::dissipation::tests::law;
use crate::multicomponent::network::wall_exchange::CircularLaminarWall;
use crate::multicomponent::wall_exchange::{SharedRecoilControl, SlidingSupport};

fn support(mass_kg: f64, velocity_m_s: f64) -> SlidingSupport {
    SlidingSupport {
        mass_kg,
        velocity_m_s,
    }
}

fn time_control() -> SharedRecoilControl {
    SharedRecoilControl {
        absolute_velocity_tolerance_m_s: 1e-10,
        absolute_energy_tolerance_j: 1e-10,
        relative_tolerance: 1e-9,
        maximum_steps: 16384,
    }
}

fn shared(n: &MixtureNetwork1d, contact: WallContact) -> NetworkWallOperator<'static> {
    uniform_operator(n, contact)
        .with_shared_branch_supports(
            vec![Some(support(0.3, 0.5)); n.branches().len()],
            time_control(),
            1e-10,
            || false,
        )
        .unwrap()
}

fn adiabatic(rate_per_s: f64) -> WallContact {
    WallContact {
        drag: WallDrag::Linear { rate_per_s },
        thermal: WallThermal::Adiabatic,
    }
}

#[test]
fn one_support_uses_total_nonuniform_gas_mass_and_original_local_heat() {
    let g = gas(2, false, 1e12);
    let mesh = Mesh1d::new(vec![0.0, 0.25, 1.0], vec![0.4; 3], vec![0.1, 0.3]).unwrap();
    for (a, b) in [(1.7_f64, 0.0_f64), (0.0, 0.4), (1.7, 0.4)] {
        let branch = MixtureEuler1d::new_on_mesh(
            g.clone(),
            mesh.clone(),
            vec![
                g.conserved(&[0.8, 0.2], 2.0, 2.0).unwrap(),
                g.conserved(&[0.4, 1.6], 2.0, 4.0).unwrap(),
            ],
            Boundaries::Periodic,
            || false,
        )
        .unwrap();
        let mut n = MixtureNetwork1d::new(vec![branch], vec![], vec![], || false).unwrap();
        let before = fields(&n);
        let contact = WallContact {
            drag: WallDrag::LinearQuadratic {
                linear_rate_per_s: a,
                quadratic_per_m: b,
                wall_velocity_m_s: 999.0,
            },
            ..lumped()
        };
        let mut wall = shared(&n, contact);
        let dt = 0.2;
        let report = wall
            .advance_for(&mut n, Time::new(dt), controls(), || false)
            .unwrap();
        // Equal initial gas velocities and acceleration coefficients stay equal.
        // The two distinct gas amounts reduce to ONE mass 0.1+0.6 against M=0.3.
        let (m, mass, slip0) = (0.7_f64, 0.3_f64, 1.5_f64);
        let mu = m * mass / (m + mass);
        let tau = dt * m / mu;
        let slip = if a == 0.0 {
            slip0 / (1.0 + b * slip0 * tau)
        } else {
            let decay = (-a * tau).exp();
            slip0 * decay / (1.0 + b * slip0 * (1.0 - decay) / a)
        };
        let u = 2.0 - mu / m * (slip0 - slip);
        let common_u = 0.5 + mu / mass * (slip0 - slip);
        let heat = 0.5 * mu * (slip0 * slip0 - slip * slip);
        let mut gas_change = 0.0;
        let mut thermal_change = 0.0;
        for (i, gas_mass) in [0.1, 0.6].into_iter().enumerate() {
            let cell = n.branches()[0].cell(i).unwrap();
            let p = g.primitive(cell, || false).unwrap();
            let hot = 2.0 + (heat * gas_mass / m) / (2.5 * gas_mass);
            let q = dt * 0.3 * (5.0 - hot) / (1.0 + dt * 0.3 * (1.0 / (2.5 * gas_mass) + 1.0));
            close(p.velocity, u, 2e-7);
            close(p.temperature, hot + q / (2.5 * gas_mass), 3e-7);
            close(
                wall_temperature(wall.branch_contacts(0).unwrap()[i]),
                5.0 - q,
                2e-7,
            );
            assert_eq!(cell.partial_densities, before.0[0][i].partial_densities);
            gas_change += mesh.volumes()[i] * (cell.energy - before.0[0][i].energy);
            thermal_change += wall_temperature(wall.branch_contacts(0).unwrap()[i]) - 5.0;
        }
        let s = wall.branch_shared_support(0).unwrap();
        let exchange = wall.accepted_branch_shared_support_exchange(0).unwrap();
        close(s.velocity_m_s, common_u, 3e-7);
        close(exchange.drag_heating_j, heat, 3e-7);
        let kinetic = mass * (s.velocity_m_s - 0.5) * (0.5 * s.velocity_m_s + 0.25);
        close(
            exchange.support_impulse_ns,
            mass * (s.velocity_m_s - 0.5),
            1e-13,
        );
        close(exchange.kinetic_energy_change_j, kinetic, 1e-13);
        close(
            exchange.gas_impulse_ns + exchange.support_impulse_ns,
            0.0,
            1e-10,
        );
        close(exchange.gas_energy_change_j + kinetic, 0.0, 1e-10);
        close(gas_change + thermal_change + kinetic, 0.0, 2e-10);
        close(report.mechanical_work_to_gas_j + kinetic, 0.0, 1e-10);
        close(report.drag_heating_j, exchange.drag_heating_j, 1e-13);
        close(report.energy_residual_j, 0.0, 2e-10);
        assert_eq!(n.time(), Time::new(0.0));
    }
}

// Independent continuum two-rate star: s'=-A*s. The two real eigenmodes
// integrate both velocity and physical heat m*a*integral(s^2 dt) analytically.
fn linear_star(dt: f64) -> ([f64; 2], f64, [f64; 2]) {
    let masses = [2.0_f64, 5.0_f64];
    let rates = [0.7_f64, 2.3_f64];
    let support_mass = 3.0;
    let initial_u = [4.0, -1.2];
    let initial_support = -0.4;
    let s = [
        initial_u[0] - initial_support,
        initial_u[1] - initial_support,
    ];
    let a = [
        [
            rates[0] * (1.0 + masses[0] / support_mass),
            rates[1] * masses[1] / support_mass,
        ],
        [
            rates[0] * masses[0] / support_mass,
            rates[1] * (1.0 + masses[1] / support_mass),
        ],
    ];
    let gap = ((a[0][0] - a[1][1]).powi(2) + 4.0 * a[0][1] * a[1][0]).sqrt();
    let lambda = [
        0.5 * (a[0][0] + a[1][1] - gap),
        0.5 * (a[0][0] + a[1][1] + gap),
    ];
    let integral = |rate: f64| -(-rate * dt).exp_m1() / rate;
    let mut u = initial_u;
    let mut heat = [0.0; 2];
    for i in 0..2 {
        let derivative = -a[i][0] * s[0] - a[i][1] * s[1];
        let c0 = (derivative + lambda[1] * s[i]) / gap;
        let c1 = s[i] - c0;
        u[i] -= rates[i] * (c0 * integral(lambda[0]) + c1 * integral(lambda[1]));
        heat[i] = masses[i]
            * rates[i]
            * (c0 * c0 * integral(2.0 * lambda[0])
                + 2.0 * c0 * c1 * integral(lambda[0] + lambda[1])
                + c1 * c1 * integral(2.0 * lambda[1]));
    }
    let support = initial_support
        - (masses[0] * (u[0] - initial_u[0]) + masses[1] * (u[1] - initial_u[1])) / support_mass;
    (u, support, heat)
}

#[test]
fn heterogeneous_rates_follow_one_shared_continuum_and_each_cell_heat() {
    let g = gas(1, false, 0.0);
    let branch = MixtureEuler1d::new(
        g.clone(),
        Length::new(1.0),
        1.0,
        vec![
            g.conserved(&[2.0], 4.0, 4.0).unwrap(),
            g.conserved(&[5.0], -1.2, 10.0).unwrap(),
        ],
        Boundaries::Periodic,
        || false,
    )
    .unwrap();
    let mut n = MixtureNetwork1d::new(vec![branch], vec![], vec![], || false).unwrap();
    let before = fields(&n);
    let mut wall = NetworkWallOperator::new(
        &n,
        vec![vec![adiabatic(0.7), adiabatic(2.3)]],
        vec![],
        || false,
    )
    .unwrap()
    .with_shared_branch_supports(
        vec![Some(support(3.0, -0.4))],
        time_control(),
        1e-10,
        || false,
    )
    .unwrap();
    let dt = 0.9;
    let report = wall
        .advance_for(&mut n, Time::new(dt), controls(), || false)
        .unwrap();
    let (u, common_u, heat) = linear_star(dt);
    for (i, mass) in [2.0, 5.0].into_iter().enumerate() {
        let old = before.0[0][i].as_cell();
        let cell = n.branches()[0].cell(i).unwrap();
        let p = g.primitive(cell, || false).unwrap();
        let old_u = old.momentum / mass;
        close(p.velocity, u[i], 3e-7);
        close(
            cell.energy - old.energy - 0.5 * mass * (p.velocity - old_u) * (p.velocity + old_u),
            heat[i],
            2e-6,
        );
        assert_eq!(cell.partial_densities, old.partial_densities);
    }
    close(
        wall.branch_shared_support(0).unwrap().velocity_m_s,
        common_u,
        3e-7,
    );
    close(report.drag_heating_j, heat.into_iter().sum(), 2e-6);
    let exchange = wall.accepted_branch_shared_support_exchange(0).unwrap();
    close(exchange.energy_residual_j, 0.0, 2e-10);
    close(exchange.momentum_residual_ns, 0.0, 1e-10);
}

#[test]
fn shared_support_carries_both_source_halves_and_complete_endpoint_receipts() {
    let g = gas(1, false, 0.0);
    let n = network(&g, &[2.0], false);
    let before = fields(&n);
    let mut reference = network(&g, &[2.0], false);
    let mut once = shared(&reference, adiabatic(3.0));
    once.advance_for(&mut reference, Time::new(0.4), controls(), || false)
        .unwrap();
    let mut wall = shared(&n, adiabatic(3.0));
    let mut trial = network(&g, &[2.0], false);
    wall.begin_trial(&mut || false).unwrap();
    wall.advance_for_trial(&mut trial, Time::new(0.2), controls(), &mut || false)
        .unwrap();
    assert_eq!(wall.branch_shared_support(0), Some(support(0.3, 0.5)));
    assert!(wall.last_report().is_none());
    assert!(wall.accepted_branch_shared_support_exchange(0).is_none());
    assert_eq!(fields(&n), before);
    wall.advance_for_trial(&mut trial, Time::new(0.2), controls(), &mut || false)
        .unwrap();
    wall.publish_trial_report();
    for i in 0..2 {
        close(
            trial.branches()[0].cell(i).unwrap().momentum,
            reference.branches()[0].cell(i).unwrap().momentum,
            3e-7,
        );
        close(
            trial.branches()[0].cell(i).unwrap().energy,
            reference.branches()[0].cell(i).unwrap().energy,
            3e-7,
        );
    }
    let s = wall.branch_shared_support(0).unwrap();
    close(
        s.velocity_m_s,
        once.branch_shared_support(0).unwrap().velocity_m_s,
        3e-7,
    );
    let receipt = wall.accepted_branch_shared_support_exchange(0).unwrap();
    close(
        receipt.support_impulse_ns,
        s.mass_kg * (s.velocity_m_s - 0.5),
        1e-14,
    );
    close(
        receipt.kinetic_energy_change_j,
        receipt.support_impulse_ns * (0.5 * s.velocity_m_s + 0.25),
        1e-14,
    );
    close(receipt.energy_residual_j, 0.0, 1e-10);
    assert_eq!(wall.last_report().unwrap().elapsed, Time::new(0.4));
    assert_eq!(trial.time(), Time::new(0.0));
}

#[test]
fn shared_source_cancellation_preserves_prior_publication_and_exact_retry() {
    let g = gas(1, false, 0.0);
    let contact = WallContact {
        drag: WallDrag::Linear { rate_per_s: 3.0 },
        ..lumped()
    };
    let make = || {
        let mut n = network(&g, &[2.0], false);
        let mut wall = shared(&n, contact);
        wall.advance_for(&mut n, Time::new(0.02), controls(), || false)
            .unwrap();
        (n, wall)
    };
    let (mut expected, mut reference) = make();
    let mut polls = 0;
    reference
        .advance_for(&mut expected, Time::new(0.04), controls(), || {
            polls += 1;
            false
        })
        .unwrap();
    assert!(polls > 20);
    for stop in [1, polls / 4, polls / 2, polls - 1, polls] {
        let (mut n, mut wall) = make();
        let before = fields(&n);
        let contacts = wall.branch_contacts(0).unwrap().to_vec();
        let original_support = wall.branch_shared_support(0);
        let receipt = wall.accepted_branch_shared_support_exchange(0);
        let report = wall.last_report();
        let mut visits = 0;
        assert_eq!(
            wall.advance_for(&mut n, Time::new(0.04), controls(), || {
                visits += 1;
                visits == stop
            }),
            Err(WallExchangeError::Cancelled)
        );
        assert_eq!(fields(&n), before);
        assert_eq!(wall.branch_contacts(0).unwrap(), contacts);
        assert_eq!(wall.branch_shared_support(0), original_support);
        assert_eq!(wall.accepted_branch_shared_support_exchange(0), receipt);
        assert_eq!(wall.last_report(), report);
        wall.advance_for(&mut n, Time::new(0.04), controls(), || false)
            .unwrap();
        assert_eq!(fields(&n), fields(&expected));
        assert_eq!(wall.branch_contacts(0), reference.branch_contacts(0));
        assert_eq!(
            wall.branch_shared_support(0),
            reference.branch_shared_support(0)
        );
        assert_eq!(
            wall.accepted_branch_shared_support_exchange(0),
            reference.accepted_branch_shared_support_exchange(0)
        );
        assert_eq!(wall.last_report(), reference.last_report());
    }
    let mut n = network(&g, &[2.0], false);
    let before = fields(&n);
    // Unequal rates retain the adaptive path; a common linear rate now has
    // an exact collective solution even for a one-attempt budget.
    let mut limited = NetworkWallOperator::new(
        &n,
        vec![vec![
            contact,
            WallContact {
                drag: WallDrag::Linear { rate_per_s: 5.0 },
                ..contact
            },
        ]],
        vec![],
        || false,
    )
    .unwrap()
    .with_shared_branch_supports(
        vec![Some(support(0.3, 0.5))],
        SharedRecoilControl {
            maximum_steps: 1,
            ..time_control()
        },
        1e-10,
        || false,
    )
    .unwrap();
    assert!(
        limited
            .advance_for(&mut n, Time::new(0.4), controls(), || false)
            .is_err()
    );
    assert_eq!(fields(&n), before);
    assert_eq!(limited.branch_shared_support(0), Some(support(0.3, 0.5)));
    assert!(limited.last_report().is_none());
}

#[test]
fn independent_supports_on_other_branches_preserve_their_original_path() {
    let g = gas(1, false, 0.0);
    let contact = adiabatic(3.0);
    let descriptor = vec![vec![None; 2], vec![Some(support(0.1, -0.2)); 2]];
    let mut reference = network(&g, &[2.0, -1.0], false);
    let mut original = uniform_operator(&reference, contact)
        .with_finite_supports(descriptor.clone(), 1e-10, || false)
        .unwrap();
    let original_report = original
        .advance_for(&mut reference, Time::new(0.1), controls(), || false)
        .unwrap();
    let mut unchanged = network(&g, &[2.0, -1.0], false);
    let mut none = uniform_operator(&unchanged, contact)
        .with_shared_branch_supports(vec![None; 2], time_control(), 1e-10, || false)
        .unwrap()
        .with_finite_supports(descriptor.clone(), 1e-10, || false)
        .unwrap();
    let report = none
        .advance_for(&mut unchanged, Time::new(0.1), controls(), || false)
        .unwrap();
    assert_eq!(fields(&unchanged), fields(&reference));
    assert_eq!(report, original_report);
    for shared_first in [false, true] {
        let mut n = network(&g, &[2.0, -1.0], false);
        let op = uniform_operator(&n, contact);
        let mut wall = if shared_first {
            op.with_shared_branch_supports(
                vec![Some(support(0.3, 0.5)), None],
                time_control(),
                1e-10,
                || false,
            )
            .unwrap()
            .with_finite_supports(descriptor.clone(), 1e-10, || false)
            .unwrap()
        } else {
            op.with_finite_supports(descriptor.clone(), 1e-10, || false)
                .unwrap()
                .with_shared_branch_supports(
                    vec![Some(support(0.3, 0.5)), None],
                    time_control(),
                    1e-10,
                    || false,
                )
                .unwrap()
        };
        wall.advance_for(&mut n, Time::new(0.1), controls(), || false)
            .unwrap();
        assert_eq!(fields(&n).0[1], fields(&reference).0[1]);
        assert_eq!(wall.branch_supports(1), original.branch_supports(1));
        assert_eq!(
            wall.accepted_branch_support_exchanges(1),
            original.accepted_branch_support_exchanges(1)
        );
        assert_eq!(wall.branch_shared_support(1), None);
        assert!(wall.accepted_branch_shared_support_exchange(0).is_some());
    }
    let n = network(&g, &[2.0, -1.0], false);
    assert!(
        uniform_operator(&n, contact)
            .with_finite_supports(descriptor.clone(), 1e-10, || false,)
            .unwrap()
            .with_shared_branch_supports(
                vec![None, Some(support(0.3, 0.5))],
                time_control(),
                1e-10,
                || false,
            )
            .is_err()
    );
    assert!(
        uniform_operator(&n, contact)
            .with_shared_branch_supports(
                vec![None, Some(support(0.3, 0.5))],
                time_control(),
                1e-10,
                || false,
            )
            .unwrap()
            .with_finite_supports(descriptor, 1e-10, || false)
            .is_err()
    );
    for descriptor in [
        vec![],
        vec![Some(support(0.0, 0.0)); 2],
        vec![Some(support(0.3, f64::NAN)); 2],
    ] {
        assert!(
            uniform_operator(&n, contact)
                .with_shared_branch_supports(descriptor, time_control(), 1e-10, || false,)
                .is_err()
        );
    }
}

#[test]
fn shared_transport_admits_actual_slip_and_rejects_unsupported_final_slip() {
    let g = gas(1, false, 0.0);
    let geometry = CircularLaminarWall {
        diameter: Length::new((0.4 / core::f64::consts::FRAC_PI_4).sqrt()),
        nusselt_number: 0.0,
        reynolds_number: [0.0, 2300.0],
        maximum_mach_number: 0.3,
    };
    let mu = 0.004;
    let rate = 32.0 * mu / geometry.diameter.value().powi(2);
    for boost in [0.0, 20.0] {
        let mut n = network(&g, &[0.4 + boost], false);
        let mut wall = uniform_operator(&n, adiabatic(999.0))
            .with_shared_branch_supports(
                vec![Some(support(0.1, 0.2 + boost))],
                time_control(),
                1e-10,
                || false,
            )
            .unwrap()
            .with_circular_laminar_transport(
                law(&g, mu, 0.0, 0.0, 0.0),
                vec![Some(geometry)],
                || false,
            )
            .unwrap();
        wall.advance_for(&mut n, Time::new(0.2), controls(), || false)
            .unwrap();
        // Gas mass is 2*0.2. The true shared-support relative rate is 5*a.
        let removed = 0.2 * (1.0 - (-5.0 * rate * 0.2).exp());
        let common_u = 0.2 + boost + 0.8 * removed;
        close(
            wall.branch_shared_support(0).unwrap().velocity_m_s,
            common_u,
            3e-7,
        );
        for i in 0..2 {
            close(
                g.primitive(n.branches()[0].cell(i).unwrap(), || false)
                    .unwrap()
                    .velocity,
                0.4 + boost - 0.2 * removed,
                3e-7,
            );
            match wall.branch_contacts(0).unwrap()[i].drag {
                WallDrag::MovingLinear {
                    rate_per_s,
                    wall_velocity_m_s,
                } => {
                    close(rate_per_s, rate, 2e-13);
                    close(
                        wall_velocity_m_s,
                        wall.branch_shared_support(0).unwrap().velocity_m_s,
                        1e-13,
                    );
                }
                _ => panic!("lost source-owned drag or common support velocity"),
            }
        }
    }
    let mut n = network(&g, &[0.4], false);
    let before = fields(&n);
    let mut domain = geometry;
    domain.reynolds_number[0] = 20.0;
    let mut wall = uniform_operator(&n, adiabatic(0.0))
        .with_shared_branch_supports(
            vec![Some(support(0.01, 0.2))],
            time_control(),
            1e-10,
            || false,
        )
        .unwrap()
        .with_circular_laminar_transport(law(&g, mu, 0.0, 0.0, 0.0), vec![Some(domain)], || false)
        .unwrap();
    // Actual common recoil crosses Re=20, although using the old support speed
    // would leave the endpoint inside the original material correlation domain.
    assert!(matches!(
        wall.advance_for(&mut n, Time::new(0.2), controls(), || false),
        Err(WallExchangeError::TransportDomain {
            quantity: "Reynolds number",
            ..
        })
    ));
    assert_eq!(fields(&n), before);
    assert_eq!(wall.branch_shared_support(0), Some(support(0.01, 0.2)));
    assert!(wall.last_report().is_none());
    wall.advance_for(&mut n, Time::new(0.01), controls(), || false)
        .unwrap();
}

#[test]
fn stiff_common_rate_branch_advances_with_one_shared_mechanical_step() {
    let g = gas(1, false, 0.0);
    let branch = MixtureEuler1d::new(
        g.clone(),
        Length::new(0.5),
        0.4,
        vec![
            g.conserved(&[1.0], 0.6, 2.0).unwrap(),
            g.conserved(&[2.0], -0.2, 4.0).unwrap(),
        ],
        Boundaries::Periodic,
        || false,
    )
    .unwrap();
    let mut n = MixtureNetwork1d::new(vec![branch], vec![], vec![], || false).unwrap();
    let before = fields(&n);
    let initial_support = support(1e-6, -0.3);
    let rate = 3.0_f64;
    let dt = 1.0_f64;
    let mut wall = uniform_operator(&n, adiabatic(rate))
        .with_shared_branch_supports(
            vec![Some(initial_support)],
            SharedRecoilControl {
                maximum_steps: 1,
                ..time_control()
            },
            1e-10,
            || false,
        )
        .unwrap();
    let report = wall
        .advance_for(&mut n, Time::new(dt), controls(), || false)
        .unwrap();
    // Independent collective/deviation modes. The support mode has rate
    // 1,800,003/s while the two gas deviations retain the physical rate 3/s.
    let masses = [0.2_f64, 0.4_f64];
    let velocities = [0.6_f64, -0.2_f64];
    let gas_mass = masses.into_iter().sum::<f64>();
    let gas_mean = (masses[0] * velocities[0] + masses[1] * velocities[1]) / gas_mass;
    let mean_slip = gas_mean - initial_support.velocity_m_s;
    let collective_rate = rate * (1.0 + gas_mass / initial_support.mass_kg);
    let center = (gas_mass * gas_mean + initial_support.mass_kg * initial_support.velocity_m_s)
        / (gas_mass + initial_support.mass_kg);
    let collective_decay = (-collective_rate * dt).exp();
    let expected_support =
        center - gas_mass / (gas_mass + initial_support.mass_kg) * mean_slip * collective_decay;
    let exposure = |a: f64| -(-a * dt).exp_m1() / a;
    let mut heat = 0.0;
    let mut gas_change = 0.0;
    for i in 0..2 {
        let deviation = velocities[i] - gas_mean;
        let expected_velocity = center
            + initial_support.mass_kg / (gas_mass + initial_support.mass_kg)
                * mean_slip
                * collective_decay
            + deviation * (-rate * dt).exp();
        let expected_heat = masses[i]
            * rate
            * (deviation * deviation * exposure(2.0 * rate)
                + 2.0 * deviation * mean_slip * exposure(rate + collective_rate)
                + mean_slip * mean_slip * exposure(2.0 * collective_rate));
        let old = before.0[0][i].as_cell();
        let cell = n.branches()[0].cell(i).unwrap();
        let velocity = g.primitive(cell, || false).unwrap().velocity;
        let energy_change = 0.2 * (cell.energy - old.energy);
        let kinetic_change =
            0.5 * masses[i] * (velocity - velocities[i]) * (velocity + velocities[i]);
        close(velocity, expected_velocity, 2e-12);
        close(energy_change - kinetic_change, expected_heat, 3e-12);
        assert_eq!(cell.partial_densities, old.partial_densities);
        heat += expected_heat;
        gas_change += energy_change;
    }
    close(
        wall.branch_shared_support(0).unwrap().velocity_m_s,
        expected_support,
        2e-12,
    );
    let exchange = wall.accepted_branch_shared_support_exchange(0).unwrap();
    close(report.drag_heating_j, heat, 3e-12);
    close(gas_change + exchange.kinetic_energy_change_j, 0.0, 3e-12);
    close(exchange.momentum_residual_ns, 0.0, 3e-12);
    assert_eq!(n.time(), Time::new(0.0));
}
