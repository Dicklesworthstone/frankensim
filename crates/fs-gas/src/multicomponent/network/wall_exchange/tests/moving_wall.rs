//! Prescribed wall drives on the actual finite-volume mixture fields.

use super::*;

fn drag(darcy: bool, wall: f64) -> WallDrag {
    if darcy {
        WallDrag::MovingDarcy {
            friction_factor: 0.6,
            hydraulic_diameter_m: 0.3,
            wall_velocity_m_s: wall,
        }
    } else {
        WallDrag::MovingLinear {
            rate_per_s: 2.0,
            wall_velocity_m_s: wall,
        }
    }
}

#[test]
fn moving_wall_network_entrainment_braking_and_reactions_close_actual_energy() {
    let g = gas(1, false, 0.0);
    let velocities = [0.0, 0.7, -0.4];
    let speeds = [0.8, 0.2, -0.7];
    for resistance in 0..3 {
        let mut n = network(&g, &velocities, false);
        let before = fields(&n);
        let mut wall = NetworkWallOperator::new(
            &n,
            speeds
                .iter()
                .map(|&speed| {
                    vec![
                        WallContact {
                            drag: if resistance == 2 {
                                WallDrag::LinearQuadratic {
                                    linear_rate_per_s: 2.0,
                                    quadratic_per_m: 1.0,
                                    wall_velocity_m_s: speed,
                                }
                            } else {
                                drag(resistance == 1, speed)
                            },
                            thermal: WallThermal::Adiabatic,
                        };
                        2
                    ]
                })
                .collect(),
            vec![],
            || false,
        )
        .unwrap();
        let dt = 0.2;
        let report = wall
            .advance_for(&mut n, Time::new(dt), controls(), || false)
            .unwrap();
        let mut work = 0.0;
        let mut heating = 0.0;
        for (b, (&u, &speed)) in velocities.iter().zip(&speeds).enumerate() {
            let slip: f64 = u - speed;
            let final_slip = if resistance == 2 {
                slip / ((2.0 * dt).exp() + 0.5 * slip.abs() * (2.0 * dt).exp_m1())
            } else if resistance == 1 {
                slip / (1.0 + slip.abs() * dt)
            } else {
                slip * (-2.0 * dt).exp()
            };
            let final_velocity = speed + final_slip;
            let temperature = 2.0 + (slip * slip - final_slip * final_slip) / 5.0;
            let mut impulse = 0.0;
            for i in 0..2 {
                let cell = n.branches()[b].cell(i).unwrap();
                let old = &before.0[b][i];
                let primitive = g.primitive(cell, || false).unwrap();
                let volume = n.branches()[b].mesh().volumes()[i];
                close(primitive.velocity, final_velocity, 3e-14);
                close(primitive.temperature, temperature, 3e-11);
                assert_eq!(cell.partial_densities, old.partial_densities);
                let gas_impulse = (cell.momentum - old.momentum) * volume;
                close(
                    (cell.energy - old.energy) * volume,
                    speed * gas_impulse,
                    3e-14,
                );
                impulse -= gas_impulse;
                work += speed * gas_impulse;
                heating += 0.5 * volume * (slip * slip - final_slip * final_slip);
            }
            close(wall.accepted_branch_impulses().unwrap()[b], impulse, 3e-14);
        }
        close(report.mechanical_work_to_gas_j, work, 3e-14);
        close(report.drag_heating_j, heating, 3e-14);
        assert!(report.drag_heating_j > 0.0);
        assert_eq!(report.heat_to_gas_j, 0.0);
        assert_eq!(report.wall_energy_change_j, 0.0);
        assert_eq!(report.external_heat_to_gas_j, 0.0);
        close(report.energy_residual_j, 0.0, 3e-14);
        assert_eq!(n.time(), Time::new(0.0));
    }
}

#[test]
fn moving_wall_late_cancellation_and_hot_source_refusal_publish_nothing() {
    let g = gas(1, false, 0.0);
    let make = || {
        let mut n = network(&g, &[0.2], false);
        let mut wall = uniform_operator(
            &n,
            WallContact {
                drag: drag(false, 0.8),
                ..lumped()
            },
        );
        wall.advance_for(&mut n, Time::new(0.02), controls(), || false)
            .unwrap();
        (n, wall)
    };
    let (mut reference, mut reference_wall) = make();
    let mut polls = 0;
    reference_wall
        .advance_for(&mut reference, Time::new(0.03), controls(), || {
            polls += 1;
            false
        })
        .unwrap();
    let (mut n, mut wall) = make();
    let before = fields(&n);
    let contacts = wall.branch_contacts(0).unwrap().to_vec();
    let accepted = wall.last_report();
    let impulses = wall.accepted_branch_impulses().unwrap().to_vec();
    let mut visits = 0;
    assert_eq!(
        wall.advance_for(&mut n, Time::new(0.03), controls(), || {
            visits += 1;
            visits >= polls
        }),
        Err(WallExchangeError::Cancelled)
    );
    assert_eq!(fields(&n), before);
    assert_eq!(wall.branch_contacts(0).unwrap(), contacts);
    assert_eq!(wall.accepted_branch_impulses().unwrap(), impulses);
    assert_eq!(wall.last_report(), accepted);
    wall.advance_for(&mut n, Time::new(0.03), controls(), || false)
        .unwrap();
    assert_eq!(fields(&n), fields(&reference));
    assert_eq!(wall.last_report(), reference_wall.last_report());

    let before = fields(&n);
    let contacts = vec![
        WallContact {
            drag: drag(false, 0.8),
            thermal: WallThermal::Adiabatic,
        },
        WallContact {
            drag: drag(false, 1e6),
            thermal: WallThermal::Adiabatic,
        },
    ];
    let mut refused =
        NetworkWallOperator::new(&n, vec![contacts.clone()], vec![], || false).unwrap();
    assert!(
        refused
            .advance_for(&mut n, Time::new(0.2), controls(), || false)
            .is_err()
    );
    assert_eq!(fields(&n), before);
    assert_eq!(refused.branch_contacts(0).unwrap(), contacts);
    assert_eq!(refused.last_report(), None);
}

#[test]
fn moving_wall_plenums_and_nonfinite_motion_remain_unsupported() {
    let g = gas(1, false, 0.0);
    let n = network(&g, &[0.0], true);
    let fixed = WallContact {
        drag: WallDrag::None,
        thermal: WallThermal::Adiabatic,
    };
    assert!(
        NetworkWallOperator::new(
            &n,
            vec![vec![fixed; 2]],
            vec![WallContact {
                drag: drag(false, 0.3),
                ..fixed
            }],
            || false
        )
        .is_err()
    );
    for speed in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(
            NetworkWallOperator::new(
                &n,
                vec![vec![
                    WallContact {
                        drag: drag(true, speed),
                        ..fixed
                    };
                    2
                ]],
                vec![fixed],
                || false
            )
            .is_err()
        );
    }
}
