use super::*;

fn topology(side: BoundarySide, pressure_pa: f64, opposite: BoundaryCondition) -> Boundaries2d {
    let outlet = BoundaryCondition::PressureOutlet { pressure_pa };
    let pair = AxisBoundaries::Ends {
        lower: if side.sign() < 0.0 { outlet } else { opposite },
        upper: if side.sign() > 0.0 { outlet } else { opposite },
    };
    if side.axis() == 0 {
        Boundaries2d {
            x: pair,
            ..PERIODIC
        }
    } else {
        Boundaries2d {
            y: pair,
            ..PERIODIC
        }
    }
}

#[test]
fn g0_pressure_outlet_uses_physical_hydrostatic_face_in_both_axes_and_gravity_directions() {
    use fs_qty::Acceleration;
    let grid = CartesianGrid::new(4, 3, 0.1, 0.2, 0.4).unwrap();
    for side in SIDES {
        for acceleration in [-8.0, 8.0] {
            for reconstruction in [Reconstruction::FirstOrder, Reconstruction::Muscl] {
                let reference = IsothermalAtmosphere {
                    acceleration: Acceleration::new(acceleration),
                    temperature: Temperature::new(2.0),
                    reference_pressure: Pressure::new(3.0),
                    reference_position: Length::new(0.0),
                };
                let closed = if side.axis() == 0 {
                    Boundaries2d {
                        x: WALLS,
                        ..PERIODIC
                    }
                } else {
                    Boundaries2d {
                        y: WALLS,
                        ..PERIODIC
                    }
                };
                let equilibrium = Euler2d::new_isothermal_atmosphere(
                    gas(0.4),
                    grid,
                    side.axis(),
                    reference,
                    closed,
                    || false,
                )
                .unwrap();
                let cell = match side {
                    BoundarySide::XMin | BoundarySide::YMin => 0,
                    BoundarySide::XMax => grid.shape[0] - 1,
                    BoundarySide::YMax => grid.shape[0] * (grid.shape[1] - 1),
                };
                let pressure = equilibrium.hydrostatic.as_ref().unwrap().trace(
                    &grid,
                    cell,
                    side.axis(),
                    side.sign(),
                )[3];
                let initial = equilibrium.states().to_vec();
                let mut flow = Euler2d::new_isothermal_atmosphere(
                    equilibrium.gas().clone(),
                    grid,
                    side.axis(),
                    reference,
                    topology(side, pressure, BoundaryCondition::Reflecting),
                    || false,
                )
                .unwrap();
                for _ in 0..3 {
                    let report = flow
                        .step(
                            Time::new(0.001),
                            StepControl {
                                reconstruction,
                                ..StepControl::default()
                            },
                            || false,
                        )
                        .unwrap();
                    check_balance(&report);
                    assert_eq!(report.halvings, 0);
                    close(report.boundary_change.mass, 0.0, 3e-12);
                    close(report.boundary_change.energy, 0.0, 3e-12);
                }
                for (actual, expected) in flow.states().iter().zip(initial) {
                    for (a, e) in actual.array().into_iter().zip(expected.array()) {
                        close(a, e, 3e-12);
                    }
                }
            }
        }
    }
}

#[test]
fn g0_planar_pressure_outlet_preserves_uniform_advection_and_shear_energy() {
    let grid = CartesianGrid::new(2, 3, 0.1, 0.2, 0.4).unwrap();
    for side in SIDES {
        for reconstruction in [Reconstruction::FirstOrder, Reconstruction::Muscl] {
            let mut velocity = [0.2, -0.3];
            velocity[side.axis()] = 0.4 * side.sign();
            let mut flow = field(
                gas(0.4),
                grid,
                topology(side, 1.0, BoundaryCondition::Transmissive),
                |_| [1.0, velocity[0], velocity[1], 1.0],
            );
            let initial = flow.states().to_vec();
            let report = flow
                .step(
                    Time::new(0.001),
                    StepControl {
                        reconstruction,
                        ..StepControl::default()
                    },
                    || false,
                )
                .unwrap();
            check_balance(&report);
            for (actual, expected) in flow.states().iter().zip(&initial) {
                for (a, e) in actual.array().into_iter().zip(expected.array()) {
                    close(a, e, 2e-11);
                }
            }
            for impulse in flow.boundary_impulses(side) {
                let mass = 0.4 * grid.face_area(side) * report.elapsed.value();
                close(impulse.mass, mass, 2e-12);
                close(
                    impulse.momentum[1 - side.axis()],
                    mass * velocity[1 - side.axis()],
                    2e-12,
                );
                close(impulse.energy, mass * (initial[0].energy + 1.0), 2e-11);
            }
        }
    }
}

#[test]
fn g1_planar_pressure_outlet_chokes_and_accounts_for_every_outgoing_face() {
    let grid = CartesianGrid::new(2, 3, 0.1, 0.2, 0.4).unwrap();
    for side in SIDES {
        for reconstruction in [Reconstruction::FirstOrder, Reconstruction::Muscl] {
            let mut velocity = [0.2, -0.3];
            velocity[side.axis()] = 0.0;
            let mut choked: Option<Vec<State>> = None;
            // All pressures, including vacuum, are beyond the sonic expansion.
            // The far-field tail leaves the source domain; the sampled face
            // remains valid and must be unchanged.
            for pressure in [0.01, 1e-12, 0.0] {
                let mut flow = field(
                    gas(0.0),
                    grid,
                    topology(side, pressure, BoundaryCondition::Reflecting),
                    |_| [1.0, velocity[0], velocity[1], 1.0],
                );
                let report = flow
                    .step(
                        Time::new(0.001),
                        StepControl {
                            reconstruction,
                            ..StepControl::default()
                        },
                        || false,
                    )
                    .unwrap();
                check_balance(&report);
                assert!(report.boundary_change.mass < 0.0);
                assert!(report.boundary_change.energy < 0.0);
                assert!(side.sign() * report.after.momentum[side.axis()] > 0.0);
                let mut mass = 0.0;
                let mut energy = 0.0;
                for impulse in flow.boundary_impulses(side) {
                    assert!(impulse.mass > 0.0);
                    close(
                        impulse.momentum[1 - side.axis()],
                        impulse.mass * velocity[1 - side.axis()],
                        2e-12,
                    );
                    mass += impulse.mass;
                    energy += impulse.energy;
                }
                close(mass, -report.boundary_change.mass, 2e-12);
                close(energy, -report.boundary_change.energy, 2e-12);
                if let Some(states) = &choked {
                    for (actual, expected) in flow.states().iter().zip(states) {
                        for (a, e) in actual.array().into_iter().zip(expected.array()) {
                            close(a, e, 2e-10);
                        }
                    }
                } else {
                    choked = Some(flow.states().to_vec());
                }
            }
        }
    }
}

#[test]
fn g4_planar_pressure_outlet_invalid_pressure_backflow_and_late_cancel_are_atomic() {
    let grid = CartesianGrid::new(2, 2, 0.1, 0.2, 0.4).unwrap();
    let source = gas(0.4);
    let state = State::from_primitive(&source, 1.0, [0.4, -0.2], 1.0).unwrap();
    let make = |pressure, side| {
        Euler2d::new(
            source.clone(),
            grid,
            vec![state; grid.len()],
            topology(side, pressure, BoundaryCondition::Transmissive),
            || false,
        )
    };
    for pressure in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(matches!(
            make(pressure, BoundarySide::XMax),
            Err(Error::Invalid(_))
        ));
    }
    let mut reversed = make(1.0, BoundarySide::XMin).unwrap();
    let before = reversed.states().to_vec();
    assert!(matches!(
        reversed.step(Time::new(0.001), StepControl::default(), || false),
        Err(Error::Invalid(_))
    ));
    assert_eq!(reversed.states(), before);
    assert_eq!(reversed.time().value(), 0.0);
    assert!(reversed.last_step().is_none());

    let mut probe = make(1.0, BoundarySide::XMax).unwrap();
    probe
        .step(Time::new(0.001), StepControl::default(), || false)
        .unwrap();
    let mut polls = 0;
    probe
        .step(Time::new(0.001), StepControl::default(), || {
            polls += 1;
            false
        })
        .unwrap();
    let mut flow = make(1.0, BoundarySide::XMax).unwrap();
    flow.step(Time::new(0.001), StepControl::default(), || false)
        .unwrap();
    let states = flow.states().to_vec();
    let loads = flow.boundary_impulses.clone();
    let report = *flow.last_step().unwrap();
    let time = flow.time();
    let mut count = 0;
    assert!(matches!(
        flow.step(Time::new(0.001), StepControl::default(), || {
            count += 1;
            count >= polls - 1
        }),
        Err(Error::Cancelled)
    ));
    assert_eq!(flow.states(), states);
    assert_eq!(flow.boundary_impulses, loads);
    assert_eq!(flow.last_step(), Some(&report));
    assert_eq!(flow.time(), time);
}


fn characteristic_topology(
    side: BoundarySide,
    reference_pressure_pa: f64,
    reference_velocity_m_s: f64,
) -> Boundaries2d {
    let outlet = BoundaryCondition::CharacteristicOutlet {
        reference_pressure_pa,
        reference_velocity_m_s,
    };
    let pair = AxisBoundaries::Ends {
        lower: if side.sign() < 0.0 { outlet } else { BoundaryCondition::Transmissive },
        upper: if side.sign() > 0.0 { outlet } else { BoundaryCondition::Transmissive },
    };
    if side.axis() == 0 {
        Boundaries2d { x: pair, ..PERIODIC }
    } else {
        Boundaries2d { y: pair, ..PERIODIC }
    }
}

#[test]
fn g0_characteristic_outlet_preserves_planar_advection_and_shear_on_every_side() {
    let grid = CartesianGrid::new(2, 3, 0.1, 0.2, 0.4).unwrap();
    for side in SIDES {
        for reconstruction in [Reconstruction::FirstOrder, Reconstruction::Muscl] {
            let mut velocity = [0.2, -0.3];
            velocity[side.axis()] = 0.4 * side.sign();
            let mut flow = field(
                gas(0.4),
                grid,
                characteristic_topology(side, 1.0, velocity[side.axis()]),
                |_| [1.0, velocity[0], velocity[1], 1.0],
            );
            let initial = flow.states().to_vec();
            let report = flow.step(
                Time::new(0.001),
                StepControl { reconstruction, ..StepControl::default() },
                || false,
            ).unwrap();
            check_balance(&report);
            assert_eq!(report.halvings, 0);
            for (actual, expected) in flow.states().iter().zip(initial) {
                for (a, e) in actual.array().into_iter().zip(expected.array()) {
                    close(a, e, 3e-11);
                }
            }
        }
    }
}

fn planar_outgoing_pulse_remainder(outlet: BoundaryCondition) -> f64 {
    let cells = 96;
    let grid = CartesianGrid::new(cells, 2, 1.0 / cells as f64, 0.5, 0.4).unwrap();
    let boundaries = Boundaries2d {
        x: AxisBoundaries::Ends {
            lower: BoundaryCondition::Transmissive,
            upper: outlet,
        },
        y: AxisBoundaries::Periodic,
    };
    let mut flow = field(gas(0.0), grid, boundaries, |x| {
        // Independent gamma=1.4 simple wave over the manufactured Cp/R=3.5 source.
        let pressure = 1.0 + 1e-3 * (-((x[0] - 0.65) / 0.055).powi(2)).exp();
        let density = pressure.powf(1.0 / 1.4);
        let sound_speed = (1.4 * pressure / density).sqrt();
        [density, 0.3 + 5.0 * (sound_speed - 1.4_f64.sqrt()), 0.2, pressure]
    });
    run(&mut flow, 0.6);
    flow.states().iter().map(|state| {
        (state.primitive(flow.gas(), || false).unwrap().pressure - 1.0).powi(2)
    }).sum::<f64>().sqrt()
}

#[test]
fn g1_characteristic_outlet_releases_a_pressure_pulse_from_the_planar_rk_owner() {
    let outgoing = planar_outgoing_pulse_remainder(BoundaryCondition::CharacteristicOutlet {
        reference_pressure_pa: 1.0,
        reference_velocity_m_s: 0.3,
    });
    let fixed = planar_outgoing_pulse_remainder(BoundaryCondition::PressureOutlet {
        pressure_pa: 1.0,
    });
    assert!(fixed > 1e-4, "fixed-pressure reflection={fixed:e}");
    assert!(outgoing < 0.05 * fixed, "characteristic={outgoing:e}, fixed={fixed:e}");
}

#[test]
fn g4_characteristic_outlet_admission_backflow_and_cancellation_are_atomic() {
    let grid = CartesianGrid::new(2, 2, 0.1, 0.1, 0.4).unwrap();
    for (pressure, velocity) in [(0.0, 0.3), (f64::NAN, 0.3), (1.0, -0.3), (1.0, f64::INFINITY)] {
        assert!(characteristic_topology(BoundarySide::XMax, pressure, velocity)
            .x.admit::<Error>().is_err());
    }
    let mut flow = field(
        gas(0.4), grid,
        characteristic_topology(BoundarySide::XMax, 1.0, 0.3),
        |_| [1.0, -0.2, 0.1, 1.0],
    );
    let initial = flow.states().to_vec();
    assert!(flow.step(Time::new(0.001), StepControl::default(), || false).is_err());
    assert_eq!(flow.states(), initial);
    assert_eq!(flow.time().value(), 0.0);
    assert!(matches!(
        flow.step(Time::new(0.001), StepControl::default(), || true),
        Err(Error::Cancelled)
    ));
    assert_eq!(flow.states(), initial);
    assert_eq!(flow.time().value(), 0.0);
}
