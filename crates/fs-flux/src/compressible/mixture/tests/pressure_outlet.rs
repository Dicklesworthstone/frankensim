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
fn g1_planar_mixture_pressure_outlets_preserve_species_shear_and_energy_on_every_side() {
    let grid = CartesianGrid::new(2, 3, 0.1, 0.2, 0.4).unwrap();
    for side in SIDES {
        for reconstruction in [Reconstruction::FirstOrder, Reconstruction::Muscl] {
            let mut choked: Option<Vec<f64>> = None;
            for pressure in [1.0, 0.8, 0.01, 1e-12, 0.0] {
                let mut velocity = [0.2, -0.3];
                velocity[side.axis()] = if pressure == 1.0 {
                    0.4 * side.sign()
                } else {
                    0.0
                };
                let opposite = if pressure == 1.0 {
                    BoundaryCondition::Transmissive
                } else {
                    BoundaryCondition::Reflecting
                };
                let mut flow = field(
                    standard(0.4),
                    grid,
                    topology(side, pressure, opposite),
                    |_| ([0.35, 0.65, 0.0], velocity, 1.0),
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
                    .unwrap_or_else(|error| panic!("side={side:?}, reconstruction={reconstruction:?}, pressure={pressure:e}: {error:?}"));
                balance(report);
                let mut mass = 0.0;
                let mut energy = 0.0;
                for impulse in flow.boundary_impulses(side) {
                    let outward = impulse.species_masses.iter().sum::<f64>();
                    assert!(outward > 0.0);
                    close(impulse.species_masses[0], 0.35 * outward, 3e-12);
                    close(impulse.species_masses[1], 0.65 * outward, 3e-12);
                    assert_eq!(impulse.species_masses[2], 0.0);
                    close(
                        impulse.momentum[1 - side.axis()],
                        velocity[1 - side.axis()] * outward,
                        3e-12,
                    );
                    mass += outward;
                    energy += impulse.energy;
                }
                for index in 0..flow.cell_count() {
                    let q = flow.cell(index).unwrap();
                    let p = q.primitive(flow.gas(), || false).unwrap();
                    assert_eq!(q.partial_densities[2], 0.0);
                    close(q.partial_densities[0] / p.density, 0.35, 3e-11);
                    close(
                        p.velocity[1 - side.axis()],
                        velocity[1 - side.axis()],
                        3e-11,
                    );
                }
                if pressure == 1.0 {
                    for (&actual, expected) in flow.states().iter().zip(initial) {
                        close(actual, expected, 3e-11);
                    }
                } else {
                    assert!(report.boundary_change.mass() < 0.0);
                    close(mass, -report.boundary_change.mass(), 3e-12);
                    close(energy, -report.boundary_change.energy, 3e-12);
                }
                if pressure < 0.1 {
                    if let Some(states) = &choked {
                        for (&actual, &expected) in flow.states().iter().zip(states) {
                            close(actual, expected, 3e-10);
                        }
                    } else {
                        choked = Some(flow.states().to_vec());
                    }
                }
            }
        }
    }
}

#[test]
fn g4_planar_mixture_pressure_outlet_refuses_invalid_data_backflow_and_reservoir_override() {
    let grid = CartesianGrid::new(2, 2, 0.1, 0.2, 0.4).unwrap();
    let initial = field(standard(0.4), grid, PERIODIC, |_| {
        ([0.35, 0.65, 0.0], [0.4, -0.2], 1.0)
    });
    let make = |pressure, side| {
        MixtureEuler2d::new(
            initial.gas().clone(),
            grid,
            initial.states().to_vec(),
            topology(side, pressure, BoundaryCondition::Transmissive),
            initial.workspace_bytes(),
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
    assert!(matches!(
        reversed.step(Time::new(0.001), StepControl::default(), || false),
        Err(Error::Invalid(_))
    ));
    assert_eq!(reversed.states(), initial.states());
    assert_eq!(reversed.time().value(), 0.0);
    assert!(reversed.last_step().is_none());
    let input = ReservoirBoundary {
        side: BoundarySide::XMax,
        partial_densities: &[0.35, 0.65, 0.0],
        velocity: [0.4, -0.2],
        pressure_pa: 1.0,
    };
    let bytes =
        MixtureEuler2d::required_workspace_bytes_with_reservoirs(initial.gas(), &grid, 1).unwrap();
    assert!(matches!(
        MixtureEuler2d::new_with_reservoirs(
            initial.gas().clone(),
            grid,
            initial.states().to_vec(),
            topology(input.side, 1.0, BoundaryCondition::Transmissive),
            &[input],
            bytes,
            || false,
        ),
        Err(Error::Invalid(_))
    ));

    let mut flow = make(1.0, BoundarySide::XMax).unwrap();
    flow.step(Time::new(0.001), StepControl::default(), || false)
        .unwrap();
    let states = flow.states().to_vec();
    let loads = flow.boundary_impulses.clone();
    let previous = *flow.last_step().unwrap();
    let time = flow.time();
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
    assert_eq!(flow.last_step(), Some(&previous));
    assert_eq!(flow.time(), time);
}


#[test]
fn g0_characteristic_outlet_preserves_planar_species_entropy_and_shear() {
    let grid = CartesianGrid::new(2, 3, 0.1, 0.2, 0.4).unwrap();
    for side in SIDES {
        for reconstruction in [Reconstruction::FirstOrder, Reconstruction::Muscl] {
            let mut velocity = [0.2, -0.3];
            velocity[side.axis()] = 0.4 * side.sign();
            let outlet = BoundaryCondition::CharacteristicOutlet {
                reference_pressure_pa: 1.0,
                reference_velocity_m_s: velocity[side.axis()],
            };
            let pair = AxisBoundaries::Ends {
                lower: if side.sign() < 0.0 { outlet } else { BoundaryCondition::Transmissive },
                upper: if side.sign() > 0.0 { outlet } else { BoundaryCondition::Transmissive },
            };
            let boundaries = if side.axis() == 0 {
                Boundaries2d { x: pair, ..PERIODIC }
            } else {
                Boundaries2d { y: pair, ..PERIODIC }
            };
            let mut flow = field(standard(0.4), grid, boundaries, |_| {
                ([0.35, 0.65, 0.0], velocity, 1.0)
            });
            let initial = flow.states().to_vec();
            let report = flow.step(
                Time::new(0.001),
                StepControl { reconstruction, ..StepControl::default() },
                || false,
            ).unwrap();
            balance(report);
            for (&actual, expected) in flow.states().iter().zip(initial) {
                close(actual, expected, 3e-11);
            }
            for impulse in flow.boundary_impulses(side) {
                let mass: f64 = impulse.species_masses.iter().sum();
                assert!(mass > 0.0);
                close(impulse.species_masses[0], 0.35 * mass, 3e-12);
                close(impulse.species_masses[1], 0.65 * mass, 3e-12);
                assert_eq!(impulse.species_masses[2], 0.0);
                close(
                    impulse.momentum[1 - side.axis()],
                    velocity[1 - side.axis()] * mass,
                    3e-12,
                );
            }
        }
    }
}
