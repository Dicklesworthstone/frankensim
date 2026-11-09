use super::*;

fn topology(side: BoundarySide, pressure_pa: f64, opposite: BoundaryCondition) -> Boundaries3d {
    let outlet = BoundaryCondition::PressureOutlet { pressure_pa };
    let pair = AxisBoundaries::Ends {
        lower: if side.sign() < 0.0 { outlet } else { opposite },
        upper: if side.sign() > 0.0 { outlet } else { opposite },
    };
    let mut axes = [AxisBoundaries::Periodic; 3];
    axes[side.axis()] = pair;
    Boundaries3d {
        x: axes[0],
        y: axes[1],
        z: axes[2],
    }
}

#[test]
fn g1_volumetric_pressure_outlets_preserve_species_both_shears_and_energy_on_every_side() {
    let grid = CartesianGrid::new(2, 3, 4, 0.1, 0.2, 0.3).unwrap();
    for side in BoundarySide::ALL {
        for reconstruction in [Reconstruction::FirstOrder, Reconstruction::Muscl] {
            let mut choked: Option<Vec<f64>> = None;
            for pressure in [1.0, 0.8, 0.01, 1e-12, 0.0] {
                let mut velocity = [0.2, -0.3, 0.1];
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
                    for tangent in tangents(side.axis()) {
                        close(
                            impulse.momentum[tangent],
                            velocity[tangent] * outward,
                            3e-12,
                        );
                    }
                    mass += outward;
                    energy += impulse.energy;
                }
                for index in 0..flow.cell_count() {
                    let q = flow.cell(index).unwrap();
                    let p = q.primitive(flow.gas(), || false).unwrap();
                    assert_eq!(q.partial_densities[2], 0.0);
                    close(q.partial_densities[0] / p.density, 0.35, 3e-11);
                    for tangent in tangents(side.axis()) {
                        close(p.velocity[tangent], velocity[tangent], 3e-11);
                    }
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
fn g1_volumetric_mixture_pressure_outlet_bounds_the_initial_interval_without_retry() {
    let grid = CartesianGrid::new(2, 2, 2, 0.1, 0.1, 0.1).unwrap();
    let mut flow = field(
        standard(0.0),
        grid,
        topology(BoundarySide::ZMin, 1e-12, BoundaryCondition::Reflecting),
        |_| ([0.35, 0.65, 0.0], [0.0; 3], 1.0),
    );
    let controls = StepControl {
        cfl: 0.1,
        maximum_halvings: 0,
        reconstruction: Reconstruction::FirstOrder,
    };
    let report = flow.step(Time::new(1.0), controls, || false).unwrap();
    let sonic = 2.0 * 1.4_f64.sqrt() / 2.4;
    let bound = controls.cfl * 0.1 / (3.0 * 2.0 * sonic);
    assert!(report.elapsed.value() <= bound * (1.0 + 2e-10));
    assert_eq!(report.halvings, 0);
    assert!(report.boundary_change.mass() < 0.0);
    balance(report);
}

#[test]
fn g4_volumetric_mixture_pressure_outlet_refuses_invalid_pressure_and_reverse_flow() {
    let grid = CartesianGrid::new(2, 2, 2, 0.1, 0.2, 0.3).unwrap();
    let initial = field(standard(0.4), grid, PERIODIC, |_| {
        ([0.35, 0.65, 0.0], [0.2, -0.1, 0.4], 1.0)
    });
    let make = |pressure| {
        MixtureEuler3d::new(
            initial.gas().clone(),
            grid,
            initial.states().to_vec(),
            topology(
                BoundarySide::ZMin,
                pressure,
                BoundaryCondition::Transmissive,
            ),
            initial.workspace_bytes(),
            || false,
        )
    };
    for pressure in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(matches!(make(pressure), Err(Error::Invalid(_))));
    }
    let input = ReservoirBoundary {
        side: BoundarySide::ZMin,
        partial_densities: &[0.35, 0.65, 0.0],
        velocity: [0.2, -0.1, 0.4],
        pressure_pa: 1.0,
    };
    let bytes =
        MixtureEuler3d::required_workspace_bytes_with_reservoirs(initial.gas(), &grid, 1).unwrap();
    assert!(matches!(
        MixtureEuler3d::new_with_reservoirs(
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
    let mut reversed = make(1.0).unwrap();
    assert!(matches!(
        reversed.step(Time::new(0.001), StepControl::default(), || false),
        Err(Error::Invalid(_))
    ));
    assert_eq!(reversed.states(), initial.states());
    assert_eq!(reversed.time().value(), 0.0);
    assert!(reversed.last_step().is_none());
    assert!(reversed.boundary_impulses.iter().all(|&value| value == 0.0));
}
