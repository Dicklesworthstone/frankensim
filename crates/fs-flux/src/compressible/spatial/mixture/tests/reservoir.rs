use super::*;

fn states(
    gas: &MixtureGas,
    grid: &CartesianGrid,
    densities: &[f64],
    velocity: [f64; 3],
    pressure: f64,
) -> Vec<f64> {
    let width = gas.species_count() + 4;
    let mut values = vec![0.0; grid.shape().iter().product::<usize>() * width];
    for row in values.chunks_exact_mut(width) {
        state_into(gas, densities, velocity, pressure, row).unwrap();
    }
    values
}

fn supplied(
    gas: MixtureGas,
    grid: CartesianGrid,
    topology: Boundaries3d,
    initial: ([f64; 3], [f64; 3], f64),
    inputs: &[ReservoirBoundary<'_>],
) -> MixtureEuler3d {
    let values = states(&gas, &grid, &initial.0, initial.1, initial.2);
    let bytes =
        MixtureEuler3d::required_workspace_bytes_with_reservoirs(&gas, &grid, inputs.len()).unwrap();
    MixtureEuler3d::new_with_reservoirs(
        gas, grid, values, topology, inputs, bytes, || false,
    )
    .unwrap()
}

fn one_open_side(side: BoundarySide) -> Boundaries3d {
    let mut axes = [AxisBoundaries::Periodic; 3];
    axes[side.axis()] = AxisBoundaries::Ends {
        lower: if side.sign() < 0.0 {
            BoundaryCondition::Transmissive
        } else {
            BoundaryCondition::Reflecting
        },
        upper: if side.sign() > 0.0 {
            BoundaryCondition::Transmissive
        } else {
            BoundaryCondition::Reflecting
        },
    };
    Boundaries3d {
        x: axes[0],
        y: axes[1],
        z: axes[2],
    }
}

#[test]
fn g0_volumetric_mixture_reservoirs_preserve_uniform_composition_and_pressure() {
    let grid = CartesianGrid::new(3, 2, 2, 0.2, 0.3, 0.4).unwrap();
    let density = [0.7, 0.3, 0.0];
    for reconstruction in [Reconstruction::FirstOrder, Reconstruction::Muscl] {
        for velocity in [[0.0; 3], [0.2, -0.15, 0.1]] {
            let inputs = BoundarySide::ALL.map(|side| ReservoirBoundary {
                side,
                partial_densities: &density,
                velocity,
                pressure_pa: 2.0,
            });
            let source = standard(0.2);
            let mut flow = supplied(
                source.clone(),
                grid,
                Boundaries3d { x: OPEN, y: OPEN, z: OPEN },
                (density, velocity, 2.0),
                &inputs,
            );
            let before = flow.states().to_vec();
            for _ in 0..3 {
                balance(
                    flow.step(
                        Time::new(0.01),
                        StepControl { reconstruction, ..StepControl::default() },
                        || false,
                    )
                    .unwrap(),
                );
                for (actual, expected) in flow.states().iter().zip(&before) {
                    close(*actual, *expected, 3e-12);
                }
                assert_eq!(flow.gas(), &source);
                for input in inputs {
                    assert_eq!(flow.reservoir(input.side), Some(input));
                }
            }
        }
    }
}

#[test]
fn g1_pressure_reservoirs_transport_species_and_physical_energy_on_all_six_sides() {
    let grid = CartesianGrid::new(2, 2, 2, 0.1, 0.1, 0.1).unwrap();
    let donor = [0.0, 1.0, 0.0];
    for reconstruction in [Reconstruction::FirstOrder, Reconstruction::Muscl] {
        let mut reference: Option<(Totals, usize, f64)> = None;
        for side in BoundarySide::ALL {
            let input = ReservoirBoundary {
                side,
                partial_densities: &donor,
                velocity: [0.0; 3],
                pressure_pa: 2.0,
            };
            let mut flow = supplied(
                standard(0.2),
                grid,
                one_open_side(side),
                ([1.0, 0.0, 0.0], [0.0; 3], 1.0),
                &[input],
            );
            for _ in 0..4 {
                let report = flow
                    .step(
                        Time::new(0.002),
                        StepControl { reconstruction, ..StepControl::default() },
                        || false,
                    )
                    .unwrap();
                balance(report);
                assert!(report.boundary_change.mass() > 0.0);
                assert!(report.boundary_change.species_masses[1] > 0.0);
                let mut outward_species = [0.0; 3];
                let mut outward_energy = 0.0;
                for face in flow.boundary_impulses(side) {
                    for (total, mass) in outward_species.iter_mut().zip(face.species_masses) {
                        *total += mass;
                    }
                    outward_energy += face.energy
                        + face.species_masses.iter().zip(flow.gas().energy_references())
                            .map(|(mass, reference)| mass * reference).sum::<f64>();
                }
                for (s, outward) in outward_species.iter().enumerate() {
                    close(
                        report.after.species_masses[s] - report.before.species_masses[s],
                        -*outward,
                        3e-12,
                    );
                }
                close(
                    report.after.absolute_energy(flow.gas()).unwrap()
                        - report.before.absolute_energy(flow.gas()).unwrap(),
                    -outward_energy,
                    3e-10,
                );
            }
            let totals = flow.totals();
            assert!(totals.species_masses[1] > 0.0);
            assert_eq!(totals.species_masses[2], 0.0);
            assert!(-side.sign() * totals.momentum[side.axis()] > 0.0);
            for index in 0..flow.cell_count() {
                let q = flow.cell(index).unwrap();
                let p = q.primitive(flow.gas(), || false).unwrap();
                assert_eq!(q.partial_densities[2], 0.0);
                // The original species molecular masses give R_A=1, R_B=1/2.
                close(
                    p.temperature,
                    p.pressure / (q.partial_densities[0] + 0.5 * q.partial_densities[1]),
                    3e-12,
                );
            }
            if let Some((expected, axis, sign)) = reference {
                for s in 0..3 {
                    close(totals.species_masses[s], expected.species_masses[s], 3e-12);
                }
                close(totals.energy, expected.energy, 3e-12);
                close(
                    -side.sign() * totals.momentum[side.axis()],
                    -sign * expected.momentum[axis],
                    3e-12,
                );
            } else {
                reference = Some((totals, side.axis(), side.sign()));
            }
        }
    }
}

#[test]
fn g1_hot_mixture_reservoir_bounds_the_first_three_axis_step_without_retry() {
    let input = ReservoirBoundary {
        side: BoundarySide::XMin,
        partial_densities: &[0.5, 0.5, 0.0],
        velocity: [0.0; 3],
        pressure_pa: 32.0,
    };
    let mut flow = supplied(
        standard(0.0),
        CartesianGrid::new(2, 2, 2, 0.1, 0.1, 0.1).unwrap(),
        one_open_side(input.side),
        ([1.0, 0.0, 0.0], [0.0; 3], 0.01),
        &[input],
    );
    let controls = StepControl {
        cfl: 0.2,
        maximum_halvings: 0,
        reconstruction: Reconstruction::FirstOrder,
    };
    let report = flow.step(Time::new(1.0), controls, || false).unwrap();
    let bound = controls.cfl * 0.1 / (3.0 * (1.4_f64 * 32.0).sqrt());
    assert!(report.elapsed.value() <= bound * (1.0 + 2e-14));
    assert_eq!(report.halvings, 0);
    assert!(report.boundary_change.species_masses[1] > 0.0);
    balance(report);
}

#[test]
fn g1_pressure_fed_mixture_and_body_force_share_species_and_work_accounts() {
    let input = ReservoirBoundary {
        side: BoundarySide::XMin,
        partial_densities: &[0.0, 1.0, 0.0],
        velocity: [0.0; 3],
        pressure_pa: 2.0,
    };
    let mut flow = supplied(
        standard(0.2),
        CartesianGrid::new(2, 2, 2, 0.1, 0.1, 0.1).unwrap(),
        one_open_side(input.side),
        ([1.0, 0.0, 0.0], [0.0; 3], 1.0),
        &[input],
    );
    let mut workspace = flow.forcing_workspace(flow.required_forcing_workspace_bytes().unwrap())
        .unwrap();
    let report = flow.step_with_forcing(
        Time::new(0.001),
        StepControl::default(),
        BodyForce3d::UniformAcceleration([0.0, 0.25, 0.0]),
        ForceControl3d { maximum_mach: 2.0 },
        &mut workspace,
        || false,
    ).unwrap();
    assert!(report.boundary_change.species_masses[1] > 0.0);
    assert!(report.after.momentum[0] > 0.0);
    assert!(report.source_change.impulse[1] > 0.0);
    assert!(report.source_change.work > 0.0);
    balance(report);
}

#[test]
fn g4_mixture_reservoir_admission_budget_and_late_cancellation_are_atomic() {
    let source = standard(0.0);
    let grid = CartesianGrid::new(2, 2, 2, 0.1, 0.1, 0.1).unwrap();
    let input = ReservoirBoundary {
        side: BoundarySide::XMin,
        partial_densities: &[0.0, 1.0, 0.0],
        velocity: [0.0; 3],
        pressure_pa: 2.0,
    };
    let topology = one_open_side(input.side);
    let values = states(&source, &grid, &[1.0, 0.0, 0.0], [0.0; 3], 1.0);
    let bytes =
        MixtureEuler3d::required_workspace_bytes_with_reservoirs(&source, &grid, 1).unwrap();
    let make = |boundaries, inputs: &[ReservoirBoundary<'_>], budget| {
        MixtureEuler3d::new_with_reservoirs(
            source.clone(), grid, values.clone(), boundaries, inputs, budget, || false,
        )
    };
    assert!(make(topology, &[input], bytes - 1).is_err());
    assert!(make(PERIODIC, &[input], MAX_MIXTURE_WORKSPACE_BYTES).is_err());
    assert!(make(
        Boundaries3d { x: WALLS, y: WALLS, z: WALLS },
        &[input],
        MAX_MIXTURE_WORKSPACE_BYTES,
    ).is_err());
    assert!(make(topology, &[input, input], MAX_MIXTURE_WORKSPACE_BYTES).is_err());
    for bad in [
        ReservoirBoundary { partial_densities: &[1.0, 0.0], ..input },
        ReservoirBoundary { partial_densities: &[-0.1, 1.1, 0.0], ..input },
        ReservoirBoundary { pressure_pa: -1.0, ..input },
        ReservoirBoundary { pressure_pa: 100.0, ..input },
        ReservoirBoundary { velocity: [f64::NAN, 0.0, 0.0], ..input },
    ] {
        assert!(make(topology, &[bad], bytes).is_err());
    }
    let mut flow = make(topology, &[input], bytes).unwrap();
    assert_eq!(flow.workspace_bytes(), bytes);
    flow.step(Time::new(0.002), StepControl::default(), || false).unwrap();
    let before = flow.states().to_vec();
    let time = flow.time();
    let previous = *flow.last_step().unwrap();
    let loads = flow.boundary_impulses.clone();
    let mut probe = make(topology, &[input], bytes).unwrap();
    probe.step(Time::new(0.002), StepControl::default(), || false).unwrap();
    let mut calls = 0;
    probe.step(Time::new(0.002), StepControl::default(), || {
        calls += 1;
        false
    }).unwrap();
    let cutoff = calls * 3 / 4;
    let mut polls = 0;
    assert_eq!(
        flow.step(Time::new(0.002), StepControl::default(), || {
            polls += 1;
            polls == cutoff
        }),
        Err(Error::Cancelled),
    );
    assert_eq!(flow.states(), before);
    assert_eq!(flow.time(), time);
    assert_eq!(flow.last_step(), Some(&previous));
    assert_eq!(flow.boundary_impulses, loads);
    assert_eq!(flow.reservoir(input.side), Some(input));
}
