use super::*;

fn unequal_gas(linear: f64) -> MixtureGas {
    gas(vec![
        source("A", R, 3.5, 0.0, 30.0),
        source("B", 2.0 * R, 4.5, linear, -60.0),
        source("absent-C", 0.5 * R, 4.5, 0.2, 7.0),
    ])
}

fn control() -> ContactEquilibriumControl {
    ContactEquilibriumControl {
        relative_tolerance: 1e-8,
        maximum_energy_defect_j: 10.0,
    }
}

fn workspace(cells: usize) -> ContactWorkspace {
    ContactWorkspace::new(
        cells,
        ContactWorkspace::required_workspace_bytes(cells).unwrap(),
        || false,
    )
    .unwrap()
}

fn contact(linear: f64, direction: usize) -> MixtureEuler2d {
    let velocity = match direction {
        0 => [0.35, 0.0],
        1 => [0.0, 0.35],
        _ => [0.35, -0.2],
    };
    field(
        unequal_gas(linear),
        CartesianGrid::new(8, 8, 0.125, 0.125, 0.3).unwrap(),
        PERIODIC,
        |x| {
            let coordinate = if direction < 2 {
                x[direction]
            } else {
                0.5 * (x[0] + x[1])
            };
            let a = if coordinate < 0.5 { 1.0 } else { 0.0 };
            ([a, 1.0 - a, 0.0], velocity, 1.0)
        },
    )
}

#[test]
fn g1_planar_contact_equilibrium_preserves_unequal_calorics_and_accounts_energy() {
    for linear in [0.0, 0.4] {
        for direction in 0..3 {
            for reconstruction in [Reconstruction::FirstOrder, Reconstruction::Muscl] {
                let mut flow = contact(linear, direction);
                let mut conservative = contact(linear, direction);
                let mut work = workspace(flow.cell_count());
                let controls = StepControl {
                    reconstruction,
                    ..StepControl::default()
                };
                conservative
                    .step(Time::new(0.01), controls, || false)
                    .unwrap();
                assert!((0..conservative.cell_count()).any(|k| {
                    (conservative
                        .cell(k)
                        .unwrap()
                        .primitive(conservative.gas(), || false)
                        .unwrap()
                        .pressure
                        - 1.0)
                        .abs()
                        > 1e-5
                }));
                let velocity = flow
                    .cell(0)
                    .unwrap()
                    .primitive(flow.gas(), || false)
                    .unwrap()
                    .velocity;
                for step in 0..3 {
                    let report = flow
                        .step_with_contact_equilibrium(
                            Time::new(0.01),
                            controls,
                            control(),
                            &mut work,
                            || false,
                        )
                        .unwrap();
                    assert!(report.correction.corrected_cells > 0);
                    assert!(report.correction.corrected_cells <= flow.cell_count());
                    assert!(report.correction.absolute_energy_j > 0.0);
                    assert!(
                        report.correction.absolute_energy_j + 1e-14
                            >= report.correction.signed_energy_j.abs()
                    );
                    close(report.corrected_energy_residual_j, 0.0, 3e-10);
                    close(
                        report.transport.conservation_residual.energy,
                        report.correction.signed_energy_j,
                        3e-10,
                    );
                    assert_eq!(report.transport.source_change, ForceChange2d::default());
                    for value in report
                        .transport
                        .conservation_residual
                        .species_masses
                        .into_iter()
                        .chain(report.transport.conservation_residual.momentum)
                    {
                        close(value, 0.0, 3e-12);
                    }
                    for k in 0..flow.cell_count() {
                        let q = flow.cell(k).unwrap();
                        let p = q.primitive(flow.gas(), || false).unwrap();
                        let [a, b, absent] = <[f64; 3]>::try_from(q.partial_densities).unwrap();
                        assert_eq!(absent, 0.0);
                        close(p.density, 1.0, 3e-11);
                        close(p.pressure, 1.0, 3e-9);
                        for axis in 0..2 {
                            close(p.velocity[axis], velocity[axis], 3e-10);
                        }
                        let temperature = 1.0 / (a + 0.5 * b);
                        close(p.temperature, temperature, 3e-9);
                        // Independent NASA polynomial integral, including both
                        // reference offsets and the complete shear kinetic term.
                        let thermal = (2.5 * a + 1.75 * b) * temperature
                            + 0.25 * linear * b * temperature.powi(2)
                            + 30.0 * a
                            - 30.0 * b;
                        let physical = q.energy
                            + a * flow.gas().energy_references()[0]
                            + b * flow.gas().energy_references()[1];
                        close(
                            physical,
                            thermal + 0.5 * velocity.iter().map(|u| u * u).sum::<f64>(),
                            3e-9,
                        );
                    }
                    if step == 0 && direction < 2 && reconstruction == Reconstruction::FirstOrder {
                        let k = if direction == 0 { 4 } else { 4 * 8 };
                        let courant: f64 = 0.35 * 0.01 / 0.125;
                        // Exact SSP-RK2 upwind mass entering the first pure-B cell.
                        close(
                            flow.cell(k).unwrap().partial_densities[0],
                            courant - 0.5 * courant.powi(2),
                            3e-12,
                        );
                    }
                }
            }
        }
    }
}

#[test]
fn g1_planar_contact_mode_leaves_pressure_gradients_shocks_and_shear_bit_identical() {
    for reconstruction in [Reconstruction::FirstOrder, Reconstruction::Muscl] {
        for case in 0..3 {
            let make = || {
                field(
                    unequal_gas(0.4),
                    CartesianGrid::new(8, 4, 0.125, 0.25, 0.3).unwrap(),
                    PERIODIC,
                    |x| {
                        let wave = (2.0 * std::f64::consts::PI * x[0]).sin();
                        let a = 0.5 + 0.2 * wave;
                        let pressure = match case {
                            0 => 1.0 + 0.1 * wave,
                            1 => {
                                if (x[0] * 8.0) as usize % 2 == 0 {
                                    1.0
                                } else {
                                    0.5
                                }
                            }
                            _ => 1.0,
                        };
                        (
                            [a, 1.0 - a, 0.0],
                            [0.3, if case == 2 { 0.1 * wave } else { 0.2 }],
                            pressure,
                        )
                    },
                )
            };
            let mut ordinary = make();
            let mut optional = make();
            let mut work = workspace(optional.cell_count());
            let controls = StepControl {
                reconstruction,
                ..StepControl::default()
            };
            let original = ordinary.step(Time::new(0.002), controls, || false).unwrap();
            let projected = optional
                .step_with_contact_equilibrium(
                    Time::new(0.002),
                    controls,
                    control(),
                    &mut work,
                    || false,
                )
                .unwrap();
            assert_eq!(
                optional.states(),
                ordinary.states(),
                "case={case:?}, order={reconstruction:?}"
            );
            assert_eq!(projected.transport, original);
            assert_eq!(projected.correction, Default::default());
        }
    }
}

#[test]
fn g1_planar_contact_mode_retains_original_physical_pressure_boundary_path() {
    for condition in [
        BoundaryCondition::Reflecting,
        BoundaryCondition::PressureOutlet { pressure_pa: 0.8 },
    ] {
        let make = || {
            field(
                unequal_gas(0.4),
                CartesianGrid::new(2, 4, 0.25, 0.25, 0.3).unwrap(),
                Boundaries2d {
                    x: AxisBoundaries::Ends {
                        lower: condition,
                        upper: condition,
                    },
                    ..PERIODIC
                },
                |x| {
                    let a = if x[1] < 0.5 { 1.0 } else { 0.0 };
                    ([a, 1.0 - a, 0.0], [0.0, 0.3], 1.0)
                },
            )
        };
        let mut ordinary = make();
        let mut optional = make();
        let mut work = workspace(optional.cell_count());
        let original = ordinary
            .step(Time::new(0.001), StepControl::default(), || false)
            .unwrap();
        let projected = optional
            .step_with_contact_equilibrium(
                Time::new(0.001),
                StepControl::default(),
                control(),
                &mut work,
                || false,
            )
            .unwrap();
        assert_eq!(optional.states(), ordinary.states());
        assert_eq!(projected.transport, original);
        assert_eq!(projected.correction, Default::default());
        assert_eq!(optional.boundary_impulses, ordinary.boundary_impulses);
    }
}

#[test]
fn g1_planar_contact_budget_and_cancellation_never_publish_a_trial() {
    let mut flow = contact(0.4, 2);
    let initial = flow.states().to_vec();
    let mut work = workspace(flow.cell_count());
    let refused = flow.step_with_contact_equilibrium(
        Time::new(0.01),
        StepControl {
            maximum_halvings: 0,
            ..StepControl::default()
        },
        ContactEquilibriumControl {
            maximum_energy_defect_j: 0.0,
            ..control()
        },
        &mut work,
        || false,
    );
    assert!(
        matches!(refused, Err(Error::StepRetries { last, .. }) if matches!(*last, Error::ContactEnergyBudget { .. }))
    );
    assert_eq!(flow.states(), initial);
    assert_eq!(flow.time().value(), 0.0);
    assert!(flow.last_step().is_none());
    assert!(flow.boundary_impulses.iter().all(|v| *v == 0.0));

    let mut polls = 0;
    flow.step_with_contact_equilibrium(
        Time::new(0.01),
        StepControl::default(),
        control(),
        &mut work,
        || {
            polls += 1;
            false
        },
    )
    .unwrap();
    for stop in [1, polls / 2, polls] {
        let mut cancelled = contact(0.4, 2);
        let mut trial_work = workspace(cancelled.cell_count());
        let mut count = 0;
        let result = cancelled.step_with_contact_equilibrium(
            Time::new(0.01),
            StepControl::default(),
            control(),
            &mut trial_work,
            || {
                count += 1;
                count == stop
            },
        );
        assert!(
            matches!(result, Err(Error::Cancelled)),
            "stop={stop}, count={count}"
        );
        assert_eq!(cancelled.states(), initial);
        assert_eq!(cancelled.time().value(), 0.0);
        assert!(cancelled.last_step().is_none());
        assert!(cancelled.boundary_impulses.iter().all(|v| *v == 0.0));
    }
}
