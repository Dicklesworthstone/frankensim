//! Moving unequal-caloric contacts, honest energy defects and atomic refusal.

use super::*;
use crate::multicomponent::pressure_equilibrium::{ContactEquilibriumControl, ContactWorkspace};

fn calorics(linear_cp: f64) -> MixtureGas {
    let mut domain = CaloricDomain::new(Temperature::new(0.001), Temperature::new(50.0));
    domain.absolute_energy_tolerance = 1e-13;
    domain.relative_energy_tolerance = 1e-14;
    MixtureGas::new(
        bank(vec![
            source("A", R, 3.5, 0.0, 1e18),
            source("B", 2.0 * R, 4.5, linear_cp, -1e18),
            source("absent", 3.0 * R, 5.0, 0.1, 7e17),
        ]),
        domain,
        || false,
    )
    .unwrap()
}

fn contact(gas: MixtureGas, velocity: f64, count: usize, boundaries: Boundaries) -> MixtureEuler1d {
    let states = (0..count)
        .map(|i| {
            let densities = if i < count / 2 {
                [1.0, 0.0, 0.0]
            } else {
                [0.0, 1.0, 0.0]
            };
            gas.conserved(&densities, velocity, 1.0).unwrap()
        })
        .collect();
    MixtureEuler1d::new(
        gas,
        Length::new(1.0 / count as f64),
        0.7,
        states,
        boundaries,
        || false,
    )
    .unwrap()
}

fn workspace(count: usize) -> ContactWorkspace {
    ContactWorkspace::new(
        count,
        ContactWorkspace::required_workspace_bytes(count).unwrap(),
        || false,
    )
    .unwrap()
}

fn control(budget: f64) -> ContactEquilibriumControl {
    ContactEquilibriumControl {
        relative_tolerance: 1e-10,
        maximum_energy_defect_j: budget,
    }
}

fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a} != {b}");
}

#[test]
fn moving_unequal_caloric_contacts_preserve_pressure_with_reported_energy_defects() {
    for linear_cp in [0.0, 0.4] {
        for velocity in [-0.7, 0.7] {
            for reconstruction in [Reconstruction::FirstOrder, Reconstruction::Muscl] {
                let mut flow = contact(calorics(linear_cp), velocity, 32, Boundaries::Periodic);
                let before = flow.totals();
                let mut scratch = workspace(flow.cell_count());
                let bytes = scratch.workspace_bytes();
                let controls = StepControl {
                    reconstruction,
                    maximum_halvings: 0,
                    ..StepControl::default()
                };
                let mut signed = 0.0;
                let mut absolute = 0.0;
                let mut steps = 0;
                // The material interface crosses several cells in both directions.
                while flow.time().value() < 0.15 {
                    let report = flow
                        .step_with_contact_equilibrium(
                            Time::new(0.15 - flow.time().value()),
                            controls,
                            control(1.0),
                            &mut scratch,
                            || false,
                        )
                        .unwrap();
                    assert_eq!(report.transport.hll_fallback_faces, 0);
                    assert_eq!(report.transport.halvings, 0);
                    assert!(!report.transport.first_order_retry);
                    assert!(report.correction.corrected_cells > 0);
                    assert!(
                        report.correction.absolute_energy_j
                            >= report.correction.signed_energy_j.abs()
                    );
                    near(report.corrected_energy_residual_j, 0.0, 3e-14);
                    for residual in &report.transport.conservation_residual.species_masses[..3] {
                        near(*residual, 0.0, 3e-14);
                    }
                    near(report.transport.conservation_residual.momentum, 0.0, 3e-14);
                    signed += report.correction.signed_energy_j;
                    absolute += report.correction.absolute_energy_j;
                    for i in 0..flow.cell_count() {
                        let q = flow.cell(i).unwrap();
                        let primitive = flow.gas().primitive(q, || false).unwrap();
                        near(primitive.pressure, 1.0, 3e-10);
                        near(primitive.velocity, velocity, 3e-10);
                        assert_eq!(q.partial_densities[2], 0.0);
                        assert!(q.partial_densities.iter().all(|rho| *rho >= 0.0));
                    }
                    steps += 1;
                    assert!(steps < 200);
                }
                assert!(absolute > 1e-4);
                let after = flow.totals();
                near(after.energy - before.energy, signed, 2e-13);
                for k in 0..3 {
                    near(after.species_masses[k], before.species_masses[k], 2e-13);
                }
                near(after.momentum, before.momentum, 2e-13);
                assert_eq!(scratch.workspace_bytes(), bytes);
            }
        }
    }
}

#[test]
fn non_equilibrium_stencils_and_physical_boundary_cells_keep_the_original_update() {
    for reconstruction in [Reconstruction::FirstOrder, Reconstruction::Muscl] {
        for pressure_gradient in [false, true] {
            let gas = calorics(0.4);
            let states: Vec<_> = (0..8)
                .map(|i| {
                    gas.conserved(
                        &[0.25 + 0.02 * i as f64, 0.75 - 0.02 * i as f64, 0.0],
                        if pressure_gradient {
                            0.7
                        } else {
                            0.05 * i as f64
                        },
                        if pressure_gradient {
                            1.0 + 0.03 * i as f64
                        } else {
                            1.0
                        },
                    )
                    .unwrap()
                })
                .collect();
            let make = || {
                MixtureEuler1d::new(
                    gas.clone(),
                    Length::new(0.125),
                    0.7,
                    states.clone(),
                    Boundaries::Periodic,
                    || false,
                )
                .unwrap()
            };
            let mut ordinary = make();
            let mut checked = make();
            let controls = StepControl {
                reconstruction,
                ..StepControl::default()
            };
            let expected = ordinary.step(Time::new(0.001), controls, || false).unwrap();
            let actual = checked
                .step_with_contact_equilibrium(
                    Time::new(0.001),
                    controls,
                    control(0.0),
                    &mut workspace(8),
                    || false,
                )
                .unwrap();
            assert_eq!(actual.correction.corrected_cells, 0);
            assert_eq!(actual.correction.absolute_energy_j, 0.0);
            assert_eq!(actual.transport, expected);
            assert_eq!(checked.states, ordinary.states);
        }
        let ends = Boundaries::Ends {
            left: Boundary::Transmissive,
            right: Boundary::Transmissive,
        };
        let mut ordinary = contact(calorics(0.4), 0.7, 2, ends);
        let mut checked = contact(calorics(0.4), 0.7, 2, ends);
        let controls = StepControl {
            reconstruction,
            ..StepControl::default()
        };
        let expected = ordinary.step(Time::new(0.001), controls, || false).unwrap();
        let actual = checked
            .step_with_contact_equilibrium(
                Time::new(0.001),
                controls,
                control(0.0),
                &mut workspace(2),
                || false,
            )
            .unwrap();
        assert_eq!(actual.correction.corrected_cells, 0);
        assert_eq!(actual.transport, expected);
        assert_eq!(checked.states, ordinary.states);
    }
}

#[test]
fn contact_budget_retries_preserve_muscl_and_refusal_is_atomic() {
    let make = || contact(calorics(0.4), 0.7, 16, Boundaries::Periodic);
    let controls = StepControl {
        reconstruction: Reconstruction::Muscl,
        maximum_halvings: 0,
        ..StepControl::default()
    };
    let mut reference = make();
    let report = reference
        .step_with_contact_equilibrium(
            Time::new(0.001),
            controls,
            control(1.0),
            &mut workspace(16),
            || false,
        )
        .unwrap();
    assert!(report.correction.absolute_energy_j > 0.0);
    let mut checked = make();
    let before = checked.states.clone();
    let mut scratch = workspace(16);
    let refused = checked.step_with_contact_equilibrium(
        Time::new(0.001),
        controls,
        control(0.0),
        &mut scratch,
        || false,
    );
    assert!(
        matches!(refused, Err(Error::StepRetries { last, .. }) if matches!(*last, Error::ContactEnergyBudget { .. }))
    );
    assert_eq!(checked.states, before);
    assert_eq!(checked.time().value(), 0.0);
    let budget = report.correction.absolute_energy_j / 4.0;
    let retried = checked
        .step_with_contact_equilibrium(
            Time::new(0.001),
            StepControl {
                maximum_halvings: 12,
                ..controls
            },
            control(budget),
            &mut scratch,
            || false,
        )
        .unwrap();
    assert!(retried.transport.halvings > 0);
    assert!(!retried.transport.first_order_retry);
    assert!(retried.correction.absolute_energy_j <= budget);
    near(retried.corrected_energy_residual_j, 0.0, 2e-14);
}

#[test]
fn late_contact_cancellation_retains_state_and_can_replay() {
    let make = || contact(calorics(0.4), -0.7, 16, Boundaries::Periodic);
    let mut reference = make();
    let mut calls = 0;
    let expected = reference
        .step_with_contact_equilibrium(
            Time::new(0.001),
            StepControl::default(),
            control(1.0),
            &mut workspace(16),
            || {
                calls += 1;
                false
            },
        )
        .unwrap();
    assert!(calls > 100);
    let mut checked = make();
    let before = checked.states.clone();
    let mut scratch = workspace(16);
    let mut cancelled_calls = 0;
    let refused = checked.step_with_contact_equilibrium(
        Time::new(0.001),
        StepControl::default(),
        control(1.0),
        &mut scratch,
        || {
            cancelled_calls += 1;
            cancelled_calls == calls
        },
    );
    assert_eq!(refused, Err(Error::Cancelled));
    assert_eq!(checked.states, before);
    assert_eq!(checked.time().value(), 0.0);
    let actual = checked
        .step_with_contact_equilibrium(
            Time::new(0.001),
            StepControl::default(),
            control(1.0),
            &mut scratch,
            || false,
        )
        .unwrap();
    assert_eq!(actual.transport.elapsed, expected.transport.elapsed);
    assert_eq!(
        actual.correction.corrected_cells,
        expected.correction.corrected_cells
    );
    // Cancellation can retain private caloric root guesses. The accepted
    // field above is bitwise untouched; a later solve agrees within the
    // original source-inversion tolerance, regardless of those warm starts.
    near(
        actual.correction.signed_energy_j,
        expected.correction.signed_energy_j,
        2e-14,
    );
    near(
        actual.correction.absolute_energy_j,
        expected.correction.absolute_energy_j,
        2e-14,
    );
    near(actual.corrected_energy_residual_j, 0.0, 2e-14);
    for (&actual, &expected) in checked.states.iter().zip(&reference.states) {
        near(actual, expected, 2e-13);
    }
}

#[test]
fn nonuniform_constant_area_contacts_work_and_variable_area_refuses_before_publication() {
    for varying_area in [false, true] {
        let gas = calorics(0.4);
        let positions = vec![0.0, 0.1, 0.3, 0.6, 1.0];
        let areas = if varying_area {
            vec![0.7, 0.9, 1.0, 0.9, 0.7]
        } else {
            vec![0.7; 5]
        };
        let volumes = positions
            .windows(2)
            .zip(areas.windows(2))
            .map(|(x, a)| (x[1] - x[0]) * 0.5 * (a[0] + a[1]))
            .collect();
        let mesh = Mesh1d::new(positions, areas, volumes).unwrap();
        let states = [
            [1.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
        ]
        .iter()
        .map(|d| gas.conserved(d, 0.7, 1.0).unwrap())
        .collect();
        let mut flow =
            MixtureEuler1d::new_on_mesh(gas, mesh, states, Boundaries::Periodic, || false).unwrap();
        let before = flow.states.clone();
        let result = flow.step_with_contact_equilibrium(
            Time::new(0.001),
            StepControl::default(),
            control(1.0),
            &mut workspace(4),
            || false,
        );
        if varying_area {
            assert!(matches!(result, Err(Error::Invalid(_))));
            assert_eq!(flow.states, before);
            assert_eq!(flow.time().value(), 0.0);
        } else {
            let report = result.unwrap();
            near(report.corrected_energy_residual_j, 0.0, 2e-14);
            for i in 0..4 {
                let p = flow
                    .gas()
                    .primitive(flow.cell(i).unwrap(), || false)
                    .unwrap();
                near(p.pressure, 1.0, 2e-12);
                near(p.velocity, 0.7, 2e-12);
            }
        }
    }
}
