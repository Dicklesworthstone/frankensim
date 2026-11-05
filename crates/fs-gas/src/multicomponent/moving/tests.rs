use super::*;
use crate::multicomponent::tests::{gas, source};
use fs_qty::{Length, Temperature};
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;

const CLOSED: Boundaries = Boundaries::Ends {
    left: Boundary::Reflecting,
    right: Boundary::Reflecting,
};

fn mixture(references: [f64; 3]) -> MixtureGas {
    gas(vec![
        source("A", R, 3.5, 0.4, references[0]),
        source("B", 2.0 * R, 4.5, 0.2, references[1]),
        source("absent", 0.5 * R, 3.0, 0.1, references[2]),
    ])
}

fn flow(count: usize, references: [f64; 3], velocity: f64) -> MixtureEuler1d {
    let gas = mixture(references);
    let states = (0..count)
        .map(|i| {
            let first = if i < count / 2 { 0.8 } else { 0.2 };
            gas.conserved(&[first, 1.0 - first, 0.0], velocity, 1.0)
                .unwrap()
        })
        .collect();
    MixtureEuler1d::new(
        gas,
        Length::new(1.0 / count as f64),
        0.7,
        states,
        CLOSED,
        || false,
    )
    .unwrap()
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}"
    );
}

#[test]
fn g0_species_ale_flux_and_stationary_step_retain_original_conventions() {
    let gas = mixture([0.0; 3]);
    let mut q = [0.0; 5];
    let p = gas.state_at(&[0.3, 0.8, 0.0], 0.4, 1.2, &mut q).unwrap();
    let mut fixed = [0.0; 5];
    gas.interface(&q, p, &q, p, &mut fixed, &mut || false)
        .unwrap();
    for speed in [-0.7, 0.0, 1.3] {
        let mut moving = [0.0; 5];
        let receipt = gas
            .moving_interface(&q, p, &q, p, speed, &mut moving, &mut || false)
            .unwrap();
        for k in 0..5 {
            close(moving[k], fixed[k] - speed * q[k], 3e-12);
        }
        assert_eq!(moving[2], 0.0);
        close(
            receipt.signal_speed,
            (p.velocity - speed).abs() + p.sound_speed,
            3e-14,
        );
        if speed == 0.0 {
            assert_eq!(moving, fixed);
        }
    }
    let mut ordinary = flow(16, [0.0; 3], 0.1);
    let mut moving = flow(16, [0.0; 3], 0.1);
    for _ in 0..4 {
        let expected = ordinary
            .step(Time::new(0.002), StepControl::default(), || false)
            .unwrap();
        let actual = moving
            .step_with_moving_walls(Time::new(0.002), StepControl::default(), [0.0; 2], || false)
            .unwrap();
        assert_eq!(actual, expected);
        assert_eq!(moving.states, ordinary.states);
        assert_eq!(moving.mesh(), ordinary.mesh());
        assert_eq!(moving.time(), ordinary.time());
        assert!(moving.moving.is_none());
    }
}

#[test]
fn g0_comoving_unequal_caloric_contacts_preserve_species_and_pressure() {
    let speed = 0.125;
    let mut field = flow(24, [0.0; 3], speed);
    let original = field.states.clone();
    let positions = field.mesh().positions().to_vec();
    let volumes = field.mesh().volumes().to_vec();
    let before = field.totals();
    for _ in 0..24 {
        let report = field
            .step_with_moving_walls(Time::new(0.002), StepControl::default(), [speed; 2], || {
                false
            })
            .unwrap();
        for k in 0..3 {
            assert_eq!(report.boundary_change.species_masses[k], 0.0);
            close(
                report.after.species_masses[k],
                before.species_masses[k],
                3e-13,
            );
        }
        assert_eq!(field.mesh().volumes(), volumes);
    }
    for i in 0..field.cell_count() {
        let q = field.cell(i).unwrap();
        let p = field.gas().primitive(q, || false).unwrap();
        close(p.pressure, 1.0, 3e-10);
        close(p.velocity, speed, 3e-10);
        assert_eq!(q.partial_densities[2], 0.0);
    }
    for (&actual, &expected) in field.states.iter().zip(&original) {
        close(actual, expected, 3e-10);
    }
    for (&x, &start) in field.mesh().positions().iter().zip(&positions) {
        close(x, start + speed * field.time().value(), 3e-14);
    }
}

#[test]
fn g1_mixture_piston_preserves_each_species_and_independent_energy_gauges() {
    let mut field = flow(32, [0.0; 3], 0.0);
    let mut shifted = flow(32, [1e18, -7e17, 3e17], 0.0);
    let before = field.totals();
    let mut work = 0.0;
    let mut impulse = 0.0;
    let mut steps = 0;
    while field.time().value() < 0.08 {
        let dt = Time::new(0.08 - field.time().value());
        let report = field
            .step_with_moving_walls(dt, StepControl::default(), [0.2, 0.0], || false)
            .unwrap();
        let other = shifted
            .step_with_moving_walls(dt, StepControl::default(), [0.2, 0.0], || false)
            .unwrap();
        assert_eq!(report, other);
        assert_eq!(field.states, shifted.states);
        assert_eq!(field.mesh(), shifted.mesh());
        for k in 0..3 {
            assert_eq!(report.boundary_change.species_masses[k], 0.0);
            close(report.conservation_residual.species_masses[k], 0.0, 3e-13);
        }
        for q in field.states.chunks_exact(5) {
            assert_eq!(q[2], 0.0, "moving pressure must not seed absent species");
        }
        assert!(report.left_boundary_change.energy > 0.0);
        assert_eq!(report.right_boundary_change.energy, 0.0);
        close(report.conservation_residual.energy, 0.0, 3e-12);
        assert!(report.maximum_cfl <= 0.40000000000001);
        work += report.boundary_change.energy;
        impulse += report.boundary_change.momentum;
        steps += 1;
        assert!(steps < 1000);
    }
    let after = field.totals();
    for k in 0..3 {
        close(after.species_masses[k], before.species_masses[k], 3e-12);
    }
    close(after.energy - before.energy, work, 3e-12);
    close(after.momentum - before.momentum, impulse, 3e-12);
    close(
        field.mesh().volumes().iter().sum(),
        0.7 * (1.0 - 0.2 * 0.08),
        3e-14,
    );
    assert!(
        field
            .gas()
            .primitive(field.cell(0).unwrap(), || false)
            .unwrap()
            .pressure
            > 1.1
    );
}

#[test]
fn g1_single_species_piston_reduces_to_the_original_frozen_owner() {
    use crate::dynamics::{CaloricDomain, Euler1d, FrozenGas};
    use fs_thermochem::mixture::FrozenIdealGasMixtureModelV1;
    use fs_thermochem::{Composition, CompositionBasis};
    let model = source("A", R, 3.5, 0.4, 0.0);
    let mixed_gas = gas(vec![model.clone()]);
    let frozen_gas = FrozenGas::new(
        FrozenIdealGasMixtureModelV1::new(
            vec![model],
            Composition::new(CompositionBasis::MoleFraction, vec![1.0]).unwrap(),
        )
        .unwrap(),
        CaloricDomain::new(Temperature::new(0.001), Temperature::new(50.0)),
        || false,
    )
    .unwrap();
    let mixed_q = mixed_gas.conserved(&[1.0], 0.0, 1.0).unwrap();
    let frozen_q = frozen_gas.conserved(1.0, 0.0, 1.0).unwrap();
    let mut mixed = MixtureEuler1d::new(
        mixed_gas,
        Length::new(0.05),
        0.7,
        vec![mixed_q; 20],
        CLOSED,
        || false,
    )
    .unwrap();
    let mut frozen = Euler1d::new(
        frozen_gas,
        Length::new(0.05),
        0.7,
        vec![frozen_q; 20],
        CLOSED,
        || false,
    )
    .unwrap();
    for _ in 0..12 {
        let a = mixed
            .step_with_moving_walls(
                Time::new(0.0005),
                StepControl::default(),
                [0.2, 0.0],
                || false,
            )
            .unwrap();
        let b = frozen
            .step_with_moving_walls(
                Time::new(0.0005),
                StepControl::default(),
                [0.2, 0.0],
                || false,
            )
            .unwrap();
        assert_eq!(a.elapsed, b.elapsed);
        assert_eq!(mixed.mesh(), frozen.mesh());
        close(a.boundary_change.energy, b.boundary_change.energy, 3e-11);
        for (i, &q) in frozen.states().iter().enumerate() {
            let m = mixed.cell(i).unwrap();
            close(m.partial_densities[0], q.density, 3e-10);
            close(m.momentum, q.momentum, 3e-10);
            close(m.energy, q.energy, 3e-10);
        }
    }
}

#[test]
fn g4_mixture_moving_geometry_and_late_cancellation_preserve_complete_state() {
    let mut counted = flow(12, [0.0; 3], 0.0);
    let mut polls = 0;
    counted
        .step_with_moving_walls(Time::new(0.003), StepControl::default(), [0.2, 0.0], || {
            polls += 1;
            false
        })
        .unwrap();
    // Five polls interrupts the shared geometry velocity preparation; its
    // dynamics cancellation must not become a retryable Geometry(Cancelled).
    for cutoff in [1, 5, polls / 2, polls] {
        let mut field = flow(12, [0.0; 3], 0.0);
        let states = field.states.clone();
        let mesh = field.mesh().clone();
        let mut seen = 0;
        assert_eq!(
            field.step_with_moving_walls(
                Time::new(0.003),
                StepControl::default(),
                [0.2, 0.0],
                || {
                    seen += 1;
                    seen >= cutoff
                }
            ),
            Err(Error::Cancelled)
        );
        assert_eq!(seen, cutoff);
        assert_eq!(field.states, states);
        assert_eq!(field.mesh(), &mesh);
        assert_eq!(field.time().value(), 0.0);
    }
    for velocity in [[f64::NAN, 0.0], [100.0, 0.0]] {
        let mut field = flow(12, [0.0; 3], 0.0);
        let states = field.states.clone();
        let mesh = field.mesh().clone();
        assert!(
            field
                .step_with_moving_walls(
                    Time::new(0.003),
                    StepControl {
                        maximum_halvings: 0,
                        ..StepControl::default()
                    },
                    velocity,
                    || false,
                )
                .is_err()
        );
        assert_eq!(field.states, states);
        assert_eq!(field.mesh(), &mesh);
        assert_eq!(field.time().value(), 0.0);
    }
    let gas = mixture([0.0; 3]);
    let q = gas.conserved(&[0.3, 0.7, 0.0], 0.0, 1.0).unwrap();
    let mut periodic = MixtureEuler1d::new(
        gas.clone(),
        Length::new(0.5),
        0.7,
        vec![q.clone(); 2],
        Boundaries::Periodic,
        || false,
    )
    .unwrap();
    assert!(
        periodic
            .step_with_moving_walls(Time::new(0.003), StepControl::default(), [0.2, 0.0], || {
                false
            })
            .is_err()
    );
    let mesh = Mesh1d::new(vec![0.0, 0.5, 1.0], vec![1.0; 3], vec![0.4; 2]).unwrap();
    let mut invalid =
        MixtureEuler1d::new_on_mesh(gas, mesh.clone(), vec![q; 2], CLOSED, || false).unwrap();
    assert!(
        invalid
            .step_with_moving_walls(Time::new(0.003), StepControl::default(), [0.2, 0.0], || {
                false
            })
            .is_err()
    );
    assert_eq!(invalid.mesh(), &mesh);
    assert_eq!(invalid.time().value(), 0.0);
}

fn contact_field(count: usize, velocity: f64) -> MixtureEuler1d {
    let source = mixture([1e18, -7e17, 3e17]).source().clone();
    let mut domain =
        crate::dynamics::CaloricDomain::new(Temperature::new(0.001), Temperature::new(50.0));
    domain.absolute_energy_tolerance = 1e-13;
    domain.relative_energy_tolerance = 1e-14;
    let gas = MixtureGas::new(source, domain, || false).unwrap();
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
        CLOSED,
        || false,
    )
    .unwrap()
}

fn contact_workspace(count: usize) -> ContactWorkspace {
    ContactWorkspace::new(
        count,
        ContactWorkspace::required_workspace_bytes(count).unwrap(),
        || false,
    )
    .unwrap()
}

fn contact_control(budget: f64) -> ContactEquilibriumControl {
    ContactEquilibriumControl {
        relative_tolerance: 1e-10,
        maximum_energy_defect_j: budget,
    }
}

#[test]
fn g0_ale_contact_pressure_survives_translation_and_volume_change_with_accounted_work() {
    for wall_velocity in [[0.4, 0.4], [-0.4, -0.4], [0.25, -0.125], [-0.25, 0.125]] {
        for reconstruction in [Reconstruction::FirstOrder, Reconstruction::Muscl] {
            let mut ordinary = contact_field(64, 0.0);
            let mut checked = contact_field(64, 0.0);
            let before = checked.totals();
            let mut workspace = contact_workspace(64);
            let bytes = workspace.workspace_bytes();
            let controls = StepControl {
                reconstruction,
                maximum_halvings: 0,
                ..StepControl::default()
            };
            let (mut wall_work, mut wall_impulse, mut defect, mut absolute) = (0.0, 0.0, 0.0, 0.0);
            for _ in 0..24 {
                ordinary
                    .step_with_moving_walls(Time::new(0.001), controls, wall_velocity, || false)
                    .unwrap();
                let report = checked
                    .step_with_moving_walls_and_contact_equilibrium(
                        Time::new(0.001),
                        controls,
                        wall_velocity,
                        contact_control(1.0),
                        &mut workspace,
                        || false,
                    )
                    .unwrap();
                assert_eq!(report.transport.halvings, 0);
                assert!(!report.transport.first_order_retry);
                assert!(
                    report.correction.absolute_energy_j >= report.correction.signed_energy_j.abs()
                );
                close(report.corrected_energy_residual_j, 0.0, 2e-12);
                close(report.transport.conservation_residual.momentum, 0.0, 2e-12);
                for k in 0..3 {
                    assert_eq!(report.transport.boundary_change.species_masses[k], 0.0);
                    close(
                        report.transport.conservation_residual.species_masses[k],
                        0.0,
                        2e-12,
                    );
                }
                for q in checked.states.chunks_exact(5) {
                    assert_eq!(q[2], 0.0);
                    assert!(q[..3].iter().all(|rho| *rho >= 0.0));
                }
                wall_work += report.transport.boundary_change.energy;
                wall_impulse += report.transport.boundary_change.momentum;
                defect += report.correction.signed_energy_j;
                absolute += report.correction.absolute_energy_j;
            }
            assert_eq!(ordinary.time(), checked.time());
            assert_eq!(ordinary.mesh(), checked.mesh());
            // The material front moves relative to the mesh. The short time
            // leaves this interior window outside the physical wall waves.
            let mut ordinary_error = 0.0_f64;
            let mut corrected_error = 0.0_f64;
            for i in 24..40 {
                let old = ordinary
                    .gas()
                    .primitive(ordinary.cell(i).unwrap(), || false)
                    .unwrap();
                let new = checked
                    .gas()
                    .primitive(checked.cell(i).unwrap(), || false)
                    .unwrap();
                ordinary_error = ordinary_error.max((old.pressure - 1.0).abs());
                corrected_error = corrected_error.max((new.pressure - 1.0).abs());
                close(new.velocity, 0.0, 3e-10);
            }
            assert!(
                ordinary_error > 1e-4,
                "baseline pressure error: {ordinary_error:e}"
            );
            assert!(
                corrected_error < 3e-10,
                "corrected pressure error: {corrected_error:e}"
            );
            assert!(absolute > 1e-7);
            assert!(wall_work.abs() > 1e-6);
            let after = checked.totals();
            close(after.energy - before.energy - wall_work, defect, 3e-12);
            close(after.momentum - before.momentum, wall_impulse, 3e-12);
            for k in 0..3 {
                close(after.species_masses[k], before.species_masses[k], 3e-12);
            }
            assert_eq!(workspace.workspace_bytes(), bytes);
        }
    }
}

#[test]
fn g0_ale_contact_zero_motion_and_resolved_acoustic_gradients_keep_native_maps() {
    let mut fixed = contact_field(16, 0.7);
    let mut moving = contact_field(16, 0.7);
    let mut fixed_workspace = contact_workspace(16);
    let mut moving_workspace = contact_workspace(16);
    for _ in 0..3 {
        let expected = fixed
            .step_with_contact_equilibrium(
                Time::new(0.001),
                StepControl::default(),
                contact_control(1.0),
                &mut fixed_workspace,
                || false,
            )
            .unwrap();
        let actual = moving
            .step_with_moving_walls_and_contact_equilibrium(
                Time::new(0.001),
                StepControl::default(),
                [0.0; 2],
                contact_control(1.0),
                &mut moving_workspace,
                || false,
            )
            .unwrap();
        assert_eq!(actual, expected);
        assert_eq!(moving.states, fixed.states);
        assert_eq!(moving.mesh(), fixed.mesh());
        assert_eq!(moving.time(), fixed.time());
        assert!(moving.moving.is_none());
    }
    for pressure_gradient in [false, true] {
        let make = || {
            let gas = mixture([0.0; 3]);
            let states = (0..12)
                .map(|i| {
                    gas.conserved(
                        &[0.2 + 0.03 * i as f64, 0.8 - 0.03 * i as f64, 0.0],
                        if pressure_gradient {
                            0.1
                        } else {
                            0.03 * i as f64
                        },
                        if pressure_gradient {
                            1.0 + 0.02 * i as f64
                        } else {
                            1.0
                        },
                    )
                    .unwrap()
                })
                .collect();
            MixtureEuler1d::new(gas, Length::new(1.0 / 12.0), 0.7, states, CLOSED, || false)
                .unwrap()
        };
        let mut ordinary = make();
        let mut checked = make();
        let expected = ordinary
            .step_with_moving_walls(
                Time::new(0.001),
                StepControl::default(),
                [0.2, -0.1],
                || false,
            )
            .unwrap();
        let actual = checked
            .step_with_moving_walls_and_contact_equilibrium(
                Time::new(0.001),
                StepControl::default(),
                [0.2, -0.1],
                contact_control(0.0),
                &mut contact_workspace(12),
                || false,
            )
            .unwrap();
        assert_eq!(actual.correction.corrected_cells, 0);
        assert_eq!(actual.correction.absolute_energy_j, 0.0);
        assert_eq!(actual.transport, expected);
        assert_eq!(checked.states, ordinary.states);
        assert_eq!(checked.mesh(), ordinary.mesh());
    }
}

#[test]
fn g4_ale_contact_budget_retry_retains_order_and_refusal_preserves_geometry() {
    let controls = StepControl {
        maximum_halvings: 0,
        reconstruction: Reconstruction::Muscl,
        ..StepControl::default()
    };
    let mut reference = contact_field(16, 0.0);
    let full = reference
        .step_with_moving_walls_and_contact_equilibrium(
            Time::new(0.002),
            controls,
            [0.25, -0.125],
            contact_control(1.0),
            &mut contact_workspace(16),
            || false,
        )
        .unwrap();
    assert!(full.correction.absolute_energy_j > 0.0);
    let mut checked = contact_field(16, 0.0);
    let before = checked.states.clone();
    let mesh = checked.mesh().clone();
    let mut workspace = contact_workspace(16);
    let refused = checked.step_with_moving_walls_and_contact_equilibrium(
        Time::new(0.002),
        controls,
        [0.25, -0.125],
        contact_control(0.0),
        &mut workspace,
        || false,
    );
    assert!(
        matches!(refused, Err(Error::StepRetries { last, .. }) if matches!(*last, Error::ContactEnergyBudget { .. }))
    );
    assert_eq!(checked.states, before);
    assert_eq!(checked.mesh(), &mesh);
    assert_eq!(checked.time().value(), 0.0);
    let budget = full.correction.absolute_energy_j / 4.0;
    let report = checked
        .step_with_moving_walls_and_contact_equilibrium(
            Time::new(0.002),
            StepControl {
                maximum_halvings: 12,
                ..controls
            },
            [0.25, -0.125],
            contact_control(budget),
            &mut workspace,
            || false,
        )
        .unwrap();
    assert!(report.transport.halvings > 0);
    assert!(!report.transport.first_order_retry);
    assert!(report.correction.absolute_energy_j <= budget);
    close(report.corrected_energy_residual_j, 0.0, 2e-13);
    close(
        checked.mesh().positions()[0],
        0.25 * report.transport.elapsed.value(),
        2e-15,
    );
    close(
        checked.mesh().positions()[16],
        1.0 - 0.125 * report.transport.elapsed.value(),
        2e-15,
    );
}

#[test]
fn g4_ale_contact_late_cancellation_is_atomic_and_reusable() {
    let mut reference = contact_field(16, 0.0);
    let mut calls = 0;
    let expected = reference
        .step_with_moving_walls_and_contact_equilibrium(
            Time::new(0.002),
            StepControl::default(),
            [0.25, -0.125],
            contact_control(1.0),
            &mut contact_workspace(16),
            || {
                calls += 1;
                false
            },
        )
        .unwrap();
    assert!(calls > 100);
    for cutoff in [1, calls / 2, calls] {
        let mut checked = contact_field(16, 0.0);
        let before = checked.states.clone();
        let mesh = checked.mesh().clone();
        let mut workspace = contact_workspace(16);
        let mut seen = 0;
        let refused = checked.step_with_moving_walls_and_contact_equilibrium(
            Time::new(0.002),
            StepControl::default(),
            [0.25, -0.125],
            contact_control(1.0),
            &mut workspace,
            || {
                seen += 1;
                seen == cutoff
            },
        );
        assert_eq!(refused, Err(Error::Cancelled));
        assert_eq!(checked.states, before);
        assert_eq!(checked.mesh(), &mesh);
        assert_eq!(checked.time().value(), 0.0);
        let actual = checked
            .step_with_moving_walls_and_contact_equilibrium(
                Time::new(0.002),
                StepControl::default(),
                [0.25, -0.125],
                contact_control(1.0),
                &mut workspace,
                || false,
            )
            .unwrap();
        assert_eq!(actual.transport.elapsed, expected.transport.elapsed);
        assert_eq!(checked.mesh(), reference.mesh());
        close(
            actual.correction.signed_energy_j,
            expected.correction.signed_energy_j,
            3e-13,
        );
        close(actual.corrected_energy_residual_j, 0.0, 3e-13);
        for (&a, &b) in checked.states.iter().zip(&reference.states) {
            close(a, b, 3e-12);
        }
    }
}
