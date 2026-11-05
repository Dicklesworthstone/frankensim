//! Actual caloric contacts driving the native inertial wall and suspension.

use super::*;

fn contact_field(count: usize, references: [f64; 3]) -> MixtureEuler1d {
    let mut domain =
        crate::dynamics::CaloricDomain::new(Temperature::new(0.001), Temperature::new(50.0));
    domain.absolute_energy_tolerance = 1e-13;
    domain.relative_energy_tolerance = 1e-14;
    let gas = MixtureGas::new(mixture(references).source().clone(), domain, || false).unwrap();
    let states = (0..count)
        .map(|i| {
            let densities = if i < count / 2 {
                [1.0, 0.0, 0.0]
            } else {
                [0.0, 1.0, 0.0]
            };
            gas.conserved(&densities, 0.0, 1.0).unwrap()
        })
        .collect();
    MixtureEuler1d::new(
        gas,
        Length::new(1.0 / count as f64),
        AREA,
        states,
        CLOSED,
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

fn equilibrium(budget: f64) -> ContactEquilibriumControl {
    ContactEquilibriumControl {
        relative_tolerance: 1e-10,
        maximum_energy_defect_j: budget,
    }
}

#[test]
fn g0_contact_piston_retains_pressure_species_impulse_and_accounted_energy() {
    for direction in [-1.0, 1.0] {
        for reconstruction in [Reconstruction::FirstOrder, Reconstruction::Muscl] {
            let mut flow = contact_field(64, [1e18, -7e17, 3e17]);
            let mut native = contact_field(64, [1e18, -7e17, 3e17]);
            let mut piston = Piston::new(0.4, direction * 0.25).unwrap();
            let original = piston;
            let mut native_piston = piston;
            let before = flow.totals();
            let mut scratch = workspace(64);
            let bytes = scratch.workspace_bytes();
            let gas_controls = StepControl {
                reconstruction,
                maximum_halvings: 0,
                ..StepControl::default()
            };
            let load = PistonLoad {
                exterior_pressure_pa: 1.0 - direction * 0.3,
                external_force_n: direction * 0.05,
            };
            let (mut work, mut impulse, mut support, mut defect, mut absolute) =
                (0.0, 0.0, 0.0, 0.0, 0.0);
            for _ in 0..24 {
                native
                    .step_with_piston(
                        &mut native_piston,
                        Time::new(0.001),
                        gas_controls,
                        load,
                        controls(),
                        || false,
                    )
                    .unwrap();
                let step = flow
                    .step_with_piston_and_contact_equilibrium(
                        &mut piston,
                        Time::new(0.001),
                        gas_controls,
                        load,
                        controls(),
                        equilibrium(1.0),
                        &mut scratch,
                        || false,
                    )
                    .unwrap();
                let gas = step.gas.transport;
                let correction = step.gas.correction;
                assert!(step.admitted(controls()));
                assert_eq!(gas.halvings, 0);
                assert!(!gas.first_order_retry);
                assert!(correction.absolute_energy_j >= correction.signed_energy_j.abs());
                close(step.gas.corrected_energy_residual_j, 0.0, 3e-12);
                close(gas.conservation_residual.momentum, 0.0, 3e-12);
                close(
                    gas.right_boundary_change.energy,
                    -step.midpoint_velocity_m_s * step.gas_impulse_n_s,
                    3e-14,
                );
                let raw_total = gas.after.energy - gas.before.energy
                    + step.kinetic_energy_change_j
                    + step.exterior_pressure_work_j
                    - step.external_force_work_j;
                close(
                    raw_total - correction.signed_energy_j,
                    step.total_energy_residual_j,
                    3e-14,
                );
                for k in 0..3 {
                    assert_eq!(gas.boundary_change.species_masses[k], 0.0);
                    close(gas.conservation_residual.species_masses[k], 0.0, 3e-12);
                }
                for cell in flow.states.chunks_exact(5) {
                    assert_eq!(cell[2], 0.0);
                    assert!(cell[..3].iter().all(|rho| *rho >= 0.0));
                }
                work += step.external_force_work_j - step.exterior_pressure_work_j;
                impulse += gas.elapsed.value()
                    * (load.external_force_n - AREA * load.exterior_pressure_pa);
                support += step.fixed_wall_impulse_n_s;
                defect += correction.signed_energy_j;
                absolute += correction.absolute_energy_j;
            }
            // The short window leaves this material front outside wall waves.
            let (mut native_error, mut corrected_error) = (0.0_f64, 0.0_f64);
            for i in 24..40 {
                let old = native
                    .gas()
                    .primitive(native.cell(i).unwrap(), || false)
                    .unwrap();
                let new = flow
                    .gas()
                    .primitive(flow.cell(i).unwrap(), || false)
                    .unwrap();
                native_error = native_error.max((old.pressure - 1.0).abs());
                corrected_error = corrected_error.max((new.pressure - 1.0).abs());
                close(new.velocity, 0.0, 3e-10);
            }
            assert!(
                native_error > 1e-4,
                "native pressure error {native_error:e}"
            );
            assert!(
                corrected_error < 3e-10,
                "corrected pressure error {corrected_error:e}"
            );
            assert!(absolute > 1e-7);
            let after = flow.totals();
            let kinetic = 0.5
                * piston.mass_kg()
                * (piston.velocity_m_s().powi(2) - original.velocity_m_s().powi(2));
            close(after.energy - before.energy + kinetic - work, defect, 3e-10);
            close(
                after.momentum - before.momentum
                    + piston.mass_kg() * (piston.velocity_m_s() - original.velocity_m_s())
                    + support,
                impulse,
                3e-10,
            );
            for k in 0..3 {
                close(after.species_masses[k], before.species_masses[k], 3e-12);
            }
            assert_eq!(scratch.workspace_bytes(), bytes);
        }
    }
}

#[test]
fn g0_final_corrected_trial_replays_and_zero_suspension_is_exact() {
    let mut coupled = contact_field(24, [0.0; 3]);
    let mut replay = contact_field(24, [0.0; 3]);
    let mut suspended = contact_field(24, [0.0; 3]);
    let mut piston = Piston::new(0.4, -0.125).unwrap();
    let mut other = piston;
    let load = PistonLoad {
        exterior_pressure_pa: 1.7,
        external_force_n: 0.05,
    };
    let step = coupled
        .step_with_piston_and_contact_equilibrium(
            &mut piston,
            Time::new(0.003),
            StepControl::default(),
            load,
            controls(),
            equilibrium(1.0),
            &mut workspace(24),
            || false,
        )
        .unwrap();
    assert!(
        step.gas_trials > 1,
        "must exercise multiple scalar proposals"
    );
    let expected = replay
        .step_with_moving_walls_and_contact_equilibrium(
            step.gas.transport.elapsed,
            StepControl::default(),
            [0.0, step.midpoint_velocity_m_s],
            equilibrium(1.0),
            &mut workspace(24),
            || false,
        )
        .unwrap();
    assert_eq!(coupled.mesh(), replay.mesh());
    assert_eq!(coupled.time(), replay.time());
    for (&actual, &expected) in coupled.states.iter().zip(&replay.states) {
        close(actual, expected, 3e-12);
    }
    close(
        step.gas.correction.signed_energy_j,
        expected.correction.signed_energy_j,
        3e-12,
    );
    close(
        step.gas.correction.absolute_energy_j,
        expected.correction.absolute_energy_j,
        3e-12,
    );
    assert_eq!(
        step.gas.correction.corrected_cells,
        expected.correction.corrected_cells
    );
    close(
        step.gas.transport.right_boundary_change.energy,
        expected.transport.right_boundary_change.energy,
        3e-12,
    );
    let zero = suspended
        .step_with_suspended_piston_and_contact_equilibrium(
            &mut other,
            Time::new(0.003),
            StepControl::default(),
            load,
            PistonSuspension::new(0.0, 0.0, -3.0).unwrap(),
            controls(),
            equilibrium(1.0),
            &mut workspace(24),
            || false,
        )
        .unwrap();
    assert_eq!(zero.piston, step);
    assert_eq!(other, piston);
    assert_eq!(suspended.states, coupled.states);
    assert_eq!(suspended.mesh(), coupled.mesh());
}

#[test]
fn g1_suspended_contact_preserves_gauges_and_separates_loss_from_numerical_defect() {
    let mut flow = contact_field(24, [0.0; 3]);
    let mut shifted = contact_field(24, [1e18, -7e17, 3e17]);
    let mut piston = Piston::new(0.7, 0.125).unwrap();
    let mut other = piston;
    let spring = PistonSuspension::new(8.0, 1.2, 0.92).unwrap();
    let load = PistonLoad {
        exterior_pressure_pa: 1.3,
        external_force_n: 0.05,
    };
    let before = flow.totals();
    let initial =
        before.energy + 0.5 * 0.7 * 0.125_f64.powi(2) + spring.stored_energy_j(1.0).unwrap();
    let (mut heat, mut work, mut defect, mut absolute) = (0.0, 0.0, 0.0, 0.0);
    let mut scratch = workspace(24);
    let mut other_scratch = workspace(24);
    for _ in 0..32 {
        let step = flow
            .step_with_suspended_piston_and_contact_equilibrium(
                &mut piston,
                Time::new(0.001),
                StepControl::default(),
                load,
                spring,
                controls(),
                equilibrium(1.0),
                &mut scratch,
                || false,
            )
            .unwrap();
        let reference = shifted
            .step_with_suspended_piston_and_contact_equilibrium(
                &mut other,
                Time::new(0.001),
                StepControl::default(),
                load,
                spring,
                controls(),
                equilibrium(1.0),
                &mut other_scratch,
                || false,
            )
            .unwrap();
        assert_eq!(step, reference);
        assert_eq!(piston, other);
        assert_eq!(flow.states, shifted.states);
        assert_eq!(flow.mesh(), shifted.mesh());
        assert!(step.total_energy_residual_j.abs() <= controls().energy_tolerance_j);
        assert!(step.damping_energy_j >= 0.0);
        heat += step.damping_energy_j;
        work += step.applied_force_work_j - step.piston.exterior_pressure_work_j;
        defect += step.piston.gas.correction.signed_energy_j;
        absolute += step.piston.gas.correction.absolute_energy_j;
    }
    let after = flow.totals();
    let final_energy = after.energy
        + 0.5 * 0.7 * piston.velocity_m_s().powi(2)
        + spring.stored_energy_j(flow.mesh().positions()[24]).unwrap();
    close(final_energy - initial + heat - work, defect, 4e-10);
    assert!(heat > 1e-7);
    assert!(absolute > 1e-7);
    for k in 0..3 {
        close(after.species_masses[k], before.species_masses[k], 3e-12);
    }
}

#[test]
fn g0_resolved_pressure_and_velocity_gradients_keep_native_piston_map() {
    for pressure_gradient in [false, true] {
        let make = || {
            let gas = mixture([0.0; 3]);
            let states = (0..12)
                .map(|i| {
                    gas.conserved(
                        &[0.25 + 0.03 * i as f64, 0.75 - 0.03 * i as f64, 0.0],
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
            MixtureEuler1d::new(gas, Length::new(1.0 / 12.0), AREA, states, CLOSED, || false)
                .unwrap()
        };
        let mut native = make();
        let mut corrected = make();
        let mut a = Piston::new(0.4, 0.125).unwrap();
        let mut b = a;
        let load = PistonLoad {
            exterior_pressure_pa: 1.5,
            external_force_n: 0.05,
        };
        let expected = native
            .step_with_piston(
                &mut a,
                Time::new(0.001),
                StepControl::default(),
                load,
                controls(),
                || false,
            )
            .unwrap();
        let actual = corrected
            .step_with_piston_and_contact_equilibrium(
                &mut b,
                Time::new(0.001),
                StepControl::default(),
                load,
                controls(),
                equilibrium(0.0),
                &mut workspace(12),
                || false,
            )
            .unwrap();
        assert_eq!(actual.gas.correction, ContactCorrection::default());
        assert_eq!(actual.map_gas(|report| report.transport), expected);
        assert_eq!(a, b);
        assert_eq!(corrected.states, native.states);
        assert_eq!(corrected.mesh(), native.mesh());
    }
}

#[test]
fn g4_piston_contact_budget_retry_and_cancellation_publish_no_partial_owners() {
    let no_retry = StepControl {
        maximum_halvings: 0,
        ..StepControl::default()
    };
    let load = PistonLoad {
        exterior_pressure_pa: 1.6,
        external_force_n: 0.05,
    };
    let duration = Time::new(0.003);
    let initial = Piston::new(0.4, -0.125).unwrap();
    let mut baseline = contact_field(16, [0.0; 3]);
    let mut baseline_piston = initial;
    let mut polls = 0;
    let full = baseline
        .step_with_piston_and_contact_equilibrium(
            &mut baseline_piston,
            duration,
            no_retry,
            load,
            controls(),
            equilibrium(1.0),
            &mut workspace(16),
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    assert!(full.gas.correction.absolute_energy_j > 0.0);
    let mut flow = contact_field(16, [0.0; 3]);
    let state = flow.states.clone();
    let mesh = flow.mesh().clone();
    let mut piston = initial;
    let mut scratch = workspace(16);
    let refused = flow.step_with_piston_and_contact_equilibrium(
        &mut piston,
        duration,
        no_retry,
        load,
        controls(),
        equilibrium(0.0),
        &mut scratch,
        || false,
    );
    assert!(matches!(refused, Err(PistonError::StepRetries { last, .. })
        if matches!(*last, PistonError::MixtureGas(Error::ContactEnergyBudget { .. }))));
    assert_eq!(flow.states, state);
    assert_eq!(flow.mesh(), &mesh);
    assert_eq!(flow.time().value(), 0.0);
    assert_eq!(piston, initial);
    let budget = full.gas.correction.absolute_energy_j / 8.0;
    let step = flow
        .step_with_piston_and_contact_equilibrium(
            &mut piston,
            duration,
            StepControl {
                maximum_halvings: 12,
                ..no_retry
            },
            load,
            controls(),
            equilibrium(budget),
            &mut scratch,
            || false,
        )
        .unwrap();
    assert!(step.gas.transport.halvings > 0);
    assert!(!step.gas.transport.first_order_retry);
    assert!(step.gas.correction.absolute_energy_j <= budget);
    assert!(step.admitted(controls()));
    for cutoff in [1, polls / 2, polls] {
        let mut flow = contact_field(16, [0.0; 3]);
        let mut piston = initial;
        let mut scratch = workspace(16);
        let mut seen = 0;
        let result = flow.step_with_piston_and_contact_equilibrium(
            &mut piston,
            duration,
            no_retry,
            load,
            controls(),
            equilibrium(1.0),
            &mut scratch,
            || {
                seen += 1;
                seen >= cutoff
            },
        );
        assert_eq!(result, Err(PistonError::Cancelled));
        assert_eq!(seen, cutoff);
        assert_eq!(flow.states, state);
        assert_eq!(flow.mesh(), &mesh);
        assert_eq!(flow.time().value(), 0.0);
        assert_eq!(piston, initial);
        let replay = flow
            .step_with_piston_and_contact_equilibrium(
                &mut piston,
                duration,
                no_retry,
                load,
                controls(),
                equilibrium(1.0),
                &mut scratch,
                || false,
            )
            .unwrap();
        close(
            replay.gas.correction.signed_energy_j,
            full.gas.correction.signed_energy_j,
            3e-12,
        );
        close(piston.velocity_m_s(), baseline_piston.velocity_m_s(), 3e-12);
        assert_eq!(flow.mesh(), baseline.mesh());
        for (&q, &expected) in flow.states.iter().zip(&baseline.states) {
            close(q, expected, 3e-12);
        }
    }
}
