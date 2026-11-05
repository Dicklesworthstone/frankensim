//! Two independent caloric contact fields in one pressure/inertia transaction.

use super::*;
use fs_qty::Temperature;

fn contact_gas(right: bool, datum: f64) -> MixtureGas {
    let models = if right {
        vec![
            source("right-A", 3.0 * R, 3.25, 0.25, datum),
            source("right-B", 0.75 * R, 5.0, 0.1, -2.0 * datum),
            source("right-absent-1", R, 3.5, 0.0, 3.0 * datum),
            source("right-absent-2", 2.0 * R, 4.0, 0.0, -4.0 * datum),
        ]
    } else {
        vec![
            source("left-A", R, 3.5, 0.4, datum),
            source("left-B", 2.0 * R, 4.5, 0.2, -2.0 * datum),
            source("left-absent", 0.5 * R, 3.0, 0.1, 3.0 * datum),
        ]
    };
    let mut domain =
        crate::dynamics::CaloricDomain::new(Temperature::new(0.001), Temperature::new(50.0));
    domain.absolute_energy_tolerance = 1e-13;
    domain.relative_energy_tolerance = 1e-14;
    MixtureGas::new(gas(models).source().clone(), domain, || false).unwrap()
}

fn contact_chamber(
    right: bool,
    length: f64,
    count: usize,
    pressure: f64,
    datum: f64,
) -> MixtureEuler1d {
    let gas = contact_gas(right, datum);
    let start = if right { 1.0 } else { 0.0 };
    let positions: Vec<_> = (0..=count)
        .map(|i| start + length * (i as f64 / count as f64))
        .collect();
    let volumes = positions.windows(2).map(|p| AREA * (p[1] - p[0])).collect();
    let mesh = Mesh1d::new(positions, vec![AREA; count + 1], volumes).unwrap();
    let states = (0..count)
        .map(|i| {
            let mut densities = vec![0.0; gas.species_count()];
            densities[usize::from(i >= count / 2)] = 1.0;
            gas.conserved(&densities, 0.0, pressure).unwrap()
        })
        .collect();
    MixtureEuler1d::new_on_mesh(gas, mesh, states, CLOSED, || false).unwrap()
}

fn contact_pair(pressure: [f64; 2], datum: [f64; 2]) -> (MixtureEuler1d, MixtureEuler1d) {
    (
        contact_chamber(false, 1.0, 64, pressure[0], datum[0]),
        contact_chamber(true, 1.0, 80, pressure[1], datum[1]),
    )
}

fn workspace(count: usize) -> ContactWorkspace {
    ContactWorkspace::new(
        count,
        ContactWorkspace::required_workspace_bytes(count).unwrap(),
        || false,
    )
    .unwrap()
}

fn equilibrium(budget: [f64; 2]) -> [ContactEquilibriumControl; 2] {
    budget.map(|maximum_energy_defect_j| ContactEquilibriumControl {
        relative_tolerance: 1e-10,
        maximum_energy_defect_j,
    })
}

#[test]
fn g0_independent_contact_pressures_survive_shared_motion_with_full_balances() {
    for direction in [-1.0, 1.0] {
        for reconstruction in [Reconstruction::FirstOrder, Reconstruction::Muscl] {
            let pressure = [1.0 + direction * 0.2, 1.0 - direction * 0.2];
            let (mut left, mut right) = contact_pair(pressure, [1e18, -2e18]);
            let (mut native_left, mut native_right) = contact_pair(pressure, [1e18, -2e18]);
            let mut piston = Piston::new(0.4, direction * 0.25).unwrap();
            let original = piston;
            let mut native_piston = piston;
            let before = [left.totals(), right.totals()];
            let (mut a, mut b) = (workspace(64), workspace(80));
            let bytes = [a.workspace_bytes(), b.workspace_bytes()];
            let gas_controls = StepControl {
                reconstruction,
                maximum_halvings: 0,
                ..StepControl::default()
            };
            let external = direction * 0.05;
            let (mut work, mut impulse, mut supports) = (0.0, 0.0, 0.0);
            let (mut signed, mut absolute) = ([0.0; 2], [0.0; 2]);
            for _ in 0..24 {
                native_left
                    .step_with_shared_piston(
                        &mut native_right,
                        &mut native_piston,
                        Time::new(0.001),
                        gas_controls,
                        external,
                        controls(),
                        || false,
                    )
                    .unwrap();
                let step = left
                    .step_with_shared_piston_and_contact_equilibrium(
                        &mut right,
                        &mut piston,
                        Time::new(0.001),
                        gas_controls,
                        external,
                        controls(),
                        equilibrium([1.0; 2]),
                        [&mut a, &mut b],
                        || false,
                    )
                    .unwrap();
                assert!(step.admitted(controls()));
                assert_eq!(step.exterior_pressure_work_j, 0.0);
                let reports = [step.gas.left, step.gas.right];
                assert_eq!(reports[0].transport.elapsed, reports[1].transport.elapsed);
                for (i, side) in reports.iter().enumerate() {
                    assert_eq!(side.transport.halvings, 0);
                    assert!(!side.transport.first_order_retry);
                    close(side.corrected_energy_residual_j, 0.0, 3e-12);
                    close(side.transport.conservation_residual.momentum, 0.0, 3e-12);
                    for mass in side.transport.boundary_change.species_masses {
                        assert_eq!(mass, 0.0);
                    }
                    for mass in side.transport.conservation_residual.species_masses {
                        close(mass, 0.0, 3e-12);
                    }
                    signed[i] += side.correction.signed_energy_j;
                    absolute[i] += side.correction.absolute_energy_j;
                }
                close(
                    reports[0].transport.right_boundary_change.energy
                        + reports[1].transport.left_boundary_change.energy,
                    -step.midpoint_velocity_m_s * step.gas_impulse_n_s,
                    3e-14,
                );
                let raw =
                    sum(reports
                        .map(|side| side.transport.after.energy - side.transport.before.energy))
                        + step.kinetic_energy_change_j
                        - step.external_force_work_j;
                let correction = sum(reports.map(|side| side.correction.signed_energy_j));
                close(raw - correction, step.total_energy_residual_j, 3e-14);
                work += step.external_force_work_j;
                impulse += external * reports[0].transport.elapsed.value();
                supports += step.fixed_wall_impulse_n_s;
                assert_eq!(left.time(), right.time());
                assert_eq!(left.mesh().positions()[64], right.mesh().positions()[0]);
            }
            for (i, flow, native) in [(0, &left, &native_left), (1, &right, &native_right)] {
                let (mut native_error, mut corrected_error) = (0.0_f64, 0.0_f64);
                for cell in flow.cell_count() / 2 - 6..flow.cell_count() / 2 + 6 {
                    let old = native
                        .gas()
                        .primitive(native.cell(cell).unwrap(), || false)
                        .unwrap();
                    let new = flow
                        .gas()
                        .primitive(flow.cell(cell).unwrap(), || false)
                        .unwrap();
                    native_error = native_error.max((old.pressure - pressure[i]).abs());
                    corrected_error = corrected_error.max((new.pressure - pressure[i]).abs());
                    close(new.velocity, 0.0, 3e-10);
                }
                assert!(
                    native_error > 1e-4,
                    "side {i} native error {native_error:e}"
                );
                assert!(
                    corrected_error < 3e-10,
                    "side {i} corrected error {corrected_error:e}"
                );
                assert!(absolute[i] > 1e-7);
                for cell in 0..flow.cell_count() {
                    assert!(
                        flow.cell(cell).unwrap().partial_densities[2..]
                            .iter()
                            .all(|rho| *rho == 0.0)
                    );
                }
                let after = flow.totals();
                for k in 0..flow.gas().species_count() {
                    close(after.species_masses[k], before[i].species_masses[k], 3e-12);
                }
            }
            let after = [left.totals(), right.totals()];
            let kinetic = 0.5
                * piston.mass_kg()
                * (piston.velocity_m_s().powi(2) - original.velocity_m_s().powi(2));
            close(
                sum([
                    after[0].energy - before[0].energy,
                    after[1].energy - before[1].energy,
                ]) + kinetic
                    - work,
                sum(signed),
                4e-10,
            );
            close(
                sum([
                    after[0].momentum - before[0].momentum,
                    after[1].momentum - before[1].momentum,
                ]) + piston.mass_kg() * (piston.velocity_m_s() - original.velocity_m_s())
                    + supports,
                impulse,
                4e-10,
            );
            close(
                left.mesh().volumes().iter().sum::<f64>()
                    + right.mesh().volumes().iter().sum::<f64>(),
                2.0 * AREA,
                4e-14,
            );
            assert_eq!([a.workspace_bytes(), b.workspace_bytes()], bytes);
        }
    }
}

#[test]
fn g0_fine_side_cfl_final_trial_replay_and_zero_suspension_are_shared() {
    let make = || {
        (
            contact_chamber(false, 1.0, 16, 1.3, 0.0),
            contact_chamber(true, 0.125, 64, 0.7, 0.0),
        )
    };
    let (mut left, mut right) = make();
    let (mut replay_left, mut replay_right) = make();
    let (mut zero_left, mut zero_right) = make();
    let mut piston = Piston::new(0.4, 0.125).unwrap();
    let mut other = piston;
    let step = left
        .step_with_shared_piston_and_contact_equilibrium(
            &mut right,
            &mut piston,
            Time::new(0.02),
            StepControl::default(),
            0.02,
            controls(),
            equilibrium([1.0; 2]),
            [&mut workspace(16), &mut workspace(64)],
            || false,
        )
        .unwrap();
    assert!(step.gas_trials > 1);
    assert!(step.gas.left.transport.elapsed.value() < 0.001);
    assert_eq!(
        step.gas.left.transport.elapsed,
        step.gas.right.transport.elapsed
    );
    for (replay, accepted, side, speeds) in [
        (
            &mut replay_left,
            &left,
            step.gas.left,
            [0.0, step.midpoint_velocity_m_s],
        ),
        (
            &mut replay_right,
            &right,
            step.gas.right,
            [step.midpoint_velocity_m_s, 0.0],
        ),
    ] {
        let expected = replay
            .step_with_moving_walls_and_contact_equilibrium(
                side.transport.elapsed,
                StepControl::default(),
                speeds,
                equilibrium([1.0; 2])[0],
                &mut workspace(replay.cell_count()),
                || false,
            )
            .unwrap();
        assert_eq!(replay.mesh(), accepted.mesh());
        assert_eq!(replay.time(), accepted.time());
        for (&actual, &reference) in accepted.states.iter().zip(&replay.states) {
            close(actual, reference, 4e-12);
        }
        close(
            side.correction.signed_energy_j,
            expected.correction.signed_energy_j,
            4e-12,
        );
        close(
            side.correction.absolute_energy_j,
            expected.correction.absolute_energy_j,
            4e-12,
        );
        assert_eq!(
            side.correction.corrected_cells,
            expected.correction.corrected_cells
        );
        close(
            side.transport.boundary_change.energy,
            expected.transport.boundary_change.energy,
            4e-12,
        );
    }
    let zero = zero_left
        .step_with_shared_suspended_piston_and_contact_equilibrium(
            &mut zero_right,
            &mut other,
            Time::new(0.02),
            StepControl::default(),
            0.02,
            PistonSuspension::new(0.0, 0.0, -3.0).unwrap(),
            controls(),
            equilibrium([1.0; 2]),
            [&mut workspace(16), &mut workspace(64)],
            || false,
        )
        .unwrap();
    assert_eq!(zero.piston, step);
    assert_eq!(other, piston);
    assert_eq!(zero_left.states, left.states);
    assert_eq!(zero_right.states, right.states);
}

#[test]
fn g1_shared_suspension_closes_both_gases_storage_loss_and_independent_datums() {
    let make = |datums: [f64; 2]| {
        (
            contact_chamber(false, 1.0, 24, 1.2, datums[0]),
            contact_chamber(true, 1.0, 40, 0.8, datums[1]),
        )
    };
    let (mut left, mut right) = make([0.0; 2]);
    let (mut shifted_left, mut shifted_right) = make([1e18, -2e18]);
    let mut piston = Piston::new(0.7, 0.125).unwrap();
    let mut other = piston;
    let spring = PistonSuspension::new(8.0, 1.2, 0.92).unwrap();
    let initial = left.totals().energy
        + right.totals().energy
        + 0.5 * 0.7 * 0.125_f64.powi(2)
        + spring.stored_energy_j(1.0).unwrap();
    let (mut heat, mut work, mut defect) = (0.0, 0.0, 0.0);
    let (mut a, mut b) = (workspace(24), workspace(40));
    let (mut c, mut d) = (workspace(24), workspace(40));
    for _ in 0..32 {
        let step = left
            .step_with_shared_suspended_piston_and_contact_equilibrium(
                &mut right,
                &mut piston,
                Time::new(0.001),
                StepControl::default(),
                0.03,
                spring,
                controls(),
                equilibrium([1.0; 2]),
                [&mut a, &mut b],
                || false,
            )
            .unwrap();
        let shifted = shifted_left
            .step_with_shared_suspended_piston_and_contact_equilibrium(
                &mut shifted_right,
                &mut other,
                Time::new(0.001),
                StepControl::default(),
                0.03,
                spring,
                controls(),
                equilibrium([1.0; 2]),
                [&mut c, &mut d],
                || false,
            )
            .unwrap();
        assert_eq!(step, shifted);
        assert_eq!(piston, other);
        assert_eq!(left.states, shifted_left.states);
        assert_eq!(right.states, shifted_right.states);
        assert!(step.total_energy_residual_j.abs() <= controls().energy_tolerance_j);
        assert_eq!(step.piston.exterior_pressure_work_j, 0.0);
        heat += step.damping_energy_j;
        work += step.applied_force_work_j;
        defect += sum([
            step.piston.gas.left.correction.signed_energy_j,
            step.piston.gas.right.correction.signed_energy_j,
        ]);
    }
    let final_energy = left.totals().energy
        + right.totals().energy
        + 0.5 * 0.7 * piston.velocity_m_s().powi(2)
        + spring.stored_energy_j(left.mesh().positions()[24]).unwrap();
    close(final_energy - initial + heat - work, defect, 5e-10);
    assert!(heat > 1e-7);
}

#[test]
fn g0_resolved_gradients_preserve_the_native_shared_suspension_map() {
    for pressure_gradient in [false, true] {
        let make = || {
            let mut pair = (
                contact_chamber(false, 1.0, 12, 1.2, 0.0),
                contact_chamber(true, 1.0, 16, 0.8, 0.0),
            );
            for flow in [&mut pair.0, &mut pair.1] {
                let n = flow.gas().species_count();
                for i in 0..flow.cell_count() {
                    let state = flow
                        .gas()
                        .conserved(
                            flow.cell(i).unwrap().partial_densities,
                            if pressure_gradient {
                                0.1
                            } else {
                                0.02 * i as f64
                            },
                            if pressure_gradient {
                                1.0 + 0.02 * i as f64
                            } else {
                                1.0
                            },
                        )
                        .unwrap();
                    let target = &mut flow.states[i * (n + 2)..(i + 1) * (n + 2)];
                    target[..n].copy_from_slice(&state.partial_densities);
                    target[n] = state.momentum;
                    target[n + 1] = state.energy;
                }
            }
            pair
        };
        let (mut left, mut right) = make();
        let (mut native_left, mut native_right) = make();
        let mut piston = Piston::new(0.4, 0.125).unwrap();
        let mut other = piston;
        let spring = PistonSuspension::new(2.0, 0.3, 1.04).unwrap();
        let expected = native_left
            .step_with_shared_suspended_piston(
                &mut native_right,
                &mut other,
                Time::new(0.001),
                StepControl::default(),
                0.05,
                spring,
                controls(),
                || false,
            )
            .unwrap();
        let actual = left
            .step_with_shared_suspended_piston_and_contact_equilibrium(
                &mut right,
                &mut piston,
                Time::new(0.001),
                StepControl::default(),
                0.05,
                spring,
                controls(),
                equilibrium([0.0; 2]),
                [&mut workspace(12), &mut workspace(16)],
                || false,
            )
            .unwrap();
        assert_eq!(
            actual.piston.gas.left.correction,
            ContactCorrection::default()
        );
        assert_eq!(
            actual.piston.gas.right.correction,
            ContactCorrection::default()
        );
        assert_eq!(
            actual.map_gas(|gas| SharedPistonGasReport {
                left: gas.left.transport,
                right: gas.right.transport,
            }),
            expected
        );
        assert_eq!(piston, other);
        assert_eq!(left.states, native_left.states);
        assert_eq!(right.states, native_right.states);
        assert_eq!(left.mesh(), native_left.mesh());
        assert_eq!(right.mesh(), native_right.mesh());
    }
}

#[test]
fn g4_each_chambers_defect_budget_and_late_cancellation_are_jointly_atomic() {
    let make = || {
        (
            contact_chamber(false, 1.0, 16, 1.2, 0.0),
            contact_chamber(true, 1.0, 24, 0.8, 0.0),
        )
    };
    let no_retry = StepControl {
        maximum_halvings: 0,
        ..StepControl::default()
    };
    let duration = Time::new(0.003);
    let initial = Piston::new(0.4, 0.125).unwrap();
    let (mut baseline_left, mut baseline_right) = make();
    let mut baseline_piston = initial;
    let mut polls = 0;
    let full = baseline_left
        .step_with_shared_piston_and_contact_equilibrium(
            &mut baseline_right,
            &mut baseline_piston,
            duration,
            no_retry,
            0.05,
            controls(),
            equilibrium([1.0; 2]),
            [&mut workspace(16), &mut workspace(24)],
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    assert!(
        full.gas.left.correction.signed_energy_j * full.gas.right.correction.signed_energy_j < 0.0
    );
    let required = [
        full.gas.left.correction.absolute_energy_j,
        full.gas.right.correction.absolute_energy_j,
    ];
    for limited in 0..2 {
        let (mut left, mut right) = make();
        let states = [left.states.clone(), right.states.clone()];
        let meshes = [left.mesh().clone(), right.mesh().clone()];
        let mut piston = initial;
        let (mut a, mut b) = (workspace(16), workspace(24));
        let mut budget = [1.0; 2];
        budget[limited] = 0.0;
        let refused = left.step_with_shared_piston_and_contact_equilibrium(
            &mut right,
            &mut piston,
            duration,
            no_retry,
            0.05,
            controls(),
            equilibrium(budget),
            [&mut a, &mut b],
            || false,
        );
        assert!(matches!(refused, Err(PistonError::StepRetries { last, .. })
            if matches!(*last, PistonError::MixtureGas(Error::ContactEnergyBudget { .. }))));
        assert_eq!(left.states, states[0]);
        assert_eq!(right.states, states[1]);
        assert_eq!(left.mesh(), &meshes[0]);
        assert_eq!(right.mesh(), &meshes[1]);
        assert_eq!([left.time().value(), right.time().value()], [0.0; 2]);
        assert_eq!(piston, initial);
        budget[limited] = required[limited] / 8.0;
        let retry = left
            .step_with_shared_piston_and_contact_equilibrium(
                &mut right,
                &mut piston,
                duration,
                StepControl {
                    maximum_halvings: 12,
                    ..no_retry
                },
                0.05,
                controls(),
                equilibrium(budget),
                [&mut a, &mut b],
                || false,
            )
            .unwrap();
        assert!(retry.admitted(controls()));
        assert!(retry.gas.left.transport.halvings > 0);
        assert_eq!(
            retry.gas.left.transport.halvings,
            retry.gas.right.transport.halvings
        );
        assert_eq!(
            retry.gas.left.transport.elapsed,
            retry.gas.right.transport.elapsed
        );
        for (i, side) in [retry.gas.left, retry.gas.right].iter().enumerate() {
            assert!(!side.transport.first_order_retry);
            assert!(side.correction.absolute_energy_j <= budget[i]);
        }
    }
    for cutoff in [1, polls / 2, polls] {
        let (mut left, mut right) = make();
        let states = [left.states.clone(), right.states.clone()];
        let meshes = [left.mesh().clone(), right.mesh().clone()];
        let mut piston = initial;
        let (mut a, mut b) = (workspace(16), workspace(24));
        let mut seen = 0;
        let cancelled = left.step_with_shared_piston_and_contact_equilibrium(
            &mut right,
            &mut piston,
            duration,
            no_retry,
            0.05,
            controls(),
            equilibrium([1.0; 2]),
            [&mut a, &mut b],
            || {
                seen += 1;
                seen >= cutoff
            },
        );
        assert_eq!(cancelled, Err(PistonError::Cancelled));
        assert_eq!(seen, cutoff);
        assert_eq!(left.states, states[0]);
        assert_eq!(right.states, states[1]);
        assert_eq!(left.mesh(), &meshes[0]);
        assert_eq!(right.mesh(), &meshes[1]);
        assert_eq!([left.time().value(), right.time().value()], [0.0; 2]);
        assert_eq!(piston, initial);
        let replay = left
            .step_with_shared_piston_and_contact_equilibrium(
                &mut right,
                &mut piston,
                duration,
                no_retry,
                0.05,
                controls(),
                equilibrium([1.0; 2]),
                [&mut a, &mut b],
                || false,
            )
            .unwrap();
        close(piston.velocity_m_s(), baseline_piston.velocity_m_s(), 4e-12);
        for (actual, expected) in [
            (replay.gas.left, full.gas.left),
            (replay.gas.right, full.gas.right),
        ] {
            close(
                actual.correction.signed_energy_j,
                expected.correction.signed_energy_j,
                4e-12,
            );
            close(
                actual.correction.absolute_energy_j,
                expected.correction.absolute_energy_j,
                4e-12,
            );
        }
        for (actual, expected) in [(&left, &baseline_left), (&right, &baseline_right)] {
            assert_eq!(actual.mesh(), expected.mesh());
            for (&a, &b) in actual.states.iter().zip(&expected.states) {
                close(a, b, 4e-12);
            }
        }
    }
}
