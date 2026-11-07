//! Independent graded-depth, force, spring and native-continuation references.

use super::*;
use fs_contact::normal_patch::{FiniteGapPlasticDrive, FiniteGapPlasticLoadControls};

fn controls() -> FiniteGapPlasticLoadControls {
    FiniteGapPlasticLoadControls {
        minimum_approach_m: -1.0e-4,
        maximum_approach_m: 1.0e-4,
        absolute_force_tolerance_n: 1.0e-14,
        relative_force_tolerance: 1.0e-12,
        approach_tolerance_m: 1.0e-18,
        maximum_equilibrium_evaluations: 128,
    }
}

fn compliance(source: &SourceDepthHardnessWear) -> f64 {
    let config = source.model().config();
    4.0 * config.grid.cell_width_m * (1.0 + 2.0_f64.sqrt()).ln()
        / (core::f64::consts::PI * config.reduced_modulus_pa)
}

fn local_slip(distance: f64) -> Vec<f64> {
    let mut slip = vec![0.0; 25];
    slip[12] = distance;
    slip
}

// Independent positive root of H0*z + gradient*z²/2 = k*p*s.
fn affine_depth(pressure: f64, distance: f64) -> f64 {
    let integral = K[0] * pressure * distance;
    2.0 * integral / (H0 + (H0 * H0 + 2.0 * GRADIENT * integral).sqrt())
}

#[test]
fn g1_integrated_load_wear_matches_affine_primitive_and_actuator_work() {
    with_cx(|cx, _| {
        let fixture = Fixture::new(GRADIENT);
        let initial = fixture
            .owner(false, FiniteGapPlasticLinearSolver::Dense, 1.0e-4, cx)
            .advance_normal(5.0e-5, cx)
            .unwrap()
            .candidate;
        let force = 0.005;
        let pressure = force / AREA;
        let loaded = initial
            .advance_normal_at_load(force, controls(), cx)
            .unwrap();
        let initial_energy = loaded.normal.receipt.normal.reversible_energy_j;
        let before = loaded.normal.candidate;
        let saved = before.clone();
        let distance = 1.0e-4;
        let result = before
            .advance_wear_integrated_at_load(&local_slip(distance), force, controls(), cx)
            .unwrap();
        let wear = &result.wear;
        let next = wear.candidate.contact_state();
        let depths = [
            affine_depth(pressure, distance),
            K[1] * pressure * distance / HB,
        ];
        assert!(depths[0] > K[0] * pressure * distance / H0);
        for (side, depth) in depths.into_iter().enumerate() {
            close(wear.recession_increment_m[side][12], depth);
            close(wear.candidate.recession_m()[side][12], depth);
            close(wear.wear[side].increment_m3, AREA * depth);
            close(wear.candidate.volumes_m3()[side], AREA * depth);
        }
        let relative = depths[0] + depths[1];
        close(result.recession_load.achieved_force_n, force);
        close(result.hardness_load.achieved_force_n, force);
        close(next.pressure_pa()[12], pressure);
        close(
            next.approach_m() - before.contact_state().approach_m(),
            relative,
        );
        close(
            wear.recession.normal_work_j + wear.hardness.normal_work_j,
            force * relative,
        );
        close(
            wear.recession.parameter_work_j + wear.hardness.parameter_work_j,
            -force * relative,
        );
        close(wear.hardness.elastic_energy_j, initial_energy);
        assert_eq!(
            next.plastic_depth_m(),
            before.contact_state().plastic_depth_m()
        );
        assert_eq!(wear.recession.response.plastic_work_j, 0.0);
        assert_eq!(wear.hardness.response.plastic_work_j, 0.0);
        close(
            wear.candidate.model().config().hardness_pa[12],
            H0 + GRADIENT * depths[0],
        );
        assert_eq!(
            wear.hardness.predecessor_state_identity,
            wear.recession.successor_state_identity
        );
        assert_eq!(
            result.hardness_load.candidate_state_identity,
            next.identity()
        );
        assert_eq!(before, saved);

        // Pressure is fixed in this one-cell elastic regime, so exact profile
        // integration composes across unequal physical slip subdivisions.
        let first = before
            .advance_wear_integrated_at_load(&local_slip(0.3 * distance), force, controls(), cx)
            .unwrap();
        let split = first
            .wear
            .candidate
            .advance_wear_integrated_at_load(&local_slip(0.7 * distance), force, controls(), cx)
            .unwrap();
        for side in 0..2 {
            close(split.wear.candidate.recession_m()[side][12], depths[side]);
            close(
                split.wear.candidate.volumes_m3()[side],
                wear.candidate.volumes_m3()[side],
            );
        }
        let unloaded = wear
            .candidate
            .advance_normal_at_load(0.0025, controls(), cx)
            .unwrap();
        assert_eq!(
            unloaded.normal.candidate.recession_m(),
            wear.candidate.recession_m()
        );
        assert_eq!(
            unloaded.normal.candidate.volumes_m3(),
            wear.candidate.volumes_m3()
        );
    });
}

#[test]
fn g1_integrated_spring_wear_resolves_force_relaxation_softening_and_storage() {
    with_cx(|cx, _| {
        let fixture = Fixture::new(GRADIENT);
        for (initial_pressure, distance, initial_approach) in
            [(0.5e6, 1.0e-4, None), (H0, 1.0e-5, Some(5.0e-5))]
        {
            let pristine = fixture.owner(false, FiniteGapPlasticLinearSolver::Dense, 1.0e-4, cx);
            let c = compliance(&pristine);
            let drive = FiniteGapPlasticDrive {
                stiffness_n_per_m: AREA / c,
                position_m: initial_approach.unwrap_or(c * initial_pressure) + c * initial_pressure,
            };
            let loaded = pristine
                .advance_normal_with_drive(drive, controls(), cx)
                .unwrap();
            let previous_spring_energy = loaded.drive.spring_energy_j;
            let previous_contact_energy = loaded.normal.receipt.normal.reversible_energy_j;
            let before = loaded.normal.candidate;
            let saved = before.clone();
            close(before.contact_state().pressure_pa()[12], initial_pressure);
            let result = before
                .advance_wear_integrated_with_drive(&local_slip(distance), drive, controls(), cx)
                .unwrap();
            let wear = &result.wear;
            let next = wear.candidate.contact_state();
            let depths = [
                affine_depth(initial_pressure, distance),
                K[1] * initial_pressure * distance / HB,
            ];
            let relative = depths[0] + depths[1];
            let cap = H0 + GRADIENT * depths[0];
            let old_plastic = before.contact_state().plastic_depth_m()[12];
            let elastic_pressure = (drive.position_m - old_plastic - relative) / (2.0 * c);
            let expected_pressure = elastic_pressure.min(cap);
            let expected_approach = drive.position_m - c * expected_pressure;
            close(next.pressure_pa()[12], expected_pressure);
            close(next.approach_m(), expected_approach);
            close(wear.candidate.model().config().hardness_pa[12], cap);
            if elastic_pressure <= cap {
                // The elastic branch retains the old irreversible depth
                // exactly; subtracting its three displacements only measures
                // floating cancellation around that zero increment.
                assert_eq!(next.plastic_depth_m()[12], old_plastic);
            } else {
                close(
                    next.plastic_depth_m()[12],
                    expected_approach - relative - c * expected_pressure,
                );
            }
            close(
                result.hardness_drive.achieved_force_n,
                AREA * expected_pressure,
            );
            close(
                result.hardness_drive.drive_force_n,
                AREA * expected_pressure,
            );
            for (side, depth) in depths.into_iter().enumerate() {
                close(wear.candidate.recession_m()[side][12], depth);
                close(wear.candidate.volumes_m3()[side], AREA * depth);
            }
            let released = previous_spring_energy - result.hardness_drive.spring_energy_j;
            close(
                wear.recession.normal_work_j + wear.hardness.normal_work_j,
                released,
            );
            let plastic_work =
                wear.recession.response.plastic_work_j + wear.hardness.response.plastic_work_j;
            close(
                wear.hardness.elastic_energy_j - previous_contact_energy + plastic_work,
                released + wear.recession.parameter_work_j + wear.hardness.parameter_work_j,
            );
            if initial_approach.is_some() {
                assert!(wear.hardness.response.plastic_work_j > 0.0);
                assert!(next.plastic_depth_m()[12] > old_plastic);
            } else {
                assert_eq!(
                    next.plastic_depth_m(),
                    before.contact_state().plastic_depth_m()
                );
            }
            assert_eq!(before, saved);
            let retained = wear
                .candidate
                .advance_normal_with_drive(drive, controls(), cx)
                .unwrap();
            assert_eq!(retained.normal.candidate, wear.candidate);
        }
    });
}

#[test]
fn g3_integrated_constant_profiles_recover_native_boundaries_in_all_lanes() {
    with_cx(|cx, _| {
        let fixture = Fixture::new(0.0);
        for solver in [
            FiniteGapPlasticLinearSolver::Dense,
            FiniteGapPlasticLinearSolver::MatrixFree {
                maximum_iterations: 256,
                maximum_corrections: 3,
            },
            FiniteGapPlasticLinearSolver::Fft {
                maximum_iterations: 256,
                maximum_corrections: 3,
            },
        ] {
            let pristine = fixture.owner(false, solver, 1.0e-4, cx);
            let force = 0.005;
            let loaded = pristine
                .advance_normal_at_load(force, controls(), cx)
                .unwrap()
                .normal
                .candidate;
            let slip = local_slip(1.0e-5);
            let held = loaded.advance_wear_integrated(&slip, cx).unwrap();
            let explicit = loaded.advance_wear(&slip, cx).unwrap();
            close(
                held.candidate.contact_state().pressure_pa()[12],
                explicit.candidate.contact_state().pressure_pa()[12],
            );
            close(
                held.candidate.recession_m()[0][12],
                explicit.candidate.recession_m()[0][12],
            );
            close(held.wear[0].increment_m3, explicit.wear[0].increment_m3);
            let controlled = loaded
                .advance_wear_integrated_at_load(&slip, force, controls(), cx)
                .unwrap();
            let relative: Vec<_> = controlled.wear.recession_increment_m[0]
                .iter()
                .zip(&controlled.wear.recession_increment_m[1])
                .map(|(a, b)| a + b)
                .collect();
            let native = loaded
                .model()
                .recede_at_load(loaded.contact_state(), &relative, force, controls(), cx)
                .unwrap();
            assert_eq!(
                controlled.wear.candidate.contact_state(),
                &native.continuation.candidate
            );
            assert_eq!(controlled.wear.recession, native.continuation.receipt);
            let zero = controlled
                .wear
                .candidate
                .advance_wear_integrated_at_load(&[0.0; 25], force, controls(), cx)
                .unwrap();
            assert_eq!(zero.wear.candidate, controlled.wear.candidate);

            let c = compliance(&pristine);
            let drive = FiniteGapPlasticDrive {
                stiffness_n_per_m: AREA / c,
                position_m: 2.0 * c * force / AREA,
            };
            let spring = pristine
                .advance_normal_with_drive(drive, controls(), cx)
                .unwrap()
                .normal
                .candidate;
            let controlled = spring
                .advance_wear_integrated_with_drive(&slip, drive, controls(), cx)
                .unwrap();
            let relative: Vec<_> = controlled.wear.recession_increment_m[0]
                .iter()
                .zip(&controlled.wear.recession_increment_m[1])
                .map(|(a, b)| a + b)
                .collect();
            let native = spring
                .model()
                .recede_with_drive(spring.contact_state(), &relative, drive, controls(), cx)
                .unwrap();
            assert_eq!(
                controlled.wear.candidate.contact_state(),
                &native.continuation.candidate
            );
            assert_eq!(controlled.wear.recession, native.continuation.receipt);
            assert_eq!(controlled.wear.hardness.response.plastic_work_j, 0.0);
        }
    });
}

#[test]
fn g4_integrated_control_refusals_preserve_both_surfaces_and_mechanical_history() {
    with_cx(|cx, gate| {
        let fixture = Fixture::new(GRADIENT);
        let pristine = fixture.owner(false, FiniteGapPlasticLinearSolver::Dense, 1.0e-4, cx);
        let force = 0.009;
        let before = pristine
            .advance_normal_at_load(force, controls(), cx)
            .unwrap()
            .normal
            .candidate;
        let saved = before.clone();
        // Both depth integrals are supported, and recession at the original
        // caps can carry the force. Only final exposure of the softer layer
        // makes the maintained force unattainable inside the original bounds.
        let error = before
            .advance_wear_integrated_at_load(&local_slip(1.0e-4), force, controls(), cx)
            .unwrap_err();
        assert!(
            matches!(error, DepthHardnessWearError::Load(_)),
            "{error:?}"
        );
        assert_eq!(before, saved);
        for slip in [
            local_slip(1.0),
            local_slip(f64::from_bits(1)),
            local_slip(-1.0),
        ] {
            assert!(
                before
                    .advance_wear_integrated_at_load(&slip, force, controls(), cx)
                    .is_err()
            );
            assert_eq!(before, saved);
        }
        assert!(
            before
                .advance_wear_integrated_at_load(&[0.0; 25], 0.005, controls(), cx)
                .is_err()
        );
        let mut outside = controls();
        outside.maximum_approach_m = 2.0e-4;
        assert!(
            before
                .advance_wear_integrated_at_load(&[0.0; 25], force, outside, cx)
                .is_err()
        );
        let c = compliance(&before);
        let mismatched = FiniteGapPlasticDrive {
            stiffness_n_per_m: AREA / c,
            position_m: 0.0,
        };
        assert!(
            before
                .advance_wear_integrated_with_drive(&[0.0; 25], mismatched, controls(), cx)
                .is_err()
        );
        assert_eq!(before, saved);
        gate.request();
        assert!(before.advance_wear_integrated(&[0.0; 25], cx).is_err());
        assert!(
            before
                .advance_wear_integrated_at_load(&[0.0; 25], force, controls(), cx)
                .is_err()
        );
        assert_eq!(before, saved);
    });
}
