//! Independent maintained-force Archard and external-work references.

use super::*;
use fs_contact::normal_patch::FiniteGapPlasticLoadControls;
use fs_contact::normal_patch::plastic_source::SourcePlasticContactError;

fn load_controls() -> FiniteGapPlasticLoadControls {
    FiniteGapPlasticLoadControls {
        minimum_approach_m: -1.0e-4,
        maximum_approach_m: 1.0e-4,
        absolute_force_tolerance_n: 1.0e-14,
        relative_force_tolerance: 1.0e-12,
        approach_tolerance_m: 1.0e-18,
        maximum_equilibrium_evaluations: 128,
    }
}

#[test]
fn maintained_load_matches_independent_archard_volume_approach_and_work() {
    with_cx(|cx, _| {
        let model = fixture([0.01, 0.03], [1.0e6, 2.0e6], [1.0e-4; 2], false, cx);
        let dented = model
            .advance_normal(&model.initial_state(), 5.0e-5, cx)
            .unwrap()
            .candidate;
        let loaded = model
            .advance_normal_at_load(&dented, 0.005, load_controls(), cx)
            .unwrap();
        let previous_energy = loaded.receipt.normal.reversible_energy_j;
        let before = loaded.candidate;
        let saved = before.clone();
        let previous = before.contact_state().contact_state();
        assert!(previous.plastic_depth_m()[12] > 0.0);

        // One loaded cell has p=F/A=0.5 MPa, below the 1 MPa pressure cap.
        // For s=1e-4 m, each body independently removes V_i=k_i*F*s/H_i.
        let result = model
            .advance_wear_at_load(&before, &[1.0e-4; 25], 0.005, load_controls(), cx)
            .unwrap();
        let wear = &result.wear;
        let next = wear.candidate.contact_state().contact_state();
        assert_eq!(wear.source_identity, model.identity());
        close(result.load.achieved_force_n, 0.005);
        close(wear.receipt.response.normal.normal_force_n, 0.005);
        close(wear.pressure_weighted_sliding_j, 0.005 * 1.0e-4);
        for side in 0..2 {
            let volume = [5.0e-15, 7.5e-15][side];
            let depth = [0.5e-6, 0.75e-6][side];
            close(wear.wear[side].increment_m3, volume);
            close(wear.candidate.volumes_m3()[side], volume);
            close(wear.recession_increment_m[side][12], depth);
            close(wear.candidate.recession_m()[side][12], depth);
            assert!(
                wear.candidate.recession_m()[side]
                    .iter()
                    .enumerate()
                    .all(|(index, depth)| index == 12 || *depth == 0.0)
            );
        }
        close(next.recession_m()[12], 1.25e-6);
        close(next.approach_m() - previous.approach_m(), 1.25e-6);
        close(wear.receipt.previous_approach_m, previous.approach_m());
        close(wear.receipt.normal_work_j, 6.25e-9);
        close(wear.receipt.parameter_work_j, -6.25e-9);
        close(wear.receipt.elastic_energy_j, previous_energy);
        assert!(wear.receipt.elastic_energy_change_j.abs() <= 1.0e-8 * 6.25e-9);
        assert_eq!(wear.receipt.response.plastic_work_j, 0.0);
        assert_eq!(next.plastic_depth_m(), previous.plastic_depth_m());
        assert_eq!(
            next.cumulative_plastic_work_j(),
            previous.cumulative_plastic_work_j()
        );
        assert_eq!(before, saved);
    });
}

#[test]
fn explicit_load_changes_retain_existing_inventories_and_allow_further_wear() {
    with_cx(|cx, _| {
        let model = fixture([0.01, 0.03], [1.0e6, 2.0e6], [1.0e-4; 2], false, cx);
        let dented = model
            .advance_normal(&model.initial_state(), 5.0e-5, cx)
            .unwrap()
            .candidate;
        let first = model
            .advance_normal_at_load(&dented, 0.005, load_controls(), cx)
            .unwrap()
            .candidate;
        let worn = model
            .advance_wear_at_load(&first, &[1.0e-4; 25], 0.005, load_controls(), cx)
            .unwrap()
            .wear
            .candidate;
        let saved = worn.clone();
        let before = worn.contact_state().contact_state();
        let changed = model
            .advance_normal_at_load(&worn, 0.0025, load_controls(), cx)
            .unwrap();
        let next = changed.candidate.contact_state().contact_state();
        close(changed.receipt.normal.normal_force_n, 0.0025);
        assert_eq!(changed.candidate.volumes_m3(), worn.volumes_m3());
        assert_eq!(changed.candidate.recession_m(), worn.recession_m());
        assert_eq!(next.recession_m(), before.recession_m());
        assert_eq!(next.plastic_depth_m(), before.plastic_depth_m());
        assert_eq!(
            next.cumulative_plastic_work_j(),
            before.cumulative_plastic_work_j()
        );
        let retained_gap = before.plastic_depth_m()[12] + before.recession_m()[12];
        close(
            next.approach_m(),
            retained_gap + 0.5 * (before.approach_m() - retained_gap),
        );

        let continued = model
            .advance_wear_at_load(
                &changed.candidate,
                &[2.0e-4; 25],
                0.0025,
                load_controls(),
                cx,
            )
            .unwrap();
        let measure = 0.005 * 1.0e-4 + 0.0025 * 2.0e-4;
        for side in 0..2 {
            let expected = [0.01, 0.03][side] * measure / [1.0e6, 2.0e6][side];
            close(continued.wear.candidate.volumes_m3()[side], expected);
            close(
                continued.wear.candidate.recession_m()[side][12],
                expected / 1.0e-8,
            );
        }
        assert_eq!(
            continued
                .wear
                .candidate
                .contact_state()
                .contact_state()
                .plastic_depth_m(),
            before.plastic_depth_m()
        );
        assert_eq!(worn, saved);
    });
}

#[test]
fn maintained_load_mismatch_late_bounds_and_cancellation_are_atomic() {
    with_cx(|cx, gate| {
        let model = fixture([0.01, 0.03], [1.0e6, 2.0e6], [1.0e-4; 2], false, cx);
        let before = model
            .advance_normal_at_load(&model.initial_state(), 0.005, load_controls(), cx)
            .unwrap()
            .candidate;
        let saved = before.clone();
        let slip = [1.0e-4; 25];
        // Changing applied load must be explicit, even when both forces lie
        // safely below the hardness plateau.
        assert!(matches!(
            model.advance_wear_at_load(&before, &slip, 0.0075, load_controls(), cx),
            Err(SourcePlasticWearError::Load(_))
        ));
        let invalid_bracket = FiniteGapPlasticLoadControls {
            maximum_approach_m: 1.01e-4,
            ..load_controls()
        };
        assert!(matches!(
            model.advance_wear_at_load(&before, &slip, 0.005, invalid_bracket, cx),
            Err(SourcePlasticWearError::Contact(
                SourcePlasticContactError::OutsideApplicability { .. }
            ))
        ));
        let other = fixture([0.02, 0.03], [1.0e6, 2.0e6], [1.0e-4; 2], false, cx);
        assert_ne!(other.identity(), model.identity());
        assert!(
            other
                .advance_normal_at_load(&before, 0.005, load_controls(), cx)
                .is_err()
        );
        assert!(
            other
                .advance_wear_at_load(&before, &slip, 0.005, load_controls(), cx)
                .is_err()
        );

        let bounded = fixture([0.01, 0.03], [1.0e6, 2.0e6], [1.0e-4, 1.0e-12], false, cx);
        let bounded_before = bounded
            .advance_normal_at_load(&bounded.initial_state(), 0.005, load_controls(), cx)
            .unwrap()
            .candidate;
        let bounded_saved = bounded_before.clone();
        assert!(matches!(
            bounded.advance_wear_at_load(&bounded_before, &slip, 0.005, load_controls(), cx),
            Err(SourcePlasticWearError::Contact(
                SourcePlasticContactError::OutsideApplicability {
                    quantity: "body B cumulative recession [m]",
                    ..
                }
            ))
        ));
        assert_eq!(bounded_before, bounded_saved);
        gate.request();
        assert!(
            model
                .advance_normal_at_load(&before, 0.005, load_controls(), cx)
                .is_err()
        );
        assert!(
            model
                .advance_wear_at_load(&before, &slip, 0.005, load_controls(), cx)
                .is_err()
        );
        assert_eq!(before, saved);
    });
}

#[test]
fn local_wear_redistributes_two_cell_pressure_while_retaining_total_load() {
    with_cx(|cx, _| {
        let model = fixture([0.01, 0.03], [1.0e6, 2.0e6], [1.0e-4; 2], true, cx);
        let before = model
            .advance_normal_at_load(&model.initial_state(), 0.005, load_controls(), cx)
            .unwrap()
            .candidate;
        let saved = before.clone();
        let previous = before.contact_state().contact_state();
        close(previous.pressure_pa()[12], 0.005 / (2.0 * 1.0e-8));
        close(previous.pressure_pa()[13], previous.pressure_pa()[12]);
        let mut slip = [0.0; 25];
        slip[12] = 1.0e-4;
        let result = model
            .advance_wear_at_load(&before, &slip, 0.005, load_controls(), cx)
            .unwrap();
        let next = result.wear.candidate.contact_state().contact_state();
        close(result.load.achieved_force_n, 0.005);
        close(
            (next.pressure_pa()[12] + next.pressure_pa()[13]) * 1.0e-8,
            0.005,
        );
        assert!(next.pressure_pa()[12] < previous.pressure_pa()[12]);
        assert!(next.pressure_pa()[13] > previous.pressure_pa()[13]);
        assert!(next.pressure_pa()[12] > 0.0 && next.pressure_pa()[13] < 1.0e6);
        // Symmetry makes the approach increase half the recession of the one
        // worn cell at fixed total force, independently of the two compliances.
        close(next.recession_m()[12], 0.625e-6);
        assert_eq!(next.recession_m()[13], 0.0);
        close(next.approach_m() - previous.approach_m(), 0.3125e-6);
        close(
            result.wear.pressure_weighted_sliding_j,
            0.005 * 1.0e-4 / 2.0,
        );
        for side in 0..2 {
            let expected = [0.01, 0.03][side] * (0.005 * 1.0e-4 / 2.0) / [1.0e6, 2.0e6][side];
            close(result.wear.wear[side].increment_m3, expected);
            assert_eq!(result.wear.candidate.recession_m()[side][13], 0.0);
        }
        assert_eq!(result.wear.receipt.response.plastic_work_j, 0.0);
        assert_eq!(next.plastic_depth_m(), previous.plastic_depth_m());
        assert_eq!(before, saved);
    });
}
