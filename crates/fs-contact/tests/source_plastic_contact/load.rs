//! Prescribed force keeps the original sources, physical bounds, and plastic history.

use super::*;
use fs_contact::normal_patch::FiniteGapPlasticLoadControls;
use fs_contact::normal_patch::plastic_source::load::SourcePlasticLoadError;

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
fn prescribed_load_retains_exact_cards_and_rejects_equal_valued_source_substitution() {
    with_cx(|cx, _| {
        let (normal, hardness) = sources([1.0e6, 2.0e6], "load-original", 0.0);
        let law = build(normal.clone(), hardness.clone(), envelope(), cx).unwrap();
        let (other_normal, other_hardness) =
            sources([1.0e6, 2.0e6], "different-original-test", 0.0);
        let other = build(other_normal, other_hardness.clone(), envelope(), cx).unwrap();
        assert_eq!(law.model().identity(), other.model().identity());
        assert_ne!(law.identity(), other.identity());
        assert_ne!(
            hardness[0].resolved().card_identity(),
            other_hardness[0].resolved().card_identity()
        );

        let before = law.initial_state();
        let saved = before.clone();
        let force = 0.005;
        let result = law
            .advance_at_load(&before, force, load_controls(), cx)
            .unwrap();
        assert_eq!(result.source_identity, law.identity());
        assert_eq!(result.load.model_identity, law.model().identity());
        assert_eq!(
            result.load.predecessor_state_identity,
            before.contact_state().identity()
        );
        assert_eq!(
            result.load.candidate_state_identity,
            result.candidate.contact_state().identity()
        );
        assert_eq!(result.load.response_identity, result.receipt.identity);
        assert_eq!(law.normal_source(), &normal);
        assert_eq!(law.hardness_sources(), &hardness);
        common::close(result.load.target_force_n, force);
        common::close(result.load.achieved_force_n, force);
        common::close(result.receipt.normal.normal_force_n, force);
        common::close(
            result.candidate.contact_state().pressure_pa()[12],
            force / 1.0e-8,
        );
        assert_eq!(result.receipt.plastic_work_j, 0.0);
        assert!(matches!(
            other.advance_at_load(&before, force, load_controls(), cx),
            Err(SourcePlasticLoadError::Contact(
                SourcePlasticContactError::Invalid(_)
            ))
        ));
        assert!(build(normal, other_hardness, envelope(), cx).is_err());
        assert_eq!(before, saved);
    });
}

#[test]
fn load_search_and_final_candidate_obey_original_approach_and_depth_bounds() {
    with_cx(|cx, _| {
        let (normal, hardness) = sources([1.0e6, 2.0e6], "load-bounds", 0.0);
        let mut approach_bound = envelope();
        approach_bound.maximum_approach_m = 1.0e-5;
        let law = build(normal.clone(), hardness.clone(), approach_bound, cx).unwrap();
        let before = law.initial_state();
        let saved = before.clone();
        assert!(matches!(
            law.advance_at_load(&before, 0.005, load_controls(), cx),
            Err(SourcePlasticLoadError::Contact(
                SourcePlasticContactError::OutsideApplicability {
                    quantity: "load-control maximum approach [m]",
                    ..
                }
            ))
        ));
        let bounded_search = FiniteGapPlasticLoadControls {
            maximum_approach_m: 1.0e-5,
            ..load_controls()
        };
        let accepted = law
            .advance_at_load(&before, 0.005, bounded_search, cx)
            .unwrap();
        common::close(accepted.receipt.normal.normal_force_n, 0.005);
        assert!(accepted.candidate.contact_state().approach_m() <= 1.0e-5);
        assert_eq!(before, saved);

        let mut shallow = envelope();
        shallow.half_space_depth_m = [1.0e-5; 2];
        let shallow_law = build(normal.clone(), hardness.clone(), shallow, cx).unwrap();
        let shallow_before = shallow_law.initial_state();
        let shallow_saved = shallow_before.clone();
        assert!(matches!(
            shallow_law.advance_at_load(&shallow_before, 0.005, load_controls(), cx),
            Err(SourcePlasticLoadError::Contact(
                SourcePlasticContactError::OutsideApplicability {
                    quantity: "patch radius / homogeneous depth",
                    ..
                }
            ))
        ));
        assert_eq!(shallow_before, shallow_saved);

        // The zero-gap cell yields before the recessed cell touches. At 0.015 N,
        // the first carries 0.01 N and the second carries 0.005 N: a unique
        // equilibrium whose permanent indentation exceeds the original bound.
        let mut gap = vec![1.0e-3; 25];
        gap[12] = 0.0;
        gap[13] = 2.0e-5;
        let mut depth_bound = envelope();
        depth_bound.maximum_plastic_depth_m = 1.0e-6;
        let depth_law = SourcePlasticContact::new(
            normal,
            hardness,
            FiniteGapGrid {
                cells_x: 5,
                cells_y: 5,
                cell_width_m: 1.0e-4,
                cell_depth_m: 1.0e-4,
            },
            gap,
            PlasticContactControls {
                maximum_active_set_iterations: 100,
                complementarity_tolerance_m: 1.0e-12,
                boundary_clearance_cells: 1,
            },
            depth_bound,
            cx,
        )
        .unwrap();
        let depth_before = depth_law.initial_state();
        let depth_saved = depth_before.clone();
        assert!(matches!(
            depth_law.advance_at_load(&depth_before, 0.015, load_controls(), cx),
            Err(SourcePlasticLoadError::Contact(
                SourcePlasticContactError::OutsideApplicability {
                    quantity: "permanent relative indentation [m]",
                    ..
                }
            ))
        ));
        assert_eq!(depth_before, depth_saved);
    });
}

#[test]
fn load_trials_cannot_heal_or_add_to_the_accepted_plastic_dent() {
    with_cx(|cx, gate| {
        let (normal, hardness) = sources([1.0e6, 2.0e6], "load-history", 0.0);
        let law = build(normal, hardness, envelope(), cx).unwrap();
        let loaded = law.advance(&law.initial_state(), 5.0e-5, cx).unwrap();
        let before = loaded.candidate;
        let saved = before.clone();
        let previous = before.contact_state();
        let dent = previous.plastic_depth_m()[12];
        assert!(dent > 0.0 && previous.cumulative_plastic_work_j() > 0.0);

        // The upper bracket trial yields, but the requested sub-cap load has
        // the unique elastic unloading solution on the original retained dent.
        let force = 0.005;
        let result = law
            .advance_at_load(&before, force, load_controls(), cx)
            .unwrap();
        let candidate = result.candidate.contact_state();
        common::close(result.receipt.normal.normal_force_n, force);
        common::close(
            candidate.approach_m(),
            dent + 0.5 * (previous.approach_m() - dent),
        );
        assert_eq!(candidate.plastic_depth_m(), previous.plastic_depth_m());
        assert_eq!(candidate.recession_m(), previous.recession_m());
        assert_eq!(
            candidate.cumulative_plastic_work_j(),
            previous.cumulative_plastic_work_j()
        );
        assert_eq!(result.receipt.plastic_work_j, 0.0);
        assert!(
            result
                .receipt
                .plastic_increment_m
                .iter()
                .all(|depth| *depth == 0.0)
        );

        let unloaded = law
            .advance_at_load(&result.candidate, 0.0, load_controls(), cx)
            .unwrap();
        assert_eq!(unloaded.receipt.normal.normal_force_n, 0.0);
        assert_eq!(
            unloaded.candidate.contact_state().approach_m(),
            load_controls().minimum_approach_m
        );
        assert_eq!(
            unloaded.candidate.contact_state().plastic_depth_m(),
            previous.plastic_depth_m()
        );
        assert_eq!(
            unloaded
                .candidate
                .contact_state()
                .cumulative_plastic_work_j(),
            previous.cumulative_plastic_work_j()
        );
        gate.request();
        assert!(
            law.advance_at_load(&before, force, load_controls(), cx)
                .is_err()
        );
        assert_eq!(before, saved);
    });
}
