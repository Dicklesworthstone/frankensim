//! G1/G3/G4: first-onset progress uses retained physical loads and native solves.
use super::*;
use fs_solid::brittle::solid::initiation::first_event::{
    TensileFirstEventControls, TensileFirstEventError,
};

fn event_controls() -> TensileFirstEventControls {
    TensileFirstEventControls {
        maximum_replays: 8,
        minimum_fraction: 1e-12,
    }
}

#[test]
fn native_first_event_progress_preserves_preload_work_and_exact_remainder() {
    with_cx(false, |cx| {
        let nodes = points(0.03);
        let geometry =
            TetLoadGeometry3::new(&nodes, &TETS, TetAssemblyBudget::standard(), 4, cx).unwrap();
        let source = card(STRENGTH, "first tensile event");
        let materials = [material(&source, &geometry, 0)];
        let virgin = TensileElasticTetSolid::try_new(
            &geometry,
            &materials,
            TetAssemblyBudget::standard(),
            cx,
        )
        .unwrap();
        let preload = 0.4 * STRENGTH / YOUNG;
        let target = prescribed_y(0.002);
        let accepted = elastic(
            virgin
                .probe(
                    &prescribed_y(preload),
                    &[],
                    IncrementSettings::default(),
                    initiation_settings(),
                    None,
                    cx,
                )
                .unwrap(),
        );
        let snapshot = accepted.solution().clone();
        let prediction = onset(
            accepted
                .probe(
                    &target,
                    &[],
                    IncrementSettings::default(),
                    initiation_settings(),
                    None,
                    cx,
                )
                .unwrap(),
        );
        let first = prediction
            .advance_before_first_event(event_controls(), cx)
            .unwrap();
        assert!(first.reached_fraction > 0.0);
        assert!(first.reached_fraction <= prediction.first_fraction_bracket()[0]);
        assert_eq!(first.remaining_fraction, 1.0 - first.reached_fraction);
        assert!(first.replays >= 1 && first.replays <= event_controls().maximum_replays);
        assert_eq!(first.prescribed_end_m, target);
        assert!(first.nodal_forces_end_n.is_empty());
        assert!(first.force_intervals_n.is_empty());
        let strain = (1.0 - first.reached_fraction) * preload + first.reached_fraction * 0.002;
        close(first.solid.solution().element_strain[0][1], strain, 2e-9);
        close(
            first.solid.solution().element_stress[0][1],
            YOUNG * strain,
            2e-9,
        );
        assert!(first.solid.solution().element_stress[0][1] < STRENGTH);
        let volume = geometry.element_volumes_m3()[0];
        close(
            first.elastic_work_j,
            0.5 * YOUNG * (strain * strain - preload * preload) * volume,
            2e-9,
        );
        assert_eq!(accepted.solution(), &snapshot);
        assert_eq!(prediction.initial_solution(), &snapshot);

        let remaining = onset(first.probe_remainder(cx).unwrap());
        assert_eq!(remaining.initial_solution(), first.solid.solution());
        let physical_fraction = (STRENGTH / YOUNG - preload) / (0.002 - preload);
        let b = remaining.first_fraction_bracket();
        let mapped = [
            first.reached_fraction + first.remaining_fraction * b[0],
            first.reached_fraction + first.remaining_fraction * b[1],
        ];
        assert!(mapped[0] <= physical_fraction + 2e-8);
        assert!(mapped[1] >= physical_fraction - 2e-8);
        let bands = remaining.resolve_bands(cx).unwrap();
        let handoff = remaining.handoff(&bands, 1e-9, cx).unwrap();
        close(
            handoff.solid.stored_energy_j(),
            first.solid.stored_energy_j(),
            2e-9,
        );
        let fracture = handoff
            .solid
            .advance(
                &first.prescribed_end_m,
                &first.nodal_forces_end_n,
                IncrementSettings::default(),
                controls(),
                cx,
            )
            .unwrap();
        assert!(fracture.solid.solution().element_states[0].damage() > 0.0);
        close(
            first.elastic_work_j + fracture.energy.material_work_j,
            fracture.solid.stored_energy_j() - accepted.stored_energy_j()
                + fracture.energy.fracture_dissipation_j,
            2e-9,
        );
    });
}

#[test]
fn original_applied_force_entries_and_support_loads_survive_event_splitting() {
    with_cx(false, |cx| {
        let nodes = points(0.03);
        let geometry =
            TetLoadGeometry3::new(&nodes, &TETS, TetAssemblyBudget::standard(), 4, cx).unwrap();
        let source = card(STRENGTH, "force-driven first event");
        let materials = [material(&source, &geometry, 0)];
        let virgin = TensileElasticTetSolid::try_new(
            &geometry,
            &materials,
            TetAssemblyBudget::standard(),
            cx,
        )
        .unwrap();
        let fixed: Vec<_> = prescribed(0.03, 0.0)
            .into_iter()
            .filter(|(d, _)| *d != 3)
            .collect();
        let peak_force = STRENGTH * geometry.element_volumes_m3()[0] / 0.03;
        let initial_forces = [(3, 0.2 * peak_force), (0, 17.0)];
        let target_forces = [(3, peak_force), (0, -31.0), (3, peak_force)];
        let accepted = elastic(
            virgin
                .probe(
                    &fixed,
                    &initial_forces,
                    IncrementSettings::default(),
                    initiation_settings(),
                    None,
                    cx,
                )
                .unwrap(),
        );
        assert_eq!(accepted.nodal_forces_n(), &initial_forces);
        let prediction = onset(
            accepted
                .probe(
                    &fixed,
                    &target_forces,
                    IncrementSettings::default(),
                    initiation_settings(),
                    None,
                    cx,
                )
                .unwrap(),
        );
        let first = prediction
            .advance_before_first_event(event_controls(), cx)
            .unwrap();
        assert_eq!(first.nodal_forces_end_n, target_forces);
        assert_eq!(first.prescribed_end_m, fixed);
        assert_eq!(first.force_intervals_n.len(), 2);
        let expected_axial =
            ((1.0 - first.reached_fraction) * 0.2 + first.reached_fraction * 2.0) * peak_force;
        let expected_support =
            (1.0 - first.reached_fraction) * 17.0 + first.reached_fraction * (-31.0);
        for (dof, values) in &first.force_intervals_n {
            let accepted_force = first
                .solid
                .nodal_forces_n()
                .iter()
                .find(|(d, _)| d == dof)
                .unwrap()
                .1;
            assert_eq!(values[0], accepted_force);
            match *dof {
                0 => {
                    close(values[0], expected_support, 2e-12);
                    assert_eq!(values[1], -31.0);
                }
                3 => {
                    close(values[0], expected_axial, 2e-12);
                    assert_eq!(values[1], 2.0 * peak_force);
                }
                _ => panic!("unexpected applied-force coordinate"),
            }
        }
        close(
            first.solid.solution().element_stress[0][0],
            expected_axial * 0.03 / geometry.element_volumes_m3()[0],
            2e-9,
        );
        let remaining = onset(first.probe_remainder(cx).unwrap());
        close(
            remaining.elastic_predictor().element_stress[0][0],
            2.0 * STRENGTH,
            2e-9,
        );
        assert_eq!(accepted.nodal_forces_n(), &initial_forces);
    });
}

#[test]
fn bounded_and_cancelled_event_requests_preserve_the_original_elastic_endpoint() {
    with_cx(false, |cx| {
        let nodes = points(0.03);
        let geometry =
            TetLoadGeometry3::new(&nodes, &TETS, TetAssemblyBudget::standard(), 4, cx).unwrap();
        let source = card(STRENGTH, "event refusal");
        let materials = [material(&source, &geometry, 0)];
        let virgin = TensileElasticTetSolid::try_new(
            &geometry,
            &materials,
            TetAssemblyBudget::standard(),
            cx,
        )
        .unwrap();
        let snapshot = virgin.solution().clone();
        let prediction = onset(
            virgin
                .probe(
                    &prescribed_y(0.002),
                    &[],
                    IncrementSettings::default(),
                    initiation_settings(),
                    None,
                    cx,
                )
                .unwrap(),
        );
        for maximum_replays in [0, 129] {
            assert!(matches!(
                prediction.advance_before_first_event(
                    TensileFirstEventControls {
                        maximum_replays,
                        ..event_controls()
                    },
                    cx,
                ),
                Err(TensileFirstEventError::Invalid(_))
            ));
        }
        assert!(matches!(
            prediction.advance_before_first_event(
                TensileFirstEventControls {
                    minimum_fraction: 0.5,
                    ..event_controls()
                },
                cx,
            ),
            Err(TensileFirstEventError::MinimumFraction { .. })
        ));
        with_cx(true, |cancelled| {
            assert!(
                prediction
                    .advance_before_first_event(event_controls(), cancelled)
                    .is_err()
            );
        });
        assert_eq!(virgin.solution(), &snapshot);
        let first = prediction
            .advance_before_first_event(event_controls(), cx)
            .unwrap();
        let repeated = prediction
            .advance_before_first_event(event_controls(), cx)
            .unwrap();
        assert_eq!(first.solid.solution(), repeated.solid.solution());
        assert_eq!(first.reached_fraction, repeated.reached_fraction);
        assert_eq!(first.force_intervals_n, repeated.force_intervals_n);

        // A deliberately unresolved numerical event can have a zero lower
        // bracket. No stationary state is mislabeled as positive progress.
        let coarse = TensileInitiationSettings {
            relative_stress_tolerance: 0.01,
            fraction_tolerance: 1.0,
            ..initiation_settings()
        };
        let near = elastic(
            virgin
                .probe(
                    &prescribed_y(0.999 * STRENGTH / YOUNG),
                    &[],
                    IncrementSettings::default(),
                    coarse,
                    None,
                    cx,
                )
                .unwrap(),
        );
        let unresolved = onset(
            near.probe(
                &prescribed_y(1.0001 * STRENGTH / YOUNG),
                &[],
                IncrementSettings::default(),
                coarse,
                None,
                cx,
            )
            .unwrap(),
        );
        assert_eq!(unresolved.first_fraction_bracket()[0], 0.0);
        assert!(matches!(
            unresolved.advance_before_first_event(event_controls(), cx),
            Err(TensileFirstEventError::NoProgress { fraction: 0.0 })
        ));
    });
}

#[test]
fn physical_hint_is_retained_when_the_original_target_is_reprobed() {
    with_cx(false, |cx| {
        let nodes = points(0.03);
        let geometry =
            TetLoadGeometry3::new(&nodes, &TETS, TetAssemblyBudget::standard(), 4, cx).unwrap();
        let source = card(STRENGTH, "first event hint");
        let materials = [material(&source, &geometry, 0)];
        let virgin = TensileElasticTetSolid::try_new(
            &geometry,
            &materials,
            TetAssemblyBudget::standard(),
            cx,
        )
        .unwrap();
        let strain = 2.0 * (1.0 - 2.0 * POISSON) * STRENGTH / YOUNG;
        let prescribed: Vec<_> = (0..12)
            .filter(|&dof| dof != 4)
            .map(|dof| {
                let length = match dof {
                    3 => 0.03,
                    7 => 0.02,
                    11 => 0.01,
                    _ => 0.0,
                };
                (dof, length * strain)
            })
            .collect();
        let hints = [Some([0.0, 1.0, 0.0])];
        let prediction = onset(
            virgin
                .probe(
                    &prescribed,
                    &[],
                    IncrementSettings::default(),
                    initiation_settings(),
                    Some(&hints),
                    cx,
                )
                .unwrap(),
        );
        let first = prediction
            .advance_before_first_event(event_controls(), cx)
            .unwrap();
        let remainder = onset(first.probe_remainder(cx).unwrap());
        let event = remainder.events()[0].as_ref().unwrap();
        assert_eq!(event.native().principal_multiplicity(), 3);
        close(event.native().normal()[1].abs(), 1.0, 2e-12);
        assert_eq!(first.prescribed_end_m, prescribed);
        assert_eq!(virgin.stored_energy_j(), 0.0);
    });
}
