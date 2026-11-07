//! Coefficient activation must become native fracture, with no synthetic path.
use super::*;
use fs_solid::brittle::solid::expansion::initiation::ThermalTensileTetActivation;

fn activated<'a>(probe: ThermalTensileTetProbe<'a>) -> ThermalTensileTetActivation<'a> {
    match probe {
        ThermalTensileTetProbe::Activation(activation) => activation,
        _ => panic!("expected tensile activation from changed coefficients"),
    }
}

#[test]
fn coefficient_activation_replays_native_release_at_held_mechanical_strain_without_healing() {
    with_cx(false, |cx| {
        let nodes = points(0.03);
        let geometry =
            TetLoadGeometry3::new(&nodes, &TETS, TetAssemblyBudget::standard(), 4, cx).unwrap();
        let card = weakening_card();
        let reference = [resolved_expansion(&card, &geometry, 300.0, 300.0)];
        let hot = [reference[0].at_temperature(&card, 400.0).unwrap()];
        let virgin = ThermalTensileElasticTetSolid::try_new(
            &geometry,
            &reference,
            TetAssemblyBudget::standard(),
            cx,
        )
        .unwrap();
        let strain = 0.75 * STRENGTH / YOUNG;
        let initial_targets = affine_targets(&nodes, [0.0, strain, 0.0]);
        let accepted = elastic(
            virgin
                .probe(
                    &reference,
                    &initial_targets,
                    &[],
                    IncrementSettings::default(),
                    initiation_settings(),
                    controls(),
                    None,
                    cx,
                )
                .unwrap(),
        );
        let saved = accepted.solid.solution().clone();
        let free = hot[0].free_linear_strain();
        // Supports follow the exact integrated free strain, holding the physical
        // MECHANICAL preload. Coefficient work is the only source of damage.
        let targets = affine_targets(&nodes, [free, free + strain, free]);
        let prediction = activated(
            accepted
                .solid
                .probe(
                    &hot,
                    &targets,
                    &[],
                    IncrementSettings::default(),
                    initiation_settings(),
                    controls(),
                    None,
                    cx,
                )
                .unwrap(),
        );
        assert_eq!(prediction.initial_solution(), &saved);
        let observation = prediction.activations()[0].as_ref().unwrap();
        assert_eq!(observation.source().material(), hot[0].band().material());
        close(observation.native().strain()[1], strain, 1e-12);
        close(
            observation.native().stress_pa()[1],
            0.5 * YOUNG * strain,
            1e-12,
        );
        assert!(observation.native().normal()[1].abs() > 1.0 - 1e-14);
        assert!(
            observation.native().principal_stresses_pa()[2]
                > observation.native().tensile_strength_pa()
        );

        let bands = prediction.resolve_bands(cx).unwrap();
        assert_eq!(bands.initial[0].expansion(), reference[0].expansion());
        assert_eq!(bands.next[0].expansion(), hot[0].expansion());
        let handoff = prediction.handoff(&bands.initial, 1e-10, cx).unwrap();
        assert_eq!(
            handoff.solid.solution().displacement_m,
            saved.displacement_m
        );
        assert_eq!(handoff.solid.solution().element_states[0].damage(), 0.0);
        assert_eq!(handoff.solid.fracture_dissipation_j(), 0.0);
        close(
            handoff.solid.stored_energy_j(),
            accepted.solid.stored_energy_j(),
            1e-12,
        );
        assert!(handoff.storage_defect_j.abs() < 1e-14);

        let replay = prediction.replay(&bands, controls(), cx).unwrap();
        let volume = geometry.element_volumes_m3()[0];
        let width = 0.02 / 3.0;
        let new_strength = 0.25 * STRENGTH;
        let new_modulus = 0.5 * YOUNG;
        let onset_strain = new_strength / new_modulus;
        let failure_strain = 2.0 * FRACTURE_WORK / (new_strength * width);
        let stress = new_strength * (failure_strain - strain) / (failure_strain - onset_strain);
        let intact_storage = 0.5 * volume * new_modulus * strain * strain;
        let storage = 0.5 * volume * stress * strain;
        let released = intact_storage - storage;
        assert!(released > 0.0);
        close(replay.solid.solution().element_stress[0][1], stress, 2e-10);
        close(replay.solid.stored_energy_j(), storage, 2e-10);
        close(
            replay.energy.parameter_work_j,
            -0.5 * accepted.solid.stored_energy_j(),
            2e-10,
        );
        close(
            replay.energy.parameter_fracture_dissipation_j,
            released,
            2e-10,
        );
        close(replay.energy.fracture_dissipation_j, released, 2e-10);
        close(replay.solid.fracture_dissipation_j(), released, 2e-10);
        assert!(replay.energy.constitutive_strain_work_j.abs() < 1e-12);
        close(
            replay.energy.mechanical_work_j,
            volume * stress * free,
            2e-10,
        );
        close(
            replay.energy.eigenstrain_work_j,
            -volume * stress * free,
            2e-10,
        );
        check_balances(&replay.energy);
        let damage = replay.solid.solution().element_states[0].damage();
        close(damage, 1.0 - stress / (new_modulus * strain), 2e-10);

        // Cooling restores stiffness/strength, while the original damage and
        // cumulative fracture release remain physical history.
        let cooled = replay
            .solid
            .advance(
                &bands.initial,
                &initial_targets,
                &[],
                IncrementSettings::default(),
                controls(),
                cx,
            )
            .unwrap();
        assert!(cooled.solid.solution().element_states[0].damage() >= damage);
        close(
            cooled.solid.solution().element_states[0].damage(),
            damage,
            2e-10,
        );
        close(cooled.solid.fracture_dissipation_j(), released, 2e-10);
        assert!(cooled.energy.fracture_dissipation_j.abs() < 1e-12);
        check_balances(&cooled.energy);
        assert_eq!(accepted.solid.solution(), &saved);
    });
}

#[test]
fn coefficient_activation_admits_the_original_request_and_refuses_source_substitution_atomically() {
    with_cx(false, |cx| {
        let nodes = points(0.03);
        let geometry =
            TetLoadGeometry3::new(&nodes, &TETS, TetAssemblyBudget::standard(), 4, cx).unwrap();
        let card = weakening_card();
        let reference = [resolved_expansion(&card, &geometry, 300.0, 300.0)];
        let hot = [reference[0].at_temperature(&card, 400.0).unwrap()];
        let virgin = ThermalTensileElasticTetSolid::try_new(
            &geometry,
            &reference,
            TetAssemblyBudget::standard(),
            cx,
        )
        .unwrap();
        let targets = affine_targets(&nodes, [0.0, 0.75 * STRENGTH / YOUNG, 0.0]);
        let accepted = elastic(
            virgin
                .probe(
                    &reference,
                    &targets,
                    &[],
                    IncrementSettings::default(),
                    initiation_settings(),
                    controls(),
                    None,
                    cx,
                )
                .unwrap(),
        );
        let saved = accepted.solid.solution().clone();

        for forces in [vec![(12, 1.0)], vec![(7, f64::NAN)]] {
            assert!(
                accepted
                    .solid
                    .probe(
                        &hot,
                        &targets,
                        &forces,
                        IncrementSettings::default(),
                        initiation_settings(),
                        controls(),
                        None,
                        cx
                    )
                    .is_err()
            );
        }
        let mut bad_targets = targets.clone();
        bad_targets[0].1 = f64::INFINITY;
        assert!(
            accepted
                .solid
                .probe(
                    &hot,
                    &bad_targets,
                    &[],
                    IncrementSettings::default(),
                    initiation_settings(),
                    controls(),
                    None,
                    cx
                )
                .is_err()
        );
        let mut bad_settings = IncrementSettings::default();
        bad_settings.max_newton_iterations = 0;
        assert!(
            accepted
                .solid
                .probe(
                    &hot,
                    &targets,
                    &[],
                    bad_settings,
                    initiation_settings(),
                    controls(),
                    None,
                    cx
                )
                .is_err()
        );
        with_cx(true, |cancelled| {
            assert!(
                accepted
                    .solid
                    .probe(
                        &hot,
                        &targets,
                        &[],
                        IncrementSettings::default(),
                        initiation_settings(),
                        controls(),
                        None,
                        cancelled
                    )
                    .is_err()
            );
        });

        let prediction = activated(
            accepted
                .solid
                .probe(
                    &hot,
                    &targets,
                    &[],
                    IncrementSettings::default(),
                    initiation_settings(),
                    controls(),
                    None,
                    cx,
                )
                .unwrap(),
        );
        let bands = prediction.resolve_bands(cx).unwrap();
        // The same temperature and coefficients cannot substitute the old x
        // plane for the actual y plane after source-bound activation.
        let mut wrong_plane = bands.clone();
        wrong_plane.next[0] = hot[0].clone();
        assert!(prediction.replay(&wrong_plane, controls(), cx).is_err());
        assert!(prediction.handoff(&bands.next, 1e-10, cx).is_err());
        let mut wrong_reference = bands.clone();
        wrong_reference.next[0] = resolved_expansion(&card, &geometry, 400.0, 400.0)
            .with_resolved_band(bands.next[0].band().clone())
            .unwrap();
        assert_eq!(wrong_reference.next[0].band(), bands.next[0].band());
        assert!(prediction.replay(&wrong_reference, controls(), cx).is_err());
        assert_eq!(accepted.solid.solution(), &saved);
        let replay = prediction.replay(&bands, controls(), cx).unwrap();
        assert!(replay.energy.parameter_fracture_dissipation_j > 0.0);
        assert_eq!(accepted.solid.solution(), &saved);
    });
}
