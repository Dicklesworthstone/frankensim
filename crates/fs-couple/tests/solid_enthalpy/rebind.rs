//! Source-card replacement reaches retained mechanics and reports both inventories.
use super::*;
use fs_solid::loads3::transient::events::MaterialRebindConfig3;

fn controls() -> MaterialRebindConfig3 {
    MaterialRebindConfig3 {
        maximum_absolute_work_j: 1.0,
        absolute_momentum_tolerance_kg_m_s: 1e-12,
        maximum_strain: 0.1,
    }
}

#[test]
fn g1_card_replacement_transfers_motion_then_continues_on_the_new_caloric_sources() {
    with_cx(|cx| {
        let fixed: Vec<_> = (0..12).filter(|d| *d != 3).collect();
        let zero = [[0.0; 3]; 4];
        let mut displacement = zero;
        let mut velocity = zero;
        displacement[1][0] = 0.001;
        velocity[1][0] = 0.02;
        let (old_card, old_caloric) =
            material_with_density(60.0, 60.0, 1e-5, 10.0, 0.0, false, 6.0);
        let (new_card, new_caloric) =
            material_with_density(180.0, 180.0, 2e-5, 20.0, 0.0, false, 12.0);
        let context = point(300.0);
        let old_thermal = initial(old_caloric.curve(), &[(30, 10_000.0)], 1.0);
        let new_thermal = initial(new_caloric.curve(), &[(30, 40_000.0)], 2.0);
        let old_model = SolidEnthalpy::try_new(
            &old_thermal,
            &NODES,
            &TETS,
            vec![source(30, &old_card, old_caloric.curve(), &context)],
            TetAssemblyBudget::standard(),
            0.1,
            config().mesh_frame,
            cx,
        )
        .unwrap();
        let new_model = SolidEnthalpy::try_new(
            &new_thermal,
            &NODES,
            &TETS,
            vec![source(30, &new_card, new_caloric.curve(), &context)],
            TetAssemblyBudget::standard(),
            0.1,
            config().mesh_frame,
            cx,
        )
        .unwrap();
        let loads = old_model
            .geometry()
            .assemble(&[], VolumeLoad3::None, &[], cx)
            .unwrap();
        let start = old_model.resolve(&old_thermal, cx).unwrap();
        let start_problem = start.problem(&fixed);
        let run = loads
            .begin_transient(
                &start_problem,
                Some(start.thermal_strains()),
                0.0,
                &displacement,
                &velocity,
                Default::default(),
                cx,
            )
            .unwrap();

        // Arrive at a nonzero accepted clock and step count before replacement.
        let old_hot = old_thermal
            .advance_to(0.01, thermal_config(), || Ok(()))
            .unwrap()
            .0;
        let new_hot = new_thermal
            .advance_to(0.01, thermal_config(), || Ok(()))
            .unwrap()
            .0;
        let before = old_model.resolve(&old_hot, cx).unwrap();
        let after = new_model.resolve(&new_hot, cx).unwrap();
        let before_problem = before.problem(&fixed);
        let after_problem = after.problem(&fixed);
        let (run, _) = before
            .advance_transient_candidate(&start, &run, &before_problem, &loads, cx)
            .unwrap();
        let q = run.state().displacement[0];
        let v = run.state().velocity[0];
        let (replacement, report) = after
            .rebind_transient_candidate(&before, &run, &after_problem, &loads, controls(), cx)
            .unwrap();
        assert_eq!(replacement.state().time_s, 0.01);
        assert_eq!(replacement.state().step, 1);
        assert_eq!(replacement.state().displacement, run.state().displacement);
        near(replacement.state().velocity[0], 0.5 * v, 1e-13);
        // Actual integrated alpha(T); source-derived new density doubles mass.
        let strain_before = 1e-5 * (10.0 + 0.005 * 100.0);
        let strain_after = 2.0 * strain_before;
        near(
            replacement.state().acceleration[0],
            -150.0 * (q - strain_after),
            1e-11,
        );
        near(report.mechanical.kinetic_work_j, -0.025 * v * v, 1e-14);
        let stored =
            |young: f64, strain: f64| young / 12.0 * ((q - strain).powi(2) + 2.0 * strain * strain);
        near(
            report.mechanical.elastic_work_j,
            stored(180.0, strain_after) - stored(60.0, strain_before),
            1e-14,
        );
        assert_eq!(report.cells[0].reference_mass_kg, [1.0, 2.0]);
        near(report.cells[0].caloric_chart_energy_j[0], 100.0, 1e-10);
        near(report.cells[0].caloric_chart_energy_j[1], 400.0, 1e-10);
        for (phase, card) in report.cells[0].phase.iter().zip([&old_card, &new_card]) {
            assert_eq!(phase.material_card_identity(), card.content_hash());
            near(phase.temperature_k(), 310.0, 1e-10);
        }
        assert_ne!(
            report.cells[0].phase[0].phase_curve_identity(),
            report.cells[0].phase[1].phase_curve_identity()
        );
        assert!(report.mechanical.reduced_momentum_residual_kg_m_s < 1e-12);

        // Ordinary heating now retains the replacement mass and source claims;
        // it does not repeat the external material-replacement event.
        let hotter = new_hot
            .advance_to(0.02, thermal_config(), || Ok(()))
            .unwrap()
            .0;
        let final_state = new_model.resolve(&hotter, cx).unwrap();
        let final_problem = final_state.problem(&fixed);
        let (continued, heating) = final_state
            .advance_transient_candidate(&after, &replacement, &final_problem, &loads, cx)
            .unwrap();
        assert_eq!(continued.state().step, 2);
        assert_eq!(continued.state().time_s, 0.02);
        assert!(heating.physical_balance_error_j.abs() < 1e-12);
        near(final_state.materials()[0].density_kg_m3, 12.0, 1e-12);
        near(
            continued
                .recover_response(config().mesh_frame, 0.1, cx)
                .unwrap()
                .stress
                .elements()[0]
                .stress_mandel_pa()[0],
            180.0 * (continued.state().displacement[0] - 2e-5 * (20.0 + 0.005 * 400.0)),
            1e-10,
        );
        assert_eq!(old_hot.cells()[&CellId(30)].mass_kg, 1.0);
        assert_eq!(new_hot.cells()[&CellId(30)].mass_kg, 2.0);
        assert_eq!(run.state().velocity[0], v);

        // The same external operation is reversible when no time has elapsed.
        let (restored, reverse) = before
            .rebind_transient_candidate(
                &after,
                &replacement,
                &before_problem,
                &loads,
                controls(),
                cx,
            )
            .unwrap();
        near(restored.state().velocity[0], v, 1e-13);
        near(
            report.mechanical.supplied_work_j + reverse.mechanical.supplied_work_j,
            0.0,
            1e-13,
        );
    });
}

#[test]
fn g4_source_replacement_refuses_unrelated_state_and_late_work_without_mutation() {
    with_cx(|cx| {
        let fixed: Vec<_> = (0..12).filter(|d| *d != 3).collect();
        let zero = [[0.0; 3]; 4];
        let mut displacement = zero;
        displacement[1][0] = 0.001;
        let context = point(300.0);
        let (old_card, old_caloric) =
            material_with_density(60.0, 60.0, 1e-5, 10.0, 0.0, false, 6.0);
        let (new_card, new_caloric) =
            material_with_density(180.0, 180.0, 2e-5, 20.0, 0.0, false, 12.0);
        let old_thermal = initial(old_caloric.curve(), &[(30, 0.0)], 1.0);
        let new_thermal = initial(new_caloric.curve(), &[(30, 0.0)], 2.0);
        let old_model = SolidEnthalpy::try_new(
            &old_thermal,
            &NODES,
            &TETS,
            vec![source(30, &old_card, old_caloric.curve(), &context)],
            TetAssemblyBudget::standard(),
            0.1,
            config().mesh_frame,
            cx,
        )
        .unwrap();
        let new_model = SolidEnthalpy::try_new(
            &new_thermal,
            &NODES,
            &TETS,
            vec![source(30, &new_card, new_caloric.curve(), &context)],
            TetAssemblyBudget::standard(),
            0.1,
            config().mesh_frame,
            cx,
        )
        .unwrap();
        let before = old_model.resolve(&old_thermal, cx).unwrap();
        let after = new_model.resolve(&new_thermal, cx).unwrap();
        let before_problem = before.problem(&fixed);
        let after_problem = after.problem(&fixed);
        let loads = old_model
            .geometry()
            .assemble(&[], VolumeLoad3::None, &[], cx)
            .unwrap();
        let run = loads
            .begin_transient(
                &before_problem,
                Some(before.thermal_strains()),
                0.0,
                &displacement,
                &zero,
                Default::default(),
                cx,
            )
            .unwrap();
        let state = run.state().clone();
        assert!(
            after
                .rebind_transient_candidate(&after, &run, &after_problem, &loads, controls(), cx)
                .is_err()
        );
        assert!(
            after
                .rebind_transient_candidate(&before, &run, &before_problem, &loads, controls(), cx)
                .is_err()
        );
        assert!(
            after
                .rebind_transient_candidate(
                    &before,
                    &run,
                    &after_problem,
                    &loads,
                    MaterialRebindConfig3 {
                        maximum_strain: 0.2,
                        ..controls()
                    },
                    cx
                )
                .is_err()
        );
        assert!(matches!(
            after.rebind_transient_candidate(
                &before,
                &run,
                &after_problem,
                &loads,
                MaterialRebindConfig3 {
                    maximum_absolute_work_j: 0.0,
                    ..controls()
                },
                cx
            ),
            Err(SolidEnthalpyError::Transient(
                fs_solid::loads3::transient::TransientError3::Invalid(_)
            ))
        ));

        let later_thermal = new_thermal
            .advance_to(0.01, thermal_config(), || Ok(()))
            .unwrap()
            .0;
        let later = new_model.resolve(&later_thermal, cx).unwrap();
        let later_problem = later.problem(&fixed);
        assert!(
            later
                .rebind_transient_candidate(&before, &run, &later_problem, &loads, controls(), cx)
                .is_err()
        );
        let unrelated_thermal = initial(new_caloric.curve(), &[(31, 0.0)], 2.0);
        let unrelated_model = SolidEnthalpy::try_new(
            &unrelated_thermal,
            &NODES,
            &TETS,
            vec![source(31, &new_card, new_caloric.curve(), &context)],
            TetAssemblyBudget::standard(),
            0.1,
            config().mesh_frame,
            cx,
        )
        .unwrap();
        let unrelated = unrelated_model.resolve(&unrelated_thermal, cx).unwrap();
        let unrelated_problem = unrelated.problem(&fixed);
        assert!(
            unrelated
                .rebind_transient_candidate(
                    &before,
                    &run,
                    &unrelated_problem,
                    &loads,
                    controls(),
                    cx
                )
                .is_err()
        );
        let gate = CancelGate::new();
        gate.request();
        ArenaPool::new(ArenaConfig::default()).scope(|arena| {
            let cancelled = Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 0,
                    kernel_id: 1,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            );
            assert!(matches!(
                after.rebind_transient_candidate(
                    &before,
                    &run,
                    &after_problem,
                    &loads,
                    controls(),
                    &cancelled
                ),
                Err(SolidEnthalpyError::Cancelled)
            ));
        });
        assert_eq!(run.state(), &state);
        assert_eq!(old_thermal.cells()[&CellId(30)].mass_kg, 1.0);
        assert_eq!(new_thermal.cells()[&CellId(30)].mass_kg, 2.0);
    });
}
