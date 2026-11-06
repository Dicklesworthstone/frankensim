//! G3/G4: adaptive finite activation uses original states and atomic publication.
use super::*;

fn adaptive_controls(tolerance: f64) -> NuclearAdaptiveControls {
    NuclearAdaptiveControls {
        step: StepConfig {
            max_dt_s: 0.5,
            energy_tolerance_j: 1e-6,
        },
        initial_dt_s: 0.5,
        minimum_dt_s: 1e-8,
        maximum_trials: 20_000,
        maximum_substeps: 20_000,
        tolerances: NuclearTimeTolerances {
            specific_energy_j_kg: tolerance,
            temperature_k: tolerance,
            nuclei: tolerance,
            nuclear_energy_j: tolerance,
            neutron_population: tolerance,
            neutron_count: tolerance,
            relative: 0.0,
        },
    }
}

fn heating() -> (
    DepositionSimulation<'static>,
    FusionNeutronActivationHeating,
) {
    let initial = material_cells();
    let state = activation_state(&[[100.0, 0.0, 0.0], [80.0, 0.0, 0.0]]);
    let source = history(&state);
    let cells = FusionNeutronActivationCellSource::bind_with_history(
        state,
        Time::new(0.0),
        bindings(),
        source,
    )
    .unwrap();
    let heating =
        FusionNeutronActivationHeating::new(&initial, cells, activation_budget()).unwrap();
    (initial, heating)
}

#[test]
fn adaptive_activation_resolves_depletion_and_material_feedback_on_original_history() {
    let (initial, heating) = heating();
    let (reference, reference_source, _) = heating
        .advance_to(&initial, 0.625, controls(1.0 / 1024.0), || Ok(()))
        .unwrap();
    let mut errors = Vec::new();
    for tolerance in [1e-2, 1e-5] {
        let (material, evolved, report) = heating
            .advance_adaptive_to(&initial, 0.625, adaptive_controls(tolerance), || Ok(()))
            .unwrap();
        let mut error = (material.energy_j() - reference.energy_j()).abs();
        for cell in [CellId(20), CellId(10)] {
            let source = evolved.source(cell).unwrap();
            let target = reference_source.source(cell).unwrap();
            let index = source.cell_index();
            for (&a, &b) in source.state().inventories()[index]
                .populations()
                .iter()
                .zip(target.state().inventories()[index].populations())
            {
                error = error.max((a - b).abs());
            }
            error = error.max(
                (source.state().neutrons().neutron_counts()[index]
                    - target.state().neutrons().neutron_counts()[index])
                    .abs(),
            );
            assert!(source.state().inventories()[index].populations()[0] < [100.0, 80.0][index]);
            assert_eq!(
                evolved
                    .last_step(cell)
                    .unwrap()
                    .activation()
                    .spatial
                    .neutrons
                    .injected_neutrons,
                0.0
            );
        }
        errors.push(error);
        assert!(report.maximum_accepted_normalized_error <= 1.0);
        assert_eq!(report.thermal.substeps, 2 * report.accepted_comparisons);
        assert_eq!(evolved.audit().substeps, report.thermal.substeps);
        if tolerance < 1e-4 {
            assert!(report.accuracy_rejections > 0);
        }
        close(
            material.energy_j() - initial.energy_j(),
            evolved.audit().thermal_credit_j,
            1e-7,
        );
        let source = evolved.source(CellId(20)).unwrap();
        let generated = 5.0 * evolved.audit().captures + 3.5 * evolved.audit().decays;
        close(
            source.state().kinetic_energy_j().unwrap()
                + evolved.cumulative_released().unwrap().total_j(),
            56.0 + generated,
            1e-7,
        );
    }
    assert!(errors[1] < errors[0] / 3.0, "refinement errors {errors:?}");
    assert!(errors[1] < 0.005, "refinement errors {errors:?}");
    assert_eq!(initial.time_s(), 0.0);
    assert_eq!(heating.time_s(), 0.0);
}

#[test]
fn adaptive_activation_refuses_exhausted_work_without_publishing_trials() {
    let (initial, heating) = heating();
    let mut controls = adaptive_controls(1e-8);
    controls.maximum_trials = 1;
    assert!(matches!(
        heating.advance_adaptive_to(&initial, 0.625, controls, || Ok(())),
        Err(NuclearAdaptiveError::TrialBudget)
    ));
    controls = adaptive_controls(1e-4);
    let mut successful_polls = 0;
    heating
        .advance_adaptive_to(&initial, 0.625, controls, || {
            successful_polls += 1;
            Ok(())
        })
        .unwrap();
    let mut polls = 0;
    let result = heating.advance_adaptive_to(&initial, 0.625, controls, || {
        polls += 1;
        if polls == successful_polls - 1 {
            Err(DepositionError::Cancelled)
        } else {
            Ok(())
        }
    });
    assert!(matches!(
        result,
        Err(NuclearAdaptiveError::Source(FissionThermalError::Thermal(
            DepositionError::Cancelled
        )))
    ));
    assert_eq!(heating.time_s(), 0.0);
    assert_eq!(initial.time_s(), 0.0);
    assert_eq!(
        heating.cumulative_released().unwrap(),
        EnergyChannels::default()
    );
    assert_eq!(
        heating.source(CellId(20)).unwrap().state().inventories()[0].populations(),
        &[100.0, 0.0, 0.0]
    );
}
