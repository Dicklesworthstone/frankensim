//! G1/G3/G4 original PR endpoints, temporal convergence and bounded atomic work.
use super::*;
use network::adaptive::{
    DiffusionNetworkAdvanceControls, RealGasDiffusionNetworkAdvanceControls,
    RealGasDiffusionNetworkAdvanceError as Error, advance_candidate,
};

fn temporal() -> RealGasDiffusionNetworkAdvanceControls {
    RealGasDiffusionNetworkAdvanceControls {
        temporal: DiffusionNetworkAdvanceControls {
            maximum_step_s: 40.0,
            minimum_step_s: 1e-5,
            maximum_evaluations: 4096,
            absolute_species_tolerance_kg: 1e-7,
            relative_species_tolerance: 0.0,
            temperature_tolerance_k: 0.1,
        },
        maximum_face_evaluations: 8192,
    }
}

fn link<'a>(model: &'a RealGasDiffusionModel<'a>) -> RealGasDiffusionNetworkLink<'a> {
    RealGasDiffusionNetworkLink {
        id: 1,
        parcels: [10, 20],
        model,
        geometry: GasDiffusionGeometry {
            area_m2: 0.02,
            ..geometry()
        },
        solve: controls(),
    }
}

#[test]
fn accepted_fine_states_and_gross_enthalpy_are_the_original_two_transactions() {
    let card = card();
    let context = point(330.0, 4e6);
    for references in [[0.0; 2], [1e19, -2e19]] {
        let [left, right] = pair(references);
        let model = RealGasDiffusionModel::new(left.model(), &[binding(&card, &context)], || false)
            .unwrap();
        let links = [link(&model)];
        let original = BTreeMap::from([(10, left), (20, right)]);
        let endpoint_controls = endpoints(&original);
        let first = network::step_candidate(
            &original,
            &links,
            20.0,
            &endpoint_controls,
            physical(),
            || false,
        )
        .unwrap();
        let second = network::step_candidate(
            first.parcels(),
            &links,
            20.0,
            &endpoint_controls,
            physical(),
            || false,
        )
        .unwrap();
        let result = advance_candidate(
            &original,
            &links,
            40.0,
            &endpoint_controls,
            physical(),
            temporal(),
            || false,
        )
        .unwrap();
        assert_eq!(result.parcels(), second.parcels());
        let report = result.report();
        assert_eq!(report.duration_s, 40.0);
        assert_eq!(report.accepted_intervals, 1);
        assert_eq!(report.accepted_substeps, 2);
        assert_eq!(report.rejected_intervals, 0);
        assert_eq!(report.evaluations, 3);
        assert_eq!(report.face_evaluations, 3);
        assert_eq!(report.pair_evaluations, 0);
        assert_eq!(report.last_substep, *second.step());
        assert!(report.maximum_species_error_ratio <= 1.0);
        assert!(report.maximum_temperature_error_ratio <= 1.0);
        assert!(report.endpoint_energy_change_j.abs() <= physical().exchange.energy_tolerance_j);
        assert!(report.total_energy_residual_j.abs() <= physical().exchange.energy_tolerance_j);
        assert!(report.entropy_change_j_k > 0.0);
        assert!(report.total_mass_residual_kg.abs() <= physical().exchange.mass_tolerance_kg);
        assert!(report.maximum_species_mass_residual_kg <= physical().exchange.mass_tolerance_kg);
        for row in &report.transfers {
            let actual: Vec<_> = first
                .step()
                .exchange
                .fluxes
                .iter()
                .chain(&second.step().exchange.fluxes)
                .filter(|flux| {
                    flux.transfer.donor == row.donor
                        && flux.transfer.receiver == row.receiver
                        && flux.transfer.species == row.species
                })
                .collect();
            assert_eq!(actual.len(), 2);
            close(
                row.mass_kg,
                actual[0].transfer.mass_kg + actual[1].transfer.mass_kg,
                1e-19,
            );
            let expected = actual[0].carried_enthalpy_j + actual[1].carried_enthalpy_j;
            close(
                row.carried_enthalpy_j,
                expected,
                2e-15 * expected.abs().max(1.0),
            );
            close(
                species_mass(&result.parcels()[&row.receiver], &row.species)
                    - species_mass(&original[&row.receiver], &row.species),
                row.mass_kg,
                physical().exchange.mass_tolerance_kg,
            );
        }
        for (&id, old) in &original {
            let new = &result.parcels()[&id];
            assert_eq!(new.volume_m3(), old.volume_m3());
            assert_eq!(
                new.model().mechanics().binary_interactions(),
                old.model().mechanics().binary_interactions()
            );
        }
    }
}

#[test]
fn tighter_species_and_pr_temperature_control_converges_to_small_native_steps() {
    let [left, right] = pair([0.0; 2]);
    let card = card();
    let context = point(330.0, 4e6);
    let model =
        RealGasDiffusionModel::new(left.model(), &[binding(&card, &context)], || false).unwrap();
    let links = [link(&model)];
    let original = BTreeMap::from([(10, left), (20, right)]);
    let endpoint_controls = endpoints(&original);
    let coarse = advance_candidate(
        &original,
        &links,
        40.0,
        &endpoint_controls,
        physical(),
        temporal(),
        || false,
    )
    .unwrap();
    let mut tight = temporal();
    tight.temporal.absolute_species_tolerance_kg = 1e-11;
    tight.temporal.temperature_tolerance_k = 1e-6;
    let fine = advance_candidate(
        &original,
        &links,
        40.0,
        &endpoint_controls,
        physical(),
        tight,
        || false,
    )
    .unwrap();
    assert!(fine.report().accepted_substeps > coarse.report().accepted_substeps);
    assert!(fine.report().rejected_intervals > 0);
    assert_eq!(fine.report().face_evaluations, fine.report().evaluations);
    assert!(fine.report().evaluations > 3 * fine.report().accepted_intervals);
    assert!(fine.report().maximum_species_error_ratio <= 1.0);
    assert!(fine.report().maximum_temperature_error_ratio <= 1.0);
    let mut reference = original.clone();
    for _ in 0..256 {
        reference = network::step_candidate(
            &reference,
            &links,
            40.0 / 256.0,
            &endpoint_controls,
            physical(),
            || false,
        )
        .unwrap()
        .into_parts()
        .0;
    }
    let mut coarse_species_error = 0.0_f64;
    let mut fine_species_error = 0.0_f64;
    let mut coarse_temperature_error = 0.0_f64;
    let mut fine_temperature_error = 0.0_f64;
    for (&id, state) in &reference {
        for species in model.species() {
            let expected = species_mass(state, species);
            coarse_species_error = coarse_species_error
                .max((species_mass(&coarse.parcels()[&id], species) - expected).abs());
            fine_species_error = fine_species_error
                .max((species_mass(&fine.parcels()[&id], species) - expected).abs());
        }
        coarse_temperature_error = coarse_temperature_error.max(
            (coarse.parcels()[&id].state().temperature_k() - state.state().temperature_k()).abs(),
        );
        fine_temperature_error = fine_temperature_error.max(
            (fine.parcels()[&id].state().temperature_k() - state.state().temperature_k()).abs(),
        );
    }
    assert!(fine_species_error < 0.25 * coarse_species_error);
    assert!(fine_temperature_error < 0.25 * coarse_temperature_error);
}

#[test]
fn rejected_partial_face_probes_and_accepted_progress_share_finite_work_caps() {
    let [left, right] = pair([0.0; 2]);
    let card = card();
    let context = point(330.0, 4e6);
    let model =
        RealGasDiffusionModel::new(left.model(), &[binding(&card, &context)], || false).unwrap();
    let active = link(&model);
    let closed = RealGasDiffusionNetworkLink {
        id: 2,
        geometry: GasDiffusionGeometry {
            area_m2: 0.0,
            ..active.geometry
        },
        ..active
    };
    let links = [
        active,
        closed,
        RealGasDiffusionNetworkLink { id: 3, ..active },
    ];
    let original = BTreeMap::from([(10, left), (20, right)]);
    let saved = original.clone();
    let endpoint_controls = endpoints(&original);
    let mut guarded = physical();
    guarded.exchange.maximum_species_fraction = 1e-12;
    let mut bounded = temporal();
    bounded.maximum_face_evaluations = 3;
    assert_eq!(
        advance_candidate(
            &original,
            &links,
            40.0,
            &endpoint_controls,
            guarded,
            bounded,
            || false,
        ),
        Err(Error::FaceWorkLimit {
            elapsed_s: 0.0,
            accepted_intervals: 0,
            face_evaluations: 3,
            maximum_face_evaluations: 3,
        }),
    );
    // The first rejected probe consumed two calls. Its retry consumes only one
    // before the cap: neither a failed probe nor its partial face list is free.
    let mut bounded = temporal();
    bounded.temporal.maximum_step_s = 1.0;
    bounded.temporal.maximum_evaluations = 3;
    assert_eq!(
        advance_candidate(
            &original,
            &[active],
            2.0,
            &endpoint_controls,
            physical(),
            bounded,
            || false,
        ),
        Err(Error::WorkLimit {
            elapsed_s: 1.0,
            accepted_intervals: 1,
            evaluations: 3,
        }),
    );
    let mut bounded = temporal();
    bounded.maximum_face_evaluations = 1;
    let idle = advance_candidate(
        &original,
        &[closed],
        40.0,
        &endpoint_controls,
        physical(),
        bounded,
        || false,
    )
    .unwrap();
    assert_eq!(idle.parcels(), &original);
    assert_eq!(idle.report().face_evaluations, 0);
    assert_eq!(original, saved);
}

#[test]
fn source_incompatibility_is_fatal_and_unresolved_accuracy_cannot_publish() {
    let [left, right] = pair([0.0; 2]);
    let card = card();
    let context = point(330.0, 4e6);
    let model =
        RealGasDiffusionModel::new(left.model(), &[binding(&card, &context)], || false).unwrap();
    let links = [link(&model)];
    let incompatible = at(real_model(source(0.35, [0.0; 2]), 0.04), 330.0, 4e6, 22.0);
    let original = BTreeMap::from([(10, left.clone()), (20, incompatible)]);
    let mut bounded = temporal();
    bounded.maximum_face_evaluations = 1;
    assert!(matches!(
        advance_candidate(
            &original,
            &links,
            40.0,
            &endpoints(&original),
            physical(),
            bounded,
            || false,
        ),
        Err(Error::RealGasStep {
            elapsed_s: 0.0,
            step_s: 40.0,
            error: network::RealGasDiffusionNetworkError::Face { .. },
        })
    ));
    let original = BTreeMap::from([(10, left), (20, right)]);
    let mut exact = temporal();
    exact.temporal.minimum_step_s = 40.0;
    exact.temporal.absolute_species_tolerance_kg = 1e-30;
    exact.temporal.temperature_tolerance_k = 1e-12;
    assert!(matches!(
        advance_candidate(
            &original,
            &links,
            40.0,
            &endpoints(&original),
            physical(),
            exact,
            || false,
        ),
        Err(Error::Accuracy {
            elapsed_s: 0.0,
            step_s: 40.0,
            ..
        })
    ));
}

#[test]
fn late_cancellation_discards_the_complete_local_candidate() {
    let [left, right] = pair([0.0; 2]);
    let card = card();
    let context = point(330.0, 4e6);
    let model =
        RealGasDiffusionModel::new(left.model(), &[binding(&card, &context)], || false).unwrap();
    let links = [link(&model)];
    let original = BTreeMap::from([(10, left), (20, right)]);
    let saved = original.clone();
    let endpoint_controls = endpoints(&original);
    let mut polls = 0usize;
    let complete = advance_candidate(
        &original,
        &links,
        40.0,
        &endpoint_controls,
        physical(),
        temporal(),
        || {
            polls += 1;
            false
        },
    )
    .unwrap();
    assert!(complete.report().accepted_substeps > 0);
    let mut seen = 0usize;
    assert_eq!(
        advance_candidate(
            &original,
            &links,
            40.0,
            &endpoint_controls,
            physical(),
            temporal(),
            || {
                seen += 1;
                seen == polls
            },
        ),
        Err(Error::Cancelled),
    );
    assert_eq!(original, saved);
}
