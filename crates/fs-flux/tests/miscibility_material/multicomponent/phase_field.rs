//! Original-source ternary demixing, endpoint qualification and atomic refusals.
//! All thermodynamic coefficients and material surfaces are synthetic fixtures.

use super::*;
use fs_flux::miscibility::multicomponent::phase_field::{
    MulticomponentPhaseField, MulticomponentPhaseFieldNetwork,
    adaptive::MulticomponentPhaseFieldTimeControl,
    material::MulticomponentPhaseFieldMaterialError as PhaseMaterialError,
};
use fs_thermochem::solution::multicomponent::stability::MulticomponentStabilityControlV1;

const INITIAL: [[f64; 3]; 2] = [[0.31, 0.29, 0.4], [0.29, 0.31, 0.4]];

fn demixing_model(order: [usize; 3]) -> MulticomponentRegularSolutionV1 {
    MulticomponentRegularSolutionV1::new(
        order.map(species).into(),
        order
            .into_iter()
            .flat_map(|i| order.map(|j| if i == j { 0.0 } else { 4.0 }))
            .collect(),
        vec![0.0; 9],
        Temperature::new(290.0),
        Temperature::new(600.0),
    )
    .unwrap()
}

fn phase_field(initial: [[f64; 3]; 2]) -> MulticomponentPhaseField {
    let transport = MaxwellStefanV1::new(
        [1, 2, 0]
            .map(|i| (species(i), MolarMass::new(MASSES[i])))
            .into(),
    )
    .unwrap();
    MulticomponentPhaseFieldNetwork::new(
        demixing_model([2, 1, 0]),
        transport,
        Temperature::new(300.0),
        AmountConcentration::new(1000.0),
        vec![MulticomponentMiscibilityCell { volume_m3: 0.001 }; 2],
        vec![MulticomponentMiscibilityFace {
            cells: [0, 1],
            area_m2: 1.0,
            distance_m: 1.0,
            diffusivity_matrix: 0,
        }],
        vec![vec![0.0, 1e-6, 1e-4, 1e-6, 0.0, 1e-4, 1e-4, 1e-4, 0.0]],
        0.01,
        16 * 1024 * 1024,
        || false,
    )
    .unwrap()
    .bind(initial.map(Vec::from).into(), || false)
    .unwrap()
}

fn phase_time() -> MulticomponentPhaseFieldTimeControl {
    MulticomponentPhaseFieldTimeControl {
        composition: time(),
        absolute_interfacial_energy_tolerance_j: 1e-9,
        relative_interfacial_energy_tolerance: 0.0,
    }
}

#[test]
fn spinodal_demixing_refreshes_original_sources_without_flashing_local_cells() {
    let card = source_card(1.0, 0.0);
    // Canonicalization preserves full model equality across independent orders.
    let source = source(&card, demixing_model([1, 0, 2]), MASSES);
    let numerical = phase_field(INITIAL);
    let field = numerical
        .resolve_liquid_field(&source, 2, 1e-12, || false)
        .unwrap();
    let before_cells = field.cells().to_vec();
    let controls = MulticomponentMiscibilitySolve::default();
    let expected = numerical.step_candidate(0.1, controls, || false).unwrap();
    let candidate = field.step_candidate(0.1, controls, || false).unwrap();
    let end = candidate.field();
    assert_eq!(end.field(), expected.field());
    assert_eq!(candidate.step(), expected.step());
    assert_eq!(end.network(), numerical.network());
    assert!(std::ptr::eq(end.source(), &source));
    assert_eq!(end.maximum_resolved_cells(), 2);
    assert_eq!(end.maximum_relative_volume_mismatch(), 1e-12);
    assert_eq!(end.field().elapsed_s(), 0.1);
    // The A-rich cell enriches further; the positive interface cost opposes
    // demixing while the complete bulk-plus-interface free energy decreases.
    assert!(end.field().cell_mole_fractions(0).unwrap()[0] > INITIAL[0][0]);
    assert!(candidate.step().bulk_free_energy_change_j < 0.0);
    assert!(candidate.step().interfacial_energy_change_j > 0.0);
    assert!(candidate.step().free_energy_change_j <= controls.free_energy_tolerance_j);
    assert_ne!(end.cells(), before_cells.as_slice());
    for (cell, resolved) in end.cells().iter().enumerate() {
        let amounts = end.field().cell_component_moles(cell).unwrap();
        let fractions = end.field().cell_mole_fractions(cell).unwrap();
        let original = source
            .homogeneous(Temperature::new(300.0), amounts, || false)
            .unwrap();
        assert_eq!(resolved, &original);
        assert_eq!(resolved.component_moles(), amounts);
        assert_eq!(resolved.thermodynamics().model(), source.thermodynamics());
        assert!(
            resolved
                .thermodynamics()
                .local_stability(MulticomponentStabilityControlV1::default(), || false)
                .unwrap()
                .smallest_eigenvalue()
                .unwrap()
                < 0.0
        );
        near(resolved.volume_m3(), 0.001);
        near(resolved.molar_concentration_mol_m3(), 1000.0);
        let material = resolved.liquid().material();
        assert_eq!(material.card_identity(), card.content_hash());
        assert_eq!(material.state_coordinate("T"), Some(300.0));
        assert_eq!(material.state_coordinate("pressure"), Some(100_000.0));
        near(material.state_coordinate("x_a").unwrap(), fractions[0]);
        near(material.state_coordinate("x_b").unwrap(), fractions[1]);
        near(
            resolved.transport_state().unwrap().density_kg_m3(),
            44.0 - 42.0 * fractions[0] - 16.0 * fractions[1],
        );
        near(
            resolved.transport_state().unwrap().specific_heat_j_kg_k(),
            1500.0 + 1000.0 * fractions[0] + 2500.0 * fractions[1],
        );
        for property in material.properties() {
            let original = source
                .reference_state()
                .material()
                .property(property.requirement().name())
                .unwrap();
            assert_eq!(
                property.answer().receipt.selected,
                original.answer().receipt.selected
            );
        }
    }
    for species in 0..3 {
        near(
            end.field().cell_component_moles(0).unwrap()[species]
                + end.field().cell_component_moles(1).unwrap()[species],
            INITIAL[0][species] + INITIAL[1][species],
        );
    }
    let (accepted, step) = candidate.into_parts();
    assert_eq!(step.duration_s, 0.1);
    assert!(std::ptr::eq(accepted.source(), &source));
    let raw = accepted
        .field()
        .advance_candidate(0.1, phase_time(), controls, || false)
        .unwrap();
    let advanced = accepted
        .advance_candidate(0.1, phase_time(), controls, || false)
        .unwrap();
    assert_eq!(advanced.field().field(), raw.field());
    assert_eq!(advanced.advance(), raw.advance());
    let (advanced, advance) = advanced.into_parts();
    assert!(advance.accepted_intervals > 0);
    assert!(advance.maximum_accepted_error_ratio <= 1.0);
    assert!(advance.free_energy_change_j <= controls.free_energy_tolerance_j);
    assert_eq!(advanced.field().elapsed_s(), 0.2);
    assert!(std::ptr::eq(advanced.source(), &source));
    for (cell, resolved) in advanced.cells().iter().enumerate() {
        assert_eq!(
            resolved,
            &source
                .homogeneous(
                    Temperature::new(300.0),
                    advanced.field().cell_component_moles(cell).unwrap(),
                    || false,
                )
                .unwrap()
        );
    }
    assert_eq!(field.field(), &numerical);
    assert_eq!(field.cells(), before_cells.as_slice());
}

#[test]
fn phase_field_source_admission_reuses_model_mass_budget_and_volume_rules() {
    let card = source_card(1.0, 0.0);
    let source = source(&card, demixing_model([0, 1, 2]), MASSES);
    let numerical = phase_field(INITIAL);
    assert!(matches!(
        numerical.resolve_liquid_field(&source, 1, 1e-12, || false),
        Err(PhaseMaterialError::CellBudget {
            required: 2,
            maximum: 1,
        })
    ));
    for allowance in [-1.0, f64::NAN, 1e-5] {
        assert!(matches!(
            numerical.resolve_liquid_field(&source, 2, allowance, || false),
            Err(PhaseMaterialError::Invalid(_))
        ));
    }
    let other_model = self::source(&card, model([0, 1, 2]), MASSES);
    assert!(matches!(
        numerical.resolve_liquid_field(&other_model, 2, 1e-12, || false),
        Err(PhaseMaterialError::Invalid(_))
    ));
    let wrong_mass = self::source(
        &card,
        demixing_model([0, 1, 2]),
        [0.002, 0.028, 0.045],
    );
    assert!(matches!(
        numerical.resolve_liquid_field(&wrong_mass, 2, 1e-12, || false),
        Err(PhaseMaterialError::Invalid(_))
    ));
    assert!(matches!(
        numerical.resolve_liquid_field(&source, 2, 1e-12, || true),
        Err(PhaseMaterialError::Cancelled)
    ));
    let incompatible_card = source_card(1.0, 1.0);
    let incompatible_source =
        self::source(&incompatible_card, demixing_model([0, 1, 2]), MASSES);
    assert!(matches!(
        numerical.resolve_liquid_field(&incompatible_source, 2, 1e-12, || false),
        Err(PhaseMaterialError::VolumeMismatch { .. })
    ));
}

#[test]
fn phase_field_late_source_volume_and_cancellation_refusals_are_atomic() {
    let controls = MulticomponentMiscibilitySolve::default();
    for (initial, upper_x, perturbation) in [
        (INITIAL, 0.31, 0.0),
        ([[0.3, 0.3, 0.4], [0.31, 0.1, 0.59]], 1.0, 1.0),
    ] {
        let numerical = phase_field(initial);
        let card = source_card(upper_x, perturbation);
        let source = source(&card, demixing_model([2, 0, 1]), MASSES);
        let field = numerical
            .resolve_liquid_field(&source, 2, 1e-12, || false)
            .unwrap();
        let before_cells = field.cells().to_vec();
        let raw = numerical.step_candidate(0.1, controls, || false).unwrap();
        if perturbation == 0.0 {
            assert!(raw.field().cell_mole_fractions(0).unwrap()[0] > upper_x);
        }
        // Numerical demixing succeeds; only final original material support
        // or its fixed-volume compatibility refuses publication afterwards.
        numerical
            .advance_candidate(0.1, phase_time(), controls, || false)
            .unwrap();
        for error in [
            field.step_candidate(0.1, controls, || false).unwrap_err(),
            field
                .advance_candidate(0.1, phase_time(), controls, || false)
                .unwrap_err(),
        ] {
            if perturbation == 0.0 {
                assert!(matches!(error, PhaseMaterialError::Source { cell: 0, .. }));
            } else {
                assert!(matches!(error, PhaseMaterialError::VolumeMismatch { .. }));
            }
        }
        assert_eq!(field.field(), &numerical);
        assert_eq!(field.cells(), before_cells.as_slice());
        assert!(std::ptr::eq(field.source(), &source));
    }

    let numerical = phase_field(INITIAL);
    let card = source_card(1.0, 0.0);
    let source = source(&card, demixing_model([1, 2, 0]), MASSES);
    let field = numerical
        .resolve_liquid_field(&source, 2, 1e-12, || false)
        .unwrap();
    let before_cells = field.cells().to_vec();
    let mut completed_polls = 0;
    let accepted = field
        .step_candidate(0.1, controls, || {
            completed_polls += 1;
            false
        })
        .unwrap();
    // The last poll follows the complete numerical and source evaluations.
    let mut polls = 0;
    assert!(matches!(
        field.step_candidate(0.1, controls, || {
            polls += 1;
            polls == completed_polls
        }),
        Err(PhaseMaterialError::Cancelled)
    ));
    assert_eq!(field.field(), &numerical);
    assert_eq!(field.cells(), before_cells.as_slice());
    let replay = field.step_candidate(0.1, controls, || false).unwrap();
    assert_eq!(replay.field().field(), accepted.field().field());
    assert_eq!(replay.field().cells(), accepted.field().cells());
    assert_eq!(replay.step(), accepted.step());
}
