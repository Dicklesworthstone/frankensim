//! Original-source endpoint qualification after ternary stirring and diffusion.

use super::*;
use fs_flux::miscibility::multicomponent::advection::{
    MulticomponentMiscibilityAdvection, MulticomponentMiscibilityAdvectionControls,
    MulticomponentMiscibilityAdvectionFace,
};
use fs_flux::miscibility::multicomponent::material::ResolvedMulticomponentMiscibilityField;

fn stirred_numerical() -> MulticomponentMiscibilityField {
    let transport = MaxwellStefanV1::new(
        [2, 0, 1]
            .map(|i| (species(i), MolarMass::new(MASSES[i])))
            .into(),
    )
    .unwrap();
    MulticomponentMiscibilityNetwork::new(
        model([1, 2, 0]),
        transport,
        Temperature::new(300.0),
        AmountConcentration::new(1000.0),
        vec![MulticomponentMiscibilityCell { volume_m3: 0.001 }; 3],
        [[0, 1], [0, 2], [1, 2]]
            .map(|cells| MulticomponentMiscibilityFace {
                cells,
                area_m2: 1.0,
                distance_m: 1.0,
                diffusivity_matrix: 0,
            })
            .into(),
        vec![vec![0.0, 1e-6, 1e-4, 1e-6, 0.0, 1e-4, 1e-4, 1e-4, 0.0]],
        16 * 1024 * 1024,
        || false,
    )
    .unwrap()
    .bind(
        vec![
            vec![0.3, 0.3, 0.4],
            vec![0.31, 0.1, 0.59],
            vec![0.305, 0.2, 0.495],
        ],
        || false,
    )
    .unwrap()
}

fn circulation(field: &MulticomponentMiscibilityField, q: f64) -> MulticomponentMiscibilityAdvection {
    // Closed circulation 0 -> 1 -> 2 -> 0, in canonical graph face order.
    field
        .network()
        .prepare_advection(
            vec![
                MulticomponentMiscibilityAdvectionFace {
                    forward_m3_s: q,
                    reverse_m3_s: 0.0,
                },
                MulticomponentMiscibilityAdvectionFace {
                    forward_m3_s: 0.0,
                    reverse_m3_s: q,
                },
                MulticomponentMiscibilityAdvectionFace {
                    forward_m3_s: q,
                    reverse_m3_s: 0.0,
                },
            ],
            1e-13,
            || false,
        )
        .unwrap()
}

fn assert_source_endpoint(
    field: &ResolvedMulticomponentMiscibilityField<'_, '_>,
    source: &MulticomponentLiquidSolutionLaw<'_>,
    card: &MaterialCard,
    expected: &MulticomponentMiscibilityField,
) {
    assert!(std::ptr::eq(field.source(), source));
    assert_eq!(field.field(), expected);
    assert_eq!(field.maximum_resolved_cells(), 3);
    assert_eq!(field.maximum_relative_volume_mismatch(), 1e-12);
    for (cell, receipt) in field.cells().iter().enumerate() {
        let amounts = expected.cell_component_moles(cell).unwrap();
        let original = source
            .homogeneous(expected.network().temperature(), amounts, || false)
            .unwrap();
        // Equality checks all original density, viscosity, Cp and conductivity
        // receipts at the actual final composition, not the initial source point.
        assert_eq!(receipt, &original);
        assert_eq!(receipt.component_moles(), amounts);
        assert_eq!(receipt.thermodynamics().model(), source.thermodynamics());
        assert_eq!(receipt.liquid().material().card_identity(), card.content_hash());
        assert_eq!(receipt.liquid().material().state_coordinate("pressure"), Some(100_000.0));
        assert_eq!(receipt.liquid().material().state_coordinate("T"), Some(300.0));
        near(receipt.volume_m3(), 0.001);
        near(receipt.molar_concentration_mol_m3(), 1000.0);
    }
}

#[test]
fn stirred_paths_refresh_the_same_original_source_after_complete_transport() {
    let card = source_card(1.0, 0.0);
    let source = source(&card, model([2, 0, 1]), MASSES);
    let numerical = stirred_numerical();
    let flow = circulation(&numerical, 5e-4);
    let field = numerical
        .resolve_liquid_field(&source, 3, 1e-12, || false)
        .unwrap();
    let original_cells = field.cells().to_vec();
    let advection = MulticomponentMiscibilityAdvectionControls::default();
    let diffusion = MulticomponentMiscibilitySolve::default();

    let expected = numerical.advect_candidate(&flow, 0.1, advection, || false).unwrap();
    let actual = field.advect_candidate(&flow, 0.1, advection, || false).unwrap();
    assert_source_endpoint(actual.field(), &source, &card, expected.field());
    assert_eq!(actual.step(), expected.step());
    let (advected, report) = actual.into_parts();
    assert_eq!(advected.field().elapsed_s(), 0.1);
    assert_eq!(report.duration_s, 0.1);
    assert_ne!(advected.cells(), original_cells.as_slice());

    let expected = numerical.transport_candidate(&flow, 0.1, advection, diffusion, || false).unwrap();
    let actual = field.transport_candidate(&flow, 0.1, advection, diffusion, || false).unwrap();
    assert_source_endpoint(actual.field(), &source, &card, expected.field());
    assert_eq!(actual.step(), expected.step());
    let (transported, report) = actual.into_parts();
    assert_eq!(transported.field().elapsed_s(), 0.1);
    assert_eq!(report.duration_s, 0.1);

    let expected = numerical.advance_transport_candidate(&flow, 0.1, advection, diffusion, time(), || false).unwrap();
    let actual = field.advance_transport_candidate(&flow, 0.1, advection, diffusion, time(), || false).unwrap();
    assert_source_endpoint(actual.field(), &source, &card, expected.field());
    assert_eq!(actual.advance(), expected.advance());
    let (advanced, report) = actual.into_parts();
    assert_eq!(advanced.field().elapsed_s(), 0.1);
    assert!(report.accepted_intervals > 0);
    assert!(report.maximum_accepted_error_ratio <= 1.0);
    assert_eq!(field.field(), &numerical);
    assert_eq!(field.cells(), original_cells.as_slice());
}

#[test]
fn late_material_and_cancellation_refusals_preserve_stirring_state() {
    let numerical = stirred_numerical();
    let flow = circulation(&numerical, 1e-9);
    let advection = MulticomponentMiscibilityAdvectionControls::default();
    let diffusion = MulticomponentMiscibilitySolve::default();
    // Tiny but nonzero circulation leaves the established uphill ternary
    // cross-diffusion response visible at the source-domain boundary.
    let expected = numerical.transport_candidate(&flow, 0.1, advection, diffusion, || false).unwrap();
    assert!(expected.field().cell_mole_fractions(1).unwrap()[0] > 0.31);
    for (maximum_x_a, perturbation) in [(0.31, 0.0), (1.0, 1.0)] {
        let card = source_card(maximum_x_a, perturbation);
        let source = source(&card, model([0, 2, 1]), MASSES);
        let field = numerical.resolve_liquid_field(&source, 3, 1e-12, || false).unwrap();
        let original_cells = field.cells().to_vec();
        for error in [
            field.transport_candidate(&flow, 0.1, advection, diffusion, || false).unwrap_err(),
            field.advance_transport_candidate(&flow, 0.1, advection, diffusion, time(), || false).unwrap_err(),
        ] {
            if perturbation == 0.0 {
                assert!(matches!(error, Error::Source { cell: 1, .. }));
            } else {
                assert!(matches!(error, Error::VolumeMismatch { .. }));
            }
        }
        assert_eq!(field.field(), &numerical);
        assert_eq!(field.cells(), original_cells.as_slice());
    }

    let card = source_card(1.0, 0.0);
    let source = source(&card, model([1, 0, 2]), MASSES);
    let field = numerical.resolve_liquid_field(&source, 3, 1e-12, || false).unwrap();
    let original_cells = field.cells().to_vec();
    let mut total_polls = 0;
    let accepted = field.transport_candidate(&flow, 0.1, advection, diffusion, || {
        total_polls += 1;
        false
    }).unwrap();
    let mut polls = 0;
    assert!(matches!(field.transport_candidate(&flow, 0.1, advection, diffusion, || {
        polls += 1;
        polls == total_polls
    }), Err(Error::Cancelled)));
    assert_eq!(field.field(), &numerical);
    assert_eq!(field.cells(), original_cells.as_slice());
    let replay = field.transport_candidate(&flow, 0.1, advection, diffusion, || false).unwrap();
    assert_eq!(replay.field().field(), accepted.field().field());
    assert_eq!(replay.field().cells(), accepted.field().cells());
    assert_eq!(replay.step(), accepted.step());
}
