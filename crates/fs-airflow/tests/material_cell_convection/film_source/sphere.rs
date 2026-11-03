//! G1/G4: original liquid bulk/surface queries on the finite-cell boundary.
use super::*;

fn sphere(cell: CellId) -> SourceConvectiveCell {
    SourceConvectiveCell {
        cell,
        reference_temperature: Temperature::new(300.0),
        characteristic_length: Length::new(0.1),
        property_temperature: ConvectionPropertyTemperature::Reservoir,
        law: SourceCellConvectionLaw::ForcedSphere { speed_m_s: 0.1 },
    }
}

#[test]
fn g1_sphere_boundary_uses_bulk_transport_and_original_wall_viscosity() {
    let chart = curve();
    let thermal = initial(&chart, [310.0, 300.0], 0.0, 0.002);
    let card = liquid_card(0.5, 0.01, 320.0);
    let law = liquid_law(&card);
    let mut queries = Vec::new();
    let result = evaluate_source_cell_convection(
        &thermal,
        &[sphere(IDS[0])],
        1,
        |binding, t| {
            queries.push((binding.cell, t.value()));
            law.at_temperature(t)
        },
        || false,
    )
    .unwrap();
    assert_eq!(queries, vec![(IDS[0], 300.0), (IDS[0], 310.0)]);
    let source = &result.materials()[0];
    let surface = source.surface.as_ref().unwrap();
    for material in [&source.material, &surface.material] {
        assert_eq!(material.material().card_identity(), card.content_hash());
        assert_eq!(material.material().query().axes()["p"], 100000.0);
        assert_eq!(
            material.material().query().axes()["solute_mass_fraction"],
            0.2
        );
    }
    assert_eq!(source.thermal.temperature_k(), 300.0);
    assert_eq!(surface.thermal.temperature_k(), 310.0);
    let sample = &result.convection().samples()[0];
    let evaluation = sample.evaluation.as_ref().unwrap();
    let ratio = 1.0_f64 / 0.9;
    let nu = 2.0
        + (0.4 * 10000.0_f64.sqrt() + 0.06 * 10000.0_f64.powf(2.0 / 3.0))
            * 8.0_f64.powf(0.4)
            * ratio.powf(0.25);
    near(evaluation.nusselt.groups()["Re"], 10000.0, 1e-9);
    near(evaluation.nusselt.groups()["Pr"], 8.0, 1e-12);
    near(
        evaluation.nusselt.groups()["mu_bulk_over_surface"],
        ratio,
        1e-12,
    );
    near(evaluation.coefficient.value.value(), nu * 0.5 / 0.1, 1e-9);
    near(
        sample.heat_into_cell_w,
        -10.0 * 0.008 * nu * 0.5 / 0.1,
        1e-9,
    );
    assert_eq!(thermal.temperature_k(IDS[0]).unwrap(), 310.0);
}

#[test]
fn g4_surface_source_and_sphere_geometry_refusals_do_not_publish_heat() {
    let chart = curve();
    let thermal = initial(&chart, [310.0, 290.0], 0.0, 0.002);
    let saved = thermal.cells().clone();
    let card = liquid_card(0.0, 0.01, 305.0);
    let law = liquid_law(&card);
    let result = evaluate_source_cell_convection(
        &thermal,
        &[sphere(IDS[0])],
        1,
        |_, t| law.at_temperature(t),
        || false,
    );
    assert!(matches!(
        result,
        Err(SourceCellConvectionError::Source {
            property_temperature_k,
            ..
        }) if property_temperature_k == 310.0
    ));
    let card = liquid_card(0.0, 0.01, 320.0);
    let law = liquid_law(&card);
    let mut calls = 0;
    let film = SourceConvectiveCell {
        property_temperature: ConvectionPropertyTemperature::FilmMean,
        ..sphere(IDS[0])
    };
    assert!(
        evaluate_source_cell_convection(
            &thermal,
            &[film],
            1,
            |_, t| {
                calls += 1;
                law.at_temperature(t)
            },
            || false
        )
        .is_err()
    );
    assert_eq!(calls, 0);
    for binding in [
        sphere(IDS[1]), // mu_bulk/mu_wall < 1 is outside the original card.
        SourceConvectiveCell {
            characteristic_length: Length::new(0.02),
            ..sphere(IDS[0])
        },
        SourceConvectiveCell {
            law: SourceCellConvectionLaw::ForcedSphere { speed_m_s: 0.0 },
            ..sphere(IDS[0])
        },
    ] {
        assert!(
            evaluate_source_cell_convection(
                &thermal,
                &[binding],
                1,
                |_, t| law.at_temperature(t),
                || false
            )
            .is_err()
        );
    }
    // A resolver cannot pass off the bulk snapshot as the requested wall state.
    assert!(
        evaluate_source_cell_convection(
            &thermal,
            &[sphere(IDS[0])],
            1,
            |_, _| law.at_temperature(Temperature::new(300.0)),
            || false
        )
        .is_err()
    );
    let queries = std::cell::Cell::new(0);
    assert!(matches!(
        evaluate_source_cell_convection(
            &thermal,
            &[sphere(IDS[0])],
            1,
            |_, t| {
                queries.set(queries.get() + 1);
                law.at_temperature(t)
            },
            || queries.get() == 1
        ),
        Err(SourceCellConvectionError::Convection(
            ConvectionSourceError::Cancelled
        ))
    ));
    assert_eq!(queries.get(), 1);
    assert_eq!(thermal.cells(), &saved);
}
