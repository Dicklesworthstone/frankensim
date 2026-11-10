//! Whole heated-tube laminar wall properties in source and native heat paths.
use super::*;
use fs_airflow::graph::thermal::medium::convection::{
    evaluate_laminar_tube_from_material,
    implicit::PrescribedMassFlowConvection,
};

const LENGTH_OVER_DIAMETER: f64 = 20.0;

fn laminar_settings() -> LiquidExteriorConvectionSettings {
    // The native helper's AREA is the complete pi D L, not a wall patch.
    let diameter = (AREA / (core::f64::consts::PI * LENGTH_OVER_DIAMETER)).sqrt();
    let mut config = settings(SourceCellConvectionLaw::Forced {
        correlation: CorrelationId::SiederTateLaminar,
        speed_m_s: 0.001 / diameter, // Re=1000 at the authored mean bulk state
        settings: ForcedConvectionSettings {
            length_over_hydraulic_diameter: Some(LENGTH_OVER_DIAMETER),
            ..Default::default()
        },
    });
    config.characteristic_length = Length::new(diameter);
    config.property_temperature = ConvectionPropertyTemperature::Reservoir;
    config
}

// Independent Sieder-Tate formula and original synthetic source algebra.
// Bulk Re/Pr/k are held at the authored arithmetic mean Tbulk=300 K.
fn reference_coefficient(surface: f64, kinematic: bool) -> (f64, f64) {
    let offset = surface - 300.0;
    let (wall_mu, wall_log_slope) = if kinematic {
        let rho = 1000.0 - 0.3 * offset;
        let nu = 1e-6 - 1e-8 * offset;
        (rho * nu, -0.3 / rho - 1e-8 / nu)
    } else {
        let mu = 0.001 * (1.0 - 0.01 * offset);
        (mu, -1e-5 / mu)
    };
    let diameter = laminar_settings().characteristic_length.value();
    let coefficient = 1.86 * (8000.0 / LENGTH_OVER_DIAMETER).cbrt()
        * (0.001 / wall_mu).powf(0.14) * 0.5 / diameter;
    (coefficient, -0.14 * coefficient * wall_log_slope)
}

#[test]
fn g1_mean_laminar_liquid_retains_bulk_state_and_complete_wall_viscosity_tangent() {
    let body = curve();
    for kinematic in [false, true] {
        let knots = if kinematic {
            vec![(280.0, 1.2e-6), (320.0, 0.8e-6)]
        } else {
            vec![(280.0, 0.0012), (320.0, 0.0008)]
        };
        let card = viscosity_card(knots, kinematic);
        let point = QueryPoint::new()
            .with_quantity("T", temperature_quantity(), 300.0).unwrap()
            .with("p", 100000.0).unwrap()
            .with("solute_mass_fraction", 0.2).unwrap();
        let liquid = LiquidTemperatureLaw::try_new(
            &card, &point, "T",
            if kinematic { LiquidViscosityBasis::Kinematic } else { LiquidViscosityBasis::Dynamic },
            MaterialPropertySelection::SingleClaimOnly,
        ).unwrap();
        let config = laminar_settings();
        let law = ImplicitLiquidConvection::try_new(
            &body, [-20000.0, 20000.0], &liquid, config,
        ).unwrap();
        let bulk = liquid.at_temperature(Temperature::new(300.0)).unwrap();
        for surface in [290.0, 300.0, 310.0] {
            let actual = law.evaluate(surface, 300.0).unwrap();
            let (expected, derivative) = reference_coefficient(surface, kinematic);
            near(actual.coefficient.value.value(), expected, 2e-12 * expected);
            near(actual.coefficient_surface_derivative_w_m2_k2, derivative, 2e-11);
            assert_eq!(actual.property_temperature_k, 300.0);
            let original_bulk = bulk.transport_state().unwrap();
            for (prepared, original) in [
                (actual.transport.density_kg_m3(), original_bulk.density_kg_m3()),
                (actual.transport.specific_heat_j_kg_k(), original_bulk.specific_heat_j_kg_k()),
                (actual.transport.thermal_conductivity_w_m_k(), original_bulk.thermal_conductivity_w_m_k()),
                (actual.transport.dynamic_viscosity_pa_s(), original_bulk.dynamic_viscosity_pa_s()),
            ] {
                near(prepared, original, 2e-14 * original);
            }
            assert_eq!(actual.nusselt.card().id, CorrelationId::SiederTateLaminar);
            assert_eq!(actual.nusselt.groups()["L_over_Dh"], LENGTH_OVER_DIAMETER);
            assert!(actual.tube_section.is_none());
            assert!(actual.turbulent_property_correction.is_none());
            let wall = liquid.at_temperature(Temperature::new(surface)).unwrap();
            let wall_mu = wall.transport_state().unwrap().dynamic_viscosity_pa_s();
            near(actual.surface_dynamic_viscosity_pa_s.unwrap(), wall_mu, 2e-17);
            let direct = evaluate_laminar_tube_from_material(
                &bulk, &wall, 0.001 / config.characteristic_length.value(),
                config.characteristic_length, LENGTH_OVER_DIAMETER,
            ).unwrap();
            near(actual.coefficient.value.value(), direct.bulk().coefficient().value.value(), 2e-12 * expected);
            check_derivative(&law, surface, 300.0);
            let flux = law.flux_state_at_surface(0.0, surface, 300.0).unwrap();
            near(flux.convection_temperature_tangent_w_m2_k,
                expected + (surface - 300.0) * derivative, 2e-11 * expected);
        }
    }
}

#[test]
fn g3_mean_laminar_mass_flow_keeps_reynolds_and_surface_partials() {
    let body = curve();
    let card = liquid_card(0.3, 0.01, 320.0);
    let liquid = liquid_law(&card);
    let config = laminar_settings();
    let diameter = config.characteristic_length.value();
    let area = core::f64::consts::PI * diameter * diameter / 4.0;
    let rate = 1000.0 * (0.001 / diameter) * area;
    let source = ImplicitLiquidConvection::try_new(
        &body, [-20000.0, 20000.0], &liquid, config,
    ).unwrap();
    let flow = PrescribedMassFlowConvection::for_liquid(&source, rate, area).unwrap();
    let response = flow.evaluate(310.0, 300.0).unwrap();
    let (expected, derivative) = reference_coefficient(310.0, false);
    near(response.convection.coefficient.value.value(), expected, 2e-12 * expected);
    near(response.convection.coefficient_surface_derivative_w_m2_k2, derivative, 2e-11);
    let slope = response.coefficient_mass_flow_derivative_w_m2_k_s_kg.unwrap();
    near(slope, expected / (3.0 * rate), 2e-12 * slope);
    let delta = 1e-4 * rate;
    let at = |mdot| PrescribedMassFlowConvection::for_liquid(&source, mdot, area)
        .unwrap().evaluate(310.0, 300.0).unwrap().convection.coefficient.value.value();
    near(slope, (at(rate + delta) - at(rate - delta)) / (2.0 * delta), 1e-8 * slope);
}

#[test]
fn g3_source_callbacks_retain_actual_wall_source_and_whole_tube_heat() {
    let body = curve();
    let heat = initial(&body, [310.0, 290.0], 0.0, 0.0);
    let card = liquid_card(0.3, 0.01, 320.0);
    let liquid = liquid_law(&card);
    let config = laminar_settings();
    let bindings = IDS.map(|cell| SourceConvectiveCell {
        cell, reference_temperature: Temperature::new(300.0),
        characteristic_length: config.characteristic_length,
        property_temperature: config.property_temperature, law: config.law,
    });
    let mut queries = Vec::new();
    let result = evaluate_source_cell_convection(
        &heat, &bindings, 2,
        |binding, temperature| {
            queries.push((binding.cell, temperature.value()));
            liquid.at_temperature(temperature)
        },
        || false,
    ).unwrap();
    assert_eq!(queries.len(), 4);
    for (source, sample) in result.materials().iter().zip(result.convection().samples()) {
        let wall = heat.temperature_k(sample.cell).unwrap();
        let retained = source.surface.as_ref().unwrap();
        assert_eq!(source.property_temperature_k, 300.0);
        assert_eq!(retained.thermal.temperature_k(), wall);
        assert_eq!(retained.material.material().card_identity(), card.content_hash());
        assert_eq!(retained.material.material().query().axes()["p"], 100000.0);
        assert_eq!(retained.material.material().query().axes()["solute_mass_fraction"], 0.2);
        assert!(queries.contains(&(sample.cell, 300.0)));
        assert!(queries.contains(&(sample.cell, wall)));
        let expected = reference_coefficient(wall, false).0;
        let evaluation = sample.evaluation.as_ref().unwrap();
        assert_eq!(evaluation.nusselt.card().id, CorrelationId::SiederTateLaminar);
        assert!(evaluation.tube_section.is_none());
        near(evaluation.coefficient.value.value(), expected, 2e-12 * expected);
        near(sample.exposed_area_m2, AREA, 1e-15);
        near(sample.heat_into_cell_w, AREA * expected * (300.0 - wall), 2e-10);
    }
    // A remaining patch cannot inherit the whole-tube mean correlation.
    let partial = initial(&body, [310.0, 290.0], 0.0, 0.002);
    assert!(evaluate_source_cell_convection(
        &partial, &bindings, 2, |_, t| liquid.at_temperature(t), || false,
    ).is_err());
    assert_eq!(heat.time_s(), 0.0);
}

#[test]
fn g0_laminar_mean_requires_bulk_convention_whole_length_and_original_source_support() {
    let body = curve();
    let card = liquid_card(0.3, 0.01, 320.0);
    let liquid = liquid_law(&card);
    let mut config = laminar_settings();
    config.property_temperature = ConvectionPropertyTemperature::FilmMean;
    assert!(ImplicitLiquidConvection::try_new(
        &body, [-20000.0, 20000.0], &liquid, config,
    ).is_err());
    for authored in [
        ForcedConvectionSettings::default(),
        ForcedConvectionSettings {
            length_over_hydraulic_diameter: Some(LENGTH_OVER_DIAMETER),
            aspect_ratio: Some(0.5), ..Default::default()
        },
    ] {
        let mut config = laminar_settings();
        if let SourceCellConvectionLaw::Forced { ref mut settings, .. } = config.law {
            *settings = authored;
        }
        assert!(ImplicitLiquidConvection::try_new(
            &body, [-20000.0, 20000.0], &liquid, config,
        ).is_err());
    }
    for (speed_factor, length_ratio) in [(2.1, LENGTH_OVER_DIAMETER), (1.0, 10_000.0)] {
        let mut config = laminar_settings();
        if let SourceCellConvectionLaw::Forced { ref mut speed_m_s, ref mut settings, .. } = config.law {
            *speed_m_s *= speed_factor;
            settings.length_over_hydraulic_diameter = Some(length_ratio);
        }
        let source = ImplicitLiquidConvection::try_new(
            &body, [-20000.0, 20000.0], &liquid, config,
        ).unwrap();
        assert!(source.evaluate(310.0, 300.0).is_err());
    }
    let kinked = viscosity_card(
        vec![(280.0, 0.0012), (300.0, 0.001), (320.0, 0.0009)], false,
    );
    let liquid = liquid_law(&kinked);
    let source = ImplicitLiquidConvection::try_new(
        &body, [-20000.0, 20000.0], &liquid, laminar_settings(),
    ).unwrap();
    assert!(source.evaluate(300.0, 300.0).is_err());
    assert!(source.evaluate(310.0, 300.0).is_ok());
    assert!(source.evaluate(330.0, 300.0).is_err());
}

fn surface_reference(body_temperature: f64) -> (f64, [f64; 2]) {
    let powers = |surface: f64| [
        AREA * reference_coefficient(surface, false).0 * (surface - 300.0),
        AREA * 0.7 * 5.670_374_419e-8 * (surface.powi(4) - 305.0_f64.powi(4)),
    ];
    let (mut lo, mut hi) = (body_temperature.min(300.0), body_temperature.max(305.0));
    for _ in 0..90 {
        let surface = 0.5 * (lo + hi);
        let conducted = AREA / 0.001 * (0.5 * (body_temperature - surface)
            + 0.001 * ((body_temperature - 300.0).powi(2) - (surface - 300.0).powi(2)));
        let p = powers(surface);
        if conducted > p[0] + p[1] { lo = surface; } else { hi = surface; }
    }
    let surface = 0.5 * (lo + hi);
    (surface, powers(surface))
}

#[test]
fn g3_mean_laminar_wall_tangent_closes_native_surface_energy_and_entropy() {
    source_surface_energy_with(laminar_settings(), surface_reference);
}
