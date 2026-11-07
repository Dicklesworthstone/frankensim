//! Finite physical sections -> sourced tangents -> the native surface/history.
use super::*;
use fs_airflow::graph::thermal::medium::convection::{
    evaluate_flat_plate_section_from_material, evaluate_laminar_tube_section_from_material,
    evaluate_turbulent_tube_section_from_material,
    implicit::LiquidExteriorConvectionError,
};

fn plate(correlation: CorrelationId, speed_m_s: f64, start_m: f64) -> SourceCellConvectionLaw {
    SourceCellConvectionLaw::ForcedPlateSection {
        correlation, speed_m_s, start: Length::new(start_m),
    }
}

#[test]
fn g1_plate_section_film_chain_matches_original_sources_and_thin_strips() {
    let body = curve();
    let card = liquid_card(0.3, 0.01, 320.0);
    let liquid = liquid_law(&card);
    for (correlation, speed) in [
        (CorrelationId::FlatPlateLaminarAverage, 0.1),
        (CorrelationId::FlatPlateTurbulentAverage, 40.0),
    ] {
        for start in [0.0, 0.025, 0.1 * (1.0 - 1e-12)] {
            let convection = ImplicitLiquidConvection::try_new(
                &body, [-20000.0, 20000.0], &liquid, settings(plate(correlation, speed, start)),
            ).unwrap();
            assert_eq!(convection.transport_sources().count(), 4);
            for surface in [290.0, 310.0] {
                let actual = convection.evaluate(surface, 300.0).unwrap();
                let film = 0.5 * (surface + 300.0);
                let original = liquid.at_temperature(Temperature::new(film)).unwrap();
                let expected = evaluate_flat_plate_section_from_material(
                    correlation, &original, speed, Length::new(start), Length::new(0.1),
                ).unwrap();
                let h = expected.coefficient().value.value();
                near(actual.coefficient.value.value(), h, 1e-11 * h);
                assert_eq!(actual.property_temperature_k, film);
                assert_eq!(actual.nusselt.groups(), expected.section().downstream().groups());
                assert!(actual.tube_section.is_none());
                let details = actual.plate_section.unwrap();
                assert_eq!(details.start.value(), start);
                assert_eq!(details.end.value(), 0.1);
                near(details.nusselt, expected.section().nusselt(), 1e-10);
                near(
                    details.nusselt * actual.transport.thermal_conductivity_w_m_k() / 0.1,
                    actual.coefficient.value.value(), 1e-10 * h,
                );
                assert_eq!(details.upstream.is_some(), start > 0.0);
                if start > 0.0 {
                    assert!(actual.coefficient.model.discrepancy_rel.is_infinite());
                    assert_ne!(details.nusselt, actual.nusselt.evidence().value);
                } else {
                    assert_eq!(details.nusselt, actual.nusselt.evidence().value);
                }
                check_derivative(&convection, surface, 300.0);
            }
        }
    }
}

#[test]
fn g1_kinematic_plate_section_keeps_density_and_viscosity_temperature_chain() {
    let card = viscosity_card(vec![(280.0, 1.2e-6), (320.0, 0.8e-6)], true);
    let point = QueryPoint::new()
        .with_quantity("T", temperature_quantity(), 300.0).unwrap()
        .with("p", 100000.0).unwrap()
        .with("solute_mass_fraction", 0.2).unwrap();
    let liquid = LiquidTemperatureLaw::try_new(
        &card, &point, "T", LiquidViscosityBasis::Kinematic,
        MaterialPropertySelection::SingleClaimOnly,
    ).unwrap();
    let body = curve();
    let convection = ImplicitLiquidConvection::try_new(
        &body, [-20000.0, 20000.0], &liquid,
        settings(plate(CorrelationId::FlatPlateLaminarAverage, 0.1, 0.025)),
    ).unwrap();
    for surface in [290.0, 310.0] {
        check_derivative(&convection, surface, 300.0);
        let original = liquid.at_temperature(Temperature::new(0.5 * (surface + 300.0)))
            .unwrap();
        let expected = evaluate_flat_plate_section_from_material(
            CorrelationId::FlatPlateLaminarAverage, &original, 0.1,
            Length::new(0.025), Length::new(0.1),
        ).unwrap();
        near(convection.evaluate(surface, 300.0).unwrap().coefficient.value.value(),
            expected.coefficient().value.value(), 1e-9);
    }
}

#[test]
fn g0_plate_sections_refuse_source_kinks_wrong_convention_and_moving_regime_bounds() {
    let body = curve();
    let card = liquid_card(0.3, 0.01, 320.0);
    let liquid = liquid_law(&card);
    let mut config = settings(plate(CorrelationId::FlatPlateLaminarAverage, 0.1, 0.025));
    config.property_temperature = ConvectionPropertyTemperature::Reservoir;
    assert!(ImplicitLiquidConvection::try_new(
        &body, [-20000.0, 20000.0], &liquid, config,
    ).is_err());
    config.property_temperature = ConvectionPropertyTemperature::FilmMean;
    let law = ImplicitLiquidConvection::try_new(
        &body, [-20000.0, 20000.0], &liquid, config,
    ).unwrap();
    assert!(law.evaluate(319.0, 330.0).is_err());

    for (speed, start) in [(10.0, 0.025), (10.0, 0.05)] {
        let law = ImplicitLiquidConvection::try_new(
            &body, [-20000.0, 20000.0], &liquid,
            settings(plate(CorrelationId::FlatPlateTurbulentAverage, speed, start)),
        ).unwrap();
        // At film 300 K the downstream Re=1e6 is admitted. Its positive
        // upstream endpoint is either below or exactly at the turbulent bound.
        let error = law.evaluate(310.0, 290.0).unwrap_err();
        if start == 0.05 {
            assert!(matches!(error, LiquidExteriorConvectionError::MissingDerivative { axis: "Re" }));
        }
    }
    let kinked_card = viscosity_card(
        vec![(280.0, 0.0012), (300.0, 0.001), (320.0, 0.0009)], false,
    );
    let kinked_liquid = liquid_law(&kinked_card);
    let law = ImplicitLiquidConvection::try_new(
        &body, [-20000.0, 20000.0], &kinked_liquid, config,
    ).unwrap();
    assert!(law.evaluate(310.0, 290.0).is_err());
}

#[test]
fn g1_fixed_reynolds_bound_keeps_independent_prandtl_and_conductivity_chain() {
    let body = curve();
    // rho and mu are constant; Cp and k vary. Re stays exactly at the laminar
    // upper boundary, so only Pr and k need local source partials.
    let card = liquid_card(0.0, 0.0, 320.0);
    let liquid = liquid_law(&card);
    let convection = ImplicitLiquidConvection::try_new(
        &body, [-20000.0, 20000.0], &liquid,
        settings(plate(CorrelationId::FlatPlateLaminarAverage, 5.0, 0.025)),
    ).unwrap();
    let actual = convection.evaluate(310.0, 300.0).unwrap();
    assert_eq!(actual.nusselt.groups()["Re"], 500000.0);
    assert!(!actual.coefficient.sensitivity.d_qoi.contains_key("Re"));
    assert!(actual.coefficient.sensitivity.d_qoi.contains_key("Pr"));
    check_derivative(&convection, 310.0, 300.0);
}

#[test]
fn g1_tube_sections_retain_entrance_geometry_bulk_source_and_zero_wall_tangent() {
    let body = curve();
    let card = liquid_card(0.3, 0.01, 320.0);
    let liquid = liquid_law(&card);
    for start in [0.0, 1.0, 2.0 * (1.0 - 1e-12)] {
        let mut config = settings(SourceCellConvectionLaw::ForcedTubeSection {
            speed_m_s: 0.1, start: Length::new(start), end: Length::new(2.0),
        });
        assert!(ImplicitLiquidConvection::try_new(
            &body, [-20000.0, 20000.0], &liquid, config,
        ).is_err());
        config.property_temperature = ConvectionPropertyTemperature::Reservoir;
        let convection = ImplicitLiquidConvection::try_new(
            &body, [-20000.0, 20000.0], &liquid, config,
        ).unwrap();
        for surface in [290.0, 300.0, 310.0] {
            let actual = convection.evaluate(surface, 300.0).unwrap();
            let original = liquid.at_temperature(Temperature::new(300.0)).unwrap();
            let expected = evaluate_turbulent_tube_section_from_material(
                &original, 0.1, Length::new(0.1), Length::new(start), Length::new(2.0),
            ).unwrap();
            near(actual.coefficient.value.value(), expected.coefficient().value.value(), 1e-9);
            assert_eq!(actual.property_temperature_k, 300.0);
            assert_eq!(actual.coefficient_surface_derivative_w_m2_k2, 0.0);
            assert_eq!(actual.nusselt.groups(), expected.section().downstream().groups());
            assert!(actual.plate_section.is_none());
            let details = actual.tube_section.unwrap();
            assert_eq!(details.diameter.value(), 0.1);
            assert_eq!(details.start.value(), start);
            assert_eq!(details.end.value(), 2.0);
            near(details.nusselt, expected.section().nusselt(), 1e-12);
            assert_eq!(details.upstream.is_some(), start > 0.0);
            if start > 0.0 {
                assert!(actual.coefficient.model.discrepancy_rel.is_infinite());
                assert_ne!(details.nusselt, actual.nusselt.evidence().value);
            }
            check_derivative(&convection, surface, 300.0);
        }
    }
}

// Independent scalar section formula, source algebra and surface bisection;
// neither the production section evaluator nor its tangent enters this oracle.
fn section_reference(body_temperature: f64, law: SourceCellConvectionLaw) -> (f64, [f64; 2]) {
    let powers = |surface: f64| {
        let offset = if matches!(law, SourceCellConvectionLaw::ForcedTubeSection { .. }
            | SourceCellConvectionLaw::ForcedLaminarTubeSection { .. }) {
            0.0
        } else { 0.5 * (surface + 300.0) - 300.0 };
        let (rho, cp, k, mu) = (
            1000.0 - 0.3 * offset, 4000.0 + 2.0 * offset,
            0.5 + 0.002 * offset, 0.001 * (1.0 - 0.01 * offset),
        );
        let pr: f64 = mu * cp / k;
        let nusselt = match law {
            SourceCellConvectionLaw::ForcedPlateSection { correlation, speed_m_s, start } => {
                let re = rho * speed_m_s * 0.1 / mu;
                let fraction = start.value() / 0.1;
                match correlation {
                    CorrelationId::FlatPlateLaminarAverage =>
                        0.664 * re.sqrt() * pr.cbrt() / (1.0 + fraction.sqrt()),
                    CorrelationId::FlatPlateTurbulentAverage => {
                        let leading = 0.037 * re.powf(0.8);
                        if fraction == 0.0 { (leading - 871.0) * pr.cbrt() }
                        else { leading * (1.0 - fraction.powf(0.8)) / (1.0 - fraction) * pr.cbrt() }
                    }
                    _ => unreachable!("plate oracle uses original plate cards only"),
                }
            }
            SourceCellConvectionLaw::ForcedTubeSection { speed_m_s, start, end } => {
                let re = rho * speed_m_s * 0.1 / mu;
                let friction = (1.8 * re.log10() - 1.5).powi(-2);
                let developed = friction / 8.0 * (re - 1000.0) * pr
                    / (1.0 + 12.7 * (friction / 8.0).sqrt() * (pr.powf(2.0 / 3.0) - 1.0));
                let (a, b) = ((start.value() / 0.1).cbrt(), (end.value() / 0.1).cbrt());
                developed * (1.0 + 1.0 / (a * a + a * b + b * b))
            }
            SourceCellConvectionLaw::ForcedLaminarTubeSection { speed_m_s, start, end } => {
                let peclet = rho * speed_m_s * 0.1 * cp / k;
                let cumulative = |x: f64| {
                    if x == 0.0 { return 0.0; }
                    let gz = peclet / (x / 0.1);
                    x * (3.66 + 0.0668 * gz / (1.0 + 0.04 * gz.powf(2.0 / 3.0)))
                };
                (cumulative(end.value()) - cumulative(start.value())) / (end.value() - start.value())
            }
            _ => unreachable!("finite-section oracle only"),
        };
        [
            AREA * (nusselt * k / 0.1) * (surface - 300.0),
            AREA * 0.7 * 5.670_374_419e-8 * (surface.powi(4) - 305.0_f64.powi(4)),
        ]
    };
    let (mut lo, mut hi) = (body_temperature.min(300.0), body_temperature.max(305.0));
    for _ in 0..90 {
        let surface = 0.5 * (lo + hi);
        let conductive = AREA / 0.001 * (0.5 * (body_temperature - surface)
            + 0.001 * ((body_temperature - 300.0).powi(2) - (surface - 300.0).powi(2)));
        let p = powers(surface);
        if conductive > p[0] + p[1] { lo = surface; } else { hi = surface; }
    }
    let surface = 0.5 * (lo + hi);
    (surface, powers(surface))
}

#[test]
fn g3_plate_sections_close_native_surface_and_backward_euler_heat_and_entropy() {
    for (correlation, speed) in [
        (CorrelationId::FlatPlateLaminarAverage, 0.1),
        (CorrelationId::FlatPlateTurbulentAverage, 40.0),
    ] {
        let law = plate(correlation, speed, 0.025);
        source_surface_energy_with(settings(law), |t| section_reference(t, law));
    }
}

#[test]
fn g3_tube_section_closes_native_surface_and_backward_euler_heat_and_entropy() {
    let law = SourceCellConvectionLaw::ForcedTubeSection {
        speed_m_s: 0.1, start: Length::new(1.0), end: Length::new(2.0),
    };
    let mut config = settings(law);
    config.property_temperature = ConvectionPropertyTemperature::Reservoir;
    source_surface_energy_with(config, |t| section_reference(t, law));
}

#[test]
fn g1_laminar_sections_retain_original_bulk_source_and_exact_heat_tangent() {
    let body = curve();
    let card = liquid_card(0.3, 0.01, 320.0);
    let liquid = liquid_law(&card);
    for start in [0.0, 1.0, 2.0 * (1.0 - 1e-12)] {
        let mut config = settings(SourceCellConvectionLaw::ForcedLaminarTubeSection {
            speed_m_s: 0.01, start: Length::new(start), end: Length::new(2.0),
        });
        assert!(ImplicitLiquidConvection::try_new(
            &body, [-20000.0, 20000.0], &liquid, config,
        ).is_err(), "Hausen bulk properties must not silently become film properties");
        config.property_temperature = ConvectionPropertyTemperature::Reservoir;
        let convection = ImplicitLiquidConvection::try_new(
            &body, [-20000.0, 20000.0], &liquid, config,
        ).unwrap();
        for surface in [290.0, 300.0, 310.0] {
            let actual = convection.evaluate(surface, 300.0).unwrap();
            let original = liquid.at_temperature(Temperature::new(300.0)).unwrap();
            let expected = evaluate_laminar_tube_section_from_material(
                &original, 0.01, Length::new(0.1), Length::new(start), Length::new(2.0),
            ).unwrap();
            near(actual.coefficient.value.value(), expected.coefficient().value.value(), 1e-9);
            assert_eq!(actual.property_temperature_k, 300.0);
            assert_eq!(actual.coefficient_surface_derivative_w_m2_k2, 0.0);
            assert_eq!(actual.nusselt.groups(), expected.section().downstream().groups());
            let details = actual.tube_section.as_ref().unwrap();
            assert_eq!(details.start.value(), start);
            assert_eq!(details.end.value(), 2.0);
            assert_eq!(details.upstream.is_some(), start > 0.0);
            assert_eq!(actual.coefficient.model.discrepancy_rel.is_infinite(), start > 0.0);
            near(details.nusselt, expected.section().nusselt(), 1e-12);
            check_derivative(&convection, surface, 300.0);
            let flux = convection.flux_state_at_surface(0.0, surface, 300.0).unwrap();
            near(flux.convection_temperature_tangent_w_m2_k, actual.coefficient.value.value(), 1e-10);
        }
    }
}

#[test]
fn g3_laminar_section_closes_native_conducting_surface_energy_and_entropy() {
    let law = SourceCellConvectionLaw::ForcedLaminarTubeSection {
        speed_m_s: 0.01, start: Length::new(1.0), end: Length::new(2.0),
    };
    let mut config = settings(law);
    config.property_temperature = ConvectionPropertyTemperature::Reservoir;
    source_surface_energy_with(config, |t| section_reference(t, law));
}

#[test]
fn g0_pointwise_exterior_refuses_unverifiable_four_wall_rectangular_geometry() {
    let body = curve();
    let card = liquid_card(0.0, 0.0, 320.0);
    let liquid = liquid_law(&card);
    let mut config = settings(SourceCellConvectionLaw::ForcedRectangularDuctSection {
        speed_m_s: 0.01, aspect_ratio: 0.5,
        start: Length::new(1.0), end: Length::new(2.0),
    });
    config.property_temperature = ConvectionPropertyTemperature::Reservoir;
    assert!(matches!(ImplicitLiquidConvection::try_new(
        &body, [-20000.0, 20000.0], &liquid, config,
    ), Err(LiquidExteriorConvectionError::Invalid(_))));
}
