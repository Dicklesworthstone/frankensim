//! G1/G3: original wall properties in developing turbulent convection and heat.
use super::*;
use fs_airflow::graph::thermal::medium::convection::implicit::PrescribedMassFlowConvection;

fn tube_settings(start: f64) -> LiquidExteriorConvectionSettings {
    let mut configuration = settings(SourceCellConvectionLaw::ForcedTubeSection {
        speed_m_s: 0.2,
        start: Length::new(start),
        end: Length::new(2.0),
    });
    configuration.property_temperature = ConvectionPropertyTemperature::Reservoir;
    configuration
}

// Independent corrected Gnielinski section and original synthetic liquid
// source algebra. No production section, material evaluator or tangent.
fn reference_coefficient(surface: f64, start: f64, kinematic: bool) -> (f64, f64) {
    let reynolds: f64 = 20_000.0;
    let prandtl: f64 = 8.0;
    let friction = (1.8 * reynolds.log10() - 1.5).powi(-2);
    let developed = friction / 8.0 * (reynolds - 1000.0) * prandtl
        / (1.0 + 12.7 * (friction / 8.0).sqrt() * (prandtl.powf(2.0 / 3.0) - 1.0));
    let a = (start / 0.1).cbrt();
    let b = 20.0_f64.cbrt();
    let baseline = developed * (1.0 + 1.0 / (a * a + a * b + b * b)) * 5.0;
    let offset = surface - 300.0;
    let rho = 1000.0 - 0.3 * offset;
    let cp = 4000.0 + 2.0 * offset;
    let k = 0.5 + 0.002 * offset;
    let (mu, logarithmic_mu_slope) = if kinematic {
        let nu = 1e-6 - 1e-8 * offset;
        (rho * nu, -0.3 / rho - 1e-8 / nu)
    } else {
        (0.001 * (1.0 - 0.01 * offset), -0.01 / (1.0 - 0.01 * offset))
    };
    let wall_prandtl = mu * cp / k;
    let coefficient = baseline * (prandtl / wall_prandtl).powf(0.11);
    let derivative = -0.11 * coefficient * (logarithmic_mu_slope + 2.0 / cp - 0.002 / k);
    (coefficient, derivative)
}

#[test]
fn g1_liquid_wall_prandtl_and_complete_tangent_reach_the_original_section() {
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
        for start in [0.0, 1.0] {
            let uncorrected = ImplicitLiquidConvection::try_new(
                &body, [-20000.0, 20000.0], &liquid, tube_settings(start),
            ).unwrap();
            let baseline = uncorrected.evaluate(310.0, 300.0).unwrap();
            assert!(baseline.surface_transport.is_none());
            assert_eq!(baseline.coefficient_surface_derivative_w_m2_k2, 0.0);
            let corrected = uncorrected.with_gnielinski_wall_correction().unwrap();
            assert!(core::ptr::eq(corrected.caloric_law(), &body as &dyn ThermalLaw));
            for surface in [290.0, 300.0, 310.0] {
                let actual = corrected.evaluate(surface, 300.0).unwrap();
                let (expected, slope) = reference_coefficient(surface, start, kinematic);
                near(actual.coefficient.value.value(), expected, 2e-12 * expected);
                near(actual.coefficient_surface_derivative_w_m2_k2, slope, 2e-11);
                assert_eq!(actual.property_temperature_k, 300.0);
                assert_eq!(actual.transport, baseline.transport);
                let wall = actual.surface_transport.unwrap();
                let original = liquid.at_temperature(Temperature::new(surface))
                    .unwrap().transport_state().unwrap();
                for (prepared, source) in [
                    (wall.density_kg_m3(), original.density_kg_m3()),
                    (wall.specific_heat_j_kg_k(), original.specific_heat_j_kg_k()),
                    (wall.thermal_conductivity_w_m_k(), original.thermal_conductivity_w_m_k()),
                    (wall.dynamic_viscosity_pa_s(), original.dynamic_viscosity_pa_s()),
                ] {
                    near(prepared, source, 2e-14 * source);
                }
                let flux = corrected.flux_state_at_surface(0.0, surface, 300.0).unwrap();
                near(flux.convection_temperature_tangent_w_m2_k,
                    expected + (surface - 300.0) * slope, 2e-11 * expected);
                let details = actual.tube_section.unwrap();
                assert_eq!(details.start.value(), start);
                assert_eq!(details.end.value(), 2.0);
                assert_eq!(details.upstream.is_some(), start > 0.0);
                near(details.nusselt * 5.0, expected, 2e-12 * expected);
                check_derivative(&corrected, surface, 300.0);
            }
        }
    }
}

#[test]
fn g3_prescribed_mass_flow_retains_wall_properties_and_both_independent_tangents() {
    let body = curve();
    let card = liquid_card(0.3, 0.01, 320.0);
    let liquid = liquid_law(&card);
    let source = ImplicitLiquidConvection::try_new(
        &body, [-20000.0, 20000.0], &liquid, tube_settings(1.0),
    ).unwrap().with_gnielinski_wall_correction().unwrap();
    let flow = PrescribedMassFlowConvection::for_liquid(&source, 2.0, 0.01).unwrap();
    for surface in [290.0, 310.0] {
        let result = flow.evaluate(surface, 300.0).unwrap();
        let (expected, derivative) = reference_coefficient(surface, 1.0, false);
        near(result.speed_m_s, 0.2, 1e-15);
        near(result.convection.coefficient.value.value(), expected, 2e-12 * expected);
        near(result.convection.coefficient_surface_derivative_w_m2_k2, derivative, 2e-11);
        let delta = 1e-4;
        let at = |rate| PrescribedMassFlowConvection::for_liquid(&source, rate, 0.01)
            .unwrap().evaluate(surface, 300.0).unwrap().convection.coefficient.value.value();
        let numerical = (at(2.0 + delta) - at(2.0 - delta)) / (2.0 * delta);
        near(result.coefficient_mass_flow_derivative_w_m2_k_s_kg.unwrap(), numerical,
            2e-7 * numerical.abs());
        let flux = flow.flux_state_at_surface(0.0, surface, 300.0).unwrap();
        near(flux.convection_temperature_tangent_w_m2_k,
            expected + (surface - 300.0) * derivative, 2e-11 * expected);
    }
}

#[test]
fn g0_wall_correction_refuses_wrong_law_missing_support_and_original_source_kinks() {
    let body = curve();
    let card = liquid_card(0.3, 0.01, 320.0);
    let liquid = liquid_law(&card);
    let mut configuration = tube_settings(1.0);
    configuration.surface_temperature_bounds_k = [281.0, 330.0];
    assert!(ImplicitLiquidConvection::try_new(
        &body, [-20000.0, 20000.0], &liquid, configuration,
    ).unwrap().with_gnielinski_wall_correction().is_err());
    let mut configuration = tube_settings(1.0);
    configuration.law = SourceCellConvectionLaw::ForcedLaminarTubeSection {
        speed_m_s: 0.001, start: Length::new(1.0), end: Length::new(2.0),
    };
    assert!(ImplicitLiquidConvection::try_new(
        &body, [-20000.0, 20000.0], &liquid, configuration,
    ).unwrap().with_gnielinski_wall_correction().is_err());

    let kinked = viscosity_card(
        vec![(280.0, 0.0012), (300.0, 0.001), (320.0, 0.0009)], false,
    );
    let liquid = liquid_law(&kinked);
    let source = ImplicitLiquidConvection::try_new(
        &body, [-20000.0, 20000.0], &liquid, tube_settings(1.0),
    ).unwrap();
    assert!(source.evaluate(300.0, 300.0).is_ok());
    let corrected = source.with_gnielinski_wall_correction().unwrap();
    assert!(corrected.evaluate(300.0, 300.0).is_err());
    assert!(corrected.evaluate(310.0, 300.0).is_ok());
    assert!(corrected.evaluate(330.0, 300.0).is_err());
    assert!(corrected.state_at_surface(30000.0, 310.0, 300.0).is_err());
}

// Independent nonlinear wall balance with a finite original conductive drop.
fn surface_reference(body_temperature: f64) -> (f64, [f64; 2]) {
    let powers = |surface: f64| [
        AREA * reference_coefficient(surface, 1.0, false).0 * (surface - 300.0),
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
fn g3_corrected_surface_closes_native_conduction_heat_storage_and_bath_entropy() {
    source_surface_energy_with_wall_correction(tube_settings(1.0), true, surface_reference);
}
