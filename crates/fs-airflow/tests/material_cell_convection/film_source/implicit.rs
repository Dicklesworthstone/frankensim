//! G1/G3: original liquid curves -> actual nonlinear surface -> native heat history.
use super::*;

#[path = "implicit/sections.rs"]
mod sections;
#[path = "implicit/gas.rs"]
mod gas;
#[path = "implicit/mass_flow.rs"]
mod mass_flow;
#[path = "implicit/turbulent_wall.rs"]
mod turbulent_wall;
#[path = "implicit/laminar_tube.rs"]
mod laminar_tube;
use fs_airflow::graph::thermal::medium::convection::{
    AidingVerticalPlateSettings, HorizontalPlateFacing, NaturalConvectionGeometry,
    evaluate_aiding_vertical_plate_from_material, evaluate_natural_from_material,
    implicit::{ImplicitLiquidConvection, LiquidExteriorConvectionSettings},
};
use fs_matdb::{ClaimSelection, SelectionPolicy};
use fs_material::liquid::{KINEMATIC_VISCOSITY_DIMS, KINEMATIC_VISCOSITY_PROPERTY};
use fs_material::thermal::{
    conduction::LayerConductionConfig,
    transport::{
        conduction::{CaloricConductionLayer, IntegratedCaloricGrayBoundary},
        isobaric::{
            IsobaricHeatCapacityLaw, IsobaricHeatCapacityQuery, conduction::IsobaricConductionLaw,
        },
    },
};
use fs_qty::Pressure;
use fs_thermal::enthalpy::conduction::{
    ConductionObserver, ConductiveReservoir,
    implicit::{
        ImplicitConductionChain, ImplicitConductionControls, ImplicitConductionLaw,
        ImplicitConductionWorkspace,
        boundary::{ImplicitBoundaryGeometry, ImplicitGrayExteriorLaw},
    },
};

fn settings(law: SourceCellConvectionLaw) -> LiquidExteriorConvectionSettings {
    LiquidExteriorConvectionSettings {
        fluid_temperature_bounds_k: [280.0, 320.0],
        surface_temperature_bounds_k: [281.0, 319.0],
        characteristic_length: Length::new(0.1),
        property_temperature: ConvectionPropertyTemperature::FilmMean,
        law,
        emissivity: 0.7,
    }
}

fn check_derivative(law: &ImplicitLiquidConvection<'_, '_>, surface: f64, bath: f64) {
    let response = law.evaluate(surface, bath).unwrap();
    let delta = 1e-3;
    let value = |temperature| {
        law.evaluate(temperature, bath)
            .unwrap()
            .coefficient
            .value
            .value()
    };
    let difference = (value(surface + delta) - value(surface - delta)) / (2.0 * delta);
    near(
        response.coefficient_surface_derivative_w_m2_k2,
        difference,
        2e-7 * difference.abs().max(1.0),
    );
}

fn mixed_settings(speed: f64, gravity: f64) -> AidingVerticalPlateSettings {
    AidingVerticalPlateSettings {
        speed_up_m_s: speed,
        gravity_m_s2: gravity,
        maximum_relative_density_departure: 0.02,
    }
}

#[test]
fn g1_mixed_film_chain_retains_signed_sources_and_pure_flow_limits() {
    let body = curve();
    for density_slope in [0.3, -0.3] {
        let card = liquid_card(density_slope, 0.01, 320.0);
        let liquid = liquid_law(&card)
            .with_density_expansion(Temperature::new(280.0), Temperature::new(320.0))
            .unwrap();
        for surface in [290.0, 310.0] {
            let speed = 0.1 * (density_slope * (surface - 300.0)).signum();
            for (speed, gravity) in [(speed, 9.81), (0.0, 9.81), (speed, 0.0)] {
                let driving = mixed_settings(speed, gravity);
                let convection = ImplicitLiquidConvection::try_new(
                    &body, [-20000.0, 20000.0], &liquid,
                    settings(SourceCellConvectionLaw::AidingVerticalPlate { settings: driving }),
                ).unwrap();
                let actual = convection.evaluate(surface, 300.0).unwrap();
                let original = liquid.at_temperature(Temperature::new(0.5 * (surface + 300.0)))
                    .unwrap();
                let expected = evaluate_aiding_vertical_plate_from_material(
                    &original, driving, original.volumetric_expansion_per_k().unwrap(),
                    Temperature::new(surface), Temperature::new(300.0), Length::new(0.1),
                ).unwrap();
                near(actual.coefficient.value.value(), expected.coefficient().value.value(), 1e-9);
                assert_eq!(actual.nusselt.groups(), expected.nusselt().groups());
                assert_eq!(actual.property_temperature_k, 0.5 * (surface + 300.0));
                check_derivative(&convection, surface, 300.0);
                let flux = convection.flux_state_at_surface(0.0, surface, 300.0).unwrap();
                near(flux.convection_temperature_tangent_w_m2_k,
                    actual.coefficient.value.value()
                        + (surface - 300.0) * actual.coefficient_surface_derivative_w_m2_k2, 1e-9);
            }
        }
    }
}

#[test]
fn g0_mixed_surface_refuses_opposing_flow_and_unadmitted_material_or_density() {
    let body = curve();
    let card = liquid_card(0.3, 0.01, 320.0);
    let unexpanded = liquid_law(&card);
    let liquid = liquid_law(&card)
        .with_density_expansion(Temperature::new(280.0), Temperature::new(320.0))
        .unwrap();
    let mut config = settings(SourceCellConvectionLaw::AidingVerticalPlate {
        settings: mixed_settings(0.1, 9.81),
    });
    assert!(ImplicitLiquidConvection::try_new(
        &body, [-20000.0, 20000.0], &unexpanded, config).is_err());
    let convection = ImplicitLiquidConvection::try_new(
        &body, [-20000.0, 20000.0], &liquid, config).unwrap();
    assert!(convection.evaluate(290.0, 300.0).is_err()); // downward buoyancy, upward stream
    assert!(convection.evaluate(319.0, 330.0).is_err()); // film leaves original source
    config.property_temperature = ConvectionPropertyTemperature::Reservoir;
    assert!(ImplicitLiquidConvection::try_new(
        &body, [-20000.0, 20000.0], &liquid, config).is_err());
    config.property_temperature = ConvectionPropertyTemperature::FilmMean;
    config.law = SourceCellConvectionLaw::AidingVerticalPlate {
        settings: AidingVerticalPlateSettings {
            maximum_relative_density_departure: 0.001,
            ..mixed_settings(0.1, 9.81)
        },
    };
    let bounded = ImplicitLiquidConvection::try_new(
        &body, [-20000.0, 20000.0], &liquid, config).unwrap();
    assert!(bounded.evaluate(310.0, 300.0).is_err());
}

#[test]
fn g1_mean_film_natural_derivatives_keep_signed_expansion_and_full_driving_difference() {
    for slope in [0.3, -0.3] {
        let card = liquid_card(slope, 0.01, 320.0);
        let liquid = liquid_law(&card)
            .with_density_expansion(Temperature::new(280.0), Temperature::new(320.0))
            .unwrap();
        let body = curve();
        for geometry in [
            NaturalConvectionGeometry::VerticalPlate,
            NaturalConvectionGeometry::HorizontalCylinder,
            NaturalConvectionGeometry::HorizontalPlate(HorizontalPlateFacing::Up),
            NaturalConvectionGeometry::HorizontalPlate(HorizontalPlateFacing::Down),
            NaturalConvectionGeometry::Sphere,
        ] {
            let convection = ImplicitLiquidConvection::try_new(
                &body,
                [-20000.0, 20000.0],
                &liquid,
                settings(SourceCellConvectionLaw::Natural {
                    geometry,
                    gravity_m_s2: 9.81,
                }),
            )
            .unwrap();
            assert!(core::ptr::eq(convection.liquid(), &liquid));
            assert_eq!(convection.transport_sources().count(), 4);
            for surface in [290.0, 310.0] {
                let result = convection.evaluate(surface, 300.0).unwrap();
                let temperature = 0.5 * (surface + 300.0);
                let original = liquid
                    .at_temperature(Temperature::new(temperature))
                    .unwrap();
                let expected = evaluate_natural_from_material(
                    geometry,
                    &original,
                    9.81,
                    original.volumetric_expansion_per_k().unwrap(),
                    Temperature::new(surface),
                    Temperature::new(300.0),
                    Length::new(0.1),
                )
                .unwrap();
                assert_eq!(result.property_temperature_k, temperature);
                near(
                    result.coefficient.value.value(),
                    expected.coefficient().value.value(),
                    1e-9,
                );
                assert_eq!(result.nusselt.card().id, expected.nusselt().card().id);
                check_derivative(&convection, surface, 300.0);
            }
        }
        let mut bulk = settings(SourceCellConvectionLaw::VerticalPlate { gravity_m_s2: 9.81 });
        bulk.property_temperature = ConvectionPropertyTemperature::Reservoir;
        let convection =
            ImplicitLiquidConvection::try_new(&body, [-20000.0, 20000.0], &liquid, bulk).unwrap();
        check_derivative(&convection, 310.0, 300.0);
    }
}

#[test]
fn g1_sieder_tate_and_whitaker_only_differentiate_original_wall_viscosity() {
    let card = liquid_card(0.3, 0.01, 320.0);
    let liquid = liquid_law(&card);
    let body = curve();
    for model in [
        SourceCellConvectionLaw::Forced {
            correlation: CorrelationId::SiederTateTurbulent,
            speed_m_s: 2.0,
            settings: ForcedConvectionSettings {
                length_over_hydraulic_diameter: Some(100.0),
                ..Default::default()
            },
        },
        SourceCellConvectionLaw::ForcedSphere { speed_m_s: 0.1 },
    ] {
        let mut config = settings(model);
        config.characteristic_length = Length::new(0.01);
        config.property_temperature = ConvectionPropertyTemperature::Reservoir;
        let law =
            ImplicitLiquidConvection::try_new(&body, [-20000.0, 20000.0], &liquid, config).unwrap();
        let response = law.evaluate(310.0, 300.0).unwrap();
        assert_eq!(response.property_temperature_k, 300.0);
        near(response.transport.thermal_conductivity_w_m_k(), 0.5, 1e-15);
        near(response.transport.specific_heat_j_kg_k(), 4000.0, 1e-11);
        near(
            response.surface_dynamic_viscosity_pa_s.unwrap(),
            0.0009,
            1e-18,
        );
        let htc = response.coefficient.value.value();
        let expected = if matches!(model, SourceCellConvectionLaw::ForcedSphere { .. }) {
            (htc - 2.0 * 0.5 / 0.01) * 0.25 * 0.00001 / 0.0009
        } else {
            htc * 0.14 * 0.00001 / 0.0009
        };
        near(
            response.coefficient_surface_derivative_w_m2_k2,
            expected,
            1e-10,
        );
        check_derivative(&law, 310.0, 300.0);
    }
}

fn viscosity_card(knots: Vec<(f64, f64)>, kinematic: bool) -> MaterialCard {
    let mut claims = ClaimSet::new();
    let viscosity_name = if kinematic {
        KINEMATIC_VISCOSITY_PROPERTY
    } else {
        DYNAMIC_VISCOSITY_PROPERTY
    };
    let viscosity_dims = if kinematic {
        KINEMATIC_VISCOSITY_DIMS
    } else {
        DYNAMIC_VISCOSITY_DIMS
    };
    for (name, dims, values) in [
        (
            DENSITY_PROPERTY,
            Density::DIMS,
            vec![(280.0, 1006.0), (320.0, 994.0)],
        ),
        (
            SPECIFIC_HEAT_CAPACITY_PROPERTY,
            SPECIFIC_HEAT_CAPACITY_DIMS,
            vec![(280.0, 3960.0), (320.0, 4040.0)],
        ),
        (
            THERMAL_CONDUCTIVITY_PROPERTY,
            THERMAL_CONDUCTIVITY_DIMS,
            vec![(280.0, 0.46), (320.0, 0.54)],
        ),
        (viscosity_name, viscosity_dims, knots),
    ] {
        let quantity = if name == THERMAL_CONDUCTIVITY_PROPERTY {
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::ThermalConductivity,
                ValueForm::Static,
            ))
        } else {
            QuantitySpec::dimensional(dims)
        };
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::with_quantity(name, quantity),
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: values,
                    dims,
                },
                validity: ValidityDomain::unconstrained()
                    .with_quantity("T", temperature_quantity(), 280.0, 320.0)
                    .with("p", 100000.0, 100000.0)
                    .with("solute_mass_fraction", 0.2, 0.2),
                interpolation: InterpolationPolicy::LinearInside,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic implicit convection derivative fixture".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic liquid".into(),
            phase: "liquid".into(),
            process: "original tangent fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}

#[test]
fn g1_kinematic_wall_viscosity_includes_density_chain_without_resampling() {
    let card = viscosity_card(vec![(280.0, 1.2e-6), (320.0, 0.8e-6)], true);
    let point = QueryPoint::new()
        .with_quantity("T", temperature_quantity(), 300.0)
        .unwrap()
        .with("p", 100000.0)
        .unwrap()
        .with("solute_mass_fraction", 0.2)
        .unwrap();
    let liquid = LiquidTemperatureLaw::try_new(
        &card,
        &point,
        "T",
        LiquidViscosityBasis::Kinematic,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let body = curve();
    let mut config = settings(SourceCellConvectionLaw::ForcedSphere { speed_m_s: 0.1 });
    config.characteristic_length = Length::new(0.01);
    config.property_temperature = ConvectionPropertyTemperature::Reservoir;
    let law =
        ImplicitLiquidConvection::try_new(&body, [-20000.0, 20000.0], &liquid, config).unwrap();
    let result = law.evaluate(310.0, 300.0).unwrap();
    let wall_mu = 997.0 * 0.9e-6;
    let wall_slope = -0.3 * 0.9e-6 - 997.0 * 1e-8;
    near(
        result.surface_dynamic_viscosity_pa_s.unwrap(),
        wall_mu,
        1e-18,
    );
    let expected = -(result.coefficient.value.value() - 100.0) * 0.25 * wall_slope / wall_mu;
    near(
        result.coefficient_surface_derivative_w_m2_k2,
        expected,
        1e-10,
    );
    check_derivative(&law, 310.0, 300.0);
}

#[test]
fn g0_true_source_kinks_and_original_property_conventions_refuse() {
    let card = viscosity_card(
        vec![(280.0, 0.0012), (300.0, 0.001), (320.0, 0.0005)],
        false,
    );
    let liquid = liquid_law(&card);
    let body = curve();
    let mut config = settings(SourceCellConvectionLaw::Forced {
        correlation: CorrelationId::FlatPlateLaminarAverage,
        speed_m_s: 0.1,
        settings: ForcedConvectionSettings::default(),
    });
    let law =
        ImplicitLiquidConvection::try_new(&body, [-20000.0, 20000.0], &liquid, config).unwrap();
    assert!(law.evaluate(310.0, 290.0).is_err()); // film lies exactly on mu kink
    assert!(law.evaluate(319.0, 330.0).is_err()); // actual film leaves source
    config.property_temperature = ConvectionPropertyTemperature::Reservoir;
    let fixed_bulk =
        ImplicitLiquidConvection::try_new(&body, [-20000.0, 20000.0], &liquid, config).unwrap();
    // A fixed source coordinate needs no invented derivative at its knot.
    assert_eq!(
        fixed_bulk
            .evaluate(310.0, 300.0)
            .unwrap()
            .coefficient_surface_derivative_w_m2_k2,
        0.0
    );
    config.law = SourceCellConvectionLaw::ForcedSphere { speed_m_s: 0.1 };
    config.property_temperature = ConvectionPropertyTemperature::FilmMean;
    assert!(
        ImplicitLiquidConvection::try_new(&body, [-20000.0, 20000.0], &liquid, config).is_err()
    );
    config.property_temperature = ConvectionPropertyTemperature::Reservoir;
    config.fluid_temperature_bounds_k = [280.0, 321.0];
    assert!(
        ImplicitLiquidConvection::try_new(&body, [-20000.0, 20000.0], &liquid, config).is_err()
    );
}

#[test]
fn g1_sphere_equilibrium_has_finite_heat_tangent_without_fabricated_coefficient_partial() {
    let card = liquid_card(0.3, 0.01, 320.0);
    let liquid = liquid_law(&card)
        .with_density_expansion(Temperature::new(280.0), Temperature::new(320.0))
        .unwrap();
    let body = curve();
    let config = settings(SourceCellConvectionLaw::Natural {
        geometry: NaturalConvectionGeometry::Sphere,
        gravity_m_s2: 9.81,
    });
    let law =
        ImplicitLiquidConvection::try_new(&body, [-20000.0, 20000.0], &liquid, config).unwrap();
    assert!(law.state_at_surface(0.0, 300.0, 300.0).is_err());
    let flux = law.flux_state_at_surface(0.0, 300.0, 300.0).unwrap();
    near(flux.coefficients.convection_w_m2_k, 10.0, 1e-14);
    near(flux.convection_temperature_tangent_w_m2_k, 10.0, 1e-14);
    let nonzero = law.evaluate(310.0, 300.0).unwrap();
    let flux = law.flux_state_at_surface(0.0, 310.0, 300.0).unwrap();
    near(
        flux.convection_temperature_tangent_w_m2_k,
        nonzero.coefficient.value.value() + 10.0 * nonzero.coefficient_surface_derivative_w_m2_k2,
        1e-11,
    );
    for (density_slope, gravity) in [(0.3, 0.0), (0.0, 9.81)] {
        let zero_card = liquid_card(density_slope, 0.01, 320.0);
        let zero_liquid = liquid_law(&zero_card)
            .with_density_expansion(Temperature::new(280.0), Temperature::new(320.0))
            .unwrap();
        let zero = ImplicitLiquidConvection::try_new(
            &body,
            [-20000.0, 20000.0],
            &zero_liquid,
            settings(SourceCellConvectionLaw::Natural {
                geometry: NaturalConvectionGeometry::Sphere,
                gravity_m_s2: gravity,
            }),
        )
        .unwrap();
        let conductive = zero.evaluate(310.0, 300.0).unwrap();
        near(conductive.coefficient.value.value(), 10.2, 1e-13);
        near(
            conductive.coefficient_surface_derivative_w_m2_k2,
            0.02,
            1e-14,
        );
        let flux = zero.flux_state_at_surface(0.0, 310.0, 300.0).unwrap();
        near(flux.convection_temperature_tangent_w_m2_k, 10.4, 1e-13);
    }
    let plate = ImplicitLiquidConvection::try_new(
        &body,
        [-20000.0, 20000.0],
        &liquid,
        settings(SourceCellConvectionLaw::VerticalPlate { gravity_m_s2: 9.81 }),
    )
    .unwrap();
    assert!(plate.flux_state_at_surface(0.0, 300.0, 300.0).is_err());
}

#[test]
fn g1_mixed_equilibrium_uses_heat_limits_without_inventing_coefficient_partials() {
    // Pr=0.72 makes the buoyant mixed-coefficient response singular at Ra=0.
    // Its contribution to signed heat still has zero first derivative there.
    let card = viscosity_card(vec![(280.0, 9e-5), (320.0, 9e-5)], false);
    let liquid = liquid_law(&card)
        .with_density_expansion(Temperature::new(280.0), Temperature::new(320.0))
        .unwrap();
    let body = curve();
    for speed in [0.01, 0.0] {
        let convection = ImplicitLiquidConvection::try_new(
            &body, [-20000.0, 20000.0], &liquid,
            settings(SourceCellConvectionLaw::AidingVerticalPlate {
                settings: mixed_settings(speed, 9.81),
            }),
        ).unwrap();
        assert!(convection.state_at_surface(0.0, 300.0, 300.0).is_err());
        let flux = convection.flux_state_at_surface(0.0, 300.0, 300.0).unwrap();
        let re: f64 = 1000.0 * speed * 0.1 / 9e-5;
        let pr: f64 = 9e-5 * 4000.0 / 0.5;
        let expected = 2.0 * 0.339 * re.sqrt() * pr.cbrt()
            / (1.0 + 0.1 * pr.powf(-0.75)).powf(2.0 / 9.0) * 0.5 / 0.1;
        near(flux.coefficients.convection_w_m2_k, expected, 1e-10);
        near(flux.convection_temperature_tangent_w_m2_k, expected, 1e-10);
        let mut previous = f64::INFINITY;
        for delta in [1e-2, 1e-4, 1e-6] {
            // q(delta)/delta = h(delta); the derivative is one-sided because
            // negative Delta T would oppose the positive stream.
            let quotient = convection.evaluate(300.0 + delta, 300.0).unwrap()
                .coefficient.value.value();
            let error = (quotient - expected).abs();
            assert!(error < previous);
            if speed == 0.0 && previous.is_finite() {
                near(error / previous, 0.01_f64.powf(0.25), 1e-4);
            }
            previous = error;
        }
        if speed != 0.0 {
            assert!(previous < 1e-5 * expected);
        }
        assert!(convection.flux_state_at_surface(30000.0, 300.0, 300.0).is_err());
        assert!(convection.flux_state_at_surface(0.0, 330.0, 330.0).is_err());
    }
    // The resting limit must still pass the ORIGINAL card's Pr admission.
    let outside = viscosity_card(vec![(280.0, 0.02), (320.0, 0.02)], false);
    let outside = liquid_law(&outside)
        .with_density_expansion(Temperature::new(280.0), Temperature::new(320.0))
        .unwrap();
    let convection = ImplicitLiquidConvection::try_new(
        &body, [-20000.0, 20000.0], &outside,
        settings(SourceCellConvectionLaw::AidingVerticalPlate {
            settings: mixed_settings(0.0, 9.81),
        }),
    ).unwrap();
    assert!(convection.flux_state_at_surface(0.0, 300.0, 300.0).is_err());
}

// Independent source algebra and absolute-temperature bisection. No production
// source derivative, material integral, surface root or Jacobian is used.
fn reference(body_temperature: f64, mixed: bool) -> (f64, [f64; 2]) {
    let powers = |surface: f64| {
        let t = 0.5 * (surface + 300.0) - 300.0;
        let (rho, cp, k, mu) = (
            1000.0 - 0.3 * t,
            4000.0 + 2.0 * t,
            0.5 + 0.002 * t,
            0.001 * (1.0 - 0.01 * t),
        );
        let re: f64 = rho * 0.1 * 0.1 / mu;
        let pr: f64 = mu * cp / k;
        let nusselt = if mixed {
            let ra = (9.81 * (0.3 / rho) * (surface - 300.0) * 0.1_f64.powi(3)
                / ((mu / rho) * (k / (rho * cp)))).abs();
            let forced = 2.0 * 0.339 * re.sqrt() * pr.cbrt()
                / (1.0 + 0.100 * pr.powf(-0.75)).powf(2.0 / 9.0);
            let natural = 4.0 * 0.503 / 3.0 * ra.powf(0.25)
                / (1.0 + 0.670 * pr.powf(-9.0 / 16.0)).powf(4.0 / 9.0);
            let exponent = 3.5 * pr.powf(0.075);
            (forced.powf(exponent) + natural.powf(exponent)).powf(1.0 / exponent)
        } else {
            0.664 * re.sqrt() * pr.cbrt()
        };
        let htc = nusselt * k / 0.1;
        [
            AREA * htc * (surface - 300.0),
            AREA * 0.7 * 5.670_374_419e-8 * (surface.powi(4) - 305.0_f64.powi(4)),
        ]
    };
    let (mut lo, mut hi) = (body_temperature.min(300.0), body_temperature.max(305.0));
    for _ in 0..90 {
        let surface = 0.5 * (lo + hi);
        let conductive = AREA / 0.001
            * (0.5 * (body_temperature - surface)
                + 0.001 * ((body_temperature - 300.0).powi(2) - (surface - 300.0).powi(2)));
        let p = powers(surface);
        if conductive > p[0] + p[1] {
            lo = surface;
        } else {
            hi = surface;
        }
    }
    let surface = 0.5 * (lo + hi);
    (surface, powers(surface))
}

#[test]
fn g3_source_surface_closes_integrated_half_cell_and_native_backward_euler_energy() {
    source_surface_energy(false);
}

#[test]
fn g3_mixed_surface_closes_integrated_half_cell_and_native_backward_euler_energy() {
    source_surface_energy(true);
}

fn source_surface_energy(mixed: bool) {
    let law = if mixed {
        SourceCellConvectionLaw::AidingVerticalPlate { settings: mixed_settings(0.1, 9.81) }
    } else {
        SourceCellConvectionLaw::Forced {
            correlation: CorrelationId::FlatPlateLaminarAverage,
            speed_m_s: 0.1,
            settings: ForcedConvectionSettings::default(),
        }
    };
    source_surface_energy_with(settings(law), |temperature| reference(temperature, mixed));
}

fn source_surface_energy_with(
    configuration: LiquidExteriorConvectionSettings,
    reference: impl Fn(f64) -> (f64, [f64; 2]),
) {
    source_surface_energy_with_wall_correction(configuration, false, reference);
}

fn source_surface_energy_with_wall_correction(
    configuration: LiquidExteriorConvectionSettings,
    corrected: bool,
    reference: impl Fn(f64) -> (f64, [f64; 2]),
) {
    let card = isobaric_liquid_card(0.3, 0.01, 320.0);
    let liquid = isobaric_liquid_law(&card)
        .with_density_expansion(Temperature::new(280.0), Temperature::new(320.0))
        .unwrap();
    let point = |t| {
        liquid
            .reference_point()
            .clone()
            .with_quantity("p", QuantitySpec::dimensional(Pressure::DIMS), 100000.0)
            .unwrap()
            .with("T", t)
            .unwrap()
    };
    let cp_key = PropertyKey::new(SPECIFIC_HEAT_CAPACITY_PROPERTY, SPECIFIC_HEAT_CAPACITY_DIMS);
    let k_key = PropertyKey::with_quantity(
        THERMAL_CONDUCTIVITY_PROPERTY,
        QuantitySpec::semantic(SemanticType::new(
            QuantityKind::ThermalConductivity,
            ValueForm::Static,
        )),
    );
    let selection = ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly);
    let caloric = IsobaricHeatCapacityLaw::try_new(
        &card,
        IsobaricHeatCapacityQuery {
            specific_heat_cp: &cp_key,
            temperature_axis: "T",
            pressure_axis: "p",
            lower: &point(280.0),
            upper: &point(320.0),
            selection,
        },
        300.0,
    )
    .unwrap();
    let conduction = IsobaricConductionLaw::try_new(&card, &caloric, &k_key, selection).unwrap();
    let layer = CaloricConductionLayer::from_isobaric(&conduction).unwrap();
    let source = ImplicitLiquidConvection::try_new(
        layer.caloric_law(),
        layer.energy_bounds_j_kg(),
        &liquid,
        configuration,
    )
    .unwrap();
    let source = if corrected {
        source.with_gnielinski_wall_correction().unwrap()
    } else {
        source
    };
    let geometry = ImplicitBoundaryGeometry {
        area_m2: AREA,
        half_distance_m: 0.001,
        surface_resistance_m2_k_w: 0.0,
    };
    let boundary = IntegratedCaloricGrayBoundary::new(
        &layer,
        geometry,
        0.0,
        0.0,
        305.0,
        LayerConductionConfig {
            max_iterations: 64,
            absolute_heat_tolerance_w: 1e-10,
            relative_heat_tolerance: 1e-12,
        },
    )
    .unwrap()
    .with_exterior_law(&source)
    .unwrap();
    let energy = |temperature: f64| {
        let delta = temperature - 300.0;
        4000.0 * delta + delta * delta
    };
    let body_temperature = 315.0;
    let result = boundary
        .evaluate(energy(body_temperature), 300.0, &mut || Ok(()))
        .unwrap();
    let (surface, powers) = reference(body_temperature);
    near(result.surface_temperature_k, surface, 2e-9);
    near(result.convection_power_w, powers[0], 1e-8);
    near(result.radiation_power_w, powers[1], 1e-8);
    near(result.surface_residual_w, 0.0, 1e-8);
    let temperature_delta = 1e-3;
    let a = reference(body_temperature + temperature_delta).1;
    let b = reference(body_temperature - temperature_delta).1;
    let derivative = ((a[0] + a[1]) - (b[0] + b[1]))
        / (energy(body_temperature + temperature_delta)
            - energy(body_temperature - temperature_delta));
    near(
        result.boundary.power_enthalpy_derivative_w_kg_j,
        derivative,
        1e-9,
    );

    let cell = CellId(0);
    let material = MaterialId(0);
    let initial = DepositionSimulation::new(
        vec![TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: cell,
                segment_id: 0,
                body_id: 0,
                material_id: material,
                mass_kg: MASS,
                specific_h_j_kg: energy(body_temperature),
                geometry: CellGeometry {
                    start_m: [0.0; 3],
                    end_m: [0.001, 0.0, 0.0],
                    volume_m3: 1e-6,
                    surface_area_m2: AREA,
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 0.0,
                    source_w: 0.0,
                    reservoir: Some(Reservoir {
                        area_m2: AREA,
                        conductance_w_k: 0.0,
                        temperature_k: 300.0,
                    }),
                },
            }),
        }],
        &[(material, layer.caloric_law())],
        DepositionBudget {
            max_cells: 1,
            max_links: 0,
            max_events: 1,
            max_substeps: 8,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 1.0,
            energy_tolerance_j: 1e-10,
        },
        || Ok(()),
    )
    .unwrap()
    .0;
    let observer = ConductionObserver::new(
        &initial,
        &[(material, conduction.conductivity_law())],
        &[],
        &[ConductiveReservoir {
            cell_id: cell,
            half_distance_m: 0.001,
            surface_resistance_m2_k_w: 0.0,
        }],
        || Ok(()),
    )
    .unwrap();
    let solver = ImplicitConductionChain::new(
        &initial,
        &observer,
        &[(material, &conduction as &dyn ImplicitConductionLaw)],
        || Ok(()),
    )
    .unwrap()
    .with_boundary_exchanges(&initial, &[(cell, &boundary)], &[], || Ok(()))
    .unwrap();
    let dt = 0.1;
    let mut workspace = ImplicitConductionWorkspace::new(1).unwrap();
    let (next, report) = solver
        .advance_to(
            &initial,
            dt,
            ImplicitConductionControls {
                maximum_dt_s: dt,
                maximum_halvings: 0,
                cell_tolerance_j: 1e-7,
                total_tolerance_j: 1e-7,
                ..Default::default()
            },
            &mut workspace,
            || Ok(()),
        )
        .unwrap();
    let temperature = next.temperature_k(cell).unwrap();
    let (_, p) = reference(temperature);
    near(
        MASS * (energy(temperature) - energy(body_temperature)) + dt * (p[0] + p[1]),
        0.0,
        1e-7,
    );
    near(report.energy.external_heat_j, -dt * (p[0] + p[1]), 1e-7);
    near(report.energy.residual_j, 0.0, 1e-7);
    near(
        report.entropy.reservoir_entropy_input_j_k,
        -dt * (p[0] / 300.0 + p[1] / 305.0),
        1e-9,
    );
    assert_eq!(report.energy.substeps, 1);
    assert_eq!(initial.time_s(), 0.0);
    assert_eq!(
        next.cells()[&cell].boundary,
        initial.cells()[&cell].boundary
    );
}
