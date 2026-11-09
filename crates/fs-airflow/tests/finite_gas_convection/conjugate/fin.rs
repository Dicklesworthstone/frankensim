//! G0/G1/G3/G4: original source fins on the existing reciprocal gas/body owner.
#[path = "fin/nonlinear.rs"]
mod nonlinear;
#[path = "fin/transient.rs"]
mod transient;
use super::*;
use fs_airflow::graph::thermal::medium::convection::finite_gas::fin::{
    FinTipExposure, SourceFinError, SourceFinLimits, StraightFinArray,
};
use fs_convection::fin::StraightFinGeometry;

fn array(tip: FinTipExposure) -> StraightFinArray {
    StraightFinArray::try_new(
        StraightFinGeometry::try_new(0.04, 0.00002, 0.042).unwrap(),
        20,
        AREA - 20.0 * 0.00002,
        tip,
    )
    .unwrap()
}
fn limits() -> SourceFinLimits {
    SourceFinLimits {
        maximum_transverse_biot: 0.1,
        maximum_relative_conductivity_variation: 1e-10,
    }
}
fn film<'a>(
    model: &'a fs_material::thermochemical_gas::ThermochemicalGasModel,
    geometry: StraightFinArray,
) -> FilmGasSurface<'a> {
    FilmGasSurface {
        model,
        law: FilmGasHeatLaw::Forced(ForcedGasSurface {
            area_m2: geometry.wetted_area_m2(),
            speed_m_s: 0.0,
            ..surface()
        }),
    }
}

#[test]
fn fin_source_retains_original_films_and_solid_with_independent_base_heat_in_both_directions() {
    let model = source::model_transport(1.0, true);
    let gas = model.evaluate(&source::point(350.0, 100000.0)).unwrap();
    let card = solid(false, 100.0);
    let curve = caloric(&card);
    let material = layer(&card, &curve, [250.0, 440.0]);
    for tip in [FinTipExposure::Adiabatic, FinTipExposure::Exposed] {
        let geometry = array(tip);
        let owner = film(&model, geometry)
            .with_prescribed_mass_flow(RATE, FLOW_AREA)
            .unwrap()
            .with_straight_fins(&material, geometry, limits())
            .unwrap();
        assert!(core::ptr::eq(owner.material(), &material));
        for temperature in [400.0, 300.0, 350.0] {
            let sample = owner
                .resolve(&gas, Temperature::new(temperature), || Ok(()))
                .unwrap();
            let source = &sample.source;
            let property = source.film.property_gas.as_ref().unwrap();
            let h = plate_reference(property, gas.state().density);
            let k = 60.0;
            let a = geometry.fin().cross_section_m2();
            let p = geometry.fin().exposed_perimeter_m();
            let length = geometry.fin().length_m();
            let m = (h * p / (k * a)).sqrt();
            let beta = if tip == FinTipExposure::Exposed {
                h / (m * k)
            } else {
                0.0
            };
            let x = m * length;
            let one = k * a * m * (x.sinh() + beta * x.cosh()) / (x.cosh() + beta * x.sinh());
            let total = h * geometry.bare_base_area_m2() + f64::from(geometry.count()) * one;
            near(sample.conductance_w_k, total, 2e-13);
            near(
                source.applied_heat_into_gas_w,
                total * (temperature - 350.0),
                2e-11,
            );
            assert!(sample.conductance_w_k < h * geometry.wetted_area_m2());
            assert_eq!(source.film.interface.gas, gas);
            assert_eq!(source.film.interface.body.value(), temperature);
            assert_eq!(property.source_packs(), gas.source_packs());
            assert_eq!(property.state().pressure, gas.state().pressure);
            near(
                property.state().temperature,
                0.5 * (temperature + 350.0),
                1e-12,
            );
            assert!(Arc::ptr_eq(
                &source.conductivity_source,
                material.source_line()
            ));
            let velocity = source.film.prescribed_mass_flow.unwrap();
            near(
                velocity.speed_m_s,
                RATE / (gas.state().density * FLOW_AREA),
                2e-12,
            );
            near(
                source.tip_temperature_k,
                350.0 + (temperature - 350.0) / (x.cosh() + beta * x.sinh()),
                2e-12,
            );
            if temperature == 350.0 {
                assert_eq!(source.applied_heat_into_gas_w, 0.0);
                assert_eq!(source.tip_temperature_k, 350.0);
            }
        }
        let wrong_film = FilmGasSurface {
            model: &model,
            law: FilmGasHeatLaw::Forced(ForcedGasSurface {
                area_m2: geometry.base_area_m2(),
                ..surface()
            }),
        };
        assert!(
            wrong_film
                .with_straight_fins(&material, geometry, limits())
                .unwrap()
                .resolve(&gas, Temperature::new(400.0), || Ok(()))
                .is_err()
        );
    }
}

fn conductivity_with_interior_value(interior: f64) -> MaterialCard {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS),
            value: PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![
                    (200.0, 1.0),
                    (350.0, 1.0),
                    (375.0, interior),
                    (400.0, 1.0),
                    (600.0, 1.0),
                ],
                dims: THERMAL_CONDUCTIVITY_DIMS,
            },
            validity: ValidityDomain::unconstrained().with("T", 200.0, 600.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::LinearInside,
            observations: vec![],
            provenance: Provenance {
                source: "manufactured interior fin conductivity extremum".into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        })
        .unwrap();
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "manufactured fin".into(),
            phase: "solid".into(),
            process: "test".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}

#[test]
fn original_interior_conductivity_peak_and_tip_domain_cannot_hide_between_valid_endpoints() {
    let model = source::model_transport(1.0, true);
    let gas = model.evaluate(&source::point(350.0, 100000.0)).unwrap();
    let card = conductivity_with_interior_value(2.0);
    let curve = caloric(&card);
    let material = layer(&card, &curve, [250.0, 440.0]);
    let geometry = StraightFinArray::try_new(
        StraightFinGeometry::try_new(0.08, 0.000002, 0.042).unwrap(),
        10,
        AREA - 10.0 * 0.000002,
        FinTipExposure::Adiabatic,
    )
    .unwrap();
    let selection = film(&model, geometry)
        .with_prescribed_mass_flow(RATE, FLOW_AREA)
        .unwrap();
    let allowed = selection
        .with_straight_fins(
            &material,
            geometry,
            SourceFinLimits {
                maximum_relative_conductivity_variation: 1.1,
                ..limits()
            },
        )
        .unwrap();
    let witness = allowed
        .resolve(&gas, Temperature::new(400.0), || Ok(()))
        .unwrap();
    let tip = witness.source.tip_temperature_k;
    assert!(tip >= 350.0 && tip < 350.01);
    let k_base = material
        .conductivity_curve()
        .conductivity_w_m_k(400.0)
        .unwrap();
    let k_tip = material
        .conductivity_curve()
        .conductivity_w_m_k(tip)
        .unwrap();
    assert!((k_base - k_tip).abs() / k_base < 0.001);
    assert_eq!(witness.source.relative_conductivity_variation, 1.0);
    let strict = selection
        .with_straight_fins(
            &material,
            geometry,
            SourceFinLimits {
                maximum_relative_conductivity_variation: 0.2,
                ..limits()
            },
        )
        .unwrap();
    assert!(matches!(
        strict.resolve(&gas, Temperature::new(400.0), || Ok(())),
        Err(SourceFinError::ConductivityVariation { actual: 1.0, .. })
    ));

    // A permissive k-variation allowance cannot hide a loss of transverse
    // isothermality inside the fin. The base and tip have k=1, while the
    // original interior source has k=0.01 and hence one hundred times Bi.
    let valley_card = conductivity_with_interior_value(0.01);
    let valley_curve = caloric(&valley_card);
    let valley_material = layer(&valley_card, &valley_curve, [250.0, 440.0]);
    let base_biot = witness.source.fin.transverse_biot();
    let maximum = 10.0 * base_biot;
    let valley = selection
        .with_straight_fins(
            &valley_material,
            geometry,
            SourceFinLimits {
                maximum_transverse_biot: maximum,
                maximum_relative_conductivity_variation: 1.1,
            },
        )
        .unwrap();
    assert!(maximum > base_biot);
    match valley.resolve(&gas, Temperature::new(400.0), || Ok(())) {
        Err(SourceFinError::TransverseBiot {
            actual,
            maximum: admitted,
        }) => {
            assert_eq!(admitted, maximum);
            assert!(actual > admitted);
            near(actual, 100.0 * base_biot, 2e-15);
        }
        other => panic!("interior source Biot must refuse: {other:?}"),
    }

    let narrow = layer(&card, &curve, [390.0, 440.0]);
    let unsupported = selection
        .with_straight_fins(
            &narrow,
            geometry,
            SourceFinLimits {
                maximum_relative_conductivity_variation: 1.1,
                ..limits()
            },
        )
        .unwrap();
    assert!(matches!(
        unsupported.resolve(&gas, Temperature::new(400.0), || Ok(())),
        Err(SourceFinError::Thermal(
            fs_thermal::enthalpy::DepositionError::MaterialDomain
        ))
    ));
    let mut polls = 0;
    let cancelled = allowed.resolve(&gas, Temperature::new(400.0), || {
        polls += 1;
        if polls == 2 {
            Err(fs_thermal::enthalpy::DepositionError::Cancelled)
        } else {
            Ok(())
        }
    });
    assert!(matches!(
        cancelled,
        Err(SourceFinError::Thermal(
            fs_thermal::enthalpy::DepositionError::Cancelled
        ))
    ));
    let retry = allowed
        .resolve(&gas, Temperature::new(400.0), || Ok(()))
        .unwrap();
    assert_eq!(retry.conductance_w_k, witness.conductance_w_k);
    assert_eq!(retry.source.fin, witness.source.fin);
}

#[test]
fn native_finned_body_closes_reciprocal_energy_with_actual_base_area_and_source_feedback() {
    scope(|cx, gate| {
        let old = IsochoricGasParcel::new(
            source::model_transport(1.0, true),
            source::point(350.0, 100000.0),
            0.001,
        )
        .unwrap();
        let saved = old.clone();
        let boundary = uniform_boundary();
        let mut environment = uniform_environment(&old, &boundary);
        environment.exchange.heat_bracket_j = [-20.0, 20.0];
        let wall = WallConvectionSolve {
            temperature_increment_bracket_k: [-2.0, 2.0],
            temperature_tolerance_k: 1e-8,
            maximum_iterations: 128,
        };
        let geometry = array(FinTipExposure::Exposed);
        let mut heat = Vec::new();
        for scale in [1.0, 100.0] {
            let card = solid(false, scale);
            let curve = caloric(&card);
            let material = layer(&card, &curve, [250.0, 440.0]);
            let owner = film(old.model(), geometry)
                .with_prescribed_mass_flow(RATE, FLOW_AREA)
                .unwrap()
                .with_straight_fins(&material, geometry, limits())
                .unwrap();
            let body = UniformEnthalpyStepInput {
                surface_area_m2: geometry.base_area_m2(),
                duration_s: 0.1,
                ..thermal(&curve, 420.0)
            };
            let next = owner
                .isochoric_body_candidate(cx, &environment, body, wall)
                .unwrap();
            let report = &next.report;
            let feedback = &report.convection.source;
            let endpoint = &feedback.endpoint_source;
            let gas = report.next_gas.state();
            assert_eq!(endpoint.base_temperature_k, next.state.temperature_k());
            assert_eq!(&endpoint.film.interface.gas, gas);
            assert_eq!(endpoint.film.interface.area_m2, geometry.wetted_area_m2());
            near(
                report.convection.coefficient_w_per_m2_k,
                feedback.source.total_conductance_w_k / geometry.base_area_m2(),
                1e-12,
            );
            near(
                report.solid_to_gas_heat_j,
                body.duration_s * endpoint.applied_heat_into_gas_w,
                2.0 * environment.exchange.source_tolerance_j,
            );
            let dh = body.mass_kg
                * (next.state.specific_enthalpy_j_kg() - body.initial.specific_enthalpy_j_kg());
            let du = old.mass_kg()
                * (gas.specific_internal_energy_j_kg()
                    - old.state().specific_internal_energy_j_kg());
            near(dh + du, 0.0, environment.exchange.energy_tolerance_j);
            near(
                du,
                report.solid_to_gas_heat_j,
                environment.exchange.energy_tolerance_j,
            );
            near(
                report.combined_energy_residual_j,
                0.0,
                environment.exchange.energy_tolerance_j,
            );
            assert!(report.entropy.convection_production_j_k >= 0.0);
            assert_eq!(report.radiation_heat_j, 0.0);
            assert!(next.state.temperature_k() < body.initial.temperature_k());
            assert!(gas.state().temperature > old.state().state().temperature);
            assert_eq!(old, saved);
            heat.push(report.solid_to_gas_heat_j);
            assert!(
                owner
                    .resolve_isochoric_coefficient(
                        cx,
                        old.state(),
                        geometry.wetted_area_m2(),
                        Temperature::new(420.0)
                    )
                    .is_err()
            );
            let radiation = LumpedThermalEnvironment {
                transport: LumpedThermalTransport::try_declared(20.0, 0.2).unwrap(),
                ..boundary.clone()
            };
            assert!(
                owner
                    .isochoric_body_candidate(
                        cx,
                        &uniform_environment(&old, &radiation),
                        body,
                        wall
                    )
                    .is_err()
            );
            if scale == 100.0 {
                gate.request();
                assert!(matches!(
                    owner.isochoric_body_candidate(cx, &environment, body, wall),
                    Err(ConvectionError::Thermal(ThermalPistonError::Cancelled))
                ));
                assert_eq!(old, saved);
            }
        }
        assert!(
            heat[1] > heat[0] * 1.2,
            "source k must change actual fin cooling: {heat:?}"
        );
    });
}
