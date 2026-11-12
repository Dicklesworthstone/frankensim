//! G1/G3/G4: original conducting surfaces and finite gas/body heat ownership.
#[path = "conjugate/specialized.rs"]
mod specialized;
#[path = "conjugate/window.rs"]
mod window;
use super::*;
use fs_airflow::graph::thermal::medium::convection::finite_gas::{
    conjugate::ConjugateGasError,
    film::{FilmGasHeatLaw, FilmGasSurface},
};
use fs_couple::finite_gas_thermal::{WallConvectionError, WallConvectionSolve};
use fs_evidence::ValidityDomain;
use fs_matdb::*;
use fs_material::{
    fluid_transport::FluidTransport,
    phase::{EnthalpyPhaseKnot, SolidLiquidPhase},
    state_point::THERMAL_CONDUCTIVITY_DIMS,
    thermal::{
        conduction::{
            LayerConductionConfig, ReservoirLayerGeometry, source::ResolvedConductionLayer,
        },
        source::ThermalMaterialPoint,
    },
};
use fs_qty::Temperature;
use std::sync::Arc;

const AREA: f64 = 0.02;
const DEPTH: f64 = 0.0005;
const RATE: f64 = 0.02;
const FLOW_AREA: f64 = 0.01;

fn config() -> LayerConductionConfig {
    LayerConductionConfig {
        max_iterations: 96,
        absolute_heat_tolerance_w: 1e-11,
        relative_heat_tolerance: 1e-13,
    }
}

fn solid(variable: bool, scale: f64) -> MaterialCard {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS),
            value: PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: [200.0, 350.0, 600.0]
                    .map(|t| (t, scale * if variable { 0.002 * t } else { 0.6 }))
                    .to_vec(),
                dims: THERMAL_CONDUCTIVITY_DIMS,
            },
            validity: ValidityDomain::unconstrained().with("T", 200.0, 600.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::LinearInside,
            observations: vec![],
            provenance: Provenance {
                source: "manufactured conducting gas wall; not measured".into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        })
        .unwrap();
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic conducting wall".into(),
            phase: "equilibrium phase".into(),
            process: "manufactured".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}

fn caloric(card: &MaterialCard) -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_new(
        card.content_hash(),
        [
            (0.0, 200.0, 0.0),
            (250000.0, 450.0, 0.0),
            (310000.0, 450.0, 1.0),
            (460000.0, 600.0, 1.0),
        ]
        .map(|(h, t, fraction)| EnthalpyPhaseKnot {
            specific_enthalpy_j_kg: h,
            temperature_k: t,
            liquid_mass_fraction: fraction,
            bulk_density_kg_m3: 1000.0,
        })
        .to_vec(),
    )
    .unwrap()
}

fn layer<'a>(
    card: &MaterialCard,
    curve: &'a EquilibriumEnthalpyPhaseCurve,
    span: [f64; 2],
) -> ResolvedConductionLayer<'a> {
    let context = QueryPoint::new().with("T", span[0]).unwrap();
    ResolvedConductionLayer::new(
        curve,
        ThermalMaterialPoint {
            material: card,
            context: &context,
            temperature_axis: "T",
            liquid_fraction_axis: None,
        },
        &PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS),
        ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
        SolidLiquidPhase::Solid,
        span,
        || Ok(()),
    )
    .unwrap()
}

fn thermal(
    curve: &EquilibriumEnthalpyPhaseCurve,
    temperature: f64,
) -> UniformEnthalpyStepInput<'_> {
    UniformEnthalpyStepInput {
        initial: curve
            .state_at_specific_enthalpy(1000.0 * (temperature - 200.0))
            .unwrap(),
        mass_kg: 0.02,
        volume_m3: 2.0 * AREA * DEPTH,
        surface_area_m2: AREA,
        duration_s: 0.005,
        ..body(curve, 0.0)
    }
}

fn plate(model: &fs_material::thermochemical_gas::ThermochemicalGasModel) -> FilmGasSurface<'_> {
    FilmGasSurface {
        model,
        law: FilmGasHeatLaw::Forced(ForcedGasSurface {
            speed_m_s: 0.0,
            ..surface()
        }),
    }
}

// Original flat-plate equation, using the solved surface film and ORIGINAL
// bulk-density velocity; no production evaluator supplies the reference.
fn plate_reference(film: &impl FluidTransport, bulk_density: f64) -> f64 {
    let p = film.transport_state().unwrap();
    let re =
        p.density_kg_m3() * (RATE / (bulk_density * FLOW_AREA)) * 0.1 / p.dynamic_viscosity_pa_s();
    0.664 * re.sqrt() * p.prandtl().cbrt() * p.thermal_conductivity_w_m_k() / 0.1
}

#[test]
fn original_radial_layer_recovers_analytic_series_heat_and_surface_in_both_directions() {
    let gas_model = source::model_transport(1.0, true);
    let bulk = gas_model.evaluate(&source::point(350.0, 100000.0)).unwrap();
    let geometry = ReservoirLayerGeometry::cylindrical(0.005, 0.006, 0.5).unwrap();
    let area = geometry.surface_area_m2();
    let film = FilmGasSurface {
        model: &gas_model,
        law: FilmGasHeatLaw::BulkForced(ForcedGasSurface {
            area_m2: area,
            characteristic_length: Length::new(0.01),
            speed_m_s: 0.1,
            correlation: CorrelationId::CircularDuctLaminarCwt,
            settings: ForcedConvectionSettings {
                length_over_hydraulic_diameter: Some(50.0),
                ..ForcedConvectionSettings::default()
            },
        }),
    };
    let g_film = 3.66 * bulk.transport_state().unwrap().thermal_conductivity_w_m_k() / 0.01 * area;
    let multiplier = core::f64::consts::TAU * 0.5 / (0.006_f64 / 0.005).ln();
    for variable in [false, true] {
        let card = solid(variable, 0.01);
        let curve = caloric(&card);
        let wall = layer(&card, &curve, [250.0, 450.0]);
        let owner = film
            .with_conducting_layer(&wall, geometry, config())
            .unwrap();
        for temperature in [400.0, 300.0, 350.0] {
            let body = thermal(&curve, temperature).initial;
            let resolved = owner.resolve(&bulk, body, || Ok(())).unwrap();
            let sample = &resolved.source;
            let response = sample.layer.as_ref().unwrap();
            let (a, b) = if variable {
                (0.0, 0.00002)
            } else {
                (0.006, 0.0)
            };
            let linear = multiplier * (a + b * temperature) + g_film;
            let difference = temperature - 350.0;
            // Stable small root of C[k(Tb)x-b*x²/2] = Gfilm[(Tb-Tg)-x].
            let drop = 2.0 * g_film * difference
                / (linear + (linear * linear - 2.0 * multiplier * b * g_film * difference).sqrt());
            let surface = temperature - drop;
            let power = g_film * (surface - 350.0);
            near(response.surface_temperature_k, surface, 2e-9);
            near(response.power_to_bath_w, power, 2e-10);
            near(
                response.solid_power_to_bath_w,
                multiplier
                    * (a * (temperature - surface)
                        + 0.5 * b * (temperature - surface) * (temperature + surface)),
                2e-10,
            );
            near(sample.applied_heat_into_gas_w, power, 2e-10);
            assert_eq!(sample.body, body);
            assert_eq!(response.film.interface.gas, bulk);
            assert_eq!(
                response.film.interface.body.value(),
                response.surface_temperature_k
            );
            assert_eq!(
                response
                    .film
                    .property_gas
                    .as_ref()
                    .unwrap()
                    .state()
                    .temperature,
                350.0
            );
            assert_eq!(
                response.film.property_gas.as_ref().unwrap().source_packs(),
                bulk.source_packs()
            );
            assert!(Arc::ptr_eq(&sample.conductivity_source, wall.source_line()));
            assert!(sample.insulated_source.is_none());
            if difference == 0.0 {
                let solid_g = multiplier * (a + b * temperature);
                near(
                    resolved.conductance_w_k,
                    solid_g * g_film / (solid_g + g_film),
                    1e-12,
                );
                assert_eq!(response.power_to_bath_w, 0.0);
            } else {
                assert!((surface - temperature).abs() > 1.0);
                assert!(resolved.conductance_w_k < g_film);
            }
        }
        // The same 450 K is not sufficient: a latent enthalpy leaves this pure-solid path.
        let molten = curve.state_at_specific_enthalpy(280000.0).unwrap();
        assert!(matches!(
            owner.resolve(&bulk, molten, || Ok(())),
            Err(ConjugateGasError::Layer(_))
        ));
        let narrow = layer(&card, &curve, [300.0, 380.0]);
        let narrow_owner = film
            .with_conducting_layer(&narrow, geometry, config())
            .unwrap();
        assert!(
            narrow_owner
                .resolve(&bulk, thermal(&curve, 380.0).initial, || Ok(()))
                .is_ok()
        );
        assert!(
            narrow_owner
                .resolve(&bulk, thermal(&curve, 400.0).initial, || Ok(()))
                .is_err()
        );
    }

    // A source-pinned h endpoint may round to a native T just outside support.
    // Inactive area still admits the participating body; no active root may
    // mask this shared initial/midpoint/final-state check.
    let card = solid(true, 0.01);
    let shifted = EquilibriumEnthalpyPhaseCurve::try_single_phase(
        card.content_hash(),
        SolidLiquidPhase::Solid,
        [(1e7, 280.0), (1e7 + 40000.0, 320.0)]
            .map(|(h, t)| EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: h,
                temperature_k: t,
                liquid_mass_fraction: 0.0,
                bulk_density_kg_m3: 1000.0,
            })
            .to_vec(),
    )
    .unwrap();
    let lower = 305.23456789;
    let pinned = layer(&card, &shifted, [lower, 320.0]);
    let h = pinned.enthalpy_range_j_kg()[0];
    let state = shifted.state_at_specific_enthalpy(h).unwrap();
    assert_eq!(pinned.temperature_at_enthalpy(h).unwrap(), lower);
    assert!(state.temperature_k() < lower);
    let insulated = FilmGasSurface {
        model: &gas_model,
        law: FilmGasHeatLaw::Forced(ForcedGasSurface {
            area_m2: 0.0,
            speed_m_s: 0.0,
            ..surface()
        }),
    }
    .with_conducting_layer(&pinned, geometry, config())
    .unwrap();
    assert!(matches!(
        insulated.resolve(&bulk, state, || Ok(())),
        Err(ConjugateGasError::Thermal(
            fs_thermal::enthalpy::DepositionError::MaterialDomain
        ))
    ));
}

#[test]
fn conducting_surface_refreshes_mass_flux_film_and_closes_original_piston_energy() {
    let card = solid(true, 1.0);
    let curve = caloric(&card);
    let wall = layer(&card, &curve, [300.0, 440.0]);
    let geometry = ReservoirLayerGeometry::planar(AREA, DEPTH).unwrap();
    let old = GasPiston::new(
        ClosedGasParcel::new(
            source::model_transport(1.0, true),
            source::point(350.0, 100000.0),
            0.001,
        )
        .unwrap(),
        1.0,
        0.01,
        0.0,
    )
    .unwrap();
    let saved = old.clone();
    let film = plate(old.gas().model());
    let stream = film.with_prescribed_mass_flow(RATE, FLOW_AREA).unwrap();
    let owner = stream
        .with_conducting_layer(&wall, geometry, config())
        .unwrap();
    let body = UniformEnthalpyStepInput {
        // The free piston's gas spring, sqrt(rho c^2 A^2 / (m V)) = 118.3 rad/s
        // (1 kg, 0.01 m2, 1 L of ideal gas at 1 bar, gamma 1.4), is admitted by fs-phs only for
        // rate*dt <= 0.04 (MIDPOINT_ONE_CENT_ANGULAR_STEP): dt <= 3.38e-4 s. 0.3 ms gives 0.0355.
        duration_s: 0.0003,
        ..thermal(&curve, 420.0)
    };
    let piston_input = PistonInput {
        duration_s: body.duration_s,
        ..input()
    };
    let heat = HeatExchangeSolve {
        heat_bracket_j: [-1.0, 1.0],
        ..exchange()
    };
    let next = owner
        .finite_body_candidate(&old, piston_input, body, controls(), heat, || false)
        .unwrap();
    let report = next.step();
    let sample = &next.thermal().report;
    let response = sample.layer.as_ref().unwrap();
    let bulk = &response.film.interface.gas;
    let property = response.film.property_gas.as_ref().unwrap();
    let velocity = response.film.prescribed_mass_flow.unwrap();
    assert_eq!(bulk, &report.piston.midpoint);
    assert_eq!(sample.body, report.thermal_midpoint);
    assert_eq!(
        response.film.interface.body.value(),
        response.surface_temperature_k
    );
    assert!(sample.body_temperature_k > response.surface_temperature_k);
    assert!((sample.body_temperature_k - response.surface_temperature_k).abs() > 0.1);
    near(
        property.state().temperature,
        0.5 * (response.surface_temperature_k + bulk.state().temperature),
        1e-12,
    );
    assert_eq!(property.state().pressure, bulk.state().pressure);
    assert_eq!(property.source_packs(), bulk.source_packs());
    assert_eq!(velocity.bulk_density_kg_m3, bulk.state().density);
    assert!((property.state().density / bulk.state().density - 1.0).abs() > 0.05);
    near(
        velocity.speed_m_s * bulk.state().density * FLOW_AREA,
        RATE,
        2e-15,
    );
    let h = plate_reference(property, bulk.state().density);
    let film_power = AREA * h * (response.surface_temperature_k - bulk.state().temperature);
    let solid_power = (AREA / DEPTH)
        * 0.001
        * (sample.body_temperature_k - response.surface_temperature_k)
        * (sample.body_temperature_k + response.surface_temperature_k);
    near(response.power_to_bath_w, film_power, 2e-10);
    near(response.solid_power_to_bath_w, solid_power, 2e-9);
    near(
        report.exchanged_heat_j,
        body.duration_s * film_power,
        heat.source_tolerance_j,
    );
    assert!((body.duration_s * response.heat_residual_w).abs() <= heat.source_tolerance_j);
    let du = old.gas().mass_kg()
        * (next.piston().gas().state().specific_internal_energy_j_kg()
            - old.gas().state().specific_internal_energy_j_kg());
    let dh = body.mass_kg
        * (next.thermal().state.specific_enthalpy_j_kg() - body.initial.specific_enthalpy_j_kg());
    let kinetic = 0.5 * old.mass_kg() * next.piston().velocity_m_s().powi(2);
    let exterior = piston_input.exterior_pressure_pa
        * (next.piston().gas().volume_m3() - old.gas().volume_m3());
    near(du + dh + kinetic + exterior, 0.0, heat.energy_tolerance_j);
    near(dh + report.exchanged_heat_j, 0.0, heat.energy_tolerance_j);
    assert_eq!(next.piston().gas().mass_kg(), old.gas().mass_kg());
    assert!(Arc::ptr_eq(&sample.conductivity_source, wall.source_line()));

    let same_numbers_other_object = curve.clone();
    assert!(
        owner
            .finite_body_candidate(
                &old,
                piston_input,
                UniformEnthalpyStepInput {
                    curve: &same_numbers_other_object,
                    ..body
                },
                controls(),
                heat,
                || false,
            )
            .is_err()
    );
    let wrong_model = source::model_transport(1.0, false);
    let wrong = plate(&wrong_model)
        .with_prescribed_mass_flow(RATE, FLOW_AREA)
        .unwrap()
        .with_conducting_layer(&wall, geometry, config())
        .unwrap();
    assert!(
        wrong
            .resolve(old.gas().state(), body.initial, || Ok(()))
            .is_err()
    );
    let mismatched = stream
        .with_conducting_layer(
            &wall,
            ReservoirLayerGeometry::planar(0.5 * AREA, DEPTH).unwrap(),
            config(),
        )
        .unwrap();
    assert!(
        mismatched
            .resolve(old.gas().state(), body.initial, || Ok(()))
            .is_err()
    );

    // Every energy midpoint is supported, including the zero-Q seed at 400 K;
    // only the proposed final body leaves this layer's 400 K source endpoint.
    let limited = layer(&card, &curve, [300.0, 400.0]);
    let limited_owner = stream
        .with_conducting_layer(&limited, geometry, config())
        .unwrap();
    let warming = UniformEnthalpyStepInput {
        internal_heat_j: 4.0,
        duration_s: body.duration_s,
        ..thermal(&curve, 399.9)
    };
    assert!(matches!(
        limited_owner.finite_body_candidate(
            &old,
            piston_input,
            warming,
            controls(),
            HeatExchangeSolve {
                heat_bracket_j: [0.0, 1.0],
                ..heat
            },
            || false,
        ),
        Err(ThermalPistonError::Interface(ConjugateGasError::Layer(_)))
    ));

    let isolated = FilmGasSurface {
        model: &wrong_model,
        law: FilmGasHeatLaw::Forced(ForcedGasSurface {
            area_m2: 0.0,
            speed_m_s: 0.0,
            ..surface()
        }),
    }
    .with_prescribed_mass_flow(RATE, FLOW_AREA)
    .unwrap()
    .with_conducting_layer(&wall, geometry, config())
    .unwrap()
    .finite_body_candidate(&old, piston_input, body, controls(), heat, || false)
    .unwrap();
    assert_eq!(isolated.step().exchanged_heat_j, 0.0);
    assert_eq!(isolated.thermal().state, body.initial);
    assert!(isolated.thermal().report.layer.is_none());
    let zero_source = isolated.thermal().report.insulated_source.as_ref().unwrap();
    assert!(zero_source.property_gas.is_none());
    assert!(zero_source.interface.evaluation.is_none());

    let mut calls = 0;
    assert!(matches!(
        owner.finite_body_candidate(&old, piston_input, body, controls(), heat, || {
            calls += 1;
            calls == 8
        },),
        Err(ThermalPistonError::Cancelled)
    ));
    assert_eq!(old, saved);
}

#[test]
fn isochoric_conjugate_heat_keeps_native_energy_entropy_biot_and_surface_ownership() {
    scope(|cx, _| {
        let card = solid(true, 1.0);
        let curve = caloric(&card);
        let wall = layer(&card, &curve, [300.0, 440.0]);
        let geometry = ReservoirLayerGeometry::planar(AREA, DEPTH).unwrap();
        let old = IsochoricGasParcel::new(
            source::model_transport(1.0, true),
            source::point(350.0, 100000.0),
            0.003,
        )
        .unwrap();
        let saved = old.clone();
        let transport = LumpedThermalTransport::from_material_card_without_radiation(
            &card,
            &PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS),
            "T",
            &[
                QueryPoint::new().with("T", 200.0).unwrap(),
                QueryPoint::new().with("T", 600.0).unwrap(),
            ],
            SelectionPolicy::SingleClaimOnly,
        )
        .unwrap();
        assert!(transport.radiation_is_disabled());
        let boundary = LumpedThermalEnvironment {
            transport,
            ..uniform_boundary()
        };
        let mut environment = uniform_environment(&old, &boundary);
        environment.exchange.heat_bracket_j = [-4.0, 4.0];
        let body = UniformEnthalpyStepInput {
            duration_s: 0.1,
            ..thermal(&curve, 420.0)
        };
        let solve = WallConvectionSolve {
            temperature_increment_bracket_k: [-2.0, 2.0],
            temperature_tolerance_k: 1e-8,
            maximum_iterations: 128,
        };
        let owner = plate(old.model())
            .with_prescribed_mass_flow(RATE, FLOW_AREA)
            .unwrap()
            .with_conducting_layer(&wall, geometry, config())
            .unwrap();
        let next = owner
            .isochoric_body_candidate(cx, &environment, body, solve)
            .unwrap();
        let report = &next.report;
        let feedback = &report.convection.source;
        let sample = &feedback.source;
        let response = sample.layer.as_ref().unwrap();
        let endpoint_sample = &feedback.endpoint_source;
        let endpoint_surface = endpoint_sample.layer.as_ref().unwrap();
        let bulk = report.next_gas.state();
        let property = response.film.property_gas.as_ref().unwrap();
        let velocity = response.film.prescribed_mass_flow.unwrap();
        assert_eq!(&response.film.interface.gas, bulk);
        assert_eq!(sample.body.value(), feedback.temperature_k);
        assert_eq!(sample.body_temperature_k, feedback.temperature_k);
        assert_eq!(
            endpoint_sample.body_temperature_k,
            next.state.temperature_k()
        );
        assert!(response.surface_temperature_k < sample.body_temperature_k);
        assert!(response.surface_temperature_k > bulk.state().temperature);
        near(
            property.state().temperature,
            0.5 * (response.surface_temperature_k + bulk.state().temperature),
            1e-12,
        );
        assert_eq!(property.source_packs(), bulk.source_packs());
        assert_eq!(property.state().pressure, bulk.state().pressure);
        assert_eq!(velocity.bulk_density_kg_m3, bulk.state().density);
        near(
            velocity.bulk_density_kg_m3,
            old.state().state().density,
            2e-12,
        );
        near(
            velocity.speed_m_s,
            RATE / (old.state().state().density * FLOW_AREA),
            2e-12,
        );
        assert!(Arc::ptr_eq(&sample.conductivity_source, wall.source_line()));
        let h = plate_reference(property, bulk.state().density);
        near(
            response.power_to_bath_w,
            AREA * h * (response.surface_temperature_k - bulk.state().temperature),
            2e-10,
        );
        near(
            report.convection.coefficient_w_per_m2_k,
            response.conductance_w_k / AREA,
            1e-12,
        );
        // The endpoint reevaluation is at the thermal owner's actual body T.
        near(
            report.solid_to_gas_heat_j,
            body.duration_s * endpoint_surface.power_to_bath_w,
            2.0 * environment.exchange.source_tolerance_j,
        );
        near(
            report.solid_to_gas_heat_j,
            body.duration_s
                * AREA
                * report.convection.coefficient_w_per_m2_k
                * (next.state.temperature_k() - bulk.state().temperature),
            environment.exchange.source_tolerance_j,
        );
        let du = old.mass_kg()
            * (bulk.specific_internal_energy_j_kg() - old.state().specific_internal_energy_j_kg());
        let dh = body.mass_kg
            * (next.state.specific_enthalpy_j_kg() - body.initial.specific_enthalpy_j_kg());
        near(du, report.solid_to_gas_heat_j, 1e-7);
        near(dh, -report.solid_to_gas_heat_j, 1e-7);
        near(report.combined_energy_residual_j, 0.0, 1e-7);
        let q = report.solid_to_gas_heat_j;
        let ts = endpoint_surface.surface_temperature_k;
        let solid_entropy = q * (1.0 / ts - 1.0 / next.state.temperature_k());
        let film_entropy = q * (1.0 / bulk.state().temperature - 1.0 / ts);
        assert!(solid_entropy >= 0.0 && film_entropy >= 0.0);
        near(
            report.entropy.convection_production_j_k,
            solid_entropy + film_entropy,
            1e-12,
        );
        near(report.entropy.balance_residual_j_k, 0.0, 1e-12);
        assert_eq!(report.radiation_heat_j, 0.0);
        assert_eq!(report.next_gas.mass_kg(), old.mass_kg());
        assert!(next.state.temperature_k() < body.initial.temperature_k());
        assert!(bulk.state().temperature > old.state().state().temperature);

        let radiating = LumpedThermalEnvironment {
            transport: LumpedThermalTransport::try_declared(0.6, 0.2).unwrap(),
            ..boundary.clone()
        };
        assert!(!radiating.transport.radiation_is_disabled());
        assert!(matches!(
            owner
                .isochoric_body_candidate(cx, &uniform_environment(&old, &radiating), body, solve,),
            Err(ConvectionError::Convection(WallConvectionError::Source(
                ConjugateGasError::Invalid(_),
            )))
        ));
        let poor_core = LumpedThermalEnvironment {
            transport: LumpedThermalTransport::try_declared(1e-6, 0.0).unwrap(),
            ..boundary
        };
        assert!(
            owner
                .isochoric_body_candidate(cx, &uniform_environment(&old, &poor_core), body, solve,)
                .is_err()
        );
        assert_eq!(old, saved);
    });
}

#[path = "conjugate/multiple.rs"]
mod multiple;

#[path = "conjugate/fin.rs"]
mod fin;
