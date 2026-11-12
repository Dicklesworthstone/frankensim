//! G1/G3/G4: original wall-property laws share one conducting/native heat owner.
use super::*;
use fs_airflow::graph::thermal::medium::convection::finite_gas::{
    conjugate::{ConjugateFilmGasSurface, ConjugateGasFilm},
    sphere::ForcedSphereGasSurface,
    tube::ForcedTubeGasSurface,
};
use fs_material::thermochemical_gas::ThermochemicalGasModel;

const TUBE_D: f64 = 0.01;
const TUBE_L: f64 = 1.0;
const SPHERE_D: f64 = 0.1;
const STREAM_AREA: f64 = 0.005;

fn tube(model: &ThermochemicalGasModel, area: f64) -> ForcedTubeGasSurface<'_> {
    ForcedTubeGasSurface {
        model,
        area_m2: area,
        diameter: Length::new(TUBE_D),
        heated_length: Length::new(TUBE_L),
        speed_m_s: 0.0,
    }
}

fn tube_flow_area() -> f64 {
    core::f64::consts::PI * TUBE_D * TUBE_D / 4.0
}

// Independent original equations use mass flux directly, not sampled velocity.
fn tube_reference(
    bulk: &impl FluidTransport,
    surface: &impl FluidTransport,
    rate: f64,
    laminar: bool,
) -> (f64, f64) {
    let p = bulk.transport_state().unwrap();
    let w = surface.transport_state().unwrap();
    let re = rate * TUBE_D / (tube_flow_area() * p.dynamic_viscosity_pa_s());
    let ratio = p.dynamic_viscosity_pa_s() / w.dynamic_viscosity_pa_s();
    let nu = if laminar {
        1.86 * (re * p.prandtl() * TUBE_D / TUBE_L).cbrt() * ratio.powf(0.14)
    } else {
        0.027 * re.powf(0.8) * p.prandtl().cbrt() * ratio.powf(0.14)
    };
    (re, nu * p.thermal_conductivity_w_m_k() / TUBE_D)
}

#[test]
fn both_original_tube_cards_sample_eliminated_surface_and_close_piston_energy() {
    let card = solid(true, 1.0);
    let curve = caloric(&card);
    let wall = layer(&card, &curve, [300.0, 440.0]);
    let geometry = ReservoirLayerGeometry::cylindrical(0.005, 0.0055, TUBE_L).unwrap();
    let area = geometry.surface_area_m2();
    let old = GasPiston::new(
        ClosedGasParcel::new(
            source::model_transport(1.0, true),
            source::point(400.0, 100000.0),
            0.001,
        )
        .unwrap(),
        1.0,
        0.01,
        0.0,
    )
    .unwrap();
    let saved = old.clone();
    // The original source bank's authored mixture must not replace trial pure A.
    let bank = source::model_transport(0.25, true);
    let authored = tube(&bank, area);
    let body = UniformEnthalpyStepInput {
        surface_area_m2: area,
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
    for laminar in [true, false] {
        let initial_re = if laminar { 1600.0 } else { 20000.0 };
        let rate = initial_re
            * old
                .gas()
                .state()
                .transport_state()
                .unwrap()
                .dynamic_viscosity_pa_s()
            * tube_flow_area()
            / TUBE_D;
        let stream = if laminar {
            authored
                .laminar()
                .with_prescribed_mass_flow(rate, tube_flow_area())
        } else {
            authored.with_prescribed_mass_flow(rate, tube_flow_area())
        }
        .unwrap();
        let owner = stream
            .with_conducting_layer(&wall, geometry, config())
            .unwrap();
        let initial = owner
            .resolve(old.gas().state(), body.initial, || Ok(()))
            .unwrap();
        let speed = initial
            .source
            .layer
            .as_ref()
            .unwrap()
            .film
            .prescribed_mass_flow
            .unwrap()
            .speed_m_s;
        let fixed = ForcedTubeGasSurface {
            speed_m_s: speed,
            ..authored
        };
        let fixed = if laminar {
            fixed
                .laminar()
                .with_conducting_layer(&wall, geometry, config())
                .unwrap()
                .resolve(old.gas().state(), body.initial, || Ok(()))
        } else {
            fixed
                .with_conducting_layer(&wall, geometry, config())
                .unwrap()
                .resolve(old.gas().state(), body.initial, || Ok(()))
        }
        .unwrap();
        near(fixed.conductance_w_k, initial.conductance_w_k, 2e-12);
        assert!(
            fixed
                .source
                .layer
                .as_ref()
                .unwrap()
                .film
                .prescribed_mass_flow
                .is_none()
        );

        let next = owner
            .finite_body_candidate(&old, piston_input, body, controls(), heat, || false)
            .unwrap();
        let report = next.step();
        let sample = &next.thermal().report;
        let response = sample.layer.as_ref().unwrap();
        let film = &response.film;
        let bulk = &film.interface.gas;
        let surface = film.surface_gas.as_ref().unwrap();
        let flow = film.prescribed_mass_flow.unwrap();
        let evaluated = film.interface.evaluation.as_ref().unwrap();
        let (re, h) = tube_reference(bulk, surface, rate, laminar);
        assert_eq!(bulk, &report.piston.midpoint);
        assert_eq!(sample.body, report.thermal_midpoint);
        assert_eq!(film.interface.body.value(), response.surface_temperature_k);
        assert_eq!(surface.state().temperature, response.surface_temperature_k);
        assert!((sample.body_temperature_k - response.surface_temperature_k) > 0.05);
        assert_eq!(surface.state().pressure, bulk.state().pressure);
        assert_eq!(surface.source_packs(), bulk.source_packs());
        assert_eq!(surface.components().len(), 1);
        assert_eq!(
            surface.components()[0].card_identity(),
            bulk.components()[0].card_identity()
        );
        assert_eq!(evaluated.transport, bulk.transport_state().unwrap());
        assert_eq!(
            evaluated.nusselt.card().id,
            if laminar {
                CorrelationId::SiederTateLaminar
            } else {
                CorrelationId::SiederTateTurbulent
            }
        );
        near(evaluated.nusselt.groups()["Re"], re, 2e-12 * re);
        near(evaluated.coefficient.value.value(), h, 2e-12 * h);
        assert_eq!(flow.bulk_density_kg_m3, bulk.state().density);
        assert!((surface.state().density / bulk.state().density - 1.0).abs() > 0.02);
        assert!((bulk.state().density / old.gas().state().state().density - 1.0).abs() > 1e-10);
        near(
            flow.speed_m_s * bulk.state().density * tube_flow_area(),
            rate,
            2e-14 * rate,
        );
        let film_power = area * h * (response.surface_temperature_k - bulk.state().temperature);
        let radial_multiplier = core::f64::consts::TAU * TUBE_L / (0.0055_f64 / 0.005).ln();
        let solid_power = radial_multiplier
            * 0.001
            * (sample.body_temperature_k - response.surface_temperature_k)
            * (sample.body_temperature_k + response.surface_temperature_k);
        near(response.power_to_bath_w, film_power, 2e-10);
        near(response.solid_power_to_bath_w, solid_power, 2e-9);
        assert!(report.conductance_w_k < area * h);
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
            * (next.thermal().state.specific_enthalpy_j_kg()
                - body.initial.specific_enthalpy_j_kg());
        let kinetic = 0.5 * old.mass_kg() * next.piston().velocity_m_s().powi(2);
        let outside = piston_input.exterior_pressure_pa
            * (next.piston().gas().volume_m3() - old.gas().volume_m3());
        near(du + dh + kinetic + outside, 0.0, heat.energy_tolerance_j);
        near(dh + report.exchanged_heat_j, 0.0, heat.energy_tolerance_j);
        assert_eq!(next.piston().gas().mass_kg(), old.gas().mass_kg());
        assert!(Arc::ptr_eq(&sample.conductivity_source, wall.source_line()));
        assert_eq!(old, saved);
    }
}

#[test]
fn whitaker_surface_and_radial_material_share_isochoric_energy_and_entropy() {
    scope(|cx, _| {
        let card = solid(true, 2.0);
        let curve = caloric(&card);
        let wall = layer(&card, &curve, [300.0, 440.0]);
        let geometry = ReservoirLayerGeometry::spherical_face(0.05, 0.048).unwrap();
        let area = geometry.surface_area_m2();
        let old = IsochoricGasParcel::new(
            source::model_transport(1.0, true),
            source::point(600.0, 100000.0),
            0.003,
        )
        .unwrap();
        let saved = old.clone();
        let authored = ForcedSphereGasSurface {
            model: old.model(),
            area_m2: area,
            diameter: Length::new(SPHERE_D),
            speed_m_s: 0.0,
        };
        let rate = 1000.0
            * old
                .state()
                .transport_state()
                .unwrap()
                .dynamic_viscosity_pa_s()
            * STREAM_AREA
            / SPHERE_D;
        let owner = authored
            .with_prescribed_mass_flow(rate, STREAM_AREA)
            .unwrap()
            .with_conducting_layer(&wall, geometry, config())
            .unwrap();
        // A hollow sphere with a 10 mm wall. The solid sphere has V/A = 16.7 mm
        // and a convection Biot h(V/A)/k_min = 7.67 * 0.0167 / 0.8 = 0.160,
        // above the lumped body's ceiling of 0.1; this shell has V/A = 8.13 mm
        // and Biot 0.078. The radial layer and film are unchanged.
        let volume = 4.0 * core::f64::consts::PI
            * ((0.5 * SPHERE_D).powi(3) - 0.04_f64.powi(3)) / 3.0;
        let body = UniformEnthalpyStepInput {
            mass_kg: 1000.0 * volume,
            volume_m3: volume,
            surface_area_m2: area,
            duration_s: 0.02,
            ..thermal(&curve, 420.0)
        };
        let initial = owner.resolve(old.state(), body.initial, || Ok(())).unwrap();
        let fixed = ForcedSphereGasSurface {
            speed_m_s: rate / (old.state().state().density * STREAM_AREA),
            ..authored
        }
        .with_conducting_layer(&wall, geometry, config())
        .unwrap()
        .resolve(old.state(), body.initial, || Ok(()))
        .unwrap();
        near(fixed.conductance_w_k, initial.conductance_w_k, 2e-12);
        assert!(
            fixed
                .source
                .layer
                .as_ref()
                .unwrap()
                .film
                .prescribed_mass_flow
                .is_none()
        );
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
        let boundary = LumpedThermalEnvironment {
            transport,
            ..uniform_boundary()
        };
        let mut environment = uniform_environment(&old, &boundary);
        environment.exchange.heat_bracket_j = [-4.0, 4.0];
        let solve = WallConvectionSolve {
            temperature_increment_bracket_k: [-2.0, 2.0],
            temperature_tolerance_k: 1e-8,
            maximum_iterations: 128,
        };
        let next = owner
            .isochoric_body_candidate(cx, &environment, body, solve)
            .unwrap();
        let report = &next.report;
        let feedback = &report.convection.source;
        let sample = &feedback.source;
        let response = sample.layer.as_ref().unwrap();
        let endpoint = feedback.endpoint_source.layer.as_ref().unwrap();
        let film = &response.film;
        let bulk = report.next_gas.state();
        let surface = film.surface_gas.as_ref().unwrap();
        let flow = film.prescribed_mass_flow.unwrap();
        assert_eq!(&film.interface.gas, bulk);
        assert_eq!(sample.body_temperature_k, feedback.temperature_k);
        assert_eq!(
            feedback.endpoint_source.body_temperature_k,
            next.state.temperature_k()
        );
        assert_eq!(surface.state().temperature, response.surface_temperature_k);
        assert!(response.surface_temperature_k > sample.body_temperature_k + 1.0);
        assert!(response.surface_temperature_k < 440.0);
        assert!(bulk.state().temperature > 440.0); // Gas need not lie on the material k(T) chart.
        assert_eq!(surface.state().pressure, bulk.state().pressure);
        assert_eq!(surface.source_packs(), bulk.source_packs());
        assert_eq!(flow.bulk_density_kg_m3, bulk.state().density);
        near(flow.bulk_density_kg_m3, old.state().state().density, 2e-12);
        assert!((surface.state().density / bulk.state().density - 1.0).abs() > 0.3);
        near(
            flow.speed_m_s * bulk.state().density * STREAM_AREA,
            rate,
            2e-14 * rate,
        );
        let p = bulk.transport_state().unwrap();
        let w = surface.transport_state().unwrap();
        let re = rate * SPHERE_D / (STREAM_AREA * p.dynamic_viscosity_pa_s());
        let ratio = p.dynamic_viscosity_pa_s() / w.dynamic_viscosity_pa_s();
        let nu = 2.0
            + (0.4 * re.sqrt() + 0.06 * re.powf(2.0 / 3.0))
                * p.prandtl().powf(0.4)
                * ratio.powf(0.25);
        let h = nu * p.thermal_conductivity_w_m_k() / SPHERE_D;
        let evaluated = film.interface.evaluation.as_ref().unwrap();
        assert_eq!(evaluated.nusselt.card().id, CorrelationId::WhitakerSphere);
        assert_eq!(evaluated.transport, p);
        near(evaluated.nusselt.groups()["Re"], re, 2e-12 * re);
        near(evaluated.coefficient.value.value(), h, 2e-12 * h);
        near(
            response.power_to_bath_w,
            area * h * (response.surface_temperature_k - bulk.state().temperature),
            2e-10,
        );
        let radial_multiplier = 4.0 * core::f64::consts::PI * 0.048 * 0.05 / (0.05 - 0.048);
        let solid_power = radial_multiplier
            * 0.002
            * (sample.body_temperature_k - response.surface_temperature_k)
            * (sample.body_temperature_k + response.surface_temperature_k);
        near(response.solid_power_to_bath_w, solid_power, 2e-9);
        near(
            report.convection.coefficient_w_per_m2_k,
            response.conductance_w_k / area,
            1e-12,
        );
        near(
            report.solid_to_gas_heat_j,
            body.duration_s * endpoint.power_to_bath_w,
            2.0 * environment.exchange.source_tolerance_j,
        );
        let q = report.solid_to_gas_heat_j;
        let du = old.mass_kg()
            * (bulk.specific_internal_energy_j_kg() - old.state().specific_internal_energy_j_kg());
        let dh = body.mass_kg
            * (next.state.specific_enthalpy_j_kg() - body.initial.specific_enthalpy_j_kg());
        near(du, q, 1e-7);
        near(dh, -q, 1e-7);
        near(report.combined_energy_residual_j, 0.0, 1e-7);
        let ts = endpoint.surface_temperature_k;
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
        assert!(q < 0.0);
        assert!(next.state.temperature_k() > body.initial.temperature_k());
        assert!(bulk.state().temperature < 600.0);
        assert!(Arc::ptr_eq(&sample.conductivity_source, wall.source_line()));
        assert_eq!(old, saved);
    });
}

#[test]
fn specialized_admission_insulation_and_legacy_film_type_remain_explicit() {
    let card = solid(true, 1.0);
    let curve = caloric(&card);
    let wall = layer(&card, &curve, [300.0, 440.0]);
    let gas_model = source::model_transport(1.0, true);
    let wrong_model = source::model_transport(1.0, false);
    let bulk = gas_model.evaluate(&source::point(400.0, 100000.0)).unwrap();
    let body = thermal(&curve, 420.0).initial;
    let flat = ReservoirLayerGeometry::planar(AREA, DEPTH).unwrap();
    let fixed_film = FilmGasSurface {
        model: &gas_model,
        law: FilmGasHeatLaw::Forced(ForcedGasSurface {
            speed_m_s: RATE / (bulk.state().density * FLOW_AREA),
            ..surface()
        }),
    };
    // Existing explicitly annotated return types must work for BOTH constructors.
    let fixed: ConjugateFilmGasSurface<'_, '_, '_, _> = fixed_film
        .with_conducting_layer(&wall, flat, config())
        .unwrap();
    let flowed: ConjugateFilmGasSurface<'_, '_, '_, _> = plate(&gas_model)
        .with_prescribed_mass_flow(RATE, FLOW_AREA)
        .unwrap()
        .with_conducting_layer(&wall, flat, config())
        .unwrap();
    let a = fixed.resolve(&bulk, body, || Ok(())).unwrap();
    let b = flowed.resolve(&bulk, body, || Ok(())).unwrap();
    near(a.conductance_w_k, b.conductance_w_k, 2e-12);
    assert!(
        a.source
            .layer
            .as_ref()
            .unwrap()
            .film
            .prescribed_mass_flow
            .is_none()
    );
    assert!(
        b.source
            .layer
            .as_ref()
            .unwrap()
            .film
            .prescribed_mass_flow
            .is_some()
    );

    let radial = ReservoirLayerGeometry::cylindrical(0.005, 0.0055, TUBE_L).unwrap();
    let area = radial.surface_area_m2();
    let rate = 1600.0 * bulk.transport_state().unwrap().dynamic_viscosity_pa_s() * tube_flow_area()
        / TUBE_D;
    let partial = tube(&gas_model, 0.5 * area)
        .laminar()
        .with_prescribed_mass_flow(rate, tube_flow_area())
        .unwrap()
        .with_conducting_layer(&wall, radial.with_area_fraction(0.5).unwrap(), config())
        .unwrap();
    assert!(matches!(
        partial.resolve(&bulk, body, || Ok(())),
        Err(ConjugateGasError::Film(_))
    ));
    let outside_regime = tube(&gas_model, area)
        .laminar()
        .with_prescribed_mass_flow(2.0 * rate, tube_flow_area())
        .unwrap()
        .with_conducting_layer(&wall, radial, config())
        .unwrap();
    assert!(matches!(
        outside_regime.resolve(&bulk, body, || Ok(())),
        Err(ConjugateGasError::Film(_))
    ));
    let wrong = tube(&wrong_model, area)
        .laminar()
        .with_prescribed_mass_flow(rate, tube_flow_area())
        .unwrap()
        .with_conducting_layer(&wall, radial, config())
        .unwrap();
    assert!(matches!(
        wrong.resolve(&bulk, body, || Ok(())),
        Err(ConjugateGasError::Film(_))
    ));
    let zero_tube = tube(&wrong_model, 0.0)
        .laminar()
        .with_prescribed_mass_flow(rate, tube_flow_area())
        .unwrap()
        .with_conducting_layer(&wall, radial, config())
        .unwrap();
    let zero = zero_tube.resolve(&bulk, body, || Ok(())).unwrap();
    assert_eq!(zero.conductance_w_k, 0.0);
    assert_eq!(zero.source.applied_heat_into_gas_w, 0.0);
    assert!(zero.source.layer.is_none());
    let insulated = zero.source.insulated_source.as_ref().unwrap();
    assert!(insulated.surface_gas.is_none());
    assert!(insulated.interface.evaluation.is_none());
    assert!(insulated.prescribed_mass_flow.is_some());

    let spherical = ReservoirLayerGeometry::spherical_face(0.05, 0.048).unwrap();
    let sphere = ForcedSphereGasSurface {
        model: &gas_model,
        area_m2: spherical.surface_area_m2(),
        diameter: Length::new(SPHERE_D),
        speed_m_s: 1.0,
    };
    let heated = sphere
        .with_conducting_layer(&wall, spherical, config())
        .unwrap();
    // Whitaker's original viscosity-ratio range is not clamped for a hotter wall.
    assert!(matches!(
        heated.resolve(&bulk, body, || Ok(())),
        Err(ConjugateGasError::Film(_))
    ));
    let zero_sphere = ForcedSphereGasSurface {
        model: &wrong_model,
        area_m2: 0.0,
        ..sphere
    }
    .with_prescribed_mass_flow(0.001, STREAM_AREA)
    .unwrap()
    .with_conducting_layer(&wall, spherical, config())
    .unwrap();
    let zero = zero_sphere.resolve(&bulk, body, || Ok(())).unwrap();
    assert_eq!(zero.conductance_w_k, 0.0);
    assert!(zero.source.layer.is_none());
    let insulated = zero.source.insulated_source.as_ref().unwrap();
    assert!(insulated.surface_gas.is_none());
    assert!(insulated.interface.evaluation.is_none());
    // Insulation skips unused gas properties but retains the original pure-phase gate.
    let latent = curve.state_at_specific_enthalpy(280000.0).unwrap();
    assert!(matches!(
        zero_tube.resolve(&bulk, latent, || Ok(())),
        Err(ConjugateGasError::Layer(_))
    ));
    assert!(matches!(
        zero_sphere.resolve(&bulk, latent, || Ok(())),
        Err(ConjugateGasError::Layer(_))
    ));
}
