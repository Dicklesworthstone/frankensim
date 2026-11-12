//! G1/G3: prescribed bulk mass flux on the same finite gas/body heat owners.
use super::*;
use fs_airflow::graph::thermal::medium::convection::{
    AidingVerticalPlateSettings,
    finite_gas::film::{
        AidingVerticalPlateGasSurface, ForcedPlateGasSurface,
        ForcedRectangularDuctSectionGasSurface, ForcedTubeSectionGasSurface,
    },
};
use fs_couple::finite_gas_thermal::WallConvectionSolve;

fn plate_reference(property: &impl FluidTransport, speed: f64, length: f64) -> f64 {
    let p = property.transport_state().unwrap();
    let re = p.density_kg_m3() * speed * length / p.dynamic_viscosity_pa_s();
    0.664 * re.sqrt() * p.prandtl().cbrt() * p.thermal_conductivity_w_m_k() / length
}

fn circle_area(diameter: f64) -> f64 {
    core::f64::consts::PI * diameter * diameter / 4.0
}

#[test]
fn g1_mass_flux_uses_original_bulk_density_even_when_properties_live_at_the_film() {
    let model = source::model_transport(0.75, true);
    let curve = source::body_curve();
    let wall = body(&curve, 260000.0).initial; // 500 K.
    let authored = FilmGasSurface {
        model: &model,
        law: FilmGasHeatLaw::Forced(ForcedGasSurface {
            speed_m_s: 0.0,
            ..surface()
        }),
    };
    let flow = authored.with_prescribed_mass_flow(0.02, 0.01).unwrap();
    assert_eq!(flow.mass_flow_kg_s(), 0.02);
    assert_eq!(flow.flow_area_m2(), 0.01);
    let mut conductances = Vec::new();
    for (temperature, pressure) in [(300.0, 100000.0), (400.0, 100000.0), (400.0, 200000.0)] {
        let bulk = model
            .evaluate(&source::point(temperature, pressure))
            .unwrap();
        let sample = flow.resolve(&bulk, wall).unwrap();
        let velocity = sample.source.prescribed_mass_flow.unwrap();
        let property = sample.source.property_gas.as_ref().unwrap();
        let evaluated = sample.source.interface.evaluation.as_ref().unwrap();
        assert_eq!(sample.source.interface.gas, bulk);
        assert_eq!(velocity.bulk_density_kg_m3, bulk.state().density);
        near(
            velocity.speed_m_s * velocity.bulk_density_kg_m3 * 0.01,
            0.02,
            2e-16,
        );
        assert_eq!(property.state().temperature, 0.5 * (temperature + 500.0));
        assert_eq!(property.state().pressure, pressure);
        assert_eq!(property.source_packs(), bulk.source_packs());
        assert_eq!(
            property.thermodynamics().receipt().declared_fraction_bits(),
            bulk.thermodynamics().receipt().declared_fraction_bits()
        );
        let expected = plate_reference(property, 0.02 / (bulk.state().density * 0.01), 0.1);
        near(
            evaluated.coefficient.value.value(),
            expected,
            2e-12 * expected,
        );
        let using_film_density =
            plate_reference(property, 0.02 / (property.state().density * 0.01), 0.1);
        assert!((expected / using_film_density - 1.0).abs() > 0.03);
        let same_flux = authored
            .with_prescribed_mass_flow(0.04, 0.02)
            .unwrap()
            .resolve(&bulk, wall)
            .unwrap();
        assert_eq!(same_flux.conductance_w_k, sample.conductance_w_k);
        conductances.push(sample.conductance_w_k);
    }
    // At fixed composition/temperatures the original ideal source doubles rho
    // with pressure. U halves; the same mass flux preserves the forced Re/h.
    near(conductances[1], conductances[2], 2e-12 * conductances[1]);

    for (rate, area) in [
        (0.0, 0.01),
        (-1.0, 0.01),
        (1.0, 0.0),
        (f64::NAN, 0.01),
        (1.0, f64::INFINITY),
        (f64::MAX, f64::MIN_POSITIVE),
        (f64::MIN_POSITIVE, f64::MAX),
    ] {
        assert!(authored.with_prescribed_mass_flow(rate, area).is_err());
    }
    for law in [
        FilmGasHeatLaw::Natural(natural_surface()),
        FilmGasHeatLaw::AidingVerticalPlate(AidingVerticalPlateGasSurface {
            area_m2: 0.02,
            characteristic_length: Length::new(0.1),
            settings: AidingVerticalPlateSettings {
                speed_up_m_s: 1.0,
                gravity_m_s2: 9.80665,
                maximum_relative_density_departure: 0.1,
            },
        }),
    ] {
        assert!(
            FilmGasSurface { model: &model, law }
                .with_prescribed_mass_flow(0.02, 0.01)
                .is_err()
        );
    }
    let bulk = model.evaluate(&source::point(400.0, 100000.0)).unwrap();
    let wrong = source::model_transport(0.75, false);
    assert!(
        FilmGasSurface {
            model: &wrong,
            ..authored
        }
        .with_prescribed_mass_flow(0.02, 0.01)
        .unwrap()
        .resolve(&bulk, wall)
        .is_err()
    );
    let insulated = FilmGasSurface {
        model: &wrong,
        law: FilmGasHeatLaw::Forced(ForcedGasSurface {
            area_m2: 0.0,
            speed_m_s: 0.0,
            ..surface()
        }),
    }
    .with_prescribed_mass_flow(0.02, 0.01)
    .unwrap()
    .resolve(&bulk, wall)
    .unwrap();
    assert_eq!(insulated.conductance_w_k, 0.0);
    assert!(insulated.source.property_gas.is_none());
    assert!(insulated.source.interface.evaluation.is_none());
    let fixed = FilmGasSurface {
        model: &model,
        law: FilmGasHeatLaw::Forced(surface()),
    }
    .resolve(&bulk, wall)
    .unwrap();
    assert!(fixed.source.prescribed_mass_flow.is_none());
}

#[test]
fn g1_g0_stream_selection_preserves_passage_area_inlet_coordinates_and_original_regimes() {
    let model = source::model_transport(1.0, true);
    let bulk = model.evaluate(&source::point(400.0, 100000.0)).unwrap();
    let transport = bulk.transport_state().unwrap();
    let curve = source::body_curve();
    let wall = body(&curve, 100000.0).initial;
    let section = ForcedTubeSectionGasSurface {
        area_m2: core::f64::consts::PI * 0.01 * (0.3 - 0.1),
        diameter: Length::new(0.01),
        start: Length::new(0.1),
        end: Length::new(0.3),
        speed_m_s: 0.0,
    };
    let laws = [
        (
            FilmGasHeatLaw::BulkForced(ForcedGasSurface {
                area_m2: 0.02,
                characteristic_length: Length::new(0.01),
                speed_m_s: 0.0,
                correlation: CorrelationId::CircularDuctHausen,
                settings: ForcedConvectionSettings {
                    length_over_hydraulic_diameter: Some(100.0),
                    ..Default::default()
                },
            }),
            1000.0,
            0.01,
            circle_area(0.01),
        ),
        (
            FilmGasHeatLaw::BulkForcedLaminarTubeSection(section),
            1000.0,
            0.01,
            circle_area(0.01),
        ),
        (
            FilmGasHeatLaw::BulkForcedTubeSection(section),
            20000.0,
            0.01,
            circle_area(0.01),
        ),
        (
            FilmGasHeatLaw::ForcedPlateSection(ForcedPlateGasSurface {
                area_m2: 0.02,
                start: Length::new(0.025),
                end: Length::new(0.1),
                speed_m_s: 0.0,
                correlation: CorrelationId::FlatPlateLaminarAverage,
            }),
            5000.0,
            0.1,
            0.01,
        ),
    ];
    for (law, re, length, area) in laws {
        let rate = re * transport.dynamic_viscosity_pa_s() * area / length;
        let surface = FilmGasSurface { model: &model, law };
        let sample = surface
            .with_prescribed_mass_flow(rate, area)
            .unwrap()
            .resolve(&bulk, wall)
            .unwrap();
        let evaluated = sample.source.interface.evaluation.as_ref().unwrap();
        near(evaluated.nusselt.groups()["Re"], re, 2e-12 * re);
        near(
            sample.source.prescribed_mass_flow.unwrap().speed_m_s,
            rate / (transport.density_kg_m3() * area),
            1e-12,
        );
        if let Some(actual) = &evaluated.tube_section {
            assert_eq!(actual.start, section.start);
            assert_eq!(actual.end, section.end);
            assert_eq!(actual.diameter, section.diameter);
            assert!(actual.upstream.is_some());
        }
        if let Some(actual) = &evaluated.section {
            assert_eq!(actual.start.value(), 0.025);
            assert_eq!(actual.end.value(), 0.1);
        }
        if !matches!(law, FilmGasHeatLaw::ForcedPlateSection(_)) {
            assert!(surface.with_prescribed_mass_flow(rate, area * 0.5).is_err());
        }
    }
    let rate = 1000.0 * transport.dynamic_viscosity_pa_s() * circle_area(0.01) / 0.01;
    assert!(
        FilmGasSurface {
            model: &model,
            law: FilmGasHeatLaw::BulkForcedLaminarTubeSection(section)
        }
        .with_prescribed_mass_flow(rate * 20.0, circle_area(0.01))
        .unwrap()
        .resolve(&bulk, wall)
        .is_err()
    );

    let model = source::model_pr072();
    let bulk = model
        .evaluate(&source::point(
            400.0,
            fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * 400.0,
        ))
        .unwrap();
    let duct = ForcedRectangularDuctSectionGasSurface {
        area_m2: 4.5 * 0.01 * 0.08,
        hydraulic_diameter: Length::new(0.01),
        aspect_ratio: 0.5,
        start: Length::new(0.04),
        end: Length::new(0.12),
        speed_m_s: 0.0,
    };
    let area = 0.01_f64.powi(2) * 1.5 * 1.5 / 2.0;
    let rate = 100.0 * bulk.transport_state().unwrap().dynamic_viscosity_pa_s() * area / 0.01;
    let surface = FilmGasSurface {
        model: &model,
        law: FilmGasHeatLaw::BulkForcedRectangularDuctSection(duct),
    };
    let sample = surface
        .with_prescribed_mass_flow(rate, area)
        .unwrap()
        .resolve(&bulk, wall)
        .unwrap();
    let evaluated = sample.source.interface.evaluation.as_ref().unwrap();
    near(evaluated.nusselt.groups()["Re"], 100.0, 1e-12);
    let actual = evaluated.rectangular_section.as_ref().unwrap();
    assert_eq!(actual.start, duct.start);
    assert_eq!(actual.end, duct.end);
    assert_eq!(actual.aspect_ratio, 0.5);
    assert!(
        surface
            .with_prescribed_mass_flow(rate, circle_area(0.01))
            .is_err()
    );
}

#[test]
fn g3_prescribed_stream_rebinds_piston_trials_and_keeps_closed_mass_and_energy() {
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
    let curve = source::body_curve();
    let thermal = UniformEnthalpyStepInput {
        // The free piston's gas spring, sqrt(rho c^2 A^2 / (m V)) = 118.3 rad/s
        // (1 kg, 0.01 m2, 1 L of ideal gas at 1 bar, gamma 1.4), is admitted by fs-phs only for
        // rate*dt <= 0.04 (MIDPOINT_ONE_CENT_ANGULAR_STEP): dt <= 3.38e-4 s. 0.3 ms gives 0.0355.
        duration_s: 0.0003,
        ..body(&curve, 260000.0)
    };
    let forcing = PistonInput {
        duration_s: thermal.duration_s,
        ..input()
    };
    let exchange = HeatExchangeSolve {
        heat_bracket_j: [-1.0, 1.0],
        ..exchange()
    };
    let surface = FilmGasSurface {
        model: old.gas().model(),
        law: FilmGasHeatLaw::Forced(surface()),
    };
    let stream = surface.with_prescribed_mass_flow(0.02, 0.01).unwrap();
    let next = stream
        .finite_body_candidate(&old, forcing, thermal, controls(), exchange, || false)
        .unwrap();
    let report = next.step();
    let sample = &next.thermal().report;
    let flow = sample.prescribed_mass_flow.unwrap();
    assert_eq!(sample.interface.gas, report.piston.midpoint);
    assert_eq!(sample.interface.body, report.thermal_midpoint);
    assert_eq!(
        flow.bulk_density_kg_m3,
        report.piston.midpoint.state().density
    );
    near(flow.speed_m_s * flow.bulk_density_kg_m3 * 0.01, 0.02, 2e-16);
    assert!((flow.bulk_density_kg_m3 / old.gas().state().state().density - 1.0).abs() > 1e-10);
    let property = sample.property_gas.as_ref().unwrap();
    let h = plate_reference(property, flow.speed_m_s, 0.1);
    near(
        report.conductance_w_k,
        sample.interface.area_m2 * h,
        2e-12 * h,
    );
    near(
        report.exchanged_heat_j,
        thermal.duration_s
            * report.conductance_w_k
            * (report.thermal_midpoint.temperature_k()
                - report.piston.midpoint.state().temperature),
        exchange.source_tolerance_j,
    );
    let gas_heat = old.gas().mass_kg()
        * (next.piston().gas().state().specific_internal_energy_j_kg()
            - old.gas().state().specific_internal_energy_j_kg());
    let body_heat = thermal.mass_kg
        * (next.thermal().state.specific_enthalpy_j_kg()
            - thermal.initial.specific_enthalpy_j_kg());
    let kinetic = 0.5 * old.mass_kg() * next.piston().velocity_m_s().powi(2);
    let exterior =
        forcing.exterior_pressure_pa * (next.piston().gas().volume_m3() - old.gas().volume_m3());
    near(gas_heat + body_heat + kinetic + exterior, 0.0, 1e-7);
    near(body_heat + report.exchanged_heat_j, 0.0, 1e-7);
    assert_eq!(next.piston().gas().mass_kg(), old.gas().mass_kg());
    assert!(matches!(
        stream.finite_body_candidate(&old, forcing, thermal, controls(), exchange, || true,),
        Err(ThermalPistonError::Cancelled)
    ));
    assert_eq!(old, saved);
}

#[test]
fn g3_fixed_volume_stream_keeps_bulk_density_while_wall_and_gas_temperatures_change() {
    scope(|cx, _| {
        let old = IsochoricGasParcel::new(
            source::model_transport(1.0, true),
            source::point(400.0, 100000.0),
            0.003,
        )
        .unwrap();
        let saved = old.clone();
        let curve = source::body_curve();
        let thermal = UniformEnthalpyStepInput {
            duration_s: 0.2,
            internal_heat_j: 0.01,
            ..body(&curve, 120000.0)
        };
        let boundary = uniform_boundary();
        let mut environment = uniform_environment(&old, &boundary);
        environment.gas_heat_j = 0.02;
        environment.exchange.heat_bracket_j = [-2.0, 2.0];
        let surface = FilmGasSurface {
            model: old.model(),
            law: FilmGasHeatLaw::Forced(surface()),
        };
        let next = surface
            .with_prescribed_mass_flow(0.02, 0.01)
            .unwrap()
            .isochoric_body_candidate(
                cx,
                &environment,
                thermal,
                WallConvectionSolve {
                    temperature_increment_bracket_k: [-5.0, 5.0],
                    temperature_tolerance_k: 1e-8,
                    maximum_iterations: 128,
                },
            )
            .unwrap();
        let report = &next.report;
        let bulk = report.next_gas.state();
        let sample = &report.convection.source.source;
        let flow = sample.prescribed_mass_flow.unwrap();
        assert_eq!(&sample.interface.gas, bulk);
        assert_eq!(flow.bulk_density_kg_m3, bulk.state().density);
        near(
            flow.speed_m_s,
            0.02 / (old.state().state().density * 0.01),
            2e-12,
        );
        assert_ne!(bulk.state().temperature, old.state().state().temperature);
        let property = sample.property_gas.as_ref().unwrap();
        let h = plate_reference(property, flow.speed_m_s, 0.1);
        near(report.convection.coefficient_w_per_m2_k, h, 2e-12 * h);
        near(
            report.solid_to_gas_heat_j,
            thermal.duration_s
                * thermal.surface_area_m2
                * h
                * (next.state.temperature_k() - bulk.state().temperature),
            environment.exchange.source_tolerance_j,
        );
        let gas_heat = old.mass_kg()
            * (bulk.specific_internal_energy_j_kg() - old.state().specific_internal_energy_j_kg());
        let body_heat = thermal.mass_kg
            * (next.state.specific_enthalpy_j_kg() - thermal.initial.specific_enthalpy_j_kg());
        near(
            gas_heat,
            report.solid_to_gas_heat_j + environment.gas_heat_j,
            1e-7,
        );
        near(
            body_heat,
            thermal.internal_heat_j + report.radiation_heat_j - report.solid_to_gas_heat_j,
            1e-7,
        );
        near(report.combined_energy_residual_j, 0.0, 1e-7);
        assert!(report.entropy.convection_production_j_k >= 0.0);
        near(report.entropy.balance_residual_j_k, 0.0, 1e-12);
        assert_eq!(report.next_gas.mass_kg(), old.mass_kg());
        assert_eq!(old, saved);
    });
}
