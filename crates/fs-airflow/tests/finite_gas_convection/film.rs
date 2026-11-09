//! Actual finite-gas heat/motion uses source-pinned film properties.
use super::*;
use fs_airflow::graph::thermal::medium::convection::finite_gas::film::{
    FilmGasHeatLaw, FilmGasSurface,
};
use fs_convection::{HorizontalPlateFacing, NaturalConvectionGeometry};
use fs_material::fluid_transport::FluidTransport;
use fs_qty::Temperature;
use fs_qty::semantic::{Composition, CompositionBasis};

#[path = "film/mixed.rs"]
mod mixed;
#[path = "film/mass_flow.rs"]
mod mass_flow;

fn natural_surface() -> NaturalGasSurface {
    NaturalGasSurface {
        geometry: NaturalConvectionGeometry::VerticalPlate,
        area_m2: 0.02,
        characteristic_length: Length::new(0.1),
        gravity_m_s2: 9.80665,
    }
}

#[test]
fn oriented_film_convection_reaches_original_gas_body_and_motion_energy() {
    let curve = source::body_curve();
    let old = piston(0.75, 400.0, 100000.0);
    let saved = old.clone();
    let thermal = body(&curve, 260000.0); // 500 K
    let mut conductances = Vec::new();
    for (geometry, expected) in [
        (
            NaturalConvectionGeometry::HorizontalPlate(HorizontalPlateFacing::Up),
            CorrelationId::HorizontalPlateUnstableLaminar,
        ),
        (
            NaturalConvectionGeometry::HorizontalPlate(HorizontalPlateFacing::Down),
            CorrelationId::HorizontalPlateStable,
        ),
        (
            NaturalConvectionGeometry::HorizontalCylinder,
            CorrelationId::ChurchillChuHorizontalCylinder,
        ),
    ] {
        let patch = NaturalGasSurface {
            geometry,
            ..natural_surface()
        };
        // Horizontal source cards cannot silently fall back to bulk properties.
        assert!(patch.resolve(old.gas().state(), thermal.initial).is_err());
        let film = FilmGasSurface {
            model: old.gas().model(),
            law: FilmGasHeatLaw::Natural(patch),
        };
        let fixed = film.resolve(old.gas().state(), thermal.initial).unwrap();
        let sample = fixed.source.interface.evaluation.as_ref().unwrap();
        assert_eq!(sample.nusselt.card().id, expected);
        conductances.push(fixed.conductance_w_k);
        let next = film
            .finite_body_candidate(&old, input(), thermal, controls(), exchange(), || false)
            .unwrap();
        let report = next.step();
        let source = &next.thermal().report;
        let property = source.property_gas.as_ref().unwrap();
        let gas = &source.interface.gas;
        near(
            property.state().temperature,
            0.5 * (gas.state().temperature + source.interface.body.temperature_k()),
            1e-12,
        );
        assert_eq!(property.state().pressure, gas.state().pressure);
        assert_eq!(
            property.thermodynamics().receipt().declared_fraction_bits(),
            gas.thermodynamics().receipt().declared_fraction_bits()
        );
        let du = old.gas().mass_kg()
            * (next.piston().gas().state().specific_internal_energy_j_kg()
                - old.gas().state().specific_internal_energy_j_kg());
        let dh = thermal.mass_kg
            * (next.thermal().state.specific_enthalpy_j_kg()
                - thermal.initial.specific_enthalpy_j_kg());
        let kinetic = 0.5 * old.mass_kg() * next.piston().velocity_m_s().powi(2);
        let work = input().exterior_pressure_pa
            * (next.piston().gas().volume_m3() - old.gas().volume_m3());
        near(du + dh + kinetic + work, 0.0, 1e-7);
        near(dh + report.exchanged_heat_j, 0.0, 1e-7);
        assert!(report.exchanged_heat_j > 0.0);
    }
    // Same source, temperatures, length and area: stability alone halves h.
    near(conductances[0], 2.0 * conductances[1], 1e-12);
    assert_eq!(old, saved);
}

#[test]
fn source_bank_retains_actual_composition_pressure_claims_and_full_thermal_contrast() {
    let bank = source::model_transport(0.75, true);
    let actual = bank
        .with_composition(
            Composition::new(CompositionBasis::MoleFraction, vec![0.25, 0.75]).unwrap(),
        )
        .unwrap();
    let bulk = actual.evaluate(&source::point(400.0, 200000.0)).unwrap();
    let curve = source::body_curve();
    let wall = body(&curve, 260000.0).initial;
    let patch = natural_surface();
    let result = FilmGasSurface {
        model: &bank,
        law: FilmGasHeatLaw::Natural(patch),
    }
    .resolve(&bulk, wall)
    .unwrap();
    let film = result.source.property_gas.as_ref().unwrap();
    assert_eq!(film.state().temperature, 450.0);
    assert_eq!(film.state().pressure, bulk.state().pressure);
    assert_eq!(
        film.thermodynamics().receipt().declared_fraction_bits(),
        bulk.thermodynamics().receipt().declared_fraction_bits()
    );
    assert_eq!(film.source_packs(), bulk.source_packs());
    for (property, original) in film.components().iter().zip(bulk.components()) {
        assert_eq!(property.card_identity(), original.card_identity());
        for (a, b) in property.properties().iter().zip(original.properties()) {
            assert_eq!(a.answer().receipt.selected, b.answer().receipt.selected);
        }
    }
    assert_eq!(result.source.interface.gas, bulk);
    assert_eq!(result.source.interface.body, wall);
    let transport = film.transport_state().unwrap();
    let rho = transport.density_kg_m3();
    let cp = transport.specific_heat_j_kg_k();
    let mu = transport.dynamic_viscosity_pa_s();
    let k = transport.thermal_conductivity_w_m_k();
    let pr = mu * cp / k;
    let ra =
        patch.gravity_m_s2 / 450.0 * (500.0 - 400.0) * 0.1_f64.powi(3) * rho * rho * cp / (mu * k);
    let nu = (0.825
        + 0.387 * ra.powf(1.0 / 6.0) / (1.0 + (0.492 / pr).powf(9.0 / 16.0)).powf(8.0 / 27.0))
    .powi(2);
    let evaluation = result.source.interface.evaluation.as_ref().unwrap();
    near(evaluation.nusselt.groups()["Ra"], ra, ra * 1e-12);
    near(result.conductance_w_k, patch.area_m2 * k * nu / 0.1, 1e-11);
    assert!(
        (result.conductance_w_k - patch.resolve(&bulk, wall).unwrap().conductance_w_k).abs() > 1e-5
    );

    // A bank's inactive species cannot be resurrected by property evaluation.
    let pure = bank
        .with_composition(Composition::new(CompositionBasis::MoleFraction, vec![1.0, 0.0]).unwrap())
        .unwrap();
    let pure_bulk = pure.evaluate(&source::point(400.0, 100000.0)).unwrap();
    let pure_film = bank
        .evaluate_temperature_from_snapshot(&pure_bulk, Temperature::new(450.0))
        .unwrap();
    assert_eq!(pure_film.components().len(), 1);
    assert_eq!(
        pure_film.thermodynamics().receipt().species(),
        pure_bulk.thermodynamics().receipt().species()
    );
    // Preserve the declared mass-basis convention as well as the physical mole fractions.
    let mass = bank
        .with_composition(
            Composition::new(CompositionBasis::MassFraction, vec![0.25, 0.75]).unwrap(),
        )
        .unwrap()
        .evaluate(&source::point(400.0, 100000.0))
        .unwrap();
    let mass_film = bank
        .evaluate_temperature_from_snapshot(&mass, Temperature::new(450.0))
        .unwrap();
    assert_eq!(
        mass_film.thermodynamics().receipt().declared_basis(),
        CompositionBasis::MassFraction
    );
    assert_eq!(
        mass_film.thermodynamics().receipt().mole_fraction_bits(),
        mass.thermodynamics().receipt().mole_fraction_bits()
    );
}

#[test]
fn film_changes_the_same_finite_gas_body_and_piston_solution_without_extra_heat() {
    let curve = source::body_curve();
    let old = piston(0.75, 400.0, 100000.0);
    let saved = old.clone();
    let thermal = body(&curve, 260000.0);
    let bulk = surface()
        .finite_body_candidate(&old, input(), thermal, controls(), exchange(), || false)
        .unwrap();
    let next = FilmGasSurface {
        model: old.gas().model(),
        law: FilmGasHeatLaw::Forced(surface()),
    }
    .finite_body_candidate(&old, input(), thermal, controls(), exchange(), || false)
    .unwrap();
    let report = next.step();
    let source = &next.thermal().report;
    let film = source.property_gas.as_ref().unwrap();
    let gas = &source.interface.gas;
    assert_eq!(*gas, report.piston.midpoint);
    assert_eq!(source.interface.body, report.thermal_midpoint);
    near(
        film.state().temperature,
        0.5 * gas.state().temperature + 0.5 * source.interface.body.temperature_k(),
        1e-12,
    );
    assert_eq!(film.state().pressure, gas.state().pressure);
    // Constant mu/k/cp fixtures give h_film/h_bulk = sqrt(T_bulk/T_film).
    let same_bulk = surface().resolve(gas, source.interface.body).unwrap();
    near(
        report.conductance_w_k / same_bulk.conductance_w_k,
        (gas.state().temperature / film.state().temperature).sqrt(),
        1e-12,
    );
    near(
        report.exchanged_heat_j,
        DT * report.conductance_w_k
            * (source.interface.body.temperature_k() - gas.state().temperature),
        exchange().source_tolerance_j,
    );
    let du = old.gas().mass_kg()
        * (next.piston().gas().state().specific_internal_energy_j_kg()
            - old.gas().state().specific_internal_energy_j_kg());
    let dh = thermal.mass_kg
        * (next.thermal().state.specific_enthalpy_j_kg()
            - thermal.initial.specific_enthalpy_j_kg());
    let dk = 0.5 * old.mass_kg() * next.piston().velocity_m_s().powi(2);
    let work =
        input().exterior_pressure_pa * (next.piston().gas().volume_m3() - old.gas().volume_m3());
    near(du + dh + dk + work, 0.0, 1e-7);
    near(dh + report.exchanged_heat_j, 0.0, 1e-7);
    assert!(next.step().exchanged_heat_j < bulk.step().exchanged_heat_j);
    assert!(
        next.piston().gas().state().state().temperature
            < bulk.piston().gas().state().state().temperature
    );
    assert_eq!(old, saved);
}

#[test]
fn film_source_support_and_wrong_bank_refuse_and_late_cancellation_is_atomic() {
    let curve = source::body_curve();
    let old = piston(0.75, 400.0, 100000.0);
    let saved = old.clone();
    let thermal = body(&curve, 260000.0);
    let bank = old.gas().model();
    assert!(
        bank.evaluate_temperature_from_snapshot(old.gas().state(), Temperature::new(1100.0))
            .is_err()
    );
    let wrong = source::model_transport(0.75, true);
    let mismatched = FilmGasSurface {
        model: &wrong,
        law: FilmGasHeatLaw::Forced(surface()),
    };
    assert!(
        mismatched
            .finite_body_candidate(&old, input(), thermal, controls(), exchange(), || false)
            .is_err()
    );
    let good = FilmGasSurface {
        model: bank,
        law: FilmGasHeatLaw::Natural(natural_surface()),
    };
    let mut count = 0;
    good.finite_body_candidate(&old, input(), thermal, controls(), exchange(), || {
        count += 1;
        false
    })
    .unwrap();
    let mut calls = 0;
    assert!(matches!(
        good.finite_body_candidate(&old, input(), thermal, controls(), exchange(), || {
            calls += 1;
            calls == count
        }),
        Err(ThermalPistonError::Cancelled)
    ));
    assert_eq!(old, saved);
    let equilibrium = body(&curve, 100000.0).initial;
    let zero = good.resolve(old.gas().state(), equilibrium).unwrap();
    assert_eq!(zero.conductance_w_k, 0.0);
    assert!(zero.source.interface.evaluation.is_none());
}

#[test]
fn sphere_film_conductive_and_buoyant_limits_exchange_signed_heat_with_original_gas_body() {
    let curve = source::body_curve();
    // Pure species A has Pr approximately 0.72, inside the sphere's source
    // domain. The fixture's helium-like mixture need not have Pr >= 0.7.
    let old = piston(1.0, 400.0, 100000.0);
    let saved = old.clone();
    for gravity in [0.0, 9.80665] {
        for initial_h in [80000.0, 120000.0] {
            // 380 K and 420 K
            let mut thermal = body(&curve, initial_h);
            let diameter = (6.0 * thermal.volume_m3 / core::f64::consts::PI).cbrt();
            let patch = NaturalGasSurface {
                geometry: NaturalConvectionGeometry::Sphere,
                area_m2: core::f64::consts::PI * diameter * diameter,
                characteristic_length: Length::new(diameter),
                gravity_m_s2: gravity,
            };
            thermal.surface_area_m2 = patch.area_m2;
            assert!(patch.resolve(old.gas().state(), thermal.initial).is_err());
            let law = FilmGasSurface {
                model: old.gas().model(),
                law: FilmGasHeatLaw::Natural(patch),
            };
            let next = law
                .finite_body_candidate(&old, input(), thermal, controls(), exchange(), || false)
                .unwrap();
            let report = next.step();
            let source = &next.thermal().report;
            let film = source.property_gas.as_ref().unwrap();
            let gas = &source.interface.gas;
            let wall = source.interface.body.temperature_k();
            let evaluation = source.interface.evaluation.as_ref().unwrap();
            assert_eq!(evaluation.nusselt.card().id, CorrelationId::ChurchillSphere);
            assert_eq!(evaluation.transport, film.transport_state().unwrap());
            assert_eq!(*gas, report.piston.midpoint);
            near(
                film.state().temperature,
                0.5 * (wall + gas.state().temperature),
                1e-12,
            );
            assert_eq!(film.state().pressure, gas.state().pressure);
            assert_eq!(film.source_packs(), gas.source_packs());
            assert_eq!(
                film.thermodynamics().receipt().declared_fraction_bits(),
                gas.thermodynamics().receipt().declared_fraction_bits()
            );
            let conductive_g =
                2.0 * evaluation.transport.thermal_conductivity_w_m_k() / diameter * patch.area_m2;
            if gravity == 0.0 {
                assert_eq!(evaluation.nusselt.groups()["Ra"], 0.0);
                assert_eq!(evaluation.nusselt.evidence().value, 2.0);
                near(report.conductance_w_k, conductive_g, 1e-12);
            } else {
                assert!(report.conductance_w_k > conductive_g);
            }
            near(
                report.exchanged_heat_j,
                DT * report.conductance_w_k * (wall - gas.state().temperature),
                exchange().source_tolerance_j,
            );
            assert!(report.exchanged_heat_j * (thermal.initial.temperature_k() - 400.0) > 0.0);
            let du = old.gas().mass_kg()
                * (next.piston().gas().state().specific_internal_energy_j_kg()
                    - old.gas().state().specific_internal_energy_j_kg());
            let dh = thermal.mass_kg
                * (next.thermal().state.specific_enthalpy_j_kg()
                    - thermal.initial.specific_enthalpy_j_kg());
            let kinetic = 0.5 * old.mass_kg() * next.piston().velocity_m_s().powi(2);
            let work = input().exterior_pressure_pa
                * (next.piston().gas().volume_m3() - old.gas().volume_m3());
            near(du + dh + kinetic + work, 0.0, 1e-7);
            near(dh + report.exchanged_heat_j, 0.0, 1e-7);
            assert_eq!(next.piston().gas().mass_kg(), old.gas().mass_kg());

            let equilibrium = law
                .resolve(old.gas().state(), body(&curve, 100000.0).initial)
                .unwrap();
            assert_eq!(equilibrium.conductance_w_k, 0.0);
            assert!(equilibrium.source.interface.evaluation.is_none());
        }
    }
    assert_eq!(old, saved);
    // The zero-gravity exception belongs only to the sphere's conductive
    // limit; it must not invent an intercept for other natural cards.
    for geometry in [
        NaturalConvectionGeometry::VerticalPlate,
        NaturalConvectionGeometry::HorizontalCylinder,
    ] {
        let law = FilmGasSurface {
            model: old.gas().model(),
            law: FilmGasHeatLaw::Natural(NaturalGasSurface {
                geometry,
                gravity_m_s2: 0.0,
                ..natural_surface()
            }),
        };
        assert!(
            law.resolve(old.gas().state(), body(&curve, 260000.0).initial)
                .is_err()
        );
    }
}

fn bulk_duct(
    correlation: CorrelationId,
    reynolds: f64,
    gas: &impl FluidTransport,
) -> ForcedGasSurface {
    let transport = gas.transport_state().unwrap();
    ForcedGasSurface {
        correlation,
        speed_m_s: reynolds * transport.dynamic_viscosity_pa_s()
            / (transport.density_kg_m3() * 0.1),
        settings: ForcedConvectionSettings {
            length_over_hydraulic_diameter: Some(20.0),
            direction: ThermalDirection::CoolingFluid,
            ..ForcedConvectionSettings::default()
        },
        ..surface()
    }
}

#[test]
fn bulk_duct_properties_and_heat_direction_match_the_original_correlations() {
    let bank = source::model_transport(1.0, true);
    let bulk = bank.evaluate(&source::point(400.0, 100000.0)).unwrap();
    let transport = bulk.transport_state().unwrap();
    let pr = transport.prandtl();
    let curve = source::body_curve();
    for (correlation, reynolds) in [
        (CorrelationId::CircularDuctHausen, 1000.0),
        (CorrelationId::DittusBoelter, 20000.0),
        (CorrelationId::Gnielinski, 20000.0),
    ] {
        let patch = bulk_duct(correlation, reynolds, &bulk);
        // A modest +/- 5% wall/bulk contrast distinguishes the references.
        // These manufactured properties do not validate an empirical gas.
        for enthalpy in [80000.0, 120000.0] {
            let wall = body(&curve, enthalpy).initial;
            let result = FilmGasSurface {
                model: &bank,
                law: FilmGasHeatLaw::BulkForced(patch),
            }
            .resolve(&bulk, wall)
            .unwrap();
            let property = result.source.property_gas.as_ref().unwrap();
            let evaluation = result.source.interface.evaluation.as_ref().unwrap();
            assert_eq!(property.state().temperature, 400.0);
            assert_eq!(property.state().pressure, bulk.state().pressure);
            assert_eq!(property.source_packs(), bulk.source_packs());
            assert_eq!(property.thermodynamics().receipt().declared_fraction_bits(),
                bulk.thermodynamics().receipt().declared_fraction_bits());
            assert_eq!(evaluation.transport, transport);
            assert_eq!(result.source.interface.gas, bulk);
            near(evaluation.nusselt.groups()["Re"], reynolds, 1e-9);
            near(evaluation.nusselt.groups()["Pr"], pr, 1e-14);
            let nu = match correlation {
                CorrelationId::CircularDuctHausen => {
                    let gz = reynolds * pr / 20.0;
                    3.66 + 0.0668 * gz / (1.0 + 0.04 * gz.powf(2.0 / 3.0))
                }
                CorrelationId::DittusBoelter => {
                    let exponent = if enthalpy > 100000.0 { 0.4 } else { 0.3 };
                    0.023 * reynolds.powf(0.8) * pr.powf(exponent)
                }
                CorrelationId::Gnielinski => {
                    let friction = (0.79 * reynolds.ln() - 1.64).powi(-2);
                    (friction / 8.0) * (reynolds - 1000.0) * pr
                        / (1.0 + 12.7 * (friction / 8.0).sqrt()
                            * (pr.powf(2.0 / 3.0) - 1.0))
                }
                _ => unreachable!("only the three declared duct laws enter this oracle"),
            };
            let expected_h = nu * transport.thermal_conductivity_w_m_k() / 0.1;
            near(evaluation.nusselt.evidence().value, nu, 1e-10);
            near(evaluation.coefficient.value.value(), expected_h, 1e-10);
            near(result.conductance_w_k, patch.area_m2 * expected_h, 1e-12);

            let film = bank.evaluate_temperature_from_snapshot(
                &bulk, Temperature::new(0.5 * (400.0 + wall.temperature_k())),
            ).unwrap();
            let mut settings = patch.settings;
            if enthalpy > 100000.0 {
                settings.direction = ThermalDirection::HeatingFluid;
            }
            let wrong_reference = evaluate_forced_from_material(
                correlation, &film, patch.speed_m_s, patch.characteristic_length, settings,
            ).unwrap();
            assert!((wrong_reference.coefficient().value.value() - expected_h).abs() > 1e-3);
            assert!(FilmGasSurface {
                model: &bank, law: FilmGasHeatLaw::Forced(patch),
            }.resolve(&bulk, wall).is_err());
        }
    }
}

#[test]
fn bulk_reference_still_checks_source_identity_and_keeps_external_film_properties() {
    let bank = source::model_transport(1.0, true);
    let bulk = bank.evaluate(&source::point(400.0, 100000.0)).unwrap();
    let curve = source::body_curve();
    let patch = bulk_duct(CorrelationId::Gnielinski, 20000.0, &bulk);
    // Equal numerical properties at 400 K do not make different source cards
    // interchangeable, including at exact thermal equilibrium.
    let wrong_bank = source::model_transport(1.0, false);
    assert!(FilmGasSurface {
        model: &wrong_bank, law: FilmGasHeatLaw::BulkForced(patch),
    }.resolve(&bulk, body(&curve, 100000.0).initial).is_err());
    let wall = body(&curve, 120000.0).initial;
    let external = FilmGasSurface {
        model: &bank, law: FilmGasHeatLaw::Forced(surface()),
    }.resolve(&bulk, wall).unwrap();
    assert_eq!(external.source.property_gas.as_ref().unwrap().state().temperature, 410.0);
    assert!(FilmGasSurface {
        model: &bank, law: FilmGasHeatLaw::BulkForced(surface()),
    }.resolve(&bulk, wall).is_err());
}

#[test]
fn bulk_duct_heat_closes_the_original_piston_body_energy_in_both_directions() {
    let curve = source::body_curve();
    let old = GasPiston::new(
        ClosedGasParcel::new(
            source::model_transport(1.0, true), source::point(400.0, 100000.0), 0.001,
        ).unwrap(),
        1.0, 0.01, 0.0,
    ).unwrap();
    let saved = old.clone();
    let patch = bulk_duct(CorrelationId::Gnielinski, 20000.0, old.gas().state());
    let law = FilmGasSurface {
        model: old.gas().model(), law: FilmGasHeatLaw::BulkForced(patch),
    };
    for enthalpy in [80000.0, 120000.0] {
        let thermal = body(&curve, enthalpy);
        let next = law.finite_body_candidate(
            &old, input(), thermal, controls(), exchange(), || false,
        ).unwrap();
        let report = next.step();
        let sample = &next.thermal().report;
        let gas = &sample.interface.gas;
        let property = sample.property_gas.as_ref().unwrap();
        assert_eq!(*gas, report.piston.midpoint);
        assert_eq!(property.state().temperature, gas.state().temperature);
        assert_eq!(sample.interface.evaluation.as_ref().unwrap().transport,
            gas.transport_state().unwrap());
        assert_eq!(report.exchanged_heat_j > 0.0, enthalpy > 100000.0);
        near(report.exchanged_heat_j,
            DT * report.conductance_w_k
                * (sample.interface.body.temperature_k() - gas.state().temperature),
            exchange().source_tolerance_j);
        let du = old.gas().mass_kg()
            * (next.piston().gas().state().specific_internal_energy_j_kg()
                - old.gas().state().specific_internal_energy_j_kg());
        let dh = thermal.mass_kg
            * (next.thermal().state.specific_enthalpy_j_kg()
                - thermal.initial.specific_enthalpy_j_kg());
        let kinetic = 0.5 * old.mass_kg() * next.piston().velocity_m_s().powi(2);
        let work = input().exterior_pressure_pa
            * (next.piston().gas().volume_m3() - old.gas().volume_m3());
        near(du + dh + kinetic + work, 0.0, exchange().energy_tolerance_j);
        near(dh + report.exchanged_heat_j, 0.0, exchange().energy_tolerance_j);
        assert_eq!(next.thermal().external_heat_j, -report.exchanged_heat_j);
    }
    assert_eq!(old, saved);
}

#[test]
fn isothermal_piston_surfaces_reject_constant_heat_flux_duct_laws() {
    let old = piston(1.0, 400.0, 100000.0);
    let saved = old.clone();
    let curve = source::body_curve();
    let thermal = body(&curve, 120000.0);
    for correlation in [
        CorrelationId::CircularDuctLaminarChf,
        CorrelationId::RectangularDuctLaminarChf,
    ] {
        let mut patch = bulk_duct(correlation, 1000.0, old.gas().state());
        patch.settings.aspect_ratio = Some(0.5);
        assert!(patch.resolve(old.gas().state(), thermal.initial).is_err());
        assert!(matches!(patch.finite_body_candidate(
            &old, input(), thermal, controls(), exchange(), || false,
        ), Err(ThermalPistonError::Interface(_))));
        for law in [FilmGasHeatLaw::Forced(patch), FilmGasHeatLaw::BulkForced(patch)] {
            let source_law = FilmGasSurface { model: old.gas().model(), law };
            assert!(source_law.resolve(old.gas().state(), thermal.initial).is_err());
            assert!(matches!(source_law.finite_body_candidate(
                &old, input(), thermal, controls(), exchange(), || false,
            ), Err(ThermalPistonError::Interface(_))));
        }
        let inactive = ForcedGasSurface { area_m2: 0.0, speed_m_s: 0.0, ..patch };
        let insulated = FilmGasSurface {
            model: old.gas().model(), law: FilmGasHeatLaw::BulkForced(inactive),
        }.resolve(old.gas().state(), thermal.initial).unwrap();
        assert_eq!(insulated.conductance_w_k, 0.0);
        assert!(insulated.source.property_gas.is_none());
    }
    assert_eq!(old, saved);
}
