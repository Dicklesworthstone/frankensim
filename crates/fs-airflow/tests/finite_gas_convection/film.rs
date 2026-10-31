//! Actual finite-gas heat/motion uses source-pinned film properties.
use super::*;
use fs_airflow::graph::thermal::medium::convection::finite_gas::film::{
    FilmGasHeatLaw, FilmGasSurface,
};
use fs_convection::{HorizontalPlateFacing, NaturalConvectionGeometry};
use fs_material::fluid_transport::FluidTransport;
use fs_qty::Temperature;
use fs_qty::semantic::{Composition, CompositionBasis};

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
