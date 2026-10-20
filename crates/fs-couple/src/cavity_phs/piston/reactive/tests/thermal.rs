//! Synthetic source cards exercise chemistry, material heat and momentum together.
use super::*;
use crate::cavity_phs::piston::reactive::thermal::ReactiveThermalPistonError;
use crate::cavity_phs::piston::thermal::layers::MaterialConductionLayers;
use crate::cavity_phs::piston::thermal::{
    ConductanceSample, HeatExchangeSolve, ThermalPistonError,
};
use fs_matdb::SelectionPolicy;
use fs_material::phase::{
    HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase, UniformEnthalpyStepInput,
};
use fs_material::state_point::SPECIFIC_HEAT_CAPACITY_DIMS;
use std::convert::Infallible;

fn point(t: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", quantity(QuantityKind::AbsoluteTemperature), t)
        .unwrap()
        .with("lot", 7.0)
        .unwrap()
}
fn wall(k: f64) -> (MaterialCard, PropertyKey, ResolvedHeatCapacityCurve) {
    let conductivity = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let cp = PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS);
    let density = PropertyKey::new("density", Density::DIMS);
    let mut claims = ClaimSet::new();
    for (key, lo, hi) in [
        (conductivity.clone(), k, 2.0 * k),
        (cp.clone(), 500.0, 500.0),
        (density.clone(), 4000.0, 4000.0),
    ] {
        let dims = key.quantity().dims();
        claims
            .insert_claim(PropertyClaim {
                key,
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![(200.0, lo), (1000.0, hi)],
                    dims,
                },
                validity: ValidityDomain::unconstrained()
                    .with_quantity(
                        "T",
                        quantity(QuantityKind::AbsoluteTemperature),
                        200.0,
                        1000.0,
                    )
                    .with("lot", 7.0, 7.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: vec![],
                provenance: provenance(),
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "thermal-lining".into(),
            phase: "solid".into(),
            process: "synthetic".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let caloric = ResolvedHeatCapacityCurve::try_new(
        &card,
        HeatCapacitySourceQuery {
            heat_capacity: &cp,
            density: &density,
            temperature_axis: "T",
            lower: &point(200.0),
            upper: &point(1000.0),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Solid,
        0.0,
        1e-5,
    )
    .unwrap();
    (card, conductivity, caloric)
}
fn thermal_body(caloric: &ResolvedHeatCapacityCurve) -> UniformEnthalpyStepInput<'_> {
    UniformEnthalpyStepInput {
        curve: caloric.curve(),
        initial: caloric
            .curve()
            .state_at_specific_enthalpy(150000.0)
            .unwrap(),
        mass_kg: 0.001,
        volume_m3: 0.001 / 4000.0,
        surface_area_m2: 0.02,
        internal_heat_j: 0.0,
        duration_s: drive().duration_s,
    }
}
fn layers<'a>(
    card: &'a MaterialCard,
    key: &'a PropertyKey,
    context: &'a QueryPoint,
) -> MaterialConductionLayers<'a> {
    MaterialConductionLayers {
        material: card,
        context,
        temperature_axis: "T",
        liquid_fraction_axis: None,
        conductivity: key,
        selection: MaterialPropertySelection::SingleClaimOnly,
        area_m2: 0.02,
        gas_distance_m: 0.0001,
        body_distance_m: 0.002,
    }
}
fn exchange() -> HeatExchangeSolve {
    HeatExchangeSolve {
        heat_bracket_j: [-1.0, 1.0],
        source_tolerance_j: 1e-10,
        energy_tolerance_j: 1e-7,
        maximum_iterations: 128,
    }
}

#[test]
fn g1_reaction_and_material_wall_close_energy_without_an_extra_chemical_heat_credit() {
    let run = |conductivity, release, dilution, density| {
        let law = reaction(release, 3.5, 1000.0);
        let (card, key, caloric) = wall(conductivity);
        let context = point(500.0);
        let inertial = body(density);
        let mut piston = MaterialGasPiston::new(
            gas(&law, [1.0 - dilution, 0.0, dilution]),
            inertial.clone(),
            0.01,
            0.0,
        )
        .unwrap();
        let initial = piston.clone();
        let mut body = thermal_body(&caloric);
        let h0 = body.initial.specific_enthalpy_j_kg();
        let mut heat_into_body = 0.0;
        for _ in 0..32 {
            let previous_tick = piston.piston().ports.next_tick();
            let expected_reaction = law
                .parcel_candidate(
                    piston.gas(),
                    drive().duration_s,
                    chemical_controls(),
                    || false,
                )
                .unwrap();
            let candidate = piston
                .reactive_thermal_layers_candidate(
                    &law,
                    drive(),
                    body,
                    chemical_controls(),
                    mechanical_controls(),
                    exchange(),
                    layers(&card, &key, &context),
                    || false,
                )
                .unwrap();
            let report = candidate.step();
            // Chemical extent does not depend on how many heat/momentum trials run.
            assert_eq!(&report.reaction, expected_reaction.step());
            assert_eq!(
                candidate.piston().piston().ports.next_tick(),
                previous_tick + 1
            );
            assert_eq!(candidate.piston().body(), &inertial);
            assert_eq!(candidate.piston().gas().mass_kg(), initial.gas().mass_kg());
            assert!(report.thermal.exchange_trials > 1);
            assert!(report.total_energy_residual_j.abs() <= exchange().energy_tolerance_j);
            assert!(report.thermal.source_residual_j.abs() <= exchange().source_tolerance_j);
            assert_eq!(
                candidate.thermal().state.material_card_identity(),
                card.content_hash()
            );
            let source = layers(&card, &key, &context)
                .resolve(
                    &report.thermal.piston.midpoint,
                    report.thermal.thermal_midpoint,
                )
                .unwrap();
            assert_eq!(candidate.thermal().report, source.source);
            assert_eq!(report.thermal.conductance_w_k, source.conductance_w_k);
            assert_eq!(
                candidate
                    .piston()
                    .gas()
                    .model()
                    .thermodynamics()
                    .mole_composition(),
                expected_reaction
                    .gas()
                    .model()
                    .thermodynamics()
                    .mole_composition()
            );
            let (next, thermal, report) = candidate.into_parts();
            heat_into_body -= report.thermal.exchanged_heat_j;
            body.initial = thermal.state;
            piston = next;
        }
        let gas_change = piston.gas().mass_kg()
            * (piston.gas().state().specific_internal_energy_j_kg()
                - initial.gas().state().specific_internal_energy_j_kg());
        let body_change = body.mass_kg * (body.initial.specific_enthalpy_j_kg() - h0);
        let kinetic = 0.5 * piston.body().inertia().mass() * piston.piston().velocity_m_s().powi(2);
        let outside =
            drive().exterior_pressure_pa * (piston.gas().volume_m3() - initial.gas().volume_m3());
        near(
            gas_change + body_change + kinetic + outside,
            0.0,
            32.0 * exchange().energy_tolerance_j,
        );
        near(
            body_change,
            heat_into_body,
            32.0 * exchange().energy_tolerance_j,
        );
        assert!(heat_into_body > 0.0);
        (piston, body.initial.temperature_k(), heat_into_body)
    };
    let (base, base_t, base_heat) = run(0.02, 4000.0, 0.0, 8000.0);
    let (conductive, conductive_t, conductive_heat) = run(2.0, 4000.0, 0.0, 8000.0);
    assert!(conductive_heat > base_heat && conductive_t > base_t);
    assert!(conductive.piston().velocity_m_s() < base.piston().velocity_m_s());
    assert!(conductive.gas().state().state().sound_speed < base.gas().state().state().sound_speed);
    let (energetic, _, energetic_heat) = run(0.02, 8000.0, 0.0, 8000.0);
    assert!(energetic.piston().velocity_m_s() > base.piston().velocity_m_s());
    assert!(energetic_heat > base_heat);
    let (diluted, _, _) = run(0.02, 4000.0, 0.5, 8000.0);
    assert_ne!(
        diluted.gas().state().state().thermal_conductivity,
        base.gas().state().state().thermal_conductivity
    );
    assert_ne!(
        diluted.piston().velocity_m_s(),
        base.piston().velocity_m_s()
    );
    let (light, _, _) = run(0.02, 4000.0, 0.0, 1000.0);
    assert!(light.piston().velocity_m_s() > base.piston().velocity_m_s());
}

#[test]
fn g3_insulation_and_absent_reagents_reduce_to_the_existing_physics_owners() {
    let law = reaction(4000.0, 3.5, 1000.0);
    let (card, key, caloric) = wall(0.0);
    let context = point(500.0);
    let body = thermal_body(&caloric);
    let before = GasPiston::new(gas(&law, [1.0, 0.0, 0.0]), 1.0, 0.01, 0.0).unwrap();
    let candidate = before
        .reactive_thermal_layers_candidate(
            &law,
            drive(),
            body,
            chemical_controls(),
            mechanical_controls(),
            exchange(),
            layers(&card, &key, &context),
            || false,
        )
        .unwrap();
    let expected = before
        .reactive_candidate(
            &law,
            drive(),
            chemical_controls(),
            mechanical_controls(),
            || false,
        )
        .unwrap();
    assert_eq!(candidate.piston(), expected.piston());
    assert_eq!(candidate.piston().ports, expected.piston().ports);
    assert_eq!(&candidate.step().reaction, &expected.step().reaction);
    assert_eq!(&candidate.step().thermal.piston, &expected.step().piston);
    assert_eq!(candidate.thermal().state, body.initial);
    assert_eq!(candidate.step().thermal.exchanged_heat_j, 0.0);

    let (card, key, caloric) = wall(0.2);
    let body = thermal_body(&caloric);
    let inert = GasPiston::new(gas(&law, [0.0, 0.0, 1.0]), 1.0, 0.01, 0.0).unwrap();
    let candidate = inert
        .reactive_thermal_layers_candidate(
            &law,
            drive(),
            body,
            chemical_controls(),
            mechanical_controls(),
            exchange(),
            layers(&card, &key, &context),
            || false,
        )
        .unwrap();
    let expected = inert
        .thermal_candidate(
            drive(),
            body,
            mechanical_controls(),
            exchange(),
            |g, b| layers(&card, &key, &context).resolve(g, b),
            || false,
        )
        .unwrap();
    assert_eq!(candidate.piston(), expected.piston());
    assert_eq!(candidate.piston().ports, expected.piston().ports);
    assert_eq!(candidate.thermal(), expected.thermal());
    assert_eq!(&candidate.step().thermal, expected.step());
    assert_eq!(candidate.step().reaction.extent_mol, 0.0);
}

#[test]
fn g4_source_refusal_wrong_clock_and_late_cancellation_preserve_all_owners() {
    let law = reaction(4000.0, 3.5, 1000.0);
    let (card, key, caloric) = wall(0.2);
    let context = point(500.0);
    let body = thermal_body(&caloric);
    let piston =
        MaterialGasPiston::new(gas(&law, [1.0, 0.0, 0.0]), super::body(8000.0), 0.01, 0.0).unwrap();
    let before = piston.clone();
    let mut polls = 0;
    piston
        .reactive_thermal_layers_candidate(
            &law,
            drive(),
            body,
            chemical_controls(),
            mechanical_controls(),
            exchange(),
            layers(&card, &key, &context),
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    assert!(polls > 4);
    for stop in [1, polls] {
        let mut calls = 0;
        let result = piston.reactive_thermal_layers_candidate(
            &law,
            drive(),
            body,
            chemical_controls(),
            mechanical_controls(),
            exchange(),
            layers(&card, &key, &context),
            || {
                calls += 1;
                calls == stop
            },
        );
        assert!(matches!(result, Err(ReactiveThermalPistonError::Cancelled)));
    }
    let bad_source = piston.reactive_thermal_candidate(
        &law,
        drive(),
        body,
        chemical_controls(),
        mechanical_controls(),
        exchange(),
        |_, _| -> Result<ConductanceSample<()>, &str> { Err("unavailable material interface") },
        || false,
    );
    assert!(matches!(
        bad_source,
        Err(ReactiveThermalPistonError::Thermal(
            ThermalPistonError::Interface(_)
        ))
    ));
    let wrong_clock = UniformEnthalpyStepInput {
        duration_s: body.duration_s * 2.0,
        ..body
    };
    assert!(matches!(
        piston.reactive_thermal_candidate(
            &law,
            drive(),
            wrong_clock,
            chemical_controls(),
            mechanical_controls(),
            exchange(),
            |_, _| -> Result<_, Infallible> {
                Ok(ConductanceSample {
                    conductance_w_k: 0.0,
                    source: (),
                })
            },
            || false
        ),
        Err(ReactiveThermalPistonError::Invalid(_))
    ));
    let other = wall(2.0).0;
    assert!(
        piston
            .reactive_thermal_layers_candidate(
                &law,
                drive(),
                body,
                chemical_controls(),
                mechanical_controls(),
                exchange(),
                layers(&other, &key, &context),
                || false
            )
            .is_err()
    );
    assert_eq!(piston, before);
    assert_eq!(piston.piston().ports, before.piston().ports);
    assert_eq!(
        body.initial,
        caloric
            .curve()
            .state_at_specific_enthalpy(150000.0)
            .unwrap()
    );
}
