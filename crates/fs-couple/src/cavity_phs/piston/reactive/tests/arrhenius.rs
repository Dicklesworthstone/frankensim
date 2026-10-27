//! Synthetic coefficient data exercise actual conservative gas/piston owners.
//! These are integration/constitutive checks, not measured combustion kinetics.
use super::*;
use crate::cavity_phs::piston::reactive::implicit::ImplicitReactionSolve;
use crate::cavity_phs::piston::reactive::material_rate::arrhenius::*;
use crate::cavity_phs::piston::reactive::material_rate::{
    ACTIVITY_RATE_SCALE_DIMS, MaterialRateScale,
};
use fs_qty::Dims;

fn rate_card(a: f64, b: f64, ea: f64, c0: f64, tmax: f64, energy_dims: Dims) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, dims, value) in [
        (PREFACTOR_PROPERTY, ACTIVITY_RATE_SCALE_DIMS, a),
        (TEMPERATURE_EXPONENT_PROPERTY, Dims([0; 6]), b),
        (ACTIVATION_ENERGY_PROPERTY, energy_dims, ea),
        (REFERENCE_TEMPERATURE_PROPERTY, Temperature::DIMS, 600.0),
        (
            REFERENCE_CONCENTRATION_PROPERTY,
            MOLAR_CONCENTRATION_DIMS,
            c0,
        ),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained()
                    .with_quantity(
                        "temperature",
                        quantity(QuantityKind::AbsoluteTemperature),
                        200.0,
                        tmax,
                    )
                    .with_quantity("pressure", quantity(QuantityKind::Pressure), 1000.0, 1e7),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: provenance(),
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic-Arrhenius-isomers".into(),
            phase: "gas".into(),
            process: "declared-normalized-kinetics".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn source(a: f64, b: f64, ea: f64, c0: f64, basis: ArrheniusRateBasis) -> MaterialRateScale {
    MaterialRateScale::modified_arrhenius(
        rate_card(a, b, ea, c0, 1800.0, MOLAR_ACTIVATION_ENERGY_DIMS),
        MaterialPropertySelection::SingleClaimOnly,
        Pressure::new(100000.0),
        basis,
    )
    .unwrap()
}
fn gas_at(
    law: &GasReaction,
    temperature: f64,
    pressure: f64,
    fractions: [f64; 3],
) -> ClosedGasParcel {
    ClosedGasParcel::new(
        law.catalog()
            .with_composition(
                Composition::new(CompositionBasis::MoleFraction, fractions.to_vec()).unwrap(),
            )
            .unwrap(),
        QueryPoint::new()
            .with_quantity(
                "temperature",
                quantity(QuantityKind::AbsoluteTemperature),
                temperature,
            )
            .unwrap()
            .with_quantity("pressure", quantity(QuantityKind::Pressure), pressure)
            .unwrap(),
        0.001,
    )
    .unwrap()
}

#[test]
fn g1_activity_arrhenius_resolves_actual_temperature_with_explicit_normalization() {
    for (a, b, ea) in [(1200.0, 0.75, 8000.0), (50.0, -0.5, -3000.0)] {
        let law = reaction(4000.0, 3.5, 1.0)
            .with_material_rate(source(a, b, ea, 1.0, ArrheniusRateBasis::PressureActivity))
            .unwrap();
        for temperature in [400.0, 600.0, 1000.0] {
            let gas = gas_at(&law, temperature, 200000.0, [0.75, 0.0, 0.25]);
            let original = gas.clone();
            let step = law
                .parcel_candidate(&gas, 1e-6, chemical_controls(), || false)
                .unwrap();
            // Standard library arithmetic is intentionally independent of the
            // deterministic log-domain operation tree used in production.
            let k = a * (temperature / 600.0_f64).powf(b) * (-ea / (R * temperature)).exp();
            near(step.step().extent_mol, 1e-6 * 0.001 * k * 1.5, 1e-15);
            let resolved = step.step().rate_material.as_ref().unwrap();
            assert_eq!(resolved.properties().len(), 4);
            assert_eq!(
                resolved.query().axes().get("temperature").copied(),
                Some(temperature)
            );
            assert_eq!(
                resolved.card_identity(),
                law.rate_material().unwrap().material().content_hash()
            );
            assert_eq!(gas, original);
            assert_eq!(
                step,
                law.parcel_candidate(&gas, 1e-6, chemical_controls(), || false)
                    .unwrap()
            );
        }
    }
}

#[test]
fn g1_concentration_basis_uses_forward_order_and_preserves_reverse_equilibrium() {
    for order in [1, 2] {
        // For 2 A <-> B, B has twice the mass and twice the molar cv, so
        // changing molecule count does not change T for the zero-offset case.
        let catalog = catalog(
            0.0,
            if order == 1 { 3.5 } else { 6.0 },
            0.028 * order as f64,
        );
        let kinetic = kinetics_with_atoms(&catalog, 1.0, order);
        let law = GasReaction::new(catalog, kinetic)
            .unwrap()
            .with_material_rate(source(
                0.1,
                0.0,
                0.0,
                2.0,
                ArrheniusRateBasis::NormalizedMolarConcentration,
            ))
            .unwrap();
        let mut rates = Vec::new();
        for pressure in [100000.0, 200000.0] {
            let gas = gas_at(&law, 600.0, pressure, [1.0, 0.0, 0.0]);
            let step = law
                .parcel_candidate(&gas, 1e-6, chemical_controls(), || false)
                .unwrap();
            let concentration = gas.mass_kg() / 0.028 / gas.volume_m3();
            let expected = 0.1 * (concentration / 2.0).powi(order as i32);
            let rate = step.step().rate.forward_progress_rate().unwrap().value();
            near(rate, expected, 1e-11);
            near(
                step.step().extent_mol,
                1e-6 * gas.volume_m3() * expected,
                1e-16,
            );
            assert_eq!(
                step.step()
                    .rate_material
                    .as_ref()
                    .unwrap()
                    .properties()
                    .len(),
                5
            );
            rates.push(rate);
        }
        near(rates[1] / rates[0], 2.0_f64.powi(order as i32), 1e-12);
        let products = gas_at(&law, 600.0, 100000.0, [0.0, 1.0, 0.0]);
        let reverse = law
            .parcel_candidate(&products, 1e-6, chemical_controls(), || false)
            .unwrap();
        let equilibrium = law
            .kinetics()
            .reverse_rate_closure()
            .equilibrium()
            .evaluate(Temperature::new(600.0))
            .unwrap()
            .equilibrium_constant()
            .unwrap()
            .value();
        let forward_scale = 0.1 * (100000.0 / (R * 600.0 * 2.0)).powi(order as i32);
        near(
            reverse.step().rate.net_progress_rate().unwrap().value(),
            -forward_scale / equilibrium,
            1e-10,
        );
        assert!(reverse.step().extent_mol < 0.0);
    }
}

#[test]
fn g1_equivalent_normalized_concentration_data_give_same_physical_trajectory() {
    // A/c0 is a conventional first-order coefficient. Changing the arbitrary
    // concentration normalization changes A as DATA, not the reaction dynamics.
    let run = |a, c0| {
        let law = reaction(4000.0, 3.5, 1.0)
            .with_material_rate(source(
                a,
                0.25,
                4000.0,
                c0,
                ArrheniusRateBasis::NormalizedMolarConcentration,
            ))
            .unwrap();
        let mut gas = gas(&law, [0.75, 0.0, 0.25]);
        for _ in 0..16 {
            gas = law
                .parcel_candidate(&gas, 1e-4, chemical_controls(), || false)
                .unwrap()
                .into_parts()
                .0;
        }
        gas
    };
    let a = run(20.0, 1.0);
    let b = run(100.0, 5.0);
    near(
        a.state().state().temperature,
        b.state().state().temperature,
        1e-7,
    );
    near(a.state().state().pressure, b.state().state().pressure, 1e-5);
    near(
        a.state().state().sound_speed,
        b.state().state().sound_speed,
        1e-7,
    );
}

#[test]
fn g3_activation_energy_data_changes_material_piston_motion_and_live_source_queries() {
    let run = |activation| {
        let law = reaction(4000.0, 3.5, 1.0)
            .with_material_rate(source(
                2000.0,
                0.5,
                activation,
                1.0,
                ArrheniusRateBasis::PressureActivity,
            ))
            .unwrap();
        let mut piston =
            MaterialGasPiston::new(gas(&law, [1.0, 0.0, 0.0]), body(8000.0), 0.01, 0.0).unwrap();
        let initial = piston.clone();
        let initial_u =
            initial.gas().mass_kg() * initial.gas().state().specific_internal_energy_j_kg();
        let mut observed_temperatures = Vec::new();
        for _ in 0..32 {
            let temperature = piston.gas().state().state().temperature;
            let candidate = piston
                .reactive_candidate(
                    &law,
                    drive(),
                    chemical_controls(),
                    mechanical_controls(),
                    || false,
                )
                .unwrap();
            let report = candidate.step();
            let source = report.reaction.rate_material.as_ref().unwrap();
            assert_eq!(
                source.query().axes().get("temperature").copied(),
                Some(temperature)
            );
            assert_eq!(
                source.card_identity(),
                law.rate_material().unwrap().material().content_hash()
            );
            assert!(
                report.total_energy_residual_j.abs() <= mechanical_controls().energy_tolerance_j
            );
            observed_temperatures.push(temperature);
            piston = candidate.into_parts().0;
        }
        assert!(observed_temperatures.last().unwrap() > &observed_temperatures[0]);
        let energy = piston.gas().mass_kg() * piston.gas().state().specific_internal_energy_j_kg()
            + 0.5 * piston.body().inertia().mass() * piston.piston().velocity_m_s().powi(2)
            + drive().exterior_pressure_pa * (piston.gas().volume_m3() - initial.gas().volume_m3());
        near(
            energy,
            initial_u,
            32.0 * mechanical_controls().energy_tolerance_j,
        );
        assert_eq!(piston.body(), initial.body());
        assert_eq!(piston.gas().mass_kg(), initial.gas().mass_kg());
        piston
    };
    let faster = run(5000.0);
    let slower = run(10000.0);
    assert!(faster.piston().velocity_m_s() > slower.piston().velocity_m_s());
    assert!(slower.piston().velocity_m_s() > 0.0);
    assert!(faster.gas().state().state().temperature > slower.gas().state().state().temperature);
    assert!(faster.gas().state().state().sound_speed > slower.gas().state().state().sound_speed);
}

#[test]
fn g1_implicit_arrhenius_samples_the_reconstructed_chemical_endpoint() {
    let law = reaction(4000.0, 3.5, 1.0)
        .with_material_rate(source(
            2000.0,
            0.5,
            5000.0,
            1.0,
            ArrheniusRateBasis::PressureActivity,
        ))
        .unwrap()
        .with_implicit(ImplicitReactionSolve {
            extent_bracket_mol: [-1.0, 1.0],
            extent_tolerance_mol: 1e-12,
            maximum_iterations: 128,
        })
        .unwrap();
    let gas = gas(&law, [0.75, 0.0, 0.25]);
    let step = law
        .parcel_candidate(&gas, drive().duration_s, chemical_controls(), || false)
        .unwrap();
    let temperature = step.gas().state().state().temperature;
    assert!(temperature > 600.0);
    near(
        step.step()
            .rate_material
            .as_ref()
            .unwrap()
            .query()
            .axes()
            .get("temperature")
            .copied()
            .unwrap(),
        temperature,
        1e-7,
    );
    let expected = 2000.0 * (temperature / 600.0).sqrt() * (-5000.0 / (R * temperature)).exp();
    near(
        step.step()
            .rate
            .reverse_rate()
            .forward_progress_rate_scale()
            .value(),
        expected,
        1e-8,
    );
    let residual = step.step().extent_mol
        - drive().duration_s
            * gas.volume_m3()
            * step.step().rate.net_progress_rate().unwrap().value();
    assert!(residual.abs() <= 1e-12);
}

#[test]
fn g4_bad_units_unrepresentable_rates_and_evolved_source_refusals_preserve_piston() {
    let base = reaction(4000.0, 3.5, 1.0);
    let piston =
        MaterialGasPiston::new(gas(&base, [1.0, 0.0, 0.0]), body(8000.0), 0.01, 0.0).unwrap();
    let original = piston.clone();
    for (ea, dims) in [
        (5000.0, Temperature::DIMS),
        (1e20, MOLAR_ACTIVATION_ENERGY_DIMS),
        (-1e20, MOLAR_ACTIVATION_ENERGY_DIMS),
    ] {
        let source = MaterialRateScale::modified_arrhenius(
            rate_card(2000.0, 0.5, ea, 1.0, 1800.0, dims),
            MaterialPropertySelection::SingleClaimOnly,
            Pressure::new(100000.0),
            ArrheniusRateBasis::PressureActivity,
        )
        .unwrap();
        let law = base.clone().with_material_rate(source).unwrap();
        assert!(
            piston
                .reactive_candidate(
                    &law,
                    drive(),
                    chemical_controls(),
                    mechanical_controls(),
                    || false
                )
                .is_err()
        );
        assert_eq!(piston, original);
    }
    let source = MaterialRateScale::modified_arrhenius(
        rate_card(
            2000.0,
            0.5,
            5000.0,
            1.0,
            600.000001,
            MOLAR_ACTIVATION_ENERGY_DIMS,
        ),
        MaterialPropertySelection::SingleClaimOnly,
        Pressure::new(100000.0),
        ArrheniusRateBasis::PressureActivity,
    )
    .unwrap();
    let law = base.clone().with_material_rate(source).unwrap();
    let heated = piston
        .reactive_candidate(
            &law,
            drive(),
            chemical_controls(),
            mechanical_controls(),
            || false,
        )
        .unwrap()
        .into_parts()
        .0;
    let accepted = heated.clone();
    assert!(heated.gas().state().state().temperature > 600.000001);
    assert!(matches!(
        heated.reactive_candidate(
            &law,
            drive(),
            chemical_controls(),
            mechanical_controls(),
            || false
        ),
        Err(GasReactionError::RateMaterial(_))
    ));
    assert_eq!(heated, accepted);
    let mut polls = 0;
    assert!(matches!(
        piston.reactive_candidate(
            &law,
            drive(),
            chemical_controls(),
            mechanical_controls(),
            || {
                polls += 1;
                polls > 3
            }
        ),
        Err(GasReactionError::Cancelled)
    ));
    assert_eq!(piston, original);
}

#[test]
fn g4_arrhenius_pins_do_not_switch_disjoint_coefficient_sources_during_heating() {
    use fs_material::state_point::{
        ScalarAdmissibility, ScalarPropertyRequirement, resolve_material_state_point,
    };
    let base = rate_card(
        2000.0,
        0.5,
        5000.0,
        1.0,
        1800.0,
        MOLAR_ACTIVATION_ENERGY_DIMS,
    );
    let mut claims = ClaimSet::new();
    for (_, claim) in base.claims().claims_ordered() {
        if claim.key.name() != PREFACTOR_PROPERTY {
            claims.insert_claim(claim.clone()).unwrap();
        }
    }
    let mut cold = base.claims_for(PREFACTOR_PROPERTY)[0].1.clone();
    cold.validity = ValidityDomain::unconstrained()
        .with_quantity(
            "temperature",
            quantity(QuantityKind::AbsoluteTemperature),
            200.0,
            600.000001,
        )
        .with_quantity("pressure", quantity(QuantityKind::Pressure), 1000.0, 1e7);
    let cold_id = claims.insert_claim(cold.clone()).unwrap();
    let mut hot = cold;
    hot.value = PropertyValue::Scalar {
        value: 5000.0,
        dims: ACTIVITY_RATE_SCALE_DIMS,
    };
    hot.validity = ValidityDomain::unconstrained()
        .with_quantity(
            "temperature",
            quantity(QuantityKind::AbsoluteTemperature),
            600.000002,
            1800.0,
        )
        .with_quantity("pressure", quantity(QuantityKind::Pressure), 1000.0, 1e7);
    claims.insert_claim(hot).unwrap();
    let card = MaterialCard::assemble(base.id().clone(), claims, vec![]).unwrap();
    assert!(
        MaterialRateScale::modified_arrhenius(
            card.clone(),
            MaterialPropertySelection::SingleClaimOnly,
            Pressure::new(100000.0),
            ArrheniusRateBasis::PressureActivity
        )
        .is_err()
    );
    assert!(
        MaterialRateScale::modified_arrhenius(
            card.clone(),
            MaterialPropertySelection::PreferObservationBacked,
            Pressure::new(100000.0),
            ArrheniusRateBasis::PressureActivity
        )
        .is_err()
    );
    let pins = [
        PREFACTOR_PROPERTY,
        TEMPERATURE_EXPONENT_PROPERTY,
        ACTIVATION_ENERGY_PROPERTY,
        REFERENCE_TEMPERATURE_PROPERTY,
    ]
    .into_iter()
    .map(|name| {
        (
            name.to_owned(),
            if name == PREFACTOR_PROPERTY {
                cold_id
            } else {
                card.claims_for(name)[0].0
            },
        )
    })
    .collect();
    let source = MaterialRateScale::modified_arrhenius(
        card.clone(),
        MaterialPropertySelection::PinnedByProperty(pins),
        Pressure::new(100000.0),
        ArrheniusRateBasis::PressureActivity,
    )
    .unwrap();
    let law = reaction(4000.0, 3.5, 1.0)
        .with_material_rate(source)
        .unwrap();
    let piston =
        MaterialGasPiston::new(gas(&law, [1.0, 0.0, 0.0]), body(8000.0), 0.01, 0.0).unwrap();
    let original = piston.clone();
    let heated = piston
        .reactive_candidate(
            &law,
            drive(),
            chemical_controls(),
            mechanical_controls(),
            || false,
        )
        .unwrap()
        .into_parts()
        .0;
    let accepted = heated.clone();
    assert!(heated.gas().state().state().temperature > 600.000002);
    // A new temperature-local selection really WOULD find the hot source.
    let hot = resolve_material_state_point(
        &card,
        heated.gas().point(),
        &[ScalarPropertyRequirement::try_new(
            PREFACTOR_PROPERTY,
            ACTIVITY_RATE_SCALE_DIMS,
            ScalarAdmissibility::StrictlyPositive,
        )
        .unwrap()],
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    assert_eq!(hot.property(PREFACTOR_PROPERTY).unwrap().value_si(), 5000.0);
    // The immutable kinetic binding must retain the original cold claim.
    assert!(matches!(
        heated.reactive_candidate(
            &law,
            drive(),
            chemical_controls(),
            mechanical_controls(),
            || false
        ),
        Err(GasReactionError::RateMaterial(_))
    ));
    assert_eq!(heated, accepted);
    // The same rule holds inside a tentative implicit endpoint search.
    let implicit = law
        .with_implicit(ImplicitReactionSolve {
            extent_bracket_mol: [-1.0, 1.0],
            extent_tolerance_mol: 1e-12,
            maximum_iterations: 128,
        })
        .unwrap();
    assert!(matches!(
        piston.reactive_candidate(
            &implicit,
            drive(),
            chemical_controls(),
            mechanical_controls(),
            || false
        ),
        Err(GasReactionError::RateMaterial(_))
    ));
    assert_eq!(piston, original);
}
