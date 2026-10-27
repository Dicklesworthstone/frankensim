//! Synthetic reaction/collider data exercise native conservative gas and motion.
//! No empirical ignition, fuel-mechanism or engine-fidelity claim is made.
mod troe;

use super::*;
use crate::cavity_phs::piston::reactive::implicit::ImplicitReactionSolve;
use crate::cavity_phs::piston::reactive::material_rate::collider::*;
use crate::cavity_phs::piston::reactive::material_rate::{
    ACTIVITY_RATE_SCALE_DIMS, ACTIVITY_RATE_SCALE_PROPERTY, MaterialRateScale, arrhenius,
};
use fs_material::state_point::{
    ScalarAdmissibility, ScalarPropertyRequirement, resolve_material_state_point,
};
use fs_qty::Dims;
use fs_thermochem::equilibrium::MassActionDirectionalRateStatusV1;

const EFF_A: &str = "collider_efficiency_A";
const EFF_B: &str = "collider_efficiency_B";
const EFF_I: &str = "collider_efficiency_I";
const C0: f64 = 20.0;

fn collider_card(scale: f64, low: f64, efficiencies: [f64; 3]) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, dims, value) in [
        (
            ACTIVITY_RATE_SCALE_PROPERTY,
            ACTIVITY_RATE_SCALE_DIMS,
            scale,
        ),
        (LOW_PRESSURE_SCALE_PROPERTY, ACTIVITY_RATE_SCALE_DIMS, low),
        (
            REFERENCE_CONCENTRATION_PROPERTY,
            arrhenius::MOLAR_CONCENTRATION_DIMS,
            C0,
        ),
        (EFF_A, Dims([0; 6]), efficiencies[0]),
        (EFF_B, Dims([0; 6]), efficiencies[1]),
        (EFF_I, Dims([0; 6]), efficiencies[2]),
        (
            arrhenius::PREFACTOR_PROPERTY,
            ACTIVITY_RATE_SCALE_DIMS,
            2000.0,
        ),
        (arrhenius::TEMPERATURE_EXPONENT_PROPERTY, Dims([0; 6]), 0.0),
        (
            arrhenius::ACTIVATION_ENERGY_PROPERTY,
            arrhenius::MOLAR_ACTIVATION_ENERGY_DIMS,
            5000.0,
        ),
        (
            arrhenius::REFERENCE_TEMPERATURE_PROPERTY,
            Temperature::DIMS,
            600.0,
        ),
    ] {
        let (value, interpolation) = if name == EFF_B {
            (
                PropertyValue::Curve {
                    abscissa: "temperature".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![(200.0, value / 3.0), (1800.0, 3.0 * value)],
                    dims,
                },
                InterpolationPolicy::LinearInside,
            )
        } else {
            (
                PropertyValue::Scalar { value, dims },
                InterpolationPolicy::ConstantWithinValidity,
            )
        };
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value,
                validity: ValidityDomain::unconstrained()
                    .with_quantity(
                        "temperature",
                        quantity(QuantityKind::AbsoluteTemperature),
                        200.0,
                        1800.0,
                    )
                    .with_quantity("pressure", quantity(QuantityKind::Pressure), 1000.0, 1e7),
                uncertainty: UncertaintyModel::Unstated,
                interpolation,
                observations: vec![],
                provenance: provenance(),
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic-reaction-colliders".into(),
            phase: "gas".into(),
            process: "declared-collision-law".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn bindings(regime: ColliderRegime) -> MaterialColliderLaw {
    // Deliberately permuted input, including the nonreactive species.
    MaterialColliderLaw::new(
        regime,
        [("I", EFF_I), ("B", EFF_B), ("A", EFF_A)]
            .into_iter()
            .map(|(id, property)| ColliderEfficiency {
                species: SpeciesId::new(id).unwrap(),
                property: property.into(),
            })
            .collect(),
    )
    .unwrap()
}
fn source(
    scale: f64,
    low: f64,
    efficiencies: [f64; 3],
    regime: ColliderRegime,
    activated: bool,
) -> MaterialRateScale {
    let card = collider_card(scale, low, efficiencies);
    let base = if activated {
        MaterialRateScale::modified_arrhenius(
            card,
            MaterialPropertySelection::SingleClaimOnly,
            Pressure::new(100000.0),
            arrhenius::ArrheniusRateBasis::PressureActivity,
        )
    } else {
        MaterialRateScale::new(
            card,
            MaterialPropertySelection::SingleClaimOnly,
            Pressure::new(100000.0),
        )
    }
    .unwrap();
    base.with_colliders(bindings(regime), MaterialPropertySelection::SingleClaimOnly)
        .unwrap()
}
fn parcel(law: &GasReaction, pressure: f64, fractions: [f64; 3]) -> ClosedGasParcel {
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
                600.0,
            )
            .unwrap()
            .with_quantity("pressure", quantity(QuantityKind::Pressure), pressure)
            .unwrap(),
        0.001,
    )
    .unwrap()
}

#[test]
fn g1_inert_collision_pressure_changes_rates_at_fixed_reactant_partial_pressure() {
    let law = reaction(0.0, 3.5, 1.0)
        .with_material_rate(source(
            100.0,
            1.0,
            [0.0, 0.0, 1.0],
            ColliderRegime::ThirdBody,
            false,
        ))
        .unwrap();
    let mut rates = Vec::new();
    for (pressure, fractions) in [(100000.0, [0.5, 0.0, 0.5]), (200000.0, [0.25, 0.0, 0.75])] {
        let before = parcel(&law, pressure, fractions);
        let step = law
            .parcel_candidate(&before, 1e-6, chemical_controls(), || false)
            .unwrap();
        let expected_q = pressure * fractions[2] / (R * 600.0 * C0);
        let expected_rate = 100.0 * 0.5 * expected_q;
        near(
            step.step().rate.kinetic_multiplier().value(),
            expected_q,
            1e-12,
        );
        near(
            step.step().rate.forward_progress_rate().unwrap().value(),
            expected_rate,
            1e-10,
        );
        near(step.step().extent_mol, expected_rate * 0.001 * 1e-6, 1e-16);
        near(
            step.step().species[2].final_mol,
            step.step().species[2].initial_mol,
            1e-14,
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
        rates.push(step.step().rate.forward_progress_rate().unwrap().value());
    }
    near(rates[1] / rates[0], 3.0, 1e-12);
    // The only change is an efficiency value in the reaction material card.
    let efficient = reaction(0.0, 3.5, 1.0)
        .with_material_rate(source(
            100.0,
            1.0,
            [0.0, 0.0, 5.0],
            ColliderRegime::ThirdBody,
            false,
        ))
        .unwrap();
    let gas = parcel(&efficient, 100000.0, [0.5, 0.0, 0.5]);
    let changed = efficient
        .parcel_candidate(&gas, 1e-6, chemical_controls(), || false)
        .unwrap();
    near(
        changed.step().rate.forward_progress_rate().unwrap().value() / rates[0],
        5.0,
        1e-12,
    );
}

#[test]
fn g1_lindemann_uses_actual_mixture_and_recovers_low_and_high_pressure_limits() {
    for low in [2e-4, 20.0, 2e10] {
        for pressure in [50000.0, 100000.0, 200000.0] {
            let law = reaction(0.0, 3.5, 1.0)
                .with_material_rate(source(
                    200.0,
                    low,
                    [1.0, 2.0, 4.0],
                    ColliderRegime::Lindemann,
                    false,
                ))
                .unwrap();
            let before = parcel(&law, pressure, [0.75, 0.0, 0.25]);
            let step = law
                .parcel_candidate(&before, 1e-6, chemical_controls(), || false)
                .unwrap();
            let q = (0.75 + 4.0 * 0.25) * pressure / (R * 600.0 * C0);
            let reduced = low * q / 200.0;
            let expected = 200.0 * reduced / (1.0 + reduced);
            let measured = step.step().rate.forward_progress_rate().unwrap().value()
                / (0.75 * pressure / 100000.0);
            near(measured, expected, 1e-9);
            assert!(step.step().rate.kinetic_multiplier().value() > 0.0);
            assert!(step.step().rate.kinetic_multiplier().value() <= 1.0);
            if low == 2e-4 {
                near(measured / (low * q), 1.0, 4e-6);
            }
            if low == 2e10 {
                near(measured / 200.0, 1.0, 2e-8);
            }
        }
    }
    // Both directional rates receive the same collision factor, so the
    // original pressure-activity equilibrium ratio is preserved.
    let law = reaction(4000.0, 3.5, 1.0)
        .with_material_rate(source(
            200.0,
            20.0,
            [1.0, 2.0, 4.0],
            ColliderRegime::Lindemann,
            false,
        ))
        .unwrap();
    let before = parcel(&law, 100000.0, [0.25, 0.25, 0.5]);
    let result = law
        .parcel_candidate(&before, 1e-6, chemical_controls(), || false)
        .unwrap();
    let rate = &result.step().rate;
    let equilibrium = rate
        .reverse_rate()
        .equilibrium()
        .equilibrium_constant()
        .unwrap()
        .value();
    near(
        rate.forward_progress_rate().unwrap().value()
            / rate.reverse_progress_rate().unwrap().value(),
        equilibrium,
        1e-10,
    );
}

#[test]
fn g0_absent_colliders_stop_chemistry_without_erasing_real_species() {
    for regime in [ColliderRegime::ThirdBody, ColliderRegime::Lindemann] {
        let law = reaction(4000.0, 3.5, 1.0)
            .with_material_rate(source(100.0, 20.0, [0.0, 0.0, 0.0], regime, false))
            .unwrap();
        let before = parcel(&law, 100000.0, [0.25, 0.25, 0.5]);
        let result = law
            .parcel_candidate(&before, drive().duration_s, chemical_controls(), || false)
            .unwrap();
        assert_eq!(result.gas(), &before);
        assert_eq!(result.step().extent_mol, 0.0);
        assert_eq!(result.step().stored_energy_change_j, 0.0);
        let rate = &result.step().rate;
        assert_eq!(rate.kinetic_multiplier().value(), 0.0);
        assert_eq!(
            rate.forward_status(),
            MassActionDirectionalRateStatusV1::ZeroKineticMultiplier
        );
        assert_eq!(
            rate.reverse_status(),
            MassActionDirectionalRateStatusV1::ZeroKineticMultiplier
        );
        assert_eq!(rate.forward_activity_product().unwrap().value(), 0.25);
        assert_eq!(rate.reverse_activity_product().unwrap().value(), 0.25);
        assert!(
            result
                .step()
                .species
                .iter()
                .all(|s| s.initial_mol == s.final_mol)
        );
    }
    let law = reaction(4000.0, 3.5, 1.0)
        .with_material_rate(source(
            1e-200,
            20.0,
            [1e-200; 3],
            ColliderRegime::ThirdBody,
            false,
        ))
        .unwrap();
    let before = parcel(&law, 100000.0, [0.25, 0.25, 0.5]);
    let original = before.clone();
    assert!(
        law.parcel_candidate(&before, drive().duration_s, chemical_controls(), || false)
            .is_err()
    );
    assert_eq!(
        before, original,
        "positive unrepresentable kinetics cannot be published as extinction"
    );
}

#[test]
fn g3_collider_material_swap_changes_conservative_piston_motion_and_hot_rates() {
    let run = |efficiency| {
        let law = reaction(4000.0, 3.5, 1.0)
            .with_material_rate(source(
                100.0,
                2000.0,
                [1.0, 2.0, efficiency],
                ColliderRegime::ThirdBody,
                true,
            ))
            .unwrap();
        let mut piston = MaterialGasPiston::new(
            parcel(&law, 100000.0, [0.75, 0.0, 0.25]),
            body(8000.0),
            0.01,
            0.0,
        )
        .unwrap();
        let initial = piston.clone();
        let energy0 =
            initial.gas().mass_kg() * initial.gas().state().specific_internal_energy_j_kg();
        for _ in 0..32 {
            let temperature = piston.gas().state().state().temperature;
            let amounts = law.inventory(piston.gas()).unwrap();
            let expected_q =
                (amounts[0] + 2.0 * temperature / 600.0 * amounts[1] + efficiency * amounts[2])
                    / (piston.gas().volume_m3() * C0);
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
            near(
                report.reaction.rate.kinetic_multiplier().value(),
                expected_q,
                1e-10,
            );
            let resolved = report.reaction.rate_material.as_ref().unwrap();
            assert_eq!(
                resolved.query().axes().get("temperature").copied(),
                Some(temperature)
            );
            near(
                resolved.property(EFF_B).unwrap().value_si(),
                2.0 * temperature / 600.0,
                1e-12,
            );
            assert_eq!(
                resolved.card_identity(),
                law.rate_material().unwrap().material().content_hash()
            );
            assert!(
                report.total_energy_residual_j.abs() <= mechanical_controls().energy_tolerance_j
            );
            piston = candidate.into_parts().0;
        }
        let energy = piston.gas().mass_kg() * piston.gas().state().specific_internal_energy_j_kg()
            + 0.5 * piston.body().inertia().mass() * piston.piston().velocity_m_s().powi(2)
            + drive().exterior_pressure_pa * (piston.gas().volume_m3() - initial.gas().volume_m3());
        near(
            energy,
            energy0,
            32.0 * mechanical_controls().energy_tolerance_j,
        );
        assert_eq!(piston.body(), initial.body());
        assert_eq!(piston.gas().mass_kg(), initial.gas().mass_kg());
        piston
    };
    let ordinary = run(1.0);
    let efficient = run(5.0);
    assert!(efficient.piston().velocity_m_s() > ordinary.piston().velocity_m_s());
    assert!(ordinary.piston().velocity_m_s() > 0.0);
    assert!(
        efficient.gas().state().state().temperature > ordinary.gas().state().state().temperature
    );
    assert!(
        efficient.gas().state().state().sound_speed > ordinary.gas().state().state().sound_speed
    );
}

#[test]
fn g1_implicit_falloff_samples_actual_endpoint_composition_and_temperature() {
    let law = reaction(4000.0, 3.5, 1.0)
        .with_material_rate(source(
            100.0,
            2000.0,
            [1.0, 2.0, 5.0],
            ColliderRegime::Lindemann,
            true,
        ))
        .unwrap()
        .with_implicit(ImplicitReactionSolve {
            extent_bracket_mol: [-1.0, 1.0],
            extent_tolerance_mol: 1e-12,
            maximum_iterations: 128,
        })
        .unwrap();
    let before = parcel(&law, 100000.0, [0.75, 0.0, 0.25]);
    let result = law
        .parcel_candidate(&before, drive().duration_s, chemical_controls(), || false)
        .unwrap();
    let after = result.gas();
    let temperature = after.state().state().temperature;
    let n = law.inventory(after).unwrap();
    let q = (n[0] + 2.0 * temperature / 600.0 * n[1] + 5.0 * n[2]) / (after.volume_m3() * C0);
    let high = 2000.0 * (-5000.0 / (R * temperature)).exp();
    let reduced = 2000.0 * q / high;
    near(
        result.step().rate.kinetic_multiplier().value(),
        reduced / (1.0 + reduced),
        1e-10,
    );
    near(
        result
            .step()
            .rate
            .reverse_rate()
            .forward_progress_rate_scale()
            .value(),
        high,
        1e-8,
    );
    near(
        result.step().rate_material.as_ref().unwrap().query().axes()["temperature"],
        temperature,
        1e-7,
    );
    let residual = result.step().extent_mol
        - drive().duration_s
            * before.volume_m3()
            * result.step().rate.net_progress_rate().unwrap().value();
    assert!(residual.abs() <= 1e-12);
    assert!(temperature > 600.0);
}

#[test]
fn g4_catalog_omissions_and_collision_source_switches_refuse_atomically() {
    let missing = MaterialColliderLaw::new(
        ColliderRegime::ThirdBody,
        bindings(ColliderRegime::ThirdBody).efficiencies()[..2].to_vec(),
    )
    .unwrap();
    let incomplete = MaterialRateScale::new(
        collider_card(100.0, 20.0, [1.0; 3]),
        MaterialPropertySelection::SingleClaimOnly,
        Pressure::new(100000.0),
    )
    .unwrap()
    .with_colliders(missing, MaterialPropertySelection::SingleClaimOnly)
    .unwrap();
    assert!(
        reaction(4000.0, 3.5, 1.0)
            .with_material_rate(incomplete)
            .is_err()
    );
    let duplicate = ColliderEfficiency {
        species: SpeciesId::new("A").unwrap(),
        property: EFF_A.into(),
    };
    assert!(
        MaterialColliderLaw::new(
            ColliderRegime::ThirdBody,
            vec![duplicate.clone(), duplicate]
        )
        .is_err()
    );

    let original_card = collider_card(1000.0, 20.0, [1.0, 2.0, 5.0]);
    let mut claims = ClaimSet::new();
    for (_, claim) in original_card.claims().claims_ordered() {
        if claim.key.name() != EFF_I {
            claims.insert_claim(claim.clone()).unwrap();
        }
    }
    let mut cold = original_card.claims_for(EFF_I)[0].1.clone();
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
        value: 10.0,
        dims: Dims([0; 6]),
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
    let card = MaterialCard::assemble(original_card.id().clone(), claims, vec![]).unwrap();
    let base = MaterialRateScale::new(
        card.clone(),
        MaterialPropertySelection::SingleClaimOnly,
        Pressure::new(100000.0),
    )
    .unwrap();
    assert!(
        base.clone()
            .with_colliders(
                bindings(ColliderRegime::ThirdBody),
                MaterialPropertySelection::SingleClaimOnly
            )
            .is_err()
    );
    let pins = [REFERENCE_CONCENTRATION_PROPERTY, EFF_A, EFF_B, EFF_I]
        .into_iter()
        .map(|name| {
            (
                name.to_owned(),
                if name == EFF_I {
                    cold_id
                } else {
                    card.claims_for(name)[0].0
                },
            )
        })
        .collect();
    let law = reaction(4000.0, 3.5, 1.0)
        .with_material_rate(
            base.with_colliders(
                bindings(ColliderRegime::ThirdBody),
                MaterialPropertySelection::PinnedByProperty(pins),
            )
            .unwrap(),
        )
        .unwrap();
    let piston = MaterialGasPiston::new(
        parcel(&law, 100000.0, [0.75, 0.0, 0.25]),
        body(8000.0),
        0.01,
        0.0,
    )
    .unwrap();
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
    let available = resolve_material_state_point(
        &card,
        heated.gas().point(),
        &[ScalarPropertyRequirement::try_new(
            EFF_I,
            Dims([0; 6]),
            ScalarAdmissibility::NonNegative,
        )
        .unwrap()],
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    assert_eq!(available.property(EFF_I).unwrap().value_si(), 10.0);
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
                polls > 8
            }
        ),
        Err(GasReactionError::Cancelled)
    ));
    assert_eq!(piston, original);
}
