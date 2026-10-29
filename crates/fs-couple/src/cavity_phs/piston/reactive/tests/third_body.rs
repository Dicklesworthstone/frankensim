//! Synthetic data exercise the native chemical, caloric and mechanical owners.
//! These are constitutive/integration checks, not measured combustion kinetics.
use super::*;
use crate::cavity_phs::piston::reactive::material_rate::{
    ACTIVITY_RATE_SCALE_DIMS, MaterialRateScale,
    arrhenius::{self, ArrheniusRateBasis},
    collider::{
        ColliderEfficiency, ColliderRegime, LOW_PRESSURE_SCALE_PROPERTY, MaterialColliderLaw,
        REFERENCE_CONCENTRATION_PROPERTY,
    },
};
use fs_qty::Dims;
use fs_thermochem::equilibrium::MassActionDirectionalRateStatusV1;

fn card(a: f64, c_ref: f64, efficiencies: [f64; 3], tmax: f64) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, dims, value) in [
        (arrhenius::PREFACTOR_PROPERTY, ACTIVITY_RATE_SCALE_DIMS, a),
        (LOW_PRESSURE_SCALE_PROPERTY, ACTIVITY_RATE_SCALE_DIMS, 10.0),
        (arrhenius::TEMPERATURE_EXPONENT_PROPERTY, Dims([0; 6]), 0.5),
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
        (
            arrhenius::REFERENCE_CONCENTRATION_PROPERTY,
            arrhenius::MOLAR_CONCENTRATION_DIMS,
            2.0,
        ),
        (
            REFERENCE_CONCENTRATION_PROPERTY,
            arrhenius::MOLAR_CONCENTRATION_DIMS,
            c_ref,
        ),
        ("collision_A", Dims([0; 6]), efficiencies[0]),
        ("collision_B", Dims([0; 6]), efficiencies[1]),
        ("collision_I", Dims([0; 6]), efficiencies[2]),
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
            chemistry: "synthetic-collider-isomers".into(),
            phase: "gas".into(),
            process: "explicit-normalized-third-body".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn colliders() -> Vec<ColliderEfficiency> {
    // Deliberately shuffled: source bindings are canonical by original species.
    ["I", "A", "B"]
        .into_iter()
        .map(|species| ColliderEfficiency {
            species: SpeciesId::new(species).unwrap(),
            property: format!("collision_{species}"),
        })
        .collect()
}
fn rate_source(
    a: f64,
    c_ref: f64,
    efficiencies: [f64; 3],
    basis: ArrheniusRateBasis,
    regime: ColliderRegime,
) -> MaterialRateScale {
    MaterialRateScale::modified_arrhenius(
        card(a, c_ref, efficiencies, 1800.0),
        MaterialPropertySelection::SingleClaimOnly,
        Pressure::new(100000.0),
        basis,
    )
    .unwrap()
    .with_colliders(
        MaterialColliderLaw::new(regime, colliders()).unwrap(),
        MaterialPropertySelection::SingleClaimOnly,
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
fn collider_pressure_composition_and_normalization_follow_independent_rate_equations() {
    for basis in [
        ArrheniusRateBasis::PressureActivity,
        ArrheniusRateBasis::NormalizedMolarConcentration,
    ] {
        for (a, c_ref) in [(50.0, 1.0), (250.0, 5.0)] {
            let law = reaction(0.0, 3.5, 1.0)
                .with_material_rate(rate_source(
                    a,
                    c_ref,
                    [1.0, 2.0, 4.0],
                    basis,
                    ColliderRegime::ThirdBody,
                ))
                .unwrap();
            for (temperature, pressure, fractions) in [
                (600.0, 100000.0, [0.75, 0.0, 0.25]),
                (600.0, 200000.0, [0.75, 0.0, 0.25]),
                (900.0, 200000.0, [0.25, 0.25, 0.5]),
            ] {
                let gas = gas_at(&law, temperature, pressure, fractions);
                let original = gas.clone();
                let candidate = law
                    .parcel_candidate(&gas, 1e-7, chemical_controls(), || false)
                    .unwrap();
                let rate = &candidate.step().rate;
                let k = a * (temperature / 600.0_f64).sqrt() * (-5000.0 / (R * temperature)).exp();
                let factor = pressure / (R * temperature * c_ref)
                    * (fractions[0] + 2.0 * fractions[1] + 4.0 * fractions[2]);
                let activity_conversion = match basis {
                    ArrheniusRateBasis::PressureActivity => 1.0,
                    ArrheniusRateBasis::NormalizedMolarConcentration => {
                        100000.0 / (R * temperature * 2.0)
                    }
                };
                let forward = k * activity_conversion * fractions[0] * pressure / 100000.0 * factor;
                // A and B have identical standard-state thermodynamics: Kp=1.
                let reverse = k * activity_conversion * fractions[1] * pressure / 100000.0 * factor;
                near(rate.kinetic_multiplier().value(), factor, 1e-11);
                near(rate.forward_progress_rate().unwrap().value(), forward, 1e-9);
                near(rate.reverse_progress_rate().unwrap().value(), reverse, 1e-9);
                near(
                    candidate.step().extent_mol,
                    1e-7 * gas.volume_m3() * (forward - reverse),
                    1e-15,
                );
                assert_eq!(candidate.step().species[2].stoichiometric_change_mol, 0.0);
                // The endpoint amounts are rebuilt from the new parcel's own
                // composition, so the inert collider keeps its amount to rounding.
                let collider = &candidate.step().species[2];
                assert!(
                    (collider.final_mol - collider.initial_mol).abs()
                        <= 4.0 * f64::EPSILON * collider.initial_mol,
                    "inert collider amount changed: {} -> {}",
                    collider.initial_mol,
                    collider.final_mol
                );
                let resolved = candidate.step().rate_material.as_ref().unwrap();
                assert_eq!(resolved.property("collision_I").unwrap().value_si(), 4.0);
                assert_eq!(resolved.query().axes().get("pressure"), Some(&pressure));
                assert_eq!(gas, original);
            }
        }
    }
}

#[test]
fn absent_effective_colliders_skip_extreme_normalization_but_still_validate_sources() {
    for regime in [ColliderRegime::ThirdBody, ColliderRegime::Lindemann] {
        for c_ref in [1.0, 1e-320] {
            let law = reaction(4000.0, 3.5, 1.0)
                .with_material_rate(rate_source(
                    50.0,
                    c_ref,
                    [0.0, 0.0, 4.0],
                    ArrheniusRateBasis::PressureActivity,
                    regime,
                ))
                .unwrap();
            let before = gas(&law, [0.75, 0.25, 0.0]);
            let original = before.clone();
            let mut polls = 0;
            let result = law
                .parcel_candidate(&before, 0.01, chemical_controls(), || {
                    polls += 1;
                    false
                })
                .unwrap();
            assert_eq!(result.gas(), &before);
            assert_eq!(result.step().extent_mol, 0.0);
            assert_eq!(result.step().stored_energy_change_j, 0.0);
            let rate = &result.step().rate;
            assert_eq!(
                rate.forward_status(),
                MassActionDirectionalRateStatusV1::ZeroKineticMultiplier
            );
            assert_eq!(
                rate.reverse_status(),
                MassActionDirectionalRateStatusV1::ZeroKineticMultiplier
            );
            assert_eq!(rate.kinetic_multiplier().value(), 0.0);
            assert!(rate.forward_activity_product().unwrap().value() > 0.0);
            assert!(rate.reverse_activity_product().unwrap().value() > 0.0);
            assert!(
                result
                    .step()
                    .species
                    .iter()
                    .all(|s| s.initial_mol == s.final_mol)
            );
            // The exact-zero fast path must still honor late cancellation.
            let mut seen = 0;
            assert!(matches!(
                law.parcel_candidate(&before, 0.01, chemical_controls(), || {
                    seen += 1;
                    seen == polls
                }),
                Err(GasReactionError::Cancelled)
            ));
            assert_eq!(before, original);
        }
        // Every claim is resolved even when the sole efficient species is
        // absent. Invalid normalization or that species' efficiency is refused.
        for (c_ref, efficiencies) in [(0.0, [0.0, 0.0, 4.0]), (1.0, [0.0, 0.0, -4.0])] {
            let law = reaction(4000.0, 3.5, 1.0)
                .with_material_rate(rate_source(
                    50.0,
                    c_ref,
                    efficiencies,
                    ArrheniusRateBasis::PressureActivity,
                    regime,
                ))
                .unwrap();
            let before = gas(&law, [0.75, 0.25, 0.0]);
            let original = before.clone();
            assert!(matches!(
                law.parcel_candidate(&before, 0.01, chemical_controls(), || false),
                Err(GasReactionError::RateMaterial(_))
            ));
            assert_eq!(before, original);
        }
    }
}
