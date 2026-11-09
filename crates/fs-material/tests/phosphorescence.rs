//! G0/G1/G3 source admission and physical quenching using manufactured claims.

#[path = "phosphorescence/annihilation.rs"]
mod annihilation;
#[path = "phosphorescence/pump.rs"]
mod pump;

use fs_evidence::ValidityDomain;
use fs_la::decay::DecayLimits;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::phosphorescence::*;
use fs_material::state_point::{
    MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
    resolve_material_state_point,
};
use fs_qty::{
    Dims, QuantitySpec, Temperature,
    semantic::{QuantityKind, SemanticType, ValueForm},
};

fn claims() -> Vec<PropertyClaim> {
    [
        (
            PHOSPHORESCENCE_ENERGY_PROPERTIES,
            PHOSPHORESCENCE_ENERGY_DIMS,
            [4e-19, 3e-19],
        ),
        (
            PHOSPHORESCENCE_RADIATIVE_RATE_PROPERTIES,
            PHOSPHORESCENCE_RATE_DIMS,
            [8.0, 0.25],
        ),
        (
            PHOSPHORESCENCE_NONRADIATIVE_RATE_PROPERTIES,
            PHOSPHORESCENCE_RATE_DIMS,
            [1.0, 0.125],
        ),
        (
            PHOSPHORESCENCE_CROSSING_RATE_PROPERTIES,
            PHOSPHORESCENCE_RATE_DIMS,
            [4.0, 0.0625],
        ),
    ]
    .into_iter()
    .flat_map(|(names, dims, values)| {
        names
            .into_iter()
            .zip(values)
            .map(move |(name, value)| PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained(),
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                uncertainty: UncertaintyModel::Unstated,
                observations: Vec::new(),
                provenance: Provenance {
                    source: format!("synthetic phosphorescence coefficient: {name}"),
                    license: "internal-test-use".into(),
                    artifact: None,
                },
            })
    })
    .collect()
}

fn assemble(claims: Vec<PropertyClaim>) -> MaterialCard {
    let mut set = ClaimSet::new();
    for claim in claims {
        set.insert_claim(claim).unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic singlet/triplet fixture".into(),
            phase: "solid".into(),
            process: "manufactured kinetic test".into(),
            revision: 0,
        },
        set,
        Vec::new(),
    )
    .unwrap()
}

fn resolve(card: &MaterialCard, point: &QueryPoint) -> PhosphorescenceStatePoint {
    PhosphorescenceStatePoint::resolve(card, point, MaterialPropertySelection::SingleClaimOnly)
        .unwrap()
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 3e-11 * expected.abs().max(1e-300),
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn g0_original_roles_and_complete_unrelated_properties_are_retained() {
    let mut input = claims();
    let mut density = input[0].clone();
    density.key = PropertyKey::new("density", Dims([-3, 1, 0, 0, 0, 0]));
    density.value = PropertyValue::Scalar {
        value: 1200.0,
        dims: density.key.dims(),
    };
    input.push(density);
    let card = assemble(input);
    let initial = resolve(&card, &QueryPoint::new());
    assert_eq!(initial.decay().excitation_energies_j(), [4e-19, 3e-19]);
    assert_eq!(initial.decay().radiative_rates_s(), [8.0, 0.25]);
    assert_eq!(initial.decay().nonradiative_rates_s(), [1.0, 0.125]);
    assert_eq!(initial.decay().crossing_rates_s(), [4.0, 0.0625]);
    let expanded = initial
        .resolved()
        .with_additional_requirements(
            &card,
            &[ScalarPropertyRequirement::try_new(
                "density",
                Dims([-3, 1, 0, 0, 0, 0]),
                ScalarAdmissibility::StrictlyPositive,
            )
            .unwrap()],
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
    let complete = PhosphorescenceStatePoint::try_from_material(expanded.clone()).unwrap();
    assert_eq!(complete.resolved(), &expanded);
    assert_eq!(complete.resolved().card_identity(), card.content_hash());
    assert_eq!(
        complete.resolved().property("density").unwrap().value_si(),
        1200.0
    );
    assert_eq!(complete.decay(), initial.decay());
    assert_ne!(complete.identity(), initial.identity());
}

#[test]
fn g0_temperature_refresh_keeps_original_claims_and_refuses_source_escape() {
    let mut input = claims();
    let reverse = input.last_mut().unwrap();
    reverse.value = PropertyValue::Curve {
        abscissa: "T".into(),
        abscissa_dims: Temperature::DIMS,
        knots: vec![(300.0, 0.0625), (600.0, 0.5625)],
        dims: PHOSPHORESCENCE_RATE_DIMS,
    };
    reverse.validity = ValidityDomain::unconstrained().with("T", 300.0, 600.0);
    reverse.interpolation = InterpolationPolicy::LinearInside;
    let card = assemble(input.clone());
    let point = |temperature| QueryPoint::new().with("T", temperature).unwrap();
    let cold = resolve(&card, &point(300.0));
    let warm = cold.at_query_point(&card, &point(450.0)).unwrap();
    assert_eq!(warm.decay().crossing_rates_s(), [4.0, 0.3125]);
    assert_eq!(cold.decay().crossing_rates_s(), [4.0, 0.0625]);
    assert_eq!(
        cold.decay().excitation_energies_j(),
        warm.decay().excitation_energies_j()
    );
    assert_ne!(cold.identity(), warm.identity());
    for (a, b) in cold
        .resolved()
        .properties()
        .iter()
        .zip(warm.resolved().properties())
    {
        assert_eq!(a.requirement(), b.requirement());
        assert_eq!(a.answer().receipt.selected, b.answer().receipt.selected);
    }
    assert!(cold.at_query_point(&card, &point(601.0)).is_err());
    assert!(cold.at_query_point(&card, &QueryPoint::new()).is_err());
    input[0].provenance.source.push_str(" different source");
    assert!(
        cold.at_query_point(&assemble(input), &point(450.0))
            .is_err()
    );
    assert_eq!(cold, resolve(&card, &point(300.0)));
}

#[test]
fn g1_sourced_oxygen_quenching_changes_lifetime_light_and_heat_together() {
    let oxygen = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AmountConcentration,
        ValueForm::Static,
    ));
    let mut input = claims();
    for claim in &mut input[6..8] {
        claim.value = PropertyValue::Scalar {
            value: 0.0,
            dims: PHOSPHORESCENCE_RATE_DIMS,
        };
    }
    // Manufactured dynamic-quenching law knr = 0.125 + kq [O2], with
    // kq = 0.5 m^3/(mol s); this is not a measured oxygen/material constant.
    // See https://goldbook.iupac.org/terms/view/S06004 for Stern--Volmer kinetics.
    input[5].value = PropertyValue::Curve {
        abscissa: "dissolved_O2".into(),
        abscissa_dims: oxygen.dims(),
        knots: vec![(0.0, 0.125), (4.0, 2.125)],
        dims: PHOSPHORESCENCE_RATE_DIMS,
    };
    input[5].validity =
        ValidityDomain::unconstrained().with_quantity("dissolved_O2", oxygen, 0.0, 4.0);
    input[5].interpolation = InterpolationPolicy::LinearInside;
    let card = assemble(input);
    let point = |c| {
        QueryPoint::new()
            .with_quantity("dissolved_O2", oxygen, c)
            .unwrap()
    };
    let unquenched = resolve(&card, &point(0.0));
    let quenched = unquenched.at_query_point(&card, &point(4.0)).unwrap();
    let h = 2.0;
    let n0 = 1e12;
    let step = |source: &PhosphorescenceStatePoint| {
        source
            .decay()
            .prepare(h, DecayLimits::default(), || false)
            .unwrap()
            .advance([0.0, n0], [0.0; 2])
            .unwrap()
    };
    let a = step(&unquenched);
    let b = step(&quenched);
    for (source, result) in [(&unquenched, &a), (&quenched, &b)] {
        let kr = 0.25;
        let knr = source.decay().nonradiative_rates_s()[1];
        let remaining = n0 * (-(kr + knr) * h).exp();
        close(result.populations[1], remaining);
        close(
            result.radiative_events[1],
            kr / (kr + knr) * (n0 - remaining),
        );
        close(
            result.nonradiative_heat_j[1],
            knr / (kr + knr) * (n0 - remaining) * 3e-19,
        );
        assert!(result.residual_j.abs() < 3e-18);
    }
    assert!(b.populations[1] < a.populations[1]);
    assert!(b.radiative_events[1] < a.radiative_events[1]);
    assert!(b.nonradiative_heat_j[1] > a.nonradiative_heat_j[1]);
    assert!(quenched.at_query_point(&card, &point(4.1)).is_err());
}

#[test]
fn g0_missing_negative_and_reversed_source_roles_refuse_without_defaults() {
    for edit in [0, 1, 2] {
        let mut input = claims();
        match edit {
            0 => {
                input.pop();
            }
            1 => {
                input[3].value = PropertyValue::Scalar {
                    value: -0.25,
                    dims: PHOSPHORESCENCE_RATE_DIMS,
                };
            }
            _ => {
                input[1].value = PropertyValue::Scalar {
                    value: 5e-19,
                    dims: PHOSPHORESCENCE_ENERGY_DIMS,
                };
            }
        }
        assert!(
            PhosphorescenceStatePoint::resolve(
                &assemble(input),
                &QueryPoint::new(),
                MaterialPropertySelection::SingleClaimOnly
            )
            .is_err()
        );
    }
    let mut input = claims();
    for claim in &mut input[2..] {
        claim.value = PropertyValue::Scalar {
            value: 0.0,
            dims: PHOSPHORESCENCE_RATE_DIMS,
        };
    }
    let stable = resolve(&assemble(input), &QueryPoint::new());
    let step = stable
        .decay()
        .prepare(100.0, DecayLimits::default(), || false)
        .unwrap()
        .advance([5.0, 7.0], [0.0; 2])
        .unwrap();
    assert_eq!(step.populations, [5.0, 7.0]);
    assert_eq!(step.radiative_events, [0.0; 2]);
}

#[test]
fn g0_adopting_a_permissive_bundle_rechecks_original_quantity_and_values() {
    for wrong_units in [false, true] {
        let mut input = claims();
        let dims = if wrong_units {
            PHOSPHORESCENCE_ENERGY_DIMS
        } else {
            PHOSPHORESCENCE_RATE_DIMS
        };
        input[3].key = PropertyKey::new(PHOSPHORESCENCE_RADIATIVE_RATE_PROPERTIES[1], dims);
        input[3].value = PropertyValue::Scalar { value: -0.25, dims };
        let requirements: Vec<_> = input
            .iter()
            .map(|claim| {
                ScalarPropertyRequirement::try_with_key(&claim.key, ScalarAdmissibility::Finite)
                    .unwrap()
            })
            .collect();
        let resolved = resolve_material_state_point(
            &assemble(input),
            &QueryPoint::new(),
            &requirements,
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
        assert!(PhosphorescenceStatePoint::try_from_material(resolved).is_err());
    }
}

#[test]
fn g0_history_curves_and_frequency_applicability_cannot_be_frozen_into_rates() {
    for history in [false, true] {
        let mut input = claims();
        let (axis, dims) = if history {
            ("elapsed_time", Dims([0, 0, 1, 0, 0, 0]))
        } else {
            ("omega", Dims([0, 0, -1, 0, 0, 0]))
        };
        let quantity = QuantitySpec::dimensional(dims);
        input[3].validity = ValidityDomain::unconstrained().with_quantity(axis, quantity, 1.0, 2.0);
        if history {
            input[3].value = PropertyValue::Curve {
                abscissa: axis.into(),
                abscissa_dims: dims,
                knots: vec![(1.0, 0.25), (2.0, 0.5)],
                dims: PHOSPHORESCENCE_RATE_DIMS,
            };
            input[3].interpolation = InterpolationPolicy::LinearInside;
        }
        let point = QueryPoint::new()
            .with_quantity(axis, quantity, 1.5)
            .unwrap();
        assert!(matches!(
            PhosphorescenceStatePoint::resolve(
                &assemble(input),
                &point,
                MaterialPropertySelection::SingleClaimOnly
            ),
            Err(PhosphorescenceError::Property { .. })
        ));
    }
}
