//! Authored kinetic fits check native source admission, not material validation.
use super::*;
use fs_material::field_storage::hysteresis::vector::VectorHysteresisControls;

fn viscosity_keys() -> Vec<PropertyKey> {
    (0..2)
        .map(|i| PropertyKey::new(format!("branch_{i}_viscosity"), MAGNETIC_VISCOSITY_DIMS))
        .collect()
}

fn viscosity_card(
    values: [f64; 2],
    mut change: impl FnMut(&mut PropertyClaim),
) -> Arc<MaterialCard> {
    let original = card(|_| {});
    let mut claims = original.claims().clone();
    for (key, value) in viscosity_keys().into_iter().zip(values) {
        let mut claim = PropertyClaim {
            key,
            value: PropertyValue::Scalar {
                value,
                dims: MAGNETIC_VISCOSITY_DIMS,
            },
            validity: ValidityDomain::unconstrained().with_quantity(
                "T",
                QuantitySpec::dimensional(Temperature::DIMS),
                280.0,
                420.0,
            ),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![],
            provenance: Provenance {
                source: "authored intrinsic switching viscosity fixture".into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        };
        change(&mut claim);
        claims.insert_claim(claim).unwrap();
    }
    Arc::new(MaterialCard::assemble(original.id().clone(), claims, vec![]).unwrap())
}

fn kinetic_source(card: Arc<MaterialCard>) -> ResolvedVectorHysteresis {
    ResolvedVectorHysteresis::from_isotropic_source_with_viscosity(
        resolve(card),
        viscosity_keys(),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}

fn controls() -> VectorHysteresisControls {
    VectorHysteresisControls {
        relative_tolerance: 1e-12,
        absolute_field_tolerance_a_m: 1e-12,
        absolute_flux_tolerance_t: 1e-12,
        work_tolerance_j_m3: 2e-11,
        maximum_branch_iterations: 100,
        maximum_inverse_iterations: 60,
        maximum_backtracks: 40,
    }
}

#[test]
fn original_kinetic_receipts_survive_temperature_refresh_and_retain_native_memory() {
    let card = viscosity_card([0.4, 0.9], |_| {});
    let original = resolve(card.clone());
    let fit = kinetic_source(card);
    assert_eq!(fit.sources(), original.sources());
    assert_eq!(fit.viscosity_sources(), viscosity_keys());
    assert_eq!(fit.law().viscosities_m_s_h(), Some([0.4, 0.9].as_slice()));
    assert_eq!(
        fit.material().properties().len(),
        original.material().properties().len() + 2
    );
    for property in original.material().properties() {
        assert_eq!(
            fit.material()
                .property(property.requirement().name())
                .unwrap(),
            property
        );
    }
    let first = fit
        .law()
        .candidate_field_timed(
            &fit.law().initial_state(),
            [4.0, 0.0, 0.0],
            0.1,
            controls(),
            &mut || false,
        )
        .unwrap();
    let warm = fit.at_query_point(&point(400.0)).unwrap();
    let advance = |source: &ResolvedVectorHysteresis| {
        source
            .law()
            .candidate_field_timed(&first.state, [4.0, 1.0, 0.0], 0.2, controls(), &mut || {
                false
            })
            .unwrap()
    };
    assert_eq!(advance(&fit), advance(&warm));
    assert!(advance(&warm).viscous_heat_j_m3 > 0.0);
    assert_eq!(
        warm.material().card_identity(),
        fit.material().card_identity()
    );
    assert_eq!(warm.material().query(), point(400.0));
    for (before, after) in fit
        .material()
        .properties()
        .iter()
        .zip(warm.material().properties())
    {
        assert_eq!(
            before.answer().receipt.selected,
            after.answer().receipt.selected
        );
    }
    fit.admit_temperature_interval("T", [280.0, 420.0]).unwrap();
    // The reversible/pinning fit alone covers [250,600]. Kinetic validity is
    // stricter and must be visible to every coupled thermal consumer.
    assert!(fit.admit_temperature_interval("T", [250.0, 420.0]).is_err());
    assert!(fit.admit_temperature_interval("T", [280.0, 450.0]).is_err());
    assert!(fit.at_query_point(&point(420.1)).is_err());
    assert_eq!(fit.material().query(), point(300.0));
}

#[test]
fn explicit_zero_profile_recovers_legacy_map_but_keeps_its_source_domain() {
    let card = viscosity_card([0.0, 0.0], |_| {});
    let scalar = resolve(card.clone());
    let original = ResolvedVectorHysteresis::from_isotropic_source(scalar).unwrap();
    let zero = kinetic_source(card);
    assert!(!zero.law().is_rate_dependent());
    assert_eq!(zero.law().viscosities_m_s_h(), None);
    assert_eq!(zero.viscosity_sources(), viscosity_keys());
    let mut state = original.law().initial_state();
    for field in [[4.0, 0.0, 0.0], [0.0, 3.0, 0.0], [-5.0, 0.0, 1.0]] {
        let legacy = original
            .law()
            .candidate_field(&state, field, controls(), &mut || false)
            .unwrap();
        let sourced = zero
            .law()
            .candidate_field(&state, field, controls(), &mut || false)
            .unwrap();
        assert_eq!(legacy, sourced);
        state = legacy.state;
    }
    assert!(zero.at_query_point(&point(450.0)).is_err());
    assert!(
        zero.admit_temperature_interval("T", [250.0, 500.0])
            .is_err()
    );
    let unconditioned = kinetic_source(viscosity_card([0.4, 0.9], |claim| {
        claim.validity = ValidityDomain::unconstrained();
    }));
    assert!(
        unconditioned
            .admit_temperature_interval("T", [280.0, 420.0])
            .is_err()
    );
}

#[test]
fn viscosity_requires_complete_scalar_roles_and_constant_nonnegative_sources() {
    let scalar = resolve(viscosity_card([0.4, 0.9], |_| {}));
    let attempt = |source, keys| {
        ResolvedVectorHysteresis::from_isotropic_source_with_viscosity(
            source,
            keys,
            MaterialPropertySelection::SingleClaimOnly,
        )
    };
    assert!(attempt(scalar.clone(), vec![]).is_err());
    assert!(attempt(scalar.clone(), vec![viscosity_keys()[0].clone()]).is_err());
    let mut wrong = viscosity_keys();
    wrong[0] = PropertyKey::new("branch_0_viscosity", MAGNETIC_PERMEABILITY_DIMS);
    assert!(attempt(scalar, wrong).is_err());
    assert!(
        attempt(
            resolve(viscosity_card([-0.1, 0.9], |_| {})),
            viscosity_keys()
        )
        .is_err()
    );
    let varying = resolve(viscosity_card([0.4, 0.9], |claim| {
        claim.value = PropertyValue::Curve {
            abscissa: "T".into(),
            abscissa_dims: Temperature::DIMS,
            knots: vec![(280.0, 0.4), (420.0, 0.9)],
            dims: MAGNETIC_VISCOSITY_DIMS,
        };
        claim.interpolation = InterpolationPolicy::LinearInside;
    }));
    assert!(attempt(varying, viscosity_keys()).is_err());
    let dynamic = viscosity_card([0.4, 0.9], |claim| {
        claim.validity = claim.validity.clone().with_quantity(
            "omega",
            QuantitySpec::dimensional(Dims([0, 0, -1, 0, 0, 0])),
            0.0,
            10.0,
        );
    });
    let source = ResolvedScalarHysteresis::resolve(
        dynamic,
        &point(300.0)
            .with_quantity(
                "omega",
                QuantitySpec::dimensional(Dims([0, 0, -1, 0, 0, 0])),
                1.0,
            )
            .unwrap(),
        roles(),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    assert!(attempt(source, viscosity_keys()).is_err());
}

#[test]
fn repeated_viscosity_roles_share_one_explicitly_pinned_original_claim() {
    let card = viscosity_card([0.4, 0.9], |_| {});
    let mut claims = card.claims().clone();
    let key = viscosity_keys()[0].clone();
    let chosen = claims
        .insert_claim(PropertyClaim {
            key: key.clone(),
            value: PropertyValue::Scalar {
                value: 0.7,
                dims: MAGNETIC_VISCOSITY_DIMS,
            },
            validity: ValidityDomain::unconstrained().with_quantity(
                "T",
                QuantitySpec::dimensional(Temperature::DIMS),
                280.0,
                420.0,
            ),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![],
            provenance: Provenance {
                source: "separately pinned intrinsic viscosity".into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        })
        .unwrap();
    let card = Arc::new(MaterialCard::assemble(card.id().clone(), claims, vec![]).unwrap());
    let source = resolve(card);
    let repeated = vec![key.clone(), key.clone()];
    assert!(
        ResolvedVectorHysteresis::from_isotropic_source_with_viscosity(
            source.clone(),
            repeated.clone(),
            MaterialPropertySelection::SingleClaimOnly,
        )
        .is_err()
    );
    let fit = ResolvedVectorHysteresis::from_isotropic_source_with_viscosity(
        source.clone(),
        repeated,
        MaterialPropertySelection::PinnedByProperty(vec![(key.name().into(), chosen)]),
    )
    .unwrap();
    assert_eq!(fit.law().viscosities_m_s_h(), Some([0.7, 0.7].as_slice()));
    assert_eq!(
        fit.material().properties().len(),
        source.material().properties().len() + 1
    );
    let warm = fit.at_query_point(&point(400.0)).unwrap();
    assert_eq!(
        warm.material()
            .property_by_key(&key)
            .unwrap()
            .answer()
            .receipt
            .selected,
        chosen
    );
    assert_eq!(warm.law(), fit.law());
}
