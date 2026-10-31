//! Named, independently sourced backstresses through the original cyclic law.
//! Synthetic verification data do not establish empirical material calibration.

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::SmallStrainLaw;
use fs_material::plastic::chaboche::source::{
    BackstressProperties, CHABOCHE_BACKSTRESS_MODULUS_PROPERTY as C,
    CHABOCHE_DYNAMIC_RECOVERY_PROPERTY as GAMMA, ResolvedChabocheBranches,
    chaboche_branch_requirements, chaboche_requirements, resolve_chaboche_branches,
    resolve_chaboche_j2,
};
use fs_material::plastic::energy::EnergyTolerances;
use fs_material::plastic::source::{
    HardeningPropertyBasis, ISOTROPIC_HARDENING_MODULUS_PROPERTY,
    KINEMATIC_HARDENING_MODULUS_PROPERTY, POST_YIELD_TANGENT_MODULUS_PROPERTY,
};
use fs_material::plastic::voce::source::{
    VOCE_SATURATION_INCREMENT_PROPERTY as Q, VOCE_SATURATION_RATE_PROPERTY as B,
};
use fs_material::state_point::{
    DENSITY_PROPERTY, MaterialPropertySelection, POISSON_RATIO_PROPERTY, ScalarAdmissibility,
    ScalarPropertyRequirement, YIELD_STRESS_PROPERTY, YOUNG_MODULUS_PROPERTY,
    resolve_material_state_point,
};
use fs_qty::semantic::{QuantityKind, SemanticType, StrainBasis, StrainComponent, ValueForm};
use fs_qty::{Density, Dims, Frequency, Pressure, QuantitySpec, Temperature};

const ISO: HardeningPropertyBasis = HardeningPropertyBasis::PlasticStrainModulus;
const MIXED: HardeningPropertyBasis = HardeningPropertyBasis::MixedPlasticStrainModuli;
const E: f64 = 1200.0;
const SY: f64 = 3.0;
const H: f64 = 40.0;
const HK: f64 = 160.0;

fn names<const N: usize>() -> [BackstressProperties; N] {
    // Descending names expose accidental property-sort/history-sort coupling.
    core::array::from_fn(|i| {
        BackstressProperties::new(
            format!("backstress_{}_modulus", N - i),
            format!("backstress_{}_recovery", N - i),
        )
    })
}

fn scalar(name: &str, dims: Dims, value: f64) -> PropertyClaim {
    PropertyClaim {
        key: PropertyKey::new(name, dims),
        value: PropertyValue::Scalar { value, dims },
        validity: ValidityDomain::unconstrained(),
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        uncertainty: UncertaintyModel::Unstated,
        observations: vec![],
        provenance: Provenance {
            source: "synthetic named cyclic branches".into(),
            license: "test-only".into(),
            artifact: None,
        },
    }
}

fn properties<const N: usize>(
    schema: &[BackstressProperties; N],
    values: &[(f64, f64); N],
) -> Vec<PropertyClaim> {
    let mut data = vec![
        scalar(YOUNG_MODULUS_PROPERTY, Pressure::DIMS, E),
        scalar(POISSON_RATIO_PROPERTY, Dims([0; 6]), 0.2),
        scalar(YIELD_STRESS_PROPERTY, Pressure::DIMS, SY),
        scalar(ISOTROPIC_HARDENING_MODULUS_PROPERTY, Pressure::DIMS, H),
        scalar(KINEMATIC_HARDENING_MODULUS_PROPERTY, Pressure::DIMS, HK),
        scalar(Q, Pressure::DIMS, 2.0),
        scalar(B, Dims([0; 6]), 80.0),
    ];
    for (names, &(c, gamma)) in schema.iter().zip(values) {
        data.push(scalar(&names.modulus, Pressure::DIMS, c));
        data.push(scalar(&names.recovery, Dims([0; 6]), gamma));
    }
    data
}

fn card(data: Vec<PropertyClaim>) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for property in data {
        claims.insert_claim(property).unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic multi-backstress solid".into(),
            phase: "solid".into(),
            process: "numerical fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}

fn resolve<const N: usize>(
    card: &MaterialCard,
    schema: [BackstressProperties; N],
) -> ResolvedChabocheBranches<N> {
    resolve_chaboche_branches(
        card,
        &QueryPoint::new(),
        MaterialPropertySelection::SingleClaimOnly,
        MIXED,
        0.04,
        schema,
    )
    .unwrap()
}

fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-10 * expected.abs().max(1.0),
        "{actual:e} != {expected:e}"
    );
}

fn cyclic<const N: usize>(values: [(f64, f64); N]) {
    let schema = names();
    let source = card(properties(&schema, &values));
    let resolved = resolve(&source, schema.clone());
    assert_eq!(resolved.properties(), &schema);
    assert_eq!(resolved.hardening_basis(), MIXED);
    assert_eq!(resolved.material().card_identity(), source.content_hash());
    assert_eq!(resolved.material().properties().len(), 7 + 2 * N);
    for (i, ((c, gamma), names)) in values.iter().zip(&schema).enumerate() {
        assert_eq!(resolved.law().branches()[i].modulus_pa(), *c);
        assert_eq!(resolved.law().branches()[i].recovery_per_strain(), *gamma);
        assert_eq!(
            resolved.modulus_sources()[i].requirement().name(),
            names.modulus
        );
        assert_eq!(
            resolved.recovery_sources()[i].requirement().name(),
            names.recovery
        );
    }
    let mut old_strain = [0.0; 6];
    let mut state = resolved.initial_state();
    for shear in [0.012, -0.014, 0.016, -0.009] {
        let strain = [0.0, 0.0, 0.0, shear, 0.0, 0.0];
        let stress = resolved.stress(&strain, &state);
        let next = resolved.update_state(&strain, &state);
        let delta = next.plastic.alpha - state.plastic.alpha;
        assert!(delta > 0.0);
        let mut shift = (2.0 / 3.0) * HK * next.plastic.plastic_strain[3];
        let mut expected_recovery = 0.0;
        for (i, &(c, gamma)) in values.iter().enumerate() {
            for component in 0..6 {
                let dp = next.plastic.plastic_strain[component]
                    - state.plastic.plastic_strain[component];
                near(
                    next.backstress_pa[i][component],
                    (state.backstress_pa[i][component] + (2.0 / 3.0) * c * dp)
                        / (1.0 + gamma * delta),
                );
            }
            shift += next.backstress_pa[i][3];
            if c > 0.0 {
                expected_recovery += 3.0 * gamma * delta * next.backstress_pa[i][3].powi(2) / c;
            } else {
                assert_eq!(next.backstress_pa[i], [0.0; 6]);
            }
        }
        near(
            3.0_f64.sqrt() * (stress[3] - shift).abs(),
            SY + H * next.plastic.alpha + 2.0 * -(-80.0 * next.plastic.alpha).exp_m1(),
        );
        assert_eq!(stress, resolved.law().stress(&strain, &state));
        assert_eq!(next, resolved.law().update_state(&strain, &state));
        assert_eq!(
            resolved.tangent(&strain, &state),
            resolved.law().tangent(&strain, &state)
        );
        let energy = resolved
            .law()
            .advance_with_energy(
                &old_strain,
                &state,
                &strain,
                EnergyTolerances {
                    energy_j_m3: 1e-9,
                    yield_pa: 1e-9,
                },
            )
            .unwrap();
        assert_eq!(energy.state, next);
        near(energy.recovery_dissipation_j_m3, expected_recovery);
        near(
            energy.energy.physical_dissipation_j_m3,
            SY * delta + expected_recovery,
        );
        old_strain = strain;
        state = next;
    }
    for property in resolved.material().properties() {
        source
            .claims()
            .verify_receipt(&property.answer().receipt)
            .unwrap();
    }
}

#[test]
fn two_and_three_independent_recovery_scales_survive_reversal_and_ordered_energy() {
    cyclic([(300.0, 800.0), (120.0, 25.0)]);
    cyclic([(300.0, 800.0), (120.0, 80.0), (40.0, 8.0)]);
    cyclic([(300.0, 800.0), (0.0, 80.0), (40.0, 0.0)]);
}

#[test]
fn equal_total_modulus_does_not_erase_distinct_recovery_or_cyclic_response() {
    let schema = names::<2>();
    let values_a = [(300.0, 20.0), (100.0, 500.0)];
    let values_b = [(300.0, 150.0), (100.0, 150.0)];
    let a = resolve(&card(properties(&schema, &values_a)), schema.clone());
    let b = resolve(&card(properties(&schema, &values_b)), schema);
    assert_eq!(
        a.law()
            .branches()
            .iter()
            .map(|x| x.modulus_pa())
            .sum::<f64>(),
        b.law()
            .branches()
            .iter()
            .map(|x| x.modulus_pa())
            .sum::<f64>()
    );
    let mut state_a = a.initial_state();
    let mut state_b = b.initial_state();
    for shear in [0.02, -0.015, 0.018] {
        let strain = [0.0, 0.0, 0.0, shear, 0.0, 0.0];
        assert!((a.stress(&strain, &state_a)[3] - b.stress(&strain, &state_b)[3]).abs() > 0.01);
        state_a = a.update_state(&strain, &state_a);
        state_b = b.update_state(&strain, &state_b);
    }
}

#[test]
fn full_bundle_retains_density_extra_roles_identity_and_branch_source_pins() {
    let schema = names::<3>();
    let mut data = properties(&schema, &[(300.0, 800.0), (120.0, 80.0), (40.0, 8.0)]);
    data.push(scalar(DENSITY_PROPERTY, Density::DIMS, 7800.0));
    data.push(scalar(
        POST_YIELD_TANGENT_MODULUS_PROPERTY,
        Pressure::DIMS,
        999.0,
    ));
    let source = card(data);
    let mut requirements = chaboche_branch_requirements(MIXED, &schema).unwrap();
    requirements.push(
        ScalarPropertyRequirement::try_new(
            DENSITY_PROPERTY,
            Density::DIMS,
            ScalarAdmissibility::StrictlyPositive,
        )
        .unwrap(),
    );
    requirements.push(
        HardeningPropertyBasis::UniaxialTangentModulus
            .requirement()
            .unwrap(),
    );
    let bundle = resolve_material_state_point(
        &source,
        &QueryPoint::new(),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let resolved =
        ResolvedChabocheBranches::try_from_material(bundle.clone(), MIXED, 0.04, schema.clone())
            .unwrap();
    assert_eq!(resolved.material(), &bundle);
    assert_eq!(resolved.backbone().material(), &bundle);
    assert_eq!(resolved.material().identity(), bundle.identity());
    assert_eq!(
        resolved
            .material()
            .property(DENSITY_PROPERTY)
            .unwrap()
            .value_si(),
        7800.0
    );
    for (i, names) in schema.iter().enumerate() {
        assert_eq!(
            resolved.modulus_sources()[i],
            bundle.property(&names.modulus).unwrap()
        );
        assert_eq!(
            resolved.recovery_sources()[i],
            bundle.property(&names.recovery).unwrap()
        );
    }
    assert!(
        ResolvedChabocheBranches::try_from_material(
            bundle,
            HardeningPropertyBasis::UniaxialTangentModulus,
            0.04,
            schema
        )
        .is_err()
    );
}

#[test]
fn one_canonical_branch_matches_the_original_single_branch_api_exactly() {
    let schema = [BackstressProperties::new(C, GAMMA)];
    assert_eq!(
        chaboche_branch_requirements(MIXED, &schema).unwrap(),
        chaboche_requirements(MIXED).unwrap()
    );
    let source = card(properties(&schema, &[(300.0, 150.0)]));
    let original = resolve_chaboche_j2(
        &source,
        &QueryPoint::new(),
        MaterialPropertySelection::SingleClaimOnly,
        MIXED,
        0.04,
    )
    .unwrap();
    let named = resolve(&source, schema);
    assert_eq!(original.law(), named.law());
    assert_eq!(original.material(), named.material());
    assert_eq!(original.modulus_source(), named.modulus_sources()[0]);
    assert_eq!(original.recovery_source(), named.recovery_sources()[0]);
    let mut state = original.initial_state();
    for shear in [0.012, -0.014, 0.016] {
        let strain = [0.0, 0.0, 0.0, shear, 0.0, 0.0];
        assert_eq!(
            original.stress(&strain, &state),
            named.stress(&strain, &state)
        );
        assert_eq!(
            original.tangent(&strain, &state),
            named.tangent(&strain, &state)
        );
        assert_eq!(
            original.update_state(&strain, &state),
            named.update_state(&strain, &state)
        );
        state = named.update_state(&strain, &state);
    }
}

#[test]
fn counts_and_aliased_branch_backbone_or_density_names_refuse() {
    assert!(chaboche_branch_requirements(ISO, &names::<0>()).is_err());
    assert!(chaboche_branch_requirements(ISO, &names::<9>()).is_err());
    let eight = names::<8>();
    let values = core::array::from_fn(|i| (10.0 + i as f64, 30.0 * i as f64));
    let admitted = resolve(&card(properties(&eight, &values)), eight);
    assert_eq!(admitted.law().branches().len(), 8);
    for reserved in [
        YOUNG_MODULUS_PROPERTY,
        POISSON_RATIO_PROPERTY,
        YIELD_STRESS_PROPERTY,
        ISOTROPIC_HARDENING_MODULUS_PROPERTY,
        KINEMATIC_HARDENING_MODULUS_PROPERTY,
        POST_YIELD_TANGENT_MODULUS_PROPERTY,
        Q,
        B,
        DENSITY_PROPERTY,
        "",
        " white_space ",
    ] {
        let schema = [BackstressProperties::new(reserved, "own_gamma")];
        assert!(
            chaboche_branch_requirements(ISO, &schema).is_err(),
            "accepted {reserved}"
        );
        let schema = [BackstressProperties::new("own_C", reserved)];
        assert!(
            chaboche_branch_requirements(MIXED, &schema).is_err(),
            "accepted {reserved}"
        );
    }
    for schema in [
        [
            BackstressProperties::new("C1", "g1"),
            BackstressProperties::new("C1", "g2"),
        ],
        [
            BackstressProperties::new("C1", "g1"),
            BackstressProperties::new("C2", "g1"),
        ],
        [
            BackstressProperties::new("C1", "g1"),
            BackstressProperties::new("g1", "g2"),
        ],
        [
            BackstressProperties::new("same", "same"),
            BackstressProperties::new("C2", "g2"),
        ],
    ] {
        assert!(chaboche_branch_requirements(MIXED, &schema).is_err());
    }
}

#[test]
fn every_branch_requires_correct_units_nonnegative_values_and_unambiguous_selection() {
    let schema = names::<2>();
    let original = properties(&schema, &[(300.0, 800.0), (120.0, 25.0)]);
    for name in [&schema[1].modulus, &schema[1].recovery] {
        let mut missing = original.clone();
        missing.retain(|p| p.key.name() != name);
        assert!(
            resolve_chaboche_branches(
                &card(missing),
                &QueryPoint::new(),
                MaterialPropertySelection::SingleClaimOnly,
                MIXED,
                0.04,
                schema.clone()
            )
            .is_err()
        );
        let dims = if *name == schema[1].modulus {
            Pressure::DIMS
        } else {
            Dims([0; 6])
        };
        for (dims, value) in [(dims, -1.0), (Frequency::DIMS, 12.0)] {
            let mut invalid = original.clone();
            *invalid.iter_mut().find(|p| p.key.name() == name).unwrap() = scalar(name, dims, value);
            assert!(
                resolve_chaboche_branches(
                    &card(invalid),
                    &QueryPoint::new(),
                    MaterialPropertySelection::SingleClaimOnly,
                    MIXED,
                    0.04,
                    schema.clone()
                )
                .is_err()
            );
        }
    }
    let selected = resolve(&card(original.clone()), schema.clone());
    let pins: Vec<_> = selected
        .material()
        .properties()
        .iter()
        .map(|p| {
            (
                p.requirement().name().to_owned(),
                p.answer().receipt.selected,
            )
        })
        .collect();
    let mut duplicate = original;
    duplicate.push(scalar(&schema[1].modulus, Pressure::DIMS, 999.0));
    let source = card(duplicate);
    assert!(
        resolve_chaboche_branches(
            &source,
            &QueryPoint::new(),
            MaterialPropertySelection::SingleClaimOnly,
            MIXED,
            0.04,
            schema.clone()
        )
        .is_err()
    );
    let pinned = resolve_chaboche_branches(
        &source,
        &QueryPoint::new(),
        MaterialPropertySelection::PinnedByProperty(pins.clone()),
        MIXED,
        0.04,
        schema.clone(),
    )
    .unwrap();
    assert_eq!(pinned.law().branches()[1].modulus_pa(), 120.0);
    for property in pinned.material().properties() {
        source
            .claims()
            .verify_receipt(&property.answer().receipt)
            .unwrap();
    }
    let mut ambiguous = pins;
    ambiguous.push(
        ambiguous
            .iter()
            .find(|(name, _)| name == &schema[1].modulus)
            .unwrap()
            .clone(),
    );
    assert!(
        resolve_chaboche_branches(
            &source,
            &QueryPoint::new(),
            MaterialPropertySelection::PinnedByProperty(ambiguous),
            MIXED,
            0.04,
            schema
        )
        .is_err()
    );
}

#[test]
fn later_branches_cannot_freeze_evolving_coordinates_and_fixed_temperature_support_is_retained() {
    let schema = names::<2>();
    let original = properties(&schema, &[(300.0, 800.0), (120.0, 25.0)]);
    for quantity in [
        QuantitySpec::semantic(SemanticType::new(QuantityKind::Stress, ValueForm::Static)),
        QuantitySpec::semantic(SemanticType::new(
            QuantityKind::Strain {
                basis: StrainBasis::Engineering,
                component: StrainComponent::Normal,
            },
            ValueForm::Static,
        )),
        QuantitySpec::dimensional(Frequency::DIMS),
        QuantitySpec::dimensional(Dims([0; 6])),
    ] {
        let mut invalid = original.clone();
        let property = invalid
            .iter_mut()
            .find(|p| p.key.name() == schema[1].recovery)
            .unwrap();
        property.validity =
            ValidityDomain::unconstrained().with_quantity("condition", quantity, 0.0, 1.0);
        let point = QueryPoint::new()
            .with_quantity("condition", quantity, 0.1)
            .unwrap();
        assert!(
            resolve_chaboche_branches(
                &card(invalid),
                &point,
                MaterialPropertySelection::SingleClaimOnly,
                MIXED,
                0.04,
                schema.clone()
            )
            .is_err()
        );
    }
    let mut data = original;
    let property = data
        .iter_mut()
        .find(|p| p.key.name() == schema[1].modulus)
        .unwrap();
    property.value = PropertyValue::Curve {
        abscissa: "T".into(),
        abscissa_dims: Temperature::DIMS,
        knots: vec![(300.0, 100.0), (400.0, 200.0)],
        dims: Pressure::DIMS,
    };
    property.validity = ValidityDomain::unconstrained().with("T", 300.0, 400.0);
    property.interpolation = InterpolationPolicy::LinearInside;
    let source = card(data);
    let at = |t| QueryPoint::new().with("T", t).unwrap();
    let fixed = resolve_chaboche_branches(
        &source,
        &at(350.0),
        MaterialPropertySelection::SingleClaimOnly,
        MIXED,
        0.04,
        schema.clone(),
    )
    .unwrap();
    near(fixed.law().branches()[1].modulus_pa(), 150.0);
    assert_eq!(
        fixed.modulus_sources()[1].source_axes(),
        &[("T".into(), Temperature::DIMS)]
    );
    assert_eq!(
        fixed.modulus_sources()[1].answer().receipt.query_point,
        vec![("T".into(), 350.0)]
    );
    assert!(
        resolve_chaboche_branches(
            &source,
            &at(401.0),
            MaterialPropertySelection::SingleClaimOnly,
            MIXED,
            0.04,
            schema
        )
        .is_err()
    );
}
