//! G0/G1 source admission and cyclic response for independently fitted C/gamma.
//! Synthetic numerical fixtures verify integration, not empirical calibration.

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::SmallStrainLaw;
use fs_material::plastic::chaboche::source::{
    CHABOCHE_BACKSTRESS_MODULUS_PROPERTY as C, CHABOCHE_DYNAMIC_RECOVERY_PROPERTY as GAMMA,
    ResolvedChabocheJ2, chaboche_requirements, resolve_chaboche_j2,
};
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
const TANGENT: HardeningPropertyBasis = HardeningPropertyBasis::UniaxialTangentModulus;

fn scalar(name: &str, dims: Dims, value: f64) -> PropertyClaim {
    PropertyClaim {
        key: PropertyKey::new(name, dims),
        value: PropertyValue::Scalar { value, dims },
        validity: ValidityDomain::unconstrained(),
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        uncertainty: UncertaintyModel::Unstated,
        observations: vec![],
        provenance: Provenance {
            source: "synthetic independently specified Chaboche fixture".into(),
            license: "test-only".into(),
            artifact: None,
        },
    }
}

fn properties() -> Vec<PropertyClaim> {
    vec![
        scalar(YOUNG_MODULUS_PROPERTY, Pressure::DIMS, 1200.0),
        scalar(POISSON_RATIO_PROPERTY, Dims([0; 6]), 0.2),
        scalar(YIELD_STRESS_PROPERTY, Pressure::DIMS, 3.0),
        scalar(ISOTROPIC_HARDENING_MODULUS_PROPERTY, Pressure::DIMS, 40.0),
        scalar(Q, Pressure::DIMS, 2.0),
        scalar(B, Dims([0; 6]), 80.0),
        scalar(KINEMATIC_HARDENING_MODULUS_PROPERTY, Pressure::DIMS, 160.0),
        scalar(C, Pressure::DIMS, 300.0),
        scalar(GAMMA, Dims([0; 6]), 150.0),
    ]
}

fn card(properties: Vec<PropertyClaim>) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for property in properties {
        claims.insert_claim(property).unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic cyclic hardening solid".into(),
            phase: "solid".into(),
            process: "numerical fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}

fn resolve(source: &MaterialCard, basis: HardeningPropertyBasis) -> ResolvedChabocheJ2 {
    resolve_chaboche_j2(
        source,
        &QueryPoint::new(),
        MaterialPropertySelection::SingleClaimOnly,
        basis,
        0.04,
    )
    .unwrap()
}

fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-10 * expected.abs().max(1.0),
        "{actual} != {expected}",
    );
}

#[test]
fn independently_sourced_recovery_drives_cyclic_history_without_changing_the_voce_split() {
    let source = card(properties());
    for (basis, hk, count) in [(ISO, 0.0, 8), (MIXED, 160.0, 9)] {
        let resolved = resolve(&source, basis);
        assert_eq!(resolved.material().card_identity(), source.content_hash());
        assert_eq!(resolved.material().properties().len(), count);
        assert_eq!(resolved.hardening_basis(), basis);
        assert_eq!(resolved.law().branches()[0].modulus_pa(), 300.0);
        assert_eq!(resolved.law().branches()[0].recovery_per_strain(), 150.0);
        assert_eq!(resolved.backbone().law().saturation_pa(), 2.0);
        assert_eq!(resolved.backbone().law().rate_per_strain(), 80.0);
        assert_eq!(resolved.backbone().backbone().hardening_modulus_pa(), 40.0);
        assert_eq!(
            resolved
                .backbone()
                .backbone()
                .kinematic_hardening_modulus_pa(),
            hk
        );
        assert_eq!(resolved.modulus_source().value_si(), 300.0);
        assert_eq!(resolved.recovery_source().value_si(), 150.0);

        let mut state = resolved.initial_state();
        for shear in [0.01, -0.012, 0.014] {
            let strain = [0.0, 0.0, 0.0, shear, 0.0, 0.0];
            let stress = resolved.stress(&strain, &state);
            let next = resolved.update_state(&strain, &state);
            let delta_alpha = next.plastic.alpha - state.plastic.alpha;
            assert!(delta_alpha > 0.0);
            // Check the independently selected branch against the implicit
            // Armstrong-Frederick evolution, including reverse plastic flow.
            let delta_plastic = next.plastic.plastic_strain[3] - state.plastic.plastic_strain[3];
            near(
                next.backstress_pa[0][3],
                (state.backstress_pa[0][3] + (2.0 / 3.0) * 300.0 * delta_plastic)
                    / (1.0 + 150.0 * delta_alpha),
            );
            let shift =
                (2.0 / 3.0) * hk * next.plastic.plastic_strain[3] + next.backstress_pa[0][3];
            near(
                3.0_f64.sqrt() * (stress[3] - shift).abs(),
                3.0 + 40.0 * next.plastic.alpha + 2.0 * (1.0 - (-80.0 * next.plastic.alpha).exp()),
            );
            assert_eq!(stress, resolved.law().stress(&strain, &state));
            assert_eq!(
                resolved.tangent(&strain, &state),
                resolved.law().tangent(&strain, &state)
            );
            assert_eq!(next, resolved.law().update_state(&strain, &state));
            state = next;
        }
        for property in resolved.material().properties() {
            source
                .claims()
                .verify_receipt(&property.answer().receipt)
                .unwrap();
        }
    }
}

#[test]
fn complete_bundles_keep_voce_other_physics_and_original_receipts_without_tangent_inference() {
    let mut all = properties();
    all.push(scalar(DENSITY_PROPERTY, Density::DIMS, 7800.0));
    all.push(scalar(
        POST_YIELD_TANGENT_MODULUS_PROPERTY,
        Pressure::DIMS,
        500.0,
    ));
    let source = card(all);
    let mut requirements = chaboche_requirements(MIXED).unwrap();
    requirements.push(
        ScalarPropertyRequirement::try_new(
            DENSITY_PROPERTY,
            Density::DIMS,
            ScalarAdmissibility::StrictlyPositive,
        )
        .unwrap(),
    );
    requirements.push(TANGENT.requirement().unwrap());
    let bundle = resolve_material_state_point(
        &source,
        &QueryPoint::new(),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let resolved = ResolvedChabocheJ2::try_from_material(bundle.clone(), MIXED, 0.04).unwrap();
    assert_eq!(resolved.material(), &bundle);
    assert_eq!(resolved.backbone().material(), &bundle);
    assert_eq!(resolved.backbone().backbone().material(), &bundle);
    assert_eq!(resolved.material().properties().len(), 11);
    assert_eq!(resolved.modulus_source(), bundle.property(C).unwrap());
    assert_eq!(resolved.recovery_source(), bundle.property(GAMMA).unwrap());
    assert_eq!(
        resolved.backbone().saturation_source(),
        bundle.property(Q).unwrap()
    );
    assert_eq!(
        resolved.backbone().rate_source(),
        bundle.property(B).unwrap()
    );
    assert_eq!(
        resolved
            .material()
            .property(DENSITY_PROPERTY)
            .unwrap()
            .value_si(),
        7800.0
    );
    assert!(ResolvedChabocheJ2::try_from_material(bundle, TANGENT, 0.04).is_err());
    assert!(chaboche_requirements(TANGENT).is_err());
    assert!(
        resolve_chaboche_j2(
            &source,
            &QueryPoint::new(),
            MaterialPropertySelection::SingleClaimOnly,
            TANGENT,
            0.04,
        )
        .is_err()
    );
}

#[test]
fn each_ambiguous_branch_parameter_needs_an_explicit_pin_and_preserves_its_receipt() {
    let original = card(properties());
    let selected = resolve(&original, MIXED);
    for (name, dims, value) in [(C, Pressure::DIMS, 900.0), (GAMMA, Dims([0; 6]), 15.0)] {
        let pins = selected
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
        let mut ambiguous = properties();
        ambiguous.push(scalar(name, dims, value));
        let source = card(ambiguous);
        assert!(
            resolve_chaboche_j2(
                &source,
                &QueryPoint::new(),
                MaterialPropertySelection::SingleClaimOnly,
                MIXED,
                0.04,
            )
            .is_err()
        );
        let resolved = resolve_chaboche_j2(
            &source,
            &QueryPoint::new(),
            MaterialPropertySelection::PinnedByProperty(pins),
            MIXED,
            0.04,
        )
        .unwrap();
        assert_eq!(resolved.law().branches()[0].modulus_pa(), 300.0);
        assert_eq!(resolved.law().branches()[0].recovery_per_strain(), 150.0);
        assert_eq!(
            resolved
                .material()
                .property(name)
                .unwrap()
                .answer()
                .receipt
                .selected,
            selected
                .material()
                .property(name)
                .unwrap()
                .answer()
                .receipt
                .selected
        );
        assert_eq!(resolved.material().card_identity(), source.content_hash());
        for property in resolved.material().properties() {
            source
                .claims()
                .verify_receipt(&property.answer().receipt)
                .unwrap();
        }
    }
}

#[test]
fn fixed_temperature_coefficients_keep_source_domains_and_refuse_extrapolation() {
    let mut all = properties();
    for (name, low, high) in [(C, 200.0, 400.0), (GAMMA, 100.0, 200.0)] {
        let property = all.iter_mut().find(|p| p.key.name() == name).unwrap();
        property.value = PropertyValue::Curve {
            abscissa: "T".into(),
            abscissa_dims: Temperature::DIMS,
            knots: vec![(300.0, low), (400.0, high)],
            dims: property.key.quantity().dims(),
        };
        property.validity = ValidityDomain::unconstrained().with("T", 300.0, 400.0);
        property.interpolation = InterpolationPolicy::LinearInside;
    }
    let source = card(all);
    let resolved = resolve_chaboche_j2(
        &source,
        &QueryPoint::new().with("T", 350.0).unwrap(),
        MaterialPropertySelection::SingleClaimOnly,
        MIXED,
        0.04,
    )
    .unwrap();
    near(resolved.law().branches()[0].modulus_pa(), 300.0);
    near(resolved.law().branches()[0].recovery_per_strain(), 150.0);
    for property in [resolved.modulus_source(), resolved.recovery_source()] {
        assert_eq!(property.source_axes(), &[("T".into(), Temperature::DIMS)]);
        assert_eq!(
            property.answer().receipt.query_point,
            vec![("T".into(), 350.0)]
        );
        source
            .claims()
            .verify_receipt(&property.answer().receipt)
            .unwrap();
    }
    assert!(
        resolve_chaboche_j2(
            &source,
            &QueryPoint::new().with("T", 401.0).unwrap(),
            MaterialPropertySelection::SingleClaimOnly,
            MIXED,
            0.04,
        )
        .is_err()
    );
}

#[test]
fn missing_wrong_quantity_and_negative_branch_data_refuse_but_explicit_limits_are_admitted() {
    for missing in [C, GAMMA, Q, B] {
        let mut all = properties();
        all.retain(|p| p.key.name() != missing);
        assert!(
            resolve_chaboche_j2(
                &card(all),
                &QueryPoint::new(),
                MaterialPropertySelection::SingleClaimOnly,
                ISO,
                0.04,
            )
            .is_err()
        );
    }
    for (name, dims, value) in [
        (C, Pressure::DIMS, -1.0),
        (C, Density::DIMS, 300.0),
        (GAMMA, Dims([0; 6]), -1.0),
        (GAMMA, Frequency::DIMS, 150.0),
    ] {
        let mut all = properties();
        *all.iter_mut().find(|p| p.key.name() == name).unwrap() = scalar(name, dims, value);
        assert!(
            resolve_chaboche_j2(
                &card(all),
                &QueryPoint::new(),
                MaterialPropertySelection::SingleClaimOnly,
                ISO,
                0.04,
            )
            .is_err()
        );
    }
    for (modulus, recovery) in [(0.0, 150.0), (300.0, 0.0), (0.0, 0.0)] {
        let mut all = properties();
        *all.iter_mut().find(|p| p.key.name() == C).unwrap() = scalar(C, Pressure::DIMS, modulus);
        *all.iter_mut().find(|p| p.key.name() == GAMMA).unwrap() =
            scalar(GAMMA, Dims([0; 6]), recovery);
        let resolved = resolve(&card(all), MIXED);
        let mut state = resolved.initial_state();
        for shear in [0.01, 0.008, -0.012] {
            let strain = [0.0, 0.0, 0.0, shear, 0.0, 0.0];
            let next = resolved.update_state(&strain, &state);
            if modulus == 0.0 {
                assert_eq!(
                    resolved.stress(&strain, &state),
                    resolved.backbone().stress(&strain, &state.plastic)
                );
                assert_eq!(
                    resolved.tangent(&strain, &state),
                    resolved.backbone().tangent(&strain, &state.plastic)
                );
                assert_eq!(
                    next.plastic,
                    resolved.backbone().update_state(&strain, &state.plastic)
                );
                assert_eq!(next.backstress_pa, [[0.0; 6]]);
            } else {
                near(
                    next.backstress_pa[0][3],
                    (2.0 / 3.0) * modulus * next.plastic.plastic_strain[3],
                );
            }
            state = next;
        }
    }
}

#[test]
fn branch_coefficients_cannot_freeze_evolving_mechanical_history_or_rate_conditions() {
    let static_quantity = |kind| QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static));
    for name in [C, GAMMA] {
        for quantity in [
            static_quantity(QuantityKind::Stress),
            static_quantity(QuantityKind::Strain {
                basis: StrainBasis::Engineering,
                component: StrainComponent::Normal,
            }),
            QuantitySpec::dimensional(Frequency::DIMS),
            QuantitySpec::dimensional(Dims([0; 6])),
        ] {
            let mut all = properties();
            let property = all.iter_mut().find(|p| p.key.name() == name).unwrap();
            property.validity =
                ValidityDomain::unconstrained().with_quantity("condition", quantity, 0.0, 1.0);
            let source = card(all);
            let point = QueryPoint::new()
                .with_quantity("condition", quantity, 0.1)
                .unwrap();
            let bundle = resolve_material_state_point(
                &source,
                &point,
                &chaboche_requirements(ISO).unwrap(),
                MaterialPropertySelection::SingleClaimOnly,
            )
            .unwrap();
            assert!(ResolvedChabocheJ2::try_from_material(bundle, ISO, 0.04).is_err());
        }
    }
    let mut all = properties();
    let property = all.iter_mut().find(|p| p.key.name() == GAMMA).unwrap();
    property.value = PropertyValue::Curve {
        abscissa: "strain".into(),
        abscissa_dims: Dims([0; 6]),
        knots: vec![(0.0, 150.0), (0.02, 100.0)],
        dims: Dims([0; 6]),
    };
    property.validity = ValidityDomain::unconstrained().with("strain", 0.0, 0.02);
    property.interpolation = InterpolationPolicy::LinearInside;
    assert!(
        resolve_chaboche_j2(
            &card(all),
            &QueryPoint::new().with("strain", 0.01).unwrap(),
            MaterialPropertySelection::SingleClaimOnly,
            ISO,
            0.04,
        )
        .is_err()
    );
}
