//! Original source knots drive native plastic response without losing meaning.
//! Synthetic coefficients verify implementation, not alloy calibration.
#![deny(unsafe_code)]

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimId, ClaimSelection, ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId,
    PropertyClaim, PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy,
    UncertaintyModel,
};
use fs_material::SmallStrainLaw;
use fs_material::plastic::energy::EnergyTolerances;
use fs_material::plastic::section::energy::NativeSectionEnergyLaw;
use fs_material::plastic::source::KINEMATIC_HARDENING_MODULUS_PROPERTY;
use fs_material::plastic::tabulated::source::{
    ACCUMULATED_PLASTIC_STRAIN_AXIS as ALPHA, CAUCHY_YIELD_CURVE_PROPERTY, HardeningCurveRequest,
    KIRCHHOFF_YIELD_CURVE_PROPERTY, ResolvedTabulatedJ2, TabulatedJ2SourceRequest,
    TabulatedKinematicHardening, YieldStressMeasure, resolve_tabulated_hardening,
    resolve_tabulated_j2, tabulated_j2_requirements,
};
use fs_material::state_point::{
    DENSITY_PROPERTY, MaterialPropertySelection, POISSON_RATIO_PROPERTY, ScalarAdmissibility,
    ScalarPropertyRequirement, YOUNG_MODULUS_PROPERTY, resolve_material_state_point,
};
use fs_qty::semantic::{QuantityKind, SemanticType, StrainBasis, StrainComponent, ValueForm};
use fs_qty::{Density, Dims, Pressure, QuantitySpec};

const TEMPERATURE: &str = "specimen-temperature";
const SINGLE: ClaimSelection = ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly);
const KNOTS: [(f64, f64); 5] = [
    (0.0, 3.0),
    (0.002, 4.0),
    (0.006, 4.0),
    (0.012, 7.0),
    (0.04, 8.4),
];

fn absolute_temperature() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}
fn context() -> QueryPoint {
    QueryPoint::new()
        .with_quantity(TEMPERATURE, absolute_temperature(), 300.0)
        .unwrap()
}
fn request(measure: YieldStressMeasure) -> HardeningCurveRequest {
    HardeningCurveRequest {
        stress_measure: measure,
        maximum_points: KNOTS.len(),
    }
}
fn scalar(name: &str, dims: Dims, value: f64) -> PropertyClaim {
    PropertyClaim {
        key: PropertyKey::new(name, dims),
        value: PropertyValue::Scalar { value, dims },
        validity: ValidityDomain::unconstrained().with_quantity(
            TEMPERATURE,
            absolute_temperature(),
            290.0,
            330.0,
        ),
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        uncertainty: UncertaintyModel::Unstated,
        observations: vec![],
        provenance: Provenance {
            source: "synthetic hardening calibration".into(),
            license: "test-only".into(),
            artifact: None,
        },
    }
}
fn curve(measure: YieldStressMeasure) -> PropertyClaim {
    let mut claim = scalar(measure.property_name(), Pressure::DIMS, 3.0);
    claim.value = PropertyValue::Curve {
        abscissa: ALPHA.into(),
        abscissa_dims: Dims::NONE,
        knots: KNOTS.to_vec(),
        dims: Pressure::DIMS,
    };
    claim.validity =
        claim
            .validity
            .with_quantity(ALPHA, QuantitySpec::dimensional(Dims::NONE), 0.0, 0.04);
    claim.interpolation = InterpolationPolicy::TabulatedOnly;
    claim.uncertainty = UncertaintyModel::HalfWidth {
        half_width: 0.05,
        confidence: 0.95,
    };
    claim
}
fn properties() -> Vec<PropertyClaim> {
    vec![
        scalar(YOUNG_MODULUS_PROPERTY, Pressure::DIMS, 1200.0),
        scalar(POISSON_RATIO_PROPERTY, Dims::NONE, 0.2),
        scalar(KINEMATIC_HARDENING_MODULUS_PROPERTY, Pressure::DIMS, 160.0),
        scalar(DENSITY_PROPERTY, Density::DIMS, 7800.0),
        curve(YieldStressMeasure::Cauchy),
    ]
}
fn card(properties: Vec<PropertyClaim>) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for property in properties {
        claims.insert_claim(property).unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic alloy".into(),
            phase: "solid".into(),
            process: "explicit tabulated flow fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}"
    );
}

#[test]
fn complete_original_knots_and_interpolation_policy_survive_preparation() {
    let original = curve(YieldStressMeasure::Cauchy);
    let source = card(vec![original.clone()]);
    let r = resolve_tabulated_hardening(
        &source,
        &context(),
        &request(YieldStressMeasure::Cauchy),
        SINGLE,
    )
    .unwrap();
    assert_eq!(r.source().card_identity(), source.content_hash());
    assert_eq!(r.source().material_state(), source.id());
    assert_eq!(r.source().context(), &context());
    assert_eq!(r.source().claim(), &original);
    assert_eq!(r.source().selection(), SINGLE);
    assert_eq!(r.source().request(), &request(YieldStressMeasure::Cauchy));
    assert_eq!(r.stress_measure(), YieldStressMeasure::Cauchy);
    assert_eq!(r.source().observations().len(), KNOTS.len());
    source
        .claims()
        .verify_receipt(&r.source().selection_answer().receipt)
        .unwrap();
    for ((point, answer), &(alpha, stress)) in r.source().observations().iter().zip(&KNOTS) {
        assert_eq!(point.axes()[ALPHA], alpha);
        assert_eq!(
            point.axis_quantities()[ALPHA],
            QuantitySpec::dimensional(Dims::NONE)
        );
        assert_eq!(answer.evidence.value.value, stress);
        source.claims().verify_receipt(&answer.receipt).unwrap();
    }
    // Independent line segment: Y(.009)=4+500*(.009-.006)=5.5 Pa.
    near(r.hardening().yield_stress_pa(0.009).unwrap(), 5.5, 2e-15);
    assert!(r.hardening().yield_stress_pa(0.040001).is_err());
    let at = context()
        .with_quantity(ALPHA, QuantitySpec::dimensional(Dims::NONE), 0.009)
        .unwrap();
    assert!(
        source
            .claims()
            .query_typed(&original.key, &at, SelectionPolicy::SingleClaimOnly)
            .is_err(),
        "constitutive interpolation must not rewrite the original TabulatedOnly claim"
    );
    assert_eq!(
        r,
        resolve_tabulated_hardening(
            &source,
            &context(),
            &request(YieldStressMeasure::Cauchy),
            SINGLE
        )
        .unwrap()
    );
}

#[test]
fn sourced_curve_and_independent_prager_drive_tensor_flow_and_native_section_energy() {
    let source = card(properties());
    let mut requirements = tabulated_j2_requirements(TabulatedKinematicHardening::Prager).unwrap();
    requirements.push(
        ScalarPropertyRequirement::try_new(
            DENSITY_PROPERTY,
            Density::DIMS,
            ScalarAdmissibility::StrictlyPositive,
        )
        .unwrap(),
    );
    let bundle = resolve_material_state_point(
        &source,
        &context(),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let prepared = resolve_tabulated_hardening(
        &source,
        &context(),
        &request(YieldStressMeasure::Cauchy),
        SINGLE,
    )
    .unwrap();
    let resolved = ResolvedTabulatedJ2::try_from_material(
        bundle.clone(),
        prepared,
        TabulatedKinematicHardening::Prager,
        0.05,
    )
    .unwrap();
    assert_eq!(resolved.material(), &bundle);
    assert_eq!(resolved.material().properties().len(), 4);
    assert_eq!(
        resolved
            .material()
            .property(DENSITY_PROPERTY)
            .unwrap()
            .value_si(),
        7800.0
    );
    for property in resolved.material().properties() {
        source
            .claims()
            .verify_receipt(&property.answer().receipt)
            .unwrap();
    }
    // Prescribe accumulated flow, then independently construct the required
    // total deviatoric strain. mu=500; sigma_eq = Y(alpha)+Hk*alpha.
    let alpha = 0.009;
    let effective = 5.5 + 160.0 * alpha;
    let axial = alpha + effective / 1500.0;
    let strain = [axial, -0.5 * axial, -0.5 * axial, 0.0, 0.0, 0.0];
    let virgin = resolved.initial_state();
    let r = resolved.response(&strain, &virgin).unwrap();
    near(r.state.alpha, alpha, 2e-16);
    near(r.stress_pa[0], 2.0 * effective / 3.0, 2e-12);
    for column in 0..6 {
        let mut plus = strain;
        let mut minus = strain;
        let h = 1e-7;
        plus[column] += h;
        minus[column] -= h;
        let sp = resolved.stress(&plus, &virgin);
        let sm = resolved.stress(&minus, &virgin);
        for row in 0..6 {
            near(
                r.tangent_pa[row][column],
                (sp[row] - sm[row]) / (2.0 * h),
                2e-5,
            );
        }
    }
    let gates = EnergyTolerances {
        energy_j_m3: 1e-10,
        yield_pa: 1e-9,
    };
    let (state, stress, energy) = resolved
        .native_energy_increment(&[0.0; 6], &virgin, &strain, gates)
        .unwrap();
    assert_eq!(state, r.state);
    assert_eq!(stress, r.stress_pa);
    near(energy.physical_dissipation_j_m3, 3.0 * alpha, 1e-14);
    near(
        energy.final_state.hardening_j_m3,
        0.001 + 0.004 + 0.00525,
        1e-14,
    );
    near(
        energy.final_state.kinematic_hardening_j_m3,
        0.5 * 160.0 * alpha * alpha,
        1e-14,
    );
    near(energy.residual_j_m3, 0.0, 1e-10);
    resolved
        .validate_committed_history(&strain, &state, gates.yield_pa)
        .unwrap();
    let simple = resolve_tabulated_j2(
        &source,
        &context(),
        &TabulatedJ2SourceRequest {
            kinematic_hardening: TabulatedKinematicHardening::None,
            strain_limit: 0.05,
            maximum_points: 5,
        },
        MaterialPropertySelection::SingleClaimOnly,
        SINGLE,
    )
    .unwrap();
    assert_eq!(simple.law().backbone().kinematic_hardening, 0.0);
    assert_eq!(simple.material().properties().len(), 2);
}

#[test]
fn source_meaning_complete_support_and_preparation_bounds_are_enforced() {
    let original = curve(YieldStressMeasure::Cauchy);
    let source = card(vec![original.clone()]);
    assert!(
        resolve_tabulated_hardening(
            &source,
            &context(),
            &HardeningCurveRequest {
                maximum_points: 4,
                ..request(YieldStressMeasure::Cauchy)
            },
            SINGLE
        )
        .is_err()
    );
    assert!(
        resolve_tabulated_hardening(
            &source,
            &context(),
            &request(YieldStressMeasure::Kirchhoff),
            SINGLE
        )
        .is_err()
    );
    let mut wrong_axis = original.clone();
    if let PropertyValue::Curve { abscissa, .. } = &mut wrong_axis.value {
        *abscissa = "engineering-total-strain".into();
    }
    let mut wrong_quantity = original.clone();
    wrong_quantity.key = PropertyKey::new(CAUCHY_YIELD_CURVE_PROPERTY, Dims::NONE);
    if let PropertyValue::Curve { dims, .. } = &mut wrong_quantity.value {
        *dims = Dims::NONE;
    }
    let mut wrong_convention = original.clone();
    wrong_convention.validity = wrong_convention.validity.with_quantity(
        ALPHA,
        QuantitySpec::semantic(SemanticType::new(
            QuantityKind::Strain {
                basis: StrainBasis::Engineering,
                component: StrainComponent::Normal,
            },
            ValueForm::Static,
        )),
        0.0,
        0.04,
    );
    let mut clipped = original.clone();
    clipped.validity =
        clipped
            .validity
            .with_quantity(ALPHA, QuantitySpec::dimensional(Dims::NONE), 0.0, 0.03);
    let mut softening = original.clone();
    if let PropertyValue::Curve { knots, .. } = &mut softening.value {
        knots[2].1 = 3.5;
    }
    for invalid in [
        wrong_axis,
        wrong_quantity,
        wrong_convention,
        clipped,
        softening,
    ] {
        assert!(
            resolve_tabulated_hardening(
                &card(vec![invalid]),
                &context(),
                &request(YieldStressMeasure::Cauchy),
                SINGLE
            )
            .is_err()
        );
    }
    let runtime = context()
        .with_quantity(
            "strain-rate",
            QuantitySpec::dimensional(Dims([0, 0, -1, 0, 0, 0])),
            1e-3,
        )
        .unwrap();
    assert!(
        resolve_tabulated_hardening(
            &source,
            &runtime,
            &request(YieldStressMeasure::Cauchy),
            SINGLE
        )
        .is_err()
    );
    let alpha_in_context = context()
        .with_quantity(ALPHA, QuantitySpec::dimensional(Dims::NONE), 0.0)
        .unwrap();
    assert!(
        resolve_tabulated_hardening(
            &source,
            &alpha_in_context,
            &request(YieldStressMeasure::Cauchy),
            SINGLE
        )
        .is_err()
    );
    let missing_temperature = QueryPoint::new();
    assert!(
        resolve_tabulated_hardening(
            &source,
            &missing_temperature,
            &request(YieldStressMeasure::Cauchy),
            SINGLE
        )
        .is_err()
    );
}

#[test]
fn source_ambiguity_needs_one_explicit_claim_pin_for_the_complete_curve() {
    let first = curve(YieldStressMeasure::Cauchy);
    let mut second = first.clone();
    second.provenance.source = "a distinct synthetic flow calibration".into();
    if let PropertyValue::Curve { knots, .. } = &mut second.value {
        for (_, y) in knots {
            *y *= 2.0;
        }
    }
    let pins = [first.content_hash(), second.content_hash()];
    let source = card(vec![first, second]);
    assert!(
        resolve_tabulated_hardening(
            &source,
            &context(),
            &request(YieldStressMeasure::Cauchy),
            SINGLE
        )
        .is_err()
    );
    for (index, pin) in pins.into_iter().enumerate() {
        let selection = ClaimSelection::Pinned(ClaimId(pin));
        let r = resolve_tabulated_hardening(
            &source,
            &context(),
            &request(YieldStressMeasure::Cauchy),
            selection,
        )
        .unwrap();
        assert_eq!(r.source().selection(), selection);
        near(
            r.hardening().yield_stress_pa(0.009).unwrap(),
            5.5 * (index + 1) as f64,
            4e-15,
        );
        assert!(
            r.source()
                .observations()
                .iter()
                .all(|(_, answer)| answer.receipt.selected == ClaimId(pin))
        );
    }
}

#[test]
fn same_material_context_and_stress_measure_are_required_for_native_binding() {
    let source = card(properties());
    let requirements = tabulated_j2_requirements(TabulatedKinematicHardening::None).unwrap();
    let material = resolve_material_state_point(
        &source,
        &context(),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let bind = |curve| {
        ResolvedTabulatedJ2::try_from_material(
            material.clone(),
            curve,
            TabulatedKinematicHardening::None,
            0.05,
        )
    };
    let hot = resolve_tabulated_hardening(
        &source,
        &context().with(TEMPERATURE, 301.0).unwrap(),
        &request(YieldStressMeasure::Cauchy),
        SINGLE,
    )
    .unwrap();
    assert!(
        bind(hot).is_err(),
        "same values at another temperature cannot transfer source context"
    );
    let mut different = properties();
    different[0].provenance.source = "another elastic calibration".into();
    let other_card = card(different);
    let other_curve = resolve_tabulated_hardening(
        &other_card,
        &context(),
        &request(YieldStressMeasure::Cauchy),
        SINGLE,
    )
    .unwrap();
    assert!(
        bind(other_curve).is_err(),
        "equal numerical curves on a different card cannot replace the source"
    );
    let mut finite_properties = properties();
    finite_properties.push(curve(YieldStressMeasure::Kirchhoff));
    let finite_source = card(finite_properties);
    let finite_material = resolve_material_state_point(
        &finite_source,
        &context(),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let kirchhoff = resolve_tabulated_hardening(
        &finite_source,
        &context(),
        &request(YieldStressMeasure::Kirchhoff),
        SINGLE,
    )
    .unwrap();
    assert_eq!(
        kirchhoff.source().claim().key.name(),
        KIRCHHOFF_YIELD_CURVE_PROPERTY
    );
    assert!(
        ResolvedTabulatedJ2::try_from_material(
            finite_material,
            kirchhoff,
            TabulatedKinematicHardening::None,
            0.05
        )
        .is_err()
    );
}
