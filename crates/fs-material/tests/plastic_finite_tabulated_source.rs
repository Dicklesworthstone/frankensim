//! Original Kirchhoff curve/model bindings drive native finite flow.
//! Synthetic source cards verify transport and mechanics, not alloy calibration.

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSelection, ClaimSet, ConstitutiveModelCard, InitialStatePolicy, InterpolationPolicy,
    LawId, LawParameter, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey, PropertyValue,
    Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::plastic::finite::FiniteJ2Material;
use fs_material::plastic::finite::tabulated::source::{
    FINITE_TABULATED_J2_LAW_ID, ResolvedTabulatedFiniteJ2,
};
use fs_material::plastic::tabulated::source::{
    ACCUMULATED_PLASTIC_STRAIN_AXIS, CAUCHY_YIELD_CURVE_PROPERTY, KIRCHHOFF_YIELD_CURVE_PROPERTY,
    YieldStressMeasure,
};
use fs_qty::{Dims, Pressure, QuantitySpec};

const SINGLE: ClaimSelection = ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly);
const MU: f64 = 80e6;
const K: f64 = 160e6;

fn provenance() -> Provenance {
    Provenance {
        source: "synthetic bounded Kirchhoff plastic-flow fixture".into(),
        license: "test-only".into(),
        artifact: None,
    }
}

fn model() -> ConstitutiveModelCard {
    ConstitutiveModelCard {
        law: LawId(FINITE_TABULATED_J2_LAW_ID.into()),
        law_version: 1,
        parameters: [
            ("shear_modulus_pa", MU, Pressure::DIMS),
            ("bulk_modulus_pa", K, Pressure::DIMS),
            ("minimum_jacobian", 0.2, Dims::NONE),
            ("maximum_elastic_log_strain", 0.6, Dims::NONE),
            ("maximum_plastic_log_strain", 2.0, Dims::NONE),
            ("maximum_accumulated_plastic_strain", 0.2, Dims::NONE),
        ]
        .into_iter()
        .map(|(name, value, dims)| (name.into(), LawParameter { value, dims }))
        .collect(),
        state_schema_version: 1,
        initial_state: InitialStatePolicy::ZeroInternalState,
        validity: ValidityDomain::unconstrained().with("T", 300.0, 400.0),
        sources: vec![],
        provenance: provenance(),
    }
}

fn card(model: ConstitutiveModelCard, property: &str) -> MaterialCard {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(property, Pressure::DIMS),
            value: PropertyValue::Curve {
                abscissa: ACCUMULATED_PLASTIC_STRAIN_AXIS.into(),
                abscissa_dims: Dims::NONE,
                knots: vec![(0.0, 2e6), (0.05, 4e6), (0.2, 6e6)],
                dims: Pressure::DIMS,
            },
            validity: ValidityDomain::unconstrained()
                .with("T", 300.0, 400.0)
                .with_quantity(
                    ACCUMULATED_PLASTIC_STRAIN_AXIS,
                    QuantitySpec::dimensional(Dims::NONE),
                    0.0,
                    0.2,
                ),
            interpolation: InterpolationPolicy::LinearInside,
            uncertainty: UncertaintyModel::Unstated,
            observations: vec![],
            provenance: provenance(),
        })
        .unwrap();
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic alloy".into(),
            phase: "solid".into(),
            process: "bounded curve fixture".into(),
            revision: 0,
        },
        claims,
        vec![model],
    )
    .unwrap()
}

fn context() -> QueryPoint {
    QueryPoint::new().with("T", 350.0).unwrap()
}
fn close(a: f64, b: f64) {
    assert!(
        (a - b).abs() <= 3e-10 * b.abs().max(1.0),
        "{a:.16e} != {b:.16e}"
    );
}

#[test]
fn selected_original_curve_drives_kirchhoff_flow_at_nonunit_volume() {
    let card = card(model(), KIRCHHOFF_YIELD_CURVE_PROPERTY);
    let q = context();
    let source =
        ResolvedTabulatedFiniteJ2::try_from_material_card_at(&card, 0, &q, SINGLE, 3).unwrap();
    assert_eq!(source.model_card(), &card.models()[0]);
    assert_eq!(source.model_index(), 0);
    assert_eq!(source.context(), &q);
    assert_eq!(
        source.hardening_source().source().card_identity(),
        card.content_hash()
    );
    assert_eq!(source.hardening_source().source().context(), &q);
    assert_eq!(
        source.hardening_source().stress_measure(),
        YieldStressMeasure::Kirchhoff
    );
    assert_eq!(source.hardening_source().source().observations().len(), 3);
    assert!(std::ptr::eq(source.law(), source.finite_j2()));

    // Prescribe a plastic endpoint on the second supplied affine segment.
    // The scalar Kirchhoff relation independently gives total deviatoric log strain.
    let (alpha, strength, mean) = (0.08, 4.4e6, 0.03_f64);
    let axial: f64 = alpha + strength / (3.0 * MU);
    let f = [
        (mean + axial).exp(),
        0.0,
        0.0,
        0.0,
        (mean - 0.5 * axial).exp(),
        0.0,
        0.0,
        0.0,
        (mean - 0.5 * axial).exp(),
    ];
    let response = source
        .finite_j2()
        .response(&f, &source.law().initial_state())
        .unwrap();
    close(response.state.alpha, alpha);
    close(response.equivalent_kirchhoff_pa, strength);
    close(response.jacobian, (3.0 * mean).exp());
    assert!(response.equivalent_kirchhoff_pa / response.jacobian < strength);
    let stored_hardening = 0.5 * 2e6 * 0.05 + 0.5 * (2e6 + 2.4e6) * 0.03;
    close(
        response.stored_energy_j_m3,
        0.5 * K * (3.0 * mean).powi(2) + strength * strength / (6.0 * MU) + stored_hardening,
    );
    close(response.plastic_dissipation_j_m3, 2e6 * alpha);
    source
        .law()
        .validate_committed(&f, &response.state, 1e-3)
        .unwrap();
}

#[test]
fn source_stress_convention_support_and_fixed_context_cannot_be_substituted() {
    let cauchy = card(model(), CAUCHY_YIELD_CURVE_PROPERTY);
    assert!(
        ResolvedTabulatedFiniteJ2::try_from_material_card_at(&cauchy, 0, &context(), SINGLE, 3)
            .is_err()
    );
    let original = card(model(), KIRCHHOFF_YIELD_CURVE_PROPERTY);
    assert!(
        ResolvedTabulatedFiniteJ2::try_from_material_card_at(&original, 1, &context(), SINGLE, 3)
            .is_err()
    );
    assert!(
        ResolvedTabulatedFiniteJ2::try_from_material_card_at(&original, 0, &context(), SINGLE, 2)
            .is_err()
    );
    assert!(
        ResolvedTabulatedFiniteJ2::try_from_material_card_at(
            &original,
            0,
            &QueryPoint::new(),
            SINGLE,
            3
        )
        .is_err()
    );
    let hot = QueryPoint::new().with("T", 401.0).unwrap();
    assert!(
        ResolvedTabulatedFiniteJ2::try_from_material_card_at(&original, 0, &hot, SINGLE, 3)
            .is_err()
    );
    let frozen_history = context()
        .with(ACCUMULATED_PLASTIC_STRAIN_AXIS, 0.0)
        .unwrap();
    assert!(
        ResolvedTabulatedFiniteJ2::try_from_material_card_at(
            &original,
            0,
            &frozen_history,
            SINGLE,
            3
        )
        .is_err()
    );

    let mut too_long = model();
    too_long
        .parameters
        .get_mut("maximum_accumulated_plastic_strain")
        .unwrap()
        .value = 0.21;
    let unsupported = card(too_long, KIRCHHOFF_YIELD_CURVE_PROPERTY);
    assert!(
        ResolvedTabulatedFiniteJ2::try_from_material_card_at(
            &unsupported,
            0,
            &context(),
            SINGLE,
            3
        )
        .is_err()
    );
}

#[test]
fn table_model_requires_its_original_typed_coefficients_and_explicit_selector() {
    for modification in 0..5 {
        let mut wrong = model();
        match modification {
            0 => {
                wrong.parameters.remove("bulk_modulus_pa");
            }
            1 => {
                wrong.parameters.get_mut("shear_modulus_pa").unwrap().dims = Dims::NONE;
            }
            2 => {
                wrong.parameters.insert(
                    "hardening_modulus_pa".into(),
                    LawParameter {
                        value: 0.0,
                        dims: Pressure::DIMS,
                    },
                );
            }
            3 => {
                wrong.law = LawId("finite-j2-hencky-kirchhoff".into());
            }
            _ => {
                wrong.initial_state = InitialStatePolicy::RequiresDeclaredState;
            }
        }
        let card = card(wrong, KIRCHHOFF_YIELD_CURVE_PROPERTY);
        assert!(
            ResolvedTabulatedFiniteJ2::try_from_material_card_at(&card, 0, &context(), SINGLE, 3)
                .is_err(),
            "case {modification}"
        );
    }
}
