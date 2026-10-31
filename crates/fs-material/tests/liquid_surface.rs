//! Source admission for isothermal momentum/capillarity. These are synthetic
//! constitutive fixtures, not measured liquids or a capillary accuracy oracle.

#[path = "liquid_surface/contact_line.rs"]
mod contact_line;

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::liquid::surface::{
    LiquidCapillarity, ResolvedLiquidSurfaceState, SURFACE_TENSION_DIMS, SURFACE_TENSION_PROPERTY,
    resolve_liquid_surface_state,
};
use fs_material::liquid::{
    DYNAMIC_VISCOSITY_DIMS, DYNAMIC_VISCOSITY_PROPERTY, KINEMATIC_VISCOSITY_DIMS,
    KINEMATIC_VISCOSITY_PROPERTY,
};
use fs_material::liquid_temperature::LiquidViscosityBasis;
use fs_material::state_point::{
    DENSITY_PROPERTY, MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
    resolve_material_state_point,
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Density, Dims, QuantitySpec};

fn point() -> QueryPoint {
    QueryPoint::new()
        .with_quantity(
            "temperature",
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::AbsoluteTemperature,
                ValueForm::Static,
            )),
            300.0,
        )
        .unwrap()
        .with("pressure", 100_000.0)
        .unwrap()
        .with("surface_environment", 7.0)
        .unwrap()
}

fn card(
    density: f64,
    viscosity: f64,
    basis: LiquidViscosityBasis,
    surface: Option<(f64, Dims)>,
) -> MaterialCard {
    let (viscosity_name, viscosity_dims) = match basis {
        LiquidViscosityBasis::Dynamic => (DYNAMIC_VISCOSITY_PROPERTY, DYNAMIC_VISCOSITY_DIMS),
        LiquidViscosityBasis::Kinematic => (KINEMATIC_VISCOSITY_PROPERTY, KINEMATIC_VISCOSITY_DIMS),
    };
    let mut rows = vec![
        (DENSITY_PROPERTY, Density::DIMS, density),
        (viscosity_name, viscosity_dims, viscosity),
    ];
    if let Some((sigma, dims)) = surface {
        rows.push((SURFACE_TENSION_PROPERTY, dims, sigma));
    }
    let mut claims = ClaimSet::new();
    for (name, dims, value) in rows {
        let mut validity = ValidityDomain::unconstrained()
            .with_quantity(
                "temperature",
                QuantitySpec::semantic(SemanticType::new(
                    QuantityKind::AbsoluteTemperature,
                    ValueForm::Static,
                )),
                290.0,
                310.0,
            )
            .with("pressure", 100_000.0, 100_000.0);
        if name == SURFACE_TENSION_PROPERTY {
            validity = validity.with("surface_environment", 7.0, 7.0);
        }
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity,
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic liquid surface fixture".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic liquid".into(),
            phase: "liquid".into(),
            process: "numerical fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}

#[test]
fn isothermal_surface_resolves_only_used_claims_and_retains_full_context() {
    let card = card(
        1_000.0,
        0.1,
        LiquidViscosityBasis::Dynamic,
        Some((0.02, SURFACE_TENSION_DIMS)),
    );
    let state = resolve_liquid_surface_state(
        &card,
        &point(),
        LiquidViscosityBasis::Dynamic,
        LiquidCapillarity::Sourced,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    assert_eq!(state.material().properties().len(), 3);
    assert_eq!(state.density_kg_m3(), 1_000.0);
    assert_eq!(state.dynamic_viscosity_pa_s(), 0.1);
    assert_eq!(state.kinematic_viscosity_m2_s(), 0.0001);
    assert_eq!(state.surface_tension_n_m(), Some(0.02));
    assert_eq!(state.capillarity(), LiquidCapillarity::Sourced);
    assert_eq!(state.material().card_identity(), card.content_hash());
    assert!(
        state
            .material()
            .property("specific_heat_capacity")
            .is_none()
    );
    assert!(state.material().property("thermal_conductivity").is_none());
    for property in state.material().properties() {
        assert_eq!(
            property.answer().receipt.axis_quantities,
            *point().axis_quantities()
        );
        assert!(
            property
                .answer()
                .receipt
                .query_point
                .contains(&("surface_environment".into(), 7.0))
        );
    }
    let pins = state
        .material()
        .properties()
        .iter()
        .map(|property| {
            (
                property.requirement().name().to_owned(),
                property.answer().receipt.selected,
            )
        })
        .collect();
    let pinned = resolve_liquid_surface_state(
        &card,
        &point(),
        LiquidViscosityBasis::Dynamic,
        LiquidCapillarity::Sourced,
        MaterialPropertySelection::PinnedByProperty(pins),
    )
    .unwrap();
    assert_eq!(pinned.surface_tension_n_m(), state.surface_tension_n_m());
}

#[test]
fn kinematic_source_uses_same_point_density_without_creating_a_mu_claim() {
    let card = card(
        800.0,
        2e-6,
        LiquidViscosityBasis::Kinematic,
        Some((0.03, SURFACE_TENSION_DIMS)),
    );
    let state = resolve_liquid_surface_state(
        &card,
        &point(),
        LiquidViscosityBasis::Kinematic,
        LiquidCapillarity::Sourced,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    assert_eq!(state.kinematic_viscosity_m2_s(), 2e-6);
    assert!((state.dynamic_viscosity_pa_s() - 0.0016).abs() < 1e-18);
    assert_eq!(state.viscosity_basis(), LiquidViscosityBasis::Kinematic);
    assert!(
        state
            .material()
            .property(KINEMATIC_VISCOSITY_PROPERTY)
            .is_some()
    );
    assert!(
        state
            .material()
            .property(DYNAMIC_VISCOSITY_PROPERTY)
            .is_none()
    );
    assert!(
        resolve_liquid_surface_state(
            &card,
            &point(),
            LiquidViscosityBasis::Dynamic,
            LiquidCapillarity::Sourced,
            MaterialPropertySelection::SingleClaimOnly
        )
        .is_err()
    );
}

#[test]
fn shear_rate_source_retains_evolution_obligation_and_original_claims() {
    for basis in [
        LiquidViscosityBasis::Dynamic,
        LiquidViscosityBasis::Kinematic,
    ] {
        let original = card(1000.0, 0.1, basis, None);
        let (viscosity, dims, scale) = match basis {
            LiquidViscosityBasis::Dynamic => {
                (DYNAMIC_VISCOSITY_PROPERTY, DYNAMIC_VISCOSITY_DIMS, 1.0)
            }
            LiquidViscosityBasis::Kinematic => (
                KINEMATIC_VISCOSITY_PROPERTY,
                KINEMATIC_VISCOSITY_DIMS,
                0.001,
            ),
        };
        let mut claims = ClaimSet::new();
        for (_, original) in original.claims().claims_ordered() {
            let mut claim = original.clone();
            if claim.key.name() == viscosity {
                claim.value = PropertyValue::Curve {
                    abscissa: "rate".into(),
                    abscissa_dims: fs_qty::Frequency::DIMS,
                    knots: vec![(0.0, 0.1 * scale), (100.0, 0.3 * scale)],
                    dims,
                };
                claim.interpolation = InterpolationPolicy::LinearInside;
                claim.validity = claim.validity.with_quantity(
                    "rate",
                    QuantitySpec::dimensional(fs_qty::Frequency::DIMS),
                    0.0,
                    100.0,
                );
            }
            claims.insert_claim(claim).unwrap();
        }
        let source = MaterialCard::assemble(original.id().clone(), claims, vec![]).unwrap();
        let point = point()
            .with_quantity(
                "rate",
                QuantitySpec::dimensional(fs_qty::Frequency::DIMS),
                50.0,
            )
            .unwrap();
        let requirements: Vec<_> = [(DENSITY_PROPERTY, Density::DIMS), (viscosity, dims)]
            .into_iter()
            .map(|(name, dims)| {
                ScalarPropertyRequirement::try_new(
                    name,
                    dims,
                    ScalarAdmissibility::StrictlyPositive,
                )
                .unwrap()
            })
            .collect();
        let bundle = resolve_material_state_point(
            &source,
            &point,
            &requirements,
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
        assert!(
            ResolvedLiquidSurfaceState::try_from_material(
                bundle.clone(),
                basis,
                LiquidCapillarity::Neglected
            )
            .is_err()
        );
        let evolving = ResolvedLiquidSurfaceState::try_for_shear_rate_evolution(
            bundle.clone(),
            basis,
            "temperature",
            "rate",
        )
        .unwrap();
        assert_eq!(evolving.rheology_axes(), Some(("temperature", "rate")));
        assert_eq!(evolving.material(), &bundle);
        assert_eq!(evolving.capillarity(), LiquidCapillarity::Neglected);
        assert!((evolving.dynamic_viscosity_pa_s() - 0.2).abs() < 1e-15);
        assert!((evolving.kinematic_viscosity_m2_s() - 0.0002).abs() < 1e-18);
        assert!(
            ResolvedLiquidSurfaceState::try_for_shear_rate_evolution(
                bundle.clone(),
                basis,
                "temperature",
                "absent",
            )
            .is_err()
        );
        assert!(
            ResolvedLiquidSurfaceState::try_for_shear_rate_evolution(
                bundle,
                basis,
                "temperature",
                "temperature",
            )
            .is_err()
        );
    }
    let source = card(1000.0, 0.1, LiquidViscosityBasis::Dynamic, None);
    let wrong = point()
        .with_quantity(
            "rate",
            QuantitySpec::dimensional(fs_qty::Temperature::DIMS),
            1.0,
        )
        .unwrap();
    let bundle = resolve_liquid_surface_state(
        &source,
        &wrong,
        LiquidViscosityBasis::Dynamic,
        LiquidCapillarity::Neglected,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    assert!(
        ResolvedLiquidSurfaceState::try_for_shear_rate_evolution(
            bundle.material().clone(),
            LiquidViscosityBasis::Dynamic,
            "temperature",
            "rate",
        )
        .is_err()
    );
}

#[test]
fn omitting_capillary_stress_is_explicit_and_does_not_require_sigma() {
    let card = card(1_000.0, 0.1, LiquidViscosityBasis::Dynamic, None);
    let state = resolve_liquid_surface_state(
        &card,
        &point(),
        LiquidViscosityBasis::Dynamic,
        LiquidCapillarity::Neglected,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    assert_eq!(state.material().properties().len(), 2);
    assert_eq!(state.surface_tension_n_m(), None);
    assert_eq!(state.capillarity(), LiquidCapillarity::Neglected);
    assert!(
        resolve_liquid_surface_state(
            &card,
            &point(),
            LiquidViscosityBasis::Dynamic,
            LiquidCapillarity::Sourced,
            MaterialPropertySelection::SingleClaimOnly
        )
        .is_err()
    );
}

#[test]
fn surface_context_units_and_positive_domain_are_not_silently_replaced() {
    let source = card(
        1_000.0,
        0.1,
        LiquidViscosityBasis::Dynamic,
        Some((0.02, SURFACE_TENSION_DIMS)),
    );
    for query in [
        point().with("temperature", 311.0).unwrap(),
        point().with("surface_environment", 8.0).unwrap(),
        point().with("pressure", 200_000.0).unwrap(),
    ] {
        assert!(
            resolve_liquid_surface_state(
                &source,
                &query,
                LiquidViscosityBasis::Dynamic,
                LiquidCapillarity::Sourced,
                MaterialPropertySelection::SingleClaimOnly
            )
            .is_err()
        );
    }
    for surface in [
        (0.0, SURFACE_TENSION_DIMS),
        (-0.01, SURFACE_TENSION_DIMS),
        (0.02, Dims([1, 1, -2, 0, 0, 0])),
    ] {
        let wrong = card(1_000.0, 0.1, LiquidViscosityBasis::Dynamic, Some(surface));
        assert!(
            resolve_liquid_surface_state(
                &wrong,
                &point(),
                LiquidViscosityBasis::Dynamic,
                LiquidCapillarity::Sourced,
                MaterialPropertySelection::SingleClaimOnly
            )
            .is_err()
        );
    }
    for (rho, viscosity, basis) in [
        (1e300, 1e100, LiquidViscosityBasis::Kinematic),
        (1e100, 1e-300, LiquidViscosityBasis::Dynamic),
    ] {
        let wrong = card(rho, viscosity, basis, None);
        assert!(
            resolve_liquid_surface_state(
                &wrong,
                &point(),
                basis,
                LiquidCapillarity::Neglected,
                MaterialPropertySelection::SingleClaimOnly
            )
            .is_err()
        );
    }
}

#[test]
fn shared_bundle_is_adopted_without_requery_or_loss_of_other_physics() {
    let original = card(
        1_000.0,
        0.1,
        LiquidViscosityBasis::Dynamic,
        Some((0.02, SURFACE_TENSION_DIMS)),
    );
    let mut claims = original.claims().clone();
    for (name, dims, value) in [
        (KINEMATIC_VISCOSITY_PROPERTY, KINEMATIC_VISCOSITY_DIMS, 2e-4),
        (
            "specific_heat_capacity",
            Dims([2, 0, -2, -1, 0, 0]),
            4_000.0,
        ),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained(),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic multiphysics observation".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let shared_card = MaterialCard::assemble(original.id().clone(), claims, vec![]).unwrap();
    let requirements: Vec<_> = [
        (DENSITY_PROPERTY, Density::DIMS),
        (DYNAMIC_VISCOSITY_PROPERTY, DYNAMIC_VISCOSITY_DIMS),
        (KINEMATIC_VISCOSITY_PROPERTY, KINEMATIC_VISCOSITY_DIMS),
        (SURFACE_TENSION_PROPERTY, SURFACE_TENSION_DIMS),
        ("specific_heat_capacity", Dims([2, 0, -2, -1, 0, 0])),
    ]
    .into_iter()
    .map(|(name, dims)| {
        ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::Finite).unwrap()
    })
    .collect();
    let bundle = resolve_material_state_point(
        &shared_card,
        &point(),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let dynamic = ResolvedLiquidSurfaceState::try_from_material(
        bundle.clone(),
        LiquidViscosityBasis::Dynamic,
        LiquidCapillarity::Sourced,
    )
    .unwrap();
    let kinematic = ResolvedLiquidSurfaceState::try_from_material(
        bundle.clone(),
        LiquidViscosityBasis::Kinematic,
        LiquidCapillarity::Sourced,
    )
    .unwrap();
    assert_eq!(dynamic.material(), &bundle);
    assert_eq!(kinematic.material(), &bundle);
    assert_eq!(dynamic.dynamic_viscosity_pa_s(), 0.1);
    assert_eq!(kinematic.dynamic_viscosity_pa_s(), 0.2);
    assert_eq!(
        kinematic
            .material()
            .property("specific_heat_capacity")
            .unwrap()
            .value_si(),
        4_000.0
    );
    // The explicitly independent nu source does not rewrite the other mu observation.
    assert_eq!(
        kinematic
            .material()
            .property(DYNAMIC_VISCOSITY_PROPERTY)
            .unwrap()
            .value_si(),
        0.1
    );

    let nonpositive = card(0.0, 0.1, LiquidViscosityBasis::Dynamic, None);
    let broad = resolve_material_state_point(
        &nonpositive,
        &point(),
        &requirements[..2],
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    assert!(
        ResolvedLiquidSurfaceState::try_from_material(
            broad,
            LiquidViscosityBasis::Dynamic,
            LiquidCapillarity::Neglected
        )
        .is_err()
    );
}

#[test]
fn typed_temperature_aliases_are_static_but_rate_and_frequency_response_cannot_be_frozen() {
    let original = card(
        1_000.0,
        0.1,
        LiquidViscosityBasis::Dynamic,
        Some((0.02, SURFACE_TENSION_DIMS)),
    );
    let absolute = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ));
    for property_name in [DYNAMIC_VISCOSITY_PROPERTY, SURFACE_TENSION_PROPERTY] {
        for (axis, quantity, lo, at, hi, admitted) in [
            ("temperature", absolute, 290.0, 300.0, 310.0, true),
            (
                "T",
                QuantitySpec::dimensional(fs_qty::Temperature::DIMS),
                290.0,
                300.0,
                310.0,
                true,
            ),
            (
                "renamed",
                QuantitySpec::dimensional(fs_qty::Temperature::DIMS),
                290.0,
                300.0,
                310.0,
                false,
            ),
            (
                "omega",
                QuantitySpec::dimensional(fs_qty::Frequency::DIMS),
                1.0,
                5.0,
                10.0,
                false,
            ),
            (
                "strain_rate",
                QuantitySpec::dimensional(fs_qty::Frequency::DIMS),
                1.0,
                5.0,
                10.0,
                false,
            ),
            (
                "strain",
                QuantitySpec::dimensional(Dims::NONE),
                0.0,
                0.05,
                0.1,
                false,
            ),
        ] {
            // Constant-valued data restricted to a dynamic band also cannot
            // become an unrestricted Newtonian/surface law during evolution.
            for as_curve in [true, false] {
                let mut claims = ClaimSet::new();
                for (_, source) in original.claims().claims_ordered() {
                    let mut claim = source.clone();
                    if claim.key.name() == property_name {
                        claim.validity = claim.validity.with_quantity(axis, quantity, lo, hi);
                        if as_curve {
                            let PropertyValue::Scalar { value, dims } = claim.value else {
                                unreachable!()
                            };
                            claim.value = PropertyValue::Curve {
                                abscissa: axis.into(),
                                abscissa_dims: quantity.dims(),
                                knots: vec![(lo, value), (hi, 2.0 * value)],
                                dims,
                            };
                            claim.interpolation = InterpolationPolicy::LinearInside;
                        }
                    }
                    claims.insert_claim(claim).unwrap();
                }
                let card = MaterialCard::assemble(original.id().clone(), claims, vec![]).unwrap();
                let query = point().with_quantity(axis, quantity, at).unwrap();
                let result = resolve_liquid_surface_state(
                    &card,
                    &query,
                    LiquidViscosityBasis::Dynamic,
                    LiquidCapillarity::Sourced,
                    MaterialPropertySelection::SingleClaimOnly,
                );
                let expected = admitted || (!as_curve && axis == "renamed");
                assert_eq!(
                    result.is_ok(),
                    expected,
                    "property={property_name}, axis={axis}, curve={as_curve}"
                );
            }
        }
    }
}
