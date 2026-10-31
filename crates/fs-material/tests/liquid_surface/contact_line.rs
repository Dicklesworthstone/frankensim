//! The same original interface bundle must serve transport and contact laws.
use super::*;
use fs_matdb::{InterfaceSystemCard, SurfaceSpec, SystemContext};
use fs_material::liquid::surface::adhesion::{
    EQUILIBRIUM_WORK_OF_ADHESION_PROPERTY, EquilibriumAngleOrigin,
};
use fs_material::liquid::surface::contact_line::*;
use fs_material::liquid::surface::{
    EQUILIBRIUM_CONTACT_ANGLE_PROPERTY, ResolvedLiquidWettingState, resolve_liquid_wetting_state,
};
use fs_material::state_point::resolve_interface_state_point;

fn liquid_at(query: &QueryPoint) -> (MaterialCard, ResolvedLiquidSurfaceState) {
    let source = card(
        1_000.0,
        0.1,
        LiquidViscosityBasis::Dynamic,
        Some((0.05, SURFACE_TENSION_DIMS)),
    );
    let liquid = resolve_liquid_surface_state(
        &source,
        query,
        LiquidViscosityBasis::Dynamic,
        LiquidCapillarity::Sourced,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    (source, liquid)
}

fn interface(
    liquid: &MaterialCard,
    origin: EquilibriumAngleOrigin,
    modify: impl Fn(&mut PropertyClaim),
) -> InterfaceSystemCard {
    let equilibrium = match origin {
        EquilibriumAngleOrigin::SourcedAngle => {
            (EQUILIBRIUM_CONTACT_ANGLE_PROPERTY, 1.0, Dims::NONE)
        }
        EquilibriumAngleOrigin::YoungDupre => (
            EQUILIBRIUM_WORK_OF_ADHESION_PROPERTY,
            0.05 * (1.0 + 1.0_f64.cos()),
            SURFACE_TENSION_DIMS,
        ),
    };
    let mut claims = ClaimSet::new();
    for (name, value, dims) in [
        equilibrium,
        (CONTACT_LINE_FRICTION_PROPERTY, 5.0, DYNAMIC_VISCOSITY_DIMS),
        (ADVANCING_CONTACT_ANGLE_PROPERTY, 1.2, Dims::NONE),
        (RECEDING_CONTACT_ANGLE_PROPERTY, 0.8, Dims::NONE),
    ] {
        let mut claim = PropertyClaim {
            key: PropertyKey::new(name, dims),
            value: PropertyValue::Scalar { value, dims },
            validity: ValidityDomain::unconstrained()
                .with_quantity(
                    "temperature",
                    QuantitySpec::semantic(SemanticType::new(
                        QuantityKind::AbsoluteTemperature,
                        ValueForm::Static,
                    )),
                    290.0,
                    310.0,
                )
                .with("surface_environment", 7.0, 7.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![],
            provenance: Provenance {
                source: "synthetic contact-line source".into(),
                license: "test-only".into(),
                artifact: None,
            },
        };
        modify(&mut claim);
        claims.insert_claim(claim).unwrap();
    }
    InterfaceSystemCard::assemble(
        SurfaceSpec {
            material: MaterialStateId {
                chemistry: "substrate".into(),
                phase: "solid".into(),
                process: "prepared fixture".into(),
                revision: 0,
            },
            texture_frame: "declared finish".into(),
        },
        SurfaceSpec {
            material: liquid.id().clone(),
            texture_frame: "liquid".into(),
        },
        SystemContext {
            medium: "declared vapor".into(),
            third_body: None,
            environment: "synthetic isothermal source".into(),
            history: "fixed clean preparation".into(),
        },
        claims,
        vec![],
    )
    .unwrap()
}

#[test]
fn both_angle_origins_keep_the_original_complete_source_bundle() {
    let (card, liquid) = liquid_at(&point());
    for origin in [
        EquilibriumAngleOrigin::SourcedAngle,
        EquilibriumAngleOrigin::YoungDupre,
    ] {
        let source = interface(&card, origin, |_| {});
        let wetting = resolve_contact_line_wetting(
            &source,
            &liquid,
            origin,
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
        let line = FrozenContactLineWetting::new(&liquid, &wetting).unwrap();
        assert!(core::ptr::eq(line.wetting(), &wetting));
        assert_eq!(line.friction_pa_s(), 5.0);
        assert_ne!(line.friction_pa_s(), liquid.dynamic_viscosity_pa_s());
        assert!((line.equilibrium_angle_radians() - 1.0).abs() < 3e-15);
        assert_eq!(line.advancing_angle_radians(), 1.2);
        assert_eq!(line.receding_angle_radians(), 0.8);
        assert_eq!(wetting.angle_origin(), origin);
        assert_eq!(wetting.context().history, "fixed clean preparation");
        assert_eq!(wetting.resolved().card_identity(), source.content_hash());
        assert_eq!(wetting.resolved().properties().len(), 4);
        for property in wetting.resolved().properties() {
            assert_eq!(
                property.answer().receipt.query_point,
                liquid.material().query_point()
            );
            assert_eq!(
                property.answer().receipt.axis_quantities,
                *point().axis_quantities()
            );
        }
        assert_eq!(
            wetting
                .resolved()
                .property(EQUILIBRIUM_CONTACT_ANGLE_PROPERTY)
                .is_some(),
            origin == EquilibriumAngleOrigin::SourcedAngle
        );
    }
}

#[test]
fn missing_coefficients_and_wrong_liquid_point_never_fall_back() {
    let (card, liquid) = liquid_at(&point());
    let source = interface(&card, EquilibriumAngleOrigin::SourcedAngle, |_| {});
    let incomplete = resolve_liquid_wetting_state(
        &source,
        &point(),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    assert!(FrozenContactLineWetting::new(&liquid, &incomplete).is_err());
    let wetting = resolve_contact_line_wetting(
        &source,
        &liquid,
        EquilibriumAngleOrigin::SourcedAngle,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let (_, different_point) = liquid_at(&point().with("extra_context", 1.0).unwrap());
    assert!(FrozenContactLineWetting::new(&different_point, &wetting).is_err());
    assert!(FrozenContactLineWetting::new(&liquid, &wetting).is_ok());
}

#[test]
fn frozen_contact_line_also_requires_a_history_independent_equilibrium_potential() {
    let quantity = QuantitySpec::dimensional(fs_qty::Time::DIMS);
    let query = point().with_quantity("exposure", quantity, 1.0).unwrap();
    let (card, liquid) = liquid_at(&query);
    for origin in [
        EquilibriumAngleOrigin::SourcedAngle,
        EquilibriumAngleOrigin::YoungDupre,
    ] {
        let source = interface(&card, origin, |claim| {
            if [
                EQUILIBRIUM_CONTACT_ANGLE_PROPERTY,
                EQUILIBRIUM_WORK_OF_ADHESION_PROPERTY,
            ]
            .contains(&claim.key.name())
            {
                claim.validity = claim
                    .validity
                    .clone()
                    .with_quantity("exposure", quantity, 0.0, 2.0);
            }
        });
        let wetting = resolve_contact_line_wetting(
            &source,
            &liquid,
            origin,
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
        assert!(FrozenContactLineWetting::new(&liquid, &wetting).is_err());
    }
}

#[test]
fn wrongly_typed_friction_and_unordered_hysteresis_are_refused() {
    let (card, liquid) = liquid_at(&point());
    let wrong = interface(&card, EquilibriumAngleOrigin::SourcedAngle, |claim| {
        if claim.key.name() == CONTACT_LINE_FRICTION_PROPERTY {
            claim.key = PropertyKey::new(CONTACT_LINE_FRICTION_PROPERTY, SURFACE_TENSION_DIMS);
            claim.value = PropertyValue::Scalar {
                value: 5.0,
                dims: SURFACE_TENSION_DIMS,
            };
        }
    });
    let requirements: Vec<_> = [
        (EQUILIBRIUM_CONTACT_ANGLE_PROPERTY, Dims::NONE),
        (CONTACT_LINE_FRICTION_PROPERTY, SURFACE_TENSION_DIMS),
        (ADVANCING_CONTACT_ANGLE_PROPERTY, Dims::NONE),
        (RECEDING_CONTACT_ANGLE_PROPERTY, Dims::NONE),
    ]
    .into_iter()
    .map(|(name, dims)| {
        ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::StrictlyPositive)
            .unwrap()
    })
    .collect();
    let bundle = resolve_interface_state_point(
        &wrong,
        &point(),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let wetting = ResolvedLiquidWettingState::try_from_resolved(&wrong, bundle).unwrap();
    assert!(FrozenContactLineWetting::new(&liquid, &wetting).is_err());
    for (name, value) in [
        (RECEDING_CONTACT_ANGLE_PROPERTY, 1.1),
        (ADVANCING_CONTACT_ANGLE_PROPERTY, 0.9),
        (ADVANCING_CONTACT_ANGLE_PROPERTY, core::f64::consts::PI),
    ] {
        let source = interface(&card, EquilibriumAngleOrigin::SourcedAngle, |claim| {
            if claim.key.name() == name {
                claim.value = PropertyValue::Scalar {
                    value,
                    dims: Dims::NONE,
                };
            }
        });
        let wetting = resolve_contact_line_wetting(
            &source,
            &liquid,
            EquilibriumAngleOrigin::SourcedAngle,
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
        assert!(
            FrozenContactLineWetting::new(&liquid, &wetting).is_err(),
            "{name}={value}"
        );
    }
}

#[test]
fn temperature_support_is_frozen_but_speed_rate_and_history_support_are_not() {
    let temperature = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ));
    for (axis, quantity, low, at, high, admitted) in [
        ("temperature", temperature, 290.0, 300.0, 310.0, true),
        (
            "speed",
            QuantitySpec::dimensional(Dims([1, 0, -1, 0, 0, 0])),
            0.0,
            1.0,
            2.0,
            false,
        ),
        (
            "rate",
            QuantitySpec::dimensional(fs_qty::Frequency::DIMS),
            0.0,
            1.0,
            2.0,
            false,
        ),
        (
            "exposure",
            QuantitySpec::dimensional(fs_qty::Time::DIMS),
            0.0,
            1.0,
            2.0,
            false,
        ),
        (
            "cycle_count",
            QuantitySpec::dimensional(Dims::NONE),
            0.0,
            1.0,
            2.0,
            false,
        ),
    ] {
        for as_curve in [false, true] {
            let query = point().with_quantity(axis, quantity, at).unwrap();
            let (card, liquid) = liquid_at(&query);
            let source = interface(&card, EquilibriumAngleOrigin::SourcedAngle, |claim| {
                if claim.key.name() == CONTACT_LINE_FRICTION_PROPERTY {
                    claim.validity = claim
                        .validity
                        .clone()
                        .with_quantity(axis, quantity, low, high);
                    if as_curve {
                        claim.value = PropertyValue::Curve {
                            abscissa: axis.into(),
                            abscissa_dims: quantity.dims(),
                            knots: vec![(low, 4.0), (high, 6.0)],
                            dims: DYNAMIC_VISCOSITY_DIMS,
                        };
                        claim.interpolation = InterpolationPolicy::LinearInside;
                    }
                }
            });
            let wetting = resolve_contact_line_wetting(
                &source,
                &liquid,
                EquilibriumAngleOrigin::SourcedAngle,
                MaterialPropertySelection::SingleClaimOnly,
            )
            .unwrap();
            let result = FrozenContactLineWetting::new(&liquid, &wetting);
            assert_eq!(result.is_ok(), admitted, "axis={axis}, curve={as_curve}");
            if let Ok(line) = result {
                assert_eq!(line.friction_pa_s(), 5.0);
            }
        }
    }
}
