//! Original junction claims follow face T while the normal owner stays frozen.
use super::*;
use core::f64::consts::PI;
use fs_contact::interface_binding::adhesive_shear::spinning::temperature::SourceThermalAdhesivePatchError;
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{QuantitySpec, Temperature};

fn curve(name: &str, axis: &str, range: [f64; 2], values: [f64; 2]) -> PropertyClaim {
    PropertyClaim {
        key: PropertyKey::new(name, Pressure::DIMS),
        value: PropertyValue::Curve {
            abscissa: axis.into(),
            abscissa_dims: Temperature::DIMS,
            knots: vec![(range[0], values[0]), (range[1], values[1])],
            dims: Pressure::DIMS,
        },
        validity: ValidityDomain::unconstrained().with(axis, range[0], range[1]),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::LinearInside,
        observations: vec![],
        provenance: Provenance {
            source: format!("manufactured thermal junction: {name}"),
            license: "CC0-1.0".into(),
            artifact: None,
        },
    }
}

fn strengths(axis: &str) -> Vec<PropertyClaim> {
    vec![
        curve(
            STATIC_JUNCTION_SHEAR_PROPERTY,
            axis,
            [280.0, 320.0],
            [3e6, 1.5e6],
        ),
        curve(
            KINETIC_JUNCTION_SHEAR_PROPERTY,
            axis,
            [280.0, 320.0],
            [2e6, 1e6],
        ),
    ]
}

fn point() -> QueryPoint {
    QueryPoint::new()
        .with("T", 300.0)
        .unwrap()
        .with_quantity(
            "p",
            QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static)),
            2e5,
        )
        .unwrap()
}

fn normal_with(claims: Vec<PropertyClaim>, point: QueryPoint) -> BoundNormalInterfaceState {
    let a = common::solid("thermal-junction-A", 1000.0, 1e7);
    let b = common::solid("thermal-junction-B", 2000.0, 2e7);
    let mut all = ClaimSet::new();
    // The normal source itself has thermal dependence: this split must retain
    // its original W and area rather than silently refreshing that other owner.
    let mut adhesion = claim(ADHESION_ENERGY_PROPERTY, 0.05, ADHESION_ENERGY_DIMS, false);
    adhesion.value = PropertyValue::Curve {
        abscissa: "T".into(),
        abscissa_dims: Temperature::DIMS,
        knots: vec![(280.0, 0.04), (320.0, 0.06)],
        dims: ADHESION_ENERGY_DIMS,
    };
    adhesion.interpolation = InterpolationPolicy::LinearInside;
    all.insert_claim(adhesion).unwrap();
    for claim in claims {
        all.insert_claim(claim).unwrap();
    }
    let card = InterfaceSystemCard::assemble(
        SurfaceSpec {
            material: a.id().clone(),
            texture_frame: "A-frame".into(),
        },
        SurfaceSpec {
            material: b.id().clone(),
            texture_frame: "B-frame".into(),
        },
        SystemContext {
            medium: "dry".into(),
            third_body: None,
            environment: "manufactured thermal test".into(),
            history: "initial".into(),
        },
        all,
        vec![],
    )
    .unwrap();
    let selection = MaterialPropertySelection::SingleClaimOnly;
    let sa = resolve_isotropic_solid_state_point(&a, &point, selection.clone()).unwrap();
    let sb = resolve_isotropic_solid_state_point(&b, &point, selection.clone()).unwrap();
    let dry = bind_dry_interface_system_card(&card, InputAuthority::SyntheticFixture).unwrap();
    let elastic = bind_isotropic_elastic_interface(&dry, &sa, &sb).unwrap();
    let normal = resolve_interface_state_point(
        &card,
        &point,
        &[ScalarPropertyRequirement::try_new(
            ADHESION_ENERGY_PROPERTY,
            ADHESION_ENERGY_DIMS,
            ScalarAdmissibility::NonNegative,
        )
        .unwrap()],
        selection,
    )
    .unwrap();
    bind_normal_interface_state(&elastic, &normal).unwrap()
}

fn jkr(normal: &BoundNormalInterfaceState) -> SourceJkrSphereContact {
    SourceJkrSphereContact::new(
        normal.clone(),
        0.01,
        1e-9,
        JkrEnvelope {
            maximum_indentation_m: 1e-5,
            half_space_depth_m: 0.1,
            layer_thickness_m: 0.1,
            maximum_patch_to_radius: 0.1,
            maximum_patch_to_depth: 0.1,
            maximum_patch_to_layer: 0.1,
            maximum_nominal_strain: 0.1,
            maximum_repulsive_pressure_to_yield: 0.1,
            maximum_process_zone_to_patch: 0.01,
            minimum_tabor_parameter: 5.0,
        },
    )
    .unwrap()
}

#[test]
fn thermal_junction_torque_and_heat_follow_selected_strength_at_zero_and_tensile_load() {
    let normal = normal_with(strengths("T"), point());
    let shear =
        SourceAdhesiveShear::resolve(&normal, MaterialPropertySelection::SingleClaimOnly).unwrap();
    let owner = jkr(&normal);
    let radius = 0.01;
    let a0 = (4.5 * PI * normal.adhesion_energy_j_per_m2() * radius * radius
        / normal.elastic().reduced_modulus_pa())
    .cbrt();
    let zero_indentation = a0 * a0 / radius
        - (2.0 * PI * normal.adhesion_energy_j_per_m2() * a0
            / normal.elastic().reduced_modulus_pa())
        .sqrt();
    let zero = owner.evaluate_attached(zero_indentation).unwrap();
    let tensile = owner.evaluate_attached(0.0).unwrap();
    assert!(zero.normal_force_n.abs() < 1e-14);
    assert!(tensile.normal_force_n < 0.0);
    let before_normal = normal.clone();
    let before_shear = shear.clone();
    for response in [zero, tensile] {
        let original = shear.sliding_patch_jkr(&normal, &response).unwrap();
        let trial = original.at_temperature("T", 320.0).unwrap();
        assert_eq!(trial.original(), original);
        assert_eq!(trial.friction_source().state_coordinate("T"), Some(320.0));
        assert_eq!(
            trial
                .original()
                .normal_state()
                .interface_state()
                .state_coordinate("T"),
            Some(300.0)
        );
        assert_eq!(
            trial.law().parameters().contact_radius_m,
            response.contact_radius_m
        );
        common::close(trial.law().parameters().static_shear_strength_pa, 1.5e6);
        common::close(trial.law().parameters().kinetic_shear_strength_pa, 1e6);
        let spin = trial.law().kinetic_patch().evaluate([0.0; 2], 5.0).unwrap();
        let moment = (2.0 / 3.0) * 1e6 * response.contact_area_m2 * response.contact_radius_m;
        assert_eq!(spin.force_n, [0.0; 2]);
        common::close(spin.spin_moment_nm, -moment);
        common::close(spin.dissipated_power_w, 5.0 * moment);
        common::close(
            spin.dissipated_power_w / original.evaluate([0.0; 2], 5.0).unwrap().dissipated_power_w,
            2.0 / 3.0,
        );
    }
    assert_eq!(normal, before_normal);
    assert_eq!(shear, before_shear);
}

#[test]
fn refresh_pins_claims_schema_and_nonthermal_coordinates_when_other_claims_become_eligible() {
    let mut claims = strengths("T");
    claims.push(curve(
        KINETIC_JUNCTION_SHEAR_PROPERTY,
        "T",
        [310.0, 320.0],
        [9e6, 9e6],
    ));
    let normal = normal_with(claims, point());
    let shear =
        SourceAdhesiveShear::resolve(&normal, MaterialPropertySelection::SingleClaimOnly).unwrap();
    let response = jkr(&normal).evaluate_attached(0.0).unwrap();
    let original = shear.sliding_patch_jkr(&normal, &response).unwrap();
    let trial = original.at_temperature("T", 315.0).unwrap();
    common::close(trial.law().parameters().kinetic_shear_strength_pa, 1.125e6);
    assert_eq!(
        trial.friction_source().card_identity(),
        shear.resolved().card_identity()
    );
    assert_eq!(
        trial.friction_source().axis_quantities(),
        shear.resolved().axis_quantities()
    );
    assert_eq!(trial.friction_source().state_coordinate("p"), Some(2e5));
    for (old, new) in shear
        .resolved()
        .properties()
        .iter()
        .zip(trial.friction_source().properties())
    {
        assert_eq!(old.requirement(), new.requirement());
        assert_eq!(old.answer().receipt.selected, new.answer().receipt.selected);
    }
}

#[test]
fn expired_pinned_claim_invalid_temperature_and_detachment_refuse_without_mutating_sources() {
    let mut claims = strengths("T");
    claims[1] = curve(
        KINETIC_JUNCTION_SHEAR_PROPERTY,
        "T",
        [280.0, 305.0],
        [2e6, 1e6],
    );
    claims.push(curve(
        KINETIC_JUNCTION_SHEAR_PROPERTY,
        "T",
        [310.0, 320.0],
        [1e6, 1e6],
    ));
    let normal = normal_with(claims, point());
    let shear =
        SourceAdhesiveShear::resolve(&normal, MaterialPropertySelection::SingleClaimOnly).unwrap();
    let owner = jkr(&normal);
    let response = owner.evaluate_attached(0.0).unwrap();
    let original = shear.sliding_patch_jkr(&normal, &response).unwrap();
    let before = (normal.clone(), shear.clone());
    assert!(matches!(
        original.at_temperature("T", 315.0),
        Err(SourceThermalAdhesivePatchError::Material(_))
    ));
    for value in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(matches!(
            original.at_temperature("T", value),
            Err(SourceThermalAdhesivePatchError::Invalid(_))
        ));
    }
    for axis in ["p", "missing", "temperature"] {
        assert!(matches!(
            original.at_temperature(axis, 300.0),
            Err(SourceThermalAdhesivePatchError::Invalid(_))
        ));
    }
    let detached = owner
        .advance_monotone(owner.detached_state(-1e-6).unwrap(), -1e-6)
        .unwrap()
        .response;
    let empty = shear.sliding_patch_jkr(&normal, &detached).unwrap();
    assert!(matches!(
        empty.at_temperature("T", 300.0),
        Err(SourceThermalAdhesivePatchError::Invalid(_))
    ));
    assert_eq!(normal, before.0);
    assert_eq!(shear, before.1);
    assert_eq!(
        original
            .at_temperature("T", 300.0)
            .unwrap()
            .friction_source(),
        shear.resolved()
    );
}

#[test]
fn semantic_face_temperature_is_owned_but_a_second_temperature_bound_is_not_frozen_silently() {
    let face_temperature = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ));
    let query = point()
        .with_quantity("face_T", face_temperature, 300.0)
        .unwrap();
    for second_temperature_bound in [false, true] {
        let mut claims = strengths("face_T");
        // The claims' face_T axis carries the quantity the query gives face_T:
        // a query axis typed differently from its claim's axis is outside the
        // claim's domain (AxisQuantityMismatch).
        for claim in &mut claims {
            claim.validity = ValidityDomain::unconstrained().with_quantity(
                "face_T",
                face_temperature,
                280.0,
                320.0,
            );
        }
        if second_temperature_bound {
            claims[0].validity = claims[0].validity.clone().with("T", 280.0, 320.0);
        }
        let normal = normal_with(claims, query.clone());
        let shear =
            SourceAdhesiveShear::resolve(&normal, MaterialPropertySelection::SingleClaimOnly)
                .unwrap();
        let response = jkr(&normal).evaluate_attached(0.0).unwrap();
        let original = shear.sliding_patch_jkr(&normal, &response).unwrap();
        if second_temperature_bound {
            assert!(matches!(
                original.at_temperature("face_T", 310.0),
                Err(SourceThermalAdhesivePatchError::Invalid(_))
            ));
        } else {
            let trial = original.at_temperature("face_T", 310.0).unwrap();
            common::close(trial.law().parameters().kinetic_shear_strength_pa, 1.25e6);
            assert_eq!(
                trial.friction_source().state_coordinate("face_T"),
                Some(310.0)
            );
            assert_eq!(trial.friction_source().state_coordinate("T"), Some(300.0));
        }
    }
}

#[test]
fn temperature_strength_crossing_refuses_the_original_uniform_junction_law() {
    let mut claims = strengths("T");
    claims[1] = curve(
        KINETIC_JUNCTION_SHEAR_PROPERTY,
        "T",
        [280.0, 320.0],
        [2e6, 2e6],
    );
    let normal = normal_with(claims, point());
    let shear =
        SourceAdhesiveShear::resolve(&normal, MaterialPropertySelection::SingleClaimOnly).unwrap();
    let response = jkr(&normal).evaluate_attached(0.0).unwrap();
    let original = shear.sliding_patch_jkr(&normal, &response).unwrap();
    assert!(matches!(
        original.at_temperature("T", 320.0),
        Err(SourceThermalAdhesivePatchError::Wrench(_))
    ));
    assert_eq!(
        original
            .at_temperature("T", 300.0)
            .unwrap()
            .friction_source(),
        shear.resolved()
    );
}
