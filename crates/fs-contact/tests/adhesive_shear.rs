//! Source identity, physical data and zero-load junction shear through real resolvers.
#[path = "source_point/common.rs"]
mod common;

use fs_contact::interface_binding::adhesive_shear::*;
use fs_contact::interface_binding::*;
use fs_contact::normal_patch::adhesive::{JkrEnvelope, SourceJkrSphereContact};
use fs_contact::normal_patch::dmt::{DmtEnvelope, SourceDmtSphereContact};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterfaceSystemCard, InterpolationPolicy, PropertyClaim, PropertyKey, PropertyValue,
    Provenance, QueryPoint, SurfaceSpec, SystemContext, UncertaintyModel,
};
use fs_material::state_point::{
    MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
    resolve_interface_state_point, resolve_isotropic_solid_state_point,
};
use fs_qty::{Dims, Pressure};
use fs_tribo::adhesive::AdhesiveContactState;
use fs_tribo::{
    ContactFrame, FrictionRegime, InputAuthority, TangentialLaw, TangentialSlip,
    TangentialSolveConfig, TangentialStepInput,
};

fn claim(name: &str, value: f64, dims: Dims, rate_axis: bool) -> PropertyClaim {
    let mut validity = ValidityDomain::unconstrained().with("T", 280.0, 320.0);
    if rate_axis {
        validity = validity.with("slip_speed", 0.0, 1.0);
    }
    PropertyClaim {
        key: PropertyKey::new(name, dims),
        value: PropertyValue::Scalar { value, dims },
        validity,
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        observations: vec![],
        provenance: Provenance {
            source: format!("synthetic junction test: {name}"),
            license: "CC0-1.0".into(),
            artifact: None,
        },
    }
}

fn normal(
    strengths: Option<(f64, f64)>,
    shear_dims: Dims,
    rate_axis: bool,
    reversed: bool,
) -> BoundNormalInterfaceState {
    let a = common::solid("junction-A", 1000.0, 1e7);
    let b = common::solid("junction-B", 2000.0, 2e7);
    let (a, b) = if reversed { (&b, &a) } else { (&a, &b) };
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(claim(
            ADHESION_ENERGY_PROPERTY,
            0.05,
            ADHESION_ENERGY_DIMS,
            false,
        ))
        .unwrap();
    if let Some((static_strength, kinetic_strength)) = strengths {
        for (name, value) in [
            (STATIC_JUNCTION_SHEAR_PROPERTY, static_strength),
            (KINETIC_JUNCTION_SHEAR_PROPERTY, kinetic_strength),
        ] {
            claims
                .insert_claim(claim(name, value, shear_dims, rate_axis))
                .unwrap();
        }
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
            environment: "synthetic-test".into(),
            history: "initial".into(),
        },
        claims,
        vec![],
    )
    .unwrap();
    let mut point = QueryPoint::new().with("T", 293.15).unwrap();
    if rate_axis {
        point = point.with("slip_speed", 0.0).unwrap();
    }
    let selection = MaterialPropertySelection::SingleClaimOnly;
    let sa = resolve_isotropic_solid_state_point(a, &point, selection.clone()).unwrap();
    let sb = resolve_isotropic_solid_state_point(b, &point, selection.clone()).unwrap();
    let dry = bind_dry_interface_system_card(&card, InputAuthority::SyntheticFixture).unwrap();
    let elastic = bind_isotropic_elastic_interface(&dry, &sa, &sb).unwrap();
    let interface = resolve_interface_state_point(
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
    bind_normal_interface_state(&elastic, &interface).unwrap()
}

fn slip() -> TangentialSlip {
    TangentialSlip::new(
        &ContactFrame::new([0.0, 0.0, 1.0]).unwrap(),
        [0.03, 0.04, 0.0],
    )
    .unwrap()
}

#[test]
fn original_interface_strengths_produce_area_shear_and_keep_receipts() {
    let normal = normal(Some((3e6, 2e6)), Pressure::DIMS, false, false);
    let source =
        SourceAdhesiveShear::resolve(&normal, MaterialPropertySelection::SingleClaimOnly).unwrap();
    assert_eq!(
        source.resolved().card_identity(),
        normal.interface_state().card_identity()
    );
    assert_eq!(
        source.resolved().query_point(),
        normal.interface_state().query_point()
    );
    assert_eq!(source.resolved().properties().len(), 2);
    let contact = AdhesiveContactState::new(2e-6, 0.0).unwrap();
    let response = source.evaluate(&normal, contact, slip()).unwrap();
    common::close(response.static_capacity_n, 6.0);
    common::close(response.sliding_capacity_n.unwrap(), 4.0);
    common::close(response.traction_n[0], -2.4);
    common::close(response.traction_n[1], -3.2);
    common::close(response.dissipated_power_w, 0.2);
    assert!(std::ptr::eq(
        response.provenance,
        normal.elastic().interface().interface().provenance()
    ));
    // This explicitly selected model has no extra Coulomb contribution.
    let loaded = source
        .evaluate(
            &normal,
            AdhesiveContactState::new(2e-6, 1e6).unwrap(),
            slip(),
        )
        .unwrap();
    assert_eq!(loaded.traction_n, response.traction_n);
}

#[test]
fn missing_strengths_and_wrong_units_do_not_fall_back_to_bulk_or_adhesion() {
    for normal in [
        normal(None, Pressure::DIMS, false, false),
        normal(Some((3e6, 2e6)), Dims::NONE, false, false),
    ] {
        assert!(
            SourceAdhesiveShear::resolve(&normal, MaterialPropertySelection::SingleClaimOnly,)
                .is_err()
        );
    }
}

#[test]
fn rate_restrictions_and_inverted_strengths_refuse() {
    for normal in [
        normal(Some((3e6, 2e6)), Pressure::DIMS, true, false),
        normal(Some((2e6, 3e6)), Pressure::DIMS, false, false),
        normal(Some((-1.0, 0.0)), Pressure::DIMS, false, false),
    ] {
        assert!(
            SourceAdhesiveShear::resolve(&normal, MaterialPropertySelection::SingleClaimOnly,)
                .is_err()
        );
    }
}

#[test]
fn reversed_surfaces_and_changed_strengths_cannot_reuse_source_binding() {
    let original = normal(Some((3e6, 2e6)), Pressure::DIMS, false, false);
    let source =
        SourceAdhesiveShear::resolve(&original, MaterialPropertySelection::SingleClaimOnly)
            .unwrap();
    let contact = AdhesiveContactState::new(2e-6, 0.0).unwrap();
    for other in [
        normal(Some((3e6, 2e6)), Pressure::DIMS, false, true),
        normal(Some((4e6, 2e6)), Pressure::DIMS, false, false),
    ] {
        assert!(matches!(
            source.evaluate(&other, contact, slip()),
            Err(AdhesiveShearError::Invalid { .. })
        ));
    }
}

#[test]
fn detached_geometry_produces_no_shear_or_frictional_heat() {
    let normal = normal(Some((3e6, 2e6)), Pressure::DIMS, false, false);
    let source =
        SourceAdhesiveShear::resolve(&normal, MaterialPropertySelection::SingleClaimOnly).unwrap();
    let response = source
        .evaluate(
            &normal,
            AdhesiveContactState::new(0.0, 0.0).unwrap(),
            slip(),
        )
        .unwrap();
    assert_eq!(response.static_capacity_n, 0.0);
    assert_eq!(response.traction_n, [0.0; 3]);
    assert_eq!(response.dissipated_power_w, 0.0);
}

#[test]
fn sourced_jkr_area_drives_shear_at_zero_and_tensile_net_force_then_detaches() {
    let normal = normal(Some((3e6, 2e6)), Pressure::DIMS, false, false);
    let shear =
        SourceAdhesiveShear::resolve(&normal, MaterialPropertySelection::SingleClaimOnly).unwrap();
    let radius = 0.01;
    let jkr = SourceJkrSphereContact::new(
        normal.clone(),
        radius,
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
    .unwrap();
    let w = normal.adhesion_energy_j_per_m2();
    let e = normal.elastic().reduced_modulus_pa();
    // Independent zero-load JKR radius: a0^3 = 9 pi W R^2 / (2 E*).
    let a0 = (4.5 * core::f64::consts::PI * w * radius * radius / e).cbrt();
    let delta0 = a0 * a0 / radius - (2.0 * core::f64::consts::PI * w * a0 / e).sqrt();
    let zero_load = jkr.evaluate_attached(delta0).unwrap();
    assert!(zero_load.normal_force_n.abs() < 1e-13 * w * radius);
    common::close(zero_load.contact_radius_m, a0);
    let tensile = jkr.evaluate_attached(0.0).unwrap();
    assert!(tensile.normal_force_n < 0.0);
    for contact in [zero_load, tensile] {
        let response = shear.evaluate_jkr(&normal, &contact, slip()).unwrap();
        common::close(response.static_capacity_n, 3e6 * contact.contact_area_m2);
        common::close(
            response.sliding_capacity_n.unwrap(),
            2e6 * contact.contact_area_m2,
        );
        common::close(
            response.dissipated_power_w,
            2e6 * contact.contact_area_m2 * 0.05,
        );
        let held = shear.bind_jkr(&normal, &contact).unwrap();
        assert_eq!(held.normal_state().identity(), normal.identity());
        assert_eq!(held.shear_source().resolved(), shear.resolved());
        let x = TangentialStepInput {
            normal_n: contact.repulsive_normal_force_n,
            slip_increment_m: [0.5 * response.static_capacity_n, 0.0],
            dt_s: 1.0,
            contact_temperature_k: Some(293.15),
            nominal_pressure_pa: None,
            compliance_m_per_n: 1.0,
            solve: TangentialSolveConfig {
                max_iterations: 80,
                maximum_force_n: 10.0,
                force_absolute_tolerance_n: 1e-14,
                force_relative_tolerance: 1e-12,
            },
        };
        let initial = held.initial_state();
        let stuck = held.advance(&initial, &x).unwrap();
        assert_eq!(stuck.regime, FrictionRegime::Sticking);
        common::close(stuck.traction_n[0], 0.5 * response.static_capacity_n);
        assert_eq!(stuck.dissipated_j, 0.0);
        let slide = held
            .advance(
                &initial,
                &TangentialStepInput {
                    slip_increment_m: [2.0 * response.static_capacity_n, 0.0],
                    ..x
                },
            )
            .unwrap();
        assert_eq!(slide.regime, FrictionRegime::Sliding);
        common::close(slide.traction_n[0], response.sliding_capacity_n.unwrap());
        assert!(slide.dissipated_j > 0.0);
        assert!(
            held.advance(
                &initial,
                &TangentialStepInput {
                    contact_temperature_k: Some(300.0),
                    ..x
                }
            )
            .is_err()
        );
        assert!(
            held.advance(
                &initial,
                &TangentialStepInput {
                    normal_n: -1.0,
                    ..x
                }
            )
            .is_err()
        );
        assert_eq!(initial, held.initial_state());
        assert_eq!(contact.authority, InputAuthority::SyntheticFixture);
        assert_eq!(
            response.provenance.authority(),
            InputAuthority::SyntheticFixture
        );
    }
    let detached = jkr
        .advance_monotone(tensile.state(), 1.1 * jkr.minimum_attached_indentation_m())
        .unwrap();
    let response = shear
        .evaluate_jkr(&normal, &detached.response, slip())
        .unwrap();
    assert_eq!(response.static_capacity_n, 0.0);
    assert_eq!(response.dissipated_power_w, 0.0);
    assert!(detached.irreversible_work_j > 0.0);
}

#[test]
fn sourced_dmt_distinguishes_contact_friction_from_attraction_across_a_gap() {
    let soft = normal(Some((3e6, 2e6)), Pressure::DIMS, false, false);
    let point = soft.interface_state().query();
    let a = common::solid("junction-A", 1000.0, 1e11);
    let b = common::solid("junction-B", 2000.0, 2e11);
    let sa =
        resolve_isotropic_solid_state_point(&a, &point, MaterialPropertySelection::SingleClaimOnly)
            .unwrap();
    let sb =
        resolve_isotropic_solid_state_point(&b, &point, MaterialPropertySelection::SingleClaimOnly)
            .unwrap();
    let elastic = bind_isotropic_elastic_interface(soft.elastic().interface(), &sa, &sb).unwrap();
    let normal = bind_normal_interface_state(&elastic, soft.interface_state()).unwrap();
    let shear =
        SourceAdhesiveShear::resolve(&normal, MaterialPropertySelection::SingleClaimOnly).unwrap();
    let dmt = SourceDmtSphereContact::new(
        normal.clone(),
        1e-6,
        1e-9,
        DmtEnvelope {
            maximum_indentation_m: 1e-9,
            half_space_depth_m: 1e-3,
            layer_thickness_m: 1e-3,
            maximum_interaction_radius_to_radius: 0.1,
            maximum_interaction_radius_to_depth: 0.1,
            maximum_interaction_radius_to_layer: 0.1,
            maximum_nominal_strain: 0.1,
            maximum_repulsive_pressure_to_yield: 0.3,
            maximum_adhesive_stress_to_yield: 0.1,
            maximum_interaction_range_to_radius: 0.01,
            maximum_tabor_parameter: 0.1,
        },
    )
    .unwrap();
    let contact = dmt.evaluate_contact_load(0.0).unwrap();
    let friction = shear.evaluate_dmt(&normal, &contact, slip()).unwrap();
    assert!(friction.dissipated_power_w > 0.0);
    common::close(
        friction.sliding_capacity_n.unwrap(),
        2e6 * contact.contact_area_m2,
    );
    let gap = dmt.evaluate(-0.5e-9).unwrap();
    assert!(gap.normal_force_n < 0.0);
    let friction = shear.evaluate_dmt(&normal, &gap, slip()).unwrap();
    assert_eq!(friction.static_capacity_n, 0.0);
    assert_eq!(friction.dissipated_power_w, 0.0);
    assert!(shear.evaluate_dmt(&soft, &contact, slip()).is_err());
}
