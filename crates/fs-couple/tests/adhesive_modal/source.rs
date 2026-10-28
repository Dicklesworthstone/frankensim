//! Actual source-backed DMT area/load changes in the native mechanical owner.
use super::*;
#[path = "../../../fs-contact/tests/source_point/common.rs"]
mod common;
use fs_contact::interface_binding::adhesive_shear::{
    KINETIC_JUNCTION_SHEAR_PROPERTY, STATIC_JUNCTION_SHEAR_PROPERTY, SourceAdhesiveContactUpdate,
    SourceAdhesiveShear,
};
use fs_contact::interface_binding::*;
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
use fs_qty::Pressure;
use fs_tribo::InputAuthority;

fn normal(temperature: f64) -> SourceDmtSphereContact {
    let a = common::solid("moving-junction-A", 1000.0, 1e11);
    let b = common::solid("moving-junction-B", 2000.0, 2e11);
    let mut claims = ClaimSet::new();
    for (name, value, dims) in [
        (ADHESION_ENERGY_PROPERTY, 0.05, ADHESION_ENERGY_DIMS),
        (STATIC_JUNCTION_SHEAR_PROPERTY, 3e6, Pressure::DIMS),
        (KINETIC_JUNCTION_SHEAR_PROPERTY, 2e6, Pressure::DIMS),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained().with("T", 280.0, 320.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source: format!("synthetic native adhesive geometry test: {name}"),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = InterfaceSystemCard::assemble(
        SurfaceSpec {
            material: a.id().clone(),
            texture_frame: "A".into(),
        },
        SurfaceSpec {
            material: b.id().clone(),
            texture_frame: "B".into(),
        },
        SystemContext {
            medium: "dry".into(),
            third_body: None,
            environment: "synthetic test".into(),
            history: "initial".into(),
        },
        claims,
        vec![],
    )
    .unwrap();
    let point = QueryPoint::new().with("T", temperature).unwrap();
    let selection = MaterialPropertySelection::SingleClaimOnly;
    let a = resolve_isotropic_solid_state_point(&a, &point, selection.clone()).unwrap();
    let b = resolve_isotropic_solid_state_point(&b, &point, selection.clone()).unwrap();
    let dry = bind_dry_interface_system_card(&card, InputAuthority::SyntheticFixture).unwrap();
    let elastic = bind_isotropic_elastic_interface(&dry, &a, &b).unwrap();
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
    SourceDmtSphereContact::new(
        bind_normal_interface_state(&elastic, &interface).unwrap(),
        0.01,
        1e-6,
        DmtEnvelope {
            maximum_indentation_m: 1e-6,
            half_space_depth_m: 0.1,
            layer_thickness_m: 0.1,
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
    .unwrap()
}

#[test]
fn real_normal_responses_update_scalar_shear_without_replacing_original_sources() {
    let normal = normal(293.15);
    let shear = SourceAdhesiveShear::resolve(
        normal.normal_state(),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let first = normal.evaluate(0.2e-6).unwrap();
    let binding =
        AdhesiveModalBinding::from_source(shear.bind_dmt(normal.normal_state(), &first).unwrap())
            .unwrap();
    let mut s = mechanics_with_binding(1.0, binding);
    let receipts = s
        .adhesive_binding()
        .unwrap()
        .source()
        .unwrap()
        .shear_source()
        .resolved()
        .properties()
        .as_ptr();
    let initial_energy = s.total_energy_j().unwrap();
    let initial_momentum = momentum(&s);
    let mut heat = 0.0;
    for indentation in [0.8e-6, -0.5e-6, 0.4e-6] {
        let response = normal.evaluate(indentation).unwrap();
        let update = SourceAdhesiveContactUpdate::from_dmt(&response).unwrap();
        let frame = *s
            .step_with_adhesive_contact(
                &[0.0, 0.0],
                FrictionDrive {
                    normal_force_n: response.repulsive_normal_force_n,
                    ..drive()
                },
                AdhesiveModalContactUpdate::Source(update),
            )
            .unwrap();
        close(frame.traction_n, -2e6 * response.contact_area_m2, 1e-12);
        assert_eq!(
            frame.regime,
            if indentation > 0.0 {
                ModalFrictionRegime::Sliding
            } else {
                ModalFrictionRegime::Released
            }
        );
        let binding = s.adhesive_binding().unwrap();
        let source = binding.source().unwrap();
        assert_eq!(binding.tangential(), source.tangential());
        assert_eq!(
            source.tangential().contact(),
            response.friction_contact_state().unwrap()
        );
        assert_eq!(
            source.normal_state().identity(),
            normal.normal_state().identity()
        );
        assert_eq!(source.temperature_kelvin(), Some(293.15));
        assert_eq!(source.shear_source().resolved(), shear.resolved());
        assert_eq!(
            source.shear_source().resolved().properties().as_ptr(),
            receipts
        );
        assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
        heat += frame.friction_dissipation_j;
        close(momentum(&s), initial_momentum, 1e-14);
    }
    close(initial_energy - s.total_energy_j().unwrap(), heat, 1e-13);
    let states = states(&s);
    let last = s.last_frame().copied();
    let binding = s.adhesive_binding().unwrap().clone();
    let rejected = normal.evaluate(0.9e-6).unwrap();
    assert!(
        s.step_with_adhesive_contact(
            &[0.0, 0.0],
            FrictionDrive {
                speed_m_s: 20.0,
                normal_force_n: rejected.repulsive_normal_force_n
            },
            AdhesiveModalContactUpdate::Source(
                SourceAdhesiveContactUpdate::from_dmt(&rejected).unwrap()
            )
        )
        .is_err()
    );
    assert_eq!(s.adhesive_binding().unwrap(), &binding);
    // Even identical numeric coefficients at another source temperature are
    // a different state, not permission to overwrite this original query.
    let foreign = self::normal(300.0).evaluate(0.4e-6).unwrap();
    let update = SourceAdhesiveContactUpdate::from_dmt(&foreign).unwrap();
    assert!(
        s.step_with_adhesive_contact(
            &[0.0, 0.0],
            FrictionDrive {
                normal_force_n: foreign.repulsive_normal_force_n,
                ..drive()
            },
            AdhesiveModalContactUpdate::Source(update)
        )
        .is_err()
    );
    assert!(
        s.step_with_adhesive_contact(&[0.0, 0.0], drive(), geometry(0.0, 0.0))
            .is_err()
    );
    assert_eq!(super::states(&s), states);
    assert_eq!(s.last_frame().copied(), last);
    assert_eq!(s.adhesive_binding().unwrap(), &binding);
}
