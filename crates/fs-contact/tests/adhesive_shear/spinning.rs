//! The original normal patch and shear claims determine torsional resistance.
use super::*;
use core::f64::consts::PI;
use fs_contact::normal_patch::maugis::{MaugisEnvelope, SourceMaugisSphereContact};

fn jkr(normal: BoundNormalInterfaceState) -> SourceJkrSphereContact {
    SourceJkrSphereContact::new(
        normal,
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

fn stiff(normal: &BoundNormalInterfaceState) -> BoundNormalInterfaceState {
    let point = normal.interface_state().query();
    let a = common::solid("junction-A", 1000.0, 1e11);
    let b = common::solid("junction-B", 2000.0, 2e11);
    let selection = MaterialPropertySelection::SingleClaimOnly;
    let sa = resolve_isotropic_solid_state_point(&a, &point, selection.clone()).unwrap();
    let sb = resolve_isotropic_solid_state_point(&b, &point, selection).unwrap();
    let elastic = bind_isotropic_elastic_interface(normal.elastic().interface(), &sa, &sb).unwrap();
    bind_normal_interface_state(&elastic, normal.interface_state()).unwrap()
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-13 * actual.abs().max(expected.abs()).max(1e-30),
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn tensile_jkr_patch_transmits_torsion_using_original_strengths_and_receipts() {
    let normal = normal(Some((3e6, 2e6)), Pressure::DIMS, false, false);
    let source =
        SourceAdhesiveShear::resolve(&normal, MaterialPropertySelection::SingleClaimOnly).unwrap();
    let owner = jkr(normal.clone());
    let response = owner.evaluate_attached(0.0).unwrap();
    assert!(response.normal_force_n < 0.0);
    let binding = source.sliding_patch_jkr(&normal, &response).unwrap();
    assert!(std::ptr::eq(binding.normal_state(), &normal));
    assert!(std::ptr::eq(binding.shear_source(), &source));
    assert_eq!(binding.contact_area_m2(), response.contact_area_m2);
    close(
        binding.patch().sliding_capacity_n(),
        2e6 * response.contact_area_m2,
    );
    close(
        binding.static_patch().sliding_capacity_n(),
        3e6 * response.contact_area_m2,
    );
    let spin = binding.evaluate([0.0; 2], 4.0).unwrap();
    let torque = (2.0 / 3.0) * 2e6 * response.contact_area_m2 * response.contact_radius_m;
    close(spin.spin_moment_nm, -torque);
    close(spin.dissipated_power_w, 4.0 * torque);
    assert_eq!(spin.force_n, [0.0; 2]);
    // Detachment removes the solid-area torque; snap energy is not friction heat.
    let detached = owner
        .advance_monotone(
            response.state(),
            1.1 * owner.minimum_attached_indentation_m(),
        )
        .unwrap();
    assert!(detached.irreversible_work_j > 0.0);
    let empty = source
        .sliding_patch_jkr(&normal, &detached.response)
        .unwrap();
    assert_eq!(
        empty
            .evaluate([0.03, 0.04], 4.0)
            .unwrap()
            .dissipated_power_w,
        0.0
    );
    assert_eq!(normal, owner.normal_state().clone());
}

#[test]
fn dmt_gap_attraction_has_no_torsional_junction_but_zero_load_contact_does() {
    let soft = normal(Some((3e6, 2e6)), Pressure::DIMS, false, false);
    let normal = stiff(&soft);
    let source =
        SourceAdhesiveShear::resolve(&normal, MaterialPropertySelection::SingleClaimOnly).unwrap();
    let owner = SourceDmtSphereContact::new(
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
    let contact = owner.evaluate_contact_load(0.0).unwrap();
    let held = source.sliding_patch_dmt(&normal, &contact).unwrap();
    assert!(held.evaluate([0.0; 2], 10.0).unwrap().dissipated_power_w > 0.0);
    for indentation in [0.0, -0.5e-9, -2e-9] {
        let gap = owner.evaluate(indentation).unwrap();
        let empty = source.sliding_patch_dmt(&normal, &gap).unwrap();
        assert_eq!(empty.contact_area_m2(), 0.0);
        let spin = empty.evaluate([0.03, 0.04], 10.0).unwrap();
        assert_eq!(spin.force_n, [0.0; 2]);
        assert_eq!(spin.spin_moment_nm, 0.0);
        assert_eq!(spin.dissipated_power_w, 0.0);
    }
    assert!(source.sliding_patch_dmt(&soft, &contact).is_err());
}

#[test]
fn maugis_solid_radius_excludes_cohesive_annulus_from_spin_friction() {
    let normal = normal(Some((3e6, 2e6)), Pressure::DIMS, false, false);
    let source =
        SourceAdhesiveShear::resolve(&normal, MaterialPropertySelection::SingleClaimOnly).unwrap();
    let radius = 0.01;
    let e = normal.elastic().reduced_modulus_pa();
    let w = normal.adhesion_energy_j_per_m2();
    let a0 = (PI * w * radius * radius / (4.0 * e / 3.0)).cbrt();
    let d0 = a0 * a0 / radius;
    let owner = SourceMaugisSphereContact::new(
        normal.clone(),
        radius,
        2.0 * d0 / PI,
        MaugisEnvelope {
            maximum_indentation_m: 20.0 * d0,
            half_space_depth_m: 0.1,
            layer_thickness_m: 0.1,
            maximum_interaction_radius_to_radius: 0.1,
            maximum_interaction_radius_to_depth: 0.1,
            maximum_interaction_radius_to_layer: 0.1,
            maximum_nominal_strain: 0.1,
            maximum_repulsive_pressure_to_yield: 0.5,
            maximum_adhesive_stress_to_yield: 0.9,
            maximum_interaction_range_to_radius: 0.01,
        },
    )
    .unwrap();
    let response = owner.evaluate_contact_load(0.0).unwrap();
    assert!(response.cohesive_radius_m > response.contact_radius_m);
    let binding = source.sliding_patch_maugis(&normal, &response).unwrap();
    assert_eq!(
        binding.patch().parameters().contact_radius_m,
        response.contact_radius_m
    );
    close(
        binding.patch().spinning_capacity_nm(),
        (2.0 / 3.0) * 2e6 * response.contact_area_m2 * response.contact_radius_m,
    );
    assert!(
        binding.patch().spinning_capacity_nm()
            < (2.0 / 3.0) * 2e6 * PI * response.cohesive_radius_m.powi(3)
    );
}

#[test]
fn incompatible_normal_sources_or_noncircular_diagnostics_cannot_relabel_a_patch() {
    let original = normal(Some((3e6, 2e6)), Pressure::DIMS, false, false);
    let source =
        SourceAdhesiveShear::resolve(&original, MaterialPropertySelection::SingleClaimOnly)
            .unwrap();
    let response = jkr(original.clone()).evaluate_attached(0.0).unwrap();
    for foreign in [
        normal(Some((3e6, 2e6)), Pressure::DIMS, false, true),
        normal(Some((4e6, 2e6)), Pressure::DIMS, false, false),
    ] {
        let wrong = jkr(foreign.clone()).evaluate_attached(0.0).unwrap();
        assert!(source.sliding_patch_jkr(&foreign, &response).is_err());
        assert!(source.sliding_patch_jkr(&original, &wrong).is_err());
    }
    let mut wrong_area = response;
    wrong_area.contact_area_m2 *= 2.0;
    assert!(source.sliding_patch_jkr(&original, &wrong_area).is_err());
    assert_eq!(source.normal_state_identity(), original.identity());
}
