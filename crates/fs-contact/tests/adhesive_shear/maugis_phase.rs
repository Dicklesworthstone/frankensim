//! Complete finite-zone phases feed the original sourced tangential owner.
use super::*;
use core::f64::consts::PI;
use fs_contact::normal_patch::maugis::{
    MaugisEnvelope, SourceMaugisSphereContact, continuation::MaugisPhase,
};

fn model(normal: BoundNormalInterfaceState, radius: f64) -> SourceMaugisSphereContact {
    let a0 = (PI * normal.adhesion_energy_j_per_m2() * radius * radius
        / (4.0 * normal.elastic().reduced_modulus_pa() / 3.0))
        .cbrt();
    let d0 = a0 * a0 / radius;
    SourceMaugisSphereContact::new(
        normal,
        radius,
        2.0 * d0 / (PI * 0.2),
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
    .unwrap()
}

#[test]
fn complete_phases_preserve_solid_area_and_exclude_noncontact_attraction_from_shear() {
    let normal = normal(Some((3e6, 2e6)), Pressure::DIMS, false, false);
    let law = model(normal.clone(), 0.01);
    let shear =
        SourceAdhesiveShear::resolve(&normal, MaterialPropertySelection::SingleClaimOnly).unwrap();
    let contact = law.evaluate_contact_energy(0.0).unwrap();
    assert!(contact.normal_force_n() < 0.0);
    let limits = law.noncontact_limits().unwrap().unwrap();
    let gap = law
        .evaluate_noncontact(
            0.5 * limits.minimum_indentation_m + 0.5 * limits.first_contact_indentation_m,
        )
        .unwrap();
    let detached = law
        .evaluate_detached(-2.0 * law.interaction_range_m())
        .unwrap();
    assert_eq!(gap.phase(), MaugisPhase::Noncontact);
    assert!(gap.normal_force_n() < 0.0 && gap.cohesive_radius_m() > 0.0);
    let solid = contact.contact_response().unwrap();
    let mut held = shear.bind_maugis_equilibrium(&law, &contact).unwrap();
    assert_eq!(held, shear.bind_maugis(&normal, solid).unwrap());
    assert_eq!(
        SourceAdhesiveContactUpdate::from_maugis_equilibrium(&law, &contact).unwrap(),
        SourceAdhesiveContactUpdate::from_maugis(solid).unwrap()
    );
    let original_receipts = held.shear_source().resolved().properties().as_ptr();
    for equilibrium in [contact, gap, detached, contact] {
        let update =
            SourceAdhesiveContactUpdate::from_maugis_equilibrium(&law, &equilibrium).unwrap();
        let state = equilibrium.friction_contact_state().unwrap();
        assert_eq!(update.normal_state_identity(), normal.identity());
        assert_eq!(update.contact(), state);
        held.update_contact(update).unwrap();
        assert_eq!(held.tangential().contact(), state);
        assert_eq!(held.normal_state(), &normal);
        assert_eq!(
            held.shear_source().resolved().properties().as_ptr(),
            original_receipts
        );
        assert_eq!(
            shear.bind_maugis_equilibrium(&law, &equilibrium).unwrap(),
            held
        );
        let response = shear
            .evaluate_maugis_equilibrium(&law, &equilibrium, slip())
            .unwrap();
        common::close(response.static_capacity_n, 3e6 * state.area_m2());
        common::close(response.sliding_capacity_n.unwrap(), 2e6 * state.area_m2());
        assert!(std::ptr::eq(
            response.provenance,
            law.normal_state()
                .elastic()
                .interface()
                .interface()
                .provenance()
        ));
        if equilibrium.phase() == MaugisPhase::Contact {
            let legacy = shear.evaluate_maugis(&normal, solid, slip()).unwrap();
            assert_eq!(response.traction_n, legacy.traction_n);
            assert_eq!(response.dissipated_power_w, legacy.dissipated_power_w);
            assert!(response.dissipated_power_w > 0.0);
        } else {
            assert_eq!(state.area_m2(), 0.0);
            assert_eq!(state.repulsive_normal_force_n(), 0.0);
            assert_eq!(response.traction_n, [0.0; 3]);
            assert_eq!(response.dissipated_power_w, 0.0);
        }
    }
}

#[test]
fn complete_equilibrium_ingress_rejects_wrong_model_and_ordered_source() {
    let normal = normal(Some((3e6, 2e6)), Pressure::DIMS, false, false);
    let law = model(normal.clone(), 0.01);
    let shear =
        SourceAdhesiveShear::resolve(&normal, MaterialPropertySelection::SingleClaimOnly).unwrap();
    let contact = law.evaluate_contact_energy(0.0).unwrap();
    let mut held = shear.bind_maugis_equilibrium(&law, &contact).unwrap();
    let accepted = held.clone();
    // Even with the identical material state and zero solid area, another
    // radius belongs to a different normal model and cannot supply its token.
    let other_geometry = model(normal, 0.02);
    let other = other_geometry
        .evaluate_detached(-2.0 * other_geometry.interaction_range_m())
        .unwrap();
    assert!(SourceAdhesiveContactUpdate::from_maugis_equilibrium(&law, &other).is_err());
    assert!(shear.bind_maugis_equilibrium(&law, &other).is_err());
    assert!(
        shear
            .evaluate_maugis_equilibrium(&law, &other, slip())
            .is_err()
    );

    let reversed = model(
        super::normal(Some((3e6, 2e6)), Pressure::DIMS, false, true),
        0.01,
    );
    let other = reversed
        .evaluate_detached(-2.0 * reversed.interaction_range_m())
        .unwrap();
    let update = SourceAdhesiveContactUpdate::from_maugis_equilibrium(&reversed, &other).unwrap();
    assert!(shear.bind_maugis_equilibrium(&reversed, &other).is_err());
    assert!(
        shear
            .evaluate_maugis_equilibrium(&reversed, &other, slip())
            .is_err()
    );
    assert!(held.tangential_for_update(update).is_err());
    assert!(held.update_contact(update).is_err());
    assert_eq!(held, accepted);
}
