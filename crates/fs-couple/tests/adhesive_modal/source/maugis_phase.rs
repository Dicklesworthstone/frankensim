//! Attractive Maugis gaps release the actual native tangential consumer.
use super::*;
use core::f64::consts::PI;
use fs_contact::normal_patch::maugis::{
    MaugisEnvelope, SourceMaugisSphereContact, continuation::MaugisPhase,
};

fn maugis(temperature: f64) -> SourceMaugisSphereContact {
    let normal = super::normal(temperature).normal_state().clone();
    let radius = 0.01;
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
fn original_maugis_contact_noncontact_and_recontact_drive_native_shear_and_work() {
    let normal = maugis(293.15);
    let shear = SourceAdhesiveShear::resolve(
        normal.normal_state(),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let contact = normal.evaluate_contact_energy(0.0).unwrap();
    assert!(contact.normal_force_n() < 0.0);
    let bounds = normal.noncontact_limits().unwrap().unwrap();
    let gap = normal
        .evaluate_noncontact(
            0.5 * bounds.minimum_indentation_m + 0.5 * bounds.first_contact_indentation_m,
        )
        .unwrap();
    let detached = normal
        .evaluate_detached(-2.0 * normal.interaction_range_m())
        .unwrap();
    assert!(gap.normal_force_n() < 0.0 && gap.cohesive_radius_m() > 0.0);
    let binding = AdhesiveModalBinding::from_source(
        shear.bind_maugis_equilibrium(&normal, &contact).unwrap(),
    )
    .unwrap();
    let mut system = mechanics_with_binding(1.0, binding);
    let receipts = system
        .adhesive_binding()
        .unwrap()
        .source()
        .unwrap()
        .shear_source()
        .resolved()
        .properties()
        .as_ptr();
    let initial_energy = system.total_energy_j().unwrap();
    let initial_momentum = momentum(&system);
    let mut heat = 0.0;
    for equilibrium in [contact, gap, detached, contact] {
        let update =
            SourceAdhesiveContactUpdate::from_maugis_equilibrium(&normal, &equilibrium).unwrap();
        let solid = update.contact();
        let before = system.total_energy_j().unwrap();
        let frame = *system
            .step_with_adhesive_contact_under_gate(
                &[0.0, 0.0],
                FrictionDrive {
                    normal_force_n: solid.repulsive_normal_force_n(),
                    ..drive()
                },
                AdhesiveModalContactUpdate::Source(update),
                &CancelGate::new_clock_free(),
            )
            .unwrap();
        if equilibrium.phase() == MaugisPhase::Contact {
            assert_eq!(frame.regime, ModalFrictionRegime::Sliding);
            close(frame.traction_n, -2e6 * solid.area_m2(), 1e-12);
            assert!(frame.friction_dissipation_j > 0.0);
        } else {
            // Both zero-area phases carry zero repulsive load; attractive
            // normal work stays with the normal owner and never becomes heat.
            assert_eq!(solid.repulsive_normal_force_n(), 0.0);
            assert_eq!(frame.regime, ModalFrictionRegime::Released);
            assert_eq!(frame.traction_n, 0.0);
            assert_eq!(frame.friction_dissipation_j, 0.0);
            close(system.total_energy_j().unwrap(), before, 1e-16);
        }
        let binding = system.adhesive_binding().unwrap();
        let source = binding.source().unwrap();
        assert_eq!(binding.tangential(), source.tangential());
        assert_eq!(source.tangential().contact(), solid);
        assert_eq!(
            source.normal_state().identity(),
            normal.normal_state().identity()
        );
        assert_eq!(source.shear_source().resolved(), shear.resolved());
        assert_eq!(
            source.shear_source().resolved().properties().as_ptr(),
            receipts
        );
        assert_eq!(source.temperature_kelvin(), Some(293.15));
        assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
        heat += frame.friction_dissipation_j;
        close(momentum(&system), initial_momentum, 1e-14);
    }
    assert_eq!(system.samples_rendered(), 4);
    close(
        initial_energy - system.total_energy_j().unwrap(),
        heat,
        1e-13,
    );
}

#[test]
fn foreign_zero_area_equilibrium_cannot_release_an_original_native_junction() {
    let normal = maugis(293.15);
    let shear = SourceAdhesiveShear::resolve(
        normal.normal_state(),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let contact = normal.evaluate_contact_energy(0.0).unwrap();
    let binding = AdhesiveModalBinding::from_source(
        shear.bind_maugis_equilibrium(&normal, &contact).unwrap(),
    )
    .unwrap();
    let mut system = mechanics_with_binding(1.0, binding);
    system
        .step(
            &[0.0, 0.0],
            FrictionDrive {
                normal_force_n: contact
                    .friction_contact_state()
                    .unwrap()
                    .repulsive_normal_force_n(),
                ..drive()
            },
        )
        .unwrap();
    let before = states(&system);
    let last = system.last_frame().copied();
    let sample = system.samples_rendered();
    let binding = system.adhesive_binding().unwrap().clone();
    let foreign = maugis(300.0);
    let bounds = foreign.noncontact_limits().unwrap().unwrap();
    let gap = foreign
        .evaluate_noncontact(
            0.5 * bounds.minimum_indentation_m + 0.5 * bounds.first_contact_indentation_m,
        )
        .unwrap();
    assert!(gap.normal_force_n() < 0.0);
    let update = SourceAdhesiveContactUpdate::from_maugis_equilibrium(&foreign, &gap).unwrap();
    assert_eq!(update.contact().area_m2(), 0.0);
    assert!(
        system
            .step_with_adhesive_contact(
                &[0.0, 0.0],
                drive(),
                AdhesiveModalContactUpdate::Source(update),
            )
            .is_err()
    );
    assert_eq!(states(&system), before);
    assert_eq!(system.last_frame().copied(), last);
    assert_eq!(system.samples_rendered(), sample);
    assert_eq!(system.adhesive_binding().unwrap(), &binding);
}
