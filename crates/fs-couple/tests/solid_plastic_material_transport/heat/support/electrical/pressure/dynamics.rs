//! Inertia and both sourced interface laws use one accepted support state.
use super::*;
use fs_couple::solid_enthalpy::plastic::dynamics::ThermoplasticDynamics;

#[test]
fn retained_motion_reaches_pressure_dependent_current_and_thermal_supports() {
    with_current_support(
        0.01,
        8e6,
        0.1,
        0.3,
        |cx, _, t, s, support, electrical, p, c| {
        let mut velocity = vec![[0.0; 3]; s.solution().displacement_m.len()];
            velocity[3][0] = 2.0;
            let moving = ThermoplasticDynamics::try_new(s.clone(), velocity, cx).unwrap();
            let next = moving
                .advance_electrical_with_pressure_contacts_candidate(
                    t,
                    c,
                    p,
                    std::slice::from_ref(support),
                    &[electrical],
                    &[],
                    2e-6,
                    &displacements(0.001, 0.0),
                    &[],
                    IncrementSettings::default(),
                    powered_controls(),
                    &[0.0; 2],
                    cx,
                )
                .unwrap();
            let r = next.report();
            let contact = &r.material.support_contacts[0];
            close(contact.compression_pa, 2e6, 1e-4);
            assert_eq!(contact.resistance.receipt().selected, support.claim);
            assert!(contact.bulk_conductance_w_k.is_some());
            assert_eq!(
                r.material.circuit_thermal.source.destinations[0].heat_fraction_to_first,
                0.3
            );
            assert!(r.material.circuit_thermal.electrical.contacts[0].heat_j > 0.0);
            assert!(next.dynamics().solid().solution().displacement_m[3][0] > 0.0);
            assert!(r.inertia.kinetic_before_j > 0.0 && r.inertia.numerical_loss_j > 0.0);
            close(r.material.physical_plastic_heat_j, 0.0, 1e-12);
            close(
                r.material.circuit_thermal.thermal.stored_energy_change_j,
                r.material.circuit_thermal.deposited_joule_j,
                1e-9,
            );
            close(
                r.material.combined_residual_j,
                0.0,
                powered_controls().coupled.total_energy_tolerance_j,
            );
            assert_eq!(moving.time_s(), 0.0);
            assert_eq!(c.charges_c(), &[0.1]);
            assert_eq!(t.time_s(), 0.0);
        },
    );
}
