//! Directional contact after thermal fracture keeps the original matter and
//! eigenstrain. Manufactured friction parameters confer no material authority.

#[path = "friction/sustained.rs"]
mod sustained;

use super::*;
use fs_contact::frictional_multicontact::{FrictionalContactConfig, FrictionalContactDeterminacy};
use fs_solid::brittle::solid::expansion::dynamics::fracture::contact::friction::{
    BrittleFaceFriction, BrittleFaceFrictionControls, EllipticFriction,
};

fn controls() -> BrittleFaceFrictionControls {
    let normal = impact_controls();
    BrittleFaceFrictionControls {
        contact: FrictionalContactConfig {
            velocity_tolerance: 1e-11,
            impulse_tolerance: 1e-11,
            complementarity_tolerance: 1e-11,
            energy_tolerance: 1e-11,
            ..FrictionalContactConfig::default()
        },
        impulse: normal.impulse,
        activation_tolerance_m: normal.activation_tolerance_m,
        tangential_tolerance_m: normal.tangential_tolerance_m,
        velocity_tolerance_mps: normal.velocity_tolerance_mps,
        energy_tolerance_j: normal.energy_tolerance_j,
        maximum_contact_rows: normal.maximum_contact_rows,
    }
}

#[test]
fn directional_thermal_fragment_impact_preserves_source_eigenstrain_and_original_mass() {
    approach::with_approaching_motion(8, |moving, cx| {
        // The existing journey has already cooled rho(T) from 6 to 3,
        // fractured, split, opened and reversed with the original rho=6 mass.
        // Two prescribed native midpoints arrive with relative tangential
        // speed: midway tangential offsets +/-t*dt/4 give final velocities
        // drift-/+t. Halving the gap first stops normal approach; the next
        // half-gap restores +/-SPEED at exactly paired touching positions.
        let t = [SPEED / 2.0, -SPEED / 4.0];
        let mut midway = moving.solid().solution().displacement_m.clone();
        let mut touching = midway.clone();
        for (element, tet) in moving.solid().geometry().tetrahedra().iter().enumerate() {
            let sign = if element == 0 { 1.0 } else { -1.0 };
            for &node in tet {
                for c in 0..2 {
                    midway[node][c] += DRIFT[c] * INTERVAL + sign * t[c] * INTERVAL / 4.0;
                    touching[node][c] += 2.0 * DRIFT[c] * INTERVAL;
                }
                midway[node][2] += sign * SPEED * INTERVAL / 2.0;
                touching[node][2] = 0.0;
            }
        }
        let halfway = moving
            .advance_midpoint(
                moving.solid().materials(),
                INTERVAL,
                &prescribed(&midway),
                &[],
                settings(),
                midpoint_controls(),
                cx,
            )
            .unwrap()
            .motion;
        let arrived = halfway
            .advance_midpoint(
                halfway.solid().materials(),
                INTERVAL,
                &prescribed(&touching),
                &[],
                settings(),
                midpoint_controls(),
                cx,
            )
            .unwrap()
            .motion;
        for (element, tet) in arrived.solid().geometry().tetrahedra().iter().enumerate() {
            let sign = if element == 0 { 1.0 } else { -1.0 };
            for &node in tet {
                let expected = [DRIFT[0] - sign * t[0], DRIFT[1] - sign * t[1], sign * SPEED];
                for c in 0..3 {
                    close(arrived.velocity_mps()[node][c], expected[c], 1e-12);
                }
            }
        }
        assert_eq!(arrived.accepted_steps(), 5);
        assert_eq!(arrived.remaining_steps(), 3);
        approach::retained_thermal_matter(moving, &arrived);
        let friction = [BrittleFaceFriction {
            tangent_direction: [1.0, 0.0, 0.0],
            vertices: [EllipticFriction {
                coefficients_mu: [0.0, 1.0],
                principal_axis: [1.0, 0.0],
            }; 3],
        }];
        let hit = arrived
            .impact_released_faces_frictional(&[0.5], &friction, &[], ORIGIN, controls(), cx)
            .unwrap();
        assert_eq!(hit.rows.len(), 3);
        assert_eq!(hit.motion.free_linear_strain_rates_per_s(), arrived.free_linear_strain_rates_per_s());
        assert!(matches!(
            hit.contact.transformed.determinacy,
            FrictionalContactDeterminacy::UniqueImpulses
        ));
        // Each original tet has mass1 despite the current source rho=3.
        // M=(I+ones)/20 gives uniform paired-face mobility 16. Physical X
        // remains exactly free; physical Y sticks with py=-t_y/8.
        let p = [0.0, SPEED / 32.0, 3.0 * SPEED / 16.0];
        for (&impulse, &velocity) in hit
            .contact
            .impulses
            .iter()
            .zip(&hit.contact.relative_velocities_after)
        {
            for c in 0..3 {
                close(impulse[c], [p[2], p[0], p[1]][c], 1e-10);
                close(velocity[c], [SPEED, SPEED, 0.0][c], 2e-9);
            }
        }
        for (element, tet) in arrived.solid().geometry().tetrahedra().iter().enumerate() {
            let sign = if element == 0 { -1.0 } else { 1.0 };
            for &node in tet {
                let is_face = hit.rows.iter().any(|row| row.nodes[element] == node);
                let response = if is_face { 8.0 } else { -12.0 };
                for c in 0..3 {
                    close(
                        hit.motion.velocity_mps()[node][c],
                        arrived.velocity_mps()[node][c] + sign * response * p[c],
                        2e-9,
                    );
                }
            }
        }
        let old = inventory(&arrived);
        let new = inventory(&hit.motion);
        let loss = 39.0 * SPEED * SPEED / 128.0;
        for c in 0..3 {
            close(new.0[c], old.0[c], 2e-10);
            close(new.1[c], old.1[c], 2e-10);
        }
        close(hit.contact.dissipated_energy, loss, 2e-10);
        close(hit.contact.relative_work, -loss, 2e-10);
        close(hit.relative_work_j, -loss, 2e-10);
        close(hit.contact_work_error_j, 0.0, 2e-10);
        close(new.2 - old.2, -loss, 2e-10);
        close(hit.motion.kinetic_energy_j(), new.2, 2e-12);
        close(hit.impulse.kinetic_change_j, -loss, 2e-10);
        close(hit.impulse.support_work_j, 0.0, 2e-12);
        close(hit.relative_work_error_j, 0.0, 2e-10);
        close(hit.energy_balance_error_j, 0.0, 2e-10);
        close(hit.absolute_contact_couple_kg_m2_s, 0.0, 2e-12);

        approach::retained_thermal_matter(&arrived, &hit.motion);
        let before = arrived.solid().solution();
        let after = hit.motion.solid().solution();
        assert_eq!(after.displacement_m, before.displacement_m);
        assert_eq!(after.element_strain, before.element_strain);
        assert_eq!(after.element_stress, before.element_stress);
        assert_eq!(after.element_states, before.element_states);
        assert_eq!(hit.motion.reference_mass_kg(), arrived.reference_mass_kg());
        assert_eq!(hit.motion.time_s().to_bits(), arrived.time_s().to_bits());
        assert_eq!(hit.motion.accepted_steps(), arrived.accepted_steps());
        assert_eq!(hit.motion.remaining_steps(), arrived.remaining_steps());
        assert_eq!(
            hit.motion.solid().stored_energy_j(),
            arrived.solid().stored_energy_j()
        );
        assert_eq!(
            hit.motion.solid().fracture_dissipation_j(),
            arrived.solid().fracture_dissipation_j()
        );
        for (a, b) in arrived
            .solid()
            .materials()
            .iter()
            .zip(hit.motion.solid().materials())
        {
            assert_eq!(
                a.free_linear_strain().to_bits(),
                b.free_linear_strain().to_bits()
            );
            assert_eq!(a.expansion(), b.expansion());
            close(b.free_linear_strain(), -0.01, 1e-15);
            close(b.band().reference_density_kg_m3().unwrap(), 3.0, 1e-15);
        }

        let gate = CancelGate::new_clock_free();
        gate.request();
        ArenaPool::new(ArenaConfig::default()).scope(|arena| {
            let cancelled = Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 598,
                    kernel_id: 598,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            );
            assert!(
                arrived
                    .impact_released_faces_frictional(
                        &[0.5],
                        &friction,
                        &[],
                        ORIGIN,
                        controls(),
                        &cancelled
                    )
                    .is_err()
            );
        });
        let retry = arrived
            .impact_released_faces_frictional(&[0.5], &friction, &[], ORIGIN, controls(), cx)
            .unwrap();
        assert_eq!(retry.contact, hit.contact);
        assert_eq!(retry.motion.velocity_mps(), hit.motion.velocity_mps());
        let continued = hit
            .motion
            .advance_midpoint(
                hit.motion.solid().materials(),
                1.0 / 4096.0,
                &[],
                &[],
                settings(),
                midpoint_controls(),
                cx,
            )
            .unwrap();
        assert_eq!(continued.motion.accepted_steps(), 6);
        assert_eq!(continued.motion.remaining_steps(), 2);
        assert!(core::ptr::eq(
            continued.motion.reference_density_kg_m3().as_ptr(),
            moving.reference_density_kg_m3().as_ptr()
        ));
        assert!(core::ptr::eq(
            continued.motion.solid().materials().as_ptr(),
            moving.solid().materials().as_ptr()
        ));
        close(continued.motion.solid().fracture_dissipation_j(), GC, 1e-13);
        close(
            continued.motion.kinetic_energy_j() + continued.motion.solid().stored_energy_j(),
            hit.motion.kinetic_energy_j() + hit.motion.solid().stored_energy_j(),
            1e-10,
        );
        close(continued.energy.total_balance_error_j, 0.0, 1e-10);
    });
}
