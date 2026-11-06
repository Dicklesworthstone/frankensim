//! Oblique post-fracture impacts against exact P1 consistent-mass oracles.
//! Friction coefficients are explicit manufactured models, not material data.

#[path = "friction/source.rs"]
mod source;
#[path = "friction/sustained.rs"]
mod sustained;
#[path = "friction/time.rs"]
mod time;

use super::*;
use fs_contact::frictional_multicontact::{FrictionalContactConfig, FrictionalContactDeterminacy};
use fs_solid::brittle::solid::dynamics::fracture::contact::friction::{
    BrittleFaceFriction, BrittleFaceFrictionControls, BrittleFaceFrictionStep, EllipticFriction,
};

fn controls() -> BrittleFaceFrictionControls {
    let base = impact_controls();
    let contact = FrictionalContactConfig {
        velocity_tolerance: 1e-11,
        impulse_tolerance: 1e-11,
        complementarity_tolerance: 1e-11,
        energy_tolerance: 1e-11,
        ..FrictionalContactConfig::default()
    };
    BrittleFaceFrictionControls {
        contact,
        impulse: base.impulse,
        activation_tolerance_m: base.activation_tolerance_m,
        tangential_tolerance_m: base.tangential_tolerance_m,
        velocity_tolerance_mps: 1e-9,
        energy_tolerance_j: base.energy_tolerance_j,
        maximum_contact_rows: base.maximum_contact_rows,
    }
}

fn law(mu: [f64; 2], tangent: [f64; 3], principal_axis: [f64; 2]) -> BrittleFaceFriction {
    BrittleFaceFriction {
        tangent_direction: tangent,
        vertices: [EllipticFriction {
            coefficients_mu: mu,
            principal_axis,
        }; 3],
    }
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    (0..3).map(|c| a[c] * b[c]).sum()
}

// Reuse the original source bundle and fracture law, preparing its geometric
// crack normal for the requested specimen orientation. No detached mass,
// constitutive substitute, or artificially failed state enters this journey.
fn with_oblique_motion(
    nodes: &[[f64; 3]; 5],
    normal: [f64; 3],
    velocities: [[f64; 3]; 2],
    run: impl FnOnce(&BrittleTetMotion<'_>, &Cx<'_>),
) {
    with_oblique_motion_budget(nodes, normal, velocities, 2, run);
}

fn with_oblique_motion_budget(
    nodes: &[[f64; 3]; 5],
    normal: [f64; 3],
    velocities: [[f64; 3]; 2],
    maximum_steps: usize,
    run: impl FnOnce(&BrittleTetMotion<'_>, &Cx<'_>),
) {
    with_cx(false, |cx| {
        let geometry = TetLoadGeometry3::new(nodes, &TETS, budget(), 8, cx).unwrap();
        let materials: Vec<_> = sources(&geometry)
            .iter()
            .enumerate()
            .map(|(element, source)| {
                let band = TetCrackBandGeometry::from_element(&geometry, element, normal).unwrap();
                ResolvedBrittleCrackBand::try_from_material(
                    source.material().clone(),
                    band.material_settings(0.03, 1e-9),
                )
                .unwrap()
            })
            .collect();
        let displacement: Vec<_> = nodes
            .iter()
            .map(|&x| normal.map(|n| n * dot(x, normal) / 128.0))
            .collect();
        let solid = BrittleTetSolid::try_new(&geometry, &materials, budget(), cx)
            .unwrap()
            .advance(
                &prescribed(&displacement),
                &[],
                settings(),
                work_controls(),
                cx,
            )
            .unwrap()
            .solid;
        assert!(solid.solution().element_states.iter().all(|s| s.failed()));
        assert_eq!(solid.stored_energy_j(), 0.0);
        let cut = topology(
            &geometry,
            &[InteriorSeam {
                id: 951,
                nodes: [0, 1, 2],
                side_a_element: 0,
            }],
            cx,
        );
        let target =
            TetLoadGeometry3::new(cut.nodes_m(), cut.tetrahedra(), budget(), 8, cx).unwrap();
        let split = solid
            .split_fracture(&cut, &target, fracture_controls(), cx)
            .unwrap();
        let common: [f64; 3] =
            core::array::from_fn(|c| (velocities[0][c] + velocities[1][c]) / 2.0);
        let mut opened = split.solid.solution().displacement_m.clone();
        let mut endpoint = opened.clone();
        let mut velocity = vec![[0.0; 3]; opened.len()];
        for (element, tet) in target.tetrahedra().iter().enumerate() {
            for &node in tet {
                velocity[node] = velocities[element];
                for c in 0..3 {
                    endpoint[node][c] += common[c] * APPROACH_DT;
                    opened[node][c] = endpoint[node][c] - velocity[node][c] * APPROACH_DT;
                }
            }
        }
        let open = split
            .solid
            .advance(&prescribed(&opened), &[], settings(), work_controls(), cx)
            .unwrap()
            .solid;
        let motion =
            BrittleTetMotion::try_from_solid(&open, &velocity, 2.0, maximum_steps, cx).unwrap();
        let motion = motion
            .advance_midpoint(
                APPROACH_DT,
                &prescribed(&endpoint),
                &[],
                settings(),
                midpoint_controls(),
                cx,
            )
            .unwrap()
            .motion;
        for (&actual, &expected) in motion.velocity_mps().iter().zip(&velocity) {
            for c in 0..3 {
                close(actual[c], expected[c], 2e-14);
            }
        }
        assert_eq!(motion.accepted_steps(), 1);
        assert_eq!(motion.remaining_steps(), maximum_steps - 1);
        run(&motion, cx);
    });
}

// For each unit-mass tet M=(I+ones)/20, M^-1=20I-4ones. Three equal
// face impulses therefore change a face velocity by 8*p and the opposing
// apex velocity by -12*p. Held Cartesian coordinates have zero increment.
fn check_nodal_response(
    before: &BrittleTetMotion<'_>,
    hit: &BrittleFaceFrictionStep<'_>,
    impulse_on_b: [f64; 3],
    held_a: [bool; 3],
) {
    for (element, tet) in before.solid().geometry().tetrahedra().iter().enumerate() {
        let sign = if element == 0 { -1.0 } else { 1.0 };
        for &node in tet {
            let is_face = hit.rows.iter().any(|row| row.nodes[element] == node);
            let response = if is_face { 8.0 } else { -12.0 };
            for c in 0..3 {
                let delta = if element == 0 && held_a[c] {
                    0.0
                } else {
                    sign * response * impulse_on_b[c]
                };
                close(
                    hit.motion.velocity_mps()[node][c],
                    before.velocity_mps()[node][c] + delta,
                    2e-9,
                );
            }
        }
    }
}

fn check_balances(
    before: &BrittleTetMotion<'_>,
    hit: &BrittleFaceFrictionStep<'_>,
    expected_loss: f64,
    support_per_row: [f64; 3],
) {
    retained(before, &hit.motion);
    assert_eq!(hit.rows.len(), 3);
    assert!(matches!(
        hit.contact.transformed.determinacy,
        FrictionalContactDeterminacy::UniqueImpulses
    ));
    let old = inventory(before);
    let new = inventory(&hit.motion);
    let mut support_work = 0.0;
    let mut support_moment = [0.0; 3];
    for row in &hit.rows {
        let node = row.nodes[0];
        let position = core::array::from_fn(|c| {
            before.solid().geometry().nodes_m()[node][c]
                + before.solid().solution().displacement_m[node][c]
                - ORIGIN[c]
        });
        support_work += dot(support_per_row, before.velocity_mps()[node]);
        let moment = cross(position, support_per_row);
        for c in 0..3 {
            support_moment[c] += moment[c];
        }
    }
    for c in 0..3 {
        close(new.0[c] - old.0[c], 3.0 * support_per_row[c], 2e-10);
        close(new.1[c] - old.1[c], support_moment[c], 2e-10);
    }
    close(hit.contact.dissipated_energy, expected_loss, 2e-10);
    close(hit.contact.relative_work, -expected_loss, 2e-10);
    close(new.2 - old.2, support_work - expected_loss, 2e-10);
    close(hit.motion.kinetic_energy_j(), new.2, 2e-12);
    close(hit.impulse.kinetic_change_j, new.2 - old.2, 2e-12);
    close(hit.impulse.support_work_j, support_work, 2e-10);
    close(hit.relative_work_error_j, 0.0, 2e-10);
    close(hit.energy_balance_error_j, 0.0, 2e-10);
    close(hit.absolute_contact_couple_kg_m2_s, 0.0, 2e-12);
}

#[test]
fn oblique_impact_has_exact_frictionless_sliding_and_sticking_mass_responses() {
    let t = [SPEED / 2.0, -SPEED / 4.0];
    let velocities = [
        [DRIFT[0] - t[0], DRIFT[1] - t[1], SPEED],
        [DRIFT[0] + t[0], DRIFT[1] + t[1], -SPEED],
    ];
    with_oblique_motion(&NODES, [0.0, 0.0, 1.0], velocities, |motion, cx| {
        let e = 0.5;
        for mu in [0.0, 0.125, 1.0] {
            let friction = [law([mu; 2], [1.0, 0.0, 0.0], [1.0, 0.0])];
            let hit = motion
                .impact_released_faces_frictional(&[e], &friction, &[], ORIGIN, controls(), cx)
                .unwrap();
            // The paired-face mobility is 40I-8ones, uniform eigenvalue 16.
            let pn = (1.0 + e) * SPEED / 8.0;
            let tangent_norm = t[0].hypot(t[1]);
            let traction = (tangent_norm / 8.0).min(mu * pn);
            let pt = t.map(|v| -traction * v / tangent_norm);
            for (&p, &v) in hit
                .contact
                .impulses
                .iter()
                .zip(&hit.contact.relative_velocities_after)
            {
                for c in 0..3 {
                    close(p[c], [pn, pt[0], pt[1]][c], 1e-10);
                    close(
                        v[c],
                        [
                            2.0 * e * SPEED,
                            2.0 * t[0] + 16.0 * pt[0],
                            2.0 * t[1] + 16.0 * pt[1],
                        ][c],
                        2e-9,
                    );
                }
            }
            let loss = 3.0 * (1.0 - e * e) * SPEED * SPEED / 8.0
                - 3.0
                    * (2.0 * (t[0] * pt[0] + t[1] * pt[1]) + 8.0 * (pt[0] * pt[0] + pt[1] * pt[1]));
            check_nodal_response(motion, &hit, [pt[0], pt[1], pn], [false; 3]);
            check_balances(motion, &hit, loss, [0.0; 3]);
            if mu == 0.0 {
                let normal = motion
                    .impact_released_faces(&[e], &[], ORIGIN, impact_controls(), cx)
                    .unwrap();
                for (a, b) in hit
                    .motion
                    .velocity_mps()
                    .iter()
                    .zip(normal.motion.velocity_mps())
                {
                    for c in 0..3 {
                        close(a[c], b[c], 2e-9);
                    }
                }
            }
            if mu == 1.0 {
                // Refusals leave this same used-clock specimen retryable.
                let mut bounded = controls();
                bounded.maximum_contact_rows = 2;
                assert!(
                    motion
                        .impact_released_faces_frictional(&[e], &friction, &[], ORIGIN, bounded, cx)
                        .is_err()
                );
                with_cx(true, |cancelled| {
                    assert!(
                        motion
                            .impact_released_faces_frictional(
                                &[e],
                                &friction,
                                &[],
                                ORIGIN,
                                controls(),
                                cancelled,
                            )
                            .is_err()
                    );
                });
                let retry = motion
                    .impact_released_faces_frictional(&[e], &friction, &[], ORIGIN, controls(), cx)
                    .unwrap();
                assert_eq!(retry.motion.velocity_mps(), hit.motion.velocity_mps());
                assert_eq!(retry.contact, hit.contact);
                retained(motion, &retry.motion);
            }
        }
    });
}

#[test]
fn directional_free_slip_and_sticking_are_covariant_under_contact_frame_rotation() {
    let velocities = [[-SPEED, -SPEED / 4.0, SPEED], [SPEED, SPEED / 4.0, -SPEED]];
    with_oblique_motion(&NODES, [0.0, 0.0, 1.0], velocities, |motion, cx| {
        let e = 0.5;
        let pn = (1.0 + e) * SPEED / 8.0;
        let py = -SPEED / 32.0;
        let expected_loss =
            3.0 * (1.0 - e * e) * SPEED * SPEED / 8.0 + 3.0 * (SPEED / 4.0).powi(2) / 8.0;
        let mut previous: Option<Vec<[f64; 3]>> = None;
        for friction in [
            law([0.0, 1.0], [1.0, 0.0, 0.0], [1.0, 0.0]),
            // Project away the normal component; the declared texture stays
            // along physical X after the contact basis rotates by 45 degrees.
            law([0.0, 1.0], [1.0, 1.0, 7.0], [2.0, -2.0]),
        ] {
            let hit = motion
                .impact_released_faces_frictional(&[e], &[friction], &[], ORIGIN, controls(), cx)
                .unwrap();
            for (row, frame) in hit.frames.iter().enumerate() {
                let p = hit.contact.impulses[row];
                let v = hit.contact.relative_velocities_after[row];
                for c in 0..3 {
                    let physical_p = p[0] * frame.normal[c]
                        + p[1] * frame.tangent1[c]
                        + p[2] * frame.tangent2[c];
                    let physical_v = v[0] * frame.normal[c]
                        + v[1] * frame.tangent1[c]
                        + v[2] * frame.tangent2[c];
                    close(physical_p, [0.0, py, pn][c], 1e-10);
                    close(physical_v, [2.0 * SPEED, 0.0, 2.0 * e * SPEED][c], 2e-9);
                }
            }
            check_nodal_response(motion, &hit, [0.0, py, pn], [false; 3]);
            check_balances(motion, &hit, expected_loss, [0.0; 3]);
            if let Some(old) = previous {
                for (a, b) in hit.motion.velocity_mps().iter().zip(old) {
                    for c in 0..3 {
                        close(a[c], b[c], 2e-9);
                    }
                }
            }
            previous = Some(hit.motion.velocity_mps().to_vec());
        }
    });
}

#[test]
fn moving_held_fragment_has_separate_friction_loss_and_support_work() {
    let a = [2.0 * SPEED, -SPEED, -SPEED / 2.0];
    let b = [a[0] + SPEED, a[1] - SPEED / 2.0, a[2] - 2.0 * SPEED];
    with_oblique_motion(&NODES, [0.0, 0.0, 1.0], [a, b], |motion, cx| {
        let held: Vec<_> = motion.solid().geometry().tetrahedra()[0]
            .iter()
            .flat_map(|&node| (0..3).map(move |c| 3 * node + c))
            .collect();
        let hit = motion
            .impact_released_faces_frictional(
                &[0.0],
                &[law([2.0; 2], [1.0, 0.0, 0.0], [1.0, 0.0])],
                &held,
                ORIGIN,
                controls(),
                cx,
            )
            .unwrap();
        // Only B is mobile: the uniform contact mobility is 8, not 16.
        let p = [-SPEED / 8.0, SPEED / 16.0, SPEED / 4.0];
        let loss = 63.0 * SPEED * SPEED / 64.0;
        check_nodal_response(motion, &hit, p, [true; 3]);
        check_balances(motion, &hit, loss, p);
        close(
            hit.impulse.support_work_j,
            -21.0 * SPEED * SPEED / 16.0,
            2e-10,
        );
        close(
            hit.impulse.kinetic_change_j,
            -147.0 * SPEED * SPEED / 64.0,
            2e-10,
        );
        for v in &hit.contact.relative_velocities_after {
            for &component in v {
                close(component, 0.0, 2e-9);
            }
        }
    });
}

#[test]
fn rotated_partial_support_solves_normal_and_directional_friction_together() {
    let q = core::f64::consts::FRAC_1_SQRT_2;
    let normal = [q, 0.0, q];
    let tangent = [q, 0.0, -q];
    let nodes = NODES.map(|x| [q * (x[0] + x[2]), x[1], q * (x[2] - x[0])]);
    let a = core::array::from_fn(|c| SPEED * (normal[c] - tangent[c]));
    let b = a.map(|v| -v);
    with_oblique_motion(&nodes, normal, [a, b], |motion, cx| {
        let held: Vec<_> = motion.solid().geometry().tetrahedra()[0]
            .iter()
            .map(|&node| 3 * node + 2)
            .collect();
        let hit = motion
            .impact_released_faces_frictional(
                &[0.5],
                &[law([0.25, 0.5], tangent, [1.0, 0.0])],
                &held,
                ORIGIN,
                controls(),
                cx,
            )
            .unwrap();
        // Constraining global Z on A leaves uniform physical mobility
        // diag(16,16,8). In the rotated [n,t1,t2] frame it is
        // [[12,4,0],[4,12,0],[0,0,16]]. Sliding gives pt=-pn/4, so
        // 12*pn+4*pt=3*SPEED => pn=3*SPEED/11. Solving the normal
        // problem first and clipping friction would incorrectly use /12.
        let pn = 3.0 * SPEED / 11.0;
        let pt = -pn / 4.0;
        let physical = core::array::from_fn(|c| pn * normal[c] + pt * tangent[c]);
        for (&p, &v) in hit
            .contact
            .impulses
            .iter()
            .zip(&hit.contact.relative_velocities_after)
        {
            for c in 0..3 {
                close(p[c], [pn, pt, 0.0][c], 1e-10);
                close(v[c], [SPEED, 2.0 * SPEED + pn, 0.0][c], 2e-9);
            }
        }
        let loss = 3.0 * (SPEED * pn + pn * pn / 8.0);
        check_nodal_response(motion, &hit, physical, [false, false, true]);
        check_balances(motion, &hit, loss, [0.0, 0.0, physical[2]]);
        close(hit.impulse.support_work_j, 15.0 * pn * SPEED / 4.0, 2e-10);
        // The moving constraint supplies more work than friction dissipates;
        // increased total kinetic energy must not be mislabeled active friction.
        assert!(hit.impulse.kinetic_change_j > 0.0);
        assert!(hit.contact.dissipated_energy > 0.0);
    });
}
