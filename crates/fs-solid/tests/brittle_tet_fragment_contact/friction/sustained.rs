//! Sustained directional sticking uses original force-driven midpoint fields.
//! All coefficients are manufactured models. Paired witnesses never authorize
//! a finite sliding trajectory; an overloaded static interface must refuse.

use super::*;
use fs_contact::sustained::solve::{SolveControls, SolveError};
use fs_solid::brittle::solid::dynamics::fracture::contact::friction::sustained::{
    BrittleFaceSustainedFrictionControls, BrittleFaceSustainedFrictionError,
    BrittleFaceSustainedFrictionStep,
};
use fs_solid::brittle::solid::dynamics::fracture::contact::sustained::{
    BrittleFaceSustainedControls, BrittleFaceSustainedError,
};
use fs_solver::NewtonKrylovConfig;

const H: f64 = 1.0 / 64.0;
const VELOCITY: [f64; 3] = [SPEED, -SPEED / 2.0, SPEED / 4.0];

pub(super) fn sustained_controls() -> BrittleFaceSustainedFrictionControls {
    BrittleFaceSustainedFrictionControls {
        contact: BrittleFaceSustainedControls {
            solve: SolveControls {
                newton: NewtonKrylovConfig {
                    absolute_tolerance: 1e-11,
                    relative_tolerance: 0.0,
                    ..NewtonKrylovConfig::default()
                },
                maximum_response_evaluations: 512,
                difference_step_n: 1e-3,
                force_residual_tolerance_n: 1e-9,
                maximum_contacts: 3,
                ..SolveControls::default()
            },
            work: midpoint_controls(),
            normal_scale_n_per_m: 4096.0,
            activation_tolerance_m: 1e-10,
            tangential_tolerance_m: 1e-10,
            velocity_tolerance_mps: 1e-8,
            force_tolerance_n: 1e-8,
            complementarity_tolerance_j: 1e-10,
            contact_work_tolerance_j: 1e-10,
            force_moment_tolerance_n_m: 1e-10,
        },
        tangent_scale_n_per_m: 4096.0,
        slip_tolerance_m: 1e-10,
    }
}

fn with_touching_motion(run: impl FnOnce(&BrittleTetMotion<'_>, &Cx<'_>)) {
    with_cx(false, |cx| {
        let geometry = TetLoadGeometry3::new(&NODES, &TETS, budget(), 8, cx).unwrap();
        let materials = sources(&geometry);
        let solid = failed(&geometry, &materials, cx);
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
        let velocities = vec![VELOCITY; target.nodes_m().len()];
        let motion =
            BrittleTetMotion::try_from_solid(&split.solid, &velocities, 2.0, 3, cx).unwrap();
        run(&motion, cx);
    });
}

fn body_forces(motion: &BrittleTetMotion<'_>, acceleration_a: [f64; 3]) -> Vec<(usize, f64)> {
    let mut forces = Vec::new();
    for (element, tet) in motion.solid().geometry().tetrahedra().iter().enumerate() {
        let sign = if element == 0 { 1.0 } else { -1.0 };
        // Integral of each shape function is V/4; rho*V=1 for each tet.
        for &node in tet {
            for (c, acceleration) in acceleration_a.iter().enumerate() {
                forces.push((3 * node + c, sign * acceleration / 4.0));
            }
        }
    }
    forces
}

fn supports(motion: &BrittleTetMotion<'_>, h: f64, held: [bool; 3]) -> Vec<(usize, f64)> {
    let mut endpoints = Vec::new();
    for &node in &motion.solid().geometry().tetrahedra()[0] {
        for c in 0..3 {
            if held[c] {
                endpoints.push((
                    3 * node + c,
                    motion.solid().solution().displacement_m[node][c]
                        + h * motion.velocity_mps()[node][c],
                ));
            }
        }
    }
    endpoints.sort_unstable_by_key(|&(dof, _)| dof);
    endpoints
}

fn audit_interval(
    before: &BrittleTetMotion<'_>,
    step: &BrittleFaceSustainedFrictionStep<'_>,
    h: f64,
    applied: &[(usize, f64)],
    support_forces: &[(usize, f64)],
) {
    // Integrate constant physical external forces at actual midpoint positions,
    // independently of both contact assembly and native mass-matrix reports.
    let mut force = vec![[0.0; 3]; before.velocity_mps().len()];
    for &(dof, value) in applied.iter().chain(support_forces) {
        force[dof / 3][dof % 3] += value;
    }
    let mut expected_p = [0.0; 3];
    let mut expected_l = [0.0; 3];
    for (node, f) in force.iter().enumerate() {
        let x = core::array::from_fn(|c| {
            before.solid().geometry().nodes_m()[node][c] - ORIGIN[c]
                + 0.5
                    * (before.solid().solution().displacement_m[node][c]
                        + step.motion.solid().solution().displacement_m[node][c])
        });
        let moment = cross(x, *f);
        for c in 0..3 {
            expected_p[c] += h * f[c];
            expected_l[c] += h * moment[c];
        }
    }
    let old = inventory(before);
    let next = inventory(&step.motion);
    for c in 0..3 {
        close(next.0[c] - old.0[c], expected_p[c], 2e-9);
        close(next.1[c] - old.1[c], expected_l[c], 2e-9);
    }
    close(step.motion.kinetic_energy_j(), next.2, 2e-12);
    close(step.energy.inertia.kinetic_change_j, next.2 - old.2, 2e-12);
    close(step.energy.total_balance_error_j, 0.0, 2e-10);
    close(step.energy.stored_after_j, 0.0, 2e-10);
    close(step.normal_contact_work_j, 0.0, 2e-10);
    close(step.tangential_contact_work_j, 0.0, 2e-10);
    assert!(step.absolute_contact_work_j <= 2e-10);
    close(step.coordinate_work_error_j, 0.0, 2e-10);
    close(step.contact_work_balance_error_j, 0.0, 2e-10);
    close(step.absolute_contact_couple_n_m, 0.0, 2e-10);
    close(
        step.native_energy.dead_load_mean_work_j,
        step.energy.dead_load_mean_work_j
            + step.normal_contact_work_j
            + step.tangential_contact_work_j,
        2e-10,
    );
    assert_eq!(step.energy.inertia.numerical_loss_j, 0.0);
    assert_eq!(step.energy.fracture_dissipation_j, 0.0);
    assert!(core::ptr::eq(
        before.solid().materials().as_ptr(),
        step.motion.solid().materials().as_ptr()
    ));
    assert!(core::ptr::eq(
        before.reference_density_kg_m3().as_ptr(),
        step.motion.reference_density_kg_m3().as_ptr()
    ));
    assert_eq!(step.motion.reference_mass_kg(), before.reference_mass_kg());
    assert_eq!(
        step.motion.solid().released_faces(),
        before.solid().released_faces()
    );
    assert_eq!(
        step.motion.solid().fracture_dissipation_j(),
        before.solid().fracture_dissipation_j()
    );
    for (old, new) in before
        .solid()
        .solution()
        .element_states
        .iter()
        .zip(&step.motion.solid().solution().element_states)
    {
        assert!(new.failed());
        assert_eq!(new.damage(), old.damage());
        assert_eq!(new.dissipated_j_m3(), old.dissipated_j_m3());
        assert!(core::ptr::eq(
            new.cohesive_point().law(),
            old.cohesive_point().law()
        ));
    }
    assert_eq!(step.motion.time_s(), before.time_s() + h);
    assert_eq!(step.motion.accepted_steps(), before.accepted_steps() + 1);
    assert_eq!(step.motion.remaining_steps() + 1, before.remaining_steps());
}

#[test]
fn directional_sticking_body_loads_retain_consistent_mass_internal_motion_and_force_units() {
    with_touching_motion(|motion, cx| {
        let a = [0.25, -0.125, 1.0];
        let forces = body_forces(motion, a);
        let friction = [law([0.5, 0.25], [1.0, 0.0, 0.0], [1.0, 0.0])];
        for h in [H, H / 2.0] {
            let step = motion
                .advance_released_face_sustained_frictional(
                    h,
                    &friction,
                    &[],
                    &forces,
                    settings(),
                    sustained_controls(),
                    cx,
                )
                .unwrap();
            assert_eq!(step.rows.len(), 3);
            // M=(I+ones)/20 and M^-1=20I-4ones. A uniform paired force
            // changes relative face acceleration by 16*f. Body acceleration
            // is -2*a, so the static contact force is a/8 for BOTH durations.
            for (force, path) in step.forces_n.iter().zip(&step.paths_m) {
                for c in 0..3 {
                    close(force[c], [a[2] / 8.0, a[0] / 8.0, a[1] / 8.0][c], 2e-8);
                    close(path[c], 0.0, 1e-10);
                }
                assert!(path[0] >= 0.0);
            }
            for (element, tet) in motion.solid().geometry().tetrahedra().iter().enumerate() {
                let sign = if element == 0 { 1.0 } else { -1.0 };
                for &node in tet {
                    let is_face = step.rows.iter().any(|row| row.nodes[element] == node);
                    for c in 0..3 {
                        let acceleration = if is_face { 0.0 } else { sign * 2.5 * a[c] };
                        close(
                            step.motion.velocity_mps()[node][c],
                            VELOCITY[c] + h * acceleration,
                            2e-8,
                        );
                        close(
                            step.motion.solid().solution().displacement_m[node][c]
                                - motion.solid().solution().displacement_m[node][c],
                            h * VELOCITY[c] + 0.5 * h * h * acceleration,
                            2e-10,
                        );
                    }
                }
            }
            // Native consistent mass retains the moving apex; rigid point
            // masses would incorrectly stop both bodies and erase this work.
            let work = 5.0 * dot(a, a) * h * h / 8.0;
            close(step.energy.dead_load_mean_work_j, work, 2e-10);
            close(step.energy.support_mean_work_j, 0.0, 2e-10);
            close(step.energy.inertia.kinetic_change_j, work, 2e-10);
            audit_interval(motion, &step, h, &forces, &[]);
        }
    });
}

#[test]
fn moving_support_transfers_static_friction_work_then_unloads_without_gluing() {
    with_touching_motion(|motion, cx| {
        let p = [1.0 / 32.0, -1.0 / 64.0, 1.0 / 8.0];
        let face = &motion.solid().released_faces()[0];
        let mut forces = Vec::new();
        let mut reactions = Vec::new();
        for vertex in 0..3 {
            for c in 0..3 {
                forces.push((3 * face.paired_nodes[1][vertex] + c, -p[c]));
                reactions.push((3 * face.paired_nodes[0][vertex] + c, p[c]));
            }
        }
        let friction = [law([0.5, 0.25], [1.0, 0.0, 0.0], [1.0, 0.0])];
        let loaded = motion
            .advance_released_face_sustained_frictional(
                H,
                &friction,
                &supports(motion, H, [true; 3]),
                &forces,
                settings(),
                sustained_controls(),
                cx,
            )
            .unwrap();
        for force in &loaded.forces_n {
            for c in 0..3 {
                close(force[c], [p[2], p[0], p[1]][c], 2e-8);
            }
        }
        for node in 0..motion.velocity_mps().len() {
            for c in 0..3 {
                close(loaded.motion.velocity_mps()[node][c], VELOCITY[c], 2e-8);
                close(
                    loaded.motion.solid().solution().displacement_m[node][c]
                        - motion.solid().solution().displacement_m[node][c],
                    H * VELOCITY[c],
                    2e-10,
                );
            }
        }
        // Sticking transmits power between the supported and loaded bodies;
        // their equal signed contributions cancel without a frictional loss.
        let work = 3.0 * dot(p, VELOCITY) * H;
        assert!(work > 0.0);
        close(loaded.energy.support_mean_work_j, work, 2e-10);
        close(loaded.energy.dead_load_mean_work_j, -work, 2e-10);
        close(loaded.energy.inertia.kinetic_change_j, 0.0, 2e-10);
        audit_interval(motion, &loaded, H, &forces, &reactions);

        // Removing the support and applying outward body loads opens the face.
        // Tangential separation is valid here because all contact forces vanish.
        let a = [0.25, -0.125, -1.0];
        let separating = body_forces(&loaded.motion, a);
        let released = loaded
            .motion
            .advance_released_face_sustained_frictional(
                H,
                &friction,
                &[],
                &separating,
                settings(),
                sustained_controls(),
                cx,
            )
            .unwrap();
        for (force, path) in released.forces_n.iter().zip(&released.paths_m) {
            assert_eq!(*force, [0.0; 3]);
            for c in 0..3 {
                close(
                    path[c],
                    [-a[2] * H * H, -a[0] * H * H, -a[1] * H * H][c],
                    2e-10,
                );
            }
        }
        for (element, tet) in motion.solid().geometry().tetrahedra().iter().enumerate() {
            let sign = if element == 0 { 1.0 } else { -1.0 };
            for &node in tet {
                for c in 0..3 {
                    close(
                        released.motion.velocity_mps()[node][c],
                        VELOCITY[c] + sign * a[c] * H,
                        2e-8,
                    );
                    close(
                        released.motion.solid().solution().displacement_m[node][c]
                            - loaded.motion.solid().solution().displacement_m[node][c],
                        VELOCITY[c] * H + 0.5 * sign * a[c] * H * H,
                        2e-10,
                    );
                }
            }
        }
        let work = dot(a, a) * H * H;
        close(released.energy.dead_load_mean_work_j, work, 2e-10);
        close(released.energy.support_mean_work_j, 0.0, 2e-10);
        close(released.energy.inertia.kinetic_change_j, work, 2e-10);
        audit_interval(&loaded.motion, &released, H, &separating, &[]);
        assert_eq!(released.motion.accepted_steps(), 2);
        assert_eq!(released.motion.remaining_steps(), 1);
    });
}

#[test]
fn rotated_partial_support_couples_static_tractions_and_refuses_sliding_transactionally() {
    let q = core::f64::consts::FRAC_1_SQRT_2;
    let n = [q, 0.0, q];
    let tangent = [q, 0.0, -q];
    let nodes = NODES.map(|x| [q * (x[0] + x[2]), x[1], q * (x[2] - x[0])]);
    with_oblique_motion(&nodes, n, [VELOCITY; 2], |motion, cx| {
        let a = core::array::from_fn(|c| n[c] + 0.25 * tangent[c]);
        let forces = body_forces(motion, a);
        let prescribed = supports(motion, H, [false, false, true]);
        let friction = [law([0.5, 0.25], tangent, [1.0, 0.0])];
        let original = motion.solid().solution().clone();
        let velocity = motion.velocity_mps().to_vec();

        // Required ft/fn=1/4 exceeds this first principal capacity. Native
        // Coulomb response has finite slip; paired-node geometry cannot follow
        // that trajectory and must not project it back into a fake stick.
        let weak = [law([0.125, 0.25], tangent, [1.0, 0.0])];
        assert!(matches!(
            motion.advance_released_face_sustained_frictional(
                H,
                &weak,
                &prescribed,
                &forces,
                settings(),
                sustained_controls(),
                cx,
            ),
            Err(BrittleFaceSustainedFrictionError::Sustained(
                BrittleFaceSustainedError::Mode { .. }
            ))
        ));
        let mut limited = sustained_controls();
        limited.contact.solve.maximum_response_evaluations = 1;
        let error = motion
            .advance_released_face_sustained_frictional(
                H,
                &friction,
                &prescribed,
                &forces,
                settings(),
                limited,
                cx,
            )
            .unwrap_err();
        let BrittleFaceSustainedFrictionError::Sustained(BrittleFaceSustainedError::Solve(error)) =
            error
        else {
            panic!("expected native force-response budget refusal");
        };
        assert!(matches!(
            *error,
            SolveError::ResponseBudget { evaluations: 1, .. }
        ));
        with_cx(true, |cancelled| {
            assert!(
                motion
                    .advance_released_face_sustained_frictional(
                        H,
                        &friction,
                        &prescribed,
                        &forces,
                        settings(),
                        sustained_controls(),
                        cancelled,
                    )
                    .is_err()
            );
        });
        assert_eq!(
            motion.solid().solution().displacement_m,
            original.displacement_m
        );
        assert_eq!(
            motion.solid().solution().element_states,
            original.element_states
        );
        assert_eq!(motion.velocity_mps(), velocity);
        assert_eq!(motion.accepted_steps(), 1);
        assert_eq!(motion.remaining_steps(), 1);

        let step = motion
            .advance_released_face_sustained_frictional(
                H,
                &friction,
                &prescribed,
                &forces,
                settings(),
                sustained_controls(),
                cx,
            )
            .unwrap();
        // Holding A's global Z gives physical uniform mobility diag(16,16,8).
        // In [n,t1,t2] it is [[12,4,0],[4,12,0],[0,0,16]]. Sticking solves
        // both rows: [fn,ft]=[1/8,1/32]. A normal-first solve gives 13/96.
        for force in &step.forces_n {
            for c in 0..3 {
                close(force[c], [1.0 / 8.0, 1.0 / 32.0, 0.0][c], 2e-8);
            }
        }
        let mut reactions = Vec::new();
        for (element, tet) in motion.solid().geometry().tetrahedra().iter().enumerate() {
            let sign = if element == 0 { 1.0 } else { -1.0 };
            for &node in tet {
                let is_face = step.rows.iter().any(|row| row.nodes[element] == node);
                if element == 0 {
                    reactions.push((
                        3 * node + 2,
                        if is_face { -a[2] / 8.0 } else { -a[2] / 4.0 },
                    ));
                }
                for c in 0..3 {
                    let acceleration = if is_face || (element == 0 && c == 2) {
                        0.0
                    } else {
                        sign * 2.5 * a[c]
                    };
                    close(
                        step.motion.velocity_mps()[node][c],
                        VELOCITY[c] + H * acceleration,
                        2e-8,
                    );
                    close(
                        step.motion.solid().solution().displacement_m[node][c]
                            - motion.solid().solution().displacement_m[node][c],
                        VELOCITY[c] * H + 0.5 * H * H * acceleration,
                        2e-10,
                    );
                }
            }
        }
        let dead = 5.0 * H * H * (2.0 * a[0] * a[0] + 2.0 * a[1] * a[1] + a[2] * a[2]) / 16.0;
        let support = -5.0 * a[2] * VELOCITY[2] * H / 8.0;
        close(step.energy.dead_load_mean_work_j, dead, 2e-10);
        close(step.energy.support_mean_work_j, support, 2e-10);
        close(step.energy.inertia.kinetic_change_j, dead + support, 2e-10);
        audit_interval(motion, &step, H, &forces, &reactions);
        let retry = motion
            .advance_released_face_sustained_frictional(
                H,
                &friction,
                &prescribed,
                &forces,
                settings(),
                sustained_controls(),
                cx,
            )
            .unwrap();
        assert_eq!(retry.motion.velocity_mps(), step.motion.velocity_mps());
        assert_eq!(
            retry.motion.solid().solution().displacement_m,
            step.motion.solid().solution().displacement_m
        );
        assert_eq!(retry.forces_n, step.forces_n);
        assert_eq!(retry.energy, step.energy);
    });
}
