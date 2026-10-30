//! G1/G3/G4: original Hertz patches share actual finite body recoil and work.

use core::f64::consts::PI;
use fs_exec::CancelGate;
use fs_mbd::tangential_friction::RigidContactBinding;
use fs_mbd::tangential_friction::hertz::multiple::*;
use fs_mbd::tangential_friction::hertz::*;
use fs_mbd::{MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3};
use fs_tribo::partial_slip::TangentFrame;
use fs_tribo::sliding_patch::contact::*;

const DT: f64 = 0.01;

fn close(a: f64, b: f64, tol: f64) {
    assert!(
        (a - b).abs() <= tol,
        "{a:.17e} != {b:.17e}; tolerance {tol:e}"
    );
}
fn close_vec(a: Vec3, b: Vec3, tol: f64) {
    close(a.x, b.x, tol);
    close(a.y, b.y, tol);
    close(a.z, b.z, tol);
}
fn array(v: Vec3) -> [f64; 3] {
    [v.x, v.y, v.z]
}
fn vector(v: [f64; 3]) -> Vec3 {
    Vec3::new(v[0], v[1], v[2])
}
fn controls() -> RigidHertzIslandControls {
    RigidHertzIslandControls {
        mechanical: RigidHertzControls {
            maximum_bodies: 3,
            geometry_tolerance_m: 2e-12,
            velocity_tolerance_m_s: 3e-9,
            angular_velocity_tolerance_rad_s: 3e-8,
            linear_momentum_tolerance_ns: 3e-11,
            angular_momentum_tolerance_nms: 3e-11,
            energy_tolerance_j: 2e-9,
            contact: HertzWrenchConfig {
                maximum_static_iterations: 96,
                maximum_iterations: 256,
                maximum_evaluations: 512,
                force_tolerance_n: 5e-8,
                moment_tolerance_nm: 2e-9,
                displacement_tolerance_m: 2e-12,
                energy_tolerance_j: 2e-10,
                relative_matrix_tolerance: 1e-13,
                maximum_force_n: 20.0,
                maximum_spin_moment_nm: 2.0,
            },
        },
        maximum_contacts: 4,
        maximum_matrix_terms: 144,
        maximum_sweeps: 96,
        maximum_evaluations: 20000,
    }
}
fn law(n: f64, a: f64, mu: f64) -> HertzWrenchContact {
    HertzWrenchContact::new(HertzWrenchParameters {
        normal_force_n: n,
        contact_radius_m: a,
        static_mu: mu + 0.1,
        kinetic_mu: mu,
    })
    .unwrap()
}
fn state(pose: Pose, properties: MassProperties, v: Vec3, w: Vec3) -> RigidBodyState {
    RigidBodyState::new(
        pose,
        v.scale(properties.mass()),
        properties
            .angular_momentum_body_checked(pose.orientation().rotate_world_to_body(w))
            .unwrap(),
    )
    .unwrap()
}

struct Fixture {
    bodies: Vec<RigidBodyState>,
    properties: Vec<MassProperties>,
    bindings: [RigidContactBinding; 2],
    frames: [TangentFrame; 2],
    laws: [HertzWrenchContact; 2],
}
impl Fixture {
    fn contacts(&self) -> [RigidHertzIslandContact<'_>; 2] {
        [0, 1].map(|i| RigidHertzIslandContact {
            patch_id: if i == 0 { "upper" } else { "lower" },
            contact: RigidHertzContact {
                binding: self.bindings[i],
                frame: self.frames[i],
                law: &self.laws[i],
            },
        })
    }
    fn prepared(&self) -> PreparedRigidHertzIsland {
        PreparedRigidHertzIsland::new(
            &self.bodies,
            &self.properties,
            &self.contacts(),
            DT,
            controls(),
            &CancelGate::new(),
        )
        .unwrap()
    }
    fn laws(&self) -> [RigidHertzIslandLaw<'_>; 2] {
        [
            RigidHertzIslandLaw {
                patch_id: "upper",
                law: &self.laws[0],
            },
            RigidHertzIslandLaw {
                patch_id: "lower",
                law: &self.laws[1],
            },
        ]
    }
    fn run(&self) -> RigidHertzIslandCandidate {
        self.prepared()
            .advance(&self.laws(), &CancelGate::new())
            .unwrap()
    }
}
fn fixture(slip: bool, spin: bool) -> Fixture {
    let masses = [2.0, 3.0, 5.0];
    let inertia = [
        Vec3::new(0.08, 0.10, 0.12),
        Vec3::new(0.15, 0.20, 0.25),
        Vec3::new(0.25, 0.30, 0.40),
    ];
    let properties: Vec<_> = (0..3)
        .map(|i| MassProperties::new(masses[i], Vec3::ZERO, inertia[i]).unwrap())
        .collect();
    let heights = [0.2, 0.0, -0.3];
    let speeds = [2.0, 0.3, -1.0];
    let spins = [5.0, 1.0, -2.0];
    let bodies = (0..3)
        .map(|i| {
            state(
                Pose::new(Vec3::new(0.0, 0.0, heights[i]), UnitQuaternion::IDENTITY).unwrap(),
                properties[i],
                Vec3::new(if slip { speeds[i] } else { 0.0 }, 0.0, 0.0),
                Vec3::new(0.0, 0.0, if spin { spins[i] } else { 0.0 }),
            )
        })
        .collect();
    Fixture {
        bodies,
        properties,
        bindings: [
            RigidContactBinding {
                body_a: 0,
                body_b: 1,
                arm_a_body: Vec3::new(0.0, 0.0, -0.1),
                arm_b_body: Vec3::new(0.0, 0.0, 0.1),
            },
            RigidContactBinding {
                body_a: 1,
                body_b: 2,
                arm_a_body: Vec3::new(0.0, 0.0, -0.1),
                arm_b_body: Vec3::new(0.0, 0.0, 0.2),
            },
        ],
        frames: [TangentFrame::new([0.0, 0.0, 1.0], [1.0, 0.0, 0.0]).unwrap(); 2],
        laws: [law(10.0, 0.03, 0.4), law(7.0, 0.02, 0.4)],
    }
}
fn translation_mobility() -> [[f64; 2]; 2] {
    [
        [
            1.0 / 2.0 + 1.0 / 3.0 + 0.01 / 0.10 + 0.01 / 0.20,
            -1.0 / 3.0 + 0.01 / 0.20,
        ],
        [
            -1.0 / 3.0 + 0.01 / 0.20,
            1.0 / 3.0 + 1.0 / 5.0 + 0.01 / 0.20 + 0.04 / 0.30,
        ],
    ]
}

#[test]
fn two_sliding_patches_share_middle_body_recoil_and_independent_exact_heat() {
    let f = fixture(true, false);
    let prepared = f.prepared();
    let candidate = prepared.advance(&f.laws(), &CancelGate::new()).unwrap();
    assert_eq!(candidate.body_indices, vec![0, 1, 2]);
    assert_eq!(
        prepared.patch_ids().collect::<Vec<_>>(),
        vec!["lower", "upper"]
    );
    let forces = [-4.0, -2.8];
    let w = translation_mobility();
    let before = [1.7, 1.3];
    let expected_p = [
        forces[0] * DT,
        (forces[1] - forces[0]) * DT,
        -forces[1] * DT,
    ];
    let expected_l = [
        -0.1 * forces[0] * DT,
        -0.1 * (forces[0] + forces[1]) * DT,
        -0.2 * forces[1] * DT,
    ];
    let mut expected_heat = 0.0;
    for i in 0..3 {
        close_vec(
            candidate.bodies_after[i].linear_momentum_world(),
            f.bodies[i]
                .linear_momentum_world()
                .add(Vec3::new(expected_p[i], 0.0, 0.0)),
            2e-14,
        );
        close_vec(
            candidate.bodies_after[i].angular_momentum_body(),
            f.bodies[i]
                .angular_momentum_body()
                .add(Vec3::new(0.0, expected_l[i], 0.0)),
            2e-14,
        );
    }
    for i in 0..2 {
        let slot = 1 - i;
        let r = &candidate.report.contacts[slot];
        let after = before[i] + DT * (w[i][0] * forces[0] + w[i][1] * forces[1]);
        let heat = -forces[i] * DT * 0.5 * (before[i] + after);
        expected_heat += heat;
        close(r.contact.force_n[0], forces[i], 2e-12);
        close(r.slip_velocity_after_m_s[0], after, 2e-14);
        close(r.dissipated_energy_j, heat, 2e-13);
        close(r.contact.spin_moment_nm, 0.0, 2e-14);
    }
    // The negative off-block term is the actual finite middle-body response.
    close(prepared.scaled_mobility_per_kg()[3], w[0][1], 2e-15);
    close(candidate.report.dissipated_energy_j, expected_heat, 3e-13);
    close(
        candidate.report.kinetic_before_j - candidate.report.kinetic_after_j,
        expected_heat,
        3e-13,
    );
    close(candidate.report.total_energy_residual_j, 0.0, 3e-13);
}

#[test]
fn two_spinning_patches_retain_opposite_shared_spin_coupling() {
    let f = fixture(false, true);
    let prepared = f.prepared();
    let candidate = f.run();
    let moments = [-3.0 * PI / 16.0 * 4.0 * 0.03, -3.0 * PI / 16.0 * 2.8 * 0.02];
    let w = [
        [1.0 / 0.12 + 1.0 / 0.25, -1.0 / 0.25],
        [-1.0 / 0.25, 1.0 / 0.25 + 1.0 / 0.40],
    ];
    let expected_l = [
        moments[0] * DT,
        (moments[1] - moments[0]) * DT,
        -moments[1] * DT,
    ];
    for i in 0..3 {
        close_vec(
            candidate.bodies_after[i].linear_momentum_world(),
            f.bodies[i].linear_momentum_world(),
            3e-14,
        );
        close_vec(
            candidate.bodies_after[i].angular_momentum_body(),
            f.bodies[i]
                .angular_momentum_body()
                .add(Vec3::new(0.0, 0.0, expected_l[i])),
            3e-14,
        );
    }
    for i in 0..2 {
        let before = [4.0, 3.0][i];
        let after = before + DT * (w[i][0] * moments[0] + w[i][1] * moments[1]);
        let r = &candidate.report.contacts[1 - i];
        close(r.contact.spin_moment_nm, moments[i], 3e-13);
        close(r.spin_after_rad_s, after, 3e-13);
        close(
            r.dissipated_energy_j,
            -moments[i] * DT * 0.5 * (before + after),
            3e-13,
        );
    }
    close(
        prepared.scaled_mobility_per_kg()[2 * 6 + 5],
        -0.02 * 0.03 / 0.25,
        3e-17,
    );
}

fn inverse(m: [f64; 9]) -> [f64; 9] {
    let [a, b, c, d, e, f, g, h, i] = m;
    let det = a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g);
    [
        (e * i - f * h) / det,
        (c * h - b * i) / det,
        (b * f - c * e) / det,
        (f * g - d * i) / det,
        (a * i - c * g) / det,
        (c * d - a * f) / det,
        (d * h - e * g) / det,
        (b * g - a * h) / det,
        (a * e - b * d) / det,
    ]
}
fn matvec(m: [f64; 9], v: Vec3) -> Vec3 {
    Vec3::new(
        m[0] * v.x + m[1] * v.y + m[2] * v.z,
        m[3] * v.x + m[4] * v.y + m[5] * v.z,
        m[6] * v.x + m[7] * v.y + m[8] * v.z,
    )
}
fn dense_mobility(f: &Fixture) -> Vec<f64> {
    // Independent J M^-1 J^T with an explicit adjugate inverse, not the
    // production principal-square-root factors or its two-body sparse Gram.
    let mut translation = vec![[Vec3::ZERO; 3]; 6];
    let mut rotation = vec![[Vec3::ZERO; 3]; 6];
    for slot in 0..2 {
        let i = 1 - slot;
        let b = f.bindings[i];
        let axes = [
            vector(f.frames[i].longitudinal()),
            vector(f.frames[i].lateral()),
            vector(f.frames[i].normal()),
        ];
        for (body, arm, sign) in [
            (b.body_a, b.arm_a_body, 1.0),
            (b.body_b, b.arm_b_body, -1.0),
        ] {
            for axis in 0..3 {
                let body_axis = f.bodies[body]
                    .pose()
                    .orientation()
                    .rotate_world_to_body(axes[axis]);
                if axis < 2 {
                    translation[3 * slot + axis][body] = axes[axis].scale(sign);
                    rotation[3 * slot + axis][body] = arm.cross(body_axis).scale(sign);
                } else {
                    rotation[3 * slot + axis][body] =
                        body_axis.scale(sign * f.laws[i].parameters().contact_radius_m);
                }
            }
        }
    }
    let mut matrix = vec![0.0; 36];
    for row in 0..6 {
        for col in 0..6 {
            for b in 0..3 {
                matrix[6 * row + col] += translation[row][b].dot(translation[col][b])
                    / f.properties[b].mass()
                    + rotation[row][b].dot(matvec(
                        inverse(f.properties[b].inertia_tensor_body().m),
                        rotation[col][b],
                    ));
            }
        }
    }
    matrix
}

#[test]
fn full_tensor_cross_patch_matrix_world_motion_and_tangent_gauges_are_covariant() {
    let mut f = fixture(true, true);
    for i in 0..3 {
        let old = f.properties[i];
        let basis = UnitQuaternion::from_axis_angle(Vec3::new(1.0, 2.0, 3.0), 0.4 + 0.1 * i as f64)
            .unwrap();
        f.properties[i] = MassProperties::from_principal_frame(
            old.mass(),
            Vec3::ZERO,
            old.principal_inertia_body(),
            basis,
        )
        .unwrap();
        f.bodies[i] = state(
            f.bodies[i].pose(),
            f.properties[i],
            Vec3::new([2.0, 0.3, -1.0][i], 0.0, 0.0),
            Vec3::new(0.0, 0.0, [5.0, 1.0, -2.0][i]),
        );
    }
    let prepared = f.prepared();
    for (&actual, reference) in prepared
        .scaled_mobility_per_kg()
        .iter()
        .zip(dense_mobility(&f))
    {
        close(actual, reference, 3e-14);
    }
    assert!(prepared.scaled_mobility_per_kg()[2].abs() > 1e-5);
    let reference = f.run();
    let q = UnitQuaternion::from_axis_angle(Vec3::new(1.0, -2.0, 0.5), 0.63).unwrap();
    let shift = Vec3::new(0.7, -0.3, 0.2);
    let boost = Vec3::new(0.2, -0.1, 0.05);
    for i in 0..3 {
        let b = f.bodies[i];
        f.bodies[i] = RigidBodyState::new(
            Pose::new(
                q.rotate_body_to_world(b.pose().position_world()).add(shift),
                q,
            )
            .unwrap(),
            q.rotate_body_to_world(b.linear_momentum_world())
                .add(boost.scale(f.properties[i].mass())),
            b.angular_momentum_body(),
        )
        .unwrap();
    }
    for i in 0..2 {
        let angle = [0.31_f64, -0.72][i];
        f.frames[i] = TangentFrame::new(
            array(q.rotate_body_to_world(Vec3::new(0.0, 0.0, 1.0))),
            array(q.rotate_body_to_world(Vec3::new(angle.cos(), angle.sin(), 0.0))),
        )
        .unwrap();
    }
    let rotated = f.run();
    for i in 0..3 {
        close_vec(
            rotated.bodies_after[i].linear_momentum_world(),
            q.rotate_body_to_world(reference.bodies_after[i].linear_momentum_world())
                .add(boost.scale(f.properties[i].mass())),
            2e-8,
        );
        close_vec(
            rotated.bodies_after[i].angular_momentum_body(),
            reference.bodies_after[i].angular_momentum_body(),
            2e-9,
        );
    }
    close(
        rotated.report.dissipated_energy_j,
        reference.report.dissipated_energy_j,
        2e-8,
    );
    for i in 0..2 {
        close_vec(
            rotated.report.contacts[i].force_impulse_world_ns,
            q.rotate_body_to_world(reference.report.contacts[i].force_impulse_world_ns),
            2e-8,
        );
        close_vec(
            rotated.report.contacts[i].spin_angular_impulse_world_nms,
            q.rotate_body_to_world(reference.report.contacts[i].spin_angular_impulse_world_nms),
            2e-9,
        );
    }
}

#[test]
fn canonical_patch_order_and_energetic_spectators_preserve_exact_native_result() {
    let mut f = fixture(true, true);
    let reference = f.run();
    let p = MassProperties::new(2.0, Vec3::ZERO, Vec3::new(0.2, 0.3, 0.4)).unwrap();
    f.properties.push(p);
    f.bodies.push(state(
        Pose::new(Vec3::new(3.0, 2.0, -1.0), UnitQuaternion::IDENTITY).unwrap(),
        p,
        Vec3::new(1e12, -2e12, 3e12),
        Vec3::new(4e11, -5e11, 6e11),
    ));
    let mut contacts = f.contacts();
    contacts.reverse();
    let prepared = PreparedRigidHertzIsland::new(
        &f.bodies,
        &f.properties,
        &contacts,
        DT,
        controls(),
        &CancelGate::new(),
    )
    .unwrap();
    let mut laws = f.laws();
    laws.reverse();
    assert_eq!(
        prepared.advance(&laws, &CancelGate::new()).unwrap(),
        reference
    );
    let saved = f.bodies.clone();
    let mut published = saved.clone();
    let report = apply_rigid_hertz_island(
        &mut published,
        &f.properties,
        &contacts,
        DT,
        controls(),
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(report, reference.report);
    assert_eq!(&published[..3], reference.bodies_after);
    assert_eq!(published[3], saved[3]);
    // Changing the caller's later container cannot alter the prepared trial.
    f.bodies[0] = RigidBodyState::new(f.bodies[0].pose(), Vec3::ZERO, Vec3::ZERO).unwrap();
    assert_eq!(
        prepared.advance(&f.laws(), &CancelGate::new()).unwrap(),
        reference
    );
}

#[test]
fn prepared_coefficient_trials_and_single_patch_reduce_to_original_physics() {
    let f = fixture(true, false);
    let prepared = f.prepared();
    let w = translation_mobility();
    for mu in [0.0, 0.2, 0.4] {
        let trials = [law(10.0, 0.03, mu), law(7.0, 0.02, mu)];
        let keyed = [
            RigidHertzIslandLaw {
                patch_id: "upper",
                law: &trials[0],
            },
            RigidHertzIslandLaw {
                patch_id: "lower",
                law: &trials[1],
            },
        ];
        let candidate = prepared.advance(&keyed, &CancelGate::new()).unwrap();
        let forces = [-10.0 * mu, -7.0 * mu];
        for i in 0..2 {
            let before = [1.7, 1.3][i];
            let after = before + DT * (w[i][0] * forces[0] + w[i][1] * forces[1]);
            close(
                candidate.report.contacts[1 - i].dissipated_energy_j,
                -forces[i] * DT * 0.5 * (before + after),
                3e-12,
            );
        }
    }
    let original = f.contacts()[0];
    let pair = PreparedRigidHertzPair::new(
        &f.bodies,
        &f.properties,
        original.contact,
        DT,
        controls().mechanical,
        &CancelGate::new(),
    )
    .unwrap();
    let expected = pair.advance(&f.laws[0], &CancelGate::new()).unwrap();
    let one = PreparedRigidHertzIsland::new(
        &f.bodies,
        &f.properties,
        &[original],
        DT,
        controls(),
        &CancelGate::new(),
    )
    .unwrap()
    .advance(&[f.laws()[0]], &CancelGate::new())
    .unwrap();
    assert_eq!(one.bodies_after, expected.bodies_after);
    close(
        one.report.dissipated_energy_j,
        expected.report.dissipated_energy_j,
        2e-14,
    );
    let rest = fixture(false, false);
    let stationary = rest.run();
    assert_eq!(stationary.bodies_after, rest.bodies);
    assert_eq!(stationary.report.dissipated_energy_j, 0.0);
    assert!(
        stationary
            .report
            .contacts
            .iter()
            .all(|r| r.contact.mode == HertzWrenchMode::Sticking)
    );
}

#[test]
fn patch_identity_geometry_rank_effort_and_work_budgets_refuse_without_publication() {
    let f = fixture(true, false);
    let saved = f.bodies.clone();
    let prepared = f.prepared();
    let changed = law(11.0, 0.03, 0.4);
    assert!(
        prepared
            .advance(
                &[
                    RigidHertzIslandLaw {
                        patch_id: "upper",
                        law: &changed
                    },
                    f.laws()[1]
                ],
                &CancelGate::new()
            )
            .is_err()
    );
    assert!(
        prepared
            .advance(&[f.laws()[0], f.laws()[0]], &CancelGate::new())
            .is_err()
    );
    let mut duplicate = f.contacts();
    duplicate[1].patch_id = "upper";
    assert!(
        PreparedRigidHertzIsland::new(
            &saved,
            &f.properties,
            &duplicate,
            DT,
            controls(),
            &CancelGate::new()
        )
        .is_err()
    );
    // Distinct names do not make duplicate kinematic rows independent. This
    // admitted SPD rung refuses the rank defect without diagonal compliance.
    duplicate = [f.contacts()[0], f.contacts()[0]];
    duplicate[1].patch_id = "another";
    assert!(
        PreparedRigidHertzIsland::new(
            &saved,
            &f.properties,
            &duplicate,
            DT,
            controls(),
            &CancelGate::new()
        )
        .is_err()
    );
    for limited in [
        RigidHertzIslandControls {
            maximum_contacts: 1,
            ..controls()
        },
        RigidHertzIslandControls {
            maximum_matrix_terms: 35,
            ..controls()
        },
        RigidHertzIslandControls {
            maximum_evaluations: 1,
            ..controls()
        },
        RigidHertzIslandControls {
            mechanical: RigidHertzControls {
                contact: HertzWrenchConfig {
                    maximum_force_n: 0.01,
                    ..controls().mechanical.contact
                },
                ..controls().mechanical
            },
            ..controls()
        },
    ] {
        let mut bodies = saved.clone();
        assert!(
            apply_rigid_hertz_island(
                &mut bodies,
                &f.properties,
                &f.contacts(),
                DT,
                limited,
                &CancelGate::new()
            )
            .is_err()
        );
        assert_eq!(bodies, saved);
    }
    assert_eq!(
        prepared.advance(&f.laws(), &CancelGate::new()).unwrap(),
        f.run()
    );
}

#[test]
fn coarse_events_late_global_angular_audit_and_cancel_leave_exact_retry() {
    let f = fixture(true, false);
    let saved = f.bodies.clone();
    let mut bodies = saved.clone();
    assert!(
        apply_rigid_hertz_island(
            &mut bodies,
            &f.properties,
            &f.contacts(),
            2.0,
            controls(),
            &CancelGate::new()
        )
        .is_err()
    );
    assert_eq!(bodies, saved);
    // Geometry admission alone is insufficient: a tolerated endpoint mismatch
    // generates an unbalanced couple caught by the final whole-island audit.
    let old = bodies[2];
    bodies[2] = RigidBodyState::new(
        Pose::new(
            old.pose().position_world().add(Vec3::new(0.0, 0.0, 1e-6)),
            old.pose().orientation(),
        )
        .unwrap(),
        old.linear_momentum_world(),
        old.angular_momentum_body(),
    )
    .unwrap();
    let shifted = bodies.clone();
    let loose_geometry = RigidHertzIslandControls {
        mechanical: RigidHertzControls {
            geometry_tolerance_m: 2e-6,
            ..controls().mechanical
        },
        ..controls()
    };
    assert!(matches!(
        apply_rigid_hertz_island(
            &mut bodies,
            &f.properties,
            &f.contacts(),
            DT,
            loose_geometry,
            &CancelGate::new()
        ),
        Err(RigidHertzIslandError::Audit(_))
    ));
    assert_eq!(bodies, shifted);
    let gate = CancelGate::new();
    gate.request();
    bodies = saved.clone();
    assert!(
        apply_rigid_hertz_island(
            &mut bodies,
            &f.properties,
            &f.contacts(),
            DT,
            controls(),
            &gate
        )
        .is_err()
    );
    assert_eq!(bodies, saved);
    let mut retry = saved.clone();
    let result = apply_rigid_hertz_island(
        &mut retry,
        &f.properties,
        &f.contacts(),
        DT,
        controls(),
        &CancelGate::new(),
    )
    .unwrap();
    let clean = f.run();
    assert_eq!(retry, clean.bodies_after);
    assert_eq!(result, clean.report);
}
