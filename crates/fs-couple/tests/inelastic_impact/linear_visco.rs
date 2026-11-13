use super::*;
use fs_couple::inelastic_impact::linear_visco::{
    LinearViscoImpactControls, LinearViscoImpactError, LinearViscoImpactReport,
    apply_linear_visco_rigid_impact,
};
use fs_material::visco::tensor::{Mandel6, PronyBranch, TensorProny};
use fs_solid::linear3::TetStaticSolveConfig;
use fs_solid::linear3::visco::{
    SolidMeshRef, TetViscoField, TetViscoMaterial, ViscoDynamicConfig, ViscoDynamicState,
    ViscoElements, ViscoStateBudget,
};

fn diagonal(value: f64) -> Mandel6 {
    core::array::from_fn(|i| core::array::from_fn(|j| if i == j { value } else { 0.0 }))
}

fn material() -> TetViscoMaterial {
    TetViscoMaterial {
        law: TensorProny::new_mandel(
            diagonal(10_000.0),
            vec![PronyBranch::new(diagonal(3_000.0), 0.02).unwrap()],
            0.02,
        )
        .unwrap()
        .with_kelvin_viscosity_mandel(diagonal(0.2))
        .unwrap(),
        // Unit tetrahedron: mass=1 kg and M=(I+11^T)/20 in each axis.
        density_kg_m3: 6.0,
        material_state_identity: [81; 32],
    }
}

fn system(law: &TetViscoMaterial, free: &[usize], cx: &Cx<'_>) -> ViscoElements {
    let held: Vec<_> = (0..12).filter(|d| !free.contains(d)).collect();
    ViscoElements::new_with_full_reference_mass(
        SolidMeshRef::Tet4 {
            nodes_m: &NODES,
            elements: &TETS,
        },
        TetViscoField::Uniform(law),
        &held,
        TetAssemblyBudget::standard(),
        ViscoStateBudget {
            max_internal_variables: 1000,
        },
        cx,
    )
    .unwrap()
}

fn config() -> ViscoDynamicConfig {
    ViscoDynamicConfig {
        maximum_steps: 20,
        maximum_free_dofs: 12,
        solve: TetStaticSolveConfig {
            maximum_iterations: 100,
            relative_residual_tolerance: 1e-13,
        },
        absolute_energy_tolerance_j: 1e-12,
        relative_energy_tolerance: 1e-10,
    }
}

fn impact_controls() -> LinearViscoImpactControls {
    LinearViscoImpactControls {
        geometry_tolerance_m: 1e-12,
        velocity_tolerance_m_s: 1e-10,
        linear_tolerance_kg_m_s: 1e-10,
        angular_tolerance_kg_m2_s: 1e-10,
        energy_tolerance_j: 1e-10,
    }
}

fn request(restitution: f64) -> InelasticRigidContact<usize> {
    InelasticRigidContact {
        solid_point: 1,
        rigid_arm_body_m: Vec3::new(0.0, 0.0, 0.2),
        normal_world: Vec3::new(2.0, 0.0, 0.0),
        activation_tolerance_m: 0.0,
        restitution,
    }
}

fn point(state: &ViscoDynamicState<'_, '_>, node: usize) -> [f64; 3] {
    core::array::from_fn(|a| {
        NODES[node][a] + state.material_history().displacement_m()[3 * node + a]
    })
}

fn check_report(report: &LinearViscoImpactReport) {
    near(norm(report.linear_momentum_residual_kg_m_s), 0.0, 2e-11);
    near(norm(report.angular_momentum_residual_kg_m2_s), 0.0, 2e-11);
    near(report.endpoint_velocity_residual_m_s, 0.0, 2e-11);
    near(report.work_energy_residual_j, 0.0, 2e-11);
    near(report.dissipation_residual_j, 0.0, 2e-11);
    near(report.absolute_linear_residual_kg_m_s, 0.0, 2e-11);
    near(report.absolute_angular_residual_kg_m2_s, 0.0, 2e-11);
    near(report.absolute_work_residual_j, 0.0, 2e-11);
}

// Independent complete-state momenta for the analytic unit tetrahedron. This
// does not consume the implementation's impulse or residual receipt.
fn momenta(
    solid: &ViscoDynamicState<'_, '_>,
    rigid: RigidBodyState,
    origin: [f64; 3],
) -> (Vec3, Vec3) {
    let origin = vector(origin);
    let mut p = rigid.linear_momentum_world();
    let mut l = rigid.pose().position_world().sub(origin).cross(p).add(
        rigid
            .pose()
            .orientation()
            .rotate_body_to_world(rigid.angular_momentum_body()),
    );
    let sums: [f64; 3] =
        core::array::from_fn(|a| (0..4).map(|n| solid.velocity_mps()[3 * n + a]).sum());
    for node in 0..4 {
        let nodal = Vec3::new(
            (solid.velocity_mps()[3 * node] + sums[0]) / 20.0,
            (solid.velocity_mps()[3 * node + 1] + sums[1]) / 20.0,
            (solid.velocity_mps()[3 * node + 2] + sums[2]) / 20.0,
        );
        p = p.add(nodal);
        l = l.add(vector(point(solid, node)).sub(origin).cross(nodal));
    }
    (p, l)
}

fn unchanged_history(before: &ViscoDynamicState<'_, '_>, after: &ViscoDynamicState<'_, '_>) {
    assert_eq!(after.time_s(), before.time_s());
    assert_eq!(after.steps(), before.steps());
    assert_eq!(after.eigenstrain_work_j(), before.eigenstrain_work_j());
    assert_eq!(after.parameter_work_j(), before.parameter_work_j());
    assert_eq!(
        after.material_history().energy(),
        before.material_history().energy()
    );
    assert_eq!(
        after.material_history().displacement_m(),
        before.material_history().displacement_m()
    );
}

#[test]
fn g3_translated_prony_contact_retains_displacement_gap_and_angular_impulse() {
    with_cx(|_, cx| {
        let law = material();
        let offset = 0.015625;
        for shift in [0.0, 2.0_f64.powi(48), -2.0_f64.powi(48)] {
            let nodes = NODES.map(|[x, y, z]| [x, y, z + shift]);
            let mesh = ViscoElements::new_with_full_reference_mass(
                SolidMeshRef::Tet4 {
                    nodes_m: &nodes,
                    elements: &TETS,
                },
                TetViscoField::Uniform(&law),
                &[],
                TetAssemblyBudget::standard(),
                ViscoStateBudget {
                    max_internal_variables: 1000,
                },
                cx,
            )
            .unwrap();
            let prepared = mesh.dynamic_operator(1e-5, config(), cx).unwrap();
            let mut history = mesh.state(cx).unwrap();
            let mut displacement = [0.0; 12];
            displacement[5] = offset;
            history
                .step_internal_force(&displacement, 0.01, cx)
                .unwrap();
            let initial = prepared.start(history, &[0.0; 12], cx).unwrap();
            let before_rigid = RigidBodyState::new(
                Pose::new(vector(nodes[1]), UnitQuaternion::IDENTITY).unwrap(),
                Vec3::new(-2.0, 0.0, 0.0),
                Vec3::ZERO,
            )
            .unwrap();
            let mut solid = initial.clone();
            let mut rigid = before_rigid;
            let report = apply_linear_visco_rigid_impact(
                &mut solid,
                &mut rigid,
                properties(),
                InelasticRigidContact {
                    solid_point: 1,
                    rigid_arm_body_m: Vec3::new(0.0, 0.0, offset),
                    normal_world: Vec3::new(1.0, 0.0, 0.0),
                    activation_tolerance_m: 0.0,
                    restitution: 0.0,
                },
                nodes[1],
                impact_controls(),
                cx,
            )
            .unwrap();
            // The local z offset remains the moment arm even when the world
            // point rounds onto the reference node. I^-1_yy = 25/12 kg^-1 m^-2.
            let impulse = 1.0 / (16.0 + 0.5 + offset * offset * 25.0 / 12.0);
            near(report.contact.normal_impulse, impulse, 2e-12);
            near(rigid.angular_momentum_body().y, offset * impulse, 2e-12);
            assert_eq!(report.gap_m, 0.0);
            unchanged_history(&initial, &solid);
            check_report(&report);

            // Independent consistent-mass momenta about the translated node.
            // Form physical arms directly from the mesh and displacement,
            // without first rounding their absolute world coordinates.
            let sums: [f64; 3] =
                core::array::from_fn(|a| (0..4).map(|n| solid.velocity_mps()[3 * n + a]).sum());
            let mut linear = rigid.linear_momentum_world();
            let mut angular = rigid.angular_momentum_body();
            for node in 0..4 {
                let p = vector(core::array::from_fn(|a| {
                    (solid.velocity_mps()[3 * node + a] + sums[a]) / 20.0
                }));
                let arm = vector(core::array::from_fn(|a| {
                    (nodes[node][a] - nodes[1][a]) + displacement[3 * node + a]
                }));
                linear = linear.add(p);
                angular = angular.add(arm.cross(p));
            }
            near(
                norm(linear.sub(before_rigid.linear_momentum_world())),
                0.0,
                2e-11,
            );
            near(norm(angular), 0.0, 2e-11);

            let mut solid = initial.clone();
            let mut rigid =
                RigidBodyState::new(before_rigid.pose(), Vec3::new(0.0, 0.0, -2.0), Vec3::ZERO)
                    .unwrap();
            let before_rigid = rigid;
            let report = apply_linear_visco_rigid_impact(
                &mut solid,
                &mut rigid,
                properties(),
                InelasticRigidContact {
                    solid_point: 1,
                    rigid_arm_body_m: Vec3::new(0.0, 0.0, 2.0 * offset),
                    normal_world: Vec3::new(0.0, 0.0, 1.0),
                    activation_tolerance_m: 0.0,
                    restitution: 0.0,
                },
                nodes[1],
                impact_controls(),
                cx,
            )
            .unwrap();
            assert_eq!(report.gap_m, offset);
            assert!(!report.contact.impacted);
            assert_eq!(report.contact.normal_impulse, 0.0);
            assert_eq!(rigid, before_rigid);
            assert_eq!(solid.velocity_mps(), initial.velocity_mps());
            unchanged_history(&initial, &solid);
            check_report(&report);
        }
    });
}

#[test]
fn native_prony_contact_uses_consistent_mass_and_full_inertia_without_resetting_history() {
    with_cx(|_, cx| {
        let law = material();
        let mesh = system(&law, &(0..12).collect::<Vec<_>>(), cx);
        let prepared = mesh.dynamic_operator(1e-5, config(), cx).unwrap();
        let mut history = mesh.state(cx).unwrap();
        let mut preload = [0.0; 12];
        preload[3] = 1e-4;
        preload[7] = -2e-5;
        history.step_internal_force(&preload, 0.01, cx).unwrap();
        let velocity: Vec<_> = (0..4).flat_map(|_| [0.2, -0.01, 0.03]).collect();
        let mut initial = prepared.start(history, &velocity, cx).unwrap();
        initial.advance(&[0.0; 12], cx).unwrap();
        assert!(initial.material_history().energy().dissipated_j > 0.0);
        assert_eq!(initial.steps(), 1);
        for restitution in [0.0, 0.35, 1.0] {
            let mut solid = initial.clone();
            // The native gap is formed exactly, (pose + arm) - (reference
            // node + displacement), so the rigid point is placed on the
            // node itself: the pose takes the reference x and the arm the
            // displacement x. point() rounds their sum, which the exact gap
            // sees open by that rounding (1.1e-16 m) at zero activation
            // tolerance. An arm along the normal moves neither the moment
            // arm nor the mobility.
            let at = point(&solid, 1);
            let displacement_x = solid.material_history().displacement_m()[3];
            let mut rigid = rigid_at([NODES[1][0], at[1], at[2]], solid.velocity_mps()[3] - 0.4);
            let before_rigid = rigid;
            let before_momenta = momenta(&solid, rigid, ORIGIN);
            let report = apply_linear_visco_rigid_impact(
                &mut solid,
                &mut rigid,
                properties(),
                InelasticRigidContact {
                    rigid_arm_body_m: Vec3::new(displacement_x, 0.0, 0.2),
                    ..request(restitution)
                },
                ORIGIN,
                impact_controls(),
                cx,
            )
            .unwrap();
            let mobility = 16.0 + 7.0 / 12.0;
            let j = (1.0 + restitution) * 0.4 / mobility;
            near(report.inverse_effective_mass_per_kg, mobility, 2e-11);
            near(report.contact.normal_impulse, j, 2e-12);
            near(report.relative_normal_speed_before_m_s, -0.4, 2e-12);
            near(
                report.relative_normal_speed_after_m_s,
                0.4 * restitution,
                2e-11,
            );
            near(
                report.contact.dissipated_energy,
                0.5 * (1.0 - restitution) * 0.4 * j,
                2e-12,
            );
            for node in 0..4 {
                near(
                    solid.velocity_mps()[3 * node] - initial.velocity_mps()[3 * node],
                    if node == 1 { -16.0 * j } else { 4.0 * j },
                    2e-11,
                );
                for a in 1..3 {
                    assert_eq!(
                        solid.velocity_mps()[3 * node + a],
                        initial.velocity_mps()[3 * node + a]
                    );
                }
            }
            let after_momenta = momenta(&solid, rigid, ORIGIN);
            near(norm(after_momenta.0.sub(before_momenta.0)), 0.0, 2e-11);
            near(norm(after_momenta.1.sub(before_momenta.1)), 0.0, 2e-11);
            assert_eq!(rigid.pose(), before_rigid.pose());
            assert_eq!(report.solid_impulse.support_work_j, 0.0);
            unchanged_history(&initial, &solid);
            check_report(&report);
            near(
                solid.external_work_j() - initial.external_work_j(),
                report.solid_impulse.external_work_j,
                2e-12,
            );

            // Equality of visible energies is insufficient: exercise the
            // retained internal Maxwell variables with the same next ramp.
            let mut old_memory = initial.material_history().clone();
            let mut new_memory = solid.material_history().clone();
            let a = old_memory.step_internal_force(&preload, 0.003, cx).unwrap();
            let b = new_memory.step_internal_force(&preload, 0.003, cx).unwrap();
            assert_eq!(a.full_force_n, b.full_force_n);
            assert_eq!(a.full_mean_force_n, b.full_mean_force_n);
            assert_eq!(a.energy, b.energy);
            let step = solid.advance(&[0.0; 12], cx).unwrap();
            assert_eq!(solid.steps(), initial.steps() + 1);
            near(step.time_s, initial.time_s() + 1e-5, 1e-16);
            assert!(step.constitutive.step_dissipated_j >= 0.0);
        }
    });
}

#[test]
fn moving_support_contact_keeps_cross_mass_reactions_and_signed_actuator_work() {
    with_cx(|_, cx| {
        let law = material();
        for fully_held in [false, true] {
            let free = if fully_held { vec![] } else { vec![3] };
            let mesh = system(&law, &free, cx);
            let prepared = mesh.dynamic_operator(1e-5, config(), cx).unwrap();
            let mut velocity = [0.0; 12];
            for (node, speed) in [0.2, 0.05, -0.1, 0.3].into_iter().enumerate() {
                velocity[3 * node] = speed;
            }
            let mut solid = prepared
                .start_with_support_motion(mesh.state(cx).unwrap(), &velocity, cx)
                .unwrap();
            let before = solid.clone();
            let mut rigid = rigid_at(point(&solid, 1), -0.35);
            let before_momenta = momenta(&solid, rigid, ORIGIN);
            let report = apply_linear_visco_rigid_impact(
                &mut solid,
                &mut rigid,
                properties(),
                request(0.5),
                ORIGIN,
                impact_controls(),
                cx,
            )
            .unwrap();
            let mobility = if fully_held { 0.0 } else { 10.0 };
            let j = 0.6 / (mobility + 7.0 / 12.0);
            near(
                report.inverse_effective_mass_per_kg,
                mobility + 7.0 / 12.0,
                2e-11,
            );
            near(report.contact.normal_impulse, j, 2e-12);
            let mut support_p = Vec3::ZERO;
            let mut support_l = Vec3::ZERO;
            for node in 0..4 {
                let reaction = if fully_held {
                    if node == 1 { j } else { 0.0 }
                } else if node == 1 {
                    0.0
                } else {
                    -0.5 * j
                };
                near(
                    report.solid_impulse.support_impulse_kg_m_s[3 * node],
                    reaction,
                    2e-12,
                );
                let p = Vec3::new(reaction, 0.0, 0.0);
                support_p = support_p.add(p);
                support_l =
                    support_l.add(vector(point(&before, node)).sub(vector(ORIGIN)).cross(p));
                for a in 0..3 {
                    let d = 3 * node + a;
                    if !free.contains(&d) {
                        assert_eq!(solid.velocity_mps()[d], velocity[d]);
                    }
                    if a > 0 {
                        assert_eq!(report.solid_impulse.support_impulse_kg_m_s[d], 0.0);
                    }
                }
            }
            // Original mass coupling gives -J/2 on each held x coordinate.
            // Their distinct prescribed rates do signed work, never heat.
            let support_work = if fully_held { 0.05 * j } else { -0.2 * j };
            near(report.solid_impulse.support_work_j, support_work, 2e-12);
            near(solid.support_work_j(), support_work, 2e-12);
            near(solid.velocity_mps()[3], 0.05 - mobility * j, 2e-11);
            let after_momenta = momenta(&solid, rigid, ORIGIN);
            near(
                norm(after_momenta.0.sub(before_momenta.0).sub(support_p)),
                0.0,
                2e-11,
            );
            near(
                norm(after_momenta.1.sub(before_momenta.1).sub(support_l)),
                0.0,
                2e-11,
            );
            unchanged_history(&before, &solid);
            check_report(&report);
            assert_eq!(solid.material_history().energy().dissipated_j, 0.0);
            if fully_held {
                assert_eq!(report.solid_impulse.iterations, 0);
                assert_eq!(report.solid_impulse.kinetic_change_j, 0.0);
            }
        }
    });
}

#[test]
fn late_lost_rigid_recoil_and_cancellation_publish_neither_owner() {
    with_cx(|gate, cx| {
        let law = material();
        let mesh = system(&law, &(0..12).collect::<Vec<_>>(), cx);
        let prepared = mesh.dynamic_operator(1e-5, config(), cx).unwrap();
        let mut history = mesh.state(cx).unwrap();
        let mut preload = [0.0; 12];
        preload[3] = 1e-4;
        history.step_internal_force(&preload, 0.01, cx).unwrap();
        let mut solid = prepared.start(history, &[0.0; 12], cx).unwrap();
        solid.advance(&[0.0; 12], cx).unwrap();
        let before = solid.clone();
        let massive = MassProperties::from_inertia_tensor(
            1e20,
            Vec3::ZERO,
            fs_ga::Mat3 {
                m: [1e20, 0.0, 0.0, 0.0, 1e20, 0.0, 0.0, 0.0, 1e20],
            },
        )
        .unwrap();
        let mut rigid = RigidBodyState::new(
            Pose::new(vector(point(&solid, 1)), UnitQuaternion::IDENTITY).unwrap(),
            Vec3::new(-1e20, 0.0, 0.0),
            Vec3::ZERO,
        )
        .unwrap();
        let before_rigid = rigid;
        let request = InelasticRigidContact {
            rigid_arm_body_m: Vec3::ZERO,
            ..request(0.5)
        };
        let error = apply_linear_visco_rigid_impact(
            &mut solid,
            &mut rigid,
            massive,
            request,
            ORIGIN,
            impact_controls(),
            cx,
        )
        .unwrap_err();
        // The small recoil cannot change the stored -1e20 momentum. The
        // solid's own staged impulse is valid, but the joint event must fail.
        let LinearViscoImpactError::Audit(report) = error else {
            panic!("expected a late joint recoil refusal, got {error:?}");
        };
        assert!(report.absolute_linear_residual_kg_m_s > 0.01);
        assert_eq!(report.rigid_impulse.linear_increment_world, Vec3::ZERO);
        assert_eq!(rigid, before_rigid);
        assert_eq!(solid.velocity_mps(), before.velocity_mps());
        assert_eq!(solid.kinetic_energy_j(), before.kinetic_energy_j());
        assert_eq!(solid.external_work_j(), before.external_work_j());
        assert_eq!(solid.support_work_j(), before.support_work_j());
        unchanged_history(&before, &solid);
        gate.request();
        assert!(matches!(
            apply_linear_visco_rigid_impact(
                &mut solid,
                &mut rigid,
                massive,
                request,
                ORIGIN,
                impact_controls(),
                cx,
            ),
            Err(LinearViscoImpactError::Cancelled)
        ));
        assert_eq!(rigid, before_rigid);
        assert_eq!(solid.velocity_mps(), before.velocity_mps());
        assert_eq!(solid.kinetic_energy_j(), before.kinetic_energy_j());
        assert_eq!(solid.external_work_j(), before.external_work_j());
        assert_eq!(solid.support_work_j(), before.support_work_j());
        unchanged_history(&before, &solid);
    });
}
