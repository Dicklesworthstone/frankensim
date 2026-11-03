//! Analytic graded-mass/full-inertia contact, with retained constitutive history.
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_couple::inelastic_impact::*;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_material::hyper::{Hyperelastic, HyperelasticModel};
use fs_material::plastic::finite::{FiniteJ2, FiniteJ2Limits, FiniteJ2Parameters};
use fs_material::visco_hyper::{BranchEnergy, FiniteStrainBranch, ViscoHyperelastic};
use fs_mbd::{MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3};
use fs_solid::hyper3::{
    HyperTetSettings,
    plastic::{
        FinitePlasticLawField, FinitePlasticStepControls, FinitePlasticTetSolid,
        dynamics::FinitePlasticTetMotion,
    },
    visco::{
        FiniteViscoStepControls, FiniteViscoTetSolid, ViscoelasticLawField,
        dynamics::{FiniteViscoTetMotion, SolidImpulseControls, TetImpulsePoint},
    },
};
use fs_solid::linear3::TetAssemblyBudget;
use fs_solver::CheckedCgConfig;

const NODES: [[f64; 3]; 4] = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
];
const TETS: [[usize; 4]; 1] = [[0, 1, 2, 3]];
const RHO: [f64; 4] = [3.0, 15.0, 3.0, 3.0];
const POINT: TetImpulsePoint = TetImpulsePoint {
    tetrahedron: 0,
    weights: [0.0, 1.0, 0.0, 0.0],
};
const ORIGIN: [f64; 3] = [0.3, -0.2, 0.4];
fn with_cx<R>(f: impl FnOnce(&CancelGate, &Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(
            &gate,
            &Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 3,
                    kernel_id: 17,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            ),
        )
    })
}
fn controls() -> InelasticImpactControls {
    InelasticImpactControls {
        solid_impulse: SolidImpulseControls {
            solve: CheckedCgConfig {
                tolerance: 1e-13,
                max_iterations: 64,
                max_corrections: 2,
            },
            maximum_workspace_components: 4096,
            linear_tolerance_kg_m_s: 1e-11,
            angular_tolerance_kg_m2_s: 1e-11,
            work_tolerance_j: 1e-11,
        },
        geometry_tolerance_m: 1e-12,
        velocity_tolerance_m_s: 1e-10,
        linear_tolerance_kg_m_s: 1e-10,
        angular_tolerance_kg_m2_s: 1e-10,
        energy_tolerance_j: 1e-10,
    }
}
fn settings() -> HyperTetSettings {
    HyperTetSettings {
        load_steps: 1,
        relative_tolerance: 1e-11,
        ..Default::default()
    }
}
fn near(a: f64, b: f64, tol: f64) {
    assert!((a - b).abs() <= tol, "{a:e} != {b:e}, tolerance {tol:e}");
}
fn norm(v: Vec3) -> f64 {
    v.x.hypot(v.y).hypot(v.z)
}
fn vector(a: [f64; 3]) -> Vec3 {
    Vec3::new(a[0], a[1], a[2])
}
fn properties() -> MassProperties {
    // Off-diagonal I_xy must affect the normal effective mass. With arm .2z
    // and direction x, rotational mobility is .04*(I^-1)_yy = 1/12 kg^-1.
    MassProperties::from_inertia_tensor(
        2.0,
        Vec3::ZERO,
        fs_ga::Mat3 {
            m: [0.5, 0.1, 0.0, 0.1, 0.5, 0.0, 0.0, 0.0, 0.4],
        },
    )
    .unwrap()
}
fn rigid_at(point: [f64; 3], speed: f64) -> RigidBodyState {
    RigidBodyState::new(
        Pose::new(
            vector(point).sub(Vec3::new(0.0, 0.0, 0.2)),
            UnitQuaternion::IDENTITY,
        )
        .unwrap(),
        Vec3::new(2.0 * speed, 0.0, 0.0),
        Vec3::ZERO,
    )
    .unwrap()
}
fn contact(e: f64) -> InelasticRigidContact {
    InelasticRigidContact {
        solid_point: POINT,
        rigid_arm_body_m: Vec3::new(0.0, 0.0, 0.2),
        normal_world: Vec3::new(2.0, 0.0, 0.0),
        activation_tolerance_m: 0.0,
        restitution: e,
    }
}
fn audit(r: &InelasticImpactReport) {
    near(norm(r.linear_momentum_residual_kg_m_s), 0.0, 2e-11);
    near(norm(r.angular_momentum_residual_kg_m2_s), 0.0, 2e-11);
    near(r.endpoint_velocity_residual_m_s, 0.0, 2e-11);
    near(r.work_energy_residual_j, 0.0, 2e-11);
    near(r.dissipation_residual_j, 0.0, 2e-11);
    near(r.absolute_linear_residual_kg_m_s, 0.0, 2e-11);
    near(r.absolute_angular_residual_kg_m2_s, 0.0, 2e-11);
    near(r.absolute_work_residual_j, 0.0, 2e-11);
}

#[test]
fn relaxed_viscoelastic_body_collides_through_full_rigid_inertia_and_moving_supports() {
    with_cx(|gate, cx| {
        let law = ViscoHyperelastic::new(
            Hyperelastic::new(
                HyperelasticModel::NeoHookean {
                    mu: 1200.0,
                    lambda: 2400.0,
                },
                4.0,
            )
            .unwrap(),
            vec![FiniteStrainBranch {
                energy: BranchEnergy::NeoHookean { mu_pa: 3600.0 },
                tau_s: 0.07,
            }],
        )
        .unwrap();
        let mut solid = FiniteViscoTetSolid::try_new(
            &NODES,
            &TETS,
            ViscoelasticLawField::Uniform(&law),
            TetAssemblyBudget::default(),
            cx,
        )
        .unwrap();
        let fixed: Vec<_> = (0..12)
            .map(|dof| (dof, if dof == 4 { 1e-4 } else { 0.0 }))
            .collect();
        for dt in [0.0, 0.007] {
            solid = solid
                .advance(
                    dt,
                    &fixed,
                    &[],
                    [0.0; 3],
                    settings(),
                    FiniteViscoStepControls {
                        equilibrium_work_tolerance_j: 1e-9,
                        maximum_endpoint_work_defect_j: 1.0,
                    },
                    cx,
                )
                .unwrap()
                .solid;
        }
        assert!(solid.element_states()[0].dissipated > 0.0);
        let motion = FiniteViscoTetMotion::try_new_with_nodal_density(
            solid,
            &RHO,
            &[[0.2, 0.0, 0.0]; 4],
            cx,
        )
        .unwrap();
        let before = motion.solid().element_states().to_vec();
        let displacement = motion.solid().displacement_m().to_vec();
        let stress = motion.solid().element_first_piola_pa().to_vec();
        let storage = motion.solid().stored_energy_j();
        let point = motion.point_kinematics(POINT, cx).unwrap().position_m;
        let rigid = rigid_at(point, -0.8);
        let held: Vec<_> = (0..12).filter(|&dof| dof != 3).collect();
        let event = visco_rigid_impact_candidate(
            &motion,
            &held,
            rigid,
            properties(),
            contact(0.5),
            ORIGIN,
            controls(),
            cx,
        )
        .unwrap();
        // Exact rho*N_1^2=.15 kg. Full rigid mobility=1/2+1/12.
        // Combined w=20/3+7/12=29/4, so J=(1+.5)/w=6/29.
        let j = 6.0 / 29.0;
        near(
            event.report.inverse_effective_mass_per_kg,
            29.0 / 4.0,
            2e-12,
        );
        near(event.report.contact.normal_impulse, j, 2e-12);
        near(event.report.relative_normal_speed_before_m_s, -1.0, 1e-14);
        near(event.report.relative_normal_speed_after_m_s, 0.5, 2e-11);
        near(event.report.contact.dissipated_energy, 3.0 / 58.0, 2e-12);
        near(event.report.support_work_j, -7.0 / 145.0, 2e-12);
        near(event.report.kinetic_change_j, -0.1, 2e-12);
        for node in 0..4 {
            near(
                event.solid.velocity_mps()[node][0],
                0.2 - if node == 1 { j / 0.15 } else { 0.0 },
                2e-12,
            );
        }
        near(event.rigid.linear_momentum_world().x, -1.6 + j, 2e-12);
        near(event.rigid.angular_momentum_body().y, 0.2 * j, 2e-12);
        near(
            properties()
                .angular_velocity_body(event.rigid.angular_momentum_body())
                .x,
            -1.0 / 58.0,
            2e-12,
        );
        let expected_solid_moment = vector(point)
            .sub(vector(ORIGIN))
            .cross(Vec3::new(-j, 0.0, 0.0));
        near(
            norm(
                vector(event.report.solid_impulse.applied.angular_kg_m2_s)
                    .sub(expected_solid_moment),
            ),
            0.0,
            2e-12,
        );
        audit(&event.report);
        assert_eq!(event.rigid.pose(), rigid.pose());
        assert_eq!(event.solid.solid().element_states(), before);
        assert_eq!(event.solid.solid().displacement_m(), displacement);
        assert_eq!(event.solid.solid().element_first_piola_pa(), stress);
        assert_eq!(event.solid.solid().stored_energy_j(), storage);
        assert_eq!(event.solid.solid().time_s(), 0.007);
        assert_eq!(event.solid.reference_mass_kg(), motion.reference_mass_kg());
        assert_eq!(
            event.solid.reference_nodal_density_kg_m3(),
            Some(RHO.as_slice())
        );
        let mut wrong_line = point;
        wrong_line[1] += 1e-5;
        assert!(matches!(
            visco_rigid_impact_candidate(
                &motion,
                &held,
                rigid_at(wrong_line, -0.8),
                properties(),
                contact(0.5),
                ORIGIN,
                controls(),
                cx
            ),
            Err(InelasticImpactError::Invalid(_))
        ));
        let mut penetrating = point;
        penetrating[0] -= 1e-5;
        assert!(matches!(
            visco_rigid_impact_candidate(
                &motion,
                &held,
                rigid_at(penetrating, -0.8),
                properties(),
                contact(0.5),
                ORIGIN,
                controls(),
                cx
            ),
            Err(InelasticImpactError::Contact(
                fs_contact::impact::NormalImpactError::InvalidGap
            ))
        ));
        gate.request();
        assert!(matches!(
            visco_rigid_impact_candidate(
                &motion,
                &held,
                rigid,
                properties(),
                contact(0.5),
                ORIGIN,
                controls(),
                cx
            ),
            Err(InelasticImpactError::Cancelled)
        ));
        assert_eq!(motion.solid().element_states(), before);
        assert_eq!(motion.solid().displacement_m(), displacement);
        assert_eq!(motion.velocity_mps(), [[0.2, 0.0, 0.0]; 4]);
    });
}

#[test]
fn freely_moving_plastic_body_uses_original_graded_inverse_and_keeps_its_return_history() {
    with_cx(|_, cx| {
        let (mu, kappa, yield_pa, hardening) = (8000.0, 16000.0, 200.0, 800.0);
        let law = FiniteJ2::new(
            FiniteJ2Parameters {
                shear_modulus_pa: mu,
                bulk_modulus_pa: kappa,
                yield_stress_pa: yield_pa,
                hardening_modulus_pa: hardening,
            },
            FiniteJ2Limits {
                minimum_jacobian: 0.25,
                maximum_elastic_log_strain: 0.5,
                maximum_plastic_log_strain: 2.0,
                maximum_accumulated_plastic_strain: 3.0,
            },
        )
        .unwrap();
        let mut solid = FinitePlasticTetSolid::try_new(
            &NODES,
            &TETS,
            FinitePlasticLawField::Uniform(&law),
            TetAssemblyBudget::default(),
            cx,
        )
        .unwrap();
        let plastic = (3.0 * mu * 0.2 - yield_pa) / (3.0 * mu + hardening);
        for strain in [0.2, plastic] {
            let stretches = [strain.exp(), (-0.5 * strain).exp(), (-0.5 * strain).exp()];
            let fixed: Vec<_> = (0..12)
                .map(|dof| (dof, (stretches[dof % 3] - 1.0) * NODES[dof / 3][dof % 3]))
                .collect();
            solid = solid
                .advance(
                    &fixed,
                    &[],
                    [0.0; 3],
                    settings(),
                    FinitePlasticStepControls {
                        equilibrium_work_tolerance_j: 1e-5,
                        maximum_endpoint_work_defect_j: 200.0,
                        yield_tolerance_pa: 1e-4,
                    },
                    cx,
                )
                .unwrap()
                .solid;
        }
        assert!(solid.element_states()[0].alpha > 0.18);
        let motion = FinitePlasticTetMotion::try_new_with_nodal_density(
            solid,
            &RHO,
            &[[0.2, 0.0, 0.0]; 4],
            cx,
        )
        .unwrap();
        let states = motion.solid().element_states().to_vec();
        let deformation = motion.solid().element_deformation().to_vec();
        let displacement = motion.solid().displacement_m().to_vec();
        let storage = motion.solid().stored_energy_j();
        let solid_point = TetImpulsePoint {
            weights: [0.0, 0.5, 0.5, 0.0],
            ..POINT
        };
        let request = |e| InelasticRigidContact {
            solid_point,
            ..contact(e)
        };
        let point = motion.point_kinematics(solid_point, cx).unwrap().position_m;
        let rigid = rigid_at(point, -0.8);
        for e in [0.0, 1.0] {
            let event = plastic_rigid_impact_candidate(
                &motion,
                &[],
                rigid,
                properties(),
                request(e),
                ORIGIN,
                controls(),
                cx,
            )
            .unwrap();
            // Independent exact M^-1 B^T for this edge material point is
            // [-304,260,548,-304]/71. Mean density would give [-4,6,6,-4].
            let response = [-304.0, 260.0, 548.0, -304.0];
            let w = 404.0 / 71.0 + 7.0 / 12.0;
            let j = (1.0 + e) / w;
            near(event.report.inverse_effective_mass_per_kg, w, 2e-11);
            near(event.report.contact.normal_impulse, j, 2e-12);
            for node in 0..4 {
                near(
                    event.solid.velocity_mps()[node][0],
                    0.2 - j * response[node] / 71.0,
                    2e-11,
                );
            }
            near(event.report.relative_normal_speed_after_m_s, e, 2e-11);
            near(event.report.kinetic_change_j, -0.5 * (1.0 - e) * j, 2e-11);
            near(event.report.support_work_j, 0.0, 0.0);
            assert!(
                event
                    .report
                    .solid_impulse
                    .support_impulses_kg_m_s
                    .is_empty()
            );
            audit(&event.report);
            assert_eq!(event.solid.solid().element_states(), states);
            assert_eq!(event.solid.solid().element_deformation(), deformation);
            assert_eq!(event.solid.solid().displacement_m(), displacement);
            assert_eq!(event.solid.solid().stored_energy_j(), storage);
            assert_eq!(
                event.solid.reference_nodal_density_kg_m3(),
                Some(RHO.as_slice())
            );
            assert_eq!(event.rigid.pose(), rigid.pose());
        }
        let separating = rigid_at(point, 0.3);
        let unchanged = plastic_rigid_impact_candidate(
            &motion,
            &[],
            separating,
            properties(),
            request(0.5),
            ORIGIN,
            controls(),
            cx,
        )
        .unwrap();
        assert!(!unchanged.report.contact.impacted);
        assert_eq!(unchanged.report.contact.normal_impulse, 0.0);
        assert_eq!(unchanged.solid.velocity_mps(), motion.velocity_mps());
        assert_eq!(unchanged.rigid, separating);
        let mut limited = controls();
        limited.solid_impulse.maximum_workspace_components = 0;
        assert!(matches!(
            plastic_rigid_impact_candidate(
                &motion,
                &[],
                rigid,
                properties(),
                request(0.5),
                ORIGIN,
                limited,
                cx
            ),
            Err(InelasticImpactError::Solid(_))
        ));
        assert_eq!(motion.solid().element_states(), states);
        assert_eq!(motion.solid().displacement_m(), displacement);
        assert_eq!(motion.velocity_mps(), [[0.2, 0.0, 0.0]; 4]);
    });
}
