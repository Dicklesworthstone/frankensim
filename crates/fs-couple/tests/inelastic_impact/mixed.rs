//! Signed quadratic contact, independently integrated mass and full rigid inertia.

use super::*;
use fs_material::hyper::mixed::VolumetricTreatment;
use fs_solid::hyper3::mixed::MixedUpSettings;
use fs_solid::hyper3::mixed::visco::{
    MixedViscoLawField, MixedViscoTetError, MixedViscoTetSolid,
    dynamics::{MixedViscoTetMotion, Tet10ImpulsePoint},
};
use fs_solid::linear3::tet10::{TET10_EDGES, promote_p1_to_p2};
use fs_solver::nonlinear::NewtonKrylovConfig;

fn mixed_settings() -> MixedUpSettings {
    MixedUpSettings {
        load_steps: 1,
        maximum_newton_iterations: 30,
        newton: NewtonKrylovConfig {
            absolute_tolerance: 2e-11,
            relative_tolerance: 0.0,
            ..NewtonKrylovConfig::default()
        },
        length_scale_m: 1.0,
        stress_scale_pa: 12.0,
        maximum_assembly_bytes: 2 * 1024 * 1024,
        maximum_factor_bytes: 2 * 1024 * 1024,
    }
}

#[test]
fn curved_quadratic_contact_preserves_pressure_memory_and_native_momentum() {
    with_cx(|_, cx| {
        let (nodes, elements) = promote_p1_to_p2(&NODES, &TETS, None, cx).unwrap();
        let law = ViscoHyperelastic::new(
            Hyperelastic::new(
                HyperelasticModel::MooneyRivlin {
                    c10: 6.0,
                    c01: 0.0,
                    kappa: 120.0,
                },
                4.0,
            )
            .unwrap(),
            vec![FiniteStrainBranch {
                energy: BranchEnergy::NeoHookean { mu_pa: 36.0 },
                tau_s: 0.07,
            }],
        )
        .unwrap();
        let new_solid = |treatment| {
            MixedViscoTetSolid::try_new(
                &nodes,
                &elements,
                MixedViscoLawField::Uniform(&law),
                treatment,
                TetAssemblyBudget::standard(),
                4,
                64,
                cx,
            )
            .unwrap()
        };
        let mut solid = new_solid(VolumetricTreatment::OriginalMaterial);
        let fixed: Vec<_> = nodes
            .iter()
            .enumerate()
            .flat_map(|(i, x)| {
                [
                    0.002 * x[0] + 0.002 * x[1] * x[2],
                    -0.001 * x[1],
                    0.001 * x[2],
                ]
                .into_iter()
                .enumerate()
                .map(move |(c, u)| (3 * i + c, u))
            })
            .collect();
        for dt in [0.0, 0.007] {
            solid = solid
                .advance(
                    dt,
                    &fixed,
                    &[],
                    [0.0; 3],
                    None,
                    mixed_settings(),
                    FiniteViscoStepControls {
                        equilibrium_work_tolerance_j: 2e-9,
                        maximum_endpoint_work_defect_j: 1.0,
                    },
                    cx,
                )
                .unwrap()
                .solid;
        }
        assert!(solid.quadrature_states().iter().all(|s| s.dissipated > 0.0));
        assert!(solid.pressure_pa().iter().any(|p| p.abs() > 0.01));
        let velocity = vec![[0.2, 0.0, 0.0]; nodes.len()];
        let motion =
            MixedViscoTetMotion::try_new(solid, &[6.0], &velocity, 128 * 1024, cx).unwrap();
        let material_point = Tet10ImpulsePoint {
            tetrahedron: 0,
            barycentric: [0.0, 0.25, 0.25, 0.5],
        };
        let point = motion
            .point_kinematics(material_point, cx)
            .unwrap()
            .position_m;
        // The quadratic yz displacement vanishes at every corner. A four-node
        // trace would place this contact 0.00025 m away in the normal direction.
        for (x, expected) in point.into_iter().zip([0.25075, 0.24975, 0.5005]) {
            near(x, expected, 2e-13);
        }
        let edge = TET10_EDGES
            .iter()
            .position(|e| e.contains(&1) && e.contains(&2))
            .unwrap();
        let free_node = elements[0][4 + edge];
        let held: Vec<_> = (0..3 * nodes.len())
            .filter(|&d| d != 3 * free_node)
            .collect();
        let request = InelasticRigidContact {
            solid_point: material_point,
            rigid_arm_body_m: Vec3::new(0.0, 0.0, 0.2),
            normal_world: Vec3::new(2.0, 0.0, 0.0),
            activation_tolerance_m: 0.0,
            restitution: 0.5,
        };
        let rigid = rigid_at(point, -0.8);
        let event = mixed_visco_rigid_impact_candidate(
            &motion,
            &held,
            rigid,
            properties(),
            request,
            ORIGIN,
            controls(),
            cx,
        )
        .unwrap();

        // Independent affine reference integrals: total mass 1 kg, M_edge,edge
        // =8/105 kg, integral rho*N_edge=1/5 kg and point N_edge=1/4.
        // Supported solid mobility is 105/128; full rigid mobility is 7/12.
        let j = 576.0 / 539.0;
        near(
            event.report.inverse_effective_mass_per_kg,
            539.0 / 384.0,
            2e-11,
        );
        near(event.report.contact.normal_impulse, j, 2e-11);
        near(event.report.relative_normal_speed_before_m_s, -1.0, 2e-13);
        near(event.report.relative_normal_speed_after_m_s, 0.5, 2e-11);
        near(event.report.contact.dissipated_energy, j / 4.0, 2e-11);
        near(event.report.support_work_j, 11.0 * j / 160.0, 2e-11);
        near(event.report.kinetic_change_j, -29.0 * j / 160.0, 2e-11);
        for (node, v) in event.solid.velocity_mps().iter().enumerate() {
            near(
                v[0],
                0.2 - if node == free_node {
                    105.0 * j / 32.0
                } else {
                    0.0
                },
                2e-11,
            );
            near(v[1], 0.0, 0.0);
            near(v[2], 0.0, 0.0);
        }
        near(event.rigid.linear_momentum_world().x, -1.6 + j, 2e-11);
        near(event.rigid.angular_momentum_body().y, 0.2 * j, 2e-11);
        let expected_torque = vector(point)
            .sub(vector(ORIGIN))
            .cross(Vec3::new(-j, 0.0, 0.0));
        near(
            norm(vector(event.report.solid_impulse.applied.angular_kg_m2_s).sub(expected_torque)),
            0.0,
            2e-11,
        );
        audit(&event.report);
        assert_eq!(event.rigid.pose(), rigid.pose());
        assert_eq!(
            event.solid.solid().displacement_m(),
            motion.solid().displacement_m()
        );
        assert_eq!(
            event.solid.solid().pressure_pa(),
            motion.solid().pressure_pa()
        );
        assert_eq!(
            event.solid.solid().quadrature_states(),
            motion.solid().quadrature_states()
        );
        assert_eq!(
            event.solid.solid().first_piola_pa(),
            motion.solid().first_piola_pa()
        );
        assert_eq!(
            event.solid.solid().remainder_energy_j(),
            motion.solid().remainder_energy_j()
        );
        assert_eq!(
            event.solid.solid().bulk_energy_j(),
            motion.solid().bulk_energy_j()
        );
        assert_eq!(event.solid.solid().time_s(), 0.007);
        assert_eq!(event.solid.reference_density_kg_m3(), &[6.0]);
        assert_eq!(motion.velocity_mps(), velocity);

        check_quadratic_manifold(&motion, &held, cx);

        let exact = MixedViscoTetMotion::try_new(
            new_solid(VolumetricTreatment::ExactIncompressible),
            &[6.0],
            &velocity,
            128 * 1024,
            cx,
        )
        .unwrap();
        assert!(matches!(
            mixed_visco_rigid_impact_candidate(
                &exact,
                &held,
                rigid_at([0.25, 0.25, 0.5], -0.8),
                properties(),
                request,
                ORIGIN,
                controls(),
                cx,
            ),
            Err(InelasticImpactError::Solid(
                MixedViscoTetError::ExactIncompressibleImpulse
            ))
        ));
        assert_eq!(exact.velocity_mps(), velocity);
    });
}

fn check_quadratic_manifold(motion: &MixedViscoTetMotion<'_>, held: &[usize], cx: &Cx<'_>) {
    use fs_couple::inelastic_impact::multiple::elliptic::{
        EllipticFriction, mixed_visco_mixed_elliptic_contact_candidate,
    };

    let points =
        [[0.0, 0.25, 0.25, 0.5], [0.0, 0.5, 0.25, 0.25]].map(|barycentric| Tet10ImpulsePoint {
            tetrahedron: 0,
            barycentric,
        });
    let contacts: [_; 2] = core::array::from_fn(|body| InelasticMixedContact {
        id: 20 + body as u64,
        point: RigidImpactPoint {
            body,
            arm_body: Vec3::new(0.0, 0.0, 0.2),
        },
        counterpart: FiniteSolidContactCounterpart::Solid(points[body]),
        normal_world: Vec3::new(1.0, 0.0, 0.0),
        tangent_world: Vec3::new(0.0, 1.0, 0.0),
        activation_tolerance_m: 1e-12,
        friction_coefficient: 0.0,
        normal_law: NormalVelocityLaw::NewtonRestitution { coefficient: 0.5 },
    });
    let bodies: [_; 2] = core::array::from_fn(|body| {
        let point = motion
            .point_kinematics(points[body], cx)
            .unwrap()
            .position_m;
        let state = rigid_at(point, -0.8 - body as f64);
        RigidBodyState::new(
            state.pose(),
            state.linear_momentum_world().add(Vec3::new(
                0.0,
                if body == 0 { 0.8 } else { -0.5 },
                0.6,
            )),
            Vec3::ZERO,
        )
        .unwrap()
    });
    let properties = [properties(); 2];
    let mut controls = mixed_controls();
    controls.solid_impulse.maximum_workspace_components = 8192;
    let event = mixed_visco_mixed_contact_candidate(
        motion,
        held,
        &bodies,
        &properties,
        &contacts,
        ORIGIN,
        controls,
        cx,
    )
    .unwrap();

    // Only the original edge-12 x component is free. Its mass is 8/105 kg
    // and these traces are 1/4 and 1/2. The two normal rows are therefore
    // [[539/384, 105/64], [105/64, 371/96]], including each full rigid tensor.
    // Solving this coupled 2x2 system against [3/2, 3] gives the values below;
    // two separate point solves cannot give these impulses.
    let matrix = &event.report.contact_mobility_per_kg;
    for (index, expected) in [
        (0, 539.0 / 384.0),
        (3, 105.0 / 64.0),
        (18, 105.0 / 64.0),
        (21, 371.0 / 96.0),
    ] {
        near(matrix[index], expected, 3e-11);
    }
    for (body, expected) in [576.0 / 1799.0, 1152.0 / 1799.0].into_iter().enumerate() {
        near(event.report.contact.impulses[body][0], expected, 2e-10);
        near(event.report.contact.impulses[body][1], 0.0, 1e-12);
        near(event.report.contact.impulses[body][2], 0.0, 1e-12);
        near(
            event.report.contact.relative_velocities_after[body][0],
            0.5 * (1 + body) as f64,
            2e-10,
        );
    }
    near(
        event.report.contact.dissipated_energy,
        720.0 / 1799.0,
        2e-10,
    );
    near(event.report.support_work_j, -162.0 / 8995.0, 2e-10);
    near(event.report.kinetic_change_j, -3762.0 / 8995.0, 2e-10);
    mixed_audit(&event.report);

    let circular = contacts.map(|contact| InelasticMixedContact {
        normal_law: NormalVelocityLaw::NewtonRestitution { coefficient: 0.0 },
        ..contact.with_friction(0.2)
    });
    let circle = mixed_visco_mixed_contact_candidate(
        motion,
        held,
        &bodies,
        &properties,
        &circular,
        ORIGIN,
        controls,
        cx,
    )
    .unwrap();
    let equal = mixed_visco_mixed_elliptic_contact_candidate(
        motion,
        held,
        &bodies,
        &properties,
        &circular.map(|c| {
            c.with_friction(EllipticFriction {
                coefficients_mu: [0.2; 2],
                principal_axis: [3.0, 4.0],
            })
        }),
        ORIGIN,
        controls,
        cx,
    )
    .unwrap();
    assert_eq!(equal.rigid_bodies, circle.rigid_bodies);
    assert_eq!(equal.solid.velocity_mps(), circle.solid.velocity_mps());
    assert_eq!(equal.report.contact.transformed, circle.report.contact);

    let free_z = mixed_visco_mixed_elliptic_contact_candidate(
        motion,
        held,
        &bodies,
        &properties,
        &circular.map(|c| {
            c.with_friction(EllipticFriction {
                coefficients_mu: [0.2, 0.0],
                principal_axis: [1.0, 0.0],
            })
        }),
        ORIGIN,
        controls,
        cx,
    )
    .unwrap();
    for impulse in &free_z.report.contact.impulses {
        near(impulse[2], 0.0, 1e-13);
        assert!(impulse[1].abs() > 1e-4);
    }
    for next in [&event.solid, &circle.solid, &free_z.solid] {
        assert_eq!(next.solid().pressure_pa(), motion.solid().pressure_pa());
        assert_eq!(
            next.solid().displacement_m(),
            motion.solid().displacement_m()
        );
        assert_eq!(
            next.solid().quadrature_states(),
            motion.solid().quadrature_states()
        );
        assert_eq!(next.solid().time_s(), motion.solid().time_s());
    }
    assert_eq!(motion.velocity_mps(), [[0.2, 0.0, 0.0]; 10]);
}
