//! Complete curved material contact and declared rigid-pair motion [G1/G3/G4].

use super::*;
use fs_couple::inelastic_impact::surface::{
    FiniteSolidSphere, FiniteSolidSurfaceControls,
    quadratic::{
        QuadraticSolidContactSurface, QuadraticSolidSurfaceControls,
        path::QuadraticSolidPathControls, scene::*,
    },
};
use fs_material::hyper::{discrete_gradient::DiscreteGradientSettings, mixed::VolumetricTreatment};
use fs_query::triangle::quadratic::QuadraticQueryControls;
use fs_solid::hyper3::mixed::{
    MixedLawField, MixedUpSettings,
    conservative::{ConservativeFormulation, ConservativeTetDefinition},
    visco::{MixedViscoLawField, MixedViscoTetSolid, dynamics::MixedViscoTetMotion},
};
use fs_solid::linear3::tet10::promote_p1_to_p2;
use fs_solver::nonlinear::NewtonKrylovConfig;

fn geometry_controls() -> QuadraticSolidSurfaceControls {
    QuadraticSolidSurfaceControls {
        surface: FiniteSolidSurfaceControls {
            activation_tolerance_m: 1e-6,
            maximum_contacts: 8,
            ..Default::default()
        },
        query: QuadraticQueryControls {
            distance_tolerance: 1e-9,
            witness_tolerance: 1e-5,
            ..Default::default()
        },
        maximum_witness_span_m: 3e-5,
    }
}

fn path_controls() -> QuadraticSolidPathControls {
    QuadraticSolidPathControls {
        initial: geometry_controls(),
        ..Default::default()
    }
}

fn sphere(id: u64, rigid_body: usize) -> FiniteSolidSphere {
    FiniteSolidSphere {
        id,
        rigid_body,
        center_body_m: Vec3::ZERO,
        radius_m: 0.125,
        friction_coefficient: 0.0,
        normal_law: NormalVelocityLaw::NewtonRestitution { coefficient: 0.25 },
    }
}

fn body(position: Vec3, momentum: Vec3) -> RigidBodyState {
    RigidBodyState::new(
        Pose::new(position, UnitQuaternion::IDENTITY).unwrap(),
        momentum,
        Vec3::ZERO,
    )
    .unwrap()
}

#[test]
fn curved_native_trace_and_rigid_pair_share_one_supported_momentum_solve() {
    with_cx(|gate, cx| {
        let (nodes, elements) = promote_p1_to_p2(&NODES, &TETS, None, cx).unwrap();
        let elastic = Hyperelastic::new(
            HyperelasticModel::MooneyRivlin {
                c10: 6.0,
                c01: 0.0,
                kappa: 120.0,
            },
            4.0,
        )
        .unwrap();
        let material = ViscoHyperelastic::new(
            elastic.clone(),
            vec![FiniteStrainBranch {
                energy: BranchEnergy::NeoHookean { mu_pa: 36.0 },
                tau_s: 0.07,
            }],
        )
        .unwrap();
        let displacement: Vec<_> = nodes
            .iter()
            .map(|p| [0.2 * p[1] * p[2], 0.0, 0.0])
            .collect();
        let prescribed: Vec<_> = displacement
            .as_flattened()
            .iter()
            .copied()
            .enumerate()
            .collect();
        let mut accepted = MixedViscoTetSolid::try_new(
            &nodes,
            &elements,
            MixedViscoLawField::Uniform(&material),
            VolumetricTreatment::OriginalMaterial,
            TetAssemblyBudget::standard(),
            4,
            64,
            cx,
        )
        .unwrap();
        let settings = MixedUpSettings {
            load_steps: 1,
            maximum_newton_iterations: 20,
            newton: NewtonKrylovConfig {
                absolute_tolerance: 2e-11,
                relative_tolerance: 0.0,
                ..Default::default()
            },
            length_scale_m: 1.0,
            stress_scale_pa: 12.0,
            maximum_assembly_bytes: 2 * 1024 * 1024,
            maximum_factor_bytes: 2 * 1024 * 1024,
        };
        for dt in [0.0, 0.007] {
            accepted = accepted
                .advance(
                    dt,
                    &prescribed,
                    &[],
                    [0.0; 3],
                    None,
                    settings,
                    FiniteViscoStepControls {
                        equilibrium_work_tolerance_j: 2e-9,
                        maximum_endpoint_work_defect_j: 1.0,
                    },
                    cx,
                )
                .unwrap()
                .solid;
        }
        let velocity = vec![[0.0; 3]; nodes.len()];
        let density = [6.0];
        let motion =
            MixedViscoTetMotion::try_new(accepted, &density, &velocity, 128 * 1024, cx).unwrap();
        let original = motion.solid().quadrature_states().to_vec();
        assert!(original.iter().any(|state| state.dissipated > 0.0));
        let pressure = motion.solid().pressure_pa().to_vec();
        let surface = QuadraticSolidContactSurface::for_mixed_visco(&motion, 4, cx).unwrap();
        let point = Vec3::new(0.0125, 0.25, 0.25);
        let n = Vec3::new(-1.0, 0.05, 0.05);
        let normal = n.scale(1.0 / norm(n));
        let gap = 2e-7;
        let masses = [2.0, 3.0];
        let properties = masses.map(|mass| {
            MassProperties::new(
                mass,
                Vec3::ZERO,
                Vec3::new(mass / 160.0, mass / 160.0, mass / 160.0),
            )
            .unwrap()
        });
        let first_center = point.add(normal.scale(0.125 + gap));
        let bodies = [
            body(first_center, normal.scale(-2.0)),
            body(
                first_center.add(normal.scale(0.25 + gap)),
                normal.scale(-6.0),
            ),
        ];
        let spheres = [sphere(11, 0), sphere(22, 1)];
        let law = SpherePairLaw {
            first: 11,
            second: 22,
            friction_coefficient: 0.0,
            restitution: 0.5,
        };
        let scene = FiniteSolidSphereScene::new(&spheres, &[law], 1).unwrap();
        let query = scene
            .mixed_visco_contacts(&surface, &motion, &bodies, geometry_controls(), cx)
            .unwrap();
        assert_eq!(query.contacts.len(), 2);
        assert_eq!(query.surface.contacts.len(), 1);
        assert_eq!(query.sphere_pair_candidates, 1);
        assert_eq!(
            query.keys[0],
            FiniteSolidSphereSceneKey::Spheres {
                first: 11,
                second: 22
            }
        );
        assert_eq!(query.contacts[0].id, 1);
        let witness = query.surface.contacts[0];
        assert_eq!(query.contacts[1].id, 2 * witness.contact.id);
        assert_eq!(query.contacts[1].counterpart, witness.contact.counterpart);
        assert_eq!(query.contacts[1].normal_world, witness.contact.normal_world);
        assert_eq!(query.contacts[1].point, witness.contact.point);
        assert!(witness.solid_point_world_m.x > 0.01);
        near(norm(witness.solid_point_world_m.sub(point)), 0.0, 3e-5);
        assert!(query.maximum_contact_gap_m <= geometry_controls().surface.activation_tolerance_m);

        let reordered = FiniteSolidSphereScene::new(&[spheres[1], spheres[0]], &[law], 1)
            .unwrap()
            .mixed_visco_contacts(&surface, &motion, &bodies, geometry_controls(), cx)
            .unwrap();
        assert_eq!(query.keys, reordered.keys);
        assert_eq!(query.contacts, reordered.contacts);

        // The same complete scene consumes an original exact-volume native
        // geometry without rebuilding its original material point coordinates.
        let conservative = ConservativeTetDefinition {
            nodes_m: &nodes,
            elements: &elements,
            material: MixedLawField::Uniform(&elastic),
            reference_density_kg_m3: &density,
            formulation: ConservativeFormulation::ExactIncompressible {
                maximum_relative_volume_residual: 1e-9,
                maximum_initial_relative_volume_rate_per_s: 1e-10,
                maximum_point_volume_work_defect: 1e-12,
            },
            budget: TetAssemblyBudget::standard(),
            quadrature_points_per_dir: 4,
            maximum_quadrature_points: 64,
            maximum_mass_assembly_bytes: 128 * 1024,
            material_settings: DiscreteGradientSettings {
                minimum_relative_strain_increment: 1e-6,
                absolute_work_tolerance_pa: 1e-10,
                relative_work_tolerance: 1e-12,
            },
        }
        .prepare_with_pressure(&displacement, &velocity, &[0.3; 4], cx)
        .unwrap();
        let conservative_query = scene
            .conservative_contacts(&surface, &conservative, &bodies, geometry_controls(), cx)
            .unwrap();
        assert_eq!(conservative_query.contacts, query.contacts);
        assert_eq!(conservative_query.keys, query.keys);

        let mut impact_controls = mixed_controls();
        impact_controls.geometry_tolerance_m = 2e-6;
        let held: Vec<_> = (0..3 * nodes.len()).collect();
        let event = mixed_visco_mixed_contact_candidate(
            &motion,
            &held,
            &bodies,
            &properties,
            &query.contacts,
            ORIGIN,
            impact_controls,
            cx,
        )
        .unwrap();
        // In (pair,solid) order, W=[[5/6,-1/2],[-1/2,1/2]]. Incoming
        // speeds are -1 and -1; independent targets are 1/2 and 1/4.
        // Solving gives J_pair=33/4, J_solid=43/4, and loss=195/32 J.
        near(event.report.contact_mobility_per_kg[0], 5.0 / 6.0, 2e-9);
        near(event.report.contact_mobility_per_kg[3], -0.5, 2e-9);
        near(event.report.contact_mobility_per_kg[21], 0.5, 2e-9);
        near(event.report.contact.impulses[0][0], 33.0 / 4.0, 1e-5);
        near(event.report.contact.impulses[1][0], 43.0 / 4.0, 1e-5);
        near(event.report.contact.dissipated_energy, 195.0 / 32.0, 1e-5);
        near(event.report.kinetic_change_j, -195.0 / 32.0, 1e-5);
        near(
            norm(
                event.rigid_bodies[0]
                    .linear_momentum_world()
                    .sub(normal.scale(0.5)),
            ),
            0.0,
            1e-4,
        );
        near(
            norm(
                event.rigid_bodies[1]
                    .linear_momentum_world()
                    .sub(normal.scale(2.25)),
            ),
            0.0,
            1e-4,
        );
        let impulse = witness
            .contact
            .normal_world
            .scale(event.report.contact.impulses[1][0]);
        let expected_torque = witness
            .solid_point_world_m
            .sub(vector(ORIGIN))
            .cross(impulse);
        assert!(norm(expected_torque) > 0.1);
        near(
            norm(vector(event.report.solid_impulse.support.linear_kg_m_s).sub(impulse)),
            0.0,
            2e-9,
        );
        near(
            norm(vector(event.report.solid_impulse.support.angular_kg_m2_s).sub(expected_torque)),
            0.0,
            2e-9,
        );
        mixed_audit(&event.report);
        assert_eq!(event.solid.solid().quadrature_states(), original);
        assert_eq!(event.solid.solid().pressure_pa(), pressure);
        assert_eq!(event.solid.solid().displacement_m(), displacement);
        assert_eq!(event.solid.solid().time_s(), 0.007);

        assert!(FiniteSolidSphereScene::new(&spheres, &[], 1).is_err());
        let limited = QuadraticSolidSurfaceControls {
            surface: FiniteSolidSurfaceControls {
                maximum_contacts: 1,
                ..geometry_controls().surface
            },
            ..geometry_controls()
        };
        assert!(matches!(
            scene.mixed_visco_contacts(&surface, &motion, &bodies, limited, cx),
            Err(QuadraticSphereSceneError::Scene(
                FiniteSolidSphereSceneError::ContactBudget {
                    maximum_contacts: 1
                }
            ))
        ));
        let penetrating = [
            bodies[0],
            body(
                first_center.add(normal.scale(0.25 - 1e-5)),
                bodies[1].linear_momentum_world(),
            ),
        ];
        assert!(matches!(
            scene.mixed_visco_contacts(&surface, &motion, &penetrating, geometry_controls(), cx),
            Err(QuadraticSphereSceneError::Scene(FiniteSolidSphereSceneError::Penetration {
                first: 11, second: 22, gap_m
            })) if gap_m < 0.0
        ));
        assert_eq!(motion.solid().quadrature_states(), original);
        gate.request();
        assert!(matches!(
            scene.mixed_visco_contacts(&surface, &motion, &bodies, geometry_controls(), cx),
            Err(QuadraticSphereSceneError::Cancelled)
        ));
    });
}

#[test]
fn complete_paths_retain_earliest_canonical_pair_and_shared_work_limits() {
    with_cx(|gate, cx| {
        let (nodes, elements) = promote_p1_to_p2(&NODES, &TETS, None, cx).unwrap();
        let surface = QuadraticSolidContactSurface::from_reference(
            &nodes,
            &elements,
            TetAssemblyBudget::standard(),
            4,
            cx,
        )
        .unwrap();
        let ids = [11, 22, 3, 4];
        let spheres = core::array::from_fn::<_, 4, _>(|i| sphere(ids[i], i));
        let laws: Vec<_> = [[3, 4], [3, 11], [3, 22], [4, 11], [4, 22], [11, 22]]
            .map(|[first, second]| SpherePairLaw {
                first,
                second,
                friction_coefficient: 0.0,
                restitution: 0.5,
            })
            .into();
        let scene = FiniteSolidSphereScene::new(&spheres, &laws, 6).unwrap();
        let start = [(-3.0, -3.0), (-1.0, -3.0), (-3.0, -5.0), (-1.0, -5.0)]
            .map(|(x, y)| body(Vec3::new(x, y, -3.0), Vec3::ZERO));
        let end = [(-1.0, -3.0), (-3.0, -3.0), (-1.0, -5.0), (-3.0, -5.0)]
            .map(|(x, y)| body(Vec3::new(x, y, -3.0), Vec3::ZERO));
        let controls = path_controls();
        let path = scene
            .affine_quadratic_sphere_paths_at_positions(
                &surface, &nodes, &nodes, &start, &end, controls, cx,
            )
            .unwrap();
        let QuadraticSphereScenePathVerdict::PossibleContact {
            key,
            fraction_start,
            fraction_end,
        } = path.verdict
        else {
            panic!("separated endpoints conceal both pair crossings");
        };
        assert_eq!(
            key,
            QuadraticSphereScenePathKey::Spheres {
                first: 3,
                second: 4
            }
        );
        assert!(fraction_start <= 7.0 / 16.0 && fraction_end >= 7.0 / 16.0);
        assert!(fraction_end - fraction_start <= controls.sweep.time.fraction_tolerance);
        assert_eq!(path.sphere_pair_candidates, 2);
        assert_eq!(path.sphere_pairs_checked, 2);
        assert_eq!(path.minimum_pair_clearance_m, None);
        let reordered = FiniteSolidSphereScene::new(
            &[spheres[3], spheres[1], spheres[0], spheres[2]],
            &laws.iter().rev().copied().collect::<Vec<_>>(),
            6,
        )
        .unwrap()
        .affine_quadratic_sphere_paths_at_positions(
            &surface, &nodes, &nodes, &start, &end, controls, cx,
        )
        .unwrap();
        assert_eq!(path, reordered);
        let mut limited = controls;
        limited.initial.surface.maximum_query_visits = path.boundary.visited_nodes
            + path.boundary.tested_triangles
            + path.sphere_pair_comparisons
            + path.sphere_pairs_checked
            - 1;
        assert!(matches!(
            scene.affine_quadratic_sphere_paths_at_positions(
                &surface, &nodes, &nodes, &start, &end, limited, cx,
            ),
            Err(QuadraticSphereSceneError::PathBudget)
        ));
        let clear = scene
            .affine_quadratic_sphere_paths_at_positions(
                &surface, &nodes, &nodes, &start, &start, controls, cx,
            )
            .unwrap();
        assert_eq!(clear.sphere_pair_candidates, 0);
        assert_eq!(clear.sphere_pairs_checked, 0);
        assert_eq!(
            clear.minimum_pair_clearance_m,
            Some(controls.initial.surface.activation_tolerance_m)
        );
        assert!(matches!(clear.verdict,
            QuadraticSphereScenePathVerdict::Clear { minimum_clearance_m }
                if minimum_clearance_m == controls.initial.surface.activation_tolerance_m));
        gate.request();
        assert!(matches!(
            scene.affine_quadratic_sphere_paths_at_positions(
                &surface, &nodes, &nodes, &start, &end, controls, cx,
            ),
            Err(QuadraticSphereSceneError::Cancelled)
        ));
    });
}

#[test]
fn declared_pair_paths_admit_exact_recession_and_refuse_one_ulp_penetration() {
    with_cx(|_, cx| {
        let (nodes, elements) = promote_p1_to_p2(&NODES, &TETS, None, cx).unwrap();
        let surface = QuadraticSolidContactSurface::from_reference(
            &nodes,
            &elements,
            TetAssemblyBudget::standard(),
            4,
            cx,
        )
        .unwrap();
        let scene = FiniteSolidSphereScene::new(
            &[sphere(7, 0), sphere(9, 1)],
            &[SpherePairLaw {
                first: 7,
                second: 9,
                friction_coefficient: 0.0,
                restitution: 0.0,
            }],
            1,
        )
        .unwrap();
        let start = [
            body(Vec3::new(-3.0, -3.0, -3.0), Vec3::ZERO),
            body(Vec3::new(-2.75, -3.0, -3.0), Vec3::ZERO),
        ];
        let receding = [start[0], body(Vec3::new(-2.5, -3.0, -3.0), Vec3::ZERO)];
        let report = scene
            .affine_quadratic_sphere_paths_at_positions(
                &surface,
                &nodes,
                &nodes,
                &start,
                &receding,
                path_controls(),
                cx,
            )
            .unwrap();
        assert!(matches!(report.verdict,
            QuadraticSphereScenePathVerdict::Clear { minimum_clearance_m } if minimum_clearance_m == 0.0));
        assert_eq!(report.minimum_pair_clearance_m, Some(0.0));
        let penetrating = [
            start[0],
            body(Vec3::new((-2.75_f64).next_down(), -3.0, -3.0), Vec3::ZERO),
        ];
        let unresolved = scene
            .affine_quadratic_sphere_paths_at_positions(
                &surface,
                &nodes,
                &nodes,
                &start,
                &penetrating,
                path_controls(),
                cx,
            )
            .unwrap();
        assert!(matches!(
            unresolved.verdict,
            QuadraticSphereScenePathVerdict::PossibleContact {
                key: QuadraticSphereScenePathKey::Spheres {
                    first: 7,
                    second: 9
                },
                ..
            }
        ));
        assert!(matches!(
            scene.affine_quadratic_sphere_paths_at_positions(
                &surface, &nodes, &nodes, &penetrating, &receding, path_controls(), cx,
            ),
            Err(QuadraticSphereSceneError::Scene(FiniteSolidSphereSceneError::Penetration {
                first: 7, second: 9, gap_m
            })) if gap_m < 0.0
        ));
    });
}
