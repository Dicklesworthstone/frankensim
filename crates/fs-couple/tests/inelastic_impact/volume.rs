//! G1/G3: simultaneous boundary impacts through the original P2/P1 pressure factor.

use super::*;
use fs_couple::inelastic_impact::multiple::elliptic::EllipticFriction;
use fs_couple::inelastic_impact::multiple::volume::*;
use fs_material::hyper::discrete_gradient::DiscreteGradientSettings;
use fs_material::hyper::mixed::VolumetricTreatment;
use fs_solid::hyper3::mixed::MixedLawField;
use fs_solid::hyper3::mixed::conservative::{ConservativeFormulation, ConservativeTetDefinition};
use fs_solid::hyper3::mixed::impulse::MixedImpulseControls;
use fs_solid::hyper3::mixed::visco::{
    MixedViscoLawField, MixedViscoTetSolid,
    dynamics::{MixedViscoTetMotion, Tet10ImpulsePoint},
};
use fs_solid::linear3::tet10::promote_p1_to_p2;

pub(super) fn volume_controls() -> VolumeContactControls {
    let base = mixed_controls();
    let mut balance = base.solid_impulse;
    balance.maximum_workspace_components = 8192;
    balance.solve.tolerance = 1e-11;
    VolumeContactControls {
        maximum_rigid_bodies: base.maximum_rigid_bodies,
        geometry_tolerance_m: base.geometry_tolerance_m,
        velocity_tolerance_m_s: base.velocity_tolerance_m_s,
        linear_tolerance_kg_m_s: base.linear_tolerance_kg_m_s,
        angular_tolerance_kg_m2_s: base.angular_tolerance_kg_m2_s,
        energy_tolerance_j: base.energy_tolerance_j,
        solid_impulse: MixedImpulseControls {
            balance,
            velocity_scale_m_s: 1.0,
            pressure_impulse_scale_pa_s: 1.0,
            work_scale_j: 1.0,
            maximum_assembly_bytes: 2 * 1024 * 1024,
            maximum_factor_bytes: 2 * 1024 * 1024,
            maximum_relative_volume_rate_per_s: 1e-10,
        },
        contact: base.contact,
    }
}

fn volume_audit<R>(report: &VolumeContactReport<R>) {
    near(norm(report.linear_momentum_residual_kg_m_s), 0.0, 2e-10);
    near(norm(report.angular_momentum_residual_kg_m2_s), 0.0, 2e-10);
    near(report.work_energy_residual_j, 0.0, 2e-10);
    near(report.dissipation_residual_j, 0.0, 2e-10);
    near(report.endpoint_velocity_residual_m_s, 0.0, 2e-10);
    near(report.absolute_linear_residual_kg_m_s, 0.0, 2e-10);
    near(report.absolute_angular_residual_kg_m2_s, 0.0, 2e-10);
    near(report.absolute_work_residual_j, 0.0, 2e-10);
    near(report.solid_impulse.pressure_work_j, 0.0, 1e-11);
    for rate in report
        .solid_impulse
        .volume_rate_before_m3_s
        .iter()
        .chain(&report.solid_impulse.volume_rate_after_m3_s)
    {
        near(*rate, 0.0, 1e-10 / 6.0);
    }
}

#[test]
fn exact_boundary_manifold_retains_pressure_response_and_native_material_work() {
    with_cx(|_, cx| {
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
        let velocity = vec![[0.2, 0.0, 0.0]; nodes.len()];
        let solid = MixedViscoTetSolid::try_new(
            &nodes,
            &elements,
            MixedViscoLawField::Uniform(&material),
            VolumetricTreatment::ExactIncompressible,
            TetAssemblyBudget::standard(),
            4,
            64,
            cx,
        )
        .unwrap();
        let motion =
            MixedViscoTetMotion::try_new(solid, &[6.0], &velocity, 128 * 1024, cx).unwrap();
        // Both witnesses lie on the physical x=0 face. Its outward normal is
        // -x and the authored tangent is +y, so the second tangent is -z.
        let points =
            [[0.5, 0.0, 0.25, 0.25], [0.25, 0.0, 0.5, 0.25]].map(|barycentric| Tet10ImpulsePoint {
                tetrahedron: 0,
                barycentric,
            });
        let contacts: [_; 2] = core::array::from_fn(|body| InelasticMixedContact {
            id: 40 + body as u64,
            point: RigidImpactPoint {
                body,
                arm_body: Vec3::new(0.0, 0.0, 0.2),
            },
            counterpart: FiniteSolidContactCounterpart::Solid(points[body]),
            normal_world: Vec3::new(-1.0, 0.0, 0.0),
            tangent_world: Vec3::new(0.0, 1.0, 0.0),
            activation_tolerance_m: 1e-12,
            friction_coefficient: 0.0,
            normal_law: NormalVelocityLaw::NewtonRestitution { coefficient: 0.5 },
        });
        let bodies: [_; 2] = core::array::from_fn(|body| {
            let x = motion
                .point_kinematics(points[body], cx)
                .unwrap()
                .position_m;
            let rigid = rigid_at(x, 1.2);
            RigidBodyState::new(
                rigid.pose(),
                Vec3::new(2.4, if body == 0 { 0.8 } else { -0.5 }, 0.6),
                Vec3::ZERO,
            )
            .unwrap()
        });
        let props = [properties(); 2];
        let controls = volume_controls();
        let event = mixed_visco_contact_candidate(
            &motion,
            &[],
            &bodies,
            &props,
            &contacts,
            None,
            ORIGIN,
            controls,
            cx,
        )
        .unwrap();
        // Independent rational elimination of the 30+4 mass/volume equations.
        // The normal block includes rigid 7/12 on its diagonal; the cross-axis
        // entry includes rigid -1/60 from the full off-diagonal inertia tensor.
        for (index, expected) in [
            (0, 6763.0 / 1088.0),
            (3, 8005.0 / 2176.0),
            (18, 8005.0 / 2176.0),
            (21, 24989.0 / 3264.0),
            (1, 2585.0 / 3264.0 - 1.0 / 60.0),
            (4, -1975.0 / 1632.0),
        ] {
            near(event.report.contact_mobility_per_kg[index], expected, 3e-10);
        }
        for (i, expected) in [4984896.0 / 28456609.0, 3180096.0 / 28456609.0]
            .into_iter()
            .enumerate()
        {
            near(event.report.contact.impulses[i][0], expected, 2e-10);
            near(event.report.contact.impulses[i][1], 0.0, 1e-12);
            near(event.report.contact.impulses[i][2], 0.0, 1e-12);
            near(
                event.report.contact.relative_velocities_after[i][0],
                0.5,
                2e-10,
            );
        }
        let loss = 2041248.0 / 28456609.0;
        near(event.report.contact.dissipated_energy, loss, 2e-10);
        near(event.report.kinetic_change_j, -loss, 2e-10);
        near(event.report.support_work_j, 0.0, 1e-12);
        assert!(
            event
                .report
                .solid_impulse
                .pressure_impulses_pa_s
                .iter()
                .any(|p| p.abs() > 1e-6)
        );
        volume_audit(&event.report);

        let circular = contacts.map(|c| InelasticMixedContact {
            normal_law: NormalVelocityLaw::NewtonRestitution { coefficient: 0.0 },
            ..c.with_friction(0.05)
        });
        let circle = mixed_visco_contact_candidate(
            &motion,
            &[],
            &bodies,
            &props,
            &circular,
            None,
            ORIGIN,
            controls,
            cx,
        )
        .unwrap();
        let equal = mixed_visco_elliptic_contact_candidate(
            &motion,
            &[],
            &bodies,
            &props,
            &circular.map(|c| {
                c.with_friction(EllipticFriction {
                    coefficients_mu: [0.05; 2],
                    principal_axis: [3.0, 4.0],
                })
            }),
            None,
            ORIGIN,
            controls,
            cx,
        )
        .unwrap();
        assert_eq!(equal.rigid_bodies, circle.rigid_bodies);
        assert_eq!(equal.solid.velocity_mps(), circle.solid.velocity_mps());
        assert_eq!(equal.report.contact.transformed, circle.report.contact);
        let free = mixed_visco_elliptic_contact_candidate(
            &motion,
            &[],
            &bodies,
            &props,
            &circular.map(|c| {
                c.with_friction(EllipticFriction {
                    coefficients_mu: [0.05, 0.0],
                    principal_axis: [1.0, 0.0],
                })
            }),
            None,
            ORIGIN,
            controls,
            cx,
        )
        .unwrap();
        for impulse in &free.report.contact.impulses {
            near(impulse[2], 0.0, 1e-13);
            assert!(impulse[1].abs() > 1e-4);
        }
        volume_audit(&circle.report);
        volume_audit(&free.report);
        for next in [&event.solid, &circle.solid, &free.solid] {
            assert_eq!(next.solid().pressure_pa(), motion.solid().pressure_pa());
            assert_eq!(
                next.solid().quadrature_states(),
                motion.solid().quadrature_states()
            );
            assert_eq!(
                next.solid().displacement_m(),
                motion.solid().displacement_m()
            );
            assert_eq!(next.solid().time_s(), motion.solid().time_s());
        }

        let density = [6.0];
        for exact in [true, false] {
            let definition = ConservativeTetDefinition {
                nodes_m: &nodes,
                elements: &elements,
                material: MixedLawField::Uniform(&elastic),
                reference_density_kg_m3: &density,
                formulation: if exact {
                    ConservativeFormulation::ExactIncompressible {
                        maximum_relative_volume_residual: 1e-9,
                        maximum_initial_relative_volume_rate_per_s: 1e-10,
                        maximum_point_volume_work_defect: 1e-12,
                    }
                } else {
                    ConservativeFormulation::FullMaterial
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
            };
            let source = definition
                .prepare_with_pressure(
                    &vec![[0.0; 3]; nodes.len()],
                    &velocity,
                    if exact { &[0.3; 4] } else { &[] },
                    cx,
                )
                .unwrap();
            let next = conservative_contact_candidate(
                &source,
                &[],
                &bodies,
                &props,
                &contacts,
                None,
                ORIGIN,
                controls,
                cx,
            )
            .unwrap();
            volume_audit(&next.report);
            if exact {
                assert!(next.report.solid_impulse.kkt_inertia.is_some());
                for (actual, expected) in next
                    .report
                    .contact_mobility_per_kg
                    .iter()
                    .zip(&event.report.contact_mobility_per_kg)
                {
                    near(*actual, *expected, 3e-10);
                }
                near(next.report.contact.dissipated_energy, loss, 2e-10);
            } else {
                assert!(next.report.solid_impulse.kkt_inertia.is_none());
                assert!(next.report.solid_impulse.pressure_impulses_pa_s.is_empty());
            }
            // The conservative owner's cumulative gate sees the contact event's
            // actual external work, while stored energy and its pressure stay fixed.
            near(
                next.solid.total_energy_j() - source.total_energy_j(),
                next.report.solid_impulse.external_work_j,
                2e-10,
            );
            assert_eq!(next.solid.stored_energy_j(), source.stored_energy_j());
            assert_eq!(
                next.solid.algorithmic_pressure_pa(),
                source.algorithmic_pressure_pa()
            );
            assert_eq!(next.solid.displacement_m(), source.displacement_m());
            assert_eq!(next.solid.time_s(), source.time_s());
            let elliptic = conservative_elliptic_contact_candidate(
                &source,
                &[],
                &bodies,
                &props,
                &contacts.map(|c| {
                    c.with_friction(EllipticFriction {
                        coefficients_mu: [0.0; 2],
                        principal_axis: [3.0, 4.0],
                    })
                }),
                None,
                ORIGIN,
                controls,
                cx,
            )
            .unwrap();
            assert_eq!(elliptic.rigid_bodies, next.rigid_bodies);
            assert_eq!(elliptic.solid.velocity_mps(), next.solid.velocity_mps());
            assert_eq!(elliptic.report.contact.transformed, next.report.contact);
            assert_eq!(source.velocity_mps(), velocity);
        }
        assert_eq!(motion.velocity_mps(), velocity);
    });
}
