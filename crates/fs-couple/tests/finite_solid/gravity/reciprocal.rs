//! Native tidal response, finite-source recoil and supported pair momentum.
#[path = "reciprocal/trajectory.rs"]
mod trajectory;

use super::*;
use fs_couple::finite_solid::gravity::reciprocal::{
    QuadrupoleSourceExtent, ReciprocalGravityKickControls,
};
use fs_mbd::gravity::{CentralGravity, QuadrupoleSourceGravity};
use fs_mbd::{MassProperties, Pose, RigidBodyState, UnitQuaternion};

fn kick_controls() -> ReciprocalGravityKickControls {
    ReciprocalGravityKickControls {
        frame: FRAME,
        gravitational_constant_m3_kg_s2: 2.0,
        duration_s: 0.02,
        balance_origin_m: [0.7, -0.8, 0.4],
        impulse: Default::default(),
    }
}
fn source_at(position: [f64; 3]) -> (RigidBodyState, MassProperties) {
    let inertia = MassProperties::new(3.0, Vec3::ZERO, Vec3::new(0.2, 0.3, 0.4)).unwrap();
    let body = RigidBodyState::new(
        Pose::new(vector(position), UnitQuaternion::IDENTITY).unwrap(),
        Vec3::new(0.15, -0.25, 0.1),
        Vec3::new(0.01, -0.02, 0.03),
    )
    .unwrap();
    (body, inertia)
}
const SOURCE_EXTENT: QuadrupoleSourceExtent = QuadrupoleSourceExtent {
    radius_m: 0.7,
    maximum_radius_ratio: 0.3,
};
fn extended_source() -> (RigidBodyState, MassProperties) {
    let (source, properties) = source_at([-3.0, 0.3, -1.0]);
    let properties = MassProperties::from_principal_frame(
        properties.mass(),
        Vec3::ZERO,
        properties.principal_inertia_body(),
        UnitQuaternion::from_axis_angle(Vec3::new(1.0, 2.0, -0.5), 0.6).unwrap(),
    )
    .unwrap();
    let pose = Pose::new(
        source.pose().position_world(),
        UnitQuaternion::from_axis_angle(Vec3::new(-2.0, 1.0, 0.5), 0.4).unwrap(),
    )
    .unwrap();
    (
        RigidBodyState::new(
            pose,
            source.linear_momentum_world(),
            source.angular_momentum_body(),
        )
        .unwrap(),
        properties,
    )
}
fn near_vector(a: [f64; 3], b: [f64; 3], tolerance: f64) {
    for c in 0..3 {
        near(a[c], b[c], tolerance);
    }
}

/// Newton's inverse-square field integrated independently on the two admitted
/// reference rules; no production field, quadrature, force or mass operation.
fn oracle(x: [[f64; 3]; 4], rho: [f64; 4], source: [f64; 3], graded: bool) -> ([[f64; 3]; 4], f64) {
    let c = kick_controls();
    let mut impulse = [[0.0; 3]; 4];
    let mut potential = 0.0;
    for q in 0..if graded { 8 } else { 4 } {
        let shape: [f64; 4] = std::array::from_fn(|a| {
            if graded {
                if q < 4 {
                    if a == q { 1.0 } else { 0.0 }
                } else if a == q - 4 {
                    0.0
                } else {
                    1.0 / 3.0
                }
            } else if a == q {
                (5.0 + 3.0 * 5.0_f64.sqrt()) / 20.0
            } else {
                (5.0 - 5.0_f64.sqrt()) / 20.0
            }
        });
        let weight = if graded {
            let density: f64 = (0..4).map(|a| shape[a] * rho[a]).sum();
            VOLUME * density * if q < 4 { 1.0 / 40.0 } else { 9.0 / 40.0 }
        } else {
            VOLUME * rho[0] / 4.0
        };
        let position: [f64; 3] = std::array::from_fn(|i| (0..4).map(|a| shape[a] * x[a][i]).sum());
        let r = vector(position).sub(vector(source));
        let radius = (r.x * r.x + r.y * r.y + r.z * r.z).sqrt();
        let acceleration = r.scale(-c.gravitational_constant_m3_kg_s2 * 3.0 / radius.powi(3));
        potential -= weight * c.gravitational_constant_m3_kg_s2 * 3.0 / radius;
        for a in 0..4 {
            for i in 0..3 {
                impulse[a][i] += c.duration_s * weight * shape[a] * components(acceleration)[i];
            }
        }
    }
    (impulse, potential)
}

/// Exact P1 consistent mass integral of rho*N_a*N_b, including graded rho.
fn mass(a: usize, b: usize, rho: [f64; 4]) -> f64 {
    VOLUME * (rho.iter().sum::<f64>() + rho[a] + rho[b]) / if a == b { 60.0 } else { 120.0 }
}

#[test]
fn extended_source_spin_is_the_discrete_potential_derivative_with_actual_pair_recoil() {
    with_cx(|_, cx| {
        let material = bind(&source_card(6.0, 1.0, false, false));
        let (source, properties) = extended_source();
        let c = kick_controls();
        for fixed in [vec![], vec![(0, 0.0), (1, 0.0), (2, 0.0)]] {
            let body = FiniteSolidBody::new(
                problem(&material, &fixed, &[]),
                &[[0.0; 3]; 4],
                &[[0.0; 3]; 4],
                0.0,
                Default::default(),
                cx,
            )
            .unwrap();
            let before = body.native_state().clone();
            let (next, recoiled, report) = body
                .reciprocal_quadrupole_gravity_kick_candidate(
                    source,
                    properties,
                    SOURCE_EXTENT,
                    c,
                    cx,
                )
                .unwrap();
            let spin = recoiled
                .angular_momentum_body()
                .sub(source.angular_momentum_body());
            assert!(spin.x.hypot(spin.y).hypot(spin.z) > 1e-7);
            let prepared = body.problem().prepare(cx).unwrap();
            let potential = |orientation| {
                let pose = Pose::new(source.pose().position_world(), orientation).unwrap();
                let gravity = QuadrupoleSourceGravity::new(
                    FRAME,
                    pose,
                    properties,
                    c.gravitational_constant_m3_kg_s2,
                    SOURCE_EXTENT,
                )
                .unwrap();
                let field = fs_couple::solid_gravity::bind_field(&gravity, FRAME);
                prepared
                    .gravity_force_and_potential(&[6.0], &[0.0; 12], &field, cx)
                    .unwrap()
                    .1
            };
            for axis in 0..3 {
                let h = 1e-5;
                let mut direction = [0.0; 3];
                direction[axis] = h;
                let plus = source
                    .pose()
                    .orientation()
                    .right_exp(vector(direction))
                    .unwrap();
                let minus = source
                    .pose()
                    .orientation()
                    .right_exp(vector(direction).scale(-1.0))
                    .unwrap();
                near(
                    components(spin)[axis] / c.duration_s,
                    -(potential(plus) - potential(minus)) / (2.0 * h),
                    2e-9,
                );
            }
            let p0 = body.momentum_about(c.balance_origin_m, cx).unwrap();
            let p1 = next.momentum_about(c.balance_origin_m, cx).unwrap();
            let source_dp = recoiled
                .linear_momentum_world()
                .sub(source.linear_momentum_world());
            let linear = vector(std::array::from_fn(|i| {
                p1.linear_kg_m_s[i] - p0.linear_kg_m_s[i]
            }))
            .add(source_dp);
            let angular = vector(std::array::from_fn(|i| {
                p1.angular_kg_m2_s[i] - p0.angular_kg_m2_s[i]
            }))
            .add(
                source
                    .pose()
                    .position_world()
                    .sub(vector(c.balance_origin_m))
                    .cross(source_dp),
            )
            .add(source.pose().orientation().rotate_body_to_world(spin));
            near_vector(
                components(linear),
                report.solid.support.linear_kg_m_s,
                3e-13,
            );
            near_vector(
                components(angular),
                report.solid.support.angular_kg_m2_s,
                3e-13,
            );
            if !fixed.is_empty() {
                assert!(
                    report
                        .solid
                        .support
                        .linear_kg_m_s
                        .iter()
                        .any(|x| x.abs() > 1e-6)
                );
            }
            near(
                report.kinetic_energy_change_j,
                next.kinetic_energy_j().unwrap() - body.kinetic_energy_j().unwrap()
                    + recoiled.kinetic_energy(properties).unwrap()
                    - source.kinetic_energy(properties).unwrap(),
                2e-14,
            );
            near(report.impulse_work_j, report.kinetic_energy_change_j, 2e-14);
            assert_eq!(body.native_state(), &before);
            assert_eq!(recoiled.pose(), source.pose());
            let (shifted, shifted_source, _) = body
                .reciprocal_quadrupole_gravity_kick_candidate(
                    source,
                    properties,
                    SOURCE_EXTENT,
                    ReciprocalGravityKickControls {
                        balance_origin_m: [-1.1, 2.3, -0.7],
                        ..c
                    },
                    cx,
                )
                .unwrap();
            assert_eq!(shifted.native_state(), next.native_state());
            assert_eq!(shifted_source, recoiled);
        }
    });
}

#[test]
fn isotropic_source_kick_retains_monopole_bits_and_extent_refusal_is_atomic() {
    with_cx(|_, cx| {
        let material = bind(&source_card(6.0, 1.0, false, false));
        let body = FiniteSolidBody::new(
            problem(&material, &[], &[]),
            &[[0.0; 3]; 4],
            &[[0.0; 3]; 4],
            0.0,
            Default::default(),
            cx,
        )
        .unwrap();
        let (source, _) = extended_source();
        let properties = MassProperties::new(3.0, Vec3::ZERO, Vec3::new(0.3, 0.3, 0.3)).unwrap();
        let c = kick_controls();
        let (monopole, monopole_source, mono) = body
            .reciprocal_gravity_kick_candidate(source, properties, c, cx)
            .unwrap();
        let (extended, extended_source, quad) = body
            .reciprocal_quadrupole_gravity_kick_candidate(source, properties, SOURCE_EXTENT, c, cx)
            .unwrap();
        assert_eq!(monopole.native_state(), extended.native_state());
        assert_eq!(monopole_source, extended_source);
        assert_eq!(mono.source, quad.source);
        assert_eq!(mono.potential_energy_j, quad.potential_energy_j);
        let before = body.native_state().clone();
        assert!(
            body.reciprocal_quadrupole_gravity_kick_candidate(
                source,
                properties,
                QuadrupoleSourceExtent {
                    radius_m: 3.0,
                    ..SOURCE_EXTENT
                },
                c,
                cx
            )
            .is_err()
        );
        assert_eq!(body.native_state(), &before);
    });
}

#[test]
fn reciprocal_tidal_kick_uses_actual_mass_and_recoils_with_origin_covariant_momentum() {
    with_cx(|_, cx| {
        let material = bind(&source_card(6.0, 1.0, false, false));
        let u = [
            [0.2, -0.1, 0.3],
            [0.35, -0.1, 0.3],
            [0.3, -0.2, 0.3],
            [0.2, -0.1, 0.5],
        ];
        let x = std::array::from_fn(|a| std::array::from_fn(|i| NODES[a][i] + u[a][i]));
        let v = [[0.2, -0.1, 0.3]; 4];
        let (source, inertia) = source_at([-2.0, 0.3, -0.7]);
        for graded in [false, true] {
            let rho = if graded {
                [3.0, 9.0, 12.0, 6.0]
            } else {
                [6.0; 4]
            };
            let p = problem(&material, &[], &[]);
            let body = if graded {
                FiniteSolidBody::new_with_reference_density_field(
                    p,
                    graded_density(rho),
                    &u,
                    &v,
                    0.4,
                    Default::default(),
                    cx,
                )
            } else {
                FiniteSolidBody::new(p, &u, &v, 0.4, Default::default(), cx)
            }
            .unwrap();
            let before = body.native_state().clone();
            let (expected, potential) =
                oracle(x, rho, components(source.pose().position_world()), graded);
            let total: [f64; 3] = std::array::from_fn(|i| expected.iter().map(|j| j[i]).sum());
            let c = kick_controls();
            let (next, recoiled, report) = body
                .reciprocal_gravity_kick_candidate(source, inertia, c, cx)
                .unwrap();
            near(report.potential_energy_j, potential, 2e-14);
            near_vector(
                components(
                    recoiled
                        .linear_momentum_world()
                        .sub(source.linear_momentum_world()),
                ),
                total.map(|j| -j),
                2e-16,
            );
            assert_eq!(
                recoiled.angular_momentum_body(),
                source.angular_momentum_body()
            );
            assert_eq!(recoiled.pose(), source.pose());
            assert_eq!(body.native_state(), &before);
            assert_eq!(next.displacement_m(), u);
            assert_eq!(next.native_state().t, before.t);
            assert_eq!(next.native_state().a, before.a);
            assert_eq!(next.native_state().history, before.history);
            let next_v = next.velocity_m_s();
            assert_ne!(next_v[0], next_v[1]); // resolved tidal loading, not a COM force
            let mut angular = Vec3::ZERO;
            let mut kinetic = 0.0;
            for a in 0..4 {
                let dp: [f64; 3] = std::array::from_fn(|i| {
                    (0..4)
                        .map(|b| mass(a, b, rho) * (next_v[b][i] - v[b][i]))
                        .sum()
                });
                near_vector(dp, expected[a], 2e-14);
                angular = angular.add(
                    vector(x[a])
                        .sub(vector(c.balance_origin_m))
                        .cross(vector(dp)),
                );
                for i in 0..3 {
                    kinetic += 0.5 * (v[a][i] + next_v[a][i]) * dp[i];
                }
            }
            let source_dp = recoiled
                .linear_momentum_world()
                .sub(source.linear_momentum_world());
            angular = angular.add(
                source
                    .pose()
                    .position_world()
                    .sub(vector(c.balance_origin_m))
                    .cross(source_dp),
            );
            near_vector(components(angular), [0.0; 3], 5e-14);
            let source_mean_v = source
                .center_of_mass_velocity_world(inertia)
                .unwrap()
                .add(recoiled.center_of_mass_velocity_world(inertia).unwrap())
                .scale(0.5);
            kinetic += source_mean_v.dot(source_dp);
            near(report.kinetic_energy_change_j, kinetic, 2e-14);
            near(report.impulse_work_j, kinetic, 2e-14);
            let shifted = ReciprocalGravityKickControls {
                balance_origin_m: [-1.3, 2.2, -0.6],
                ..c
            };
            let (shifted_body, shifted_source, shifted_report) = body
                .reciprocal_gravity_kick_candidate(source, inertia, shifted, cx)
                .unwrap();
            assert_eq!(next.native_state(), shifted_body.native_state());
            assert_eq!(recoiled, shifted_source);
            let expected_shift = vector(shifted.balance_origin_m)
                .sub(vector(c.balance_origin_m))
                .cross(vector(total));
            near_vector(
                std::array::from_fn(|i| {
                    report.solid.applied.angular_kg_m2_s[i]
                        - shifted_report.solid.applied.angular_kg_m2_s[i]
                }),
                components(expected_shift),
                2e-14,
            );

            // Force-only execution must retain the original full-tangent path's bits.
            let prepared = body.problem().prepare(cx).unwrap();
            let field = CentralGravity::new(FRAME, source.pose().position_world(), 6.0).unwrap();
            let field = fs_couple::solid_gravity::bind_field(&field, FRAME);
            let (full, light) = if graded {
                (
                    prepared
                        .gravity_load_with_nodal_density(&rho, u.as_flattened(), &field, cx)
                        .unwrap(),
                    prepared
                        .gravity_force_and_potential_with_nodal_density(
                            &rho,
                            u.as_flattened(),
                            &field,
                            cx,
                        )
                        .unwrap(),
                )
            } else {
                (
                    prepared
                        .gravity_load(&[6.0], u.as_flattened(), &field, cx)
                        .unwrap(),
                    prepared
                        .gravity_force_and_potential(&[6.0], u.as_flattened(), &field, cx)
                        .unwrap(),
                )
            };
            assert_eq!(full.force_n, light.0);
            assert_eq!(full.potential_energy_j, light.1);
        }
    });
}

#[test]
fn anchored_solid_keeps_support_impulse_and_source_recoil_separate() {
    with_cx(|_, cx| {
        let material = bind(&source_card(6.0, 1.0, false, false));
        let fixed = [(0, 0.0), (1, 0.0), (2, 0.0)];
        let body = FiniteSolidBody::new(
            problem(&material, &fixed, &[]),
            &[[0.0; 3]; 4],
            &[[0.0; 3]; 4],
            0.0,
            Default::default(),
            cx,
        )
        .unwrap();
        let (source, inertia) = source_at([-2.0, 0.3, -0.7]);
        let c = kick_controls();
        let (expected, _) = oracle(
            NODES,
            [6.0; 4],
            components(source.pose().position_world()),
            false,
        );
        let free: [f64; 3] = std::array::from_fn(|i| (1..4).map(|a| expected[a][i]).sum());
        let support: [f64; 3] = std::array::from_fn(|i| free[i] / 4.0 - expected[0][i]);
        let (next, recoiled, report) = body
            .reciprocal_gravity_kick_candidate(source, inertia, c, cx)
            .unwrap();
        assert_eq!(next.velocity_m_s()[0], [0.0; 3]);
        for a in 1..4 {
            near_vector(
                next.velocity_m_s()[a],
                std::array::from_fn(|i| 20.0 * (expected[a][i] - free[i] / 4.0)),
                2e-14,
            );
        }
        near_vector(report.solid.support.linear_kg_m_s, support, 2e-14);
        let actual_solid = next.momentum_about(c.balance_origin_m, cx).unwrap();
        let actual_pair = vector(actual_solid.linear_kg_m_s).add(
            recoiled
                .linear_momentum_world()
                .sub(source.linear_momentum_world()),
        );
        near_vector(components(actual_pair), support, 2e-14);
        let actual_angular = vector(actual_solid.angular_kg_m2_s).add(
            source
                .pose()
                .position_world()
                .sub(vector(c.balance_origin_m))
                .cross(
                    recoiled
                        .linear_momentum_world()
                        .sub(source.linear_momentum_world()),
                ),
        );
        near_vector(
            components(actual_angular),
            components(
                vector(NODES[0])
                    .sub(vector(c.balance_origin_m))
                    .cross(vector(support)),
            ),
            2e-14,
        );
        assert!(support.iter().any(|x| x.abs() > 1e-6));
        near_vector(report.linear_residual_kg_m_s, [0.0; 3], 2e-14);
        near_vector(report.angular_residual_kg_m2_s, [0.0; 3], 2e-14);
    });
}

#[test]
fn singular_field_unresolved_recoil_and_cancellation_publish_neither_owner() {
    with_cx(|gate, cx| {
        let material = bind(&source_card(6.0, 1.0, false, false));
        let body = FiniteSolidBody::new_with_reference_density_field(
            problem(&material, &[], &[]),
            graded_density([3.0, 9.0, 12.0, 6.0]),
            &[[0.0; 3]; 4],
            &[[0.0; 3]; 4],
            0.4,
            Default::default(),
            cx,
        )
        .unwrap();
        let before = body.native_state().clone();
        let (source, inertia) = source_at(NODES[0]); // admitted nodal rule samples this singularity
        assert!(
            body.reciprocal_gravity_kick_candidate(source, inertia, kick_controls(), cx)
                .is_err()
        );
        let (source, inertia) = source_at([-2.0, 0.3, -0.7]);
        let unresolved = RigidBodyState::new(
            source.pose(),
            Vec3::new(1e20, 1e20, 1e20),
            source.angular_momentum_body(),
        )
        .unwrap();
        assert!(matches!(
            body.reciprocal_gravity_kick_candidate(unresolved, inertia, kick_controls(), cx),
            Err(FiniteSolidError::ImpulseBalance { .. })
        ));
        assert_eq!(body.native_state(), &before);
        gate.request();
        assert!(
            body.reciprocal_gravity_kick_candidate(source, inertia, kick_controls(), cx)
                .is_err()
        );
        assert_eq!(body.native_state(), &before);
    });
}
