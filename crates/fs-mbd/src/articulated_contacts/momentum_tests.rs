//! Native impulse and contact-solve oracles for simultaneous momentum transfer.

use super::*;
use crate::articulated::{JointModel, Link, SpatialInertia, free_floating_momentum};
use fs_ga::{Mat3, Se3, So3};

fn body(mass: f64) -> ArticulatedModel {
    let inertia = SpatialInertia::new(
        mass,
        Vec3::default(),
        Mat3 {
            m: [mass, 0.0, 0.0, 0.0, mass, 0.0, 0.0, 0.0, mass],
        },
    )
    .unwrap();
    ArticulatedModel::new(vec![Link::new(
        "body",
        None,
        Se3::identity(),
        JointModel::FIXED,
        inertia,
    )])
    .unwrap()
}

fn system(
    model: &ArticulatedModel,
    position: Vec3,
    velocity: Vec3,
) -> ArticulatedContactSystem<'_> {
    ArticulatedContactSystem {
        model,
        base: FreeFloatingBaseState::new(
            Se3::from_parts(So3::identity(), position).unwrap(),
            Twist::new(Vec3::default(), velocity),
        ),
        joint_position: &[],
        joint_velocity: &[],
    }
}

/// Each group has opposing large tangential impulses with a small remainder.
/// At y=+1 the resultant is Jx=1, Lz=-1; at y=-1 it is Jx=2, Lz=2.
/// The unit body's unique endpoint is vx=3, wz=1, with kinetic gain 5 J.
/// Supports at the two midpoint point speeds (1 and 2 m/s) supply exactly 5 J.
#[test]
fn simultaneous_point_resultants_retain_small_force_and_couple() {
    let model = body(1.0);
    let large = 2.0_f64.powi(56);
    let impulses = [large, 1.0, -large, large, 2.0, -large];
    for boost in [0.0, 2.0_f64.powi(30)] {
        let state = system(
            &model,
            Vec3::new(2.0, -3.0, 4.0),
            Vec3::new(0.0, 0.0, boost),
        );
        let systems = [state];
        let workspace = EventWorkspace::new(&systems, &mut || false).unwrap();
        let initial = free_floating_momentum(&model, state.base, &[], &[]).unwrap();
        for order in [[0, 1, 2, 3, 4, 5], [5, 4, 3, 2, 1, 0], [0, 3, 1, 4, 2, 5]] {
            let mut prepared = Vec::new();
            let mut speeds = Vec::new();
            for index in order {
                let (arm_y, support_speed) = if index < 3 { (1.0, 1.0) } else { (-1.0, 2.0) };
                let point = LinkContactPoint {
                    system_index: 0,
                    link_index: 0,
                    point_body: Vec3::new(0.0, arm_y, 0.0),
                };
                let axis = Vec3::new(1.0, 0.0, 0.0);
                let (first, world_point) = prepare_point(point, axis, &workspace.before).unwrap();
                prepared.push(Prepared {
                    first,
                    second: None,
                    axis_world: axis,
                    world_point,
                    plane_speed: support_speed,
                    minimum_speed: None,
                });
                speeds.push(support_speed);
            }
            let ordered_impulses = order.map(|index| impulses[index]);
            let law = || EventResponse {
                impulses: &ordered_impulses,
                predicted_speeds: &speeds,
                stored_energy_change: 0.0,
                dissipated_energy: 0.0,
                velocity_tolerance: 1e-10,
                energy_tolerance: 1e-10,
            };
            let mut polls = 0;
            let result = workspace
                .clone()
                .apply(&systems, &prepared, law(), &mut || {
                    polls += 1;
                    false
                })
                .unwrap();
            let response = &result.responses[0];
            assert_eq!(
                response.base_after.twist_body.linear,
                Vec3::new(3.0, 0.0, boost)
            );
            assert_eq!(
                response.base_after.twist_body.angular,
                Vec3::new(0.0, 0.0, 1.0)
            );
            assert_eq!(
                response.linear_momentum_change_world,
                Vec3::new(3.0, 0.0, 0.0)
            );
            assert_eq!(
                response.angular_momentum_change_world,
                Vec3::new(0.0, 12.0, 10.0)
            );
            assert_eq!(result.kinetic_energy_change, 5.0);
            assert_eq!(result.plane_work, 5.0);
            assert_eq!(result.energy_balance_residual, 0.0);
            assert_eq!(result.linear_momentum_residual, Vec3::default());
            assert_eq!(result.angular_momentum_residual, Vec3::default());
            assert_eq!(result.endpoint_velocity_residual, 0.0);
            assert_eq!(
                response.base_after.world_from_base,
                state.base.world_from_base
            );
            if boost != 0.0 {
                assert_eq!(
                    response.momentum_before.kinetic_energy,
                    response.momentum_after.kinetic_energy
                );
            }

            // Refuse during resultant accumulation and at final admission;
            // neither path publishes a partial native endpoint.
            for stop in [1, polls / 2, polls] {
                let mut calls = 0;
                let refused = workspace
                    .clone()
                    .apply(&systems, &prepared, law(), &mut || {
                        calls += 1;
                        calls == stop
                    });
                assert!(matches!(refused, Err(ArticulatedContactError::Cancelled)));
                assert_eq!(
                    free_floating_momentum(&model, state.base, &[], &[]).unwrap(),
                    initial
                );
            }
        }
    }
}

/// Three genuine normal-contact solves share one event. Elastic moving-wall
/// reflections supply +2^55 and -2^55 J; the middle inelastic impact receives
/// 1 J, keeps 0.5 J and dissipates 0.5 J. Both residuals and signed work must
/// retain those small quantities, independent of system/contact order or a
/// common transverse boost that changes no contact speed.
#[test]
fn moving_support_work_cancellation_retains_one_joule_through_the_contact_solver() {
    let model = body(2.0_f64.powi(54));
    let slow = 2.0_f64.powi(-27);
    let large = 2.0_f64.powi(55);
    for boost in [0.0, 2.0_f64.powi(20)] {
        for order in [[0, 1, 2], [2, 1, 0], [0, 2, 1]] {
            let velocities = [0.0, 0.0, 2.0];
            let support_speeds = [1.0, slow, 1.0];
            let signs = [1.0, 1.0, -1.0];
            let restitution = [1.0, 0.0, 1.0];
            let systems: Vec<_> = order
                .iter()
                .map(|&index| {
                    system(
                        &model,
                        Vec3::new(3.0 * index as f64, 2.0, 0.0),
                        Vec3::new(velocities[index], boost, 0.0),
                    )
                })
                .collect();
            let contacts: Vec<_> = order
                .iter()
                .enumerate()
                .map(|(system_index, &index)| ArticulatedNormalContact {
                    point: LinkContactPoint {
                        system_index,
                        link_index: 0,
                        point_body: Vec3::default(),
                    },
                    other: ContactCounterpart::PrescribedPlane {
                        point_world: Vec3::new(3.0 * index as f64, 2.0, 0.0),
                        velocity_world: Vec3::new(support_speeds[index], boost, 0.0),
                    },
                    normal_world: Vec3::new(signs[index], 0.0, 0.0),
                    activation_tolerance: 1e-10,
                    restitution: restitution[index],
                })
                .collect();
            let result = resolve_articulated_contacts(
                &systems,
                &contacts,
                MultiImpactConfig::default(),
                || false,
            )
            .unwrap();
            for (slot, index) in order.into_iter().enumerate() {
                assert_eq!(
                    result.normal_impulses[slot],
                    [large, 2.0_f64.powi(27), large][index]
                );
                assert_eq!(
                    result.responses[slot].kinetic_energy_change,
                    [large, 0.5, -large][index]
                );
                assert_eq!(
                    result.responses[slot].base_after.twist_body.linear,
                    Vec3::new([2.0, slow, 0.0][index], boost, 0.0)
                );
            }
            assert_eq!(result.plane_work, 1.0);
            assert_eq!(result.dissipated_energy, 0.5);
            assert_eq!(result.energy_balance_residual, 0.0);
            assert_eq!(result.linear_momentum_residual, Vec3::default());
            assert_eq!(result.angular_momentum_residual, Vec3::default());
        }
    }
}
