//! G0/G1/G3: native impulses and trajectories across independent contact islands.
#[path = "frictional_islands/connected.rs"]
mod connected;

use fs_contact::frictional_multicontact::{FrictionalContactConfig, NormalVelocityLaw};
use fs_exec::CancelGate;
use fs_mbd::frictional_impact::islands::{IslandContactError, apply_rigid_contact_islands};
use fs_mbd::frictional_impact::step::elliptic::step_rigid_elliptic_contact_islands;
use fs_mbd::frictional_impact::step::{
    RigidContactStepControls, RigidContactStepError, step_rigid_contact_islands,
};
use fs_mbd::frictional_impact::{
    RigidFrictionalContact, RigidFrictionalImpactError, RigidImpactControls,
    RigidImpactCounterpart, RigidImpactPoint, apply_rigid_frictional_impacts,
};
use fs_mbd::{Gravity, MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3, Wrench};

fn controls() -> RigidImpactControls {
    RigidImpactControls {
        maximum_bodies: 128,
        maximum_contacts: 64,
        geometry_tolerance_m: 1e-8,
        velocity_tolerance_mps: 2e-8,
        linear_momentum_tolerance_ns: 1e-9,
        angular_momentum_tolerance_nms: 1e-9,
        energy_tolerance_j: 2e-8,
        contact: FrictionalContactConfig::default(),
    }
}

fn step_controls() -> RigidContactStepControls {
    RigidContactStepControls {
        impact: controls(),
        drift_energy_tolerance_j: 1e-9,
        drift_angular_momentum_tolerance_nms: 1e-9,
        work_energy_tolerance_j: 2e-8,
        linear_momentum_tolerance_ns: 1e-9,
        angular_momentum_tolerance_nms: 1e-9,
    }
}

fn property() -> MassProperties {
    MassProperties::new(1.0, Vec3::ZERO, Vec3::new(1.0, 1.0, 1.0)).unwrap()
}

fn body(position: Vec3, momentum: Vec3) -> RigidBodyState {
    RigidBodyState::new(
        Pose::new(position, UnitQuaternion::IDENTITY).unwrap(),
        momentum,
        Vec3::ZERO,
    )
    .unwrap()
}

fn plane(index: usize, id: u64) -> RigidFrictionalContact {
    RigidFrictionalContact {
        id,
        point: RigidImpactPoint {
            body: index,
            arm_body: Vec3::ZERO,
        },
        other: RigidImpactCounterpart::PrescribedPlane {
            point_world: Vec3::ZERO,
            velocity_world: Vec3::ZERO,
        },
        normal_world: Vec3::new(0.0, 0.0, 1.0),
        tangent_world: Vec3::new(1.0, 0.0, 0.0),
        activation_tolerance_m: 1e-8,
        friction_coefficient: 0.25,
        normal_law: NormalVelocityLaw::NewtonRestitution { coefficient: 0.0 },
    }
}

fn close(a: f64, b: f64) {
    assert!((a - b).abs() < 2e-8, "{a:e} != {b:e}");
}

#[test]
fn g1_many_floor_impacts_match_original_native_solves_and_keep_canonical_order() {
    let n = 32;
    let properties = vec![property(); n];
    let initial: Vec<_> = (0..n)
        .map(|i| {
            body(
                Vec3::new(3.0 * i as f64, 0.0, 1.0),
                Vec3::new(1.0, 0.0, -2.0),
            )
        })
        .collect();
    let contacts: Vec<_> = (0..n)
        .map(|i| RigidFrictionalContact {
            point: RigidImpactPoint {
                body: i,
                arm_body: Vec3::new(0.0, 0.0, -1.0),
            },
            normal_law: NormalVelocityLaw::NewtonRestitution { coefficient: 0.5 },
            ..plane(i, (n - i) as u64)
        })
        .collect();
    // This was a total-scene limit even though the exact matrix is block diagonal.
    assert!(
        apply_rigid_frictional_impacts(
            &mut initial.clone(),
            &properties,
            &contacts,
            controls(),
            &CancelGate::new(),
        )
        .is_err()
    );

    let mut reference = initial.clone();
    let mut reference_loss = 0.0;
    for i in 0..n {
        let mut contact = contacts[i];
        contact.point.body = 0;
        let receipt = apply_rigid_frictional_impacts(
            &mut reference[i..=i],
            &properties[i..=i],
            &[contact],
            controls(),
            &CancelGate::new(),
        )
        .unwrap();
        reference_loss += receipt.contact.dissipated_energy;
    }
    let mut original_report = None;
    for reverse in [false, true] {
        let mut input = contacts.clone();
        if reverse {
            input.reverse();
        }
        let mut states = initial.clone();
        let report = apply_rigid_contact_islands(
            &mut states,
            &properties,
            &input,
            controls(),
            &CancelGate::new(),
        )
        .unwrap();
        assert_eq!(states, reference);
        assert_eq!(
            report.contact.islands.len(),
            n,
            "a shared floor must not join bodies"
        );
        assert_eq!(report.contact_ids, (1..=n as u64).collect::<Vec<_>>());
        for (slot, island) in report.contact.islands.iter().enumerate() {
            assert_eq!(island.body_indices, [n - slot - 1]);
            assert_eq!(island.contact_indices.len(), 1);
            assert_eq!(
                input[island.contact_indices[0]].id,
                report.contact_ids[slot]
            );
            assert_eq!(
                island.impact.contact.impulses[0],
                report.contact.impulses[slot]
            );
        }
        close(report.contact.dissipated_energy, reference_loss);
        close(reference_loss, 1.75 * n as f64);
        close(report.total_energy_residual_j, 0.0);
        if let Some((impulses, velocities)) = &original_report {
            assert_eq!(&report.contact.impulses, impulses);
            assert_eq!(&report.contact.relative_velocities_after, velocities);
        } else {
            original_report = Some((
                report.contact.impulses,
                report.contact.relative_velocities_after,
            ));
        }
    }
}

fn joined_scene() -> (Vec<RigidBodyState>, Vec<RigidFrictionalContact>) {
    let states = vec![
        body(Vec3::ZERO, Vec3::new(0.0, 0.0, -2.0)),
        body(Vec3::ZERO, Vec3::new(0.0, 0.0, 1.0)),
        body(Vec3::new(4.0, 0.0, 0.0), Vec3::new(0.0, 0.0, -1.0)),
    ];
    let contacts = vec![
        plane(0, 40),
        plane(2, 20),
        RigidFrictionalContact {
            other: RigidImpactCounterpart::Body(RigidImpactPoint {
                body: 1,
                arm_body: Vec3::ZERO,
            }),
            ..plane(0, 10)
        },
    ];
    (states, contacts)
}

#[test]
fn g3_shared_bodies_join_contacts_even_when_global_ids_interleave_islands() {
    let (initial, contacts) = joined_scene();
    let properties = vec![property(); initial.len()];
    let mut expected = initial.clone();
    let native = apply_rigid_frictional_impacts(
        &mut expected,
        &properties,
        &contacts,
        controls(),
        &CancelGate::new(),
    )
    .unwrap();
    let mut states = initial;
    let report = apply_rigid_contact_islands(
        &mut states,
        &properties,
        &contacts,
        controls(),
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(report.contact_ids, [10, 20, 40]);
    assert_eq!(report.contact.islands.len(), 2);
    assert_eq!(report.contact.islands[0].body_indices, [0, 1]);
    assert_eq!(report.contact.islands[0].contact_indices, [2, 0]);
    assert_eq!(report.contact.islands[1].body_indices, [2]);
    for (actual, expected) in states.iter().zip(expected) {
        close(
            actual.linear_momentum_world().z,
            expected.linear_momentum_world().z,
        );
    }
    for (actual, expected) in report.contact.impulses.iter().zip(native.contact.impulses) {
        for i in 0..3 {
            close(actual[i], expected[i]);
        }
    }
    close(
        report.contact.dissipated_energy,
        native.contact.dissipated_energy,
    );
}

#[test]
fn g0_connected_island_limit_and_late_native_failure_publish_nothing() {
    let (initial, contacts) = joined_scene();
    let properties = vec![property(); initial.len()];
    let mut limit = controls();
    limit.contact.maximum_contacts = 1;
    let mut states = initial.clone();
    assert!(matches!(
        apply_rigid_contact_islands(
            &mut states,
            &properties,
            &contacts,
            limit,
            &CancelGate::new()
        ),
        Err(RigidFrictionalImpactError::Contact(
            IslandContactError::ContactLimit {
                first_contact_id: 10,
                contacts: 2,
                maximum: 1,
            }
        ))
    ));
    assert_eq!(states, initial);

    // The first independent impact is admissible; a later material-law refusal
    // must discard it rather than partially committing the scene.
    let contacts = [
        plane(0, 10),
        RigidFrictionalContact {
            friction_coefficient: f64::NAN,
            ..plane(2, 20)
        },
    ];
    assert!(matches!(
        apply_rigid_contact_islands(
            &mut states,
            &properties,
            &contacts,
            controls(),
            &CancelGate::new()
        ),
        Err(RigidFrictionalImpactError::Contact(
            IslandContactError::Solve {
                first_contact_id: Some(20),
                ..
            }
        ))
    ));
    assert_eq!(states, initial);
}

#[test]
fn g0_global_angular_allowance_does_not_grow_with_the_number_of_islands() {
    // Tangential forces acting across each small positive normal gap produce
    // a real gap cross impulse couple. Each one is inside the chosen numerical
    // allowance; their sum is not. No contact point may be silently relocated.
    let gap = 1e-8;
    let initial: Vec<_> = (0..2)
        .flat_map(|i| {
            [
                body(
                    Vec3::new(4.0 * i as f64, 0.0, gap),
                    Vec3::new(1.0, 0.0, -1.0),
                ),
                body(Vec3::new(4.0 * i as f64, 0.0, 0.0), Vec3::ZERO),
            ]
        })
        .collect();
    let properties = vec![property(); initial.len()];
    let contacts: Vec<_> = (0..2)
        .map(|i| RigidFrictionalContact {
            other: RigidImpactCounterpart::Body(RigidImpactPoint {
                body: 2 * i + 1,
                arm_body: Vec3::ZERO,
            }),
            friction_coefficient: 1.0,
            ..plane(2 * i, i as u64)
        })
        .collect();
    let mut config = controls();
    config.angular_momentum_tolerance_nms = 6e-9;
    for i in 0..2 {
        let mut states = initial[2 * i..2 * i + 2].to_vec();
        let mut contact = contacts[i];
        contact.point.body = 0;
        contact.other = RigidImpactCounterpart::Body(RigidImpactPoint {
            body: 1,
            arm_body: Vec3::ZERO,
        });
        let report = apply_rigid_frictional_impacts(
            &mut states,
            &properties[..2],
            &[contact],
            config,
            &CancelGate::new(),
        )
        .unwrap();
        close(report.angular_momentum_residual_nms.y, -0.5 * gap);
    }
    let mut states = initial.clone();
    let error = apply_rigid_contact_islands(
        &mut states,
        &properties,
        &contacts,
        config,
        &CancelGate::new(),
    )
    .unwrap_err();
    let RigidFrictionalImpactError::Audit(report) = error else {
        panic!("wrong refusal: {error:?}")
    };
    assert!(report.angular_momentum_residual_nms.y.abs() > config.angular_momentum_tolerance_nms);
    assert_eq!(report.contact.islands.len(), 2);
    assert_eq!(states, initial);
}

#[test]
fn g0_empty_scene_global_caps_duplicate_ids_and_cancellation() {
    let initial = [body(Vec3::ZERO, Vec3::new(1.0, 0.0, -1.0))];
    let properties = [property()];
    let mut states = initial;
    let report = apply_rigid_contact_islands(
        &mut states,
        &properties,
        &[],
        controls(),
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(states, initial);
    assert!(report.contact.islands.is_empty());
    assert_eq!(report.kinetic_before_j, report.kinetic_after_j);
    assert_eq!(report.contact.dissipated_energy, 0.0);
    let mut config = controls();
    config.maximum_contacts = 1;
    assert!(
        apply_rigid_contact_islands(
            &mut states,
            &properties,
            &[plane(0, 1), plane(0, 2)],
            config,
            &CancelGate::new(),
        )
        .is_err()
    );
    assert!(
        apply_rigid_contact_islands(
            &mut states,
            &properties,
            &[plane(0, 1), plane(0, 1)],
            controls(),
            &CancelGate::new(),
        )
        .is_err()
    );
    let gate = CancelGate::new();
    gate.request();
    assert!(matches!(
        apply_rigid_contact_islands(&mut states, &properties, &[plane(0, 1)], controls(), &gate),
        Err(RigidFrictionalImpactError::Cancelled)
    ));
    assert_eq!(states, initial);
}

#[test]
fn g1_directional_friction_uses_the_same_independent_native_islands() {
    use fs_contact::frictional_multicontact::elliptic::EllipticFriction;
    let n = 24;
    let properties = vec![property(); n];
    let mut states: Vec<_> = (0..n)
        .map(|i| {
            body(
                Vec3::new(3.0 * i as f64, 0.0, 0.0),
                Vec3::new(1.0, 1.0, -1.0),
            )
        })
        .collect();
    let contacts: Vec<_> = (0..n)
        .map(|i| {
            let original = plane(i, i as u64);
            RigidFrictionalContact {
                friction_coefficient: EllipticFriction {
                    coefficients_mu: [0.2, 0.0],
                    principal_axis: [1.0, 0.0],
                },
                id: original.id,
                point: original.point,
                other: original.other,
                normal_world: original.normal_world,
                tangent_world: original.tangent_world,
                activation_tolerance_m: original.activation_tolerance_m,
                normal_law: original.normal_law,
            }
        })
        .collect();
    let report = step_rigid_elliptic_contact_islands(
        &mut states,
        &properties,
        &vec![Wrench::ZERO; n],
        Gravity::ZERO,
        &contacts,
        0.1,
        step_controls(),
        |path| {
            // Every supplied point stays on the plane: zero angular motion,
            // zero normal velocity and a constant tangential free drift.
            assert!(path.post_contact.iter().all(|body| {
                body.angular_momentum_body() == Vec3::ZERO && body.linear_momentum_world().z == 0.0
            }));
            Ok::<_, ()>(())
        },
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(report.impact.contact.islands.len(), n);
    for (i, state) in states.into_iter().enumerate() {
        close(state.linear_momentum_world().x, 0.8);
        close(state.linear_momentum_world().y, 1.0);
        close(state.linear_momentum_world().z, 0.0);
        close(state.pose().position_world().x, 3.0 * i as f64 + 0.08);
        close(state.pose().position_world().y, 0.1);
        close(state.pose().position_world().z, 0.0);
    }
    close(report.impact.contact.dissipated_energy, n as f64 * 0.68);
    close(report.total_energy_residual_j, 0.0);
}

#[test]
fn g1_island_force_contact_drift_and_late_path_cancellation_are_atomic() {
    let n = 24;
    let properties = vec![property(); n];
    let initial: Vec<_> = (0..n)
        .map(|i| {
            body(
                Vec3::new(3.0 * i as f64, 0.0, 0.0),
                Vec3::new(0.0, 0.0, -1.0),
            )
        })
        .collect();
    let contacts: Vec<_> = (0..n).map(|i| plane(i, i as u64)).collect();
    let controls = step_controls();
    for cancel_at_path in [false, true] {
        let gate = CancelGate::new();
        let mut states = initial.clone();
        let result = step_rigid_contact_islands(
            &mut states,
            &properties,
            &vec![Wrench::ZERO; n],
            Gravity::ZERO,
            &contacts,
            0.1,
            controls,
            |path| {
                // All post-impact velocities vanish and this exact drift is
                // constant: the plane gap is zero throughout the whole path.
                assert_eq!(path.post_contact, path.endpoint);
                assert!(
                    path.post_contact
                        .iter()
                        .all(|b| b.linear_momentum_world() == Vec3::ZERO)
                );
                if cancel_at_path {
                    gate.request();
                }
                Ok::<_, ()>(())
            },
            &gate,
        );
        if cancel_at_path {
            assert!(matches!(result, Err(RigidContactStepError::Cancelled)));
            assert_eq!(states, initial);
        } else {
            let report = result.unwrap();
            assert_eq!(report.impact.contact.islands.len(), n);
            close(report.total_energy_residual_j, 0.0);
            close(report.impact.contact.dissipated_energy, 0.5 * n as f64);
        }
    }
}
