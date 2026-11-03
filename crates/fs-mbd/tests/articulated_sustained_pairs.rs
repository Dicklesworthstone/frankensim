//! G1/G3/G4: finite sustained exchange, native support work and reciprocal geometry.

use fs_contact::sustained::{NaturalMapScales, solve::SolveControls};
use fs_ga::{Mat3, Se3, So3, Twist, Vec3, Wrench};
use fs_mbd::articulated::contact_geometry::ContactGeometryError;
use fs_mbd::articulated::motion::prescribed::{
    PrescribedBaseMotionError, PrescribedBaseMotionState, PrescribedBasePath,
};
use fs_mbd::articulated::motion::{FreeFloatingMotionLoads, FreeFloatingMotionState};
use fs_mbd::articulated::{
    ArticulatedModel, BaseState, FreeFloatingBaseState, JointModel, Link, SpatialInertia,
};
use fs_mbd::articulated_contacts::step::sustained::{
    ArticulatedSustainedContact, ArticulatedSustainedContactPath, ArticulatedSustainedControls,
    ArticulatedSustainedError, ArticulatedSustainedMotion, ArticulatedSustainedStepCandidate,
    step_sustained_articulated_contacts,
};
use fs_mbd::articulated_contacts::trajectory::mixed::{
    MixedArticulatedTrajectoryState as State, MixedArticulatedTrajectorySystem as System,
};
use fs_mbd::articulated_contacts::{ContactCounterpart, LinkContactPoint};
use std::cell::Cell;

fn inertia() -> SpatialInertia {
    SpatialInertia::new(2.0, Vec3::default(), Mat3::identity()).unwrap()
}
fn body() -> ArticulatedModel {
    ArticulatedModel::new(vec![Link::new(
        "finite body",
        None,
        Se3::identity(),
        JointModel::FIXED,
        inertia(),
    )])
    .unwrap()
}
fn cartesian() -> ArticulatedModel {
    ArticulatedModel::new(vec![
        Link::new(
            "x carriage",
            None,
            Se3::identity(),
            JointModel::prismatic(Vec3::new(1.0, 0.0, 0.0), None).unwrap(),
            inertia(),
        ),
        Link::new(
            "y carriage",
            Some(0),
            Se3::identity(),
            JointModel::prismatic(Vec3::new(0.0, 1.0, 0.0), None).unwrap(),
            inertia(),
        ),
        Link::new(
            "z shoe",
            Some(1),
            Se3::identity(),
            JointModel::prismatic(Vec3::new(0.0, 0.0, 1.0), None).unwrap(),
            inertia(),
        ),
    ])
    .unwrap()
}
fn pose(position: Vec3) -> Se3 {
    Se3::from_parts(So3::identity(), position).unwrap()
}
fn free(position: Vec3, velocity: Vec3) -> State {
    State::Free(FreeFloatingMotionState {
        base: FreeFloatingBaseState::new(pose(position), Twist::new(Vec3::default(), velocity)),
        joint_position: vec![],
        joint_velocity: vec![],
    })
}
fn supported(position: Vec3, velocity: Vec3, dofs: usize) -> State {
    State::Prescribed(PrescribedBaseMotionState {
        time_s: 0.0,
        base: BaseState::prescribed(
            pose(position),
            Twist::new(Vec3::default(), velocity),
            Twist::zero(),
        ),
        joint_position: vec![0.0; dofs],
        joint_velocity: vec![0.0; dofs],
    })
}
fn controls() -> ArticulatedSustainedControls {
    ArticulatedSustainedControls {
        duration: 0.125,
        angular_origin_world: Vec3::new(3.0, -2.0, 1.0),
        solve: SolveControls {
            difference_step_n: 1.0,
            maximum_response_evaluations: 512,
            ..SolveControls::default()
        },
        loaded_gap_tolerance_m: 1e-9,
        force_tolerance_n: 1e-8,
        normal_velocity_tolerance_m_s: 1e-8,
        complementarity_tolerance_j: 1e-8,
        energy_tolerance_j: 1e-8,
        linear_momentum_tolerance_ns: 1e-8,
        angular_momentum_tolerance_nms: 1e-8,
    }
}
fn pair(id: u64, a: usize, b: usize, link: usize, mu: f64) -> ArticulatedSustainedContact {
    ArticulatedSustainedContact {
        id,
        point: LinkContactPoint {
            system_index: a,
            link_index: link,
            point_body: Vec3::new(0.0, -1.0, 0.0),
        },
        other: ContactCounterpart::Body(LinkContactPoint {
            system_index: b,
            link_index: link,
            point_body: Vec3::new(0.0, 1.0, 0.0),
        }),
        normal_world: Vec3::new(0.0, 1.0, 0.0),
        tangent_world: Vec3::new(1.0, 0.0, 0.0),
        friction_coefficient: mu,
        scales: NaturalMapScales {
            normal_n_per_m: 128.0,
            tangent_n_per_m: 128.0,
        },
    }
}
fn support_path(
    _: usize,
    state: &PrescribedBaseMotionState,
    duration: f64,
) -> Result<PrescribedBasePath, PrescribedBaseMotionError> {
    let at = |elapsed| {
        BaseState::prescribed(
            pose(
                state.base.world_from_base.translation()
                    + state.base.twist_body.linear.scale(elapsed),
            ),
            state.base.twist_body,
            Twist::zero(),
        )
    };
    Ok(PrescribedBasePath {
        midpoint: at(duration * 0.5),
        endpoint: at(duration),
    })
}
fn manufactured_translation(
    path: &ArticulatedSustainedContactPath<'_>,
) -> Result<(), &'static str> {
    // These manufactured bodies have fixed orientation, Cartesian joints and
    // constant original loads. Material trajectories are quadratic. The loaded
    // fixtures have identical initial witness velocities and accelerations, so
    // witnesses coincide throughout; the unloaded fixture separates outward.
    // No other surfaces are present. Those family facts, not stage sampling,
    // establish complete-path/frame/contact-mode applicability for THESE cases.
    for state in path.endpoint {
        let twist = match state {
            State::Free(state) => state.base.twist_body,
            State::Prescribed(state) => state.base.twist_body,
        };
        if twist.angular.norm() > 1e-12 {
            return Err("left the manufactured translating family");
        }
    }
    Ok(())
}
fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-9 * expected.abs().max(1.0),
        "{actual:e} != {expected:e}"
    );
}
fn vector(actual: Vec3, expected: Vec3) {
    close(actual.x, expected.x);
    close(actual.y, expected.y);
    close(actual.z, expected.z);
}
fn audit<P>(step: &ArticulatedSustainedStepCandidate<P>) {
    close(step.kinetic_work_residual_j, 0.0);
    close(step.contact_work_quadrature_defect_j, 0.0);
    vector(step.linear_momentum_residual_world, Vec3::default());
    vector(
        step.angular_momentum_residual_about_origin_world,
        Vec3::default(),
    );
    for pair in &step.pair_exchanges {
        vector(
            pair.first_linear_impulse_world + pair.second_linear_impulse_world,
            Vec3::default(),
        );
        vector(
            pair.first_angular_impulse_about_origin_world
                + pair.second_angular_impulse_about_origin_world,
            Vec3::default(),
        );
        vector(pair.linear_impulse_residual_world, Vec3::default());
        vector(pair.angular_impulse_residual_world, Vec3::default());
        close(pair.first_work_j + pair.second_work_j, pair.relative_work_j);
    }
    assert!(step.solve.paths_m.iter().all(|path| path[0] >= 0.0));
}

#[test]
fn free_finite_bodies_share_push_and_reciprocal_spatial_impulses() {
    let models = [body(), body()];
    let states = [
        free(Vec3::new(2.0, 1.0, 0.0), Vec3::default()),
        free(Vec3::new(2.0, -1.0, 0.0), Vec3::default()),
    ];
    let forces = [
        [Wrench::default()],
        [Wrench::new(Vec3::default(), Vec3::new(0.0, 16.0, 0.0))],
    ];
    for order in [[0, 1], [1, 0]] {
        let systems: Vec<_> = order
            .iter()
            .map(|&i| System {
                model: &models[i],
                state: &states[i],
                loads: FreeFloatingMotionLoads {
                    generalized_force: &[],
                    external_body_wrench: &forces[i],
                    gravity_world: Vec3::default(),
                },
            })
            .collect();
        let upper = order.iter().position(|&i| i == 0).unwrap();
        let lower = 1 - upper;
        let step = step_sustained_articulated_contacts(
            &systems,
            &[pair(40, upper, lower, 0, 0.0)],
            &[[0.0; 3]],
            0.0,
            controls(),
            |_, _, _| panic!("free bodies have no prescribed path"),
            manufactured_translation,
            || false,
        )
        .unwrap();
        // Two genuine masses of 2 kg share 16 N: a=4 m/s², N=8 N.
        close(step.solve.forces_n[0][0], 8.0);
        assert!(step.solve.nonlinear.iterations > 0);
        for (index, after) in step.states_after.iter().enumerate() {
            let State::Free(after) = after else { panic!() };
            let State::Free(before) = &states[order[index]] else {
                panic!()
            };
            vector(after.base.twist_body.linear, Vec3::new(0.0, 0.5, 0.0));
            vector(
                after.base.world_from_base.translation()
                    - before.base.world_from_base.translation(),
                Vec3::new(0.0, 0.03125, 0.0),
            );
        }
        vector(step.linear_momentum_change_world, Vec3::new(0.0, 2.0, 0.0));
        vector(step.external_linear_impulse_world, Vec3::new(0.0, 2.0, 0.0));
        vector(
            step.external_angular_impulse_about_origin_world,
            Vec3::new(2.0, 0.0, -2.0),
        );
        close(step.kinetic_energy_change_j, 0.5);
        close(step.body_wrench_work_j, 0.5);
        close(step.frictional_dissipation_j, 0.0);
        close(step.plane_work_j, 0.0);
        assert_eq!(step.pair_exchanges.len(), 1);
        let exchange = &step.pair_exchanges[0];
        vector(
            exchange.first_linear_impulse_world,
            Vec3::new(0.0, 1.0, 0.0),
        );
        vector(
            exchange.first_angular_impulse_about_origin_world,
            Vec3::new(1.0, 0.0, -1.0),
        );
        close(exchange.first_work_j, 0.25);
        close(exchange.second_work_j, -0.25);
        close(exchange.relative_work_j, 0.0);
        audit(&step);
    }
}

#[test]
fn moving_three_body_column_preserves_support_or_plane_work_and_internal_exchange() {
    let models = [body(), body(), body()];
    let zero = [Wrench::default()];
    for prescribed_floor in [true, false] {
        let states = [
            free(Vec3::new(2.0, 3.0, 0.0), Vec3::new(0.0, 1.0, 0.0)),
            free(Vec3::new(2.0, 1.0, 0.0), Vec3::new(0.0, 1.0, 0.0)),
            if prescribed_floor {
                supported(Vec3::new(2.0, -1.0, 0.0), Vec3::new(0.0, 1.0, 0.0), 0)
            } else {
                free(Vec3::new(2.0, -1.0, 0.0), Vec3::new(0.0, 1.0, 0.0))
            },
        ];
        let systems: Vec<_> = models
            .iter()
            .zip(&states)
            .map(|(model, state)| System {
                model,
                state,
                loads: FreeFloatingMotionLoads {
                    generalized_force: &[],
                    external_body_wrench: &zero,
                    gravity_world: Vec3::new(0.0, -10.0, 0.0),
                },
            })
            .collect();
        let mut contacts = vec![pair(90, 0, 1, 0, 0.0), pair(10, 1, 2, 0, 0.0)];
        if !prescribed_floor {
            let mut floor = pair(50, 2, 0, 0, 0.0);
            floor.other = ContactCounterpart::PrescribedPlane {
                point_world: Vec3::new(2.0, -2.0, 0.0),
                velocity_world: Vec3::new(0.0, 1.0, 0.0),
            };
            contacts.push(floor);
        }
        for reversed in [false, true] {
            if reversed {
                contacts.reverse();
            }
            let warm: Vec<_> = contacts
                .iter()
                .map(|contact| {
                    [
                        match contact.id {
                            90 => 20.0,
                            10 => 40.0,
                            50 => 60.0,
                            _ => unreachable!(),
                        },
                        0.0,
                        0.0,
                    ]
                })
                .collect();
            let calls = Cell::new(0);
            let step = step_sustained_articulated_contacts(
                &systems,
                &contacts,
                &warm,
                0.0,
                controls(),
                |index, state, duration| {
                    assert!(prescribed_floor);
                    assert_eq!(index, 2);
                    calls.set(calls.get() + 1);
                    support_path(index, state, duration)
                },
                manufactured_translation,
                || false,
            )
            .unwrap();
            assert_eq!(calls.get(), usize::from(prescribed_floor));
            assert_eq!(
                step.pair_exchanges
                    .iter()
                    .map(|p| p.contact_id)
                    .collect::<Vec<_>>(),
                vec![10, 90]
            );
            close(step.kinetic_energy_change_j, 0.0);
            close(step.gravity_work_j, -7.5);
            close(
                step.support_work_j,
                if prescribed_floor { 7.5 } else { 0.0 },
            );
            close(step.plane_work_j, if prescribed_floor { 0.0 } else { 7.5 });
            close(step.frictional_dissipation_j, 0.0);
            vector(step.linear_momentum_change_world, Vec3::default());
            vector(step.external_linear_impulse_world, Vec3::default());
            vector(
                step.external_angular_impulse_about_origin_world,
                Vec3::default(),
            );
            close(step.pair_exchanges[0].first_work_j, 5.0);
            close(step.pair_exchanges[1].first_work_j, 2.5);
            if prescribed_floor {
                let ArticulatedSustainedMotion::Prescribed(native) = &step.motion[2] else {
                    panic!()
                };
                vector(
                    native.support_linear_impulse_world,
                    Vec3::new(0.0, 7.5, 0.0),
                );
            }
            audit(&step);
        }
    }
}

#[test]
fn finite_guided_bodies_transfer_static_friction_without_heat_or_pseudo_mass() {
    let models = [cartesian(), cartesian()];
    let states = [
        supported(Vec3::new(2.0, 1.0, 0.0), Vec3::default(), 3),
        supported(Vec3::new(2.0, -1.0, 0.0), Vec3::default(), 3),
    ];
    let forces = [
        [Wrench::default(); 3],
        [
            Wrench::default(),
            Wrench::default(),
            Wrench::new(Vec3::default(), Vec3::new(12.0, 16.0, 0.0)),
        ],
    ];
    let systems: Vec<_> = (0..2)
        .map(|i| System {
            model: &models[i],
            state: &states[i],
            loads: FreeFloatingMotionLoads {
                generalized_force: &[0.0; 3],
                external_body_wrench: &forces[i],
                gravity_world: Vec3::default(),
            },
        })
        .collect();
    let step = step_sustained_articulated_contacts(
        &systems,
        &[pair(70, 0, 1, 2, 1.0)],
        &[[8.0, 6.0, 0.0]],
        0.0,
        controls(),
        support_path,
        manufactured_translation,
        || false,
    )
    .unwrap();
    // Native Cartesian mobilities are diag(1/6,1/4,1/2) for EACH finite tree.
    // The common acceleration is (1,2,0), requiring static traction (6,8,0).
    close(step.solve.forces_n[0][0], 8.0);
    close(step.solve.forces_n[0][1], 6.0);
    close(step.solve.paths_m[0][1], 0.0);
    for after in &step.states_after {
        let State::Prescribed(after) = after else {
            panic!()
        };
        close(after.joint_velocity[0], 0.125);
        close(after.joint_velocity[1], 0.25);
        close(after.joint_velocity[2], 0.0);
    }
    vector(step.linear_momentum_change_world, Vec3::new(1.5, 2.0, 0.0));
    vector(step.external_linear_impulse_world, Vec3::new(1.5, 2.0, 0.0));
    close(step.kinetic_energy_change_j, 0.34375);
    close(step.body_wrench_work_j, 0.34375);
    close(step.support_work_j, 0.0);
    close(step.frictional_dissipation_j, 0.0);
    close(step.pair_exchanges[0].first_work_j, 0.171875);
    close(step.pair_exchanges[0].second_work_j, -0.171875);
    audit(&step);
}

#[test]
fn inactive_finite_pair_separates_without_adhesion() {
    let models = [body(), body()];
    let states = [
        free(Vec3::new(2.0, 1.0, 0.0), Vec3::new(0.0, 0.5, 0.0)),
        free(Vec3::new(2.0, -1.0, 0.0), Vec3::new(0.0, -0.5, 0.0)),
    ];
    let zero = [Wrench::default()];
    let systems: Vec<_> = (0..2)
        .map(|i| System {
            model: &models[i],
            state: &states[i],
            loads: FreeFloatingMotionLoads {
                generalized_force: &[],
                external_body_wrench: &zero,
                gravity_world: Vec3::default(),
            },
        })
        .collect();
    let step = step_sustained_articulated_contacts(
        &systems,
        &[pair(40, 0, 1, 0, 0.5)],
        &[[0.0; 3]],
        0.0,
        controls(),
        |_, _, _| panic!("free bodies have no prescribed path"),
        manufactured_translation,
        || false,
    )
    .unwrap();
    assert_eq!(step.solve.forces_n[0], [0.0; 3]);
    close(step.solve.paths_m[0][0], 0.125);
    close(step.kinetic_energy_change_j, 0.0);
    close(step.pair_exchanges[0].relative_work_j, 0.0);
    vector(step.external_linear_impulse_world, Vec3::default());
    audit(&step);
}

#[test]
fn sliding_material_witnesses_refuse_the_artificial_normal_force_couple() {
    let models = [cartesian(), cartesian()];
    let mut upper = supported(Vec3::new(2.0, 1.0, 0.0), Vec3::default(), 3);
    let State::Prescribed(ref mut state) = upper else {
        panic!()
    };
    state.joint_velocity[0] = 1.0;
    let states = [
        upper,
        supported(Vec3::new(2.0, -1.0, 0.0), Vec3::default(), 3),
    ];
    let saved = states.clone();
    let forces = [
        [Wrench::default(); 3],
        [
            Wrench::default(),
            Wrench::default(),
            Wrench::new(Vec3::default(), Vec3::new(0.0, 16.0, 0.0)),
        ],
    ];
    let systems: Vec<_> = (0..2)
        .map(|i| System {
            model: &models[i],
            state: &states[i],
            loads: FreeFloatingMotionLoads {
                generalized_force: &[0.0; 3],
                external_body_wrench: &forces[i],
                gravity_world: Vec3::default(),
            },
        })
        .collect();
    let error = step_sustained_articulated_contacts(
        &systems,
        &[pair(70, 0, 1, 2, 0.0)],
        &[[8.0, 0.0, 0.0]],
        0.0,
        controls(),
        support_path,
        |_| -> Result<(), &'static str> { panic!("nonreciprocal geometry reached path admission") },
        || false,
    )
    .err()
    .unwrap();
    assert!(matches!(
        error,
        ArticulatedSustainedError::Geometry(ContactGeometryError::NonCollinearContact { .. })
    ));
    assert_eq!(states, saved);
}

#[test]
fn tiny_normal_clearance_cannot_hide_a_large_relative_tangent_couple() {
    let models = [cartesian(), cartesian()];
    let gap = 2.0_f64.powi(-40);
    let states = [
        supported(Vec3::new(2.0, 1.0 + gap, 0.0), Vec3::default(), 3),
        supported(Vec3::new(2.0, -1.0, 0.0), Vec3::default(), 3),
    ];
    let forces = [
        [Wrench::default(); 3],
        [
            Wrench::default(),
            Wrench::default(),
            Wrench::new(Vec3::default(), Vec3::new(12.0, 16.0, 0.0)),
        ],
    ];
    let systems: Vec<_> = (0..2)
        .map(|i| System {
            model: &models[i],
            state: &states[i],
            loads: FreeFloatingMotionLoads {
                generalized_force: &[0.0; 3],
                external_body_wrench: &forces[i],
                gravity_world: Vec3::default(),
            },
        })
        .collect();
    let mut c = controls();
    // This residual is resolved by the declared normal-force allowance, but
    // the same gap times tangential force exceeds the independently tighter
    // angular-impulse allowance. A law-only check would miss that distinction.
    c.solve.newton.absolute_tolerance = 1e-8;
    c.angular_momentum_tolerance_nms = 1e-14;
    let error = step_sustained_articulated_contacts(
        &systems,
        &[pair(70, 0, 1, 2, 1.0)],
        &[[8.0, 6.0, 0.0]],
        0.0,
        c,
        support_path,
        |_| -> Result<(), &'static str> { panic!("artificial couple reached path admission") },
        || false,
    )
    .err()
    .unwrap();
    let ArticulatedSustainedError::Reciprocity {
        contact_id,
        angular_impulse_world,
    } = error
    else {
        panic!("expected reciprocal-moment refusal, got {error:?}");
    };
    assert_eq!(contact_id, 70);
    assert!(angular_impulse_world.z < -c.angular_momentum_tolerance_nms);
    close(angular_impulse_world.z / gap, -0.75);
}
