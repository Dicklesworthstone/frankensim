//! Independent reductions and momentum transfer through rigid closed linkages.
use fs_contact::frictional_multicontact::{FrictionalContactConfig, NormalVelocityLaw};
use fs_ga::{Mat3, Se3, So3, So3Tangent, Twist, Vec3, Wrench};
use fs_mbd::articulated::motion::spatial_loads::implicit::distance::{
    DistanceLoop, step_free_floating_distance_loops,
};
use fs_mbd::articulated::motion::spatial_loads::{MaterialPoint, SpatialMotionLoads};
use fs_mbd::articulated::motion::{FreeFloatingMotionLoads, FreeFloatingMotionState};
use fs_mbd::articulated::{
    ArticulatedModel, FreeFloatingBaseState, JointModel, Link, SpatialInertia,
};
use fs_mbd::articulated_contacts::friction::constrained::{
    ConstrainedArticulatedContactResult, ConstrainedImpactConfig, ImpactDistanceConstraint,
    MAX_IMPACT_DISTANCE_CONSTRAINTS, resolve_constrained_articulated_frictional_contacts,
};
use fs_mbd::articulated_contacts::friction::{
    ArticulatedFrictionalContact, ArticulatedFrictionalContactResult,
    resolve_articulated_frictional_contacts,
};
use fs_mbd::articulated_contacts::{
    ArticulatedContactSystem, ContactCounterpart, LinkContactPoint,
};

fn close(a: f64, b: f64) {
    assert!((a - b).abs() < 3e-10, "{a:.16e} != {b:.16e}");
}
fn vector(a: Vec3, b: Vec3) {
    assert!((a - b).norm() < 3e-10, "{a:?} != {b:?}");
}
fn pose(position: Vec3) -> Se3 {
    Se3::from_parts(So3::identity(), position).unwrap()
}
fn inertia(mass: f64, center: Vec3) -> SpatialInertia {
    SpatialInertia::new(
        mass,
        center,
        Mat3 {
            m: [1.2, 0.11, 0.04, 0.11, 1.4, 0.08, 0.04, 0.08, 1.6],
        },
    )
    .unwrap()
}
fn body(mass: f64) -> ArticulatedModel {
    ArticulatedModel::new(vec![Link::new(
        "body",
        None,
        Se3::identity(),
        JointModel::FIXED,
        inertia(mass, Vec3::default()),
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
        base: FreeFloatingBaseState::new(pose(position), Twist::new(Vec3::default(), velocity)),
        joint_position: &[],
        joint_velocity: &[],
    }
}
fn point(system_index: usize, link_index: usize) -> LinkContactPoint {
    LinkContactPoint {
        system_index,
        link_index,
        point_body: Vec3::default(),
    }
}
fn rod(a: LinkContactPoint, b: LinkContactPoint, length_m: f64) -> ImpactDistanceConstraint {
    ImpactDistanceConstraint { a, b, length_m }
}
fn contact(
    id: u64,
    point: LinkContactPoint,
    world: Vec3,
    normal: Vec3,
    tangent: Vec3,
    mu: f64,
    e: f64,
) -> ArticulatedFrictionalContact {
    ArticulatedFrictionalContact {
        id,
        point,
        other: ContactCounterpart::PrescribedPlane {
            point_world: world,
            velocity_world: Vec3::default(),
        },
        normal_world: normal,
        tangent_world: tangent,
        activation_tolerance: 1e-8,
        friction_coefficient: mu,
        normal_law: NormalVelocityLaw::NewtonRestitution { coefficient: e },
    }
}
fn solve(
    states: &[ArticulatedContactSystem<'_>],
    contacts: &[ArticulatedFrictionalContact],
    rods: &[ImpactDistanceConstraint],
) -> ConstrainedArticulatedContactResult {
    resolve_constrained_articulated_frictional_contacts(
        states,
        contacts,
        rods,
        FrictionalContactConfig::default(),
        ConstrainedImpactConfig::default(),
        || false,
    )
    .unwrap()
}
fn audit(result: &ConstrainedArticulatedContactResult) {
    let r = &result.impact;
    close(r.energy_balance_residual, 0.);
    close(r.endpoint_velocity_residual, 0.);
    vector(r.linear_momentum_residual, Vec3::default());
    vector(r.angular_momentum_residual, Vec3::default());
    assert!(result.velocity_after_mps.iter().all(|v| v.abs() < 1e-9));
    assert!(result.length_residual_m.iter().all(|v| v.abs() < 1e-9));
    let change: f64 = r
        .responses
        .iter()
        .map(|p| p.momentum_after.kinetic_energy - p.momentum_before.kinetic_energy)
        .sum();
    close(
        change - r.plane_work + r.contact.dissipated_energy,
        r.energy_balance_residual,
    );
    let work: f64 = result
        .distance_impulses_ns
        .iter()
        .zip(&result.velocity_before_mps)
        .zip(&result.velocity_after_mps)
        .map(|((j, a), b)| 0.5 * j * (a + b))
        .sum();
    close(result.constraint_work_j, work);
}

#[test]
fn cross_tree_rod_transmits_impact_and_keeps_independent_tangential_mobility() {
    let models = [body(2.), body(3.)];
    let states = [
        system(&models[0], Vec3::default(), Vec3::new(-1., 2., 0.)),
        system(&models[1], Vec3::new(1., 0., 0.), Vec3::new(-1., 2., 0.)),
    ];
    let rods = [rod(point(0, 0), point(1, 0), 1.)];
    for e in [0., 0.5, 1.] {
        let c = contact(
            7,
            point(0, 0),
            Vec3::default(),
            Vec3::new(1., 0., 0.),
            Vec3::new(0., 1., 0.),
            0.2,
            e,
        );
        let r = solve(&states, &[c], &rods);
        audit(&r);
        let jn = 5. * (1. + e);
        let jt = -0.2 * jn;
        close(r.impact.delassus[0], 0.2);
        close(r.impact.delassus[4], 0.5);
        close(r.impact.contact.impulses[0][0], jn);
        close(r.impact.contact.impulses[0][1], jt);
        close(r.distance_impulses_ns[0], 3. * (1. + e));
        vector(
            r.impact.responses[0].base_after.twist_body.linear,
            Vec3::new(e, 2. + jt / 2., 0.),
        );
        vector(
            r.impact.responses[1].base_after.twist_body.linear,
            Vec3::new(e, 2., 0.),
        );
        close(
            r.impact.contact.dissipated_energy,
            2.5 * (1. - e * e) - 2. * jt - jt * jt / 4.,
        );
        close(r.constraint_work_j, 0.);
    }
}

fn slider(fixed: bool) -> ArticulatedModel {
    ArticulatedModel::new(vec![
        Link::new(
            "carrier",
            None,
            Se3::identity(),
            JointModel::FIXED,
            inertia(3., Vec3::new(0.1, -0.04, 0.08)),
        ),
        Link::new(
            "slider",
            Some(0),
            pose(Vec3::new(if fixed { 1.2 } else { 1. }, 0., 0.)),
            if fixed {
                JointModel::FIXED
            } else {
                JointModel::prismatic(Vec3::new(1., 0., 0.), None).unwrap()
            },
            inertia(2., Vec3::new(-0.1, 0.07, -0.02)),
        ),
    ])
    .unwrap()
}
fn compare(a: &ArticulatedFrictionalContactResult, b: &ArticulatedFrictionalContactResult) {
    for (x, y) in a.delassus.iter().zip(&b.delassus) {
        close(*x, *y);
    }
    for (x, y) in a
        .contact
        .impulses
        .iter()
        .flatten()
        .zip(b.contact.impulses.iter().flatten())
    {
        close(*x, *y);
    }
    for (x, y) in a.responses.iter().zip(&b.responses) {
        vector(
            x.base_after.twist_body.linear,
            y.base_after.twist_body.linear,
        );
        vector(
            x.base_after.twist_body.angular,
            y.base_after.twist_body.angular,
        );
    }
    close(a.contact.dissipated_energy, b.contact.dissipated_energy);
}
#[test]
fn locked_articulated_slider_matches_independent_fixed_joint_reduction_with_full_inertia() {
    let model = slider(false);
    let reduced = slider(true);
    let rods = [rod(point(0, 0), point(0, 1), 1.2)];
    for spin in [Vec3::default(), Vec3::new(0.3, -0.4, 0.2)] {
        let rotation = So3::exp(So3Tangent::new(spin)).unwrap();
        let scene = Se3::from_parts(rotation, Vec3::new(1., -2., 0.5)).unwrap();
        let base = FreeFloatingBaseState::new(
            scene,
            Twist::new(Vec3::new(0.2, -0.15, 0.1), Vec3::new(-2., 1.4, 0.2)),
        );
        let original = ArticulatedContactSystem {
            model: &model,
            base,
            joint_position: &[0.2],
            joint_velocity: &[0.],
        };
        let independent = ArticulatedContactSystem {
            model: &reduced,
            base,
            joint_position: &[],
            joint_velocity: &[],
        };
        let endpoint = LinkContactPoint {
            point_body: Vec3::new(0.13, -0.2, 0.17),
            ..point(0, 1)
        };
        let world = scene
            .transform_point(Vec3::new(1.2, 0., 0.) + endpoint.point_body)
            .unwrap();
        for mu in [0., 0.4] {
            let c = contact(
                1,
                endpoint,
                // A small positive witness gap stays within activation while
                // avoiding roundoff-negative penetration under scene rotation.
                world - rotation.rotate(Vec3::new(1e-10, 0., 0.)).unwrap(),
                rotation.rotate(Vec3::new(1., 0., 0.)).unwrap(),
                rotation.rotate(Vec3::new(0., 1., 0.)).unwrap(),
                mu,
                0.35,
            );
            let actual = solve(&[original], &[c], &rods);
            audit(&actual);
            let expected = resolve_articulated_frictional_contacts(
                &[independent],
                &[c],
                FrictionalContactConfig::default(),
                || false,
            )
            .unwrap();
            compare(&actual.impact, &expected);
            close(actual.impact.responses[0].joint_velocity_after[0], 0.);
        }
    }
}

#[test]
fn simultaneous_shared_rods_and_redundant_contacts_keep_unique_velocity_and_order_parity() {
    let models = [body(1.), body(2.), body(3.)];
    let states: Vec<_> = (0..3)
        .map(|i| {
            system(
                &models[i],
                Vec3::new(i as f64, 0., 0.),
                Vec3::new(-1., 0., 0.),
            )
        })
        .collect();
    let rods = [
        rod(point(0, 0), point(1, 0), 1.),
        rod(point(1, 0), point(2, 0), 1.),
    ];
    let contacts = [
        contact(
            11,
            point(0, 0),
            Vec3::default(),
            Vec3::new(1., 0., 0.),
            Vec3::new(0., 1., 0.),
            0.,
            0.,
        ),
        contact(
            5,
            point(2, 0),
            Vec3::new(2., 0., 0.),
            Vec3::new(1., 0., 0.),
            Vec3::new(0., 1., 0.),
            0.,
            0.,
        ),
    ];
    let reference = solve(&states, &contacts, &rods);
    audit(&reference);
    close(
        reference.impact.contact.impulses.iter().map(|p| p[0]).sum(),
        6.,
    );
    close(reference.impact.contact.dissipated_energy, 3.);
    for r in &reference.impact.responses {
        vector(r.base_after.twist_body.linear, Vec3::default());
    }
    let permuted = solve(&states, &[contacts[1], contacts[0]], &[rods[1], rods[0]]);
    audit(&permuted);
    compare(&reference.impact, &permuted.impact);
    close(
        reference.distance_impulses_ns[0],
        permuted.distance_impulses_ns[1],
    );
    close(
        reference.distance_impulses_ns[1],
        permuted.distance_impulses_ns[0],
    );
}

#[test]
fn constrained_contact_continues_the_original_rigid_loop_motion() {
    let model = slider(false);
    let mut state = FreeFloatingMotionState {
        base: FreeFloatingBaseState::new(
            Se3::identity(),
            Twist::new(Vec3::default(), Vec3::new(-1., 0., 0.)),
        ),
        joint_position: vec![0.2],
        joint_velocity: vec![0.],
    };
    let links = [
        MaterialPoint {
            link: 0,
            point_body: Vec3::default(),
        },
        MaterialPoint {
            link: 1,
            point_body: Vec3::default(),
        },
    ];
    let loops = [DistanceLoop {
        a: links[0],
        b: links[1],
        length_m: 1.2,
    }];
    let zero = Wrench::new(Vec3::default(), Vec3::default());
    let loads = SpatialMotionLoads {
        held: FreeFloatingMotionLoads {
            generalized_force: &[0.],
            external_body_wrench: &[zero; 2],
            gravity_world: Vec3::default(),
        },
        world_points: &[],
        connections: &[],
    };
    state = step_free_floating_distance_loops(
        &model,
        &state,
        loads,
        &loops,
        0.1,
        Default::default(),
        || false,
    )
    .unwrap()
    .state_after;
    let supplied = ArticulatedContactSystem {
        model: &model,
        base: state.base,
        joint_position: &state.joint_position,
        joint_velocity: &state.joint_velocity,
    };
    let endpoint = LinkContactPoint {
        point_body: Vec3::new(0.1, -0.2, 0.15),
        ..point(0, 1)
    };
    let world = state
        .base
        .world_from_base
        .transform_point(Vec3::new(1.2, 0., 0.) + endpoint.point_body)
        .unwrap();
    let c = contact(
        1,
        endpoint,
        world,
        Vec3::new(1., 0., 0.),
        Vec3::new(0., 1., 0.),
        0.2,
        0.4,
    );
    let impact = solve(&[supplied], &[c], &[rod(point(0, 0), point(0, 1), 1.2)]);
    audit(&impact);
    state.base = impact.impact.responses[0].base_after;
    state.joint_velocity = impact.impact.responses[0].joint_velocity_after.clone();
    let next = step_free_floating_distance_loops(
        &model,
        &state,
        loads,
        &loops,
        0.02,
        Default::default(),
        || false,
    )
    .unwrap();
    close(next.state_after.joint_position[0], 0.2);
    close(next.state_after.joint_velocity[0], 0.);
    assert!(next.position_residual_after_m[0].abs() < 1e-9);
    assert!(next.velocity_residual_after_mps[0].abs() < 1e-9);
    assert!(next.linear_momentum_residual_world.norm() < 1e-9);
    assert!(next.angular_momentum_residual_world.norm() < 1e-9);
}

#[test]
fn empty_constraint_path_is_the_original_contact_owner() {
    let model = body(2.);
    let state = system(&model, Vec3::default(), Vec3::new(-1., 0.5, 0.));
    let c = contact(
        2,
        point(0, 0),
        Vec3::default(),
        Vec3::new(1., 0., 0.),
        Vec3::new(0., 1., 0.),
        0.3,
        0.2,
    );
    let plain = resolve_articulated_frictional_contacts(
        &[state],
        &[c],
        FrictionalContactConfig::default(),
        || false,
    )
    .unwrap();
    let delegated = solve(&[state], &[c], &[]);
    assert_eq!(delegated.impact.delassus, plain.delassus);
    assert_eq!(delegated.impact.contact, plain.contact);
    assert_eq!(
        delegated.impact.responses[0].base_after,
        plain.responses[0].base_after
    );
    assert!(delegated.distance_impulses_ns.is_empty());
    close(delegated.constraint_work_j, 0.);
}

#[test]
fn inconsistent_singular_and_late_cancelled_events_never_publish_a_partial_state() {
    let models = [body(2.), body(3.)];
    let states = [
        system(&models[0], Vec3::default(), Vec3::new(-1., 0.5, 0.)),
        system(&models[1], Vec3::new(1., 0., 0.), Vec3::new(-1., 0.5, 0.)),
    ];
    let original = [states[0].base, states[1].base];
    let c = contact(
        1,
        point(0, 0),
        Vec3::default(),
        Vec3::new(1., 0., 0.),
        Vec3::new(0., 1., 0.),
        0.2,
        0.2,
    );
    let d = rod(point(0, 0), point(1, 0), 1.);
    let refuse = |states: &[ArticulatedContactSystem<'_>], rods: &[ImpactDistanceConstraint]| {
        resolve_constrained_articulated_frictional_contacts(
            states,
            &[c],
            rods,
            FrictionalContactConfig::default(),
            ConstrainedImpactConfig::default(),
            || false,
        )
        .is_err()
    };
    assert!(refuse(&states, &[d, d]));
    assert!(refuse(
        &states,
        &[ImpactDistanceConstraint { length_m: 1.1, ..d }]
    ));
    assert!(refuse(
        &states,
        &vec![d; MAX_IMPACT_DISTANCE_CONSTRAINTS + 1]
    ));
    let mut inconsistent = states;
    inconsistent[1].base.twist_body.linear.x = 0.;
    assert!(refuse(&inconsistent, &[d]));
    let mut count = 0;
    let expected = resolve_constrained_articulated_frictional_contacts(
        &states,
        &[c],
        &[d],
        FrictionalContactConfig::default(),
        ConstrainedImpactConfig::default(),
        || {
            count += 1;
            false
        },
    )
    .unwrap();
    assert!(count > 10);
    for stop in [1, count / 2, count] {
        let mut polls = 0;
        assert!(
            resolve_constrained_articulated_frictional_contacts(
                &states,
                &[c],
                &[d],
                FrictionalContactConfig::default(),
                ConstrainedImpactConfig::default(),
                || {
                    polls += 1;
                    polls == stop
                }
            )
            .is_err()
        );
        assert_eq!([states[0].base, states[1].base], original);
        let retry = solve(&states, &[c], &[d]);
        assert_eq!(retry.impact.contact, expected.impact.contact);
        assert_eq!(retry.distance_impulses_ns, expected.distance_impulses_ns);
        for (a, b) in retry
            .impact
            .responses
            .iter()
            .zip(&expected.impact.responses)
        {
            assert_eq!(a.base_after, b.base_after);
        }
    }
}
