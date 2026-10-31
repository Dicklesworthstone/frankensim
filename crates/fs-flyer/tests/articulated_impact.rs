//! Full dry/fluid impulse physics through the original global normal law.
use fs_contact::multicontact::MultiImpactConfig;
use fs_flyer::addedmass::Strip;
use fs_flyer::addedmass::articulated::{LinkedStrip, assemble_articulated_strips};
use fs_flyer::partitioned::MemoryState;
use fs_flyer::partitioned::articulated::ArticulatedState;
use fs_flyer::partitioned::articulated::impact::{
    ArticulatedNormalContact, ContactCounterpart, LinkContactPoint, resolve,
};
use fs_ga::{Mat3, Se3, So3, Twist, Vec3, Wrench};
use fs_mbd::articulated::motion::FreeFloatingMotionState;
use fs_mbd::articulated::{
    ArticulatedModel, BaseState, FreeFloatingBaseState, JointLimits, JointModel, Link,
    SpatialInertia, apply_free_floating_impulses, forward_kinematics, free_floating_momentum,
};
use std::cell::Cell;

fn v(x: f64, y: f64, z: f64) -> Vec3 {
    Vec3::new(x, y, z)
}
fn close(a: f64, b: f64) {
    assert!((a - b).abs() < 3e-10 * (1.0 + b.abs()), "{a:e} != {b:e}");
}
fn vector_close(a: Vec3, b: Vec3) {
    assert!((a - b).norm() < 3e-10 * (1.0 + b.norm()), "{a:?} != {b:?}");
}
fn inertia(m: f64, c: Vec3, d: [f64; 3]) -> SpatialInertia {
    SpatialInertia::new(
        m,
        c,
        Mat3 {
            m: [d[0], 0.0, 0.0, 0.0, d[1], 0.0, 0.0, 0.0, d[2]],
        },
    )
    .unwrap()
}
fn translation(x: f64) -> Se3 {
    Se3::from_parts(So3::identity(), v(x, 0.0, 0.0)).unwrap()
}
fn point(link: usize, p: Vec3) -> LinkContactPoint {
    LinkContactPoint {
        system_index: 0,
        link_index: link,
        point_body: p,
    }
}
fn plane(
    point: LinkContactPoint,
    position: Vec3,
    normal: Vec3,
    speed: Vec3,
    restitution: f64,
) -> ArticulatedNormalContact {
    ArticulatedNormalContact {
        point,
        other: ContactCounterpart::PrescribedPlane {
            point_world: position,
            velocity_world: speed,
        },
        normal_world: normal,
        activation_tolerance: 0.0,
        restitution,
    }
}
fn fixture() -> (
    ArticulatedModel,
    ArticulatedState,
    Vec<ArticulatedNormalContact>,
    Strip,
) {
    let limits = JointLimits::new(-2.0, 2.0, 100.0, 1e-8).unwrap();
    let model = ArticulatedModel::new(vec![
        Link::new(
            "carrier",
            None,
            Se3::identity(),
            JointModel::FIXED,
            inertia(4.0, v(0.1, -0.05, 0.0), [2.0, 2.5, 3.0]),
        ),
        Link::new(
            "hinge",
            Some(0),
            translation(0.5),
            JointModel::revolute(v(0.0, 1.0, 0.0), Some(limits)).unwrap(),
            inertia(1.8, v(0.2, 0.0, 0.05), [0.2, 0.3, 0.4]),
        ),
        Link::new(
            "second hinge",
            Some(1),
            translation(0.5),
            JointModel::revolute(v(0.0, 0.0, 1.0), Some(limits)).unwrap(),
            inertia(1.2, v(0.1, -0.1, 0.1), [0.15, 0.2, 0.25]),
        ),
    ])
    .unwrap();
    let state = ArticulatedState {
        motion: FreeFloatingMotionState {
            base: FreeFloatingBaseState::new(
                Se3::identity(),
                Twist::new(v(0.1, -0.1, 0.2), v(-1.0, 0.3, -0.2)),
            ),
            joint_position: vec![0.0; 2],
            joint_velocity: vec![1.0, 0.4],
        },
        memory: vec![MemoryState {
            tau_s: 0.1,
            value: 0.7,
        }],
    };
    let contacts = vec![
        ArticulatedNormalContact {
            point: point(2, v(0.5, 0.2, 0.0)),
            other: ContactCounterpart::Body(point(0, v(1.5, 0.2, 0.0))),
            normal_world: v(0.0, 0.0, 1.0),
            activation_tolerance: 0.0,
            restitution: 0.0,
        },
        plane(
            point(1, v(0.0, 0.1, 0.2)),
            v(0.5, 0.1, 0.2),
            v(1.0, 0.0, 0.0),
            v(0.2, 0.0, 0.0),
            0.0,
        ),
    ];
    let strip = Strip {
        name: "moving fluid plate",
        chord_m: 1.2,
        span_m: 1.6,
        position_m: [0.3, -0.1, 0.1],
        normal: [0.0, 0.0, 1.0],
        span_axis: [0.0, 1.0, 0.0],
        control_coord: None,
        control_gain: 0.0,
        control_angular_gain: [0.0; 3],
    };
    (model, state, contacts, strip)
}

fn total(
    model: &ArticulatedModel,
    state: &ArticulatedState,
    strip: &Strip,
    rho: f64,
) -> (f64, Vec3, Vec3) {
    let motion = &state.motion;
    let dry = free_floating_momentum(
        model,
        motion.base,
        &motion.joint_position,
        &motion.joint_velocity,
    )
    .unwrap();
    let kin = forward_kinematics(
        model,
        BaseState::prescribed(
            motion.base.world_from_base,
            motion.base.twist_body,
            Twist::zero(),
        ),
        &motion.joint_position,
        &motion.joint_velocity,
    )
    .unwrap();
    // Independent physical local strip energy and momentum, not the prepared
    // generalized matrix, Delassus solve or returned work accounting.
    let speed = kin.body_twist[2];
    let r = v(
        strip.position_m[0],
        strip.position_m[1],
        strip.position_m[2],
    );
    let normal = v(0.0, 0.0, 1.0);
    let span = v(0.0, 1.0, 0.0);
    let chord = span.cross(normal);
    let m = rho * std::f64::consts::PI * strip.chord_m.powi(2) * strip.span_m / 4.0;
    let pitch = m * strip.chord_m.powi(2) / 32.0;
    let rotary = m * strip.span_m.powi(2) / 12.0;
    let vn = normal.dot(speed.linear + speed.angular.cross(r));
    let wp = span.dot(speed.angular);
    let wr = chord.dot(speed.angular);
    let force = normal.scale(m * vn);
    let local = Wrench::new(
        r.cross(force) + span.scale(pitch * wp) + chord.scale(rotary * wr),
        force,
    );
    let world = local.transform_by(&kin.world_from_link[2]);
    (
        dry.kinetic_energy + 0.5 * (m * vn * vn + pitch * wp * wp + rotary * wr * wr),
        dry.linear_momentum_world + world.force,
        dry.angular_momentum_world + world.torque,
    )
}

#[test]
fn simultaneous_self_and_moving_plane_impacts_close_full_fluid_energy_and_momentum() {
    let (model, state, contacts, strip) = fixture();
    let saved = state.clone();
    let rho = 1.1;
    let attached = [LinkedStrip::new(2, strip.clone()).unwrap()];
    let before = total(&model, &state, &strip, rho);
    let event = resolve(
        &model,
        &state,
        &contacts,
        MultiImpactConfig::default(),
        |map| assemble_articulated_strips(rho, &attached, map),
        || false,
    )
    .unwrap();
    let after = total(&model, &event.state, &strip, rho);
    let response = &event.generalized;
    let impulse = &response.law.normal_impulses;
    assert!(impulse.iter().all(|x| *x > 0.0), "{impulse:?}");
    assert!(response.delassus[1].abs() > 1e-4);
    let external = v(impulse[1], 0.0, 0.0);
    vector_close(after.1 - before.1, external);
    vector_close(after.2 - before.2, v(0.5, 0.1, 0.2).cross(external));
    close(response.kinetic_before, before.0);
    close(response.kinetic_after, after.0);
    close(response.kinetic_change, after.0 - before.0);
    close(response.support_work, 0.2 * impulse[1]);
    close(
        after.0 - before.0 - response.support_work + response.law.dissipated_energy,
        0.0,
    );
    for speed in &event.relative_normal_speeds_after {
        close(*speed, 0.0);
    }
    assert!(event.linear_momentum_residual.norm() < 1e-10);
    assert!(event.angular_momentum_residual.norm() < 1e-10);
    assert_eq!(
        event.state.motion.base.world_from_base,
        state.motion.base.world_from_base
    );
    assert_eq!(
        event.state.motion.joint_position,
        state.motion.joint_position
    );
    assert_eq!(event.state.memory, state.memory);
    assert_eq!(state, saved);
    let dry = resolve(
        &model,
        &state,
        &contacts,
        MultiImpactConfig::default(),
        |map| assemble_articulated_strips(rho, &[], map),
        || false,
    )
    .unwrap();
    assert!(
        (dry.generalized.law.normal_impulses[0] - impulse[0]).abs() > 1e-3,
        "fluid inertia must change the simultaneous impulse"
    );
}

#[test]
fn zero_fluid_reset_matches_original_native_aba_impulse_application() {
    let (model, state, contacts, _) = fixture();
    let event = resolve(
        &model,
        &state,
        &contacts,
        MultiImpactConfig::default(),
        |map| assemble_articulated_strips(1.1, &[], map),
        || false,
    )
    .unwrap();
    let impulse = &event.generalized.law.normal_impulses;
    // At this zero-angle event all original material frames have world axes.
    // Build action/reaction at the declared points independently of row prep.
    let mut loads = [Wrench::default(); 3];
    let first = v(0.0, 0.0, impulse[0]);
    let second = v(impulse[1], 0.0, 0.0);
    loads[2] = Wrench::new(v(0.5, 0.2, 0.0).cross(first), first);
    loads[0] = Wrench::new(v(1.5, 0.2, 0.0).cross(first.scale(-1.0)), first.scale(-1.0));
    loads[1] = Wrench::new(v(0.0, 0.1, 0.2).cross(second), second);
    let native = apply_free_floating_impulses(
        &model,
        state.motion.base,
        &state.motion.joint_position,
        &state.motion.joint_velocity,
        &[0.0; 2],
        &loads,
    )
    .unwrap();
    for (actual, expected) in event
        .state
        .motion
        .base
        .twist_body
        .to_array()
        .into_iter()
        .zip(native.base_after.twist_body.to_array())
    {
        close(actual, expected);
    }
    for (actual, expected) in event
        .state
        .motion
        .joint_velocity
        .iter()
        .zip(native.joint_velocity_after)
    {
        close(*actual, expected);
    }
    close(
        event.generalized.kinetic_after,
        native.momentum_after.kinetic_energy,
    );
    // The tiny native effort limits do not constrain an impulse, which has
    // different units. Native final velocity limits still apply.
    assert!(impulse.iter().any(|x| *x > 1e-8));
}

fn bar() -> (
    ArticulatedModel,
    ArticulatedState,
    Vec<ArticulatedNormalContact>,
) {
    let model = ArticulatedModel::new(vec![Link::new(
        "bar",
        None,
        Se3::identity(),
        JointModel::FIXED,
        inertia(2.0, Vec3::default(), [0.5; 3]),
    )])
    .unwrap();
    let state = ArticulatedState {
        motion: FreeFloatingMotionState {
            base: FreeFloatingBaseState::new(
                Se3::identity(),
                Twist::new(v(0.0, -0.75, 0.0), v(0.0, 0.0, -0.25)),
            ),
            joint_position: vec![],
            joint_velocity: vec![],
        },
        memory: vec![],
    };
    let contacts = [-1.0, 1.0]
        .map(|x| {
            plane(
                point(0, v(x, 0.0, 0.0)),
                v(x, 0.0, 0.0),
                v(0.0, 0.0, 1.0),
                Vec3::default(),
                1.0,
            )
        })
        .to_vec();
    (model, state, contacts)
}

#[test]
fn incompatible_simultaneous_rebound_refuses_instead_of_injecting_global_energy() {
    let (model, state, mut contacts) = bar();
    let saved = state.clone();
    // Independent W=[[2.5,-1.5],[-1.5,2.5]], initial speeds[-1,.5].
    // Elastic pairwise targets would inject 0.109375J globally.
    let error = resolve(
        &model,
        &state,
        &contacts,
        MultiImpactConfig::default(),
        |map| assemble_articulated_strips(1.0, &[], map),
        || false,
    )
    .unwrap_err();
    assert_eq!(error.code, "articulated-impact-law");
    assert!(error.message.contains("EnergyIncrease"));
    assert_eq!(state, saved);
    for contact in &mut contacts {
        contact.restitution = 0.0;
    }
    let out = resolve(
        &model,
        &state,
        &contacts,
        MultiImpactConfig::default(),
        |map| assemble_articulated_strips(1.0, &[], map),
        || false,
    )
    .unwrap();
    for (got, want) in out.generalized.delassus.iter().zip([2.5, -1.5, -1.5, 2.5]) {
        close(*got, want);
    }
    close(out.generalized.law.normal_impulses[0], 0.4375);
    close(out.generalized.law.normal_impulses[1], 0.0625);
    close(out.generalized.law.dissipated_energy, 0.203125);
    vector_close(out.state.motion.base.twist_body.linear, Vec3::default());
    vector_close(out.state.motion.base.twist_body.angular, Vec3::default());
}

#[test]
fn geometry_load_ownership_caps_and_cancellation_refuse_atomically() {
    let (model, state, mut contacts) = bar();
    let saved = state.clone();
    for contact in &mut contacts {
        contact.restitution = 0.0;
    }
    let cancel = Cell::new(false);
    let error = resolve(
        &model,
        &state,
        &contacts,
        MultiImpactConfig::default(),
        |map| {
            cancel.set(true);
            assemble_articulated_strips(1.0, &[], map)
        },
        || cancel.get(),
    )
    .unwrap_err();
    assert_eq!(error.code, "articulated-impact-cancelled");
    let error = resolve(
        &model,
        &state,
        &contacts,
        MultiImpactConfig::default(),
        |map| {
            let mut fluid = assemble_articulated_strips(1.0, &[], map)?;
            fluid.q_rigid_nonaccel[0] = 1.0;
            Ok(fluid)
        },
        || false,
    )
    .unwrap_err();
    assert_eq!(error.code, "coupled-load-ownership");
    let mut penetrating = contacts[0];
    penetrating.other = ContactCounterpart::PrescribedPlane {
        point_world: v(-1.0, 0.0, 1e-12),
        velocity_world: Vec3::default(),
    };
    penetrating.activation_tolerance = 1e-3;
    assert!(
        resolve(
            &model,
            &state,
            &[penetrating],
            MultiImpactConfig::default(),
            |map| assemble_articulated_strips(1.0, &[], map),
            || false
        )
        .is_err()
    );
    assert!(
        resolve(
            &model,
            &state,
            &[contacts[0]; 2],
            MultiImpactConfig::default(),
            |map| assemble_articulated_strips(1.0, &[], map),
            || false
        )
        .is_err()
    );
    let mut inactive = contacts[0];
    inactive.other = ContactCounterpart::PrescribedPlane {
        point_world: v(-1.0, 0.0, -0.1),
        velocity_world: Vec3::default(),
    };
    let out = resolve(
        &model,
        &state,
        &[inactive; 64],
        MultiImpactConfig::default(),
        |map| assemble_articulated_strips(1.0, &[], map),
        || false,
    )
    .unwrap();
    assert_eq!(out.state, state);
    assert!(out.generalized.law.active_indices.is_empty());
    assert!(
        resolve(
            &model,
            &state,
            &[inactive; 65],
            MultiImpactConfig::default(),
            |map| assemble_articulated_strips(1.0, &[], map),
            || false
        )
        .is_err()
    );
    let (tree, tree_state, mut self_contacts, _) = fixture();
    self_contacts[0].other = ContactCounterpart::Body(point(0, v(1.6, 0.2, 0.0)));
    let error = resolve(
        &tree,
        &tree_state,
        &self_contacts,
        MultiImpactConfig::default(),
        |map| assemble_articulated_strips(1.0, &[], map),
        || false,
    )
    .unwrap_err();
    assert_eq!(error.code, "articulated-impact-invalid");
    assert!(error.message.contains("noncollinear"));
    assert_eq!(state, saved);
}
