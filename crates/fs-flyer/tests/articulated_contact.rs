//! Sustained dry/fluid contact: independent force, motion and work oracles.
use fs_contact::frictional_multicontact::{
    FrictionalContactConfig, FrictionalContactDeterminacy, NormalVelocityLaw,
};
use fs_flyer::Refusal;
use fs_flyer::addedmass::{
    Strip,
    articulated::{LinkedStrip, assemble_articulated_strips},
};
use fs_flyer::partitioned::MemoryState;
use fs_flyer::partitioned::articulated::contact::{
    ArticulatedFrictionalContact, ContactPath, ContactStepConfig, step_with_contacts,
};
use fs_flyer::partitioned::articulated::{
    ArticulatedLoads, ArticulatedState, ArticulatedStepConfig, JointImpedance,
};
use fs_ga::{Mat3, Se3, So3, So3Tangent, Twist, Vec3, Wrench};
use fs_mbd::articulated::contact_geometry::{ContactCounterpart, LinkContactPoint};
use fs_mbd::articulated::generalized::GeneralizedKinematics;
use fs_mbd::articulated::motion::FreeFloatingMotionState;
use fs_mbd::articulated::{
    ArticulatedModel, FreeFloatingBaseState, JointModel, Link, SpatialInertia,
};
use std::cell::Cell;

const RHO: f64 = 1.2;
const MASS: f64 = 2.0;
const G: f64 = 9.8;
const GAP: f64 = 1e-6;
fn v(x: f64, y: f64, z: f64) -> Vec3 {
    Vec3::new(x, y, z)
}
fn close(a: f64, b: f64) {
    assert!((a - b).abs() < 3e-9 * (1.0 + b.abs()), "{a:e} != {b:e}");
}
fn vec_close(a: Vec3, b: Vec3) {
    close((a - b).norm(), 0.0);
}
fn refused(message: &str) -> Refusal {
    Refusal {
        code: "path-falsifier",
        message: message.into(),
        ranked_repairs: vec![],
    }
}
fn inertia(m: f64, j: f64) -> SpatialInertia {
    SpatialInertia::new(
        m,
        Vec3::default(),
        Mat3 {
            m: [j, 0., 0., 0., j, 0., 0., 0., j],
        },
    )
    .unwrap()
}
fn body() -> ArticulatedModel {
    ArticulatedModel::new(vec![Link::new(
        "body",
        None,
        Se3::identity(),
        JointModel::FIXED,
        inertia(MASS, 5.),
    )])
    .unwrap()
}
fn initial(pose: Se3, linear: Vec3) -> ArticulatedState {
    ArticulatedState {
        motion: FreeFloatingMotionState {
            base: FreeFloatingBaseState::new(pose, Twist::new(Vec3::default(), linear)),
            joint_position: vec![],
            joint_velocity: vec![],
        },
        memory: vec![MemoryState {
            tau_s: 0.2,
            value: 0.7,
        }],
    }
}
fn strip() -> LinkedStrip {
    LinkedStrip::new(
        0,
        Strip {
            name: "transverse plate",
            chord_m: 1.,
            span_m: 2.,
            position_m: [0.; 3],
            normal: [1., 0., 0.],
            span_axis: [0., 1., 0.],
            control_coord: None,
            control_gain: 0.,
            control_angular_gain: [0.; 3],
        },
    )
    .unwrap()
}
fn added_mass() -> f64 {
    RHO * std::f64::consts::PI * 2. / 4.
}
fn loads(
    plate: &LinkedStrip,
    pose: Se3,
    force: Vec3,
    gravity: f64,
    map: &GeneralizedKinematics,
) -> Result<ArticulatedLoads, Refusal> {
    Ok(ArticulatedLoads {
        fluid: assemble_articulated_strips(RHO, std::slice::from_ref(plate), map)?,
        world_wrenches: vec![Wrench::new(
            Vec3::default(),
            pose.rotation().rotate(force).unwrap(),
        )],
        gravity_world: pose.rotation().rotate(v(0., 0., -gravity)).unwrap(),
        joint_effort: vec![],
    })
}
fn config() -> ContactStepConfig<'static> {
    ContactStepConfig {
        smooth: ArticulatedStepConfig {
            max_iterations: 64,
            linear_tolerance_mps: 1e-12,
            angular_tolerance_rad_s: 1e-12,
            joint_velocity_tolerances: &[],
        },
        contact: FrictionalContactConfig::default(),
    }
}
fn plane(
    id: u64,
    point: Vec3,
    pose: Se3,
    speed: f64,
    time: f64,
    mu: f64,
) -> ArticulatedFrictionalContact {
    let normal = pose.rotation().rotate(v(0., 0., 1.)).unwrap();
    let conveyor = pose.rotation().rotate(v(speed, 0., 0.)).unwrap();
    ArticulatedFrictionalContact {
        id,
        point: LinkContactPoint {
            system_index: 0,
            link_index: 0,
            point_body: point,
        },
        other: ContactCounterpart::PrescribedPlane {
            point_world: pose.transform_point(point).unwrap() - normal.scale(GAP)
                + conveyor.scale(time),
            velocity_world: conveyor,
        },
        normal_world: normal,
        tangent_world: pose.rotation().rotate(v(1., 0., 0.)).unwrap(),
        activation_tolerance: 2. * GAP,
        friction_coefficient: mu,
        normal_law: NormalVelocityLaw::Target { speed: 0. },
    }
}

// Complete interval lower bound for a root-fixed point on the ORIGINAL Lie
// midpoint path. |x''| <= |omega| (|v| + |omega| |r|); planes translate affinely.
// This is a conservative analytic bound, not a no-missed-contact claim from samples.
fn admit_root_planes(path: &ContactPath<'_>) -> Result<Vec<f64>, Refusal> {
    let u = path.midpoint.motion.base.twist_body;
    let pose = path.initial.motion.base.world_from_base;
    path.contacts
        .iter()
        .map(|c| {
            assert_eq!(c.point.link_index, 0);
            let ContactCounterpart::PrescribedPlane {
                point_world,
                velocity_world,
            } = c.other
            else {
                panic!("plane fixture")
            };
            let normal = c.normal_world.scale(1. / c.normal_world.norm());
            let point = pose.transform_point(c.point.point_body).unwrap();
            let point_speed = pose
                .rotation()
                .rotate(u.linear + u.angular.cross(c.point.point_body))
                .unwrap();
            let gap = normal.dot(point - point_world);
            let end_linear = gap + path.duration_s * normal.dot(point_speed - velocity_world);
            let acceleration =
                u.angular.norm() * (u.linear.norm() + u.angular.norm() * c.point.point_body.norm());
            let lower =
                gap.min(end_linear) - 0.5 * acceleration * path.duration_s * path.duration_s;
            if lower < 0.0 {
                Err(refused("complete Lie path cannot admit plane separation"))
            } else {
                Ok(lower)
            }
        })
        .collect()
}

#[test]
fn sustained_sliding_and_sticking_reuse_every_returned_state_with_physical_work() {
    let model = body();
    let plate = strip();
    let dt = 0.02;
    let mu = 0.3;
    let mx = MASS + added_mass();
    for (speed, drive, conveyor) in [(3., 0., 0.5), (0., 1., 0.)] {
        let mut state = initial(Se3::identity(), v(speed, 0., 0.));
        let mut expected_x = 0.;
        let mut expected_v = speed;
        for index in 0..20 {
            let saved = state.clone();
            let c = plane(
                3,
                Vec3::default(),
                Se3::identity(),
                conveyor,
                index as f64 * dt,
                mu,
            );
            let out = step_with_contacts(
                &model,
                &[],
                &state,
                &[0.2],
                index as f64 * dt,
                dt,
                &[c],
                config(),
                |_, _, map| loads(&plate, Se3::identity(), v(drive, 0., 0.), G, map),
                admit_root_planes,
                || false,
            )
            .unwrap_or_else(|error| {
                panic!(
                    "speed={speed} drive={drive} index={index} input_twist={:?}: {error:?}",
                    state.motion.base.twist_body,
                )
            });
            let pn = MASS * G * dt;
            let pt = if speed == 0. { -drive * dt } else { -mu * pn };
            let next_v = expected_v + (drive * dt + pt) / mx;
            let midpoint = 0.5 * (expected_v + next_v);
            expected_x += dt * midpoint;
            close(out.reaction_solve.impulses[0][0], pn);
            close(out.reaction_solve.impulses[0][1], pt);
            close(out.step.state.motion.base.twist_body.linear.x, next_v);
            close(out.step.state.motion.base.twist_body.linear.z, 0.);
            close(
                out.step.state.motion.base.world_from_base.translation().x,
                expected_x,
            );
            close(out.friction_work_j, pt * (midpoint - conveyor));
            close(out.support_work_j, pt * conveyor);
            close(out.reaction_work_j, pt * midpoint);
            close(out.normal_work_j, 0.);
            close(
                0.5 * mx * (next_v * next_v - expected_v * expected_v),
                drive * dt * midpoint + out.reaction_work_j,
            );
            assert!(out.endpoint_gaps[0] >= 0. && out.path_admission[0] >= 0.);
            assert!(out.initial_normal_speeds[0] >= -out.normal_speed_roundoff_mps[0]);
            assert!(out.normal_speed_roundoff_mps[0] < 1e-10);
            assert_eq!(state, saved);
            // The auxiliary midpoint projection removes a free normal load
            // response even in perfect sticking: that loss is not contact heat.
            if speed == 0. {
                assert!(out.reaction_solve.dissipated_energy > 0.);
                close(out.friction_work_j, 0.);
            }
            expected_v = next_v;
            state = out.step.state;
        }
        close(state.memory[0].value, 0.2 + 0.5 * (-20. * dt / 0.2).exp());
    }
}

#[test]
fn conveyor_motion_is_frame_covariant_and_upward_load_releases_support() {
    let model = body();
    let plate = strip();
    let dt = 0.02;
    let mu = 0.3;
    let mx = MASS + added_mass();
    for pose in [
        Se3::identity(),
        Se3::from_parts(
            So3::exp(So3Tangent::new(v(0.2, -0.1, 0.3))).unwrap(),
            v(2., -1., 0.7),
        )
        .unwrap(),
    ] {
        let state = initial(pose, v(3., 0., 0.));
        let c = plane(9, Vec3::default(), pose, 0.5, 0., mu);
        let out = step_with_contacts(
            &model,
            &[],
            &state,
            &[0.2],
            0.,
            dt,
            &[c],
            config(),
            |_, _, map| loads(&plate, pose, Vec3::default(), G, map),
            admit_root_planes,
            || false,
        )
        .unwrap();
        let pt = -mu * MASS * G * dt;
        let end = 3. + pt / mx;
        let midpoint = 0.5 * (3. + end);
        vec_close(out.step.state.motion.base.twist_body.linear, v(end, 0., 0.));
        vec_close(
            out.step.state.motion.base.world_from_base.translation(),
            pose.translation() + pose.rotation().rotate(v(dt * midpoint, 0., 0.)).unwrap(),
        );
        close(out.friction_work_j, pt * (midpoint - 0.5));
        close(out.support_work_j, pt * 0.5);
        let state = initial(pose, Vec3::default());
        let c = plane(9, Vec3::default(), pose, 0., 0., mu);
        let released = step_with_contacts(
            &model,
            &[],
            &state,
            &[0.2],
            0.,
            dt,
            &[c],
            config(),
            |_, _, map| loads(&plate, pose, v(0., 0., 2. * MASS * G), G, map),
            admit_root_planes,
            || false,
        )
        .unwrap();
        close(released.reaction_solve.impulses[0][0], 0.);
        vec_close(
            released.step.state.motion.base.twist_body.linear,
            v(0., 0., G * dt),
        );
        vec_close(
            released
                .step
                .state
                .motion
                .base
                .world_from_base
                .translation(),
            pose.translation()
                + pose
                    .rotation()
                    .rotate(v(0., 0., 0.5 * G * dt * dt))
                    .unwrap(),
        );
        assert!(released.endpoint_gaps[0] > GAP);
        close(released.reaction_work_j, 0.);
    }
}

#[test]
fn redundant_supports_retain_determined_motion_and_canonical_permutation() {
    let model = body();
    let plate = strip();
    let state = initial(Se3::identity(), v(3., 0., 0.));
    let dt = 0.02;
    let mu = 0.2;
    let mut contacts: Vec<_> = [(17, -1., -1.), (4, -1., 1.), (31, 1., -1.), (9, 1., 1.)]
        .into_iter()
        .map(|(id, x, y)| plane(id, v(x, y, 0.), Se3::identity(), 0., 0., mu))
        .collect();
    let run = |cs: &[ArticulatedFrictionalContact]| {
        step_with_contacts(
            &model,
            &[],
            &state,
            &[0.2],
            0.,
            dt,
            cs,
            config(),
            |_, _, map| loads(&plate, Se3::identity(), Vec3::default(), G, map),
            admit_root_planes,
            || false,
        )
        .unwrap()
    };
    let out = run(&contacts);
    contacts.reverse();
    let reversed = run(&contacts);
    assert_eq!(out.step.state, reversed.step.state);
    assert_eq!(out.reaction_solve, reversed.reaction_solve);
    assert_eq!(out.contact_ids, vec![4, 9, 17, 31]);
    assert!(matches!(
        out.reaction_solve.determinacy,
        FrictionalContactDeterminacy::UniqueVelocity(_)
    ));
    close(
        out.reaction_solve.impulses.iter().map(|p| p[0]).sum(),
        MASS * G * dt,
    );
    close(
        out.reaction_solve.impulses.iter().map(|p| p[1]).sum(),
        -mu * MASS * G * dt,
    );
    close(
        out.step.state.motion.base.twist_body.linear.x,
        3. - mu * MASS * G * dt / (MASS + added_mass()),
    );
    // The translated midpoint COM gives the support impulse a nonzero
    // ORBITAL moment about the fixed world origin, even with zero body spin.
    let normal_impulse = MASS * G * dt;
    let midpoint_speed = 3.0 - 0.5 * mu * normal_impulse / (MASS + added_mass());
    let midpoint_origin = v(0.5 * dt * midpoint_speed, 0.0, 0.0);
    vec_close(
        out.external_angular_impulse_world,
        midpoint_origin.cross(v(-mu * normal_impulse, 0.0, normal_impulse)),
    );
    assert!(out.endpoint_gaps.iter().all(|g| *g >= 0.));
}

#[test]
fn impacts_unresolved_friction_stops_curved_penetration_and_late_cancellation_refuse_atomically() {
    let model = body();
    let plate = strip();
    let c = plane(1, Vec3::default(), Se3::identity(), 0., 0., 0.5);
    let state = initial(Se3::identity(), v(0., 0., -0.1));
    let saved = state.clone();
    let evaluations = Cell::new(0);
    let error = step_with_contacts(
        &model,
        &[],
        &state,
        &[0.2],
        0.,
        0.02,
        &[c],
        config(),
        |_, _, map| {
            evaluations.set(evaluations.get() + 1);
            loads(&plate, Se3::identity(), Vec3::default(), G, map)
        },
        admit_root_planes,
        || false,
    )
    .unwrap_err();
    assert_eq!(error.code, "articulated-contact-impact-required");
    assert_eq!(evaluations.get(), 1);
    assert_eq!(state, saved);
    // Neither a loose physical solve tolerance nor a tiny absolute speed
    // converts a resolved incoming impact into sustained support.
    for (approach, gravity, tolerance) in [(0.01, G, 0.25), (1e-14, 0.0, 1e-9)] {
        let state = initial(Se3::identity(), v(0.0, 0.0, -approach));
        let mut cfg = config();
        cfg.contact.velocity_tolerance = tolerance;
        let error = step_with_contacts(
            &model,
            &[],
            &state,
            &[0.2],
            0.0,
            0.02,
            &[c],
            cfg,
            |_, _, map| loads(&plate, Se3::identity(), Vec3::default(), gravity, map),
            admit_root_planes,
            || false,
        )
        .unwrap_err();
        assert_eq!(error.code, "articulated-contact-impact-required");
    }
    let state = initial(Se3::identity(), v(0.1, 0., 0.));
    let saved = state.clone();
    let error = step_with_contacts(
        &model,
        &[],
        &state,
        &[0.2],
        0.,
        0.1,
        &[c],
        config(),
        |_, _, map| loads(&plate, Se3::identity(), Vec3::default(), G, map),
        admit_root_planes,
        || false,
    )
    .unwrap_err();
    assert!(error.message.contains("stopping"));
    assert_eq!(state, saved);
    let mut state = initial(
        Se3::from_parts(So3::identity(), v(0., 0., 0.1)).unwrap(),
        Vec3::default(),
    );
    state.motion.base.twist_body.angular = v(0., 1., 0.);
    let mut curved = plane(8, v(1., 0., 0.), Se3::identity(), 0., 0., 0.);
    curved.other = ContactCounterpart::PrescribedPlane {
        point_world: v(1., 0., 0.),
        velocity_world: Vec3::default(),
    };
    curved.activation_tolerance = 0.;
    let observed = Cell::new(false);
    let error = step_with_contacts(
        &model,
        &[],
        &state,
        &[0.2],
        0.,
        0.2,
        &[curved],
        config(),
        |_, _, map| loads(&plate, Se3::identity(), Vec3::default(), 0., map),
        |_| {
            observed.set(true);
            Ok(())
        },
        || false,
    )
    .unwrap_err();
    assert!(error.message.contains("penetrates"));
    assert!(!observed.get());
    // Positive initial, midpoint and endpoint gaps do not admit a path that
    // makes two complete turns through the plane between those observations.
    state.motion.base.world_from_base = Se3::identity();
    curved.point.point_body = v(0., 0., 1.);
    curved.other = ContactCounterpart::PrescribedPlane {
        point_world: v(0., 0., 0.5),
        velocity_world: Vec3::default(),
    };
    let error = step_with_contacts(
        &model,
        &[],
        &state,
        &[0.2],
        0.,
        4. * std::f64::consts::PI,
        &[curved],
        config(),
        |_, _, map| loads(&plate, Se3::identity(), Vec3::default(), 0., map),
        |path| {
            let inside = path.sample(0.25, || false)?;
            assert!(
                inside.world_from_link[0]
                    .transform_point(curved.point.point_body)
                    .unwrap()
                    .z
                    < 0.5
            );
            Err::<(), _>(refused("a true interior path point penetrates"))
        },
        || false,
    )
    .unwrap_err();
    assert_eq!(error.code, "path-falsifier");
    let state = initial(Se3::identity(), v(3., 0., 0.));
    let saved = state.clone();
    let stop = Cell::new(false);
    let error = step_with_contacts(
        &model,
        &[],
        &state,
        &[0.2],
        0.,
        0.02,
        &[c],
        config(),
        |_, _, map| loads(&plate, Se3::identity(), Vec3::default(), G, map),
        |path| {
            let receipt = admit_root_planes(path)?;
            stop.set(true);
            Ok(receipt)
        },
        || stop.get(),
    )
    .unwrap_err();
    assert_eq!(error.code, "articulated-step-cancelled");
    assert_eq!(state, saved);
}

#[test]
fn finite_gap_tangential_self_reaction_cannot_create_an_internal_couple() {
    // Large, positive physical rotary inertias keep this fixed-feature fixture
    // inside the native collinearity arithmetic envelope while making no change
    // to the analytic artificial torque gap*P_t that MUST be rejected.
    let model = ArticulatedModel::new(vec![
        Link::new(
            "base",
            None,
            Se3::identity(),
            JointModel::FIXED,
            inertia(2., 1e12),
        ),
        Link::new(
            "normal slider",
            Some(0),
            Se3::identity(),
            JointModel::prismatic(v(0., 0., 1.), None).unwrap(),
            inertia(1., 1e12),
        ),
        Link::new(
            "tangent slider",
            Some(1),
            Se3::identity(),
            JointModel::prismatic(v(1., 0., 0.), None).unwrap(),
            inertia(1., 1e12),
        ),
    ])
    .unwrap();
    let plate = strip();
    let mut state = initial(Se3::identity(), Vec3::default());
    state.motion.joint_position = vec![0.1, 0.];
    state.motion.joint_velocity = vec![0., 0.];
    let saved = state.clone();
    let c = ArticulatedFrictionalContact {
        id: 1,
        point: LinkContactPoint {
            system_index: 0,
            link_index: 2,
            point_body: Vec3::default(),
        },
        other: ContactCounterpart::Body(LinkContactPoint {
            system_index: 0,
            link_index: 0,
            point_body: Vec3::default(),
        }),
        normal_world: v(0., 0., 1.),
        tangent_world: v(1., 0., 0.),
        activation_tolerance: 0.2,
        friction_coefficient: 0.5,
        normal_law: NormalVelocityLaw::Target { speed: 0. },
    };
    let impedance = vec![
        JointImpedance {
            damping: 0.,
            stiffness: 0.,
            rest_position: 0.
        };
        2
    ];
    let mut cfg = config();
    cfg.smooth.joint_velocity_tolerances = &[1e-12, 1e-12];
    let geometry_called = Cell::new(false);
    let error = step_with_contacts(
        &model,
        &impedance,
        &state,
        &[0.2],
        0.,
        0.01,
        &[c],
        cfg,
        |_, _, map| {
            Ok(ArticulatedLoads {
                fluid: assemble_articulated_strips(RHO, std::slice::from_ref(&plate), map)?,
                world_wrenches: vec![Wrench::default(); 3],
                gravity_world: Vec3::default(),
                joint_effort: vec![-10., 2.],
            })
        },
        |_| {
            geometry_called.set(true);
            Ok(())
        },
        || false,
    )
    .unwrap_err();
    assert!(error.message.contains("couple"), "{}", error.message);
    assert!(!geometry_called.get());
    assert_eq!(state, saved);
}
