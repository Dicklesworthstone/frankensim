//! Physical coupled travel-stop events, recoil, support and release.
use fs_contact::multicontact::MultiImpactConfig;
use fs_flyer::Refusal;
use fs_flyer::addedmass::Strip;
use fs_flyer::addedmass::articulated::{LinkedStrip, assemble_articulated_strips};
use fs_flyer::partitioned::MemoryState;
use fs_flyer::partitioned::articulated::stops::{JointTravel, StopStepConfig, step_with_stops};
use fs_flyer::partitioned::articulated::{
    ArticulatedLoads, ArticulatedState, ArticulatedStepConfig, JointImpedance,
};
use fs_ga::{Mat3, Se3, Twist, Vec3, Wrench};
use fs_mbd::articulated::generalized::GeneralizedKinematics;
use fs_mbd::articulated::motion::FreeFloatingMotionState;
use fs_mbd::articulated::{
    ArticulatedModel, FreeFloatingBaseState, JointLimits, JointModel, Link, SpatialInertia,
};

const RHO: f64 = 1.2;
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * (1.0 + b.abs()),
        "{a:e} != {b:e}"
    );
}
fn inertia(value: f64) -> SpatialInertia {
    SpatialInertia::new(
        2.0,
        Vec3::default(),
        Mat3 {
            m: [value, 0.0, 0.0, 0.0, value, 0.0, 0.0, 0.0, value],
        },
    )
    .unwrap()
}
fn model(bounds: &[[f64; 2]]) -> ArticulatedModel {
    let mut links = vec![Link::new(
        "carrier",
        None,
        Se3::identity(),
        JointModel::FIXED,
        inertia(5.0),
    )];
    for (i, bound) in bounds.iter().enumerate() {
        let joint = JointModel::revolute(
            Vec3::new(0.0, 1.0, 0.0),
            Some(JointLimits::new(bound[0], bound[1], 1000.0, 1e9).unwrap()),
        )
        .unwrap();
        links.push(Link::new(
            format!("rotor {i}"),
            Some(i),
            Se3::identity(),
            joint,
            inertia(2.0 + i as f64),
        ));
    }
    ArticulatedModel::new(links).unwrap()
}
fn strips(n: usize) -> Vec<LinkedStrip> {
    (0..n)
        .map(|i| {
            LinkedStrip::new(
                i + 1,
                Strip {
                    name: "plate",
                    chord_m: 1.1,
                    span_m: 2.0,
                    position_m: [0.0; 3],
                    normal: [0.0, 0.0, 1.0],
                    span_axis: [0.0, 1.0, 0.0],
                    control_coord: None,
                    control_gain: 0.0,
                    control_angular_gain: [0.0; 3],
                },
            )
            .unwrap()
        })
        .collect()
}
fn fluid_pitch() -> f64 {
    RHO * std::f64::consts::PI * 1.1_f64.powi(4) * 2.0 / 128.0
}
fn state(q: &[f64], v: &[f64], omega: f64) -> ArticulatedState {
    ArticulatedState {
        motion: FreeFloatingMotionState {
            base: FreeFloatingBaseState::new(
                Se3::identity(),
                Twist::new(Vec3::new(0.0, omega, 0.0), Vec3::default()),
            ),
            joint_position: q.to_vec(),
            joint_velocity: v.to_vec(),
        },
        memory: vec![MemoryState {
            tau_s: 0.02,
            value: 0.4,
        }],
    }
}
fn loads(
    plates: &[LinkedStrip],
    links: usize,
    effort: &[f64],
    map: &GeneralizedKinematics,
) -> Result<ArticulatedLoads, Refusal> {
    Ok(ArticulatedLoads {
        fluid: assemble_articulated_strips(RHO, plates, map)?,
        world_wrenches: vec![Wrench::default(); links],
        gravity_world: Vec3::default(),
        joint_effort: effort.to_vec(),
    })
}
fn impedance(n: usize) -> Vec<JointImpedance> {
    vec![
        JointImpedance {
            damping: 0.0,
            stiffness: 0.0,
            rest_position: 0.0
        };
        n
    ]
}
fn travel(joint: usize, lower: f64, upper: f64, restitution: f64) -> JointTravel {
    JointTravel {
        joint,
        lower,
        upper,
        activation_tolerance: 1e-9,
        restitution,
    }
}
fn config(tolerances: &[f64]) -> StopStepConfig<'_> {
    let mut contact = MultiImpactConfig::default();
    contact.solver.gradient_tolerance = 1e-12;
    StopStepConfig {
        smooth: ArticulatedStepConfig {
            max_iterations: 64,
            linear_tolerance_mps: 1e-12,
            angular_tolerance_rad_s: 1e-12,
            joint_velocity_tolerances: tolerances,
        },
        maximum_step_s: 1.0,
        minimum_interval_s: 1e-13,
        maximum_attempts: 2048,
        maximum_events: 16,
        contact,
    }
}

#[test]
fn localized_rebound_uses_complete_mass_and_preserves_internal_momentum() {
    let model = model(&[[-0.4, 0.3]]);
    let plates = strips(1);
    let k = 2.0 + fluid_pitch();
    let ratio = k / (5.0 + k);
    let initial = state(&[0.0], &[2.0], -2.0 * ratio);
    let saved = initial.clone();
    let result = step_with_stops(
        &model,
        &impedance(1),
        &initial,
        &[1.0],
        0.0,
        0.25,
        &[travel(0, -0.4, 0.3, 0.5)],
        config(&[1e-12]),
        |_, _, map| loads(&plates, model.link_count(), &[0.0], map),
        || false,
    )
    .unwrap();
    assert_eq!(result.events.len(), 1);
    close(result.events[0].time_s, 0.15, 2e-9);
    close(result.state.motion.joint_position[0], 0.2, 3e-9);
    close(result.state.motion.joint_velocity[0], -1.0, 2e-12);
    close(result.state.motion.base.twist_body.angular.y, ratio, 2e-12);
    let schur = 5.0 * k / (5.0 + k);
    close(
        result.impact_dissipation_j,
        0.5 * schur * 4.0 * (1.0 - 0.25),
        3e-12,
    );
    close(result.events[0].energy_balance_residual_j, 0.0, 3e-12);
    close(
        (5.0 + k) * result.state.motion.base.twist_body.angular.y
            + k * result.state.motion.joint_velocity[0],
        0.0,
        3e-12,
    );
    close(
        result.state.memory[0].value,
        1.0 - 0.6 * (-0.25_f64 / 0.02).exp(),
        2e-14,
    );
    assert_eq!(initial, saved);
}

#[test]
fn simultaneous_stops_share_one_coupled_impulse_and_locked_momentum() {
    let model = model(&[[-1.0, 1.0], [-1.0, 1.0]]);
    let plates = strips(2);
    let initial = state(&[0.0, 0.0], &[1.0, 2.0], 0.0);
    let k1 = 2.0 + fluid_pitch();
    let k2 = 3.0 + fluid_pitch();
    let total = 5.0 + k1 + k2;
    let momentum = k1 * 1.0 + k2 * 3.0;
    let energy = 0.5 * k1 + 0.5 * k2 * 9.0;
    let result = step_with_stops(
        &model,
        &impedance(2),
        &initial,
        &[1.0],
        0.0,
        0.3,
        &[travel(0, -0.5, 0.2, 0.0), travel(1, -0.6, 0.4, 0.0)],
        config(&[1e-12; 2]),
        |_, _, map| loads(&plates, model.link_count(), &[0.0; 2], map),
        || false,
    )
    .unwrap();
    assert_eq!(result.events.len(), 1);
    assert_eq!(result.events[0].normal_impulses.len(), 2);
    assert!(result.events[0].normal_impulses.iter().all(|j| *j > 0.1));
    close(result.events[0].time_s, 0.2, 2e-9);
    close(result.state.motion.joint_velocity[0], 0.0, 0.0);
    close(result.state.motion.joint_velocity[1], 0.0, 0.0);
    close(
        result.state.motion.base.twist_body.angular.y,
        momentum / total,
        2e-12,
    );
    close(
        result.impact_dissipation_j,
        energy - 0.5 * momentum * momentum / total,
        4e-12,
    );
    close(result.state.motion.joint_position[0], 0.2, 2e-9);
    close(result.state.motion.joint_position[1], 0.4, 2e-9);
}

#[test]
fn stiff_stop_holds_releases_and_recontacts_with_full_fluid_recoil() {
    let model = model(&[[-0.4, 0.3]]);
    let plates = strips(1);
    let impedance = [JointImpedance {
        damping: 0.0,
        stiffness: 1e6,
        rest_position: 0.6,
    }];
    let stop = [travel(0, -0.4, 0.3, 0.0)];
    let initial = state(&[0.3], &[0.0], 0.0);
    let held = step_with_stops(
        &model,
        &impedance,
        &initial,
        &[1.0],
        0.0,
        0.1,
        &stop,
        config(&[1e-12]),
        |_, _, map| loads(&plates, model.link_count(), &[0.0], map),
        || false,
    )
    .unwrap();
    assert!(held.events.is_empty());
    assert_eq!(held.state.motion.joint_position[0], 0.3);
    assert_eq!(held.state.motion.joint_velocity[0], 0.0);
    close(held.reaction_impulses[0][1], 30_000.0, 3e-12);
    close(held.reaction_work_j, 0.0, 1e-12);
    let mut resolved = config(&[1e-12]);
    resolved.maximum_step_s = 1e-4;
    let released = step_with_stops(
        &model,
        &impedance,
        &held.state,
        &[1.0],
        0.1,
        0.0002,
        &stop,
        resolved,
        |_, _, map| loads(&plates, model.link_count(), &[-400_000.0], map),
        || false,
    )
    .unwrap();
    assert!(released.state.motion.joint_position[0] < 0.3 - 1e-4);
    assert!(released.state.motion.joint_velocity[0] < -1.0);
    close(released.reaction_impulses[0][1], 0.0, 0.0);
    let returned = step_with_stops(
        &model,
        &impedance,
        &released.state,
        &[1.0],
        0.1002,
        0.005,
        &stop,
        resolved,
        |_, _, map| loads(&plates, model.link_count(), &[0.0], map),
        || false,
    )
    .unwrap();
    assert_eq!(returned.events.len(), 1);
    assert!(returned.impact_dissipation_j > 1.0);
    close(returned.state.motion.joint_position[0], 0.3, 2e-9);
    close(returned.state.motion.joint_velocity[0], 0.0, 0.0);
    let k = 2.0 + fluid_pitch();
    close(
        (5.0 + k) * returned.state.motion.base.twist_body.angular.y
            + k * returned.state.motion.joint_velocity[0],
        0.0,
        2e-9,
    );
    close(
        returned.state.memory[0].value,
        1.0 - 0.6 * (-0.1052_f64 / 0.02).exp(),
        4e-14,
    );
}

#[test]
fn event_budgets_and_late_load_refusals_publish_no_partial_motion_or_memory() {
    let model = model(&[[-0.4, 0.3]]);
    let plates = strips(1);
    let initial = state(&[0.0], &[2.0], 0.0);
    let saved = initial.clone();
    let stops = [travel(0, -0.4, 0.3, 0.5)];
    let mut no_events = config(&[1e-12]);
    no_events.maximum_events = 0;
    let error = step_with_stops(
        &model,
        &impedance(1),
        &initial,
        &[1.0],
        0.0,
        0.3,
        &stops,
        no_events,
        |_, _, map| loads(&plates, model.link_count(), &[0.0], map),
        || false,
    )
    .unwrap_err();
    assert_eq!(error.code, "joint-stop-unresolved");
    assert_eq!(initial, saved);
    let error = step_with_stops(
        &model,
        &impedance(1),
        &initial,
        &[1.0],
        0.0,
        0.3,
        &stops,
        config(&[1e-12]),
        |time, _, map| {
            if time > 0.16 {
                return Err(Refusal {
                    code: "fixture-load-refused",
                    message: "late source refusal".into(),
                    ranked_repairs: vec![],
                });
            }
            loads(&plates, model.link_count(), &[0.0], map)
        },
        || false,
    )
    .unwrap_err();
    assert_eq!(error.code, "fixture-load-refused");
    assert_eq!(initial, saved);
    let error = step_with_stops(
        &model,
        &impedance(1),
        &initial,
        &[1.0],
        0.0,
        0.3,
        &stops,
        config(&[1e-12]),
        |_, _, map| loads(&plates, model.link_count(), &[0.0], map),
        || true,
    )
    .unwrap_err();
    assert_eq!(error.code, "articulated-step-cancelled");
    assert_eq!(initial, saved);
    let invalid = [travel(0, 0.1, 0.3, 0.5)];
    assert!(
        step_with_stops(
            &model,
            &impedance(1),
            &initial,
            &[1.0],
            0.0,
            0.3,
            &invalid,
            config(&[1e-12]),
            |_, _, map| loads(&plates, model.link_count(), &[0.0], map),
            || false
        )
        .is_err()
    );
}

#[test]
fn wide_trials_detect_interior_crossing_and_departure_then_return() {
    let model = model(&[[-2.0, 2.0]]);
    let plates = strips(1);
    let k = 2.0 + fluid_pitch();
    let effective = 5.0 * k / (5.0 + k);
    let ratio = k / (5.0 + k);
    let stop = [travel(0, -0.7, 0.3, 0.0)];
    // q=4t-10t² has an admissible final value at .32 s despite crossing
    // the upper stop at .1 s and returning through it at .3 s.
    let initial = state(&[0.0], &[4.0], -4.0 * ratio);
    let crossed = step_with_stops(
        &model,
        &impedance(1),
        &initial,
        &[1.0],
        0.0,
        0.32,
        &stop,
        config(&[1e-12]),
        |_, _, map| loads(&plates, model.link_count(), &[-20.0 * effective], map),
        || false,
    )
    .unwrap();
    assert_eq!(crossed.events.len(), 1);
    close(crossed.events[0].time_s, 0.1, 3e-9);
    close(crossed.state.motion.joint_position[0], -0.184, 3e-8);
    close(crossed.state.motion.joint_velocity[0], -4.4, 3e-8);
    close(
        crossed.state.motion.base.twist_body.angular.y,
        4.4 * ratio,
        3e-8,
    );
    close(crossed.impact_dissipation_j, 2.0 * effective, 4e-8);

    // Start AT the stop and genuinely depart. A .25 s trial must find
    // the .2 s closing return before supporting the last .05 s.
    let departing = state(&[0.3], &[-2.0], 2.0 * ratio);
    let returned = step_with_stops(
        &model,
        &impedance(1),
        &departing,
        &[1.0],
        0.0,
        0.25,
        &stop,
        config(&[1e-12]),
        |_, _, map| loads(&plates, model.link_count(), &[20.0 * effective], map),
        || false,
    )
    .unwrap();
    assert_eq!(returned.events.len(), 1);
    close(returned.events[0].time_s, 0.2, 3e-9);
    close(returned.state.motion.joint_position[0], 0.3, 2e-9);
    close(returned.state.motion.joint_velocity[0], 0.0, 0.0);
    close(returned.state.motion.base.twist_body.angular.y, 0.0, 3e-11);
    close(
        returned.reaction_impulses[0][1],
        20.0 * effective * 0.05,
        3e-8,
    );
    close(returned.impact_dissipation_j, 2.0 * effective, 4e-8);
    close(
        returned.state.memory[0].value,
        1.0 - 0.6 * (-0.25_f64 / 0.02).exp(),
        4e-14,
    );
}

#[test]
fn simultaneous_support_and_release_retain_forced_base_and_joint_coupling() {
    let model = model(&[[-1.0, 1.0], [-1.0, 1.0]]);
    let plates = strips(2);
    let k1 = 2.0 + fluid_pitch();
    let k2 = 3.0 + fluid_pitch();
    let initial = state(&[0.2, 0.4], &[0.0; 2], 0.0);
    let stops = [travel(0, -0.5, 0.2, 0.0), travel(1, -0.6, 0.4, 0.0)];
    let torque = 15.0;
    let dt = 0.01;
    for effort2 in [6.0, 2.0] {
        let (alpha, acceleration2, reaction1, reaction2) = if effort2 == 6.0 {
            let a = torque / (5.0 + k1 + k2);
            (a, 0.0, 20.0 - (k1 + k2) * a, effort2 - k2 * a)
        } else {
            let a = (torque - effort2) / (5.0 + k1);
            (a, effort2 / k2 - a, 20.0 - k1 * a - effort2, 0.0)
        };
        let result = step_with_stops(
            &model,
            &impedance(2),
            &initial,
            &[1.0],
            0.0,
            dt,
            &stops,
            config(&[1e-12; 2]),
            |_, _, map| {
                let mut result = loads(&plates, model.link_count(), &[20.0, effort2], map)?;
                result.world_wrenches[0].torque.y = torque;
                Ok(result)
            },
            || false,
        )
        .unwrap();
        assert!(result.events.is_empty());
        close(
            result.state.motion.base.twist_body.angular.y,
            dt * alpha,
            4e-12,
        );
        assert_eq!(result.state.motion.joint_position[0], 0.2);
        assert_eq!(result.state.motion.joint_velocity[0], 0.0);
        close(
            result.state.motion.joint_position[1],
            0.4 + 0.5 * dt * dt * acceleration2,
            4e-12,
        );
        close(
            result.state.motion.joint_velocity[1],
            dt * acceleration2,
            4e-12,
        );
        close(result.reaction_impulses[0][1], dt * reaction1, 4e-12);
        close(result.reaction_impulses[1][1], dt * reaction2, 4e-12);
        close(result.reaction_work_j, 0.0, 4e-13);
        let momentum = (5.0 + k1 + k2) * result.state.motion.base.twist_body.angular.y
            + (k1 + k2) * result.state.motion.joint_velocity[0]
            + k2 * result.state.motion.joint_velocity[1];
        close(momentum, torque * dt, 4e-12);
    }
}
