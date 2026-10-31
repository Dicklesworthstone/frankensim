//! Independent multi-tree impact oracles with physical attached fluid inertia.
use fs_contact::multicontact::MultiImpactConfig;
use fs_flyer::Refusal;
use fs_flyer::addedmass::Strip;
use fs_flyer::addedmass::articulated::{LinkedStrip, assemble_articulated_strips};
use fs_flyer::partitioned::MemoryState;
use fs_flyer::partitioned::articulated::ArticulatedState;
use fs_flyer::partitioned::articulated::impact::{
    ArticulatedImpactSystem, ArticulatedNormalContact, ContactCounterpart, LinkContactPoint,
    resolve_many,
};
use fs_ga::{Mat3, Se3, So3, Twist, Vec3};
use fs_mbd::articulated::motion::FreeFloatingMotionState;
use fs_mbd::articulated::{
    ArticulatedModel, FreeFloatingBaseState, JointLimits, JointModel, Link, SpatialInertia,
};

fn v(x: f64, y: f64, z: f64) -> Vec3 {
    Vec3::new(x, y, z)
}
fn close(a: f64, b: f64) {
    assert!((a - b).abs() < 3e-11 * (1.0 + b.abs()), "{a:e} != {b:e}");
}
fn inertia(mass: f64, spin: f64) -> SpatialInertia {
    SpatialInertia::new(
        mass,
        Vec3::default(),
        Mat3 {
            m: [spin, 0.0, 0.0, 0.0, spin, 0.0, 0.0, 0.0, spin],
        },
    )
    .unwrap()
}
fn rotor(rate_limit: f64) -> ArticulatedModel {
    ArticulatedModel::new(vec![
        Link::new(
            "carrier",
            None,
            Se3::identity(),
            JointModel::FIXED,
            inertia(3.0, 5.0),
        ),
        Link::new(
            "rotor",
            Some(0),
            Se3::identity(),
            JointModel::revolute(
                v(0.0, 1.0, 0.0),
                Some(JointLimits::new(-1.0, 1.0, rate_limit, 1000.0).unwrap()),
            )
            .unwrap(),
            inertia(2.0, 2.0),
        ),
    ])
    .unwrap()
}
fn state(x: f64, base_velocity: Vec3, rates: &[f64]) -> ArticulatedState {
    ArticulatedState {
        motion: FreeFloatingMotionState {
            base: FreeFloatingBaseState::new(
                Se3::from_parts(So3::identity(), v(x, 0.0, 0.0)).unwrap(),
                Twist::new(Vec3::default(), base_velocity),
            ),
            joint_position: vec![0.0; rates.len()],
            joint_velocity: rates.to_vec(),
        },
        memory: vec![MemoryState {
            tau_s: 0.2,
            value: 0.6,
        }],
    }
}
fn plate(link: usize, normal: [f64; 3]) -> LinkedStrip {
    LinkedStrip::new(
        link,
        Strip {
            name: "independent attached plate",
            chord_m: 1.2,
            span_m: 1.6,
            position_m: [0.0; 3],
            normal,
            span_axis: [0.0, 1.0, 0.0],
            control_coord: None,
            control_gain: 0.0,
            control_angular_gain: [0.0; 3],
        },
    )
    .unwrap()
}
fn fluid_mass(rho: f64) -> f64 {
    rho * std::f64::consts::PI * 1.2_f64.powi(2) * 1.6 / 4.0
}
fn point(system: usize, link: usize, x: f64) -> LinkContactPoint {
    LinkContactPoint {
        system_index: system,
        link_index: link,
        point_body: v(x, 0.0, 0.0),
    }
}
fn pair(
    first: LinkContactPoint,
    second: LinkContactPoint,
    normal: Vec3,
    restitution: f64,
) -> ArticulatedNormalContact {
    ArticulatedNormalContact {
        point: first,
        other: ContactCounterpart::Body(second),
        normal_world: normal,
        activation_tolerance: 0.0,
        restitution,
    }
}

#[test]
fn separate_fluid_rotors_exchange_impulse_with_exact_joint_and_orbital_momentum() {
    let models = [rotor(100.0), rotor(100.0), rotor(100.0)];
    let states = [
        state(0.0, Vec3::default(), &[1.0]),
        state(2.0, Vec3::default(), &[0.0]),
        state(7.0, v(0.2, 0.0, 0.0), &[0.3]),
    ];
    let saved = states.clone();
    let rho = [1.0, 1.4, 0.8];
    let plates = [plate(1, [0.0, 0.0, 1.0])];
    let systems: Vec<_> = models
        .iter()
        .zip(&states)
        .map(|(model, state)| ArticulatedImpactSystem { model, state })
        .collect();
    let contact = pair(point(0, 1, 1.0), point(1, 1, -1.0), v(0.0, 0.0, 1.0), 0.5);
    let event = resolve_many(
        &systems,
        &[contact],
        MultiImpactConfig::default(),
        |i, map| assemble_articulated_strips(rho[i], &plates, map),
        || false,
    )
    .unwrap();
    let m = [5.0 + fluid_mass(rho[0]), 5.0 + fluid_mass(rho[1])];
    let k = [
        2.0 + fluid_mass(rho[0]) * 1.2_f64.powi(2) / 32.0,
        2.0 + fluid_mass(rho[1]) * 1.2_f64.powi(2) / 32.0,
    ];
    let mobility = 1.0 / m[0] + 1.0 / k[0] + 1.0 / m[1] + 1.0 / k[1];
    let impulse = 1.5 / mobility;
    close(event.delassus[0], mobility);
    close(event.law.normal_impulses[0], impulse);
    close(
        event.mechanics.states[0].motion.base.twist_body.linear.z,
        impulse / m[0],
    );
    close(
        event.mechanics.states[1].motion.base.twist_body.linear.z,
        -impulse / m[1],
    );
    close(
        event.mechanics.states[0].motion.joint_velocity[0],
        1.0 - impulse / k[0],
    );
    close(
        event.mechanics.states[1].motion.joint_velocity[0],
        -impulse / k[1],
    );
    close(event.mechanics.relative_speeds_after[0], 0.5);
    close(event.law.dissipated_energy, 0.5 * (1.0 - 0.25) / mobility);
    let a = &event.mechanics.states[0].motion;
    let b = &event.mechanics.states[1].motion;
    close(
        m[0] * a.base.twist_body.linear.z + m[1] * b.base.twist_body.linear.z,
        0.0,
    );
    close(
        k[0] * a.joint_velocity[0] + k[1] * b.joint_velocity[0]
            - 2.0 * m[1] * b.base.twist_body.linear.z,
        k[0],
    );
    close(
        a.base.twist_body.angular.norm() + b.base.twist_body.angular.norm(),
        0.0,
    );
    assert_eq!(event.mechanics.states[2], states[2]);
    assert_eq!(states, saved);

    // Reindexing the independent systems changes no physical response.
    let reversed = [systems[1], systems[0], systems[2]];
    let swapped = pair(point(1, 1, 1.0), point(0, 1, -1.0), v(0.0, 0.0, 1.0), 0.5);
    let other = resolve_many(
        &reversed,
        &[swapped],
        MultiImpactConfig::default(),
        |i, map| assemble_articulated_strips(rho[[1, 0, 2][i]], &plates, map),
        || false,
    )
    .unwrap();
    for (left, right) in event.mechanics.velocities_after.iter().zip([
        &other.mechanics.velocities_after[1],
        &other.mechanics.velocities_after[0],
        &other.mechanics.velocities_after[2],
    ]) {
        for (a, b) in left.iter().zip(right) {
            close(*a, *b);
        }
    }
}

#[test]
fn one_global_three_body_solve_retains_negative_cross_contact_mobility() {
    let dry_mass = [2.0, 3.0, 4.0];
    let rho = [0.5, 1.0, 1.5];
    let models: Vec<_> = dry_mass
        .iter()
        .map(|m| {
            ArticulatedModel::new(vec![Link::new(
                "body",
                None,
                Se3::identity(),
                JointModel::FIXED,
                inertia(*m, 2.0),
            )])
            .unwrap()
        })
        .collect();
    let initial = [2.0, 0.0, -1.0];
    let states: Vec<_> = initial
        .iter()
        .map(|u| state(0.0, v(*u, 0.0, 0.0), &[]))
        .collect();
    let systems: Vec<_> = models
        .iter()
        .zip(&states)
        .map(|(model, state)| ArticulatedImpactSystem { model, state })
        .collect();
    let plates = [plate(0, [1.0, 0.0, 0.0])];
    let contacts = [
        pair(point(1, 0, 0.0), point(0, 0, 0.0), v(1.0, 0.0, 0.0), 0.0),
        pair(point(2, 0, 0.0), point(1, 0, 0.0), v(1.0, 0.0, 0.0), 0.0),
    ];
    let event = resolve_many(
        &systems,
        &contacts,
        MultiImpactConfig::default(),
        |i, map| assemble_articulated_strips(rho[i], &plates, map),
        || false,
    )
    .unwrap();
    let m: [f64; 3] = std::array::from_fn(|i| dry_mass[i] + fluid_mass(rho[i]));
    let final_speed = (2.0 * m[0] - m[2]) / m.iter().sum::<f64>();
    close(event.delassus[1], -1.0 / m[1]);
    close(event.delassus[2], -1.0 / m[1]);
    close(event.law.normal_impulses[0], m[0] * (2.0 - final_speed));
    close(event.law.normal_impulses[1], m[2] * (final_speed + 1.0));
    for state in &event.mechanics.states {
        close(state.motion.base.twist_body.linear.x, final_speed);
    }
    let before = 0.5 * (4.0 * m[0] + m[2]);
    let after = 0.5 * m.iter().sum::<f64>() * final_speed.powi(2);
    close(event.law.dissipated_energy, before - after);
    close(event.mechanics.kinetic_before, before);
    close(event.mechanics.kinetic_after, after);
    let swapped = resolve_many(
        &systems,
        &[contacts[1], contacts[0]],
        MultiImpactConfig::default(),
        |i, map| assemble_articulated_strips(rho[i], &plates, map),
        || false,
    )
    .unwrap();
    for state in &swapped.mechanics.states {
        close(state.motion.base.twist_body.linear.x, final_speed);
    }
}

#[test]
fn multi_tree_failures_are_atomic_including_late_native_rate_admission() {
    let models = [rotor(100.0), rotor(0.1)];
    let states = [
        state(0.0, Vec3::default(), &[1.0]),
        state(2.0, Vec3::default(), &[0.0]),
    ];
    let saved = states.clone();
    let systems: Vec<_> = models
        .iter()
        .zip(&states)
        .map(|(model, state)| ArticulatedImpactSystem { model, state })
        .collect();
    let contacts = [pair(
        point(0, 1, 1.0),
        point(1, 1, -1.0),
        v(0.0, 0.0, 1.0),
        0.5,
    )];
    let plates = [plate(1, [0.0, 0.0, 1.0])];
    let failure = resolve_many(
        &systems,
        &contacts,
        MultiImpactConfig::default(),
        |_, map| assemble_articulated_strips(1.0, &plates, map),
        || false,
    )
    .unwrap_err();
    assert_eq!(failure.code, "articulated-impact-invalid");
    assert!(failure.message.starts_with("joint speed "), "{failure:?}");
    assert!(
        failure
            .message
            .contains("for link 1 exceeds absolute limit 0.1"),
        "{failure:?}"
    );
    let failure = resolve_many(
        &systems,
        &contacts,
        MultiImpactConfig::default(),
        |i, map| {
            if i == 1 {
                Err(Refusal {
                    code: "fixture-fluid-refused",
                    message: "late fluid refusal".into(),
                    ranked_repairs: vec![],
                })
            } else {
                assemble_articulated_strips(1.0, &plates, map)
            }
        },
        || false,
    )
    .unwrap_err();
    assert_eq!(failure.code, "fixture-fluid-refused");
    let failure = resolve_many(
        &systems,
        &contacts,
        MultiImpactConfig::default(),
        |_, map| assemble_articulated_strips(1.0, &plates, map),
        || true,
    )
    .unwrap_err();
    assert_eq!(failure.code, "articulated-impact-cancelled");
    let failure = resolve_many(
        &vec![systems[0]; 129],
        &[],
        MultiImpactConfig::default(),
        |_, _| panic!("oversized event must refuse before fluid evaluation"),
        || false,
    )
    .unwrap_err();
    assert_eq!(failure.code, "articulated-impact-invalid");
    assert_eq!(states, saved);
}
