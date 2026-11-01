//! Independent sustained-contact oracles across complete dry/fluid trees.
use fs_contact::frictional_multicontact::{FrictionalContactConfig, NormalVelocityLaw};
use fs_flyer::Refusal;
use fs_flyer::addedmass::{
    Strip,
    articulated::{LinkedStrip, assemble_articulated_strips},
};
use fs_flyer::partitioned::MemoryState;
use fs_flyer::partitioned::articulated::contact::ArticulatedFrictionalContact;
use fs_flyer::partitioned::articulated::contact::systems::{
    ContactSystem, ContactSystemsPath, step_with_contacts,
};
use fs_flyer::partitioned::articulated::{
    ArticulatedLoads, ArticulatedStage, ArticulatedState, ArticulatedStepConfig, step,
};
use fs_ga::{Mat3, Se3, So3, Twist, Vec3, Wrench};
use fs_mbd::articulated::contact_geometry::{ContactCounterpart, LinkContactPoint};
use fs_mbd::articulated::motion::FreeFloatingMotionState;
use fs_mbd::articulated::{
    ArticulatedModel, FreeFloatingBaseState, JointModel, Link, SpatialInertia,
};
use std::cell::Cell;

const RHO: f64 = 1.2;
const G: f64 = 9.8;
const GAP: f64 = 1.0 / 16_777_216.0;
const HELD: [f64; 1] = [0.2];
fn v(x: f64, y: f64, z: f64) -> Vec3 {
    Vec3::new(x, y, z)
}
fn close(a: f64, b: f64) {
    assert!((a - b).abs() < 3e-9 * (1.0 + b.abs()), "{a:e} != {b:e}");
}
fn refusal(message: &str) -> Refusal {
    Refusal {
        code: "path-falsifier",
        message: message.into(),
        ranked_repairs: vec![],
    }
}
fn model(mass: f64) -> ArticulatedModel {
    ArticulatedModel::new(vec![Link::new(
        "body",
        None,
        Se3::identity(),
        JointModel::FIXED,
        SpatialInertia::new(
            mass,
            Vec3::default(),
            Mat3 {
                m: [4., 0., 0., 0., 4., 0., 0., 0., 4.],
            },
        )
        .unwrap(),
    )])
    .unwrap()
}
fn state(position: Vec3, linear: Vec3) -> ArticulatedState {
    ArticulatedState {
        motion: FreeFloatingMotionState {
            base: FreeFloatingBaseState::new(
                Se3::from_parts(So3::identity(), position).unwrap(),
                Twist::new(Vec3::default(), linear),
            ),
            joint_position: vec![],
            joint_velocity: vec![],
        },
        memory: vec![MemoryState {
            tau_s: 0.2,
            value: 0.7,
        }],
    }
}
fn strip(span: f64, normal: Vec3) -> LinkedStrip {
    LinkedStrip::new(
        0,
        Strip {
            name: "attached plate",
            chord_m: 1.,
            span_m: span,
            position_m: [0.; 3],
            normal: [normal.x, normal.y, normal.z],
            span_axis: if normal.x == 0. {
                [1., 0., 0.]
            } else {
                [0., 1., 0.]
            },
            control_coord: None,
            control_gain: 0.,
            control_angular_gain: [0.; 3],
        },
    )
    .unwrap()
}
fn fluid_mass(span: f64) -> f64 {
    RHO * std::f64::consts::PI * span / 4.
}
fn smooth() -> ArticulatedStepConfig<'static> {
    ArticulatedStepConfig {
        max_iterations: 64,
        linear_tolerance_mps: 1e-12,
        angular_tolerance_rad_s: 1e-12,
        joint_velocity_tolerances: &[],
    }
}
fn systems<'a>(
    models: &'a [ArticulatedModel],
    states: &'a [ArticulatedState],
) -> Vec<ContactSystem<'a>> {
    models
        .iter()
        .zip(states)
        .map(|(model, state)| ContactSystem {
            model,
            state,
            impedance: &[],
            memory_inputs: &HELD,
            smooth: smooth(),
        })
        .collect()
}
fn loads(
    stage: ArticulatedStage<'_>,
    plate: &LinkedStrip,
    force: Vec3,
    g: f64,
) -> Result<ArticulatedLoads, Refusal> {
    Ok(ArticulatedLoads {
        fluid: assemble_articulated_strips(RHO, std::slice::from_ref(plate), stage.kinematics)?,
        world_wrenches: vec![Wrench::new(Vec3::default(), force)],
        gravity_world: v(0., 0., -g),
        joint_effort: vec![],
    })
}
fn contact(
    id: u64,
    first: usize,
    second: usize,
    arm: f64,
    mu: f64,
) -> ArticulatedFrictionalContact {
    ArticulatedFrictionalContact {
        id,
        point: LinkContactPoint {
            system_index: first,
            link_index: 0,
            point_body: v(0., 0., -arm),
        },
        other: ContactCounterpart::Body(LinkContactPoint {
            system_index: second,
            link_index: 0,
            point_body: v(0., 0., arm),
        }),
        normal_world: v(0., 0., 1.),
        tangent_world: v(1., 0., 0.),
        activation_tolerance: 2. * GAP,
        friction_coefficient: mu,
        normal_law: NormalVelocityLaw::Target { speed: 0. },
    }
}

// Fixed world normals and nonrotating root translations make every gap affine
// over the COMPLETE numerical path. Its two endpoint gaps give the exact lower
// bound in real arithmetic; sampling is only an API consistency check here.
fn admit_affine(path: &ContactSystemsPath<'_>) -> Result<Vec<[f64; 2]>, Refusal> {
    for midpoint in path.midpoints {
        if midpoint.motion.base.twist_body.angular.norm() != 0.
            || !midpoint.motion.joint_position.is_empty()
        {
            return Err(refusal("fixture is no longer a pure translation"));
        }
    }
    let samples = path.sample(0.375, || false)?;
    assert_eq!(samples.len(), path.systems.len());
    path.contacts
        .iter()
        .map(|c| {
            let ContactCounterpart::Body(other) = c.other else {
                panic!("body-pair fixture")
            };
            let point = |s: &ArticulatedState, p: Vec3| {
                s.motion.base.world_from_base.transform_point(p).unwrap()
            };
            let gap = |first: &ArticulatedState, second: &ArticulatedState| {
                c.normal_world
                    .dot(point(first, c.point.point_body) - point(second, other.point_body))
            };
            let before = gap(
                path.systems[c.point.system_index].state,
                path.systems[other.system_index].state,
            );
            let after = gap(
                path.endpoints[c.point.system_index],
                path.endpoints[other.system_index],
            );
            if before < 0. || after < 0. {
                Err(refusal("affine gap crosses contact"))
            } else {
                Ok([before, after])
            }
        })
        .collect()
}

#[test]
fn repeated_two_tree_support_uses_both_fluid_masses_and_physical_work() {
    let dry = [1., 2.];
    let spans = [0.5, 3.];
    let masses = [dry[0] + fluid_mass(spans[0]), dry[1] + fluid_mass(spans[1])];
    let models = [model(dry[0]), model(dry[1])];
    let plates = [
        strip(spans[0], v(0., 0., 1.)),
        strip(spans[1], v(0., 0., 1.)),
    ];
    let mut states = vec![
        state(v(2., -3., 1. + GAP), v(0., 0., -0.125)),
        state(v(2., -3., 0.), v(0., 0., -0.125)),
    ];
    let acceleration = -G * (dry[0] + dry[1]) / (masses[0] + masses[1]);
    let reaction = dry[0] * G + masses[0] * acceleration;
    assert!(reaction > 0.);
    let dt = 0.02;
    let c = contact(7, 0, 1, 0.5, 0.);
    let mut expected_v = -0.125;
    for tick in 0..20 {
        let before = states.clone();
        let out = step_with_contacts(
            &systems(&models, &states),
            tick as f64 * dt,
            dt,
            &[c],
            FrictionalContactConfig::default(),
            |time, stages| {
                close(time, (tick as f64 + 0.5) * dt);
                stages
                    .iter()
                    .enumerate()
                    .map(|(i, stage)| {
                        close(
                            stage.state.memory[0].value,
                            states[i].memory[0].advanced(HELD[0], 0.5 * dt).value,
                        );
                        loads(*stage, &plates[i], Vec3::default(), G)
                    })
                    .collect()
            },
            admit_affine,
            || false,
        )
        .unwrap();
        close(out.reaction_solve.impulses[0][0], dt * reaction);
        close(out.reaction_work_j, 0.);
        close(out.normal_work_j, 0.);
        close(out.friction_work_j, 0.);
        close(out.support_work_j, 0.);
        close(out.external_linear_impulse_world.norm(), 0.);
        close(out.external_angular_impulse_world.norm(), 0.);
        let next_v = expected_v + dt * acceleration;
        let mut energy_change = 0.;
        let mut momentum_change = 0.;
        for (i, step) in out.steps.iter().enumerate() {
            let next = &step.state;
            close(next.motion.base.twist_body.linear.z, next_v);
            close(next.motion.base.twist_body.angular.norm(), 0.);
            let dz = next.motion.base.world_from_base.translation().z
                - before[i].motion.base.world_from_base.translation().z;
            close(dz, 0.5 * dt * (expected_v + next_v));
            energy_change +=
                0.5 * masses[i] * (next_v * next_v - expected_v * expected_v) + dry[i] * G * dz;
            momentum_change += masses[i] * (next_v - expected_v);
            close(
                next.memory[0].value,
                before[i].memory[0].advanced(HELD[0], dt).value,
            );
        }
        close(energy_change, 0.);
        close(momentum_change, -dt * G * (dry[0] + dry[1]));
        assert!(out.path_admission[0].iter().all(|g| *g >= 0.));
        assert!(out.initial_normal_speeds[0] >= -out.normal_speed_roundoff_mps[0]);
        assert_eq!(states, before);
        states = out.steps.into_iter().map(|s| s.state).collect();
        expected_v = next_v;
    }
}

#[test]
fn three_tree_stack_has_simultaneous_reactions_and_system_reindexing_invariance() {
    let dry = [1., 2., 3.];
    let spans = [0.25, 2., 6.];
    let masses: Vec<_> = dry
        .iter()
        .zip(spans)
        .map(|(m, s)| m + fluid_mass(s))
        .collect();
    let models = [model(dry[0]), model(dry[1]), model(dry[2])];
    let plates = spans.map(|s| strip(s, v(0., 0., 1.)));
    let states = [
        state(v(2., -3., 2. + 2. * GAP), Vec3::default()),
        state(v(2., -3., 1. + GAP), Vec3::default()),
        state(v(2., -3., 0.), Vec3::default()),
    ];
    let dt = 0.025;
    let a = -G * dry.iter().sum::<f64>() / masses.iter().sum::<f64>();
    let expected = [
        dt * (dry[0] * G + masses[0] * a),
        dt * ((dry[0] + dry[1]) * G + (masses[0] + masses[1]) * a),
    ];
    for order in [[0, 1, 2], [2, 0, 1]] {
        let index = |physical| order.iter().position(|i| *i == physical).unwrap();
        let inputs: Vec<_> = order
            .iter()
            .map(|&i| ContactSystem {
                model: &models[i],
                state: &states[i],
                impedance: &[],
                memory_inputs: &HELD,
                smooth: smooth(),
            })
            .collect();
        // Deliberately reverse identity order relative to the vertical chain.
        let contacts = [
            contact(91, index(0), index(1), 0.5, 0.),
            contact(7, index(1), index(2), 0.5, 0.),
        ];
        let out = step_with_contacts(
            &inputs,
            0.,
            dt,
            &contacts,
            FrictionalContactConfig::default(),
            |_, stages| {
                stages
                    .iter()
                    .enumerate()
                    .map(|(i, s)| loads(*s, &plates[order[i]], Vec3::default(), G))
                    .collect()
            },
            admit_affine,
            || false,
        )
        .unwrap();
        assert_eq!(out.contact_ids, [7, 91]);
        close(out.reaction_solve.impulses[0][0], expected[1]);
        close(out.reaction_solve.impulses[1][0], expected[0]);
        for physical in 0..3 {
            let next = &out.steps[index(physical)].state;
            close(next.motion.base.twist_body.linear.z, dt * a);
            close(
                next.motion.base.world_from_base.translation().z,
                states[physical].motion.base.world_from_base.translation().z + 0.5 * dt * dt * a,
            );
        }
        close(out.reaction_work_j, 0.);
        close(out.external_angular_impulse_world.norm(), 0.);
    }
}

#[test]
fn sliding_pair_uses_relative_full_mass_and_refuses_an_internal_finite_gap_couple() {
    let dry = [2., 4.];
    let spans = [1., 3.];
    let masses = [dry[0] + fluid_mass(spans[0]), dry[1] + fluid_mass(spans[1])];
    let models = [model(dry[0]), model(dry[1])];
    let plates = spans.map(|s| strip(s, v(1., 0., 0.)));
    let dt = 1. / 32.;
    let pn = 4. * dt;
    let pt = -0.25 * pn;
    let speeds = [3., -1.];
    let next = [speeds[0] + pt / masses[0], speeds[1] - pt / masses[1]];
    let middle = [0.5 * (speeds[0] + next[0]), 0.5 * (speeds[1] + next[1])];
    for gap in [0., GAP] {
        // Original material points coincide at the accepted midpoint in the
        // admissible case. The positive-gap twin has a real artificial couple.
        let states = [
            state(
                v(2. - 0.5 * dt * middle[0], -3., 1. + gap),
                v(speeds[0], 0., 0.),
            ),
            state(v(2. - 0.5 * dt * middle[1], -3., 1.), v(speeds[1], 0., 0.)),
        ];
        let result = step_with_contacts(
            &systems(&models, &states),
            0.,
            dt,
            &[contact(5, 0, 1, 0., 0.25)],
            FrictionalContactConfig::default(),
            |_, stages| {
                stages
                    .iter()
                    .enumerate()
                    .map(|(i, stage)| {
                        loads(
                            *stage,
                            &plates[i],
                            v(0., 0., if i == 0 { -4. } else { 4. }),
                            0.,
                        )
                    })
                    .collect()
            },
            admit_affine,
            || false,
        );
        if gap != 0. {
            let error = result.unwrap_err();
            assert_eq!(error.code, "articulated-contact-unresolved");
            assert!(
                error
                    .message
                    .contains("spatial impulse or finite-gap couple")
            );
            continue;
        }
        let out = result.unwrap();
        close(out.reaction_solve.impulses[0][0], pn);
        close(out.reaction_solve.impulses[0][1], pt);
        let mut kinetic_change = 0.;
        let mut momentum_change = 0.;
        for i in 0..2 {
            close(out.steps[i].state.motion.base.twist_body.linear.x, next[i]);
            close(out.steps[i].state.motion.base.twist_body.linear.z, 0.);
            kinetic_change += 0.5 * masses[i] * (next[i] * next[i] - speeds[i] * speeds[i]);
            momentum_change += masses[i] * (next[i] - speeds[i]);
        }
        close(momentum_change, 0.);
        close(kinetic_change, pt * (middle[0] - middle[1]));
        close(out.friction_work_j, kinetic_change);
        close(out.reaction_work_j, kinetic_change);
        close(out.normal_work_j, 0.);
        close(out.support_work_j, 0.);
        close(out.external_linear_impulse_world.norm(), 0.);
        close(out.external_angular_impulse_world.norm(), 0.);
    }
}

#[test]
fn empty_contacts_preserve_smooth_motion_and_synchronized_cross_tree_spring_loads() {
    let dry = [2., 4.];
    let spans = [1., 3.];
    let masses = [dry[0] + fluid_mass(spans[0]), dry[1] + fluid_mass(spans[1])];
    let models = [model(dry[0]), model(dry[1])];
    let plates = spans.map(|s| strip(s, v(1., 0., 0.)));
    let states = [
        state(v(2., 0., 0.), v(0.25, 0., 0.)),
        state(v(0., 0., 0.), v(-0.125, 0., 0.)),
    ];
    let dt = 1. / 32.;
    let h = 0.5 * dt;
    let r0 = 1.; // x0-x1 minus spring rest length1.
    let relative_speed = 0.375;
    for stiffness in [0., 32.] {
        let inverse_mass = 1. / masses[0] + 1. / masses[1];
        let middle_extension = (r0 + h * relative_speed) / (1. + h * h * stiffness * inverse_mass);
        let impulse = -dt * stiffness * middle_extension;
        let out = step_with_contacts(
            &systems(&models, &states),
            0.,
            dt,
            &[],
            FrictionalContactConfig::default(),
            |time, stages| {
                close(time, h);
                assert_eq!(stages.len(), 2);
                let extension = stages[0].state.motion.base.world_from_base.translation().x
                    - stages[1].state.motion.base.world_from_base.translation().x
                    - 1.;
                stages
                    .iter()
                    .enumerate()
                    .map(|(i, stage)| {
                        close(
                            stage.state.memory[0].value,
                            states[i].memory[0].advanced(HELD[0], h).value,
                        );
                        loads(
                            *stage,
                            &plates[i],
                            v(
                                if i == 0 {
                                    -stiffness * extension
                                } else {
                                    stiffness * extension
                                },
                                0.,
                                0.,
                            ),
                            0.,
                        )
                    })
                    .collect()
            },
            admit_affine,
            || false,
        )
        .unwrap();
        assert!(out.reaction_solve.impulses.is_empty());
        assert!(out.path_admission.is_empty());
        let mut kinetic_change = 0.;
        for i in 0..2 {
            let original = states[i].motion.base.twist_body.linear.x;
            let expected = original + if i == 0 { impulse } else { -impulse } / masses[i];
            let next = &out.steps[i].state;
            close(next.motion.base.twist_body.linear.x, expected);
            close(
                next.motion.base.world_from_base.translation().x,
                states[i].motion.base.world_from_base.translation().x + h * (original + expected),
            );
            kinetic_change += 0.5 * masses[i] * (expected * expected - original * original);
            if stiffness == 0. {
                let independent = step(
                    &models[i],
                    &[],
                    &states[i],
                    &HELD,
                    0.,
                    dt,
                    smooth(),
                    |_, state, kinematics| {
                        loads(
                            ArticulatedStage { state, kinematics },
                            &plates[i],
                            Vec3::default(),
                            0.,
                        )
                    },
                    || false,
                )
                .unwrap();
                assert_eq!(next, &independent.state);
            }
        }
        let extension_after = out.steps[0]
            .state
            .motion
            .base
            .world_from_base
            .translation()
            .x
            - out.steps[1]
                .state
                .motion
                .base
                .world_from_base
                .translation()
                .x
            - 1.;
        close(
            kinetic_change + 0.5 * stiffness * (extension_after * extension_after - r0 * r0),
            0.,
        );
    }
}

#[test]
fn resolved_impacts_strict_geometry_and_late_path_failure_never_publish_partial_states() {
    let models = [model(1.), model(2.)];
    let plates = [strip(0.5, v(0., 0., 1.)), strip(3., v(0., 0., 1.))];
    let states = [
        state(v(2., -3., 1. + GAP), Vec3::default()),
        state(v(2., -3., 0.), Vec3::default()),
    ];
    let saved = states.clone();
    let c = contact(1, 0, 1, 0.5, 0.);
    for cancel_at_end in [false, true] {
        let cancel = Cell::new(false);
        let observed = Cell::new(false);
        let result = step_with_contacts(
            &systems(&models, &states),
            0.,
            0.02,
            &[c],
            FrictionalContactConfig::default(),
            |_, stages| {
                stages
                    .iter()
                    .enumerate()
                    .map(|(i, stage)| loads(*stage, &plates[i], Vec3::default(), G))
                    .collect()
            },
            |path| {
                let admission = admit_affine(path)?;
                assert_eq!(path.endpoints.len(), 2);
                observed.set(true);
                if cancel_at_end {
                    cancel.set(true);
                    Ok(admission)
                } else {
                    Err(refusal("entire candidate rejected by geometry owner"))
                }
            },
            || cancel.get(),
        );
        assert!(observed.get());
        assert_eq!(
            result.unwrap_err().code,
            if cancel_at_end {
                "articulated-step-cancelled"
            } else {
                "path-falsifier"
            }
        );
        assert_eq!(states, saved);
    }
    let mut incoming = states.clone();
    incoming[0].motion.base.twist_body.linear.z = -1e-14;
    let mut loose = FrictionalContactConfig::default();
    loose.velocity_tolerance = 0.25;
    let err = step_with_contacts(
        &systems(&models, &incoming),
        0.,
        0.02,
        &[c],
        loose,
        |_, stages| {
            stages
                .iter()
                .enumerate()
                .map(|(i, stage)| loads(*stage, &plates[i], Vec3::default(), 0.))
                .collect()
        },
        admit_affine,
        || false,
    )
    .unwrap_err();
    assert_eq!(err.code, "articulated-contact-impact-required");
    let mut penetrating = states.clone();
    penetrating[0].motion.base.world_from_base =
        Se3::from_parts(So3::identity(), v(2., -3., 1. - GAP)).unwrap();
    let err = step_with_contacts(
        &systems(&models, &penetrating),
        0.,
        0.02,
        &[c],
        FrictionalContactConfig::default(),
        |_, _| panic!("invalid actual gap must refuse before loads"),
        admit_affine,
        || false,
    )
    .unwrap_err();
    assert_eq!(err.code, "articulated-contact-unresolved");
    assert_eq!(states, saved);
}
