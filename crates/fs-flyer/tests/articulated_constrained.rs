//! Independent rigid-link reductions of simultaneous impacts with attached fluid.
//! The oracles use physical impulse and spatial-momentum balance, not the
//! production Schur factor or its projected mass matrix.

use fs_contact::frictional_multicontact::{
    FrictionalContactConfig, FrictionalContactDeterminacy, NormalVelocityLaw,
};
use fs_flyer::Refusal;
use fs_flyer::addedmass::articulated::{LinkedStrip, assemble_articulated_strips};
use fs_flyer::addedmass::{AeroGeneralizedLoads, Strip};
use fs_flyer::partitioned::MemoryState;
use fs_flyer::partitioned::articulated::ArticulatedState;
use fs_flyer::partitioned::articulated::impact::friction::ArticulatedFrictionalContact;
use fs_flyer::partitioned::articulated::impact::friction::constrained::{
    ConstrainedArticulatedImpact, ConstrainedImpactConfig, ImpactBilateralConstraint,
    ImpactDistanceConstraint, ImpactPointConstraint, ImpactWeldConstraint, resolve,
};
use fs_flyer::partitioned::articulated::impact::{
    ArticulatedImpactSystem, ContactCounterpart, LinkContactPoint,
};
use fs_ga::{Mat3, Se3, So3, So3Tangent, Twist, Vec3};
use fs_mbd::articulated::generalized::GeneralizedKinematics;
use fs_mbd::articulated::motion::FreeFloatingMotionState;
use fs_mbd::articulated::{
    ArticulatedModel, FreeFloatingBaseState, JointModel, Link, SpatialInertia,
};
use fs_mbd::articulated_contacts::ArticulatedContactSystem;
use fs_mbd::articulated_contacts::friction::constrained::bilateral::resolve_bilateral_articulated_frictional_contacts;
use std::cell::Cell;

#[derive(Clone, Copy)]
struct Body {
    mass: f64,
    rotary: f64,
    density: f64,
}
impl Body {
    // Unit-chord, two-metre-span strip, normal z and span y, centred at the COM.
    fn fluid_mass(self) -> f64 {
        self.density * std::f64::consts::PI / 2.
    }
    fn normal_mass(self) -> f64 {
        self.mass + self.fluid_mass()
    }
    fn pitch_inertia(self) -> f64 {
        self.rotary + self.fluid_mass() / 32.
    }
    fn model(self) -> ArticulatedModel {
        ArticulatedModel::new(vec![Link::new(
            "body",
            None,
            Se3::identity(),
            JointModel::FIXED,
            SpatialInertia::new(
                self.mass,
                Vec3::default(),
                Mat3 {
                    m: [
                        self.rotary,
                        0.,
                        0.,
                        0.,
                        self.rotary,
                        0.,
                        0.,
                        0.,
                        self.rotary,
                    ],
                },
            )
            .unwrap(),
        )])
        .unwrap()
    }
    fn fluid(self, map: &GeneralizedKinematics) -> Result<AeroGeneralizedLoads, Refusal> {
        if self.density == 0. {
            return assemble_articulated_strips(1., &[], map);
        }
        let strip = LinkedStrip::new(
            0,
            Strip {
                name: "physical strip",
                chord_m: 1.,
                span_m: 2.,
                position_m: [0.; 3],
                normal: [0., 0., 1.],
                span_axis: [0., 1., 0.],
                control_coord: None,
                control_gain: 0.,
                control_angular_gain: [0.; 3],
            },
        )?;
        assemble_articulated_strips(self.density, &[strip], map)
    }
}
fn v(x: f64, y: f64, z: f64) -> Vec3 {
    Vec3::new(x, y, z)
}
fn close(a: f64, b: f64) {
    assert!((a - b).abs() <= 4e-9 * (1. + b.abs()), "{a:e} != {b:e}");
}
fn vector(a: Vec3, b: Vec3) {
    close((a - b).norm(), 0.);
}
fn frames() -> [Se3; 2] {
    [
        Se3::identity(),
        Se3::from_parts(
            So3::exp(So3Tangent::new(v(0.31, -0.23, 0.17))).unwrap(),
            v(1.4, -0.7, 2.3),
        )
        .unwrap(),
    ]
}
fn state(frame: Se3, position: Vec3, velocity: Vec3) -> ArticulatedState {
    let pose = frame
        .compose(Se3::from_parts(So3::identity(), position).unwrap())
        .unwrap();
    ArticulatedState {
        motion: FreeFloatingMotionState {
            base: FreeFloatingBaseState::new(pose, Twist::new(Vec3::default(), velocity)),
            joint_position: Vec::new(),
            joint_velocity: Vec::new(),
        },
        memory: vec![MemoryState {
            tau_s: 0.3,
            value: 0.7,
        }],
    }
}
fn point(system_index: usize, point_body: Vec3) -> LinkContactPoint {
    LinkContactPoint {
        system_index,
        link_index: 0,
        point_body,
    }
}
fn rod(a: usize, b: usize) -> ImpactBilateralConstraint {
    ImpactBilateralConstraint::Distance(ImpactDistanceConstraint {
        a: point(a, Vec3::default()),
        b: point(b, Vec3::default()),
        length_m: 1.,
    })
}
fn plane(
    id: u64,
    system: usize,
    s: &ArticulatedState,
    mu: f64,
    e: f64,
) -> ArticulatedFrictionalContact {
    let pose = s.motion.base.world_from_base;
    ArticulatedFrictionalContact {
        id,
        point: point(system, Vec3::default()),
        other: ContactCounterpart::PrescribedPlane {
            point_world: pose.translation(),
            velocity_world: Vec3::default(),
        },
        normal_world: pose.rotation().rotate(v(0., 0., 1.)).unwrap(),
        tangent_world: pose.rotation().rotate(v(1., 0., 0.)).unwrap(),
        // The rotated native material point differs from pose.translation() by
        // a few ulps; this explicit activation band retains the same plane.
        activation_tolerance: 1e-12,
        friction_coefficient: mu,
        normal_law: NormalVelocityLaw::NewtonRestitution { coefficient: e },
    }
}
fn systems<'a>(
    models: &'a [ArticulatedModel],
    states: &'a [ArticulatedState],
) -> Vec<ArticulatedImpactSystem<'a>> {
    models
        .iter()
        .zip(states)
        .map(|(model, state)| ArticulatedImpactSystem { model, state })
        .collect()
}
fn run(
    bodies: &[Body],
    states: &[ArticulatedState],
    contacts: &[ArticulatedFrictionalContact],
    constraints: &[ImpactBilateralConstraint],
    cfg: ConstrainedImpactConfig,
) -> ConstrainedArticulatedImpact {
    let models: Vec<_> = bodies.iter().map(|b| b.model()).collect();
    resolve(
        &systems(&models, states),
        contacts,
        constraints,
        FrictionalContactConfig::default(),
        cfg,
        |i, map| bodies[i].fluid(map),
        || false,
    )
    .unwrap()
}

// Independent strip spatial moments: normal translation, span-distributed
// normal rotation about x, and section pitch about y. No assembled matrix is read.
fn spatial(bodies: &[Body], states: &[ArticulatedState]) -> (Vec3, Vec3, f64) {
    let (mut linear, mut angular, mut kinetic) = (Vec3::default(), Vec3::default(), 0.);
    for (&body, s) in bodies.iter().zip(states) {
        let pose = s.motion.base.world_from_base;
        let u = s.motion.base.twist_body;
        let p = v(
            body.mass * u.linear.x,
            body.mass * u.linear.y,
            body.normal_mass() * u.linear.z,
        );
        let l = v(
            (body.rotary + body.fluid_mass() / 3.) * u.angular.x,
            body.pitch_inertia() * u.angular.y,
            body.rotary * u.angular.z,
        );
        let world_p = pose.rotation().rotate(p).unwrap();
        linear = linear + world_p;
        angular = angular + pose.rotation().rotate(l).unwrap() + pose.translation().cross(world_p);
        kinetic += 0.5 * (u.linear.dot(p) + u.angular.dot(l));
    }
    (linear, angular, kinetic)
}
fn audit(
    bodies: &[Body],
    before: &[ArticulatedState],
    contacts: &[ArticulatedFrictionalContact],
    result: &ConstrainedArticulatedImpact,
) {
    let m = &result.impact.mechanics;
    let (p0, l0, t0) = spatial(bodies, before);
    let (p1, l1, t1) = spatial(bodies, &m.states);
    let (mut external_p, mut external_l) = (Vec3::default(), Vec3::default());
    for (&id, impulse) in result
        .impact
        .contact_ids
        .iter()
        .zip(&result.impact.contact.impulses)
    {
        let c = contacts.iter().find(|c| c.id == id).unwrap();
        let force = c.normal_world.scale(impulse[0])
            + c.tangent_world.scale(impulse[1])
            + c.normal_world.cross(c.tangent_world).scale(impulse[2]);
        let pos = before[c.point.system_index]
            .motion
            .base
            .world_from_base
            .transform_point(c.point.point_body)
            .unwrap();
        external_p = external_p + force;
        external_l = external_l + pos.cross(force);
    }
    vector(p1 - p0, external_p);
    vector(l1 - l0, external_l);
    close(t0, m.kinetic_before);
    close(t1, m.kinetic_after);
    close(t1 - t0, m.kinetic_change);
    close(
        t1 - t0 + result.impact.contact.dissipated_energy,
        result.constraint_work_j,
    );
    close(m.support_work, 0.);
    close(m.energy_balance_residual, 0.);
    vector(m.linear_momentum_residual, Vec3::default());
    vector(m.angular_momentum_residual, Vec3::default());
    close(m.endpoint_velocity_residual, 0.);
    assert!(result.velocity_after.iter().all(|x| x.abs() <= 1e-9));
    for (old, new) in before.iter().zip(&m.states) {
        assert_eq!(
            old.motion.base.world_from_base,
            new.motion.base.world_from_base
        );
        assert_eq!(old.motion.joint_position, new.motion.joint_position);
        assert_eq!(old.memory, new.memory);
    }
}
fn compare_native(
    bodies: &[Body],
    states: &[ArticulatedState],
    contacts: &[ArticulatedFrictionalContact],
    constraints: &[ImpactBilateralConstraint],
    result: &ConstrainedArticulatedImpact,
) {
    assert!(bodies.iter().all(|b| b.density == 0.));
    let models: Vec<_> = bodies.iter().map(|b| b.model()).collect();
    let native: Vec<_> = models
        .iter()
        .zip(states)
        .map(|(model, state)| ArticulatedContactSystem {
            model,
            base: state.motion.base,
            joint_position: &state.motion.joint_position,
            joint_velocity: &state.motion.joint_velocity,
        })
        .collect();
    let expected = resolve_bilateral_articulated_frictional_contacts(
        &native,
        contacts,
        constraints,
        FrictionalContactConfig::default(),
        ConstrainedImpactConfig::default(),
        || false,
    )
    .unwrap();
    for (a, b) in expected.impact.delassus.iter().zip(&result.impact.delassus) {
        close(*a, *b);
    }
    for (a, b) in expected
        .impact
        .contact
        .impulses
        .iter()
        .flatten()
        .zip(result.impact.contact.impulses.iter().flatten())
    {
        close(*a, *b);
    }
    for (a, b) in expected
        .bilateral_impulses_ns
        .iter()
        .zip(&result.bilateral_impulses)
    {
        close(*a, *b);
    }
    for (a, b) in expected
        .impact
        .responses
        .iter()
        .zip(&result.impact.mechanics.states)
    {
        vector(
            a.base_after.twist_body.linear,
            b.motion.base.twist_body.linear,
        );
        vector(
            a.base_after.twist_body.angular,
            b.motion.base.twist_body.angular,
        );
    }
    close(
        expected.impact.contact.dissipated_energy,
        result.impact.contact.dissipated_energy,
    );
}

#[test]
fn rigid_distance_transmits_full_fluid_normal_impulse_without_locking_tangential_motion() {
    for density in [0., 1.2] {
        let bodies = [
            Body {
                mass: 2.,
                rotary: 0.7,
                density,
            },
            Body {
                mass: 3.,
                rotary: 1.1,
                density: 0.6 * density,
            },
        ];
        for frame in frames() {
            let states = [
                state(frame, Vec3::default(), v(2., 0., -1.)),
                state(frame, v(0., 0., 1.), v(2., 0., -1.)),
            ];
            for e in [0., 0.4] {
                let contacts = [plane(7, 0, &states[0], 0.2, e)];
                let constraints = [rod(0, 1)];
                let result = run(
                    &bodies,
                    &states,
                    &contacts,
                    &constraints,
                    ConstrainedImpactConfig::default(),
                );
                let mass = bodies[0].normal_mass() + bodies[1].normal_mass();
                let pn = mass * (1. + e);
                let pt = -0.2 * pn;
                close(result.impact.contact.impulses[0][0], pn);
                close(result.impact.contact.impulses[0][1], pt);
                close(result.impact.delassus[0], 1. / mass);
                close(result.impact.delassus[4], 1. / bodies[0].mass);
                close(
                    result.bilateral_impulses[0],
                    bodies[1].normal_mass() * (1. + e),
                );
                vector(
                    result.impact.mechanics.states[0]
                        .motion
                        .base
                        .twist_body
                        .linear,
                    v(2. + pt / bodies[0].mass, 0., e),
                );
                vector(
                    result.impact.mechanics.states[1]
                        .motion
                        .base
                        .twist_body
                        .linear,
                    v(2., 0., e),
                );
                close(
                    result.impact.contact.dissipated_energy,
                    0.5 * mass * (1. - e * e) - 2. * pt - pt * pt / (2. * bodies[0].mass),
                );
                close(result.constraint_work_j, 0.);
                audit(&bodies, &states, &contacts, &result);
                if density == 0. {
                    compare_native(&bodies, &states, &contacts, &constraints, &result);
                }
            }
        }
    }
}

#[test]
fn point_and_weld_impacts_match_distinct_physical_rotor_reductions() {
    let (a, b) = (0.6, 0.9);
    for density in [0., 1.2] {
        let bodies = [
            Body {
                mass: 2.,
                rotary: 0.7,
                density,
            },
            Body {
                mass: 3.,
                rotary: 1.1,
                density: 0.6 * density,
            },
        ];
        let (ma, mb, ja, jb) = (
            bodies[0].normal_mass(),
            bodies[1].normal_mass(),
            bodies[0].pitch_inertia(),
            bodies[1].pitch_inertia(),
        );
        for frame in frames() {
            let states = [
                state(frame, v(-a, 0., 0.), v(0., 0., -1.)),
                state(frame, v(b, 0., 0.), v(0., 0., -1.)),
            ];
            let contacts = [plane(9, 0, &states[0], 0., 0.)];
            for weld in [false, true] {
                let pa = point(0, v(a, 0., 0.));
                let pb = point(1, v(-b, 0., 0.));
                let constraint = if weld {
                    ImpactBilateralConstraint::Weld(ImpactWeldConstraint { a: pa, b: pb })
                } else {
                    ImpactBilateralConstraint::PointCoincidence(ImpactPointConstraint {
                        a: pa,
                        b: pb,
                    })
                };
                let result = run(
                    &bodies,
                    &states,
                    &contacts,
                    &[constraint],
                    ConstrainedImpactConfig::default(),
                );
                // Linear/angular impulse balance plus coincidence (and equal
                // angular rate for the weld) gives these scalar closed forms.
                let (lambda, wa, wb, vb, pure_couple) = if weld {
                    let w = mb * (a + b) / (ja + jb + mb * (a + b) * (a + b));
                    let lambda = mb * (1. - (a + b) * w);
                    (lambda, w, w, -(a + b) * w, a * lambda - ja * w)
                } else {
                    let lambda = 1. / (1. / mb + a * a / ja + b * b / jb);
                    (
                        lambda,
                        a * lambda / ja,
                        b * lambda / jb,
                        -1. + lambda / mb,
                        0.,
                    )
                };
                let pn = ma + lambda;
                close(result.impact.contact.impulses[0][0], pn);
                close(result.impact.delassus[0], 1. / pn);
                close(result.bilateral_impulses[2], lambda);
                vector(
                    result.impact.mechanics.states[0]
                        .motion
                        .base
                        .twist_body
                        .linear,
                    Vec3::default(),
                );
                vector(
                    result.impact.mechanics.states[1]
                        .motion
                        .base
                        .twist_body
                        .linear,
                    v(0., 0., vb),
                );
                vector(
                    result.impact.mechanics.states[0]
                        .motion
                        .base
                        .twist_body
                        .angular,
                    v(0., wa, 0.),
                );
                vector(
                    result.impact.mechanics.states[1]
                        .motion
                        .base
                        .twist_body
                        .angular,
                    v(0., wb, 0.),
                );
                if weld {
                    close(result.bilateral_impulses[4], pure_couple);
                }
                close(result.impact.contact.dissipated_energy, 0.5 * pn);
                audit(&bodies, &states, &contacts, &result);
                if density == 0. {
                    compare_native(&bodies, &states, &contacts, &[constraint], &result);
                }
            }
        }
    }
}

#[test]
fn toleranced_point_gap_retains_its_geometric_couple_and_world_angular_momentum() {
    let bodies = [
        Body {
            mass: 2.,
            rotary: 0.7,
            density: 1.2,
        },
        Body {
            mass: 3.,
            rotary: 1.1,
            density: 0.8,
        },
    ];
    let gap = 5e-5;
    let cfg = ConstrainedImpactConfig {
        position_tolerance_m: 1e-4,
        ..Default::default()
    };
    let lambda = 1. / (1. / bodies[1].normal_mass() + gap * gap / bodies[0].pitch_inertia());
    let spin = gap * lambda / bodies[0].pitch_inertia();
    for frame in frames() {
        let states = [
            state(frame, Vec3::default(), v(0., 0., -1.)),
            state(frame, v(gap, 0., 0.), v(0., 0., -1.)),
        ];
        let contacts = [plane(3, 0, &states[0], 0., 0.)];
        let constraints = [ImpactBilateralConstraint::PointCoincidence(
            ImpactPointConstraint {
                a: point(0, Vec3::default()),
                b: point(1, Vec3::default()),
            },
        )];
        let result = run(&bodies, &states, &contacts, &constraints, cfg);
        close(result.position_residual_m[0], gap);
        close(result.bilateral_impulses[2], lambda);
        close(
            result.impact.contact.impulses[0][0],
            bodies[0].normal_mass() + lambda,
        );
        vector(
            result.impact.mechanics.states[0]
                .motion
                .base
                .twist_body
                .angular,
            v(0., spin, 0.),
        );
        vector(
            result.impact.mechanics.states[1]
                .motion
                .base
                .twist_body
                .angular,
            Vec3::default(),
        );
        close(
            result.impact.mechanics.states[1]
                .motion
                .base
                .twist_body
                .linear
                .z,
            -gap * spin,
        );
        // Merely applying equal/opposite endpoint forces would leave a torque
        // gap*lambda; the A-side geometric couple must cancel it exactly.
        assert!(gap * lambda > 1e-4);
        audit(&bodies, &states, &contacts, &result);
    }
}

#[test]
fn linked_chain_with_redundant_contacts_has_determined_motion_under_all_reindexings() {
    let bodies = [
        Body {
            mass: 1.,
            rotary: 0.7,
            density: 1.2,
        },
        Body {
            mass: 2.,
            rotary: 1.1,
            density: 0.7,
        },
        Body {
            mass: 3.,
            rotary: 1.3,
            density: 1.8,
        },
    ];
    let mass: f64 = bodies.iter().map(|b| b.normal_mass()).sum();
    for reverse_systems in [false, true] {
        let index = |i| if reverse_systems { 2 - i } else { i };
        let ordered: Vec<_> = (0..3).map(|i| bodies[index(i)]).collect();
        let states: Vec<_> = (0..3)
            .map(|i| state(Se3::identity(), v(0., 0., index(i) as f64), v(0., 0., -1.)))
            .collect();
        for reverse_rows in [false, true] {
            let mut constraints = [rod(index(0), index(1)), rod(index(1), index(2))];
            let mut contacts = [
                plane(11, index(0), &states[index(0)], 0., 0.),
                plane(5, index(2), &states[index(2)], 0., 0.),
            ];
            if reverse_rows {
                constraints.reverse();
                contacts.reverse();
            }
            let result = run(
                &ordered,
                &states,
                &contacts,
                &constraints,
                ConstrainedImpactConfig::default(),
            );
            assert_eq!(result.impact.contact_ids, vec![5, 11]);
            assert!(matches!(
                result.impact.contact.determinacy,
                FrictionalContactDeterminacy::UniqueVelocity(_)
            ));
            close(
                result.impact.contact.impulses.iter().map(|p| p[0]).sum(),
                mass,
            );
            close(result.impact.contact.dissipated_energy, 0.5 * mass);
            for s in &result.impact.mechanics.states {
                vector(s.motion.base.twist_body.linear, Vec3::default());
                vector(s.motion.base.twist_body.angular, Vec3::default());
            }
            audit(&ordered, &states, &contacts, &result);
        }
    }
}

#[test]
fn inconsistent_or_dependent_closures_and_late_failures_publish_no_state() {
    let bodies = [
        Body {
            mass: 2.,
            rotary: 0.7,
            density: 1.2,
        },
        Body {
            mass: 3.,
            rotary: 1.1,
            density: 0.8,
        },
    ];
    let models: Vec<_> = bodies.iter().map(|b| b.model()).collect();
    let states = [
        state(Se3::identity(), Vec3::default(), v(0., 0., -1.)),
        state(Se3::identity(), v(0., 0., 1.), v(0., 0., -1.)),
    ];
    let saved = states.clone();
    let contacts = [plane(1, 0, &states[0], 0., 0.)];
    let constraints = [rod(0, 1)];
    let cfg = ConstrainedImpactConfig::default();
    let run_bad = |s: &[ArticulatedState], c: &[ImpactBilateralConstraint]| {
        resolve(
            &systems(&models, s),
            &contacts,
            c,
            FrictionalContactConfig::default(),
            cfg,
            |i, map| bodies[i].fluid(map),
            || false,
        )
    };
    assert!(run_bad(&states, &[constraints[0], constraints[0]]).is_err());
    let wrong_length = ImpactBilateralConstraint::Distance(ImpactDistanceConstraint {
        a: point(0, Vec3::default()),
        b: point(1, Vec3::default()),
        length_m: 1.1,
    });
    assert!(run_bad(&states, &[wrong_length]).is_err());
    let mut bad_rate = states.clone();
    bad_rate[1].motion.base.twist_body.linear.z += 0.1;
    let retained_bad_rate = bad_rate.clone();
    assert!(run_bad(&bad_rate, &constraints).is_err());
    assert_eq!(bad_rate, retained_bad_rate);
    let calls = Cell::new(0);
    let error = resolve(
        &systems(&models, &states),
        &contacts,
        &constraints,
        FrictionalContactConfig::default(),
        cfg,
        |i, map| {
            calls.set(calls.get() + 1);
            if i == 1 {
                Err(Refusal {
                    code: "late-fluid",
                    message: "injected second-tree refusal".into(),
                    ranked_repairs: Vec::new(),
                })
            } else {
                bodies[i].fluid(map)
            }
        },
        || false,
    )
    .unwrap_err();
    assert_eq!(error.code, "late-fluid");
    assert_eq!(calls.get(), 2);
    let polls = Cell::new(0);
    resolve(
        &systems(&models, &states),
        &contacts,
        &constraints,
        FrictionalContactConfig::default(),
        cfg,
        |i, map| bodies[i].fluid(map),
        || {
            polls.set(polls.get() + 1);
            false
        },
    )
    .unwrap();
    let final_poll = polls.get();
    assert!(final_poll > 1);
    polls.set(0);
    let error = resolve(
        &systems(&models, &states),
        &contacts,
        &constraints,
        FrictionalContactConfig::default(),
        cfg,
        |i, map| bodies[i].fluid(map),
        || {
            polls.set(polls.get() + 1);
            polls.get() == final_poll
        },
    )
    .unwrap_err();
    assert!(error.code.contains("cancel") || error.message.contains("cancel"));
    assert_eq!(states, saved);
}
