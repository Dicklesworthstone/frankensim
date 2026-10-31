//! G1/G3/G4: complete fluid inertia in coupled normal/tangent impacts.
use fs_contact::frictional_multicontact::{
    FrictionalContactConfig, FrictionalContactDeterminacy, NormalVelocityLaw,
};
use fs_flyer::Refusal;
use fs_flyer::addedmass::articulated::{LinkedStrip, assemble_articulated_strips};
use fs_flyer::addedmass::{AeroGeneralizedLoads, Strip};
use fs_flyer::partitioned::MemoryState;
use fs_flyer::partitioned::articulated::ArticulatedState;
use fs_flyer::partitioned::articulated::impact::friction::{ArticulatedFrictionalContact, resolve};
use fs_flyer::partitioned::articulated::impact::{
    ArticulatedImpactSystem, ContactCounterpart, LinkContactPoint,
};
use fs_ga::{Mat3, Se3, So3, So3Tangent, Twist, Vec3};
use fs_mbd::articulated::generalized::GeneralizedKinematics;
use fs_mbd::articulated::motion::FreeFloatingMotionState;
use fs_mbd::articulated::{
    ArticulatedModel, FreeFloatingBaseState, JointLimits, JointModel, Link, SpatialInertia,
};
use fs_mbd::articulated_contacts::ArticulatedContactSystem;
use fs_mbd::articulated_contacts::friction::resolve_articulated_frictional_contacts;
use std::cell::Cell;

const RHO: f64 = 1.2;
fn v(x: f64, y: f64, z: f64) -> Vec3 {
    Vec3::new(x, y, z)
}
fn close(a: f64, b: f64) {
    assert!((a - b).abs() <= 3e-9 * (1.0 + b.abs()), "{a:e} != {b:e}");
}
fn vec_close(a: Vec3, b: Vec3) {
    close((a - b).norm(), 0.0);
}
fn inertia(mass: f64, j: f64) -> SpatialInertia {
    SpatialInertia::new(
        mass,
        Vec3::default(),
        Mat3 {
            m: [j, 0., 0., 0., j, 0., 0., 0., j],
        },
    )
    .unwrap()
}
fn body(mass: f64, j: f64) -> ArticulatedModel {
    ArticulatedModel::new(vec![Link::new(
        "body",
        None,
        Se3::identity(),
        JointModel::FIXED,
        inertia(mass, j),
    )])
    .unwrap()
}
fn hinge() -> ArticulatedModel {
    ArticulatedModel::new(vec![
        Link::new(
            "carrier",
            None,
            Se3::identity(),
            JointModel::FIXED,
            inertia(2., 5.),
        ),
        Link::new(
            "moving plate",
            Some(0),
            Se3::identity(),
            JointModel::revolute(
                v(0., 1., 0.),
                Some(JointLimits::new(-1., 1., 1000., 1e-8).unwrap()),
            )
            .unwrap(),
            inertia(3., 2.),
        ),
    ])
    .unwrap()
}
fn state(pose: Se3, velocity: Vec3, rates: Vec<f64>) -> ArticulatedState {
    ArticulatedState {
        motion: FreeFloatingMotionState {
            base: FreeFloatingBaseState::new(pose, Twist::new(Vec3::default(), velocity)),
            joint_position: vec![0.; rates.len()],
            joint_velocity: rates,
        },
        memory: vec![MemoryState {
            tau_s: 0.2,
            value: 0.7,
        }],
    }
}
fn fluid(
    rho: f64,
    link: usize,
    map: &GeneralizedKinematics,
) -> Result<AeroGeneralizedLoads, Refusal> {
    if rho == 0.0 {
        return assemble_articulated_strips(RHO, &[], map);
    }
    let strip = LinkedStrip::new(
        link,
        Strip {
            name: "plate",
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
    assemble_articulated_strips(rho, &[strip], map)
}
fn fluid_mass(rho: f64) -> f64 {
    rho * std::f64::consts::PI * 2. / 4.
}
fn fluid_pitch(rho: f64) -> f64 {
    rho * std::f64::consts::PI * 2. / 128.
}
fn point(system: usize, link: usize, p: Vec3) -> LinkContactPoint {
    LinkContactPoint {
        system_index: system,
        link_index: link,
        point_body: p,
    }
}
fn plane(
    id: u64,
    point: LinkContactPoint,
    pose: Se3,
    speed: Vec3,
    mu: f64,
) -> ArticulatedFrictionalContact {
    ArticulatedFrictionalContact {
        id,
        point,
        other: ContactCounterpart::PrescribedPlane {
            point_world: pose.transform_point(point.point_body).unwrap(),
            velocity_world: pose.rotation().rotate(speed).unwrap(),
        },
        normal_world: pose.rotation().rotate(v(0., 0., 1.)).unwrap(),
        tangent_world: pose.rotation().rotate(v(1., 0., 0.)).unwrap(),
        activation_tolerance: 0.,
        friction_coefficient: mu,
        normal_law: NormalVelocityLaw::NewtonRestitution { coefficient: 0. },
    }
}
fn energy(s: &ArticulatedState, m: f64, j: f64, rho: f64) -> f64 {
    let u = s.motion.base.twist_body;
    let a = fluid_mass(rho);
    let p = fluid_pitch(rho);
    0.5 * (m * (u.linear.x * u.linear.x + u.linear.y * u.linear.y)
        + (m + a) * u.linear.z * u.linear.z
        + (j + a / 3.) * u.angular.x * u.angular.x
        + (j + p) * u.angular.y * u.angular.y
        + j * u.angular.z * u.angular.z)
}

#[test]
fn oblique_hinge_impact_has_joint_recoil_and_normal_tangent_fluid_coupling() {
    let model = hinge();
    let p = v(0.6, 0., 0.2);
    let mu = 0.25;
    let total_mass = 5.;
    let mz = total_mass + fluid_mass(RHO);
    let jy = 2. + fluid_pitch(RHO);
    let wnn = 1. / mz + p.x * p.x / jy;
    let wnt = -p.x * p.z / jy;
    let wtt = 1. / total_mass + p.z * p.z / jy;
    let pn = 2. / (wnn - mu * wnt);
    let pt = -mu * pn;
    for pose in [
        Se3::identity(),
        Se3::from_parts(
            So3::exp(So3Tangent::new(v(0.2, -0.1, 0.3))).unwrap(),
            v(2., -1., 0.7),
        )
        .unwrap(),
    ] {
        let initial = state(pose, v(3.42, 0., -1.56), vec![0.4]);
        let saved = initial.clone();
        let c = plane(7, point(0, 1, p), pose, v(0.5, 0., 0.2), mu);
        let out = resolve(
            &[ArticulatedImpactSystem {
                model: &model,
                state: &initial,
            }],
            &[c],
            FrictionalContactConfig::default(),
            |_, map| fluid(RHO, 1, map),
            || false,
        )
        .unwrap();
        close(out.contact.impulses[0][0], pn);
        close(out.contact.impulses[0][1], pt);
        close(out.contact.impulses[0][2], 0.);
        let after = &out.mechanics.states[0];
        vec_close(after.motion.base.twist_body.angular, Vec3::default());
        close(
            after.motion.base.twist_body.linear.x,
            3.42 + pt / total_mass,
        );
        close(after.motion.base.twist_body.linear.z, -1.56 + pn / mz);
        close(
            after.motion.joint_velocity[0],
            0.4 + (p.z * pt - p.x * pn) / jy,
        );
        close(out.mechanics.relative_speeds_after[0], 0.);
        close(
            out.mechanics.relative_speeds_after[1],
            3. + (wnt - mu * wtt) * pn,
        );
        close(out.mechanics.support_work, 0.5 * pt + 0.2 * pn);
        let kinetic = |s: &ArticulatedState| {
            let u = s.motion.base.twist_body;
            0.5 * (total_mass * u.linear.x * u.linear.x
                + mz * u.linear.z * u.linear.z
                + 5. * u.angular.y * u.angular.y
                + jy * (u.angular.y + s.motion.joint_velocity[0]).powi(2))
        };
        close(out.mechanics.kinetic_before, kinetic(&initial));
        close(out.mechanics.kinetic_after, kinetic(after));
        close(
            kinetic(after) - kinetic(&initial) - out.mechanics.support_work,
            -out.contact.dissipated_energy,
        );
        vec_close(out.mechanics.linear_momentum_residual, Vec3::default());
        vec_close(out.mechanics.angular_momentum_residual, Vec3::default());
        assert_eq!(initial, saved);
        assert_eq!(after.memory, initial.memory);
        assert_eq!(after.motion.joint_position, initial.motion.joint_position);
        assert_eq!(after.motion.base.world_from_base, pose);
        assert!((pn - 2. / wnn).abs() > 0.1); // sequential normal then tangent would be wrong
    }
}

#[test]
fn separate_bodies_exchange_frictional_momentum_with_each_fluids_inertia() {
    let a = body(2., 0.8);
    let b = body(3., 1.3);
    let poses = [
        Se3::from_parts(So3::identity(), v(0., 0., 1.)).unwrap(),
        Se3::from_parts(So3::identity(), v(0., 0., -1.)).unwrap(),
    ];
    let states = [
        state(poses[0], v(1., 0., -1.), vec![]),
        state(poses[1], Vec3::default(), vec![]),
    ];
    let systems = [
        ArticulatedImpactSystem {
            model: &a,
            state: &states[0],
        },
        ArticulatedImpactSystem {
            model: &b,
            state: &states[1],
        },
    ];
    let c = ArticulatedFrictionalContact {
        other: ContactCounterpart::Body(point(1, 0, v(0., 0., 1.))),
        ..plane(
            13,
            point(0, 0, v(0., 0., -1.)),
            poses[0],
            Vec3::default(),
            0.1,
        )
    };
    for rho in [0., RHO] {
        let out = resolve(
            &systems,
            &[c],
            FrictionalContactConfig::default(),
            |_, map| fluid(rho, 0, map),
            || false,
        )
        .unwrap();
        let mz = [2. + fluid_mass(rho), 3. + fluid_mass(rho)];
        let jy = [0.8 + fluid_pitch(rho), 1.3 + fluid_pitch(rho)];
        let pn = 1. / (1. / mz[0] + 1. / mz[1]);
        let pt = -0.1 * pn;
        close(out.contact.impulses[0][0], pn);
        close(out.contact.impulses[0][1], pt);
        let u = out.mechanics.states[0].motion.base.twist_body;
        let w = out.mechanics.states[1].motion.base.twist_body;
        vec_close(u.linear, v(1. + pt / 2., 0., -1. + pn / mz[0]));
        vec_close(w.linear, v(-pt / 3., 0., -pn / mz[1]));
        vec_close(u.angular, v(0., -pt / jy[0], 0.));
        vec_close(w.angular, v(0., -pt / jy[1], 0.));
        let before = energy(&states[0], 2., 0.8, rho) + energy(&states[1], 3., 1.3, rho);
        let after = energy(&out.mechanics.states[0], 2., 0.8, rho)
            + energy(&out.mechanics.states[1], 3., 1.3, rho);
        close(out.mechanics.kinetic_before, before);
        close(out.mechanics.kinetic_after, after);
        close(before - after, out.contact.dissipated_energy);
        close(out.mechanics.support_work, 0.);
        vec_close(out.mechanics.linear_momentum_residual, Vec3::default());
        vec_close(out.mechanics.angular_momentum_residual, Vec3::default());
        // Geometric activation alone cannot authorize a tangential couple
        // between two separated material points: its net angular impulse has
        // no external source. Keep that physical refusal at finite gap.
        let mut separated = states[1].clone();
        separated.motion.base.world_from_base =
            Se3::from_parts(So3::identity(), v(0., 0., -1.001)).unwrap();
        let separated_systems = [
            systems[0],
            ArticulatedImpactSystem {
                model: &b,
                state: &separated,
            },
        ];
        let separated_contact = ArticulatedFrictionalContact {
            activation_tolerance: 0.01,
            ..c
        };
        let refused = resolve(
            &separated_systems,
            &[separated_contact],
            FrictionalContactConfig::default(),
            |_, map| fluid(rho, 0, map),
            || false,
        )
        .unwrap_err();
        assert_eq!(refused.code, "articulated-impact-unresolved");
        assert!(refused.message.contains("spatial impulse"));
        if rho == 0. {
            let native =
                [(&a, &states[0]), (&b, &states[1])].map(|(model, s)| ArticulatedContactSystem {
                    model,
                    base: s.motion.base,
                    joint_position: &s.motion.joint_position,
                    joint_velocity: &s.motion.joint_velocity,
                });
            let dry = resolve_articulated_frictional_contacts(
                &native,
                &[c],
                FrictionalContactConfig::default(),
                || false,
            )
            .unwrap();
            for (got, expected) in out.mechanics.states.iter().zip(dry.responses) {
                vec_close(
                    got.motion.base.twist_body.linear,
                    expected.base_after.twist_body.linear,
                );
                vec_close(
                    got.motion.base.twist_body.angular,
                    expected.base_after.twist_body.angular,
                );
            }
        }
    }
}

#[test]
fn redundant_four_foot_contact_determines_fluid_loaded_motion_without_inventing_local_loads() {
    let model = body(3., 2.);
    let pose = Se3::from_parts(So3::identity(), v(0., 0., 1.)).unwrap();
    let initial = state(pose, v(3., 0., -1.), vec![]);
    let contacts = [(31, 1., 1.), (4, -1., 1.), (17, -1., -1.), (9, 1., -1.)]
        .map(|(id, x, y)| plane(id, point(0, 0, v(x, y, -1.)), pose, Vec3::default(), 0.1));
    let systems = [ArticulatedImpactSystem {
        model: &model,
        state: &initial,
    }];
    let out = resolve(
        &systems,
        &contacts,
        FrictionalContactConfig::default(),
        |_, map| fluid(RHO, 0, map),
        || false,
    )
    .unwrap();
    let FrictionalContactDeterminacy::UniqueVelocity(witness) = &out.contact.determinacy else {
        panic!("four-foot pressure gauge expected")
    };
    assert_eq!(witness.dimension, 12);
    assert_eq!(witness.rank, 6);
    let pn = 3. + fluid_mass(RHO);
    let vx = 3. - 0.1 * pn / 3.;
    let after = &out.mechanics.states[0];
    vec_close(after.motion.base.twist_body.linear, v(vx, 0., 0.));
    vec_close(after.motion.base.twist_body.angular, Vec3::default());
    close(out.contact.impulses.iter().map(|p| p[0]).sum(), pn);
    close(out.contact.impulses.iter().map(|p| p[1]).sum(), -0.1 * pn);
    close(
        out.contact.dissipated_energy,
        energy(&initial, 3., 2., RHO) - energy(after, 3., 2., RHO),
    );
    let mut reversed = contacts;
    reversed.reverse();
    let same = resolve(
        &systems,
        &reversed,
        FrictionalContactConfig::default(),
        |_, map| fluid(RHO, 0, map),
        || false,
    )
    .unwrap();
    assert_eq!(out.contact_ids, vec![4, 9, 17, 31]);
    assert_eq!(out.contact.impulses, same.contact.impulses);
    assert_eq!(out.mechanics.states, same.mechanics.states);
}

#[test]
fn ambiguous_support_work_duplicate_ids_and_late_cancellation_refuse_atomically() {
    let model = body(3., 2.);
    let initial = state(Se3::identity(), v(3., 0., -1.), vec![]);
    let saved = initial.clone();
    let systems = [ArticulatedImpactSystem {
        model: &model,
        state: &initial,
    }];
    let c = plane(
        1,
        point(0, 0, Vec3::default()),
        Se3::identity(),
        Vec3::default(),
        0.1,
    );
    let other = ArticulatedFrictionalContact {
        id: 2,
        other: ContactCounterpart::PrescribedPlane {
            point_world: Vec3::default(),
            velocity_world: v(0.01, 0., 0.),
        },
        ..c
    };
    assert_eq!(
        resolve(
            &systems,
            &[c, other],
            FrictionalContactConfig::default(),
            |_, map| fluid(RHO, 0, map),
            || false
        )
        .unwrap_err()
        .code,
        "articulated-friction-support-work"
    );
    assert_eq!(
        resolve(
            &systems,
            &[c, c],
            FrictionalContactConfig::default(),
            |_, map| fluid(RHO, 0, map),
            || false
        )
        .unwrap_err()
        .code,
        "articulated-friction-invalid"
    );
    let probes = Cell::new(0);
    resolve(
        &systems,
        &[c],
        FrictionalContactConfig::default(),
        |_, map| fluid(RHO, 0, map),
        || {
            probes.set(probes.get() + 1);
            false
        },
    )
    .unwrap();
    let mut seen = 0;
    assert!(
        resolve(
            &systems,
            &[c],
            FrictionalContactConfig::default(),
            |_, map| fluid(RHO, 0, map),
            || {
                seen += 1;
                seen == probes.get()
            }
        )
        .is_err()
    );
    assert_eq!(initial, saved);
}
