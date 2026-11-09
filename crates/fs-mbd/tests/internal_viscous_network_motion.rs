//! G1/G3 moving passive graphs: continuous modes, material frames and actual heat.

#[path = "internal_viscous_network_motion/driven.rs"]
mod driven;

#[path = "internal_viscous_network_motion/relaxation.rs"]
mod relaxation;

use fs_mbd::spatial_loads::internal_impulse::{
    InternalWrenchImpulseControls,
    viscous::{
        BushingCoordinateScales, PassiveBushingMatrix, ViscousConnectorControls,
        network::{ViscousNetworkLimits, ViscousNetworkLink},
    },
};
use fs_mbd::spatial_loads::internal_motion::{
    HeldExternalLoads, HeldSpatialLoads,
    viscous::{
        ViscousMotionLaw,
        network::{
            ViscousNetworkMotionControls, ViscousNetworkMotionError, step_viscous_network,
            step_viscous_network_with_loads, step_viscous_network_with_spatial_loads,
        },
        step_viscous_connector, step_viscous_connector_with_spatial_loads,
    },
};
use fs_mbd::spatial_loads::{
    PreparedSpatialLoads, WorldPointForce, midpoint::SpatialMidpointConfig,
};
use fs_mbd::{
    Gravity, MassProperties, Pose, RigidBodyIntegrator, RigidBodyState, UnitQuaternion, Vec3,
    Wrench,
};

const X: Vec3 = Vec3::new(1., 0., 0.);
fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:e} != {b:e}, tolerance {tolerance:e}"
    );
}
fn vector(a: Vec3, b: Vec3, tolerance: f64) {
    near(a.x, b.x, tolerance);
    near(a.y, b.y, tolerance);
    near(a.z, b.z, tolerance);
}
fn norm(a: Vec3) -> f64 {
    a.x.hypot(a.y).hypot(a.z)
}
fn controls() -> ViscousNetworkMotionControls {
    ViscousNetworkMotionControls {
        source: ViscousConnectorControls {
            impulse: InternalWrenchImpulseControls {
                angular_origin_world_m: Vec3::new(0.2, -0.3, 0.1),
                linear_tolerance_ns: 3e-10,
                angular_tolerance_nms: 3e-10,
                work_tolerance_j: 3e-10,
            },
            linear_speed_tolerance_m_s: 3e-10,
            angular_speed_tolerance_rad_s: 3e-10,
        },
        source_substeps: 1,
        limits: ViscousNetworkLimits {
            maximum_bodies: 32,
            maximum_connectors: 64,
            maximum_pair_relaxations: 100_000,
            maximum_live_bytes: 1_000_000,
        },
    }
}
fn matrix(c: [[f64; 6]; 6]) -> PassiveBushingMatrix {
    PassiveBushingMatrix::try_new(
        c,
        BushingCoordinateScales {
            translation_m: 1.,
            rotation_rad: 1.,
        },
    )
    .unwrap()
}
fn axial(bodies: [usize; 2], coefficient: f64) -> ViscousNetworkLink {
    let mut c = [[0.; 6]; 6];
    c[0][0] = coefficient;
    ViscousNetworkLink {
        bodies,
        arms_body_m: [Vec3::ZERO; 2],
        damping: matrix(c),
    }
}
fn state(p: MassProperties, x: Vec3, q: UnitQuaternion, v: Vec3, w: Vec3) -> RigidBodyState {
    RigidBodyState::new(
        Pose::new(x, q).unwrap(),
        v.scale(p.mass()),
        p.angular_momentum_body_checked(w).unwrap(),
    )
    .unwrap()
}
fn line() -> (
    [RigidBodyState; 3],
    [MassProperties; 3],
    [ViscousNetworkLink; 2],
) {
    let p = [MassProperties::new(2., Vec3::ZERO, Vec3::new(2., 3., 4.)).unwrap(); 3];
    let s = core::array::from_fn(|i| {
        state(
            p[i],
            X.scale(i as f64 - 1.),
            UnitQuaternion::IDENTITY,
            X.scale([1.2, -0.8, 0.2][i]),
            Vec3::ZERO,
        )
    });
    (s, p, [axial([0, 1], 0.6), axial([1, 2], 0.6)])
}

#[test]
fn single_edge_motion_is_the_existing_native_pair_without_new_physics() {
    let (s, p, links) = line();
    let pair = ViscousMotionLaw {
        arms_body_m: links[0].arms_body_m,
        damping: links[0].damping.clone(),
    };
    let expected = step_viscous_connector(
        [s[0], s[1]],
        [p[0], p[1]],
        &pair,
        0.3,
        controls().source,
        || false,
    )
    .unwrap();
    let result =
        step_viscous_network(&s[..2], &p[..2], &links[..1], 0.3, controls(), || false).unwrap();
    assert_eq!(result.states_after(), expected.states_after);
    assert_eq!(result.dissipated_energy_j, expected.dissipated_energy_j);
    assert_eq!(
        result.connector_dissipation_j,
        [expected.dissipated_energy_j]
    );
    assert_eq!(
        result.kinetic_energy_change_j,
        expected.kinetic_energy_change_j
    );
    assert_eq!(result.pair_relaxations, 2);
}

#[test]
fn shared_chain_positions_velocities_and_heat_refine_to_both_continuous_graph_modes() {
    let (s, p, links) = line();
    let t: f64 = 0.8;
    // Path graph modes [1,0,-1] and [1,-2,1], both amplitude 1/2.
    // Translational rigid mode is 0.2 m/s; decay rates are c/m and 3c/m.
    let e1 = (-0.3 * t).exp();
    let e3 = (-0.9 * t).exp();
    let v = [
        0.2 + 0.5 * e1 + 0.5 * e3,
        0.2 - e3,
        0.2 - 0.5 * e1 + 0.5 * e3,
    ];
    let d1 = 0.5 * -(-0.3 * t).exp_m1() / 0.3;
    let d3 = 0.5 * -(-0.9 * t).exp_m1() / 0.9;
    let x = [
        -1. + 0.2 * t + d1 + d3,
        0.2 * t - 2. * d3,
        1. + 0.2 * t - d1 + d3,
    ];
    let loss = 0.5 * -(-0.6 * t).exp_m1() + 1.5 * -(-1.8 * t).exp_m1();
    let mut errors = [0.; 3];
    for (rung, steps) in [4, 8, 16].into_iter().enumerate() {
        let mut states = s.to_vec();
        let mut heat = 0.;
        for _ in 0..steps {
            let out =
                step_viscous_network(&states, &p, &links, t / steps as f64, controls(), || false)
                    .unwrap();
            heat += out.dissipated_energy_j;
            near(
                out.connector_dissipation_j.iter().sum(),
                out.dissipated_energy_j,
                2e-14,
            );
            near(out.work_heat_residual_j, 0., 2e-13);
            vector(out.linear_momentum_residual_ns, Vec3::ZERO, 3e-14);
            vector(out.angular_momentum_residual_nms, Vec3::ZERO, 3e-14);
            assert_eq!(out.pair_relaxations, 6);
            states = out.into_states_after();
        }
        let mut error = (heat - loss).abs();
        for i in 0..3 {
            error += (states[i].pose().position_world().x - x[i]).abs();
            error += (states[i].center_of_mass_velocity_world(p[i]).unwrap().x - v[i]).abs();
        }
        near(
            states
                .iter()
                .map(|s| s.pose().position_world().x)
                .sum::<f64>()
                / 3.,
            0.2 * t,
            3e-14,
        );
        errors[rung] = error;
    }
    for pair in errors.windows(2) {
        assert!((3.8..4.2).contains(&(pair[0] / pair[1])), "{errors:?}");
    }
}

fn full() -> (
    [RigidBodyState; 3],
    [MassProperties; 3],
    [ViscousNetworkLink; 3],
) {
    let properties = core::array::from_fn(|i| {
        MassProperties::from_principal_frame(
            2. + i as f64,
            Vec3::ZERO,
            Vec3::new(1. + i as f64, 2. + i as f64, 3. + i as f64),
            UnitQuaternion::from_axis_angle(Vec3::new(1., 2., -0.5), 0.4 + 0.1 * i as f64).unwrap(),
        )
        .unwrap()
    });
    let states = core::array::from_fn(|i| {
        let f = i as f64;
        state(
            properties[i],
            Vec3::new(f - 0.7, 0.3 - 0.2 * f, -0.1 + 0.4 * f),
            UnitQuaternion::from_axis_angle(Vec3::new(1., 2. - f, -0.5), 0.2 + 0.17 * f).unwrap(),
            Vec3::new(0.5 - 0.3 * f, -0.2 + 0.4 * f, 0.1 - 0.2 * f),
            Vec3::new(0.4 - 0.2 * f, 0.3 - 0.1 * f, -0.2 + 0.3 * f),
        )
    });
    let links = core::array::from_fn(|edge| {
        let d = [0.6, 0.8, 0.5, 0.2, 0.3, 0.4];
        let b = [0.1, -0.2, 0.05, 0.12, 0.03, -0.14];
        let c = core::array::from_fn(|i| {
            core::array::from_fn(|j| b[i] * b[j] + if i == j { d[i] } else { 0. })
        });
        ViscousNetworkLink {
            bodies: [edge, (edge + 1) % 3],
            arms_body_m: [Vec3::new(0.2, -0.1, 0.15), Vec3::new(-0.15, 0.3, 0.1)],
            damping: matrix(c),
        }
    });
    (states, properties, links)
}

// Independent continuous Newton-Euler equations, not native pair/source maps.
// Full nonprincipal inertia comes from the existing mass-property owner only.
#[derive(Clone, Copy)]
struct Reference {
    x: [Vec3; 3],
    p: [Vec3; 3],
    h: [Vec3; 3],
    q: [[f64; 4]; 3],
    heat: [f64; 3],
}
impl Reference {
    fn from_states(s: [RigidBodyState; 3]) -> Self {
        Self {
            x: s.map(|s| s.pose().position_world()),
            p: s.map(|s| s.linear_momentum_world()),
            h: s.map(|s| s.angular_momentum_body()),
            q: s.map(|s| s.pose().orientation().components()),
            heat: [0.; 3],
        }
    }
    fn shifted(self, rhs: Self, h: f64) -> Self {
        Self {
            x: core::array::from_fn(|i| self.x[i].add(rhs.x[i].scale(h))),
            p: core::array::from_fn(|i| self.p[i].add(rhs.p[i].scale(h))),
            h: core::array::from_fn(|i| self.h[i].add(rhs.h[i].scale(h))),
            q: core::array::from_fn(|i| core::array::from_fn(|j| self.q[i][j] + h * rhs.q[i][j])),
            heat: core::array::from_fn(|i| self.heat[i] + h * rhs.heat[i]),
        }
    }
}
fn derivative(
    s: Reference,
    properties: [MassProperties; 3],
    links: &[ViscousNetworkLink; 3],
) -> Reference {
    let q =
        s.q.map(|q| UnitQuaternion::new(q[0], q[1], q[2], q[3]).unwrap());
    let w = core::array::from_fn::<_, 3, _>(|i| {
        properties[i].angular_velocity_body_checked(s.h[i]).unwrap()
    });
    let world_w = core::array::from_fn::<_, 3, _>(|i| q[i].rotate_body_to_world(w[i]));
    let velocity = core::array::from_fn::<_, 3, _>(|i| s.p[i].scale(1. / properties[i].mass()));
    let mut out = Reference {
        x: velocity,
        p: [Vec3::ZERO; 3],
        h: core::array::from_fn(|i| w[i].cross(s.h[i]).scale(-1.)),
        q: core::array::from_fn(|i| {
            let [a, b, c, d] = s.q[i];
            let o = w[i];
            [
                -0.5 * (b * o.x + c * o.y + d * o.z),
                0.5 * (a * o.x + c * o.z - d * o.y),
                0.5 * (a * o.y + d * o.x - b * o.z),
                0.5 * (a * o.z + b * o.y - c * o.x),
            ]
        }),
        heat: [0.; 3],
    };
    for (edge, link) in links.iter().enumerate() {
        let [a, b] = link.bodies;
        let arms = [
            q[a].rotate_body_to_world(link.arms_body_m[0]),
            q[b].rotate_body_to_world(link.arms_body_m[1]),
        ];
        let separation = s.x[b].add(arms[1]).sub(s.x[a].add(arms[0]));
        let slip = q[a].rotate_world_to_body(
            velocity[b].add(world_w[b].cross(arms[1])).sub(
                velocity[a]
                    .add(world_w[a].cross(arms[0]))
                    .add(world_w[a].cross(separation)),
            ),
        );
        let spin = q[a].rotate_world_to_body(world_w[b].sub(world_w[a]));
        let twist = [slip.x, slip.y, slip.z, spin.x, spin.y, spin.z];
        let c = link.damping.native_coefficients();
        let load: [f64; 6] =
            core::array::from_fn(|i| -(0..6).map(|j| c[i][j] * twist[j]).sum::<f64>());
        let force = q[a].rotate_body_to_world(Vec3::new(load[0], load[1], load[2]));
        let couple = q[a].rotate_body_to_world(Vec3::new(load[3], load[4], load[5]));
        out.p[a] = out.p[a].sub(force);
        out.p[b] = out.p[b].add(force);
        out.h[a] = out.h[a].add(
            q[a].rotate_world_to_body(
                arms[0]
                    .cross(force)
                    .scale(-1.)
                    .sub(couple)
                    .sub(separation.cross(force)),
            ),
        );
        out.h[b] = out.h[b].add(q[b].rotate_world_to_body(arms[1].cross(force).add(couple)));
        out.heat[edge] = -(0..6).map(|i| twist[i] * load[i]).sum::<f64>();
    }
    out
}
fn reference(
    initial: [RigidBodyState; 3],
    properties: [MassProperties; 3],
    links: &[ViscousNetworkLink; 3],
    duration: f64,
) -> Reference {
    let mut s = Reference::from_states(initial);
    let h = duration / 2048.;
    for _ in 0..2048 {
        let k1 = derivative(s, properties, links);
        let k2 = derivative(s.shifted(k1, 0.5 * h), properties, links);
        let k3 = derivative(s.shifted(k2, 0.5 * h), properties, links);
        let k4 = derivative(s.shifted(k3, h), properties, links);
        s = s
            .shifted(k1, h / 6.)
            .shifted(k2, h / 3.)
            .shifted(k3, h / 3.)
            .shifted(k4, h / 6.);
    }
    s
}

#[test]
fn moving_cycle_with_full_coupled_spatial_losses_matches_continuous_newton_euler_flow() {
    let (initial, properties, links) = full();
    let duration = 0.3;
    let expected = reference(initial, properties, &links, duration);
    let mut errors = [0.; 3];
    for (rung, steps) in [2, 4, 8].into_iter().enumerate() {
        let mut states = initial.to_vec();
        let mut heat = [0.; 3];
        for _ in 0..steps {
            let result = step_viscous_network(
                &states,
                &properties,
                &links,
                duration / steps as f64,
                controls(),
                || false,
            )
            .unwrap();
            for i in 0..3 {
                heat[i] += result.connector_dissipation_j[i];
            }
            near(result.work_heat_residual_j, 0., 4e-11);
            vector(result.linear_momentum_residual_ns, Vec3::ZERO, 3e-12);
            vector(result.angular_momentum_residual_nms, Vec3::ZERO, 3e-11);
            assert!(result.connector_dissipation_j.iter().all(|h| *h >= 0.));
            states = result.into_states_after();
        }
        let mut error = 0.;
        for i in 0..3 {
            error += norm(states[i].pose().position_world().sub(expected.x[i]));
            error += norm(states[i].linear_momentum_world().sub(expected.p[i]));
            error += norm(states[i].angular_momentum_body().sub(expected.h[i]));
            error += (heat[i] - expected.heat[i]).abs();
            for (a, b) in states[i]
                .pose()
                .orientation()
                .components()
                .into_iter()
                .zip(expected.q[i])
            {
                error += (a - b).abs();
            }
        }
        errors[rung] = error;
    }
    for pair in errors.windows(2) {
        assert!((3.4..4.5).contains(&(pair[0] / pair[1])), "{errors:?}");
    }
}

#[test]
fn isolated_bodies_and_zero_laws_keep_the_native_loaded_motion() {
    let (states, properties, _) = full();
    let links = [axial([0, 1], 0.)];
    let wrenches = [Wrench {
        force_world: Vec3::new(0.2, -0.1, 0.3),
        torque_body: Vec3::new(0.03, -0.02, 0.01),
    }; 3];
    let gravity = Gravity::new(Vec3::new(0., -9.81, 0.)).unwrap();
    let dt = 0.04;
    let result = step_viscous_network_with_loads(
        &states,
        &properties,
        &links,
        HeldExternalLoads {
            wrenches: &wrenches,
            gravity,
        },
        dt,
        controls(),
        || false,
    )
    .unwrap();
    for i in 0..3 {
        let expected = RigidBodyIntegrator::new(gravity)
            .step_cayley(states[i], properties[i], wrenches[i], dt)
            .unwrap();
        assert_eq!(result.states_after()[i], expected.state_after);
    }
    assert_eq!(result.pair_relaxations, 0);
    assert_eq!(result.connector_dissipation_j, [0.]);
    assert_eq!(result.dissipated_energy_j, 0.);
    near(
        result.mechanical_energy_change_j,
        result.external_loads.applied_work_j,
        1e-12,
    );
    near(
        result.gravitational_potential_energy_change_j,
        -result.external_loads.gravity_work_j,
        0.,
    );
}

#[test]
fn gravity_translates_the_entire_shared_graph_without_becoming_viscous_heat() {
    let (s, p, links) = line();
    let dt = 0.1;
    let gravity = Gravity::new(Vec3::new(0., -9.81, 0.)).unwrap();
    let free = step_viscous_network(&s, &p, &links, dt, controls(), || false).unwrap();
    let falling = step_viscous_network_with_loads(
        &s,
        &p,
        &links,
        HeldExternalLoads {
            wrenches: &[Wrench::ZERO; 3],
            gravity,
        },
        dt,
        controls(),
        || false,
    )
    .unwrap();
    near(falling.dissipated_energy_j, free.dissipated_energy_j, 3e-14);
    for i in 0..3 {
        vector(
            falling.states_after()[i].linear_momentum_world(),
            free.states_after()[i]
                .linear_momentum_world()
                .add(gravity.acceleration_world().scale(dt * p[i].mass())),
            3e-14,
        );
        vector(
            falling.states_after()[i].pose().position_world(),
            free.states_after()[i]
                .pose()
                .position_world()
                .add(gravity.acceleration_world().scale(0.5 * dt * dt)),
            3e-14,
        );
    }
    near(falling.work_heat_residual_j, 0., 2e-13);
}

#[test]
fn rotating_material_point_loads_reuse_pair_work_and_advance_an_isolated_body() {
    let (s, p, links) = full();
    let prepared: [_; 3] = core::array::from_fn(|i| {
        PreparedSpatialLoads::new(
            Wrench::ZERO,
            Vec3::new(0.03, -0.02, 0.01),
            &[WorldPointForce {
                arm_body: Vec3::new(0.1, -0.2, 0.3),
                force_world: Vec3::new(0.2, 0.1 * i as f64, -0.1),
            }],
            || false,
        )
        .unwrap()
    });
    let midpoint = SpatialMidpointConfig {
        maximum_iterations: 24,
        maximum_backtracks: 12,
        angular_momentum_tolerance_nms: 1e-13,
        relative_tolerance: 1e-13,
        maximum_cayley_increment: 0.5,
    };
    let dt = 0.04;
    let out = step_viscous_network_with_spatial_loads(
        &s,
        &p,
        &links[..1],
        HeldSpatialLoads {
            bodies: &prepared,
            gravity: Gravity::ZERO,
            midpoint,
        },
        dt,
        controls(),
        || false,
    )
    .unwrap();
    let pair_law = ViscousMotionLaw {
        arms_body_m: links[0].arms_body_m,
        damping: links[0].damping.clone(),
    };
    let pair = step_viscous_connector_with_spatial_loads(
        [s[0], s[1]],
        [p[0], p[1]],
        &pair_law,
        HeldSpatialLoads {
            bodies: &prepared[..2],
            gravity: Gravity::ZERO,
            midpoint,
        },
        dt,
        controls().source,
        || false,
    )
    .unwrap();
    assert_eq!(&out.states_after()[..2], &pair.states_after);
    assert_eq!(out.dissipated_energy_j, pair.dissipated_energy_j);
    let isolated = RigidBodyIntegrator::new(Gravity::ZERO)
        .step_spatial_midpoint(s[2], p[2], &prepared[2], dt, midpoint, || false)
        .unwrap();
    assert_eq!(out.states_after()[2], isolated.motion.step.state_after);
    near(
        out.external_loads.applied_force_work_j,
        pair.external_loads.applied_force_work_j + isolated.motion.force_work_j,
        2e-14,
    );
    near(
        out.external_loads.applied_torque_work_j,
        pair.external_loads.applied_torque_work_j + isolated.motion.couple_work_j,
        2e-14,
    );
    near(out.work_heat_residual_j, 0., 3e-12);
}

#[test]
fn complete_work_and_memory_caps_and_late_cancellation_are_atomic() {
    let (s, p, links) = full();
    let saved_links = links.clone();
    let mut polls = 0;
    let accepted = step_viscous_network(&s, &p, &links, 0.07, controls(), || {
        polls += 1;
        false
    })
    .unwrap();
    let mut exact = controls();
    exact.limits.maximum_live_bytes = accepted.planned_live_bytes;
    exact.limits.maximum_pair_relaxations = accepted.pair_relaxations;
    assert_eq!(
        step_viscous_network(&s, &p, &links, 0.07, exact, || false).unwrap(),
        accepted
    );
    let mut insufficient = exact;
    insufficient.limits.maximum_live_bytes -= 1;
    assert!(step_viscous_network(&s, &p, &links, 0.07, insufficient, || false).is_err());
    insufficient = exact;
    insufficient.limits.maximum_pair_relaxations -= 1;
    assert!(step_viscous_network(&s, &p, &links, 0.07, insufficient, || false).is_err());
    let mut later = 0;
    assert_eq!(
        step_viscous_network(&s, &p, &links, 0.07, controls(), || {
            later += 1;
            later == polls
        }),
        Err(ViscousNetworkMotionError::Cancelled)
    );
    for duration in [0., -0.1, f64::NAN, f64::from_bits(1)] {
        assert!(step_viscous_network(&s, &p, &links, duration, controls(), || false).is_err());
    }
    assert!(step_viscous_network(&s, &p[..2], &links, 0.07, controls(), || false).is_err());
    insufficient = controls();
    insufficient.source_substeps = 0;
    assert!(step_viscous_network(&s, &p, &links, 0.07, insufficient, || false).is_err());
    assert_eq!(links, saved_links);
    assert_eq!(
        step_viscous_network(&s, &p, &links, 0.07, controls(), || false).unwrap(),
        accepted
    );
}
