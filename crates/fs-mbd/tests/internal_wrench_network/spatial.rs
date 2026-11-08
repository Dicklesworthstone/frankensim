//! G1/G3/G4: physical material-point loads on complete finite actuator networks.

use super::*;
use fs_mbd::spatial_loads::{
    PreparedSpatialLoads, WorldPointForce,
    internal_motion::{
        HeldSpatialLoads, InternalWrenchStepError,
        network::step_internal_wrench_network_with_spatial_loads, step_internal_wrench_with_loads,
        step_internal_wrench_with_spatial_loads,
    },
    midpoint::{SpatialMidpointConfig, SpatialMidpointError},
};

fn midpoint_config() -> SpatialMidpointConfig {
    SpatialMidpointConfig {
        maximum_iterations: 12,
        maximum_backtracks: 8,
        angular_momentum_tolerance_nms: 1e-14,
        relative_tolerance: 1e-13,
        maximum_cayley_increment: 1.0,
    }
}

struct Sources {
    held: [Wrench; 3],
    world_couples: [Vec3; 3],
    points: [[WorldPointForce; 2]; 3],
    gravity: Gravity,
}
impl Sources {
    fn new() -> Self {
        Self {
            held: core::array::from_fn(|i| Wrench {
                force_world: [X.scale(0.3), Y.scale(-0.2), Z.scale(0.4)][i],
                torque_body: [
                    Vec3::new(0.1, -0.2, 0.3),
                    Vec3::new(-0.2, 0.4, 0.1),
                    Vec3::new(0.3, 0.1, -0.2),
                ][i],
            }),
            world_couples: [
                Vec3::new(0.2, -0.1, 0.4),
                Vec3::new(-0.1, 0.2, -0.3),
                Vec3::new(0.3, 0.1, 0.2),
            ],
            points: core::array::from_fn(|i| {
                [
                    WorldPointForce {
                        arm_body: Vec3::new(0.3 + i as f64 * 0.1, -0.2, 0.1),
                        force_world: Vec3::new(-0.4, 0.6, 0.2),
                    },
                    WorldPointForce {
                        arm_body: Vec3::new(-0.1, 0.4, -0.2 - i as f64 * 0.1),
                        force_world: Vec3::new(0.2, -0.1, -0.3),
                    },
                ]
            }),
            gravity: Gravity::new(Vec3::new(0.2, -0.4, 0.1)).unwrap(),
        }
    }

    fn prepared(&self) -> [PreparedSpatialLoads; 3] {
        core::array::from_fn(|i| {
            PreparedSpatialLoads::new(self.held[i], self.world_couples[i], &self.points[i], || {
                false
            })
            .unwrap()
        })
    }

    // Direct physical point geometry, independent of the compiled moment bank
    // and native effective body torque. Returns force/couple/gravity work and impulses.
    fn endpoint_account(
        &self,
        before: &[RigidBodyState],
        after: &[RigidBodyState],
        time: f64,
        origin: Vec3,
    ) -> (f64, f64, f64, Vec3, Vec3) {
        let p = properties();
        let (mut work, mut couple_work, mut gravity_work, mut impulse, mut angular) =
            (0., 0., 0., Vec3::ZERO, Vec3::ZERO);
        for i in 0..before.len() {
            let (a, b) = (before[i].pose(), after[i].pose());
            let dx = b.position_world().sub(a.position_world());
            let mut force = self.held[i]
                .force_world
                .add(self.gravity.acceleration_world().scale(p[i].mass()));
            let mut moment = self.world_couples[i].add(
                a.orientation()
                    .rotate_body_to_world(self.held[i].torque_body)
                    .add(
                        b.orientation()
                            .rotate_body_to_world(self.held[i].torque_body),
                    )
                    .scale(0.5),
            );
            work += self.held[i].force_world.dot(dx);
            // Independent finite Cayley increment from the represented relative
            // quaternion, not the receipt's angular-velocity or torque fields.
            let [w0, x0, y0, z0] = a.orientation().components();
            let [w1, x1, y1, z1] = b.orientation().components();
            let (v0, v1) = (Vec3::new(x0, y0, z0), Vec3::new(x1, y1, z1));
            let c = v1
                .scale(w0)
                .sub(v0.scale(w1))
                .sub(v0.cross(v1))
                .scale(2. / (w0 * w1 + v0.dot(v1)));
            couple_work += self.held[i].torque_body.dot(c)
                + self.world_couples[i].dot(rotate([w0, x0, y0, z0], c));
            gravity_work += p[i].mass() * self.gravity.acceleration_world().dot(dx);
            for point in self.points[i] {
                let r0 = a.orientation().rotate_body_to_world(point.arm_body);
                let r1 = b.orientation().rotate_body_to_world(point.arm_body);
                force = force.add(point.force_world);
                moment = moment.add(r0.add(r1).scale(0.5).cross(point.force_world));
                work += point.force_world.dot(dx.add(r1.sub(r0)));
            }
            let radius = a
                .position_world()
                .add(b.position_world())
                .scale(0.5)
                .sub(origin);
            impulse = impulse.add(force.scale(time));
            angular = angular.add(moment.add(radius.cross(force)).scale(time));
        }
        (work, couple_work, gravity_work, impulse, angular)
    }
}

#[test]
fn g1_rotating_point_loads_supply_actual_work_and_impulse_to_pair_and_network() {
    let sources = Sources::new();
    let loads = sources.prepared();
    let (s, p, links) = (initial(), properties(), drives());
    let origin = controls().balance.angular_origin_world_m;
    let environment = HeldSpatialLoads {
        bodies: &loads[..2],
        gravity: sources.gravity,
        midpoint: midpoint_config(),
    };
    let duration = 0.08;
    let pair = step_internal_wrench_with_spatial_loads(
        [s[0], s[1]],
        [p[0], p[1]],
        links[0].drive,
        environment,
        duration,
        controls().balance,
        || false,
    )
    .unwrap();
    let network = step_internal_wrench_network_with_spatial_loads(
        &s[..2],
        &p[..2],
        &links[..1],
        environment,
        duration,
        controls(),
        || false,
    )
    .unwrap();
    for i in 0..2 {
        assert!(error(pair.states_after[i], network.states_after[i]) < 2e-12);
    }
    near(pair.actuator_work_j, network.actuator_work_j, 2e-12);
    let (work, couple_work, gravity_work, impulse, angular) =
        sources.endpoint_account(&s[..2], &pair.states_after, duration, origin);
    near(pair.external_loads.applied_force_work_j, work, 2e-13);
    near(
        pair.external_loads.applied_torque_work_j,
        couple_work,
        2e-13,
    );
    near(pair.external_loads.gravity_work_j, gravity_work, 2e-13);
    vector(pair.external_loads.linear_impulse_world_ns, impulse, 2e-13);
    vector(
        pair.external_loads.angular_impulse_world_nms,
        angular,
        2e-13,
    );
    let old = totals(&s[..2], origin);
    let new = totals(&pair.states_after, origin);
    vector(new.0.sub(old.0), impulse, 2e-12);
    vector(new.1.sub(old.1), angular, 2e-12);
    let expected_kinetic = kinetic(&pair.states_after) - kinetic(&s[..2]);
    near(
        expected_kinetic,
        work + gravity_work + couple_work + pair.actuator_work_j,
        2e-12,
    );
    // The load really changes the rotational trajectory; a COM-only surrogate
    // omits both rotating point moments and the world-held free couples.
    let com = step_internal_wrench_with_loads(
        [s[0], s[1]],
        [p[0], p[1]],
        links[0].drive,
        HeldExternalLoads {
            wrenches: &sources.held[..2],
            gravity: sources.gravity,
        },
        duration,
        controls().balance,
        || false,
    )
    .unwrap();
    assert!(
        pair.states_after[0]
            .angular_momentum_body()
            .sub(com.states_after[0].angular_momentum_body())
            .norm_squared()
            > 1e-5
    );
    let shift = Vec3::new(5., -3., 2.);
    let mut shifted_controls = controls().balance;
    shifted_controls.angular_origin_world_m = origin.add(shift);
    let shifted = step_internal_wrench_with_spatial_loads(
        [s[0], s[1]],
        [p[0], p[1]],
        links[0].drive,
        environment,
        duration,
        shifted_controls,
        || false,
    )
    .unwrap();
    assert_eq!(shifted.states_after, pair.states_after);
    vector(
        shifted.external_loads.angular_impulse_world_nms,
        angular.sub(shift.cross(impulse)),
        3e-13,
    );
}

#[test]
fn g1_g3_spatial_actuator_cycle_refines_to_independent_newton_euler_and_power_integrals() {
    let sources = Sources::new();
    let loads = sources.prepared();
    let (initial, p, links) = (initial(), properties(), drives());
    // Original Newton/Euler equations plus each actuator's power and distinct
    // material-force, free-couple and gravity integrals. No native inertia
    // inversion, compiled moments, integrator or work receipt is used here.
    let derivative = |y: [f64; 45]| {
        let body: Ode = core::array::from_fn(|i| y[i]);
        let base = rhs(body, &links);
        let mut out: [f64; 45] = core::array::from_fn(|i| if i < 42 { base[i] } else { 0. });
        for i in 0..3 {
            let n = 13 * i;
            let q = [y[n + 6], y[n + 7], y[n + 8], y[n + 9]];
            let omega_body = omega(triple(&body, n + 10), i);
            let omega_world = rotate(q, omega_body);
            let velocity = triple(&body, n + 3).scale(1. / p[i].mass());
            let gravity = sources.gravity.acceleration_world().scale(p[i].mass());
            let mut force = sources.held[i].force_world.add(gravity);
            let mut moment = sources.world_couples[i];
            out[42] += sources.held[i].force_world.dot(velocity);
            out[43] += sources.world_couples[i].dot(omega_world)
                + sources.held[i].torque_body.dot(omega_body);
            out[44] += gravity.dot(velocity);
            for point in sources.points[i] {
                let arm = rotate(q, point.arm_body);
                force = force.add(point.force_world);
                moment = moment.add(arm.cross(point.force_world));
                out[42] += point.force_world.dot(velocity.add(omega_world.cross(arm)));
            }
            let torque =
                rotate([q[0], -q[1], -q[2], -q[3]], moment).add(sources.held[i].torque_body);
            for (j, f) in [force.x, force.y, force.z].into_iter().enumerate() {
                out[n + 3 + j] += f;
            }
            for (j, t) in [torque.x, torque.y, torque.z].into_iter().enumerate() {
                out[n + 10 + j] += t;
            }
        }
        out
    };
    let time = 0.4;
    let integrate = |count| {
        let start = ode_initial(initial);
        let mut y: [f64; 45] = core::array::from_fn(|i| if i < 42 { start[i] } else { 0. });
        let dt = time / f64::from(count);
        for _ in 0..count {
            let a = derivative(y);
            let b = derivative(core::array::from_fn(|i| y[i] + 0.5 * dt * a[i]));
            let c = derivative(core::array::from_fn(|i| y[i] + 0.5 * dt * b[i]));
            let d = derivative(core::array::from_fn(|i| y[i] + dt * c[i]));
            y = core::array::from_fn(|i| y[i] + dt * (a[i] + 2. * b[i] + 2. * c[i] + d[i]) / 6.);
        }
        y
    };
    let coarse = integrate(1024);
    let reference = integrate(2048);
    for i in 0..45 {
        near(coarse[i], reference[i], 3e-12);
    }
    let target = core::array::from_fn::<_, 3, _>(|i| {
        state(
            Vec3::new(
                reference[13 * i],
                reference[13 * i + 1],
                reference[13 * i + 2],
            ),
            UnitQuaternion::new(
                reference[13 * i + 6],
                reference[13 * i + 7],
                reference[13 * i + 8],
                reference[13 * i + 9],
            )
            .unwrap(),
            Vec3::new(
                reference[13 * i + 3],
                reference[13 * i + 4],
                reference[13 * i + 5],
            ),
            Vec3::new(
                reference[13 * i + 10],
                reference[13 * i + 11],
                reference[13 * i + 12],
            ),
        )
    });
    let (mut motion_errors, mut work_errors) = (Vec::new(), Vec::new());
    for count in [20, 40, 80] {
        let mut s = initial.to_vec();
        let mut works = [0.; 6];
        let dt = time / f64::from(count);
        for _ in 0..count {
            let out = step_internal_wrench_network_with_spatial_loads(
                &s,
                &p,
                &links,
                HeldSpatialLoads {
                    bodies: &loads,
                    gravity: sources.gravity,
                    midpoint: midpoint_config(),
                },
                dt,
                controls(),
                || false,
            )
            .unwrap();
            for stage in &out.half_impulses {
                for (i, connector) in stage.connectors.iter().enumerate() {
                    works[i] += connector.connector_work_j;
                }
            }
            works[3] += out.external_loads.applied_force_work_j;
            works[4] += out.external_loads.applied_torque_work_j;
            works[5] += out.external_loads.gravity_work_j;
            s = out.states_after;
        }
        motion_errors.push((0..3).map(|i| error(s[i], target[i])).fold(0., f64::max));
        work_errors.push(
            (0..6)
                .map(|i| (works[i] - reference[39 + i]).abs())
                .fold(0., f64::max),
        );
    }
    for values in [&motion_errors, &work_errors] {
        assert!(values[2] < 1e-4, "{values:?}");
        assert!(
            values[1] < 0.28 * values[0] && values[2] < 0.28 * values[1],
            "{values:?}"
        );
    }
}

#[test]
fn g3_g4_spatial_zero_loads_preserve_motion_and_native_refusals_publish_no_partial_network() {
    let (s, p, links) = (initial(), properties(), drives());
    let zero = [PreparedSpatialLoads::ZERO; 3];
    let environment = HeldSpatialLoads {
        bodies: &zero,
        gravity: Gravity::ZERO,
        midpoint: midpoint_config(),
    };
    let spatial = step_internal_wrench_network_with_spatial_loads(
        &s,
        &p,
        &links,
        environment,
        0.05,
        controls(),
        || false,
    )
    .unwrap();
    let native = step_internal_wrench_network(&s, &p, &links, 0.05, controls(), || false).unwrap();
    assert_eq!(spatial, native);
    let pair = step_internal_wrench_with_spatial_loads(
        [s[0], s[1]],
        [p[0], p[1]],
        links[0].drive,
        HeldSpatialLoads {
            bodies: &zero[..2],
            ..environment
        },
        0.05,
        controls().balance,
        || false,
    )
    .unwrap();
    assert_eq!(
        pair,
        step_internal_wrench(
            [s[0], s[1]],
            [p[0], p[1]],
            links[0].drive,
            0.05,
            controls().balance,
            || false
        )
        .unwrap()
    );
    let sources = Sources::new();
    let loads = sources.prepared();
    let environment = HeldSpatialLoads {
        bodies: &loads,
        gravity: sources.gravity,
        midpoint: midpoint_config(),
    };
    let mut count = 0;
    let accepted = step_internal_wrench_network_with_spatial_loads(
        &s,
        &p,
        &links,
        environment,
        0.05,
        controls(),
        || {
            count += 1;
            false
        },
    )
    .unwrap();
    for stop in [count / 2, count] {
        let mut calls = 0;
        assert_eq!(
            step_internal_wrench_network_with_spatial_loads(
                &s,
                &p,
                &links,
                environment,
                0.05,
                controls(),
                || {
                    calls += 1;
                    calls == stop
                }
            ),
            Err(InternalWrenchNetworkStepError::Cancelled)
        );
        assert_eq!(s, initial());
    }
    let mut limited = environment;
    limited.midpoint.maximum_cayley_increment = 1e-14;
    assert!(matches!(
        step_internal_wrench_network_with_spatial_loads(
            &s,
            &p,
            &links,
            limited,
            0.05,
            controls(),
            || false
        ),
        Err(InternalWrenchNetworkStepError::Motion(
            InternalWrenchStepError::SpatialMidpoint(SpatialMidpointError::IncrementLimit)
        ))
    ));
    assert_eq!(
        accepted,
        step_internal_wrench_network_with_spatial_loads(
            &s,
            &p,
            &links,
            environment,
            0.05,
            controls(),
            || false
        )
        .unwrap()
    );
}
