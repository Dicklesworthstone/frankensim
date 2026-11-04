//! Simultaneous internal actuation: recoil, source work and finite trajectories.

use fs_mbd::spatial_loads::internal_impulse::{
    InternalWrenchImpulse, InternalWrenchImpulseControls, InternalWrenchImpulseError,
    apply_internal_wrench_impulse,
    network::{
        InternalWrenchLink, InternalWrenchNetworkControls, InternalWrenchNetworkError,
        apply_internal_wrench_network,
    },
};
use fs_mbd::spatial_loads::internal_motion::{
    PrescribedInternalWrench,
    network::{
        InternalWrenchNetworkStepError, PrescribedInternalLink, step_internal_wrench_network,
    },
    step_internal_wrench,
};
use fs_mbd::{
    Gravity, MassProperties, Pose, RigidBodyIntegrator, RigidBodyState, UnitQuaternion, Vec3,
    Wrench,
};

const X: Vec3 = Vec3::new(1., 0., 0.);
const Y: Vec3 = Vec3::new(0., 1., 0.);
const Z: Vec3 = Vec3::new(0., 0., 1.);
fn near(a: f64, b: f64, t: f64) {
    assert!((a - b).abs() <= t, "{a:e} != {b:e}, tolerance {t:e}");
}
fn norm(v: Vec3) -> f64 {
    v.x.hypot(v.y).hypot(v.z)
}
fn vector(a: Vec3, b: Vec3, t: f64) {
    assert!(norm(a.sub(b)) <= t, "{a:?} != {b:?}, tolerance {t:e}");
}
fn state(x: Vec3, q: UnitQuaternion, p: Vec3, h: Vec3) -> RigidBodyState {
    RigidBodyState::new(Pose::new(x, q).unwrap(), p, h).unwrap()
}
fn controls() -> InternalWrenchNetworkControls {
    InternalWrenchNetworkControls {
        balance: InternalWrenchImpulseControls {
            angular_origin_world_m: Vec3::new(0.3, -0.2, 0.1),
            linear_tolerance_ns: 3e-10,
            angular_tolerance_nms: 3e-10,
            work_tolerance_j: 3e-10,
        },
        maximum_bodies: 32,
        maximum_connectors: 64,
        maximum_live_bytes: 1_000_000,
    }
}
fn properties() -> [MassProperties; 3] {
    core::array::from_fn(|i| {
        MassProperties::from_principal_frame(
            [2., 3., 5.][i],
            Vec3::ZERO,
            Vec3::new(2. + i as f64, 3. + i as f64, 4. + i as f64),
            UnitQuaternion::from_axis_angle(Z, core::f64::consts::FRAC_PI_4).unwrap(),
        )
        .unwrap()
    })
}
// Explicit inverse of I=[[a,-.5,0],[-.5,a,0],[0,0,4+i]], independent of the
// production full-inertia and native rotational stepping implementations.
fn omega(h: Vec3, i: usize) -> Vec3 {
    let a = 2.5 + i as f64;
    let det = (2. + i as f64) * (3. + i as f64);
    Vec3::new(
        (a * h.x + 0.5 * h.y) / det,
        (0.5 * h.x + a * h.y) / det,
        h.z / (4. + i as f64),
    )
}
fn kinetic(s: &[RigidBodyState]) -> f64 {
    s.iter()
        .enumerate()
        .map(|(i, s)| {
            s.linear_momentum_world().dot(s.linear_momentum_world()) / (2. * [2., 3., 5.][i])
                + 0.5
                    * s.angular_momentum_body()
                        .dot(omega(s.angular_momentum_body(), i))
        })
        .sum()
}
fn totals(s: &[RigidBodyState], origin: Vec3) -> (Vec3, Vec3) {
    s.iter().fold((Vec3::ZERO, Vec3::ZERO), |(p, h), s| {
        (
            p.add(s.linear_momentum_world()),
            h.add(
                s.pose()
                    .position_world()
                    .sub(origin)
                    .cross(s.linear_momentum_world()),
            )
            .add(
                s.pose()
                    .orientation()
                    .rotate_body_to_world(s.angular_momentum_body()),
            ),
        )
    })
}
fn error(a: RigidBodyState, b: RigidBodyState) -> f64 {
    let mut e = norm(a.pose().position_world().sub(b.pose().position_world()))
        .max(norm(
            a.linear_momentum_world().sub(b.linear_momentum_world()),
        ))
        .max(norm(
            a.angular_momentum_body().sub(b.angular_momentum_body()),
        ));
    for axis in [X, Y, Z] {
        e = e.max(norm(
            a.pose()
                .orientation()
                .rotate_body_to_world(axis)
                .sub(b.pose().orientation().rotate_body_to_world(axis)),
        ));
    }
    e
}
fn impulse(bodies: [usize; 2], force: Vec3) -> InternalWrenchLink {
    InternalWrenchLink {
        bodies,
        impulse: InternalWrenchImpulse {
            arms_body_m: [Vec3::ZERO; 2],
            impulse_on_second_world_ns: force,
            couple_on_second_world_nms: Vec3::ZERO,
        },
    }
}

#[test]
fn g1_shared_body_resultants_retain_cross_energy_and_allocate_simultaneous_source_work() {
    let p = properties();
    let v = [0.4, -0.3, 0.1];
    let s = core::array::from_fn::<_, 3, _>(|i| {
        state(
            X.scale([-1., 0., 2.][i]),
            UnitQuaternion::IDENTITY,
            X.scale(p[i].mass() * v[i]),
            Vec3::ZERO,
        )
    });
    let links = [impulse([0, 1], X.scale(1.2)), impulse([2, 1], X.scale(0.7))];
    let out = apply_internal_wrench_network(&s, &p, &links, controls(), || false).unwrap();
    let dp = [-1.2, 1.9, -0.7];
    let mid = core::array::from_fn::<_, 3, _>(|i| v[i] + 0.5 * dp[i] / p[i].mass());
    let expected = (0..3)
        .map(|i| v[i] * dp[i] + 0.5 * dp[i] * dp[i] / p[i].mass())
        .sum::<f64>();
    near(out.kinetic_energy_change_j, expected, 3e-15);
    near(out.connector_work_j, expected, 3e-15);
    near(
        out.connectors[0].connector_work_j,
        1.2 * (mid[1] - mid[0]),
        3e-15,
    );
    near(
        out.connectors[1].connector_work_j,
        0.7 * (mid[1] - mid[2]),
        3e-15,
    );
    for i in 0..3 {
        assert_eq!(out.bodies[i].state_before, s[i]);
        vector(
            out.bodies[i].state_after.linear_momentum_world(),
            s[i].linear_momentum_world().add(X.scale(dp[i])),
            3e-15,
        );
    }
    let isolated = links
        .iter()
        .map(|l| {
            apply_internal_wrench_impulse(
                l.bodies.map(|i| s[i]),
                l.bodies.map(|i| p[i]),
                l.impulse,
                controls().balance,
                || false,
            )
            .unwrap()
            .connector_work_j
        })
        .sum::<f64>();
    // This is the cross term an independent-pair energy sum would lose.
    near(
        out.connector_work_j - isolated,
        1.2 * 0.7 / p[1].mass(),
        3e-15,
    );
    let reordered =
        apply_internal_wrench_network(&s, &p, &[links[1], links[0]], controls(), || false).unwrap();
    for i in 0..3 {
        assert!(error(out.bodies[i].state_after, reordered.bodies[i].state_after) < 2e-14);
    }
    near(
        out.connectors[0].connector_work_j,
        reordered.connectors[1].connector_work_j,
        3e-15,
    );
    let one =
        apply_internal_wrench_network(&s[..2], &p[..2], &links[..1], controls(), || false).unwrap();
    let pair = apply_internal_wrench_impulse(
        [s[0], s[1]],
        [p[0], p[1]],
        links[0].impulse,
        controls().balance,
        || false,
    )
    .unwrap();
    assert_eq!(one.bodies.as_slice(), &pair.bodies);
    near(one.connector_work_j, pair.connector_work_j, 3e-15);
}

fn initial() -> [RigidBodyState; 3] {
    core::array::from_fn(|i| {
        state(
            [
                Vec3::new(-0.8, 0.3, -0.2),
                Vec3::new(0.4, -0.5, 0.3),
                Vec3::new(0.9, 0.2, -0.4),
            ][i],
            UnitQuaternion::from_axis_angle(
                Vec3::new(1., 2. - i as f64, -0.5),
                0.2 + 0.17 * i as f64,
            )
            .unwrap(),
            Vec3::new(0.2, -0.1, 0.15).scale(1. + 0.3 * i as f64),
            Vec3::new(0.3, -0.2, 0.1).scale(1. + 0.2 * i as f64),
        )
    })
}
fn drives() -> [PrescribedInternalLink; 3] {
    core::array::from_fn(|i| PrescribedInternalLink {
        bodies: [i, (i + 1) % 3],
        drive: PrescribedInternalWrench {
            arms_body_m: [
                Vec3::new(0.2, -0.1, 0.15).scale(1. + i as f64 * 0.1),
                Vec3::new(-0.1, 0.3, 0.2).scale(0.7 + i as f64 * 0.2),
            ],
            force_on_second_world_n: [
                Vec3::new(0.6, -0.4, 0.3),
                Vec3::new(-0.2, 0.5, 0.1),
                Vec3::new(0.3, 0.2, -0.4),
            ][i],
            couple_on_second_world_nm: [
                Vec3::new(0.1, 0.2, -0.15),
                Vec3::new(-0.2, 0.1, 0.05),
                Vec3::new(0.15, -0.1, 0.2),
            ][i],
        },
    })
}
fn impulses(d: &[PrescribedInternalLink; 3], h: f64) -> [InternalWrenchLink; 3] {
    d.map(|l| InternalWrenchLink {
        bodies: l.bodies,
        impulse: InternalWrenchImpulse {
            arms_body_m: l.drive.arms_body_m,
            impulse_on_second_world_ns: l.drive.force_on_second_world_n.scale(h),
            couple_on_second_world_nms: l.drive.couple_on_second_world_nm.scale(h),
        },
    })
}
fn compose(a: UnitQuaternion, b: UnitQuaternion) -> UnitQuaternion {
    let [w, x, y, z] = a.components();
    let [s, u, v, t] = b.components();
    UnitQuaternion::new(
        w * s - x * u - y * v - z * t,
        w * u + x * s + y * t - z * v,
        w * v - x * t + y * s + z * u,
        w * t + x * v - y * u + z * s,
    )
    .unwrap()
}

#[test]
fn g1_g3_full_tensor_cycle_matches_spatial_wrench_algebra_and_frame_boost_covariance() {
    let p = properties();
    let s = initial();
    let links = impulses(&drives(), 0.2);
    let out = apply_internal_wrench_network(&s, &p, &links, controls(), || false).unwrap();
    let mut force = [Vec3::ZERO; 3];
    let mut torque = [Vec3::ZERO; 3];
    for l in links {
        let [a, b] = l.bodies;
        let j = l.impulse.impulse_on_second_world_ns;
        let rb = s[b]
            .pose()
            .orientation()
            .rotate_body_to_world(l.impulse.arms_body_m[1]);
        // Independent spatial transport: receiver central wrench (J, C+rB×J)
        // fixes donor central torque -(xB-xA+rB)×J-C, whatever A's reference arm.
        force[a] = force[a].sub(j);
        force[b] = force[b].add(j);
        torque[a] = torque[a]
            .sub(
                s[b].pose()
                    .position_world()
                    .sub(s[a].pose().position_world())
                    .add(rb)
                    .cross(j),
            )
            .sub(l.impulse.couple_on_second_world_nms);
        torque[b] = torque[b]
            .add(rb.cross(j))
            .add(l.impulse.couple_on_second_world_nms);
    }
    let after: Vec<_> = out.bodies.iter().map(|r| r.state_after).collect();
    for i in 0..3 {
        vector(
            after[i].linear_momentum_world(),
            s[i].linear_momentum_world().add(force[i]),
            2e-13,
        );
        vector(
            after[i].angular_momentum_body(),
            s[i].angular_momentum_body()
                .add(s[i].pose().orientation().rotate_world_to_body(torque[i])),
            2e-13,
        );
    }
    let before = totals(&s, controls().balance.angular_origin_world_m);
    let total = totals(&after, controls().balance.angular_origin_world_m);
    vector(total.0, before.0, 2e-13);
    vector(total.1, before.1, 3e-13);
    near(kinetic(&after) - kinetic(&s), out.connector_work_j, 3e-13);
    let q = UnitQuaternion::from_axis_angle(Vec3::new(1., 2., -1.), 0.7).unwrap();
    let shift = Vec3::new(3., -2., 0.4);
    let boost = Vec3::new(0.5, -0.4, 0.2);
    let moved = core::array::from_fn::<_, 3, _>(|i| {
        state(
            q.rotate_body_to_world(s[i].pose().position_world())
                .add(shift),
            compose(q, s[i].pose().orientation()),
            q.rotate_body_to_world(s[i].linear_momentum_world())
                .add(boost.scale(p[i].mass())),
            s[i].angular_momentum_body(),
        )
    });
    let transformed = links.map(|l| InternalWrenchLink {
        bodies: l.bodies,
        impulse: InternalWrenchImpulse {
            arms_body_m: l.impulse.arms_body_m,
            impulse_on_second_world_ns: q
                .rotate_body_to_world(l.impulse.impulse_on_second_world_ns),
            couple_on_second_world_nms: q
                .rotate_body_to_world(l.impulse.couple_on_second_world_nms),
        },
    });
    let mut c = controls();
    c.balance.angular_origin_world_m = q
        .rotate_body_to_world(c.balance.angular_origin_world_m)
        .add(shift);
    let moved = apply_internal_wrench_network(&moved, &p, &transformed, c, || false).unwrap();
    for i in 0..3 {
        vector(
            moved.bodies[i].state_after.linear_momentum_world(),
            q.rotate_body_to_world(after[i].linear_momentum_world())
                .add(boost.scale(p[i].mass())),
            2e-13,
        );
        vector(
            moved.bodies[i].state_after.angular_momentum_body(),
            after[i].angular_momentum_body(),
            2e-13,
        );
        near(
            moved.connectors[i].connector_work_j,
            out.connectors[i].connector_work_j,
            3e-13,
        );
    }
    // Relabeling bodies changes no physical result or source attribution.
    let order = [2, 0, 1];
    let inverse = [1, 2, 0];
    let permuted = links.map(|l| InternalWrenchLink {
        bodies: l.bodies.map(|i| inverse[i]),
        ..l
    });
    let reordered = apply_internal_wrench_network(
        &order.map(|i| s[i]),
        &order.map(|i| p[i]),
        &permuted,
        controls(),
        || false,
    )
    .unwrap();
    for i in 0..3 {
        assert!(error(reordered.bodies[inverse[i]].state_after, after[i]) < 2e-13);
    }
}

#[test]
fn g1_finite_actuator_chain_has_exact_unequal_mass_recoil_torque_and_source_work() {
    let p = properties();
    let x = [-1., 0., 2.];
    let v = [0.2, -0.1, 0.3];
    let w = [0.1, -0.2, 0.05];
    let angles = [0.2, -0.1, 0.3];
    let s = core::array::from_fn::<_, 3, _>(|i| {
        state(
            X.scale(x[i]),
            UnitQuaternion::from_axis_angle(Z, angles[i]).unwrap(),
            X.scale(p[i].mass() * v[i]),
            Z.scale((4. + i as f64) * w[i]),
        )
    });
    let f = [0.8, -0.3];
    let c = [0.25, -0.15];
    let links = core::array::from_fn::<_, 2, _>(|i| PrescribedInternalLink {
        bodies: [i, i + 1],
        drive: PrescribedInternalWrench {
            arms_body_m: [Vec3::ZERO; 2],
            force_on_second_world_n: X.scale(f[i]),
            couple_on_second_world_nm: Z.scale(c[i]),
        },
    });
    let force = [-f[0], f[0] - f[1], f[1]];
    let torque = [-c[0], c[0] - c[1], c[1]];
    let time = 0.3;
    let expected = core::array::from_fn::<_, 3, _>(|i| {
        state(
            X.scale(x[i] + v[i] * time + 0.5 * force[i] / p[i].mass() * time * time),
            UnitQuaternion::from_axis_angle(
                Z,
                angles[i] + w[i] * time + 0.5 * torque[i] / (4. + i as f64) * time * time,
            )
            .unwrap(),
            X.scale(p[i].mass() * v[i] + time * force[i]),
            Z.scale((4. + i as f64) * w[i] + time * torque[i]),
        )
    });
    for count in [1, 6] {
        let h = time / f64::from(count);
        let mut current = s.to_vec();
        let mut works = [0.; 2];
        for _ in 0..count {
            let out = step_internal_wrench_network(&current, &p, &links, h, controls(), || false)
                .unwrap();
            for i in 0..3 {
                assert_eq!(
                    out.free_steps[i],
                    RigidBodyIntegrator::new(Gravity::ZERO)
                        .step(
                            out.half_impulses[0].bodies[i].state_after,
                            p[i],
                            Wrench::ZERO,
                            h
                        )
                        .unwrap()
                );
            }
            for stage in &out.half_impulses {
                for (i, r) in stage.connectors.iter().enumerate() {
                    works[i] += r.connector_work_j;
                }
            }
            current = out.states_after;
        }
        for i in 0..3 {
            assert!(error(current[i], expected[i]) < 3e-13);
        }
        for i in 0..2 {
            let dx = core::array::from_fn::<_, 3, _>(|j| {
                v[j] * time + 0.5 * force[j] / p[j].mass() * time * time
            });
            let da = core::array::from_fn::<_, 3, _>(|j| {
                w[j] * time + 0.5 * torque[j] / (4. + j as f64) * time * time
            });
            near(
                works[i],
                f[i] * (dx[i + 1] - dx[i]) + c[i] * (da[i + 1] - da[i]),
                3e-13,
            );
        }
        near(works.iter().sum(), kinetic(&expected) - kinetic(&s), 3e-13);
    }
    let pair = step_internal_wrench(
        [s[0], s[1]],
        [p[0], p[1]],
        links[0].drive,
        time,
        controls().balance,
        || false,
    )
    .unwrap();
    let one =
        step_internal_wrench_network(&s[..2], &p[..2], &links[..1], time, controls(), || false)
            .unwrap();
    for i in 0..2 {
        assert!(error(one.states_after[i], pair.states_after[i]) < 2e-14);
    }
    near(one.actuator_work_j, pair.actuator_work_j, 3e-14);
}

// Independent continuous 3-body Newton/Euler ODE, plus each connector's physical
// work integral. It does not call native impulse, inertia inversion or stepping.
type Ode = [f64; 42];
fn triple(y: &Ode, i: usize) -> Vec3 {
    Vec3::new(y[i], y[i + 1], y[i + 2])
}
fn put(y: &mut Ode, i: usize, v: Vec3) {
    y[i] = v.x;
    y[i + 1] = v.y;
    y[i + 2] = v.z;
}
fn rotate(q: [f64; 4], v: Vec3) -> Vec3 {
    let length = q.iter().map(|v| v * v).sum::<f64>().sqrt();
    let [w, x, y, z] = q.map(|v| v / length);
    let u = Vec3::new(x, y, z);
    v.add(u.cross(u.cross(v).add(v.scale(w))).scale(2.))
}
fn rhs(y: Ode, links: &[PrescribedInternalLink; 3]) -> Ode {
    let mut out = [0.; 42];
    let mut f = [Vec3::ZERO; 3];
    let mut tau = [Vec3::ZERO; 3];
    let q =
        core::array::from_fn::<_, 3, _>(|i| core::array::from_fn::<_, 4, _>(|j| y[13 * i + 6 + j]));
    let spin = core::array::from_fn::<_, 3, _>(|i| omega(triple(&y, 13 * i + 10), i));
    for (k, l) in links.iter().enumerate() {
        let [a, b] = l.bodies;
        let force = l.drive.force_on_second_world_n;
        let rb = rotate(q[b], l.drive.arms_body_m[1]);
        let cb = l.drive.couple_on_second_world_nm.add(rb.cross(force));
        let ca = l.drive.couple_on_second_world_nm.scale(-1.).sub(
            triple(&y, 13 * b)
                .sub(triple(&y, 13 * a))
                .add(rb)
                .cross(force),
        );
        f[a] = f[a].sub(force);
        f[b] = f[b].add(force);
        tau[a] = tau[a].add(ca);
        tau[b] = tau[b].add(cb);
        out[39 + k] = force.dot(
            triple(&y, 13 * b + 3)
                .scale(1. / [2., 3., 5.][b])
                .sub(triple(&y, 13 * a + 3).scale(1. / [2., 3., 5.][a])),
        ) + ca.dot(rotate(q[a], spin[a]))
            + cb.dot(rotate(q[b], spin[b]));
    }
    for i in 0..3 {
        let n = 13 * i;
        let [s, qx, qy, qz] = q[i];
        let w = spin[i];
        put(&mut out, n, triple(&y, n + 3).scale(1. / [2., 3., 5.][i]));
        put(&mut out, n + 3, f[i]);
        out[n + 6] = -0.5 * (qx * w.x + qy * w.y + qz * w.z);
        out[n + 7] = 0.5 * (s * w.x + qy * w.z - qz * w.y);
        out[n + 8] = 0.5 * (s * w.y + qz * w.x - qx * w.z);
        out[n + 9] = 0.5 * (s * w.z + qx * w.y - qy * w.x);
        put(
            &mut out,
            n + 10,
            rotate([s, -qx, -qy, -qz], tau[i]).sub(w.cross(triple(&y, n + 10))),
        );
    }
    out
}
fn ode_initial(s: [RigidBodyState; 3]) -> Ode {
    let mut y = [0.; 42];
    for i in 0..3 {
        put(&mut y, 13 * i, s[i].pose().position_world());
        put(&mut y, 13 * i + 3, s[i].linear_momentum_world());
        for (j, x) in s[i]
            .pose()
            .orientation()
            .components()
            .into_iter()
            .enumerate()
        {
            y[13 * i + 6 + j] = x;
        }
        put(&mut y, 13 * i + 10, s[i].angular_momentum_body());
    }
    y
}
fn reverse(s: RigidBodyState) -> RigidBodyState {
    state(
        s.pose().position_world(),
        s.pose().orientation(),
        s.linear_momentum_world().scale(-1.),
        s.angular_momentum_body().scale(-1.),
    )
}

#[test]
fn g1_rotating_full_tensor_cycle_refines_to_independent_continuous_trajectory_and_work() {
    let p = properties();
    let initial = initial();
    let links = drives();
    let time = 0.6;
    let integrate = |count| {
        let mut y = ode_initial(initial);
        let h = time / f64::from(count);
        for _ in 0..count {
            let a = rhs(y, &links);
            let b = rhs(core::array::from_fn(|i| y[i] + h * a[i] / 2.), &links);
            let c = rhs(core::array::from_fn(|i| y[i] + h * b[i] / 2.), &links);
            let d = rhs(core::array::from_fn(|i| y[i] + h * c[i]), &links);
            y = core::array::from_fn(|i| y[i] + h * (a[i] + 2. * b[i] + 2. * c[i] + d[i]) / 6.);
        }
        y
    };
    let reference = integrate(4096);
    let fine = integrate(8192);
    for i in 0..42 {
        near(reference[i], fine[i], 2e-12);
    }
    let expected = core::array::from_fn::<_, 3, _>(|i| {
        state(
            triple(&fine, 13 * i),
            UnitQuaternion::new(
                fine[13 * i + 6],
                fine[13 * i + 7],
                fine[13 * i + 8],
                fine[13 * i + 9],
            )
            .unwrap(),
            triple(&fine, 13 * i + 3),
            triple(&fine, 13 * i + 10),
        )
    });
    near(
        fine[39..].iter().sum(),
        kinetic(&expected) - kinetic(&initial),
        4e-12,
    );
    let before = totals(&initial, controls().balance.angular_origin_world_m);
    let mut errors = Vec::new();
    let mut work_errors = Vec::new();
    for count in [20, 40, 80] {
        let h = time / f64::from(count);
        let mut s = initial.to_vec();
        let mut work = [0.; 3];
        for _ in 0..count {
            let out =
                step_internal_wrench_network(&s, &p, &links, h, controls(), || false).unwrap();
            for stage in &out.half_impulses {
                for (i, r) in stage.connectors.iter().enumerate() {
                    work[i] += r.connector_work_j;
                }
            }
            s = out.states_after;
        }
        errors.push(
            s.iter()
                .zip(expected)
                .map(|(&a, b)| error(a, b))
                .fold(0., f64::max),
        );
        work_errors.push(
            (0..3)
                .map(|i| (work[i] - fine[39 + i]).abs())
                .fold(0., f64::max),
        );
        let after = totals(&s, controls().balance.angular_origin_world_m);
        vector(after.0, before.0, 2e-11);
        vector(after.1, before.1, 3e-10);
        near(kinetic(&s) - kinetic(&initial), work.iter().sum(), 3e-10);
        s = s.into_iter().map(reverse).collect();
        for _ in 0..count {
            s = step_internal_wrench_network(&s, &p, &links, h, controls(), || false)
                .unwrap()
                .states_after;
        }
        for i in 0..3 {
            assert!(error(reverse(s[i]), initial[i]) < 3e-9);
        }
    }
    assert!(errors[2] < 3e-5, "{errors:?}");
    for series in [&errors, &work_errors] {
        for pair in series.windows(2) {
            assert!((3.6..4.4).contains(&(pair[0] / pair[1])), "{series:?}");
        }
    }
}

#[test]
fn g4_whole_network_budgets_swallowed_updates_and_every_cancellation_refuse_atomically() {
    let s = initial();
    let p = properties();
    let drives = drives();
    let links = impulses(&drives, 0.1);
    let mut polls = 0;
    let out = apply_internal_wrench_network(&s, &p, &links, controls(), || {
        polls += 1;
        false
    })
    .unwrap();
    let mut exact = controls();
    exact.maximum_live_bytes = out.planned_live_bytes;
    assert_eq!(
        out,
        apply_internal_wrench_network(&s, &p, &links, exact, || false).unwrap()
    );
    exact.maximum_live_bytes -= 1;
    assert!(matches!(
        apply_internal_wrench_network(&s, &p, &links, exact, || false),
        Err(InternalWrenchNetworkError::Limit { .. })
    ));
    for stop in 1..=polls {
        let mut n = 0;
        assert!(matches!(
            apply_internal_wrench_network(&s, &p, &links, controls(), || {
                n += 1;
                n == stop
            }),
            Err(InternalWrenchNetworkError::Cancelled)
        ));
    }
    assert_eq!(
        out,
        apply_internal_wrench_network(&s, &p, &links, controls(), || false).unwrap()
    );
    let mut polls = 0;
    let out = step_internal_wrench_network(&s, &p, &drives, 0.05, controls(), || {
        polls += 1;
        false
    })
    .unwrap();
    let mut exact = controls();
    exact.maximum_live_bytes = out.planned_live_bytes;
    assert_eq!(
        out,
        step_internal_wrench_network(&s, &p, &drives, 0.05, exact, || false).unwrap()
    );
    exact.maximum_live_bytes -= 1;
    assert!(matches!(
        step_internal_wrench_network(&s, &p, &drives, 0.05, exact, || false),
        Err(InternalWrenchNetworkStepError::Network(
            InternalWrenchNetworkError::Limit { .. }
        ))
    ));
    for stop in 1..=polls {
        let mut n = 0;
        assert!(
            matches!(
                step_internal_wrench_network(&s, &p, &drives, 0.05, controls(), || {
                    n += 1;
                    n == stop
                }),
                Err(InternalWrenchNetworkStepError::Cancelled)
            ),
            "stop={stop}"
        );
    }
    assert_eq!(
        out,
        step_internal_wrench_network(&s, &p, &drives, 0.05, controls(), || false).unwrap()
    );
    let mut count = controls();
    count.maximum_bodies = 2;
    assert!(apply_internal_wrench_network(&s, &p, &links, count, || false).is_err());
    count = controls();
    count.maximum_connectors = 2;
    assert!(apply_internal_wrench_network(&s, &p, &links, count, || false).is_err());
    let mut invalid = links;
    invalid[2].bodies = [2, 2];
    assert!(apply_internal_wrench_network(&s, &p, &invalid, controls(), || false).is_err());
    let mut late = drives;
    late[2].drive.force_on_second_world_n.x = f64::INFINITY;
    assert!(step_internal_wrench_network(&s, &p, &late, 0.1, controls(), || false).is_err());
    assert_eq!(s, initial());
    let empty = step_internal_wrench_network(&[], &[], &[], 0.1, controls(), || false).unwrap();
    assert_eq!(empty.planned_live_bytes, 0);
    assert!(empty.states_after.is_empty());
    let free = step_internal_wrench_network(&s, &p, &[], 0.05, controls(), || false).unwrap();
    for i in 0..3 {
        assert_eq!(
            free.states_after[i],
            RigidBodyIntegrator::new(Gravity::ZERO)
                .step(s[i], p[i], Wrench::ZERO, 0.05)
                .unwrap()
                .state_after
        );
    }
    let huge = [
        state(
            X.scale(-1.),
            UnitQuaternion::IDENTITY,
            X.scale(1e20),
            Vec3::ZERO,
        ),
        state(X, UnitQuaternion::IDENTITY, X.scale(1e20), Vec3::ZERO),
    ];
    let unit = MassProperties::new(1., Vec3::ZERO, Vec3::new(1., 1., 1.)).unwrap();
    assert!(matches!(
        apply_internal_wrench_network(&huge, &[unit; 2], &[impulse([0, 1], X)], controls(), || {
            false
        }),
        Err(InternalWrenchNetworkError::Impulse(
            InternalWrenchImpulseError::Unresolved { .. }
        ))
    ));
    // Exactly opposed simultaneous sources on each body cancel physically;
    // they are distinct from two undelivered opposite BODY resultants.
    let cancelled = apply_internal_wrench_network(
        &huge,
        &[unit; 2],
        &[impulse([0, 1], X), impulse([0, 1], X.scale(-1.))],
        controls(),
        || false,
    )
    .unwrap();
    for i in 0..2 {
        assert_eq!(cancelled.bodies[i].state_after, huge[i]);
    }
}
