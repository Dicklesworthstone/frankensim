//! Physical graph sources: shared-body recoil, stiffness, null modes and refinement.

use fs_mbd::spatial_loads::internal_impulse::{
    InternalWrenchImpulseControls,
    network::InternalWrenchNetworkError,
    viscous::{
        BushingCoordinateScales, PassiveBushingMatrix, PreparedViscousConnector,
        ViscousConnectorControls,
        network::{
            PreparedViscousNetwork, ViscousNetworkError, ViscousNetworkLimits, ViscousNetworkLink,
        },
    },
};
use fs_mbd::{MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3};

const X: Vec3 = Vec3::new(1., 0., 0.);
const Y: Vec3 = Vec3::new(0., 1., 0.);
const Z: Vec3 = Vec3::new(0., 0., 1.);
fn a(v: Vec3) -> [f64; 3] {
    [v.x, v.y, v.z]
}
fn v(a: [f64; 3]) -> Vec3 {
    Vec3::new(a[0], a[1], a[2])
}
fn near(a: f64, b: f64, tol: f64) {
    assert!((a - b).abs() <= tol, "{a:e} != {b:e}, tolerance {tol:e}");
}
fn vector(a: Vec3, b: Vec3, tol: f64) {
    for (a, b) in a_components(a).into_iter().zip(a_components(b)) {
        near(a, b, tol);
    }
}
fn a_components(v: Vec3) -> [f64; 3] {
    [v.x, v.y, v.z]
}
fn controls() -> ViscousConnectorControls {
    ViscousConnectorControls {
        impulse: InternalWrenchImpulseControls {
            angular_origin_world_m: Vec3::ZERO,
            linear_tolerance_ns: 3e-10,
            angular_tolerance_nms: 3e-10,
            work_tolerance_j: 3e-10,
        },
        linear_speed_tolerance_m_s: 3e-10,
        angular_speed_tolerance_rad_s: 3e-10,
    }
}
fn limits() -> ViscousNetworkLimits {
    ViscousNetworkLimits {
        maximum_bodies: 32,
        maximum_connectors: 64,
        maximum_pair_relaxations: 100_000,
        maximum_live_bytes: 1_000_000,
    }
}
fn law(c: [[f64; 6]; 6]) -> PassiveBushingMatrix {
    PassiveBushingMatrix::try_new(
        c,
        BushingCoordinateScales {
            translation_m: 1.,
            rotation_rad: 1.,
        },
    )
    .unwrap()
}
fn diagonal(c: [f64; 6]) -> PassiveBushingMatrix {
    law(core::array::from_fn(|i| {
        core::array::from_fn(|j| if i == j { c[i] } else { 0. })
    }))
}
fn coupled() -> PassiveBushingMatrix {
    let d = [0.6, 0.8, 0.5, 0.2, 0.3, 0.4];
    let b = [0.1, -0.2, 0.05, 0.02, 0.03, -0.04];
    law(core::array::from_fn(|i| {
        core::array::from_fn(|j| b[i] * b[j] + if i == j { d[i] } else { 0. })
    }))
}
fn state(x: Vec3, q: UnitQuaternion, p: Vec3, h: Vec3) -> RigidBodyState {
    RigidBodyState::new(Pose::new(x, q).unwrap(), p, h).unwrap()
}
fn line_states(velocity: [f64; 3]) -> [RigidBodyState; 3] {
    core::array::from_fn(|i| {
        state(
            X.scale(i as f64 - 1.),
            UnitQuaternion::IDENTITY,
            X.scale(2. * velocity[i]),
            Vec3::ZERO,
        )
    })
}
fn line_properties() -> [MassProperties; 3] {
    [MassProperties::new(2., Vec3::ZERO, Vec3::new(2., 3., 4.)).unwrap(); 3]
}
fn axial_link(bodies: [usize; 2], coefficient: f64) -> ViscousNetworkLink {
    ViscousNetworkLink {
        bodies,
        arms_body_m: [Vec3::ZERO; 2],
        damping: diagonal([coefficient, 0., 0., 0., 0., 0.]),
    }
}

#[test]
fn one_active_connector_reuses_the_exact_native_pair_and_original_decay_rate() {
    let properties = [
        MassProperties::new(2., Vec3::ZERO, Vec3::new(2., 3., 4.)).unwrap(),
        MassProperties::new(3., Vec3::ZERO, Vec3::new(3., 4., 5.)).unwrap(),
    ];
    let states = [
        state(
            Vec3::ZERO,
            UnitQuaternion::IDENTITY,
            X.scale(0.8),
            Vec3::ZERO,
        ),
        state(
            Vec3::ZERO,
            UnitQuaternion::IDENTITY,
            X.scale(-0.6),
            Vec3::ZERO,
        ),
    ];
    let link = axial_link([0, 1], 6.);
    let pair = PreparedViscousConnector::new(
        states.map(|s| s.pose()),
        properties,
        link.arms_body_m,
        link.damping.clone(),
        0.3,
        || false,
    )
    .unwrap();
    near(pair.maximum_decay_rate_per_s(), 5., 2e-14);
    let expected = pair.relax(states, controls(), || false).unwrap();
    let network = PreparedViscousNetwork::new(
        &states.map(|s| s.pose()),
        &properties,
        &[link.clone()],
        0.3,
        1,
        limits(),
        || false,
    )
    .unwrap();
    let actual = network.relax(&states, controls(), || false).unwrap();
    assert_eq!(
        actual.states_after,
        expected.impulse.bodies.map(|b| b.state_after)
    );
    assert_eq!(
        actual.connector_dissipation_j,
        [expected.dissipated_energy_j]
    );
    assert_eq!(actual.dissipated_energy_j, expected.dissipated_energy_j);
    assert_eq!(actual.pair_relaxations, 1);
    near(actual.kinetic_heat_residual_j, 0., 2e-15);
    near(actual.connector_heat_residual_j, 0., 2e-15);
    let no_time = PreparedViscousConnector::new(
        states.map(|s| s.pose()),
        properties,
        link.arms_body_m,
        link.damping,
        0.,
        || false,
    )
    .unwrap();
    // Zero means preparation was intentionally skipped, not a physical rate
    // claim about this nonzero constitutive matrix.
    assert_eq!(no_time.maximum_decay_rate_per_s(), 0.);
    let no_law = PreparedViscousConnector::new(
        states.map(|s| s.pose()),
        properties,
        [Vec3::ZERO; 2],
        diagonal([0.; 6]),
        1.,
        || false,
    )
    .unwrap();
    assert_eq!(no_law.maximum_decay_rate_per_s(), 0.);
}

#[test]
fn shared_three_body_chain_refines_quadratically_to_independent_graph_modes() {
    let states = line_states([1.2, -0.8, 0.2]);
    let properties = line_properties();
    let links = [axial_link([0, 1], 0.6), axial_link([1, 2], 0.6)];
    let duration = 0.8;
    // Equal masses: the path-graph Laplacian has modes [1,0,-1] and
    // [1,-2,1], with rates c/m and 3c/m. No production eigensolver here.
    let first = 0.5_f64 * (-0.3_f64 * duration).exp();
    let second = 0.5_f64 * (-0.9_f64 * duration).exp();
    let expected = [
        0.2 + first + second,
        0.2 - 2. * second,
        0.2 - first + second,
    ];
    let heat =
        2. * (0.25 * -(-0.6_f64 * duration).exp_m1() + 0.75 * -(-1.8_f64 * duration).exp_m1());
    let mut errors = Vec::new();
    for n in [4, 8, 16] {
        let prepared = PreparedViscousNetwork::new(
            &states.map(|s| s.pose()),
            &properties,
            &links,
            duration,
            n,
            limits(),
            || false,
        )
        .unwrap();
        let out = prepared.relax(&states, controls(), || false).unwrap();
        let error = out
            .states_after
            .iter()
            .enumerate()
            .map(|(i, s)| {
                (s.center_of_mass_velocity_world(properties[i]).unwrap().x - expected[i]).abs()
            })
            .fold(0., f64::max);
        assert!(error > 1e-9);
        errors.push(error);
        assert_eq!(out.pair_relaxations, 3 * n);
        assert_eq!(prepared.active_connectors(), 2);
        assert_eq!(prepared.pair_relaxations(), 3 * n);
        vector(out.linear_momentum_change_world_ns, Vec3::ZERO, 2e-14);
        vector(out.angular_momentum_change_world_nms, Vec3::ZERO, 2e-14);
        near(out.kinetic_heat_residual_j, 0., 2e-13);
        near(
            out.connector_dissipation_j.iter().sum(),
            out.dissipated_energy_j,
            2e-15,
        );
        near(out.dissipated_energy_j, heat, 10. * error);
    }
    for pair in errors.windows(2) {
        assert!((3.8..4.2).contains(&(pair[0] / pair[1])), "{errors:?}");
    }
}

fn full_properties() -> [MassProperties; 3] {
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
fn full_states() -> [RigidBodyState; 3] {
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
fn cycle_links() -> [ViscousNetworkLink; 3] {
    core::array::from_fn(|i| ViscousNetworkLink {
        bodies: [i, (i + 1) % 3],
        arms_body_m: [
            Vec3::new(0.2, -0.1, 0.15).scale(1. + 0.1 * i as f64),
            Vec3::new(-0.1, 0.3, 0.2).scale(0.7 + 0.2 * i as f64),
        ],
        damping: coupled(),
    })
}
// Explicit inverse of I=[[a,-.5,0],[-.5,a,0],[0,0,4+i]] from the stated
// principal tensor rotated 45 degrees. Independent of native inertia actions.
fn omega(h: Vec3, i: usize) -> Vec3 {
    let x = 2.5 + i as f64;
    let det = (2. + i as f64) * (3. + i as f64);
    Vec3::new(
        (x * h.x + 0.5 * h.y) / det,
        (0.5 * h.x + x * h.y) / det,
        h.z / (4. + i as f64),
    )
}
fn initial_ode(states: &[RigidBodyState; 3]) -> [f64; 19] {
    let mut y = [0.; 19];
    for i in 0..3 {
        y[6 * i..6 * i + 3].copy_from_slice(&a(states[i].linear_momentum_world()));
        y[6 * i + 3..6 * i + 6].copy_from_slice(&a(states[i].angular_momentum_body()));
    }
    y
}
fn rhs(y: [f64; 19], poses: &[Pose; 3], links: &[ViscousNetworkLink; 3]) -> [f64; 19] {
    let velocity: [Vec3; 3] = core::array::from_fn(|i| {
        v([y[6 * i], y[6 * i + 1], y[6 * i + 2]]).scale(1. / [2., 3., 5.][i])
    });
    let spin: [Vec3; 3] = core::array::from_fn(|i| {
        poses[i]
            .orientation()
            .rotate_body_to_world(omega(v([y[6 * i + 3], y[6 * i + 4], y[6 * i + 5]]), i))
    });
    let mut out = [0.; 19];
    for link in links {
        let [a, b] = link.bodies;
        let qa = poses[a].orientation();
        let rb = poses[b]
            .orientation()
            .rotate_body_to_world(link.arms_body_m[1]);
        let separation = poses[b]
            .position_world()
            .sub(poses[a].position_world())
            .add(rb);
        // Transport A's velocity to B's material attachment. Its own reference
        // attachment arm cancels from this spatial twist/wrench pairing.
        let translation = qa.rotate_world_to_body(
            velocity[b]
                .sub(velocity[a])
                .add(spin[b].cross(rb))
                .sub(spin[a].cross(separation)),
        );
        let rotation = qa.rotate_world_to_body(spin[b].sub(spin[a]));
        let slip = [
            translation.x,
            translation.y,
            translation.z,
            rotation.x,
            rotation.y,
            rotation.z,
        ];
        let c = link.damping.native_coefficients();
        let effort: [f64; 6] =
            core::array::from_fn(|i| -(0..6).map(|j| c[i][j] * slip[j]).sum::<f64>());
        let force = qa.rotate_body_to_world(v([effort[0], effort[1], effort[2]]));
        let couple = qa.rotate_body_to_world(v([effort[3], effort[4], effort[5]]));
        let moments = [
            separation.cross(force).add(couple).scale(-1.),
            rb.cross(force).add(couple),
        ];
        for (side, body) in [a, b].into_iter().enumerate() {
            let torque = poses[body]
                .orientation()
                .rotate_world_to_body(moments[side]);
            for axis in 0..3 {
                out[6 * body + axis] +=
                    a_components(force)[axis] * if side == 0 { -1. } else { 1. };
                out[6 * body + 3 + axis] += a_components(torque)[axis];
            }
        }
        out[18] -= (0..6).map(|i| effort[i] * slip[i]).sum::<f64>();
    }
    out
}
fn reference(
    states: &[RigidBodyState; 3],
    links: &[ViscousNetworkLink; 3],
    duration: f64,
) -> [f64; 19] {
    let poses = states.map(|s| s.pose());
    let mut y = initial_ode(states);
    let dt = duration / 4096.;
    for _ in 0..4096 {
        let k1 = rhs(y, &poses, links);
        let k2 = rhs(
            core::array::from_fn(|i| y[i] + 0.5 * dt * k1[i]),
            &poses,
            links,
        );
        let k3 = rhs(
            core::array::from_fn(|i| y[i] + 0.5 * dt * k2[i]),
            &poses,
            links,
        );
        let k4 = rhs(core::array::from_fn(|i| y[i] + dt * k3[i]), &poses, links);
        for i in 0..19 {
            y[i] += dt * (k1[i] + 2. * k2[i] + 2. * k3[i] + k4[i]) / 6.;
        }
    }
    y
}
fn totals(states: &[RigidBodyState], origin: Vec3) -> (Vec3, Vec3) {
    states.iter().fold((Vec3::ZERO, Vec3::ZERO), |(p, h), s| {
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

#[test]
fn full_tensor_shared_cycle_converges_to_independent_force_torque_and_heat_equations() {
    let states = full_states();
    let properties = full_properties();
    let links = cycle_links();
    let expected = reference(&states, &links, 0.3);
    let mut errors = Vec::new();
    let before = totals(&states, Vec3::ZERO);
    for n in [4, 8, 16] {
        let plan = PreparedViscousNetwork::new(
            &states.map(|s| s.pose()),
            &properties,
            &links,
            0.3,
            n,
            limits(),
            || false,
        )
        .unwrap();
        let out = plan.relax(&states, controls(), || false).unwrap();
        let mut error = (out.dissipated_energy_j - expected[18]).abs();
        for i in 0..3 {
            assert_eq!(out.states_after[i].pose(), states[i].pose());
            for (axis, value) in a(out.states_after[i].linear_momentum_world())
                .into_iter()
                .enumerate()
            {
                error = error.max((value - expected[6 * i + axis]).abs());
            }
            for (axis, value) in a(out.states_after[i].angular_momentum_body())
                .into_iter()
                .enumerate()
            {
                error = error.max((value - expected[6 * i + 3 + axis]).abs());
            }
        }
        errors.push(error);
        let after = totals(&out.states_after, Vec3::ZERO);
        vector(after.0, before.0, 3e-13);
        vector(after.1, before.1, 3e-13);
        near(out.kinetic_heat_residual_j, 0., 3e-13);
        near(out.connector_heat_residual_j, 0., 3e-13);
        assert!(out.connector_dissipation_j.iter().all(|h| *h > 0.));
        assert_eq!(out.pair_relaxations, 5 * n);
    }
    assert!(errors[0] > 1e-8 && errors[2] < 2e-5, "{errors:?}");
    for pair in errors.windows(2) {
        assert!((3.7..4.3).contains(&(pair[0] / pair[1])), "{errors:?}");
    }
}

#[test]
fn common_rigid_motion_and_isolated_body_are_exact_free_modes() {
    let masses = [1., 2., 4., 8.];
    let omega = Z.scale(0.25);
    let translation = Vec3::new(0.5, 0.125, 0.25);
    let properties =
        masses.map(|m| MassProperties::new(m, Vec3::ZERO, Vec3::new(4., 8., 8.)).unwrap());
    let positions = [X.scale(-2.), Y.scale(4.), X.scale(2.), Z.scale(7.)];
    let states = core::array::from_fn::<_, 4, _>(|i| {
        state(
            positions[i],
            UnitQuaternion::IDENTITY,
            translation.add(omega.cross(positions[i])).scale(masses[i]),
            Z.scale(2.),
        )
    });
    let links = core::array::from_fn::<_, 3, _>(|i| ViscousNetworkLink {
        bodies: [i, (i + 1) % 3],
        arms_body_m: [Vec3::new(0.125, -0.25, 0.), Vec3::new(0.25, 0.125, 0.)],
        damping: diagonal([0.5, 1., 0.25, 0.25, 0.5, 1.]),
    });
    let prepared = PreparedViscousNetwork::new(
        &states.map(|s| s.pose()),
        &properties,
        &links,
        1e200,
        2,
        limits(),
        || false,
    )
    .unwrap();
    let out = prepared.relax(&states, controls(), || false).unwrap();
    assert_eq!(out.states_after, states);
    assert_eq!(out.dissipated_energy_j, 0.);
    assert_eq!(out.connector_dissipation_j, [0.; 3]);
    assert_eq!(out.kinetic_energy_change_j, 0.);
}

#[test]
fn shifted_world_and_common_boost_preserve_the_same_spatial_exchange() {
    let properties = [MassProperties::new(2., Vec3::ZERO, Vec3::new(4., 8., 8.)).unwrap(); 3];
    let positions = [
        Vec3::new(-8., 0., 4.),
        Vec3::new(0., 8., -4.),
        Vec3::new(8., -8., 0.),
    ];
    let states = core::array::from_fn::<_, 3, _>(|i| {
        state(
            positions[i],
            UnitQuaternion::IDENTITY,
            Vec3::new(0.125, -0.25, 0.5).scale(1. + i as f64),
            Vec3::new(-0.25, 0.125, 0.5).scale(0.5 + i as f64),
        )
    });
    let links = cycle_links();
    let shift = Vec3::new(2_f64.powi(54), -2_f64.powi(54), 2_f64.powi(54));
    let boost = Vec3::new(0.5, -0.25, 0.125);
    let shifted = core::array::from_fn::<_, 3, _>(|i| {
        state(
            states[i].pose().position_world().add(shift),
            states[i].pose().orientation(),
            states[i]
                .linear_momentum_world()
                .add(boost.scale(properties[i].mass())),
            states[i].angular_momentum_body(),
        )
    });
    let base = PreparedViscousNetwork::new(
        &states.map(|s| s.pose()),
        &properties,
        &links,
        0.1,
        4,
        limits(),
        || false,
    )
    .unwrap()
    .relax(&states, controls(), || false)
    .unwrap();
    let mut changed = controls();
    changed.impulse.angular_origin_world_m = shift;
    let moved = PreparedViscousNetwork::new(
        &shifted.map(|s| s.pose()),
        &properties,
        &links,
        0.1,
        4,
        limits(),
        || false,
    )
    .unwrap()
    .relax(&shifted, changed, || false)
    .unwrap();
    for i in 0..3 {
        vector(
            moved.states_after[i].linear_momentum_world(),
            base.states_after[i]
                .linear_momentum_world()
                .add(boost.scale(properties[i].mass())),
            3e-13,
        );
        vector(
            moved.states_after[i].angular_momentum_body(),
            base.states_after[i].angular_momentum_body(),
            3e-13,
        );
    }
    near(moved.dissipated_energy_j, base.dissipated_energy_j, 3e-13);
    near(
        moved.kinetic_energy_change_j,
        base.kinetic_energy_change_j,
        3e-13,
    );
}

#[test]
fn inactive_edges_empty_network_and_zero_duration_keep_native_identity() {
    let states = line_states([0.4, -0.3, 0.7]);
    let properties = line_properties();
    let links = [
        axial_link([0, 2], 0.),
        axial_link([0, 1], 0.7),
        axial_link([1, 2], 0.),
    ];
    let prepared = PreparedViscousNetwork::new(
        &states.map(|s| s.pose()),
        &properties,
        &links,
        0.3,
        4,
        limits(),
        || false,
    )
    .unwrap();
    let out = prepared.relax(&states, controls(), || false).unwrap();
    assert_eq!(prepared.active_connectors(), 1);
    assert_eq!(out.pair_relaxations, 4);
    assert_eq!(out.states_after[2], states[2]);
    assert_eq!(out.connector_dissipation_j[0], 0.);
    assert_eq!(out.connector_dissipation_j[2], 0.);
    let pair = PreparedViscousConnector::new(
        [states[0].pose(), states[1].pose()],
        [properties[0], properties[1]],
        [Vec3::ZERO; 2],
        links[1].damping.clone(),
        0.3,
        || false,
    )
    .unwrap()
    .relax([states[0], states[1]], controls(), || false)
    .unwrap();
    for i in 0..2 {
        vector(
            out.states_after[i].linear_momentum_world(),
            pair.impulse.bodies[i].state_after.linear_momentum_world(),
            2e-15,
        );
    }
    near(out.dissipated_energy_j, pair.dissipated_energy_j, 2e-15);
    let mut idle = limits();
    idle.maximum_pair_relaxations = 0;
    for edge_set in [&links[..], &[][..]] {
        let p = PreparedViscousNetwork::new(
            &states.map(|s| s.pose()),
            &properties,
            edge_set,
            0.,
            8,
            idle,
            || false,
        )
        .unwrap();
        let r = p.relax(&states, controls(), || false).unwrap();
        assert_eq!(r.states_after, states);
        assert_eq!(r.dissipated_energy_j, 0.);
        assert_eq!(r.pair_relaxations, 0);
    }
    let empty = PreparedViscousNetwork::new(
        &[],
        &[],
        &[],
        0.,
        1,
        ViscousNetworkLimits {
            maximum_bodies: 0,
            maximum_connectors: 0,
            maximum_pair_relaxations: 0,
            maximum_live_bytes: 0,
        },
        || false,
    )
    .unwrap()
    .relax(&[], controls(), || false)
    .unwrap();
    assert!(empty.states_after.is_empty());
    assert!(empty.connector_dissipation_j.is_empty());
    assert_eq!(empty.planned_live_bytes, 0);
}

#[test]
fn budgets_stale_bindings_mid_sweep_cancellation_and_late_physical_refusal_are_atomic() {
    let states = full_states();
    let properties = full_properties();
    let links = cycle_links();
    let poses = states.map(|s| s.pose());
    let prepared =
        PreparedViscousNetwork::new(&poses, &properties, &links, 0.3, 4, limits(), || false)
            .unwrap();
    let reference = prepared.relax(&states, controls(), || false).unwrap();
    for (label, bad) in [
        (
            "bodies",
            ViscousNetworkLimits {
                maximum_bodies: 2,
                ..limits()
            },
        ),
        (
            "connectors",
            ViscousNetworkLimits {
                maximum_connectors: 2,
                ..limits()
            },
        ),
        (
            "pair relaxations",
            ViscousNetworkLimits {
                maximum_pair_relaxations: 19,
                ..limits()
            },
        ),
        (
            "live bytes",
            ViscousNetworkLimits {
                maximum_live_bytes: prepared.planned_live_bytes() - 1,
                ..limits()
            },
        ),
    ] {
        assert!(matches!(
            PreparedViscousNetwork::new(&poses,&properties,&links,0.3,4,bad,||false),
            Err(ViscousNetworkError::Resource(InternalWrenchNetworkError::Limit {resource,..})) if resource==label
        ));
    }
    assert!(
        PreparedViscousNetwork::new(&poses, &properties, &links, 0.3, 0, limits(), || false)
            .is_err()
    );
    // An odd subnormal full interval can round upward on halving: two
    // rounded halves would then evolve noncentral edges for a different time.
    assert!(matches!(
        PreparedViscousNetwork::new(
            &poses,
            &properties,
            &links[..2],
            f64::from_bits(3),
            1,
            limits(),
            || false,
        ),
        Err(ViscousNetworkError::Invalid(
            "positive source subdivision is unrepresentable"
        ))
    ));
    let mut invalid = links.clone();
    invalid[1].bodies = [1, 1];
    assert!(
        PreparedViscousNetwork::new(&poses, &properties, &invalid, 0.3, 4, limits(), || false)
            .is_err()
    );
    let mut changed = states;
    changed[1] = state(
        poses[1].position_world().add(X),
        poses[1].orientation(),
        states[1].linear_momentum_world(),
        states[1].angular_momentum_body(),
    );
    assert!(matches!(
        prepared.relax(&changed, controls(), || false),
        Err(ViscousNetworkError::Invalid(_))
    ));
    let mut polls = 0;
    assert_eq!(
        prepared.relax(&states, controls(), || {
            polls += 1;
            polls == 50
        }),
        Err(ViscousNetworkError::Cancelled)
    );
    assert_eq!(states, full_states());
    assert_eq!(
        prepared.relax(&states, controls(), || false).unwrap(),
        reference
    );
    let mut excessive = states;
    excessive[2] = state(
        poses[2].position_world(),
        poses[2].orientation(),
        X.scale(1e30),
        states[2].angular_momentum_body(),
    );
    assert!(matches!(
        prepared.relax(&excessive, controls(), || false),
        Err(ViscousNetworkError::Connector { connector: 1, .. })
    ));
    assert_eq!(excessive[0], states[0]);
    assert_eq!(excessive[1], states[1]);
    assert_eq!(
        prepared.relax(&states, controls(), || false).unwrap(),
        reference
    );
}
