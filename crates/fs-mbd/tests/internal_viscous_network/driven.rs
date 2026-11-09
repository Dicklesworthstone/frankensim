//! Simultaneous graph actuation, original constitutive loss and native recoil.
use super::*;
use fs_mbd::spatial_loads::internal_impulse::viscous::{
    driven::PreparedDrivenViscousConnector,
    network::driven::{
        DrivenViscousNetworkError, DrivenViscousNetworkLimits, PreparedDrivenViscousNetwork,
    },
};

fn driven_limits() -> DrivenViscousNetworkLimits {
    DrivenViscousNetworkLimits {
        maximum_bodies: 8,
        maximum_connectors: 8,
        maximum_live_bytes: 64_000_000,
        maximum_preparation_terms: 20_000_000_000,
    }
}

fn integral_decay(rate: f64, dt: f64) -> f64 {
    -(-rate * dt).exp_m1() / rate
}

#[test]
fn g1_simultaneous_cycle_retains_terminal_slip_and_separate_link_heat() {
    let states = line_states([1.0, 0.0, -1.0]);
    let properties = line_properties();
    let c = 0.7;
    let links = [
        axial_link([0, 1], c),
        axial_link([1, 2], c),
        axial_link([0, 2], c),
    ];
    let dt = 0.35;
    for force in [0.4, -0.2, 0.0] {
        let efforts = [
            [force, 0.0, 0.0, 0.0, 0.0, 0.0],
            [force, 0.0, 0.0, 0.0, 0.0, 0.0],
            [2.0 * force, 0.0, 0.0, 0.0, 0.0, 0.0],
        ];
        let mut source = PreparedDrivenViscousNetwork::new(
            &states.map(|s| s.pose()),
            &properties,
            &links,
            dt,
            driven_limits(),
            || false,
        )
        .unwrap();
        let out = source
            .advance(&states, &efforts, controls(), || false)
            .unwrap();
        // Complete K3 graph mode v=(-u,0,u), with u'=(3/m)*(f-c*u).
        let rate = 1.5 * c;
        let terminal = force / c;
        let amplitude = -1.0 - terminal;
        let u = terminal + amplitude * (-rate * dt).exp();
        let first = terminal * dt + amplitude * integral_decay(rate, dt);
        let second = terminal * terminal * dt
            + 2.0 * terminal * amplitude * integral_decay(rate, dt)
            + amplitude * amplitude * integral_decay(2.0 * rate, dt);
        for i in 0..3 {
            near(
                out.states_after[i]
                    .center_of_mass_velocity_world(properties[i])
                    .unwrap()
                    .x,
                [-u, 0.0, u][i],
                3e-12,
            );
            near(
                out.connector_actuator_work_j[i],
                [1.0, 1.0, 4.0][i] * force * first,
                3e-12,
            );
            near(
                out.connector_dissipation_j[i],
                [1.0, 1.0, 4.0][i] * c * second,
                3e-12,
            );
        }
        near(out.kinetic_balance_residual_j, 0.0, 4e-12);
        near(out.connector_balance_residual_j, 0.0, 4e-12);
        near(out.heat_partition_residual_j, 0.0, 4e-12);
        near(out.actuator_partition_residual_j, 0.0, 4e-12);
        vector(out.linear_momentum_change_world_ns, Vec3::ZERO, 3e-13);
        vector(out.angular_momentum_change_world_nms, Vec3::ZERO, 3e-13);
        if force > 0.0 {
            assert!(out.actuator_work_j < 0.0);
        }
    }
}

#[test]
fn g1_parallel_drives_can_supply_and_recover_work_on_one_trajectory() {
    let states = line_states([0.0, 1.0, 0.0]);
    let properties = line_properties();
    let links = [axial_link([0, 1], 0.0), axial_link([0, 1], 0.0)];
    let mut source = PreparedDrivenViscousNetwork::new(
        &states.map(|s| s.pose()),
        &properties,
        &links,
        0.4,
        driven_limits(),
        || false,
    )
    .unwrap();
    let out = source
        .advance(
            &states,
            &[
                [2.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                [-1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            ],
            controls(),
            || false,
        )
        .unwrap();
    // Reduced mass is one kilogram: s(t)=1+t, integral s dt=.48.
    near(out.connector_actuator_work_j[0], 0.96, 2e-13);
    near(out.connector_actuator_work_j[1], -0.48, 2e-13);
    near(out.kinetic_energy_change_j, 0.48, 2e-13);
    assert_eq!(out.connector_dissipation_j, [0.0, 0.0]);
    assert_eq!(out.states_after[2], states[2]);
    near(out.states_after[0].linear_momentum_world().x, -0.4, 2e-13);
    near(out.states_after[1].linear_momentum_world().x, 2.4, 2e-13);
}

#[test]
fn g3_single_edge_matches_original_full_tensor_driven_pair() {
    let states = full_states();
    let properties = full_properties();
    let link = cycle_links()[0].clone();
    let efforts = [0.2, -0.1, 0.3, 0.04, -0.03, 0.02];
    let mut pair = PreparedDrivenViscousConnector::new(
        [states[0].pose(), states[1].pose()],
        [properties[0], properties[1]],
        link.arms_body_m,
        link.damping.clone(),
        0.2,
        || false,
    )
    .unwrap();
    let expected = pair
        .advance([states[0], states[1]], efforts, controls(), || false)
        .unwrap();
    let mut graph = PreparedDrivenViscousNetwork::new(
        &states.map(|s| s.pose()),
        &properties,
        &[link],
        0.2,
        driven_limits(),
        || false,
    )
    .unwrap();
    let out = graph
        .advance(&states, &[efforts], controls(), || false)
        .unwrap();
    for i in 0..2 {
        vector(
            out.states_after[i].linear_momentum_world(),
            expected.impulse.bodies[i]
                .state_after
                .linear_momentum_world(),
            4e-11,
        );
        vector(
            out.states_after[i].angular_momentum_body(),
            expected.impulse.bodies[i]
                .state_after
                .angular_momentum_body(),
            4e-11,
        );
    }
    near(out.actuator_work_j, expected.actuator_work_j, 4e-11);
    near(out.dissipated_energy_j, expected.dissipated_energy_j, 4e-11);
    assert_eq!(out.states_after[2], states[2]);
}

// Independent original native-coordinate ODE: actual masses, full inverse
// inertias, physical arms and transported force/couple. No production graph
// matrix, body kinetic transform or source observer appears in this reference.
fn driven_rhs(
    y: [f64; 24],
    poses: &[Pose; 3],
    links: &[ViscousNetworkLink; 3],
    efforts: &[[f64; 6]; 3],
) -> [f64; 24] {
    let velocity: [Vec3; 3] = core::array::from_fn(|i| {
        v([y[6 * i], y[6 * i + 1], y[6 * i + 2]]).scale(1.0 / [2.0, 3.0, 5.0][i])
    });
    let spin: [Vec3; 3] = core::array::from_fn(|i| {
        poses[i]
            .orientation()
            .rotate_body_to_world(omega(v([y[6 * i + 3], y[6 * i + 4], y[6 * i + 5]]), i))
    });
    let mut out = [0.0; 24];
    for (e, link) in links.iter().enumerate() {
        let [first, second] = link.bodies;
        let qa = poses[first].orientation();
        let rb = poses[second]
            .orientation()
            .rotate_body_to_world(link.arms_body_m[1]);
        let arm = poses[second]
            .position_world()
            .sub(poses[first].position_world())
            .add(rb);
        let linear = qa.rotate_world_to_body(
            velocity[second]
                .sub(velocity[first])
                .add(spin[second].cross(rb))
                .sub(spin[first].cross(arm)),
        );
        let angular = qa.rotate_world_to_body(spin[second].sub(spin[first]));
        let s = [
            linear.x, linear.y, linear.z, angular.x, angular.y, angular.z,
        ];
        let c = link.damping.native_coefficients();
        let viscous: [f64; 6] = core::array::from_fn(|i| (0..6).map(|j| c[i][j] * s[j]).sum());
        let f = qa.rotate_body_to_world(v(core::array::from_fn(|i| efforts[e][i] - viscous[i])));
        let couple = qa.rotate_body_to_world(v(core::array::from_fn(|i| {
            efforts[e][i + 3] - viscous[i + 3]
        })));
        let moments = [
            arm.cross(f).add(couple).scale(-1.0),
            rb.cross(f).add(couple),
        ];
        for (side, body) in [first, second].into_iter().enumerate() {
            let moment = poses[body]
                .orientation()
                .rotate_world_to_body(moments[side]);
            for i in 0..3 {
                out[6 * body + i] += a(f)[i] * if side == 0 { -1.0 } else { 1.0 };
                out[6 * body + 3 + i] += a(moment)[i];
            }
        }
        out[18 + e] = (0..6).map(|i| viscous[i] * s[i]).sum();
        out[21 + e] = (0..6).map(|i| efforts[e][i] * s[i]).sum();
    }
    out
}

#[test]
fn g1_full_tensor_cycle_matches_independent_ode_and_link_reordering() {
    let states = full_states();
    let properties = full_properties();
    let poses = states.map(|s| s.pose());
    let links = cycle_links();
    let efforts = [
        [0.2, -0.1, 0.3, 0.04, -0.03, 0.02],
        [-0.3, 0.2, 0.1, -0.02, 0.05, 0.01],
        [0.1, 0.15, -0.2, 0.02, 0.01, -0.04],
    ];
    let mut y = [0.0; 24];
    y[..18].copy_from_slice(&initial_ode(&states)[..18]);
    let dt = 0.2 / 1024.0;
    for _ in 0..1024 {
        let k1 = driven_rhs(y, &poses, &links, &efforts);
        let k2 = driven_rhs(
            core::array::from_fn(|i| y[i] + 0.5 * dt * k1[i]),
            &poses,
            &links,
            &efforts,
        );
        let k3 = driven_rhs(
            core::array::from_fn(|i| y[i] + 0.5 * dt * k2[i]),
            &poses,
            &links,
            &efforts,
        );
        let k4 = driven_rhs(
            core::array::from_fn(|i| y[i] + dt * k3[i]),
            &poses,
            &links,
            &efforts,
        );
        for i in 0..24 {
            y[i] += dt * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0;
        }
    }
    let mut source = PreparedDrivenViscousNetwork::new(
        &poses,
        &properties,
        &links,
        0.2,
        driven_limits(),
        || false,
    )
    .unwrap();
    let out = source
        .advance(&states, &efforts, controls(), || false)
        .unwrap();
    for i in 0..3 {
        vector(
            out.states_after[i].linear_momentum_world(),
            v(core::array::from_fn(|j| y[6 * i + j])),
            3e-10,
        );
        vector(
            out.states_after[i].angular_momentum_body(),
            v(core::array::from_fn(|j| y[6 * i + 3 + j])),
            3e-10,
        );
        near(out.connector_dissipation_j[i], y[18 + i], 3e-10);
        near(out.connector_actuator_work_j[i], y[21 + i], 3e-10);
        assert!(out.connector_dissipation_j[i] > 0.0);
    }
    let order = [2, 0, 1];
    let mut reordered = PreparedDrivenViscousNetwork::new(
        &poses,
        &properties,
        &order.map(|i| links[i].clone()),
        0.2,
        driven_limits(),
        || false,
    )
    .unwrap();
    let other = reordered
        .advance(&states, &order.map(|i| efforts[i]), controls(), || false)
        .unwrap();
    for i in 0..3 {
        vector(
            other.states_after[i].linear_momentum_world(),
            out.states_after[i].linear_momentum_world(),
            3e-11,
        );
        vector(
            other.states_after[i].angular_momentum_body(),
            out.states_after[i].angular_momentum_body(),
            3e-11,
        );
        near(
            other.connector_dissipation_j[i],
            out.connector_dissipation_j[order[i]],
            3e-11,
        );
        near(
            other.connector_actuator_work_j[i],
            out.connector_actuator_work_j[order[i]],
            3e-11,
        );
    }
}

#[test]
fn g1_weak_parallel_loss_survives_beside_a_stronger_connector() {
    let states = line_states([0.5, -0.5, 0.25]);
    let properties = line_properties();
    let links = [axial_link([0, 1], 1.0), axial_link([0, 1], 1e-24)];
    let mut source = PreparedDrivenViscousNetwork::new(
        &states.map(|s| s.pose()),
        &properties,
        &links,
        0.3,
        driven_limits(),
        || false,
    )
    .unwrap();
    let out = source
        .advance(&states, &[[0.0; 6]; 2], controls(), || false)
        .unwrap();
    let expected = 0.5 * -(-0.6_f64).exp_m1();
    near(out.connector_dissipation_j[0], expected, 3e-13);
    assert!(out.connector_dissipation_j[1] > 0.0);
    near(out.connector_dissipation_j[1] / 1e-24, expected, 3e-13);
    assert_eq!(out.states_after[2], states[2]);
}

#[test]
fn g0_source_refusals_are_atomic_and_work_and_storage_caps_cover_all_links() {
    let states = line_states([1.0, 0.0, -1.0]);
    let properties = line_properties();
    let links = [axial_link([0, 1], 0.7), axial_link([1, 2], 0.7)];
    let poses = states.map(|s| s.pose());
    let mut source = PreparedDrivenViscousNetwork::new(
        &poses,
        &properties,
        &links,
        0.1,
        driven_limits(),
        || false,
    )
    .unwrap();
    let force = [[0.2, 0.0, 0.0, 0.0, 0.0, 0.0]; 2];
    let mut polls = 0;
    let expected = source
        .advance(&states, &force, controls(), || {
            polls += 1;
            false
        })
        .unwrap();
    let mut seen = 0;
    assert_eq!(
        source.advance(&states, &force, controls(), || {
            seen += 1;
            seen == polls
        }),
        Err(DrivenViscousNetworkError::Cancelled)
    );
    assert_eq!(
        source
            .advance(&states, &force, controls(), || false)
            .unwrap(),
        expected
    );
    let bad = [
        DrivenViscousNetworkLimits {
            maximum_bodies: 2,
            ..driven_limits()
        },
        DrivenViscousNetworkLimits {
            maximum_connectors: 1,
            ..driven_limits()
        },
        DrivenViscousNetworkLimits {
            maximum_live_bytes: source.planned_live_bytes() - 1,
            ..driven_limits()
        },
        DrivenViscousNetworkLimits {
            maximum_preparation_terms: source.planned_preparation_terms() - 1,
            ..driven_limits()
        },
    ];
    for limits in bad {
        assert!(
            PreparedDrivenViscousNetwork::new(&poses, &properties, &links, 0.1, limits, || false)
                .is_err()
        );
    }
    let mut invalid = links.clone();
    invalid[1].bodies = [1, 1];
    assert!(
        PreparedDrivenViscousNetwork::new(
            &poses,
            &properties,
            &invalid,
            0.1,
            driven_limits(),
            || false
        )
        .is_err()
    );
    assert!(
        source
            .advance(&states, &force[..1], controls(), || false)
            .is_err()
    );
    let mut invalid_force = force;
    invalid_force[1][2] = f64::NAN;
    assert!(
        source
            .advance(&states, &invalid_force, controls(), || false)
            .is_err()
    );
    let mut changed = states;
    changed[0] = state(
        Vec3::new(0.2, 0.0, 0.0),
        UnitQuaternion::IDENTITY,
        changed[0].linear_momentum_world(),
        Vec3::ZERO,
    );
    assert!(
        source
            .advance(&changed, &force, controls(), || false)
            .is_err()
    );
    let mut zero = PreparedDrivenViscousNetwork::new(
        &poses,
        &properties,
        &links,
        0.0,
        DrivenViscousNetworkLimits {
            maximum_preparation_terms: 0,
            ..driven_limits()
        },
        || false,
    )
    .unwrap();
    let out = zero.advance(&states, &force, controls(), || false).unwrap();
    assert_eq!(out.states_after, states);
    assert_eq!(out.connector_dissipation_j, [0.0, 0.0]);
    assert_eq!(out.connector_actuator_work_j, [0.0, 0.0]);
}
