//! Moving driven graphs reuse the original native body and source owners.
use super::*;
use fs_mbd::spatial_loads::internal_motion::viscous::{
    driven::{HeldInternalEffort, step_driven_viscous_connector_with_loads},
    network::driven::{
        DrivenViscousNetworkLimits, DrivenViscousNetworkMotionControls,
        DrivenViscousNetworkMotionError, step_driven_viscous_network,
        step_driven_viscous_network_with_loads,
    },
};
fn driven_controls() -> DrivenViscousNetworkMotionControls {
    DrivenViscousNetworkMotionControls {
        source: controls().source,
        limits: DrivenViscousNetworkLimits {
            maximum_bodies: 8,
            maximum_connectors: 8,
            maximum_live_bytes: 64_000_000,
            maximum_preparation_terms: 40_000_000_000,
        },
    }
}

#[test]
fn g3_moving_single_edge_matches_the_original_driven_pair_with_external_loads() {
    let (states, properties, links) = full();
    let law = ViscousMotionLaw {
        arms_body_m: links[0].arms_body_m,
        damping: links[0].damping.clone(),
    };
    let efforts = HeldInternalEffort::World([0.2, -0.1, 0.3, 0.04, 0.02, -0.03]);
    let wrenches = [
        Wrench {
            force_world: Vec3::new(0.1, 0.2, -0.1),
            torque_body: Vec3::new(0.01, 0.0, 0.02),
        },
        Wrench {
            force_world: Vec3::new(-0.1, 0.05, 0.2),
            torque_body: Vec3::new(0.0, 0.03, -0.01),
        },
    ];
    let loads = HeldExternalLoads {
        wrenches: &wrenches,
        gravity: Gravity::new(Vec3::new(0.0, 0.0, -0.2)).unwrap(),
    };
    let expected = step_driven_viscous_connector_with_loads(
        [states[0], states[1]],
        [properties[0], properties[1]],
        &law,
        efforts,
        loads,
        0.002,
        controls().source,
        || false,
    )
    .unwrap();
    let out = step_driven_viscous_network_with_loads(
        &states[..2],
        &properties[..2],
        &links[..1],
        &[efforts],
        loads,
        0.002,
        driven_controls(),
        || false,
    )
    .unwrap();
    for i in 0..2 {
        vector(
            out.states_after()[i].pose().position_world(),
            expected.states_after[i].pose().position_world(),
            3e-11,
        );
        vector(
            out.states_after()[i].linear_momentum_world(),
            expected.states_after[i].linear_momentum_world(),
            3e-11,
        );
        vector(
            out.states_after()[i].angular_momentum_body(),
            expected.states_after[i].angular_momentum_body(),
            3e-11,
        );
    }
    near(out.actuator_work_j, expected.actuator_work_j, 3e-11);
    near(out.dissipated_energy_j, expected.dissipated_energy_j, 3e-11);
    near(
        out.external_loads.applied_work_j,
        expected.external_loads.applied_work_j,
        3e-11,
    );
    near(
        out.gravitational_potential_energy_change_j,
        expected.gravitational_potential_energy_change_j,
        3e-11,
    );
}

#[test]
fn g1_driven_cycle_positions_refine_quadratically_without_graph_source_splitting() {
    let (_, properties, _) = line();
    let states = core::array::from_fn::<_, 3, _>(|i| {
        state(
            properties[i],
            X.scale(i as f64 - 1.0),
            UnitQuaternion::IDENTITY,
            X.scale([1.0, 0.0, -1.0][i]),
            Vec3::ZERO,
        )
    });
    let links = [axial([0, 1], 0.7), axial([1, 2], 0.7), axial([0, 2], 0.7)];
    let force = 0.4;
    let efforts = [
        HeldInternalEffort::MaterialFirst([force, 0.0, 0.0, 0.0, 0.0, 0.0]),
        HeldInternalEffort::MaterialFirst([force, 0.0, 0.0, 0.0, 0.0, 0.0]),
        HeldInternalEffort::MaterialFirst([2.0 * force, 0.0, 0.0, 0.0, 0.0, 0.0]),
    ];
    let t = 0.5;
    let rate = 1.05_f64;
    let terminal = force / 0.7;
    let amplitude = -1.0 - terminal;
    let slip = terminal + amplitude * (-rate * t).exp();
    let integral = terminal * t + amplitude * -(-rate * t).exp_m1() / rate;
    let second = terminal * terminal * t
        + 2.0 * terminal * amplitude * -(-rate * t).exp_m1() / rate
        + amplitude * amplitude * -(-2.0 * rate * t).exp_m1() / (2.0 * rate);
    let mut errors = [0.0; 3];
    for (rung, steps) in [2, 4, 8].into_iter().enumerate() {
        let mut current = states.to_vec();
        let mut heat = 0.0;
        let mut work = 0.0;
        for _ in 0..steps {
            let out = step_driven_viscous_network(
                &current,
                &properties,
                &links,
                &efforts,
                t / steps as f64,
                driven_controls(),
                || false,
            )
            .unwrap();
            heat += out.dissipated_energy_j;
            work += out.actuator_work_j;
            near(out.work_heat_residual_j, 0.0, 5e-12);
            assert!(out.connector_dissipation_j.iter().all(|x| *x >= 0.0));
            current = out.into_states_after();
        }
        near(heat, 6.0 * 0.7 * second, 2e-11);
        near(work, 6.0 * force * integral, 2e-11);
        for i in 0..3 {
            near(
                current[i]
                    .center_of_mass_velocity_world(properties[i])
                    .unwrap()
                    .x,
                [-slip, 0.0, slip][i],
                2e-11,
            );
            errors[rung] += (current[i].pose().position_world().x
                - [-1.0 - integral, 0.0, 1.0 + integral][i])
                .abs();
        }
    }
    for pair in errors.windows(2) {
        assert!((3.9..4.1).contains(&(pair[0] / pair[1])), "{errors:?}");
    }
}

#[test]
fn g0_graph_motion_bounds_and_late_cancellation_preserve_all_original_states() {
    let (states, properties, links) = line();
    let original = states;
    let efforts = [HeldInternalEffort::World([0.2, 0.0, 0.0, 0.0, 0.0, 0.0]); 2];
    let mut count = 0;
    let expected = step_driven_viscous_network(
        &states,
        &properties,
        &links,
        &efforts,
        0.02,
        driven_controls(),
        || {
            count += 1;
            false
        },
    )
    .unwrap();
    let mut seen = 0;
    assert_eq!(
        step_driven_viscous_network(
            &states,
            &properties,
            &links,
            &efforts,
            0.02,
            driven_controls(),
            || {
                seen += 1;
                seen == count
            },
        ),
        Err(DrivenViscousNetworkMotionError::Cancelled)
    );
    assert_eq!(states, original);
    let retry = step_driven_viscous_network(
        &states,
        &properties,
        &links,
        &efforts,
        0.02,
        driven_controls(),
        || false,
    )
    .unwrap();
    assert_eq!(retry, expected);
    for limit in [
        DrivenViscousNetworkLimits {
            maximum_live_bytes: expected.planned_live_bytes - 1,
            ..driven_controls().limits
        },
        DrivenViscousNetworkLimits {
            maximum_preparation_terms: expected.planned_preparation_terms - 1,
            ..driven_controls().limits
        },
    ] {
        assert!(
            step_driven_viscous_network(
                &states,
                &properties,
                &links,
                &efforts,
                0.02,
                DrivenViscousNetworkMotionControls {
                    limits: limit,
                    ..driven_controls()
                },
                || false
            )
            .is_err()
        );
    }
}
