//! G1/G3/G4: continuous dry friction through native bodies, joints and supports.

use super::*;
use fs_contact::tangential::impact::multiple::CoupledFrictionImpulseControls;
use fs_contact::tangential::impact::multiple::relaxation::CoupledFrictionFlowControls;
use fs_contact::tangential::impact::multiple::relaxation::stiction::{
    CoupledStictionFlowControls, MAX_STICTION_MODES_PER_PHASE,
};
use fs_mbd::articulated_contacts::tangential::flow::{
    ArticulatedDryFrictionControls, ArticulatedDryFrictionResult, ArticulatedDryFrictionSite,
    advance_articulated_dry_friction,
};

fn controls() -> ArticulatedDryFrictionControls {
    ArticulatedDryFrictionControls {
        flow: CoupledStictionFlowControls {
            flow: CoupledFrictionFlowControls {
                projection: CoupledFrictionImpulseControls {
                    maximum_iterations: 128,
                    impulse_tolerance_ns: 2e-11,
                    velocity_tolerance_m_s: 2e-11,
                    energy_tolerance_j: 2e-10,
                    relative_matrix_tolerance: 1e-12,
                    maximum_impulse_ns: 100.0,
                },
                maximum_segments: 64,
                acceleration_tolerance_m_s2: 2e-11,
                maximum_force_n: 100.0,
            },
            maximum_mode_evaluations: MAX_STICTION_MODES_PER_PHASE * 64,
        },
        maximum_duration_s: 10.0,
        maximum_slip_speed_m_s: 20.0,
        geometry_tolerance_m: 1e-10,
        velocity_tolerance_m_s: 2e-11,
        energy_tolerance_j: 2e-10,
    }
}

fn site(law: &FrictionLaw, arm: Vec3) -> ArticulatedDryFrictionSite<'_> {
    ArticulatedDryFrictionSite {
        id: 10,
        contact: contact(arm, 0.0),
        normal_force_n: 1.0,
        law,
        drive_force_n: 0.0,
    }
}

fn audit(out: &ArticulatedDryFrictionResult) {
    close(out.energy_balance_residual_j, 0.0);
    close(out.endpoint_velocity_residual_m_s, 0.0);
    close(out.path.flow.energy_residual_j, 0.0);
    close(
        out.kinetic_energy_change_j + out.path.flow.dissipated_j,
        out.support_work_j + out.path.flow.total_drive_work_j,
    );
    close_vector(out.linear_momentum_residual, Vec3::default());
    close_vector(out.angular_momentum_residual, Vec3::default());
    assert!(out.path.flow.dissipated_j >= 0.0);
    for i in 0..out.path.flow.count {
        assert!(out.path.flow.friction_work_j[i] <= 0.0);
        close(
            out.relative_speeds_after_m_s[i],
            out.path.flow.flow_velocity_after_m_s[i],
        );
    }
}

#[test]
fn off_centre_sliding_arrests_without_midpoint_velocity_reflection() {
    let model = body(2.0, 1.0);
    let initial = system(
        &model,
        Vec3::default(),
        Vec3::new(3.0, 0.0, 0.0),
        Vec3::default(),
    );
    let law = FrictionLaw::Coulomb {
        static_mu: 4.0,
        kinetic_mu: 2.0,
    };
    let port = site(&law, Vec3::new(0.0, -1.0, 0.0));
    for duration in [0.25_f64, 1.0, 2.0] {
        let out =
            advance_articulated_dry_friction(&[initial], &[port], duration, controls(), || false)
                .unwrap();
        // Actual point mobility includes rotation: W=1/m+r²/I=3/2.
        // dv/dt=-3 until t=1, then static force is zero for the rest of the phase.
        let sliding_time = duration.min(1.0);
        let impulse = -2.0 * sliding_time;
        let travel = 3.0 * sliding_time - 1.5 * sliding_time * sliding_time;
        close(out.mobility_per_kg[0], 1.5);
        close(out.path.flow.friction_impulse_ns[0], impulse);
        close(out.path.flow.displacement_m[0], travel);
        close(out.path.flow.dissipated_j, 2.0 * travel);
        close(out.relative_speeds_after_m_s[0], 3.0 - 3.0 * sliding_time);
        close(
            out.responses[0].base_after.twist_body.linear.x,
            3.0 + impulse / 2.0,
        );
        close(out.responses[0].base_after.twist_body.angular.z, impulse);
        assert_eq!(
            out.responses[0].base_after.world_from_base,
            initial.base.world_from_base
        );
        if duration >= 1.0 {
            assert_eq!(out.path.flow.flow_velocity_after_m_s[0], 0.0);
            assert_eq!(out.path.flow.stops, 1);
        }
        audit(&out);
    }
}

#[test]
fn native_driven_stiction_distinguishes_breakaway_from_arbitrarily_small_slip() {
    let model = body(1.0, 1.0);
    let law = FrictionLaw::Coulomb {
        static_mu: 2.0,
        kinetic_mu: 1.0,
    };
    let mut port = site(&law, Vec3::default());
    port.drive_force_n = 1.5;
    let held = system(&model, Vec3::default(), Vec3::default(), Vec3::default());
    let rest =
        advance_articulated_dry_friction(&[held], &[port], 2.0, controls(), || false).unwrap();
    assert_eq!(rest.relative_speeds_after_m_s[0], 0.0);
    close(rest.path.flow.friction_impulse_ns[0], -3.0);
    assert_eq!(rest.path.flow.dissipated_j, 0.0);
    assert_eq!(rest.path.flow.total_drive_work_j, 0.0);
    audit(&rest);

    let v = 1e-15;
    let moving = system(
        &model,
        Vec3::default(),
        Vec3::new(v, 0.0, 0.0),
        Vec3::default(),
    );
    let out =
        advance_articulated_dry_friction(&[moving], &[port], 2.0, controls(), || false).unwrap();
    close(out.relative_speeds_after_m_s[0], 1.0 + v);
    close(out.path.flow.displacement_m[0], 1.0 + 2.0 * v);
    close(out.path.flow.friction_impulse_ns[0], -2.0);
    close(out.path.flow.total_drive_work_j, 1.5 + 3.0 * v);
    assert_eq!(out.path.flow.minimum_speed_m_s[0], v);
    assert_eq!(out.path.flow.stops, 0);
    audit(&out);

    port.drive_force_n = 2.5;
    let released =
        advance_articulated_dry_friction(&[held], &[port], 1.0, controls(), || false).unwrap();
    close(released.relative_speeds_after_m_s[0], 1.5);
    close(released.path.flow.distance_m[0], 0.75);
    close(released.path.flow.dissipated_j, 0.75);
    assert_eq!(released.path.diagnostics.breakaways, 1);
    audit(&released);
}

#[test]
fn simultaneous_ports_retain_static_recoil_and_continuous_work_under_permutation() {
    // m=2/3, I=2 and arms y=+/-1 give W=[[2,1],[1,2]].
    let model = body(2.0 / 3.0, 2.0);
    let initial = system(
        &model,
        Vec3::default(),
        Vec3::new(1.0, 0.0, 0.0),
        Vec3::new(0.0, 0.0, 1.0),
    );
    let laws = [
        FrictionLaw::Coulomb {
            static_mu: 2.0,
            kinetic_mu: 0.5,
        },
        FrictionLaw::Coulomb {
            static_mu: 2.0,
            kinetic_mu: 1.0,
        },
    ];
    let ports = [
        site(&laws[0], Vec3::new(0.0, 1.0, 0.0)),
        ArticulatedDryFrictionSite {
            id: 20,
            ..site(&laws[1], Vec3::new(0.0, -1.0, 0.0))
        },
    ];
    let mut baseline = None;
    for ordered in [ports, [ports[1], ports[0]]] {
        let out = advance_articulated_dry_friction(&[initial], &ordered, 1.0, controls(), || false)
            .unwrap();
        assert_eq!(out.contact_ids, [10, 20]);
        for (actual, expected) in out.mobility_per_kg.iter().zip([2.0, 1.0, 1.0, 2.0]) {
            close(*actual, expected);
        }
        // Holding the first point requires f1=+1/2 while f2=-1 slides.
        close(out.path.flow.friction_impulse_ns[0], 0.5);
        close(out.path.flow.friction_impulse_ns[1], -1.0);
        close(out.relative_speeds_after_m_s[0], 0.0);
        close(out.relative_speeds_after_m_s[1], 0.5);
        close(out.path.flow.distance_m[0], 0.0);
        close(out.path.flow.distance_m[1], 1.25);
        close(out.path.flow.dissipated_j, 1.25);
        close(out.responses[0].base_after.twist_body.linear.x, 0.25);
        close(out.responses[0].base_after.twist_body.angular.z, 0.25);
        if let Some(previous) = &baseline {
            assert_eq!(&out.path, previous);
        } else {
            baseline = Some(out.path.clone());
        }
        audit(&out);
    }
    let stopped =
        advance_articulated_dry_friction(&[initial], &ports, 2.0, controls(), || false).unwrap();
    close(stopped.path.flow.dissipated_j, 4.0 / 3.0);
    close(stopped.path.flow.friction_impulse_ns[0], 2.0 / 3.0);
    close(stopped.path.flow.friction_impulse_ns[1], -4.0 / 3.0);
    assert_eq!(stopped.path.flow.flow_velocity_after_m_s[..2], [0.0, 0.0]);
    audit(&stopped);
}

#[test]
fn finite_counterpart_recoils_and_moving_support_work_has_its_own_source() {
    let models = [body(2.0, 1.0), body(3.0, 1.0)];
    let law = FrictionLaw::Coulomb {
        static_mu: 2.0,
        kinetic_mu: 1.0,
    };
    for boost in [0.0, 7.0] {
        let systems = [
            system(
                &models[0],
                Vec3::default(),
                Vec3::new(3.0 + boost, 0.0, 0.0),
                Vec3::default(),
            ),
            system(
                &models[1],
                Vec3::default(),
                Vec3::new(-2.0 + boost, 0.0, 0.0),
                Vec3::default(),
            ),
        ];
        let mut port = site(&law, Vec3::default());
        port.contact.other = ContactCounterpart::Body(LinkContactPoint {
            system_index: 1,
            link_index: 0,
            point_body: Vec3::default(),
        });
        let out = advance_articulated_dry_friction(&systems, &[port], 10.0, controls(), || false)
            .unwrap();
        close(out.mobility_per_kg[0], 5.0 / 6.0);
        close(out.path.flow.friction_impulse_ns[0], -6.0);
        close(out.path.flow.distance_m[0], 15.0);
        close(out.path.flow.dissipated_j, 15.0);
        for response in &out.responses {
            close(response.base_after.twist_body.linear.x, boost);
        }
        close(out.support_work_j, 0.0);
        audit(&out);
    }

    let model = body(1.0, 1.0);
    let initial = system(&model, Vec3::default(), Vec3::default(), Vec3::default());
    let mut port = site(&law, Vec3::default());
    port.contact = contact(Vec3::default(), 2.0);
    let out =
        advance_articulated_dry_friction(&[initial], &[port], 3.0, controls(), || false).unwrap();
    close(out.responses[0].base_after.twist_body.linear.x, 2.0);
    close(out.path.flow.displacement_m[0], -2.0);
    close(out.path.flow.dissipated_j, 2.0);
    close(out.support_work_j, 4.0);
    close(out.kinetic_energy_change_j, 2.0);
    audit(&out);
}

#[test]
fn friction_arrest_uses_original_prismatic_and_revolute_joint_response() {
    let law = FrictionLaw::Coulomb {
        static_mu: 2.0,
        kinetic_mu: 1.0,
    };
    let model = super::super::friction::two_link(
        JointModel::prismatic(Vec3::new(1.0, 0.0, 0.0), None).unwrap(),
    );
    let initial = ArticulatedContactSystem {
        model: &model,
        base: FreeFloatingBaseState::new(
            Se3::identity(),
            Twist::new(Vec3::default(), Vec3::new(0.3, 0.0, 0.0)),
        ),
        joint_position: &[0.0],
        joint_velocity: &[0.7],
    };
    let mut port = site(&law, Vec3::default());
    port.contact.point.link_index = 1;
    let out =
        advance_articulated_dry_friction(&[initial], &[port], 2.0, controls(), || false).unwrap();
    close(out.mobility_per_kg[0], 1.0);
    close(out.responses[0].base_after.twist_body.linear.x, 0.3);
    close(out.responses[0].joint_velocity_after[0], -0.3);
    close(out.path.flow.dissipated_j, 0.5);
    audit(&out);

    let model = super::super::friction::two_link(
        JointModel::revolute(Vec3::new(0.0, 0.0, 1.0), None).unwrap(),
    );
    let initial = ArticulatedContactSystem {
        model: &model,
        base: FreeFloatingBaseState::new(
            Se3::identity(),
            Twist::new(Vec3::default(), Vec3::new(0.0, 0.1, 0.0)),
        ),
        joint_position: &[0.0],
        joint_velocity: &[2.0],
    };
    let mut port = site(&law, Vec3::new(1.0, 0.0, 0.0));
    port.contact.point.link_index = 1;
    port.contact.tangent_world = Vec3::new(0.0, 1.0, 0.0);
    let out =
        advance_articulated_dry_friction(&[initial], &[port], 2.0, controls(), || false).unwrap();
    close(out.mobility_per_kg[0], 7.0 / 3.0);
    close(out.responses[0].base_after.twist_body.linear.y, -0.2);
    close(out.responses[0].joint_velocity_after[0], 0.2);
    close(out.path.flow.friction_impulse_ns[0], -0.9);
    close(out.path.flow.dissipated_j, 0.945);
    audit(&out);
}

#[test]
fn refusals_and_cancellation_preserve_native_inputs_and_retry() {
    let model = body(2.0, 1.0);
    let initial = system(
        &model,
        Vec3::default(),
        Vec3::new(3.0, 0.0, 0.0),
        Vec3::default(),
    );
    let law = FrictionLaw::Coulomb {
        static_mu: 4.0,
        kinetic_mu: 2.0,
    };
    let port = site(&law, Vec3::new(0.0, -1.0, 0.0));
    let mut polls = 0;
    let baseline = advance_articulated_dry_friction(&[initial], &[port], 2.0, controls(), || {
        polls += 1;
        false
    })
    .unwrap();
    for stop_at in [1, polls / 2, polls] {
        let mut count = 0;
        assert!(
            advance_articulated_dry_friction(&[initial], &[port], 2.0, controls(), || {
                count += 1;
                count == stop_at
            })
            .is_err()
        );
        assert_eq!(initial.base.twist_body.linear.x, 3.0);
    }
    let retry =
        advance_articulated_dry_friction(&[initial], &[port], 2.0, controls(), || false).unwrap();
    assert_eq!(retry.path, baseline.path);
    assert_eq!(
        retry.responses[0].base_after,
        baseline.responses[0].base_after
    );
    assert!(
        advance_articulated_dry_friction(&[initial], &[port, port], 2.0, controls(), || false)
            .is_err()
    );
    let mut config = controls();
    config.flow.flow.maximum_segments = 1;
    assert!(advance_articulated_dry_friction(&[initial], &[port], 2.0, config, || false).is_err());
    let mut other = port;
    other.id += 1;
    assert!(
        advance_articulated_dry_friction(&[initial], &[port, other], 2.0, controls(), || false)
            .is_err(),
        "dependent mobility must not choose arbitrary reactions"
    );
    let rate = FrictionLaw::Stribeck {
        static_mu: 2.0,
        kinetic_mu: 1.0,
        characteristic_speed: 0.1,
        viscous_per_speed: 0.0,
    };
    assert!(
        advance_articulated_dry_friction(
            &[initial],
            &[ArticulatedDryFrictionSite { law: &rate, ..port }],
            2.0,
            controls(),
            || false
        )
        .is_err()
    );
    let identity =
        advance_articulated_dry_friction(&[initial], &[port], 0.0, controls(), || false).unwrap();
    assert_eq!(identity.responses[0].base_after, initial.base);
    assert_eq!(identity.path.flow.dissipated_j, 0.0);
    let empty =
        advance_articulated_dry_friction(&[initial], &[], 2.0, controls(), || false).unwrap();
    assert_eq!(empty.responses[0].base_after, initial.base);
    assert!(empty.contact_ids.is_empty());
    audit(&empty);
}
