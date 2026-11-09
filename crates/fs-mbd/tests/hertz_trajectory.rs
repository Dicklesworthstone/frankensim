//! Actual Hertz stopping times, finite recoil, travel and physical heat.

use core::f64::consts::PI;
use fs_exec::CancelGate;
use fs_mbd::tangential_friction::RigidContactBinding;
use fs_mbd::tangential_friction::hertz::trajectory::*;
use fs_mbd::tangential_friction::hertz::{RigidHertzContact, RigidHertzControls};
use fs_mbd::{MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3};
use fs_tribo::partial_slip::TangentFrame;
use fs_tribo::sliding_patch::contact::relaxation::HertzRelaxationError;
use fs_tribo::sliding_patch::contact::{
    HertzWrenchConfig, HertzWrenchContact, HertzWrenchParameters,
};

const RADIUS: f64 = 0.03;

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 3e-12 * expected.abs().max(1e-2),
        "actual {actual:.17e}, expected {expected:.17e}"
    );
}
fn close_vec(actual: Vec3, expected: Vec3) {
    close(actual.x, expected.x);
    close(actual.y, expected.y);
    close(actual.z, expected.z);
}
fn controls() -> RigidHertzControls {
    RigidHertzControls {
        maximum_bodies: 4,
        geometry_tolerance_m: 2e-12,
        velocity_tolerance_m_s: 2e-11,
        angular_velocity_tolerance_rad_s: 2e-10,
        linear_momentum_tolerance_ns: 2e-11,
        angular_momentum_tolerance_nms: 2e-11,
        energy_tolerance_j: 2e-11,
        contact: HertzWrenchConfig {
            maximum_static_iterations: 1,
            maximum_iterations: 1,
            maximum_evaluations: 1,
            force_tolerance_n: 2e-12,
            moment_tolerance_nm: 2e-12,
            displacement_tolerance_m: 2e-12,
            energy_tolerance_j: 2e-12,
            relative_matrix_tolerance: 1e-13,
            maximum_force_n: 20.0,
            maximum_spin_moment_nm: 2.0,
        },
    }
}
fn law() -> HertzWrenchContact {
    HertzWrenchContact::new(HertzWrenchParameters {
        normal_force_n: 10.0,
        contact_radius_m: RADIUS,
        static_mu: 0.8,
        kinetic_mu: 0.4,
    })
    .unwrap()
}
fn state(pose: Pose, mass: MassProperties, velocity: Vec3, spin: f64) -> RigidBodyState {
    RigidBodyState::new(
        pose,
        velocity.scale(mass.mass()),
        Vec3::new(0.0, 0.0, spin * mass.principal_inertia_body().z),
    )
    .unwrap()
}
struct Fixture {
    bodies: [RigidBodyState; 2],
    properties: [MassProperties; 2],
    binding: RigidContactBinding,
    frame: TangentFrame,
}
impl Fixture {
    fn contact<'a>(&self, law: &'a HertzWrenchContact) -> RigidHertzContact<'a> {
        RigidHertzContact {
            binding: self.binding,
            frame: self.frame,
            law,
        }
    }
}
fn axial(slip: f64, spin: f64) -> Fixture {
    let properties = [
        MassProperties::new(2.0, Vec3::ZERO, Vec3::new(0.08, 0.10, 0.12)).unwrap(),
        MassProperties::new(3.0, Vec3::ZERO, Vec3::new(0.15, 0.20, 0.25)).unwrap(),
    ];
    let poses = [
        Pose::new(Vec3::new(0.0, 0.0, 0.2), UnitQuaternion::IDENTITY).unwrap(),
        Pose::new(Vec3::new(0.0, 0.0, -0.15), UnitQuaternion::IDENTITY).unwrap(),
    ];
    Fixture {
        bodies: [
            state(poses[0], properties[0], Vec3::new(slip, 0.0, 0.0), spin),
            state(poses[1], properties[1], Vec3::ZERO, 0.0),
        ],
        properties,
        binding: RigidContactBinding {
            body_a: 0,
            body_b: 1,
            arm_a_body: Vec3::new(0.0, 0.0, -0.2),
            arm_b_body: Vec3::new(0.0, 0.0, 0.15),
        },
        frame: TangentFrame::new([0.0, 0.0, 1.0], [1.0, 0.0, 0.0]).unwrap(),
    }
}
fn combined() -> Fixture {
    let beta = 0.2;
    let arm = Vec3::new(0.0, -beta * RADIUS, 0.0);
    let pose = Pose::new(arm.scale(-1.0), UnitQuaternion::IDENTITY).unwrap();
    let properties = [2.0, 3.0].map(|mass| {
        // Hertz rim traction has F=4 M/a. These original inertias give
        // W*[F,0,M/a] parallel to [1,0,1], including its off-diagonal terms.
        let inertia = mass * RADIUS * RADIUS * (1.0 + 3.0 * beta - 4.0 * beta * beta) / 4.0;
        MassProperties::new(mass, Vec3::ZERO, Vec3::new(inertia, inertia, inertia)).unwrap()
    });
    Fixture {
        bodies: [
            state(
                pose,
                properties[0],
                Vec3::new((1.0 - beta) * RADIUS * 3.0, 0.0, 0.0),
                3.0,
            ),
            state(pose, properties[1], Vec3::ZERO, 0.0),
        ],
        properties,
        binding: RigidContactBinding {
            body_a: 0,
            body_b: 1,
            arm_a_body: arm,
            arm_b_body: arm,
        },
        frame: TangentFrame::new([0.0, 0.0, 1.0], [1.0, 0.0, 0.0]).unwrap(),
    }
}

#[test]
fn sliding_crosses_stop_with_correct_point_recoil_and_zero_dwell_work() {
    let f = axial(1.0, 0.0);
    let law = law();
    let mobility = 1.0 / 2.0 + 1.0 / 3.0 + 0.2 * 0.2 / 0.10 + 0.15 * 0.15 / 0.20;
    let stop = 1.0 / (4.0 * mobility);
    let mut accepted = None;
    for duration in [1.5 * stop, 5.0 * stop] {
        let mut bodies = f.bodies;
        let report = apply_rigid_hertz_trajectory(
            &mut bodies,
            &f.properties,
            f.contact(&law),
            duration,
            controls(),
            &CancelGate::new(),
        )
        .unwrap();
        close(report.relaxation.kinetic_phase_end_s.unwrap(), stop);
        close(report.relaxation.stationary_duration_s, duration - stop);
        close(report.relaxation.slip_increment_m[0], 0.5 * stop);
        close(report.relaxation.dissipated_energy_j, 0.5 / mobility);
        close(report.kinetic_change_j, -0.5 / mobility);
        close(report.slip_velocity_after_m_s[0], 0.0);
        close(bodies[0].linear_momentum_world().x, 2.0 - 4.0 * stop);
        close(bodies[1].linear_momentum_world().x, 4.0 * stop);
        close(bodies[0].angular_momentum_body().y, 0.2 * 4.0 * stop);
        close(bodies[1].angular_momentum_body().y, 0.15 * 4.0 * stop);
        close_vec(report.linear_momentum_change_ns, Vec3::ZERO);
        close_vec(report.angular_momentum_change_nms, Vec3::ZERO);
        for i in 0..2 {
            assert_eq!(bodies[i].pose(), f.bodies[i].pose());
        }
        if let Some(previous) = accepted {
            assert_eq!(bodies, previous);
        }
        accepted = Some(bodies);
    }
}

#[test]
fn spin_shares_angular_momentum_and_heats_only_relative_rotation() {
    let f = axial(0.0, 0.3);
    let law = law();
    let ia = 0.12;
    let ib = 0.25;
    let inverse_inertia = 1.0 / ia + 1.0 / ib;
    let moment = (3.0 * PI / 16.0) * 4.0 * RADIUS;
    let stop = 0.3 / (inverse_inertia * moment);
    let mut bodies = f.bodies;
    let report = apply_rigid_hertz_trajectory(
        &mut bodies,
        &f.properties,
        f.contact(&law),
        5.0 * stop,
        controls(),
        &CancelGate::new(),
    )
    .unwrap();
    let common = ia * 0.3 / (ia + ib);
    let heat = 0.5 * ia * ib / (ia + ib) * 0.3 * 0.3;
    close(report.relaxation.kinetic_phase_end_s.unwrap(), stop);
    close(report.relaxation.spin_increment_rad, 0.5 * 0.3 * stop);
    close(report.relaxation.dissipated_energy_j, heat);
    close(report.kinetic_change_j, -heat);
    close(bodies[0].angular_momentum_body().z / ia, common);
    close(bodies[1].angular_momentum_body().z / ib, common);
    close(
        bodies[0].angular_momentum_body().z + bodies[1].angular_momentum_body().z,
        ia * 0.3,
    );
    assert_eq!(report.force_impulse_world_ns, Vec3::ZERO);
    close_vec(
        report.spin_angular_impulse_world_nms,
        Vec3::new(0.0, 0.0, -moment * stop),
    );
    close(report.total_energy_residual_j, 0.0);
}

#[test]
fn combined_motion_uses_native_cross_response_and_rotates_the_tangent_gauge() {
    let f = combined();
    let law = law();
    let force = 4.0 * 9.0 * PI / 32.0;
    let moment = 4.0 * RADIUS * 9.0 * PI / 128.0;
    let inverse_mass = 1.0 / 2.0 + 1.0 / 3.0;
    let inverse_inertia = 1.0 / f.properties[0].principal_inertia_body().z
        + 1.0 / f.properties[1].principal_inertia_body().z;
    let arm = 0.2 * RADIUS;
    let deceleration =
        (inverse_mass + arm * arm * inverse_inertia) * force + arm * inverse_inertia * moment;
    let stop = RADIUS * 3.0 / deceleration;
    let heat = 0.5 * (force * RADIUS * 3.0 + moment * 3.0) * stop;
    let mut original: Option<[RigidBodyState; 2]> = None;
    for frame in [
        f.frame,
        TangentFrame::new([0.0, 0.0, 1.0], [0.0, 1.0, 0.0]).unwrap(),
    ] {
        let candidate = prepare_rigid_hertz_trajectory(
            &f.bodies,
            &f.properties,
            RigidHertzContact {
                frame,
                ..f.contact(&law)
            },
            3.0 * stop,
            controls(),
            &CancelGate::new(),
        )
        .unwrap();
        let report = candidate.report;
        close(report.relaxation.kinetic_phase_end_s.unwrap(), stop);
        close(report.relaxation.dissipated_energy_j, heat);
        close(report.kinetic_change_j, -heat);
        close_vec(
            report.force_impulse_world_ns,
            Vec3::new(-force * stop, 0.0, 0.0),
        );
        close_vec(
            report.spin_angular_impulse_world_nms,
            Vec3::new(0.0, 0.0, -moment * stop),
        );
        assert!(
            report.scaled_mobility_per_kg[0][2]
                .abs()
                .max(report.scaled_mobility_per_kg[1][2].abs())
                > 0.1
        );
        close(report.spin_after_rad_s, 0.0);
        close_vec(report.angular_momentum_change_nms, Vec3::ZERO);
        if let Some(previous) = original {
            for i in 0..2 {
                close_vec(
                    candidate.bodies_after[i].linear_momentum_world(),
                    previous[i].linear_momentum_world(),
                );
                close_vec(
                    candidate.bodies_after[i].angular_momentum_body(),
                    previous[i].angular_momentum_body(),
                );
            }
        }
        original = Some(candidate.bodies_after);
    }
}

#[test]
fn split_native_interval_matches_complete_stop_and_preserves_unbound_bodies() {
    let f = combined();
    let law = law();
    let mut whole = f.bodies;
    let total = apply_rigid_hertz_trajectory(
        &mut whole,
        &f.properties,
        f.contact(&law),
        1.0,
        controls(),
        &CancelGate::new(),
    )
    .unwrap();
    let first_dt = 0.4 * total.relaxation.kinetic_phase_end_s.unwrap();
    let untouched = f.bodies[1];
    let mut split = [f.bodies[0], untouched, f.bodies[1]];
    let properties = [f.properties[0], f.properties[1], f.properties[1]];
    let contact = RigidHertzContact {
        binding: RigidContactBinding {
            body_b: 2,
            ..f.binding
        },
        ..f.contact(&law)
    };
    let first = apply_rigid_hertz_trajectory(
        &mut split,
        &properties,
        contact,
        first_dt,
        controls(),
        &CancelGate::new(),
    )
    .unwrap();
    let second = apply_rigid_hertz_trajectory(
        &mut split,
        &properties,
        contact,
        1.0 - first_dt,
        controls(),
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(split[1], untouched);
    close(
        first.relaxation.dissipated_energy_j + second.relaxation.dissipated_energy_j,
        total.relaxation.dissipated_energy_j,
    );
    close(
        first.relaxation.slip_increment_m[0] + second.relaxation.slip_increment_m[0],
        total.relaxation.slip_increment_m[0],
    );
    close(
        first.relaxation.spin_increment_rad + second.relaxation.spin_increment_rad,
        total.relaxation.spin_increment_rad,
    );
    for (i, index) in [(0, 0), (1, 2)] {
        close_vec(
            split[index].linear_momentum_world(),
            whole[i].linear_momentum_world(),
        );
        close_vec(
            split[index].angular_momentum_body(),
            whole[i].angular_momentum_body(),
        );
    }
}

#[test]
fn curved_motion_or_cancelled_phase_leaves_the_native_pair_unchanged() {
    let f = axial(0.1, 3.0);
    let law = law();
    let mut bodies = f.bodies;
    assert!(matches!(
        apply_rigid_hertz_trajectory(
            &mut bodies,
            &f.properties,
            f.contact(&law),
            1.0,
            controls(),
            &CancelGate::new()
        ),
        Err(RigidHertzTrajectoryError::Relaxation(
            HertzRelaxationError::CurvedTrajectory { .. }
        ))
    ));
    assert_eq!(bodies, f.bodies);
    let f = combined();
    let mut bodies = f.bodies;
    let gate = CancelGate::new();
    gate.request();
    assert!(
        apply_rigid_hertz_trajectory(
            &mut bodies,
            &f.properties,
            f.contact(&law),
            1.0,
            controls(),
            &gate
        )
        .is_err()
    );
    assert_eq!(bodies, f.bodies);
    let mut limited = controls();
    limited.contact.maximum_spin_moment_nm = 0.0;
    assert!(matches!(
        apply_rigid_hertz_trajectory(
            &mut bodies,
            &f.properties,
            f.contact(&law),
            1.0,
            limited,
            &CancelGate::new()
        ),
        Err(RigidHertzTrajectoryError::Relaxation(
            HertzRelaxationError::EffortBudget
        ))
    ));
    assert_eq!(bodies, f.bodies);
    apply_rigid_hertz_trajectory(
        &mut bodies,
        &f.properties,
        f.contact(&law),
        1.0,
        controls(),
        &CancelGate::new(),
    )
    .unwrap();
}

#[test]
fn prepared_coefficient_trials_reuse_original_pair_and_refuse_changed_patch_geometry() {
    let f = axial(1.0, 0.0);
    let original = law();
    let gate = CancelGate::new();
    let prepared = PreparedRigidHertzTrajectory::new(
        &f.bodies,
        &f.properties,
        f.contact(&original),
        0.02,
        controls(),
        &gate,
    )
    .unwrap();
    for coefficient in [0.2, 0.6, 0.2] {
        let candidate_law = HertzWrenchContact::new(HertzWrenchParameters {
            kinetic_mu: coefficient,
            ..original.parameters()
        })
        .unwrap();
        let candidate = prepared.advance(&candidate_law, &gate).unwrap();
        let direct = prepare_rigid_hertz_trajectory(
            &f.bodies,
            &f.properties,
            f.contact(&candidate_law),
            0.02,
            controls(),
            &gate,
        )
        .unwrap();
        assert_eq!(candidate, direct);
    }
    for parameters in [
        HertzWrenchParameters {
            normal_force_n: 11.0,
            ..original.parameters()
        },
        HertzWrenchParameters {
            contact_radius_m: RADIUS * 1.1,
            ..original.parameters()
        },
    ] {
        assert!(
            prepared
                .advance(&HertzWrenchContact::new(parameters).unwrap(), &gate)
                .is_err()
        );
    }
    assert_eq!(prepared.contact_points_world(), [Vec3::ZERO; 2]);
    assert_eq!(f.bodies, axial(1.0, 0.0).bodies);
}
