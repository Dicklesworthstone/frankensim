//! G1/G3: curved original Hertz traction, finite recoil and physical stopping.
use fs_exec::CancelGate;
use fs_mbd::tangential_friction::RigidContactBinding;
use fs_mbd::tangential_friction::hertz::transient::*;
use fs_mbd::tangential_friction::hertz::{RigidHertzContact, RigidHertzControls};
use fs_mbd::{MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3};
use fs_tribo::partial_slip::TangentFrame;
use fs_tribo::sliding_patch::contact::{
    HertzWrenchConfig, HertzWrenchContact, HertzWrenchParameters,
};
use fs_tribo::sliding_patch::{HertzSlidingPatch, HertzSlidingPatchParameters};

const R: f64 = 0.03;
fn law() -> HertzWrenchContact {
    HertzWrenchContact::new(HertzWrenchParameters {
        normal_force_n: 10.0,
        contact_radius_m: R,
        static_mu: 0.8,
        kinetic_mu: 0.4,
    })
    .unwrap()
}
fn controls() -> RigidHertzControls {
    RigidHertzControls {
        maximum_bodies: 4,
        geometry_tolerance_m: 2e-12,
        velocity_tolerance_m_s: 2e-11,
        angular_velocity_tolerance_rad_s: 2e-10,
        linear_momentum_tolerance_ns: 2e-11,
        angular_momentum_tolerance_nms: 2e-12,
        energy_tolerance_j: 2e-9,
        contact: HertzWrenchConfig {
            maximum_static_iterations: 64,
            maximum_iterations: 128,
            maximum_evaluations: 256,
            force_tolerance_n: 2e-12,
            moment_tolerance_nm: 2e-12,
            displacement_tolerance_m: 2e-12,
            energy_tolerance_j: 1e-10,
            relative_matrix_tolerance: 1e-13,
            maximum_force_n: 20.0,
            maximum_spin_moment_nm: 2.0,
        },
    }
}
fn transient() -> HertzTransientControls {
    HertzTransientControls {
        maximum_step_s: 0.005,
        minimum_step_s: 1e-16,
        maximum_attempts: 50_000,
        maximum_evaluations: 200_000,
        scaled_rate_tolerance_m_s: 1e-9,
        scaled_impulse_tolerance_ns: 1e-10,
        scaled_travel_tolerance_m: 1e-10,
        heat_tolerance_j: 1e-11,
        stopping_time_tolerance_s: 1e-8,
    }
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:.17e} vs {b:.17e}, tolerance {tolerance:.3e}"
    );
}
fn close_vec(a: Vec3, b: Vec3, tolerance: f64) {
    close(a.x, b.x, tolerance);
    close(a.y, b.y, tolerance);
    close(a.z, b.z, tolerance);
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
    // Independent planar point/free-couple response: each body contributes
    // diag(1/m,1/m,0) + [-r_y,r_x,a] [-r_y,r_x,a]^T / I_z.
    fn mobility(&self) -> [[f64; 3]; 3] {
        let mut w = [[0.0; 3]; 3];
        for (index, arm) in [self.binding.arm_a_body, self.binding.arm_b_body]
            .into_iter()
            .enumerate()
        {
            let v = [-arm.y, arm.x, R];
            let p = self.properties[index];
            for i in 0..3 {
                for j in 0..3 {
                    w[i][j] += v[i] * v[j] / p.principal_inertia_body().z;
                }
            }
            w[0][0] += p.mass().recip();
            w[1][1] += p.mass().recip();
        }
        w
    }
}
fn fixture() -> Fixture {
    let arms = [
        Vec3::new(0.2, -0.3, 0.0).scale(R),
        Vec3::new(-0.15, 0.1, 0.0).scale(R),
    ];
    let properties = [2.0, 3.0].map(|mass| {
        MassProperties::new(
            mass,
            Vec3::ZERO,
            Vec3::new(0.7, 0.9, 1.1).scale(mass * R * R),
        )
        .unwrap()
    });
    let poses = arms.map(|arm| Pose::new(arm.scale(-1.0), UnitQuaternion::IDENTITY).unwrap());
    let spin = 3.0;
    let omega = Vec3::new(0.0, 0.0, spin);
    let com = Vec3::new(0.8, -0.3, 0.0)
        .scale(R * spin)
        .sub(omega.cross(arms[0]));
    Fixture {
        bodies: [
            RigidBodyState::new(
                poses[0],
                com.scale(properties[0].mass()),
                Vec3::new(0.0, 0.0, spin * properties[0].principal_inertia_body().z),
            )
            .unwrap(),
            RigidBodyState::new(poses[1], Vec3::ZERO, Vec3::ZERO).unwrap(),
        ],
        properties,
        binding: RigidContactBinding {
            body_a: 0,
            body_b: 1,
            arm_a_body: arms[0],
            arm_b_body: arms[1],
        },
        frame: TangentFrame::new([0.0, 0.0, 1.0], [1.0, 0.0, 0.0]).unwrap(),
    }
}
fn kinetic(bodies: [RigidBodyState; 2], properties: [MassProperties; 2]) -> f64 {
    (0..2)
        .map(|i| {
            let p = bodies[i].linear_momentum_world();
            let l = bodies[i].angular_momentum_body();
            let inertia = properties[i].principal_inertia_body();
            0.5 * (p.dot(p) / properties[i].mass()
                + l.x * l.x / inertia.x
                + l.y * l.y / inertia.y
                + l.z * l.z / inertia.z)
        })
        .sum()
}
fn solve(w: [[f64; 3]; 3], q: [f64; 3]) -> [f64; 3] {
    // Cofactor reference is independent of production's metric/Cholesky solve.
    let a = w[0][0];
    let b = w[0][1];
    let c = w[0][2];
    let d = w[1][1];
    let e = w[1][2];
    let f = w[2][2];
    let cof = [
        [d * f - e * e, c * e - b * f, b * e - c * d],
        [c * e - b * f, a * f - c * c, b * c - a * e],
        [b * e - c * d, b * c - a * e, a * d - b * b],
    ];
    let determinant = a * cof[0][0] + b * cof[0][1] + c * cof[0][2];
    std::array::from_fn(|i| (0..3).map(|j| cof[i][j] * q[j]).sum::<f64>() / determinant)
}

// Fixed-step RK4 of the original Hertz law with independently assembled W.
// Additional coordinates integrate travel and physical power, not numerical loss.
fn reference(w: [[f64; 3]; 3], duration: f64) -> [f64; 7] {
    let patch = HertzSlidingPatch::new(HertzSlidingPatchParameters {
        normal_force_n: 10.0,
        contact_radius_m: R,
        friction_coefficient: 0.4,
    })
    .unwrap();
    let rhs = |y: [f64; 7]| {
        let wrench = patch.evaluate([y[0], y[1]], y[2] / R).unwrap();
        let f = [
            wrench.force_n[0],
            wrench.force_n[1],
            wrench.spin_moment_nm / R,
        ];
        [
            (0..3).map(|j| w[0][j] * f[j]).sum(),
            (0..3).map(|j| w[1][j] * f[j]).sum(),
            (0..3).map(|j| w[2][j] * f[j]).sum(),
            y[0],
            y[1],
            y[2],
            wrench.dissipated_power_w,
        ]
    };
    let mut y = [0.8 * R * 3.0, -0.3 * R * 3.0, R * 3.0, 0.0, 0.0, 0.0, 0.0];
    let h = duration / 4096.0;
    for _ in 0..4096 {
        let a = rhs(y);
        let b = rhs(std::array::from_fn(|i| y[i] + 0.5 * h * a[i]));
        let c = rhs(std::array::from_fn(|i| y[i] + 0.5 * h * b[i]));
        let d = rhs(std::array::from_fn(|i| y[i] + h * c[i]));
        y = std::array::from_fn(|i| y[i] + h * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0);
    }
    y
}

#[test]
fn curved_native_path_matches_independent_response_and_rk4_work_travel() {
    let f = fixture();
    let law = law();
    let mut bodies = f.bodies;
    let report = apply_rigid_hertz_transient(
        &mut bodies,
        &f.properties,
        f.contact(&law),
        0.008,
        controls(),
        transient(),
        &CancelGate::new(),
    )
    .unwrap();
    assert!(!report.transient.exact_ray && report.transient.stopping_time_bracket_s.is_none());
    let w = f.mobility();
    assert!(w[0][1] != 0.0 && w[0][2] != 0.0 && w[1][2] != 0.0);
    for i in 0..3 {
        for j in 0..3 {
            close(report.scaled_mobility_per_kg[i][j], w[i][j], 2e-15);
        }
    }
    let expected = reference(w, 0.008);
    for i in 0..2 {
        close(report.slip_velocity_m_s[1][i], expected[i], 3e-9);
        close(report.transient.slip_increment_m[i], expected[3 + i], 3e-10);
    }
    close(R * report.spin_rad_s[1], expected[2], 3e-9);
    close(R * report.transient.spin_increment_rad, expected[5], 3e-10);
    close(report.dissipated_energy_j, expected[6], 3e-10);
    close(
        kinetic(bodies, f.properties) - kinetic(f.bodies, f.properties),
        -expected[6],
        3e-10,
    );
    assert_eq!(report.audit.normal_velocity_bound_m_s, 0.0);
    let initial = [0.8 * R * 3.0, -0.3 * R * 3.0, R * 3.0];
    let expected_impulse = solve(w, std::array::from_fn(|i| expected[i] - initial[i]));
    close_vec(
        report.force_impulse_world_ns,
        Vec3::new(expected_impulse[0], expected_impulse[1], 0.0),
        3e-9,
    );
    close(
        report.spin_angular_impulse_world_nms.z,
        R * expected_impulse[2],
        3e-10,
    );
    assert_ne!(bodies[1].linear_momentum_world(), Vec3::ZERO);
    for i in 0..2 {
        assert_eq!(bodies[i].pose(), f.bodies[i].pose());
    }
}

#[test]
fn curved_arrest_loses_only_relative_energy_and_stationary_dwell_adds_nothing() {
    let f = fixture();
    let law = law();
    let gate = CancelGate::new();
    let q = [0.8 * R * 3.0, -0.3 * R * 3.0, R * 3.0];
    let response = solve(f.mobility(), q);
    let expected_heat = 0.5 * (0..3).map(|i| q[i] * response[i]).sum::<f64>();
    let mut previous = None;
    for duration in [1.0, 10.0] {
        let mut bodies = f.bodies;
        let report = apply_rigid_hertz_transient(
            &mut bodies,
            &f.properties,
            f.contact(&law),
            duration,
            controls(),
            transient(),
            &gate,
        )
        .unwrap();
        let [lo, hi] = report.transient.stopping_time_bracket_s.unwrap();
        assert!(!report.transient.exact_ray && lo < hi && hi < duration);
        assert!(hi - lo <= transient().stopping_time_tolerance_s);
        assert!(report.transient.stopping_tail_heat_j > 0.0);
        close(report.dissipated_energy_j, expected_heat, 2e-10);
        close(report.kinetic_change_j, -expected_heat, 2e-10);
        for velocity in report.slip_velocity_m_s[1] {
            close(velocity, 0.0, controls().velocity_tolerance_m_s);
        }
        close(
            report.spin_rad_s[1],
            0.0,
            controls().angular_velocity_tolerance_rad_s,
        );
        close_vec(
            report.force_impulse_world_ns,
            Vec3::new(-response[0], -response[1], 0.0),
            2e-10,
        );
        close(
            report.spin_angular_impulse_world_nms.z,
            -R * response[2],
            2e-11,
        );
        assert!(kinetic(bodies, f.properties) > 0.0);
        close_vec(report.audit.linear_momentum_change_ns, Vec3::ZERO, 2e-12);
        close_vec(report.audit.angular_momentum_change_nms, Vec3::ZERO, 2e-12);
        if let Some((previous_bodies, previous_heat)) = previous {
            assert_eq!(bodies, previous_bodies);
            assert_eq!(report.dissipated_energy_j, previous_heat);
        }
        previous = Some((bodies, report.dissipated_energy_j));
    }
}

#[test]
fn tangent_gauge_and_world_rotation_preserve_native_recoil_and_heat() {
    let f = fixture();
    let law = law();
    let gate = CancelGate::new();
    let base = prepare_rigid_hertz_transient(
        &f.bodies,
        &f.properties,
        f.contact(&law),
        0.008,
        controls(),
        transient(),
        &gate,
    )
    .unwrap();
    let changed_gauge = prepare_rigid_hertz_transient(
        &f.bodies,
        &f.properties,
        RigidHertzContact {
            frame: TangentFrame::new([0.0, 0.0, 1.0], [0.0, 1.0, 0.0]).unwrap(),
            ..f.contact(&law)
        },
        0.008,
        controls(),
        transient(),
        &gate,
    )
    .unwrap();
    let rotation = UnitQuaternion::from_axis_angle(Vec3::new(0.0, 1.0, 0.0), 0.7).unwrap();
    let rotated = f.bodies.map(|body| {
        RigidBodyState::new(
            Pose::new(
                rotation.rotate_body_to_world(body.pose().position_world()),
                rotation,
            )
            .unwrap(),
            rotation.rotate_body_to_world(body.linear_momentum_world()),
            body.angular_momentum_body(),
        )
        .unwrap()
    });
    let array = |v: Vec3| [v.x, v.y, v.z];
    let result = prepare_rigid_hertz_transient(
        &rotated,
        &f.properties,
        RigidHertzContact {
            frame: TangentFrame::new(
                array(rotation.rotate_body_to_world(Vec3::new(0.0, 0.0, 1.0))),
                array(rotation.rotate_body_to_world(Vec3::new(1.0, 0.0, 0.0))),
            )
            .unwrap(),
            ..f.contact(&law)
        },
        0.008,
        controls(),
        transient(),
        &gate,
    )
    .unwrap();
    for i in 0..2 {
        close_vec(
            changed_gauge.bodies_after[i].linear_momentum_world(),
            base.bodies_after[i].linear_momentum_world(),
            2e-10,
        );
        close_vec(
            changed_gauge.bodies_after[i].angular_momentum_body(),
            base.bodies_after[i].angular_momentum_body(),
            2e-11,
        );
        close_vec(
            result.bodies_after[i].linear_momentum_world(),
            rotation.rotate_body_to_world(base.bodies_after[i].linear_momentum_world()),
            2e-10,
        );
        close_vec(
            result.bodies_after[i].angular_momentum_body(),
            base.bodies_after[i].angular_momentum_body(),
            2e-11,
        );
    }
    close(
        result.report.dissipated_energy_j,
        base.report.dissipated_energy_j,
        2e-11,
    );
    close(
        changed_gauge.report.dissipated_energy_j,
        base.report.dissipated_energy_j,
        2e-11,
    );
}

#[test]
fn prepared_trials_and_late_budget_normal_cancellation_refusals_preserve_bodies() {
    let f = fixture();
    let law = law();
    let gate = CancelGate::new();
    let prepared = PreparedRigidHertzTransient::new(
        &f.bodies,
        &f.properties,
        f.contact(&law),
        0.008,
        controls(),
        transient(),
        &gate,
    )
    .unwrap();
    let first = prepared.advance(&law, &gate).unwrap();
    let changed = HertzWrenchContact::new(HertzWrenchParameters {
        kinetic_mu: 0.2,
        ..law.parameters()
    })
    .unwrap();
    let second = prepared.advance(&changed, &gate).unwrap();
    assert!(second.report.dissipated_energy_j < first.report.dissipated_energy_j);
    assert_eq!(first, prepared.advance(&law, &gate).unwrap());
    assert_eq!(prepared.contact_points_world(), [Vec3::ZERO; 2]);
    let mut bodies = [f.bodies[0], f.bodies[1], f.bodies[1]];
    let properties = [f.properties[0], f.properties[1], f.properties[1]];
    let saved = bodies;
    let mut limited = transient();
    limited.maximum_evaluations = 2;
    assert!(
        apply_rigid_hertz_transient(
            &mut bodies,
            &properties,
            f.contact(&law),
            0.008,
            controls(),
            limited,
            &gate
        )
        .is_err()
    );
    assert_eq!(bodies, saved);
    gate.request();
    assert!(
        apply_rigid_hertz_transient(
            &mut bodies,
            &properties,
            f.contact(&law),
            0.008,
            controls(),
            transient(),
            &gate
        )
        .is_err()
    );
    assert_eq!(bodies, saved);
    let gate = CancelGate::new();
    apply_rigid_hertz_transient(
        &mut bodies,
        &properties,
        f.contact(&law),
        0.008,
        controls(),
        transient(),
        &gate,
    )
    .unwrap();
    assert_eq!(bodies[2], saved[2]);
    let mut lifted = f.bodies;
    let arm = f.binding.arm_a_body.add(Vec3::new(0.0, 0.0, R));
    lifted[0] = RigidBodyState::new(
        Pose::new(arm.scale(-1.0), UnitQuaternion::IDENTITY).unwrap(),
        f.bodies[0].linear_momentum_world(),
        f.bodies[0].angular_momentum_body(),
    )
    .unwrap();
    let before = lifted;
    assert!(matches!(
        apply_rigid_hertz_transient(
            &mut lifted,
            &f.properties,
            RigidHertzContact {
                binding: RigidContactBinding {
                    arm_a_body: arm,
                    ..f.binding
                },
                ..f.contact(&law)
            },
            0.008,
            controls(),
            transient(),
            &gate
        ),
        Err(RigidHertzTransientError::Invalid(_))
    ));
    assert_eq!(lifted, before);
}
