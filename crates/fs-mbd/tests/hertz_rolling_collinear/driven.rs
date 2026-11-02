use super::*;
use fs_mbd::Wrench;

fn driven(f: &Fixture, external: [Wrench; 2], dt: f64) -> RigidHertzRollingCollinearCandidate {
    PreparedRigidHertzRollingCollinear::new_driven(
        &f.bodies,
        &f.properties,
        &f.history,
        f.contact(),
        external,
        dt,
        controls(),
        &CancelGate::new(),
    )
    .unwrap()
    .advance(&f.hertz, &CancelGate::new())
    .unwrap()
}

// Independent sphere inverse inertia: a body torque changes both rolling
// and point speed; a COM force compensates that point-speed contribution.
fn acceleration(f: &Fixture, shear: f64, rolling: f64) -> [Wrench; 2] {
    let m = f.properties[0].mass();
    let inertia = 0.4 * m * R * R;
    [
        Wrench {
            force_world: Vec3::new(m * (shear + rolling), 0.0, 0.0),
            torque_body: Vec3::new(0.0, inertia * rolling / R, 0.0),
        },
        Wrench::ZERO,
    ]
}

fn balance(f: &Fixture, c: &RigidHertzRollingCollinearCandidate) {
    let r = &c.report;
    close(
        r.kinetic_after_j - r.kinetic_before_j + r.surface_heat_j + r.bulk_heat_j,
        r.external_work_j,
        2e-11,
    );
    close(r.total_energy_residual_j, 0.0, 2e-11);
    close_vec(r.linear_momentum_residual_ns, Vec3::ZERO, 2e-11);
    close_vec(r.angular_momentum_residual_nms, Vec3::ZERO, 2e-11);
    assert!(r.surface_heat_j >= 0.0 && r.bulk_heat_j >= 0.0);
    for (before, after) in f.bodies.iter().zip(c.bodies_after) {
        assert_eq!(before.pose(), after.pose());
    }
}

#[test]
fn zero_drive_reduces_exactly_and_linear_drive_has_independent_work_and_heat() {
    let f = fixture(0.0, [0.3, 0.2], 0.2, 0.04);
    assert_eq!(f.run(0.01), driven(&f, [Wrench::ZERO; 2], 0.01));
    let dt = 0.03;
    let external = acceleration(&f, 0.7, 0.2);
    let c = driven(&f, external, dt);
    balance(&f, &c);
    let [[a, b], _] = f.response();
    let force = -0.2 * f.patch.normal_load_n();
    let travel = 0.2 * dt + 0.5 * (a * force + 0.7) * dt * dt;
    close(c.report.phases[0].shear_displacement_m, travel, 2e-15);
    close(
        c.report.phases[0].shear_velocity_after_m_s,
        0.2 + (a * force + 0.7) * dt,
        2e-15,
    );
    close(
        c.report.phases[0].rolling_rate_after_rad_s,
        0.04 + (b * force + 0.2) * dt / R,
        2e-14,
    );
    close(c.report.surface_heat_j, -force * travel, 2e-15);
    close(c.report.bulk_heat_j, 0.0, 0.0);
    // With zero rolling viscosity all body accelerations are constant, so
    // endpoint mean COM/angular velocities give exact continuous external work.
    let mut work = 0.0;
    for i in 0..2 {
        let old = f.bodies[i]
            .point_kinematics(f.properties[i], Vec3::ZERO)
            .unwrap();
        let new = c.bodies_after[i]
            .point_kinematics(f.properties[i], Vec3::ZERO)
            .unwrap();
        work += 0.5
            * dt
            * (external[i].force_world.dot(
                old.center_of_mass_velocity_world
                    .add(new.center_of_mass_velocity_world),
            ) + external[i]
                .torque_body
                .dot(old.angular_velocity_body.add(new.angular_velocity_body)));
    }
    close(c.report.external_work_j, work, 2e-15);
}

#[test]
fn drive_finds_an_interior_stop_even_when_unresolved_endpoint_is_sliding_positive() {
    let f = fixture(0.01, [0.002, 0.001], 0.00002, -0.02);
    let [[a, b], [_b, c]] = f.response();
    let force = -0.001 * f.patch.normal_load_n();
    let gamma = f.tau * f.patch.normal_load_n() / R;
    let du = 0.0005 - a * force;
    let dr = -b * force;
    let lambda = gamma * c;
    let unresolved =
        |t: f64| 0.00002 + 0.0005 * t - b * gamma * (-0.002) * (-(-lambda * t).exp_m1()) / lambda;
    assert!(unresolved(5.0) > 0.0 && unresolved(0.1) < 0.0);
    let mut lo = 0.0;
    let mut hi = 0.1;
    for _ in 0..80 {
        let mid = (lo + hi) * 0.5;
        if unresolved(mid) > 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let result = driven(&f, acceleration(&f, du, dr), 5.0);
    balance(&f, &result);
    assert!(result.report.phases[0].reached_stop);
    close(result.report.phases[0].duration_s, (lo + hi) * 0.5, 2e-13);
    assert!(result.history.is_stopped());
    assert_eq!(
        result.report.phases[1].mode,
        RigidHertzRollingCollinearMode::Sticking
    );
}

#[test]
fn motor_drive_crosses_static_capacity_and_split_windows_preserve_the_event() {
    let mut f = fixture(0.01, [0.00005, 0.00004], 0.0, 0.0);
    let external = acceleration(&f, 0.0, 0.002);
    let [[a, b], [_b, c]] = f.response();
    let gamma = f.tau * f.patch.normal_load_n() / R;
    let decay = gamma * (c - b * b / a);
    let equilibrium = 0.002 / decay;
    let limit = a * 0.00005 * f.patch.normal_load_n() / (-b * gamma);
    let exit = -(1.0 - limit / equilibrium).ln() / decay;
    let whole = driven(&f, external, 1.2);
    balance(&f, &whole);
    assert_eq!(whole.report.phases.len(), 2);
    assert!(whole.report.phases[0].reached_static_limit);
    close(whole.report.phases[0].duration_s, exit, 2e-13);
    close(whole.report.phases[0].surface_heat_j, 0.0, 0.0);
    assert_eq!(
        whole.report.phases[1].mode,
        RigidHertzRollingCollinearMode::SlidingPositive
    );
    assert!(!whole.history.is_stopped());
    let split = whole.report.phases[0].duration_s;
    let first = driven(&f, external, split);
    f.bodies = first.bodies_after;
    f.history = first.history;
    let second = driven(&f, external, 1.2 - split);
    balance(&f, &second);
    for i in 0..2 {
        close_vec(
            second.bodies_after[i].linear_momentum_world(),
            whole.bodies_after[i].linear_momentum_world(),
            2e-12,
        );
        close_vec(
            second.bodies_after[i].angular_momentum_body(),
            whole.bodies_after[i].angular_momentum_body(),
            2e-12,
        );
    }
    close(
        first.report.surface_heat_j + second.report.surface_heat_j,
        whole.report.surface_heat_j,
        2e-13,
    );
    close(
        first.report.bulk_heat_j + second.report.bulk_heat_j,
        whole.report.bulk_heat_j,
        2e-13,
    );
    close(
        first.report.external_work_j + second.report.external_work_j,
        whole.report.external_work_j,
        2e-13,
    );
}

#[test]
fn common_acceleration_changes_external_work_without_changing_friction_heat() {
    let f = fixture(0.01, [0.6, 0.1], 0.004, 0.01);
    let dt = 0.4;
    let free = f.run(dt);
    let external = f.properties.map(|p| Wrench {
        force_world: Vec3::new(0.0, p.mass() * 0.3, 0.0),
        torque_body: Vec3::ZERO,
    });
    let result = driven(&f, external, dt);
    balance(&f, &result);
    close(
        result.report.surface_heat_j,
        free.report.surface_heat_j,
        2e-16,
    );
    close(result.report.bulk_heat_j, free.report.bulk_heat_j, 2e-16);
    close(
        result.report.external_work_j,
        0.5 * (f.properties[0].mass() + f.properties[1].mass()) * (0.3 * dt).powi(2),
        2e-14,
    );
    let invalid = [
        Wrench {
            force_world: Vec3::new(0.0, 0.0, 1.0),
            torque_body: Vec3::ZERO,
        },
        Wrench::ZERO,
    ];
    assert!(
        PreparedRigidHertzRollingCollinear::new_driven(
            &f.bodies,
            &f.properties,
            &f.history,
            f.contact(),
            invalid,
            dt,
            controls(),
            &CancelGate::new(),
        )
        .is_err()
    );
    let gate = CancelGate::new();
    gate.request();
    assert!(
        PreparedRigidHertzRollingCollinear::new_driven(
            &f.bodies,
            &f.properties,
            &f.history,
            f.contact(),
            external,
            dt,
            controls(),
            &gate,
        )
        .is_err()
    );
}

#[test]
fn driven_tiny_slip_and_equal_static_kinetic_departure_retain_physical_events() {
    for tau in [0.0, 0.01] {
        let tiny = fixture(tau, [0.8, 0.1], 1e-20, 0.0);
        let force = -0.1 * tiny.patch.normal_load_n();
        let drive = 0.05;
        let rolling_drive = -tiny.response()[1][0] * force;
        let expected_acceleration = tiny.response()[0][0] * force + drive;
        let result = driven(&tiny, acceleration(&tiny, drive, rolling_drive), 0.001);
        assert!(result.report.phases[0].reached_stop);
        close(
            result.report.phases[0].duration_s,
            -1e-20 / expected_acceleration,
            1e-33,
        );
        let expected_distance = -0.5e-40 / expected_acceleration;
        close(
            result.report.surface_heat_j,
            -force * expected_distance,
            1e-52,
        );
    }
    let equal = fixture(0.01, [0.00005; 2], 0.0, 0.0);
    let result = driven(&equal, acceleration(&equal, 0.0, 0.002), 1.2);
    balance(&equal, &result);
    assert!(result.report.phases[0].reached_static_limit);
    assert_eq!(
        result.report.phases[1].mode,
        RigidHertzRollingCollinearMode::SlidingPositive
    );
    assert!(result.report.phases[1].shear_velocity_after_m_s > 0.0);
}
