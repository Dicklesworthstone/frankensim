use super::*;
use fs_mbd::UnitQuaternion;

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e} (tol={tolerance:e})"
    );
}
fn near_vec(actual: Vec3, expected: Vec3, tolerance: f64) {
    for (a, b) in [
        (actual.x, expected.x),
        (actual.y, expected.y),
        (actual.z, expected.z),
    ] {
        near(a, b, tolerance);
    }
}
fn fixture() -> (Vec<Cell>, Body) {
    let rotation = UnitQuaternion::from_axis_angle(Vec3::new(1.0, -2.0, 0.7), 0.63).unwrap();
    let axes = [
        Vec3::new(1.0, 0.0, 0.0),
        Vec3::new(0.0, 1.0, 0.0),
        Vec3::new(0.0, 0.0, 1.0),
    ]
    .map(|axis| rotation.rotate_body_to_world(axis));
    let body = Body {
        velocity: Vec3::new(0.11, -0.07, 0.09),
        omega: Vec3::new(-0.13, 0.17, 0.03),
        inverse_sqrt_mass: 1.0 / 2.3_f64.sqrt(),
        inverse_sqrt_inertia: [0.31_f64, 0.42, 0.55].map(|i| 1.0 / i.sqrt()),
        principal_axes: axes,
    };
    let cells = [
        (
            0.27,
            Vec3::new(-0.7, 0.2, 0.4),
            Vec3::new(-0.03, 0.05, 0.07),
        ),
        (
            0.53,
            Vec3::new(0.3, -0.4, 0.1),
            Vec3::new(0.09, 0.02, -0.01),
        ),
        (
            0.18,
            Vec3::new(0.4, 0.5, -0.3),
            Vec3::new(-0.08, -0.09, 0.06),
        ),
        (
            0.46,
            Vec3::new(-0.1, 0.1, -0.6),
            Vec3::new(0.12, -0.04, -0.05),
        ),
        (0.11, Vec3::new(0.2, 0.2, 0.3), Vec3::new(-0.04, 0.01, 0.08)),
    ]
    .map(|(mass, arm, velocity)| Cell {
        mass,
        arm,
        velocity,
    })
    .to_vec();
    (cells, body)
}
fn inverse_inertia(body: Body, angular: Vec3) -> Vec3 {
    (0..3).fold(Vec3::ZERO, |sum, k| {
        sum.add(
            body.principal_axes[k]
                .scale(body.principal_axes[k].dot(angular) * body.inverse_sqrt_inertia[k].powi(2)),
        )
    })
}
fn after(cells: &[Cell], body: Body, plan: &Plan) -> (Vec<Cell>, Body) {
    let mut endpoint = body;
    endpoint.velocity = body
        .velocity
        .add(plan.body_linear.scale(body.inverse_sqrt_mass.powi(2)));
    endpoint.omega = body.omega.add(inverse_inertia(body, plan.body_angular));
    let liquid = cells
        .iter()
        .zip(&plan.impulses)
        .map(|(c, j)| Cell {
            velocity: c.velocity.add(j.scale(1.0 / c.mass)),
            ..*c
        })
        .collect();
    (liquid, endpoint)
}
fn accounts(cells: &[Cell], body: Body) -> (Vec3, Vec3, f64) {
    let mass = 1.0 / body.inverse_sqrt_mass.powi(2);
    let mut p = body.velocity.scale(mass);
    let mut l = Vec3::ZERO;
    let mut energy = 0.5 * mass * body.velocity.dot(body.velocity);
    for k in 0..3 {
        let omega = body.omega.dot(body.principal_axes[k]);
        let moment = 1.0 / body.inverse_sqrt_inertia[k].powi(2);
        l = l.add(body.principal_axes[k].scale(moment * omega));
        energy += 0.5 * moment * omega * omega;
    }
    for cell in cells {
        p = p.add(cell.velocity.scale(cell.mass));
        l = l.add(cell.arm.cross(cell.velocity.scale(cell.mass)));
        energy += 0.5 * cell.mass * cell.velocity.dot(cell.velocity);
    }
    (p, l, energy)
}

#[test]
fn two_finite_masses_match_analytic_relative_decay() {
    let (_, mut body) = fixture();
    body.omega = Vec3::ZERO;
    let cell = Cell {
        mass: 0.7,
        arm: Vec3::ZERO,
        velocity: Vec3::new(-0.04, 0.02, 0.01),
    };
    let mass = 1.0 / body.inverse_sqrt_mass.powi(2);
    for z in [0.0, 1e-12, 0.002, 0.7, 30.0, f64::INFINITY] {
        let plan = solve(&[cell], body, z, &CancelGate::new_clock_free()).unwrap();
        let factor =
            mass * cell.mass / (mass + cell.mass) * -(-z * (1.0 + cell.mass / mass)).exp_m1();
        let expected = body.velocity.sub(cell.velocity).scale(factor);
        near_vec(plan.impulses[0], expected, 6e-17);
        near_vec(plan.body_linear, expected.scale(-1.0), 6e-17);
        near_vec(plan.body_angular, Vec3::ZERO, 0.0);
        let (cells, end) = after(&[cell], body, &plan);
        let slip = end.velocity.sub(cells[0].velocity);
        near_vec(
            slip,
            body.velocity
                .sub(cell.velocity)
                .scale((-z * (1.0 + cell.mass / mass)).exp()),
            8e-17,
        );
    }
}

/// Independent physical force/torque RHS, with direct world-frame inverse inertia.
fn rhs(state: &[Vec3], cells: &[Cell], body: Body) -> Vec<Vec3> {
    let mut out = vec![Vec3::ZERO; state.len()];
    let mut force = Vec3::ZERO;
    let mut torque = Vec3::ZERO;
    for (i, cell) in cells.iter().enumerate() {
        let acceleration = state[0].add(state[1].cross(cell.arm)).sub(state[i + 2]);
        let f = acceleration.scale(cell.mass);
        force = force.sub(f);
        torque = torque.sub(cell.arm.cross(f));
        out[i + 2] = acceleration;
    }
    out[0] = force.scale(body.inverse_sqrt_mass.powi(2));
    out[1] = inverse_inertia(body, torque);
    out
}
fn shifted(state: &[Vec3], rate: &[Vec3], h: f64) -> Vec<Vec3> {
    state
        .iter()
        .zip(rate)
        .map(|(a, b)| a.add(b.scale(h)))
        .collect()
}
#[test]
fn full_spatial_metric_matches_independent_force_torque_rk4() {
    let (cells, body) = fixture();
    let z = 1.7;
    let mut reference = vec![body.velocity, body.omega];
    reference.extend(cells.iter().map(|c| c.velocity));
    let h = z / 4096.0;
    for _ in 0..4096 {
        let k1 = rhs(&reference, &cells, body);
        let k2 = rhs(&shifted(&reference, &k1, h / 2.0), &cells, body);
        let k3 = rhs(&shifted(&reference, &k2, h / 2.0), &cells, body);
        let k4 = rhs(&shifted(&reference, &k3, h), &cells, body);
        for i in 0..reference.len() {
            reference[i] = reference[i].add(
                k1[i]
                    .add(k2[i].scale(2.0))
                    .add(k3[i].scale(2.0))
                    .add(k4[i])
                    .scale(h / 6.0),
            );
        }
    }
    let plan = solve(&cells, body, z, &CancelGate::new_clock_free()).unwrap();
    let (liquid, end) = after(&cells, body, &plan);
    near_vec(end.velocity, reference[0], 4e-14);
    near_vec(end.omega, reference[1], 4e-14);
    for (i, c) in liquid.iter().enumerate() {
        near_vec(c.velocity, reference[i + 2], 4e-14);
    }
}
#[test]
fn exact_source_has_semigroup_and_galilean_covariance() {
    let (cells, body) = fixture();
    let gate = CancelGate::new_clock_free();
    let full = solve(&cells, body, 0.83, &gate).unwrap();
    let (full_cells, full_body) = after(&cells, body, &full);
    let first = solve(&cells, body, 0.31, &gate).unwrap();
    let (middle_cells, middle_body) = after(&cells, body, &first);
    let second = solve(&middle_cells, middle_body, 0.52, &gate).unwrap();
    let (twice_cells, twice_body) = after(&middle_cells, middle_body, &second);
    near_vec(full_body.velocity, twice_body.velocity, 5e-16);
    near_vec(full_body.omega, twice_body.omega, 5e-16);
    for (a, b) in full_cells.iter().zip(twice_cells) {
        near_vec(a.velocity, b.velocity, 5e-16);
    }
    let boost = Vec3::new(30.0, -20.0, 10.0);
    let boosted: Vec<Cell> = cells
        .iter()
        .map(|c| Cell {
            velocity: c.velocity.add(boost),
            ..*c
        })
        .collect();
    let boosted_body = Body {
        velocity: body.velocity.add(boost),
        ..body
    };
    let translated = solve(&boosted, boosted_body, 0.83, &gate).unwrap();
    for (a, b) in full.impulses.iter().zip(translated.impulses) {
        near_vec(*a, b, 8e-15);
    }
}
#[test]
fn stiff_limit_relaxes_all_slip_and_conserves_spatial_momentum() {
    let (cells, body) = fixture();
    let initial = accounts(&cells, body);
    for exposure in [0.1, 2.0, 1e8, f64::INFINITY] {
        let plan = solve(&cells, body, exposure, &CancelGate::new_clock_free()).unwrap();
        let (end_cells, end_body) = after(&cells, body, &plan);
        let end = accounts(&end_cells, end_body);
        near_vec(initial.0, end.0, 8e-16);
        near_vec(initial.1, end.1, 8e-16);
        assert!(end.2 < initial.2);
        if exposure > 100.0 {
            for c in end_cells {
                near_vec(c.velocity, end_body.point_velocity(c.arm), 6e-16);
            }
        }
    }
}
#[test]
fn rank_deficient_stencil_keeps_unsupported_spin_and_fluid_modes() {
    let (_, mut body) = fixture();
    body.principal_axes = [
        Vec3::new(1.0, 0.0, 0.0),
        Vec3::new(0.0, 1.0, 0.0),
        Vec3::new(0.0, 0.0, 1.0),
    ];
    body.velocity = Vec3::ZERO;
    body.omega = Vec3::new(0.0, 0.0, 0.17);
    let cells = [
        Cell {
            mass: 0.5,
            arm: Vec3::new(0.0, 0.0, 0.3),
            velocity: Vec3::new(0.02, 0.0, 0.0),
        },
        Cell {
            mass: 0.5,
            arm: Vec3::new(0.0, 0.0, 0.3),
            velocity: Vec3::new(-0.02, 0.0, 0.0),
        },
    ];
    let plan = solve(&cells, body, 0.7, &CancelGate::new_clock_free()).unwrap();
    let (end_cells, end) = after(&cells, body, &plan);
    near_vec(end.omega, body.omega, 0.0);
    near_vec(end.velocity, Vec3::ZERO, 0.0);
    near(end_cells[0].velocity.x, 0.02 * (-0.7_f64).exp(), 5e-18);
    near(end_cells[1].velocity.x, -0.02 * (-0.7_f64).exp(), 5e-18);
}
#[test]
fn small_exposure_response_and_null_metric_have_correct_limits() {
    for c in [0.0, 1e-18, 0.2, 5.0, 1e12] {
        let z = 1e-10 / (1.0 + c);
        let (b, f) = weights(z, c).unwrap();
        near(b / z, 1.0, 6e-11);
        near(f / (0.5 * z * z), 1.0, 7e-11);
    }
    for z in [0.2, 1.0, 50.0] {
        let (_, f) = weights(z, 0.0).unwrap();
        near(f, 1.0 - (-z).exp() * (1.0 + z), 3e-16);
    }
    let (cells, body) = fixture();
    let gate = CancelGate::new_clock_free();
    gate.request();
    assert!(matches!(
        solve(&cells, body, 1.0, &gate),
        Err(LiquidMomentumError::Cancelled)
    ));
}
