//! G0/G3 native panel/wake constraints, circulation convention and similarity.

use fs_bem::BemError;
use fs_bem::panel2d::{Airfoil2d, PreparedPanel2d, naca4_symmetric, solve};

fn vortex_velocity(point: [f64; 2], source: [f64; 2], gamma: f64, core: f64) -> [f64; 2] {
    let dx = point[0] - source[0];
    let dy = point[1] - source[1];
    let factor = gamma / (std::f64::consts::TAU * (dx * dx + dy * dy + core * core));
    [factor * dy, -factor * dx]
}

fn assert_close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * (1.0 + expected.abs()),
        "actual={actual:.17e}, expected={expected:.17e}"
    );
}

#[test]
fn g0_retained_uniform_solve_matches_existing_hess_smith() {
    let foil = naca4_symmetric(0.12, 64).unwrap();
    let prepared = PreparedPanel2d::new(&foil).unwrap();
    assert_eq!(prepared.nodes(), foil.nodes());
    for alpha in [-0.13, 0.0, 0.09] {
        let uniform = [fs_math::det::cos(alpha), fs_math::det::sin(alpha)];
        let state = prepared
            .solve_incident(&vec![uniform; prepared.midpoints().len()])
            .unwrap();
        let existing = solve(&foil, alpha).unwrap();
        for (got, want) in state.sources.iter().zip(&existing.sources) {
            assert_eq!(got.to_bits(), want.to_bits());
        }
        assert_eq!(state.sheet_density.to_bits(), existing.gamma.to_bits());
        for (&got, &want) in state.tangential_velocity.iter().zip(&existing.vt) {
            assert_close(got, want, 3e-12);
        }
        assert!(state.max_normal_residual < 2e-12);
        assert!(state.kutta_residual.abs() < 2e-12);
    }
}

#[test]
fn g0_kelvin_newborn_satisfies_normal_flow_and_kutta_at_each_position() {
    let foil = naca4_symmetric(0.12, 64).unwrap();
    let prepared = PreparedPanel2d::new(&foil).unwrap();
    let core = 0.03;
    let previous_bound = 0.14;
    let old_wake = [([1.7, 0.06], -0.12), ([2.5, -0.08], -0.02)];
    let old_incident = |point| {
        let mut velocity = [0.9968, 0.0799];
        for (position, gamma) in old_wake {
            let induced = vortex_velocity(point, position, gamma, core);
            velocity[0] += induced[0];
            velocity[1] += induced[1];
        }
        velocity
    };
    let incident: Vec<_> = prepared
        .midpoints()
        .iter()
        .map(|&p| old_incident(p))
        .collect();
    for newborn_position in [[1.04, 0.0], [1.10, 0.025]] {
        let unit_newborn: Vec<_> = prepared
            .midpoints()
            .iter()
            .map(|&point| vortex_velocity(point, newborn_position, 1.0, core))
            .collect();
        let (state, newborn) = prepared
            .solve_kelvin(&incident, &unit_newborn, previous_bound)
            .unwrap();
        assert_close(state.bound_circulation + newborn, previous_bound, 2e-12);
        let total_circulation = state.bound_circulation
            + newborn
            + old_wake.iter().map(|(_, gamma)| gamma).sum::<f64>();
        assert!(total_circulation.abs() < 2e-12);
        assert!(state.max_normal_residual < 2e-12);
        assert!(state.kutta_residual.abs() < 2e-12);

        // An ordinary RHS solve with the resulting newborn independently
        // checks the linear-combination solution of the augmented system.
        let combined: Vec<_> = incident
            .iter()
            .zip(&unit_newborn)
            .map(|(old, unit)| [old[0] + newborn * unit[0], old[1] + newborn * unit[1]])
            .collect();
        let resolved = prepared.solve_incident(&combined).unwrap();
        assert_close(state.sheet_density, resolved.sheet_density, 2e-12);
        for (&got, &want) in state.sources.iter().zip(&resolved.sources) {
            assert_close(got, want, 2e-12);
        }

        // Check the actual exterior field, not only the retained matrix rows.
        // The small normal offset takes the exterior limit of the source jump.
        let mut exterior_kutta = 0.0;
        let n = foil.nodes().len();
        for (i, &midpoint) in prepared.midpoints().iter().enumerate() {
            let a = foil.nodes()[i];
            let b = foil.nodes()[(i + 1) % n];
            let dx = b[0] - a[0];
            let dy = b[1] - a[1];
            let length = (dx * dx + dy * dy).sqrt();
            let tangent = [dx / length, dy / length];
            let normal = [-tangent[1], tangent[0]];
            let point = [
                midpoint[0] + 1e-9 * normal[0],
                midpoint[1] + 1e-9 * normal[1],
            ];
            let body = prepared.body_velocity(&state, point).unwrap();
            let old = old_incident(point);
            let new = vortex_velocity(point, newborn_position, newborn, core);
            let velocity = [body[0] + old[0] + new[0], body[1] + old[1] + new[1]];
            let normal_velocity = velocity[0] * normal[0] + velocity[1] * normal[1];
            let tangential_velocity = velocity[0] * tangent[0] + velocity[1] * tangent[1];
            assert!(
                normal_velocity.abs() < 5e-6,
                "panel={i}, normal={normal_velocity}"
            );
            assert_close(tangential_velocity, state.tangential_velocity[i], 5e-6);
            if i == 0 || i == n - 1 {
                exterior_kutta += tangential_velocity;
            }
        }
        assert!(exterior_kutta.abs() < 5e-6);
    }
}

#[test]
fn g0_body_field_circulation_is_clockwise_density_times_perimeter() {
    let foil = naca4_symmetric(0.12, 48).unwrap();
    let prepared = PreparedPanel2d::new(&foil).unwrap();
    let incident = vec![[0.994, 0.11]; foil.nodes().len()];
    let state = prepared.solve_incident(&incident).unwrap();
    let mut circulation = 0.0;
    let count = 64;
    let dtheta = std::f64::consts::TAU / f64::from(count);
    let radius = 4.0;
    for i in 0..count {
        let theta = dtheta * f64::from(i);
        let (sin, cos) = theta.sin_cos();
        let point = [0.5 + radius * cos, radius * sin];
        let velocity = prepared.body_velocity(&state, point).unwrap();
        circulation += (velocity[0] * sin - velocity[1] * cos) * radius * dtheta;
    }
    assert!(state.bound_circulation > 0.0);
    assert_close(circulation, state.bound_circulation, 2e-11);
}

#[test]
fn g3_incident_solve_and_exterior_field_respect_rigid_motion_and_scale() {
    let foil = naca4_symmetric(0.12, 48).unwrap();
    let prepared = PreparedPanel2d::new(&foil).unwrap();
    let incident: Vec<_> = prepared
        .midpoints()
        .iter()
        .map(|&p| {
            let wake = vortex_velocity(p, [1.3, 0.1], -0.15, 0.04);
            [0.97 + wake[0], 0.11 + wake[1]]
        })
        .collect();
    let original = prepared.solve_incident(&incident).unwrap();
    let scale = 3.5;
    let (sin, cos) = 0.37f64.sin_cos();
    let rotate = |p: [f64; 2]| [cos * p[0] - sin * p[1], sin * p[0] + cos * p[1]];
    let transform = |p| {
        let rotated = rotate(p);
        [2.0 + scale * rotated[0], -1.0 + scale * rotated[1]]
    };
    let transformed = Airfoil2d::new(foil.nodes().iter().map(|&p| transform(p)).collect()).unwrap();
    let moved = PreparedPanel2d::new(&transformed).unwrap();
    let moved_incident: Vec<_> = incident.iter().map(|&velocity| rotate(velocity)).collect();
    let state = moved.solve_incident(&moved_incident).unwrap();
    assert_close(state.sheet_density, original.sheet_density, 3e-11);
    assert_close(
        state.bound_circulation,
        scale * original.bound_circulation,
        3e-11,
    );
    for (&got, &want) in state.sources.iter().zip(&original.sources) {
        assert_close(got, want, 3e-11);
    }
    let point = [1.4, 0.3];
    let expected = rotate(prepared.body_velocity(&original, point).unwrap());
    let actual = moved.body_velocity(&state, transform(point)).unwrap();
    assert_close(actual[0], expected[0], 3e-11);
    assert_close(actual[1], expected[1], 3e-11);
}

#[test]
fn invalid_incidence_and_singular_kelvin_coupling_are_refused() {
    let foil = naca4_symmetric(0.12, 32).unwrap();
    let prepared = PreparedPanel2d::new(&foil).unwrap();
    let incident = vec![[0.97, 0.15]; foil.nodes().len()];
    assert!(matches!(
        prepared.solve_incident(&incident[1..]),
        Err(BemError::VectorLength { .. })
    ));
    let mut bad_incident = incident.clone();
    bad_incident[2][1] = f64::NAN;
    assert!(matches!(
        prepared.solve_incident(&bad_incident),
        Err(BemError::InvalidScalar { .. })
    ));
    let state = prepared.solve_incident(&incident).unwrap();
    // Linearity constructs a response Γ₁ = -1, the singular Kelvin pivot.
    let unit_response: Vec<_> = incident
        .iter()
        .map(|u| {
            [
                -u[0] / state.bound_circulation,
                -u[1] / state.bound_circulation,
            ]
        })
        .collect();
    assert!(matches!(
        prepared.solve_kelvin(&incident, &unit_response, 0.0),
        Err(BemError::LinearSolve {
            stage: "2D panel Kelvin coupling",
            ..
        })
    ));
    assert!(prepared.body_velocity(&state, [f64::NAN, 0.0]).is_err());
}
