use super::*;

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}"
    );
}

fn triangle() -> [Vector; 3] {
    [[0.0; 3], [2.0, 0.0, 0.0], [0.0, 3.0, 0.0]]
}

fn check_balance(q: [Vector; 3], v: [Vector; 3], response: TriangleSurfaceViscosityResponse) {
    let mut resultant = [0.0; 3];
    let mut torque = [0.0; 3];
    let mut power = 0.0;
    for i in 0..3 {
        let moment = cross(q[i], response.forces_n[i]);
        power += dot(v[i], response.forces_n[i]);
        for d in 0..3 {
            resultant[d] += response.forces_n[i][d];
            torque[d] += moment[d];
        }
    }
    for d in 0..3 {
        close(resultant[d], 0.0, 2e-13);
        close(torque[d], 0.0, 2e-13);
    }
    assert!(response.dissipation_w >= 0.0);
    close(power, -response.dissipation_w, 2e-13);
}

#[test]
fn analytical_shear_and_dilation_separate_the_two_surface_viscosities() {
    let q = triangle();
    let viscosity = SurfaceViscosity {
        shear_n_s_m: 0.7,
        dilatational_n_s_m: 1.1,
    };
    let shear_rate = 0.4;
    let v = q.map(|p| [shear_rate * p[1], 0.0, 0.0]);
    let shear = triangle_surface_viscosity(q, v, viscosity).unwrap();
    close(
        shear.dissipation_w,
        3.0 * viscosity.shear_n_s_m * shear_rate * shear_rate,
        2e-15,
    );
    close(
        shear.forces_n[1][1],
        -1.5 * viscosity.shear_n_s_m * shear_rate,
        2e-15,
    );
    close(
        shear.forces_n[2][0],
        -viscosity.shear_n_s_m * shear_rate,
        2e-15,
    );
    check_balance(q, v, shear);
    let expansion_rate = -0.3;
    let v = q.map(|p| p.map(|x| expansion_rate * x));
    let dilation = triangle_surface_viscosity(q, v, viscosity).unwrap();
    close(
        dilation.dissipation_w,
        4.0 * 3.0 * viscosity.dilatational_n_s_m * expansion_rate * expansion_rate,
        3e-15,
    );
    close(
        dilation.forces_n[1][0],
        -3.0 * viscosity.dilatational_n_s_m * expansion_rate,
        3e-15,
    );
    check_balance(q, v, dilation);
}

#[test]
fn rigid_motion_is_not_damped_and_superposed_spin_does_not_change_the_stress() {
    let q = [[0.2, -0.1, 0.7], [1.3, 0.4, -0.2], [-0.3, 1.6, 0.8]];
    let omega = [0.4, -0.8, 1.3];
    let rigid = q.map(|p| {
        let spin = cross(omega, p);
        [spin[0] + 0.7, spin[1] - 0.2, spin[2] + 0.1]
    });
    let viscosity = SurfaceViscosity {
        shear_n_s_m: 0.8,
        dilatational_n_s_m: 0.6,
    };
    let response = triangle_surface_viscosity(q, rigid, viscosity).unwrap();
    close(response.dissipation_w, 0.0, 1e-29);
    for force in response.forces_n.into_iter().flatten() {
        close(force, 0.0, 3e-15);
    }
    let v = [[0.2, -0.4, 0.8], [0.1, 0.3, -0.2], [0.7, -0.5, 0.4]];
    let original = triangle_surface_viscosity(q, v, viscosity).unwrap();
    let superposed = core::array::from_fn(|i| core::array::from_fn(|d| v[i][d] + rigid[i][d]));
    let moved = triangle_surface_viscosity(q, superposed, viscosity).unwrap();
    close(original.dissipation_w, moved.dissipation_w, 3e-15);
    for i in 0..3 {
        for d in 0..3 {
            close(original.forces_n[i][d], moved.forces_n[i][d], 3e-15);
        }
    }
    check_balance(q, v, original);
}

#[test]
fn changing_geometry_and_velocity_tangent_matches_the_same_objective_force() {
    let q = [[0.2, -0.1, 0.7], [1.3, 0.4, -0.2], [-0.3, 1.6, 0.8]];
    let v = [[0.2, -0.4, 0.8], [0.1, 0.3, -0.2], [0.7, -0.5, 0.4]];
    let dq = [[0.1, 0.4, -0.2], [-0.3, 0.1, 0.2], [0.4, -0.1, 0.3]];
    let dv = [[0.2, 0.1, -0.4], [-0.1, 0.3, 0.2], [0.2, 0.5, -0.1]];
    let viscosity = SurfaceViscosity {
        shear_n_s_m: 0.7,
        dilatational_n_s_m: 1.1,
    };
    for (dx, du) in [(dq, [[0.0; 3]; 3]), ([[0.0; 3]; 3], dv), (dq, dv)] {
        let exact = triangle_surface_viscosity_tangent(q, v, viscosity, dx, du).unwrap();
        let epsilon = 2e-6;
        let shift = |base: [Vector; 3], direction: [Vector; 3], sign: f64| {
            core::array::from_fn(|i| {
                core::array::from_fn(|d| base[i][d] + sign * epsilon * direction[i][d])
            })
        };
        let plus =
            triangle_surface_viscosity(shift(q, dx, 1.0), shift(v, du, 1.0), viscosity).unwrap();
        let minus =
            triangle_surface_viscosity(shift(q, dx, -1.0), shift(v, du, -1.0), viscosity).unwrap();
        for i in 0..3 {
            for d in 0..3 {
                close(
                    exact[i][d],
                    (plus.forces_n[i][d] - minus.forces_n[i][d]) / (2.0 * epsilon),
                    2e-9,
                );
            }
        }
    }
}

#[test]
fn patch_assembly_preserves_power_and_winding_independence_with_nonuniform_coefficients() {
    let q = [[0.0; 3], [2.0, 0.0, 0.0], [0.0, 3.0, 0.0], [2.0, 3.0, 0.0]];
    let v = q.map(|p| [0.2 * p[0] + 0.4 * p[1], -0.1 * p[1], 0.3 * p[0]]);
    let coefficients = [
        SurfaceViscosity {
            shear_n_s_m: 0.7,
            dilatational_n_s_m: 1.1,
        },
        SurfaceViscosity {
            shear_n_s_m: 1.4,
            dilatational_n_s_m: 0.3,
        },
    ];
    let response = assemble_surface_viscosity(
        &q,
        &v,
        &[[0, 1, 2], [1, 3, 2]],
        SurfaceViscosityField::PerFace(&coefficients),
    )
    .unwrap();
    let reverse = assemble_surface_viscosity(
        &q,
        &v,
        &[[0, 2, 1], [1, 2, 3]],
        SurfaceViscosityField::PerFace(&coefficients),
    )
    .unwrap();
    close(response.dissipation_w, reverse.dissipation_w, 3e-15);
    let power: f64 = v
        .iter()
        .zip(&response.forces_n)
        .map(|(&v, &f)| dot(v, f))
        .sum();
    close(power, -response.dissipation_w, 3e-15);
    for (a, b) in response
        .forces_n
        .iter()
        .flatten()
        .zip(reverse.forces_n.iter().flatten())
    {
        close(*a, *b, 3e-15);
    }
}

#[test]
fn malformed_surface_viscosities_and_nonfinite_kinematics_are_refused() {
    let q = triangle();
    for coefficient in [-1.0, f64::INFINITY, f64::NAN] {
        let bad = SurfaceViscosity {
            shear_n_s_m: coefficient,
            dilatational_n_s_m: 0.0,
        };
        assert_eq!(
            triangle_surface_viscosity(q, q, bad),
            Err(CapillarySurfaceError::InvalidSurfaceViscosity)
        );
    }
    let mut v = q;
    v[1][0] = f64::NAN;
    assert_eq!(
        triangle_surface_viscosity(q, v, SurfaceViscosity::ZERO),
        Err(CapillarySurfaceError::NonFiniteVelocity)
    );
    assert_eq!(
        assemble_surface_viscosity(&q, &q, &[[0, 1, 2]], SurfaceViscosityField::PerFace(&[])),
        Err(CapillarySurfaceError::LengthMismatch)
    );
}
