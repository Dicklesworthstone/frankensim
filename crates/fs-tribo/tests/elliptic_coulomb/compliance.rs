use super::*;
use fs_tribo::tangential::elliptic::SymmetricTangentialCompliance;

fn response(matrix: [[f64; 2]; 2]) -> SymmetricTangentialCompliance {
    SymmetricTangentialCompliance::new(matrix).unwrap()
}

fn rotated_matrix(c: [[f64; 2]; 2], angle: f64) -> [[f64; 2]; 2] {
    let u = rotate([1.0, 0.0], angle);
    let v = rotate([0.0, 1.0], angle);
    let entry = |i, j| {
        c[0][0] * u[i] * u[j] + c[0][1] * (u[i] * v[j] + v[i] * u[j]) + c[1][1] * v[i] * v[j]
    };
    let xy = entry(0, 1);
    [[entry(0, 0), xy], [xy, entry(1, 1)]]
}

#[test]
fn direct_full_response_matches_independent_ellipse_energy_minimization() {
    let l = law();
    for c in [
        [[0.008, 0.003], [0.003, 0.016]],
        [[0.012, -0.006], [-0.006, 0.02]],
        [[0.02, 0.009], [0.009, 0.006]],
    ] {
        for free in [[0.15, 0.04], [-0.13, 0.07], [0.09, -0.12]] {
            let x = input(free, 0.0);
            let step = l
                .advance_with_compliance(&l.initial_state(), &x, response(c))
                .unwrap();
            assert_eq!(step.regime, FrictionRegime::Sliding);
            close2(step.traction_n, angular_projection(free, c), 4e-7);
            for i in 0..2 {
                close(
                    step.slip_increment_m[i]
                        + c[i][0] * step.traction_n[0]
                        + c[i][1] * step.traction_n[1],
                    free[i],
                    1e-15,
                );
            }
            let u = step.slip_increment_m;
            close(
                step.dissipated_j,
                (36.0 * u[0] * u[0] + 4.0 * u[1] * u[1]).sqrt(),
                2e-11,
            );
            assert!(step.iterations <= x.solve.max_iterations);
            assert!(
                step.compatibility_residual_m
                    .iter()
                    .all(|r| r.abs() < 1e-15)
            );
            let wrapped = advance_with_matrix_compliance(
                &l,
                &l.initial_state(),
                &x,
                c,
                matrix_controls(),
                || false,
            )
            .unwrap();
            close2(step.traction_n, wrapped.contact.traction_n, 5e-8);
        }
    }
}

#[test]
fn matrix_static_kinetic_branches_and_scalar_fast_paths_retain_original_owners() {
    let l = law();
    let c = response([[0.01, 0.004], [0.004, 0.02]]);
    let x = input([0.07, 0.028], 0.0);
    let resting = l
        .advance_with_compliance(&l.initial_state(), &x, c)
        .unwrap();
    assert_eq!(resting.regime, FrictionRegime::Sticking);
    close2(resting.traction_n, [7.0, 0.0], 2e-14);
    assert_eq!(resting.slip_increment_m, [0.0; 2]);
    let moving = l
        .advance_kinetic_with_compliance(&l.initial_state(), &x, c)
        .unwrap();
    assert_eq!(moving.regime, FrictionRegime::Sliding);
    assert!(moving.dissipated_j > 0.0);
    let stopped = l.advance_with_compliance(&moving.state, &x, c).unwrap();
    assert_eq!(stopped.regime, FrictionRegime::Sticking);
    assert_eq!(stopped.state.slip_m(), moving.state.slip_m());
    for angle in [0.0, 0.73] {
        let l = EllipticCoulomb::new([0.8, 0.4], [0.6, 0.2], rotate([1.0, 0.0], angle)).unwrap();
        for scalar in [0.0, 0.005, 0.02] {
            for free in [[0.0, 0.0], [0.02, 0.03], [0.15, -0.04]] {
                let scalar_input = input(free, scalar);
                let matrix_input = input(free, 0.0);
                let c = response([[scalar, 0.0], [0.0, scalar]]);
                assert_eq!(
                    l.advance(&l.initial_state(), &scalar_input).unwrap(),
                    l.advance_with_compliance(&l.initial_state(), &matrix_input, c)
                        .unwrap()
                );
                assert_eq!(
                    l.advance_kinetic(&l.initial_state(), &scalar_input)
                        .unwrap(),
                    l.advance_kinetic_with_compliance(&l.initial_state(), &matrix_input, c)
                        .unwrap()
                );
            }
        }
    }
}

#[test]
fn full_response_rotates_jointly_with_texture_and_body_coordinates() {
    let l = law();
    let c = [[0.008, 0.003], [0.003, 0.016]];
    let free = [0.16, 0.09];
    let original = l
        .advance_with_compliance(&l.initial_state(), &input(free, 0.0), response(c))
        .unwrap();
    for angle in [0.073, 0.37, 1.12, -0.82] {
        let rotated =
            EllipticCoulomb::new([0.8, 0.4], [0.6, 0.2], rotate([1.0, 0.0], angle)).unwrap();
        let step = rotated
            .advance_with_compliance(
                &rotated.initial_state(),
                &input(rotate(free, angle), 0.0),
                response(rotated_matrix(c, angle)),
            )
            .unwrap();
        close2(rotate(step.traction_n, -angle), original.traction_n, 2e-11);
        close2(
            rotate(step.slip_increment_m, -angle),
            original.slip_increment_m,
            2e-13,
        );
        close(step.dissipated_j, original.dissipated_j, 2e-12);
        assert!(
            step.compatibility_residual_m
                .iter()
                .all(|r| r.abs() < 2e-16)
        );
    }
}

#[test]
fn one_axis_friction_preserves_cross_coupled_free_motion_and_exact_zero_heat() {
    let c = [[0.02, 0.006], [0.006, 0.01]];
    for angle in [0.0, 0.292, 0.803] {
        let l = EllipticCoulomb::new([0.0, 0.4], [0.0, 0.2], rotate([1.0, 0.0], angle)).unwrap();
        for (free, force, slip, heat) in [
            ([0.1, 0.01], [0.0, 1.0], [0.094, 0.0], 0.0),
            ([0.1, 0.05], [0.0, 2.0], [0.088, 0.03], 0.06),
        ] {
            let step = l
                .advance_with_compliance(
                    &l.initial_state(),
                    &input(rotate(free, angle), 0.0),
                    response(rotated_matrix(c, angle)),
                )
                .unwrap();
            close2(rotate(step.traction_n, -angle), force, 2e-14);
            close2(rotate(step.slip_increment_m, -angle), slip, 2e-16);
            close(step.dissipated_j, heat, 2e-16);
            if heat == 0.0 {
                assert_eq!(step.dissipated_j, 0.0);
            }
        }
    }
    let l = EllipticCoulomb::new([0.0, 0.4], [0.0, 0.2], [1.0, 0.0]).unwrap();
    let singular = response([[0.0, 0.0], [0.0, 0.01]]);
    let static_step = l
        .advance_with_compliance(&l.initial_state(), &input([0.0, 0.03], 0.0), singular)
        .unwrap();
    assert_eq!(static_step.regime, FrictionRegime::Sticking);
    close2(static_step.traction_n, [0.0, 3.0], 1e-14);
    let moving = l
        .advance_kinetic_with_compliance(&l.initial_state(), &input([0.1, 0.03], 0.0), singular)
        .unwrap();
    close2(moving.traction_n, [0.0, 2.0], 1e-14);
    close2(moving.slip_increment_m, [0.1, 0.01], 1e-15);
}

#[test]
fn singular_sliding_is_resolved_but_rank_deficient_static_and_invalid_inputs_refuse() {
    let l = law();
    let state = l.initial_state();
    let saved = state.clone();
    let c = [[0.01, 0.0], [0.0, 0.0]];
    let free = [0.02, 0.03];
    let step = l
        .advance_with_compliance(&state, &input(free, 0.0), response(c))
        .unwrap();
    close2(step.traction_n, angular_projection(free, c), 3e-7);
    assert!(
        l.advance_with_compliance(&state, &input([0.01, 0.0], 0.0), response(c))
            .is_err()
    );
    assert!(
        l.advance_with_compliance(&state, &input(free, 0.01), response(c))
            .is_err()
    );
    for bad in [
        [[1.0, 0.1], [0.2, 1.0]],
        [[1.0, 2.0], [2.0, 1.0]],
        [[-1.0, 0.0], [0.0, 1.0]],
        [[1.0, f64::NAN], [f64::NAN, 1.0]],
        [[1e308, 1e-308], [1e-308, 0.0]],
        [[1e308, 1e-308], [1e-308, 1e100]],
    ] {
        assert!(SymmetricTangentialCompliance::new(bad).is_err());
    }
    let mut x = input([0.15, 0.04], 0.0);
    x.solve.max_iterations = 1;
    assert!(
        l.advance_with_compliance(&state, &x, response([[0.008, 0.003], [0.003, 0.016]]))
            .is_err()
    );
    let mut x = input([0.07, 0.028], 0.0);
    x.solve.maximum_force_n = 6.5;
    assert!(
        l.advance_with_compliance(&state, &x, response([[0.01, 0.004], [0.004, 0.02]]))
            .is_err()
    );
    let foreign = EllipticCoulomb::new([0.8, 0.4], [0.5, 0.2], [1.0, 0.0]).unwrap();
    assert!(
        foreign
            .advance_with_compliance(&state, &input(free, 0.0), response(c))
            .is_err()
    );
    assert_eq!(state, saved);
    assert_eq!(
        step,
        l.advance_with_compliance(&saved, &input(free, 0.0), response(c))
            .unwrap()
    );
}

#[test]
fn rank_one_tangent_line_has_one_static_reaction_and_no_slip_heat() {
    // The line Fx + Fy = 10 is tangent to (Fx/8)^2 + (Fy/6)^2 = 1.
    // Lagrange multipliers give the unique reaction (6.4, 3.6), without
    // selecting an arbitrary null-space component or adding compliance.
    let l = EllipticCoulomb::new([0.8, 0.6], [0.4, 0.3], [1.0, 0.0]).unwrap();
    let state = l.initial_state();
    for scale in [0.000_976_562_5, 1.0, 1024.0] {
        for sign in [-1.0, 1.0] {
            let c = response([[scale, scale], [scale, scale]]);
            let mut x = input([sign * 10.0 * scale; 2], 0.0);
            let result = l.advance_with_compliance(&state, &x, c).unwrap();
            assert_eq!(result.regime, FrictionRegime::Sticking);
            close2(result.traction_n, [sign * 6.4, sign * 3.6], 2e-15);
            assert_eq!(result.slip_increment_m, [0.0; 2]);
            assert_eq!(result.state.slip_m(), state.slip_m());
            assert_eq!(result.dissipated_j, 0.0);
            assert_eq!(result.closure_residual_j, 0.0);
            assert_eq!(result.iterations, 0);
            assert!(
                result
                    .compatibility_residual_m
                    .iter()
                    .all(|r| r.abs() <= scale * 2e-15)
            );

            // A moving phase retains the weaker kinetic surface and must slide.
            let kinetic = l.advance_kinetic_with_compliance(&state, &x, c).unwrap();
            assert_eq!(kinetic.regime, FrictionRegime::Sliding);
            assert!(kinetic.dissipated_j > 0.0);
            close2(kinetic.traction_n, [sign * 3.2, sign * 1.8], 2e-11);

            // An effort budget cannot change the chosen physical branch.
            x.solve.maximum_force_n = 7.0;
            assert!(l.advance_with_compliance(&state, &x, c).is_err());
        }
    }
    assert_eq!(state, l.initial_state());
}

#[test]
fn rank_one_static_boundary_does_not_capture_neighboring_interior_or_sliding_loads() {
    let l = law();
    let c = response([[1.0, 0.0], [0.0, 0.0]]);
    let state = l.initial_state();
    let boundary = l
        .advance_with_compliance(&state, &input([8.0, 0.0], 0.0), c)
        .unwrap();
    assert_eq!(boundary.regime, FrictionRegime::Sticking);
    assert_eq!(boundary.traction_n, [8.0, 0.0]);
    // Even an interior point one representable number away is nonunique.
    for demand in [4.0, 8.0_f64.next_down()] {
        let mut x = input([demand, 0.0], 0.0);
        x.solve.force_absolute_tolerance_n = 1.0;
        assert!(l.advance_with_compliance(&state, &x, c).is_err());
    }
    let moving = l
        .advance_with_compliance(&state, &input([8.0_f64.next_up(), 0.0], 0.0), c)
        .unwrap();
    assert_eq!(moving.regime, FrictionRegime::Sliding);
    assert!(moving.slip_increment_m[0] > 0.0);
    assert!(moving.dissipated_j > 0.0);
}

#[test]
fn rank_one_boundary_respects_material_axis_and_signed_coordinate_changes() {
    for axis in [[1.0, 0.0], [0.0, 1.0], [-1.0, 0.0], [0.0, -1.0]] {
        let l = EllipticCoulomb::new([0.8, 0.6], [0.4, 0.3], axis).unwrap();
        for sign in [-1.0, 1.0] {
            let c = response([[1.0, sign], [sign, 1.0]]);
            let step = l
                .advance_with_compliance(&l.initial_state(), &input([10.0, sign * 10.0], 0.0), c)
                .unwrap();
            let expected = if axis[0] != 0.0 {
                [6.4, sign * 3.6]
            } else {
                [3.6, sign * 6.4]
            };
            close2(step.traction_n, expected, 2e-15);
            assert_eq!(step.regime, FrictionRegime::Sticking);
            assert_eq!(step.dissipated_j, 0.0);
            assert_eq!(step.state.principal_axis(), axis);
        }
    }
}
