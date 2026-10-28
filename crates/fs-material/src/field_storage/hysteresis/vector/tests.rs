use super::*;

fn branch(chi: f64, hc: f64, zs: f64) -> HysteresisBranch {
    HysteresisBranch {
        susceptibility_h_m: chi,
        pinning_field_a_m: hc,
        saturation_polarization_t: zs,
    }
}
fn controls() -> VectorHysteresisControls {
    VectorHysteresisControls {
        relative_tolerance: 1e-12,
        absolute_field_tolerance_a_m: 1e-12,
        absolute_flux_tolerance_t: 1e-12,
        work_tolerance_j_m3: 2e-11,
        maximum_branch_iterations: 100,
        maximum_inverse_iterations: 60,
        maximum_backtracks: 40,
    }
}
fn law() -> VectorHysteresis {
    VectorHysteresis::new(0.01, 30.0, vec![branch(0.2, 2.0, 1.0)]).unwrap()
}
fn near(a: f64, b: f64, tol: f64) {
    assert!((a - b).abs() <= tol, "{a:.16e} != {b:.16e}; tol {tol:e}");
}
fn vector_near(a: Vector, b: Vector, tol: f64) {
    for i in 0..3 {
        near(a[i], b[i], tol);
    }
}
fn step(law: &VectorHysteresis, state: &VectorHysteresisState, h: Vector) -> VectorHysteresisStep {
    law.candidate_field(state, h, controls(), &mut || false)
        .unwrap()
}

#[test]
fn collinear_history_reduces_to_original_scalar_states_and_physical_heat() {
    let branches = vec![
        branch(0.2, 2.0, 1.0),
        branch(0.07, 4.0, 0.4),
        branch(0.1, 0.0, 0.6),
    ];
    let scalar = ScalarHysteresis::new(0.01, 30.0, branches.clone()).unwrap();
    let vector = VectorHysteresis::new(0.01, 30.0, branches).unwrap();
    let mut a = scalar.initial_state();
    let mut b = vector.initial_state();
    for h in [
        0.0, 1.0, 3.0, 12.0, 2.0, -3.0, -12.0, 4.0, 0.0, 8.0, 1.0, 8.0,
    ] {
        let x = scalar.candidate_field(&a, h).unwrap();
        let y = step(&vector, &b, [h, 0.0, 0.0]);
        near(y.state.flux_density_t()[0], x.state.flux_density_t(), 2e-12);
        near(
            y.state.stored_energy_j_m3(),
            x.state.stored_energy_j_m3(),
            2e-12,
        );
        near(y.irreversible_heat_j_m3, x.irreversible_heat_j_m3, 2e-12);
        for (z, &s) in y
            .state
            .polarizations_t()
            .iter()
            .zip(x.state.polarizations_t())
        {
            vector_near(*z, [s, 0.0, 0.0], 2e-12);
        }
        near(
            y.energy_change_j_m3 + y.irreversible_heat_j_m3,
            x.magnetic_work_j_m3,
            3e-12,
        );
        assert!(y.numerical_dissipation_j_m3 >= 0.0);
        near(
            y.endpoint_work_j_m3,
            y.energy_change_j_m3 + y.irreversible_heat_j_m3 + y.numerical_dissipation_j_m3,
            2e-11,
        );
        a = x.state;
        b = y.state;
    }
}

#[test]
fn saturated_rotating_return_map_matches_independent_spherical_minimization() {
    let material = law();
    let old = step(&material, &material.initial_state(), [5.0, 0.0, 0.0]).state;
    vector_near(old.polarizations_t()[0], [0.6, 0.0, 0.0], 1e-14);
    let h = [3.0, 8.0, 0.0];
    let result = step(&material, &old, h);
    let z = result.state.polarizations_t()[0];
    near(norm(z), 1.0, 1e-14);
    // Isotropy keeps the minimizer in the plane spanned by old z and H.
    // Independently minimize the original objective on that unit circle,
    // first by exhaustive angular bracketing, then golden-section refinement.
    let objective = |theta: f64| {
        let z = [theta.cos(), theta.sin(), 0.0];
        2.5 - dot(h, z) + 2.0 * norm(sub(z, old.polarizations_t()[0]))
    };
    let count = 8192;
    let width = 2.0 * core::f64::consts::PI / count as f64;
    let best = (0..count)
        .min_by(|&i, &j| objective(i as f64 * width).total_cmp(&objective(j as f64 * width)))
        .unwrap();
    let (mut lo, mut hi) = ((best as f64 - 1.0) * width, (best as f64 + 1.0) * width);
    let ratio = (5.0_f64.sqrt() - 1.0) / 2.0;
    for _ in 0..80 {
        let a = hi - ratio * (hi - lo);
        let b = lo + ratio * (hi - lo);
        if objective(a) < objective(b) {
            hi = b;
        } else {
            lo = a;
        }
    }
    let theta = f64::midpoint(lo, hi);
    vector_near(z, [theta.cos(), theta.sin(), 0.0], 3e-8);
    // The commonly tempting unconstrained-update-then-project rule is wrong.
    let projected = [0.6 / 1.8_f64.sqrt(), 1.2 / 1.8_f64.sqrt(), 0.0];
    assert!(norm(sub(z, projected)) > 0.01);
    assert!(result.maximum_constitutive_residual_a_m < 2e-11);
}

fn rotation() -> Matrix {
    let a = [
        1.0 / 14.0_f64.sqrt(),
        2.0 / 14.0_f64.sqrt(),
        3.0 / 14.0_f64.sqrt(),
    ];
    let (s, c) = 0.73_f64.sin_cos();
    let skew = [[0.0, -a[2], a[1]], [a[2], 0.0, -a[0]], [-a[1], a[0], 0.0]];
    std::array::from_fn(|i| {
        std::array::from_fn(|j| {
            (if i == j { c } else { 0.0 }) + (1.0 - c) * a[i] * a[j] + s * skew[i][j]
        })
    })
}

#[test]
fn whole_history_is_rotation_equivariant_including_heat_and_tensor_tangent() {
    let material = VectorHysteresis::new(
        0.01,
        30.0,
        vec![branch(0.2, 2.0, 1.0), branch(0.1, 0.0, 0.5)],
    )
    .unwrap();
    let r = rotation();
    let mut original = material.initial_state();
    let mut rotated = material.initial_state();
    for h in [
        [4.0, 0.0, 0.0],
        [3.0, 5.0, 1.0],
        [-2.0, 7.0, 3.0],
        [-8.0, -3.0, 2.0],
        [0.0, 0.0, 0.0],
        [2.0, -4.0, 6.0],
    ] {
        let a = step(&material, &original, h);
        let b = step(&material, &rotated, mv(r, h));
        vector_near(
            b.state.flux_density_t(),
            mv(r, a.state.flux_density_t()),
            2e-10,
        );
        for (&x, &y) in a
            .state
            .polarizations_t()
            .iter()
            .zip(b.state.polarizations_t())
        {
            vector_near(y, mv(r, x), 2e-10);
        }
        near(a.irreversible_heat_j_m3, b.irreversible_heat_j_m3, 2e-10);
        near(
            a.numerical_dissipation_j_m3,
            b.numerical_dissipation_j_m3,
            2e-10,
        );
        near(
            a.state.stored_energy_j_m3(),
            b.state.stored_energy_j_m3(),
            2e-10,
        );
        for i in 0..3 {
            for j in 0..3 {
                let expected = (0..3)
                    .flat_map(|k| {
                        (0..3)
                            .map(move |l| r[i][k] * a.incremental_permeability_h_m[k][l] * r[j][l])
                    })
                    .sum();
                near(b.incremental_permeability_h_m[i][j], expected, 2e-10);
            }
        }
        original = a.state;
        rotated = b.state;
    }
}

#[test]
fn analytic_rotating_loop_has_nonzero_loss_and_converges_under_history_refinement() {
    let material = VectorHysteresis::new(0.01, 30.0, vec![branch(0.2, 2.0, 10.0)]).unwrap();
    let amplitude: f64 = 5.0;
    let pinning: f64 = 2.0;
    let chi = 0.2;
    // Independent steady continuous solution: H-z/chi is tangent to the
    // polarization circle. Therefore radius=chi*sqrt(Hamp²-hc²),
    // phase lag=asin(hc/Hamp), and Q/cycle=2*pi*hc*radius.
    let radius = chi * (amplitude * amplitude - pinning * pinning).sqrt();
    let heat_exact = 2.0 * core::f64::consts::PI * pinning * radius;
    let mut errors = vec![];
    let mut numerical = vec![];
    for n in [64, 128, 256] {
        let mut state = material.initial_state();
        let (mut heat, mut num, mut work) = (0.0, 0.0, 0.0);
        let mut first_energy = 0.0;
        for i in 0..=8 * n {
            let angle = 2.0 * core::f64::consts::PI * i as f64 / n as f64;
            let result = step(
                &material,
                &state,
                [amplitude * angle.cos(), amplitude * angle.sin(), 0.0],
            );
            if i == 7 * n {
                first_energy = result.state.stored_energy_j_m3();
            }
            if i > 7 * n {
                heat += result.irreversible_heat_j_m3;
                num += result.numerical_dissipation_j_m3;
                work += result.endpoint_work_j_m3;
            }
            state = result.state;
        }
        let z = state.polarizations_t()[0];
        assert!(z[1] < -0.1); // persistent phase lag, not independent scalar play.
        let delta = 2.0 * core::f64::consts::PI / n as f64;
        let discrete_lag = (pinning * (delta / 2.0).cos() / amplitude).asin();
        let discrete_radius = chi
            * ((amplitude * amplitude - pinning * pinning * (delta / 2.0).cos().powi(2)).sqrt()
                - pinning * (delta / 2.0).sin());
        vector_near(
            z,
            [
                discrete_radius * discrete_lag.cos(),
                -discrete_radius * discrete_lag.sin(),
                0.0,
            ],
            2e-8,
        );
        near(
            work,
            heat + num + state.stored_energy_j_m3() - first_energy,
            2e-10,
        );
        assert!(heat > 10.0 && num > 0.0);
        errors.push((heat - heat_exact).abs());
        numerical.push(num);
    }
    for pair in errors.windows(2) {
        assert!(
            (0.8..1.2).contains(&(pair[0] / pair[1]).log2()),
            "{errors:?}"
        );
    }
    for pair in numerical.windows(2) {
        assert!((1.8..2.2).contains(&(pair[0] / pair[1])), "{numerical:?}");
    }
}

#[test]
fn consistent_active_set_tangent_matches_finite_differences_and_retains_spd_background() {
    let material = law();
    let old = step(&material, &material.initial_state(), [5.0, 0.0, 0.0]).state;
    for h in [
        [3.0, 0.0, 0.0],
        [4.0, 2.0, 1.0],
        [3.0, 8.0, 0.0],
        [10.0, 0.0, 0.0],
        [-4.0, -3.0, -2.0],
    ] {
        let result = step(&material, &old, h);
        let jac = result.incremental_permeability_h_m;
        for j in 0..3 {
            let mut hp = h;
            let mut hm = h;
            hp[j] += 1e-5;
            hm[j] -= 1e-5;
            let plus = step(&material, &old, hp).state.flux_density_t();
            let minus = step(&material, &old, hm).state.flux_density_t();
            for i in 0..3 {
                near((plus[i] - minus[i]) / 2e-5, jac[i][j], 2e-7);
                near(jac[i][j], jac[j][i], 2e-12);
            }
        }
        for v in [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [1.0, -2.0, 3.0]] {
            assert!(dot(v, mv(jac, v)) >= 0.01 * dot(v, v) - 1e-12);
        }
    }
}

#[test]
fn flux_driven_newton_recovers_the_actual_vector_history_and_transfers() {
    let material = VectorHysteresis::new(
        0.02,
        30.0,
        vec![branch(0.2, 2.0, 1.0), branch(0.07, 3.0, 0.4)],
    )
    .unwrap();
    let mut state = material.initial_state();
    for h in [
        [4.0, 0.0, 0.0],
        [3.0, 5.0, 1.0],
        [-2.0, 7.0, 3.0],
        [-8.0, -3.0, 2.0],
        [0.0, 0.0, 0.0],
        [2.0, -4.0, 6.0],
    ] {
        let direct = step(&material, &state, h);
        let inverse = material
            .candidate_flux_density(
                &state,
                direct.state.flux_density_t(),
                controls(),
                &mut || false,
            )
            .unwrap();
        vector_near(inverse.state.field_a_m(), h, 3e-9);
        vector_near(
            inverse.state.flux_density_t(),
            direct.state.flux_density_t(),
            3e-11,
        );
        for (&a, &b) in direct
            .state
            .polarizations_t()
            .iter()
            .zip(inverse.state.polarizations_t())
        {
            vector_near(a, b, 3e-10);
        }
        near(
            inverse.irreversible_heat_j_m3,
            direct.irreversible_heat_j_m3,
            3e-9,
        );
        near(inverse.energy_change_j_m3, direct.energy_change_j_m3, 3e-9);
        assert!(inverse.inverse_residual_t < 3e-11);
        state = direct.state;
    }
}

#[test]
fn reversible_and_empty_branch_limits_separate_endpoint_remainder_from_heat() {
    for branches in [vec![], vec![branch(0.2, 0.0, 1.0)]] {
        let reversible_slope = 0.01 + branches.iter().map(|b| b.susceptibility_h_m).sum::<f64>();
        let material = VectorHysteresis::new(0.01, 30.0, branches).unwrap();
        let initial = material.initial_state();
        let zero = step(&material, &initial, [0.0; 3]);
        for i in 0..3 {
            near(
                zero.incremental_permeability_h_m[i][i],
                reversible_slope,
                1e-15,
            );
        }
        let a = step(&material, &initial, [7.0, 8.0, 3.0]);
        let b = step(&material, &a.state, [0.0; 3]);
        assert_eq!(a.irreversible_heat_j_m3, 0.0);
        assert_eq!(b.irreversible_heat_j_m3, 0.0);
        assert!(a.numerical_dissipation_j_m3 > 0.0 && b.numerical_dissipation_j_m3 > 0.0);
        near(a.energy_change_j_m3 + b.energy_change_j_m3, 0.0, 1e-12);
        vector_near(b.state.flux_density_t(), [0.0; 3], 1e-13);
        let hold = step(&material, &a.state, a.state.field_a_m());
        near(hold.irreversible_heat_j_m3, 0.0, 1e-13);
        near(hold.energy_change_j_m3, 0.0, 1e-13);
    }
}

#[test]
fn original_parameters_domains_work_budgets_and_late_cancellation_are_transactional() {
    let material = law();
    let old = step(&material, &material.initial_state(), [5.0, 0.0, 0.0]).state;
    let before = old.clone();
    let field = [3.0, 8.0, 0.0];
    let other = VectorHysteresis::new(0.02, 30.0, vec![branch(0.2, 2.0, 1.0)]).unwrap();
    assert!(
        other
            .candidate_field(&old, field, controls(), &mut || false)
            .is_err()
    );
    assert!(
        material
            .candidate_field(&old, [31.0, 0.0, 0.0], controls(), &mut || false)
            .is_err()
    );
    assert!(
        material
            .candidate_field(&old, [f64::NAN, 0.0, 0.0], controls(), &mut || false)
            .is_err()
    );
    assert!(
        material
            .candidate_flux_density(&old, [f64::MAX; 3], controls(), &mut || false)
            .is_err()
    );
    assert!(
        material
            .candidate_flux_density(&old, [100.0, 0.0, 0.0], controls(), &mut || false)
            .is_err()
    );
    let mut exhausted = controls();
    exhausted.maximum_branch_iterations = 0;
    assert!(matches!(
        material.candidate_field(&old, field, exhausted, &mut || false),
        Err(VectorHysteresisError::BranchNonconvergence { .. })
    ));
    let actual = step(&material, &old, field);
    exhausted = controls();
    exhausted.maximum_inverse_iterations = 0;
    assert!(matches!(
        material.candidate_flux_density(
            &old,
            actual.state.flux_density_t(),
            exhausted,
            &mut || false
        ),
        Err(VectorHysteresisError::FluxNonconvergence { .. })
    ));
    let mut calls = 0;
    material
        .candidate_flux_density(&old, actual.state.flux_density_t(), controls(), &mut || {
            calls += 1;
            false
        })
        .unwrap();
    assert!(calls > 10);
    let mut count = 0;
    assert!(matches!(
        material.candidate_flux_density(
            &old,
            actual.state.flux_density_t(),
            controls(),
            &mut || {
                count += 1;
                count == calls - 1
            }
        ),
        Err(VectorHysteresisError::Cancelled)
    ));
    assert_eq!(old, before);
    assert_eq!(step(&material, &old, field), actual);
}
