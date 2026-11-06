use super::*;

fn geometry() -> JournalGeometry {
    JournalGeometry { inner_radius_m: 0.1, outer_radius_m: 0.11, length_m: 0.2 }
}
fn solve() -> JournalShearSolve {
    JournalShearSolve { maximum_iterations: 96,
        absolute_angular_rate_tolerance_per_s: 1e-13,
        relative_angular_rate_tolerance: 1e-13 }
}
fn point(rate: f64, mu: f64) -> JournalViscosityKnot {
    JournalViscosityKnot { shear_rate_per_s: rate, dynamic_viscosity_pa_s: mu }
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance * a.abs().max(b.abs()).max(1e-30), "{a:e} != {b:e}");
}

#[test]
fn g1_newtonian_source_recovers_exact_finite_gap_torque_and_zero_speed() {
    let g = geometry();
    let journal = GeneralizedNewtonianJournal::new(g,
        vec![point(0., 0.2), point(100., 0.2)], solve()).unwrap();
    let zero = journal.response(0.).unwrap();
    assert_eq!(zero.shear_rate_per_s, [0.; 2]);
    assert_eq!(zero.dissipated_power_w, 0.);
    for omega in [-2., 2.] {
        let response = journal.response(omega).unwrap();
        let torque = -4. * core::f64::consts::PI * 0.2 * g.length_m * omega
            / (g.inner_radius_m.recip().powi(2) - g.outer_radius_m.recip().powi(2));
        close(response.resisting_torque_nm, torque, 3e-15);
        close(response.dissipated_power_w, -torque * omega, 3e-15);
        close(response.shear_rate_per_s[0], torque.abs()
            / (2. * core::f64::consts::PI * g.length_m * g.inner_radius_m.powi(2) * 0.2), 3e-15);
        assert_eq!(response.iterations, 0);
    }
}

#[test]
fn g1_g3_power_law_torque_and_wall_stress_use_actual_cylinder_radii() {
    // The positive zero-shear plateau is never reached by these solutions.
    let g = geometry();
    for k in [0.01, 0.03] {
        let journal = GeneralizedNewtonianJournal::new(g,
            vec![point(0., 5. * k), point(5., 5. * k), point(100., 100. * k)], solve()).unwrap();
        for omega in [-2., -1., 1., 2.] {
            let response = journal.response(omega).unwrap();
            let inner = omega.abs() / (1. - g.inner_radius_m / g.outer_radius_m);
            let outer = inner * g.inner_radius_m / g.outer_radius_m;
            let torque = -omega.signum() * 2. * core::f64::consts::PI
                * g.length_m * g.inner_radius_m.powi(2) * k * inner.powi(2);
            close(response.shear_rate_per_s[0], inner, 3e-12);
            close(response.shear_rate_per_s[1], outer, 3e-12);
            close(response.resisting_torque_nm, torque, 5e-12);
            close(response.dissipated_power_w, -torque * omega, 5e-12);
            assert!(response.shear_rate_per_s[1] > 5.);
        }
    }
}

#[test]
fn g1_piecewise_shear_thinning_matches_independent_radial_quadrature() {
    let g = JournalGeometry { inner_radius_m: 0.1, outer_radius_m: 0.3, length_m: 0.4 };
    let knots = vec![point(0., 0.2), point(3., 0.2), point(8., 0.14), point(30., 0.14)];
    let journal = GeneralizedNewtonianJournal::new(g, knots.clone(), solve()).unwrap();
    let inner_stress = 2.0;
    // Independent radial midpoint rule and direct μ(gamma)*gamma inversion.
    let mut omega = 0.;
    let n = 32768;
    let dr = (g.outer_radius_m - g.inner_radius_m) / f64::from(n);
    for i in 0..n {
        let radius = g.inner_radius_m + (f64::from(i) + 0.5) * dr;
        let stress = inner_stress * (g.inner_radius_m / radius).powi(2);
        let mut low = 0.;
        let mut high = 30.;
        for _ in 0..60 {
            let rate = 0.5 * (low + high);
            let mu = if rate < 3. { 0.2 } else if rate < 8. {
                0.2 - 0.012 * (rate - 3.)
            } else { 0.14 };
            if rate * mu < stress { low = rate; } else { high = rate; }
        }
        omega += 0.5 * (low + high) * dr / radius;
    }
    let response = journal.response(omega).unwrap();
    let torque = 2. * core::f64::consts::PI * g.length_m
        * g.inner_radius_m.powi(2) * inner_stress;
    close(-response.resisting_torque_nm, torque, 3e-9);
    assert!(response.shear_rate_per_s[0] > 8. && response.shear_rate_per_s[1] < 3.);
}

#[test]
fn g3_almost_constant_source_keeps_the_newtonian_limit() {
    let g = geometry();
    let constant = GeneralizedNewtonianJournal::new(g,
        vec![point(0., 0.2), point(100., 0.2)], solve()).unwrap();
    let slight = GeneralizedNewtonianJournal::new(g,
        vec![point(0., 0.2), point(100., 0.2 * (1. + 1e-12))], solve()).unwrap();
    close(slight.response(2.).unwrap().resisting_torque_nm,
        constant.response(2.).unwrap().resisting_torque_nm, 2e-12);
}

#[test]
fn g1_g3_shear_thinning_retains_positive_wall_rates_near_rest() {
    let g = geometry();
    let controls = JournalShearSolve { absolute_angular_rate_tolerance_per_s: 1e-40,
        relative_angular_rate_tolerance: 1e-13, ..solve() };
    let thinning = GeneralizedNewtonianJournal::new(g,
        vec![point(0., 1.), point(100., 0.75)], controls).unwrap();
    let limit = GeneralizedNewtonianJournal::new(g,
        vec![point(0., 1.), point(100., 1.)], controls).unwrap();
    for omega in [-1e-20, 1e-20] {
        let response = thinning.response(omega).unwrap();
        let expected = limit.response(omega).unwrap();
        // This affine law tends to mu(0)=1; neither positive wall rate may be
        // obtained by subtracting a large inverse increment from gamma=100.
        assert!(response.shear_rate_per_s[0] > response.shear_rate_per_s[1]);
        assert!(response.shear_rate_per_s[1] > 0.0);
        for wall in 0..2 {
            close(response.shear_rate_per_s[wall], expected.shear_rate_per_s[wall], 2e-12);
        }
        close(response.resisting_torque_nm, expected.resisting_torque_nm, 2e-12);
        close(response.dissipated_power_w, expected.dissipated_power_w, 2e-12);
    }
}

#[test]
fn g1_g3_narrow_clearance_preserves_nonlinear_wall_stress_torque_and_heat() {
    let inner_radius_m = 0.1_f64;
    for outer_radius_m in [0.1001, 0.1000001, 0.1000000001,
        f64::from_bits(inner_radius_m.to_bits() + 1)] {
        let g = JournalGeometry { inner_radius_m, outer_radius_m, length_m: 0.2 };
        let controls = JournalShearSolve { absolute_angular_rate_tolerance_per_s: 1e-30,
            relative_angular_rate_tolerance: 2e-14, ..solve() };
        let k = 0.01;
        let inner = 20.0;
        let journal = GeneralizedNewtonianJournal::new(g,
            vec![point(0., 5. * k), point(5., 5. * k), point(100., 100. * k)],
            controls).unwrap();
        // For mu=k*gamma, gamma(r)=gamma_i*ri/r. Keep the actual clearance
        // in this closed form instead of evaluating 1-ri/ro near one.
        let speed = inner * ((outer_radius_m - inner_radius_m) / outer_radius_m);
        for omega in [-speed, speed] {
            let response = journal.response(omega).unwrap();
            let torque = -omega.signum() * 2. * core::f64::consts::PI
                * g.length_m * inner_radius_m.powi(2) * k * inner.powi(2);
            close(response.shear_rate_per_s[0], inner, 2e-12);
            close(response.shear_rate_per_s[1], inner * inner_radius_m / outer_radius_m, 2e-12);
            close(response.resisting_torque_nm, torque, 2e-12);
            close(response.dissipated_power_w, -torque * omega, 2e-12);
            assert!(response.dissipated_power_w > 0.0);
            assert!(response.angular_rate_defect_per_s.abs()
                <= controls.absolute_angular_rate_tolerance_per_s
                    + controls.relative_angular_rate_tolerance * speed);
        }
    }
}

#[test]
fn g1_narrow_gap_crosses_a_constitutive_kink_without_losing_work() {
    // Below gamma=16 the active source is mu=gamma/64; above it mu=1/4.
    // Binary-exact source knots and stress increments isolate radial-gap error.
    for gap in [2.0_f64.powi(-24), 2.0_f64.powi(-34), 2.0_f64.powi(-44)] {
        let g = JournalGeometry { inner_radius_m: 0.125,
            outer_radius_m: 0.125 + gap, length_m: 0.25 };
        let journal = GeneralizedNewtonianJournal::new(g,
            vec![point(0., 0.125), point(8., 0.125), point(16., 0.25), point(64., 0.25)],
            JournalShearSolve { absolute_angular_rate_tolerance_per_s: 1e-30,
                relative_angular_rate_tolerance: 2e-14, ..solve() }).unwrap();
        let inner = 16.0 + 64.0 * gap;
        let radius_ratio = g.inner_radius_m / g.outer_radius_m;
        let delta = (inner - 16.0) / 16.0;
        let root = (1.0 + delta).sqrt();
        // Independent closed form: half the constant-viscosity rate change,
        // plus the power-law rate change. Rationalize 1-q*sqrt(1+delta).
        let speed = 0.5 * (inner - 16.0) + 16.0
            * (gap / g.outer_radius_m - radius_ratio * delta / (1.0 + root));
        for omega in [-speed, speed] {
            let response = journal.response(omega).unwrap();
            let torque = -omega.signum() * 2. * core::f64::consts::PI
                * g.length_m * g.inner_radius_m.powi(2) * 0.25 * inner;
            assert!(response.shear_rate_per_s[0] > 16.0);
            assert!(response.shear_rate_per_s[1] < 16.0);
            close(response.shear_rate_per_s[0], inner, 2e-12);
            close(response.shear_rate_per_s[1], 16.0 * root * radius_ratio, 2e-12);
            close(response.resisting_torque_nm, torque, 2e-12);
            close(response.dissipated_power_w, -torque * omega, 2e-12);
        }
    }
}

#[test]
fn g1_exact_tangent_resolves_wide_source_support_with_a_small_solve_budget() {
    let g = geometry();
    let controls = JournalShearSolve { maximum_iterations: 12,
        absolute_angular_rate_tolerance_per_s: 1e-24,
        relative_angular_rate_tolerance: 1e-13 };
    let journal = GeneralizedNewtonianJournal::new(g,
        vec![point(0., 0.05), point(5., 0.05), point(1e6, 1e4)], controls).unwrap();
    for omega in [-2.0_f64, 2.0, -1e-9, 1e-9] {
        let response = journal.response(omega).unwrap();
        let (inner, stress) = if omega.abs() > 1e-6 {
            let inner = omega.abs() * g.outer_radius_m / (g.outer_radius_m - g.inner_radius_m);
            (inner, 0.01 * inner.powi(2))
        } else {
            let inner = 2. * omega.abs() * g.outer_radius_m.powi(2)
                / ((g.outer_radius_m - g.inner_radius_m) * (g.outer_radius_m + g.inner_radius_m));
            (inner, 0.05 * inner)
        };
        let torque = -omega.signum() * 2. * core::f64::consts::PI
            * g.length_m * g.inner_radius_m.powi(2) * stress;
        close(response.shear_rate_per_s[0], inner, 2e-12);
        close(response.resisting_torque_nm, torque, 2e-12);
        close(response.dissipated_power_w, -torque * omega, 2e-12);
        assert!(response.iterations <= controls.maximum_iterations);
    }
}

#[test]
fn g0_source_domain_nonmonotonicity_and_solver_exhaustion_refuse() {
    let good = vec![point(0., 0.05), point(5., 0.05), point(100., 1.)];
    let model = GeneralizedNewtonianJournal::new(geometry(), good.clone(), solve()).unwrap();
    assert!(model.response(100.).is_err());
    assert!(model.response(f64::NAN).is_err());
    for bad in [vec![point(0., 0.), point(100., 1.)],
        vec![point(1., 0.1), point(100., 1.)],
        vec![point(0., 1.), point(10., 0.1)],
        vec![point(0., 1.), point(10., 1.), point(10., 2.)]] {
        assert!(GeneralizedNewtonianJournal::new(geometry(), bad, solve()).is_err());
    }
    let mut exhausted = solve();
    exhausted.maximum_iterations = 1;
    assert!(GeneralizedNewtonianJournal::new(geometry(), good, exhausted).unwrap().response(2.).is_err());
}
