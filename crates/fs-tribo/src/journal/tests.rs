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
