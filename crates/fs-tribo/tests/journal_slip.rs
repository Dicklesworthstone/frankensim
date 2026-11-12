//! G1/G3: cylindrical wall slip, source rheology and separate physical heat.
use fs_tribo::journal::{
    GeneralizedNewtonianJournal, JournalGeometry, JournalShearSolve, JournalViscosityKnot,
};
use fs_tribo::journal::slip::{JournalWallSlipLengths, NavierSlipJournal};

fn geometry() -> JournalGeometry {
    JournalGeometry { inner_radius_m: 0.1, outer_radius_m: 0.14, length_m: 0.2 }
}

fn source(g: JournalGeometry, knots: &[(f64, f64)]) -> GeneralizedNewtonianJournal {
    GeneralizedNewtonianJournal::new(g, knots.iter().map(|&(rate, mu)| JournalViscosityKnot {
        shear_rate_per_s: rate, dynamic_viscosity_pa_s: mu,
    }).collect(), JournalShearSolve { maximum_iterations: 96,
        absolute_angular_rate_tolerance_per_s: 1e-13,
        relative_angular_rate_tolerance: 1e-13 }).unwrap()
}

fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance * a.abs().max(b.abs()).max(1e-30), "{a:e} != {b:e}");
}

#[test]
fn g1_curved_newtonian_slip_matches_exact_torque_and_three_heat_integrals() {
    let g = geometry();
    let mu = 0.2;
    let slip = JournalWallSlipLengths { inner_m: 0.003, outer_m: 0.007 };
    let journal = NavierSlipJournal::new(source(g, &[(0., mu), (1000., mu)]), slip).unwrap();
    let resistance = g.inner_radius_m.recip().powi(2) - g.outer_radius_m.recip().powi(2)
        + 2. * slip.inner_m / g.inner_radius_m.powi(3)
        + 2. * slip.outer_m / g.outer_radius_m.powi(3);
    for omega in [-4., 4.] {
        let result = journal.response(omega).unwrap();
        let torque = 4. * core::f64::consts::PI * mu * g.length_m * omega.abs() / resistance;
        close(result.shear.resisting_torque_nm, -omega.signum() * torque, 4e-15);
        close(result.shear.dissipated_power_w, torque * omega.abs(), 4e-15);
        let scale = torque.powi(2) / (2. * core::f64::consts::PI * mu * g.length_m);
        close(result.bulk_dissipated_power_w, 0.5 * scale
            * (g.inner_radius_m.recip().powi(2) - g.outer_radius_m.recip().powi(2)), 5e-15);
        close(result.wall_dissipated_power_w[0], scale * slip.inner_m / g.inner_radius_m.powi(3), 5e-15);
        close(result.wall_dissipated_power_w[1], scale * slip.outer_m / g.outer_radius_m.powi(3), 5e-15);
        assert_eq!(result.shear.iterations, 0);
        assert!(result.shear.resisting_torque_nm.abs()
            < journal.journal().response(omega).unwrap().resisting_torque_nm.abs());
        assert!(result.wall_slip_velocity_m_per_s[0] * omega < 0.);
        assert!(result.wall_slip_velocity_m_per_s[1] * omega > 0.);
        assert!(result.power_defect_w.abs() <= 1e-14 * result.shear.dissipated_power_w);
    }
}

#[test]
fn g1_g3_nonlinear_wall_viscosities_recover_power_law_couette_slip() {
    let g = geometry();
    let k = 0.03;
    let slip = JournalWallSlipLengths { inner_m: 0.01, outer_m: 0.025 };
    // The positive zero-shear plateau is outside this solution's rate range.
    let journal = NavierSlipJournal::new(source(g, &[(0., 5. * k), (5., 5. * k), (100., 100. * k)]), slip).unwrap();
    let inner_rate = 20.;
    let outer_rate = inner_rate * g.inner_radius_m / g.outer_radius_m;
    let bulk_speed = inner_rate * (1. - g.inner_radius_m / g.outer_radius_m);
    let omega = bulk_speed + slip.inner_m * inner_rate / g.inner_radius_m
        + slip.outer_m * outer_rate / g.outer_radius_m;
    let torque = 2. * core::f64::consts::PI * g.length_m
        * g.inner_radius_m.powi(2) * k * inner_rate.powi(2);
    let forward = journal.response(omega).unwrap();
    let reverse = journal.response(-omega).unwrap();
    close(-forward.shear.resisting_torque_nm, torque, 1e-12);
    close(forward.shear.shear_rate_per_s[0], inner_rate, 1e-12);
    close(forward.shear.shear_rate_per_s[1], outer_rate, 1e-12);
    close(forward.bulk_dissipated_power_w, torque * bulk_speed, 1e-12);
    close(forward.wall_dissipated_power_w[0], torque * slip.inner_m * inner_rate / g.inner_radius_m, 1e-12);
    close(forward.wall_dissipated_power_w[1], torque * slip.outer_m * outer_rate / g.outer_radius_m, 1e-12);
    assert_eq!(forward.shear.resisting_torque_nm, -reverse.shear.resisting_torque_nm);
    assert_eq!(forward.bulk_dissipated_power_w, reverse.bulk_dissipated_power_w);
    assert_eq!(forward.wall_dissipated_power_w, reverse.wall_dissipated_power_w);
    assert_eq!(forward.wall_slip_velocity_m_per_s.map(|v| -v), reverse.wall_slip_velocity_m_per_s);
}

#[test]
fn g1_g3_thinning_source_crosses_knots_and_matches_independent_volume_dissipation() {
    let g = JournalGeometry { outer_radius_m: 0.3, ..geometry() };
    let slip = JournalWallSlipLengths { inner_m: 0.01, outer_m: 0.04 };
    let journal = NavierSlipJournal::new(source(g, &[(0., 0.2), (3., 0.2), (8., 0.14), (30., 0.14)]), slip).unwrap();
    let tau_i = 1.2;
    let inverse = |stress: f64| {
        let mut low = 0.;
        let mut high = 30.;
        for _ in 0..60 {
            let rate = 0.5 * (low + high);
            let mu = if rate < 3. { 0.2 } else if rate < 8. { 0.2 - 0.012 * (rate - 3.) } else { 0.14 };
            if rate * mu < stress { low = rate; } else { high = rate; }
        }
        0.5 * (low + high)
    };
    // Composite Simpson on the three radial pieces between the walls and the
    // radii of the 8/s and 3/s source knots, where gamma(r)/r is smooth. With
    // the 20,000-cell midpoint rule this replaces, the library's inner-wall
    // loss at the reference speed was 9.8e-9 below the loss at the chosen
    // stress, beyond the 6e-9 checked below. By a replay of both computations
    // the library's bulk speed matches the closed form to an ulp, the midpoint
    // sum left it 8.1e-9 low (its h^2 error on the steep inner-wall
    // integrand), and 2,000 intervals per piece leave about 4e-13.
    let knot_radius = |stress: f64| g.inner_radius_m * (tau_i / stress).sqrt();
    let radii = [g.inner_radius_m, knot_radius(8. * 0.14), knot_radius(3. * 0.2), g.outer_radius_m];
    let mut bulk_speed = 0.;
    let mut bulk_heat = 0.;
    for piece in radii.windows(2) {
        let n = 2_000_u32;
        let dr = (piece[1] - piece[0]) / f64::from(n);
        for i in 0..=n {
            let r = piece[0] + f64::from(i) * dr;
            let simpson = if i == 0 || i == n { 1. } else if i % 2 == 1 { 4. } else { 2. };
            let weight = simpson * dr / 3.;
            let stress = tau_i * (g.inner_radius_m / r).powi(2);
            let rate = inverse(stress);
            bulk_speed += weight * rate / r;
            bulk_heat += weight * 2. * core::f64::consts::PI * g.length_m * stress * rate * r;
        }
    }
    let rates = [inverse(tau_i), inverse(tau_i * (g.inner_radius_m / g.outer_radius_m).powi(2))];
    let omega = bulk_speed + slip.inner_m * rates[0] / g.inner_radius_m
        + slip.outer_m * rates[1] / g.outer_radius_m;
    let result = journal.response(omega).unwrap();
    let torque = 2. * core::f64::consts::PI * g.length_m * g.inner_radius_m.powi(2) * tau_i;
    close(-result.shear.resisting_torque_nm, torque, 6e-9);
    close(result.bulk_dissipated_power_w, bulk_heat, 6e-9);
    assert!(result.shear.shear_rate_per_s[0] > 8. && result.shear.shear_rate_per_s[1] < 3.);
    close(result.wall_dissipated_power_w[0], torque * slip.inner_m * rates[0] / g.inner_radius_m, 6e-9);
    close(result.wall_dissipated_power_w[1], torque * slip.outer_m * rates[1] / g.outer_radius_m, 6e-9);
    assert!(result.power_defect_w.abs() < 1e-12 * result.shear.dissipated_power_w);
}

#[test]
fn g3_zero_slip_retains_original_response_and_rest_has_no_heat() {
    let journal = source(geometry(), &[(0., 0.2), (10., 0.18), (100., 0.18)]);
    let slipping = NavierSlipJournal::new(journal.clone(), JournalWallSlipLengths::default()).unwrap();
    for omega in [-2., 0., 2.] {
        let result = slipping.response(omega).unwrap();
        assert_eq!(result.shear, journal.response(omega).unwrap());
        assert_eq!(result.wall_dissipated_power_w, [0.; 2]);
        assert_eq!(result.wall_slip_velocity_m_per_s, [0.; 2]);
    }
    let slipping = NavierSlipJournal::new(journal, JournalWallSlipLengths { inner_m: 0.01, outer_m: 0.02 }).unwrap();
    let rest = slipping.response(0.).unwrap();
    assert_eq!(rest.shear.resisting_torque_nm, 0.);
    assert_eq!(rest.bulk_dissipated_power_w, 0.);
    assert_eq!(rest.wall_dissipated_power_w, [0.; 2]);
    assert_eq!(rest.power_defect_w, 0.);
}

#[test]
fn g2_slip_preserves_material_rate_support_and_refuses_invalid_boundary_data() {
    let g = geometry();
    let journal = source(g, &[(0., 0.2), (10., 0.2)]);
    for value in [-1., f64::NAN, f64::INFINITY] {
        assert!(NavierSlipJournal::new(journal.clone(), JournalWallSlipLengths { inner_m: value, outer_m: 0. }).is_err());
        assert!(NavierSlipJournal::new(journal.clone(), JournalWallSlipLengths { inner_m: 0., outer_m: value }).is_err());
    }
    let slipping = NavierSlipJournal::new(journal.clone(), JournalWallSlipLengths { inner_m: 0.1, outer_m: 0.2 }).unwrap();
    assert!(journal.response(5.).is_err());
    assert!(slipping.response(5.).unwrap().shear.shear_rate_per_s[0] < 10.);
    assert!(slipping.response(100.).is_err());
    assert!(slipping.response(f64::NAN).is_err());
}

#[test]
fn g3_representable_narrow_gap_retains_bulk_heat_beside_interface_losses() {
    let inner = 1.0_f64;
    let outer = f64::from_bits(inner.to_bits() + 1);
    let g = JournalGeometry { inner_radius_m: inner, outer_radius_m: outer, length_m: 0.2 };
    let gap = outer - inner;
    let slipping = NavierSlipJournal::new(source(g, &[(0., 0.2), (100., 0.2)]),
        JournalWallSlipLengths { inner_m: 2. * gap, outer_m: 3. * gap }).unwrap();
    let result = slipping.response(12. * gap).unwrap();
    close(result.shear.shear_rate_per_s[0], 2., 2e-15);
    close(result.bulk_dissipated_power_w / result.shear.dissipated_power_w, 1. / 6., 3e-15);
    close(result.wall_dissipated_power_w[0] / result.shear.dissipated_power_w, 2. / 6., 3e-15);
    close(result.wall_dissipated_power_w[1] / result.shear.dissipated_power_w, 3. / 6., 3e-15);
}
