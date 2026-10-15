//! Local transport acts on the actual kernels; reference profiles are analytic.
use fs_lbm::{ThermalLbm, d3q19::ThermalLbm3};
use fs_lbm::thermal::sources::{ThermalSourceLimits, transport::RelaxationField};

fn limits(n: usize) -> ThermalSourceLimits {
    ThermalSourceLimits { max_cells: n, temperature_range: [-1.0, 2.0], balance_tolerance: 1e-11 }
}
fn near(a: f64, b: f64, tol: f64) { assert!((a-b).abs() <= tol, "{a} != {b}"); }

#[test]
fn constant_fields_preserve_original_source_step_in_both_dimensions() {
    let mut a = ThermalLbm::slab(3, 6, 0.8, 0.9, 1e-6);
    let mut b = ThermalLbm::slab(3, 6, 0.8, 0.9, 1e-6);
    a.gpop[10][0] += 0.01; b.gpop[10][0] += 0.01;
    let mut heat = vec![0.001; 24]; heat[..3].fill(0.0); heat[21..].fill(0.0);
    let forces = [[1e-6, 0.0]; 24];
    let expected = a.step_with_sources(&heat, &forces, limits(24)).unwrap();
    let actual = b.step_with_sources_and_transport(&heat, &forces,
        RelaxationField { momentum: &[0.8;24], thermal: &[0.9;24] }, limits(24)).unwrap();
    assert_eq!(expected, actual); assert_eq!(a.grid.f, b.grid.f); assert_eq!(a.gpop, b.gpop);
    assert_eq!(a.grid.tau, b.grid.tau);

    let mut a = ThermalLbm3::slab(3, 1, 6, 0.8, 0.9, 1e-6);
    let mut b = ThermalLbm3::slab(3, 1, 6, 0.8, 0.9, 1e-6);
    a.gpop[10][0] += 0.01; b.gpop[10][0] += 0.01;
    let forces = [[1e-6, 0.0, 0.0]; 24];
    let expected = a.step_with_sources(&heat, &forces, limits(24)).unwrap();
    let actual = b.step_with_sources_and_transport(&heat, &forces,
        RelaxationField { momentum: &[0.8;24], thermal: &[0.9;24] }, limits(24)).unwrap();
    assert_eq!(expected, actual); assert_eq!(a.grid.f, b.grid.f); assert_eq!(a.gpop, b.gpop);
    assert_eq!(a.grid.tau, b.grid.tau);
}

#[test]
fn local_momentum_uses_existing_collision_without_leaving_stale_rates() {
    let mut a = ThermalLbm::slab(3, 6, 0.8, 0.9, 0.0);
    let mut direct = ThermalLbm::slab(3, 6, 0.8, 0.9, 0.0);
    a.grid.f[10][0] -= 0.002; a.grid.f[10][1] += 0.001; a.grid.f[10][3] += 0.001;
    direct.grid.f.clone_from(&a.grid.f);
    let before_mass: f64 = a.grid.f.iter().flatten().sum();
    let mut momentum = vec![0.8;24]; momentum[10] = 1.2;
    direct.grid.tau.clone_from(&momentum);
    let expected = direct.step_with_sources(&[0.0;24], &[[0.0;2];24], limits(24)).unwrap();
    let actual = a.step_with_sources_and_transport(&[0.0;24], &[[0.0;2];24],
        RelaxationField { momentum: &momentum, thermal: &[0.9;24] }, limits(24)).unwrap();
    assert_eq!(a.grid.f, direct.grid.f); assert_eq!(a.gpop, direct.gpop); assert_eq!(actual, expected);
    assert_eq!(a.grid.tau, vec![0.8;24]); assert_eq!(a.tau_g, 0.9);
    near(a.grid.f.iter().flatten().sum(), before_mass, 1e-12);
    // A subsequent legacy step uses its original setup, not one stale half of
    // the last transport field. Both history and background rates are retained.
    direct.grid.tau.fill(0.8);
    a.step(); direct.step(); assert_eq!(a.grid.f, direct.grid.f); assert_eq!(a.gpop, direct.gpop);
}

fn profile_error(n: usize, dimension: usize) -> f64 {
    let cells = n+2;
    let mut a = ThermalLbm::slab(1, n, 0.8, 0.9, 0.0);
    let mut b = ThermalLbm3::slab(1, 1, n, 0.8, 0.9, 0.0);
    let momentum = vec![0.8;cells];
    let thermal: Vec<_> = (0..cells).map(|i| {
        let y = ((i as f64 - 0.5) / n as f64).clamp(0.0,1.0);
        let alpha = 0.08*(1.0+y);
        0.5 + alpha / if dimension == 2 { 1.0/3.0 } else { 0.25 }
    }).collect();
    let rates = RelaxationField { momentum: &momentum, thermal: &thermal };
    let heat = vec![0.0;cells]; let f2 = vec![[0.0;2];cells]; let f3 = vec![[0.0;3];cells];
    let mut last_bottom = 0.0;
    for _ in 0..30*n*n {
        let r = if dimension == 2 {
            a.step_with_sources_and_transport(&heat, &f2, rates, limits(cells)).unwrap()
        } else { b.step_with_sources_and_transport(&heat, &f3, rates, limits(cells)).unwrap() };
        assert!(r.defect.abs() < 1e-11); last_bottom = r.bottom_wall;
    }
    // k(y) proportional to 1+y: integral(dy/k) gives a logarithmic profile.
    let mut error = 0.0_f64;
    for i in 1..=n {
        let y = (i as f64-0.5)/n as f64;
        let exact = 1.0-(1.0+y).ln()/2.0_f64.ln();
        let actual = if dimension == 2 { a.temperature(0,i) } else { b.temperature(0,0,i) };
        error = error.max((actual-exact).abs());
    }
    near(last_bottom*n as f64, 0.08/2.0_f64.ln(), 1e-4);
    error
}

#[test]
fn variable_conductivity_converges_to_the_independent_flux_profile() {
    for d in [2,3] {
        let errors: Vec<_> = [8,16,32].map(|n| profile_error(n,d)).to_vec();
        assert!(errors[0] < 0.003);
        assert!(errors[1] < 0.3*errors[0]); assert!(errors[2] < 0.3*errors[1]);
    }
}

#[test]
fn malformed_rates_and_failed_candidates_do_not_change_any_state() {
    let mut a = ThermalLbm3::slab(2, 1, 6, 0.8, 0.9, 0.0);
    let f = a.grid.f.clone(); let g = a.gpop.clone(); let loads = a.grid.fext.clone();
    for bad in [vec![0.9;15], vec![0.9;17], vec![0.5;16], vec![f64::NAN;16], vec![f64::INFINITY;16]] {
        for thermal in [true,false] {
            let good = vec![0.8;16];
            let rates = if thermal { RelaxationField { momentum: &good, thermal: &bad } }
                else { RelaxationField { momentum: &bad, thermal: &good } };
            assert!(a.step_with_sources_and_transport(&[0.0;16], &[[0.0;3];16], rates, limits(16)).is_err());
        }
    }
    let mut hot = vec![0.0;16]; hot[6]=10.0;
    assert!(a.step_with_sources_and_transport(&hot, &[[0.0;3];16],
        RelaxationField { momentum: &[1.2;16], thermal: &[1.3;16] }, limits(16)).is_err());
    assert_eq!(a.grid.f,f); assert_eq!(a.gpop,g); assert_eq!(a.grid.fext,loads);
    assert_eq!(a.grid.tau,vec![0.8;16]); assert_eq!(a.tau_g,0.9);
}
