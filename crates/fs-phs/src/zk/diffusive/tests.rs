use super::*;

fn model(poles: usize, error: f64) -> Result<ZkDiffusive, PhsError> {
    // tau = 1, so the test band is exactly viscous shear rv in [1,30].
    ZkDiffusive::new(1.0, 1.0, 1.0, 1.4, 0.71, poles, [1.0, 900.0], error)
}

#[test]
fn zk_pole_ladder_resolves_complex_storage_and_real_loss() {
    let mut previous_complex = f64::INFINITY;
    let mut previous_real = f64::INFINITY;
    for poles in [2, 4, 8, 16] {
        let model = model(poles, 0.95).unwrap();
        assert!(model.max_complex_error() < previous_complex);
        assert!(model.max_real_loss_error() < previous_real);
        previous_complex = model.max_complex_error();
        previous_real = model.max_real_loss_error();
        for i in 0..129 {
            let w = if i == 128 {
                900.0
            } else {
                det::pow(900.0, f64::from(i) / 128.0)
            };
            let (z, y) = model.response_at(w).unwrap();
            assert!(z.re > 0.0 && z.im > 0.0 && y.re > 0.0 && y.im > 0.0);
        }
        println!("poles={poles} complex={previous_complex:e} real={previous_real:e}");
    }
    assert!(previous_complex < 1e-4 && previous_real < 1e-4);
    assert!(model(4, 1e-4).is_err());
    assert!(model(16, 1e-4).is_ok());
}

#[test]
fn pole_storage_keeps_poiseuille_inertance_and_isothermal_compliance() {
    for poles in [1, 2, 4, 8, 16, 32, MAX_ZK_DIFFUSIVE_POLES] {
        let model =
            ZkDiffusive::new(0.001, 1.2, 1.8e-5, 1.4, 0.71, poles, [0.1, 1.0], 0.01).unwrap();
        let viscous_storage: f64 = model
            .series_terms()
            .iter()
            .map(|t| t.gain_per_s / t.rate_per_s)
            .sum();
        let thermal_storage: f64 = model
            .shunt_terms()
            .iter()
            .map(|t| t.gain_per_s / t.rate_per_s)
            .sum();
        assert!((1.0 + viscous_storage - 4.0 / 3.0).abs() < 2e-12);
        assert!((1.0 + thermal_storage - 1.4).abs() < 2e-12);
        assert!((model.series_drag_per_s() - 120.0).abs() < 1e-12);
        assert!(
            model
                .series_terms()
                .iter()
                .chain(model.shunt_terms())
                .all(|t| t.gain_per_s > 0.0 && t.rate_per_s > 0.0)
        );
    }
}

#[test]
fn material_diffusion_scale_changes_both_memory_and_pressure_response() {
    let a = ZkDiffusive::new(0.001, 1.2, 1.8e-5, 1.4, 0.71, 16, [10.0, 5000.0], 1e-4).unwrap();
    let b = ZkDiffusive::new(0.001, 0.6, 1.8e-5, 1.4, 0.71, 16, [10.0, 5000.0], 1e-4).unwrap();
    assert_eq!(b.series_drag_per_s(), 2.0 * a.series_drag_per_s());
    for (a, b) in a.series_terms().iter().zip(b.series_terms()) {
        assert_eq!(b.rate_per_s, 2.0 * a.rate_per_s);
        assert_eq!(b.gain_per_s, 2.0 * a.gain_per_s);
    }
    let (az, ay) = a.response_at(1000.0).unwrap();
    let (bz, by) = b.response_at(1000.0).unwrap();
    assert!((az - bz).abs() > 10.0 && (ay - by).abs() > 10.0);
    let isothermal =
        ZkDiffusive::new(0.001, 1.2, 1.8e-5, 1.0, 0.71, 16, [10.0, 5000.0], 1e-4).unwrap();
    assert!(isothermal.shunt_terms().is_empty());
    assert_eq!(
        isothermal.response_at(1000.0).unwrap().1,
        C64::new(0.0, 1000.0)
    );
}

#[test]
fn invalid_or_unresolved_diffusion_never_becomes_lossless() {
    for count in [0, MAX_ZK_DIFFUSIVE_POLES + 1] {
        assert!(model(count, 0.1).is_err());
    }
    for limit in [0.0, 1.0, f64::NAN, f64::INFINITY] {
        assert!(model(16, limit).is_err());
    }
    for band in [[0.0, 1.0], [1.0, 1.0], [2.0, 1.0], [1.0, f64::INFINITY]] {
        assert!(ZkDiffusive::new(1.0, 1.0, 1.0, 1.4, 0.71, 16, band, 0.1).is_err());
    }
    for physical in [
        [0.0, 1.0, 1.0, 1.4, 0.71],
        [1.0, 0.0, 1.0, 1.4, 0.71],
        [1.0, 1.0, 0.0, 1.4, 0.71],
        [1.0, 1.0, 1.0, 0.99, 0.71],
        [1.0, 1.0, 1.0, 1.4, 0.0],
        [f64::NAN, 1.0, 1.0, 1.4, 0.71],
    ] {
        let [a, rho, mu, gamma, pr] = physical;
        assert!(ZkDiffusive::new(a, rho, mu, gamma, pr, 16, [1.0, 10.0], 0.1).is_err());
    }
    let model = model(16, 1e-4).unwrap();
    assert!(model.response_at(0.0).is_err());
    assert!(model.response_at(901.0).is_err());
    assert!(model.response_at(f64::NAN).is_err());
    assert!(ZkDiffusive::new(1.0, 1.0, 1.0, 1.4, 0.71, 16, [1e-18, 1e-15], 0.1).is_err()); // existing exact-Bessel shear domain
}
