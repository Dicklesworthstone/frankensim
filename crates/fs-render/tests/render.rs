//! Battery for the spectral path-tracing core (fs-render). The load-bearing
//! tests are the analytic ones: the low-discrepancy sequence hits known values,
//! the Lambertian FURNACE returns exactly `albedo·radiance` (energy
//! conservation), MIS weights sum to 1 (no energy at strategy boundaries) and
//! integrate unbiasedly, and hero-wavelength integration is exact on a constant
//! spectrum and accurate on a linear one.

use fs_render::{
    Lambertian, balance_heuristic, cosine_sample_hemisphere, halton, mis_integrate_unit,
    mis_weight_sum, power_heuristic, radical_inverse, spectral_integral,
};

#[test]
fn the_radical_inverse_hits_known_values() {
    assert!((radical_inverse(2, 1) - 0.5).abs() < 1e-15);
    assert!((radical_inverse(2, 2) - 0.25).abs() < 1e-15);
    assert!((radical_inverse(2, 3) - 0.75).abs() < 1e-15);
    assert!((radical_inverse(3, 1) - 1.0 / 3.0).abs() < 1e-15);
    // halton dim 0 is the base-2 radical inverse.
    assert!((halton(0, 5) - radical_inverse(2, 5)).abs() < 1e-15);
}

#[test]
fn cosine_samples_are_unit_vectors_with_the_right_pdf() {
    for i in 1..50u64 {
        let (d, pdf) = cosine_sample_hemisphere(radical_inverse(2, i), radical_inverse(3, i));
        let len2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
        assert!((len2 - 1.0).abs() < 1e-12, "not unit: {len2}");
        assert!(d[2] >= 0.0); // upper hemisphere
        assert!((pdf - d[2] / std::f64::consts::PI).abs() < 1e-15); // cosθ/π
    }
}

#[test]
fn the_furnace_test_conserves_energy_exactly() {
    // closed uniform-emission scene: reflected radiance == albedo · incident.
    let lam = Lambertian { albedo: 0.7 };
    let est = lam.furnace_radiance(2.0, 64);
    assert!((est - 1.4).abs() < 1e-12, "furnace {est} != 1.4");
    // holds across albedos, with zero variance (cosine importance sampling).
    for &a in &[0.1, 0.5, 0.9] {
        let r = Lambertian { albedo: a }.furnace_radiance(1.0, 16);
        assert!((r - a).abs() < 1e-12, "albedo {a}: {r}");
    }
}

#[test]
fn mis_weights_sum_to_one() {
    // no energy lost or gained at the strategy boundary.
    for &(pf, pg) in &[(1.0, 2.0), (0.3, 0.7), (5.0, 0.01), (1.0, 1.0)] {
        assert!(
            (mis_weight_sum(pf, pg) - 1.0).abs() < 1e-12,
            "sum for ({pf},{pg})"
        );
    }
    // heuristics favor the higher-pdf strategy; power is sharper than balance.
    assert!(balance_heuristic(1, 3.0, 1, 1.0) > 0.5);
    assert!(power_heuristic(1, 3.0, 1, 1.0) > balance_heuristic(1, 3.0, 1, 1.0));
}

#[test]
fn mis_integration_is_unbiased() {
    // ∫₀¹ x dx = 1/2, combining uniform + linear-importance strategies.
    let est = mis_integrate_unit(|x| x, 2048);
    assert!((est - 0.5).abs() < 5e-3, "MIS estimate {est}");
    // ∫₀¹ x² dx = 1/3.
    let est2 = mis_integrate_unit(|x| x * x, 2048);
    assert!((est2 - 1.0 / 3.0).abs() < 5e-3, "MIS estimate {est2}");
}

#[test]
fn hero_wavelength_integration_is_exact_on_a_constant_and_accurate_on_a_ramp() {
    // constant spectrum: ∫ 2 dλ over [400,700] = 600, exact.
    let c = spectral_integral(|_| 2.0, 400.0, 700.0, 16);
    assert!((c - 600.0).abs() < 1e-9, "const spectral {c}");
    // linear spectrum S(λ)=λ: ∫ = (700²-400²)/2 = 165000.
    let ramp = spectral_integral(|l| l, 400.0, 700.0, 128);
    assert!(
        (ramp - 165_000.0).abs() / 165_000.0 < 0.01,
        "ramp spectral {ramp}"
    );
}

#[test]
fn rendering_primitives_are_deterministic() {
    let a = Lambertian { albedo: 0.6 }.furnace_radiance(1.0, 32);
    let b = Lambertian { albedo: 0.6 }.furnace_radiance(1.0, 32);
    assert_eq!(a.to_bits(), b.to_bits());
    assert_eq!(
        mis_integrate_unit(|x| x, 256).to_bits(),
        mis_integrate_unit(|x| x, 256).to_bits()
    );
}
