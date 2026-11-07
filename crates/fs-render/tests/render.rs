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
fn g0_mis_preserves_power_across_density_scales_and_sample_counts() {
    // Common changes of measure cannot change a strategy's relative mass.
    // Include subnormal densities, underflowing squares, overflowing squares,
    // and count products that exceed binary64 despite finite input pdfs.
    for scale in [f64::from_bits(16), 1e-300, 1e-170, 1.0, 1e170, 1e300] {
        let (pf, pg) = (2.0 * scale, 5.0 * scale);
        for (heuristic, expected) in [
            (balance_heuristic as fn(u32, f64, u32, f64) -> f64, 6.0 / 41.0),
            (power_heuristic, 36.0 / 1261.0),
        ] {
            let forward = heuristic(3, pf, 7, pg);
            let reverse = heuristic(7, pg, 3, pf);
            assert!((forward - expected).abs() < 4e-16, "scale={scale}: {forward}");
            assert!((forward + reverse - 1.0).abs() < 4e-16);
        }
    }
    for heuristic in [balance_heuristic, power_heuristic] {
        assert_eq!(heuristic(u32::MAX, f64::MAX, u32::MAX, f64::MAX), 0.5);
        assert_eq!(heuristic(1, f64::MAX, 1, f64::MAX), 0.5);
        assert_eq!(heuristic(1, f64::from_bits(1), 1, f64::from_bits(1)), 0.5);
    }
    // A small strategy can still have a representable weight. In particular,
    // neither its pdf nor its squared mass may be dropped before normalization.
    let weak_balance = balance_heuristic(u32::MAX, 1e-20, 2, 1e308);
    let expected = (f64::from(u32::MAX) * 1e-20 / 1e308) * 0.5;
    assert!(weak_balance > 0.0);
    assert!((weak_balance - expected).abs() <= f64::from_bits(1));
    assert!((power_heuristic(1, 1e-200, 1, 1e-100) / 1e-200 - 1.0).abs() < 1e-14);
}

#[test]
fn g0_mis_inactive_proposals_and_invalid_densities_are_explicit() {
    for heuristic in [balance_heuristic, power_heuristic] {
        assert_eq!(heuristic(1, 0.0, 1, 0.0), 0.0);
        assert_eq!(heuristic(0, f64::NAN, 0, f64::INFINITY), 0.0);
        assert_eq!(heuristic(0, f64::NAN, 1, 1e-300), 0.0);
        assert_eq!(heuristic(1, 1e-300, 0, f64::INFINITY), 1.0);
        for invalid in [-1.0, f64::NAN, f64::INFINITY] {
            assert!(heuristic(1, invalid, 1, 1.0).is_nan());
            assert!(heuristic(1, 1.0, 1, invalid).is_nan());
        }
    }
}

#[test]
fn g3_mis_constant_radiance_integral_is_invariant_to_change_of_measure() {
    // Integrate f(x)=1/width on [0,width] with uniform and linearly weighted
    // proposals. This independent unit-power oracle exercises both MIS legs;
    // neither a very sharp nor a very diffuse proposal may erase the energy.
    for heuristic in [balance_heuristic, power_heuristic] {
        for width in [1e-200, 1.0, 1e200] {
            let uniform_pdf = 1.0 / width;
            let mut estimate = 0.0;
            let count = 4096;
            for i in 1..=count {
                let uniform = radical_inverse(2, i);
                let linear = radical_inverse(3, i).sqrt();
                estimate += heuristic(1, uniform_pdf, 1, 2.0 * uniform / width);
                estimate += heuristic(1, 2.0 * linear / width, 1, uniform_pdf) / (2.0 * linear);
            }
            estimate /= count as f64;
            assert!((estimate - 1.0).abs() < 1e-3, "width={width}: {estimate}");
        }
    }
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
