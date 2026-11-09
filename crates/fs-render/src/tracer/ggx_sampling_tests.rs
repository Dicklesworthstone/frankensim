//! G1/G3 controls for the reflection-only bounded VNDF. The finite-direction
//! angular integrals do not use the sampling warp, and the normal-incidence
//! variance oracle has a closed form independent of both implementations.

use super::*;

fn local_direction(cosine: f64, azimuth: f64) -> [f64; 3] {
    let sine = (1.0 - cosine * cosine).sqrt();
    [sine * det::cos(azimuth), sine * det::sin(azimuth), cosine]
}

fn half_vector(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    let h = core::array::from_fn::<_, 3, _>(|i| a[i] + b[i]);
    let length = h.iter().map(|x| x * x).sum::<f64>().sqrt();
    h.map(|x| x / length)
}

// Independent slope-form single-scattering GGX f*cos_i for unit Fresnel.
fn projected_brdf(ax: f64, ay: f64, wo: [f64; 3], wi: [f64; 3]) -> f64 {
    let h = half_vector(wo, wi);
    let slope_x = h[0] / (ax * h[2]);
    let slope_y = h[1] / (ay * h[2]);
    let denominator = 1.0 + slope_x * slope_x + slope_y * slope_y;
    let d = 1.0 / (PI * ax * ay * h[2].powi(4) * denominator * denominator);
    let lambda = |v: [f64; 3]| {
        let tangent2 = ((ax * v[0]).powi(2) + (ay * v[1]).powi(2)) / v[2].powi(2);
        0.5 * ((1.0 + tangent2).sqrt() - 1.0)
    };
    d / (4.0 * wo[2] * (1.0 + lambda(wo) + lambda(wi)))
}

#[test]
fn unit_roughness_produces_uniform_upper_hemisphere_without_rejection() {
    const SIDE: u32 = 32;
    for cosine in [0.01, 0.4, 1.0] {
        let local = local_direction(cosine, 0.7);
        for n in [Vec3::new(0.0, 0.0, 1.0), Vec3::new(0.0, 0.0, -1.0),
            Vec3::new(0.6, 0.0, 0.8)] {
            let wo = to_world_all_sphere(n, local);
            for y in 0..SIDE {
                for x in 0..SIDE {
                    let u = (f64::from(x) + 0.5) / f64::from(SIDE);
                    let v = (f64::from(y) + 0.5) / f64::from(SIDE);
                    let (wi, pdf) = sample_ggx_vndf_reflection(n, wo, 1.0, u, v)
                        .expect("unit-roughness cap is exactly the contributing hemisphere");
                    assert!(n.dot(wi) > 0.0);
                    assert!((wi.norm() - 1.0).abs() < 3e-14);
                    assert!((pdf * 2.0 * PI - 1.0).abs() < 3e-14);
                    assert!((n.dot(wi) - (1.0 - u)).abs() < 3e-14);
                }
            }
        }
    }
}

#[test]
fn anisotropic_sample_mass_moments_and_power_match_angular_quadrature() {
    const SIDE: u32 = 256;
    let count = f64::from(SIDE * SIDE);
    for (ax, ay) in [(0.3, 0.7), (0.7, 0.3), (0.6, 0.6), (1.0, 0.4)] {
        for cosine in [0.15, 0.65, 1.0] {
            let wo = local_direction(cosine, 0.43);
            let mut sampled = [0.0; 5];
            let mut angular = [0.0; 5];
            for y in 0..SIDE {
                for x in 0..SIDE {
                    let u = (f64::from(x) + 0.5) / f64::from(SIDE);
                    let v = (f64::from(y) + 0.5) / f64::from(SIDE);
                    if let Some((wi, pdf)) = sample_ggx_vndf_reflection_aniso(wo, ax, ay, u, v) {
                        let evaluated = ggx_vndf_reflection_pdf_aniso(ax, ay, wo, half_vector(wo, wi));
                        assert!((pdf / evaluated - 1.0).abs() < 3e-12);
                        assert!((wi.iter().map(|x| x*x).sum::<f64>() - 1.0).abs() < 3e-14);
                        let values = [1.0, wi[0], wi[1], wi[2], projected_brdf(ax, ay, wo, wi) / pdf];
                        for i in 0..5 { sampled[i] += values[i] / count; }
                    }
                    let wi = local_direction(u, 2.0 * PI * v);
                    let pdf = ggx_vndf_reflection_pdf_aniso(ax, ay, wo, half_vector(wo, wi));
                    let values = [pdf, pdf * wi[0], pdf * wi[1], pdf * wi[2], projected_brdf(ax, ay, wo, wi)];
                    for i in 0..5 { angular[i] += values[i] * 2.0 * PI / count; }
                }
            }
            assert!(sampled[0] > 0.0 && sampled[0] <= 1.0);
            assert!(sampled[4] > 0.0 && sampled[4] < 1.0);
            for i in 0..5 {
                assert!((sampled[i] - angular[i]).abs() < 1.5e-3,
                    "alpha=({ax},{ay}), mu={cosine}, moment={i}: sample={} integral={}",
                    sampled[i], angular[i]);
            }
        }
    }
}

#[test]
fn bounded_proposal_preserves_physical_scalar_and_polarized_normalization() {
    for (ax, ay) in [(0.0001, 0.0001), (0.2, 0.7), (0.7, 0.2), (1.0, 1.0)] {
        for cosine in [0.001, 0.2, 1.0] {
            let wo = local_direction(cosine, 0.29);
            let masking_o = smith_g1_aniso(ax, ay, wo);
            let mass = ggx_reflection_cap_mass(ax, ay, cosine, masking_o);
            assert!((0.5..=1.0).contains(&mass));
            for (u, v) in [(0.0, 0.0), (0.31, 0.73), (0.99, 0.41)] {
                let Some((wi, pdf)) = sample_ggx_vndf_reflection_aniso(wo, ax, ay, u, v)
                    else { continue; };
                let weighted = mass * crate::smith_correlated_weight(
                    masking_o, smith_g1_aniso(ax, ay, wi),
                );
                let expected = projected_brdf(ax, ay, wo, wi) / pdf;
                assert!((weighted / expected - 1.0).abs() < 2e-10);
                // The reverse density belongs to the reversed incident view;
                // matching forward/reverse PDFs is not a reciprocity condition.
                let reverse = ggx_vndf_reflection_pdf_aniso(ax, ay, wi, half_vector(wo, wi));
                let reverse_mass = ggx_reflection_cap_mass(ax, ay, wi[2], smith_g1_aniso(ax, ay, wi));
                let forward_brdf = pdf * weighted / wi[2];
                let reverse_brdf = reverse * reverse_mass * crate::smith_correlated_weight(
                    smith_g1_aniso(ax, ay, wi), masking_o,
                ) / cosine;
                assert!((forward_brdf / reverse_brdf - 1.0).abs() < 2e-10);
            }
        }
    }
}

#[test]
fn rough_white_furnace_retains_analytic_mean_with_lower_equal_count_variance() {
    // At alpha=1 and normal view the old VNDF is uniform on the sphere;
    // half of its draws are rejected. The new cap is the upper hemisphere.
    // With unit Fresnel, f*cos_i = mu/(2*pi*(1+mu)). Thus the mean is
    // 1-ln(2), and the new second moment is 3/2-2*ln(2); the old second
    // moment is twice that. This is the unchanged SINGLE-scattered energy.
    const SAMPLES: u32 = 1 << 15;
    let n = Vec3::new(0.0, 0.0, 1.0);
    let mut sum = 0.0;
    let mut squares = 0.0;
    let mut old_sum = 0.0;
    let mut old_squares = 0.0;
    let mut old_accepted = 0;
    for i in 0..SAMPLES {
        let u = (f64::from(i) + 0.5) / f64::from(SAMPLES);
        let (wi, pdf) = sample_ggx_vndf_reflection(n, n, 1.0, u, 0.37).unwrap();
        let value = wi.z / (2.0 * PI * (1.0 + wi.z)) / pdf;
        sum += value;
        squares += value * value;
        // Original VNDF at normal incidence is the original NDF warp.
        let mu = 1.0 - 2.0 * u;
        if mu > 0.0 {
            old_accepted += 1;
            let value = 2.0 * mu / (1.0 + mu);
            old_sum += value;
            old_squares += value * value;
        }
    }
    let count = f64::from(SAMPLES);
    let mean = sum / count;
    let old_mean = old_sum / count;
    let variance = squares / count - mean * mean;
    let old_variance = old_squares / count - old_mean * old_mean;
    let expected_mean = 1.0 - core::f64::consts::LN_2;
    let expected_second = 1.5 - 2.0 * core::f64::consts::LN_2;
    assert_eq!(old_accepted, SAMPLES / 2);
    assert!((mean - expected_mean).abs() < 2e-10);
    assert!((old_mean - expected_mean).abs() < 2e-10);
    assert!((squares / count - expected_second).abs() < 2e-10);
    assert!(variance < 0.15 * old_variance);
    eprintln!("bounded GGX furnace: mean={mean:.12}, old_mean={old_mean:.12}, variance={variance:.12}, old_variance={old_variance:.12}");
}

#[test]
fn cap_admission_keeps_rejections_and_does_not_resample() {
    let wo = local_direction(0.1, 0.0);
    assert!(sample_ggx_vndf_reflection_aniso(wo, 0.1, 0.1, 0.8, 0.5).is_none());
    for (u, v) in [(-0.1, 0.5), (1.0, 0.5), (0.5, f64::NAN), (0.5, 1.0)] {
        assert!(sample_ggx_vndf_reflection_aniso(wo, 0.3, 0.7, u, v).is_none());
    }
    for a in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(sample_ggx_vndf_reflection_aniso(wo, a, 0.7, 0.3, 0.7).is_none());
    }
}
