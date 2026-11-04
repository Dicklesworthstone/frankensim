use super::*;
use crate::dielectric::{
    DielectricEvent, DielectricSurface, evaluate_rough_dielectric, sample_rough_dielectric,
};

fn n() -> Vec3 {
    Vec3::new(0.0, 0.0, 1.0)
}
fn x() -> Vec3 {
    Vec3::new(1.0, 0.0, 0.0)
}
fn normalized(v: Vec3) -> Vec3 {
    v.scale(1.0 / v.norm())
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1.0e-8),
        "actual={actual:.17e}, expected={expected:.17e}, relative tolerance={tolerance}"
    );
}
fn surface(ax: f64, ay: f64) -> DielectricSurface {
    DielectricSurface::try_anisotropic_rough(ax, ay, x()).unwrap()
}

#[test]
fn equal_slopes_preserve_the_original_isotropic_bits_and_null_measure() {
    for alpha in [1.0e-4, 0.2, 1.0] {
        let isotropic = DielectricSurface::try_rough(alpha).unwrap();
        let equal =
            DielectricSurface::try_anisotropic_rough(alpha, alpha, Vec3::new(0.4, -0.3, 0.8))
                .unwrap();
        assert_eq!(equal, isotropic);
        let interface = equal.interface(550.0, false).unwrap();
        for eta in [1.0, 1.5] {
            for wo in [n(), normalized(Vec3::new(0.3, 0.2, 0.8))] {
                for branch in [0.0, 0.5, 0.999] {
                    let old = sample_rough_dielectric(n(), wo, 1.0, eta, alpha, 0.2, 0.7, branch)
                        .unwrap();
                    let new = interface
                        .sample_rough(n(), wo, 1.0, eta, alpha, 0.2, 0.7, branch)
                        .unwrap();
                    assert_eq!(old, new);
                    if let Some(sample) = new {
                        assert_eq!(sample.delta, eta == 1.0);
                        let a =
                            evaluate_rough_dielectric(n(), wo, sample.direction, 1.0, eta, alpha)
                                .unwrap();
                        let b = interface
                            .evaluate_rough(n(), wo, sample.direction, 1.0, eta, alpha)
                            .unwrap();
                        assert_eq!(a.value.to_bits(), b.value.to_bits());
                        assert_eq!(a.pdf.to_bits(), b.pdf.to_bits());
                    }
                }
            }
        }
    }
}

// Independently evaluated Walter BSDF and Smith VNDF density. This oracle uses
// tangent-plane slopes, scalar Snell/Fresnel and the explicit half-vector
// Jacobian; it never invokes the production D, G1 or sampler.
fn oracle(wo: Vec3, wi: Vec3, eta_i: f64, eta_t: f64, ax: f64, ay: f64) -> (f64, f64) {
    let reflection = wi.z > 0.0;
    let eta = if reflection { 1.0 } else { eta_t / eta_i };
    let sum = Vec3::new(wo.x + eta * wi.x, wo.y + eta * wi.y, wo.z + eta * wi.z);
    let mut m = normalized(sum);
    if m.z < 0.0 {
        m = m.scale(-1.0);
    }
    let a = wo.dot(m);
    let b = wi.dot(m);
    if a <= 0.0 || m.z <= 0.0 || (reflection && b <= 0.0) || (!reflection && b >= 0.0) {
        return (0.0, 0.0);
    }
    let sin_t2 = (eta_i / eta_t).powi(2) * (1.0 - a * a);
    let f = if sin_t2 >= 1.0 {
        1.0
    } else {
        let ct = (1.0 - sin_t2).sqrt();
        let rs = (eta_i * a - eta_t * ct) / (eta_i * a + eta_t * ct);
        let rp = (eta_t * a - eta_i * ct) / (eta_t * a + eta_i * ct);
        0.5 * (rs * rs + rp * rp)
    };
    let slope_denominator = 1.0 + (m.x / (ax * m.z)).powi(2) + (m.y / (ay * m.z)).powi(2);
    let d = 1.0 / (PI * ax * ay * m.z.powi(4) * slope_denominator.powi(2));
    let g1 = |w: Vec3| {
        2.0 / (1.0 + (1.0 + ((ax * w.x).powi(2) + (ay * w.y).powi(2)) / w.z.powi(2)).sqrt())
    };
    let visible = d * g1(wo) * a / wo.z;
    if reflection {
        (
            f * d * g1(wo) * g1(wi) / (4.0 * wo.z * wi.z),
            visible * f / (4.0 * a),
        )
    } else {
        let den = a + eta * b;
        (
            (1.0 - f) * d * g1(wo) * g1(wi) * (a * b).abs() / (wo.z * wi.z.abs() * den * den),
            visible * (1.0 - f) * eta * eta * b.abs() / (den * den),
        )
    }
}

#[test]
fn reflected_and_refracted_values_match_independent_anisotropic_walter_equations() {
    for (ax, ay) in [(0.08, 0.6), (0.7, 0.15), (1.0, 0.4)] {
        let surface = surface(ax, ay);
        let interface = surface.interface(610.0, false).unwrap();
        for (eta_i, eta_t) in [(1.0, 1.5), (1.6, 1.0)] {
            for wo in [n(), normalized(Vec3::new(0.4, -0.2, 0.6))] {
                for wi in [
                    normalized(Vec3::new(-0.2, 0.1, 0.8)),
                    normalized(Vec3::new(0.5, -0.3, 0.4)),
                    normalized(Vec3::new(-0.1, 0.15, -0.8)),
                    normalized(Vec3::new(-0.6, 0.15, -0.7)),
                ] {
                    let actual = interface
                        .evaluate_rough(n(), wo, wi, eta_i, eta_t, ax)
                        .unwrap();
                    let (value, pdf) = oracle(wo, wi, eta_i, eta_t, ax, ay);
                    close(actual.value, value, 2.0e-11);
                    close(actual.pdf, pdf, 2.0e-11);
                }
            }
        }
    }
    let glass = surface(0.08, 0.6);
    let interface = glass.interface(550.0, false).unwrap();
    let a = interface
        .evaluate_rough(
            n(),
            n(),
            normalized(Vec3::new(0.4, 0.0, 1.0)),
            1.0,
            1.5,
            0.08,
        )
        .unwrap();
    let b = interface
        .evaluate_rough(
            n(),
            n(),
            normalized(Vec3::new(0.0, 0.4, 1.0)),
            1.0,
            1.5,
            0.08,
        )
        .unwrap();
    assert!(
        b.value > 15.0 * a.value,
        "authored axes must change the actual reflected lobe"
    );
}

#[test]
fn both_hemispheres_keep_sampling_density_reciprocity_and_radiance_eta_factors() {
    let surface = surface(0.18, 0.72);
    let mut reflected = 0;
    let mut transmitted = 0;
    for (eta_i, eta_t) in [(1.0, 1.5), (1.5, 1.0)] {
        let interface = surface.interface(550.0, eta_i > eta_t).unwrap();
        let wo = normalized(Vec3::new(0.35, -0.2, 0.8));
        for i in 0..240 {
            let u = (i as f64 + 0.5) / 240.0;
            let v = ((i * 83 % 240) as f64 + 0.5) / 240.0;
            let branch = if i % 3 == 0 { 0.0 } else { 0.99 };
            let Some(sample) = interface
                .sample_rough(n(), wo, eta_i, eta_t, 0.18, u, v, branch)
                .unwrap()
            else {
                continue;
            };
            let eval = interface
                .evaluate_rough(n(), wo, sample.direction, eta_i, eta_t, 0.18)
                .unwrap();
            // Public replay renormalizes the original input once more than
            // the sampler's already admitted frame, allowing roundoff only.
            close(sample.value, eval.value, 2.0e-12);
            close(sample.pdf, eval.pdf, 2.0e-12);
            let transmission = sample.event == DielectricEvent::Transmission;
            let reverse = if transmission {
                transmitted += 1;
                surface
                    .interface(550.0, eta_i < eta_t)
                    .unwrap()
                    .evaluate_rough(n().scale(-1.0), sample.direction, wo, eta_t, eta_i, 0.18)
                    .unwrap()
            } else {
                reflected += 1;
                interface
                    .evaluate_rough(n(), sample.direction, wo, eta_i, eta_t, 0.18)
                    .unwrap()
            };
            let factor = if transmission {
                (eta_t / eta_i).powi(2)
            } else {
                1.0
            };
            close(reverse.value, sample.value * factor, 2.0e-11);
            assert!(reverse.pdf > 0.0 && reverse.pdf.is_finite());
            let wi = sample.direction;
            let g1 = 2.0
                / (1.0
                    + (1.0 + ((0.18 * wi.x).powi(2) + (0.72 * wi.y).powi(2)) / wi.z.powi(2))
                        .sqrt());
            close(sample.radiance_weight * factor, g1, 2.0e-11);
        }
    }
    assert!(reflected > 50 && transmitted > 50);
}

fn radical(mut index: usize, base: usize) -> f64 {
    let mut value = 0.0;
    let mut weight = 1.0 / base as f64;
    while index > 0 {
        value += weight * (index % base) as f64;
        index /= base;
        weight /= base as f64;
    }
    value
}

#[test]
fn vndf_density_mass_and_furnace_power_match_independent_sphere_quadrature() {
    let surface = surface(0.48, 0.82);
    let interface = surface.interface(550.0, false).unwrap();
    let wo = normalized(Vec3::new(0.5, 0.3, 0.7));
    let (eta_i, eta_t) = (1.0_f64, 1.8_f64);
    let mut pdf_mass = 0.0;
    let mut power = 0.0;
    let (nz, na) = (320, 640);
    let domega = 4.0 * PI / (nz * na) as f64;
    for iz in 0..nz {
        let z = -1.0 + 2.0 * (iz as f64 + 0.5) / nz as f64;
        let r = (1.0 - z * z).sqrt();
        for ia in 0..na {
            let angle = 2.0 * PI * (ia as f64 + 0.5) / na as f64;
            let wi = Vec3::new(r * angle.cos(), r * angle.sin(), z);
            let (value, pdf) = oracle(wo, wi, eta_i, eta_t, 0.48, 0.82);
            pdf_mass += pdf * domega;
            power += value
                * z.abs()
                * domega
                * if z < 0.0 {
                    (eta_t / eta_i).powi(2)
                } else {
                    1.0
                };
        }
    }
    let samples = 65_536;
    let mut accepted = 0;
    let mut sampled_power = 0.0;
    for i in 1..=samples {
        if let Some(sample) = interface
            .sample_rough(
                n(),
                wo,
                eta_i,
                eta_t,
                0.48,
                radical(i, 2),
                radical(i, 3),
                radical(i, 5),
            )
            .unwrap()
        {
            accepted += 1;
            sampled_power += sample.radiance_weight
                * if sample.event == DielectricEvent::Transmission {
                    (eta_t / eta_i).powi(2)
                } else {
                    1.0
                };
        }
    }
    close(accepted as f64 / samples as f64, pdf_mass, 2.0e-3);
    close(sampled_power / samples as f64, power, 2.0e-3);
    assert!(
        power > 0.8 && power < 1.0,
        "single-scatter furnace power={power}"
    );
    assert!(
        pdf_mass > power && pdf_mass < 1.0,
        "unresampled rejected facets retain probability mass"
    );
}

fn rotate(v: Vec3) -> Vec3 {
    // Exact proper cyclic rotation avoids a second rotation-matrix oracle.
    Vec3::new(v.z, v.x, v.y)
}

#[test]
fn authored_axes_follow_rigid_rotation_swap_and_boundary_reversal() {
    let surface = surface(0.12, 0.65);
    let rotated = surface.with_tangent_axis(rotate(x())).unwrap();
    let swapped =
        DielectricSurface::try_anisotropic_rough(0.65, 0.12, Vec3::new(0.0, 1.0, 0.0)).unwrap();
    let wo = normalized(Vec3::new(0.4, 0.2, 0.7));
    for i in 0..100 {
        let u = (i as f64 + 0.5) / 100.0;
        for branch in [0.0, 0.9] {
            let a = surface
                .interface(550.0, false)
                .unwrap()
                .sample_rough(n(), wo, 1.0, 1.5, 0.12, u, 0.37, branch)
                .unwrap();
            let b = rotated
                .interface(550.0, false)
                .unwrap()
                .sample_rough(rotate(n()), rotate(wo), 1.0, 1.5, 0.12, u, 0.37, branch)
                .unwrap();
            assert_eq!(a.is_some(), b.is_some());
            if let (Some(a), Some(b)) = (a, b) {
                assert!((rotate(a.direction).dot(b.direction) - 1.0).abs() < 2.0e-14);
                close(a.radiance_weight, b.radiance_weight, 2.0e-12);
                close(a.pdf, b.pdf, 2.0e-12);
                let same = swapped
                    .interface(550.0, false)
                    .unwrap()
                    .evaluate_rough(n(), wo, a.direction, 1.0, 1.5, 0.65)
                    .unwrap();
                close(a.value, same.value, 2.0e-12);
                close(a.pdf, same.pdf, 2.0e-12);
            }
        }
    }
    let axes = surface.anisotropy().unwrap();
    let (t, b) = axes.tangent_frame(n()).unwrap();
    let (rt, rb) = axes.tangent_frame(n().scale(-1.0)).unwrap();
    assert_eq!(t, rt);
    close(b.dot(rb), -1.0, 1.0e-15);
}

#[test]
fn absorbing_coating_retains_survival_and_reversed_order_on_anisotropic_facets() {
    use crate::conductor::{
        CauchyCoatingLayer, CoherentPeriodicCoating, ConductorDataStatus, ConductorSource,
    };
    use crate::dielectric::{CauchyIor, DielectricCoating};
    let layers = [
        CauchyCoatingLayer::try_new(CauchyIor::try_constant(1.3).unwrap(), 0.02, 120.0e-9).unwrap(),
        CauchyCoatingLayer::try_new(CauchyIor::try_constant(2.0).unwrap(), 0.04, 80.0e-9).unwrap(),
    ];
    let source = ConductorSource::try_new(
        fs_blake3::ContentHash([87; 32]),
        ConductorDataStatus::Representative,
    )
    .unwrap();
    let coating =
        DielectricCoating::try_new(CoherentPeriodicCoating::try_new(&layers, 1, source).unwrap())
            .unwrap();
    let surface = surface(0.2, 0.6).with_coating(coating);
    for (eta_i, eta_t, reverse) in [(1.0, 1.5, false), (1.5, 1.0, true)] {
        let interface = surface.interface(550.0, reverse).unwrap();
        let wo = normalized(Vec3::new(0.2, -0.1, 0.9));
        for branch in [0.0, 0.999] {
            let Some(sample) = interface
                .sample_rough(n(), wo, eta_i, eta_t, 0.2, 0.3, 0.2, branch)
                .unwrap()
            else {
                panic!("supported facet");
            };
            let transmission = sample.event == DielectricEvent::Transmission;
            let ratio = if transmission { eta_t / eta_i } else { 1.0 };
            let mut m = normalized(Vec3::new(
                wo.x + ratio * sample.direction.x,
                wo.y + ratio * sample.direction.y,
                wo.z + ratio * sample.direction.z,
            ));
            if m.z < 0.0 {
                m = m.scale(-1.0);
            }
            let fresnel = interface.fresnel(wo.dot(m), eta_i, eta_t).unwrap();
            let survival = fresnel.reflectance + fresnel.transmittance;
            assert!(survival > 0.0 && survival < 0.99);
            let wi = sample.direction;
            let mask = 2.0
                / (1.0
                    + (1.0 + ((0.2 * wi.x).powi(2) + (0.6 * wi.y).powi(2)) / wi.z.powi(2)).sqrt());
            close(
                sample.radiance_weight * ratio * ratio,
                survival * mask,
                2.0e-11,
            );
        }
    }
}

#[test]
fn declared_axes_and_slope_admission_never_fall_back_to_isotropy() {
    for (ax, ay) in [(0.0, 0.2), (0.2, 0.0), (1.01, 0.2), (0.2, f64::NAN)] {
        assert_eq!(
            DielectricSurface::try_anisotropic_rough(ax, ay, x()),
            Err(DielectricError::InvalidRoughness)
        );
    }
    assert_eq!(
        DielectricSurface::try_anisotropic_rough(0.2, 0.6, Vec3::new(0.0, 0.0, 0.0)),
        Err(DielectricError::InvalidDirection)
    );
    let parallel = DielectricSurface::try_anisotropic_rough(0.2, 0.6, n()).unwrap();
    assert!(matches!(
        parallel.interface(550.0, false).unwrap().sample_rough(
            n(),
            n(),
            1.0,
            1.5,
            0.2,
            0.2,
            0.3,
            0.8
        ),
        Err(DielectricError::InvalidDirection)
    ));
    let surface = surface(0.2, 0.6);
    let interface = surface.interface(550.0, false).unwrap();
    assert!(matches!(
        interface.evaluate_rough(n(), n(), n(), 1.0, 1.5, 0.6),
        Err(DielectricError::InvalidRoughness)
    ));
    assert!(matches!(
        interface.sample_rough(n(), n(), 1.0, 1.5, 0.2, 1.0, 0.2, 0.8),
        Err(DielectricError::InvalidInterface)
    ));
    let null = interface
        .sample_rough(n(), n(), 1.5, 1.5, 0.2, 0.2, 0.3, 0.4)
        .unwrap()
        .unwrap();
    assert!(null.delta);
    assert_eq!(null.direction, n().scale(-1.0));
    assert_eq!(null.radiance_weight, 1.0);
}

#[test]
fn grazing_total_reflection_and_weak_index_contrasts_remain_finite() {
    for (ax, ay) in [(1.0e-4, 1.0), (1.0, 1.0e-4)] {
        let surface = surface(ax, ay);
        let interface = surface.interface(550.0, false).unwrap();
        for cosine in [1.0, 0.1, 1.0e-7] {
            let wo = Vec3::new((1.0_f64 - cosine * cosine).sqrt(), 0.0, cosine);
            for (eta_i, eta_t) in [(1.5, 1.0), (1.0, f64::from_bits(1.0_f64.to_bits() + 1))] {
                for u in [0.0, 0.2, 0.9] {
                    if let Some(sample) = interface
                        .sample_rough(n(), wo, eta_i, eta_t, ax, u, 0.3, 0.99)
                        .unwrap()
                    {
                        assert!(
                            [sample.value, sample.pdf, sample.radiance_weight]
                                .iter()
                                .all(|v| v.is_finite() && *v >= 0.0)
                        );
                        let reverse_normal = if sample.event == DielectricEvent::Transmission {
                            n().scale(-1.0)
                        } else {
                            n()
                        };
                        assert!(reverse_normal.dot(sample.direction) > 0.0);
                    }
                }
            }
        }
    }
}
