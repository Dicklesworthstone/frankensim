//! G1/G3 independent field and transport checks for transmitting lossy films.

use fs_geom::Vec3;
use fs_math::c64::C64;
use fs_render::conductor::{
    CauchyCoatingLayer, CoatingIorSample, CoatingLayer, CoherentAbsorbingCoating,
    CoherentDielectricCoating, CoherentPeriodicCoating, ConductorDataStatus, ConductorSource,
};
use fs_render::dielectric::{
    CauchyIor, DielectricCoating, DielectricCoatingStack, DielectricEvent, DielectricSurface,
};
use fs_render::polarization::Mueller;

fn source() -> ConductorSource {
    ConductorSource::try_new(
        fs_blake3::hash_domain("org.frankensim.test.absorbing-glass", b"declared-film-data"),
        ConductorDataStatus::Representative,
    ).unwrap()
}

fn layer(n: f64, k: f64, thickness_nm: f64) -> CauchyCoatingLayer {
    CauchyCoatingLayer::try_new(CauchyIor::try_constant(n).unwrap(), k, thickness_nm * 1.0e-9).unwrap()
}

fn surface(layers: &[CauchyCoatingLayer], repetitions: usize) -> DielectricSurface {
    DielectricSurface::SMOOTH.with_coating(DielectricCoating::try_new(
        CoherentPeriodicCoating::try_new(layers, repetitions, source()).unwrap(),
    ).unwrap())
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!((actual - expected).abs() <= tolerance,
        "actual={actual:.17e}, expected={expected:.17e}, tolerance={tolerance:.3e}");
}

fn terms(matrix: Mueller) -> [f64; 4] {
    [matrix.0[0][0], matrix.0[0][1], matrix.0[2][2], matrix.0[3][2]]
}

fn matrix_terms([s, p]: [C64; 2]) -> [f64; 4] {
    let cross = s * p.conj();
    [0.5 * (s.norm_sq() + p.norm_sq()), 0.5 * (s.norm_sq() - p.norm_sq()), cross.re, cross.im]
}

// Independent interface-scattering/Airy recursion, not the production
// projective characteristic-field propagation. Every exponential decays in
// its actual physical layer; only the two surrounding media are lossless.
fn airy(layers: &[(f64, f64, f64)], wavelength: f64, cosine: f64, eta_i: f64, eta_t: f64)
    -> ([f64; 4], [f64; 4])
{
    let transverse = eta_i * eta_i * (1.0 - cosine * cosine);
    let mut indices = vec![C64::from_re(eta_i)];
    indices.extend(layers.iter().map(|&(n, k, _)| C64::new(n, k)));
    indices.push(C64::from_re(eta_t));
    let q: Vec<_> = indices.iter().map(|&n| (n * n - C64::from_re(transverse)).sqrt()).collect();
    let mut reflected = [C64::ZERO; 2];
    let mut transmitted = [C64::ZERO; 2];
    for pol in 0..2 {
        let y: Vec<_> = indices.iter().zip(&q).map(|(&n, &q)| if pol == 0 { q } else { q / (n*n) }).collect();
        let last = layers.len();
        let mut r = (y[last] - y[last+1]) / (y[last] + y[last+1]);
        let mut t = y[last].scale(2.0) / (y[last] + y[last+1]);
        for j in (0..last).rev() {
            let delta = q[j+1].scale(core::f64::consts::TAU * layers[j].2 / wavelength);
            let attenuation = (-delta.im).exp();
            let phase = C64::new(delta.re.cos() * attenuation, delta.re.sin() * attenuation);
            let interface_r = (y[j] - y[j+1]) / (y[j] + y[j+1]);
            let denominator = C64::ONE + interface_r * r * phase * phase;
            t = (y[j].scale(2.0) / (y[j] + y[j+1])) * t * phase / denominator;
            r = (interface_r + r * phase * phase) / denominator;
        }
        reflected[pol] = r;
        transmitted[pol] = if y[last+1].re > 0.0 {
            t.scale((y[last+1].re / y[0].re).sqrt())
        } else { C64::ZERO };
    }
    (matrix_terms(reflected), matrix_terms(transmitted))
}

#[test]
fn g1_absorbing_multilayers_match_independent_airy_fields() {
    for k in [0.0, 1.0e-6, 0.03, 0.7, 8.0] {
        let laws = [(1.35, 0.0, 65.0), (2.1, k, 43.0), (1.6, 0.2*k, 97.0)];
        let compact = laws.map(|(n, k, d)| layer(n, k, d));
        let coated = surface(&compact, 1);
        for wavelength in [380.0, 455.0, 550.0, 685.0, 780.0] {
            for cosine in [0.03, 0.3, 0.7, 1.0] {
                for (eta_i, eta_t) in [(1.0, 1.7), (1.6, 1.0)] {
                    let response = coated.interface(wavelength, false).unwrap()
                        .fresnel_polarized(cosine, eta_i, eta_t).unwrap();
                    let expected = airy(&laws, wavelength, cosine, eta_i, eta_t);
                    for (actual, target) in terms(response.reflection).into_iter().zip(expected.0) {
                        near(actual, target, 3.0e-12);
                    }
                    for (actual, target) in terms(response.transmission).into_iter().zip(expected.1) {
                        near(actual, target, 3.0e-12);
                    }
                    assert!(response.scalar.absorption_residual() >= -4.0e-15);
                    if k > 0.0 { assert!(response.scalar.absorption_residual() > 0.0); }
                }
            }
        }
    }
}

#[test]
fn g3_absorbing_stack_reversal_keeps_reciprocal_transmission() {
    let coated = surface(&[layer(1.35, 0.02, 90.0), layer(2.1, 0.6, 33.0)], 3);
    for cosine in [0.2_f64, 0.6, 1.0] {
        let forward = coated.interface(540.0, false).unwrap()
            .fresnel_polarized(cosine, 1.0, 1.7).unwrap();
        let exit_cosine = (1.0 - (1.0-cosine*cosine)/1.7_f64.powi(2)).sqrt();
        let reverse = coated.interface(540.0, true).unwrap()
            .fresnel_polarized(exit_cosine, 1.7, 1.0).unwrap();
        for (a,b) in terms(forward.transmission).into_iter().zip(terms(reverse.transmission)) {
            near(a, b, 2.0e-12);
        }
        assert!((forward.scalar.reflectance - reverse.scalar.reflectance).abs() > 1.0e-3);
        assert!(forward.scalar.absorption_residual() > 0.1);
        assert!(reverse.scalar.absorption_residual() > 0.1);
    }
}

#[test]
fn g1_absorbing_equal_index_and_tir_interfaces_do_not_become_lossless() {
    let coated = surface(&[layer(1.5, 0.12, 300.0)], 1);
    let normal = Vec3::new(0.0, 0.0, 1.0);
    let same = coated.interface(550.0, false).unwrap().fresnel(1.0, 1.5, 1.5).unwrap();
    assert!(same.reflectance > 0.0 && same.transmittance < 0.6);
    let sample = coated.interface(550.0, false).unwrap().sample_smooth(normal, normal, 1.5, 1.5, 0.9).unwrap();
    assert_eq!(sample.event, DielectricEvent::Transmission);
    near(sample.radiance_weight, same.reflectance + same.transmittance, 3.0e-15);
    let tir = coated.interface(550.0, true).unwrap().fresnel(0.5, 1.6, 1.0).unwrap();
    assert!(tir.total_internal_reflection);
    assert_eq!(tir.transmittance, 0.0);
    assert!(tir.reflectance > 0.0 && tir.reflectance < 0.95);
    let wo = Vec3::new(0.75_f64.sqrt(), 0.0, 0.5);
    let sample = coated.interface(550.0, true).unwrap().sample_smooth(normal, wo, 1.6, 1.0, 0.7).unwrap();
    assert_eq!(sample.event, DielectricEvent::Reflection);
    assert_eq!(sample.probability, 1.0);
    near(sample.radiance_weight, tir.reflectance, 3.0e-15);
}

#[test]
fn g1_lossy_smooth_and_ggx_sampling_preserve_survival_power() {
    let coated = surface(&[layer(1.8, 0.25, 130.0)], 1);
    let interface = coated.interface(550.0, false).unwrap();
    let expected = interface.fresnel_polarized(1.0, 1.0, 1.5).unwrap();
    let survival = expected.scalar.reflectance + expected.scalar.transmittance;
    assert!(survival < 0.8);
    let normal = Vec3::new(0.0, 0.0, 1.0);
    let mut power = 0.0;
    for u in [0.0, 1.0_f64.next_down()] {
        let sample = interface.sample_smooth(normal, normal, 1.0, 1.5, u).unwrap();
        let jacobian = if sample.event == DielectricEvent::Transmission { 1.5_f64.powi(2) } else { 1.0 };
        near(sample.radiance_weight * jacobian, survival, 3.0e-15);
        power += sample.probability * sample.radiance_weight * jacobian;
        let polarized = interface.sample_smooth_polarized(normal, normal, 1.0, 1.5, u).unwrap();
        near(polarized.radiance_weight.0[0][0], sample.radiance_weight, 3.0e-15);
        let rough = interface.sample_rough(normal, normal, 1.0, 1.5, 0.2, 0.0, 0.0, u).unwrap().unwrap();
        assert_eq!(rough.event, sample.event);
        near(rough.radiance_weight * jacobian, survival, 5.0e-15);
        let evaluation = interface.evaluate_rough(normal, normal, rough.direction, 1.0, 1.5, 0.2).unwrap();
        near(rough.pdf, evaluation.pdf, 0.0);
        near(rough.value, evaluation.value, 0.0);
    }
    near(power, survival, 3.0e-15);
    near(power + expected.scalar.absorption_residual(), 1.0, 3.0e-15);
}

#[test]
fn g1_weak_loss_opaque_limit_and_zero_thickness_keep_physical_channels() {
    for k in [1.0e-6, 1.0e-10, 1.0e-12] {
        let coated = surface(&[layer(1.5, k, 100.0)], 1);
        let result = coated.interface(550.0, false).unwrap().fresnel(1.0, 1.5, 1.5).unwrap();
        let expected = airy(&[(1.5,k,100.0)], 550.0, 1.0, 1.5, 1.5);
        near(result.reflectance, expected.0[0], 1.0e-24);
        near(result.transmittance, expected.1[0], 5.0e-15);
        assert!(result.absorption_residual() > 0.0);
    }
    let opaque = surface(&[layer(1.8, 2.0, 100_000.0)], 16);
    let result = opaque.interface(550.0, false).unwrap().fresnel(0.8, 1.0, 1.5).unwrap();
    assert_eq!(result.transmittance, 0.0);
    let (reference, _) = airy(&[(1.8,2.0,100_000.0)], 550.0, 0.8, 1.0, 1.5);
    near(result.reflectance, reference[0], 2.0e-14);
    assert!(result.absorption_residual() > 0.1);
    let absent = surface(&[layer(1.8, 2.0, 0.0)], 1);
    assert_eq!(absent.interface(550.0, false).unwrap().fresnel(0.8, 1.0, 1.5).unwrap(),
        fs_render::dielectric::fresnel_dielectric(0.8, 1.0, 1.5).unwrap());
}

#[test]
fn g3_original_spectral_mixed_stack_preserves_tables_and_work_bound() {
    let transparent = CoherentDielectricCoating::try_new(CauchyIor::try_constant(1.3).unwrap(), 70e-9, source()).unwrap();
    let samples = core::array::from_fn(|i| CoatingIorSample::try_new(
        380.0+50.0*i as f64, 1.7+0.025*i as f64, 0.01+0.035*i as f64).unwrap());
    let absorbing = CoherentAbsorbingCoating::try_new(samples, 80e-9, source()).unwrap();
    let layers = [CoatingLayer::Dielectric(&transparent), CoatingLayer::Absorbing(&absorbing)];
    let stack = DielectricCoatingStack::try_mixed(&layers).unwrap();
    for sample in samples {
        let wavelength = sample.wavelength_nm();
        let result = stack.fresnel_polarized(wavelength, false, 0.65, 1.0, 1.6).unwrap();
        let expected = airy(&[(1.3,0.0,70.0),(sample.eta(),sample.k(),80.0)], wavelength, 0.65, 1.0, 1.6);
        for (a,b) in terms(result.reflection).into_iter().zip(expected.0) { near(a,b,2.0e-13); }
        for (a,b) in terms(result.transmission).into_iter().zip(expected.1) { near(a,b,2.0e-13); }
    }
    assert_eq!(absorbing.samples(), samples);
    assert_eq!(stack.layer_count(), 2);
    assert!(DielectricCoatingStack::try_mixed(&[layers[0];65]).is_err());
    assert!(stack.fresnel(379.0, false, 0.8, 1.0, 1.5).is_err());
    assert!(stack.fresnel(550.0, false, 0.0, 1.0, 1.5).is_err());
    assert!(CoatingIorSample::try_new(550.0,1.5,-1e-9).is_err());
    assert!(CauchyCoatingLayer::try_new(CauchyIor::try_constant(1.5).unwrap(),-1e-9,80e-9).is_err());
}
