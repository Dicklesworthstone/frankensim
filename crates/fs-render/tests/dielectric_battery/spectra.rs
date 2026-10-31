//! G1/G3: source-shaped spectra traverse actual closed refractive geometry.
use super::*;
use fs_render::dielectric::{DielectricIor, DielectricIorSample, TabulatedIor};
use fs_render::spectral::xyz_of_spectrum;

fn sampled_material(indices: [f64; 9], absorption: [f64; 9]) -> Material {
    let absorption_samples: [(f64, f64); 9] =
        core::array::from_fn(|i| (380.0 + 50.0 * i as f64, absorption[i]));
    Material::Dielectric {
        glass: DielectricGlass::with_ior(
            DielectricIor::Tabulated(
                TabulatedIor::try_new(core::array::from_fn(|i| {
                    DielectricIorSample::try_new(380.0 + 50.0 * i as f64, indices[i]).unwrap()
                }))
                .unwrap(),
            ),
            BeerLambertAbsorption::try_from_spectral_extinction(&absorption_samples).unwrap(),
            GlassProvenance::Custom,
        ),
        surface: DielectricSurface::SMOOTH,
    }
}

// Independent regular-grid interpolant; no production optical evaluator is
// used by the radiance oracle below.
fn interpolated(values: [f64; 9], wavelength_nm: f64) -> f64 {
    let position = (wavelength_nm - 380.0) / 50.0;
    let index = (position.floor() as usize).min(7);
    let fraction = position - index as f64;
    values[index] * (1.0 - fraction) + values[index + 1] * fraction
}

#[test]
fn g1_sampled_dispersion_and_absorption_match_independent_slab_radiance() {
    // Deliberately nonmonotonic dispersion and an absorption peak. A Cauchy
    // replacement or RGB tint cannot represent these chosen source responses.
    let indices = [1.45, 1.49, 1.55, 1.70, 1.60, 1.53, 1.50, 1.48, 1.46];
    let absorption = [0.1, 0.2, 0.5, 1.4, 2.8, 1.2, 0.4, 0.2, 0.1];
    let glass = sampled_material(indices, absorption);
    let emission = lift_rgb([1.0; 3]);
    with_cx(|cx| {
        for cosine in [1.0_f64, 0.8] {
            let mut scene = scene_with_boundaries(&[(0.0, -0.2, glass, false)], -1.0, 4.0);
            scene.camera.forward = Vec3::new((1.0 - cosine * cosine).sqrt(), 0.0, -cosine);
            // Three vertices retain only direct transmission: front, back,
            // emitter. Internal round trips cannot reach the light in this
            // depth budget. Exit cancels entry's radiance eta-squared factor.
            let expected = xyz_of_spectrum(|wavelength| {
                let n = interpolated(indices, wavelength);
                let sigma = interpolated(absorption, wavelength);
                let cos_t = (1.0 - (1.0 - cosine * cosine) / (n * n)).sqrt();
                let rs = (cosine - n * cos_t) / (cosine + n * cos_t);
                let rp = (n * cosine - cos_t) / (n * cosine + cos_t);
                let reflection = 0.5 * (rs * rs + rp * rp);
                3.0 * emission.eval(wavelength)
                    * (1.0 - reflection).powi(2)
                    * (-sigma * 0.2 / cos_t).exp()
            });
            let mut settings = settings(DirectStrategy::Mis, 8192, 3);
            settings.sampler = Sampler::OwenSobol;
            let reference = render(&scene, cx, &settings).unwrap();
            for (channel, want) in expected.into_iter().enumerate() {
                let mean = reference.xyz[0][channel] / f64::from(reference.spp_done);
                assert!(
                    (mean / want - 1.0).abs() < 0.008,
                    "cosine={cosine}, channel={channel}: rendered {mean}, analytic {want}"
                );
            }
            for strategy in [
                DirectStrategy::BsdfOnly,
                DirectStrategy::NeeOnly,
                DirectStrategy::PowerMis,
            ] {
                settings.strategy = strategy;
                assert_film_bits_eq(
                    &render(&scene, cx, &settings).unwrap(),
                    &reference,
                    "sampled delta transmission cannot gain a light-sampling competitor",
                );
            }
        }
    });
}

#[test]
fn g3_constant_tables_match_analytic_transport_and_every_knot_rekeys_identity() {
    let sampled = sampled_material([1.5; 9], [0.7; 9]);
    let analytic = material(1.5, 0.7, DielectricSurface::SMOOTH);
    with_cx(|cx| {
        let settings = settings(DirectStrategy::Mis, 256, 5);
        let render_one = |material| {
            render(
                &scene_with_boundaries(&[(0.0, -0.2, material, false)], -1.0, 4.0),
                cx,
                &settings,
            )
            .unwrap()
        };
        assert_film_bits_eq(
            &render_one(sampled),
            &render_one(analytic),
            "constant tables preserve nondispersive analytic transport",
        );
    });
    // Model identity records the selected law, even at equivalent values.
    let identity = sampled.content_identity();
    assert_ne!(identity, analytic.content_identity());
    for i in 0..9 {
        let mut index = [1.5; 9];
        index[i] = 1.51;
        let mut absorption = [0.7; 9];
        absorption[i] = 0.71;
        assert_ne!(
            identity,
            sampled_material(index, [0.7; 9]).content_identity()
        );
        assert_ne!(
            identity,
            sampled_material([1.5; 9], absorption).content_identity()
        );
    }
}
