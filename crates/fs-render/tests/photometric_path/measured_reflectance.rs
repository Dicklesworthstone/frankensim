//! Independent spectral and energy references for passive measured diffuse walls.
//! The tables are declared reference materials, not measured-fixture certificates.

use super::*;
use fs_render::camera::{Aperture, CameraProjection, PhysicalCamera};
use fs_render::photometry::area::{EmissionSidedness, PhotometricAreaEmitter, PhotometricEmitters};
use fs_render::photometry::image::LuminanceJob;
use fs_render::photometry::medium::{HomogeneousMedium, ScatterSettings};
use fs_render::photometry::path::{IlluminanceSampling, PathIlluminanceJob};
use fs_render::photometry::{IlluminanceProbe, PhotometricScene};
use fs_render::spectral::SpectralReflectance;

use super::spectral_reference::cie_reference_xyz;

#[test]
fn g1_measured_diffuse_direct_image_matches_independent_photopic_reflection() {
    with_cx(|cx, _| {
        let mut bins = [0.0; SPECTRAL_BINS];
        bins[13] = 1.0; // 447.5 nm: equal radiant power, not equal lumens.
        bins[54] = 1.0; // 652.5 nm.
        let sources = [isotropic(
            [0.0, 0.0, 3.0],
            90.0,
            SourceSpectrum::try_from_bins(bins).unwrap(),
        )];
        let camera = PhysicalCamera::try_look_at(
            Point3::new(0.0, 0.0, 2.0),
            Point3::new(0.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            CameraProjection::try_vertical_fov(PI / 3.0).unwrap(),
            2.0,
            Aperture::try_circular(0.0).unwrap(),
        )
        .unwrap();
        let blue_y = cie_reference_xyz(447.5)[1];
        let red_y = cie_reference_xyz(652.5)[1];
        for (blue_rho, red_rho) in [(0.0, 0.0), (1.0, 1.0), (0.8, 0.05), (0.05, 0.8)] {
            let reflectance = if blue_rho == red_rho {
                SpectralReflectance::constant(blue_rho).unwrap()
            } else {
                SpectralReflectance::table(&[
                    (380.0, blue_rho),
                    (447.5, blue_rho),
                    (652.5, red_rho),
                    (780.0, red_rho),
                ])
                .unwrap()
            };
            let mut scene = receiver_scene();
            scene.primitives[0].material = Material::SpectralDiffuse { reflectance };
            let direct = PhotometricScene::try_new(
                &scene.primitives,
                &sources,
                HomogeneousMedium::vacuum(),
                TransportOptions::default(),
                cx,
            )
            .unwrap();
            let mut image =
                LuminanceJob::new(&direct, &camera, 1, 1, 8.0, ScatterSettings::default()).unwrap();
            assert_eq!(image.advance(cx, 1).unwrap(), 1);
            assert!(image.is_complete());
            let value = image.results()[0];
            // I/r² = 90 cd/(3 m)² = 10 lx. Independent original-CIE
            // weights and the declared table values give rho_eff*E/pi.
            // No production reflectance evaluator supplies the reference.
            let expected = (10.0 / PI) * (blue_y * blue_rho + red_y * red_rho) / (blue_y + red_y);
            close(
                value.surface_cd_m2,
                expected,
                3.0e-12,
                "measured spectral rho*E/pi",
            );
            close(
                value.luminance_cd_m2,
                expected,
                3.0e-12,
                "direct measured-wall luminance",
            );
            assert_eq!(value.single_scatter.luminance_cd_m2, 0.0);
            assert!(
                (0.0..=(10.0 / PI) * (1.0 + 3.0e-12)).contains(&value.surface_cd_m2),
                "a passive diffuse wall cannot return more than E/pi"
            );
            assert_eq!(value.numerical.kind, fs_evidence::NumericalKind::NoClaim);
            assert!(value.total_error_bound_cd_m2.is_none());

            // A measured passive wall supplies no LTE emission. This also
            // exercises the native path admission for exact black and white.
            let unlit = render_pixel(&scene, &[], PathMedium::Vacuum, settings(1, 4), cx);
            assert_eq!(unlit.xyz, [0.0; 3]);
        }
    });
}

#[test]
fn g1_measured_diffuse_interreflection_retains_each_wavelength() {
    with_cx(|cx, _| {
        let emitted = 8.0;
        let mut scene = empty_scene([0.0; 3], [0.0, 0.0, 1.0]);
        let mut enclosure = super::area::diffuse_emitting_cube(0.0);
        enclosure.material = Material::SpectralDiffuse {
            reflectance: SpectralReflectance::table(&[(380.0, 0.0), (780.0, 1.0)]).unwrap(),
        };
        scene.primitives.push(enclosure);
        // The retained original CIE data make every 5 nm bin carry equal
        // photopic source power. All 24 m² of inward-facing wall emits Le.
        assert_eq!(SPECTRAL_BINS, 80);
        let spectrum = SourceSpectrum::try_from_bins(core::array::from_fn(|bin| {
            1.0 / cie_reference_xyz(382.5 + 5.0 * bin as f64)[1]
        }))
        .unwrap();
        let surfaces = [PhotometricAreaEmitter::try_new(
            0,
            24.0 * PI * emitted,
            1.0,
            spectrum,
            EmissionSidedness::Front,
        )
        .unwrap()];
        let emitters = PhotometricEmitters {
            points: &[],
            surfaces: &surfaces,
        };
        // At the 80 represented midpoints, rho_b=(b+1/2)/80.
        // Every point sees a complete hemisphere of the same wall spectrum,
        // so depth D contributes Le*mean(sum_{k=0}^{D-1} rho_b^k).
        // mean(rho)=1/2, mean(rho²)=1/3-1/(12*80²). Averaging rho
        // before the second reflection would incorrectly substitute 1/4.
        let second_moment = 1.0 / 3.0 - 1.0 / (12.0 * 80.0 * 80.0);
        let three_orders = 1.0 + 0.5 + second_moment;
        for (depth, factor) in [(1, 1.0), (3, three_orders)] {
            let mut image = PathLuminanceJob::new_with_emitters(
                &scene,
                emitters,
                PathMedium::Vacuum,
                TransportOptions::default(),
                Settings {
                    strategy: DirectStrategy::BsdfOnly,
                    ..settings(depth, 8192)
                },
                cx,
            )
            .unwrap();
            image.advance(cx, usize::MAX).unwrap();
            let value = image.pixel_estimate(0).unwrap();
            close(
                value.luminance_cd_m2,
                emitted * factor,
                0.005,
                "wavelength-resolved diffuse reflection series",
            );
            assert_eq!(value.numerical.kind, fs_evidence::NumericalKind::NoClaim);
            assert!(value.total_error_bound_cd_m2.is_none());
        }

        // Probe MIS includes terminal-wall NEE, hence depth two includes
        // orders 0,1,2: the same three terms as the hit-only camera above.
        let probes =
            [
                IlluminanceProbe::try_new(Point3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 0.0, 1.0))
                    .unwrap(),
            ];
        let mut probe = PathIlluminanceJob::new_with_emitters(
            &scene,
            emitters,
            PathMedium::Vacuum,
            &probes,
            TransportOptions::default(),
            IlluminanceSampling {
                samples_per_probe: 65_536,
                max_depth: 2,
                sampler: Sampler::OwenSobol,
                seed: 173,
            },
            cx,
        )
        .unwrap();
        probe.advance(cx, usize::MAX).unwrap();
        let value = probe.probe_estimate(0).unwrap();
        close(
            value.illuminance_lux,
            PI * emitted * three_orders,
            0.01,
            "measured-wall cavity illuminance in lux",
        );
        assert_eq!(value.direct_lux, 0.0);
        assert_eq!(value.numerical.kind, fs_evidence::NumericalKind::NoClaim);
        assert!(value.total_error_bound_lux.is_none());
    });
}
