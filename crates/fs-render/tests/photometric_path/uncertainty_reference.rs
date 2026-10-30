//! Independent known-population checks of the reported IID sampling error.
//! These fixed regression tolerances are not confidence or physical-error bounds.

use super::*;
use fs_render::photometry::area::{EmissionSidedness, PhotometricAreaEmitter, PhotometricEmitters};

const PANEL_LUMINANCE_CD_M2: f64 = 40.0;

fn equal_photopic_bins() -> SourceSpectrum {
    // The original CIE CSV, not the production observer, supplies V(lambda).
    // S_b proportional to 1/V_b makes each represented spectral bin carry
    // equal photopic power. A hit therefore returns the same Y for every
    // sampled spectral packet; only its geometric hit/miss remains random.
    assert_eq!(SPECTRAL_BINS, 80);
    SourceSpectrum::try_from_bins(core::array::from_fn(|bin| {
        1.0 / super::spectral_reference::cie_reference_xyz(382.5 + 5.0 * bin as f64)[1]
    }))
    .unwrap()
}

fn quarter_pixel_panels(two_emitters: bool) -> (Scene, Vec<PhotometricAreaEmitter>) {
    let mut scene = empty_scene([0.0; 3], [0.0, 0.0, 1.0]);
    scene.camera.half_tan = 0.5;
    // At z=2 m this square pinhole pixel covers [-1,1] x [-1,1] m.
    // Uniform raster jitter makes each 0.5 m x 2 m strip occupy exactly
    // one quarter of its samples, irrespective of ray-length normalization.
    // The two strips are disjoint and both face the camera (negative Z).
    let spectrum = equal_photopic_bins();
    let mut emitters = Vec::new();
    let strips = [(-1.0, -0.5, 1.0), (0.5, 1.0, 3.0)];
    for &(left, right, multiplier) in strips.iter().take(if two_emitters { 2 } else { 1 }) {
        let index = scene.primitives.len();
        scene.primitives.push(Primitive {
            shape: Shape::Mesh(TriMesh::new(
                vec![
                    [left, -1.0, 2.0],
                    [right, -1.0, 2.0],
                    [right, 1.0, 2.0],
                    [left, 1.0, 2.0],
                ],
                vec![[0, 2, 1], [0, 3, 2]],
            )),
            material: Material::Lambertian {
                reflectance: LiftedSpectrum { c: [0.0; 3] },
            },
            emission: None,
        });
        // Each panel has area 1 m^2. A front-emitting Lambertian surface
        // with luminance L emits pi*A*L lumens.
        emitters.push(
            PhotometricAreaEmitter::try_new(
                index,
                PI * multiplier * PANEL_LUMINANCE_CD_M2,
                1.0,
                spectrum.clone(),
                EmissionSidedness::Front,
            )
            .unwrap(),
        );
    }
    (scene, emitters)
}

#[test]
fn g1_iid_luminance_error_matches_independent_population_moments() {
    with_cx(|cx, _| {
        for two_emitters in [false, true] {
            let (scene, surfaces) = quarter_pixel_panels(two_emitters);
            let emitters = PhotometricEmitters {
                points: &[],
                surfaces: &surfaces,
            };
            let policy = Settings {
                sampler: Sampler::Iid,
                strategy: DirectStrategy::BsdfOnly,
                ..settings(1, 16_384)
            };
            let mut job = PathLuminanceJob::new_with_emitters(
                &scene,
                emitters,
                PathMedium::Vacuum,
                TransportOptions::default(),
                policy,
                cx,
            )
            .unwrap();
            assert!(job.pixel_estimate(0).is_none());
            assert_eq!(job.advance(cx, 1).unwrap(), 1);
            assert!(
                job.pixel_estimate(0)
                    .unwrap()
                    .estimated_standard_error_cd_m2
                    .is_none()
            );

            let l = PANEL_LUMINANCE_CD_M2;
            // One strip: X is {0,L} with probabilities {3/4,1/4}.
            // Two strips: X is {0,L,3L} with probabilities {1/2,1/4,1/4}.
            // Hence the population moments below use neither sampled values
            // nor the renderer's own estimated standard error as a reference.
            let (mean, variance) = if two_emitters {
                (l, 1.5 * l * l)
            } else {
                (l / 4.0, 3.0 * l * l / 16.0)
            };
            let mut reported_errors = Vec::new();
            for samples in [4_096_u32, 16_384] {
                let remaining = u64::from(samples) - job.completed_path_samples();
                assert_eq!(
                    job.advance(cx, remaining as usize).unwrap(),
                    remaining as usize
                );
                let measured = job.pixel_estimate(0).unwrap();
                assert_eq!(measured.samples, samples);
                let population_error = (variance / f64::from(samples)).sqrt();
                assert!(
                    (measured.luminance_cd_m2 - mean).abs() <= 6.0 * population_error,
                    "two_emitters={two_emitters}, N={samples}: mean={}, expected={mean}",
                    measured.luminance_cd_m2,
                );
                let reported = measured.estimated_standard_error_cd_m2.unwrap();
                close(
                    reported,
                    population_error,
                    0.06,
                    "known-population IID standard error",
                );
                assert_eq!(measured.numerical.kind, fs_evidence::NumericalKind::NoClaim);
                assert!(measured.total_error_bound_cd_m2.is_none());
                reported_errors.push(reported);
            }
            assert!(job.is_complete());
            close(
                reported_errors[1] / reported_errors[0],
                0.5,
                0.1,
                "four times the samples halves IID standard error",
            );

            // The same finite-source scene under correlated Sobol sampling
            // must not inherit the IID standard-error interpretation.
            let mut sobol = PathLuminanceJob::new_with_emitters(
                &scene,
                emitters,
                PathMedium::Vacuum,
                TransportOptions::default(),
                Settings {
                    sampler: Sampler::OwenSobol,
                    spp: 256,
                    ..policy
                },
                cx,
            )
            .unwrap();
            sobol.advance(cx, usize::MAX).unwrap();
            let measured = sobol.pixel_estimate(0).unwrap();
            assert_eq!(measured.samples, 256);
            assert!(measured.estimated_standard_error_cd_m2.is_none());
            assert_eq!(measured.numerical.kind, fs_evidence::NumericalKind::NoClaim);
            assert!(measured.total_error_bound_cd_m2.is_none());
        }
    });
}
