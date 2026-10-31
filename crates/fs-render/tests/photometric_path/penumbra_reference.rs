//! Finite-source penumbra from an intervening opaque mesh, with a geometric oracle.
//! A clipped emitting polygon supplies the reference without production visibility.

use super::*;
use fs_render::photometry::IlluminanceProbe;
use fs_render::photometry::area::{EmissionSidedness, PhotometricAreaEmitter, PhotometricEmitters};
use fs_render::photometry::path::{IlluminanceSampling, PathIlluminanceJob};
use fs_render::spectral::SpectralReflectance;

const LUMINANCE: f64 = 50.0;
const RHO: f64 = 0.4;

fn black_quad(corners: [[f64; 3]; 4]) -> Primitive {
    Primitive {
        shape: Shape::Mesh(TriMesh::new(corners.to_vec(), vec![[0, 2, 1], [0, 3, 2]])),
        material: Material::SpectralDiffuse {
            reflectance: SpectralReflectance::constant(0.0).unwrap(),
        },
        emission: None,
    }
}

fn reference_lux(x: f64) -> f64 {
    // Receiver (x,0,0), source point (u,v,2): at the z=1 blocker plane the
    // connection crosses ((x+u)/2,v/2,1). The blocker covers x<=0 across
    // every source connection in this test, so visibility is exactly u>-x.
    // Integrate the remaining rectangle; a visible-area fraction alone is
    // insufficient because source distance and projected angle vary over it.
    let left = (-x).max(-1.0);
    if left >= 1.0 {
        return 0.0;
    }
    super::area::projected_polygon_lux(
        LUMINANCE,
        [
            [left, -1.0, 2.0],
            [1.0, -1.0, 2.0],
            [1.0, 1.0, 2.0],
            [left, 1.0, 2.0],
        ],
        [x, 0.0, 0.0],
        [0.0, 0.0, 1.0],
    )
}

fn check_estimate(actual: f64, expected: f64, estimated_se: f64, context: &str) {
    if expected == 0.0 {
        assert_eq!(actual, 0.0, "{context}: complete geometric shadow");
    } else {
        // This fixed acceptance threshold cannot expand with estimated noise.
        close(actual, expected, 0.05, context);
    }
    // Supplemental IID consistency check, not a confidence interval or any
    // bound on numerical, fixture-calibration or physical-model uncertainty.
    assert!(estimated_se.is_finite() && estimated_se >= 0.0);
    assert!(
        (actual - expected).abs() <= 6.0 * estimated_se + 1.0e-10 * expected.max(1.0),
        "{context}: observed={actual}, expected={expected}, estimated_se={estimated_se}"
    );
}

#[test]
fn g1_opaque_mesh_penumbra_matches_clipped_source_lux_and_reflection() {
    with_cx(|cx, _| {
        let mut scene = empty_scene([0.0, 0.0, 0.5], [0.0, 0.0, -1.0]);
        let mut receiver = floor(0.0, 10.0, LiftedSpectrum { c: [0.0; 3] });
        receiver.material = Material::SpectralDiffuse {
            reflectance: SpectralReflectance::constant(RHO).unwrap(),
        };
        scene.primitives.push(receiver);
        scene.primitives.push(black_quad([
            [-1.0, -1.0, 2.0],
            [1.0, -1.0, 2.0],
            [1.0, 1.0, 2.0],
            [-1.0, 1.0, 2.0],
        ]));
        scene.primitives.push(black_quad([
            [-4.0, -2.0, 1.0],
            [0.0, -2.0, 1.0],
            [0.0, 2.0, 1.0],
            [-4.0, 2.0, 1.0],
        ]));
        // Flatten photopic bin power using original CIE data. This removes
        // spectral packet noise so the test resolves source/shadow sampling.
        let spectrum = SourceSpectrum::try_from_bins(core::array::from_fn(|bin| {
            1.0 / super::spectral_reference::cie_reference_xyz(382.5 + 5.0 * bin as f64)[1]
        }))
        .unwrap();
        let surfaces = [PhotometricAreaEmitter::try_new(
            1,
            4.0 * PI * LUMINANCE,
            1.0,
            spectrum,
            EmissionSidedness::Front,
        )
        .unwrap()];
        let emitters = PhotometricEmitters {
            points: &[],
            surfaces: &surfaces,
        };
        let positions = [-1.25, -0.5, 0.0, 0.5, 1.25];
        let probes = positions.map(|x| {
            IlluminanceProbe::try_new(Point3::new(x, 0.0, 0.0), Vec3::new(0.0, 0.0, 1.0)).unwrap()
        });
        let mut probes_job = PathIlluminanceJob::new_with_emitters(
            &scene,
            emitters,
            PathMedium::Vacuum,
            &probes,
            TransportOptions::default(),
            IlluminanceSampling {
                samples_per_probe: 32_768,
                max_depth: 1,
                sampler: Sampler::Iid,
                seed: 173,
            },
            cx,
        )
        .unwrap();
        probes_job.advance(cx, usize::MAX).unwrap();
        assert!(probes_job.is_complete());
        for (x, estimate) in positions.into_iter().zip(probes_job.estimates()) {
            let estimate = estimate.unwrap();
            assert_eq!(estimate.direct_lux, 0.0);
            assert_eq!(estimate.samples, 32_768);
            check_estimate(
                estimate.illuminance_lux,
                reference_lux(x),
                estimate.estimated_standard_error_lux.unwrap(),
                &format!("geometric penumbra probe at x={x}"),
            );
            assert_eq!(estimate.numerical.kind, fs_evidence::NumericalKind::NoClaim);
            assert!(estimate.total_error_bound_lux.is_none());
        }
        drop(probes_job);

        // Black source/blocker materials and a planar receiver leave no
        // higher-order reflection contribution. Depth two therefore includes
        // the complete rho*E/pi solution, including hit-only source discovery.
        // Compare every strategy in the penumbra, and MIS in dark/clear regions.
        for x in [-1.25, 0.0, 1.25] {
            scene.camera.eye.x = x;
            for strategy in [
                DirectStrategy::NeeOnly,
                DirectStrategy::Mis,
                DirectStrategy::PowerMis,
                DirectStrategy::BsdfOnly,
            ] {
                if x != 0.0 && strategy != DirectStrategy::Mis {
                    continue;
                }
                let samples = if x < -1.0 {
                    64
                } else if strategy == DirectStrategy::BsdfOnly {
                    65_536
                } else {
                    32_768
                };
                let mut image = PathLuminanceJob::new_with_emitters(
                    &scene,
                    emitters,
                    PathMedium::Vacuum,
                    TransportOptions::default(),
                    Settings {
                        strategy,
                        sampler: Sampler::Iid,
                        ..settings(2, samples)
                    },
                    cx,
                )
                .unwrap();
                image.advance(cx, usize::MAX).unwrap();
                assert!(image.is_complete());
                let estimate = image.pixel_estimate(0).unwrap();
                assert_eq!(estimate.samples, samples);
                check_estimate(
                    estimate.luminance_cd_m2,
                    RHO * reference_lux(x) / PI,
                    estimate.estimated_standard_error_cd_m2.unwrap(),
                    &format!("{strategy:?} reflected penumbra at x={x}"),
                );
                assert_eq!(estimate.numerical.kind, fs_evidence::NumericalKind::NoClaim);
                assert!(estimate.total_error_bound_cd_m2.is_none());
            }
        }
    });
}
