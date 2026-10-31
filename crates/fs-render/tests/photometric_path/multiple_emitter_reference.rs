//! Independent colored multi-emitter illumination through source-choice MIS.
//! Solid angles and original CIE data supply the reference, not another render.

use super::*;
use fs_render::photometry::IlluminanceProbe;
use fs_render::photometry::area::{EmissionSidedness, PhotometricAreaEmitter, PhotometricEmitters};
use fs_render::photometry::path::{IlluminanceSampling, PathIlluminanceJob};
use fs_render::spectral::SpectralReflectance;

use super::area::projected_polygon_lux;
use super::spectral_reference::cie_reference_xyz;

fn statistical_match(actual: f64, expected: f64, estimated_se: f64, context: &str) {
    assert!(estimated_se.is_finite() && estimated_se >= 0.0);
    // A fixed regression rule using estimated sampling variation, not a
    // confidence interval, numerical enclosure or physical accuracy claim.
    let tolerance = 6.0 * estimated_se + 2.0e-11 * expected.abs().max(1.0);
    assert!(
        actual.is_finite() && (actual - expected).abs() <= tolerance,
        "{context}: observed={actual}, expected={expected}, estimated_se={estimated_se}"
    );
}

#[test]
fn g1_colored_finite_emitters_match_independent_reflected_xyz_and_lux() {
    with_cx(|cx, _| {
        let mut scene = empty_scene([0.0, 0.0, 0.5], [0.0, 0.0, -1.0]);
        let mut receiver = floor(0.0, 10.0, LiftedSpectrum { c: [0.0; 3] });
        receiver.material = Material::SpectralDiffuse {
            reflectance: SpectralReflectance::table(&[
                (380.0, 0.65),
                (540.0, 0.65),
                (620.0, 0.2),
                (780.0, 0.2),
            ])
            .unwrap(),
        };
        scene.primitives.push(receiver);

        // Corners, independently known m², full-output cd/m², dimmer,
        // first spectral bin, and the declared wall rho over that band.
        // Positive polygon winding supplies the solid-angle reference;
        // the actual emitter triangles face downward toward the receiver.
        let panels = [
            (
                [
                    [-2.5, -1.2, 2.0],
                    [-0.1, -1.2, 2.0],
                    [-0.1, 1.2, 2.0],
                    [-2.5, 1.2, 2.0],
                ],
                5.76,
                80.0,
                0.5,
                12,
                0.65,
            ),
            (
                [
                    [0.2, -0.8, 1.5],
                    [2.2, -0.8, 1.5],
                    [2.2, 0.8, 1.5],
                    [0.2, 0.8, 1.5],
                ],
                3.2,
                160.0,
                0.75,
                48,
                0.2,
            ),
        ];
        let probe_geometry = [
            ([0.0; 3], [0.0, 0.0, 1.0]),
            ([0.35, 0.15, 0.0], [0.2, 0.0, 0.96_f64.sqrt()]),
        ];
        let probes = probe_geometry.map(|(p, n)| {
            IlluminanceProbe::try_new(Point3::new(p[0], p[1], p[2]), Vec3::new(n[0], n[1], n[2]))
                .unwrap()
        });
        let mut emitters = Vec::new();
        let mut expected_xyz = [0.0; 3];
        let mut expected_lux = [0.0; 2];
        for (index, (corners, area, full_luminance, dimmer, first_bin, rho)) in
            panels.into_iter().enumerate()
        {
            let mut bins = [0.0; SPECTRAL_BINS];
            let mut color = [0.0; 3];
            // Disjoint 100 nm bands: 440--540 and 620--720 nm. Equal
            // photopic power per bin makes each four-lane calibrated packet
            // carry exactly one active bin per source, avoiding rare-band
            // luminance noise while retaining very different CIE colors.
            for (bin, power) in bins.iter_mut().enumerate().skip(first_bin).take(20) {
                let observer = cie_reference_xyz(382.5 + 5.0 * bin as f64);
                *power = 1.0 / observer[1];
                for channel in 0..3 {
                    color[channel] += observer[channel] / observer[1] / 20.0;
                }
            }
            scene.primitives.push(Primitive {
                shape: Shape::Mesh(TriMesh::new(corners.to_vec(), vec![[0, 2, 1], [0, 3, 2]])),
                // Black sources make one wall reflection the complete
                // solution: no untraced interreflection tail is being fitted.
                material: Material::SpectralDiffuse {
                    reflectance: SpectralReflectance::constant(0.0).unwrap(),
                },
                emission: None,
            });
            emitters.push(
                PhotometricAreaEmitter::try_new(
                    index + 1,
                    PI * area * full_luminance,
                    dimmer,
                    SourceSpectrum::try_from_bins(bins).unwrap(),
                    EmissionSidedness::Front,
                )
                .unwrap(),
            );
            let luminance = full_luminance * dimmer;
            let incident = projected_polygon_lux(luminance, corners, [0.0; 3], [0.0, 0.0, 1.0]);
            // Le=Phi/(pi*A), E=Le*projected_solid_angle, Lr=rho*E/pi.
            // Original CIE ratios convert each source's photopic contribution
            // into XYZ; no production spectrum or reflectance evaluator is
            // used in this reference. Emission and reflection remain colored.
            for channel in 0..3 {
                expected_xyz[channel] += incident * rho / PI * color[channel];
            }
            for (expected, (p, n)) in expected_lux.iter_mut().zip(probe_geometry) {
                *expected += projected_polygon_lux(luminance, corners, p, n);
            }
        }

        for reverse in [false, true] {
            if reverse {
                // Reorder bindings only, retaining their primitive identities.
                // Different CDF addresses must still have the same expectation.
                emitters.reverse();
            }
            let sources = PhotometricEmitters {
                points: &[],
                surfaces: &emitters,
            };
            for strategy in [
                DirectStrategy::NeeOnly,
                DirectStrategy::Mis,
                DirectStrategy::PowerMis,
                DirectStrategy::BsdfOnly,
            ] {
                let mut image = PathLuminanceJob::new_with_emitters(
                    &scene,
                    sources,
                    PathMedium::Vacuum,
                    TransportOptions::default(),
                    Settings {
                        width: 32,
                        sampler: Sampler::Iid,
                        strategy,
                        ..settings(2, 1024)
                    },
                    cx,
                )
                .unwrap();
                image.advance(cx, usize::MAX).unwrap();
                assert!(image.is_complete());
                // half_tan=0 makes all pixels observe exactly the same point.
                // Independent IID pixel streams give replicated means, from
                // which sampling SE is estimated separately for X, Y and Z.
                let pixels: Vec<_> = image
                    .estimates()
                    .map(|value| {
                        let value = value.unwrap();
                        assert_eq!(value.samples, 1024);
                        assert_eq!(value.numerical.kind, fs_evidence::NumericalKind::NoClaim);
                        assert!(value.total_error_bound_cd_m2.is_none());
                        value.xyz
                    })
                    .collect();
                let n = pixels.len() as f64;
                for channel in 0..3 {
                    let mean = pixels.iter().map(|xyz| xyz[channel]).sum::<f64>() / n;
                    let se = (pixels
                        .iter()
                        .map(|xyz| (xyz[channel] - mean).powi(2))
                        .sum::<f64>()
                        / (n * (n - 1.0)))
                        .sqrt();
                    statistical_match(
                        mean,
                        expected_xyz[channel],
                        se,
                        &format!("{strategy:?}, reversed={reverse}, XYZ channel {channel}"),
                    );
                    if channel == 1 {
                        close(mean, expected_xyz[1], 0.03, "absolute reflected luminance");
                    }
                }
            }
            let mut job = PathIlluminanceJob::new_with_emitters(
                &scene,
                sources,
                PathMedium::Vacuum,
                &probes,
                TransportOptions::default(),
                IlluminanceSampling {
                    samples_per_probe: 16_384,
                    max_depth: 1,
                    sampler: Sampler::Iid,
                    seed: 173,
                },
                cx,
            )
            .unwrap();
            job.advance(cx, usize::MAX).unwrap();
            assert!(job.is_complete());
            for (value, expected) in job.estimates().zip(expected_lux) {
                let value = value.unwrap();
                assert_eq!(value.direct_lux, 0.0);
                statistical_match(
                    value.illuminance_lux,
                    expected,
                    value.estimated_standard_error_lux.unwrap(),
                    "incident lux from unequal colored finite emitters",
                );
                close(
                    value.illuminance_lux,
                    expected,
                    0.03,
                    "absolute incident lux",
                );
                assert_eq!(value.numerical.kind, fs_evidence::NumericalKind::NoClaim);
                assert!(value.total_error_bound_lux.is_none());
            }
        }
    });
}
