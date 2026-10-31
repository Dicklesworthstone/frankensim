//! CIE 171:2006 section 5.2, diffuse point-source analytical subcase only.
//!
//! The 3 m height and 0.5 m probe grid are publicly described on printed page
//! 17 of Dau Design & Consulting's 2007 AGi32 validation report:
//! <https://www.agi32.com/Downloads/TechnicalDocs/Report%20on%20AGI%2032%20validation%20of%20CIE%20171_Compiled_070620.pdf>.
//! Lumen Pro's validation description independently states I(theta) = 1000 cos(theta) cd:
//! <https://lumenpro.eu/cie-1712006-validation/>.
//! This is original code using those factual parameters and a derived formula;
//! no report text, diagrams or rounded reference tables are redistributed.
//! It tests neither the separate asymmetric photometry nor a file-format parser,
//! and establishes no full CIE 171 conformance or measured-fixture accuracy.

use super::*;
use fs_evidence::NumericalKind;
use fs_render::photometry::medium::HomogeneousMedium;
use fs_render::photometry::path::{IlluminanceSampling, PathIlluminanceJob};
use fs_render::photometry::{IlluminanceProbe, PhotometricScene};

#[test]
fn g1_cie171_5_2_diffuse_point_source_matches_derived_illuminance() {
    with_cx(|cx, _| {
        // Horizontal positions A--J relative to the source's vertical projection.
        let positions = [
            [0.0, 0.0],
            [0.5, 0.0],
            [1.0, 0.0],
            [1.5, 0.0],
            [0.5, 0.5],
            [1.0, 0.5],
            [1.5, 0.5],
            [1.0, 1.0],
            [1.5, 1.0],
            [1.5, 1.5],
        ];
        let probes = positions.map(|[x, y]| {
            IlluminanceProbe::try_new(Point3::new(x, y, 0.0), Vec3::new(0.0, 0.0, 1.0)).unwrap()
        });
        // Integral of I0*cos(theta) over the downward hemisphere is pi*I0 lm.
        let sources = [PhotometricSource::try_new(
            RigidTransform::try_new([1.0, 0.0, 0.0, 0.0], [0.0, 0.0, 3.0]).unwrap(),
            1000.0 * PI,
            1.0,
            AngularDistribution::cosine_power(1.0, FRAC_PI_2).unwrap(),
            SourceSpectrum::equal_energy(),
        )
        .unwrap()];
        // Both cosines are h/r: E = I0*h^2/r^4, in lux. Compute this directly,
        // without the production source evaluator or rounded tabulated values.
        let expected = positions.map(|[x, y]| {
            let distance_squared = 9.0 + x * x + y * y;
            9000.0 / (distance_squared * distance_squared)
        });
        let scene = empty_scene([0.0; 3], [0.0, 0.0, 1.0]);
        let direct = PhotometricScene::try_new(
            &scene.primitives,
            &sources,
            HomogeneousMedium::vacuum(),
            TransportOptions::default(),
            cx,
        )
        .unwrap();
        for (probe, reference) in probes.iter().zip(expected) {
            let value = direct.illuminance(cx, probe).unwrap();
            close(value.lux, reference, 2.0e-13, "CIE 171 5.2 direct lux");
            assert_eq!(value.numerical.kind, NumericalKind::NoClaim);
            assert!(value.total_error_bound_lux.is_none());
        }

        for sampler in [Sampler::Iid, Sampler::OwenSobol] {
            for samples_per_probe in [1, 4] {
                let mut job = PathIlluminanceJob::new(
                    &scene,
                    &sources,
                    PathMedium::Vacuum,
                    &probes,
                    TransportOptions::default(),
                    IlluminanceSampling {
                        samples_per_probe,
                        max_depth: 2,
                        sampler,
                        seed: 173,
                    },
                    cx,
                )
                .unwrap();
                job.advance(cx, usize::MAX).unwrap();
                assert!(job.is_complete());
                assert_eq!(job.probe_count(), probes.len());
                for (index, reference) in expected.into_iter().enumerate() {
                    let value = job.probe_estimate(index).unwrap();
                    close(
                        value.direct_lux,
                        reference,
                        2.0e-13,
                        "CIE 171 5.2 path direct lux",
                    );
                    close(
                        value.illuminance_lux,
                        reference,
                        2.0e-13,
                        "CIE 171 5.2 path total lux",
                    );
                    assert_eq!(value.continuous_lux, 0.0);
                    assert_eq!(value.samples, samples_per_probe);
                    assert_eq!(value.numerical.kind, NumericalKind::NoClaim);
                    assert!(value.total_error_bound_lux.is_none());
                    // The deterministic direct term has no sampling variance;
                    // zero IID standard error still is not a total error bound.
                    assert_eq!(
                        value.estimated_standard_error_lux,
                        if sampler == Sampler::Iid && samples_per_probe >= 2 {
                            Some(0.0)
                        } else {
                            None
                        }
                    );
                }
            }
        }
    });
}
