//! Point-fixture energy injected into a closed, nonuniformly lit diffuse cavity.
//! The whole-wall flux oracle needs neither a radiosity solve nor emissive walls.

use super::*;
use fs_render::photometry::IlluminanceProbe;
use fs_render::photometry::path::{IlluminanceSampling, PathIlluminanceJob};
use fs_render::spectral::SpectralReflectance;

// Positions are strictly interior; their offsets make the wall field nonuniform.
const FIXTURES: [([f64; 3], f64); 2] = [([0.2, -0.15, 0.25], 600.0), ([-0.3, 0.1, -0.2], 400.0)];

// Gauss--Legendre nodes and weights on [-1,1], tensorized on each 2x2 m wall.
// These are geometric integration weights in square metres, not ray PDFs.
const GAUSS_6: [(f64, f64); 6] = [
    (-0.932_469_514_203_152_1, 0.171_324_492_379_170_4),
    (-0.661_209_386_466_264_5, 0.360_761_573_048_138_6),
    (-0.238_619_186_083_196_9, 0.467_913_934_572_691),
    (0.238_619_186_083_196_9, 0.467_913_934_572_691),
    (0.661_209_386_466_264_5, 0.360_761_573_048_138_6),
    (0.932_469_514_203_152_1, 0.171_324_492_379_170_4),
];
const GAUSS_8: [(f64, f64); 8] = [
    (-0.960_289_856_497_536_3, 0.101_228_536_290_376_3),
    (-0.796_666_477_413_626_7, 0.222_381_034_453_374_5),
    (-0.525_532_409_916_329, 0.313_706_645_877_887_3),
    (-0.183_434_642_495_649_8, 0.362_683_783_378_362),
    (0.183_434_642_495_649_8, 0.362_683_783_378_362),
    (0.525_532_409_916_329, 0.313_706_645_877_887_3),
    (0.796_666_477_413_626_7, 0.222_381_034_453_374_5),
    (0.960_289_856_497_536_3, 0.101_228_536_290_376_3),
];

fn wall_quadrature(rule: &[(f64, f64)]) -> (Vec<IlluminanceProbe>, Vec<f64>) {
    let mut probes = Vec::with_capacity(6 * rule.len() * rule.len());
    let mut areas = Vec::with_capacity(probes.capacity());
    for axis in 0..3 {
        for side in [-1.0, 1.0] {
            let mut normal = [0.0; 3];
            normal[axis] = -side;
            for &(u, wu) in rule {
                for &(v, wv) in rule {
                    let mut point = [0.0; 3];
                    point[axis] = side;
                    point[(axis + 1) % 3] = u;
                    point[(axis + 2) % 3] = v;
                    probes.push(
                        IlluminanceProbe::try_new(
                            Point3::new(point[0], point[1], point[2]),
                            Vec3::new(normal[0], normal[1], normal[2]),
                        )
                        .unwrap(),
                    );
                    areas.push(wu * wv);
                }
            }
        }
    }
    (probes, areas)
}

fn direct_reference_lux(probe: IlluminanceProbe) -> f64 {
    let p = probe.position_m();
    let n = probe.normal();
    FIXTURES
        .iter()
        .map(|&(source, flux)| {
            let d = [source[0] - p.x, source[1] - p.y, source[2] - p.z];
            let r2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
            let projected = n.x * d[0] + n.y * d[1] + n.z * d[2];
            assert!(projected > 0.0);
            (flux / (4.0 * PI)) * projected / (r2 * r2.sqrt())
        })
        .sum()
}

fn integrated_flux(
    scene: &Scene,
    sources: &[PhotometricSource],
    rule: &[(f64, f64)],
    depth: u32,
    samples: u32,
    cx: &Cx<'_>,
) -> (f64, f64) {
    let (probes, areas) = wall_quadrature(rule);
    let mut job = PathIlluminanceJob::new(
        scene,
        sources,
        PathMedium::Vacuum,
        &probes,
        TransportOptions::default(),
        IlluminanceSampling {
            samples_per_probe: samples,
            max_depth: depth,
            sampler: Sampler::Iid,
            seed: 0x504f_494e_545f_454e,
        },
        cx,
    )
    .unwrap();
    job.advance(cx, usize::MAX).unwrap();
    assert!(job.is_complete());
    let mut direct_lm = 0.0;
    let mut continuous_lm = 0.0;
    for ((probe, area), value) in probes.iter().zip(areas).zip(job.estimates()) {
        let value = value.unwrap();
        close(
            value.direct_lux,
            direct_reference_lux(*probe),
            3.0e-12,
            "off-center point-fixture direct illuminance",
        );
        assert_eq!(value.samples, samples);
        assert_eq!(value.numerical.kind, fs_evidence::NumericalKind::NoClaim);
        assert!(value.total_error_bound_lux.is_none());
        assert!(value.estimated_standard_error_lux.is_some());
        direct_lm += area * value.direct_lux;
        continuous_lm += area * value.continuous_lux;
    }
    (direct_lm, continuous_lm)
}

#[test]
fn g1_point_fixture_multiple_interreflection_conserves_enclosed_flux() {
    with_cx(|cx, _| {
        let emitted_lm = 1000.0;
        // Flatten the photopic power using the original CIE table, so spectral
        // packet sampling cannot hide geometry/transport errors in noise.
        // This declared synthetic spectrum is not a fixture calibration claim.
        let spectrum = SourceSpectrum::try_from_bins(core::array::from_fn(|bin| {
            1.0 / super::spectral_reference::cie_reference_xyz(382.5 + 5.0 * bin as f64)[1]
        }))
        .unwrap();
        let sources = FIXTURES.map(|(position, lumens)| {
            PhotometricSource::try_new(
                RigidTransform::try_new([0.0, 0.0, 0.0, 1.0], position).unwrap(),
                lumens,
                1.0,
                AngularDistribution::isotropic(),
                spectrum.clone(),
            )
            .unwrap()
        });
        // A closed surface intercepts every emitted lumen. Independent direct
        // quadrature is 1000.10882094 lm (6x6) and 1000.01505287 lm (8x8).
        // Its residual is distinct from path sampling and depth truncation.
        let mut direct_residuals = Vec::new();
        for rule in [&GAUSS_6[..], &GAUSS_8[..]] {
            let (probes, areas) = wall_quadrature(rule);
            let integrated: f64 = probes
                .into_iter()
                .zip(areas)
                .map(|(probe, area)| area * direct_reference_lux(probe))
                .sum();
            direct_residuals.push((integrated - emitted_lm).abs());
        }
        assert!(direct_residuals[0] < 0.2);
        assert!(direct_residuals[1] < 0.03);
        assert!(direct_residuals[1] < direct_residuals[0] / 5.0);

        let mut scene = empty_scene([0.0; 3], [0.0, 0.0, 1.0]);
        let mut walls = super::area::diffuse_emitting_cube(0.0);
        // No area-emitter binding: all light must originate at the point
        // fixtures and reach other walls through actual diffuse transport.
        walls.material = Material::SpectralDiffuse {
            reflectance: SpectralReflectance::constant(0.0).unwrap(),
        };
        scene.primitives.push(walls);
        let (direct, continuous) = integrated_flux(&scene, &sources, &GAUSS_8, 1, 2, cx);
        close(
            direct,
            emitted_lm,
            3.0e-5,
            "black cavity intercepts emitted lumens",
        );
        assert_eq!(continuous, 0.0, "black walls cannot create reflected light");

        let rho = 0.6_f64;
        scene.primitives[0].material = Material::SpectralDiffuse {
            reflectance: SpectralReflectance::constant(rho).unwrap(),
        };
        let mut depth_four = None;
        for (depth, rule) in [
            (1, &GAUSS_8[..]),
            (4, &GAUSS_6[..]),
            (4, &GAUSS_8[..]),
            (12, &GAUSS_8[..]),
        ] {
            let (direct, continuous) = integrated_flux(&scene, &sources, rule, depth, 256, cx);
            // The receiver consumes no vertex. Each subsequent wall reflects
            // rho of its incoming power, and the closed cavity loses no ray.
            // Therefore orders 0..D contribute Phi * sum_{k=0}^D rho^k.
            let remainder_fraction = rho.powi(depth as i32 + 1);
            let expected = emitted_lm * (1.0 - remainder_fraction) / (1.0 - rho);
            close(
                direct + continuous,
                expected,
                0.025,
                "all-wall incident flux series",
            );
            close(
                continuous,
                expected - emitted_lm,
                0.025,
                "reflected point-fixture flux",
            );
            let absorbed = (1.0 - rho) * (direct + continuous);
            close(
                absorbed + emitted_lm * remainder_fraction,
                emitted_lm,
                0.025,
                "absorbed power plus the known finite-depth residual",
            );
            if depth == 4 {
                if let Some(coarse) = depth_four {
                    close(
                        direct + continuous,
                        coarse,
                        0.02,
                        "whole-wall quadrature refinement",
                    );
                }
                depth_four = Some(direct + continuous);
                // Infinite-bounce illumination would violate this finite-depth
                // expectation by 8.4%, well beyond the regression tolerance.
                assert!(absorbed < 0.96 * emitted_lm);
            }
            if depth == 12 {
                close(
                    direct + continuous,
                    emitted_lm / (1.0 - rho),
                    0.03,
                    "approach to the independent infinite-bounce solution",
                );
            }
        }
        // Fixed quadrature and IID regression tolerances are not numerical
        // enclosures or confidence bounds. The geometric tail above is valid
        // only for this closed, uniform-gray, nonabsorbing reference model.
    });
}
