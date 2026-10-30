//! Explicit, host-specific qualification of reflected illuminance probe work.

use super::*;
use fs_render::photometry::IlluminanceProbe;
use fs_render::photometry::path::{IlluminanceSampling, PathIlluminanceJob};

/// This measures one fixed probe, one reflecting surface and 300 point sources,
/// not general scene throughput. Execute explicitly with `--ignored --nocapture`.
#[test]
#[ignore = "host-specific calibrated illuminance probe throughput qualification"]
fn photometric_probe_300_fixture_throughput() {
    with_cx(|cx, _| {
        let mut scene = empty_scene([0.0; 3], [0.0, 0.0, 1.0]);
        scene
            .primitives
            .push(floor(2.0, 100.0, LiftedSpectrum { c: [0.0; 3] }));
        let sources: Vec<_> = (0..300)
            .map(|_| isotropic([0.0, 0.0, 1.0], 1.0, SourceSpectrum::equal_energy()))
            .collect();
        let probes =
            [
                IlluminanceProbe::try_new(Point3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 0.0, 1.0))
                    .unwrap(),
            ];
        let sampling = IlluminanceSampling {
            samples_per_probe: 8192,
            max_depth: 1,
            sampler: Sampler::OwenSobol,
            seed: 173,
        };
        let mut job = PathIlluminanceJob::new(
            &scene,
            &sources,
            PathMedium::Vacuum,
            &probes,
            TransportOptions::default(),
            sampling,
            cx,
        )
        .unwrap();

        let start = std::time::Instant::now();
        let accepted = job.advance(cx, usize::MAX).unwrap();
        let elapsed = start.elapsed().as_secs_f64();
        assert!(job.is_complete());
        assert_eq!(accepted, 8192);
        assert_eq!(job.completed_path_samples(), 8192);
        let estimate = job.probe_estimate(0).unwrap();
        assert_eq!(estimate.samples, 8192);

        // The sources are 1 m above the upward-facing probe: direct E = I.
        let intensity_cd = sources.len() as f64;
        close(
            estimate.direct_lux,
            intensity_cd,
            1.0e-12,
            "300-fixture exact direct probe illuminance",
        );

        // For an infinite ceiling at z=2, source z=1 and rho=0.5, integrate
        // rho*I/pi * 1/(r^2+1)^(3/2) * 4/(r^2+4)^2 over ceiling area:
        // E = 8*rho*I * integral_1^infinity dt/[t^2*(t^2+3)^2]
        //   = rho*I * (1 - 4*pi/(9*sqrt(3))).
        // The finite square contains the radius-100 disk. Its missing positive
        // contribution is below 8*rho*I/(5*100^5) = 2.4e-8 lx. This geometric
        // truncation estimate is not a bound on the sampled result's total error.
        let reflected_lux = 0.5 * intensity_cd * (1.0 - 4.0 * PI / (9.0 * 3.0_f64.sqrt()));
        close(
            estimate.continuous_lux,
            reflected_lux,
            0.02,
            "300-fixture reflected probe illuminance independent of direct term",
        );
        assert_eq!(
            estimate.illuminance_lux.to_bits(),
            (estimate.direct_lux + estimate.continuous_lux).to_bits()
        );
        assert!(estimate.estimated_standard_error_lux.is_none());
        assert!(estimate.total_error_bound_lux.is_none());

        // Most cosine-sampled paths hit the large finite ceiling. Count these
        // connections approximately: a few directions can miss the square.
        println!(
            "photometric probe: fixtures={} probes=1 spp=8192 depth=1 paths={accepted} elapsed_s={elapsed:.6} paths_per_s={:.3} approximate_fixture_vertex_connections_per_s={:.3} direct_lux={:.9} reflected_lux={:.9} reflected_reference_lux={reflected_lux:.9}",
            sources.len(),
            accepted as f64 / elapsed,
            (accepted * sources.len()) as f64 / elapsed,
            estimate.direct_lux,
            estimate.continuous_lux,
        );
    });
}
