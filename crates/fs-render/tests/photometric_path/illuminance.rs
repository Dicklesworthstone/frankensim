//! Hemispherical illuminance references independent of receiver material.

use super::*;
use fs_blake3::hash_domain;
use fs_render::instances::{GeometryInstance, SharedGeometry};
use fs_render::photometry::path::{
    IlluminanceSampling, PathIlluminanceEstimate, PathIlluminanceJob,
};
use fs_render::photometry::{CalculationGrid, IlluminanceProbe};

fn probe(point: [f64; 3], normal: [f64; 3]) -> IlluminanceProbe {
    IlluminanceProbe::try_new(
        Point3::new(point[0], point[1], point[2]),
        Vec3::new(normal[0], normal[1], normal[2]),
    )
    .unwrap()
}

fn sampling(depth: u32, count: u32) -> IlluminanceSampling {
    IlluminanceSampling {
        samples_per_probe: count,
        max_depth: depth,
        sampler: Sampler::Iid,
        seed: 173,
    }
}

fn measured(job: &mut PathIlluminanceJob<'_, '_>, cx: &Cx<'_>) -> PathIlluminanceEstimate {
    assert!(job.probe_estimate(0).is_none());
    job.advance(cx, usize::MAX).unwrap();
    assert!(job.is_complete());
    let estimate = job.probe_estimate(0).unwrap();
    assert_eq!(estimate.samples, job.samples_per_probe());
    assert_eq!(
        estimate.illuminance_lux.to_bits(),
        (estimate.direct_lux + estimate.continuous_lux).to_bits()
    );
    assert!(estimate.total_error_bound_lux.is_none());
    assert_eq!(estimate.numerical.kind, fs_evidence::NumericalKind::NoClaim);
    estimate
}

fn upward_fixture(candela: f64) -> PhotometricSource {
    PhotometricSource::try_new(
        RigidTransform::try_new([0.0, 0.0, 0.0, 1.0], [0.0, 0.0, 1.0]).unwrap(),
        2.0 * PI * candela,
        1.0,
        AngularDistribution::uniform_cone(FRAC_PI_2).unwrap(),
        SourceSpectrum::equal_energy(),
    )
    .unwrap()
}

fn illuminance_bits(estimate: PathIlluminanceEstimate) -> ([u64; 3], u32, Option<u64>) {
    (
        [
            estimate.illuminance_lux.to_bits(),
            estimate.direct_lux.to_bits(),
            estimate.continuous_lux.to_bits(),
        ],
        estimate.samples,
        estimate.estimated_standard_error_lux.map(f64::to_bits),
    )
}

#[test]
fn g1_exact_direct_fixture_is_not_counted_again_by_hemisphere_paths() {
    with_cx(|cx, _| {
        let scene = empty_scene([0.0; 3], [0.0, 0.0, 1.0]);
        let sources = [isotropic(
            [0.0, 0.0, 3.0],
            90.0,
            SourceSpectrum::equal_energy(),
        )];
        let probes = [probe([0.0; 3], [0.0, 0.0, 1.0])];
        for sampler in [Sampler::Iid, Sampler::OwenSobol] {
            let mut policy = sampling(4, 32);
            policy.sampler = sampler;
            let mut job = PathIlluminanceJob::new(
                &scene,
                &sources,
                PathMedium::Vacuum,
                &probes,
                TransportOptions::default(),
                policy,
                cx,
            )
            .unwrap();
            assert_eq!(job.probe_count(), 1);
            let result = measured(&mut job, cx);
            close(result.direct_lux, 10.0, 2.0e-13, "90 cd at 3 m");
            assert_eq!(result.continuous_lux, 0.0);
            assert_eq!(
                result.illuminance_lux.to_bits(),
                result.direct_lux.to_bits()
            );
            assert_eq!(
                result.estimated_standard_error_lux,
                if sampler == Sampler::Iid {
                    Some(0.0)
                } else {
                    None
                }
            );
        }
    });
}

#[test]
fn g1_reflected_ceiling_matches_irradiance_integral_without_receiver_albedo() {
    with_cx(|cx, _| {
        let sources = [upward_fixture(100.0)];
        let probes = [probe([0.0; 3], [0.0, 0.0, 1.0])];
        // At ceiling point (x,y,2): E_light=100/(r²+1)^(3/2).
        // Its rho=0.5 diffuse radiance is E_light/(2*pi); the area-to-probe
        // geometry is 4/(r²+4)^2. There is NO receiver BRDF in illuminance.
        let cells = 256;
        let width = 8.0 / f64::from(cells);
        let mut integral = 0.0;
        for ix in 0..cells {
            for iy in 0..cells {
                let x = -4.0 + (f64::from(ix) + 0.5) * width;
                let y = -4.0 + (f64::from(iy) + 0.5) * width;
                let r2 = x * x + y * y;
                integral += width * width / ((r2 + 1.0) * (r2 + 1.0).sqrt() * (r2 + 4.0).powi(2));
            }
        }
        let expected = 0.5 * 100.0 * 4.0 * integral / PI;
        let mut first_bits = None;
        for receiver_coefficient in [-2.0, 2.0] {
            let mut scene = empty_scene([0.0, 0.0, 0.5], [0.0, 0.0, -1.0]);
            scene.primitives.push(floor(
                0.0,
                10.0,
                LiftedSpectrum {
                    c: [0.0, 0.0, receiver_coefficient],
                },
            ));
            scene
                .primitives
                .push(floor(2.0, 4.0, LiftedSpectrum { c: [0.0; 3] }));
            let mut job = PathIlluminanceJob::new(
                &scene,
                &sources,
                PathMedium::Vacuum,
                &probes,
                TransportOptions::default(),
                sampling(1, 65_536),
                cx,
            )
            .unwrap();
            let result = measured(&mut job, cx);
            assert_eq!(result.direct_lux, 0.0);
            close(
                result.continuous_lux,
                expected,
                0.03,
                "one-reflection ceiling irradiance",
            );
            if let Some(bits) = first_bits {
                assert_eq!(
                    result.illuminance_lux.to_bits(),
                    bits,
                    "a virtual irradiance receiver has no reflectance factor"
                );
            } else {
                first_bits = Some(result.illuminance_lux.to_bits());
            }
        }
    });
}

#[test]
fn g1_direct_term_preserves_wavelength_dependent_grid_extinction() {
    with_cx(|cx, _| {
        let scene = empty_scene([0.0; 3], [0.0, 0.0, 1.0]);
        let density = [2.0];
        let medium = AbsorbingGrid::try_new(
            cx,
            VolumeGrid::new([1; 3], &density, [-1.0, -1.0, 1.0], [2.0, 2.0, 1.0]),
            MassAbsorptionSpectrum::table(&[(380.0, 0.1), (780.0, 1.1)]).unwrap(),
            None,
            medium_budget(),
        )
        .unwrap();
        let probes = [probe([0.0; 3], [0.0, 0.0, 1.0])];
        let mut observed = Vec::new();
        for bin in [12, 60] {
            let mut bins = [0.0; SPECTRAL_BINS];
            bins[bin] = 1.0;
            let sources = [isotropic(
                [0.0, 0.0, 3.0],
                90.0,
                SourceSpectrum::try_from_bins(bins).unwrap(),
            )];
            let mut job = PathIlluminanceJob::new(
                &scene,
                &sources,
                PathMedium::Absorbing(&medium),
                &probes,
                TransportOptions::default(),
                sampling(2, 16),
                cx,
            )
            .unwrap();
            let result = measured(&mut job, cx);
            let kappa = 0.1 + (wavelength_nm(bin) - 380.0) / 400.0;
            close(
                result.illuminance_lux,
                10.0 * (-2.0 * kappa).exp(),
                3.0e-12,
                "spectrally filtered direct illuminance",
            );
            assert_eq!(result.continuous_lux, 0.0);
            observed.push(result.illuminance_lux);
        }
        assert!(observed[0] > 3.0 * observed[1]);
    });
}

#[test]
fn g1_single_scattered_illuminance_matches_independent_volume_integral() {
    with_cx(|cx, _| {
        let scene = empty_scene([0.0; 3], [0.0, 0.0, 1.0]);
        let sources = [isotropic(
            [0.0, 0.0, 2.0],
            100.0,
            SourceSpectrum::equal_energy(),
        )];
        let probes = [probe([0.0; 3], [0.0, 0.0, 1.0])];
        let density = [1.0];
        let medium = ScatteringGrid::try_new(
            cx,
            AbsorbingGrid::try_new(
                cx,
                VolumeGrid::new([1; 3], &density, [-1.0, -1.0, 0.5], [2.0, 2.0, 0.5]),
                MassAbsorptionSpectrum::constant(0.2).unwrap(),
                None,
                medium_budget(),
            )
            .unwrap(),
            MassScatteringSpectrum::constant(0.3).unwrap(),
            PhaseFunction::Isotropic,
        )
        .unwrap();
        // Integral over the actual box, rather than over sampled ray paths:
        // I*sigma_s*p * exp[-sigma_t*(ell_light+ell_probe)]
        //               * z/(r_probe^3*r_light^2) dV.
        // Both exits are through horizontal faces because both endpoints lie
        // on the central axis outside the slab.
        let cells = 64;
        let dx = 2.0 / f64::from(cells);
        let dz = 0.5 / f64::from(cells);
        let mut integral = 0.0;
        for ix in 0..cells {
            for iy in 0..cells {
                for iz in 0..cells {
                    let x = -1.0 + (f64::from(ix) + 0.5) * dx;
                    let y = -1.0 + (f64::from(iy) + 0.5) * dx;
                    let z = 0.5 + (f64::from(iz) + 0.5) * dz;
                    let rp2 = x * x + y * y + z * z;
                    let rl2 = x * x + y * y + (2.0 - z) * (2.0 - z);
                    let lp = rp2.sqrt() * (z - 0.5) / z;
                    let ll = rl2.sqrt() * (1.0 - z) / (2.0 - z);
                    integral +=
                        dx * dx * dz * (-0.5 * (lp + ll)).exp() * z / (rp2 * rp2.sqrt() * rl2);
                }
            }
        }
        let expected_scatter = 100.0 * 0.3 * integral / (4.0 * PI);
        let mut job = PathIlluminanceJob::new(
            &scene,
            &sources,
            PathMedium::Scattering(&medium),
            &probes,
            TransportOptions::default(),
            sampling(1, 131_072),
            cx,
        )
        .unwrap();
        let result = measured(&mut job, cx);
        close(
            result.direct_lux,
            25.0 * (-0.25_f64).exp(),
            2.0e-13,
            "unscattered slab direct term",
        );
        close(
            result.continuous_lux,
            expected_scatter,
            0.04,
            "single-scattering volume irradiance",
        );
        assert!(result.continuous_lux > 0.0);
    });
}

#[test]
fn g1_thermal_volume_emission_reaches_a_source_free_illuminance_probe() {
    with_cx(|cx, _| {
        let scene = empty_scene([0.0; 3], [0.0, 0.0, 1.0]);
        let probes = [probe([0.0; 3], [0.0, 0.0, 1.0])];
        let density = [1.0];
        let temperatures = [2000.0];
        // Independent Planck integral in W/(m² sr nm), then photopic units.
        let n_lambda = 4096;
        let dl = 400.0 / f64::from(n_lambda);
        let mut blackbody_cd_m2 = 0.0;
        for i in 0..n_lambda {
            let nm = 380.0 + (f64::from(i) + 0.5) * dl;
            let wavelength = nm * 1.0e-9;
            let h = 6.626_070_15e-34;
            let c = 299_792_458.0;
            let k = 1.380_649e-23;
            let planck = 2.0 * h * c * c
                / (wavelength.powi(5) * (h * c / (wavelength * k * 2000.0)).exp_m1())
                * 1.0e-9;
            blackbody_cd_m2 +=
                683.0 * super::spectral_reference::cie_reference_xyz(nm)[1] * planck * dl;
        }
        // E = 2*pi*Bv * integral_0^1 mu [1-exp(-tau/mu)] dmu.
        // The finite side faces affect only a grazing angular region whose
        // cosine-weighted measure is below 1e-10 for this aspect ratio.
        let angular: f64 = (0..8192)
            .map(|i| {
                let mu = (f64::from(i) + 0.5) / 8192.0;
                mu * -(-0.4 / mu).exp_m1() / 8192.0
            })
            .sum();
        let expected = 2.0 * PI * blackbody_cd_m2 * angular;
        let options = TransportOptions {
            max_distance_m: 1.0e8,
            ..TransportOptions::default()
        };
        // Equal optical depth must also survive a receiver-adjacent cell
        // thinner than the geometric ray epsilon. Transport starts at the
        // physical receiver, even when its numerical ray origin is beyond it.
        for (start, thickness) in [(1.0, 1.0), (0.0, 1.0e-8)] {
            let medium = AbsorbingGrid::try_new(
                cx,
                VolumeGrid::new(
                    [1; 3],
                    &density,
                    [-1.0e6, -1.0e6, start],
                    [2.0e6, 2.0e6, thickness],
                ),
                MassAbsorptionSpectrum::constant(0.4 / thickness).unwrap(),
                Some(&temperatures),
                medium_budget(),
            )
            .unwrap();
            let mut job = PathIlluminanceJob::new(
                &scene,
                &[],
                PathMedium::Absorbing(&medium),
                &probes,
                options,
                sampling(1, 32_768),
                cx,
            )
            .unwrap();
            let result = measured(&mut job, cx);
            assert_eq!(result.direct_lux, 0.0);
            close(
                result.continuous_lux,
                expected,
                0.025,
                "hemispherical thermal-medium illuminance, including thin adjacent cells",
            );
        }
    });
}

#[test]
fn g4_g5_probe_grid_progress_cancellation_and_iid_error_units_replay() {
    with_cx(|cx, gate| {
        let mut scene = empty_scene([0.0; 3], [0.0, 0.0, 1.0]);
        scene
            .primitives
            .push(floor(2.0, 4.0, LiftedSpectrum { c: [0.0; 3] }));
        let sources = [upward_fixture(100.0)];
        let grid = CalculationGrid::try_new(
            Point3::new(-0.5, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            2,
            1,
        )
        .unwrap();
        let probes = [grid.probe(0).unwrap(), grid.probe(1).unwrap()];
        let policy = sampling(1, 64);
        let mut reference = PathIlluminanceJob::for_grid(
            &scene,
            &sources,
            PathMedium::Vacuum,
            grid,
            TransportOptions::default(),
            policy,
            cx,
        )
        .unwrap();
        assert_eq!(reference.probe_count(), 2);
        reference.advance(cx, usize::MAX).unwrap();
        let reference_bits: Vec<_> = reference
            .estimates()
            .map(|estimate| illuminance_bits(estimate.unwrap()))
            .collect();
        let mut job = PathIlluminanceJob::new(
            &scene,
            &sources,
            PathMedium::Vacuum,
            &probes,
            TransportOptions::default(),
            policy,
            cx,
        )
        .unwrap();
        assert!(job.probe_estimate(0).is_none());
        assert!(job.probe_estimate(2).is_none());
        job.advance(cx, 1).unwrap();
        let first = job.probe_estimate(0).unwrap();
        assert_eq!(first.samples, 1);
        assert!(first.estimated_standard_error_lux.is_none());
        assert!(job.probe_estimate(1).is_none());
        job.advance(cx, 2).unwrap();
        let second = job.probe_estimate(0).unwrap();
        assert_eq!(second.samples, 2);
        close(
            second.estimated_standard_error_lux.unwrap(),
            (second.continuous_lux - first.continuous_lux).abs(),
            1.0e-12,
            "two-sample IID standard error in lux",
        );
        let prefix: Vec<_> = job
            .estimates()
            .map(|estimate| illuminance_bits(estimate.unwrap()))
            .collect();
        gate.request();
        assert!(matches!(
            job.advance(cx, 8),
            Err(PhotometryError::Cancelled)
        ));
        assert_eq!(job.completed_path_samples(), 3);
        assert_eq!(
            job.estimates()
                .map(|estimate| illuminance_bits(estimate.unwrap()))
                .collect::<Vec<_>>(),
            prefix
        );
        with_cx(|resumed, _| {
            while !job.is_complete() {
                job.advance(resumed, 7).unwrap();
            }
        });
        assert_eq!(
            job.estimates()
                .map(|estimate| illuminance_bits(estimate.unwrap()))
                .collect::<Vec<_>>(),
            reference_bits
        );
    });
}

#[test]
fn g3_near_south_pole_probe_ignores_the_unused_scene_camera() {
    with_cx(|cx, _| {
        let mut scene = empty_scene([f64::NAN; 3], [0.0; 3]);
        scene.camera.up = Vec3::new(0.0, 0.0, 0.0);
        scene.camera.half_tan = f64::NAN;
        scene
            .primitives
            .push(floor(-1.0, 4.0, LiftedSpectrum { c: [0.0; 3] }));
        let sources = [isotropic([0.0; 3], 90.0, SourceSpectrum::equal_energy())];
        let probes = [probe([0.0, 0.0, 3.0], [0.0003, 0.0, -1.0])];
        let normal_length = (1.0_f64 + 0.0003_f64.powi(2)).sqrt();
        // Symmetry cancels the receiving normal's tiny x component. At the
        // floor the fixture is 1 m above and the receiver 4 m above; integrate
        // the reflected irradiance independently in area coordinates.
        let cells = 256;
        let width = 8.0 / f64::from(cells);
        let mut integral = 0.0;
        for ix in 0..cells {
            for iy in 0..cells {
                let x = -4.0 + (f64::from(ix) + 0.5) * width;
                let y = -4.0 + (f64::from(iy) + 0.5) * width;
                let r2 = x * x + y * y;
                integral += width * width / ((r2 + 1.0) * (r2 + 1.0).sqrt() * (r2 + 16.0).powi(2));
            }
        }
        let reflected = 0.5 * 90.0 * 16.0 * integral / (PI * normal_length);
        let mut job = PathIlluminanceJob::new(
            &scene,
            &sources,
            PathMedium::Vacuum,
            &probes,
            TransportOptions::default(),
            sampling(1, 32_768),
            cx,
        )
        .unwrap();
        let result = measured(&mut job, cx);
        close(
            result.direct_lux,
            10.0 / normal_length,
            2.0e-13,
            "near-south-pole receiving cosine",
        );
        close(
            result.continuous_lux,
            reflected,
            0.03,
            "near-south-pole reflected irradiance",
        );
    });
}

#[test]
fn g3_distant_mesh_instances_match_native_probes_and_refuse_short_ranges() {
    with_cx(|cx, _| {
        let sources = [upward_fixture(1.0e12)];
        let probes = [probe([0.0; 3], [0.0, 0.0, 1.0])];
        let reflectance = LiftedSpectrum { c: [0.0; 3] };
        let native = floor(12_000.0, 1000.0, reflectance);
        let Shape::Mesh(local) = floor(0.0, 1000.0, reflectance).shape else {
            unreachable!()
        };
        let instance = Primitive {
            shape: Shape::Instance(
                GeometryInstance::try_new(
                    1,
                    hash_domain("org.frankensim.test.photometric-probe", b"distant-ceiling"),
                    SharedGeometry::mesh(local),
                    RigidTransform::try_new([0.0, 0.0, 0.0, 1.0], [0.0, 0.0, 12_000.0]).unwrap(),
                )
                .unwrap(),
            ),
            material: Material::Lambertian { reflectance },
            emission: None,
        };
        let mut reference = None;
        for primitive in [native, instance] {
            let mut scene = empty_scene([0.0; 3], [0.0, 0.0, 1.0]);
            scene.primitives.push(primitive);
            let options = TransportOptions {
                max_distance_m: 20_000.0,
                ..TransportOptions::default()
            };
            let mut job = PathIlluminanceJob::new(
                &scene,
                &sources,
                PathMedium::Vacuum,
                &probes,
                options,
                sampling(1, 4096),
                cx,
            )
            .unwrap();
            let result = measured(&mut job, cx);
            assert_eq!(result.direct_lux, 0.0);
            assert!(
                result.continuous_lux > 1.0,
                "distant instances must not silently disappear beyond the legacy 10km limit"
            );
            if let Some(expected) = reference {
                close(
                    result.continuous_lux,
                    expected,
                    1.0e-10,
                    "native and instanced distant ceiling",
                );
            } else {
                reference = Some(result.continuous_lux);
            }
            let short = TransportOptions {
                max_distance_m: 11_000.0,
                ..options
            };
            let mut refused = PathIlluminanceJob::new(
                &scene,
                &sources,
                PathMedium::Vacuum,
                &probes,
                short,
                sampling(1, 4096),
                cx,
            )
            .unwrap();
            assert!(matches!(
                refused.advance(cx, usize::MAX),
                Err(PhotometryError::Geometry(TracerError::FixturePathRange))
            ));
            assert!(refused.completed_path_samples() < 4096);
        }
    });
}
