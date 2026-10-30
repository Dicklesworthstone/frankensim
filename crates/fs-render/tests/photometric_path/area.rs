//! Independent radiometric references for calibrated native-mesh emitters.

use super::*;
use fs_render::photometry::IlluminanceProbe;
use fs_render::photometry::area::{EmissionSidedness, PhotometricAreaEmitter, PhotometricEmitters};
use fs_render::photometry::path::{
    IlluminanceSampling, PathIlluminanceEstimate, PathIlluminanceJob,
};

fn panel(z: f64, a: f64, b: f64, positive_z: bool) -> Primitive {
    Primitive {
        shape: Shape::Mesh(TriMesh::new(
            vec![[-a, -b, z], [a, -b, z], [a, b, z], [-a, b, z]],
            if positive_z {
                vec![[0, 1, 2], [0, 2, 3]]
            } else {
                vec![[0, 2, 1], [0, 3, 2]]
            },
        )),
        material: Material::Lambertian {
            reflectance: LiftedSpectrum { c: [0.0; 3] },
        },
        emission: None,
    }
}

fn emitter(index: usize, flux_lm: f64, sidedness: EmissionSidedness) -> PhotometricAreaEmitter {
    PhotometricAreaEmitter::try_new(
        index,
        flux_lm,
        1.0,
        SourceSpectrum::equal_energy(),
        sidedness,
    )
    .unwrap()
}

fn area_pixel(
    scene: &Scene,
    points: &[PhotometricSource],
    surfaces: &[PhotometricAreaEmitter],
    medium: PathMedium<'_, '_>,
    policy: Settings,
    cx: &Cx<'_>,
) -> PathLuminanceEstimate {
    let mut job = PathLuminanceJob::new_with_emitters(
        scene,
        PhotometricEmitters { points, surfaces },
        medium,
        TransportOptions::default(),
        policy,
        cx,
    )
    .unwrap();
    job.advance(cx, usize::MAX).unwrap();
    assert!(job.is_complete());
    let result = job.pixel_estimate(0).unwrap();
    assert!(result.total_error_bound_cd_m2.is_none());
    result
}

fn area_probe(
    scene: &Scene,
    points: &[PhotometricSource],
    surfaces: &[PhotometricAreaEmitter],
    probe: IlluminanceProbe,
    samples: u32,
    cx: &Cx<'_>,
) -> PathIlluminanceEstimate {
    let probes = [probe];
    let mut job = PathIlluminanceJob::new_with_emitters(
        scene,
        PhotometricEmitters { points, surfaces },
        PathMedium::Vacuum,
        &probes,
        TransportOptions::default(),
        IlluminanceSampling {
            samples_per_probe: samples,
            max_depth: 1,
            sampler: Sampler::OwenSobol,
            seed: 173,
        },
        cx,
    )
    .unwrap();
    job.advance(cx, usize::MAX).unwrap();
    assert!(job.is_complete());
    let result = job.probe_estimate(0).unwrap();
    assert!(result.total_error_bound_lux.is_none());
    assert!(result.estimated_standard_error_lux.is_none());
    result
}

fn axial_rectangle_lux(luminance: f64, a: f64, b: f64, z: f64) -> f64 {
    // Exact integral of L*z²/(x²+y²+z²)² over [-a,a] x [-b,b].
    let az = (a * a + z * z).sqrt();
    let bz = (b * b + z * z).sqrt();
    2.0 * luminance * (a / az * (b / az).atan() + b / bz * (a / bz).atan())
}

fn projected_polygon_lux(
    luminance: f64,
    corners: [[f64; 3]; 4],
    receiver: [f64; 3],
    normal: [f64; 3],
) -> f64 {
    // Stokes' theorem on the unit sphere gives
    // integral_Omega omega dOmega = (1/2) sum_edges theta * unit(u cross v).
    // Dotting this vector with the receiving normal integrates the cosine.
    // The polygon must be wholly above the horizon, with positive winding
    // as seen from the receiver. No renderer PDF, ray or quadrature is used.
    let dot = |a: [f64; 3], b: [f64; 3]| a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    let normal_length = dot(normal, normal).sqrt();
    let normal = normal.map(|value| value / normal_length);
    let directions = corners.map(|corner| {
        let delta = std::array::from_fn(|axis| corner[axis] - receiver[axis]);
        let distance = dot(delta, delta).sqrt();
        let direction = delta.map(|value| value / distance);
        assert!(
            dot(normal, direction) > 0.0,
            "reference requires no horizon clipping"
        );
        direction
    });
    let mut projection = 0.0;
    for edge in 0..4 {
        let u = directions[edge];
        let v = directions[(edge + 1) % 4];
        let cross = [
            u[1] * v[2] - u[2] * v[1],
            u[2] * v[0] - u[0] * v[2],
            u[0] * v[1] - u[1] * v[0],
        ];
        let sine = dot(cross, cross).sqrt();
        projection += dot(normal, cross) * sine.atan2(dot(u, v)) / (2.0 * sine);
    }
    assert!(projection > 0.0 && projection <= PI);
    luminance * projection
}

fn diffuse_emitting_cube(coefficient: f64) -> Primitive {
    // Exactly represented planar enclosure: 6 faces * 4 m², inward winding.
    // Unlike a tessellated sphere, this reference has no curvature error.
    Primitive {
        shape: Shape::Mesh(TriMesh::new(
            vec![
                [-1.0, -1.0, -1.0],
                [1.0, -1.0, -1.0],
                [1.0, 1.0, -1.0],
                [-1.0, 1.0, -1.0],
                [-1.0, -1.0, 1.0],
                [1.0, -1.0, 1.0],
                [1.0, 1.0, 1.0],
                [-1.0, 1.0, 1.0],
            ],
            vec![
                [0, 1, 2],
                [0, 2, 3],
                [4, 6, 5],
                [4, 7, 6],
                [0, 3, 7],
                [0, 7, 4],
                [1, 5, 6],
                [1, 6, 2],
                [0, 4, 5],
                [0, 5, 1],
                [3, 2, 6],
                [3, 6, 7],
            ],
        )),
        material: Material::Lambertian {
            reflectance: LiftedSpectrum {
                c: [0.0, 0.0, coefficient],
            },
        },
        emission: None,
    }
}

#[test]
fn g1_off_axis_inclined_receiver_matches_exact_projected_solid_angle() {
    with_cx(|cx, _| {
        let mut scene = empty_scene([0.0; 3], [0.0, 0.0, 1.0]);
        scene.primitives.push(panel(2.0, 1.0, 0.5, false));
        let sources = [emitter(0, 100.0 * PI, EmissionSidedness::Front)];
        let corners = [
            [-1.0, -0.5, 2.0],
            [1.0, -0.5, 2.0],
            [1.0, 0.5, 2.0],
            [-1.0, 0.5, 2.0],
        ];
        // Near-field, off-axis and far-field views exercise both geometric
        // cosines. Flux is 100*pi lm over 2 m², hence L = 50 cd/m².
        for (point, normal) in [
            ([0.7, 0.1, 1.7], [0.1, 0.1, 1.0]),
            ([0.6, -0.3, 0.3], [0.35, 0.25, 1.0]),
            ([-1.2, 0.7, -3.0], [0.3, -0.4, 1.0]),
        ] {
            let expected = projected_polygon_lux(50.0, corners, point, normal);
            let probe = IlluminanceProbe::try_new(
                Point3::new(point[0], point[1], point[2]),
                Vec3::new(normal[0], normal[1], normal[2]),
            )
            .unwrap();
            let result = area_probe(&scene, &[], &sources, probe, 32_768, cx);
            close(
                result.illuminance_lux,
                expected,
                0.02,
                "projected solid angle",
            );
            assert_eq!(result.direct_lux, 0.0);
            // Any passive aperture covers at most the receiving hemisphere.
            assert!(result.illuminance_lux <= 1.02 * PI * 50.0);
        }
    });
}

#[test]
fn g1_diffuse_enclosure_matches_reflection_series_and_absorbed_power() {
    with_cx(|cx, _| {
        let emitted = 20.0;
        let flux = 24.0 * PI * emitted;
        let sources = [emitter(0, flux, EmissionSidedness::Front)];
        // c=0 and c=3/4 represent spectrally constant rho=1/2 and 4/5.
        // Every wall sees a complete hemisphere of the same radiance, so
        // L_infinity = Le/(1-rho), independent of position or mesh sampling.
        for (coefficient, rho) in [(0.0, 0.5_f64), (0.75, 0.8_f64)] {
            let mut scene = empty_scene([0.0; 3], [0.0, 0.0, 1.0]);
            scene.primitives.push(diffuse_emitting_cube(coefficient));
            for depth in [1, 4, 12] {
                let mut policy = settings(depth, 8192);
                policy.strategy = DirectStrategy::BsdfOnly;
                let value = area_pixel(&scene, &[], &sources, PathMedium::Vacuum, policy, cx);
                // A hit-only depth D observes emission after 0..D-1 bounces.
                let tail_fraction = rho.powi(depth as i32);
                let expected = emitted * (1.0 - tail_fraction) / (1.0 - rho);
                close(
                    value.luminance_cd_m2,
                    expected,
                    0.002,
                    "diffuse reflection series",
                );
                // Each wall absorbs (1-rho)*pi*L. Missing absorption equals
                // the analytically known omitted tail, not sampling noise.
                let absorbed = (1.0 - rho) * 24.0 * PI * value.luminance_cd_m2;
                close(
                    absorbed + flux * tail_fraction,
                    flux,
                    0.002,
                    "enclosure power balance",
                );
                assert!(value.total_error_bound_cd_m2.is_none());
            }
        }
    });
}

#[test]
fn g1_calibrated_multiple_scattering_preserves_enclosure_energy_balance() {
    with_cx(|cx, _| {
        let emitted = 20.0;
        let rho = 0.5;
        let flux = 24.0 * PI * emitted;
        let mut scene = empty_scene([0.0; 3], [0.0, 0.0, 1.0]);
        scene.primitives.push(diffuse_emitting_cube(0.0));
        let sources = [emitter(0, flux, EmissionSidedness::Front)];
        let density = [1.0];
        let medium = ScatteringGrid::try_new(
            cx,
            AbsorbingGrid::try_new(
                cx,
                VolumeGrid::new([1; 3], &density, [-1.0; 3], [2.0; 3]),
                MassAbsorptionSpectrum::constant(0.0).unwrap(),
                None,
                medium_budget(),
            )
            .unwrap(),
            MassScatteringSpectrum::constant(0.4).unwrap(),
            PhaseFunction::henyey_greenstein(0.6).unwrap(),
        )
        .unwrap();
        // A constant isotropic radiance has equal in- and out-scattering for
        // any normalized reciprocal phase law. With no volume absorption,
        // L=Le+rho*L remains the exact all-order solution (40 cd/m²), and
        // (1-rho)*A*pi*L equals the enclosure's emitted luminous flux.
        // Scattering consumes path vertices too: the vacuum rho^depth tail
        // is NOT a bound here. Two finite depths are regression convergence
        // evidence only; neither the test tolerance nor sampling uncertainty
        // certifies the untraced remainder or a measured haze model.
        let mut previous = None;
        for depth in [32, 64] {
            let value = area_pixel(
                &scene,
                &[],
                &sources,
                PathMedium::Scattering(&medium),
                settings(depth, 32_768),
                cx,
            );
            close(
                value.luminance_cd_m2,
                emitted / (1.0 - rho),
                0.03,
                "conservative scattering in a calibrated diffuse enclosure",
            );
            close(
                (1.0 - rho) * 24.0 * PI * value.luminance_cd_m2,
                flux,
                0.03,
                "emitted versus absorbed lumens with multiple scattering",
            );
            if let Some(previous) = previous {
                close(
                    value.luminance_cd_m2,
                    previous,
                    0.005,
                    "finite-depth convergence in this enclosure",
                );
            }
            previous = Some(value.luminance_cd_m2);
            assert!(value.total_error_bound_cd_m2.is_none());
            assert!(value.estimated_standard_error_cd_m2.is_none());
        }
    });
}

#[test]
fn g1_diffuse_enclosure_mis_and_lux_preserve_terminal_reflection_order() {
    with_cx(|cx, _| {
        let emitted = 20.0;
        let rho = 0.5_f64;
        let mut scene = empty_scene([0.0; 3], [0.0, 0.0, 1.0]);
        scene.primitives.push(diffuse_emitting_cube(0.0));
        let sources = [emitter(0, 24.0 * PI * emitted, EmissionSidedness::Front)];
        let probes = [
            IlluminanceProbe::try_new(Point3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 0.0, 1.0))
                .unwrap(),
            IlluminanceProbe::try_new(Point3::new(-0.3, 0.4, 0.1), Vec3::new(1.0, 0.0, 0.0))
                .unwrap(),
        ];
        for depth in [1, 4] {
            // NEE at the terminal wall includes one more reflected emission
            // than a hit-only path: orders 0..D, in both camera and probe APIs.
            let expected = emitted * (1.0 - rho.powi(depth as i32 + 1)) / (1.0 - rho);
            let value = area_pixel(
                &scene,
                &[],
                &sources,
                PathMedium::Vacuum,
                settings(depth, 32_768),
                cx,
            );
            close(
                value.luminance_cd_m2,
                expected,
                0.025,
                "terminal MIS reflection order",
            );
            let mut job = PathIlluminanceJob::new_with_emitters(
                &scene,
                PhotometricEmitters {
                    points: &[],
                    surfaces: &sources,
                },
                PathMedium::Vacuum,
                &probes,
                TransportOptions::default(),
                IlluminanceSampling {
                    samples_per_probe: 32_768,
                    max_depth: depth,
                    sampler: Sampler::OwenSobol,
                    seed: 173,
                },
                cx,
            )
            .unwrap();
            job.advance(cx, usize::MAX).unwrap();
            for value in job.estimates().map(Option::unwrap) {
                close(
                    value.illuminance_lux,
                    PI * expected,
                    0.025,
                    "uniform-cavity illuminance",
                );
                assert_eq!(value.direct_lux, 0.0);
                assert!(value.estimated_standard_error_lux.is_none());
                assert!(value.total_error_bound_lux.is_none());
            }
        }
        // The 2--2.5% bands here qualify fixed deterministic regression samples.
        // They are not confidence intervals, calibration errors or certificates.
    });
}

#[test]
fn g1_area_source_luminance_uses_actual_area_sides_and_dimmer() {
    with_cx(|cx, _| {
        let mut scene = empty_scene([0.0, 0.0, 2.0], [0.0, 0.0, -1.0]);
        scene.primitives.push(panel(0.0, 1.0, 0.5, true));
        // Area=2 m². Phi=80*pi lm gives L=40 cd/m² on the front side.
        for (sidedness, dimmer, expected) in [
            (EmissionSidedness::Front, 1.0, 40.0),
            (EmissionSidedness::Both, 1.0, 20.0),
            (EmissionSidedness::Front, 0.25, 10.0),
        ] {
            let sources = [PhotometricAreaEmitter::try_new(
                0,
                80.0 * PI,
                dimmer,
                SourceSpectrum::equal_energy(),
                sidedness,
            )
            .unwrap()];
            let front = area_pixel(
                &scene,
                &[],
                &sources,
                PathMedium::Vacuum,
                settings(1, 8192),
                cx,
            );
            close(
                front.luminance_cd_m2,
                expected,
                0.002,
                "emitted surface luminance",
            );
            scene.camera.eye.z = -2.0;
            scene.camera.forward.z = 1.0;
            let back = area_pixel(
                &scene,
                &[],
                &sources,
                PathMedium::Vacuum,
                settings(1, 8192),
                cx,
            );
            if sidedness == EmissionSidedness::Front {
                assert_eq!(back.luminance_cd_m2, 0.0);
            } else {
                close(
                    back.luminance_cd_m2,
                    expected,
                    0.002,
                    "equal two-sided emission",
                );
            }
            scene.camera.eye.z = 2.0;
            scene.camera.forward.z = -1.0;
        }
    });
}

#[test]
fn g1_area_probe_matches_exact_rectangular_illuminance() {
    with_cx(|cx, _| {
        let mut scene = empty_scene([0.0; 3], [0.0, 0.0, 1.0]);
        scene.primitives.push(panel(2.0, 1.0, 0.5, false));
        let sources = [emitter(0, 100.0 * PI, EmissionSidedness::Front)];
        let probe = IlluminanceProbe::try_new(Point3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 0.0, 1.0))
            .unwrap();
        let result = area_probe(&scene, &[], &sources, probe, 32_768, cx);
        assert_eq!(
            result.direct_lux, 0.0,
            "delta fixtures have their own exact component"
        );
        close(
            result.continuous_lux,
            axial_rectangle_lux(50.0, 1.0, 0.5, 2.0),
            0.02,
            "area NEE plus hemisphere hits must integrate once",
        );
    });
}

#[test]
fn g1_polygonal_disk_uses_mesh_area_and_analytic_disk_brackets() {
    with_cx(|cx, _| {
        let count = 8;
        let radius = 1.0;
        let z = 2.0;
        let mut vertices = vec![[0.0; 3]];
        for index in 0..count {
            let angle = 2.0 * PI * f64::from(index) / f64::from(count);
            vertices.push([radius * angle.cos(), radius * angle.sin(), 0.0]);
        }
        let triangles: Vec<_> = (0..count)
            .map(|index| [0, index + 1, (index + 1) % count + 1])
            .collect();
        let area = 0.5 * f64::from(count) * radius * radius * (2.0 * PI / f64::from(count)).sin();
        let mut scene = empty_scene([0.0, 0.0, z], [0.0, 0.0, -1.0]);
        scene.primitives.push(Primitive {
            shape: Shape::Mesh(TriMesh::new(vertices, triangles)),
            material: Material::Lambertian {
                reflectance: LiftedSpectrum { c: [0.0; 3] },
            },
            emission: None,
        });
        let sources = [emitter(0, 300.0, EmissionSidedness::Front)];
        let luminance = 300.0 / (PI * area);
        let camera = area_pixel(
            &scene,
            &[],
            &sources,
            PathMedium::Vacuum,
            settings(1, 8192),
            cx,
        );
        close(
            camera.luminance_cd_m2,
            luminance,
            0.002,
            "polygon uses its actual emitting area",
        );
        let probe =
            IlluminanceProbe::try_new(Point3::new(0.0, 0.0, z), Vec3::new(0.0, 0.0, -1.0)).unwrap();
        let value = area_probe(&scene, &[], &sources, probe, 32_768, cx).illuminance_lux;
        let inner = radius * (PI / f64::from(count)).cos();
        let lower = PI * luminance * inner * inner / (inner * inner + z * z);
        let upper = PI * luminance * radius * radius / (radius * radius + z * z);
        // Geometry bounds are independent of the Monte Carlo estimator. The
        // explicit 1% allowance is a numerical qualification tolerance.
        assert!(
            value >= lower * 0.99 && value <= upper * 1.01,
            "polygonal disk: measured {value}, geometric bracket [{lower}, {upper}]"
        );
    });
}

#[test]
fn g1_area_nee_and_bsdf_hits_agree_without_double_counting() {
    with_cx(|cx, _| {
        let mut scene = receiver_scene();
        scene.camera.eye.z = 0.5;
        scene.primitives.push(panel(2.0, 1.0, 0.5, false));
        let sources = [emitter(1, 100.0 * PI, EmissionSidedness::Front)];
        let expected = 0.5 * axial_rectangle_lux(50.0, 1.0, 0.5, 2.0) / PI;
        for strategy in [
            DirectStrategy::Mis,
            DirectStrategy::PowerMis,
            DirectStrategy::NeeOnly,
            DirectStrategy::BsdfOnly,
        ] {
            let mut policy = settings(2, 65_536);
            policy.strategy = strategy;
            let result = area_pixel(&scene, &[], &sources, PathMedium::Vacuum, policy, cx);
            close(
                result.luminance_cd_m2,
                expected,
                0.04,
                "independent area-light sampling strategies",
            );
        }
    });
}

#[test]
fn g1_narrowband_area_emission_retains_spectral_absorption() {
    with_cx(|cx, _| {
        let mut scene = empty_scene([0.0; 3], [0.0, 0.0, 1.0]);
        scene.primitives.push(panel(2.0, 1.0, 0.5, false));
        let density = [1.0];
        let medium = AbsorbingGrid::try_new(
            cx,
            VolumeGrid::new([1; 3], &density, [-2.0, -2.0, 0.0], [4.0, 4.0, 2.0]),
            MassAbsorptionSpectrum::table(&[(380.0, 0.1), (780.0, 1.1)]).unwrap(),
            None,
            medium_budget(),
        )
        .unwrap();
        let mut observed = Vec::new();
        for bin in [12, 60] {
            let mut bins = [0.0; SPECTRAL_BINS];
            bins[bin] = 1.0;
            let sources = [PhotometricAreaEmitter::try_new(
                0,
                100.0 * PI,
                1.0,
                SourceSpectrum::try_from_bins(bins).unwrap(),
                EmissionSidedness::Front,
            )
            .unwrap()];
            let result = area_pixel(
                &scene,
                &[],
                &sources,
                PathMedium::Absorbing(&medium),
                settings(1, 65_536),
                cx,
            );
            let kappa = 0.1 + (wavelength_nm(bin) - 380.0) / 400.0;
            close(
                result.luminance_cd_m2,
                50.0 * (-2.0 * kappa).exp(),
                0.002,
                "spectral source-hit radiance through physical absorption",
            );
            observed.push(result.luminance_cd_m2);
        }
        assert!(observed[0] > 3.0 * observed[1]);
    });
}

#[test]
fn g3_faces_of_the_emitting_primitive_can_occlude_other_faces() {
    with_cx(|cx, _| {
        let mut scene = empty_scene([0.0; 3], [0.0, 0.0, 1.0]);
        // The near square faces away from the receiver. The far square faces
        // it, but every connection is blocked by the near square. Both belong
        // to the same emitting primitive; excluding that whole object is wrong.
        scene.primitives.push(Primitive {
            shape: Shape::Mesh(TriMesh::new(
                vec![
                    [-2.0, -2.0, 1.0],
                    [2.0, -2.0, 1.0],
                    [2.0, 2.0, 1.0],
                    [-2.0, 2.0, 1.0],
                    [-1.0, -1.0, 2.0],
                    [1.0, -1.0, 2.0],
                    [1.0, 1.0, 2.0],
                    [-1.0, 1.0, 2.0],
                ],
                vec![[0, 1, 2], [0, 2, 3], [4, 6, 5], [4, 7, 6]],
            )),
            material: Material::Lambertian {
                reflectance: LiftedSpectrum { c: [0.0; 3] },
            },
            emission: None,
        });
        let sources = [emitter(0, 1000.0, EmissionSidedness::Front)];
        let probe = IlluminanceProbe::try_new(Point3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 0.0, 1.0))
            .unwrap();
        let result = area_probe(&scene, &[], &sources, probe, 4096, cx);
        assert_eq!(result.illuminance_lux, 0.0);
    });
}

#[test]
fn g1_off_axis_area_source_single_scattering_matches_volume_integral() {
    with_cx(|cx, _| {
        let mut scene = empty_scene([0.0; 3], [0.0, 0.0, 1.0]);
        // An off-axis panel is invisible to the pinhole; the measured pixel
        // receives only single-scattered light from its finite emitting area.
        scene.primitives.push(Primitive {
            shape: Shape::Mesh(TriMesh::new(
                vec![
                    [0.7, -0.3, 2.0],
                    [1.3, -0.3, 2.0],
                    [1.3, 0.3, 2.0],
                    [0.7, 0.3, 2.0],
                ],
                vec![[0, 2, 1], [0, 3, 2]],
            )),
            material: Material::Lambertian {
                reflectance: LiftedSpectrum { c: [0.0; 3] },
            },
            emission: None,
        });
        let sources = [emitter(0, 60.0 * PI * 0.36, EmissionSidedness::Front)];
        let density = [1.0];
        let medium = ScatteringGrid::try_new(
            cx,
            AbsorbingGrid::try_new(
                cx,
                VolumeGrid::new([1; 3], &density, [-4.0, -4.0, 0.0], [8.0, 8.0, 1.0]),
                MassAbsorptionSpectrum::constant(0.2).unwrap(),
                None,
                medium_budget(),
            )
            .unwrap(),
            MassScatteringSpectrum::constant(0.6).unwrap(),
            PhaseFunction::Isotropic,
        )
        .unwrap();
        let (nz, nx) = (96, 48);
        let dz = 1.0 / f64::from(nz);
        let dx = 0.6 / f64::from(nx);
        let mut integral = 0.0;
        for iz in 0..nz {
            let z = (f64::from(iz) + 0.5) * dz;
            let height = 2.0 - z;
            for ix in 0..nx {
                let x = 0.7 + (f64::from(ix) + 0.5) * dx;
                for iy in 0..nx {
                    let y = -0.3 + (f64::from(iy) + 0.5) * dx;
                    let distance = (x * x + y * y + height * height).sqrt();
                    let medium_distance = z + distance * (1.0 - z) / height;
                    integral +=
                        dz * dx * dx * height / distance.powi(3) * (-0.8 * medium_distance).exp();
                }
            }
        }
        let expected = 60.0 * 0.6 * integral / (4.0 * PI);
        let result = area_pixel(
            &scene,
            &[],
            &sources,
            PathMedium::Scattering(&medium),
            settings(1, 65_536),
            cx,
        );
        assert!(expected > 0.01);
        close(
            result.luminance_cd_m2,
            expected,
            0.04,
            "finite-area single-scatter integral with both physical extinction legs",
        );
    });
}

#[test]
fn g3_point_and_area_emission_add_linearly_in_the_shared_surface_estimator() {
    with_cx(|cx, _| {
        let mut scene = receiver_scene();
        scene.camera.eye.z = 0.5;
        scene.primitives.push(panel(2.0, 1.0, 0.5, false));
        let points = [isotropic(
            [1.5, 0.0, 1.0],
            100.0,
            SourceSpectrum::equal_energy(),
        )];
        let surfaces = [emitter(1, 100.0 * PI, EmissionSidedness::Front)];
        let policy = settings(1, 16_384);
        let point_only = area_pixel(&scene, &points, &[], PathMedium::Vacuum, policy, cx);
        let area_only = area_pixel(&scene, &[], &surfaces, PathMedium::Vacuum, policy, cx);
        let together = area_pixel(&scene, &points, &surfaces, PathMedium::Vacuum, policy, cx);
        close(
            together.luminance_cd_m2,
            point_only.luminance_cd_m2 + area_only.luminance_cd_m2,
            5.0e-13,
            "independent emission components preserve their sample streams",
        );
    });
}

#[test]
fn g1_noncoplanar_faces_of_one_emitter_illuminate_each_other() {
    with_cx(|cx, _| {
        let mut scene = empty_scene([0.0, 0.0, 0.5], [0.0, 0.0, -1.0]);
        // One primitive: a 4 m² horizontal patch and a 4 m² vertical patch.
        // Both emit; the horizontal material also reflects the vertical patch.
        scene.primitives.push(Primitive {
            shape: Shape::Mesh(TriMesh::new(
                vec![
                    [-1.0, -1.0, 0.0],
                    [1.0, -1.0, 0.0],
                    [1.0, 1.0, 0.0],
                    [-1.0, 1.0, 0.0],
                    [1.0, -1.0, 0.0],
                    [1.0, 1.0, 0.0],
                    [1.0, 1.0, 2.0],
                    [1.0, -1.0, 2.0],
                ],
                vec![[0, 1, 2], [0, 2, 3], [4, 6, 5], [4, 7, 6]],
            )),
            material: Material::Lambertian {
                reflectance: LiftedSpectrum { c: [0.0; 3] },
            },
            emission: None,
        });
        let sources = [emitter(0, 400.0 * PI, EmissionSidedness::Front)];
        let value = area_pixel(
            &scene,
            &[],
            &sources,
            PathMedium::Vacuum,
            settings(1, 65_536),
            cx,
        );
        // Integral of L*z/(1+y²+z²)² over y∈[-1,1], z∈[0,2].
        let root_five = 5.0_f64.sqrt();
        let wall_illuminance = 50.0 * (PI / 4.0 - (1.0 / root_five).atan() / root_five);
        let reflected = 0.5 * wall_illuminance / PI;
        close(
            value.luminance_cd_m2 - 50.0,
            reflected,
            0.03,
            "emission must retain reflection from noncoplanar faces of the same source",
        );
    });
}

#[test]
fn g3_g5_rigid_area_emitters_preserve_mesh_equivalence_and_parallel_replay() {
    use fs_render::instances::{GeometryInstance, SharedGeometry};
    use fs_render::tracer::RenderExecutionConfig;

    with_cx(|cx, _| {
        let half_angle = 0.2_f64;
        let placement = RigidTransform::try_new(
            [half_angle.sin(), 0.0, 0.0, half_angle.cos()],
            [0.375, -0.25, 2.0],
        )
        .unwrap();
        let Shape::Mesh(local_mesh) = panel(0.0, 1.0, 0.5, false).shape else {
            unreachable!()
        };
        let world_vertices = local_mesh
            .vertices
            .iter()
            .map(|&[x, y, z]| {
                let point = placement.transform_point(Point3::new(x, y, z));
                [point.x, point.y, point.z]
            })
            .collect();
        let explicit_mesh = TriMesh::new(world_vertices, local_mesh.triangles.clone());
        let instance = GeometryInstance::try_new(
            91,
            fs_blake3::hash_domain(
                "org.frankensim.test.photometric-area",
                b"rotated-translated-panel",
            ),
            SharedGeometry::mesh(local_mesh),
            placement,
        )
        .unwrap();
        let execution =
            RenderExecutionConfig::try_new(4, 2, 4, 64 << 20, fs_exec::RunId(491)).unwrap();
        let probes = [(0.0, 0.0), (0.3, 0.1), (-0.2, 0.2)].map(|(x, y)| {
            IlluminanceProbe::try_new(Point3::new(x, y, 0.0), Vec3::new(0.0, 0.0, 1.0)).unwrap()
        });
        let surfaces = [emitter(1, 100.0 * PI, EmissionSidedness::Front)];
        let emitters = PhotometricEmitters {
            points: &[],
            surfaces: &surfaces,
        };
        let image_policy = Settings {
            width: 3,
            height: 2,
            sampler: Sampler::Iid,
            ..settings(2, 256)
        };
        let probe_policy = IlluminanceSampling {
            samples_per_probe: 512,
            max_depth: 1,
            sampler: Sampler::Iid,
            seed: 173,
        };
        let mut explicit_image = Vec::new();
        let mut explicit_probes = Vec::new();
        for shape in [Shape::Mesh(explicit_mesh), Shape::Instance(instance)] {
            let mut scene = receiver_scene();
            scene.camera.eye.z = 0.5;
            scene.camera.half_tan = 0.25;
            scene.primitives.push(Primitive {
                shape,
                material: Material::Lambertian {
                    reflectance: LiftedSpectrum { c: [0.0; 3] },
                },
                emission: None,
            });
            let mut serial = PathLuminanceJob::new_with_emitters(
                &scene,
                emitters,
                PathMedium::Vacuum,
                TransportOptions::default(),
                image_policy,
                cx,
            )
            .unwrap();
            let mut parallel = PathLuminanceJob::new_with_emitters(
                &scene,
                emitters,
                PathMedium::Vacuum,
                TransportOptions::default(),
                image_policy,
                cx,
            )
            .unwrap();
            serial.advance(cx, usize::MAX).unwrap();
            assert_eq!(
                parallel.advance_with_execution(cx, 17, &execution).unwrap(),
                17
            );
            parallel
                .advance_with_execution(cx, usize::MAX, &execution)
                .unwrap();
            let pixel_bits = |value: PathLuminanceEstimate| {
                (
                    value.xyz.map(f64::to_bits),
                    value.luminance_cd_m2.to_bits(),
                    value.samples,
                    value.estimated_standard_error_cd_m2.map(f64::to_bits),
                    value.total_error_bound_cd_m2.map(f64::to_bits),
                )
            };
            assert!(serial.is_complete() && parallel.is_complete());
            assert_eq!(
                serial
                    .estimates()
                    .map(|value| pixel_bits(value.unwrap()))
                    .collect::<Vec<_>>(),
                parallel
                    .estimates()
                    .map(|value| pixel_bits(value.unwrap()))
                    .collect::<Vec<_>>()
            );
            let image: Vec<_> = serial
                .estimates()
                .map(|value| value.unwrap().luminance_cd_m2)
                .collect();
            assert!(image.iter().all(|&value| value > 0.0));
            if explicit_image.is_empty() {
                explicit_image = image;
            } else {
                for (&observed, &expected) in image.iter().zip(&explicit_image) {
                    close(
                        observed,
                        expected,
                        1.0e-10,
                        "placed area emitter versus explicit mesh",
                    );
                }
            }

            let mut serial = PathIlluminanceJob::new_with_emitters(
                &scene,
                emitters,
                PathMedium::Vacuum,
                &probes,
                TransportOptions::default(),
                probe_policy,
                cx,
            )
            .unwrap();
            let mut parallel = PathIlluminanceJob::new_with_emitters(
                &scene,
                emitters,
                PathMedium::Vacuum,
                &probes,
                TransportOptions::default(),
                probe_policy,
                cx,
            )
            .unwrap();
            serial.advance(cx, usize::MAX).unwrap();
            assert_eq!(
                parallel.advance_with_execution(cx, 11, &execution).unwrap(),
                11
            );
            parallel
                .advance_with_execution(cx, usize::MAX, &execution)
                .unwrap();
            let probe_bits = |value: PathIlluminanceEstimate| {
                (
                    [
                        value.direct_lux,
                        value.continuous_lux,
                        value.illuminance_lux,
                    ]
                    .map(f64::to_bits),
                    value.samples,
                    value.estimated_standard_error_lux.map(f64::to_bits),
                    value.total_error_bound_lux.map(f64::to_bits),
                )
            };
            assert!(serial.is_complete() && parallel.is_complete());
            assert_eq!(
                serial
                    .estimates()
                    .map(|value| probe_bits(value.unwrap()))
                    .collect::<Vec<_>>(),
                parallel
                    .estimates()
                    .map(|value| probe_bits(value.unwrap()))
                    .collect::<Vec<_>>()
            );
            let values: Vec<_> = serial
                .estimates()
                .map(|value| {
                    let value = value.unwrap();
                    assert_eq!(value.direct_lux, 0.0);
                    assert!(value.continuous_lux > 0.0);
                    value.illuminance_lux
                })
                .collect();
            if explicit_probes.is_empty() {
                explicit_probes = values;
            } else {
                for (&observed, &expected) in values.iter().zip(&explicit_probes) {
                    close(
                        observed,
                        expected,
                        1.0e-10,
                        "placed area probe versus explicit mesh",
                    );
                }
            }
        }
    });
}
