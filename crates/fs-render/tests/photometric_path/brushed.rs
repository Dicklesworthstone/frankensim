//! G1/G3 calibrated anisotropic reflection, using independent angular integrals.

use super::*;
use fs_render::conductor::{
    ConductorDataStatus, ConductorIorSample, ConductorOptics, ConductorSource,
};
use fs_render::instances::{GeometryInstance, SharedGeometry};
use fs_render::photometry::area::{EmissionSidedness, PhotometricAreaEmitter, PhotometricEmitters};
use fs_render::surface_detail::{BrushingPattern, SurfaceDetail};

fn optics() -> ConductorOptics {
    // A declared synthetic, spectrally constant n=2, k=3 conductor. This makes
    // the independent angular reference independent of wavelength quadrature.
    ConductorOptics::try_new(
        std::array::from_fn(|i| {
            ConductorIorSample::try_new(380.0 + 50.0 * i as f64, 2.0, 3.0).unwrap()
        }),
        ConductorSource::try_new(
            fs_blake3::hash_domain("org.frankensim.test.photometric-brushed", b"n2-k3"),
            ConductorDataStatus::Representative,
        )
        .unwrap(),
    )
    .unwrap()
}

fn brushed_scene(alpha: f64, ratio: f64) -> Scene {
    let mut scene = empty_scene([1.0, 0.0, 1.0], [0.0, 0.0, -1.0]);
    let mut receiver = floor(0.0, 10.0, LiftedSpectrum { c: [0.0; 3] });
    receiver.material = Material::BrushedConductor {
        optics: optics(),
        detail: SurfaceDetail::try_new(alpha, ratio, BrushingPattern::Circular, 0.0, None).unwrap(),
    };
    scene.primitives.push(receiver);
    scene
}

fn brdf_normal_view(wx: f64, wy: f64, wz: f64, alpha_x: f64, alpha_y: f64) -> f64 {
    // Tangent-space Trowbridge-Reitz distribution and Smith G1 at wo=(0,0,1).
    // This angular reference does not call the renderer's BSDF or sampler.
    // https://pbr-book.org/4ed/Reflection_Models/Roughness_Using_Microfacet_Theory
    let half_length = (2.0 + 2.0 * wz).sqrt();
    let hx = wx / half_length;
    let hy = wy / half_length;
    let hz = (1.0 + wz) / half_length;
    let ellipse = (hx / alpha_x).powi(2) + (hy / alpha_y).powi(2) + hz * hz;
    let distribution = 1.0 / (PI * alpha_x * alpha_y * ellipse * ellipse);
    let smith =
        2.0 / (1.0 + (1.0 + ((alpha_x * wx).powi(2) + (alpha_y * wy).powi(2)) / (wz * wz)).sqrt());
    distribution * smith * fresnel_n2_k3_reference(hz) / (4.0 * wz)
}

fn fresnel_n2_k3_reference(cosine: f64) -> f64 {
    // Real-valued unpolarized complex-index Fresnel, independently expanded
    // from n=2,k=3. The production path uses complex arithmetic instead.
    let sine2 = 1.0 - cosine * cosine;
    let t0 = 4.0 - 9.0 - sine2;
    let a2_b2 = (t0 * t0 + 4.0 * 4.0 * 9.0).sqrt();
    let a = (0.5 * (a2_b2 + t0)).sqrt();
    let t1 = a2_b2 + cosine * cosine;
    let t2 = 2.0 * cosine * a;
    let rs = (t1 - t2) / (t1 + t2);
    let t3 = cosine * cosine * a2_b2 + sine2 * sine2;
    let t4 = t2 * sine2;
    let rp = rs * (t3 - t4) / (t3 + t4);
    0.5 * (rs + rp)
}

fn tilted_radiance_bsdf_reference(wo: [f64; 3], wi: [f64; 3], angle: f64) -> f64 {
    // Radial grain at the +X axis: T=(1,0,0). An explicit rotation about
    // that fixed axis gives B=(0,cos(theta),-sin(theta)) and
    // Ns=(0,sin(theta),cos(theta)), without using production frame helpers.
    let (sine, cosine) = angle.sin_cos();
    let to_local = |v: [f64; 3]| {
        [
            v[0],
            cosine * v[1] - sine * v[2],
            sine * v[1] + cosine * v[2],
        ]
    };
    let local_o = to_local(wo);
    let local_i = to_local(wi);
    if wo[2] <= 0.0 || wi[2] <= 0.0 || local_o[2] <= 0.0 || local_i[2] <= 0.0 {
        return 0.0;
    }
    let sum = std::array::from_fn::<_, 3, _>(|axis| wo[axis] + wi[axis]);
    let length = sum.iter().map(|value| value * value).sum::<f64>().sqrt();
    let half = sum.map(|value| value / length);
    let local_half = to_local(half);
    let ellipse =
        (local_half[0] / 0.12).powi(2) + (local_half[1] / 0.36).powi(2) + local_half[2].powi(2);
    let distribution = 1.0 / (PI * 0.12 * 0.36 * ellipse * ellipse);
    let masking = |v: [f64; 3]| {
        2.0 / (1.0 + (1.0 + ((0.12 * v[0]).powi(2) + (0.36 * v[1]).powi(2)) / v[2].powi(2)).sqrt())
    };
    let half_cosine = wo
        .iter()
        .zip(half)
        .map(|(left, right)| left * right)
        .sum::<f64>();
    let shading_brdf = distribution
        * masking(local_o)
        * masking(local_i)
        * fresnel_n2_k3_reference(half_cosine.clamp(0.0, 1.0))
        / (4.0 * local_o[2] * local_i[2]);
    // The radiance kernel uses geometric incident lux. Convert the shading
    // projected measure to that measure: f_R = f_s * (Ns.wi)/(Ng.wi).
    // Only f_s is reciprocal; f_R must not be compared by an unweighted swap.
    // This decorative normal model does not resolve measured microgeometry.
    shading_brdf * local_i[2] / wi[2]
}

#[test]
fn g1_brushed_point_luminance_retains_anisotropy_and_body_frame() {
    with_cx(|cx, _| {
        let placement = RigidTransform::try_new(
            [
                0.0,
                0.0,
                core::f64::consts::FRAC_1_SQRT_2,
                core::f64::consts::FRAC_1_SQRT_2,
            ],
            [7.0, -3.0, 4.0],
        )
        .unwrap();
        let mut local_values = Vec::new();
        for placed in [false, true] {
            let mut scene = brushed_scene(0.12, 3.0);
            if placed {
                let Shape::Mesh(mesh) = scene.primitives.remove(0).shape else {
                    unreachable!()
                };
                let material = brushed_scene(0.12, 3.0).primitives[0].material;
                scene.primitives.push(Primitive {
                    shape: Shape::Instance(
                        GeometryInstance::try_new(
                            211,
                            fs_blake3::hash_domain(
                                "org.frankensim.test.photometric-brushed",
                                b"placed-floor",
                            ),
                            SharedGeometry::mesh(mesh),
                            placement,
                        )
                        .unwrap(),
                    ),
                    material,
                    emission: None,
                });
                scene.camera.eye = placement.transform_point(scene.camera.eye);
                scene.camera.forward = placement.transform_vector(scene.camera.forward);
                scene.camera.up = placement.transform_vector(scene.camera.up);
            }
            // At local (1,0,0), circular grooves run along +Y and across -X.
            for (case, (position, wx, wy)) in
                [([1.0, 1.5, 2.0], 0.6, 0.0), ([2.5, 0.0, 2.0], 0.0, -0.6)]
                    .into_iter()
                    .enumerate()
            {
                let mut position = Point3::new(position[0], position[1], position[2]);
                if placed {
                    position = placement.transform_point(position);
                }
                let sources = [isotropic(
                    [position.x, position.y, position.z],
                    100.0,
                    SourceSpectrum::equal_energy(),
                )];
                let result =
                    render_pixel(&scene, &sources, PathMedium::Vacuum, settings(1, 8192), cx);
                let expected =
                    (100.0 / 2.5_f64.powi(2)) * 0.8 * brdf_normal_view(wx, wy, 0.8, 0.12, 0.36);
                close(
                    result.luminance_cd_m2,
                    expected,
                    0.002,
                    "anisotropic absolute point-source luminance",
                );
                if placed {
                    close(
                        result.luminance_cd_m2,
                        local_values[case],
                        1.0e-11,
                        "body-locked brushing under rigid placement",
                    );
                } else {
                    local_values.push(result.luminance_cd_m2);
                }
            }
        }
        assert!(
            local_values[1] > 10.0 * local_values[0],
            "isotropic fallback must not replace the brushed frame"
        );
    });
}

#[test]
fn g1_tilted_brushed_point_luminance_matches_independent_radiance_reference() {
    with_cx(|cx, _| {
        // Flatten photopic weight using the independent published CIE table,
        // so this constant-n,k test measures geometry/BRDF without wavelength
        // sampling variation. The source still enters through native lumens.
        let spectrum = SourceSpectrum::try_from_bins(std::array::from_fn(|bin| {
            1.0 / super::spectral_reference::cie_reference_xyz(382.5 + 5.0 * bin as f64)[1]
        }))
        .unwrap();
        let target = [0.5, 0.0, 0.0];
        for view_cosine in [0.1_f64, 0.02] {
            let view = [0.0, (1.0 - view_cosine * view_cosine).sqrt(), view_cosine];
            let light = [0.0, -0.6, 0.8];
            let mut recovered_shading_brdfs = Vec::new();
            for (wo, wi) in [(view, light), (light, view)] {
                let mut scene = empty_scene(
                    std::array::from_fn(|axis| target[axis] + 2.0 * wo[axis]),
                    wo.map(|component| -component),
                );
                // View directions lie in YZ, so +X is the orthogonal camera up.
                scene.camera.up = Vec3::new(1.0, 0.0, 0.0);
                let mut receiver = floor(0.0, 10.0, LiftedSpectrum { c: [0.0; 3] });
                receiver.material = Material::BrushedConductor {
                    optics: optics(),
                    detail: SurfaceDetail::try_new(0.12, 3.0, BrushingPattern::Radial, 0.35, None)
                        .unwrap(),
                };
                scene.primitives.push(receiver);
                let sources = [isotropic(
                    std::array::from_fn(|axis| target[axis] + 3.0 * wi[axis]),
                    100.0,
                    spectrum.clone(),
                )];
                let incident_lux = (100.0 / 9.0) * wi[2];
                let expected_bsdf = tilted_radiance_bsdf_reference(wo, wi, 0.35);
                let measured =
                    render_pixel(&scene, &sources, PathMedium::Vacuum, settings(1, 64), cx);
                close(
                    measured.luminance_cd_m2,
                    incident_lux * expected_bsdf,
                    1.0e-9,
                    "tilted normal uses geometric incident lux and the shading radiance projection",
                );
                assert_eq!(measured.numerical.kind, fs_evidence::NumericalKind::NoClaim);
                assert!(measured.total_error_bound_cd_m2.is_none());
                let shading_incident_cosine = 0.35_f64.sin() * wi[1] + 0.35_f64.cos() * wi[2];
                recovered_shading_brdfs.push(
                    measured.luminance_cd_m2 / incident_lux * wi[2] / shading_incident_cosine,
                );
            }
            close(
                recovered_shading_brdfs[0],
                recovered_shading_brdfs[1],
                2.0e-9,
                "removing the radiance projection recovers the reciprocal local GGX law",
            );
        }
    });
}

fn panel_reference(side: usize) -> f64 {
    let mut sum = 0.0;
    for y in 0..side {
        for x in 0..side {
            let dx = 0.3 + (x as f64 + 0.5) / side as f64;
            let dy = -0.4 + 0.8 * (y as f64 + 0.5) / side as f64;
            let r2 = dx * dx + dy * dy + 4.0;
            let r = r2.sqrt();
            let cosine = 2.0 / r;
            sum +=
                50.0 * brdf_normal_view(dy / r, -dx / r, cosine, 0.2, 0.5) * cosine * cosine / r2;
        }
    }
    sum * 0.8 / (side * side) as f64
}

#[test]
fn g1_brushed_area_strategies_match_independent_irradiance_integral() {
    with_cx(|cx, _| {
        let mut scene = brushed_scene(0.2, 2.5);
        scene.primitives.push(Primitive {
            shape: Shape::Mesh(TriMesh::new(
                vec![
                    [1.3, -0.4, 2.0],
                    [2.3, -0.4, 2.0],
                    [2.3, 0.4, 2.0],
                    [1.3, 0.4, 2.0],
                ],
                vec![[0, 2, 1], [0, 3, 2]],
            )),
            material: Material::Lambertian {
                reflectance: LiftedSpectrum { c: [0.0; 3] },
            },
            emission: None,
        });
        let surfaces = [PhotometricAreaEmitter::try_new(
            1,
            40.0 * PI,
            1.0,
            SourceSpectrum::equal_energy(),
            EmissionSidedness::Front,
        )
        .unwrap()];
        let expected = panel_reference(128);
        close(
            panel_reference(64),
            expected,
            0.0002,
            "reference panel quadrature convergence",
        );
        for strategy in [
            DirectStrategy::NeeOnly,
            DirectStrategy::Mis,
            DirectStrategy::PowerMis,
            DirectStrategy::BsdfOnly,
        ] {
            let mut job = PathLuminanceJob::new_with_emitters(
                &scene,
                PhotometricEmitters {
                    points: &[],
                    surfaces: &surfaces,
                },
                PathMedium::Vacuum,
                TransportOptions::default(),
                Settings {
                    strategy,
                    sampler: Sampler::Iid,
                    ..settings(2, 32768)
                },
                cx,
            )
            .unwrap();
            job.advance(cx, usize::MAX).unwrap();
            let result = job.pixel_estimate(0).unwrap();
            let error = result.estimated_standard_error_cd_m2.unwrap();
            assert!(
                (result.luminance_cd_m2 - expected).abs() < 6.0 * error + 0.002 * expected,
                "{strategy:?}: luminance={}, reference={expected}, estimated SE={error}",
                result.luminance_cd_m2
            );
            assert!(result.total_error_bound_cd_m2.is_none());
        }
    });
}
