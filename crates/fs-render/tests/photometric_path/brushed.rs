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
    // Real-valued unpolarized complex-index Fresnel, independently expanded
    // from n=2,k=3. The production path uses complex arithmetic instead.
    let sine2 = 1.0 - hz * hz;
    let t0 = 4.0 - 9.0 - sine2;
    let a2_b2 = (t0 * t0 + 4.0 * 4.0 * 9.0).sqrt();
    let a = (0.5 * (a2_b2 + t0)).sqrt();
    let t1 = a2_b2 + hz * hz;
    let t2 = 2.0 * hz * a;
    let rs = (t1 - t2) / (t1 + t2);
    let t3 = hz * hz * a2_b2 + sine2 * sine2;
    let t4 = t2 * sine2;
    let rp = rs * (t3 - t4) / (t3 + t4);
    distribution * smith * 0.5 * (rs + rp) / (4.0 * wz)
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
