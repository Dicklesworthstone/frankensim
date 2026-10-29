//! G1/G3: physical hot mesh faces illuminate the original surface integrator.

use super::*;
use fs_render::tracer::{Film, render_range};

fn quad(z: f64, x: f64, half_width: f64, material: Material) -> Primitive {
    Primitive {
        shape: Shape::Mesh(TriMesh::new(
            vec![
                [x - half_width, -half_width, z],
                [x + half_width, -half_width, z],
                [x + half_width, half_width, z],
                [x - half_width, half_width, z],
            ],
            vec![[0, 1, 2], [0, 2, 3]],
        )),
        material,
        emission: None,
    }
}

fn illuminated_scene(temperature_k: f64, emissivity: f64) -> Scene {
    let mut s = scene(Material::Lambertian {
        reflectance: lift_rgb([0.4; 3]),
    });
    s.camera.eye.z = 0.5;
    s.primitives.push(quad(
        2.0,
        0.0,
        1.0,
        Material::GrayDiffuseThermal {
            surface: surface(temperature_k, emissivity),
        },
    ));
    s
}

// Independent midpoint area quadrature of cos(theta_s) cos(theta_l) / r².
// The receiver is at the origin and the square occupies [-1,1]² at z=2.
fn projected_solid_angle() -> f64 {
    let n = 256;
    let dx = 2.0 / f64::from(n);
    let mut result = 0.0;
    for i in 0..n {
        let x = -1.0 + (f64::from(i) + 0.5) * dx;
        for j in 0..n {
            let y = -1.0 + (f64::from(j) + 0.5) * dx;
            let r2 = x * x + y * y + 4.0;
            result += 4.0 / (r2 * r2) * dx * dx;
        }
    }
    result
}

#[test]
fn g1_hot_mesh_illumination_matches_area_quadrature_across_estimators() {
    with_cx(|cx| {
        let body = surface(2200.0, 0.7);
        let reflectance = lift_rgb([0.4; 3]);
        let geometry = projected_solid_angle() / core::f64::consts::PI;
        let expected = xyz_of_spectrum(|lambda| {
            reflectance.eval(lambda) * body.radiance_w_m2_sr_nm(lambda).unwrap() * geometry
        });
        for strategy in [
            DirectStrategy::NeeOnly,
            DirectStrategy::BsdfOnly,
            DirectStrategy::Mis,
            DirectStrategy::PowerMis,
        ] {
            let mut settings = settings(strategy);
            settings.sampler = Sampler::OwenSobolFullPath;
            settings.spp = 16_384;
            let image = render(&illuminated_scene(2200.0, 0.7), cx, &settings).unwrap();
            for (channel, want) in expected.into_iter().enumerate() {
                let got = image.xyz[0][channel] / f64::from(image.spp_done);
                assert!(
                    (got / want - 1.0).abs() < 0.015,
                    "{strategy:?}: {got} vs {want}"
                );
            }
        }
    });
}

#[test]
fn g1_final_bounce_keeps_complete_thermal_next_event_weight() {
    with_cx(|cx| {
        let mut settings = settings(DirectStrategy::Mis);
        settings.sampler = Sampler::OwenSobolFullPath;
        settings.spp = 8192;
        settings.max_depth = 1;
        let image = render(&illuminated_scene(2200.0, 1.0), cx, &settings).unwrap();
        let body = surface(2200.0, 1.0);
        let reflectance = lift_rgb([0.4; 3]);
        let geometry = projected_solid_angle() / core::f64::consts::PI;
        let expected = xyz_of_spectrum(|lambda| {
            reflectance.eval(lambda) * body.radiance_w_m2_sr_nm(lambda).unwrap() * geometry
        });
        for (channel, want) in expected.into_iter().enumerate() {
            let got = image.xyz[0][channel] / f64::from(image.spp_done);
            assert!((got / want - 1.0).abs() < 0.01, "{got} vs {want}");
        }
    });
}

#[test]
fn g3_thermal_light_obeys_occlusion_and_physical_emissivity() {
    with_cx(|cx| {
        let mut settings = settings(DirectStrategy::NeeOnly);
        settings.spp = 1024;
        settings.max_depth = 1;
        let high = render(&illuminated_scene(2200.0, 0.8), cx, &settings).unwrap();
        let low = render(&illuminated_scene(2200.0, 0.2), cx, &settings).unwrap();
        for channel in 0..3 {
            assert!(high.xyz[0][channel] > 0.0);
            assert!((high.xyz[0][channel] / low.xyz[0][channel] - 4.0).abs() < 1.0e-11);
        }
        let mut blocked = illuminated_scene(2200.0, 0.8);
        blocked.primitives.push(quad(
            1.0,
            0.0,
            2.0,
            Material::Lambertian {
                reflectance: lift_rgb([0.0; 3]),
            },
        ));
        let image = render(&blocked, cx, &settings).unwrap();
        assert_eq!(image.xyz, vec![[0.0; 3]]);
    });
}

#[test]
fn g3_same_thermal_body_occludes_faces_separated_by_micrometres() {
    with_cx(|cx| {
        let mut s = illuminated_scene(2200.0, 1.0);
        // Both layers belong to one primitive. The rear layer is closer than
        // the former 8*RAY_EPS visibility allowance, but remains fully hidden.
        s.primitives[1].shape = Shape::Mesh(TriMesh::new(
            vec![
                [-1.0, -1.0, 2.0],
                [1.0, -1.0, 2.0],
                [1.0, 1.0, 2.0],
                [-1.0, 1.0, 2.0],
                [-1.0, -1.0, 2.0 + 4.0e-6],
                [1.0, -1.0, 2.0 + 4.0e-6],
                [1.0, 1.0, 2.0 + 4.0e-6],
                [-1.0, 1.0, 2.0 + 4.0e-6],
            ],
            vec![[0, 1, 2], [0, 2, 3], [4, 5, 6], [4, 6, 7]],
        ));
        let mut settings = settings(DirectStrategy::NeeOnly);
        settings.sampler = Sampler::OwenSobolFullPath;
        settings.spp = 8192;
        settings.max_depth = 1;
        let image = render(&s, cx, &settings).unwrap();
        let body = surface(2200.0, 1.0);
        let reflectance = lift_rgb([0.4; 3]);
        let geometry = projected_solid_angle() / core::f64::consts::PI;
        let expected = xyz_of_spectrum(|lambda| {
            reflectance.eval(lambda) * body.radiance_w_m2_sr_nm(lambda).unwrap() * geometry
        });
        for (channel, want) in expected.into_iter().enumerate() {
            let got = image.xyz[0][channel] / f64::from(image.spp_done);
            assert!((got / want - 1.0).abs() < 0.015, "{got} vs {want}");
        }
    });
}

#[test]
fn g5_thermal_lights_replay_across_primitive_order_and_progressive_ranges() {
    with_cx(|cx| {
        let mut s = illuminated_scene(2200.0, 0.7);
        s.primitives.push(quad(
            3.0,
            2.0,
            0.4,
            Material::GrayDiffuseThermal {
                surface: surface(3100.0, 0.3),
            },
        ));
        let mut settings = settings(DirectStrategy::PowerMis);
        settings.sampler = Sampler::OwenSobolFullPath;
        settings.spp = 256;
        let original = render(&s, cx, &settings).unwrap();
        s.primitives.reverse();
        let reordered = render(&s, cx, &settings).unwrap();
        assert_eq!(original, reordered);
        let mut resumed = Film::new(1, 1);
        render_range(&s, cx, &settings, &mut resumed, 0, 73).unwrap();
        render_range(&s, cx, &settings, &mut resumed, 73, 256).unwrap();
        assert_eq!(original, resumed);
    });
}
