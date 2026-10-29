//! G1/G3: physical thermal spectra reach actual tracer pixels and reflection.
use asupersync::types::Budget;
use fs_blake3::ContentHash;
use fs_exec::{CancelGate, Cx, ExecMode, StreamKey};
use fs_geom::{Point3, Vec3};
use fs_render::charts::TriMesh;
use fs_render::dielectric::{DielectricGlass, DielectricSurface};
use fs_render::lighting::EnvironmentMap;
use fs_render::spectral::{lift_rgb, xyz_of_spectrum};
use fs_render::thermal_surface::GrayDiffuseThermalSurface;
use fs_render::tracer::{
    Camera, DirectStrategy, Material, Primitive, Sampler, Scene, Settings, Shape, TracerError,
    render,
};

#[path = "thermal_surface/lighting.rs"]
mod lighting;

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        f(&Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 40,
                kernel_id: 1,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}

fn surface(t: f64, e: f64) -> GrayDiffuseThermalSurface {
    GrayDiffuseThermalSurface::try_new(t, e, ContentHash([40; 32])).unwrap()
}

fn scene(material: Material) -> Scene {
    Scene {
        primitives: vec![Primitive {
            shape: Shape::Mesh(TriMesh::new(
                vec![
                    [-10.0, -10.0, 0.0],
                    [10.0, -10.0, 0.0],
                    [10.0, 10.0, 0.0],
                    [-10.0, 10.0, 0.0],
                ],
                vec![[0, 1, 2], [0, 2, 3]],
            )),
            material,
            emission: None,
        }],
        lights: Vec::new(),
        environment: None,
        camera: Camera {
            eye: Point3::new(0.0, 0.0, 1.0),
            forward: Vec3::new(0.0, 0.0, -1.0),
            up: Vec3::new(0.0, 1.0, 0.0),
            half_tan: 0.0,
        },
    }
}

fn settings(strategy: DirectStrategy) -> Settings {
    Settings {
        width: 1,
        height: 1,
        spp: 4096,
        max_depth: 2,
        sampler: Sampler::OwenSobol,
        strategy,
        seed: 40,
    }
}

#[test]
fn g1_hot_surface_lights_a_dark_scene_at_its_physical_radiance() {
    with_cx(|cx| {
        let body = surface(1800.0, 0.4);
        let expected = xyz_of_spectrum(|l| body.radiance_w_m2_sr_nm(l).unwrap());
        for strategy in [
            DirectStrategy::BsdfOnly,
            DirectStrategy::NeeOnly,
            DirectStrategy::Mis,
            DirectStrategy::PowerMis,
        ] {
            let image = render(
                &scene(Material::GrayDiffuseThermal { surface: body }),
                cx,
                &settings(strategy),
            )
            .unwrap();
            for (got, want) in image.xyz[0].iter().zip(expected) {
                let mean = got / f64::from(image.spp_done);
                assert!(
                    (mean / want - 1.0).abs() < 0.005,
                    "{strategy:?}: {mean} vs {want}"
                );
            }
        }
        let hotter = surface(2200.0, 0.4);
        let a = render(
            &scene(Material::GrayDiffuseThermal { surface: body }),
            cx,
            &settings(DirectStrategy::Mis),
        )
        .unwrap();
        let b = render(
            &scene(Material::GrayDiffuseThermal { surface: hotter }),
            cx,
            &settings(DirectStrategy::Mis),
        )
        .unwrap();
        assert!(b.xyz[0][1] > a.xyz[0][1]);
        assert!(b.xyz[0][2] / b.xyz[0][0] > a.xyz[0][2] / a.xyz[0][0]);
    });
}

#[test]
fn g3_thermal_surfaces_keep_reciprocal_reflection_and_have_one_emission_owner() {
    with_cx(|cx| {
        // At 1 K the visible emission underflows to zero. The original camera
        // path must still reflect the exterior field with exactly 1-epsilon.
        let make = |e| {
            let mut s = scene(Material::GrayDiffuseThermal {
                surface: surface(1.0, e),
            });
            s.environment =
                Some(EnvironmentMap::try_from_linear_srgb(1, 1, vec![[1.0; 3]], 0.0).unwrap());
            s
        };
        let a = render(&make(0.25), cx, &settings(DirectStrategy::BsdfOnly)).unwrap();
        let b = render(&make(0.75), cx, &settings(DirectStrategy::BsdfOnly)).unwrap();
        for channel in 0..3 {
            assert!(a.xyz[0][channel] > 0.0);
            assert!((a.xyz[0][channel] / b.xyz[0][channel] - 3.0).abs() < 1.0e-10);
        }
        let mut bad = make(0.5);
        bad.primitives[0].emission = Some((lift_rgb([1.0; 3]), 1.0));
        assert_eq!(
            render(&bad, cx, &settings(DirectStrategy::Mis)).unwrap_err(),
            TracerError::DuplicateThermalEmission
        );
    });
}

#[cfg(feature = "material-state")]
#[path = "thermal_surface/material.rs"]
mod material;

#[test]
fn g0_thermal_body_can_be_viewed_through_glass_but_not_embedded_in_it() {
    with_cx(|cx| {
        let with_glass = |bottom: f64| {
            let mut s = scene(Material::GrayDiffuseThermal {
                surface: surface(1800.0, 0.4),
            });
            s.primitives.push(Primitive {
                shape: Shape::Mesh(TriMesh::new(
                    vec![
                        [-20.0, -20.0, bottom],
                        [20.0, -20.0, bottom],
                        [20.0, 20.0, bottom],
                        [-20.0, 20.0, bottom],
                        [-20.0, -20.0, 0.5],
                        [20.0, -20.0, 0.5],
                        [20.0, 20.0, 0.5],
                        [-20.0, 20.0, 0.5],
                    ],
                    vec![
                        [0, 2, 1],
                        [0, 3, 2],
                        [4, 5, 6],
                        [4, 6, 7],
                        [0, 1, 5],
                        [0, 5, 4],
                        [1, 2, 6],
                        [1, 6, 5],
                        [2, 3, 7],
                        [2, 7, 6],
                        [3, 0, 4],
                        [3, 4, 7],
                    ],
                )),
                material: Material::Dielectric {
                    glass: DielectricGlass::representative_borosilicate(),
                    surface: DielectricSurface::SMOOTH,
                },
                emission: None,
            });
            s
        };
        let settings = Settings {
            spp: 64,
            max_depth: 4,
            ..settings(DirectStrategy::BsdfOnly)
        };
        let exterior = render(&with_glass(0.25), cx, &settings).unwrap();
        assert!(exterior.xyz[0].iter().all(|x| x.is_finite() && *x > 0.0));
        assert_eq!(
            render(&with_glass(-0.5), cx, &settings).unwrap_err(),
            TracerError::ThermalSurfaceMediumUnsupported
        );
    });
}
