//! Physical ambient fields coexist with closed, volume-displacing glass.

use super::*;
use fs_render::dielectric::{BeerLambertAbsorption, CauchyIor, GlassProvenance};
use fs_render::tracer::{ScatteringFilm, render_scattering, render_scattering_range};
use fs_render::volumes::transport::scattering::{
    MassScatteringSpectrum, PhaseFunction, ScatteringGrid,
};

fn glass(ior: CauchyIor, extinction: f64) -> Material {
    Material::Dielectric {
        glass: DielectricGlass::new(
            ior,
            BeerLambertAbsorption::try_constant(extinction).unwrap(),
            GlassProvenance::Custom,
        ),
        surface: DielectricSurface::SMOOTH,
    }
}

fn slab(material: Material) -> Primitive {
    // A complete outward-oriented closed box; its large side faces do not
    // intersect the normal/oblique direct transmission fixtures.
    Primitive {
        shape: Shape::Mesh(TriMesh::new(
            vec![
                [-4.0, -4.0, -0.2],
                [4.0, -4.0, -0.2],
                [4.0, 4.0, -0.2],
                [-4.0, 4.0, -0.2],
                [-4.0, -4.0, 0.0],
                [4.0, -4.0, 0.0],
                [4.0, 4.0, 0.0],
                [-4.0, 4.0, 0.0],
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
        material,
        emission: None,
    }
}

fn transmission_scene(material: Material, cosine: f64) -> (Scene, (LiftedSpectrum, f64)) {
    let mut scene = empty_scene();
    scene.camera.eye = Point3::new(0.0, 0.0, 1.0);
    scene.camera.forward = Vec3::new((1.0 - cosine * cosine).sqrt(), 0.0, -cosine);
    scene.primitives.push(slab(material));
    let emission = (lift_rgb([1.0; 3]), 3.0);
    let mut light = quad(
        -1.0,
        4.0,
        Material::Lambertian {
            reflectance: lift_rgb([0.0; 3]),
        },
    );
    light.emission = Some(emission);
    scene.primitives.push(light);
    scene.lights.push(RectLight {
        corner: Point3::new(-4.0, -4.0, -1.0),
        edge_u: Vec3::new(8.0, 0.0, 0.0),
        edge_v: Vec3::new(0.0, 8.0, 0.0),
        prim: 1,
        emission,
    });
    (scene, emission)
}

fn ambient<'a>(
    cx: &Cx<'_>,
    density: &'a [f64],
    absorption: MassAbsorptionSpectrum,
) -> AbsorbingGrid<'a> {
    AbsorbingGrid::try_new(
        cx,
        VolumeGrid::new([1; 3], density, [-8.0, -8.0, -1.0], [16.0, 16.0, 2.0]),
        absorption,
        None,
        budget(),
    )
    .unwrap()
}

// Independent unpolarized Fresnel formula and Snell cosine. The parallel
// exit has identical reflectance and cancels the two radiance eta² factors.
fn interface(index: f64, cosine: f64) -> (f64, f64) {
    let transmitted_cosine = (1.0 - (1.0 - cosine * cosine) / (index * index)).sqrt();
    let rs = (cosine - index * transmitted_cosine) / (cosine + index * transmitted_cosine);
    let rp = (index * cosine - transmitted_cosine) / (index * cosine + transmitted_cosine);
    (0.5 * (rs * rs + rp * rp), transmitted_cosine)
}

#[test]
fn g1_dispersive_refraction_preserves_ambient_transport_on_every_split_suffix() {
    with_cx(|cx, _| {
        let density = [1.0];
        let absorption =
            MassAbsorptionSpectrum::table(&[(380.0, 0.1), (580.0, 0.55), (780.0, 1.0)]).unwrap();
        let medium = ambient(cx, &density, absorption);
        let ior = CauchyIor::try_new(1.4, 0.012, 0.0).unwrap();
        for cosine in [1.0, 0.8] {
            let (scene, (emission, scale)) = transmission_scene(glass(ior, 0.6), cosine);
            let mut settings = settings(DirectStrategy::Mis);
            settings.max_depth = 3;
            settings.spp = 4096;
            let expected = xyz_of_spectrum(|lambda| {
                let um = lambda * 1.0e-3;
                let n = 1.4 + 0.012 / (um * um);
                let (fresnel, transmitted_cosine) = interface(n, cosine);
                let ambient_length = 1.8 / cosine;
                emission.eval(lambda)
                    * scale
                    * (1.0 - fresnel)
                    * (1.0 - fresnel)
                    * (-0.6 * 0.2 / transmitted_cosine).exp()
                    * (-absorption.at(lambda).unwrap() * ambient_length).exp()
            });
            let reference = render_absorbing(&scene, &medium, cx, &settings).unwrap();
            assert_mean(&reference, expected, 0.006);
            for strategy in [
                DirectStrategy::BsdfOnly,
                DirectStrategy::NeeOnly,
                DirectStrategy::PowerMis,
            ] {
                settings.strategy = strategy;
                assert_eq!(
                    render_absorbing(&scene, &medium, cx, &settings).unwrap(),
                    reference,
                    "delta transmission must not acquire an unevaluated manifold competitor"
                );
            }
        }
    });
}

#[test]
fn g1_forced_fresnel_reflection_and_transmission_both_keep_the_ambient_field() {
    with_cx(|cx, _| {
        let density = [1.0];
        let coefficient = 0.7;
        let medium = ambient(
            cx,
            &density,
            MassAbsorptionSpectrum::constant(coefficient).unwrap(),
        );
        let (mut scene, (emission, scale)) = environment_scene();
        scene.camera.eye.z = 1.0;
        scene
            .primitives
            .push(slab(glass(CauchyIor::try_constant(1.5).unwrap(), 0.0)));
        let mut settings = settings(DirectStrategy::Mis);
        settings.max_depth = 3;
        settings.spp = 2048;
        // This depth includes primary reflection and two-interface transmission;
        // internally reflected suffixes have not reached an emitter yet.
        let (r, _) = interface(1.5, 1.0);
        let response =
            r * (-2.0 * coefficient).exp() + (1.0 - r) * (1.0 - r) * (-1.8 * coefficient).exp();
        assert_mean(
            &render_absorbing(&scene, &medium, cx, &settings).unwrap(),
            xyz_of_spectrum(|lambda| emission.eval(lambda) * scale * response),
            0.005,
        );
    });
}

#[test]
fn g1_dielectric_interior_replaces_instead_of_double_counting_ambient_extinction() {
    with_cx(|cx, _| {
        let density = [1.0];
        let sigma = 0.7;
        let medium = ambient(
            cx,
            &density,
            MassAbsorptionSpectrum::constant(sigma).unwrap(),
        );
        let (mut scene, _) =
            transmission_scene(glass(CauchyIor::try_constant(1.0).unwrap(), sigma), 0.8);
        let mut settings = settings(DirectStrategy::BsdfOnly);
        settings.max_depth = 3;
        settings.spp = 128;
        let with_glass = render_absorbing(&scene, &medium, cx, &settings).unwrap();
        scene.primitives.remove(0);
        scene.lights[0].prim = 0;
        let continuous = render_absorbing(&scene, &medium, cx, &settings).unwrap();
        for channel in 0..3 {
            assert!(
                (with_glass.xyz[0][channel] / continuous.xyz[0][channel] - 1.0).abs() < 2.0e-13,
                "null-index glass with equal extinction must reproduce the displaced field: channel {channel}, glass {}, continuous {}, ratio {}",
                with_glass.xyz[0][channel],
                continuous.xyz[0][channel],
                with_glass.xyz[0][channel] / continuous.xyz[0][channel],
            );
        }
    });
}

fn fog_source_scene(with_glass: bool) -> Scene {
    let mut scene = empty_scene();
    scene.camera.eye.z = 1.0;
    if with_glass {
        scene
            .primitives
            .push(slab(glass(CauchyIor::try_constant(1.5).unwrap(), 0.0)));
    }
    let black = Material::Lambertian {
        reflectance: lift_rgb([0.0; 3]),
    };
    scene.primitives.push(quad(-1.0, 4.0, black));
    let index = scene.primitives.len();
    let emission = (lift_rgb([1.0; 3]), 3.0);
    scene.primitives.push(Primitive {
        shape: Shape::Mesh(TriMesh::new(
            vec![
                [1.0, -0.5, -0.7],
                [1.0, 0.5, -0.7],
                [1.0, 0.5, -0.3],
                [1.0, -0.5, -0.3],
            ],
            vec![[0, 1, 2], [0, 2, 3]],
        )),
        material: black,
        emission: Some(emission),
    });
    scene.lights.push(RectLight {
        corner: Point3::new(1.0, -0.5, -0.7),
        edge_u: Vec3::new(0.0, 1.0, 0.0),
        edge_v: Vec3::new(0.0, 0.0, 0.4),
        prim: index,
        emission,
    });
    scene
}

#[test]
fn g1_scattering_behind_glass_receives_an_off_axis_light() {
    with_cx(|cx, _| {
        let density = [1.0];
        let medium = ScatteringGrid::try_new(
            cx,
            AbsorbingGrid::try_new(
                cx,
                VolumeGrid::new([1; 3], &density, [-2.0, -2.0, -0.9], [4.0, 4.0, 0.6]),
                MassAbsorptionSpectrum::constant(0.0).unwrap(),
                None,
                budget(),
            )
            .unwrap(),
            MassScatteringSpectrum::constant(0.8).unwrap(),
            PhaseFunction::Isotropic,
        )
        .unwrap();
        let mut settings = settings(DirectStrategy::Mis);
        settings.spp = 16_384;
        settings.max_depth = 1;
        let bare = render_scattering(&fog_source_scene(false), &medium, cx, &settings).unwrap();
        settings.max_depth = 3;
        let filtered = render_scattering(&fog_source_scene(true), &medium, cx, &settings).unwrap();
        let (r, _) = interface(1.5, 1.0);
        for channel in 0..3 {
            assert!(bare.xyz[0][channel] > 0.0);
            let ratio = filtered.xyz[0][channel] / bare.xyz[0][channel];
            assert!(
                (ratio / ((1.0 - r) * (1.0 - r)) - 1.0).abs() < 0.025,
                "scattering suffix lost the medium or glass throughput: {ratio}"
            );
        }
    });
}

#[test]
fn g1_rough_glass_light_connections_agree_with_bsdf_transport_through_ambient_absorption() {
    with_cx(|cx, _| {
        let density = [1.0];
        let medium = AbsorbingGrid::try_new(
            cx,
            VolumeGrid::new([1; 3], &density, [-8.0, -8.0, -0.9], [16.0, 16.0, 0.6]),
            MassAbsorptionSpectrum::constant(2.0).unwrap(),
            None,
            budget(),
        )
        .unwrap();
        let material = Material::Dielectric {
            glass: DielectricGlass::new(
                CauchyIor::try_constant(1.5).unwrap(),
                BeerLambertAbsorption::CLEAR,
                GlassProvenance::Custom,
            ),
            surface: DielectricSurface::try_rough(0.25).unwrap(),
        };
        let (scene, _) = transmission_scene(material, 1.0);
        let mut settings = settings(DirectStrategy::BsdfOnly);
        settings.spp = 32_768;
        settings.max_depth = 8;
        let reference = render_absorbing(&scene, &medium, cx, &settings).unwrap();
        // Exterior NEE from the second glass boundary must read the ambient
        // field even though the incident segment was inside the glass stack.
        for strategy in [DirectStrategy::Mis, DirectStrategy::PowerMis] {
            settings.strategy = strategy;
            let film = render_absorbing(&scene, &medium, cx, &settings).unwrap();
            for channel in 0..3 {
                assert!(reference.xyz[0][channel] > 0.0);
                let ratio = film.xyz[0][channel] / reference.xyz[0][channel];
                assert!(
                    (ratio - 1.0).abs() < 0.04,
                    "rough-interface NEE omitted ambient shadow extinction: {ratio}"
                );
            }
        }
    });
}

#[test]
fn g1_lossless_glass_in_conservative_scattering_preserves_uniform_radiance() {
    with_cx(|cx, _| {
        let density = [1.0];
        let medium = ScatteringGrid::try_new(
            cx,
            AbsorbingGrid::try_new(
                cx,
                VolumeGrid::new([1; 3], &density, [-3.0; 3], [6.0; 3]),
                MassAbsorptionSpectrum::constant(0.0).unwrap(),
                None,
                budget(),
            )
            .unwrap(),
            MassScatteringSpectrum::constant(0.15).unwrap(),
            PhaseFunction::Isotropic,
        )
        .unwrap();
        let (mut scene, (emission, scale)) = environment_scene();
        scene.camera.eye.z = 1.0;
        scene
            .primitives
            .push(slab(glass(CauchyIor::try_constant(1.5).unwrap(), 0.0)));
        for strategy in [
            DirectStrategy::BsdfOnly,
            DirectStrategy::NeeOnly,
            DirectStrategy::Mis,
            DirectStrategy::PowerMis,
        ] {
            let mut settings = settings(strategy);
            settings.spp = 8192;
            settings.max_depth = 32;
            assert_mean(
                &render_scattering(&scene, &medium, cx, &settings).unwrap(),
                xyz_of_spectrum(|lambda| emission.eval(lambda) * scale),
                0.035,
            );
        }
    });
}

#[test]
fn g3_refractive_volume_progressive_and_vacuum_paths_replay() {
    with_cx(|cx, _| {
        let density = [1.0];
        let medium = ScatteringGrid::try_new(
            cx,
            ambient(cx, &density, MassAbsorptionSpectrum::constant(0.4).unwrap()),
            MassScatteringSpectrum::constant(0.2).unwrap(),
            PhaseFunction::Rayleigh,
        )
        .unwrap();
        let (scene, _) = transmission_scene(
            glass(CauchyIor::try_new(1.4, 0.012, 0.0).unwrap(), 0.1),
            0.8,
        );
        let mut settings = settings(DirectStrategy::PowerMis);
        settings.max_depth = 6;
        settings.spp = 64;
        let complete = render_scattering(&scene, &medium, cx, &settings).unwrap();
        let mut film = ScatteringFilm::try_new(&medium, &settings).unwrap();
        render_scattering_range(&scene, &medium, cx, &settings, &mut film, 0, 19).unwrap();
        render_scattering_range(&scene, &medium, cx, &settings, &mut film, 19, 64).unwrap();
        assert_eq!(film.film(), &complete);
        let vacuum = ambient(cx, &density, MassAbsorptionSpectrum::constant(0.0).unwrap());
        assert_eq!(
            render_absorbing(&scene, &vacuum, cx, &settings).unwrap(),
            render(&scene, cx, &settings).unwrap()
        );
    });
}
