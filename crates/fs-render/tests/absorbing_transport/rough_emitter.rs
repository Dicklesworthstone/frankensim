//! G1/G3/G5: physical volume sources connect to the ambient side of rough glass.

use super::*;
use fs_exec::RunId;
use fs_render::dielectric::{BeerLambertAbsorption, CauchyIor, GlassProvenance};
use fs_render::tracer::{
    RenderExecutionConfig, ScatteringFilm, render_scattering, render_scattering_range,
    render_scattering_range_with_execution, render_scattering_with_execution,
};
use fs_render::volumes::transport::scattering::{
    MassScatteringSpectrum, PhaseFunction, ScatteringGrid,
};

const TEMPERATURE_K: f64 = 2400.0;
const ABSORPTION: f64 = 0.7;
const ALPHA: f64 = 0.25;

fn glass_scene(a: f64, b: f64, glass_absorption: f64) -> Scene {
    let mut scene = empty_scene();
    scene.camera.eye.z = 0.5;
    // Complete outward-oriented slab. First-interface transmitted directions
    // from the normal camera reach its bottom, away from the side faces.
    scene.primitives.push(Primitive {
        shape: Shape::Mesh(TriMesh::new(
            vec![
                [-8.0, -8.0, -0.2],
                [8.0, -8.0, -0.2],
                [8.0, 8.0, -0.2],
                [-8.0, 8.0, -0.2],
                [-8.0, -8.0, 0.0],
                [8.0, -8.0, 0.0],
                [8.0, 8.0, 0.0],
                [-8.0, 8.0, 0.0],
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
            glass: DielectricGlass::new(
                CauchyIor::try_new(a, b, 0.0).unwrap(),
                BeerLambertAbsorption::try_constant(glass_absorption).unwrap(),
                GlassProvenance::Custom,
            ),
            surface: DielectricSurface::try_rough(ALPHA).unwrap(),
        },
        emission: None,
    });
    scene
}

fn source<'field>(
    cx: &Cx<'_>,
    density: &'field [f64],
    temperature: &'field [f64],
    z_min: f64,
    thickness: f64,
) -> AbsorbingGrid<'field> {
    AbsorbingGrid::try_new(
        cx,
        VolumeGrid::new([1; 3], density, [-1.0, -1.0, z_min], [2.0, 2.0, thickness]),
        MassAbsorptionSpectrum::constant(ABSORPTION).unwrap(),
        Some(temperature),
        budget(),
    )
    .unwrap()
}

// Independent normal-view GGX reflection of [-1,1]^2 × [1,2]. Integrate the
// source's near face: dω = dA/r³, the homogeneous chord emits B*(1-exp(-κL)),
// and f_r*cos_i = F*D*G1(cos_i)/4. This uses no production BSDF, cell traversal,
// proposal, or MIS function. Radiance remains in Planck's SI spectral units.
fn reflected_factor(index: f64, bins: u32) -> f64 {
    let spacing = 2.0 / f64::from(bins);
    let alpha2 = ALPHA * ALPHA;
    let mut integral = 0.0;
    for ix in 0..bins {
        for iy in 0..bins {
            let x = -1.0 + (f64::from(ix) + 0.5) * spacing;
            let y = -1.0 + (f64::from(iy) + 0.5) * spacing;
            let distance2 = x * x + y * y + 1.0;
            let distance = distance2.sqrt();
            let cosine = 1.0 / distance;
            let half_cosine2 = 0.5 * (1.0 + cosine);
            let half_cosine = half_cosine2.sqrt();
            let transmitted_cosine = (1.0 - (1.0 - half_cosine2) / (index * index)).sqrt();
            let rs = (half_cosine - index * transmitted_cosine)
                / (half_cosine + index * transmitted_cosine);
            let rp = (index * half_cosine - transmitted_cosine)
                / (index * half_cosine + transmitted_cosine);
            let fresnel = 0.5 * (rs * rs + rp * rp);
            let denominator = 1.0 + (alpha2 - 1.0) * half_cosine2;
            let distribution = alpha2 / (core::f64::consts::PI * denominator * denominator);
            let masking =
                2.0 * cosine / (cosine + (alpha2 + (1.0 - alpha2) * cosine * cosine).sqrt());
            let exit_distance = (2.0 * distance)
                .min(distance / x.abs())
                .min(distance / y.abs());
            let emitted_fraction = -(-ABSORPTION * (exit_distance - distance)).exp_m1();
            integral += fresnel * distribution * masking * emitted_fraction * spacing * spacing
                / (4.0 * distance2 * distance);
        }
    }
    integral
}

fn reflected_reference(a: f64, b: f64) -> [f64; 3] {
    xyz_of_spectrum(|lambda| {
        let micrometres = lambda * 1.0e-3;
        planck(TEMPERATURE_K, lambda) * reflected_factor(a + b / (micrometres * micrometres), 96)
    })
}

#[test]
fn g1_terminal_rough_reflection_matches_independent_source_integral() {
    with_cx(|cx, _| {
        let density = [1.0];
        let temperature = [TEMPERATURE_K];
        let medium = source(cx, &density, &temperature, 1.0, 1.0);
        let scene = glass_scene(1.7, 0.0, 0.0);
        assert!((reflected_factor(1.7, 48) / reflected_factor(1.7, 96) - 1.0).abs() < 0.001);
        let mut settings = settings(DirectStrategy::NeeOnly);
        settings.spp = 16_384;
        settings.max_depth = 1;
        let direct = render_absorbing(&scene, &medium, cx, &settings).unwrap();
        assert_mean(&direct, reflected_reference(1.7, 0.0), 0.03);
        for strategy in [DirectStrategy::Mis, DirectStrategy::PowerMis] {
            settings.strategy = strategy;
            assert_eq!(
                render_absorbing(&scene, &medium, cx, &settings).unwrap(),
                direct,
                "the last rough vertex has no traced source continuation"
            );
        }
        settings.strategy = DirectStrategy::BsdfOnly;
        assert_eq!(
            render_absorbing(&scene, &medium, cx, &settings)
                .unwrap()
                .xyz,
            vec![[0.0; 3]],
            "the source is behind the camera and no continuation is traced"
        );
    });
}

#[test]
fn g1_dispersive_rough_source_mis_preserves_each_wavelengths_energy() {
    with_cx(|cx, _| {
        let density = [1.0];
        let temperature = [TEMPERATURE_K];
        let medium = source(cx, &density, &temperature, 1.0, 1.0);
        let scene = glass_scene(1.7, 0.04, 0.0);
        let expected = reflected_reference(1.7, 0.04);
        for strategy in [
            DirectStrategy::BsdfOnly,
            DirectStrategy::NeeOnly,
            DirectStrategy::Mis,
            DirectStrategy::PowerMis,
        ] {
            let mut settings = settings(strategy);
            settings.spp = 32_768;
            settings.max_depth = 2;
            assert_mean(
                &render_absorbing(&scene, &medium, cx, &settings).unwrap(),
                expected,
                0.055,
            );
        }
    });
}

#[test]
fn g1_rough_exit_source_connections_match_transmitted_path_radiance() {
    with_cx(|cx, _| {
        let density = [1.0];
        let temperature = [TEMPERATURE_K];
        let medium = source(cx, &density, &temperature, -2.0, 1.0);
        let scene = glass_scene(1.7, 0.04, 0.8);
        // Only transmission through both slab interfaces reaches the source.
        // Depth three admits its escaped source segment; a terminal direct
        // connection reaches the same integral already at the second vertex.
        let mut settings = settings(DirectStrategy::BsdfOnly);
        settings.spp = 32_768;
        settings.max_depth = 3;
        let reference = render_absorbing(&scene, &medium, cx, &settings).unwrap();
        assert!(reference.xyz[0].iter().all(|value| *value > 0.0));
        let expected = reference.xyz[0].map(|value| value / f64::from(settings.spp));
        for strategy in [
            DirectStrategy::NeeOnly,
            DirectStrategy::Mis,
            DirectStrategy::PowerMis,
        ] {
            settings.strategy = strategy;
            assert_mean(
                &render_absorbing(&scene, &medium, cx, &settings).unwrap(),
                expected,
                0.045,
            );
        }
        settings.max_depth = 2;
        settings.strategy = DirectStrategy::NeeOnly;
        let terminal = render_absorbing(&scene, &medium, cx, &settings).unwrap();
        assert_mean(&terminal, expected, 0.045);
        for strategy in [DirectStrategy::Mis, DirectStrategy::PowerMis] {
            settings.strategy = strategy;
            assert_eq!(
                render_absorbing(&scene, &medium, cx, &settings).unwrap(),
                terminal
            );
        }
    });
}

#[test]
fn g3_hot_grid_displaced_by_glass_has_no_interior_source_connections() {
    with_cx(|cx, _| {
        let density = [1.0];
        let temperature = [TEMPERATURE_K];
        // The entire proposed emitting cell is inside the closed glass body.
        // Neither a connection into glass nor an internal reflected segment
        // may treat that ambient field as a source in the dielectric interior.
        let medium = source(cx, &density, &temperature, -0.18, 0.16);
        let scene = glass_scene(1.7, 0.04, 0.0);
        for strategy in [
            DirectStrategy::BsdfOnly,
            DirectStrategy::NeeOnly,
            DirectStrategy::Mis,
            DirectStrategy::PowerMis,
        ] {
            let mut settings = settings(strategy);
            settings.spp = 128;
            settings.max_depth = 5;
            assert_eq!(
                render_absorbing(&scene, &medium, cx, &settings)
                    .unwrap()
                    .xyz,
                vec![[0.0; 3]],
                "ambient source emission leaked into the displaced dielectric volume"
            );
        }
    });
}

#[test]
fn g5_rough_source_paths_replay_across_scattering_tiles_and_sample_ranges() {
    with_cx(|cx, _| {
        let density = [1.0];
        let temperature = [TEMPERATURE_K];
        let medium = ScatteringGrid::try_new(
            cx,
            source(cx, &density, &temperature, 1.0, 1.0),
            MassScatteringSpectrum::constant(0.4).unwrap(),
            PhaseFunction::Rayleigh,
        )
        .unwrap();
        let mut scene = glass_scene(1.7, 0.04, 0.2);
        scene.camera.half_tan = 0.2;
        let execution = RenderExecutionConfig::try_new(2, 1, 3, 8 << 20, RunId(972)).unwrap();
        for sampler in [Sampler::Iid, Sampler::OwenSobolFullPath] {
            let mut settings = settings(DirectStrategy::PowerMis);
            settings.width = 3;
            settings.height = 2;
            settings.spp = 32;
            settings.max_depth = 5;
            settings.sampler = sampler;
            let serial = render_scattering(&scene, &medium, cx, &settings).unwrap();
            assert!(serial.xyz.iter().any(|pixel| pixel[1] > 0.0));
            assert_eq!(
                render_scattering_with_execution(&scene, &medium, cx, &settings, &execution)
                    .unwrap()
                    .film,
                serial
            );
            let mut progressive = ScatteringFilm::try_new(&medium, &settings).unwrap();
            render_scattering_range(&scene, &medium, cx, &settings, &mut progressive, 0, 11)
                .unwrap();
            render_scattering_range_with_execution(
                &scene,
                &medium,
                cx,
                &settings,
                &mut progressive,
                11,
                settings.spp,
                &execution,
            )
            .unwrap();
            assert_eq!(*progressive.film(), serial);
        }
    });
}
