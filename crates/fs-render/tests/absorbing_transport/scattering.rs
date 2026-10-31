//! Independent transport integrals and conservation checks for coupled paths.

use super::*;
use fs_render::tracer::{ScatteringFilm, render_scattering, render_scattering_range};
use fs_render::volumes::transport::scattering::{
    MassScatteringSpectrum, PhaseFunction, ScatteringGrid,
};

#[path = "scattering/execution.rs"]
mod execution;

fn physical_grid<'a>(
    cx: &Cx<'_>,
    density: &'a [f64],
    origin: [f64; 3],
    cell: [f64; 3],
    absorption: f64,
    scattering: f64,
    phase: PhaseFunction,
) -> ScatteringGrid<'a> {
    ScatteringGrid::try_new(
        cx,
        AbsorbingGrid::try_new(
            cx,
            VolumeGrid::new([1; 3], density, origin, cell),
            MassAbsorptionSpectrum::constant(absorption).unwrap(),
            None,
            budget(),
        )
        .unwrap(),
        MassScatteringSpectrum::constant(scattering).unwrap(),
        phase,
    )
    .unwrap()
}

fn backlit_slab() -> (Scene, (LiftedSpectrum, f64)) {
    let mut scene = empty_scene();
    scene.camera.eye = Point3::new(0.0, 0.0, 0.0);
    scene.camera.forward = Vec3::new(0.0, 0.0, 1.0);
    let black = Material::Lambertian {
        reflectance: lift_rgb([0.0; 3]),
    };
    scene.primitives.push(quad(1.0, 4.0, black));
    let emission = (lift_rgb([1.0; 3]), 3.0);
    let mut emitter = quad(-1.0, 0.5, black);
    emitter.emission = Some(emission);
    scene.primitives.push(emitter);
    scene.lights.push(RectLight {
        corner: Point3::new(-0.5, -0.5, -1.0),
        edge_u: Vec3::new(1.0, 0.0, 0.0),
        edge_v: Vec3::new(0.0, 1.0, 0.0),
        prim: 1,
        emission,
    });
    (scene, emission)
}

// Independent midpoint quadrature over position and source area. Camera
// traverses z=0..1, source is a unit square at z=-1. Both beam extinctions,
// physical phase value and projected emitter area appear explicitly.
fn single_scatter_integral(absorption: f64, scattering: f64, g: f64) -> f64 {
    let sigma = absorption + scattering;
    let n = 96;
    let step = 1.0 / f64::from(n);
    let mut result = 0.0;
    for iz in 0..n {
        let z = (f64::from(iz) + 0.5) * step;
        for ix in 0..n {
            let x = (f64::from(ix) + 0.5) * step - 0.5;
            for iy in 0..n {
                let y = (f64::from(iy) + 0.5) * step - 0.5;
                let distance = (x * x + y * y + (1.0 + z) * (1.0 + z)).sqrt();
                let cosine = -(1.0 + z) / distance;
                let phase = (1.0 - g * g)
                    / (4.0 * core::f64::consts::PI * (1.0 + g * g - 2.0 * g * cosine).powf(1.5));
                let shadow_length = z * distance / (1.0 + z);
                result += (-sigma * (z + shadow_length)).exp() * scattering * phase * (1.0 + z)
                    / (distance * distance * distance)
                    * step
                    * step
                    * step;
            }
        }
    }
    result
}

#[test]
fn g1_single_scattering_matches_independent_source_integral_and_terminal_mis() {
    with_cx(|cx, _| {
        let density = [1.0];
        let (scene, (emission, scale)) = backlit_slab();
        for g in [0.0, 0.6, -0.6] {
            let medium = physical_grid(
                cx,
                &density,
                [-2.0, -2.0, 0.0],
                [4.0, 4.0, 1.0],
                0.3,
                0.7,
                PhaseFunction::henyey_greenstein(g).unwrap(),
            );
            let mut settings = settings(DirectStrategy::NeeOnly);
            settings.max_depth = 1;
            settings.spp = 16_384;
            let film = render_scattering(&scene, &medium, cx, &settings).unwrap();
            let integral = single_scatter_integral(0.3, 0.7, g);
            assert_mean(
                &film,
                xyz_of_spectrum(|lambda| emission.eval(lambda) * scale * integral),
                0.018,
            );
            for strategy in [DirectStrategy::Mis, DirectStrategy::PowerMis] {
                settings.strategy = strategy;
                assert_eq!(
                    render_scattering(&scene, &medium, cx, &settings).unwrap(),
                    film,
                    "a terminal volume vertex has no competing continuation"
                );
            }
        }
    });
}

#[test]
fn g1_physical_thermal_mesh_light_scattering_preserves_planck_radiance() {
    with_cx(|cx, _| {
        let density = [1.0];
        let medium = physical_grid(
            cx,
            &density,
            [-2.0, -2.0, 0.0],
            [4.0, 4.0, 1.0],
            0.3,
            0.7,
            PhaseFunction::Isotropic,
        );
        let (mut scene, _) = backlit_slab();
        scene.lights.clear();
        scene.primitives[1].emission = None;
        scene.primitives[1].material = Material::GrayDiffuseThermal {
            surface: GrayDiffuseThermalSurface::try_new(2000.0, 0.7, ContentHash([71; 32]))
                .unwrap(),
        };
        let mut settings = settings(DirectStrategy::Mis);
        settings.max_depth = 1;
        settings.spp = 16_384;
        let film = render_scattering(&scene, &medium, cx, &settings).unwrap();
        let integral = single_scatter_integral(0.3, 0.7, 0.0);
        assert_mean(
            &film,
            xyz_of_spectrum(|lambda| 0.7 * planck(2000.0, lambda) * integral),
            0.025,
        );
    });
}

#[test]
fn g1_surface_direct_lighting_uses_total_extinction_on_shadow_segments() {
    with_cx(|cx, _| {
        let density = [2.0];
        let medium = physical_grid(
            cx,
            &density,
            [-2.0, -2.0, 1.0],
            [4.0, 4.0, 0.5],
            0.3,
            0.4,
            PhaseFunction::Rayleigh,
        );
        let (scene, reflectance, (emission, scale)) = lit_receiver();
        let mut settings = settings(DirectStrategy::Mis);
        settings.max_depth = 1;
        let film = render_scattering(&scene, &medium, cx, &settings).unwrap();
        let integral = shadow_square_integral(1.4);
        assert_mean(
            &film,
            xyz_of_spectrum(|lambda| {
                reflectance.eval(lambda) * emission.eval(lambda) * scale * integral
                    / core::f64::consts::PI
            }),
            0.006,
        );
    });
}

#[test]
fn g1_multiple_spectral_scattering_conserves_uniform_environment_radiance() {
    with_cx(|cx, _| {
        let density = [1.0];
        let medium = ScatteringGrid::try_new(
            cx,
            AbsorbingGrid::try_new(
                cx,
                VolumeGrid::new([1; 3], &density, [-1.0; 3], [2.0; 3]),
                MassAbsorptionSpectrum::constant(0.0).unwrap(),
                None,
                budget(),
            )
            .unwrap(),
            MassScatteringSpectrum::table(&[(380.0, 0.05), (580.0, 1.2), (780.0, 0.1)]).unwrap(),
            PhaseFunction::HenyeyGreenstein { g: 0.35 },
        )
        .unwrap();
        let (mut scene, (emission, scale)) = environment_scene();
        scene.camera.eye = Point3::new(0.0, 0.0, 0.0);
        let expected = xyz_of_spectrum(|lambda| emission.eval(lambda) * scale);
        for strategy in [
            DirectStrategy::BsdfOnly,
            DirectStrategy::NeeOnly,
            DirectStrategy::Mis,
            DirectStrategy::PowerMis,
        ] {
            let mut settings = settings(strategy);
            settings.max_depth = 64;
            settings.spp = 32_768;
            assert_mean(
                &render_scattering(&scene, &medium, cx, &settings).unwrap(),
                expected,
                0.035,
            );
        }
    });
}

#[test]
fn g1_spectral_phase_direct_light_matches_wavelength_resolved_source_integrals() {
    with_cx(|cx, _| {
        let density = [1.0];
        let (scene, (emission, scale)) = backlit_slab();
        let medium = physical_grid(
            cx,
            &density,
            [-2.0, -2.0, 0.0],
            [4.0, 4.0, 1.0],
            0.3,
            0.7,
            PhaseFunction::spectral_table(&[
                (380.0, PhaseFunction::henyey_greenstein(0.6).unwrap()),
                (780.0, PhaseFunction::henyey_greenstein(-0.6).unwrap()),
            ])
            .unwrap(),
        );
        // Linearity of the transport integral in the angular law supplies
        // an independent reference at EVERY wavelength, not a gray average.
        let blue = single_scatter_integral(0.3, 0.7, 0.6);
        let red = single_scatter_integral(0.3, 0.7, -0.6);
        let expected = xyz_of_spectrum(|lambda| {
            let t = (lambda - 380.0) / 400.0;
            emission.eval(lambda) * scale * ((1.0 - t) * blue + t * red)
        });
        let mut settings = settings(DirectStrategy::NeeOnly);
        settings.max_depth = 1;
        settings.spp = 16_384;
        let film = render_scattering(&scene, &medium, cx, &settings).unwrap();
        assert_mean(&film, expected, 0.018);
        for strategy in [DirectStrategy::Mis, DirectStrategy::PowerMis] {
            settings.strategy = strategy;
            assert_eq!(
                render_scattering(&scene, &medium, cx, &settings).unwrap(),
                film
            );
        }
    });
}

#[test]
fn g1_spectral_phase_multiple_scattering_preserves_the_environment_furnace() {
    with_cx(|cx, _| {
        let density = [1.0];
        let medium = physical_grid(
            cx,
            &density,
            [-1.0; 3],
            [2.0; 3],
            0.0,
            0.4,
            PhaseFunction::spectral_table(&[
                (380.0, PhaseFunction::henyey_greenstein(0.65).unwrap()),
                (780.0, PhaseFunction::henyey_greenstein(-0.6).unwrap()),
            ])
            .unwrap(),
        );
        let (mut scene, (emission, scale)) = environment_scene();
        scene.camera.eye = Point3::new(0.0, 0.0, 0.0);
        let expected = xyz_of_spectrum(|lambda| emission.eval(lambda) * scale);
        for strategy in [
            DirectStrategy::BsdfOnly,
            DirectStrategy::NeeOnly,
            DirectStrategy::Mis,
            DirectStrategy::PowerMis,
        ] {
            let mut settings = settings(strategy);
            settings.max_depth = 32;
            settings.spp = 32_768;
            assert_mean(
                &render_scattering(&scene, &medium, cx, &settings).unwrap(),
                expected,
                0.035,
            );
        }
    });
}

#[test]
fn g3_zero_scattering_reuses_absorption_and_vacuum_exactly() {
    with_cx(|cx, _| {
        let density = [1.0];
        let (scene, _) = environment_scene();
        for absorption in [0.0, 0.7] {
            let medium = physical_grid(
                cx,
                &density,
                [-1.0; 3],
                [2.0; 3],
                absorption,
                0.0,
                PhaseFunction::HenyeyGreenstein { g: 0.9 },
            );
            for strategy in [
                DirectStrategy::BsdfOnly,
                DirectStrategy::NeeOnly,
                DirectStrategy::Mis,
                DirectStrategy::PowerMis,
            ] {
                let mut settings = settings(strategy);
                settings.spp = 128;
                let expected =
                    render_absorbing(&scene, medium.absorption_grid(), cx, &settings).unwrap();
                assert_eq!(
                    render_scattering(&scene, &medium, cx, &settings).unwrap(),
                    expected
                );
                if absorption == 0.0 {
                    assert_eq!(expected, render(&scene, cx, &settings).unwrap());
                }
            }
        }
    });
}

#[test]
fn g3_progressive_scattering_replays_and_rejects_changed_phase_atomically() {
    with_cx(|cx, gate| {
        let density = [1.0];
        let make = |phase| physical_grid(cx, &density, [-1.0; 3], [2.0; 3], 0.2, 0.7, phase);
        let medium = make(PhaseFunction::Isotropic);
        let (scene, _) = environment_scene();
        let mut settings = settings(DirectStrategy::PowerMis);
        settings.spp = 128;
        settings.max_depth = 8;
        let full = render_scattering(&scene, &medium, cx, &settings).unwrap();
        let mut film = ScatteringFilm::try_new(&medium, &settings).unwrap();
        render_scattering_range(&scene, &medium, cx, &settings, &mut film, 0, 37).unwrap();
        let checkpoint = film.clone();
        render_scattering_range(&scene, &medium, cx, &settings, &mut film, 37, 128).unwrap();
        assert_eq!(film.film(), &full);
        let mut film = checkpoint.clone();
        let changed = make(PhaseFunction::Rayleigh);
        assert_eq!(
            render_scattering_range(&scene, &changed, cx, &settings, &mut film, 37, 128),
            Err(TracerError::InvalidInput)
        );
        assert_eq!(film, checkpoint);
        gate.request();
        assert_eq!(
            render_scattering_range(&scene, &medium, cx, &settings, &mut film, 37, 128),
            Err(TracerError::Cancelled)
        );
        assert_eq!(film, checkpoint);
    });
}
