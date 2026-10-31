//! Conditional medium/surface depth tails do not certify sampled results.

use super::*;
use fs_render::photometry::IlluminanceProbe;
use fs_render::photometry::area::{EmissionSidedness, PhotometricAreaEmitter, PhotometricEmitters};
use fs_render::photometry::path::{
    DiffuseTransportRemainder, IlluminanceSampling, PathIlluminanceJob,
};
use fs_render::spectral::SpectralReflectance;

fn cavity() -> (Scene, [PhotometricAreaEmitter; 1]) {
    let mut scene = empty_scene([0.13, -0.17, 0.09], [0.0, 0.0, 1.0]);
    let mut wall = super::area::diffuse_emitting_cube(0.0);
    wall.material = Material::SpectralDiffuse {
        reflectance: SpectralReflectance::constant(0.25).unwrap(),
    };
    scene.primitives.push(wall);
    let emitter = PhotometricAreaEmitter::try_new(
        0,
        24.0 * PI * 8.0,
        1.0,
        SourceSpectrum::equal_energy(),
        EmissionSidedness::Front,
    )
    .unwrap();
    (scene, [emitter])
}

fn grid<'a>(
    cx: &Cx<'_>,
    density: &'a [f64],
    absorption: MassAbsorptionSpectrum,
    temperature: Option<&'a [f64]>,
) -> AbsorbingGrid<'a> {
    AbsorbingGrid::try_new(
        cx,
        VolumeGrid::new([1, 1, 2], density, [-1.0; 3], [2.0, 2.0, 1.0]),
        absorption,
        temperature,
        medium_budget(),
    )
    .unwrap()
}

fn check_both(
    cx: &Cx<'_>,
    scene: &Scene,
    surfaces: &[PhotometricAreaEmitter],
    medium: PathMedium<'_, '_>,
) -> Option<DiffuseTransportRemainder> {
    let emitters = PhotometricEmitters {
        points: &[],
        surfaces,
    };
    let mut image = PathLuminanceJob::new_with_emitters(
        scene,
        emitters,
        medium,
        TransportOptions::default(),
        settings(3, 2),
        cx,
    )
    .unwrap();
    let probes = [
        IlluminanceProbe::try_new(Point3::new(0.0, 0.0, -0.4), Vec3::new(0.6, 0.0, 0.8)).unwrap(),
    ];
    let mut probe = PathIlluminanceJob::new_with_emitters(
        scene,
        emitters,
        medium,
        &probes,
        TransportOptions::default(),
        IlluminanceSampling {
            samples_per_probe: 2,
            max_depth: 3,
            sampler: Sampler::OwenSobol,
            seed: 173,
        },
        cx,
    )
    .unwrap();
    let remainder = image.diffuse_transport_remainder();
    assert_eq!(probe.diffuse_transport_remainder(), remainder);
    if !matches!(medium, PathMedium::Vacuum) {
        // The earlier explicit-vacuum API keeps its narrower admission contract.
        assert!(image.diffuse_vacuum_remainder().is_none());
        assert!(probe.diffuse_vacuum_remainder().is_none());
    }
    for budget in [1, usize::MAX] {
        image.advance(cx, budget).unwrap();
        probe.advance(cx, budget).unwrap();
        assert_eq!(image.diffuse_transport_remainder(), remainder);
        assert_eq!(probe.diffuse_transport_remainder(), remainder);
        let pixel = image.pixel_estimate(0).unwrap();
        let value = probe.probe_estimate(0).unwrap();
        assert_eq!(pixel.numerical.kind, fs_evidence::NumericalKind::NoClaim);
        assert_eq!(value.numerical.kind, fs_evidence::NumericalKind::NoClaim);
        assert!(pixel.total_error_bound_cd_m2.is_none());
        assert!(value.total_error_bound_lux.is_none());
    }
    remainder
}

#[test]
fn g1_nonemitting_absorption_preserves_the_diffuse_vacuum_tail() {
    with_cx(|cx, _| {
        let (scene, surfaces) = cavity();
        let vacuum = check_both(cx, &scene, &surfaces, PathMedium::Vacuum).unwrap();
        assert_eq!(vacuum.max_depth(), 3);
        assert_eq!(vacuum.max_surface_reflectance(), 0.25);
        assert_eq!(vacuum.max_scattering_gain(), 0.0);
        assert_eq!(vacuum.contraction_upper_bound(), 0.25);
        close(
            vacuum.source_luminance_ceiling_cd_m2(),
            8.0,
            3e-12,
            "area luminance",
        );
        // Absorption can remove paths but cannot add a continuation order.
        let tail = 8.0 * 0.25_f64.powi(3) / (1.0 - 0.25);
        close(
            vacuum.luminance_bound_cd_m2(),
            tail,
            3e-12,
            "absorbing tail",
        );
        for density in [[0.0, 0.0], [0.25, 2.0], [200.0, 0.0]] {
            let medium = grid(
                cx,
                &density,
                MassAbsorptionSpectrum::table(&[(380.0, 0.1), (550.0, 2.0), (780.0, 0.3)]).unwrap(),
                None,
            );
            assert_eq!(
                check_both(cx, &scene, &surfaces, PathMedium::Absorbing(&medium)),
                Some(vacuum)
            );
        }
    });
}

#[test]
fn g1_colored_scattering_tail_uses_combined_knots_and_all_density_layers() {
    with_cx(|cx, _| {
        let (scene, surfaces) = cavity();
        let phase =
            PhaseFunction::tabulated(&[(-1.0, 0.2), (-0.3, 0.6), (0.5, 1.9), (1.0, 0.3)]).unwrap();
        // The albedo maximum is at an absorption-only knot, not at either
        // band endpoint or any scattering knot: ks(500)=14/11, ka(500)=1/4.
        let omega = 56.0 / 67.0;
        let nominal_tail = 8.0 * omega * omega * omega / (1.0 - omega);
        let mut previous = None;
        for density in [[0.5, 0.5], [0.25, 2.0], [30.0, 0.01]] {
            let medium = ScatteringGrid::try_new(
                cx,
                grid(
                    cx,
                    &density,
                    MassAbsorptionSpectrum::table(&[
                        (380.0, 3.0),
                        (500.0, 0.25),
                        (620.0, 3.0),
                        (780.0, 3.0),
                    ])
                    .unwrap(),
                    None,
                ),
                MassScatteringSpectrum::table(&[(380.0, 1.0), (600.0, 1.5), (780.0, 1.0)]).unwrap(),
                phase.clone(),
            )
            .unwrap();
            let tail = check_both(cx, &scene, &surfaces, PathMedium::Scattering(&medium)).unwrap();
            assert_eq!(tail.max_depth(), 3);
            assert_eq!(tail.max_surface_reflectance(), 0.25);
            assert!(tail.max_scattering_gain() >= omega);
            close(
                tail.max_scattering_gain(),
                omega,
                3e-12,
                "combined-knot albedo",
            );
            assert_eq!(tail.contraction_upper_bound(), tail.max_scattering_gain());
            close(
                tail.luminance_bound_cd_m2(),
                nominal_tail,
                3e-12,
                "scattering tail",
            );
            close(
                tail.illuminance_bound_lux(),
                PI * nominal_tail,
                3e-12,
                "hemispheric tail",
            );
            assert!(tail.luminance_bound_cd_m2() + 1e-12 >= nominal_tail);
            if let Some(previous) = previous {
                assert_eq!(tail, previous);
            }
            previous = Some(tail);
        }
        // Check the retained phase interpolant independently by trapezoids;
        // the bound must not interpret normalization roundoff as absorption.
        let PhaseFunction::Tabulated(table) = phase else {
            unreachable!()
        };
        let mass: f64 = table
            .knots_per_sr()
            .windows(2)
            .map(|pair| PI * (pair[1].0 - pair[0].0) * (pair[0].1 + pair[1].1))
            .sum();
        close(mass, 1.0, 3e-14, "normalized retained phase mass");
        assert!(previous.unwrap().max_scattering_gain() + 2e-15 >= omega * mass);
    });
}

#[test]
fn g0_conservative_scattering_and_thermal_sources_have_no_medium_tail_claim() {
    with_cx(|cx, _| {
        let (scene, surfaces) = cavity();
        let density = [0.25, 2.0];
        let temperatures = [1800.0, 2200.0];
        for (absorption, temperature) in [(0.0, None), (0.25, Some(temperatures.as_slice()))] {
            let absorbing = grid(
                cx,
                &density,
                MassAbsorptionSpectrum::constant(absorption).unwrap(),
                temperature,
            );
            if temperature.is_some() {
                assert!(
                    check_both(cx, &scene, &surfaces, PathMedium::Absorbing(&absorbing)).is_none()
                );
            }
            let medium = ScatteringGrid::try_new(
                cx,
                absorbing,
                MassScatteringSpectrum::constant(0.75).unwrap(),
                PhaseFunction::Isotropic,
            )
            .unwrap();
            assert!(check_both(cx, &scene, &surfaces, PathMedium::Scattering(&medium)).is_none());
        }
    });
}
