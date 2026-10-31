//! Conditional model tails are separate from sampling and numerical accuracy.

use super::*;
use fs_render::photometry::IlluminanceProbe;
use fs_render::photometry::area::{
    EmissionSidedness, GoboTransmission, PhotometricAreaEmitter, PhotometricEmitters,
    ThinLensProjector,
};
use fs_render::photometry::path::{
    DiffuseVacuumRemainder, IlluminanceSampling, PathIlluminanceJob,
};
use fs_render::spectral::SpectralReflectance;

fn diffuse(rho: f64) -> Material {
    Material::SpectralDiffuse {
        reflectance: SpectralReflectance::constant(rho).unwrap(),
    }
}

fn cavity(rho: f64) -> Scene {
    let mut scene = empty_scene([0.0; 3], [0.0, 0.0, 1.0]);
    let mut wall = super::area::diffuse_emitting_cube(0.0);
    wall.material = diffuse(rho);
    scene.primitives.push(wall);
    scene
}

fn cavity_emitter() -> PhotometricAreaEmitter {
    // Original CIE rows flatten photopic power independently of production
    // observer lookup. The inward cube has exactly 24 m² and emits Le=8 cd/m².
    let spectrum = SourceSpectrum::try_from_bins(core::array::from_fn(|bin| {
        1.0 / super::spectral_reference::cie_reference_xyz(382.5 + 5.0 * bin as f64)[1]
    }))
    .unwrap();
    PhotometricAreaEmitter::try_new(0, 24.0 * PI * 8.0, 1.0, spectrum, EmissionSidedness::Front)
        .unwrap()
}

fn bound(
    scene: &Scene,
    emitters: PhotometricEmitters<'_>,
    medium: PathMedium<'_, '_>,
    cx: &Cx<'_>,
) -> Option<DiffuseVacuumRemainder> {
    PathLuminanceJob::new_with_emitters(
        scene,
        emitters,
        medium,
        TransportOptions::default(),
        settings(3, 1),
        cx,
    )
    .unwrap()
    .diffuse_vacuum_remainder()
}

#[test]
fn g1_diffuse_cavity_tail_matches_geometric_series_without_certifying_samples() {
    with_cx(|cx, _| {
        let scene = cavity(0.5);
        let surfaces = [cavity_emitter()];
        let emitters = PhotometricEmitters {
            points: &[],
            surfaces: &surfaces,
        };
        let probes =
            [
                IlluminanceProbe::try_new(Point3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 0.0, 1.0))
                    .unwrap(),
            ];
        for depth in [1, 4, 12] {
            let mut image = PathLuminanceJob::new_with_emitters(
                &scene,
                emitters,
                PathMedium::Vacuum,
                TransportOptions::default(),
                Settings {
                    strategy: DirectStrategy::BsdfOnly,
                    ..settings(depth, 8192)
                },
                cx,
            )
            .unwrap();
            let remainder = image.diffuse_vacuum_remainder().unwrap();
            let tail = 8.0 * 0.5_f64.powi(depth as i32) / (1.0 - 0.5);
            assert_eq!(remainder.max_depth(), depth);
            assert_eq!(remainder.max_reflectance(), 0.5);
            close(
                remainder.source_luminance_ceiling_cd_m2(),
                8.0,
                3.0e-12,
                "source ceiling",
            );
            close(
                remainder.luminance_bound_cd_m2(),
                tail,
                3.0e-12,
                "geometric tail",
            );
            // Allow only normalization roundoff relative to nominal Le=8.
            assert!(remainder.luminance_bound_cd_m2() + 1.0e-12 >= tail);
            assert!(remainder.illuminance_bound_lux() + 1.0e-12 >= PI * tail);
            image.advance(cx, 1).unwrap();
            assert_eq!(image.diffuse_vacuum_remainder(), Some(remainder));
            image.advance(cx, usize::MAX).unwrap();
            let estimate = image.pixel_estimate(0).unwrap();
            // A separate sampling-regression tolerance, never the model tail.
            close(
                estimate.luminance_cd_m2,
                16.0 - tail,
                0.005,
                "finite-depth cavity mean",
            );
            assert_eq!(estimate.numerical.kind, fs_evidence::NumericalKind::NoClaim);
            assert!(estimate.total_error_bound_cd_m2.is_none());
            assert_eq!(image.diffuse_vacuum_remainder(), Some(remainder));
            let mut probe = PathIlluminanceJob::new_with_emitters(
                &scene,
                emitters,
                PathMedium::Vacuum,
                &probes,
                TransportOptions::default(),
                IlluminanceSampling {
                    samples_per_probe: 2,
                    max_depth: depth,
                    sampler: Sampler::Iid,
                    seed: 173,
                },
                cx,
            )
            .unwrap();
            assert_eq!(probe.diffuse_vacuum_remainder(), Some(remainder));
            probe.advance(cx, 1).unwrap();
            assert_eq!(probe.diffuse_vacuum_remainder(), Some(remainder));
            let estimate = probe.probe_estimate(0).unwrap();
            assert_eq!(estimate.numerical.kind, fs_evidence::NumericalKind::NoClaim);
            assert!(estimate.total_error_bound_lux.is_none());
        }
    });
}

#[test]
fn g0_global_reflectance_includes_unlit_colored_geometry_and_exact_black() {
    with_cx(|cx, _| {
        let surfaces = [cavity_emitter()];
        let emitters = PhotometricEmitters {
            points: &[],
            surfaces: &surfaces,
        };
        let mut scene = cavity(0.25);
        let mut hidden = floor(3.0, 1.0, LiftedSpectrum { c: [0.0; 3] });
        hidden.material = Material::SpectralDiffuse {
            reflectance: SpectralReflectance::table(&[(380.0, 0.1), (500.0, 0.75), (780.0, 0.3)])
                .unwrap(),
        };
        scene.primitives.push(hidden);
        let remainder = bound(&scene, emitters, PathMedium::Vacuum, cx).unwrap();
        assert_eq!(remainder.max_reflectance(), 0.75);
        close(
            remainder.luminance_bound_cd_m2(),
            13.5,
            3.0e-12,
            "global spectral contraction",
        );
        let black = cavity(0.0);
        let remainder = bound(&black, emitters, PathMedium::Vacuum, cx).unwrap();
        assert_eq!(remainder.luminance_bound_cd_m2(), 0.0);
        assert_eq!(remainder.illuminance_bound_lux(), 0.0);
    });
}

#[test]
fn g0_model_remainder_refuses_noncontracting_or_unsupported_transport() {
    with_cx(|cx, _| {
        let surfaces = [cavity_emitter()];
        let emitters = PhotometricEmitters {
            points: &[],
            surfaces: &surfaces,
        };
        for material in [
            diffuse(1.0),
            Material::Lambertian {
                reflectance: LiftedSpectrum { c: [0.0; 3] },
            },
            Material::Ggx {
                reflectance: LiftedSpectrum { c: [0.0; 3] },
                alpha: 0.3,
            },
        ] {
            let mut scene = cavity(0.5);
            scene.primitives[0].material = material;
            assert!(bound(&scene, emitters, PathMedium::Vacuum, cx).is_none());
        }
        let scene = cavity(0.5);
        for intensity in [0.0, 10.0] {
            let points = [isotropic(
                [0.5; 3],
                intensity,
                SourceSpectrum::equal_energy(),
            )];
            let with_point = PhotometricEmitters {
                points: &points,
                surfaces: &surfaces,
            };
            assert_eq!(
                bound(&scene, with_point, PathMedium::Vacuum, cx).is_some(),
                intensity == 0.0
            );
        }
        // Even a zero-density physical grid is outside the explicit-vacuum model.
        let density = [0.0];
        let medium = AbsorbingGrid::try_new(
            cx,
            VolumeGrid::new([1; 3], &density, [-2.0; 3], [4.0; 3]),
            MassAbsorptionSpectrum::constant(0.1).unwrap(),
            None,
            medium_budget(),
        )
        .unwrap();
        assert!(bound(&scene, emitters, PathMedium::Absorbing(&medium), cx).is_none());
    });
}

#[test]
fn g0_inactive_projector_binding_does_not_acquire_a_diffuse_tail_bound() {
    with_cx(|cx, _| {
        let mut scene = empty_scene([0.0, 0.0, 3.0], [0.0, 0.0, -1.0]);
        let mut pupil = floor(0.0, 0.05, LiftedSpectrum { c: [0.0; 3] });
        pupil.material = diffuse(0.0);
        scene.primitives.push(pupil);
        let gobo = GoboTransmission::try_new(0.04, 0.02, 1, 1, vec![1.0]).unwrap();
        let projector =
            ThinLensProjector::try_new(RigidTransform::identity(), 0.1, 0.08, gobo).unwrap();
        let surfaces = [PhotometricAreaEmitter::try_new(
            0,
            1.0,
            0.0,
            SourceSpectrum::equal_energy(),
            EmissionSidedness::Front,
        )
        .unwrap()
        .with_projector(projector)
        .unwrap()];
        let emitters = PhotometricEmitters {
            points: &[],
            surfaces: &surfaces,
        };
        assert!(bound(&scene, emitters, PathMedium::Vacuum, cx).is_none());
    });
}
