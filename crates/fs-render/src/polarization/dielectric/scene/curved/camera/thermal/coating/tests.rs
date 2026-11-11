//! Native polarized thermal radiance, geometry frames and exact optical binding.

use super::CoatingThermalEmitter;
use crate::charts::{Ray, TriMesh};
use crate::conductor::{
    CoatingIorSample, CoatingLayer, CoherentAbsorbingCoating, ConductorDataStatus,
    ConductorIorSample, ConductorOptics, ConductorSource, ConductorSurface,
};
use crate::dielectric::{
    BeerLambertAbsorption, DielectricGlass, DielectricIorSample, GlassProvenance, TabulatedIor,
};
use crate::gradient_index::GradientIndexOptions;
use crate::gradient_index::event::{BoundarySide, GradientEventOptions};
use crate::instances::{GeometryInstance, InstanceScene, RigidTransform, SharedGeometry};
use crate::polarization::dielectric::scene::curved::camera::{
    CameraAnalyzer, CameraRadianceError, CurvedCameraScene, EmissionSides, RadianceSpectrum,
    SurfaceEmitter,
};
use crate::polarization::dielectric::scene::curved::{
    CurvedPolarizedScene, CurvedSceneError, CurvedTraceSettings, SceneMedium,
};
use crate::polarization::dielectric::scene::{Interaction, SceneError, SurfaceBinding};
use crate::thermal_surface::CoatingThermalEmission;
use core::f64::consts::{FRAC_PI_2, FRAC_PI_4};
use fs_blake3::{ContentHash, hash_domain};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_geom::{Point3, Vec3};

const Y: Vec3 = Vec3::new(0.0, 1.0, 0.0);
const VISIBLE: (f64, f64) = (380.0, 780.0);

fn with_cx(f: impl FnOnce(&Cx<'_>, &CancelGate)) {
    let gate = CancelGate::new_clock_free();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 97,
                kernel_id: 29,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx, &gate);
    });
}

fn settings() -> CurvedTraceSettings {
    CurvedTraceSettings {
        max_surface_events: 4,
        max_distance_m: 4.0,
        max_total_steps: 16_384,
        max_query_work: 1_000_000,
        event: GradientEventOptions {
            integration: GradientIndexOptions {
                position_tolerance_m: 1e-11,
                direction_tolerance: 1e-12,
                ..GradientIndexOptions::default()
            },
            surface_tolerance_m: 1e-8,
            ..GradientEventOptions::default()
        },
    }
}

fn patch(cosine: f64, azimuth: f64) -> GeometryInstance {
    // Rotate the original body-local +z normal; the view stays along world z.
    let (s, c) = (0.5 * cosine.acos()).sin_cos();
    let placement =
        RigidTransform::try_new([-azimuth.sin() * s, azimuth.cos() * s, 0.0, c], [0.0; 3]).unwrap();
    GeometryInstance::try_new(
        9,
        hash_domain("native-coating-thermal-test", b"oriented-face"),
        SharedGeometry::mesh(TriMesh::new(
            vec![[-2.0, -2.0, 0.0], [2.0, -2.0, 0.0], [0.0, 4.0, 0.0]],
            vec![[0, 1, 2]],
        )),
        placement,
    )
    .unwrap()
}

fn optical_source(tag: u8) -> ConductorSource {
    ConductorSource::try_new(ContentHash([tag; 32]), ConductorDataStatus::Representative).unwrap()
}

fn metal(bounds: (f64, f64), tag: u8) -> ConductorOptics {
    let gap = (bounds.1 - bounds.0) / 8.0;
    ConductorOptics::try_new_in_band(
        core::array::from_fn(|i| {
            ConductorIorSample::try_new_in_band(bounds.0 + gap * i as f64, 0.3, 3.0).unwrap()
        }),
        gap,
        optical_source(tag),
    )
    .unwrap()
}

fn medium(bounds: (f64, f64), index: f64, absorption: f64) -> DielectricGlass {
    let gap = (bounds.1 - bounds.0) / 8.0;
    let phase = TabulatedIor::try_new_in_band(
        core::array::from_fn(|i| {
            DielectricIorSample::try_new_in_band(bounds.0 + gap * i as f64, index).unwrap()
        }),
        gap,
    )
    .unwrap();
    DielectricGlass::new(
        phase,
        BeerLambertAbsorption::try_constant(absorption).unwrap(),
        GlassProvenance::Custom,
    )
}

fn smooth() -> Vec<SurfaceBinding> {
    vec![SurfaceBinding {
        object_id: 9,
        interaction: Interaction::SmoothConductor {
            optics: 0,
            medium: 0,
        },
    }]
}

fn bare(optics: ConductorOptics, temperature: f64, tag: u8) -> CoatingThermalEmission {
    CoatingThermalEmission::try_opaque(&[], &[], optics, temperature, ContentHash([tag; 32]))
        .unwrap()
}

fn ray(origin_z: f64) -> Ray {
    Ray {
        origin: Point3::new(0.0, 0.0, origin_z),
        dir: Vec3::new(0.0, 0.0, -origin_z),
    }
}

fn observe(
    cx: &Cx<'_>,
    camera: &CurvedCameraScene<'_, '_>,
    origin_z: f64,
    wavelength: f64,
    analyzer: Option<CameraAnalyzer>,
) -> f64 {
    let sample = match analyzer {
        Some(analyzer) => camera.sample_ray_analyzed_with_samples(
            cx,
            ray(origin_z),
            wavelength,
            0,
            analyzer,
            settings(),
            |_| panic!("smooth thermal conductor needs no sampled surface direction"),
        ),
        None => {
            camera.sample_ray_with_samples(cx, ray(origin_z), wavelength, 0, settings(), |_| {
                panic!("smooth thermal conductor needs no sampled surface direction")
            })
        }
    }
    .unwrap();
    assert_eq!(sample.path().surface_events(), 1);
    sample.radiance_w_m2_sr_nm()
}

fn near(actual: f64, expected: f64) {
    near_relative(actual, expected, 2e-9);
}

fn near_relative(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs().max(f64::MIN_POSITIVE),
        "actual={actual:.17e}, expected={expected:.17e}, relative={relative:.3e}",
    );
}

fn planck_si(temperature: f64, wavelength_nm: f64) -> f64 {
    let wavelength = wavelength_nm * 1e-9;
    let (h, c, k) = (6.626_070_15e-34, 299_792_458.0, 1.380_649e-23);
    2.0 * h * c * c / (wavelength.powi(5) * (h * c / (wavelength * k * temperature)).exp_m1())
        * 1e-9
}

fn bare_radiance_sp(temperature: f64, wavelength: f64, cosine: f64) -> [f64; 2] {
    // Independent complex Fresnel powers for the declared eta=0.3, k=3 metal.
    // q=sqrt((eta+i*k)^2-sin(theta)^2) is the passive transmitted normal mode.
    let (eta, k) = (0.3_f64, 3.0_f64);
    let real = eta * eta - k * k;
    let imaginary = 2.0 * eta * k;
    let q_real = real - (1.0 - cosine * cosine);
    let modulus = q_real.hypot(imaginary);
    let a = (0.5 * (modulus + q_real)).sqrt();
    let b = (0.5 * (modulus - q_real)).sqrt();
    let rs = ((cosine - a).powi(2) + b * b) / ((cosine + a).powi(2) + b * b);
    let rp = ((real * cosine - a).powi(2) + (imaginary * cosine - b).powi(2))
        / ((real * cosine + a).powi(2) + (imaginary * cosine + b).powi(2));
    let per_polarization = 0.5 * planck_si(temperature, wavelength);
    [per_polarization * (1.0 - rs), per_polarization * (1.0 - rp)]
}

#[test]
fn hot_bare_metal_resolves_s_p_and_body_rotated_detector_frames() {
    with_cx(|cx, _| {
        let host = medium(VISIBLE, 1.0, 0.0);
        let media = [SceneMedium::Homogeneous(&host)];
        let metals = [metal(VISIBLE, 101)];
        let emitter =
            CoatingThermalEmitter::try_new(9, BoundarySide::Positive, bare(metals[0], 2400.0, 102))
                .unwrap();
        for (cosine, azimuth) in [(1.0, 0.0), (0.6, 0.0), (0.6, FRAC_PI_4), (0.6, FRAC_PI_2)] {
            let geometry = InstanceScene::try_new(vec![patch(cosine, azimuth)]).unwrap();
            let native = CurvedPolarizedScene::try_new_with_conductors(
                cx,
                &geometry,
                &media,
                &metals,
                smooth(),
                100_000,
            )
            .unwrap();
            let camera = CurvedCameraScene::try_new(cx, &native, &[], None)
                .unwrap()
                .with_coating_thermal_emitters(cx, core::slice::from_ref(&emitter))
                .unwrap();
            for wavelength in [450.0, 650.0] {
                let [ls, lp] = bare_radiance_sp(2400.0, wavelength, cosine);
                near(observe(cx, &camera, 1.0, wavelength, None), ls + lp);
                for angle in [0.0, FRAC_PI_4, FRAC_PI_2] {
                    let analyzer = CameraAnalyzer::try_linear(Y, angle).unwrap();
                    let relative = azimuth - angle;
                    near(
                        observe(cx, &camera, 1.0, wavelength, Some(analyzer)),
                        ls * relative.cos().powi(2) + lp * relative.sin().powi(2),
                    );
                }
                for right_handed in [false, true] {
                    near(
                        observe(
                            cx,
                            &camera,
                            1.0,
                            wavelength,
                            Some(CameraAnalyzer::try_circular(Y, right_handed).unwrap()),
                        ),
                        0.5 * (ls + lp),
                    );
                }
            }
        }
    });
}

#[test]
fn nonisothermal_absorbing_film_keeps_original_layer_emission_and_beer_path() {
    with_cx(|cx, _| {
        let cosine = 0.6;
        let geometry = InstanceScene::try_new(vec![patch(cosine, 0.0)]).unwrap();
        let host = medium(VISIBLE, 1.0, 0.4);
        let media = [SceneMedium::Homogeneous(&host)];
        let bulk = metal(VISIBLE, 103);
        let film = CoherentAbsorbingCoating::try_new(
            core::array::from_fn(|i| {
                CoatingIorSample::try_new(380.0 + 50.0 * i as f64, 1.6, 0.35).unwrap()
            }),
            90e-9,
            optical_source(104),
        )
        .unwrap();
        let metals = [bulk.with_absorbing_coating(film)];
        let native = CurvedPolarizedScene::try_new_with_conductors(
            cx,
            &geometry,
            &media,
            &metals,
            smooth(),
            100_000,
        )
        .unwrap();
        // The camera ray stops at its contact with the patch, not on the patch.
        // Arrival is admitted once the witness lies within surface_tolerance_m
        // of the face (gradient_index/event.rs:519). Along this ray, at approach
        // cosine `cosine`, the stop is therefore at most surface_tolerance_m / cosine
        // short of the face. The host's Beer factor covers the path up to the stop,
        // so the radiance can exceed state.i * exp(-0.4) by at most
        // expm1(0.4 * surface_tolerance_m / cosine). That is 6.67e-9 relative here;
        // R-c measured 3.31e-9, from a stop 8.27e-9 m short. The file's 2e-9
        // allowance for everything else is kept on top of it.
        let contact_bound_m = settings().event.surface_tolerance_m / cosine;
        let tolerance = 2e-9 + (0.4 * contact_bound_m).exp_m1();
        for (film_temperature, substrate_temperature) in [(1200.0, 2500.0), (2500.0, 1200.0)] {
            let source = CoatingThermalEmission::try_opaque(
                &[CoatingLayer::Absorbing(&film)],
                &[film_temperature],
                bulk,
                substrate_temperature,
                ContentHash([105; 32]),
            )
            .unwrap();
            let emitter =
                CoatingThermalEmitter::try_new(9, BoundarySide::Positive, source.clone()).unwrap();
            let camera = CurvedCameraScene::try_new(cx, &native, &[], None)
                .unwrap()
                .with_coating_thermal_emitters(cx, &[emitter])
                .unwrap();
            for wavelength in [450.0, 513.25, 650.0] {
                let original = source.evaluate(wavelength, cosine).unwrap();
                let state = original.stokes();
                let attenuation = (-0.4_f64).exp();
                near_relative(
                    observe(cx, &camera, 1.0, wavelength, None),
                    state.i * attenuation,
                    tolerance,
                );
                for (angle, sign) in [(0.0, 1.0), (FRAC_PI_2, -1.0)] {
                    near_relative(
                        observe(
                            cx,
                            &camera,
                            1.0,
                            wavelength,
                            Some(CameraAnalyzer::try_linear(Y, angle).unwrap()),
                        ),
                        0.5 * (state.i + sign * state.q) * attenuation,
                        tolerance,
                    );
                }
                assert!(
                    original.film_radiance_w_m2_sr_nm()[0]
                        .iter()
                        .all(|value| *value > 0.0)
                );
                assert!(
                    original
                        .substrate_radiance_w_m2_sr_nm()
                        .unwrap()
                        .iter()
                        .all(|value| *value > 0.0)
                );
            }
        }
    });
}

#[test]
fn physical_faces_remain_independent_and_incompatible_bindings_refuse() {
    with_cx(|cx, gate| {
        let geometry = InstanceScene::try_new(vec![patch(1.0, 0.0)]).unwrap();
        let host = medium(VISIBLE, 1.0, 0.0);
        let media = [SceneMedium::Homogeneous(&host)];
        let metals = [metal(VISIBLE, 106)];
        let native = CurvedPolarizedScene::try_new_with_conductors(
            cx,
            &geometry,
            &media,
            &metals,
            smooth(),
            100_000,
        )
        .unwrap();
        let front_source = bare(metals[0], 2400.0, 107);
        let back_source = bare(metals[0], 1400.0, 108);
        let front = CoatingThermalEmitter::try_new(9, BoundarySide::Positive, front_source.clone())
            .unwrap();
        let back =
            CoatingThermalEmitter::try_new(9, BoundarySide::Negative, back_source.clone()).unwrap();
        assert_eq!(front.object_id(), 9);
        assert_eq!(front.side(), BoundarySide::Positive);
        assert_eq!(front.source().identity(), front_source.identity());
        let empty = || CurvedCameraScene::try_new(cx, &native, &[], None).unwrap();
        let camera = empty()
            .with_coating_thermal_emitters(cx, &[back.clone(), front.clone()])
            .unwrap();
        for (origin, source) in [(1.0, &front_source), (-1.0, &back_source)] {
            near(
                observe(cx, &camera, origin, 550.0, None),
                source.evaluate(550.0, 1.0).unwrap().stokes().i,
            );
        }
        let one_face = empty()
            .with_coating_thermal_emitters(cx, core::slice::from_ref(&front))
            .unwrap();
        assert_eq!(observe(cx, &one_face, -1.0, 550.0, None), 0.0);
        let cleared = camera.with_coating_thermal_emitters(cx, &[]).unwrap();
        assert_eq!(observe(cx, &cleared, 1.0, 550.0, None), 0.0);
        assert!(
            CoatingThermalEmitter::try_new(9, BoundarySide::OnSurface, front_source.clone())
                .is_err()
        );
        // Identical n+k numbers with a different original optical source are
        // still a different binding; this is not a tolerance-based Fresnel fit.
        let wrong_source = CoatingThermalEmitter::try_new(
            9,
            BoundarySide::Positive,
            bare(metal(VISIBLE, 109), 2400.0, 107),
        )
        .unwrap();
        let unknown =
            CoatingThermalEmitter::try_new(10, BoundarySide::Positive, front_source.clone())
                .unwrap();
        for emitters in [
            vec![front.clone(), front.clone()],
            vec![wrong_source],
            vec![unknown],
        ] {
            assert!(
                empty()
                    .with_coating_thermal_emitters(cx, &emitters)
                    .is_err()
            );
        }
        let calibrated = [SurfaceEmitter::new(
            9,
            RadianceSpectrum::try_constant(1.0).unwrap(),
            EmissionSides::Both,
        )];
        assert!(
            CurvedCameraScene::try_new(cx, &native, &calibrated, None)
                .unwrap()
                .with_coating_thermal_emitters(cx, core::slice::from_ref(&front))
                .is_err()
        );
        let rough = CurvedPolarizedScene::try_new_with_conductors(
            cx,
            &geometry,
            &media,
            &metals,
            vec![SurfaceBinding {
                object_id: 9,
                interaction: Interaction::RoughConductor {
                    optics: 0,
                    medium: 0,
                    surface: ConductorSurface::try_rough(0.2).unwrap(),
                },
            }],
            100_000,
        )
        .unwrap();
        assert!(
            CurvedCameraScene::try_new(cx, &rough, &[], None)
                .unwrap()
                .with_coating_thermal_emitters(cx, core::slice::from_ref(&front))
                .is_err()
        );
        let refractive_host = medium(VISIBLE, 1.2, 0.0);
        let refractive_media = [SceneMedium::Homogeneous(&refractive_host)];
        let refractive = CurvedPolarizedScene::try_new_with_conductors(
            cx,
            &geometry,
            &refractive_media,
            &metals,
            smooth(),
            100_000,
        )
        .unwrap();
        let in_material = CurvedCameraScene::try_new(cx, &refractive, &[], None)
            .unwrap()
            .with_coating_thermal_emitters(cx, core::slice::from_ref(&front))
            .unwrap();
        assert!(matches!(
            in_material.sample_ray_with_samples(cx, ray(1.0), 550.0, 0, settings(), |_| {
                panic!("smooth conductor")
            }),
            Err(CameraRadianceError::UnsupportedSurface { object_id: 9 })
        ));
        let cancelled = empty();
        gate.request();
        assert!(matches!(
            cancelled.with_coating_thermal_emitters(cx, &[front]),
            Err(CameraRadianceError::Scene(CurvedSceneError::Scene(
                SceneError::Cancelled
            )))
        ));
    });
}

#[test]
fn room_temperature_infrared_uses_its_original_supported_planck_band() {
    with_cx(|cx, _| {
        let bounds = (8000.0, 16_000.0);
        let cosine = 0.6;
        let geometry = InstanceScene::try_new(vec![patch(cosine, 0.0)]).unwrap();
        let host = medium(bounds, 1.0, 0.0);
        let media = [SceneMedium::Homogeneous(&host)];
        let metals = [metal(bounds, 110)];
        let native = CurvedPolarizedScene::try_new_with_conductors(
            cx,
            &geometry,
            &media,
            &metals,
            smooth(),
            100_000,
        )
        .unwrap();
        let emitter =
            CoatingThermalEmitter::try_new(9, BoundarySide::Positive, bare(metals[0], 300.0, 111))
                .unwrap();
        let camera = CurvedCameraScene::try_new(cx, &native, &[], None)
            .unwrap()
            .with_coating_thermal_emitters(cx, &[emitter])
            .unwrap();
        for wavelength in [9750.25, 14_000.0] {
            let [ls, lp] = bare_radiance_sp(300.0, wavelength, cosine);
            near(observe(cx, &camera, 1.0, wavelength, None), ls + lp);
            near(
                observe(
                    cx,
                    &camera,
                    1.0,
                    wavelength,
                    Some(CameraAnalyzer::try_linear(Y, 0.0).unwrap()),
                ),
                ls,
            );
            near(
                observe(
                    cx,
                    &camera,
                    1.0,
                    wavelength,
                    Some(CameraAnalyzer::try_linear(Y, FRAC_PI_2).unwrap()),
                ),
                lp,
            );
        }
        assert!(
            camera
                .sample_ray_with_samples(cx, ray(1.0), 550.0, 0, settings(), |_| {
                    panic!("unsupported wavelength must not sample a surface")
                })
                .is_err()
        );
    });
}

#[test]
fn positive_thermal_channel_contraction_retains_weak_signal_and_exact_extinction() {
    use crate::polarization::Stokes;
    use crate::polarization::dielectric::scene::curved::camera::emission::thermal_polarization_weight;
    use crate::polarization::dielectric::transport::{PolarizationFrame, PolarizedBeam};

    let host = medium(VISIBLE, 1.0, 0.0);
    let weak = 1e-20;
    assert_eq!(0.5 * ((weak + 1.0) + (weak - 1.0)), 0.0);
    for azimuth in [0.0_f64, FRAC_PI_4] {
        let axis = Vec3::new(-azimuth.sin(), azimuth.cos(), 0.0);
        let camera_frame = PolarizationFrame::try_new(Vec3::new(0.0, 0.0, -1.0), axis).unwrap();
        let source_frame = PolarizationFrame::try_new(Vec3::new(0.0, 0.0, 1.0), axis).unwrap();
        for (pol, q) in [(0, 0.5), (1, -0.5)] {
            let beam = PolarizedBeam::try_new(
                Point3::new(0.0, 0.0, 0.0),
                camera_frame,
                550.0,
                host,
                Stokes {
                    i: 0.5,
                    q,
                    u: 0.0,
                    v: 0.0,
                },
            )
            .unwrap();
            for channels in [[weak, 1.0], [1.0, weak], [0.0, 1.0], [1.0, 0.0]] {
                let weight = thermal_polarization_weight(beam, source_frame, channels).unwrap();
                let measured = beam.power().i * (channels[0] + channels[1]) * weight;
                if channels[pol] == 0.0 {
                    assert_eq!(measured, 0.0);
                } else {
                    near(measured, channels[pol]);
                }
            }
        }
    }
}
