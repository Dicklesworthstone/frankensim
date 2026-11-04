//! Original mixed spectral films in real native sampling and source paths.

use super::*;
use crate::conductor::{
    CoatingIorSample, CoatingLayer, CoherentAbsorbingCoating, ConductorDataStatus, ConductorSource,
};
use crate::dielectric::{DielectricError, SellmeierIor};
use crate::polarization::dielectric::scene::curved::bundle::CurvedBundleSettings;
use crate::polarization::dielectric::scene::curved::{
    CurvedPolarizedScene, CurvedTraceSettings, SceneMedium,
};

#[cfg(feature = "photometry")]
mod camera;

fn films() -> (CoherentDielectricCoating, CoherentAbsorbingCoating) {
    let source = ConductorSource::try_new(
        hash_domain("native-spectral-stack", b"original-laws"),
        ConductorDataStatus::Representative,
    )
    .unwrap();
    let silica = SellmeierIor::try_new(
        [0.696_166_3, 0.407_942_6, 0.897_479_4],
        [
            0.068_404_3_f64.powi(2),
            0.116_241_4_f64.powi(2),
            9.896_161_f64.powi(2),
        ],
    )
    .unwrap();
    let eta = [2.1, 2.05, 1.9, 1.85, 1.8, 1.9, 2.15, 2.0, 1.95];
    let k = [0.01, 0.025, 0.07, 0.12, 0.08, 0.16, 0.24, 0.11, 0.04];
    (
        CoherentDielectricCoating::try_new(silica, 73e-9, source).unwrap(),
        CoherentAbsorbingCoating::try_new(
            core::array::from_fn(|i| {
                CoatingIorSample::try_new(380.0 + 50.0 * i as f64, eta[i], k[i]).unwrap()
            }),
            83e-9,
            source,
        )
        .unwrap(),
    )
}

// Independent Python cmath recursion of interface amplitudes and exp(i*k*d),
// with normal transmitted Poynting factor n_exit/n_incident. Film phase uses
// the original Sellmeier law; the lossy eta/k uses its original interval.
// Columns: wavelength, outside R, inside R, reciprocal T, for bulk n=1,1.6.
const POWERS: [(f64, f64, f64, f64); 3] = [
    (
        430.0,
        0.016685242321602747,
        0.01512526197252545,
        0.92649839129717337,
    ),
    (
        550.0,
        0.0034498132599879048,
        0.019783434805665863,
        0.82311688411774631,
    ),
    (
        730.0,
        0.069337179319487074,
        0.092307200658427971,
        0.79443157349579197,
    ),
];

fn curved_settings() -> CurvedTraceSettings {
    CurvedTraceSettings {
        max_distance_m: 10.0,
        max_surface_events: 12,
        max_total_steps: 16_384,
        max_query_work: 1_000_000,
        ..CurvedTraceSettings::default()
    }
}

fn close(a: f64, b: f64, relative: f64) {
    assert!(
        (a - b).abs() <= relative * b.abs().max(1e-15),
        "{a:.17e} != {b:.17e}"
    );
}

#[test]
fn g1_original_stack_sampling_retains_both_media_beer_and_reverse_bundle_density() {
    with_cx(|cx| {
        let (transparent, absorbing) = films();
        let layers = [
            CoatingLayer::Dielectric(&transparent),
            CoatingLayer::Absorbing(&absorbing),
        ];
        let stack = DielectricCoatingStack::try_mixed(&layers).unwrap();
        let geometry = InstanceScene::try_new(vec![
            plane(1, 0.0, true),
            plane(2, 2.0, true),
            plane(3, -2.0, true),
        ])
        .unwrap();
        let media = [glass(1.0, 0.1), glass(1.6, 0.3)];
        let native_media = [
            SceneMedium::Homogeneous(&media[0]),
            SceneMedium::Homogeneous(&media[1]),
        ];
        let surfaces = vec![boundary(1, 0, 1), detector(2), detector(3)];
        let stacks = [DielectricStackBinding {
            object_id: 1,
            stack,
        }];
        let straight = PolarizedScene::try_new(&geometry, &media, surfaces.clone())
            .unwrap()
            .with_dielectric_coating_stacks(cx, &stacks)
            .unwrap();
        let curved = CurvedPolarizedScene::try_new(cx, &geometry, &native_media, surfaces, 100_000)
            .unwrap()
            .with_dielectric_coating_stacks(cx, &stacks)
            .unwrap();
        for (wavelength, r_out, r_in, transmission) in POWERS {
            for initial in 0..2 {
                let (reflection, opposite_r, sign) = if initial == 0 {
                    (r_out, r_in, 1.0)
                } else {
                    (r_in, r_out, -1.0)
                };
                let survival = reflection + transmission;
                let p_reflection = reflection / survival;
                let direction = Vec3::new(0.0, 0.0, -sign);
                let frame =
                    PolarizationFrame::try_new(direction, Vec3::new(1.0, 0.0, 0.0)).unwrap();
                let input = PolarizedBeam::try_new(
                    Point3::new(0.0, 0.0, sign),
                    frame,
                    wavelength,
                    media[initial],
                    Stokes::try_new(1.0, 0.3, 0.4, 0.2).unwrap(),
                )
                .unwrap();
                // R is the smaller branch for these independently evaluated
                // spectra. Each midpoint remains away from a selection threshold.
                for (draw, event) in [
                    (0.5 * p_reflection, DielectricEvent::Reflection),
                    (0.5 + 0.5 * p_reflection, DielectricEvent::Transmission),
                ] {
                    let reflected = event == DielectricEvent::Reflection;
                    let target = if reflected { initial } else { 1 - initial };
                    let endpoint = if reflected != (initial == 1) { 2 } else { 3 };
                    let expected =
                        survival * (-[0.1_f64, 0.3][initial] - 2.0 * [0.1, 0.3][target]).exp();
                    let samples = |_| [0.2, 0.7, draw];
                    let a = straight
                        .trace_with_samples(cx, input, initial, settings(), samples)
                        .unwrap();
                    let b = curved
                        .trace_with_samples(cx, input, initial, curved_settings(), samples)
                        .unwrap();
                    assert_eq!(a.termination(), PathTermination::Detected);
                    assert_eq!(b.termination(), PathTermination::Detected);
                    assert_eq!(a.medium_index(), target);
                    assert_eq!(b.medium_index(), target);
                    assert_eq!(a.last_hit().unwrap().object_id, endpoint);
                    assert_eq!(b.last_contact().unwrap().object_id, endpoint);
                    for (power, tolerance) in [(a.beam().power(), 2e-11), (b.beam().power(), 2e-6)]
                    {
                        close(power.i, expected, tolerance);
                        close(power.q, 0.3 * expected, tolerance);
                        let phase_sign = if reflected { -1.0 } else { 1.0 };
                        close(power.u, phase_sign * 0.4 * expected, tolerance);
                        close(power.v, phase_sign * 0.2 * expected, tolerance);
                    }
                    close(b.traveled_m(), 3.0, 2e-6);
                    close(
                        b.optical_path_length_m(),
                        [1.0_f64, 1.6][initial] + 2.0 * [1.0, 1.6][target],
                        2e-6,
                    );

                    let bundle = curved
                        .launch_bundle(
                            input.ray(),
                            frame.horizontal(),
                            wavelength,
                            input.power(),
                            initial,
                        )
                        .unwrap();
                    let branch = curved
                        .trace_bundle_with_branches(
                            cx,
                            bundle,
                            initial,
                            CurvedBundleSettings {
                                trace: curved_settings(),
                                ..CurvedBundleSettings::default()
                            },
                            |_, _| event,
                        )
                        .unwrap();
                    let physical_power = if reflected { reflection } else { transmission };
                    close(
                        branch.bundle().interface_log_power().exp(),
                        physical_power,
                        3e-11,
                    );
                    let reverse_pdf = if reflected {
                        p_reflection
                    } else {
                        transmission / (opposite_r + transmission)
                    };
                    close(
                        branch.bundle().reverse_branch_log_probability().exp(),
                        reverse_pdf,
                        3e-11,
                    );
                }
            }
        }
    });
}

#[test]
fn g3_original_stack_ggx_proposals_keep_polarized_reciprocity_and_frames() {
    let (transparent, absorbing) = films();
    let layers = [
        CoatingLayer::Dielectric(&transparent),
        CoatingLayer::Absorbing(&absorbing),
    ];
    let stack = DielectricCoatingStack::try_mixed(&layers).unwrap();
    let normal = Vec3::new(0.0, 0.0, 1.0);
    let surface = DielectricSurface::try_rough(0.2).unwrap();
    let mut events = [0_u32; 2];
    for (eta_i, eta_t, sign) in [(1.0, 1.6, -1.0), (1.6, 1.0, 1.0)] {
        let input = PolarizedBeam::try_new(
            Point3::new(0.0, 0.0, 0.0),
            PolarizationFrame::try_new(
                Vec3::new(0.6 * sign, 0.0, 0.8 * sign),
                Vec3::new(0.0, 1.0, 0.0),
            )
            .unwrap(),
            550.0,
            glass(eta_i, 0.0),
            Stokes::try_new(2.0, 0.4, 0.7, 0.9).unwrap(),
        )
        .unwrap();
        for draw in [[0.0, 0.0, 0.0], [0.17, 0.31, 0.0], [0.17, 0.31, 0.999]] {
            let Some(sample) = input
                .sample_stack_interface(normal, glass(eta_t, 0.0), surface, stack, draw)
                .unwrap()
            else {
                continue;
            };
            assert!(!sample.delta);
            let transmitted = sample.event == DielectricEvent::Transmission;
            events[usize::from(transmitted)] += 1;
            let wi = sample.beam.frame().direction();
            let source =
                PolarizationFrame::try_new(wi.scale(-1.0), Vec3::new(0.0, 1.0, 0.0)).unwrap();
            let response = input
                .evaluate_stack_dielectric_connection(
                    normal,
                    glass(eta_t, 0.0),
                    surface,
                    stack,
                    source,
                    Stokes::unpolarized(1.0),
                )
                .unwrap();
            close(response.pdf, sample.pdf, 2e-10);
            let factor = sample.beam.power().i / input.power().i * sample.pdf
                / normal.dot(wi).abs()
                * if transmitted {
                    (eta_i / eta_t).powi(2)
                } else {
                    1.0
                };
            close(response.factor, factor, 2e-10);
            let state = Stokes::try_new(1.0, 0.1, 0.4, -0.6).unwrap();
            let reference = input
                .evaluate_stack_dielectric_connection(
                    normal,
                    glass(eta_t, 0.0),
                    surface,
                    stack,
                    source,
                    state,
                )
                .unwrap();
            let rotated_camera = input.with_reference_axis(Vec3::new(0.3, 0.7, 0.9)).unwrap();
            let rotated_source =
                PolarizationFrame::try_new(source.direction(), Vec3::new(0.8, 0.4, 0.2)).unwrap();
            let rotated = rotated_camera
                .evaluate_stack_dielectric_connection(
                    normal,
                    glass(eta_t, 0.0),
                    surface,
                    stack,
                    rotated_source,
                    source.express_in(state, rotated_source),
                )
                .unwrap();
            close(rotated.factor, reference.factor, 2e-10);
            assert_eq!(rotated.pdf, reference.pdf);
        }
    }
    assert!(events.into_iter().all(|count| count > 0));
}

#[test]
fn g0_stack_binding_preserves_disabled_paths_and_original_domain_refusals() {
    with_cx(|cx| {
        let (transparent, absorbing) = films();
        let layers = [
            CoatingLayer::Dielectric(&transparent),
            CoatingLayer::Absorbing(&absorbing),
        ];
        let stack = DielectricCoatingStack::try_mixed(&layers).unwrap();
        let geometry =
            InstanceScene::try_new(vec![plane(1, 0.0, true), plane(2, -2.0, true)]).unwrap();
        let media = [glass(1.0, 0.1), glass(1.6, 0.3)];
        let scene = || {
            PolarizedScene::try_new(&geometry, &media, vec![boundary(1, 0, 1), detector(2)])
                .unwrap()
        };
        let expected = scene()
            .trace(cx, beam(1.0, media[0]), 0, settings(), 91, 3)
            .unwrap();
        let actual = scene()
            .with_dielectric_coating_stacks(cx, &[])
            .unwrap()
            .trace(cx, beam(1.0, media[0]), 0, settings(), 91, 3)
            .unwrap();
        assert_eq!(actual, expected);
        for bindings in [
            vec![DielectricStackBinding {
                object_id: 99,
                stack,
            }],
            vec![DielectricStackBinding {
                object_id: 2,
                stack,
            }],
            vec![
                DielectricStackBinding {
                    object_id: 1,
                    stack
                };
                2
            ],
        ] {
            assert!(matches!(
                scene().with_dielectric_coating_stacks(cx, &bindings),
                Err(SceneError::InvalidBindings)
            ));
        }
        let old_films = [transparent];
        let old = PolarizedScene::try_new_with_optics(
            &geometry,
            &media,
            &[],
            &old_films,
            vec![
                SurfaceBinding {
                    object_id: 1,
                    interaction: Interaction::CoatedDielectric {
                        outside: 0,
                        inside: 1,
                        coating: 0,
                    },
                },
                detector(2),
            ],
        )
        .unwrap();
        assert!(matches!(
            old.with_dielectric_coating_stacks(
                cx,
                &[DielectricStackBinding {
                    object_id: 1,
                    stack
                }]
            ),
            Err(SceneError::InvalidBindings)
        ));
        let empty = DielectricCoatingStack::try_new_in_band(&[], (500.0, 600.0)).unwrap();
        assert!(matches!(
            empty.interface(430.0, false),
            Err(DielectricError::InvalidWavelength)
        ));
        let restricted = scene()
            .with_dielectric_coating_stacks(
                cx,
                &[DielectricStackBinding {
                    object_id: 1,
                    stack: empty,
                }],
            )
            .unwrap();
        let input = beam(1.0, media[0]);
        let outside_band = PolarizedBeam::try_new(
            input.ray().origin,
            input.frame(),
            430.0,
            input.medium(),
            input.power(),
        )
        .unwrap();
        assert!(matches!(
            restricted.trace_with_samples(cx, outside_band, 0, settings(), transmit),
            Err(SceneError::Beam(BeamError::Dielectric(
                DielectricError::InvalidWavelength
            )))
        ));
        let null = empty
            .interface(550.0, false)
            .unwrap()
            .sample_rough(
                Vec3::new(0.0, 0.0, 1.0),
                Vec3::new(0.0, 0.0, 1.0),
                1.6,
                1.6,
                0.2,
                0.2,
                0.3,
                0.4,
            )
            .unwrap()
            .unwrap();
        assert!(null.delta);
        assert_eq!(null.radiance_weight, 1.0);
        assert!(matches!(
            stack.interface(550.0, false).unwrap().sample_rough(
                Vec3::new(0.0, 0.0, 1.0),
                Vec3::new(0.0, 0.0, 1.0),
                1.6,
                1.6,
                0.2,
                0.2,
                0.3,
                0.4
            ),
            Err(DielectricError::CoatedNullRoughInterface)
        ));
    });
}
