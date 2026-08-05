//! G0/G3/G4/G5 acceptance battery for deterministic Euler modal synthesis.

use core::f64::consts::TAU;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_blake3::{ContentHash, hash_domain};
use fs_euler_disc_e2e::{
    MAX_MODAL_SPATIAL_PARTICIPATION, ModalComponentValues, ModalCouplingClass, ModalDriveFrame,
    ModalPresetAuthority, ModalSpatialParticipation, ModalSynthesisBudget, ModalSynthesisError,
    ModalSynthesisModel, ModalSynthesisModelInput, RepresentativeDiscMaterial,
    representative_modal_preset,
};
use fs_evidence::cinematic_sound::{SoundModalComponent, SoundMode, SoundModeParticipation};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_math::det;

fn with_cx<R>(cancelled: bool, operation: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    if cancelled {
        gate.request();
    }
    let pool = ArenaPool::new(ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 0x4d4f_4441_4c53_594e,
                kernel_id: 0x4555_4c45,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        operation(&cx)
    })
}

fn identity(label: &str) -> ContentHash {
    hash_domain("org.frankensim.test.modal-synthesis.v1", label.as_bytes())
}

fn generous_budget(total_frames: u64) -> ModalSynthesisBudget {
    ModalSynthesisBudget {
        maximum_total_sample_frames: total_frames,
        maximum_chunk_sample_frames: usize::try_from(total_frames).unwrap(),
        maximum_abs_displacement_m: 10.0,
        maximum_abs_velocity_m_per_s: 100_000.0,
        maximum_mode_energy_j: 1.0e12,
        maximum_total_energy_j: 1.0e13,
        maximum_abs_output_fs: 1.0e9,
    }
}

fn mode(
    mode_id: u32,
    component: SoundModalComponent,
    frequency_hz: f64,
    damping_ratio: f64,
) -> SoundMode {
    SoundMode {
        mode_id,
        component,
        frequency_hz,
        damping_ratio,
        modal_mass_kg: 0.2,
        source_participation: match component {
            SoundModalComponent::Disc => SoundModeParticipation {
                disc: 1.0,
                glass_plate: 0.0,
                base_assembly: 0.0,
            },
            SoundModalComponent::GlassPlate => SoundModeParticipation {
                disc: 0.0,
                glass_plate: 1.0,
                base_assembly: 0.0,
            },
            SoundModalComponent::BaseAssembly => SoundModeParticipation {
                disc: 0.0,
                glass_plate: 0.0,
                base_assembly: 1.0,
            },
        },
        radiation_gain_fs_s_per_m: 0.1,
        material_identity: identity("material"),
        base_identity: identity("base"),
    }
}

fn build(modes: Vec<SoundMode>, budget: ModalSynthesisBudget) -> ModalSynthesisModel {
    with_cx(false, |cx| {
        ModalSynthesisModel::try_new(
            ModalSynthesisModelInput {
                sample_rate_hz: 48_000,
                modes,
                budget,
            },
            cx,
        )
        .unwrap()
    })
}

fn drive(force: ModalComponentValues, impulse: ModalComponentValues) -> ModalDriveFrame {
    ModalDriveFrame {
        localized_generalized_force_n: force,
        distributed_generalized_force_n: ModalComponentValues::ZERO,
        localized_boundary_impulse_n_s: impulse,
        distributed_boundary_impulse_n_s: ModalComponentValues::ZERO,
    }
}

fn assert_close(actual: f64, expected: f64, relative: f64, absolute: f64) {
    let tolerance = absolute.max(relative * expected.abs().max(actual.abs()));
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:.17e} expected={expected:.17e} tolerance={tolerance:.17e}"
    );
}

#[derive(Clone, Copy)]
struct Complex {
    re: f64,
    im: f64,
}

impl Complex {
    fn add(self, other: Self) -> Self {
        Self {
            re: self.re + other.re,
            im: self.im + other.im,
        }
    }

    fn sub(self, other: Self) -> Self {
        Self {
            re: self.re - other.re,
            im: self.im - other.im,
        }
    }

    fn mul(self, other: Self) -> Self {
        Self {
            re: self.re * other.re - self.im * other.im,
            im: self.re * other.im + self.im * other.re,
        }
    }

    fn scale(self, value: f64) -> Self {
        Self {
            re: self.re * value,
            im: self.im * value,
        }
    }

    fn div(self, other: Self) -> Self {
        let denominator = other.re * other.re + other.im * other.im;
        Self {
            re: (self.re * other.re + self.im * other.im) / denominator,
            im: (self.im * other.re - self.re * other.im) / denominator,
        }
    }

    fn magnitude(self) -> f64 {
        det::sqrt(self.re * self.re + self.im * self.im)
    }
}

#[test]
fn g0_impulse_step_and_harmonic_match_independent_oracles() {
    let frequency_hz = 750.0;
    let damping = 0.05;
    let modal_mass = 0.2;
    let impulse_n_s = 1.0e-3;
    let model = build(
        vec![mode(1, SoundModalComponent::Disc, frequency_hz, damping)],
        generous_budget(60_000),
    );
    let initial = with_cx(false, |cx| model.initial_checkpoint(cx).unwrap());
    let impulse_frame = drive(
        ModalComponentValues::ZERO,
        ModalComponentValues {
            disc: impulse_n_s,
            ..ModalComponentValues::ZERO
        },
    );
    let impulse = with_cx(false, |cx| {
        model
            .synthesize_chunk(
                &initial,
                &[impulse_frame],
                ModalSpatialParticipation::Declared,
                cx,
            )
            .unwrap()
    });
    let dt = 1.0 / 48_000.0;
    let omega = TAU * frequency_hz;
    let damped = omega * det::sqrt((1.0 - damping) * (1.0 + damping));
    let decay = det::exp(-damping * omega * dt);
    let expected_q = impulse_n_s / modal_mass * decay * det::sin(damped * dt) / damped;
    let expected_v = impulse_n_s / modal_mass
        * decay
        * (det::cos(damped * dt) - damping * omega / damped * det::sin(damped * dt));
    let state = impulse.successor.states()[0];
    assert_close(state.displacement_m, expected_q, 2.0e-13, 1.0e-18);
    assert_close(state.velocity_m_per_s, expected_v, 2.0e-13, 1.0e-18);
    assert_close(
        impulse.diagnostics.maximum_total_modal_energy_j,
        impulse_n_s * impulse_n_s / (2.0 * modal_mass),
        2.0e-15,
        1.0e-20,
    );

    let step_force_n = 0.2;
    let step_count = 257;
    let step_drive = vec![
        drive(
            ModalComponentValues {
                disc: step_force_n,
                ..ModalComponentValues::ZERO
            },
            ModalComponentValues::ZERO,
        );
        step_count
    ];
    let step = with_cx(false, |cx| {
        model
            .synthesize_chunk(
                &initial,
                &step_drive,
                ModalSpatialParticipation::Declared,
                cx,
            )
            .unwrap()
    });
    let time = step_count as f64 * dt;
    let decay = det::exp(-damping * omega * time);
    let transient = det::cos(damped * time) + damping * omega / damped * det::sin(damped * time);
    let expected_q = step_force_n / (modal_mass * omega * omega) * (1.0 - decay * transient);
    let expected_v = step_force_n / modal_mass * decay * det::sin(damped * time) / damped;
    let state = step.successor.states()[0];
    assert_close(state.displacement_m, expected_q, 2.0e-12, 1.0e-17);
    assert_close(state.velocity_m_per_s, expected_v, 2.0e-12, 1.0e-15);
    assert_close(model.quality_factor(1).unwrap(), 10.0, 0.0, 1.0e-15);

    let natural_hz = 800.0;
    let drive_hz = 600.0;
    let damping = 0.08;
    let force_peak = 0.1;
    let mut harmonic_mode = mode(1, SoundModalComponent::Disc, natural_hz, damping);
    harmonic_mode.modal_mass_kg = 0.25;
    harmonic_mode.radiation_gain_fs_s_per_m = 0.2;
    let harmonic_model = build(vec![harmonic_mode], generous_budget(48_000));
    let harmonic_initial = with_cx(false, |cx| harmonic_model.initial_checkpoint(cx).unwrap());
    let phase = TAU * drive_hz / 48_000.0;
    let harmonic_drive: Vec<_> = (0..48_000)
        .map(|index| {
            drive(
                ModalComponentValues {
                    disc: force_peak * det::cos(phase * index as f64),
                    ..ModalComponentValues::ZERO
                },
                ModalComponentValues::ZERO,
            )
        })
        .collect();
    let harmonic = with_cx(false, |cx| {
        harmonic_model
            .synthesize_chunk(
                &harmonic_initial,
                &harmonic_drive,
                ModalSpatialParticipation::Declared,
                cx,
            )
            .unwrap()
    });

    // Independent exact discrete ZOH transfer at z = exp(i phase).
    let omega = TAU * natural_hz;
    let nu = omega * det::sqrt((1.0 - damping) * (1.0 + damping));
    let decay = det::exp(-damping * omega / 48_000.0);
    let s = decay * det::sin(nu / 48_000.0) / nu;
    let c = decay * det::cos(nu / 48_000.0);
    let a00 = c + damping * omega * s;
    let a01 = s;
    let a10 = -omega * omega * s;
    let a11 = c - damping * omega * s;
    let b0 = (1.0 - a00) / (omega * omega * harmonic_mode.modal_mass_kg);
    let b1 = s / harmonic_mode.modal_mass_kg;
    let z = Complex {
        re: det::cos(phase),
        im: det::sin(phase),
    };
    let z_minus_a00 = z.sub(Complex { re: a00, im: 0.0 });
    let z_minus_a11 = z.sub(Complex { re: a11, im: 0.0 });
    let determinant = z_minus_a00.mul(z_minus_a11).sub(Complex {
        re: a01 * a10,
        im: 0.0,
    });
    let velocity_numerator = z_minus_a00.scale(b1).add(Complex {
        re: a10 * b0,
        im: 0.0,
    });
    let expected_peak = velocity_numerator
        .div(determinant)
        .scale(force_peak * harmonic_mode.radiation_gain_fs_s_per_m)
        .magnitude();
    let analysis_start = 40_000;
    let analysis = &harmonic.mixed_samples_fs[analysis_start..];
    let mut measured_re = 0.0;
    let mut measured_im = 0.0;
    for (offset, sample) in analysis.iter().enumerate() {
        let published_index = analysis_start + offset + 1;
        measured_re += sample * det::cos(phase * published_index as f64);
        measured_im -= sample * det::sin(phase * published_index as f64);
    }
    let scale = 2.0 / analysis.len() as f64;
    let measured_peak = det::sqrt(
        (measured_re * scale) * (measured_re * scale)
            + (measured_im * scale) * (measured_im * scale),
    );
    assert_close(measured_peak, expected_peak, 2.0e-10, 1.0e-14);
}

#[test]
fn g0_zero_critical_overdamped_and_long_decay_energy_are_stable() {
    for damping in [0.0, 1.0, 8.0] {
        let model = build(
            vec![mode(1, SoundModalComponent::Disc, 350.0, damping)],
            generous_budget(20_000),
        );
        let initial = with_cx(false, |cx| model.initial_checkpoint(cx).unwrap());
        let mut frames = vec![ModalDriveFrame::default(); 20_000];
        frames[0].localized_boundary_impulse_n_s.disc = 1.0e-3;
        let result = with_cx(false, |cx| {
            model
                .synthesize_chunk(&initial, &frames, ModalSpatialParticipation::Declared, cx)
                .unwrap()
        });
        assert!(
            result
                .mixed_samples_fs
                .iter()
                .all(|value| value.is_finite())
        );
        assert!(
            result
                .total_modal_energy_j
                .iter()
                .all(|value| value.is_finite() && *value >= 0.0)
        );
        let first = result.total_modal_energy_j[0];
        if damping == 0.0 {
            let (minimum, maximum) = result
                .total_modal_energy_j
                .iter()
                .copied()
                .fold((f64::INFINITY, 0.0_f64), |(minimum, maximum), value| {
                    (minimum.min(value), maximum.max(value))
                });
            assert!((maximum - minimum) <= 2.0e-10 * first);
            assert!(model.quality_factor(1).unwrap().is_infinite());
        } else {
            for pair in result.total_modal_energy_j.windows(2) {
                assert!(pair[1] <= pair[0] + 2.0e-13 * first);
            }
            assert!(result.total_modal_energy_j.last().unwrap() < &first);
        }
    }
}

#[test]
fn g3_nyquist_nonfinite_and_state_bounds_refuse_atomically() {
    let guard_hz: f64 = 48_000.0 * 0.5 * 0.9;
    let predecessor = f64::from_bits(guard_hz.to_bits() - 1);
    assert!(
        with_cx(false, |cx| {
            ModalSynthesisModel::try_new(
                ModalSynthesisModelInput {
                    sample_rate_hz: 48_000,
                    modes: vec![mode(1, SoundModalComponent::Disc, predecessor, 0.02)],
                    budget: generous_budget(2),
                },
                cx,
            )
        })
        .is_ok()
    );
    assert_eq!(
        with_cx(false, |cx| ModalSynthesisModel::try_new(
            ModalSynthesisModelInput {
                sample_rate_hz: 48_000,
                modes: vec![mode(1, SoundModalComponent::Disc, guard_hz, 0.02)],
                budget: generous_budget(2),
            },
            cx,
        ))
        .unwrap_err(),
        ModalSynthesisError::InvalidMode {
            mode_id: 1,
            field: "frequency_hz",
        }
    );

    let mut bounded = generous_budget(4);
    bounded.maximum_abs_velocity_m_per_s = 1.0e-4;
    let model = build(
        vec![mode(1, SoundModalComponent::Disc, 500.0, 0.02)],
        bounded,
    );
    let checkpoint = with_cx(false, |cx| model.initial_checkpoint(cx).unwrap());
    let original = checkpoint.clone();
    let excessive = drive(
        ModalComponentValues::ZERO,
        ModalComponentValues {
            disc: 1.0,
            ..ModalComponentValues::ZERO
        },
    );
    assert!(matches!(
        with_cx(false, |cx| model.synthesize_chunk(
            &checkpoint,
            &[excessive],
            ModalSpatialParticipation::Declared,
            cx,
        )),
        Err(ModalSynthesisError::LimitExceeded {
            field: "modal velocity",
            ..
        })
    ));
    assert_eq!(checkpoint, original);

    // A highly overdamped transition can hide a large instantaneous kick by
    // the frame's end. The left-boundary state must still obey the budget.
    let mut kick_bounded = generous_budget(1);
    kick_bounded.maximum_abs_velocity_m_per_s = 1.0e-4;
    let kick_model = build(
        vec![mode(1, SoundModalComponent::Disc, predecessor, 16.0)],
        kick_bounded,
    );
    let kick_checkpoint = with_cx(false, |cx| kick_model.initial_checkpoint(cx).unwrap());
    assert!(matches!(
        with_cx(false, |cx| kick_model.synthesize_chunk(
            &kick_checkpoint,
            &[drive(
                ModalComponentValues::ZERO,
                ModalComponentValues {
                    disc: 1.0e-3,
                    ..ModalComponentValues::ZERO
                },
            )],
            ModalSpatialParticipation::Declared,
            cx,
        )),
        Err(ModalSynthesisError::LimitExceeded {
            sample_frame: 0,
            mode_id: Some(1),
            field: "modal velocity",
            ..
        })
    ));

    let mut aggregate_bounded = generous_budget(1);
    aggregate_bounded.maximum_mode_energy_j = 1.0;
    aggregate_bounded.maximum_total_energy_j = 1.0e-7;
    let aggregate_model = build(
        vec![mode(1, SoundModalComponent::Disc, 500.0, 0.02)],
        aggregate_bounded,
    );
    let aggregate_checkpoint = with_cx(false, |cx| aggregate_model.initial_checkpoint(cx).unwrap());
    assert!(matches!(
        with_cx(false, |cx| aggregate_model.synthesize_chunk(
            &aggregate_checkpoint,
            &[drive(
                ModalComponentValues::ZERO,
                ModalComponentValues {
                    disc: 1.0e-3,
                    ..ModalComponentValues::ZERO
                },
            )],
            ModalSpatialParticipation::Declared,
            cx,
        )),
        Err(ModalSynthesisError::LimitExceeded {
            sample_frame: 0,
            mode_id: None,
            field: "total modal energy",
            ..
        })
    ));

    let invalid = drive(
        ModalComponentValues {
            disc: f64::NAN,
            ..ModalComponentValues::ZERO
        },
        ModalComponentValues::ZERO,
    );
    assert_eq!(
        with_cx(false, |cx| model.synthesize_chunk(
            &checkpoint,
            &[invalid],
            ModalSpatialParticipation::Declared,
            cx,
        )),
        Err(ModalSynthesisError::NonFiniteDrive {
            frame: 0,
            field: "localized disc force",
        })
    );
}

#[test]
fn g3_mode_permutation_duplicate_degeneracy_and_zero_coupling_are_explicit() {
    let first = mode(7, SoundModalComponent::GlassPlate, 1_200.0, 0.03);
    let second = mode(2, SoundModalComponent::Disc, 500.0, 0.01);
    let model_a = build(vec![first, second], generous_budget(32));
    let model_b = build(vec![second, first], generous_budget(32));
    assert_eq!(model_a.identity(), model_b.identity());
    assert_eq!(model_a.modes(), model_b.modes());
    let frames = vec![
        drive(
            ModalComponentValues {
                disc: 0.1,
                glass_plate: -0.2,
                base_assembly: 0.0,
            },
            ModalComponentValues::ZERO,
        );
        32
    ];
    let run = |model: &ModalSynthesisModel| {
        let checkpoint = with_cx(false, |cx| model.initial_checkpoint(cx).unwrap());
        with_cx(false, |cx| {
            model
                .synthesize_chunk(
                    &checkpoint,
                    &frames,
                    ModalSpatialParticipation::Declared,
                    cx,
                )
                .unwrap()
        })
    };
    assert_eq!(run(&model_a), run(&model_b));

    assert_eq!(
        with_cx(false, |cx| ModalSynthesisModel::try_new(
            ModalSynthesisModelInput {
                sample_rate_hz: 48_000,
                modes: vec![second, second],
                budget: generous_budget(2),
            },
            cx,
        ))
        .unwrap_err(),
        ModalSynthesisError::DuplicateModeId(2)
    );

    let mut degenerate = second;
    degenerate.mode_id = 3;
    let degenerate_model = build(vec![second, degenerate], generous_budget(2));
    assert_eq!(degenerate_model.modes().len(), 2);
    let checkpoint = with_cx(false, |cx| degenerate_model.initial_checkpoint(cx).unwrap());
    let excited = with_cx(false, |cx| {
        degenerate_model
            .synthesize_chunk(
                &checkpoint,
                &[drive(
                    ModalComponentValues::ZERO,
                    ModalComponentValues {
                        disc: 1.0e-3,
                        ..ModalComponentValues::ZERO
                    },
                )],
                ModalSpatialParticipation::Declared,
                cx,
            )
            .unwrap()
    });
    assert_ne!(excited.mixed_samples_fs[0], 0.0);

    let mut uncoupled = second;
    uncoupled.source_participation = SoundModeParticipation {
        disc: 0.0,
        glass_plate: 0.0,
        base_assembly: 0.0,
    };
    let silent_model = build(vec![uncoupled], generous_budget(4));
    let silent_checkpoint = with_cx(false, |cx| silent_model.initial_checkpoint(cx).unwrap());
    let silent = with_cx(false, |cx| {
        silent_model
            .synthesize_chunk(
                &silent_checkpoint,
                &[drive(
                    ModalComponentValues {
                        disc: 1.0,
                        glass_plate: 1.0,
                        base_assembly: 1.0,
                    },
                    ModalComponentValues {
                        disc: 1.0,
                        glass_plate: 1.0,
                        base_assembly: 1.0,
                    },
                )],
                ModalSpatialParticipation::Declared,
                cx,
            )
            .unwrap()
    });
    assert_eq!(silent.mixed_samples_fs, vec![0.0]);
    assert_eq!(silent.total_modal_energy_j, vec![0.0]);
}

#[test]
fn g3_declared_cross_routing_and_source_location_factors_are_observable() {
    let glass = mode(1, SoundModalComponent::GlassPlate, 900.0, 0.02);
    let independent = build(vec![glass], generous_budget(4));
    assert_eq!(independent.coupling(), ModalCouplingClass::Independent);
    let initial = with_cx(false, |cx| independent.initial_checkpoint(cx).unwrap());
    let disc_impulse = drive(
        ModalComponentValues::ZERO,
        ModalComponentValues {
            disc: 1.0e-3,
            ..ModalComponentValues::ZERO
        },
    );
    let silent = with_cx(false, |cx| {
        independent
            .synthesize_chunk(
                &initial,
                &[disc_impulse],
                ModalSpatialParticipation::Declared,
                cx,
            )
            .unwrap()
    });
    assert_eq!(silent.mixed_samples_fs, vec![0.0]);

    let mut cross_mode = glass;
    cross_mode.source_participation.disc = 0.5;
    let coupled = build(vec![cross_mode], generous_budget(4));
    assert_eq!(
        coupled.coupling(),
        ModalCouplingClass::DeclaredCrossParticipation
    );
    let initial = with_cx(false, |cx| coupled.initial_checkpoint(cx).unwrap());
    let positive = with_cx(false, |cx| {
        coupled
            .synthesize_chunk(
                &initial,
                &[disc_impulse],
                ModalSpatialParticipation::PerFrameModeFactors(&[1.0]),
                cx,
            )
            .unwrap()
    });
    let negative = with_cx(false, |cx| {
        coupled
            .synthesize_chunk(
                &initial,
                &[disc_impulse],
                ModalSpatialParticipation::PerFrameModeFactors(&[-1.0]),
                cx,
            )
            .unwrap()
    });
    assert_ne!(positive.stem_frames[0].glass_plate_fs, 0.0);
    assert_eq!(
        positive.mixed_samples_fs[0].to_bits(),
        (-negative.mixed_samples_fs[0]).to_bits()
    );
    assert_eq!(
        with_cx(false, |cx| coupled.synthesize_chunk(
            &initial,
            &[disc_impulse],
            ModalSpatialParticipation::PerFrameModeFactors(
                &[MAX_MODAL_SPATIAL_PARTICIPATION + 1.0]
            ),
            cx,
        )),
        Err(ModalSynthesisError::InvalidSpatialParticipation {
            frame: 0,
            mode_id: 1,
        })
    );
}

#[test]
fn g3_spatial_factor_modulates_localized_but_not_distributed_drive() {
    let model = build(
        vec![mode(1, SoundModalComponent::Disc, 700.0, 0.02)],
        generous_budget(1),
    );
    let initial = with_cx(false, |cx| model.initial_checkpoint(cx).unwrap());
    let force = ModalComponentValues {
        disc: 0.25,
        ..ModalComponentValues::ZERO
    };
    let impulse = ModalComponentValues {
        disc: 1.0e-5,
        ..ModalComponentValues::ZERO
    };

    let localized = drive(force, impulse);
    let localized_suppressed = with_cx(false, |cx| {
        model
            .synthesize_chunk(
                &initial,
                &[localized],
                ModalSpatialParticipation::PerFrameModeFactors(&[0.0]),
                cx,
            )
            .unwrap()
    });
    assert_eq!(localized_suppressed.mixed_samples_fs, vec![0.0]);
    assert_eq!(localized_suppressed.total_modal_energy_j, vec![0.0]);

    let distributed = ModalDriveFrame {
        localized_generalized_force_n: ModalComponentValues::ZERO,
        distributed_generalized_force_n: force,
        localized_boundary_impulse_n_s: ModalComponentValues::ZERO,
        distributed_boundary_impulse_n_s: impulse,
    };
    let distributed_declared = with_cx(false, |cx| {
        model
            .synthesize_chunk(
                &initial,
                &[distributed],
                ModalSpatialParticipation::Declared,
                cx,
            )
            .unwrap()
    });
    let distributed_with_zero_factor = with_cx(false, |cx| {
        model
            .synthesize_chunk(
                &initial,
                &[distributed],
                ModalSpatialParticipation::PerFrameModeFactors(&[0.0]),
                cx,
            )
            .unwrap()
    });
    assert_ne!(distributed_declared.mixed_samples_fs, vec![0.0]);
    assert_eq!(distributed_with_zero_factor, distributed_declared);

    let combined = ModalDriveFrame {
        localized_generalized_force_n: force,
        localized_boundary_impulse_n_s: impulse,
        ..distributed
    };
    let combined_with_zero_factor = with_cx(false, |cx| {
        model
            .synthesize_chunk(
                &initial,
                &[combined],
                ModalSpatialParticipation::PerFrameModeFactors(&[0.0]),
                cx,
            )
            .unwrap()
    });
    assert_eq!(combined_with_zero_factor, distributed_declared);

    let preparticipated_with_distributed = with_cx(false, |cx| {
        model
            .synthesize_chunk(
                &initial,
                &[distributed],
                ModalSpatialParticipation::PreparticipatedLocalizedDrive {
                    generalized_force_n: &[force.disc],
                    boundary_impulse_n_s: &[impulse.disc],
                },
                cx,
            )
            .unwrap()
    });
    let combined_declared = with_cx(false, |cx| {
        model
            .synthesize_chunk(
                &initial,
                &[combined],
                ModalSpatialParticipation::Declared,
                cx,
            )
            .unwrap()
    });
    assert_eq!(preparticipated_with_distributed, combined_declared);
}

#[test]
fn g0_representative_presets_are_distinct_complete_and_honestly_labeled() {
    let tungsten = representative_modal_preset(RepresentativeDiscMaterial::Tungsten);
    let stainless = representative_modal_preset(RepresentativeDiscMaterial::StainlessSteel);
    assert_ne!(tungsten.identity(), stainless.identity());
    assert_eq!(
        tungsten.authority(),
        ModalPresetAuthority::RepresentativeUncalibrated
    );
    assert!(
        tungsten
            .disclosure()
            .contains("not measured eigenfrequencies")
    );
    for preset in [&tungsten, &stainless] {
        assert!(preset.modes().iter().all(|mode| mode.damping_ratio > 0.0));
        for component in [
            SoundModalComponent::Disc,
            SoundModalComponent::GlassPlate,
            SoundModalComponent::BaseAssembly,
        ] {
            assert!(
                preset
                    .modes()
                    .iter()
                    .any(|mode| mode.component == component)
            );
        }
    }
    let tungsten_model = build(tungsten.modes().to_vec(), generous_budget(256));
    let stainless_model = build(stainless.modes().to_vec(), generous_budget(256));
    assert_ne!(tungsten_model.identity(), stainless_model.identity());
    let frames = {
        let mut values = vec![ModalDriveFrame::default(); 256];
        values[0].localized_boundary_impulse_n_s.disc = 1.0e-5;
        values
    };
    let render = |model: &ModalSynthesisModel| {
        let initial = with_cx(false, |cx| model.initial_checkpoint(cx).unwrap());
        with_cx(false, |cx| {
            model
                .synthesize_chunk(&initial, &frames, ModalSpatialParticipation::Declared, cx)
                .unwrap()
        })
    };
    assert_ne!(
        render(&tungsten_model).mixed_samples_fs,
        render(&stainless_model).mixed_samples_fs
    );
}

#[test]
fn g4_precancellation_publishes_no_chunk_or_successor() {
    let model = build(
        vec![mode(1, SoundModalComponent::Disc, 600.0, 0.02)],
        generous_budget(4),
    );
    let initial = with_cx(false, |cx| model.initial_checkpoint(cx).unwrap());
    let before = initial.clone();
    let result = with_cx(true, |cx| {
        model.synthesize_chunk(
            &initial,
            &[ModalDriveFrame::default()],
            ModalSpatialParticipation::Declared,
            cx,
        )
    });
    assert_eq!(result, Err(ModalSynthesisError::Cancelled));
    assert_eq!(initial, before);
}

#[test]
fn g5_split_resume_and_replay_are_bit_exact_and_wrong_model_refuses() {
    let mut modes = vec![
        mode(3, SoundModalComponent::GlassPlate, 1_100.0, 0.03),
        mode(1, SoundModalComponent::Disc, 430.0, 0.015),
    ];
    let model = build(modes.clone(), generous_budget(1_024));
    let frames: Vec<_> = (0..1_024)
        .map(|index| {
            let signed = (index % 13) as f64 - 6.0;
            let impulse = if index % 137 == 0 { 1.0e-5 } else { 0.0 };
            drive(
                ModalComponentValues {
                    disc: signed * 1.0e-3,
                    glass_plate: -signed * 5.0e-4,
                    base_assembly: 0.0,
                },
                ModalComponentValues {
                    disc: impulse,
                    glass_plate: -0.5 * impulse,
                    base_assembly: 0.0,
                },
            )
        })
        .collect();
    let initial = with_cx(false, |cx| model.initial_checkpoint(cx).unwrap());
    let one_shot = with_cx(false, |cx| {
        model
            .synthesize_chunk(&initial, &frames, ModalSpatialParticipation::Declared, cx)
            .unwrap()
    });
    let replay = with_cx(false, |cx| {
        model
            .synthesize_chunk(&initial, &frames, ModalSpatialParticipation::Declared, cx)
            .unwrap()
    });
    assert_eq!(one_shot, replay);

    let split = 377;
    let first = with_cx(false, |cx| {
        model
            .synthesize_chunk(
                &initial,
                &frames[..split],
                ModalSpatialParticipation::Declared,
                cx,
            )
            .unwrap()
    });
    let second = with_cx(false, |cx| {
        model
            .synthesize_chunk(
                &first.successor,
                &frames[split..],
                ModalSpatialParticipation::Declared,
                cx,
            )
            .unwrap()
    });
    let mut joined_samples = first.mixed_samples_fs;
    joined_samples.extend_from_slice(&second.mixed_samples_fs);
    let mut joined_stems = first.stem_frames;
    joined_stems.extend_from_slice(&second.stem_frames);
    let mut joined_energy = first.total_modal_energy_j;
    joined_energy.extend_from_slice(&second.total_modal_energy_j);
    assert_eq!(joined_samples, one_shot.mixed_samples_fs);
    assert_eq!(joined_stems, one_shot.stem_frames);
    assert_eq!(joined_energy, one_shot.total_modal_energy_j);
    assert_eq!(second.successor, one_shot.successor);

    modes[0].frequency_hz += 1.0;
    let other = build(modes, generous_budget(1_024));
    assert_eq!(
        with_cx(false, |cx| other.synthesize_chunk(
            &initial,
            &[ModalDriveFrame::default()],
            ModalSpatialParticipation::Declared,
            cx,
        )),
        Err(ModalSynthesisError::CheckpointIdentityMismatch)
    );
}
