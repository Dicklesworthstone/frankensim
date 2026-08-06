//! G0/G1/G3/G4/G5 integration evidence for the strict Euler cinematic WAV
//! boundary, deterministic dry mixer, meters, and typed sound manifest.

use core::f64::consts::TAU;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_blake3::{ContentHash, hash_domain};
use fs_euler_disc_e2e::{
    AudioArtifactBudget, AudioArtifactError, AudioArtifactRole, AudioDryMixSpec, AudioMasterSource,
    AudioSignalPath, ModalStemFrame, SoundWavArtifact, StemGainPan, StereoSample, WavMetadata,
    WavSampleEncoding, decode_stereo_wav, encode_stereo_wav, measure_audio, mix_dry_modal_stems,
    verify_wav_against_manifest,
};
use fs_evidence::{
    cinematic::{CinematicClock, CinematicClockDomain, CinematicDeliverableError, SoundAuthority},
    cinematic_config::{CinematicComponentRef, CinematicComponentRole},
    cinematic_sound::{
        ListenerFrame, ListenerPose, SOUND_MASTER_SAMPLE_RATE_HZ, SOUND_SYNTHESIS_SCHEMA_VERSION,
        SoundAmplitudeReference, SoundChannelLayout, SoundExcitationChannel,
        SoundExcitationControl, SoundModalComponent, SoundMode, SoundModeParticipation,
        SoundModelAssumption, SoundRoomResponse, SoundSynthesisConfig, SoundSynthesisInput,
        SoundTerminalPolicy, SoundTrajectoryDisposition,
    },
};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};

const PCM24_ZERO_FRAME_WAV: [u8; 44] = [
    b'R', b'I', b'F', b'F', 36, 0, 0, 0, b'W', b'A', b'V', b'E', b'f', b'm', b't', b' ', 16, 0, 0,
    0, 1, 0, 2, 0, 0x80, 0xbb, 0, 0, 0x00, 0x65, 0x04, 0, 6, 0, 24, 0, b'd', b'a', b't', b'a', 0,
    0, 0, 0,
];

const FLOAT32_ZERO_FRAME_WAV: [u8; 58] = [
    b'R', b'I', b'F', b'F', 50, 0, 0, 0, b'W', b'A', b'V', b'E', b'f', b'm', b't', b' ', 18, 0, 0,
    0, 3, 0, 2, 0, 0x80, 0xbb, 0, 0, 0x00, 0xdc, 0x05, 0, 8, 0, 32, 0, 0, 0, b'f', b'a', b'c',
    b't', 4, 0, 0, 0, 0, 0, 0, 0, b'd', b'a', b't', b'a', 0, 0, 0, 0,
];

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
                seed: 0x4155_4449_4f41_5254,
                kernel_id: 0x5741_565f,
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
    hash_domain("org.frankensim.test.audio-artifact.v1", label.as_bytes())
}

fn component(role: CinematicComponentRole, label: &str) -> CinematicComponentRef {
    CinematicComponentRef::try_new(role, identity(label), 1).unwrap()
}

fn sound_configuration(video_frames: i64, headroom_db: f64) -> SoundSynthesisConfig {
    let video_clock =
        CinematicClock::try_new(CinematicClockDomain::Video, 24, 1, 0, video_frames).unwrap();
    let audio_clock = CinematicClock::try_new(
        CinematicClockDomain::Audio,
        SOUND_MASTER_SAMPLE_RATE_HZ,
        1,
        0,
        video_frames * 2_000,
    )
    .unwrap();
    let mode = SoundMode {
        mode_id: 1,
        component: SoundModalComponent::Disc,
        frequency_hz: 440.0,
        damping_ratio: 0.02,
        modal_mass_kg: 0.2,
        source_participation: SoundModeParticipation {
            disc: 1.0,
            glass_plate: 0.0,
            base_assembly: 0.0,
        },
        radiation_gain_fs_s_per_m: 0.1,
        material_identity: identity("material"),
        base_identity: identity("base"),
    };
    SoundSynthesisConfig::try_admit(SoundSynthesisInput {
        schema_version: SOUND_SYNTHESIS_SCHEMA_VERSION,
        authority: SoundAuthority::PhysicallyInformed,
        trajectory: component(CinematicComponentRole::Trajectory, "trajectory"),
        excitation: component(CinematicComponentRole::AudioExcitation, "excitation"),
        sound_model: component(CinematicComponentRole::SoundModel, "sound-model"),
        microphone: component(CinematicComponentRole::Microphone, "microphone"),
        room: component(CinematicComponentRole::Room, "room"),
        timeline: component(CinematicComponentRole::Timeline, "timeline"),
        video_clock,
        audio_clock,
        channel_layout: SoundChannelLayout::Stereo,
        listener: ListenerPose {
            frame: ListenerFrame::AnimatedCamera,
            position_m: [0.0, 0.0, 0.0],
            forward: [0.0, 0.0, -1.0],
            up: [0.0, 1.0, 0.0],
        },
        excitation_controls: vec![SoundExcitationControl {
            channel: SoundExcitationChannel::ContactNormalForce,
            target_component: SoundModalComponent::Disc,
            source_scale: 1.0,
        }],
        modes: vec![mode],
        room_response: SoundRoomResponse::Dry,
        amplitude_reference: SoundAmplitudeReference::DigitalFullScale { headroom_db },
        trajectory_disposition: SoundTrajectoryDisposition::HorizonCensored,
        terminal_policy: SoundTerminalPolicy::FadeAtLastAccepted {
            fade_sample_frames: 240,
        },
        resampler_identity: identity("resampler"),
        resampler_version: 1,
        filter_identity: identity("filter"),
        filter_version: 1,
        assumptions: vec![
            SoundModelAssumption::LinearModalSuperposition,
            SoundModelAssumption::TimeInvariantDamping,
            SoundModelAssumption::DeclaredExcitationCompleteness,
            SoundModelAssumption::DeclaredRoomResponse,
        ],
        calibration: None,
    })
    .unwrap()
}

fn read_u32(bytes: &[u8], offset: usize) -> u32 {
    u32::from_le_bytes(bytes[offset..offset + 4].try_into().unwrap())
}

fn write_riff_size(bytes: &mut [u8]) {
    let riff_size = u32::try_from(bytes.len() - 8).unwrap();
    bytes[4..8].copy_from_slice(&riff_size.to_le_bytes());
}

fn assert_close(actual: f64, expected: f64, tolerance: f64, label: &str) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{label}: actual={actual:.16e} expected={expected:.16e} tolerance={tolerance:.3e}",
    );
}

#[test]
fn g1_riff_known_answers_cover_zero_and_one_frame_for_both_encodings() {
    with_cx(false, |cx| {
        let metadata = WavMetadata::default();

        let (pcm_zero, pcm_zero_receipt) = encode_stereo_wav(
            &[],
            SOUND_MASTER_SAMPLE_RATE_HZ,
            WavSampleEncoding::Pcm24,
            &metadata,
            AudioArtifactBudget::DEFAULT,
            cx,
        )
        .unwrap();
        assert_eq!(pcm_zero, PCM24_ZERO_FRAME_WAV);
        assert_eq!(pcm_zero_receipt.byte_len(), 44);
        assert_eq!(pcm_zero_receipt.sample_frame_count(), 0);
        assert!(
            decode_stereo_wav(&pcm_zero, AudioArtifactBudget::DEFAULT, cx)
                .unwrap()
                .samples
                .is_empty()
        );

        let (float_zero, float_zero_receipt) = encode_stereo_wav(
            &[],
            SOUND_MASTER_SAMPLE_RATE_HZ,
            WavSampleEncoding::Float32,
            &metadata,
            AudioArtifactBudget::DEFAULT,
            cx,
        )
        .unwrap();
        assert_eq!(float_zero, FLOAT32_ZERO_FRAME_WAV);
        assert_eq!(float_zero_receipt.byte_len(), 58);
        assert_eq!(float_zero_receipt.sample_frame_count(), 0);

        let pcm_frame = [StereoSample {
            left_fs: -1.0,
            right_fs: 1.0,
        }];
        let (pcm_one, pcm_one_receipt) = encode_stereo_wav(
            &pcm_frame,
            SOUND_MASTER_SAMPLE_RATE_HZ,
            WavSampleEncoding::Pcm24,
            &metadata,
            AudioArtifactBudget::DEFAULT,
            cx,
        )
        .unwrap();
        let mut expected_pcm_one = PCM24_ZERO_FRAME_WAV.to_vec();
        expected_pcm_one[4..8].copy_from_slice(&42_u32.to_le_bytes());
        expected_pcm_one[40..44].copy_from_slice(&6_u32.to_le_bytes());
        expected_pcm_one.extend_from_slice(&[0x00, 0x00, 0x80, 0xff, 0xff, 0x7f]);
        assert_eq!(pcm_one, expected_pcm_one);
        assert_eq!(pcm_one_receipt.byte_len(), 50);

        let float_frame = [StereoSample {
            left_fs: 0.5,
            right_fs: -0.25,
        }];
        let (float_one, float_one_receipt) = encode_stereo_wav(
            &float_frame,
            SOUND_MASTER_SAMPLE_RATE_HZ,
            WavSampleEncoding::Float32,
            &metadata,
            AudioArtifactBudget::DEFAULT,
            cx,
        )
        .unwrap();
        let mut expected_float_one = FLOAT32_ZERO_FRAME_WAV.to_vec();
        expected_float_one[4..8].copy_from_slice(&58_u32.to_le_bytes());
        expected_float_one[46..50].copy_from_slice(&1_u32.to_le_bytes());
        expected_float_one[54..58].copy_from_slice(&8_u32.to_le_bytes());
        expected_float_one.extend_from_slice(&0.5_f32.to_bits().to_le_bytes());
        expected_float_one.extend_from_slice(&(-0.25_f32).to_bits().to_le_bytes());
        assert_eq!(float_one, expected_float_one);
        assert_eq!(float_one_receipt.byte_len(), 66);
    });
}

#[test]
fn g0_pcm24_extrema_sign_extension_and_ties_are_exact() {
    with_cx(false, |cx| {
        let lsb = 1.0 / 8_388_608.0;
        let frames = [
            StereoSample {
                left_fs: -1.0,
                right_fs: 1.0,
            },
            StereoSample {
                left_fs: lsb,
                right_fs: -lsb,
            },
            StereoSample {
                left_fs: 0.5 * lsb,
                right_fs: 1.5 * lsb,
            },
        ];
        let (bytes, _) = encode_stereo_wav(
            &frames,
            SOUND_MASTER_SAMPLE_RATE_HZ,
            WavSampleEncoding::Pcm24,
            &WavMetadata::default(),
            AudioArtifactBudget::DEFAULT,
            cx,
        )
        .unwrap();
        assert_eq!(
            &bytes[44..],
            &[
                0x00, 0x00, 0x80, 0xff, 0xff, 0x7f, 0x01, 0x00, 0x00, 0xff, 0xff, 0xff, 0x00, 0x00,
                0x00, 0x02, 0x00, 0x00,
            ],
            "packed PCM must be signed little-endian and round half to even",
        );
        let decoded = decode_stereo_wav(&bytes, AudioArtifactBudget::DEFAULT, cx).unwrap();
        assert_eq!(decoded.samples[0].left_fs, -1.0);
        assert_eq!(decoded.samples[0].right_fs, 1.0 - lsb);
        assert_eq!(decoded.samples[1], frames[1]);
        assert_eq!(decoded.samples[2].left_fs, 0.0);
        assert_eq!(decoded.samples[2].right_fs, 2.0 * lsb);
    });
}

#[test]
fn g0_metadata_layout_padding_and_validation_are_canonical() {
    with_cx(false, |cx| {
        let metadata = WavMetadata::try_new(Some("ab".to_owned())).unwrap();
        let (bytes, receipt) = encode_stereo_wav(
            &[],
            SOUND_MASTER_SAMPLE_RATE_HZ,
            WavSampleEncoding::Pcm24,
            &metadata,
            AudioArtifactBudget::DEFAULT,
            cx,
        )
        .unwrap();
        assert_eq!(bytes.len(), 68);
        assert_eq!(receipt.byte_len(), 68);
        assert_eq!(&bytes[36..40], b"LIST");
        assert_eq!(read_u32(&bytes, 40), 16);
        assert_eq!(&bytes[44..48], b"INFO");
        assert_eq!(&bytes[48..52], b"ICMT");
        assert_eq!(read_u32(&bytes, 52), 3);
        assert_eq!(&bytes[56..60], &[b'a', b'b', 0, 0]);
        assert_eq!(&bytes[60..64], b"data");
        let decoded = decode_stereo_wav(&bytes, AudioArtifactBudget::DEFAULT, cx).unwrap();
        assert_eq!(decoded.metadata, metadata);

        let mut nonzero_pad = bytes.clone();
        nonzero_pad[59] = 1;
        assert_eq!(
            decode_stereo_wav(&nonzero_pad, AudioArtifactBudget::DEFAULT, cx),
            Err(AudioArtifactError::MalformedWav("nonzero chunk pad byte")),
        );
        assert_eq!(
            WavMetadata::try_new(Some(String::new())),
            Err(AudioArtifactError::InvalidMetadata),
        );
        assert_eq!(
            WavMetadata::try_new(Some("embedded\0nul".to_owned())),
            Err(AudioArtifactError::InvalidMetadata),
        );
        assert_eq!(
            WavMetadata::try_new(Some("non-ASCII café".to_owned())),
            Err(AudioArtifactError::InvalidMetadata),
        );
        assert_eq!(
            WavMetadata::try_new(Some("tab\tcontrol".to_owned())),
            Err(AudioArtifactError::InvalidMetadata),
        );
        assert_eq!(
            WavMetadata::try_new(Some("line one\nline two".to_owned()))
                .unwrap()
                .comment(),
            Some("line one\nline two"),
        );
        assert_eq!(
            WavMetadata::try_new(Some("x".repeat(4_097))),
            Err(AudioArtifactError::InvalidMetadata),
        );
    });
}

#[test]
fn g3_reader_refuses_every_truncation_unknown_features_and_trailing_bytes() {
    with_cx(false, |cx| {
        let metadata = WavMetadata::try_new(Some("prefix-sentinel".to_owned())).unwrap();
        let (complete, _) = encode_stereo_wav(
            &[StereoSample {
                left_fs: 0.25,
                right_fs: -0.5,
            }],
            SOUND_MASTER_SAMPLE_RATE_HZ,
            WavSampleEncoding::Float32,
            &metadata,
            AudioArtifactBudget::DEFAULT,
            cx,
        )
        .unwrap();
        for prefix_len in 0..complete.len() {
            assert!(
                decode_stereo_wav(&complete[..prefix_len], AudioArtifactBudget::DEFAULT, cx)
                    .is_err(),
                "truncated prefix {prefix_len}/{} was accepted",
                complete.len(),
            );
        }
        assert!(decode_stereo_wav(&complete, AudioArtifactBudget::DEFAULT, cx).is_ok());

        let mut non_riff = complete.clone();
        non_riff[0..4].copy_from_slice(b"RF64");
        assert_eq!(
            decode_stereo_wav(&non_riff, AudioArtifactBudget::DEFAULT, cx),
            Err(AudioArtifactError::UnsupportedWav("non-RIFF container")),
        );
        let mut non_wave = complete.clone();
        non_wave[8..12].copy_from_slice(b"AVI ");
        assert_eq!(
            decode_stereo_wav(&non_wave, AudioArtifactBudget::DEFAULT, cx),
            Err(AudioArtifactError::UnsupportedWav("non-WAVE RIFF form")),
        );

        let mut leading_junk = complete.clone();
        leading_junk.splice(12..12, [b'J', b'U', b'N', b'K', 0, 0, 0, 0]);
        write_riff_size(&mut leading_junk);
        assert_eq!(
            decode_stereo_wav(&leading_junk, AudioArtifactBudget::DEFAULT, cx),
            Err(AudioArtifactError::UnsupportedWav(
                "noncanonical chunk before fmt"
            )),
        );

        let mut missing_fact = complete.clone();
        missing_fact.drain(38..50);
        write_riff_size(&mut missing_fact);
        assert_eq!(
            decode_stereo_wav(&missing_fact, AudioArtifactBudget::DEFAULT, cx),
            Err(AudioArtifactError::UnsupportedWav(
                "noncanonical chunk before float fact"
            )),
        );

        let mut reordered_fact = complete.clone();
        let fact: Vec<_> = reordered_fact.drain(38..50).collect();
        reordered_fact.extend_from_slice(&fact);
        write_riff_size(&mut reordered_fact);
        assert_eq!(
            decode_stereo_wav(&reordered_fact, AudioArtifactBudget::DEFAULT, cx),
            Err(AudioArtifactError::UnsupportedWav(
                "noncanonical chunk before float fact"
            )),
        );

        let mut extended_fact = complete.clone();
        extended_fact[42..46].copy_from_slice(&8_u32.to_le_bytes());
        extended_fact.splice(50..50, [0, 0, 0, 0]);
        write_riff_size(&mut extended_fact);
        assert_eq!(
            decode_stereo_wav(&extended_fact, AudioArtifactBudget::DEFAULT, cx),
            Err(AudioArtifactError::UnsupportedWav(
                "extended float fact chunk"
            )),
        );

        let mut unsupported_tag = complete.clone();
        unsupported_tag[20..22].copy_from_slice(&2_u16.to_le_bytes());
        assert_eq!(
            decode_stereo_wav(&unsupported_tag, AudioArtifactBudget::DEFAULT, cx),
            Err(AudioArtifactError::UnsupportedWav("audio format tag")),
        );
        let mut mono = complete.clone();
        mono[22..24].copy_from_slice(&1_u16.to_le_bytes());
        assert_eq!(
            decode_stereo_wav(&mono, AudioArtifactBudget::DEFAULT, cx),
            Err(AudioArtifactError::UnsupportedWav(
                "non-stereo channel layout"
            )),
        );

        let data_offset = complete
            .windows(4)
            .position(|window| window == b"data")
            .unwrap();
        let mut unknown_chunk = complete.clone();
        unknown_chunk[data_offset..data_offset + 4].copy_from_slice(b"JUNK");
        assert_eq!(
            decode_stereo_wav(&unknown_chunk, AudioArtifactBudget::DEFAULT, cx),
            Err(AudioArtifactError::UnsupportedWav(
                "noncanonical or unknown chunk"
            )),
        );

        let mut truncated_tail = complete.clone();
        truncated_tail.push(0);
        write_riff_size(&mut truncated_tail);
        assert_eq!(
            decode_stereo_wav(&truncated_tail, AudioArtifactBudget::DEFAULT, cx),
            Err(AudioArtifactError::MalformedWav("truncated chunk header")),
        );

        let mut trailing_chunk = complete.clone();
        trailing_chunk.extend_from_slice(b"JUNK");
        trailing_chunk.extend_from_slice(&0_u32.to_le_bytes());
        write_riff_size(&mut trailing_chunk);
        assert_eq!(
            decode_stereo_wav(&trailing_chunk, AudioArtifactBudget::DEFAULT, cx),
            Err(AudioArtifactError::UnsupportedWav("chunk after data")),
        );

        let mut float_nan = complete.clone();
        let payload = data_offset + 8;
        float_nan[payload..payload + 4].copy_from_slice(&f32::NAN.to_bits().to_le_bytes());
        assert!(matches!(
            decode_stereo_wav(&float_nan, AudioArtifactBudget::DEFAULT, cx),
            Err(AudioArtifactError::NonFiniteSample {
                frame: 0,
                channel: "left"
            })
        ));
    });
}

#[test]
fn g5_codec_replay_is_byte_and_receipt_deterministic() {
    with_cx(false, |cx| {
        let samples = [
            StereoSample {
                left_fs: -0.75,
                right_fs: 0.125,
            },
            StereoSample {
                left_fs: 0.25,
                right_fs: -0.5,
            },
        ];
        let metadata = WavMetadata::try_new(Some("deterministic replay".to_owned())).unwrap();
        for encoding in [WavSampleEncoding::Pcm24, WavSampleEncoding::Float32] {
            let first = encode_stereo_wav(
                &samples,
                SOUND_MASTER_SAMPLE_RATE_HZ,
                encoding,
                &metadata,
                AudioArtifactBudget::DEFAULT,
                cx,
            )
            .unwrap();
            let second = encode_stereo_wav(
                &samples,
                SOUND_MASTER_SAMPLE_RATE_HZ,
                encoding,
                &metadata,
                AudioArtifactBudget::DEFAULT,
                cx,
            )
            .unwrap();
            assert_eq!(first, second, "{encoding:?} replay drifted");
            let decoded = decode_stereo_wav(&first.0, AudioArtifactBudget::DEFAULT, cx).unwrap();
            assert_eq!(decoded.receipt, first.1);
            assert_eq!(decoded.metadata, metadata);
        }
    });
}

#[test]
fn g0_g3_dry_mix_obeys_pan_gain_order_and_fails_closed() {
    with_cx(false, |cx| {
        let stems = [ModalStemFrame {
            disc_fs: 0.2,
            glass_plate_fs: 0.4,
            base_assembly_fs: 0.6,
        }];
        let spec = AudioDryMixSpec {
            disc: StemGainPan {
                gain_db: 0.0,
                pan: -1.0,
            },
            glass_plate: StemGainPan {
                gain_db: 0.0,
                pan: 1.0,
            },
            base_assembly: StemGainPan {
                gain_db: 0.0,
                pan: 0.0,
            },
            master_gain_db: -6.020_599_913_279_624,
        };
        let mixed = mix_dry_modal_stems(&stems, spec, AudioArtifactBudget::DEFAULT, cx).unwrap();
        let centre = 0.6 / 2.0_f64.sqrt();
        assert_close(mixed[0].left_fs, 0.5 * (0.2 + centre), 1.0e-15, "left mix");
        assert_close(
            mixed[0].right_fs,
            0.5 * (0.4 + centre),
            1.0e-15,
            "right mix",
        );
        assert_ne!(spec.identity(), AudioDryMixSpec::NEUTRAL.identity());

        let invalid_pan = AudioDryMixSpec {
            disc: StemGainPan {
                gain_db: 0.0,
                pan: 1.000_001,
            },
            ..AudioDryMixSpec::NEUTRAL
        };
        assert_eq!(
            mix_dry_modal_stems(&stems, invalid_pan, AudioArtifactBudget::DEFAULT, cx),
            Err(AudioArtifactError::InvalidMix("disc stem")),
        );

        let quiet_stem = [ModalStemFrame {
            disc_fs: 1.0e-6,
            glass_plate_fs: 0.0,
            base_assembly_fs: 0.0,
        }];
        let explicit_mastering = AudioDryMixSpec {
            master_gain_db: 120.0,
            ..AudioDryMixSpec::NEUTRAL
        };
        let mastered = mix_dry_modal_stems(
            &quiet_stem,
            explicit_mastering,
            AudioArtifactBudget::DEFAULT,
            cx,
        )
        .unwrap();
        assert_close(
            mastered[0].left_fs,
            1.0 / 2.0_f64.sqrt(),
            1.0e-12,
            "explicit positive master gain",
        );
        let excessive_mastering = AudioDryMixSpec {
            master_gain_db: 120.000_001,
            ..AudioDryMixSpec::NEUTRAL
        };
        assert_eq!(
            mix_dry_modal_stems(
                &quiet_stem,
                excessive_mastering,
                AudioArtifactBudget::DEFAULT,
                cx,
            ),
            Err(AudioArtifactError::InvalidMix("master gain")),
        );

        let nonfinite = [ModalStemFrame {
            disc_fs: f64::NAN,
            glass_plate_fs: 0.0,
            base_assembly_fs: 0.0,
        }];
        assert!(matches!(
            mix_dry_modal_stems(
                &nonfinite,
                AudioDryMixSpec::NEUTRAL,
                AudioArtifactBudget::DEFAULT,
                cx,
            ),
            Err(AudioArtifactError::NonFiniteSample {
                frame: 0,
                channel: "disc"
            })
        ));

        let clipping_stems = [ModalStemFrame {
            disc_fs: 1.0,
            glass_plate_fs: 1.0,
            base_assembly_fs: 1.0,
        }];
        let clipping = mix_dry_modal_stems(
            &clipping_stems,
            AudioDryMixSpec::NEUTRAL,
            AudioArtifactBudget::DEFAULT,
            cx,
        )
        .unwrap();
        assert!(clipping[0].left_fs > 1.0 && clipping[0].right_fs > 1.0);
        assert!(matches!(
            encode_stereo_wav(
                &clipping,
                SOUND_MASTER_SAMPLE_RATE_HZ,
                WavSampleEncoding::Float32,
                &WavMetadata::default(),
                AudioArtifactBudget::DEFAULT,
                cx,
            ),
            Err(AudioArtifactError::SampleOutOfRange { frame: 0, .. })
        ));
    });
}

#[test]
fn g0_g3_budgets_nonfinite_samples_and_precancellation_are_typed() {
    let sample = [StereoSample {
        left_fs: 0.25,
        right_fs: -0.25,
    }];
    with_cx(false, |cx| {
        let frame_budget = AudioArtifactBudget {
            maximum_sample_frames: 0,
            ..AudioArtifactBudget::DEFAULT
        };
        assert!(matches!(
            encode_stereo_wav(
                &sample,
                SOUND_MASTER_SAMPLE_RATE_HZ,
                WavSampleEncoding::Pcm24,
                &WavMetadata::default(),
                frame_budget,
                cx,
            ),
            Err(AudioArtifactError::BudgetExceeded {
                artifact: "WAV sample frames",
                requested: 1,
                limit: 0,
            })
        ));
        let work_budget = AudioArtifactBudget {
            maximum_work_items: 3,
            ..AudioArtifactBudget::DEFAULT
        };
        assert!(matches!(
            encode_stereo_wav(
                &sample,
                SOUND_MASTER_SAMPLE_RATE_HZ,
                WavSampleEncoding::Pcm24,
                &WavMetadata::default(),
                work_budget,
                cx,
            ),
            Err(AudioArtifactError::BudgetExceeded {
                artifact: "WAV encoding and hashing work",
                requested: 54,
                limit: 3,
            })
        ));
        let mix_work_budget = AudioArtifactBudget {
            maximum_work_items: 23,
            ..AudioArtifactBudget::DEFAULT
        };
        assert!(matches!(
            mix_dry_modal_stems(
                &[ModalStemFrame {
                    disc_fs: 0.1,
                    glass_plate_fs: 0.0,
                    base_assembly_fs: 0.0,
                }],
                AudioDryMixSpec::NEUTRAL,
                mix_work_budget,
                cx,
            ),
            Err(AudioArtifactError::BudgetExceeded {
                artifact: "dry mixing work",
                requested: 24,
                limit: 23,
            })
        ));
        let meter_work_budget = AudioArtifactBudget {
            maximum_work_items: 127,
            ..AudioArtifactBudget::DEFAULT
        };
        assert!(matches!(
            measure_audio(&sample, meter_work_budget, cx),
            Err(AudioArtifactError::BudgetExceeded {
                artifact: "audio metering work",
                requested: 128,
                limit: 127,
            })
        ));
        let byte_budget = AudioArtifactBudget {
            maximum_wav_bytes: 44,
            ..AudioArtifactBudget::DEFAULT
        };
        assert!(matches!(
            encode_stereo_wav(
                &sample,
                SOUND_MASTER_SAMPLE_RATE_HZ,
                WavSampleEncoding::Pcm24,
                &WavMetadata::default(),
                byte_budget,
                cx,
            ),
            Err(AudioArtifactError::BudgetExceeded {
                artifact: "complete WAV bytes",
                requested: 50,
                limit: 44,
            })
        ));
        assert_eq!(
            encode_stereo_wav(
                &[],
                44_100,
                WavSampleEncoding::Pcm24,
                &WavMetadata::default(),
                AudioArtifactBudget::DEFAULT,
                cx,
            ),
            Err(AudioArtifactError::InvalidSampleRate(44_100)),
        );
        let metadata = WavMetadata::try_new(Some("ab".to_owned())).unwrap();
        let metadata_budget = AudioArtifactBudget {
            maximum_metadata_bytes: 1,
            ..AudioArtifactBudget::DEFAULT
        };
        assert!(matches!(
            encode_stereo_wav(
                &[],
                SOUND_MASTER_SAMPLE_RATE_HZ,
                WavSampleEncoding::Pcm24,
                &metadata,
                metadata_budget,
                cx,
            ),
            Err(AudioArtifactError::BudgetExceeded {
                artifact: "WAV metadata bytes",
                requested: 2,
                limit: 1,
            })
        ));
        let invalid_budget = AudioArtifactBudget {
            maximum_wav_bytes: 43,
            ..AudioArtifactBudget::DEFAULT
        };
        assert_eq!(
            measure_audio(&[], invalid_budget, cx),
            Err(AudioArtifactError::InvalidBudget("maximum_wav_bytes")),
        );
        assert!(matches!(
            encode_stereo_wav(
                &[StereoSample {
                    left_fs: f64::INFINITY,
                    right_fs: 0.0,
                }],
                SOUND_MASTER_SAMPLE_RATE_HZ,
                WavSampleEncoding::Float32,
                &WavMetadata::default(),
                AudioArtifactBudget::DEFAULT,
                cx,
            ),
            Err(AudioArtifactError::NonFiniteSample {
                frame: 0,
                channel: "left"
            })
        ));
    });

    with_cx(true, |cancelled_cx| {
        assert_eq!(
            encode_stereo_wav(
                &sample,
                SOUND_MASTER_SAMPLE_RATE_HZ,
                WavSampleEncoding::Pcm24,
                &WavMetadata::default(),
                AudioArtifactBudget::DEFAULT,
                cancelled_cx,
            ),
            Err(AudioArtifactError::Cancelled),
        );
        assert_eq!(
            measure_audio(&sample, AudioArtifactBudget::DEFAULT, cancelled_cx),
            Err(AudioArtifactError::Cancelled),
        );
        let stems = [ModalStemFrame {
            disc_fs: 0.1,
            glass_plate_fs: 0.0,
            base_assembly_fs: 0.0,
        }];
        assert_eq!(
            mix_dry_modal_stems(
                &stems,
                AudioDryMixSpec::NEUTRAL,
                AudioArtifactBudget::DEFAULT,
                cancelled_cx,
            ),
            Err(AudioArtifactError::Cancelled),
        );
    });
}

#[test]
fn g1_meters_cover_silence_short_vectors_and_a_known_stereo_tone() {
    with_cx(false, |cx| {
        let empty = measure_audio(&[], AudioArtifactBudget::DEFAULT, cx).unwrap();
        assert_eq!(empty.sample_peak_fs, 0.0);
        assert_eq!(empty.true_peak_estimate_fs, 0.0);
        assert_eq!(empty.stereo_rms_fs, 0.0);
        assert_eq!(empty.dc_left_fs, 0.0);
        assert_eq!(empty.dc_right_fs, 0.0);
        assert_eq!(empty.integrated_loudness_lufs, None);
        assert_eq!(empty.loudness_block_count, 0);

        let short = [StereoSample {
            left_fs: 0.25,
            right_fs: -0.25,
        }];
        let short_meters = measure_audio(&short, AudioArtifactBudget::DEFAULT, cx).unwrap();
        assert_eq!(short_meters.sample_peak_fs, 0.25);
        assert_eq!(short_meters.true_peak_estimate_fs, 0.25);
        assert_eq!(short_meters.stereo_rms_fs, 0.25);
        assert_eq!(short_meters.dc_left_fs, 0.25);
        assert_eq!(short_meters.dc_right_fs, -0.25);
        assert_eq!(short_meters.integrated_loudness_lufs, None);

        let mut plateau = vec![StereoSample::default(); 33];
        plateau[15].left_fs = 1.0;
        plateau[16].left_fs = 1.0;
        let plateau_meters = measure_audio(&plateau, AudioArtifactBudget::DEFAULT, cx).unwrap();
        assert_eq!(plateau_meters.sample_peak_fs, 1.0);
        assert_close(
            plateau_meters.true_peak_estimate_fs,
            1.264_683_493_325_069_2,
            1.0e-12,
            "Lanczos-8 plateau intersample peak",
        );

        let tone: Vec<_> = (0..SOUND_MASTER_SAMPLE_RATE_HZ)
            .map(|frame| {
                let phase =
                    TAU * 1_000.0 * f64::from(frame) / f64::from(SOUND_MASTER_SAMPLE_RATE_HZ);
                let value = 0.1 * phase.sin();
                StereoSample {
                    left_fs: value,
                    right_fs: value,
                }
            })
            .collect();
        let tone_meters = measure_audio(&tone, AudioArtifactBudget::DEFAULT, cx).unwrap();
        assert_close(tone_meters.sample_peak_fs, 0.1, 1.0e-14, "tone sample peak");
        assert_close(
            tone_meters.stereo_rms_fs,
            0.1 / 2.0_f64.sqrt(),
            1.0e-12,
            "tone RMS",
        );
        assert_close(tone_meters.dc_left_fs, 0.0, 1.0e-15, "tone left DC");
        assert_close(tone_meters.dc_right_fs, 0.0, 1.0e-15, "tone right DC");
        assert_eq!(tone_meters.loudness_block_count, 7);
        assert_eq!(tone_meters.absolute_gated_block_count, 7);
        assert_eq!(tone_meters.relative_gated_block_count, 7);
        let loudness = tone_meters.integrated_loudness_lufs.unwrap();
        assert_close(loudness, -20.0, 0.15, "1 kHz stereo programme loudness");
        assert!(
            tone_meters.true_peak_estimate_fs >= tone_meters.sample_peak_fs
                && tone_meters.true_peak_estimate_fs <= 0.101,
            "unexpected 4x true-peak estimate {}",
            tone_meters.true_peak_estimate_fs,
        );

        let gated_tone: Vec<_> = (0..8 * SOUND_MASTER_SAMPLE_RATE_HZ)
            .map(|frame| {
                let amplitude = if frame < 4 * SOUND_MASTER_SAMPLE_RATE_HZ {
                    0.1
                } else {
                    0.000_1
                };
                let phase = TAU * 997.0 * f64::from(frame) / f64::from(SOUND_MASTER_SAMPLE_RATE_HZ);
                StereoSample {
                    left_fs: amplitude * phase.sin(),
                    right_fs: 0.0,
                }
            })
            .collect();
        let gated_meters = measure_audio(&gated_tone, AudioArtifactBudget::DEFAULT, cx).unwrap();
        assert_eq!(gated_meters.loudness_block_count, 77);
        assert_eq!(gated_meters.absolute_gated_block_count, 41);
        assert_eq!(gated_meters.relative_gated_block_count, 40);
        assert_close(
            gated_meters.integrated_loudness_lufs.unwrap(),
            -23.176_3,
            1.0e-3,
            "997 Hz mixed-level mono programme loudness",
        );
    });
}

#[test]
fn g3_g5_sound_artifact_binds_eight_second_dry_master_and_spatialized_path() {
    with_cx(false, |cx| {
        const VIDEO_FRAMES: i64 = 8 * 24;
        const AUDIO_FRAMES: usize = 8 * SOUND_MASTER_SAMPLE_RATE_HZ as usize;
        let configuration = sound_configuration(VIDEO_FRAMES, 6.0);
        let dry: Vec<_> = (0..AUDIO_FRAMES)
            .map(|frame| ModalStemFrame {
                disc_fs: 0.05
                    * (TAU * 440.0 * frame as f64 / f64::from(SOUND_MASTER_SAMPLE_RATE_HZ)).sin(),
                glass_plate_fs: 0.0,
                base_assembly_fs: 0.0,
            })
            .collect();
        let metadata =
            WavMetadata::try_new(Some("Euler disc canonical dry master".to_owned())).unwrap();
        let artifact = SoundWavArtifact::try_build(
            &configuration,
            AudioMasterSource::DryModalStems {
                frames: &dry,
                mix: AudioDryMixSpec::NEUTRAL,
                source_synthesis: configuration.receipt(),
            },
            WavSampleEncoding::Pcm24,
            metadata.clone(),
            AudioArtifactBudget::DEFAULT,
            cx,
        )
        .unwrap();
        let decoded = artifact.verify(AudioArtifactBudget::DEFAULT, cx).unwrap();
        assert_eq!(decoded.samples.len(), AUDIO_FRAMES);
        assert_eq!(decoded.metadata, metadata);
        assert_eq!(
            artifact.manifest().wav().sample_frame_count(),
            AUDIO_FRAMES as u64
        );
        assert_eq!(
            artifact.manifest().wav().encoding(),
            WavSampleEncoding::Pcm24
        );
        assert_eq!(
            artifact.manifest().role(),
            AudioArtifactRole::QuantizedPcm24Derivative,
        );
        assert_eq!(
            artifact.manifest().authority(),
            SoundAuthority::PhysicallyInformed
        );
        assert_eq!(
            artifact.manifest().channel_layout().path(),
            AudioSignalPath::CanonicalDryStereo,
        );
        assert_eq!(
            artifact.manifest().mix_identity(),
            Some(AudioDryMixSpec::NEUTRAL.identity())
        );
        assert_eq!(artifact.manifest().synthesis(), configuration.receipt());
        let json = artifact.manifest().to_manifest_json();
        assert!(json.contains("\"sample_frames\":384000"));
        assert!(json.contains("\"audio_frames_per_video_frame\":2000"));
        assert!(json.contains("\"artifact_role\":\"quantized-pcm24-derivative\""));
        assert!(json.contains("\"calibrated_acoustic_prediction\":false"));

        let mut mutated = artifact.wav_bytes().to_vec();
        let last = mutated.len() - 1;
        mutated[last] ^= 1;
        assert!(matches!(
            verify_wav_against_manifest(
                artifact.manifest(),
                &mutated,
                AudioArtifactBudget::DEFAULT,
                cx,
            ),
            Err(AudioArtifactError::WavIdentityMismatch { .. })
        ));
        drop(decoded);
        drop(artifact);
        drop(dry);

        let short_configuration = sound_configuration(1, 6.0);
        let short_spatialized = vec![
            StereoSample {
                left_fs: 0.01,
                right_fs: -0.02,
            };
            2_000
        ];
        assert!(matches!(
            SoundWavArtifact::try_build(
                &short_configuration,
                AudioMasterSource::SpatializedStereo {
                    frames: &short_spatialized,
                    spatialization_identity: identity("spatializer"),
                    source_synthesis: short_configuration.receipt(),
                },
                WavSampleEncoding::Float32,
                WavMetadata::default(),
                AudioArtifactBudget::DEFAULT,
                cx,
            ),
            Err(AudioArtifactError::OutsideCinematicDeliverable(
                CinematicDeliverableError::FrameCountOutOfRange {
                    got: 1,
                    minimum: 192,
                    maximum: 288,
                }
            ))
        ));

        let spatialized = vec![
            StereoSample {
                left_fs: 0.01,
                right_fs: -0.02,
            };
            AUDIO_FRAMES
        ];
        let spatialized_artifact = SoundWavArtifact::try_build(
            &configuration,
            AudioMasterSource::SpatializedStereo {
                frames: &spatialized,
                spatialization_identity: identity("spatializer"),
                source_synthesis: configuration.receipt(),
            },
            WavSampleEncoding::Float32,
            WavMetadata::default(),
            AudioArtifactBudget::DEFAULT,
            cx,
        )
        .unwrap();
        assert!(matches!(
            spatialized_artifact.manifest().channel_layout().path(),
            AudioSignalPath::SpatializedStereo { spatialization_identity }
                if spatialization_identity == identity("spatializer")
        ));
        assert_eq!(spatialized_artifact.manifest().mix_identity(), None);
        assert_eq!(
            spatialized_artifact.manifest().role(),
            AudioArtifactRole::AuthoritativeFloat32Master,
        );
        assert!(
            spatialized_artifact
                .manifest()
                .to_manifest_json()
                .contains("\"artifact_role\":\"authoritative-float32-master\"")
        );
        let spatialized_decoded = spatialized_artifact
            .verify(AudioArtifactBudget::DEFAULT, cx)
            .unwrap();
        assert_eq!(spatialized_decoded.samples.len(), AUDIO_FRAMES);
        assert_eq!(
            spatialized_decoded.samples[0],
            StereoSample {
                left_fs: f64::from(0.01_f32),
                right_fs: f64::from(-0.02_f32),
            },
            "float32 storage semantics must be visible after verification",
        );

        let mut shifted_input = configuration.input().clone();
        shifted_input.video_clock =
            CinematicClock::try_new(CinematicClockDomain::Video, 24, 1, -48, 144).unwrap();
        shifted_input.audio_clock = CinematicClock::try_new(
            CinematicClockDomain::Audio,
            SOUND_MASTER_SAMPLE_RATE_HZ,
            1,
            -96_000,
            288_000,
        )
        .unwrap();
        let shifted_configuration = SoundSynthesisConfig::try_admit(shifted_input).unwrap();
        let shifted_artifact = SoundWavArtifact::try_build(
            &shifted_configuration,
            AudioMasterSource::SpatializedStereo {
                frames: &spatialized,
                spatialization_identity: identity("spatializer"),
                source_synthesis: shifted_configuration.receipt(),
            },
            WavSampleEncoding::Float32,
            WavMetadata::default(),
            AudioArtifactBudget::DEFAULT,
            cx,
        )
        .unwrap();
        assert_eq!(
            shifted_artifact.manifest().wav().sample_frame_count(),
            AUDIO_FRAMES as u64,
            "deliverable admission must use clock differences, not zero origins",
        );
        assert_eq!(shifted_artifact.manifest().video_ticks(), (-48, 144));
        assert_eq!(
            shifted_artifact.manifest().audio_ticks(),
            (-96_000, 288_000),
        );
        assert_eq!(
            shifted_artifact.manifest().audio_frames_per_video_frame(),
            2_000,
            "24 fps and 48 kHz must bind an exact integral master-clock ratio",
        );
        assert_eq!(
            shifted_artifact.manifest().synthesis(),
            shifted_configuration.receipt(),
        );
        drop(shifted_artifact);

        let mismatched_configuration = sound_configuration(VIDEO_FRAMES, 5.0);
        assert!(matches!(
            SoundWavArtifact::try_build(
                &configuration,
                AudioMasterSource::SpatializedStereo {
                    frames: &spatialized,
                    spatialization_identity: identity("spatializer"),
                    source_synthesis: mismatched_configuration.receipt(),
                },
                WavSampleEncoding::Float32,
                WavMetadata::default(),
                AudioArtifactBudget::DEFAULT,
                cx,
            ),
            Err(AudioArtifactError::SourceSynthesisMismatch { expected, actual })
                if expected == configuration.receipt()
                    && actual == mismatched_configuration.receipt()
        ));

        assert!(matches!(
            SoundWavArtifact::try_build(
                &configuration,
                AudioMasterSource::SpatializedStereo {
                    frames: &spatialized[..AUDIO_FRAMES - 1],
                    spatialization_identity: identity("spatializer"),
                    source_synthesis: configuration.receipt(),
                },
                WavSampleEncoding::Float32,
                WavMetadata::default(),
                AudioArtifactBudget::DEFAULT,
                cx,
            ),
            Err(AudioArtifactError::SampleCountMismatch {
                expected: 384_000,
                actual: 383_999,
            })
        ));

        drop(spatialized_decoded);
        drop(spatialized_artifact);

        let over_headroom = vec![
            StereoSample {
                left_fs: 0.51,
                right_fs: 0.51,
            };
            AUDIO_FRAMES
        ];
        assert!(matches!(
            SoundWavArtifact::try_build(
                &configuration,
                AudioMasterSource::SpatializedStereo {
                    frames: &over_headroom,
                    spatialization_identity: identity("spatializer"),
                    source_synthesis: configuration.receipt(),
                },
                WavSampleEncoding::Float32,
                WavMetadata::default(),
                AudioArtifactBudget::DEFAULT,
                cx,
            ),
            Err(AudioArtifactError::HeadroomExceeded {
                observed_peak_fs,
                allowed_peak_fs,
            }) if observed_peak_fs >= 0.51 && allowed_peak_fs < 0.51
        ));
        drop(over_headroom);

        // This impulse is exactly admissible before encoding. PCM24 rounds it
        // upward, so only the mandatory decoded-domain recheck can refuse it.
        let quantization_configuration = sound_configuration(VIDEO_FRAMES, 3.0);
        let allowed_peak = 0.707_945_784_384_137_9_f64;
        let mut quantization_edge = vec![StereoSample::default(); AUDIO_FRAMES];
        quantization_edge[AUDIO_FRAMES / 2] = StereoSample {
            left_fs: allowed_peak,
            right_fs: allowed_peak,
        };
        let pre_encode_meters =
            measure_audio(&quantization_edge, AudioArtifactBudget::DEFAULT, cx).unwrap();
        let pre_encode_peak = pre_encode_meters
            .sample_peak_fs
            .max(pre_encode_meters.true_peak_estimate_fs);
        assert!(
            pre_encode_peak <= allowed_peak,
            "source-domain peak {pre_encode_peak} must pass the exact {allowed_peak} threshold",
        );

        let (quantized_bytes, _) = encode_stereo_wav(
            &quantization_edge,
            SOUND_MASTER_SAMPLE_RATE_HZ,
            WavSampleEncoding::Pcm24,
            &WavMetadata::default(),
            AudioArtifactBudget::DEFAULT,
            cx,
        )
        .unwrap();
        let quantized =
            decode_stereo_wav(&quantized_bytes, AudioArtifactBudget::DEFAULT, cx).unwrap();
        let quantized_meters =
            measure_audio(&quantized.samples, AudioArtifactBudget::DEFAULT, cx).unwrap();
        let expected_quantized_peak = 5_938_680.0 / 8_388_608.0;
        assert_eq!(quantized_meters.sample_peak_fs, expected_quantized_peak);
        assert_eq!(
            quantized_meters.true_peak_estimate_fs,
            expected_quantized_peak,
        );
        assert!(expected_quantized_peak > allowed_peak);

        assert!(matches!(
            SoundWavArtifact::try_build(
                &quantization_configuration,
                AudioMasterSource::SpatializedStereo {
                    frames: &quantization_edge,
                    spatialization_identity: identity("spatializer"),
                    source_synthesis: quantization_configuration.receipt(),
                },
                WavSampleEncoding::Pcm24,
                WavMetadata::default(),
                AudioArtifactBudget::DEFAULT,
                cx,
            ),
            Err(AudioArtifactError::HeadroomExceeded {
                observed_peak_fs,
                allowed_peak_fs,
            }) if observed_peak_fs > allowed_peak_fs
                && (allowed_peak_fs - allowed_peak).abs() <= 1.0e-15
        ));
    });
}
