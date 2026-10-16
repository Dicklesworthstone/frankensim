//! Deterministic stereo WAV artifacts for the Euler-disc cinematic pipeline.
//!
//! The implementation intentionally covers one small, replayable RIFF/WAVE
//! subset: 48 kHz stereo PCM24 or IEEE-float32, an optional bounded `ICMT`
//! comment, deterministic dry-stem mixing, and a strict reader for exactly the
//! emitted layout.  It is not a general media library.  No limiter, automatic
//! gain, dither, compression, room response, or acoustic calibration is hidden
//! in this boundary.

use core::fmt;

use fs_blake3::{ContentHash, DomainHasher};
use fs_evidence::{
    cinematic::{CinematicDeliverableContract, CinematicDeliverableError, SoundAuthority},
    cinematic_sound::{
        SOUND_MASTER_SAMPLE_RATE_HZ, SoundAmplitudeReference, SoundChannelLayout,
        SoundSynthesisConfig, SoundSynthesisReceipt, SoundTrajectoryDisposition,
    },
};
use fs_exec::Cx;
use fs_math::det;

use crate::modal_synthesis::ModalStemFrame;

/// Version of the sound-artifact manifest and canonical mixing policy.
///
/// V2 separates synthesis, presentation, and final-artifact authority so an
/// artistic spatialization cannot inherit a physically-informed synthesis
/// label in the standalone manifest.
pub const AUDIO_ARTIFACT_SCHEMA_VERSION: u16 = 2;
/// Version of the strict RIFF/WAVE subset emitted by this module.
pub const EULER_WAV_CODEC_VERSION: u16 = 1;
/// Maximum ASCII bytes in the optional canonical `ICMT` value.
pub const MAX_WAV_COMMENT_BYTES: usize = 4_096;
/// Cancellation polling interval for sample-domain work.
pub const AUDIO_ARTIFACT_CANCELLATION_POLL_FRAMES: usize = 1_024;
/// Cancellation polling interval for byte-domain hashing.
pub const AUDIO_ARTIFACT_CANCELLATION_POLL_BYTES: usize = 64 * 1_024;
/// Four-times oversampling used by the declared intersample-peak estimate.
pub const AUDIO_TRUE_PEAK_OVERSAMPLE_FACTOR: usize = 4;
/// Largest explicit post-mix presentation gain admitted by the dry mixer.
///
/// Physically informed but uncalibrated modal velocities can sit many orders
/// of magnitude below digital full scale. This bound permits their explicit
/// normalization without changing the mechanics-to-force mapping. Independent
/// sample and intersample headroom gates still refuse any over-range artifact;
/// this is never an SPL or acoustic-pressure claim.
pub const MAX_AUDIO_MASTER_GAIN_DB: f64 = 180.0;

const WAV_IDENTITY_DOMAIN: &str = "org.frankensim.euler-cinematic.wav.v1";
const WAV_METADATA_IDENTITY_DOMAIN: &str = "org.frankensim.euler-cinematic.wav-metadata.v1";
const CHANNEL_RECEIPT_IDENTITY_DOMAIN: &str =
    "org.frankensim.euler-cinematic.audio-channel-receipt.v1";
const DRY_MIX_IDENTITY_DOMAIN: &str = "org.frankensim.euler-cinematic.dry-mix.v1";
const MASTER_SOURCE_IDENTITY_DOMAIN: &str = "org.frankensim.euler-cinematic.master-source.v1";
const AUDIO_MANIFEST_IDENTITY_DOMAIN: &str = "org.frankensim.euler-cinematic.audio-manifest.v1";
const RIFF_HEADER_BYTES: u64 = 12;
const FMT_CHUNK_BYTES: u64 = 24;
const FACT_CHUNK_BYTES: u64 = 12;
const DATA_CHUNK_HEADER_BYTES: u64 = 8;
const MAX_CANONICAL_WAV_BYTES: u64 = u32::MAX as u64 + 8;
const LOUDNESS_BLOCK_FRAMES: usize = 19_200;
const LOUDNESS_HOP_FRAMES: usize = 4_800;
const LOUDNESS_ABSOLUTE_GATE_LUFS: f64 = -70.0;
const LOUDNESS_RELATIVE_GATE_LU: f64 = -10.0;
const LOUDNESS_OFFSET_LUFS: f64 = -0.691;
const LN_10: f64 = core::f64::consts::LN_10;
const PEAK_INTERPOLATOR_RADIUS: i64 = 8;

/// Supported sample representation in the emitted WAV data chunk.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(u16)]
pub enum WavSampleEncoding {
    /// Signed packed 24-bit little-endian PCM.
    Pcm24 = 1,
    /// IEEE-754 binary32 little-endian samples (`WAVE_FORMAT_IEEE_FLOAT`).
    Float32 = 3,
}

impl WavSampleEncoding {
    /// Stable manifest spelling.
    #[must_use]
    pub const fn code(self) -> &'static str {
        match self {
            Self::Pcm24 => "pcm-s24le",
            Self::Float32 => "ieee-f32le",
        }
    }

    const fn bytes_per_sample(self) -> u16 {
        match self {
            Self::Pcm24 => 3,
            Self::Float32 => 4,
        }
    }

    const fn bits_per_sample(self) -> u16 {
        self.bytes_per_sample() * 8
    }
}

/// Authority role of one supported encoding in the frozen cinematic contract.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum AudioArtifactRole {
    /// The `euler-disc-v1` authoritative stereo float32 WAV master.
    AuthoritativeFloat32Master,
    /// A deterministic PCM24 derivative that cannot replace the float master.
    QuantizedPcm24Derivative,
}

impl AudioArtifactRole {
    /// Stable manifest spelling.
    #[must_use]
    pub const fn code(self) -> &'static str {
        match self {
            Self::AuthoritativeFloat32Master => "authoritative-float32-master",
            Self::QuantizedPcm24Derivative => "quantized-pcm24-derivative",
        }
    }

    const fn for_encoding(encoding: WavSampleEncoding) -> Self {
        match encoding {
            WavSampleEncoding::Pcm24 => Self::QuantizedPcm24Derivative,
            WavSampleEncoding::Float32 => Self::AuthoritativeFloat32Master,
        }
    }
}

/// One stereo sample frame in digital-full-scale coordinates.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct StereoSample {
    /// Left channel sample.
    pub left_fs: f64,
    /// Right channel sample.
    pub right_fs: f64,
}

/// Explicit origin of the stereo master.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum AudioSignalPath {
    /// Static, camera-independent pan/gain of the three dry modal stems.
    CanonicalDryStereo,
    /// Stereo samples supplied by a separately identified spatializer.
    SpatializedStereo {
        /// Nonzero identity of the spatialization transform and its inputs.
        spatialization_identity: ContentHash,
    },
}

impl AudioSignalPath {
    /// Stable manifest spelling.
    #[must_use]
    pub const fn code(self) -> &'static str {
        match self {
            Self::CanonicalDryStereo => "canonical-dry-stereo",
            Self::SpatializedStereo { .. } => "spatialized-stereo",
        }
    }

    const fn presentation_authority(self) -> AudioPresentationAuthority {
        match self {
            Self::CanonicalDryStereo => AudioPresentationAuthority::SynthesisBoundCanonicalDryMix,
            // An identity alone cannot prove the authority of an upstream
            // spatializer. The artifact boundary therefore assigns the
            // conservative presentation class admitted for arbitrary supplied
            // stereo samples.
            Self::SpatializedStereo { .. } => AudioPresentationAuthority::Artistic,
        }
    }
}

/// Authority of the channel-presentation transform applied after synthesis.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum AudioPresentationAuthority {
    /// The deterministic canonical dry mix introduces no independent acoustic
    /// or spatial claim; final authority remains bounded by synthesis.
    SynthesisBoundCanonicalDryMix,
    /// Spatial presentation is selected for communication or aesthetics and
    /// is not a calibrated acoustic prediction.
    Artistic,
}

impl AudioPresentationAuthority {
    /// Stable manifest spelling.
    #[must_use]
    pub const fn code(self) -> &'static str {
        match self {
            Self::SynthesisBoundCanonicalDryMix => "synthesis-bound-canonical-dry-mix",
            Self::Artistic => "artistic",
        }
    }
}

/// Immutable receipt that prevents dry and camera-relative signals from being
/// confused at the WAV boundary.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct AudioChannelLayoutReceipt {
    identity: ContentHash,
    layout: SoundChannelLayout,
    path: AudioSignalPath,
}

impl AudioChannelLayoutReceipt {
    fn try_new(path: AudioSignalPath) -> Result<Self, AudioArtifactError> {
        if let AudioSignalPath::SpatializedStereo {
            spatialization_identity,
        } = path
            && is_zero(spatialization_identity)
        {
            return Err(AudioArtifactError::InvalidIdentity(
                "spatialization identity",
            ));
        }
        let layout = SoundChannelLayout::Stereo;
        let mut hasher = DomainHasher::new(CHANNEL_RECEIPT_IDENTITY_DOMAIN);
        hasher.update(&[layout.channels()]);
        match path {
            AudioSignalPath::CanonicalDryStereo => hasher.update(&[1]),
            AudioSignalPath::SpatializedStereo {
                spatialization_identity,
            } => {
                hasher.update(&[2]);
                hasher.update(spatialization_identity.as_bytes());
            }
        }
        Ok(Self {
            identity: hasher.finalize(),
            layout,
            path,
        })
    }

    /// Complete layout/path identity.
    #[must_use]
    pub const fn identity(self) -> ContentHash {
        self.identity
    }

    /// The v1 master is always stereo.
    #[must_use]
    pub const fn layout(self) -> SoundChannelLayout {
        self.layout
    }

    /// Dry or separately spatialized provenance.
    #[must_use]
    pub const fn path(self) -> AudioSignalPath {
        self.path
    }
}

/// Gain and equal-power pan for one mono stem.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StemGainPan {
    /// Gain in decibels; no implicit normalization follows it.
    pub gain_db: f64,
    /// Pan in `[-1, 1]`, from hard left through centre to hard right.
    pub pan: f64,
}

/// Complete deterministic dry mix.  Stem order is fixed as disc, glass, base.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AudioDryMixSpec {
    /// Disc-body stem treatment.
    pub disc: StemGainPan,
    /// Glass-plate stem treatment.
    pub glass_plate: StemGainPan,
    /// Base-assembly stem treatment.
    pub base_assembly: StemGainPan,
    /// Final explicit gain after stem summation, in decibels.
    pub master_gain_db: f64,
}

impl AudioDryMixSpec {
    /// A neutral, centred mix.  The equal-power centre law applies `-3.01 dB`
    /// to each output channel while preserving summed stereo power.
    pub const NEUTRAL: Self = Self {
        disc: StemGainPan {
            gain_db: 0.0,
            pan: 0.0,
        },
        glass_plate: StemGainPan {
            gain_db: 0.0,
            pan: 0.0,
        },
        base_assembly: StemGainPan {
            gain_db: 0.0,
            pan: 0.0,
        },
        master_gain_db: 0.0,
    };

    fn validate(self) -> Result<(), AudioArtifactError> {
        for (field, stem) in [
            ("disc stem", self.disc),
            ("glass stem", self.glass_plate),
            ("base stem", self.base_assembly),
        ] {
            if !stem.gain_db.is_finite() || !(-120.0..=24.0).contains(&stem.gain_db) {
                return Err(AudioArtifactError::InvalidMix(field));
            }
            if !stem.pan.is_finite() || !(-1.0..=1.0).contains(&stem.pan) {
                return Err(AudioArtifactError::InvalidMix(field));
            }
        }
        if !self.master_gain_db.is_finite()
            || !(-MAX_AUDIO_MASTER_GAIN_DB..=MAX_AUDIO_MASTER_GAIN_DB)
                .contains(&self.master_gain_db)
        {
            return Err(AudioArtifactError::InvalidMix("master gain"));
        }
        Ok(())
    }

    /// Content identity of all gains, pans, ordering, and the pan law.
    #[must_use]
    pub fn identity(self) -> ContentHash {
        let mut hasher = DomainHasher::new(DRY_MIX_IDENTITY_DOMAIN);
        hasher.update(&AUDIO_ARTIFACT_SCHEMA_VERSION.to_le_bytes());
        for stem in [self.disc, self.glass_plate, self.base_assembly] {
            hasher.update(&stem.gain_db.to_bits().to_le_bytes());
            hasher.update(&stem.pan.to_bits().to_le_bytes());
        }
        hasher.update(&self.master_gain_db.to_bits().to_le_bytes());
        hasher.update(b"equal-power-sqrt-v1");
        hasher.finalize()
    }
}

/// Master samples supplied to the artifact boundary.
#[derive(Debug, Clone, Copy)]
pub enum AudioMasterSource<'a> {
    /// Camera-independent component stems from modal synthesis.
    DryModalStems {
        /// One frame of disc/glass/base output per master sample frame.
        frames: &'a [ModalStemFrame],
        /// Explicit deterministic mix.
        mix: AudioDryMixSpec,
        /// Receipt asserted by the caller for the synthesis that produced the
        /// stems. It must exactly match the artifact configuration.
        source_synthesis: SoundSynthesisReceipt,
    },
    /// Stereo samples produced by an optional upstream spatializer.
    SpatializedStereo {
        /// Already-spatialized stereo frames. Because this boundary receives
        /// only samples plus an identity, their presentation is conservatively
        /// classified as artistic even when synthesis is physically informed.
        frames: &'a [StereoSample],
        /// Identity of that transform and its declared listener/room inputs.
        spatialization_identity: ContentHash,
        /// Receipt asserted by the caller for the synthesis input consumed by
        /// the spatializer. It must exactly match the artifact configuration.
        source_synthesis: SoundSynthesisReceipt,
    },
}

/// Optional bounded ASCII metadata emitted as one canonical `LIST/INFO/ICMT` chunk.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct WavMetadata {
    comment: Option<String>,
}

impl WavMetadata {
    /// Construct metadata. V1 accepts printable ASCII plus line feed. RIFF
    /// readers otherwise disagree about unmarked text encodings, so non-ASCII,
    /// empty strings, NULs, and other controls are refused.
    pub fn try_new(comment: Option<String>) -> Result<Self, AudioArtifactError> {
        if let Some(value) = &comment
            && (value.is_empty()
                || value.len() > MAX_WAV_COMMENT_BYTES
                || !value
                    .bytes()
                    .all(|byte| byte == b'\n' || (0x20..=0x7e).contains(&byte)))
        {
            return Err(AudioArtifactError::InvalidMetadata);
        }
        Ok(Self { comment })
    }

    /// Optional comment without its canonical terminating NUL.
    #[must_use]
    pub fn comment(&self) -> Option<&str> {
        self.comment.as_deref()
    }

    fn identity(&self) -> ContentHash {
        let mut hasher = DomainHasher::new(WAV_METADATA_IDENTITY_DOMAIN);
        match &self.comment {
            None => hasher.update(&[0]),
            Some(comment) => {
                hasher.update(&[1]);
                hasher.update(&(comment.len() as u64).to_le_bytes());
                hasher.update(comment.as_bytes());
            }
        }
        hasher.finalize()
    }
}

/// Caller-controlled resource ceilings under the 32-bit RIFF hard limit.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AudioArtifactBudget {
    /// Maximum stereo sample frames.
    pub maximum_sample_frames: u64,
    /// Maximum complete WAV bytes.
    pub maximum_wav_bytes: u64,
    /// Maximum metadata ASCII bytes.
    pub maximum_metadata_bytes: usize,
    /// Maximum deterministic sample/byte work estimate for one public
    /// operation. The high-level builder preflights its complete multi-stage
    /// transaction.
    pub maximum_work_items: u64,
}

impl AudioArtifactBudget {
    /// Fifteen minutes at 48 kHz, with room for the declared peak/loudness pass.
    pub const DEFAULT: Self = Self {
        maximum_sample_frames: 15 * 60 * SOUND_MASTER_SAMPLE_RATE_HZ as u64,
        maximum_wav_bytes: 512 * 1_024 * 1_024,
        maximum_metadata_bytes: MAX_WAV_COMMENT_BYTES,
        maximum_work_items: 10_000_000_000,
    };
}

impl Default for AudioArtifactBudget {
    fn default() -> Self {
        Self::DEFAULT
    }
}

/// Exact and explicitly estimated level diagnostics over decoded WAV samples.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AudioMeters {
    /// Exact maximum absolute stored sample.
    pub sample_peak_fs: f64,
    /// Four-times Lanczos-8 windowed-sinc intersample estimate under
    /// half-sample-even boundary extension. This is not dBTP, the
    /// continuous-time supremum, or a BS.1770 true-peak certificate.
    pub true_peak_estimate_fs: f64,
    /// Population RMS across both channels.
    pub stereo_rms_fs: f64,
    /// Mean left-channel sample.
    pub dc_left_fs: f64,
    /// Mean right-channel sample.
    pub dc_right_fs: f64,
    /// K-weighted, absolute/relative-gated integrated programme loudness.
    /// `None` means silence, no gated blocks, or a programme shorter than 400 ms.
    pub integrated_loudness_lufs: Option<f64>,
    /// Complete 400 ms blocks before gating.
    pub loudness_block_count: u64,
    /// Blocks surviving the absolute `-70 LUFS` gate.
    pub absolute_gated_block_count: u64,
    /// Blocks surviving both absolute and relative gates.
    pub relative_gated_block_count: u64,
}

/// Receipt for the strict WAV bytes themselves.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct WavCodecReceipt {
    wav_identity: ContentHash,
    metadata_identity: ContentHash,
    byte_len: u64,
    sample_frame_count: u64,
    sample_rate_hz: u32,
    encoding: WavSampleEncoding,
}

impl WavCodecReceipt {
    /// Domain-separated identity of every emitted byte.
    #[must_use]
    pub const fn wav_identity(self) -> ContentHash {
        self.wav_identity
    }

    /// Identity of the optional metadata value.
    #[must_use]
    pub const fn metadata_identity(self) -> ContentHash {
        self.metadata_identity
    }

    /// Exact complete byte length.
    #[must_use]
    pub const fn byte_len(self) -> u64 {
        self.byte_len
    }

    /// Exact number of interleaved stereo frames.
    #[must_use]
    pub const fn sample_frame_count(self) -> u64 {
        self.sample_frame_count
    }

    /// Declared sample rate.
    #[must_use]
    pub const fn sample_rate_hz(self) -> u32 {
        self.sample_rate_hz
    }

    /// Stored sample representation.
    #[must_use]
    pub const fn encoding(self) -> WavSampleEncoding {
        self.encoding
    }
}

/// Strictly decoded WAV subset.
#[derive(Debug, Clone, PartialEq)]
pub struct DecodedStereoWav {
    /// Decoded and format-quantized frames.
    pub samples: Vec<StereoSample>,
    /// Optional canonical metadata.
    pub metadata: WavMetadata,
    /// Byte-level receipt.
    pub receipt: WavCodecReceipt,
}

/// Immutable artifact manifest binding source, synthesis, layout, bytes, and meters.
#[derive(Debug, Clone, PartialEq)]
pub struct AudioArtifactManifest {
    identity: ContentHash,
    synthesis: SoundSynthesisReceipt,
    authority: SoundAuthority,
    channel_layout: AudioChannelLayoutReceipt,
    source_signal_identity: ContentHash,
    mix_identity: Option<ContentHash>,
    wav: WavCodecReceipt,
    role: AudioArtifactRole,
    meters: AudioMeters,
    video_start_tick: i64,
    video_end_tick_exclusive: i64,
    audio_start_tick: i64,
    audio_end_tick_exclusive: i64,
    audio_frames_per_video_frame: u32,
    admitted_headroom_db: f64,
}

impl AudioArtifactManifest {
    /// Complete manifest identity.
    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }

    /// Frozen sound configuration receipt.
    #[must_use]
    pub const fn synthesis(&self) -> SoundSynthesisReceipt {
        self.synthesis
    }

    /// Authority of the mechanics-driven synthesis before channel presentation.
    #[must_use]
    pub const fn synthesis_authority(&self) -> SoundAuthority {
        self.authority
    }

    /// Compatibility accessor for [`Self::synthesis_authority`].
    ///
    /// This is not the final artifact authority when presentation is artistic.
    #[must_use]
    pub const fn authority(&self) -> SoundAuthority {
        self.synthesis_authority()
    }

    /// Authority of the dry-mix or spatial-presentation stage.
    #[must_use]
    pub const fn presentation_authority(&self) -> AudioPresentationAuthority {
        self.channel_layout.path.presentation_authority()
    }

    /// Conservative authority of the samples as finally presented.
    #[must_use]
    pub const fn artifact_authority(&self) -> SoundAuthority {
        match self.presentation_authority() {
            AudioPresentationAuthority::SynthesisBoundCanonicalDryMix => self.synthesis_authority(),
            AudioPresentationAuthority::Artistic => SoundAuthority::Artistic,
        }
    }

    /// Explicit dry/spatialized stereo receipt.
    #[must_use]
    pub const fn channel_layout(&self) -> AudioChannelLayoutReceipt {
        self.channel_layout
    }

    /// Identity of pre-encoding source frames and their path declaration.
    #[must_use]
    pub const fn source_signal_identity(&self) -> ContentHash {
        self.source_signal_identity
    }

    /// Dry-mix identity, absent for already-spatialized input.
    #[must_use]
    pub const fn mix_identity(&self) -> Option<ContentHash> {
        self.mix_identity
    }

    /// Byte-level WAV receipt.
    #[must_use]
    pub const fn wav(&self) -> WavCodecReceipt {
        self.wav
    }

    /// Whether these bytes are the frozen float master or a PCM derivative.
    #[must_use]
    pub const fn role(&self) -> AudioArtifactRole {
        self.role
    }

    /// Reproducible meters over decoded/stored samples.
    #[must_use]
    pub const fn meters(&self) -> AudioMeters {
        self.meters
    }

    /// Exact video-clock interval bound by the sound configuration.
    #[must_use]
    pub const fn video_ticks(&self) -> (i64, i64) {
        (self.video_start_tick, self.video_end_tick_exclusive)
    }

    /// Exact audio-clock interval bound by the sound configuration.
    #[must_use]
    pub const fn audio_ticks(&self) -> (i64, i64) {
        (self.audio_start_tick, self.audio_end_tick_exclusive)
    }

    /// Exact integral master-clock ratio.
    #[must_use]
    pub const fn audio_frames_per_video_frame(&self) -> u32 {
        self.audio_frames_per_video_frame
    }

    /// Configured peak headroom applied as a refusal threshold.
    #[must_use]
    pub const fn admitted_headroom_db(&self) -> f64 {
        self.admitted_headroom_db
    }

    /// Deterministic JSON view generated from authoritative typed fields.
    #[must_use]
    pub fn to_manifest_json(&self) -> String {
        let loudness = self
            .meters
            .integrated_loudness_lufs
            .map_or_else(|| "null".to_string(), |value| value.to_string());
        let mix_identity = self.mix_identity.map_or_else(
            || "null".to_string(),
            |identity| format!("\"{}\"", identity.to_hex()),
        );
        let spatialization_identity = match self.channel_layout.path {
            AudioSignalPath::CanonicalDryStereo => "null".to_string(),
            AudioSignalPath::SpatializedStereo {
                spatialization_identity,
            } => format!("\"{}\"", spatialization_identity.to_hex()),
        };
        format!(
            concat!(
                "{{\"schema_version\":{},\"codec_version\":{},",
                "\"identity\":\"{}\",\"synthesis_schema_version\":{},",
                "\"synthesis_configuration_identity\":\"{}\",",
                "\"source_trajectory_identity\":\"{}\",",
                "\"excitation_identity\":\"{}\",",
                "\"sound_model_identity\":\"{}\",",
                "\"timeline_identity\":\"{}\",",
                "\"synthesis_authority\":\"{}\",",
                "\"presentation_authority\":\"{}\",",
                "\"artifact_authority\":\"{}\",",
                "\"channel_receipt_identity\":\"{}\",",
                "\"signal_path\":\"{}\",",
                "\"spatialization_identity\":{},",
                "\"source_signal_identity\":\"{}\",",
                "\"mix_identity\":{},",
                "\"channel_layout\":\"stereo\",\"sample_rate_hz\":{},",
                "\"sample_frames\":{},\"encoding\":\"{}\",",
                "\"artifact_role\":\"{}\",",
                "\"wav_bytes\":{},\"wav_identity\":\"{}\",",
                "\"metadata_identity\":\"{}\",",
                "\"sample_peak_fs\":{},\"true_peak_estimate_fs\":{},",
                "\"peak_estimator\":\"lanczos8-windowed-sinc-4x-half-sample-even-v1\",",
                "\"stereo_rms_fs\":{},\"dc_left_fs\":{},\"dc_right_fs\":{},",
                "\"integrated_loudness_lufs\":{},",
                "\"loudness_block_count\":{},",
                "\"absolute_gated_block_count\":{},",
                "\"relative_gated_block_count\":{},",
                "\"video_start_tick\":{},\"video_end_tick_exclusive\":{},",
                "\"audio_start_tick\":{},\"audio_end_tick_exclusive\":{},",
                "\"audio_frames_per_video_frame\":{},",
                "\"admitted_headroom_db\":{},",
                "\"calibrated_acoustic_prediction\":{}}}"
            ),
            AUDIO_ARTIFACT_SCHEMA_VERSION,
            EULER_WAV_CODEC_VERSION,
            self.identity.to_hex(),
            self.synthesis.schema_version,
            self.synthesis.configuration_identity.to_hex(),
            self.synthesis.trajectory_identity.to_hex(),
            self.synthesis.excitation_identity.to_hex(),
            self.synthesis.sound_model_identity.to_hex(),
            self.synthesis.timeline_identity.to_hex(),
            self.synthesis_authority().code(),
            self.presentation_authority().code(),
            self.artifact_authority().code(),
            self.channel_layout.identity.to_hex(),
            self.channel_layout.path.code(),
            spatialization_identity,
            self.source_signal_identity.to_hex(),
            mix_identity,
            self.wav.sample_rate_hz,
            self.wav.sample_frame_count,
            self.wav.encoding.code(),
            self.role.code(),
            self.wav.byte_len,
            self.wav.wav_identity.to_hex(),
            self.wav.metadata_identity.to_hex(),
            self.meters.sample_peak_fs,
            self.meters.true_peak_estimate_fs,
            self.meters.stereo_rms_fs,
            self.meters.dc_left_fs,
            self.meters.dc_right_fs,
            loudness,
            self.meters.loudness_block_count,
            self.meters.absolute_gated_block_count,
            self.meters.relative_gated_block_count,
            self.video_start_tick,
            self.video_end_tick_exclusive,
            self.audio_start_tick,
            self.audio_end_tick_exclusive,
            self.audio_frames_per_video_frame,
            self.admitted_headroom_db,
            self.artifact_authority() == SoundAuthority::Calibrated,
        )
    }
}

/// Complete canonical WAV plus its typed manifest.
#[derive(Debug, Clone, PartialEq)]
pub struct SoundWavArtifact {
    wav_bytes: Vec<u8>,
    manifest: AudioArtifactManifest,
}

impl SoundWavArtifact {
    /// Build a complete artifact without hidden normalization or limiting.
    pub fn try_build(
        configuration: &SoundSynthesisConfig,
        source: AudioMasterSource<'_>,
        encoding: WavSampleEncoding,
        metadata: WavMetadata,
        budget: AudioArtifactBudget,
        cx: &Cx<'_>,
    ) -> Result<Self, AudioArtifactError> {
        validate_budget(budget)?;
        checkpoint(cx)?;
        if configuration.input().trajectory_disposition
            == SoundTrajectoryDisposition::NumericalRefusal
        {
            return Err(AudioArtifactError::WaveformForbiddenByConfiguration);
        }
        if configuration.input().channel_layout != SoundChannelLayout::Stereo {
            return Err(AudioArtifactError::InvalidChannelLayout);
        }
        let expected_frames = expected_audio_frames(configuration)?;
        if expected_frames > budget.maximum_sample_frames {
            return Err(AudioArtifactError::BudgetExceeded {
                artifact: "master sample frames",
                requested: expected_frames,
                limit: budget.maximum_sample_frames,
            });
        }
        // The complete transaction performs source hashing, mixing/copying,
        // pre-encode peak admission, encoding and hashing, decoding and
        // hashing, and final metering. This conservative aggregate preflight
        // happens before any frame-sized allocation or traversal.
        check_work(expected_frames, 320, budget, "complete audio artifact work")?;
        let video_frames = expected_video_frames(configuration)?;
        CinematicDeliverableContract::euler_disc_v1()
            .validate_timeline(video_frames, expected_frames)
            .map_err(AudioArtifactError::OutsideCinematicDeliverable)?;
        let (master, channel_layout, source_signal_identity, mix_identity) = match source {
            AudioMasterSource::DryModalStems {
                frames,
                mix,
                source_synthesis,
            } => {
                validate_source_synthesis(configuration, source_synthesis)?;
                if frames.len() as u64 != expected_frames {
                    return Err(AudioArtifactError::SampleCountMismatch {
                        expected: expected_frames,
                        actual: frames.len() as u64,
                    });
                }
                let source_identity = hash_dry_source(frames, mix, source_synthesis, cx)?;
                let master = mix_dry_modal_stems(frames, mix, budget, cx)?;
                (
                    master,
                    AudioChannelLayoutReceipt::try_new(AudioSignalPath::CanonicalDryStereo)?,
                    source_identity,
                    Some(mix.identity()),
                )
            }
            AudioMasterSource::SpatializedStereo {
                frames,
                spatialization_identity,
                source_synthesis,
            } => {
                validate_source_synthesis(configuration, source_synthesis)?;
                if frames.len() as u64 != expected_frames {
                    return Err(AudioArtifactError::SampleCountMismatch {
                        expected: expected_frames,
                        actual: frames.len() as u64,
                    });
                }
                let path = AudioSignalPath::SpatializedStereo {
                    spatialization_identity,
                };
                let layout = AudioChannelLayoutReceipt::try_new(path)?;
                let source_identity = hash_stereo_source(frames, path, source_synthesis, cx)?;
                let master = copy_and_validate_stereo(frames, budget, cx)?;
                (master, layout, source_identity, None)
            }
        };
        let admitted_headroom_db = amplitude_headroom_db(configuration.input().amplitude_reference);
        let allowed_peak_fs = det::exp(-admitted_headroom_db * LN_10 / 20.0);
        let (sample_peak_fs, true_peak_estimate_fs) = measure_peaks(&master, cx)?;
        let observed_peak_fs = sample_peak_fs.max(true_peak_estimate_fs);
        if observed_peak_fs > allowed_peak_fs {
            return Err(AudioArtifactError::HeadroomExceeded {
                observed_peak_fs,
                allowed_peak_fs,
            });
        }
        let (wav_bytes, wav_receipt) = encode_stereo_wav(
            &master,
            SOUND_MASTER_SAMPLE_RATE_HZ,
            encoding,
            &metadata,
            budget,
            cx,
        )?;
        drop(master);
        let decoded = decode_stereo_wav(&wav_bytes, budget, cx)?;
        if decoded.receipt != wav_receipt || decoded.metadata != metadata {
            return Err(AudioArtifactError::NonCanonicalWav);
        }
        let meters = measure_audio(&decoded.samples, budget, cx)?;
        let decoded_peak_fs = meters.sample_peak_fs.max(meters.true_peak_estimate_fs);
        if decoded_peak_fs > allowed_peak_fs {
            return Err(AudioArtifactError::HeadroomExceeded {
                observed_peak_fs: decoded_peak_fs,
                allowed_peak_fs,
            });
        }
        let manifest = build_manifest(
            configuration,
            channel_layout,
            source_signal_identity,
            mix_identity,
            wav_receipt,
            meters,
            admitted_headroom_db,
        )?;
        Ok(Self {
            wav_bytes,
            manifest,
        })
    }

    /// Exact canonical WAV bytes.
    #[must_use]
    pub fn wav_bytes(&self) -> &[u8] {
        &self.wav_bytes
    }

    /// Typed authoritative manifest.
    #[must_use]
    pub const fn manifest(&self) -> &AudioArtifactManifest {
        &self.manifest
    }

    /// Decode and independently re-check byte identity, sample count, metadata,
    /// meters, and immutable manifest identity.
    pub fn verify(
        &self,
        budget: AudioArtifactBudget,
        cx: &Cx<'_>,
    ) -> Result<DecodedStereoWav, AudioArtifactError> {
        verify_wav_against_manifest(&self.manifest, &self.wav_bytes, budget, cx)
    }
}

/// Precise fail-closed errors for mixing, codec, metering, and verification.
#[derive(Debug, Clone, PartialEq)]
pub enum AudioArtifactError {
    /// Execution scope requested cancellation before atomic publication.
    Cancelled,
    /// One explicit budget is zero or exceeds a hard ceiling.
    InvalidBudget(&'static str),
    /// A checked work, frame, byte, or metadata ceiling was exceeded.
    BudgetExceeded {
        /// Named resource.
        artifact: &'static str,
        /// Requested amount.
        requested: u64,
        /// Admitted maximum.
        limit: u64,
    },
    /// Fallible allocation refused the requested result.
    Capacity {
        /// Named allocation.
        artifact: &'static str,
        /// Requested elements or bytes.
        requested: u64,
    },
    /// The sound configuration does not admit the fixed stereo master.
    InvalidChannelLayout,
    /// The caller-declared source synthesis does not match the artifact
    /// configuration. This catches accidental cross-configuration relabeling;
    /// it is not an independent proof that the samples were synthesized by it.
    SourceSynthesisMismatch {
        /// Receipt expected by the artifact configuration.
        expected: SoundSynthesisReceipt,
        /// Receipt declared for the supplied source samples.
        actual: SoundSynthesisReceipt,
    },
    /// A required identity was all zero.
    InvalidIdentity(&'static str),
    /// One gain or pan violates the declared finite range.
    InvalidMix(&'static str),
    /// Input frame count disagrees with the exact sound clock.
    SampleCountMismatch {
        /// Exact expected frames.
        expected: u64,
        /// Supplied or decoded frames.
        actual: u64,
    },
    /// One input sample is NaN or infinite.
    NonFiniteSample {
        /// Zero-based frame index.
        frame: u64,
        /// Stable channel/stem name.
        channel: &'static str,
    },
    /// A finite sample lies outside digital full scale.
    SampleOutOfRange {
        /// Zero-based frame index.
        frame: u64,
        /// Stable channel name.
        channel: &'static str,
    },
    /// Sample or intersample peak violates configured headroom. Samples are not changed.
    HeadroomExceeded {
        /// Largest observed linear peak.
        observed_peak_fs: f64,
        /// Maximum admitted linear peak.
        allowed_peak_fs: f64,
    },
    /// Numerical-refusal configurations explicitly forbid a waveform.
    WaveformForbiddenByConfiguration,
    /// Exact frame/sample count lies outside the frozen 8--12 second master.
    OutsideCinematicDeliverable(CinematicDeliverableError),
    /// Optional metadata is empty, too large, non-ASCII, or contains a
    /// disallowed control byte.
    InvalidMetadata,
    /// WAV sample rate is not the frozen 48 kHz master.
    InvalidSampleRate(u32),
    /// Structural RIFF/WAVE corruption or unsupported ordering.
    MalformedWav(&'static str),
    /// A valid RIFF feature lies outside the intentionally narrow subset.
    UnsupportedWav(&'static str),
    /// Decoded bytes do not reproduce the canonical writer receipt.
    NonCanonicalWav,
    /// WAV identity differs from the manifest.
    WavIdentityMismatch {
        /// Manifest identity.
        expected: ContentHash,
        /// Identity of supplied bytes.
        actual: ContentHash,
    },
    /// Manifest fields no longer produce the stored manifest identity.
    ManifestIdentityMismatch,
    /// Recomputed decoded-sample meters differ from the manifest.
    MeterMismatch,
}

impl fmt::Display for AudioArtifactError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "{self:?}")
    }
}

impl std::error::Error for AudioArtifactError {}

/// Mix the three canonical modal stems into stereo with a frozen equal-power
/// pan law.  Output is transactional and no gain is inferred from the signal.
pub fn mix_dry_modal_stems(
    frames: &[ModalStemFrame],
    spec: AudioDryMixSpec,
    budget: AudioArtifactBudget,
    cx: &Cx<'_>,
) -> Result<Vec<StereoSample>, AudioArtifactError> {
    validate_budget(budget)?;
    spec.validate()?;
    let frame_count = frames.len() as u64;
    if frame_count > budget.maximum_sample_frames {
        return Err(AudioArtifactError::BudgetExceeded {
            artifact: "dry mix sample frames",
            requested: frame_count,
            limit: budget.maximum_sample_frames,
        });
    }
    check_work(frame_count, 24, budget, "dry mixing work")?;
    let disc = stem_coefficients(spec.disc);
    let glass = stem_coefficients(spec.glass_plate);
    let base = stem_coefficients(spec.base_assembly);
    let master = db_to_linear(spec.master_gain_db);
    let mut output = Vec::new();
    reserve_exact(&mut output, frames.len(), "dry stereo master")?;
    for (frame_index, frame) in frames.iter().copied().enumerate() {
        if frame_index % AUDIO_ARTIFACT_CANCELLATION_POLL_FRAMES == 0 {
            checkpoint(cx)?;
        }
        for (channel, value) in [
            ("disc", frame.disc_fs),
            ("glass-plate", frame.glass_plate_fs),
            ("base-assembly", frame.base_assembly_fs),
        ] {
            if !value.is_finite() {
                return Err(AudioArtifactError::NonFiniteSample {
                    frame: frame_index as u64,
                    channel,
                });
            }
        }
        let mut left = CompensatedSum::new();
        left.add(frame.disc_fs * disc.0);
        left.add(frame.glass_plate_fs * glass.0);
        left.add(frame.base_assembly_fs * base.0);
        let mut right = CompensatedSum::new();
        right.add(frame.disc_fs * disc.1);
        right.add(frame.glass_plate_fs * glass.1);
        right.add(frame.base_assembly_fs * base.1);
        let sample = StereoSample {
            left_fs: left.total() * master,
            right_fs: right.total() * master,
        };
        if !sample.left_fs.is_finite() || !sample.right_fs.is_finite() {
            return Err(AudioArtifactError::NonFiniteSample {
                frame: frame_index as u64,
                channel: "mixed stereo",
            });
        }
        output.push(sample);
    }
    checkpoint(cx)?;
    Ok(output)
}

/// Encode the canonical stereo WAV subset.  Zero-frame files are valid.  The
/// returned vector is published only after the complete file has been built.
pub fn encode_stereo_wav(
    samples: &[StereoSample],
    sample_rate_hz: u32,
    encoding: WavSampleEncoding,
    metadata: &WavMetadata,
    budget: AudioArtifactBudget,
    cx: &Cx<'_>,
) -> Result<(Vec<u8>, WavCodecReceipt), AudioArtifactError> {
    validate_budget(budget)?;
    if sample_rate_hz != SOUND_MASTER_SAMPLE_RATE_HZ {
        return Err(AudioArtifactError::InvalidSampleRate(sample_rate_hz));
    }
    if samples.len() as u64 > budget.maximum_sample_frames {
        return Err(AudioArtifactError::BudgetExceeded {
            artifact: "WAV sample frames",
            requested: samples.len() as u64,
            limit: budget.maximum_sample_frames,
        });
    }
    let comment_len = metadata.comment().map_or(0, str::len);
    if comment_len > budget.maximum_metadata_bytes {
        return Err(AudioArtifactError::BudgetExceeded {
            artifact: "WAV metadata bytes",
            requested: comment_len as u64,
            limit: budget.maximum_metadata_bytes as u64,
        });
    }
    checkpoint(cx)?;

    let block_align = encoding
        .bytes_per_sample()
        .checked_mul(2)
        .ok_or(AudioArtifactError::MalformedWav("block alignment overflow"))?;
    let data_bytes = (samples.len() as u64)
        .checked_mul(u64::from(block_align))
        .ok_or(AudioArtifactError::BudgetExceeded {
            artifact: "WAV data bytes",
            requested: u64::MAX,
            limit: budget.maximum_wav_bytes,
        })?;
    if data_bytes > u64::from(u32::MAX) {
        return Err(AudioArtifactError::BudgetExceeded {
            artifact: "RIFF data chunk bytes",
            requested: data_bytes,
            limit: u64::from(u32::MAX),
        });
    }
    let fmt_bytes = match encoding {
        WavSampleEncoding::Pcm24 => FMT_CHUNK_BYTES,
        // WAVEFORMATEX for non-PCM includes a zero cbSize field.
        WavSampleEncoding::Float32 => FMT_CHUNK_BYTES + 2,
    };
    let fact_bytes = if encoding == WavSampleEncoding::Float32 {
        FACT_CHUNK_BYTES
    } else {
        0
    };
    let metadata_bytes = metadata_chunk_total_bytes(metadata)?;
    let total_bytes = RIFF_HEADER_BYTES
        .checked_add(fmt_bytes)
        .and_then(|value| value.checked_add(fact_bytes))
        .and_then(|value| value.checked_add(metadata_bytes))
        .and_then(|value| value.checked_add(DATA_CHUNK_HEADER_BYTES))
        .and_then(|value| value.checked_add(data_bytes))
        .ok_or(AudioArtifactError::BudgetExceeded {
            artifact: "complete WAV bytes",
            requested: u64::MAX,
            limit: budget.maximum_wav_bytes,
        })?;
    if total_bytes > budget.maximum_wav_bytes || total_bytes > MAX_CANONICAL_WAV_BYTES {
        return Err(AudioArtifactError::BudgetExceeded {
            artifact: "complete WAV bytes",
            requested: total_bytes,
            limit: budget.maximum_wav_bytes.min(MAX_CANONICAL_WAV_BYTES),
        });
    }
    check_combined_work(
        &[(samples.len() as u64, 4)],
        total_bytes,
        budget,
        "WAV encoding and hashing work",
    )?;
    let riff_size =
        u32::try_from(total_bytes - 8).map_err(|_| AudioArtifactError::BudgetExceeded {
            artifact: "RIFF size field",
            requested: total_bytes - 8,
            limit: u64::from(u32::MAX),
        })?;
    let frame_count_u32 =
        u32::try_from(samples.len()).map_err(|_| AudioArtifactError::BudgetExceeded {
            artifact: "WAV fact sample frames",
            requested: samples.len() as u64,
            limit: u64::from(u32::MAX),
        })?;
    let byte_rate = sample_rate_hz
        .checked_mul(u32::from(block_align))
        .ok_or(AudioArtifactError::MalformedWav("byte rate overflow"))?;
    let capacity = usize::try_from(total_bytes).map_err(|_| AudioArtifactError::Capacity {
        artifact: "complete WAV bytes",
        requested: total_bytes,
    })?;
    let mut bytes = Vec::new();
    reserve_exact(&mut bytes, capacity, "complete WAV bytes")?;
    bytes.extend_from_slice(b"RIFF");
    bytes.extend_from_slice(&riff_size.to_le_bytes());
    bytes.extend_from_slice(b"WAVE");
    bytes.extend_from_slice(b"fmt ");
    match encoding {
        WavSampleEncoding::Pcm24 => bytes.extend_from_slice(&16_u32.to_le_bytes()),
        WavSampleEncoding::Float32 => bytes.extend_from_slice(&18_u32.to_le_bytes()),
    }
    bytes.extend_from_slice(&(encoding as u16).to_le_bytes());
    bytes.extend_from_slice(&2_u16.to_le_bytes());
    bytes.extend_from_slice(&sample_rate_hz.to_le_bytes());
    bytes.extend_from_slice(&byte_rate.to_le_bytes());
    bytes.extend_from_slice(&block_align.to_le_bytes());
    bytes.extend_from_slice(&encoding.bits_per_sample().to_le_bytes());
    if encoding == WavSampleEncoding::Float32 {
        bytes.extend_from_slice(&0_u16.to_le_bytes());
        bytes.extend_from_slice(b"fact");
        bytes.extend_from_slice(&4_u32.to_le_bytes());
        bytes.extend_from_slice(&frame_count_u32.to_le_bytes());
    }
    append_metadata_chunk(&mut bytes, metadata);
    bytes.extend_from_slice(b"data");
    bytes.extend_from_slice(&(data_bytes as u32).to_le_bytes());
    for (frame_index, sample) in samples.iter().copied().enumerate() {
        if frame_index % AUDIO_ARTIFACT_CANCELLATION_POLL_FRAMES == 0 {
            checkpoint(cx)?;
        }
        encode_sample(
            &mut bytes,
            sample.left_fs,
            encoding,
            frame_index as u64,
            "left",
        )?;
        encode_sample(
            &mut bytes,
            sample.right_fs,
            encoding,
            frame_index as u64,
            "right",
        )?;
    }
    checkpoint(cx)?;
    if bytes.len() != capacity {
        return Err(AudioArtifactError::NonCanonicalWav);
    }
    let wav_identity = hash_bytes_cancellable(
        WAV_IDENTITY_DOMAIN,
        &bytes,
        budget,
        "WAV byte hashing work",
        cx,
    )?;
    let receipt = WavCodecReceipt {
        wav_identity,
        metadata_identity: metadata.identity(),
        byte_len: total_bytes,
        sample_frame_count: samples.len() as u64,
        sample_rate_hz,
        encoding,
    };
    Ok((bytes, receipt))
}

/// Decode exactly the subset emitted by [`encode_stereo_wav`].  Valid RIFF
/// features outside that subset return `UnsupportedWav` instead of being guessed.
pub fn decode_stereo_wav(
    bytes: &[u8],
    budget: AudioArtifactBudget,
    cx: &Cx<'_>,
) -> Result<DecodedStereoWav, AudioArtifactError> {
    validate_budget(budget)?;
    checkpoint(cx)?;
    if bytes.len() as u64 > budget.maximum_wav_bytes {
        return Err(AudioArtifactError::BudgetExceeded {
            artifact: "decoded WAV bytes",
            requested: bytes.len() as u64,
            limit: budget.maximum_wav_bytes,
        });
    }
    check_combined_work(
        &[],
        bytes.len() as u64,
        budget,
        "WAV decoding and hashing work",
    )?;
    if bytes.len() < RIFF_HEADER_BYTES as usize {
        return Err(AudioArtifactError::MalformedWav("truncated RIFF header"));
    }
    if bytes.get(0..4) != Some(b"RIFF") {
        return Err(AudioArtifactError::UnsupportedWav("non-RIFF container"));
    }
    if bytes.get(8..12) != Some(b"WAVE") {
        return Err(AudioArtifactError::UnsupportedWav("non-WAVE RIFF form"));
    }
    let declared_riff_size = read_u32(bytes, 4)? as u64;
    if declared_riff_size.checked_add(8) != Some(bytes.len() as u64) {
        return Err(AudioArtifactError::MalformedWav(
            "RIFF size or trailing bytes",
        ));
    }
    let mut position = 12_usize;
    let (fmt_id, fmt) = take_chunk(bytes, &mut position)?;
    if &fmt_id != b"fmt " {
        return Err(AudioArtifactError::UnsupportedWav(
            "noncanonical chunk before fmt",
        ));
    }
    let format_tag = read_u16(fmt, 0)?;
    let encoding = match format_tag {
        1 => WavSampleEncoding::Pcm24,
        3 => WavSampleEncoding::Float32,
        _ => return Err(AudioArtifactError::UnsupportedWav("audio format tag")),
    };
    let expected_fmt_len = match encoding {
        WavSampleEncoding::Pcm24 => 16,
        WavSampleEncoding::Float32 => 18,
    };
    if fmt.len() != expected_fmt_len {
        return Err(AudioArtifactError::UnsupportedWav("fmt extension"));
    }
    if encoding == WavSampleEncoding::Float32 && read_u16(fmt, 16)? != 0 {
        return Err(AudioArtifactError::UnsupportedWav("nonempty fmt extension"));
    }
    let channels = read_u16(fmt, 2)?;
    let sample_rate_hz = read_u32(fmt, 4)?;
    let byte_rate = read_u32(fmt, 8)?;
    let block_align = read_u16(fmt, 12)?;
    let bits_per_sample = read_u16(fmt, 14)?;
    let expected_block_align = encoding.bytes_per_sample() * 2;
    if channels != 2 {
        return Err(AudioArtifactError::UnsupportedWav(
            "non-stereo channel layout",
        ));
    }
    if sample_rate_hz != SOUND_MASTER_SAMPLE_RATE_HZ {
        return Err(AudioArtifactError::InvalidSampleRate(sample_rate_hz));
    }
    if block_align != expected_block_align
        || bits_per_sample != encoding.bits_per_sample()
        || byte_rate != sample_rate_hz * u32::from(expected_block_align)
    {
        return Err(AudioArtifactError::MalformedWav(
            "inconsistent fmt rates or widths",
        ));
    }
    let fact_frame_count = if encoding == WavSampleEncoding::Float32 {
        let (fact_id, fact) = take_chunk(bytes, &mut position)?;
        if &fact_id != b"fact" {
            return Err(AudioArtifactError::UnsupportedWav(
                "noncanonical chunk before float fact",
            ));
        }
        if fact.len() < 4 {
            return Err(AudioArtifactError::MalformedWav("float fact chunk size"));
        }
        if fact.len() > 4 {
            return Err(AudioArtifactError::UnsupportedWav(
                "extended float fact chunk",
            ));
        }
        Some(u64::from(read_u32(fact, 0)?))
    } else {
        None
    };
    let (mut next_id, mut next_payload) = take_chunk(bytes, &mut position)?;
    let metadata = if &next_id == b"LIST" {
        let parsed = parse_metadata_chunk(next_payload, budget)?;
        (next_id, next_payload) = take_chunk(bytes, &mut position)?;
        parsed
    } else {
        WavMetadata::default()
    };
    if &next_id != b"data" {
        return Err(AudioArtifactError::UnsupportedWav(
            "noncanonical or unknown chunk",
        ));
    }
    if position != bytes.len() {
        // A complete following chunk is a valid RIFF feature outside the
        // canonical subset. `take_chunk` preserves MalformedWav for a
        // truncated header/payload or nonzero/missing padding.
        let mut trailing_position = position;
        let _ = take_chunk(bytes, &mut trailing_position)?;
        return Err(AudioArtifactError::UnsupportedWav("chunk after data"));
    }
    if next_payload.len() % usize::from(block_align) != 0 {
        return Err(AudioArtifactError::MalformedWav(
            "partial stereo sample frame",
        ));
    }
    let sample_frames = next_payload.len() / usize::from(block_align);
    if sample_frames as u64 > budget.maximum_sample_frames {
        return Err(AudioArtifactError::BudgetExceeded {
            artifact: "decoded WAV sample frames",
            requested: sample_frames as u64,
            limit: budget.maximum_sample_frames,
        });
    }
    if fact_frame_count.is_some_and(|count| count != sample_frames as u64) {
        return Err(AudioArtifactError::MalformedWav(
            "fact sample count mismatch",
        ));
    }
    check_combined_work(
        &[(sample_frames as u64, 4)],
        bytes.len() as u64,
        budget,
        "WAV decoding and hashing work",
    )?;
    let mut samples = Vec::new();
    reserve_exact(&mut samples, sample_frames, "decoded WAV samples")?;
    for frame_index in 0..sample_frames {
        if frame_index % AUDIO_ARTIFACT_CANCELLATION_POLL_FRAMES == 0 {
            checkpoint(cx)?;
        }
        let offset = frame_index * usize::from(block_align);
        let left = decode_sample(next_payload, offset, encoding, frame_index as u64, "left")?;
        let right = decode_sample(
            next_payload,
            offset + usize::from(encoding.bytes_per_sample()),
            encoding,
            frame_index as u64,
            "right",
        )?;
        samples.push(StereoSample {
            left_fs: left,
            right_fs: right,
        });
    }
    checkpoint(cx)?;
    let receipt = WavCodecReceipt {
        wav_identity: hash_bytes_cancellable(
            WAV_IDENTITY_DOMAIN,
            bytes,
            budget,
            "WAV byte hashing work",
            cx,
        )?,
        metadata_identity: metadata.identity(),
        byte_len: bytes.len() as u64,
        sample_frame_count: sample_frames as u64,
        sample_rate_hz,
        encoding,
    };
    Ok(DecodedStereoWav {
        samples,
        metadata,
        receipt,
    })
}

/// Recompute final sample-domain metrics. The K-weighted loudness diagnostic
/// follows the BS.1770 stereo coefficient/gating recipe at the frozen 48 kHz
/// rate, but is not an external standards conformance certificate.
pub fn measure_audio(
    samples: &[StereoSample],
    budget: AudioArtifactBudget,
    cx: &Cx<'_>,
) -> Result<AudioMeters, AudioArtifactError> {
    validate_budget(budget)?;
    if samples.len() as u64 > budget.maximum_sample_frames {
        return Err(AudioArtifactError::BudgetExceeded {
            artifact: "meter sample frames",
            requested: samples.len() as u64,
            limit: budget.maximum_sample_frames,
        });
    }
    check_work(samples.len() as u64, 128, budget, "audio metering work")?;
    let (sample_peak_fs, true_peak_estimate_fs) = measure_peaks(samples, cx)?;
    let mut square_sum = CompensatedSum::new();
    let mut left_sum = CompensatedSum::new();
    let mut right_sum = CompensatedSum::new();
    for (frame_index, sample) in samples.iter().copied().enumerate() {
        if frame_index % AUDIO_ARTIFACT_CANCELLATION_POLL_FRAMES == 0 {
            checkpoint(cx)?;
        }
        validate_stereo_sample(sample, frame_index as u64)?;
        square_sum.add(sample.left_fs * sample.left_fs);
        square_sum.add(sample.right_fs * sample.right_fs);
        left_sum.add(sample.left_fs);
        right_sum.add(sample.right_fs);
    }
    let count = samples.len() as f64;
    let (stereo_rms_fs, dc_left_fs, dc_right_fs) = if samples.is_empty() {
        (0.0, 0.0, 0.0)
    } else {
        (
            det::sqrt((square_sum.total() / (2.0 * count)).max(0.0)),
            left_sum.total() / count,
            right_sum.total() / count,
        )
    };
    let loudness = measure_programme_loudness(samples, cx)?;
    checkpoint(cx)?;
    Ok(AudioMeters {
        sample_peak_fs,
        true_peak_estimate_fs,
        stereo_rms_fs,
        dc_left_fs,
        dc_right_fs,
        integrated_loudness_lufs: loudness.integrated_lufs,
        loudness_block_count: loudness.total_blocks,
        absolute_gated_block_count: loudness.absolute_blocks,
        relative_gated_block_count: loudness.relative_blocks,
    })
}

/// Verify supplied WAV bytes against an independently retained typed manifest.
pub fn verify_wav_against_manifest(
    manifest: &AudioArtifactManifest,
    wav_bytes: &[u8],
    budget: AudioArtifactBudget,
    cx: &Cx<'_>,
) -> Result<DecodedStereoWav, AudioArtifactError> {
    validate_budget(budget)?;
    checkpoint(cx)?;
    if wav_bytes.len() as u64 > budget.maximum_wav_bytes {
        return Err(AudioArtifactError::BudgetExceeded {
            artifact: "verified WAV bytes",
            requested: wav_bytes.len() as u64,
            limit: budget.maximum_wav_bytes,
        });
    }
    if manifest.wav.sample_frame_count > budget.maximum_sample_frames {
        return Err(AudioArtifactError::BudgetExceeded {
            artifact: "verified WAV sample frames",
            requested: manifest.wav.sample_frame_count,
            limit: budget.maximum_sample_frames,
        });
    }
    // Verification hashes once before parsing and once while producing the
    // decoder receipt, then performs decoding and complete metering.
    check_combined_work(
        &[(manifest.wav.sample_frame_count, 132)],
        (wav_bytes.len() as u64).checked_mul(2).unwrap_or(u64::MAX),
        budget,
        "WAV verification work",
    )?;
    if manifest_identity(manifest) != manifest.identity {
        return Err(AudioArtifactError::ManifestIdentityMismatch);
    }
    let actual_identity = hash_bytes_cancellable(
        WAV_IDENTITY_DOMAIN,
        wav_bytes,
        budget,
        "WAV verification hashing work",
        cx,
    )?;
    if actual_identity != manifest.wav.wav_identity {
        return Err(AudioArtifactError::WavIdentityMismatch {
            expected: manifest.wav.wav_identity,
            actual: actual_identity,
        });
    }
    let decoded = decode_stereo_wav(wav_bytes, budget, cx)?;
    if decoded.receipt != manifest.wav {
        return Err(AudioArtifactError::NonCanonicalWav);
    }
    if manifest.role != AudioArtifactRole::for_encoding(decoded.receipt.encoding) {
        return Err(AudioArtifactError::NonCanonicalWav);
    }
    let meters = measure_audio(&decoded.samples, budget, cx)?;
    if meters != manifest.meters {
        return Err(AudioArtifactError::MeterMismatch);
    }
    let allowed_peak_fs = det::exp(-manifest.admitted_headroom_db * LN_10 / 20.0);
    let observed_peak_fs = meters.sample_peak_fs.max(meters.true_peak_estimate_fs);
    if observed_peak_fs > allowed_peak_fs {
        return Err(AudioArtifactError::HeadroomExceeded {
            observed_peak_fs,
            allowed_peak_fs,
        });
    }
    checkpoint(cx)?;
    Ok(decoded)
}

fn validate_budget(budget: AudioArtifactBudget) -> Result<(), AudioArtifactError> {
    if budget.maximum_wav_bytes < 44 {
        return Err(AudioArtifactError::InvalidBudget("maximum_wav_bytes"));
    }
    if budget.maximum_wav_bytes > MAX_CANONICAL_WAV_BYTES {
        return Err(AudioArtifactError::InvalidBudget(
            "maximum_wav_bytes hard ceiling",
        ));
    }
    if budget.maximum_metadata_bytes > MAX_WAV_COMMENT_BYTES {
        return Err(AudioArtifactError::InvalidBudget(
            "maximum_metadata_bytes hard ceiling",
        ));
    }
    if budget.maximum_work_items == 0 {
        return Err(AudioArtifactError::InvalidBudget("maximum_work_items"));
    }
    Ok(())
}

fn check_work(
    frames: u64,
    items_per_frame: u64,
    budget: AudioArtifactBudget,
    artifact: &'static str,
) -> Result<(), AudioArtifactError> {
    check_combined_work(&[(frames, items_per_frame)], 0, budget, artifact)
}

fn check_combined_work(
    counted_terms: &[(u64, u64)],
    additional_items: u64,
    budget: AudioArtifactBudget,
    artifact: &'static str,
) -> Result<(), AudioArtifactError> {
    let mut requested = additional_items;
    for &(count, items_per_count) in counted_terms {
        requested = requested
            .checked_add(count.checked_mul(items_per_count).unwrap_or(u64::MAX))
            .unwrap_or(u64::MAX);
    }
    if requested > budget.maximum_work_items {
        Err(AudioArtifactError::BudgetExceeded {
            artifact,
            requested,
            limit: budget.maximum_work_items,
        })
    } else {
        Ok(())
    }
}

fn reserve_exact<T>(
    values: &mut Vec<T>,
    requested: usize,
    artifact: &'static str,
) -> Result<(), AudioArtifactError> {
    values
        .try_reserve_exact(requested)
        .map_err(|_| AudioArtifactError::Capacity {
            artifact,
            requested: requested as u64,
        })
}

fn expected_audio_frames(configuration: &SoundSynthesisConfig) -> Result<u64, AudioArtifactError> {
    let clock = configuration.input().audio_clock;
    let difference = i128::from(clock.end_tick_exclusive()) - i128::from(clock.start_tick());
    if difference < 0 {
        return Err(AudioArtifactError::MalformedWav(
            "negative configured audio horizon",
        ));
    }
    u64::try_from(difference).map_err(|_| AudioArtifactError::BudgetExceeded {
        artifact: "configured audio horizon",
        requested: u64::MAX,
        limit: u64::MAX,
    })
}

fn expected_video_frames(configuration: &SoundSynthesisConfig) -> Result<u32, AudioArtifactError> {
    let contract = CinematicDeliverableContract::euler_disc_v1();
    let clock = configuration.input().video_clock;
    let difference = i128::from(clock.end_tick_exclusive()) - i128::from(clock.start_tick());
    if difference < 0 {
        return Err(AudioArtifactError::MalformedWav(
            "negative configured video horizon",
        ));
    }
    u32::try_from(difference).map_err(|_| {
        AudioArtifactError::OutsideCinematicDeliverable(
            CinematicDeliverableError::FrameCountOutOfRange {
                got: u32::MAX,
                minimum: contract.minimum_frame_count(),
                maximum: contract.maximum_frame_count(),
            },
        )
    })
}

fn amplitude_headroom_db(reference: SoundAmplitudeReference) -> f64 {
    match reference {
        SoundAmplitudeReference::DigitalFullScale { headroom_db }
        | SoundAmplitudeReference::CalibratedPressure { headroom_db, .. } => headroom_db,
    }
}

fn validate_source_synthesis(
    configuration: &SoundSynthesisConfig,
    source_synthesis: SoundSynthesisReceipt,
) -> Result<(), AudioArtifactError> {
    let expected = configuration.receipt();
    if source_synthesis != expected {
        return Err(AudioArtifactError::SourceSynthesisMismatch {
            expected,
            actual: source_synthesis,
        });
    }
    Ok(())
}

fn stem_coefficients(stem: StemGainPan) -> (f64, f64) {
    let gain = db_to_linear(stem.gain_db);
    let left = gain * det::sqrt((1.0 - stem.pan) * 0.5);
    let right = gain * det::sqrt((1.0 + stem.pan) * 0.5);
    (left, right)
}

fn db_to_linear(db: f64) -> f64 {
    det::exp(db * LN_10 / 20.0)
}

fn copy_and_validate_stereo(
    frames: &[StereoSample],
    budget: AudioArtifactBudget,
    cx: &Cx<'_>,
) -> Result<Vec<StereoSample>, AudioArtifactError> {
    check_work(frames.len() as u64, 2, budget, "spatialized copy work")?;
    let mut output = Vec::new();
    reserve_exact(&mut output, frames.len(), "spatialized stereo master")?;
    for (frame_index, sample) in frames.iter().copied().enumerate() {
        if frame_index % AUDIO_ARTIFACT_CANCELLATION_POLL_FRAMES == 0 {
            checkpoint(cx)?;
        }
        validate_stereo_sample(sample, frame_index as u64)?;
        output.push(sample);
    }
    checkpoint(cx)?;
    Ok(output)
}

fn validate_stereo_sample(sample: StereoSample, frame: u64) -> Result<(), AudioArtifactError> {
    for (channel, value) in [("left", sample.left_fs), ("right", sample.right_fs)] {
        if !value.is_finite() {
            return Err(AudioArtifactError::NonFiniteSample { frame, channel });
        }
    }
    Ok(())
}

fn hash_dry_source(
    frames: &[ModalStemFrame],
    mix: AudioDryMixSpec,
    source_synthesis: SoundSynthesisReceipt,
    cx: &Cx<'_>,
) -> Result<ContentHash, AudioArtifactError> {
    mix.validate()?;
    let mut hasher = DomainHasher::new(MASTER_SOURCE_IDENTITY_DOMAIN);
    hasher.update(&[1]);
    hash_synthesis_receipt(&mut hasher, source_synthesis);
    hasher.update(mix.identity().as_bytes());
    hasher.update(&(frames.len() as u64).to_le_bytes());
    for (frame_index, frame) in frames.iter().copied().enumerate() {
        if frame_index % AUDIO_ARTIFACT_CANCELLATION_POLL_FRAMES == 0 {
            checkpoint(cx)?;
        }
        for (channel, value) in [
            ("disc", frame.disc_fs),
            ("glass-plate", frame.glass_plate_fs),
            ("base-assembly", frame.base_assembly_fs),
        ] {
            if !value.is_finite() {
                return Err(AudioArtifactError::NonFiniteSample {
                    frame: frame_index as u64,
                    channel,
                });
            }
            hasher.update(&value.to_bits().to_le_bytes());
        }
    }
    checkpoint(cx)?;
    Ok(hasher.finalize())
}

fn hash_stereo_source(
    frames: &[StereoSample],
    path: AudioSignalPath,
    source_synthesis: SoundSynthesisReceipt,
    cx: &Cx<'_>,
) -> Result<ContentHash, AudioArtifactError> {
    let mut hasher = DomainHasher::new(MASTER_SOURCE_IDENTITY_DOMAIN);
    hasher.update(&[2]);
    hash_synthesis_receipt(&mut hasher, source_synthesis);
    if let AudioSignalPath::SpatializedStereo {
        spatialization_identity,
    } = path
    {
        if is_zero(spatialization_identity) {
            return Err(AudioArtifactError::InvalidIdentity(
                "spatialization identity",
            ));
        }
        hasher.update(spatialization_identity.as_bytes());
    } else {
        return Err(AudioArtifactError::InvalidChannelLayout);
    }
    hasher.update(&(frames.len() as u64).to_le_bytes());
    for (frame_index, sample) in frames.iter().copied().enumerate() {
        if frame_index % AUDIO_ARTIFACT_CANCELLATION_POLL_FRAMES == 0 {
            checkpoint(cx)?;
        }
        validate_stereo_sample(sample, frame_index as u64)?;
        hasher.update(&sample.left_fs.to_bits().to_le_bytes());
        hasher.update(&sample.right_fs.to_bits().to_le_bytes());
    }
    checkpoint(cx)?;
    Ok(hasher.finalize())
}

fn metadata_chunk_total_bytes(metadata: &WavMetadata) -> Result<u64, AudioArtifactError> {
    let Some(comment) = metadata.comment() else {
        return Ok(0);
    };
    let text_bytes = (comment.len() as u64)
        .checked_add(1)
        .ok_or(AudioArtifactError::InvalidMetadata)?;
    let inner_padded = text_bytes
        .checked_add(text_bytes & 1)
        .ok_or(AudioArtifactError::InvalidMetadata)?;
    // LIST header + INFO form tag + ICMT header + NUL-terminated text/pad.
    8_u64
        .checked_add(4)
        .and_then(|value| value.checked_add(8))
        .and_then(|value| value.checked_add(inner_padded))
        .ok_or(AudioArtifactError::InvalidMetadata)
}

fn append_metadata_chunk(bytes: &mut Vec<u8>, metadata: &WavMetadata) {
    let Some(comment) = metadata.comment() else {
        return;
    };
    let data_len = comment.len() + 1;
    let inner_total = 8 + data_len + (data_len & 1);
    let list_payload_len = 4 + inner_total;
    bytes.extend_from_slice(b"LIST");
    bytes.extend_from_slice(&(list_payload_len as u32).to_le_bytes());
    bytes.extend_from_slice(b"INFO");
    bytes.extend_from_slice(b"ICMT");
    bytes.extend_from_slice(&(data_len as u32).to_le_bytes());
    bytes.extend_from_slice(comment.as_bytes());
    bytes.push(0);
    if data_len & 1 == 1 {
        bytes.push(0);
    }
}

fn encode_sample(
    bytes: &mut Vec<u8>,
    value: f64,
    encoding: WavSampleEncoding,
    frame: u64,
    channel: &'static str,
) -> Result<(), AudioArtifactError> {
    if !value.is_finite() {
        return Err(AudioArtifactError::NonFiniteSample { frame, channel });
    }
    if !(-1.0..=1.0).contains(&value) {
        return Err(AudioArtifactError::SampleOutOfRange { frame, channel });
    }
    match encoding {
        WavSampleEncoding::Pcm24 => {
            let quantized = quantize_pcm24(value);
            let encoded = quantized.to_le_bytes();
            bytes.extend_from_slice(&encoded[..3]);
        }
        WavSampleEncoding::Float32 => {
            let converted = value as f32;
            if !converted.is_finite() {
                return Err(AudioArtifactError::NonFiniteSample { frame, channel });
            }
            bytes.extend_from_slice(&converted.to_bits().to_le_bytes());
        }
    }
    Ok(())
}

fn quantize_pcm24(value: f64) -> i32 {
    const SCALE: f64 = 8_388_608.0;
    const MINIMUM: i32 = -8_388_608;
    const MAXIMUM: i32 = 8_388_607;
    if value <= -1.0 {
        MINIMUM
    } else if value >= 1.0 {
        MAXIMUM
    } else {
        (value * SCALE)
            .round_ties_even()
            .clamp(f64::from(MINIMUM), f64::from(MAXIMUM)) as i32
    }
}

fn decode_sample(
    bytes: &[u8],
    offset: usize,
    encoding: WavSampleEncoding,
    frame: u64,
    channel: &'static str,
) -> Result<f64, AudioArtifactError> {
    match encoding {
        WavSampleEncoding::Pcm24 => {
            let sample = bytes
                .get(offset..offset + 3)
                .ok_or(AudioArtifactError::MalformedWav("truncated PCM24 sample"))?;
            let packed =
                u32::from(sample[0]) | (u32::from(sample[1]) << 8) | (u32::from(sample[2]) << 16);
            let signed = if packed & 0x0080_0000 != 0 {
                (packed | 0xff00_0000) as i32
            } else {
                packed as i32
            };
            Ok(f64::from(signed) / 8_388_608.0)
        }
        WavSampleEncoding::Float32 => {
            let sample = bytes
                .get(offset..offset + 4)
                .ok_or(AudioArtifactError::MalformedWav("truncated float32 sample"))?;
            let bits = u32::from_le_bytes([sample[0], sample[1], sample[2], sample[3]]);
            let value = f32::from_bits(bits);
            if !value.is_finite() {
                return Err(AudioArtifactError::NonFiniteSample { frame, channel });
            }
            if !(-1.0..=1.0).contains(&value) {
                return Err(AudioArtifactError::SampleOutOfRange { frame, channel });
            }
            Ok(f64::from(value))
        }
    }
}

fn take_chunk<'a>(
    bytes: &'a [u8],
    position: &mut usize,
) -> Result<([u8; 4], &'a [u8]), AudioArtifactError> {
    let header = bytes
        .get(*position..position.saturating_add(8))
        .ok_or(AudioArtifactError::MalformedWav("truncated chunk header"))?;
    let id = [header[0], header[1], header[2], header[3]];
    let size = u32::from_le_bytes([header[4], header[5], header[6], header[7]]) as usize;
    let payload_start = position
        .checked_add(8)
        .ok_or(AudioArtifactError::MalformedWav("chunk offset overflow"))?;
    let payload_end = payload_start
        .checked_add(size)
        .ok_or(AudioArtifactError::MalformedWav("chunk size overflow"))?;
    let payload = bytes
        .get(payload_start..payload_end)
        .ok_or(AudioArtifactError::MalformedWav("truncated chunk payload"))?;
    let padded_end = payload_end
        .checked_add(size & 1)
        .ok_or(AudioArtifactError::MalformedWav("chunk padding overflow"))?;
    if size & 1 == 1 && bytes.get(payload_end).copied() != Some(0) {
        return Err(AudioArtifactError::MalformedWav("nonzero chunk pad byte"));
    }
    if padded_end > bytes.len() {
        return Err(AudioArtifactError::MalformedWav("missing chunk pad byte"));
    }
    *position = padded_end;
    Ok((id, payload))
}

fn parse_metadata_chunk(
    payload: &[u8],
    budget: AudioArtifactBudget,
) -> Result<WavMetadata, AudioArtifactError> {
    if payload.get(0..4) != Some(b"INFO") {
        return Err(AudioArtifactError::UnsupportedWav("non-INFO LIST chunk"));
    }
    let mut position = 4_usize;
    let (id, comment_payload) = take_chunk(payload, &mut position)?;
    if &id != b"ICMT" || position != payload.len() {
        return Err(AudioArtifactError::UnsupportedWav(
            "noncanonical INFO metadata",
        ));
    }
    if comment_payload.is_empty() || comment_payload.last().copied() != Some(0) {
        return Err(AudioArtifactError::MalformedWav(
            "ICMT is not NUL terminated",
        ));
    }
    let text = &comment_payload[..comment_payload.len() - 1];
    if text.len() > budget.maximum_metadata_bytes || text.len() > MAX_WAV_COMMENT_BYTES {
        return Err(AudioArtifactError::BudgetExceeded {
            artifact: "decoded WAV metadata bytes",
            requested: text.len() as u64,
            limit: budget.maximum_metadata_bytes.min(MAX_WAV_COMMENT_BYTES) as u64,
        });
    }
    if text.contains(&0) {
        return Err(AudioArtifactError::InvalidMetadata);
    }
    let decoded = core::str::from_utf8(text).map_err(|_| AudioArtifactError::InvalidMetadata)?;
    let mut comment = String::new();
    comment
        .try_reserve_exact(decoded.len())
        .map_err(|_| AudioArtifactError::Capacity {
            artifact: "decoded WAV metadata",
            requested: decoded.len() as u64,
        })?;
    comment.push_str(decoded);
    WavMetadata::try_new(Some(comment))
}

fn read_u16(bytes: &[u8], offset: usize) -> Result<u16, AudioArtifactError> {
    let value = bytes
        .get(offset..offset + 2)
        .ok_or(AudioArtifactError::MalformedWav("truncated u16 field"))?;
    Ok(u16::from_le_bytes([value[0], value[1]]))
}

fn read_u32(bytes: &[u8], offset: usize) -> Result<u32, AudioArtifactError> {
    let value = bytes
        .get(offset..offset + 4)
        .ok_or(AudioArtifactError::MalformedWav("truncated u32 field"))?;
    Ok(u32::from_le_bytes([value[0], value[1], value[2], value[3]]))
}

fn measure_peaks(samples: &[StereoSample], cx: &Cx<'_>) -> Result<(f64, f64), AudioArtifactError> {
    let mut sample_peak = 0.0_f64;
    for (frame_index, sample) in samples.iter().copied().enumerate() {
        if frame_index % AUDIO_ARTIFACT_CANCELLATION_POLL_FRAMES == 0 {
            checkpoint(cx)?;
        }
        validate_stereo_sample(sample, frame_index as u64)?;
        sample_peak = sample_peak
            .max(sample.left_fs.abs())
            .max(sample.right_fs.abs());
    }
    if samples.len() < 2 {
        return Ok((sample_peak, sample_peak));
    }
    let coefficients = peak_interpolator_coefficients();
    let mut estimate = sample_peak;
    for interval in 0..samples.len() - 1 {
        if interval % AUDIO_ARTIFACT_CANCELLATION_POLL_FRAMES == 0 {
            checkpoint(cx)?;
        }
        for phase in 0..AUDIO_TRUE_PEAK_OVERSAMPLE_FACTOR - 1 {
            let mut left = CompensatedSum::new();
            let mut right = CompensatedSum::new();
            for (tap, coefficient) in coefficients[phase].iter().copied().enumerate() {
                let offset = tap as i64 - (PEAK_INTERPOLATOR_RADIUS - 1);
                let source = reflect_half_sample_even(
                    interval as i128 + i128::from(offset),
                    samples.len() as u64,
                ) as usize;
                left.add(coefficient * samples[source].left_fs);
                right.add(coefficient * samples[source].right_fs);
            }
            estimate = estimate.max(left.total().abs()).max(right.total().abs());
        }
    }
    checkpoint(cx)?;
    Ok((sample_peak, estimate))
}

fn peak_interpolator_coefficients() -> [[f64; 16]; 3] {
    let mut phases = [[0.0_f64; 16]; 3];
    for (phase_index, coefficients) in phases.iter_mut().enumerate() {
        let fraction = (phase_index + 1) as f64 / AUDIO_TRUE_PEAK_OVERSAMPLE_FACTOR as f64;
        let mut normalization = CompensatedSum::new();
        for (tap, coefficient) in coefficients.iter_mut().enumerate() {
            let offset = tap as i64 - (PEAK_INTERPOLATOR_RADIUS - 1);
            let distance = offset as f64 - fraction;
            let value = sinc(distance) * sinc(distance / PEAK_INTERPOLATOR_RADIUS as f64);
            *coefficient = value;
            normalization.add(value);
        }
        let normalization = normalization.total();
        for coefficient in coefficients {
            *coefficient /= normalization;
        }
    }
    phases
}

fn sinc(value: f64) -> f64 {
    if value == 0.0 {
        1.0
    } else {
        let argument = core::f64::consts::PI * value;
        det::sin(argument) / argument
    }
}

fn reflect_half_sample_even(index: i128, length: u64) -> u64 {
    debug_assert!(length > 0);
    let length = i128::from(length);
    let period = 2 * length;
    let reduced = index.rem_euclid(period);
    let reflected = if reduced < length {
        reduced
    } else {
        period - 1 - reduced
    };
    reflected as u64
}

#[derive(Debug, Clone, Copy)]
struct LoudnessResult {
    integrated_lufs: Option<f64>,
    total_blocks: u64,
    absolute_blocks: u64,
    relative_blocks: u64,
}

fn measure_programme_loudness(
    samples: &[StereoSample],
    cx: &Cx<'_>,
) -> Result<LoudnessResult, AudioArtifactError> {
    if samples.len() < LOUDNESS_BLOCK_FRAMES {
        return Ok(LoudnessResult {
            integrated_lufs: None,
            total_blocks: 0,
            absolute_blocks: 0,
            relative_blocks: 0,
        });
    }
    let mut left_shelf = KWeightingShelf::new();
    let mut left_high_pass = KWeightingHighPass::new();
    let mut right_shelf = KWeightingShelf::new();
    let mut right_high_pass = KWeightingHighPass::new();
    let mut energies = Vec::new();
    reserve_exact(&mut energies, samples.len(), "K-weighted sample energies")?;
    for (frame_index, sample) in samples.iter().copied().enumerate() {
        if frame_index % AUDIO_ARTIFACT_CANCELLATION_POLL_FRAMES == 0 {
            checkpoint(cx)?;
        }
        let left = left_high_pass.process(left_shelf.process(sample.left_fs));
        let right = right_high_pass.process(right_shelf.process(sample.right_fs));
        let energy = left.mul_add(left, right * right);
        if !energy.is_finite() {
            return Err(AudioArtifactError::NonFiniteSample {
                frame: frame_index as u64,
                channel: "K-weighted energy",
            });
        }
        energies.push(energy);
    }
    let block_count = 1 + (samples.len() - LOUDNESS_BLOCK_FRAMES) / LOUDNESS_HOP_FRAMES;
    let mut blocks = Vec::new();
    reserve_exact(&mut blocks, block_count, "programme loudness blocks")?;
    let mut block_sum_items = 0_usize;
    for block_index in 0..block_count {
        let start = block_index * LOUDNESS_HOP_FRAMES;
        let end = start + LOUDNESS_BLOCK_FRAMES;
        let mut sum = CompensatedSum::new();
        for value in &energies[start..end] {
            if block_sum_items % AUDIO_ARTIFACT_CANCELLATION_POLL_FRAMES == 0 {
                checkpoint(cx)?;
            }
            sum.add(*value);
            block_sum_items += 1;
        }
        let mean_energy = (sum.total() / LOUDNESS_BLOCK_FRAMES as f64).max(0.0);
        blocks.push(mean_energy);
    }
    let mut absolute_energy = Vec::new();
    reserve_exact(
        &mut absolute_energy,
        blocks.len(),
        "absolute-gated loudness blocks",
    )?;
    for (block_index, energy) in blocks.iter().copied().enumerate() {
        if block_index % AUDIO_ARTIFACT_CANCELLATION_POLL_FRAMES == 0 {
            checkpoint(cx)?;
        }
        if energy > 0.0 && loudness_from_energy(energy) > LOUDNESS_ABSOLUTE_GATE_LUFS {
            absolute_energy.push(energy);
        }
    }
    if absolute_energy.is_empty() {
        return Ok(LoudnessResult {
            integrated_lufs: None,
            total_blocks: block_count as u64,
            absolute_blocks: 0,
            relative_blocks: 0,
        });
    }
    let absolute_mean = compensated_mean(&absolute_energy, cx)?;
    let relative_gate = (loudness_from_energy(absolute_mean) + LOUDNESS_RELATIVE_GATE_LU)
        .max(LOUDNESS_ABSOLUTE_GATE_LUFS);
    let mut final_sum = CompensatedSum::new();
    let mut final_count = 0_u64;
    for (block_index, energy) in absolute_energy.iter().copied().enumerate() {
        if block_index % AUDIO_ARTIFACT_CANCELLATION_POLL_FRAMES == 0 {
            checkpoint(cx)?;
        }
        if loudness_from_energy(energy) > relative_gate {
            final_sum.add(energy);
            final_count += 1;
        }
    }
    let integrated_lufs = if final_count == 0 {
        None
    } else {
        Some(loudness_from_energy(final_sum.total() / final_count as f64))
    };
    checkpoint(cx)?;
    Ok(LoudnessResult {
        integrated_lufs,
        total_blocks: block_count as u64,
        absolute_blocks: absolute_energy.len() as u64,
        relative_blocks: final_count,
    })
}

fn loudness_from_energy(energy: f64) -> f64 {
    LOUDNESS_OFFSET_LUFS + 10.0 * det::ln(energy) / LN_10
}

fn compensated_mean(values: &[f64], cx: &Cx<'_>) -> Result<f64, AudioArtifactError> {
    let mut sum = CompensatedSum::new();
    for (index, value) in values.iter().enumerate() {
        if index % AUDIO_ARTIFACT_CANCELLATION_POLL_FRAMES == 0 {
            checkpoint(cx)?;
        }
        sum.add(*value);
    }
    Ok(sum.total() / values.len() as f64)
}

#[derive(Debug, Clone, Copy)]
struct Biquad {
    b0: f64,
    b1: f64,
    b2: f64,
    a1: f64,
    a2: f64,
    x1: f64,
    x2: f64,
    y1: f64,
    y2: f64,
}

impl Biquad {
    fn process(&mut self, input: f64) -> f64 {
        let output = self.b0.mul_add(
            input,
            self.b1.mul_add(
                self.x1,
                self.b2
                    .mul_add(self.x2, (-self.a1).mul_add(self.y1, -self.a2 * self.y2)),
            ),
        );
        self.x2 = self.x1;
        self.x1 = input;
        self.y2 = self.y1;
        self.y1 = output;
        output
    }
}

struct KWeightingShelf;

impl KWeightingShelf {
    fn new() -> Biquad {
        Biquad {
            b0: 1.535_124_859_586_97,
            b1: -2.691_696_189_406_38,
            b2: 1.198_392_810_852_85,
            a1: -1.690_659_293_182_41,
            a2: 0.732_480_774_215_85,
            x1: 0.0,
            x2: 0.0,
            y1: 0.0,
            y2: 0.0,
        }
    }
}

struct KWeightingHighPass;

impl KWeightingHighPass {
    fn new() -> Biquad {
        Biquad {
            b0: 1.0,
            b1: -2.0,
            b2: 1.0,
            a1: -1.990_047_454_833_98,
            a2: 0.990_072_250_366_21,
            x1: 0.0,
            x2: 0.0,
            y1: 0.0,
            y2: 0.0,
        }
    }
}

fn build_manifest(
    configuration: &SoundSynthesisConfig,
    channel_layout: AudioChannelLayoutReceipt,
    source_signal_identity: ContentHash,
    mix_identity: Option<ContentHash>,
    wav: WavCodecReceipt,
    meters: AudioMeters,
    admitted_headroom_db: f64,
) -> Result<AudioArtifactManifest, AudioArtifactError> {
    let input = configuration.input();
    let video_frames = i128::from(input.video_clock.end_tick_exclusive())
        - i128::from(input.video_clock.start_tick());
    let audio_frames = i128::from(input.audio_clock.end_tick_exclusive())
        - i128::from(input.audio_clock.start_tick());
    if video_frames < 0
        || audio_frames < 0
        || audio_frames != video_frames * 2_000
        || wav.sample_frame_count != audio_frames as u64
    {
        return Err(AudioArtifactError::SampleCountMismatch {
            expected: audio_frames.max(0) as u64,
            actual: wav.sample_frame_count,
        });
    }
    let mut manifest = AudioArtifactManifest {
        identity: zero_hash(),
        synthesis: configuration.receipt(),
        authority: configuration.authority(),
        channel_layout,
        source_signal_identity,
        mix_identity,
        wav,
        role: AudioArtifactRole::for_encoding(wav.encoding),
        meters,
        video_start_tick: input.video_clock.start_tick(),
        video_end_tick_exclusive: input.video_clock.end_tick_exclusive(),
        audio_start_tick: input.audio_clock.start_tick(),
        audio_end_tick_exclusive: input.audio_clock.end_tick_exclusive(),
        audio_frames_per_video_frame: 2_000,
        admitted_headroom_db,
    };
    manifest.identity = manifest_identity(&manifest);
    Ok(manifest)
}

fn manifest_identity(manifest: &AudioArtifactManifest) -> ContentHash {
    let mut hasher = DomainHasher::new(AUDIO_MANIFEST_IDENTITY_DOMAIN);
    hasher.update(&AUDIO_ARTIFACT_SCHEMA_VERSION.to_le_bytes());
    hash_synthesis_receipt(&mut hasher, manifest.synthesis);
    hasher.update(manifest.authority.code().as_bytes());
    hasher.update(manifest.channel_layout.identity.as_bytes());
    hasher.update(manifest.source_signal_identity.as_bytes());
    match manifest.mix_identity {
        None => hasher.update(&[0]),
        Some(identity) => {
            hasher.update(&[1]);
            hasher.update(identity.as_bytes());
        }
    }
    hash_wav_receipt(&mut hasher, manifest.wav);
    hasher.update(manifest.role.code().as_bytes());
    hash_meters(&mut hasher, manifest.meters);
    hasher.update(&manifest.video_start_tick.to_le_bytes());
    hasher.update(&manifest.video_end_tick_exclusive.to_le_bytes());
    hasher.update(&manifest.audio_start_tick.to_le_bytes());
    hasher.update(&manifest.audio_end_tick_exclusive.to_le_bytes());
    hasher.update(&manifest.audio_frames_per_video_frame.to_le_bytes());
    hasher.update(&manifest.admitted_headroom_db.to_bits().to_le_bytes());
    hasher.finalize()
}

fn hash_synthesis_receipt(hasher: &mut DomainHasher, receipt: SoundSynthesisReceipt) {
    hasher.update(&receipt.schema_version.to_le_bytes());
    hasher.update(receipt.configuration_identity.as_bytes());
    hasher.update(receipt.authority.code().as_bytes());
    hasher.update(receipt.trajectory_identity.as_bytes());
    hasher.update(receipt.excitation_identity.as_bytes());
    hasher.update(receipt.sound_model_identity.as_bytes());
    hasher.update(receipt.timeline_identity.as_bytes());
}

fn hash_wav_receipt(hasher: &mut DomainHasher, receipt: WavCodecReceipt) {
    hasher.update(receipt.wav_identity.as_bytes());
    hasher.update(receipt.metadata_identity.as_bytes());
    hasher.update(&receipt.byte_len.to_le_bytes());
    hasher.update(&receipt.sample_frame_count.to_le_bytes());
    hasher.update(&receipt.sample_rate_hz.to_le_bytes());
    hasher.update(&(receipt.encoding as u16).to_le_bytes());
}

fn hash_meters(hasher: &mut DomainHasher, meters: AudioMeters) {
    hasher.update(&meters.sample_peak_fs.to_bits().to_le_bytes());
    hasher.update(&meters.true_peak_estimate_fs.to_bits().to_le_bytes());
    hasher.update(&meters.stereo_rms_fs.to_bits().to_le_bytes());
    hasher.update(&meters.dc_left_fs.to_bits().to_le_bytes());
    hasher.update(&meters.dc_right_fs.to_bits().to_le_bytes());
    match meters.integrated_loudness_lufs {
        None => hasher.update(&[0]),
        Some(value) => {
            hasher.update(&[1]);
            hasher.update(&value.to_bits().to_le_bytes());
        }
    }
    hasher.update(&meters.loudness_block_count.to_le_bytes());
    hasher.update(&meters.absolute_gated_block_count.to_le_bytes());
    hasher.update(&meters.relative_gated_block_count.to_le_bytes());
}

fn hash_bytes_cancellable(
    domain: &str,
    bytes: &[u8],
    budget: AudioArtifactBudget,
    artifact: &'static str,
    cx: &Cx<'_>,
) -> Result<ContentHash, AudioArtifactError> {
    check_combined_work(&[], bytes.len() as u64, budget, artifact)?;
    hash_bytes_with_checkpoint(domain, bytes, &mut || checkpoint(cx))
}

fn hash_bytes_with_checkpoint(
    domain: &str,
    bytes: &[u8],
    checkpoint_fn: &mut impl FnMut() -> Result<(), AudioArtifactError>,
) -> Result<ContentHash, AudioArtifactError> {
    let mut hasher = DomainHasher::new(domain);
    for chunk in bytes.chunks(AUDIO_ARTIFACT_CANCELLATION_POLL_BYTES) {
        checkpoint_fn()?;
        hasher.update(chunk);
    }
    checkpoint_fn()?;
    Ok(hasher.finalize())
}

fn zero_hash() -> ContentHash {
    ContentHash([0; 32])
}

fn is_zero(identity: ContentHash) -> bool {
    identity.as_bytes().iter().all(|byte| *byte == 0)
}

fn checkpoint(cx: &Cx<'_>) -> Result<(), AudioArtifactError> {
    cx.checkpoint().map_err(|_| AudioArtifactError::Cancelled)
}

#[derive(Debug, Clone, Copy)]
struct CompensatedSum {
    sum: f64,
    correction: f64,
}

impl CompensatedSum {
    const fn new() -> Self {
        Self {
            sum: 0.0,
            correction: 0.0,
        }
    }

    fn add(&mut self, value: f64) {
        let next = self.sum + value;
        if self.sum.abs() >= value.abs() {
            self.correction += (self.sum - next) + value;
        } else {
            self.correction += (value - next) + self.sum;
        }
        self.sum = next;
    }

    fn total(self) -> f64 {
        self.sum + self.correction
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn g4_g5_chunked_hash_is_replay_stable_and_stops_at_cancelled_boundary() {
        let bytes = vec![0x5a; AUDIO_ARTIFACT_CANCELLATION_POLL_BYTES * 2 + 1];
        let mut expected_hasher = DomainHasher::new(WAV_IDENTITY_DOMAIN);
        expected_hasher.update(&bytes);
        let expected = expected_hasher.finalize();

        let mut completed_checkpoints = 0_usize;
        let actual = hash_bytes_with_checkpoint(WAV_IDENTITY_DOMAIN, &bytes, &mut || {
            completed_checkpoints += 1;
            Ok(())
        })
        .expect("uninterrupted chunked hash must complete");
        assert_eq!(actual, expected, "chunking must not change BLAKE3 identity");
        assert_eq!(completed_checkpoints, 4, "three chunks plus final poll");

        let mut cancelled_checkpoints = 0_usize;
        let cancelled = hash_bytes_with_checkpoint(WAV_IDENTITY_DOMAIN, &bytes, &mut || {
            cancelled_checkpoints += 1;
            if cancelled_checkpoints == 2 {
                Err(AudioArtifactError::Cancelled)
            } else {
                Ok(())
            }
        });
        assert_eq!(cancelled, Err(AudioArtifactError::Cancelled));
        assert_eq!(cancelled_checkpoints, 2);
    }
}
