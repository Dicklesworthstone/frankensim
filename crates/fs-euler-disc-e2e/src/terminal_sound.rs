//! Honest terminal-sound policy and the optional artistic singularity chirp
//! (bead frankensim-h7xu5.7.7).
//!
//! The motivating disc makes an accelerating terminal sound associated with
//! finite-time-singularity models, but the reduced runner terminates at an
//! inclination threshold and does not resolve the final elastic / contact /
//! thin-air-film acoustics. Extrapolating it into audible infinity would be
//! misleading, so the audible ending is a POLICY with typed authority:
//!
//! - [`TerminalSoundPolicy::TaperAtLastSupported`] (default): the existing
//!   deterministic taper at the last supported sample; zero synthesized
//!   samples, and the output equals the supported mix exactly.
//! - [`TerminalSoundPolicy::ArtisticChirp`]: a frozen deterministic chirp
//!   law driven by the final supported trend, admitted ONLY where the
//!   disposition supports it, with every synthesized sample marked in the
//!   receipt.
//!
//! Disposition gating is the honesty core:
//! - `TerminalInclination` (an observed physical stop): the chirp is
//!   admissible as declared artistry on top of a supported terminal event.
//! - `HorizonReached` (censored data): a "stop" sound would invent an
//!   ending the simulation never observed. The chirp refuses unless the
//!   caller declares [`ChirpRequest::explicit_edit_intent`], and even then
//!   the receipt records the intent.
//! - `NumericalRefusal`: never decorated. Taper only.
//!
//! No-claims: a chirped ending is presentation, not physics. The receipt's
//! `synthesized_sample_frames` is the exact count of samples that are NOT
//! backed by supported trajectory time; downstream mixes must carry it.

use crate::coupled_runner::CoupledTerminal;
use fs_evidence::cinematic_sound::SOUND_MASTER_SAMPLE_RATE_HZ;

/// Frozen chirp frequency cap: comfortably below the 24 kHz Nyquist of the
/// 48 kHz master, so the band limit holds by construction (no partial is
/// ever synthesized above it; nothing to alias).
pub const CHIRP_FREQUENCY_CAP_HZ: f64 = 16_000.0;
/// Frozen maximum chirp duration.
pub const CHIRP_MAX_DURATION_S: f64 = 1.25;
/// Frozen minimum admitted chirp start frequency.
pub const CHIRP_MIN_START_HZ: f64 = 40.0;
/// Frozen raised-cosine crossfade from the supported stem into the chirp.
pub const CHIRP_CROSSFADE_SAMPLE_FRAMES: u32 = 480; // 10 ms at 48 kHz
/// Frozen terminal release so the chirp never ends on a click.
pub const CHIRP_RELEASE_SAMPLE_FRAMES: u32 = 240; // 5 ms at 48 kHz
/// Frozen peak amplitude of the chirp relative to full scale.
pub const CHIRP_PEAK_FS: f64 = 0.5;

/// The three declared terminal policies.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TerminalSoundPolicy {
    /// End at the last supported sample with the existing deterministic
    /// taper. Always admissible; synthesizes nothing.
    TaperAtLastSupported,
    /// Bounded extrapolation under an explicit parametric model is a
    /// DECLARED policy without an admitted model: no parametric terminal
    /// model has passed validation, so selecting it refuses with
    /// [`TerminalSoundError::NoAdmittedExtrapolationModel`] rather than
    /// synthesizing under an unvalidated law.
    BoundedExtrapolation,
    /// The artistic singularity chirp under the frozen law.
    ArtisticChirp,
}

/// Caller request for the artistic chirp.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ChirpRequest {
    /// Instantaneous drive frequency at the last supported sample, Hz.
    /// Derived from the final supported trend by the caller; this module
    /// does not inspect the trajectory.
    pub final_supported_frequency_hz: f64,
    /// Fractional growth of the drive frequency per second measured over
    /// the final supported window (>= 0; zero trend means a flat chirp).
    pub terminal_growth_per_s: f64,
    /// Chirp peak amplitude in the SUPPORTED STEM'S OWN units (full scale
    /// for an fs-domain stem, pascals for a pressure stem: the caller owns
    /// the unit conversion, e.g. `CHIRP_PEAK_FS / gain_fs_per_pa`).
    /// Downstream loudness/peak gates remain the clipping authority.
    pub peak_amplitude: f64,
    /// Explicit edit intent, required to chirp censored (horizon) data.
    pub explicit_edit_intent: bool,
}

/// Typed refusals of the terminal-sound policy boundary.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TerminalSoundError {
    /// BoundedExtrapolation has no admitted parametric model.
    NoAdmittedExtrapolationModel,
    /// A chirp on censored data without explicit edit intent invents an
    /// ending the simulation never observed.
    CensoredWithoutEditIntent,
    /// A numerical refusal is never decorated with an artistic ending.
    RefusalNeverDecorated,
    /// The request's trend fields are outside the admitted domain.
    InvalidTrend,
    /// The supported stem is too short to carry the declared crossfade.
    StemTooShortForCrossfade,
}

impl core::fmt::Display for TerminalSoundError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::NoAdmittedExtrapolationModel => f.write_str(
                "bounded extrapolation is declared but no parametric terminal model is admitted",
            ),
            Self::CensoredWithoutEditIntent => f.write_str(
                "horizon-censored trajectory: an artistic stop sound requires explicit edit intent",
            ),
            Self::RefusalNeverDecorated => {
                f.write_str("numerical-refusal trajectory: terminal sound is taper-only")
            }
            Self::InvalidTrend => f.write_str(
                "chirp trend must be finite, start frequency in the admitted band, growth >= 0",
            ),
            Self::StemTooShortForCrossfade => {
                f.write_str("supported stem shorter than the declared crossfade")
            }
        }
    }
}

impl std::error::Error for TerminalSoundError {}

/// Authority receipt for one terminal-sound decision. Every synthesized
/// sample beyond supported trajectory time is counted here; a mix that
/// drops this record has laundered artistry into physics.
#[derive(Debug, Clone, PartialEq)]
pub struct TerminalSoundReceipt {
    /// Which policy actually ran.
    pub policy: TerminalSoundPolicy,
    /// The trajectory disposition the policy was admitted against.
    pub disposition: TerminalDispositionClass,
    /// Sample frames of the SUPPORTED stem consumed by the crossfade.
    pub crossfade_sample_frames: u32,
    /// Exact count of synthesized sample frames not backed by supported
    /// trajectory time (zero under taper).
    pub synthesized_sample_frames: u64,
    /// Chirp law parameters when the chirp ran.
    pub chirp: Option<ChirpLawReceipt>,
    /// Whether explicit edit intent was declared (censored chirps only).
    pub explicit_edit_intent: bool,
}

/// Coarse disposition class retained in the receipt.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TerminalDispositionClass {
    /// Observed physical stop.
    TerminalInclination,
    /// Horizon elapsed with the run still live (censored).
    HorizonCensored,
    /// The trajectory refused to continue.
    NumericalRefusal,
}

impl TerminalDispositionClass {
    /// Classify the coupled-runner terminal vocabulary.
    #[must_use]
    pub const fn from_coupled(terminal: &CoupledTerminal) -> Self {
        match terminal {
            CoupledTerminal::TerminalInclination => Self::TerminalInclination,
            CoupledTerminal::HorizonReached => Self::HorizonCensored,
            CoupledTerminal::NumericalRefusal { .. } => Self::NumericalRefusal,
        }
    }
}

/// Frozen chirp-law record: everything needed to regenerate the synthetic
/// tail bit-exactly.
#[derive(Debug, Clone, PartialEq)]
pub struct ChirpLawReceipt {
    /// Start frequency, Hz (clamped into the admitted band).
    pub start_hz: f64,
    /// Exponential growth rate per second actually applied.
    pub growth_per_s: f64,
    /// Frequency cap actually reached or respected.
    pub cap_hz: f64,
    /// Sample frame at which the instantaneous frequency first hit the cap
    /// (`None` when the cap was never reached inside the duration).
    pub cap_crossing_frame: Option<u64>,
    /// Total chirp duration in sample frames.
    pub duration_sample_frames: u64,
}

/// Apply the terminal-sound policy to a supported mono stem at the master
/// rate. Returns the full output (supported stem, possibly crossfaded into
/// a synthetic tail) plus the authority receipt.
///
/// Under [`TerminalSoundPolicy::TaperAtLastSupported`] the returned samples
/// are EXACTLY the input (the production taper already lives in the
/// listening-master policy); the receipt records zero synthesized frames —
/// disabling the artistic layer yields the supported mix bit-exactly.
///
/// # Errors
/// Typed refusals per the disposition-gating doctrine above.
pub fn apply_terminal_policy(
    supported_stem: &[f64],
    disposition: TerminalDispositionClass,
    policy: TerminalSoundPolicy,
    request: Option<ChirpRequest>,
) -> Result<(Vec<f64>, TerminalSoundReceipt), TerminalSoundError> {
    match policy {
        TerminalSoundPolicy::TaperAtLastSupported => Ok((
            supported_stem.to_vec(),
            TerminalSoundReceipt {
                policy,
                disposition,
                crossfade_sample_frames: 0,
                synthesized_sample_frames: 0,
                chirp: None,
                explicit_edit_intent: false,
            },
        )),
        TerminalSoundPolicy::BoundedExtrapolation => {
            Err(TerminalSoundError::NoAdmittedExtrapolationModel)
        }
        TerminalSoundPolicy::ArtisticChirp => {
            let request = request.ok_or(TerminalSoundError::InvalidTrend)?;
            match disposition {
                TerminalDispositionClass::NumericalRefusal => {
                    return Err(TerminalSoundError::RefusalNeverDecorated);
                }
                TerminalDispositionClass::HorizonCensored if !request.explicit_edit_intent => {
                    return Err(TerminalSoundError::CensoredWithoutEditIntent);
                }
                _ => {}
            }
            if !(request.final_supported_frequency_hz.is_finite()
                && request.terminal_growth_per_s.is_finite()
                && request.terminal_growth_per_s >= 0.0
                && request.final_supported_frequency_hz >= CHIRP_MIN_START_HZ
                && request.final_supported_frequency_hz <= CHIRP_FREQUENCY_CAP_HZ
                && request.peak_amplitude.is_finite()
                && request.peak_amplitude > 0.0)
            {
                return Err(TerminalSoundError::InvalidTrend);
            }
            let crossfade = CHIRP_CROSSFADE_SAMPLE_FRAMES as usize;
            if supported_stem.len() < crossfade {
                return Err(TerminalSoundError::StemTooShortForCrossfade);
            }
            let (tail, law) = synthesize_chirp(&request);
            let mut output = supported_stem.to_vec();
            let splice = output.len() - crossfade;
            // Raised-cosine crossfade: the supported stem eases out while
            // the chirp eases in over the same frames, so the first
            // difference at the splice stays bounded (no click).
            for index in 0..crossfade {
                let phase = core::f64::consts::PI * (index as f64) / (crossfade as f64);
                let fade_out = 0.5 * (1.0 + phase.cos());
                let fade_in = 1.0 - fade_out;
                output[splice + index] = output[splice + index] * fade_out + tail[index] * fade_in;
            }
            output.extend_from_slice(&tail[crossfade..]);
            let synthesized = (tail.len() - crossfade) as u64;
            Ok((
                output,
                TerminalSoundReceipt {
                    policy,
                    disposition,
                    crossfade_sample_frames: CHIRP_CROSSFADE_SAMPLE_FRAMES,
                    synthesized_sample_frames: synthesized,
                    chirp: Some(law),
                    explicit_edit_intent: request.explicit_edit_intent,
                },
            ))
        }
    }
}

/// Deterministic band-limited chirp under the frozen law: exponential
/// frequency growth from the final supported trend, hard cap below
/// Nyquist, exponential amplitude decay with a raised-cosine release.
fn synthesize_chirp(request: &ChirpRequest) -> (Vec<f64>, ChirpLawReceipt) {
    let rate = f64::from(SOUND_MASTER_SAMPLE_RATE_HZ);
    let duration_frames = (CHIRP_MAX_DURATION_S * rate) as u64;
    let dt = rate.recip();
    let start_hz = request.final_supported_frequency_hz;
    let growth = request.terminal_growth_per_s;
    // Amplitude decays so the chirp fades toward the singularity instead
    // of screaming: half-life tied to the chirp duration.
    let amplitude_decay_per_s = 4.0 / CHIRP_MAX_DURATION_S;
    let peak = request.peak_amplitude;
    let release_start = duration_frames.saturating_sub(u64::from(CHIRP_RELEASE_SAMPLE_FRAMES));
    let mut samples = Vec::with_capacity(usize::try_from(duration_frames).unwrap_or(0));
    let mut phase = 0.0f64;
    let mut cap_crossing_frame = None;
    for frame in 0..duration_frames {
        let time_s = frame as f64 * dt;
        let mut frequency = start_hz * (growth * time_s).exp();
        if frequency >= CHIRP_FREQUENCY_CAP_HZ {
            frequency = CHIRP_FREQUENCY_CAP_HZ;
            if cap_crossing_frame.is_none() {
                cap_crossing_frame = Some(frame);
            }
        }
        phase += 2.0 * core::f64::consts::PI * frequency * dt;
        let mut amplitude = peak * (-amplitude_decay_per_s * time_s).exp();
        if frame >= release_start {
            let into =
                (frame - release_start) as f64 / f64::from(CHIRP_RELEASE_SAMPLE_FRAMES.max(1));
            amplitude *= 0.5 * (1.0 + (core::f64::consts::PI * into).cos());
        }
        samples.push(amplitude * phase.sin());
    }
    (
        samples,
        ChirpLawReceipt {
            start_hz,
            growth_per_s: growth,
            cap_hz: CHIRP_FREQUENCY_CAP_HZ,
            cap_crossing_frame,
            duration_sample_frames: duration_frames,
        },
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    fn stem(frames: usize) -> Vec<f64> {
        // A quiet supported stem with nonzero content so bit-exactness is
        // a real assertion.
        (0..frames)
            .map(|index| 0.25 * ((index as f64) * 0.01).sin())
            .collect()
    }

    fn chirp_request() -> ChirpRequest {
        ChirpRequest {
            final_supported_frequency_hz: 220.0,
            terminal_growth_per_s: 6.0,
            peak_amplitude: CHIRP_PEAK_FS,
            explicit_edit_intent: false,
        }
    }

    #[test]
    fn taper_policy_is_bit_exact_and_synthesizes_nothing() {
        let input = stem(4_800);
        for disposition in [
            TerminalDispositionClass::TerminalInclination,
            TerminalDispositionClass::HorizonCensored,
            TerminalDispositionClass::NumericalRefusal,
        ] {
            let (output, receipt) = apply_terminal_policy(
                &input,
                disposition,
                TerminalSoundPolicy::TaperAtLastSupported,
                None,
            )
            .expect("taper always admits");
            assert_eq!(
                output, input,
                "disabling artistry yields the supported mix exactly"
            );
            assert_eq!(receipt.synthesized_sample_frames, 0);
            assert!(receipt.chirp.is_none());
        }
    }

    #[test]
    fn disposition_gating_is_the_honesty_core() {
        let input = stem(4_800);
        // Refusal trajectories are never decorated.
        assert_eq!(
            apply_terminal_policy(
                &input,
                TerminalDispositionClass::NumericalRefusal,
                TerminalSoundPolicy::ArtisticChirp,
                Some(chirp_request()),
            )
            .expect_err("refusal never decorated"),
            TerminalSoundError::RefusalNeverDecorated
        );
        // Censored data refuses without explicit edit intent...
        assert_eq!(
            apply_terminal_policy(
                &input,
                TerminalDispositionClass::HorizonCensored,
                TerminalSoundPolicy::ArtisticChirp,
                Some(chirp_request()),
            )
            .expect_err("censored without intent refuses"),
            TerminalSoundError::CensoredWithoutEditIntent
        );
        // ...and admits WITH it, recording the intent in the receipt.
        let mut intended = chirp_request();
        intended.explicit_edit_intent = true;
        let (_, receipt) = apply_terminal_policy(
            &input,
            TerminalDispositionClass::HorizonCensored,
            TerminalSoundPolicy::ArtisticChirp,
            Some(intended),
        )
        .expect("explicit intent admits");
        assert!(receipt.explicit_edit_intent);
        // Bounded extrapolation refuses: declared policy, no admitted model.
        assert_eq!(
            apply_terminal_policy(
                &input,
                TerminalDispositionClass::TerminalInclination,
                TerminalSoundPolicy::BoundedExtrapolation,
                None,
            )
            .expect_err("no admitted model"),
            TerminalSoundError::NoAdmittedExtrapolationModel
        );
    }

    #[test]
    fn chirp_is_deterministic_capped_and_click_free() {
        let input = stem(4_800);
        let run = || {
            apply_terminal_policy(
                &input,
                TerminalDispositionClass::TerminalInclination,
                TerminalSoundPolicy::ArtisticChirp,
                Some(chirp_request()),
            )
            .expect("terminal inclination admits the chirp")
        };
        let (first, receipt) = run();
        let (again, _) = run();
        assert_eq!(first, again, "chirp output is deterministic");
        let law = receipt.chirp.expect("law recorded");
        // Extreme trend hits the cap and records the crossing frame.
        assert!(
            law.cap_crossing_frame.is_some(),
            "growth 6/s from 220 Hz must cap"
        );
        assert!(law.cap_hz <= CHIRP_FREQUENCY_CAP_HZ);
        // The receipt counts exactly the frames beyond the supported stem.
        assert_eq!(
            first.len() as u64,
            input.len() as u64 + receipt.synthesized_sample_frames,
        );
        // Click-freedom at the splice and everywhere: the first difference
        // stays below a generous bound for band-limited content at 48 kHz.
        let max_step = first
            .windows(2)
            .map(|pair| (pair[1] - pair[0]).abs())
            .fold(0.0f64, f64::max);
        assert!(
            max_step < 1.2,
            "first difference bounded, no click: {max_step}"
        );
        // The chirp ends in silence (release window).
        assert!(first.last().copied().expect("nonempty").abs() < 1e-6);
    }

    #[test]
    fn zero_trend_admits_and_never_caps() {
        let input = stem(4_800);
        let request = ChirpRequest {
            final_supported_frequency_hz: 440.0,
            terminal_growth_per_s: 0.0,
            peak_amplitude: CHIRP_PEAK_FS,
            explicit_edit_intent: false,
        };
        let (_, receipt) = apply_terminal_policy(
            &input,
            TerminalDispositionClass::TerminalInclination,
            TerminalSoundPolicy::ArtisticChirp,
            Some(request),
        )
        .expect("zero trend admits (flat chirp)");
        let law = receipt.chirp.expect("law");
        assert_eq!(law.cap_crossing_frame, None, "flat 440 Hz never caps");
    }

    #[test]
    fn invalid_trends_and_short_stems_refuse() {
        let input = stem(4_800);
        for bad in [
            ChirpRequest {
                final_supported_frequency_hz: f64::NAN,
                ..chirp_request()
            },
            ChirpRequest {
                final_supported_frequency_hz: 10.0, // below admitted band
                ..chirp_request()
            },
            ChirpRequest {
                final_supported_frequency_hz: 30_000.0, // above cap
                ..chirp_request()
            },
            ChirpRequest {
                terminal_growth_per_s: -1.0,
                ..chirp_request()
            },
        ] {
            assert_eq!(
                apply_terminal_policy(
                    &input,
                    TerminalDispositionClass::TerminalInclination,
                    TerminalSoundPolicy::ArtisticChirp,
                    Some(bad),
                )
                .expect_err("invalid trend refuses"),
                TerminalSoundError::InvalidTrend
            );
        }
        // Missing final state: no request at all.
        assert_eq!(
            apply_terminal_policy(
                &input,
                TerminalDispositionClass::TerminalInclination,
                TerminalSoundPolicy::ArtisticChirp,
                None,
            )
            .expect_err("missing request refuses"),
            TerminalSoundError::InvalidTrend
        );
        // Stem shorter than the crossfade refuses.
        assert_eq!(
            apply_terminal_policy(
                &stem(100),
                TerminalDispositionClass::TerminalInclination,
                TerminalSoundPolicy::ArtisticChirp,
                Some(chirp_request()),
            )
            .expect_err("short stem refuses"),
            TerminalSoundError::StemTooShortForCrossfade
        );
    }

    #[test]
    fn disposition_classification_matches_the_coupled_vocabulary() {
        use crate::coupled_runner::CoupledNumericalRefusalReason;
        assert_eq!(
            TerminalDispositionClass::from_coupled(&CoupledTerminal::TerminalInclination),
            TerminalDispositionClass::TerminalInclination
        );
        assert_eq!(
            TerminalDispositionClass::from_coupled(&CoupledTerminal::HorizonReached),
            TerminalDispositionClass::HorizonCensored
        );
        assert_eq!(
            TerminalDispositionClass::from_coupled(&CoupledTerminal::NumericalRefusal {
                reason: CoupledNumericalRefusalReason::ReimpactLimitExceeded,
            }),
            TerminalDispositionClass::NumericalRefusal
        );
    }
}
