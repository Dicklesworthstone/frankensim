//! Typed bridge from one acoustic pressure artifact to a cinematic stem.
//!
//! The bridge is deliberately a boundary, not a room/mastering or PCM path:
//! it consumes already sampled time-domain pascals without altering values.
//! Source, clock, observer, bandwidth, uncertainty, and authority remain
//! bound to the returned stem so a caller cannot silently promote a synthetic
//! waveform into a calibrated presentation input.

use fs_blake3::{ContentHash, DomainHasher, hash_domain};

/// Semantic version of the pressure-to-cinematic-stem boundary.
pub const CINEMATIC_STEM_SEMANTICS_VERSION: u32 = 1;

const OBSERVER_IDENTITY_DOMAIN: &str = "org.frankensim.fs-couple.acoustic-observer.v1";
const PRESSURE_PAYLOAD_IDENTITY_DOMAIN: &str =
    "org.frankensim.fs-couple.acoustic-pressure-payload.v1";
const STEM_IDENTITY_DOMAIN: &str = "org.frankensim.fs-couple.cinematic-pressure-stem.v1";

/// Opaque proof that calibration authority was admitted by this boundary's
/// producer path. It has no public constructor, so callers can request a
/// calibrated stem tier but cannot mint a calibrated source artifact.
///
/// ```compile_fail
/// use fs_couple::cinematic_stem::AcousticAuthority;
///
/// // Calibration is not a caller-mintable unit enum variant.
/// let _forged = AcousticAuthority::Calibrated;
/// ```
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct CalibratedAuthority {
    _sealed: (),
}

impl CalibratedAuthority {
    const INTERNAL: Self = Self { _sealed: () };
}

/// Authority retained by an acoustic pressure artifact.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum AcousticAuthority {
    /// Analytic or synthetic fixture data. Useful for integration, never a
    /// calibration claim.
    Synthetic,
    /// A source with declared limitations but no calibration authority.
    Estimated,
    /// A source whose calibration authority was admitted upstream.
    Calibrated(CalibratedAuthority),
}

impl AcousticAuthority {
    const fn identity_tag(self) -> u8 {
        match self {
            Self::Synthetic => 0,
            Self::Estimated => 1,
            Self::Calibrated(_) => 2,
        }
    }
}

/// A consumer's minimum acoustic authority; this is intentionally distinct
/// from the producer-carried [`AcousticAuthority`]. A caller may require the
/// calibrated tier without being able to manufacture it on an artifact.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CinematicAuthorityRequirement {
    /// Accept a synthetic, estimated, or calibrated artifact.
    Synthetic,
    /// Accept an estimated or calibrated artifact.
    Estimated,
    /// Accept only an upstream-admitted calibrated artifact.
    Calibrated,
}

impl CinematicAuthorityRequirement {
    const fn is_satisfied_by(self, authority: AcousticAuthority) -> bool {
        matches!(
            (self, authority),
            (Self::Synthetic, _)
                | (
                    Self::Estimated,
                    AcousticAuthority::Estimated | AcousticAuthority::Calibrated(_)
                )
                | (Self::Calibrated, AcousticAuthority::Calibrated(_))
        )
    }
}

/// Declared representation of the supplied acoustic values.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PressureRepresentation {
    /// Real time-domain pressure samples in pascals.
    TimeDomainPascals,
    /// Frequency-domain phasors reported with peak amplitude.
    FrequencyPhasorPeak,
    /// Frequency-domain phasors reported with RMS amplitude.
    FrequencyPhasorRms,
}

/// Immutable observation geometry for one pressure artifact.
#[derive(Clone, Debug, PartialEq)]
pub struct AcousticObserver {
    position_m: [f64; 3],
    frame_identity: ContentHash,
    identity: ContentHash,
}

impl AcousticObserver {
    /// Construct an observation geometry in metres and bind it to its frame.
    ///
    /// # Errors
    /// Refuses non-finite coordinates or a missing frame identity.
    pub fn try_new(
        position_m: [f64; 3],
        frame_identity: ContentHash,
    ) -> Result<Self, CinematicStemError> {
        if frame_identity == ContentHash([0; 32]) {
            return Err(CinematicStemError::MissingIdentity {
                field: "observer frame",
            });
        }
        if position_m.iter().any(|coordinate| !coordinate.is_finite()) {
            return Err(CinematicStemError::NonFinite {
                field: "observer position",
            });
        }
        let position_m = position_m.map(canonical_zero);
        let mut payload = Vec::with_capacity(4 + 32 + 3 * 8);
        payload.extend_from_slice(&CINEMATIC_STEM_SEMANTICS_VERSION.to_le_bytes());
        payload.extend_from_slice(frame_identity.as_bytes());
        for coordinate in position_m {
            payload.extend_from_slice(&coordinate.to_bits().to_le_bytes());
        }
        Ok(Self {
            position_m,
            frame_identity,
            identity: hash_domain(OBSERVER_IDENTITY_DOMAIN, &payload),
        })
    }

    /// Position in the declared frame, in metres.
    #[must_use]
    pub const fn position_m(&self) -> [f64; 3] {
        self.position_m
    }

    /// Immutable frame identity.
    #[must_use]
    pub const fn frame_identity(&self) -> ContentHash {
        self.frame_identity
    }

    /// Content identity of the complete observation geometry.
    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }
}

/// Typed refusal at the acoustic-to-cinematic boundary.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum CinematicStemError {
    /// A required content binding was the all-zero sentinel.
    MissingIdentity {
        /// Missing binding.
        field: &'static str,
    },
    /// A supplied scalar or sample was non-finite.
    NonFinite {
        /// Invalid field.
        field: &'static str,
    },
    /// A positive, bounded physical value was required.
    InvalidValue {
        /// Invalid field.
        field: &'static str,
    },
    /// No time-domain samples were supplied.
    EmptySamples,
    /// The bridge accepts waveform samples, not ambiguous phasors.
    UnsupportedRepresentation {
        /// Representation that cannot be converted without an explicit
        /// phase/amplitude convention transform.
        observed: PressureRepresentation,
    },
    /// The expected cinematic clock differs from the source clock.
    ClockMismatch,
    /// The expected listener geometry differs from the source observer.
    ObserverMismatch,
    /// The source authority cannot satisfy the requested stem tier.
    InsufficientAuthority {
        /// Authority carried by the source artifact.
        observed: AcousticAuthority,
        /// Minimum authority required by the consumer.
        required: CinematicAuthorityRequirement,
    },
}

impl core::fmt::Display for CinematicStemError {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::MissingIdentity { field } => {
                write!(formatter, "cinematic stem requires {field} identity")
            }
            Self::NonFinite { field } => write!(formatter, "cinematic stem {field} must be finite"),
            Self::InvalidValue { field } => write!(
                formatter,
                "cinematic stem {field} is outside its admitted range"
            ),
            Self::EmptySamples => {
                formatter.write_str("cinematic stem requires at least one pressure sample")
            }
            Self::UnsupportedRepresentation { observed } => write!(
                formatter,
                "cinematic stem refuses {observed:?}; only time-domain pascal samples are admitted"
            ),
            Self::ClockMismatch => formatter
                .write_str("cinematic stem source clock does not match the requested clock"),
            Self::ObserverMismatch => formatter
                .write_str("cinematic stem source observer does not match the requested observer"),
            Self::InsufficientAuthority { observed, required } => write!(
                formatter,
                "cinematic stem source authority {observed:?} is below required {required:?}"
            ),
        }
    }
}

impl std::error::Error for CinematicStemError {}

/// A validated time-domain acoustic artifact. Its samples are physical Pa,
/// not normalized PCM or a mastering-ready signal.
#[derive(Clone, Debug, PartialEq)]
pub struct AcousticPressureArtifact {
    pressure_pa: Vec<f64>,
    sample_rate_hz: u32,
    clock_identity: ContentHash,
    observer: AcousticObserver,
    bandwidth_hz: f64,
    uncertainty_pa: f64,
    authority: AcousticAuthority,
    source_identity: ContentHash,
    pressure_identity: ContentHash,
}

impl AcousticPressureArtifact {
    /// Admit synthetic or analytic time-domain pressure samples in Pa.
    ///
    /// This factory is the public lower-authority producer seam. It fixes the
    /// representation to time-domain pascals and the authority to
    /// [`AcousticAuthority::Synthetic`]; callers cannot select either tag.
    #[allow(clippy::too_many_arguments)]
    pub fn synthetic(
        pressure_pa: Vec<f64>,
        sample_rate_hz: u32,
        clock_identity: ContentHash,
        observer: AcousticObserver,
        bandwidth_hz: f64,
        uncertainty_pa: f64,
        source_identity: ContentHash,
    ) -> Result<Self, CinematicStemError> {
        Self::try_new(
            pressure_pa,
            sample_rate_hz,
            clock_identity,
            observer,
            bandwidth_hz,
            uncertainty_pa,
            AcousticAuthority::Synthetic,
            source_identity,
            PressureRepresentation::TimeDomainPascals,
        )
    }

    /// Admit estimated time-domain pressure samples in Pa.
    ///
    /// This factory preserves the estimated tier for declared-but-uncalibrated
    /// sources. It cannot construct calibrated authority or admit a phasor
    /// representation.
    #[allow(clippy::too_many_arguments)]
    pub fn estimated(
        pressure_pa: Vec<f64>,
        sample_rate_hz: u32,
        clock_identity: ContentHash,
        observer: AcousticObserver,
        bandwidth_hz: f64,
        uncertainty_pa: f64,
        source_identity: ContentHash,
    ) -> Result<Self, CinematicStemError> {
        Self::try_new(
            pressure_pa,
            sample_rate_hz,
            clock_identity,
            observer,
            bandwidth_hz,
            uncertainty_pa,
            AcousticAuthority::Estimated,
            source_identity,
            PressureRepresentation::TimeDomainPascals,
        )
    }

    /// Shared admission implementation for fixed-authority producer paths.
    /// Public callers cannot select an authority or representation tag. No
    /// public calibrated upstream producer adapter is currently exposed from
    /// this crate.
    ///
    /// The input is accepted only when it names time-domain pressure in Pa;
    /// peak/RMS phasors require their own explicit inverse transform first.
    #[allow(clippy::too_many_arguments)]
    fn try_new(
        pressure_pa: Vec<f64>,
        sample_rate_hz: u32,
        clock_identity: ContentHash,
        observer: AcousticObserver,
        bandwidth_hz: f64,
        uncertainty_pa: f64,
        authority: AcousticAuthority,
        source_identity: ContentHash,
        representation: PressureRepresentation,
    ) -> Result<Self, CinematicStemError> {
        if representation != PressureRepresentation::TimeDomainPascals {
            return Err(CinematicStemError::UnsupportedRepresentation {
                observed: representation,
            });
        }
        if pressure_pa.is_empty() {
            return Err(CinematicStemError::EmptySamples);
        }
        if pressure_pa.iter().any(|sample| !sample.is_finite()) {
            return Err(CinematicStemError::NonFinite {
                field: "pressure samples",
            });
        }
        if sample_rate_hz == 0 {
            return Err(CinematicStemError::InvalidValue {
                field: "sample rate",
            });
        }
        let nyquist_hz = f64::from(sample_rate_hz) * 0.5;
        if !bandwidth_hz.is_finite() || !(bandwidth_hz > 0.0 && bandwidth_hz <= nyquist_hz) {
            return Err(CinematicStemError::InvalidValue { field: "bandwidth" });
        }
        if !uncertainty_pa.is_finite() || uncertainty_pa < 0.0 {
            return Err(CinematicStemError::InvalidValue {
                field: "uncertainty",
            });
        }
        if clock_identity == ContentHash([0; 32]) {
            return Err(CinematicStemError::MissingIdentity { field: "clock" });
        }
        if source_identity == ContentHash([0; 32]) {
            return Err(CinematicStemError::MissingIdentity { field: "source" });
        }
        let pressure_identity = pressure_identity(sample_rate_hz, &pressure_pa);
        Ok(Self {
            pressure_pa,
            sample_rate_hz,
            clock_identity,
            observer,
            bandwidth_hz: canonical_zero(bandwidth_hz),
            uncertainty_pa: canonical_zero(uncertainty_pa),
            authority,
            source_identity,
            pressure_identity,
        })
    }

    /// ```compile_fail
    /// use fs_couple::cinematic_stem::AcousticPressureArtifact;
    ///
    /// // Raw samples may enter only through an admitted producer path.
    /// let _forged = AcousticPressureArtifact::try_new;
    /// ```

    /// Consume this artifact into one compatible cinematic stem without
    /// resampling, gain staging, peak normalization, or value rewriting.
    pub fn into_cinematic_stem(
        self,
        request: &CinematicStemRequest,
    ) -> Result<CinematicPressureStem, CinematicStemError> {
        if self.clock_identity != request.clock_identity {
            return Err(CinematicStemError::ClockMismatch);
        }
        if self.observer.identity() != request.observer_identity {
            return Err(CinematicStemError::ObserverMismatch);
        }
        if !request.minimum_authority.is_satisfied_by(self.authority) {
            return Err(CinematicStemError::InsufficientAuthority {
                observed: self.authority,
                required: request.minimum_authority,
            });
        }
        let identity = stem_identity(
            self.sample_rate_hz,
            self.clock_identity,
            self.observer.identity(),
            self.bandwidth_hz,
            self.uncertainty_pa,
            self.authority,
            self.source_identity,
            self.pressure_identity,
        );
        Ok(CinematicPressureStem {
            pressure_pa: self.pressure_pa,
            sample_rate_hz: self.sample_rate_hz,
            clock_identity: self.clock_identity,
            observer: self.observer,
            bandwidth_hz: self.bandwidth_hz,
            uncertainty_pa: self.uncertainty_pa,
            authority: self.authority,
            source_identity: self.source_identity,
            pressure_identity: self.pressure_identity,
            identity,
        })
    }
}

/// Compatibility requirements imposed by a cinematic consumer.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct CinematicStemRequest {
    clock_identity: ContentHash,
    observer_identity: ContentHash,
    minimum_authority: CinematicAuthorityRequirement,
}

impl CinematicStemRequest {
    /// Bind a consumer to one clock, observer, and minimum authority tier.
    pub fn try_new(
        clock_identity: ContentHash,
        observer_identity: ContentHash,
        minimum_authority: CinematicAuthorityRequirement,
    ) -> Result<Self, CinematicStemError> {
        if clock_identity == ContentHash([0; 32]) {
            return Err(CinematicStemError::MissingIdentity { field: "clock" });
        }
        if observer_identity == ContentHash([0; 32]) {
            return Err(CinematicStemError::MissingIdentity { field: "observer" });
        }
        Ok(Self {
            clock_identity,
            observer_identity,
            minimum_authority,
        })
    }
}

/// A cinematic-ready, still-physical pressure stem. It deliberately does not
/// encode PCM or apply room/mastering transforms.
#[derive(Clone, Debug, PartialEq)]
pub struct CinematicPressureStem {
    pressure_pa: Vec<f64>,
    sample_rate_hz: u32,
    clock_identity: ContentHash,
    observer: AcousticObserver,
    bandwidth_hz: f64,
    uncertainty_pa: f64,
    authority: AcousticAuthority,
    source_identity: ContentHash,
    pressure_identity: ContentHash,
    identity: ContentHash,
}

impl CinematicPressureStem {
    /// Unmodified time-domain Pa samples.
    #[must_use]
    pub fn pressure_pa(&self) -> &[f64] {
        &self.pressure_pa
    }

    /// Sample rate in Hz.
    #[must_use]
    pub const fn sample_rate_hz(&self) -> u32 {
        self.sample_rate_hz
    }

    /// Bound acoustic clock identity.
    #[must_use]
    pub const fn clock_identity(&self) -> ContentHash {
        self.clock_identity
    }

    /// Bound observer geometry.
    #[must_use]
    pub const fn observer(&self) -> &AcousticObserver {
        &self.observer
    }

    /// Declared maximum represented frequency in Hz.
    #[must_use]
    pub const fn bandwidth_hz(&self) -> f64 {
        self.bandwidth_hz
    }

    /// Declared amplitude uncertainty in Pa.
    #[must_use]
    pub const fn uncertainty_pa(&self) -> f64 {
        self.uncertainty_pa
    }

    /// Authority carried through the bridge without promotion.
    #[must_use]
    pub const fn authority(&self) -> AcousticAuthority {
        self.authority
    }

    /// Upstream source identity.
    #[must_use]
    pub const fn source_identity(&self) -> ContentHash {
        self.source_identity
    }

    /// Content identity of the exact retained Pa sample sequence.
    #[must_use]
    pub const fn pressure_identity(&self) -> ContentHash {
        self.pressure_identity
    }

    /// Identity of the bound cinematic-stem metadata.
    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }
}

fn stem_identity(
    sample_rate_hz: u32,
    clock_identity: ContentHash,
    observer_identity: ContentHash,
    bandwidth_hz: f64,
    uncertainty_pa: f64,
    authority: AcousticAuthority,
    source_identity: ContentHash,
    pressure_identity: ContentHash,
) -> ContentHash {
    let mut payload = Vec::with_capacity(4 + 4 + 32 * 4 + 8 + 8 + 1);
    payload.extend_from_slice(&CINEMATIC_STEM_SEMANTICS_VERSION.to_le_bytes());
    payload.extend_from_slice(&sample_rate_hz.to_le_bytes());
    payload.extend_from_slice(clock_identity.as_bytes());
    payload.extend_from_slice(observer_identity.as_bytes());
    payload.extend_from_slice(&bandwidth_hz.to_bits().to_le_bytes());
    payload.extend_from_slice(&uncertainty_pa.to_bits().to_le_bytes());
    payload.push(authority.identity_tag());
    payload.extend_from_slice(source_identity.as_bytes());
    payload.extend_from_slice(pressure_identity.as_bytes());
    hash_domain(STEM_IDENTITY_DOMAIN, &payload)
}

fn pressure_identity(sample_rate_hz: u32, samples_pa: &[f64]) -> ContentHash {
    let mut hasher = DomainHasher::new(PRESSURE_PAYLOAD_IDENTITY_DOMAIN);
    hasher.update(&CINEMATIC_STEM_SEMANTICS_VERSION.to_le_bytes());
    hasher.update(&sample_rate_hz.to_le_bytes());
    for sample in samples_pa {
        hasher.update(&sample.to_bits().to_le_bytes());
    }
    hasher.finalize()
}

fn canonical_zero(value: f64) -> f64 {
    if value == 0.0 { 0.0 } else { value }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn identity(label: &[u8]) -> ContentHash {
        hash_domain("org.frankensim.fs-couple.cinematic-stem.test", label)
    }

    fn observer() -> AcousticObserver {
        AcousticObserver::try_new([1.0, -0.0, 0.5], identity(b"frame")).expect("observer")
    }

    fn synthetic_artifact() -> AcousticPressureArtifact {
        AcousticPressureArtifact::synthetic(
            vec![-0.0, 0.25, -0.5, 0.125],
            48_000,
            identity(b"clock"),
            observer(),
            20_000.0,
            0.01,
            identity(b"source"),
        )
        .expect("synthetic fixture")
    }

    fn estimated_artifact() -> AcousticPressureArtifact {
        AcousticPressureArtifact::estimated(
            vec![-0.0, 0.25, -0.5, 0.125],
            48_000,
            identity(b"clock"),
            observer(),
            20_000.0,
            0.01,
            identity(b"estimated-source"),
        )
        .expect("estimated fixture")
    }

    #[test]
    fn time_domain_pascal_artifact_crosses_preview_without_rewriting_samples() {
        let artifact = synthetic_artifact();
        let request = CinematicStemRequest::try_new(
            identity(b"clock"),
            observer().identity(),
            CinematicAuthorityRequirement::Synthetic,
        )
        .expect("request");
        let stem = artifact.into_cinematic_stem(&request).expect("bridge");
        assert_eq!(
            stem.pressure_pa()
                .iter()
                .map(|sample| sample.to_bits())
                .collect::<Vec<_>>(),
            [
                (-0.0_f64).to_bits(),
                0.25_f64.to_bits(),
                (-0.5_f64).to_bits(),
                0.125_f64.to_bits()
            ],
            "the bridge must not normalize, gain-stage, or resample pressure"
        );
        assert_eq!(stem.sample_rate_hz(), 48_000);
        assert_eq!(stem.authority(), AcousticAuthority::Synthetic);
        assert_eq!(stem.observer().position_m()[1].to_bits(), 0.0_f64.to_bits());
    }

    #[test]
    fn authority_clock_observer_and_phasor_promotion_attacks_refuse() {
        let request = CinematicStemRequest::try_new(
            identity(b"clock"),
            observer().identity(),
            CinematicAuthorityRequirement::Calibrated,
        )
        .expect("calibrated request");
        assert!(matches!(
            synthetic_artifact().into_cinematic_stem(&request),
            Err(CinematicStemError::InsufficientAuthority {
                observed: AcousticAuthority::Synthetic,
                required: CinematicAuthorityRequirement::Calibrated,
            })
        ));
        let calibrated = AcousticPressureArtifact::try_new(
            vec![-0.0, 0.25, -0.5, 0.125],
            48_000,
            identity(b"clock"),
            observer(),
            20_000.0,
            0.01,
            AcousticAuthority::Calibrated(CalibratedAuthority::INTERNAL),
            identity(b"calibrated-source"),
            PressureRepresentation::TimeDomainPascals,
        )
        .expect("calibrated artifact")
        .into_cinematic_stem(&request)
        .expect("calibrated source crosses its requested tier");
        assert!(matches!(
            calibrated.authority(),
            AcousticAuthority::Calibrated(_)
        ));
        let wrong_clock = CinematicStemRequest::try_new(
            identity(b"other-clock"),
            observer().identity(),
            CinematicAuthorityRequirement::Synthetic,
        )
        .expect("request");
        assert_eq!(
            synthetic_artifact().into_cinematic_stem(&wrong_clock),
            Err(CinematicStemError::ClockMismatch)
        );
        let other_observer =
            AcousticObserver::try_new([2.0, 0.0, 0.5], identity(b"frame")).expect("other observer");
        let wrong_observer = CinematicStemRequest::try_new(
            identity(b"clock"),
            other_observer.identity(),
            CinematicAuthorityRequirement::Synthetic,
        )
        .expect("request");
        assert_eq!(
            synthetic_artifact().into_cinematic_stem(&wrong_observer),
            Err(CinematicStemError::ObserverMismatch)
        );
        for representation in [
            PressureRepresentation::FrequencyPhasorPeak,
            PressureRepresentation::FrequencyPhasorRms,
        ] {
            assert_eq!(
                AcousticPressureArtifact::try_new(
                    vec![1.0],
                    48_000,
                    identity(b"clock"),
                    observer(),
                    20_000.0,
                    0.01,
                    AcousticAuthority::Calibrated(CalibratedAuthority::INTERNAL),
                    identity(b"source"),
                    representation,
                ),
                Err(CinematicStemError::UnsupportedRepresentation {
                    observed: representation,
                })
            );
        }
    }

    #[test]
    fn estimated_artifact_crosses_estimated_tier_without_calibration_promotion() {
        let estimated_request = CinematicStemRequest::try_new(
            identity(b"clock"),
            observer().identity(),
            CinematicAuthorityRequirement::Estimated,
        )
        .expect("estimated request");
        let stem = estimated_artifact()
            .into_cinematic_stem(&estimated_request)
            .expect("estimated source crosses its requested tier");
        assert_eq!(stem.authority(), AcousticAuthority::Estimated);

        let calibrated_request = CinematicStemRequest::try_new(
            identity(b"clock"),
            observer().identity(),
            CinematicAuthorityRequirement::Calibrated,
        )
        .expect("calibrated request");
        assert!(matches!(
            estimated_artifact().into_cinematic_stem(&calibrated_request),
            Err(CinematicStemError::InsufficientAuthority {
                observed: AcousticAuthority::Estimated,
                required: CinematicAuthorityRequirement::Calibrated,
            })
        ));
    }

    #[test]
    fn metadata_identity_is_deterministic_and_source_sensitive() {
        let request = CinematicStemRequest::try_new(
            identity(b"clock"),
            observer().identity(),
            CinematicAuthorityRequirement::Synthetic,
        )
        .expect("request");
        let first = synthetic_artifact()
            .into_cinematic_stem(&request)
            .expect("first");
        let second = synthetic_artifact()
            .into_cinematic_stem(&request)
            .expect("second");
        assert_eq!(first.identity(), second.identity());

        let changed_source = AcousticPressureArtifact::try_new(
            vec![-0.0, 0.25, -0.5, 0.125],
            48_000,
            identity(b"clock"),
            observer(),
            20_000.0,
            0.01,
            AcousticAuthority::Synthetic,
            identity(b"another-source"),
            PressureRepresentation::TimeDomainPascals,
        )
        .expect("changed source")
        .into_cinematic_stem(&request)
        .expect("bridge");
        assert_ne!(first.identity(), changed_source.identity());

        let changed_samples = AcousticPressureArtifact::try_new(
            vec![-0.0, 0.25, -0.25, 0.125],
            48_000,
            identity(b"clock"),
            observer(),
            20_000.0,
            0.01,
            AcousticAuthority::Synthetic,
            identity(b"source"),
            PressureRepresentation::TimeDomainPascals,
        )
        .expect("changed samples")
        .into_cinematic_stem(&request)
        .expect("bridge");
        assert_ne!(
            first.pressure_identity(),
            changed_samples.pressure_identity()
        );
        assert_ne!(first.identity(), changed_samples.identity());
    }
}
