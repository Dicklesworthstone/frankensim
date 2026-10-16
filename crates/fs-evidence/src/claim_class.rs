//! Closed E09 claim-to-evidence taxonomy shared by governance and evidence.
//!
//! The table rows and capability mappings remain in `fs-govern`. These enums
//! live at the lower evidence layer so typed refusals can name the exact
//! doctrine class without introducing an upward dependency.

/// One claim family in the closed v1 certificate-regime doctrine.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ClaimClass {
    /// A root or event time over one declared finite domain.
    RootOrEventTime,
    /// Collision, exclusion, or reachability over a declared short horizon.
    ShortHorizonReachability,
    /// Conservation or balance over a declared discrete operator and boundary.
    ConservedQuantity,
    /// Local stability over a stated model, state, and parameter domain.
    LocalStability,
    /// A long-horizon mean load or other statistical observable.
    LongHorizonMeanLoad,
    /// A turbulent or broadband distributional/spectral observable.
    BroadbandSpectrum,
    /// Reliability over a declared population of duty cycles.
    DutyCycleReliability,
    /// One exact trajectory far beyond an admitted chaotic predictability horizon.
    ExactLongChaoticTrajectory,
}

impl ClaimClass {
    /// Every v1 claim family in canonical schema order.
    pub const ALL: [Self; 8] = [
        Self::RootOrEventTime,
        Self::ShortHorizonReachability,
        Self::ConservedQuantity,
        Self::LocalStability,
        Self::LongHorizonMeanLoad,
        Self::BroadbandSpectrum,
        Self::DutyCycleReliability,
        Self::ExactLongChaoticTrajectory,
    ];

    /// Stable machine code.
    #[must_use]
    pub const fn code(self) -> &'static str {
        match self {
            Self::RootOrEventTime => "root-or-event-time",
            Self::ShortHorizonReachability => "short-horizon-reachability",
            Self::ConservedQuantity => "conserved-quantity",
            Self::LocalStability => "local-stability",
            Self::LongHorizonMeanLoad => "long-horizon-mean-load",
            Self::BroadbandSpectrum => "broadband-spectrum",
            Self::DutyCycleReliability => "duty-cycle-reliability",
            Self::ExactLongChaoticTrajectory => "exact-long-chaotic-trajectory",
        }
    }

    /// Human-readable claim name.
    #[must_use]
    pub const fn title(self) -> &'static str {
        match self {
            Self::RootOrEventTime => "Root or event time",
            Self::ShortHorizonReachability => "Short-horizon collision or reachability",
            Self::ConservedQuantity => "Conserved quantity",
            Self::LocalStability => "Local stability",
            Self::LongHorizonMeanLoad => "Long-horizon mean load",
            Self::BroadbandSpectrum => "Turbulent or broadband spectrum",
            Self::DutyCycleReliability => "Reliability over duty cycles",
            Self::ExactLongChaoticTrajectory => "Exact long chaotic trajectory",
        }
    }

    /// Stable index into the closed schema-v1 doctrine table.
    #[must_use]
    pub const fn canonical_index(self) -> usize {
        match self {
            Self::RootOrEventTime => 0,
            Self::ShortHorizonReachability => 1,
            Self::ConservedQuantity => 2,
            Self::LocalStability => 3,
            Self::LongHorizonMeanLoad => 4,
            Self::BroadbandSpectrum => 5,
            Self::DutyCycleReliability => 6,
            Self::ExactLongChaoticTrajectory => 7,
        }
    }

    /// The only evidence regime admitted by doctrine v1 for this claim class.
    #[must_use]
    pub const fn required_evidence(self) -> EvidenceRegime {
        match self {
            Self::RootOrEventTime => EvidenceRegime::IntervalRootOrTaylorEnclosure,
            Self::ShortHorizonReachability => EvidenceRegime::ValidatedReachabilityTube,
            Self::ConservedQuantity => EvidenceRegime::DiscreteBalanceCertificate,
            Self::LocalStability => EvidenceRegime::SpectralOrLyapunovCertificate,
            Self::LongHorizonMeanLoad => EvidenceRegime::StatisticalObservableWithModelEvidence,
            Self::BroadbandSpectrum => EvidenceRegime::DistributionalSpectralValidation,
            Self::DutyCycleReliability => EvidenceRegime::SequentialRareEventStatistics,
            Self::ExactLongChaoticTrajectory => EvidenceRegime::NoUsefulBound,
        }
    }
}

/// Evidence object or explicit no-useful-bound outcome selected by doctrine.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum EvidenceRegime {
    /// Interval root isolation or a Taylor enclosure on a finite domain.
    IntervalRootOrTaylorEnclosure,
    /// A validated finite-horizon tube.
    ValidatedReachabilityTube,
    /// A discrete conservation or balance certificate.
    DiscreteBalanceCertificate,
    /// A spectral or Lyapunov certificate on a stated local domain.
    SpectralOrLyapunovCertificate,
    /// A statistical observable with sampling and model-form evidence.
    StatisticalObservableWithModelEvidence,
    /// Distributional and spectral validation against external observations.
    DistributionalSpectralValidation,
    /// Sequential or rare-event statistical evidence.
    SequentialRareEventStatistics,
    /// Honest refusal to claim a useful bound for the requested proposition.
    NoUsefulBound,
}

impl EvidenceRegime {
    /// Stable machine code.
    #[must_use]
    pub const fn code(self) -> &'static str {
        match self {
            Self::IntervalRootOrTaylorEnclosure => "interval-root-or-taylor-enclosure",
            Self::ValidatedReachabilityTube => "validated-reachability-tube",
            Self::DiscreteBalanceCertificate => "discrete-balance-certificate",
            Self::SpectralOrLyapunovCertificate => "spectral-or-lyapunov-certificate",
            Self::StatisticalObservableWithModelEvidence => {
                "statistical-observable-with-model-evidence"
            }
            Self::DistributionalSpectralValidation => "distributional-spectral-validation",
            Self::SequentialRareEventStatistics => "sequential-rare-event-statistics",
            Self::NoUsefulBound => "no-useful-bound",
        }
    }

    /// Human-readable evidence name.
    #[must_use]
    pub const fn title(self) -> &'static str {
        match self {
            Self::IntervalRootOrTaylorEnclosure => "Interval root isolation or Taylor enclosure",
            Self::ValidatedReachabilityTube => "Validated finite-horizon tube",
            Self::DiscreteBalanceCertificate => "Discrete conservation or balance certificate",
            Self::SpectralOrLyapunovCertificate => {
                "Spectral or Lyapunov certificate in a stated domain"
            }
            Self::StatisticalObservableWithModelEvidence => {
                "Statistical observable with sampling and model evidence"
            }
            Self::DistributionalSpectralValidation => "Distributional and spectral validation",
            Self::SequentialRareEventStatistics => "Sequential or rare-event statistical evidence",
            Self::NoUsefulBound => "NoUsefulBound",
        }
    }
}
