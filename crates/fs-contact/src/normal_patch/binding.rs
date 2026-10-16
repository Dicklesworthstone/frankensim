//! Common conservative-potential view of the existing source-bound owners.
//! No power-law fit is made to a finite-gap table, and no authority is upgraded.
use super::{
    FiniteGapResponseCurve, FiniteGapResponseCurveError, NormalPatchError,
    source::SourcePointContact,
};
use fs_blake3::ContentHash;

/// Frozen physical normal law; the point coefficient is material/card derived,
/// while the finite-gap curve retains its explicit interpolation estimate.
#[derive(Debug, Clone, PartialEq)]
pub enum NormalLawSource {
    /// Hertz or Hunt--Crossley point owner, including its source envelope.
    HertzPoint(Box<SourcePointContact>),
    /// Whole-footprint force curve from retained finite-gap half-space solves.
    /// Its force and energy already integrate the physical surface-cell areas.
    FiniteGap(Box<FiniteGapResponseCurve>),
}

/// Original source-owner refusal, preserved across the potential view.
#[derive(Debug, Clone, PartialEq)]
pub enum NormalLawSourceError {
    /// Source-bound analytic point law refused.
    Point(NormalPatchError),
    /// Finite-gap curve refused; extrapolation is never substituted.
    FiniteGap(FiniteGapResponseCurveError),
}

impl core::fmt::Display for NormalLawSourceError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Point(e) => core::fmt::Display::fmt(e, f),
            Self::FiniteGap(e) => core::fmt::Display::fmt(e, f),
        }
    }
}
impl std::error::Error for NormalLawSourceError {}

impl From<SourcePointContact> for NormalLawSource {
    fn from(source: SourcePointContact) -> Self {
        Self::HertzPoint(Box::new(source))
    }
}
impl From<FiniteGapResponseCurve> for NormalLawSource {
    fn from(source: FiniteGapResponseCurve) -> Self {
        Self::FiniteGap(Box::new(source))
    }
}

impl NormalLawSource {
    /// Original point limiting-response receipt (including geometry/envelope),
    /// or complete curve identity, without minting a new receipt.
    #[must_use]
    pub fn identity(&self) -> ContentHash {
        match self {
            Self::HertzPoint(s) => s.limiting_response().receipt_id,
            Self::FiniteGap(s) => s.identity,
        }
    }

    /// Conservative `(U [J], F [N], dF/d(delta) [N/m])` at an approach [m].
    /// The point owner is evaluated at zero rate because internal loss is a
    /// separate port, not stored energy. The table owner integrates its own
    /// interpolation segment, so its derivative and energy cannot drift apart.
    pub fn potential(&self, delta_m: f64) -> Result<(f64, f64, f64), NormalLawSourceError> {
        match self {
            Self::HertzPoint(source) => {
                let r = source
                    .evaluate(
                        "normal-law-source/potential",
                        delta_m,
                        0.0,
                        source.envelope().maximum_step_s,
                    )
                    .map_err(NormalLawSourceError::Point)?;
                Ok((r.reversible_energy_j, r.normal_force_n, r.tangent_n_per_m))
            }
            Self::FiniteGap(source) => {
                let r = source
                    .evaluate(delta_m)
                    .map_err(NormalLawSourceError::FiniteGap)?;
                Ok((
                    r.reversible_energy_j,
                    r.normal_force_n,
                    r.normal_tangent_n_per_m,
                ))
            }
        }
    }
}
