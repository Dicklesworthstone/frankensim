//! Common conservative-potential view of the existing source-bound owners.
//! No power-law fit is made to a finite-gap table, and no authority is upgraded.
use super::{
    FiniteGapResponseCurve, FiniteGapResponseCurveError, NormalPatchError,
    adhesive::{JkrError, SourceJkrSphereContact},
    dmt::{DmtError, SourceDmtSphereContact},
    maugis::{
        MaugisError, continuation::MaugisContinuationError, mechanical::MaugisAttachedPotential,
    },
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
    /// Conservative DMT contact and its explicitly selected finite-range attraction.
    DmtSphere(Box<SourceDmtSphereContact>),
    /// Explicit attached Maugis branch. Its limit refuses; no snap or detached
    /// branch is silently substituted by a mechanical solver.
    MaugisAttached(Box<MaugisAttachedPotential>),
    /// Explicit stable attached JKR branch, including negative approach and
    /// tensile force. Crossing its fold refuses instead of choosing a snap.
    JkrAttached(Box<SourceJkrSphereContact>),
}

/// Original source-owner refusal, preserved across the potential view.
#[derive(Debug, Clone, PartialEq)]
pub enum NormalLawSourceError {
    /// Source-bound analytic point law refused.
    Point(NormalPatchError),
    /// Finite-gap curve refused; extrapolation is never substituted.
    FiniteGap(FiniteGapResponseCurveError),
    /// Original adhesive owner refused the approach or smooth tangent.
    Dmt(DmtError),
    /// Original finite-zone branch, quadrature or independent energy refusal.
    Maugis(MaugisContinuationError),
    /// Original short-range adhesive source, branch or tangent refused.
    Jkr(JkrError),
}

impl core::fmt::Display for NormalLawSourceError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Point(e) => core::fmt::Display::fmt(e, f),
            Self::FiniteGap(e) => core::fmt::Display::fmt(e, f),
            Self::Dmt(e) => core::fmt::Display::fmt(e, f),
            Self::Maugis(e) => core::fmt::Display::fmt(e, f),
            Self::Jkr(e) => core::fmt::Display::fmt(e, f),
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

impl From<SourceDmtSphereContact> for NormalLawSource {
    fn from(source: SourceDmtSphereContact) -> Self {
        Self::DmtSphere(Box::new(source))
    }
}
impl From<MaugisAttachedPotential> for NormalLawSource {
    fn from(source: MaugisAttachedPotential) -> Self {
        Self::MaugisAttached(Box::new(source))
    }
}
impl From<SourceJkrSphereContact> for NormalLawSource {
    /// Select the attached branch explicitly; no history or snap is inferred.
    fn from(source: SourceJkrSphereContact) -> Self {
        Self::JkrAttached(Box::new(source))
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
            Self::DmtSphere(s) => s.identity(),
            Self::MaugisAttached(s) => s.identity(),
            Self::JkrAttached(s) => s.identity(),
        }
    }

    /// Whether signed approach belongs to this owner's declared potential.
    /// Its branch/domain limits still apply; negative approach is not clipped.
    /// Unilateral Hertz/finite-gap owners retain their existing clipping policy.
    #[must_use]
    pub const fn supports_signed_approach(&self) -> bool {
        matches!(
            self,
            Self::DmtSphere(_) | Self::MaugisAttached(_) | Self::JkrAttached(_)
        )
    }

    /// Continuous energy and force, including adhesive touch/range corners.
    /// A discrete-gradient step needs these values, not a nonexistent smooth
    /// tangent at a piecewise-law corner.
    pub fn energy_and_force(&self, delta_m: f64) -> Result<(f64, f64), NormalLawSourceError> {
        if let Self::MaugisAttached(source) = self {
            let r = source
                .evaluate(delta_m)
                .map_err(NormalLawSourceError::Maugis)?;
            Ok((r.potential_j(), r.normal_force_n()))
        } else if let Self::DmtSphere(source) = self {
            let r = source
                .evaluate(delta_m)
                .map_err(NormalLawSourceError::Dmt)?;
            Ok((r.reversible_energy_j, r.normal_force_n))
        } else if let Self::JkrAttached(source) = self {
            let r = source
                .evaluate_attached(delta_m)
                .map_err(NormalLawSourceError::Jkr)?;
            Ok((r.reversible_energy_j, r.normal_force_n))
        } else {
            self.potential(delta_m).map(|(u, f, _)| (u, f))
        }
    }

    /// Cancellation-stable scalar Gonzalez numerator when the original owner
    /// supplies it. `None` preserves the native endpoint-subtraction fallback.
    pub fn midpoint_energy_remainder(
        &self,
        from_m: f64,
        to_m: f64,
        midpoint_m: f64,
    ) -> Result<Option<f64>, NormalLawSourceError> {
        match self {
            Self::MaugisAttached(source) => source
                .midpoint_energy_remainder(from_m, to_m, midpoint_m)
                .map(Some)
                .map_err(NormalLawSourceError::Maugis),
            Self::DmtSphere(source) => source
                .midpoint_energy_remainder(from_m, to_m, midpoint_m)
                .map(Some)
                .map_err(NormalLawSourceError::Dmt),
            _ => Ok(None),
        }
    }

    /// Mean effort on the original straight approach path when the owner
    /// supplies one. This permits native componentwise work and reciprocal
    /// force projection without a global energy correction on other bodies.
    /// Unilateral owners include the inactive part of a SIGNED affine path:
    /// clipping both endpoints before this call would overcount onset work.
    /// Point-contact loss and actual rate/duration admission remain separate.
    pub fn mean_force_between(
        &self,
        from_m: f64,
        to_m: f64,
    ) -> Result<Option<f64>, NormalLawSourceError> {
        match self {
            Self::JkrAttached(source) => source
                .attached_path_average(from_m, to_m)
                .map(Some)
                .map_err(NormalLawSourceError::Jkr),
            Self::MaugisAttached(source) => source
                .path(from_m, to_m)
                .map(|p| Some(p.mean_force_n))
                .map_err(NormalLawSourceError::Maugis),
            Self::DmtSphere(source) => {
                let force = if from_m == to_m {
                    source
                        .evaluate(from_m)
                        .map_err(NormalLawSourceError::Dmt)?
                        .normal_force_n
                } else {
                    source
                        .reversible_work_between(from_m, to_m)
                        .map_err(NormalLawSourceError::Dmt)?
                        / (to_m - from_m)
                };
                if !force.is_finite() {
                    return Err(NormalLawSourceError::Dmt(DmtError::Unrepresentable {
                        field: "mechanical mean effort",
                    }));
                }
                Ok(Some(force))
            }
            Self::HertzPoint(source) => {
                let (from, to, active_fraction) = unilateral_path(from_m, to_m).ok_or(
                    NormalLawSourceError::Point(NormalPatchError::InvalidInput {
                        field: "finite signed normal approach path",
                    }),
                )?;
                let active_force = source
                    .conservative_path_average(from, to)
                    .map_err(NormalLawSourceError::Point)?
                    .mean_force_n;
                let force = active_force * active_fraction;
                if !force.is_finite()
                    || (active_force > 0.0 && active_fraction > 0.0 && force == 0.0)
                {
                    return Err(NormalLawSourceError::Point(NormalPatchError::Overflow {
                        field: "mechanical mean effort",
                    }));
                }
                Ok(Some(force))
            }
            Self::FiniteGap(source) => {
                let (from, to, active_fraction) = unilateral_path(from_m, to_m).ok_or(
                    NormalLawSourceError::FiniteGap(FiniteGapResponseCurveError::InvalidInput {
                        field: "finite signed normal approach path",
                    }),
                )?;
                let active_force = source
                    .mean_force_between(from, to)
                    .map_err(NormalLawSourceError::FiniteGap)?;
                let force = active_force * active_fraction;
                if !force.is_finite()
                    || (active_force > 0.0 && active_fraction > 0.0 && force == 0.0)
                {
                    return Err(NormalLawSourceError::FiniteGap(
                        FiniteGapResponseCurveError::InvalidInput {
                            field: "mechanical mean effort",
                        },
                    ));
                }
                Ok(Some(force))
            }
        }
    }

    /// Conservative `(U [J], F [N], dF/d(delta) [N/m])` at an approach [m].
    /// The point owner is evaluated at zero rate because internal loss is a
    /// separate port, not stored energy. The table owner integrates its own
    /// interpolation segment, so its derivative and energy cannot drift apart.
    pub fn potential(&self, delta_m: f64) -> Result<(f64, f64, f64), NormalLawSourceError> {
        match self {
            Self::JkrAttached(source) => {
                let r = source
                    .evaluate_attached(delta_m)
                    .map_err(NormalLawSourceError::Jkr)?;
                let tangent = r
                    .tangent_n_per_m
                    .ok_or(NormalLawSourceError::Jkr(JkrError::SingularTangent))?;
                Ok((r.reversible_energy_j, r.normal_force_n, tangent))
            }
            Self::MaugisAttached(source) => {
                let r = source
                    .evaluate(delta_m)
                    .map_err(NormalLawSourceError::Maugis)?;
                let tangent = r.tangent_n_per_m().ok_or_else(|| {
                    NormalLawSourceError::Maugis(MaugisError::EndpointTangent.into())
                })?;
                Ok((r.potential_j(), r.normal_force_n(), tangent))
            }
            Self::DmtSphere(source) => {
                let r = source
                    .evaluate(delta_m)
                    .map_err(NormalLawSourceError::Dmt)?;
                let tangent = r
                    .tangent_n_per_m
                    .ok_or(NormalLawSourceError::Dmt(DmtError::CornerTangent))?;
                Ok((r.reversible_energy_j, r.normal_force_n, tangent))
            }
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

// Preserve the fraction of the ORIGINAL affine coordinate spent in contact.
// The storage source itself still receives only its admitted indentation.
fn unilateral_path(from_m: f64, to_m: f64) -> Option<(f64, f64, f64)> {
    if !from_m.is_finite() || !to_m.is_finite() {
        return None;
    }
    let lo = from_m.min(to_m);
    let hi = from_m.max(to_m);
    let fraction = if hi <= 0.0 {
        0.0
    } else if lo >= 0.0 {
        1.0
    } else {
        let span = hi - lo;
        if !span.is_finite() {
            return None;
        }
        hi / span
    };
    (fraction.is_finite() && (hi <= 0.0 || fraction > 0.0)).then_some((
        from_m.max(0.0),
        to_m.max(0.0),
        fraction,
    ))
}
