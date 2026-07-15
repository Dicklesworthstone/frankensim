//! Codimensional thickness gating (bead rjnd, E1 query upgrades,
//! part 5).
//!
//! Shells and rods are represented by a midsurface/centerline plus a
//! certified thickness radius. Contact between two such bodies is a
//! statement about the EFFECTIVE gap
//! `dist(midsurfaces) - (t_a + t_b)`, so the whole primitive is an
//! outward-rounded interval subtraction plus a fail-closed verdict:
//!
//! - `ProvenClear`: the effective gap's certified lower bound is
//!   positive — the offset bodies are disjoint.
//! - `ProvenContact`: the certified upper bound is negative — the
//!   offset bodies must interpenetrate.
//! - `Unresolved`: the bracket straddles zero; NOTHING is claimed.
//!
//! The caller owns the midsurface-distance certificate. Soundness is
//! directional: `ProvenClear` needs `lo` to be a true lower bound on
//! the midsurface distance (a convex-hull separation lower bound
//! qualifies, since hulls only shrink distances), while
//! `ProvenContact` needs `hi` to be realized between actual
//! midsurface points (support-map witnesses qualify for our
//! primitives). [`codim_gap_from_separation`] composes with
//! [`crate::ConvexSeparation`] under exactly that reading.

use crate::{ContactInflation, ConvexSeparation, QueryError};

/// A certified codimensional thickness radius (half-thickness for a
/// shell, cross-section radius for a rod), measured outward from the
/// midsurface/centerline.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CodimThickness {
    radius: f64,
}

impl CodimThickness {
    /// Validated construction. Zero is legal (a degenerate midsurface
    /// body: the effective gap IS the midsurface distance).
    ///
    /// # Errors
    /// [`QueryError::CodimInvalidThickness`] for non-finite or
    /// negative radii.
    pub fn new(radius: f64) -> Result<CodimThickness, QueryError> {
        if radius.is_finite() && radius >= 0.0 {
            Ok(CodimThickness { radius })
        } else {
            Err(QueryError::CodimInvalidThickness {
                thickness_bits: radius.to_bits(),
            })
        }
    }

    /// The certified radius.
    #[must_use]
    pub fn radius(&self) -> f64 {
        self.radius
    }
}

/// The fail-closed contact verdict for one codimensional pair.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CodimVerdict {
    /// The offset bodies are proven disjoint (`effective lo > 0`).
    ProvenClear,
    /// The offset bodies are proven to interpenetrate
    /// (`effective hi < 0`).
    ProvenContact,
    /// The bracket straddles zero: no claim either way.
    Unresolved,
}

/// A certified effective-gap bracket with its verdict.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CodimGap {
    /// Certified lower bound on `dist(midsurfaces) - (t_a + t_b)`.
    pub lo: f64,
    /// Certified upper bound on the same quantity.
    pub hi: f64,
    /// The fail-closed verdict derived from the bracket alone.
    pub verdict: CodimVerdict,
}

/// Effective gap from a caller-certified midsurface distance
/// enclosure `[distance_lo, distance_hi]` and two thickness radii.
///
/// The subtraction is rounded against each claim: the combined
/// thickness is rounded UP before the lower bound and DOWN before the
/// upper bound, so the returned bracket always contains the true
/// effective gap whenever the inputs hold.
///
/// # Errors
/// [`QueryError::CodimInvalidDistance`] for a non-finite, inverted,
/// or negative distance enclosure (set distances are nonnegative; a
/// negative claim is a caller bug, not a contact proof).
pub fn codim_gap(
    distance_lo: f64,
    distance_hi: f64,
    a: CodimThickness,
    b: CodimThickness,
) -> Result<CodimGap, QueryError> {
    if !(distance_lo.is_finite() && distance_hi.is_finite())
        || distance_lo > distance_hi
        || distance_lo < 0.0
    {
        return Err(QueryError::CodimInvalidDistance {
            reason: "midsurface distance enclosure must be finite, ordered, and nonnegative",
        });
    }
    let thickness_up = (a.radius + b.radius).next_up();
    let thickness_dn = (a.radius + b.radius).next_down();
    let lo = (distance_lo - thickness_up).next_down();
    let hi = (distance_hi - thickness_dn).next_up();
    let verdict = if lo > 0.0 {
        CodimVerdict::ProvenClear
    } else if hi < 0.0 {
        CodimVerdict::ProvenContact
    } else {
        CodimVerdict::Unresolved
    };
    Ok(CodimGap { lo, hi, verdict })
}

/// Effective codimensional gap with representation/motion uncertainty on
/// both midsurface inputs.
///
/// The ordinary distance enclosure is validated first. The composed
/// proof-bearing radius then lowers its nonnegative lower endpoint and raises
/// its upper endpoint before thickness is subtracted. Exact-zero inflations
/// return the bit-identical ordinary result.
///
/// # Errors
/// The refusals from [`codim_gap`], plus
/// [`QueryError::InvalidContactInflation`] if the two radii cannot be composed
/// or applied without finite outward-rounded endpoints.
pub fn codim_gap_with_inflation(
    distance_lo: f64,
    distance_hi: f64,
    a: CodimThickness,
    b: CodimThickness,
    inflation_a: ContactInflation,
    inflation_b: ContactInflation,
) -> Result<CodimGap, QueryError> {
    let nominal = codim_gap(distance_lo, distance_hi, a, b)?;
    let inflation = inflation_a.compose(inflation_b)?;
    if inflation.radius() == 0.0 {
        return Ok(nominal);
    }
    let widened_lo = inflation.deflate_nonnegative(distance_lo)?;
    let widened_hi = inflation.inflate_upper(distance_hi)?;
    codim_gap(widened_lo, widened_hi, a, b)
}

/// Compose with a convex midsurface-hull separation.
///
/// `separation.lo` is always a sound lower bound for the midsurface
/// distance (hulls only shrink distances), so `ProvenClear` is
/// unconditional. `ProvenContact` additionally requires the
/// separation's witness points to lie IN the midsurface sets (true
/// for the shipped support maps; a caller obligation for custom
/// ones) — when `witnesses_in_sets` is false, a contact-side bracket
/// is DOWNGRADED to [`CodimVerdict::Unresolved`] instead of claiming
/// contact from hull evidence.
///
/// # Errors
/// [`QueryError::CodimInvalidDistance`] as for [`codim_gap`].
pub fn codim_gap_from_separation(
    separation: &ConvexSeparation,
    a: CodimThickness,
    b: CodimThickness,
    witnesses_in_sets: bool,
) -> Result<CodimGap, QueryError> {
    let mut gap = codim_gap(separation.lo, separation.hi, a, b)?;
    if gap.verdict == CodimVerdict::ProvenContact && !witnesses_in_sets {
        gap.verdict = CodimVerdict::Unresolved;
    }
    Ok(gap)
}
