//! Sparse energy deposition from a caller-supplied radiation response.
//!
//! Coefficients describe the fraction of energy leaving each nuclear source
//! that actually deposits in each receiver, plus energy escaping the modeled
//! domain. They must describe a resolved transport problem for the declared
//! geometry, materials and emitted spectrum. They may be supplied independently
//! or prepared by the explicitly scoped planar photon [`transport`] owner.
//! Attenuation and particle counts are not equated with heating. No nuclear
//! data, dose model or automatic geometry update is supplied.

use std::collections::{BTreeMap, BTreeSet};
use std::sync::Arc;

/// Energy-dependent deposition on an explicit emitted-radiation group grid.
pub mod grouped;

/// Evaluated decay spectra reconciled to the original outward energy account.
pub mod decay;

/// Conservative planar photon transport that prepares energy-group responses.
pub mod transport;

/// Fixed admission tolerance for the sum of a supplied response row.
/// Coefficients are never renormalized; admitted defects are charged in joules.
pub const ROW_CLOSURE_TOLERANCE: f64 = 16.0 * f64::EPSILON;

/// One supplied deposited-energy fraction, not an interaction probability.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RadiationDeposit {
    /// Physical receiver identity, independent of any thermal representation.
    pub receiver_id: u64,
    /// Fraction of this source's emitted radiation energy thermalized here.
    pub fraction: f64,
}

/// Complete, explicitly closed response for one source spectrum and location.
#[derive(Clone, Debug, PartialEq)]
pub struct RadiationSourceResponse {
    /// Physical emitting source identity.
    pub source_id: u64,
    /// Unique receiving cells, possibly empty when all energy escapes.
    pub deposits: Vec<RadiationDeposit>,
    /// Fraction escaping the modeled domain after any transport interactions.
    pub escape_fraction: f64,
}

/// Allocation and application-work limits for an immutable sparse response.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct RadiationResponseLimits {
    /// Maximum number of source rows.
    pub max_sources: usize,
    /// Maximum number of distinct receiver identities across all rows.
    pub max_receivers: usize,
    /// Maximum number of explicitly authored source-to-receiver entries.
    pub max_transfers: usize,
}

/// Shared immutable response data. Cloning does not copy its sparse rows.
#[derive(Clone, Debug)]
pub struct RadiationResponse(Arc<ResponseData>);

#[derive(Clone, Debug)]
struct ResponseData {
    data_id: String,
    rows: Vec<RadiationSourceResponse>,
    receiver_ids: Vec<u64>,
    closure_defects: Vec<(Sum, f64)>,
}

/// Accepted response in joules, distinct from the original source accounts.
#[derive(Clone, Debug, PartialEq)]
pub struct RadiationDeposition {
    /// Radiation heat to credit exactly once to each physical receiver.
    pub receiver_energy_j: BTreeMap<u64, f64>,
    /// Total radiation energy supplied by the emitting sources.
    pub emitted_j: f64,
    /// Sum of the returned receiver energies.
    pub deposited_j: f64,
    /// Energy remaining outside the modeled receiving domain.
    pub escaped_j: f64,
    /// Signed `deposited + escaped - emitted` scalar-account residual.
    pub residual_j: f64,
    /// Noncancelling row-closure, product, reduction, scalar-export and final
    /// residual defects. This accounts for arithmetic, not transport-model error.
    pub roundoff_j: f64,
}

/// Refusal without partially accepted receiver credits or mutable history.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum RadiationError {
    /// Invalid coefficients, dataset identity, source coverage or energy input.
    InvalidInput(&'static str),
    /// Explicit sparse-work or energy-accounting bound exceeded.
    Budget(&'static str),
    /// A required positive energy product or aggregate is not representable.
    Arithmetic,
    /// Caller requested cooperative cancellation.
    Cancelled,
}

impl core::fmt::Display for RadiationError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::InvalidInput(detail) => write!(f, "invalid radiation response: {detail}"),
            Self::Budget(detail) => write!(f, "radiation response budget exceeded: {detail}"),
            Self::Arithmetic => f.write_str("radiation response arithmetic is not representable"),
            Self::Cancelled => f.write_str("radiation response cancelled"),
        }
    }
}
impl std::error::Error for RadiationError {}

impl RadiationResponse {
    /// Admit and canonically sort sparse rows without changing coefficients.
    /// The identity must name the supplied response and its applicability; it
    /// does not establish transport accuracy or validate physical source data.
    pub fn new(
        data_id: String,
        mut rows: Vec<RadiationSourceResponse>,
        limits: RadiationResponseLimits,
    ) -> Result<Self, RadiationError> {
        if data_id.trim().is_empty() {
            return Err(RadiationError::InvalidInput(
                "a response data identity is required",
            ));
        }
        if rows.len() > limits.max_sources {
            return Err(RadiationError::Budget("source count"));
        }
        let mut receivers = BTreeSet::new();
        let mut transfer_count = 0_usize;
        for row in &rows {
            transfer_count = transfer_count
                .checked_add(row.deposits.len())
                .ok_or(RadiationError::Budget("transfer count"))?;
            if transfer_count > limits.max_transfers {
                return Err(RadiationError::Budget("transfer count"));
            }
            for deposit in &row.deposits {
                receivers.insert(deposit.receiver_id);
                if receivers.len() > limits.max_receivers {
                    return Err(RadiationError::Budget("receiver count"));
                }
            }
        }
        rows.sort_by_key(|row| row.source_id);
        if rows
            .windows(2)
            .any(|pair| pair[0].source_id == pair[1].source_id)
        {
            return Err(RadiationError::InvalidInput("duplicate source identity"));
        }
        let mut closure_defects = Vec::with_capacity(rows.len());
        for row in &mut rows {
            if !fraction(row.escape_fraction) {
                return Err(RadiationError::InvalidInput(
                    "escape fraction must be in [0, 1]",
                ));
            }
            row.deposits.sort_by_key(|deposit| deposit.receiver_id);
            if row
                .deposits
                .windows(2)
                .any(|pair| pair[0].receiver_id == pair[1].receiver_id)
            {
                return Err(RadiationError::InvalidInput(
                    "duplicate receiver within source row",
                ));
            }
            let mut closure = Sum::default();
            let mut error = ErrorBudget::new(ROW_CLOSURE_TOLERANCE);
            for deposit in &row.deposits {
                if !fraction(deposit.fraction) {
                    return Err(RadiationError::InvalidInput(
                        "deposition fraction must be in [0, 1]",
                    ));
                }
                error.charge(closure.add(deposit.fraction)?)?;
            }
            error.charge(closure.add(row.escape_fraction)?)?;
            error.charge(closure.add(-1.0)?)?;
            let accumulation_error = error.used;
            error.charge(closure.hi.abs())?;
            error.charge(closure.lo.abs())?;
            closure_defects.push((closure, accumulation_error));
        }
        Ok(Self(Arc::new(ResponseData {
            data_id,
            rows,
            receiver_ids: receivers.into_iter().collect(),
            closure_defects,
        })))
    }

    /// Caller-supplied dataset and applicability identity.
    #[must_use]
    pub fn data_id(&self) -> &str {
        &self.0.data_id
    }

    /// Canonical source rows; supplied fractions remain unchanged.
    #[must_use]
    pub fn rows(&self) -> &[RadiationSourceResponse] {
        &self.0.rows
    }

    // Generated transport responses retain their measured preparation defects
    // in the same per-emitted-joule audit as supplied row closure. They cannot
    // silently discard a solve residual or endpoint-rounding adjustment.
    pub(super) fn add_preparation_roundoff(&mut self, fractions: &[f64]) -> Result<(), RadiationError> {
        if fractions.len() != self.0.rows.len() || !fractions.iter().copied().all(nonnegative) {
            return Err(RadiationError::InvalidInput("exact nonnegative preparation defect coverage required"));
        }
        let data = Arc::make_mut(&mut self.0);
        for ((_, uncertainty), &fraction) in data.closure_defects.iter_mut().zip(fractions) {
            let mut error = ErrorBudget::new(f64::MAX);
            error.charge(*uncertainty)?;
            error.charge(fraction)?;
            *uncertainty = error.used;
        }
        Ok(())
    }

    // Preserve the complete original coefficient audit when generated rows are
    // conservatively combined, including defects that could cancel in a sum.
    pub(super) fn row_roundoff(&self, index: usize) -> Option<f64> {
        let (closure, uncertainty) = self.0.closure_defects.get(index)?;
        let mut error = ErrorBudget::new(f64::MAX);
        for defect in [closure.hi.abs(), closure.lo.abs(), *uncertainty] {
            error.charge(defect).ok()?;
        }
        Some(error.used)
    }

    /// Unique physical emitting identities in ascending order.
    pub fn source_ids(&self) -> impl ExactSizeIterator<Item = u64> + '_ {
        self.0.rows.iter().map(|row| row.source_id)
    }

    /// Unique physical receiving identities in ascending order.
    pub fn receiver_ids(&self) -> impl ExactSizeIterator<Item = u64> + '_ {
        self.0.receiver_ids.iter().copied()
    }

    /// Distribute externally produced radiation energy using this exact response.
    ///
    /// Input may be unordered but must contain every source exactly once,
    /// including zero emitters. Sources, routes and receivers are reduced in
    /// ascending identity order independently of input order. Work/storage are
    /// sparse: O((sources + transfers + receivers) log(sources + receivers)).
    /// Every input, route and receiver polls the checkpoint. Callback errors are
    /// returned unchanged; refusal exposes no partial result and changes no state.
    pub fn apply(
        &self,
        emitted_by_source_j: &[(u64, f64)],
        energy_tolerance_j: f64,
        mut checkpoint: impl FnMut() -> Result<(), RadiationError>,
    ) -> Result<RadiationDeposition, RadiationError> {
        if !nonnegative(energy_tolerance_j) {
            return Err(RadiationError::InvalidInput(
                "energy tolerance must be finite and nonnegative",
            ));
        }
        if emitted_by_source_j.len() != self.0.rows.len() {
            return Err(RadiationError::InvalidInput(
                "source coverage must be exact",
            ));
        }
        checkpoint()?;
        let mut emissions = BTreeMap::new();
        for &(id, energy) in emitted_by_source_j {
            checkpoint()?;
            if !nonnegative(energy) {
                return Err(RadiationError::InvalidInput(
                    "emitted energy must be finite and nonnegative",
                ));
            }
            if emissions.insert(id, energy).is_some() {
                return Err(RadiationError::InvalidInput(
                    "duplicate emitted source identity",
                ));
            }
        }
        if !self.source_ids().eq(emissions.keys().copied()) {
            return Err(RadiationError::InvalidInput(
                "source coverage must be exact",
            ));
        }
        let mut receivers: BTreeMap<_, Sum> = BTreeMap::new();
        for id in self.receiver_ids() {
            checkpoint()?;
            receivers.insert(id, Sum::default());
        }
        let mut emitted = Sum::default();
        let mut escaped = Sum::default();
        let mut errors = ErrorBudget::new(energy_tolerance_j);
        for (row, (closure, uncertainty)) in self.0.rows.iter().zip(&self.0.closure_defects) {
            checkpoint()?;
            let energy = emissions[&row.source_id];
            errors.charge(emitted.add(energy)?)?;
            for defect in [closure.hi.abs(), closure.lo.abs(), *uncertainty] {
                errors.charge(scaled_bound(energy, defect)?)?;
            }
            for deposit in &row.deposits {
                checkpoint()?;
                let (heat, defect) = product(energy, deposit.fraction)?;
                errors.charge(defect)?;
                errors.charge(
                    receivers
                        .get_mut(&deposit.receiver_id)
                        .expect("admitted receiver")
                        .add(heat)?,
                )?;
            }
            let (export, defect) = product(energy, row.escape_fraction)?;
            errors.charge(defect)?;
            errors.charge(escaped.add(export)?)?;
        }
        let mut receiver_energy_j = BTreeMap::new();
        let mut deposited = Sum::default();
        for (id, energy) in receivers {
            checkpoint()?;
            let energy = energy.export(&mut errors)?;
            errors.charge(deposited.add(energy)?)?;
            receiver_energy_j.insert(id, energy);
        }
        let emitted_j = emitted.export(&mut errors)?;
        let deposited_j = deposited.export(&mut errors)?;
        let escaped_j = escaped.export(&mut errors)?;
        let mut residual = Sum {
            hi: -emitted_j,
            lo: 0.0,
        };
        errors.charge(residual.add(deposited_j)?)?;
        errors.charge(residual.add(escaped_j)?)?;
        let residual_j = residual.export(&mut errors)?;
        errors.charge(residual_j.abs())?;
        checkpoint()?;
        Ok(RadiationDeposition {
            receiver_energy_j,
            emitted_j,
            deposited_j,
            escaped_j,
            residual_j,
            roundoff_j: errors.used,
        })
    }
}

fn nonnegative(value: f64) -> bool {
    value.is_finite() && value >= 0.0
}
fn fraction(value: f64) -> bool {
    nonnegative(value) && value <= 1.0
}
fn finite(value: f64) -> Result<f64, RadiationError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(RadiationError::Arithmetic)
    }
}

// Error-free two-sum, followed by a two-component reduction. A discarded third
// component is reported instead of allowing cancellation to hide local loss.
fn two_sum(a: f64, b: f64) -> Result<(f64, f64), RadiationError> {
    let sum = finite(a + b)?;
    let v = sum - a;
    Ok((sum, finite((a - (sum - v)) + (b - v))?))
}

#[derive(Clone, Copy, Debug, Default)]
struct Sum {
    hi: f64,
    lo: f64,
}
impl Sum {
    fn add(&mut self, value: f64) -> Result<f64, RadiationError> {
        let (s, e) = two_sum(self.hi, value)?;
        let (t, u) = two_sum(self.lo, e)?;
        let (hi, lo) = two_sum(s, t)?;
        let (tail, lost) = two_sum(lo, u)?;
        let (hi, lo) = two_sum(hi, tail)?;
        *self = Self { hi, lo };
        Ok(lost.abs())
    }
    fn export(self, errors: &mut ErrorBudget) -> Result<f64, RadiationError> {
        let (value, error) = two_sum(self.hi, self.lo)?;
        errors.charge(error.abs())?;
        Ok(value)
    }
}

struct ErrorBudget {
    tolerance: f64,
    used: f64,
}
impl ErrorBudget {
    const fn new(tolerance: f64) -> Self {
        Self {
            tolerance,
            used: 0.0,
        }
    }
    fn charge(&mut self, error: f64) -> Result<(), RadiationError> {
        let (sum, remainder) = two_sum(self.used, error)?;
        // Round this nonnegative bound upwards, including tiny lost charges.
        self.used = if remainder > 0.0 {
            finite(sum.next_up())?
        } else {
            sum
        };
        if self.used > self.tolerance {
            return Err(RadiationError::Budget("energy-accounting defect"));
        }
        Ok(())
    }
}

// Product residuals are exact with FMA except at the underflow boundary. There
// integer significands detect exact products; otherwise one extra subnormal
// quantum conservatively bounds the unrepresentable part of the residual.
fn product(a: f64, b: f64) -> Result<(f64, f64), RadiationError> {
    let value = finite(a * b)?;
    if a > 0.0 && b > 0.0 && value == 0.0 {
        return Err(RadiationError::Arithmetic);
    }
    let mut error = a.mul_add(b, -value).abs();
    if value > 0.0
        && value < 2.0 * f64::MIN_POSITIVE / f64::EPSILON
        && error < f64::MIN_POSITIVE
        && !exact_product(a, b, value)
    {
        error = error.next_up();
    }
    Ok((value, error))
}

fn exact_product(a: f64, b: f64, product: f64) -> bool {
    fn parts(value: f64) -> (u128, i32) {
        let bits = value.to_bits();
        let field = ((bits >> 52) & 0x7ff) as i32;
        let mantissa = (bits & ((1_u64 << 52) - 1)) as u128;
        if field == 0 {
            (mantissa, -1074)
        } else {
            (mantissa | (1_u128 << 52), field - 1075)
        }
    }
    let (a, ae) = parts(a);
    let (b, be) = parts(b);
    let (p, pe) = parts(product);
    let ab = a * b;
    let az = ab.trailing_zeros();
    let pz = p.trailing_zeros();
    (ab >> az) == (p >> pz) && ae + be + az as i32 == pe + pz as i32
}

// A bound that underflows is rounded up instead of being mistaken for zero.
fn scaled_bound(a: f64, b: f64) -> Result<f64, RadiationError> {
    if a == 0.0 || b == 0.0 {
        return Ok(0.0);
    }
    if a * b == 0.0 {
        return Ok(f64::from_bits(1));
    }
    let (value, residual) = product(a, b)?;
    let (total, remainder) = two_sum(value, residual)?;
    if remainder > 0.0 {
        finite(total.next_up())
    } else {
        Ok(total)
    }
}

#[cfg(test)]
#[path = "radiation/tests.rs"]
mod tests;
