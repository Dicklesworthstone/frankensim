//! Cumulative acquisition allowance, separate from numerical Cx budgets.
//!
//! Reserve a named attempt, persist the account, THEN consume its single-use
//! permit to do I/O. Failed, cancelled or abandoned attempts are never refunded.
//! Numerical application of an acquired reading does not reserve another attempt.

use fs_blake3::{ContentHash, DomainHasher};
use fs_evidence::color_leaf_identity_reason;
use fs_exec::{AdmittedBudget, Cx};
use crate::{AcquisitionCost, AcquisitionCostSpec, MAX_CAMPAIGN_SENSORS, OedError};
use super::{AcquiredMeasurement, AcquisitionError, MeasurementProvider};
use crate::plan::measurement::{MeasurementError, MeasurementRequest, MeasurementStudy};

const MAGIC: &[u8; 8] = b"FSOEDBD\0";
const VERSION: u32 = 1;
const DOMAIN: &str = "org.frankensim.sensorforge.acquisition-account.v1";
const HEADER: usize = 12 + crate::ACQUISITION_COST_SPEC_ENCODED_LEN + 8 + 4 + 4;
const MAX_FRAME: usize = HEADER + MAX_CAMPAIGN_SENSORS * (2 + 256 + 32 + 8) + 32;

/// No reservation or provider call occurs when admission fails.
#[derive(Debug, Clone, PartialEq)]
pub enum AcquisitionLimitError {
    /// Invalid identity, cap, allocation or saved transport.
    Invalid(&'static str),
    /// Scientific request admission or numerical execution refused.
    Measurement(MeasurementError),
    /// An allowance cannot silently reinterpret another resource or cost role.
    Schema { expected: AcquisitionCostSpec, actual: AcquisitionCostSpec },
    /// An earlier attempt already consumes this identity and its charge.
    DuplicateAttempt { id: String },
    /// No more named attempts can be reserved.
    AttemptLimit { maximum: usize },
    /// The conservative next cost would exceed the allowance.
    CostLimit { maximum: f64, committed_upper: f64, requested: f64 },
}
impl From<MeasurementError> for AcquisitionLimitError {
    fn from(error: MeasurementError) -> Self { Self::Measurement(error) }
}
impl core::fmt::Display for AcquisitionLimitError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Invalid(reason) => f.write_str(reason),
            Self::Measurement(error) => write!(f, "{error}"),
            Self::Schema { expected, actual } => write!(f, "acquisition cost schema {actual:?} differs from allowance {expected:?}"),
            Self::DuplicateAttempt { id } => write!(f, "acquisition attempt {id} was already reserved"),
            Self::AttemptLimit { maximum } => write!(f, "acquisition attempt limit {maximum} reached"),
            Self::CostLimit { maximum, committed_upper, requested } => write!(f, "acquisition cost {requested} exceeds allowance {maximum} after committed upper bound {committed_upper}"),
        }
    }
}
impl std::error::Error for AcquisitionLimitError {}

/// One budget commitment. Its presence does NOT assert successful physical I/O.
#[derive(Debug, Clone, PartialEq)]
pub struct AcquisitionCharge { attempt: String, request: ContentHash, cost: f64 }
impl AcquisitionCharge {
    #[must_use]
    pub fn attempt_id(&self) -> &str { &self.attempt }
    #[must_use]
    pub const fn request_identity(&self) -> ContentHash { self.request }
    /// Declared cost under the account's exact cost role/schema.
    #[must_use]
    pub const fn cost(&self) -> f64 { self.cost }
}

/// A finite allowance for one explicit additive acquisition resource.
/// Intentionally not Clone: copying the owner must not casually renew its budget.
/// This is application accounting, not authentication or a distributed lock.
#[derive(Debug, PartialEq)]
pub struct AcquisitionAccount {
    limit: AcquisitionCost,
    max_attempts: usize,
    charges: Vec<AcquisitionCharge>,
    committed_upper: f64,
}

/// A reserved attempt usable for at most one provider call. Intentionally not
/// Clone and not serializable: after a crash, an attempt of unknown I/O status
/// cannot silently be reissued. Recover its saved reading or explicitly reserve
/// a differently named attempt, consuming another allowance.
#[derive(Debug)]
#[must_use = "a dropped reservation stays charged; persist its account before I/O"]
pub struct ReservedAcquisition { attempt: String, request: MeasurementRequest }
impl ReservedAcquisition {
    #[must_use]
    pub fn attempt_id(&self) -> &str { &self.attempt }
    #[must_use]
    pub const fn request(&self) -> &MeasurementRequest { &self.request }
    /// Consume the permit and use the ordinary acquisition path. Staleness and
    /// cancellation are rechecked BEFORE the provider; a returned reading is
    /// retained in AcquisitionError::After on post-I/O refusal. No refund occurs.
    pub fn acquire<P: MeasurementProvider>(self, study: &MeasurementStudy,
        provider: &mut P, cx: &Cx<'_>) -> Result<AcquiredMeasurement, AcquisitionError<P::Error>> {
        study.acquire(&self.request, provider, cx)
    }
}

impl AcquisitionAccount {
    /// A positive finite allowance in the SAME exact role/schema as sensor costs.
    /// A zero attempt cap disables acquisition. At most 4,096 named attempts.
    pub fn new(limit: AcquisitionCost, max_attempts: usize) -> Result<Self, AcquisitionLimitError> {
        if max_attempts > MAX_CAMPAIGN_SENSORS { return Err(invalid("acquisition attempt cap exceeds 4096")); }
        Ok(Self { limit, max_attempts, charges: Vec::new(), committed_upper: 0.0 })
    }
    #[must_use]
    pub const fn limit(&self) -> AcquisitionCost { self.limit }
    #[must_use]
    pub const fn max_attempts(&self) -> usize { self.max_attempts }
    #[must_use]
    pub fn charges(&self) -> &[AcquisitionCharge] { &self.charges }
    /// Conservative upper sum, in the allowance's declared coherent units.
    /// Not an assertion of invoiced or physically consumed resource.
    #[must_use]
    pub const fn committed_cost_upper(&self) -> f64 { self.committed_upper }

    /// Reserve BEFORE external I/O. Duplicate IDs always refuse, even if their
    /// earlier permit was abandoned or I/O failed. An explicit new ID is required
    /// for a deliberate retry. No scientific state is changed or outcome invented.
    /// Save to_checkpoint() durably before consuming the returned permit when
    /// restart-safe spending matters; an unsaved debit cannot survive a crash.
    pub fn reserve(&mut self, study: &MeasurementStudy, request: &MeasurementRequest,
        attempt_id: impl Into<String>, cx: &Cx<'_>) -> Result<ReservedAcquisition, AcquisitionLimitError> {
        study.check_current(request)?;
        let attempt = attempt_id.into();
        identity(&attempt)?;
        let cost = request.acquisition_cost();
        if cost.spec() != self.limit.spec() {
            return Err(AcquisitionLimitError::Schema { expected: self.limit.spec(), actual: cost.spec() });
        }
        let mut scope = admit(cx, self.charges.len() as u64 + 8)?;
        let next = self.prepare(&attempt, cost.value(), &mut scope, cx)?;
        let retained_id = copy(&attempt)?;
        let permit = ReservedAcquisition { attempt, request: request.clone() };
        self.charges.try_reserve(1).map_err(|_| invalid("acquisition charge allocation refused"))?;
        poll(&mut scope, cx)?;
        // All fallible work precedes this one-way commitment. Never refund it
        // based on numerical/provider failure or whether a permit is dropped.
        self.charges.push(AcquisitionCharge { attempt: retained_id, request: request.identity(), cost: cost.value() });
        self.committed_upper = next;
        Ok(permit)
    }

    fn prepare(&self, attempt: &str, cost: f64, scope: &mut AdmittedBudget<'_>, cx: &Cx<'_>) -> Result<f64, AcquisitionLimitError> {
        for (index, charge) in self.charges.iter().enumerate() {
            if index % 32 == 0 { poll(scope, cx)?; }
            charge_work(scope, 1)?;
            if charge.attempt == attempt {
                return Err(AcquisitionLimitError::DuplicateAttempt { id: copy(attempt)? });
            }
        }
        if self.charges.len() >= self.max_attempts {
            return Err(AcquisitionLimitError::AttemptLimit { maximum: self.max_attempts });
        }
        let next = add_upper(self.committed_upper, cost);
        if !next.is_finite() || next > self.limit.value() {
            return Err(AcquisitionLimitError::CostLimit { maximum: self.limit.value(), committed_upper: self.committed_upper, requested: cost });
        }
        Ok(next)
    }

    /// Encode allowance and all commitments, including attempts with NO reading.
    /// The caller owns durable storage and must not roll back to an older account.
    pub fn to_checkpoint(&self, cx: &Cx<'_>) -> Result<Vec<u8>, AcquisitionLimitError> {
        let bound = HEADER + self.charges.iter().map(|c| 42 + c.attempt.len()).sum::<usize>() + 32;
        let mut scope = admit(cx, 4 * bound as u64 + 8)?;
        let mut bytes = Vec::new();
        bytes.try_reserve_exact(bound).map_err(|_| invalid("acquisition checkpoint allocation refused"))?;
        bytes.extend_from_slice(MAGIC); bytes.extend_from_slice(&VERSION.to_le_bytes());
        bytes.extend_from_slice(&self.limit.spec().canonical_bytes());
        bytes.extend_from_slice(&self.limit.value().to_le_bytes());
        bytes.extend_from_slice(&(self.max_attempts as u32).to_le_bytes());
        bytes.extend_from_slice(&(self.charges.len() as u32).to_le_bytes());
        for charge in &self.charges {
            poll(&mut scope, cx)?;
            bytes.extend_from_slice(&(charge.attempt.len() as u16).to_le_bytes());
            bytes.extend_from_slice(charge.attempt.as_bytes());
            bytes.extend_from_slice(charge.request.as_bytes());
            bytes.extend_from_slice(&charge.cost.to_le_bytes());
        }
        let digest = checksum(&bytes, &mut scope, cx)?;
        bytes.extend_from_slice(digest.as_bytes());
        charge_work(&mut scope, bound as u64)?;
        poll(&mut scope, cx)?;
        Ok(bytes)
    }

    /// Recover commitments WITHOUT calling a provider or trusting a stored sum.
    /// The caller must supply the original allowance/cap; mismatches refuse.
    /// Duplicate attempts, nonfinite costs and over-budget histories refuse.
    /// A checksum is not authentication or protection against storage rollback.
    pub fn from_checkpoint(limit: AcquisitionCost, max_attempts: usize, bytes: &[u8], cx: &Cx<'_>) -> Result<Self, AcquisitionLimitError> {
        let mut account = Self::new(limit, max_attempts)?;
        if !(HEADER + 32..=MAX_FRAME).contains(&bytes.len()) { return Err(invalid("acquisition account frame length invalid")); }
        let payload_len = bytes.len() - 32;
        let mut reader = Reader { bytes: &bytes[..payload_len], position: 0 };
        if reader.take(8)? != MAGIC || reader.u32()? != VERSION { return Err(invalid("unknown acquisition account type or version")); }
        if reader.take(crate::ACQUISITION_COST_SPEC_ENCODED_LEN)? != limit.spec().canonical_bytes().as_slice()
            || reader.take(8)? != limit.value().to_le_bytes().as_slice() || reader.u32()? as usize != max_attempts {
            return Err(invalid("saved acquisition allowance differs from supplied allowance"));
        }
        let n = reader.u32()? as usize;
        if n > max_attempts || payload_len < HEADER + 43 * n { return Err(invalid("invalid acquisition account count")); }
        let work = 4 * bytes.len() as u64 + (n as u64) * (n as u64 + 1) / 2 + 16;
        let mut scope = admit(cx, work)?;
        if checksum(reader.bytes, &mut scope, cx)?.as_bytes().as_slice() != &bytes[payload_len..] {
            return Err(invalid("acquisition account checksum mismatch"));
        }
        account.charges.try_reserve_exact(n).map_err(|_| invalid("acquisition account allocation refused"))?;
        for _ in 0..n {
            poll(&mut scope, cx)?;
            let length = u16::from_le_bytes(reader.array()?) as usize;
            if length == 0 || length > 256 { return Err(invalid("invalid acquisition attempt identity length")); }
            let attempt = copy(core::str::from_utf8(reader.take(length)?).map_err(|_| invalid("acquisition attempt is not UTF-8"))?)?;
            identity(&attempt)?;
            let request = ContentHash(reader.array()?);
            let cost = f64::from_le_bytes(reader.array()?);
            if !cost.is_finite() || cost <= 0.0 { return Err(invalid("invalid saved acquisition cost")); }
            let next = account.prepare(&attempt, cost, &mut scope, cx)?;
            account.charges.push(AcquisitionCharge { attempt, request, cost });
            account.committed_upper = next;
        }
        if reader.position != payload_len { return Err(invalid("trailing acquisition account data")); }
        charge_work(&mut scope, payload_len as u64)?;
        poll(&mut scope, cx)?;
        Ok(account)
    }
}

// Error-free TwoSum recovers the round-to-nearest error for finite positive
// operands. Advance one representable value only when the rounded sum is low.
// Repeating this gives an upper bound without ever losing a tiny positive cost.
fn add_upper(left: f64, right: f64) -> f64 {
    let sum = left + right;
    if !sum.is_finite() { return sum; }
    let part = sum - left;
    let error = (left - (sum - part)) + (right - part);
    if error > 0.0 { f64::from_bits(sum.to_bits() + 1) } else { sum }
}
fn identity(value: &str) -> Result<(), AcquisitionLimitError> {
    if value.len() > 256 || color_leaf_identity_reason(value).is_some() { return Err(invalid("invalid acquisition attempt identity")); }
    Ok(())
}
fn copy(value: &str) -> Result<String, AcquisitionLimitError> {
    let mut result = String::new();
    result.try_reserve_exact(value.len()).map_err(|_| invalid("acquisition identity allocation refused"))?;
    result.push_str(value); Ok(result)
}
fn admit<'s>(cx: &Cx<'s>, work: u64) -> Result<AdmittedBudget<'s>, AcquisitionLimitError> {
    let mut scope = AdmittedBudget::admit_ambient(cx, work).map_err(execution)?;
    poll(&mut scope, cx)?; Ok(scope)
}
fn poll(scope: &mut AdmittedBudget<'_>, cx: &Cx<'_>) -> Result<(), AcquisitionLimitError> {
    scope.checkpoint("acquisition allowance", cx).map_err(execution)
}
fn charge_work(scope: &mut AdmittedBudget<'_>, work: u64) -> Result<(), AcquisitionLimitError> {
    scope.charge_cost("acquisition allowance", work).map_err(execution)
}
fn execution(error: fs_exec::BudgetRefusal) -> AcquisitionLimitError {
    AcquisitionLimitError::Measurement(MeasurementError::Oed(OedError::BudgetRefused(error)))
}
fn invalid(reason: &'static str) -> AcquisitionLimitError { AcquisitionLimitError::Invalid(reason) }
fn checksum(bytes: &[u8], scope: &mut AdmittedBudget<'_>, cx: &Cx<'_>) -> Result<ContentHash, AcquisitionLimitError> {
    let mut hash = DomainHasher::new(DOMAIN);
    for chunk in bytes.chunks(1024) { poll(scope, cx)?; hash.update(chunk); charge_work(scope, chunk.len() as u64)?; }
    Ok(hash.finalize())
}
struct Reader<'a> { bytes: &'a [u8], position: usize }
impl<'a> Reader<'a> {
    fn take(&mut self, n: usize) -> Result<&'a [u8], AcquisitionLimitError> {
        let end = self.position.checked_add(n).ok_or_else(|| invalid("acquisition frame offset overflow"))?;
        let result = self.bytes.get(self.position..end).ok_or_else(|| invalid("truncated acquisition account"))?;
        self.position = end; Ok(result)
    }
    fn array<const N: usize>(&mut self) -> Result<[u8; N], AcquisitionLimitError> {
        self.take(N)?.try_into().map_err(|_| invalid("invalid acquisition account field"))
    }
    fn u32(&mut self) -> Result<u32, AcquisitionLimitError> { Ok(u32::from_le_bytes(self.array()?)) }
}
