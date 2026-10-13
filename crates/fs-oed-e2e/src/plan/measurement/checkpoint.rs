//! Storage-neutral recovery of actual measurement studies (plan P7/P9).
//!
//! Store the returned bytes in the caller's durable store. Recovery replays the
//! original declarations and acquired readings through the SAME recommendation
//! and assimilation code. Cached posteriors are never decoded as authority.
//! One work/byte/poll accountant covers transport, initialization and all replay.
//! Same-build numerical replay is required; a checksum is not authentication.

use fs_blake3::DomainHasher;
use fs_exec::Cx;
use fs_qty::{QUANTITY_SPEC_ENCODED_LEN, QtyAny, QuantitySpec, semantic::SemanticQty};

use super::{MeasurementAlternative, MeasurementDecision, MeasurementError, MeasurementStudy, reserve};
use crate::{
    AcquisitionCost, AcquisitionCostSpec, ByteShape, CampaignBytePlan, CampaignProgress,
    CampaignWorkPlan, MAX_CAMPAIGN_CANDIDATES, MAX_CAMPAIGN_SENSORS, MAX_CANDIDATE_NAME_BYTES,
    ObjectiveValue, OedError, canonicalize_zero,
};

/// Exact transport version. Unknown versions refuse rather than guessing.
pub const MEASUREMENT_CHECKPOINT_VERSION: u32 = 1;
/// Bound checked before hashing, allocation or numerical reconstruction.
pub const MAX_MEASUREMENT_CHECKPOINT_BYTES: usize = 2 * 1024 * 1024;
const MAGIC: &[u8; 8] = b"FSOEDST\0";
const DOMAIN: &str = "org.frankensim.sensorforge.measurement-checkpoint.v1";
const HEADER: usize = 24;
const FIXED: usize = HEADER + 2 * QUANTITY_SPEC_ENCODED_LEN + 1 + 8 + 64;

impl MeasurementStudy {
    /// Encode initial typed declarations and every accepted actual reading.
    /// No I/O occurs. Cancellation/allocation failure leaves the study unchanged.
    /// The caller must durably store the COMPLETE frame before relying on it.
    pub fn to_checkpoint(&self, cx: &Cx<'_>) -> Result<Vec<u8>, MeasurementError> {
        let n = self.alternatives.len();
        let count = self.updates.len();
        // Bounds follow the sealed menu/reading caps, including 256-byte IDs.
        let bound = FIXED + n * (2 + MAX_CANDIDATE_NAME_BYTES + 32) + count * (32 + 8 + 2 + 256);
        let mut scope = Scope::new(n, 0, bound, cx)?;
        let mut writer = Writer::new(bound)?;
        writer.put(MAGIC)?;
        writer.u32(MEASUREMENT_CHECKPOINT_VERSION)?;
        writer.u32(n as u32)?;
        writer.u32(self.max_readings as u32)?;
        writer.u32(count as u32)?;
        writer.put(&self.alternatives[0].prior.spec().canonical_bytes())?;
        writer.put(&self.alternatives[0].cost.spec().canonical_bytes())?;
        writer.number(self.threshold.value())?;
        for alternative in &self.alternatives {
            scope.poll(cx)?;
            writer.text(&alternative.name)?;
            for value in [alternative.prior.value(), alternative.variance, alternative.noise, alternative.cost.value()] {
                writer.number(value)?;
            }
        }
        for update in &self.updates {
            scope.poll(cx)?;
            writer.put(update.request.identity.as_bytes())?;
            writer.number(update.outcome.value.value())?;
            writer.text(&update.outcome.instrument)?;
        }
        writer.put(self.history.as_bytes())?;
        scope.transport(writer.bytes.len())?;
        let digest = checksum(&writer.bytes, DOMAIN, &mut scope, cx)?;
        writer.put(&digest)?;
        scope.poll(cx)?;
        Ok(writer.bytes)
    }

    /// Restore the full study, rechecking every selected request and the final
    /// input/output history. This never calls a MeasurementProvider. Only a
    /// complete, numerically replayed study is returned; bytes survive refusal.
    ///
    /// Cost grows with accepted history and menu size. One cumulative accountant
    /// and the existing action-evaluation cap apply to the ENTIRE recovery, not
    /// a fresh budget for each row. A valid large checkpoint may exceed that cap.
    /// This is same-build/compatible-profile replay, not automatic migration.
    pub fn from_checkpoint(bytes: &[u8], cx: &Cx<'_>) -> Result<Self, MeasurementError> {
        let (n, max_readings, count) = header(bytes)?;
        let rounds = 2 * count + 1; // Initialization plus recommendation and update per row.
        let mut scope = Scope::new(n, rounds, bytes.len(), cx)?;
        let payload = verified_payload(bytes, DOMAIN, &mut scope, cx)?;
        let mut reader = Reader { bytes: payload, offset: HEADER };
        let objective_spec = reader.spec()?;
        let cost_role = reader.take(1)?[0];
        let cost_token = reader.take(QUANTITY_SPEC_ENCODED_LEN)?;
        let cost_spec = match cost_role {
            0 if cost_token.iter().all(|byte| *byte == 0) => AcquisitionCostSpec::RelativeWeight,
            1 => AcquisitionCostSpec::Quantity(QuantitySpec::from_canonical_bytes(cost_token)
                .map_err(|_| invalid("invalid checkpoint acquisition-cost schema"))?),
            _ => return Err(invalid("invalid checkpoint acquisition-cost role")),
        };
        let threshold = reader.number()?;
        let mut alternatives = reserve(n)?;
        for _ in 0..n {
            scope.poll(cx)?;
            let name = reader.text(MAX_CANDIDATE_NAME_BYTES)?;
            if alternatives.last().is_some_and(|previous: &MeasurementAlternative| previous.name >= name) {
                return Err(invalid("checkpoint alternatives are not in canonical unique order"));
            }
            let prior = objective(objective_spec, reader.number()?)?;
            let variance = reader.number()?;
            let noise = reader.number()?;
            let cost = cost(cost_spec, reader.number()?)?;
            let variance_dims = prior.spec().variance_dims()?;
            alternatives.push(MeasurementAlternative::new(name, prior,
                QtyAny::new(variance, variance_dims), QtyAny::new(noise, variance_dims), cost)?);
        }
        // No serialized mean, covariance, EVPI, score or Color is installed.
        let threshold = alternatives[0].prior.spec().decision_value(threshold)?;
        let mut study = Self::new_planned(alternatives, threshold, max_readings,
            scope.plan, scope.shape, &mut scope.progress, cx)?;
        for _ in 0..count {
            scope.poll(cx)?;
            let expected_request = reader.array::<32>()?;
            let value = objective(objective_spec, reader.number()?)?;
            let instrument = reader.text(256)?;
            let MeasurementDecision::Acquire(request) = study.recommend_planned(
                scope.plan, scope.shape, &mut scope.progress, cx)? else {
                return Err(invalid("checkpoint contains a reading after the replayed study stopped"));
            };
            if request.identity.as_bytes() != &expected_request {
                return Err(invalid("checkpoint request differs from numerical replay"));
            }
            let outcome = request.reading(value, instrument)?;
            study.assimilate_planned(&request, &outcome, scope.plan, scope.shape, &mut scope.progress, cx)?;
        }
        let expected_history = reader.array::<32>()?;
        if reader.offset != payload.len() {
            return Err(invalid("trailing checkpoint data"));
        }
        if study.history.as_bytes() != &expected_history {
            return Err(invalid("checkpoint history differs from numerical replay"));
        }
        scope.transport(payload.len())?;
        scope.poll(cx)?;
        Ok(study)
    }
}

fn header(bytes: &[u8]) -> Result<(usize, usize, usize), MeasurementError> {
    if bytes.len() < FIXED || bytes.len() > MAX_MEASUREMENT_CHECKPOINT_BYTES {
        return Err(invalid("checkpoint length is outside the bounded transport"));
    }
    let mut reader = Reader { bytes, offset: 0 };
    if reader.take(8)? != MAGIC || reader.u32()? != MEASUREMENT_CHECKPOINT_VERSION {
        return Err(invalid("unknown measurement checkpoint type or version"));
    }
    let n = reader.u32()? as usize;
    let max_readings = reader.u32()? as usize;
    let count = reader.u32()? as usize;
    if n == 0 || n > MAX_CAMPAIGN_CANDIDATES || max_readings > MAX_CAMPAIGN_SENSORS || count > max_readings {
        return Err(invalid("checkpoint menu or history count exceeds its bounds"));
    }
    // Require enough bytes for even the shortest legal rows before allocation.
    if bytes.len() < FIXED + n * (2 + 1 + 32) + count * (32 + 8 + 2 + 1) {
        return Err(invalid("truncated measurement checkpoint"));
    }
    Ok((n, max_readings, count))
}

struct Scope<'s> {
    plan: CampaignWorkPlan,
    shape: ByteShape,
    progress: CampaignProgress<'s>,
}
impl<'s> Scope<'s> {
    fn new(n: usize, rounds: usize, wire_bytes: usize, cx: &Cx<'s>) -> Result<Self, MeasurementError> {
        if wire_bytes > MAX_MEASUREMENT_CHECKPOINT_BYTES {
            return Err(invalid("checkpoint exceeds the transport byte cap"));
        }
        let mut plan = CampaignWorkPlan::checked(n, rounds)?;
        let shape = ByteShape { n: n as u128, sum_name: (n * MAX_CANDIDATE_NAME_BYTES) as u128,
            max_name: MAX_CANDIDATE_NAME_BYTES as u128 };
        let mut bytes = CampaignBytePlan::checked(n, rounds, n * MAX_CANDIDATE_NAME_BYTES, MAX_CANDIDATE_NAME_BYTES)?;
        // Covers encode/decode, bounded copying, frame hash, and fixed fields.
        // Dimensions/counts are already capped; still use checked arithmetic.
        let extra = 4 * wire_bytes as u128 + 4096;
        plan.admitted_work_units = plan.admitted_work_units.checked_add(extra)
            .ok_or_else(|| invalid("checkpoint work size overflow"))?;
        bytes.admitted_byte_units = bytes.admitted_byte_units.checked_add(extra)
            .ok_or_else(|| invalid("checkpoint byte size overflow"))?;
        let cost = u64::try_from(plan.admitted_work_units).map_err(|_| invalid("checkpoint work size overflow"))?;
        let progress = CampaignProgress::admit(cx, 0, cost, bytes)?;
        let mut scope = Self { plan, shape, progress };
        scope.poll(cx)?;
        Ok(scope)
    }
    fn poll(&mut self, cx: &Cx<'_>) -> Result<(), MeasurementError> {
        self.progress.checkpoint(cx, self.plan, "measurement checkpoint")?;
        Ok(())
    }
    fn transport(&mut self, bytes: usize) -> Result<(), MeasurementError> {
        self.progress.advance(bytes as u128);
        self.progress.charge_bytes("measurement checkpoint transport", bytes as u128)?;
        Ok(())
    }
}

fn checksum(bytes: &[u8], domain: &str, scope: &mut Scope<'_>, cx: &Cx<'_>) -> Result<[u8; 32], MeasurementError> {
    let mut hash = DomainHasher::new(domain);
    for chunk in bytes.chunks(1024) {
        scope.poll(cx)?;
        hash.update(chunk);
        scope.transport(chunk.len())?;
    }
    Ok(*hash.finalize().as_bytes())
}

fn verified_payload<'a>(bytes: &'a [u8], domain: &str, scope: &mut Scope<'_>, cx: &Cx<'_>) -> Result<&'a [u8], MeasurementError> {
    let split = bytes.len().checked_sub(32).ok_or_else(|| invalid("truncated checkpoint checksum"))?;
    let (payload, digest) = bytes.split_at(split);
    if checksum(payload, domain, scope, cx)?.as_slice() != digest {
        return Err(invalid("checkpoint checksum mismatch"));
    }
    Ok(payload)
}

fn objective(spec: QuantitySpec, value: f64) -> Result<ObjectiveValue, MeasurementError> {
    let quantity = QtyAny::new(value, spec.dims());
    match spec.semantic_type() {
        Some(kind) => Ok(ObjectiveValue::semantic(SemanticQty::new(quantity, kind)
            .map_err(|_| invalid("checkpoint objective violates its semantic schema"))?)),
        None => Ok(ObjectiveValue::dimensional(quantity)?),
    }
}

fn cost(spec: AcquisitionCostSpec, value: f64) -> Result<AcquisitionCost, MeasurementError> {
    match spec {
        AcquisitionCostSpec::RelativeWeight => Ok(AcquisitionCost::relative(value)?),
        AcquisitionCostSpec::Quantity(spec) => {
            let quantity = QtyAny::new(value, spec.dims());
            match spec.semantic_type() {
                Some(kind) => Ok(AcquisitionCost::semantic(SemanticQty::new(quantity, kind)
                    .map_err(|_| invalid("checkpoint cost violates its semantic schema"))?)?),
                None => Ok(AcquisitionCost::dimensional(quantity)?),
            }
        }
    }
}

struct Writer { bytes: Vec<u8>, bound: usize }
impl Writer {
    fn new(bound: usize) -> Result<Self, MeasurementError> {
        Ok(Self { bytes: reserve(bound)?, bound })
    }
    fn put(&mut self, bytes: &[u8]) -> Result<(), MeasurementError> {
        if bytes.len() > self.bound.saturating_sub(self.bytes.len()) {
            return Err(invalid("checkpoint encoding exceeds admitted storage"));
        }
        self.bytes.extend_from_slice(bytes);
        Ok(())
    }
    fn u32(&mut self, value: u32) -> Result<(), MeasurementError> { self.put(&value.to_le_bytes()) }
    fn number(&mut self, value: f64) -> Result<(), MeasurementError> {
        self.put(&canonicalize_zero(value).to_le_bytes())
    }
    fn text(&mut self, value: &str) -> Result<(), MeasurementError> {
        let len = u16::try_from(value.len()).map_err(|_| invalid("checkpoint identity too long"))?;
        self.put(&len.to_le_bytes())?;
        self.put(value.as_bytes())
    }
}

struct Reader<'a> { bytes: &'a [u8], offset: usize }
impl<'a> Reader<'a> {
    fn take(&mut self, count: usize) -> Result<&'a [u8], MeasurementError> {
        let end = self.offset.checked_add(count).ok_or_else(|| invalid("checkpoint offset overflow"))?;
        let value = self.bytes.get(self.offset..end).ok_or_else(|| invalid("truncated checkpoint field"))?;
        self.offset = end;
        Ok(value)
    }
    fn array<const N: usize>(&mut self) -> Result<[u8; N], MeasurementError> {
        self.take(N)?.try_into().map_err(|_| invalid("checkpoint field length mismatch"))
    }
    fn u32(&mut self) -> Result<u32, MeasurementError> { Ok(u32::from_le_bytes(self.array()?)) }
    fn number(&mut self) -> Result<f64, MeasurementError> {
        let value = f64::from_le_bytes(self.array()?);
        if !value.is_finite() || value.to_bits() == (-0.0_f64).to_bits() {
            return Err(invalid("checkpoint number is not finite canonical data"));
        }
        Ok(value)
    }
    fn spec(&mut self) -> Result<QuantitySpec, MeasurementError> {
        QuantitySpec::from_canonical_bytes(self.take(QUANTITY_SPEC_ENCODED_LEN)?)
            .map_err(|_| invalid("invalid checkpoint objective schema"))
    }
    fn text(&mut self, max: usize) -> Result<String, MeasurementError> {
        let len = u16::from_le_bytes(self.array()?) as usize;
        if len == 0 || len > max { return Err(invalid("checkpoint identity length exceeds its bounds")); }
        let text = core::str::from_utf8(self.take(len)?).map_err(|_| invalid("checkpoint identity is not UTF-8"))?;
        let mut value = String::new();
        value.try_reserve_exact(len).map_err(|_| OedError::OutputAllocationRefused { what: "checkpoint identity" })?;
        value.push_str(text);
        Ok(value)
    }
}

fn invalid(reason: &'static str) -> MeasurementError { MeasurementError::Invalid(reason) }

const ACQUIRED_MAGIC: &[u8; 8] = b"FSOEDAC\0";
const ACQUIRED_DOMAIN: &str = "org.frankensim.sensorforge.acquired-measurement.v1";
const ACQUIRED_FIXED: usize = 8 + 4 + 64 + QUANTITY_SPEC_ENCODED_LEN + 8 + 2 + 32;

pub(in crate::plan) fn acquisition_to_checkpoint(
    acquired: &crate::plan::acquisition::AcquiredMeasurement, cx: &Cx<'_>,
) -> Result<Vec<u8>, MeasurementError> {
    let request = acquired.request();
    let outcome = acquired.outcome();
    if outcome.request_identity() != request.identity() {
        return Err(MeasurementError::RequestMismatch);
    }
    let bound = ACQUIRED_FIXED + outcome.instrument().len();
    let mut scope = Scope::new(1, 0, bound, cx)?;
    let mut writer = Writer::new(bound)?;
    writer.put(ACQUIRED_MAGIC)?;
    writer.u32(MEASUREMENT_CHECKPOINT_VERSION)?;
    writer.put(request.history.as_bytes())?;
    writer.put(request.identity().as_bytes())?;
    writer.put(&outcome.value().spec().canonical_bytes())?;
    writer.number(outcome.value().value())?;
    writer.text(outcome.instrument())?;
    scope.transport(writer.bytes.len())?;
    let digest = checksum(&writer.bytes, ACQUIRED_DOMAIN, &mut scope, cx)?;
    writer.put(&digest)?;
    scope.poll(cx)?;
    Ok(writer.bytes)
}

pub(in crate::plan) fn acquisition_from_checkpoint(
    study: &MeasurementStudy, bytes: &[u8], cx: &Cx<'_>,
) -> Result<(super::MeasurementRequest, super::MeasurementOutcome), MeasurementError> {
    if !(ACQUIRED_FIXED + 1..=ACQUIRED_FIXED + 256).contains(&bytes.len()) {
        return Err(invalid("acquisition checkpoint length is outside its bounds"));
    }
    let mut reader = Reader { bytes, offset: 0 };
    if reader.take(8)? != ACQUIRED_MAGIC || reader.u32()? != MEASUREMENT_CHECKPOINT_VERSION {
        return Err(invalid("unknown acquired measurement type or version"));
    }
    let mut scope = Scope::new(study.alternatives.len(), 1, bytes.len(), cx)?;
    let payload = verified_payload(bytes, ACQUIRED_DOMAIN, &mut scope, cx)?;
    let mut reader = Reader { bytes: payload, offset: 12 };
    let source_history = reader.array::<32>()?;
    let expected_request = reader.array::<32>()?;
    if study.history.as_bytes() != &source_history { return Err(MeasurementError::StaleRequest); }
    let spec = reader.spec()?;
    let value = objective(spec, reader.number()?)?;
    let instrument = reader.text(256)?;
    if reader.offset != payload.len() { return Err(invalid("trailing acquisition checkpoint data")); }
    scope.transport(payload.len())?;
    let MeasurementDecision::Acquire(request) = study.recommend_planned(
        scope.plan, scope.shape, &mut scope.progress, cx)? else {
        return Err(MeasurementError::StaleRequest);
    };
    if request.identity().as_bytes() != &expected_request { return Err(MeasurementError::RequestMismatch); }
    let outcome = request.reading(value, instrument)?;
    scope.poll(cx)?;
    Ok((request, outcome))
}
