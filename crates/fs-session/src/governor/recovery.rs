//! Strict binary codecs and restart-only recovery for durable session receipts.
//!
//! Event JSON remains human-readable audit evidence; this codec is the typed,
//! bounded reconstruction source. Every decoder consumes its complete input
//! and refuses unknown schema versions or trailing bytes.

use super::*;
use crate::token::{MAX_CAPABILITY_OP_BYTES, MAX_CAPABILITY_OPS, MAX_LEDGER_SCOPE_BYTES};
use fs_ledger::EventRow;

pub(super) const TERMINAL_SCHEMA_VERSION: u32 = 1;
pub(super) const KIND_OPEN: &str = "session-open";
pub(super) const KIND_METER: &str = "meter-report";
pub(super) const KIND_PRESSURE: &str = "pressure-action";
pub(super) const KIND_SUBMISSION: &str = "submission";
pub(super) const KIND_PAUSE_ACK: &str = "pause-acknowledgement";
pub(super) const KIND_RESUME_ACTIVATION: &str = "resume-activation";
pub(crate) const KIND_PROGRAM_RISK_REPORT: &str = "program-risk-report";

const MAX_CODEC_STRING_BYTES: usize = MAX_IDEMPOTENCY_INPUT_BYTES;

#[derive(Default)]
struct Encoder {
    bytes: Vec<u8>,
}

impl Encoder {
    fn schema(&mut self) {
        self.u32(TERMINAL_SCHEMA_VERSION);
    }

    fn u8(&mut self, value: u8) {
        self.bytes.push(value);
    }

    fn u32(&mut self, value: u32) {
        self.bytes.extend_from_slice(&value.to_le_bytes());
    }

    fn u64(&mut self, value: u64) {
        self.bytes.extend_from_slice(&value.to_le_bytes());
    }

    fn i64(&mut self, value: i64) {
        self.bytes.extend_from_slice(&value.to_le_bytes());
    }

    fn hash(&mut self, value: fs_blake3::ContentHash) {
        self.bytes.extend_from_slice(value.as_bytes());
    }

    fn string(&mut self, value: &str) {
        self.u64(u64::try_from(value.len()).expect("bounded codec string fits u64"));
        self.bytes.extend_from_slice(value.as_bytes());
    }

    fn optional_hash(&mut self, value: Option<fs_blake3::ContentHash>) {
        match value {
            Some(value) => {
                self.u8(1);
                self.hash(value);
            }
            None => self.u8(0),
        }
    }

    fn optional_u64(&mut self, value: Option<u64>) {
        match value {
            Some(value) => {
                self.u8(1);
                self.u64(value);
            }
            None => self.u8(0),
        }
    }

    fn optional_i64(&mut self, value: Option<i64>) {
        match value {
            Some(value) => {
                self.u8(1);
                self.i64(value);
            }
            None => self.u8(0),
        }
    }

    fn finish(self) -> Vec<u8> {
        self.bytes
    }
}

struct Decoder<'a> {
    bytes: &'a [u8],
    offset: usize,
    kind: &'static str,
    authority: fs_blake3::ContentHash,
}

impl<'a> Decoder<'a> {
    fn new(
        bytes: &'a [u8],
        kind: &'static str,
        authority: fs_blake3::ContentHash,
    ) -> Result<Self, SessionError> {
        let mut decoder = Self {
            bytes,
            offset: 0,
            kind,
            authority,
        };
        let found = decoder.u32()?;
        if found > TERMINAL_SCHEMA_VERSION {
            return Err(SessionError::UnsupportedTerminalSchema {
                found,
                supported: TERMINAL_SCHEMA_VERSION,
            });
        }
        if found != TERMINAL_SCHEMA_VERSION {
            return Err(decoder.corrupt(format!(
                "terminal codec schema v{found} is not the required v{TERMINAL_SCHEMA_VERSION}"
            )));
        }
        Ok(decoder)
    }

    fn corrupt(&self, detail: impl Into<String>) -> SessionError {
        SessionError::TerminalCorrupt {
            kind: self.kind,
            authority: self.authority,
            detail: detail.into(),
        }
    }

    fn take<const N: usize>(&mut self) -> Result<[u8; N], SessionError> {
        let end = self
            .offset
            .checked_add(N)
            .ok_or_else(|| self.corrupt("codec offset overflow"))?;
        let slice = self
            .bytes
            .get(self.offset..end)
            .ok_or_else(|| self.corrupt(format!("truncated fixed-width {N}-byte field")))?;
        self.offset = end;
        slice
            .try_into()
            .map_err(|_| self.corrupt(format!("invalid fixed-width {N}-byte field")))
    }

    fn u8(&mut self) -> Result<u8, SessionError> {
        Ok(self.take::<1>()?[0])
    }

    fn u32(&mut self) -> Result<u32, SessionError> {
        Ok(u32::from_le_bytes(self.take()?))
    }

    fn u64(&mut self) -> Result<u64, SessionError> {
        Ok(u64::from_le_bytes(self.take()?))
    }

    fn i64(&mut self) -> Result<i64, SessionError> {
        Ok(i64::from_le_bytes(self.take()?))
    }

    fn hash(&mut self) -> Result<fs_blake3::ContentHash, SessionError> {
        Ok(fs_blake3::ContentHash(self.take()?))
    }

    fn string(&mut self, max_bytes: usize) -> Result<String, SessionError> {
        let len = usize::try_from(self.u64()?)
            .map_err(|_| self.corrupt("string length does not fit usize"))?;
        if len > max_bytes {
            return Err(self.corrupt(format!(
                "string length {len} exceeds the {max_bytes}-byte codec bound"
            )));
        }
        let end = self
            .offset
            .checked_add(len)
            .ok_or_else(|| self.corrupt("string offset overflow"))?;
        let value = self
            .bytes
            .get(self.offset..end)
            .ok_or_else(|| self.corrupt("truncated string field"))?;
        self.offset = end;
        std::str::from_utf8(value)
            .map(str::to_owned)
            .map_err(|_| self.corrupt("string field is not UTF-8"))
    }

    fn optional_hash(&mut self) -> Result<Option<fs_blake3::ContentHash>, SessionError> {
        match self.u8()? {
            0 => Ok(None),
            1 => self.hash().map(Some),
            value => Err(self.corrupt(format!("invalid optional-hash tag {value}"))),
        }
    }

    fn optional_u64(&mut self) -> Result<Option<u64>, SessionError> {
        match self.u8()? {
            0 => Ok(None),
            1 => self.u64().map(Some),
            value => Err(self.corrupt(format!("invalid optional-u64 tag {value}"))),
        }
    }

    fn optional_i64(&mut self) -> Result<Option<i64>, SessionError> {
        match self.u8()? {
            0 => Ok(None),
            1 => self.i64().map(Some),
            value => Err(self.corrupt(format!("invalid optional-i64 tag {value}"))),
        }
    }

    fn finish(self) -> Result<(), SessionError> {
        if self.offset == self.bytes.len() {
            Ok(())
        } else {
            Err(self.corrupt(format!(
                "{} trailing byte(s) after the terminal envelope",
                self.bytes.len() - self.offset
            )))
        }
    }
}

fn encode_evidence(encoder: &mut Encoder, evidence: &RetainedEvidence) {
    encoder.string(&evidence.preview);
    encoder.u64(u64::try_from(evidence.byte_len).expect("bounded evidence length fits u64"));
    encoder.hash(evidence.digest);
}

fn decode_evidence(decoder: &mut Decoder<'_>) -> Result<RetainedEvidence, SessionError> {
    let preview = decoder.string(MAX_RETAINED_EVIDENCE_BYTES)?;
    let byte_len = usize::try_from(decoder.u64()?)
        .map_err(|_| decoder.corrupt("evidence length does not fit usize"))?;
    let digest = decoder.hash()?;
    if preview.len() > byte_len || preview.len() > MAX_RETAINED_EVIDENCE_BYTES {
        return Err(decoder.corrupt("evidence preview/length envelope is inconsistent"));
    }
    Ok(RetainedEvidence {
        preview,
        byte_len,
        digest,
    })
}

fn encode_snapshot(encoder: &mut Encoder, snapshot: MeterSnapshot) {
    encoder.u64(snapshot.core_s.to_bits());
    encoder.u64(snapshot.mem_peak_bytes);
    encoder.u64(snapshot.wall_s.to_bits());
    encoder.u32(snapshot.throttled);
    encoder.u32(snapshot.paused);
}

fn decode_snapshot(decoder: &mut Decoder<'_>) -> Result<MeterSnapshot, SessionError> {
    let snapshot = MeterSnapshot {
        core_s: f64::from_bits(decoder.u64()?),
        mem_peak_bytes: decoder.u64()?,
        wall_s: f64::from_bits(decoder.u64()?),
        throttled: decoder.u32()?,
        paused: decoder.u32()?,
    };
    validate_resource("recovered core-seconds meter", snapshot.core_s)?;
    validate_resource("recovered wall-seconds meter", snapshot.wall_s)?;
    Ok(snapshot)
}

fn encode_enforcement(encoder: &mut Encoder, enforcement: &Enforcement) {
    match enforcement {
        Enforcement::Ok => encoder.u8(0),
        Enforcement::Throttled {
            resource,
            used,
            granted,
        } => {
            encoder.u8(1);
            encoder.string(resource);
            encoder.u64(used.to_bits());
            encoder.u64(granted.to_bits());
        }
        Enforcement::Paused {
            resource,
            used,
            granted,
            resume_hint,
        } => {
            encoder.u8(2);
            encoder.string(resource);
            encoder.u64(used.to_bits());
            encoder.u64(granted.to_bits());
            encoder.string(resume_hint);
        }
    }
}

fn decode_resource(decoder: &mut Decoder<'_>) -> Result<&'static str, SessionError> {
    match decoder.string(32)?.as_str() {
        "core-seconds" => Ok("core-seconds"),
        "memory-bytes" => Ok("memory-bytes"),
        "wall-seconds" => Ok("wall-seconds"),
        other => Err(decoder.corrupt(format!("unknown enforcement resource {other:?}"))),
    }
}

fn decode_enforcement(decoder: &mut Decoder<'_>) -> Result<Enforcement, SessionError> {
    match decoder.u8()? {
        0 => Ok(Enforcement::Ok),
        1 => {
            let resource = decode_resource(decoder)?;
            let used = f64::from_bits(decoder.u64()?);
            let granted = f64::from_bits(decoder.u64()?);
            validate_resource("recovered enforcement used", used)?;
            validate_resource("recovered enforcement grant", granted)?;
            Ok(Enforcement::Throttled {
                resource,
                used,
                granted,
            })
        }
        2 => {
            let resource = decode_resource(decoder)?;
            let used = f64::from_bits(decoder.u64()?);
            let granted = f64::from_bits(decoder.u64()?);
            validate_resource("recovered enforcement used", used)?;
            validate_resource("recovered enforcement grant", granted)?;
            let resume_hint = decoder.string(MAX_CODEC_STRING_BYTES)?;
            Ok(Enforcement::Paused {
                resource,
                used,
                granted,
                resume_hint,
            })
        }
        value => Err(decoder.corrupt(format!("unknown enforcement tag {value}"))),
    }
}

fn encode_meter_receipt(encoder: &mut Encoder, receipt: &MeterReceipt) {
    encoder.hash(receipt.report_id.content_hash);
    encoder.u64(receipt.commit_ordinal);
    encoder.u64(receipt.delta.core_s.to_bits());
    encoder.u64(receipt.delta.mem_peak_bytes);
    encoder.u64(receipt.delta.wall_s.to_bits());
    encode_snapshot(encoder, receipt.before);
    encode_snapshot(encoder, receipt.after);
    encode_enforcement(encoder, &receipt.enforcement);
    encoder.hash(receipt.content_hash);
}

fn decode_meter_receipt(
    decoder: &mut Decoder<'_>,
    report_id: MeterReportId,
) -> Result<MeterReceipt, SessionError> {
    let stored_report = decoder.hash()?;
    if stored_report != report_id.content_hash {
        return Err(decoder.corrupt("meter receipt names a different report authority"));
    }
    let commit_ordinal = decoder.u64()?;
    let delta = Charge {
        core_s: f64::from_bits(decoder.u64()?),
        mem_peak_bytes: decoder.u64()?,
        wall_s: f64::from_bits(decoder.u64()?),
    };
    validate_resource("recovered core-seconds charge", delta.core_s)?;
    validate_resource("recovered wall-seconds charge", delta.wall_s)?;
    let before = decode_snapshot(decoder)?;
    let after = decode_snapshot(decoder)?;
    let enforcement = decode_enforcement(decoder)?;
    let content_hash = decoder.hash()?;
    let receipt = MeterReceipt {
        report_id,
        commit_ordinal,
        delta,
        before,
        after,
        enforcement,
        content_hash,
    };
    if meter_receipt_hash(
        report_id,
        commit_ordinal,
        delta,
        before,
        after,
        &receipt.enforcement,
    ) != content_hash
    {
        return Err(decoder.corrupt("meter receipt content hash does not verify"));
    }
    Ok(receipt)
}

pub(super) fn encode_open_payload(
    token: &CapabilityToken,
    gate_binding: Option<fs_blake3::ContentHash>,
) -> Vec<u8> {
    let mut encoder = Encoder::default();
    encoder.schema();
    encoder.u64(token.session.0);
    encoder.u64(u64::try_from(token.ops.len()).expect("bounded operator count fits u64"));
    for operation in &token.ops {
        encoder.string(operation);
    }
    encoder.u64(token.core_s.to_bits());
    encoder.u64(token.mem_bytes);
    encoder.u64(token.wall_s.to_bits());
    encoder.u64(token.cores);
    encoder.string(&token.ledger_scope);
    encoder.optional_hash(gate_binding);
    encoder.finish()
}

fn decode_open_payload(
    bytes: &[u8],
    authority: fs_blake3::ContentHash,
) -> Result<(CapabilityToken, Option<fs_blake3::ContentHash>), SessionError> {
    let mut decoder = Decoder::new(bytes, KIND_OPEN, authority)?;
    let session = SessionId(decoder.u64()?);
    let op_count = usize::try_from(decoder.u64()?)
        .map_err(|_| decoder.corrupt("operator count does not fit usize"))?;
    if op_count > MAX_CAPABILITY_OPS {
        return Err(decoder.corrupt(format!(
            "operator count {op_count} exceeds the {MAX_CAPABILITY_OPS}-item limit"
        )));
    }
    let mut ops = Vec::with_capacity(op_count);
    for _ in 0..op_count {
        ops.push(decoder.string(MAX_CAPABILITY_OP_BYTES)?);
    }
    let token = CapabilityToken {
        session,
        ops,
        core_s: f64::from_bits(decoder.u64()?),
        mem_bytes: decoder.u64()?,
        wall_s: f64::from_bits(decoder.u64()?),
        cores: decoder.u64()?,
        ledger_scope: decoder.string(MAX_LEDGER_SCOPE_BYTES)?,
    };
    let gate_binding = decoder.optional_hash()?;
    decoder.finish()?;
    token.validate_operator_grants()?;
    CapabilityToken::validate_ledger_scope(&token.ledger_scope)?;
    validate_resource("recovered core-seconds grant", token.core_s)?;
    validate_resource("recovered wall-seconds grant", token.wall_s)?;
    Ok((token, gate_binding))
}

pub(super) fn encode_open_receipt(receipt: &SessionOpenReceipt) -> Vec<u8> {
    let mut encoder = Encoder::default();
    encoder.schema();
    encoder.hash(receipt.open_id.content_hash);
    encoder.hash(receipt.token_digest);
    encoder.optional_hash(receipt.gate_identity);
    encoder.string(&receipt.permit.ledger_scope);
    encoder.hash(receipt.content_hash);
    encoder.finish()
}

fn decode_open_receipt(
    bytes: &[u8],
    governor_id: fs_blake3::ContentHash,
    open_id: SessionOpenId,
) -> Result<SessionOpenReceipt, SessionError> {
    let mut decoder = Decoder::new(bytes, KIND_OPEN, open_id.content_hash)?;
    if decoder.hash()? != open_id.content_hash {
        return Err(decoder.corrupt("open receipt names a different authority"));
    }
    let token_digest = decoder.hash()?;
    let gate_identity = decoder.optional_hash()?;
    let ledger_scope = decoder.string(MAX_LEDGER_SCOPE_BYTES)?;
    let content_hash = decoder.hash()?;
    decoder.finish()?;
    let receipt = SessionOpenReceipt {
        open_id,
        token_digest,
        gate_identity,
        permit: ScopeFlushPermit {
            governor_id,
            ledger_scope: ledger_scope.clone(),
        },
        content_hash,
    };
    if session_open_receipt_hash(open_id, token_digest, gate_identity, &ledger_scope)
        != content_hash
    {
        return Err(SessionError::TerminalCorrupt {
            kind: KIND_OPEN,
            authority: open_id.content_hash,
            detail: "open receipt content hash does not verify".to_string(),
        });
    }
    Ok(receipt)
}

pub(super) fn encode_meter_payload(delta: Charge) -> Vec<u8> {
    let mut encoder = Encoder::default();
    encoder.schema();
    encoder.u64(delta.core_s.to_bits());
    encoder.u64(delta.mem_peak_bytes);
    encoder.u64(delta.wall_s.to_bits());
    encoder.finish()
}

fn decode_meter_payload(
    bytes: &[u8],
    authority: fs_blake3::ContentHash,
) -> Result<Charge, SessionError> {
    let mut decoder = Decoder::new(bytes, KIND_METER, authority)?;
    let delta = Charge {
        core_s: f64::from_bits(decoder.u64()?),
        mem_peak_bytes: decoder.u64()?,
        wall_s: f64::from_bits(decoder.u64()?),
    };
    decoder.finish()?;
    validate_resource("recovered core-seconds charge", delta.core_s)?;
    validate_resource("recovered wall-seconds charge", delta.wall_s)?;
    Ok(delta)
}

pub(super) fn encode_meter_terminal_receipt(receipt: &MeterReceipt) -> Vec<u8> {
    let mut encoder = Encoder::default();
    encoder.schema();
    encode_meter_receipt(&mut encoder, receipt);
    encoder.finish()
}

fn decode_meter_terminal_receipt(
    bytes: &[u8],
    report_id: MeterReportId,
) -> Result<MeterReceipt, SessionError> {
    let mut decoder = Decoder::new(bytes, KIND_METER, report_id.content_hash)?;
    let receipt = decode_meter_receipt(&mut decoder, report_id)?;
    decoder.finish()?;
    Ok(receipt)
}

fn encode_pause_request(encoder: &mut Encoder, request: PauseRequestId) {
    encoder.hash(request.governor_id);
    encoder.u64(request.session.0);
    encoder.u64(request.gate_generation);
    encoder.i64(request.requested_ordinal);
}

fn decode_pause_request(decoder: &mut Decoder<'_>) -> Result<PauseRequestId, SessionError> {
    Ok(PauseRequestId {
        governor_id: decoder.hash()?,
        session: SessionId(decoder.u64()?),
        gate_generation: decoder.u64()?,
        requested_ordinal: decoder.i64()?,
    })
}

fn encode_pressure_action(encoder: &mut Encoder, action: PressureActionId) {
    encoder.hash(action.governor_id);
    encoder.u64(action.session.0);
    encoder.hash(action.session_open);
    encoder.u64(action.generation);
    encoder.hash(action.content_hash);
}

fn decode_pressure_action(decoder: &mut Decoder<'_>) -> Result<PressureActionId, SessionError> {
    Ok(PressureActionId {
        governor_id: decoder.hash()?,
        session: SessionId(decoder.u64()?),
        session_open: decoder.hash()?,
        generation: decoder.u64()?,
        content_hash: decoder.hash()?,
    })
}

fn encode_degradation_event(encoder: &mut Encoder, event: &DegradationEvent) {
    encoder.u64(event.session.0);
    encoder.u8(match event.step {
        DegradationStep::SpillColdArenas => 0,
        DegradationStep::CoarsenAdaptively => 1,
        DegradationStep::PauseSerializeResume => 2,
    });
    encoder.u8(event.pressure_level);
    encoder.u8(match event.phase {
        StepPhase::Declared => 0,
        StepPhase::Requested => 1,
        StepPhase::Complete => 2,
    });
    encoder.string(&event.attribution);
    encoder.i64(event.ordinal);
    encoder.optional_i64(event.requested_ordinal);
    match &event.checkpoint {
        Some(evidence) => {
            encoder.u8(1);
            encode_evidence(encoder, evidence);
        }
        None => encoder.u8(0),
    }
    encoder.optional_u64(event.gate_generation);
    match event.pause_request_id {
        Some(request) => {
            encoder.u8(1);
            encode_pause_request(encoder, request);
        }
        None => encoder.u8(0),
    }
    match event.pressure_action_id {
        Some(action) => {
            encoder.u8(1);
            encode_pressure_action(encoder, action);
        }
        None => encoder.u8(0),
    }
}

fn decode_degradation_event(decoder: &mut Decoder<'_>) -> Result<DegradationEvent, SessionError> {
    let session = SessionId(decoder.u64()?);
    let step = match decoder.u8()? {
        0 => DegradationStep::SpillColdArenas,
        1 => DegradationStep::CoarsenAdaptively,
        2 => DegradationStep::PauseSerializeResume,
        value => return Err(decoder.corrupt(format!("unknown degradation-step tag {value}"))),
    };
    let pressure_level = decoder.u8()?;
    if !(1..=3).contains(&pressure_level) {
        return Err(decoder.corrupt(format!("pressure level {pressure_level} is outside 1..=3")));
    }
    let phase = match decoder.u8()? {
        0 => StepPhase::Declared,
        1 => StepPhase::Requested,
        2 => StepPhase::Complete,
        value => return Err(decoder.corrupt(format!("unknown degradation-phase tag {value}"))),
    };
    let attribution = decoder.string(MAX_RETAINED_EVIDENCE_BYTES)?;
    let ordinal = decoder.i64()?;
    let requested_ordinal = decoder.optional_i64()?;
    let checkpoint = match decoder.u8()? {
        0 => None,
        1 => Some(decode_evidence(decoder)?),
        value => return Err(decoder.corrupt(format!("invalid checkpoint tag {value}"))),
    };
    let gate_generation = decoder.optional_u64()?;
    let pause_request_id = match decoder.u8()? {
        0 => None,
        1 => Some(decode_pause_request(decoder)?),
        value => return Err(decoder.corrupt(format!("invalid pause-request tag {value}"))),
    };
    let pressure_action_id = match decoder.u8()? {
        0 => None,
        1 => Some(decode_pressure_action(decoder)?),
        value => return Err(decoder.corrupt(format!("invalid pressure-action tag {value}"))),
    };
    Ok(DegradationEvent {
        session,
        step,
        pressure_level,
        phase,
        attribution,
        ordinal,
        requested_ordinal,
        checkpoint,
        gate_generation,
        pause_request_id,
        pressure_action_id,
    })
}

pub(super) fn encode_pressure_payload(level: u8) -> Vec<u8> {
    let mut encoder = Encoder::default();
    encoder.schema();
    encoder.u8(level);
    encoder.finish()
}

fn decode_pressure_payload(
    bytes: &[u8],
    authority: fs_blake3::ContentHash,
) -> Result<u8, SessionError> {
    let mut decoder = Decoder::new(bytes, KIND_PRESSURE, authority)?;
    let level = decoder.u8()?;
    decoder.finish()?;
    if !(1..=3).contains(&level) {
        return Err(SessionError::TerminalCorrupt {
            kind: KIND_PRESSURE,
            authority,
            detail: format!("pressure level {level} is outside 1..=3"),
        });
    }
    Ok(level)
}

pub(super) fn encode_pressure_terminal_receipt(receipt: &PressureReceipt) -> Vec<u8> {
    let mut encoder = Encoder::default();
    encoder.schema();
    encode_pressure_action(&mut encoder, receipt.action_id);
    encoder.u8(receipt.level);
    encoder.u64(u64::try_from(receipt.events.len()).expect("bounded event count fits u64"));
    for event in &receipt.events {
        encode_degradation_event(&mut encoder, event);
    }
    encoder.hash(receipt.content_hash);
    encoder.finish()
}

fn decode_pressure_terminal_receipt(
    bytes: &[u8],
    action_id: PressureActionId,
) -> Result<PressureReceipt, SessionError> {
    let mut decoder = Decoder::new(bytes, KIND_PRESSURE, action_id.content_hash)?;
    if decode_pressure_action(&mut decoder)? != action_id {
        return Err(decoder.corrupt("pressure receipt names a different action authority"));
    }
    let level = decoder.u8()?;
    let event_count = usize::try_from(decoder.u64()?)
        .map_err(|_| decoder.corrupt("pressure event count does not fit usize"))?;
    if event_count == 0 || event_count > LADDER.len() {
        return Err(decoder.corrupt(format!(
            "pressure receipt event count {event_count} is outside 1..={} ",
            LADDER.len()
        )));
    }
    let mut events = Vec::with_capacity(event_count);
    for _ in 0..event_count {
        events.push(decode_degradation_event(&mut decoder)?);
    }
    let content_hash = decoder.hash()?;
    decoder.finish()?;
    let canonical_events = events.iter().enumerate().all(|(index, event)| {
        let expected_step = LADDER[index];
        let is_pause = expected_step == DegradationStep::PauseSerializeResume;
        let expected_phase = if is_pause {
            StepPhase::Requested
        } else {
            StepPhase::Declared
        };
        let expected_request = is_pause.then_some(PauseRequestId {
            governor_id: action_id.governor_id,
            session: action_id.session,
            gate_generation: action_id.generation,
            requested_ordinal: event.ordinal,
        });
        event.session == action_id.session
            && event.step == expected_step
            && event.pressure_level == level
            && event.phase == expected_phase
            && event.attribution == degradation_attribution(expected_step)
            && event.requested_ordinal.is_none()
            && event.checkpoint.is_none()
            && event.gate_generation == is_pause.then_some(action_id.generation)
            && event.pause_request_id == expected_request
            && event.pressure_action_id == Some(action_id)
    });
    let dense_ordinals = events.windows(2).all(|pair| {
        pair[0]
            .ordinal
            .checked_add(1)
            .is_some_and(|next| next == pair[1].ordinal)
    });
    if level != u8::try_from(event_count).expect("ladder length fits u8")
        || !canonical_events
        || !dense_ordinals
        || pressure_receipt_hash(action_id, level, &events) != content_hash
    {
        return Err(SessionError::TerminalCorrupt {
            kind: KIND_PRESSURE,
            authority: action_id.content_hash,
            detail:
                "pressure receipt authority, level, event group, or content hash does not verify"
                    .to_string(),
        });
    }
    Ok(PressureReceipt {
        action_id,
        level,
        events,
        content_hash,
    })
}

pub(super) fn encode_submission_payload(request_id: SubmissionRequestId) -> Vec<u8> {
    let mut encoder = Encoder::default();
    encoder.schema();
    encoder.hash(request_id.key_hash);
    encoder.hash(request_id.request_hash);
    encoder.finish()
}

fn decode_submission_payload(
    bytes: &[u8],
    request_id: SubmissionRequestId,
) -> Result<(), SessionError> {
    let mut decoder = Decoder::new(bytes, KIND_SUBMISSION, request_id.content_hash)?;
    let key_hash = decoder.hash()?;
    let request_hash = decoder.hash()?;
    decoder.finish()?;
    if key_hash != request_id.key_hash || request_hash != request_id.request_hash {
        return Err(SessionError::MutationConflict {
            kind: KIND_SUBMISSION,
            id: request_id.content_hash,
        });
    }
    Ok(())
}

pub(super) fn encode_submission_terminal_receipt(
    state: &IdemState,
) -> Result<Vec<u8>, SessionError> {
    let mut encoder = Encoder::default();
    encoder.schema();
    match state {
        IdemState::Done {
            admission_ordinal,
            receipt,
            charge,
            meter_receipt,
            ..
        } => {
            encoder.u8(0);
            encoder.u64(*admission_ordinal);
            encoder.hash(receipt.0);
            encoder.u64(charge.core_s.to_bits());
            encoder.u64(charge.mem_peak_bytes);
            encoder.u64(charge.wall_s.to_bits());
            encode_meter_receipt(&mut encoder, meter_receipt);
        }
        IdemState::Failed {
            admission_ordinal,
            receipt,
            evidence,
            ..
        } => {
            encoder.u8(1);
            encoder.u64(*admission_ordinal);
            encoder.hash(receipt.0);
            encode_evidence(&mut encoder, evidence);
        }
        IdemState::Pending { .. } => {
            return Err(SessionError::Persistence {
                what: "pending submission has no terminal receipt to encode".to_string(),
            });
        }
    }
    Ok(encoder.finish())
}

enum RecoveredSubmission {
    Done {
        admission_ordinal: u64,
        receipt: SubmissionReceipt,
        charge: Charge,
        meter_receipt: MeterReceipt,
    },
    Failed {
        admission_ordinal: u64,
        receipt: SubmissionReceipt,
        evidence: RetainedEvidence,
    },
}

fn decode_submission_terminal_receipt(
    bytes: &[u8],
    request_id: SubmissionRequestId,
    ledger_scope: &str,
) -> Result<RecoveredSubmission, SessionError> {
    let mut decoder = Decoder::new(bytes, KIND_SUBMISSION, request_id.content_hash)?;
    let tag = decoder.u8()?;
    let admission_ordinal = decoder.u64()?;
    if admission_ordinal == 0 || admission_ordinal > i64::MAX as u64 {
        return Err(decoder.corrupt(format!(
            "submission admission ordinal {admission_ordinal} is outside 1..=i64::MAX"
        )));
    }
    let receipt = SubmissionReceipt(decoder.hash()?);
    let recovered = match tag {
        0 => {
            let charge = Charge {
                core_s: f64::from_bits(decoder.u64()?),
                mem_peak_bytes: decoder.u64()?,
                wall_s: f64::from_bits(decoder.u64()?),
            };
            validate_resource("recovered submission core-seconds charge", charge.core_s)?;
            validate_resource("recovered submission wall-seconds charge", charge.wall_s)?;
            let report_id = Governor::submission_meter_report_id(request_id);
            let meter_receipt = decode_meter_receipt(&mut decoder, report_id)?;
            if !same_charge(charge, meter_receipt.delta)
                || !receipt.matches_success(
                    request_id,
                    ledger_scope,
                    admission_ordinal,
                    charge,
                    &meter_receipt,
                )
            {
                return Err(decoder.corrupt(
                    "successful submission receipt, charge, or causal meter receipt does not verify",
                ));
            }
            RecoveredSubmission::Done {
                admission_ordinal,
                receipt,
                charge,
                meter_receipt,
            }
        }
        1 => {
            let evidence = decode_evidence(&mut decoder)?;
            if !receipt.matches_failure(request_id, ledger_scope, admission_ordinal, &evidence) {
                return Err(decoder.corrupt("failed submission receipt does not verify"));
            }
            RecoveredSubmission::Failed {
                admission_ordinal,
                receipt,
                evidence,
            }
        }
        value => return Err(decoder.corrupt(format!("unknown submission terminal tag {value}"))),
    };
    decoder.finish()?;
    Ok(recovered)
}

pub(super) fn encode_pause_ack_payload(evidence: &RetainedEvidence) -> Vec<u8> {
    let mut encoder = Encoder::default();
    encoder.schema();
    encode_evidence(&mut encoder, evidence);
    encoder.finish()
}

fn decode_pause_ack_payload(
    bytes: &[u8],
    request_id: PauseRequestId,
) -> Result<RetainedEvidence, SessionError> {
    let authority = pause_ack_authority(request_id);
    let mut decoder = Decoder::new(bytes, KIND_PAUSE_ACK, authority)?;
    let evidence = decode_evidence(&mut decoder)?;
    decoder.finish()?;
    Ok(evidence)
}

pub(super) fn pause_ack_authority(request_id: PauseRequestId) -> fs_blake3::ContentHash {
    let mut encoder = Encoder::default();
    encoder.hash(request_id.governor_id);
    encoder.u64(request_id.session.0);
    encoder.u64(request_id.gate_generation);
    encoder.i64(request_id.requested_ordinal);
    fs_blake3::hash_domain(PAUSE_ACKNOWLEDGEMENT_ID_DOMAIN, &encoder.finish())
}

pub(super) fn encode_pause_ack_terminal_receipt(acknowledgement: &PauseAcknowledgement) -> Vec<u8> {
    let mut encoder = Encoder::default();
    encoder.schema();
    encode_pause_request(&mut encoder, acknowledgement.request_id);
    encode_degradation_event(&mut encoder, &acknowledgement.event);
    encoder.u64(acknowledgement.resume_generation);
    encoder.hash(acknowledgement.gate_binding);
    encoder.hash(acknowledgement.content_hash);
    encoder.finish()
}

fn decode_pause_ack_terminal_receipt(
    bytes: &[u8],
    request_id: PauseRequestId,
    resume_gate: Arc<CancelGate>,
) -> Result<PauseAcknowledgement, SessionError> {
    let authority = pause_ack_authority(request_id);
    let mut decoder = Decoder::new(bytes, KIND_PAUSE_ACK, authority)?;
    if decode_pause_request(&mut decoder)? != request_id {
        return Err(decoder.corrupt("pause acknowledgement names a different request"));
    }
    let event = decode_degradation_event(&mut decoder)?;
    let resume_generation = decoder.u64()?;
    let gate_binding = decoder.hash()?;
    let content_hash = decoder.hash()?;
    decoder.finish()?;
    let expected_resume_generation = request_id.gate_generation.checked_add(1);
    let expected_attribution = event.checkpoint.as_ref().map(|evidence| {
        format!(
            "pause complete: checkpoint evidence ({} bytes, digest {}) acknowledges \
             the request at ordinal {} and rotates gate generation {} to {resume_generation}",
            evidence.byte_len,
            evidence.digest,
            request_id.requested_ordinal,
            request_id.gate_generation,
        )
    });
    let canonical_action = event.pressure_action_id.is_some_and(|action| {
        action.governor_id == request_id.governor_id
            && action.session == request_id.session
            && action.generation == request_id.gate_generation
    });
    if event.session != request_id.session
        || event.step != DegradationStep::PauseSerializeResume
        || event.pressure_level != 3
        || event.pause_request_id != Some(request_id)
        || event.phase != StepPhase::Complete
        || event.checkpoint.is_none()
        || event.attribution != expected_attribution.unwrap_or_default()
        || event.ordinal <= request_id.requested_ordinal
        || event.requested_ordinal != Some(request_id.requested_ordinal)
        || event.gate_generation != Some(request_id.gate_generation)
        || !canonical_action
        || expected_resume_generation != Some(resume_generation)
        || resumed_gate_binding(request_id, resume_generation) != gate_binding
        || pause_acknowledgement_hash(request_id, &event, resume_generation, gate_binding)
            != content_hash
    {
        return Err(SessionError::TerminalCorrupt {
            kind: KIND_PAUSE_ACK,
            authority,
            detail:
                "pause acknowledgement event, generation, binding, or receipt hash does not verify"
                    .to_string(),
        });
    }
    Ok(PauseAcknowledgement {
        request_id,
        event,
        resume_gate,
        resume_generation,
        gate_binding,
        content_hash,
    })
}

pub(super) fn encode_activation_payload(acknowledgement: &PauseAcknowledgement) -> Vec<u8> {
    encode_activation_payload_parts(
        acknowledgement.content_hash,
        acknowledgement.resume_generation,
        acknowledgement.gate_binding,
    )
}

pub(super) fn encode_activation_payload_parts(
    acknowledgement_hash: fs_blake3::ContentHash,
    resume_generation: u64,
    gate_binding: fs_blake3::ContentHash,
) -> Vec<u8> {
    let mut encoder = Encoder::default();
    encoder.schema();
    encoder.hash(acknowledgement_hash);
    encoder.u64(resume_generation);
    encoder.hash(gate_binding);
    encoder.finish()
}

fn decode_activation_payload(
    bytes: &[u8],
    activation_id: ResumeActivationId,
) -> Result<(fs_blake3::ContentHash, u64, fs_blake3::ContentHash), SessionError> {
    let mut decoder = Decoder::new(bytes, KIND_RESUME_ACTIVATION, activation_id.content_hash)?;
    let acknowledgement_hash = decoder.hash()?;
    let resume_generation = decoder.u64()?;
    let gate_binding = decoder.hash()?;
    decoder.finish()?;
    Ok((acknowledgement_hash, resume_generation, gate_binding))
}

pub(super) fn encode_activation_terminal_receipt(receipt: ResumeActivationReceipt) -> Vec<u8> {
    let mut encoder = Encoder::default();
    encoder.schema();
    encoder.hash(receipt.activation_id.content_hash);
    encoder.hash(receipt.acknowledgement_hash);
    encoder.hash(receipt.gate_binding);
    encoder.hash(receipt.content_hash);
    encoder.finish()
}

fn decode_activation_terminal_receipt(
    bytes: &[u8],
    activation_id: ResumeActivationId,
) -> Result<ResumeActivationReceipt, SessionError> {
    let mut decoder = Decoder::new(bytes, KIND_RESUME_ACTIVATION, activation_id.content_hash)?;
    if decoder.hash()? != activation_id.content_hash {
        return Err(decoder.corrupt("activation receipt names a different authority"));
    }
    let acknowledgement_hash = decoder.hash()?;
    let gate_binding = decoder.hash()?;
    let content_hash = decoder.hash()?;
    decoder.finish()?;
    let receipt = resume_activation_receipt(activation_id, acknowledgement_hash, gate_binding);
    if receipt.content_hash != content_hash {
        return Err(SessionError::TerminalCorrupt {
            kind: KIND_RESUME_ACTIVATION,
            authority: activation_id.content_hash,
            detail: "activation receipt content hash does not verify".to_string(),
        });
    }
    Ok(receipt)
}

fn ledger_error(context: &'static str, error: fs_ledger::LedgerError) -> SessionError {
    SessionError::Persistence {
        what: format!("{context}: {error}"),
    }
}

fn same_snapshot(left: MeterSnapshot, right: MeterSnapshot) -> bool {
    left.core_s.to_bits() == right.core_s.to_bits()
        && left.mem_peak_bytes == right.mem_peak_bytes
        && left.wall_s.to_bits() == right.wall_s.to_bits()
        && left.throttled == right.throttled
        && left.paused == right.paused
}

fn cached_pause_acknowledgement(
    inner: &Inner,
    request_id: PauseRequestId,
    expected: &PauseAcknowledgement,
    resume_gate: &Arc<CancelGate>,
) -> Result<Option<PauseAcknowledgement>, SessionError> {
    let Some(replay) = inner.pause_acknowledgements.get(&request_id) else {
        return Ok(None);
    };
    let completed = inner.completed_pause.get(&request_id.session.0).copied();
    let token = inner
        .tokens
        .get(&request_id.session.0)
        .ok_or(SessionError::UnknownSession {
            id: request_id.session.0,
        })?;
    let event = inner
        .scopes
        .get(&token.ledger_scope)
        .and_then(|scope| scope.events.get(replay.completion_event_index))
        .ok_or_else(|| SessionError::TerminalCorrupt {
            kind: KIND_PAUSE_ACK,
            authority: pause_ack_authority(request_id),
            detail: "recovered pause acknowledgement lost its completion event".to_string(),
        })?;
    let current_gate =
        inner
            .gates
            .get(&request_id.session.0)
            .cloned()
            .ok_or(SessionError::UngatedSession {
                id: request_id.session.0,
            })?;
    if event.as_ref() != &expected.event
        || !completed.is_some_and(|completed| {
            completed.request_id == request_id
                && completed.resume_generation == replay.resume_generation
                && completed.gate_binding == replay.gate_binding
                && completed.acknowledgement_hash == replay.content_hash
        })
        || inner.gate_generations.get(&request_id.session.0) != Some(&replay.resume_generation)
        || replay.resume_generation != expected.resume_generation
        || replay.gate_binding != expected.gate_binding
        || replay.content_hash != expected.content_hash
        || !Arc::ptr_eq(&current_gate, resume_gate)
    {
        return Err(SessionError::PauseAcknowledgementConflict {
            id: request_id.session.0,
            requested_ordinal: request_id.requested_ordinal,
        });
    }
    Ok(Some(PauseAcknowledgement {
        request_id,
        event: event.as_ref().clone(),
        resume_gate: current_gate,
        resume_generation: replay.resume_generation,
        gate_binding: replay.gate_binding,
        content_hash: replay.content_hash,
    }))
}

fn cached_resume_activation(
    inner: &Inner,
    activation_id: ResumeActivationId,
    acknowledgement: &PauseAcknowledgement,
) -> Result<Option<ResumeActivationReceipt>, SessionError> {
    let Some(receipt) = inner.resume_activations.get(&activation_id) else {
        return Ok(None);
    };
    let request_id = acknowledgement.request_id;
    let completed = inner
        .completed_pause
        .get(&request_id.session.0)
        .copied()
        .ok_or(SessionError::ResumeAcknowledgementMismatch {
            id: request_id.session.0,
        })?;
    let current_gate =
        inner
            .gates
            .get(&request_id.session.0)
            .ok_or(SessionError::UngatedSession {
                id: request_id.session.0,
            })?;
    if completed.request_id != request_id
        || completed.acknowledgement_hash != acknowledgement.content_hash
        || completed.gate_binding != acknowledgement.gate_binding
        || inner.gate_phases.get(&request_id.session.0) != Some(&GatePhase::Running)
        || !Arc::ptr_eq(current_gate, &acknowledgement.resume_gate)
        || receipt.acknowledgement_hash != acknowledgement.content_hash
        || receipt.gate_binding != acknowledgement.gate_binding
    {
        return Err(SessionError::ResumeAcknowledgementMismatch {
            id: request_id.session.0,
        });
    }
    Ok(Some(*receipt))
}

impl Governor {
    pub(crate) fn recovery_ledger(
        &self,
        ledger: &fs_ledger::Ledger,
    ) -> Result<fs_ledger::LedgerInstanceId, SessionError> {
        let expected = self.durable_sink.ok_or_else(|| SessionError::Persistence {
            what: "durable recovery requires Governor::new_durable".to_string(),
        })?;
        let attempted = ledger
            .checked_instance_id()
            .map_err(|error| ledger_error("recovery ledger identity validation", error))?;
        if attempted != expected {
            return Err(SessionError::RecoveryLedgerMismatch {
                expected,
                attempted,
            });
        }
        Ok(attempted)
    }

    #[allow(clippy::too_many_arguments)] // Exact expected durable claim envelope.
    pub(crate) fn validate_recovered_terminal(
        &self,
        terminal: &fs_ledger::session_registry::StoredSessionTerminal,
        ledger_instance_id: fs_ledger::LedgerInstanceId,
        kind: &'static str,
        authority: fs_blake3::ContentHash,
        session: SessionId,
        ledger_scope: &str,
        session_open: fs_blake3::ContentHash,
        generation: u64,
        causal_ordinal: Option<u64>,
        payload: &[u8],
        expected_events: &[EventRow<'_>],
    ) -> Result<(), SessionError> {
        let claim = &terminal.claim;
        if claim.governor_hash != self.id {
            return Err(SessionError::RecoveryGovernorMismatch {
                expected: self.id,
                found: claim.governor_hash,
            });
        }
        let structural_match = claim.authority == authority
            && claim.ledger_instance_id == ledger_instance_id
            && claim.kind == kind
            && claim.session == session.0
            && claim.ledger_scope == ledger_scope
            && claim.session_open_hash == session_open
            && claim.generation == generation
            && claim.causal_ordinal == causal_ordinal;
        if !structural_match {
            return Err(SessionError::TerminalCorrupt {
                kind,
                authority,
                detail: "stored claim authority, ledger, session, scope, open identity, generation, or causal ordinal does not match the typed request"
                    .to_string(),
            });
        }
        if claim.payload != payload {
            return Err(SessionError::MutationConflict {
                kind,
                id: authority,
            });
        }
        let events_hash =
            fs_ledger::session_registry::session_terminal_events_hash(expected_events);
        if terminal.event_count != expected_events.len() || terminal.events_hash != events_hash {
            return Err(SessionError::TerminalCorrupt {
                kind,
                authority,
                detail:
                    "typed receipt does not reproduce the authenticated ordered audit-event group"
                        .to_string(),
            });
        }
        Ok(())
    }

    /// Reconstruct and install one already-durable session open without
    /// dirtying a flush cursor or changing audit rows.
    ///
    /// The supplied gate is a fresh process-local binding for the persisted
    /// semantic gate identity. Repeating recovery in this governor requires
    /// the same [`Arc`] identity, just like ordinary open replay.
    ///
    /// # Errors
    /// Missing/corrupt/conflicting terminals, foreign ledgers/governors,
    /// altered tokens/gating, and ordinary open validation failures fail
    /// closed.
    #[allow(clippy::too_many_lines)] // Validation and no-dirty installation are one recovery transaction.
    pub fn recover_open(
        &self,
        ledger: &fs_ledger::Ledger,
        open_id: SessionOpenId,
        token: CapabilityToken,
        gate: Option<Arc<CancelGate>>,
    ) -> Result<SessionOpenReceipt, SessionError> {
        if open_id.governor_id != self.id || open_id.session != token.session {
            return Err(SessionError::MutationAuthorityMismatch {
                kind: KIND_OPEN,
                id: open_id.content_hash,
            });
        }
        let ledger_instance_id = self.recovery_ledger(ledger)?;
        let terminal = ledger
            .session_terminal(&open_id.content_hash)
            .map_err(|error| ledger_error("recover open terminal", error))?
            .ok_or(SessionError::RecoveryRequired {
                kind: KIND_OPEN,
                authority: open_id.content_hash,
            })?;
        let payload =
            encode_open_payload(&token, gate.as_ref().map(|_| session_gate_binding(open_id)));
        let (stored_token, stored_gate_binding) =
            decode_open_payload(&terminal.claim.payload, open_id.content_hash)?;
        let expected_gate_binding = gate.as_ref().map(|_| session_gate_binding(open_id));
        if stored_token != token || stored_gate_binding != expected_gate_binding {
            return Err(SessionError::MutationConflict {
                kind: KIND_OPEN,
                id: open_id.content_hash,
            });
        }
        let stored_receipt = decode_open_receipt(&terminal.receipt, self.id, open_id)?;
        let reconstructed_token_digest = capability_token_identity(&token);
        if stored_receipt.token_digest != reconstructed_token_digest
            || stored_receipt.gate_identity != expected_gate_binding
        {
            return Err(SessionError::TerminalCorrupt {
                kind: KIND_OPEN,
                authority: open_id.content_hash,
                detail: "open terminal token digest or gate binding does not verify".to_string(),
            });
        }
        let reconstructed_receipt = SessionOpenReceipt {
            open_id,
            token_digest: reconstructed_token_digest,
            gate_identity: expected_gate_binding,
            permit: ScopeFlushPermit {
                governor_id: self.id,
                ledger_scope: token.ledger_scope.clone(),
            },
            content_hash: session_open_receipt_hash(
                open_id,
                reconstructed_token_digest,
                expected_gate_binding,
                &token.ledger_scope,
            ),
        };
        if reconstructed_receipt != stored_receipt {
            return Err(SessionError::TerminalCorrupt {
                kind: KIND_OPEN,
                authority: open_id.content_hash,
                detail: "reconstructed open receipt differs from durable receipt".to_string(),
            });
        }
        let expected_event =
            buffered_open_receipt(&token.ledger_scope, open_id, &stored_receipt, &token);
        let expected_rows = [expected_event.as_row()];
        self.validate_recovered_terminal(
            &terminal,
            ledger_instance_id,
            KIND_OPEN,
            open_id.content_hash,
            token.session,
            &token.ledger_scope,
            stored_receipt.content_hash,
            0,
            None,
            &payload,
            &expected_rows,
        )?;

        let receipt = self.register_session(open_id, token, gate, true)?;
        if receipt != stored_receipt {
            return Err(SessionError::TerminalCorrupt {
                kind: KIND_OPEN,
                authority: open_id.content_hash,
                detail: "reconstructed open receipt differs from durable receipt".to_string(),
            });
        }
        let mut inner = self.inner.lock().expect("governor lock");
        let scope = inner
            .scopes
            .get_mut(receipt.flush_permit().ledger_scope())
            .expect("recovered open registered its scope");
        scope.dirty_open_receipts.remove(&open_id);
        scope.sink.get_or_insert(ledger_instance_id);
        self.mark_durable_claim_recovered(&mut inner, open_id.content_hash);
        Ok(receipt)
    }

    /// Reconstruct one contiguous durable meter commit exactly once.
    ///
    /// # Errors
    /// Missing/corrupt/conflicting terminals, causal gaps, foreign authority,
    /// capacity, and meter-transition mismatch fail closed.
    #[allow(clippy::too_many_lines)] // Full causal verification precedes one lock-held install.
    pub fn recover_meter(
        &self,
        ledger: &fs_ledger::Ledger,
        report_id: MeterReportId,
        delta: Charge,
    ) -> Result<MeterReceipt, SessionError> {
        if report_id.governor_id != self.id {
            return Err(SessionError::MutationAuthorityMismatch {
                kind: KIND_METER,
                id: report_id.content_hash,
            });
        }
        let ledger_instance_id = self.recovery_ledger(ledger)?;
        {
            let inner = self.inner.lock().expect("governor lock");
            if let Some(receipt) = inner.meter_reports.get(&report_id) {
                return if same_charge(receipt.delta, delta) {
                    Ok(receipt.clone())
                } else {
                    Err(SessionError::MutationConflict {
                        kind: KIND_METER,
                        id: report_id.content_hash,
                    })
                };
            }
        }
        let terminal = ledger
            .session_terminal(&report_id.content_hash)
            .map_err(|error| ledger_error("recover meter terminal", error))?
            .ok_or(SessionError::RecoveryRequired {
                kind: KIND_METER,
                authority: report_id.content_hash,
            })?;
        let payload = encode_meter_payload(delta);
        if !same_charge(
            decode_meter_payload(&terminal.claim.payload, report_id.content_hash)?,
            delta,
        ) {
            return Err(SessionError::MutationConflict {
                kind: KIND_METER,
                id: report_id.content_hash,
            });
        }
        let receipt = decode_meter_terminal_receipt(&terminal.receipt, report_id)?;
        let ledger_scope = terminal.claim.ledger_scope.clone();
        let expected_event = buffered_meter_receipt(&ledger_scope, report_id, &receipt)?;
        let expected_rows = [expected_event.as_row()];
        self.validate_recovered_terminal(
            &terminal,
            ledger_instance_id,
            KIND_METER,
            report_id.content_hash,
            report_id.session,
            &ledger_scope,
            report_id.session_open,
            report_id.generation,
            Some(receipt.commit_ordinal),
            &payload,
            &expected_rows,
        )?;

        let mut inner = self.inner.lock().expect("governor lock");
        if let Some(existing) = inner.meter_reports.get(&report_id) {
            return if existing == &receipt && same_charge(existing.delta, delta) {
                Ok(existing.clone())
            } else {
                Err(SessionError::MutationConflict {
                    kind: KIND_METER,
                    id: report_id.content_hash,
                })
            };
        }
        let token = inner.tokens.get(&report_id.session.0).cloned().ok_or(
            SessionError::UnknownSession {
                id: report_id.session.0,
            },
        )?;
        if token.ledger_scope != ledger_scope
            || Self::current_open_identity(&inner, report_id.session)? != report_id.session_open
        {
            return Err(SessionError::MutationAuthorityMismatch {
                kind: KIND_METER,
                id: report_id.content_hash,
            });
        }
        let current_generation = inner
            .gate_generations
            .get(&report_id.session.0)
            .copied()
            .unwrap_or(0);
        if report_id.generation > current_generation {
            return Err(SessionError::StaleMutationGeneration {
                kind: KIND_METER,
                id: report_id.session.0,
                supplied: report_id.generation,
                current: current_generation,
            });
        }
        let expected_ordinal =
            inner
                .next_meter_commit_ordinal
                .checked_add(1)
                .ok_or(SessionError::LimitExceeded {
                    resource: "meter_commit_ordinal",
                    limit: i64::MAX as usize,
                    observed_at_least: usize::MAX,
                })?;
        if receipt.commit_ordinal != expected_ordinal {
            return Err(SessionError::RecoveryCausalGap {
                session: report_id.session.0,
                expected: expected_ordinal,
                found: receipt.commit_ordinal,
            });
        }
        let before = inner
            .meters
            .get(&report_id.session.0)
            .cloned()
            .unwrap_or_default();
        if !same_snapshot(before.snapshot(), receipt.before) {
            return Err(SessionError::TerminalCorrupt {
                kind: KIND_METER,
                authority: report_id.content_hash,
                detail: "meter receipt before-state does not continue the recovered meter"
                    .to_string(),
            });
        }
        let (after, enforcement) = meter_transition(&token, &before, delta)?;
        if !same_snapshot(after.snapshot(), receipt.after) || enforcement != receipt.enforcement {
            return Err(SessionError::TerminalCorrupt {
                kind: KIND_METER,
                authority: report_id.content_hash,
                detail: "meter receipt after-state or enforcement does not recompute".to_string(),
            });
        }
        let reports = inner
            .meter_report_ids
            .get(&report_id.session.0)
            .expect("recovered open created meter index");
        if reports.len() >= MAX_METER_REPORTS_PER_SESSION {
            return Err(SessionError::LimitExceeded {
                resource: "meter_reports_per_session",
                limit: MAX_METER_REPORTS_PER_SESSION,
                observed_at_least: reports.len().saturating_add(1),
            });
        }
        ensure_retained_capacity(&inner, &ledger_scope, MAX_METER_RECEIPT_RETAINED_BYTES)?;
        inner.next_meter_commit_ordinal = receipt.commit_ordinal;
        inner.meters.insert(report_id.session.0, after);
        inner
            .meter_report_ids
            .get_mut(&report_id.session.0)
            .expect("recovered open created meter index")
            .insert(report_id);
        inner.meter_reports.insert(report_id, receipt.clone());
        inner
            .scopes
            .get_mut(&ledger_scope)
            .expect("recovered scope")
            .sink
            .get_or_insert(ledger_instance_id);
        commit_retained_bytes(&mut inner, &ledger_scope, MAX_METER_RECEIPT_RETAINED_BYTES);
        self.mark_durable_claim_recovered(&mut inner, report_id.content_hash);
        Ok(receipt)
    }

    /// Reconstruct one durable submission terminal without accepting or
    /// invoking caller work. Missing terminals are Indeterminate by design.
    ///
    /// # Errors
    /// Pending/missing, corrupt, conflicting, causally discontinuous, foreign,
    /// and over-cap terminals fail closed; no closure is accepted by this API.
    #[allow(clippy::too_many_lines)] // Terminal decode, causal meter install, and idempotency install are indivisible.
    pub fn recover_submission(
        &self,
        ledger: &fs_ledger::Ledger,
        request_id: SubmissionRequestId,
        canonical_program: &str,
    ) -> Result<SubmitOutcome, SessionError> {
        if request_id.governor_id != self.id {
            return Err(SessionError::MutationAuthorityMismatch {
                kind: KIND_SUBMISSION,
                id: request_id.content_hash,
            });
        }
        let supplied_request_hash = bounded_request_digest(
            "idempotency_program_text_bytes",
            IDEMPOTENCY_PROGRAM_DOMAIN,
            canonical_program,
        )?;
        if supplied_request_hash != request_id.request_hash {
            return Err(SessionError::MutationConflict {
                kind: KIND_SUBMISSION,
                id: request_id.content_hash,
            });
        }
        let ledger_instance_id = self.recovery_ledger(ledger)?;
        {
            let inner = self.inner.lock().expect("governor lock");
            match inner.idempotency.get(&request_id) {
                Some(IdemState::Done {
                    admission_ordinal,
                    receipt,
                    meter_receipt,
                    ..
                }) => {
                    return Ok(SubmitOutcome::Duplicate {
                        admission_ordinal: *admission_ordinal,
                        receipt: *receipt,
                        enforcement: meter_receipt.enforcement.clone(),
                        meter_receipt: meter_receipt.clone(),
                    });
                }
                Some(IdemState::Failed {
                    admission_ordinal,
                    receipt,
                    evidence,
                    ..
                }) => {
                    return Ok(SubmitOutcome::Failed {
                        admission_ordinal: *admission_ordinal,
                        receipt: *receipt,
                        evidence: evidence.as_ref().clone(),
                    });
                }
                Some(IdemState::Pending { .. }) => return Ok(SubmitOutcome::InFlight),
                None => {}
            }
        }

        let payload = encode_submission_payload(request_id);
        let terminal = match ledger
            .session_terminal(&request_id.content_hash)
            .map_err(|error| ledger_error("recover submission terminal", error))?
        {
            Some(terminal) => terminal,
            None => {
                let expected_scope = {
                    let inner = self.inner.lock().expect("governor lock");
                    let token = inner.tokens.get(&request_id.session.0).ok_or(
                        SessionError::UnknownSession {
                            id: request_id.session.0,
                        },
                    )?;
                    if Self::current_open_identity(&inner, request_id.session)?
                        != request_id.session_open
                    {
                        return Err(SessionError::MutationAuthorityMismatch {
                            kind: KIND_SUBMISSION,
                            id: request_id.content_hash,
                        });
                    }
                    let current_generation = inner
                        .gate_generations
                        .get(&request_id.session.0)
                        .copied()
                        .unwrap_or(0);
                    if request_id.generation > current_generation {
                        return Err(SessionError::StaleMutationGeneration {
                            kind: KIND_SUBMISSION,
                            id: request_id.session.0,
                            supplied: request_id.generation,
                            current: current_generation,
                        });
                    }
                    token.ledger_scope.clone()
                };
                if let Some(claim) = ledger
                    .session_mutation_claim(&request_id.content_hash)
                    .map_err(|error| ledger_error("recover submission claim", error))?
                {
                    if claim.payload != payload {
                        return Err(SessionError::MutationConflict {
                            kind: KIND_SUBMISSION,
                            id: request_id.content_hash,
                        });
                    }
                    if claim.authority != request_id.content_hash
                        || claim.ledger_instance_id != ledger_instance_id
                        || claim.governor_hash != self.id
                        || claim.session_open_hash != request_id.session_open
                        || claim.schema_version
                            != fs_ledger::session_registry::SESSION_REGISTRY_ROW_SCHEMA_VERSION
                        || claim.kind != KIND_SUBMISSION
                        || claim.session != request_id.session.0
                        || claim.ledger_scope != expected_scope
                        || claim.generation != request_id.generation
                    {
                        return Err(SessionError::TerminalCorrupt {
                            kind: KIND_SUBMISSION,
                            authority: request_id.content_hash,
                            detail:
                                "Pending submission claim envelope differs from its typed authority"
                                    .to_string(),
                        });
                    }
                }
                return Err(SessionError::IndeterminateMutation {
                    kind: KIND_SUBMISSION,
                    authority: request_id.content_hash,
                });
            }
        };
        decode_submission_payload(&terminal.claim.payload, request_id)?;
        let ledger_scope = terminal.claim.ledger_scope.clone();
        let recovered =
            decode_submission_terminal_receipt(&terminal.receipt, request_id, &ledger_scope)?;
        let (state, event) = match &recovered {
            RecoveredSubmission::Done {
                admission_ordinal,
                receipt,
                charge,
                meter_receipt,
            } => {
                let state = IdemState::Done {
                    admission_ordinal: *admission_ordinal,
                    receipt: *receipt,
                    charge: *charge,
                    meter_receipt: meter_receipt.clone(),
                    durable_permit: None,
                };
                let (event, _) = buffered_submission_success(&ledger_scope, request_id, &state)?;
                (state, event)
            }
            RecoveredSubmission::Failed {
                admission_ordinal,
                receipt,
                evidence,
            } => {
                let state = IdemState::Failed {
                    admission_ordinal: *admission_ordinal,
                    receipt: *receipt,
                    evidence: Arc::new(evidence.clone()),
                    durable_permit: None,
                };
                let (event, _) = buffered_submission_failure(&ledger_scope, request_id, &state)?;
                (state, event)
            }
        };
        let expected_rows = [event.as_row()];
        let admission_ordinal = match &recovered {
            RecoveredSubmission::Done {
                admission_ordinal, ..
            }
            | RecoveredSubmission::Failed {
                admission_ordinal, ..
            } => *admission_ordinal,
        };
        if let Some(claim_ordinal) = terminal.claim.causal_ordinal
            && claim_ordinal != admission_ordinal
        {
            return Err(SessionError::TerminalCorrupt {
                kind: KIND_SUBMISSION,
                authority: request_id.content_hash,
                detail: format!(
                    "submission claim admission ordinal {claim_ordinal} disagrees with terminal ordinal {admission_ordinal}"
                ),
            });
        }
        self.validate_recovered_terminal(
            &terminal,
            ledger_instance_id,
            KIND_SUBMISSION,
            request_id.content_hash,
            request_id.session,
            &ledger_scope,
            request_id.session_open,
            request_id.generation,
            terminal.claim.causal_ordinal,
            &payload,
            &expected_rows,
        )?;

        let mut inner = self.inner.lock().expect("governor lock");
        match inner.idempotency.get(&request_id) {
            Some(IdemState::Done {
                admission_ordinal,
                receipt,
                meter_receipt,
                ..
            }) => {
                return Ok(SubmitOutcome::Duplicate {
                    admission_ordinal: *admission_ordinal,
                    receipt: *receipt,
                    enforcement: meter_receipt.enforcement.clone(),
                    meter_receipt: meter_receipt.clone(),
                });
            }
            Some(IdemState::Failed {
                admission_ordinal,
                receipt,
                evidence,
                ..
            }) => {
                return Ok(SubmitOutcome::Failed {
                    admission_ordinal: *admission_ordinal,
                    receipt: *receipt,
                    evidence: evidence.as_ref().clone(),
                });
            }
            Some(IdemState::Pending { .. }) => return Ok(SubmitOutcome::InFlight),
            None => {}
        }
        let token = inner.tokens.get(&request_id.session.0).cloned().ok_or(
            SessionError::UnknownSession {
                id: request_id.session.0,
            },
        )?;
        if token.ledger_scope != ledger_scope
            || Self::current_open_identity(&inner, request_id.session)? != request_id.session_open
        {
            return Err(SessionError::MutationAuthorityMismatch {
                kind: KIND_SUBMISSION,
                id: request_id.content_hash,
            });
        }
        let current_generation = inner
            .gate_generations
            .get(&request_id.session.0)
            .copied()
            .unwrap_or(0);
        if request_id.generation > current_generation {
            return Err(SessionError::StaleMutationGeneration {
                kind: KIND_SUBMISSION,
                id: request_id.session.0,
                supplied: request_id.generation,
                current: current_generation,
            });
        }
        let key_index = inner
            .idempotency_keys
            .get(&request_id.session.0)
            .expect("recovered open created idempotency index");
        if let Some(existing) = key_index.get(&request_id.key_hash)
            && existing != &request_id
        {
            return Err(SessionError::MutationConflict {
                kind: KIND_SUBMISSION,
                id: request_id.content_hash,
            });
        }
        if key_index.len() >= MAX_IDEMPOTENCY_KEYS_PER_SESSION {
            return Err(SessionError::LimitExceeded {
                resource: "idempotency_keys_per_session",
                limit: MAX_IDEMPOTENCY_KEYS_PER_SESSION,
                observed_at_least: key_index.len().saturating_add(1),
            });
        }
        let (admission_ordinal, terminal_bytes, meter_receipt) = match &state {
            IdemState::Done {
                admission_ordinal,
                meter_receipt,
                ..
            } => (
                *admission_ordinal,
                enforcement_retained_bytes(&meter_receipt.enforcement),
                Some(meter_receipt.clone()),
            ),
            IdemState::Failed {
                admission_ordinal,
                evidence,
                ..
            } => (*admission_ordinal, evidence.preview.len(), None),
            IdemState::Pending { .. } => unreachable!("decoded terminal is never pending"),
        };
        if let Some(existing) = inner.submission_admissions.get(&admission_ordinal)
            && existing != &request_id
        {
            return Err(SessionError::TerminalCorrupt {
                kind: KIND_SUBMISSION,
                authority: request_id.content_hash,
                detail: format!(
                    "admission ordinal {admission_ordinal} is already owned by submission {}",
                    existing.content_hash
                ),
            });
        }
        let retained_bytes = SUBMISSION_REQUEST_RETAINED_BYTES
            .checked_add(terminal_bytes)
            .and_then(|bytes| {
                bytes.checked_add(if meter_receipt.is_some() {
                    MAX_METER_RECEIPT_RETAINED_BYTES
                } else {
                    0
                })
            })
            .ok_or(SessionError::LimitExceeded {
                resource: "retained_bytes_per_scope",
                limit: MAX_RETAINED_BYTES_PER_SCOPE,
                observed_at_least: usize::MAX,
            })?;
        ensure_retained_capacity(&inner, &ledger_scope, retained_bytes)?;
        if let Some(meter_receipt) = &meter_receipt {
            let report_id = meter_receipt.report_id;
            let report_count = inner
                .meter_report_ids
                .get(&request_id.session.0)
                .map_or(0, BTreeSet::len);
            if report_count >= MAX_METER_REPORTS_PER_SESSION {
                return Err(SessionError::LimitExceeded {
                    resource: "meter_reports_per_session",
                    limit: MAX_METER_REPORTS_PER_SESSION,
                    observed_at_least: report_count.saturating_add(1),
                });
            }
            let expected_meter_ordinal = inner.next_meter_commit_ordinal.checked_add(1).ok_or(
                SessionError::LimitExceeded {
                    resource: "meter_commit_ordinal",
                    limit: i64::MAX as usize,
                    observed_at_least: usize::MAX,
                },
            )?;
            if meter_receipt.commit_ordinal != expected_meter_ordinal {
                return Err(SessionError::RecoveryCausalGap {
                    session: request_id.session.0,
                    expected: expected_meter_ordinal,
                    found: meter_receipt.commit_ordinal,
                });
            }
            let before = inner
                .meters
                .get(&request_id.session.0)
                .cloned()
                .unwrap_or_default();
            if !same_snapshot(before.snapshot(), meter_receipt.before) {
                return Err(SessionError::TerminalCorrupt {
                    kind: KIND_SUBMISSION,
                    authority: request_id.content_hash,
                    detail: "submission meter before-state is not contiguous".to_string(),
                });
            }
            let (after, enforcement) = meter_transition(&token, &before, meter_receipt.delta)?;
            if !same_snapshot(after.snapshot(), meter_receipt.after)
                || enforcement != meter_receipt.enforcement
            {
                return Err(SessionError::TerminalCorrupt {
                    kind: KIND_SUBMISSION,
                    authority: request_id.content_hash,
                    detail: "submission meter after-state or enforcement does not recompute"
                        .to_string(),
                });
            }
            inner.next_meter_commit_ordinal = meter_receipt.commit_ordinal;
            inner.meters.insert(request_id.session.0, after);
            inner
                .meter_report_ids
                .get_mut(&request_id.session.0)
                .expect("recovered open created meter index")
                .insert(report_id);
            inner.meter_reports.insert(report_id, meter_receipt.clone());
        }
        inner.next_submission_ordinal = inner.next_submission_ordinal.max(admission_ordinal);
        inner
            .submission_admissions
            .insert(admission_ordinal, request_id);
        inner
            .idempotency_keys
            .get_mut(&request_id.session.0)
            .expect("recovered open created idempotency index")
            .insert(request_id.key_hash, request_id);
        inner.idempotency.insert(request_id, state);
        inner
            .scopes
            .get_mut(&ledger_scope)
            .expect("recovered scope")
            .sink
            .get_or_insert(ledger_instance_id);
        commit_retained_bytes(&mut inner, &ledger_scope, retained_bytes);
        self.mark_durable_claim_recovered(&mut inner, request_id.content_hash);

        match recovered {
            RecoveredSubmission::Done {
                admission_ordinal,
                receipt,
                meter_receipt,
                ..
            } => Ok(SubmitOutcome::Duplicate {
                admission_ordinal,
                receipt,
                enforcement: meter_receipt.enforcement.clone(),
                meter_receipt,
            }),
            RecoveredSubmission::Failed {
                admission_ordinal,
                receipt,
                evidence,
            } => Ok(SubmitOutcome::Failed {
                admission_ordinal,
                receipt,
                evidence,
            }),
        }
    }

    /// Reconstruct one durable declared pressure action and its indivisible
    /// initial event group.
    ///
    /// # Errors
    /// Missing/corrupt/conflicting terminals, causal event gaps, foreign or
    /// stale authority, gate mismatch, and capacity exhaustion fail closed.
    #[allow(clippy::too_many_lines)] // Full ladder/gate/reservation state is one recovery transition.
    pub fn recover_pressure(
        &self,
        ledger: &fs_ledger::Ledger,
        action_id: PressureActionId,
        level: u8,
    ) -> Result<PressureReceipt, SessionError> {
        if !(1..=3).contains(&level) {
            return Err(SessionError::InvalidPressureLevel { level });
        }
        if action_id.governor_id != self.id {
            return Err(SessionError::MutationAuthorityMismatch {
                kind: KIND_PRESSURE,
                id: action_id.content_hash,
            });
        }
        let ledger_instance_id = self.recovery_ledger(ledger)?;
        {
            let inner = self.inner.lock().expect("governor lock");
            if let Some(replay) = inner.pressure_actions.get(&action_id) {
                if replay.level != level {
                    return Err(SessionError::MutationConflict {
                        kind: KIND_PRESSURE,
                        id: action_id.content_hash,
                    });
                }
                let token =
                    inner
                        .tokens
                        .get(&action_id.session.0)
                        .ok_or(SessionError::UnknownSession {
                            id: action_id.session.0,
                        })?;
                let scope = inner
                    .scopes
                    .get(&token.ledger_scope)
                    .expect("registered scope");
                return Ok(PressureReceipt {
                    action_id,
                    level,
                    events: scope.events[replay.event_start..replay.event_start + replay.event_len]
                        .iter()
                        .map(|event| event.as_ref().clone())
                        .collect(),
                    content_hash: replay.content_hash,
                });
            }
        }
        let terminal = ledger
            .session_terminal(&action_id.content_hash)
            .map_err(|error| ledger_error("recover pressure terminal", error))?
            .ok_or(SessionError::RecoveryRequired {
                kind: KIND_PRESSURE,
                authority: action_id.content_hash,
            })?;
        let payload = encode_pressure_payload(level);
        if decode_pressure_payload(&terminal.claim.payload, action_id.content_hash)? != level {
            return Err(SessionError::MutationConflict {
                kind: KIND_PRESSURE,
                id: action_id.content_hash,
            });
        }
        let receipt = decode_pressure_terminal_receipt(&terminal.receipt, action_id)?;
        if receipt.level != level {
            return Err(SessionError::MutationConflict {
                kind: KIND_PRESSURE,
                id: action_id.content_hash,
            });
        }
        let ledger_scope = terminal.claim.ledger_scope.clone();
        let expected_events: Vec<_> = receipt
            .events
            .iter()
            .map(|event| buffered_degradation_event(&ledger_scope, event, receipt.content_hash))
            .collect::<Result<_, _>>()?;
        let expected_rows: Vec<_> = expected_events
            .iter()
            .map(BufferedLedgerEvent::as_row)
            .collect();
        let first_ordinal = receipt
            .events
            .first()
            .expect("decoded pressure receipt is nonempty")
            .ordinal;
        let causal_ordinal =
            u64::try_from(first_ordinal).map_err(|_| SessionError::TerminalCorrupt {
                kind: KIND_PRESSURE,
                authority: action_id.content_hash,
                detail: "pressure action first ordinal is negative".to_string(),
            })?;
        self.validate_recovered_terminal(
            &terminal,
            ledger_instance_id,
            KIND_PRESSURE,
            action_id.content_hash,
            action_id.session,
            &ledger_scope,
            action_id.session_open,
            action_id.generation,
            Some(causal_ordinal),
            &payload,
            &expected_rows,
        )?;

        let mut inner = self.inner.lock().expect("governor lock");
        if let Some(replay) = inner.pressure_actions.get(&action_id) {
            if replay.level != level || replay.content_hash != receipt.content_hash {
                return Err(SessionError::MutationConflict {
                    kind: KIND_PRESSURE,
                    id: action_id.content_hash,
                });
            }
            let token =
                inner
                    .tokens
                    .get(&action_id.session.0)
                    .ok_or(SessionError::UnknownSession {
                        id: action_id.session.0,
                    })?;
            let scope = inner
                .scopes
                .get(&token.ledger_scope)
                .expect("registered scope");
            return Ok(PressureReceipt {
                action_id,
                level,
                events: scope.events[replay.event_start..replay.event_start + replay.event_len]
                    .iter()
                    .map(|event| event.as_ref().clone())
                    .collect(),
                content_hash: replay.content_hash,
            });
        }
        let token = inner.tokens.get(&action_id.session.0).cloned().ok_or(
            SessionError::UnknownSession {
                id: action_id.session.0,
            },
        )?;
        if token.ledger_scope != ledger_scope
            || Self::current_open_identity(&inner, action_id.session)? != action_id.session_open
        {
            return Err(SessionError::MutationAuthorityMismatch {
                kind: KIND_PRESSURE,
                id: action_id.content_hash,
            });
        }
        let current_generation = inner
            .gate_generations
            .get(&action_id.session.0)
            .copied()
            .unwrap_or(0);
        if current_generation != action_id.generation {
            return Err(SessionError::StaleMutationGeneration {
                kind: KIND_PRESSURE,
                id: action_id.session.0,
                supplied: action_id.generation,
                current: current_generation,
            });
        }
        let action_count = inner
            .pressure_action_ids
            .get(&action_id.session.0)
            .map_or(0, BTreeSet::len);
        if action_count >= MAX_PRESSURE_ACTIONS_PER_SESSION {
            return Err(SessionError::LimitExceeded {
                resource: "pressure_actions_per_session",
                limit: MAX_PRESSURE_ACTIONS_PER_SESSION,
                observed_at_least: action_count.saturating_add(1),
            });
        }
        let expected_first =
            inner
                .next_ordinal
                .checked_add(1)
                .ok_or(SessionError::LimitExceeded {
                    resource: "degradation_ordinal",
                    limit: i64::MAX as usize,
                    observed_at_least: usize::MAX,
                })?;
        if first_ordinal != expected_first
            || receipt.events.iter().enumerate().any(|(index, event)| {
                first_ordinal.checked_add(i64::try_from(index).expect("bounded index fits i64"))
                    != Some(event.ordinal)
            })
        {
            return Err(SessionError::TerminalCorrupt {
                kind: KIND_PRESSURE,
                authority: action_id.content_hash,
                detail: "pressure event ordinals do not continue the recovered dense prefix"
                    .to_string(),
            });
        }
        let immediate_bytes = receipt.events.iter().try_fold(0usize, |total, event| {
            total
                .checked_add(degradation_event_retained_bytes(event)?)
                .ok_or(SessionError::LimitExceeded {
                    resource: "retained_bytes_per_scope",
                    limit: MAX_RETAINED_BYTES_PER_SCOPE,
                    observed_at_least: usize::MAX,
                })
        })?;
        let is_pause = level == 3;
        let scope = inner.scopes.get(&ledger_scope).expect("recovered scope");
        let required_event_capacity = scope
            .events
            .len()
            .checked_add(scope.reserved_pause_completions)
            .and_then(|count| count.checked_add(receipt.events.len()))
            .and_then(|count| count.checked_add(usize::from(is_pause)))
            .ok_or(SessionError::LimitExceeded {
                resource: "degradation_events_per_scope",
                limit: MAX_DEGRADATION_EVENTS_PER_SCOPE,
                observed_at_least: usize::MAX,
            })?;
        if required_event_capacity > MAX_DEGRADATION_EVENTS_PER_SCOPE {
            return Err(SessionError::LimitExceeded {
                resource: "degradation_events_per_scope",
                limit: MAX_DEGRADATION_EVENTS_PER_SCOPE,
                observed_at_least: required_event_capacity,
            });
        }
        let reserved_completion_bytes = if is_pause {
            MAX_PAUSE_COMPLETION_RETAINED_BYTES
        } else {
            0
        };
        let retained_bytes = immediate_bytes
            .checked_add(reserved_completion_bytes)
            .and_then(|bytes| bytes.checked_add(PRESSURE_ACTION_RETAINED_BYTES))
            .ok_or(SessionError::LimitExceeded {
                resource: "retained_bytes_per_scope",
                limit: MAX_RETAINED_BYTES_PER_SCOPE,
                observed_at_least: usize::MAX,
            })?;
        ensure_retained_capacity(&inner, &ledger_scope, retained_bytes)?;
        if is_pause {
            inner
                .next_ordinal
                .checked_add(i64::try_from(receipt.events.len()).expect("bounded events fit i64"))
                .and_then(|ordinal| ordinal.checked_add(1))
                .ok_or(SessionError::LimitExceeded {
                    resource: "degradation_ordinal",
                    limit: i64::MAX as usize,
                    observed_at_least: usize::MAX,
                })?;
            if inner.pending_pause.contains_key(&action_id.session.0) {
                return Err(SessionError::PauseAlreadyPending {
                    id: action_id.session.0,
                    requested_ordinal: inner.pending_pause[&action_id.session.0]
                        .request_id
                        .requested_ordinal,
                });
            }
            let request = receipt
                .events
                .last()
                .and_then(|event| event.pause_request_id)
                .ok_or_else(|| SessionError::TerminalCorrupt {
                    kind: KIND_PRESSURE,
                    authority: action_id.content_hash,
                    detail: "level-3 pressure terminal lacks pause request authority".to_string(),
                })?;
            let gate = inner.gates.get(&action_id.session.0).cloned().ok_or(
                SessionError::UngatedSession {
                    id: action_id.session.0,
                },
            )?;
            if inner.gate_phases.get(&action_id.session.0) != Some(&GatePhase::Running) {
                return Err(SessionError::SessionGateDraining {
                    id: action_id.session.0,
                    generation: action_id.generation,
                });
            }
            gate.request();
            inner.pending_pause.insert(
                action_id.session.0,
                PendingPause {
                    request_id: request,
                    pressure_action_id: action_id,
                    reserved_retained_bytes: reserved_completion_bytes,
                },
            );
            inner.reserved_pause_ordinals = inner.reserved_pause_ordinals.checked_add(1).ok_or(
                SessionError::LimitExceeded {
                    resource: "pause_ordinal_reservations",
                    limit: usize::MAX,
                    observed_at_least: usize::MAX,
                },
            )?;
            inner
                .scopes
                .get_mut(&ledger_scope)
                .expect("recovered scope")
                .reserved_pause_completions += 1;
        }
        let event_start = inner
            .scopes
            .get(&ledger_scope)
            .expect("recovered scope")
            .events
            .len();
        inner
            .scopes
            .get_mut(&ledger_scope)
            .expect("recovered scope")
            .events
            .extend(receipt.events.iter().cloned().map(Arc::new));
        inner.next_ordinal = receipt
            .events
            .last()
            .expect("decoded pressure receipt is nonempty")
            .ordinal;
        inner.pressure_actions.insert(
            action_id,
            PressureReplay {
                level,
                event_start,
                event_len: receipt.events.len(),
                content_hash: receipt.content_hash,
            },
        );
        inner
            .pressure_action_ids
            .get_mut(&action_id.session.0)
            .expect("recovered open created pressure index")
            .insert(action_id);
        inner
            .scopes
            .get_mut(&ledger_scope)
            .expect("recovered scope")
            .sink
            .get_or_insert(ledger_instance_id);
        commit_retained_bytes(&mut inner, &ledger_scope, retained_bytes);
        self.mark_durable_claim_recovered(&mut inner, action_id.content_hash);
        Ok(receipt)
    }

    /// Reconstruct one durable pause acknowledgement onto a caller-supplied
    /// fresh process-local gate.
    ///
    /// # Errors
    /// Missing/corrupt/conflicting terminals, noncontiguous lifecycle state,
    /// pending work, stale authority, or a requested replacement gate fail
    /// closed.
    #[allow(clippy::too_many_lines)] // Checkpoint, event, reservation, and gate rotation are one transition.
    pub fn recover_pause_acknowledgement(
        &self,
        ledger: &fs_ledger::Ledger,
        request_id: PauseRequestId,
        checkpoint: &fs_ledger::session_registry::SolverCheckpointReceipt,
        resume_gate: Arc<CancelGate>,
    ) -> Result<PauseAcknowledgement, SessionError> {
        if request_id.governor_id != self.id {
            return Err(SessionError::PauseRequestMismatch {
                id: request_id.session.0,
                requested_ordinal: request_id.requested_ordinal,
            });
        }
        let authority = pause_ack_authority(request_id);
        let checkpoint_refusal = |reason| SessionError::PauseCheckpointMismatch {
            id: request_id.session.0,
            requested_ordinal: request_id.requested_ordinal,
            reason,
        };
        if checkpoint.session() != request_id.session.0 {
            return Err(checkpoint_refusal("cross-session"));
        }
        if checkpoint.gate_generation() != request_id.gate_generation {
            return Err(checkpoint_refusal("stale-generation"));
        }
        if checkpoint.pause_authority() != authority {
            return Err(checkpoint_refusal("foreign-pause-request"));
        }
        ledger
            .verify_solver_checkpoint_receipt(checkpoint)
            .map_err(|_| checkpoint_refusal("unverified-ledger-receipt"))?;
        let ledger_instance_id = self.recovery_ledger(ledger)?;
        let terminal = ledger
            .session_terminal(&authority)
            .map_err(|error| ledger_error("recover pause acknowledgement terminal", error))?
            .ok_or(SessionError::RecoveryRequired {
                kind: KIND_PAUSE_ACK,
                authority,
            })?;
        let evidence = RetainedEvidence::capture(&checkpoint.content_hash().to_hex());
        let payload = encode_pause_ack_payload(&evidence);
        if decode_pause_ack_payload(&terminal.claim.payload, request_id)? != evidence {
            return Err(SessionError::PauseAcknowledgementConflict {
                id: request_id.session.0,
                requested_ordinal: request_id.requested_ordinal,
            });
        }
        let acknowledgement = decode_pause_ack_terminal_receipt(
            &terminal.receipt,
            request_id,
            Arc::clone(&resume_gate),
        )?;
        if acknowledgement.event.checkpoint.as_ref() != Some(&evidence) {
            return Err(SessionError::TerminalCorrupt {
                kind: KIND_PAUSE_ACK,
                authority,
                detail: "pause receipt checkpoint evidence differs from its payload".to_string(),
            });
        }
        let (ledger_scope, session_open, action_receipt_hash) = {
            let inner = self.inner.lock().expect("governor lock");
            if let Some(replayed) =
                cached_pause_acknowledgement(&inner, request_id, &acknowledgement, &resume_gate)?
            {
                return Ok(replayed);
            }
            if resume_gate.is_requested() {
                return Err(SessionError::ResumeGateAlreadyRequested {
                    id: request_id.session.0,
                    generation: request_id.gate_generation.saturating_add(1),
                });
            }
            let token =
                inner
                    .tokens
                    .get(&request_id.session.0)
                    .ok_or(SessionError::UnknownSession {
                        id: request_id.session.0,
                    })?;
            let session_open = Self::current_open_identity(&inner, request_id.session)?;
            let action_id = acknowledgement.event.pressure_action_id.ok_or_else(|| {
                SessionError::TerminalCorrupt {
                    kind: KIND_PAUSE_ACK,
                    authority,
                    detail: "pause completion lacks originating pressure action".to_string(),
                }
            })?;
            let action_receipt_hash = inner
                .pressure_actions
                .get(&action_id)
                .ok_or_else(|| SessionError::RecoveryRequired {
                    kind: KIND_PRESSURE,
                    authority: action_id.content_hash,
                })?
                .content_hash;
            (
                token.ledger_scope.clone(),
                session_open,
                action_receipt_hash,
            )
        };
        if let Some(pending) = ledger
            .pending_session_mutation(
                self.id,
                session_open,
                KIND_SUBMISSION,
                request_id.session.0,
                &ledger_scope,
                request_id.gate_generation,
            )
            .map_err(|error| ledger_error("recover pause Pending-claim probe", error))?
        {
            return Err(SessionError::IndeterminateMutation {
                kind: KIND_SUBMISSION,
                authority: pending.authority,
            });
        }
        let expected_event =
            buffered_degradation_event(&ledger_scope, &acknowledgement.event, action_receipt_hash)?;
        let expected_rows = [expected_event.as_row()];
        let causal_ordinal = u64::try_from(acknowledgement.event.ordinal).map_err(|_| {
            SessionError::TerminalCorrupt {
                kind: KIND_PAUSE_ACK,
                authority,
                detail: "pause completion ordinal is negative".to_string(),
            }
        })?;
        self.validate_recovered_terminal(
            &terminal,
            ledger_instance_id,
            KIND_PAUSE_ACK,
            authority,
            request_id.session,
            &ledger_scope,
            session_open,
            acknowledgement.resume_generation,
            Some(causal_ordinal),
            &payload,
            &expected_rows,
        )?;

        let mut inner = self.inner.lock().expect("governor lock");
        if let Some(replayed) =
            cached_pause_acknowledgement(&inner, request_id, &acknowledgement, &resume_gate)?
        {
            return Ok(replayed);
        }
        if resume_gate.is_requested() {
            return Err(SessionError::ResumeGateAlreadyRequested {
                id: request_id.session.0,
                generation: request_id.gate_generation.saturating_add(1),
            });
        }
        if inner
            .pending_submissions
            .get(&request_id.session.0)
            .copied()
            .unwrap_or(0)
            != 0
        {
            return Err(SessionError::PauseDrainPending {
                id: request_id.session.0,
                pending_submissions: inner.pending_submissions[&request_id.session.0],
            });
        }
        let pending = inner
            .pending_pause
            .get(&request_id.session.0)
            .copied()
            .filter(|pending| pending.request_id == request_id)
            .ok_or(SessionError::PauseRequestMismatch {
                id: request_id.session.0,
                requested_ordinal: request_id.requested_ordinal,
            })?;
        if acknowledgement.event.pressure_action_id != Some(pending.pressure_action_id) {
            return Err(SessionError::TerminalCorrupt {
                kind: KIND_PAUSE_ACK,
                authority,
                detail: "pause completion names a different originating pressure action"
                    .to_string(),
            });
        }
        let old_gate =
            inner
                .gates
                .get(&request_id.session.0)
                .ok_or(SessionError::UngatedSession {
                    id: request_id.session.0,
                })?;
        if inner.gate_generations.get(&request_id.session.0) != Some(&request_id.gate_generation)
            || inner.gate_phases.get(&request_id.session.0) != Some(&GatePhase::Running)
            || !old_gate.is_requested()
        {
            return Err(SessionError::PauseRequestMismatch {
                id: request_id.session.0,
                requested_ordinal: request_id.requested_ordinal,
            });
        }
        if inner.next_ordinal.checked_add(1) != Some(acknowledgement.event.ordinal)
            || inner.reserved_pause_ordinals == 0
            || inner
                .scopes
                .get(&ledger_scope)
                .expect("recovered scope")
                .reserved_pause_completions
                == 0
        {
            return Err(SessionError::TerminalCorrupt {
                kind: KIND_PAUSE_ACK,
                authority,
                detail: "pause completion does not consume the reserved next event/ordinal"
                    .to_string(),
            });
        }
        let event_bytes = degradation_event_retained_bytes(&acknowledgement.event)?;
        if event_bytes > pending.reserved_retained_bytes {
            return Err(SessionError::TerminalCorrupt {
                kind: KIND_PAUSE_ACK,
                authority,
                detail: "pause completion exceeds its recovered retained-byte reservation"
                    .to_string(),
            });
        }
        inner.pending_pause.remove(&request_id.session.0);
        inner.reserved_pause_ordinals -= 1;
        inner
            .gates
            .insert(request_id.session.0, Arc::clone(&resume_gate));
        inner
            .gate_generations
            .insert(request_id.session.0, acknowledgement.resume_generation);
        inner
            .gate_phases
            .insert(request_id.session.0, GatePhase::ReadyToResume);
        let event_index = inner
            .scopes
            .get(&ledger_scope)
            .expect("recovered scope")
            .events
            .len();
        inner.completed_pause.insert(
            request_id.session.0,
            CompletedPause {
                request_id,
                checkpoint_byte_len: evidence.byte_len,
                checkpoint_digest: evidence.digest,
                completion_event_index: event_index,
                completion_ordinal: acknowledgement.event.ordinal,
                resume_generation: acknowledgement.resume_generation,
                gate_binding: acknowledgement.gate_binding,
                acknowledgement_hash: acknowledgement.content_hash,
            },
        );
        inner.pause_acknowledgements.insert(
            request_id,
            PauseAcknowledgementReplay {
                completion_event_index: event_index,
                resume_generation: acknowledgement.resume_generation,
                gate_binding: acknowledgement.gate_binding,
                content_hash: acknowledgement.content_hash,
            },
        );
        inner.next_ordinal = acknowledgement.event.ordinal;
        {
            let scope = inner
                .scopes
                .get_mut(&ledger_scope)
                .expect("recovered scope");
            scope.reserved_pause_completions -= 1;
            scope.events.push(Arc::new(acknowledgement.event.clone()));
            scope.sink.get_or_insert(ledger_instance_id);
        }
        release_retained_bytes(
            &mut inner,
            &ledger_scope,
            pending.reserved_retained_bytes - event_bytes,
        );
        self.mark_durable_claim_recovered(&mut inner, authority);
        Ok(acknowledgement)
    }

    /// Reconstruct the zero-event durable terminal that proves resumed workers
    /// adopted an acknowledged gate generation.
    ///
    /// # Errors
    /// Missing/corrupt/conflicting terminals, stale acknowledgement/gate, or
    /// foreign recovery identity fail closed.
    pub fn recover_resume_activation(
        &self,
        ledger: &fs_ledger::Ledger,
        acknowledgement: &PauseAcknowledgement,
    ) -> Result<ResumeActivationReceipt, SessionError> {
        let request_id = acknowledgement.request_id;
        if request_id.governor_id != self.id {
            return Err(SessionError::ResumeAcknowledgementMismatch {
                id: request_id.session.0,
            });
        }
        let ledger_instance_id = self.recovery_ledger(ledger)?;
        let (session_open, ledger_scope) =
            {
                let inner = self.inner.lock().expect("governor lock");
                let token = inner.tokens.get(&request_id.session.0).ok_or(
                    SessionError::UnknownSession {
                        id: request_id.session.0,
                    },
                )?;
                (
                    Self::current_open_identity(&inner, request_id.session)?,
                    token.ledger_scope.clone(),
                )
            };
        let activation_id = resume_activation_id(
            self.id,
            request_id.session,
            session_open,
            acknowledgement.content_hash,
            acknowledgement.resume_generation,
        );
        {
            let inner = self.inner.lock().expect("governor lock");
            if let Some(replayed) =
                cached_resume_activation(&inner, activation_id, acknowledgement)?
            {
                return Ok(replayed);
            }
        }
        let terminal = ledger
            .session_terminal(&activation_id.content_hash)
            .map_err(|error| ledger_error("recover resume activation terminal", error))?
            .ok_or(SessionError::RecoveryRequired {
                kind: KIND_RESUME_ACTIVATION,
                authority: activation_id.content_hash,
            })?;
        let payload = encode_activation_payload(acknowledgement);
        let decoded_payload = decode_activation_payload(&terminal.claim.payload, activation_id)?;
        if decoded_payload
            != (
                acknowledgement.content_hash,
                acknowledgement.resume_generation,
                acknowledgement.gate_binding,
            )
        {
            return Err(SessionError::MutationConflict {
                kind: KIND_RESUME_ACTIVATION,
                id: activation_id.content_hash,
            });
        }
        let receipt = decode_activation_terminal_receipt(&terminal.receipt, activation_id)?;
        if receipt.acknowledgement_hash != acknowledgement.content_hash
            || receipt.gate_binding != acknowledgement.gate_binding
        {
            return Err(SessionError::TerminalCorrupt {
                kind: KIND_RESUME_ACTIVATION,
                authority: activation_id.content_hash,
                detail: "activation receipt differs from acknowledgement".to_string(),
            });
        }
        self.validate_recovered_terminal(
            &terminal,
            ledger_instance_id,
            KIND_RESUME_ACTIVATION,
            activation_id.content_hash,
            request_id.session,
            &ledger_scope,
            session_open,
            acknowledgement.resume_generation,
            Some(u64::try_from(acknowledgement.event.ordinal).map_err(|_| {
                SessionError::TerminalCorrupt {
                    kind: KIND_RESUME_ACTIVATION,
                    authority: activation_id.content_hash,
                    detail: "activation causal ordinal is negative".to_string(),
                }
            })?),
            &payload,
            &[],
        )?;
        let mut inner = self.inner.lock().expect("governor lock");
        if let Some(replayed) = cached_resume_activation(&inner, activation_id, acknowledgement)? {
            return Ok(replayed);
        }
        let completed = inner
            .completed_pause
            .get(&request_id.session.0)
            .copied()
            .ok_or(SessionError::ResumeAcknowledgementMismatch {
                id: request_id.session.0,
            })?;
        let current_gate = inner.gates.get(&request_id.session.0).cloned().ok_or(
            SessionError::UngatedSession {
                id: request_id.session.0,
            },
        )?;
        if completed.request_id != request_id
            || completed.acknowledgement_hash != acknowledgement.content_hash
            || completed.gate_binding != acknowledgement.gate_binding
            || inner.gate_phases.get(&request_id.session.0) != Some(&GatePhase::ReadyToResume)
            || !Arc::ptr_eq(&current_gate, &acknowledgement.resume_gate)
            || current_gate.is_requested()
        {
            return Err(SessionError::ResumeAcknowledgementMismatch {
                id: request_id.session.0,
            });
        }
        inner
            .gate_phases
            .insert(request_id.session.0, GatePhase::Running);
        inner.resume_activations.insert(activation_id, receipt);
        inner
            .scopes
            .get_mut(&ledger_scope)
            .expect("recovered scope")
            .sink
            .get_or_insert(ledger_instance_id);
        self.mark_durable_claim_recovered(&mut inner, activation_id.content_hash);
        Ok(receipt)
    }
}
