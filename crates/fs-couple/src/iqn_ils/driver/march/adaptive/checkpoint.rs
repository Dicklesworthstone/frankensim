//! Portable bytes for an adaptive attempt-boundary checkpoint.
//!
//! The caller binds the physical model, mesh/order, initial data, forcing,
//! estimator, state codec and executable semantics through `model_binding`.
//! The runtime separately binds every numerical control and mandatory endpoint.
//! Payload encoding/admission remains with the state owner: hashing opaque bytes
//! does not prove that they describe a valid physical state. The digest detects
//! corruption, NOT authenticity; a malicious writer can recompute it.
//!
//! Restore performs no physics and commits only after every check and the state
//! decoder succeeds. File durability/atomic publication belong to the consumer.
//! Bitwise continuation requires the same admitted producer/toolchain profile;
//! portable little-endian encoding alone is not a cross-ISA replay guarantee.

use fs_blake3::{DomainHasher, hash_domain};
use super::AdaptiveEvolution;
use super::super::CouplingMethod;

const MAGIC: &[u8; 8] = b"FSCPAD01";
const DOMAIN: &str = "org.frankensim.couple.adaptive-checkpoint.v1";
const MODEL_DOMAIN: &str = "org.frankensim.couple.adaptive-checkpoint.model.v1";
const POLICY_DOMAIN: &str = "org.frankensim.couple.adaptive-checkpoint.policy.v1";
/// Hard ceiling; callers can and should declare a smaller per-checkpoint cap.
pub const MAX_CHECKPOINT_BYTES: usize = 256 * 1024 * 1024;
/// Bound on the caller's complete, canonical model description.
pub const MAX_MODEL_BINDING_BYTES: usize = 64 * 1024;

/// A refused checkpoint; the destination evolution remains unchanged.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CheckpointError {
    /// Byte count, allocation or platform integer range exceeded its limit.
    Limit,
    /// Unsupported version, truncated/trailing data or inconsistent cursor.
    Format,
    /// The content digest does not match the retained bytes.
    Integrity,
    /// Empty/oversized model binding, or a different physical/codec identity.
    Model,
    /// Coupling controls, time policy or mandatory schedule differ.
    Policy,
    /// The state owner refused its bounded payload.
    State(String),
}
impl core::fmt::Display for CheckpointError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "adaptive checkpoint refused: {self:?}")
    }
}
impl std::error::Error for CheckpointError {}

fn binding(bytes: &[u8]) -> Result<[u8; 32], CheckpointError> {
    if bytes.is_empty() || bytes.len() > MAX_MODEL_BINDING_BYTES {
        return Err(CheckpointError::Model);
    }
    Ok(*hash_domain(MODEL_DOMAIN, bytes).as_bytes())
}
fn integer(out: &mut Vec<u8>, value: usize) -> Result<(), CheckpointError> {
    out.extend_from_slice(&u64::try_from(value).map_err(|_| CheckpointError::Limit)?.to_le_bytes());
    Ok(())
}
fn real(out: &mut Vec<u8>, value: f64) { out.extend_from_slice(&value.to_bits().to_le_bytes()); }

impl<S: Clone> AdaptiveEvolution<S> {
    fn checkpoint_policy(&self) -> [u8; 32] {
        let mut h = DomainHasher::new(POLICY_DOMAIN);
        let c = &self.base.controls;
        h.update(&(c.max_evaluations as u64).to_le_bytes());
        h.update(&c.relaxation.to_bits().to_le_bytes());
        match c.method {
            CouplingMethod::RelaxedPicard => h.update(&[0]),
            CouplingMethod::IqnIls(config) => {
                h.update(&[1]);
                h.update(&(config.max_history as u64).to_le_bytes());
                h.update(&config.relative_rank_tolerance.to_bits().to_le_bytes());
            }
        }
        h.update(&(c.interfaces.len() as u64).to_le_bytes());
        for rule in &c.interfaces {
            for value in [rule.scale, rule.absolute_tolerance, rule.relative_tolerance] {
                h.update(&value.to_bits().to_le_bytes());
            }
        }
        h.update(&(c.balances.len() as u64).to_le_bytes());
        for rule in &c.balances {
            h.update(&(rule.name.len() as u64).to_le_bytes());
            h.update(rule.name.as_bytes());
            h.update(&rule.absolute_tolerance.to_bits().to_le_bytes());
        }
        for value in [self.settings.initial_step_s, self.settings.minimum_step_s,
            self.settings.maximum_step_s] { h.update(&value.to_bits().to_le_bytes()); }
        h.update(&self.settings.method_order.to_le_bytes());
        h.update(&(self.base.times_s.len() as u64).to_le_bytes());
        for value in &self.base.times_s { h.update(&value.to_bits().to_le_bytes()); }
        *h.finalize().as_bytes()
    }

    /// Encode the complete runtime continuation and caller-encoded physical state.
    ///
    /// `state_payload` must encode exactly `self.state()` and any consumer-owned
    /// accounting needed on restart. `model_binding` must identify all inputs
    /// not owned by this runtime, including the distance callback and codec.
    /// The byte cap includes the header, interface, payload and digest. No
    /// allocation occurs until the full output length has been admitted.
    pub fn checkpoint_bytes(&self, model_binding: &[u8], state_payload: &[u8], max_bytes: usize)
        -> Result<Vec<u8>, CheckpointError>
    {
        let model = binding(model_binding)?;
        let length = self.base.interface.len().checked_mul(8)
            .and_then(|n| n.checked_add(152))
            .and_then(|n| n.checked_add(state_payload.len())).ok_or(CheckpointError::Limit)?;
        if length > max_bytes.min(MAX_CHECKPOINT_BYTES) { return Err(CheckpointError::Limit); }
        let mut out = Vec::new();
        out.try_reserve_exact(length).map_err(|_| CheckpointError::Limit)?;
        out.extend_from_slice(MAGIC);
        out.extend_from_slice(&model);
        out.extend_from_slice(&self.checkpoint_policy());
        integer(&mut out, self.base.next)?;
        real(&mut out, self.time_s);
        real(&mut out, self.next_step_s);
        integer(&mut out, self.accepted_steps)?;
        integer(&mut out, self.base.interface.len())?;
        for &value in &self.base.interface { real(&mut out, value); }
        integer(&mut out, state_payload.len())?;
        out.extend_from_slice(state_payload);
        let digest = hash_domain(DOMAIN, &out);
        out.extend_from_slice(digest.as_bytes());
        debug_assert_eq!(out.len(), length);
        Ok(out)
    }

    /// Restore into an independently constructed evolution with the same policy.
    ///
    /// Length, digest, model/policy binding, finite interface values and cursor
    /// consistency are checked BEFORE invoking `decode_state`. The decoder must
    /// admit physical field shape/values and consume its complete payload. The
    /// old state is replaced only after it succeeds. No producer is evaluated.
    /// Accepted time, the possibly shrunken next duration, mandatory-endpoint
    /// cursor and accepted count are restored together, including after rejection.
    pub fn restore_checkpoint<F>(&mut self, bytes: &[u8], model_binding: &[u8], max_bytes: usize,
        decode_state: F) -> Result<(), CheckpointError>
    where F: FnOnce(&[u8]) -> Result<S, String>
    {
        if bytes.len() > max_bytes.min(MAX_CHECKPOINT_BYTES) { return Err(CheckpointError::Limit); }
        if bytes.len() < 152 || &bytes[..8] != MAGIC { return Err(CheckpointError::Format); }
        let split = bytes.len() - 32;
        if hash_domain(DOMAIN, &bytes[..split]).as_bytes().as_slice() != &bytes[split..] {
            return Err(CheckpointError::Integrity);
        }
        let mut input = Reader(&bytes[8..split]);
        if input.take(32)? != binding(model_binding)? { return Err(CheckpointError::Model); }
        if input.take(32)? != self.checkpoint_policy() { return Err(CheckpointError::Policy); }
        let cursor = input.integer()?;
        let time = input.real()?;
        let next_step = input.real()?;
        let accepted = input.integer()?;
        let dimension = input.integer()?;
        if cursor >= self.base.times_s.len() || dimension != self.base.interface.len()
            || !time.is_finite() || !next_step.is_finite()
            || next_step < self.settings.minimum_step_s || next_step > self.settings.maximum_step_s
        { return Err(CheckpointError::Format); }
        let lower = self.base.times_s[cursor];
        let complete = cursor + 1 == self.base.times_s.len();
        if time < lower || (complete && time != lower)
            || (!complete && time >= self.base.times_s[cursor + 1])
            || accepted < cursor || (time > lower && accepted == cursor)
            || (time == self.base.times_s[0] && accepted != 0)
        { return Err(CheckpointError::Format); }
        let mut interface = Vec::new();
        // The expected dimension comes from the already admitted destination,
        // not from an unchecked length in the file. Check bytes before allocation.
        let raw = input.take(dimension.checked_mul(8).ok_or(CheckpointError::Limit)?)?;
        interface.try_reserve_exact(dimension).map_err(|_| CheckpointError::Limit)?;
        for chunk in raw.chunks_exact(8) {
            let value = f64::from_le_bytes(chunk.try_into().map_err(|_| CheckpointError::Format)?);
            if !value.is_finite() { return Err(CheckpointError::Format); }
            interface.push(value);
        }
        let payload_len = input.integer()?;
        let payload = input.take(payload_len)?;
        if !input.0.is_empty() { return Err(CheckpointError::Format); }
        let state = decode_state(payload).map_err(CheckpointError::State)?;
        self.base.state = state;
        self.base.interface = interface;
        self.base.next = cursor;
        self.time_s = time;
        self.next_step_s = next_step;
        self.accepted_steps = accepted;
        Ok(())
    }
}

struct Reader<'a>(&'a [u8]);
impl<'a> Reader<'a> {
    fn take(&mut self, count: usize) -> Result<&'a [u8], CheckpointError> {
        if count > self.0.len() { return Err(CheckpointError::Format); }
        let (value, rest) = self.0.split_at(count); self.0 = rest; Ok(value)
    }
    fn word(&mut self) -> Result<[u8; 8], CheckpointError> {
        self.take(8)?.try_into().map_err(|_| CheckpointError::Format)
    }
    fn integer(&mut self) -> Result<usize, CheckpointError> {
        usize::try_from(u64::from_le_bytes(self.word()?)).map_err(|_| CheckpointError::Limit)
    }
    fn real(&mut self) -> Result<f64, CheckpointError> { Ok(f64::from_le_bytes(self.word()?)) }
}

#[cfg(test)]
mod tests;
