//! State ring buffer + record/replay (bead wf-root-guzez.4.2.2, E3.2-ii).
//! First CODE consumer of the E0.9a-frozen input-trace contract
//! (data/wright-flyer/replay-identity-schema-v1.json): `InputTraceV1`
//! carries {schema id, end_tick_exclusive, ordered resolved applied
//! events}, the trace EXTENT is part of the hash domain (two event-free
//! runs stopped at different ticks differ), and
//! `InputTraceId = H("fs-flyer/applied-input-trace/v1", canonical bytes)`.
//!
//! Replay verification is LOCALIZING: a mismatch names the first divergent
//! tick, not just "hashes differ" — the seed of the E3.5 structured-
//! checkpoint program.

use crate::Refusal;
use crate::spine::{Loads, RigidBody, SixDofState, step, tick_digest};
use fs_blake3::hash_domain;

/// The frozen input-trace hash domain (replay-identity-schema-v1).
pub const INPUT_TRACE_DOMAIN: &str = "fs-flyer/applied-input-trace/v1";
/// The v1 input-trace schema id carried inside the preimage.
pub const INPUT_TRACE_SCHEMA_ID: &str = "org.frankensim.wright-flyer.input-trace.v1";
/// Ring-buffer capacity cap (slots).
pub const MAX_RING_CAPACITY: usize = 16_384;
/// Control-channel count cap.
pub const MAX_CHANNELS: usize = 16;

/// One resolved applied event (the CANONICAL trace element — no
/// acquisition-clock fields, per the frozen Round-5 boundary).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AppliedEvent {
    /// Control channel index.
    pub channel: u32,
    /// Simulation tick at which the value applies.
    pub applied_tick: u32,
    /// Ordinal within the tick (canonical determinism for same-tick events).
    pub ordinal_within_tick: u32,
    /// Quantized control value (1/4096 grid upstream).
    pub quantized_value: f64,
}

/// The closed applied-input trace (InputTraceV1 shape).
#[derive(Clone, Debug, PartialEq)]
pub struct InputTrace {
    /// Trace extent: one past the last simulated tick (part of the DOMAIN).
    pub end_tick_exclusive: u32,
    /// Ordered resolved applied events.
    pub events: Vec<AppliedEvent>,
}

impl InputTrace {
    /// Validate ordering (applied_tick ascending; ordinals ascending within
    /// a tick and dense from 0) and extent.
    ///
    /// # Errors
    /// `trace-order-invalid`, `trace-extent-invalid`.
    pub fn admit(&self) -> Result<(), Refusal> {
        let mut prev: Option<(u32, u32)> = None;
        for e in &self.events {
            if e.applied_tick >= self.end_tick_exclusive {
                return Err(Refusal {
                    code: "trace-extent-invalid",
                    message: format!(
                        "event at tick {} outside the trace extent {}",
                        e.applied_tick, self.end_tick_exclusive
                    ),
                    ranked_repairs: vec!["extend end_tick_exclusive or drop the event".into()],
                });
            }
            if !e.quantized_value.is_finite() {
                return Err(Refusal {
                    code: "trace-order-invalid",
                    message: format!("non-finite value on channel {}", e.channel),
                    ranked_repairs: vec!["quantized values are finite by construction".into()],
                });
            }
            if let Some((pt, po)) = prev {
                let ok = e.applied_tick > pt
                    || (e.applied_tick == pt && e.ordinal_within_tick == po + 1);
                if !ok {
                    return Err(Refusal {
                        code: "trace-order-invalid",
                        message: format!(
                            "event ({}, ord {}) after ({pt}, ord {po}) breaks canonical order",
                            e.applied_tick, e.ordinal_within_tick
                        ),
                        ranked_repairs: vec![
                            "events sort by (applied_tick, ordinal); ordinals are dense".into(),
                        ],
                    });
                }
            } else if e.ordinal_within_tick != 0 && self.events.first().map(|f| f.applied_tick)
                == Some(e.applied_tick)
            {
                return Err(Refusal {
                    code: "trace-order-invalid",
                    message: "first event of a tick must carry ordinal 0".into(),
                    ranked_repairs: vec!["re-run the input scheduler".into()],
                });
            }
            prev = Some((e.applied_tick, e.ordinal_within_tick));
        }
        Ok(())
    }

    /// The frozen `InputTraceId`: domain-separated hash of the canonical
    /// bytes {schema id, end_tick_exclusive, events in order}. The EXTENT
    /// is in the domain by construction.
    #[must_use]
    pub fn trace_id(&self) -> String {
        let mut payload = Vec::with_capacity(16 + 20 * self.events.len());
        payload.extend_from_slice(INPUT_TRACE_SCHEMA_ID.as_bytes());
        payload.extend_from_slice(&self.end_tick_exclusive.to_le_bytes());
        for e in &self.events {
            payload.extend_from_slice(&e.channel.to_le_bytes());
            payload.extend_from_slice(&e.applied_tick.to_le_bytes());
            payload.extend_from_slice(&e.ordinal_within_tick.to_le_bytes());
            payload.extend_from_slice(&e.quantized_value.to_bits().to_le_bytes());
        }
        hash_domain(INPUT_TRACE_DOMAIN, &payload).to_hex()
    }
}

/// Fixed-capacity tick-indexed state ring (the presentation/checkpoint
/// window). Slots are overwritten oldest-first; reads outside the live
/// window refuse.
#[derive(Clone, Debug)]
pub struct StateRing {
    slots: Vec<(u32, SixDofState)>,
    next: usize,
    count: usize,
}

impl StateRing {
    /// Create a ring with `capacity` slots.
    ///
    /// # Errors
    /// `ring-capacity-invalid` (zero or above [`MAX_RING_CAPACITY`],
    /// tested at cap and cap+1).
    pub fn new(capacity: usize) -> Result<StateRing, Refusal> {
        if capacity == 0 || capacity > MAX_RING_CAPACITY {
            return Err(Refusal {
                code: "ring-capacity-invalid",
                message: format!("capacity {capacity} outside [1, {MAX_RING_CAPACITY}]"),
                ranked_repairs: vec!["a few seconds of ticks is the intended window".into()],
            });
        }
        Ok(StateRing { slots: Vec::with_capacity(capacity), next: 0, count: 0 })
    }

    /// Push a tick's state (overwrites the oldest slot when full).
    pub fn push(&mut self, tick: u32, state: SixDofState) {
        if self.slots.len() < self.slots.capacity() {
            self.slots.push((tick, state));
        } else {
            self.slots[self.next] = (tick, state);
        }
        self.next = (self.next + 1) % self.slots.capacity();
        self.count += 1;
    }

    /// Fetch a tick's state if it is still inside the live window.
    #[must_use]
    pub fn get(&self, tick: u32) -> Option<&SixDofState> {
        self.slots.iter().find(|(t, _)| *t == tick).map(|(_, s)| s)
    }

    /// Number of pushes ever made (not the live count).
    #[must_use]
    pub fn pushes(&self) -> usize {
        self.count
    }
}

/// A recorded run: everything needed to reproduce the tick stream
/// bit-exactly on the same artifact (E3.2 DONE-WHEN clause).
#[derive(Clone, Debug, PartialEq)]
pub struct RunRecord {
    /// Initial state at tick 0.
    pub initial: SixDofState,
    /// Rigid-body properties.
    pub body: RigidBody,
    /// Fixed timestep [s].
    pub dt_s: f64,
    /// The closed applied-input trace.
    pub trace: InputTrace,
    /// The frozen trace id (recorded at close).
    pub input_trace_id: String,
    /// Per-tick state digests (the determinism receipt).
    pub tick_digests: Vec<String>,
}

/// Run the spine over an input trace: events set control-channel values at
/// the START of their applied tick (ordinal order); `loads` maps
/// (t, state, controls) to generalized loads. Returns the final state, the
/// digest trace, and a ring holding the trailing window.
///
/// # Errors
/// Trace admission refusals; `channel-outside-domain`; spine refusals.
pub fn run_recorded<F>(
    body: &RigidBody,
    initial: &SixDofState,
    dt_s: f64,
    trace: &InputTrace,
    ring_capacity: usize,
    mut loads: F,
) -> Result<(SixDofState, Vec<String>, StateRing), Refusal>
where
    F: FnMut(f64, &SixDofState, &[f64]) -> Loads,
{
    trace.admit()?;
    if trace.events.iter().any(|e| e.channel as usize >= MAX_CHANNELS) {
        return Err(Refusal {
            code: "channel-outside-domain",
            message: format!("channel above cap {MAX_CHANNELS}"),
            ranked_repairs: vec!["the 1903 Flyer has 2 control channels".into()],
        });
    }
    let mut ring = StateRing::new(ring_capacity)?;
    let mut controls = [0.0f64; MAX_CHANNELS];
    let mut state = *initial;
    let mut digests = Vec::with_capacity(trace.end_tick_exclusive as usize);
    let mut ev = trace.events.iter().peekable();
    for tick in 0..trace.end_tick_exclusive {
        while let Some(e) = ev.peek() {
            if e.applied_tick == tick {
                controls[e.channel as usize] = e.quantized_value;
                ev.next();
            } else {
                break;
            }
        }
        let t = f64::from(tick) * dt_s;
        let c = controls;
        state = step(body, &state, t, dt_s, |time, s| loads(time, s, &c))?;
        digests.push(tick_digest(tick, &state));
        ring.push(tick, state);
    }
    Ok((state, digests, ring))
}

/// Record a run (producing the closed trace id + digest receipt), then the
/// caller can [`verify_replay`] it on the same artifact.
///
/// # Errors
/// As [`run_recorded`].
pub fn record<F>(
    body: &RigidBody,
    initial: &SixDofState,
    dt_s: f64,
    trace: InputTrace,
    loads: F,
) -> Result<RunRecord, Refusal>
where
    F: FnMut(f64, &SixDofState, &[f64]) -> Loads,
{
    let (_, digests, _) = run_recorded(body, initial, dt_s, &trace, 64, loads)?;
    let input_trace_id = trace.trace_id();
    Ok(RunRecord { initial: *initial, body: *body, dt_s, trace, input_trace_id, tick_digests: digests })
}

/// Replay a record and verify bit-identity, LOCALIZING any divergence.
///
/// # Errors
/// `replay-trace-id-mismatch` (the record's id does not match its own
/// trace — tampering); `replay-digest-mismatch` naming the FIRST divergent
/// tick with expected/observed digests; length mismatches are typed too.
pub fn verify_replay<F>(record: &RunRecord, loads: F) -> Result<(), Refusal>
where
    F: FnMut(f64, &SixDofState, &[f64]) -> Loads,
{
    if record.trace.trace_id() != record.input_trace_id {
        return Err(Refusal {
            code: "replay-trace-id-mismatch",
            message: "the record's input trace does not hash to its recorded id".into(),
            ranked_repairs: vec!["the record was tampered with or corrupted; restore it".into()],
        });
    }
    let (_, digests, _) =
        run_recorded(&record.body, &record.initial, record.dt_s, &record.trace, 64, loads)?;
    if digests.len() != record.tick_digests.len() {
        return Err(Refusal {
            code: "replay-digest-mismatch",
            message: format!(
                "tick count diverged: replay {} vs record {}",
                digests.len(),
                record.tick_digests.len()
            ),
            ranked_repairs: vec!["the trace extent defines the count; check the record".into()],
        });
    }
    for (tick, (got, want)) in digests.iter().zip(&record.tick_digests).enumerate() {
        if got != want {
            return Err(Refusal {
                code: "replay-digest-mismatch",
                message: format!(
                    "first divergence at tick {tick}: expected {want}, observed {got}"
                ),
                ranked_repairs: vec![
                    "same-artifact replay must be bit-identical; a divergence here is a \
                     determinism defect or a changed load model"
                        .into(),
                ],
            });
        }
    }
    Ok(())
}
