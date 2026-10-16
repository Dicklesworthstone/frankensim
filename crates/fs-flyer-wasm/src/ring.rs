//! Versioned worker ABI + leased FieldSourceSnapshotV1 ring (bead
//! wf-root-guzez.6.1, E5.0). Plan Round-5: the ONLY state channel from
//! the sim worker to consumers is a leased ring of ≥3 snapshot slots
//! with a seqlock header per slot:
//!
//!   abi_version · payload_layout_hash · payload_bytes · run_epoch ·
//!   run_anchor_digest_prefix · model_id_prefix · tick · slot · sequence
//!
//! Writer: begin (sequence → odd) → fill → end (sequence → even).
//! Reader: lease the newest EVEN slot, copy, revalidate sequence —
//! any change means a torn read and the lease is refused, counted,
//! and retried by the caller. A run restart bumps `run_epoch`, which
//! invalidates every outstanding lease (the ABA guard: a slot reused
//! for a new run can never satisfy an old lease). Drop counters record
//! every unread overwrite — starvation is visible, never silent.
//!
//! This module is the PROTOCOL in host-testable form (deterministic
//! single-threaded interleavings drive the battery). The cross-thread
//! memory-ordering execution of the same protocol is browser-lane
//! scope (E0.6b/E6.4) and is recorded as such — the state machine,
//! validation rules, and refusal surface here are what those builds
//! reuse verbatim.

use crate::Refusal;
use fs_blake3::hash_domain;

/// ABI version (bump = renegotiate; mismatch is a typed refusal).
pub const ABI_VERSION: u32 = 1;

/// Slot-count window (≥3 guarantees writer progress under one lease).
pub const MIN_SLOTS: usize = 3;
/// Slot-count cap.
pub const MAX_SLOTS: usize = 8;

/// Payload float-count cap (absurd-input guard).
pub const MAX_PAYLOAD_F64: usize = 4096;

/// The FieldSourceSnapshotV1 payload layout descriptor (versioned; its
/// hash prefix rides in every header).
pub const PAYLOAD_LAYOUT_V1: &str =
    "org.frankensim.wf.field-source-snapshot.v1:[pos3,quat4,vel3,omega3,dc,dw,omega_prop,ou...]";

/// Layout hash prefix (u64 of the blake3 of the descriptor).
#[must_use]
pub fn payload_layout_hash() -> u64 {
    let hex = hash_domain(
        "org.frankensim.wf.snapshot-layout.v1",
        PAYLOAD_LAYOUT_V1.as_bytes(),
    )
    .to_hex();
    u64::from_str_radix(&hex[..16], 16).unwrap_or(0)
}

/// One slot's header (the ABI words).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct SlotHeader {
    /// ABI version.
    pub abi_version: u32,
    /// Payload layout identity prefix.
    pub payload_layout_hash: u64,
    /// Payload length in f64 words.
    pub payload_len: u32,
    /// Run epoch (bumped on restart — the ABA guard).
    pub run_epoch: u64,
    /// Tick-0 digest prefix (binds the ring to the equilibrated run).
    pub run_anchor_digest_prefix: u64,
    /// Model identity prefix.
    pub model_id_prefix: u64,
    /// Physics tick this snapshot is OF.
    pub tick: u64,
    /// Slot index.
    pub slot: u32,
    /// Seqlock sequence (even = stable, odd = being written).
    pub sequence: u64,
}

struct Slot {
    header: SlotHeader,
    payload: Vec<f64>,
}

/// A reader lease: slot + the sequence/epoch it validated against.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Lease {
    /// Leased slot.
    pub slot: u32,
    /// Sequence observed at acquisition (even).
    pub sequence: u64,
    /// Epoch observed at acquisition.
    pub run_epoch: u64,
}

/// Protocol counters (receipts — starvation and tearing are visible).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct RingCounters {
    /// Snapshots published.
    pub published: u64,
    /// Slots overwritten before any reader saw them.
    pub dropped_unread: u64,
    /// Torn-read refusals handed to readers.
    pub torn_reads: u64,
    /// Stale-epoch lease refusals.
    pub stale_epoch_refusals: u64,
    /// Writer skips because the target slot was leased.
    pub writer_skips: u64,
}

/// The ring.
pub struct SnapshotRing {
    slots: Vec<Slot>,
    read_flags: Vec<bool>,
    leased: Vec<bool>,
    newest: Option<u32>,
    run_epoch: u64,
    run_anchor_digest_prefix: u64,
    model_id_prefix: u64,
    payload_len: u32,
    counters: RingCounters,
}

impl SnapshotRing {
    /// Create a ring bound to a run identity.
    ///
    /// # Errors
    /// `ring-config-invalid` (slots outside [MIN, MAX] — both edges AND
    /// one past each; payload length 0 or above the cap).
    pub fn new(
        n_slots: usize,
        payload_len: usize,
        run_anchor_digest_prefix: u64,
        model_id_prefix: u64,
    ) -> Result<SnapshotRing, Refusal> {
        if !(MIN_SLOTS..=MAX_SLOTS).contains(&n_slots)
            || payload_len == 0
            || payload_len > MAX_PAYLOAD_F64
        {
            return Err(Refusal {
                code: "ring-config-invalid",
                message: format!("{n_slots} slots, payload {payload_len}"),
                ranked_repairs: vec![format!(
                    "slots in [{MIN_SLOTS}, {MAX_SLOTS}]; payload in [1, {MAX_PAYLOAD_F64}]"
                )],
            });
        }
        let layout = payload_layout_hash();
        let slots = (0..n_slots)
            .map(|i| Slot {
                header: SlotHeader {
                    abi_version: ABI_VERSION,
                    payload_layout_hash: layout,
                    payload_len: payload_len as u32,
                    run_epoch: 1,
                    run_anchor_digest_prefix,
                    model_id_prefix,
                    tick: 0,
                    slot: i as u32,
                    sequence: 0,
                },
                payload: vec![0.0; payload_len],
            })
            .collect();
        Ok(SnapshotRing {
            slots,
            read_flags: vec![true; n_slots], // nothing to drop initially
            leased: vec![false; n_slots],
            newest: None,
            run_epoch: 1,
            run_anchor_digest_prefix,
            model_id_prefix,
            payload_len: payload_len as u32,
            counters: RingCounters::default(),
        })
    }

    /// Publish a snapshot for `tick`. Chooses the oldest unleased slot;
    /// counts drops of unread slots and skips of leased ones.
    ///
    /// # Errors
    /// `ring-publish-invalid` (payload length mismatch, non-finite
    /// entries, non-monotonic tick).
    pub fn publish(&mut self, tick: u64, payload: &[f64]) -> Result<u32, Refusal> {
        if payload.len() != self.payload_len as usize || payload.iter().any(|v| !v.is_finite()) {
            return Err(Refusal {
                code: "ring-publish-invalid",
                message: format!(
                    "payload len {} (want {}) or non-finite",
                    payload.len(),
                    self.payload_len
                ),
                ranked_repairs: vec!["publish the declared layout exactly".into()],
            });
        }
        if let Some(n) = self.newest {
            let newest_tick = self.slots[n as usize].header.tick;
            if tick <= newest_tick {
                return Err(Refusal {
                    code: "ring-publish-invalid",
                    message: format!("tick {tick} not after newest {newest_tick}"),
                    ranked_repairs: vec!["ticks are strictly monotone within an epoch".into()],
                });
            }
        }
        // Oldest unleased slot (deterministic scan by tick, then index).
        let mut target: Option<usize> = None;
        for (i, s) in self.slots.iter().enumerate() {
            if self.leased[i] {
                self.counters.writer_skips += 1;
                continue;
            }
            match target {
                None => target = Some(i),
                Some(t) => {
                    if s.header.tick < self.slots[t].header.tick {
                        target = Some(i);
                    }
                }
            }
        }
        let Some(i) = target else {
            return Err(Refusal {
                code: "ring-publish-invalid",
                message: "all slots leased (protocol violation with >=3 slots)".into(),
                ranked_repairs: vec!["a reader is hoarding leases".into()],
            });
        };
        if !self.read_flags[i] {
            self.counters.dropped_unread += 1;
        }
        // Seqlock write: odd, mutate, even.
        let s = &mut self.slots[i];
        s.header.sequence += 1; // odd = writing
        s.payload.copy_from_slice(payload);
        s.header.tick = tick;
        s.header.run_epoch = self.run_epoch;
        s.header.run_anchor_digest_prefix = self.run_anchor_digest_prefix;
        s.header.sequence += 1; // even = stable
        self.read_flags[i] = false;
        self.newest = Some(i as u32);
        self.counters.published += 1;
        Ok(i as u32)
    }

    /// Acquire a lease on the newest stable slot.
    ///
    /// # Errors
    /// `ring-empty` (nothing published this epoch);
    /// `ring-abi-mismatch` (caller's expected version differs).
    pub fn acquire(&mut self, expected_abi: u32) -> Result<Lease, Refusal> {
        if expected_abi != ABI_VERSION {
            return Err(Refusal {
                code: "ring-abi-mismatch",
                message: format!("caller expects ABI {expected_abi}, ring is {ABI_VERSION}"),
                ranked_repairs: vec!["renegotiate the worker ABI version".into()],
            });
        }
        let Some(n) = self.newest else {
            return Err(Refusal {
                code: "ring-empty",
                message: "no snapshot published this epoch".into(),
                ranked_repairs: vec!["wait for the first publish".into()],
            });
        };
        let i = n as usize;
        let h = self.slots[i].header;
        debug_assert!(h.sequence % 2 == 0);
        self.leased[i] = true;
        self.read_flags[i] = true;
        Ok(Lease {
            slot: n,
            sequence: h.sequence,
            run_epoch: h.run_epoch,
        })
    }

    /// Read under a lease: copies the payload out, REVALIDATING the
    /// sequence and epoch — a torn or stale read is a typed refusal and
    /// is counted.
    ///
    /// # Errors
    /// `ring-lease-torn` (sequence moved — writer touched the slot);
    /// `ring-epoch-stale` (run restarted since the lease — ABA guard).
    pub fn read(&mut self, lease: &Lease, out: &mut [f64]) -> Result<SlotHeader, Refusal> {
        let i = lease.slot as usize;
        let h = self.slots[i].header;
        if lease.run_epoch != self.run_epoch || h.run_epoch != lease.run_epoch {
            self.counters.stale_epoch_refusals += 1;
            return Err(Refusal {
                code: "ring-epoch-stale",
                message: format!(
                    "lease epoch {} vs ring epoch {} — the run restarted",
                    lease.run_epoch, self.run_epoch
                ),
                ranked_repairs: vec!["drop the lease and re-acquire".into()],
            });
        }
        if h.sequence != lease.sequence {
            self.counters.torn_reads += 1;
            return Err(Refusal {
                code: "ring-lease-torn",
                message: format!(
                    "slot {} sequence {} vs leased {} — writer reused the slot",
                    lease.slot, h.sequence, lease.sequence
                ),
                ranked_repairs: vec!["re-acquire the newest slot".into()],
            });
        }
        if out.len() != self.payload_len as usize {
            return Err(Refusal {
                code: "ring-publish-invalid",
                message: format!(
                    "reader buffer {} vs payload {}",
                    out.len(),
                    self.payload_len
                ),
                ranked_repairs: vec!["size the reader buffer to the declared layout".into()],
            });
        }
        out.copy_from_slice(&self.slots[i].payload);
        Ok(h)
    }

    /// Release a lease (idempotent).
    pub fn release(&mut self, lease: &Lease) {
        let i = lease.slot as usize;
        if i < self.leased.len() {
            self.leased[i] = false;
        }
    }

    /// Restart the run: epoch bump + new anchor. Every outstanding
    /// lease becomes stale (ABA guard); slots are reset unpublished.
    pub fn restart(&mut self, run_anchor_digest_prefix: u64) {
        self.run_epoch += 1;
        self.run_anchor_digest_prefix = run_anchor_digest_prefix;
        self.newest = None;
        for (i, s) in self.slots.iter_mut().enumerate() {
            s.header.run_epoch = self.run_epoch;
            s.header.tick = 0;
            self.read_flags[i] = true;
        }
    }

    /// Counters snapshot (receipts).
    #[must_use]
    pub fn counters(&self) -> RingCounters {
        self.counters
    }

    /// Current epoch.
    #[must_use]
    pub fn epoch(&self) -> u64 {
        self.run_epoch
    }

    /// Simulate a writer overrun for the battery's torn-read twin: two
    /// full publishes advance the leased slot's sequence, so the held
    /// lease MUST fail revalidation. (Deterministic interleaving driver
    /// — the cross-thread execution is browser-lane scope.)
    #[doc(hidden)]
    pub fn force_slot_rewrite(&mut self, slot: u32, tick: u64) {
        let i = slot as usize;
        self.leased[i] = false; // simulate the misbehaving writer path
        let s = &mut self.slots[i];
        s.header.sequence += 1;
        s.header.tick = tick;
        s.header.sequence += 1;
    }
}
