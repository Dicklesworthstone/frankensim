//! Run-scoped operation memory lease (bead wf9.16).
//!
//! One operation — a TilePool run and everything it allocates — obtains a
//! single [`OperationMemoryLease`] and every constituent mechanism charges
//! it: executor root metadata before worker launch, and every tile arena's
//! chunks while the operation holds them. The lease and [`ArenaPool`]'s
//! process-wide `limit_bytes` are DIFFERENT ledgers with different
//! lifetimes: the pool counts OS-reserved bytes (in-use + free-listed,
//! across operations), the lease counts one operation's live set. A chunk
//! recycled from the pool free list charges the acquiring operation's lease
//! exactly while held and never twice; free-list inventory belongs to no
//! operation. Both gates must admit; a refusal names whichever refused.
//!
//! Receipts are deterministic in structure and, for a fixed tile plan, in
//! their `requested` totals; `peak_bytes` is a CONSERVATIVE logical
//! high-water (reservation attempts count when entered), matching the
//! wf9.15 accounting doctrine. Thread stacks and allocator overhead are
//! explicitly NOT claimed (CONTRACT no-claims).
//!
//! [`ArenaPool`]: crate::ArenaPool

use core::fmt;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::{Arc, Mutex};

/// One refused lease reservation (the FIRST refusal is retained verbatim in
/// the receipt; later refusals only count).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LeaseRefusal {
    /// Component that requested the bytes (e.g. `"tilepool-root-metadata"`,
    /// `"arena-chunk"`).
    pub what: &'static str,
    /// Bytes the component asked for.
    pub requested_bytes: u64,
    /// Lease bytes in use at refusal time.
    pub used_bytes: u64,
    /// The lease limit in force.
    pub limit_bytes: u64,
}

impl fmt::Display for LeaseRefusal {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "operation memory lease refused {} B for `{}` with {} B of the {} B lease in use",
            self.requested_bytes, self.what, self.used_bytes, self.limit_bytes
        )
    }
}

/// Deterministic lease accounting snapshot.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LeaseReceipt {
    /// The limit in force (`None` = unbounded legacy wrapper).
    pub limit_bytes: Option<u64>,
    /// Cumulative bytes of granted reservations.
    pub requested_bytes: u64,
    /// Conservative logical high-water of concurrently held bytes.
    pub peak_bytes: u64,
    /// Bytes still held when the snapshot was taken.
    pub used_bytes: u64,
    /// Number of refused reservations.
    pub refusals: u64,
    /// The first refusal, verbatim.
    pub first_refusal: Option<LeaseRefusal>,
}

impl LeaseReceipt {
    /// Canonical JSON object (deterministic field order).
    #[must_use]
    pub fn to_json(&self) -> String {
        use std::fmt::Write as _;
        let mut out = String::from("{\"schema\":\"fs-alloc-operation-lease-v1\"");
        match self.limit_bytes {
            Some(limit) => {
                let _ = write!(out, ",\"limit_bytes\":{limit}");
            }
            None => out.push_str(",\"limit_bytes\":null"),
        }
        let _ = write!(
            out,
            ",\"requested_bytes\":{},\"peak_bytes\":{},\"used_bytes\":{},\"refusals\":{}",
            self.requested_bytes, self.peak_bytes, self.used_bytes, self.refusals
        );
        match &self.first_refusal {
            Some(refusal) => {
                let _ = write!(
                    out,
                    ",\"first_refusal\":{{\"what\":\"{}\",\"requested_bytes\":{},\"used_bytes\":{},\"limit_bytes\":{}}}",
                    refusal.what, refusal.requested_bytes, refusal.used_bytes, refusal.limit_bytes
                );
            }
            None => out.push_str(",\"first_refusal\":null"),
        }
        out.push('}');
        out
    }
}

struct LeaseShared {
    limit_bytes: Option<u64>,
    used_bytes: AtomicU64,
    peak_bytes: AtomicU64,
    requested_bytes: AtomicU64,
    refusals: AtomicU64,
    first_refusal: Mutex<Option<LeaseRefusal>>,
}

/// Cloneable run-scoped memory lease with atomic reserve/release. Clones
/// share one ledger; the value is cheap to hand to every worker.
#[derive(Clone)]
pub struct OperationMemoryLease {
    shared: Arc<LeaseShared>,
}

impl fmt::Debug for OperationMemoryLease {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("OperationMemoryLease")
            .field("receipt", &self.receipt().to_json())
            .finish_non_exhaustive()
    }
}

impl OperationMemoryLease {
    /// A lease enforcing a hard byte limit.
    #[must_use]
    pub fn bounded(limit_bytes: u64) -> Self {
        Self::with_limit(Some(limit_bytes))
    }

    /// The legacy-wrapper lease: accounts but never refuses.
    #[must_use]
    pub fn unbounded() -> Self {
        Self::with_limit(None)
    }

    fn with_limit(limit_bytes: Option<u64>) -> Self {
        OperationMemoryLease {
            shared: Arc::new(LeaseShared {
                limit_bytes,
                used_bytes: AtomicU64::new(0),
                peak_bytes: AtomicU64::new(0),
                requested_bytes: AtomicU64::new(0),
                refusals: AtomicU64::new(0),
                first_refusal: Mutex::new(None),
            }),
        }
    }

    /// The limit in force.
    #[must_use]
    pub fn limit_bytes(&self) -> Option<u64> {
        self.shared.limit_bytes
    }

    /// Atomically reserve `bytes` for `what`, returning a guard that
    /// releases on drop (panic containment: an unwinding holder releases
    /// its charge on the way out).
    ///
    /// # Errors
    /// [`LeaseRefusal`] when the reservation would exceed the limit; the
    /// refusal is also recorded in the receipt.
    pub fn reserve(
        &self,
        what: &'static str,
        bytes: u64,
    ) -> Result<LeaseCharge, LeaseRefusal> {
        self.try_reserve_raw(what, bytes)?;
        Ok(LeaseCharge {
            lease: self.clone(),
            bytes,
        })
    }

    /// Reserve without a guard; the caller owns the matching
    /// [`OperationMemoryLease::release_raw`]. The arena integration uses
    /// this because chunk lifetime is managed by `Arena::drop`.
    ///
    /// # Errors
    /// As [`OperationMemoryLease::reserve`].
    pub fn try_reserve_raw(&self, what: &'static str, bytes: u64) -> Result<(), LeaseRefusal> {
        loop {
            let used = self.shared.used_bytes.load(Ordering::Acquire);
            let Some(next) = used.checked_add(bytes) else {
                return Err(self.record_refusal(what, bytes, used));
            };
            if self
                .shared
                .limit_bytes
                .is_some_and(|limit| next > limit)
            {
                return Err(self.record_refusal(what, bytes, used));
            }
            if self
                .shared
                .used_bytes
                .compare_exchange_weak(used, next, Ordering::AcqRel, Ordering::Acquire)
                .is_ok()
            {
                self.shared.peak_bytes.fetch_max(next, Ordering::AcqRel);
                self.shared
                    .requested_bytes
                    .fetch_add(bytes, Ordering::Relaxed);
                return Ok(());
            }
        }
    }

    /// Release bytes previously reserved through
    /// [`OperationMemoryLease::try_reserve_raw`].
    pub fn release_raw(&self, bytes: u64) {
        let previous = self.shared.used_bytes.fetch_sub(bytes, Ordering::AcqRel);
        debug_assert!(
            previous >= bytes,
            "operation lease released more than it reserved"
        );
    }

    fn record_refusal(&self, what: &'static str, bytes: u64, used: u64) -> LeaseRefusal {
        let refusal = LeaseRefusal {
            what,
            requested_bytes: bytes,
            used_bytes: used,
            limit_bytes: self.shared.limit_bytes.unwrap_or(u64::MAX),
        };
        self.shared.refusals.fetch_add(1, Ordering::Relaxed);
        let mut first = self
            .shared
            .first_refusal
            .lock()
            .expect("fs-alloc lease refusal record poisoned");
        if first.is_none() {
            *first = Some(refusal.clone());
        }
        refusal
    }

    /// Deterministic accounting snapshot.
    #[must_use]
    pub fn receipt(&self) -> LeaseReceipt {
        LeaseReceipt {
            limit_bytes: self.shared.limit_bytes,
            requested_bytes: self.shared.requested_bytes.load(Ordering::Acquire),
            peak_bytes: self.shared.peak_bytes.load(Ordering::Acquire),
            used_bytes: self.shared.used_bytes.load(Ordering::Acquire),
            refusals: self.shared.refusals.load(Ordering::Acquire),
            first_refusal: self
                .shared
                .first_refusal
                .lock()
                .expect("fs-alloc lease refusal record poisoned")
                .clone(),
        }
    }
}

/// RAII lease charge: releases its bytes on drop (including unwinds).
#[derive(Debug)]
pub struct LeaseCharge {
    lease: OperationMemoryLease,
    bytes: u64,
}

impl LeaseCharge {
    /// Bytes held by this charge.
    #[must_use]
    pub fn bytes(&self) -> u64 {
        self.bytes
    }
}

impl Drop for LeaseCharge {
    fn drop(&mut self) {
        self.lease.release_raw(self.bytes);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reserve_release_and_receipt_account_exactly() {
        let lease = OperationMemoryLease::bounded(1000);
        let a = lease.reserve("root", 600).expect("fits");
        assert_eq!(lease.receipt().used_bytes, 600);
        let refusal = lease.reserve("chunk", 500).expect_err("over limit");
        assert_eq!(refusal.what, "chunk");
        assert_eq!(refusal.used_bytes, 600);
        assert_eq!(refusal.limit_bytes, 1000);
        let b = lease.reserve("chunk", 400).expect("exactly fits");
        drop(a);
        drop(b);
        let receipt = lease.receipt();
        assert_eq!(receipt.used_bytes, 0);
        assert_eq!(receipt.requested_bytes, 1000);
        assert_eq!(receipt.peak_bytes, 1000);
        assert_eq!(receipt.refusals, 1);
        let first = receipt.first_refusal.as_ref().expect("recorded");
        assert_eq!(first.what, "chunk");
        assert_eq!(first.requested_bytes, 500);
        assert!(receipt.to_json().contains("\"refusals\":1"));
    }

    #[test]
    fn unbounded_lease_accounts_but_never_refuses() {
        let lease = OperationMemoryLease::unbounded();
        let charge = lease
            .reserve("huge", u64::MAX / 2)
            .expect("unbounded admits");
        assert_eq!(lease.receipt().peak_bytes, u64::MAX / 2);
        drop(charge);
        assert_eq!(lease.receipt().used_bytes, 0);
        assert_eq!(lease.receipt().refusals, 0);
    }

    #[test]
    fn charges_release_on_unwind() {
        let lease = OperationMemoryLease::bounded(100);
        let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            let _held = lease.reserve("tile", 80).expect("fits");
            panic!("tile body panicked");
        }));
        assert!(result.is_err());
        assert_eq!(
            lease.receipt().used_bytes,
            0,
            "an unwinding holder must release its charge"
        );
    }

    #[test]
    fn concurrent_reservations_never_exceed_the_limit() {
        let lease = OperationMemoryLease::bounded(64);
        std::thread::scope(|s| {
            for _ in 0..8 {
                let lease = lease.clone();
                s.spawn(move || {
                    for _ in 0..200 {
                        if let Ok(charge) = lease.reserve("hammer", 16) {
                            assert!(lease.receipt().used_bytes <= 64);
                            drop(charge);
                        }
                    }
                });
            }
        });
        let receipt = lease.receipt();
        assert_eq!(receipt.used_bytes, 0);
        assert!(receipt.peak_bytes <= 64);
    }
}
