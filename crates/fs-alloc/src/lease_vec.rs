//! Lease-admitted owned storage (bead wf9.16.1).
//!
//! [`LeasedVec`] is the sealed output-storage building block: a growable
//! owned buffer whose backing bytes are ADMITTED through an
//! [`OperationMemoryLease`] BEFORE the allocator is asked for them —
//! admission precedes allocation, so there is never an over-limit live
//! interval. The charge is an RAII [`LeaseCharge`] held inside the value:
//! it releases when the vector drops, including unwinds, with no manual
//! bookkeeping. Growth admits the new capacity in full before reallocating
//! (both buffers are live across the move, and the lease sees exactly
//! that), then drops the old charge.
//!
//! Receipts therefore track PAYLOAD bytes exactly; allocator bookkeeping
//! beyond the `try_reserve_exact`-sized buffer is an explicit no-claim.
//!
//! The detached-empty state ([`LeasedVec::detached_empty`]) exists for
//! contextless fold identities (`Reduce::identity()` has no lease): it owns
//! no storage, charges nothing, and REFUSES growth — an identity can be
//! merged from, never grown into.

use core::fmt;

use crate::arena::AllocError;
use crate::lease::{LeaseCharge, OperationMemoryLease};

/// A growable owned buffer whose capacity is lease-admitted before every
/// allocation. See the module docs for the admission discipline.
pub struct LeasedVec<T> {
    values: Vec<T>,
    /// RAII charge for the current buffer (`None` only when detached or
    /// zero-capacity).
    charge: Option<LeaseCharge>,
    /// `None` only for the detached-empty fold identity.
    lease: Option<OperationMemoryLease>,
    what: &'static str,
}

impl<T> LeasedVec<T> {
    /// Admit `capacity` elements against `lease`, then allocate exactly
    /// that capacity.
    ///
    /// # Errors
    /// [`AllocError::LeaseExhausted`] when the lease refuses;
    /// [`AllocError::OutOfMemory`] when the allocator refuses (the charge
    /// guard releases before returning); [`AllocError::LayoutOverflow`]
    /// when the byte size is unrepresentable.
    pub fn with_capacity(
        lease: &OperationMemoryLease,
        what: &'static str,
        capacity: usize,
    ) -> Result<Self, AllocError> {
        let (values, charge) = admitted_buffer::<T>(lease, what, capacity)?;
        Ok(LeasedVec {
            values,
            charge,
            lease: Some(lease.clone()),
            what,
        })
    }

    /// The contextless fold identity: owns nothing, charges nothing,
    /// refuses growth.
    #[must_use]
    pub fn detached_empty(what: &'static str) -> Self {
        LeasedVec {
            values: Vec::new(),
            charge: None,
            lease: None,
            what,
        }
    }

    /// Whether this is the detached fold identity.
    #[must_use]
    pub fn is_detached(&self) -> bool {
        self.lease.is_none()
    }

    /// Push within admitted capacity, or admit a geometric growth step
    /// first (admission precedes allocation; the old buffer's charge drops
    /// only after the move succeeds).
    ///
    /// # Errors
    /// As [`LeasedVec::with_capacity`]; growth on a detached identity is
    /// refused as `LeaseExhausted` with a zero limit.
    pub fn push(&mut self, value: T) -> Result<(), AllocError> {
        if self.values.len() == self.values.capacity() {
            let target = self.values.capacity().max(4).saturating_mul(2);
            self.regrow(target)?;
        }
        self.values.push(value);
        Ok(())
    }

    /// Move every element of `other` into `self`, admitting more capacity
    /// if needed. `other`'s charge releases when it drops here — the peak
    /// (both live) is the honest concurrent cost of a concatenating merge.
    ///
    /// # Errors
    /// As [`LeasedVec::push`].
    pub fn append(&mut self, mut other: LeasedVec<T>) -> Result<(), AllocError> {
        let needed = self.values.len().checked_add(other.values.len()).ok_or(
            AllocError::LayoutOverflow {
                site: self.what,
                len: usize::MAX,
                elem_bytes: size_of::<T>(),
            },
        )?;
        if needed > self.values.capacity() {
            self.regrow(needed)?;
        }
        self.values.append(&mut other.values);
        Ok(())
    }

    fn regrow(&mut self, target_capacity: usize) -> Result<(), AllocError> {
        let Some(lease) = self.lease.clone() else {
            return Err(AllocError::LeaseExhausted {
                site: self.what,
                requested_bytes: payload_bytes::<T>(self.what, target_capacity).unwrap_or(u64::MAX),
                used_bytes: 0,
                limit_bytes: 0,
            });
        };
        let (mut replacement, new_charge) =
            admitted_buffer::<T>(&lease, self.what, target_capacity)?;
        replacement.append(&mut self.values);
        self.values = replacement;
        // The old buffer is gone after the move; dropping its charge is the
        // release.
        self.charge = new_charge;
        Ok(())
    }

    /// Elements as a slice.
    #[must_use]
    pub fn as_slice(&self) -> &[T] {
        &self.values
    }

    /// Element count.
    #[must_use]
    pub fn len(&self) -> usize {
        self.values.len()
    }

    /// Whether no elements are stored.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.values.is_empty()
    }

    /// Admitted capacity in elements.
    #[must_use]
    pub fn capacity(&self) -> usize {
        self.values.capacity()
    }
}

/// Admit `capacity * size_of::<T>()` bytes, then allocate exactly that
/// buffer. The returned guard owns the admission; dropping it (on any later
/// failure or at the end of the buffer's life) is the release.
fn admitted_buffer<T>(
    lease: &OperationMemoryLease,
    what: &'static str,
    capacity: usize,
) -> Result<(Vec<T>, Option<LeaseCharge>), AllocError> {
    let bytes = payload_bytes::<T>(what, capacity)?;
    if bytes == 0 {
        return Ok((Vec::new(), None));
    }
    let charge = lease
        .reserve(what, bytes)
        .map_err(|refusal| AllocError::LeaseExhausted {
            site: what,
            requested_bytes: refusal.requested_bytes,
            used_bytes: refusal.used_bytes,
            limit_bytes: refusal.limit_bytes,
        })?;
    let mut values = Vec::new();
    if values.try_reserve_exact(capacity).is_err() {
        // `charge` drops here: admission released before the error returns.
        return Err(AllocError::OutOfMemory {
            site: what,
            requested_bytes: usize::try_from(bytes).unwrap_or(usize::MAX),
        });
    }
    Ok((values, Some(charge)))
}

impl<T: fmt::Debug> fmt::Debug for LeasedVec<T> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("LeasedVec")
            .field("len", &self.values.len())
            .field(
                "charged_bytes",
                &self.charge.as_ref().map_or(0, LeaseCharge::bytes),
            )
            .field("detached", &self.lease.is_none())
            .finish_non_exhaustive()
    }
}

impl<T: PartialEq> PartialEq for LeasedVec<T> {
    fn eq(&self, other: &Self) -> bool {
        self.values == other.values
    }
}

fn payload_bytes<T>(what: &'static str, capacity: usize) -> Result<u64, AllocError> {
    let bytes = capacity
        .checked_mul(size_of::<T>())
        .ok_or(AllocError::LayoutOverflow {
            site: what,
            len: capacity,
            elem_bytes: size_of::<T>(),
        })?;
    Ok(bytes as u64)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn admission_precedes_allocation_and_drop_releases() {
        let lease = OperationMemoryLease::bounded(1024);
        {
            let v = LeasedVec::<u64>::with_capacity(&lease, "t/payload", 64).expect("512 B fits");
            assert_eq!(lease.receipt().used_bytes, 512);
            assert_eq!(v.capacity(), 64);
        }
        assert_eq!(lease.receipt().used_bytes, 0);
        let refusal = LeasedVec::<u64>::with_capacity(&lease, "t/too-big", 200)
            .expect_err("1600 B over a 1024 B lease");
        assert!(
            matches!(
                refusal,
                AllocError::LeaseExhausted {
                    site: "t/too-big",
                    ..
                }
            ),
            "{refusal:?}"
        );
        assert_eq!(lease.receipt().used_bytes, 0, "refusal charges nothing");
    }

    #[test]
    fn growth_admits_before_reallocating_and_peaks_both_buffers() {
        let lease = OperationMemoryLease::bounded(4096);
        let mut v = LeasedVec::<u64>::with_capacity(&lease, "t/grow", 4).expect("fits");
        for i in 0..16u64 {
            v.push(i).expect("bounded growth");
        }
        assert_eq!(v.as_slice(), (0..16).collect::<Vec<_>>().as_slice());
        // Growth peaked with old+new live: 4→8: 32+64; 8→16: 64+128 = 192.
        assert!(lease.receipt().peak_bytes >= 192);
        drop(v);
        assert_eq!(lease.receipt().used_bytes, 0);
    }

    #[test]
    fn detached_identity_charges_nothing_and_refuses_growth() {
        let mut identity = LeasedVec::<u64>::detached_empty("t/identity");
        assert!(identity.is_detached());
        let error = identity.push(1).expect_err("identities never grow");
        assert!(matches!(
            error,
            AllocError::LeaseExhausted { limit_bytes: 0, .. }
        ));
    }

    #[test]
    fn append_merges_and_releases_the_source_charge() {
        let lease = OperationMemoryLease::bounded(8192);
        let mut left = LeasedVec::<u64>::with_capacity(&lease, "t/left", 8).expect("fits");
        let mut right = LeasedVec::<u64>::with_capacity(&lease, "t/right", 8).expect("fits");
        for i in 0..8u64 {
            left.push(i).expect("in capacity");
            right.push(100 + i).expect("in capacity");
        }
        assert_eq!(lease.receipt().used_bytes, 128, "two 64 B buffers live");
        left.append(right).expect("merge admits");
        assert_eq!(left.len(), 16);
        // Exactly the merged buffer remains charged: the source's 64 B and
        // the left's original 64 B both released; the admitted 16-element
        // replacement (128 B) is the whole live set.
        assert_eq!(lease.receipt().used_bytes, 128);
        // The merge's honest concurrent peak: left(64) + right(64) + the
        // replacement buffer(128) all live across the move.
        assert!(lease.receipt().peak_bytes >= 256);
        drop(left);
        assert_eq!(lease.receipt().used_bytes, 0);
    }

    #[test]
    fn unwind_releases_the_charge() {
        let lease = OperationMemoryLease::bounded(4096);
        let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            let _v = LeasedVec::<u8>::with_capacity(&lease, "t/unwind", 128).expect("fits");
            panic!("holder unwinds");
        }));
        assert!(result.is_err());
        assert_eq!(lease.receipt().used_bytes, 0);
    }

    #[test]
    fn zero_capacity_charges_nothing() {
        let lease = OperationMemoryLease::bounded(16);
        let v = LeasedVec::<u64>::with_capacity(&lease, "t/zero", 0).expect("free");
        assert_eq!(lease.receipt().used_bytes, 0);
        assert_eq!(v.capacity(), 0);
    }
}
