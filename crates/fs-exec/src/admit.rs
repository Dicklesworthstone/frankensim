//! Sealed lease-admitted output contract (bead wf9.16.1).
//!
//! The operation memory lease can only enforce a full-run live-set ceiling
//! if kernel OUTPUTS cannot smuggle heap past it: `size_of::<K::Out>()`
//! sees a `Vec`'s three words, never its payload. [`AdmittedStorage`] is
//! the sealed inductive contract that closes that hole — a type implements
//! it exactly when every byte it owns is either INLINE (visible to
//! `size_of`, charged by the pool's root-metadata slot accounting) or held
//! by lease-admitted storage ([`fs_alloc::LeasedVec`], whose buffer was
//! admitted before allocation). The trait is sealed: downstream crates
//! cannot declare a heap-bearing type admitted, so a hostile custom output
//! FAILS TO COMPILE on the leased production path (see the `compile_fail`
//! battery below).
//!
//! `Vec<T>` deliberately implements [`Reduce`] but NOT [`AdmittedStorage`]:
//! it remains available to the legacy unleased entries and is exactly the
//! bypass the leased bound refuses.
//!
//! ```compile_fail
//! // A hostile heap-bearing output cannot enter the admitted contract:
//! fn requires_admitted<T: fs_exec::AdmittedStorage>() {}
//! requires_admitted::<Vec<u8>>(); // ERROR: Vec is not admitted storage
//! ```
//!
//! ```compile_fail
//! // Sealing: downstream types cannot self-declare admission.
//! struct Smuggler(String);
//! impl fs_exec::AdmittedStorage for Smuggler {} // ERROR: sealed trait
//! ```

use crate::kernel::Reduce;
use fs_alloc::LeasedVec;

mod sealed {
    /// Sealing token: only fs-exec names the admitted set.
    pub trait Sealed {}
}

/// Every owned byte is inline or lease-admitted (see module docs). Sealed.
pub trait AdmittedStorage: Send + sealed::Sealed {}

macro_rules! admit_inline {
    ($($ty:ty),* $(,)?) => {
        $(
            impl sealed::Sealed for $ty {}
            impl AdmittedStorage for $ty {}
        )*
    };
}

admit_inline!(
    (),
    bool,
    char,
    u8,
    u16,
    u32,
    u64,
    u128,
    usize,
    i8,
    i16,
    i32,
    i64,
    i128,
    isize,
    f32,
    f64,
);

impl<E: AdmittedStorage, const N: usize> sealed::Sealed for [E; N] {}
impl<E: AdmittedStorage, const N: usize> AdmittedStorage for [E; N] {}

impl<A: AdmittedStorage, B: AdmittedStorage> sealed::Sealed for (A, B) {}
impl<A: AdmittedStorage, B: AdmittedStorage> AdmittedStorage for (A, B) {}
impl<A: AdmittedStorage, B: AdmittedStorage, C: AdmittedStorage> sealed::Sealed for (A, B, C) {}
impl<A: AdmittedStorage, B: AdmittedStorage, C: AdmittedStorage> AdmittedStorage for (A, B, C) {}
impl<A: AdmittedStorage, B: AdmittedStorage, C: AdmittedStorage, D: AdmittedStorage> sealed::Sealed
    for (A, B, C, D)
{
}
impl<A: AdmittedStorage, B: AdmittedStorage, C: AdmittedStorage, D: AdmittedStorage> AdmittedStorage
    for (A, B, C, D)
{
}

// The inductive step: a LeasedVec's buffer was admitted before allocation,
// and its elements are admitted by induction — nested outputs have
// explicit bounded storage all the way down.
impl<E: AdmittedStorage> sealed::Sealed for LeasedVec<E> {}
impl<E: AdmittedStorage> AdmittedStorage for LeasedVec<E> {}

/// Outputs eligible for the leased production path: mergeable AND
/// admission-visible. Blanket-derived; never implemented directly.
pub trait LeaseAdmittedOut: Reduce + AdmittedStorage {}
impl<T: Reduce + AdmittedStorage> LeaseAdmittedOut for T {}

/// Concatenating fold wrapper over lease-admitted element storage: the
/// list-shaped kernel output for the leased path (where `Vec<T>` is
/// refused).
///
/// `merge` re-admits capacity when the concatenation outgrows the left
/// side. `Reduce::merge` cannot return an error, so an admission refusal
/// there PANICS with a structured message — the pool contains fold-time
/// unwinds as [`crate::RunError::ReductionPanicked`], the outcome documented
/// for merge-time refusals; every lease charge releases during the unwind.
#[derive(Debug, PartialEq)]
pub struct Concat<E: AdmittedStorage>(pub LeasedVec<E>);

impl<E: AdmittedStorage> sealed::Sealed for Concat<E> {}
impl<E: AdmittedStorage> AdmittedStorage for Concat<E> {}

impl<E: AdmittedStorage> Reduce for Concat<E> {
    fn identity() -> Self {
        Concat(LeasedVec::detached_empty("fs-exec/concat-identity"))
    }

    fn merge(mut self, other: Self) -> Self {
        if self.0.is_detached() && self.0.is_empty() {
            return other;
        }
        match self.0.append(other.0) {
            Ok(()) => self,
            Err(error) => panic!(
                "leased concatenating fold refused mid-merge (reported as \
                 ReductionPanicked; charges release with the unwind): {error}"
            ),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_alloc::OperationMemoryLease;

    fn assert_admitted<T: LeaseAdmittedOut>() {}

    #[test]
    fn the_admitted_surface_compiles() {
        assert_admitted::<u64>();
        assert_admitted::<f64>();
        assert_admitted::<()>();
        assert_admitted::<Concat<u64>>();
        assert_admitted::<Concat<[f64; 4]>>();
        assert_admitted::<Concat<LeasedVec<u64>>>();
    }

    #[test]
    fn concat_identity_folds_and_merges_release_exactly() {
        let lease = OperationMemoryLease::bounded(1 << 16);
        let mut a = LeasedVec::with_capacity(&lease, "t/a", 4).expect("fits");
        let mut b = LeasedVec::with_capacity(&lease, "t/b", 4).expect("fits");
        for i in 0..4u64 {
            a.push(i).expect("in capacity");
            b.push(10 + i).expect("in capacity");
        }
        let merged = Concat::identity().merge(Concat(a)).merge(Concat(b));
        assert_eq!(merged.0.as_slice(), &[0, 1, 2, 3, 10, 11, 12, 13]);
        drop(merged);
        assert_eq!(lease.receipt().used_bytes, 0);
    }

    #[test]
    fn merge_admission_refusal_panics_for_pool_containment() {
        // Both sides fit alone; their concatenation cannot be admitted.
        let lease = OperationMemoryLease::bounded(96);
        let mut a = LeasedVec::with_capacity(&lease, "t/a", 6).expect("48 B fits");
        let mut b = LeasedVec::with_capacity(&lease, "t/b", 6).expect("48 B more fits");
        for i in 0..6u64 {
            a.push(i).expect("in capacity");
            b.push(i).expect("in capacity");
        }
        let unwound = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            let _ = Concat(a).merge(Concat(b));
        }));
        assert!(unwound.is_err(), "refused merge must unwind for the pool");
        assert_eq!(
            lease.receipt().used_bytes,
            0,
            "every charge releases with the unwind"
        );
    }
}
