//! Software-prefetch capsule (bead fz2.2): one total, safe hint —
//! `read_ahead(slice, idx)` — behind which each ISA issues its L1
//! read-prefetch (`prefetcht0` on x86-64, `prfm pldl1keep` on
//! aarch64; a documented no-op elsewhere). Registered capsule;
//! SAFETY.md beside this file.
//!
//! Prefetch is architecturally side-effect-free on memory state: it
//! never faults, never writes, and never changes any value — it can
//! only move cache lines. Hence the façade is total (out-of-range
//! indices are a no-op, not a panic) and the hint can NEVER affect
//! results, only timing (determinism P2 by construction).
//!
//! WHERE IT PAYS: data-DEPENDENT access (gather/scatter through index
//! arrays — SpMV columns, unstructured-mesh neighbors) that hardware
//! prefetchers cannot predict. Sequential streams are already covered
//! by hardware; the sweep harness measures the distance per machine
//! and the tuner records it.
#![allow(unsafe_code)] // registered capsule — see SAFETY.md beside this file

/// Hint the cache that `slice[idx]` will be READ soon. Total: an
/// out-of-range `idx` does nothing. Never affects results — only
/// residency (timing).
#[inline]
pub fn read_ahead<T>(slice: &[T], idx: usize) {
    if idx >= slice.len() {
        return;
    }
    // A pointer inside a live slice; the instruction dereferences
    // NOTHING architecturally.
    let p = unsafe { slice.as_ptr().add(idx) };
    imp::l1_read(p.cast());
}

#[cfg(target_arch = "x86_64")]
mod imp {
    /// # Safety-free by contract
    /// `prefetcht0` has no architectural effect on memory or flags; a
    /// stale/invalid address cannot fault through it.
    #[inline]
    pub fn l1_read(p: *const u8) {
        // SAFETY: prefetch never dereferences; any address is allowed
        // architecturally (invalid ones are simply ignored by the CPU).
        unsafe {
            core::arch::x86_64::_mm_prefetch::<{ core::arch::x86_64::_MM_HINT_T0 }>(p.cast());
        }
    }
}

#[cfg(target_arch = "aarch64")]
mod imp {
    #[inline]
    pub fn l1_read(p: *const u8) {
        // SAFETY: PRFM PLDL1KEEP has no architectural effect on memory
        // or flags and cannot fault; the register holds a plain address.
        unsafe {
            core::arch::asm!(
                "prfm pldl1keep, [{p}]",
                p = in(reg) p,
                options(nostack, preserves_flags)
            );
        }
    }
}

#[cfg(not(any(target_arch = "x86_64", target_arch = "aarch64")))]
mod imp {
    #[inline]
    pub fn l1_read(_p: *const u8) {
        // Documented no-op: no prefetch instruction wired for this ISA.
    }
}

#[cfg(test)]
mod tests {
    use super::read_ahead;

    #[test]
    fn total_and_effect_free_on_values() {
        let v: Vec<u64> = (0..1024).collect();
        for i in [0usize, 1, 511, 1023, 1024, 10_000] {
            read_ahead(&v, i); // in and out of range: total
        }
        read_ahead::<u64>(&[], 0);
        // The hint can never change data.
        assert!(v.iter().enumerate().all(|(i, &x)| x == i as u64));
    }
}
