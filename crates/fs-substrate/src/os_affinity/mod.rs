//! OS thread-affinity + page-placement queries (bead fz2.2): raw Linux
//! syscalls via `asm!` — no libc (Decalogue P1) — behind a safe façade.
//!
//! CAPSULE: registered in unsafe-capsules.json; SAFETY.md beside this
//! file. Linux x86-64 only — every other target gets the same API
//! returning [`OsAffinityError::Unsupported`] (an honest refusal, never
//! a silent no-op: a caller that thinks it pinned but didn't would
//! ledger fake locality claims). macOS P/E pinning needs QoS APIs
//! outside the dependency policy — the documented no-claim in
//! `bandwidth.rs` stands; the portable Apple mechanism is work-stealing
//! (fs-la xlvx dispenser, fs-exec pool).
//!
//! Pinning here is ADVISORY infrastructure for A/B locality harnesses
//! and CCD-aware scheduling; correctness of results never depends on
//! it (determinism discipline P2 — placement changes timing, not bits).

#![allow(unsafe_code)] // capsule: registered in unsafe-capsules.json, SAFETY.md beside

use core::fmt;

/// Why an affinity operation could not be performed.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OsAffinityError {
    /// The target OS/arch has no supported implementation.
    Unsupported(&'static str),
    /// The syscall itself failed (negative errno returned).
    Syscall {
        /// Which call failed.
        call: &'static str,
        /// The (positive) errno value.
        errno: i64,
    },
    /// A caller-side argument problem (empty CPU list, oversized ids).
    BadArgument(&'static str),
}

impl fmt::Display for OsAffinityError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            OsAffinityError::Unsupported(why) => {
                write!(f, "os-affinity unsupported on this target: {why}")
            }
            OsAffinityError::Syscall { call, errno } => {
                write!(f, "{call} failed with errno {errno}")
            }
            OsAffinityError::BadArgument(what) => write!(f, "bad argument: {what}"),
        }
    }
}

impl core::error::Error for OsAffinityError {}

/// Pin the CURRENT thread to the given logical CPUs.
///
/// # Errors
/// [`OsAffinityError`] — unsupported target, empty/oversized CPU list,
/// or the raw `sched_setaffinity` failure with errno.
pub fn pin_current_thread(cpus: &[u32]) -> Result<(), OsAffinityError> {
    imp::pin_current_thread(cpus)
}

/// The logical CPU the current thread is executing on right now.
///
/// # Errors
/// [`OsAffinityError`] — unsupported target or `getcpu` failure.
pub fn current_cpu() -> Result<u32, OsAffinityError> {
    imp::current_cpu()
}

/// NUMA node of each page backing `buf` (first-touch AUDIT: the
/// `move_pages` QUERY form — nodes=NULL moves nothing). Pages not yet
/// faulted in report `-ENOENT` (-2): audit AFTER touching.
///
/// # Errors
/// [`OsAffinityError`] — unsupported target or syscall failure.
pub fn page_nodes(buf: &[u8], page_size: usize) -> Result<Vec<i32>, OsAffinityError> {
    imp::page_nodes(buf, page_size)
}

/// Advise transparent hugepages for the page-aligned INTERIOR of
/// `buf` (`madvise(MADV_HUGEPAGE)`); returns the advised byte count
/// (0 when the aligned interior is empty). Only effective when the
/// kernel's THP mode is `madvise` or `always` — read the mode via
/// [`crate::affinity::thp_mode`] and ledger it with any measurement.
///
/// # Errors
/// [`OsAffinityError`] — unsupported target or syscall failure.
pub fn advise_hugepages(buf: &[u8], page_size: usize) -> Result<usize, OsAffinityError> {
    imp::advise_hugepages(buf.as_ptr() as usize, buf.len(), page_size)
}

/// [`advise_hugepages`] for a `u64` buffer (no unsafe byte-casting at
/// call sites): advise BEFORE first write — THP coalesces on fault,
/// and `vec![0u64; n]` maps lazily, so advising the fresh zeroed
/// allocation and then writing it faults huge pages.
///
/// # Errors
/// [`OsAffinityError`] — unsupported target or syscall failure.
pub fn advise_hugepages_words(buf: &[u64], page_size: usize) -> Result<usize, OsAffinityError> {
    imp::advise_hugepages(buf.as_ptr() as usize, buf.len() * 8, page_size)
}

#[cfg(all(target_os = "linux", target_arch = "x86_64"))]
mod imp {
    use super::OsAffinityError;
    use core::arch::asm;

    const SYS_SCHED_SETAFFINITY: u64 = 203;
    const SYS_MOVE_PAGES: u64 = 279;
    const SYS_GETCPU: u64 = 309;

    /// One raw Linux syscall (up to 6 args), x86-64 convention.
    ///
    /// SAFETY: caller guarantees the argument registers hold values
    /// valid for the specific syscall number (pointers live and sized
    /// as the kernel expects for the duration of the call).
    unsafe fn syscall6(nr: u64, a1: u64, a2: u64, a3: u64, a4: u64, a5: u64, a6: u64) -> i64 {
        let ret: i64;
        // SAFETY: the Linux x86-64 syscall ABI clobbers only rcx/r11
        // (declared); all pointer arguments obey the caller contract.
        unsafe {
            asm!(
                "syscall",
                inlateout("rax") nr => ret,
                in("rdi") a1,
                in("rsi") a2,
                in("rdx") a3,
                in("r10") a4,
                in("r8") a5,
                in("r9") a6,
                lateout("rcx") _,
                lateout("r11") _,
                options(nostack)
            );
        }
        ret
    }

    pub fn pin_current_thread(cpus: &[u32]) -> Result<(), OsAffinityError> {
        if cpus.is_empty() {
            return Err(OsAffinityError::BadArgument("empty CPU list"));
        }
        // 1024-bit cpu_set_t: sixteen u64 words.
        let mut mask = [0u64; 16];
        for &c in cpus {
            let (word, bit) = ((c / 64) as usize, c % 64);
            if word >= mask.len() {
                return Err(OsAffinityError::BadArgument("cpu id >= 1024"));
            }
            mask[word] |= 1u64 << bit;
        }
        // SAFETY: pid 0 = current thread; the mask pointer/len name a
        // live stack array for the (synchronous) call's duration.
        let ret = unsafe {
            syscall6(
                SYS_SCHED_SETAFFINITY,
                0,
                core::mem::size_of_val(&mask) as u64,
                mask.as_ptr() as u64,
                0,
                0,
                0,
            )
        };
        if ret < 0 {
            return Err(OsAffinityError::Syscall {
                call: "sched_setaffinity",
                errno: -ret,
            });
        }
        Ok(())
    }

    pub fn current_cpu() -> Result<u32, OsAffinityError> {
        let mut cpu: u32 = 0;
        let mut node: u32 = 0;
        // SAFETY: both out-pointers name live stack slots; tcache=NULL
        // is the documented modern form.
        let ret = unsafe {
            syscall6(
                SYS_GETCPU,
                core::ptr::from_mut(&mut cpu) as u64,
                core::ptr::from_mut(&mut node) as u64,
                0,
                0,
                0,
                0,
            )
        };
        if ret < 0 {
            return Err(OsAffinityError::Syscall {
                call: "getcpu",
                errno: -ret,
            });
        }
        Ok(cpu)
    }

    pub fn advise_hugepages(
        base: usize,
        len: usize,
        page_size: usize,
    ) -> Result<usize, OsAffinityError> {
        const SYS_MADVISE: u64 = 28;
        const MADV_HUGEPAGE: u64 = 14;
        if len == 0 || page_size == 0 {
            return Ok(0);
        }
        let start = base.next_multiple_of(page_size);
        let end = (base + len) / page_size * page_size;
        if start >= end {
            return Ok(0);
        }
        // SAFETY: [start, end) is a page-aligned sub-range of the live
        // caller-owned buffer; MADV_HUGEPAGE is purely advisory and
        // never invalidates or discards memory contents.
        let ret = unsafe {
            syscall6(
                SYS_MADVISE,
                start as u64,
                (end - start) as u64,
                MADV_HUGEPAGE,
                0,
                0,
                0,
            )
        };
        if ret < 0 {
            return Err(OsAffinityError::Syscall {
                call: "madvise",
                errno: -ret,
            });
        }
        Ok(end - start)
    }

    pub fn page_nodes(buf: &[u8], page_size: usize) -> Result<Vec<i32>, OsAffinityError> {
        if buf.is_empty() || page_size == 0 {
            return Ok(Vec::new());
        }
        let base = buf.as_ptr() as usize;
        let first = base / page_size * page_size;
        let count = (base + buf.len()).div_ceil(page_size) - first / page_size;
        let pages: Vec<u64> = (0..count).map(|i| (first + i * page_size) as u64).collect();
        let mut status = vec![i32::MIN; count];
        // SAFETY: query form — nodes=NULL moves nothing; pages/status
        // are live, correctly sized arrays for the synchronous call.
        let ret = unsafe {
            syscall6(
                SYS_MOVE_PAGES,
                0,
                count as u64,
                pages.as_ptr() as u64,
                0,
                status.as_mut_ptr() as u64,
                0,
            )
        };
        if ret < 0 {
            return Err(OsAffinityError::Syscall {
                call: "move_pages",
                errno: -ret,
            });
        }
        Ok(status)
    }
}

#[cfg(not(all(target_os = "linux", target_arch = "x86_64")))]
mod imp {
    use super::OsAffinityError;

    const WHY: &str = "raw-syscall affinity capsule implements Linux x86-64 only \
                       (macOS pinning needs QoS APIs outside the dependency policy)";

    pub fn pin_current_thread(_cpus: &[u32]) -> Result<(), OsAffinityError> {
        Err(OsAffinityError::Unsupported(WHY))
    }

    pub fn current_cpu() -> Result<u32, OsAffinityError> {
        Err(OsAffinityError::Unsupported(WHY))
    }

    pub fn page_nodes(_buf: &[u8], _page_size: usize) -> Result<Vec<i32>, OsAffinityError> {
        Err(OsAffinityError::Unsupported(WHY))
    }

    pub fn advise_hugepages(
        _base: usize,
        _len: usize,
        _page_size: usize,
    ) -> Result<usize, OsAffinityError> {
        Err(OsAffinityError::Unsupported(WHY))
    }
}
