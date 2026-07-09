//! SME2 streaming-mode GEMM prototype (bead wf9.3, feature
//! `frontier-sme2`, [F]): outer-product accumulation on the ZA tile.
//! EXPLORATORY tier — never load-bearing, never in the [`crate::ops`]
//! table, runtime-capability-gated (NEVER compile-time assumed), and
//! inert wherever the hardware or OS support is absent: NEON remains
//! the committed path.
//!
//! Shape: fixed 16×16 f32 microkernel (requires SVL = 512 bits — the
//! Apple M4 class; other streaming vector lengths report unavailable
//! rather than guessing). Packing/blocking integration with fs-la-gemm
//! is the xdgf-side successor; this module carries its own panel
//! packing so only the microkernel would swap.
//!
//! Determinism: per-element accumulation order equals the scalar
//! twin's k-order with fused multiply-add — the equivalence battery
//! measures whether that yields BITWISE equality on this hardware and
//! ledgers the answer; NO cross-ISA determinism-mode claim is made
//! until the G5 report characterizes it (the bead's explicit
//! non-goal).

#![allow(unsafe_code)]

use std::sync::OnceLock;

/// Tile dimension of the fixed-shape prototype (SVL 512 = 16 f32).
pub const TILE: usize = 16;

/// Runtime capability probe: SME2 present AND the streaming vector
/// length matches the prototype's fixed shape. macOS: `sysctl`
/// subprocess (the fs-substrate P1 pattern — no FFI); Linux/aarch64:
/// `/proc/cpuinfo` features; elsewhere: false.
#[must_use]
pub fn sme2_available() -> bool {
    static AVAILABLE: OnceLock<bool> = OnceLock::new();
    *AVAILABLE.get_or_init(|| probe_os() && streaming_vl_bytes() == Some(64))
}

fn probe_os() -> bool {
    #[cfg(target_os = "macos")]
    {
        let out = std::process::Command::new("sysctl")
            .args(["-n", "hw.optional.arm.FEAT_SME2"])
            .output();
        matches!(out, Ok(o) if String::from_utf8_lossy(&o.stdout).trim() == "1")
    }
    #[cfg(all(target_os = "linux", target_arch = "aarch64"))]
    {
        std::fs::read_to_string("/proc/cpuinfo")
            .map(|s| {
                s.lines()
                    .any(|l| l.starts_with("Features") && l.contains(" sme2"))
            })
            .unwrap_or(false)
    }
    #[cfg(not(any(target_os = "macos", all(target_os = "linux", target_arch = "aarch64"))))]
    {
        false
    }
}

/// The streaming vector length in bytes (rdsvl — readable without
/// entering streaming mode), or None when the OS probe already said no.
#[must_use]
pub fn streaming_vl_bytes() -> Option<u64> {
    if !probe_os() {
        return None;
    }
    let svl: u64;
    // SAFETY: `rdsvl` only reads the streaming vector length register;
    // it is architecturally defined whenever FEAT_SME is present,
    // which `probe_os` just confirmed. No state is modified.
    unsafe {
        core::arch::asm!(
            ".arch armv9.2-a+sme2",
            "rdsvl {out}, #1",
            out = out(reg) svl,
            options(nomem, nostack, pure)
        );
    }
    Some(svl)
}

/// The streaming-mode 16×16 GEMM microkernel: C[16×16] = Σ_k
/// a_col(k) ⊗ b_row(k). `a_panel` is column-major per k (k×16),
/// `b_panel` row-major per k (k×16), `c` row-major 16×16 with row
/// stride 16.
///
/// # Panics
/// If SME2 is unavailable (callers gate on [`sme2_available`]) or the
/// panel lengths disagree with `k`.
pub fn gemm_tile_f32(a_panel: &[f32], b_panel: &[f32], c: &mut [f32], k: usize) {
    assert!(sme2_available(), "SME2 tile kernel needs the capability");
    assert!(k >= 1, "empty accumulation");
    assert_eq!(a_panel.len(), k * TILE, "a panel is k x 16");
    assert_eq!(b_panel.len(), k * TILE, "b panel is k x 16");
    assert_eq!(c.len(), TILE * TILE, "c is 16 x 16");
    // SAFETY: one self-contained streaming-mode region. `smstart`
    // enters streaming SVE + ZA and `smstop` leaves before the block
    // ends, so no Rust FP/SIMD code interleaves with streaming state.
    // All v0–v31 are declared clobbered (z-registers alias them);
    // rustc allocates neither predicate nor ZA registers, and the
    // block touches only p0/p1 and za0.s beyond the declared operands.
    // Pointer arithmetic stays inside the asserted panel bounds:
    // exactly k loads of 64 bytes from each panel and 16 stores of 64
    // bytes to c. w12 is the mandated ZA slice-index register class.
    unsafe {
        core::arch::asm!(
            ".arch armv9.2-a+sme2",
            "smstart",
            "ptrue p0.s",
            "ptrue p1.s",
            "zero {{za}}",
            // k-loop: za0.s += a_col ⊗ b_row.
            "2:",
            "ld1w {{z0.s}}, p0/z, [{a}]",
            "ld1w {{z1.s}}, p1/z, [{b}]",
            "fmopa za0.s, p0/m, p1/m, z0.s, z1.s",
            "add {a}, {a}, #64",
            "add {b}, {b}, #64",
            "subs {k}, {k}, #1",
            "b.ne 2b",
            // Row extraction: horizontal ZA slices to memory.
            "mov w12, #0",
            "3:",
            "st1w {{za0h.s[w12, 0]}}, p0, [{c}]",
            "add {c}, {c}, #64",
            "add w12, w12, #1",
            "cmp w12, #16",
            "b.ne 3b",
            "smstop",
            a = inout(reg) a_panel.as_ptr() => _,
            b = inout(reg) b_panel.as_ptr() => _,
            c = inout(reg) c.as_mut_ptr() => _,
            k = inout(reg) k => _,
            out("w12") _,
            out("v0") _, out("v1") _, out("v2") _, out("v3") _,
            out("v4") _, out("v5") _, out("v6") _, out("v7") _,
            out("v8") _, out("v9") _, out("v10") _, out("v11") _,
            out("v12") _, out("v13") _, out("v14") _, out("v15") _,
            out("v16") _, out("v17") _, out("v18") _, out("v19") _,
            out("v20") _, out("v21") _, out("v22") _, out("v23") _,
            out("v24") _, out("v25") _, out("v26") _, out("v27") _,
            out("v28") _, out("v29") _, out("v30") _, out("v31") _,
            options(nostack)
        );
    }
}

/// The scalar twin: identical panel layout, identical per-element
/// k-order, fused multiply-add, and overwrite semantics — the G0
/// equivalence reference.
pub fn gemm_tile_f32_scalar(a_panel: &[f32], b_panel: &[f32], c: &mut [f32], k: usize) {
    assert_eq!(a_panel.len(), k * TILE);
    assert_eq!(b_panel.len(), k * TILE);
    assert_eq!(c.len(), TILE * TILE);
    c.fill(0.0);
    for kk in 0..k {
        for i in 0..TILE {
            let a = a_panel[kk * TILE + i];
            for j in 0..TILE {
                let b = b_panel[kk * TILE + j];
                c[i * TILE + j] = a.mul_add(b, c[i * TILE + j]);
            }
        }
    }
}
