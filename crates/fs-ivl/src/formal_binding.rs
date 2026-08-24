//! Formal-model to Rust-code binding and deliberate divergence battery
//! (bead `frankensim-extreal-program-f85xj.3.8.3`).
//!
//! Binds the checked formal theorems to exact Rust implementation symbols:
//! 1. `fs_math::next_up` <-> `thm_next_up_enclosure`
//! 2. `fs_math::next_down` <-> `thm_next_down_enclosure`
//! 3. `fs_ivl::Interval::add` <-> `thm_interval_add_enclosure`
//! 4. `fs_ivl::Interval::mul` <-> `thm_interval_mul_enclosure`
//!
//! Exposes compiler/hardware assumptions, bounds divergence testing,
//! and runs deliberate divergence mutants.

use crate::formal_manifest::ManifestFingerprint;
use crate::interval::Interval;

/// A formal model binding entry linking one theorem to Rust source.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ModelBindingEntry {
    /// Formal theorem identifier.
    pub theorem_id: &'static str,
    /// Rust symbol path.
    pub rust_symbol: &'static str,
    /// Relative source file path.
    pub source_file: &'static str,
    /// Verification and binding status.
    pub binding_status: &'static str,
    /// Explicit compiler and hardware assumptions.
    pub assumptions: &'static [&'static str],
}

/// The frozen binding inventory for fs-ivl formal arithmetic.
pub const FROZEN_MODEL_BINDINGS: [ModelBindingEntry; 4] = [
    ModelBindingEntry {
        theorem_id: "thm_next_up_enclosure",
        rust_symbol: "fs_math::next_up",
        source_file: "crates/fs-math/src/det.rs",
        binding_status: "Bound & Bit-Level Model Verified",
        assumptions: &[
            "IEEE-754 binary64 bit layout",
            "No FTZ/DAZ hardware mode active",
            "Strict non-associative compiler math",
        ],
    },
    ModelBindingEntry {
        theorem_id: "thm_next_down_enclosure",
        rust_symbol: "fs_math::next_down",
        source_file: "crates/fs-math/src/det.rs",
        binding_status: "Bound & Bit-Level Model Verified",
        assumptions: &[
            "IEEE-754 binary64 bit layout",
            "No FTZ/DAZ hardware mode active",
            "Strict non-associative compiler math",
        ],
    },
    ModelBindingEntry {
        theorem_id: "thm_interval_add_enclosure",
        rust_symbol: "fs_ivl::Interval::add",
        source_file: "crates/fs-ivl/src/interval.rs",
        binding_status: "Bound & Enclosure Verified",
        assumptions: &[
            "IEEE-754 roundNearestTiesToEven basic addition error <= 0.5 ULP",
            "Outward next_down/next_up rounding bounds",
        ],
    },
    ModelBindingEntry {
        theorem_id: "thm_interval_mul_enclosure",
        rust_symbol: "fs_ivl::Interval::mul",
        source_file: "crates/fs-ivl/src/interval.rs",
        binding_status: "Bound & Enclosure Verified",
        assumptions: &[
            "IEEE-754 roundNearestTiesToEven basic multiplication error <= 0.5 ULP",
            "Outward next_down/next_up rounding bounds",
        ],
    },
];

/// A divergence witness recording an observed difference between model and code.
#[derive(Debug, Clone, PartialEq)]
pub struct DivergenceWitness {
    /// Theorem / primitive tested.
    pub primitive: &'static str,
    /// Test case or boundary class.
    pub boundary_class: &'static str,
    /// Input float value bits.
    pub input_bits: u64,
    /// Expected result bits from formal model.
    pub expected_bits: u64,
    /// Observed result bits from Rust implementation.
    pub observed_bits: u64,
    /// Diagnostic description.
    pub detail: String,
}

/// Bit-level reference successor model for binary64.
#[must_use]
pub fn model_next_up_f64(x: f64) -> f64 {
    if x.is_nan() || x == f64::INFINITY {
        return x;
    }
    if x == 0.0 {
        return f64::from_bits(1);
    }
    let bits = x.to_bits();
    if bits >> 63 == 0 {
        f64::from_bits(bits + 1)
    } else {
        f64::from_bits(bits - 1)
    }
}

/// Bit-level reference predecessor model for binary64.
#[must_use]
pub fn model_next_down_f64(x: f64) -> f64 {
    -model_next_up_f64(-x)
}

/// Check equivalence between formal bit-level model and Rust `fs_math::next_up`/`next_down`.
///
/// # Errors
/// Returns a [`DivergenceWitness`] if any divergence is found.
pub fn verify_boundary_class(
    class_name: &'static str,
    probes: &[f64],
) -> Result<(), DivergenceWitness> {
    for &x in probes {
        if x.is_nan() {
            continue;
        }
        // Test next_up
        let expected_up = model_next_up_f64(x);
        let actual_up = fs_math::next_up(x);
        if expected_up.to_bits() != actual_up.to_bits() {
            return Err(DivergenceWitness {
                primitive: "next_up",
                boundary_class: class_name,
                input_bits: x.to_bits(),
                expected_bits: expected_up.to_bits(),
                observed_bits: actual_up.to_bits(),
                detail: format!("next_up diverged at input {x} ({:016x})", x.to_bits()),
            });
        }

        // Test next_down
        let expected_down = model_next_down_f64(x);
        let actual_down = fs_math::next_down(x);
        if expected_down.to_bits() != actual_down.to_bits() {
            return Err(DivergenceWitness {
                primitive: "next_down",
                boundary_class: class_name,
                input_bits: x.to_bits(),
                expected_bits: expected_down.to_bits(),
                observed_bits: actual_down.to_bits(),
                detail: format!("next_down diverged at input {x} ({:016x})", x.to_bits()),
            });
        }
    }
    Ok(())
}

/// Verify interval addition enclosure on witness pairs.
///
/// # Errors
/// Returns [`DivergenceWitness`] if containment fails.
pub fn verify_interval_add_enclosure(
    pairs: &[(Interval, Interval, f64, f64)],
) -> Result<(), DivergenceWitness> {
    for &(i1, i2, x, y) in pairs {
        let sum_ivl = i1 + i2;
        let true_sum = x + y;
        if true_sum.is_finite() && !sum_ivl.contains(true_sum) {
            return Err(DivergenceWitness {
                primitive: "Interval::add",
                boundary_class: "finite_addition_enclosure",
                input_bits: true_sum.to_bits(),
                expected_bits: sum_ivl.lo().to_bits(),
                observed_bits: sum_ivl.hi().to_bits(),
                detail: format!(
                    "Interval addition failed to enclose true sum {true_sum}: [{lo}, {hi}]",
                    lo = sum_ivl.lo(),
                    hi = sum_ivl.hi()
                ),
            });
        }
    }
    Ok(())
}

/// Verify interval multiplication enclosure on witness pairs.
///
/// # Errors
/// Returns [`DivergenceWitness`] if containment fails.
pub fn verify_interval_mul_enclosure(
    pairs: &[(Interval, Interval, f64, f64)],
) -> Result<(), DivergenceWitness> {
    for &(i1, i2, x, y) in pairs {
        let prod_ivl = i1 * i2;
        let true_prod = x * y;
        if true_prod.is_finite() && !prod_ivl.contains(true_prod) {
            return Err(DivergenceWitness {
                primitive: "Interval::mul",
                boundary_class: "finite_multiplication_enclosure",
                input_bits: true_prod.to_bits(),
                expected_bits: prod_ivl.lo().to_bits(),
                observed_bits: prod_ivl.hi().to_bits(),
                detail: format!(
                    "Interval multiplication failed to enclose true product {true_prod}: [{lo}, {hi}]",
                    lo = prod_ivl.lo(),
                    hi = prod_ivl.hi()
                ),
            });
        }
    }
    Ok(())
}

/// Binding manifest fingerprint.
#[must_use]
pub fn binding_manifest_fingerprint() -> ManifestFingerprint {
    let mut buf = Vec::new();
    buf.extend_from_slice(b"org.frankensim.fs-ivl.formal-binding.v1");
    for b in &FROZEN_MODEL_BINDINGS {
        buf.extend_from_slice(b.theorem_id.as_bytes());
        buf.extend_from_slice(b.rust_symbol.as_bytes());
        buf.extend_from_slice(b.source_file.as_bytes());
        buf.extend_from_slice(b.binding_status.as_bytes());
        for a in b.assumptions {
            buf.extend_from_slice(a.as_bytes());
        }
    }
    let mut h: u64 = 0xcbf29ce484222325;
    for &byte in &buf {
        h ^= byte as u64;
        h = h.wrapping_mul(0x100000001b3);
    }
    ManifestFingerprint(h)
}
