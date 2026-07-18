//! Versioned certificate-plugin registry (bead
//! frankensim-checker-semantic-plugins-9e8n): independently checkable,
//! solver-free witness families the standalone checker can RECOMPUTE —
//! closing the gap where schema-v7 could authenticate a certificate's
//! integrity and origin but not its mathematics.
//!
//! Each plugin owns a closed canonical witness schema (magic + version +
//! bounded fields), explicit resource bounds, a semantic recheck that
//! recomputes the claim from the witness data with outward-rounded
//! native arithmetic, failure localization (the exact op index or matrix
//! row), and a no-claim boundary. Unknown families and versions are
//! explicit [`PluginVerdict::CapabilityRefused`] — never a generic Pass.
//!
//! v1 families:
//! - `interval-enclosure-chain` — a tape of interval operations over
//!   exact f64-bit endpoints; the checker replays the tape with outward
//!   rounding and verifies the claimed final enclosure CONTAINS the
//!   recomputed one.
//! - `linear-residual-linf` — a bounded dense (A, x, b) with a claimed
//!   ∞-norm residual bound; the checker recomputes every row residual as
//!   an outward interval and verifies the whole enclosure sits inside
//!   [−bound, bound].
//!
//! No-claims: a refuted witness refutes THE CERTIFICATE ("this witness
//! proves the claim"), not necessarily the underlying mathematical fact;
//! plugin verdicts supply the INDEPENDENT SEMANTIC VERIFICATION axis and
//! deliberately do not restate integrity ([`crate::IntegrityStatus`]) or
//! origin ([`crate::OriginStatus`]) authority; no solver or geometry
//! dependency enters the crate — recomputation uses native IEEE
//! arithmetic only, so WASM/standalone builds stay green.

use core::fmt;
use std::collections::BTreeMap;

/// Plugin protocol version stamped into every registry decision.
pub const PLUGIN_PROTOCOL_VERSION: u32 = 1;
/// Maximum accepted witness transport bytes.
pub const MAX_WITNESS_BYTES: usize = 256 * 1024;
/// Maximum interval-tape operations.
pub const MAX_CHAIN_OPS: usize = 4096;
/// Maximum interval-tape inputs.
pub const MAX_CHAIN_INPUTS: usize = 1024;
/// Maximum dense matrix dimension per side.
pub const MAX_RESIDUAL_DIM: usize = 64;

const CHAIN_MAGIC: &[u8; 4] = b"FSW1";
const RESIDUAL_MAGIC: &[u8; 4] = b"FSW2";

/// FNV-1a 64 witness content identity (native, dependency-free).
#[must_use]
pub fn witness_identity(bytes: &[u8]) -> u64 {
    let mut h: u64 = 0xcbf2_9ce4_8422_2325;
    for b in bytes {
        h ^= u64::from(*b);
        h = h.wrapping_mul(0x0000_0100_0000_01b3);
    }
    h
}

/// The typed plugin verdict. Exactly one variant grants semantic
/// verification; every other outcome names its axis and location.
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum PluginVerdict {
    /// The witness recomputes and proves its claim.
    SemanticallyVerified {
        /// Recheck detail (what was recomputed).
        detail: String,
    },
    /// Structurally valid witness whose mathematics do NOT check.
    SemanticallyRefuted {
        /// The exact failing location (op index / matrix row).
        location: String,
        /// What failed.
        detail: String,
    },
    /// The witness bytes violate the closed schema or its bounds.
    WitnessMalformed {
        /// Byte offset of the refusal.
        offset: usize,
        /// What was malformed.
        detail: String,
    },
    /// No registered plugin owns (family, version): an explicit
    /// capability refusal, never a generic Pass.
    CapabilityRefused {
        /// The requested family.
        family: String,
        /// The requested version.
        version: u32,
        /// Why nothing can check it.
        detail: String,
    },
}

impl PluginVerdict {
    /// True only for [`PluginVerdict::SemanticallyVerified`].
    #[must_use]
    pub const fn is_verified(&self) -> bool {
        matches!(self, Self::SemanticallyVerified { .. })
    }
}

impl fmt::Display for PluginVerdict {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::SemanticallyVerified { detail } => write!(f, "verified: {detail}"),
            Self::SemanticallyRefuted { location, detail } => {
                write!(f, "REFUTED at {location}: {detail}")
            }
            Self::WitnessMalformed { offset, detail } => {
                write!(f, "malformed at byte {offset}: {detail}")
            }
            Self::CapabilityRefused {
                family,
                version,
                detail,
            } => write!(f, "capability refused for {family}@v{version}: {detail}"),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum BuiltinChecker {
    IntervalEnclosureChain,
    LinearResidualLinf,
}

/// The versioned plugin registry: a CLOSED set per protocol version.
#[derive(Debug, Clone)]
pub struct PluginRegistry {
    entries: BTreeMap<(String, u32), BuiltinChecker>,
}

impl PluginRegistry {
    /// The v1 registry: exactly the two solver-free launch families.
    #[must_use]
    pub fn v1() -> PluginRegistry {
        let mut entries = BTreeMap::new();
        entries.insert(
            ("interval-enclosure-chain".to_owned(), 1),
            BuiltinChecker::IntervalEnclosureChain,
        );
        entries.insert(
            ("linear-residual-linf".to_owned(), 1),
            BuiltinChecker::LinearResidualLinf,
        );
        PluginRegistry { entries }
    }

    /// The registered (family, version) keys, deterministic order.
    #[must_use]
    pub fn families(&self) -> Vec<(String, u32)> {
        self.entries.keys().cloned().collect()
    }

    /// Check one witness under (family, version). Unknown keys refuse
    /// explicitly; nothing outside the closed registry can Pass.
    #[must_use]
    pub fn check(&self, family: &str, version: u32, witness: &[u8]) -> PluginVerdict {
        if witness.len() > MAX_WITNESS_BYTES {
            return PluginVerdict::WitnessMalformed {
                offset: 0,
                detail: format!(
                    "witness {} bytes above the {MAX_WITNESS_BYTES}-byte bound",
                    witness.len()
                ),
            };
        }
        match self.entries.get(&(family.to_owned(), version)) {
            None => PluginVerdict::CapabilityRefused {
                family: family.to_owned(),
                version,
                detail: "no registered plugin owns this family/version; independent \
                         semantic verification is unavailable"
                    .to_owned(),
            },
            Some(BuiltinChecker::IntervalEnclosureChain) => check_chain(witness),
            Some(BuiltinChecker::LinearResidualLinf) => check_residual(witness),
        }
    }
}

// ---------------------------------------------------------------- reader

struct Reader<'a> {
    bytes: &'a [u8],
    pos: usize,
}

impl<'a> Reader<'a> {
    fn take(&mut self, n: usize) -> Result<&'a [u8], PluginVerdict> {
        let end = self
            .pos
            .checked_add(n)
            .filter(|&e| e <= self.bytes.len())
            .ok_or(PluginVerdict::WitnessMalformed {
                offset: self.pos,
                detail: "truncated witness".to_owned(),
            })?;
        let out = &self.bytes[self.pos..end];
        self.pos = end;
        Ok(out)
    }

    fn u16(&mut self) -> Result<u16, PluginVerdict> {
        let b = self.take(2)?;
        Ok(u16::from_be_bytes([b[0], b[1]]))
    }

    fn u32(&mut self) -> Result<u32, PluginVerdict> {
        let b = self.take(4)?;
        Ok(u32::from_be_bytes([b[0], b[1], b[2], b[3]]))
    }

    fn f64_bits(&mut self) -> Result<f64, PluginVerdict> {
        let b = self.take(8)?;
        let mut arr = [0u8; 8];
        arr.copy_from_slice(b);
        Ok(f64::from_bits(u64::from_be_bytes(arr)))
    }

    fn finite(&mut self, what: &str) -> Result<f64, PluginVerdict> {
        let offset = self.pos;
        let v = self.f64_bits()?;
        if !v.is_finite() {
            return Err(PluginVerdict::WitnessMalformed {
                offset,
                detail: format!("{what} is not finite"),
            });
        }
        Ok(v)
    }

    fn done(&self) -> Result<(), PluginVerdict> {
        if self.pos == self.bytes.len() {
            Ok(())
        } else {
            Err(PluginVerdict::WitnessMalformed {
                offset: self.pos,
                detail: "trailing bytes after the closed schema".to_owned(),
            })
        }
    }
}

fn interval(lo: f64, hi: f64, offset: usize, what: &str) -> Result<(f64, f64), PluginVerdict> {
    if lo > hi {
        return Err(PluginVerdict::WitnessMalformed {
            offset,
            detail: format!("{what} interval inverted"),
        });
    }
    Ok((lo, hi))
}

// -------------------------------------------- family 1: interval chain

fn outward_add(a: (f64, f64), b: (f64, f64)) -> (f64, f64) {
    ((a.0 + b.0).next_down(), (a.1 + b.1).next_up())
}

fn outward_mul(a: (f64, f64), b: (f64, f64)) -> (f64, f64) {
    let products = [a.0 * b.0, a.0 * b.1, a.1 * b.0, a.1 * b.1];
    let lo = products.iter().copied().fold(f64::INFINITY, f64::min);
    let hi = products.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    (lo.next_down(), hi.next_up())
}

fn check_chain(witness: &[u8]) -> PluginVerdict {
    match check_chain_inner(witness) {
        Ok(verdict) | Err(verdict) => verdict,
    }
}

fn check_chain_inner(witness: &[u8]) -> Result<PluginVerdict, PluginVerdict> {
    let mut r = Reader {
        bytes: witness,
        pos: 0,
    };
    if r.take(4)? != CHAIN_MAGIC {
        return Ok(PluginVerdict::WitnessMalformed {
            offset: 0,
            detail: "bad interval-chain magic".to_owned(),
        });
    }
    let proto = r.u32()?;
    if proto != PLUGIN_PROTOCOL_VERSION {
        return Ok(PluginVerdict::WitnessMalformed {
            offset: 4,
            detail: format!("unknown chain protocol {proto}"),
        });
    }
    let n_inputs = r.u16()? as usize;
    if n_inputs == 0 || n_inputs > MAX_CHAIN_INPUTS {
        return Ok(PluginVerdict::WitnessMalformed {
            offset: r.pos,
            detail: format!("input count {n_inputs} outside 1..={MAX_CHAIN_INPUTS}"),
        });
    }
    let mut tape: Vec<(f64, f64)> = Vec::with_capacity(n_inputs);
    for i in 0..n_inputs {
        let offset = r.pos;
        let lo = r.finite(&format!("input {i} lo"))?;
        let hi = r.finite(&format!("input {i} hi"))?;
        tape.push(interval(lo, hi, offset, &format!("input {i}"))?);
    }
    let n_ops = r.u16()? as usize;
    if n_ops == 0 || n_ops > MAX_CHAIN_OPS {
        return Ok(PluginVerdict::WitnessMalformed {
            offset: r.pos,
            detail: format!("op count {n_ops} outside 1..={MAX_CHAIN_OPS}"),
        });
    }
    for op_index in 0..n_ops {
        let offset = r.pos;
        let opcode = r.take(1)?[0];
        let a = r.u16()? as usize;
        let b = r.u16()? as usize;
        let limit = tape.len();
        let fetch = |k: usize| -> Result<(f64, f64), PluginVerdict> {
            tape.get(k).copied().ok_or(PluginVerdict::WitnessMalformed {
                offset,
                detail: format!("op {op_index} operand {k} out of range (tape {limit})"),
            })
        };
        let value = match opcode {
            1 => outward_add(fetch(a)?, fetch(b)?),
            2 => outward_mul(fetch(a)?, fetch(b)?),
            3 => {
                let x = fetch(a)?;
                (-x.1, -x.0)
            }
            other => {
                return Ok(PluginVerdict::WitnessMalformed {
                    offset,
                    detail: format!("op {op_index} unknown opcode {other}"),
                });
            }
        };
        if !value.0.is_finite() || !value.1.is_finite() {
            return Ok(PluginVerdict::SemanticallyRefuted {
                location: format!("op {op_index}"),
                detail: "recomputed enclosure escaped the finite range".to_owned(),
            });
        }
        tape.push(value);
    }
    let claim_offset = r.pos;
    let claimed_lo = r.finite("claimed lo")?;
    let claimed_hi = r.finite("claimed hi")?;
    let claimed = interval(claimed_lo, claimed_hi, claim_offset, "claimed enclosure")?;
    r.done()?;
    let recomputed = *tape.last().expect("tape has at least the inputs");
    if claimed.0 <= recomputed.0 && recomputed.1 <= claimed.1 {
        Ok(PluginVerdict::SemanticallyVerified {
            detail: format!(
                "replayed {n_ops} interval ops over {n_inputs} inputs; recomputed \
                 [{:e}, {:e}] within claimed [{claimed_lo:e}, {claimed_hi:e}]",
                recomputed.0, recomputed.1
            ),
        })
    } else {
        Ok(PluginVerdict::SemanticallyRefuted {
            location: format!("final op {}", n_ops - 1),
            detail: format!(
                "claimed enclosure [{claimed_lo:e}, {claimed_hi:e}] does not contain the \
                 recomputed [{:e}, {:e}]",
                recomputed.0, recomputed.1
            ),
        })
    }
}

// ---------------------------------------- family 2: linear residual ∞

fn check_residual(witness: &[u8]) -> PluginVerdict {
    match check_residual_inner(witness) {
        Ok(verdict) | Err(verdict) => verdict,
    }
}

fn check_residual_inner(witness: &[u8]) -> Result<PluginVerdict, PluginVerdict> {
    let mut r = Reader {
        bytes: witness,
        pos: 0,
    };
    if r.take(4)? != RESIDUAL_MAGIC {
        return Ok(PluginVerdict::WitnessMalformed {
            offset: 0,
            detail: "bad linear-residual magic".to_owned(),
        });
    }
    let proto = r.u32()?;
    if proto != PLUGIN_PROTOCOL_VERSION {
        return Ok(PluginVerdict::WitnessMalformed {
            offset: 4,
            detail: format!("unknown residual protocol {proto}"),
        });
    }
    let rows = r.u16()? as usize;
    let cols = r.u16()? as usize;
    if rows == 0 || cols == 0 || rows > MAX_RESIDUAL_DIM || cols > MAX_RESIDUAL_DIM {
        return Ok(PluginVerdict::WitnessMalformed {
            offset: r.pos,
            detail: format!("dimensions {rows}x{cols} outside 1..={MAX_RESIDUAL_DIM} per side"),
        });
    }
    let mut a = Vec::with_capacity(rows * cols);
    for k in 0..rows * cols {
        a.push(r.finite(&format!("A[{k}]"))?);
    }
    let mut x = Vec::with_capacity(cols);
    for k in 0..cols {
        x.push(r.finite(&format!("x[{k}]"))?);
    }
    let mut b = Vec::with_capacity(rows);
    for k in 0..rows {
        b.push(r.finite(&format!("b[{k}]"))?);
    }
    let bound_offset = r.pos;
    let bound = r.finite("claimed bound")?;
    if bound < 0.0 {
        return Ok(PluginVerdict::WitnessMalformed {
            offset: bound_offset,
            detail: "claimed residual bound is negative".to_owned(),
        });
    }
    r.done()?;

    for i in 0..rows {
        // Outward enclosure of r_i = b_i − Σ_j A_ij·x_j.
        let mut acc = (b[i], b[i]);
        for j in 0..cols {
            let product = a[i * cols + j] * x[j];
            let outward = (product.next_down(), product.next_up());
            acc = outward_add(acc, (-outward.1, -outward.0));
        }
        if !(acc.0 >= -bound && acc.1 <= bound) {
            return Ok(PluginVerdict::SemanticallyRefuted {
                location: format!("row {i}"),
                detail: format!(
                    "residual enclosure [{:e}, {:e}] is not provably inside \
                     [-{bound:e}, {bound:e}] — this witness cannot certify the bound",
                    acc.0, acc.1
                ),
            });
        }
    }
    Ok(PluginVerdict::SemanticallyVerified {
        detail: format!(
            "recomputed all {rows} row residuals in outward interval arithmetic; \
             every enclosure inside [-{bound:e}, {bound:e}]"
        ),
    })
}

// ---------------------------------------------------------- producers

/// Producer helper: encode an interval-chain witness (the exact closed
/// schema the checker replays). Producers and packages share this
/// canonical form; the standalone checker only ever reads it.
#[must_use]
pub fn encode_chain_witness(
    inputs: &[(f64, f64)],
    ops: &[(u8, u16, u16)],
    claimed: (f64, f64),
) -> Vec<u8> {
    let mut out = Vec::new();
    out.extend_from_slice(CHAIN_MAGIC);
    out.extend_from_slice(&PLUGIN_PROTOCOL_VERSION.to_be_bytes());
    out.extend_from_slice(&(inputs.len() as u16).to_be_bytes());
    for (lo, hi) in inputs {
        out.extend_from_slice(&lo.to_bits().to_be_bytes());
        out.extend_from_slice(&hi.to_bits().to_be_bytes());
    }
    out.extend_from_slice(&(ops.len() as u16).to_be_bytes());
    for (opcode, a, b) in ops {
        out.push(*opcode);
        out.extend_from_slice(&a.to_be_bytes());
        out.extend_from_slice(&b.to_be_bytes());
    }
    out.extend_from_slice(&claimed.0.to_bits().to_be_bytes());
    out.extend_from_slice(&claimed.1.to_bits().to_be_bytes());
    out
}

/// Producer helper: encode a linear-residual witness.
#[must_use]
pub fn encode_residual_witness(
    rows: usize,
    cols: usize,
    a: &[f64],
    x: &[f64],
    b: &[f64],
    claimed_bound: f64,
) -> Vec<u8> {
    let mut out = Vec::new();
    out.extend_from_slice(RESIDUAL_MAGIC);
    out.extend_from_slice(&PLUGIN_PROTOCOL_VERSION.to_be_bytes());
    out.extend_from_slice(&(rows as u16).to_be_bytes());
    out.extend_from_slice(&(cols as u16).to_be_bytes());
    for v in a.iter().chain(x).chain(b) {
        out.extend_from_slice(&v.to_bits().to_be_bytes());
    }
    out.extend_from_slice(&claimed_bound.to_bits().to_be_bytes());
    out
}
