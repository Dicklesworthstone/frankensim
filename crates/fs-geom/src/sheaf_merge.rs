//! SHEAF-ADJUDICATED THREE-WAY MERGE (addendum Proposal 10, bead
//! lmp4.12 — THE CROWN JEWEL; [M], behind the `sheaf-merge` feature
//! until its Gauntlet tier + kill metric are green): the sheaf
//! machinery built for watertightness is, UNMODIFIED, a merge-conflict
//! classifier. Auto-resolution happens EXACTLY when it is
//! mathematically licensed — no more, no less:
//!
//! - the union mismatch's COBOUNDARY component admits a canonical
//!   least-squares reconciliation → applied automatically, with the
//!   reconciled state's own certificate RE-VERIFIED before the merge is
//!   reported resolved (Sev-0: a passing certificate is never attached
//!   over a state that then fails watertightness);
//! - the HARMONIC component is a genuine obstruction class — no local
//!   adjustment can resolve it → a STRUCTURAL-CONFLICT object localized
//!   to its supporting interface cells, carrying both parents'
//!   provenance;
//! - NON-GEOMETRIC collisions (load cases, material assignments) are
//!   TYPE-LEVEL conflicts caught before any decomposition runs;
//! - trust is CONDITIONED on the complex's spectral gap (Proposal 5,
//!   risk R5): merges in degraded-gap regions are flagged
//!   low-confidence.

use crate::sheaf_repair::{SheafSkeleton, apply_gauge, hodge_decompose};
use std::collections::BTreeMap;

/// Energy fraction below which a Hodge component is treated as absent.
pub const COMPONENT_FLOOR: f64 = 1e-9;

/// One branch's edits: a mismatch cochain plus non-geometric keyed
/// assignments (load cases, materials — the typed layer's inputs).
#[derive(Debug, Clone)]
pub struct BranchState {
    /// Provenance label (commit root, agent id…).
    pub provenance: String,
    /// The branch's interface mismatch cochain.
    pub mismatch: Vec<f64>,
    /// Non-geometric keyed assignments.
    pub assignments: BTreeMap<String, String>,
}

/// A structural (harmonic) conflict: no local fix exists.
#[derive(Debug, Clone, PartialEq)]
pub struct StructuralConflict {
    /// Supporting interface cells (patch pairs) with magnitudes,
    /// strongest first.
    pub cells: Vec<((usize, usize), f64)>,
    /// Parent provenances (X, Y).
    pub parents: (String, String),
}

/// A type-level conflict: both branches edited the same non-geometric
/// key to different values.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TypeConflict {
    /// The colliding key.
    pub key: String,
    /// X's value.
    pub x_value: String,
    /// Y's value.
    pub y_value: String,
}

/// Merge trust, conditioned on the spectral gap (Proposal 5 / R5).
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Confidence {
    /// The complex is well-coupled at the given gap.
    Normal {
        /// The algebraic-connectivity gap.
        gap: f64,
    },
    /// Degraded gap: the harmonic/coboundary split is less trustworthy
    /// here — treat the merge as provisional.
    LowGap {
        /// The measured gap.
        gap: f64,
        /// The threshold it fell below.
        threshold: f64,
    },
}

/// The reconciliation certificate attached to auto-resolved merges —
/// RE-VERIFIED against the reconciled state, never assumed.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MergeCertificate {
    /// Worst post-reconciliation interface mismatch.
    pub post_norm: f64,
    /// The tolerance it passed.
    pub tol: f64,
}

/// The merge verdict.
#[derive(Debug, Clone, PartialEq)]
pub enum MergeOutcome {
    /// A boundary fast path fired (X == Y, or one branch unchanged):
    /// no decomposition needed.
    Trivial {
        /// Which fast path.
        reason: &'static str,
        /// The merged cochain.
        merged: Vec<f64>,
    },
    /// The mismatch was (numerically) pure coboundary and the canonical
    /// least-squares reconciliation REACHED a certificate-passing state.
    Resolved {
        /// The reconciled cochain.
        merged: Vec<f64>,
        /// The gauge offsets applied.
        gauge: Vec<f64>,
        /// The re-verified certificate.
        certificate: MergeCertificate,
        /// Gap-conditioned trust.
        confidence: Confidence,
    },
    /// Genuine conflicts: structural (harmonic) and/or type-level.
    Conflicted {
        /// Harmonic obstructions (empty if only type-level).
        structural: Vec<StructuralConflict>,
        /// Keyed-assignment collisions.
        type_conflicts: Vec<TypeConflict>,
        /// Gap-conditioned trust.
        confidence: Confidence,
    },
    /// The Sev-0 guard: reconciliation could NOT reach a
    /// certificate-passing state (e.g. a coexact residue) — escalated
    /// unresolved rather than certified falsely.
    EscalatedUnresolved {
        /// The post-reconciliation norm that failed.
        post_norm: f64,
        /// The tolerance it failed against.
        tol: f64,
        /// What is left (exact, coexact, harmonic fractions).
        fractions: (f64, f64, f64),
    },
}

fn norm_inf(v: &[f64]) -> f64 {
    v.iter().fold(0.0f64, |a, &b| a.max(b.abs()))
}

/// Dense symmetric eigenvalues by cyclic Jacobi (small matrices —
/// patch counts, not DOF counts).
#[allow(clippy::needless_range_loop)] // rotations touch (k,p),(k,q) pairs
fn jacobi_eigenvalues(mut a: Vec<Vec<f64>>) -> Vec<f64> {
    let n = a.len();
    for _sweep in 0..64 {
        let mut off = 0.0f64;
        for p in 0..n {
            for q in (p + 1)..n {
                off += a[p][q] * a[p][q];
            }
        }
        if off < 1e-24 {
            break;
        }
        for p in 0..n {
            for q in (p + 1)..n {
                if a[p][q].abs() < 1e-300 {
                    continue;
                }
                let theta = (a[q][q] - a[p][p]) / (2.0 * a[p][q]);
                let t = theta.signum() / (theta.abs() + (theta * theta + 1.0).sqrt());
                let c = 1.0 / (t * t + 1.0).sqrt();
                let s = t * c;
                for k in 0..n {
                    let (akp, akq) = (a[k][p], a[k][q]);
                    a[k][p] = c * akp - s * akq;
                    a[k][q] = s * akp + c * akq;
                }
                for k in 0..n {
                    let (apk, aqk) = (a[p][k], a[q][k]);
                    a[p][k] = c * apk - s * aqk;
                    a[q][k] = s * apk + c * aqk;
                }
            }
        }
    }
    let mut ev: Vec<f64> = (0..n).map(|i| a[i][i]).collect();
    ev.sort_by(f64::total_cmp);
    ev
}

/// The spectral gap (algebraic connectivity λ₂) of the weighted
/// patch-adjacency Laplacian — the Proposal-5 trust signal. Weights
/// default to 1 (e.g. interface sample counts belong here).
#[must_use]
pub fn spectral_gap(skeleton: &SheafSkeleton, weights: Option<&[f64]>) -> f64 {
    let n = skeleton.n_patches;
    let mut lap = vec![vec![0.0f64; n]; n];
    for (k, &(u, v)) in skeleton.edges.iter().enumerate() {
        let w = weights.map_or(1.0, |ws| ws[k]);
        lap[u][u] += w;
        lap[v][v] += w;
        lap[u][v] -= w;
        lap[v][u] -= w;
    }
    let ev = jacobi_eigenvalues(lap);
    // λ₁ ≈ 0 (kernel = components); the gap is the next eigenvalue.
    ev.iter().copied().find(|&e| e > 1e-9).unwrap_or(0.0)
}

/// Detect keyed-assignment collisions (the typed layer's cheap half —
/// coupling-graph legality of the merged assignment set is fs-iface's
/// contract at its own layer).
#[must_use]
pub fn type_conflicts(x: &BranchState, y: &BranchState) -> Vec<TypeConflict> {
    let mut out = Vec::new();
    for (k, xv) in &x.assignments {
        if let Some(yv) = y.assignments.get(k)
            && xv != yv
        {
            out.push(TypeConflict {
                key: k.clone(),
                x_value: xv.clone(),
                y_value: yv.clone(),
            });
        }
    }
    out
}

/// The three-way merge. `base` is the common ancestor's mismatch
/// cochain; `tol` is the watertightness tolerance the reconciled state
/// must PASS to be reported resolved; `gap_threshold` conditions trust.
#[must_use]
#[allow(clippy::too_many_lines)]
pub fn three_way_merge(
    skeleton: &SheafSkeleton,
    base: &[f64],
    x: &BranchState,
    y: &BranchState,
    weights: Option<&[f64]>,
    tol: f64,
    gap_threshold: f64,
) -> MergeOutcome {
    assert_eq!(base.len(), skeleton.edges.len(), "base cochain length");
    assert_eq!(x.mismatch.len(), base.len(), "X cochain length");
    assert_eq!(y.mismatch.len(), base.len(), "Y cochain length");
    // Type-level conflicts are caught BEFORE any decomposition.
    let tc = type_conflicts(x, y);
    let gap = spectral_gap(skeleton, weights);
    let confidence = if gap < gap_threshold {
        Confidence::LowGap {
            gap,
            threshold: gap_threshold,
        }
    } else {
        Confidence::Normal { gap }
    };
    if !tc.is_empty() {
        return MergeOutcome::Conflicted {
            structural: Vec::new(),
            type_conflicts: tc,
            confidence,
        };
    }
    // Boundary fast paths: no decomposition, no false ceremony.
    let bits = |v: &[f64]| -> Vec<u64> { v.iter().map(|f| f.to_bits()).collect() };
    if bits(&x.mismatch) == bits(&y.mismatch) {
        return MergeOutcome::Trivial {
            reason: "branches identical",
            merged: x.mismatch.clone(),
        };
    }
    if bits(&x.mismatch) == bits(base) {
        return MergeOutcome::Trivial {
            reason: "X unchanged from base",
            merged: y.mismatch.clone(),
        };
    }
    if bits(&y.mismatch) == bits(base) {
        return MergeOutcome::Trivial {
            reason: "Y unchanged from base",
            merged: x.mismatch.clone(),
        };
    }
    // The naive union of edits at the cochain level: X + Y − B.
    let union: Vec<f64> = x
        .mismatch
        .iter()
        .zip(&y.mismatch)
        .zip(base)
        .map(|((a, b), c)| a + b - c)
        .collect();
    let split = hodge_decompose(skeleton, &union);
    // Coboundary reconciliation FIRST: the canonical least-squares
    // gauge. Auto-resolution is licensed exactly when the reconciled
    // state passes its own certificate — a harmonic remnant BELOW the
    // watertightness tolerance is not an obstruction.
    let merged = apply_gauge(skeleton, &union, &split.potential);
    // Sev-0 RE-VERIFICATION: the reconciled state's own certificate,
    // checked — never attached on faith.
    let post_norm = norm_inf(&merged);
    if post_norm <= tol {
        return MergeOutcome::Resolved {
            merged,
            gauge: split.potential,
            certificate: MergeCertificate { post_norm, tol },
            confidence,
        };
    }
    // Verification failed: classify the failure. A dominant harmonic
    // residue is a genuine obstruction class (no local fix exists);
    // anything else (coexact circulation) escalates unresolved.
    let harmonic_norm = norm_inf(&split.harmonic);
    let coexact_norm = norm_inf(&split.coexact);
    if harmonic_norm > tol && harmonic_norm >= coexact_norm {
        let mut cells: Vec<((usize, usize), f64)> = skeleton
            .edges
            .iter()
            .zip(&split.harmonic)
            .filter(|(_, h)| h.abs() > tol)
            .map(|(&e, &h)| (e, h.abs()))
            .collect();
        cells.sort_by(|a, b| b.1.total_cmp(&a.1));
        return MergeOutcome::Conflicted {
            structural: vec![StructuralConflict {
                cells,
                parents: (x.provenance.clone(), y.provenance.clone()),
            }],
            type_conflicts: Vec::new(),
            confidence,
        };
    }
    MergeOutcome::EscalatedUnresolved {
        post_norm,
        tol,
        fractions: split.fractions,
    }
}

/// The KILL-CRITERION harness (Proposal 10): run seeded random
/// three-way merges and measure the harmonic-conflict rate. If more
/// than ~25% of realistic merges are structural, agents are colliding
/// and merge-based concurrency is the wrong model.
#[must_use]
pub fn harmonic_conflict_rate(
    skeleton: &SheafSkeleton,
    trials: usize,
    edit_scale: f64,
    seed: u64,
) -> f64 {
    let mut state = seed;
    let mut lcg = move || {
        state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((state >> 11) as f64) / (1u64 << 53) as f64 - 0.5
    };
    let n_edges = skeleton.edges.len();
    let n_patches = skeleton.n_patches;
    let mut conflicts = 0usize;
    for _ in 0..trials {
        // Realistic edits: each branch re-gauges a random subset of
        // patches (coboundary-style work) — collisions arise when the
        // edits interact around cycles.
        let mut edit = |scale: f64| -> Vec<f64> {
            let offsets: Vec<f64> = (0..n_patches).map(|_| scale * lcg()).collect();
            let mut m = skeleton.d0(&offsets);
            // A small amount of independent interface noise.
            for v in &mut m {
                *v += 0.01 * scale * lcg();
            }
            m
        };
        let base = vec![0.0f64; n_edges];
        let x = BranchState {
            provenance: "trial-x".to_string(),
            mismatch: edit(edit_scale),
            assignments: BTreeMap::new(),
        };
        let y = BranchState {
            provenance: "trial-y".to_string(),
            mismatch: edit(edit_scale),
            assignments: BTreeMap::new(),
        };
        let out = three_way_merge(skeleton, &base, &x, &y, None, 0.05 * edit_scale, 1e-6);
        if matches!(out, MergeOutcome::Conflicted { ref structural, .. } if !structural.is_empty())
        {
            conflicts += 1;
        }
    }
    #[allow(clippy::cast_precision_loss)]
    {
        conflicts as f64 / trials as f64
    }
}
