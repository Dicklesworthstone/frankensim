//! SHEAF CERTIFICATES (plan §7.3, Bet 11): the cellular-sheaf
//! watertightness certificate for multi-representation models. A model
//! whose patches live in different charts is a cellular sheaf over the
//! patch-adjacency complex: stalks are per-patch sample spaces,
//! restriction maps are selection operators onto shared interface
//! samples, and GLOBAL CONSISTENCY is the existence of a section. The
//! watertightness certificate is an INTERVAL-VERIFIED bound
//! `‖δs‖∞ ≤ tol` on the interface mismatch cocycle; when the coboundary
//! cannot be driven below tolerance, the H¹ obstruction is reported WITH
//! THE OFFENDING INTERFACE CELLS ATTACHED — exactly the diagnostic an
//! agent needs to fix a leaky model.
//!
//! The construction is finite linear algebra: δ⁰ and δ¹ assemble as
//! sparse matrices with entries in {−1, 0, +1} (restrictions are point
//! samplers), so `δ¹·δ⁰ = 0` holds BITWISE — small-integer f64
//! arithmetic is exact. The least-squares section solve (per-patch gauge
//! offsets over the adjacency Laplacian) splits the mismatch into a
//! reconcilable coboundary part and the structural residual — the same
//! split Proposal 10's merge semantics reuses unmodified.

use crate::{Aabb, Chart, ChartSample, Point3};
use fs_evidence::{
    Evidence, ModelEvidence, NumericalCertificate, NumericalKind, ProvenanceHash,
    SensitivitySummary, StatisticalCertificate,
};
use fs_exec::Cx;
use fs_ivl::Interval;
use fs_sparse::{Coo, Csr};
use std::fmt::Write as _;

/// Samples drawn per pairwise interface.
pub const SAMPLES_PER_INTERFACE: usize = 32;

/// Zero-band half-width as a fraction of the overlap-box diagonal:
/// a point belongs to the shared surface region when BOTH charts place
/// it within this band of their zero set.
pub const BAND_FRACTION: f64 = 0.05;

/// One shared interface sample.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct InterfaceSample {
    /// The world-space point.
    pub point: Point3,
    /// Enclosures of the two charts' signed distances at the point
    /// (outward bounds from each chart's own error certificate).
    pub values: [Interval; 2],
}

/// One interface (edge of the adjacency complex).
#[derive(Debug, Clone)]
pub struct Interface {
    /// Patch indices (u < v; the edge is oriented u → v).
    pub patches: (usize, usize),
    /// Shared samples.
    pub samples: Vec<InterfaceSample>,
}

/// A triple junction (2-cell): three patches with a common overlap.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TripleCell {
    /// Patch indices (sorted).
    pub patches: (usize, usize, usize),
    /// Sample count shared by all three pairwise interfaces.
    pub samples: usize,
}

/// The patch-adjacency sheaf complex.
#[derive(Debug)]
pub struct SheafComplex {
    /// Patch count.
    pub n_patches: usize,
    /// Pairwise interfaces (sorted by patch pair — deterministic).
    pub interfaces: Vec<Interface>,
    /// Triple junctions.
    pub triples: Vec<TripleCell>,
}

/// One interface's assessed mismatch (verdict bits are predicate-sound;
/// reported magnitudes are midpoint±width reconstructions).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct InterfaceBound {
    /// The patch pair.
    pub patches: (usize, usize),
    /// Every sample's |mismatch| enclosure lies inside `[0, tol]`.
    pub all_within: bool,
    /// Some sample's enclosure lies ENTIRELY above tol (proven leak).
    pub proven_leak: bool,
    /// Reported worst lower bound.
    pub lo_report: f64,
    /// Reported worst upper bound.
    pub hi_report: f64,
}

/// The certificate verdict.
#[derive(Debug, Clone, PartialEq)]
pub enum SheafVerdict {
    /// `‖δs‖∞ ≤ tol` with the enclosure upper bound as the margin —
    /// per-interface margins attached.
    Pass {
        /// Upper bound of the worst interface mismatch.
        worst_mismatch: f64,
        /// Per-interface (patch pair, mismatch upper bound).
        margins: Vec<((usize, usize), f64)>,
    },
    /// The H¹ obstruction: some interface's mismatch enclosure lies
    /// ENTIRELY above tolerance — a proven leak, localized.
    Fail {
        /// Offending interfaces: (patch pair, mismatch lower bound).
        obstruction: Vec<((usize, usize), f64)>,
        /// The reconcilable (coboundary) share of the raw mismatch in
        /// [0, 1]: how much a re-gauge would fix (the merge-semantics
        /// split — near 1 means gauge drift, near 0 means structure).
        coboundary_share: f64,
    },
    /// Enclosures straddle the tolerance: no sound claim either way
    /// (tighten chart certificates or the band and re-run).
    Unknown {
        /// Straddling interfaces: (patch pair, lower, upper).
        straddling: Vec<((usize, usize), f64, f64)>,
    },
}

fn fnv(bytes: &[u8]) -> u64 {
    let mut h = 0xcbf2_9ce4_8422_2325u64;
    for &b in bytes {
        h ^= u64::from(b);
        h = h.wrapping_mul(0x0000_0100_0000_01B3);
    }
    h
}

/// Geometry-derived deterministic seed: identical boxes give identical
/// samples regardless of patch indexing (re-index invariance is exact).
fn box_seed(b: &Aabb) -> u64 {
    let mut bytes = Vec::with_capacity(48);
    for v in [b.min.x, b.min.y, b.min.z, b.max.x, b.max.y, b.max.z] {
        bytes.extend_from_slice(&v.to_bits().to_le_bytes());
    }
    fnv(&bytes)
}

fn lcg(state: &mut u64) -> f64 {
    *state = state
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    ((*state >> 11) as f64) / (1u64 << 53) as f64
}

/// The overlap box of two supports, if any.
fn overlap(a: &Aabb, b: &Aabb) -> Option<Aabb> {
    let min = Point3::new(
        a.min.x.max(b.min.x),
        a.min.y.max(b.min.y),
        a.min.z.max(b.min.z),
    );
    let max = Point3::new(
        a.max.x.min(b.max.x),
        a.max.y.min(b.max.y),
        a.max.z.min(b.max.z),
    );
    (min.x < max.x && min.y < max.y && min.z < max.z).then(|| Aabb::new(min, max))
}

/// Outward enclosure of a chart sample's signed distance from its own
/// error certificate (honest: NoClaim charts poison the interface into
/// an infinite interval, which can only yield Unknown).
fn sample_interval(s: &ChartSample) -> Interval {
    match s.error.kind {
        NumericalKind::Exact | NumericalKind::Enclosure => Interval::new(s.error.lo, s.error.hi),
        NumericalKind::Estimate => {
            // Estimates are widened by their own band again (documented
            // conservatism: estimates carry no rigor).
            let w = (s.error.hi - s.error.lo).abs();
            Interval::new(s.error.lo - w, s.error.hi + w)
        }
        NumericalKind::NoClaim => Interval::new(f64::NEG_INFINITY, f64::INFINITY),
    }
}

impl SheafComplex {
    /// Build the complex from charts: interface discovery via
    /// support-overlap, shared-surface samples via the zero band of BOTH
    /// charts, triple junctions via triple overlaps.
    pub fn from_charts(charts: &[&dyn Chart], cx: &Cx<'_>) -> SheafComplex {
        let n = charts.len();
        let mut interfaces = Vec::new();
        for u in 0..n {
            for v in (u + 1)..n {
                let Some(shared) = overlap(&charts[u].support(), &charts[v].support()) else {
                    continue;
                };
                let diag = shared.max.delta_from(shared.min).norm();
                let band = BAND_FRACTION * diag;
                let mut state = box_seed(&shared);
                let mut samples = Vec::new();
                // Rejection-sample the shared zero band (bounded draws).
                for _ in 0..SAMPLES_PER_INTERFACE * 64 {
                    if samples.len() >= SAMPLES_PER_INTERFACE {
                        break;
                    }
                    let p = Point3::new(
                        shared.min.x + lcg(&mut state) * (shared.max.x - shared.min.x),
                        shared.min.y + lcg(&mut state) * (shared.max.y - shared.min.y),
                        shared.min.z + lcg(&mut state) * (shared.max.z - shared.min.z),
                    );
                    let su = charts[u].eval(p, cx);
                    let sv = charts[v].eval(p, cx);
                    if su.signed_distance.abs() <= band && sv.signed_distance.abs() <= band {
                        samples.push(InterfaceSample {
                            point: p,
                            values: [sample_interval(&su), sample_interval(&sv)],
                        });
                    }
                }
                if !samples.is_empty() {
                    interfaces.push(Interface {
                        patches: (u, v),
                        samples,
                    });
                }
            }
        }
        // Triple junctions: patch triples whose three pairwise interfaces
        // all exist (2-cells of the adjacency complex).
        let mut triples = Vec::new();
        let has_edge =
            |a: usize, b: usize| interfaces.iter().any(|i| i.patches == (a.min(b), a.max(b)));
        for a in 0..n {
            for b in (a + 1)..n {
                for c in (b + 1)..n {
                    if has_edge(a, b) && has_edge(b, c) && has_edge(a, c) {
                        let count = interfaces
                            .iter()
                            .filter(|i| {
                                let (p, q) = i.patches;
                                (p == a || p == b) && (q == b || q == c)
                            })
                            .map(|i| i.samples.len())
                            .min()
                            .unwrap_or(0);
                        triples.push(TripleCell {
                            patches: (a, b, c),
                            samples: count,
                        });
                    }
                }
            }
        }
        SheafComplex {
            n_patches: n,
            interfaces,
            triples,
        }
    }

    /// Assemble δ⁰ (edges × patches) with ±1 entries: one row per
    /// interface SAMPLE, `+1` on patch v's slot, `−1` on patch u's.
    /// (Per-sample rows: the stalk of an edge is its sample space.)
    #[must_use]
    pub fn delta0(&self) -> Csr {
        let rows: usize = self.interfaces.iter().map(|i| i.samples.len()).sum();
        let mut coo = Coo::new(rows, self.n_patches);
        let mut r = 0usize;
        for iface in &self.interfaces {
            for _ in &iface.samples {
                coo.push(r, iface.patches.0, -1.0);
                coo.push(r, iface.patches.1, 1.0);
                r += 1;
            }
        }
        coo.assemble()
    }

    /// Assemble δ¹ (triples × edges) with ±1 entries per the oriented
    /// triangle boundary: for triple (a,b,c) with edges e_ab, e_bc, e_ac:
    /// `+e_ab + e_bc − e_ac` (edge-level stalks: one column per edge).
    #[must_use]
    pub fn delta1(&self) -> Csr {
        let edge_index = |a: usize, b: usize| {
            self.interfaces
                .iter()
                .position(|i| i.patches == (a.min(b), a.max(b)))
                .expect("triple implies edges")
        };
        let mut coo = Coo::new(self.triples.len(), self.interfaces.len());
        for (t, triple) in self.triples.iter().enumerate() {
            let (a, b, c) = triple.patches;
            coo.push(t, edge_index(a, b), 1.0);
            coo.push(t, edge_index(b, c), 1.0);
            coo.push(t, edge_index(a, c), -1.0);
        }
        coo.assemble()
    }

    /// Edge-level δ⁰ (edges × patches, one row per INTERFACE): the
    /// companion of [`Self::delta1`] for the bitwise δδ = 0 identity.
    #[must_use]
    pub fn delta0_edges(&self) -> Csr {
        let mut coo = Coo::new(self.interfaces.len(), self.n_patches);
        for (r, iface) in self.interfaces.iter().enumerate() {
            coo.push(r, iface.patches.0, -1.0);
            coo.push(r, iface.patches.1, 1.0);
        }
        coo.assemble()
    }

    /// Per-interface mismatch assessment. The VERDICT bits come from
    /// fs-ivl's sound predicates (`encloses`/`contains` — no bound
    /// extraction); the reported magnitudes are midpoint±width/2
    /// reconstructions (within an ulp, for humans and logs).
    #[must_use]
    pub fn mismatch_bounds(&self, tol: f64) -> Vec<InterfaceBound> {
        self.interfaces
            .iter()
            .map(|iface| {
                let ok_band = Interval::new(0.0, tol);
                let mut all_within = true;
                let mut proven_leak = false;
                let mut lo = 0.0f64;
                let mut hi = 0.0f64;
                for s in &iface.samples {
                    let d = (s.values[1] - s.values[0]).abs();
                    let within = ok_band.encloses(d);
                    all_within &= within;
                    // |mismatch| enclosure entirely above tol: a proven
                    // violation (sound: the true value is inside d).
                    if !within && !d.contains(tol) {
                        proven_leak = true;
                    }
                    let (m, w) = (d.midpoint(), d.width());
                    lo = lo.max(m - 0.5 * w);
                    hi = hi.max(m + 0.5 * w);
                }
                InterfaceBound {
                    patches: iface.patches,
                    all_within,
                    proven_leak,
                    lo_report: lo,
                    hi_report: hi,
                }
            })
            .collect()
    }

    /// Least-squares section solve: per-patch gauge offsets minimizing
    /// the mean-square mismatch (graph-Laplacian normal equations, gauge
    /// fixed by pinning patch 0; deterministic Gauss–Seidel). Returns
    /// (offsets, raw ms mismatch, residual ms mismatch).
    #[must_use]
    pub fn section_solve(&self) -> (Vec<f64>, f64, f64) {
        let n = self.n_patches;
        let mut offsets = vec![0.0f64; n];
        // Edge means of the midpoint mismatch.
        let edges: Vec<((usize, usize), f64, usize)> = self
            .interfaces
            .iter()
            .map(|iface| {
                let mut sum = 0.0;
                for s in &iface.samples {
                    sum += s.values[1].midpoint() - s.values[0].midpoint();
                }
                (iface.patches, sum, s_len(iface))
            })
            .collect();
        let raw_ms = mean_square(&edges, &offsets);
        for _ in 0..200 {
            for p in 1..n {
                // Optimal c_p given the rest: weighted average balance.
                let mut num = 0.0f64;
                let mut den = 0.0f64;
                for ((u, v), sum, count) in &edges {
                    #[allow(clippy::cast_precision_loss)]
                    let w = *count as f64;
                    if *u == p {
                        num += sum + w * offsets[*v];
                        den += w;
                    } else if *v == p {
                        num += -sum + w * offsets[*u];
                        den += w;
                    }
                }
                if den > 0.0 {
                    offsets[p] = num / den;
                }
            }
        }
        let residual_ms = mean_square(&edges, &offsets);
        (offsets, raw_ms, residual_ms)
    }

    /// The watertightness certificate: interval-verified verdict as
    /// Evidence (enclosure numerics; content-addressed provenance).
    /// PASS requires every sample's enclosure INSIDE `[0, tol]` (sound);
    /// FAIL requires a proven-above-tolerance interface; anything else
    /// is an honest Unknown.
    #[must_use]
    pub fn watertightness(&self, tol: f64) -> Evidence<SheafVerdict> {
        let bounds = self.mismatch_bounds(tol);
        let worst_hi = bounds.iter().map(|b| b.hi_report).fold(0.0f64, f64::max);
        let worst_lo = bounds.iter().map(|b| b.lo_report).fold(0.0f64, f64::max);
        // PASS requires at least one DISCOVERED interface whose samples all lie
        // inside [0, tol]. An empty `bounds` means NO interface was found (charts
        // are disjoint/gapped/near-tangent, or the geometry is empty) — the
        // interface-agreement check gathered NO evidence, so the honest verdict
        // is `Unknown`, not a positive `worst_mismatch = 0` PASS. `all()` on an
        // empty set is vacuously true; guarding on non-emptiness closes the
        // vacuous-truth false certificate (bead obnw).
        let all_pass = !bounds.is_empty() && bounds.iter().all(|b| b.all_within);
        let obstruction: Vec<((usize, usize), f64)> = bounds
            .iter()
            .filter(|b| b.proven_leak)
            .map(|b| (b.patches, b.lo_report))
            .collect();
        let verdict = if all_pass {
            SheafVerdict::Pass {
                worst_mismatch: worst_hi,
                margins: bounds.iter().map(|b| (b.patches, b.hi_report)).collect(),
            }
        } else if !obstruction.is_empty() {
            let (_, raw, residual) = self.section_solve();
            let share = if raw > 0.0 {
                (1.0 - residual / raw).clamp(0.0, 1.0)
            } else {
                0.0
            };
            SheafVerdict::Fail {
                obstruction,
                coboundary_share: share,
            }
        } else {
            SheafVerdict::Unknown {
                straddling: bounds
                    .iter()
                    .filter(|b| !b.all_within)
                    .map(|b| (b.patches, b.lo_report, b.hi_report))
                    .collect(),
            }
        };
        let mut canon = format!(
            "sheaf-watertightness;patches={};interfaces={};tol={tol}",
            self.n_patches,
            self.interfaces.len()
        );
        for b in &bounds {
            let _ = write!(
                canon,
                ";{}-{}:{}:{}",
                b.patches.0, b.patches.1, b.lo_report, b.hi_report
            );
        }
        Evidence {
            qoi: worst_hi,
            numerical: NumericalCertificate::enclosure(worst_lo, worst_hi),
            statistical: StatisticalCertificate::None,
            model: ModelEvidence::none(),
            sensitivity: SensitivitySummary::default(),
            provenance: ProvenanceHash::of_bytes(canon.as_bytes()),
            adjoint_ref: None,
            value: verdict,
        }
    }
}

fn s_len(iface: &Interface) -> usize {
    iface.samples.len()
}

fn mean_square(edges: &[((usize, usize), f64, usize)], offsets: &[f64]) -> f64 {
    let mut num = 0.0f64;
    let mut den = 0.0f64;
    for ((u, v), sum, count) in edges {
        #[allow(clippy::cast_precision_loss)]
        let w = *count as f64;
        // Mean mismatch on the edge after gauging.
        let gauged = sum / w + offsets[*v] - offsets[*u];
        num += w * gauged * gauged;
        den += w;
    }
    if den > 0.0 { num / den } else { 0.0 }
}

/// THE INDEPENDENT FALSIFIER (registry pairing: watertightness →
/// ray-parity): a different algorithm on a different code path. March a
/// segment through the union-of-charts model counting sign changes of
/// the min-SDF; a closed (watertight) model yields an EVEN count on
/// segments with both endpoints outside. Returns the violating ray
/// index, if any.
pub fn ray_parity_falsifier(
    charts: &[&dyn Chart],
    rays: &[(Point3, Point3)],
    steps: usize,
    cx: &Cx<'_>,
) -> Option<usize> {
    for (ri, (a, b)) in rays.iter().enumerate() {
        let mut crossings = 0usize;
        let mut prev_sign = None;
        for k in 0..=steps {
            #[allow(clippy::cast_precision_loss)]
            let t = k as f64 / steps as f64;
            let p = Point3::new(
                a.x + t * (b.x - a.x),
                a.y + t * (b.y - a.y),
                a.z + t * (b.z - a.z),
            );
            let d = charts
                .iter()
                .map(|c| c.eval(p, cx).signed_distance)
                .fold(f64::INFINITY, f64::min);
            let sign = d < 0.0;
            if let Some(ps) = prev_sign
                && ps != sign
            {
                crossings += 1;
            }
            prev_sign = Some(sign);
        }
        if !crossings.is_multiple_of(2) {
            return Some(ri);
        }
    }
    None
}
