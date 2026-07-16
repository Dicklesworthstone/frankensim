//! SHEAF REPAIR (patch Rev L, bead wqd.14; [M] — behind the
//! `sheaf-repair` feature until certifier trials pass): upgrade the sheaf
//! machinery from diagnosis to explicit GAUGE-CORRECTION PLANNING. The current routine
//! sequentially fits the interface mismatch to the patch-coboundary image,
//! then fits that residual to the retained triangle-coboundary image, and
//! retains the final remainder. The fixed-iteration results are Hodge-inspired
//! diagnostics, not a per-result certified orthogonal decomposition. Each
//! output has an INTERPRETATION CONTRACT:
//!
//! - EXACT (`δ⁰c`): algebraically removable from the sampled mismatch by a
//!   patch 0-cochain — a candidate chart/gauge adjustment bounded by each
//!   chart's declared error budget;
//! - COEXACT (`δ¹ᵀw`): circulation-like inconsistency around retained
//!   triple cells. Converter orientation/trace errors are one hypothesis, but
//!   chart/model, junction, sampling, and numerical errors can produce the
//!   same algebraic signature; the decomposition alone does not assign cause;
//! - HARMONIC (the remainder): the part left by the current deterministic
//!   patch-potential and triple-junction projections. Because those numerical
//!   solves have no per-result convergence certificate, a generic remainder is
//!   only a candidate. Calling it H¹ or ruling out gauge repair additionally
//!   requires a retained closed, non-exact witness. It does not by itself prove
//!   that geometry topology must change.
//!
//! Repairs are PROPOSALS. `apply_gauge` only corrects the retained mismatch
//! cochain; it does not mutate or re-evaluate a chart, publish geometry, or
//! prove that a chart-level edit realizes the algebraic correction. Only the
//! algebraic gauge proposal currently has a
//! directly evaluated post-repair seam norm; other proposal kinds retain an
//! unavailable (`+∞`) prediction rather than comparing unlike quantities.
//! Optional Rep-Router reroute costs remain cost estimates, and repairs apply
//! only under an explicit budget.

use crate::router::{CostOracle, RoutePlan, RoutePlanError, RouteRequest, Router};
use crate::sheaf::{
    SHEAF_MAX_CHARTS, SHEAF_MAX_PAIR_CANDIDATES, SHEAF_MAX_RETAINED_INTERFACE_SAMPLES,
    SHEAF_MAX_TRIPLE_CANDIDATES, SheafComplex,
};
use fs_exec::{AdmittedBudget, BudgetConsumption, BudgetRefusal, Cx};
use std::fmt::Write as _;

/// The complex skeleton the decomposition runs over (extractable from a
/// [`SheafComplex`] or built directly for controlled fixtures).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SheafSkeleton {
    /// Patch count.
    pub n_patches: usize,
    /// Interfaces as (u, v) with u < v (edge k orients u → v).
    pub edges: Vec<(usize, usize)>,
    /// Triple junctions (a, b, c) sorted; boundary = +e_ab + e_bc − e_ac.
    pub triangles: Vec<(usize, usize, usize)>,
}

/// Structurally admitted repair skeleton.
///
/// Unlike [`SheafSkeleton`], this type cannot be assembled with unchecked
/// public fields. Its canonical incidence is validated once and then retained
/// immutably, so the fallible incidence operators below do not need to rescan
/// topology or rely on indexing assertions.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AdmittedSheafSkeleton {
    n_patches: usize,
    edges: Vec<(usize, usize)>,
    triangles: Vec<(usize, usize, usize)>,
}

/// Failure to extract a repair skeleton from a public complex.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SheafSkeletonError {
    /// The public complex violates ordering, incidence, range, sample, or
    /// sampling-domain invariants required by the incidence operators.
    MalformedComplex,
    /// At least one patch is required by the repair algebra.
    EmptyComplex,
    /// A validation-work or retained skeleton cardinality exceeds its
    /// defensive ceiling.
    WorkLimit {
        /// Stable validation stage.
        stage: &'static str,
        /// Caller-supplied cardinality.
        requested: usize,
        /// Defensive ceiling.
        cap: usize,
    },
    /// Edges must be strictly increasing canonical pairs in range.
    InvalidEdge {
        /// Position of the first invalid edge.
        index: usize,
    },
    /// Triangles must be strictly increasing canonical triples in range and
    /// every oriented boundary edge must be retained.
    InvalidTriangle {
        /// Position of the first invalid triangle.
        index: usize,
    },
    /// A cochain does not match the admitted incidence space.
    CochainLength {
        /// Stable cochain role.
        role: &'static str,
        /// Required number of scalars.
        expected: usize,
        /// Supplied number of scalars.
        actual: usize,
    },
    /// A cochain scalar is not finite.
    NonFiniteCochain {
        /// Stable cochain role.
        role: &'static str,
        /// First non-finite scalar.
        index: usize,
    },
    /// Finite inputs overflowed an incidence arithmetic operation.
    NumericalOverflow {
        /// Stable arithmetic stage.
        stage: &'static str,
    },
    /// A bounded output allocation could not be reserved.
    ResourceExhausted {
        /// Stable allocation stage.
        stage: &'static str,
    },
}

impl core::fmt::Display for SheafSkeletonError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::MalformedComplex => write!(
                f,
                "cannot extract a skeleton from a malformed sheaf complex"
            ),
            Self::EmptyComplex => write!(f, "repair skeleton requires at least one patch"),
            Self::WorkLimit {
                stage,
                requested,
                cap,
            } => write!(
                f,
                "repair skeleton stage {stage} requests {requested} items above cap {cap}"
            ),
            Self::InvalidEdge { index } => write!(
                f,
                "repair skeleton edge {index} is non-canonical, duplicated, or out of range"
            ),
            Self::InvalidTriangle { index } => write!(
                f,
                "repair skeleton triangle {index} is non-canonical, duplicated, out of range, or missing a boundary edge"
            ),
            Self::CochainLength {
                role,
                expected,
                actual,
            } => write!(
                f,
                "repair {role} cochain requires {expected} scalars, got {actual}"
            ),
            Self::NonFiniteCochain { role, index } => {
                write!(f, "repair {role} cochain scalar {index} is not finite")
            }
            Self::NumericalOverflow { stage } => {
                write!(f, "repair incidence arithmetic overflowed during {stage}")
            }
            Self::ResourceExhausted { stage } => {
                write!(f, "repair incidence could not reserve storage for {stage}")
            }
        }
    }
}

impl std::error::Error for SheafSkeletonError {}

fn validate_finite_cochain(
    values: &[f64],
    expected: usize,
    role: &'static str,
) -> Result<(), SheafSkeletonError> {
    if values.len() != expected {
        return Err(SheafSkeletonError::CochainLength {
            role,
            expected,
            actual: values.len(),
        });
    }
    if let Some(index) = values.iter().position(|value| !value.is_finite()) {
        return Err(SheafSkeletonError::NonFiniteCochain { role, index });
    }
    Ok(())
}

fn zeroed_output(len: usize, stage: &'static str) -> Result<Vec<f64>, SheafSkeletonError> {
    let mut values = Vec::new();
    values
        .try_reserve_exact(len)
        .map_err(|_| SheafSkeletonError::ResourceExhausted { stage })?;
    values.resize(len, 0.0);
    Ok(values)
}

fn validate_skeleton_incidence(
    n_patches: usize,
    edges: &[(usize, usize)],
    triangles: &[(usize, usize, usize)],
) -> Result<(), SheafSkeletonError> {
    if n_patches == 0 {
        return Err(SheafSkeletonError::EmptyComplex);
    }
    for (stage, requested, cap) in [
        ("patches", n_patches, SHEAF_MAX_CHARTS),
        ("edges", edges.len(), SHEAF_MAX_PAIR_CANDIDATES),
        ("triangles", triangles.len(), SHEAF_MAX_TRIPLE_CANDIDATES),
    ] {
        if requested > cap {
            return Err(SheafSkeletonError::WorkLimit {
                stage,
                requested,
                cap,
            });
        }
    }

    let mut previous_edge = None;
    for (index, &(u, v)) in edges.iter().enumerate() {
        if u >= v || v >= n_patches || previous_edge.is_some_and(|previous| previous >= (u, v)) {
            return Err(SheafSkeletonError::InvalidEdge { index });
        }
        previous_edge = Some((u, v));
    }

    let mut previous_triangle = None;
    for (index, &(a, b, c)) in triangles.iter().enumerate() {
        if a >= b
            || b >= c
            || c >= n_patches
            || previous_triangle.is_some_and(|previous| previous >= (a, b, c))
            || edges.binary_search(&(a, b)).is_err()
            || edges.binary_search(&(a, c)).is_err()
            || edges.binary_search(&(b, c)).is_err()
        {
            return Err(SheafSkeletonError::InvalidTriangle { index });
        }
        previous_triangle = Some((a, b, c));
    }
    Ok(())
}

fn validate_raw_skeleton_shape(skeleton: &SheafSkeleton) -> Result<(), SheafSkeletonError> {
    if skeleton.n_patches == 0 {
        return Err(SheafSkeletonError::EmptyComplex);
    }
    for (stage, requested, cap) in [
        ("patches", skeleton.n_patches, SHEAF_MAX_CHARTS),
        ("edges", skeleton.edges.len(), SHEAF_MAX_PAIR_CANDIDATES),
        (
            "triangles",
            skeleton.triangles.len(),
            SHEAF_MAX_TRIPLE_CANDIDATES,
        ),
    ] {
        if requested > cap {
            return Err(SheafSkeletonError::WorkLimit {
                stage,
                requested,
                cap,
            });
        }
    }

    for (index, &(u, v)) in skeleton.edges.iter().enumerate() {
        if u >= v || v >= skeleton.n_patches {
            return Err(SheafSkeletonError::InvalidEdge { index });
        }
    }
    for (index, &(a, b, c)) in skeleton.triangles.iter().enumerate() {
        if a >= b || b >= c || c >= skeleton.n_patches {
            return Err(SheafSkeletonError::InvalidTriangle { index });
        }
    }
    Ok(())
}

fn validate_raw_skeleton_cross_structure(
    skeleton: &SheafSkeleton,
) -> Result<(), SheafSkeletonError> {
    let mut indexed_edges = Vec::new();
    indexed_edges
        .try_reserve_exact(skeleton.edges.len())
        .map_err(|_| SheafSkeletonError::ResourceExhausted {
            stage: "raw-validation-edges",
        })?;
    indexed_edges.extend(
        skeleton
            .edges
            .iter()
            .copied()
            .enumerate()
            .map(|(index, edge)| (edge, index)),
    );
    indexed_edges.sort_unstable();
    if let Some(index) = indexed_edges
        .windows(2)
        .filter(|pair| pair[0].0 == pair[1].0)
        .map(|pair| pair[1].1)
        .min()
    {
        return Err(SheafSkeletonError::InvalidEdge { index });
    }

    let mut indexed_triangles = Vec::new();
    indexed_triangles
        .try_reserve_exact(skeleton.triangles.len())
        .map_err(|_| SheafSkeletonError::ResourceExhausted {
            stage: "raw-validation-triangles",
        })?;
    indexed_triangles.extend(
        skeleton
            .triangles
            .iter()
            .copied()
            .enumerate()
            .map(|(index, triangle)| (triangle, index)),
    );
    indexed_triangles.sort_unstable();
    if let Some(index) = indexed_triangles
        .windows(2)
        .filter(|pair| pair[0].0 == pair[1].0)
        .map(|pair| pair[1].1)
        .min()
    {
        return Err(SheafSkeletonError::InvalidTriangle { index });
    }

    for (index, &(a, b, c)) in skeleton.triangles.iter().enumerate() {
        for boundary in [(a, b), (a, c), (b, c)] {
            if indexed_edges
                .binary_search_by_key(&boundary, |(edge, _)| *edge)
                .is_err()
            {
                return Err(SheafSkeletonError::InvalidTriangle { index });
            }
        }
    }
    Ok(())
}

impl AdmittedSheafSkeleton {
    /// Validate and seal canonical repair incidence supplied by a caller.
    ///
    /// Validation order is deterministic: cardinalities, edges in input
    /// order, then triangles in input order. Edges and triangles must already
    /// be strictly lexicographically ordered, making duplicate and orientation
    /// errors unambiguous and every later lookup bounded.
    pub fn try_new(
        n_patches: usize,
        edges: Vec<(usize, usize)>,
        triangles: Vec<(usize, usize, usize)>,
    ) -> Result<Self, SheafSkeletonError> {
        validate_skeleton_incidence(n_patches, &edges, &triangles)?;

        Ok(Self {
            n_patches,
            edges,
            triangles,
        })
    }

    /// Validate and seal a raw repair skeleton without copying its storage.
    pub fn admit(raw: SheafSkeleton) -> Result<Self, SheafSkeletonError> {
        Self::try_new(raw.n_patches, raw.edges, raw.triangles)
    }

    /// Structurally validate a public complex and retain its canonical edge
    /// incidence. Candidate clique triples remain omitted because the base
    /// builder has not verified a common triple overlap.
    pub fn of(complex: &SheafComplex) -> Result<Self, SheafSkeletonError> {
        let raw = SheafSkeleton::of(complex)?;
        Self::try_new(raw.n_patches, raw.edges, raw.triangles)
    }

    /// Number of retained patches.
    #[must_use]
    pub const fn n_patches(&self) -> usize {
        self.n_patches
    }

    /// Canonically ordered retained interfaces.
    #[must_use]
    pub fn edges(&self) -> &[(usize, usize)] {
        &self.edges
    }

    /// Canonically ordered retained triangles.
    #[must_use]
    pub fn triangles(&self) -> &[(usize, usize, usize)] {
        &self.triangles
    }

    /// Apply `delta^0` to a finite vertex cochain.
    pub fn d0(&self, c: &[f64]) -> Result<Vec<f64>, SheafSkeletonError> {
        validate_finite_cochain(c, self.n_patches, "vertex")?;
        let mut out = zeroed_output(self.edges.len(), "d0-output")?;
        for (value, &(u, v)) in out.iter_mut().zip(&self.edges) {
            *value = c[v] - c[u];
            if !value.is_finite() {
                return Err(SheafSkeletonError::NumericalOverflow { stage: "d0" });
            }
        }
        Ok(out)
    }

    /// Apply `delta^0` transpose to a finite edge cochain.
    pub fn d0t(&self, m: &[f64]) -> Result<Vec<f64>, SheafSkeletonError> {
        validate_finite_cochain(m, self.edges.len(), "edge")?;
        let mut out = zeroed_output(self.n_patches, "d0t-output")?;
        for (k, &(u, v)) in self.edges.iter().enumerate() {
            out[u] -= m[k];
            out[v] += m[k];
            if !(out[u].is_finite() && out[v].is_finite()) {
                return Err(SheafSkeletonError::NumericalOverflow { stage: "d0t" });
            }
        }
        Ok(out)
    }

    /// Apply `delta^1` to a finite edge cochain.
    pub fn d1(&self, m: &[f64]) -> Result<Vec<f64>, SheafSkeletonError> {
        validate_finite_cochain(m, self.edges.len(), "edge")?;
        let mut out = zeroed_output(self.triangles.len(), "d1-output")?;
        for (triangle, (value, &(a, b, c))) in out.iter_mut().zip(&self.triangles).enumerate() {
            let eab = self
                .edges
                .binary_search(&(a, b))
                .map_err(|_| SheafSkeletonError::InvalidTriangle { index: triangle })?;
            let ebc = self
                .edges
                .binary_search(&(b, c))
                .map_err(|_| SheafSkeletonError::InvalidTriangle { index: triangle })?;
            let eac = self
                .edges
                .binary_search(&(a, c))
                .map_err(|_| SheafSkeletonError::InvalidTriangle { index: triangle })?;
            *value = (m[eab] + m[ebc]) - m[eac];
            if !value.is_finite() {
                return Err(SheafSkeletonError::NumericalOverflow { stage: "d1" });
            }
        }
        Ok(out)
    }

    /// Apply `delta^1` transpose to a finite triangle cochain.
    pub fn d1t(&self, w: &[f64]) -> Result<Vec<f64>, SheafSkeletonError> {
        validate_finite_cochain(w, self.triangles.len(), "triangle")?;
        let mut out = zeroed_output(self.edges.len(), "d1t-output")?;
        for (triangle, &(a, b, c)) in self.triangles.iter().enumerate() {
            let eab = self
                .edges
                .binary_search(&(a, b))
                .map_err(|_| SheafSkeletonError::InvalidTriangle { index: triangle })?;
            let ebc = self
                .edges
                .binary_search(&(b, c))
                .map_err(|_| SheafSkeletonError::InvalidTriangle { index: triangle })?;
            let eac = self
                .edges
                .binary_search(&(a, c))
                .map_err(|_| SheafSkeletonError::InvalidTriangle { index: triangle })?;
            out[eab] += w[triangle];
            out[ebc] += w[triangle];
            out[eac] -= w[triangle];
            if !(out[eab].is_finite() && out[ebc].is_finite() && out[eac].is_finite()) {
                return Err(SheafSkeletonError::NumericalOverflow { stage: "d1t" });
            }
        }
        Ok(out)
    }
}

impl SheafSkeleton {
    /// Structurally validate and extract caller-supplied adjacency for
    /// diagnostic repair algebra. This accepts a raw public complex and does not
    /// authenticate chart-sampling origin or confer topology authority. The base
    /// builder's `TripleCell`s are pairwise-interface clique completions rather
    /// than verified common triple overlaps, so they are deliberately omitted.
    ///
    /// # Errors
    /// Returns a deterministic empty, work-limit, structural, or allocation
    /// refusal before publishing a skeleton.
    pub fn of(complex: &SheafComplex) -> Result<SheafSkeleton, SheafSkeletonError> {
        if complex.n_patches == 0 {
            return Err(SheafSkeletonError::EmptyComplex);
        }
        for (stage, requested, cap) in [
            ("patches", complex.n_patches, SHEAF_MAX_CHARTS),
            (
                "interfaces",
                complex.interfaces.len(),
                SHEAF_MAX_PAIR_CANDIDATES,
            ),
            (
                "triples",
                complex.triples.len(),
                SHEAF_MAX_TRIPLE_CANDIDATES,
            ),
        ] {
            if requested > cap {
                return Err(SheafSkeletonError::WorkLimit {
                    stage,
                    requested,
                    cap,
                });
            }
        }
        let mut samples = 0usize;
        for interface in &complex.interfaces {
            samples = samples.checked_add(interface.samples.len()).ok_or(
                SheafSkeletonError::WorkLimit {
                    stage: "interface-samples",
                    requested: usize::MAX,
                    cap: SHEAF_MAX_RETAINED_INTERFACE_SAMPLES,
                },
            )?;
            if samples > SHEAF_MAX_RETAINED_INTERFACE_SAMPLES {
                return Err(SheafSkeletonError::WorkLimit {
                    stage: "interface-samples",
                    requested: samples,
                    cap: SHEAF_MAX_RETAINED_INTERFACE_SAMPLES,
                });
            }
        }
        if !complex.structure_is_valid() {
            return Err(SheafSkeletonError::MalformedComplex);
        }
        let mut edges = Vec::new();
        edges
            .try_reserve_exact(complex.interfaces.len())
            .map_err(|_| SheafSkeletonError::ResourceExhausted {
                stage: "raw-complex-edges",
            })?;
        edges.extend(complex.interfaces.iter().map(|interface| interface.patches));
        Ok(SheafSkeleton {
            n_patches: complex.n_patches,
            edges,
            triangles: Vec::new(),
        })
    }

    fn edge_index(&self, a: usize, b: usize) -> Option<usize> {
        let key = (a.min(b), a.max(b));
        self.edges.iter().position(|&e| e == key)
    }

    /// Apply δ⁰ to a vertex cochain: `(δ⁰c)_e = c_v − c_u`.
    ///
    /// # Errors
    /// Returns a deterministic structural, cardinality, finiteness,
    /// allocation, or arithmetic refusal without indexing unchecked input.
    pub fn d0(&self, c: &[f64]) -> Result<Vec<f64>, SheafSkeletonError> {
        validate_raw_skeleton_shape(self)?;
        validate_finite_cochain(c, self.n_patches, "vertex")?;
        validate_raw_skeleton_cross_structure(self)?;
        self.d0_validated(c)
    }

    /// Apply δ⁰ᵀ to an edge cochain.
    ///
    /// # Errors
    /// Returns a deterministic structural, cardinality, finiteness,
    /// allocation, or arithmetic refusal without indexing unchecked input.
    pub fn d0t(&self, m: &[f64]) -> Result<Vec<f64>, SheafSkeletonError> {
        validate_raw_skeleton_shape(self)?;
        validate_finite_cochain(m, self.edges.len(), "edge")?;
        validate_raw_skeleton_cross_structure(self)?;
        self.d0t_validated(m)
    }

    /// Apply δ¹ to an edge cochain: signed sum around each triangle.
    ///
    /// # Errors
    /// Returns a deterministic structural, cardinality, finiteness,
    /// allocation, or arithmetic refusal without indexing unchecked input.
    pub fn d1(&self, m: &[f64]) -> Result<Vec<f64>, SheafSkeletonError> {
        validate_raw_skeleton_shape(self)?;
        validate_finite_cochain(m, self.edges.len(), "edge")?;
        validate_raw_skeleton_cross_structure(self)?;
        self.d1_validated(m)
    }

    /// Apply δ¹ᵀ to a triangle cochain.
    ///
    /// # Errors
    /// Returns a deterministic structural, cardinality, finiteness,
    /// allocation, or arithmetic refusal without indexing unchecked input.
    pub fn d1t(&self, w: &[f64]) -> Result<Vec<f64>, SheafSkeletonError> {
        validate_raw_skeleton_shape(self)?;
        validate_finite_cochain(w, self.triangles.len(), "triangle")?;
        validate_raw_skeleton_cross_structure(self)?;
        self.d1t_validated(w)
    }

    fn d0_validated(&self, c: &[f64]) -> Result<Vec<f64>, SheafSkeletonError> {
        let mut out = zeroed_output(self.edges.len(), "raw-d0-output")?;
        for (value, &(u, v)) in out.iter_mut().zip(&self.edges) {
            *value = c[v] - c[u];
            if !value.is_finite() {
                return Err(SheafSkeletonError::NumericalOverflow { stage: "d0" });
            }
        }
        Ok(out)
    }

    fn d0t_validated(&self, m: &[f64]) -> Result<Vec<f64>, SheafSkeletonError> {
        let mut out = zeroed_output(self.n_patches, "raw-d0t-output")?;
        for (k, &(u, v)) in self.edges.iter().enumerate() {
            out[u] -= m[k];
            out[v] += m[k];
            if !(out[u].is_finite() && out[v].is_finite()) {
                return Err(SheafSkeletonError::NumericalOverflow { stage: "d0t" });
            }
        }
        Ok(out)
    }

    fn d1_validated(&self, m: &[f64]) -> Result<Vec<f64>, SheafSkeletonError> {
        let mut out = zeroed_output(self.triangles.len(), "raw-d1-output")?;
        for (triangle, (value, &(a, b, c))) in out.iter_mut().zip(&self.triangles).enumerate() {
            let eab = self
                .edge_index(a, b)
                .ok_or(SheafSkeletonError::InvalidTriangle { index: triangle })?;
            let ebc = self
                .edge_index(b, c)
                .ok_or(SheafSkeletonError::InvalidTriangle { index: triangle })?;
            let eac = self
                .edge_index(a, c)
                .ok_or(SheafSkeletonError::InvalidTriangle { index: triangle })?;
            *value = (m[eab] + m[ebc]) - m[eac];
            if !value.is_finite() {
                return Err(SheafSkeletonError::NumericalOverflow { stage: "d1" });
            }
        }
        Ok(out)
    }

    fn d1t_validated(&self, w: &[f64]) -> Result<Vec<f64>, SheafSkeletonError> {
        let mut out = zeroed_output(self.edges.len(), "raw-d1t-output")?;
        for (t, &(a, b, c)) in self.triangles.iter().enumerate() {
            let eab = self
                .edge_index(a, b)
                .ok_or(SheafSkeletonError::InvalidTriangle { index: t })?;
            let ebc = self
                .edge_index(b, c)
                .ok_or(SheafSkeletonError::InvalidTriangle { index: t })?;
            let eac = self
                .edge_index(a, c)
                .ok_or(SheafSkeletonError::InvalidTriangle { index: t })?;
            out[eab] += w[t];
            out[ebc] += w[t];
            out[eac] -= w[t];
            if !(out[eab].is_finite() && out[ebc].is_finite() && out[eac].is_finite()) {
                return Err(SheafSkeletonError::NumericalOverflow { stage: "d1t" });
            }
        }
        Ok(out)
    }
}

/// The Hodge-inspired sequential diagnostic split of an edge mismatch cochain.
#[derive(Debug, Clone, PartialEq)]
pub struct HodgeSplit {
    /// The fitted exact (coboundary) component `δ⁰c`.
    pub exact: Vec<f64>,
    /// The vertex potential `c` (gauge offsets; `c[0]` pinned to 0).
    pub potential: Vec<f64>,
    /// The fitted coexact component `δ¹ᵀw` of the first residual.
    pub coexact: Vec<f64>,
    /// The remainder retained after both fixed-iteration fits.
    pub harmonic: Vec<f64>,
    /// Separate squared-norm ratios (exact, coexact, remainder) over ‖m‖².
    /// Without certified orthogonality these diagnostic ratios need not sum to
    /// one.
    pub fractions: (f64, f64, f64),
}

/// Explicit resource envelope for the admitted Hodge diagnostic.
///
/// `sweeps` is the exact deterministic sweep count used by each non-empty
/// least-squares stage. The remaining fields are hard admission ceilings, not
/// post-hoc observations. Deadline and capability budgets remain carried by
/// the caller's [`Cx`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SheafRepairBudget {
    /// Deterministic Gauss-Seidel sweeps per non-empty projection stage.
    pub sweeps: usize,
    /// Maximum incidence-operator applications across the entire split.
    pub max_operator_evaluations: usize,
    /// Maximum scalar, graph, comparison, and formatting work items.
    pub max_work_items: usize,
    /// Conservative maximum number of simultaneously live/retained scalar
    /// slots admitted for the split.
    pub max_scalar_slots: usize,
    /// Maximum scalar, graph, comparison, or formatting work items between
    /// cancellation checkpoints.
    pub poll_stride: usize,
}

/// Measured consumption retained with one admitted diagnostic.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SheafRepairUsage {
    /// Projection sweeps completed across exact and coexact stages.
    pub completed_sweeps: usize,
    /// Incidence-operator applications performed.
    pub operator_evaluations: usize,
    /// Scalar, graph, comparison, and formatting work items completed.
    pub work_items: usize,
    /// Conservative work-item envelope admitted before work began.
    pub admitted_work_items: usize,
    /// Conservative scalar-slot envelope admitted before work began.
    pub admitted_scalar_slots: usize,
    /// Enforced ambient deadline, poll, and cost-budget consumption.
    pub ambient_budget: BudgetConsumption,
}

/// A bounded diagnostic plus its resource consumption.
#[derive(Debug, Clone, PartialEq)]
pub struct BoundedHodgeSplit {
    /// The Hodge-inspired diagnostic payload.
    pub split: HodgeSplit,
    /// Exact caller-admitted resource envelope used for this run.
    pub budget: SheafRepairBudget,
    /// Exact measured work and admitted memory envelope.
    pub usage: SheafRepairUsage,
}

/// Additional output-allocation limits for bounded repair planning.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SheafRepairPlanBudget {
    /// Decomposition, scalar-work, cancellation-poll, and ambient budget.
    pub repair: SheafRepairBudget,
    /// Maximum cumulative bytes requested for plan-only retained storage.
    pub max_plan_bytes: usize,
    /// Maximum UTF-8 bytes retained across every proposal action.
    pub max_action_bytes: usize,
    /// Maximum proposals that may be published.
    pub max_proposals: usize,
    /// Maximum retained harmonic-support interfaces.
    pub max_harmonic_support: usize,
}

/// Measured consumption retained with one bounded repair plan.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SheafRepairPlanUsage {
    /// Shared decomposition and whole-plan execution accounting.
    pub repair: SheafRepairUsage,
    /// Conservative plan-only byte envelope admitted before work began.
    pub plan_memory_envelope: usize,
    /// Cumulative plan-only capacity bytes actually requested.
    pub reserved_plan_bytes: usize,
    /// UTF-8 action bytes retained by the published proposals.
    pub action_bytes: usize,
    /// Proposals retained in the published plan.
    pub proposals: usize,
    /// Interfaces retained in the harmonic-support diagnostic.
    pub harmonic_support: usize,
}

/// A repair plan published only after all bounded work completes.
#[derive(Debug, Clone, PartialEq)]
pub struct BoundedRepairPlan {
    /// Complete repair proposal payload.
    pub plan: RepairPlan,
    /// Exact caller-admitted resource envelope.
    pub budget: SheafRepairPlanBudget,
    /// Enforced and measured resource consumption.
    pub usage: SheafRepairPlanUsage,
}

/// Structured refusal from the admitted, cancellation-aware repair path.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SheafRepairError {
    /// The admitted skeleton or cochain failed a total incidence operation.
    Skeleton(SheafSkeletonError),
    /// One per-patch gauge budget is non-finite or negative.
    InvalidGaugeBudget {
        /// First invalid budget in caller order.
        index: usize,
    },
    /// A budget field that must be positive was zero.
    InvalidBudget {
        /// Stable budget field name.
        field: &'static str,
    },
    /// The deterministic operator schedule exceeds the admitted work cap.
    WorkBudgetExceeded {
        /// Conservative required operator-application envelope.
        required: u128,
        /// Caller-admitted ceiling.
        cap: usize,
    },
    /// The conservative scalar/graph/string work schedule exceeds its cap.
    WorkItemBudgetExceeded {
        /// Stable preflight or execution stage.
        stage: &'static str,
        /// Conservative or next required work-item count.
        required: u128,
        /// Caller-admitted ceiling.
        cap: usize,
    },
    /// The conservative live/retained scalar envelope exceeds the admitted
    /// memory cap.
    MemoryBudgetExceeded {
        /// Conservative admitted scalar-slot envelope.
        required: u128,
        /// Caller-admitted ceiling.
        cap: usize,
    },
    /// The conservative plan-only allocation envelope exceeds its byte cap.
    PlanMemoryBudgetExceeded {
        /// Conservative required plan-only bytes.
        required: u128,
        /// Caller-admitted byte ceiling.
        cap: usize,
    },
    /// A retained plan output exceeds its explicit cardinality/byte cap.
    OutputBudgetExceeded {
        /// Stable output resource name.
        resource: &'static str,
        /// Required output cardinality or bytes.
        required: u128,
        /// Caller-admitted ceiling.
        cap: usize,
    },
    /// Checked admission arithmetic exceeded `u128`.
    BudgetArithmeticOverflow {
        /// Stable preflight stage.
        stage: &'static str,
    },
    /// Caller cancellation was observed at a bounded checkpoint.
    Cancelled {
        /// Stable execution stage.
        stage: &'static str,
        /// Sweeps fully completed before cancellation.
        completed_sweeps: usize,
        /// Operator applications completed before cancellation.
        operator_evaluations: usize,
        /// Scalar/graph/string work items completed before cancellation.
        work_items: usize,
    },
    /// The ambient deadline, poll, or cost budget refused execution.
    AmbientBudgetRefused {
        /// First refusal latched by the ambient accountant.
        refusal: BudgetRefusal,
        /// Sweeps fully completed before refusal.
        completed_sweeps: usize,
        /// Operator applications completed before refusal.
        operator_evaluations: usize,
        /// Scalar/graph/string work items completed before refusal.
        work_items: usize,
    },
    /// Finite diagnostic arithmetic overflowed.
    NumericalOverflow {
        /// Stable arithmetic stage.
        stage: &'static str,
    },
}

impl core::fmt::Display for SheafRepairError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Skeleton(source) => write!(f, "{source}"),
            Self::InvalidGaugeBudget { index } => write!(
                f,
                "sheaf repair gauge budget {index} must be finite and non-negative"
            ),
            Self::InvalidBudget { field } => {
                write!(f, "sheaf repair budget field {field} must be positive")
            }
            Self::WorkBudgetExceeded { required, cap } => write!(
                f,
                "sheaf repair operator envelope requires {required} evaluations above cap {cap}"
            ),
            Self::WorkItemBudgetExceeded {
                stage,
                required,
                cap,
            } => write!(
                f,
                "sheaf repair work envelope requires {required} items above cap {cap} during {stage}"
            ),
            Self::MemoryBudgetExceeded { required, cap } => write!(
                f,
                "sheaf repair scalar envelope requires {required} slots above cap {cap}"
            ),
            Self::PlanMemoryBudgetExceeded { required, cap } => write!(
                f,
                "sheaf repair plan envelope requires {required} bytes above cap {cap}"
            ),
            Self::OutputBudgetExceeded {
                resource,
                required,
                cap,
            } => write!(
                f,
                "sheaf repair output {resource} requires {required} above cap {cap}"
            ),
            Self::BudgetArithmeticOverflow { stage } => {
                write!(
                    f,
                    "sheaf repair budget arithmetic overflowed during {stage}"
                )
            }
            Self::Cancelled {
                stage,
                completed_sweeps,
                operator_evaluations,
                work_items,
            } => write!(
                f,
                "sheaf repair cancelled during {stage} after {completed_sweeps} sweeps, {operator_evaluations} operator evaluations, and {work_items} work items"
            ),
            Self::AmbientBudgetRefused {
                refusal,
                completed_sweeps,
                operator_evaluations,
                work_items,
            } => write!(
                f,
                "sheaf repair ambient budget refused after {completed_sweeps} sweeps, {operator_evaluations} operator evaluations, and {work_items} work items: {refusal}"
            ),
            Self::NumericalOverflow { stage } => {
                write!(f, "sheaf repair arithmetic overflowed during {stage}")
            }
        }
    }
}

impl std::error::Error for SheafRepairError {}

impl From<SheafSkeletonError> for SheafRepairError {
    fn from(source: SheafSkeletonError) -> Self {
        Self::Skeleton(source)
    }
}

pub(super) struct RepairAccountant<'a, 'cx> {
    cx: &'a Cx<'cx>,
    budget: SheafRepairBudget,
    ambient: AdmittedBudget<'cx>,
    operator_evaluations: usize,
    completed_sweeps: usize,
    work_items: usize,
    plan_memory_cap: usize,
    reserved_plan_bytes: usize,
    action_bytes_cap: usize,
    action_bytes: usize,
}

impl<'a, 'cx> RepairAccountant<'a, 'cx> {
    pub(super) fn new(
        cx: &'a Cx<'cx>,
        budget: SheafRepairBudget,
        planned_cost: u64,
        plan_memory_cap: usize,
        action_bytes_cap: usize,
    ) -> Result<Self, SheafRepairError> {
        let ambient = AdmittedBudget::admit_ambient(cx, planned_cost).map_err(|refusal| {
            SheafRepairError::AmbientBudgetRefused {
                refusal,
                completed_sweeps: 0,
                operator_evaluations: 0,
                work_items: 0,
            }
        })?;
        Ok(Self {
            cx,
            budget,
            ambient,
            operator_evaluations: 0,
            completed_sweeps: 0,
            work_items: 0,
            plan_memory_cap,
            reserved_plan_bytes: 0,
            action_bytes_cap,
            action_bytes: 0,
        })
    }

    fn map_refusal(&self, stage: &'static str, refusal: BudgetRefusal) -> SheafRepairError {
        if matches!(refusal, BudgetRefusal::Cancelled { .. }) {
            SheafRepairError::Cancelled {
                stage,
                completed_sweeps: self.completed_sweeps,
                operator_evaluations: self.operator_evaluations,
                work_items: self.work_items,
            }
        } else {
            SheafRepairError::AmbientBudgetRefused {
                refusal,
                completed_sweeps: self.completed_sweeps,
                operator_evaluations: self.operator_evaluations,
                work_items: self.work_items,
            }
        }
    }

    pub(super) fn checkpoint(&mut self, stage: &'static str) -> Result<(), SheafRepairError> {
        let result = self.ambient.checkpoint(stage, self.cx);
        result.map_err(|refusal| self.map_refusal(stage, refusal))
    }

    fn begin_operator(&mut self, stage: &'static str) -> Result<(), SheafRepairError> {
        if self.operator_evaluations >= self.budget.max_operator_evaluations {
            return Err(SheafRepairError::WorkBudgetExceeded {
                required: self.operator_evaluations as u128 + 1,
                cap: self.budget.max_operator_evaluations,
            });
        }
        self.checkpoint(stage)?;
        Ok(())
    }

    pub(super) fn consume_item(&mut self, stage: &'static str) -> Result<(), SheafRepairError> {
        if self.work_items >= self.budget.max_work_items {
            return Err(SheafRepairError::WorkItemBudgetExceeded {
                stage,
                required: self.work_items as u128 + 1,
                cap: self.budget.max_work_items,
            });
        }
        if self.work_items.is_multiple_of(self.budget.poll_stride) {
            self.checkpoint(stage)?;
        }
        let next =
            self.work_items
                .checked_add(1)
                .ok_or(SheafRepairError::BudgetArithmeticOverflow {
                    stage: "work-items",
                })?;
        let charge = self.ambient.charge_cost(stage, 1);
        let charge = charge.map_err(|refusal| self.map_refusal(stage, refusal));
        charge?;
        self.work_items = next;
        Ok(())
    }

    fn finish_operator(&mut self, stage: &'static str) -> Result<(), SheafRepairError> {
        self.operator_evaluations = self.operator_evaluations.checked_add(1).ok_or(
            SheafRepairError::BudgetArithmeticOverflow {
                stage: "operator-evaluations",
            },
        )?;
        let charge = self.ambient.charge_cost(stage, 1);
        let charge = charge.map_err(|refusal| self.map_refusal(stage, refusal));
        charge?;
        Ok(())
    }

    fn complete_sweep(&mut self, stage: &'static str) -> Result<(), SheafRepairError> {
        self.completed_sweeps = self.completed_sweeps.checked_add(1).ok_or(
            SheafRepairError::BudgetArithmeticOverflow {
                stage: "completed-sweeps",
            },
        )?;
        self.checkpoint(stage)
    }

    pub(super) fn reserve_plan_bytes(
        &mut self,
        stage: &'static str,
        bytes: usize,
    ) -> Result<(), SheafRepairError> {
        let required = self.reserved_plan_bytes.checked_add(bytes).ok_or(
            SheafRepairError::BudgetArithmeticOverflow {
                stage: "plan-reserved-bytes",
            },
        )?;
        if required > self.plan_memory_cap {
            return Err(SheafRepairError::PlanMemoryBudgetExceeded {
                required: required as u128,
                cap: self.plan_memory_cap,
            });
        }
        self.checkpoint(stage)?;
        self.reserved_plan_bytes = required;
        Ok(())
    }

    pub(super) fn release_plan_bytes(
        &mut self,
        stage: &'static str,
        bytes: usize,
    ) -> Result<(), SheafRepairError> {
        let remaining = self.reserved_plan_bytes.checked_sub(bytes).ok_or(
            SheafRepairError::BudgetArithmeticOverflow {
                stage: "plan-released-bytes",
            },
        )?;
        self.checkpoint(stage)?;
        self.reserved_plan_bytes = remaining;
        Ok(())
    }

    fn retain_action_bytes(
        &mut self,
        stage: &'static str,
        bytes: usize,
    ) -> Result<(), SheafRepairError> {
        let required = self.action_bytes.checked_add(bytes).ok_or(
            SheafRepairError::BudgetArithmeticOverflow {
                stage: "action-bytes",
            },
        )?;
        if required > self.action_bytes_cap {
            return Err(SheafRepairError::OutputBudgetExceeded {
                resource: "action-bytes",
                required: required as u128,
                cap: self.action_bytes_cap,
            });
        }
        self.reserve_plan_bytes(stage, bytes)?;
        self.action_bytes = required;
        Ok(())
    }

    pub(super) fn usage(
        &self,
        admitted_scalar_slots: usize,
        admitted_work_items: usize,
    ) -> SheafRepairUsage {
        SheafRepairUsage {
            completed_sweeps: self.completed_sweeps,
            operator_evaluations: self.operator_evaluations,
            work_items: self.work_items,
            admitted_work_items,
            admitted_scalar_slots,
            ambient_budget: self.ambient.consumption(),
        }
    }

    pub(super) const fn reserved_plan_bytes(&self) -> usize {
        self.reserved_plan_bytes
    }
}

fn checked_norm2(
    values: &[f64],
    stage: &'static str,
    accountant: &mut RepairAccountant<'_, '_>,
) -> Result<f64, SheafRepairError> {
    let mut total = 0.0f64;
    for value in values {
        accountant.consume_item(stage)?;
        let square = value * value;
        total += square;
        if !(square.is_finite() && total.is_finite()) {
            return Err(SheafRepairError::NumericalOverflow { stage });
        }
    }
    Ok(total)
}

pub(super) fn zeroed_output_bounded(
    len: usize,
    stage: &'static str,
    accountant: &mut RepairAccountant<'_, '_>,
) -> Result<Vec<f64>, SheafRepairError> {
    accountant.checkpoint(stage)?;
    let mut output = Vec::new();
    output
        .try_reserve_exact(len)
        .map_err(|_| SheafSkeletonError::ResourceExhausted { stage })?;
    for _ in 0..len {
        accountant.consume_item(stage)?;
        output.push(0.0);
    }
    accountant.checkpoint(stage)?;
    Ok(output)
}

fn checked_scalar_envelope(skeleton: &AdmittedSheafSkeleton) -> Result<u128, SheafRepairError> {
    let dimensions = (skeleton.n_patches as u128)
        .checked_add(skeleton.edges.len() as u128)
        .and_then(|value| value.checked_add(skeleton.triangles.len() as u128))
        .ok_or(SheafRepairError::BudgetArithmeticOverflow {
            stage: "scalar-envelope-dimensions",
        })?;
    dimensions
        .checked_mul(6)
        .ok_or(SheafRepairError::BudgetArithmeticOverflow {
            stage: "scalar-envelope",
        })
}

fn projection_operator_evaluations(
    unknowns: usize,
    sweeps: usize,
    pin_first: bool,
) -> Result<u128, SheafRepairError> {
    let active = unknowns.saturating_sub(usize::from(pin_first && unknowns > 0)) as u128;
    let sweep_work = active
        .checked_mul(sweeps as u128)
        .and_then(|value| value.checked_mul(2))
        .ok_or(SheafRepairError::BudgetArithmeticOverflow {
            stage: "projection-operator-evaluations",
        })?;
    1u128
        .checked_add(unknowns as u128)
        .and_then(|value| value.checked_add(sweep_work))
        .and_then(|value| value.checked_add(1))
        .ok_or(SheafRepairError::BudgetArithmeticOverflow {
            stage: "projection-operator-evaluations",
        })
}

fn checked_operator_schedule(
    skeleton: &AdmittedSheafSkeleton,
    sweeps: usize,
) -> Result<u128, SheafRepairError> {
    let exact = projection_operator_evaluations(skeleton.n_patches, sweeps, true)?;
    if skeleton.triangles.is_empty() {
        return Ok(exact);
    }
    exact
        .checked_add(projection_operator_evaluations(
            skeleton.triangles.len(),
            sweeps,
            false,
        )?)
        .ok_or(SheafRepairError::BudgetArithmeticOverflow {
            stage: "hodge-operator-schedule",
        })
}

fn checked_hodge_work_envelope(
    skeleton: &AdmittedSheafSkeleton,
    required_operators: u128,
) -> Result<u128, SheafRepairError> {
    let span = (skeleton.n_patches as u128)
        .checked_add(skeleton.edges.len() as u128)
        .and_then(|value| value.checked_add(skeleton.triangles.len() as u128))
        .and_then(|value| value.checked_add(1))
        .ok_or(SheafRepairError::BudgetArithmeticOverflow {
            stage: "hodge-work-span",
        })?;
    let edge_search_steps = if skeleton.edges.is_empty() {
        0
    } else {
        usize::BITS - skeleton.edges.len().leading_zeros()
    } as u128;
    let per_span = edge_search_steps
        .checked_mul(3)
        .and_then(|value| value.checked_add(4))
        .ok_or(SheafRepairError::BudgetArithmeticOverflow {
            stage: "hodge-work-search-factor",
        })?;
    required_operators
        .checked_add(16)
        .and_then(|value| value.checked_mul(span))
        .and_then(|value| value.checked_mul(per_span))
        .ok_or(SheafRepairError::BudgetArithmeticOverflow {
            stage: "hodge-work-envelope",
        })
}

#[derive(Clone, Copy)]
pub(super) struct RepairAdmission {
    pub(super) scalar_slots: usize,
    pub(super) operator_evaluations: usize,
    pub(super) work_items: usize,
}

pub(super) fn admit_repair_budget(
    skeleton: &AdmittedSheafSkeleton,
    budget: SheafRepairBudget,
) -> Result<RepairAdmission, SheafRepairError> {
    if budget.sweeps == 0 {
        return Err(SheafRepairError::InvalidBudget { field: "sweeps" });
    }
    if budget.poll_stride == 0 {
        return Err(SheafRepairError::InvalidBudget {
            field: "poll_stride",
        });
    }
    let scalar_slots = checked_scalar_envelope(skeleton)?;
    if scalar_slots > budget.max_scalar_slots as u128 {
        return Err(SheafRepairError::MemoryBudgetExceeded {
            required: scalar_slots,
            cap: budget.max_scalar_slots,
        });
    }
    let operator_evaluations = checked_operator_schedule(skeleton, budget.sweeps)?;
    if operator_evaluations > budget.max_operator_evaluations as u128 {
        return Err(SheafRepairError::WorkBudgetExceeded {
            required: operator_evaluations,
            cap: budget.max_operator_evaluations,
        });
    }
    let work_items = checked_hodge_work_envelope(skeleton, operator_evaluations)?;
    if work_items > budget.max_work_items as u128 {
        return Err(SheafRepairError::WorkItemBudgetExceeded {
            stage: "hodge-work-preflight",
            required: work_items,
            cap: budget.max_work_items,
        });
    }
    Ok(RepairAdmission {
        scalar_slots: usize::try_from(scalar_slots).map_err(|_| {
            SheafRepairError::BudgetArithmeticOverflow {
                stage: "scalar-envelope-publication",
            }
        })?,
        operator_evaluations: usize::try_from(operator_evaluations).map_err(|_| {
            SheafRepairError::BudgetArithmeticOverflow {
                stage: "operator-envelope-publication",
            }
        })?,
        work_items: usize::try_from(work_items).map_err(|_| {
            SheafRepairError::BudgetArithmeticOverflow {
                stage: "work-envelope-publication",
            }
        })?,
    })
}

pub(super) fn planned_cost(admission: RepairAdmission) -> Result<u64, SheafRepairError> {
    let cost = (admission.work_items as u128)
        .checked_add(admission.operator_evaluations as u128)
        .ok_or(SheafRepairError::BudgetArithmeticOverflow {
            stage: "ambient-cost-plan",
        })?;
    u64::try_from(cost).map_err(|_| SheafRepairError::BudgetArithmeticOverflow {
        stage: "ambient-cost-plan",
    })
}

fn checked_difference(
    left: &[f64],
    right: &[f64],
    stage: &'static str,
    accountant: &mut RepairAccountant<'_, '_>,
) -> Result<Vec<f64>, SheafRepairError> {
    if left.len() != right.len() {
        return Err(SheafRepairError::Skeleton(
            SheafSkeletonError::CochainLength {
                role: stage,
                expected: left.len(),
                actual: right.len(),
            },
        ));
    }
    let mut output = zeroed_output_bounded(left.len(), stage, accountant)?;
    for ((value, a), b) in output.iter_mut().zip(left).zip(right) {
        accountant.consume_item(stage)?;
        *value = a - b;
        if !value.is_finite() {
            return Err(SheafRepairError::NumericalOverflow { stage });
        }
    }
    Ok(output)
}

pub(super) fn validate_bounded_cochain(
    values: &[f64],
    expected: usize,
    role: &'static str,
    stage: &'static str,
    accountant: &mut RepairAccountant<'_, '_>,
) -> Result<(), SheafRepairError> {
    if values.len() != expected {
        return Err(SheafSkeletonError::CochainLength {
            role,
            expected,
            actual: values.len(),
        }
        .into());
    }
    for (index, value) in values.iter().enumerate() {
        accountant.consume_item(stage)?;
        if !value.is_finite() {
            return Err(SheafSkeletonError::NonFiniteCochain { role, index }.into());
        }
    }
    Ok(())
}

fn validate_cochain_length(
    values: &[f64],
    expected: usize,
    role: &'static str,
) -> Result<(), SheafRepairError> {
    if values.len() != expected {
        return Err(SheafSkeletonError::CochainLength {
            role,
            expected,
            actual: values.len(),
        }
        .into());
    }
    Ok(())
}

pub(super) fn bounded_d0(
    skeleton: &AdmittedSheafSkeleton,
    values: &[f64],
    stage: &'static str,
    accountant: &mut RepairAccountant<'_, '_>,
) -> Result<Vec<f64>, SheafRepairError> {
    validate_cochain_length(values, skeleton.n_patches, "vertex")?;
    accountant.begin_operator(stage)?;
    let mut output = zeroed_output_bounded(skeleton.edges.len(), "bounded-d0-output", accountant)?;
    for (value, &(u, v)) in output.iter_mut().zip(&skeleton.edges) {
        accountant.consume_item(stage)?;
        *value = values[v] - values[u];
        if !value.is_finite() {
            return Err(SheafRepairError::NumericalOverflow { stage });
        }
    }
    accountant.finish_operator(stage)?;
    Ok(output)
}

fn bounded_d0t(
    skeleton: &AdmittedSheafSkeleton,
    values: &[f64],
    stage: &'static str,
    accountant: &mut RepairAccountant<'_, '_>,
) -> Result<Vec<f64>, SheafRepairError> {
    validate_cochain_length(values, skeleton.edges.len(), "edge")?;
    accountant.begin_operator(stage)?;
    let mut output = zeroed_output_bounded(skeleton.n_patches, "bounded-d0t-output", accountant)?;
    for (edge, &(u, v)) in skeleton.edges.iter().enumerate() {
        accountant.consume_item(stage)?;
        output[u] -= values[edge];
        output[v] += values[edge];
        if !(output[u].is_finite() && output[v].is_finite()) {
            return Err(SheafRepairError::NumericalOverflow { stage });
        }
    }
    accountant.finish_operator(stage)?;
    Ok(output)
}

fn bounded_edge_index(
    skeleton: &AdmittedSheafSkeleton,
    target: (usize, usize),
    triangle: usize,
    stage: &'static str,
    accountant: &mut RepairAccountant<'_, '_>,
) -> Result<usize, SheafRepairError> {
    let mut lower = 0usize;
    let mut upper = skeleton.edges.len();
    while lower < upper {
        accountant.consume_item(stage)?;
        let middle = lower + (upper - lower) / 2;
        match skeleton.edges[middle].cmp(&target) {
            core::cmp::Ordering::Less => lower = middle + 1,
            core::cmp::Ordering::Equal => return Ok(middle),
            core::cmp::Ordering::Greater => upper = middle,
        }
    }
    Err(SheafSkeletonError::InvalidTriangle { index: triangle }.into())
}

fn bounded_d1(
    skeleton: &AdmittedSheafSkeleton,
    values: &[f64],
    stage: &'static str,
    accountant: &mut RepairAccountant<'_, '_>,
) -> Result<Vec<f64>, SheafRepairError> {
    validate_cochain_length(values, skeleton.edges.len(), "edge")?;
    accountant.begin_operator(stage)?;
    let mut output =
        zeroed_output_bounded(skeleton.triangles.len(), "bounded-d1-output", accountant)?;
    for (triangle, (value, &(a, b, c))) in output.iter_mut().zip(&skeleton.triangles).enumerate() {
        accountant.consume_item(stage)?;
        let eab = bounded_edge_index(skeleton, (a, b), triangle, stage, accountant)?;
        let ebc = bounded_edge_index(skeleton, (b, c), triangle, stage, accountant)?;
        let eac = bounded_edge_index(skeleton, (a, c), triangle, stage, accountant)?;
        *value = (values[eab] + values[ebc]) - values[eac];
        if !value.is_finite() {
            return Err(SheafRepairError::NumericalOverflow { stage });
        }
    }
    accountant.finish_operator(stage)?;
    Ok(output)
}

fn bounded_d1t(
    skeleton: &AdmittedSheafSkeleton,
    values: &[f64],
    stage: &'static str,
    accountant: &mut RepairAccountant<'_, '_>,
) -> Result<Vec<f64>, SheafRepairError> {
    validate_cochain_length(values, skeleton.triangles.len(), "triangle")?;
    accountant.begin_operator(stage)?;
    let mut output = zeroed_output_bounded(skeleton.edges.len(), "bounded-d1t-output", accountant)?;
    for (triangle, &(a, b, c)) in skeleton.triangles.iter().enumerate() {
        accountant.consume_item(stage)?;
        let eab = bounded_edge_index(skeleton, (a, b), triangle, stage, accountant)?;
        let ebc = bounded_edge_index(skeleton, (b, c), triangle, stage, accountant)?;
        let eac = bounded_edge_index(skeleton, (a, c), triangle, stage, accountant)?;
        output[eab] += values[triangle];
        output[ebc] += values[triangle];
        output[eac] -= values[triangle];
        if !(output[eab].is_finite() && output[ebc].is_finite() && output[eac].is_finite()) {
            return Err(SheafRepairError::NumericalOverflow { stage });
        }
    }
    accountant.finish_operator(stage)?;
    Ok(output)
}

#[derive(Clone, Copy)]
enum ProjectionKind {
    Exact,
    Coexact,
}

fn apply_projection(
    kind: ProjectionKind,
    skeleton: &AdmittedSheafSkeleton,
    values: &[f64],
    stage: &'static str,
    accountant: &mut RepairAccountant<'_, '_>,
) -> Result<Vec<f64>, SheafRepairError> {
    match kind {
        ProjectionKind::Exact => bounded_d0(skeleton, values, stage, accountant),
        ProjectionKind::Coexact => bounded_d1t(skeleton, values, stage, accountant),
    }
}

fn apply_projection_transpose(
    kind: ProjectionKind,
    skeleton: &AdmittedSheafSkeleton,
    values: &[f64],
    stage: &'static str,
    accountant: &mut RepairAccountant<'_, '_>,
) -> Result<Vec<f64>, SheafRepairError> {
    match kind {
        ProjectionKind::Exact => bounded_d0t(skeleton, values, stage, accountant),
        ProjectionKind::Coexact => bounded_d1(skeleton, values, stage, accountant),
    }
}

fn least_squares_bounded(
    skeleton: &AdmittedSheafSkeleton,
    m: &[f64],
    n_unknowns: usize,
    kind: ProjectionKind,
    pin_first: bool,
    stage: &'static str,
    accountant: &mut RepairAccountant<'_, '_>,
) -> Result<Vec<f64>, SheafRepairError> {
    let mut x = zeroed_output_bounded(n_unknowns, "least-squares-solution", accountant)?;
    let rhs = apply_projection_transpose(kind, skeleton, m, stage, accountant)?;
    let mut diag = zeroed_output_bounded(n_unknowns, "least-squares-diagonal", accountant)?;
    for (i, diagonal) in diag.iter_mut().enumerate() {
        let mut basis = zeroed_output_bounded(n_unknowns, "least-squares-basis", accountant)?;
        basis[i] = 1.0;
        let image = apply_projection(kind, skeleton, &basis, stage, accountant)?;
        *diagonal = checked_norm2(&image, "least-squares-diagonal", accountant)?;
    }
    for _ in 0..accountant.budget.sweeps {
        for i in 0..n_unknowns {
            if (pin_first && i == 0) || diag[i] <= 0.0 {
                continue;
            }
            accountant.consume_item(stage)?;
            let image = apply_projection(kind, skeleton, &x, stage, accountant)?;
            let normal_image =
                apply_projection_transpose(kind, skeleton, &image, stage, accountant)?;
            let gradient = normal_image[i] - rhs[i];
            let step = gradient / diag[i];
            let next = x[i] - step;
            if !(gradient.is_finite() && step.is_finite() && next.is_finite()) {
                return Err(SheafRepairError::NumericalOverflow { stage });
            }
            x[i] = next;
        }
        accountant.complete_sweep(stage)?;
    }
    Ok(x)
}

pub(super) fn hodge_decompose_accounted(
    skeleton: &AdmittedSheafSkeleton,
    mismatch: &[f64],
    accountant: &mut RepairAccountant<'_, '_>,
) -> Result<HodgeSplit, SheafRepairError> {
    validate_bounded_cochain(
        mismatch,
        skeleton.edges.len(),
        "edge",
        "mismatch-validation",
        accountant,
    )?;

    let potential = least_squares_bounded(
        skeleton,
        mismatch,
        skeleton.n_patches,
        ProjectionKind::Exact,
        true,
        "exact-projection",
        accountant,
    )?;
    let exact = bounded_d0(skeleton, &potential, "exact-publication", accountant)?;
    let first_residual = checked_difference(mismatch, &exact, "exact-residual", accountant)?;

    let coexact = if skeleton.triangles.is_empty() {
        zeroed_output_bounded(mismatch.len(), "empty-coexact", accountant)?
    } else {
        let triangle_potential = least_squares_bounded(
            skeleton,
            &first_residual,
            skeleton.triangles.len(),
            ProjectionKind::Coexact,
            false,
            "coexact-projection",
            accountant,
        )?;
        bounded_d1t(
            skeleton,
            &triangle_potential,
            "coexact-publication",
            accountant,
        )?
    };
    let harmonic = checked_difference(&first_residual, &coexact, "harmonic-residual", accountant)?;
    let total = checked_norm2(mismatch, "input-norm", accountant)?.max(f64::MIN_POSITIVE);
    let fractions = (
        checked_norm2(&exact, "exact-norm", accountant)? / total,
        checked_norm2(&coexact, "coexact-norm", accountant)? / total,
        checked_norm2(&harmonic, "harmonic-norm", accountant)? / total,
    );
    if [fractions.0, fractions.1, fractions.2]
        .into_iter()
        .any(|fraction| !fraction.is_finite())
    {
        return Err(SheafRepairError::NumericalOverflow {
            stage: "component-fractions",
        });
    }
    Ok(HodgeSplit {
        exact,
        potential,
        coexact,
        harmonic,
        fractions,
    })
}

/// Run the fixed-iteration Hodge-inspired diagnostic over sealed incidence
/// under one explicit resource envelope and caller cancellation context.
///
/// This function retains the same no-claim boundary as [`hodge_decompose`]: a
/// successful return proves that the declared deterministic arithmetic
/// completed within budget. It does not certify convergence, orthogonality, or
/// a topological interpretation of the remainder.
pub fn hodge_decompose_bounded(
    skeleton: &AdmittedSheafSkeleton,
    mismatch: &[f64],
    budget: SheafRepairBudget,
    cx: &Cx<'_>,
) -> Result<BoundedHodgeSplit, SheafRepairError> {
    validate_cochain_length(mismatch, skeleton.edges.len(), "edge")?;
    let admission = admit_repair_budget(skeleton, budget)?;
    let mut accountant = RepairAccountant::new(cx, budget, planned_cost(admission)?, 0, 0)?;
    accountant.checkpoint("admission")?;
    let split = hodge_decompose_accounted(skeleton, mismatch, &mut accountant)?;
    accountant.checkpoint("publication")?;
    let usage = accountant.usage(admission.scalar_slots, admission.work_items);
    Ok(BoundedHodgeSplit {
        split,
        budget,
        usage,
    })
}

fn apply_raw_projection(
    kind: ProjectionKind,
    skeleton: &SheafSkeleton,
    values: &[f64],
) -> Result<Vec<f64>, SheafRepairError> {
    match kind {
        ProjectionKind::Exact => skeleton
            .d0_validated(values)
            .map_err(SheafRepairError::from),
        ProjectionKind::Coexact => skeleton
            .d1t_validated(values)
            .map_err(SheafRepairError::from),
    }
}

fn apply_raw_projection_transpose(
    kind: ProjectionKind,
    skeleton: &SheafSkeleton,
    values: &[f64],
) -> Result<Vec<f64>, SheafRepairError> {
    match kind {
        ProjectionKind::Exact => skeleton
            .d0t_validated(values)
            .map_err(SheafRepairError::from),
        ProjectionKind::Coexact => skeleton
            .d1_validated(values)
            .map_err(SheafRepairError::from),
    }
}

fn checked_raw_norm2(values: &[f64], stage: &'static str) -> Result<f64, SheafRepairError> {
    let mut total = 0.0f64;
    for value in values {
        let square = value * value;
        total += square;
        if !(square.is_finite() && total.is_finite()) {
            return Err(SheafRepairError::NumericalOverflow { stage });
        }
    }
    Ok(total)
}

/// Least squares `min ‖m − A x‖²` via 400 deterministic Gauss–Seidel
/// sweeps on the normal equations. The raw compatibility path is capped by
/// validated skeleton cardinalities and every allocation/arithmetic step is
/// fallible; new authority paths should use the explicit bounded API.
fn least_squares_raw(
    skeleton: &SheafSkeleton,
    m: &[f64],
    n_unknowns: usize,
    kind: ProjectionKind,
    pin_first: bool,
    stage: &'static str,
) -> Result<Vec<f64>, SheafRepairError> {
    let mut x = zeroed_output(n_unknowns, "raw-least-squares-solution")?;
    let rhs = apply_raw_projection_transpose(kind, skeleton, m)?;
    let mut diag = zeroed_output(n_unknowns, "raw-least-squares-diagonal")?;
    for (index, diagonal) in diag.iter_mut().enumerate() {
        let mut basis = zeroed_output(n_unknowns, "raw-least-squares-basis")?;
        basis[index] = 1.0;
        let image = apply_raw_projection(kind, skeleton, &basis)?;
        *diagonal = checked_raw_norm2(&image, "raw-least-squares-diagonal")?;
    }
    for _ in 0..400 {
        for index in 0..n_unknowns {
            if (pin_first && index == 0) || diag[index] <= 0.0 {
                continue;
            }
            let image = apply_raw_projection(kind, skeleton, &x)?;
            let normal_image = apply_raw_projection_transpose(kind, skeleton, &image)?;
            let gradient = normal_image[index] - rhs[index];
            let step = gradient / diag[index];
            let next = x[index] - step;
            if !(gradient.is_finite() && step.is_finite() && next.is_finite()) {
                return Err(SheafRepairError::NumericalOverflow { stage });
            }
            x[index] = next;
        }
    }
    Ok(x)
}

fn hodge_decompose_raw_validated(
    skeleton: &SheafSkeleton,
    mismatch: &[f64],
) -> Result<HodgeSplit, SheafRepairError> {
    let potential = least_squares_raw(
        skeleton,
        mismatch,
        skeleton.n_patches,
        ProjectionKind::Exact,
        true,
        "raw-exact-projection",
    )?;
    let exact = skeleton.d0_validated(&potential)?;
    let mut first_residual = zeroed_output(mismatch.len(), "raw-exact-residual")?;
    for ((value, input), projected) in first_residual.iter_mut().zip(mismatch).zip(&exact) {
        *value = input - projected;
        if !value.is_finite() {
            return Err(SheafRepairError::NumericalOverflow {
                stage: "raw-exact-residual",
            });
        }
    }
    let coexact = if skeleton.triangles.is_empty() {
        zeroed_output(mismatch.len(), "raw-empty-coexact")?
    } else {
        let triangle_potential = least_squares_raw(
            skeleton,
            &first_residual,
            skeleton.triangles.len(),
            ProjectionKind::Coexact,
            false,
            "raw-coexact-projection",
        )?;
        skeleton.d1t_validated(&triangle_potential)?
    };
    let mut harmonic = zeroed_output(mismatch.len(), "raw-harmonic-residual")?;
    for ((value, residual), projected) in harmonic.iter_mut().zip(&first_residual).zip(&coexact) {
        *value = residual - projected;
        if !value.is_finite() {
            return Err(SheafRepairError::NumericalOverflow {
                stage: "raw-harmonic-residual",
            });
        }
    }
    let total = checked_raw_norm2(mismatch, "raw-input-norm")?.max(f64::MIN_POSITIVE);
    let fractions = (
        checked_raw_norm2(&exact, "raw-exact-norm")? / total,
        checked_raw_norm2(&coexact, "raw-coexact-norm")? / total,
        checked_raw_norm2(&harmonic, "raw-harmonic-norm")? / total,
    );
    if [fractions.0, fractions.1, fractions.2]
        .into_iter()
        .any(|fraction| !fraction.is_finite())
    {
        return Err(SheafRepairError::NumericalOverflow {
            stage: "raw-component-fractions",
        });
    }
    Ok(HodgeSplit {
        exact,
        potential,
        coexact,
        harmonic,
        fractions,
    })
}

/// Sequentially fit an edge cochain over a skeleton. A retained fixture checks
/// the first fit against an independent dense reference, but this fixed-count
/// solver returns no convergence or orthogonality certificate. Consumers must
/// verify residual identities such as `d0t(remainder) ≈ 0` and
/// `d1(remainder) ≈ 0` before assigning stronger meaning to a result.
///
/// # Errors
/// Returns a typed structural, cardinality, finiteness, allocation, or
/// arithmetic refusal. No partial decomposition is published on failure.
pub fn hodge_decompose(
    skeleton: &SheafSkeleton,
    mismatch: &[f64],
) -> Result<HodgeSplit, SheafRepairError> {
    validate_raw_skeleton_shape(skeleton)?;
    validate_finite_cochain(mismatch, skeleton.edges.len(), "mismatch")?;
    validate_raw_skeleton_cross_structure(skeleton)?;
    hodge_decompose_raw_validated(skeleton, mismatch)
}

/// One ranked repair proposal (the agent-facing format).
#[derive(Debug, Clone, PartialEq)]
pub struct RepairProposal {
    /// What to do, concretely.
    pub action: String,
    /// Expected post-repair worst interface mismatch. `+∞` means this proposal
    /// has no comparable constructive seam-norm prediction yet.
    pub expected_post_norm: f64,
    /// Cost estimate in seconds (router-modeled where applicable).
    pub cost_s: f64,
}

/// The repair verdict for one model.
#[derive(Debug, Clone, PartialEq)]
pub struct RepairPlan {
    /// The decomposition driving the plan.
    pub split: HodgeSplit,
    /// Ranked proposals (best first).
    pub proposals: Vec<RepairProposal>,
    /// Gauge offsets the eligible exact-component step would apply (per patch).
    pub gauge: Vec<f64>,
    /// True when the exact-component gauge step fits EVERY patch budget.
    /// This does not claim the complete repair is automatic when coexact or
    /// retained harmonic components remain.
    pub gauge_step_eligible: bool,
    /// Interfaces in the retained harmonic support with their magnitudes.
    /// This is not a graph-theoretic minimal cut-set.
    pub harmonic_support: Vec<((usize, usize), f64)>,
    /// Structured reason an optional router alternative could not be planned.
    /// `None` means no reroute was requested or a proposal was produced.
    pub reroute_error: Option<RoutePlanError>,
}

fn checked_plan_work_envelope(
    skeleton: &AdmittedSheafSkeleton,
    reroute_edges: usize,
    hodge_work_items: usize,
    action_bytes: usize,
) -> Result<u128, SheafRepairError> {
    let span = (skeleton.n_patches as u128)
        .checked_add(skeleton.edges.len() as u128)
        .and_then(|value| value.checked_add(skeleton.triangles.len() as u128))
        .and_then(|value| value.checked_add(reroute_edges as u128))
        .and_then(|value| value.checked_add(1))
        .ok_or(SheafRepairError::BudgetArithmeticOverflow {
            stage: "plan-work-span",
        })?;
    let plan_work = span
        .checked_mul(span)
        .and_then(|value| value.checked_mul(8))
        .and_then(|value| value.checked_add(128))
        .ok_or(SheafRepairError::BudgetArithmeticOverflow {
            stage: "plan-work-envelope",
        })?;
    let formatting_work = (action_bytes as u128).checked_mul(2).ok_or(
        SheafRepairError::BudgetArithmeticOverflow {
            stage: "plan-action-work-envelope",
        },
    )?;
    (hodge_work_items as u128)
        .checked_add(plan_work)
        .and_then(|value| value.checked_add(formatting_work))
        .ok_or(SheafRepairError::BudgetArithmeticOverflow {
            stage: "whole-plan-work-envelope",
        })
}

fn checked_plan_memory_envelope(
    skeleton: &AdmittedSheafSkeleton,
    budget: SheafRepairPlanBudget,
) -> Result<u128, SheafRepairError> {
    let patches = skeleton.n_patches as u128;
    let component_bytes = patches
        .checked_mul(core::mem::size_of::<usize>() as u128)
        .and_then(|value| {
            patches
                .checked_mul(4 * core::mem::size_of::<f64>() as u128)
                .and_then(|floats| value.checked_add(floats))
        })
        .ok_or(SheafRepairError::BudgetArithmeticOverflow {
            stage: "plan-component-memory",
        })?;
    let support_bytes = (budget.max_harmonic_support as u128)
        .checked_mul(core::mem::size_of::<((usize, usize), f64)>() as u128)
        .ok_or(SheafRepairError::BudgetArithmeticOverflow {
            stage: "plan-support-memory",
        })?;
    let proposal_bytes = (budget.max_proposals as u128)
        .checked_mul(core::mem::size_of::<RepairProposal>() as u128)
        .ok_or(SheafRepairError::BudgetArithmeticOverflow {
            stage: "plan-proposal-memory",
        })?;
    component_bytes
        .checked_add(support_bytes)
        .and_then(|value| value.checked_add(proposal_bytes))
        .and_then(|value| value.checked_add(budget.max_action_bytes as u128))
        .ok_or(SheafRepairError::BudgetArithmeticOverflow {
            stage: "plan-memory-envelope",
        })
}

fn plan_vec_with_capacity<T>(
    capacity: usize,
    stage: &'static str,
    accountant: &mut RepairAccountant<'_, '_>,
) -> Result<Vec<T>, SheafRepairError> {
    let bytes = capacity.checked_mul(core::mem::size_of::<T>()).ok_or(
        SheafRepairError::BudgetArithmeticOverflow {
            stage: "plan-vector-capacity",
        },
    )?;
    accountant.reserve_plan_bytes(stage, bytes)?;
    let mut output = Vec::new();
    output
        .try_reserve_exact(capacity)
        .map_err(|_| SheafSkeletonError::ResourceExhausted { stage })?;
    accountant.checkpoint(stage)?;
    Ok(output)
}

#[derive(Default)]
struct ActionByteCounter {
    bytes: usize,
}

impl core::fmt::Write for ActionByteCounter {
    fn write_str(&mut self, value: &str) -> core::fmt::Result {
        self.bytes = self
            .bytes
            .checked_add(value.len())
            .ok_or(core::fmt::Error)?;
        Ok(())
    }
}

struct AccountedActionWriter<'a, 'cx, 'clock, W: core::fmt::Write + ?Sized> {
    inner: &'a mut W,
    accountant: &'a mut RepairAccountant<'cx, 'clock>,
    stage: &'static str,
    bytes_written: usize,
    base_action_bytes: usize,
    byte_cap: usize,
    error: Option<SheafRepairError>,
}

impl<W: core::fmt::Write + ?Sized> core::fmt::Write for AccountedActionWriter<'_, '_, '_, W> {
    fn write_str(&mut self, value: &str) -> core::fmt::Result {
        for character in value.chars() {
            let next = match self.bytes_written.checked_add(character.len_utf8()) {
                Some(next) => next,
                None => {
                    self.error = Some(SheafRepairError::BudgetArithmeticOverflow {
                        stage: "action-length",
                    });
                    return Err(core::fmt::Error);
                }
            };
            if next > self.byte_cap {
                let required = match self.base_action_bytes.checked_add(next) {
                    Some(required) => required,
                    None => {
                        self.error = Some(SheafRepairError::BudgetArithmeticOverflow {
                            stage: "action-length",
                        });
                        return Err(core::fmt::Error);
                    }
                };
                let cap = match self.base_action_bytes.checked_add(self.byte_cap) {
                    Some(cap) => cap,
                    None => {
                        self.error = Some(SheafRepairError::BudgetArithmeticOverflow {
                            stage: "action-length",
                        });
                        return Err(core::fmt::Error);
                    }
                };
                self.error = Some(SheafRepairError::OutputBudgetExceeded {
                    resource: "action-bytes",
                    required: required as u128,
                    cap,
                });
                return Err(core::fmt::Error);
            }
            for _ in 0..character.len_utf8() {
                if let Err(error) = self.accountant.consume_item(self.stage) {
                    self.error = Some(error);
                    return Err(core::fmt::Error);
                }
            }
            self.inner.write_char(character)?;
            self.bytes_written = next;
        }
        Ok(())
    }
}

fn try_build_action(
    stage: &'static str,
    accountant: &mut RepairAccountant<'_, '_>,
    render: impl Fn(&mut dyn core::fmt::Write) -> core::fmt::Result,
) -> Result<String, SheafRepairError> {
    accountant.checkpoint(stage)?;
    let mut counter = ActionByteCounter::default();
    let action_base = accountant.action_bytes;
    let action_remaining = accountant.action_bytes_cap.checked_sub(action_base).ok_or(
        SheafRepairError::BudgetArithmeticOverflow {
            stage: "action-remaining-capacity",
        },
    )?;
    let count_result;
    let count_error;
    {
        let mut writer = AccountedActionWriter {
            inner: &mut counter,
            accountant,
            stage,
            bytes_written: 0,
            base_action_bytes: action_base,
            byte_cap: action_remaining,
            error: None,
        };
        count_result = render(&mut writer);
        count_error = writer.error;
    }
    if let Some(error) = count_error {
        return Err(error);
    }
    count_result.map_err(|_| SheafRepairError::BudgetArithmeticOverflow {
        stage: "action-length",
    })?;
    accountant.retain_action_bytes(stage, counter.bytes)?;
    let mut action = String::new();
    action
        .try_reserve_exact(counter.bytes)
        .map_err(|_| SheafSkeletonError::ResourceExhausted { stage })?;
    let render_result;
    let render_error;
    {
        let action_base = accountant.action_bytes.checked_sub(counter.bytes).ok_or(
            SheafRepairError::BudgetArithmeticOverflow {
                stage: "action-render-base",
            },
        )?;
        let mut writer = AccountedActionWriter {
            inner: &mut action,
            accountant,
            stage,
            bytes_written: 0,
            base_action_bytes: action_base,
            byte_cap: counter.bytes,
            error: None,
        };
        render_result = render(&mut writer);
        render_error = writer.error;
    }
    if let Some(error) = render_error {
        return Err(error);
    }
    render_result.map_err(|_| SheafRepairError::NumericalOverflow { stage })?;
    if action.len() != counter.bytes {
        return Err(SheafRepairError::NumericalOverflow { stage });
    }
    accountant.checkpoint(stage)?;
    Ok(action)
}

fn bounded_find_root(
    parent: &[usize],
    mut patch: usize,
    accountant: &mut RepairAccountant<'_, '_>,
) -> Result<usize, SheafRepairError> {
    for _ in 0..=parent.len() {
        accountant.consume_item("gauge-component-find")?;
        let next = parent[patch];
        if next == patch {
            return Ok(patch);
        }
        patch = next;
    }
    Err(SheafRepairError::NumericalOverflow {
        stage: "gauge-component-find",
    })
}

fn bounded_gauge_representative(
    skeleton: &AdmittedSheafSkeleton,
    potential: &[f64],
    budgets: &[f64],
    accountant: &mut RepairAccountant<'_, '_>,
) -> Result<(Vec<f64>, bool), SheafRepairError> {
    let mut parent =
        plan_vec_with_capacity::<usize>(skeleton.n_patches, "gauge-component-parent", accountant)?;
    for patch in 0..skeleton.n_patches {
        accountant.consume_item("gauge-component-parent")?;
        parent.push(patch);
    }
    for &(u, v) in &skeleton.edges {
        accountant.consume_item("gauge-component-union")?;
        let u_root = bounded_find_root(&parent, u, accountant)?;
        let v_root = bounded_find_root(&parent, v, accountant)?;
        if u_root != v_root {
            let (low, high) = if u_root < v_root {
                (u_root, v_root)
            } else {
                (v_root, u_root)
            };
            parent[high] = low;
        }
    }

    let mut lower =
        plan_vec_with_capacity::<f64>(skeleton.n_patches, "gauge-component-lower", accountant)?;
    let mut upper =
        plan_vec_with_capacity::<f64>(skeleton.n_patches, "gauge-component-upper", accountant)?;
    for _ in 0..skeleton.n_patches {
        accountant.consume_item("gauge-component-bounds")?;
        lower.push(f64::NEG_INFINITY);
        upper.push(f64::INFINITY);
    }
    for patch in 0..skeleton.n_patches {
        accountant.consume_item("gauge-component-bounds")?;
        let root = bounded_find_root(&parent, patch, accountant)?;
        lower[root] = lower[root].max(-budgets[patch] - potential[patch]);
        upper[root] = upper[root].min(budgets[patch] - potential[patch]);
        if !(lower[root].is_finite() && upper[root].is_finite()) {
            return Err(SheafRepairError::NumericalOverflow {
                stage: "gauge-component-bounds",
            });
        }
    }

    let mut shifts =
        plan_vec_with_capacity::<f64>(skeleton.n_patches, "gauge-component-shifts", accountant)?;
    let mut feasible = true;
    for patch in 0..skeleton.n_patches {
        accountant.consume_item("gauge-component-shifts")?;
        let shift = if parent[patch] == patch {
            if lower[patch] > upper[patch] {
                feasible = false;
                0.0
            } else {
                f64::midpoint(lower[patch], upper[patch])
            }
        } else {
            0.0
        };
        if !shift.is_finite() {
            return Err(SheafRepairError::NumericalOverflow {
                stage: "gauge-component-shifts",
            });
        }
        shifts.push(shift);
    }

    if feasible {
        for patch in 0..skeleton.n_patches {
            accountant.consume_item("gauge-representative-check")?;
            let root = bounded_find_root(&parent, patch, accountant)?;
            let shifted = potential[patch] + shifts[root];
            if !shifted.is_finite() || shifted.abs() > budgets[patch] {
                feasible = false;
                break;
            }
        }
    }
    let mut gauge =
        plan_vec_with_capacity::<f64>(skeleton.n_patches, "gauge-representative", accountant)?;
    for patch in 0..skeleton.n_patches {
        accountant.consume_item("gauge-representative")?;
        let shifted = if feasible {
            let root = bounded_find_root(&parent, patch, accountant)?;
            potential[patch] + shifts[root]
        } else {
            potential[patch]
        };
        if !shifted.is_finite() {
            return Err(SheafRepairError::NumericalOverflow {
                stage: "gauge-representative",
            });
        }
        gauge.push(shifted);
    }
    Ok((gauge, feasible))
}

fn bounded_expected_after_gauge(
    skeleton: &AdmittedSheafSkeleton,
    mismatch: &[f64],
    gauge: &[f64],
    accountant: &mut RepairAccountant<'_, '_>,
) -> Result<f64, SheafRepairError> {
    let mut expected = 0.0f64;
    for (edge, &(u, v)) in skeleton.edges.iter().enumerate() {
        accountant.consume_item("gauge-post-norm")?;
        let correction = gauge[v] - gauge[u];
        let residual = mismatch[edge] - correction;
        if !(correction.is_finite() && residual.is_finite()) {
            return Err(SheafRepairError::NumericalOverflow {
                stage: "gauge-post-norm",
            });
        }
        expected = expected.max(residual.abs());
    }
    Ok(expected)
}

fn bounded_gauge_proposal(
    gauge: &[f64],
    eligible: bool,
    expected: f64,
    accountant: &mut RepairAccountant<'_, '_>,
) -> Result<RepairProposal, SheafRepairError> {
    let mut worst = 0usize;
    for index in 0..gauge.len() {
        accountant.consume_item("gauge-proposal-rank")?;
        if gauge[index].abs().total_cmp(&gauge[worst].abs()).is_ge() {
            worst = index;
        }
    }
    for _ in gauge {
        accountant.consume_item("gauge-proposal-format")?;
    }
    let action = try_build_action("gauge-proposal-action", accountant, |writer| {
        write!(
            writer,
            "project patch P{worst} gauge by {:+.3e} (exact-component section projection; offsets per patch: [",
            gauge[worst]
        )?;
        for (index, offset) in gauge.iter().enumerate() {
            if index > 0 {
                writer.write_str(", ")?;
            }
            write!(writer, "{offset:+.3e}")?;
        }
        writer.write_str("])")?;
        if !eligible {
            writer.write_str(" — EXCEEDS a patch budget; needs acceptance")?;
        }
        Ok(())
    })?;
    Ok(RepairProposal {
        action,
        expected_post_norm: expected,
        cost_s: 0.001,
    })
}

fn bounded_coexact_proposal(
    skeleton: &AdmittedSheafSkeleton,
    mismatch: &[f64],
    accountant: &mut RepairAccountant<'_, '_>,
) -> Result<RepairProposal, SheafRepairError> {
    let mut worst: Option<((usize, usize, usize), f64)> = None;
    for (triangle, &(a, b, c)) in skeleton.triangles.iter().enumerate() {
        accountant.consume_item("coexact-proposal-localize")?;
        let eab = bounded_edge_index(
            skeleton,
            (a, b),
            triangle,
            "coexact-proposal-localize",
            accountant,
        )?;
        let ebc = bounded_edge_index(
            skeleton,
            (b, c),
            triangle,
            "coexact-proposal-localize",
            accountant,
        )?;
        let eac = bounded_edge_index(
            skeleton,
            (a, c),
            triangle,
            "coexact-proposal-localize",
            accountant,
        )?;
        let circulation = (mismatch[eab] + mismatch[ebc]) - mismatch[eac];
        if !circulation.is_finite() {
            return Err(SheafRepairError::NumericalOverflow {
                stage: "coexact-proposal-localize",
            });
        }
        let replace = match worst.as_ref() {
            None => true,
            Some((_, value)) => circulation.abs().total_cmp(&value.abs()).is_ge(),
        };
        if replace {
            worst = Some(((a, b, c), circulation));
        }
    }
    let worst_triangle = worst.map(|(triangle, _)| triangle);
    let action = try_build_action("coexact-proposal-action", accountant, |writer| {
        write!(
            writer,
            "coexact circulation candidate around retained triangle {worst_triangle:?}: inspect chart/model/junction/sampling evidence and converter orientation/trace conventions; algebra alone does not assign cause"
        )
    })?;
    Ok(RepairProposal {
        action,
        expected_post_norm: f64::INFINITY,
        cost_s: 0.0,
    })
}

fn bounded_harmonic_support(
    skeleton: &AdmittedSheafSkeleton,
    split: &HodgeSplit,
    cap: usize,
    accountant: &mut RepairAccountant<'_, '_>,
) -> Result<Vec<((usize, usize), f64)>, SheafRepairError> {
    if split.fractions.2 <= COMPONENT_FLOOR {
        return plan_vec_with_capacity(0, "harmonic-support", accountant);
    }
    let mut scale = 0.0f64;
    for value in &split.harmonic {
        accountant.consume_item("harmonic-support-scale")?;
        scale = scale.max(value.abs());
    }
    let floor = scale * COMPONENT_FLOOR.sqrt();
    let mut required = 0usize;
    for value in &split.harmonic {
        accountant.consume_item("harmonic-support-count")?;
        if value.abs() > floor {
            required =
                required
                    .checked_add(1)
                    .ok_or(SheafRepairError::BudgetArithmeticOverflow {
                        stage: "harmonic-support-count",
                    })?;
        }
    }
    if required > cap {
        return Err(SheafRepairError::OutputBudgetExceeded {
            resource: "harmonic-support",
            required: required as u128,
            cap,
        });
    }
    let mut support = plan_vec_with_capacity(required, "harmonic-support", accountant)?;
    for (&edge, &value) in skeleton.edges.iter().zip(&split.harmonic) {
        accountant.consume_item("harmonic-support-retain")?;
        if value.abs() > floor {
            support.push((edge, value.abs()));
        }
    }
    for index in 1..support.len() {
        let mut cursor = index;
        while cursor > 0 {
            accountant.consume_item("harmonic-support-sort")?;
            if support[cursor - 1].1.total_cmp(&support[cursor].1).is_ge() {
                break;
            }
            support.swap(cursor - 1, cursor);
            cursor -= 1;
        }
    }
    Ok(support)
}

fn bounded_harmonic_proposal(
    support: &[((usize, usize), f64)],
    accountant: &mut RepairAccountant<'_, '_>,
) -> Result<RepairProposal, SheafRepairError> {
    for _ in support {
        accountant.consume_item("harmonic-proposal-format")?;
    }
    let action = try_build_action("harmonic-proposal-action", accountant, |writer| {
        writer.write_str(
            "retained harmonic remainder after deterministic gauge projection; no generic exactness or topology claim; inspect interface support [",
        )?;
        for (index, ((u, v), _)) in support.iter().enumerate() {
            if index > 0 {
                writer.write_str(", ")?;
            }
            write!(writer, "({u}, {v})")?;
        }
        writer.write_str("]")
    })?;
    Ok(RepairProposal {
        action,
        expected_post_norm: f64::INFINITY,
        cost_s: f64::INFINITY,
    })
}

fn bounded_reroute_proposal(
    route: &RoutePlan,
    accountant: &mut RepairAccountant<'_, '_>,
) -> Result<RepairProposal, SheafRepairError> {
    for _ in route.edges() {
        accountant.consume_item("reroute-proposal-format")?;
    }
    let request = route.request();
    let action = try_build_action("reroute-proposal-action", accountant, |writer| {
        write!(
            writer,
            "reroute worst patch {} -> {} via [",
            request.from, request.to
        )?;
        for (index, edge) in route.edges().iter().enumerate() {
            if index > 0 {
                writer.write_str(", ")?;
            }
            writer.write_str(edge)?;
        }
        writer.write_str("] (router-planned alternative chart)")
    })?;
    Ok(RepairProposal {
        action,
        expected_post_norm: f64::INFINITY,
        cost_s: route.predicted_cost_s(),
    })
}

fn push_bounded_proposal(
    proposals: &mut Vec<RepairProposal>,
    proposal: RepairProposal,
    cap: usize,
    accountant: &mut RepairAccountant<'_, '_>,
) -> Result<(), SheafRepairError> {
    accountant.consume_item("proposal-retain")?;
    let required =
        proposals
            .len()
            .checked_add(1)
            .ok_or(SheafRepairError::BudgetArithmeticOverflow {
                stage: "proposal-count",
            })?;
    if required > cap {
        return Err(SheafRepairError::OutputBudgetExceeded {
            resource: "proposals",
            required: required as u128,
            cap,
        });
    }
    proposals.push(proposal);
    Ok(())
}

fn sort_bounded_proposals(
    proposals: &mut [RepairProposal],
    accountant: &mut RepairAccountant<'_, '_>,
) -> Result<(), SheafRepairError> {
    for index in 1..proposals.len() {
        let mut cursor = index;
        while cursor > 0 {
            accountant.consume_item("proposal-sort")?;
            let ordering = proposals[cursor]
                .expected_post_norm
                .total_cmp(&proposals[cursor - 1].expected_post_norm)
                .then(
                    proposals[cursor]
                        .cost_s
                        .total_cmp(&proposals[cursor - 1].cost_s),
                );
            if !ordering.is_lt() {
                break;
            }
            proposals.swap(cursor - 1, cursor);
            cursor -= 1;
        }
    }
    Ok(())
}

/// Build a complete repair plan over admitted incidence under one explicit
/// decomposition, work, memory, output, deadline, and cancellation envelope.
///
/// The optional route is already admitted by the Rep Router; this function
/// never performs an unbounded live graph search. A refusal publishes no
/// partial plan. The returned proposals retain the same diagnostic no-claim
/// boundaries as [`plan_repair`].
pub fn plan_repair_bounded(
    skeleton: &AdmittedSheafSkeleton,
    mismatch: &[f64],
    gauge_budgets: &[f64],
    reroute: Option<&RoutePlan>,
    budget: SheafRepairPlanBudget,
    cx: &Cx<'_>,
) -> Result<BoundedRepairPlan, SheafRepairError> {
    validate_cochain_length(mismatch, skeleton.edges.len(), "mismatch")?;
    if gauge_budgets.len() != skeleton.n_patches {
        return Err(SheafSkeletonError::CochainLength {
            role: "gauge-budget",
            expected: skeleton.n_patches,
            actual: gauge_budgets.len(),
        }
        .into());
    }
    let hodge_admission = admit_repair_budget(skeleton, budget.repair)?;
    let work_items = checked_plan_work_envelope(
        skeleton,
        reroute.map_or(0, |route| route.edges().len()),
        hodge_admission.work_items,
        budget.max_action_bytes,
    )?;
    if work_items > budget.repair.max_work_items as u128 {
        return Err(SheafRepairError::WorkItemBudgetExceeded {
            stage: "plan-work-preflight",
            required: work_items,
            cap: budget.repair.max_work_items,
        });
    }
    let plan_memory_envelope = checked_plan_memory_envelope(skeleton, budget)?;
    if plan_memory_envelope > budget.max_plan_bytes as u128 {
        return Err(SheafRepairError::PlanMemoryBudgetExceeded {
            required: plan_memory_envelope,
            cap: budget.max_plan_bytes,
        });
    }
    let admission = RepairAdmission {
        scalar_slots: hodge_admission.scalar_slots,
        operator_evaluations: hodge_admission.operator_evaluations,
        work_items: usize::try_from(work_items).map_err(|_| {
            SheafRepairError::BudgetArithmeticOverflow {
                stage: "plan-work-publication",
            }
        })?,
    };
    let plan_memory_envelope = usize::try_from(plan_memory_envelope).map_err(|_| {
        SheafRepairError::BudgetArithmeticOverflow {
            stage: "plan-memory-publication",
        }
    })?;
    let mut accountant = RepairAccountant::new(
        cx,
        budget.repair,
        planned_cost(admission)?,
        budget.max_plan_bytes,
        budget.max_action_bytes,
    )?;
    accountant.checkpoint("plan-admission")?;
    validate_bounded_cochain(
        mismatch,
        skeleton.edges.len(),
        "mismatch",
        "mismatch-validation",
        &mut accountant,
    )?;
    for (index, gauge_budget) in gauge_budgets.iter().enumerate() {
        accountant.consume_item("gauge-budget-validation")?;
        if !gauge_budget.is_finite() || *gauge_budget < 0.0 {
            return Err(SheafRepairError::InvalidGaugeBudget { index });
        }
    }

    let split = hodge_decompose_accounted(skeleton, mismatch, &mut accountant)?;
    let (gauge, gauge_feasible) =
        bounded_gauge_representative(skeleton, &split.potential, gauge_budgets, &mut accountant)?;
    let expected_after_gauge =
        bounded_expected_after_gauge(skeleton, mismatch, &gauge, &mut accountant)?;
    let gauge_step_eligible = split.fractions.0 > COMPONENT_FLOOR && gauge_feasible;
    let harmonic_support = bounded_harmonic_support(
        skeleton,
        &split,
        budget.max_harmonic_support,
        &mut accountant,
    )?;
    let mut required_proposals = 0usize;
    for needed in [
        split.fractions.0 > COMPONENT_FLOOR,
        split.fractions.1 > COMPONENT_FLOOR,
        !harmonic_support.is_empty(),
        reroute.is_some(),
    ] {
        accountant.consume_item("proposal-count")?;
        if needed {
            required_proposals = required_proposals.checked_add(1).ok_or(
                SheafRepairError::BudgetArithmeticOverflow {
                    stage: "proposal-count",
                },
            )?;
        }
    }
    if required_proposals > budget.max_proposals {
        return Err(SheafRepairError::OutputBudgetExceeded {
            resource: "proposals",
            required: required_proposals as u128,
            cap: budget.max_proposals,
        });
    }
    let mut proposals = plan_vec_with_capacity::<RepairProposal>(
        required_proposals,
        "repair-proposals",
        &mut accountant,
    )?;
    if split.fractions.0 > COMPONENT_FLOOR {
        let proposal = bounded_gauge_proposal(
            &gauge,
            gauge_step_eligible,
            expected_after_gauge,
            &mut accountant,
        )?;
        push_bounded_proposal(
            &mut proposals,
            proposal,
            budget.max_proposals,
            &mut accountant,
        )?;
    }
    if split.fractions.1 > COMPONENT_FLOOR {
        let proposal = bounded_coexact_proposal(skeleton, mismatch, &mut accountant)?;
        push_bounded_proposal(
            &mut proposals,
            proposal,
            budget.max_proposals,
            &mut accountant,
        )?;
    }
    if !harmonic_support.is_empty() {
        let proposal = bounded_harmonic_proposal(&harmonic_support, &mut accountant)?;
        push_bounded_proposal(
            &mut proposals,
            proposal,
            budget.max_proposals,
            &mut accountant,
        )?;
    }
    if let Some(route) = reroute {
        let proposal = bounded_reroute_proposal(route, &mut accountant)?;
        push_bounded_proposal(
            &mut proposals,
            proposal,
            budget.max_proposals,
            &mut accountant,
        )?;
    }
    sort_bounded_proposals(&mut proposals, &mut accountant)?;
    accountant.checkpoint("plan-publication")?;
    let repair_usage = accountant.usage(admission.scalar_slots, admission.work_items);
    let usage = SheafRepairPlanUsage {
        repair: repair_usage,
        plan_memory_envelope,
        reserved_plan_bytes: accountant.reserved_plan_bytes,
        action_bytes: accountant.action_bytes,
        proposals: proposals.len(),
        harmonic_support: harmonic_support.len(),
    };
    Ok(BoundedRepairPlan {
        plan: RepairPlan {
            split,
            proposals,
            gauge,
            gauge_step_eligible,
            harmonic_support,
            reroute_error: None,
        },
        budget,
        usage,
    })
}

/// Threshold below which a component is treated as absent (fractions).
pub const COMPONENT_FLOOR: f64 = 1e-6;

/// Choose a deterministic maximum-slack midpoint of the feasible constant-shift
/// interval independently on each connected component (or its finite boundary
/// when the interval is half-unbounded). Adding such a constant leaves `δ⁰c`
/// unchanged mathematically; the returned gauge is the representative that the
/// planner will actually apply.
fn gauge_representative_within_budgets(
    skeleton: &SheafSkeleton,
    potential: &[f64],
    budgets: &[f64],
) -> Option<Vec<f64>> {
    if potential.len() != skeleton.n_patches
        || budgets.len() != skeleton.n_patches
        || potential.iter().any(|value| !value.is_finite())
        || budgets
            .iter()
            .any(|budget| budget.is_nan() || *budget < 0.0)
    {
        return None;
    }

    let mut adjacency = vec![Vec::new(); skeleton.n_patches];
    for &(u, v) in &skeleton.edges {
        adjacency[u].push(v);
        adjacency[v].push(u);
    }

    let mut gauge = potential.to_vec();
    let mut seen = vec![false; skeleton.n_patches];
    for root in 0..skeleton.n_patches {
        if seen[root] {
            continue;
        }
        seen[root] = true;
        let mut stack = vec![root];
        let mut component = Vec::new();
        while let Some(patch) = stack.pop() {
            component.push(patch);
            for &neighbor in &adjacency[patch] {
                if !seen[neighbor] {
                    seen[neighbor] = true;
                    stack.push(neighbor);
                }
            }
        }

        let mut lower = f64::NEG_INFINITY;
        let mut upper = f64::INFINITY;
        for &patch in &component {
            let budget = budgets[patch];
            if budget.is_finite() {
                lower = lower.max(-budget - potential[patch]);
                upper = upper.min(budget - potential[patch]);
            }
        }
        if lower > upper {
            return None;
        }
        let shift = match (lower.is_finite(), upper.is_finite()) {
            (true, true) => f64::midpoint(lower, upper),
            (true, false) => lower,
            (false, true) => upper,
            (false, false) => 0.0,
        };
        if !shift.is_finite() {
            return None;
        }
        for patch in component {
            let shifted = potential[patch] + shift;
            if !shifted.is_finite() || shifted.abs() > budgets[patch] {
                return None;
            }
            gauge[patch] = shifted;
        }
    }
    Some(gauge)
}

/// Build the repair plan: decompose, interpret, rank. `budgets` is each
/// patch's declared error budget — the exact-component gauge repair is
/// only auto-appliable when |offset| stays within it for EVERY patch
/// (a repair must never silently distort geometry beyond budget).
/// `reroute` optionally consults the Rep Router for a conversion-based
/// alternative for the worst-offending patch.
///
/// Shape, topology, and scalar validation completes before decomposition or
/// proposal allocation. A refusal returns no partial plan.
///
/// # Errors
/// Returns [`SheafRepairError`] for malformed raw incidence, wrong cochain or
/// budget cardinality, non-finite mismatch/gauge budgets, or finite arithmetic
/// overflow.
#[must_use]
pub fn plan_repair(
    skeleton: &SheafSkeleton,
    mismatch: &[f64],
    budgets: &[f64],
    reroute: Option<(&Router, &dyn CostOracle, &RouteRequest)>,
) -> Result<RepairPlan, SheafRepairError> {
    validate_raw_skeleton_shape(skeleton)?;
    validate_finite_cochain(mismatch, skeleton.edges.len(), "mismatch")?;
    if budgets.len() != skeleton.n_patches {
        return Err(SheafSkeletonError::CochainLength {
            role: "gauge-budget",
            expected: skeleton.n_patches,
            actual: budgets.len(),
        }
        .into());
    }
    if let Some(index) = budgets
        .iter()
        .position(|budget| !budget.is_finite() || *budget < 0.0)
    {
        return Err(SheafRepairError::InvalidGaugeBudget { index });
    }
    validate_raw_skeleton_cross_structure(skeleton)?;

    let split = hodge_decompose_raw_validated(skeleton, mismatch)?;
    if split
        .exact
        .iter()
        .chain(&split.potential)
        .chain(&split.coexact)
        .chain(&split.harmonic)
        .any(|value| !value.is_finite())
        || [split.fractions.0, split.fractions.1, split.fractions.2]
            .into_iter()
            .any(|value| !value.is_finite())
    {
        return Err(SheafRepairError::NumericalOverflow {
            stage: "repair-decomposition",
        });
    }
    let feasible_gauge = gauge_representative_within_budgets(skeleton, &split.potential, budgets);
    let gauge_step_is_feasible = feasible_gauge.is_some();
    let gauge = feasible_gauge.unwrap_or_else(|| split.potential.clone());
    let residual_after_exact = try_apply_gauge(skeleton, mismatch, &gauge)?;
    let expected_after_gauge = residual_after_exact
        .iter()
        .fold(0.0f64, |a, &b| a.max(b.abs()));
    let gauge_step_eligible = split.fractions.0 > COMPONENT_FLOOR && gauge_step_is_feasible;
    let mut proposals = Vec::new();
    if split.fractions.0 > COMPONENT_FLOOR {
        proposals.push(gauge_proposal(
            &gauge,
            gauge_step_eligible,
            expected_after_gauge,
        ));
    }
    if split.fractions.1 > COMPONENT_FLOOR {
        proposals.push(coexact_proposal(skeleton, mismatch)?);
    }
    // First require the whole retained component to be significant relative
    // to the input mismatch. Otherwise scaling a localization threshold by the
    // remainder's own maximum guarantees that even roundoff residue promotes
    // itself into scary-looking support and a +inf proposal. Once admitted,
    // use a within-component relative amplitude floor to localize it without a
    // unit-dependent absolute threshold. The raw split always retains the
    // sub-floor remainder for diagnosis.
    let mut harmonic_support: Vec<((usize, usize), f64)> = if split.fractions.2 > COMPONENT_FLOOR {
        let harmonic_scale = split
            .harmonic
            .iter()
            .fold(0.0f64, |scale, value| scale.max(value.abs()));
        let support_floor = harmonic_scale * COMPONENT_FLOOR.sqrt();
        skeleton
            .edges
            .iter()
            .zip(&split.harmonic)
            .filter(|(_, h)| h.abs() > support_floor)
            .map(|(&e, &h)| (e, h.abs()))
            .collect()
    } else {
        Vec::new()
    };
    harmonic_support.sort_by(|a, b| b.1.total_cmp(&a.1));
    if !harmonic_support.is_empty() {
        proposals.push(RepairProposal {
            action: format!(
                "retained harmonic remainder after deterministic gauge projection; no \
                 generic exactness or topology claim; inspect interface support {:?}",
                harmonic_support.iter().map(|(e, _)| *e).collect::<Vec<_>>()
            ),
            expected_post_norm: f64::INFINITY,
            cost_s: f64::INFINITY,
        });
    }
    let mut reroute_error = None;
    if let Some((router, oracle, req)) = reroute {
        match router.plan(req, oracle) {
            Ok(route) => proposals.push(RepairProposal {
                action: format!(
                    "reroute worst patch {} -> {} via [{}] (router-planned alternative chart)",
                    req.from,
                    req.to,
                    route.edges().join(", ")
                ),
                expected_post_norm: f64::INFINITY,
                cost_s: route.predicted_cost_s(),
            }),
            Err(error) => reroute_error = Some(error),
        }
    }
    proposals.sort_by(|a, b| {
        a.expected_post_norm
            .total_cmp(&b.expected_post_norm)
            .then(a.cost_s.total_cmp(&b.cost_s))
    });
    Ok(RepairPlan {
        gauge,
        split,
        proposals,
        gauge_step_eligible,
        harmonic_support,
        reroute_error,
    })
}

/// The exact-component proposal: the concrete per-patch gauge projection.
fn gauge_proposal(gauge: &[f64], gauge_step_eligible: bool, expected: f64) -> RepairProposal {
    let worst = gauge
        .iter()
        .enumerate()
        .max_by(|a, b| a.1.abs().total_cmp(&b.1.abs()))
        .map_or(0, |(i, _)| i);
    let mut action = format!(
        "project patch P{worst} gauge by {:+.3e} (exact-component section \
         projection; offsets per patch: [",
        gauge[worst]
    );
    for (i, off) in gauge.iter().enumerate() {
        if i > 0 {
            action.push_str(", ");
        }
        let _ = write!(action, "{off:+.3e}");
    }
    action.push_str("])");
    if !gauge_step_eligible {
        let _ = write!(action, " — EXCEEDS a patch budget; needs acceptance");
    }
    RepairProposal {
        action,
        expected_post_norm: expected,
        cost_s: 0.001, // local gauge arithmetic
    }
}

/// The coexact-component proposal: a non-causal diagnostic localized to the
/// retained triangle with the largest circulation residual.
fn coexact_proposal(
    skeleton: &SheafSkeleton,
    mismatch: &[f64],
) -> Result<RepairProposal, SheafRepairError> {
    let d1m = skeleton.d1_validated(mismatch)?;
    let worst_tri = skeleton
        .triangles
        .iter()
        .enumerate()
        .max_by(|a, b| d1m[a.0].abs().total_cmp(&d1m[b.0].abs()))
        .map(|(_, t)| *t);
    Ok(RepairProposal {
        action: format!(
            "coexact circulation candidate around retained triangle {worst_tri:?}: inspect \
             chart/model/junction/sampling evidence and converter orientation/trace \
             conventions; algebra alone does not assign cause"
        ),
        expected_post_norm: f64::INFINITY,
        cost_s: 0.0,
    })
}

/// Apply one algebraic gauge correction to an edge cochain:
/// `m ← m − δ⁰c`. Re-planning a converged repaired model can yield a zero
/// follow-up gauge; applying the same nonzero gauge twice is not idempotent.
/// This does not mutate or re-evaluate any source chart.
///
/// # Errors
/// Returns [`SheafRepairError`] before allocation for malformed incidence,
/// wrong cochain cardinality, or non-finite input, and during construction if
/// output reservation or finite subtraction fails.
#[must_use]
pub fn try_apply_gauge(
    skeleton: &SheafSkeleton,
    mismatch: &[f64],
    gauge: &[f64],
) -> Result<Vec<f64>, SheafRepairError> {
    validate_raw_skeleton_shape(skeleton)?;
    validate_finite_cochain(mismatch, skeleton.edges.len(), "mismatch")?;
    validate_finite_cochain(gauge, skeleton.n_patches, "gauge")?;
    validate_raw_skeleton_cross_structure(skeleton)?;

    let mut repaired = zeroed_output(skeleton.edges.len(), "apply-gauge-output")?;
    for (edge, (value, &(u, v))) in repaired.iter_mut().zip(&skeleton.edges).enumerate() {
        let correction = gauge[v] - gauge[u];
        *value = mismatch[edge] - correction;
        if !(correction.is_finite() && value.is_finite()) {
            return Err(SheafRepairError::NumericalOverflow {
                stage: "apply-gauge",
            });
        }
    }
    Ok(repaired)
}

/// Typed compatibility name for [`try_apply_gauge`]. Refusals are never
/// collapsed into a sentinel cochain.
///
/// # Errors
/// Returns the same typed structural, cardinality, finiteness, allocation, or
/// arithmetic refusal as [`try_apply_gauge`].
pub fn apply_gauge(
    skeleton: &SheafSkeleton,
    mismatch: &[f64],
    gauge: &[f64],
) -> Result<Vec<f64>, SheafRepairError> {
    try_apply_gauge(skeleton, mismatch, gauge)
}
