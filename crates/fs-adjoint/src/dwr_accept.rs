//! DWR GOAL-ORIENTED ACCEPT TEST (addendum Proposal 9, bead lmp4.4;
//! [F] — behind the `dwr-accept` feature): dual-weighted-residual
//! estimates target the QUERY's actual quantity of interest, which
//! sharpens the accept test enormously — but DWR constants are NOT
//! guaranteed, so a DWR-only accept carries ESTIMATED color. Promotion
//! to VERIFIED additionally requires a typed proof that the supplied dual is
//! the exact dual of THIS query. The current v0 API cannot express that proof,
//! so even an independently reverified Cauchy–Schwarz energy-product bracket
//! remains an ESTIMATED diagnostic and never promotes or vetoes acceptance.
//! False certification is worse than temporarily losing the promotion path.

use fs_evidence::Color;
use fs_verify::estimator::{EstimatorFamily, verify};
use fs_verify::fem1d::{MAX_FEM1D_POLY_COEFFICIENTS, MmsProblem, gauss5};
use fs_verify::interval::up;

/// Maximum mesh nodes admitted by the in-process rigorous bracket verifier.
pub const MAX_BRACKET_MESH_NODES: usize = 1_000_000;
/// Maximum coarse mesh/candidate nodes admitted by the DWR execution path.
pub const MAX_DWR_MESH_NODES: usize = 1_000_000;
/// Maximum manufactured-solution polynomial coefficients admitted by DWR.
pub const MAX_DWR_POLY_COEFFICIENTS: usize = MAX_FEM1D_POLY_COEFFICIENTS;
/// Maximum conservative scalar work units admitted by one DWR execution.
pub const MAX_DWR_WORK_UNITS: usize = 100_000_000;
const MAX_DWR_REFINED_NODES: usize = MAX_DWR_MESH_NODES * 2 - 1;

/// A QoI query: what the caller actually asked.
#[derive(Debug, Clone, PartialEq)]
pub struct DwrQuery {
    /// The quantity of interest (provenance label).
    pub qoi: String,
    /// The tolerance the answer must meet.
    pub tolerance: f64,
}

/// An independently reverified primal/dual energy-product diagnostic.
///
/// The fields are sealed: safe downstream code cannot assert that an arbitrary
/// number is certified. A bracket can only be created by
/// [`Bracket::cauchy_schwarz`], which reruns the equilibrated-flux verifier on
/// the exact problem/candidate pairs. It is not yet a QoI-error certificate:
/// the v0 query type does not bind the dual problem to the requested functional.
///
/// ```compile_fail
/// use fs_adjoint::dwr_accept::Bracket;
///
/// let forged = Bracket {
///     bound: 0.0,
///     source: "caller assertion".to_string(),
/// };
/// ```
#[derive(Debug, Clone, PartialEq)]
pub struct Bracket {
    /// The outward-rounded product of the two energy-error upper bounds.
    bound: f64,
    /// Where the bound came from (audit trail).
    source: String,
}

/// Why a rigorous bracket could not be issued.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum BracketError {
    /// A problem/candidate pair is malformed or exceeds the verifier envelope.
    InvalidInput {
        /// `primal` or `dual`.
        factor: &'static str,
        /// Stable refusal reason.
        reason: &'static str,
    },
    /// The independently rerun verifier panicked; no authority escaped.
    VerifierPanicked {
        /// `primal` or `dual`.
        factor: &'static str,
    },
    /// The verifier did not return a complete finite equilibrated certificate.
    VerifierRefused {
        /// `primal` or `dual`.
        factor: &'static str,
    },
    /// The outward-rounded product is not a finite usable QoI bound.
    ProductOverflow,
}

impl core::fmt::Display for BracketError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::InvalidInput { factor, reason } => {
                write!(f, "{factor} bracket input refused: {reason}")
            }
            Self::VerifierPanicked { factor } => {
                write!(f, "{factor} equilibrated verifier panicked")
            }
            Self::VerifierRefused { factor } => {
                write!(
                    f,
                    "{factor} equilibrated verifier produced no finite certificate"
                )
            }
            Self::ProductOverflow => f.write_str("Cauchy-Schwarz bracket product is not finite"),
        }
    }
}

impl std::error::Error for BracketError {}

fn verify_factor(
    factor: &'static str,
    problem: &MmsProblem,
    candidate: &[f64],
) -> Result<fs_verify::estimator::VerifierReport, BracketError> {
    if !(2..=MAX_BRACKET_MESH_NODES).contains(&problem.mesh().len()) {
        return Err(BracketError::InvalidInput {
            factor,
            reason: "mesh node count is outside 2..=MAX_BRACKET_MESH_NODES",
        });
    }
    if candidate.len() != problem.mesh().len() {
        return Err(BracketError::InvalidInput {
            factor,
            reason: "candidate length differs from mesh length",
        });
    }
    if !candidate.iter().all(|value| value.is_finite()) {
        return Err(BracketError::InvalidInput {
            factor,
            reason: "candidate contains a non-finite value",
        });
    }
    if !problem.mesh().iter().all(|value| value.is_finite())
        || !problem.mesh().windows(2).all(|pair| pair[0] < pair[1])
    {
        return Err(BracketError::InvalidInput {
            factor,
            reason: "mesh coordinates must be finite and strictly increasing",
        });
    }

    let report = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        verify(problem, candidate, f64::MAX)
    }))
    .map_err(|_| BracketError::VerifierPanicked { factor })?;
    let color_matches = matches!(
        &report.color,
        Some(Color::Verified { lo, hi })
            if lo.to_bits() == 0.0_f64.to_bits()
                && hi.to_bits() == report.bound.hi.to_bits()
    );
    if report.family != EstimatorFamily::EquilibratedFlux.id()
        || report.tolerance.to_bits() != f64::MAX.to_bits()
        || !report.accept
        || !color_matches
        || report.bound.lo.is_nan()
        || !report.bound.hi.is_finite()
        || report.bound.lo < 0.0
        || report.bound.hi < report.bound.lo
    {
        return Err(BracketError::VerifierRefused { factor });
    }
    Ok(report)
}

impl Bracket {
    /// The Cauchy–Schwarz bracket from two equilibrated energy-norm
    /// enclosures: `|a(e_u, e_z)| ≤ ‖e_u‖_E · ‖e_z‖_E`, outward-rounded.
    /// Both reports are independently recomputed here; public report fields are
    /// never accepted as authority.
    ///
    /// # Errors
    /// [`BracketError`] when either input is malformed, verification panics or
    /// refuses, or the product overflows.
    pub fn cauchy_schwarz(
        primal_problem: &MmsProblem,
        primal_candidate: &[f64],
        dual_problem: &MmsProblem,
        dual_candidate: &[f64],
    ) -> Result<Bracket, BracketError> {
        let primal = verify_factor("primal", primal_problem, primal_candidate)?;
        let dual = verify_factor("dual", dual_problem, dual_candidate)?;
        let raw_bound = primal.bound.hi * dual.bound.hi;
        let bound = up(raw_bound);
        if !bound.is_finite() || bound < 0.0 {
            return Err(BracketError::ProductOverflow);
        }
        Ok(Bracket {
            bound,
            source: format!(
                "cauchy-schwarz(equilibrated primal {:.3e} flux {:016x} x equilibrated dual {:.3e} flux {:016x})",
                primal.bound.hi, primal.flux_hash, dual.bound.hi, dual.flux_hash
            ),
        })
    }

    /// Outward-rounded energy-error product. This is diagnostic until a typed
    /// QoI-dual relation is verified.
    #[must_use]
    pub fn bound(&self) -> f64 {
        self.bound
    }

    /// Deterministic verifier audit label.
    #[must_use]
    pub fn source(&self) -> &str {
        &self.source
    }
}

/// The colored accept outcome.
#[derive(Debug, Clone, PartialEq)]
pub struct AcceptOutcome {
    /// Was the query discharged?
    pub accepted: bool,
    /// The color the answer carries.
    pub color: Color,
    /// True when malformed public inputs prevented an accept/reject decision.
    /// False includes both a valid accept and a valid over-tolerance rejection.
    pub refused: bool,
    /// The audit trail.
    pub audit: String,
}

fn malformed_refusal(estimator: &str, audit: String) -> AcceptOutcome {
    AcceptOutcome {
        accepted: false,
        color: Color::Estimated {
            estimator: estimator.to_string(),
            dispersion: f64::INFINITY,
        },
        refused: true,
        audit,
    }
}

/// The accept test. Color logic (mechanical, auditable):
/// - no acceptance path → rejected (estimated color on the estimate);
/// - DWR-only accept (`|η| ≤ tol`, no valid bracket) → ESTIMATED;
/// - a sealed energy-product bracket is retained in the audit, but cannot
///   promote or veto because the v0 query does not bind its dual relation.
#[must_use]
pub fn accept(query: &DwrQuery, dwr_abs: f64, bracket: Option<&Bracket>) -> AcceptOutcome {
    if !query.tolerance.is_finite() || query.tolerance <= 0.0 {
        return malformed_refusal(
            "dwr-invalid-tolerance",
            format!(
                "REFUSED: tolerance must be finite and positive, got {:.3e}",
                query.tolerance
            ),
        );
    }
    if let Some(b) = bracket
        && (!b.bound.is_finite() || b.bound < 0.0)
    {
        return malformed_refusal(
            "dwr-invalid-guaranteed-bracket",
            format!(
                "REFUSED: guaranteed bracket from {} has invalid bound {:.3e}",
                b.source, b.bound
            ),
        );
    }
    let dwr_is_usable = dwr_abs.is_finite() && dwr_abs >= 0.0;
    if !dwr_is_usable {
        return malformed_refusal(
            "dwr-invalid-estimate",
            format!(
                "REFUSED: DWR absolute estimate must be finite and non-negative, got {dwr_abs:.3e}"
            ),
        );
    }
    if dwr_abs <= query.tolerance {
        let bracket_note = bracket.map_or_else(
            || "no energy-product diagnostic".to_string(),
            |value| {
                format!(
                    "energy-product diagnostic {:.3e} from {}; QoI-dual relation unverified",
                    value.bound, value.source
                )
            },
        );
        return AcceptOutcome {
            accepted: true,
            color: Color::Estimated {
                estimator: if bracket.is_some() {
                    "dwr-with-unbound-energy-diagnostic".to_string()
                } else {
                    "dwr-unbracketed".to_string()
                },
                dispersion: dwr_abs,
            },
            refused: false,
            audit: format!(
                "estimated-only accept: dwr {:.3e} <= tol {:.3e}; {bracket_note}",
                dwr_abs, query.tolerance,
            ),
        };
    }
    let bracket_note = bracket.map_or_else(
        || "no energy-product diagnostic".to_string(),
        |value| {
            format!(
                "energy-product diagnostic {:.3e} from {}; QoI-dual relation unverified",
                value.bound, value.source
            )
        },
    );
    AcceptOutcome {
        accepted: false,
        color: Color::Estimated {
            estimator: "dwr-rejected".to_string(),
            dispersion: dwr_abs,
        },
        refused: false,
        audit: format!(
            "rejected: dwr {:.3e} > tol {:.3e}; {bracket_note}",
            dwr_abs, query.tolerance,
        ),
    }
}

/// The 1-D reference DWR estimator for integral QoIs
/// `J(u) = ∫_{w_lo}^{w_hi} u dx` over an fs-verify problem: the dual
/// `−z″ = 1_{[w_lo, w_hi]}` solves by P1 FEM on the ONCE-REFINED mesh
/// (the enriched dual), and the estimate is the dual-weighted residual
/// `η = r(z_f − I_h z_f)` with per-COARSE-element indicators.
#[derive(Debug, Clone)]
pub struct DwrOutput {
    /// `J(u_h)`.
    pub j_primal: f64,
    /// The signed estimate `η ≈ J(u) − J(u_h)`.
    pub eta: f64,
    /// Per-coarse-element |indicator| (refinement guidance).
    pub indicators: Vec<f64>,
}

/// Why the public DWR execution path refused an input or derived state.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum DwrError {
    /// Coarse mesh node count is outside the bounded execution envelope.
    MeshNodeCount {
        /// Supplied node count.
        count: usize,
        /// Required minimum.
        minimum: usize,
        /// Admitted maximum.
        maximum: usize,
    },
    /// Candidate node count exceeds the bounded execution envelope.
    CandidateNodeCount {
        /// Supplied value count.
        count: usize,
        /// Admitted maximum.
        maximum: usize,
    },
    /// Candidate and mesh shapes differ.
    CandidateLengthMismatch {
        /// Mesh node count.
        mesh_nodes: usize,
        /// Candidate value count.
        candidate_values: usize,
    },
    /// Manufactured-solution polynomial size is outside the bounded envelope.
    PolynomialCoefficientCount {
        /// Supplied coefficient count.
        count: usize,
        /// Required minimum.
        minimum: usize,
        /// Admitted maximum.
        maximum: usize,
    },
    /// A manufactured-solution coefficient is NaN or infinite.
    NonFinitePolynomialCoefficient {
        /// Coefficient index.
        index: usize,
    },
    /// A candidate nodal value is NaN or infinite.
    NonFiniteCandidate {
        /// Candidate index.
        index: usize,
    },
    /// A mesh coordinate is NaN or infinite.
    NonFiniteMeshNode {
        /// Mesh-node index.
        index: usize,
    },
    /// A mesh cell is not strictly increasing.
    NonIncreasingMeshCell {
        /// Left node / cell index.
        cell: usize,
    },
    /// Subtracting a cell's finite endpoints produced a non-finite width.
    NonFiniteCellWidth {
        /// Left node / cell index.
        cell: usize,
    },
    /// The once-refined midpoint is not strictly inside its coarse cell.
    NonInteriorMidpoint {
        /// Coarse cell index.
        cell: usize,
    },
    /// A coarse or refined cell has a non-finite reciprocal width.
    NonFiniteReciprocal {
        /// Coarse cell index.
        cell: usize,
        /// `None` for the coarse cell, `Some(0|1)` for a refined half.
        refined_half: Option<u8>,
    },
    /// The QoI integration window is non-finite or inverted.
    InvalidQoiWindow {
        /// Stable refusal reason.
        reason: &'static str,
    },
    /// Computing `2 * nodes - 1` overflowed or exceeded the refined cap.
    RefinedMeshSizeOverflow {
        /// Supplied coarse node count.
        mesh_nodes: usize,
    },
    /// The mesh/polynomial cross-product exceeds the bounded execution budget.
    WorkBudgetExceeded {
        /// Supplied coarse node count.
        mesh_nodes: usize,
        /// Supplied manufactured-solution coefficient count.
        polynomial_coefficients: usize,
        /// Estimated work, or `None` when the estimate itself overflowed.
        estimated_work: Option<usize>,
        /// Admitted maximum.
        maximum: usize,
    },
    /// Tridiagonal storage lengths are inconsistent.
    LinearSystemShape,
    /// Assembly or elimination produced a non-finite linear-system value.
    NonFiniteLinearSystem {
        /// Stable component/stage name.
        component: &'static str,
        /// Row index.
        index: usize,
    },
    /// Thomas elimination encountered a zero or non-finite pivot.
    InvalidLinearPivot {
        /// Pivot row.
        row: usize,
    },
    /// A derived quadrature, slope, residual, or output value is non-finite.
    NonFiniteDerived {
        /// Stable derived quantity name.
        quantity: &'static str,
        /// Cell/coefficient index when applicable.
        index: Option<usize>,
    },
}

impl core::fmt::Display for DwrError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::MeshNodeCount {
                count,
                minimum,
                maximum,
            } => write!(
                f,
                "DWR mesh has {count} nodes; expected {minimum}..={maximum}"
            ),
            Self::CandidateNodeCount { count, maximum } => {
                write!(f, "DWR candidate has {count} values; maximum is {maximum}")
            }
            Self::CandidateLengthMismatch {
                mesh_nodes,
                candidate_values,
            } => write!(
                f,
                "DWR candidate length {candidate_values} differs from mesh length {mesh_nodes}"
            ),
            Self::PolynomialCoefficientCount {
                count,
                minimum,
                maximum,
            } => write!(
                f,
                "DWR polynomial has {count} coefficients; expected {minimum}..={maximum}"
            ),
            Self::NonFinitePolynomialCoefficient { index } => {
                write!(f, "DWR polynomial coefficient {index} is non-finite")
            }
            Self::NonFiniteCandidate { index } => {
                write!(f, "DWR candidate value {index} is non-finite")
            }
            Self::NonFiniteMeshNode { index } => {
                write!(f, "DWR mesh node {index} is non-finite")
            }
            Self::NonIncreasingMeshCell { cell } => {
                write!(f, "DWR mesh cell {cell} is not strictly increasing")
            }
            Self::NonFiniteCellWidth { cell } => {
                write!(f, "DWR mesh cell {cell} has a non-finite width")
            }
            Self::NonInteriorMidpoint { cell } => write!(
                f,
                "DWR mesh cell {cell} has no representable strictly interior midpoint"
            ),
            Self::NonFiniteReciprocal { cell, refined_half } => match refined_half {
                Some(half) => write!(
                    f,
                    "DWR refined half {half} of coarse cell {cell} has a non-finite reciprocal width"
                ),
                None => write!(
                    f,
                    "DWR coarse cell {cell} has a non-finite reciprocal width"
                ),
            },
            Self::InvalidQoiWindow { reason } => {
                write!(f, "DWR QoI window refused: {reason}")
            }
            Self::RefinedMeshSizeOverflow { mesh_nodes } => write!(
                f,
                "DWR refined mesh size overflowed its bound for {mesh_nodes} coarse nodes"
            ),
            Self::WorkBudgetExceeded {
                mesh_nodes,
                polynomial_coefficients,
                estimated_work,
                maximum,
            } => match estimated_work {
                Some(work) => write!(
                    f,
                    "DWR work estimate {work} for {mesh_nodes} mesh nodes x {polynomial_coefficients} coefficients exceeds {maximum}"
                ),
                None => write!(
                    f,
                    "DWR work estimate overflowed for {mesh_nodes} mesh nodes x {polynomial_coefficients} coefficients (maximum {maximum})"
                ),
            },
            Self::LinearSystemShape => {
                f.write_str("DWR tridiagonal linear-system shapes are inconsistent")
            }
            Self::NonFiniteLinearSystem { component, index } => write!(
                f,
                "DWR linear-system {component} is non-finite at row {index}"
            ),
            Self::InvalidLinearPivot { row } => {
                write!(f, "DWR Thomas pivot {row} is zero or non-finite")
            }
            Self::NonFiniteDerived { quantity, index } => match index {
                Some(index) => write!(f, "DWR derived {quantity} is non-finite at index {index}"),
                None => write!(f, "DWR derived {quantity} is non-finite"),
            },
        }
    }
}

impl std::error::Error for DwrError {}

fn refined_node_count(mesh_nodes: usize) -> Result<usize, DwrError> {
    mesh_nodes
        .checked_mul(2)
        .and_then(|count| count.checked_sub(1))
        .filter(|&count| count <= MAX_DWR_REFINED_NODES)
        .ok_or(DwrError::RefinedMeshSizeOverflow { mesh_nodes })
}

fn dwr_work_units(mesh_nodes: usize, polynomial_coefficients: usize) -> Option<usize> {
    // Per coarse cell: 5 primal-QoI points, 10 dual-load points, and 10
    // residual points whose forcing evaluation is a Horner walk over the
    // polynomial. This deliberately overestimates the twice-derived forcing.
    let cells = mesh_nodes.checked_sub(1)?;
    let per_cell = polynomial_coefficients.checked_mul(10)?.checked_add(15)?;
    cells.checked_mul(per_cell)
}

fn validate_dwr_inputs(
    problem: &MmsProblem,
    candidate: &[f64],
    w_lo: f64,
    w_hi: f64,
) -> Result<usize, DwrError> {
    let mesh = problem.mesh();
    let coefficients = problem.exact_solution().coefficients();
    if !(2..=MAX_DWR_MESH_NODES).contains(&mesh.len()) {
        return Err(DwrError::MeshNodeCount {
            count: mesh.len(),
            minimum: 2,
            maximum: MAX_DWR_MESH_NODES,
        });
    }
    if candidate.len() > MAX_DWR_MESH_NODES {
        return Err(DwrError::CandidateNodeCount {
            count: candidate.len(),
            maximum: MAX_DWR_MESH_NODES,
        });
    }
    if candidate.len() != mesh.len() {
        return Err(DwrError::CandidateLengthMismatch {
            mesh_nodes: mesh.len(),
            candidate_values: candidate.len(),
        });
    }
    if !(1..=MAX_DWR_POLY_COEFFICIENTS).contains(&coefficients.len()) {
        return Err(DwrError::PolynomialCoefficientCount {
            count: coefficients.len(),
            minimum: 1,
            maximum: MAX_DWR_POLY_COEFFICIENTS,
        });
    }
    let estimated_work = dwr_work_units(mesh.len(), coefficients.len());
    if estimated_work.is_none_or(|work| work > MAX_DWR_WORK_UNITS) {
        return Err(DwrError::WorkBudgetExceeded {
            mesh_nodes: mesh.len(),
            polynomial_coefficients: coefficients.len(),
            estimated_work,
            maximum: MAX_DWR_WORK_UNITS,
        });
    }
    if !w_lo.is_finite() || !w_hi.is_finite() {
        return Err(DwrError::InvalidQoiWindow {
            reason: "endpoints must be finite",
        });
    }
    if w_lo >= w_hi {
        return Err(DwrError::InvalidQoiWindow {
            reason: "lower endpoint must be strictly below upper endpoint",
        });
    }
    if let Some(index) = coefficients.iter().position(|value| !value.is_finite()) {
        return Err(DwrError::NonFinitePolynomialCoefficient { index });
    }
    if let Some(index) = candidate.iter().position(|value| !value.is_finite()) {
        return Err(DwrError::NonFiniteCandidate { index });
    }
    if let Some(index) = mesh.iter().position(|value| !value.is_finite()) {
        return Err(DwrError::NonFiniteMeshNode { index });
    }
    for (cell, nodes) in mesh.windows(2).enumerate() {
        let (x0, x1) = (nodes[0], nodes[1]);
        if x0 >= x1 {
            return Err(DwrError::NonIncreasingMeshCell { cell });
        }
        let width = x1 - x0;
        if !width.is_finite() {
            return Err(DwrError::NonFiniteCellWidth { cell });
        }
        if !(1.0 / width).is_finite() {
            return Err(DwrError::NonFiniteReciprocal {
                cell,
                refined_half: None,
            });
        }
        let midpoint = f64::midpoint(x0, x1);
        if !(x0 < midpoint && midpoint < x1) {
            return Err(DwrError::NonInteriorMidpoint { cell });
        }
        for (half, half_width) in [(0_u8, midpoint - x0), (1_u8, x1 - midpoint)] {
            if !half_width.is_finite() || !(1.0 / half_width).is_finite() {
                return Err(DwrError::NonFiniteReciprocal {
                    cell,
                    refined_half: Some(half),
                });
            }
        }
    }
    refined_node_count(mesh.len())
}

fn non_finite_derived(quantity: &'static str, index: Option<usize>) -> DwrError {
    DwrError::NonFiniteDerived { quantity, index }
}

fn thomas_solve(sub: &[f64], diag: &[f64], sup: &[f64], rhs: &mut [f64]) -> Result<(), DwrError> {
    let n = rhs.len();
    if sub.len() != n || diag.len() != n || sup.len() != n {
        return Err(DwrError::LinearSystemShape);
    }
    if n == 0 {
        return Ok(());
    }
    for (component, values) in [
        ("subdiagonal", sub),
        ("diagonal", diag),
        ("superdiagonal", sup),
    ] {
        if let Some(index) = values.iter().position(|value| !value.is_finite()) {
            return Err(DwrError::NonFiniteLinearSystem { component, index });
        }
    }
    if let Some(index) = rhs.iter().position(|value| !value.is_finite()) {
        return Err(DwrError::NonFiniteLinearSystem {
            component: "right-hand side",
            index,
        });
    }
    let mut c = vec![0.0f64; n];
    let mut d = diag[0];
    if !d.is_finite() || d == 0.0 {
        return Err(DwrError::InvalidLinearPivot { row: 0 });
    }
    if n > 1 {
        c[0] = sup[0] / d;
        if !c[0].is_finite() {
            return Err(DwrError::NonFiniteLinearSystem {
                component: "forward coefficient",
                index: 0,
            });
        }
    }
    rhs[0] /= d;
    if !rhs[0].is_finite() {
        return Err(DwrError::NonFiniteLinearSystem {
            component: "forward right-hand side",
            index: 0,
        });
    }
    for i in 1..n {
        d = diag[i] - sub[i] * c[i - 1];
        if !d.is_finite() || d == 0.0 {
            return Err(DwrError::InvalidLinearPivot { row: i });
        }
        if i < n - 1 {
            c[i] = sup[i] / d;
            if !c[i].is_finite() {
                return Err(DwrError::NonFiniteLinearSystem {
                    component: "forward coefficient",
                    index: i,
                });
            }
        }
        rhs[i] = (rhs[i] - sub[i] * rhs[i - 1]) / d;
        if !rhs[i].is_finite() {
            return Err(DwrError::NonFiniteLinearSystem {
                component: "forward right-hand side",
                index: i,
            });
        }
    }
    for i in (0..n - 1).rev() {
        rhs[i] -= c[i] * rhs[i + 1];
        if !rhs[i].is_finite() {
            return Err(DwrError::NonFiniteLinearSystem {
                component: "back-substitution",
                index: i,
            });
        }
    }
    Ok(())
}

/// P1 FEM solve of `−z″ = w` (zero Dirichlet BC) on `mesh`, with
/// `w = 1` on `[w_lo, w_hi]` — deterministic Thomas solve.
fn dual_solve(mesh: &[f64], w_lo: f64, w_hi: f64) -> Result<Vec<f64>, DwrError> {
    let n = mesh.len();
    let free = n.saturating_sub(2);
    if free == 0 {
        return Ok(vec![0.0; n]);
    }
    let mut sub = vec![0.0f64; free];
    let mut diag = vec![0.0f64; free];
    let mut sup = vec![0.0f64; free];
    let mut rhs = vec![0.0f64; free];
    for e in 0..n - 1 {
        let h = mesh[e + 1] - mesh[e];
        let k = 1.0 / h;
        if !h.is_finite() || h <= 0.0 || !k.is_finite() {
            return Err(DwrError::NonFiniteReciprocal {
                cell: e / 2,
                refined_half: Some((e % 2) as u8),
            });
        }
        for (a, b, v) in [(e, e, k), (e + 1, e + 1, k), (e, e + 1, -k), (e + 1, e, -k)] {
            if a >= 1 && a <= free && b >= 1 && b <= free {
                let (i, j) = (a - 1, b - 1);
                if i == j {
                    diag[i] += v;
                } else if j == i + 1 {
                    sup[i] += v;
                } else {
                    sub[i] += v;
                }
            }
        }
        // Load: ∫ w φ_a over the element (Gauss).
        for (gx, gw) in gauss5(mesh[e], mesh[e + 1]) {
            if !gx.is_finite() || !gw.is_finite() {
                return Err(non_finite_derived("dual quadrature", Some(e)));
            }
            let w = f64::from(u8::from(gx >= w_lo && gx <= w_hi));
            let xi = (gx - mesh[e]) / h;
            if !xi.is_finite() {
                return Err(non_finite_derived("dual reference coordinate", Some(e)));
            }
            for (node, shape) in [(e, 1.0 - xi), (e + 1, xi)] {
                if node >= 1 && node <= free {
                    let contribution = gw * w * shape;
                    let updated = rhs[node - 1] + contribution;
                    if !contribution.is_finite() || !updated.is_finite() {
                        return Err(DwrError::NonFiniteLinearSystem {
                            component: "assembled right-hand side",
                            index: node - 1,
                        });
                    }
                    rhs[node - 1] = updated;
                }
            }
        }
    }
    thomas_solve(&sub, &diag, &sup, &mut rhs)?;
    let mut z = vec![0.0f64; n];
    z[1..=free].copy_from_slice(&rhs);
    Ok(z)
}

fn refine(mesh: &[f64], refined_nodes: usize) -> Result<Vec<f64>, DwrError> {
    if refined_node_count(mesh.len())? != refined_nodes {
        return Err(DwrError::RefinedMeshSizeOverflow {
            mesh_nodes: mesh.len(),
        });
    }
    let mut out = Vec::with_capacity(refined_nodes);
    for e in 0..mesh.len() - 1 {
        let midpoint = f64::midpoint(mesh[e], mesh[e + 1]);
        if !(mesh[e] < midpoint && midpoint < mesh[e + 1]) {
            return Err(DwrError::NonInteriorMidpoint { cell: e });
        }
        out.push(mesh[e]);
        out.push(midpoint);
    }
    let Some(&last) = mesh.last() else {
        return Err(DwrError::MeshNodeCount {
            count: 0,
            minimum: 2,
            maximum: MAX_DWR_MESH_NODES,
        });
    };
    out.push(last);
    Ok(out)
}

/// Run the 1-D goal-oriented estimate (see [`DwrOutput`]).
///
/// # Errors
/// [`DwrError`] when public inputs exceed the bounded execution envelope or any
/// input/derived numerical value cannot support a finite Estimated diagnostic.
pub fn dwr_integral_qoi(
    problem: &MmsProblem,
    candidate: &[f64],
    w_lo: f64,
    w_hi: f64,
) -> Result<DwrOutput, DwrError> {
    let refined_nodes = validate_dwr_inputs(problem, candidate, w_lo, w_hi)?;
    let mesh = problem.mesh();
    let f = problem.forcing();
    if let Some(index) = f.coefficients().iter().position(|value| !value.is_finite()) {
        return Err(non_finite_derived("forcing coefficient", Some(index)));
    }
    // J(u_h): the P1 interpolant integrated over the window.
    let mut j_primal = 0.0f64;
    for e in 0..mesh.len() - 1 {
        let h = mesh[e + 1] - mesh[e];
        for (gx, gw) in gauss5(mesh[e], mesh[e + 1]) {
            if !gx.is_finite() || !gw.is_finite() {
                return Err(non_finite_derived("primal quadrature", Some(e)));
            }
            if gx >= w_lo && gx <= w_hi {
                let xi = (gx - mesh[e]) / h;
                let interpolated = (1.0 - xi) * candidate[e] + xi * candidate[e + 1];
                let contribution = gw * interpolated;
                let updated = j_primal + contribution;
                if !xi.is_finite()
                    || !interpolated.is_finite()
                    || !contribution.is_finite()
                    || !updated.is_finite()
                {
                    return Err(non_finite_derived("primal QoI", Some(e)));
                }
                j_primal = updated;
            }
        }
    }
    // Enriched dual on the refined mesh.
    let fine = refine(mesh, refined_nodes)?;
    let z = dual_solve(&fine, w_lo, w_hi)?;
    // Coarse-node interpolant of z, subtracted (Galerkin orthogonality
    // makes the coarse part vanish; the fine remainder drives η).
    let mut eta = 0.0f64;
    let mut indicators = vec![0.0f64; mesh.len() - 1];
    for e in 0..mesh.len() - 1 {
        let (x0, x1) = (mesh[e], mesh[e + 1]);
        let slope = (candidate[e + 1] - candidate[e]) / (x1 - x0);
        if !slope.is_finite() {
            return Err(non_finite_derived("primal slope", Some(e)));
        }
        let (z0, z1) = (z[2 * e], z[2 * e + 2]);
        let mut local = 0.0f64;
        // Two fine halves of the coarse element.
        for half in 0..2usize {
            let (fa, fb) = (fine[2 * e + half], fine[2 * e + half + 1]);
            let (za, zb) = (z[2 * e + half], z[2 * e + half + 1]);
            let zslope = (zb - za) / (fb - fa);
            // Coarse interpolant of z on this fine piece.
            let islope = (z1 - z0) / (x1 - x0);
            if !zslope.is_finite() || !islope.is_finite() {
                return Err(non_finite_derived("dual slope", Some(e)));
            }
            for (gx, gw) in gauss5(fa, fb) {
                let xi_f = (gx - fa) / (fb - fa);
                let zf = (1.0 - xi_f) * za + xi_f * zb;
                let zi = z0 + (gx - x0) * islope;
                // r(v) = ∫ f v − ∫ u_h′ v′ with v = z_f − I_h z_f.
                let forcing = f.eval(gx);
                let contribution = gw * (forcing * (zf - zi) - slope * (zslope - islope));
                let updated = local + contribution;
                if !gx.is_finite()
                    || !gw.is_finite()
                    || !xi_f.is_finite()
                    || !zf.is_finite()
                    || !zi.is_finite()
                    || !forcing.is_finite()
                    || !contribution.is_finite()
                    || !updated.is_finite()
                {
                    return Err(non_finite_derived("dual-weighted residual", Some(e)));
                }
                local = updated;
            }
        }
        let updated_eta = eta + local;
        if !updated_eta.is_finite() {
            return Err(non_finite_derived("global DWR estimate", Some(e)));
        }
        eta = updated_eta;
        indicators[e] = local.abs();
    }
    if !j_primal.is_finite()
        || !eta.is_finite()
        || indicators.iter().any(|indicator| !indicator.is_finite())
    {
        return Err(non_finite_derived("DWR output", None));
    }
    Ok(DwrOutput {
        j_primal,
        eta,
        indicators,
    })
}

#[cfg(test)]
mod execution_tests {
    use super::*;

    #[test]
    fn zero_interior_linear_and_dual_systems_are_total() {
        let mut rhs = Vec::new();
        thomas_solve(&[], &[], &[], &mut rhs).expect("empty system is solved");
        let dual = dual_solve(&[0.0, 1.0], 0.0, 1.0).expect("boundary-only dual");
        assert_eq!(dual, vec![0.0, 0.0]);
    }

    #[test]
    fn refined_size_arithmetic_is_checked() {
        assert!(matches!(
            refined_node_count(usize::MAX),
            Err(DwrError::RefinedMeshSizeOverflow {
                mesh_nodes: usize::MAX
            })
        ));
        assert!(dwr_work_units(MAX_DWR_MESH_NODES, 1).expect("bounded work") <= MAX_DWR_WORK_UNITS);
        assert!(
            dwr_work_units(2, MAX_DWR_POLY_COEFFICIENTS).expect("bounded work")
                <= MAX_DWR_WORK_UNITS
        );
    }
}
