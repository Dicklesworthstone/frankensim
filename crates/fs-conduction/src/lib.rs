//! fs-conduction — STEADY heat conduction on tetrahedral complexes
//! (bead `frankensim-extreal-program-f85xj.5.1`). Layer: L3 FLUX.
//!
//! The strong form solved here is
//!
//! ```text
//!   −∇·(k(T, x) ∇T) = f            in Ω
//!                 T = T_D          on Γ_D
//!         (−k∇T)·n  = q_n          on Γ_N
//!         (−k∇T)·n  = h (T − T_ref) on Γ_R
//! ```
//!
//! with `k` a symmetric positive-definite conductivity TENSOR that may
//! depend on temperature, `f` a volumetric source in W/m³, `q_n` the
//! prescribed OUTWARD heat flux in W/m² (positive = leaving Ω), and
//! `h`/`T_ref` the convective transfer coefficient (W/(m²·K)) and its
//! declared reference temperature (K).
//!
//! Discretization: continuous P₁ Lagrange (the FEEC 0-form / Whitney
//! vertex-hat space on `fs-feec`'s tet complexes). Element geometry —
//! signed volumes and the constant barycentric gradients ∇λ_a — comes
//! from [`fs_feec::element_geometry`]; nothing here re-derives it.
//! Assembly stages triplets into `fs-sparse`'s [`fs_sparse::Coo`], whose
//! canonicalization makes the matrix a pure function of the staged
//! multiset, so tile interleaving cannot move a bit.
//!
//! What travels WITH a solve: the `fs-matdb` property-usage receipts
//! that produced every conductivity number (see [`material`]), the typed
//! residual claim of every linear solve RE-MEASURED as a Euclidean
//! residual by this crate (see [`solve::LinearSolveEvidence`]), and an
//! energy-balance closure over the assembled operator.
//!
//! # No-claim boundaries (short form; `CONTRACT.md` is authoritative)
//!
//! - STEADY ONLY. There is no time derivative, no heat capacity, and no
//!   transient staging in this crate.
//! - NO radiation and NO convection PHYSICS. The Robin row is a
//!   convective boundary *coupling*: `h` is an input, never a computed
//!   correlation. Surface-to-ambient radiation and enclosure exchange
//!   are separate beads that couple in through the same row.
//! - Observed convergence orders are OBSERVED. The G1 battery fits and
//!   gates slopes on the fixture ladders it runs; it is not a proof of
//!   the order for arbitrary meshes, coefficients, or data.
//! - A converged residual is not an error bound. This crate reports
//!   residuals and an energy-balance closure; it does not ship a
//!   certified bound on `‖T − T_h‖`.

pub mod adjoint;
pub mod assemble;
pub mod bc;
pub mod field;
pub mod fixtures;
pub mod material;
pub mod mesh;
pub mod solve;

use core::fmt;

pub use assemble::{ASSEMBLY_TILE, AssembledSystem, DofMap, assemble_jacobian, assemble_operator};
pub use bc::{ThermalBc, ThermalBoundary, ThermalBoundaryBuilder};
pub use field::ScalarField;
pub use material::{
    CONDUCTIVITY_DIMS, ConductivityModel, ConductivityTable, ProvenanceClass, TemperatureSpan,
};
pub use mesh::{BoundaryFace, ConductionMesh};
pub use solve::{
    ConductionProblem, ConductionReport, ConductionSolution, ConductionSolver, ConductionState,
    EnergyBalance, InitialGuess, LineSearch, LinearConfig, LinearSolveEvidence, Nonlinearity,
    SolveConfig, StopReason, StopRule, solve,
};

/// Crate version, re-exported for provenance stamping.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

/// SI exponents (m, kg, s, K, A, mol) of a temperature.
pub const TEMPERATURE_DIMS: fs_qty::Dims = fs_qty::Dims([0, 0, 0, 1, 0, 0]);
/// SI exponents of a heat flux density, W/m².
pub const HEAT_FLUX_DIMS: fs_qty::Dims = fs_qty::Dims([0, 1, -3, 0, 0, 0]);
/// SI exponents of a convective transfer coefficient, W/(m²·K).
pub const HTC_DIMS: fs_qty::Dims = fs_qty::Dims([0, 1, -3, -1, 0, 0]);
/// SI exponents of a volumetric heat source, W/m³.
pub const VOLUMETRIC_SOURCE_DIMS: fs_qty::Dims = fs_qty::Dims([-1, 1, -3, 0, 0, 0]);

/// Everything this crate refuses, as typed values. Refusals name the
/// rule, the offending quantity, and (where one exists) the fix — no
/// panics cross the crate boundary and no path silently repairs a
/// physically meaningless input.
#[derive(Debug, Clone, PartialEq)]
pub enum ConductionError {
    /// A cancellation request was observed at a tile boundary. The stage
    /// name and the tile index that was ABOUT to run are recorded, so a
    /// driver can report where the drain happened.
    Cancelled {
        /// Which stage drained (`"assemble-elements"`, `"newton-step"`, …).
        stage: &'static str,
        /// The next unprocessed tile / iteration index.
        at: usize,
    },
    /// Mesh input is structurally unusable.
    Mesh {
        /// What is wrong.
        what: String,
        /// How to fix it.
        fix: String,
    },
    /// A tetrahedron is degenerate (zero or near-zero signed volume).
    DegenerateElement {
        /// Element index.
        element: usize,
        /// Its signed volume.
        signed_volume: f64,
    },
    /// A supplied number is not finite where the model requires one.
    NonFinite {
        /// Which field.
        field: &'static str,
        /// Raw bits of the offending value.
        bits: u64,
    },
    /// A dimensioned input carries the wrong SI exponents.
    Dimensions {
        /// What was being read.
        context: String,
        /// SI exponents the contract demands.
        expected: [i8; 6],
        /// SI exponents supplied.
        found: [i8; 6],
    },
    /// A boundary region overlaps one already declared.
    OverlappingRegion {
        /// The region being declared.
        region: String,
        /// The region that already owns the face.
        owner: String,
        /// The contested boundary-face slot.
        face: usize,
    },
    /// A region name was declared twice.
    DuplicateRegion {
        /// The repeated name.
        region: String,
    },
    /// Boundary faces carry no declared condition and no explicit
    /// adiabatic remainder was requested.
    UntaggedBoundary {
        /// How many boundary faces are untagged.
        count: usize,
        /// The lowest untagged boundary-face slot.
        first: usize,
    },
    /// The problem has no free degree of freedom (every vertex is
    /// Dirichlet-pinned), so there is nothing to solve.
    NoFreeDofs,
    /// A pure-Neumann/Robin-free problem: with no Dirichlet row and no
    /// Robin row the steady operator is singular up to a constant.
    SingularPureNeumann,
    /// A conductivity model is not admissible.
    Conductivity {
        /// Diagnosis.
        what: String,
    },
    /// A temperature left the span the conductivity table was sampled
    /// over. Extrapolation is never implicit.
    OutsideTemperatureSpan {
        /// The temperature asked for.
        temperature: f64,
        /// Lowest sampled temperature.
        low: f64,
        /// Highest sampled temperature.
        high: f64,
    },
    /// An `fs-matdb` query refused (unknown property, out of validity,
    /// ambiguous selection, …). The upstream refusal is carried verbatim.
    MaterialQuery {
        /// The property name asked for.
        property: String,
        /// The temperature the query was made at.
        temperature: f64,
        /// The upstream refusal, rendered.
        upstream: String,
    },
    /// A `fs-scenario` boundary row cannot be lowered into a thermal
    /// condition this crate solves.
    ScenarioRow {
        /// The region the row is attached to.
        region: String,
        /// Why it cannot be lowered.
        what: String,
        /// What to supply instead.
        fix: String,
    },
    /// A field's nodal array has the wrong length.
    FieldLength {
        /// Which field.
        field: &'static str,
        /// Expected length.
        expected: usize,
        /// Supplied length.
        found: usize,
    },
    /// A configuration value is outside its admissible range.
    Config {
        /// The parameter name.
        parameter: &'static str,
        /// Why it is inadmissible.
        what: String,
    },
    /// The nonlinear iteration exhausted its budget without meeting the
    /// declared stop rule.
    NotConverged {
        /// Iterations performed.
        iterations: usize,
        /// Final residual 2-norm on free dofs.
        residual: f64,
        /// The stop rule's threshold at that scale.
        threshold: f64,
    },
    /// The globalization could not find an admissible step.
    LineSearchFailed {
        /// Newton iteration index.
        iteration: usize,
        /// Backtracks attempted.
        backtracks: usize,
        /// Smallest step length tried.
        smallest_step: f64,
    },
    /// A linear solve failed to reach its tolerance in the recomputed
    /// EUCLIDEAN residual (not the solver's recurrence estimate).
    LinearSolveFailed {
        /// Nonlinear iteration index.
        iteration: usize,
        /// Krylov iterations performed.
        krylov_iterations: usize,
        /// Recomputed `‖b − Ax‖₂/‖b‖₂`.
        true_relative_residual: f64,
        /// The tolerance it was asked for.
        tolerance: f64,
    },
    /// A resumable snapshot refused.
    Snapshot {
        /// The rendered upstream envelope/codec refusal.
        upstream: String,
    },
}

impl fmt::Display for ConductionError {
    #[allow(clippy::too_many_lines)]
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ConductionError::Cancelled { stage, at } => {
                write!(f, "cancelled during {stage} at tile/iteration {at}")
            }
            ConductionError::Mesh { what, fix } => write!(f, "mesh refusal: {what}; fix: {fix}"),
            ConductionError::DegenerateElement {
                element,
                signed_volume,
            } => write!(
                f,
                "tet {element} is degenerate (signed volume {signed_volume:e}); \
                 a conduction operator cannot be assembled on it"
            ),
            ConductionError::NonFinite { field, bits } => {
                write!(f, "{field} is not finite (bits 0x{bits:016x})")
            }
            ConductionError::Dimensions {
                context,
                expected,
                found,
            } => write!(
                f,
                "{context}: dimensions {found:?} do not match the required {expected:?} \
                 (m, kg, s, K, A, mol)"
            ),
            ConductionError::OverlappingRegion {
                region,
                owner,
                face,
            } => write!(
                f,
                "boundary region {region:?} claims face {face}, already owned by {owner:?}; \
                 regions must partition the boundary"
            ),
            ConductionError::DuplicateRegion { region } => {
                write!(f, "boundary region {region:?} declared twice")
            }
            ConductionError::UntaggedBoundary { count, first } => write!(
                f,
                "{count} boundary faces carry no condition (first is slot {first}); \
                 declare them or call adiabatic_remainder() to say so explicitly"
            ),
            ConductionError::NoFreeDofs => write!(
                f,
                "every vertex is Dirichlet-pinned: there is no unknown to solve for"
            ),
            ConductionError::SingularPureNeumann => write!(
                f,
                "no Dirichlet and no Robin row: the steady conduction operator is singular \
                 (T + const is also a solution); pin a temperature or add a convective row"
            ),
            ConductionError::Conductivity { what } => {
                write!(f, "inadmissible conductivity model: {what}")
            }
            ConductionError::OutsideTemperatureSpan {
                temperature,
                low,
                high,
            } => write!(
                f,
                "temperature {temperature} K is outside the sampled conductivity span \
                 [{low}, {high}] K; this crate never extrapolates material data"
            ),
            ConductionError::MaterialQuery {
                property,
                temperature,
                upstream,
            } => write!(
                f,
                "fs-matdb refused {property:?} at T = {temperature} K: {upstream}"
            ),
            ConductionError::ScenarioRow { region, what, fix } => write!(
                f,
                "scenario row on {region:?} cannot be lowered: {what}; fix: {fix}"
            ),
            ConductionError::FieldLength {
                field,
                expected,
                found,
            } => write!(f, "{field} has {found} nodal values, expected {expected}"),
            ConductionError::Config { parameter, what } => {
                write!(f, "configuration parameter {parameter}: {what}")
            }
            ConductionError::NotConverged {
                iterations,
                residual,
                threshold,
            } => write!(
                f,
                "nonlinear iteration did not converge in {iterations} iterations: \
                 residual {residual:e} > threshold {threshold:e}"
            ),
            ConductionError::LineSearchFailed {
                iteration,
                backtracks,
                smallest_step,
            } => write!(
                f,
                "Armijo backtracking failed at Newton iteration {iteration} after \
                 {backtracks} backtracks (smallest step {smallest_step:e})"
            ),
            ConductionError::LinearSolveFailed {
                iteration,
                krylov_iterations,
                true_relative_residual,
                tolerance,
            } => write!(
                f,
                "linear solve at nonlinear iteration {iteration} reached recomputed \
                 Euclidean relative residual {true_relative_residual:e} after \
                 {krylov_iterations} Krylov iterations, tolerance was {tolerance:e}"
            ),
            ConductionError::Snapshot { upstream } => {
                write!(f, "solver-state snapshot refused: {upstream}")
            }
        }
    }
}

impl std::error::Error for ConductionError {}

impl ConductionError {
    /// A stable machine-readable rule slug for structured logs.
    #[must_use]
    pub const fn rule(&self) -> &'static str {
        match self {
            ConductionError::Cancelled { .. } => "conduction-cancelled",
            ConductionError::Mesh { .. } => "conduction-mesh",
            ConductionError::DegenerateElement { .. } => "conduction-degenerate-element",
            ConductionError::NonFinite { .. } => "conduction-non-finite",
            ConductionError::Dimensions { .. } => "conduction-dimensions",
            ConductionError::OverlappingRegion { .. } => "conduction-overlapping-region",
            ConductionError::DuplicateRegion { .. } => "conduction-duplicate-region",
            ConductionError::UntaggedBoundary { .. } => "conduction-untagged-boundary",
            ConductionError::NoFreeDofs => "conduction-no-free-dofs",
            ConductionError::SingularPureNeumann => "conduction-singular-pure-neumann",
            ConductionError::Conductivity { .. } => "conduction-conductivity",
            ConductionError::OutsideTemperatureSpan { .. } => "conduction-outside-span",
            ConductionError::MaterialQuery { .. } => "conduction-material-query",
            ConductionError::ScenarioRow { .. } => "conduction-scenario-row",
            ConductionError::FieldLength { .. } => "conduction-field-length",
            ConductionError::Config { .. } => "conduction-config",
            ConductionError::NotConverged { .. } => "conduction-not-converged",
            ConductionError::LineSearchFailed { .. } => "conduction-line-search",
            ConductionError::LinearSolveFailed { .. } => "conduction-linear-solve",
            ConductionError::Snapshot { .. } => "conduction-snapshot",
        }
    }
}

pub(crate) fn require_finite(field: &'static str, value: f64) -> Result<f64, ConductionError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(ConductionError::NonFinite {
            field,
            bits: value.to_bits(),
        })
    }
}
