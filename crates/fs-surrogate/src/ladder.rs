//! THE ABSTRACTION LADDER (addendum Proposal A, bead knh1.4; [F] —
//! behind the `abstraction-ladder` feature): reduced-order VIEWS at
//! multiple abstraction levels whose fidelity to the level below is a
//! QUANTIFIED, leak-alarmed estimator — the operator reasons at the concept
//! level until the estimator says the abstraction is no longer faithful
//! HERE, and drills down only there.
//!
//! THE BEACHHEAD: REDUCED-BASIS for the affine-parametric elliptic
//! family `−(a(x;μ) u′)′ = f` with `a = 1 + μ·χ` — offline snapshots +
//! energy-orthonormal basis + online k×k Galerkin, and the TEXTBOOK
//! a-posteriori bound `‖u−u_rb‖_a ≤ ‖r‖_{V′}/√α_LB` with the residual
//! dual norm assembled offline and the coercivity floor exact for the
//! affine family. The compliance QoI estimator retains both the classic
//! squared energy estimator and the floating reduced solve's computable
//! Galerkin defect.
//!
//! EVIDENCE HONESTY (bead y6yv): every color this module emits is
//! ESTIMATED. The RB bound is real mathematics over the reals, but it
//! is COMPUTED here in round-to-nearest f64 with no outward rounding,
//! no Riesz-solve error control, and no linear-solve certificate — a
//! textbook bound evaluated in floating point is NOT an executable
//! enclosure, so it cannot mint `Color::Verified`. Likewise level 0 is
//! the DECLARED FE truth semantics, but its Thomas solve is plain
//! floating point: a point value under declared semantics, not a
//! verified interval. The certified-ladder DESTINATION (outward-rounded
//! residual/Riesz/solve certificates upgrading these colors to
//! Verified) is recorded in the CONTRACT as future scope; the type
//! system keeps today's estimates from impersonating it.

use core::{convert::Infallible, fmt, mem::size_of, ops::ControlFlow};
use fs_alloc::{LeaseReceipt, OperationMemoryLease};
use fs_evidence::{Color, NumericalCertificate, validate_color_payload};
use fs_exec::{
    Budget, CancelGate, Cancelled, RunError, RunId, RunReport, TileKernel, TilePlan, TilePool,
};
use std::sync::Mutex;

/// Largest full-order dimension a ladder will allocate (n vectors of n
/// f64 during training: 1<<20 keeps worst-case snapshots well bounded).
pub const MAX_TRUTH_NODES: usize = 1 << 20;
/// Largest RB basis size per rung.
pub const MAX_RB_DIM: usize = 512;
/// Largest concept grid / probe count.
pub const MAX_CONCEPT_POINTS: usize = 4096;
/// Largest number of RB rungs in one ladder.
pub const MAX_RB_LEVELS: usize = 16;
/// Maximum aggregate basis/scratch/dense storage estimated for one
/// training plan (512 MiB).
pub const MAX_TRAINING_BYTES: usize = 512 * 1024 * 1024;
/// Maximum scalar work estimate admitted for one training plan.
pub const MAX_TRAINING_WORK_UNITS: usize = 1_000_000_000;
/// Maximum points on either axis of an RB-coverage battery.
pub const MAX_COVERAGE_AXIS_POINTS: usize = 4096;
/// Maximum Cartesian-product size of an RB-coverage battery.
pub const MAX_COVERAGE_QUERIES: usize = 1_000_000;
/// Maximum conservative scalar work admitted for one synchronous coverage
/// battery, including every possible RB descent and one truth fallback per
/// parameter plus every tolerance comparison.
pub const MAX_COVERAGE_WORK_UNITS: usize = 1_000_000_000;
/// Maximum conservative live scratch admitted for one production coverage
/// run, excluding the immutable ladder and the caller-owned exact battery
/// plan (512 MiB).
pub const MAX_COVERAGE_MEMORY_BYTES: usize = 512 * 1024 * 1024;
/// Smaller work ceiling for the compatibility synchronous wrapper. Production
/// batteries use [`rb_coverage_scoped`] so they have explicit cancellation,
/// drain, progress, and memory-accounting semantics.
pub const MAX_SYNCHRONOUS_COVERAGE_WORK_UNITS: usize = 10_000_000;
/// Maximum logical scalar work between cooperative cancellation checkpoints
/// inside one coverage parameter tile.
pub const RB_COVERAGE_CHECKPOINT_WORK_UNITS: u64 = 256;

const DEFAULT_CONCEPT_GRID_POINTS: usize = 5;
const DEFAULT_CONCEPT_PROBES: usize = 7;

const CONCEPT_ESTIMATOR_ID: &str = "fs-surrogate.concept-cross-rung-v1";
const RB_ESTIMATOR_ID: &str = "fs-surrogate.rb-a-posteriori-f64-v1";
const TRUTH_ESTIMATOR_ID: &str = "fs-surrogate.fe-truth-f64-no-certificate-v1";
const COVERAGE_ESTIMATOR_ID: &str = "fs-surrogate.rb-coverage-f64-v2";
const COVERAGE_KERNEL_NAME: &str = "fs-surrogate/rb-coverage-v2";

/// A structured ladder refusal (Decalogue P10): every invalid input or
/// non-finite derived quantity is a named error, never a panic, NaN, or
/// silently wrong color.
#[derive(Debug, Clone, PartialEq)]
pub enum SurrogateError {
    /// Full-order dimension outside `1..=MAX_TRUTH_NODES`.
    InvalidDimension {
        /// The rejected n.
        n: usize,
    },
    /// A μ-range that is non-finite, empty, or reversed.
    InvalidRange {
        /// The rejected range.
        lo: f64,
        /// Upper end.
        hi: f64,
    },
    /// A parameter that breaks coercivity (`1 + μ ≤ 0`) or is
    /// non-finite — the bilinear form is no longer elliptic and every
    /// bound downstream is meaningless.
    InvalidCoercivity {
        /// The rejected μ.
        mu: f64,
    },
    /// A query parameter outside the rung's declared training range —
    /// estimator constants were calibrated on the range, so extrapolated
    /// ladder claims are refused.
    OutOfRange {
        /// The rejected μ.
        mu: f64,
        /// The declared range.
        range: (f64, f64),
    },
    /// A tolerance that is non-finite or non-positive.
    InvalidTolerance {
        /// The rejected tolerance.
        tol: f64,
    },
    /// Basis/grid/probe geometry outside its operation's allowed bounds.
    InvalidGeometry {
        /// Which knob.
        what: &'static str,
        /// The rejected value.
        got: usize,
    },
    /// A vector supplied to a public numerical operation has the wrong shape.
    InvalidVectorLength {
        /// Operation/operand identity.
        what: &'static str,
        /// Required length.
        expected: usize,
        /// Supplied length.
        got: usize,
    },
    /// An input vector contains a NaN or infinity.
    NonFiniteInput {
        /// Operation/operand identity.
        what: &'static str,
        /// First invalid element.
        index: usize,
    },
    /// A requested start rung does not exist.
    InvalidLevel {
        /// Requested rung.
        requested: usize,
        /// Highest declared rung.
        top: usize,
    },
    /// A rung is not bound to the ladder's truth-space/range identity.
    FamilyMismatch {
        /// Mismatched component.
        what: &'static str,
    },
    /// An RB ladder has too many rungs.
    TooManyRbLevels {
        /// Requested rung count.
        got: usize,
        /// Configured cap.
        max: usize,
    },
    /// One requested RB dimension is empty, exceeds its cap, or exceeds
    /// the truth-space dimension.
    InvalidRbDimension {
        /// Zero-based rung position.
        rung: usize,
        /// Requested basis dimension.
        dimension: usize,
        /// Maximum admitted at this rung.
        maximum: usize,
    },
    /// Rungs must become strictly coarser as their index increases.
    NonDecreasingFidelity {
        /// Zero-based position of the offending rung.
        rung: usize,
        /// Previous (finer) requested dimension.
        previous: usize,
        /// Current requested dimension.
        current: usize,
    },
    /// Orthogonalization retained a rung that is not strictly coarser
    /// than the prior rung, despite decreasing requested dimensions.
    NonDecreasingRetainedFidelity {
        /// Zero-based position of the offending rung.
        rung: usize,
        /// Prior retained basis dimension.
        previous: usize,
        /// Current retained basis dimension.
        current: usize,
        /// Current requested basis dimension.
        requested: usize,
    },
    /// Floating-point interpolation could not represent a strictly
    /// increasing grid.
    UnrepresentableGrid {
        /// Grid identity.
        what: &'static str,
        /// Index of the point that failed strict ordering.
        index: usize,
        /// Previous representable point.
        previous: f64,
        /// Candidate next point.
        next: f64,
    },
    /// Checked integer arithmetic overflowed while sizing work.
    BudgetArithmeticOverflow {
        /// Budget being computed.
        resource: &'static str,
    },
    /// A preflighted training resource exceeds its hard cap.
    BudgetExceeded {
        /// Resource identity.
        resource: &'static str,
        /// Requested amount.
        requested: usize,
        /// Configured cap.
        limit: usize,
    },
    /// A coverage axis is empty.
    EmptyCoverageBattery {
        /// Empty axis (`mu` or `tolerance`).
        axis: &'static str,
    },
    /// A coverage axis exceeds its independent cap.
    CoverageAxisTooLarge {
        /// Oversized axis.
        axis: &'static str,
        /// Requested count.
        requested: usize,
        /// Configured cap.
        limit: usize,
    },
    /// The coverage Cartesian product exceeds its cap.
    CoverageProductTooLarge {
        /// Requested query count.
        requested: usize,
        /// Configured cap.
        limit: usize,
    },
    /// Gram–Schmidt retained no basis vector (degenerate snapshots).
    EmptyBasis,
    /// A dense solve hit a zero or non-finite pivot.
    SingularSystem {
        /// Pivot column.
        column: usize,
    },
    /// A derived quantity (bound, QoI, solution entry) left the finite
    /// domain.
    NonFiniteDerived {
        /// What was being derived.
        what: &'static str,
    },
    /// A quantity that must be nonnegative was computed as negative.
    NegativeEnergy {
        /// Energy identity.
        what: &'static str,
        /// Rejected value.
        value: f64,
    },
    /// An emitted evidence color violated the shared payload grammar.
    InvalidEvidencePayload {
        /// Shared validator diagnostic.
        reason: String,
    },
    /// A production coverage run failed outside ordinary cooperative
    /// cancellation after all launched workers were drained.
    CoverageExecutionFailed {
        /// Structured executor or admission diagnostic.
        reason: String,
    },
}

impl fmt::Display for SurrogateError {
    #[allow(clippy::too_many_lines)] // exhaustive, one-to-one structured refusal rendering
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            SurrogateError::InvalidDimension { n } => {
                write!(f, "full-order dimension {n} outside 1..={MAX_TRUTH_NODES}")
            }
            SurrogateError::InvalidRange { lo, hi } => {
                write!(f, "mu range [{lo}, {hi}] must be finite with lo < hi")
            }
            SurrogateError::InvalidCoercivity { mu } => write!(
                f,
                "mu = {mu} breaks coercivity (need finite mu with 1 + mu > 0)"
            ),
            SurrogateError::OutOfRange { mu, range } => write!(
                f,
                "mu = {mu} outside the rung's declared training range [{}, {}] — ladder \
                 estimators are never extrapolated",
                range.0, range.1
            ),
            SurrogateError::InvalidTolerance { tol } => {
                write!(f, "tolerance {tol} must be finite and positive")
            }
            SurrogateError::InvalidGeometry { what, got } => {
                write!(f, "{what} = {got} is outside its allowed bounds")
            }
            SurrogateError::InvalidVectorLength {
                what,
                expected,
                got,
            } => write!(
                f,
                "{what} length {got} does not match required length {expected}"
            ),
            SurrogateError::NonFiniteInput { what, index } => {
                write!(f, "{what} contains a non-finite value at index {index}")
            }
            SurrogateError::InvalidLevel { requested, top } => {
                write!(
                    f,
                    "ladder level {requested} exceeds declared top level {top}"
                )
            }
            SurrogateError::FamilyMismatch { what } => {
                write!(f, "{what} does not match the ladder family identity")
            }
            SurrogateError::TooManyRbLevels { got, max } => {
                write!(f, "RB rung count {got} exceeds the cap {max}")
            }
            SurrogateError::InvalidRbDimension {
                rung,
                dimension,
                maximum,
            } => write!(
                f,
                "RB rung {rung} dimension {dimension} is outside 1..={maximum}"
            ),
            SurrogateError::NonDecreasingFidelity {
                rung,
                previous,
                current,
            } => write!(
                f,
                "RB rung {rung} dimension {current} must be strictly below prior dimension {previous}"
            ),
            SurrogateError::NonDecreasingRetainedFidelity {
                rung,
                previous,
                current,
                requested,
            } => write!(
                f,
                "RB rung {rung} requested dimension {requested} retained {current}, which is not \
                 strictly below the prior retained dimension {previous}"
            ),
            SurrogateError::UnrepresentableGrid {
                what,
                index,
                previous,
                next,
            } => write!(
                f,
                "{what} point {index} = {next} is not representably above prior point {previous}"
            ),
            SurrogateError::BudgetArithmeticOverflow { resource } => {
                write!(f, "{resource} budget arithmetic overflowed")
            }
            SurrogateError::BudgetExceeded {
                resource,
                requested,
                limit,
            } => write!(f, "{resource} request {requested} exceeds hard cap {limit}"),
            SurrogateError::EmptyCoverageBattery { axis } => {
                write!(f, "RB coverage {axis} axis must not be empty")
            }
            SurrogateError::CoverageAxisTooLarge {
                axis,
                requested,
                limit,
            } => write!(
                f,
                "RB coverage {axis} axis count {requested} exceeds cap {limit}"
            ),
            SurrogateError::CoverageProductTooLarge { requested, limit } => write!(
                f,
                "RB coverage Cartesian product {requested} exceeds cap {limit}"
            ),
            SurrogateError::EmptyBasis => {
                write!(
                    f,
                    "Gram–Schmidt retained no basis vector (degenerate snapshots)"
                )
            }
            SurrogateError::SingularSystem { column } => {
                write!(
                    f,
                    "linear solve hit an inadmissible/non-finite pivot at column {column}"
                )
            }
            SurrogateError::NonFiniteDerived { what } => {
                write!(
                    f,
                    "{what} left the finite domain — refusing to color a NaN/Inf"
                )
            }
            SurrogateError::NegativeEnergy { what, value } => {
                write!(f, "{what} computed inadmissible negative energy {value}")
            }
            SurrogateError::InvalidEvidencePayload { reason } => {
                write!(f, "ladder evidence payload is invalid: {reason}")
            }
            SurrogateError::CoverageExecutionFailed { reason } => {
                write!(f, "RB coverage execution failed after drain: {reason}")
            }
        }
    }
}

impl core::error::Error for SurrogateError {}

/// Internal arithmetic result that keeps cooperative cancellation distinct
/// from a scientific/numerical refusal. Public synchronous wrappers use the
/// uninhabited [`Infallible`] cancellation type; the scoped coverage kernel
/// uses [`Cancelled`].
#[derive(Debug)]
enum PolledError<C> {
    Cancelled(C),
    Numerical(SurrogateError),
}

impl<C> From<SurrogateError> for PolledError<C> {
    fn from(error: SurrogateError) -> Self {
        Self::Numerical(error)
    }
}

/// Sealed polling seam shared by the compatibility arithmetic and the
/// production coverage kernel. Implementations are private so callers cannot
/// weaken the fixed logical checkpoint protocol.
trait WorkPoll {
    type Cancel;

    fn work(&mut self, units: u64) -> Result<(), Self::Cancel>;

    fn checkpoint(&mut self) -> Result<(), Self::Cancel> {
        self.work(0)
    }
}

struct NoWorkPoll;

impl WorkPoll for NoWorkPoll {
    type Cancel = Infallible;

    fn work(&mut self, _units: u64) -> Result<(), Self::Cancel> {
        Ok(())
    }
}

fn poll_work<P: WorkPoll>(poll: &mut P, units: u64) -> Result<(), PolledError<P::Cancel>> {
    poll.work(units).map_err(PolledError::Cancelled)
}

fn poll_checkpoint<P: WorkPoll>(poll: &mut P) -> Result<(), PolledError<P::Cancel>> {
    poll.checkpoint().map_err(PolledError::Cancelled)
}

fn finish_unpolled<T>(
    result: Result<T, PolledError<Infallible>>,
) -> Result<T, SurrogateError> {
    match result {
        Ok(value) => Ok(value),
        Err(PolledError::Numerical(error)) => Err(error),
        Err(PolledError::Cancelled(never)) => match never {},
    }
}

fn filled_vec_polled<P: WorkPoll>(
    len: usize,
    value: f64,
    poll: &mut P,
) -> Result<Vec<f64>, PolledError<P::Cancel>> {
    poll_checkpoint(poll)?;
    let mut values = Vec::with_capacity(len);
    while values.len() < len {
        let old_len = values.len();
        let next_len = old_len
            .saturating_add(RB_COVERAGE_CHECKPOINT_WORK_UNITS as usize)
            .min(len);
        values.resize(next_len, value);
        poll_work(poll, (next_len - old_len) as u64)?;
    }
    poll_checkpoint(poll)?;
    Ok(values)
}

fn clone_slice_polled<P: WorkPoll>(
    values: &[f64],
    poll: &mut P,
) -> Result<Vec<f64>, PolledError<P::Cancel>> {
    poll_checkpoint(poll)?;
    let mut clone = Vec::with_capacity(values.len());
    for chunk in values.chunks(RB_COVERAGE_CHECKPOINT_WORK_UNITS as usize) {
        clone.extend_from_slice(chunk);
        poll_work(poll, chunk.len() as u64)?;
    }
    poll_checkpoint(poll)?;
    Ok(clone)
}

fn require_coercive(mu: f64) -> Result<(), SurrogateError> {
    if !mu.is_finite() || 1.0 + mu <= 0.0 {
        return Err(SurrogateError::InvalidCoercivity { mu });
    }
    Ok(())
}

fn require_range(lo: f64, hi: f64) -> Result<(), SurrogateError> {
    if !(lo.is_finite() && hi.is_finite() && lo < hi) {
        return Err(SurrogateError::InvalidRange { lo, hi });
    }
    require_coercive(lo)?;
    require_coercive(hi)?;
    let width = hi - lo;
    if !width.is_finite() || width <= 0.0 {
        return Err(SurrogateError::InvalidRange { lo, hi });
    }
    let mid = f64::midpoint(lo, hi);
    if !(mid.is_finite() && lo < mid && mid < hi) {
        return Err(SurrogateError::UnrepresentableGrid {
            what: "mu range",
            index: 1,
            previous: lo,
            next: mid,
        });
    }
    Ok(())
}

fn require_query_mu(mu: f64, range: (f64, f64)) -> Result<(), SurrogateError> {
    require_coercive(mu)?;
    if mu < range.0 || mu > range.1 {
        return Err(SurrogateError::OutOfRange { mu, range });
    }
    Ok(())
}

fn validate_vector_polled<P: WorkPoll>(
    values: &[f64],
    expected: usize,
    what: &'static str,
    poll: &mut P,
) -> Result<(), PolledError<P::Cancel>> {
    if values.len() != expected {
        return Err(SurrogateError::InvalidVectorLength {
            what,
            expected,
            got: values.len(),
        }
        .into());
    }
    for (index, value) in values.iter().enumerate() {
        if !value.is_finite() {
            return Err(SurrogateError::NonFiniteInput { what, index }.into());
        }
        poll_work(poll, 1)?;
    }
    Ok(())
}

fn strictly_spaced_grid(
    range: (f64, f64),
    points: usize,
    what: &'static str,
) -> Result<Vec<f64>, SurrogateError> {
    require_range(range.0, range.1)?;
    let mut grid = Vec::with_capacity(points);
    for index in 0..points {
        let value = if points == 1 {
            f64::midpoint(range.0, range.1)
        } else if index == 0 {
            range.0
        } else if index + 1 == points {
            range.1
        } else {
            #[allow(clippy::cast_precision_loss)]
            let t = index as f64 / (points - 1) as f64;
            range.0 + (range.1 - range.0) * t
        };
        if !value.is_finite() {
            return Err(SurrogateError::NonFiniteDerived { what });
        }
        if let Some(&previous) = grid.last()
            && value <= previous
        {
            return Err(SurrogateError::UnrepresentableGrid {
                what,
                index,
                previous,
                next: value,
            });
        }
        grid.push(value);
    }
    Ok(grid)
}

fn strictly_spaced_probe_grid(
    range: (f64, f64),
    points: usize,
) -> Result<Vec<f64>, SurrogateError> {
    require_range(range.0, range.1)?;
    let denominator = points
        .checked_mul(2)
        .ok_or(SurrogateError::BudgetArithmeticOverflow {
            resource: "concept probe grid",
        })?;
    let mut grid = Vec::with_capacity(points);
    for index in 0..points {
        #[allow(clippy::cast_precision_loss)]
        let t = (2 * index + 1) as f64 / denominator as f64;
        let value = range.0 + (range.1 - range.0) * t;
        let previous = grid.last().copied().unwrap_or(range.0);
        if !(value.is_finite() && previous < value && value < range.1) {
            return Err(SurrogateError::UnrepresentableGrid {
                what: "concept probe grid",
                index,
                previous,
                next: value,
            });
        }
        grid.push(value);
    }
    Ok(grid)
}

fn estimated_color(estimator: &'static str, dispersion: f64) -> Result<Color, SurrogateError> {
    let color = Color::Estimated {
        estimator: estimator.to_string(),
        dispersion,
    };
    validate_color_payload(&color).map_err(|error| SurrogateError::InvalidEvidencePayload {
        reason: error.to_string(),
    })?;
    Ok(color)
}

/// Immutable identity binding all rungs to one truth space and μ range.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct LadderFamily {
    truth_nodes: usize,
    mu_lo_bits: u64,
    mu_hi_bits: u64,
}

impl LadderFamily {
    fn new(truth: &TruthModel, range: (f64, f64)) -> Self {
        LadderFamily {
            truth_nodes: truth.n(),
            mu_lo_bits: range.0.to_bits(),
            mu_hi_bits: range.1.to_bits(),
        }
    }

    /// Full-order interior-node count.
    #[must_use]
    pub fn truth_nodes(self) -> usize {
        self.truth_nodes
    }

    /// Exact floating-point training range encoded by this identity.
    #[must_use]
    pub fn mu_range(self) -> (f64, f64) {
        (
            f64::from_bits(self.mu_lo_bits),
            f64::from_bits(self.mu_hi_bits),
        )
    }
}

/// The full-order "truth" model: P1 finite elements on a uniform grid
/// for `−(a u′)′ = 1`, `u(0) = u(1) = 0`, `a(x;μ) = 1 + μ·χ_{[½,1]}`.
/// The FE model IS level 0's DECLARED semantics — discretization
/// honesty lives in the CONTRACT; its floating-point solve carries no
/// numeric certificate (bead y6yv), so level-0 answers are Estimated.
///
/// Dimensions are SEALED: a `TruthModel` exists only through
/// [`TruthModel::new`], which bounds `n` before any allocation.
#[derive(Debug, Clone)]
pub struct TruthModel {
    n: usize,
}

impl TruthModel {
    /// A bounded truth model.
    ///
    /// # Errors
    /// [`SurrogateError::InvalidDimension`] outside `1..=MAX_TRUTH_NODES`.
    pub fn new(n: usize) -> Result<TruthModel, SurrogateError> {
        if n == 0 || n > MAX_TRUTH_NODES {
            return Err(SurrogateError::InvalidDimension { n });
        }
        Ok(TruthModel { n })
    }

    /// Interior node count.
    #[must_use]
    pub fn n(&self) -> usize {
        self.n
    }

    /// Assemble and solve at parameter `mu` (Thomas algorithm).
    /// Returns interior nodal values.
    ///
    /// # Errors
    /// [`SurrogateError::InvalidCoercivity`] for non-finite or
    /// non-coercive `mu`; [`SurrogateError::NonFiniteDerived`] if the
    /// solve leaves the finite domain.
    pub fn solve(&self, mu: f64) -> Result<Vec<f64>, SurrogateError> {
        let mut poll = NoWorkPoll;
        finish_unpolled(self.solve_polled(mu, &mut poll))
    }

    fn solve_polled<P: WorkPoll>(
        &self,
        mu: f64,
        poll: &mut P,
    ) -> Result<Vec<f64>, PolledError<P::Cancel>> {
        require_coercive(mu)?;
        let n = self.n;
        #[allow(clippy::cast_precision_loss)]
        let h = 1.0 / (n as f64 + 1.0);
        let mut diag = filled_vec_polled(n, 0.0, poll)?;
        let mut off = filled_vec_polled(n.saturating_sub(1), 0.0, poll)?;
        let rhs = filled_vec_polled(n, h, poll)?;
        for i in 0..=n {
            #[allow(clippy::cast_precision_loss)]
            let mid = (i as f64 + 0.5) * h;
            let coefficient = if mid >= 0.5 { 1.0 + mu } else { 1.0 };
            let a = coefficient / h;
            if !(mid.is_finite() && coefficient.is_finite() && coefficient > 0.0 && a.is_finite()) {
                return Err(SurrogateError::NonFiniteDerived {
                    what: "full-order stiffness assembly",
                }
                .into());
            }
            if i < n {
                diag[i] += a;
                if !diag[i].is_finite() {
                    return Err(SurrogateError::NonFiniteDerived {
                        what: "full-order stiffness diagonal",
                    }
                    .into());
                }
            }
            if i > 0 {
                diag[i - 1] += a;
                if !diag[i - 1].is_finite() {
                    return Err(SurrogateError::NonFiniteDerived {
                        what: "full-order stiffness diagonal",
                    }
                    .into());
                }
            }
            if i > 0 && i < n {
                off[i - 1] -= a;
                if !off[i - 1].is_finite() {
                    return Err(SurrogateError::NonFiniteDerived {
                        what: "full-order stiffness off-diagonal",
                    }
                    .into());
                }
            }
            poll_work(poll, 1)?;
        }
        // Thomas solve (diagonally dominant for coercive a; pivots
        // stay positive, but the finiteness of the RESULT is still
        // checked — refusal beats a colored NaN).
        let mut c = clone_slice_polled(&off, poll)?;
        let mut d = rhs;
        let first_pivot = diag[0];
        if !first_pivot.is_finite() || first_pivot <= 0.0 {
            return Err(SurrogateError::SingularSystem { column: 0 }.into());
        }
        if let Some(c0) = c.first_mut() {
            *c0 /= first_pivot;
            if !c0.is_finite() {
                return Err(SurrogateError::NonFiniteDerived {
                    what: "full-order Thomas coefficient",
                }
                .into());
            }
        }
        d[0] /= first_pivot;
        if !d[0].is_finite() {
            return Err(SurrogateError::NonFiniteDerived {
                what: "full-order Thomas right-hand side",
            }
            .into());
        }
        for i in 1..n {
            let m = diag[i] - off[i - 1] * c[i - 1];
            if !m.is_finite() || m <= 0.0 {
                return Err(SurrogateError::SingularSystem { column: i }.into());
            }
            if i < n - 1 {
                c[i] = off[i] / m;
                if !c[i].is_finite() {
                    return Err(SurrogateError::NonFiniteDerived {
                        what: "full-order Thomas coefficient",
                    }
                    .into());
                }
            }
            d[i] = (d[i] - off[i - 1] * d[i - 1]) / m;
            if !d[i].is_finite() {
                return Err(SurrogateError::NonFiniteDerived {
                    what: "full-order Thomas right-hand side",
                }
                .into());
            }
            poll_work(poll, 1)?;
        }
        let mut u = d;
        for i in (0..n.saturating_sub(1)).rev() {
            u[i] -= c[i] * u[i + 1];
            if !u[i].is_finite() {
                return Err(SurrogateError::NonFiniteDerived {
                    what: "full-order back-substitution",
                }
                .into());
            }
            poll_work(poll, 1)?;
        }
        for value in &u {
            if !value.is_finite() {
                return Err(SurrogateError::NonFiniteDerived {
                    what: "full-order solution",
                }
                .into());
            }
            poll_work(poll, 1)?;
        }
        poll_checkpoint(poll)?;
        Ok(u)
    }

    /// The energy inner product `a(u, v; mu)`.
    ///
    /// # Errors
    /// Refuses invalid μ, shape mismatches, non-finite operands, and any
    /// non-finite derived term or accumulation.
    pub fn energy(&self, u: &[f64], v: &[f64], mu: f64) -> Result<f64, SurrogateError> {
        let mut poll = NoWorkPoll;
        finish_unpolled(self.energy_polled(u, v, mu, &mut poll))
    }

    fn energy_polled<P: WorkPoll>(
        &self,
        u: &[f64],
        v: &[f64],
        mu: f64,
        poll: &mut P,
    ) -> Result<f64, PolledError<P::Cancel>> {
        require_coercive(mu)?;
        let n = self.n;
        validate_vector_polled(u, n, "energy left operand", poll)?;
        validate_vector_polled(v, n, "energy right operand", poll)?;
        #[allow(clippy::cast_precision_loss)]
        let h = 1.0 / (n as f64 + 1.0);
        let mut acc = 0.0f64;
        for i in 0..=n {
            #[allow(clippy::cast_precision_loss)]
            let mid = (i as f64 + 0.5) * h;
            let a = if mid >= 0.5 { 1.0 + mu } else { 1.0 };
            let du = if i == 0 {
                u[0]
            } else if i == n {
                -u[n - 1]
            } else {
                u[i] - u[i - 1]
            };
            let dv = if i == 0 {
                v[0]
            } else if i == n {
                -v[n - 1]
            } else {
                v[i] - v[i - 1]
            };
            let term = a * du * dv / h;
            if !term.is_finite() {
                return Err(SurrogateError::NonFiniteDerived {
                    what: "energy element contribution",
                }
                .into());
            }
            acc += term;
            if !acc.is_finite() {
                return Err(SurrogateError::NonFiniteDerived {
                    what: "energy accumulation",
                }
                .into());
            }
            poll_work(poll, 1)?;
        }
        poll_checkpoint(poll)?;
        Ok(acc)
    }

    /// The compliance QoI `∫ f u = h·Σu` (f = 1).
    ///
    /// # Errors
    /// Refuses a shape mismatch, a non-finite input, or non-finite sum.
    pub fn compliance(&self, u: &[f64]) -> Result<f64, SurrogateError> {
        let mut poll = NoWorkPoll;
        finish_unpolled(self.compliance_polled(u, &mut poll))
    }

    fn compliance_polled<P: WorkPoll>(
        &self,
        u: &[f64],
        poll: &mut P,
    ) -> Result<f64, PolledError<P::Cancel>> {
        validate_vector_polled(u, self.n, "compliance operand", poll)?;
        #[allow(clippy::cast_precision_loss)]
        let h = 1.0 / (self.n as f64 + 1.0);
        let mut sum = 0.0;
        for value in u {
            sum += value;
            if !sum.is_finite() {
                return Err(SurrogateError::NonFiniteDerived {
                    what: "compliance accumulation",
                }
                .into());
            }
            poll_work(poll, 1)?;
        }
        let compliance = h * sum;
        if !compliance.is_finite() {
            return Err(SurrogateError::NonFiniteDerived {
                what: "compliance result",
            }
            .into());
        }
        poll_checkpoint(poll)?;
        Ok(compliance)
    }
}

/// One RB rung: an energy-orthonormal basis with the offline residual
/// machinery for the online dual-norm bound. The bound is honest
/// mathematics evaluated in f64 — its answers are ESTIMATED (module
/// docs).
#[derive(Debug, Clone)]
pub struct RbLevel {
    truth: TruthModel,
    basis: Vec<Vec<f64>>,
    mu_range: (f64, f64),
    family: LadderFamily,
}

impl RbLevel {
    /// OFFLINE: snapshots at `k` training parameters spread over
    /// `mu_range`, Gram–Schmidt in the reference (μ = mid) energy.
    ///
    /// # Errors
    /// Geometry/range refusals; [`SurrogateError::EmptyBasis`] when no
    /// snapshot survives orthonormalization.
    pub fn train(
        truth: &TruthModel,
        mu_range: (f64, f64),
        k: usize,
    ) -> Result<RbLevel, SurrogateError> {
        require_range(mu_range.0, mu_range.1)?;
        validate_training_plan(truth, &[k], false)?;
        let training_grid = strictly_spaced_grid(mu_range, k, "RB training grid")?;
        let mid = f64::midpoint(mu_range.0, mu_range.1);
        let mut basis: Vec<Vec<f64>> = Vec::with_capacity(k);
        for mu in training_grid {
            let mut snap = truth.solve(mu)?;
            for b in &basis {
                let proj = truth.energy(&snap, b, mid)?;
                for (s, bi) in snap.iter_mut().zip(b) {
                    *s -= proj * bi;
                    if !s.is_finite() {
                        return Err(SurrogateError::NonFiniteDerived {
                            what: "RB orthogonalization",
                        });
                    }
                }
            }
            let norm_squared = truth.energy(&snap, &snap, mid)?;
            if norm_squared < 0.0 {
                return Err(SurrogateError::NegativeEnergy {
                    what: "RB snapshot norm",
                    value: norm_squared,
                });
            }
            let norm = norm_squared.sqrt();
            if !norm.is_finite() {
                return Err(SurrogateError::NonFiniteDerived {
                    what: "RB snapshot norm",
                });
            }
            if norm > 1e-10 {
                for s in &mut snap {
                    *s /= norm;
                    if !s.is_finite() {
                        return Err(SurrogateError::NonFiniteDerived {
                            what: "RB basis normalization",
                        });
                    }
                }
                basis.push(snap);
            }
        }
        if basis.is_empty() {
            return Err(SurrogateError::EmptyBasis);
        }
        Ok(RbLevel {
            truth: truth.clone(),
            basis,
            mu_range,
            family: LadderFamily::new(truth, mu_range),
        })
    }

    /// Basis size.
    #[must_use]
    pub fn dim(&self) -> usize {
        self.basis.len()
    }

    /// Declared training/query range.
    #[must_use]
    pub fn mu_range(&self) -> (f64, f64) {
        self.mu_range
    }

    /// Immutable truth-space/range identity.
    #[must_use]
    pub fn family(&self) -> LadderFamily {
        self.family
    }

    /// ONLINE: Galerkin solve in the basis + the energy and QoI bounds
    /// (f64-evaluated — ESTIMATED authority). Returns (`u_rb` in nodal
    /// form, compliance, `energy_bound`, `qoi_bound`).
    ///
    /// # Errors
    /// [`SurrogateError::OutOfRange`] outside the declared training
    /// range (bounds are never extrapolated);
    /// [`SurrogateError::SingularSystem`] /
    /// [`SurrogateError::NonFiniteDerived`] from the dense solve.
    pub fn query(&self, mu: f64) -> Result<(Vec<f64>, f64, f64, f64), SurrogateError> {
        let mut poll = NoWorkPoll;
        finish_unpolled(self.query_polled(mu, &mut poll))
    }

    fn query_polled<P: WorkPoll>(
        &self,
        mu: f64,
        poll: &mut P,
    ) -> Result<(Vec<f64>, f64, f64, f64), PolledError<P::Cancel>> {
        require_query_mu(mu, self.mu_range)?;
        let k = self.basis.len();
        poll_checkpoint(poll)?;
        let mut a = Vec::with_capacity(k);
        for _ in 0..k {
            a.push(filled_vec_polled(k, 0.0, poll)?);
        }
        let mut b = filled_vec_polled(k, 0.0, poll)?;
        for i in 0..k {
            for (j, aij) in a[i].iter_mut().enumerate() {
                *aij = self
                    .truth
                    .energy_polled(&self.basis[j], &self.basis[i], mu, poll)?;
            }
            b[i] = self.truth.compliance_polled(&self.basis[i], poll)?;
            poll_checkpoint(poll)?;
        }
        let coef = solve_dense_polled(&mut a, &mut b, poll)?;
        let n = self.truth.n();
        let mut u = filled_vec_polled(n, 0.0, poll)?;
        for (c, basis_vec) in coef.iter().zip(&self.basis) {
            for (ui, bi) in u.iter_mut().zip(basis_vec) {
                *ui += c * bi;
                if !ui.is_finite() {
                    return Err(SurrogateError::NonFiniteDerived {
                        what: "RB solution reconstruction",
                    }
                    .into());
                }
                poll_work(poll, 1)?;
            }
        }
        let s_rb = self.truth.compliance_polled(&u, poll)?;
        let riesz = self.residual_riesz_polled(&u, mu, poll)?;
        let dual_energy = self.truth.energy_polled(&riesz, &riesz, 0.0, poll)?;
        if dual_energy < 0.0 {
            return Err(SurrogateError::NegativeEnergy {
                what: "residual dual norm",
                value: dual_energy,
            }
            .into());
        }
        let dual_norm = dual_energy.sqrt();
        let alpha_lb = 1.0f64.min(1.0 + mu);
        if !(dual_norm.is_finite() && alpha_lb.is_finite() && alpha_lb > 0.0) {
            return Err(SurrogateError::NonFiniteDerived {
                what: "RB residual norm or coercivity floor",
            }
            .into());
        }
        let energy_bound = dual_norm / alpha_lb.sqrt();
        // Exact Galerkin orthogonality would make the compliance error equal
        // the energy error squared. Floating elimination leaves a computable
        // defect r(u_rb) = f(u_rb) - a(u_rb,u_rb), which must be retained.
        let reduced_energy = self.truth.energy_polled(&u, &u, mu, poll)?;
        let galerkin_defect = (s_rb - reduced_energy).abs();
        let qoi_bound = energy_bound.mul_add(energy_bound, galerkin_defect);
        if !(energy_bound.is_finite()
            && galerkin_defect.is_finite()
            && qoi_bound.is_finite()
            && qoi_bound >= 0.0)
        {
            return Err(SurrogateError::NonFiniteDerived {
                what: "rb bound arithmetic",
            }
            .into());
        }
        poll_checkpoint(poll)?;
        Ok((u, s_rb, energy_bound, qoi_bound))
    }

    /// Solve the reference-Laplacian Riesz problem for the residual.
    fn residual_riesz_polled<P: WorkPoll>(
        &self,
        u_rb: &[f64],
        mu: f64,
        poll: &mut P,
    ) -> Result<Vec<f64>, PolledError<P::Cancel>> {
        require_query_mu(mu, self.mu_range)?;
        let n = self.truth.n();
        validate_vector_polled(u_rb, n, "RB residual state", poll)?;
        #[allow(clippy::cast_precision_loss)]
        let h = 1.0 / (n as f64 + 1.0);
        let mut rhs = filled_vec_polled(n, h, poll)?;
        for (i, ri) in rhs.iter_mut().enumerate() {
            let mut acc = 0.0f64;
            for e in [i, i + 1] {
                #[allow(clippy::cast_precision_loss)]
                let mid = (e as f64 + 0.5) * h;
                let a = if mid >= 0.5 { 1.0 + mu } else { 1.0 };
                let du = if e == 0 {
                    u_rb[0]
                } else if e == n {
                    -u_rb[n - 1]
                } else {
                    u_rb[e] - u_rb[e - 1]
                };
                let dphi = if e == i { 1.0 } else { -1.0 };
                let term = a * du * dphi / h;
                if !term.is_finite() {
                    return Err(SurrogateError::NonFiniteDerived {
                        what: "RB residual element contribution",
                    }
                    .into());
                }
                acc += term;
                if !acc.is_finite() {
                    return Err(SurrogateError::NonFiniteDerived {
                        what: "RB residual accumulation",
                    }
                    .into());
                }
            }
            *ri -= acc;
            if !ri.is_finite() {
                return Err(SurrogateError::NonFiniteDerived {
                    what: "RB residual right-hand side",
                }
                .into());
            }
            poll_work(poll, 1)?;
        }
        let diag = filled_vec_polled(n, 2.0 / h, poll)?;
        let off = filled_vec_polled(n.saturating_sub(1), -1.0 / h, poll)?;
        let mut c = clone_slice_polled(&off, poll)?;
        let mut d = rhs;
        let first_pivot = diag[0];
        if !first_pivot.is_finite() || first_pivot <= 0.0 {
            return Err(SurrogateError::SingularSystem { column: 0 }.into());
        }
        if let Some(c0) = c.first_mut() {
            *c0 /= first_pivot;
            if !c0.is_finite() {
                return Err(SurrogateError::NonFiniteDerived {
                    what: "Riesz Thomas coefficient",
                }
                .into());
            }
        }
        d[0] /= first_pivot;
        if !d[0].is_finite() {
            return Err(SurrogateError::NonFiniteDerived {
                what: "Riesz Thomas right-hand side",
            }
            .into());
        }
        for i in 1..n {
            let m = diag[i] - off[i - 1] * c[i - 1];
            if !m.is_finite() || m <= 0.0 {
                return Err(SurrogateError::SingularSystem { column: i }.into());
            }
            if i < n - 1 {
                c[i] = off[i] / m;
                if !c[i].is_finite() {
                    return Err(SurrogateError::NonFiniteDerived {
                        what: "Riesz Thomas coefficient",
                    }
                    .into());
                }
            }
            d[i] = (d[i] - off[i - 1] * d[i - 1]) / m;
            if !d[i].is_finite() {
                return Err(SurrogateError::NonFiniteDerived {
                    what: "Riesz Thomas right-hand side",
                }
                .into());
            }
            poll_work(poll, 1)?;
        }
        for i in (0..n.saturating_sub(1)).rev() {
            d[i] -= c[i] * d[i + 1];
            if !d[i].is_finite() {
                return Err(SurrogateError::NonFiniteDerived {
                    what: "Riesz back-substitution",
                }
                .into());
            }
            poll_work(poll, 1)?;
        }
        poll_checkpoint(poll)?;
        Ok(d)
    }
}

fn solve_dense(a: &mut [Vec<f64>], b: &mut [f64]) -> Result<Vec<f64>, SurrogateError> {
    let mut poll = NoWorkPoll;
    finish_unpolled(solve_dense_polled(a, b, &mut poll))
}

fn solve_dense_polled<P: WorkPoll>(
    a: &mut [Vec<f64>],
    b: &mut [f64],
    poll: &mut P,
) -> Result<Vec<f64>, PolledError<P::Cancel>> {
    let n = b.len();
    if n == 0 {
        return Err(SurrogateError::InvalidGeometry {
            what: "dense system dimension",
            got: 0,
        }
        .into());
    }
    if a.len() != n {
        return Err(SurrogateError::InvalidVectorLength {
            what: "dense matrix row count",
            expected: n,
            got: a.len(),
        }
        .into());
    }
    validate_vector_polled(b, n, "dense right-hand side", poll)?;
    for row in a.iter() {
        validate_vector_polled(row, n, "dense matrix row", poll)?;
    }
    for col in 0..n {
        let mut piv = col;
        for r in col + 1..n {
            if a[r][col].abs() > a[piv][col].abs() {
                piv = r;
            }
            poll_work(poll, 1)?;
        }
        a.swap(col, piv);
        b.swap(col, piv);
        let p = a[col][col];
        if p == 0.0 || !p.is_finite() {
            return Err(SurrogateError::SingularSystem { column: col }.into());
        }
        let pivot_rhs = b[col];
        for (r, rhs) in b.iter_mut().enumerate().skip(col + 1) {
            let (rows_before, rows_from_r) = a.split_at_mut(r);
            let pivot_row = &rows_before[col];
            let row = &mut rows_from_r[0];
            let f = row[col] / p;
            if !f.is_finite() {
                return Err(SurrogateError::NonFiniteDerived {
                    what: "dense elimination factor",
                }
                .into());
            }
            for (arc, pc) in row[col..].iter_mut().zip(&pivot_row[col..]) {
                *arc -= f * pc;
                if !arc.is_finite() {
                    return Err(SurrogateError::NonFiniteDerived {
                        what: "dense elimination matrix update",
                    }
                    .into());
                }
                poll_work(poll, 1)?;
            }
            *rhs -= f * pivot_rhs;
            if !rhs.is_finite() {
                return Err(SurrogateError::NonFiniteDerived {
                    what: "dense elimination right-hand side",
                }
                .into());
            }
            poll_work(poll, 1)?;
        }
        poll_checkpoint(poll)?;
    }
    let mut x = filled_vec_polled(n, 0.0, poll)?;
    for r in (0..n).rev() {
        let mut acc = b[r];
        for c in r + 1..n {
            acc -= a[r][c] * x[c];
            if !acc.is_finite() {
                return Err(SurrogateError::NonFiniteDerived {
                    what: "dense back-substitution accumulation",
                }
                .into());
            }
            poll_work(poll, 1)?;
        }
        let pivot = a[r][r];
        if pivot == 0.0 || !pivot.is_finite() {
            return Err(SurrogateError::SingularSystem { column: r }.into());
        }
        x[r] = acc / pivot;
        if !x[r].is_finite() {
            return Err(SurrogateError::NonFiniteDerived {
                what: "dense back-substitution",
            }
            .into());
        }
        poll_checkpoint(poll)?;
    }
    poll_checkpoint(poll)?;
    Ok(x)
}

/// A NON-RB concept level: QoI lookup by linear interpolation over a
/// training grid, total dispersion calibrated by CROSS-RUNG PROBES as
/// `|concept − lower RB| + lower RB QoI estimator`.
#[derive(Debug, Clone)]
pub struct ConceptLevel {
    grid: Vec<(f64, f64)>,
    dispersion: f64,
    mu_range: (f64, f64),
    family: LadderFamily,
}

impl ConceptLevel {
    /// Build from a training grid and calibrate against a lower rung.
    ///
    /// # Errors
    /// Geometry refusals (`grid_points < 2` would divide by zero;
    /// `probes < 1` calibrates nothing; both capped) and every rung
    /// query refusal.
    pub fn train(
        rb: &RbLevel,
        grid_points: usize,
        probes: usize,
    ) -> Result<ConceptLevel, SurrogateError> {
        if !(2..=MAX_CONCEPT_POINTS).contains(&grid_points) {
            return Err(SurrogateError::InvalidGeometry {
                what: "concept grid points",
                got: grid_points,
            });
        }
        if probes == 0 || probes > MAX_CONCEPT_POINTS {
            return Err(SurrogateError::InvalidGeometry {
                what: "concept probes",
                got: probes,
            });
        }
        let evaluations =
            grid_points
                .checked_add(probes)
                .ok_or(SurrogateError::BudgetArithmeticOverflow {
                    resource: "concept training work",
                })?;
        let work_units = checked_budget_product(
            "concept training work",
            &[rb.truth.n(), rb.dim(), rb.dim(), evaluations],
        )?;
        if work_units > MAX_TRAINING_WORK_UNITS {
            return Err(SurrogateError::BudgetExceeded {
                resource: "concept training work",
                requested: work_units,
                limit: MAX_TRAINING_WORK_UNITS,
            });
        }
        let range = rb.mu_range();
        let training_grid = strictly_spaced_grid(range, grid_points, "concept training grid")?;
        let probe_grid = strictly_spaced_probe_grid(range, probes)?;
        let mut grid = Vec::with_capacity(grid_points);
        for mu in training_grid {
            grid.push((mu, rb.query(mu)?.1));
        }
        let mut level = ConceptLevel {
            grid,
            dispersion: f64::INFINITY,
            mu_range: range,
            family: rb.family(),
        };
        let mut disp = 0.0f64;
        for mu in probe_grid {
            let (_, rb_value, _, rb_qoi_estimator) = rb.query(mu)?;
            let discrepancy = (level.lookup(mu)? - rb_value).abs();
            if !discrepancy.is_finite() {
                return Err(SurrogateError::NonFiniteDerived {
                    what: "concept discrepancy probe",
                });
            }
            let total_dispersion = discrepancy + rb_qoi_estimator;
            if !(rb_qoi_estimator.is_finite()
                && rb_qoi_estimator >= 0.0
                && total_dispersion.is_finite())
            {
                return Err(SurrogateError::NonFiniteDerived {
                    what: "concept total dispersion",
                });
            }
            disp = disp.max(total_dispersion);
        }
        level.dispersion = disp;
        Ok(level)
    }

    /// Linear interpolation lookup (grid length ≥ 2 by construction).
    ///
    /// # Errors
    /// Refuses non-finite, non-coercive, or out-of-range queries before
    /// inspecting the lookup grid; also refuses every non-finite
    /// interpolation result.
    pub fn lookup(&self, mu: f64) -> Result<f64, SurrogateError> {
        require_query_mu(mu, self.mu_range)?;
        let g = &self.grid;
        let pos = g.partition_point(|(m, _)| *m < mu).clamp(1, g.len() - 1);
        let (m0, v0) = g[pos - 1];
        let (m1, v1) = g[pos];
        let width = m1 - m0;
        if !(width.is_finite() && width > 0.0) {
            return Err(SurrogateError::UnrepresentableGrid {
                what: "concept lookup grid",
                index: pos,
                previous: m0,
                next: m1,
            });
        }
        let t = (mu - m0) / width;
        if !(t.is_finite() && (0.0..=1.0).contains(&t)) {
            return Err(SurrogateError::NonFiniteDerived {
                what: "concept interpolation coordinate",
            });
        }
        let value = v0 + t * (v1 - v0);
        if !value.is_finite() {
            return Err(SurrogateError::NonFiniteDerived {
                what: "concept interpolation result",
            });
        }
        Ok(value)
    }

    /// Calibrated maximum cross-rung discrepancy over the probe grid.
    #[must_use]
    pub fn dispersion(&self) -> f64 {
        self.dispersion
    }

    /// Declared training/query range.
    #[must_use]
    pub fn mu_range(&self) -> (f64, f64) {
        self.mu_range
    }

    /// Immutable truth-space/range identity inherited from the RB rung.
    #[must_use]
    pub fn family(&self) -> LadderFamily {
        self.family
    }
}

/// A ladder answer: the value, its (always ESTIMATED — module docs)
/// color, and the drill-down forensics.
#[derive(Debug, Clone)]
pub struct LadderAnswer {
    value: f64,
    color: Color,
    level_used: usize,
    leaks: Vec<usize>,
}

impl LadderAnswer {
    /// QoI value.
    #[must_use]
    pub fn value(&self) -> f64 {
        self.value
    }

    /// Evidence color. Every color emitted by this module is Estimated;
    /// level 0 carries infinite dispersion because it has no numerical
    /// enclosure.
    #[must_use]
    pub fn color(&self) -> &Color {
        &self.color
    }

    /// Level that ultimately answered (0 = full order).
    #[must_use]
    pub fn level_used(&self) -> usize {
        self.level_used
    }

    /// Levels whose estimator exceeded tolerance during descent.
    #[must_use]
    pub fn leaks(&self) -> &[usize] {
        &self.leaks
    }
}

/// The assembled ladder: level 0 = full order (declared truth), higher
/// = more abstract. `at_level(k)` starts at rung k and DESCENDS
/// AUTOMATICALLY on leak — invisible until it leaks.
#[derive(Debug, Clone)]
pub struct Ladder {
    truth: TruthModel,
    rb_levels: Vec<RbLevel>,
    concept: Option<ConceptLevel>,
    family: LadderFamily,
}

/// A view of the ladder starting at a rung.
pub struct LevelView<'a> {
    ladder: &'a Ladder,
    start: usize,
}

fn checked_budget_product(
    resource: &'static str,
    factors: &[usize],
) -> Result<usize, SurrogateError> {
    factors.iter().try_fold(1usize, |product, factor| {
        product
            .checked_mul(*factor)
            .ok_or(SurrogateError::BudgetArithmeticOverflow { resource })
    })
}

#[allow(clippy::too_many_lines)] // one atomic preflight; no allocation may precede any check
fn validate_training_plan(
    truth: &TruthModel,
    rb_dims: &[usize],
    concept: bool,
) -> Result<(), SurrogateError> {
    if rb_dims.is_empty() {
        return Err(SurrogateError::InvalidGeometry {
            what: "rb rung count",
            got: 0,
        });
    }
    if rb_dims.len() > MAX_RB_LEVELS {
        return Err(SurrogateError::TooManyRbLevels {
            got: rb_dims.len(),
            max: MAX_RB_LEVELS,
        });
    }

    let maximum = MAX_RB_DIM.min(truth.n());
    let mut previous = None;
    let mut memory_bytes = 0usize;
    let mut largest_dimension = 0usize;
    let mut work_units = 0usize;
    for (rung, &dimension) in rb_dims.iter().enumerate() {
        if dimension == 0 || dimension > maximum {
            return Err(SurrogateError::InvalidRbDimension {
                rung,
                dimension,
                maximum,
            });
        }
        if let Some(previous) = previous
            && dimension >= previous
        {
            return Err(SurrogateError::NonDecreasingFidelity {
                rung,
                previous,
                current: dimension,
            });
        }
        previous = Some(dimension);
        largest_dimension = largest_dimension.max(dimension);

        let rung_bytes = checked_budget_product(
            "aggregate ladder training memory",
            &[truth.n(), dimension, size_of::<f64>()],
        )?;
        memory_bytes = memory_bytes.checked_add(rung_bytes).ok_or(
            SurrogateError::BudgetArithmeticOverflow {
                resource: "aggregate ladder training memory",
            },
        )?;
        let rung_work =
            checked_budget_product("RB training work", &[truth.n(), dimension, dimension])?;
        work_units =
            work_units
                .checked_add(rung_work)
                .ok_or(SurrogateError::BudgetArithmeticOverflow {
                    resource: "RB training work",
                })?;
    }

    // Conservative peak scratch: one state plus the tridiagonal solve's
    // diagonal/off-diagonal/right-hand-side/work vectors. Dense online
    // training queries additionally hold one k×k matrix.
    let scratch_bytes = checked_budget_product(
        "aggregate ladder training memory",
        &[truth.n(), 6, size_of::<f64>()],
    )?;
    let dense_bytes = checked_budget_product(
        "aggregate ladder training memory",
        &[largest_dimension, largest_dimension, size_of::<f64>()],
    )?;
    memory_bytes = memory_bytes
        .checked_add(scratch_bytes)
        .and_then(|bytes| bytes.checked_add(dense_bytes))
        .ok_or(SurrogateError::BudgetArithmeticOverflow {
            resource: "aggregate ladder training memory",
        })?;

    if concept {
        let coarsest = *rb_dims.last().ok_or(SurrogateError::InvalidGeometry {
            what: "rb rung count",
            got: 0,
        })?;
        let concept_queries = DEFAULT_CONCEPT_GRID_POINTS
            .checked_add(DEFAULT_CONCEPT_PROBES)
            .ok_or(SurrogateError::BudgetArithmeticOverflow {
                resource: "concept training work",
            })?;
        let concept_work = checked_budget_product(
            "concept training work",
            &[truth.n(), coarsest, coarsest, concept_queries],
        )?;
        work_units = work_units.checked_add(concept_work).ok_or(
            SurrogateError::BudgetArithmeticOverflow {
                resource: "aggregate ladder training work",
            },
        )?;
    }

    if memory_bytes > MAX_TRAINING_BYTES {
        return Err(SurrogateError::BudgetExceeded {
            resource: "aggregate ladder training memory",
            requested: memory_bytes,
            limit: MAX_TRAINING_BYTES,
        });
    }
    if work_units > MAX_TRAINING_WORK_UNITS {
        return Err(SurrogateError::BudgetExceeded {
            resource: "aggregate ladder training work",
            requested: work_units,
            limit: MAX_TRAINING_WORK_UNITS,
        });
    }
    Ok(())
}

impl Ladder {
    /// Build: truth + RB rungs of the given basis sizes (descending
    /// fidelity) + a concept rung calibrated against the last RB rung.
    ///
    /// # Errors
    /// Every constructor/training refusal, plus
    /// [`SurrogateError::InvalidGeometry`] for an empty `rb_dims`.
    pub fn build(
        n: usize,
        mu_range: (f64, f64),
        rb_dims: &[usize],
        concept: bool,
    ) -> Result<Ladder, SurrogateError> {
        require_range(mu_range.0, mu_range.1)?;
        let truth = TruthModel::new(n)?;
        validate_training_plan(&truth, rb_dims, concept)?;
        let family = LadderFamily::new(&truth, mu_range);
        let mut rb_levels = Vec::with_capacity(rb_dims.len());
        let mut previous_retained = None;
        for (rung, &k) in rb_dims.iter().enumerate() {
            let level = RbLevel::train(&truth, mu_range, k)?;
            if level.family() != family {
                return Err(SurrogateError::FamilyMismatch { what: "RB rung" });
            }
            if let Some(previous) = previous_retained
                && level.dim() >= previous
            {
                return Err(SurrogateError::NonDecreasingRetainedFidelity {
                    rung,
                    previous,
                    current: level.dim(),
                    requested: k,
                });
            }
            previous_retained = Some(level.dim());
            rb_levels.push(level);
        }
        let concept = if concept {
            let rb = rb_levels.last().ok_or(SurrogateError::InvalidGeometry {
                what: "rb rung count",
                got: 0,
            })?;
            let level =
                ConceptLevel::train(rb, DEFAULT_CONCEPT_GRID_POINTS, DEFAULT_CONCEPT_PROBES)?;
            if level.family() != family {
                return Err(SurrogateError::FamilyMismatch {
                    what: "concept rung",
                });
            }
            Some(level)
        } else {
            None
        };
        Ok(Ladder {
            truth,
            rb_levels,
            concept,
            family,
        })
    }

    /// The view starting at rung `k` (0 = full order, 1.. = RB rungs,
    /// rb_levels.len()+1 = concept rung if present).
    ///
    /// # Errors
    /// [`SurrogateError::InvalidLevel`] when `k` exceeds [`Self::top`].
    pub fn at_level(&self, k: usize) -> Result<LevelView<'_>, SurrogateError> {
        if k > self.top() {
            return Err(SurrogateError::InvalidLevel {
                requested: k,
                top: self.top(),
            });
        }
        Ok(LevelView {
            ladder: self,
            start: k,
        })
    }

    /// Highest rung index.
    #[must_use]
    pub fn top(&self) -> usize {
        self.rb_levels.len() + usize::from(self.concept.is_some())
    }

    /// Number of RB rungs (excluding truth and concept).
    #[must_use]
    pub fn rb_level_count(&self) -> usize {
        self.rb_levels.len()
    }

    /// Immutable RB rungs, finest first.
    #[must_use]
    pub fn rb_levels(&self) -> &[RbLevel] {
        &self.rb_levels
    }

    /// Optional immutable concept rung.
    #[must_use]
    pub fn concept(&self) -> Option<&ConceptLevel> {
        self.concept.as_ref()
    }

    /// Immutable truth-space/range identity shared by every rung.
    #[must_use]
    pub fn family(&self) -> LadderFamily {
        self.family
    }

    /// Declared μ range shared by every ladder query.
    #[must_use]
    pub fn mu_range(&self) -> (f64, f64) {
        self.family.mu_range()
    }
}

impl LevelView<'_> {
    /// Query with AUTOMATIC ESTIMATOR-GUARDED DESCENT: each rung answers only
    /// if its estimator meets `tol`; a leaking rung is recorded and the query
    /// drills down. Level 0 answers unconditionally (it IS the declared
    /// truth semantics — and its color still says ESTIMATED, because a
    /// floating solve carries no certificate; bead y6yv).
    ///
    /// # Errors
    /// [`SurrogateError::InvalidTolerance`] /
    /// [`SurrogateError::InvalidCoercivity`] /
    /// [`SurrogateError::OutOfRange`] before any rung runs; rung solve
    /// refusals propagate.
    pub fn query(&self, mu: f64, tol: f64) -> Result<LadderAnswer, SurrogateError> {
        if !tol.is_finite() || tol <= 0.0 {
            return Err(SurrogateError::InvalidTolerance { tol });
        }
        require_query_mu(mu, self.ladder.mu_range())?;
        let mut leaks = Vec::new();
        let mut level = self.start;
        loop {
            if level > self.ladder.rb_levels.len() {
                let c = self
                    .ladder
                    .concept
                    .as_ref()
                    .ok_or(SurrogateError::InvalidLevel {
                        requested: level,
                        top: self.ladder.top(),
                    })?;
                if c.dispersion() > tol {
                    leaks.push(level);
                    level -= 1;
                    continue;
                }
                let concept_value = c.lookup(mu)?;
                let rb = self
                    .ladder
                    .rb_levels
                    .last()
                    .ok_or(SurrogateError::InvalidGeometry {
                        what: "concept lower RB rung",
                        got: 0,
                    })?;
                let (_, rb_value, _, rb_qoi_estimator) = rb.query(mu)?;
                let local_cross_rung = (concept_value - rb_value).abs();
                let local_dispersion = local_cross_rung + rb_qoi_estimator;
                if !(local_cross_rung.is_finite()
                    && rb_qoi_estimator.is_finite()
                    && rb_qoi_estimator >= 0.0
                    && local_dispersion.is_finite())
                {
                    return Err(SurrogateError::NonFiniteDerived {
                        what: "concept query-local dispersion",
                    });
                }
                let dispersion = c.dispersion().max(local_dispersion);
                if dispersion <= tol {
                    return Ok(LadderAnswer {
                        value: concept_value,
                        color: estimated_color(CONCEPT_ESTIMATOR_ID, dispersion)?,
                        level_used: level,
                        leaks,
                    });
                }
                leaks.push(level);
                level -= 1;
                continue;
            }
            if level >= 1 {
                let rb = &self.ladder.rb_levels[level - 1];
                let (_, s_rb, _, qoi_bound) = rb.query(mu)?;
                if qoi_bound <= tol {
                    return Ok(LadderAnswer {
                        value: s_rb,
                        // DEMOTED from Verified (bead y6yv): the
                        // a-posteriori bound is evaluated in
                        // round-to-nearest f64 with no outward rounding
                        // — honest mathematics, not an executable
                        // enclosure. The bound travels as dispersion.
                        color: estimated_color(RB_ESTIMATOR_ID, qoi_bound)?,
                        level_used: level,
                        leaks,
                    });
                }
                leaks.push(level);
                level -= 1;
                continue;
            }
            let u = self.ladder.truth.solve(mu)?;
            let s = self.ladder.truth.compliance(&u)?;
            return Ok(LadderAnswer {
                value: s,
                // DEMOTED from Verified{[s,s]} (bead y6yv): a plain
                // floating solve is the DECLARED truth semantics, not a
                // verified enclosure. Infinite dispersion is the
                // explicit no-spread-claim payload.
                color: estimated_color(TRUTH_ESTIMATOR_ID, f64::INFINITY)?,
                level_used: 0,
                leaks,
            });
        }
    }
}

fn conservative_coverage_query_work(ladder: &Ladder) -> Result<usize, SurrogateError> {
    let truth_nodes = ladder.truth.n();
    let mut work = checked_budget_product("RB coverage truth fallback", &[truth_nodes, 8])?;
    for rb in &ladder.rb_levels {
        let dimension = rb.dim();
        let matrix = checked_budget_product(
            "RB coverage reduced matrix",
            &[truth_nodes, dimension, dimension],
        )?;
        let projections = checked_budget_product(
            "RB coverage reduced projections",
            &[truth_nodes, dimension, 2],
        )?;
        let residual = checked_budget_product("RB coverage residual", &[truth_nodes, 8])?;
        let dense = checked_budget_product(
            "RB coverage dense solve",
            &[dimension, dimension, dimension],
        )?;
        work = work
            .checked_add(matrix)
            .and_then(|value| value.checked_add(projections))
            .and_then(|value| value.checked_add(residual))
            .and_then(|value| value.checked_add(dense))
            .ok_or(SurrogateError::BudgetArithmeticOverflow {
                resource: "RB coverage query work",
            })?;
    }
    Ok(work)
}

/// Exact, validated coverage battery bound to one immutable ladder family and
/// rung shape. Floating inputs are retained as raw bits so signed zero and
/// every finite representable value replay exactly; order and duplicates are
/// semantic and are never canonicalized away.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RbCoveragePlan {
    family: LadderFamily,
    rb_dimensions: Vec<usize>,
    mu_bits: Vec<u64>,
    tolerance_bits: Vec<u64>,
    strictest_tolerance_bits: u64,
    total_queries: usize,
    per_mu_work_units: usize,
    declared_work_units: usize,
}

impl RbCoveragePlan {
    /// Validate and retain an exact production coverage battery. The complete
    /// Cartesian shape and conservative work are admitted before input bytes
    /// are copied into the plan.
    ///
    /// # Errors
    /// Every bounded-battery, parameter, tolerance, arithmetic, or work-cap
    /// refusal used by [`rb_coverage`].
    pub fn new(ladder: &Ladder, mus: &[f64], tolerances: &[f64]) -> Result<Self, SurrogateError> {
        Self::prepare(
            ladder,
            mus,
            tolerances,
            MAX_COVERAGE_WORK_UNITS,
            "RB coverage total work",
        )
    }

    fn prepare(
        ladder: &Ladder,
        mus: &[f64],
        tolerances: &[f64],
        work_limit: usize,
        work_resource: &'static str,
    ) -> Result<Self, SurrogateError> {
        if mus.is_empty() {
            return Err(SurrogateError::EmptyCoverageBattery { axis: "mu" });
        }
        if tolerances.is_empty() {
            return Err(SurrogateError::EmptyCoverageBattery { axis: "tolerance" });
        }
        for (axis, count) in [("mu", mus.len()), ("tolerance", tolerances.len())] {
            if count > MAX_COVERAGE_AXIS_POINTS {
                return Err(SurrogateError::CoverageAxisTooLarge {
                    axis,
                    requested: count,
                    limit: MAX_COVERAGE_AXIS_POINTS,
                });
            }
        }
        let total_queries = mus.len().checked_mul(tolerances.len()).ok_or(
            SurrogateError::BudgetArithmeticOverflow {
                resource: "RB coverage Cartesian product",
            },
        )?;
        if total_queries > MAX_COVERAGE_QUERIES {
            return Err(SurrogateError::CoverageProductTooLarge {
                requested: total_queries,
                limit: MAX_COVERAGE_QUERIES,
            });
        }
        let per_mu_work_units = conservative_coverage_query_work(ladder)?;
        let solver_work = per_mu_work_units.checked_mul(mus.len()).ok_or(
            SurrogateError::BudgetArithmeticOverflow {
                resource: "RB coverage total work",
            },
        )?;
        let declared_work_units = solver_work.checked_add(total_queries).ok_or(
            SurrogateError::BudgetArithmeticOverflow {
                resource: "RB coverage total work",
            },
        )?;
        if declared_work_units > work_limit {
            return Err(SurrogateError::BudgetExceeded {
                resource: work_resource,
                requested: declared_work_units,
                limit: work_limit,
            });
        }
        for &mu in mus {
            require_query_mu(mu, ladder.mu_range())?;
        }
        for &tolerance in tolerances {
            if !tolerance.is_finite() || tolerance <= 0.0 {
                return Err(SurrogateError::InvalidTolerance { tol: tolerance });
            }
        }

        let strictest_tolerance = tolerances
            .iter()
            .copied()
            .fold(f64::INFINITY, f64::min);
        Ok(Self {
            family: ladder.family(),
            rb_dimensions: ladder.rb_levels().iter().map(RbLevel::dim).collect(),
            mu_bits: mus.iter().map(|value| value.to_bits()).collect(),
            tolerance_bits: tolerances.iter().map(|value| value.to_bits()).collect(),
            strictest_tolerance_bits: strictest_tolerance.to_bits(),
            total_queries,
            per_mu_work_units,
            declared_work_units,
        })
    }

    /// Ladder family captured when this battery was validated.
    #[must_use]
    pub const fn family(&self) -> LadderFamily {
        self.family
    }

    /// Exact retained basis dimensions, finest rung first.
    #[must_use]
    pub fn rb_dimensions(&self) -> &[usize] {
        &self.rb_dimensions
    }

    /// Parameter-axis values as exact IEEE-754 bit patterns.
    #[must_use]
    pub fn mu_bits(&self) -> &[u64] {
        &self.mu_bits
    }

    /// Tolerance-axis values as exact IEEE-754 bit patterns.
    #[must_use]
    pub fn tolerance_bits(&self) -> &[u64] {
        &self.tolerance_bits
    }

    /// Cartesian query count.
    #[must_use]
    pub const fn total_queries(&self) -> usize {
        self.total_queries
    }

    /// Conservative admitted work for one parameter tile.
    #[must_use]
    pub const fn per_mu_work_units(&self) -> usize {
        self.per_mu_work_units
    }

    /// Conservative admitted work for the complete battery.
    #[must_use]
    pub const fn declared_work_units(&self) -> usize {
        self.declared_work_units
    }

    fn mu(&self, index: usize) -> f64 {
        f64::from_bits(self.mu_bits[index])
    }

    fn tolerance(&self, index: usize) -> f64 {
        f64::from_bits(self.tolerance_bits[index])
    }

    fn strictest_tolerance(&self) -> f64 {
        f64::from_bits(self.strictest_tolerance_bits)
    }

    fn matches(&self, ladder: &Ladder) -> bool {
        self.family == ladder.family()
            && self.rb_dimensions.len() == ladder.rb_levels().len()
            && self
                .rb_dimensions
                .iter()
                .zip(ladder.rb_levels())
                .all(|(dimension, level)| *dimension == level.dim())
    }
}

/// Replayable progress for the first unfinished parameter in the canonical
/// prefix. These are operational facts only; they carry no coverage fraction
/// and cannot be promoted into scientific authority.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RbCoverageParameterProgress {
    parameter_index: usize,
    rungs_completed: usize,
    truth_fallback_started: bool,
    truth_fallback_completed: bool,
    tolerances_classified: usize,
    logical_work_units: u64,
}

impl RbCoverageParameterProgress {
    /// Zero-based parameter tile in the exact plan order.
    #[must_use]
    pub const fn parameter_index(&self) -> usize {
        self.parameter_index
    }

    /// Fully evaluated RB rungs for this parameter.
    #[must_use]
    pub const fn rungs_completed(&self) -> usize {
        self.rungs_completed
    }

    /// Whether the single admitted truth fallback was entered.
    #[must_use]
    pub const fn truth_fallback_started(&self) -> bool {
        self.truth_fallback_started
    }

    /// Whether the single admitted truth fallback completed.
    #[must_use]
    pub const fn truth_fallback_completed(&self) -> bool {
        self.truth_fallback_completed
    }

    /// Tolerances classified locally but not published as a fraction.
    #[must_use]
    pub const fn tolerances_classified(&self) -> usize {
        self.tolerances_classified
    }

    /// Completed logical scalar operations observed by the fixed-stride
    /// checkpoint meter.
    #[must_use]
    pub const fn logical_work_units(&self) -> u64 {
        self.logical_work_units
    }
}

/// Semantic receipt retained after every drained production run. It contains
/// only the longest contiguous completed parameter prefix plus the first
/// unfinished parameter; out-of-order scheduler completions are excluded and
/// remain visible only in [`RunReport`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RbCoverageProgressReceipt {
    plan: RbCoveragePlan,
    run: RunId,
    tile_budget: Budget,
    declared_memory_bytes: usize,
    completed_parameter_prefix: usize,
    rb_queries_in_prefix: usize,
    truth_fallbacks_in_prefix: usize,
    logical_work_units_in_prefix: u64,
    current_parameter: Option<RbCoverageParameterProgress>,
    finalized: bool,
}

impl RbCoverageProgressReceipt {
    /// Exact battery and ladder shape supporting this receipt.
    #[must_use]
    pub const fn plan(&self) -> &RbCoveragePlan {
        &self.plan
    }

    /// Caller-declared logical run identity.
    #[must_use]
    pub const fn run(&self) -> RunId {
        self.run
    }

    /// Per-tile execution budget stamped into every logical parameter tile.
    #[must_use]
    pub const fn tile_budget(&self) -> Budget {
        self.tile_budget
    }

    /// Conservative live scratch reservation held through drain/finalization.
    #[must_use]
    pub const fn declared_memory_bytes(&self) -> usize {
        self.declared_memory_bytes
    }

    /// Longest contiguous fully committed parameter prefix.
    #[must_use]
    pub const fn completed_parameter_prefix(&self) -> usize {
        self.completed_parameter_prefix
    }

    /// RB queries completed inside the retained prefix.
    #[must_use]
    pub const fn rb_queries_in_prefix(&self) -> usize {
        self.rb_queries_in_prefix
    }

    /// Truth fallbacks completed inside the retained prefix.
    #[must_use]
    pub const fn truth_fallbacks_in_prefix(&self) -> usize {
        self.truth_fallbacks_in_prefix
    }

    /// Logical scalar work retained by the completed prefix.
    #[must_use]
    pub const fn logical_work_units_in_prefix(&self) -> u64 {
        self.logical_work_units_in_prefix
    }

    /// Progress for the first unfinished parameter, if any. Later completed
    /// tiles are deliberately not semantic evidence.
    #[must_use]
    pub const fn current_parameter(&self) -> Option<&RbCoverageParameterProgress> {
        self.current_parameter.as_ref()
    }

    /// True only after every parameter drained successfully and the final
    /// ambient cancellation checkpoint admitted publication.
    #[must_use]
    pub const fn is_finalized(&self) -> bool {
        self.finalized
    }
}

/// Scientific outcome of one drained production coverage run. The incomplete
/// variant has no fraction field by construction.
#[derive(Debug, Clone, PartialEq)]
pub enum RbCoverageOutcome {
    /// Every parameter completed and final publication was admitted.
    Complete {
        /// Fraction of the exact Cartesian battery served by an RB rung.
        coverage: f64,
        /// Number of covered parameter/tolerance pairs.
        covered_queries: usize,
        /// Estimated-only authority: RB bounds remain round-to-nearest f64.
        authority: Color,
        /// Finalized replay receipt.
        receipt: RbCoverageProgressReceipt,
    },
    /// Cancellation was observed, all workers drained, and no fraction was
    /// published. The no-claim certificate is absorbing downstream.
    Incomplete {
        /// Explicit refusal to make a numerical coverage claim.
        no_claim: NumericalCertificate,
        /// Deterministic retained prefix/progress receipt.
        receipt: RbCoverageProgressReceipt,
    },
}

impl RbCoverageOutcome {
    /// Semantic progress/finalization receipt for either outcome.
    #[must_use]
    pub const fn receipt(&self) -> &RbCoverageProgressReceipt {
        match self {
            Self::Complete { receipt, .. } | Self::Incomplete { receipt, .. } => receipt,
        }
    }
}

/// Production coverage result plus schedule/latency and memory-accounting
/// diagnostics. Measured fields never enter the semantic progress receipt.
#[derive(Debug, Clone, PartialEq)]
pub struct RbCoverageRun {
    outcome: RbCoverageOutcome,
    report: RunReport,
    memory: LeaseReceipt,
}

impl RbCoverageRun {
    /// Scientific complete/incomplete outcome.
    #[must_use]
    pub const fn outcome(&self) -> &RbCoverageOutcome {
        &self.outcome
    }

    /// Executor diagnostics, including measured cancellation latency.
    #[must_use]
    pub const fn report(&self) -> &RunReport {
        &self.report
    }

    /// Canonical operation-memory accounting snapshot after scratch release.
    #[must_use]
    pub const fn memory_receipt(&self) -> &LeaseReceipt {
        &self.memory
    }
}

#[derive(Debug, Clone)]
struct CoverageParameterResult {
    covered_queries: usize,
    progress: RbCoverageParameterProgress,
}

#[derive(Debug)]
struct CoverageSlot {
    progress: RbCoverageParameterProgress,
    result: Option<CoverageParameterResult>,
    error: Option<SurrogateError>,
}

impl CoverageSlot {
    fn new(parameter_index: usize) -> Self {
        Self {
            progress: RbCoverageParameterProgress {
                parameter_index,
                rungs_completed: 0,
                truth_fallback_started: false,
                truth_fallback_completed: false,
                tolerances_classified: 0,
                logical_work_units: 0,
            },
            result: None,
            error: None,
        }
    }
}

struct CoverageWorkPoll<'a, Caps> {
    task_cx: &'a asupersync::Cx<Caps>,
    tile_cx: &'a fs_exec::Cx<'a>,
    gate: &'a CancelGate,
    slot: &'a Mutex<CoverageSlot>,
    progress: RbCoverageParameterProgress,
    work_since_checkpoint: u64,
}

impl<Caps> CoverageWorkPoll<'_, Caps> {
    fn sync_progress(&self) {
        let mut slot = self
            .slot
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        slot.progress = self.progress.clone();
    }

    fn record_error(&self, error: SurrogateError) {
        let mut slot = self
            .slot
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        slot.progress = self.progress.clone();
        slot.error = Some(error);
    }

    fn commit(&self, covered_queries: usize) {
        let mut slot = self
            .slot
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        slot.progress = self.progress.clone();
        slot.result = Some(CoverageParameterResult {
            covered_queries,
            progress: self.progress.clone(),
        });
    }

    fn observe_cancellation(&self) -> Result<(), Cancelled> {
        self.sync_progress();
        if self.tile_cx.checkpoint().is_err() {
            return Err(Cancelled);
        }
        if self.task_cx.checkpoint().is_err() {
            self.gate.request();
            return Err(Cancelled);
        }
        Ok(())
    }
}

impl<Caps> WorkPoll for CoverageWorkPoll<'_, Caps> {
    type Cancel = Cancelled;

    fn work(&mut self, units: u64) -> Result<(), Self::Cancel> {
        self.progress.logical_work_units =
            self.progress.logical_work_units.saturating_add(units);
        self.work_since_checkpoint = self.work_since_checkpoint.saturating_add(units);
        if self.work_since_checkpoint >= RB_COVERAGE_CHECKPOINT_WORK_UNITS {
            self.work_since_checkpoint %= RB_COVERAGE_CHECKPOINT_WORK_UNITS;
            self.observe_cancellation()?;
        }
        Ok(())
    }

    fn checkpoint(&mut self) -> Result<(), Self::Cancel> {
        self.work_since_checkpoint = 0;
        self.observe_cancellation()
    }
}

enum CoverageComputeError {
    Cancelled(Cancelled),
    Numerical(SurrogateError),
}

impl From<PolledError<Cancelled>> for CoverageComputeError {
    fn from(error: PolledError<Cancelled>) -> Self {
        match error {
            PolledError::Cancelled(cancelled) => Self::Cancelled(cancelled),
            PolledError::Numerical(error) => Self::Numerical(error),
        }
    }
}

struct CoverageKernel<'a, Caps> {
    ladder: &'a Ladder,
    plan: &'a RbCoveragePlan,
    task_cx: &'a asupersync::Cx<Caps>,
    gate: &'a CancelGate,
    slots: &'a [Mutex<CoverageSlot>],
}

impl<Caps: Send + Sync + 'static> TileKernel for CoverageKernel<'_, Caps> {
    type Out = ();

    fn tiles(&self) -> TilePlan {
        TilePlan::new(COVERAGE_KERNEL_NAME, self.plan.mu_bits.len() as u64)
    }

    fn run(&self, tile: u64, tile_cx: &fs_exec::Cx<'_>) -> ControlFlow<Cancelled, Self::Out> {
        let parameter_index = tile as usize;
        let slot = &self.slots[parameter_index];
        let mut poll = CoverageWorkPoll {
            task_cx: self.task_cx,
            tile_cx,
            gate: self.gate,
            slot,
            progress: RbCoverageParameterProgress {
                parameter_index,
                rungs_completed: 0,
                truth_fallback_started: false,
                truth_fallback_completed: false,
                tolerances_classified: 0,
                logical_work_units: 0,
            },
            work_since_checkpoint: 0,
        };

        let evaluated = (|| -> Result<usize, CoverageComputeError> {
            poll.checkpoint().map_err(CoverageComputeError::Cancelled)?;
            let mu = self.plan.mu(parameter_index);
            let mut best_bound = f64::INFINITY;
            for rb in self.ladder.rb_levels.iter().rev() {
                poll.checkpoint().map_err(CoverageComputeError::Cancelled)?;
                let (state, _, _, qoi_bound) = rb.query_polled(mu, &mut poll)?;
                drop(state);
                best_bound = best_bound.min(qoi_bound);
                poll.progress.rungs_completed += 1;
                poll.checkpoint().map_err(CoverageComputeError::Cancelled)?;
                if best_bound <= self.plan.strictest_tolerance() {
                    break;
                }
            }
            if best_bound > self.plan.strictest_tolerance() {
                poll.checkpoint().map_err(CoverageComputeError::Cancelled)?;
                poll.progress.truth_fallback_started = true;
                let truth_state = self.ladder.truth.solve_polled(mu, &mut poll)?;
                let _ = self
                    .ladder
                    .truth
                    .compliance_polled(&truth_state, &mut poll)?;
                poll.progress.truth_fallback_completed = true;
                poll.checkpoint().map_err(CoverageComputeError::Cancelled)?;
            }

            let mut covered_queries = 0usize;
            for tolerance_index in 0..self.plan.tolerance_bits.len() {
                if best_bound <= self.plan.tolerance(tolerance_index) {
                    covered_queries += 1;
                }
                poll.progress.tolerances_classified += 1;
                poll.work(1).map_err(CoverageComputeError::Cancelled)?;
            }
            poll.checkpoint().map_err(CoverageComputeError::Cancelled)?;
            Ok(covered_queries)
        })();

        match evaluated {
            Ok(covered_queries) => {
                poll.commit(covered_queries);
                ControlFlow::Continue(())
            }
            Err(CoverageComputeError::Cancelled(cancelled)) => {
                poll.sync_progress();
                ControlFlow::Break(cancelled)
            }
            Err(CoverageComputeError::Numerical(error)) => {
                poll.record_error(error);
                self.gate.request();
                ControlFlow::Break(Cancelled)
            }
        }
    }
}

fn conservative_coverage_memory_bytes(
    ladder: &Ladder,
    plan: &RbCoveragePlan,
    workers: usize,
) -> Result<usize, SurrogateError> {
    let truth_nodes = ladder.truth.n();
    let largest_dimension = ladder
        .rb_levels
        .iter()
        .map(RbLevel::dim)
        .max()
        .unwrap_or(0);
    let truth_scratch = checked_budget_product(
        "RB coverage live scratch memory",
        &[truth_nodes, 8, size_of::<f64>()],
    )?;
    let dense_values = checked_budget_product(
        "RB coverage live scratch memory",
        &[largest_dimension, largest_dimension, size_of::<f64>()],
    )?;
    let dense_rows = checked_budget_product(
        "RB coverage live scratch memory",
        &[largest_dimension, size_of::<Vec<f64>>()],
    )?;
    let reduced_vectors = checked_budget_product(
        "RB coverage live scratch memory",
        &[largest_dimension, 3, size_of::<f64>()],
    )?;
    let per_tile = truth_scratch
        .checked_add(dense_values)
        .and_then(|bytes| bytes.checked_add(dense_rows))
        .and_then(|bytes| bytes.checked_add(reduced_vectors))
        .ok_or(SurrogateError::BudgetArithmeticOverflow {
            resource: "RB coverage live scratch memory",
        })?;
    let active_tiles = workers.max(1).min(plan.mu_bits.len());
    let tile_scratch = per_tile.checked_mul(active_tiles).ok_or(
        SurrogateError::BudgetArithmeticOverflow {
            resource: "RB coverage live scratch memory",
        },
    )?;
    let slots = checked_budget_product(
        "RB coverage live scratch memory",
        &[plan.mu_bits.len(), size_of::<Mutex<CoverageSlot>>()],
    )?;
    let receipt_payload = plan
        .mu_bits
        .len()
        .checked_add(plan.tolerance_bits.len())
        .and_then(|count| count.checked_add(plan.rb_dimensions.len()))
        .and_then(|count| count.checked_mul(size_of::<u64>()))
        .ok_or(SurrogateError::BudgetArithmeticOverflow {
            resource: "RB coverage live scratch memory",
        })?;
    let declared = tile_scratch
        .checked_add(slots)
        .and_then(|bytes| bytes.checked_add(receipt_payload))
        .ok_or(SurrogateError::BudgetArithmeticOverflow {
            resource: "RB coverage live scratch memory",
        })?;
    if declared > MAX_COVERAGE_MEMORY_BYTES {
        return Err(SurrogateError::BudgetExceeded {
            resource: "RB coverage live scratch memory",
            requested: declared,
            limit: MAX_COVERAGE_MEMORY_BYTES,
        });
    }
    Ok(declared)
}

fn progress_receipt(
    plan: &RbCoveragePlan,
    slots: &[Mutex<CoverageSlot>],
    run: RunId,
    tile_budget: Budget,
    declared_memory_bytes: usize,
    finalized: bool,
) -> (RbCoverageProgressReceipt, usize) {
    let mut completed_parameter_prefix = 0usize;
    let mut covered_queries = 0usize;
    let mut rb_queries_in_prefix = 0usize;
    let mut truth_fallbacks_in_prefix = 0usize;
    let mut logical_work_units_in_prefix = 0u64;
    for slot in slots {
        let slot = slot
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        let Some(result) = &slot.result else {
            break;
        };
        completed_parameter_prefix += 1;
        covered_queries += result.covered_queries;
        rb_queries_in_prefix += result.progress.rungs_completed;
        truth_fallbacks_in_prefix += usize::from(result.progress.truth_fallback_completed);
        logical_work_units_in_prefix = logical_work_units_in_prefix
            .saturating_add(result.progress.logical_work_units);
    }
    let current_parameter = slots.get(completed_parameter_prefix).map(|slot| {
        slot.lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner)
            .progress
            .clone()
    });
    (
        RbCoverageProgressReceipt {
            plan: plan.clone(),
            run,
            tile_budget,
            declared_memory_bytes,
            completed_parameter_prefix,
            rb_queries_in_prefix,
            truth_fallbacks_in_prefix,
            logical_work_units_in_prefix,
            current_parameter,
            finalized,
        },
        covered_queries,
    )
}

/// Execute an exact coverage plan as one logical parameter tile per μ under a
/// live asupersync task and the L0 throughput executor. Cancellation is
/// request → drain → finalize: every launched worker joins before this returns,
/// and an observed request yields [`RbCoverageOutcome::Incomplete`] without a
/// fraction. A worst-case scratch charge is reserved on `memory_lease` before
/// slots or workers are created and released only after finalization.
///
/// # Errors
/// A plan/ladder mismatch, static/lease memory refusal, numerical failure, tile
/// panic/refusal, or worker-launch failure. Ordinary cooperative cancellation
/// is an outcome, not an error.
#[allow(clippy::too_many_arguments)]
pub fn rb_coverage_scoped<Caps: Send + Sync + 'static>(
    task_cx: &asupersync::Cx<Caps>,
    pool: &TilePool,
    gate: &CancelGate,
    run: RunId,
    tile_budget: Budget,
    memory_lease: &OperationMemoryLease,
    ladder: &Ladder,
    plan: &RbCoveragePlan,
) -> Result<RbCoverageRun, SurrogateError> {
    if !plan.matches(ladder) {
        return Err(SurrogateError::FamilyMismatch {
            what: "RB coverage plan",
        });
    }
    let declared_memory_bytes = conservative_coverage_memory_bytes(ladder, plan, pool.workers())?;
    if memory_lease.limit_bytes().is_none() {
        return Err(SurrogateError::CoverageExecutionFailed {
            reason: "production coverage requires a bounded operation memory lease".to_string(),
        });
    }
    let scratch_bytes = u64::try_from(declared_memory_bytes).map_err(|_| {
        SurrogateError::BudgetArithmeticOverflow {
            resource: "RB coverage live scratch memory",
        }
    })?;
    let scratch_charge = memory_lease
        .reserve("fs-surrogate-rb-coverage-scratch", scratch_bytes)
        .map_err(|error| SurrogateError::CoverageExecutionFailed {
            reason: error.to_string(),
        })?;
    let slots: Vec<Mutex<CoverageSlot>> = (0..plan.mu_bits.len())
        .map(|parameter_index| Mutex::new(CoverageSlot::new(parameter_index)))
        .collect();
    let kernel = CoverageKernel {
        ladder,
        plan,
        task_cx,
        gate,
        slots: &slots,
    };
    let (run_result, report) = pool.run_scoped(
        task_cx,
        &kernel,
        gate,
        run,
        tile_budget,
        memory_lease,
    );

    for slot in &slots {
        let slot = slot
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        if let Some(error) = &slot.error {
            return Err(error.clone());
        }
    }
    if let Err(error) = &run_result
        && !matches!(error, RunError::Cancelled { .. })
    {
        return Err(SurrogateError::CoverageExecutionFailed {
            reason: error.to_string(),
        });
    }

    let all_parameters_complete = slots.iter().all(|slot| {
        slot.lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner)
            .result
            .is_some()
    });
    let mut finalized = run_result.is_ok() && all_parameters_complete && !gate.is_requested();
    if finalized && task_cx.checkpoint().is_err() {
        gate.request();
        finalized = false;
    }
    let (receipt, covered_queries) = progress_receipt(
        plan,
        &slots,
        run,
        tile_budget,
        declared_memory_bytes,
        finalized,
    );
    let outcome = if finalized {
        #[allow(clippy::cast_precision_loss)]
        let coverage = covered_queries as f64 / plan.total_queries as f64;
        if !coverage.is_finite() {
            return Err(SurrogateError::NonFiniteDerived {
                what: "RB coverage fraction",
            });
        }
        RbCoverageOutcome::Complete {
            coverage,
            covered_queries,
            authority: estimated_color(COVERAGE_ESTIMATOR_ID, f64::INFINITY)?,
            receipt,
        }
    } else {
        RbCoverageOutcome::Incomplete {
            no_claim: NumericalCertificate::no_claim(),
            receipt,
        }
    };
    drop(scratch_charge);
    let memory = memory_lease.receipt();
    Ok(RbCoverageRun {
        outcome,
        report,
        memory,
    })
}

/// THE KILL MEASUREMENT (Proposal A): the fraction of a query battery
/// answerable at an RB rung (f64 estimator ≤ tol) WITHOUT drilling to full
/// order. Below 0.2 the beachhead is too narrow — park certification.
///
/// # Errors
/// Propagates every query refusal (an invalid battery point refuses the whole
/// measurement rather than silently skewing the fraction). This compatibility
/// wrapper has the smaller [`MAX_SYNCHRONOUS_COVERAGE_WORK_UNITS`] ceiling and
/// no interruption claim. Production-scale batteries use
/// [`rb_coverage_scoped`]. Each parameter's RB rungs are evaluated at most
/// once; their bounds classify every tolerance without repeating a solve.
#[allow(clippy::cast_precision_loss)]
pub fn rb_coverage(ladder: &Ladder, mus: &[f64], tols: &[f64]) -> Result<f64, SurrogateError> {
    let plan = RbCoveragePlan::prepare(
        ladder,
        mus,
        tols,
        MAX_SYNCHRONOUS_COVERAGE_WORK_UNITS,
        "synchronous RB coverage total work",
    )?;
    let mut covered = 0usize;
    for parameter_index in 0..plan.mu_bits.len() {
        let mu = plan.mu(parameter_index);
        let mut best_bound = f64::INFINITY;
        for rb in ladder.rb_levels.iter().rev() {
            let (_, _, _, qoi_bound) = rb.query(mu)?;
            best_bound = best_bound.min(qoi_bound);
            if best_bound <= plan.strictest_tolerance() {
                break;
            }
        }
        if best_bound > plan.strictest_tolerance() {
            let truth_state = ladder.truth.solve(mu)?;
            let _ = ladder.truth.compliance(&truth_state)?;
        }
        for tolerance_index in 0..plan.tolerance_bits.len() {
            if best_bound <= plan.tolerance(tolerance_index) {
                covered += 1;
            }
        }
    }
    let coverage = covered as f64 / plan.total_queries as f64;
    if !coverage.is_finite() {
        return Err(SurrogateError::NonFiniteDerived {
            what: "RB coverage fraction",
        });
    }
    Ok(coverage)
}

#[cfg(test)]
mod tests {
    use super::{
        NoWorkPoll, PolledError, RB_COVERAGE_CHECKPOINT_WORK_UNITS, RbLevel, SurrogateError,
        TruthModel, WorkPoll, solve_dense, solve_dense_polled,
    };
    use core::convert::Infallible;

    #[derive(Default)]
    struct RecordingPoll {
        calls: usize,
        largest_work_step: u64,
    }

    impl WorkPoll for RecordingPoll {
        type Cancel = Infallible;

        fn work(&mut self, units: u64) -> Result<(), Self::Cancel> {
            self.calls += 1;
            self.largest_work_step = self.largest_work_step.max(units);
            Ok(())
        }
    }

    struct CancelAfterWork {
        completed: u64,
        cancel_at: u64,
    }

    impl WorkPoll for CancelAfterWork {
        type Cancel = ();

        fn work(&mut self, units: u64) -> Result<(), Self::Cancel> {
            self.completed = self.completed.saturating_add(units);
            if self.completed >= self.cancel_at {
                Err(())
            } else {
                Ok(())
            }
        }
    }

    #[test]
    fn dense_solver_refuses_bad_shapes_inputs_and_pivots() {
        let mut missing_row = vec![vec![1.0]];
        let mut rhs = vec![1.0, 2.0];
        assert!(matches!(
            solve_dense(&mut missing_row, &mut rhs),
            Err(SurrogateError::InvalidVectorLength { .. })
        ));

        let mut non_finite = vec![vec![f64::NAN]];
        let mut rhs = vec![1.0];
        assert!(matches!(
            solve_dense(&mut non_finite, &mut rhs),
            Err(SurrogateError::NonFiniteInput { .. })
        ));

        let mut singular = vec![vec![1.0, 1.0], vec![2.0, 2.0]];
        let mut rhs = vec![1.0, 2.0];
        assert!(matches!(
            solve_dense(&mut singular, &mut rhs),
            Err(SurrogateError::SingularSystem { column: 1 })
        ));
    }

    #[test]
    fn numerical_coverage_path_never_hides_more_than_one_checkpoint_stride() {
        let truth = TruthModel::new(512).expect("bounded truth");
        let rb = RbLevel::train(&truth, (0.0, 4.0), 4).expect("bounded RB rung");
        let mut poll = RecordingPoll::default();
        let state = truth
            .solve_polled(1.5, &mut poll)
            .expect("recorded truth solve");
        truth
            .energy_polled(&state, &state, 1.5, &mut poll)
            .expect("recorded energy");
        truth
            .compliance_polled(&state, &mut poll)
            .expect("recorded compliance");
        rb.query_polled(1.5, &mut poll)
            .expect("recorded reduced query");
        let mut matrix = vec![vec![4.0, 1.0], vec![1.0, 3.0]];
        let mut rhs = vec![1.0, 2.0];
        solve_dense_polled(&mut matrix, &mut rhs, &mut poll).expect("recorded dense solve");
        assert!(poll.calls > 512, "long loops must reach the shared poll seam");
        assert!(
            poll.largest_work_step <= RB_COVERAGE_CHECKPOINT_WORK_UNITS,
            "one unobserved logical step exceeded the cancellation stride: {}",
            poll.largest_work_step
        );
    }

    #[test]
    fn deterministic_work_injection_interrupts_truth_and_reduced_queries() {
        let truth = TruthModel::new(512).expect("bounded truth");
        let mut truth_cancel = CancelAfterWork {
            completed: 0,
            cancel_at: RB_COVERAGE_CHECKPOINT_WORK_UNITS + 17,
        };
        assert!(matches!(
            truth.solve_polled(2.0, &mut truth_cancel),
            Err(PolledError::Cancelled(()))
        ));

        let rb = RbLevel::train(&truth, (0.0, 4.0), 4).expect("bounded RB rung");
        let mut rb_cancel = CancelAfterWork {
            completed: 0,
            cancel_at: 3 * RB_COVERAGE_CHECKPOINT_WORK_UNITS + 5,
        };
        assert!(matches!(
            rb.query_polled(2.0, &mut rb_cancel),
            Err(PolledError::Cancelled(()))
        ));

        let mut no_poll = NoWorkPoll;
        let unbroken = rb
            .query_polled(2.0, &mut no_poll)
            .expect("sealed no-op polling path remains usable");
        let public = rb.query(2.0).expect("public compatibility query");
        assert_eq!(unbroken.1.to_bits(), public.1.to_bits());
        assert_eq!(unbroken.2.to_bits(), public.2.to_bits());
        assert_eq!(unbroken.3.to_bits(), public.3.to_bits());
    }
}
