//! Tiled, cancellable, resumable execution of exact component-by-component
//! lattice construction (bead 6ys.20, execution tranche over the admission
//! tranche in [`crate::cbc`]).
//!
//! The executor performs byte-identically the same arithmetic in the same
//! logical order as [`crate::qmc::Lattice::cbc`] — points ascending within a
//! candidate, candidates ascending within a prefix, exact lowest-candidate
//! tie resolution — so the chosen generator vector is invariant under tile
//! shape and pause/resume splits by construction, not by averaging. Tiling
//! changes only where cancellation and allowance checks may observe the
//! computation, never the bytes it produces.
//!
//! Work accounting debits the SAME conservative per-unit schedule the
//! admission estimate integrates (limb charges at the admitted widths,
//! scalar charges at the declared per-primitive constants), so the running
//! total is monotone, tile-shape independent, and bounded by the admitted
//! `work_units` for every admitted problem. A run-scoped allowance slices
//! that admitted total across `run` calls: exhaustion finalizes at a tile
//! boundary with a replayable state and a named boundary class.
//!
//! Cancellation is request → drain → finalize: the poll is observed at tile
//! boundaries only, the current tile always completes, and the returned
//! state never contains a half-committed generator component (`prefix()`
//! only ever grows by whole chosen components).
//!
//! NO-CLAIM: this tranche does not yet serialize state for cross-process
//! pause/migrate/fork (the state lives in the executor value), does not
//! produce the per-prefix minimality certificate, and does not parallelize
//! candidate scoring. Those remain later 6ys.20 tranches. `korobov_error_sq`
//! stays a diagnostic f64 owned by [`crate::qmc::Lattice`].

use crate::cbc::{CbcAdmission, CbcProblem};
use crate::qmc::{ExactNat, Lattice, exact_kernel_numerator, gcd, lattice_residue};

/// Version of the executor semantics (tile classes, boundary names, debit
/// schedule binding, and cancellation protocol).
pub const CBC_EXECUTOR_SCHEMA_VERSION: u32 = 1;

// The scalar per-primitive constants of the admission schedule (v3). These
// mirror `crate::cbc` and are asserted against the estimate in the battery:
// an executor claiming an admission receipt must debit the same schedule.
const SCALAR_UNITS_PER_LATTICE_VISIT: u128 = 20;
const SCALAR_UNITS_PER_FACTOR: u128 = 1;
const SCALAR_UNITS_PER_FACTOR_LIMB: u128 = 8;
const SCALAR_UNITS_PER_GCD_STEP: u128 = 3;
const SCALAR_UNITS_PER_CANDIDATE: u128 = 16;
const SCALAR_UNITS_PER_DIMENSION: u128 = 8;

/// Cancellation verdict returned by a poll.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CbcControl {
    /// Keep executing.
    Continue,
    /// Request cancellation: drain the current tile, then finalize.
    Cancel,
}

/// The executor's cancellation source. Layer L1 owns no workspace `Cx`;
/// drivers adapt theirs onto this single-method boundary.
pub trait CbcPoll {
    /// Observed at every tile boundary; never inside a tile.
    fn poll(&mut self) -> CbcControl;
}

impl<F: FnMut() -> CbcControl> CbcPoll for F {
    fn poll(&mut self) -> CbcControl {
        self()
    }
}

/// Tile shape: how many candidates and lattice points may be processed
/// between consecutive poll/allowance observations.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CbcTileShape {
    candidate_block: u32,
    point_block: u32,
}

impl CbcTileShape {
    /// Validate a tile shape (both blocks must be at least one).
    ///
    /// # Errors
    /// [`CbcExecError::InvalidTileShape`] when either block is zero.
    pub const fn new(candidate_block: u32, point_block: u32) -> Result<Self, CbcExecError> {
        if candidate_block == 0 || point_block == 0 {
            return Err(CbcExecError::InvalidTileShape {
                candidate_block,
                point_block,
            });
        }
        Ok(Self {
            candidate_block,
            point_block,
        })
    }

    /// Candidates per tile.
    #[must_use]
    pub const fn candidate_block(self) -> u32 {
        self.candidate_block
    }

    /// Lattice points per tile.
    #[must_use]
    pub const fn point_block(self) -> u32 {
        self.point_block
    }
}

/// The tile-boundary class at which a run stopped.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CbcBoundary {
    /// Before any work of a `run` call (zero allowance).
    Entry,
    /// Between lattice-point blocks inside one accumulation or update pass.
    PointBlock,
    /// Between candidate blocks inside one prefix scan.
    CandidateBlock,
    /// Between prefixes (a whole generator component was just committed).
    Prefix,
}

/// Why a `run` call returned without completing the construction.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CbcRunStatus {
    /// Every generator component is chosen; `into_lattice` succeeds.
    Completed,
    /// The poll requested cancellation; the current tile drained and the
    /// state finalized at the named boundary. Resumable.
    Cancelled(CbcBoundary),
    /// The run-scoped work allowance was exhausted at the named boundary.
    /// Resumable.
    AllowanceExhausted(CbcBoundary),
}

/// Executor refusals. Every variant is fail-closed and leaves the state
/// unchanged (construction) or replayable (runtime).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CbcExecError {
    /// A tile block was zero.
    InvalidTileShape {
        /// Requested candidates per tile.
        candidate_block: u32,
        /// Requested points per tile.
        point_block: u32,
    },
    /// The executor's conservative debits exceeded the admitted work bound —
    /// a schedule-conformance invariant breach, never a normal outcome.
    ScheduleOverrun {
        /// Units debited so far.
        spent: u128,
        /// Units the admission covered.
        admitted: u128,
    },
    /// `run` was called after completion.
    AlreadyComplete,
}

/// One in-flight candidate accumulation (points ascending).
#[derive(Debug, Clone)]
struct ScanAccum {
    score: ExactNat,
    next_point: u32,
}

/// The resumable phase cursor. `z` only ever grows by whole components.
#[derive(Debug, Clone)]
enum Phase {
    /// First-component product initialization (candidate 1, points ascending).
    Init { next_point: u32 },
    /// Scanning candidates for the next component.
    Scan {
        candidate: u32,
        accum: Option<ScanAccum>,
        best: Option<(ExactNat, u32)>,
    },
    /// Folding the chosen candidate into the prefix products.
    Update { chosen: u32, next_point: u32 },
    /// All components chosen.
    Done,
}

/// Tiled exact-CBC executor. See the module docs for the determinism,
/// accounting, and cancellation contracts.
#[derive(Debug)]
pub struct CbcExecutor {
    problem: CbcProblem,
    admitted_work_units: u128,
    // Debit-schedule widths snapshotted from the admission estimate.
    kernel_factor_limbs: u128,
    max_source_product_limbs: u128,
    score_capacity_limbs: u128,
    product_capacity_limbs: u128,
    max_score_limbs: u128,
    gcd_step_upper_bound: u128,
    products: Vec<ExactNat>,
    z: Vec<u32>,
    phase: Phase,
    work_spent: u128,
}

impl CbcExecutor {
    /// Build an executor from an admission receipt. Allocates the product
    /// table at the admitted exact limb capacity so in-envelope arithmetic
    /// never reallocates (the admission contract's exact-capacity storage
    /// requirement).
    #[must_use]
    pub fn new(admission: CbcAdmission) -> Self {
        let problem = admission.problem();
        let estimate = admission.estimate();
        let point_count = usize::try_from(problem.point_count())
            .expect("admission target bounds proved the point count fits usize");
        let product_capacity = usize::try_from(estimate.product_capacity_limbs())
            .expect("admission target bounds proved the product capacity fits usize");
        let mut products = vec![ExactNat::one(); point_count];
        for product in &mut products {
            product.reserve_exact_limbs(product_capacity);
        }
        let gcd_step_upper_bound = u128::from(ceil_log2(problem.point_count())) * 2 + 1;
        Self {
            problem,
            admitted_work_units: estimate.work_units(),
            kernel_factor_limbs: estimate.kernel_factor_limbs(),
            max_source_product_limbs: estimate.max_source_product_limbs(),
            score_capacity_limbs: estimate.score_capacity_limbs(),
            product_capacity_limbs: estimate.product_capacity_limbs(),
            max_score_limbs: estimate.max_score_limbs(),
            gcd_step_upper_bound,
            products,
            z: Vec::with_capacity(problem.dimension()),
            phase: Phase::Init { next_point: 0 },
            work_spent: 0,
        }
    }

    /// The admitted problem.
    #[must_use]
    pub const fn problem(&self) -> CbcProblem {
        self.problem
    }

    /// Whole generator components committed so far (never half-committed).
    #[must_use]
    pub fn prefix(&self) -> &[u32] {
        &self.z
    }

    /// Conservative schedule units debited so far.
    #[must_use]
    pub const fn work_spent(&self) -> u128 {
        self.work_spent
    }

    /// Whether construction is complete.
    #[must_use]
    pub const fn is_complete(&self) -> bool {
        matches!(self.phase, Phase::Done)
    }

    /// Consume the executor; `Some` exactly when complete.
    #[must_use]
    pub fn into_lattice(self) -> Option<Lattice> {
        if matches!(self.phase, Phase::Done) {
            Some(Lattice {
                n: self.problem.point_count(),
                z: self.z,
            })
        } else {
            None
        }
    }

    /// Execute tiles until completion, cancellation, or allowance
    /// exhaustion. `allowance` is a run-scoped slice of the admitted work
    /// budget in the same units; zero performs no work.
    ///
    /// # Errors
    /// [`CbcExecError::AlreadyComplete`] when called after completion, or
    /// [`CbcExecError::ScheduleOverrun`] if debits ever exceed the admitted
    /// bound (an invariant breach, never a normal outcome).
    pub fn run(
        &mut self,
        poll: &mut dyn CbcPoll,
        tile: CbcTileShape,
        allowance: u128,
    ) -> Result<CbcRunStatus, CbcExecError> {
        if matches!(self.phase, Phase::Done) {
            return Err(CbcExecError::AlreadyComplete);
        }
        let mut remaining = allowance;
        if remaining == 0 {
            return Ok(CbcRunStatus::AllowanceExhausted(CbcBoundary::Entry));
        }
        loop {
            let boundary = self.execute_tile(tile, &mut remaining)?;
            if matches!(self.phase, Phase::Done) {
                return Ok(CbcRunStatus::Completed);
            }
            if remaining == 0 {
                return Ok(CbcRunStatus::AllowanceExhausted(boundary));
            }
            if matches!(poll.poll(), CbcControl::Cancel) {
                return Ok(CbcRunStatus::Cancelled(boundary));
            }
        }
    }

    /// Execute exactly one tile (or less at a phase edge) and return the
    /// boundary reached. Debits saturate the run allowance: a tile always
    /// completes once started, so `remaining` reaching zero is observed at
    /// the boundary, never inside the tile.
    fn execute_tile(
        &mut self,
        tile: CbcTileShape,
        remaining: &mut u128,
    ) -> Result<CbcBoundary, CbcExecError> {
        let n = self.problem.point_count();
        let dimension = self.problem.dimension();
        // Charges at admitted widths (the estimate's schedule, distributed).
        let visit_limb_units = self
            .max_source_product_limbs
            .checked_mul(self.kernel_factor_limbs)
            .and_then(|units| {
                self.max_source_product_limbs
                    .checked_mul(self.score_capacity_limbs)
                    .and_then(|carry| units.checked_add(carry))
            })
            .expect("admission proved per-visit limb charges fit u128");
        let visit_scalar_units = SCALAR_UNITS_PER_LATTICE_VISIT
            + SCALAR_UNITS_PER_FACTOR
            + self.kernel_factor_limbs * SCALAR_UNITS_PER_FACTOR_LIMB;
        let visit_units = visit_limb_units + visit_scalar_units;
        let candidate_control_units = SCALAR_UNITS_PER_CANDIDATE
            + self.gcd_step_upper_bound * SCALAR_UNITS_PER_GCD_STEP
            + 2 * self.score_capacity_limbs // score zero-fill + normalization
            + self.max_score_limbs; // comparison
        let update_visit_units = visit_units + 2 * self.product_capacity_limbs;
        let prefix_control_units = SCALAR_UNITS_PER_DIMENSION + 1; // + z push

        match &mut self.phase {
            Phase::Init { next_point } => {
                // Initialization charges one unit per point plus the update
                // visits themselves (the estimate's `+ points + dimension`
                // tail distributes here and at each z push).
                let end = (*next_point).saturating_add(tile.point_block).min(n);
                for point in *next_point..end {
                    let point_index =
                        usize::try_from(point).expect("admission proved point indices fit usize");
                    let residue = lattice_residue(point_index, 1, n);
                    self.products[point_index]
                        .mul_assign_factor(exact_kernel_numerator(n, residue));
                    debit(
                        &mut self.work_spent,
                        self.admitted_work_units,
                        remaining,
                        update_visit_units + 1,
                    )?;
                }
                *next_point = end;
                if end == n {
                    self.z.push(1);
                    debit(
                        &mut self.work_spent,
                        self.admitted_work_units,
                        remaining,
                        prefix_control_units,
                    )?;
                    self.phase = if dimension == 1 {
                        Phase::Done
                    } else {
                        Phase::Scan {
                            candidate: 1,
                            accum: None,
                            best: None,
                        }
                    };
                    Ok(CbcBoundary::Prefix)
                } else {
                    Ok(CbcBoundary::PointBlock)
                }
            }
            Phase::Scan {
                candidate,
                accum,
                best,
            } => {
                let mut candidates_in_tile = 0_u32;
                loop {
                    if *candidate == n {
                        let (_, chosen) = best
                            .take()
                            .expect("candidate 1 is coprime to every admitted n");
                        self.phase = Phase::Update {
                            chosen,
                            next_point: 0,
                        };
                        return Ok(CbcBoundary::CandidateBlock);
                    }
                    if accum.is_none() {
                        if candidates_in_tile == tile.candidate_block {
                            return Ok(CbcBoundary::CandidateBlock);
                        }
                        candidates_in_tile += 1;
                        debit(
                            &mut self.work_spent,
                            self.admitted_work_units,
                            remaining,
                            candidate_control_units,
                        )?;
                        if gcd(*candidate, n) != 1 {
                            *candidate += 1;
                            continue;
                        }
                        let mut score = ExactNat::zero();
                        score.reserve_exact_limbs(
                            usize::try_from(self.score_capacity_limbs)
                                .expect("admission target bounds proved score capacity fits usize"),
                        );
                        *accum = Some(ScanAccum {
                            score,
                            next_point: 0,
                        });
                    }
                    let running = accum.as_mut().expect("accumulator was just installed");
                    let end = running.next_point.saturating_add(tile.point_block).min(n);
                    for point in running.next_point..end {
                        let point_index = usize::try_from(point)
                            .expect("admission proved point indices fit usize");
                        let residue = lattice_residue(point_index, *candidate, n);
                        running.score.add_mul_factor(
                            &self.products[point_index],
                            exact_kernel_numerator(n, residue),
                        );
                        debit(
                            &mut self.work_spent,
                            self.admitted_work_units,
                            remaining,
                            visit_units,
                        )?;
                    }
                    running.next_point = end;
                    if end < n {
                        return Ok(CbcBoundary::PointBlock);
                    }
                    let finished = accum.take().expect("accumulator finished this candidate");
                    let mut score = finished.score;
                    score.normalize();
                    let replace = match &*best {
                        None => true,
                        Some((best_score, best_candidate)) => {
                            match score.magnitude_cmp(best_score) {
                                core::cmp::Ordering::Less => true,
                                core::cmp::Ordering::Equal => *candidate < *best_candidate,
                                core::cmp::Ordering::Greater => false,
                            }
                        }
                    };
                    if replace {
                        *best = Some((score, *candidate));
                    }
                    *candidate += 1;
                    if *remaining == 0 {
                        return Ok(CbcBoundary::CandidateBlock);
                    }
                }
            }
            Phase::Update { chosen, next_point } => {
                let chosen = *chosen;
                let end = (*next_point).saturating_add(tile.point_block).min(n);
                for point in *next_point..end {
                    let point_index =
                        usize::try_from(point).expect("admission proved point indices fit usize");
                    let residue = lattice_residue(point_index, chosen, n);
                    self.products[point_index]
                        .mul_assign_factor(exact_kernel_numerator(n, residue));
                    debit(
                        &mut self.work_spent,
                        self.admitted_work_units,
                        remaining,
                        update_visit_units,
                    )?;
                }
                *next_point = end;
                if end == n {
                    self.z.push(chosen);
                    debit(
                        &mut self.work_spent,
                        self.admitted_work_units,
                        remaining,
                        prefix_control_units,
                    )?;
                    self.phase = if self.z.len() == dimension {
                        Phase::Done
                    } else {
                        Phase::Scan {
                            candidate: 1,
                            accum: None,
                            best: None,
                        }
                    };
                    Ok(CbcBoundary::Prefix)
                } else {
                    Ok(CbcBoundary::PointBlock)
                }
            }
            Phase::Done => Ok(CbcBoundary::Prefix),
        }
    }
}

/// Debit schedule units against the admitted bound and the run allowance
/// (saturating: a started tile always completes). A free function so phase
/// bindings and the accounting fields can be borrowed disjointly.
fn debit(
    work_spent: &mut u128,
    admitted: u128,
    remaining: &mut u128,
    units: u128,
) -> Result<(), CbcExecError> {
    *work_spent = work_spent
        .checked_add(units)
        .ok_or(CbcExecError::ScheduleOverrun {
            spent: u128::MAX,
            admitted,
        })?;
    if *work_spent > admitted {
        return Err(CbcExecError::ScheduleOverrun {
            spent: *work_spent,
            admitted,
        });
    }
    *remaining = remaining.saturating_sub(units);
    Ok(())
}

fn ceil_log2(value: u32) -> u32 {
    debug_assert!(value >= 1);
    32 - value.saturating_sub(1).leading_zeros()
}
