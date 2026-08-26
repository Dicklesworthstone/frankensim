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
//! pause/migrate/fork (the state lives in the executor value) or parallelize
//! candidate scoring. `korobov_error_sq` stays a diagnostic f64 owned by
//! [`crate::qmc::Lattice`].

use crate::cbc::{CbcAdmission, CbcExecutionMode, CbcExecutionSchedule, CbcProblem};
use crate::cbc_cert::{ADMISSIBLE_RULE_UNITS, CbcPrefixCertificate, TIE_RULE_LOWEST_CANDIDATE};
use crate::cbc_limb::LimbCursor;
use crate::qmc::{ExactNat, Lattice, exact_kernel_numerator, gcd, lattice_residue};

/// Version of the executor semantics (tile classes, boundary names, debit
/// schedule binding, and cancellation protocol). v3 adds the limb-block
/// microstep dimension with persistent limb cursors (bead .20.6).
pub const CBC_EXECUTOR_SCHEMA_VERSION: u32 = 3;

/// Default mutation cells per tile for the two-argument constructor.
/// Callers partitioning below this use [`CbcTileShape::with_limbs`].
pub const DEFAULT_LIMB_BLOCK: u32 = 64;

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
    limb_block: u32,
}

impl CbcTileShape {
    /// Validate a tile shape (both outer blocks must be at least one);
    /// the limb-microstep block defaults to [`DEFAULT_LIMB_BLOCK`].
    ///
    /// # Errors
    /// [`CbcExecError::InvalidTileShape`] when either block is zero.
    pub const fn new(candidate_block: u32, point_block: u32) -> Result<Self, CbcExecError> {
        Self::with_limbs(candidate_block, point_block, DEFAULT_LIMB_BLOCK)
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

    /// Validate a tile shape including the limb-microstep block. All three
    /// blocks must be at least one; each bounds the work between
    /// consecutive poll/allowance observations at its own granularity.
    ///
    /// # Errors
    /// [`CbcExecError::InvalidTileShape`] when an outer block is zero or
    /// [`CbcExecError::InvalidLimbBlock`] when `limb_block` is zero.
    pub const fn with_limbs(
        candidate_block: u32,
        point_block: u32,
        limb_block: u32,
    ) -> Result<Self, CbcExecError> {
        if candidate_block == 0 || point_block == 0 {
            return Err(CbcExecError::InvalidTileShape {
                candidate_block,
                point_block,
            });
        }
        if limb_block == 0 {
            return Err(CbcExecError::InvalidLimbBlock { limb_block });
        }
        Ok(Self {
            candidate_block,
            point_block,
            limb_block,
        })
    }

    /// Mutation cells per tile (the microstep budget).
    #[must_use]
    pub const fn limb_block(self) -> u32 {
        self.limb_block
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
    /// Between mutation cells inside one exact operation's microprogram:
    /// the persisted limb cursor resumes exactly where it paused.
    LimbBlock,
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

/// Executor refusals. Admission-authority and storage-ceiling refusals happen
/// before the rejected operation mutates state. `ScheduleOverrun` is an
/// invariant breach whose transaction/retry semantics remain a separate
/// atomic-execution ratchet; callers must not retry it as a replayable pause.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CbcExecError {
    /// The sealed receipt's schema, target layout, schedule, or covered budget
    /// no longer matches the current admission authority.
    AdmissionAuthorityMismatch,
    /// A tile block was zero.
    InvalidTileShape {
        /// Requested candidates per tile.
        candidate_block: u32,
        /// Requested points per tile.
        point_block: u32,
    },
    /// A limb-microstep block was zero.
    InvalidLimbBlock {
        /// Requested mutation cells per tile.
        limb_block: u32,
    },
    /// The executor's conservative debits exceeded the admitted work bound —
    /// a schedule-conformance invariant breach, never a normal outcome.
    ScheduleOverrun {
        /// Units debited so far.
        spent: u128,
        /// Units the admission covered.
        admitted: u128,
    },
    /// Exact arithmetic requested more limbs than the admission-owned storage
    /// schema permits. Refused before the overflowing arithmetic mutates
    /// executor state.
    StorageScheduleOverrun {
        /// Limbs required by the next exact operation.
        required_limbs: usize,
        /// Limbs admitted for this storage class.
        admitted_limbs: usize,
    },
    /// `run` was called after completion.
    AlreadyComplete,
    /// A fallible-storage admission refused an allocation or a capacity
    /// reservation before any mutation. Carries the full diagnostic surface
    /// required by the fallible-storage authority: storage class, phase,
    /// cursor, requested/admitted/observed counts, and ranked remediations.
    Storage(CbcStorageRefusal),
    /// `enable_certificates` was called after work had already been debited
    /// (certificates must cover every scanned component or none).
    CertificatesAfterStart,
    /// Certificate production was requested from a construction-only receipt.
    CertificatesNotAdmitted,
}

/// Which admitted storage class refused a fallible allocation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CbcStorageClass {
    /// The point-indexed product owner array itself.
    ProductOwnerArray,
    /// One product's limb payload (`cursor` is the point index).
    ProductLimbs,
    /// The chosen-prefix register (`z`).
    ZPrefix,
    /// The retained-certificate record array.
    CertificateRecords,
    /// One certificate's prefix word payload.
    CertificatePrefixScratch,
    /// One certificate's winning/runner-up score limbs.
    CertificateScoreScratch,
    /// One certificate's tie-class words.
    CertificateTieScratch,
    /// A candidate scan's running-score accumulator.
    ScoreAccumulator,
    /// The shared product-multiply scratch buffer.
    MultiplyScratch,
}

/// Where in the resumable control flow a storage refusal occurred.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CbcPhaseKind {
    /// Executor construction from an admitted receipt.
    Construction,
    /// First-component product initialization.
    Init,
    /// Candidate scanning.
    Scan,
    /// Chosen-product update folding.
    Update,
    /// Certified prefix-record emission.
    CertificateEmission,
}

/// Ranked remediations for a storage refusal (best first). Static text
/// keeps the error `Copy` so refusal paths never allocate.
pub const RANKED_STORAGE_REMEDIATIONS: [&str; 2] = [
    "re-admit with a larger memory budget covering this class",
    "resume with a smaller problem inside the same receipt authority",
];

/// A typed, allocation-free storage refusal diagnostic.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CbcStorageRefusal {
    /// Which admitted storage class refused.
    pub class: CbcStorageClass,
    /// Control-flow phase at refusal.
    pub phase: CbcPhaseKind,
    /// Point index, candidate, or record slot the refusal occurred at.
    pub cursor: usize,
    /// Elements requested by the refused reservation.
    pub requested: usize,
    /// Elements admitted for this class by the sealed receipt.
    pub admitted: usize,
    /// Allocator-reported capacity evidence observed at refusal.
    pub observed: usize,
}

/// One pending microprogram's logical position. Copy-only: exposing
/// progress never leaks a mutable internal buffer.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CbcMicroCursor {
    /// Human-readable operation class (charge-class identity).
    pub op: &'static str,
    /// Source limb row.
    pub src_pos: usize,
    /// Factor limb column.
    pub factor_pos: usize,
    /// Destination limb cell.
    pub dst_pos: usize,
    /// Live carry register.
    pub carry: u64,
}

/// Scan-phase microprogram progress receipt.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CbcMicroProgress {
    /// The candidate currently being scored, when one is open.
    pub candidate: Option<u32>,
    /// Candidates charged against the current tile's block budget.
    pub scan_tile_candidates: u32,
    /// The partially advanced exact operation, when one is open.
    pub pending: Option<CbcMicroCursor>,
}

/// Runtime observation separating the admission's requested product payload
/// from allocator-reported capacity. The latter is evidence only, never an
/// admitted upper bound.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CbcStorageObservation {
    requested_product_limbs: usize,
    maximum_product_length_limbs: usize,
    minimum_observed_product_capacity_limbs: usize,
    maximum_observed_product_capacity_limbs: usize,
    requested_certificate_records: usize,
    retained_certificate_records: usize,
    observed_certificate_record_capacity: usize,
    retained_certificate_prefix_words: usize,
    observed_certificate_prefix_capacity_words: usize,
    retained_certificate_score_limbs: usize,
    observed_certificate_score_capacity_limbs: usize,
    retained_certificate_tie_words: usize,
    observed_certificate_tie_capacity_words: usize,
}

impl CbcStorageObservation {
    /// Per-product logical capacity sealed by admission.
    #[must_use]
    pub const fn requested_product_limbs(self) -> usize {
        self.requested_product_limbs
    }

    /// Largest logical product length currently retained.
    #[must_use]
    pub const fn maximum_product_length_limbs(self) -> usize {
        self.maximum_product_length_limbs
    }

    /// Smallest allocator-reported capacity across resident products.
    #[must_use]
    pub const fn minimum_observed_product_capacity_limbs(self) -> usize {
        self.minimum_observed_product_capacity_limbs
    }

    /// Largest allocator-reported capacity across resident products.
    #[must_use]
    pub const fn maximum_observed_product_capacity_limbs(self) -> usize {
        self.maximum_observed_product_capacity_limbs
    }

    /// Certificate record slots requested by a certified receipt.
    #[must_use]
    pub const fn requested_certificate_records(self) -> usize {
        self.requested_certificate_records
    }

    /// Certificate records retained so far.
    #[must_use]
    pub const fn retained_certificate_records(self) -> usize {
        self.retained_certificate_records
    }

    /// Allocator-reported capacity of the certificate owner array.
    #[must_use]
    pub const fn observed_certificate_record_capacity(self) -> usize {
        self.observed_certificate_record_capacity
    }

    /// Logical prefix words retained across all certificates.
    #[must_use]
    pub const fn retained_certificate_prefix_words(self) -> usize {
        self.retained_certificate_prefix_words
    }

    /// Allocator-reported prefix capacities across all certificates.
    #[must_use]
    pub const fn observed_certificate_prefix_capacity_words(self) -> usize {
        self.observed_certificate_prefix_capacity_words
    }

    /// Logical winning/runner-up limbs retained across all certificates.
    #[must_use]
    pub const fn retained_certificate_score_limbs(self) -> usize {
        self.retained_certificate_score_limbs
    }

    /// Allocator-reported score capacities across all certificates.
    #[must_use]
    pub const fn observed_certificate_score_capacity_limbs(self) -> usize {
        self.observed_certificate_score_capacity_limbs
    }

    /// Logical tie-class words retained across all certificates.
    #[must_use]
    pub const fn retained_certificate_tie_words(self) -> usize {
        self.retained_certificate_tie_words
    }

    /// Allocator-reported tie-class capacities across all certificates.
    #[must_use]
    pub const fn observed_certificate_tie_capacity_words(self) -> usize {
        self.observed_certificate_tie_capacity_words
    }
}

/// One in-flight candidate accumulation (points ascending).
#[derive(Debug, Clone)]
struct ScanAccum {
    score: ExactNat,
    next_point: u32,
    /// Persisted limb cursor for a partially advanced exact operation.
    /// `None` between operations; `Some` only while a point visit is
    /// mid-microprogram (its visit units were already prepaid).
    micro: Option<LimbCursor>,
    /// The lattice point whose operation `micro` belongs to. Valid
    /// whenever `micro` is `Some`; resuming without this would replay a
    /// cursor against the wrong point (bead .20.6 root cause).
    micro_point: u32,
}

/// The resumable phase cursor. `z` only ever grows by whole components.
#[derive(Debug, Clone)]
enum Phase {
    /// First-component product initialization (candidate 1, points ascending).
    Init { next_point: u32 },
    /// Scanning candidates for the next component.
    Scan {
        candidate: u32,
        /// Candidates charged against the current tile's block budget.
        /// Persistent across run() calls so resumed work consumes the
        /// same observation envelope as fresh work.
        tile_candidates: u32,
        accum: Option<ScanAccum>,
        best: Option<(ExactNat, u32)>,
        runner_up: Option<(ExactNat, u32)>,
        tie_class: Vec<u32>,
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
    schedule: CbcExecutionSchedule,
    score_capacity_limbs: usize,
    product_capacity_limbs: usize,
    admissible_candidates_per_prefix: usize,
    products: Vec<ExactNat>,
    z: Vec<u32>,
    /// Admitted multiply-scratch copy of one product's limbs. Lives inside
    /// the already-charged update-phase overlap envelope; reserving it at
    /// construction keeps every arithmetic transition allocation-free.
    scratch: Vec<u32>,
    phase: Phase,
    work_spent: u128,
    certifying: bool,
    certificates_admitted: bool,
    certificates: Vec<CbcPrefixCertificate>,
}

impl CbcExecutor {
    /// Build an executor from a current admission receipt. Arithmetic lengths
    /// stay inside its requested-capacity schedule; allocator rounding remains
    /// outside the receipt's memory claim.
    ///
    /// # Errors
    /// [`CbcExecError::AdmissionAuthorityMismatch`] if any sealed schema,
    /// schedule, layout, or budget field is stale.
    pub fn new(admission: CbcAdmission) -> Result<Self, CbcExecError> {
        if !admission.has_current_authority() {
            return Err(CbcExecError::AdmissionAuthorityMismatch);
        }
        let problem = admission.problem();
        let estimate = admission.estimate();
        let schedule = admission.execution_schedule();
        let point_count = usize::try_from(problem.point_count())
            .expect("admission target bounds proved the point count fits usize");
        let product_capacity = usize::try_from(estimate.product_capacity_limbs())
            .expect("admission target bounds proved the product capacity fits usize");
        let score_capacity = usize::try_from(estimate.score_capacity_limbs())
            .expect("admission target bounds proved the score capacity fits usize");
        let admissible_candidates_per_prefix =
            usize::try_from(estimate.admissible_candidates_per_prefix())
                .expect("admission target bounds proved the unit-group size fits usize");
        let certificate_capacity = problem.dimension().saturating_sub(1);
        let certificates_admitted = matches!(admission.mode(), CbcExecutionMode::Certified);

        // Fallible-storage authority: every construction allocation is a
        // checked reservation refused with a typed diagnostic before any
        // executor state is published.
        let refuse = |class, cursor, requested, admitted, observed| {
            CbcExecError::Storage(CbcStorageRefusal {
                class,
                phase: CbcPhaseKind::Construction,
                cursor,
                requested,
                admitted,
                observed,
            })
        };

        let mut products = Vec::new();
        products.try_reserve_exact(point_count).map_err(|_| {
            refuse(
                CbcStorageClass::ProductOwnerArray,
                0,
                point_count,
                point_count,
                products.capacity(),
            )
        })?;
        for point_index in 0..point_count {
            let product = ExactNat::try_one_with_capacity(product_capacity).map_err(|_| {
                refuse(
                    CbcStorageClass::ProductLimbs,
                    point_index,
                    product_capacity,
                    product_capacity,
                    0,
                )
            })?;
            products.push(product);
        }

        let mut z = Vec::new();
        z.try_reserve_exact(problem.dimension()).map_err(|_| {
            refuse(
                CbcStorageClass::ZPrefix,
                0,
                problem.dimension(),
                problem.dimension(),
                z.capacity(),
            )
        })?;

        let mut scratch = Vec::new();
        scratch.try_reserve_exact(product_capacity).map_err(|_| {
            refuse(
                CbcStorageClass::MultiplyScratch,
                0,
                product_capacity,
                product_capacity,
                scratch.capacity(),
            )
        })?;

        let mut certificates = Vec::new();
        if certificates_admitted {
            certificates
                .try_reserve_exact(certificate_capacity)
                .map_err(|_| {
                    refuse(
                        CbcStorageClass::CertificateRecords,
                        0,
                        certificate_capacity,
                        certificate_capacity,
                        certificates.capacity(),
                    )
                })?;
        }

        Ok(Self {
            problem,
            admitted_work_units: estimate.work_units(),
            schedule,
            score_capacity_limbs: score_capacity,
            product_capacity_limbs: product_capacity,
            admissible_candidates_per_prefix,
            products,
            z,
            scratch,
            phase: Phase::Init { next_point: 0 },
            work_spent: 0,
            certifying: false,
            certificates_admitted,
            certificates,
        })
    }

    /// Enable per-prefix certificate production for every SCANNED component
    /// (the theorem-fixed first component is the [F] ratchet's business).
    /// The receipt must have been produced for
    /// [`CbcExecutionMode::Certified`], whose schema-v4 envelope covers the
    /// retained records, score/tie storage, and emission debits.
    ///
    /// # Errors
    /// [`CbcExecError::CertificatesNotAdmitted`] for a construction-only
    /// receipt, or [`CbcExecError::CertificatesAfterStart`] once any work was
    /// debited.
    pub fn enable_certificates(&mut self) -> Result<(), CbcExecError> {
        if !self.certificates_admitted {
            return Err(CbcExecError::CertificatesNotAdmitted);
        }
        if self.work_spent != 0 {
            return Err(CbcExecError::CertificatesAfterStart);
        }
        self.certifying = true;
        Ok(())
    }

    /// Certificates emitted so far (one per committed scanned component,
    /// in commit order; empty unless enabled).
    #[must_use]
    pub fn certificates(&self) -> &[CbcPrefixCertificate] {
        &self.certificates
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

    /// Logical microprogram progress: the persisted cursor of any partially
    /// advanced exact operation plus the scan tile's candidate counter.
    /// Exposes positions only — never references to internal buffers.
    #[must_use]
    pub fn micro_progress(&self) -> Option<CbcMicroProgress> {
        let Phase::Scan {
            candidate,
            tile_candidates,
            accum,
            ..
        } = &self.phase
        else {
            return None;
        };
        let candidate = u32::try_from(*candidate).ok();
        let pending = accum
            .as_ref()
            .and_then(|running| running.micro)
            .map(|cursor| CbcMicroCursor {
                op: match cursor.kind {
                    crate::cbc_limb::ExactOpKind::ZeroFill => "zero-fill",
                    crate::cbc_limb::ExactOpKind::AddMultiplyCell => "multiply-add",
                    crate::cbc_limb::ExactOpKind::AddMultiplyDrain => "carry-drain",
                },
                src_pos: cursor.src_pos,
                factor_pos: cursor.factor_pos,
                dst_pos: cursor.dst_pos,
                carry: cursor.carry,
            });
        Some(CbcMicroProgress {
            candidate,
            scan_tile_candidates: *tile_candidates,
            pending,
        })
    }

    /// Observe logical product lengths and allocator-reported capacities.
    /// Only the logical length is constrained by the admission ceiling.
    #[must_use]
    pub fn storage_observation(&self) -> CbcStorageObservation {
        let mut maximum_length = 0;
        let mut minimum_capacity = usize::MAX;
        let mut maximum_capacity = 0;
        for product in &self.products {
            maximum_length = maximum_length.max(product.limbs().len());
            let capacity = product.capacity_limbs();
            minimum_capacity = minimum_capacity.min(capacity);
            maximum_capacity = maximum_capacity.max(capacity);
        }
        let mut retained_prefix_words = 0_usize;
        let mut observed_prefix_capacity_words = 0_usize;
        let mut retained_score_limbs = 0_usize;
        let mut observed_score_capacity_limbs = 0_usize;
        let mut retained_tie_words = 0_usize;
        let mut observed_tie_capacity_words = 0_usize;
        for certificate in &self.certificates {
            retained_prefix_words = retained_prefix_words
                .checked_add(certificate.prefix.len())
                .expect("admission proved retained prefix accounting fits usize");
            observed_prefix_capacity_words = observed_prefix_capacity_words
                .checked_add(certificate.prefix.capacity())
                .expect("observed prefix capacity accounting fits usize");
            retained_score_limbs = retained_score_limbs
                .checked_add(certificate.winning_score_limbs.len())
                .expect("admission proved retained score accounting fits usize");
            observed_score_capacity_limbs = observed_score_capacity_limbs
                .checked_add(certificate.winning_score_limbs.capacity())
                .expect("observed score capacity accounting fits usize");
            if let Some((runner_limbs, _)) = &certificate.runner_up {
                retained_score_limbs = retained_score_limbs
                    .checked_add(runner_limbs.len())
                    .expect("admission proved retained runner accounting fits usize");
                observed_score_capacity_limbs = observed_score_capacity_limbs
                    .checked_add(runner_limbs.capacity())
                    .expect("observed runner capacity accounting fits usize");
            }
            retained_tie_words = retained_tie_words
                .checked_add(certificate.tie_class.len())
                .expect("admission proved retained tie accounting fits usize");
            observed_tie_capacity_words = observed_tie_capacity_words
                .checked_add(certificate.tie_class.capacity())
                .expect("observed tie capacity accounting fits usize");
        }
        CbcStorageObservation {
            requested_product_limbs: self.product_capacity_limbs,
            maximum_product_length_limbs: maximum_length,
            minimum_observed_product_capacity_limbs: minimum_capacity,
            maximum_observed_product_capacity_limbs: maximum_capacity,
            requested_certificate_records: if self.certificates_admitted {
                self.problem.dimension().saturating_sub(1)
            } else {
                0
            },
            retained_certificate_records: self.certificates.len(),
            observed_certificate_record_capacity: self.certificates.capacity(),
            retained_certificate_prefix_words: retained_prefix_words,
            observed_certificate_prefix_capacity_words: observed_prefix_capacity_words,
            retained_certificate_score_limbs: retained_score_limbs,
            observed_certificate_score_capacity_limbs: observed_score_capacity_limbs,
            retained_certificate_tie_words: retained_tie_words,
            observed_certificate_tie_capacity_words: observed_tie_capacity_words,
        }
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
        let charges = TileCharges::sealed(self, tile);
        // Dispatch by phase with each borrow released before the mutable method call.
        if let Phase::Init { next_point } = &self.phase {
            let cursor = *next_point;
            return self.run_init_tile(cursor, charges.point_block, remaining);
        }
        if let Phase::Update { chosen, next_point } = &self.phase {
            let (fold, cursor) = (*chosen, *next_point);
            return self.run_update_tile(fold, cursor, charges.point_block, remaining);
        }
        if matches!(self.phase, Phase::Done) {
            return Ok(CbcBoundary::Prefix);
        }
        let Phase::Scan {
            candidate,
            tile_candidates,
            accum,
            best,
            runner_up,
            tie_class,
        } = &mut self.phase
        else {
            unreachable!("every executor phase was handled above");
        };
        loop {
            if *candidate == charges.n {
                let (winning_score, chosen) = best
                    .take()
                    .expect("candidate 1 is coprime to every admitted n");
                if charges.certifying {
                    let prefix_len = self
                        .z
                        .len()
                        .checked_add(1)
                        .expect("admission proved the certificate prefix length fits usize");
                    let certificate_units = self
                        .schedule
                        .certificate_prefix_units(prefix_len)
                        .expect("admission proved certificate charges fit u128");
                    // Build the full record from borrowed state first: a
                    // storage refusal here leaves every scan field intact,
                    // so the retry re-emits once and never double-charges
                    // the certificate units or loses the tie class.
                    let runner_borrowed = match &*runner_up {
                        Some((score, who)) => Some((score, *who)),
                        None => None,
                    };
                    let certificate = build_prefix_certificate(
                        &self.z,
                        &winning_score,
                        chosen,
                        runner_borrowed,
                        tie_class,
                        prefix_len,
                        charges.n,
                    )?;
                    debit(
                        &mut self.work_spent,
                        charges.admitted_work_units,
                        remaining,
                        certificate_units,
                    )?;
                    debug_assert!(
                        self.certificates.len() < self.certificates.capacity(),
                        "the record array was reserved for one certificate per scanned component"
                    );
                    let _ = runner_up.take();
                    tie_class.clear();
                    self.certificates.push(certificate);
                }
                self.phase = Phase::Update {
                    chosen,
                    next_point: 0,
                };
                return Ok(CbcBoundary::CandidateBlock);
            }
            match advance_scan_candidate(
                candidate,
                accum,
                tile_candidates,
                &self.products,
                &charges,
                &mut self.work_spent,
                remaining,
            )? {
                AdvanceScan::Boundary(boundary) => return Ok(boundary),
                AdvanceScan::Accumulating => continue,
                AdvanceScan::CandidateFinished => {}
            }
            let finished = accum.take().expect("accumulator finished this candidate");
            let mut score = finished.score;
            score.normalize();
            apply_scan_verdict(
                best,
                runner_up,
                tie_class,
                charges.certifying,
                score,
                *candidate,
            );
            *candidate += 1;
            if *remaining == 0 {
                return Ok(CbcBoundary::CandidateBlock);
            }
        }
    }

    /// Initialize first-component products for one tile of points, then
    /// seal the prefix and enter the first scan on completion.
    fn run_init_tile(
        &mut self,
        next_point: u32,
        point_block: u32,
        remaining: &mut u128,
    ) -> Result<CbcBoundary, CbcExecError> {
        let Self {
            products,
            scratch,
            z,
            phase,
            work_spent,
            admitted_work_units,
            schedule,
            product_capacity_limbs,
            admissible_candidates_per_prefix,
            certifying,
            problem,
            ..
        } = self;
        let n = problem.point_count();
        let dimension = problem.dimension();
        // Initialization charges one unit per point plus the update
        // visits themselves (the estimate's `+ points + dimension`
        // tail distributes here and at each z push).
        let end = next_point.saturating_add(point_block).min(n);
        for point in next_point..end {
            let point_index =
                usize::try_from(point).expect("admission proved point indices fit usize");
            let residue = lattice_residue(point_index, 1, n);
            products[point_index]
                .mul_assign_factor_scratch(
                    exact_kernel_numerator(n, residue),
                    *product_capacity_limbs,
                    scratch,
                )
                .map_err(|required_limbs| CbcExecError::StorageScheduleOverrun {
                    required_limbs,
                    admitted_limbs: *product_capacity_limbs,
                })?;
            debit(
                work_spent,
                *admitted_work_units,
                remaining,
                schedule.initialization_visit_units(),
            )?;
        }
        if end == n {
            debug_assert!(
                z.len() < z.capacity(),
                "the prefix register was reserved for the admitted dimension"
            );
            z.push(1);
            debit(
                work_spent,
                *admitted_work_units,
                remaining,
                schedule.prefix_control_units(),
            )?;
            *phase = if dimension == 1 {
                Phase::Done
            } else {
                Phase::Scan {
                    candidate: 1,
                    tile_candidates: 0,
                    accum: None,
                    best: None,
                    runner_up: None,
                    tie_class: try_reserve_tie_class(
                        *certifying,
                        *admissible_candidates_per_prefix,
                        1,
                    )?,
                }
            };
            Ok(CbcBoundary::Prefix)
        } else {
            *phase = Phase::Init { next_point: end };
            Ok(CbcBoundary::PointBlock)
        }
    }

    /// Fold the chosen candidate into the prefix products for one tile of
    /// points, then seal the extended prefix or continue the fold.
    fn run_update_tile(
        &mut self,
        chosen: u32,
        next_point: u32,
        point_block: u32,
        remaining: &mut u128,
    ) -> Result<CbcBoundary, CbcExecError> {
        let Self {
            products,
            scratch,
            z,
            phase,
            work_spent,
            admitted_work_units,
            schedule,
            product_capacity_limbs,
            admissible_candidates_per_prefix,
            certifying,
            problem,
            ..
        } = self;
        let n = problem.point_count();
        let dimension = problem.dimension();
        let end = next_point.saturating_add(point_block).min(n);
        for point in next_point..end {
            let point_index =
                usize::try_from(point).expect("admission proved point indices fit usize");
            let residue = lattice_residue(point_index, chosen, n);
            products[point_index]
                .mul_assign_factor_scratch(
                    exact_kernel_numerator(n, residue),
                    *product_capacity_limbs,
                    scratch,
                )
                .map_err(|required_limbs| CbcExecError::StorageScheduleOverrun {
                    required_limbs,
                    admitted_limbs: *product_capacity_limbs,
                })?;
            debit(
                work_spent,
                *admitted_work_units,
                remaining,
                schedule.product_update_visit_units(),
            )?;
        }
        if end == n {
            debug_assert!(
                z.len() < z.capacity(),
                "the prefix register was reserved for the admitted dimension"
            );
            z.push(chosen);
            debit(
                work_spent,
                *admitted_work_units,
                remaining,
                schedule.prefix_control_units(),
            )?;
            *phase = if z.len() == dimension {
                Phase::Done
            } else {
                Phase::Scan {
                    candidate: 1,
                    tile_candidates: 0,
                    accum: None,
                    best: None,
                    runner_up: None,
                    tie_class: try_reserve_tie_class(
                        *certifying,
                        *admissible_candidates_per_prefix,
                        1,
                    )?,
                }
            };
            Ok(CbcBoundary::Prefix)
        } else {
            *phase = Phase::Update {
                chosen,
                next_point: end,
            };
            Ok(CbcBoundary::PointBlock)
        }
    }
}

/// Fallibly reserve the admitted tie-class scan scratch. A free function so
/// phase-destructured executor fields stay disjoint from the reservation.
fn try_reserve_tie_class(
    certifying: bool,
    admissible_candidates_per_prefix: usize,
    candidate: u32,
) -> Result<Vec<u32>, CbcExecError> {
    let mut tie = Vec::new();
    if certifying {
        tie.try_reserve_exact(admissible_candidates_per_prefix)
            .map_err(|_| {
                CbcExecError::Storage(CbcStorageRefusal {
                    class: CbcStorageClass::CertificateTieScratch,
                    phase: CbcPhaseKind::Scan,
                    cursor: usize::try_from(candidate)
                        .expect("admission proved candidates fit usize"),
                    requested: admissible_candidates_per_prefix,
                    admitted: admissible_candidates_per_prefix,
                    observed: tie.capacity(),
                })
            })?;
    }
    Ok(tie)
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

/// Per-tile charge constants copied once from the sealed admission
/// authority. A snapshot value, so the scan-tile helpers can be free
/// functions borrowing executor fields disjointly.
#[derive(Debug, Clone, Copy)]
struct TileCharges {
    n: u32,
    certifying: bool,
    score_capacity_limbs: usize,
    admitted_work_units: u128,
    visit_units: u128,
    candidate_control_units: u128,
    point_block: u32,
    candidate_block: u32,
    limb_block: u32,
}

impl TileCharges {
    fn sealed(executor: &CbcExecutor, tile: CbcTileShape) -> Self {
        // Every charge comes from the sealed admission authority; this
        // module owns no mirrored limb-width or scalar-debit constants.
        Self {
            n: executor.problem.point_count(),
            certifying: executor.certifying,
            score_capacity_limbs: executor.score_capacity_limbs,
            admitted_work_units: executor.admitted_work_units,
            visit_units: executor.schedule.candidate_visit_units(),
            candidate_control_units: executor
                .schedule
                .candidate_control_units()
                .checked_add(if executor.certifying {
                    executor.schedule.certificate_candidate_units()
                } else {
                    0
                })
                .expect("admission proved candidate charges fit u128"),
            point_block: tile.point_block,
            candidate_block: tile.candidate_block,
            limb_block: tile.limb_block,
        }
    }
}

/// Outcome of comparing a finished candidate score against the committed
/// best. Candidates ascend, so ties never displace the class minimum.
enum Verdict {
    /// Strictly below every committed score: new winner, and a displaced
    /// best becomes the runner-up.
    NewBest,
    /// Equal to the committed winner: joins the tie class.
    Tie,
    /// Strictly above the winner: may displace the runner-up.
    Above,
}

/// Outcome of one scan-loop step after the candidate-completion check.
enum AdvanceScan {
    /// A tile boundary was reached; publish it to the caller.
    Boundary(CbcBoundary),
    /// The candidate advanced (possibly skipped as non-coprime) but its
    /// score is not complete.
    Accumulating,
    /// The candidate's exact score is complete and still needs the
    /// winner/runner-up verdict applied by the caller.
    CandidateFinished,
}

/// Advance the scan loop one step: install the candidate's accumulator
/// when absent (enforcing the per-tile candidate block boundary, debiting
/// the control charge, skipping non-coprime candidates, and reserving the
/// exact score capacity), then accumulate one point block of its running
/// score with a per-visit debit. Returns the finished flag when the whole
/// candidate score is complete.
fn advance_scan_candidate(
    candidate: &mut u32,
    accum: &mut Option<ScanAccum>,
    candidates_in_tile: &mut u32,
    products: &[ExactNat],
    charges: &TileCharges,
    work_spent: &mut u128,
    remaining: &mut u128,
) -> Result<AdvanceScan, CbcExecError> {
    if accum.is_none() {
        if *candidates_in_tile == charges.candidate_block {
            // Tile envelope consumed; the next tile starts a fresh budget.
            *candidates_in_tile = 0;
            return Ok(AdvanceScan::Boundary(CbcBoundary::CandidateBlock));
        }
        *candidates_in_tile += 1;
        let coprime = gcd(*candidate, charges.n) == 1;
        // Fallible reservation precedes the control debit: a storage
        // refusal leaves no mutation behind, so the retry re-executes this
        // candidate identically and never double-charges the schedule.
        let mut score =
            ExactNat::try_zero_with_capacity(charges.score_capacity_limbs).map_err(|_| {
                CbcExecError::Storage(CbcStorageRefusal {
                    class: CbcStorageClass::ScoreAccumulator,
                    phase: CbcPhaseKind::Scan,
                    cursor: usize::try_from(*candidate)
                        .expect("admission proved candidates fit usize"),
                    requested: charges.score_capacity_limbs,
                    admitted: charges.score_capacity_limbs,
                    observed: 0,
                })
            })?;
        if coprime {
            *accum = Some(ScanAccum {
                score,
                next_point: 0,
                micro: None,
                micro_point: 0,
            });
        }
        debit(
            work_spent,
            charges.admitted_work_units,
            remaining,
            charges.candidate_control_units,
        )?;
        if !coprime {
            *candidate += 1;
            return Ok(AdvanceScan::Accumulating);
        }
    }
    let n = charges.n;
    let running = accum.as_mut().expect("accumulator was just installed");
    let end = running
        .next_point
        .saturating_add(charges.point_block)
        .min(n);
    // Resume a paused microprogram at ITS point, not at the block start:
    // the persisted cursor is meaningless against any other operand.
    let mut point = if running.micro.is_some() {
        running.micro_point
    } else {
        running.next_point
    };
    while point < end {
        let point_index = usize::try_from(point).expect("admission proved point indices fit usize");

        // Begin (and prepay) the visit when no microprogram is in flight;
        // resume a persisted cursor without re-debiting.
        if running.micro.is_none() {
            let src = products[point_index].limbs();
            let factor = exact_kernel_numerator(n, residue_of(point_index, *candidate, n));
            let (factor_words, factor_len) = crate::cbc_limb::factor_limbs_u32(factor);
            let required = src
                .len()
                .checked_add(factor_len)
                .and_then(|length| length.checked_add(1))
                .expect("exact CBC required limb count overflow");
            if required
                > charges
                    .score_capacity_limbs
                    .max(running.score.limbs().len())
            {
                return Err(CbcExecError::StorageScheduleOverrun {
                    required_limbs: required,
                    admitted_limbs: charges.score_capacity_limbs,
                });
            }
            debit(
                work_spent,
                charges.admitted_work_units,
                remaining,
                charges.visit_units,
            )?;
            running.micro = Some(crate::cbc_limb::LimbCursor::begin_add_multiply(
                src.len(),
                factor_len,
                running.score.limbs().len(),
                required,
                running.score.capacity_limbs(),
            ));
            running.micro_point = point;
        }

        // Drive up to limb_block mutation cells of the persisted microprogram.
        let src = products[point_index].limbs();
        let factor = exact_kernel_numerator(n, residue_of(point_index, *candidate, n));
        let (factor_words, factor_len) = crate::cbc_limb::factor_limbs_u32(factor);
        let mut cursor = running.micro.expect("cursor installed above");
        let outcome = crate::cbc_limb::step_add_multiply(
            running.score.limbs_mut(),
            src,
            &factor_words,
            factor_len,
            &mut cursor,
            usize::try_from(charges.limb_block).expect("admission proved blocks fit usize"),
        );
        match outcome {
            crate::cbc_limb::StepOutcome::Refused => {
                // Carry propagation hit the admitted capacity ceiling:
                // fail closed with the score-class overrun identity.
                return Err(CbcExecError::StorageScheduleOverrun {
                    required_limbs: cursor.limit + 1,
                    admitted_limbs: charges.score_capacity_limbs,
                });
            }
            crate::cbc_limb::StepOutcome::Advanced => {
                // Budget exhausted mid-operation: persist cursor AND its
                // owning point, then observe a poll.
                running.micro = Some(cursor);
                running.micro_point = point;
                return Ok(AdvanceScan::Boundary(CbcBoundary::LimbBlock));
            }
            crate::cbc_limb::StepOutcome::Complete => {
                running.micro = None;
                point += 1;
            }
        }
    }
    running.next_point = end;
    if end < n {
        return Ok(AdvanceScan::Boundary(CbcBoundary::PointBlock));
    }
    Ok(AdvanceScan::CandidateFinished)
}

/// Residue helper mirroring the historical call form so the micro driver
/// stays readable; pure function of the same sealed inputs.
fn residue_of(point_index: usize, candidate: u32, n: u32) -> u32 {
    lattice_residue(point_index, candidate, n)
}

/// Compare a finished candidate score against the committed best and
/// update the winner/runner-up/tie-class state. Candidates ascend, so a
/// displaced best is the smallest score strictly above the new winner and
/// ties never displace the class minimum.
fn apply_scan_verdict(
    best: &mut Option<(ExactNat, u32)>,
    runner_up: &mut Option<(ExactNat, u32)>,
    tie_class: &mut Vec<u32>,
    certifying: bool,
    score: ExactNat,
    candidate: u32,
) {
    let verdict = match &*best {
        None => Verdict::NewBest,
        Some((best_score, _)) => match score.magnitude_cmp(best_score) {
            core::cmp::Ordering::Less => Verdict::NewBest,
            core::cmp::Ordering::Equal => Verdict::Tie,
            core::cmp::Ordering::Greater => Verdict::Above,
        },
    };
    match verdict {
        Verdict::NewBest => {
            // Candidates ascend, so a displaced best is the
            // smallest score strictly above the new winner.
            let displaced = best.replace((score, candidate));
            if certifying {
                *runner_up = displaced;
                tie_class.clear();
                tie_class.push(candidate);
            }
        }
        Verdict::Tie => {
            // Ascending order keeps the committed winner the
            // class minimum without re-comparison.
            if certifying {
                tie_class.push(candidate);
            }
        }
        Verdict::Above => {
            if certifying {
                let replace_runner = match &*runner_up {
                    None => true,
                    Some((runner_score, _)) => {
                        score.magnitude_cmp(runner_score) == core::cmp::Ordering::Less
                    }
                };
                if replace_runner {
                    *runner_up = Some((score, candidate));
                }
            }
        }
    }
}

/// Assemble the certified prefix record for a completed scan. Pure
/// construction over borrowed scan state: every payload is fallibly cloned
/// into freshly reserved storage, so a refusal leaves the caller's fields
/// untouched and no partial record is ever minted. The caller debits the
/// certificate charge only after this succeeds, then publishes the record.
fn build_prefix_certificate(
    z: &[u32],
    winning_score: &ExactNat,
    chosen: u32,
    runner_up: Option<(&ExactNat, u32)>,
    tie_class: &[u32],
    prefix_len: usize,
    n: u32,
) -> Result<CbcPrefixCertificate, CbcExecError> {
    let refuse = |class: CbcStorageClass, requested: usize| {
        CbcExecError::Storage(CbcStorageRefusal {
            class,
            phase: CbcPhaseKind::CertificateEmission,
            cursor: prefix_len,
            requested,
            admitted: requested,
            observed: 0,
        })
    };
    let mut prefix = Vec::new();
    prefix
        .try_reserve_exact(prefix_len)
        .map_err(|_| refuse(CbcStorageClass::CertificatePrefixScratch, prefix_len))?;
    prefix.extend_from_slice(z);
    prefix.push(chosen);
    let winning_score_limbs = winning_score.try_clone_limbs().map_err(|_| {
        refuse(
            CbcStorageClass::CertificateScoreScratch,
            winning_score.limbs().len(),
        )
    })?;
    let runner_up_limbs = match runner_up {
        Some((score, who)) => {
            let limbs = score.try_clone_limbs().map_err(|_| {
                refuse(
                    CbcStorageClass::CertificateScoreScratch,
                    score.limbs().len(),
                )
            })?;
            Some((limbs, who))
        }
        None => None,
    };
    let mut retained_tie_class = Vec::new();
    retained_tie_class
        .try_reserve_exact(tie_class.len())
        .map_err(|_| refuse(CbcStorageClass::CertificateTieScratch, tie_class.len()))?;
    retained_tie_class.extend_from_slice(tie_class);
    let denominator_exponent =
        u32::try_from(prefix.len()).expect("admitted dimensions fit u32 exponents");
    Ok(CbcPrefixCertificate {
        point_count: n,
        prefix,
        winning_score_limbs,
        tie_class: retained_tie_class,
        runner_up: runner_up_limbs,
        denominator_exponent,
        tie_rule: TIE_RULE_LOWEST_CANDIDATE,
        admissible_rule: ADMISSIBLE_RULE_UNITS,
    })
}

#[cfg(test)]
mod debit_schedule_tests {
    use super::*;
    use crate::cbc::{CbcBudget, CbcExecutionMode};

    fn executor(n: u32, dimension: usize, mode: CbcExecutionMode) -> CbcExecutor {
        let problem = CbcProblem::new(n, dimension).expect("debit fixture is structural");
        let admission = problem
            .admit_for(mode, CbcBudget::UNBOUNDED)
            .expect("debit fixture admits");
        let mut executor = CbcExecutor::new(admission).expect("fresh authority admits executor");
        if mode == CbcExecutionMode::Certified {
            executor
                .enable_certificates()
                .expect("certified fixture enables evidence before work");
        }
        executor
    }

    fn tile(candidate_block: u32, point_block: u32) -> CbcTileShape {
        CbcTileShape::new(candidate_block, point_block).expect("test tile is nonzero")
    }

    fn execute_delta(executor: &mut CbcExecutor, tile: CbcTileShape) -> (CbcBoundary, u128) {
        let before = executor.work_spent;
        let mut remaining = u128::MAX;
        let boundary = executor
            .execute_tile(tile, &mut remaining)
            .expect("isolated admitted debit class executes");
        (
            boundary,
            executor
                .work_spent
                .checked_sub(before)
                .expect("work debit is monotone"),
        )
    }

    fn scan_phase(
        candidate: u32,
        accum: Option<ScanAccum>,
        best: Option<(ExactNat, u32)>,
        tie_class: Vec<u32>,
    ) -> Phase {
        Phase::Scan {
            candidate,
            tile_candidates: 0,
            accum,
            best,
            runner_up: None,
            tie_class,
        }
    }

    /// G0: hard-pinned runtime boundary deltas independently witness every
    /// schema-v4 debit class. These values are intentionally literal rather
    /// than recomputed through `CbcExecutionSchedule` or `CbcEstimate`, so a
    /// compensating producer/consumer formula drift cannot self-certify by
    /// preserving only the aggregate total.
    #[test]
    fn g0_runtime_boundary_deltas_pin_every_debit_class() {
        // n=3,d=1: three initialization visits at 40 units apiece, with the
        // final point followed by the one 9-unit prefix commit: 40+40+49=129.
        let mut dimension_one = executor(3, 1, CbcExecutionMode::Construction);
        assert_eq!(
            execute_delta(&mut dimension_one, tile(1, 1)),
            (CbcBoundary::PointBlock, 40)
        );
        assert_eq!(
            execute_delta(&mut dimension_one, tile(1, 1)),
            (CbcBoundary::PointBlock, 40)
        );
        assert_eq!(
            execute_delta(&mut dimension_one, tile(1, 1)),
            (CbcBoundary::Prefix, 49),
            "the terminal tile is one 40-unit initialization plus 9-unit prefix control"
        );
        assert!(dimension_one.is_complete());
        assert_eq!(dimension_one.work_spent(), 129);

        // n=8,d=3 construction-mode class witnesses. Installing an already
        // admitted in-flight phase isolates the exact production debit site;
        // no schedule accessor participates in the expected constants.
        let mut construction = executor(8, 3, CbcExecutionMode::Construction);
        let mut score = ExactNat::zero();
        score.reserve_exact_limbs(construction.score_capacity_limbs);
        construction.phase = scan_phase(
            1,
            Some(ScanAccum {
                score,
                next_point: 0,
                micro: None,
                micro_point: 0,
            }),
            None,
            Vec::new(),
        );
        assert_eq!(
            execute_delta(&mut construction, tile(1, 1)),
            (CbcBoundary::PointBlock, 33),
            "one candidate/point exact visit"
        );

        construction.phase = scan_phase(2, None, None, Vec::new());
        assert_eq!(
            execute_delta(&mut construction, tile(1, 1)),
            (CbcBoundary::CandidateBlock, 44),
            "non-coprime candidate isolates base candidate control"
        );

        construction.phase = Phase::Update {
            chosen: 1,
            next_point: 0,
        };
        assert_eq!(
            execute_delta(&mut construction, tile(1, 1)),
            (CbcBoundary::PointBlock, 39),
            "one retained-product update"
        );
        construction.phase = Phase::Update {
            chosen: 1,
            next_point: 7,
        };
        assert_eq!(
            execute_delta(&mut construction, tile(1, 1)),
            (CbcBoundary::Prefix, 48),
            "the final update point is 39 units plus 9-unit prefix control"
        );

        // Certified candidate control adds exactly one certificate unit.
        let mut certified = executor(8, 3, CbcExecutionMode::Certified);
        certified.phase = scan_phase(2, None, None, Vec::new());
        assert_eq!(
            execute_delta(&mut certified, tile(1, 1)),
            (CbcBoundary::CandidateBlock, 45),
            "certified non-coprime control is base 44 plus one evidence unit"
        );

        // Certificate publication is prefix-length sensitive: the first and
        // second scanned components retain 2 and 3 prefix words, producing
        // 14- and 15-unit emissions independently of candidate/update work.
        let mut first_certificate = executor(8, 3, CbcExecutionMode::Certified);
        first_certificate.z.push(1);
        first_certificate.phase = scan_phase(8, None, Some((ExactNat::one(), 1)), vec![1]);
        assert_eq!(
            execute_delta(&mut first_certificate, tile(1, 1)),
            (CbcBoundary::CandidateBlock, 14)
        );
        assert_eq!(first_certificate.certificates.len(), 1);

        let mut second_certificate = executor(8, 3, CbcExecutionMode::Certified);
        second_certificate.z.extend_from_slice(&[1, 1]);
        second_certificate.phase = scan_phase(8, None, Some((ExactNat::one(), 1)), vec![1]);
        assert_eq!(
            execute_delta(&mut second_certificate, tile(1, 1)),
            (CbcBoundary::CandidateBlock, 15)
        );
        assert_eq!(second_certificate.certificates.len(), 1);
    }
}
