//! Bounded, tile-backed population evaluation for generation-oriented engines.
//!
//! The evaluator is deliberately separate from a particular optimizer: NSGA,
//! CMA, and future population methods share the same semantic identities,
//! admission checks, cancellation boundary, and all-or-nothing publication.
//! A callback must be deterministic for replay to be a scientific claim.
//! The output limit bounds retained semantic result bytes. The executor lease
//! bounds its root metadata, tile arenas, and retained result payload.
//! Arbitrary allocation performed inside a user callback remains outside this
//! small adapter's claim.

use core::ops::ControlFlow;
use fs_alloc::{LeasedVec, OperationMemoryLease};
use fs_exec::{
    Budget, CancelGate, Concat, Cx, RunError, RunId, TileFailure, TileKernel, TilePlan, TilePool,
    TilePoolCompletionDisposition, TilePoolCompletionWitness, TilePoolCompletionWitnessError,
};
use std::collections::BTreeSet;
use std::sync::{Arc, Mutex};

const KERNEL: &str = "fs-dfo/population-evaluation-v1";
const GENERATION_IDENTITY_DOMAIN: &str = "fs-dfo/population-generation-v2";

fn generation_identity_root(
    provenance: PopulationProvenance,
    completion: &TilePoolCompletionWitness,
    evaluation_count: usize,
) -> [u8; 32] {
    let mut preimage = [0_u8; 84];
    preimage[0..4].copy_from_slice(&provenance.schema_version.to_le_bytes());
    preimage[4..12].copy_from_slice(&provenance.generation.to_le_bytes());
    preimage[12..20].copy_from_slice(&provenance.run.0.to_le_bytes());
    preimage[20..28].copy_from_slice(
        &u64::try_from(provenance.individuals)
            .unwrap_or(u64::MAX)
            .to_le_bytes(),
    );
    preimage[28..36].copy_from_slice(
        &u64::try_from(provenance.objective_dimension)
            .unwrap_or(u64::MAX)
            .to_le_bytes(),
    );
    preimage[36..44].copy_from_slice(&provenance.tiles.to_le_bytes());
    preimage[44..52].copy_from_slice(
        &u64::try_from(evaluation_count)
            .unwrap_or(u64::MAX)
            .to_le_bytes(),
    );
    preimage[52..].copy_from_slice(&completion.root_bytes());
    *fs_blake3::hash_domain(GENERATION_IDENTITY_DOMAIN, &preimage).as_bytes()
}

/// One immutable candidate with a caller-stable semantic identity.
#[derive(Debug, Clone, PartialEq)]
pub struct PopulationCandidate {
    /// Identity stable across worker counts and retries.
    pub identity: u64,
    /// Decision coordinates supplied to the objective callback.
    pub decision: Vec<f64>,
}

/// One accepted objective result, kept in semantic candidate order.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PopulationEvaluation<'a> {
    identity: u64,
    objectives: &'a [f64],
    stream_key: u128,
}

impl PopulationEvaluation<'_> {
    /// The input candidate identity.
    #[must_use]
    pub const fn identity(&self) -> u64 {
        self.identity
    }

    /// Objective values in the plan's declared objective order.
    #[must_use]
    pub const fn objectives(&self) -> &[f64] {
        self.objectives
    }

    /// The actual TilePool stream key that produced this result.
    #[must_use]
    pub const fn stream_key(&self) -> u128 {
        self.stream_key
    }
}

/// Checked envelope for one population evaluation request.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PopulationLimits {
    /// Maximum population cardinality admitted by the caller.
    pub max_individuals: usize,
    /// Maximum logical evaluations admitted by the caller.
    pub max_work_units: u64,
    /// Maximum retained semantic result bytes (identity plus objective scalars).
    pub max_output_bytes: u64,
}

/// Replay-relevant request provenance retained with a published generation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PopulationProvenance {
    /// Schema for this semantic envelope.
    pub schema_version: u32,
    /// Logical generation being evaluated.
    pub generation: u64,
    /// Caller-ledgered executor run identity.
    pub run: RunId,
    /// Number of candidates in the complete generation.
    pub individuals: usize,
    /// Objective dimension required of every accepted callback result.
    pub objective_dimension: usize,
    /// Number of contiguous logical tiles in the request.
    pub tiles: u64,
}

/// A complete, atomically publishable generation.
#[derive(Debug, Clone, PartialEq)]
pub struct PopulationGeneration {
    inner: Arc<PopulationGenerationInner>,
}

#[derive(Debug, PartialEq)]
struct PopulationGenerationInner {
    generation: u64,
    // Each key is copied directly from `Cx::stream_key()` in the TilePool.
    // Callers cannot construct this lease-backed state or substitute a seed.
    evaluations: Concat<PopulationRow>,
    provenance: PopulationProvenance,
    completion: TilePoolCompletionWitness,
    identity_root: [u8; 32],
}

impl PopulationGeneration {
    #[must_use]
    pub fn generation(&self) -> u64 {
        self.inner.generation
    }

    /// Accepted evaluations in semantic candidate order.
    #[must_use]
    pub fn evaluations(&self) -> impl ExactSizeIterator<Item = PopulationEvaluation<'_>> + '_ {
        self.inner.evaluations.0.as_slice().iter().map(
            |(_position, identity, (objectives, _actual_dimension, _non_finite), stream_key)| {
                PopulationEvaluation {
                    identity: *identity,
                    objectives: objectives.as_slice(),
                    stream_key: *stream_key,
                }
            },
        )
    }

    #[must_use]
    pub fn evaluation_count(&self) -> usize {
        self.inner.evaluations.0.len()
    }

    #[must_use]
    pub fn provenance(&self) -> PopulationProvenance {
        self.inner.provenance
    }

    /// Executor-minted completion evidence for this exact complete generation.
    #[must_use]
    pub fn completion_witness(&self) -> &TilePoolCompletionWitness {
        &self.inner.completion
    }

    /// Immutable identity binding provenance, retained rows, and completion.
    #[must_use]
    pub fn identity_root(&self) -> [u8; 32] {
        self.inner.identity_root
    }

    fn verify(&self) -> Result<(), PopulationPublishError> {
        self.inner
            .completion
            .verify()
            .map_err(PopulationPublishError::CompletionWitness)?;
        if self.inner.completion.disposition() != TilePoolCompletionDisposition::Completed
            || self.inner.completion.cancellation_requested()
            || self.inner.completion.declared_run() != self.inner.provenance.run
            || self.inner.completion.planned_tiles() != self.inner.provenance.tiles
            || self.inner.completion.completed_tiles() != self.inner.provenance.tiles
            || generation_identity_root(
                self.inner.provenance,
                &self.inner.completion,
                self.evaluation_count(),
            ) != self.inner.identity_root
        {
            return Err(PopulationPublishError::CompletionMismatch);
        }
        Ok(())
    }
}

/// Resumable publisher state. It contains only previously committed work.
#[derive(Debug, Clone, PartialEq)]
pub struct PopulationCheckpoint {
    committed: Option<PopulationGeneration>,
    committed_work_units: u64,
    committed_identity_root: Option<[u8; 32]>,
}

impl PopulationCheckpoint {
    #[must_use]
    pub fn committed(&self) -> Option<&PopulationGeneration> {
        self.committed.as_ref()
    }
    #[must_use]
    pub const fn committed_work_units(&self) -> u64 {
        self.committed_work_units
    }

    /// Identity of the exact immutable generation retained by this checkpoint.
    #[must_use]
    pub const fn committed_identity_root(&self) -> Option<[u8; 32]> {
        self.committed_identity_root
    }

    fn verify(&self) -> Result<(), PopulationPublishError> {
        match (&self.committed, self.committed_identity_root) {
            (None, None) if self.committed_work_units == 0 => Ok(()),
            (Some(generation), Some(identity_root))
                if generation.identity_root() == identity_root
                    && self.committed_work_units
                        >= u64::try_from(generation.evaluation_count()).unwrap_or(u64::MAX) =>
            {
                generation.verify()
            }
            _ => Err(PopulationPublishError::CheckpointMismatch),
        }
    }
}

/// Input or envelope refusal before any tile is launched.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum PopulationRefusal {
    /// Tile width is zero, so no bounded scheduling plan exists.
    ZeroTileWidth,
    /// At least one objective component is required for a population result.
    ZeroObjectiveDimension,
    /// The supplied population exceeds the caller's cardinality cap.
    PopulationLimit { requested: usize, maximum: usize },
    /// The required complete evaluation count exceeds the work cap.
    WorkLimit { requested: u64, maximum: u64 },
    /// The explicit executor cost quota cannot admit the requested work.
    BudgetWorkLimit { requested: u64, maximum: u64 },
    /// Result storage cannot be represented or exceeds the output cap.
    OutputLimit { requested: u64, maximum: u64 },
    /// Candidate identities must be unique within a generation.
    DuplicateIdentity { identity: u64 },
    /// Candidate coordinates are not finite.
    NonFiniteDecision { identity: u64 },
    /// The retained provenance does not describe this exact request plan.
    ProvenanceMismatch,
}

/// Evaluation or publication failure. A failed evaluation never publishes.
#[derive(Debug)]
pub enum PopulationPublishError {
    /// The request was refused before executor entry.
    Refused(PopulationRefusal),
    /// A callback returned an objective vector of a different dimension.
    ObjectiveDimension {
        identity: u64,
        expected: usize,
        actual: usize,
    },
    /// A callback returned a non-finite objective.
    NonFiniteObjective { identity: u64 },
    /// The fs-exec run drained or failed; its structured evidence is preserved.
    Executor(RunError),
    /// Another complete generation was committed while this request ran.
    GenerationConflict { expected: u64, actual: Option<u64> },
    /// A resume would exceed the publisher's caller-declared cumulative work cap.
    CumulativeWorkLimit { requested: u64, maximum: u64 },
    /// Executor completion evidence failed its self-verifier.
    CompletionWitness(TilePoolCompletionWitnessError),
    /// Completion evidence does not describe this complete population result.
    CompletionMismatch,
    /// A checkpoint no longer binds its retained generation identity.
    CheckpointMismatch,
}

/// Runs an admitted population as contiguous fs-exec tiles.
#[derive(Clone, Copy)]
pub struct PopulationEvaluator<'pool> {
    pool: &'pool TilePool,
    tile_width: usize,
}

impl<'pool> PopulationEvaluator<'pool> {
    /// Bind an evaluator to a production tile pool and bounded tile width.
    #[must_use]
    pub const fn new(pool: &'pool TilePool, tile_width: usize) -> Self {
        Self { pool, tile_width }
    }

    /// Evaluate a complete population. Cancellation is request-drain-finalize
    /// through `fs-exec`; no partial result is returned from this method.
    ///
    /// The caller supplies the exact executor budget and operation lease. This
    /// adapter preflights declared work against an explicit cost quota and
    /// forwards the full budget to every tile. The synchronous TilePool seam
    /// has no ambient clock, so deadline and poll-quota enforcement remain the
    /// responsibility of an enclosing scoped executor; callback-owned heap is
    /// likewise outside this adapter's memory claim.
    pub fn evaluate<F>(
        &self,
        candidates: &[PopulationCandidate],
        objective_dimension: usize,
        limits: PopulationLimits,
        provenance: PopulationProvenance,
        gate: &CancelGate,
        budget: Budget,
        lease: &OperationMemoryLease,
        objective: F,
    ) -> Result<PopulationGeneration, PopulationPublishError>
    where
        F: Fn(&[f64]) -> Vec<f64> + Sync,
    {
        self.preflight(candidates, objective_dimension, limits, provenance, budget)?;
        let kernel = PopulationKernel {
            candidates,
            tile_width: self.tile_width,
            objective_dimension,
            objective: &objective,
        };
        let witnessed = self
            .pool
            .run_declared_leased_budgeted_witnessed(&kernel, gate, provenance.run, budget, lease)
            .map_err(PopulationPublishError::CompletionWitness)?;
        witnessed
            .verify_bundle()
            .map_err(PopulationPublishError::CompletionWitness)?;
        let (result, _report, completion) = witnessed.into_parts();
        let evaluations = result.map_err(PopulationPublishError::Executor)?;
        if completion.disposition() != TilePoolCompletionDisposition::Completed
            || completion.cancellation_requested()
            || completion.declared_run() != provenance.run
            || completion.planned_tiles() != provenance.tiles
            || completion.completed_tiles() != provenance.tiles
        {
            return Err(PopulationPublishError::CompletionMismatch);
        }
        for (expected_position, row) in evaluations.0.as_slice().iter().enumerate() {
            let (position, identity, (objectives, actual_dimension, non_finite), _stream_key) = row;
            if *position != expected_position || *identity != candidates[expected_position].identity
            {
                return Err(PopulationPublishError::Refused(
                    PopulationRefusal::ProvenanceMismatch,
                ));
            }
            if *actual_dimension != objective_dimension {
                return Err(PopulationPublishError::ObjectiveDimension {
                    identity: *identity,
                    expected: objective_dimension,
                    actual: *actual_dimension,
                });
            }
            if *non_finite || !objectives.as_slice().iter().all(|value| value.is_finite()) {
                return Err(PopulationPublishError::NonFiniteObjective {
                    identity: *identity,
                });
            }
        }
        Ok(PopulationGeneration {
            inner: Arc::new(PopulationGenerationInner {
                generation: provenance.generation,
                evaluations,
                provenance,
                identity_root: generation_identity_root(provenance, &completion, candidates.len()),
                completion,
            }),
        })
    }

    fn preflight(
        &self,
        candidates: &[PopulationCandidate],
        objective_dimension: usize,
        limits: PopulationLimits,
        provenance: PopulationProvenance,
        budget: Budget,
    ) -> Result<(), PopulationPublishError> {
        if self.tile_width == 0 {
            return Err(PopulationPublishError::Refused(
                PopulationRefusal::ZeroTileWidth,
            ));
        }
        if objective_dimension == 0 {
            return Err(PopulationPublishError::Refused(
                PopulationRefusal::ZeroObjectiveDimension,
            ));
        }
        if candidates.len() > limits.max_individuals {
            return Err(PopulationPublishError::Refused(
                PopulationRefusal::PopulationLimit {
                    requested: candidates.len(),
                    maximum: limits.max_individuals,
                },
            ));
        }
        let work = u64::try_from(candidates.len()).map_err(|_| {
            PopulationPublishError::Refused(PopulationRefusal::WorkLimit {
                requested: u64::MAX,
                maximum: limits.max_work_units,
            })
        })?;
        if work > limits.max_work_units {
            return Err(PopulationPublishError::Refused(
                PopulationRefusal::WorkLimit {
                    requested: work,
                    maximum: limits.max_work_units,
                },
            ));
        }
        if let Some(maximum) = budget.cost_quota
            && work > maximum
        {
            return Err(PopulationPublishError::Refused(
                PopulationRefusal::BudgetWorkLimit {
                    requested: work,
                    maximum,
                },
            ));
        }
        let per_result = u64::try_from(objective_dimension)
            .ok()
            .and_then(|n| n.checked_mul(8))
            .and_then(|n| n.checked_add(8))
            .ok_or(PopulationPublishError::Refused(
                PopulationRefusal::OutputLimit {
                    requested: u64::MAX,
                    maximum: limits.max_output_bytes,
                },
            ))?;
        let output = work
            .checked_mul(per_result)
            .ok_or(PopulationPublishError::Refused(
                PopulationRefusal::OutputLimit {
                    requested: u64::MAX,
                    maximum: limits.max_output_bytes,
                },
            ))?;
        if output > limits.max_output_bytes {
            return Err(PopulationPublishError::Refused(
                PopulationRefusal::OutputLimit {
                    requested: output,
                    maximum: limits.max_output_bytes,
                },
            ));
        }
        let expected_tiles = if candidates.is_empty() {
            0
        } else {
            candidates.len().div_ceil(self.tile_width) as u64
        };
        if provenance.individuals != candidates.len()
            || provenance.objective_dimension != objective_dimension
            || provenance.tiles != expected_tiles
        {
            return Err(PopulationPublishError::Refused(
                PopulationRefusal::ProvenanceMismatch,
            ));
        }
        let mut identities = BTreeSet::new();
        for candidate in candidates {
            if !identities.insert(candidate.identity) {
                return Err(PopulationPublishError::Refused(
                    PopulationRefusal::DuplicateIdentity {
                        identity: candidate.identity,
                    },
                ));
            }
            if !candidate.decision.iter().all(|value| value.is_finite()) {
                return Err(PopulationPublishError::Refused(
                    PopulationRefusal::NonFiniteDecision {
                        identity: candidate.identity,
                    },
                ));
            }
        }
        Ok(())
    }
}

/// Holds one all-or-nothing generation slot for pause/resume/fork callers.
#[derive(Debug)]
pub struct PopulationPublisher {
    state: Mutex<PopulationCheckpoint>,
    max_committed_work_units: u64,
}

impl PopulationPublisher {
    /// Start with no committed generation and an explicit cumulative work cap.
    #[must_use]
    pub fn new(max_committed_work_units: u64) -> Self {
        Self {
            state: Mutex::new(PopulationCheckpoint {
                committed: None,
                committed_work_units: 0,
                committed_identity_root: None,
            }),
            max_committed_work_units,
        }
    }

    /// Resume a prior immutable checkpoint under an explicit cumulative cap.
    pub fn from_checkpoint(
        checkpoint: PopulationCheckpoint,
        max_committed_work_units: u64,
    ) -> Result<Self, PopulationPublishError> {
        checkpoint.verify()?;
        if checkpoint.committed_work_units > max_committed_work_units {
            return Err(PopulationPublishError::CumulativeWorkLimit {
                requested: checkpoint.committed_work_units,
                maximum: max_committed_work_units,
            });
        }
        Ok(Self {
            state: Mutex::new(checkpoint),
            max_committed_work_units,
        })
    }

    /// Copy the complete committed state; cancelled attempts do not appear.
    #[must_use]
    pub fn checkpoint(&self) -> PopulationCheckpoint {
        match self.state.lock() {
            Ok(state) => state.clone(),
            Err(poisoned) => poisoned.into_inner().clone(),
        }
    }

    /// Atomically install a complete next generation, rejecting stale forks.
    ///
    /// Only an executor-complete generation can enter the publication lock.
    /// Cancellation authority is consumed by the witnessed evaluator; callers
    /// cannot race a fresh, mutable gate into this immutable state transition.
    pub fn publish(&self, generation: PopulationGeneration) -> Result<(), PopulationPublishError> {
        generation.verify()?;
        let mut state = match self.state.lock() {
            Ok(state) => state,
            Err(poisoned) => poisoned.into_inner(),
        };
        let actual = state
            .committed
            .as_ref()
            .map(PopulationGeneration::generation);
        let expected = actual.map_or(0, |previous| previous.saturating_add(1));
        if generation.generation() != expected {
            return Err(PopulationPublishError::GenerationConflict { expected, actual });
        }
        let generation_work = u64::try_from(generation.evaluation_count())
            .map_err(|_| PopulationPublishError::GenerationConflict { expected, actual })?;
        let committed_work_units = state
            .committed_work_units
            .checked_add(generation_work)
            .ok_or(PopulationPublishError::GenerationConflict { expected, actual })?;
        if committed_work_units > self.max_committed_work_units {
            return Err(PopulationPublishError::CumulativeWorkLimit {
                requested: committed_work_units,
                maximum: self.max_committed_work_units,
            });
        }
        state.committed_work_units = committed_work_units;
        state.committed_identity_root = Some(generation.identity_root());
        state.committed = Some(generation);
        Ok(())
    }
}

// Position, identity, lease-backed objectives with their callback outcome,
// and the actual TilePool stream key. Every owned payload is admission-visible.
type PopulationRow = (usize, u64, (LeasedVec<f64>, usize, bool), u128);

struct PopulationKernel<'a, F> {
    candidates: &'a [PopulationCandidate],
    tile_width: usize,
    objective_dimension: usize,
    objective: &'a F,
}

impl<F> TileKernel for PopulationKernel<'_, F>
where
    F: Fn(&[f64]) -> Vec<f64> + Sync,
{
    type Out = Concat<PopulationRow>;

    fn tiles(&self) -> TilePlan {
        let tiles = if self.candidates.is_empty() {
            0
        } else {
            self.candidates.len().div_ceil(self.tile_width) as u64
        };
        TilePlan::new(KERNEL, tiles)
    }

    fn run(&self, tile: u64, cx: &Cx<'_>) -> ControlFlow<fs_exec::Cancelled, Self::Out> {
        if cx.checkpoint().is_err() {
            return ControlFlow::Break(fs_exec::Cancelled);
        }
        let start = usize::try_from(tile)
            .ok()
            .and_then(|tile| tile.checked_mul(self.tile_width))
            .unwrap_or(self.candidates.len());
        let end = start
            .saturating_add(self.tile_width)
            .min(self.candidates.len());
        let Some(lease) = cx.lease() else {
            return ControlFlow::Break(fs_exec::Cancelled);
        };
        let mut rows = match LeasedVec::with_capacity(
            lease,
            "fs-dfo/population-tile-records",
            end.saturating_sub(start),
        ) {
            Ok(rows) => rows,
            Err(error) => return ControlFlow::Break(cx.refuse(TileFailure::Allocation(error))),
        };
        for (position, candidate) in self.candidates[start..end].iter().enumerate() {
            if cx.checkpoint().is_err() {
                return ControlFlow::Break(fs_exec::Cancelled);
            }
            let callback_objectives = (self.objective)(&candidate.decision);
            let actual_dimension = callback_objectives.len();
            let non_finite = !callback_objectives.iter().all(|value| value.is_finite());
            let retained_dimension = if actual_dimension == self.objective_dimension && !non_finite
            {
                actual_dimension
            } else {
                0
            };
            let mut objectives = match LeasedVec::with_capacity(
                lease,
                "fs-dfo/population-objectives",
                retained_dimension,
            ) {
                Ok(objectives) => objectives,
                Err(error) => return ControlFlow::Break(cx.refuse(TileFailure::Allocation(error))),
            };
            for objective in callback_objectives.into_iter().take(retained_dimension) {
                if let Err(error) = objectives.push(objective) {
                    return ControlFlow::Break(cx.refuse(TileFailure::Allocation(error)));
                }
            }
            if let Err(error) = rows.push((
                start + position,
                candidate.identity,
                (objectives, actual_dimension, non_finite),
                cx.stream_key().key128(),
            )) {
                return ControlFlow::Break(cx.refuse(TileFailure::Allocation(error)));
            }
        }
        ControlFlow::Continue(Concat(rows))
    }
}
