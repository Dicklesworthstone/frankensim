//! Bounded, tile-backed population evaluation for generation-oriented engines.
//!
//! The evaluator is deliberately separate from a particular optimizer: NSGA,
//! CMA, and future population methods share the same semantic identities,
//! admission checks, cancellation boundary, and all-or-nothing publication.
//! A callback must be deterministic for replay to be a scientific claim.
//! The output limit bounds retained semantic result bytes; arbitrary allocation
//! performed inside a user callback remains outside this small adapter's claim.

use core::ops::ControlFlow;
use fs_exec::{CancelGate, Cx, Reduce, RunError, RunId, TileKernel, TilePlan, TilePool};
use std::collections::BTreeSet;
use std::sync::Mutex;

const KERNEL: &str = "fs-dfo/population-evaluation-v1";

/// One immutable candidate with a caller-stable semantic identity.
#[derive(Debug, Clone, PartialEq)]
pub struct PopulationCandidate {
    /// Identity stable across worker counts and retries.
    pub identity: u64,
    /// Decision coordinates supplied to the objective callback.
    pub decision: Vec<f64>,
}

/// One accepted objective result, kept in semantic candidate order.
#[derive(Debug, Clone, PartialEq)]
pub struct PopulationEvaluation {
    /// The input candidate identity.
    pub identity: u64,
    /// Objective values in the plan's declared objective order.
    pub objectives: Vec<f64>,
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
    /// Seed used by the pool's logical tile streams.
    pub seed: u64,
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
    /// Semantic generation ordinal.
    pub generation: u64,
    /// Results in the original candidate order.
    pub evaluations: Vec<PopulationEvaluation>,
    /// Exact request identity fields for replay and resume.
    pub provenance: PopulationProvenance,
}

/// Resumable publisher state. It contains only previously committed work.
#[derive(Debug, Clone, PartialEq)]
pub struct PopulationCheckpoint {
    /// Last fully committed generation, if any.
    pub committed: Option<PopulationGeneration>,
    /// Total work committed by successful publication only.
    pub committed_work_units: u64,
}

/// Input or envelope refusal before any tile is launched.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum PopulationRefusal {
    /// Tile width is zero, so no bounded scheduling plan exists.
    ZeroTileWidth,
    /// The supplied population exceeds the caller's cardinality cap.
    PopulationLimit { requested: usize, maximum: usize },
    /// The required complete evaluation count exceeds the work cap.
    WorkLimit { requested: u64, maximum: u64 },
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
    pub fn evaluate<F>(
        &self,
        candidates: &[PopulationCandidate],
        objective_dimension: usize,
        limits: PopulationLimits,
        provenance: PopulationProvenance,
        gate: &CancelGate,
        objective: F,
    ) -> Result<PopulationGeneration, PopulationPublishError>
    where
        F: Fn(&[f64]) -> Vec<f64> + Sync,
    {
        self.preflight(candidates, objective_dimension, limits, provenance)?;
        let kernel = PopulationKernel {
            candidates,
            tile_width: self.tile_width,
            objective: &objective,
        };
        let (result, _) = self.pool.run_declared(&kernel, gate, provenance.run);
        let mut evaluations = result.map_err(PopulationPublishError::Executor)?;
        evaluations.sort_unstable_by_key(|row| row.position);
        let mut accepted = Vec::with_capacity(evaluations.len());
        for row in evaluations {
            if row.objectives.len() != objective_dimension {
                return Err(PopulationPublishError::ObjectiveDimension {
                    identity: row.identity,
                    expected: objective_dimension,
                    actual: row.objectives.len(),
                });
            }
            if !row.objectives.iter().all(|value| value.is_finite()) {
                return Err(PopulationPublishError::NonFiniteObjective {
                    identity: row.identity,
                });
            }
            accepted.push(PopulationEvaluation {
                identity: row.identity,
                objectives: row.objectives,
            });
        }
        Ok(PopulationGeneration {
            generation: provenance.generation,
            evaluations: accepted,
            provenance,
        })
    }

    fn preflight(
        &self,
        candidates: &[PopulationCandidate],
        objective_dimension: usize,
        limits: PopulationLimits,
        provenance: PopulationProvenance,
    ) -> Result<(), PopulationPublishError> {
        if self.tile_width == 0 {
            return Err(PopulationPublishError::Refused(
                PopulationRefusal::ZeroTileWidth,
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
#[derive(Debug, Default)]
pub struct PopulationPublisher {
    state: Mutex<PopulationCheckpoint>,
}

impl PopulationPublisher {
    /// Start with no committed generation.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Resume from a prior immutable checkpoint.
    #[must_use]
    pub fn from_checkpoint(checkpoint: PopulationCheckpoint) -> Self {
        Self {
            state: Mutex::new(checkpoint),
        }
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
    pub fn publish(&self, generation: PopulationGeneration) -> Result<(), PopulationPublishError> {
        let mut state = match self.state.lock() {
            Ok(state) => state,
            Err(poisoned) => poisoned.into_inner(),
        };
        let actual = state.committed.as_ref().map(|previous| previous.generation);
        let expected = actual.map_or(0, |previous| previous.saturating_add(1));
        if generation.generation != expected {
            return Err(PopulationPublishError::GenerationConflict { expected, actual });
        }
        let generation_work = u64::try_from(generation.evaluations.len())
            .map_err(|_| PopulationPublishError::GenerationConflict { expected, actual })?;
        state.committed_work_units = state
            .committed_work_units
            .checked_add(generation_work)
            .ok_or(PopulationPublishError::GenerationConflict { expected, actual })?;
        state.committed = Some(generation);
        Ok(())
    }
}

#[derive(Debug)]
struct TileRows(Vec<TileRow>);

impl Reduce for TileRows {
    fn identity() -> Self {
        Self(Vec::new())
    }
    fn merge(mut self, mut other: Self) -> Self {
        self.0.append(&mut other.0);
        self
    }
}

#[derive(Debug)]
struct TileRow {
    position: usize,
    identity: u64,
    objectives: Vec<f64>,
}

struct PopulationKernel<'a, F> {
    candidates: &'a [PopulationCandidate],
    tile_width: usize,
    objective: &'a F,
}

impl<F> TileKernel for PopulationKernel<'_, F>
where
    F: Fn(&[f64]) -> Vec<f64> + Sync,
{
    type Out = TileRows;

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
        let mut rows = Vec::with_capacity(end.saturating_sub(start));
        for (position, candidate) in self.candidates[start..end].iter().enumerate() {
            if cx.checkpoint().is_err() {
                return ControlFlow::Break(fs_exec::Cancelled);
            }
            let objectives = (self.objective)(&candidate.decision);
            rows.push(TileRow {
                position: start + position,
                identity: candidate.identity,
                objectives,
            });
        }
        ControlFlow::Continue(TileRows(rows))
    }
}
