//! fs-cmaes-viz-wasm — packed browser boundary for fs-dfo's CMA families.
//! Layer: L6.
//!
//! This crate owns transport, not optimization mathematics. Every ask/tell
//! transition delegates to [`fs_dfo::CmaOptimizer`].
//!
//! Contracts every entry inherits (fs-flyer-wasm pattern):
//!
//! - **Typed refusal.** Every output is a schema-versioned numeric success or
//!   refusal packet. Malformed inputs never become optimizer calls.
//! - **Determinism.** Seeded sampling and tie ordering inherit fs-dfo's Philox
//!   and total-order contracts; the adapter adds no randomness.
//! - **Representation honesty.** Snapshots expose only diagnostics the selected
//!   owner actually stores: dense diagonal/spectrum range, separable variances,
//!   or bounded-memory direction norms.
//!
//! No-claims: this is a synchronous browser transport. It does not add BIPOP,
//! cancellation, parallel evaluation, generic objective functions, or dense
//! diagnostic matrices for limited-memory families. Its built-in objectives
//! are the explicitly owner-composed G1 walking experiment in [`g1_walking`]
//! and KUKA household manipulation experiment in [`manipulation`].

// This crate's protocol is deliberately binary64-only. Every integer-to-word
// cast is downstream of exact safe-integer admission or a smaller browser cap;
// every word-to-integer cast is guarded by the inverse finite/range/integrality
// checks. `split_u64` intentionally takes the low and high 32-bit halves.
#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::cast_sign_loss
)]

use fs_dfo::{
    CmaAdmission, CmaAsk, CmaComplexityOrder, CmaConfig, CmaFamily, CmaFamilyError, CmaOptimizer,
    CmaShapeSnapshot, CmaSnapshot, admit_cma,
};
use fs_mbd::robot_models::G1_POLICY_DIMENSION;

#[cfg(feature = "g1-learned")]
pub mod g1_learned;
pub mod g1_walking;
pub mod manipulation;

use g1_walking::{
    G1_LINK_COUNT, G1_LINK_POSE_WORDS, G1_PUSH_END_S, G1_PUSH_PEAK_FORCE_N, G1_PUSH_START_S,
    G1_TERRAIN_AMPLITUDE_M, G1_TERRAIN_WAVENUMBER_RAD_PER_M, G1Challenge, G1Task, G1TraceSample,
    G1WalkingConfig, G1WalkingError, G1WalkingEvaluator, G1WalkingReceipt,
};
use manipulation::{
    ARM_JOINT_COUNT, ARM_LINK_COUNT, ARM_LINK_POSE_WORDS, ARM_POLICY_DIMENSION, ARM_POLICY_KNOTS,
    LIFT_TARGET_M, MAX_DECLARED_OBSTACLES, MIN_GRIPPER_WIDTH_M, ManipulationConfig,
    ManipulationError, ManipulationEvaluator, ManipulationReceipt, ManipulationTask,
    ManipulationTraceSample, OPEN_GRIPPER_WIDTH_M, ObstacleBox, PLACEMENT_TOLERANCE_M,
    manipulation_curriculum_mean, manipulation_max_population,
};

/// Kernel identity returned by the browser capability probe.
pub const KERNEL_VERSION: &str = "fs-cmaes-viz-wasm 0.6.21";
/// Exact binary64 word identifying schema-2 packets (`"CMA2"`).
pub const PACKET_MAGIC: u32 = 0x434d_4132;
/// Packed ask/tell ABI schema.
pub const PACKET_SCHEMA_VERSION: u32 = 2;

/// Input packet containing one optimizer configuration.
pub const PACKET_KIND_CONFIG: u32 = 0;
/// Output packet containing admission plus the initial/current snapshot.
pub const PACKET_KIND_ADMISSION: u32 = 1;
/// Output packet containing one complete candidate population.
pub const PACKET_KIND_ASK: u32 = 2;
/// Input packet containing objective values for the pending population.
pub const PACKET_KIND_TELL: u32 = 3;
/// Output packet containing a post-tell snapshot.
pub const PACKET_KIND_SNAPSHOT: u32 = 4;

/// Successful packet status.
pub const PACKET_STATUS_OK: u32 = 0;
/// Typed-refusal packet status.
pub const PACKET_STATUS_REFUSAL: u32 = 1;

/// Full CMA is cubic and deliberately capped at a browser-honest dimension.
pub const FULL_DIMENSION_LIMIT: usize = 256;
/// Linear and limited-memory families may admit large browser workloads.
pub const SCALABLE_DIMENSION_LIMIT: usize = 100_000;

/// Exact binary64 word identifying G1 curriculum packets (`"G1W7"`).
pub const G1_PACKET_MAGIC: u32 = 0x4731_5737;
/// Packed G1 objective/trace ABI schema.
///
/// Schema 8 makes the walking config self-describing: eleven fixed words
/// followed by a keep-out box count and seven words per box. The receipt
/// gains the maximum body penetration the guard measured, so the browser can
/// report the number the kernel already computed instead of re-deriving it.
pub const G1_PACKET_SCHEMA_VERSION: u32 = 8;
/// Input packet containing fixed walking-experiment controls.
pub const G1_PACKET_KIND_CONFIG: u32 = 0;
/// Output packet describing an admitted evaluator.
pub const G1_PACKET_KIND_ADMISSION: u32 = 1;
/// Output packet containing one candidate receipt.
pub const G1_PACKET_KIND_EVALUATION: u32 = 2;
/// Output packet containing one candidate receipt plus owner-derived poses.
pub const G1_PACKET_KIND_TRACE: u32 = 3;
/// Output packet containing objectives for a complete candidate population.
pub const G1_PACKET_KIND_POPULATION: u32 = 4;

/// Exact binary64 word identifying household-arm packets (`"ARM1"`).
pub const ARM_PACKET_MAGIC: u32 = 0x4152_4d31;
/// Packed household-manipulation objective/trace ABI schema.
///
/// Schema 3 (0.6.14) makes the config packet self-describing: twelve fixed
/// words followed by seven per caller-declared keep-out box. The three added
/// fixed words are an object-mass override and the two Coulomb coefficients;
/// all default to zero, which selects the owner's preset values and
/// reproduces the schema-2 receipts exactly.
pub const ARM_PACKET_SCHEMA_VERSION: u32 = 4;
/// Input packet containing fixed manipulation-experiment controls.
pub const ARM_PACKET_KIND_CONFIG: u32 = 0;
/// Output packet describing an admitted manipulation evaluator and scene.
pub const ARM_PACKET_KIND_ADMISSION: u32 = 1;
/// Output packet containing one manipulation candidate receipt.
pub const ARM_PACKET_KIND_EVALUATION: u32 = 2;
/// Output packet containing one receipt plus owner-derived object/link poses.
pub const ARM_PACKET_KIND_TRACE: u32 = 3;
/// Output packet containing objectives for a complete policy population.
pub const ARM_PACKET_KIND_POPULATION: u32 = 4;

const CONFIG_FIXED_WORDS: usize = 12;
const ASK_FIXED_WORDS: usize = 9;
const TELL_FIXED_WORDS: usize = 6;
const SNAPSHOT_FIXED_WORDS: usize = 31;
const REFUSAL_WORDS: usize = 7;
const MAX_SAFE_INTEGER: f64 = 9_007_199_254_740_991.0;
const BROWSER_PACKET_FIXED_FLOATS: usize = 2 * ASK_FIXED_WORDS + SNAPSHOT_FIXED_WORDS;
const MAX_BROWSER_LIVE_FLOATS: usize = 16 * 1024 * 1024;
const G1_CONFIG_FIXED_WORDS: usize = 12;
const G1_OBSTACLE_WORDS: usize = 8;
/// Largest caller-declared keep-out roster the walking owner admits.
pub const G1_MAX_OBSTACLES: usize = 64;
const G1_ADMISSION_WORDS: usize = 21;
const G1_RECEIPT_WORDS: usize = 29;
const G1_REFUSAL_WORDS: usize = 7;
const G1_TRACE_SAMPLE_WORDS: usize = 3 + G1_LINK_COUNT * G1_LINK_POSE_WORDS;
const G1_MAX_POPULATION: usize = 64;
const ARM_CONFIG_FIXED_WORDS: usize = 12;
const ARM_OBSTACLE_WORDS: usize = 8;
const ARM_ADMISSION_WORDS: usize = 40;
const ARM_RECEIPT_WORDS: usize = 22;
const ARM_REFUSAL_WORDS: usize = 7;
const ARM_TRACE_SAMPLE_WORDS: usize =
    4 + ARM_LINK_POSE_WORDS + ARM_LINK_COUNT * ARM_LINK_POSE_WORDS;

/// Stable numeric refusal codes for schema 2.
#[repr(u32)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PackedRefusalCode {
    /// Packet prefix, length, or integral field is malformed.
    MalformedPacket = 1,
    /// Packet magic or schema is unknown.
    SchemaMismatch = 2,
    /// Family selector is not registered.
    FamilyUnknown = 3,
    /// Dimension is empty or exceeds the scalable adapter limit.
    DimensionInvalid = 4,
    /// Full CMA exceeds its honest cubic browser limit.
    FullDimensionLimit = 5,
    /// Population is below four or unrepresentable.
    PopulationInvalid = 6,
    /// Memory is zero when explicit, inapplicable, or unrepresentable.
    MemoryInvalid = 7,
    /// No complete generation fits, or budget is unrepresentable.
    BudgetInvalid = 8,
    /// Seed words are not exact unsigned 32-bit integers.
    SeedInvalid = 9,
    /// Initial step size is non-finite or non-positive.
    SigmaInvalid = 10,
    /// Initial mean contains a non-finite coordinate.
    MeanInvalid = 11,
    /// Checked optimizer shape arithmetic overflowed.
    ShapeOverflow = 12,
    /// Exact Philox counter admission failed.
    RandomCounterOverflow = 13,
    /// Dense eigensolver admission refused full CMA.
    DenseWorkRefused = 14,
    /// Conservative browser float envelope exceeds the adapter cap.
    BrowserMemoryRefused = 15,
    /// An ask is already outstanding.
    AskAlreadyPending = 16,
    /// No complete population remains in the admitted budget.
    BudgetExhausted = 17,
    /// Tell has no outstanding ask.
    NoPendingAsk = 18,
    /// Tell names a different generation.
    GenerationMismatch = 19,
    /// Tell has the wrong objective count.
    ObjectiveCount = 20,
    /// One supplied objective is non-finite.
    NonFiniteObjective = 21,
    /// The owner refused an internal batch identity.
    BatchMismatch = 22,
    /// An owner numerical update left the finite domain.
    NumericalFailure = 23,
}

#[derive(Debug, Clone, Copy)]
struct PackedRefusal {
    code: PackedRefusalCode,
    generation: Option<u64>,
}

impl PackedRefusal {
    const fn new(code: PackedRefusalCode) -> Self {
        Self {
            code,
            generation: None,
        }
    }

    const fn at_generation(code: PackedRefusalCode, generation: u64) -> Self {
        Self {
            code,
            generation: Some(generation),
        }
    }
}

/// Stateful packed boundary over exactly one fs-dfo optimizer.
///
/// Invalid configurations also produce a session value; its `receipt`, `ask`,
/// and `tell` methods return typed refusal packets instead of throwing across
/// wasm-bindgen.
pub struct PackedCmaSession {
    optimizer: Option<CmaOptimizer>,
    admission: Option<CmaAdmission>,
    pending: Option<CmaAsk>,
    evaluations: usize,
    creation_refusal: Option<PackedRefusal>,
}

impl PackedCmaSession {
    /// Parse and admit a schema-2 configuration packet.
    #[must_use]
    pub fn new(config_packet: &[f64]) -> Self {
        let config = match parse_config(config_packet) {
            Ok(config) => config,
            Err(refusal) => return Self::refused(refusal),
        };
        let admission = match admit_cma(&config) {
            Ok(admission) => admission,
            Err(error) => return Self::refused(owner_refusal(&error)),
        };
        let owner_floats = admission
            .complexity
            .persistent_scalars
            .checked_add(admission.complexity.pending_generation_scalars)
            .and_then(|value| value.checked_add(admission.complexity.update_workspace_scalars));
        // During an ask return, the owner retains its pending population while
        // Rust flattens that population into a packet and wasm-bindgen copies
        // the packet into a JavaScript Float64Array. Include both transport
        // copies plus a conservative snapshot payload and packet headers
        // instead of pretending the owner allocation is the whole browser
        // live set.
        let browser_live_floats = admission
            .dimension
            .checked_mul(admission.population_size)
            .and_then(|candidate_floats| candidate_floats.checked_mul(2))
            .and_then(|transport_floats| {
                admission
                    .dimension
                    .checked_mul(3)
                    .and_then(|snapshot_floats| {
                        snapshot_floats.checked_add(admission.complexity.memory_capacity)
                    })
                    .and_then(|snapshot_floats| {
                        snapshot_floats.checked_add(BROWSER_PACKET_FIXED_FLOATS)
                    })
                    .and_then(|snapshot_floats| transport_floats.checked_add(snapshot_floats))
            })
            .and_then(|transport_and_snapshot| {
                owner_floats.and_then(|owner| owner.checked_add(transport_and_snapshot))
            });
        if browser_live_floats.is_none_or(|value| value > MAX_BROWSER_LIVE_FLOATS) {
            return Self::refused(PackedRefusal::new(PackedRefusalCode::BrowserMemoryRefused));
        }
        match CmaOptimizer::new(config) {
            Ok(optimizer) => Self {
                optimizer: Some(optimizer),
                admission: Some(admission),
                pending: None,
                evaluations: 0,
                creation_refusal: None,
            },
            Err(error) => Self::refused(owner_refusal(&error)),
        }
    }

    const fn refused(refusal: PackedRefusal) -> Self {
        Self {
            optimizer: None,
            admission: None,
            pending: None,
            evaluations: 0,
            creation_refusal: Some(refusal),
        }
    }

    /// Return admission, initial state, and exact complexity diagnostics.
    #[must_use]
    pub fn receipt_packet(&self) -> Vec<f64> {
        if let Some(refusal) = self.creation_refusal {
            return refusal_packet(PACKET_KIND_ADMISSION, refusal);
        }
        let (Some(optimizer), Some(admission)) = (self.optimizer.as_ref(), self.admission) else {
            return refusal_packet(
                PACKET_KIND_ADMISSION,
                PackedRefusal::new(PackedRefusalCode::NumericalFailure),
            );
        };
        snapshot_packet(PACKET_KIND_ADMISSION, admission, &optimizer.snapshot())
    }

    /// Ask for one complete population and retain its opaque owner token.
    #[must_use]
    pub fn ask_packet(&mut self) -> Vec<f64> {
        if let Some(refusal) = self.creation_refusal {
            return refusal_packet(PACKET_KIND_ASK, refusal);
        }
        if let Some(batch) = &self.pending {
            return refusal_packet(
                PACKET_KIND_ASK,
                PackedRefusal::at_generation(
                    PackedRefusalCode::AskAlreadyPending,
                    batch.generation(),
                ),
            );
        }
        let (Some(optimizer), Some(admission)) = (self.optimizer.as_mut(), self.admission) else {
            return refusal_packet(
                PACKET_KIND_ASK,
                PackedRefusal::new(PackedRefusalCode::NumericalFailure),
            );
        };
        match optimizer.ask() {
            Ok(batch) => {
                let packet = ask_packet(admission, self.evaluations, &batch);
                self.pending = Some(batch);
                packet
            }
            Err(error) => refusal_packet(PACKET_KIND_ASK, owner_refusal(&error)),
        }
    }

    /// Tell objective values for the pending population.
    ///
    /// Every parse or owner refusal leaves the pending batch intact, so callers
    /// can repair and retry the same generation.
    #[must_use]
    pub fn tell_packet(&mut self, objective_packet: &[f64]) -> Vec<f64> {
        if let Some(refusal) = self.creation_refusal {
            return refusal_packet(PACKET_KIND_TELL, refusal);
        }
        let Some(admission) = self.admission else {
            return refusal_packet(
                PACKET_KIND_TELL,
                PackedRefusal::new(PackedRefusalCode::NumericalFailure),
            );
        };
        let Some(batch) = self.pending.as_ref() else {
            return refusal_packet(
                PACKET_KIND_TELL,
                PackedRefusal::new(PackedRefusalCode::NoPendingAsk),
            );
        };
        let expected_generation = batch.generation();
        let (generation, objectives) =
            match parse_tell(objective_packet, admission, expected_generation) {
                Ok(parsed) => parsed,
                Err(refusal) => return refusal_packet(PACKET_KIND_TELL, refusal),
            };
        if generation != batch.generation() {
            return refusal_packet(
                PACKET_KIND_TELL,
                PackedRefusal::at_generation(
                    PackedRefusalCode::GenerationMismatch,
                    batch.generation(),
                ),
            );
        }
        let Some(optimizer) = self.optimizer.as_mut() else {
            return refusal_packet(
                PACKET_KIND_TELL,
                PackedRefusal::new(PackedRefusalCode::NumericalFailure),
            );
        };
        match optimizer.tell(batch, objectives) {
            Ok(snapshot) => {
                self.pending = None;
                self.evaluations = snapshot.evaluations;
                snapshot_packet(PACKET_KIND_SNAPSHOT, admission, &snapshot)
            }
            Err(error) => refusal_packet(PACKET_KIND_TELL, owner_refusal(&error)),
        }
    }
}

/// Stable refusal codes for the owner-composed G1 walking boundary.
#[repr(u32)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum G1PackedRefusalCode {
    /// Packet prefix, length, or integral field is malformed.
    MalformedPacket = 1,
    /// Packet magic or schema is unknown.
    SchemaMismatch = 2,
    /// A fixed rollout control is outside the browser experiment domain.
    InvalidConfig = 3,
    /// Policy coordinate count does not match the catalog-owned 5,040-D map.
    ParameterCount = 4,
    /// A policy coordinate is NaN or infinite.
    NonFiniteParameter = 5,
    /// The articulated-body owner refused the rollout.
    Robot = 6,
    /// The policy owner refused an observation or output.
    Policy = 7,
    /// The normal-contact owner refused the rollout.
    Contact = 8,
    /// The friction owner refused the rollout.
    Friction = 9,
    /// The Lie-group time owner refused the rollout.
    Time = 10,
    /// The geometry owner refused the rollout.
    Geometry = 11,
    /// A sphere/plane request returned an impossible receipt kind.
    UnexpectedContactReceipt = 12,
    /// A completed rollout produced a non-finite objective.
    NonFiniteObjective = 13,
    /// A flat population is empty, too large, or not candidate-aligned.
    PopulationInvalid = 14,
    /// Checked packet-size arithmetic overflowed.
    ShapeOverflow = 15,
}

#[derive(Debug, Clone, Copy)]
struct G1PackedRefusal {
    code: G1PackedRefusalCode,
    detail: Option<usize>,
}

impl G1PackedRefusal {
    const fn new(code: G1PackedRefusalCode) -> Self {
        Self { code, detail: None }
    }

    const fn with_detail(code: G1PackedRefusalCode, detail: usize) -> Self {
        Self {
            code,
            detail: Some(detail),
        }
    }
}

/// Reusable browser boundary for the fixed G1 walking experiment.
///
/// The source-bound model and constitutive owners are constructed once. Each
/// call then evaluates either one admitted 5,040-D policy or one flat complete
/// population without rebuilding the catalog.
pub struct PackedG1WalkingEvaluator {
    evaluator: Option<G1WalkingEvaluator>,
    creation_refusal: Option<G1PackedRefusal>,
}

impl PackedG1WalkingEvaluator {
    /// Parse and admit one packed G1 experiment configuration.
    #[must_use]
    pub fn new(config_packet: &[f64]) -> Self {
        let config = match parse_g1_config(config_packet) {
            Ok(config) => config,
            Err(refusal) => return Self::refused(refusal),
        };
        match G1WalkingEvaluator::new(config) {
            Ok(evaluator) => Self {
                evaluator: Some(evaluator),
                creation_refusal: None,
            },
            Err(error) => Self::refused(g1_owner_refusal(&error)),
        }
    }

    const fn refused(refusal: G1PackedRefusal) -> Self {
        Self {
            evaluator: None,
            creation_refusal: Some(refusal),
        }
    }

    /// Return exact fixed controls and render-layout dimensions.
    #[must_use]
    pub fn receipt_packet(&self) -> Vec<f64> {
        if let Some(refusal) = self.creation_refusal {
            return g1_refusal_packet(G1_PACKET_KIND_ADMISSION, refusal);
        }
        let Some(evaluator) = self.evaluator.as_ref() else {
            return g1_refusal_packet(
                G1_PACKET_KIND_ADMISSION,
                G1PackedRefusal::new(G1PackedRefusalCode::Robot),
            );
        };
        let config = evaluator.config();
        let mut packet = g1_success_header(G1_PACKET_KIND_ADMISSION, G1_ADMISSION_WORDS);
        packet.extend_from_slice(&[
            G1_POLICY_DIMENSION as f64,
            G1_LINK_COUNT as f64,
            G1_LINK_POSE_WORDS as f64,
            G1_TRACE_SAMPLE_WORDS as f64,
            config.step_s,
            config.duration_s,
            config.target_forward_speed_m_per_s,
            config.gait_frequency_hz,
            config.trace_stride as f64,
            f64::from(config.task as u32),
            f64::from(config.challenge as u32),
            G1_TERRAIN_AMPLITUDE_M,
            G1_TERRAIN_WAVENUMBER_RAD_PER_M,
            G1_PUSH_START_S,
            G1_PUSH_END_S,
            G1_PUSH_PEAK_FORCE_N,
        ]);
        debug_assert_eq!(packet.len(), G1_ADMISSION_WORDS);
        packet
    }

    /// Evaluate one 5,040-D policy without retaining a trajectory.
    #[must_use]
    pub fn evaluate_packet(&self, parameters: &[f64]) -> Vec<f64> {
        if let Some(refusal) = self.creation_refusal {
            return g1_refusal_packet(G1_PACKET_KIND_EVALUATION, refusal);
        }
        let Some(evaluator) = self.evaluator.as_ref() else {
            return g1_refusal_packet(
                G1_PACKET_KIND_EVALUATION,
                G1PackedRefusal::new(G1PackedRefusalCode::Robot),
            );
        };
        match evaluator.evaluate(parameters) {
            Ok(receipt) => g1_receipt_packet(G1_PACKET_KIND_EVALUATION, &receipt, false),
            Err(error) => g1_refusal_packet(G1_PACKET_KIND_EVALUATION, g1_owner_refusal(&error)),
        }
    }

    /// Evaluate a flat row-major population and return one objective per row.
    #[must_use]
    pub fn evaluate_population_packet(&self, parameters: &[f64]) -> Vec<f64> {
        if let Some(refusal) = self.creation_refusal {
            return g1_refusal_packet(G1_PACKET_KIND_POPULATION, refusal);
        }
        let Some(evaluator) = self.evaluator.as_ref() else {
            return g1_refusal_packet(
                G1_PACKET_KIND_POPULATION,
                G1PackedRefusal::new(G1PackedRefusalCode::Robot),
            );
        };
        if parameters.is_empty() || !parameters.len().is_multiple_of(G1_POLICY_DIMENSION) {
            return g1_refusal_packet(
                G1_PACKET_KIND_POPULATION,
                G1PackedRefusal::new(G1PackedRefusalCode::PopulationInvalid),
            );
        }
        let population = parameters.len() / G1_POLICY_DIMENSION;
        if population > G1_MAX_POPULATION {
            return g1_refusal_packet(
                G1_PACKET_KIND_POPULATION,
                G1PackedRefusal::new(G1PackedRefusalCode::PopulationInvalid),
            );
        }
        let Some(total_words) = 6usize.checked_add(population) else {
            return g1_refusal_packet(
                G1_PACKET_KIND_POPULATION,
                G1PackedRefusal::new(G1PackedRefusalCode::ShapeOverflow),
            );
        };
        let mut packet = g1_success_header(G1_PACKET_KIND_POPULATION, total_words);
        packet.push(population as f64);
        let (policies, remainder) = parameters.as_chunks::<G1_POLICY_DIMENSION>();
        debug_assert!(remainder.is_empty());
        for (candidate, policy) in policies.iter().enumerate() {
            match evaluator.evaluate(policy) {
                Ok(receipt) => packet.push(receipt.objective),
                Err(error) => {
                    let refusal = g1_owner_refusal(&error);
                    return g1_refusal_packet(
                        G1_PACKET_KIND_POPULATION,
                        G1PackedRefusal::with_detail(refusal.code, candidate),
                    );
                }
            }
        }
        debug_assert_eq!(packet.len(), total_words);
        packet
    }

    /// Evaluate one policy and retain decimated owner-derived link poses.
    #[must_use]
    pub fn trace_packet(&self, parameters: &[f64]) -> Vec<f64> {
        if let Some(refusal) = self.creation_refusal {
            return g1_refusal_packet(G1_PACKET_KIND_TRACE, refusal);
        }
        let Some(evaluator) = self.evaluator.as_ref() else {
            return g1_refusal_packet(
                G1_PACKET_KIND_TRACE,
                G1PackedRefusal::new(G1PackedRefusalCode::Robot),
            );
        };
        match evaluator.trace(parameters) {
            Ok(receipt) => g1_receipt_packet(G1_PACKET_KIND_TRACE, &receipt, true),
            Err(error) => g1_refusal_packet(G1_PACKET_KIND_TRACE, g1_owner_refusal(&error)),
        }
    }
}

/// Stable refusal codes for the owner-composed household-arm boundary.
#[repr(u32)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ArmPackedRefusalCode {
    /// Packet prefix, length, or integral field is malformed.
    MalformedPacket = 1,
    /// Packet magic or schema is unknown.
    SchemaMismatch = 2,
    /// A fixed rollout control or task selector is outside the browser domain.
    InvalidConfig = 3,
    /// Policy coordinate count does not match the declared 128-D knot map.
    ParameterCount = 4,
    /// A policy coordinate is NaN or infinite.
    NonFiniteParameter = 5,
    /// The source-bound articulated owner refused the rollout.
    Robot = 6,
    /// A canonical Lie-group operation refused the rollout.
    Geometry = 7,
    /// The compliant normal-contact owner refused the rollout.
    Contact = 8,
    /// The dry-friction owner refused the rollout.
    Friction = 9,
    /// A sphere/plane request returned an impossible receipt kind.
    UnexpectedContactReceipt = 10,
    /// A completed rollout produced a non-finite objective or receipt.
    NonFiniteObjective = 11,
    /// A flat population is empty, too large, or not candidate-aligned.
    PopulationInvalid = 12,
    /// Checked packet-size arithmetic overflowed.
    ShapeOverflow = 13,
    /// The certified convex-query owner refused an envelope or query.
    Query = 14,
}

#[derive(Debug, Clone, Copy)]
struct ArmPackedRefusal {
    code: ArmPackedRefusalCode,
    detail: Option<usize>,
}

impl ArmPackedRefusal {
    const fn new(code: ArmPackedRefusalCode) -> Self {
        Self { code, detail: None }
    }

    const fn with_detail(code: ArmPackedRefusalCode, detail: usize) -> Self {
        Self {
            code,
            detail: Some(detail),
        }
    }
}

/// Reusable browser boundary for the fixed KUKA household experiment.
///
/// The pinned articulated model and constitutive owners are constructed once;
/// calls then evaluate complete 128-D policies without rebuilding the catalog.
pub struct PackedManipulationEvaluator {
    evaluator: Option<ManipulationEvaluator>,
    task: Option<ManipulationTask>,
    creation_refusal: Option<ArmPackedRefusal>,
}

impl PackedManipulationEvaluator {
    /// Parse and admit one packed household-manipulation configuration.
    #[must_use]
    pub fn new(config_packet: &[f64]) -> Self {
        let config = match parse_arm_config(config_packet) {
            Ok(config) => config,
            Err(refusal) => return Self::refused(refusal),
        };
        let task = config.task;
        match ManipulationEvaluator::new(config) {
            Ok(evaluator) => Self {
                evaluator: Some(evaluator),
                task: Some(task),
                creation_refusal: None,
            },
            Err(error) => Self::refused(arm_owner_refusal(&error)),
        }
    }

    const fn refused(refusal: ArmPackedRefusal) -> Self {
        Self {
            evaluator: None,
            task: None,
            creation_refusal: Some(refusal),
        }
    }

    /// Return exact controls, scene geometry, and render-layout dimensions.
    #[must_use]
    pub fn receipt_packet(&self) -> Vec<f64> {
        if let Some(refusal) = self.creation_refusal {
            return arm_refusal_packet(ARM_PACKET_KIND_ADMISSION, refusal);
        }
        let Some(evaluator) = self.evaluator.as_ref() else {
            return arm_refusal_packet(
                ARM_PACKET_KIND_ADMISSION,
                ArmPackedRefusal::new(ArmPackedRefusalCode::Robot),
            );
        };
        let config = evaluator.config();
        let scene = evaluator.scene();
        let mut packet = arm_success_header(ARM_PACKET_KIND_ADMISSION, ARM_ADMISSION_WORDS);
        packet.extend_from_slice(&[
            ARM_POLICY_DIMENSION as f64,
            ARM_JOINT_COUNT as f64,
            ARM_POLICY_KNOTS as f64,
            ARM_LINK_COUNT as f64,
            ARM_LINK_POSE_WORDS as f64,
            ARM_TRACE_SAMPLE_WORDS as f64,
            config.step_s,
            config.duration_s,
            config.trace_stride as f64,
            f64::from(config.task as u32),
            MIN_GRIPPER_WIDTH_M,
            OPEN_GRIPPER_WIDTH_M,
            PLACEMENT_TOLERANCE_M,
            LIFT_TARGET_M,
            scene.object_mass_kg,
            scene.object_dimensions_m.x,
            scene.object_dimensions_m.y,
            scene.object_dimensions_m.z,
            scene.grasp_half_width_m,
            scene.initial_object_position_m.x,
            scene.initial_object_position_m.y,
            scene.initial_object_position_m.z,
            scene.goal_object_position_m.x,
            scene.goal_object_position_m.y,
            scene.goal_object_position_m.z,
            scene.support_height_m,
            scene.obstacle_center_m.x,
            scene.obstacle_center_m.y,
            scene.obstacle_center_m.z,
            scene.obstacle_half_extents_m.x,
            scene.obstacle_half_extents_m.y,
            scene.obstacle_half_extents_m.z,
            // Schema 3: the effective interface the rollout actually ran with,
            // so a caller override can never be read back as the owner default.
            scene.static_friction_mu,
            scene.kinetic_friction_mu,
            scene.extra_obstacle_count as f64,
        ]);
        debug_assert_eq!(packet.len(), ARM_ADMISSION_WORDS);
        packet
    }

    /// Return the disclosed source-feasible 128-D curriculum mean.
    #[must_use]
    pub fn curriculum_policy_mean(&self) -> Vec<f64> {
        self.task
            .map_or_else(Vec::new, |task| manipulation_curriculum_mean(task).to_vec())
    }

    /// Evaluate one 128-D policy without retaining a trajectory.
    #[must_use]
    pub fn evaluate_packet(&self, parameters: &[f64]) -> Vec<f64> {
        if let Some(refusal) = self.creation_refusal {
            return arm_refusal_packet(ARM_PACKET_KIND_EVALUATION, refusal);
        }
        let Some(evaluator) = self.evaluator.as_ref() else {
            return arm_refusal_packet(
                ARM_PACKET_KIND_EVALUATION,
                ArmPackedRefusal::new(ArmPackedRefusalCode::Robot),
            );
        };
        match evaluator.evaluate(parameters) {
            Ok(receipt) => arm_receipt_packet(ARM_PACKET_KIND_EVALUATION, &receipt, false),
            Err(error) => arm_refusal_packet(ARM_PACKET_KIND_EVALUATION, arm_owner_refusal(&error)),
        }
    }

    /// Evaluate a flat row-major population and return one objective per row.
    #[must_use]
    pub fn evaluate_population_packet(&self, parameters: &[f64]) -> Vec<f64> {
        if let Some(refusal) = self.creation_refusal {
            return arm_refusal_packet(ARM_PACKET_KIND_POPULATION, refusal);
        }
        let Some(evaluator) = self.evaluator.as_ref() else {
            return arm_refusal_packet(
                ARM_PACKET_KIND_POPULATION,
                ArmPackedRefusal::new(ArmPackedRefusalCode::Robot),
            );
        };
        if parameters.is_empty() || !parameters.len().is_multiple_of(ARM_POLICY_DIMENSION) {
            return arm_refusal_packet(
                ARM_PACKET_KIND_POPULATION,
                ArmPackedRefusal::new(ArmPackedRefusalCode::PopulationInvalid),
            );
        }
        let population = parameters.len() / ARM_POLICY_DIMENSION;
        if population > manipulation_max_population() {
            return arm_refusal_packet(
                ARM_PACKET_KIND_POPULATION,
                ArmPackedRefusal::new(ArmPackedRefusalCode::PopulationInvalid),
            );
        }
        let Some(total_words) = 6usize.checked_add(population) else {
            return arm_refusal_packet(
                ARM_PACKET_KIND_POPULATION,
                ArmPackedRefusal::new(ArmPackedRefusalCode::ShapeOverflow),
            );
        };
        let mut packet = arm_success_header(ARM_PACKET_KIND_POPULATION, total_words);
        packet.push(population as f64);
        let (policies, remainder) = parameters.as_chunks::<ARM_POLICY_DIMENSION>();
        debug_assert!(remainder.is_empty());
        for (candidate, policy) in policies.iter().enumerate() {
            match evaluator.evaluate(policy) {
                Ok(receipt) => packet.push(receipt.objective),
                Err(error) => {
                    let refusal = arm_owner_refusal(&error);
                    return arm_refusal_packet(
                        ARM_PACKET_KIND_POPULATION,
                        ArmPackedRefusal::with_detail(refusal.code, candidate),
                    );
                }
            }
        }
        debug_assert_eq!(packet.len(), total_words);
        packet
    }

    /// Evaluate one policy and retain decimated owner-derived poses.
    #[must_use]
    pub fn trace_packet(&self, parameters: &[f64]) -> Vec<f64> {
        if let Some(refusal) = self.creation_refusal {
            return arm_refusal_packet(ARM_PACKET_KIND_TRACE, refusal);
        }
        let Some(evaluator) = self.evaluator.as_ref() else {
            return arm_refusal_packet(
                ARM_PACKET_KIND_TRACE,
                ArmPackedRefusal::new(ArmPackedRefusalCode::Robot),
            );
        };
        match evaluator.trace(parameters) {
            Ok(receipt) => arm_receipt_packet(ARM_PACKET_KIND_TRACE, &receipt, true),
            Err(error) => arm_refusal_packet(ARM_PACKET_KIND_TRACE, arm_owner_refusal(&error)),
        }
    }
}

/// Parse a schema-8 G1 walking config packet.
///
/// Layout: `[magic, schema, kind, wordCount, step_s, duration_s,
/// target_forward_speed, gait_frequency, trace_stride, task, challenge,
/// obstacleCount]` followed by `obstacleCount` groups of
/// `[cx, cy, cz, hx, hy, hz, yaw_rad, role]` in the owner world frame (z up,
/// yaw about +Z), where role 0 is keep-out (nothing may enter) and role 1 is
/// support (things rest on it but may not sink through). `wordCount` is
/// self-describing and must equal the packet length. An empty roster leaves
/// the rollout identical to schema 7.
fn parse_g1_config(packet: &[f64]) -> Result<G1WalkingConfig, G1PackedRefusal> {
    if packet.len() < G1_CONFIG_FIXED_WORDS {
        return Err(G1PackedRefusal::new(G1PackedRefusalCode::MalformedPacket));
    }
    if exact_u32(packet[0]) != Some(G1_PACKET_MAGIC)
        || exact_u32(packet[1]) != Some(G1_PACKET_SCHEMA_VERSION)
    {
        return Err(G1PackedRefusal::new(G1PackedRefusalCode::SchemaMismatch));
    }
    if exact_u32(packet[2]) != Some(G1_PACKET_KIND_CONFIG)
        || exact_usize(packet[3]) != Some(packet.len())
    {
        return Err(G1PackedRefusal::new(G1PackedRefusalCode::MalformedPacket));
    }
    let obstacle_count = exact_usize(packet[11])
        .filter(|value| *value <= G1_MAX_OBSTACLES)
        .ok_or(G1PackedRefusal::new(G1PackedRefusalCode::InvalidConfig))?;
    if packet.len() != G1_CONFIG_FIXED_WORDS + obstacle_count * G1_OBSTACLE_WORDS {
        return Err(G1PackedRefusal::new(G1PackedRefusalCode::MalformedPacket));
    }
    let mut obstacles = Vec::with_capacity(obstacle_count);
    for index in 0..obstacle_count {
        let base = G1_CONFIG_FIXED_WORDS + index * G1_OBSTACLE_WORDS;
        let words = &packet[base..base + G1_OBSTACLE_WORDS];
        if words.iter().any(|value| !value.is_finite()) {
            return Err(G1PackedRefusal::new(G1PackedRefusalCode::InvalidConfig));
        }
        let role = match exact_u32(words[7]) {
            Some(0) => fs_scene::BodyRole::KeepOut,
            Some(1) => fs_scene::BodyRole::Support,
            _ => return Err(G1PackedRefusal::new(G1PackedRefusalCode::InvalidConfig)),
        };
        obstacles.push(crate::g1_walking::ObstacleBox {
            center_m: [words[0], words[1], words[2]],
            half_extents_m: [words[3], words[4], words[5]],
            yaw_rad: words[6],
            role,
        });
    }
    let trace_stride = exact_usize(packet[8])
        .filter(|value| (1..=1_000).contains(value))
        .ok_or(G1PackedRefusal::new(G1PackedRefusalCode::InvalidConfig))?;
    let task = match exact_u32(packet[9]) {
        Some(0) => G1Task::Balance,
        Some(1) => G1Task::Stepping,
        Some(2) => G1Task::Walking,
        _ => return Err(G1PackedRefusal::new(G1PackedRefusalCode::InvalidConfig)),
    };
    let challenge = match exact_u32(packet[10]) {
        Some(0) => G1Challenge::Flat,
        Some(1) => G1Challenge::TerrainAndPush,
        _ => return Err(G1PackedRefusal::new(G1PackedRefusalCode::InvalidConfig)),
    };
    let config = G1WalkingConfig {
        obstacles,
        task,
        challenge,
        step_s: packet[4],
        duration_s: packet[5],
        target_forward_speed_m_per_s: packet[6],
        gait_frequency_hz: packet[7],
        trace_stride,
    };
    if !(config.step_s.is_finite()
        && (1.0 / 480.0..=1.0 / 30.0).contains(&config.step_s)
        && config.duration_s.is_finite()
        && (config.step_s..=4.0).contains(&config.duration_s)
        && config.target_forward_speed_m_per_s.is_finite()
        && (0.0..=2.0).contains(&config.target_forward_speed_m_per_s)
        && config.gait_frequency_hz.is_finite()
        && (0.25..=4.0).contains(&config.gait_frequency_hz))
    {
        return Err(G1PackedRefusal::new(G1PackedRefusalCode::InvalidConfig));
    }
    Ok(config)
}

fn g1_owner_refusal(error: &G1WalkingError) -> G1PackedRefusal {
    let code = match error {
        G1WalkingError::InvalidConfig { .. } => G1PackedRefusalCode::InvalidConfig,
        G1WalkingError::Robot(_) => G1PackedRefusalCode::Robot,
        G1WalkingError::Policy(fs_mbd::robot_models::G1PolicyError::ParameterCount { .. }) => {
            G1PackedRefusalCode::ParameterCount
        }
        G1WalkingError::Policy(fs_mbd::robot_models::G1PolicyError::NonFiniteParameter {
            ..
        }) => G1PackedRefusalCode::NonFiniteParameter,
        G1WalkingError::Policy(_) => G1PackedRefusalCode::Policy,
        G1WalkingError::Contact(_) => G1PackedRefusalCode::Contact,
        G1WalkingError::Friction(_) => G1PackedRefusalCode::Friction,
        G1WalkingError::Time(_) => G1PackedRefusalCode::Time,
        G1WalkingError::Geometry(_) => G1PackedRefusalCode::Geometry,
        G1WalkingError::UnexpectedContactReceipt => G1PackedRefusalCode::UnexpectedContactReceipt,
        G1WalkingError::NonFiniteObjective => G1PackedRefusalCode::NonFiniteObjective,
    };
    G1PackedRefusal::new(code)
}

fn g1_success_header(kind: u32, total_words: usize) -> Vec<f64> {
    vec![
        f64::from(G1_PACKET_MAGIC),
        f64::from(G1_PACKET_SCHEMA_VERSION),
        f64::from(PACKET_STATUS_OK),
        f64::from(kind),
        total_words as f64,
    ]
}

fn g1_refusal_packet(kind: u32, refusal: G1PackedRefusal) -> Vec<f64> {
    vec![
        f64::from(G1_PACKET_MAGIC),
        f64::from(G1_PACKET_SCHEMA_VERSION),
        f64::from(PACKET_STATUS_REFUSAL),
        f64::from(kind),
        G1_REFUSAL_WORDS as f64,
        f64::from(refusal.code as u32),
        refusal.detail.map_or(f64::NAN, |value| value as f64),
    ]
}

fn g1_receipt_packet(kind: u32, receipt: &G1WalkingReceipt, include_trace: bool) -> Vec<f64> {
    let trace_words = if include_trace {
        receipt
            .trace
            .len()
            .checked_mul(G1_TRACE_SAMPLE_WORDS)
            .and_then(|words| words.checked_add(1))
    } else {
        Some(0)
    };
    let Some(total_words) = trace_words.and_then(|words| G1_RECEIPT_WORDS.checked_add(words))
    else {
        return g1_refusal_packet(
            kind,
            G1PackedRefusal::new(G1PackedRefusalCode::ShapeOverflow),
        );
    };
    let mut packet = g1_success_header(kind, total_words);
    packet.extend_from_slice(&[
        receipt.objective,
        receipt.distance_m,
        receipt.speed_error_integral,
        receipt.actuator_work_j,
        receipt.slip_integral,
        receipt.posture_integral,
        receipt.joint_limit_integral,
        receipt.impact_integral,
        receipt.backward_distance_m,
        receipt.lateral_error_integral,
        receipt.heading_error_integral,
        receipt.contact_schedule_mismatch_integral,
        receipt.swing_clearance_error_integral,
        receipt.single_support_s,
        receipt.double_support_s,
        receipt.flight_s,
        receipt.push_impulse_n_s,
        receipt.recovery_time_s,
        receipt.minimum_base_height_m,
        receipt.maximum_tilt_sine,
        receipt.maximum_abs_terrain_height_m,
        receipt.completed_steps as f64,
        f64::from(receipt.termination_reason as u32),
        // Schema 8: the deepest body-sphere penetration the obstacle guard
        // measured. Zero with an empty roster; the browser reports this
        // instead of re-deriving contact from the rendered poses.
        receipt.maximum_body_penetration_m,
    ]);
    if include_trace {
        packet.push(receipt.trace.len() as f64);
        for sample in &receipt.trace {
            append_g1_trace_sample(&mut packet, sample);
        }
    }
    debug_assert_eq!(packet.len(), total_words);
    packet
}

fn append_g1_trace_sample(packet: &mut Vec<f64>, sample: &G1TraceSample) {
    packet.extend_from_slice(&[
        sample.time_s,
        f64::from(u8::from(sample.foot_contact[0])),
        f64::from(u8::from(sample.foot_contact[1])),
    ]);
    for pose in &sample.link_pose {
        packet.extend_from_slice(pose);
    }
}

/// Parse a schema-3 household-arm config packet.
///
/// Layout: `[magic, schema, kind, wordCount, step_s, duration_s,
/// trace_stride, task, object_mass_kg, static_mu, kinetic_mu, obstacleCount]`
/// followed by `obstacleCount` groups of
/// `[cx, cy, cz, hx, hy, hz, yaw_rad]` in world coordinates. `wordCount` is
/// self-describing and must equal the packet length. A zero in any of the
/// three override words selects the owner's preset value.
fn parse_arm_config(packet: &[f64]) -> Result<ManipulationConfig, ArmPackedRefusal> {
    if packet.len() < ARM_CONFIG_FIXED_WORDS {
        return Err(ArmPackedRefusal::new(ArmPackedRefusalCode::MalformedPacket));
    }
    if exact_u32(packet[0]) != Some(ARM_PACKET_MAGIC)
        || exact_u32(packet[1]) != Some(ARM_PACKET_SCHEMA_VERSION)
    {
        return Err(ArmPackedRefusal::new(ArmPackedRefusalCode::SchemaMismatch));
    }
    if exact_u32(packet[2]) != Some(ARM_PACKET_KIND_CONFIG)
        || exact_usize(packet[3]) != Some(packet.len())
    {
        return Err(ArmPackedRefusal::new(ArmPackedRefusalCode::MalformedPacket));
    }
    let obstacle_count = exact_usize(packet[11])
        .filter(|value| *value <= MAX_DECLARED_OBSTACLES)
        .ok_or(ArmPackedRefusal::new(ArmPackedRefusalCode::InvalidConfig))?;
    if packet.len() != ARM_CONFIG_FIXED_WORDS + obstacle_count * ARM_OBSTACLE_WORDS {
        return Err(ArmPackedRefusal::new(ArmPackedRefusalCode::MalformedPacket));
    }
    let trace_stride = exact_usize(packet[6])
        .filter(|value| (1..=1_000).contains(value))
        .ok_or(ArmPackedRefusal::new(ArmPackedRefusalCode::InvalidConfig))?;
    let task = match exact_u32(packet[7]) {
        Some(0) => ManipulationTask::KitchenMug,
        Some(1) => ManipulationTask::LivingRoomRemote,
        Some(2) => ManipulationTask::BackyardTrowel,
        _ => {
            return Err(ArmPackedRefusal::new(ArmPackedRefusalCode::InvalidConfig));
        }
    };
    // Zero selects the owner preset; any other value must be finite and is
    // range-checked by the owner's own validate_config.
    let optional_override = |value: f64| -> Result<Option<f64>, ArmPackedRefusal> {
        if !value.is_finite() || value < 0.0 {
            return Err(ArmPackedRefusal::new(ArmPackedRefusalCode::InvalidConfig));
        }
        Ok(if value == 0.0 { None } else { Some(value) })
    };
    let object_mass_kg = optional_override(packet[8])?;
    let static_mu = optional_override(packet[9])?;
    let kinetic_mu = optional_override(packet[10])?;
    let mut obstacles = Vec::with_capacity(obstacle_count);
    for index in 0..obstacle_count {
        let base = ARM_CONFIG_FIXED_WORDS + index * ARM_OBSTACLE_WORDS;
        let words = &packet[base..base + ARM_OBSTACLE_WORDS];
        if words.iter().any(|value| !value.is_finite()) {
            return Err(ArmPackedRefusal::new(ArmPackedRefusalCode::InvalidConfig));
        }
        let role = match exact_u32(words[7]) {
            Some(0) => fs_scene::BodyRole::KeepOut,
            Some(1) => fs_scene::BodyRole::Support,
            _ => return Err(ArmPackedRefusal::new(ArmPackedRefusalCode::InvalidConfig)),
        };
        obstacles.push(ObstacleBox {
            center_m: fs_ga::Vec3::new(words[0], words[1], words[2]),
            half_extents_m: fs_ga::Vec3::new(words[3], words[4], words[5]),
            yaw_rad: words[6],
            role,
        });
    }
    let config = ManipulationConfig {
        task,
        step_s: packet[4],
        duration_s: packet[5],
        trace_stride,
        object_mass_kg,
        static_mu,
        kinetic_mu,
        obstacles,
    };
    if !(config.step_s.is_finite()
        && (1.0 / 240.0..=1.0 / 45.0).contains(&config.step_s)
        && config.duration_s.is_finite()
        && (3.0..=6.0).contains(&config.duration_s)
        && config.duration_s >= config.step_s)
    {
        return Err(ArmPackedRefusal::new(ArmPackedRefusalCode::InvalidConfig));
    }
    Ok(config)
}

fn arm_owner_refusal(error: &ManipulationError) -> ArmPackedRefusal {
    let code = match error {
        ManipulationError::InvalidConfig { .. } => ArmPackedRefusalCode::InvalidConfig,
        ManipulationError::ParameterCount { .. } => ArmPackedRefusalCode::ParameterCount,
        ManipulationError::NonFiniteParameter { .. } => ArmPackedRefusalCode::NonFiniteParameter,
        ManipulationError::Robot(_) => ArmPackedRefusalCode::Robot,
        ManipulationError::Geometry(_) => ArmPackedRefusalCode::Geometry,
        ManipulationError::Contact(_) => ArmPackedRefusalCode::Contact,
        ManipulationError::Friction(_) => ArmPackedRefusalCode::Friction,
        ManipulationError::Query(_) => ArmPackedRefusalCode::Query,
        ManipulationError::UnexpectedContactReceipt => {
            ArmPackedRefusalCode::UnexpectedContactReceipt
        }
        ManipulationError::NonFiniteObjective => ArmPackedRefusalCode::NonFiniteObjective,
    };
    ArmPackedRefusal::new(code)
}

fn arm_success_header(kind: u32, total_words: usize) -> Vec<f64> {
    vec![
        f64::from(ARM_PACKET_MAGIC),
        f64::from(ARM_PACKET_SCHEMA_VERSION),
        f64::from(PACKET_STATUS_OK),
        f64::from(kind),
        total_words as f64,
    ]
}

fn arm_refusal_packet(kind: u32, refusal: ArmPackedRefusal) -> Vec<f64> {
    vec![
        f64::from(ARM_PACKET_MAGIC),
        f64::from(ARM_PACKET_SCHEMA_VERSION),
        f64::from(PACKET_STATUS_REFUSAL),
        f64::from(kind),
        ARM_REFUSAL_WORDS as f64,
        f64::from(refusal.code as u32),
        refusal.detail.map_or(f64::NAN, |value| value as f64),
    ]
}

fn arm_receipt_packet(kind: u32, receipt: &ManipulationReceipt, include_trace: bool) -> Vec<f64> {
    let trace_words = if include_trace {
        receipt
            .trace
            .len()
            .checked_mul(ARM_TRACE_SAMPLE_WORDS)
            .and_then(|words| words.checked_add(1))
    } else {
        Some(0)
    };
    let Some(total_words) = trace_words.and_then(|words| ARM_RECEIPT_WORDS.checked_add(words))
    else {
        return arm_refusal_packet(
            kind,
            ArmPackedRefusal::new(ArmPackedRefusalCode::ShapeOverflow),
        );
    };
    let mut packet = arm_success_header(kind, total_words);
    packet.extend_from_slice(&[
        receipt.objective,
        receipt.final_object_error_m,
        receipt.minimum_reach_error_m,
        receipt.maximum_lift_m,
        receipt.actuator_work_j,
        receipt.collision_risk_integral,
        receipt.minimum_certified_clearance_m,
        receipt.possible_collision_time_s,
        receipt.collision_query_iterations as f64,
        receipt.control_limit_integral,
        receipt.first_grasp_time_s,
        receipt.grasp_duration_s,
        receipt.peak_grip_force_n,
        f64::from(u8::from(receipt.ever_grasped)),
        f64::from(u8::from(receipt.released_after_transport)),
        f64::from(u8::from(receipt.placed)),
        receipt.completed_steps as f64,
    ]);
    if include_trace {
        packet.push(receipt.trace.len() as f64);
        for sample in &receipt.trace {
            append_arm_trace_sample(&mut packet, sample);
        }
    }
    debug_assert_eq!(packet.len(), total_words);
    packet
}

fn append_arm_trace_sample(packet: &mut Vec<f64>, sample: &ManipulationTraceSample) {
    packet.extend_from_slice(&[
        sample.time_s,
        sample.gripper_width_m,
        sample.grip_normal_force_n,
        f64::from(u8::from(sample.grasped)),
    ]);
    packet.extend_from_slice(&sample.object_pose);
    for pose in &sample.link_pose {
        packet.extend_from_slice(pose);
    }
}

/// Kernel identity probe used after module instantiation.
#[must_use]
pub const fn kernel_version() -> &'static str {
    KERNEL_VERSION
}

fn parse_config(packet: &[f64]) -> Result<CmaConfig, PackedRefusal> {
    if packet.len() < CONFIG_FIXED_WORDS {
        return Err(PackedRefusal::new(PackedRefusalCode::MalformedPacket));
    }
    validate_input_header(packet, PACKET_KIND_CONFIG)?;
    let family = match exact_u32(packet[4]) {
        Some(0) => CmaFamily::Full,
        Some(1) => CmaFamily::Separable,
        Some(2) => CmaFamily::LmCma,
        Some(3) => CmaFamily::LmMa,
        _ => return Err(PackedRefusal::new(PackedRefusalCode::FamilyUnknown)),
    };
    let dimension = exact_usize(packet[5])
        .filter(|&value| (1..=SCALABLE_DIMENSION_LIMIT).contains(&value))
        .ok_or(PackedRefusal::new(PackedRefusalCode::DimensionInvalid))?;
    if family == CmaFamily::Full && dimension > FULL_DIMENSION_LIMIT {
        return Err(PackedRefusal::new(PackedRefusalCode::FullDimensionLimit));
    }
    let expected_words = CONFIG_FIXED_WORDS
        .checked_add(dimension)
        .ok_or(PackedRefusal::new(PackedRefusalCode::MalformedPacket))?;
    if packet.len() != expected_words {
        return Err(PackedRefusal::new(PackedRefusalCode::MalformedPacket));
    }
    let population = optional_usize(packet[6], PackedRefusalCode::PopulationInvalid)?;
    if population.is_some_and(|value| value < 4) {
        return Err(PackedRefusal::new(PackedRefusalCode::PopulationInvalid));
    }
    let memory = optional_usize(packet[7], PackedRefusalCode::MemoryInvalid)?;
    let max_evaluations = exact_usize(packet[8])
        .filter(|&value| value > 0)
        .ok_or(PackedRefusal::new(PackedRefusalCode::BudgetInvalid))?;
    let seed_low =
        exact_u32(packet[9]).ok_or(PackedRefusal::new(PackedRefusalCode::SeedInvalid))?;
    let seed_high =
        exact_u32(packet[10]).ok_or(PackedRefusal::new(PackedRefusalCode::SeedInvalid))?;
    let sigma = packet[11];
    if !sigma.is_finite() || sigma <= 0.0 {
        return Err(PackedRefusal::new(PackedRefusalCode::SigmaInvalid));
    }
    let mean = packet[CONFIG_FIXED_WORDS..].to_vec();
    if mean.iter().any(|value| !value.is_finite()) {
        return Err(PackedRefusal::new(PackedRefusalCode::MeanInvalid));
    }
    Ok(CmaConfig {
        family,
        mean,
        sigma,
        max_evaluations,
        seed: u64::from(seed_low) | (u64::from(seed_high) << 32),
        population_size: population,
        memory,
    })
}

fn parse_tell(
    packet: &[f64],
    admission: CmaAdmission,
    expected_generation: u64,
) -> Result<(u64, &[f64]), PackedRefusal> {
    if packet.len() < TELL_FIXED_WORDS {
        return Err(PackedRefusal::new(PackedRefusalCode::MalformedPacket));
    }
    validate_input_header(packet, PACKET_KIND_TELL)?;
    let generation = exact_u64(packet[4]).ok_or(PackedRefusal::at_generation(
        PackedRefusalCode::GenerationMismatch,
        expected_generation,
    ))?;
    let population = exact_usize(packet[5]).ok_or(PackedRefusal::at_generation(
        PackedRefusalCode::ObjectiveCount,
        expected_generation,
    ))?;
    if population != admission.population_size
        || packet.len() != TELL_FIXED_WORDS + admission.population_size
    {
        return Err(PackedRefusal::at_generation(
            PackedRefusalCode::ObjectiveCount,
            expected_generation,
        ));
    }
    let objectives = &packet[TELL_FIXED_WORDS..];
    if objectives.iter().any(|value| !value.is_finite()) {
        return Err(PackedRefusal::at_generation(
            PackedRefusalCode::NonFiniteObjective,
            expected_generation,
        ));
    }
    Ok((generation, objectives))
}

fn validate_input_header(packet: &[f64], kind: u32) -> Result<(), PackedRefusal> {
    if exact_u32(packet[0]) != Some(PACKET_MAGIC)
        || exact_u32(packet[1]) != Some(PACKET_SCHEMA_VERSION)
    {
        return Err(PackedRefusal::new(PackedRefusalCode::SchemaMismatch));
    }
    if exact_u32(packet[2]) != Some(kind) || exact_usize(packet[3]) != Some(packet.len()) {
        return Err(PackedRefusal::new(PackedRefusalCode::MalformedPacket));
    }
    Ok(())
}

fn optional_usize(value: f64, code: PackedRefusalCode) -> Result<Option<usize>, PackedRefusal> {
    let parsed = exact_usize(value).ok_or(PackedRefusal::new(code))?;
    Ok((parsed != 0).then_some(parsed))
}

fn exact_u32(value: f64) -> Option<u32> {
    (value.is_finite() && (0.0..=f64::from(u32::MAX)).contains(&value) && value.fract() == 0.0)
        .then_some(value as u32)
}

fn exact_u64(value: f64) -> Option<u64> {
    (value.is_finite() && (0.0..=MAX_SAFE_INTEGER).contains(&value) && value.fract() == 0.0)
        .then_some(value as u64)
}

fn exact_usize(value: f64) -> Option<usize> {
    exact_u32(value).and_then(|parsed| usize::try_from(parsed).ok())
}

fn owner_refusal(error: &CmaFamilyError) -> PackedRefusal {
    let code = match error {
        CmaFamilyError::EmptyMean => PackedRefusalCode::DimensionInvalid,
        CmaFamilyError::NonFiniteMean { .. } => PackedRefusalCode::MeanInvalid,
        CmaFamilyError::InvalidSigma { .. } => PackedRefusalCode::SigmaInvalid,
        CmaFamilyError::InvalidPopulation { .. } => PackedRefusalCode::PopulationInvalid,
        CmaFamilyError::BudgetTooSmall { .. } => PackedRefusalCode::BudgetInvalid,
        CmaFamilyError::InvalidMemory { .. } | CmaFamilyError::MemoryNotApplicable { .. } => {
            PackedRefusalCode::MemoryInvalid
        }
        CmaFamilyError::ShapeOverflow { .. } => PackedRefusalCode::ShapeOverflow,
        CmaFamilyError::RandomCounterOverflow => PackedRefusalCode::RandomCounterOverflow,
        CmaFamilyError::DenseEigensolver(_) => PackedRefusalCode::DenseWorkRefused,
        CmaFamilyError::AskAlreadyPending { .. } => PackedRefusalCode::AskAlreadyPending,
        CmaFamilyError::BudgetExhausted { .. } => PackedRefusalCode::BudgetExhausted,
        CmaFamilyError::NoPendingAsk => PackedRefusalCode::NoPendingAsk,
        CmaFamilyError::GenerationMismatch { .. } => PackedRefusalCode::GenerationMismatch,
        CmaFamilyError::BatchMismatch => PackedRefusalCode::BatchMismatch,
        CmaFamilyError::ObjectiveCount { .. } => PackedRefusalCode::ObjectiveCount,
        CmaFamilyError::NonFiniteObjective { .. } => PackedRefusalCode::NonFiniteObjective,
        _ => PackedRefusalCode::NumericalFailure,
    };
    let generation = match error {
        CmaFamilyError::AskAlreadyPending { generation }
        | CmaFamilyError::GenerationMismatch {
            expected: generation,
            ..
        } => Some(*generation),
        _ => None,
    };
    PackedRefusal { code, generation }
}

fn success_header(kind: u32, total_words: usize) -> Vec<f64> {
    vec![
        f64::from(PACKET_MAGIC),
        f64::from(PACKET_SCHEMA_VERSION),
        f64::from(PACKET_STATUS_OK),
        f64::from(kind),
        total_words as f64,
    ]
}

fn refusal_packet(kind: u32, refusal: PackedRefusal) -> Vec<f64> {
    vec![
        f64::from(PACKET_MAGIC),
        f64::from(PACKET_SCHEMA_VERSION),
        f64::from(PACKET_STATUS_REFUSAL),
        f64::from(kind),
        REFUSAL_WORDS as f64,
        f64::from(refusal.code as u32),
        refusal.generation.map_or(f64::NAN, |value| value as f64),
    ]
}

fn ask_packet(admission: CmaAdmission, evaluations: usize, batch: &CmaAsk) -> Vec<f64> {
    let candidate_words = admission.dimension * admission.population_size;
    let total_words = ASK_FIXED_WORDS + candidate_words;
    let mut packet = success_header(PACKET_KIND_ASK, total_words);
    packet.extend_from_slice(&[
        batch.generation() as f64,
        evaluations as f64,
        admission.dimension as f64,
        admission.population_size as f64,
    ]);
    for candidate in batch.candidates() {
        packet.extend_from_slice(candidate);
    }
    debug_assert_eq!(packet.len(), total_words);
    packet
}

fn snapshot_packet(kind: u32, admission: CmaAdmission, snapshot: &CmaSnapshot) -> Vec<f64> {
    let (shape_kind, shape_payload) = shape_payload(&snapshot.shape);
    let total_words = SNAPSHOT_FIXED_WORDS + 2 * admission.dimension + shape_payload.len();
    let mut packet = success_header(kind, total_words);
    let (normal_low, normal_high) = split_u64(admission.normal_stream_blocks);
    let best = snapshot.best.as_ref();
    packet.extend_from_slice(&[
        f64::from(family_id(snapshot.family)),
        admission.dimension as f64,
        snapshot.generation as f64,
        snapshot.evaluations as f64,
        snapshot.sigma,
        admission.population_size as f64,
        admission.parent_count as f64,
        admission.max_generations as f64,
        admission.admitted_evaluations as f64,
        f64::from(admission.stream_semantics_version),
        f64::from(admission.stream_kernel),
        f64::from(normal_low),
        f64::from(normal_high),
        f64::from(complexity_order_id(
            snapshot.complexity.sampling_per_candidate,
        )),
        f64::from(complexity_order_id(
            snapshot.complexity.update_per_generation,
        )),
        snapshot.complexity.persistent_scalars as f64,
        snapshot.complexity.pending_generation_scalars as f64,
        snapshot.complexity.update_workspace_scalars as f64,
        snapshot.complexity.dense_matrix_entries as f64,
        snapshot.complexity.memory_capacity as f64,
        f64::from(u8::from(best.is_some())),
        best.map_or(f64::NAN, |value| value.objective),
        best.map_or(f64::NAN, |value| value.generation as f64),
        best.map_or(f64::NAN, |value| value.candidate as f64),
        f64::from(shape_kind),
        shape_payload.len() as f64,
    ]);
    packet.extend_from_slice(&snapshot.mean);
    if let Some(best) = best {
        packet.extend_from_slice(&best.point);
    } else {
        packet.resize(packet.len() + admission.dimension, f64::NAN);
    }
    packet.extend_from_slice(&shape_payload);
    debug_assert_eq!(packet.len(), total_words);
    packet
}

fn shape_payload(shape: &CmaShapeSnapshot) -> (u32, Vec<f64>) {
    match shape {
        CmaShapeSnapshot::Full {
            diagonal,
            min_eigenvalue,
            max_eigenvalue,
            negative_weight_count,
        } => {
            let mut payload = Vec::with_capacity(3 + diagonal.len());
            payload.extend_from_slice(&[
                *negative_weight_count as f64,
                *min_eigenvalue,
                *max_eigenvalue,
            ]);
            payload.extend_from_slice(diagonal);
            (0, payload)
        }
        CmaShapeSnapshot::Diagonal {
            variances,
            negative_weight_count,
        } => {
            let mut payload = Vec::with_capacity(1 + variances.len());
            payload.push(*negative_weight_count as f64);
            payload.extend_from_slice(variances);
            (1, payload)
        }
        CmaShapeSnapshot::LimitedMemory {
            vectors,
            capacity,
            direction_norms,
        } => {
            let mut payload = Vec::with_capacity(2 + direction_norms.len());
            payload.extend_from_slice(&[*vectors as f64, *capacity as f64]);
            payload.extend_from_slice(direction_norms);
            (2, payload)
        }
    }
}

const fn family_id(family: CmaFamily) -> u32 {
    match family {
        CmaFamily::Full => 0,
        CmaFamily::Separable => 1,
        CmaFamily::LmCma => 2,
        CmaFamily::LmMa => 3,
    }
}

const fn complexity_order_id(order: CmaComplexityOrder) -> u32 {
    match order {
        CmaComplexityOrder::Linear => 0,
        CmaComplexityOrder::MemoryLinear => 1,
        CmaComplexityOrder::MemoryQuadratic => 4,
        CmaComplexityOrder::Quadratic => 2,
        CmaComplexityOrder::Cubic => 3,
    }
}

const fn split_u64(value: u64) -> (u32, u32) {
    (value as u32, (value >> 32) as u32)
}

#[cfg(target_arch = "wasm32")]
mod schema_two_wasm {
    use super::{PackedCmaSession, PackedG1WalkingEvaluator, PackedManipulationEvaluator};
    use wasm_bindgen::prelude::wasm_bindgen;

    /// Stateful schema-2 browser session. Construction never throws; inspect
    /// `receipt()` for admission or a typed refusal packet.
    #[wasm_bindgen]
    pub struct CmaesVizSession {
        inner: PackedCmaSession,
    }

    #[wasm_bindgen]
    impl CmaesVizSession {
        /// Create a session from one packed configuration.
        #[wasm_bindgen(constructor)]
        #[must_use]
        pub fn new(config: &[f64]) -> Self {
            Self {
                inner: PackedCmaSession::new(config),
            }
        }

        /// Return admission and the current compact snapshot.
        #[must_use]
        pub fn receipt(&self) -> Vec<f64> {
            self.inner.receipt_packet()
        }

        /// Return one complete row-major candidate population.
        #[must_use]
        pub fn ask(&mut self) -> Vec<f64> {
            self.inner.ask_packet()
        }

        /// Tell one packed objective payload and return the updated snapshot.
        #[must_use]
        pub fn tell(&mut self, objectives: &[f64]) -> Vec<f64> {
            self.inner.tell_packet(objectives)
        }
    }

    /// Stateful browser evaluator for the owner-composed G1 walking problem.
    /// Construction never throws; inspect `receipt()` for admission or a typed
    /// refusal packet.
    #[wasm_bindgen]
    pub struct G1WalkingVizEvaluator {
        inner: PackedG1WalkingEvaluator,
    }

    #[wasm_bindgen]
    impl G1WalkingVizEvaluator {
        /// Create an evaluator from one packed experiment configuration.
        #[wasm_bindgen(constructor)]
        #[must_use]
        pub fn new(config: &[f64]) -> Self {
            Self {
                inner: PackedG1WalkingEvaluator::new(config),
            }
        }

        /// Return admitted controls and exact render-layout dimensions.
        #[must_use]
        pub fn receipt(&self) -> Vec<f64> {
            self.inner.receipt_packet()
        }

        /// Return the disclosed sparse 5,040-D stabilizing curriculum mean.
        #[must_use]
        pub fn stabilizing_policy_mean(&self) -> Vec<f64> {
            super::g1_walking::g1_stabilizing_policy_mean().to_vec()
        }

        /// Return the disclosed sparse 5,040-D walking curriculum mean.
        #[must_use]
        pub fn walking_curriculum_mean(&self) -> Vec<f64> {
            super::g1_walking::g1_walking_curriculum_mean().to_vec()
        }

        /// Evaluate one 5,040-D policy without retaining link poses.
        #[must_use]
        pub fn evaluate(&self, parameters: &[f64]) -> Vec<f64> {
            self.inner.evaluate_packet(parameters)
        }

        /// Evaluate a flat complete population in one boundary call.
        #[must_use]
        pub fn evaluate_population(&self, parameters: &[f64]) -> Vec<f64> {
            self.inner.evaluate_population_packet(parameters)
        }

        /// Evaluate one policy and return decimated owner-derived link poses.
        #[must_use]
        pub fn trace(&self, parameters: &[f64]) -> Vec<f64> {
            self.inner.trace_packet(parameters)
        }
    }

    /// Stateful browser evaluator for the owner-composed household-arm problem.
    /// Construction never throws; inspect `receipt()` for admission or a typed
    /// refusal packet.
    #[wasm_bindgen]
    pub struct HouseholdManipulationVizEvaluator {
        inner: PackedManipulationEvaluator,
    }

    #[wasm_bindgen]
    impl HouseholdManipulationVizEvaluator {
        /// Create an evaluator from one packed experiment configuration.
        #[wasm_bindgen(constructor)]
        #[must_use]
        pub fn new(config: &[f64]) -> Self {
            Self {
                inner: PackedManipulationEvaluator::new(config),
            }
        }

        /// Return admitted controls, scene data, and render-layout dimensions.
        #[must_use]
        pub fn receipt(&self) -> Vec<f64> {
            self.inner.receipt_packet()
        }

        /// Return the disclosed source-feasible 128-D curriculum mean.
        #[must_use]
        pub fn curriculum_policy_mean(&self) -> Vec<f64> {
            self.inner.curriculum_policy_mean()
        }

        /// Evaluate one 128-D policy without retaining object/link poses.
        #[must_use]
        pub fn evaluate(&self, parameters: &[f64]) -> Vec<f64> {
            self.inner.evaluate_packet(parameters)
        }

        /// Evaluate a flat complete population in one boundary call.
        #[must_use]
        pub fn evaluate_population(&self, parameters: &[f64]) -> Vec<f64> {
            self.inner.evaluate_population_packet(parameters)
        }

        /// Evaluate one policy and return decimated owner-derived poses.
        #[must_use]
        pub fn trace(&self, parameters: &[f64]) -> Vec<f64> {
            self.inner.trace_packet(parameters)
        }
    }

    /// Kernel identity probe after module instantiation.
    #[wasm_bindgen]
    #[must_use]
    pub fn cmaes_viz_kernel_version() -> String {
        super::KERNEL_VERSION.to_string()
    }
}

#[cfg(test)]
#[allow(clippy::float_cmp)] // packet selector/count words require exact equality
mod schema_two_tests {
    use super::*;

    fn config_packet(
        family: CmaFamily,
        mean: &[f64],
        sigma: f64,
        population: usize,
        memory: usize,
        budget: usize,
        seed: u64,
    ) -> Vec<f64> {
        let total_words = CONFIG_FIXED_WORDS + mean.len();
        let (seed_low, seed_high) = split_u64(seed);
        let mut packet = vec![
            f64::from(PACKET_MAGIC),
            f64::from(PACKET_SCHEMA_VERSION),
            f64::from(PACKET_KIND_CONFIG),
            total_words as f64,
            f64::from(family_id(family)),
            mean.len() as f64,
            population as f64,
            memory as f64,
            budget as f64,
            f64::from(seed_low),
            f64::from(seed_high),
            sigma,
        ];
        packet.extend_from_slice(mean);
        packet
    }

    fn tell_packet(generation: u64, objectives: &[f64]) -> Vec<f64> {
        let total_words = TELL_FIXED_WORDS + objectives.len();
        let mut packet = vec![
            f64::from(PACKET_MAGIC),
            f64::from(PACKET_SCHEMA_VERSION),
            f64::from(PACKET_KIND_TELL),
            total_words as f64,
            generation as f64,
            objectives.len() as f64,
        ];
        packet.extend_from_slice(objectives);
        packet
    }

    fn g1_config_packet(duration_s: f64, trace_stride: usize) -> Vec<f64> {
        g1_config_packet_with(duration_s, trace_stride, &[])
    }

    /// Schema-8 walking config with an explicit keep-out roster. An empty
    /// roster is the schema-7 equivalent.
    fn g1_config_packet_with(
        duration_s: f64,
        trace_stride: usize,
        obstacles: &[crate::g1_walking::ObstacleBox],
    ) -> Vec<f64> {
        let mut packet = vec![
            f64::from(G1_PACKET_MAGIC),
            f64::from(G1_PACKET_SCHEMA_VERSION),
            f64::from(G1_PACKET_KIND_CONFIG),
            (G1_CONFIG_FIXED_WORDS + obstacles.len() * G1_OBSTACLE_WORDS) as f64,
            1.0 / 120.0,
            duration_s,
            0.65,
            1.55,
            trace_stride as f64,
            f64::from(G1Task::Walking as u32),
            f64::from(G1Challenge::Flat as u32),
            obstacles.len() as f64,
        ];
        for obstacle in obstacles {
            packet.extend_from_slice(&[
                obstacle.center_m[0],
                obstacle.center_m[1],
                obstacle.center_m[2],
                obstacle.half_extents_m[0],
                obstacle.half_extents_m[1],
                obstacle.half_extents_m[2],
                obstacle.yaw_rad,
                match obstacle.role {
                    fs_scene::BodyRole::KeepOut => 0.0,
                    fs_scene::BodyRole::Support => 1.0,
                },
            ]);
        }
        packet
    }

    fn arm_config_packet(task: ManipulationTask, duration_s: f64, trace_stride: usize) -> Vec<f64> {
        arm_config_packet_with(task, duration_s, trace_stride, 0.0, 0.0, 0.0, &[])
    }

    /// Schema-3 config packet with explicit overrides and keep-out boxes.
    /// Zero overrides plus an empty roster is the schema-2 equivalent.
    fn arm_config_packet_with(
        task: ManipulationTask,
        duration_s: f64,
        trace_stride: usize,
        object_mass_kg: f64,
        static_mu: f64,
        kinetic_mu: f64,
        obstacles: &[ObstacleBox],
    ) -> Vec<f64> {
        let mut packet = vec![
            f64::from(ARM_PACKET_MAGIC),
            f64::from(ARM_PACKET_SCHEMA_VERSION),
            f64::from(ARM_PACKET_KIND_CONFIG),
            (ARM_CONFIG_FIXED_WORDS + obstacles.len() * ARM_OBSTACLE_WORDS) as f64,
            1.0 / 90.0,
            duration_s,
            trace_stride as f64,
            f64::from(task as u32),
            object_mass_kg,
            static_mu,
            kinetic_mu,
            obstacles.len() as f64,
        ];
        for obstacle in obstacles {
            packet.extend_from_slice(&[
                obstacle.center_m.x,
                obstacle.center_m.y,
                obstacle.center_m.z,
                obstacle.half_extents_m.x,
                obstacle.half_extents_m.y,
                obstacle.half_extents_m.z,
                obstacle.yaw_rad,
                match obstacle.role {
                    fs_scene::BodyRole::KeepOut => 0.0,
                    fs_scene::BodyRole::Support => 1.0,
                },
            ]);
        }
        packet
    }

    fn assert_g1_success(packet: &[f64], kind: u32) {
        assert_eq!(packet[0], f64::from(G1_PACKET_MAGIC));
        assert_eq!(packet[1], f64::from(G1_PACKET_SCHEMA_VERSION));
        assert_eq!(packet[2], f64::from(PACKET_STATUS_OK));
        assert_eq!(packet[3], f64::from(kind));
        assert_eq!(packet[4], packet.len() as f64);
    }

    fn assert_arm_packet_ok(packet: &[f64], kind: u32) {
        assert_eq!(packet[0], f64::from(ARM_PACKET_MAGIC));
        assert_eq!(packet[1], f64::from(ARM_PACKET_SCHEMA_VERSION));
        assert_eq!(packet[2], f64::from(PACKET_STATUS_OK));
        assert_eq!(packet[3], f64::from(kind));
        assert_eq!(packet[4], packet.len() as f64);
    }

    fn sphere_objectives(ask: &[f64]) -> Vec<f64> {
        let dimension = ask[7] as usize;
        let population = ask[8] as usize;
        ask[ASK_FIXED_WORDS..]
            .chunks_exact(dimension)
            .take(population)
            .map(|point| point.iter().map(|value| value * value).sum())
            .collect()
    }

    fn assert_success(packet: &[f64], kind: u32) {
        assert_eq!(packet[0], f64::from(PACKET_MAGIC));
        assert_eq!(packet[1], f64::from(PACKET_SCHEMA_VERSION));
        assert_eq!(packet[2], f64::from(PACKET_STATUS_OK));
        assert_eq!(packet[3], f64::from(kind));
        assert_eq!(packet[4], packet.len() as f64);
    }

    fn assert_word_identical(left: &[f64], right: &[f64]) {
        assert_eq!(left.len(), right.len());
        for (index, (&left_word, &right_word)) in left.iter().zip(right).enumerate() {
            assert_eq!(
                left_word.to_bits(),
                right_word.to_bits(),
                "packed word {index} differs"
            );
        }
    }

    fn refusal_code(packet: &[f64]) -> u32 {
        assert_eq!(packet.len(), REFUSAL_WORDS);
        assert_eq!(packet[2], f64::from(PACKET_STATUS_REFUSAL));
        packet[5] as u32
    }

    #[test]
    fn all_families_dispatch_to_owner_shapes_and_complexities() {
        for family in [
            CmaFamily::Full,
            CmaFamily::Separable,
            CmaFamily::LmCma,
            CmaFamily::LmMa,
        ] {
            let memory = usize::from(matches!(family, CmaFamily::LmCma | CmaFamily::LmMa)) * 3;
            let config = config_packet(family, &[1.0; 8], 0.5, 8, memory, 16, 17);
            let mut session = PackedCmaSession::new(&config);
            let receipt = session.receipt_packet();
            assert_success(&receipt, PACKET_KIND_ADMISSION);
            assert_eq!(receipt[5], f64::from(family_id(family)));
            assert_eq!(receipt[23] > 0.0, family == CmaFamily::Full);
            assert_eq!(receipt[24], memory as f64);
            let (sampling_order, update_order) = match family {
                CmaFamily::Full => (2.0, 3.0),
                CmaFamily::Separable => (0.0, 0.0),
                CmaFamily::LmCma => (1.0, 4.0),
                CmaFamily::LmMa => (1.0, 1.0),
            };
            assert_eq!(receipt[18], sampling_order);
            assert_eq!(receipt[19], update_order);

            let ask = session.ask_packet();
            assert_success(&ask, PACKET_KIND_ASK);
            let objectives = sphere_objectives(&ask);
            let snapshot = session.tell_packet(&tell_packet(0, &objectives));
            assert_success(&snapshot, PACKET_KIND_SNAPSHOT);
            let expected_shape = match family {
                CmaFamily::Full => 0.0,
                CmaFamily::Separable => 1.0,
                CmaFamily::LmCma | CmaFamily::LmMa => 2.0,
            };
            assert_eq!(snapshot[29], expected_shape);
            assert_eq!(snapshot[7], 1.0);
            assert_eq!(snapshot[8], 8.0);
        }
    }

    #[test]
    fn packed_boundary_matches_direct_owner_for_every_family() {
        for family in [
            CmaFamily::Full,
            CmaFamily::Separable,
            CmaFamily::LmCma,
            CmaFamily::LmMa,
        ] {
            let memory = usize::from(matches!(family, CmaFamily::LmCma | CmaFamily::LmMa)) * 4;
            let mean = vec![0.75; 10];
            let owner_config = CmaConfig {
                family,
                mean: mean.clone(),
                sigma: 0.4,
                max_evaluations: 18,
                seed: 0x0123_4567_89AB_CDEF,
                population_size: Some(9),
                memory: (memory != 0).then_some(memory),
            };
            let admission = admit_cma(&owner_config).expect("direct owner admission");
            let mut owner = CmaOptimizer::new(owner_config).expect("direct owner construction");
            let mut packed = PackedCmaSession::new(&config_packet(
                family,
                &mean,
                0.4,
                9,
                memory,
                18,
                0x0123_4567_89AB_CDEF,
            ));

            assert_word_identical(
                &packed.receipt_packet(),
                &snapshot_packet(PACKET_KIND_ADMISSION, admission, &owner.snapshot()),
            );
            let owner_ask = owner.ask().expect("direct owner ask");
            let packed_ask = packed.ask_packet();
            assert_word_identical(&packed_ask, &ask_packet(admission, 0, &owner_ask));
            let objectives = sphere_objectives(&packed_ask);
            let owner_snapshot = owner
                .tell(&owner_ask, &objectives)
                .expect("direct owner tell");
            assert_word_identical(
                &packed.tell_packet(&tell_packet(0, &objectives)),
                &snapshot_packet(PACKET_KIND_SNAPSHOT, admission, &owner_snapshot),
            );
        }
    }

    #[test]
    fn zero_selects_reference_population_and_dimension_based_memory_defaults() {
        let mean = vec![0.0; 100];
        let default_population = 17;

        let separable = PackedCmaSession::new(&config_packet(
            CmaFamily::Separable,
            &mean,
            0.5,
            0,
            0,
            default_population,
            3,
        ));
        let separable_receipt = separable.receipt_packet();
        assert_success(&separable_receipt, PACKET_KIND_ADMISSION);
        assert_eq!(separable_receipt[10], default_population as f64);

        for family in [CmaFamily::LmCma, CmaFamily::LmMa] {
            let limited_memory =
                PackedCmaSession::new(&config_packet(family, &mean, 0.5, 4, 0, 4, 3));
            let receipt = limited_memory.receipt_packet();
            assert_success(&receipt, PACKET_KIND_ADMISSION);
            assert_eq!(receipt[10], 4.0);
            assert_eq!(receipt[24], default_population as f64);
        }
    }

    #[test]
    fn seeded_sessions_replay_ask_and_tell_packets_bit_for_bit() {
        for family in [
            CmaFamily::Full,
            CmaFamily::Separable,
            CmaFamily::LmCma,
            CmaFamily::LmMa,
        ] {
            let memory = usize::from(matches!(family, CmaFamily::LmCma | CmaFamily::LmMa)) * 4;
            let config = config_packet(
                family,
                &[1.25; 12],
                0.75,
                10,
                memory,
                30,
                0xFEDC_BA98_7654_3210,
            );
            let mut left = PackedCmaSession::new(&config);
            let mut right = PackedCmaSession::new(&config);
            for generation in 0..3 {
                let left_ask = left.ask_packet();
                let right_ask = right.ask_packet();
                assert_eq!(left_ask, right_ask);
                let objectives = sphere_objectives(&left_ask);
                let tell = tell_packet(generation, &objectives);
                assert_eq!(left.tell_packet(&tell), right.tell_packet(&tell));
            }
        }
    }

    #[test]
    fn malformed_tell_is_retryable_and_budget_is_exact() {
        let config = config_packet(CmaFamily::Separable, &[1.0; 5], 0.5, 10, 0, 25, 9);
        let mut session = PackedCmaSession::new(&config);
        let receipt = session.receipt_packet();
        assert_eq!(receipt[13], 20.0);

        let ask = session.ask_packet();
        assert_eq!(
            refusal_code(&session.ask_packet()),
            PackedRefusalCode::AskAlreadyPending as u32
        );
        let objectives = sphere_objectives(&ask);
        let wrong_count = tell_packet(9, &objectives[..9]);
        let wrong_count_refusal = session.tell_packet(&wrong_count);
        assert_eq!(
            refusal_code(&wrong_count_refusal),
            PackedRefusalCode::ObjectiveCount as u32
        );
        assert_eq!(
            wrong_count_refusal[6], 0.0,
            "repair metadata must name the pending owner generation"
        );
        let mut non_integral_count = tell_packet(9, &objectives);
        non_integral_count[5] = 9.5;
        let non_integral_count_refusal = session.tell_packet(&non_integral_count);
        assert_eq!(
            refusal_code(&non_integral_count_refusal),
            PackedRefusalCode::ObjectiveCount as u32
        );
        assert_eq!(
            non_integral_count_refusal[6], 0.0,
            "malformed count metadata must name the pending owner generation"
        );
        let wrong_generation = tell_packet(1, &objectives);
        assert_eq!(
            refusal_code(&session.tell_packet(&wrong_generation)),
            PackedRefusalCode::GenerationMismatch as u32
        );
        assert_success(
            &session.tell_packet(&tell_packet(0, &objectives)),
            PACKET_KIND_SNAPSHOT,
        );

        let second = session.ask_packet();
        let second_objectives = sphere_objectives(&second);
        assert_success(
            &session.tell_packet(&tell_packet(1, &second_objectives)),
            PACKET_KIND_SNAPSHOT,
        );
        assert_eq!(
            refusal_code(&session.ask_packet()),
            PackedRefusalCode::BudgetExhausted as u32
        );
        assert_eq!(
            refusal_code(&session.tell_packet(&tell_packet(2, &second_objectives))),
            PackedRefusalCode::NoPendingAsk as u32
        );
    }

    fn assert_scalable_5040d_generation(family: CmaFamily) {
        const N: usize = 5_040;
        let memory = usize::from(matches!(family, CmaFamily::LmCma | CmaFamily::LmMa)) * 4;
        let config = config_packet(family, &vec![0.25; N], 0.1, 16, memory, 16, 29);
        let mut session = PackedCmaSession::new(&config);
        let receipt = session.receipt_packet();
        assert_success(&receipt, PACKET_KIND_ADMISSION);
        assert_eq!(receipt[6], N as f64);
        assert_eq!(receipt[23], 0.0);
        assert!(receipt[20] < 20.0 * N as f64);

        let ask = session.ask_packet();
        assert_success(&ask, PACKET_KIND_ASK);
        assert_eq!(ask.len(), ASK_FIXED_WORDS + 16 * N);
        let objectives = sphere_objectives(&ask);
        let snapshot = session.tell_packet(&tell_packet(0, &objectives));
        assert_success(&snapshot, PACKET_KIND_SNAPSHOT);
        assert_eq!(snapshot[8], 16.0);
    }

    #[test]
    fn separable_executes_a_real_5040d_generation_without_dense_state() {
        assert_scalable_5040d_generation(CmaFamily::Separable);
    }

    #[test]
    fn lm_cma_executes_a_real_5040d_generation_without_dense_state() {
        assert_scalable_5040d_generation(CmaFamily::LmCma);
    }

    #[test]
    fn lm_cma_packed_5040d_plateau_remains_finite_for_40_generations() {
        const N: usize = 5_040;
        const POPULATION: usize = 16;
        const GENERATIONS: usize = 40;
        let config = config_packet(
            CmaFamily::LmCma,
            &vec![0.0; N],
            0.01,
            POPULATION,
            12,
            POPULATION * GENERATIONS,
            0x4731_5040,
        );
        let mut session = PackedCmaSession::new(&config);
        assert_success(&session.receipt_packet(), PACKET_KIND_ADMISSION);
        for generation in 0..GENERATIONS {
            let ask = session.ask_packet();
            assert_success(&ask, PACKET_KIND_ASK);
            let maximum_coordinate = ask[ASK_FIXED_WORDS..]
                .iter()
                .map(|value| value.abs())
                .fold(0.0, f64::max);
            assert!(
                maximum_coordinate.is_finite() && maximum_coordinate < 1.0e6,
                "packed LM-CMA escaped its finite search scale: {maximum_coordinate:e}"
            );
            let snapshot = session.tell_packet(&tell_packet(generation as u64, &[1.0; POPULATION]));
            assert_success(&snapshot, PACKET_KIND_SNAPSHOT);
            let shape_start = SNAPSHOT_FIXED_WORDS + 2 * N;
            let stored_vectors = snapshot[shape_start] as usize;
            assert_eq!(snapshot[30] as usize, stored_vectors + 2);
            assert!(
                snapshot[shape_start + 2..shape_start + 2 + stored_vectors]
                    .iter()
                    .all(|value| value.is_finite())
            );
        }
    }

    #[test]
    fn lm_ma_executes_a_real_5040d_generation_without_dense_state() {
        assert_scalable_5040d_generation(CmaFamily::LmMa);
    }

    #[test]
    fn full_large_dimension_and_nonfinite_objectives_are_typed_refusals() {
        let config = config_packet(
            CmaFamily::Full,
            &vec![0.0; FULL_DIMENSION_LIMIT + 1],
            1.0,
            8,
            0,
            8,
            1,
        );
        let session = PackedCmaSession::new(&config);
        assert_eq!(
            refusal_code(&session.receipt_packet()),
            PackedRefusalCode::FullDimensionLimit as u32
        );

        let config = config_packet(CmaFamily::LmMa, &[0.0; 16], 1.0, 8, 4, 8, 1);
        let mut session = PackedCmaSession::new(&config);
        let _ask = session.ask_packet();
        let mut objectives = vec![0.0; 8];
        objectives[3] = f64::NAN;
        assert_eq!(
            refusal_code(&session.tell_packet(&tell_packet(0, &objectives))),
            PackedRefusalCode::NonFiniteObjective as u32
        );
        objectives[3] = 0.0;
        assert_success(
            &session.tell_packet(&tell_packet(0, &objectives)),
            PACKET_KIND_SNAPSHOT,
        );

        let mut nonportable_budget =
            config_packet(CmaFamily::Separable, &[0.0; 2], 1.0, 4, 0, 4, 1);
        nonportable_budget[8] = f64::from(u32::MAX) + 1.0;
        let session = PackedCmaSession::new(&nonportable_budget);
        assert_eq!(
            refusal_code(&session.receipt_packet()),
            PackedRefusalCode::BudgetInvalid as u32,
            "native and wasm32 admission must share the same integer domain"
        );

        let extreme_population =
            config_packet(CmaFamily::LmMa, &vec![0.0; 100_000], 1.0, 221, 1, 221, 1);
        let session = PackedCmaSession::new(&extreme_population);
        assert_eq!(
            refusal_code(&session.receipt_packet()),
            PackedRefusalCode::BrowserMemoryRefused as u32,
            "admission must include the packed Rust and JavaScript ask copies"
        );
    }

    #[test]
    fn g1_packets_expose_owner_receipts_and_link_pose_traces() {
        let evaluator = PackedG1WalkingEvaluator::new(&g1_config_packet(0.10, 3));
        let admission = evaluator.receipt_packet();
        assert_g1_success(&admission, G1_PACKET_KIND_ADMISSION);
        assert_eq!(admission[5], G1_POLICY_DIMENSION as f64);
        assert_eq!(admission[6], G1_LINK_COUNT as f64);
        assert_eq!(admission[7], G1_LINK_POSE_WORDS as f64);
        assert_eq!(admission[8], G1_TRACE_SAMPLE_WORDS as f64);

        let parameters = vec![0.0; G1_POLICY_DIMENSION];
        let evaluation = evaluator.evaluate_packet(&parameters);
        assert_g1_success(&evaluation, G1_PACKET_KIND_EVALUATION);
        assert!(evaluation[5].is_finite());
        assert!(evaluation[6].is_finite());
        assert!((0.0..=6.0).contains(&evaluation[14]));

        let trace = evaluator.trace_packet(&parameters);
        assert_g1_success(&trace, G1_PACKET_KIND_TRACE);
        let sample_count = trace[G1_RECEIPT_WORDS] as usize;
        assert!(sample_count >= 2);
        assert_eq!(
            trace.len(),
            G1_RECEIPT_WORDS + 1 + sample_count * G1_TRACE_SAMPLE_WORDS
        );
        assert_eq!(trace[G1_RECEIPT_WORDS + 2], 1.0);
        assert_eq!(trace[G1_RECEIPT_WORDS + 3], 1.0);
    }

    #[test]
    fn lm_ma_ask_g1_population_tell_completes_one_real_5040d_generation() {
        let mut optimizer = PackedCmaSession::new(&config_packet(
            CmaFamily::LmMa,
            &vec![0.0; G1_POLICY_DIMENSION],
            0.02,
            4,
            4,
            4,
            0xC0FF_EE11,
        ));
        assert_success(&optimizer.receipt_packet(), PACKET_KIND_ADMISSION);
        let ask = optimizer.ask_packet();
        assert_success(&ask, PACKET_KIND_ASK);

        let evaluator = PackedG1WalkingEvaluator::new(&g1_config_packet(0.10, 3));
        let objectives = evaluator.evaluate_population_packet(&ask[ASK_FIXED_WORDS..]);
        assert_g1_success(&objectives, G1_PACKET_KIND_POPULATION);
        assert_eq!(objectives[5], 4.0);
        assert!(
            objectives[6..]
                .iter()
                .all(|objective| objective.is_finite())
        );

        let snapshot = optimizer.tell_packet(&tell_packet(0, &objectives[6..]));
        assert_success(&snapshot, PACKET_KIND_SNAPSHOT);
        assert_eq!(snapshot[7], 1.0);
        assert_eq!(snapshot[8], 4.0);
    }

    #[test]
    fn g1_population_refusals_name_the_failed_candidate() {
        let evaluator = PackedG1WalkingEvaluator::new(&g1_config_packet(0.10, 3));
        let mut parameters = vec![0.0; 2 * G1_POLICY_DIMENSION];
        parameters[G1_POLICY_DIMENSION + 17] = f64::NAN;
        let refusal = evaluator.evaluate_population_packet(&parameters);
        assert_eq!(refusal[0], f64::from(G1_PACKET_MAGIC));
        assert_eq!(refusal[2], f64::from(PACKET_STATUS_REFUSAL));
        assert_eq!(
            refusal[5],
            f64::from(G1PackedRefusalCode::NonFiniteParameter as u32)
        );
        assert_eq!(refusal[6], 1.0);
    }

    #[test]
    fn arm_packets_expose_self_describing_scenes_and_honest_owner_outcomes() {
        for task in [
            ManipulationTask::KitchenMug,
            ManipulationTask::LivingRoomRemote,
            ManipulationTask::BackyardTrowel,
        ] {
            let evaluator = PackedManipulationEvaluator::new(&arm_config_packet(task, 6.0, 3));
            let admission = evaluator.receipt_packet();
            assert_arm_packet_ok(&admission, ARM_PACKET_KIND_ADMISSION);
            assert_eq!(admission.len(), ARM_ADMISSION_WORDS);
            assert_eq!(admission[5], ARM_POLICY_DIMENSION as f64);
            assert_eq!(admission[6], ARM_JOINT_COUNT as f64);
            assert_eq!(admission[7], ARM_POLICY_KNOTS as f64);
            assert_eq!(admission[8], ARM_LINK_COUNT as f64);
            assert_eq!(admission[9], ARM_LINK_POSE_WORDS as f64);
            assert_eq!(admission[10], ARM_TRACE_SAMPLE_WORDS as f64);
            assert_eq!(admission[14], f64::from(task as u32));
            assert!(admission[15] > 0.0);
            assert!(admission[16..=32].iter().all(|word| word.is_finite()));

            let mean = evaluator.curriculum_policy_mean();
            assert_eq!(mean.len(), ARM_POLICY_DIMENSION);
            let evaluation = evaluator.evaluate_packet(&mean);
            assert_arm_packet_ok(&evaluation, ARM_PACKET_KIND_EVALUATION);
            assert_eq!(evaluation.len(), ARM_RECEIPT_WORDS);
            assert!(evaluation[5].is_finite());
            assert!(evaluation[6] <= PLACEMENT_TOLERANCE_M);
            assert!(evaluation[8] >= LIFT_TARGET_M);
            assert_eq!(&evaluation[18..=19], &[1.0, 1.0]);
            match task {
                ManipulationTask::KitchenMug | ManipulationTask::LivingRoomRemote => {
                    assert_eq!(evaluation[10], 0.0);
                    assert!(evaluation[11] >= manipulation::PLACEMENT_CLEARANCE_M);
                    assert_eq!(evaluation[12], 0.0);
                    assert_eq!(evaluation[20], 1.0);
                }
                ManipulationTask::BackyardTrowel => {
                    assert!(evaluation[10] > 0.0);
                    assert!(evaluation[11] < manipulation::PLACEMENT_CLEARANCE_M);
                    assert!(evaluation[12] > 0.0);
                    assert_eq!(evaluation[20], 0.0);
                }
            }

            let trace = evaluator.trace_packet(&mean);
            assert_arm_packet_ok(&trace, ARM_PACKET_KIND_TRACE);
            let sample_count = trace[ARM_RECEIPT_WORDS] as usize;
            assert!(sample_count >= 50);
            assert_eq!(
                trace.len(),
                ARM_RECEIPT_WORDS + 1 + sample_count * ARM_TRACE_SAMPLE_WORDS
            );
            let first_sample = ARM_RECEIPT_WORDS + 1;
            assert!(trace[first_sample].is_finite());
            assert!(
                (MIN_GRIPPER_WIDTH_M..=OPEN_GRIPPER_WIDTH_M).contains(&trace[first_sample + 1])
            );
            assert!(
                trace[first_sample + 4..first_sample + ARM_TRACE_SAMPLE_WORDS]
                    .iter()
                    .all(|word| word.is_finite())
            );
        }
    }

    #[test]
    fn every_cma_family_completes_a_real_128d_arm_generation() {
        let evaluator = PackedManipulationEvaluator::new(&arm_config_packet(
            ManipulationTask::KitchenMug,
            4.0,
            6,
        ));
        let mean = evaluator.curriculum_policy_mean();
        for family in [
            CmaFamily::Full,
            CmaFamily::Separable,
            CmaFamily::LmCma,
            CmaFamily::LmMa,
        ] {
            let memory = usize::from(matches!(family, CmaFamily::LmCma | CmaFamily::LmMa)) * 4;
            let mut optimizer = PackedCmaSession::new(&config_packet(
                family,
                &mean,
                0.001,
                4,
                memory,
                4,
                0x4152_4d31,
            ));
            assert_success(&optimizer.receipt_packet(), PACKET_KIND_ADMISSION);
            let ask = optimizer.ask_packet();
            assert_success(&ask, PACKET_KIND_ASK);
            let objectives = evaluator.evaluate_population_packet(&ask[ASK_FIXED_WORDS..]);
            assert_arm_packet_ok(&objectives, ARM_PACKET_KIND_POPULATION);
            assert_eq!(objectives[5], 4.0);
            assert!(
                objectives[6..]
                    .iter()
                    .all(|objective| objective.is_finite())
            );
            let snapshot = optimizer.tell_packet(&tell_packet(0, &objectives[6..]));
            assert_success(&snapshot, PACKET_KIND_SNAPSHOT);
            assert_eq!(snapshot[7], 1.0);
            assert_eq!(snapshot[8], 4.0);
        }
    }

    #[test]
    fn arm_population_refusals_name_the_failed_candidate() {
        let evaluator = PackedManipulationEvaluator::new(&arm_config_packet(
            ManipulationTask::KitchenMug,
            4.0,
            3,
        ));
        let mut parameters = vec![0.0; 2 * ARM_POLICY_DIMENSION];
        parameters[ARM_POLICY_DIMENSION + 17] = f64::NAN;
        let refusal = evaluator.evaluate_population_packet(&parameters);
        assert_eq!(refusal[0], f64::from(ARM_PACKET_MAGIC));
        assert_eq!(refusal[2], f64::from(PACKET_STATUS_REFUSAL));
        assert_eq!(
            refusal[5],
            f64::from(ArmPackedRefusalCode::NonFiniteParameter as u32)
        );
        assert_eq!(refusal[6], 1.0);
    }

    /// The walking owner's body-vs-obstacle guard has always been
    /// implemented; schema 8 is the first packet that can reach it. An empty
    /// roster must leave the rollout exactly as schema 7 produced it, and a
    /// wall across the walking axis must terminate it as `BodyObstacle` with
    /// a measured penetration the receipt now reports.
    #[test]
    fn g1_declared_obstacles_stop_the_body_and_report_penetration() {
        let clear = PackedG1WalkingEvaluator::new(&g1_config_packet(1.5, 12));
        let mean = crate::g1_walking::g1_walking_curriculum_mean().to_vec();
        let clear_receipt = clear.evaluate_packet(&mean);
        assert_g1_success(&clear_receipt, G1_PACKET_KIND_EVALUATION);
        assert_eq!(clear_receipt.len(), G1_RECEIPT_WORDS);
        // With no declared geometry the guard never fires, whatever else
        // ends the rollout at this step size.
        assert_ne!(
            clear_receipt[27],
            f64::from(crate::g1_walking::G1TerminationReason::BodyObstacle as u32),
            "an empty roster cannot terminate on an obstacle"
        );
        assert_eq!(clear_receipt[28], 0.0, "no roster means no measured penetration");

        // A wall straddling the start pose on the walking axis: the body
        // spheres are inside it from the first step.
        let wall = crate::g1_walking::ObstacleBox {
            center_m: [0.1, 0.0, 0.6],
            half_extents_m: [0.15, 1.0, 0.6],
            yaw_rad: 0.0,
            role: fs_scene::BodyRole::KeepOut,
        };
        let blocked = PackedG1WalkingEvaluator::new(&g1_config_packet_with(1.5, 12, &[wall]));
        let blocked_receipt = blocked.evaluate_packet(&mean);
        assert_g1_success(&blocked_receipt, G1_PACKET_KIND_EVALUATION);
        assert_eq!(
            blocked_receipt[27],
            f64::from(crate::g1_walking::G1TerminationReason::BodyObstacle as u32),
            "a wall on the walking axis must terminate the rollout"
        );
        assert!(
            blocked_receipt[28] > 0.0,
            "the guard must report the penetration it measured"
        );
        assert!(
            blocked_receipt[26] < clear_receipt[26],
            "termination must cut the completed step count"
        );

        // The same wall behind the robot leaves the rollout untouched.
        let behind = crate::g1_walking::ObstacleBox {
            center_m: [-2.0, 0.0, 0.6],
            ..wall
        };
        let unaffected = PackedG1WalkingEvaluator::new(&g1_config_packet_with(1.5, 12, &[behind]));
        let unaffected_receipt = unaffected.evaluate_packet(&mean);
        assert_eq!(
            unaffected_receipt[5].to_bits(),
            clear_receipt[5].to_bits(),
            "a box off the path must not perturb the objective"
        );
    }

    /// Schema-8 refuses a malformed or oversized keep-out roster by name.
    #[test]
    fn g1_schema_eight_refuses_bad_rosters() {
        let refusal_code = |packet: Vec<f64>| -> f64 {
            let evaluator = PackedG1WalkingEvaluator::new(&packet);
            let receipt = evaluator.receipt_packet();
            assert_eq!(receipt[2], f64::from(PACKET_STATUS_REFUSAL));
            receipt[5]
        };

        let mut wrong_count = g1_config_packet(1.5, 12);
        wrong_count[3] = 77.0;
        assert_eq!(
            refusal_code(wrong_count),
            f64::from(G1PackedRefusalCode::MalformedPacket as u32)
        );

        let too_many: Vec<crate::g1_walking::ObstacleBox> = (0..=G1_MAX_OBSTACLES)
            .map(|index| crate::g1_walking::ObstacleBox {
                center_m: [index as f64 * 0.01, 0.0, 0.5],
                half_extents_m: [0.05, 0.05, 0.05],
                yaw_rad: 0.0,
                role: fs_scene::BodyRole::KeepOut,
            })
            .collect();
        assert_eq!(
            refusal_code(g1_config_packet_with(1.5, 12, &too_many)),
            f64::from(G1PackedRefusalCode::InvalidConfig as u32)
        );

        let degenerate = crate::g1_walking::ObstacleBox {
            center_m: [0.5, 0.0, 0.5],
            half_extents_m: [0.0, 0.05, 0.05],
            yaw_rad: 0.0,
            role: fs_scene::BodyRole::KeepOut,
        };
        assert_eq!(
            refusal_code(g1_config_packet_with(1.5, 12, &[degenerate])),
            f64::from(G1PackedRefusalCode::InvalidConfig as u32)
        );

        let non_finite = crate::g1_walking::ObstacleBox {
            center_m: [0.5, 0.0, f64::NAN],
            half_extents_m: [0.05, 0.05, 0.05],
            yaw_rad: 0.0,
            role: fs_scene::BodyRole::KeepOut,
        };
        assert_eq!(
            refusal_code(g1_config_packet_with(1.5, 12, &[non_finite])),
            f64::from(G1PackedRefusalCode::InvalidConfig as u32)
        );
    }

    /// Schema-3 words are additive: a packet whose three override words are
    /// zero and whose keep-out roster is empty must reproduce the schema-2
    /// rollout bit for bit, and must report the owner's declared interface.
    #[test]
    fn schema_three_defaults_reproduce_the_preset_owner_rollout() {
        for task in [
            ManipulationTask::KitchenMug,
            ManipulationTask::LivingRoomRemote,
            ManipulationTask::BackyardTrowel,
        ] {
            let packed =
                PackedManipulationEvaluator::new(&arm_config_packet(task, 6.0, 3));
            let admission = packed.receipt_packet();
            assert_arm_packet_ok(&admission, ARM_PACKET_KIND_ADMISSION);
            assert_eq!(admission.len(), ARM_ADMISSION_WORDS);
            assert_eq!(admission[37], 0.82, "default static mu for {task:?}");
            assert_eq!(admission[38], 0.68, "default kinetic mu for {task:?}");
            assert_eq!(admission[39], 0.0, "no declared boxes for {task:?}");

            let mean = packed.curriculum_policy_mean();
            let packed_receipt = packed.evaluate_packet(&mean);
            assert_arm_packet_ok(&packed_receipt, ARM_PACKET_KIND_EVALUATION);

            // The same rollout through the internal API with an untouched
            // default config: identical bits, so the packet layer adds nothing.
            let direct = ManipulationEvaluator::new(ManipulationConfig {
                task,
                ..ManipulationConfig::default()
            })
            .expect("preset config admits");
            let direct_receipt = direct.evaluate(&mean).expect("preset rollout");
            assert_eq!(
                packed_receipt[5].to_bits(),
                direct_receipt.objective.to_bits(),
                "objective drift for {task:?}"
            );
            assert_eq!(packed_receipt[10], direct_receipt.collision_risk_integral);
            assert_eq!(packed_receipt[17], direct_receipt.peak_grip_force_n);
            assert_eq!(packed_receipt[20], f64::from(u8::from(direct_receipt.placed)));
        }
    }

    /// A declared work surface may be worked ON but not THROUGH.
    ///
    /// This is the arm's half of the lesson the humanoid learned by shipping a
    /// robot buried in the floor: a surface the renderer draws is not a body
    /// until the owner is told about it. Declaring the counter as a keep-out
    /// volume is not the fix either -- the flange envelope is a 46 mm sphere
    /// standing in for a slim gripper, so closing on an object that rests on
    /// the counter unavoidably overlaps it, and a keep-out counter refuses
    /// every rollout that does the task. Support is the role that makes the
    /// surface real without making the task impossible.
    #[test]
    fn a_work_surface_may_be_worked_on_but_not_through() {
        let baseline = PackedManipulationEvaluator::new(&arm_config_packet(
            ManipulationTask::LivingRoomRemote,
            6.0,
            3,
        ));
        let admission = baseline.receipt_packet();
        let mean = baseline.curriculum_policy_mean();
        let clear = baseline.evaluate_packet(&mean);
        assert_arm_packet_ok(&clear, ARM_PACKET_KIND_EVALUATION);
        assert_eq!(clear[20], 1.0, "the undeclared baseline places the object");

        // The counter the browser actually ships: a 0.90 x 0.90 m top, 90 mm
        // thick, centred 0.80 m out along -x with its top face at the admitted
        // support height. Sized against this owner -- the original 1.40 x 1.65 m
        // slab, whose near edge sat 150 mm from the arm's own base axis, refused
        // all three tasks because the arm's link envelopes must sweep through
        // where it was drawn.
        let support_height_m = admission[30];
        let counter = |role, lift_m: f64| ObstacleBox {
            center_m: fs_ga::Vec3::new(-0.8, 0.0, support_height_m - 0.045 + lift_m),
            half_extents_m: fs_ga::Vec3::new(0.45, 0.45, 0.045),
            yaw_rad: 0.0,
            role,
        };
        let with_counter = |body| {
            PackedManipulationEvaluator::new(&arm_config_packet_with(
                ManipulationTask::LivingRoomRemote,
                6.0,
                3,
                0.0,
                0.0,
                0.0,
                &[body],
            ))
            .evaluate_packet(&mean)
        };

        // A hard body is a motion filter, not a penalty: the owner refuses any
        // step that drives a link further into one. So the visible consequence
        // of getting the role wrong is that the task stops working.
        let as_keep_out = with_counter(counter(fs_scene::BodyRole::KeepOut, 0.0));
        assert_arm_packet_ok(&as_keep_out, ARM_PACKET_KIND_EVALUATION);
        assert_eq!(
            as_keep_out[20], 0.0,
            "a keep-out surface blocks the grasp that has to reach onto it"
        );

        // Declared as the surface it is, the same geometry leaves the task
        // exactly as it was.
        let as_support = with_counter(counter(fs_scene::BodyRole::Support, 0.0));
        assert_arm_packet_ok(&as_support, ARM_PACKET_KIND_EVALUATION);
        assert_eq!(
            as_support[20], 1.0,
            "a support surface must leave a rollout that works on it placeable"
        );

        // Support is not a blanket exemption. Raise the same surface 150 mm and
        // doing the task would mean reaching far below its top face, which is
        // exactly what the skin refuses.
        let raised = with_counter(counter(fs_scene::BodyRole::Support, 0.15));
        assert_arm_packet_ok(&raised, ARM_PACKET_KIND_EVALUATION);
        assert_eq!(
            raised[20], 0.0,
            "sinking past the support skin must still be refused"
        );
    }

    /// A declared keep-out box straddling the transport path must be scored
    /// as a hard link constraint; the identical box moved out of the workspace
    /// must leave the rollout exactly as it was.
    #[test]
    fn declared_obstacles_block_the_transport_path_only_where_they_stand() {
        let baseline = PackedManipulationEvaluator::new(&arm_config_packet(
            ManipulationTask::KitchenMug,
            6.0,
            3,
        ));
        let admission = baseline.receipt_packet();
        let mean = baseline.curriculum_policy_mean();
        let clear = baseline.evaluate_packet(&mean);
        assert_arm_packet_ok(&clear, ARM_PACKET_KIND_EVALUATION);
        assert_eq!(clear[10], 0.0, "preset mug path is collision free");
        assert_eq!(clear[20], 1.0, "preset mug rollout places the object");

        // Midpoint of the source-derived grasp and place stations.
        let midpoint = fs_ga::Vec3::new(
            0.5 * (admission[24] + admission[27]),
            0.5 * (admission[25] + admission[28]),
            0.5 * (admission[26] + admission[29]),
        );
        let blocker = ObstacleBox {
            center_m: midpoint,
            half_extents_m: fs_ga::Vec3::new(0.09, 0.09, 0.09),
            yaw_rad: 0.3,
            role: fs_scene::BodyRole::KeepOut,
        };
        let blocked = PackedManipulationEvaluator::new(&arm_config_packet_with(
            ManipulationTask::KitchenMug,
            6.0,
            3,
            0.0,
            0.0,
            0.0,
            &[blocker],
        ));
        let blocked_admission = blocked.receipt_packet();
        assert_eq!(blocked_admission[39], 1.0, "admission echoes the roster size");
        let blocked_receipt = blocked.evaluate_packet(&mean);
        assert_arm_packet_ok(&blocked_receipt, ARM_PACKET_KIND_EVALUATION);
        assert!(
            blocked_receipt[10] > 0.0,
            "a box on the path must accrue collision risk"
        );
        assert_eq!(blocked_receipt[20], 0.0, "a blocked rollout cannot be placed");

        // The same box, 1.5 m away in +x, is outside every link envelope.
        let distant = ObstacleBox {
            center_m: fs_ga::Vec3::new(midpoint.x + 1.5, midpoint.y, midpoint.z),
            ..blocker
        };
        let unaffected = PackedManipulationEvaluator::new(&arm_config_packet_with(
            ManipulationTask::KitchenMug,
            6.0,
            3,
            0.0,
            0.0,
            0.0,
            &[distant],
        ));
        let unaffected_receipt = unaffected.evaluate_packet(&mean);
        assert_eq!(
            unaffected_receipt[5].to_bits(),
            clear[5].to_bits(),
            "a distant box must not perturb the objective"
        );
        assert_eq!(unaffected_receipt[20], 1.0);
    }

    /// Mass and friction overrides reach the physics and are echoed back.
    #[test]
    fn mass_and_friction_overrides_change_the_rollout_and_the_admission() {
        let baseline = PackedManipulationEvaluator::new(&arm_config_packet(
            ManipulationTask::KitchenMug,
            6.0,
            3,
        ));
        let mean = baseline.curriculum_policy_mean();
        let baseline_admission = baseline.receipt_packet();
        let baseline_receipt = baseline.evaluate_packet(&mean);
        let preset_mass = baseline_admission[19];

        let heavy = PackedManipulationEvaluator::new(&arm_config_packet_with(
            ManipulationTask::KitchenMug,
            6.0,
            3,
            2.0,
            0.0,
            0.0,
            &[],
        ));
        let heavy_admission = heavy.receipt_packet();
        assert_eq!(heavy_admission[19], 2.0, "admission reports the effective mass");
        assert!((preset_mass - 2.0).abs() > 1.0);
        let heavy_receipt = heavy.evaluate_packet(&mean);
        assert_ne!(
            heavy_receipt[5].to_bits(),
            baseline_receipt[5].to_bits(),
            "a 6x heavier object must change the objective"
        );

        let slippery = PackedManipulationEvaluator::new(&arm_config_packet_with(
            ManipulationTask::KitchenMug,
            6.0,
            3,
            0.0,
            0.3,
            0.2,
            &[],
        ));
        let slippery_admission = slippery.receipt_packet();
        assert_eq!(slippery_admission[37], 0.3);
        assert_eq!(slippery_admission[38], 0.2);
        let slippery_receipt = slippery.evaluate_packet(&mean);
        assert_ne!(
            slippery_receipt[5].to_bits(),
            baseline_receipt[5].to_bits(),
            "a low-friction interface must change the objective"
        );
    }

    /// Every schema-3 envelope violation is refused by name, not clamped.
    #[test]
    fn schema_three_refuses_malformed_and_out_of_range_declarations() {
        let refusal_code = |packet: Vec<f64>| -> f64 {
            let evaluator = PackedManipulationEvaluator::new(&packet);
            let admission = evaluator.receipt_packet();
            assert_eq!(admission[2], f64::from(PACKET_STATUS_REFUSAL));
            admission[5]
        };

        // wordCount that disagrees with the packet length.
        let mut wrong_count = arm_config_packet(ManipulationTask::KitchenMug, 6.0, 3);
        wrong_count[3] = 99.0;
        assert_eq!(
            refusal_code(wrong_count),
            f64::from(ArmPackedRefusalCode::MalformedPacket as u32)
        );

        // A roster one box past the declared cap.
        let too_many: Vec<ObstacleBox> = (0..=MAX_DECLARED_OBSTACLES)
            .map(|index| ObstacleBox {
                center_m: fs_ga::Vec3::new(index as f64 * 0.01, 0.0, 1.0),
                half_extents_m: fs_ga::Vec3::new(0.05, 0.05, 0.05),
                yaw_rad: 0.0,
                role: fs_scene::BodyRole::KeepOut,
            })
            .collect();
        assert_eq!(
            refusal_code(arm_config_packet_with(
                ManipulationTask::KitchenMug,
                6.0,
                3,
                0.0,
                0.0,
                0.0,
                &too_many
            )),
            f64::from(ArmPackedRefusalCode::InvalidConfig as u32)
        );

        // Kinetic friction above static friction is not a physical interface.
        assert_eq!(
            refusal_code(arm_config_packet_with(
                ManipulationTask::KitchenMug,
                6.0,
                3,
                0.0,
                0.4,
                0.9,
                &[]
            )),
            f64::from(ArmPackedRefusalCode::InvalidConfig as u32)
        );

        // Negative mass.
        assert_eq!(
            refusal_code(arm_config_packet_with(
                ManipulationTask::KitchenMug,
                6.0,
                3,
                -1.0,
                0.0,
                0.0,
                &[]
            )),
            f64::from(ArmPackedRefusalCode::InvalidConfig as u32)
        );

        // Non-finite obstacle extent.
        assert_eq!(
            refusal_code(arm_config_packet_with(
                ManipulationTask::KitchenMug,
                6.0,
                3,
                0.0,
                0.0,
                0.0,
                &[ObstacleBox {
                    center_m: fs_ga::Vec3::new(0.3, 0.0, 0.8),
                    half_extents_m: fs_ga::Vec3::new(f64::NAN, 0.05, 0.05),
                    yaw_rad: 0.0,
                    role: fs_scene::BodyRole::KeepOut,
                }]
            )),
            f64::from(ArmPackedRefusalCode::InvalidConfig as u32)
        );

        // A degenerate (zero-thickness) box is refused rather than admitted.
        assert_eq!(
            refusal_code(arm_config_packet_with(
                ManipulationTask::KitchenMug,
                6.0,
                3,
                0.0,
                0.0,
                0.0,
                &[ObstacleBox {
                    center_m: fs_ga::Vec3::new(0.3, 0.0, 0.8),
                    half_extents_m: fs_ga::Vec3::new(0.0, 0.05, 0.05),
                    yaw_rad: 0.0,
                    role: fs_scene::BodyRole::KeepOut,
                }]
            )),
            f64::from(ArmPackedRefusalCode::InvalidConfig as u32)
        );
    }
}
