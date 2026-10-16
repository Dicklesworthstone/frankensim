//! G0/G3/G5 full-lifecycle replay for the public bounded hypervolume archive.
//!
//! A capacity-two, two-objective dyadic fixture exercises insertion, the
//! earliest-index policy under an exact three-way contribution tie, duplicate
//! and dominated no-ops, domination purges, a successful insertion whose
//! candidate is immediately evicted as the unique least contributor, and a
//! final two-member purge. Every candidate has a distinct decision payload.
//! The retained trajectory binds every public before/after member order and
//! objective/decision bit, public hypervolume bit, independently computed
//! denominator-256 area and exclusive contribution. Branch reason, purge, and
//! capacity-eviction records are explicitly ORACLE EXPECTATIONS reconstructed
//! from public before/candidate data; `HvArchive::insert` exposes no internal
//! branch receipt, so this target never presents those projections as observed
//! production facts. For each over-capacity projected pre-state, a separate
//! retained cross-check calls public production `hypervolume` for the full set
//! and every leave-one-out set and compares the exact binary64 results with the
//! integer denominator-256 oracle.
//!
//! The oracle below is deliberately integer-only: it enumerates the 16 by 16
//! reference grid and implements minimization dominance itself. It calls
//! neither production `hypervolume`/`dominates` nor `HvArchive`. A separate
//! disclosed `StreamKey` flips one bit in one retained after-checkpoint
//! objective. The stale payload, identity-consistently resealed payload,
//! retained-reference gate, first impossible transition, stable red evidence,
//! and local merge gate all refuse that corruption.
//!
//! This target does not claim arbitrary dimensions, capacities, or fronts;
//! malformed, non-finite, or outside-reference policy; Monte Carlo
//! contribution eviction; optimizer convergence or archive quality; sealed
//! public state; allocation, `Cx`, or cancellation behavior; cross-ISA
//! authority; authenticated persistence; performance; or authority over the
//! unobservable internal purge/rejection/eviction branch path.

#![deny(unsafe_code)]

use fs_blake3::{ContentHash, hash_domain};
use fs_dfo::{HvArchive, Individual, hypervolume};
use fs_obs::ident::{IdentityBuilder, ReplayIdentity};
use fs_obs::{Emitter, Event, EventKind, Severity};
use fs_rand::StreamKey;
use std::panic::catch_unwind;

const SUITE: &str = "fs-dfo/hv-archive-study-replay";
const CASE: &str = "capacity-two-dyadic-full-lifecycle";
const RED_CASE: &str = "seeded-retained-checkpoint-corruption";

const FIXTURE_IDENTITY_KIND: &str = "fs-dfo-hv-archive-fixture-v2";
const RESULT_IDENTITY_KIND: &str = "fs-dfo-hv-archive-result-v2";
const FIXTURE_DIGEST_DOMAIN: &str = "frankensim.fs-dfo.hv-archive-fixture.v2";
const RESULT_DIGEST_DOMAIN: &str = "frankensim.fs-dfo.hv-archive-result.v2";
const EVENT_DIGEST_DOMAIN: &str = "frankensim.fs-dfo.hv-archive-event.v2";
const SUPPLIED_FIXTURE_DIGEST_TRIPWIRE_DOMAIN: &str =
    "frankensim.fs-dfo.hv-archive-supplied-fixture-digest-tripwire.v2";

const CAPACITY: usize = 2;
const OBJECTIVES: usize = 2;
const DECISION_COORDINATES: usize = 2;
const GRID_DENOMINATOR: u16 = 16;
const AREA_DENOMINATOR: u16 = GRID_DENOMINATOR * GRID_DENOMINATOR;
const REFERENCE_GRID: [u16; OBJECTIVES] = [GRID_DENOMINATOR; OBJECTIVES];
const REFERENCE: [f64; OBJECTIVES] = [1.0, 1.0];
const CANDIDATE_COUNT: usize = 8;

const MUTATION_SEED: u64 = 0x4856_A2C4_FA11_0051;
const MUTATION_KERNEL: u32 = 0xA251;
const MUTATION_TILE: u32 = 0;
const MUTATION_BIT_BASE: u32 = 0;
const MUTATION_BIT_COUNT: u64 = 8;
const MUTABLE_STATE_CELLS: [(usize, usize, usize); 4] = [
    (2, 0, 0), // after C: B.f0
    (5, 1, 1), // after D: D.f1
    (6, 0, 1), // after E self-evicts: C.f1
    (7, 0, 0), // after F: F.f0
];

const NO_CLAIMS: &str = "arbitrary-dimensions-capacities-fronts;malformed-nonfinite-outside-reference-policy;mc-contribution-eviction;optimizer-convergence-archive-quality;sealed-public-state;allocation-Cx-cancellation;cross-ISA-authority;authenticated-persistence;performance;internal-purge-rejection-eviction-branch-authority";
const ORACLE_BRANCH_PROVENANCE: &str = "integer-oracle-expectation-not-production-observation";

const _: () = assert!(CAPACITY == 2);
const _: () = assert!(OBJECTIVES == 2);
const _: () = assert!(AREA_DENOMINATOR == 256);
const _: () = assert!(CANDIDATE_COUNT == 8);
const _: () = assert!(MUTATION_BIT_BASE + MUTATION_BIT_COUNT as u32 <= 52);

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct CandidateSpec {
    label: &'static str,
    x_bits: [u64; DECISION_COORDINATES],
    f_grid: [u16; OBJECTIVES],
}

impl CandidateSpec {
    fn individual(self) -> Individual {
        Individual {
            x: self.x_bits.into_iter().map(f64::from_bits).collect(),
            f: self.f_grid.into_iter().map(grid_value).collect(),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct IndividualBits {
    candidate_index: usize,
    x: Vec<u64>,
    f: Vec<u64>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct Checkpoint {
    members: Vec<IndividualBits>,
    hv_bits: u64,
    area_units: u16,
    exclusive_contribution_units: Vec<u16>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum OracleRejectionKind {
    DuplicateObjective,
    DominatedByMember,
}

impl OracleRejectionKind {
    const fn name(self) -> &'static str {
        match self {
            Self::DuplicateObjective => "duplicate-objective",
            Self::DominatedByMember => "dominated-by-member",
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct OracleRejectionExpectation {
    kind: OracleRejectionKind,
    blocker_index: usize,
    blocker: IndividualBits,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct OraclePurgedMemberExpectation {
    prior_index: usize,
    member: IndividualBits,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum OracleCapacityEvictionReason {
    LeastExclusiveHypervolumeContributor,
}

impl OracleCapacityEvictionReason {
    const fn name(self) -> &'static str {
        match self {
            Self::LeastExclusiveHypervolumeContributor => "least-exclusive-hypervolume-contributor",
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct OracleCapacityEvictionExpectation {
    reason: OracleCapacityEvictionReason,
    selected_index: usize,
    selected: IndividualBits,
    contribution_units: u16,
    tied_indices: Vec<usize>,
    candidate_evicted: bool,
}

/// Public production-hypervolume recomputation over one ORACLE-EXPECTED
/// over-capacity pre-state. This observes the public arithmetic API, not the
/// unexposed control-flow path taken inside `HvArchive::insert`.
#[derive(Debug, Clone, PartialEq, Eq)]
struct PublicHypervolumeCrossCheck {
    oracle_full_area_units: u16,
    oracle_leave_one_out_area_units: Vec<u16>,
    oracle_exclusive_contribution_area_units: Vec<u16>,
    full_hv_bits: u64,
    leave_one_out_hv_bits: Vec<u64>,
    exclusive_contribution_hv_bits: Vec<u64>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct TransitionRecord {
    ordinal: usize,
    candidate: IndividualBits,
    /// Directly observed public `HvArchive::insert` return value.
    insert_returned: bool,
    /// Directly observed public state and `HvArchive::hv` before insertion.
    before: Checkpoint,
    /// Integer-oracle EXPECTATION, not an observed intermediate production state.
    oracle_expected_pre_capacity: Checkpoint,
    /// Integer-oracle EXPECTATION, not an observed production purge receipt.
    oracle_expected_purged: Vec<OraclePurgedMemberExpectation>,
    /// Integer-oracle EXPECTATION, not an observed production rejection receipt.
    oracle_expected_rejection: Option<OracleRejectionExpectation>,
    /// Integer-oracle EXPECTATION, not an observed production eviction receipt.
    oracle_expected_capacity_eviction: Option<OracleCapacityEvictionExpectation>,
    /// Direct public `hypervolume` recomputation over the oracle-expected
    /// over-capacity set; not an observation of `insert`'s internal branch.
    public_hypervolume_cross_check: Option<PublicHypervolumeCrossCheck>,
    /// Directly observed public state and `HvArchive::hv` after insertion.
    after: Checkpoint,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct StudyRecord {
    initial: Checkpoint,
    transitions: Vec<TransitionRecord>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct StudyRun {
    fixture: ReplayIdentity,
    fixture_digest: ContentHash,
    record: StudyRecord,
    result: ReplayIdentity,
    result_digest: ContentHash,
}

#[derive(Debug, Clone, PartialEq, Eq)]
enum AdmissionError {
    SuppliedFixtureDigestMismatch {
        declared: [u8; 32],
        computed: [u8; 32],
    },
    RetainedFixtureIdentityMismatch {
        expected: [u8; 32],
        found: [u8; 32],
    },
    ResultPayloadIdentityMismatch {
        declared: [u8; 32],
        computed: [u8; 32],
    },
    ReferenceIdentityMismatch {
        expected: [u8; 32],
        found: [u8; 32],
    },
    SemanticInconsistency(String),
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct MutationCoordinate {
    transition_index: usize,
    member_index: usize,
    objective_index: usize,
    mantissa_bit: u32,
    selector_draws: u64,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct Mutation {
    seed: u64,
    kernel: u32,
    tile: u32,
    coordinate: MutationCoordinate,
    before: u64,
    after: u64,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct SeededCorruption {
    run: StudyRun,
    mutation: Mutation,
    stale_error: AdmissionError,
    reference_error: AdmissionError,
    semantic_error: AdmissionError,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct ExpectedRejection {
    kind: OracleRejectionKind,
    blocker_index: usize,
    blocker_candidate: usize,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct ExpectedCapacityEviction {
    selected_index: usize,
    selected_candidate: usize,
    contribution_units: u16,
    tied_indices: Vec<usize>,
    candidate_evicted: bool,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct BranchExpectation {
    candidate_index: usize,
    insert_returned: bool,
    rejection: Option<ExpectedRejection>,
    purged_prior_indices: Vec<usize>,
    purged_candidates: Vec<usize>,
    pre_capacity_candidates: Vec<usize>,
    capacity_eviction: Option<ExpectedCapacityEviction>,
    after_candidates: Vec<usize>,
    after_area_units: u16,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct OracleStep {
    insert_returned: bool,
    expected_pre_capacity: Checkpoint,
    expected_purged: Vec<OraclePurgedMemberExpectation>,
    expected_rejection: Option<OracleRejectionExpectation>,
    expected_capacity_eviction: Option<OracleCapacityEvictionExpectation>,
    after: Checkpoint,
}

fn usize_u64(value: usize) -> u64 {
    u64::try_from(value).expect("fixed HvArchive fixture cardinality fits u64")
}

fn digest_bytes(digest: ContentHash) -> [u8; 32] {
    *digest.as_bytes()
}

fn bit_hex_json(bits: &[u64]) -> String {
    format!(
        "[{}]",
        bits.iter()
            .map(|bits| format!("\"0x{bits:016x}\""))
            .collect::<Vec<_>>()
            .join(",")
    )
}

fn grid_value(units: u16) -> f64 {
    f64::from(units) / f64::from(GRID_DENOMINATOR)
}

fn area_value(units: u16) -> f64 {
    f64::from(units) / f64::from(AREA_DENOMINATOR)
}

fn candidate_specs() -> [CandidateSpec; CANDIDATE_COUNT] {
    [
        CandidateSpec {
            label: "A",
            x_bits: [0.125f64.to_bits(), (-0.125f64).to_bits()],
            f_grid: [4, 12],
        },
        CandidateSpec {
            label: "B",
            x_bits: [0.250f64.to_bits(), (-0.250f64).to_bits()],
            f_grid: [8, 8],
        },
        CandidateSpec {
            label: "C",
            x_bits: [0.375f64.to_bits(), (-0.375f64).to_bits()],
            f_grid: [12, 4],
        },
        CandidateSpec {
            label: "B-duplicate",
            x_bits: [0.500f64.to_bits(), (-0.500f64).to_bits()],
            f_grid: [8, 8],
        },
        CandidateSpec {
            label: "dominated-7/8",
            x_bits: [0.625f64.to_bits(), (-0.625f64).to_bits()],
            f_grid: [14, 14],
        },
        CandidateSpec {
            label: "D",
            x_bits: [0.750f64.to_bits(), (-0.750f64).to_bits()],
            f_grid: [6, 6],
        },
        CandidateSpec {
            label: "E",
            x_bits: [0.875f64.to_bits(), (-0.875f64).to_bits()],
            f_grid: [1, 15],
        },
        CandidateSpec {
            label: "F",
            x_bits: [1.000f64.to_bits(), (-1.000f64).to_bits()],
            f_grid: [2, 2],
        },
    ]
}

#[allow(clippy::too_many_lines)] // The literal table is the bound branch contract.
fn branch_expectations() -> Vec<BranchExpectation> {
    vec![
        BranchExpectation {
            candidate_index: 0,
            insert_returned: true,
            rejection: None,
            purged_prior_indices: vec![],
            purged_candidates: vec![],
            pre_capacity_candidates: vec![0],
            capacity_eviction: None,
            after_candidates: vec![0],
            after_area_units: 48,
        },
        BranchExpectation {
            candidate_index: 1,
            insert_returned: true,
            rejection: None,
            purged_prior_indices: vec![],
            purged_candidates: vec![],
            pre_capacity_candidates: vec![0, 1],
            capacity_eviction: None,
            after_candidates: vec![0, 1],
            after_area_units: 80,
        },
        BranchExpectation {
            candidate_index: 2,
            insert_returned: true,
            rejection: None,
            purged_prior_indices: vec![],
            purged_candidates: vec![],
            pre_capacity_candidates: vec![0, 1, 2],
            capacity_eviction: Some(ExpectedCapacityEviction {
                selected_index: 0,
                selected_candidate: 0,
                contribution_units: 16,
                tied_indices: vec![0, 1, 2],
                candidate_evicted: false,
            }),
            after_candidates: vec![1, 2],
            after_area_units: 80,
        },
        BranchExpectation {
            candidate_index: 3,
            insert_returned: false,
            rejection: Some(ExpectedRejection {
                kind: OracleRejectionKind::DuplicateObjective,
                blocker_index: 0,
                blocker_candidate: 1,
            }),
            purged_prior_indices: vec![],
            purged_candidates: vec![],
            pre_capacity_candidates: vec![1, 2],
            capacity_eviction: None,
            after_candidates: vec![1, 2],
            after_area_units: 80,
        },
        BranchExpectation {
            candidate_index: 4,
            insert_returned: false,
            rejection: Some(ExpectedRejection {
                kind: OracleRejectionKind::DominatedByMember,
                blocker_index: 0,
                blocker_candidate: 1,
            }),
            purged_prior_indices: vec![],
            purged_candidates: vec![],
            pre_capacity_candidates: vec![1, 2],
            capacity_eviction: None,
            after_candidates: vec![1, 2],
            after_area_units: 80,
        },
        BranchExpectation {
            candidate_index: 5,
            insert_returned: true,
            rejection: None,
            purged_prior_indices: vec![0],
            purged_candidates: vec![1],
            pre_capacity_candidates: vec![2, 5],
            capacity_eviction: None,
            after_candidates: vec![2, 5],
            after_area_units: 108,
        },
        BranchExpectation {
            candidate_index: 6,
            insert_returned: true,
            rejection: None,
            purged_prior_indices: vec![],
            purged_candidates: vec![],
            pre_capacity_candidates: vec![2, 5, 6],
            capacity_eviction: Some(ExpectedCapacityEviction {
                selected_index: 2,
                selected_candidate: 6,
                contribution_units: 5,
                tied_indices: vec![2],
                candidate_evicted: true,
            }),
            after_candidates: vec![2, 5],
            after_area_units: 108,
        },
        BranchExpectation {
            candidate_index: 7,
            insert_returned: true,
            rejection: None,
            purged_prior_indices: vec![0, 1],
            purged_candidates: vec![2, 5],
            pre_capacity_candidates: vec![7],
            capacity_eviction: None,
            after_candidates: vec![7],
            after_area_units: 196,
        },
    ]
}

fn member_from_spec(candidate_index: usize) -> IndividualBits {
    let spec = candidate_specs()[candidate_index];
    IndividualBits {
        candidate_index,
        x: spec.x_bits.to_vec(),
        f: spec
            .f_grid
            .into_iter()
            .map(grid_value)
            .map(f64::to_bits)
            .collect(),
    }
}

fn individual_bits(individual: &Individual) -> IndividualBits {
    let x: Vec<u64> = individual.x.iter().copied().map(f64::to_bits).collect();
    let candidate_index = candidate_specs()
        .iter()
        .position(|spec| spec.x_bits.as_slice() == x.as_slice())
        .expect("every public archive member comes from the bound candidate sequence");
    IndividualBits {
        candidate_index,
        x,
        f: individual.f.iter().copied().map(f64::to_bits).collect(),
    }
}

fn grid_point(member: &IndividualBits) -> Result<[u16; OBJECTIVES], String> {
    if member.f.len() != OBJECTIVES {
        return Err(format!(
            "candidate[{}]-objective-count:{}!={OBJECTIVES}",
            member.candidate_index,
            member.f.len()
        ));
    }
    let mut point = [0u16; OBJECTIVES];
    for (objective, &bits) in member.f.iter().enumerate() {
        let Some(units) = (0..=GRID_DENOMINATOR).find(|&units| grid_value(units).to_bits() == bits)
        else {
            return Err(format!(
                "candidate[{}].f[{objective}]=0x{bits:016x}-is-not-denominator-{GRID_DENOMINATOR}",
                member.candidate_index
            ));
        };
        point[objective] = units;
    }
    Ok(point)
}

/// Integer minimization dominance. This is independent of production
/// `fs_dfo::dominates`.
fn oracle_dominates(left: [u16; OBJECTIVES], right: [u16; OBJECTIVES]) -> bool {
    let no_worse = left
        .iter()
        .zip(right)
        .all(|(&left_coordinate, right_coordinate)| left_coordinate <= right_coordinate);
    let strictly_better = left
        .iter()
        .zip(right)
        .any(|(&left_coordinate, right_coordinate)| left_coordinate < right_coordinate);
    no_worse && strictly_better
}

/// Exact denominator-256 union area by unit-cell enumeration. This is
/// independent of production `fs_dfo::hypervolume`.
fn oracle_area(members: &[IndividualBits]) -> Result<u16, String> {
    let points: Vec<[u16; OBJECTIVES]> =
        members.iter().map(grid_point).collect::<Result<_, _>>()?;
    let mut area = 0u16;
    for x in 0..REFERENCE_GRID[0] {
        for y in 0..REFERENCE_GRID[1] {
            if points.iter().any(|point| point[0] <= x && point[1] <= y) {
                area = area.checked_add(1).expect("16 by 16 oracle area fits u16");
            }
        }
    }
    Ok(area)
}

fn oracle_metrics(members: &[IndividualBits]) -> Result<(u16, Vec<u16>), String> {
    let full = oracle_area(members)?;
    let mut contributions = Vec::with_capacity(members.len());
    for drop_index in 0..members.len() {
        let rest: Vec<IndividualBits> = members
            .iter()
            .enumerate()
            .filter(|(index, _)| *index != drop_index)
            .map(|(_, member)| member.clone())
            .collect();
        let without = oracle_area(&rest)?;
        contributions.push(
            full.checked_sub(without)
                .expect("removing a rectangle cannot increase its union area"),
        );
    }
    Ok((full, contributions))
}

fn oracle_checkpoint(members: Vec<IndividualBits>) -> Result<Checkpoint, String> {
    let (area_units, exclusive_contribution_units) = oracle_metrics(&members)?;
    Ok(Checkpoint {
        members,
        hv_bits: area_value(area_units).to_bits(),
        area_units,
        exclusive_contribution_units,
    })
}

fn public_checkpoint(archive: &HvArchive) -> Checkpoint {
    let members: Vec<IndividualBits> = archive.members.iter().map(individual_bits).collect();
    let mut checkpoint = oracle_checkpoint(members)
        .expect("the fixed public archive state stays on the denominator-16 grid");
    checkpoint.hv_bits = archive.hv().to_bits();
    checkpoint
}

/// Recompute the exact capacity-contribution inputs through the PUBLIC
/// production hypervolume function and compare every binary64 result against
/// the independent integer oracle. This proves arithmetic parity for the
/// projected set; it does not observe which internal `insert` branch ran.
fn public_hypervolume_cross_check(
    oracle_expected_pre_capacity: &Checkpoint,
) -> Result<PublicHypervolumeCrossCheck, String> {
    if oracle_expected_pre_capacity.members.len() <= CAPACITY {
        return Err(format!(
            "public-hypervolume-cross-check-needs-over-capacity-set:{}<={CAPACITY}",
            oracle_expected_pre_capacity.members.len()
        ));
    }
    if oracle_expected_pre_capacity
        .exclusive_contribution_units
        .len()
        != oracle_expected_pre_capacity.members.len()
    {
        return Err(format!(
            "public-hypervolume-cross-check-contribution-count:{}!={}",
            oracle_expected_pre_capacity
                .exclusive_contribution_units
                .len(),
            oracle_expected_pre_capacity.members.len()
        ));
    }

    let oracle_full_area_units = oracle_area(&oracle_expected_pre_capacity.members)?;
    if oracle_full_area_units != oracle_expected_pre_capacity.area_units {
        return Err(format!(
            "public-hypervolume-cross-check-oracle-full:{}!={}",
            oracle_full_area_units, oracle_expected_pre_capacity.area_units
        ));
    }
    let front: Vec<Vec<f64>> = oracle_expected_pre_capacity
        .members
        .iter()
        .map(|member| member.f.iter().copied().map(f64::from_bits).collect())
        .collect();
    let full = hypervolume(&front, &REFERENCE);
    let oracle_full = area_value(oracle_full_area_units);
    if full.to_bits() != oracle_full.to_bits() {
        return Err(format!(
            "public-hypervolume-full:0x{:016x}!=oracle-0x{:016x}",
            full.to_bits(),
            oracle_full.to_bits()
        ));
    }

    let mut leave_one_out_hv_bits = Vec::with_capacity(front.len());
    let mut exclusive_contribution_hv_bits = Vec::with_capacity(front.len());
    let mut oracle_leave_one_out_area_units = Vec::with_capacity(front.len());
    let mut oracle_exclusive_contribution_area_units = Vec::with_capacity(front.len());
    for drop_index in 0..front.len() {
        let rest: Vec<Vec<f64>> = front
            .iter()
            .enumerate()
            .filter(|(index, _)| *index != drop_index)
            .map(|(_, point)| point.clone())
            .collect();
        let leave_one_out = hypervolume(&rest, &REFERENCE);
        let oracle_rest: Vec<IndividualBits> = oracle_expected_pre_capacity
            .members
            .iter()
            .enumerate()
            .filter(|(index, _)| *index != drop_index)
            .map(|(_, member)| member.clone())
            .collect();
        let oracle_leave_one_out_units = oracle_area(&oracle_rest)?;
        let oracle_contribution_units = oracle_full_area_units
            .checked_sub(oracle_leave_one_out_units)
            .expect("an exclusive contribution cannot exceed the oracle union area");
        if oracle_contribution_units
            != oracle_expected_pre_capacity.exclusive_contribution_units[drop_index]
        {
            return Err(format!(
                "public-hypervolume-cross-check-oracle-contribution[{drop_index}]:{}!={}",
                oracle_contribution_units,
                oracle_expected_pre_capacity.exclusive_contribution_units[drop_index]
            ));
        }
        let oracle_leave_one_out = area_value(oracle_leave_one_out_units);
        if leave_one_out.to_bits() != oracle_leave_one_out.to_bits() {
            return Err(format!(
                "public-hypervolume-leave-one-out[{drop_index}]:0x{:016x}!=oracle-0x{:016x}",
                leave_one_out.to_bits(),
                oracle_leave_one_out.to_bits()
            ));
        }
        let contribution = full - leave_one_out;
        let oracle_contribution = area_value(oracle_contribution_units);
        if contribution.to_bits() != oracle_contribution.to_bits() {
            return Err(format!(
                "public-hypervolume-contribution[{drop_index}]:0x{:016x}!=oracle-0x{:016x}",
                contribution.to_bits(),
                oracle_contribution.to_bits()
            ));
        }
        leave_one_out_hv_bits.push(leave_one_out.to_bits());
        exclusive_contribution_hv_bits.push(contribution.to_bits());
        oracle_leave_one_out_area_units.push(oracle_leave_one_out_units);
        oracle_exclusive_contribution_area_units.push(oracle_contribution_units);
    }

    Ok(PublicHypervolumeCrossCheck {
        oracle_full_area_units,
        oracle_leave_one_out_area_units,
        oracle_exclusive_contribution_area_units,
        full_hv_bits: full.to_bits(),
        leave_one_out_hv_bits,
        exclusive_contribution_hv_bits,
    })
}

fn oracle_step(
    current: &[IndividualBits],
    candidate: &IndividualBits,
) -> Result<OracleStep, String> {
    let candidate_point = grid_point(candidate)?;
    for (blocker_index, blocker) in current.iter().enumerate() {
        let blocker_point = grid_point(blocker)?;
        let kind = if blocker_point == candidate_point {
            Some(OracleRejectionKind::DuplicateObjective)
        } else if oracle_dominates(blocker_point, candidate_point) {
            Some(OracleRejectionKind::DominatedByMember)
        } else {
            None
        };
        if let Some(kind) = kind {
            let checkpoint = oracle_checkpoint(current.to_vec())?;
            return Ok(OracleStep {
                insert_returned: false,
                expected_pre_capacity: checkpoint.clone(),
                expected_purged: Vec::new(),
                expected_rejection: Some(OracleRejectionExpectation {
                    kind,
                    blocker_index,
                    blocker: blocker.clone(),
                }),
                expected_capacity_eviction: None,
                after: checkpoint,
            });
        }
    }

    let mut expected_purged = Vec::new();
    let mut expected_pre_capacity_members = Vec::with_capacity(current.len() + 1);
    for (prior_index, member) in current.iter().enumerate() {
        if oracle_dominates(candidate_point, grid_point(member)?) {
            expected_purged.push(OraclePurgedMemberExpectation {
                prior_index,
                member: member.clone(),
            });
        } else {
            expected_pre_capacity_members.push(member.clone());
        }
    }
    expected_pre_capacity_members.push(candidate.clone());
    let expected_pre_capacity = oracle_checkpoint(expected_pre_capacity_members.clone())?;

    let (expected_capacity_eviction, after_members) =
        if expected_pre_capacity_members.len() > CAPACITY {
            let minimum = *expected_pre_capacity
                .exclusive_contribution_units
                .iter()
                .min()
                .expect("an over-capacity archive is nonempty");
            let tied_indices: Vec<usize> = expected_pre_capacity
                .exclusive_contribution_units
                .iter()
                .enumerate()
                .filter_map(|(index, &contribution)| (contribution == minimum).then_some(index))
                .collect();
            let selected_index = tied_indices[0];
            let selected = expected_pre_capacity_members[selected_index].clone();
            let mut after = expected_pre_capacity_members;
            after.remove(selected_index);
            (
                Some(OracleCapacityEvictionExpectation {
                    reason: OracleCapacityEvictionReason::LeastExclusiveHypervolumeContributor,
                    selected_index,
                    candidate_evicted: selected.candidate_index == candidate.candidate_index,
                    selected,
                    contribution_units: minimum,
                    tied_indices,
                }),
                after,
            )
        } else {
            (None, expected_pre_capacity_members)
        };

    Ok(OracleStep {
        insert_returned: true,
        expected_pre_capacity,
        expected_purged,
        expected_rejection: None,
        expected_capacity_eviction,
        after: oracle_checkpoint(after_members)?,
    })
}

fn member_candidates(checkpoint: &Checkpoint) -> Vec<usize> {
    checkpoint
        .members
        .iter()
        .map(|member| member.candidate_index)
        .collect()
}

fn expected_rejection(rejection: &Option<OracleRejectionExpectation>) -> Option<ExpectedRejection> {
    rejection.as_ref().map(|rejection| ExpectedRejection {
        kind: rejection.kind,
        blocker_index: rejection.blocker_index,
        blocker_candidate: rejection.blocker.candidate_index,
    })
}

fn expected_capacity_eviction(
    eviction: &Option<OracleCapacityEvictionExpectation>,
) -> Option<ExpectedCapacityEviction> {
    eviction.as_ref().map(|eviction| ExpectedCapacityEviction {
        selected_index: eviction.selected_index,
        selected_candidate: eviction.selected.candidate_index,
        contribution_units: eviction.contribution_units,
        tied_indices: eviction.tied_indices.clone(),
        candidate_evicted: eviction.candidate_evicted,
    })
}

fn branch_mismatch(
    ordinal: usize,
    candidate: &IndividualBits,
    step: &OracleStep,
    expectation: &BranchExpectation,
) -> Option<String> {
    if expectation.purged_prior_indices.len() != expectation.purged_candidates.len() {
        return Some(format!(
            "integer-oracle-branch[{ordinal}].purge-pair-count:{}!={}",
            expectation.purged_prior_indices.len(),
            expectation.purged_candidates.len()
        ));
    }
    let actual = BranchExpectation {
        candidate_index: candidate.candidate_index,
        insert_returned: step.insert_returned,
        rejection: expected_rejection(&step.expected_rejection),
        purged_prior_indices: step
            .expected_purged
            .iter()
            .map(|purged| purged.prior_index)
            .collect(),
        purged_candidates: step
            .expected_purged
            .iter()
            .map(|purged| purged.member.candidate_index)
            .collect(),
        pre_capacity_candidates: member_candidates(&step.expected_pre_capacity),
        capacity_eviction: expected_capacity_eviction(&step.expected_capacity_eviction),
        after_candidates: member_candidates(&step.after),
        after_area_units: step.after.area_units,
    };
    (actual != *expectation).then(|| {
        format!("integer-oracle-branch[{ordinal}]:actual={actual:?};expected={expectation:?}")
    })
}

fn selected_mutation_coordinate() -> MutationCoordinate {
    let mut selector = StreamKey {
        seed: MUTATION_SEED,
        kernel: MUTATION_KERNEL,
        tile: MUTATION_TILE,
    }
    .stream();
    let coordinate_slot =
        usize::try_from(selector.next_below(usize_u64(MUTABLE_STATE_CELLS.len())))
            .expect("mutation coordinate slot fits usize");
    let (transition_index, member_index, objective_index) = MUTABLE_STATE_CELLS[coordinate_slot];
    let mantissa_bit = MUTATION_BIT_BASE
        + u32::try_from(selector.next_below(MUTATION_BIT_COUNT))
            .expect("selected mantissa bit fits u32");
    MutationCoordinate {
        transition_index,
        member_index,
        objective_index,
        mantissa_bit,
        selector_draws: selector.index(),
    }
}

fn bind_candidate(
    mut builder: IdentityBuilder,
    role: &str,
    ordinal: usize,
    candidate: &IndividualBits,
) -> IdentityBuilder {
    builder = builder
        .str("individual-role", role)
        .u64("individual-ordinal", usize_u64(ordinal))
        .u64(
            "individual-candidate-index",
            usize_u64(candidate.candidate_index),
        )
        .u64("individual-x-count", usize_u64(candidate.x.len()))
        .u64("individual-f-count", usize_u64(candidate.f.len()));
    for (coordinate, &bits) in candidate.x.iter().enumerate() {
        builder = builder
            .u64("individual-x-index", usize_u64(coordinate))
            .f64_bits("individual-x", f64::from_bits(bits));
    }
    for (objective, &bits) in candidate.f.iter().enumerate() {
        builder = builder
            .u64("individual-f-index", usize_u64(objective))
            .f64_bits("individual-f", f64::from_bits(bits));
    }
    builder
}

fn bind_checkpoint(
    mut builder: IdentityBuilder,
    role: &str,
    checkpoint: &Checkpoint,
) -> IdentityBuilder {
    builder = builder
        .str("checkpoint-role", role)
        .u64(
            "checkpoint-member-count",
            usize_u64(checkpoint.members.len()),
        )
        .f64_bits("checkpoint-hypervolume", f64::from_bits(checkpoint.hv_bits))
        .u64("checkpoint-area-denominator", u64::from(AREA_DENOMINATOR))
        .u64("checkpoint-area-units", u64::from(checkpoint.area_units))
        .u64(
            "checkpoint-contribution-count",
            usize_u64(checkpoint.exclusive_contribution_units.len()),
        );
    for (member_index, member) in checkpoint.members.iter().enumerate() {
        builder = bind_candidate(builder, "checkpoint-member", member_index, member);
    }
    for (member_index, &contribution) in checkpoint.exclusive_contribution_units.iter().enumerate()
    {
        builder = builder
            .u64("contribution-member-index", usize_u64(member_index))
            .u64("contribution-area-units", u64::from(contribution));
    }
    builder
}

fn bind_public_hypervolume_cross_check(
    mut builder: IdentityBuilder,
    cross_check: &PublicHypervolumeCrossCheck,
) -> IdentityBuilder {
    builder = builder
        .str(
            "public-hypervolume-cross-check-provenance",
            "public-function-recomputation-over-oracle-expected-pre-state-not-insert-branch-observation",
        )
        .u64(
            "public-hypervolume-cross-check-oracle-full-area-units",
            u64::from(cross_check.oracle_full_area_units),
        )
        .u64(
            "public-hypervolume-cross-check-oracle-leave-one-out-count",
            usize_u64(cross_check.oracle_leave_one_out_area_units.len()),
        )
        .u64(
            "public-hypervolume-cross-check-oracle-contribution-count",
            usize_u64(
                cross_check
                    .oracle_exclusive_contribution_area_units
                    .len(),
            ),
        )
        .f64_bits(
            "public-hypervolume-cross-check-full",
            f64::from_bits(cross_check.full_hv_bits),
        )
        .u64(
            "public-hypervolume-cross-check-leave-one-out-count",
            usize_u64(cross_check.leave_one_out_hv_bits.len()),
        )
        .u64(
            "public-hypervolume-cross-check-contribution-count",
            usize_u64(cross_check.exclusive_contribution_hv_bits.len()),
        );
    for (member_index, &units) in cross_check
        .oracle_leave_one_out_area_units
        .iter()
        .enumerate()
    {
        builder = builder
            .u64(
                "public-hypervolume-cross-check-oracle-leave-one-out-member-index",
                usize_u64(member_index),
            )
            .u64(
                "public-hypervolume-cross-check-oracle-leave-one-out-area-units",
                u64::from(units),
            );
    }
    for (member_index, &units) in cross_check
        .oracle_exclusive_contribution_area_units
        .iter()
        .enumerate()
    {
        builder = builder
            .u64(
                "public-hypervolume-cross-check-oracle-contribution-member-index",
                usize_u64(member_index),
            )
            .u64(
                "public-hypervolume-cross-check-oracle-contribution-area-units",
                u64::from(units),
            );
    }
    for (member_index, &bits) in cross_check.leave_one_out_hv_bits.iter().enumerate() {
        builder = builder
            .u64(
                "public-hypervolume-cross-check-leave-one-out-member-index",
                usize_u64(member_index),
            )
            .f64_bits(
                "public-hypervolume-cross-check-leave-one-out",
                f64::from_bits(bits),
            );
    }
    for (member_index, &bits) in cross_check
        .exclusive_contribution_hv_bits
        .iter()
        .enumerate()
    {
        builder = builder
            .u64(
                "public-hypervolume-cross-check-contribution-member-index",
                usize_u64(member_index),
            )
            .f64_bits(
                "public-hypervolume-cross-check-contribution",
                f64::from_bits(bits),
            );
    }
    builder
}

fn fixture_identity() -> ReplayIdentity {
    fixture_identity_for_capacity(CAPACITY)
}

#[allow(clippy::too_many_lines)] // Every fixed input and branch expectation is identity material.
fn fixture_identity_for_capacity(declared_capacity: usize) -> ReplayIdentity {
    let coordinate = selected_mutation_coordinate();
    let specs = candidate_specs();
    let expectations = branch_expectations();
    let mut builder = IdentityBuilder::new(FIXTURE_IDENTITY_KIND)
        .str("algorithm", "fs_dfo::HvArchive")
        .str("algorithm-randomness", "none")
        .str("objective-semantics", "two-objective-minimization")
        .str("archive-order", "survivor-order-then-candidate-append")
        .str(
            "oracle-expected-capacity-eviction-rule",
            "least-exclusive-hypervolume-contribution-strict-less-earliest-index-tie",
        )
        .str(
            "independent-oracle",
            "integer-denominator-16-cell-union-and-leave-one-out-area",
        )
        .str("oracle-branch-provenance", ORACLE_BRANCH_PROVENANCE)
        .flag("internal-branch-authority-observed", false)
        .str(
            "public-hypervolume-cross-check",
            "full-and-every-leave-one-out-for-each-oracle-expected-over-capacity-pre-state",
        )
        .u64("capacity", usize_u64(declared_capacity))
        .u64("objective-count", usize_u64(OBJECTIVES))
        .u64("decision-coordinate-count", usize_u64(DECISION_COORDINATES))
        .u64("grid-denominator", u64::from(GRID_DENOMINATOR))
        .u64("area-denominator", u64::from(AREA_DENOMINATOR))
        .u64("candidate-count", usize_u64(specs.len()))
        .u64("branch-expectation-count", usize_u64(expectations.len()))
        .str("fs-dfo-version", fs_dfo::VERSION)
        .str("fs-obs-version", fs_obs::VERSION)
        .str("fs-rand-version", fs_rand::VERSION)
        .u64(
            "fs-rand-stream-semantics-version",
            u64::from(fs_rand::STREAM_SEMANTICS_VERSION),
        )
        .str(
            "fs-rand-stream-position-domain",
            fs_rand::STREAM_POSITION_IDENTITY_DOMAIN,
        )
        .u64("mutation-seed", MUTATION_SEED)
        .u64("mutation-kernel", u64::from(MUTATION_KERNEL))
        .u64("mutation-tile", u64::from(MUTATION_TILE))
        .u64("mutation-bit-base", u64::from(MUTATION_BIT_BASE))
        .u64("mutation-bit-count", MUTATION_BIT_COUNT)
        .u64(
            "selected-mutation-transition",
            usize_u64(coordinate.transition_index),
        )
        .u64(
            "selected-mutation-member",
            usize_u64(coordinate.member_index),
        )
        .u64(
            "selected-mutation-objective",
            usize_u64(coordinate.objective_index),
        )
        .u64(
            "selected-mutation-mantissa-bit",
            u64::from(coordinate.mantissa_bit),
        )
        .u64("mutation-selector-draws", coordinate.selector_draws)
        .str("fixture-digest-domain", FIXTURE_DIGEST_DOMAIN)
        .str(
            "supplied-fixture-digest-tripwire-domain",
            SUPPLIED_FIXTURE_DIGEST_TRIPWIRE_DOMAIN,
        )
        .u64(
            "wrong-fixture-tripwire-declared-capacity",
            usize_u64(CAPACITY + 1),
        )
        .str("result-digest-domain", RESULT_DIGEST_DOMAIN)
        .str("event-digest-domain", EVENT_DIGEST_DOMAIN)
        .str("no-claims", NO_CLAIMS);
    for (objective, reference) in REFERENCE.into_iter().enumerate() {
        builder = builder
            .u64("reference-objective-index", usize_u64(objective))
            .f64_bits("reference-coordinate", reference)
            .u64(
                "reference-grid-coordinate",
                u64::from(REFERENCE_GRID[objective]),
            );
    }
    for (candidate_index, spec) in specs.into_iter().enumerate() {
        builder = builder
            .u64("candidate-index", usize_u64(candidate_index))
            .str("candidate-label", spec.label);
        builder = bind_candidate(
            builder,
            "ordered-fixture-candidate",
            candidate_index,
            &member_from_spec(candidate_index),
        );
    }
    for (ordinal, expectation) in expectations.iter().enumerate() {
        builder = builder
            .u64("expectation-ordinal", usize_u64(ordinal))
            .u64(
                "expectation-candidate-index",
                usize_u64(expectation.candidate_index),
            )
            .flag("expectation-insert-returned", expectation.insert_returned)
            .u64(
                "expectation-after-area-units",
                u64::from(expectation.after_area_units),
            );
        if let Some(rejection) = &expectation.rejection {
            builder = builder
                .str("expectation-rejection", rejection.kind.name())
                .u64(
                    "expectation-rejection-blocker-index",
                    usize_u64(rejection.blocker_index),
                )
                .u64(
                    "expectation-rejection-blocker-candidate",
                    usize_u64(rejection.blocker_candidate),
                );
        } else {
            builder = builder.str("expectation-rejection", "none");
        }
        assert_eq!(
            expectation.purged_prior_indices.len(),
            expectation.purged_candidates.len(),
            "branch expectation {ordinal} must pair every oracle-expected purge index with one candidate"
        );
        builder = builder
            .u64(
                "expectation-purged-prior-index-count",
                usize_u64(expectation.purged_prior_indices.len()),
            )
            .u64(
                "expectation-purged-candidate-count",
                usize_u64(expectation.purged_candidates.len()),
            )
            .u64(
                "expectation-pre-capacity-count",
                usize_u64(expectation.pre_capacity_candidates.len()),
            )
            .u64(
                "expectation-after-count",
                usize_u64(expectation.after_candidates.len()),
            );
        for (&prior_index, &candidate_index) in expectation
            .purged_prior_indices
            .iter()
            .zip(&expectation.purged_candidates)
        {
            builder = builder
                .u64("expectation-purged-prior-index", usize_u64(prior_index))
                .u64(
                    "expectation-purged-candidate-index",
                    usize_u64(candidate_index),
                );
        }
        for &candidate_index in &expectation.pre_capacity_candidates {
            builder = builder.u64(
                "expectation-pre-capacity-candidate-index",
                usize_u64(candidate_index),
            );
        }
        if let Some(eviction) = &expectation.capacity_eviction {
            builder = builder
                .str(
                    "expectation-capacity-eviction",
                    OracleCapacityEvictionReason::LeastExclusiveHypervolumeContributor.name(),
                )
                .u64(
                    "expectation-eviction-selected-index",
                    usize_u64(eviction.selected_index),
                )
                .u64(
                    "expectation-eviction-selected-candidate",
                    usize_u64(eviction.selected_candidate),
                )
                .u64(
                    "expectation-eviction-contribution-units",
                    u64::from(eviction.contribution_units),
                )
                .flag(
                    "expectation-eviction-candidate-evicted",
                    eviction.candidate_evicted,
                )
                .u64(
                    "expectation-eviction-tie-count",
                    usize_u64(eviction.tied_indices.len()),
                );
            for &tie_index in &eviction.tied_indices {
                builder = builder.u64("expectation-eviction-tie-index", usize_u64(tie_index));
            }
        } else {
            builder = builder.str("expectation-capacity-eviction", "none");
        }
        for &candidate_index in &expectation.after_candidates {
            builder = builder.u64(
                "expectation-after-candidate-index",
                usize_u64(candidate_index),
            );
        }
    }
    for (catalog_index, &(transition, member, objective)) in MUTABLE_STATE_CELLS.iter().enumerate()
    {
        builder = builder
            .u64("mutation-catalog-index", usize_u64(catalog_index))
            .u64("mutation-catalog-transition", usize_u64(transition))
            .u64("mutation-catalog-member", usize_u64(member))
            .u64("mutation-catalog-objective", usize_u64(objective));
    }
    builder.finish()
}

fn fixture_digest(fixture: &ReplayIdentity) -> ContentHash {
    hash_domain(FIXTURE_DIGEST_DOMAIN, fixture.canonical_bytes())
}

#[allow(clippy::too_many_lines)] // Complete trajectory identity is intentionally explicit.
fn result_identity(
    fixture: &ReplayIdentity,
    strong_fixture: ContentHash,
    record: &StudyRecord,
) -> ReplayIdentity {
    let mut builder = IdentityBuilder::new(RESULT_IDENTITY_KIND)
        .child("fixture-compatibility-root", fixture)
        .bytes("fixture-canonical-bytes", fixture.canonical_bytes())
        .bytes("fixture-blake3", strong_fixture.as_bytes())
        .str("oracle-branch-provenance", ORACLE_BRANCH_PROVENANCE)
        .flag("internal-branch-authority-observed", false)
        .u64("transition-count", usize_u64(record.transitions.len()));
    builder = bind_checkpoint(builder, "initial", &record.initial);
    for transition in &record.transitions {
        builder = builder
            .u64("transition-ordinal", usize_u64(transition.ordinal))
            .flag("transition-insert-returned", transition.insert_returned);
        builder = bind_candidate(
            builder,
            "transition-candidate",
            transition.ordinal,
            &transition.candidate,
        );
        builder = bind_checkpoint(builder, "transition-before", &transition.before);
        builder = bind_checkpoint(
            builder,
            "transition-oracle-expected-pre-capacity",
            &transition.oracle_expected_pre_capacity,
        );
        builder = builder.u64(
            "transition-oracle-expected-purged-count",
            usize_u64(transition.oracle_expected_purged.len()),
        );
        for purged in &transition.oracle_expected_purged {
            builder = builder.u64(
                "transition-oracle-expected-purged-prior-index",
                usize_u64(purged.prior_index),
            );
            builder = bind_candidate(
                builder,
                "transition-oracle-expected-purged-member",
                purged.prior_index,
                &purged.member,
            );
        }
        if let Some(rejection) = &transition.oracle_expected_rejection {
            builder = builder
                .str(
                    "transition-oracle-expected-rejection",
                    rejection.kind.name(),
                )
                .u64(
                    "transition-oracle-expected-rejection-blocker-index",
                    usize_u64(rejection.blocker_index),
                );
            builder = bind_candidate(
                builder,
                "transition-oracle-expected-rejection-blocker",
                rejection.blocker_index,
                &rejection.blocker,
            );
        } else {
            builder = builder.str("transition-oracle-expected-rejection", "none");
        }
        if let Some(eviction) = &transition.oracle_expected_capacity_eviction {
            builder = builder
                .str(
                    "transition-oracle-expected-capacity-eviction",
                    eviction.reason.name(),
                )
                .u64(
                    "transition-oracle-expected-eviction-selected-index",
                    usize_u64(eviction.selected_index),
                )
                .u64(
                    "transition-oracle-expected-eviction-contribution-units",
                    u64::from(eviction.contribution_units),
                )
                .flag(
                    "transition-oracle-expected-eviction-candidate-evicted",
                    eviction.candidate_evicted,
                )
                .u64(
                    "transition-oracle-expected-eviction-tie-count",
                    usize_u64(eviction.tied_indices.len()),
                );
            builder = bind_candidate(
                builder,
                "transition-oracle-expected-evicted-member",
                eviction.selected_index,
                &eviction.selected,
            );
            for &tie_index in &eviction.tied_indices {
                builder = builder.u64(
                    "transition-oracle-expected-eviction-tie-index",
                    usize_u64(tie_index),
                );
            }
        } else {
            builder = builder.str("transition-oracle-expected-capacity-eviction", "none");
        }
        if let Some(cross_check) = &transition.public_hypervolume_cross_check {
            builder = builder.str("transition-public-hypervolume-cross-check", "present");
            builder = bind_public_hypervolume_cross_check(builder, cross_check);
        } else {
            builder = builder.str("transition-public-hypervolume-cross-check", "none");
        }
        builder = bind_checkpoint(builder, "transition-after", &transition.after);
    }
    builder.finish()
}

fn result_digest(result: &ReplayIdentity) -> ContentHash {
    hash_domain(RESULT_DIGEST_DOMAIN, result.canonical_bytes())
}

fn event_digest(event: &Event) -> ContentHash {
    hash_domain(
        EVENT_DIGEST_DOMAIN,
        event.content_identity().canonical_bytes(),
    )
}

fn run_study() -> StudyRun {
    let fixture = fixture_identity();
    let fixture_digest = fixture_digest(&fixture);
    let mut archive = HvArchive::new(CAPACITY, REFERENCE.to_vec());
    let initial = public_checkpoint(&archive);
    let mut transitions = Vec::with_capacity(CANDIDATE_COUNT);
    for (ordinal, spec) in candidate_specs().into_iter().enumerate() {
        let before = public_checkpoint(&archive);
        let candidate = member_from_spec(ordinal);
        let oracle = oracle_step(&before.members, &candidate)
            .expect("fixed archive transition stays inside the integer oracle domain");
        let public_hypervolume_cross_check =
            (oracle.expected_pre_capacity.members.len() > CAPACITY)
                .then(|| {
                    public_hypervolume_cross_check(&oracle.expected_pre_capacity).expect(
                        "public production hypervolume must match every exact integer-oracle capacity contribution",
                    )
                });
        let insert_returned = archive.insert(spec.individual());
        let after = public_checkpoint(&archive);
        transitions.push(TransitionRecord {
            ordinal,
            candidate,
            insert_returned,
            before,
            oracle_expected_pre_capacity: oracle.expected_pre_capacity,
            oracle_expected_purged: oracle.expected_purged,
            oracle_expected_rejection: oracle.expected_rejection,
            oracle_expected_capacity_eviction: oracle.expected_capacity_eviction,
            public_hypervolume_cross_check,
            after,
        });
    }
    let record = StudyRecord {
        initial,
        transitions,
    };
    let result = result_identity(&fixture, fixture_digest, &record);
    let result_digest = result_digest(&result);
    StudyRun {
        fixture,
        fixture_digest,
        record,
        result,
        result_digest,
    }
}

fn member_mismatch(
    path: &str,
    found: &IndividualBits,
    expected: &IndividualBits,
) -> Option<String> {
    if found.candidate_index != expected.candidate_index {
        return Some(format!(
            "integer-oracle-{path}.candidate-index:{}!={}",
            found.candidate_index, expected.candidate_index
        ));
    }
    if found.x.len() != expected.x.len() {
        return Some(format!(
            "integer-oracle-{path}.x-count:{}!={}",
            found.x.len(),
            expected.x.len()
        ));
    }
    for (coordinate, (&found_bits, &expected_bits)) in found.x.iter().zip(&expected.x).enumerate() {
        if found_bits != expected_bits {
            return Some(format!(
                "integer-oracle-{path}.x[{coordinate}]:0x{found_bits:016x}!=0x{expected_bits:016x}"
            ));
        }
    }
    if found.f.len() != expected.f.len() {
        return Some(format!(
            "integer-oracle-{path}.f-count:{}!={}",
            found.f.len(),
            expected.f.len()
        ));
    }
    for (objective, (&found_bits, &expected_bits)) in found.f.iter().zip(&expected.f).enumerate() {
        if found_bits != expected_bits {
            return Some(format!(
                "integer-oracle-{path}.f[{objective}]:0x{found_bits:016x}!=0x{expected_bits:016x}"
            ));
        }
    }
    None
}

fn checkpoint_mismatch(path: &str, found: &Checkpoint, expected: &Checkpoint) -> Option<String> {
    if found.members.len() != expected.members.len() {
        return Some(format!(
            "integer-oracle-{path}.member-count:{}!={}",
            found.members.len(),
            expected.members.len()
        ));
    }
    for (member_index, (found_member, expected_member)) in
        found.members.iter().zip(&expected.members).enumerate()
    {
        if let Some(mismatch) = member_mismatch(
            &format!("{path}.members[{member_index}]"),
            found_member,
            expected_member,
        ) {
            return Some(mismatch);
        }
    }
    if found.hv_bits != expected.hv_bits {
        return Some(format!(
            "integer-oracle-{path}.hv:0x{:016x}!=0x{:016x}",
            found.hv_bits, expected.hv_bits
        ));
    }
    if found.area_units != expected.area_units {
        return Some(format!(
            "integer-oracle-{path}.area-units:{}!={}",
            found.area_units, expected.area_units
        ));
    }
    if found.exclusive_contribution_units != expected.exclusive_contribution_units {
        return Some(format!(
            "integer-oracle-{path}.contributions:{:?}!={:?}",
            found.exclusive_contribution_units, expected.exclusive_contribution_units
        ));
    }
    None
}

#[allow(clippy::too_many_lines)] // First-failure ordering is part of the red receipt.
fn semantic_mismatch(record: &StudyRecord) -> Option<String> {
    if record.transitions.len() != CANDIDATE_COUNT {
        return Some(format!(
            "integer-oracle-transition-count:{}!={CANDIDATE_COUNT}",
            record.transitions.len()
        ));
    }
    let mut current = oracle_checkpoint(Vec::new()).expect("empty integer archive is valid");
    if let Some(mismatch) = checkpoint_mismatch("initial", &record.initial, &current) {
        return Some(mismatch);
    }
    let expectations = branch_expectations();
    for (ordinal, transition) in record.transitions.iter().enumerate() {
        if transition.ordinal != ordinal {
            return Some(format!(
                "integer-oracle-transition[{ordinal}].ordinal:{}!={ordinal}",
                transition.ordinal
            ));
        }
        let expected_candidate = member_from_spec(ordinal);
        if let Some(mismatch) = member_mismatch(
            &format!("transition[{ordinal}].candidate"),
            &transition.candidate,
            &expected_candidate,
        ) {
            return Some(mismatch);
        }
        if let Some(mismatch) = checkpoint_mismatch(
            &format!("transition[{ordinal}].before"),
            &transition.before,
            &current,
        ) {
            return Some(mismatch);
        }
        let step = match oracle_step(&current.members, &expected_candidate) {
            Ok(step) => step,
            Err(error) => {
                return Some(format!(
                    "integer-oracle-transition[{ordinal}]-construction:{error}"
                ));
            }
        };
        if let Some(mismatch) =
            branch_mismatch(ordinal, &expected_candidate, &step, &expectations[ordinal])
        {
            return Some(mismatch);
        }
        if transition.insert_returned != step.insert_returned {
            return Some(format!(
                "integer-oracle-transition[{ordinal}].insert-returned:{}!={}",
                transition.insert_returned, step.insert_returned
            ));
        }
        if let Some(mismatch) = checkpoint_mismatch(
            &format!("transition[{ordinal}].oracle-expected-pre-capacity"),
            &transition.oracle_expected_pre_capacity,
            &step.expected_pre_capacity,
        ) {
            return Some(mismatch);
        }
        if transition.oracle_expected_purged != step.expected_purged {
            return Some(format!(
                "integer-oracle-transition[{ordinal}].oracle-expected-purged:{:?}!={:?}",
                transition.oracle_expected_purged, step.expected_purged
            ));
        }
        if transition.oracle_expected_rejection != step.expected_rejection {
            return Some(format!(
                "integer-oracle-transition[{ordinal}].oracle-expected-rejection:{:?}!={:?}",
                transition.oracle_expected_rejection, step.expected_rejection
            ));
        }
        if transition.oracle_expected_capacity_eviction != step.expected_capacity_eviction {
            return Some(format!(
                "integer-oracle-transition[{ordinal}].oracle-expected-capacity-eviction:{:?}!={:?}",
                transition.oracle_expected_capacity_eviction, step.expected_capacity_eviction
            ));
        }
        let expected_public_cross_check = if step.expected_pre_capacity.members.len() > CAPACITY {
            match public_hypervolume_cross_check(&step.expected_pre_capacity) {
                Ok(cross_check) => Some(cross_check),
                Err(error) => {
                    return Some(format!(
                        "public-hypervolume-cross-check-transition[{ordinal}]:{error}"
                    ));
                }
            }
        } else {
            None
        };
        if transition.public_hypervolume_cross_check != expected_public_cross_check {
            return Some(format!(
                "public-hypervolume-cross-check-transition[{ordinal}]:{:?}!={:?}",
                transition.public_hypervolume_cross_check, expected_public_cross_check
            ));
        }
        if let Some(mismatch) = checkpoint_mismatch(
            &format!("transition[{ordinal}].after"),
            &transition.after,
            &step.after,
        ) {
            return Some(mismatch);
        }
        current = step.after;
    }
    None
}

fn validate_payload(run: &StudyRun) -> Result<(), AdmissionError> {
    let computed_fixture_digest = fixture_digest(&run.fixture);
    if computed_fixture_digest != run.fixture_digest {
        return Err(AdmissionError::SuppliedFixtureDigestMismatch {
            declared: digest_bytes(run.fixture_digest),
            computed: digest_bytes(computed_fixture_digest),
        });
    }
    let expected_fixture = fixture_identity();
    let expected_fixture_digest = fixture_digest(&expected_fixture);
    if run.fixture.canonical_bytes() != expected_fixture.canonical_bytes() {
        return Err(AdmissionError::RetainedFixtureIdentityMismatch {
            expected: digest_bytes(expected_fixture_digest),
            found: digest_bytes(computed_fixture_digest),
        });
    }
    let computed_result = result_identity(&run.fixture, run.fixture_digest, &run.record);
    let computed_result_digest = result_digest(&computed_result);
    if run.result.canonical_bytes() != computed_result.canonical_bytes()
        || run.result_digest != computed_result_digest
    {
        return Err(AdmissionError::ResultPayloadIdentityMismatch {
            declared: digest_bytes(run.result_digest),
            computed: digest_bytes(computed_result_digest),
        });
    }
    Ok(())
}

fn validate_semantics(run: &StudyRun) -> Result<(), AdmissionError> {
    match semantic_mismatch(&run.record) {
        Some(mismatch) => Err(AdmissionError::SemanticInconsistency(mismatch)),
        None => Ok(()),
    }
}

fn admit_reference(run: &StudyRun, reference: &StudyRun) -> Result<(), AdmissionError> {
    validate_payload(run)?;
    if run.result.canonical_bytes() == reference.result.canonical_bytes()
        && run.result_digest == reference.result_digest
    {
        Ok(())
    } else {
        Err(AdmissionError::ReferenceIdentityMismatch {
            expected: digest_bytes(reference.result_digest),
            found: digest_bytes(run.result_digest),
        })
    }
}

fn reseal(run: &mut StudyRun) {
    run.result = result_identity(&run.fixture, run.fixture_digest, &run.record);
    run.result_digest = result_digest(&run.result);
}

fn exact_state_bit_delta(reference: &StudyRun, mutant: &StudyRun, mutation: Mutation) -> bool {
    let coordinate = mutation.coordinate;
    let Some(mask) = 1u64.checked_shl(coordinate.mantissa_bit) else {
        return false;
    };
    if reference.fixture != mutant.fixture
        || reference.fixture_digest != mutant.fixture_digest
        || coordinate != selected_mutation_coordinate()
        || mutation.before ^ mutation.after != mask
    {
        return false;
    }
    let Some(reference_bits) = reference
        .record
        .transitions
        .get(coordinate.transition_index)
        .and_then(|transition| transition.after.members.get(coordinate.member_index))
        .and_then(|member| member.f.get(coordinate.objective_index))
        .copied()
    else {
        return false;
    };
    let Some(mutant_bits) = mutant
        .record
        .transitions
        .get(coordinate.transition_index)
        .and_then(|transition| transition.after.members.get(coordinate.member_index))
        .and_then(|member| member.f.get(coordinate.objective_index))
        .copied()
    else {
        return false;
    };
    if reference_bits != mutation.before || mutant_bits != mutation.after {
        return false;
    }
    let mut expected = reference.record.clone();
    expected.transitions[coordinate.transition_index]
        .after
        .members[coordinate.member_index]
        .f[coordinate.objective_index] = mutation.after;
    expected == mutant.record
}

fn seeded_corruption(reference: &StudyRun) -> SeededCorruption {
    let coordinate = selected_mutation_coordinate();
    let mut run = reference.clone();
    let target = &mut run.record.transitions[coordinate.transition_index]
        .after
        .members[coordinate.member_index]
        .f[coordinate.objective_index];
    let before = *target;
    let after = before ^ (1u64 << coordinate.mantissa_bit);
    *target = after;
    let stale_error =
        validate_payload(&run).expect_err("unsealed checkpoint mutation must be stale");
    reseal(&mut run);
    let reference_error = admit_reference(&run, reference)
        .expect_err("resealed checkpoint mutation must miss the retained reference");
    let semantic_error = validate_semantics(&run)
        .expect_err("resealed checkpoint mutation must violate the integer lifecycle oracle");
    SeededCorruption {
        run,
        mutation: Mutation {
            seed: MUTATION_SEED,
            kernel: MUTATION_KERNEL,
            tile: MUTATION_TILE,
            coordinate,
            before,
            after,
        },
        stale_error,
        reference_error,
        semantic_error,
    }
}

fn green_receipt(run: &StudyRun) -> Event {
    let tie_cross_check = run.record.transitions[2]
        .public_hypervolume_cross_check
        .as_ref()
        .expect("tie pre-state retains the public hypervolume cross-check");
    let self_eviction_cross_check = run.record.transitions[6]
        .public_hypervolume_cross_check
        .as_ref()
        .expect("self-eviction pre-state retains the public hypervolume cross-check");
    let final_checkpoint = &run
        .record
        .transitions
        .last()
        .expect("fixed lifecycle has transitions")
        .after;
    let mut emitter = Emitter::new(SUITE, CASE);
    emitter.emit(
        Severity::Info,
        EventKind::Custom {
            name: "hv-archive-full-lifecycle-replay-receipt".to_string(),
            json: format!(
                concat!(
                    "{{\"fixture_identity\":\"{}\",\"fixture_blake3\":\"{}\",",
                    "\"result_identity\":\"{}\",\"result_blake3\":\"{}\",",
                    "\"algorithm\":\"fs_dfo::HvArchive\",\"algorithm_seed\":null,",
                    "\"capacity\":{},\"reference_bits\":[\"0x{:016x}\",\"0x{:016x}\"],",
                    "\"grid_denominator\":{},\"area_denominator\":{},",
                    "\"branch_metadata_provenance\":\"{}\",",
                    "\"internal_branch_authority_observed\":false,",
                    "\"public_hypervolume_cross_checks\":[",
                    "{{\"transition\":2,\"oracle_full_area_units\":{},",
                    "\"oracle_leave_one_out_area_units\":{:?},",
                    "\"oracle_exclusive_contribution_area_units\":{:?},",
                    "\"full_hv_bits\":\"0x{:016x}\",",
                    "\"leave_one_out_hv_bits\":{},\"exclusive_contribution_hv_bits\":{}}},",
                    "{{\"transition\":6,\"oracle_full_area_units\":{},",
                    "\"oracle_leave_one_out_area_units\":{:?},",
                    "\"oracle_exclusive_contribution_area_units\":{:?},",
                    "\"full_hv_bits\":\"0x{:016x}\",",
                    "\"leave_one_out_hv_bits\":{},\"exclusive_contribution_hv_bits\":{}}}],",
                    "\"transition_count\":{},\"final_member_candidates\":{:?},",
                    "\"final_hv_bits\":\"0x{:016x}\",\"final_area_units\":{},",
                    "\"mutation_seed\":\"0x{:016x}\",",
                    "\"versions\":{{\"fs_dfo\":\"{}\",\"fs_obs\":\"{}\",",
                    "\"fs_rand\":\"{}\",\"stream_semantics\":{}}},",
                    "\"no_claims\":[\"arbitrary-dimensions-capacities-fronts\",",
                    "\"malformed-nonfinite-outside-reference-policy\",",
                    "\"mc-contribution-eviction\",\"optimizer-convergence-archive-quality\",",
                    "\"sealed-public-state\",\"allocation-Cx-cancellation\",",
                    "\"cross-ISA-authority\",\"authenticated-persistence\",",
                    "\"performance\",",
                    "\"internal-purge-rejection-eviction-branch-authority\"]}}"
                ),
                run.fixture.hex(),
                run.fixture_digest.to_hex(),
                run.result.hex(),
                run.result_digest.to_hex(),
                CAPACITY,
                REFERENCE[0].to_bits(),
                REFERENCE[1].to_bits(),
                GRID_DENOMINATOR,
                AREA_DENOMINATOR,
                ORACLE_BRANCH_PROVENANCE,
                tie_cross_check.oracle_full_area_units,
                &tie_cross_check.oracle_leave_one_out_area_units,
                &tie_cross_check.oracle_exclusive_contribution_area_units,
                tie_cross_check.full_hv_bits,
                bit_hex_json(&tie_cross_check.leave_one_out_hv_bits),
                bit_hex_json(&tie_cross_check.exclusive_contribution_hv_bits),
                self_eviction_cross_check.oracle_full_area_units,
                &self_eviction_cross_check.oracle_leave_one_out_area_units,
                &self_eviction_cross_check.oracle_exclusive_contribution_area_units,
                self_eviction_cross_check.full_hv_bits,
                bit_hex_json(&self_eviction_cross_check.leave_one_out_hv_bits),
                bit_hex_json(&self_eviction_cross_check.exclusive_contribution_hv_bits),
                run.record.transitions.len(),
                member_candidates(final_checkpoint),
                final_checkpoint.hv_bits,
                final_checkpoint.area_units,
                MUTATION_SEED,
                fs_dfo::VERSION,
                fs_obs::VERSION,
                fs_rand::VERSION,
                fs_rand::STREAM_SEMANTICS_VERSION,
            ),
        },
        None,
    )
}

fn green_verdict(run: &StudyRun) -> Event {
    let mut emitter = Emitter::new(SUITE, format!("{CASE}/verdict"));
    emitter.emit(
        Severity::Info,
        EventKind::ConformanceCase {
            suite: SUITE.to_string(),
            case: CASE.to_string(),
            pass: true,
            detail: format!(
                "fixture={}; result={}; blake3={}; transitions={}; public-outcomes=tie-removes-A,self-insertion-E-returns-true-with-state-no-op,final-removes-C-D; oracle-expected-contributions=tie:[16,16,16]/256,self:[8,54,5]/256; public-hypervolume-cross-checks=2; branch-metadata-provenance={ORACLE_BRANCH_PROVENANCE}; internal-branch-authority-observed=false; no-claim=internal-purge-rejection-eviction-branch-authority",
                run.fixture.hex(),
                run.result.hex(),
                run.result_digest.to_hex(),
                run.record.transitions.len(),
            ),
            seed: 0,
        },
        None,
    )
}

fn corruption_event(reference: &StudyRun, corruption: &SeededCorruption) -> Event {
    let mutation = corruption.mutation;
    let coordinate = mutation.coordinate;
    let detail = format!(
        "reference={}; mutant={}; seed=0x{:016x}; kernel=0x{:04x}; tile={}; selector_draws={}; target=transition[{}].after.members[{}].f[{}]; mantissa_bit={}; before=0x{:016x}; after=0x{:016x}; stale={:?}; reference_gate={:?}; first_semantic_mismatch={:?}; branch-metadata-provenance={ORACLE_BRANCH_PROVENANCE}; internal-branch-authority-observed=false; no-claim=internal-purge-rejection-eviction-branch-authority",
        reference.result_digest.to_hex(),
        corruption.run.result_digest.to_hex(),
        mutation.seed,
        mutation.kernel,
        mutation.tile,
        coordinate.selector_draws,
        coordinate.transition_index,
        coordinate.member_index,
        coordinate.objective_index,
        coordinate.mantissa_bit,
        mutation.before,
        mutation.after,
        corruption.stale_error,
        corruption.reference_error,
        corruption.semantic_error,
    );
    let mut emitter = Emitter::new(SUITE, RED_CASE);
    emitter.emit(
        Severity::Error,
        EventKind::ConformanceCase {
            suite: SUITE.to_string(),
            case: RED_CASE.to_string(),
            pass: false,
            detail,
            seed: MUTATION_SEED,
        },
        None,
    )
}

fn assert_mergeable(event: &Event) {
    let EventKind::ConformanceCase {
        case, pass, detail, ..
    } = &event.kind
    else {
        panic!("merge gate accepts only ConformanceCase evidence");
    };
    assert!(*pass, "merge gate refused {case}: {detail}");
}

fn assert_event_pair(first: &Event, second: &Event, label: &str) {
    assert_eq!(
        first.content_identity().canonical_bytes(),
        second.content_identity().canonical_bytes(),
        "{label} content must replay byte-for-byte"
    );
    assert_eq!(event_digest(first), event_digest(second));
    for event in [first, second] {
        fs_obs::lint_failure_record(event).expect("HvArchive evidence retains replay inputs");
        fs_obs::validate_line(&event.to_jsonl()).expect("HvArchive evidence is fs-obs wire-valid");
        let receipt = event.content_identity_receipt();
        event
            .admit_content_identity(&receipt)
            .expect("HvArchive evidence content identity admits exactly");
    }
}

#[test]
#[allow(clippy::too_many_lines)] // One causal test spans the lifecycle and every refusal gate.
fn hv_archive_full_lifecycle_replays_and_seeded_failure_is_refused() {
    let specs = candidate_specs();
    for left in 0..specs.len() {
        for right in (left + 1)..specs.len() {
            assert_ne!(
                specs[left].x_bits, specs[right].x_bits,
                "every candidate needs a distinct decision payload"
            );
        }
    }

    let original = run_study();
    let replay = run_study();
    assert_eq!(validate_payload(&original), Ok(()));
    assert_eq!(validate_payload(&replay), Ok(()));
    assert_eq!(validate_semantics(&original), Ok(()));
    assert_eq!(validate_semantics(&replay), Ok(()));
    assert_eq!(admit_reference(&original, &replay), Ok(()));
    assert_eq!(admit_reference(&replay, &original), Ok(()));

    let mut wrong_supplied_digest = original.clone();
    wrong_supplied_digest.fixture_digest = hash_domain(
        SUPPLIED_FIXTURE_DIGEST_TRIPWIRE_DOMAIN,
        wrong_supplied_digest.fixture.canonical_bytes(),
    );
    let supplied_digest_error = validate_payload(&wrong_supplied_digest)
        .expect_err("wrong supplied fixture digest must fail before retained-fixture admission");
    assert!(matches!(
        supplied_digest_error,
        AdmissionError::SuppliedFixtureDigestMismatch { declared, computed }
            if declared == *wrong_supplied_digest.fixture_digest.as_bytes()
                && computed == *original.fixture_digest.as_bytes()
                && declared != computed
    ));

    let mut wrong_retained_fixture = original.clone();
    wrong_retained_fixture.fixture = fixture_identity_for_capacity(CAPACITY + 1);
    wrong_retained_fixture.fixture_digest = fixture_digest(&wrong_retained_fixture.fixture);
    reseal(&mut wrong_retained_fixture);
    let retained_fixture_error = validate_payload(&wrong_retained_fixture)
        .expect_err("self-consistent wrong fixture must miss the retained canonical fixture");
    assert!(matches!(
        retained_fixture_error,
        AdmissionError::RetainedFixtureIdentityMismatch { expected, found }
            if expected == *original.fixture_digest.as_bytes()
                && found == *wrong_retained_fixture.fixture_digest.as_bytes()
                && expected != found
    ));

    assert_eq!(original.record, replay.record);
    assert_eq!(original.fixture, replay.fixture);
    assert_eq!(original.fixture_digest, replay.fixture_digest);
    assert_eq!(original.result, replay.result);
    assert_eq!(original.result_digest, replay.result_digest);
    assert_eq!(
        original.result.canonical_bytes(),
        replay.result.canonical_bytes(),
        "complete archive trajectory frames must replay byte-for-byte"
    );

    let transitions = &original.record.transitions;
    assert_eq!(
        transitions
            .iter()
            .filter(|transition| transition.public_hypervolume_cross_check.is_some())
            .map(|transition| transition.ordinal)
            .collect::<Vec<_>>(),
        vec![2, 6],
        "every and only oracle-expected over-capacity pre-state needs a public hypervolume cross-check"
    );
    let tie_eviction = transitions[2]
        .oracle_expected_capacity_eviction
        .as_ref()
        .expect("the integer oracle expects C to trigger the exact contribution tie");
    assert_eq!(
        transitions[2]
            .oracle_expected_pre_capacity
            .exclusive_contribution_units,
        vec![16, 16, 16]
    );
    assert_eq!(tie_eviction.selected_index, 0);
    assert_eq!(tie_eviction.selected.candidate_index, 0);
    assert_eq!(tie_eviction.tied_indices, vec![0, 1, 2]);
    let tie_public_cross_check = transitions[2]
        .public_hypervolume_cross_check
        .as_ref()
        .expect("the three-member tie pre-state requires a public hypervolume cross-check");
    assert_eq!(tie_public_cross_check.oracle_full_area_units, 96);
    assert_eq!(
        tie_public_cross_check.oracle_leave_one_out_area_units,
        vec![80, 80, 80]
    );
    assert_eq!(
        tie_public_cross_check.oracle_exclusive_contribution_area_units,
        vec![16, 16, 16]
    );
    assert_eq!(
        tie_public_cross_check.full_hv_bits,
        area_value(96).to_bits()
    );
    assert_eq!(
        tie_public_cross_check.leave_one_out_hv_bits,
        vec![area_value(80).to_bits(); 3]
    );
    assert_eq!(
        tie_public_cross_check.exclusive_contribution_hv_bits,
        vec![area_value(16).to_bits(); 3]
    );

    for ordinal in [3usize, 4] {
        assert!(!transitions[ordinal].insert_returned);
        assert_eq!(transitions[ordinal].before, transitions[ordinal].after);
    }
    assert_eq!(
        transitions[5]
            .oracle_expected_purged
            .iter()
            .map(|purged| purged.member.candidate_index)
            .collect::<Vec<_>>(),
        vec![1]
    );
    let self_eviction = transitions[6]
        .oracle_expected_capacity_eviction
        .as_ref()
        .expect("the integer oracle expects E to be the unique least contributor");
    assert!(transitions[6].insert_returned);
    assert!(self_eviction.candidate_evicted);
    assert_eq!(self_eviction.selected.candidate_index, 6);
    assert_eq!(self_eviction.contribution_units, 5);
    assert_eq!(self_eviction.tied_indices, vec![2]);
    let self_eviction_public_cross_check = transitions[6]
        .public_hypervolume_cross_check
        .as_ref()
        .expect(
            "the three-member self-eviction pre-state requires a public hypervolume cross-check",
        );
    assert_eq!(self_eviction_public_cross_check.oracle_full_area_units, 113);
    assert_eq!(
        self_eviction_public_cross_check.oracle_leave_one_out_area_units,
        vec![105, 59, 108]
    );
    assert_eq!(
        self_eviction_public_cross_check.oracle_exclusive_contribution_area_units,
        vec![8, 54, 5]
    );
    assert_eq!(
        self_eviction_public_cross_check.full_hv_bits,
        area_value(113).to_bits()
    );
    assert_eq!(
        self_eviction_public_cross_check.leave_one_out_hv_bits,
        vec![
            area_value(105).to_bits(),
            area_value(59).to_bits(),
            area_value(108).to_bits(),
        ]
    );
    assert_eq!(
        self_eviction_public_cross_check.exclusive_contribution_hv_bits,
        vec![
            area_value(8).to_bits(),
            area_value(54).to_bits(),
            area_value(5).to_bits(),
        ]
    );
    assert_eq!(transitions[6].before, transitions[6].after);
    assert_eq!(
        transitions[7]
            .oracle_expected_purged
            .iter()
            .map(|purged| purged.member.candidate_index)
            .collect::<Vec<_>>(),
        vec![2, 5]
    );
    assert_eq!(member_candidates(&transitions[7].after), vec![7]);
    assert_eq!(transitions[7].after.area_units, 196);

    let first_receipt = green_receipt(&original);
    let second_receipt = green_receipt(&replay);
    assert_event_pair(&first_receipt, &second_receipt, "green HvArchive receipt");
    println!("{}", first_receipt.to_jsonl());

    let first_green = green_verdict(&original);
    let second_green = green_verdict(&replay);
    assert_event_pair(&first_green, &second_green, "green HvArchive verdict");
    assert_mergeable(&first_green);
    assert_mergeable(&second_green);
    println!("{}", first_green.to_jsonl());

    let first = seeded_corruption(&original);
    let second = seeded_corruption(&replay);
    assert_eq!(
        first, second,
        "seeded retained-state corruption must replay exactly"
    );
    assert!(
        exact_state_bit_delta(&original, &first.run, first.mutation),
        "mutation must alter exactly one retained checkpoint objective bit"
    );
    assert_eq!(
        validate_payload(&first.run),
        Ok(()),
        "resealed state mutation must be internally identity-consistent"
    );
    assert!(f64::from_bits(first.mutation.after).is_finite());
    assert!(matches!(
        &first.stale_error,
        AdmissionError::ResultPayloadIdentityMismatch { declared, computed }
            if declared == original.result_digest.as_bytes()
                && computed == first.run.result_digest.as_bytes()
    ));
    assert!(matches!(
        &first.reference_error,
        AdmissionError::ReferenceIdentityMismatch { expected, found }
            if expected == original.result_digest.as_bytes()
                && found == first.run.result_digest.as_bytes()
    ));
    let coordinate = first.mutation.coordinate;
    let expected_mismatch_path = format!(
        "transition[{}].after.members[{}].f[{}]",
        coordinate.transition_index, coordinate.member_index, coordinate.objective_index
    );
    assert!(matches!(
        &first.semantic_error,
        AdmissionError::SemanticInconsistency(mismatch)
            if mismatch.contains(&expected_mismatch_path)
                && mismatch.starts_with("integer-oracle-")
    ));

    let first_red = corruption_event(&original, &first);
    let second_red = corruption_event(&replay, &second);
    assert_event_pair(&first_red, &second_red, "red HvArchive evidence");
    println!("{}", first_red.to_jsonl());

    let panic = catch_unwind(|| assert_mergeable(&first_red))
        .expect_err("merge gate must refuse the seeded checkpoint corruption");
    let message = panic
        .downcast_ref::<String>()
        .map(String::as_str)
        .or_else(|| panic.downcast_ref::<&str>().copied())
        .expect("merge-gate panic carries text");
    assert!(message.contains(RED_CASE));
    assert!(message.contains(&format!("0x{MUTATION_SEED:016x}")));
    assert!(message.contains(&expected_mismatch_path));
    assert!(message.contains("ReferenceIdentityMismatch"));
    assert!(message.contains("SemanticInconsistency"));
}
