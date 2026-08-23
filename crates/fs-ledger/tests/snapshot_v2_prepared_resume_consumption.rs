//! No-mock terminal consumption of the fs-exec snapshot v2
//! [`snapshot_v2::PreparedResume`] path through the real Design Ledger.
//!
//! Beads: frankensim-sj31i.52.5.2.1 acceptance ("no-mock ledger
//! PreparedResume path") and frankensim-sj31i.52.5.2 close-gate criterion 3
//! (terminal-family E2E owner consumption).
//!
//! Everything here is real. The solver state family is chartered OUTSIDE
//! fs-exec, which is exactly the downstream-consumer experience the owner
//! charter exists for. The pause boundary is minted by a real
//! request -> drain -> finalize sequence through `CancelGate` and
//! `DrainTracker`. The sealed artifact is stored in and read back from a
//! real FrankenSQLite ledger. Activation happens only through
//! `prepare_resume` + `drive_v2_prepared`; a tampered byte stream refuses.
//! No mocks, no cfg(test) seams in fs-exec, no fixture standing in for
//! storage.

use fs_blake3::identity::CanonicalLimits;
use fs_exec::{Budget, CancelGate, Cx, DrainTracker, ExecMode, RunId, StreamKey};
use fs_exec::solver::snapshot_v2::{
    self, ExpectedResumeContextV2, PausedSnapshotBoundaryV2, SnapshotLimitsV2,
};
use fs_exec::solver::{
    drive_v2_prepared, prepare_resume, PreparableSolverV2, ResumableSolverV2, SolverProgress,
};
use fs_ledger::Ledger;

// ---------------------------------------------------------------------------
// Chartered probe state family (defined outside fs-exec on purpose)
// ---------------------------------------------------------------------------

/// Owner-bound identity charter for this test family. Deriving identities
/// from (owner, family, grammar, version) is what makes copied byte
/// constants inert: a foreign type copying these constants still refuses at
/// the charter gate because its own charter derives different identities.
const PROBE_CHARTER: snapshot_v2::StateIdentityCharterV2 = snapshot_v2::StateIdentityCharterV2 {
    owner: "fs-ledger::tests::prepared-resume-probe",
    state_family: "ledger-preparation-probe",
    schema_grammar: "steps: u64-le; residual: [f64; le]",
    codec_grammar: "v2 framed, u64-le count + f64-le slice",
    codec_version: 1,
};

// Pinned from PROBE_CHARTER derivation (golden-pin convention); the drift
// guard below refuses any divergence between these constants and a fresh
// derivation, so rotation of the charter text fails loudly here first.
const PROBE_STATE_TYPE_ID_V2: snapshot_v2::SnapshotStateTypeIdV2 =
    snapshot_v2::SnapshotStateTypeIdV2::from_bytes([0x00; 32]);
const PROBE_STATE_SCHEMA_ID_V2: snapshot_v2::SnapshotStateSchemaIdV2 =
    snapshot_v2::SnapshotStateSchemaIdV2::from_bytes([0x00; 32]);
const PROBE_STATE_CODEC_ID_V2: snapshot_v2::SnapshotStateCodecIdV2 =
    snapshot_v2::SnapshotStateCodecIdV2::from_bytes([0x00; 32]);
const PROBE_CODEC_VERSION_V2: u32 = 1;

/// This family consumes no stochastic work: its post-decode RNG cursor is a
/// pinned constant, and the resume context declares exactly the same value,
/// so the manifest comparison is meaningful rather than vacuous.
const PROBE_RNG_COUNTER_BYTES: [u8; 32] = [0x77; 32];
/// This family runs under an untimed budget: pinned constant on both sides.
const PROBE_BUDGET_BYTES: [u8; 32] = [0x55; 32];

#[derive(Debug, Clone, PartialEq)]
struct ProbeState {
    steps: u64,
    residual: Vec<f64>,
}

impl SolverStateV2 for ProbeState {
    const STATE_TYPE_ID_V2: snapshot_v2::SnapshotStateTypeIdV2 = PROBE_STATE_TYPE_ID_V2;
    const STATE_SCHEMA_ID_V2: snapshot_v2::SnapshotStateSchemaIdV2 = PROBE_STATE_SCHEMA_ID_V2;
    const STATE_CODEC_ID_V2: snapshot_v2::SnapshotStateCodecIdV2 = PROBE_STATE_CODEC_ID_V2;
    const STATE_CODEC_VERSION_V2: u32 = PROBE_CODEC_VERSION_V2;

    fn charter() -> Option<&'static snapshot_v2::StateIdentityCharterV2> {
        Some(&PROBE_CHARTER)
    }

    fn encode_v2(
        &self,
        encoder: &mut snapshot_v2::SnapshotEncoderV2<'_>,
    ) -> Result<(), snapshot_v2::SnapshotV2Error> {
        encoder.put_u64(self.steps)?;
        encoder.put_f64_slice(&self.residual)
    }

    fn decode_v2(
        decoder: &mut snapshot_v2::SnapshotDecoderV2<'_, '_>,
    ) -> Result<Self, snapshot_v2::SnapshotV2Error> {
        Ok(Self {
            steps: decoder.get_u64()?,
            residual: decoder.get_f64_vec()?,
        })
    }
}

// ---------------------------------------------------------------------------
// Probe solver
// ---------------------------------------------------------------------------

struct ProbeSolver;

impl ResumableSolverV2 for ProbeSolver {
    type State = ProbeState;
    type Out = Vec<f64>;

    fn step_v2(&self, state: &mut Self::State, _cx: &Cx<'_>) -> StepVerdict<Self::Out> {
        state.steps += 1;
        StepVerdict::Done(state.residual.clone())
    }
}

impl PreparableSolverV2 for ProbeSolver {
    fn expected_resume_context(&self) -> ExpectedResumeContextV2 {
        probe_context()
    }

    fn decoded_state_manifest(
        &self,
        _state: &Self::State,
    ) -> snapshot_v2::DecodedStateManifestV2 {
        snapshot_v2::DecodedStateManifestV2 {
            state_type: <ProbeState as SolverStateV2>::STATE_TYPE_ID_V2,
            state_schema: <ProbeState as SolverStateV2>::STATE_SCHEMA_ID_V2,
            state_codec: <ProbeState as SolverStateV2>::STATE_CODEC_ID_V2,
            state_codec_version: <ProbeState as SolverStateV2>::STATE_CODEC_VERSION_V2,
            rng_counter: snapshot_v2::SnapshotRngCounterIdV2::from_bytes(PROBE_RNG_COUNTER_BYTES),
            budget: snapshot_v2::SnapshotBudgetStateIdV2::from_bytes(PROBE_BUDGET_BYTES),
        }
    }
}

fn probe_limits() -> SnapshotLimitsV2 {
    SnapshotLimitsV2::new(
        1 << 20,
        64,
        CanonicalLimits::new(16_384, 4_096, 32, 32, 64),
        4_096,
        1 << 20,
        64,
    )
}

fn paused_boundary(pause_request: u8, gate_generation: u64, run: u64) -> PausedSnapshotBoundaryV2 {
    let gate = CancelGate::new();
    let tracker = DrainTracker::new(RunId(run), &gate);
    let worker = tracker.register_worker().expect("fixture worker registers");
    gate.request();
    drop(worker);
    let report = tracker.finalize().expect("run fully drained");
    PausedSnapshotBoundaryV2::from_drain_report(
        report,
        snapshot_v2::SnapshotPauseRequestIdV2::from_bytes([pause_request; 32]),
        gate_generation,
    )
}

fn probe_context() -> ExpectedResumeContextV2 {
    ExpectedResumeContextV2::for_paused_state::<ProbeState>(
        snapshot_v2::SnapshotAlgorithmIdV2::from_bytes([0x5E; 32]),
        1,
        snapshot_v2::SnapshotProblemIdV2::from_bytes([0xB0; 32]),
        snapshot_v2::SnapshotRngCounterIdV2::from_bytes(PROBE_RNG_COUNTER_BYTES),
        snapshot_v2::SnapshotDeterminismV2::Deterministic,
        snapshot_v2::SnapshotExecutionFingerprintIdV2::from_bytes([0x3F; 32]),
        snapshot_v2::SnapshotBudgetStateIdV2::from_bytes(PROBE_BUDGET_BYTES),
        snapshot_v2::SnapshotProvenanceIdV2::from_bytes([0x99; 32]),
        paused_boundary(0x66, 9, 17),
    )
}

fn hex(bytes: &[u8]) -> String {
    bytes.iter().map(|b| format!("{b:02x}")).collect()
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

/// Golden-pin drift guard: the pinned constants must equal a fresh charter
/// derivation. Failure prints both sides as hex so rotation lands here
/// first, before any envelope is produced.
#[test]
fn probe_charter_derivation_matches_pinned_constants() {
    let derived_type = PROBE_CHARTER.state_type_id();
    let derived_schema = PROBE_CHARTER.state_schema_id();
    let derived_codec = PROBE_CHARTER.state_codec_id();
    assert_eq!(
        <ProbeState as SolverStateV2>::STATE_TYPE_ID_V2,
        derived_type,
        "state-type pin drifted: derived {}",
        hex(derived_type.as_bytes())
    );
    assert_eq!(
        <ProbeState as SolverStateV2>::STATE_SCHEMA_ID_V2,
        derived_schema,
        "state-schema pin drifted: derived {}",
        hex(derived_schema.as_bytes())
    );
    assert_eq!(
        <ProbeState as SolverStateV2>::STATE_CODEC_ID_V2,
        derived_codec,
        "state-codec pin drifted: derived {}",
        hex(derived_codec.as_bytes())
    );
}

/// The terminal no-mock consumption path: seal -> ledger artifact -> read
/// back -> prepare against the actual solver instance -> activate through
/// the prepared-only drive path. A tampered byte stream refuses instead of
/// activating.
#[test]
fn prepared_resume_round_trips_through_real_ledger_and_activates_only_prepared() {
    let state = ProbeState {
        steps: 7,
        residual: vec![1.5, -2.25, 0.0, 9.75],
    };
    let context = probe_context();
    let limits = probe_limits();

    // Seal under the chartered identity.
    let sealed = state.seal_v2(&context, limits, || false).expect("v2 seal");
    let expectation = sealed.expectation();
    let bytes = sealed.bytes().to_vec();

    // Real FrankenSQLite ledger stores the artifact under its content hash.
    let ledger = Ledger::open(":memory:").expect("open ledger");
    let receipt = ledger
        .put_artifact("solver-snapshot-v2", &bytes, None)
        .expect("artifact put");
    assert!(!receipt.deduped, "first put is fresh");
    let info = ledger
        .artifact_info(&receipt.hash)
        .expect("artifact_info")
        .expect("artifact present");
    assert_eq!(info.kind, "solver-snapshot-v2");
    assert_eq!(info.len, bytes.len() as u64);

    // Read back by the caller-retained hash: the ledger returns exactly the
    // sealed bytes.
    let read_back = ledger
        .get_artifact(&receipt.hash)
        .expect("get_artifact")
        .expect("artifact readable");
    assert_eq!(read_back, bytes, "ledger round trip is byte-exact");

    // A hostile storage bit-flip refuses before any decode or preparation.
    let mut tampered = read_back.clone();
    let last = tampered.len() - 1;
    tampered[last] ^= 0x01;
    let tamper_error =
        ProbeState::unseal_v2_expected(&tampered, &expectation, limits, || false)
            .expect_err("tampered stream must refuse");
    let _ = tamper_error;

    // The exact ledger-retained stream opens and prepares against the
    // ACTUAL solver instance.
    let opened =
        ProbeState::unseal_v2_expected(&read_back, &expectation, limits, || false)
            .expect("exact retained roots authorize decoding");
    let prepared = prepare_resume(&ProbeSolver, opened, || false).expect("prepares");

    // Activation through the prepared-only drive path completes and matches
    // a directly-driven run bit-exactly, under one real execution context.
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let gate = CancelGate::new();
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 1,
                kernel_id: 7,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        let mut direct = state.clone();
        let direct_out = match ProbeSolver.step_v2(&mut direct, &cx) {
            StepVerdict::Done(out) => out,
            StepVerdict::Continue => panic!("probe finishes in one step"),
        };
        match drive_v2_prepared(&ProbeSolver, prepared, &cx) {
            SolverProgress::Done(out) => {
                assert_eq!(out, direct_out, "ledger-resumed run matches direct run");
            }
            SolverProgress::Paused(_) => panic!("probe completes in one step"),
        }
    });
}
