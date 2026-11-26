//! Unit proof for `crate::stage_driver` via a non-thermal 3-stage pipeline.

use fs_blake3::identity::ContentId;
use fs_blake3::{ContentHash, hash_bytes};
use fs_exec::CancelGate;
use fs_exec::solver::{
    LegacySnapshotExpectationV1, LegacySnapshotLimitsV1, LegacySnapshotV1Adapter,
    LegacySolverStateV1, codec,
};
use fs_ledger::{EdgeRole, Ledger, LedgerError, OpOutcome};

use super::*;
use crate::solve::{
    EvidenceWork, InvocationWorkLedger, SolveDriverState, SolveOutcome, SolveRefusal, SolveRunId,
    SolveRunStatus, resume_identity,
};

const DEMO_TYPE_ID: u64 = 0x6673_636c_6964_656d; // b"fsclidem"
const DEMO_SCHEMA: PipelineSchema = PipelineSchema {
    driver_version: 1,
    stage_ir_schema: "frankensim.cli.demo-stage.v1",
    run_receipt_schema: "frankensim.cli.demo-run-receipt.v1",
    progress_schema: "frankensim.cli.demo-progress.v1",
    refusal_schema: "frankensim.cli.demo-refusal.v1",
    source_artifact_kind: "demo-source",
    stage_state_kind: "demo-stage-state",
    stage_receipt_kind: "demo-stage-receipt",
    run_receipt_kind: "demo-run-receipt",
    refusal_artifact_kind: "demo-refusal",
    cores: 1,
    terminal_verb: "demo.terminal",
    ledger_scope_prefix: "demo",
    session_open_prefix: "demo-open",
    meter_prefix: "demo-meter",
    resume_command_prefix: "frankensim demo --resume",
    budget_resume_hint: "raise budget and re-run",
    incompatible_run_hint: "re-run with matching demo driver version",
    unknown_run_hint: "pass a known demo run id",
    run_noun: "demo run",
    terminal_fallback_stage: "package",
};
const DEMO_SCHEMA_V2: PipelineSchema = PipelineSchema {
    driver_version: 2,
    ..DEMO_SCHEMA
};

struct DemoState(SolveDriverState);
impl LegacySolverStateV1 for DemoState {
    const SCHEMA_VERSION_V1: u32 = 1;
    const TYPE_ID_V1: u64 = DEMO_TYPE_ID;
    fn encode_v1(&self, enc: &mut codec::Enc) {
        self.0.encode_v1(enc);
    }
    fn decode_v1(dec: &mut codec::Dec<'_>) -> Result<Self, codec::CodecError> {
        SolveDriverState::decode_v1(dec).map(Self)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum DemoStage {
    Parse,
    Compute,
    Package,
}

impl StageDescriptor for DemoStage {
    fn all() -> &'static [Self] {
        &[Self::Parse, Self::Compute, Self::Package]
    }
    fn ordinal(self) -> u32 {
        self as u32
    }
    fn name(self) -> &'static str {
        ["parse", "compute", "package"][self as usize]
    }
    fn verb(self) -> &'static str {
        ["demo.parse", "demo.compute", "demo.package"][self as usize]
    }
    fn gap_dependency(self) -> Option<&'static str> {
        None
    }
}

#[derive(Debug, Default)]
struct DemoCtx {
    parsed: Option<u64>,
    computed: Option<(u64, ContentHash)>,
}

struct DemoPipe<'a> {
    src: &'a [u8],
    gate: CancelGate,
    cancel_after: Option<DemoStage>,
    cancel_pre_pub: Option<DemoStage>,
    gap_on: Option<DemoStage>,
}

impl<'a> DemoPipe<'a> {
    fn new(src: &'a [u8]) -> Self {
        Self {
            src,
            gate: CancelGate::new_clock_free(),
            cancel_after: None,
            cancel_pre_pub: None,
            gap_on: None,
        }
    }
}

impl StagePipeline for DemoPipe<'_> {
    type Stage = DemoStage;
    type Context = DemoCtx;
    type StageCommit = (Option<u64>, Option<(u64, ContentHash)>);

    fn schema(&self) -> &'static PipelineSchema {
        &DEMO_SCHEMA
    }
    fn canonical_source(&self) -> &[u8] {
        self.src
    }
    fn seal_state(&self, run: SolveRunId, state: &SolveDriverState) -> Vec<u8> {
        let env = DemoState(state.clone());
        LegacySnapshotV1Adapter::<DemoState>::seal(&env, envelope_provenance(run))
    }
    fn declaration_gap(&self, stage: Self::Stage) -> Option<(&'static str, String, String)> {
        (self.gap_on == Some(stage)).then(|| {
            let msg = format!("stage `{}` input undeclared", stage.name());
            ("cli-demo-undeclared", msg, "declare demo input".to_string())
        })
    }
    fn is_cancelled(&self) -> bool {
        self.gate.is_requested()
    }
    fn pre_publish_checkpoint(&self, stage: Self::Stage) -> Result<(), ()> {
        if self.cancel_pre_pub == Some(stage) {
            self.gate.request();
        }
        (!self.gate.is_requested()).then_some(()).ok_or(())
    }
    fn execute_stage(
        &mut self,
        stage: Self::Stage,
        _: &StageSessionView<'_>,
        _: &SolveDriverState,
        ctx: &mut Self::Context,
    ) -> Result<StageExecution<Self::StageCommit>, SolveRefusal> {
        let (receipt_json, side_artifacts, on_commit) = match stage {
            DemoStage::Parse => {
                let s = std::str::from_utf8(self.src).unwrap().trim();
                let n: u64 = s.parse().unwrap();
                (
                    format!("{{\"stage\":\"parse\",\"n\":{n}}}"),
                    Vec::new(),
                    (Some(n), None),
                )
            }
            DemoStage::Compute => {
                let sq = ctx.parsed.unwrap().pow(2);
                let bytes = format!("value={sq}").into_bytes();
                let artifact = hash_bytes(&bytes);
                let json = format!(
                    "{{\"stage\":\"compute\",\"sq\":{sq},\"out\":\"{}\"}}",
                    artifact.to_hex()
                );
                let side = RetainedSideArtifact {
                    kind: "demo-computed",
                    artifact,
                    bytes,
                };
                (json, vec![side], (None, Some((sq, artifact))))
            }
            DemoStage::Package => {
                let (sq, art) = ctx.computed.unwrap();
                let json = format!(
                    "{{\"stage\":\"package\",\"sq\":{sq},\"side\":\"{}\"}}",
                    art.to_hex()
                );
                (json, Vec::new(), (None, None))
            }
        };
        Ok(StageExecution {
            receipt_json,
            side_artifacts,
            on_commit,
        })
    }
    fn commit_stage(&mut self, ctx: &mut Self::Context, commit: Self::StageCommit) {
        if let Some(n) = commit.0 {
            ctx.parsed = Some(n);
            if self.cancel_after == Some(DemoStage::Parse) {
                self.gate.request();
            }
        }
        if let Some(comp) = commit.1 {
            ctx.computed = Some(comp);
            if self.cancel_after == Some(DemoStage::Compute) {
                self.gate.request();
            }
        }
    }
    fn link_stage_in_edges(
        &self,
        ledger: &Ledger,
        op: i64,
        stage: Self::Stage,
        before: &SolveDriverState,
        ctx: &Self::Context,
    ) -> Result<(), LedgerError> {
        if stage == DemoStage::Package {
            ledger.link(op, &before.completed[1].receipt, EdgeRole::In)?;
            ledger.link(op, &ctx.computed.unwrap().1, EdgeRole::In)?;
        }
        Ok(())
    }
}

struct DemoVerified {
    state: SolveDriverState,
    artifact: ContentHash,
    ctx: DemoCtx,
    src: Vec<u8>,
}
struct DemoAttestor(&'static PipelineSchema);

fn find_edge(ledger: &Ledger, op: i64, role: EdgeRole, kind: &str) -> Option<ContentHash> {
    let edges = ledger.op_artifact_edges_bounded(op, 64).unwrap().edges;
    edges.into_iter().find_map(|e| {
        let info = ledger.artifact_info(&e.artifact).unwrap();
        (e.role == role && info.is_some_and(|i| i.kind == kind)).then_some(e.artifact)
    })
}

impl ResumeAttestor for DemoAttestor {
    type Stage = DemoStage;
    type Verified = DemoVerified;

    fn schema(&self) -> &'static PipelineSchema {
        self.0
    }
    fn decode_state(
        &self,
        ledger: &Ledger,
        run: SolveRunId,
        artifact: ContentHash,
        _: EvidenceWork<'_>,
    ) -> Result<SolveDriverState, SolveRefusal> {
        let missing = || resume_identity("missing");
        let bytes = ledger
            .get_artifact(&artifact)
            .unwrap()
            .ok_or_else(missing)?;
        let cid = ContentId::parse_slice(artifact.as_bytes()).unwrap();
        let exp = LegacySnapshotExpectationV1::new(cid, DEMO_TYPE_ID, 1, envelope_provenance(run));
        let lim = LegacySnapshotLimitsV1::new(1024 * 1024, 64 * 1024);
        LegacySnapshotV1Adapter::<DemoState>::open_expected(&bytes, exp, lim, || false)
            .map(|o| o.into_parts().0.0)
            .map_err(|e| resume_identity(format!("bad demo envelope: {e:?}")))
    }
    fn validate_candidate(
        &mut self,
        ledger: &Ledger,
        run: SolveRunId,
        mut state: SolveDriverState,
        artifact: ContentHash,
        discovery_op: i64,
        work: EvidenceWork<'_>,
    ) -> Result<Self::Verified, SolveRefusal> {
        validate_state_shape_for::<DemoStage>(&state, run, self.0.cores)?;
        let sk = self.0.source_artifact_kind;
        let src_h = find_edge(ledger, state.completed[0].op_id, EdgeRole::In, sk)
            .ok_or_else(|| resume_identity("missing source"))?;
        let src = ledger.get_artifact(&src_h).unwrap().unwrap();
        let mut ctx = DemoCtx::default();
        let mut prev = None;
        let mut resolved_consumption = None;
        for (idx, done) in state.completed.iter().enumerate() {
            let edges = ledger
                .op_artifact_edges_bounded(done.op_id, 64)
                .unwrap()
                .edges;
            assert_eq!(
                ledger.op_artifact_edge_seal(done.op_id).unwrap(),
                Some(edges.len())
            );
            let cp_h =
                find_edge(ledger, done.op_id, EdgeRole::Out, self.0.stage_state_kind).unwrap();
            let cp = self.decode_state(ledger, run, cp_h, work)?;
            validate_checkpoint_prefix_for(&cp, &state, idx, prev.as_ref(), self.0.cores)?;
            let stage_consumption = resolve_stage_consumption(
                ledger,
                done.op_id,
                done.ordinal,
                &cp,
                resolved_consumption,
                self.0.cores,
            )?;
            if idx == 0 {
                ctx.parsed = Some(std::str::from_utf8(&src).unwrap().trim().parse().unwrap());
            } else if idx == 1 {
                let sq = ctx.parsed.unwrap().pow(2);
                ctx.computed = Some((sq, hash_bytes(format!("value={sq}").as_bytes())));
            }
            if idx + 1 == state.completed.len() {
                assert_eq!((discovery_op, cp_h), (done.op_id, artifact));
            }
            prev = Some(cp);
            resolved_consumption = Some(stage_consumption);
        }
        if let Some((consumed_core_s, consumed_wall_s)) = resolved_consumption {
            state.consumed_core_s = consumed_core_s;
            state.consumed_wall_s = consumed_wall_s;
        }
        Ok(DemoVerified {
            state,
            artifact,
            ctx,
            src,
        })
    }
    fn completed_len(v: &Self::Verified) -> usize {
        v.state.completed.len()
    }
}

fn demo_id(src: &[u8]) -> (ContentHash, SolveRunId) {
    let h = hash_bytes(src);
    let r = derive_run_id_bytes(
        "org.frankensim.fs-cli.demo.v1",
        h,
        "2026-07",
        "t",
        7,
        hash_bytes(b""),
        1,
    );
    (h, SolveRunId::from_content_hash(r))
}

fn open_sess<'a>(
    ledger: &'a Ledger,
    clock: &'a mut dyn FnMut() -> f64,
    run: SolveRunId,
    h: ContentHash,
    wall_s: f64,
    strip: Option<&str>,
) -> DurableStageSession<'a> {
    let mut tok =
        derive_pipeline_capability_token::<DemoStage>(run, wall_s, 1 << 20, &DEMO_SCHEMA).unwrap();
    if let Some(v) = strip {
        tok.ops.retain(|op| op != v);
    }
    let exp = SessionExplicits {
        seed: 7u64.to_le_bytes(),
        versions_json: "{\"constellation\":\"2026-07\"}".into(),
        budget_json: "{\"wall_s\":10.0}".into(),
        capability_json: "{\"ops\":[\"demo.*\"]}".into(),
    };
    DurableStageSession::open(ledger, clock, run, h, &DEMO_SCHEMA, tok, exp).unwrap()
}

fn run_demo(
    ledger: &Ledger,
    src: &[u8],
    wall_s: f64,
    clock: &mut dyn FnMut() -> f64,
    mut pipe: DemoPipe<'_>,
    strip: Option<&str>,
) -> (Result<SolveOutcome, SolveRefusal>, Vec<String>) {
    let (h, run) = demo_id(src);
    let mut sess = open_sess(ledger, clock, run, h, wall_s, strip);
    let st = SolveDriverState {
        run: *run.as_bytes(),
        project: *h.as_bytes(),
        consumed_core_s: 0.0,
        consumed_wall_s: 0.0,
        completed: Vec::new(),
    };
    let res = sess.drive(&mut pipe, st, DemoCtx::default(), None);
    let mut prog = Vec::new();
    sess.drain_progress(&mut prog);
    (res, prog)
}

fn attest(
    ledger: &Ledger,
    run: SolveRunId,
    s: &'static PipelineSchema,
) -> Result<DemoVerified, SolveRefusal> {
    let gate = CancelGate::new_clock_free();
    let inv = InvocationWorkLedger::default();
    let work = EvidenceWork::new(&gate, None, &inv);
    discover_and_attest_resume(ledger, run, work, &mut DemoAttestor(s))
}

fn resume_demo(
    ledger: &Ledger,
    run: SolveRunId,
    wall_s: f64,
    clk: &mut dyn FnMut() -> f64,
    cancel_after: Option<DemoStage>,
) -> Result<SolveOutcome, SolveRefusal> {
    let v = attest(ledger, run, &DEMO_SCHEMA)?;
    let mut sess = open_sess(ledger, clk, run, hash_bytes(&v.src), wall_s, None);
    sess.restore_budget(&v.state, v.artifact)?;
    let mut pipe = DemoPipe {
        cancel_after,
        ..DemoPipe::new(&v.src)
    };
    sess.drive(&mut pipe, v.state, v.ctx, Some(v.artifact))
}

fn step_clock(step: f64) -> impl FnMut() -> f64 {
    let mut t = 0.0;
    move || {
        t += step;
        t
    }
}

fn seq_clock(times: &'static [f64]) -> impl FnMut() -> f64 {
    let mut i = 0usize;
    move || {
        let v = times[i];
        i += 1;
        v
    }
}

fn mem() -> Ledger {
    Ledger::open(":memory:").unwrap()
}

fn err_demo(
    l: &Ledger,
    src: &[u8],
    clk: &mut dyn FnMut() -> f64,
    pipe: DemoPipe<'_>,
    strip: Option<&str>,
) -> SolveRefusal {
    run_demo(l, src, 10.0, clk, pipe, strip)
        .0
        .expect_err("refusal")
}

fn assert_err_op(l: &Ledger, op: Option<i64>) {
    let row = l.op(op.unwrap()).unwrap().unwrap();
    assert_eq!(row.outcome.as_deref(), Some(OpOutcome::Error.as_str()));
}

#[test]
fn non_thermal_pipeline_runs_and_resumes_after_kill_at_every_stage() {
    let (src, lf) = (b"12".as_slice(), mem());
    let (_, run) = demo_id(src);
    let (fresh, prog) = run_demo(
        &lf,
        src,
        10.0,
        &mut step_clock(0.5),
        DemoPipe::new(src),
        None,
    );
    let fresh = fresh.expect("fresh");
    assert_eq!(fresh.status, SolveRunStatus::Completed);
    assert_eq!(
        (fresh.stages.len(), fresh.prior_stages, prog.len()),
        (3, 0, 3)
    );
    assert!(fresh.run_receipt.is_some());
    assert!(prog[0].contains("\"stage\":\"parse\"") && prog[2].contains("\"stage\":\"package\""));

    // Kill after stage 0, then stage 1, then finish stage 2
    let (ledger, mut clk) = (mem(), step_clock(0.5));
    let mut p0 = DemoPipe::new(src);
    p0.cancel_after = Some(DemoStage::Parse);
    let e0 = err_demo(&ledger, src, &mut clk, p0, None);
    assert_eq!(
        (e0.code, e0.stage),
        ("cli-solve-cancelled", Some("compute"))
    );
    assert!(e0.fix.contains("frankensim demo --resume"));

    let e1 = resume_demo(&ledger, run, 10.0, &mut clk, Some(DemoStage::Compute)).unwrap_err();
    assert_eq!(
        (e1.code, e1.stage),
        ("cli-solve-cancelled", Some("package"))
    );
    let done = resume_demo(&ledger, run, 10.0, &mut clk, None).expect("resume finishes");
    assert_eq!(
        (done.status, done.prior_stages, done.stages.len()),
        (SolveRunStatus::Completed, 2, 1)
    );
    assert_eq!(done.stages[0].receipt, fresh.stages[2].receipt);

    // Pre-publish cancellation leaves stage uncommitted and resumable
    let (lp, mut cp) = (mem(), step_clock(0.5));
    let mut pp = DemoPipe::new(src);
    pp.cancel_pre_pub = Some(DemoStage::Compute);
    assert_eq!(
        err_demo(&lp, src, &mut cp, pp, None).code,
        "cli-solve-cancelled"
    );
    let dp = resume_demo(&lp, run, 10.0, &mut cp, None).expect("resume pre");
    assert_eq!((dp.prior_stages, dp.stages.len()), (1, 2));
    assert_eq!(dp.stages[1].receipt, fresh.stages[2].receipt);
}

#[test]
fn budget_warnings_stop_and_resume_budget_refusal_are_enforced() {
    let (ledger, src) = (mem(), b"9".as_slice());
    let (_, run) = demo_id(src);
    let mut clk = seq_clock(&[0.0, 0.6, 0.6, 1.1]);
    let (res, prog) = run_demo(&ledger, src, 1.0, &mut clk, DemoPipe::new(src), None);
    let out = res.expect("honest partial stop");
    assert!(matches!(out.status, SolveRunStatus::BudgetExceeded { .. }) && out.stages.len() == 2);
    assert!(
        prog.iter()
            .any(|l| l.contains("\"status\":\"warning\"") && l.contains("\"fraction\":0.5"))
    );
    let err = resume_demo(&ledger, run, 1.0, &mut step_clock(0.1), None).unwrap_err();
    assert_eq!(err.code, "cli-solve-resume-budget");
}

#[test]
fn capability_declaration_gap_and_discovery_ambiguity_or_version_refuse() {
    let src = b"5";
    let (h, run) = demo_id(src);
    let mut clk = step_clock(0.1);

    let lc = mem();
    let ec = err_demo(&lc, src, &mut clk, DemoPipe::new(src), Some("demo.compute"));
    assert_eq!(
        (ec.code, ec.stage),
        ("cli-solve-capability", Some("compute"))
    );
    assert_err_op(&lc, ec.recorded_op);

    let lg = mem();
    let mut pg = DemoPipe::new(src);
    pg.gap_on = Some(DemoStage::Parse);
    let eg = err_demo(&lg, src, &mut clk, pg, None);
    assert_eq!(eg.code, "cli-demo-undeclared");
    assert!(eg.recorded_op.is_some());

    let ld = mem();
    for _ in 0..2 {
        let mut p = DemoPipe::new(src);
        p.cancel_after = Some(DemoStage::Parse);
        let _ = run_demo(&ld, src, 10.0, &mut clk, p, None);
    }
    let bad_v = attest(&ld, run, &DEMO_SCHEMA_V2).err().unwrap();
    assert_eq!(bad_v.code, "cli-solve-run-incompatible");
    let unk = SolveRunId::from_content_hash(hash_bytes(b"unknown"));
    assert_eq!(
        attest(&ld, unk, &DEMO_SCHEMA).err().unwrap().code,
        "cli-solve-unknown-run"
    );
    let amb = attest(&ld, run, &DEMO_SCHEMA).err().unwrap();
    assert_eq!(amb.code, "cli-solve-resume-identity");
    assert!(
        amb.what
            .contains("competing independently valid 1-stage checkpoints")
    );
    assert!(
        format_stage_ir(&DEMO_SCHEMA, DemoStage::Parse, run, h).contains("\"driver_version\":1")
    );

    let lk = mem();
    let mut bad_clk = seq_clock(&[2.0, 1.0]);
    let ek = err_demo(&lk, src, &mut bad_clk, DemoPipe::new(src), None);
    assert_eq!(ek.code, "cli-solve-budget");
    assert_err_op(&lk, ek.recorded_op);
}

fn copy_stage(
    src: &Ledger,
    dst: &Ledger,
    op_id: i64,
    ord: u32,
    metrics: Option<(f64, f64)>,
    cp_override: Option<Vec<u8>>,
) {
    let row = src.op(op_id).unwrap().unwrap();
    let fx = FiveExplicits {
        seed: &row.seed,
        versions: &row.versions,
        budget: &row.budget,
        capability: &row.capability,
    };
    let op = dst
        .begin_op(row.session.as_deref(), &row.ir, &fx, row.t_start)
        .unwrap();
    assert_eq!(op, op_id);
    let edges = src.op_artifact_edges_bounded(op_id, 64).unwrap().edges;
    for edge in &edges {
        let kind = src.artifact_info(&edge.artifact).unwrap().unwrap().kind;
        let bytes = match (&cp_override, kind == DEMO_SCHEMA.stage_state_kind) {
            (Some(ov), true) => ov.clone(),
            _ => src.get_artifact(&edge.artifact).unwrap().unwrap(),
        };
        let put = dst.put_artifact(&kind, &bytes, None).unwrap();
        dst.link(op, &put.hash, edge.role).unwrap();
    }
    dst.seal_op_artifact_edges(op, edges.len()).unwrap();
    let t_end = stage_end_t(ord);
    if let Some((c, w)) = metrics {
        dst.record_metric(op, t_end, DRIVER_CONSUMED_CORE_S_METRIC, c)
            .unwrap();
        dst.record_metric(op, t_end, DRIVER_CONSUMED_WALL_S_METRIC, w)
            .unwrap();
    }
    dst.finish_op(op, OpOutcome::Ok, None, t_end).unwrap();
}

#[test]
fn consumption_metrics_replay_clean_and_tamper_or_legacy_enforced() {
    let src = b"11".as_slice();
    let (src_h, run) = demo_id(src);
    let (la, lb) = (mem(), mem());
    let (oa, _) = run_demo(
        &la,
        src,
        10.0,
        &mut step_clock(0.25),
        DemoPipe::new(src),
        None,
    );
    let (ob, _) = run_demo(
        &lb,
        src,
        10.0,
        &mut step_clock(1.75),
        DemoPipe::new(src),
        None,
    );
    assert_eq!(oa.unwrap().run_receipt, ob.unwrap().run_receipt);
    let v = la.replay_verdict(MAIN_BRANCH, &lb, MAIN_BRANCH).unwrap();
    assert!(v.is_replay_clean() && v.compared == 4, "{v:?}");

    // Tamper 1: missing consumption metrics on deterministic 0.0 checkpoint refuses.
    let l1 = mem();
    copy_stage(&la, &l1, 1, 0, None, None);
    let e1 = attest(&l1, run, &DEMO_SCHEMA).err().unwrap();
    assert!(e1.code == "cli-solve-resume-identity" && e1.what.contains("lacks retained budget"));

    // Tamper 2: regressing consumption metrics across stages refuses.
    let l2 = mem();
    copy_stage(&la, &l2, 1, 0, Some((0.5, 0.5)), None);
    copy_stage(&la, &l2, 2, 1, Some((0.25, 0.25)), None);
    let e2 = attest(&l2, run, &DEMO_SCHEMA).err().unwrap();
    assert!(e2.code == "cli-solve-resume-identity" && e2.what.contains("regresses"));

    // Tamper 3: single-core core_s != wall_s mismatch refuses.
    let l3 = mem();
    copy_stage(&la, &l3, 1, 0, Some((0.5, 0.4)), None);
    assert_eq!(
        attest(&l3, run, &DEMO_SCHEMA).err().unwrap().code,
        "cli-solve-resume-identity"
    );

    // Legacy v1 checkpoint with positive sealed consumption and no metrics restores.
    let l4 = mem();
    let rc_h = find_edge(&la, 1, EdgeRole::Out, DEMO_SCHEMA.stage_receipt_kind).unwrap();
    let legacy = SolveDriverState {
        run: *run.as_bytes(),
        project: *src_h.as_bytes(),
        consumed_core_s: 0.75,
        consumed_wall_s: 0.75,
        completed: vec![CompletedStage {
            ordinal: 0,
            op_id: 1,
            receipt: rc_h,
        }],
    };
    let sealed =
        LegacySnapshotV1Adapter::<DemoState>::seal(&DemoState(legacy), envelope_provenance(run));
    copy_stage(&la, &l4, 1, 0, None, Some(sealed));
    let ok = attest(&l4, run, &DEMO_SCHEMA).expect("legacy restores");
    assert!(
        (ok.state.consumed_wall_s - 0.75).abs() < 1e-12
            && (ok.state.consumed_core_s - 0.75).abs() < 1e-12
    );
}
