//! fs-detaudit real-subject acceptance (bead 6nb.6): the audit runs
//! against a REAL fs-exec non-associative pooled reduction, a REAL
//! fs-rand logical stream, and the REAL PV/OED demo study, across the
//! full host worker matrix.

use core::ops::ControlFlow;
use fs_detaudit::{AuditConfig, StageHash, StagedTrace, Subject, WorkerMatrix, audit, fnv1a64};
use fs_exec::{
    Budget, CancelGate, Cancelled, Cx, ExecMode, PoolConfig, StreamKey, TileKernel, TilePlan,
    TilePool,
};
use fs_oed_e2e::{demo_candidates, run_campaign};
use fs_substrate::affinity::CcdTopology;

fn stage(label: &str, hash: u64) -> StageHash {
    StageHash {
        label: label.to_owned(),
        hash,
    }
}

/// Non-associative float reduction over the deterministic tile pool —
/// the fs-exec acceptance subject.
struct FloatKernel {
    tiles: u64,
}

impl TileKernel for FloatKernel {
    type Out = f64;

    fn tiles(&self) -> TilePlan {
        TilePlan::new("detaudit/float", self.tiles)
    }

    fn run(&self, tile: u64, cx: &Cx<'_>) -> ControlFlow<Cancelled, f64> {
        if cx.checkpoint().is_err() {
            return ControlFlow::Break(Cancelled);
        }
        let xs = cx
            .arena()
            .alloc_slice_with(fs_alloc::Site::named("detaudit/float"), 33, |i| {
                1.0 / ((tile * 33 + i as u64 + 1) as f64)
            })
            .expect("arena");
        ControlFlow::Continue(xs.iter().sum())
    }
}

fn exec_reduction_subject() -> Subject {
    Subject {
        name: "fs-exec/pooled-float-reduction",
        run: Box::new(|workers| {
            let pool = TilePool::new(PoolConfig::new(
                workers,
                CcdTopology::APPLE_M_CLASS,
                0xD3_7A_0D17,
            ));
            let out = pool
                .run(&FloatKernel { tiles: 257 })
                .expect("deterministic pool run");
            StagedTrace {
                stages: vec![
                    stage("plan", fnv1a64(b"detaudit/float:257")),
                    stage("reduction", fnv1a64(&out.to_bits().to_le_bytes())),
                ],
            }
        }),
    }
}

fn rand_stream_subject() -> Subject {
    Subject {
        name: "fs-rand/logical-stream",
        run: Box::new(|_workers| {
            // The logical stream is keyed by (seed, kernel, tile), never
            // worker-owned: any worker count reproduces the exact draws.
            let mut bytes = Vec::with_capacity(256 * 8);
            let mut stream = fs_rand::StreamKey {
                seed: 0x5EED_CA5E,
                kernel: 7,
                tile: 3,
            }
            .stream();
            for _ in 0..256 {
                bytes.extend_from_slice(&stream.next_u64().to_le_bytes());
            }
            StagedTrace {
                stages: vec![
                    stage("key", fnv1a64(b"philox:0x5EEDCA5E:7:3")),
                    stage("draws", fnv1a64(&bytes)),
                ],
            }
        }),
    }
}

#[test]
fn fs_exec_reduction_audits_bit_identical_across_the_host_matrix() {
    let report = audit(&exec_reduction_subject(), &AuditConfig::host_default(3));
    for line in report.json_lines() {
        println!("{line}");
    }
    assert!(
        report.identical(),
        "fs-exec deterministic pool reduction must be bit-identical across \
         the worker matrix: {:?}",
        report.divergences
    );
}

#[test]
fn fs_rand_stream_audits_bit_identical_across_the_host_matrix() {
    let report = audit(&rand_stream_subject(), &AuditConfig::host_default(2));
    for line in report.json_lines() {
        println!("{line}");
    }
    assert!(report.identical(), "logical streams are worker-independent");
}

/// The PV/OED demo study as an audit subject. The campaign is logically
/// scheduled under Cx (streams keyed by identity, never by worker), so
/// the worker knob is a no-op by construction — the audit verifies
/// run-to-run replay determinism of the full study artifact.
fn pv_study_subject() -> Subject {
    Subject {
        name: "fs-oed-e2e/pv-demo-study",
        run: Box::new(|_workers| {
            let stream = StreamKey {
                seed: 0x6f65_642d_6532_6501,
                kernel_id: 1,
                tile: 0,
                iteration: 0,
            };
            let gate = CancelGate::new();
            let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
            let report = pool.scope(|arena| {
                let clock = fs_exec::VirtualClock::new();
                let cx = Cx::new(
                    &gate,
                    arena,
                    stream,
                    Budget::INFINITE,
                    ExecMode::Deterministic,
                )
                .with_time_source(&clock);
                let candidates = demo_candidates().expect("demo candidates admit");
                let threshold =
                    fs_oed_e2e::ObjectiveValue::dimensionless(0.02).expect("threshold admits");
                run_campaign(&candidates, threshold, 6, &cx).expect("demo campaign runs")
            });
            // Exact-bits projection of the study artifact.
            let mut bytes = Vec::new();
            for p in report.placements() {
                bytes.extend_from_slice(p.as_bytes());
                bytes.push(0);
            }
            bytes.extend_from_slice(report.chosen_design().as_bytes());
            bytes.push(0);
            for (name, weight) in report.allocation() {
                bytes.extend_from_slice(name.as_bytes());
                bytes.extend_from_slice(&weight.to_bits().to_le_bytes());
            }
            let evpi_bytes: Vec<u8> = report
                .evpi_trace()
                .flat_map(|v| v.value().to_bits().to_le_bytes())
                .collect();
            StagedTrace {
                stages: vec![
                    stage("placements", fnv1a64(&bytes)),
                    stage("evpi-trace", fnv1a64(&evpi_bytes)),
                ],
            }
        }),
    }
}

/// Per-ISA ledger emitter for the cross-ISA report: run on EACH host
/// (`cargo test -p fs-detaudit --test real_subjects emit_isa_ledger --
/// --ignored --nocapture`), collect the `detaudit_ledger` JSON lines into
/// one file per ISA, then classify with the `cross-isa-report` bin.
#[test]
#[ignore = "ISA-ledger emitter: run per host to feed the cross-ISA report (bead 6nb.6)"]
fn emit_isa_ledger() {
    let emit = |artifact: &str, hash: u64, value_bits: Option<u64>| {
        let value = value_bits.map_or_else(|| "null".to_owned(), |v| format!("\"{v:016x}\""));
        println!(
            "{{\"detaudit_ledger\":\"{artifact}\",\"hash\":\"{hash:016x}\",\"value_bits\":{value}}}"
        );
    };
    let exec_trace = (exec_reduction_subject().run)(4);
    emit(
        "exec/pooled-float-reduction-257",
        exec_trace.stages[1].hash,
        None,
    );
    let rand_trace = (rand_stream_subject().run)(1);
    emit("rand/philox-256-draws", rand_trace.stages[1].hash, None);
    let pv_trace = (pv_study_subject().run)(1);
    emit("oed/pv-demo-placements", pv_trace.stages[0].hash, None);
    emit("oed/pv-demo-evpi-trace", pv_trace.stages[1].hash, None);
    // Platform-libm probes: the genuine cross-ISA divergence candidates,
    // declared under the libm-ULP envelope in the report policy.
    let atan2 = 0.7f64.atan2(-1.3);
    emit(
        "libm/platform-atan2",
        fnv1a64(&atan2.to_bits().to_le_bytes()),
        Some(atan2.to_bits()),
    );
    let tan = 1.234f64.tan();
    emit(
        "libm/platform-tan",
        fnv1a64(&tan.to_bits().to_le_bytes()),
        Some(tan.to_bits()),
    );
}

#[test]
fn pv_demo_study_audits_replay_deterministic() {
    let report = audit(
        &pv_study_subject(),
        &AuditConfig {
            matrix: WorkerMatrix::explicit(vec![1, 4]),
            repeats: 3,
        },
    );
    for line in report.json_lines() {
        println!("{line}");
    }
    assert!(
        report.identical(),
        "the PV demo study must replay bit-identically: {:?}",
        report.divergences
    );
}
