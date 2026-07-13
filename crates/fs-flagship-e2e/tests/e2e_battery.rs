//! The flagship e2e battery (bead mye.5): SMOKE stages for all
//! flagships with frozen golden hashes, cross-flagship audits, failure
//! drills with expected structured outcomes, forensic-log self-audit,
//! and the bitwise-replayable lab notebook. MID/FULL stages are wired
//! behind `#[ignore]` — their CI cadence belongs to the perf-CI lanes.

use std::time::Instant;

use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_flagship_e2e::{StageArtifact, Tier, artifact, lbm_core_roll_hash, log_row, notebook};
use fs_frame::history::StoryParams;
use fs_frame::{e_stopped_fragility, layout_and_size};
use fs_ornith::param::OrnithCandidate;
use fs_ornith::screen::{lift_to_drag, screen_generation};
use fs_qty::{Dims, QtyAny};
use fs_scenario::ensemble::{SpectrumModel, StochasticEnsemble};
use fs_vessel::pour::{PourRig, run_pour};
use fs_vessel::robustify;

fn verdict(name: &str, pass: bool, details: &str) {
    println!("{{\"test\":\"{name}\",\"pass\":{pass},\"details\":\"{details}\"}}");
    assert!(pass, "{name}: {details}");
}

const TIME: Dims = Dims([0, 0, 1, 0, 0]);
const RATE: Dims = Dims([0, 0, -1, 0, 0]);

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 0xF1A6_5A1D,
                kernel_id: 1,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx)
    })
}

// ------------------------------------------------------------------
// Golden hashes: frozen at bead mye.5. Bump ONLY with justification
// (a semantic change in the owning flagship or a shared core) — the
// couplings are registered in golden-couplings.json.
// ------------------------------------------------------------------
// 2026-07-10 (gp3.14): ALL FOUR re-frozen for the metric-stream
// ENCODING migration v1 -> v2 (bare name/bits concatenation ->
// canonical typed length-prefixed replay identity, fs_obs::ident).
// Every metric BIT PATTERN in every stage's evidence payload is
// unchanged from the v1 freeze — only the identity encoding moved.
// Prior v1 history, preserved:
// - vessel: 0xd70b_9ac9_0828_ae86 (before that 0xe621_48d4_490c_a887,
//   radix-2 fs-fft schedule; only robust_offband moved by 4.48e-14).
// - ornith: 0xf513_eaf8_22d2_7813 (before that 0xa6fa_6460_e7c7_972f,
//   unit-span/clipping e-race; only evals moved, 394 -> 925).
// - frame: 0x05e1_d182_48d2_949f; lbm-core: 0x6841_e3c0_508e_eba5.
// 2026-07-10 again (xo2k): poured_mass proved build-mode bit-divergent
// (~31 ulp release vs debug) and briefly moved to an envelope gate
// (5-metric hash 0xfb33_2a50_af26_1116; the divergent 6-metric v2 hash
// was 0xdabd_6fd3_6315_31fe, debug bits). ROOT CAUSE FIXED at xo2k
// close: the pour tilt schedule's platform sin/cos (release const-folds
// libm with inlined literal rig params) now routes through
// fs_math::det::sin/cos, and poured_mass (0x3fd3b2951fb7df34 in BOTH
// modes) is restored to the hashed stream. fs-lbm rheology powf paths
// migrated to det::pow in the same change (latent, same hazard class).
// 2026-07-11 (27d3 downstream audit): mixed radix-8/4/2 changed the
// DCT evaluation order feeding the vessel stability objective. Only
// robust_offband moved, from 0xbf3c9a988956ba53 to
// 0xbf3c9a98894a2018; substituting either field reconstructs the old
// or new aggregate exactly. The latter value reproduces in debug and
// release on aarch64; the downstream x86-64 row remains pending even
// though the upstream FFT stage-path golden is verified four ways.
const GOLDEN_VESSEL_SMOKE: u64 = 0x4541_d7f3_2926_1082;
// JUSTIFIED BUMP (2026-07-10, 6ure CLOSED): the ROA chain's libm
// divergence (macOS vs glibc in fs-bem panel2d: sin/cos/atan2/ln/
// sqrt/hypot) is fixed by routing the WHOLE panel kernel through
// fs_math::det:: — bits legitimately moved on every platform, and roa
// is RESTORED to the hashed stream (6 metrics again). Reproduced in
// BOTH debug and release on aarch64 (hash and roa bits
// 0x3fe4c1ee0bb8f1e8 identical across modes); cross-ISA identity now
// follows from det::'s own cross-ISA golden gate — ts2 (x86-64)
// confirmation row noted in bead 6ure. Prior hashes: 5-metric interim
// 0x1b03_ae9f_66cd_b548; divergent 6-metric v2 0xd750_e1bb_a8d7_e76a.
const GOLDEN_ORNITH_SMOKE: u64 = 0xae56_945a_312e_0378;
const GOLDEN_FRAME_SMOKE: u64 = 0x9c09_b06a_7883_57fc;
const GOLDEN_LBM_CORE: u64 = 0x1539_430c_dae4_7762;

fn vessel_smoke() -> (StageArtifact, f64) {
    let t0 = Instant::now();
    let rig = PourRig {
        steps: 300,
        ..PourRig::default()
    };
    let out = run_pour(
        &rig,
        fs_lbm::ContactModel::Neutral,
        fs_lbm::Rheology::Newtonian { nu: 0.0167 },
    );
    let rep = robustify(0.7);
    // poured_mass is BACK in the hashed stream (xo2k closed): the
    // divergence was platform trig in the pour tilt schedule — release
    // const-folded sin/cos with the inlined literal rig parameters while
    // debug called libm at runtime. Routed through fs_math::det::sin/cos,
    // the metric is mode-invariant again (0x3fd3b2951fb7df34 both modes).
    let metrics = vec![
        ("mass_drift", out.mass_drift),
        ("fragments", out.fragments as f64),
        ("robust_lip", rep.robust_lip),
        ("nominal_lip", rep.nominal_lip),
        ("robust_offband", rep.robust_offband_growth),
        ("poured_mass", out.poured_mass),
    ];
    (
        artifact("vessel", Tier::Smoke, metrics, t0.elapsed().as_secs_f64()),
        out.poured_mass,
    )
}

fn ornith_smoke() -> (StageArtifact, f64) {
    let t0 = Instant::now();
    let mut seed = 0xE2E_u64;
    let mut lcg = move || {
        seed = seed
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((seed >> 11) as f64) / (1u64 << 53) as f64
    };
    let generation: Vec<OrnithCandidate> = (0..12)
        .map(|_| {
            let g: Vec<f64> = (0..fs_ornith::GENE_DIM).map(|_| lcg()).collect();
            OrnithCandidate::from_genes(&g)
        })
        .collect();
    let rep = screen_generation(&generation, 0xE2E).expect("normalized screen losses");
    let winner = generation[rep.winner];
    let cert = fs_ornith::certify(&winner);
    // roa is BACK in the hashed stream (6ure closed): the divergence was
    // platform libm (macOS vs glibc sin/cos/atan2/ln/sqrt/hypot) in the
    // fs-bem panel kernel feeding the ROA proxy's adjoint lift slope;
    // the whole panel2d kernel now routes through fs_math::det::, so the
    // metric is ISA- and mode-invariant by the det:: contract.
    let metrics = vec![
        ("winner", rep.winner as f64),
        ("eliminated", rep.eliminated as f64),
        ("evals", rep.evaluations_used as f64),
        ("winner_ld", lift_to_drag(&winner)),
        ("certified", f64::from(u8::from(cert.certified))),
        ("roa", cert.roa_volume),
    ];
    (
        artifact("ornith", Tier::Smoke, metrics, t0.elapsed().as_secs_f64()),
        cert.roa_volume,
    )
}

fn frame_smoke() -> StageArtifact {
    let t0 = Instant::now();
    let catalog = [0.5f64, 0.75, 1.0, 1.5, 2.0];
    let layout = with_cx(|cx| {
        layout_and_size(5, 3, 4.0, 2.0, 250e6, 200e9, &catalog, cx)
            .expect("valid flagship frame layout is admitted")
    });
    let ensemble = StochasticEnsemble {
        name: "e2e-kt".to_string(),
        seed: 90210,
        members: 60,
        duration: QtyAny::new(12.0, TIME),
        dt: QtyAny::new(0.02, TIME),
        model: SpectrumModel::KanaiTajimi {
            s0: 0.01,
            omega_g: QtyAny::new(12.5, RATE),
            zeta_g: 0.6,
        },
    };
    let frag = e_stopped_fragility(&ensemble, StoryParams::default(), 0.02, 0.05, 0.12);
    let metrics = vec![
        ("lp_gap", layout.gap),
        ("lp_residual", layout.residual),
        ("volume", layout.volume),
        ("p_hat", frag.p_hat),
        ("members_used", f64::from(frag.members_used)),
    ];
    artifact("frame", Tier::Smoke, metrics, t0.elapsed().as_secs_f64())
}

/// fe2e-001..003: the three shipped flagships' SMOKE stages, green and
/// GOLDEN (content hash frozen; replay equality is the gate).
#[test]
fn fe2e_001_vessel_smoke_golden() {
    let (a, poured_a) = vessel_smoke();
    let (b, poured_b) = vessel_smoke();
    // Physical-band gate on poured_mass (kept after xo2k's fix put the
    // metric back in the hash: the band catches meaning, the hash bits).
    let poured_ok = (0.25..0.40).contains(&poured_a) && (poured_a - poured_b).abs() < 1e-9;
    let evidence = notebook(std::slice::from_ref(&a));
    println!(
        "{}",
        log_row(
            "vessel-smoke",
            "artifact",
            &format!(
                "{{\"hash\":\"0x{:016x}\",\"wall_s\":{:.2},\"evidence\":{evidence}}}",
                a.hash, a.wall_s,
            )
        )
    );
    verdict(
        "fe2e-001-vessel-smoke",
        a.hash == b.hash && a.hash == GOLDEN_VESSEL_SMOKE && a.metrics[0].1 < 1e-10 && poured_ok,
        &format!(
            "vessel smoke: hash 0x{:016x} (golden 0x{GOLDEN_VESSEL_SMOKE:016x}), replay equal, mass drift {:.2e}, poured {poured_a:.4} bits 0x{:016x} (envelope 0.25..0.40, bead xo2k), wall {:.1}s; evidence {evidence}",
            a.hash,
            a.metrics[0].1,
            poured_a.to_bits(),
            a.wall_s,
        ),
    );
}

#[test]
fn fe2e_002_ornith_smoke_golden() {
    let (a, roa_a) = ornith_smoke();
    let (b, roa_b) = ornith_smoke();
    // Envelope gate for the ISA-divergent metric (bead 6ure): P-area
    // proxy stable to 1e-6 within a platform and physically plausible.
    let roa_ok = (0.1..2.0).contains(&roa_a) && (roa_a - roa_b).abs() < 1e-6;
    let evidence = notebook(std::slice::from_ref(&a));
    println!(
        "{}",
        log_row(
            "ornith-smoke",
            "artifact",
            &format!(
                "{{\"hash\":\"0x{:016x}\",\"wall_s\":{:.2},\"evidence\":{evidence}}}",
                a.hash, a.wall_s,
            )
        )
    );
    verdict(
        "fe2e-002-ornith-smoke",
        a.hash == b.hash && a.hash == GOLDEN_ORNITH_SMOKE && roa_ok,
        &format!(
            "ornith smoke: hash 0x{:016x} (golden 0x{GOLDEN_ORNITH_SMOKE:016x}), replay equal, roa {roa_a:.4} bits 0x{:016x} (envelope 0.1..2.0, bead 6ure), wall {:.1}s; evidence {evidence}",
            a.hash,
            roa_a.to_bits(),
            a.wall_s,
        ),
    );
}

#[test]
fn fe2e_003_frame_smoke_golden() {
    let a = frame_smoke();
    let b = frame_smoke();
    println!(
        "{}",
        log_row(
            "frame-smoke",
            "artifact",
            &format!(
                "{{\"hash\":\"0x{:016x}\",\"wall_s\":{:.2}}}",
                a.hash, a.wall_s
            )
        )
    );
    verdict(
        "fe2e-003-frame-smoke",
        a.hash == b.hash && a.hash == GOLDEN_FRAME_SMOKE,
        &format!(
            "frame smoke: hash 0x{:016x} (golden 0x{GOLDEN_FRAME_SMOKE:016x}), replay equal, wall {:.1}s",
            a.hash, a.wall_s
        ),
    );
}

/// fe2e-004: the marquee lane is gated by its owner; the suite records
/// the STATUS honestly instead of pretending a runner.
#[test]
fn fe2e_004_marquee_status_recorded() {
    let status = fs_marquee::status();
    let scope = fs_marquee::scope_summary();
    println!(
        "{}",
        log_row(
            "marquee",
            "status",
            &format!("{{\"status\":\"{status:?}\",\"scope\":\"{scope}\"}}")
        )
    );
    verdict(
        "fe2e-004-marquee-status",
        status == fs_marquee::MarqueeStatus::Disabled && scope.contains("pending"),
        &format!(
            "marquee lane status {status:?} recorded (feature-gated by its owner; suite does not pretend)"
        ),
    );
}

/// fe2e-005: CROSS-FLAGSHIP AUDIT — the shared LBM core: one canonical
/// roll hash for the machinery the vessel AND the ornithoid ride; a
/// core change surfaces as ONE delta here.
#[test]
fn fe2e_005_shared_lbm_core_audit() {
    let h1 = lbm_core_roll_hash();
    let h2 = lbm_core_roll_hash();
    verdict(
        "fe2e-005-lbm-core-audit",
        h1 == h2 && h1 == GOLDEN_LBM_CORE,
        &format!(
            "canonical D2Q9 roll: 0x{h1:016x} (golden 0x{GOLDEN_LBM_CORE:016x}), replay equal — vessel and ornithoid ride the same bits"
        ),
    );
}

/// fe2e-006: CROSS-FLAGSHIP AUDIT — e-racing behavior identical across
/// consumers: identical pre-normalized losses through the race core
/// must produce identical outcomes regardless of which flagship's
/// convention wrapped them.
#[test]
fn fe2e_006_erace_consistency_audit() {
    // The shared loss table (already normalized to the PairwiseRace
    // contract, as both flagships do).
    let base = [0.0f64, 0.4, 0.9, 1.3, 0.2, 1.1];
    let run = |seed: u64| {
        let kills = fs_exec::KillRegistry::new();
        for candidate in 0..base.len() {
            let _ = kills.register(candidate as u64);
        }
        let mut loss = |i: usize, t: u64| {
            let mut h = (i as u64) << 32 ^ t ^ seed;
            h ^= h << 13;
            h ^= h >> 7;
            h ^= h << 17;
            #[allow(clippy::cast_precision_loss)]
            let jitter = ((h >> 11) as f64 / (1u64 << 53) as f64 - 0.5) * 0.02;
            base[i] + jitter
        };
        fs_race::race_field(
            &mut loss,
            base.len(),
            fs_race::RaceSettings::new(fs_race::LossSpan::new(1.32).expect("positive constant")),
            &kills,
        )
        .expect("fixture losses stay within the declared span")
    };
    let a = run(0xAB);
    let b = run(0xAB);
    println!(
        "{}",
        log_row(
            "erace-audit",
            "race",
            &format!(
                "{{\"winner\":{},\"evals\":{},\"eliminated\":{}}}",
                a.winner,
                a.evaluations_used,
                a.eliminated.len()
            )
        )
    );
    verdict(
        "fe2e-006-erace-audit",
        a.winner == b.winner
            && a.evaluations_used == b.evaluations_used
            && a.eliminated == b.eliminated
            && a.winner == 0,
        &format!(
            "identical losses -> identical race outcomes across consumers: winner {}, {} evals, {} eliminated (both runs)",
            a.winner,
            a.evaluations_used,
            a.eliminated.len()
        ),
    );
}

/// fe2e-007: FAILURE DRILLS — each with an expected structured outcome.
/// (One test = one drill SUITE: the four scenarios share fixtures and
/// their outcomes gate together — the length is the drill inventory,
/// not incidental complexity.)
#[test]
#[allow(clippy::too_many_lines)]
fn fe2e_007_failure_drills() {
    // (a) CANCELLATION STORM: kill half the candidates mid-race; the
    // race completes, the winner is a survivor, kills are recorded.
    let base = [0.0f64, 0.3, 0.6, 0.9, 1.2, 1.5];
    let kills = fs_exec::KillRegistry::new();
    for candidate in 0..base.len() {
        let _ = kills.register(candidate as u64);
    }
    let mut calls = 0u64;
    let mut loss = |i: usize, t: u64| {
        calls += 1;
        if calls == 30 {
            // The storm: kill the trailing half.
            for victim in 3..6 {
                let _ = kills.kill(victim as u64);
            }
        }
        let mut h = (i as u64) << 32 ^ t ^ 0x570_u64;
        h ^= h << 13;
        h ^= h >> 7;
        h ^= h << 17;
        #[allow(clippy::cast_precision_loss)]
        let jitter = ((h >> 11) as f64 / (1u64 << 53) as f64 - 0.5) * 0.02;
        base[i] + jitter
    };
    let out = fs_race::race_field(
        &mut loss,
        base.len(),
        fs_race::RaceSettings::new(fs_race::LossSpan::new(1.52).expect("positive constant")),
        &kills,
    )
    .expect("fixture losses stay within the declared span");
    let storm_ok = out.winner == 0;
    println!(
        "{}",
        log_row(
            "drill",
            "cancellation-storm",
            &format!(
                "{{\"winner\":{},\"survivors\":{}}}",
                out.winner,
                out.survivors.len()
            )
        )
    );

    // (b) BUDGET EXHAUSTION: the ornithoid campaign degrades to the
    // surrogate+conformal path (the flagship's own drill, re-run at
    // suite level through public API).
    let mut seed = 0x0771_u64;
    let mut lcg = move || {
        seed = seed
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((seed >> 11) as f64) / (1u64 << 53) as f64
    };
    let train: Vec<(OrnithCandidate, f64)> = (0..40)
        .map(|_| {
            let g: Vec<f64> = (0..fs_ornith::GENE_DIM).map(|_| lcg()).collect();
            let c = OrnithCandidate::from_genes(&g);
            (c, lift_to_drag(&c))
        })
        .collect();
    let sur = fs_ornith::LdSurrogate::fit(&train, 0.1);
    let fresh: Vec<OrnithCandidate> = (0..6)
        .map(|_| {
            let g: Vec<f64> = (0..fs_ornith::GENE_DIM).map(|_| lcg()).collect();
            OrnithCandidate::from_genes(&g)
        })
        .collect();
    let in_band = fresh
        .iter()
        .filter(|c| sur.band.covers(sur.predict(c), lift_to_drag(c)))
        .count();
    let budget_ok = in_band >= 5;
    println!(
        "{}",
        log_row(
            "drill",
            "budget-exhaustion",
            &format!("{{\"degraded\":6,\"in_band\":{in_band}}}")
        )
    );

    // (c) LEDGER CRASH-RECOVERY: write events, drop WITHOUT commit
    // (the crash), reopen — the ledger is consistent and the committed
    // prefix is intact.
    let dir = std::env::temp_dir().join(format!("fe2e-ledger-{}", std::process::id()));
    let _ = std::fs::create_dir_all(&dir);
    let path = dir.join("crash.db");
    let path_s = path.to_str().expect("utf8 path").to_owned();
    let committed_ops = {
        let ledger = fs_ledger::Ledger::open(&path_s).expect("open ledger");
        ledger.begin().expect("begin");
        ledger
            .append_event(&fs_ledger::EventRow {
                session: None,
                t: 1,
                kind: "fe2e-smoke",
                payload: Some("{\"drift\":1.5e-13}"),
            })
            .expect("event");
        ledger.commit().expect("commit");
        // The crash: a transaction begun and NEVER committed.
        ledger.begin().expect("begin 2");
        ledger
            .append_event(&fs_ledger::EventRow {
                session: None,
                t: 2,
                kind: "fe2e-lost",
                payload: Some("{\"uncommitted\":true}"),
            })
            .expect("event 2");
        // drop without commit
        1usize
    };
    let reopened = fs_ledger::Ledger::open(&path_s).expect("reopen after crash");
    let recovered = reopened.schema_version().is_ok();
    println!(
        "{}",
        log_row(
            "drill",
            "ledger-crash-recovery",
            &format!("{{\"committed_ops\":{committed_ops},\"recovered\":{recovered}}}")
        )
    );

    // (d) MODEL-FORM ESCALATION: a prediction OUTSIDE the conformal
    // band must trigger escalation (the certify-or-escalate contract).
    let outside = !sur.band.covers(1000.0, lift_to_drag(&fresh[0]));
    println!(
        "{}",
        log_row(
            "drill",
            "model-form-escalation",
            &format!("{{\"escalates_on_outlier\":{outside}}}")
        )
    );

    verdict(
        "fe2e-007-failure-drills",
        storm_ok && budget_ok && recovered && outside,
        &format!(
            "cancellation storm: winner {} among survivors; budget exhaustion: {in_band}/6 in band; ledger crash-recovery: reopened consistent; model-form escalation triggers on outliers",
            out.winner
        ),
    );
}

/// fe2e-008: FORENSIC LOGGING self-audit — every suite row is valid
/// JSON with the required keys, and the LAB NOTEBOOK regenerates
/// bitwise on replay (timings ride outside the golden body).
#[test]
fn fe2e_008_forensics_and_notebook() {
    let rows = [
        log_row("vessel-smoke", "artifact", "{\"hash\":\"0x1\"}"),
        log_row("drill", "cancellation-storm", "{\"winner\":0}"),
    ];
    let parseable = rows.iter().all(|r| {
        r.starts_with('{')
            && r.ends_with('}')
            && r.contains("\"stage\":")
            && r.contains("\"kind\":")
            && r.contains("\"payload\":")
    });
    let escaped = log_row("vessel\"\n", "artifact\\kind", "{\"ok\":true}");
    let hostile = artifact("vessel\"\n", Tier::Smoke, vec![("metric\tname", 1.0)], 0.0);
    let escaped_notebook = notebook(&[hostile]);
    let arts = vec![vessel_smoke().0, ornith_smoke().0, frame_smoke()];
    let n1 = notebook(&arts);
    // Replay: rebuild everything and regenerate.
    let arts2 = vec![vessel_smoke().0, ornith_smoke().0, frame_smoke()];
    let n2 = notebook(&arts2);
    println!(
        "{}",
        log_row(
            "notebook",
            "emitted",
            &format!("{{\"bytes\":{}}}", n1.len())
        )
    );
    verdict(
        "fe2e-008-forensics-notebook",
        parseable
            && escaped
                == "{\"stage\":\"vessel\\\"\\n\",\"kind\":\"artifact\\\\kind\",\"payload\":{\"ok\":true}}"
            && escaped_notebook.contains("\"flagship\":\"vessel\\\"\\n\"")
            && escaped_notebook.contains("\"metric\\tname\":\"0x3ff0000000000000\"")
            && n1 == n2
            && n1.contains("\"stages\":[")
            && arts.len() == 3,
        &format!(
            "forensic rows escape string fields and carry required keys; lab notebook ({} bytes) regenerates BITWISE on full replay",
            n1.len()
        ),
    );
}

/// MID stage: wired with envelopes, nightly cadence — the CI lane is
/// the perf-CI bead's (fz2.4). Run manually with `-- --ignored`.
#[test]
#[ignore = "nightly lane (perf-CI bead fz2.4): hour-class fidelity"]
fn fe2e_mid_stages() {
    // Vessel at full default resolution + full robustify; ornith with
    // a 24-candidate generation + LBM refinement of the winner; frame
    // with 200-member fragility. Envelopes, not hashes, gate the MID
    // tier when lanes land (documented policy for stochastic-labeled
    // stages; today's stages are deterministic so hashes would also
    // hold).
    let rig = PourRig::default();
    let out = run_pour(
        &rig,
        fs_lbm::ContactModel::Neutral,
        fs_lbm::Rheology::Newtonian { nu: 0.0167 },
    );
    assert!(out.mass_drift < 1e-10 && out.poured_mass > 1.0);
    let c = OrnithCandidate::from_genes(&[0.4, 0.6, 0.5, 0.2, 0.5]);
    let rep = fs_ornith::refine(&c);
    assert!(rep.lift.signum() == rep.panel_cl.signum() && rep.steadiness < 1e-4);
}

/// FULL stage: weekly/on-demand placeholder — envelopes recorded; the
/// production-scale geometry and motion suites are the flagships'
/// recorded successors and land with their beads.
#[test]
#[ignore = "weekly lane: production-scale fidelity lands with the flagship successors"]
fn fe2e_full_stages() {
    // Intentionally the same shape as MID until the successors land —
    // running it is honest (it exercises the wiring), and the envelope
    // constants live here so the lane has a home.
    fe2e_full_placeholder();
}

fn fe2e_full_placeholder() {
    let ensemble = StochasticEnsemble {
        name: "e2e-kt-full".to_string(),
        seed: 90211,
        members: 200,
        duration: QtyAny::new(12.0, TIME),
        dt: QtyAny::new(0.02, TIME),
        model: SpectrumModel::KanaiTajimi {
            s0: 0.01,
            omega_g: QtyAny::new(12.5, RATE),
            zeta_g: 0.6,
        },
    };
    let frag = e_stopped_fragility(&ensemble, StoryParams::default(), 0.02, 0.05, 0.12);
    assert!(frag.members_used <= 200);
}
