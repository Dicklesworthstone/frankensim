//! The flagship e2e battery (bead mye.5): SMOKE stages for all
//! flagships with frozen golden hashes, cross-flagship audits, failure
//! drills with expected structured outcomes, forensic-log self-audit,
//! and the bitwise-replayable lab notebook. MID/FULL stages are wired
//! behind `#[ignore]` — their CI cadence belongs to the perf-CI lanes.

use std::fmt::Write as _;
use std::time::Instant;

use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_flagship_e2e::{StageArtifact, Tier, artifact, lbm_core_roll_hash, log_row, notebook};
use fs_frame::history::StoryParams;
use fs_frame::{e_stopped_fragility, layout_and_size};
use fs_ornith::param::OrnithCandidate;
use fs_ornith::screen::{lift_to_drag, screen_generation};
use fs_qty::{Dims, QtyAny};
use fs_scenario::ensemble::{SpectrumModel, StochasticEnsemble};
use fs_surrogate::{Decision, certify_or_escalate};
use fs_vessel::pour::{PourRig, run_pour};
use fs_vessel::robustify;

const SUITE: &str = "fs-flagship-e2e/e2e-battery";
const FIXED_INPUT_SEED: u64 = 0;
const ORNITH_INPUT_SEED: u64 = 0xE2E;
const FRAME_INPUT_SEED: u64 = 90_210;
const FRAME_EXECUTION_SEED: u64 = 0xF1A6_5A1D;
const ERACE_INPUT_SEED: u64 = 0xAB;
const CANCELLATION_INPUT_SEED: u64 = 0x570;
const SURROGATE_INPUT_SEED: u64 = 0x0771;

fn verdict(name: &str, pass: bool, details: &str, seed: u64) {
    let mut emitter = fs_obs::Emitter::new(SUITE, name);
    let event = emitter.emit(
        if pass {
            fs_obs::Severity::Info
        } else {
            fs_obs::Severity::Error
        },
        fs_obs::EventKind::ConformanceCase {
            suite: SUITE.to_string(),
            case: name.to_string(),
            pass,
            detail: details.to_string(),
            seed,
        },
        None,
    );
    fs_obs::lint_failure_record(&event)
        .expect("flagship e2e verdict must carry replayable failure evidence");
    let line = event.to_jsonl();
    fs_obs::validate_line(&line).expect("flagship e2e verdict must use the fs-obs wire schema");
    println!("{line}");
    assert!(pass, "{name}: {details}");
}

fn forensic_row(stage: &str, kind: &str, json: String) {
    let mut emitter = fs_obs::Emitter::new(SUITE, stage);
    let line = emitter
        .emit(
            fs_obs::Severity::Info,
            fs_obs::EventKind::Custom {
                name: kind.to_string(),
                json,
            },
            None,
        )
        .to_jsonl();
    fs_obs::validate_line(&line)
        .expect("flagship forensic companion must use the fs-obs wire schema");
    println!("{line}");
}

fn json_string(value: &str) -> String {
    let mut out = String::with_capacity(value.len() + 2);
    out.push('"');
    for ch in value.chars() {
        match ch {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\u{0008}' => out.push_str("\\b"),
            '\u{000c}' => out.push_str("\\f"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            ch if ch <= '\u{001f}' => {
                let _ = write!(out, "\\u{:04x}", u32::from(ch));
            }
            ch => out.push(ch),
        }
    }
    out.push('"');
    out
}

const TIME: Dims = Dims([0, 0, 1, 0, 0, 0]);
const RATE: Dims = Dims([0, 0, -1, 0, 0, 0]);

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: FRAME_EXECUTION_SEED,
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
    let mut seed = ORNITH_INPUT_SEED;
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
    let rep = screen_generation(&generation, ORNITH_INPUT_SEED).expect("normalized screen losses");
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
        seed: FRAME_INPUT_SEED,
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
    forensic_row(
        "vessel-smoke",
        "artifact",
        format!(
            "{{\"hash\":\"0x{:016x}\",\"wall_s\":{:.2},\"evidence\":{evidence},\
             \"input_seed\":{FIXED_INPUT_SEED}}}",
            a.hash, a.wall_s,
        ),
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
        FIXED_INPUT_SEED,
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
    forensic_row(
        "ornith-smoke",
        "artifact",
        format!(
            "{{\"hash\":\"0x{:016x}\",\"wall_s\":{:.2},\"evidence\":{evidence},\
             \"input_seed\":{ORNITH_INPUT_SEED}}}",
            a.hash, a.wall_s,
        ),
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
        ORNITH_INPUT_SEED,
    );
}

#[test]
fn fe2e_003_frame_smoke_golden() {
    let a = frame_smoke();
    let b = frame_smoke();
    forensic_row(
        "frame-smoke",
        "artifact",
        format!(
            "{{\"hash\":\"0x{:016x}\",\"wall_s\":{:.2},\
             \"input_seed\":{FRAME_INPUT_SEED},\
             \"execution_seed\":{FRAME_EXECUTION_SEED}}}",
            a.hash, a.wall_s
        ),
    );
    verdict(
        "fe2e-003-frame-smoke",
        a.hash == b.hash && a.hash == GOLDEN_FRAME_SMOKE,
        &format!(
            "frame smoke: hash 0x{:016x} (golden 0x{GOLDEN_FRAME_SMOKE:016x}), replay equal, wall {:.1}s",
            a.hash, a.wall_s
        ),
        FRAME_INPUT_SEED,
    );
}

/// fe2e-004: the marquee lane is gated by its owner; the suite records
/// the STATUS honestly instead of pretending a runner.
#[test]
fn fe2e_004_marquee_status_recorded() {
    let status = fs_marquee::status();
    let scope = fs_marquee::scope_summary();
    let status_json = json_string(&format!("{status:?}"));
    let scope_json = json_string(&scope);
    forensic_row(
        "marquee",
        "status",
        format!(
            "{{\"status\":{status_json},\"scope\":{scope_json},\
             \"input_seed\":{FIXED_INPUT_SEED}}}"
        ),
    );
    verdict(
        "fe2e-004-marquee-status",
        status == fs_marquee::MarqueeStatus::Disabled && scope.contains("pending"),
        &format!(
            "marquee lane status {status:?} recorded (feature-gated by its owner; suite does not pretend)"
        ),
        FIXED_INPUT_SEED,
    );
}

/// fe2e-005: SHARED-CORE AUDIT — the uniform-tau D2Q9 collide/stream
/// path both flagships ride, pinned by one roll hash. Its coverage is
/// exactly the fixture's (see `lbm_core_roll_hash`): rheology, free
/// surface, interior-obstacle bounce-back, non-periodic inlet/outlet
/// handling, momentum-exchange variants and D3Q19 are NOT pinned here.
#[test]
fn fe2e_005_shared_lbm_core_audit() {
    let h1 = lbm_core_roll_hash();
    let h2 = lbm_core_roll_hash();
    verdict(
        "fe2e-005-lbm-core-audit",
        h1 == h2 && h1 == GOLDEN_LBM_CORE,
        &format!(
            "canonical uniform-tau D2Q9 roll (periodic-x, wall-bounded y, 50 plain steps): \
             0x{h1:016x} (golden 0x{GOLDEN_LBM_CORE:016x}), replay equal — one shared audit \
             point for the collide/stream core; rheology, free surface, interior-obstacle \
             bounce-back, inlet/outlet columns, momentum-exchange variants and D3Q19 are \
             NOT covered by this hash"
        ),
        FIXED_INPUT_SEED,
    );
}

/// The ornithoid's declared race span (`fs_ornith::screen`): its
/// normalized base losses lie in `[0, 1.5]` and the jitter has total
/// width `0.02`. Reconstructed here so the suite can drive the SHARED
/// race core under the consumer's own convention and compare.
const ORNITH_DECLARED_SPAN: f64 = 1.52;
/// The ceiling the ornithoid normalizes its base losses onto.
const ORNITH_NORMALIZED_CEILING: f64 = 1.5;

/// Drive the shared `fs_race` core over `generation` under the
/// ORNITHOID's normalization convention, reconstructed at suite level
/// from its documented contract (`−L/D`, shifted to zero and scaled onto
/// `[0, ceiling]`, plus the ±0.01 jitter, raced under `span`).
///
/// This is the audit's independent side: `fs_ornith::screen_generation`
/// is the consumer, this is the core driven the way the consumer says it
/// drives it. If the consumer's normalization or declared span moves,
/// the two disagree — which is the drift the audit exists to catch.
fn erace_core_under_ornith_convention(
    generation: &[OrnithCandidate],
    seed: u64,
    span: f64,
    ceiling: f64,
) -> fs_race::RaceOutcome {
    let base: Vec<f64> = generation.iter().map(|c| -lift_to_drag(c)).collect();
    let lo = base.iter().fold(f64::INFINITY, |m, &v| m.min(v));
    let hi = base.iter().fold(f64::NEG_INFINITY, |m, &v| m.max(v));
    let scale = ceiling / (hi - lo).max(1e-9);
    let kills = fs_exec::KillRegistry::new();
    for candidate in 0..generation.len() {
        let _ = kills.register(candidate as u64);
    }
    let mut loss = |i: usize, t: u64| {
        let mut h = (i as u64) << 32 ^ t ^ seed;
        h ^= h << 13;
        h ^= h >> 7;
        h ^= h << 17;
        #[allow(clippy::cast_precision_loss)]
        let jitter = ((h >> 11) as f64 / (1u64 << 53) as f64 - 0.5) * 0.02;
        (base[i] - lo).mul_add(scale, jitter)
    };
    fs_race::race_field(
        &mut loss,
        generation.len(),
        fs_race::RaceSettings::new(fs_race::LossSpan::new(span).expect("positive constant")),
        &kills,
    )
    .expect("fixture losses stay within the declared span")
}

/// The vessel's screening lip fixture (`fs-vessel` vsl-005), rebuilt
/// here so the audit drives the vessel wrapper on its own home ground.
const VESSEL_SCREEN_LIPS: [f64; 5] = [0.6, 1.0, 1.6, 2.2, 2.8];
/// The seed the vessel's screening jitter hashes with.
const VESSEL_SCREEN_SEED: u64 = 0x7E55;

/// The VESSEL's declared racing convention, as three numbers the suite
/// can perturb. `declared()` reads them off `fs-vessel`'s public
/// constants, so a change in the library reaches this audit.
#[derive(Debug, Clone, Copy)]
struct VesselConvention {
    /// Multiplier applied to `(base + jitter)` before racing.
    scale: f64,
    /// Total width of the per-observation validator jitter.
    jitter_width: f64,
    /// Slack folded into the declared support on top of the fixture
    /// spread. The shipped convention sets this EQUAL to `jitter_width`
    /// — that equality is the soundness of the declaration, and drifting
    /// the two apart is the drift class this audit exists to catch.
    span_slack: f64,
}

impl VesselConvention {
    /// The convention `fs-vessel` actually ships.
    fn declared() -> Self {
        VesselConvention {
            scale: fs_vessel::race::SCREEN_SCALE,
            jitter_width: fs_vessel::race::SCREEN_JITTER_WIDTH,
            span_slack: fs_vessel::race::SCREEN_JITTER_WIDTH,
        }
    }
}

/// Drive the shared `fs_race` core over `base` under the VESSEL's
/// normalization convention, reconstructed at suite level from its
/// documented contract (`SCREEN_SCALE × (base + jitter)`, raced under
/// `SCREEN_SCALE × (fixture spread + jitter width)`).
///
/// The audit's independent side for the vessel, exactly as
/// [`erace_core_under_ornith_convention`] is for the ornithoid.
fn erace_core_under_vessel_convention(
    base: &[f64],
    seed: u64,
    conv: VesselConvention,
) -> Result<fs_race::RaceOutcome, fs_race::RaceError> {
    let base_span = base.iter().copied().fold(f64::NEG_INFINITY, f64::max)
        - base.iter().copied().fold(f64::INFINITY, f64::min);
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
        let jitter = ((h >> 11) as f64 / (1u64 << 53) as f64 - 0.5) * conv.jitter_width;
        (base[i] + jitter) * conv.scale
    };
    fs_race::race_field(
        &mut loss,
        base.len(),
        fs_race::RaceSettings::new(
            fs_race::LossSpan::new(conv.scale * (base_span + conv.span_slack))
                .expect("positive derived span"),
        ),
        &kills,
    )
}

/// `(winner, evaluations_used, eliminated)` — the comparable part of a
/// race outcome, whichever consumer produced it.
type RaceTriple = (usize, u64, usize);

/// Does the vessel's public screening wrapper agree with the shared
/// race core driven under `conv`? Returns both triples so callers can
/// report the disagreement, not just its existence.
fn vessel_wrapper_agrees_with_core(
    conv: VesselConvention,
) -> (bool, RaceTriple, Result<RaceTriple, fs_race::RaceError>) {
    let wrapper = fs_vessel::race::screen_lips(&VESSEL_SCREEN_LIPS, VESSEL_SCREEN_SEED)
        .expect("the vessel lip fixture admits a verdict");
    let a = (wrapper.winner, wrapper.evaluations_used, wrapper.eliminated);
    let core = erace_core_under_vessel_convention(&wrapper.losses, VESSEL_SCREEN_SEED, conv)
        .map(|out| (out.winner, out.evaluations_used, out.eliminated.len()));
    (core.as_ref().is_ok_and(|b| &a == b), a, core)
}

/// BOTH consumers' declared conventions over ONE shared loss table.
///
/// The table is the ornithoid's OWN screening losses (`−L/D` per
/// candidate), taken straight off `fs_ornith::screen_generation`'s
/// report, so the two consumers provably race the same numbers. The
/// ornithoid side is its public wrapper; the vessel side is driven
/// under `conv` so a drift can be seeded — at `VesselConvention::declared()`
/// this reproduces the shipped `fs_vessel::race::race_base_losses`
/// wrapper bit-for-bit, which fe2e-006 asserts rather than assumes.
fn shared_table_consumer_outcomes(
    conv: VesselConvention,
) -> (RaceTriple, Result<RaceTriple, fs_race::RaceError>) {
    let generation = erace_audit_generation();
    let ornith =
        screen_generation(&generation, ORNITH_INPUT_SEED).expect("normalized screen losses");
    let vessel = erace_core_under_vessel_convention(&ornith.losses, ORNITH_INPUT_SEED, conv)
        .map(|out| (out.winner, out.evaluations_used, out.eliminated.len()));
    (
        (ornith.winner, ornith.evaluations_used, ornith.eliminated),
        vessel,
    )
}

/// The cross-consumer claim, stated precisely: over one shared loss
/// table the two declared conventions must make the same SELECTION —
/// same winner, same number of eliminated candidates. That pins the
/// whole survivor set only when the elimination is TOTAL (n−1 of n),
/// which the audit records as `cross_consumer_survivor_set_pinned`
/// rather than assuming. Evaluation counts are convention-dependent and
/// are NOT claimed equal; they are reported.
fn selections_agree(a: RaceTriple, b: RaceTriple) -> bool {
    a.0 == b.0 && a.2 == b.2
}

/// The fe2e-006 fixture generation (the ornith smoke generation, rebuilt
/// so the audit and its falsification test share one input).
fn erace_audit_generation() -> Vec<OrnithCandidate> {
    let mut seed = ORNITH_INPUT_SEED;
    let mut lcg = move || {
        seed = seed
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((seed >> 11) as f64) / (1u64 << 53) as f64
    };
    (0..12)
        .map(|_| {
            let g: Vec<f64> = (0..fs_ornith::GENE_DIM).map(|_| lcg()).collect();
            OrnithCandidate::from_genes(&g)
        })
        .collect()
}

/// Does the ornithoid's public screening wrapper agree with the shared
/// race core driven under `span`/`ceiling`? Returns the two outcomes so
/// callers can report the disagreement, not just its existence.
fn ornith_wrapper_agrees_with_core(
    span: f64,
    ceiling: f64,
) -> (bool, (usize, u64, usize), (usize, u64, usize)) {
    let generation = erace_audit_generation();
    let wrapper =
        screen_generation(&generation, ORNITH_INPUT_SEED).expect("normalized screen losses");
    let core = erace_core_under_ornith_convention(&generation, ORNITH_INPUT_SEED, span, ceiling);
    let a = (wrapper.winner, wrapper.evaluations_used, wrapper.eliminated);
    let b = (core.winner, core.evaluations_used, core.eliminated.len());
    (a == b, a, b)
}

/// fe2e-006: CROSS-CONSUMER E-RACE AUDIT — four claims, all measured:
///
/// 1. REPLAY DETERMINISM of the shared race core on a fixed normalized
///    loss table under a fixed declared span (one closure, one seed, run
///    twice — a same-process replay check, and labelled as one).
/// 2. CONSUMER/CORE AGREEMENT, ornithoid: `fs_ornith::screen_generation`
///    must produce the same winner, evaluation count and elimination
///    count as the shared core driven under the ornithoid's own declared
///    span and normalization.
/// 3. CONSUMER/CORE AGREEMENT, vessel: `fs_vessel::race::screen_lips`
///    must likewise agree with the shared core driven under the VESSEL's
///    declared convention (`SCREEN_SCALE × (base + jitter)`, support
///    `SCREEN_SCALE × (spread + jitter width)`).
/// 4. CROSS-CONSUMER SELECTION over ONE SHARED LOSS TABLE: both public
///    wrappers race the ornithoid's own `−L/D` table, each under its own
///    declared convention, and must reach the same SELECTION — same
///    winner, same elimination count. Whether that also pins the whole
///    survivor set (it does only when the elimination is total) is
///    recorded, not assumed.
///
/// What claim 4 does NOT assert: equal evaluation counts. The two
/// conventions differ in noise-to-span ratio (the ornithoid normalizes
/// onto a fixed ceiling with ±0.01 jitter; the vessel scales by 200 with
/// a data-derived support and ±5e-5 jitter), so the SAME table costs a
/// different number of observations under each. Those counts are
/// measured and reported side by side, not equated — asserting they were
/// equal would be the same overstatement this case was filed for
/// (bead `frankensim-extreal-program-f85xj.2.31`).
///
/// The vessel's convention became auditable only when it moved out of
/// `crates/fs-vessel/tests/battery.rs` into the `fs_vessel::race`
/// library surface: a convention that exists only in a test cannot be
/// driven by anyone else, which is exactly why the original case could
/// claim a cross-flagship audit while invoking neither flagship.
#[test]
fn fe2e_006_erace_cross_consumer_audit() {
    // (1) Same-process replay determinism of the core.
    let base = [0.0f64, 0.4, 0.9, 1.3, 0.2, 1.1];
    let replay = |seed: u64| {
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
    let a = replay(ERACE_INPUT_SEED);
    let b = replay(ERACE_INPUT_SEED);
    let replay_equal = a.winner == b.winner
        && a.evaluations_used == b.evaluations_used
        && a.eliminated == b.eliminated
        && a.winner == 0;

    // (2) The ornithoid consumer vs the core under the ornithoid's own
    // declared span. This is what actually catches a normalization drift.
    let (ornith_agrees, ornith_wrapper, ornith_core) =
        ornith_wrapper_agrees_with_core(ORNITH_DECLARED_SPAN, ORNITH_NORMALIZED_CEILING);

    // (3) The vessel consumer vs the core under the VESSEL's own declared
    // convention — the half that did not exist until `fs_vessel::race`
    // took ownership of the 200x normalization and its data-derived span.
    let declared = VesselConvention::declared();
    let (vessel_agrees, vessel_wrapper, vessel_core) = vessel_wrapper_agrees_with_core(declared);
    let vessel_report = fs_vessel::race::screen_lips(&VESSEL_SCREEN_LIPS, VESSEL_SCREEN_SEED)
        .expect("the vessel lip fixture admits a verdict");

    // (4) CROSS-CONSUMER SELECTION over one shared loss table. The
    // reconstruction used for the vessel side must be the shipped wrapper
    // — asserted, not assumed — and then both conventions must select the
    // same candidate and eliminate the same number of rivals.
    let (shared_ornith, shared_vessel_core) = shared_table_consumer_outcomes(declared);
    let shared_table = screen_generation(&erace_audit_generation(), ORNITH_INPUT_SEED)
        .expect("normalized screen losses")
        .losses;
    let shared_vessel_wrapper = fs_vessel::race::race_base_losses(&shared_table, ORNITH_INPUT_SEED)
        .expect("the shared loss table admits a vessel verdict");
    let shared_vessel = (
        shared_vessel_wrapper.winner,
        shared_vessel_wrapper.evaluations_used,
        shared_vessel_wrapper.eliminated,
    );
    let reconstruction_is_the_wrapper = shared_vessel_core.as_ref() == Ok(&shared_vessel);
    let cross_agrees = selections_agree(shared_ornith, shared_vessel);
    // Equal winner + equal elimination COUNT pins the whole survivor set
    // only when the elimination is total (n-1 of n). Recorded, not
    // assumed: with a partial field, agreeing counts would leave the
    // survivor sets themselves unaudited.
    let survivor_set_pinned =
        shared_ornith.2 == shared_table.len() - 1 && shared_vessel.2 == shared_table.len() - 1;

    forensic_row(
        "erace-audit",
        "race",
        format!(
            "{{\"replay_winner\":{},\"replay_evals\":{},\"replay_eliminated\":{},\
             \"ornith_wrapper\":[{},{},{}],\"ornith_core\":[{},{},{}],\
             \"consumer_core_agree\":{ornith_agrees},\
             \"ornith_declared_span\":{ORNITH_DECLARED_SPAN},\
             \"ornith_normalized_ceiling\":{ORNITH_NORMALIZED_CEILING},\
             \"vessel_wrapper\":[{},{},{}],\"vessel_core\":{},\
             \"vessel_consumer_core_agree\":{vessel_agrees},\
             \"vessel_declared_scale\":{},\"vessel_declared_jitter_width\":{},\
             \"vessel_declared_span\":{},\
             \"shared_table_ornith\":[{},{},{}],\"shared_table_vessel\":[{},{},{}],\
             \"shared_table_len\":{},\"vessel_reconstruction_is_wrapper\":\
             {reconstruction_is_the_wrapper},\
             \"cross_consumer_selection_agree\":{cross_agrees},\
             \"cross_consumer_survivor_set_pinned\":{survivor_set_pinned},\
             \"cross_consumer_evals_claimed_equal\":false,\
             \"input_seed\":{ERACE_INPUT_SEED},\"ornith_input_seed\":{ORNITH_INPUT_SEED},\
             \"vessel_input_seed\":{VESSEL_SCREEN_SEED}}}",
            a.winner,
            a.evaluations_used,
            a.eliminated.len(),
            ornith_wrapper.0,
            ornith_wrapper.1,
            ornith_wrapper.2,
            ornith_core.0,
            ornith_core.1,
            ornith_core.2,
            vessel_wrapper.0,
            vessel_wrapper.1,
            vessel_wrapper.2,
            match &vessel_core {
                Ok((w, e, el)) => format!("[{w},{e},{el}]"),
                Err(err) => json_string(&format!("refused: {err}")),
            },
            declared.scale,
            declared.jitter_width,
            vessel_report.declared_span,
            shared_ornith.0,
            shared_ornith.1,
            shared_ornith.2,
            shared_vessel.0,
            shared_vessel.1,
            shared_vessel.2,
            shared_table.len(),
        ),
    );
    verdict(
        "fe2e-006-erace-cross-consumer",
        replay_equal
            && ornith_agrees
            && vessel_agrees
            && reconstruction_is_the_wrapper
            && cross_agrees,
        &format!(
            "race core replays identically on a fixed normalized loss table (span 1.32): \
             winner {}, {} evals, {} eliminated in both runs. BOTH public consumers agree with \
             the shared core under their OWN declared conventions: fs_ornith::screen_generation \
             {ornith_wrapper:?} vs core {ornith_core:?} at span {ORNITH_DECLARED_SPAN}/ceiling \
             {ORNITH_NORMALIZED_CEILING}; fs_vessel::race::screen_lips {vessel_wrapper:?} vs core \
             {vessel_core:?} at scale {} / jitter width {} / declared span {:.5} (suite \
             reconstruction reproduces the shipped wrapper: {reconstruction_is_the_wrapper}). \
             CROSS-CONSUMER, one shared {}-candidate loss table (the ornithoid's own -L/D scores) \
             raced under both conventions: ornith {shared_ornith:?} vs vessel {shared_vessel:?} \
             -- selection {} (winner {} vs {}; {} vs {} eliminated of {}; elimination total on \
             both sides, so the survivor set is pinned to the winner alone: \
             {survivor_set_pinned}). Evaluation counts are convention-dependent and are NOT \
             claimed equal: {} (ornith ceiling normalization, +-0.01 jitter) vs {} (vessel 200x \
             scaling, data-derived support, +-5e-5 jitter)",
            a.winner,
            a.evaluations_used,
            a.eliminated.len(),
            declared.scale,
            declared.jitter_width,
            vessel_report.declared_span,
            shared_table.len(),
            if cross_agrees { "AGREES" } else { "DIVERGES" },
            shared_ornith.0,
            shared_vessel.0,
            shared_ornith.2,
            shared_vessel.2,
            shared_table.len(),
            shared_ornith.1,
            shared_vessel.1,
        ),
        ERACE_INPUT_SEED,
    );
}

/// REGRESSION for bead `frankensim-extreal-program-f85xj.2.31`: the
/// consumer/core agreement check must be FALSIFIABLE — it has to fail
/// when the consumer's convention and the core's declared span diverge,
/// which is precisely the drift class fe2e-006 names.
///
/// The pre-fix audit could not fail that way: it ran one closure twice
/// with one seed, so it agreed with itself under ANY convention. Both
/// halves are asserted here — the blind shape stays green under a
/// drifted span, the audit's actual check does not.
#[test]
fn fe2e_006_consumer_core_agreement_is_falsifiable() {
    // The audit's own configuration agrees.
    let (agrees, wrapper, core) =
        ornith_wrapper_agrees_with_core(ORNITH_DECLARED_SPAN, ORNITH_NORMALIZED_CEILING);
    assert!(agrees, "wrapper {wrapper:?} vs core {core:?}");

    // Seed the drift the bead describes: the consumer normalizes onto a
    // different ceiling / declares a different span than the core is driven
    // with. The agreement check must SEE it.
    let drifted_ceiling = ornith_wrapper_agrees_with_core(ORNITH_DECLARED_SPAN, 1.0);
    assert!(
        !drifted_ceiling.0,
        "a changed normalization ceiling must break consumer/core agreement: \
         wrapper {:?} vs core {:?}",
        drifted_ceiling.1, drifted_ceiling.2
    );
    let drifted_span = ornith_wrapper_agrees_with_core(3.0, ORNITH_NORMALIZED_CEILING);
    assert!(
        !drifted_span.0,
        "a changed declared span must break consumer/core agreement: \
         wrapper {:?} vs core {:?}",
        drifted_span.1, drifted_span.2
    );

    // …while the pre-fix shape (one closure, one seed, twice) is blind to
    // both — it compares an implementation with itself.
    let generation = erace_audit_generation();
    for (span, ceiling) in [
        (ORNITH_DECLARED_SPAN, ORNITH_NORMALIZED_CEILING),
        (ORNITH_DECLARED_SPAN, 1.0),
        (3.0, ORNITH_NORMALIZED_CEILING),
    ] {
        let x = erace_core_under_ornith_convention(&generation, ERACE_INPUT_SEED, span, ceiling);
        let y = erace_core_under_ornith_convention(&generation, ERACE_INPUT_SEED, span, ceiling);
        assert_eq!(
            (x.winner, x.evaluations_used, x.eliminated.len()),
            (y.winner, y.evaluations_used, y.eliminated.len()),
            "self-replay agrees under span {span} / ceiling {ceiling} — which is exactly why \
             self-replay cannot audit a consumer's convention"
        );
    }
}

/// REGRESSION for the VESSEL half of bead
/// `frankensim-extreal-program-f85xj.2.31`: the vessel consumer/core
/// agreement must be falsifiable by a drift in the vessel's own declared
/// convention, and the suite-level reconstruction the falsifier perturbs
/// must be the SHIPPED wrapper at the declared settings (otherwise the
/// drills falsify a strawman).
#[test]
fn fe2e_006_vessel_consumer_core_agreement_is_falsifiable() {
    let declared = VesselConvention::declared();
    let (agrees, wrapper, core) = vessel_wrapper_agrees_with_core(declared);
    assert!(agrees, "wrapper {wrapper:?} vs core {core:?}");

    // The soundness of the vessel's declaration is `span_slack ==
    // jitter_width`: the declared support is the fixture spread plus the
    // full width a paired difference of jitters can reach. Drift the two
    // apart and the e-process runs at a different scale — the agreement
    // check must SEE it.
    for slack in [1e-3f64, 1e-2, 1e-1] {
        let drifted = vessel_wrapper_agrees_with_core(VesselConvention {
            span_slack: slack,
            ..declared
        });
        assert!(
            !drifted.0,
            "a declared support drifted to slack {slack} must break vessel consumer/core \
             agreement: wrapper {:?} vs core {:?}",
            drifted.1, drifted.2
        );
    }

    // A noisier validator under the SAME declared support is the unsound
    // direction: the race refuses (support breach), which is a
    // disagreement the check also has to register rather than swallow.
    let noisier = vessel_wrapper_agrees_with_core(VesselConvention {
        jitter_width: 1e-2,
        ..declared
    });
    assert!(
        !noisier.0,
        "a validator noisier than the declared support must break agreement: wrapper {:?} vs \
         core {:?}",
        noisier.1, noisier.2
    );
    assert!(
        matches!(noisier.2, Err(fs_race::RaceError::PairwiseInput { .. })),
        "the breach must be a STRUCTURED refusal, not a verdict: {:?}",
        noisier.2
    );

    // Honest negative result, recorded so nobody reads it as blindness: a
    // PURE RESCALE is not a drift this (or any) audit can see, because
    // scaling every loss AND the declared support by the same factor is an
    // exact invariance of the pairwise e-process. What the audit catches is
    // a normalization/declared-support MISMATCH, not the multiplier.
    for scale in [20.0f64, 2000.0] {
        let rescaled = vessel_wrapper_agrees_with_core(VesselConvention { scale, ..declared });
        assert!(
            rescaled.0,
            "a pure rescale to {scale} must be an exact invariance, not a drift: wrapper {:?} vs \
             core {:?}",
            rescaled.1, rescaled.2
        );
    }
}

/// REGRESSION: the CROSS-CONSUMER selection claim of fe2e-006 must fail
/// when one consumer's normalization convention drifts. A cross-consumer
/// audit that stays green under a drifted convention is the original
/// defect wearing a new name.
#[test]
fn fe2e_006_cross_consumer_selection_is_falsifiable() {
    let declared = VesselConvention::declared();
    let (ornith, vessel) = shared_table_consumer_outcomes(declared);
    let vessel = vessel.expect("the shared loss table admits a vessel verdict");
    assert!(
        selections_agree(ornith, vessel),
        "ornith {ornith:?} vs vessel {vessel:?}"
    );
    // The honest limit of the claim, asserted so it cannot quietly become
    // an equality: the two conventions pay DIFFERENT evaluation counts for
    // the same table.
    assert_ne!(
        ornith.1, vessel.1,
        "the two conventions are not cost-equivalent on this table; if they ever become so, the \
         reported wording must change with the evidence"
    );

    // Drift the vessel's declared support away from its jitter width. At
    // slack 10.0 the weakened e-process fails to eliminate one rival
    // inside the round budget, so the two consumers no longer make the
    // same selection — exactly the failure the audit claims to provide.
    let drifted = shared_table_consumer_outcomes(VesselConvention {
        span_slack: 10.0,
        ..declared
    })
    .1
    .expect("a looser support still yields a verdict");
    assert!(
        !selections_agree(ornith, drifted),
        "a drifted vessel declared support must break cross-consumer selection agreement: \
         ornith {ornith:?} vs drifted vessel {drifted:?}"
    );
}

// ------------------------------------------------------------------
// fe2e-007 drill mechanics. Every field a drill row publishes is
// computed here from state the drill actually moved: a budget counter
// that runs out, a ledger read back after reopen, a certify-or-escalate
// decision on the fitted band. The high-fidelity lane and the recovery
// fault are injected so the accounting can be falsified in a regression
// test without paying for (or corrupting) the real thing.
// ------------------------------------------------------------------

/// Candidates in the fe2e-007 mini-campaign.
const CAMPAIGN_CANDIDATES: usize = 6;
/// LBM refinements the campaign can fund before the budget is exhausted.
const LBM_REFINE_BUDGET: usize = 1;

/// The fe2e-007 surrogate fixture: 40 training samples fit the L/D
/// surrogate + conformal band, then `CAMPAIGN_CANDIDATES` fresh
/// candidates form the mini-campaign. One LCG, one seed, so the drills
/// and their regression tests see identical inputs.
fn surrogate_and_campaign_fixture() -> (fs_ornith::LdSurrogate, Vec<OrnithCandidate>) {
    let mut seed = SURROGATE_INPUT_SEED;
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
    let surrogate = fs_ornith::LdSurrogate::fit(&train, 0.1);
    let campaign: Vec<OrnithCandidate> = (0..CAMPAIGN_CANDIDATES)
        .map(|_| {
            let g: Vec<f64> = (0..fs_ornith::GENE_DIM).map(|_| lcg()).collect();
            OrnithCandidate::from_genes(&g)
        })
        .collect();
    (surrogate, campaign)
}

/// Measured outcome of the budget-exhaustion drill.
#[derive(Debug, Clone, Copy)]
struct BudgetExhaustion {
    /// Candidates served by the funded high-fidelity lane.
    funded: usize,
    /// Candidates that degraded to the surrogate + conformal path.
    degraded: usize,
    /// Degraded candidates whose measured L/D fell inside the band.
    in_band: usize,
    /// Every funded answer came back finite.
    funded_lift_finite: bool,
}

/// Run the campaign under a finite LBM refinement budget. `fund` is the
/// high-fidelity lane; it returns `None` when it cannot answer.
fn run_budget_exhaustion_drill(
    surrogate: &fs_ornith::LdSurrogate,
    campaign: &[OrnithCandidate],
    budget: usize,
    fund: &mut dyn FnMut(&OrnithCandidate) -> Option<f64>,
) -> BudgetExhaustion {
    let mut remaining = budget;
    let mut out = BudgetExhaustion {
        funded: 0,
        degraded: 0,
        in_band: 0,
        funded_lift_finite: true,
    };
    for candidate in campaign {
        if remaining > 0 {
            remaining -= 1;
            out.funded += 1;
            let answer = fund(candidate);
            out.funded_lift_finite &= answer.is_some_and(f64::is_finite);
        } else {
            out.degraded += 1;
            let predicted = surrogate.predict(candidate);
            if surrogate.band.covers(predicted, lift_to_drag(candidate)) {
                out.in_band += 1;
            }
        }
    }
    out
}

/// The gate fe2e-007 applies to the budget drill: the budget must have
/// actually run out (funded lane exercised AND degradation observed) and
/// the degraded estimates must mostly land inside their conformal band.
fn budget_exhaustion_gate(out: &BudgetExhaustion) -> bool {
    out.funded == LBM_REFINE_BUDGET
        && out.degraded == CAMPAIGN_CANDIDATES - LBM_REFINE_BUDGET
        && out.funded_lift_finite
        && out.in_band >= 4
}

/// A seeded recovery fault for the ledger crash-recovery drill.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum LedgerFault {
    /// Honest crash: the second transaction is dropped without commit.
    None,
    /// A recovery that loses the committed prefix.
    CommittedPrefixLost,
    /// A recovery that replays the uncommitted transaction.
    UncommittedReplayed,
}

/// Bytes committed before the crash — read back verbatim after reopen.
const LEDGER_COMMITTED_PAYLOAD: &[u8] = b"{\"stage\":\"fe2e-smoke\",\"drift\":1.5e-13}";
/// Bytes written INSIDE the uncommitted transaction — must not survive.
const LEDGER_UNCOMMITTED_PAYLOAD: &[u8] = b"{\"stage\":\"fe2e-lost\",\"uncommitted\":true}";

/// Measured outcome of the ledger crash-recovery drill.
// Each flag is a SEPARATE measured fact about the reopened ledger, and the
// drill row publishes them individually so a failure names itself. Collapsing
// them into one enum would re-hide exactly what bead
// `frankensim-extreal-program-f85xj.2.33` says must be visible.
#[allow(clippy::struct_excessive_bools)]
#[derive(Debug, Clone)]
struct LedgerRecovery {
    /// The reopened ledger reports a readable schema version. This alone
    /// was the entire pre-fix gate, and it survives both seeded faults.
    schema_version_readable: bool,
    /// Rows in `events` after reopen (the committed prefix is exactly 1).
    committed_events: u64,
    /// The committed artifact materializes with byte-identical content.
    committed_payload_readback: bool,
    /// The artifact written in the uncommitted transaction is absent.
    uncommitted_absent: bool,
    /// Stored artifacts still hash to their recorded identities.
    artifact_integrity_clean: bool,
}

impl LedgerRecovery {
    /// Crash recovery held: the committed prefix survived intact and the
    /// uncommitted write did not.
    fn recovered(&self) -> bool {
        self.schema_version_readable
            && self.committed_events == 1
            && self.committed_payload_readback
            && self.uncommitted_absent
            && self.artifact_integrity_clean
    }
}

/// Delete a ledger database and any sidecar files sharing its name.
fn remove_ledger_files(path: &std::path::Path) {
    let Some(dir) = path.parent() else { return };
    let Some(name) = path.file_name().and_then(std::ffi::OsStr::to_str) else {
        return;
    };
    let Ok(entries) = std::fs::read_dir(dir) else {
        return;
    };
    for entry in entries.flatten() {
        if entry
            .file_name()
            .to_str()
            .is_some_and(|f| f.starts_with(name))
        {
            let _ = std::fs::remove_file(entry.path());
        }
    }
}

/// Commit a prefix, crash mid-transaction, reopen, and READ THE PREFIX
/// BACK. `fault` seeds a broken recovery so the drill's own gate can be
/// falsified.
fn run_ledger_crash_recovery_drill(path: &str, fault: LedgerFault) -> LedgerRecovery {
    remove_ledger_files(std::path::Path::new(path));
    let (committed, lost) = {
        let ledger = fs_ledger::Ledger::open(path).expect("open ledger");
        ledger.begin().expect("begin");
        ledger
            .append_event(&fs_ledger::EventRow {
                session: None,
                t: 1,
                kind: "fe2e-smoke",
                payload: Some("{\"drift\":1.5e-13}"),
            })
            .expect("event");
        let committed = ledger
            .put_artifact("fe2e-drill", LEDGER_COMMITTED_PAYLOAD, None)
            .expect("committed artifact");
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
        let lost = ledger
            .put_artifact("fe2e-drill", LEDGER_UNCOMMITTED_PAYLOAD, None)
            .expect("uncommitted artifact");
        if fault == LedgerFault::UncommittedReplayed {
            ledger
                .commit()
                .expect("seeded fault: replay the uncommitted transaction");
        }
        // drop without commit
        (committed.hash, lost.hash)
    };
    if fault == LedgerFault::CommittedPrefixLost {
        remove_ledger_files(std::path::Path::new(path));
    }
    let reopened = fs_ledger::Ledger::open(path).expect("reopen after crash");
    LedgerRecovery {
        schema_version_readable: reopened.schema_version().is_ok(),
        committed_events: reopened.table_count("events").unwrap_or(u64::MAX),
        committed_payload_readback: reopened
            .get_artifact(&committed)
            .ok()
            .flatten()
            .is_some_and(|bytes| bytes == LEDGER_COMMITTED_PAYLOAD),
        uncommitted_absent: matches!(reopened.get_artifact(&lost), Ok(None)),
        artifact_integrity_clean: reopened
            .verify_artifact_integrity()
            .is_ok_and(|report| report.is_clean()),
    }
}

/// Measured outcome of the model-form escalation drill.
#[derive(Debug, Clone)]
struct ModelFormEscalation {
    /// The fitted conformal band half-width (the measured state the
    /// decision is taken on).
    band_half_width: f64,
    /// A decision tolerance TIGHTER than the band.
    tight_tolerance: f64,
    /// A decision tolerance the band satisfies.
    loose_tolerance: f64,
    /// `certify_or_escalate` escalated at the tight tolerance.
    escalated_when_band_too_wide: bool,
    /// The escalated query was actually served by the funded lane.
    served_by_high_fidelity: bool,
    /// …and the policy still serves the surrogate when the band fits.
    serves_surrogate_when_band_fits: bool,
    /// Campaign candidates whose measured L/D fell outside their band.
    conformal_violations: usize,
    /// The reason string the escalation carried.
    escalation_reason: String,
}

/// Take a REAL certify-or-escalate decision on the fitted band and act
/// on it: an escalated query is served by `serve` (the funded lane), not
/// by the surrogate.
fn run_model_form_escalation_drill(
    surrogate: &fs_ornith::LdSurrogate,
    campaign: &[OrnithCandidate],
    serve: &mut dyn FnMut(&OrnithCandidate) -> Option<f64>,
) -> ModelFormEscalation {
    let band_half_width = surrogate.band.half_width;
    let tight_tolerance = band_half_width * 0.5;
    let loose_tolerance = band_half_width * 2.0;
    let tight = certify_or_escalate(&surrogate.band, true, tight_tolerance);
    let loose = certify_or_escalate(&surrogate.band, true, loose_tolerance);
    let escalated_when_band_too_wide = matches!(tight, Decision::Escalate { .. });
    let escalation_reason = match &tight {
        Decision::Escalate { reason } => reason.clone(),
        Decision::UseSurrogate { .. } => "not-escalated".to_string(),
    };
    let served_by_high_fidelity = escalated_when_band_too_wide
        && campaign
            .first()
            .is_some_and(|c| serve(c).is_some_and(f64::is_finite));
    ModelFormEscalation {
        band_half_width,
        tight_tolerance,
        loose_tolerance,
        escalated_when_band_too_wide,
        served_by_high_fidelity,
        serves_surrogate_when_band_fits: matches!(loose, Decision::UseSurrogate { .. }),
        conformal_violations: campaign
            .iter()
            .filter(|c| !surrogate.band.covers(surrogate.predict(c), lift_to_drag(c)))
            .count(),
        escalation_reason,
    }
}

/// The gate fe2e-007 applies to the escalation drill.
fn model_form_escalation_gate(out: &ModelFormEscalation) -> bool {
    out.escalated_when_band_too_wide
        && out.served_by_high_fidelity
        && out.serves_surrogate_when_band_fits
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
        let mut h = (i as u64) << 32 ^ t ^ CANCELLATION_INPUT_SEED;
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
    forensic_row(
        "drill",
        "cancellation-storm",
        format!(
            "{{\"winner\":{},\"survivors\":{},\
             \"input_seed\":{CANCELLATION_INPUT_SEED}}}",
            out.winner,
            out.survivors.len()
        ),
    );

    // (b) BUDGET EXHAUSTION: a REAL budget counter. The campaign funds
    // `LBM_REFINE_BUDGET` high-fidelity refinements; once the counter is
    // spent, every remaining candidate DEGRADES to the surrogate+conformal
    // path. Both the funded and the degraded counts are measured here — the
    // row used to carry a hardcoded "degraded":6 with no budget, no counter
    // and no funded call anywhere in the drill (bead
    // `frankensim-extreal-program-f85xj.2.32`).
    let (sur, campaign) = surrogate_and_campaign_fixture();
    // The funded lane is the ornithoid's real LBM refinement.
    let mut fund = |c: &OrnithCandidate| {
        let report = fs_ornith::refine(c);
        (report.lift.is_finite() && report.drag.is_finite() && report.steadiness.is_finite())
            .then_some(report.lift)
    };
    let exhaustion = run_budget_exhaustion_drill(&sur, &campaign, LBM_REFINE_BUDGET, &mut fund);
    // The drill is only evidence if the budget actually ran out mid-campaign
    // (funded lane exercised AND degradation observed), and if the degraded
    // estimates stay inside their conformal band.
    let budget_ok = budget_exhaustion_gate(&exhaustion);
    let in_band = exhaustion.in_band;
    forensic_row(
        "drill",
        "budget-exhaustion",
        format!(
            "{{\"lbm_refine_budget\":{LBM_REFINE_BUDGET},\"candidates\":{CAMPAIGN_CANDIDATES},\
             \"funded\":{},\"degraded\":{},\"in_band\":{in_band},\"funded_lift_finite\":{},\
             \"input_seed\":{SURROGATE_INPUT_SEED}}}",
            exhaustion.funded, exhaustion.degraded, exhaustion.funded_lift_finite,
        ),
    );

    // (c) LEDGER CRASH-RECOVERY: commit a prefix, begin a second
    // transaction, drop the handle WITHOUT committing (the crash), reopen —
    // and READ THE COMMITTED PREFIX BACK. `schema_version().is_ok()` alone
    // could not tell a recovered prefix from a lost one (bead
    // `frankensim-extreal-program-f85xj.2.33`).
    let dir = std::env::temp_dir().join(format!("fe2e-ledger-{}", std::process::id()));
    let _ = std::fs::create_dir_all(&dir);
    let path = dir.join("crash.db");
    let path_s = path.to_str().expect("utf8 path").to_owned();
    let recovery = run_ledger_crash_recovery_drill(&path_s, LedgerFault::None);
    let recovered = recovery.recovered();
    forensic_row(
        "drill",
        "ledger-crash-recovery",
        format!(
            "{{\"schema_version_readable\":{},\"committed_events\":{},\
             \"committed_payload_readback\":{},\"uncommitted_absent\":{},\
             \"artifact_integrity_clean\":{},\"recovered\":{recovered},\
             \"input_seed\":{FIXED_INPUT_SEED}}}",
            recovery.schema_version_readable,
            recovery.committed_events,
            recovery.committed_payload_readback,
            recovery.uncommitted_absent,
            recovery.artifact_integrity_clean,
        ),
    );

    // (d) MODEL-FORM ESCALATION: a REAL certify-or-escalate decision on
    // MEASURED state — the fitted conformal band's half-width against the
    // caller's decision tolerance — and the escalation is then SPENT on the
    // funded high-fidelity lane. The row used to assert
    // "escalates_on_outlier":true from a single `!covers(1000.0, …)`
    // predicate with no escalation machinery anywhere in the drill (bead
    // `frankensim-extreal-program-f85xj.2.32`).
    let mut serve = |c: &OrnithCandidate| {
        let report = fs_ornith::refine(c);
        (report.lift.is_finite() && report.drag.is_finite() && report.steadiness.is_finite())
            .then_some(report.lift)
    };
    let escalation = run_model_form_escalation_drill(&sur, &campaign, &mut serve);
    let escalation_ok = model_form_escalation_gate(&escalation);
    forensic_row(
        "drill",
        "model-form-escalation",
        format!(
            "{{\"band_half_width\":{:e},\"tight_tolerance\":{:e},\"loose_tolerance\":{:e},\
             \"escalated_when_band_too_wide\":{},\"served_by_high_fidelity\":{},\
             \"serves_surrogate_when_band_fits\":{},\"conformal_violations\":{},\"of\":{},\
             \"escalation_reason\":{},\"input_seed\":{SURROGATE_INPUT_SEED}}}",
            escalation.band_half_width,
            escalation.tight_tolerance,
            escalation.loose_tolerance,
            escalation.escalated_when_band_too_wide,
            escalation.served_by_high_fidelity,
            escalation.serves_surrogate_when_band_fits,
            escalation.conformal_violations,
            CAMPAIGN_CANDIDATES,
            json_string(&escalation.escalation_reason),
        ),
    );

    verdict(
        "fe2e-007-failure-drills",
        storm_ok && budget_ok && recovered && escalation_ok,
        &format!(
            "cancellation storm (input seed {CANCELLATION_INPUT_SEED:#x}): winner {} \
             among survivors; budget exhaustion (input seed {SURROGATE_INPUT_SEED:#x}): \
             {} funded LBM refinements then {} candidates degraded to surrogate+conformal, \
             {in_band}/{} degraded estimates inside the band; ledger crash-recovery: \
             reopened with {} committed event(s), committed payload read back byte-exact, \
             uncommitted write absent, artifact integrity clean; model-form escalation: \
             certify_or_escalate escalated at tolerance {:e} < band {:e} and the query was \
             served by the funded LBM lane, while tolerance {:e} kept the surrogate; \
             composite aggregate seed is zero",
            out.winner,
            exhaustion.funded,
            exhaustion.degraded,
            exhaustion.degraded,
            recovery.committed_events,
            escalation.tight_tolerance,
            escalation.band_half_width,
            escalation.loose_tolerance,
        ),
        FIXED_INPUT_SEED,
    );
}

/// REGRESSION for bead `frankensim-extreal-program-f85xj.2.32` (b): the
/// budget-exhaustion row must be MEASURED. The pre-fix drill had no
/// budget variable, no counter and no funded call — its `"degraded":6`
/// was a literal in the format string, so seeding the bead's fault
/// ("make degradation impossible") left the row and the verdict
/// unchanged.
///
/// The funded lane is stubbed here so the ACCOUNTING can be falsified
/// without paying for six LBM refinements; fe2e-007 itself spends the
/// real one.
#[test]
fn fe2e_007_budget_drill_counts_real_degradation() {
    let (sur, campaign) = surrogate_and_campaign_fixture();
    let calls = std::cell::Cell::new(0usize);
    let mut stub = |c: &OrnithCandidate| {
        calls.set(calls.get() + 1);
        Some(lift_to_drag(c))
    };

    // The shipped configuration: the counter runs out mid-campaign.
    let real = run_budget_exhaustion_drill(&sur, &campaign, LBM_REFINE_BUDGET, &mut stub);
    assert_eq!(real.funded, LBM_REFINE_BUDGET);
    assert_eq!(real.degraded, CAMPAIGN_CANDIDATES - LBM_REFINE_BUDGET);
    assert_eq!(
        calls.get(),
        LBM_REFINE_BUDGET,
        "the funded lane must be spent"
    );
    assert!(budget_exhaustion_gate(&real), "{real:?}");

    // SEEDED FAULT: fund every candidate, so nothing can degrade. A drill
    // that reports its outcome must now say `degraded == 0` and FAIL its
    // gate — the pre-fix row would still have read "degraded":6.
    calls.set(0);
    let unlimited = run_budget_exhaustion_drill(&sur, &campaign, CAMPAIGN_CANDIDATES, &mut stub);
    assert_eq!(unlimited.funded, CAMPAIGN_CANDIDATES);
    assert_eq!(unlimited.degraded, 0);
    assert_eq!(calls.get(), CAMPAIGN_CANDIDATES);
    assert!(
        !budget_exhaustion_gate(&unlimited),
        "no degradation must not pass the budget-exhaustion gate: {unlimited:?}"
    );

    // SEEDED FAULT: the funded lane cannot answer. The drill must not
    // report a healthy funded lane.
    let mut dead = |_: &OrnithCandidate| None;
    let broken = run_budget_exhaustion_drill(&sur, &campaign, LBM_REFINE_BUDGET, &mut dead);
    assert!(!broken.funded_lift_finite);
    assert!(!budget_exhaustion_gate(&broken), "{broken:?}");
}

/// REGRESSION for bead `frankensim-extreal-program-f85xj.2.33`: the
/// ledger crash-recovery drill must read the committed prefix BACK. The
/// pre-fix gate was `reopened.schema_version().is_ok()`, which stays
/// green under a recovery that discards the committed prefix AND under
/// one that replays the uncommitted transaction — both asserted here.
#[test]
fn fe2e_007_ledger_drill_detects_a_broken_recovery() {
    let dir = std::env::temp_dir().join(format!("fe2e-ledger-fault-{}", std::process::id()));
    let _ = std::fs::create_dir_all(&dir);
    let path = |tag: &str| {
        dir.join(format!("{tag}.db"))
            .to_str()
            .expect("utf8 path")
            .to_owned()
    };

    let healthy = run_ledger_crash_recovery_drill(&path("healthy"), LedgerFault::None);
    assert!(healthy.recovered(), "{healthy:?}");
    assert_eq!(healthy.committed_events, 1);
    assert!(healthy.committed_payload_readback);
    assert!(healthy.uncommitted_absent);

    let lost = run_ledger_crash_recovery_drill(&path("lost"), LedgerFault::CommittedPrefixLost);
    assert!(
        lost.schema_version_readable,
        "the pre-fix gate reads Ok even here — that is the defect"
    );
    assert_eq!(lost.committed_events, 0);
    assert!(!lost.committed_payload_readback);
    assert!(
        !lost.recovered(),
        "a lost committed prefix must fail: {lost:?}"
    );

    let replayed =
        run_ledger_crash_recovery_drill(&path("replayed"), LedgerFault::UncommittedReplayed);
    assert!(
        replayed.schema_version_readable,
        "the pre-fix gate reads Ok even here — that is the defect"
    );
    assert_eq!(replayed.committed_events, 2);
    assert!(
        !replayed.uncommitted_absent,
        "the uncommitted artifact survived the replay fault"
    );
    assert!(
        !replayed.recovered(),
        "a replayed uncommitted transaction must fail: {replayed:?}"
    );
}

/// REGRESSION for bead `frankensim-extreal-program-f85xj.2.32` (d): the
/// escalation row must come from a real decision AND a real serve. The
/// pre-fix row asserted `"escalates_on_outlier":true` from a single
/// `!band.covers(1000.0, …)` predicate — no policy, no escalation, no
/// high-fidelity call.
#[test]
fn fe2e_007_escalation_drill_takes_a_real_decision_and_spends_it() {
    let (sur, campaign) = surrogate_and_campaign_fixture();
    let calls = std::cell::Cell::new(0usize);
    let mut stub = |c: &OrnithCandidate| {
        calls.set(calls.get() + 1);
        Some(lift_to_drag(c))
    };
    let out = run_model_form_escalation_drill(&sur, &campaign, &mut stub);
    assert!(sur.band.half_width.is_finite() && sur.band.half_width > 0.0);
    assert!(out.escalated_when_band_too_wide, "{out:?}");
    assert!(
        out.escalation_reason.contains("decision tolerance")
            || out.escalation_reason.contains("tolerance")
            || !out.escalation_reason.is_empty(),
        "escalation must carry a reason: {out:?}"
    );
    assert_eq!(
        calls.get(),
        1,
        "the escalated query must reach the funded lane"
    );
    assert!(out.served_by_high_fidelity);
    assert!(out.serves_surrogate_when_band_fits);
    assert!(model_form_escalation_gate(&out));

    // SEEDED FAULT: the funded lane cannot answer the escalated query.
    // The drill must not report the escalation as served.
    let mut dead = |_: &OrnithCandidate| None;
    let unserved = run_model_form_escalation_drill(&sur, &campaign, &mut dead);
    assert!(unserved.escalated_when_band_too_wide);
    assert!(!unserved.served_by_high_fidelity);
    assert!(!model_form_escalation_gate(&unserved), "{unserved:?}");
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
    forensic_row(
        "notebook",
        "emitted",
        format!(
            "{{\"bytes\":{},\"aggregate_input_seed\":{FIXED_INPUT_SEED},\
             \"vessel_input_seed\":{FIXED_INPUT_SEED},\
             \"ornith_input_seed\":{ORNITH_INPUT_SEED},\
             \"frame_input_seed\":{FRAME_INPUT_SEED},\
             \"frame_execution_seed\":{FRAME_EXECUTION_SEED}}}",
            n1.len()
        ),
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
            "constructed log-row fixtures escape string fields and carry required keys; \
             lab notebook ({} bytes) regenerates BITWISE on full replay using vessel input \
             seed zero, ornith input seed {ORNITH_INPUT_SEED:#x}, frame input seed \
             {FRAME_INPUT_SEED}, and frame execution seed {FRAME_EXECUTION_SEED:#x}; \
             composite aggregate seed is zero",
            n1.len()
        ),
        FIXED_INPUT_SEED,
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
