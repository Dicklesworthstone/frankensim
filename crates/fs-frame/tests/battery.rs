//! fs-frame conformance battery (bead mye.3, smoke tier): layout
//! diagnostics, sizing code rows, fiber-hinge hysteresis, e-stopped
//! fragility with verified coverage and ledgered savings, CVaR mass
//! minimization, and the replay/drill gates.

use fs_evidence::Color;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_frame::cvar::{RobustError, cvar, empirical_cvar};
use fs_frame::history::{StoryFrame, StoryParams, peak_drift};
use fs_frame::{LayoutError, cvar_mass_min, e_stopped_fragility, ensemble_cvar, layout_and_size};
use fs_qty::{Dims, QtyAny};
use fs_scenario::ensemble::{SpectrumModel, StochasticEnsemble};
use fs_truss::ground::TrussConstructionError;

fn verdict(name: &str, pass: bool, details: &str) {
    println!("{{\"test\":\"{name}\",\"pass\":{pass},\"details\":\"{details}\"}}");
    assert!(pass, "{name}: {details}");
}

const TIME: Dims = Dims([0, 0, 1, 0, 0, 0]);
const RATE: Dims = Dims([0, 0, -1, 0, 0, 0]);

fn with_cx<R>(cancelled: bool, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    if cancelled {
        gate.request();
    }
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 0xF2A4_E001,
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

fn kt_ensemble(members: u32, s0: f64, seed: u64) -> StochasticEnsemble {
    StochasticEnsemble {
        name: "kt-suite".to_string(),
        seed,
        members,
        duration: QtyAny::new(12.0, TIME),
        dt: QtyAny::new(0.02, TIME),
        model: SpectrumModel::KanaiTajimi {
            s0,
            omega_g: QtyAny::new(12.5, RATE),
            zeta_g: 0.6,
        },
    }
}

/// frame-001/002: layout LP diagnostics + sizing
/// audit — the fs-truss composition run end-to-end on the flagship's
/// cantilever fixture. The Michell catalogue row stays LEDGERED
/// PENDING (fs-truss contract).
#[test]
fn frame_001_layout_and_sizing() {
    let catalog: Vec<f64> = (1..=20).map(|k| 2e-4 * f64::from(k)).collect();
    let report = with_cx(false, |cx| {
        layout_and_size(5, 3, 4.0, 2.0, 250e6, 200e9, &catalog, cx)
            .expect("valid smoke-tier frame layout is admitted")
    });
    verdict(
        "frame-001-lp-diagnostics",
        report.gap < 1e-3 && report.residual < 1e-3 && report.volume > 0.0,
        &format!(
            "objective separation {:.2e}, eq residual {:.2e}, approximate volume {:.4e}",
            report.gap, report.residual, report.volume
        ),
    );
    let Color::Verified { lo, hi } = report.optimality_color else {
        panic!("frame layout optimum must carry the fs-truss certificate");
    };
    verdict(
        "frame-001-outward-optimum-bounds",
        lo.is_finite() && hi.is_finite() && lo > 0.0 && lo <= hi,
        &format!("physical optimum volume in outward interval [{lo:.8e}, {hi:.8e}]"),
    );
    verdict(
        "frame-002-sizing-code-rows",
        report.audit.all_pass
            && report.audit.eq_residual < 1e-6
            && !report.audit.members.is_empty(),
        &format!(
            "{} members sized (pruned {}), post-prune eq residual {:.2e}, all code rows pass",
            report.audit.members.len(),
            report.audit.pruned,
            report.audit.eq_residual
        ),
    );
}

/// G4: a pre-cancelled construction scope refuses the frame before publishing
/// a partial ground structure or LP.
#[test]
fn frame_002_pre_cancelled_layout_is_refused() {
    let catalog: Vec<f64> = (1..=20).map(|k| 2e-4 * f64::from(k)).collect();
    let result = with_cx(true, |cx| {
        layout_and_size(5, 3, 4.0, 2.0, 250e6, 200e9, &catalog, cx)
    });
    verdict(
        "frame-002-pre-cancel-refusal",
        matches!(
            result,
            Err(LayoutError::Construction(
                TrussConstructionError::Cancelled { .. }
            ))
        ),
        "pre-cancelled construction returns a structured refusal",
    );
}

/// frame-003: the fiber-hinge story model — elastic runs conserve
/// sanity (no drift ratcheting, Newmark stable over 10× duration) and
/// a yielding cyclic run DISSIPATES energy through the true
/// Mander/Menegotto–Pinto fibers (hysteresis with area, not a spring).
#[test]
fn frame_003_hysteresis() {
    let params = StoryParams::default();
    // Elastic: tiny sinusoidal shaking, long duration.
    let dt = 0.02;
    let n = 6000usize; // 120 s — 10× the study duration
    let ag_small: Vec<f64> = (0..u32::try_from(n).expect("small"))
        .map(|i| 0.02 * (0.8 * f64::from(i) * dt * std::f64::consts::TAU).sin())
        .collect();
    let mut frame = StoryFrame::new(params);
    let drifts = frame.run(&ag_small, dt);
    let early = peak_drift(&drifts[..n / 10], params.h);
    let late = peak_drift(&drifts[n - n / 10..], params.h);
    verdict(
        "frame-003-elastic-stability",
        late < 3.0 * early.max(1e-12) && drifts.iter().all(|d| d.is_finite()),
        &format!("peak drift early {early:.3e} late {late:.3e} over 120 s (no ratcheting)"),
    );
    // Yielding cycle: strong shaking; energy dissipated = ∮ V dx > 0.
    let ag_big: Vec<f64> = (0..600)
        .map(|i| 3.0 * (1.2 * f64::from(i) * dt * std::f64::consts::TAU).sin())
        .collect();
    let mut frame2 = StoryFrame::new(params);
    let mut dissipated = 0.0f64;
    let mut x_prev = 0.0f64;
    let mut peak = 0.0f64;
    for chunk in ag_big.chunks(1) {
        let d = frame2.run(chunk, dt);
        let x = d[0];
        let (v, _) = frame2.restoring(x);
        dissipated += v * (x - x_prev);
        x_prev = x;
        peak = peak.max((x / params.h).abs());
    }
    verdict(
        "frame-003-hysteretic-dissipation",
        dissipated > 0.0 && peak > 0.002,
        &format!("cyclic work {dissipated:.4e} > 0 at peak drift ratio {peak:.4}"),
    );
}

/// frame-004: e-stopped fragility — the study stops itself when the
/// confidence sequence is decision-grade; the interval at the STOP
/// time covers the fixed-N reference (anytime validity in action);
/// the savings vs fixed-N are measured and ledgered.
#[test]
fn frame_004_e_stopped_fragility() {
    // s0 tuned so the 2e-2 threshold DISCRIMINATES (at 0.05 every
    // member exceeded and the coverage gate was toothless).
    let ensemble = kt_ensemble(200, 0.01, 90210);
    let params = StoryParams::default();
    // Reference exceedance from the FULL fixed-N suite.
    let dt = ensemble.dt.value;
    let mut exceed = 0u32;
    // The plan's own Appendix C threshold: exceeds(peak-drift, 2e-2).
    let limit = 2e-2;
    for member in 0..ensemble.members {
        let real = ensemble.realize(member).expect("realizes");
        let mut frame = StoryFrame::new(params);
        let drifts = frame.run(&real.values, dt);
        if peak_drift(&drifts, params.h) > limit {
            exceed += 1;
        }
    }
    let p_ref = f64::from(exceed) / f64::from(ensemble.members);
    let report = e_stopped_fragility(&ensemble, params, limit, 0.05, 0.12);
    let covered = (report.p_hat - p_ref).abs() <= report.radius;
    verdict(
        "frame-004-e-stop-coverage",
        covered && report.members_used <= 200 && exceed > 0 && exceed < 200,
        &format!(
            "p_ref {p_ref:.3} in CS [{:.3} ± {:.3}] at stop after {} members ({} exceedances)",
            report.p_hat, report.radius, report.members_used, report.exceedances
        ),
    );
    verdict(
        "frame-004-e-stop-savings",
        report.stopped_early && report.members_used < 200,
        &format!(
            "LEDGER: e-stop consumed {}/200 members ({}% saving vs fixed-N); MLMC levels {}",
            report.members_used,
            100 * (200 - report.members_used) / 200,
            report.mlmc.levels.len()
        ),
    );
}

/// frame-005: CVaR-constrained mass minimization — CVaR decreases
/// monotonically in the section scale (spot check), the bisection
/// lands on a feasible minimal design, and the catalog snap preserves
/// feasibility under the independent re-check.
#[test]
fn frame_005_cvar_mass_min() {
    let ensemble = kt_ensemble(48, 0.05, 777);
    let params = StoryParams::default();
    let catalog = [0.5f64, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0];
    // Bracket the limit between the catalog extremes: this makes the
    // bisection MEANINGFUL whatever the absolute drift scale, and the
    // monotonicity gate is the physics check (bigger sections, less
    // tail drift).
    let c_weak = ensemble_cvar(&ensemble, params, 0.25, 0.9);
    let c_strong = ensemble_cvar(&ensemble, params, 4.0, 0.9);
    verdict(
        "frame-005-cvar-monotone",
        c_weak > c_strong && c_strong > 0.0,
        &format!("CVaR90: scale 0.25 -> {c_weak:.4}, scale 4.0 -> {c_strong:.4}"),
    );
    let limit = (c_weak * c_strong).sqrt();
    let design = cvar_mass_min(&ensemble, params, 0.9, limit, &catalog);
    verdict(
        "frame-005-cvar-design",
        design.cvar_snapped <= limit
            && design.scale_snapped >= design.scale_star
            && catalog.contains(&design.scale_snapped)
            && design.iters > 0,
        &format!(
            "limit {limit:.4}: scale* {:.3} -> snapped {:.3}; CVaR {:.4} -> {:.4}; mass {:.2}; {} bisections",
            design.scale_star,
            design.scale_snapped,
            design.cvar_star,
            design.cvar_snapped,
            design.mass,
            design.iters
        ),
    );
}

/// frame-006: replay determinism + the drills — a bitwise-identical
/// rerun, the budget-exhaustion path (honest indecision, no early
/// stop), and the infeasible-constraint diagnostic firing.
#[test]
fn frame_006_replay_and_drills() {
    let ensemble = kt_ensemble(40, 0.05, 4242);
    let params = StoryParams::default();
    let a = e_stopped_fragility(&ensemble, params, 2e-2, 0.05, 0.03);
    let b = e_stopped_fragility(&ensemble, params, 2e-2, 0.05, 0.03);
    verdict(
        "frame-006-replay-determinism",
        a.p_hat.to_bits() == b.p_hat.to_bits()
            && a.members_used == b.members_used
            && a.radius.to_bits() == b.radius.to_bits(),
        &format!(
            "two runs bitwise identical: p {:.4}, {} members",
            a.p_hat, a.members_used
        ),
    );
    // Budget exhaustion: a margin no smoke budget can reach.
    let tiny = kt_ensemble(12, 0.05, 5150);
    let exhausted = e_stopped_fragility(&tiny, params, 2e-2, 0.05, 1e-4);
    verdict(
        "frame-006-budget-exhaustion-drill",
        !exhausted.stopped_early && exhausted.members_used == 12 && exhausted.radius > 1e-4,
        &format!(
            "budget exhausted honestly: {}/12 members, radius {:.3} (indecision REPORTED)",
            exhausted.members_used, exhausted.radius
        ),
    );
    // Infeasible CVaR study: limit unreachable even at the max scale.
    let infeasible = std::panic::catch_unwind(|| {
        let ens = kt_ensemble(12, 0.05, 31337);
        cvar_mass_min(&ens, params, 0.9, 1e-9, &[0.5, 1.0])
    });
    verdict(
        "frame-006-infeasible-drill",
        infeasible.is_err(),
        "infeasible CVaR limit fires the diagnostic instead of returning a fake design",
    );
    let extreme_samples = [-f64::MAX, 0.0, f64::MAX];
    let frame_report = empirical_cvar(&extreme_samples, 0.25).expect("valid extreme samples");
    let canonical_report =
        fs_robust::empirical_cvar(&extreme_samples, 0.25).expect("canonical extreme samples");
    verdict(
        "frame-006-canonical-cvar-parity",
        frame_report == canonical_report
            && cvar(&extreme_samples, 0.25)
                .is_ok_and(|value| value.to_bits() == canonical_report.cvar().to_bits()),
        "frame report and scalar CVaR surfaces are exact canonical fs-robust re-exports",
    );
    let empty_cvar = empirical_cvar(&[], 0.9);
    verdict(
        "frame-006-empty-cvar-drill",
        matches!(empty_cvar, Err(RobustError::EmptySamples)),
        "empty CVaR samples return a structured refusal instead of fake zero risk",
    );
    let bad_beta = empirical_cvar(&[1.0, 2.0], 1.0);
    verdict(
        "frame-006-bad-beta-drill",
        matches!(bad_beta, Err(RobustError::BadAlpha { alpha }) if alpha.to_bits() == 1.0f64.to_bits()),
        "invalid CVaR beta returns a structured refusal before quantile indexing",
    );
    let nan_beta = empirical_cvar(&[1.0, 2.0], f64::NAN);
    verdict(
        "frame-006-nan-beta-drill",
        matches!(nan_beta, Err(RobustError::BadAlpha { alpha }) if alpha.is_nan()),
        "non-finite CVaR beta returns a structured refusal before quantile indexing",
    );
    let bad_loss = empirical_cvar(&[1.0, f64::NAN], 0.9);
    verdict(
        "frame-006-nonfinite-cvar-drill",
        matches!(bad_loss, Err(RobustError::BadSample { value }) if value.is_nan()),
        "non-finite CVaR losses return a structured refusal before tail aggregation",
    );
}
