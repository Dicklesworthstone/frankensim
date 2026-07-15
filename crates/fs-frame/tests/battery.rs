//! fs-frame conformance battery (bead mye.3, smoke tier): layout
//! diagnostics, sizing code rows, fiber-hinge hysteresis, e-stopped
//! fragility with verified coverage and ledgered savings, CVaR mass
//! minimization, and the replay/drill gates.

use fs_evidence::Color;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_frame::cvar::{RobustError, cvar, empirical_cvar};
use fs_frame::history::{
    GroundMotion, HistoryError, HistoryLimits, StoryFrame, StoryParams, peak_drift,
};
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
const EL_CENTRO_AT2: &[u8] = include_bytes!("data/el_centro_1940_ns.at2");
const EL_CENTRO_PROVENANCE: &str = include_str!("data/EL_CENTRO_PROVENANCE.md");
const STANDARD_GRAVITY_M_S2: f64 = 9.806_65;

struct At2Fixture {
    samples_g: Vec<f64>,
    dt_s: f64,
}

fn parse_pinned_el_centro(text: &str) -> Result<At2Fixture, String> {
    let mut lines = text.lines();
    let title = lines.next().ok_or("missing AT2 title")?;
    let layout = lines.next().ok_or("missing AT2 layout")?;
    let units = lines.next().ok_or("missing AT2 units")?;
    let card = lines.next().ok_or("missing AT2 NPTS/DT card")?;
    if title != "Data for El Centro 1940 North South Component (Peknold Version)" {
        return Err("unexpected AT2 title".to_string());
    }
    if layout
        != "1559 points at equal spacing of 0.02 sec Points are listed in the format of 8F10.5, i.e., 8 points across in a row with 5 decimal places"
    {
        return Err("unexpected AT2 layout declaration".to_string());
    }
    if units != "The units are (g)" || card != "NPTS=  1559, DT= .02000 SEC" {
        return Err("unexpected AT2 units or NPTS/DT card".to_string());
    }

    let mut samples_g = Vec::with_capacity(1_559);
    for token in lines.flat_map(|line| line.split_ascii_whitespace()) {
        let value = token
            .parse::<f64>()
            .map_err(|_| "non-numeric AT2 sample".to_string())?;
        if !value.is_finite() {
            return Err("non-finite AT2 sample".to_string());
        }
        samples_g.push(value);
    }
    if samples_g.len() != 1_559 {
        return Err(format!(
            "AT2 sample count {} does not match 1559",
            samples_g.len()
        ));
    }
    Ok(At2Fixture {
        samples_g,
        dt_s: 0.02,
    })
}

fn fnv1a64(bytes: &[u8]) -> u64 {
    bytes.iter().fold(0xcbf2_9ce4_8422_2325, |hash, byte| {
        (hash ^ u64::from(*byte)).wrapping_mul(0x0000_0100_0000_01b3)
    })
}

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

/// G4/G5: the recorded-motion integration surface binds units and work limits,
/// retains the displacement/restoring-shear pair, reproduces the legacy Newmark
/// path bit-for-bit, and publishes no fiber state when a sample fails to
/// converge. This is the response artifact required before a cited El Centro
/// envelope can be claimed; this fixture remains synthetic.
#[test]
#[allow(clippy::too_many_lines)] // One lifecycle proves parity, rollback, and admission.
fn frame_003_checked_response_is_bounded_atomic_and_replayable() {
    let params = StoryParams::default();
    let dt_s = 0.02;
    let acceleration_m_s2: Vec<f64> = (0..240)
        .map(|sample| 0.35 * (1.1 * f64::from(sample) * dt_s * std::f64::consts::TAU).sin())
        .collect();
    let limits = HistoryLimits::new(acceleration_m_s2.len(), 30, 1e-12, 1e-2);

    let mut legacy = StoryFrame::new(params);
    let legacy_displacement = legacy.run(&acceleration_m_s2, dt_s);
    let mut checked = StoryFrame::new(params);
    let response = checked
        .run_checked(GroundMotion::new(&acceleration_m_s2, dt_s), limits)
        .expect("admitted synthetic record converges");
    let same_displacement = legacy_displacement
        .iter()
        .zip(&response.displacement_m)
        .all(|(legacy, checked)| legacy.to_bits() == checked.to_bits());
    let peak_displacement = response
        .displacement_m
        .iter()
        .fold(0.0f64, |peak, value| peak.max(value.abs()));
    let peak_shear = response
        .restoring_shear_n
        .iter()
        .fold(0.0f64, |peak, value| peak.max(value.abs()));
    verdict(
        "frame-003-checked-response",
        same_displacement
            && response.displacement_m.len() == acceleration_m_s2.len()
            && response.restoring_shear_n.len() == acceleration_m_s2.len()
            && response
                .restoring_shear_n
                .iter()
                .all(|value| value.is_finite())
            && response.peak_abs_displacement_m.to_bits() == peak_displacement.to_bits()
            && response.peak_abs_restoring_shear_n.to_bits() == peak_shear.to_bits()
            && response.max_abs_equilibrium_residual_n < limits.equilibrium_tolerance_n
            && (1..=limits.max_newton_iterations).contains(&response.max_newton_iterations_used),
        &format!(
            "{} displacement/shear pairs; peak |x| {:.4e} m, peak |V_restore| {:.4e} N, \
             max Newton corrections {}",
            response.displacement_m.len(),
            response.peak_abs_displacement_m,
            response.peak_abs_restoring_shear_n,
            response.max_newton_iterations_used
        ),
    );

    let mut refused = StoryFrame::new(params);
    let state_before_refusal = (
        refused.x.to_bits(),
        refused.v.to_bits(),
        refused.a.to_bits(),
    );
    let refusal = refused.run_checked(
        GroundMotion::new(&[1e-9, 3.0], dt_s),
        HistoryLimits::new(2, 1, 1e-12, 1e-2),
    );
    let public_state_preserved = state_before_refusal
        == (
            refused.x.to_bits(),
            refused.v.to_bits(),
            refused.a.to_bits(),
        );
    let mut pristine = StoryFrame::new(params);
    let probe = [0.1, -0.1, 0.0];
    let after_refusal = refused
        .run_checked(
            GroundMotion::new(&probe, dt_s),
            HistoryLimits::new(probe.len(), 30, 1e-12, 1e-2),
        )
        .expect("frame remains usable after refusal");
    let pristine_response = pristine
        .run_checked(
            GroundMotion::new(&probe, dt_s),
            HistoryLimits::new(probe.len(), 30, 1e-12, 1e-2),
        )
        .expect("pristine probe converges");
    let no_state_published = after_refusal
        .displacement_m
        .iter()
        .zip(&pristine_response.displacement_m)
        .all(|(after, pristine)| after.to_bits() == pristine.to_bits())
        && after_refusal
            .restoring_shear_n
            .iter()
            .zip(&pristine_response.restoring_shear_n)
            .all(|(after, pristine)| after.to_bits() == pristine.to_bits());
    verdict(
        "frame-003-checked-response-refusal",
        matches!(
            refusal,
            Err(HistoryError::NewtonDidNotConverge {
                sample: 1,
                iterations: 1,
                ..
            })
        ) && public_state_preserved
            && no_state_published,
        "one-correction budget refuses sample 1 and rolls back sample 0's staged commit",
    );

    let invalid = StoryFrame::new(params).run_checked(
        GroundMotion::new(&[0.0, f64::NAN], dt_s),
        HistoryLimits::new(2, 30, 1e-12, 1e-2),
    );
    let over_budget = StoryFrame::new(params).run_checked(
        GroundMotion::new(&[0.0], dt_s),
        HistoryLimits::new(0, 30, 1e-12, 1e-2),
    );
    let invalid_geometry = StoryFrame::new(StoryParams {
        h: f64::MAX,
        lp: 2.0,
        ..params
    })
    .run_checked(
        GroundMotion::new(&[0.0], dt_s),
        HistoryLimits::new(1, 30, 1e-12, 1e-2),
    );
    verdict(
        "frame-003-checked-response-admission",
        matches!(
            invalid,
            Err(HistoryError::NonFiniteAcceleration { sample: 1, .. })
        ) && matches!(
            over_budget,
            Err(HistoryError::SampleLimitExceeded {
                samples: 1,
                max_samples: 0
            })
        ) && matches!(
            invalid_geometry,
            Err(HistoryError::InvalidStoryParameter { .. })
        ),
        "non-finite acceleration, sample excess, and derived-geometry overflow fail admission",
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

/// G0/G3: physically different fs-scenario realization payloads must not be
/// accepted merely because they share the same `Vec<f64>` representation.
#[test]
fn frame_007_refuses_non_ground_motion_and_empty_ensembles() {
    let params = StoryParams::default();
    let mut ensemble = kt_ensemble(0, 0.05, 7);
    let empty_fragility =
        std::panic::catch_unwind(|| e_stopped_fragility(&ensemble, params, 2e-2, 0.05, 0.1));
    assert!(
        empty_fragility.is_err(),
        "a zero-member study must refuse before confidence-sequence finalization"
    );

    ensemble.members = 1;
    ensemble.model = SpectrumModel::CarreauBand {
        eta_zero: [
            QtyAny::new(1.0, Dims([-1, 1, -1, 0, 0, 0])),
            QtyAny::new(2.0, Dims([-1, 1, -1, 0, 0, 0])),
        ],
        eta_inf: [
            QtyAny::new(0.1, Dims([-1, 1, -1, 0, 0, 0])),
            QtyAny::new(0.2, Dims([-1, 1, -1, 0, 0, 0])),
        ],
        lambda: [QtyAny::new(0.5, TIME), QtyAny::new(1.0, TIME)],
        n: [0.3, 0.8],
    };
    let material_as_motion =
        std::panic::catch_unwind(|| ensemble_cvar(&ensemble, params, 1.0, 0.9));
    assert!(
        material_as_motion.is_err(),
        "Carreau viscosity parameters must never become a structural acceleration history"
    );
}

/// G2 input slice: pin and lint the redistributable PEER/OpenSees Peknold
/// record, then drive it through the bounded frame path. The committed
/// OpenSees roof-displacement envelope remains diagnostic because its
/// distributed-column model is materially different from this concentrated-
/// hinge fixture; no public total support-reaction oracle exists there.
#[test]
#[allow(clippy::too_many_lines)] // Integrity, real-motion response, and authority boundary.
fn frame_008_el_centro_record_and_response_diagnostic() {
    let text = std::str::from_utf8(EL_CENTRO_AT2).expect("pinned AT2 is UTF-8");
    let fixture = parse_pinned_el_centro(text).expect("pinned AT2 contract parses");
    let peak_g = fixture
        .samples_g
        .iter()
        .fold(0.0f64, |peak, sample| peak.max(sample.abs()));
    let provenance_linted = [
        "fd3958fde8a3d0a350321b3fbdd3f415ee16e2a2",
        "e096f0458ae7565ac8e9fac80eda018e4715a357ea4c14504ce1cfa6556c53a6",
        "e9ae5a8c2163c28a2e52f29f8b363e2f19a3eac400e61ff73c547f8eb309b6ab",
        "b5b6e08cec39d6b1826a9d27778c3aa1f76adfb115d7ba66fb6eaf1bc7c6dd02",
        "b1f6f933b71e543e85deea342aca54d3389e0435fd6bccdf9e0da595490efdef",
        "Committed peak roof displacement: `2.47158 in = 0.062778132 m`",
        "public total-base-reaction oracle is unavailable",
        "BSD 3-Clause License",
    ]
    .iter()
    .all(|required| EL_CENTRO_PROVENANCE.contains(required));
    verdict(
        "frame-008-el-centro-record-integrity",
        EL_CENTRO_AT2.len() == 16_032
            && fnv1a64(EL_CENTRO_AT2) == 0x9eda_012e_e82c_b084
            && fixture.samples_g.len() == 1_559
            && fixture.dt_s.to_bits() == 0.02f64.to_bits()
            && peak_g.to_bits() == 0.31882f64.to_bits()
            && provenance_linted,
        "PEER/OpenSees Peknold record: 1559 samples, dt 0.02 s, PGA 0.31882 g, provenance/license pinned",
    );

    // GroundMotion uses step-end samples. This diagnostic explicitly designates
    // source row 0 as its initial forcing, then advances from index 1. That
    // timing convention is part of this artifact and is not OpenSees parity.
    let acceleration_m_s2: Vec<f64> = fixture.samples_g[1..]
        .iter()
        .map(|sample_g| sample_g * STANDARD_GRAVITY_M_S2)
        .collect();
    let params = StoryParams::default();
    let mut frame = StoryFrame::new(params);
    let limits = HistoryLimits::new(acceleration_m_s2.len(), 100, 1e-10, 1.0);
    let response = frame
        .run_checked(GroundMotion::new(&acceleration_m_s2, fixture.dt_s), limits)
        .expect("pinned record converges within its explicit response budget");
    verdict(
        "frame-008-el-centro-response-artifact",
        response.displacement_m.len() == 1_558
            && response.restoring_shear_n.len() == 1_558
            && response.peak_abs_displacement_m.is_finite()
            && response.peak_abs_displacement_m > 0.0
            && response.peak_abs_restoring_shear_n.is_finite()
            && response.peak_abs_restoring_shear_n > 0.0
            && response.max_abs_equilibrium_residual_n < limits.equilibrium_tolerance_n
            && response.max_newton_iterations_used <= limits.max_newton_iterations,
        &format!(
            "diagnostic-only concentrated-hinge response: peak |x| {:.8e} m, \
             peak |V_restore| {:.8e} N, residual {:.3e} N, max Newton {}",
            response.peak_abs_displacement_m,
            response.peak_abs_restoring_shear_n,
            response.max_abs_equilibrium_residual_n,
            response.max_newton_iterations_used
        ),
    );

    // The five-decimal OpenSees output defines a rounding envelope around the
    // external model's peak. It is pinned and linted, never applied to the
    // non-equivalent FrankenSim model.
    let external_peak_m = 2.47158 * 0.0254;
    let external_lo_m = 2.471575 * 0.0254;
    let external_hi_m = 2.471585 * 0.0254;
    let external_drift_ratio = external_peak_m / (144.0 * 0.0254);
    verdict(
        "frame-008-el-centro-reference-envelope-lint",
        external_lo_m < external_peak_m
            && external_peak_m < external_hi_m
            && (external_peak_m - 0.062_778_132).abs() < 1e-15
            && (external_drift_ratio - 0.017_163_75).abs() < 1e-15,
        "OpenSees Example 3 displacement rounding envelope pinned; model equivalence=false; support-reaction oracle=unavailable",
    );
    println!(
        "{{\"test\":\"frame-008-el-centro-model-comparison\",\"pass\":null,\
         \"authority\":\"diagnostic_only\",\"model_equivalence\":false,\
         \"frankensim_peak_displacement_m\":{:.9e},\
         \"opensees_peak_displacement_m\":{external_peak_m:.9e},\
         \"published_total_base_reaction\":null}}",
        response.peak_abs_displacement_m
    );
}
