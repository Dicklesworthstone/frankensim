//! fs-vessel conformance battery (bead mye.4, smoke tier): the
//! Orr–Sommerfeld objective path validated against the Orszag
//! reference, tilt-schedule pours with strict mass ledgers, the
//! contact-line bracketing band as a first-class output, the Carreau
//! band sweep, CVaR-vs-nominal robustification, the e-raced candidate
//! screen, and the same-bytes render deliverable.

use fs_lbm::freesurface::ContactModel;
use fs_lbm::rheology::Rheology;
use fs_vessel::pour::{PourRig, render_pour, run_pour};
use fs_vessel::race::{ScreenError, race_base_losses, screen_lips, screening_losses};
use fs_vessel::robust::{RobustError, cvar, empirical_cvar, robustify};
use fs_vessel::stability::{VesselProfile, growth_objective};

/// The battery's screening seed (the constant the inlined convention
/// hashed with before the wrapper was extracted).
const VESSEL_SCREEN_SEED: u64 = 0x7E55;

fn verdict(name: &str, pass: bool, details: &str) {
    println!("{{\"test\":\"{name}\",\"pass\":{pass},\"details\":\"{details}\"}}");
    assert!(pass, "{name}: {details}");
}

fn newtonian() -> Rheology {
    Rheology::Newtonian { nu: 0.0167 }
}

/// vsl-001: the objective path reproduces the physics it is built on —
/// growth responds to the vessel knobs along the thickness flank of
/// the U-shaped film-Re proxy (wide lip → thicker film → higher Re →
/// worse growth than the moderate lip), and the plane-Poiseuille
/// machinery underneath crosses instability between Re 5000 and 6000
/// (the Orszag bracket, re-gated through this crate's call path).
#[test]
fn vsl_001_stability_objective() {
    let narrow = growth_objective(&VesselProfile::carafe(0.6), 1.0, 1.0, 4, 4);
    let wide = growth_objective(&VesselProfile::carafe(2.4), 1.0, 1.0, 4, 4);
    // Direct Orszag bracket through the same dependency.
    let below = fs_cheb::orr_sommerfeld::max_growth(5000.0, 1.020_56, 32).expect("eig");
    let above = fs_cheb::orr_sommerfeld::max_growth(6000.0, 1.020_56, 32).expect("eig");
    verdict(
        "vsl-001-stability-objective",
        narrow < wide && below < 0.0 && above > 0.0,
        &format!(
            "growth: narrow lip {narrow:.5} < wide lip {wide:.5}; Orszag bracket through the objective's dependency: sigma(5000) = {below:.2e} < 0 < sigma(6000) = {above:.2e}"
        ),
    );
}

/// vsl-002: the pour POURS and the ledger is strict — mass crosses
/// the lip under the tilt schedule while total tracked mass drifts
/// below 1e-10 relative at every step (the make-or-break audit,
/// through the rotating-gravity moving frame).
#[test]
fn vsl_002_tilt_pour_mass() {
    let rig = PourRig::default();
    let out = run_pour(&rig, ContactModel::Neutral, newtonian());
    verdict(
        "vsl-002-tilt-pour-mass",
        out.mass_drift < 1e-10 && out.poured_mass > 1.0,
        &format!(
            "mass drift {:.2e} over the tilt schedule; poured mass {:.2} crossed the lip; {} fragments",
            out.mass_drift, out.poured_mass, out.fragments
        ),
    );
}

/// vsl-003: the CONTACT-LINE BRACKET is a first-class output — the
/// same pour under Neutral vs Wetting models yields a REPORTED
/// dribble/pour sensitivity band (the plan's honest handling of the
/// genuinely open problem: the deliverable says how wrong it might
/// be).
#[test]
fn vsl_003_contact_bracket() {
    let rig = PourRig::default();
    let a = run_pour(&rig, ContactModel::Neutral, newtonian());
    let b = run_pour(&rig, ContactModel::Wetting, newtonian());
    let pour_band = (a.poured_mass - b.poured_mass).abs();
    let dribble_band = a.dribble_cells.abs_diff(b.dribble_cells);
    verdict(
        "vsl-003-contact-bracket",
        a.mass_drift < 1e-10 && b.mass_drift < 1e-10,
        &format!(
            "BRACKET (reported, never hidden): poured {:.2} vs {:.2} (band {pour_band:.2}); dribble {} vs {} cells (band {dribble_band}); both ledgers strict",
            a.poured_mass, b.poured_mass, a.dribble_cells, b.dribble_cells
        ),
    );
}

/// vsl-004: the Carreau viscosity band — pours across the fluid
/// family all keep the strict ledger, and the validator score
/// (poured mass) responds to the fluid (the band sensitivity is
/// real, not decorative).
#[test]
fn vsl_004_carreau_band() {
    let rig = PourRig {
        steps: 500,
        ..PourRig::default()
    };
    let fluids = [
        Rheology::Carreau {
            nu0: 0.05,
            nu_inf: 0.005,
            lambda: 5.0,
            n: 0.6,
        },
        Rheology::Carreau {
            nu0: 0.02,
            nu_inf: 0.004,
            lambda: 2.0,
            n: 0.8,
        },
        newtonian(),
    ];
    let mut poured = Vec::new();
    let mut all_strict = true;
    for law in fluids {
        let out = run_pour(&rig, ContactModel::Neutral, law);
        if out.mass_drift >= 1e-10 {
            all_strict = false;
        }
        poured.push(out.poured_mass);
    }
    let spread = poured.iter().fold(f64::NEG_INFINITY, |m, &v| m.max(v))
        - poured.iter().fold(f64::INFINITY, |m, &v| m.min(v));
    verdict(
        "vsl-004-carreau-band",
        all_strict && spread > 0.01,
        &format!(
            "band pours {poured:?}; spread {spread:.3} (the family responds); all ledgers strict"
        ),
    );
}

/// vsl-005: ROBUSTIFICATION — the CVaR-over-band lip beats the
/// nominal-only lip on the off-nominal fluids (the flagship's central
/// claim), and the e-raced candidate screen eliminates dominated lips
/// with its evidence ledgered.
#[test]
fn vsl_005_cvar_and_race() {
    let report = robustify(0.7);
    verdict(
        "vsl-005-cvar-beats-nominal",
        report.robust_offband_growth <= report.nominal_offband_growth
            && report.robust_lip < report.nominal_lip,
        &format!(
            "off-band worst growth: robust {:.5} (lip {:.2}) vs nominal {:.5} (lip {:.2}) — the CVaR design serves the FAMILY",
            report.robust_offband_growth,
            report.robust_lip,
            report.nominal_offband_growth,
            report.nominal_lip
        ),
    );
    // e-raced screen: candidate lips race on the noisy validator
    // proxy (growth + deterministic per-observation jitter), through
    // the crate's PUBLIC wrapper — the wrapper owns the vessel's
    // declared LossSpan convention (losses scaled by SCREEN_SCALE,
    // support = scaled fixture spread + jitter width), so an outside
    // auditor can drive the same convention. It used to live only here
    // in the test file, which is why the fs-flagship-e2e cross-consumer
    // audit could not reach it (bead f85xj.2.31).
    let lips = [0.6f64, 1.0, 1.6, 2.2, 2.8];
    let out = screen_lips(&lips, VESSEL_SCREEN_SEED).expect("fixture screen admits a verdict");
    let expected = out
        .losses
        .iter()
        .enumerate()
        .min_by(|a, b| a.1.total_cmp(b.1))
        .map(|(i, _)| i)
        .expect("nonempty");
    verdict(
        "vsl-005-e-raced-screen",
        out.winner == expected && out.eliminated > 0,
        &format!(
            "race winner lip {} = deterministic argmin lip {} (the U-bottom design); {} dominated candidates eliminated ({} evals vs fixed-N {}) under the vessel's DECLARED span {:.5}",
            lips[out.winner],
            lips[expected],
            out.eliminated,
            out.evaluations_used,
            out.fixed_n_equivalent,
            out.declared_span,
        ),
    );
}

/// vsl-005b: the extracted racing wrapper is EXACTLY the convention the
/// battery used to inline — same declared span, same normalized losses,
/// same verdict — and it refuses structurally rather than racing on
/// evidence it cannot support.
#[test]
fn vsl_005_race_wrapper_owns_the_declared_convention() {
    let lips = [0.6f64, 1.0, 1.6, 2.2, 2.8];
    let base = screening_losses(&lips);
    let base_span = base.iter().copied().fold(f64::NEG_INFINITY, f64::max)
        - base.iter().copied().fold(f64::INFINITY, f64::min);
    let inline_span = 200.0 * (base_span + 1e-4);
    let out = screen_lips(&lips, VESSEL_SCREEN_SEED).expect("fixture screen admits a verdict");
    verdict(
        "vsl-005-declared-span-is-data-derived",
        out.declared_span.to_bits() == inline_span.to_bits()
            && out
                .losses
                .iter()
                .zip(&base)
                .all(|(a, b)| a.to_bits() == b.to_bits()),
        &format!(
            "wrapper declares span {:.6} = 200 x (fixture spread {base_span:.6} + jitter width 1e-4), bit-for-bit with the inlined convention",
            out.declared_span
        ),
    );
    // Driving the same table through the table-level entry point is the
    // same race: this is the surface a cross-consumer audit uses.
    let replay = race_base_losses(&base, VESSEL_SCREEN_SEED).expect("same table, same verdict");
    verdict(
        "vsl-005-table-entry-point-parity",
        replay == out,
        "race_base_losses over the screen's own table reproduces screen_lips exactly",
    );
    verdict(
        "vsl-005-screen-refusal-drill",
        matches!(
            race_base_losses(&[1.0], VESSEL_SCREEN_SEED),
            Err(ScreenError::TooFewCandidates { count: 1 })
        ) && matches!(
            race_base_losses(&[1.0, f64::NAN], VESSEL_SCREEN_SEED),
            Err(ScreenError::NonFiniteLoss { candidate: 1, .. })
        ) && matches!(
            race_base_losses(&[-f64::MAX, f64::MAX], VESSEL_SCREEN_SEED),
            Err(ScreenError::InvalidSpan { .. })
        ),
        "degenerate loss tables return structured refusals instead of a forged declared span",
    );
}

#[test]
fn vsl_005_cvar_rejects_invalid_risk_inputs() {
    let extreme_samples = [-f64::MAX, 0.0, f64::MAX];
    let vessel_report = empirical_cvar(&extreme_samples, 0.25).expect("valid extreme samples");
    let canonical_report =
        fs_robust::empirical_cvar(&extreme_samples, 0.25).expect("canonical extreme samples");
    verdict(
        "vsl-005-canonical-cvar-parity",
        vessel_report == canonical_report
            && cvar(&extreme_samples, 0.25)
                .is_ok_and(|value| value.to_bits() == canonical_report.cvar().to_bits()),
        "vessel report and scalar CVaR surfaces are exact canonical fs-robust re-exports",
    );
    verdict(
        "vsl-005-empty-cvar-drill",
        matches!(empirical_cvar(&[], 0.7), Err(RobustError::EmptySamples)),
        "empty CVaR losses return a structured refusal instead of fake zero risk",
    );
    verdict(
        "vsl-005-bad-beta-drill",
        matches!(
            empirical_cvar(&[1.0, 2.0], 0.0),
            Err(RobustError::BadAlpha { alpha }) if alpha == 0.0
        ),
        "invalid CVaR beta returns a structured refusal before quantile indexing",
    );
    verdict(
        "vsl-005-nan-beta-drill",
        matches!(
            empirical_cvar(&[1.0, 2.0], f64::NAN),
            Err(RobustError::BadAlpha { alpha }) if alpha.is_nan()
        ),
        "non-finite CVaR beta returns a structured refusal before quantile indexing",
    );
    verdict(
        "vsl-005-nonfinite-cvar-drill",
        matches!(
            empirical_cvar(&[1.0, f64::INFINITY], 0.7),
            Err(RobustError::BadSample { value }) if value.is_infinite()
        ),
        "non-finite CVaR losses return a structured refusal before tail aggregation",
    );
}

/// vsl-006: the DELIVERABLE — the pour rendered from the simulation's
/// own mass buffer (zero-copy borrow), bitwise-replayable, with the
/// poured jet visible as structure right of the lip.
#[test]
fn vsl_006_render_same_bytes() {
    let rig = PourRig::default();
    let out = run_pour(&rig, ContactModel::Neutral, newtonian());
    let img1 = render_pour(&out, rig.nx, rig.ny, 24);
    let img2 = render_pour(&out, rig.nx, rig.ny, 24);
    let bitwise = img1
        .iter()
        .zip(&img2)
        .all(|(a, b)| a.to_bits() == b.to_bits());
    let spread = img1.iter().fold(f64::NEG_INFINITY, |m, &v| m.max(v))
        - img1.iter().fold(f64::INFINITY, |m, &v| m.min(v));
    verdict(
        "vsl-006-render-same-bytes",
        bitwise && spread > 0.5,
        &format!(
            "render bound to the sim's own buffer: bitwise replay, transmittance range {spread:.3} (the pour is visible) — the marketing shot IS the physics"
        ),
    );
}
