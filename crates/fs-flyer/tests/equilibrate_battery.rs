//! V-11c battery (beads wf-root-guzez.5.17.1/.2, E4.6d): the full
//! HeldOnRailEquilibrated closure executes with envelope receipts;
//! canonical-anchor refusal at ±1 tick; the settle falsifier; digest
//! sensitivity to seed and branch (and NOTHING intent-shaped — no such
//! field exists); pause/resume-during-prehistory bitwise equivalence;
//! within-branch alternate-start convergence inside declared bands;
//! PrelaunchBranchMismatch both directions; caps; golden.
//! Repro: cargo test -p fs-flyer --test equilibrate_battery

use fs_atmo::ou::{OuMode, STATIONARY_ANCHOR_TICK, StationaryOuPath};
use fs_flyer::aircraft::wright_openloop_v1;
use fs_flyer::canardmech::{CANARD_MECH_V1, MechState};
use fs_flyer::equilibrate::{
    EquilibrationSpec, MAX_PREROLL_TICKS, PrelaunchBranch, SETTLE_TOL_RAD, admit_for_branch,
    equilibrate,
};

fn jlog(case: &str, payload: &str) {
    println!("{{\"suite\":\"fs-flyer-v11c\",\"case\":\"{case}\",{payload}}}");
}

fn ou_modes() -> Vec<OuMode> {
    (0..8)
        .map(|i| {
            OuMode::from_correlation_time(0.4 + 0.05 * f64::from(i), 1.5 + 0.4 * f64::from(i))
                .unwrap()
        })
        .collect()
}

fn spec() -> EquilibrationSpec {
    EquilibrationSpec {
        branch: PrelaunchBranch::HeldOnRailEquilibrated,
        seed: 1903,
        ou_modes: ou_modes(),
        anchor_tick: STATIONARY_ANCHOR_TICK,
        preroll_ticks: 960,
        rho_kg_m3: 1.294,
        trim_start: [13.0, 0.06, 0.1, 45.0],
    }
}

#[test]
fn full_closure_executes_with_envelope_receipts() {
    let d = wright_openloop_v1();
    let t0 = equilibrate(&d, &CANARD_MECH_V1, &spec()).unwrap();
    // Envelope: mech settled onto the trim deflection, trim residuals
    // are the pinned envelope, OU amplitudes finite and stationary-scale.
    assert!(
        (t0.mech.delta_rad - t0.trim.delta_canard_rad).abs() < SETTLE_TOL_RAD,
        "mech must sit on the trim deflection: {} vs {}",
        t0.mech.delta_rad,
        t0.trim.delta_canard_rad
    );
    assert!(t0.mech.rate_rad_s.abs() < 1e-3);
    assert!(t0.settle_ticks > 0 && t0.settle_ticks < 960);
    assert_eq!(t0.ou_amplitudes.len(), 8);
    for (i, a) in t0.ou_amplitudes.iter().enumerate() {
        assert!(a.is_finite() && a.abs() < 10.0, "OU amp {i}: {a}");
    }
    jlog(
        "closure",
        &format!(
            "\"settle_ticks\":{},\"trim_v\":{},\"digest\":\"{}\"",
            t0.settle_ticks,
            t0.trim.v_mps,
            &t0.digest[..16]
        ),
    );
}

#[test]
fn canonical_anchor_refused_at_plus_minus_one() {
    let d = wright_openloop_v1();
    for bad in [STATIONARY_ANCHOR_TICK - 1, STATIONARY_ANCHOR_TICK + 1, 0] {
        let mut s = spec();
        s.anchor_tick = bad;
        assert_eq!(
            equilibrate(&d, &CANARD_MECH_V1, &s).unwrap_err().code,
            "prelaunch-anchor-invalid",
            "anchor {bad} must refuse"
        );
    }
    jlog("anchor", "\"canonical_only\":true");
}

#[test]
fn settle_falsifier_and_preroll_caps() {
    let d = wright_openloop_v1();
    // Too small a budget: the typed not-settled refusal fires.
    let mut short = spec();
    short.preroll_ticks = 3;
    assert_eq!(
        equilibrate(&d, &CANARD_MECH_V1, &short).unwrap_err().code,
        "prelaunch-not-settled"
    );
    // Caps at cap AND cap+1.
    let mut at = spec();
    at.preroll_ticks = MAX_PREROLL_TICKS;
    assert!(equilibrate(&d, &CANARD_MECH_V1, &at).is_ok());
    let mut over = spec();
    over.preroll_ticks = MAX_PREROLL_TICKS + 1;
    assert_eq!(
        equilibrate(&d, &CANARD_MECH_V1, &over).unwrap_err().code,
        "prelaunch-spec-invalid"
    );
    jlog("settle", "\"falsifier_and_caps\":true");
}

#[test]
fn digest_is_seed_and_branch_sensitive_and_repeatable() {
    let d = wright_openloop_v1();
    let a = equilibrate(&d, &CANARD_MECH_V1, &spec()).unwrap();
    let b = equilibrate(&d, &CANARD_MECH_V1, &spec()).unwrap();
    assert_eq!(a.digest, b.digest, "bitwise repeat");
    let mut seeded = spec();
    seeded.seed = 1904;
    let c = equilibrate(&d, &CANARD_MECH_V1, &seeded).unwrap();
    assert_ne!(a.digest, c.digest, "seed sensitivity");
    let mut branched = spec();
    branched.branch = PrelaunchBranch::FreeAirEquilibrated;
    let e = equilibrate(&d, &CANARD_MECH_V1, &branched).unwrap();
    assert_ne!(a.digest, e.digest, "branch sensitivity");
    jlog("digest", "\"seed_branch_sensitive\":true");
}

#[test]
fn prehistory_pause_resume_is_bitwise_equivalent() {
    // The OU suffix property carried through the closure: advancing the
    // prehistory straight vs pausing at −1000 and resuming must land
    // the SAME tick-0 amplitudes, hence the same digest.
    let straight = {
        let mut ou = StationaryOuPath::stationary_at_anchor(1903, ou_modes()).unwrap();
        ou.advance_to(0).unwrap();
        ou.amplitudes().to_vec()
    };
    let paused = {
        let mut ou = StationaryOuPath::stationary_at_anchor(1903, ou_modes()).unwrap();
        ou.advance_to(-1000).unwrap();
        let cp = ou.checkpoint();
        let mut resumed = StationaryOuPath::resume(cp, ou_modes()).unwrap();
        resumed.advance_to(0).unwrap();
        resumed.amplitudes().to_vec()
    };
    for (i, (s, p)) in straight.iter().zip(&paused).enumerate() {
        assert_eq!(
            s.to_bits(),
            p.to_bits(),
            "amplitude {i} must be bitwise equal"
        );
    }
    // And the full closure digest is therefore identical too.
    let d = wright_openloop_v1();
    let a = equilibrate(&d, &CANARD_MECH_V1, &spec()).unwrap();
    let b = equilibrate(&d, &CANARD_MECH_V1, &spec()).unwrap();
    assert_eq!(a.digest, b.digest);
    jlog("pause-resume", "\"bitwise\":true");
}

#[test]
fn within_branch_alternate_starts_converge_inside_bands() {
    // Leaf ii: different mechanism initial states settle to the SAME
    // tick-0 envelope within the declared closure bands (never bitwise —
    // the canonical start defines the digest).
    let d = wright_openloop_v1();
    let canonical = equilibrate(&d, &CANARD_MECH_V1, &spec()).unwrap();
    for start in [
        MechState {
            delta_rad: -0.15,
            rate_rad_s: 0.0,
        },
        MechState {
            delta_rad: 0.20,
            rate_rad_s: 0.5,
        },
    ] {
        let (settled, ticks) = fs_flyer::equilibrate::settle_for_test(
            &CANARD_MECH_V1,
            canonical.trim.delta_canard_rad,
            start,
            960,
        )
        .unwrap();
        // Per-observable bands (declared): deflection 2e-4 rad, rate 2e-3.
        assert!(
            (settled.delta_rad - canonical.mech.delta_rad).abs() < 2.0e-4,
            "alternate start deflection outside the band: {} vs {}",
            settled.delta_rad,
            canonical.mech.delta_rad
        );
        assert!(settled.rate_rad_s.abs() < 2.0e-3);
        jlog(
            "alternate-start",
            &format!(
                "\"start_delta\":{},\"settled\":{},\"ticks\":{ticks}",
                start.delta_rad, settled.delta_rad
            ),
        );
    }
}

#[test]
fn branch_mismatch_refuses_both_directions() {
    let d = wright_openloop_v1();
    let rail = equilibrate(&d, &CANARD_MECH_V1, &spec()).unwrap();
    let mut fspec = spec();
    fspec.branch = PrelaunchBranch::FreeAirEquilibrated;
    let free = equilibrate(&d, &CANARD_MECH_V1, &fspec).unwrap();
    assert!(admit_for_branch(&rail, PrelaunchBranch::HeldOnRailEquilibrated).is_ok());
    assert!(admit_for_branch(&free, PrelaunchBranch::FreeAirEquilibrated).is_ok());
    let e1 = admit_for_branch(&rail, PrelaunchBranch::FreeAirEquilibrated).unwrap_err();
    assert_eq!(e1.code, "PrelaunchBranchMismatch");
    let e2 = admit_for_branch(&free, PrelaunchBranch::HeldOnRailEquilibrated).unwrap_err();
    assert_eq!(e2.code, "PrelaunchBranchMismatch");
    jlog("branch-mismatch", "\"both_directions\":true");
}

#[test]
fn golden_digest() {
    let d = wright_openloop_v1();
    let t0 = equilibrate(&d, &CANARD_MECH_V1, &spec()).unwrap();
    jlog("golden", &format!("\"digest\":\"{}\"", t0.digest));
    assert_eq!(
        t0.digest, "e90ff63a8bf7f1614885e9700d9dfedf97e837111af51cf21facec0b678ddf74",
        "tick-0 golden moved — determinism regression or an intentional \
         closure change requiring the golden-bump protocol"
    );
}
