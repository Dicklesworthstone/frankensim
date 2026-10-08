//! Full-tier battery (plan §15.2 step 3–4 at building scale): multi-story,
//! multi-bay RC moment frames of force-based distributed-plasticity fiber
//! members, gravity-preloaded and modally damped, under Kanai–Tajimi
//! records — modal plausibility, energy conservation of the integrator,
//! metamorphic stiffening, e-stopped building fragility against the
//! fixed-N reference, multi-fidelity MLMC, and bitwise replay.

use fs_frame::history::StoryParams;
use fs_frame::{BuildingFrame, BuildingSpec, building_fragility};
use fs_qty::{Dims, QtyAny};
use fs_scenario::ensemble::{SpectrumModel, StochasticEnsemble};

const TIME: Dims = Dims([0, 0, 1, 0, 0, 0]);
const RATE: Dims = Dims([0, 0, -1, 0, 0, 0]);

fn verdict(name: &str, pass: bool, details: &str) {
    println!("{{\"test\":\"{name}\",\"pass\":{pass},\"details\":\"{details}\"}}");
    assert!(pass, "{name}: {details}");
}

fn kt_ensemble(members: u32, s0: f64, duration: f64, seed: u64) -> StochasticEnsemble {
    StochasticEnsemble {
        name: "kt-building".to_string(),
        seed,
        members,
        duration: QtyAny::new(duration, TIME),
        dt: QtyAny::new(0.02, TIME),
        model: SpectrumModel::KanaiTajimi {
            s0,
            omega_g: QtyAny::new(12.5, RATE),
            zeta_g: 0.6,
        },
    }
}

/// building-001: the gravity-loaded 3-story, 2-bay frame has a fundamental
/// period in the RC moment-frame band with well-separated higher modes; a
/// strong record yields inelastic drifts, absorbs energy, and the Newmark
/// energy ledger closes; stiffer columns reduce the peak drift
/// (metamorphic).
#[test]
fn building_001_modes_energy_and_metamorphic_stiffening() {
    let spec = BuildingSpec::default();
    let base = BuildingFrame::new(spec).expect("default building builds");
    let t = base.periods();
    let ensemble = kt_ensemble(1, 0.02, 8.0, 7);
    let real = ensemble.realize(0).expect("realizes");
    let resp = base
        .clone()
        .run(&real.values, 0.02)
        .expect("history converges");
    let stiff =
        BuildingFrame::new(BuildingSpec { scale: 1.6, ..spec }).expect("stiff building builds");
    let resp_stiff = stiff
        .clone()
        .run(&real.values, 0.02)
        .expect("history converges");
    verdict(
        "building-001-modes",
        t.len() >= 3 && t[0] > 0.3 && t[0] < 2.0 && t[0] > 2.5 * t[1] && t[1] > t[2],
        &format!("periods {:.3} / {:.3} / {:.3} s", t[0], t[1], t[2]),
    );
    verdict(
        "building-001-energy",
        resp.energy_balance_error < 1e-6
            && resp.internal_work > 0.0
            && resp.peak_drift_ratio > 1e-4,
        &format!(
            "energy balance {:.2e}; internal work {:.3e} J; drift profile {:?}",
            resp.energy_balance_error, resp.internal_work, resp.peak_interstory_drift
        ),
    );
    verdict(
        "building-001-metamorphic",
        stiff.periods()[0] < t[0] && resp_stiff.peak_drift_ratio < resp.peak_drift_ratio,
        &format!(
            "column scale 1.6: T1 {:.3} -> {:.3} s, peak drift {:.4} -> {:.4}",
            t[0],
            stiff.periods()[0],
            resp.peak_drift_ratio,
            resp_stiff.peak_drift_ratio
        ),
    );
}

/// building-002: the e-stopped BUILDING fragility's confidence sequence
/// covers the fixed-N exceedance frequency of the full ensemble, and the
/// multi-fidelity MLMC (story model → building correction) reports both
/// levels with a finite estimator variance.
#[test]
fn building_002_e_stopped_fragility_and_multifidelity_mlmc() {
    // A 2-story, 1-bay frame keeps 64 building histories affordable in a
    // debug test; with p ≈ 0.6 the confidence sequence reaches the 0.3
    // margin before the suite is exhausted, so the stop is a real e-stop.
    let members = 64;
    let ensemble = kt_ensemble(members, 0.012, 5.0, 4242);
    let spec = BuildingSpec {
        stories: 2,
        bays: 1,
        ..BuildingSpec::default()
    };
    let limit = 1.0e-2;
    let report = building_fragility(&ensemble, spec, StoryParams::default(), limit, 0.05, 0.3)
        .expect("fragility study converges");
    // Fixed-N reference: the consumed members plus the rest of the suite.
    let base = BuildingFrame::new(spec).expect("builds");
    let mut drifts = report.peak_drifts.clone();
    for m in report.members_used..members {
        let real = ensemble.realize(m).expect("realizes");
        drifts.push(
            base.clone()
                .run(&real.values, 0.02)
                .expect("converges")
                .peak_drift_ratio,
        );
    }
    let exceed = drifts.iter().filter(|d| **d > limit).count();
    let p_ref = exceed as f64 / f64::from(members);
    verdict(
        "building-002-coverage",
        report.stopped_early
            && (report.p_hat - p_ref).abs() <= report.radius
            && exceed > 0
            && exceed < members as usize,
        &format!(
            "p_ref {p_ref:.3} in CS [{:.3} +/- {:.3}] at the e-stop after {}/{members} members ({} exceedances)",
            report.p_hat, report.radius, report.members_used, report.exceedances
        ),
    );
    verdict(
        "building-002-mlmc",
        report.mlmc.levels.len() == 2
            && report.mlmc.estimate.is_finite()
            && report.mlmc.estimate > 0.0
            && report.mlmc.estimator_variance.is_finite(),
        &format!(
            "MLMC E[peak drift] {:.4} (var {:.2e}) over story->building levels",
            report.mlmc.estimate, report.mlmc.estimator_variance
        ),
    );
}

/// building-003: bitwise replay of a building time history, and refusal
/// of degenerate specs.
#[test]
fn building_003_replay_and_refusals() {
    let spec = BuildingSpec {
        stories: 2,
        bays: 1,
        ..BuildingSpec::default()
    };
    let ensemble = kt_ensemble(1, 0.01, 3.0, 11);
    let real = ensemble.realize(0).expect("realizes");
    let run = || {
        BuildingFrame::new(spec)
            .expect("builds")
            .run(&real.values, 0.02)
            .expect("converges")
            .roof_drift
            .iter()
            .map(|x| x.to_bits())
            .collect::<Vec<u64>>()
    };
    let (a, b) = (run(), run());
    verdict(
        "building-003-replay",
        a == b,
        &format!("{} roof samples bit-identical", a.len()),
    );
    let bad = BuildingFrame::new(BuildingSpec {
        stories: 0,
        ..BuildingSpec::default()
    });
    verdict(
        "building-003-refusal",
        bad.is_err(),
        "zero-story spec refused",
    );
}
