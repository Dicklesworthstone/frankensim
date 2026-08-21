//! E4.7-iii battery (bead wf-root-guzez.5.18.3): multipole convergence
//! vs the EXACT retained segments (order check, not a vibe); the
//! registered WakeCoreEvolutionMode enters identity AND physically
//! feeds the pruning bounds; hybrid = near-exact + far-multipole
//! agrees with the pre-aggregation exact field; audited pruning with
//! the EXECUTED adversarial bound-deflation falsifier
//! (pruning-certificate-failed TERMINATES, no mutation); V-10 receipt
//! vs the fs-wakeref dense reference (shape-class, honest scale note,
//! Tier A/B KPI delta REPORTED); caps at cap AND cap+1; determinism
//! golden.
//! Repro: cargo test -p fs-vpm --test farfield3d_battery

use fs_vpm::coarsen3d::coarsen_oldest;
use fs_vpm::farfield3d::{
    FarField, MAX_NU_EFF, WakeCoreEvolutionMode, emit_v10_receipt, hybrid_velocity,
};
use fs_vpm::filament3d::{FilamentWake, WakeRateCertificate};
use fs_wakeref::{Fixture, RefereeCase, run_case, wright_geometry_v1};

fn jlog(case: &str, payload: &str) {
    println!("{{\"suite\":\"fs-vpm-farfield3d\",\"case\":\"{case}\",{payload}}}");
}

const DT: f64 = 1.0 / 120.0;
const V_MPS: f64 = 13.0;

fn line(n: usize) -> Vec<[f64; 3]> {
    (0..=n)
        .map(|i| [0.0, i as f64 - n as f64 / 2.0, 0.0])
        .collect()
}

/// Step-started elliptic-loading wake convecting downstream at V.
fn step_wake(rows: usize) -> FilamentWake {
    let cert = WakeRateCertificate {
        shed_hz: 120.0,
        n_stations: 8,
        max_rows: rows + 4,
    };
    let mut wake = FilamentWake::new(cert, line(8)).unwrap();
    let g: Vec<f64> = (0..8)
        .map(|s| {
            let y = (s as f64 - 3.5) / 4.0;
            5.0 * (1.0 - y * y).max(0.0)
        })
        .collect();
    for _ in 0..rows {
        wake.shed(&g, [-V_MPS * DT, 0.0, -0.002]).unwrap();
    }
    wake
}

fn mag(v: [f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

#[test]
fn multipole_converges_to_exact_with_distance() {
    let mut wake = step_wake(64);
    let far = FarField::aggregate(&mut wake, 32, 8, WakeCoreEvolutionMode::Frozen, DT).unwrap();
    let cell = &far.cells[0];
    // Relative expansion error at standoff d and 2d: the truncation is
    // quadrupole-order, so doubling the standoff must shrink the
    // relative error by well over 4x (cubic-class decay).
    let dir = [0.3, 0.6, 0.74];
    let base = 3.0 * cell.effective_radius_m();
    let mut errs = Vec::new();
    for k in 0..2 {
        let d = base * (1 << k) as f64;
        let p = [
            cell.centroid[0] + d * dir[0],
            cell.centroid[1] + d * dir[1],
            cell.centroid[2] + d * dir[2],
        ];
        let e = cell.exact(p);
        let a = cell.eval(p);
        let rel = mag([a[0] - e[0], a[1] - e[1], a[2] - e[2]]) / mag(e).max(1e-300);
        errs.push(rel);
    }
    assert!(errs[0] < 5e-2, "usable at 3R: {}", errs[0]);
    assert!(
        errs[0] / errs[1] > 4.0,
        "superquadratic decay: {} -> {}",
        errs[0],
        errs[1]
    );
    jlog(
        "convergence",
        &format!("\"rel_3r\":{},\"rel_6r\":{}", errs[0], errs[1]),
    );
}

#[test]
fn core_mode_enters_identity_and_feeds_bounds() {
    let frozen = {
        let mut w = step_wake(64);
        FarField::aggregate(&mut w, 32, 8, WakeCoreEvolutionMode::Frozen, DT).unwrap()
    };
    let spreading = {
        let mut w = step_wake(64);
        FarField::aggregate(
            &mut w,
            32,
            8,
            WakeCoreEvolutionMode::CoreSpreading { nu_eff_m2ps: 0.05 },
            DT,
        )
        .unwrap()
    };
    // Identity: the mode is hashed, and the spread physically changes
    // the cell state (core2), so the digests MUST differ.
    assert_ne!(frozen.digest(), spreading.digest(), "mode enters identity");
    // Physical coupling: spreading inflates the effective radius,
    // which TIGHTENS prunability (larger bound at the same probe).
    let p = [30.0, 0.0, 0.0];
    let bf = frozen.cells[0].contribution_bound(p);
    let bs = spreading.cells[0].contribution_bound(p);
    assert!(
        bs > bf,
        "spread radius must enlarge the bound: {bs} vs {bf}"
    );
    // Admission at cap AND beyond-cap.
    assert!(
        WakeCoreEvolutionMode::CoreSpreading {
            nu_eff_m2ps: MAX_NU_EFF
        }
        .admit()
        .is_ok()
    );
    for bad in [MAX_NU_EFF * 1.0000000001, 0.0, -1.0, f64::NAN] {
        assert_eq!(
            WakeCoreEvolutionMode::CoreSpreading { nu_eff_m2ps: bad }
                .admit()
                .unwrap_err()
                .code,
            "core-mode-invalid",
            "nu {bad}"
        );
    }
    jlog(
        "core-mode",
        &format!("\"bound_frozen\":{bf},\"bound_spreading\":{bs}"),
    );
}

#[test]
fn aggregation_caps_and_hybrid_matches_exact() {
    // Caps: exact multiple admits; ±1 row refuses; all rows refuses.
    let mut w = step_wake(64);
    assert_eq!(
        FarField::aggregate(&mut w, 33, 8, WakeCoreEvolutionMode::Frozen, DT)
            .unwrap_err()
            .code,
        "farfield-invalid"
    );
    assert_eq!(
        FarField::aggregate(&mut w, 31, 8, WakeCoreEvolutionMode::Frozen, DT)
            .unwrap_err()
            .code,
        "farfield-invalid"
    );
    assert_eq!(
        FarField::aggregate(&mut w, 64, 8, WakeCoreEvolutionMode::Frozen, DT)
            .unwrap_err()
            .code,
        "farfield-invalid"
    );
    assert_eq!(
        FarField::aggregate(&mut w, 32, 8, WakeCoreEvolutionMode::Frozen, 0.0)
            .unwrap_err()
            .code,
        "farfield-invalid"
    );
    // Hybrid vs the pre-aggregation EXACT field: the near probe sits
    // by the lifting line, the aggregated rows are >14 m downstream,
    // so the hybrid must reproduce the exact field to expansion error.
    let full = step_wake(64);
    let mut near = full.clone();
    let far = FarField::aggregate(&mut near, 32, 8, WakeCoreEvolutionMode::Frozen, DT).unwrap();
    assert_eq!(near.rows.len(), 32);
    assert_eq!(far.cells.len(), 4);
    let probes = [[0.3, 0.0, 0.05], [0.1, 2.0, 0.2], [-1.0, -1.5, 0.4]];
    let mut worst = 0.0f64;
    for p in probes {
        let e = full.induced_velocity(p);
        let h = hybrid_velocity(&near, &far, p);
        let rel = mag([h[0] - e[0], h[1] - e[1], h[2] - e[2]]) / mag(e).max(1e-300);
        worst = worst.max(rel);
    }
    // Measured 2026-08-21: worst 1.39e-2 — the cells span the full
    // 8 m wake width, so the nearest cell sits at r/d ~ 0.6 from the
    // worst probe and the quadrupole truncation is the honest limit.
    // Bound holds ~3.5x headroom; re-measure if the fixture moves.
    assert!(worst < 5e-2, "hybrid vs exact rel {worst}");
    jlog("hybrid", &format!("\"worst_rel\":{worst}"));
}

#[test]
fn audited_prune_lawful_and_certificate_falsifier() {
    let probes = [[0.3, 0.0, 0.05], [0.0, 2.0, 0.1]];
    // Lawful: with the wake's oldest cells ~20+ m downstream, a small
    // tolerance prunes the distant cells and the audit certifies each.
    let mut near = step_wake(240);
    let far_full =
        FarField::aggregate(&mut near, 120, 10, WakeCoreEvolutionMode::Frozen, DT).unwrap();
    let before: Vec<[f64; 3]> = probes
        .iter()
        .map(|p| hybrid_velocity(&near, &far_full, *p))
        .collect();
    let mut far = far_full.clone();
    // Measured 2026-08-21: farthest-cell bounds ~0.039, nearest ~0.6;
    // exact spot-checks of pruned cells ~1e-3. 0.05 prunes the far
    // tail while the near cells stay (bound >= tol).
    let tol = 5e-2;
    let receipt = far.prune_audited(&probes, tol).unwrap();
    assert!(
        !receipt.pruned.is_empty(),
        "the scenario must actually prune"
    );
    assert!(receipt.kept + receipt.pruned.len() == far_full.cells.len());
    // Every audit row carries its evidence, and every exact spot-check
    // sits below tol AND below its own bound (per-item oracles).
    for r in &receipt.pruned {
        assert!(
            r.exact_worst < tol,
            "cell {}: {}",
            r.cell_index,
            r.exact_worst
        );
        assert!(
            r.exact_worst <= r.bound_worst,
            "bound must actually bound: {} > {}",
            r.exact_worst,
            r.bound_worst
        );
    }
    // Whole-field effect stays bounded by the certified budget.
    for (i, p) in probes.iter().enumerate() {
        let after = hybrid_velocity(&near, &far, *p);
        let d = mag([
            after[0] - before[i][0],
            after[1] - before[i][1],
            after[2] - before[i][2],
        ]);
        assert!(
            d <= tol * receipt.pruned.len() as f64,
            "probe {i}: prune moved the field {d}"
        );
    }
    // FALSIFIER (executed): adversarially deflated bounds forge
    // candidacy for cells the audit then catches red-handed.
    let mut forged = far_full.clone();
    let n_before = forged.cells.len();
    let err = forged.prune_audited_scaled(&probes, tol, 1e-9).unwrap_err();
    assert_eq!(err.code, "pruning-certificate-failed");
    assert_eq!(
        forged.cells.len(),
        n_before,
        "certificate failure must TERMINATE with no mutation"
    );
    // Refusals: no probes; tol at zero, negative, NaN.
    let mut f2 = far_full.clone();
    assert_eq!(
        f2.prune_audited(&[], tol).unwrap_err().code,
        "prune-invalid"
    );
    for bad in [0.0, -1e-3, f64::NAN] {
        assert_eq!(
            f2.prune_audited(&probes, bad).unwrap_err().code,
            "prune-invalid",
            "tol {bad}"
        );
    }
    jlog(
        "prune",
        &format!(
            "\"pruned\":{},\"kept\":{},\"falsifier_code\":\"{}\"",
            receipt.pruned.len(),
            receipt.kept,
            err.code
        ),
    );
}

#[test]
fn v10_receipt_vs_wakeref_dense_reference() {
    // fs-vpm side: wake-induced downwash buildup at a fixed probe
    // while the step wake grows, tick by tick.
    let cert = WakeRateCertificate {
        shed_hz: 120.0,
        n_stations: 8,
        max_rows: 244,
    };
    let mut wake = FilamentWake::new(cert, line(8)).unwrap();
    let g: Vec<f64> = (0..8)
        .map(|s| {
            let y = (s as f64 - 3.5) / 4.0;
            5.0 * (1.0 - y * y).max(0.0)
        })
        .collect();
    let probe = [0.3, 0.0, 0.05];
    let mut buildup = Vec::with_capacity(240);
    for _ in 0..240 {
        wake.shed(&g, [-V_MPS * DT, 0.0, -0.002]).unwrap();
        buildup.push(mag(wake.induced_velocity(probe)));
    }
    // Referee side: dense UVLM step case, free air (same dt, 2 s).
    let series = run_case(
        &wright_geometry_v1(),
        &RefereeCase {
            fixture: Fixture::Step,
            ground_z_m: None,
            v_mps: V_MPS,
            alpha0_rad: 0.05,
            rho_kg_m3: 1.294,
            convection: 1.0,
            dt_s: DT,
            n_steps: 240,
        },
    )
    .unwrap();
    // Overlap window from tick 3 (the referee's declared
    // apparent-mass spike decays by then), each normalized by its own
    // terminal value — shape class ONLY.
    let hyb: Vec<f64> = buildup[3..].iter().map(|v| v / buildup[239]).collect();
    // The Wagner-class buildup lives on the CANARD (the stepped
    // surface); the wing series instead starts ABOVE steady (no wake
    // downwash yet) and decays — measured 1.023 at tick 3.
    let referee: Vec<f64> = series.canard_lift_n[3..]
        .iter()
        .map(|v| v / series.canard_lift_n[239])
        .collect();
    // Tier A/B KPI: probe speed at the end of the window with the
    // full exact wake (A) vs coarsen + multipole + audited prune (B).
    let tier_a = mag(wake.induced_velocity(probe));
    let mut b_wake = wake.clone();
    coarsen_oldest(&mut b_wake, 60).unwrap();
    let mut far =
        FarField::aggregate(&mut b_wake, 120, 10, WakeCoreEvolutionMode::Frozen, DT).unwrap();
    far.prune_audited(&[probe], 2e-4).unwrap();
    let tier_b = mag(hybrid_velocity(&b_wake, &far, probe));
    let receipt = emit_v10_receipt(
        WakeCoreEvolutionMode::Frozen,
        &hyb,
        &referee,
        tier_a,
        tier_b,
    )
    .unwrap();
    jlog(
        "v10",
        &format!(
            "\"shape_rms\":{},\"kpi_a\":{},\"kpi_b\":{},\"kpi_delta\":{},\"digest\":\"{}\"",
            receipt.shape_rms,
            receipt.tier_a_kpi_mps,
            receipt.tier_b_kpi_mps,
            receipt.kpi_delta_mps,
            receipt.receipt_digest
        ),
    );
    // Both are Wagner-class buildups ending at 1 by construction.
    assert!(receipt.terminal_delta < 1e-12);
    // Non-vacuity: the referee's Wagner deficiency is visible and the
    // hybrid buildup genuinely rises to its terminal. (The filament
    // probe sits in the immediate near field, so ITS buildup is fast
    // by physics — the shape window compares the slow tail.)
    // Measured 2026-08-21: canard buildup starts ~0.93 (the shallow
    // Wagner deficiency fs-wakeref's CONTRACT declares).
    assert!(referee[0] < 0.99, "referee Wagner start: {}", referee[0]);
    assert!(hyb[0] < 0.999, "hybrid buildup must rise: {}", hyb[0]);
    // Measured 2026-08-21 shape RMS 0.139 (near-field downwash vs
    // dense-UVLM lift buildup, declared shape-class); ~3.5x headroom.
    assert!(
        receipt.shape_rms < 0.5,
        "shape-class agreement: {}",
        receipt.shape_rms
    );
    // KPI delta is REPORTED, bounded loosely, never forced to vanish.
    assert!(receipt.kpi_delta_mps.abs() / tier_a.max(1e-12) < 0.2);
    // Refusal: mismatched overlap.
    assert_eq!(
        emit_v10_receipt(WakeCoreEvolutionMode::Frozen, &hyb, &referee[1..], 1.0, 1.0)
            .unwrap_err()
            .code,
        "v10-invalid"
    );
}

#[test]
fn determinism_golden() {
    let run = || {
        let mut w = step_wake(240);
        let mut far =
            FarField::aggregate(&mut w, 120, 10, WakeCoreEvolutionMode::Frozen, DT).unwrap();
        far.prune_audited(&[[0.3, 0.0, 0.05]], 2e-4).unwrap();
        far.digest()
    };
    let a = run();
    assert_eq!(a, run(), "bit-identical twice");
    jlog("golden", &format!("\"digest\":\"{a}\""));
    assert_eq!(
        a, "d55b9ef65b16d42e5bc85f838e8ba79a6d84932dba4211cd2beb4f7ace8f863f",
        "far-field golden moved — determinism regression or an \
         intentional multipole/prune change requiring the golden-bump protocol"
    );
}
