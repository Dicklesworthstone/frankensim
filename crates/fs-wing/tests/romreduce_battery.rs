//! E4.3b2-ii battery (bead wf-root-guzez.5.8.2): the shared-basis
//! block rational-Krylov ladder EXECUTED over a 9-anchor spread of the
//! frozen grid — smallest passing order wins; MIMO transfer (incl. the
//! hinge channel) checked on HELD-OUT frequencies; reduced poles
//! strictly stable (causality for a rational LTI); the FORBIDDEN
//! per-point-basis interpolation FALSIFIER (shared basis interpolates
//! the convection axis accurately, per-point bases do not); the
//! ladder-exhausted refusal executed via an adversarial tolerance;
//! determinism golden.
//! Repro: cargo test -p fs-wing --test romreduce_battery --release

use fs_wing::images::CertifiedGround;
use fs_wing::prescribedwake::{WakeOperatingPoint, frozen_grid_v1};
use fs_wing::rom::{A1Lti, assemble_a1_lti, wright_a1_layout_v1};
use fs_wing::romreduce::{
    HELD_OUT_W, LADDER, SHIFTS_0, TRANSFER_TOL, project, reduce_shared, reduce_shared_with_tol,
    small_eigenvalues, transfer_of,
};

fn jlog(case: &str, payload: &str) {
    println!("{{\"suite\":\"fs-wing-romreduce\",\"case\":\"{case}\",{payload}}}");
}

fn ground() -> CertifiedGround {
    CertifiedGround {
        z_m: 3.0,
        certificate_slope: 0.000606,
        certificate_rms_m: 0.801,
    }
}

const V: f64 = 13.0;
const ROWS: usize = 120;

/// Nine anchors spanning the grid (corners + center classes).
fn anchor_points() -> Vec<WakeOperatingPoint> {
    let grid = frozen_grid_v1();
    [0usize, 17, 40, 71, 100, 130, 143, 200, 287]
        .iter()
        .map(|&i| grid.points[i])
        .collect()
}

fn anchors() -> Vec<A1Lti> {
    let layout = wright_a1_layout_v1();
    anchor_points()
        .iter()
        .map(|p| assemble_a1_lti(&layout, p, &ground(), V, ROWS).unwrap())
        .collect()
}

#[test]
fn ladder_finds_a_small_passing_order_and_is_deterministic() {
    let sys = anchors();
    let refs: Vec<&A1Lti> = sys.iter().collect();
    let red = reduce_shared(&refs).unwrap();
    assert!(LADDER.contains(&red.order));
    assert!(
        red.order < 24,
        "reduced below the FOM order 24: {}",
        red.order
    );
    // Every rung BEFORE the winner failed (smallest-passing law).
    for rung in &red.ladder[..red.ladder.len() - 1] {
        assert!(!rung.passed, "{rung:?} passed before the winner");
    }
    let last = red.ladder.last().unwrap();
    assert!(last.passed && last.worst_rel_err < TRANSFER_TOL);
    let again = reduce_shared(&refs).unwrap();
    assert_eq!(red.digest, again.digest, "bit-identical twice");
    jlog(
        "ladder",
        &format!(
            "\"order\":{},\"worst_rel_err\":{},\"digest\":\"{}\"",
            red.order, last.worst_rel_err, red.digest
        ),
    );
    assert_eq!(
        red.digest, "7bf1914228e5486e0d7c84fbe187eb2d90b8dc951eb0c00944e47c08e9960a27",
        "reduction golden moved — determinism regression or an \
         intentional reduction change requiring the golden-bump protocol"
    );
}

#[test]
fn heldout_transfer_matches_on_every_channel_incl_hinge() {
    // Per-channel oracle at held-out frequencies for a NON-anchor grid
    // point projected through the shared basis (the scheduling claim's
    // precondition; full scheduling is 5.8.3).
    let sys = anchors();
    let refs: Vec<&A1Lti> = sys.iter().collect();
    let red = reduce_shared(&refs).unwrap();
    let layout = wright_a1_layout_v1();
    let grid = frozen_grid_v1();
    let holdout = assemble_a1_lti(&layout, &grid.points[60], &ground(), V, ROWS).unwrap();
    let r = project(&holdout, &red.basis, red.order);
    let mut worst = (0usize, 0.0f64);
    for &w in &HELD_OUT_W {
        let gf = transfer_of(
            &holdout.a,
            &holdout.b,
            &holdout.c,
            &holdout.d,
            holdout.order,
            w,
        )
        .unwrap();
        let gr = transfer_of(&r.a, &r.b, &r.c, &r.d, r.order, w).unwrap();
        for ch in 0..6 {
            let mag = gf[ch].0.hypot(gf[ch].1);
            if mag > 1e-9 {
                let err = (gf[ch].0 - gr[ch].0).hypot(gf[ch].1 - gr[ch].1) / mag;
                if err > worst.1 {
                    worst = (ch, err);
                }
            }
        }
    }
    // The hinge channels are ch 4/5 (output 2) — included above; the
    // non-anchor point is allowed a modestly looser band than the
    // anchor tolerance (declared: scheduling interpolation is 5.8.3).
    assert!(
        worst.1 < 5.0 * TRANSFER_TOL,
        "held-out point worst: {worst:?}"
    );
    jlog(
        "heldout",
        &format!("\"worst_ch\":{},\"worst_rel\":{}", worst.0, worst.1),
    );
}

#[test]
fn reduced_poles_are_strictly_stable() {
    let sys = anchors();
    let refs: Vec<&A1Lti> = sys.iter().collect();
    let red = reduce_shared(&refs).unwrap();
    for lti in &sys {
        let r = project(lti, &red.basis, red.order);
        let eigs = small_eigenvalues(&r.a, r.order).unwrap();
        assert_eq!(eigs.len(), r.order);
        for (re, im) in &eigs {
            assert!(
                *re < 0.0,
                "unstable reduced pole {re}+{im}i at {:?}",
                lti.point
            );
        }
    }
    jlog(
        "poles",
        &format!("\"systems\":{},\"all_stable\":true", sys.len()),
    );
}

#[test]
fn per_point_basis_interpolation_is_the_executed_falsifier() {
    // The FORBIDDEN scheme: reduce two convection twins independently
    // and interpolate their reduced matrices — the coordinates are
    // incompatible and the midpoint transfer is garbage; the SHARED
    // basis midpoint stays accurate against the directly-assembled
    // mid-convection truth.
    let layout = wright_a1_layout_v1();
    let grid = frozen_grid_v1();
    let base = grid.points[1]; // convection axis is the innermost pair
    let mut p085 = base;
    p085.convection = 0.85;
    let mut p100 = base;
    p100.convection = 1.0;
    let mut pmid = base;
    pmid.convection = 0.925;
    let s085 = assemble_a1_lti(&layout, &p085, &ground(), V, ROWS).unwrap();
    let s100 = assemble_a1_lti(&layout, &p100, &ground(), V, ROWS).unwrap();
    let truth = assemble_a1_lti(&layout, &pmid, &ground(), V, ROWS).unwrap();
    let w = 3.7;
    let gt = transfer_of(&truth.a, &truth.b, &truth.c, &truth.d, truth.order, w).unwrap();
    let err_of = |ga: &[(f64, f64); 6]| -> f64 {
        let mut worst = 0.0f64;
        for ch in 0..6 {
            let mag = gt[ch].0.hypot(gt[ch].1);
            if mag > 1e-9 {
                worst = worst.max((gt[ch].0 - ga[ch].0).hypot(gt[ch].1 - ga[ch].1) / mag);
            }
        }
        worst
    };
    // Shared scheme.
    let refs = [&s085, &s100];
    let shared = reduce_shared(&refs).unwrap();
    let ra = project(&s085, &shared.basis, shared.order);
    let rb = project(&s100, &shared.basis, shared.order);
    let mix = |x: &[f64], y: &[f64]| -> Vec<f64> {
        x.iter().zip(y.iter()).map(|(a, b)| 0.5 * (a + b)).collect()
    };
    let am = mix(&ra.a, &rb.a);
    let bm = mix(&ra.b, &rb.b);
    let cm = mix(&ra.c, &rb.c);
    let g_shared = transfer_of(&am, &bm, &cm, &ra.d, shared.order, w).unwrap();
    let shared_err = err_of(&g_shared);
    // Forbidden per-point scheme. A per-point basis is defined only up
    // to an orthogonal rotation WITHIN its own span — nothing anchors
    // the coordinates across operating points. Inject that legal
    // freedom (the swap-injection pattern): rotate one basis in-span,
    // verify the rotated basis is EQUALLY VALID (identical transfer at
    // its own point — the liveness check), then show the interpolation
    // is garbage. The shared basis has no such freedom by construction.
    let red_a = reduce_shared(&[&s085]).unwrap();
    let red_b = reduce_shared(&[&s100]).unwrap();
    let order = red_a.order.min(red_b.order).min(shared.order);
    let n = s100.order;
    // Deterministic in-span rotation of basis B: Givens cascade over
    // adjacent column pairs at a fixed angle.
    let mut basis_b_rot = red_b.basis.clone();
    let (cth, sth) = (0.36_f64.cos(), 0.36_f64.sin());
    for pair in 0..order.saturating_sub(1) {
        for i in 0..n {
            let u = basis_b_rot[pair * n + i];
            let v = basis_b_rot[(pair + 1) * n + i];
            basis_b_rot[pair * n + i] = cth * u + sth * v;
            basis_b_rot[(pair + 1) * n + i] = -sth * u + cth * v;
        }
    }
    let pa = project(&s085, &red_a.basis, order);
    let pb_plain = project(&s100, &red_b.basis, order);
    let pb_rot = project(&s100, &basis_b_rot, order);
    // Liveness: the rotated basis is an equally valid per-point ROM.
    let g_plain =
        transfer_of(&pb_plain.a, &pb_plain.b, &pb_plain.c, &pb_plain.d, order, w).unwrap();
    let g_rot = transfer_of(&pb_rot.a, &pb_rot.b, &pb_rot.c, &pb_rot.d, order, w).unwrap();
    for ch in 0..6 {
        let mag = g_plain[ch].0.hypot(g_plain[ch].1).max(1e-9);
        assert!(
            (g_plain[ch].0 - g_rot[ch].0).hypot(g_plain[ch].1 - g_rot[ch].1) / mag < 1e-8,
            "rotation must not change the per-point transfer (equally valid basis)"
        );
    }
    let am2 = mix(&pa.a, &pb_rot.a);
    let bm2 = mix(&pa.b, &pb_rot.b);
    let cm2 = mix(&pa.c, &pb_rot.c);
    let g_pp = transfer_of(&am2, &bm2, &cm2, &pa.d, order, w).unwrap();
    let pp_err = err_of(&g_pp);
    assert!(
        shared_err < 0.05,
        "shared-basis midpoint accurate: {shared_err}"
    );
    assert!(
        pp_err > 10.0 * shared_err.max(1e-4),
        "per-point interpolation must be garbage under the legal in-span rotation: {pp_err} vs {shared_err}"
    );
    jlog(
        "falsifier",
        &format!("\"shared_err\":{shared_err},\"per_point_err\":{pp_err}"),
    );
}

#[test]
fn ladder_exhaustion_refuses_and_forbids_a1() {
    // Adversarial tolerance: no order can pass — the typed refusal is
    // the plan's A1-forbidden outcome, executed.
    let sys = anchors();
    let refs: Vec<&A1Lti> = sys.iter().collect();
    let err = reduce_shared_with_tol(&refs, 1.0e-15).unwrap_err();
    assert_eq!(err.code, "rom-ladder-exhausted");
    assert!(err.ranked_repairs.iter().any(|r| r.contains("FORBIDDEN")));
    // Empty anchors refuse too.
    assert_eq!(reduce_shared(&[]).unwrap_err().code, "rom-anchors-empty");
    jlog("refusal", "\"ladder_exhausted_typed\":true");
}

#[test]
fn shifts_never_collide_with_held_out_frequencies() {
    for s in SHIFTS_0 {
        for w in HELD_OUT_W {
            assert!((s - w).abs() > 1e-9, "shift {s} collides with held-out {w}");
        }
    }
}
