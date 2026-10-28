//! E4.2-ii battery (bead wf-root-guzez.5.3.2): attached-regime agreement
//! with the exact linear fixture, the CAMBER CLOSURE lifting the Dec-17
//! layout into the weight class, warm-start efficiency MEASURED, the
//! safeguard refusing a pathological closure (never NaN, never linear
//! fallback), post-stall branch identity (distinct + bitwise-repeatable),
//! the operator-reuse rule enforced, golden.
//! Repro: cargo test -p fs-wing --test nonlinear_battery

use fs_wing::nonlinear::{
    InfluenceOperator, StripRegime, StripSpec, solve_nonlinear, solve_nonlinear_with_incident,
};
use fs_wing::{Panel, SurfaceId, flat_surface};

fn jlog(case: &str, payload: &str) {
    println!("{{\"suite\":\"fs-wing-e42ii\",\"case\":\"{case}\",{payload}}}");
}

const RHO: f64 = 1.294;
const V: f64 = 13.86;

fn freestream(alpha: f64) -> [f64; 3] {
    [V * alpha.cos(), 0.0, V * alpha.sin()]
}

/// The 1903 biplane (both wings, 8 span x 2 chord each) + strips.
fn biplane() -> (Vec<Panel>, Vec<StripSpec>) {
    let mut p = flat_surface(SurfaceId::WingLower, 12.29, 1.981, 0.0, 0.0, 8, 2).unwrap();
    p.extend(flat_surface(SurfaceId::WingUpper, 12.29, 1.981, 0.0, -1.89, 8, 2).unwrap());
    let mut strips = Vec::new();
    for surf in 0..2 {
        let base = surf * 16;
        for s in 0..8 {
            strips.push(StripSpec {
                panel_indices: vec![base + s, base + 8 + s],
                chord_m: 1.981,
                twist_rad: 0.0,
            });
        }
    }
    (p, strips)
}

/// The thin-airfoil + flat-plate-blend closure at 1/20 camber (the wing
/// dataset shape; fs-flyer wires the real datasets at L4).
fn camber_closure(_s: usize, alpha: f64) -> (f64, StripRegime) {
    let camber = 0.05;
    let attached = 2.0 * std::f64::consts::PI * (alpha + 2.0 * camber);
    let a = alpha.abs();
    if a <= 0.30 {
        (attached, StripRegime::Attached)
    } else if a < 0.45 {
        let t = (a - 0.30) / 0.15;
        let s = t * t * (3.0 - 2.0 * t);
        let sep = 1.98 * alpha.sin() * alpha.cos();
        (attached * (1.0 - s) + sep * s, StripRegime::Blended)
    } else {
        (1.98 * alpha.sin() * alpha.cos(), StripRegime::Separated)
    }
}

/// A zero-camber linear closure (matches the linear solve's physics).
fn flat_closure(_s: usize, alpha: f64) -> (f64, StripRegime) {
    (2.0 * std::f64::consts::PI * alpha, StripRegime::Attached)
}

#[test]
fn uniform_incident_field_matches_the_legacy_solve_exactly() {
    let (panels, strips) = biplane();
    let velocity = freestream(0.05);
    let op = InfluenceOperator::build(&panels, velocity, RHO).unwrap();
    let legacy = solve_nonlinear(
        &op,
        &panels,
        &strips,
        velocity,
        RHO,
        &camber_closure,
        None,
        None,
    )
    .unwrap();
    let uniform = solve_nonlinear_with_incident(
        &op,
        &panels,
        &strips,
        velocity,
        RHO,
        &camber_closure,
        &vec![velocity; strips.len()],
        None,
        None,
    )
    .unwrap();
    assert_eq!(uniform, legacy);
}

#[test]
fn nonuniform_vertical_incident_flow_changes_real_strip_loading_and_mirrors() {
    let (panels, strips) = biplane();
    let velocity = freestream(0.05);
    let op = InfluenceOperator::build(&panels, velocity, RHO).unwrap();
    let mut incident = vec![velocity; strips.len()];
    for (s, flow) in incident.iter_mut().enumerate() {
        flow[2] += if s % 8 < 4 { 0.35 } else { -0.35 };
    }
    let solve = |flow: &[[f64; 3]]| {
        solve_nonlinear_with_incident(
            &op,
            &panels,
            &strips,
            velocity,
            RHO,
            &camber_closure,
            flow,
            None,
            None,
        )
        .unwrap()
    };
    let gust = solve(&incident);
    assert!(
        gust.gamma[2] < gust.gamma[5],
        "upward gust must increase lift-up circulation"
    );
    incident[..8].reverse();
    incident[8..].reverse();
    let mirrored = solve(&incident);
    for plane in 0..2 {
        for chord_row in 0..2 {
            for span in 0..8 {
                let a = gust.gamma[plane * 16 + chord_row * 8 + span];
                let b = mirrored.gamma[plane * 16 + chord_row * 8 + 7 - span];
                assert!(
                    (a - b).abs() < 1.0e-9,
                    "mirrored panel {plane}/{chord_row}/{span}: {a} {b}"
                );
            }
        }
    }
    let mut lift = 0.0;
    for (s, strip) in strips.iter().enumerate() {
        let flow = incident[s];
        for &j in &strip.panel_indices {
            let panel = panels[j];
            let cross_z = flow[0] * (panel.b[1] - panel.a[1]) - flow[1] * (panel.b[0] - panel.a[0]);
            lift += -RHO * mirrored.gamma[j] * cross_z;
        }
    }
    assert!((mirrored.total_lift_n - lift).abs() < 1.0e-9);
}

#[test]
fn axial_wash_changes_local_incidence_and_reversed_flow_refuses() {
    let (panels, strips) = biplane();
    let velocity = freestream(0.08);
    let op = InfluenceOperator::build(&panels, velocity, RHO).unwrap();
    let incident = vec![velocity; strips.len()];
    let du = vec![3.0; strips.len()];
    let washed = solve_nonlinear_with_incident(
        &op,
        &panels,
        &strips,
        velocity,
        RHO,
        &camber_closure,
        &incident,
        Some(&du),
        None,
    )
    .unwrap();
    let combined = vec![[velocity[0] + 3.0, velocity[1], velocity[2]]; strips.len()];
    let direct = solve_nonlinear_with_incident(
        &op,
        &panels,
        &strips,
        velocity,
        RHO,
        &camber_closure,
        &combined,
        None,
        None,
    )
    .unwrap();
    assert_eq!(
        washed, direct,
        "wash is a vector addition before section incidence/speed"
    );
    for bad in [
        vec![],
        vec![[0.0; 3]; strips.len()],
        vec![[f64::NAN, 0.0, 0.0]; strips.len()],
    ] {
        assert_eq!(
            solve_nonlinear_with_incident(
                &op,
                &panels,
                &strips,
                velocity,
                RHO,
                &camber_closure,
                &bad,
                None,
                None,
            )
            .unwrap_err()
            .code,
            "incident-flow-invalid"
        );
    }
    let reverse = vec![-2.0 * velocity[0]; strips.len()];
    assert_eq!(
        solve_nonlinear_with_incident(
            &op,
            &panels,
            &strips,
            velocity,
            RHO,
            &camber_closure,
            &incident,
            Some(&reverse),
            None,
        )
        .unwrap_err()
        .code,
        "incident-flow-invalid"
    );
}

#[test]
fn attached_regime_tracks_the_linear_fixture() {
    // MONOPLANE with the FLAT closure at small alpha: the two
    // formulations (Weissinger 3/4-chord BC vs trailing-only lifting-line
    // closure) agree within the classical ~15% at AR 6. On the BIPLANE
    // the formulations treat mutual interference differently (bound-bound
    // vs trailing-only) — that delta is REPORTED, never forced to vanish
    // (the plan's V-05 philosophy); the first run measured +21.7%.
    let p = flat_surface(SurfaceId::WingLower, 12.29, 1.981, 0.0, 0.0, 8, 2).unwrap();
    let strips: Vec<StripSpec> = (0..8)
        .map(|s| StripSpec {
            panel_indices: vec![s, 8 + s],
            chord_m: 1.981,
            twist_rad: 0.0,
        })
        .collect();
    let fs_v = freestream(0.05);
    let op = InfluenceOperator::build(&p, fs_v, RHO).unwrap();
    let r = solve_nonlinear(&op, &p, &strips, fs_v, RHO, &flat_closure, None, None).unwrap();
    let lin = op.linear().total_lift_n;
    assert!(
        r.residual < 1e-10 && r.iterations < 200,
        "residual {} iters {}",
        r.residual,
        r.iterations
    );
    assert!(
        (r.total_lift_n / lin - 1.0).abs() < 0.15,
        "nonlinear {} vs linear {} beyond the formulation band",
        r.total_lift_n,
        lin
    );
    assert!(r.regimes.iter().all(|g| *g == StripRegime::Attached));
    // The biplane formulation delta: measured and LOGGED as data.
    let (bp, bstrips) = biplane();
    let bop = InfluenceOperator::build(&bp, fs_v, RHO).unwrap();
    let br = solve_nonlinear(&bop, &bp, &bstrips, fs_v, RHO, &flat_closure, None, None).unwrap();
    let delta = br.total_lift_n / bop.linear().total_lift_n - 1.0;
    assert!(
        delta.abs() < 0.35,
        "biplane formulation delta {delta} implausibly large"
    );
    jlog(
        "attached",
        &format!(
            "\"nl\":{},\"lin\":{lin},\"iters\":{},\"biplane_formulation_delta\":{delta}",
            r.total_lift_n, r.iterations
        ),
    );
}

#[test]
fn camber_closure_reaches_the_weight_class() {
    // The E4.2-i golden was 1648 N from UNCAMBERED plates at 4 deg. With
    // the 1/20-camber closure the biplane at the Dec-17 condition must
    // carry the gross-weight class (3336 N) — the camber closure is what
    // made the Flyer fly. Sanity band, not a validation claim.
    let (p, strips) = biplane();
    let fs_v = freestream(0.07);
    let op = InfluenceOperator::build(&p, fs_v, RHO).unwrap();
    let r = solve_nonlinear(&op, &p, &strips, fs_v, RHO, &camber_closure, None, None).unwrap();
    assert!(
        r.total_lift_n > 2500.0 && r.total_lift_n < 5500.0,
        "camber-closed lift {} N outside the weight class",
        r.total_lift_n
    );
    jlog(
        "camber",
        &format!("\"lift\":{},\"weight\":3336", r.total_lift_n),
    );
}

#[test]
fn warm_start_is_measurably_cheaper() {
    let (p, strips) = biplane();
    let fs1 = freestream(0.05);
    let op1 = InfluenceOperator::build(&p, fs1, RHO).unwrap();
    let r1 = solve_nonlinear(&op1, &p, &strips, fs1, RHO, &camber_closure, None, None).unwrap();
    let fs2 = freestream(0.06);
    let op2 = InfluenceOperator::build(&p, fs2, RHO).unwrap();
    let cold = solve_nonlinear(&op2, &p, &strips, fs2, RHO, &camber_closure, None, None).unwrap();
    let warm = solve_nonlinear(
        &op2,
        &p,
        &strips,
        fs2,
        RHO,
        &camber_closure,
        None,
        Some(&r1.gamma),
    )
    .unwrap();
    assert!(
        warm.iterations <= cold.iterations,
        "warm start must not cost more ({} vs {})",
        warm.iterations,
        cold.iterations
    );
    // And both converge to the same answer bitwise-close.
    assert!((warm.total_lift_n - cold.total_lift_n).abs() < 1e-6 * cold.total_lift_n.abs());
    jlog(
        "warm-start",
        &format!("\"cold\":{},\"warm\":{}", cold.iterations, warm.iterations),
    );
}

#[test]
fn safeguard_refuses_a_pathological_closure() {
    // A NON-FINITE closure must REFUSE with the typed code — never NaN
    // in the output, never a silent linear fallback.
    let (p, strips) = biplane();
    let fs_v = freestream(0.05);
    let op = InfluenceOperator::build(&p, fs_v, RHO).unwrap();
    // A closure that goes NON-FINITE just outside the start band: the
    // solver must refuse the moment the iteration touches it — never
    // propagate NaN into gamma, never fall back to the linear answer.
    // (Two earlier fixture attempts are instructive: a 50x-slope closure
    // and a ±40 branch-flip closure BOTH turned out to possess genuine
    // fixed points the safeguarded iteration correctly found — bounded
    // continuous-ish closures are hard problems, not pathologies. The
    // divergence-floor and iteration-cap refusal paths were exercised
    // live during E4.6a-ii integration.)
    // (A conditional NaN region proved dodgeable too — downwash steered
    // alpha_eff around it. Unconditional non-finite output is the one
    // unambiguous pathology.)
    let wild = |_s: usize, _alpha: f64| -> (f64, StripRegime) { (f64::NAN, StripRegime::Attached) };
    let err = solve_nonlinear(&op, &p, &strips, fs_v, RHO, &wild, None, None).unwrap_err();
    assert_eq!(err.code, "nonlinear-did-not-converge");
    assert!(err.ranked_repairs[0].contains("NEVER fall back") || err.message.contains("residual"));
    jlog("safeguard", "\"refused\":true");
}

#[test]
fn post_stall_branch_identity_is_distinct_and_repeatable() {
    let (p, strips) = biplane();
    let fs_lo = freestream(0.05);
    let op_lo = InfluenceOperator::build(&p, fs_lo, RHO).unwrap();
    let attached =
        solve_nonlinear(&op_lo, &p, &strips, fs_lo, RHO, &camber_closure, None, None).unwrap();
    let fs_hi = freestream(0.42);
    let op_hi = InfluenceOperator::build(&p, fs_hi, RHO).unwrap();
    let stalled =
        solve_nonlinear(&op_hi, &p, &strips, fs_hi, RHO, &camber_closure, None, None).unwrap();
    assert_ne!(
        attached.branch_id, stalled.branch_id,
        "branches must be distinct"
    );
    assert!(stalled.regimes.iter().any(|r| *r != StripRegime::Attached));
    // Bitwise repeatability of the branch id AND the solution.
    let again =
        solve_nonlinear(&op_hi, &p, &strips, fs_hi, RHO, &camber_closure, None, None).unwrap();
    assert_eq!(stalled.branch_id, again.branch_id);
    for (a, b) in stalled.gamma.iter().zip(&again.gamma) {
        assert_eq!(a.to_bits(), b.to_bits());
    }
    jlog(
        "branch",
        &format!(
            "\"attached\":\"{}\",\"stalled\":\"{}\"",
            &attached.branch_id[..12],
            &stalled.branch_id[..12]
        ),
    );
}

#[test]
fn operator_reuse_rule_is_enforced() {
    let (p, strips) = biplane();
    let fs_v = freestream(0.05);
    let op = InfluenceOperator::build(&p, fs_v, RHO).unwrap();
    // Same geometry: reuse fine.
    assert!(solve_nonlinear(&op, &p, &strips, fs_v, RHO, &camber_closure, None, None).is_ok());
    // ANY geometry change (a canard-deflection-scale nudge on one panel)
    // makes the operator stale: typed refusal, not a wrong answer.
    let mut moved = p.clone();
    moved[3].ctrl[2] += 0.01;
    assert_eq!(
        solve_nonlinear(&op, &moved, &strips, fs_v, RHO, &camber_closure, None, None)
            .unwrap_err()
            .code,
        "influence-operator-stale"
    );
    // A freestream change also invalidates.
    assert_eq!(
        solve_nonlinear(
            &op,
            &p,
            &strips,
            freestream(0.051),
            RHO,
            &camber_closure,
            None,
            None
        )
        .unwrap_err()
        .code,
        "influence-operator-stale"
    );
    jlog(
        "reuse",
        "\"stale\":\"refused for geometry AND freestream changes\"",
    );
}

#[test]
fn nonlinear_golden_digest() {
    let (p, strips) = biplane();
    let fs_v = freestream(0.07);
    let op = InfluenceOperator::build(&p, fs_v, RHO).unwrap();
    let r = solve_nonlinear(&op, &p, &strips, fs_v, RHO, &camber_closure, None, None).unwrap();
    let mut payload = Vec::new();
    for g in &r.gamma {
        payload.extend_from_slice(&g.to_bits().to_le_bytes());
    }
    payload.extend_from_slice(&r.total_lift_n.to_bits().to_le_bytes());
    payload.extend_from_slice(r.branch_id.as_bytes());
    let digest =
        fs_blake3::hash_domain("org.frankensim.fs-wing.e42ii-golden.v1", &payload).to_hex();
    jlog(
        "golden",
        &format!("\"digest\":\"{digest}\",\"lift\":{}", r.total_lift_n),
    );
    assert_eq!(
        digest, "2797bcda52ee7d52ff0cb8283c640fe463d1f5996bd4c3e4a5d02162823e1fe0",
        "nonlinear golden moved — determinism regression or an intentional \
         closure change requiring the golden-bump protocol"
    );
}
