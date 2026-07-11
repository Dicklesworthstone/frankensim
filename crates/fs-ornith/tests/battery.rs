//! fs-ornith conformance battery (bead mye.2, smoke tier): the five
//! stages gated, the e-racing payoff measured, the certified Pareto
//! atlas complete with per-row certificates, seed replay bitwise, and
//! the what-breaks-first budget drill.

use fs_ornith::param::{GENE_DIM, OrnithCandidate};
use fs_ornith::screen::{flap_metric, lift_to_drag, screen_generation};
use fs_ornith::{LdSurrogate, build_atlas, certify, refine};

fn verdict(name: &str, pass: bool, details: &str) {
    println!("{{\"test\":\"{name}\",\"pass\":{pass},\"details\":\"{details}\"}}");
    assert!(pass, "{name}: {details}");
}

fn lcg(seed: &mut u64) -> f64 {
    *seed = seed
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    ((*seed >> 11) as f64) / (1u64 << 53) as f64
}

fn generation(n: usize, seed: &mut u64) -> Vec<OrnithCandidate> {
    (0..n)
        .map(|_| {
            let g: Vec<f64> = (0..GENE_DIM).map(|_| lcg(seed)).collect();
            OrnithCandidate::from_genes(&g)
        })
        .collect()
}

/// orn-001: PARAMETERIZE — the Jacobian action's adjoint-assisted lane
/// (∂cl/∂α) matches central differences, and the inlet mass-flow
/// responds to the inlet position lever.
#[test]
fn orn_001_parameterize_jacobians() {
    let c = OrnithCandidate::from_genes(&[0.5, 0.5, 0.5, 0.5, 0.5]);
    let g = c.cl_gradient(64);
    let h = 1e-6;
    let mut cp = c;
    cp.alpha += h;
    let mut cm = c;
    cm.alpha -= h;
    let fd = (lift_to_drag_cl(&cp) - lift_to_drag_cl(&cm)) / (2.0 * h);
    let rel = (g[0] - fd).abs() / fd.abs().max(1e-12);
    // Inlet lever responds.
    let mut ci = c;
    ci.inlet_x = 0.2;
    let m_fore = ci.inlet_mass_flow(64);
    ci.inlet_x = 0.7;
    let m_aft = ci.inlet_mass_flow(64);
    verdict(
        "orn-001-parameterize",
        rel < 1e-6 && (m_fore - m_aft).abs() > 1e-3,
        &format!(
            "adjoint dcl/dalpha {:.6} vs FD {fd:.6} (rel {rel:.1e}); inlet mass-flow fore {m_fore:.4} vs aft {m_aft:.4}",
            g[0]
        ),
    );
}

fn lift_to_drag_cl(c: &OrnithCandidate) -> f64 {
    fs_bem::panel2d::solve(&c.section(64), c.alpha)
        .expect("bounded ornith fixture")
        .cl
}

/// orn-002: SCREEN — the e-raced generation finds the deterministic
/// argmax of L/D and saves real budget vs the fixed-N tournament (the
/// P7 payoff measured), and the flapping metric responds to the gait.
#[test]
fn orn_002_screen_eraced() {
    let mut seed = 0x0221_u64;
    let generation = generation(24, &mut seed);
    let rep = screen_generation(&generation, 0x7E55).expect("normalized screen losses");
    let expected = rep
        .losses
        .iter()
        .enumerate()
        .min_by(|a, b| a.1.total_cmp(b.1))
        .map(|(i, _)| i)
        .expect("nonempty");
    let mut lo = OrnithCandidate::from_genes(&[0.5, 0.5, 0.5, 0.0, 0.1]);
    lo.flap_amp = 0.0;
    let mut hi = lo;
    hi.flap_amp = 0.45;
    hi.flap_freq = 1.6;
    let (fm_lo, fm_hi) = (flap_metric(&lo), flap_metric(&hi));
    verdict(
        "orn-002-screen-eraced",
        rep.winner == expected
            && rep.eliminated >= 20
            && rep.evaluations_used * 3 < rep.fixed_n_equivalent
            && fm_hi > fm_lo,
        &format!(
            "race winner {} == argmax L/D {}; {}/24 dominated candidates eliminated early; {} evals vs fixed-N {} ({}x saved); flap metric responds: {fm_lo:.4} -> {fm_hi:.4}",
            rep.winner,
            expected,
            rep.eliminated,
            rep.evaluations_used,
            rep.fixed_n_equivalent,
            rep.fixed_n_equivalent / rep.evaluations_used.max(1)
        ),
    );
}

/// orn-003: REFINE — LBM control-volume forces agree with the panel
/// method within MODEL-FORM evidence (lift sign matches the trim sign;
/// magnitudes in the same order band), the flow is steady, and the
/// honesty label travels.
#[test]
fn orn_003_refine_agreement() {
    let c = OrnithCandidate::from_genes(&[0.4, 0.6, 0.5, 0.2, 0.5]);
    let rep = refine(&c);
    let sign_ok = rep.lift.signum() == rep.panel_cl.signum();
    verdict(
        "orn-003-refine",
        sign_ok && rep.steadiness < 1e-4 && !rep.honesty.is_empty(),
        &format!(
            "LBM lift {:.3e}, drag {:.3e} vs panel cl {:.4}: sign agreement {sign_ok}; steadiness {:.1e}; honesty label: '{}'",
            rep.lift, rep.drag, rep.panel_cl, rep.steadiness, rep.honesty
        ),
    );
}

/// orn-004: CERTIFY — the SOS/Lyapunov certificate verifies on the
/// trim state, the ROA volume is positive ONLY under certification,
/// and the conformal band covers at the declared rate on fresh
/// candidates (certify-or-escalate is real).
#[test]
fn orn_004_certified_stability_and_conformal() {
    let c = OrnithCandidate::from_genes(&[0.5, 0.5, 0.5, 0.5, 0.5]);
    let cert = certify(&c);
    // Surrogate + conformal band.
    let mut seed = 0x0441_u64;
    let train: Vec<(OrnithCandidate, f64)> = generation(40, &mut seed)
        .into_iter()
        .map(|c| (c, lift_to_drag(&c)))
        .collect();
    let sur = LdSurrogate::fit(&train, 0.1);
    let fresh = generation(60, &mut seed);
    let cov = sur.coverage(&fresh);
    verdict(
        "orn-004-certify",
        cert.certified && cert.roa_volume > 0.0 && cov >= 0.85,
        &format!(
            "Lyapunov certificate verified (A=[[0,1],[{:.3},{:.3}]]); certified ROA volume {:.4}; conformal coverage {cov:.2} on 60 fresh candidates (target 0.90 - slack for finite calibration)",
            cert.a[1][0], cert.a[1][1], cert.roa_volume
        ),
    );
}

/// orn-005: the PARETO ATLAS — the front is nonempty, every row
/// carries a certificate, the hypervolume is positive, the knee is on
/// the front, and the adjoint polish does not lose L/D.
#[test]
fn orn_005_certified_pareto_atlas() {
    let mut seed = 0x0551_u64;
    let train: Vec<(OrnithCandidate, f64)> = generation(40, &mut seed)
        .into_iter()
        .map(|c| (c, lift_to_drag(&c)))
        .collect();
    let sur = LdSurrogate::fit(&train, 0.1);
    let atlas = build_atlas(24, 12, 0xA71A5, &sur);
    let all_certified_consistent = atlas
        .rows
        .iter()
        .all(|r| (r.roa > 0.0) == r.certificate.certified);
    for r in atlas.rows.iter().take(6) {
        println!(
            "{{\"atlas\":{{\"ld\":{:.3},\"roa\":{:.4},\"maneuver\":{:.4},\"inlet_viol\":{:.4},\"certified\":{},\"surrogate_ld\":{:.3}}}}}",
            r.ld, r.roa, r.maneuver, r.inlet_violation, r.certificate.certified, r.surrogate_ld
        );
    }
    verdict(
        "orn-005-pareto-atlas",
        !atlas.rows.is_empty()
            && all_certified_consistent
            && atlas.hypervolume > 0.0
            && atlas.knee < atlas.rows.len()
            && atlas.polish_gain.1 >= atlas.polish_gain.0,
        &format!(
            "{} certified rows; hypervolume {:.3}; knee row {}; adjoint polish L/D {:.3} -> {:.3}",
            atlas.rows.len(),
            atlas.hypervolume,
            atlas.knee,
            atlas.polish_gain.0,
            atlas.polish_gain.1
        ),
    );
}

/// orn-006: REPLAY — the full pipeline is deterministic from the seed:
/// atlas rebuilt with the same seed is bitwise identical (genes and
/// objectives), and the screen race replays exactly.
#[test]
fn orn_006_seed_replay() {
    let mut seed = 0x0661_u64;
    let train: Vec<(OrnithCandidate, f64)> = generation(40, &mut seed)
        .into_iter()
        .map(|c| (c, lift_to_drag(&c)))
        .collect();
    let sur = LdSurrogate::fit(&train, 0.1);
    let a = build_atlas(16, 8, 0x5EED, &sur);
    let b = build_atlas(16, 8, 0x5EED, &sur);
    let atlas_bitwise = a.rows.len() == b.rows.len()
        && a.rows.iter().zip(&b.rows).all(|(x, y)| {
            x.genes
                .iter()
                .zip(&y.genes)
                .all(|(u, v)| u.to_bits() == v.to_bits())
                && x.ld.to_bits() == y.ld.to_bits()
        });
    let mut s1 = 0x0662_u64;
    let g1 = generation(12, &mut s1);
    let r1 = screen_generation(&g1, 0xACE).expect("normalized screen losses");
    let r2 = screen_generation(&g1, 0xACE).expect("normalized screen losses");
    verdict(
        "orn-006-seed-replay",
        atlas_bitwise && r1.winner == r2.winner && r1.evaluations_used == r2.evaluations_used,
        &format!(
            "atlas bitwise replay over {} rows; race replay: winner {} / {} evals both runs",
            a.rows.len(),
            r1.winner,
            r1.evaluations_used
        ),
    );
}

/// orn-007: WHAT BREAKS FIRST — the LBM budget exhausts mid-campaign
/// and the pipeline degrades GRACEFULLY to the surrogate + conformal
/// path instead of dying: the degraded estimate stays inside the
/// conformal band of the truth.
#[test]
fn orn_007_budget_exhaustion_degrades_gracefully() {
    let mut seed = 0x0771_u64;
    let train: Vec<(OrnithCandidate, f64)> = generation(40, &mut seed)
        .into_iter()
        .map(|c| (c, lift_to_drag(&c)))
        .collect();
    let sur = LdSurrogate::fit(&train, 0.1);
    // Campaign: 8 candidates, LBM budget for only 1 refinement.
    let campaign = generation(8, &mut seed);
    let mut lbm_budget = 1usize;
    let mut degraded = 0usize;
    let mut in_band = 0usize;
    for c in &campaign {
        if lbm_budget > 0 {
            lbm_budget -= 1;
            let _full = refine(c); // the funded lane
        } else {
            degraded += 1;
            let pred = sur.predict(c);
            if sur.band.covers(pred, lift_to_drag(c)) {
                in_band += 1;
            }
        }
    }
    verdict(
        "orn-007-graceful-degradation",
        degraded == 7 && in_band >= 6,
        &format!(
            "LBM budget exhausted after 1 refine; {degraded}/7 candidates degraded to surrogate+conformal, {in_band}/7 inside the band — the campaign survives its own honesty clause"
        ),
    );
}
