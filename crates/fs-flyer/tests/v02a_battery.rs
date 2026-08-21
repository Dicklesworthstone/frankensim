//! V-02a GATE battery (bead wf-root-guzez.5.13.3, E4.6a-iii): the
//! open-loop model meets the A4 claims at the dossier-permitted level
//! (unstable mode present, time-to-double inside the order band and
//! REPORTED, derivative signs), the deficient uncoupled-canard baseline
//! scores WORSE (anti-vacuity, executed), the eigensolver is verified on
//! analytic fixtures, Iyy insensitivity of the instability sign is
//! executed across the declared ±25%, determinism, golden.
//! Repro: cargo test -p fs-flyer --test v02a_battery

use fs_flyer::aircraft::wright_openloop_v1;
use fs_flyer::longitudinal::{IYY_KG_M2, T2_BAND_S, eig4, linearize, v02a_gate};

fn jlog(case: &str, payload: &str) {
    println!("{{\"suite\":\"fs-flyer-v02a\",\"case\":\"{case}\",{payload}}}");
}

const RHO: f64 = 1.294;
const START: [f64; 4] = [13.0, 0.06, 0.1, 45.0];

#[test]
fn v02a_gate_passes_on_the_integrated_model() {
    let d = wright_openloop_v1();
    let t = d.trim(RHO, START).unwrap();
    let rep = linearize(&d, &t, RHO).unwrap();
    let poles: Vec<String> = rep
        .poles
        .iter()
        .map(|p| format!("[{},{}]", p.re, p.im))
        .collect();
    jlog(
        "gate-receipt",
        &format!(
            "\"poles\":[{}],\"t2_s\":{},\"m_alpha\":{},\"m_q\":{},\"derivs\":{:?}",
            poles.join(","),
            rep.time_to_double_s,
            rep.m_alpha_nm_per_rad,
            rep.m_q_nm_s_per_rad,
            rep.derivatives
        ),
    );
    let v = v02a_gate(&rep);
    assert!(
        v.unstable_mode_present,
        "A4: an unstable longitudinal mode must exist"
    );
    assert!(
        v.m_alpha_sign_ok,
        "A4: M_alpha must be positive (statically unstable): {}",
        rep.m_alpha_nm_per_rad
    );
    assert!(
        v.m_q_sign_ok,
        "A4: M_q must be negative (pitch damping): {}",
        rep.m_q_nm_s_per_rad
    );
    assert!(
        v.t2_in_band,
        "A4: time-to-double {} s outside the order band {:?}",
        rep.time_to_double_s, T2_BAND_S
    );
    assert!(v.passes());
}

#[test]
fn deficient_uncoupled_canard_scores_worse() {
    // ANTI-VACUITY (executed, plan DONE-WHEN): a model with canard-wing
    // interference effectively removed (canard displaced 8 m below the
    // wing system — same arms, same areas, mutual induction dead) must
    // score WORSE against the A4 claims than the coupled model. If it
    // scored the same, the gate could not distinguish coupling fidelity
    // and would be vacuous.
    let d = wright_openloop_v1();
    let t = d.trim(RHO, START).unwrap();
    let coupled = v02a_gate(&linearize(&d, &t, RHO).unwrap());
    let mut deficient = wright_openloop_v1();
    deficient.canard_z_m = [8.35, 9.05];
    let td = deficient.trim(RHO, START).unwrap();
    let unc = v02a_gate(&linearize(&deficient, &td, RHO).unwrap());
    jlog(
        "anti-vacuity",
        &format!(
            "\"coupled_score\":{},\"deficient_score\":{}",
            coupled.score, unc.score
        ),
    );
    assert!(
        unc.score > coupled.score + 0.05,
        "the deficient baseline must score WORSE: {} vs {}",
        unc.score,
        coupled.score
    );
}

#[test]
fn instability_sign_is_iyy_insensitive_across_declared_uncertainty() {
    // The anchors file declares Iyy Estimated +/-25% and claims the pole
    // SIGN is insensitive to it — execute that claim by rescaling the
    // q-dot row (Iyy enters the A matrix only there).
    let d = wright_openloop_v1();
    let t = d.trim(RHO, START).unwrap();
    let rep = linearize(&d, &t, RHO).unwrap();
    for scale in [0.75, 1.25] {
        let mut a = rep.a;
        for c in 0..4 {
            a[2][c] *= 1.0 / scale; // Iyy' = scale*Iyy => row / scale
        }
        // Rebuild poles through the public path: a scaled-Iyy report is
        // just the same matrix with the q-row rescaled; reuse eig via a
        // fresh linearize call is not possible (Iyy is a const), so we
        // check the DOMINANT sign through the trace/determinant route:
        // an unstable real mode survives iff p(0) = det(-A)... simplest
        // executable check: power iteration on the matrix exponential
        // surrogate (I + h A)^n from a fixed start grows.
        let h = 0.001;
        let mut x = [1.0, 0.5, 0.2, 0.1];
        let mut norm0 = 0.0;
        for step in 0..40_000 {
            let mut y = [0.0f64; 4];
            for (i, yi) in y.iter_mut().enumerate() {
                *yi = x[i] + h * (0..4).map(|j| a[i][j] * x[j]).sum::<f64>();
            }
            x = y;
            if step == 0 {
                norm0 = x.iter().map(|v| v * v).sum::<f64>().sqrt();
            }
        }
        let norm1 = x.iter().map(|v| v * v).sum::<f64>().sqrt();
        assert!(
            norm1 > norm0 * 10.0,
            "instability must survive Iyy scale {scale}: growth {norm1:.3e} vs {norm0:.3e}"
        );
    }
    jlog("iyy-insensitivity", &format!("\"iyy_nominal\":{IYY_KG_M2}"));
}

#[test]
fn eigensolver_matches_analytic_fixtures() {
    // Diagonal fixture: exact real eigenvalues in canonical order.
    let a = [
        [2.0, 0.0, 0.0, 0.0],
        [0.0, -1.0, 0.0, 0.0],
        [0.0, 0.0, 0.5, 0.0],
        [0.0, 0.0, 0.0, -3.0],
    ];
    let p = eig4(&a);
    let expect = [2.0, 0.5, -1.0, -3.0];
    for (pole, e) in p.iter().zip(expect) {
        assert!(
            (pole.re - e).abs() < 1e-8 && pole.im.abs() < 1e-8,
            "diag fixture: {pole:?} vs {e}"
        );
    }
    // Complex-pair fixture: rotation block (sigma +/- i*omega) plus two
    // reals — the classic unstable-oscillation shape.
    let (sg, om) = (0.3, 2.0);
    let b = [
        [sg, om, 0.0, 0.0],
        [-om, sg, 0.0, 0.0],
        [0.0, 0.0, -0.7, 0.0],
        [0.0, 0.0, 0.0, -1.4],
    ];
    let pb = eig4(&b);
    assert!((pb[0].re - sg).abs() < 1e-8 && (pb[0].im.abs() - om).abs() < 1e-8);
    assert!((pb[1].re - sg).abs() < 1e-8);
    assert!((pb[2].re + 0.7).abs() < 1e-8);
    assert!((pb[3].re + 1.4).abs() < 1e-8);
    jlog("eig-fixtures", "\"analytic_match\":true");
}

#[test]
fn poles_are_deterministic() {
    let d = wright_openloop_v1();
    let t = d.trim(RHO, START).unwrap();
    let a = linearize(&d, &t, RHO).unwrap();
    let b = linearize(&d, &t, RHO).unwrap();
    for (pa, pb) in a.poles.iter().zip(b.poles.iter()) {
        assert_eq!(pa.re.to_bits(), pb.re.to_bits());
        assert_eq!(pa.im.to_bits(), pb.im.to_bits());
    }
    jlog("determinism", "\"bitwise\":true");
}

#[test]
fn v02a_golden_digest() {
    let d = wright_openloop_v1();
    let t = d.trim(RHO, START).unwrap();
    let rep = linearize(&d, &t, RHO).unwrap();
    let mut payload = Vec::new();
    for p in &rep.poles {
        payload.extend_from_slice(&p.re.to_bits().to_le_bytes());
        payload.extend_from_slice(&p.im.to_bits().to_le_bytes());
    }
    payload.extend_from_slice(&rep.m_alpha_nm_per_rad.to_bits().to_le_bytes());
    payload.extend_from_slice(&rep.m_q_nm_s_per_rad.to_bits().to_le_bytes());
    let digest =
        fs_blake3::hash_domain("org.frankensim.fs-flyer.v02a-golden.v1", &payload).to_hex();
    jlog(
        "golden",
        &format!("\"digest\":\"{digest}\",\"t2\":{}", rep.time_to_double_s),
    );
    assert_eq!(
        digest, "9d97aadd85219571537fbdf6c01ac76e65048104fc94eebdf5a32e50d373f270",
        "V-02a pole golden moved — determinism regression or an \
         intentional model change requiring the golden-bump protocol"
    );
}
