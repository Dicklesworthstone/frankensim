//! V-03 core battery (bead wf-root-guzez.5.12.1, E4.5-i): the E1.6
//! HOLDOUT static anchor (the 1903 reconstruction at 350 rpm vs the
//! Wrights'/LFST bench numbers, factor-2 trend band — the reconstruction
//! is Estimated by the dossier's own rule), CT monotone in J, the eta
//! curve's rise-and-fall with a peak in the holdout band's neighborhood,
//! Prandtl limits, convergence receipts, spin-up toward the 350-rpm
//! operating class, caps, golden.
//! Repro: cargo test -p fs-airscrew --test bemt_battery

use fs_airscrew::{
    BladeStation, MAX_STATIONS, Rotor, bemt_solve, engine_torque_at_prop_nm, rotor_spinup_step,
};

fn jlog(case: &str, payload: &str) {
    println!("{{\"suite\":\"fs-airscrew-v03\",\"case\":\"{case}\",{payload}}}");
}

const RHO: f64 = 1.294; // the E1.8 cold-day density

/// The 1903 reconstruction: the 1911 NTRS planform TREND scaled to the
/// sparse 1903 numbers (E1.6's declared reconstruction path #2 —
/// Estimated ceiling; tip width 8 in, straight tip, R = 1.2954 m).
fn rotor_1903() -> Rotor {
    let deg = std::f64::consts::PI / 180.0;
    Rotor {
        radius_m: 1.2954,
        n_blades: 2,
        camber_ratio: 0.04,
        stations: vec![
            BladeStation {
                r_over_r: 0.30,
                chord_m: 0.13,
                beta_rad: 40.0 * deg,
            },
            BladeStation {
                r_over_r: 0.45,
                chord_m: 0.17,
                beta_rad: 30.0 * deg,
            },
            BladeStation {
                r_over_r: 0.60,
                chord_m: 0.20,
                beta_rad: 23.0 * deg,
            },
            BladeStation {
                r_over_r: 0.75,
                chord_m: 0.21,
                beta_rad: 18.5 * deg,
            },
            BladeStation {
                r_over_r: 0.88,
                chord_m: 0.20,
                beta_rad: 15.5 * deg,
            },
            BladeStation {
                r_over_r: 0.96,
                chord_m: 0.16,
                beta_rad: 14.0 * deg,
            },
        ],
    }
}

const OMEGA_350: f64 = 350.0 / 60.0 * 2.0 * std::f64::consts::PI;

#[test]
fn static_anchor_holdout_trend() {
    // HOLDOUT (E1.6 partition): the Wrights measured 132-136 lb/pair at
    // 350 rpm (67±1 lb per prop = 298 N); the LFST repro read 285 N. The
    // reconstruction carries an Estimated ceiling: the permitted claim is
    // the factor-2 TREND band around the anchor.
    let s = bemt_solve(&rotor_1903(), RHO, 0.0, OMEGA_350).unwrap();
    let anchor_n = 285.0;
    assert!(s.thrust_n > 0.0, "static thrust must be positive");
    assert!(
        s.thrust_n / anchor_n > 0.5 && s.thrust_n / anchor_n < 2.0,
        "static thrust {} N outside the factor-2 band around {anchor_n} N",
        s.thrust_n
    );
    assert!((s.j - 0.0).abs() < 1e-12, "J = 0 exactly at the bench");
    for r in &s.stations {
        assert!(
            r.iterations < 120 && r.w_mps > 0.0,
            "receipt at r/R {}",
            r.r_over_r
        );
    }
    jlog(
        "static-anchor",
        &format!("\"thrust_n\":{},\"anchor_n\":{anchor_n}", s.thrust_n),
    );
}

#[test]
fn ct_monotone_and_eta_rises_then_falls() {
    let rotor = rotor_1903();
    let n = OMEGA_350 / (2.0 * std::f64::consts::PI);
    let d = 2.0 * rotor.radius_m;
    let mut prev_ct = f64::INFINITY;
    let mut etas = Vec::new();
    for k in 1..=9 {
        let j = 0.1 * f64::from(k);
        let v = j * n * d;
        let s = bemt_solve(&rotor, RHO, v, OMEGA_350).unwrap();
        assert!(s.ct < prev_ct, "CT must fall with J (at J = {j})");
        prev_ct = s.ct;
        let cp = s.cq * 2.0 * std::f64::consts::PI;
        etas.push((j, if cp > 0.0 { j * s.ct / cp } else { 0.0 }));
    }
    let peak = etas
        .iter()
        .cloned()
        .fold((0.0, 0.0), |a, b| if b.1 > a.1 { b } else { a });
    assert!(
        peak.1 > 0.5 && peak.1 < 0.95,
        "eta peak {} implausible",
        peak.1
    );
    assert!(
        peak.0 >= 0.5 && peak.0 <= 0.9,
        "peak J {} outside the Wright class",
        peak.0
    );
    let last = etas.last().unwrap();
    assert!(last.1 < peak.1, "eta must fall past the peak");
    jlog(
        "eta",
        &format!("\"peak_j\":{},\"peak_eta\":{}", peak.0, peak.1),
    );
}

#[test]
fn prandtl_limits_and_receipts() {
    let s = bemt_solve(&rotor_1903(), RHO, 10.7, OMEGA_350).unwrap();
    let mut prev_f = 0.0;
    for (i, r) in s.stations.iter().enumerate() {
        assert!(
            r.prandtl_f > 0.0 && r.prandtl_f <= 1.0,
            "F in (0,1] at {}",
            r.r_over_r
        );
        if i > 0 && r.r_over_r > 0.6 {
            assert!(r.prandtl_f <= prev_f + 0.05, "F must fall toward the tip");
        }
        prev_f = r.prandtl_f;
    }
    let tip = s.stations.last().unwrap();
    assert!(
        tip.prandtl_f < 0.85,
        "tip loss must bite at r/R 0.96 (F = {})",
        tip.prandtl_f
    );
    jlog("prandtl", &format!("\"tip_f\":{}", tip.prandtl_f));
}

#[test]
fn spinup_reaches_the_operating_class() {
    // From rest torque exceeds prop torque -> Omega climbs; equilibrium
    // where the declared engine curve crosses the prop's CQ. The 1903
    // operating class is 330-380 prop rpm.
    let rotor = rotor_1903();
    let mut omega = 5.0f64;
    let dt = 1.0 / 120.0;
    for _ in 0..2400 {
        let q_prop = bemt_solve(&rotor, RHO, 0.0, omega).unwrap().torque_nm;
        let q_eng = engine_torque_at_prop_nm(omega);
        omega = rotor_spinup_step(2.2, omega, q_eng, q_prop, dt).unwrap();
    }
    let rpm = omega * 60.0 / (2.0 * std::f64::consts::PI);
    assert!(
        rpm > 250.0 && rpm < 500.0,
        "equilibrium {rpm} rpm outside the operating class"
    );
    jlog("spinup", &format!("\"equilibrium_rpm\":{rpm}"));
}

#[test]
fn refusals_at_cap_and_cap_plus_one() {
    let mut r = rotor_1903();
    let st = r.stations[2];
    while r.stations.len() < MAX_STATIONS {
        let mut s2 = st;
        s2.r_over_r = r.stations.last().unwrap().r_over_r + 0.0004;
        r.stations.push(s2);
    }
    assert!(
        bemt_solve(&r, RHO, 5.0, OMEGA_350).is_ok(),
        "AT the cap solves"
    );
    let mut over = r.clone();
    let mut s2 = st;
    s2.r_over_r = over.stations.last().unwrap().r_over_r + 0.0004;
    over.stations.push(s2);
    assert_eq!(
        bemt_solve(&over, RHO, 5.0, OMEGA_350).unwrap_err().code,
        "rotor-invalid"
    );
    // Operating-point gates.
    assert_eq!(
        bemt_solve(&rotor_1903(), RHO, -1.0, OMEGA_350)
            .unwrap_err()
            .code,
        "operating-point-invalid"
    );
    assert_eq!(
        bemt_solve(&rotor_1903(), RHO, 5.0, 0.0).unwrap_err().code,
        "operating-point-invalid"
    );
    // Descending stations refuse.
    let mut bad = rotor_1903();
    bad.stations.swap(1, 2);
    assert_eq!(
        bemt_solve(&bad, RHO, 5.0, OMEGA_350).unwrap_err().code,
        "rotor-invalid"
    );
    jlog(
        "refusals",
        "\"gates\":\"stations cap/cap+1, V<0, omega=0, ordering\"",
    );
}

#[test]
fn bemt_golden_digest() {
    // The rail state: 350 rpm into the 10.73 m/s headwind (J ~ 0.71).
    let s = bemt_solve(&rotor_1903(), RHO, 10.73, OMEGA_350).unwrap();
    let mut payload = Vec::new();
    for v in [s.thrust_n, s.torque_nm, s.ct, s.cq, s.j] {
        payload.extend_from_slice(&v.to_bits().to_le_bytes());
    }
    for r in &s.stations {
        payload.extend_from_slice(&r.w_mps.to_bits().to_le_bytes());
    }
    let digest =
        fs_blake3::hash_domain("org.frankensim.fs-airscrew.v03-golden.v1", &payload).to_hex();
    jlog(
        "golden",
        &format!(
            "\"digest\":\"{digest}\",\"thrust\":{},\"j\":{}",
            s.thrust_n, s.j
        ),
    );
    assert_eq!(
        digest, "2c13d4e7efcc5091fa86b755255348e70620a8eead2a6ad6708ea6da6d8a3b91",
        "BEMT golden moved — determinism regression or an intentional \
         kernel change requiring the golden-bump protocol"
    );
}
