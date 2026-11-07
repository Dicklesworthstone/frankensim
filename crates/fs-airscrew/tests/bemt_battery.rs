//! V-03 core battery (bead wf-root-guzez.5.12.1, E4.5-i): the E1.6
//! HOLDOUT static anchor (the 1903 reconstruction at 350 rpm vs the
//! Wrights'/LFST bench numbers, factor-2 trend band — the reconstruction
//! is Estimated by the dossier's own rule), CT monotone in J, the eta
//! curve's rise-and-fall with a peak in the holdout band's neighborhood,
//! Prandtl limits, convergence receipts, spin-up toward the 350-rpm
//! operating class, caps, golden.
//! Repro: cargo test -p fs-airscrew --test bemt_battery

use fs_airscrew::drive::ShaftDrive;
use fs_airscrew::{
    BladeStation, MAX_STATIONS, Rotor, bemt_solve, engine_torque_at_prop_nm, rotor_spinup_step,
    rotor_spinup_step_with_work,
};

#[test]
fn shaft_drive_closes_source_transmission_and_bearing_power() {
    let drive = ShaftDrive::new(10000.0, 200.0, 4.0, 0.8, 2.0).unwrap();
    assert_eq!(drive.rated_prop_speed_rad_s(), 50.0);
    // Hand-calculated mechanical source, delivered torque and separate losses.
    // The last point is externally driven: bearing load exceeds engine torque.
    for (omega, command, source, torque, transmission, bearing) in [
        (0.0, 0.5, 0.0, 80.0, 0.0, 0.0),
        (25.0, 0.5, 2500.0, 30.0, 500.0, 1250.0),
        (50.0, 1.0, 10000.0, 60.0, 2000.0, 5000.0),
        (100.0, 1.0, 10000.0, -120.0, 2000.0, 20000.0),
        (75.0, 0.0, 0.0, -150.0, 0.0, 11250.0),
    ] {
        let point = drive.evaluate(omega, command).unwrap();
        assert_eq!(point.engine_speed_rad_s, 4.0 * omega);
        assert_eq!(point.engine_power_w, source);
        assert_eq!(point.prop_torque_nm, torque);
        assert_eq!(point.prop_power_w, torque * omega);
        assert!((point.transmission_heat_w - transmission).abs() < 1e-11);
        assert_eq!(point.bearing_heat_w, bearing);
        assert!(point.power_residual_w.abs() < 1e-11);
        assert!(point.transmission_heat_w >= 0.0);
        assert!(point.bearing_heat_w >= 0.0);
        assert!(point.engine_power_w <= command * drive.rated_power_w());
    }
}

#[test]
fn shaft_drive_has_continuous_torque_and_power_limits_and_correct_gearing() {
    let drive = ShaftDrive::new(1024.0, 64.0, 4.0, 1.0, 0.0).unwrap();
    for command in [0.0, 0.25, 1.0] {
        for omega in [0.0, 8.0, 16.0 - 1e-8, 16.0, 16.0 + 1e-8, 32.0, 64.0] {
            let point = drive.evaluate(omega, command).unwrap();
            let expected_torque = if omega <= 16.0 {
                command * 64.0
            } else {
                command * 1024.0 / omega
            };
            assert!((point.prop_torque_nm - expected_torque).abs() < 1e-13);
            assert!(point.prop_torque_nm <= command * 64.0);
            assert!(point.prop_power_w <= command * 1024.0 + 1e-12);
        }
    }
    let slow = drive.evaluate(8.0, 0.75).unwrap();
    let fast = ShaftDrive::new(1024.0, 64.0, 2.0, 1.0, 0.0)
        .unwrap()
        .evaluate(16.0, 0.75)
        .unwrap();
    assert_eq!(slow.engine_speed_rad_s, fast.engine_speed_rad_s);
    assert_eq!(slow.engine_torque_nm, fast.engine_torque_nm);
    assert_eq!(slow.engine_power_w, fast.engine_power_w);
    assert_eq!(slow.prop_power_w, fast.prop_power_w);
    assert_eq!(slow.prop_torque_nm, 2.0 * fast.prop_torque_nm);
}

#[test]
fn reference_drive_preserves_original_per_prop_torque_arithmetic() {
    let drive = ShaftDrive::wright_reference();
    let rated = 1025.0 / 60.0 * 2.0 * core::f64::consts::PI * (8.0 / 23.0);
    for omega in [0.0, 10.0, rated - 1e-8, rated, rated + 1e-8, 40.0, 120.0] {
        let original = if omega <= rated {
            8948.0 / rated / 2.0
        } else {
            8948.0 / omega / 2.0
        };
        assert_eq!(engine_torque_at_prop_nm(omega), original);
        assert_eq!(
            drive.evaluate(omega, 1.0).unwrap().prop_torque_nm / 2.0,
            original
        );
    }
}

#[test]
fn shaft_drive_refuses_unrepresentable_and_unsupported_inputs() {
    for parameters in [
        [0.0, 200.0, 4.0, 0.8, 2.0],
        [10000.0, -1.0, 4.0, 0.8, 2.0],
        [10000.0, 200.0, 0.0, 0.8, 2.0],
        [10000.0, 200.0, 4.0, 1.01, 2.0],
        [10000.0, 200.0, 4.0, 0.8, -1.0],
        [f64::NAN, 200.0, 4.0, 0.8, 2.0],
        [10000.0, f64::INFINITY, 4.0, 0.8, 2.0],
        [f64::MAX, f64::MIN_POSITIVE, 4.0, 0.8, 2.0],
    ] {
        let [power, speed, ratio, efficiency, drag] = parameters;
        assert_eq!(
            ShaftDrive::new(power, speed, ratio, efficiency, drag)
                .unwrap_err()
                .code,
            "shaft-drive-invalid"
        );
    }
    let drive = ShaftDrive::new(10000.0, 200.0, 4.0, 0.8, 2.0).unwrap();
    for (omega, command) in [
        (-1.0, 1.0),
        (40.0, -0.01),
        (40.0, 1.01),
        (f64::NAN, 1.0),
        (40.0, f64::NAN),
        (f64::MAX, 0.0),
    ] {
        assert_eq!(
            drive.evaluate(omega, command).unwrap_err().code,
            "shaft-drive-invalid"
        );
    }
    let tiny = ShaftDrive::new(1e-300, 1.0, 1.0, 1.0, 0.0).unwrap();
    assert_eq!(
        tiny.evaluate(1.0, 1e-100).unwrap_err().code,
        "shaft-drive-invalid"
    );
    assert_eq!(
        tiny.evaluate(1e-100, 1.0).unwrap_err().code,
        "shaft-drive-invalid"
    );
}

#[test]
fn held_torque_power_accounts_for_spinup_regeneration_and_reversal() {
    // These cases have exact binary-representable endpoints. Integrating a
    // held torque gives angle = omega_0 dt + (Q_engine-Q_load) dt²/(2I).
    for (inertia, omega, engine, load, dt) in [
        (2.0, 0.0, 8.0, 0.0, 0.5),
        (2.0, 10.0, -6.0, 2.0, 0.5),
        (2.0, 4.0, 0.0, -4.0, 0.5),
        (2.0, 1.0, -8.0, 0.0, 1.0),
        (2.0, -3.0, -8.0, -2.0, 0.5),
        (2.0, 10.0, 4.0, 4.0, 0.5),
    ] {
        let step = rotor_spinup_step_with_work(inertia, omega, engine, load, dt).unwrap();
        let acceleration = (engine - load) / inertia;
        let angle = omega * dt + 0.5 * acceleration * dt * dt;
        let endpoint = omega + acceleration * dt;
        let kinetic_change = 0.5 * inertia * (endpoint * endpoint - omega * omega);
        assert_eq!(step.rotation_rad, angle);
        assert_eq!(step.omega_after_rad_s, endpoint);
        assert_eq!(step.engine_work_j, engine * angle);
        assert_eq!(step.load_work_j, load * angle);
        assert_eq!(step.mean_engine_power_w * dt, step.engine_work_j);
        assert_eq!(step.mean_load_power_w * dt, step.load_work_j);
        assert_eq!(step.kinetic_energy_change_j, kinetic_change);
        assert_eq!(step.work_energy_residual_j, 0.0);
        assert_eq!(
            rotor_spinup_step(inertia, omega, engine, load, dt).unwrap(),
            endpoint
        );
    }
    let launch = rotor_spinup_step_with_work(2.0, 0.0, 8.0, 0.0, 0.5).unwrap();
    assert_eq!(launch.mean_engine_power_w, 8.0);
    assert_eq!(launch.engine_work_j, 4.0);
    let braking = rotor_spinup_step_with_work(2.0, 10.0, -6.0, 2.0, 0.5).unwrap();
    assert!(
        braking.engine_work_j < 0.0,
        "regeneration is not clipped to zero"
    );
    assert!(braking.load_work_j > 0.0);
}

#[test]
fn rotor_power_rejects_nonfinite_torques_inertia_and_unrepresentable_motion() {
    for inputs in [
        [f64::INFINITY, 1.0, 1.0, 0.0, 1.0],
        [1.0, 1.0, f64::NAN, 0.0, 1.0],
        [1.0, 1.0, 1.0, f64::INFINITY, 1.0],
        [1.0, 1.0, 1.0, 0.0, f64::INFINITY],
        [1.0, f64::INFINITY, 1.0, 0.0, 1.0],
        [0.0, 1.0, 1.0, 0.0, 1.0],
        [1.0, 1.0, 1.0, 0.0, 0.0],
        [f64::MIN_POSITIVE, 1.0, f64::MAX, 0.0, 1.0],
        [1.0, 1e200, 1e200, 0.0, 1.0],
    ] {
        let [inertia, omega, engine, load, dt] = inputs;
        assert_eq!(
            rotor_spinup_step_with_work(inertia, omega, engine, load, dt)
                .unwrap_err()
                .code,
            "rotor-dynamics-invalid"
        );
        assert!(rotor_spinup_step(inertia, omega, engine, load, dt).is_err());
    }
}

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
    // Golden bump 2026-10-28: the quadratic root is rationalized and
    // loads, momentum admission and station receipts now share the exact
    // accepted induction. The section law and Prandtl factors are unchanged.
    assert_eq!(
        digest, "8e5cc0361fd22d945c38b247b67ce89bd19db1c5b3cfc16ca3e8bec2ab3800e7",
        "BEMT golden moved — determinism regression or an intentional \
         kernel change requiring the golden-bump protocol"
    );
}
