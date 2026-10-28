//! E4.6b-i battery (bead wf-root-guzez.5.14.1): pilot-command sign
//! chain, stop engagement at cap AND one ulp past, regularized friction
//! (finite at zero rate; sub-Coulomb hold), per-step dissipation
//! oracles, backdrivability, dt-refinement convergence of the midpoint
//! integrator (executed, not assumed), prior-band admission at both
//! edges, determinism, golden.
//! Repro: cargo test -p fs-flyer --test canardmech_battery

use fs_flyer::canardmech::{CANARD_MECH_V1, CanardMechanism, HINGE_AXIS_PRIOR_PCT, MechState};

fn jlog(case: &str, payload: &str) {
    println!("{{\"suite\":\"fs-flyer-e46bi\",\"case\":\"{case}\",{payload}}}");
}

const REST: MechState = MechState {
    delta_rad: 0.0,
    rate_rad_s: 0.0,
};

#[test]
fn pilot_pull_back_drives_positive_deflection() {
    // control-signs-v1: positive pilot force (pull back) -> positive
    // canard command -> nose-up. The mechanism must carry the sign.
    let m = CANARD_MECH_V1;
    let mut st = REST;
    for _ in 0..200 {
        st = m.step(st, 0.0, 40.0, 0.005).unwrap().0;
    }
    assert!(
        st.delta_rad > 0.05,
        "pull-back must drive + deflection: {}",
        st.delta_rad
    );
    jlog("sign-chain", &format!("\"delta\":{}", st.delta_rad));
}

#[test]
fn stops_engage_at_cap_and_one_ulp_past() {
    let m = CANARD_MECH_V1;
    // AT the stop, zero rate: no stop torque (the stop is a strict
    // penetration spring).
    let at = MechState {
        delta_rad: m.stop_rad,
        rate_rad_s: 0.0,
    };
    let (_, r_at) = m.step(at, 0.0, 0.0, 0.001).unwrap();
    assert!(
        r_at.stop_torque_nm.abs() < 1e-15,
        "no torque exactly at the stop: {}",
        r_at.stop_torque_nm
    );
    // One ulp past: torque engages, negative (restoring).
    let past = MechState {
        delta_rad: m.stop_rad.next_up(),
        rate_rad_s: 0.0,
    };
    let (_, r_past) = m.step(past, 0.0, 0.0, 0.001).unwrap();
    assert!(
        r_past.stop_torque_nm < 0.0,
        "stop must engage one ulp past: {}",
        r_past.stop_torque_nm
    );
    // Dynamic: drive hard into the stop; deflection must stay bounded
    // near the stop (penetration < 2 deg) and settle.
    let mut st = REST;
    for _ in 0..4000 {
        st = m.step(st, 60.0, 100.0, 0.001).unwrap().0;
    }
    assert!(
        st.delta_rad < m.stop_rad + 0.035,
        "stop containment failed: {} vs stop {}",
        st.delta_rad,
        m.stop_rad
    );
    assert!(st.rate_rad_s.abs() < 0.05, "must settle on the stop");
    jlog(
        "stops",
        &format!("\"settled\":{},\"stop\":{}", st.delta_rad, m.stop_rad),
    );
}

#[test]
fn friction_is_regular_at_zero_rate_and_holds_sub_coulomb_loads() {
    let m = CANARD_MECH_V1;
    // Zero rate: finite torque, no NaN anywhere.
    let (st, r) = m.step(REST, 0.5, 0.0, 0.002).unwrap();
    assert!(st.delta_rad.is_finite() && st.rate_rad_s.is_finite());
    assert!(r.net_torque_nm.is_finite());
    // A load below the Coulomb level barely moves the surface (the
    // regularized stiction creep must stay tiny over a full second).
    let mut st = REST;
    for _ in 0..500 {
        st = m.step(st, 0.5 * m.coulomb_nm, 0.0, 0.002).unwrap().0;
    }
    assert!(
        st.delta_rad.abs() < 0.02,
        "sub-Coulomb load must essentially hold: {}",
        st.delta_rad
    );
    // And a load ABOVE Coulomb moves it decisively (liveness twin).
    let mut st2 = REST;
    for _ in 0..500 {
        st2 = m.step(st2, 3.0 * m.coulomb_nm, 0.0, 0.002).unwrap().0;
    }
    assert!(
        st2.delta_rad > 10.0 * st.delta_rad.abs().max(1e-6),
        "super-Coulomb liveness: {} vs {}",
        st2.delta_rad,
        st.delta_rad
    );
    jlog(
        "friction",
        &format!("\"creep\":{},\"drive\":{}", st.delta_rad, st2.delta_rad),
    );
}

#[test]
fn per_step_dissipation_is_nonnegative_under_oscillatory_drive() {
    let m = CANARD_MECH_V1;
    let mut st = REST;
    let mut worst = f64::INFINITY;
    for i in 0..2000 {
        let drive = 25.0 * (0.05 * f64::from(i)).sin();
        let (next, r) = m.step(st, drive, 0.0, 0.002).unwrap();
        // PER-STEP oracle, never totals: friction can only dissipate.
        assert!(
            r.friction_dissipation_j >= -1e-12,
            "step {i}: negative friction dissipation {}",
            r.friction_dissipation_j
        );
        worst = worst.min(r.friction_dissipation_j);
        st = next;
    }
    jlog("dissipation", &format!("\"worst_step_j\":{worst}"));
}

#[test]
fn released_surface_is_backdrivable() {
    // Released pilot (force 0), constant aero hinge moment: the surface
    // must accelerate away — the mechanism is backdrivable (Orville's
    // overcontrol account REQUIRES aero moments to drive the surface).
    let m = CANARD_MECH_V1;
    let mut st = REST;
    for _ in 0..300 {
        st = m.step(st, 8.0, 0.0, 0.002).unwrap().0;
    }
    assert!(
        st.delta_rad > 0.02,
        "aero moment must backdrive the released surface: {}",
        st.delta_rad
    );
    jlog("backdrive", &format!("\"delta\":{}", st.delta_rad));
}

#[test]
fn midpoint_refinement_converges_second_order() {
    // Same 0.2 s trajectory at dt and dt/2: the end-state error against
    // a fine reference must shrink ~4x (executed convergence, not an
    // assumed property of the integrator).
    let m = CANARD_MECH_V1;
    let run = |dt: f64| -> f64 {
        let n = (0.2 / dt).round() as usize;
        let mut st = REST;
        for i in 0..n {
            let drive = 15.0 * (10.0 * dt * i as f64).sin();
            st = m.step(st, drive, 0.0, dt).unwrap().0;
        }
        st.delta_rad
    };
    let fine = run(0.000_125);
    let e1 = (run(0.002) - fine).abs();
    let e2 = (run(0.001) - fine).abs();
    // The oscillatory drive crosses the stiction knee repeatedly, which
    // costs local order at the tanh regularization — CONVERGENCE is the
    // claim here (measured ratio ~2.1), not clean 2nd order.
    assert!(
        e2 < e1 / 1.8,
        "midpoint refinement must converge through the knee: e1 {e1:e}, e2 {e2:e}"
    );
    // On a SMOOTH trajectory (constant drive, rate held well above the
    // regularization velocity after startup) the order must be clean.
    let smooth = |dt: f64| -> f64 {
        let n = (0.1 / dt).round() as usize;
        let mut st = REST;
        for _ in 0..n {
            st = m.step(st, 20.0, 0.0, dt).unwrap().0;
        }
        st.delta_rad
    };
    let sf = smooth(0.000_125);
    let s1 = (smooth(0.002) - sf).abs();
    let s2 = (smooth(0.001) - sf).abs();
    assert!(
        s2 < s1 / 3.0,
        "smooth-regime refinement must be ~2nd order: s1 {s1:e}, s2 {s2:e}"
    );
    jlog(
        "refinement",
        &format!("\"knee\":[{e1:e},{e2:e}],\"smooth\":[{s1:e},{s2:e}]"),
    );
}

#[test]
fn prior_band_admission_at_both_edges() {
    let mk = |pct: f64| -> CanardMechanism {
        CanardMechanism {
            hinge_axis_pct_chord: pct,
            ..CANARD_MECH_V1
        }
    };
    assert!(mk(HINGE_AXIS_PRIOR_PCT.0).admit().is_ok(), "lo edge admits");
    assert!(mk(HINGE_AXIS_PRIOR_PCT.1).admit().is_ok(), "hi edge admits");
    assert_eq!(
        mk(HINGE_AXIS_PRIOR_PCT.0.next_down())
            .admit()
            .unwrap_err()
            .code,
        "hinge-axis-outside-prior"
    );
    assert_eq!(
        mk(HINGE_AXIS_PRIOR_PCT.1.next_up())
            .admit()
            .unwrap_err()
            .code,
        "hinge-axis-outside-prior"
    );
    // Mechanism-parameter refusals.
    let bad = CanardMechanism {
        inertia_kg_m2: 0.0,
        ..CANARD_MECH_V1
    };
    assert_eq!(bad.admit().unwrap_err().code, "canard-mech-invalid");
    assert_eq!(
        CANARD_MECH_V1
            .step(REST, f64::NAN, 0.0, 0.001)
            .unwrap_err()
            .code,
        "mech-state-invalid"
    );
    assert_eq!(
        CANARD_MECH_V1.step(REST, 0.0, 0.0, 0.0).unwrap_err().code,
        "mech-state-invalid"
    );
    jlog("priors", "\"band_edges_and_ulps\":true");
}

#[test]
fn determinism_and_golden() {
    let m = CANARD_MECH_V1;
    let run = || -> Vec<u64> {
        let mut st = REST;
        let mut out = Vec::new();
        for i in 0..400 {
            let drive = 20.0 * (0.03 * f64::from(i)).sin();
            let pilot = 30.0 * (0.011 * f64::from(i)).cos();
            st = m.step(st, drive, pilot, 0.004).unwrap().0;
            out.push(st.delta_rad.to_bits());
        }
        out
    };
    let a = run();
    assert_eq!(a, run(), "bitwise repeat");
    let mut payload = Vec::new();
    for b in &a {
        payload.extend_from_slice(&b.to_le_bytes());
    }
    let digest =
        fs_blake3::hash_domain("org.frankensim.fs-flyer.e46bi-golden.v1", &payload).to_hex();
    jlog("golden", &format!("\"digest\":\"{digest}\""));
    assert_eq!(
        // 2026-10-28: converged passive friction/stop solve replaces twelve
        // unchecked Picard substitutions; the physical trajectory changes.
        digest,
        "4d12178e5ff2df74427cf7d4b1ca6b97c65bc6d84901d64b40185275460db4a9",
        "canard-mechanism golden moved — determinism regression or an \
         intentional model change requiring the golden-bump protocol"
    );
}

fn mechanical_energy(m: CanardMechanism, state: MechState) -> f64 {
    let penetration = (state.delta_rad.abs() - m.stop_rad).max(0.0);
    0.5 * m.inertia_kg_m2 * state.rate_rad_s.powi(2)
        + 0.5 * m.stop_stiffness_nm_per_rad * penetration.powi(2)
}

#[test]
fn stiff_friction_solves_the_positive_sub_coulomb_drift() {
    // The old twelve Picard substitutions reversed the motion under this
    // positive load and reported negative dissipation. Independently solve
    // the scalar momentum balance with the standard-library tanh.
    let m = CanardMechanism {
        friction_reg_rad_s: 1e-5,
        ..CANARD_MECH_V1
    };
    let dt = 0.005;
    let drive = 1.0;
    let (next, receipt) = m.step(REST, drive, 0.0, dt).unwrap();
    let mut lower = 0.0;
    let mut upper = drive * dt / (2.0 * m.inertia_kg_m2);
    for _ in 0..100 {
        let middle = 0.5 * (lower + upper);
        let residual = 2.0 * m.inertia_kg_m2 * middle / dt
            + m.coulomb_nm * (middle / m.friction_reg_rad_s).tanh()
            + m.viscous_nm_s * middle
            - drive;
        if residual < 0.0 {
            lower = middle;
        } else {
            upper = middle;
        }
    }
    let expected_rate = lower + upper;
    assert!(next.delta_rad > 0.0 && next.rate_rad_s > 0.0);
    assert!((next.rate_rad_s - expected_rate).abs() < 1e-13);
    assert!(receipt.friction_dissipation_j > 0.0);
    assert_eq!(receipt.stop_dissipation_j, 0.0);
    let defect =
        mechanical_energy(m, next) + receipt.friction_dissipation_j - drive * next.delta_rad;
    assert!(
        defect.abs() < 1e-15,
        "stiff friction work defect {defect:e}"
    );
}

#[test]
fn fully_engaged_stop_retains_the_linear_midpoint_solution() {
    let m = CanardMechanism {
        coulomb_nm: 0.0,
        ..CANARD_MECH_V1
    };
    let before = MechState {
        delta_rad: m.stop_rad + 0.04,
        rate_rad_s: 0.2,
    };
    let dt = 0.01;
    let drive = 3.0;
    let expected_midpoint = (2.0 * m.inertia_kg_m2 * before.rate_rad_s / dt + drive
        - m.stop_stiffness_nm_per_rad * (before.delta_rad - m.stop_rad))
        / (2.0 * m.inertia_kg_m2 / dt
            + m.viscous_nm_s
            + m.stop_damping_nm_s
            + 0.5 * dt * m.stop_stiffness_nm_per_rad);
    let (after, receipt) = m.step(before, drive, 0.0, dt).unwrap();
    assert!(after.delta_rad > m.stop_rad);
    assert!((after.rate_rad_s - (2.0 * expected_midpoint - before.rate_rad_s)).abs() < 1e-13);
    assert!((after.delta_rad - before.delta_rad - dt * expected_midpoint).abs() < 1e-15);
    let defect = mechanical_energy(m, after) - mechanical_energy(m, before)
        + receipt.friction_dissipation_j
        + receipt.stop_dissipation_j
        - receipt.drive_work_j;
    assert!(defect.abs() < 1e-11, "engaged-stop work defect {defect:e}");
}

#[test]
fn stop_entry_and_release_close_mechanical_work_on_both_sides() {
    let m = CANARD_MECH_V1;
    let dt = 0.02;
    for sign in [-1.0, 1.0] {
        for (offset, rate) in [(-0.005, 1.0), (0.01, -1.0)] {
            let before = MechState {
                delta_rad: sign * (m.stop_rad + offset),
                rate_rad_s: sign * rate,
            };
            let (after, receipt) = m.step(before, 0.0, 0.0, dt).unwrap();
            let penetration_before = (before.delta_rad.abs() - m.stop_rad).max(0.0);
            let penetration_after = (after.delta_rad.abs() - m.stop_rad).max(0.0);
            if offset < 0.0 {
                assert!(penetration_after > 0.0, "entry must cross the stop");
            } else {
                assert!(
                    after.delta_rad.abs() < m.stop_rad,
                    "release must cross the stop"
                );
            }
            // Independent path-length oracle: during either one-sided crossing
            // the distance spent in contact is exactly the nonzero penetration.
            let contact_distance = penetration_before + penetration_after;
            let midpoint_rate = 0.5 * (before.rate_rad_s + after.rate_rad_s);
            let expected_stop_loss = m.stop_damping_nm_s * midpoint_rate.abs() * contact_distance;
            assert!((receipt.stop_dissipation_j - expected_stop_loss).abs() < 1e-12);
            assert!(receipt.friction_dissipation_j >= 0.0 && receipt.stop_dissipation_j >= 0.0);
            let defect = mechanical_energy(m, after) - mechanical_energy(m, before)
                + receipt.friction_dissipation_j
                + receipt.stop_dissipation_j;
            assert!(defect.abs() < 1e-11, "crossing work defect {defect:e}");
        }
    }
}

#[test]
fn stiff_unforced_hinge_dissipates_without_invalid_steps() {
    let m = CanardMechanism {
        friction_reg_rad_s: 1e-5,
        ..CANARD_MECH_V1
    };
    let mut state = MechState {
        delta_rad: m.stop_rad + 0.01,
        rate_rad_s: 0.4,
    };
    for _ in 0..100 {
        let (next, receipt) = m.step(state, 0.0, 0.0, 0.05).unwrap();
        let before = mechanical_energy(m, state);
        let after = mechanical_energy(m, next);
        assert!(after <= before + 1e-11);
        assert!(receipt.friction_dissipation_j >= 0.0 && receipt.stop_dissipation_j >= 0.0);
        let defect = after - before + receipt.friction_dissipation_j + receipt.stop_dissipation_j;
        assert!(defect.abs() < 1e-11, "passive work defect {defect:e}");
        state = next;
    }
    assert_eq!(
        m.step(REST, f64::MAX, f64::MAX, 0.1).unwrap_err().code,
        "mech-step-non-finite"
    );
}
