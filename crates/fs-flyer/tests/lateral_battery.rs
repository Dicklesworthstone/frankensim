//! Lateral build-up battery (bead frankensim-4pa2k): the reduced
//! lateral tier under the admitted ReducedAeroelasticWarp envelope.
//!
//! Executed laws:
//! - adverse-yaw SIGN law: decoupled (1901) — the induced-drag yaw
//!   component OPPOSES the warp roll command on >99% of commanded
//!   ticks; coupled (1902+ linkage) turns the NET yaw proverse on
//!   >99% of commanded ticks (the E7.4b attribution, now engine-side);
//! - decomposition sum law: induced + rudder + profile == net exactly;
//! - zero command → zero lateral drift (no phantom dynamics);
//! - roll-cap crossing refuses with the typed envelope code (never a
//!   silent saturation);
//! - determinism: identical runs produce bit-identical state rows.
//!
//! Repro: cargo test -p fs-flyer --test lateral_battery

use fs_flyer::lateral::{
    FLYER_ROLL_INERTIA_KG_M2, LateralModel, LateralState, PHI_CAP_RAD, R_CAP_RAD_S, RudderLinkage,
};

fn jlog(case: &str, payload: &str) {
    println!("{{\"suite\":\"fs-flyer-lateral\",\"case\":\"{case}\",{payload}}}");
}

const DT: f64 = 1.0 / 120.0;
const RHO: f64 = 1.294;
const U_TRIM: f64 = 12.0;

fn run_ticks(
    model: &LateralModel,
    twist_per_tick: f64,
    ticks: usize,
) -> Vec<Result<fs_flyer::lateral::LateralTick, fs_flyer::Refusal>> {
    let mut state = LateralState::default();
    (0..ticks)
        .map(|_| model.step(&mut state, twist_per_tick, U_TRIM, RHO, DT))
        .collect()
}

#[test]
fn decoupled_induced_yaw_opposes_the_warp_command() {
    let model = LateralModel::wright_v1(RudderLinkage::Decoupled);
    let rows = run_ticks(&model, 0.10, 400);
    let mut commanded = 0usize;
    let mut adverse = 0usize;
    for row in &rows {
        let tick = row.as_ref().expect("decoupled run stays inside the band");
        if tick.loaded_twist_rad > 1e-9 {
            commanded += 1;
            if tick.yaw.induced_drag_yaw_nm < 0.0 {
                adverse += 1;
            }
        }
    }
    assert!(commanded > 300, "commanded ticks {commanded}");
    let share = adverse as f64 / commanded as f64;
    jlog(
        "adverse_sign",
        &format!("\"commanded\":{commanded},\"adverse\":{adverse},\"share\":{share:.4}"),
    );
    assert!(
        share > 0.99,
        "induced-drag yaw must oppose the warp roll command (share {share})"
    );
}

#[test]
fn linked_rudder_turns_net_yaw_proverse() {
    // The declared linkage gain exceeds |induced| at the reference q so
    // the NET moment follows the command (proverse), mirroring the
    // closed E7.4b attribution for the 1902+ twin.
    let gain = 42.0 * (12.0_f64.powi(2) / 2.0 * RHO).max(76.8) / 76.8 + 6.0;
    let model = LateralModel::wright_v1(RudderLinkage::Linked {
        gain_nm_per_rad: gain,
    });
    let rows = run_ticks(&model, 0.10, 400);
    let mut commanded = 0usize;
    let mut proverse = 0usize;
    for row in &rows {
        let tick = row.as_ref().expect("linked run stays inside the band");
        if tick.loaded_twist_rad > 1e-9 {
            commanded += 1;
            if tick.yaw.net() > 0.0 {
                proverse += 1;
            }
        }
    }
    let share = proverse as f64 / commanded as f64;
    jlog(
        "proverse_net",
        &format!("\"commanded\":{commanded},\"proverse\":{proverse},\"share\":{share:.4}"),
    );
    assert!(share > 0.99, "net yaw proverse share {share}");
}

#[test]
fn decomposition_sums_to_net_exactly() {
    let model = LateralModel::wright_v1(RudderLinkage::Linked {
        gain_nm_per_rad: 30.0,
    });
    for row in run_ticks(&model, -0.07, 250) {
        let tick = row.expect("in-band run");
        let sum = tick.yaw.induced_drag_yaw_nm + tick.yaw.rudder_yaw_nm + tick.yaw.profile_yaw_nm;
        assert_eq!(
            sum,
            tick.yaw.net(),
            "the published order is the float sum order"
        );
    }
}

#[test]
fn zero_command_produces_no_lateral_drift() {
    let model = LateralModel::wright_v1(RudderLinkage::Decoupled);
    let rows = run_ticks(&model, 0.0, 600);
    for row in rows {
        let tick = row.expect("zero-command run is trivially in-band");
        assert_eq!(tick.state, LateralState::default(), "no phantom dynamics");
    }
}

#[test]
fn sustained_full_warp_refuses_at_the_roll_cap() {
    let model = LateralModel::wright_v1(RudderLinkage::Decoupled);
    let mut state = LateralState::default();
    let mut refused = None;
    for i in 0..100_000 {
        match model.step(&mut state, 0.35, U_TRIM, RHO, DT) {
            Ok(_) => {}
            Err(refusal) => {
                refused = Some((i, refusal));
                break;
            }
        }
    }
    let (tick_index, refusal) = refused.expect("sustained warp must reach the cap");
    assert_eq!(refusal.code, "lateral-envelope-exceeded");
    assert!(
        state.phi_rad.abs() <= PHI_CAP_RAD + DT,
        "refusal fires at the boundary, not beyond it"
    );
    jlog(
        "cap_crossing",
        &format!("\"tick\":{tick_index},\"phi\":{:.4}", state.phi_rad),
    );
}

#[test]
fn yaw_rate_cap_refuses_without_mutating_state() {
    let model = LateralModel::wright_v1(RudderLinkage::Linked {
        gain_nm_per_rad: 1.0e9,
    });
    let mut state = LateralState {
        r_rad_s: R_CAP_RAD_S - 0.001,
        ..LateralState::default()
    };
    let before = state;
    let refusal = model
        .step(&mut state, 0.01, U_TRIM, RHO, DT)
        .expect_err("a yaw-rate cap crossing must fail closed");
    assert_eq!(refusal.code, "lateral-envelope-exceeded");
    assert_eq!(
        state, before,
        "a refused step must not clamp or partially commit"
    );
}

#[test]
fn identical_runs_are_bit_identical() {
    let model = LateralModel::wright_v1(RudderLinkage::Decoupled);
    let to_bits = |rows: &[Result<fs_flyer::lateral::LateralTick, fs_flyer::Refusal>]| -> Vec<u64> {
        rows.iter()
            .map(|r| {
                let t = r.as_ref().expect("in-band");
                t.state.phi_rad.to_bits()
                    ^ t.state.p_rad_s.to_bits()
                    ^ t.state.psi_rad.to_bits()
                    ^ t.state.r_rad_s.to_bits()
            })
            .collect()
    };
    let a = to_bits(&run_ticks(&model, 0.05, 500));
    let b = to_bits(&run_ticks(&model, 0.05, 500));
    assert_eq!(a, b, "deterministic integration is bit-stable");
}

#[test]
fn roll_damping_limits_the_uncommanded_rate_growth() {
    // After the command stops, damping must drive |p| down, not up.
    let model = LateralModel::wright_v1(RudderLinkage::Decoupled);
    let mut state = LateralState::default();
    for _ in 0..240 {
        model
            .step(&mut state, 0.20, U_TRIM, RHO, DT)
            .expect("in-band ramp");
    }
    let p_after_command = state.p_rad_s.abs();
    for _ in 0..360 {
        model
            .step(&mut state, 0.0, U_TRIM, RHO, DT)
            .expect("coasting stays in-band");
    }
    jlog(
        "damping",
        &format!(
            "\"p_cmd\":{p_after_command:.5},\"p_coast\":{:.5}",
            state.p_rad_s.abs()
        ),
    );
    assert!(
        state.p_rad_s.abs() < p_after_command,
        "roll damping must decay the rate once the command stops"
    );
    let _ = FLYER_ROLL_INERTIA_KG_M2;
}
