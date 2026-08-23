//! E4.6c-iii battery (bead wf-root-guzez.5.16.3): the with-vs-without
//! stabilization comparison on the REAL unstable airframe loop
//! (executed), per-tick authority clamp at the cap AND one ulp past,
//! HUD-flag liveness (active even at zero output), the historical-
//! calibration isolation refusal (executed), determinism, golden.
//! Repro: cargo test -p fs-flyer --test assist_battery

use fs_flyer::aircraft::wright_openloop_v1;
use fs_flyer::assist::{ASSIST_V1, AssistSpec, MAX_AUTHORITY_FRAC, historical_calibration_admit};
use fs_flyer::canardmech::{CANARD_MECH_V1, MechState};
use fs_flyer::longitudinal::{IYY_KG_M2, linearize};
use fs_flyer::perception::perception_v1;
use fs_flyer::pilot::{PilotWrightModel, pack_cues};

fn jlog(case: &str, payload: &str) {
    println!("{{\"suite\":\"fs-flyer-e46ciii\",\"case\":\"{case}\",{payload}}}");
}

const DT: f64 = 1.0 / 120.0;

#[test]
fn with_vs_without_on_the_real_unstable_loop() {
    // The whole point of the assist: the piloted loop that DIVERGES via
    // overcontrol (H-02c record) becomes bounded with the SAS+assist
    // engaged. Same airframe, same pilot member, same latency chain.
    let d = wright_openloop_v1();
    let rho = 1.294;
    let trim = d.trim(rho, [13.0, 0.06, 0.1, 45.0]).unwrap();
    let rep = linearize(&d, &trim, rho).unwrap();
    let a = rep.a;
    let h = 0.004;
    let bp = d
        .force_buildup(
            trim.v_mps,
            trim.alpha_rad,
            trim.delta_canard_rad + h,
            trim.omega_prop_rad_s,
            0.0,
            rho,
        )
        .unwrap();
    let bm = d
        .force_buildup(
            trim.v_mps,
            trim.alpha_rad,
            trim.delta_canard_rad - h,
            trim.omega_prop_rad_s,
            0.0,
            rho,
        )
        .unwrap();
    let m = d.gross_mass_kg;
    let b = [
        (bp.force_n[0] - bm.force_n[0]) / (2.0 * h) / m,
        (bp.force_n[2] - bm.force_n[2]) / (2.0 * h) / m,
        (bp.moment_y_nm - bm.moment_y_nm) / (2.0 * h) / IYY_KG_M2,
        0.0,
    ];
    let horizon = (15.0 / DT) as usize;
    let escape = 0.35f64;
    let travel = CANARD_MECH_V1.stop_rad;
    let sim = |assist: Option<AssistSpec>| -> (f64, f64, bool) {
        let mut perc = perception_v1(5);
        for c in &mut perc.cues {
            c.remnant_sigma = 0.0;
        }
        perc.cues[0].delay_ticks = 8;
        perc.cues[1].delay_ticks = 5;
        let mut pilot = PilotWrightModel::new(2, 9).unwrap();
        pilot.gains.remnant_sigma_force_n = 0.0;
        pilot.gains.remnant_sigma_warp = 0.0;
        let mut ps = perc.init().unwrap();
        let mut st = pilot.init().unwrap();
        let mut mech = MechState {
            delta_rad: 0.0,
            rate_rad_s: 0.0,
        };
        let mut x = [0.0, 0.0, 0.0, 0.05f64];
        let mut max_theta = 0.0f64;
        let mut flag_always = true;
        for step in 0..horizon {
            let theta = x[3];
            max_theta = max_theta.max(theta.abs());
            if theta.abs() > escape {
                return (step as f64 * DT, max_theta, flag_always);
            }
            let wdot: f64 = (0..4).map(|j| a[1][j] * x[j]).sum::<f64>();
            let cues = perc
                .step(&mut ps, pack_cues(theta, x[2], wdot, 0.0, 0.0, 0.0))
                .unwrap();
            let cmd = pilot
                .step(&mut st, &cues, mech.delta_rad, mech.rate_rad_s, 0.0, 0.0)
                .unwrap();
            mech = CANARD_MECH_V1
                .step(mech, 0.0, cmd.lever_force_n, DT)
                .unwrap()
                .0;
            let mut dc = mech.delta_rad;
            if let Some(sp) = &assist {
                let out = sp.apply(x[2], theta, travel).unwrap();
                flag_always &= out.active;
                dc += out.dc_assist_rad;
            }
            for i in 0..4 {
                let mut dx: f64 = (0..4).map(|j| a[i][j] * x[j]).sum::<f64>();
                dx += b[i] * dc;
                x[i] += DT * dx;
            }
        }
        (15.0, max_theta, flag_always)
    };
    let without = sim(None);
    let with = sim(Some(ASSIST_V1));
    jlog(
        "with-vs-without",
        &format!(
            "\"without_t_escape_s\":{},\"with_t_escape_s\":{},\"with_max_theta\":{},\"flag_always\":{}",
            without.0, with.0, with.1, with.2
        ),
    );
    assert!(without.0 < 15.0, "the unassisted loop must diverge (H-02c)");
    assert!(
        with.0 >= 15.0 && with.1 < 0.2,
        "the assisted loop must stay bounded: t {}, max {}",
        with.0,
        with.1
    );
    assert!(with.2, "the HUD flag must be set EVERY assisted tick");
}

#[test]
fn authority_clamps_exactly_at_the_fraction_of_travel() {
    let travel = 0.5236f64;
    let authority = ASSIST_V1.authority_frac * travel;
    // A huge rate must clamp EXACTLY at the authority with the flag.
    let big = ASSIST_V1.apply(50.0, 0.0, travel).unwrap();
    assert_eq!(big.dc_assist_rad, -authority);
    assert!(big.clamped);
    // Just inside: unclamped (per-item check of the boundary).
    let inside_q = (authority - 1e-9) / ASSIST_V1.sas_rate_gain;
    let inside = ASSIST_V1.apply(inside_q, 0.0, travel).unwrap();
    assert!(!inside.clamped);
    assert!((inside.dc_assist_rad + authority - 1e-9).abs() < 1e-12);
    jlog("authority", &format!("\"authority_rad\":{authority}"));
}

#[test]
fn flag_is_active_even_at_zero_output() {
    // Visibility law: if the system is engaged, the flag is ON — even
    // when the current correction happens to be zero.
    let out = ASSIST_V1.apply(0.0, 0.0, 0.5).unwrap();
    assert_eq!(out.dc_assist_rad, 0.0);
    assert!(out.active, "engaged system must always flag");
    assert!(!out.clamped);
    jlog("flag", "\"active_at_zero_output\":true");
}

#[test]
fn historical_calibration_refuses_assist_tagged_specs() {
    // The ISOLATION law, executed: anything carrying the calibration-
    // subset tag is refused by the historical-calibration admission.
    let err = historical_calibration_admit(&ASSIST_V1).unwrap_err();
    assert_eq!(err.code, "assist-in-historical-calibration");
    assert!(err.message.contains("CalibrationSubsetTag"));
    jlog("isolation", "\"historical_refusal\":true");
}

#[test]
fn caps_at_cap_and_cap_plus_one() {
    let mk = |frac: f64| AssistSpec {
        authority_frac: frac,
        ..ASSIST_V1
    };
    assert!(mk(MAX_AUTHORITY_FRAC).admit().is_ok(), "cap admits");
    assert_eq!(
        mk(MAX_AUTHORITY_FRAC.next_up()).admit().unwrap_err().code,
        "assist-spec-invalid",
        "cap+1 refuses"
    );
    assert_eq!(mk(0.0).admit().unwrap_err().code, "assist-spec-invalid");
    let neg = AssistSpec {
        sas_rate_gain: -0.1,
        ..ASSIST_V1
    };
    assert_eq!(neg.admit().unwrap_err().code, "assist-spec-invalid");
    assert_eq!(
        ASSIST_V1.apply(f64::NAN, 0.0, 0.5).unwrap_err().code,
        "assist-input-invalid"
    );
    assert_eq!(
        ASSIST_V1.apply(0.0, 0.0, 0.0).unwrap_err().code,
        "assist-input-invalid"
    );
    jlog("caps", "\"cap_and_cap_plus_one\":true");
}

#[test]
fn determinism_and_golden() {
    let mut payload = Vec::new();
    for i in 0..200 {
        let q = 0.3 * fs_math::det::sin(0.1 * f64::from(i));
        let th = 0.1 * fs_math::det::cos(0.07 * f64::from(i));
        let out = ASSIST_V1.apply(q, th, 0.5236).unwrap();
        let out2 = ASSIST_V1.apply(q, th, 0.5236).unwrap();
        assert_eq!(out, out2, "bitwise repeat");
        payload.extend_from_slice(&out.dc_assist_rad.to_bits().to_le_bytes());
    }
    let digest =
        fs_blake3::hash_domain("org.frankensim.fs-flyer.e46ciii-golden.v1", &payload).to_hex();
    jlog("golden", &format!("\"digest\":\"{digest}\""));
    assert_eq!(
        digest, "27ffdac75e9db02bf4bb810fd2df7abb65a3684f2925f2ebe3f3ccfad87b94ff",
        "assist golden moved — determinism regression or an intentional \
         model change requiring the golden-bump protocol"
    );
}
