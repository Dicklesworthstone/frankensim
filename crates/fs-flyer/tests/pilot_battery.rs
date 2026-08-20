//! V-02c1 battery + H-02c record (bead wf-root-guzez.5.16.2, E4.6c-ii):
//! the GENERIC pilot-vehicle mechanism on K/s and K/s² controlled
//! elements (stability, damped-frequency band, delay margin with a live
//! divergence falsifier), saturation at the limits with receipt flags,
//! registered-family caps at cap AND cap+1, deterministic remnant, and
//! the H-02c closed-loop record against the REAL linearized airframe +
//! canard mechanism (compatibility at Estimated — a receipt, never a
//! validation claim). Repro: cargo test -p fs-flyer --test pilot_battery

use fs_flyer::aircraft::wright_openloop_v1;
use fs_flyer::canardmech::{CANARD_MECH_V1, MechState};
use fs_flyer::longitudinal::{IYY_KG_M2, linearize};
use fs_flyer::perception::{N_CUES, PerceptionModelSpec, perception_v1};
use fs_flyer::pilot::{FAMILY_SIZE, PilotWrightModel, pack_cues, pilot_family_v1};

fn jlog(case: &str, payload: &str) {
    println!("{{\"suite\":\"fs-flyer-v02c1\",\"case\":\"{case}\",{payload}}}");
}

const DT: f64 = 1.0 / 120.0;

fn quiet_perception(seed: u64) -> PerceptionModelSpec {
    let mut p = perception_v1(seed);
    for c in &mut p.cues {
        c.remnant_sigma = 0.0;
    }
    p
}

fn quiet_pilot(member: u32) -> PilotWrightModel {
    let mut m = PilotWrightModel::new(member, 9).unwrap();
    m.gains.remnant_sigma_force_n = 0.0;
    m.gains.remnant_sigma_warp = 0.0;
    m
}

/// Closed loop of the pilot's LATERAL channel around φ'' = k_el·u
/// (K/s², the crossover-model workhorse). Returns the |φ| trace.
fn run_double_integrator(pilot: &PilotWrightModel, k_el: f64, seconds: f64) -> Vec<f64> {
    let perc = quiet_perception(3);
    let mut ps = perc.init().unwrap();
    let mut st = pilot.init().unwrap();
    let (mut phi, mut p) = (0.1f64, 0.0f64);
    let mut trace = Vec::new();
    for _ in 0..((seconds / DT) as usize) {
        let cues = perc
            .step(&mut ps, pack_cues(0.0, 0.0, 0.0, phi, p, 0.0))
            .unwrap();
        let cmd = pilot.step(&mut st, &cues, 0.0, 0.0, 0.0, 0.0).unwrap();
        let u = cmd.warp_cmd_rad;
        p += DT * k_el * u;
        phi += DT * p;
        trace.push(phi.abs());
        if !phi.is_finite() || phi.abs() > 50.0 {
            trace.push(f64::INFINITY);
            break;
        }
    }
    trace
}

#[test]
fn double_integrator_stabilizes_with_damped_frequency_in_band() {
    // V-02c1 crossover clause: on K/s² the pilot's rate lead makes the
    // loop stable; the damped oscillation frequency sits in the
    // human-crossover class.
    let pilot = quiet_pilot(0);
    let trace = run_double_integrator(&pilot, 2.0, 40.0);
    let tail = &trace[trace.len() - 240..];
    assert!(
        tail.iter().all(|v| v.is_finite() && *v < 0.02),
        "K/s^2 loop must settle: tail max {:?}",
        tail.iter().cloned().fold(0.0f64, f64::max)
    );
    // Damped frequency from the first few peaks of |phi|.
    let mut peaks = Vec::new();
    for i in 1..(trace.len() - 1) {
        if trace[i] > trace[i - 1] && trace[i] > trace[i + 1] && trace[i] > 0.005 {
            peaks.push(i as f64 * DT);
        }
    }
    assert!(peaks.len() >= 2, "need visible oscillation to measure");
    let period = 2.0 * (peaks[1] - peaks[0]); // |phi| peaks twice per cycle
    let omega_d = core::f64::consts::TAU / period;
    assert!(
        (0.5..8.0).contains(&omega_d),
        "damped frequency {omega_d} rad/s outside the crossover class"
    );
    jlog(
        "k-over-s2",
        &format!("\"omega_d\":{omega_d},\"peaks\":{}", peaks.len()),
    );
}

#[test]
fn delay_margin_exists_and_its_loss_diverges() {
    // Positive effective phase margin, executed as DELAY margin: +100 ms
    // of extra reaction still converges; +400 ms diverges (the falsifier
    // proving the margin measurement is live).
    let mut margin = quiet_pilot(0);
    margin.gains.reaction_ticks += 12; // +100 ms
    let t1 = run_double_integrator(&margin, 2.0, 40.0);
    let tail_ok = t1[t1.len() - 240..]
        .iter()
        .all(|v| v.is_finite() && *v < 0.05);
    assert!(tail_ok, "+100 ms must remain inside the margin");
    let mut beyond = quiet_pilot(0);
    beyond.gains.reaction_ticks += 84; // +700 ms
    let t2 = run_double_integrator(&beyond, 2.0, 40.0);
    let diverged =
        t2.last().copied().unwrap_or(f64::INFINITY) > 1.0 || t2.iter().any(|v| !v.is_finite());
    assert!(diverged, "+700 ms must exceed the margin (live falsifier)");
    jlog(
        "delay-margin",
        "\"plus_100ms\":true,\"plus_700ms_diverges\":true",
    );
}

#[test]
fn single_integrator_settles_without_oscillation_class_failure() {
    // K/s: the easy element; settles fast and stays settled.
    let pilot = quiet_pilot(0);
    let perc = quiet_perception(3);
    let mut ps = perc.init().unwrap();
    let mut st = pilot.init().unwrap();
    let mut phi = 0.1f64;
    let mut trace = Vec::new();
    for _ in 0..(20.0 / DT) as usize {
        let cues = perc
            .step(&mut ps, pack_cues(0.0, 0.0, 0.0, phi, 0.0, 0.0))
            .unwrap();
        let cmd = pilot.step(&mut st, &cues, 0.0, 0.0, 0.0, 0.0).unwrap();
        phi += DT * 4.0 * cmd.warp_cmd_rad;
        trace.push(phi.abs());
    }
    assert!(
        trace[trace.len() - 120..].iter().all(|v| *v < 0.01),
        "K/s loop must settle"
    );
    jlog("k-over-s", "\"settled\":true");
}

#[test]
fn saturation_clamps_at_the_limits_with_receipt_flags() {
    let pilot = quiet_pilot(0);
    let perc = quiet_perception(3);
    let mut ps = perc.init().unwrap();
    let mut st = pilot.init().unwrap();
    // A huge roll error must saturate the warp channel EXACTLY at the
    // dossier limit and set the flag; the force channel likewise via a
    // huge lever error.
    let mut last = None;
    for _ in 0..120 {
        let cues = perc
            .step(&mut ps, pack_cues(0.0, 0.0, 0.0, 5.0, 0.0, 0.0))
            .unwrap();
        last = Some(pilot.step(&mut st, &cues, 3.0, 0.0, 0.0, 0.0).unwrap());
    }
    let cmd = last.unwrap();
    assert_eq!(
        cmd.warp_cmd_rad, -pilot.gains.warp_limit_rad,
        "warp clamped at the limit"
    );
    assert!(cmd.saturated[1], "warp saturation flag");
    assert_eq!(
        cmd.lever_force_n, -pilot.gains.force_limit_n,
        "force clamped at the limit"
    );
    assert!(cmd.saturated[0], "force saturation flag");
    // Tiny errors: no flags.
    let mut st2 = pilot.init().unwrap();
    let mut ps2 = quiet_perception(3).init().unwrap();
    let cues = perc
        .step(&mut ps2, pack_cues(0.0, 0.0, 0.0, 1e-4, 0.0, 0.0))
        .unwrap();
    let small = pilot.step(&mut st2, &cues, 0.0, 0.0, 0.0, 0.0).unwrap();
    assert!(!small.saturated[0] && !small.saturated[1]);
    jlog("saturation", "\"flags_and_exact_clamp\":true");
}

#[test]
fn family_caps_at_cap_and_cap_plus_one() {
    assert!(
        pilot_family_v1(FAMILY_SIZE - 1).is_ok(),
        "last member admits"
    );
    assert_eq!(
        pilot_family_v1(FAMILY_SIZE).unwrap_err().code,
        "pilot-member-invalid",
        "cap+1 refuses"
    );
    let mut bad = PilotWrightModel::new(0, 1).unwrap();
    bad.gains.reaction_ticks = fs_flyer::pilot::MAX_REACTION_TICKS + 1;
    assert_eq!(bad.admit().unwrap_err().code, "pilot-gains-invalid");
    bad.gains.reaction_ticks = fs_flyer::pilot::MAX_REACTION_TICKS;
    assert!(bad.admit().is_ok(), "reaction cap admits");
    jlog("family", "\"cap_and_cap_plus_one\":true");
}

#[test]
fn remnant_is_deterministic_and_seeded() {
    let run = |seed: u64| -> Vec<u64> {
        let pilot = PilotWrightModel::new(0, seed).unwrap();
        let perc = quiet_perception(3);
        let mut ps = perc.init().unwrap();
        let mut st = pilot.init().unwrap();
        let mut out = Vec::new();
        for _ in 0..240 {
            let cues = perc
                .step(&mut ps, pack_cues(0.01, 0.0, 0.0, 0.01, 0.0, 0.0))
                .unwrap();
            let cmd = pilot.step(&mut st, &cues, 0.0, 0.0, 0.0, 0.0).unwrap();
            out.push(cmd.lever_force_n.to_bits());
            out.push(cmd.warp_cmd_rad.to_bits());
        }
        out
    };
    assert_eq!(run(9), run(9), "bitwise repeat");
    assert_ne!(run(9), run(10), "seed liveness");
    jlog("remnant", "\"seeded_deterministic\":true");
}

#[test]
fn h02c_closed_loop_record_against_the_real_airframe() {
    // H-02c: compatibility RECORD at Estimated. Close the loop:
    // linearized 4-state airframe (E4.6a-iii) + FD control column +
    // canard mechanism (E4.6b-i) + perception (E4.6c-i) + pilot. The
    // receipt records bounded-ness, max |theta|, undulation period, and
    // the divergence-delay ratio vs open loop. Culick's account is
    // 'barely stabilizable' — the gate here is the RECEIPT plus the
    // weak liveness claim that the best member delays divergence vs
    // open loop.
    let d = wright_openloop_v1();
    let rho = 1.294;
    let trim = d.trim(rho, [13.0, 0.06, 0.1, 45.0]).unwrap();
    let rep = linearize(&d, &trim, rho).unwrap();
    let a = rep.a;
    // Control column via central FD on the canard setting.
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
    // Simulate: x = (u, w, q, theta) perturbations; dc from the
    // mechanism driven by pilot force; horizon 15 s.
    let horizon = (15.0 / DT) as usize;
    let escape = 0.35f64;
    // TRAINED-PILOT perception member (declared modeling choice for the
    // H-02c record): the standard perception delays model casual
    // sensing and would DOUBLE-COUNT latency in series with the pilot's
    // reaction delay — the crossover literature's effective pilot delay
    // (tau_e ~ 0.15-0.3 s) covers sensing+reaction combined. Total here:
    // 67 ms attitude sensing + 100 ms reaction + 80 ms neuromuscular.
    let mut fast = quiet_perception(5);
    fast.cues[0].delay_ticks = 8;
    fast.cues[1].delay_ticks = 5;
    fast.cues[2].delay_ticks = 4;
    fast.cues[0].filter_tau_s = 0.03;
    fast.cues[1].filter_tau_s = 0.03;
    let sim = |member: Option<u32>| -> (f64, f64, Option<f64>) {
        // returns (time_to_escape or horizon, max_theta, undulation period)
        let perc = fast;
        let mut ps = perc.init().unwrap();
        let pilot = member.map(quiet_pilot);
        let mut pst = pilot.as_ref().map(|p| p.init().unwrap());
        let mut mech = MechState {
            delta_rad: 0.0,
            rate_rad_s: 0.0,
        };
        let mut x = [0.0, 0.0, 0.0, 0.05f64];
        let mut max_theta = 0.0f64;
        let mut crossings: Vec<f64> = Vec::new();
        let mut prev_sign = 1.0f64;
        for step in 0..horizon {
            let theta = x[3];
            max_theta = max_theta.max(theta.abs());
            if theta.abs() > escape {
                return (step as f64 * DT, max_theta, undulation(&crossings));
            }
            if theta.signum() != prev_sign && theta.abs() > 0.005 {
                crossings.push(step as f64 * DT);
                prev_sign = theta.signum();
            }
            // Perceived cues from the physics state (w-dot proxy row 1).
            let wdot: f64 = (0..4).map(|j| a[1][j] * x[j]).sum::<f64>();
            let cues = perc
                .step(&mut ps, pack_cues(theta, x[2], wdot, 0.0, 0.0, 0.0))
                .unwrap();
            let dc = if let (Some(p), Some(st)) = (pilot.as_ref(), pst.as_mut()) {
                let cmd = p
                    .step(st, &cues, mech.delta_rad, mech.rate_rad_s, 0.0, 0.0)
                    .unwrap();
                mech = CANARD_MECH_V1
                    .step(mech, 0.0, cmd.lever_force_n, DT)
                    .unwrap()
                    .0;
                mech.delta_rad
            } else {
                0.0
            };
            for i in 0..4 {
                let mut dx: f64 = (0..4).map(|j| a[i][j] * x[j]).sum::<f64>();
                dx += b[i] * dc;
                x[i] += DT * dx;
            }
        }
        (15.0, max_theta, undulation(&crossings))
    };
    // Diagnostic probe: IDEAL static feedback dc = -k*theta - kq*q with
    // no delay, no mechanism — separates sign errors from phase loss.
    {
        let mut x = [0.0, 0.0, 0.0, 0.05f64];
        let mut maxth = 0.0f64;
        let mut esc = 15.0;
        for step in 0..horizon {
            let dc = -2.4 * x[3] - 1.8 * x[2];
            for i in 0..4 {
                let mut dx: f64 = (0..4).map(|j| a[i][j] * x[j]).sum::<f64>();
                dx += b[i] * dc;
                x[i] += DT * dx;
            }
            maxth = maxth.max(x[3].abs());
            if x[3].abs() > escape {
                esc = step as f64 * DT;
                break;
            }
        }
        jlog(
            "h02c-ideal-probe",
            &format!("\"t_escape_s\":{esc},\"max_theta\":{maxth}"),
        );
        // 'Barely stabilizable by a trained pilot' (Culick): with ZERO
        // latency the static gains hold the aircraft — the airframe is
        // stabilizable in principle; latency is what makes it hard.
        assert!(
            maxth < 0.1,
            "ideal (zero-latency) feedback must stabilize: max theta {maxth}"
        );
    }
    let open = sim(None);
    let mut best: Option<(u32, (f64, f64, Option<f64>))> = None;
    for member in 0..3u32 {
        let r = sim(Some(member));
        jlog(
            "h02c-member",
            &format!(
                "\"member\":{member},\"t_escape_s\":{},\"max_theta\":{},\"undulation_s\":{:?}",
                r.0, r.1, r.2
            ),
        );
        if best.as_ref().is_none_or(|(_, b)| r.0 > b.0) {
            best = Some((member, r));
        }
    }
    let (bm_id, br) = best.unwrap();
    // The OVERCONTROL TENDENCY (H-02c's licensed qualitative claim):
    // with realistic latency the piloted loop becomes a growing
    // pilot-involved oscillation that escapes FASTER than the open-loop
    // divergence — Orville's 'rise suddenly ... then dart for the
    // ground'; the Dec-17 flights ended within 12-59 s. Recorded, never
    // forced to look like stabilization.
    let overcontrol_present = br.0 < 15.0;
    // Sign liveness via the HOSTILE REVERSED-GAIN TWIN: a pilot pushing
    // the wrong way must be strictly worse than the best real member —
    // proving the real members' inputs are corrective in sign even
    // though latency defeats them.
    let reversed = {
        let mut m = quiet_pilot(bm_id);
        m.gains.k_theta = -m.gains.k_theta;
        m.gains.k_q = -m.gains.k_q;
        let perc = fast;
        let mut ps = perc.init().unwrap();
        let mut st = m.init().unwrap();
        let mut mech = MechState {
            delta_rad: 0.0,
            rate_rad_s: 0.0,
        };
        let mut x = [0.0, 0.0, 0.0, 0.05f64];
        let mut esc = 15.0;
        for step in 0..horizon {
            if x[3].abs() > escape {
                esc = step as f64 * DT;
                break;
            }
            let wdot: f64 = (0..4).map(|j| a[1][j] * x[j]).sum::<f64>();
            let cues = perc
                .step(&mut ps, pack_cues(x[3], x[2], wdot, 0.0, 0.0, 0.0))
                .unwrap();
            let cmd = m
                .step(&mut st, &cues, mech.delta_rad, mech.rate_rad_s, 0.0, 0.0)
                .unwrap();
            mech = CANARD_MECH_V1
                .step(mech, 0.0, cmd.lever_force_n, DT)
                .unwrap()
                .0;
            for i in 0..4 {
                let mut dx: f64 = (0..4).map(|j| a[i][j] * x[j]).sum::<f64>();
                dx += b[i] * mech.delta_rad;
                x[i] += DT * dx;
            }
        }
        esc
    };
    jlog(
        "h02c-record",
        &format!(
            "\"claim_level\":\"Estimated\",\"open_t_escape_s\":{},\"best_member\":{bm_id},\"best_t_escape_s\":{},\"best_undulation_s\":{:?},\"overcontrol_tendency\":{overcontrol_present},\"reversed_twin_t_escape_s\":{reversed},\"ideal_zero_latency_stabilizes\":true",
            open.0, br.0, br.2
        ),
    );
    assert!(
        overcontrol_present,
        "the H-02c record requires the endpoint"
    );
    assert!(
        reversed < br.0,
        "the reversed-gain twin must be strictly worse: {reversed} vs {}",
        br.0
    );
}

fn undulation(crossings: &[f64]) -> Option<f64> {
    if crossings.len() < 3 {
        return None;
    }
    // Mean full period from successive same-direction crossings.
    let mut periods = Vec::new();
    for i in 2..crossings.len() {
        periods.push(crossings[i] - crossings[i - 2]);
    }
    Some(periods.iter().sum::<f64>() / periods.len() as f64)
}

#[test]
fn golden_digest() {
    let pilot = PilotWrightModel::new(0, 42).unwrap();
    let perc = perception_v1(42);
    let mut ps = perc.init().unwrap();
    let mut st = pilot.init().unwrap();
    let mut payload = Vec::new();
    for tick in 0..240 {
        let t = tick as f64 * DT;
        let raw: [f64; N_CUES] = pack_cues(
            0.05 * (1.5 * t).sin(),
            0.07 * (1.5 * t).cos(),
            0.0,
            0.02 * t.sin(),
            0.02 * t.cos(),
            0.0,
        );
        let cues = perc.step(&mut ps, raw).unwrap();
        let cmd = pilot.step(&mut st, &cues, 0.0, 0.0, 0.0, 0.0).unwrap();
        payload.extend_from_slice(&cmd.lever_force_n.to_bits().to_le_bytes());
        payload.extend_from_slice(&cmd.warp_cmd_rad.to_bits().to_le_bytes());
    }
    let digest =
        fs_blake3::hash_domain("org.frankensim.fs-flyer.v02c1-golden.v1", &payload).to_hex();
    jlog("golden", &format!("\"digest\":\"{digest}\""));
    assert_eq!(
        digest, "3c222c45d2877d2fadf457d02cecd2b09b3235d23cb88869ed78b46bf317f1ba",
        "pilot golden moved — determinism regression or an intentional \
         model change requiring the golden-bump protocol"
    );
}
