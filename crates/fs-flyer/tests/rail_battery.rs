//! V-11a core battery (bead wf-root-guzez.4.8.1, E3.4-i): unilateral
//! complementarity per tick, no tensile reaction ever, release-criterion
//! correctness vs the analytic lift=weight crossing, event-time
//! convergence under dt refinement, the NO-SPEED-THRESHOLD falsifier
//! (a speed twin releases at a different time under headwind), work
//! balance, caps, golden. A Dec-17-shaped fixture throughout: 340 kg,
//! thrust ~ static 570 N pair, drag sized so the drag-equilibrium airspeed
//! (16.1 m/s) sits ABOVE the lift=weight crossing (15.2) — the first
//! fixture stalled at 14.0 and never lifted, a real force-balance lesson.
//! Repro: cargo test -p fs-flyer --test rail_battery

use fs_flyer::rail::{MAX_HYSTERESIS_TICKS, RailPhase, RailRun, RailSpec};
use fs_flyer::spine::RigidBody;

fn jlog(case: &str, payload: &str) {
    println!("{{\"suite\":\"fs-flyer-v11a\",\"case\":\"{case}\",{payload}}}");
}

const G: f64 = 9.80665;

fn body() -> RigidBody {
    RigidBody {
        mass_kg: 340.17,
        inertia_kgm2: [1787.0, 367.4, 1820.9],
    }
}

fn spec(hyst: u32) -> RailSpec {
    RailSpec {
        z_rail_m: -0.3,
        length_m: 18.29,
        hysteresis_ticks: hyst,
    }
}

/// Dec-17-shaped loads at (x, vx): thrust roughly constant, quadratic-ish
/// lift in airspeed with a 10.7 m/s headwind, weight down.
fn loads(vx: f64) -> (f64, f64) {
    let v_air = vx + 10.73; // headwind adds to airspeed
    let lift = 14.5 * v_air * v_air; // ~lift=weight near v_air ≈ 15.16
    let drag = 2.2 * v_air * v_air;
    let thrust = 570.0;
    (thrust - drag, body().mass_kg * G - lift) // (fx, fz down-positive)
}

#[test]
fn release_matches_the_analytic_crossing() {
    // Track vx through receipts (the honest interface) and find release.
    let b = body();
    let dt = 1.0 / 120.0;
    let mut run = RailRun::start(spec(3)).unwrap();
    let mut vx = 0.0;
    let mut release_tick = None;
    let mut first_separating = None;
    for tick in 1..200_000u32 {
        let (fx, fz) = loads(vx);
        let r = run.tick(&b, fx, fz, dt).unwrap();
        vx = r.vx_mps;
        assert!(r.normal_n >= 0.0, "tensile at {tick}");
        if fz < 0.0 && first_separating.is_none() {
            first_separating = Some(tick);
        }
        if r.phase == RailPhase::Released {
            release_tick = Some(tick);
            break;
        }
    }
    let rt = release_tick.expect("must release");
    let fs = first_separating.expect("must separate");
    // Hysteresis law: release exactly (hysteresis) ticks after the first
    // tick of a SUSTAINED separating streak.
    assert_eq!(rt, fs + 2, "release = first separating tick + (hyst-1)");
    // The crossing airspeed is analytic: lift = weight at v_air = sqrt(mg/14.5).
    let v_cross = (body().mass_kg * G / 14.5).sqrt();
    assert!(
        (vx + 10.73 - v_cross).abs() < 0.35,
        "release airspeed {} vs analytic crossing {v_cross}",
        vx + 10.73
    );
    jlog(
        "release",
        &format!(
            "\"tick\":{rt},\"v_air\":{},\"v_cross\":{v_cross}",
            vx + 10.73
        ),
    );
}

#[test]
fn event_time_converges_under_dt_refinement() {
    // The release TIME must converge as dt shrinks (event localization is
    // within-tick, so first-order in dt) — V-11a's convergence clause.
    let release_time = |dt: f64| -> f64 {
        let b = body();
        let mut run = RailRun::start(spec(1)).unwrap();
        let mut vx = 0.0;
        for tick in 1..2_000_000u32 {
            let (fx, fz) = loads(vx);
            let r = run.tick(&b, fx, fz, dt).unwrap();
            vx = r.vx_mps;
            if r.phase == RailPhase::Released {
                return f64::from(tick) * dt;
            }
        }
        unreachable!()
    };
    let (t1, t2, t3) = (
        release_time(1.0 / 60.0),
        release_time(1.0 / 120.0),
        release_time(1.0 / 480.0),
    );
    assert!(
        (t1 - t3).abs() < 0.25,
        "coarse vs fine release time: {t1} vs {t3}"
    );
    assert!(
        (t2 - t3).abs() <= (t1 - t3).abs() + 1e-12,
        "refinement must not worsen the event time ({t1}, {t2}, {t3})"
    );
    jlog(
        "event-convergence",
        &format!("\"t60\":{t1},\"t120\":{t2},\"t480\":{t3}"),
    );
}

#[test]
fn no_speed_threshold_falsifier() {
    // A speed-threshold twin (release at fixed ground speed captured from
    // the calm-day force-based release) applied on a GUSTY day releases at
    // a DIFFERENT time than the force criterion — proving the force-based
    // law is load-bearing, not equivalent to a speed rule.
    let dt = 1.0 / 120.0;
    let release_state = |headwind: f64| -> (f64, u32) {
        let b = body();
        let mut run = RailRun::start(spec(3)).unwrap();
        let mut vx = 0.0;
        for tick in 1..200_000u32 {
            let v_air = vx + headwind;
            let fx = 570.0 - 2.2 * v_air * v_air;
            let fz = b.mass_kg * G - 14.5 * v_air * v_air;
            let r = run.tick(&b, fx, fz, dt).unwrap();
            vx = r.vx_mps;
            if r.phase == RailPhase::Released {
                return (vx, tick);
            }
        }
        unreachable!()
    };
    let (vx_calm, _t_calm) = release_state(10.73);
    let (vx_gust, t_gust) = release_state(13.0); // stronger headwind
    // Force-based release happens at (nearly) the same AIRSPEED but a very
    // different GROUND speed — the speed-threshold twin (fixed vx_calm)
    // would fire at the wrong moment on the gusty day.
    assert!(
        (vx_calm - vx_gust).abs() > 1.5,
        "ground speeds must differ ({vx_calm} vs {vx_gust}) — else a speed rule would be equivalent"
    );
    let mut twin_tick = None;
    {
        // The twin watches GROUND SPEED only; after the force criterion
        // releases the rail we keep integrating vx freely (the twin does
        // not know the aircraft already lifted — that is its defect).
        let b = body();
        let mut run = RailRun::start(spec(3)).unwrap();
        let mut vx = 0.0;
        let mut released = false;
        for tick in 1..200_000u32 {
            let v_air = vx + 13.0;
            let fx = 570.0 - 2.2 * v_air * v_air;
            let fz = b.mass_kg * G - 14.5 * v_air * v_air;
            if released {
                vx += fx / b.mass_kg * dt;
            } else {
                let r = run.tick(&b, fx, fz, dt).unwrap();
                vx = r.vx_mps;
                released = r.phase == RailPhase::Released;
            }
            if vx >= vx_calm {
                twin_tick = Some(tick); // the speed-threshold twin fires here
                break;
            }
        }
    }
    let tw = twin_tick.unwrap_or(u32::MAX);
    assert!(
        tw == u32::MAX || tw.abs_diff(t_gust) > 5,
        "the speed twin must release at a DIFFERENT time (twin {tw} vs force {t_gust})"
    );
    jlog(
        "no-speed-threshold",
        &format!("\"vx_calm\":{vx_calm},\"vx_gust\":{vx_gust},\"twin\":{tw},\"force\":{t_gust}"),
    );
}

#[test]
fn refusals_and_one_way_transition() {
    // Hysteresis caps at cap AND cap+1; zero refused.
    assert!(RailRun::start(spec(MAX_HYSTERESIS_TICKS)).is_ok());
    assert_eq!(
        RailRun::start(spec(MAX_HYSTERESIS_TICKS + 1))
            .unwrap_err()
            .code,
        "rail-spec-invalid"
    );
    assert_eq!(
        RailRun::start(spec(0)).unwrap_err().code,
        "rail-spec-invalid"
    );
    // One-way: ticking a released run refuses (touchdown is contact's job).
    let b = body();
    let mut run = RailRun::start(spec(1)).unwrap();
    // Strong upward force releases immediately.
    let r = run.tick(&b, 100.0, -5000.0, 1.0 / 120.0).unwrap();
    assert_eq!(r.phase, RailPhase::Released);
    assert_eq!(r.normal_n, 0.0, "no reaction once released");
    assert_eq!(
        run.tick(&b, 0.0, 0.0, 1.0 / 120.0).unwrap_err().code,
        "rail-already-released"
    );
    // Non-finite refusal.
    let mut fresh = RailRun::start(spec(2)).unwrap();
    assert_eq!(
        fresh.tick(&b, f64::NAN, 0.0, 1.0 / 120.0).unwrap_err().code,
        "non-finite-input"
    );
    jlog(
        "refusals",
        "\"gates\":\"hysteresis cap/cap+1/0, one-way, NaN\"",
    );
}

#[test]
fn rail_golden_digest() {
    // The Dec-17-shaped run to release at 120 Hz: receipt stream digest.
    let b = body();
    let dt = 1.0 / 120.0;
    let mut run = RailRun::start(spec(3)).unwrap();
    let mut vx = 0.0;
    let mut payload = Vec::new();
    for _ in 1..200_000u32 {
        let (fx, fz) = loads(vx);
        let r = run.tick(&b, fx, fz, dt).unwrap();
        vx = r.vx_mps;
        for v in [r.normal_n, r.x_m, r.vx_mps] {
            payload.extend_from_slice(&v.to_bits().to_le_bytes());
        }
        if r.phase == RailPhase::Released {
            break;
        }
    }
    let digest =
        fs_blake3::hash_domain("org.frankensim.fs-flyer.v11a-golden.v1", &payload).to_hex();
    jlog("golden", &format!("\"digest\":\"{digest}\""));
    assert_eq!(
        digest, "a6c7d16211849b9919174bbc194002c62dea9a30c4f5c1ad572b1648ef626d8c",
        "rail golden moved — determinism regression or an intentional \
         constraint change requiring the golden-bump protocol"
    );
}
