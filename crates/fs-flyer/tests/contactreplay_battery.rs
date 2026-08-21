//! E10.5 certified-contact replay pass (bead wf-root-guzez.11.9).
//! Recorded flights' motion is fed through fs-contact's certified
//! spacetime CCD against the terrain to certify that the real-time
//! model missed NO skid/ground penetration. This battery IS the
//! offline referee pass (standing infra, dev-deps only — fs-flyer's
//! production cone is untouched): it records a REAL lifecycle run,
//! replays every airborne tick-pair through certified_ccd, and
//! publishes the receipt as JSONL with a pinned digest.
//!
//! The FALSIFIER is the ct-005 class: a deliberately tunneling trace
//! crosses a terrain ridge INSIDE one tick-pair — the naive endpoint
//! check is EXECUTED and provably blind (both endpoint boxes are
//! disjoint from the ridge), and the certified pass CATCHES it.
//! Repro: cargo test -p fs-flyer --test contactreplay_battery

use asupersync::types::Budget;
use fs_contact::{CcdVerdict, SpacetimeBody, certified_ccd};
use fs_exec::{CancelGate, Cx, ExecMode, StreamKey};
use fs_flyer::simloop::{Phase, PilotMode, ScenarioSpec, SimLoop, TerminalEvent};
use fs_ga::Motor;
use fs_geom::{Aabb, Point3};
use fs_ivl::Interval;
use fs_motion::{CertifiedMotorTube, ScrewParams, screw_tube};

fn jlog(case: &str, payload: &str) {
    println!("{{\"suite\":\"fs-flyer-contactreplay\",\"case\":\"{case}\",{payload}}}");
}

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 0x1903,
                kernel_id: 105,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx)
    })
}

const DT: f64 = 1.0 / 120.0;
/// Sample cap for one replay.
const MAX_SAMPLES: usize = 4_096;
/// CCD time tolerance [s].
const TIME_TOL: f64 = 1.0e-4;
/// CCD subwindow budget per tick-pair.
const CCD_BUDGET: usize = 512;

#[derive(Clone, Copy, Debug)]
struct MotionSample {
    tick: u64,
    x_m: f64,
    h_m: f64,
}

/// The skid proxy box CENTERED at the recorded reference point.
fn skid_box(x: f64, h: f64) -> Aabb {
    Aabb::new(
        Point3::new(x - 1.5, -0.6, h - 0.12),
        Point3::new(x + 1.5, 0.6, h + 0.12),
    )
}

fn translation_tube(axis: [f64; 3], speed: f64, domain: Interval) -> CertifiedMotorTube {
    screw_tube(
        &ScrewParams {
            axis,
            center: [0.0, 0.0, 0.0],
            omega: 0.0,
            axial_velocity: speed,
            base_pose: Motor::identity(),
        },
        domain,
        4,
        8,
    )
    .expect("analytic translation tube")
}

#[derive(Clone, Debug, PartialEq)]
enum PairVerdict {
    Clear { min_gap: f64 },
    Possible { first_window_lo: f64 },
}

#[derive(Clone, Debug)]
struct ContactReplayReceiptV1 {
    scenario: &'static str,
    pairs: Vec<(u64, PairVerdict)>,
    realtime_contact_tick: Option<u64>,
    verdict: &'static str,
    first_missed_tick: Option<u64>,
    min_gap_worst: f64,
    receipt_digest: String,
}

/// Replay a recorded trace against the terrain through certified CCD.
/// `terrain` are STATIC boxes. Refusals are plain panics here (this
/// is the referee harness; the caps are still asserted both sides).
fn replay_certify(
    scenario: &'static str,
    samples: &[MotionSample],
    terrain: &[Aabb],
    realtime_contact_tick: Option<u64>,
) -> Result<ContactReplayReceiptV1, &'static str> {
    if samples.len() < 2 {
        return Err("replay-trace-too-short");
    }
    if samples.len() > MAX_SAMPLES {
        return Err("replay-trace-too-long");
    }
    if samples.windows(2).any(|w| w[1].tick != w[0].tick + 1) {
        return Err("replay-trace-nonmonotonic");
    }
    let mut pairs = Vec::new();
    let mut min_gap_worst = f64::INFINITY;
    let mut first_missed = None;
    with_cx(|cx| {
        for w in samples.windows(2) {
            let (s0, s1) = (w[0], w[1]);
            let t0 = s0.tick as f64 * DT;
            let t1 = s1.tick as f64 * DT;
            let window = Interval::new(t0, t1);
            let dx = s1.x_m - s0.x_m;
            let dh = s1.h_m - s0.h_m;
            let dist = (dx * dx + dh * dh).sqrt();
            let (axis, speed) = if dist > 1e-12 {
                ([dx / dist, 0.0, dh / dist], dist / DT)
            } else {
                ([1.0, 0.0, 0.0], 0.0)
            };
            let craft_tube = translation_tube(axis, speed, window);
            // Fold the start pose into the support (the ct-001 idiom);
            // the tube's translation carries it across the window. The
            // tube acts from its domain START, so offset by -t0*v.
            let support = {
                let b = skid_box(s0.x_m, s0.h_m);
                Aabb::new(
                    Point3::new(
                        b.min.x - axis[0] * speed * t0,
                        b.min.y,
                        b.min.z - axis[2] * speed * t0,
                    ),
                    Point3::new(
                        b.max.x - axis[0] * speed * t0,
                        b.max.y,
                        b.max.z - axis[2] * speed * t0,
                    ),
                )
            };
            let static_tube = translation_tube([1.0, 0.0, 0.0], 0.0, window);
            let craft = SpacetimeBody::new(support, &craft_tube).expect("craft body");
            let mut pair_verdict = PairVerdict::Clear {
                min_gap: f64::INFINITY,
            };
            for ground in terrain {
                let ground_body = SpacetimeBody::new(*ground, &static_tube).expect("terrain body");
                let report = certified_ccd(&craft, &ground_body, window, TIME_TOL, CCD_BUDGET, cx)
                    .expect("ccd within budget");
                match (report.verdict, &mut pair_verdict) {
                    (CcdVerdict::ClearWindow { min_gap }, PairVerdict::Clear { min_gap: mg }) => {
                        *mg = mg.min(min_gap);
                    }
                    (CcdVerdict::ClearWindow { .. }, PairVerdict::Possible { .. }) => {}
                    (CcdVerdict::PossibleContact { windows }, v) => {
                        let lo = windows.first().map_or(t0, Interval::lo);
                        if matches!(v, PairVerdict::Clear { .. }) {
                            *v = PairVerdict::Possible {
                                first_window_lo: lo,
                            };
                        }
                    }
                }
            }
            if let PairVerdict::Clear { min_gap } = pair_verdict {
                min_gap_worst = min_gap_worst.min(min_gap);
            } else if first_missed.is_none() && realtime_contact_tick.is_none_or(|ct| s0.tick < ct)
            {
                first_missed = Some(s0.tick);
            }
            pairs.push((s0.tick, pair_verdict));
        }
    });
    let verdict = if first_missed.is_none() {
        "CERTIFIED"
    } else {
        "MISSED_CONTACT"
    };
    let mut b = scenario.as_bytes().to_vec();
    for (tick, v) in &pairs {
        b.extend_from_slice(&tick.to_le_bytes());
        match v {
            PairVerdict::Clear { min_gap } => {
                b.push(1);
                b.extend_from_slice(&min_gap.to_bits().to_le_bytes());
            }
            PairVerdict::Possible { first_window_lo } => {
                b.push(0);
                b.extend_from_slice(&first_window_lo.to_bits().to_le_bytes());
            }
        }
    }
    let receipt_digest = fs_blake3::hash_domain("org.frankensim.wf.contact-replay.v1", &b).to_hex();
    Ok(ContactReplayReceiptV1 {
        scenario,
        pairs,
        realtime_contact_tick,
        verdict,
        first_missed_tick: first_missed,
        min_gap_worst,
        receipt_digest,
    })
}

/// Record a REAL lifecycle run's motion trace (airborne portion).
fn record_reference_flight() -> (Vec<MotionSample>, Option<u64>) {
    let spec = ScenarioSpec {
        seed: 1903,
        rho_kg_m3: 1.294,
        headwind_mps: 11.0,
        pilot_mode: PilotMode::Historical(0),
        assist: None,
        catapult: None,
        rail_length_m: 18.3,
        max_ticks: 1_200,
    };
    let mut sim = SimLoop::init(spec).expect("reference scenario inits");
    let mut airborne = Vec::new();
    let mut contact_tick = None;
    loop {
        match sim.step(None) {
            Err(_) => break,
            Ok(out) => match out.phase {
                Phase::Airborne => airborne.push(MotionSample {
                    tick: out.tick,
                    x_m: out.x_m,
                    h_m: out.h_m,
                }),
                Phase::Ended(TerminalEvent::GroundContact) => {
                    contact_tick = Some(out.tick);
                    break;
                }
                Phase::Ended(_) => break,
                Phase::OnRail => {}
            },
        }
    }
    (airborne, contact_tick)
}

/// The certified flat: sand from z = -5 up to 0, wide enough for the
/// whole flight.
fn flat_ground() -> Aabb {
    Aabb::new(
        Point3::new(-50.0, -30.0, -5.0),
        Point3::new(400.0, 30.0, 0.0),
    )
}

#[test]
fn reference_flight_replay_is_certified() {
    let (trace, contact_tick) = record_reference_flight();
    assert!(
        trace.len() >= 30,
        "non-vacuity: the reference flight must actually fly ({} airborne ticks)",
        trace.len()
    );
    let receipt =
        replay_certify("dec17-ref-seed1903", &trace, &[flat_ground()], contact_tick).unwrap();
    assert_eq!(
        receipt.verdict, "CERTIFIED",
        "first missed: {:?}",
        receipt.first_missed_tick
    );
    // Per-pair oracles: every airborne pair PROVEN clear with a
    // positive certified gap (never a totals-only claim).
    let mut clear = 0;
    for (tick, v) in &receipt.pairs {
        match v {
            PairVerdict::Clear { min_gap } => {
                assert!(*min_gap > 0.0, "tick {tick}: certified gap {min_gap}");
                clear += 1;
            }
            PairVerdict::Possible { .. } => {
                // Only permissible at/after the real-time contact tick.
                assert!(
                    contact_tick.is_some_and(|ct| *tick >= ct),
                    "tick {tick}: possible contact before the model's own contact"
                );
            }
        }
    }
    assert!(clear >= 30);
    // Determinism: bit-identical receipt twice.
    let again =
        replay_certify("dec17-ref-seed1903", &trace, &[flat_ground()], contact_tick).unwrap();
    assert_eq!(
        again.receipt_digest, receipt.receipt_digest,
        "bit-identical twice"
    );
    jlog(
        "reference",
        &format!(
            "\"airborne_pairs\":{},\"min_gap_worst\":{},\"verdict\":\"{}\",\"receipt_digest\":\"{}\"",
            receipt.pairs.len(),
            receipt.min_gap_worst,
            receipt.verdict,
            receipt.receipt_digest
        ),
    );
}

#[test]
fn tunneling_falsifier_is_caught_where_endpoint_sampling_is_blind() {
    // A ridge the real-time model never saw: x in [20.9, 21.1], up to
    // z = 1.5. The doctored trace crosses it INSIDE one tick-pair at
    // 480 m/s with both endpoints clear.
    let ridge = Aabb::new(Point3::new(20.9, -30.0, 0.0), Point3::new(21.1, 30.0, 1.5));
    let trace = vec![
        MotionSample {
            tick: 0,
            x_m: 19.0,
            h_m: 1.0,
        },
        MotionSample {
            tick: 1,
            x_m: 23.0,
            h_m: 1.0,
        },
    ];
    // The EXECUTED naive endpoint check: both endpoint skid boxes are
    // disjoint from the ridge in x — endpoint sampling is blind.
    for s in &trace {
        let b = skid_box(s.x_m, s.h_m);
        let disjoint_x = b.max.x < ridge.min.x || ridge.max.x < b.min.x;
        assert!(
            disjoint_x,
            "the trap requires clear endpoints (x {})",
            s.x_m
        );
    }
    // The certified pass CATCHES the crossing.
    let receipt =
        replay_certify("tunneling-falsifier", &trace, &[flat_ground(), ridge], None).unwrap();
    assert_eq!(receipt.verdict, "MISSED_CONTACT");
    assert_eq!(receipt.first_missed_tick, Some(0));
    match &receipt.pairs[0].1 {
        PairVerdict::Possible { first_window_lo } => {
            // The true crossing starts near x=20.9-1.5(half-length) →
            // t* ≈ (19.4→20.9 gap 0.4 m at 480 m/s) ≈ 0.00083 s; the
            // certified window must contain crossing onset in [t0, t1].
            assert!(
                *first_window_lo >= 0.0 && *first_window_lo <= DT,
                "window lo {first_window_lo}"
            );
        }
        other => panic!("the falsifier must be caught, got {other:?}"),
    }
    jlog(
        "tunneling",
        &format!(
            "\"verdict\":\"{}\",\"first_missed_tick\":0",
            receipt.verdict
        ),
    );
}

#[test]
fn replay_caps_and_refusals() {
    let mk = |n: usize| -> Vec<MotionSample> {
        (0..n)
            .map(|i| MotionSample {
                tick: i as u64,
                x_m: i as f64 * 0.1,
                h_m: 3.0,
            })
            .collect()
    };
    // AT the sample cap admits; one more refuses.
    assert!(replay_certify("caps", &mk(MAX_SAMPLES), &[flat_ground()], None).is_ok());
    assert_eq!(
        replay_certify("caps", &mk(MAX_SAMPLES + 1), &[flat_ground()], None).unwrap_err(),
        "replay-trace-too-long"
    );
    assert_eq!(
        replay_certify("caps", &mk(1), &[flat_ground()], None).unwrap_err(),
        "replay-trace-too-short"
    );
    let mut gap = mk(4);
    gap[2].tick = 5;
    assert_eq!(
        replay_certify("caps", &gap, &[flat_ground()], None).unwrap_err(),
        "replay-trace-nonmonotonic"
    );
    jlog("caps", &format!("\"max_samples\":{MAX_SAMPLES}"));
}
