//! E3.5 checkpoint battery (bead wf-root-guzez.4.10). Runs the digest
//! infrastructure ON THE E3.2 SPINE (DONE-WHEN clause) and executes the
//! localization falsifiers: a load-model perturbation localizes to
//! generalized-loads (the causal root) even though integrator-state also
//! diverges the same tick; a state-only perturbation localizes to
//! integrator-state; an absent channel is data, not a wildcard.
//! Repro: cargo test -p fs-flyer --test checkpoint_battery

use fs_flyer::checkpoint::{
    CheckpointBuilder, TickCheckpoint, f64_payload, first_divergence, subsystem_digest,
    verify_streams,
};
use fs_flyer::spine::{Loads, RigidBody, SixDofState, step};

fn jlog(case: &str, payload: &str) {
    println!("{{\"suite\":\"fs-flyer-checkpoint\",\"case\":\"{case}\",{payload}}}");
}

const DT: f64 = 1.0 / 120.0;

fn body() -> RigidBody {
    RigidBody {
        mass_kg: 340.17,
        inertia_kgm2: [1787.0, 367.4, 1820.9],
    }
}

fn rest() -> SixDofState {
    SixDofState {
        pos_m: [0.0; 3],
        vel_mps: [0.0; 3],
        quat: [1.0, 0.0, 0.0, 0.0],
        omega_body: [0.01, -0.02, 0.03],
    }
}

/// Run the spine with structured checkpoints: per tick, digest the
/// generalized loads (as computed at tick start) and the integrator state
/// (after the step). `perturb_loads`/`perturb_state` inject a small relative
/// fault at `fault_tick` into the respective subsystem.
fn run_checkpointed(
    ticks: u32,
    fault_tick: Option<u32>,
    perturb_loads: bool,
    perturb_state: bool,
) -> Vec<TickCheckpoint> {
    let b = body();
    let mut s = rest();
    let mut out = Vec::with_capacity(ticks as usize);
    for tick in 0..ticks {
        let t = f64::from(tick) * DT;
        // The load model (deterministic stand-in for aero+propulsion).
        let mut loads_now = Loads {
            force_n: [5.0 * (2.0 * t).sin(), 1.0, 3336.0 - 30.0 * s.vel_mps[2]],
            moment_nm: [2.0, 8.0 * t.cos(), -s.omega_body[2]],
        };
        if perturb_loads && fault_tick == Some(tick) {
            // A small RELATIVE fault (1e-9): large enough to survive the
            // integrator's rounding, small enough that only digests see it.
            loads_now.moment_nm[1] *= 1.0 + 1.0e-9;
        }
        // Step with the (possibly perturbed) load model held for this tick.
        s = step(&b, &s, t, DT, |_, _| loads_now).unwrap();
        if perturb_state && fault_tick == Some(tick) {
            // Perturb a component the load model READS (force_n[2] uses
            // vel_mps[2]) so the cross-tick causal chain is observable.
            s.vel_mps[2] += 1.0e-9;
        }
        let mut cp = CheckpointBuilder::new(tick);
        cp.record(
            "generalized-loads",
            &f64_payload(&[
                loads_now.force_n[0],
                loads_now.force_n[1],
                loads_now.force_n[2],
                loads_now.moment_nm[0],
                loads_now.moment_nm[1],
                loads_now.moment_nm[2],
            ]),
        )
        .unwrap();
        let state_payload: Vec<f64> = s
            .pos_m
            .iter()
            .chain(s.vel_mps.iter())
            .chain(s.quat.iter())
            .chain(s.omega_body.iter())
            .copied()
            .collect();
        cp.record("integrator-state", &f64_payload(&state_payload))
            .unwrap();
        out.push(cp.finish());
    }
    out
}

#[test]
fn clean_runs_verify_and_absence_is_data() {
    let a = run_checkpointed(240, None, false, false);
    let b = run_checkpointed(240, None, false, false);
    verify_streams(&a, &b).expect("identical runs must verify");
    // Absent channels (atmosphere etc. not yet landed) are None in BOTH —
    // and a run that suddenly PRODUCES one must diverge.
    let mut with_atmo = a.clone();
    let extra = subsystem_digest(0, "atmosphere", &f64_payload(&[1.294]));
    with_atmo[0].digests[0] = Some(extra);
    let d = first_divergence(&a, &with_atmo).expect("absence vs presence must diverge");
    assert_eq!(d.subsystem, "atmosphere");
    assert_eq!(d.expected, "absent");
    assert_eq!(d.tick, 0);
    jlog(
        "clean+absence",
        "\"identical\":true,\"absence_is_data\":true",
    );
}

#[test]
fn loads_fault_localizes_to_the_causal_root() {
    // DONE-WHEN falsifier: perturb the LOAD MODEL (1e-9 relative) at tick 100.
    // Both generalized-loads and integrator-state diverge from tick 100 on;
    // localization must name generalized-loads (the causal root), the exact
    // tick, and both digests — a whole-run hash could do none of this.
    let clean = run_checkpointed(240, None, false, false);
    let faulty = run_checkpointed(240, Some(100), true, false);
    let refusal = verify_streams(&clean, &faulty).unwrap_err();
    assert_eq!(refusal.code, "checkpoint-diverged");
    let d = first_divergence(&clean, &faulty).unwrap();
    assert_eq!(d.tick, 100, "must localize to the fault tick");
    assert_eq!(
        d.subsystem, "generalized-loads",
        "must name the CAUSAL ROOT"
    );
    assert_ne!(d.expected, d.observed);
    // Confirm the downstream consequence is real (integrator also moved at
    // tick 100) so the causal-order discrimination is non-vacuous.
    let idx_state = 5;
    assert_ne!(
        clean[100].digests[idx_state], faulty[100].digests[idx_state],
        "integrator state must also diverge at the fault tick (consequence)"
    );
    jlog(
        "loads-fault",
        &format!("\"tick\":{},\"subsystem\":\"{}\"", d.tick, d.subsystem),
    );
}

#[test]
fn state_fault_localizes_to_the_integrator() {
    // Perturb the INTEGRATOR STATE directly (post-loads) at tick 100: the
    // loads channel matches at tick 100 (it was computed pre-fault), so
    // localization must name integrator-state — proving the two faults are
    // DISCRIMINATED, which is the entire point of structured checkpoints.
    let clean = run_checkpointed(240, None, false, false);
    let faulty = run_checkpointed(240, Some(100), false, true);
    let d = first_divergence(&clean, &faulty).unwrap();
    assert_eq!(d.tick, 100);
    assert_eq!(d.subsystem, "integrator-state");
    // At tick 101 the loads DO diverge (state feeds the next tick's loads)
    // — the causal order across ticks is visible in the stream.
    assert_ne!(clean[101].digests[4], faulty[101].digests[4]);
    jlog(
        "state-fault",
        &format!("\"tick\":{},\"subsystem\":\"{}\"", d.tick, d.subsystem),
    );
}

#[test]
fn builder_and_digest_laws() {
    // Unknown subsystem refuses with the registry named.
    let mut cp = CheckpointBuilder::new(7);
    let refusal = cp.record("warp-drive", b"x").unwrap_err();
    assert_eq!(refusal.code, "subsystem-unknown");
    assert!(refusal.ranked_repairs[0].contains("atmosphere"));
    // Duplicate record refuses.
    cp.record("propulsion", b"a").unwrap();
    assert_eq!(
        cp.record("propulsion", b"a").unwrap_err().code,
        "subsystem-duplicate"
    );
    // The digest binds tick, channel, AND bytes (same bytes elsewhere differ).
    let base = subsystem_digest(3, "propulsion", b"abc");
    assert_ne!(
        base,
        subsystem_digest(4, "propulsion", b"abc"),
        "tick binds"
    );
    assert_ne!(
        base,
        subsystem_digest(3, "circulation", b"abc"),
        "channel binds"
    );
    assert_ne!(
        base,
        subsystem_digest(3, "propulsion", b"abd"),
        "bytes bind"
    );
    // Length-mismatch refusal.
    let a = run_checkpointed(10, None, false, false);
    let b = run_checkpointed(11, None, false, false);
    assert_eq!(
        verify_streams(&a, &b).unwrap_err().code,
        "checkpoint-stream-length-mismatch"
    );
    jlog(
        "laws",
        "\"binds\":\"tick+channel+bytes\",\"gates\":\"unknown, duplicate, length\"",
    );
}
