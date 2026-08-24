//! E6.1-i battery (bead wf-root-guzez.7.1.1): CheckpointStateV2
//! bit-identity — checkpoint mid-run, restore into a FRESH SimLoop,
//! march to terminal, and the chained digest EQUALS the uninterrupted
//! run's (rail-phase and airborne-phase checkpoints; FixedControls AND
//! Historical mode with its remnant/delay rings); the tampered-byte
//! hostile twin refuses; terminal checkpoints refuse; size caps at cap
//! AND cap+1.
//! Repro: cargo test -p fs-flyer --test checkpoint_battery

use fs_flyer::simloop::{Phase, PilotMode, SimLoop, dec17_scenario};

fn jlog(case: &str, payload: &str) {
    println!("{{\"suite\":\"fs-flyer-e61\",\"case\":\"{case}\",{payload}}}");
}

/// Run to terminal, checkpointing at `ckpt_tick`; return
/// (uninterrupted digest, checkpoint bytes, terminal tick).
fn run_with_checkpoint(mode: PilotMode, ckpt_tick: u64) -> (String, Vec<u8>, u64) {
    let mut sim = SimLoop::init(dec17_scenario(1903, mode)).unwrap();
    let mut ckpt = None;
    let end = loop {
        let out = sim.step(None).unwrap();
        if out.tick == ckpt_tick {
            ckpt = Some(sim.save_checkpoint().unwrap());
        }
        if let Phase::Ended(_) = out.phase {
            break out.tick;
        }
    };
    (
        sim.digest_hex(),
        ckpt.expect("checkpoint tick reached"),
        end,
    )
}

fn march_restored(mode: PilotMode, bytes: &[u8]) -> String {
    let mut sim = SimLoop::restore_checkpoint(dec17_scenario(1903, mode), bytes).unwrap();
    loop {
        let out = sim.step(None).unwrap();
        if let Phase::Ended(_) = out.phase {
            break;
        }
    }
    sim.digest_hex()
}

#[test]
fn fixed_mode_bit_identity_from_rail_and_airborne_checkpoints() {
    // Rail-phase checkpoint (tick 300) — the fixed run lifts at ~626.
    let (digest, ckpt_rail, _) = run_with_checkpoint(PilotMode::FixedControls, 300);
    assert_eq!(
        march_restored(PilotMode::FixedControls, &ckpt_rail),
        digest,
        "rail checkpoint bit-identity"
    );
    // Airborne-phase checkpoint.
    let (digest2, ckpt_air, end) = run_with_checkpoint(PilotMode::FixedControls, 800);
    assert_eq!(digest2, digest, "same scenario, same digest");
    assert_eq!(
        march_restored(PilotMode::FixedControls, &ckpt_air),
        digest,
        "airborne checkpoint bit-identity"
    );
    jlog(
        "fixed",
        &format!("\"end_tick\":{end},\"ckpt_bytes\":{}", ckpt_air.len()),
    );
}

#[test]
fn historical_mode_bit_identity_with_delay_rings() {
    // Historical(3): perception + pilot delay rings and the remnant
    // stream all live inside the checkpoint.
    let (digest, ckpt, end) = run_with_checkpoint(PilotMode::Historical(3), 900);
    assert_eq!(
        march_restored(PilotMode::Historical(3), &ckpt),
        digest,
        "historical checkpoint bit-identity"
    );
    jlog(
        "historical",
        &format!("\"end_tick\":{end},\"ckpt_bytes\":{}", ckpt.len()),
    );
}

#[test]
fn tamper_terminal_and_caps_refuse() {
    // SimLoop has no Debug; extract the code via match.
    let code_of = |r: Result<SimLoop, fs_flyer::Refusal>| -> &'static str {
        match r {
            Ok(_) => "OK",
            Err(e) => e.code,
        }
    };
    let mut sim = SimLoop::init(dec17_scenario(7, PilotMode::FixedControls)).unwrap();
    for _ in 0..10 {
        sim.step(None).unwrap();
    }
    let good = sim.save_checkpoint().unwrap();
    assert_eq!(u32::from_le_bytes(good[..4].try_into().unwrap()), 2);
    // Tamper twin: flip one payload byte — the embedded digest refuses.
    let mut bad = good.clone();
    bad[20] ^= 0x01;
    assert_eq!(
        code_of(SimLoop::restore_checkpoint(
            dec17_scenario(7, PilotMode::FixedControls),
            &bad
        )),
        "checkpoint-tampered"
    );
    // Truncation is malformed-or-tampered, never a silent accept.
    let code = code_of(SimLoop::restore_checkpoint(
        dec17_scenario(7, PilotMode::FixedControls),
        &good[..good.len() - 40],
    ));
    assert!(code == "checkpoint-tampered" || code == "checkpoint-malformed");
    // Pilot-presence mismatch: restore a fixed-mode checkpoint into a
    // historical scenario — refused, never silently rewired.
    assert!(
        SimLoop::restore_checkpoint(dec17_scenario(7, PilotMode::Historical(3)), &good).is_err()
    );
    // Size caps at cap AND cap+1.
    let cap = 1 << 20;
    let huge_ok = vec![0u8; cap]; // wrong digest but inside the cap
    assert_ne!(
        code_of(SimLoop::restore_checkpoint(
            dec17_scenario(7, PilotMode::FixedControls),
            &huge_ok
        )),
        "checkpoint-too-large"
    );
    let huge = vec![0u8; cap + 1];
    assert_eq!(
        code_of(SimLoop::restore_checkpoint(
            dec17_scenario(7, PilotMode::FixedControls),
            &huge
        )),
        "checkpoint-too-large"
    );
    // A terminal run refuses to checkpoint.
    let mut sim = SimLoop::init({
        let mut s = dec17_scenario(7, PilotMode::FixedControls);
        s.max_ticks = 3;
        s
    })
    .unwrap();
    for _ in 0..3 {
        sim.step(None).unwrap();
    }
    assert_eq!(
        sim.save_checkpoint().unwrap_err().code,
        "checkpoint-after-terminal"
    );
    jlog(
        "refusals",
        "\"tamper_truncate_mismatch_caps_terminal\":true",
    );
}
