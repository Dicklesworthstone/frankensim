//! E0.9c backward-playback battery (bead wf-root-guzez.1.9.3): loads the
//! committed archive fixture through the verifying loader, replays the
//! archived generation on the current kernel, and drives the hostile twins
//! (corruption, mirror divergence, malformed envelope, wrong-kernel replay).
//! JSONL receipts on stdout.

use fs_flyer_wasm::archive::{
    ArchiveTarget, HELLO_ENVELOPE_SCHEMA, parse_hello_envelope, replay_generation,
    verify_dual_publication, verify_target_bytes,
};
use std::fs;
use std::path::PathBuf;

fn fixture_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../data/wright-flyer/archive-fixture")
}

/// Strict targets-manifest parser: `path|size|blake3` per line, no tolerance.
fn load_targets() -> Vec<ArchiveTarget> {
    let text = fs::read_to_string(fixture_root().join("targets.v1.txt")).unwrap();
    text.trim_end_matches('\n')
        .split('\n')
        .map(|line| {
            let parts: Vec<&str> = line.split('|').collect();
            assert_eq!(parts.len(), 3, "malformed targets line {line:?}");
            ArchiveTarget {
                path: parts[0].to_string(),
                size_bytes: parts[1].parse().unwrap(),
                blake3_hex: parts[2].to_string(),
            }
        })
        .collect()
}

#[test]
fn backward_playback_of_generation_0000() {
    let targets = load_targets();
    assert_eq!(targets.len(), 1, "fixture archives exactly one generation");
    let target = &targets[0];
    let bytes = fs::read(fixture_root().join(&target.path)).unwrap();
    verify_target_bytes(target, &bytes).unwrap_or_else(|r| {
        panic!(
            "fixture failed content verification: {} — {}",
            r.code, r.message
        )
    });
    // Dual-publication read-back: the local store stands in for R2+mirror.
    verify_dual_publication(target, &bytes, &bytes).unwrap();
    let env = parse_hello_envelope(std::str::from_utf8(&bytes).unwrap()).unwrap();
    assert_eq!(env.generation, 0);
    let digest = replay_generation(&env).expect("archived generation must replay old-exact");
    println!(
        "{{\"suite\":\"wf-archive\",\"case\":\"backward-playback\",\"generation\":0,\
         \"digest\":\"{digest}\",\"verdict\":\"OLD-EXACT\"}}"
    );
}

#[test]
fn corruption_twin_refuses_before_parse() {
    let targets = load_targets();
    let target = &targets[0];
    let mut bytes = fs::read(fixture_root().join(&target.path)).unwrap();
    let mid = bytes.len() / 2;
    bytes[mid] ^= 0x01; // single flipped bit, size preserved
    let refusal = verify_target_bytes(target, &bytes).unwrap_err();
    assert_eq!(refusal.code, "archive-content-digest-mismatch");
    // Truncation is caught by the SIZE check, before any hashing.
    let truncated = &bytes[..bytes.len() - 1];
    assert_eq!(
        verify_target_bytes(target, truncated).unwrap_err().code,
        "archive-size-mismatch"
    );
    println!("{{\"suite\":\"wf-archive\",\"case\":\"corruption-twin\",\"verdict\":\"REFUSED\"}}");
}

#[test]
fn mirror_divergence_twin_refuses() {
    let targets = load_targets();
    let target = &targets[0];
    let primary = fs::read(fixture_root().join(&target.path)).unwrap();
    let mut mirror = primary.clone();
    let last = mirror.len() - 1;
    mirror[last] ^= 0x01;
    let refusal = verify_dual_publication(target, &primary, &mirror).unwrap_err();
    // The tampered mirror fails its own content check first (fail-closed).
    assert_eq!(refusal.code, "archive-content-digest-mismatch");
    println!("{{\"suite\":\"wf-archive\",\"case\":\"mirror-divergence\",\"verdict\":\"REFUSED\"}}");
}

#[test]
fn malformed_envelope_twins_refuse() {
    let good =
        fs::read_to_string(fixture_root().join("generation-0000/hello-envelope.v1.txt")).unwrap();
    // Reordered keys (same content, wrong canonical order).
    let mut lines: Vec<&str> = good.trim_end().split('\n').collect();
    lines.swap(2, 3);
    let reordered = lines.join("\n");
    assert_eq!(
        parse_hello_envelope(&reordered).unwrap_err().code,
        "archive-envelope-malformed"
    );
    // Extra trailing line.
    let extended = format!("{good}extra=1\n");
    assert_eq!(
        parse_hello_envelope(&extended).unwrap_err().code,
        "archive-envelope-malformed"
    );
    // Foreign schema id.
    let foreign = good.replace(HELLO_ENVELOPE_SCHEMA, "org.example.other.v1");
    assert_eq!(
        parse_hello_envelope(&foreign).unwrap_err().code,
        "archive-envelope-malformed"
    );
    println!(
        "{{\"suite\":\"wf-archive\",\"case\":\"malformed-envelope\",\"verdict\":\"REFUSED\"}}"
    );
}

#[test]
fn wrong_kernel_replay_twin_refuses() {
    // An archived digest the current kernel cannot reproduce must be a typed
    // refusal, not a pass — this is the old-exact contract's teeth.
    let good =
        fs::read_to_string(fixture_root().join("generation-0000/hello-envelope.v1.txt")).unwrap();
    let tampered = good.replace("steps=480", "steps=481");
    let env = parse_hello_envelope(&tampered).unwrap();
    let refusal = replay_generation(&env).unwrap_err();
    assert_eq!(refusal.code, "archive-replay-digest-mismatch");
    println!(
        "{{\"suite\":\"wf-archive\",\"case\":\"wrong-kernel-replay\",\"verdict\":\"REFUSED\"}}"
    );
}
