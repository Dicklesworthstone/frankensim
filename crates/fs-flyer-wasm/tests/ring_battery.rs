//! V-18 protocol battery (bead wf-root-guzez.6.1, E5.0): lease state
//! machine, ABI-version refusal, torn-read falsifier (deterministic
//! interleaving), restart/slot-reuse ABA twins, drop counters, slot
//! caps at BOTH edges and one past each, noninterference (reads never
//! mutate published state — bitwise), golden.
//! Repro: cargo test -p fs-flyer-wasm --test ring_battery

use fs_flyer_wasm::ring::{
    ABI_VERSION, Lease, MAX_SLOTS, MIN_SLOTS, RingCounters, SnapshotRing, payload_layout_hash,
};

fn jlog(case: &str, payload: &str) {
    println!("{{\"suite\":\"fs-flyer-wasm-v18\",\"case\":\"{case}\",{payload}}}");
}

fn payload(tick: u64) -> Vec<f64> {
    (0..16)
        .map(|i| tick as f64 + f64::from(i) * 0.001)
        .collect()
}

fn ring() -> SnapshotRing {
    SnapshotRing::new(3, 16, 0xA1C0FFEE, 0x1903).unwrap()
}

#[test]
fn lease_lifecycle_and_header_words() {
    let mut r = ring();
    r.publish(1, &payload(1)).unwrap();
    let lease = r.acquire(ABI_VERSION).unwrap();
    let mut buf = vec![0.0; 16];
    let h = r.read(&lease, &mut buf).unwrap();
    // Every ABI word present and bound to the run.
    assert_eq!(h.abi_version, ABI_VERSION);
    assert_eq!(h.payload_layout_hash, payload_layout_hash());
    assert_eq!(h.payload_len, 16);
    assert_eq!(h.run_epoch, 1);
    assert_eq!(h.run_anchor_digest_prefix, 0xA1C0FFEE);
    assert_eq!(h.model_id_prefix, 0x1903);
    assert_eq!(h.tick, 1);
    assert_eq!(h.sequence % 2, 0, "leased slot must be stable (even)");
    assert_eq!(buf, payload(1));
    r.release(&lease);
    jlog(
        "lifecycle",
        &format!("\"layout\":{}", h.payload_layout_hash),
    );
}

#[test]
fn abi_version_mismatch_refuses() {
    let mut r = ring();
    r.publish(1, &payload(1)).unwrap();
    assert_eq!(
        r.acquire(ABI_VERSION + 1).unwrap_err().code,
        "ring-abi-mismatch"
    );
    assert_eq!(r.acquire(0).unwrap_err().code, "ring-abi-mismatch");
    jlog("abi", "\"mismatch_refused\":true");
}

#[test]
fn torn_read_falsifier() {
    // Deterministic interleaving: reader leases slot, the writer path
    // rewrites it (sequence advances by 2), the held lease MUST fail
    // revalidation with the typed torn code — and the counter records it.
    let mut r = ring();
    r.publish(1, &payload(1)).unwrap();
    let lease = r.acquire(ABI_VERSION).unwrap();
    r.force_slot_rewrite(lease.slot, 99);
    let mut buf = vec![0.0; 16];
    let err = r.read(&lease, &mut buf).unwrap_err();
    assert_eq!(err.code, "ring-lease-torn");
    assert_eq!(r.counters().torn_reads, 1);
    // Re-acquire heals: the new lease reads cleanly.
    let lease2 = r.acquire(ABI_VERSION).unwrap();
    assert!(r.read(&lease2, &mut buf).is_ok());
    jlog("torn", "\"falsifier_and_recovery\":true");
}

#[test]
fn restart_epoch_bump_is_the_aba_guard() {
    // The ABA twin: a lease taken before restart, with the SAME slot and
    // even sequence after restart, must still refuse via the epoch.
    let mut r = ring();
    r.publish(1, &payload(1)).unwrap();
    let old = r.acquire(ABI_VERSION).unwrap();
    r.restart(0xB00);
    // New run publishes into the same slot with fresh ticks.
    r.publish(1, &payload(1)).unwrap();
    let mut buf = vec![0.0; 16];
    let err = r.read(&old, &mut buf).unwrap_err();
    assert_eq!(err.code, "ring-epoch-stale");
    assert_eq!(r.counters().stale_epoch_refusals, 1);
    assert_eq!(r.epoch(), 2);
    // A fresh lease under the new epoch works and carries the new anchor.
    let fresh = r.acquire(ABI_VERSION).unwrap();
    let h = r.read(&fresh, &mut buf).unwrap();
    assert_eq!(h.run_anchor_digest_prefix, 0xB00);
    jlog("aba", "\"epoch_guard\":true");
}

#[test]
fn drop_counters_and_writer_skip_of_leased_slots() {
    let mut r = ring();
    // Publish 5 ticks with no reader: 3 slots -> 2 unread drops.
    for t in 1..=5 {
        r.publish(t, &payload(t)).unwrap();
    }
    assert_eq!(r.counters().dropped_unread, 2, "starvation is counted");
    // Lease the newest; further publishes must SKIP that slot.
    let lease = r.acquire(ABI_VERSION).unwrap();
    for t in 6..=9 {
        r.publish(t, &payload(t)).unwrap();
    }
    assert!(r.counters().writer_skips >= 4, "leased slot skipped");
    let mut buf = vec![0.0; 16];
    // The lease still reads ITS tick (5) — the writer never touched it.
    let h = r.read(&lease, &mut buf).unwrap();
    assert_eq!(h.tick, 5);
    assert_eq!(buf, payload(5));
    r.release(&lease);
    jlog(
        "drops",
        &format!(
            "\"dropped\":{},\"skips\":{}",
            r.counters().dropped_unread,
            r.counters().writer_skips
        ),
    );
}

#[test]
fn publish_validation_and_monotone_ticks() {
    let mut r = ring();
    r.publish(5, &payload(5)).unwrap();
    assert_eq!(
        r.publish(5, &payload(5)).unwrap_err().code,
        "ring-publish-invalid",
        "equal tick refused"
    );
    assert_eq!(
        r.publish(4, &payload(4)).unwrap_err().code,
        "ring-publish-invalid",
        "backwards tick refused"
    );
    let mut bad = payload(6);
    bad[3] = f64::NAN;
    assert_eq!(r.publish(6, &bad).unwrap_err().code, "ring-publish-invalid");
    assert_eq!(
        r.publish(6, &payload(6)[..15]).unwrap_err().code,
        "ring-publish-invalid",
        "short payload refused"
    );
    jlog("publish", "\"validation\":true");
}

#[test]
fn slot_caps_at_both_edges_and_one_past() {
    assert!(SnapshotRing::new(MIN_SLOTS, 8, 0, 0).is_ok());
    assert!(SnapshotRing::new(MAX_SLOTS, 8, 0, 0).is_ok());
    let code = |r: Result<SnapshotRing, fs_flyer_wasm::Refusal>| -> &'static str {
        match r {
            Ok(_) => "OK",
            Err(e) => e.code,
        }
    };
    assert_eq!(
        code(SnapshotRing::new(MIN_SLOTS - 1, 8, 0, 0)),
        "ring-config-invalid"
    );
    assert_eq!(
        code(SnapshotRing::new(MAX_SLOTS + 1, 8, 0, 0)),
        "ring-config-invalid"
    );
    assert_eq!(code(SnapshotRing::new(3, 0, 0, 0)), "ring-config-invalid");
    assert_eq!(
        code(SnapshotRing::new(
            3,
            fs_flyer_wasm::ring::MAX_PAYLOAD_F64 + 1,
            0,
            0
        )),
        "ring-config-invalid"
    );
    jlog("caps", "\"both_edges_and_past\":true");
}

#[test]
fn reads_never_mutate_published_state() {
    // V-18 noninterference at the protocol level: any pattern of
    // acquire/read/release leaves the published payload and header
    // words BITWISE unchanged.
    let mut r = ring();
    r.publish(1, &payload(1)).unwrap();
    r.publish(2, &payload(2)).unwrap();
    let before: Vec<u64> = {
        let lease = r.acquire(ABI_VERSION).unwrap();
        let mut buf = vec![0.0; 16];
        let h = r.read(&lease, &mut buf).unwrap();
        r.release(&lease);
        let mut v: Vec<u64> = buf.iter().map(|x| x.to_bits()).collect();
        v.push(h.sequence);
        v.push(h.tick);
        v
    };
    // Hammer the read path with varied query load.
    for _ in 0..250 {
        let lease = r.acquire(ABI_VERSION).unwrap();
        let mut buf = vec![0.0; 16];
        let _ = r.read(&lease, &mut buf).unwrap();
        r.release(&lease);
    }
    let after: Vec<u64> = {
        let lease = r.acquire(ABI_VERSION).unwrap();
        let mut buf = vec![0.0; 16];
        let h = r.read(&lease, &mut buf).unwrap();
        r.release(&lease);
        let mut v: Vec<u64> = buf.iter().map(|x| x.to_bits()).collect();
        v.push(h.sequence);
        v.push(h.tick);
        v
    };
    assert_eq!(before, after, "250 reads must not move a single bit");
    jlog("noninterference", "\"bitwise_after_250_reads\":true");
}

#[test]
fn golden_digest() {
    let mut r = ring();
    let mut transcript = Vec::new();
    for t in 1..=12u64 {
        let slot = r.publish(t, &payload(t)).unwrap();
        transcript.push(u64::from(slot));
        if t % 3 == 0 {
            let lease = r.acquire(ABI_VERSION).unwrap();
            let mut buf = vec![0.0; 16];
            let h = r.read(&lease, &mut buf).unwrap();
            transcript.push(h.sequence);
            transcript.push(h.tick);
            r.release(&lease);
        }
    }
    let c: RingCounters = r.counters();
    for v in [c.published, c.dropped_unread, c.torn_reads, c.writer_skips] {
        transcript.push(v);
    }
    let mut bytes = Vec::new();
    for v in &transcript {
        bytes.extend_from_slice(&v.to_le_bytes());
    }
    let digest =
        fs_blake3::hash_domain("org.frankensim.fs-flyer-wasm.v18-golden.v1", &bytes).to_hex();
    jlog("golden", &format!("\"digest\":\"{digest}\""));
    assert_eq!(
        digest, "ae52b2ffaa7747665f590bdd9fe9a0fe72eea9bc46d935118d5c929b13774539",
        "ring-protocol golden moved — determinism regression or an \
         intentional ABI change requiring the golden-bump protocol"
    );
}

// Keep Lease in the public test surface.
#[allow(dead_code)]
fn _anchor(l: &Lease) -> u32 {
    l.slot
}
