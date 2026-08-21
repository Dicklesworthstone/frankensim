//! E10.2b battery (bead wf-root-guzez.11.4): a SMALL PILOT CAMPAIGN
//! runs shard-resumably through the REAL lifecycle engine with a
//! deterministic partition-independent merge digest; refusal
//! accounting exact (a below-floor headwind flight is a RECORDED
//! refusal row); the surrogate-in-finals hostile twin refuses; shard
//! gap/overlap/mismatch refusals; caps at cap AND cap+1; no scoring
//! fields exist on the receipt.
//! Repro: cargo test -p fs-flyer --test hcampaign_battery

use fs_flyer::hcampaign::{
    CampaignRow, FlightSpec, HistoricalCampaignIntentManifestV1, MAX_FLIGHTS, MAX_MEMBERS,
    MemberSpec, RunOutcome, ShardResult, execute_shard, merge_shards,
};

fn jlog(case: &str, payload: &str) {
    println!("{{\"suite\":\"fs-flyer-hcampaign\",\"case\":\"{case}\",{payload}}}");
}

/// The pilot campaign: 2 exact members + 1 surrogate, 2 flights (one
/// healthy Dec-17 class, one whose headwind is below the aero floor
/// so its init REFUSES — the accounting fixture). Short tick budgets
/// keep the battery fast while every run is the REAL engine.
fn pilot_intent() -> HistoricalCampaignIntentManifestV1 {
    HistoricalCampaignIntentManifestV1 {
        label: "pilot-campaign-v1",
        members: vec![
            MemberSpec {
                member: 0,
                seed: 1903,
                surrogate: false,
            },
            MemberSpec {
                member: 1,
                seed: 1904,
                surrogate: false,
            },
            MemberSpec {
                member: 7,
                seed: 9999,
                surrogate: true,
            },
        ],
        flights: vec![
            FlightSpec {
                flight: 0,
                headwind_mps: 11.0,
                rho_kg_m3: 1.294,
                max_ticks: 240,
            },
            FlightSpec {
                flight: 1,
                headwind_mps: 0.2, // below the admission floor: refuses
                rho_kg_m3: 1.294,
                max_ticks: 240,
            },
        ],
    }
}

#[test]
fn pilot_campaign_is_shard_resumable_with_partition_independent_merge() {
    let intent = pilot_intent();
    let n = intent.work_units().len();
    assert_eq!(n, 4, "2 exact members x 2 flights (surrogate excluded)");
    // Partition A: one shard.
    let whole = execute_shard(&intent, 0, n).unwrap();
    let merged_a = merge_shards(&intent, &[whole.clone()]).unwrap();
    // Partition B: three shards (resume shape: 1 + 2 + 1).
    let s1 = execute_shard(&intent, 0, 1).unwrap();
    let s2 = execute_shard(&intent, 1, 3).unwrap();
    let s3 = execute_shard(&intent, 3, 4).unwrap();
    let merged_b = merge_shards(&intent, &[s3, s1, s2]).unwrap();
    assert_eq!(
        merged_a.merge_digest, merged_b.merge_digest,
        "merge digest is partition-independent (shards are resumable)"
    );
    // Refusal accounting: flight 1 refused for BOTH exact members,
    // flight 0 completed for both — per-row, never totals-only.
    assert_eq!(merged_a.rows.len(), 4);
    assert_eq!(merged_a.completed, 2);
    assert_eq!(merged_a.refused, 2);
    for row in &merged_a.rows {
        match (row.flight, &row.outcome) {
            (0, RunOutcome::Completed { terminal, digest }) => {
                assert!(!terminal.is_empty());
                assert_eq!(digest.len(), 64, "lifecycle digest carried");
            }
            (1, RunOutcome::Refused { code }) => {
                assert!(!code.is_empty(), "typed refusal recorded");
            }
            other => panic!("row misfiled: member {} {:?}", row.member, other),
        }
    }
    // Member-index merge: member-major order.
    let order: Vec<(u32, u32)> = merged_a.rows.iter().map(|r| (r.member, r.flight)).collect();
    assert_eq!(order, vec![(0, 0), (0, 1), (1, 0), (1, 1)]);
    jlog(
        "pilot",
        &format!(
            "\"merge_digest\":\"{}\",\"completed\":{},\"refused\":{}",
            merged_a.merge_digest, merged_a.completed, merged_a.refused
        ),
    );
}

#[test]
fn surrogate_in_finals_hostile_twin_refuses() {
    let intent = pilot_intent();
    let whole = execute_shard(&intent, 0, 4).unwrap();
    // Forge a finals row for the surrogate member 7 (the attack).
    let mut forged = whole.clone();
    forged.rows.push(CampaignRow {
        member: 7,
        flight: 0,
        outcome: RunOutcome::Completed {
            terminal: "max-ticks",
            digest: "f".repeat(64),
        },
    });
    let err = merge_shards(&intent, &[forged]).unwrap_err();
    assert_eq!(err.code, "surrogate-in-finals");
    // The executor itself NEVER runs surrogates: no work unit names
    // member 7.
    assert!(intent.work_units().iter().all(|(m, _)| m.member != 7));
    jlog("surrogate-twin", &format!("\"code\":\"{}\"", err.code));
}

#[test]
fn shard_coverage_and_identity_refusals() {
    let intent = pilot_intent();
    let s1 = execute_shard(&intent, 0, 2).unwrap();
    // Gap: missing [2, 4).
    let gap = merge_shards(&intent, &[s1.clone()]).unwrap_err();
    assert_eq!(gap.code, "campaign-shards-incomplete");
    // Overlap: [0,2) twice.
    let overlap = merge_shards(&intent, &[s1.clone(), s1.clone()]).unwrap_err();
    assert_eq!(overlap.code, "campaign-shards-incomplete");
    // RESUME: supply the missing slice and the merge completes.
    let s2 = execute_shard(&intent, 2, 4).unwrap();
    assert!(merge_shards(&intent, &[s1.clone(), s2]).is_ok());
    // Mismatched intent: a shard from a different campaign refuses.
    let mut other = pilot_intent();
    other.label = "another-campaign";
    let foreign = execute_shard(&other, 0, 4).unwrap();
    let err = merge_shards(&intent, &[foreign]).unwrap_err();
    assert_eq!(err.code, "campaign-shards-mismatched");
    // Shard range: AT the end admits (tested above); one past refuses;
    // empty refuses.
    assert_eq!(
        execute_shard(&intent, 0, 5).unwrap_err().code,
        "campaign-shard-invalid"
    );
    assert_eq!(
        execute_shard(&intent, 2, 2).unwrap_err().code,
        "campaign-shard-invalid"
    );
    jlog(
        "shards",
        "\"gap\":\"refused\",\"overlap\":\"refused\",\"resume\":\"ok\"",
    );
}

#[test]
fn intent_caps_and_validity() {
    let mk = |nm: usize, nf: usize| HistoricalCampaignIntentManifestV1 {
        label: "caps",
        members: (0..nm)
            .map(|i| MemberSpec {
                member: i as u32,
                seed: i as u64,
                surrogate: false,
            })
            .collect(),
        flights: (0..nf)
            .map(|i| FlightSpec {
                flight: i as u32,
                headwind_mps: 11.0,
                rho_kg_m3: 1.294,
                max_ticks: 10,
            })
            .collect(),
    };
    assert!(mk(MAX_MEMBERS, MAX_FLIGHTS).admit().is_ok(), "AT caps");
    assert_eq!(
        mk(MAX_MEMBERS + 1, 1).admit().unwrap_err().code,
        "campaign-intent-invalid"
    );
    assert_eq!(
        mk(1, MAX_FLIGHTS + 1).admit().unwrap_err().code,
        "campaign-intent-invalid"
    );
    assert_eq!(
        mk(0, 1).admit().unwrap_err().code,
        "campaign-intent-invalid"
    );
    // All-surrogate campaigns refuse (nothing would execute).
    let mut all_surrogate = mk(2, 1);
    for m in &mut all_surrogate.members {
        m.surrogate = true;
    }
    assert_eq!(
        all_surrogate.admit().unwrap_err().code,
        "campaign-intent-invalid"
    );
    // Duplicate member ids refuse.
    let mut dup = mk(2, 1);
    dup.members[1].member = 0;
    assert_eq!(dup.admit().unwrap_err().code, "campaign-intent-invalid");
    // Intent digest moves with a member seed (identity is real).
    let a = mk(2, 1).admit().unwrap();
    let mut edited = mk(2, 1);
    edited.members[0].seed = 42;
    assert_ne!(a, edited.admit().unwrap());
    jlog(
        "caps",
        &format!("\"members\":{MAX_MEMBERS},\"flights\":{MAX_FLIGHTS}"),
    );
}
