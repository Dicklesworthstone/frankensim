//! Tombstone-ledger conformance (addendum Proposal E, the lmp4.13 bead).
//! Acceptance: tombstones append automatically on falsification kills +
//! abandoned branches (above threshold); both indexes retrieve
//! near-neighbors; π-space collisions unify dimensionally-equivalent
//! deaths while PRECISION holds (raw-similar but π-different work is NOT
//! blocked); the orchestrator gate blocks re-exploration unless a
//! VALIDATED distinguisher is cited (free text refused, accepted ones
//! accumulate); the re-exploration-rate metric and ledger persistence
//! round out the loop.

use fs_ledger::tombstone::{Descriptor, ExplorationVerdict, TombstoneIndex, pi_distance};
use fs_qty::{Dims, QtyAny};
use std::collections::BTreeMap;

fn verdict(case: &str, detail: &str) {
    println!(
        "{{\"suite\":\"fs-ledger/tombstone\",\"case\":\"{case}\",\"verdict\":\"pass\",\
         \"detail\":\"{detail}\"}}"
    );
}

/// A bracket-in-crossflow descriptor: (ρ, V, L, μ) → one π group (Re).
fn bracket(name: &str, rho: f64, v: f64, l: f64, mu: f64) -> Descriptor {
    let mut params = BTreeMap::new();
    params.insert(
        "density".to_string(),
        QtyAny::new(rho, Dims([-3, 1, 0, 0, 0, 0])),
    );
    params.insert(
        "velocity".to_string(),
        QtyAny::new(v, Dims([1, 0, -1, 0, 0, 0])),
    );
    params.insert(
        "length".to_string(),
        QtyAny::new(l, Dims([1, 0, 0, 0, 0, 0])),
    );
    params.insert(
        "viscosity".to_string(),
        QtyAny::new(mu, Dims([-1, 1, -1, 0, 0, 0])),
    );
    Descriptor {
        name: name.to_string(),
        params,
    }
}

#[test]
fn tb_001_automatic_appends() {
    let mut index = TombstoneIndex::new();
    // Falsification kill: ALWAYS appends, carrying the falsifier JSON.
    let idx = index.record_falsification_kill(
        bracket("aluminum bracket crossflow", 2700.0, 1.1, 0.08, 1.8e-5),
        "{\"kind\":\"tombstone\",\"class\":\"conservation\"}",
        vec!["estimated".to_string()],
        420.0,
        "2026-07-07",
        "agent:CloudyFinch",
    );
    assert_eq!(index.len(), 1);
    let t = index.get(idx).expect("present");
    assert!(
        t.evidence.contains("conservation"),
        "falsifier evidence carried"
    );
    // Abandoned branch: appends ONLY above the cost threshold.
    let below = index.record_abandoned_branch(
        bracket("cheap probe", 1000.0, 0.5, 0.01, 1e-3),
        3.7,
        5.0,   // spent
        100.0, // threshold
        "2026-07-07",
        "agent:x",
    );
    assert!(below.is_none(), "cheap failures are noise, not memory");
    let above = index.record_abandoned_branch(
        bracket("expensive dead end", 1000.0, 0.5, 0.01, 1e-3),
        3.7,
        900.0,
        100.0,
        "2026-07-07",
        "agent:x",
    );
    assert!(above.is_some());
    assert_eq!(index.len(), 2);
    let json = index.get(above.expect("idx")).expect("row").to_json();
    assert!(json.contains("\"kind\":\"tombstone\"") && json.contains("900"));
    verdict(
        "tb-001",
        "kills always append with falsifier evidence; branches append only above threshold",
    );
}

#[test]
fn tb_002_pi_space_unifies_equivalent_deaths_with_precision() {
    // THE proposal example: aluminum at Re 2.0e5 vs steel at Re 2.1e5 —
    // different raw parameters, the SAME death in π-space.
    let aluminum = bracket("aluminum bracket", 1.225, 24.0, 0.12, 1.81e-5); // Re ~ 1.99e5
    let steel = bracket("steel bracket", 1.225, 25.4, 0.119, 1.81e-5); // Re ~ 2.05e5
    let sig_a = aluminum.pi_signature().expect("sig");
    let sig_s = steel.pi_signature().expect("sig");
    let d = pi_distance(&sig_a, &sig_s).expect("same structure");
    assert!(d < 0.1, "±5% Reynolds is the same death: {d} decades");
    let mut index = TombstoneIndex::new();
    index.record_falsification_kill(
        aluminum,
        "{}",
        vec!["estimated".to_string()],
        100.0,
        "2026-07-07",
        "agent:a",
    );
    assert_eq!(
        index.pi_neighbors(&steel),
        vec![0],
        "the steel variant collides with the aluminum tombstone"
    );
    // PRECISION (review-round-3 hardening): raw-similar but π-DIFFERENT
    // work must NOT be blocked. Same raw length/density, but 40x slower —
    // Re drops decades; genuinely novel.
    let creeping = bracket("slow bracket", 1.225, 0.6, 0.12, 1.81e-5); // Re ~ 5e3
    assert!(
        index.pi_neighbors(&creeping).is_empty(),
        "π-different exploration is NOT suppressed"
    );
    assert!(
        matches!(
            index.pre_exploration_check(&creeping),
            ExplorationVerdict::Clear
        ),
        "novel work must fund"
    );
    verdict(
        "tb-002",
        "±5% Re collides across materials; decades-different Re stays fundable \
         (precision fixture)",
    );
}

#[test]
fn tb_003_embedding_index_recall_and_precision() {
    let mut index = TombstoneIndex::new();
    index.record_falsification_kill(
        bracket("lattice infill wing spar", 2700.0, 10.0, 0.4, 1.8e-5),
        "{}",
        vec!["estimated".to_string()],
        50.0,
        "2026-07-07",
        "agent:a",
    );
    // Near-duplicate descriptor (same tokens, same decades): recalled.
    let near = bracket("wing spar lattice infill", 2650.0, 11.0, 0.42, 1.8e-5);
    assert_eq!(
        index.embed_neighbors(&near),
        vec![0],
        "token+decade twin recalled"
    );
    // Distinct problem (different tokens AND decades): not matched.
    let far = bracket("heat sink fin array", 8000.0, 0.02, 0.003, 1.0e-3);
    assert!(
        index.embed_neighbors(&far).is_empty(),
        "distinct work not matched"
    );
    verdict(
        "tb-003",
        "embedding recalls token/decade twins, passes distinct work",
    );
}

#[test]
fn tb_004_gate_blocks_and_validates_distinguishers() {
    let mut index = TombstoneIndex::new();
    index.record_falsification_kill(
        bracket("bracket v1", 1.225, 24.0, 0.12, 1.81e-5),
        "{}",
        vec!["estimated".to_string()],
        100.0,
        "2026-07-07",
        "agent:a",
    );
    // Re-exploring the same neighborhood: BLOCKED via π-space.
    let retry = bracket("bracket v2", 1.225, 24.5, 0.118, 1.81e-5);
    let blocked = index.pre_exploration_check(&retry);
    let neighbor = match blocked {
        ExplorationVerdict::Blocked { ref neighbors, via } => {
            assert_eq!(via, "pi-space");
            neighbors[0]
        }
        ExplorationVerdict::Clear => panic!("must block the re-run"),
    };
    // Free-text / unknown-parameter distinguishers are REFUSED.
    let bogus = index.fund_with_distinguisher(&retry, neighbor, "vibes");
    assert!(bogus.is_err(), "arbitrary text is not a distinguisher");
    // A named parameter that BARELY differs is refused with the delta.
    let same = index.fund_with_distinguisher(&retry, neighbor, "velocity");
    let refusal = same.expect_err("2% velocity is the same death");
    assert!(
        refusal.what.contains("decades"),
        "teaches the threshold: {}",
        refusal.what
    );
    // A genuinely different named parameter funds — and is LOGGED.
    let mut novel = retry.clone();
    novel.params.insert(
        "velocity".to_string(),
        QtyAny::new(90.0, Dims([1, 0, -1, 0, 0, 0])),
    );
    index
        .fund_with_distinguisher(&novel, neighbor, "velocity")
        .expect("3.7x velocity is a real distinguisher");
    let tomb = index.get(neighbor).expect("tombstone");
    assert_eq!(tomb.distinguishers.len(), 1, "distinguishers accumulate");
    assert!(tomb.distinguishers[0].contains("velocity=90"));
    verdict(
        "tb-004",
        "gate blocks; bogus + too-close distinguishers refused with teaching; real ones \
         fund and accumulate on the tombstone",
    );
}

#[test]
fn tb_005_metrics_and_ledger_persistence() {
    let mut index = TombstoneIndex::new();
    index.record_falsification_kill(
        bracket("dead end", 1.225, 24.0, 0.12, 1.81e-5),
        "{}",
        vec!["verified".to_string(), "estimated".to_string()],
        77.0,
        "2026-07-07",
        "agent:a",
    );
    // Two clear checks, one blocked.
    let _ = index.pre_exploration_check(&bracket("fresh 1", 1000.0, 0.001, 0.5, 8.0));
    let _ = index.pre_exploration_check(&bracket("fresh 2", 900.0, 3.0, 2.0, 1e-3));
    let _ = index.pre_exploration_check(&bracket("retry", 1.225, 24.1, 0.12, 1.81e-5));
    let (clear, blocked, funded, rate) = index.re_exploration_rate();
    assert_eq!((clear, blocked, funded), (2, 1, 0));
    assert!(
        (rate - 1.0 / 3.0).abs() < 1e-12,
        "rate = blocked/(clear+blocked)"
    );
    // Ledger persistence: tombstone events written and payload-complete.
    let dir = std::env::temp_dir().join(format!("fs-tombstone-{}", std::process::id()));
    std::fs::create_dir_all(&dir).expect("mkdir");
    let ledger =
        fs_ledger::Ledger::open(dir.join("t.led").to_str().expect("utf8")).expect("ledger");
    index.flush_to_ledger(&ledger).expect("flush");
    assert_eq!(ledger.table_count("events").expect("count"), 1);
    let json = index.get(0).expect("row").to_json();
    for needle in ["dead end", "77", "verified", "agent:a", "2026-07-07"] {
        assert!(json.contains(needle), "payload carries {needle:?}");
    }
    let _ = std::fs::remove_dir_all(&dir);
    verdict(
        "tb-005",
        "re-exploration rate exact; tombstone events persisted with full payloads",
    );
}

/// Bead sj31i.30: π comparison demands SEMANTIC basis compatibility.
/// Identical dimension matrices and group values under different role
/// names make no collision claim; verified crosswalk receipts restore
/// comparability for genuine aliases; partial/many-to-one/dims-changing
/// crosswalks and stale receipts refuse.
#[test]
fn tb_010_semantic_basis_compatibility_gates_pi_comparison() {
    use fs_ledger::tombstone::{
        BasisCrosswalk, PI_COMPARE_POLICY_VERSION, pi_distance_crosswalked,
    };

    assert_eq!(PI_COMPARE_POLICY_VERSION, 2);

    // A "speed"-named clone of the bracket basis: same dims matrix,
    // same group values — dimension-only agreement.
    let renamed = |name: &str, rho: f64, v: f64, l: f64, mu: f64| -> Descriptor {
        let mut params = BTreeMap::new();
        params.insert(
            "rho".to_string(),
            QtyAny::new(rho, Dims([-3, 1, 0, 0, 0, 0])),
        );
        params.insert(
            "speed".to_string(),
            QtyAny::new(v, Dims([1, 0, -1, 0, 0, 0])),
        );
        params.insert(
            "chord".to_string(),
            QtyAny::new(l, Dims([1, 0, 0, 0, 0, 0])),
        );
        params.insert(
            "mu".to_string(),
            QtyAny::new(mu, Dims([-1, 1, -1, 0, 0, 0])),
        );
        Descriptor {
            name: name.to_string(),
            params,
        }
    };

    let a = bracket("bracket", 1.225, 24.0, 0.12, 1.81e-5)
        .pi_signature()
        .expect("sig");
    let b = renamed("alien schema", 1.225, 24.0, 0.12, 1.81e-5)
        .pi_signature()
        .expect("sig");

    // Identical exponents and values — but incompatible semantics.
    assert_eq!(
        pi_distance(&a, &b),
        None,
        "a matching dimension matrix alone never makes π-neighbors"
    );

    // A dimension-only twin never becomes a PRIMARY neighbor in the
    // index either; with distinct name tokens the gate stays Clear.
    let mut index = TombstoneIndex::new();
    index.record_falsification_kill(
        bracket("bracket crossflow study", 1.225, 24.0, 0.12, 1.81e-5),
        "{}",
        vec!["estimated".to_string()],
        100.0,
        "2026-07-07",
        "agent:a",
    );
    let probe = renamed("unrelated duct jig", 1.225, 24.0, 0.12, 1.81e-5);
    assert!(
        index.pi_neighbors(&probe).is_empty(),
        "no π suppression across incompatible bases"
    );
    assert!(matches!(
        index.pre_exploration_check(&probe),
        ExplorationVerdict::Clear
    ));

    // The tombstone evidence payload binds the comparison policy.
    let json = index.get(0).expect("row").to_json();
    assert!(
        json.contains("\"pi_policy\":2"),
        "distance policy bound into evidence: {json}"
    );

    // An explicit, verified crosswalk restores comparability.
    let mut map = BTreeMap::new();
    map.insert("chord".to_string(), "length".to_string());
    map.insert("mu".to_string(), "viscosity".to_string());
    map.insert("rho".to_string(), "density".to_string());
    map.insert("speed".to_string(), "velocity".to_string());
    let crosswalk = BasisCrosswalk {
        from_schema: "alien.v1".to_string(),
        to_schema: "bracket.v1".to_string(),
        map: map.clone(),
    };
    let receipt = crosswalk
        .verify(&b, &a)
        .expect("total bijective dims-preserving alias");
    let d = pi_distance_crosswalked(&b, &a, &receipt).expect("crosswalked aliases are comparable");
    assert!(d < 1e-12, "identical physics through the alias: {d}");

    // Stale receipt: minted for (b, a), presented with a DIFFERENT
    // from-signature — refuses, no claim.
    let drifted = renamed("alien schema", 1.225, 240.0, 0.12, 1.81e-5)
        .pi_signature()
        .expect("sig");
    assert!(
        pi_distance_crosswalked(&drifted, &a, &receipt).is_some(),
        "same basis, different values: still bound (values are the distance)"
    );
    let mut drifted_basis = renamed("alien schema", 1.225, 24.0, 0.12, 1.81e-5);
    drifted_basis.params.insert(
        "spin".to_string(),
        QtyAny::new(1.0, Dims([0, 0, -1, 0, 0, 0])),
    );
    let drifted_basis = drifted_basis.pi_signature().expect("sig");
    assert_eq!(
        pi_distance_crosswalked(&drifted_basis, &a, &receipt),
        None,
        "a receipt does not bind a basis that changed since minting"
    );

    // Partial map refuses.
    let mut partial = map.clone();
    partial.remove("speed");
    assert!(
        BasisCrosswalk {
            from_schema: "alien.v1".into(),
            to_schema: "bracket.v1".into(),
            map: partial,
        }
        .verify(&b, &a)
        .is_err(),
        "partial crosswalks refuse"
    );

    // Many-to-one refuses.
    let mut collapsing = map.clone();
    collapsing.insert("speed".to_string(), "density".to_string());
    assert!(
        BasisCrosswalk {
            from_schema: "alien.v1".into(),
            to_schema: "bracket.v1".into(),
            map: collapsing,
        }
        .verify(&b, &a)
        .is_err(),
        "many-to-one crosswalks refuse"
    );

    // A dims-changing alias refuses even when names line up.
    let mut swapped = map;
    swapped.insert("speed".to_string(), "length".to_string());
    swapped.insert("chord".to_string(), "velocity".to_string());
    assert!(
        BasisCrosswalk {
            from_schema: "alien.v1".into(),
            to_schema: "bracket.v1".into(),
            map: swapped,
        }
        .verify(&b, &a)
        .is_err(),
        "an alias that changes coordinate dimensions refuses"
    );

    verdict(
        "tb-010",
        "dimension-only matches make no claim; verified crosswalks restore comparability; \
         partial/many-to-one/dims-changing/stale receipts refuse",
    );
}

/// Bead sj31i.30 metric laws: reordered parameter insertion is
/// canonicalized, distance is symmetric, satisfies the triangle
/// inequality, and shifts metamorphically by exactly the rescaling
/// decade.
#[test]
fn tb_011_pi_distance_metric_laws_hold_on_compatible_bases() {
    let a = bracket("a", 1.225, 10.0, 0.1, 1.8e-5)
        .pi_signature()
        .expect("sig");
    let b = bracket("b", 1.225, 20.0, 0.1, 1.8e-5)
        .pi_signature()
        .expect("sig");
    let c = bracket("c", 1.225, 80.0, 0.1, 1.8e-5)
        .pi_signature()
        .expect("sig");

    // Reordered insertion canonicalizes to the same signature.
    let mut reordered_params = BTreeMap::new();
    for (key, qty) in bracket("a", 1.225, 10.0, 0.1, 1.8e-5)
        .params
        .into_iter()
        .rev()
    {
        reordered_params.insert(key, qty);
    }
    let reordered = Descriptor {
        name: "a".to_string(),
        params: reordered_params,
    }
    .pi_signature()
    .expect("sig");
    assert_eq!(a, reordered, "insertion order cannot change the basis");

    let dab = pi_distance(&a, &b).expect("compatible");
    let dba = pi_distance(&b, &a).expect("compatible");
    let dbc = pi_distance(&b, &c).expect("compatible");
    let dac = pi_distance(&a, &c).expect("compatible");
    assert_eq!(dab.to_bits(), dba.to_bits(), "symmetry");
    assert!(dac <= dab + dbc + 1e-12, "triangle inequality");

    // Metamorphic: velocity x10 shifts the Reynolds group by one decade.
    let shifted = bracket("s", 1.225, 100.0, 0.1, 1.8e-5)
        .pi_signature()
        .expect("sig");
    let base = bracket("s", 1.225, 10.0, 0.1, 1.8e-5)
        .pi_signature()
        .expect("sig");
    let d = pi_distance(&base, &shifted).expect("compatible");
    assert!((d - 1.0).abs() < 1e-9, "x10 velocity is one decade: {d}");

    verdict(
        "tb-011",
        "canonical ordering, symmetry, triangle, and decade metamorphic all hold",
    );
}
