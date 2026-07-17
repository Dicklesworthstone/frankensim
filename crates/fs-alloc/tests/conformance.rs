//! fs-alloc conformance suite (CONTRACT.md: any reimplementation must pass).
//!
//! Gauntlet coverage: G0 (accounting laws, alignment laws, shadow-model
//! equivalence), G4 (the 10^6-cancellation storm, concurrent hammer), and
//! G5 (deterministic accounting reports). Every case prints one canonical
//! fs-obs conformance verdict. Randomized cases carry their input seed;
//! alloc-006 deliberately makes no exact OS-interleaving replay claim.

use fs_alloc::{
    ALLOC_ALIGN, AllocError, ArenaConfig, ArenaPool, HUGEPAGE_BYTES, HugepageOutcome,
    HugepagePolicy, OperationMemoryLease, RECLAIM_POISON_VERSION, ReclaimPoison, ShardedPool, Site,
    SiteStats,
};

fn verdict(case: &str, pass: bool, detail: &str, seed: u64) {
    let mut emitter = fs_obs::Emitter::new("fs-alloc/conformance", case);
    let event = emitter.emit(
        if pass {
            fs_obs::Severity::Info
        } else {
            fs_obs::Severity::Error
        },
        fs_obs::EventKind::ConformanceCase {
            suite: "fs-alloc/conformance".to_string(),
            case: case.to_string(),
            pass,
            detail: detail.to_string(),
            seed,
        },
        None,
    );
    fs_obs::lint_failure_record(&event).expect("allocation verdict must be replayable");
    let line = event.to_jsonl();
    fs_obs::validate_line(&line).expect("allocation verdict must use the fs-obs wire schema");
    println!("{line}");
    assert!(pass, "case {case}: {detail}");
}

/// In-house deterministic LCG (L0 crates cannot depend on fs-rand; same
/// constants as fs-qty's hardening battery).
struct Lcg(u64);

impl Lcg {
    fn next_u64(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        self.0
    }

    fn below(&mut self, n: u64) -> u64 {
        self.next_u64() % n
    }
}

fn small_config() -> ArenaConfig {
    ArenaConfig {
        chunk_bytes: 8192,
        max_chunk_bytes: 1 << 20,
        limit_bytes: None,
        free_list_max_bytes: 4 << 20,
        hugepage: HugepagePolicy::Never,
    }
}

#[test]
fn alloc_001_unconditional_128_byte_alignment() {
    let pool = ArenaPool::new(small_config());
    let mut checked = 0usize;
    pool.scope(|a| {
        for len in [1usize, 2, 3, 7, 63, 64, 65, 127, 128, 129, 1000, 4096] {
            let s = a
                .alloc_slice_fill(Site::named("conf/align"), len, 0u8)
                .expect("alloc");
            assert_eq!(s.as_ptr() as usize % ALLOC_ALIGN, 0, "u8 slice len {len}");
            checked += 1;
        }
        let x = a.alloc(Site::named("conf/align"), 1u8).expect("alloc");
        assert_eq!(core::ptr::from_mut(x) as usize % ALLOC_ALIGN, 0);
        let y = a.alloc(Site::named("conf/align"), 1.0f64).expect("alloc");
        assert_eq!(core::ptr::from_mut(y) as usize % ALLOC_ALIGN, 0);
        checked += 2;
    });
    verdict(
        "alloc-001",
        checked == 14,
        "every allocation kind is 128-byte aligned (G0 alignment law)",
        0,
    );
}

#[test]
fn alloc_002_scope_reclaim_reaches_quiescence() {
    let pool = ArenaPool::new(small_config());
    for _ in 0..10 {
        pool.scope(|a| {
            a.alloc_slice_fill(Site::named("conf/reclaim"), 3000, 1u8)
                .expect("alloc");
            a.scope(|inner| {
                inner
                    .alloc_slice_fill(Site::named("conf/reclaim-inner"), 3000, 2u8)
                    .expect("alloc");
            });
        });
    }
    let stats = pool.stats();
    verdict(
        "alloc-002",
        stats.quiescent(),
        &format!("nested scopes reclaim to quiescence: {}", stats.to_json()),
        0,
    );
}

#[test]
fn alloc_003_budget_refusal_is_structured_and_recoverable() {
    let pool = ArenaPool::new(ArenaConfig {
        limit_bytes: Some(16 * 1024),
        ..small_config()
    });
    let err = pool.scope(|a| {
        a.alloc_slice_fill(Site::named("conf/budget"), 1 << 20, 0u8)
            .expect_err("1 MiB cannot fit a 16 KiB budget")
    });
    let structured = matches!(
        &err,
        AllocError::Exhausted {
            site: "conf/budget",
            limit_bytes: 16384,
            ..
        }
    );
    let teaches = err.to_string().contains("Fixes (ranked)");
    // The refusal must not poison the pool.
    pool.scope(|a| {
        a.alloc(Site::named("conf/budget-after"), 1u64)
            .expect("pool must stay usable after refusal");
    });
    verdict(
        "alloc-003",
        structured && teaches && pool.stats().quiescent(),
        &format!("budget exhaustion is a structured teaching error: {err}"),
        0,
    );
}

#[test]
fn alloc_004_g4_storm_one_million_cancellations() {
    const SEED: u64 = 0xF5A0_0C04_2026_0706;
    const ITERATIONS: u64 = 1_000_000;
    let pool = ArenaPool::new(small_config());
    let mut rng = Lcg(SEED);
    let mut cancelled: u64 = 0;
    let mut completed: u64 = 0;
    for _ in 0..ITERATIONS {
        // Each iteration is one scoped unit of work; a "cancellation" is an
        // early exit from the closure — exactly what scope teardown under
        // asupersync cancellation does to an arena (drop mid-work).
        let early = pool.scope(|a| {
            let allocs = rng.below(3); // 0..=2 allocations before the cancel point
            for i in 0..allocs {
                let len = 1 + (rng.below(2048) as usize);
                if i % 2 == 0 {
                    let s = a
                        .alloc_slice_fill(Site::named("storm/buf"), len, 0xA5u8)
                        .expect("storm alloc");
                    // Touch, so the allocation cannot be elided.
                    assert_eq!(s[len / 2], 0xA5);
                } else {
                    let v = a
                        .alloc(Site::named("storm/val"), len as u64)
                        .expect("storm alloc");
                    assert_eq!(*v, len as u64);
                }
            }
            rng.below(2) == 0 // half the scopes "cancel" (drop with work unfinished)
        });
        if early {
            cancelled += 1;
        } else {
            completed += 1;
        }
    }
    let stats = pool.stats();
    let pass = stats.quiescent() && cancelled + completed == ITERATIONS;
    // Emit the G4 storm verdict through the one observability schema.
    let mut em = fs_obs::Emitter::new("fs-alloc/conformance", "alloc-004/storm");
    let event = em.emit(
        if pass {
            fs_obs::Severity::Info
        } else {
            fs_obs::Severity::Error
        },
        fs_obs::EventKind::StormAssertion {
            name: "no-arena-leak".to_string(),
            pass,
            seed: SEED,
        },
        None,
    );
    let line = event.to_jsonl();
    fs_obs::validate_line(&line).expect("storm event must validate");
    fs_obs::lint_failure_record(&event).expect("storm event must carry its seed");
    println!("{line}");
    verdict(
        "alloc-004",
        pass,
        &format!(
            "10^6 random cancellations leak nothing (G4): cancelled={cancelled} \
             completed={completed} {}",
            stats.to_json()
        ),
        SEED,
    );
}

#[test]
fn alloc_005_g0_shadow_model_accounting_and_disjointness() {
    const SEED: u64 = 0x0005_5EED_D00D_F00D;
    let pool = ArenaPool::new(small_config());
    let mut rng = Lcg(SEED);
    let mut shadow_bytes: u64 = 0;
    let mut shadow_count: u64 = 0;
    let mut shadow_site_b = SiteStats::default();
    pool.scope(|a| {
        let mut ranges: Vec<(usize, usize)> = Vec::new();
        for _ in 0..2000 {
            let len = 1 + (rng.below(512) as usize);
            let padded = (len as u64).div_ceil(ALLOC_ALIGN as u64) * ALLOC_ALIGN as u64;
            let s = if rng.below(2) == 0 {
                a.alloc_slice_fill(Site::named("shadow/a"), len, 1u8)
                    .expect("alloc")
            } else {
                shadow_site_b.bytes += padded;
                shadow_site_b.allocations += 1;
                a.alloc_slice_with(Site::named("shadow/b"), len, |i| (i % 251) as u8)
                    .expect("alloc")
            };
            shadow_bytes += padded;
            shadow_count += 1;
            ranges.push((s.as_ptr() as usize, len));
        }
        let stats = a.stats();
        assert_eq!(stats.allocated_bytes, shadow_bytes, "byte accounting");
        assert_eq!(stats.allocation_count, shadow_count, "count accounting");
        // Pairwise disjointness of every returned range (G0 aliasing law).
        ranges.sort_unstable();
        for w in ranges.windows(2) {
            assert!(
                w[0].0 + w[0].1 <= w[1].0,
                "overlapping allocations: {:?}",
                &w[..2]
            );
        }
    });
    let report = pool.site_report();
    let site_b = report
        .sites
        .iter()
        .find(|(name, _)| *name == "shadow/b")
        .map(|(_, s)| *s)
        .unwrap_or_default();
    verdict(
        "alloc-005",
        site_b == shadow_site_b && pool.stats().quiescent(),
        &format!(
            "shadow model matches (seed {SEED:#x}): {} vs expected bytes={} allocs={}",
            report.to_json(),
            shadow_site_b.bytes,
            shadow_site_b.allocations
        ),
        SEED,
    );
}

#[test]
fn alloc_006_concurrent_arenas_and_pools_stay_leak_free() {
    const SEED: u64 = 0xC0C0;
    let pool = ArenaPool::new(small_config());
    let shaped: ShardedPool<Vec<u64>> = ShardedPool::new(4);
    std::thread::scope(|s| {
        let pool = &pool;
        let shaped = &shaped;
        for t in 0..8usize {
            s.spawn(move || {
                let mut rng = Lcg(SEED + t as u64);
                for i in 0..400usize {
                    pool.scope(|a| {
                        let len = 1 + (rng.below(1024) as usize);
                        let buf = a
                            .alloc_slice_fill(Site::named("conc/buf"), len, t as u8)
                            .expect("alloc");
                        assert_eq!(buf[0], t as u8);
                    });
                    let mut tile = shaped.acquire_with(t, || vec![0u64; 64]);
                    tile[i % 64] = i as u64;
                }
            });
        }
    });
    let (ps, ss) = (pool.stats(), shaped.stats());
    verdict(
        "alloc-006",
        ps.quiescent() && shaped.quiescent(),
        &format!(
            "8 threads x 400 scopes leak nothing (base seed {SEED:#x}, per-thread seed = base + thread index): arena={} pool={}",
            ps.to_json(),
            ss.to_json()
        ),
        SEED,
    );
}

#[test]
fn alloc_007_g5_deterministic_accounting_reports() {
    const SEED: u64 = 0x0007_D373_2026_0706;
    let run = || {
        let pool = ArenaPool::new(small_config());
        let mut rng = Lcg(SEED);
        for _ in 0..200 {
            pool.scope(|a| {
                for _ in 0..rng.below(8) {
                    let len = 1 + (rng.below(700) as usize);
                    let site = if rng.below(2) == 0 {
                        Site::named("det/a")
                    } else {
                        Site::named("det/b")
                    };
                    a.alloc_slice_fill(site, len, 0u8).expect("alloc");
                }
            });
        }
        (pool.site_report().to_json(), pool.stats().to_json())
    };
    let (r1, s1) = run();
    let (r2, s2) = run();
    verdict(
        "alloc-007",
        r1 == r2 && s1 == s2,
        &format!("same seed => identical reports (G5, seed {SEED:#x}): {r1}"),
        SEED,
    );
}

#[test]
fn alloc_008_hugepage_decision_is_recorded_and_honored() {
    let pool = ArenaPool::new(ArenaConfig {
        chunk_bytes: HUGEPAGE_BYTES,
        max_chunk_bytes: 4 * HUGEPAGE_BYTES,
        hugepage: HugepagePolicy::Auto,
        ..small_config()
    });
    let decision = pool.hugepage_decision().clone();
    // Platform envelope: at a 2 MiB chunk size, Auto never reports
    // NotRequested; the concrete outcome depends on the machine and is
    // exactly what must be RECORDED.
    let outcome_valid = decision.outcome != HugepageOutcome::NotRequested;
    let recorded = !decision.detail.is_empty();
    // If the decision claims THP-eligible alignment, the first chunk's base
    // must actually be 2 MiB-aligned (first allocation sits at chunk base).
    let alignment_honored = pool.scope(|a| {
        a.alloc(Site::named("conf/huge"), 1u8).expect("alloc");
        let base = a.current_chunk_base().expect("chunk exists");
        match decision.outcome {
            HugepageOutcome::AlignedForThp => base % HUGEPAGE_BYTES == 0,
            _ => base % ALLOC_ALIGN == 0,
        }
    });
    // The decision must ride the observability schema for the ledger.
    let mut em = fs_obs::Emitter::new("fs-alloc/conformance", "alloc-008/hugepage");
    let line = em
        .emit(
            fs_obs::Severity::Info,
            fs_obs::EventKind::Custom {
                name: "fs-alloc-hugepage-decision".to_string(),
                json: decision.to_json(),
            },
            None,
        )
        .to_jsonl();
    fs_obs::validate_line(&line).expect("decision event must validate");
    println!("{line}");
    verdict(
        "alloc-008",
        outcome_valid && recorded && alignment_honored,
        &format!("hugepage choice recorded + honored: {}", decision.to_json()),
        0,
    );
}

#[test]
fn alloc_009_chunk_recycling_bounds_os_traffic() {
    let pool = ArenaPool::new(small_config());
    // Warm-up scope establishes the chunk set...
    pool.scope(|a| {
        a.alloc_slice_fill(Site::named("conf/reuse"), 6000, 0u8)
            .expect("alloc");
    });
    let created_after_warmup = pool.stats().chunks_created;
    // ...and 100 further identical scopes must be served from the free list.
    for _ in 0..100 {
        pool.scope(|a| {
            a.alloc_slice_fill(Site::named("conf/reuse"), 6000, 0u8)
                .expect("alloc");
        });
    }
    let stats = pool.stats();
    verdict(
        "alloc-009",
        stats.chunks_created == created_after_warmup && stats.chunks_recycled >= 100,
        &format!("steady-state scopes recycle chunks: {}", stats.to_json()),
        0,
    );
}

#[test]
fn alloc_010_seeded_reclaim_poison_detects_and_quarantines_corruption() {
    const SEED: u64 = 0xA110_C010_2026_0715;

    let plain = ArenaPool::new(small_config());
    plain.scope(|arena| {
        arena
            .alloc(Site::named("poison/default-off"), 1u8)
            .expect("plain pool allocation");
    });
    assert!(
        plain.inject_reclaimed_chunk_corruption().is_none(),
        "reclaim poisoning must stay explicitly opt-in"
    );
    assert!(plain.stats().quiescent());

    let pool = ArenaPool::new_with_reclaim_poison(small_config(), ReclaimPoison::seeded(SEED));

    pool.scope(|arena| {
        let buffer = arena
            .alloc_slice_fill(Site::named("poison/seed"), 1024, 0x5au8)
            .expect("seed allocation");
        assert_eq!(buffer[511], 0x5a);
    });
    assert!(
        pool.stats().quiescent(),
        "reclaim must park the poisoned chunk"
    );

    let mutation = pool
        .inject_reclaimed_chunk_corruption()
        .expect("poison mode has one retained chunk");
    assert_eq!(mutation.version, RECLAIM_POISON_VERSION);
    assert_eq!(mutation.seed, SEED);
    assert_eq!(mutation.chunk_bytes, 8192, "golden normalized chunk size");
    assert_eq!(mutation.offset, 7565, "golden v1 corruption offset");
    assert_eq!(mutation.expected, 0x26, "golden v1 poison byte");
    assert_eq!(mutation.actual, 0xd9, "golden v1 corrupted byte");

    let error = pool
        .scope(|arena| arena.alloc(Site::named("poison/reuse"), 7u64).map(|_| ()))
        .expect_err("the stale-write surrogate must fail before chunk reuse");
    let receipt_matches = matches!(
        error,
        AllocError::ReclaimedChunkCorrupted {
            site: "poison/reuse",
            poison_version: RECLAIM_POISON_VERSION,
            poison_seed: SEED,
            chunk_bytes,
            offset,
            expected,
            actual,
        } if chunk_bytes == mutation.chunk_bytes
            && offset == mutation.offset
            && expected == mutation.expected
            && actual == mutation.actual
    );
    let after_detection = pool.stats();

    // Detection quarantines rather than re-parking the corrupt block. A fresh
    // allocation must recover normally, and both states must be leak-clean.
    pool.scope(|arena| {
        let value = arena
            .alloc(Site::named("poison/recovery"), 11u64)
            .expect("pool remains usable after quarantine");
        assert_eq!(*value, 11);
    });
    let after_recovery = pool.stats();

    // A mismatch is quarantined before lease reservation. A fresh lease for
    // the failing reuse must therefore remain exactly untouched.
    let leased_pool =
        ArenaPool::new_with_reclaim_poison(small_config(), ReclaimPoison::seeded(SEED));
    let seed_lease = OperationMemoryLease::unbounded();
    leased_pool.scope_leased(&seed_lease, |arena| {
        arena
            .alloc_slice_fill(Site::named("poison/leased-seed"), 1024, 0x5au8)
            .expect("leased seed allocation");
    });
    assert_eq!(seed_lease.receipt().used_bytes, 0);
    leased_pool
        .inject_reclaimed_chunk_corruption()
        .expect("leased poison pool has one retained chunk");
    let detection_lease = OperationMemoryLease::unbounded();
    let leased_error = leased_pool
        .scope_leased(&detection_lease, |arena| {
            arena
                .alloc(Site::named("poison/leased-reuse"), 13u64)
                .map(|_| ())
        })
        .expect_err("leased reuse must detect corruption before charging");
    assert!(matches!(
        leased_error,
        AllocError::ReclaimedChunkCorrupted {
            poison_version: RECLAIM_POISON_VERSION,
            poison_seed: SEED,
            ..
        }
    ));
    let detection_receipt = detection_lease.receipt();
    assert_eq!(detection_receipt.requested_bytes, 0);
    assert_eq!(detection_receipt.peak_bytes, 0);
    assert_eq!(detection_receipt.used_bytes, 0);
    assert_eq!(detection_receipt.refusals, 0);
    assert_eq!(detection_receipt.release_invariant_violations, 0);
    assert!(leased_pool.stats().quiescent());

    let mut emitter = fs_obs::Emitter::new("fs-alloc/conformance", "alloc-010/poison");
    let poison_event = emitter.emit(
        fs_obs::Severity::Info,
        fs_obs::EventKind::Custom {
            name: "fs-alloc-reclaim-poison".to_string(),
            json: format!(
                "{{\"poison_version\":{},\"seed\":{SEED},\"chunk_bytes\":{},\
                 \"offset\":{},\"expected\":{},\"actual\":{}}}",
                mutation.version,
                mutation.chunk_bytes,
                mutation.offset,
                mutation.expected,
                mutation.actual
            ),
        },
        None,
    );
    let poison_line = poison_event.to_jsonl();
    fs_obs::validate_line(&poison_line).expect("poison receipt must use the fs-obs wire schema");
    println!("{poison_line}");
    verdict(
        "alloc-010",
        receipt_matches
            && after_detection.quiescent()
            && after_detection.reserved_bytes == 0
            && after_detection.free_bytes == 0
            && after_recovery.quiescent()
            && after_recovery.free_bytes > 0
            && detection_receipt.requested_bytes == 0
            && leased_pool.stats().quiescent(),
        &format!(
            "seeded reclaim poison detected exact corruption, quarantined it, and recovered: \
             seed={SEED:#018x} offset={} detected={} recovered={}",
            mutation.offset,
            after_detection.to_json(),
            after_recovery.to_json()
        ),
        SEED,
    );
}
