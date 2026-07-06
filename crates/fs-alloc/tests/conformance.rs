//! fs-alloc conformance suite (CONTRACT.md: any reimplementation must pass).
//!
//! Gauntlet coverage: G0 (accounting laws, alignment laws, shadow-model
//! equivalence), G4 (the 10^6-cancellation storm, concurrent hammer), and
//! G5 (deterministic accounting reports). Every case prints one JSON-line
//! verdict; randomized cases carry their seed so failures replay from the
//! log alone.

use fs_alloc::{
    ALLOC_ALIGN, AllocError, ArenaConfig, ArenaPool, HUGEPAGE_BYTES, HugepageOutcome,
    HugepagePolicy, ShardedPool, Site, SiteStats,
};

fn verdict(case: &str, pass: bool, detail: &str) {
    println!(
        "{{\"suite\":\"fs-alloc/conformance\",\"case\":\"{case}\",\"verdict\":\"{}\",\
         \"detail\":\"{detail}\"}}",
        if pass { "pass" } else { "fail" }
    );
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
    );
}

#[test]
fn alloc_006_concurrent_arenas_and_pools_stay_leak_free() {
    let pool = ArenaPool::new(small_config());
    let shaped: ShardedPool<Vec<u64>> = ShardedPool::new(4);
    std::thread::scope(|s| {
        let pool = &pool;
        let shaped = &shaped;
        for t in 0..8usize {
            s.spawn(move || {
                let mut rng = Lcg(0xC0C0 + t as u64);
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
            "8 threads x 400 scopes leak nothing: arena={} pool={}",
            ps.to_json(),
            ss.to_json()
        ),
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
    );
}
