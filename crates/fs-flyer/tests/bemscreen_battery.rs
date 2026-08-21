//! E4.8 battery (bead wf-root-guzez.5.20): the screening preset's
//! deviation vs full-resolution fs-bem is BOUNDED on the pinned
//! sphere fixture (analytic 1.5 surface-speed ratio referee); the
//! panel budget refuses at the budget AND one under what the preset
//! needs; the FULL cache key is load-bearing per FIELD (six per-field
//! miss oracles); the stale-key hostile twin refuses; commits cancel
//! stale jobs; caps at cap AND cap+1; determinism golden.
//! Repro: cargo test -p fs-flyer --test bemscreen_battery

use fs_flyer::bemscreen::{
    DesignCommitCache, InterferenceCacheKey, InterferenceResult, MAX_CACHE_ENTRIES,
    SCREEN_SOLVER_VERSION, schematic_preview, screening_solve,
};

fn jlog(case: &str, payload: &str) {
    println!("{{\"suite\":\"fs-flyer-bemscreen\",\"case\":\"{case}\",{payload}}}");
}

#[test]
fn preset_deviation_bounded_vs_full_resolution_on_the_pinned_sphere() {
    // The sphere's analytic max surface-speed ratio is exactly 1.5.
    let coarse = screening_solve(1, 320, 13.0).unwrap();
    let full = screening_solve(2, 320, 13.0).unwrap();
    assert_eq!(coarse.n_panels, 80);
    assert_eq!(full.n_panels, 320);
    let coarse_dev = (coarse.max_speed_ratio - 1.5).abs();
    let full_dev = (full.max_speed_ratio - 1.5).abs();
    // The honest resolution ladder: full is closer to the referee…
    assert!(
        full_dev < coarse_dev,
        "full {full_dev} must beat coarse {coarse_dev}"
    );
    // …and the coarse preset's deviation is BOUNDED (measured
    // 2026-08-21; ~5x headroom — the DONE-WHEN bound).
    assert!(coarse_dev < 0.25, "coarse deviation {coarse_dev}");
    assert!(full_dev < 0.10, "full deviation {full_dev}");
    // Determinism: bit-identical screening twice.
    let again = screening_solve(1, 320, 13.0).unwrap();
    assert_eq!(
        again.result_digest, coarse.result_digest,
        "bit-identical twice"
    );
    jlog(
        "deviation",
        &format!(
            "\"coarse_dev\":{coarse_dev},\"full_dev\":{full_dev},\"digest\":\"{}\"",
            coarse.result_digest
        ),
    );
}

#[test]
fn panel_budget_refuses_on_exhaustion() {
    // AT the exact need admits; one panel fewer refuses.
    assert!(screening_solve(1, 80, 13.0).is_ok(), "AT the need");
    let err = screening_solve(1, 79, 13.0).unwrap_err();
    assert_eq!(err.code, "bem-screen-budget-exhausted");
    // The hard tier cap also binds: subdiv 3 (1280) admits at the
    // tier cap; subdiv 4 (5120) refuses regardless of the request.
    let err = screening_solve(4, 100_000, 13.0).unwrap_err();
    assert_eq!(err.code, "bem-screen-budget-exhausted");
    jlog("budget", &format!("\"code\":\"{}\"", err.code));
}

fn base_key() -> InterferenceCacheKey {
    InterferenceCacheKey {
        geometry_digest: "geom-v1".into(),
        operating_grid_digest: "grid-v1".into(),
        panel_preset: 1,
        ground_mode: "free-air",
        solver_version: SCREEN_SOLVER_VERSION,
        coefficient_convention: "wf-coeff-v1",
    }
}

#[test]
fn full_cache_key_is_load_bearing_per_field() {
    let mut cache = DesignCommitCache::default();
    let job = cache.commit();
    let key = base_key();
    let result = InterferenceResult {
        key_digest: key.digest(),
        interference_factor: 0.83,
    };
    cache.deliver(job, &key, result.clone()).unwrap();
    // Exact key: hit, verbatim.
    assert_eq!(cache.get(&key).unwrap().interference_factor, 0.83);
    // EVERY field change is a miss (six per-field oracles — reuse
    // under a drifted key is the failure the full key exists to stop).
    let variants: Vec<InterferenceCacheKey> = vec![
        InterferenceCacheKey {
            geometry_digest: "geom-v2".into(),
            ..base_key()
        },
        InterferenceCacheKey {
            operating_grid_digest: "grid-v2".into(),
            ..base_key()
        },
        InterferenceCacheKey {
            panel_preset: 2,
            ..base_key()
        },
        InterferenceCacheKey {
            ground_mode: "image-flat",
            ..base_key()
        },
        InterferenceCacheKey {
            solver_version: "fs-bem-exterior-gmres-v2",
            ..base_key()
        },
        InterferenceCacheKey {
            coefficient_convention: "wf-coeff-v2",
            ..base_key()
        },
    ];
    for (i, k) in variants.iter().enumerate() {
        assert!(cache.get(k).is_none(), "field {i} must be load-bearing");
    }
    jlog("full-key", "\"per_field_misses\":6");
}

#[test]
fn stale_key_reuse_and_stale_jobs_refuse() {
    let mut cache = DesignCommitCache::default();
    let job = cache.commit();
    // HOSTILE TWIN: a result derived under the free-air key filed
    // under the image-flat key — refused, never adapted.
    let derived_under = base_key();
    let filed_under = InterferenceCacheKey {
        ground_mode: "image-flat",
        ..base_key()
    };
    let stale = InterferenceResult {
        key_digest: derived_under.digest(),
        interference_factor: 0.83,
    };
    let err = cache.deliver(job, &filed_under, stale).unwrap_err();
    assert_eq!(err.code, "bem-cache-key-mismatch");
    // Commit cancels stale jobs: a job minted before the new commit
    // refuses to deliver.
    let old_job = cache.commit();
    let _newer = cache.commit();
    let ok_result = InterferenceResult {
        key_digest: base_key().digest(),
        interference_factor: 0.83,
    };
    let err = cache.deliver(old_job, &base_key(), ok_result).unwrap_err();
    assert_eq!(err.code, "bem-job-stale");
    jlog(
        "hostile-twins",
        "\"stale_key\":\"refused\",\"stale_job\":\"refused\"",
    );
}

#[test]
fn cache_caps_and_schematic_preview_is_numberless() {
    let mut cache = DesignCommitCache::default();
    // Fill to the cap (fresh job per delivery epoch discipline: one
    // commit, many deliveries at that epoch).
    let job = cache.commit();
    for i in 0..MAX_CACHE_ENTRIES {
        let key = InterferenceCacheKey {
            geometry_digest: format!("geom-{i}"),
            ..base_key()
        };
        let r = InterferenceResult {
            key_digest: key.digest(),
            interference_factor: 0.5,
        };
        cache.deliver(job, &key, r).unwrap();
    }
    // Cap+1 refuses.
    let key = InterferenceCacheKey {
        geometry_digest: "geom-overflow".into(),
        ..base_key()
    };
    let r = InterferenceResult {
        key_digest: key.digest(),
        interference_factor: 0.5,
    };
    let err = cache.deliver(job, &key, r).unwrap_err();
    assert_eq!(err.code, "bem-cache-full");
    // The slider preview carries proportions and a label — the type
    // has no coefficient field to leak solver numbers through.
    let p = schematic_preview(0.154, 0.0);
    assert_eq!(p.label, "schematic preview — commit to derive");
    assert_eq!(p.gap_over_span, 0.154);
    jlog("caps", &format!("\"max_entries\":{MAX_CACHE_ENTRIES}"));
}
