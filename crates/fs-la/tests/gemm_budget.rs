//! wf9.15 battery: the pool-GEMM memory envelope. Preflight refusal
//! before any allocation or C mutation, ledgered requested/limit
//! bytes, cancellation during init, and success equivalence between
//! budgeted and unbudgeted paths.

use fs_exec::{CancelGate, PoolConfig, TilePool};
use fs_la::gemm_f64_parallel_with_pool_budgeted;
use fs_la::{GemmMemoryEnvelope, GemmRunError, gemm_f64_parallel_with_pool};

fn fixtures(m: usize, n: usize, k: usize) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let a: Vec<f64> = (0..m * k).map(|i| ((i as f64) * 0.7).sin()).collect();
    let b: Vec<f64> = (0..k * n).map(|i| ((i as f64) * 1.3).cos()).collect();
    let c: Vec<f64> = (0..m * n).map(|i| (i as f64) * 0.01 - 3.0).collect();
    (a, b, c)
}

fn pool() -> TilePool {
    TilePool::new(PoolConfig::for_host(4, 0xB0D6))
}

#[test]
fn tiny_envelope_refuses_before_touching_c() {
    let (m, n, k) = (96usize, 64, 48);
    let (a, b, c0) = fixtures(m, n, k);
    let mut c = c0.clone();
    let pool = pool();
    let gate = CancelGate::new();
    let err = gemm_f64_parallel_with_pool_budgeted(
        m,
        n,
        k,
        1.25,
        &a,
        &b,
        0.5,
        &mut c,
        &pool,
        32,
        64,
        &gate,
        fs_exec::RunId::default(),
        GemmMemoryEnvelope { limit_bytes: 64 },
    )
    .expect_err("64 bytes cannot admit this plan");
    match err {
        GemmRunError::MemoryRefused {
            what,
            requested_bytes,
            limit_bytes,
            report,
        } => {
            assert_eq!(what, "preflight-envelope");
            assert_eq!(limit_bytes, 64);
            assert!(requested_bytes > 64);
            assert_eq!(report.completed_tiles, 0);
            assert!(report.pool_runs.is_empty(), "nothing was dispatched");
            assert_eq!(report.memory.requested_bytes, requested_bytes);
            // The plan decomposition is ledgered and self-consistent.
            let sum = report.memory.staging_bytes
                + report.memory.b_pack_bytes
                + report.memory.band_metadata_bytes
                + report.memory.pool_run_bytes
                + report.memory.arena_bytes;
            assert_eq!(report.memory.requested_bytes, sum);
            assert_eq!(report.memory.staging_bytes, (m * n * 8) as u128);
        }
        other => panic!("expected MemoryRefused, got {other:?}"),
    }
    assert!(
        c.iter().zip(&c0).all(|(x, y)| x.to_bits() == y.to_bits()),
        "refusal must leave C bitwise unchanged"
    );
}

#[test]
fn roomy_envelope_is_bitwise_equivalent_to_unbudgeted() {
    let (m, n, k) = (96usize, 64, 48);
    let (a, b, c0) = fixtures(m, n, k);
    let pool = pool();
    let gate = CancelGate::new();
    let mut c_unbudgeted = c0.clone();
    let unbudgeted = gemm_f64_parallel_with_pool(
        m,
        n,
        k,
        1.25,
        &a,
        &b,
        0.5,
        &mut c_unbudgeted,
        &pool,
        32,
        64,
        &gate,
    )
    .expect("unbudgeted runs");
    let mut c_budgeted = c0.clone();
    let budgeted = gemm_f64_parallel_with_pool_budgeted(
        m,
        n,
        k,
        1.25,
        &a,
        &b,
        0.5,
        &mut c_budgeted,
        &pool,
        32,
        64,
        &gate,
        fs_exec::RunId::default(),
        GemmMemoryEnvelope {
            limit_bytes: 64 << 20,
        },
    )
    .expect("roomy envelope admits");
    assert!(
        c_unbudgeted
            .iter()
            .zip(&c_budgeted)
            .all(|(x, y)| x.to_bits() == y.to_bits()),
        "envelope admission must not change bits"
    );
    assert_eq!(budgeted.completed_tiles, budgeted.total_tiles);
    assert_eq!(budgeted.completed_tiles, unbudgeted.completed_tiles);
    // The admitted plan is ledgered on BOTH paths (unbudgeted =
    // unbounded envelope), with the caller's limit recorded verbatim.
    assert_eq!(budgeted.memory.limit_bytes, 64 << 20);
    assert_eq!(unbudgeted.memory.limit_bytes, u64::MAX);
    assert_eq!(
        budgeted.memory.requested_bytes,
        unbudgeted.memory.requested_bytes
    );
    assert!(budgeted.memory.requested_bytes > 0);
}

#[test]
fn cancellation_during_init_leaves_c_unchanged_and_ledgers_the_plan() {
    let (m, n, k) = (96usize, 64, 48);
    let (a, b, c0) = fixtures(m, n, k);
    let mut c = c0.clone();
    let pool = pool();
    let gate = CancelGate::new();
    gate.request(); // pre-cancelled: trips at the first bounded poll
    let err = gemm_f64_parallel_with_pool_budgeted(
        m,
        n,
        k,
        1.25,
        &a,
        &b,
        0.5,
        &mut c,
        &pool,
        32,
        64,
        &gate,
        fs_exec::RunId::default(),
        GemmMemoryEnvelope {
            limit_bytes: 64 << 20,
        },
    )
    .expect_err("pre-requested gate cancels");
    match err {
        GemmRunError::Cancelled(cancelled) => {
            assert_eq!(cancelled.report.completed_tiles, 0);
            // The preflight plan was computed and ledgered before the
            // cancellation was observed.
            assert!(cancelled.report.memory.requested_bytes > 0);
        }
        other => panic!("expected Cancelled, got {other:?}"),
    }
    assert!(
        c.iter().zip(&c0).all(|(x, y)| x.to_bits() == y.to_bits()),
        "cancellation during init must leave C bitwise unchanged"
    );
}

#[test]
fn no_product_charges_only_transactional_staging() {
    let (m, n, k) = (19usize, 11, 7);
    let (a, b, c0) = fixtures(m, n, k);
    let mut c = c0.clone();
    let pool = pool();
    let staging_bytes = (m * n * core::mem::size_of::<f64>()) as u64;
    let report = gemm_f64_parallel_with_pool_budgeted(
        m,
        n,
        k,
        0.0,
        &a,
        &b,
        0.5,
        &mut c,
        &pool,
        32,
        64,
        &CancelGate::new(),
        fs_exec::RunId(9),
        GemmMemoryEnvelope {
            limit_bytes: staging_bytes,
        },
    )
    .expect("alpha-zero work needs only C staging");
    assert_eq!(report.total_tiles, 0);
    assert!(report.pool_runs.is_empty());
    assert_eq!(report.memory.staging_bytes, u128::from(staging_bytes));
    assert_eq!(report.memory.requested_bytes, u128::from(staging_bytes));
    assert_eq!(report.memory.peak_used_bytes, u128::from(staging_bytes));
    assert_eq!(report.memory.refused_bytes, 0);
    assert_eq!(report.memory.b_pack_bytes, 0);
    assert_eq!(report.memory.band_metadata_bytes, 0);
    assert_eq!(report.memory.pool_run_bytes, 0);
    assert_eq!(report.memory.arena_bytes, 0);
    assert_eq!(report.memory.active_arena_workers, 0);
    assert!(
        c.iter()
            .zip(&c0)
            .all(|(actual, before)| actual.to_bits() == (before * 0.5).to_bits())
    );
}

#[test]
fn arena_plan_uses_real_chunk_size_and_active_workers() {
    let (m, n, k) = (16usize, 8, 5);
    let (a, b, mut c) = fixtures(m, n, k);
    let mut config = PoolConfig::for_host(8, 0xA11C);
    config.arena.chunk_bytes = 1 << 20;
    let pool = TilePool::new(config);
    let report = gemm_f64_parallel_with_pool_budgeted(
        m,
        n,
        k,
        1.0,
        &a,
        &b,
        0.0,
        &mut c,
        &pool,
        32,
        64,
        &CancelGate::new(),
        fs_exec::RunId(10),
        GemmMemoryEnvelope {
            limit_bytes: 16 << 20,
        },
    )
    .expect("one M band is admitted");
    let expected = pool
        .arena_pool()
        .reservation_bytes_for_slice::<f64>(
            fs_alloc::Site::named("fs-la/gemm-a-micro-panel"),
            8 * 256,
        )
        .expect("fixed A panel reservation");
    assert_eq!(report.memory.arena_bytes_per_worker, expected as u64);
    assert_eq!(report.memory.active_arena_workers, 1);
    assert_eq!(report.memory.arena_bytes, expected as u128);
    assert_ne!(report.memory.arena_bytes_per_worker, u64::MAX);
    assert_eq!(report.memory.refused_bytes, 0);
    assert!(report.memory.peak_used_bytes <= report.memory.requested_bytes);
}
