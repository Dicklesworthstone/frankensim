//! Real TilePool → fs-rand v2 integration: no narrowing, no generation folding,
//! and no dependence on scheduling or the pool's prior work.
use core::ops::ControlFlow;
use fs_exec::{CancelGate, Cancelled, Cx, PoolConfig, RunId, TileKernel, TilePlan, TilePool};
use fs_rand::{StreamKeyV2, StreamV2};
use fs_substrate::affinity::CcdTopology;

const SEED: u64 = 0xdef0_1234_5678_9abc;
const TILES: u64 = 31;

struct RandomTiles;

impl TileKernel for RandomTiles {
    type Out = Vec<(u64, StreamKeyV2, u64)>;

    fn tiles(&self) -> TilePlan {
        TilePlan::new("rand/real-pool-v2", TILES)
    }

    fn run(
        &self,
        tile: u64,
        cx: &Cx<'_>,
    ) -> ControlFlow<Cancelled, Self::Out> {
        if cx.checkpoint().is_err() {
            return ControlFlow::Break(Cancelled);
        }
        let exec_key = cx.stream_key();
        let key = StreamKeyV2::from_exec_parts_v2(
            exec_key.seed,
            exec_key.kernel_id,
            exec_key.tile,
            exec_key.iteration,
        );
        let first = key.stream().next_u64();
        ControlFlow::Continue(vec![(tile, key, first)])
    }
}

fn pool(workers: usize) -> TilePool {
    TilePool::new(PoolConfig::new(workers, CcdTopology::APPLE_M_CLASS, SEED))
}

fn declared(pool: &TilePool, generation: u64) -> Vec<(u64, StreamKeyV2, u64)> {
    pool.run_declared(&RandomTiles, &CancelGate::new(), RunId(generation))
        .0
        .expect("real pool run")
}

#[test]
fn real_pool_full_width_bridge_is_lossless_and_schedule_independent() {
    let reference_pool = pool(1);
    let kernel_id = TilePlan::new("rand/real-pool-v2", TILES).kernel_id();
    // A real hashed kernel identity is wider than the legacy u32 slot.
    assert!(kernel_id > u64::from(u32::MAX));
    for generation in [0, 1, 1_u64 << 32, u64::MAX] {
        let expected = declared(&reference_pool, generation);
        assert_eq!(expected.len(), TILES as usize);
        for (tile, (ordinal, key, draw)) in expected.iter().enumerate() {
            assert_eq!(*ordinal, tile as u64);
            assert_eq!(key.seed, SEED);
            assert_eq!(key.kernel, kernel_id);
            assert_eq!(key.tile, tile as u64);
            assert_eq!(key.iteration, generation);
            assert_eq!(*draw, StreamV2::at(*key, 0)[0]);
            let mut stream = key.stream();
            assert_eq!(stream.next_u64(), *draw);
            let retained = stream.checkpoint().to_canonical_le_bytes();
            let mut resumed = StreamV2::resume_retained(&retained).expect("valid v2 replay");
            assert_eq!(resumed.next_u64(), stream.next_u64());
        }

        for workers in [2, 4, 7] {
            let reused = pool(workers);
            assert_eq!(declared(&reused, generation), expected);
            // An unrelated generation cannot perturb this run's identity.
            let _ = declared(&reused, generation.wrapping_add(1));
            assert_eq!(declared(&reused, generation), expected);
        }
    }
}

#[test]
fn generations_cannot_alias_at_width_boundaries() {
    let pool = pool(3);
    let mut keys = std::collections::HashSet::new();
    for generation in [0, 1, 1_u64 << 32, (1_u64 << 32) | 5, u64::MAX] {
        let rows = declared(&pool, generation);
        for (_, key, draw) in rows {
            assert!(keys.insert(key), "distinct logical work must have distinct identity");
            assert_eq!(draw, StreamV2::at(key, 0)[0]);
        }
    }
    assert_eq!(keys.len(), 5 * TILES as usize);
}
