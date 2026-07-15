//! G4 seeded tile-fault battery.
//!
//! Each retained seed selects one logical tile and one numbered touch. The
//! typed failure must retain that provenance, drain siblings, reclaim tile
//! arenas, and leave the same pool reusable. This is a bounded seeded storm,
//! not exhaustive schedule model checking.

use core::ops::ControlFlow;

use fs_alloc::Site;
use fs_exec::{
    Cancelled, Cx, PoolConfig, RunError, TILE_FAULT_PLAN_VERSION, TileFailure, TileFaultPlan,
    TileKernel, TilePlan, TilePool,
};
use fs_substrate::affinity::CcdTopology;

const TILES: u64 = 41;
const TOUCHES_PER_TILE: u32 = 3;
const SEEDS: [u64; 16] = [
    0xF404_0000,
    0xF404_0001,
    0xF404_0002,
    0xF404_0003,
    0xF404_0004,
    0xF404_0005,
    0xF404_0006,
    0xF404_0007,
    0xF404_0008,
    0xF404_0009,
    0xF404_000a,
    0xF404_000b,
    0xF404_000c,
    0xF404_000d,
    0xF404_000e,
    0xF404_000f,
];

struct TouchKernel {
    fault: Option<TileFaultPlan>,
}

impl TileKernel for TouchKernel {
    type Out = u64;

    fn tiles(&self) -> TilePlan {
        TilePlan::new("g4/tile-fault-storm", TILES)
    }

    fn run(&self, tile: u64, cx: &Cx<'_>) -> ControlFlow<Cancelled, u64> {
        for touch in 1..=TOUCHES_PER_TILE {
            if cx.checkpoint().is_err() {
                return ControlFlow::Break(Cancelled);
            }
            let bytes = cx
                .arena()
                .alloc_slice_fill(
                    Site::named("g4/tile-fault-touch"),
                    64,
                    (tile as u8) ^ (touch as u8),
                )
                .expect("bounded G4 arena traffic");
            assert_eq!(bytes.len(), 64);

            if let Some(failure) = self.fault.and_then(|plan| plan.failure_at(tile, touch)) {
                return ControlFlow::Break(cx.refuse(failure));
            }
        }
        ControlFlow::Continue(tile)
    }
}

#[test]
fn g4_seeded_faults_are_structured_drained_and_replayable() {
    let pool = TilePool::new(PoolConfig::new(4, CcdTopology::APPLE_M_CLASS, 0xF404_F00D));

    for seed in SEEDS {
        let plan = TileFaultPlan::seeded(seed, TILES, TOUCHES_PER_TILE).expect("valid plan");
        let result = pool.run(&TouchKernel { fault: Some(plan) });
        match result {
            Err(RunError::TileFailed {
                kernel,
                tile,
                failure: TileFailure::InjectedFault { plan_seed, touch },
                ..
            }) => {
                assert_eq!(kernel, "g4/tile-fault-storm");
                assert_eq!(tile, plan.tile());
                assert_eq!(plan_seed, plan.seed());
                assert_eq!(touch, plan.touch());
            }
            other => panic!("seed {seed:#018x}: expected typed tile failure, got {other:?}"),
        }

        assert!(
            pool.arena_pool().stats().quiescent(),
            "seed {seed:#018x}: fault drain leaked an arena"
        );
        println!(
            "{{\"suite\":\"fs-exec/fault-storm\",\"plan_version\":{},\"seed\":\"{:#018x}\",\"tile\":{},\"touch\":{},\"verdict\":\"pass\"}}",
            TILE_FAULT_PLAN_VERSION,
            plan.seed(),
            plan.tile(),
            plan.touch(),
        );

        let healthy = pool
            .run(&TouchKernel { fault: None })
            .expect("pool remains reusable after injected fault");
        assert_eq!(healthy, TILES * (TILES - 1) / 2);
        assert!(pool.arena_pool().stats().quiescent());
    }
}
