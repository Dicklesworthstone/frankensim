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

const SUITE: &str = "fs-exec/fault-storm";
const TILES: u64 = 41;
const TOUCHES_PER_TILE: u32 = 3;
const SEEDED_POOL_EXECUTION_ROOT: u64 = 0xF404_F00D;
const EARLY_POOL_EXECUTION_ROOT: u64 = 0xF404_DA1A;
const EARLY_SEED: u64 = 0xF404_001a;
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

fn verdict(emitter: &mut fs_obs::Emitter, case: &str, pass: bool, detail: &str, seed: u64) {
    let event = emitter.emit(
        if pass {
            fs_obs::Severity::Info
        } else {
            fs_obs::Severity::Error
        },
        fs_obs::EventKind::ConformanceCase {
            suite: SUITE.to_string(),
            case: case.to_string(),
            pass,
            detail: detail.to_string(),
            seed,
        },
        None,
    );
    fs_obs::lint_failure_record(&event).expect("fault-storm verdict must be replayable");
    let line = event.to_jsonl();
    fs_obs::validate_line(&line).expect("fault-storm verdict must use the fs-obs wire schema");
    println!("{line}");
    assert!(pass, "case {case}: {detail}");
}

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
    let pool = TilePool::new(PoolConfig::new(
        4,
        CcdTopology::APPLE_M_CLASS,
        SEEDED_POOL_EXECUTION_ROOT,
    ));
    let mut emitter = fs_obs::Emitter::new(SUITE, "seeded-fault");

    for seed in SEEDS {
        let plan = TileFaultPlan::seeded(seed, TILES, TOUCHES_PER_TILE).expect("valid plan");
        let result = pool.run(&TouchKernel { fault: Some(plan) });
        let completed = match result {
            Err(RunError::TileFailed {
                kernel,
                tile,
                failure:
                    TileFailure::InjectedFault {
                        plan_version,
                        plan_seed,
                        tiles,
                        touches_per_tile,
                        touch,
                    },
                completed,
            }) => {
                assert_eq!(kernel, "g4/tile-fault-storm");
                assert_eq!(tile, plan.tile());
                assert_eq!(plan_version, plan.version());
                assert_eq!(plan_seed, plan.seed());
                assert_eq!(tiles, plan.tiles());
                assert_eq!(touches_per_tile, plan.touches_per_tile());
                assert_eq!(touch, plan.touch());
                completed
            }
            other => panic!("seed {seed:#018x}: expected typed tile failure, got {other:?}"),
        };

        assert!(
            pool.arena_pool().stats().quiescent(),
            "seed {seed:#018x}: fault drain leaked an arena"
        );
        assert!(completed < TILES, "the refusing tile cannot be completed");

        let healthy = pool
            .run(&TouchKernel { fault: None })
            .expect("pool remains reusable after injected fault");
        assert_eq!(healthy, TILES * (TILES - 1) / 2);
        assert!(pool.arena_pool().stats().quiescent());

        verdict(
            &mut emitter,
            "seeded-fault",
            true,
            &format!(
                "plan_version={}, plan_seed={:#018x}, tiles={}, touches_per_tile={}, tile={}, \
                 touch={}, completed_before_drain={}, pool_execution_root={:#018x}",
                plan.version(),
                plan.seed(),
                plan.tiles(),
                plan.touches_per_tile(),
                plan.tile(),
                plan.touch(),
                completed,
                SEEDED_POOL_EXECUTION_ROOT,
            ),
            plan.seed(),
        );
    }
}

#[test]
fn g4_early_fault_stops_single_worker_before_its_next_claim() {
    let pool = TilePool::new(PoolConfig::new(
        1,
        CcdTopology::APPLE_M_CLASS,
        EARLY_POOL_EXECUTION_ROOT,
    ));
    let plan = TileFaultPlan::seeded(EARLY_SEED, TILES, TOUCHES_PER_TILE).expect("valid plan");
    assert_eq!(
        (plan.tile(), plan.touch()),
        (0, 1),
        "golden early-fault plan"
    );

    let completed = match pool.run(&TouchKernel { fault: Some(plan) }) {
        Err(RunError::TileFailed {
            tile: 0,
            failure:
                TileFailure::InjectedFault {
                    plan_version: TILE_FAULT_PLAN_VERSION,
                    plan_seed: EARLY_SEED,
                    tiles: TILES,
                    touches_per_tile: TOUCHES_PER_TILE,
                    touch: 1,
                },
            completed,
            ..
        }) => completed,
        other => panic!("expected the golden early fault, got {other:?}"),
    };
    assert_eq!(
        completed, 0,
        "a one-worker pool must stop before claiming tile 1 after tile 0 refuses"
    );
    assert!(pool.arena_pool().stats().quiescent());
    let mut emitter = fs_obs::Emitter::new(SUITE, "early-drain");
    verdict(
        &mut emitter,
        "early-drain",
        true,
        &format!(
            "plan_version={}, plan_seed={:#018x}, tiles={}, touches_per_tile={}, tile=0, \
             touch=1, completed_before_drain=0, pool_execution_root={:#018x}",
            plan.version(),
            plan.seed(),
            plan.tiles(),
            plan.touches_per_tile(),
            EARLY_POOL_EXECUTION_ROOT,
        ),
        plan.seed(),
    );
}
