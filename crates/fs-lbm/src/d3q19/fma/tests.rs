//! Bitwise pins of the FMA capsule against the portable codegen. On an x86-64
//! FMA host the dispatcher runs the capsule and the direct call of
//! `collide_cell3_with_projection` from this (feature-less) test code runs
//! baseline codegen; elsewhere both sides are portable and the pins hold
//! trivially.

use std::sync::Mutex;
use std::sync::atomic::Ordering;

use super::super::{
    BoundaryGrid3, BoundarySpec3, CollisionError3, CollisionModel3, FaceBoundary3, Q3,
    collide_cell3_with_projection, equilibrium3,
};
use super::{PORTABLE_ONLY, collide_cell3};

// The pins share PORTABLE_ONLY; run them one at a time so neither compares
// the portable path with itself by accident.
static SERIAL: Mutex<()> = Mutex::new(());

const MODELS: [CollisionModel3; 5] = [
    CollisionModel3::Bgk { tau: 0.503 },
    CollisionModel3::Bgk { tau: 0.8 },
    CollisionModel3::CentralMoment {
        second_order_rate: 1.23,
        higher_order_rate: 1.71,
    },
    CollisionModel3::ReducedCumulant {
        second_order_rate: 1.988_071_571,
        third_order_rate: 1.0,
        fourth_order_rate: 1.0,
    },
    CollisionModel3::ReducedCumulant {
        second_order_rate: 1.23,
        third_order_rate: 1.35,
        fourth_order_rate: 1.72,
    },
];

/// xorshift64*: deterministic uniform samples in [-1, 1).
struct Samples(u64);

impl Samples {
    fn next(&mut self) -> f64 {
        self.0 ^= self.0 >> 12;
        self.0 ^= self.0 << 25;
        self.0 ^= self.0 >> 27;
        let bits = self.0.wrapping_mul(0x2545_f491_4f6c_dd1d) >> 11;
        2.0 * (bits as f64 / (1u64 << 53) as f64) - 1.0
    }
}

fn bits(result: Result<[f64; Q3], CollisionError3>) -> Result<[u64; Q3], CollisionError3> {
    result.map(|post| post.map(f64::to_bits))
}

#[test]
fn capsule_collision_is_bitwise_the_portable_collision() {
    let _serial = SERIAL
        .lock()
        .unwrap_or_else(std::sync::PoisonError::into_inner);
    let mut samples = Samples(0x9e37_79b9_7f4a_7c15);
    for model in MODELS {
        for state in 0..4000 {
            let rho = 1.0 + 0.1 * samples.next();
            let u = [
                0.12 * samples.next(),
                0.12 * samples.next(),
                0.12 * samples.next(),
            ];
            let eq = equilibrium3(rho, u);
            let populations: [f64; Q3] =
                core::array::from_fn(|q| eq[q] * 0.05f64.mul_add(samples.next(), 1.0));
            let force = if state % 2 == 0 {
                [0.0; 3]
            } else {
                [
                    1e-4 * samples.next(),
                    1e-4 * samples.next(),
                    1e-4 * samples.next(),
                ]
            };
            assert_eq!(
                bits(collide_cell3(populations, model, force)),
                bits(collide_cell3_with_projection(
                    populations,
                    model,
                    force,
                    false
                )),
                "{model:?}, state {state}"
            );
        }
        // Refusals are identical too.
        let mut empty = [0.0; Q3];
        empty[0] = -1.0;
        let mut poisoned = equilibrium3(1.0, [0.01, 0.0, 0.0]);
        poisoned[7] = f64::NAN;
        for populations in [empty, poisoned] {
            let capsule = collide_cell3(populations, model, [0.0; 3]);
            let portable = collide_cell3_with_projection(populations, model, [0.0; 3], false);
            assert!(capsule.is_err(), "{model:?}: {capsule:?}");
            assert_eq!(format!("{capsule:?}"), format!("{portable:?}"), "{model:?}");
        }
    }
}

/// A 16 x 16 x 4 grid with a solid cylinder, run for `steps` steps.
fn run_grid(model: CollisionModel3, open: bool, steps: usize) -> Vec<[u64; Q3]> {
    let periodic = FaceBoundary3::Periodic;
    let (spec, force) = if open {
        let spec = BoundarySpec3::new([
            FaceBoundary3::Velocity {
                velocity: [0.02, 0.0, 0.0],
            },
            FaceBoundary3::Pressure { density: 1.0 },
            periodic,
            periodic,
            periodic,
            periodic,
        ]);
        (spec, [0.0; 3])
    } else {
        (BoundarySpec3::new([periodic; 6]), [2e-4, 3e-5, 0.0])
    };
    let mut grid = BoundaryGrid3::with_collision_model(16, 16, 4, model, force, spec);
    grid.voxelize_sdf(|p| ((p[0] - 7.3).powi(2) + (p[1] - 8.1).powi(2)).sqrt() - 2.6);
    for _ in 0..steps {
        grid.step();
    }
    let mut state = Vec::new();
    for z in 0..4 {
        for y in 0..16 {
            for x in 0..16 {
                if !grid.is_solid(x, y, z) {
                    state.push(grid.populations(x, y, z).map(f64::to_bits));
                }
            }
        }
    }
    state
}

#[test]
fn capsule_grid_run_is_bitwise_the_portable_grid_run() {
    let _serial = SERIAL
        .lock()
        .unwrap_or_else(std::sync::PoisonError::into_inner);
    for model in MODELS {
        for open in [false, true] {
            PORTABLE_ONLY.store(true, Ordering::SeqCst);
            let portable = run_grid(model, open, 60);
            PORTABLE_ONLY.store(false, Ordering::SeqCst);
            let dispatched = run_grid(model, open, 60);
            assert!(
                portable == dispatched,
                "{model:?}, open {open}: the grid state differs between the paths"
            );
        }
    }
}
