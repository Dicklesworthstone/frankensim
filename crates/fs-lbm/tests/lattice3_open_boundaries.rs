//! G0/G1/G3/G4/G5: native open velocity, force frames, exchange and publication.

use fs_exec::{CancelGate, Cx, KernelRunner, RunError, RunReport, TileKernel, TilePool};
use fs_lbm::d3q19::{BoundarySpec3, CollisionModel3, Face3, FaceBoundary3};
use fs_lbm::d3q27::{CollisionConfig27, CollisionModel27, D3q27, Equilibrium27, Rates27};
use fs_lbm::lattice3::open_boundaries::OpenBoundaryGrid3;
use fs_lbm::lattice3::{D3q19, EquilibriumLattice3, Grid3, GridBudget3, GridError3, Lattice3};
use std::sync::atomic::{AtomicUsize, Ordering};

const SEED: u64 = 0x0B03_D327;

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            fs_exec::StreamKey {
                seed: SEED,
                kernel_id: 227,
                tile: 0,
                iteration: 0,
            },
            fs_exec::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        );
        f(&cx)
    })
}

fn budget() -> GridBudget3 {
    GridBudget3 {
        maximum_cells: 512,
        maximum_mach: 0.2,
    }
}

fn config(model: CollisionModel27) -> CollisionConfig27 {
    CollisionConfig27 {
        model,
        maximum_mach: 0.2,
    }
}

fn bgk() -> CollisionConfig27 {
    config(CollisionModel27::Bgk {
        tau: 0.9,
        equilibrium: Equilibrium27::Factorized,
    })
}

fn boundaries(axis: usize, velocity: [f64; 3], density: f64) -> BoundarySpec3 {
    let mut faces = [FaceBoundary3::Periodic; 6];
    faces[2 * axis] = FaceBoundary3::Velocity { velocity };
    faces[2 * axis + 1] = FaceBoundary3::Pressure { density };
    BoundarySpec3::new(faces)
}

fn address(p: [usize; 3], d: [usize; 3]) -> usize {
    (p[2] * d[1] + p[1]) * d[0] + p[0]
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}"
    );
}

fn totals<L: EquilibriumLattice3<Q>, const Q: usize>(
    grid: &OpenBoundaryGrid3<L, Q>,
) -> (f64, [f64; 3]) {
    let mut mass = 0.0;
    let mut momentum = [0.0; 3];
    for cell in 0..grid.populations().len() {
        let (rho, velocity) = grid.macroscopic(cell).unwrap();
        mass += rho;
        for a in 0..3 {
            momentum[a] += rho * velocity[a];
        }
    }
    (mass, momentum)
}

fn nonequilibrium_stress<L: EquilibriumLattice3<Q>, const Q: usize>(
    grid: &OpenBoundaryGrid3<L, Q>,
    cell: usize,
) -> [[f64; 3]; 3] {
    let (rho, u) = grid.macroscopic(cell).unwrap();
    let force = grid.force_density();
    core::array::from_fn(|a| {
        core::array::from_fn(|b| {
            grid.populations()[cell]
                .iter()
                .zip(L::E)
                .map(|(f, e)| f * f64::from(e[a] * e[b]))
                .sum::<f64>()
                + 0.5 * (u[a] * force[b] + force[a] * u[b])
                - rho * (u[a] * u[b] + if a == b { L::CS2 } else { 0.0 })
        })
    })
}

fn exercise<L: EquilibriumLattice3<Q>, const Q: usize>(collision: L::Collision) {
    let dims = [4, 5, 6];
    let rho = 1.3;
    let velocity = [0.021, -0.017, 0.013];
    with_cx(|cx| {
        for axis in 0..3 {
            let mut periodic = [true; 3];
            periodic[axis] = false;
            for force in [[0.0; 3], [0.00013, -0.00009, 0.00007]] {
                let original = Grid3::<L, Q>::try_new_equilibrium(
                    dims,
                    collision,
                    force,
                    periodic,
                    budget(),
                    |_| false,
                    |_| (rho, velocity),
                )
                .unwrap();
                let mut flow =
                    OpenBoundaryGrid3::try_new(original, cx, boundaries(axis, velocity, rho))
                        .unwrap();
                let (mass_before, momentum_before) = totals(&flow);
                let receipt = flow.step_serial(cx, None, [0.5, -0.75, 1.25]).unwrap();
                let (mass_after, momentum_after) = totals(&flow);
                near(mass_after - mass_before, receipt.mass_change, 3e-12);
                for a in 0..3 {
                    near(
                        momentum_after[a] - momentum_before[a],
                        flow.populations().len() as f64 * force[a] + receipt.impulse[a],
                        3e-12,
                    );
                }
                let area = dims.iter().product::<usize>() / dims[axis];
                assert_eq!(receipt.reconstructed_cells, 2 * area);
                assert_eq!(
                    receipt.open_links,
                    2 * area * (0..Q).filter(|&q| L::E[q][axis] > 0).count()
                );
                // All Cartesian components are physical velocities, including
                // tangential components and the accepted half body-force shift.
                for cell in 0..flow.populations().len() {
                    let p = [
                        cell % dims[0],
                        cell / dims[0] % dims[1],
                        cell / (dims[0] * dims[1]),
                    ];
                    let (density, u) = flow.macroscopic(cell).unwrap();
                    if p[axis] == 0 {
                        for a in 0..3 {
                            near(u[a], velocity[a], 5e-14);
                        }
                    }
                    if p[axis] + 1 == dims[axis] {
                        near(density, rho, 5e-14);
                    }
                }
                let mut p = [1; 3];
                p[axis] = 0;
                let inlet = address(p, dims);
                p[axis] = 1;
                let donor = address(p, dims);
                let inlet_stress = nonequilibrium_stress(&flow, inlet);
                let donor_stress = nonequilibrium_stress(&flow, donor);
                for a in 0..3 {
                    for b in 0..3 {
                        near(inlet_stress[a][b], donor_stress[a][b], 8e-14);
                    }
                }
                if force == [0.0; 3] {
                    for _ in 0..5 {
                        flow.step_serial(cx, None, [0.0; 3]).unwrap();
                    }
                    for cell in 0..flow.populations().len() {
                        let (density, u) = flow.macroscopic(cell).unwrap();
                        near(density, rho, 5e-13);
                        for a in 0..3 {
                            near(u[a], velocity[a], 5e-13);
                        }
                    }
                }
            }
        }
    });
}

#[test]
fn native_collision_families_keep_physical_velocity_stress_and_reservoir_balance() {
    for collision in [
        CollisionModel3::Bgk { tau: 0.9 },
        CollisionModel3::Trt {
            tau: 0.9,
            magic: 3.0 / 16.0,
        },
        CollisionModel3::CentralMoment {
            second_order_rate: 1.2,
            higher_order_rate: 1.6,
        },
        CollisionModel3::CentralMomentBulk {
            shear_rate: 1.2,
            bulk_rate: 0.9,
            higher_order_rate: 1.6,
        },
        CollisionModel3::ReducedCumulant {
            second_order_rate: 1.2,
            third_order_rate: 1.1,
            fourth_order_rate: 0.7,
        },
    ] {
        exercise::<D3q19, 19>(collision);
    }
    for equilibrium in [Equilibrium27::Quadratic, Equilibrium27::Factorized] {
        exercise::<D3q27, 27>(config(CollisionModel27::Bgk {
            tau: 0.9,
            equilibrium,
        }));
        exercise::<D3q27, 27>(config(CollisionModel27::Trt {
            tau: 0.9,
            magic: 3.0 / 16.0,
            equilibrium,
        }));
    }
    for model in [
        CollisionModel27::CentralMoment {
            rates: Rates27 { omega: [1.2; 10] },
        },
        CollisionModel27::Cumulant {
            rates: Rates27 { omega: [1.2; 10] },
        },
    ] {
        exercise::<D3q27, 27>(config(model));
    }
}

fn fixture(cx: &Cx<'_>) -> OpenBoundaryGrid3<D3q27, 27> {
    let grid = Grid3::try_new_equilibrium(
        [8, 8, 5],
        bgk(),
        [0.0001, -0.0002, 0.0003],
        [true, true, false],
        budget(),
        |_| false,
        |_| (1.2, [0.0; 3]),
    )
    .unwrap();
    OpenBoundaryGrid3::try_new(grid, cx, boundaries(2, [0.02, -0.01, 0.03], 1.2)).unwrap()
}

fn bits<L: EquilibriumLattice3<Q>, const Q: usize>(grid: &OpenBoundaryGrid3<L, Q>) -> Vec<u64> {
    grid.populations()
        .iter()
        .flatten()
        .map(|v| v.to_bits())
        .collect()
}

#[test]
fn pooled_open_step_and_receipt_match_serial_and_cancel_before_publication() {
    struct CancelAfterStream {
        pool: TilePool,
        calls: AtomicUsize,
    }
    impl KernelRunner for CancelAfterStream {
        fn workers(&self) -> usize {
            self.pool.workers()
        }
        fn run_with_gate<K: TileKernel>(
            &self,
            kernel: &K,
            gate: &CancelGate,
        ) -> (Result<K::Out, RunError>, RunReport) {
            let result = self.pool.run_with_gate(kernel, gate);
            if self.calls.fetch_add(1, Ordering::AcqRel) == 1 {
                gate.request();
            }
            result
        }
    }
    with_cx(|cx| {
        let mut serial = fixture(cx);
        let targets = boundaries(2, [0.03, -0.02, 0.015], 1.202);
        assert_ne!(serial.boundaries(), targets);
        let expected = serial
            .step_serial_with_open_boundaries(cx, targets, None, [1.25, -0.5, 2.0])
            .unwrap();
        assert_eq!(serial.boundaries(), targets);
        for workers in [1, 2, 7] {
            let mut parallel = fixture(cx);
            let actual = parallel
                .step_pooled_with_open_boundaries(
                    &TilePool::for_host(workers, SEED),
                    &CancelGate::new_clock_free(),
                    targets,
                    None,
                    [1.25, -0.5, 2.0],
                )
                .unwrap();
            assert_eq!(bits(&parallel), bits(&serial));
            assert_eq!(parallel.boundaries(), targets);
            assert_eq!(actual, expected);
        }
        let mut cancelled = fixture(cx);
        let before = bits(&cancelled);
        let previous_targets = cancelled.boundaries();
        let runner = CancelAfterStream {
            pool: TilePool::for_host(2, SEED),
            calls: AtomicUsize::new(0),
        };
        assert!(matches!(
            cancelled.step_pooled_with_open_boundaries(
                &runner,
                &CancelGate::new_clock_free(),
                targets,
                None,
                [1.25, -0.5, 2.0],
            ),
            Err(GridError3::Cancelled)
        ));
        assert_eq!(bits(&cancelled), before);
        assert_eq!(cancelled.completed_steps(), 0);
        assert_eq!(cancelled.boundaries(), previous_targets);
        let retry = cancelled
            .step_serial_with_open_boundaries(cx, targets, None, [1.25, -0.5, 2.0])
            .unwrap();
        assert_eq!(bits(&cancelled), bits(&serial));
        assert_eq!(cancelled.boundaries(), targets);
        assert_eq!(retry, expected);
        let mut shifted = fixture(cx);
        let other_origin = [2.25, 1.5, -1.0];
        let actual = shifted
            .step_serial_with_open_boundaries(cx, targets, None, other_origin)
            .unwrap();
        let shift = [1.0, 2.0, -3.0];
        for a in 0..3 {
            let (b, c) = ((a + 1) % 3, (a + 2) % 3);
            near(
                actual.angular_impulse[a],
                expected.angular_impulse[a] - shift[b] * expected.impulse[c]
                    + shift[c] * expected.impulse[b],
                3e-13,
            );
        }
    });
}

#[test]
fn wall_open_seams_keep_original_wall_and_nonopen_populations() {
    with_cx(|cx| {
        let make = || {
            Grid3::<D3q27, 27>::try_new_equilibrium(
                [5, 3, 5],
                bgk(),
                [0.0; 3],
                [false, true, false],
                budget(),
                |_| false,
                |_| (1.1, [0.005, -0.003, 0.007]),
            )
            .unwrap()
        };
        let mut closed = make();
        let mut faces = [FaceBoundary3::Periodic; 6];
        faces[0] = FaceBoundary3::stationary_wall();
        faces[1] = FaceBoundary3::stationary_wall();
        faces[4] = FaceBoundary3::Velocity {
            velocity: [0.0, 0.01, 0.025],
        };
        faces[5] = FaceBoundary3::Pressure { density: 1.1 };
        let mut open = OpenBoundaryGrid3::try_new(make(), cx, BoundarySpec3::new(faces)).unwrap();
        closed.step_serial(cx).unwrap();
        open.step_serial(cx, None, [0.0; 3]).unwrap();
        for x in [0, 4] {
            let cell = address([x, 1, 0], [5, 3, 5]);
            for (q, e) in D3q27::E.iter().enumerate() {
                let crosses_wall = (x == 0 && e[0] > 0) || (x == 4 && e[0] < 0);
                if e[2] <= 0 || crosses_wall {
                    assert_eq!(
                        open.populations()[cell][q].to_bits(),
                        closed.populations()[cell][q].to_bits()
                    );
                }
            }
        }
        let (_, actual) = open.macroscopic(address([2, 1, 0], [5, 3, 5])).unwrap();
        for (observed, requested) in actual.into_iter().zip([0.0, 0.01, 0.025]) {
            near(observed, requested, 5e-14);
        }
    });
}

#[test]
fn invalid_topology_targets_and_late_cell_mach_refusal_are_atomic() {
    with_cx(|cx| {
        let make = |dims, solid: Option<[usize; 3]>| {
            Grid3::<D3q27, 27>::try_new_equilibrium(
                dims,
                bgk(),
                [0.0; 3],
                [true, true, false],
                budget(),
                |p| Some(p) == solid,
                |_| (1.0, [0.0; 3]),
            )
            .unwrap()
        };
        let spec = boundaries(2, [0.0, 0.0, 0.03], 1.0);
        assert!(OpenBoundaryGrid3::try_new(make([3, 3, 2], None), cx, spec).is_err());
        assert!(OpenBoundaryGrid3::try_new(make([3, 3, 5], Some([1, 1, 1])), cx, spec).is_err());
        assert!(
            OpenBoundaryGrid3::try_new(make([3, 3, 5], None), cx, BoundarySpec3::periodic())
                .is_err()
        );
        let mut grid = OpenBoundaryGrid3::try_new(make([8, 8, 5], None), cx, spec).unwrap();
        let before = bits(&grid);
        let gate = CancelGate::new_clock_free();
        for target in [[f64::NAN, 0.0, 0.0], [0.09; 3]] {
            assert!(
                grid.rebind_open_boundaries(boundaries(2, target, 1.0), &gate)
                    .is_err()
            );
        }
        assert!(
            grid.rebind_open_boundaries(BoundarySpec3::periodic(), &gate)
                .is_err()
        );
        assert_eq!(grid.boundaries(), spec);
        assert_eq!(bits(&grid), before);
        let mut both = [FaceBoundary3::Periodic; 6];
        both[4] = FaceBoundary3::Velocity {
            velocity: [0.0; 3],
        };
        both[5] = both[4];
        let mut late =
            OpenBoundaryGrid3::try_new(make([8, 8, 5], None), cx, BoundarySpec3::new(both))
                .unwrap();
        let mut laws = vec![bgk(); 320];
        laws[319].maximum_mach = 0.01; // Rest collision passes; the final face target does not.
        let previous_targets = late.boundaries();
        both[4] = FaceBoundary3::Velocity { velocity: [0.0, 0.0, 0.03] };
        both[5] = both[4];
        let targets = BoundarySpec3::new(both);
        let late_before = bits(&late);
        assert!(matches!(
            late.step_serial_with_open_boundaries(cx, targets, Some(&laws), [0.0; 3]),
            Err(GridError3::Mach { cell: 319, .. })
        ));
        assert_eq!(bits(&late), late_before);
        assert_eq!(late.completed_steps(), 0);
        assert_eq!(late.boundaries(), previous_targets);
        late.step_serial_with_open_boundaries(cx, targets, None, [0.0; 3]).unwrap();
        assert_eq!(late.completed_steps(), 1);
        assert_eq!(late.boundaries(), targets);
        let cancelled = CancelGate::new_clock_free();
        cancelled.request();
        assert!(matches!(
            grid.rebind_open_boundaries(boundaries(2, [0.0; 3], 1.0), &cancelled),
            Err(GridError3::Cancelled)
        ));
        assert_eq!(grid.boundaries(), spec);
        assert_eq!(bits(&grid), before);
        assert_eq!(grid.completed_steps(), 0);
    });
}
