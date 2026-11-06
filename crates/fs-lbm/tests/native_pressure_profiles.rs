//! Spatial pressure on native D3Q19/D3Q27 reservoirs: hydrostatics, torque,
//! original wall ownership, pooled determinism, and atomic field updates.

use fs_exec::{CancelGate, Cx, TilePool};
use fs_lbm::d3q19::{BoundaryProfileError3, BoundarySpec3, CollisionModel3, Face3, FaceBoundary3};
use fs_lbm::d3q27::{CollisionConfig27, CollisionModel27, D3q27, Equilibrium27, Rates27};
use fs_lbm::lattice3::open_boundaries::OpenBoundaryGrid3;
use fs_lbm::lattice3::{D3q19, EquilibriumLattice3, Grid3, GridBudget3, GridError3, Lattice3};

const SEED: u64 = 0x5052_3350;

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        f(&Cx::new(
            &gate,
            arena,
            fs_exec::StreamKey {
                seed: SEED,
                kernel_id: 327,
                tile: 0,
                iteration: 0,
            },
            fs_exec::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        ))
    })
}

fn budget() -> GridBudget3 {
    GridBudget3 {
        maximum_cells: 1024,
        maximum_mach: 0.2,
    }
}

fn config(model: CollisionModel27) -> CollisionConfig27 {
    CollisionConfig27 {
        model,
        maximum_mach: 0.2,
    }
}

fn reservoirs(axis: usize, density: f64) -> BoundarySpec3 {
    let mut faces = [FaceBoundary3::stationary_wall(); 6];
    faces[2 * axis] = FaceBoundary3::Pressure { density };
    faces[2 * axis + 1] = FaceBoundary3::Pressure { density };
    BoundarySpec3::new(faces)
}

fn coordinates(cell: usize, dims: [usize; 3]) -> [usize; 3] {
    [
        cell % dims[0],
        cell / dims[0] % dims[1],
        cell / (dims[0] * dims[1]),
    ]
}

fn address(p: [usize; 3], dims: [usize; 3]) -> usize {
    (p[2] * dims[1] + p[1]) * dims[0] + p[0]
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}"
    );
}

fn bits<L: EquilibriumLattice3<Q>, const Q: usize>(grid: &OpenBoundaryGrid3<L, Q>) -> Vec<u64> {
    grid.populations()
        .iter()
        .flatten()
        .map(|value| value.to_bits())
        .collect()
}

fn hydrostatic<L: EquilibriumLattice3<Q>, const Q: usize>(collision: L::Collision) {
    with_cx(|cx| {
        let dims = [4, 5, 6];
        let gate = CancelGate::new_clock_free();
        for axis in 0..3 {
            let mut force = [0.0; 3];
            force[(axis + 1) % 3] = 1.7e-5;
            force[(axis + 2) % 3] = -2.3e-5;
            // Exact static balance for a uniform FORCE DENSITY: grad(p)=F,
            // with p=cs²*rho. Both transverse directions vary independently.
            let density = |point: [usize; 3]| {
                1.2 + (0..3)
                    .map(|a| force[a] * (point[a] as f64 + 0.5) / L::CS2)
                    .sum::<f64>()
            };
            let original = Grid3::<L, Q>::try_new_equilibrium(
                dims,
                collision,
                force,
                [false; 3],
                budget(),
                |_| false,
                |point| (density(point), [0.0; 3]),
            )
            .unwrap();
            let mut grid = OpenBoundaryGrid3::try_new(original, cx, reservoirs(axis, 1.2)).unwrap();
            for side in 0..2 {
                let face = Face3::ALL[2 * axis + side];
                let report = grid
                    .set_pressure_profile(face, &gate, |point| Ok(density(point)))
                    .unwrap();
                assert_eq!(report.face, face);
                assert_eq!(
                    report.sampled_cells,
                    dims.iter().product::<usize>() / dims[axis]
                );
                assert!(report.maximum_density_lu > report.minimum_density_lu);
            }
            for _ in 0..64 {
                grid.step_serial(cx, None, [0.0; 3]).unwrap();
            }
            for cell in 0..grid.populations().len() {
                let (rho, velocity) = grid.macroscopic(cell).unwrap();
                near(rho, density(coordinates(cell, dims)), 2e-11);
                for value in velocity {
                    near(value, 0.0, 2e-11);
                }
            }
        }
    });
}

#[test]
fn transverse_pressure_profiles_preserve_three_dimensional_hydrostatics() {
    hydrostatic::<D3q19, 19>(CollisionModel3::Bgk { tau: 0.9 });
    for equilibrium in [Equilibrium27::Quadratic, Equilibrium27::Factorized] {
        hydrostatic::<D3q27, 27>(config(CollisionModel27::Bgk {
            tau: 0.9,
            equilibrium,
        }));
        hydrostatic::<D3q27, 27>(config(CollisionModel27::Trt {
            tau: 0.9,
            magic: 3.0 / 16.0,
            equilibrium,
        }));
    }
}

fn extensive<L: EquilibriumLattice3<Q>, const Q: usize>(
    grid: &OpenBoundaryGrid3<L, Q>,
    origin: [f64; 3],
) -> (f64, [f64; 3], [f64; 3]) {
    let mut mass = 0.0;
    let mut momentum = [0.0; 3];
    let mut angular = [0.0; 3];
    for (cell, populations) in grid.populations().iter().enumerate() {
        let point = coordinates(cell, grid.dimensions());
        let arm: [f64; 3] = core::array::from_fn(|a| point[a] as f64 + 0.5 - origin[a]);
        let mut local = [0.0; 3];
        mass += populations.iter().sum::<f64>();
        for (q, &population) in populations.iter().enumerate() {
            for a in 0..3 {
                local[a] += population * f64::from(L::E[q][a]);
            }
        }
        for a in 0..3 {
            let (b, c) = ((a + 1) % 3, (a + 2) % 3);
            momentum[a] += local[a];
            angular[a] += arm[b] * local[c] - arm[c] * local[b];
        }
    }
    (mass, momentum, angular)
}

fn pressure_patch<L: EquilibriumLattice3<Q>, const Q: usize>(collision: L::Collision) {
    with_cx(|cx| {
        let origin = [4.0; 3];
        for axis in 0..3 {
            let original = Grid3::<L, Q>::try_new_equilibrium(
                [8; 3],
                collision,
                [0.0; 3],
                [false; 3],
                budget(),
                |_| false,
                |_| (1.0, [0.0; 3]),
            )
            .unwrap();
            let mut grid = OpenBoundaryGrid3::try_new(original, cx, reservoirs(axis, 1.0)).unwrap();
            grid.set_pressure_profile(
                Face3::ALL[2 * axis],
                &CancelGate::new_clock_free(),
                |point| {
                    Ok(
                        if point[(axis + 1) % 3] == 3 && point[(axis + 2) % 3] == 4 {
                            1.002
                        } else {
                            1.0
                        },
                    )
                },
            )
            .unwrap();
            let mut before = extensive(&grid, origin);
            // The localized perturbation cannot reach a transverse wall in
            // two native nearest-neighbor ticks. The untouched uniform wall
            // pressures cancel, so the open receipt owns the entire change.
            for _ in 0..2 {
                let receipt = grid.step_serial(cx, None, origin).unwrap();
                let after = extensive(&grid, origin);
                near(after.0 - before.0, receipt.mass_change, 3e-11);
                for a in 0..3 {
                    near(after.1[a] - before.1[a], receipt.impulse[a], 3e-11);
                    near(after.2[a] - before.2[a], receipt.angular_impulse[a], 3e-11);
                }
                before = after;
            }
            assert!(
                before.1[axis] > 1e-5,
                "pressure patch must drive inward momentum"
            );
            assert!(
                before.2.iter().map(|v| v * v).sum::<f64>().sqrt() > 1e-6,
                "off-axis pressure patch must impart a nonzero moment"
            );
            // Before a wall can respond, symmetric lateral spreading retains
            // the patch's center of pressure. This checks the physical moment
            // arm independently of the reservoir's own angular receipt.
            let mut center = [0.5; 3];
            center[(axis + 1) % 3] = 3.5;
            center[(axis + 2) % 3] = 4.5;
            let arm: [f64; 3] = core::array::from_fn(|a| center[a] - origin[a]);
            for a in 0..3 {
                let (b, c) = ((a + 1) % 3, (a + 2) % 3);
                near(
                    before.2[a],
                    arm[b] * before.1[c] - arm[c] * before.1[b],
                    3e-11,
                );
            }
        }
    });
}

#[test]
fn localized_pressure_transmits_linear_and_angular_momentum_in_native_collision_families() {
    pressure_patch::<D3q19, 19>(CollisionModel3::Bgk { tau: 0.9 });
    for model in [
        CollisionModel27::Bgk {
            tau: 0.9,
            equilibrium: Equilibrium27::Quadratic,
        },
        CollisionModel27::Trt {
            tau: 0.9,
            magic: 3.0 / 16.0,
            equilibrium: Equilibrium27::Factorized,
        },
        CollisionModel27::CentralMoment {
            rates: Rates27 { omega: [1.2; 10] },
        },
        CollisionModel27::Cumulant {
            rates: Rates27 { omega: [1.2; 10] },
        },
    ] {
        pressure_patch::<D3q27, 27>(config(model));
    }
}

fn fixture(cx: &Cx<'_>, solid: impl Fn([usize; 3]) -> bool) -> OpenBoundaryGrid3<D3q27, 27> {
    let original = Grid3::try_new_equilibrium(
        [6, 7, 5],
        config(CollisionModel27::Bgk {
            tau: 0.9,
            equilibrium: Equilibrium27::Factorized,
        }),
        [2e-5, -1e-5, 3e-5],
        [false; 3],
        budget(),
        solid,
        |_| (1.0, [0.002, -0.004, 0.006]),
    )
    .unwrap();
    OpenBoundaryGrid3::try_new(original, cx, reservoirs(2, 1.0)).unwrap()
}

fn field(point: [usize; 3]) -> f64 {
    1.01 + 0.0001 * point[0] as f64 + 0.0002 * point[1] as f64
}

fn install(grid: &mut OpenBoundaryGrid3<D3q27, 27>) {
    grid.set_pressure_profile(Face3::ZMin, &CancelGate::new_clock_free(), |point| {
        Ok(field(point))
    })
    .unwrap();
}

#[test]
fn pressure_profiles_keep_wall_owned_populations_and_pooled_exchange() {
    with_cx(|cx| {
        let origin = [0.6, 0.7, -0.4];
        let mut ordinary = fixture(cx, |_| false);
        let mut profiled = fixture(cx, |_| false);
        install(&mut profiled);
        ordinary.step_serial(cx, None, origin).unwrap();
        profiled.step_serial(cx, None, origin).unwrap();
        let dims = profiled.dimensions();
        let mut checked = 0;
        let mut changed = 0;
        for y in 0..dims[1] {
            for x in 0..dims[0] {
                if x != 0 && x + 1 != dims[0] && y != 0 && y + 1 != dims[1] {
                    continue;
                }
                let point = [x, y, 0];
                let cell = address(point, dims);
                for q in 0..27 {
                    let direction = <D3q27 as Lattice3<27>>::E[q];
                    let crosses_wall = (0..2).any(|a| {
                        let source = point[a] as isize - direction[a] as isize;
                        source < 0 || source >= dims[a] as isize
                    });
                    let open_owned = direction[2] > 0 && !crosses_wall;
                    if !open_owned {
                        assert_eq!(
                            profiled.populations()[cell][q].to_bits(),
                            ordinary.populations()[cell][q].to_bits()
                        );
                        checked += 1;
                    } else if profiled.populations()[cell][q].to_bits()
                        != ordinary.populations()[cell][q].to_bits()
                    {
                        changed += 1;
                    }
                }
            }
        }
        assert!(checked > 0 && changed > 0);
        near(
            profiled.macroscopic(address([2, 3, 0], dims)).unwrap().0,
            field([2, 3, 0]),
            2e-13,
        );
        assert!(
            (profiled.macroscopic(0).unwrap().0 - field([0, 0, 0])).abs() > 1e-5,
            "a wall seam must not silently become a complete pressure face"
        );

        let mut serial = fixture(cx, |_| false);
        install(&mut serial);
        let pointer = serial.pressure_profile(Face3::ZMin).unwrap().as_ptr();
        let receipts: Vec<_> = (0..4)
            .map(|_| serial.step_serial(cx, None, origin).unwrap())
            .collect();
        assert_eq!(
            serial.pressure_profile(Face3::ZMin).unwrap().as_ptr(),
            pointer
        );
        for workers in [1, 3] {
            let mut pooled = fixture(cx, |_| false);
            install(&mut pooled);
            let pool = TilePool::for_host(workers, SEED);
            let gate = CancelGate::new_clock_free();
            for expected in &receipts {
                assert_eq!(
                    &pooled.step_pooled(&pool, &gate, None, origin).unwrap(),
                    expected
                );
            }
            assert_eq!(bits(&pooled), bits(&serial));
        }
    });
}

#[test]
fn native_profile_updates_are_atomic_and_reuse_face_storage() {
    with_cx(|cx| {
        let gate = CancelGate::new_clock_free();
        let mut grid = fixture(cx, |_| false);
        let before = bits(&grid);
        let mut pointers = Vec::new();
        for offset in [0.0, 0.001, 0.002] {
            let mut sampled = Vec::new();
            let report = grid
                .set_pressure_profile(Face3::ZMin, &gate, |point| {
                    sampled.push(address(point, [6, 7, 5]));
                    Ok(field(point) + offset)
                })
                .unwrap();
            assert_eq!(report.sampled_cells, 42);
            assert!(sampled.windows(2).all(|pair| pair[0] < pair[1]));
            pointers.push(grid.pressure_profile(Face3::ZMin).unwrap().as_ptr());
        }
        assert_ne!(pointers[0], pointers[1]);
        assert_eq!(pointers[0], pointers[2]);
        let held = grid.pressure_profile(Face3::ZMin).unwrap().to_vec();
        assert_eq!(
            grid.set_pressure_profile(Face3::ZMin, &gate, |point| {
                Ok(if point == [5, 6, 0] { f64::NAN } else { 1.1 })
            }),
            Err(BoundaryProfileError3::InvalidDensity { cell: [5, 6, 0] })
        );
        let reason = BoundaryProfileError3::Source("pressure trace unavailable".into());
        assert_eq!(
            grid.set_pressure_profile(Face3::ZMin, &gate, |point| {
                if point == [5, 6, 0] {
                    Err(reason.clone())
                } else {
                    Ok(1.1)
                }
            }),
            Err(reason)
        );
        let cancelled = CancelGate::new_clock_free();
        assert_eq!(
            grid.set_pressure_profile(Face3::ZMin, &cancelled, |point| {
                if point == [5, 6, 0] {
                    cancelled.request();
                }
                Ok(1.1)
            }),
            Err(BoundaryProfileError3::Cancelled)
        );
        assert_eq!(grid.pressure_profile(Face3::ZMin).unwrap(), held);
        assert_eq!(bits(&grid), before);
        assert_eq!(grid.completed_steps(), 0);

        let old = grid.boundaries();
        let mut faces = Face3::ALL.map(|face| old.face(face));
        faces[4] = FaceBoundary3::Pressure { density: 1.004 };
        assert!(matches!(
            grid.rebind_open_boundaries(BoundarySpec3::new(faces), &gate),
            Err(GridError3::Input(_))
        ));
        assert_eq!(grid.boundaries(), old);
        grid.rebind_open_boundaries(old, &gate).unwrap();
        grid.set_pressure_profile(Face3::ZMin, &gate, |point| Ok(field(point)))
            .unwrap();
        assert_eq!(
            grid.pressure_profile(Face3::ZMin).unwrap().as_ptr(),
            pointers[1]
        );
        grid.clear_pressure_profile(Face3::ZMin);
        assert!(grid.pressure_profile(Face3::ZMin).is_none());
        grid.rebind_open_boundaries(BoundarySpec3::new(faces), &gate)
            .unwrap();
        let mut control = fixture(cx, |_| false);
        control
            .rebind_open_boundaries(BoundarySpec3::new(faces), &gate)
            .unwrap();
        assert_eq!(
            grid.step_serial(cx, None, [0.0; 3]).unwrap(),
            control.step_serial(cx, None, [0.0; 3]).unwrap()
        );
        assert_eq!(bits(&grid), bits(&control));

        assert_eq!(
            grid.set_pressure_profile(Face3::XMin, &gate, |_| panic!("wall sampled")),
            Err(BoundaryProfileError3::NotPressureFace)
        );
        let mut blocked = fixture(cx, |point| point == [2, 3, 0]);
        let report = blocked
            .set_pressure_profile(Face3::ZMin, &gate, |point| {
                assert_ne!(point, [2, 3, 0]);
                Ok(field(point))
            })
            .unwrap();
        assert_eq!(report.sampled_cells, 41);
        assert_eq!(blocked.pressure_profile(Face3::ZMin).unwrap().len(), 42);
        let mut empty = fixture(cx, |point| point[2] == 0);
        assert_eq!(
            empty.set_pressure_profile(Face3::ZMin, &gate, |_| panic!("solid sampled")),
            Err(BoundaryProfileError3::EmptyFace)
        );
        assert!(empty.pressure_profile(Face3::ZMin).is_none());
    });
}
