//! G0/G1 resolved-curl, physical-scaling and native-boundary regressions.
//! The Taylor--Green decay fixture is the exact 2-D solution extruded in z;
//! the genuinely 3-D fixture checks spatial curl, not an invented decay law.

use fs_exec::{CancelGate, Cx};
use fs_lbm::d3q19::{BoundaryGrid3, BoundarySpec3, CollisionModel3, FaceBoundary3, equilibrium3};
use fs_lbm::d3q27::{CollisionConfig27, CollisionModel27, D3q27, Equilibrium27, Rates27};
use fs_lbm::kinematics::{FlowKinematics3, KinematicsError, KinematicsWorkspace3};
use fs_lbm::lattice3::sources::{ForcedGrid3, ForcedLattice3};
use fs_lbm::lattice3::{D3q19, Grid3, GridBudget3, Lattice3};
use fs_lbm::nsf3::{
    NsfCharacteristicAxis, NsfCharacteristicBoundary, NsfGrid3, NsfSpec, WallThermal,
};
use fs_material::gas::GasSpec;

const TAU: f64 = core::f64::consts::TAU;

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} != {expected} (tol {tolerance})"
    );
}

fn workspace(dims: [usize; 3]) -> KinematicsWorkspace3 {
    KinematicsWorkspace3::new(dims, KinematicsWorkspace3::required_bytes(dims).unwrap()).unwrap()
}

fn moving_wall_curl<L: Lattice3<Q> + ForcedLattice3<Q>, const Q: usize>(collision: L::Collision) {
    let gate = CancelGate::new_clock_free();
    for normal in 0..3 {
        let tangent = (normal + 1) % 3;
        for height in [1, 7] {
            let mut dims = [3; 3];
            dims[normal] = height;
            let cells = dims.into_iter().product();
            for (linear, quadratic) in [(0.0, 0.0), (0.002, 0.0), (0.001, 0.0003)] {
                for spatial_force in [false, true] {
                    let force = |p: [usize; 3]| {
                        let mut f = [0.0; 3];
                        if spatial_force {
                            f[tangent] = 0.0004 * (p[normal] as f64 + 0.5);
                        }
                        f
                    };
                    let grid = Grid3::<L, Q>::try_new(
                        dims,
                        collision,
                        [0.0; 3],
                        core::array::from_fn(|axis| axis != normal),
                        GridBudget3 {
                            maximum_cells: cells,
                            maximum_mach: 0.2,
                        },
                        |_| false,
                        |p| {
                            let h = p[normal] as f64 + 0.5;
                            let rho = 1.4 + 0.002 * h;
                            let mut u = force(p).map(|f| -0.5 * f / rho);
                            u[tangent] += 0.01 + linear * h + quadratic * h * (height as f64 - h);
                            equilibrium::<L, Q>(rho, u)
                        },
                    )
                    .unwrap();
                    let probe =
                        with_cx(|cx| grid.wall_momentum_probe(cx, [0.0; 3], |_| true).unwrap());
                    let velocities: Vec<_> = probe
                        .wall_links()
                        .map(|link| {
                            let mut u = [0.0; 3];
                            u[tangent] = 0.01
                                + if link.exterior_faces[2 * normal + 1] {
                                    linear * height as f64
                                } else {
                                    0.0
                                };
                            u
                        })
                        .collect();
                    let mut w = workspace(dims);
                    let field = if spatial_force {
                        let forces = (0..cells)
                            .map(|i| {
                                force([i % dims[0], i / dims[0] % dims[1], i / dims[0] / dims[1]])
                            })
                            .collect();
                        let driven = with_cx(|cx| ForcedGrid3::try_new(grid, cx, forces).unwrap());
                        driven
                            .kinematics_with_moving_walls(&gate, &probe, &velocities, &mut w)
                            .unwrap()
                            .to_vec()
                    } else {
                        grid.kinematics_with_moving_walls(&gate, &probe, &velocities, &mut w)
                            .unwrap()
                            .to_vec()
                    };
                    for (i, sample) in field.into_iter().enumerate() {
                        let p = [i % dims[0], i / dims[0] % dims[1], i / dims[0] / dims[1]];
                        let mut gradient = [[0.0; 3]; 3];
                        gradient[tangent][normal] =
                            linear + quadratic * (height as f64 - 2.0 * (p[normal] as f64 + 0.5));
                        for (a, e) in sample
                            .unwrap()
                            .gradient()
                            .into_iter()
                            .flatten()
                            .zip(gradient.into_iter().flatten())
                        {
                            near(a, e, 3e-16);
                        }
                    }
                }
            }
        }
    }
}

#[test]
fn native_moving_wall_curl_preserves_comotion_couette_and_poiseuille_on_both_lattices() {
    moving_wall_curl::<D3q19, 19>(CollisionModel3::Bgk { tau: 0.8 });
    moving_wall_curl::<D3q27, 27>(CollisionConfig27 {
        model: CollisionModel27::Bgk {
            tau: 0.8,
            equilibrium: Equilibrium27::Factorized,
        },
        maximum_mach: 0.2,
    });
}

#[test]
fn spatial_moving_wall_curl_uses_native_axial_links_and_refuses_invalid_motion() {
    let dims = [12, 4, 3];
    let gate = CancelGate::new_clock_free();
    let k = TAU / dims[0] as f64;
    let mut grid = Grid3::<D3q19, 19>::try_new(
        dims,
        CollisionModel3::Bgk { tau: 0.8 },
        [0.0; 3],
        [true, false, true],
        GridBudget3 {
            maximum_cells: 144,
            maximum_mach: 0.2,
        },
        |_| false,
        |p| {
            equilibrium3(
                1.0,
                [
                    0.0,
                    0.0,
                    0.004 * (k * (p[0] as f64 + 0.5)).sin() * (p[1] as f64 + 0.5),
                ],
            )
        },
    )
    .unwrap();
    let probe = with_cx(|cx| {
        grid.wall_momentum_probe(cx, [0.0; 3], |link| link.exterior_faces[3])
            .unwrap()
    });
    // Use each physical link midpoint, including diagonal links. The native
    // axial link is at the correct x for the normal derivative; averaging all
    // wall links would incorrectly smooth this spatially varying boundary.
    let velocities: Vec<_> = probe
        .wall_links()
        .map(|link| [0.0, 0.0, 0.016 * (k * link.position[0]).sin()])
        .collect();
    let mut w = workspace(dims);
    let before = grid.populations().to_vec();
    let field = grid
        .kinematics_with_moving_walls(&gate, &probe, &velocities, &mut w)
        .unwrap();
    let pointer = field.as_ptr();
    for (i, sample) in field.iter().enumerate() {
        let x = i % dims[0];
        let y = i / dims[0] % dims[1];
        let omega = sample.unwrap().vorticity();
        near(omega[0], 0.004 * (k * (x as f64 + 0.5)).sin(), 1e-16);
        near(
            omega[1],
            -0.004 * k.sin() * (k * (x as f64 + 0.5)).cos() * (y as f64 + 0.5),
            1e-16,
        );
        near(omega[2], 0.0, 1e-16);
    }
    assert_eq!(before, grid.populations());
    assert_eq!(grid.completed_steps(), 0);
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(
        grid.kinematics_with_moving_walls(&cancelled, &probe, &velocities, &mut w)
            .is_err()
    );
    for invalid in [
        vec![[0.0; 3]; velocities.len() - 1],
        vec![[0.0, 0.01, 0.0]; velocities.len()],
        vec![[f64::NAN, 0.0, 0.0]; velocities.len()],
    ] {
        assert!(
            grid.kinematics_with_moving_walls(&gate, &probe, &invalid, &mut w)
                .is_err()
        );
    }
    assert!(
        grid.kinematics_with_moving_walls(&gate, &probe, &velocities, &mut workspace([6, 8, 3]))
            .is_err()
    );
    assert_eq!(
        pointer,
        grid.kinematics_with_moving_walls(&gate, &probe, &velocities, &mut w)
            .unwrap()
            .as_ptr()
    );
    assert_eq!(before, grid.populations());
    let zero = vec![[0.0; 3]; velocities.len()];
    let stationary = grid.kinematics(&gate, &mut w).unwrap().to_vec();
    assert_eq!(
        stationary,
        grid.kinematics_with_moving_walls(&gate, &probe, &zero, &mut w)
            .unwrap()
    );
    // The same admitted drive must be usable by the actual momentum owner.
    with_cx(|cx| {
        grid.step_serial_with_moving_wall_momentum(cx, &probe, &velocities, None)
            .unwrap()
    });
    assert_eq!(grid.completed_steps(), 1);
}

#[test]
fn finite_wall_curl_uses_accepted_endpoint_instead_of_last_streaming_midpoint() {
    use fs_lbm::d3q19::Face3;
    use fs_lbm::lattice3::finite_wall::{FinitePlanarWall3, FiniteWallBudget3, PlanarWallState3};
    let dims = [4, 5, 6];
    let gate = CancelGate::new_clock_free();
    let mut grid = Grid3::<D3q19, 19>::try_new(
        dims,
        CollisionModel3::Bgk { tau: 0.8 },
        [0.0; 3],
        [true, false, true],
        GridBudget3 {
            maximum_cells: 120,
            maximum_mach: 0.2,
        },
        |_| false,
        |_| equilibrium3(1.0, [0.0; 3]),
    )
    .unwrap();
    with_cx(|cx| {
        let mut wall = FinitePlanarWall3::new(
            &grid,
            cx,
            Face3::YMax,
            [2.0, 5.0, 3.0],
            PlanarWallState3 {
                mass: 4.0,
                axial_inertia: 8.0,
                velocity: [0.01, -0.006, 0.001],
            },
            |_| true,
        )
        .unwrap();
        let step = grid
            .step_serial_with_finite_wall(
                cx,
                &mut wall,
                [0.003, 0.001, -0.0004],
                FiniteWallBudget3 {
                    maximum_relative_impulse_defect: 1e-10,
                    maximum_relative_work_defect: 1e-10,
                },
                None,
            )
            .unwrap();
        assert_ne!(step.after.velocity, step.velocity_midpoint);
        let probe = grid
            .wall_momentum_probe(cx, [0.0; 3], |link| link.exterior_faces[3])
            .unwrap();
        let drives = |v: [f64; 3]| {
            probe
                .wall_links()
                .map(|link| {
                    let [x, _, z] = link.fluid_cell.map(|p| p as f64 + 0.5);
                    [v[1] + v[2] * (z - 3.0), 0.0, v[0] - v[2] * (x - 2.0)]
                })
                .collect::<Vec<_>>()
        };
        let mut w = workspace(dims);
        let expected = grid
            .kinematics_with_moving_walls(&gate, &probe, &drives(step.after.velocity), &mut w)
            .unwrap()
            .to_vec();
        let midpoint = grid
            .kinematics_with_moving_walls(&gate, &probe, &drives(step.velocity_midpoint), &mut w)
            .unwrap()
            .to_vec();
        assert_ne!(expected, midpoint);
        let before = grid.populations().to_vec();
        assert_eq!(
            grid.kinematics_with_finite_wall(&gate, &mut wall, &mut w)
                .unwrap(),
            expected
        );
        assert_eq!(wall.state(), step.after);
        assert_eq!(grid.populations(), before);
        assert_eq!(grid.completed_steps(), 1);
    });
}

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    let arenas = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    arenas.scope(|arena| {
        f(&Cx::new(
            &gate,
            arena,
            fs_exec::StreamKey {
                seed: 0xC071,
                kernel_id: 3,
                tile: 0,
                iteration: 0,
            },
            fs_exec::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        ))
    })
}

fn equilibrium<L: Lattice3<Q>, const Q: usize>(rho: f64, u: [f64; 3]) -> [f64; Q] {
    let u2: f64 = u.iter().map(|v| v * v).sum();
    core::array::from_fn(|q| {
        let eu: f64 = (0..3).map(|a| f64::from(L::E[q][a]) * u[a]).sum();
        L::W[q] * rho * (1.0 + eu / L::CS2 + 0.5 * eu * eu / L::CS2.powi(2) - 0.5 * u2 / L::CS2)
    })
}

#[test]
fn rotation_shear_strain_stretching_and_units_are_distinct() {
    let rotating =
        FlowKinematics3::from_gradient([[0.0, -0.3, 0.0], [0.3, 0.0, 0.0], [0.0; 3]]).unwrap();
    assert_eq!(rotating.vorticity(), [0.0, 0.0, 0.6]);
    near(rotating.enstrophy_density(), 0.18, 1e-15);
    near(rotating.q_criterion(), 0.09, 1e-15);
    assert_eq!(rotating.vortex_stretching(), [0.0; 3]);
    let strain =
        FlowKinematics3::from_gradient([[0.3, 0.0, 0.0], [0.0, -0.3, 0.0], [0.0; 3]]).unwrap();
    assert_eq!(strain.vorticity(), [0.0; 3]);
    near(strain.q_criterion(), -0.09, 1e-15);
    let shear = FlowKinematics3::from_gradient([[0.0, 0.4, 0.0], [0.0; 3], [0.0; 3]]).unwrap();
    assert_eq!(shear.vorticity(), [0.0, 0.0, -0.4]);
    assert_eq!(shear.q_criterion(), 0.0);
    let stretched =
        FlowKinematics3::from_gradient([[-0.1, -0.3, 0.0], [0.3, -0.1, 0.0], [0.0, 0.0, 0.2]])
            .unwrap();
    assert_eq!(stretched.divergence(), 0.0);
    near(stretched.vortex_stretching()[2], 0.12, 1e-15);
    let si = stretched.to_physical(0.25).unwrap();
    near(si.vorticity()[2], 2.4, 1e-15);
    near(
        si.enstrophy_density(),
        16.0 * stretched.enstrophy_density(),
        1e-14,
    );
    near(si.vortex_stretching()[2], 1.92, 1e-14);
    for dt in [0.0, -1.0, f64::INFINITY, f64::NAN] {
        assert!(stretched.to_physical(dt).is_err());
    }
    assert!(FlowKinematics3::from_gradient([[f64::MAX; 3]; 3]).is_err());
}

#[test]
fn moving_and_stationary_halfway_wall_curl_is_exact_in_every_orientation() {
    let gate = CancelGate::new_clock_free();
    for normal in 0..3 {
        let flow = (normal + 1) % 3;
        for quadratic in [false, true] {
            let mut faces = [FaceBoundary3::Periodic; 6];
            faces[2 * normal] = FaceBoundary3::Wall { velocity: [0.0; 3] };
            let mut wall = [0.0; 3];
            if !quadratic {
                wall[flow] = 0.024;
            }
            faces[2 * normal + 1] = FaceBoundary3::Wall { velocity: wall };
            let mut grid = BoundaryGrid3::new(8, 8, 8, 0.8, [0.0; 3], BoundarySpec3::new(faces));
            let force = |p: [usize; 3]| {
                let mut f = [0.0; 3];
                f[flow] = 0.0002 * (p[normal] as f64 + 0.5);
                f
            };
            grid.set_force_field(force);
            for z in 0..8 {
                for y in 0..8 {
                    for x in 0..8 {
                        let p = [x, y, z];
                        let h = p[normal] as f64 + 0.5;
                        let mut raw = force(p).map(|f| -0.5 * f);
                        raw[flow] += if quadratic {
                            0.001 * h * (8.0 - h)
                        } else {
                            0.003 * h
                        };
                        grid.set_populations(x, y, z, &equilibrium3(1.0, raw));
                    }
                }
            }
            let mut w = workspace([8; 3]);
            let field = grid.kinematics(&gate, &mut w).unwrap();
            for (cell, actual) in field.iter().enumerate() {
                let p = [cell % 8, cell / 8 % 8, cell / 64];
                let mut gradient = [[0.0; 3]; 3];
                gradient[flow][normal] = if quadratic {
                    0.001 * (8.0 - 2.0 * (p[normal] as f64 + 0.5))
                } else {
                    0.003
                };
                let expected = FlowKinematics3::from_gradient(gradient).unwrap();
                for (a, e) in actual
                    .unwrap()
                    .vorticity()
                    .into_iter()
                    .zip(expected.vorticity())
                {
                    near(a, e, 2e-15);
                }
                near(actual.unwrap().divergence(), 0.0, 2e-15);
            }
        }
    }
}

fn taylor_green_curl_error<L: Lattice3<Q>, const Q: usize>(
    n: usize,
    collision: L::Collision,
) -> f64 {
    let amplitude = 0.01;
    let k = TAU / n as f64;
    let grid = Grid3::<L, Q>::try_new(
        [n; 3],
        collision,
        [0.0; 3],
        [true; 3],
        GridBudget3 {
            maximum_cells: n.pow(3),
            maximum_mach: 0.2,
        },
        |_| false,
        |p| {
            let [x, y, z] = p.map(|q| k * (q as f64 + 0.5));
            equilibrium::<L, Q>(
                1.0,
                [
                    amplitude * x.sin() * y.cos() * z.cos(),
                    -amplitude * x.cos() * y.sin() * z.cos(),
                    0.0,
                ],
            )
        },
    )
    .unwrap();
    let mut w = workspace([n; 3]);
    let gate = CancelGate::new_clock_free();
    let field = grid.kinematics(&gate, &mut w).unwrap();
    let mut error2 = 0.0;
    let mut exact2 = 0.0;
    for (i, sample) in field.iter().enumerate() {
        let [x, y, z] = [i % n, i / n % n, i / n / n].map(|q| k * (q as f64 + 0.5));
        let exact = [
            -amplitude * k * x.cos() * y.sin() * z.sin(),
            -amplitude * k * x.sin() * y.cos() * z.sin(),
            2.0 * amplitude * k * x.sin() * y.sin() * z.cos(),
        ];
        for (a, e) in sample.unwrap().vorticity().into_iter().zip(exact) {
            error2 += (a - e).powi(2);
            exact2 += e * e;
        }
        near(sample.unwrap().divergence(), 0.0, 1e-15);
    }
    (error2 / exact2).sqrt()
}

#[test]
fn three_dimensional_taylor_green_curl_converges_at_second_order_on_both_lattices() {
    let d27 = CollisionConfig27 {
        model: CollisionModel27::Cumulant {
            rates: Rates27 { omega: [1.25; 10] },
        },
        maximum_mach: 0.2,
    };
    for errors in [
        [8, 16, 32]
            .map(|n| taylor_green_curl_error::<D3q19, 19>(n, CollisionModel3::Bgk { tau: 0.8 })),
        [8, 16, 32].map(|n| taylor_green_curl_error::<D3q27, 27>(n, d27)),
    ] {
        for pair in errors.windows(2) {
            assert!(
                pair[0] / pair[1] > 3.8,
                "second-order curl ladder: {errors:?}"
            );
        }
    }
}

#[test]
fn native_step_diffuses_vorticity_with_analytic_taylor_green_enstrophy_decay() {
    let n = 24usize;
    let k = TAU / n as f64;
    let amplitude = 0.003;
    let nu = 0.1;
    let dims = [n, n, 4];
    let mut grid = Grid3::<D3q27, 27>::try_new(
        dims,
        CollisionConfig27 {
            model: CollisionModel27::Bgk {
                tau: 0.5 + 3.0 * nu,
                equilibrium: Equilibrium27::Factorized,
            },
            maximum_mach: 0.2,
        },
        [0.0; 3],
        [true; 3],
        GridBudget3 {
            maximum_cells: n * n * 4,
            maximum_mach: 0.2,
        },
        |_| false,
        |p| {
            let x = k * (p[0] as f64 + 0.5);
            let y = k * (p[1] as f64 + 0.5);
            let rho = 1.0 + 0.75 * amplitude * amplitude * ((2.0 * x).cos() + (2.0 * y).cos());
            let mut populations = equilibrium::<D3q27, 27>(
                rho,
                [
                    amplitude * x.sin() * y.cos(),
                    -amplitude * x.cos() * y.sin(),
                    0.0,
                ],
            );
            // Initialize the analytic viscous stress as well as velocity.
            // For this vortex S_xy=0 and S_xx=-S_yy=A k cos(x) cos(y).
            // The first Chapman--Enskog correction is
            // f_neq=-rho*tau*w_q*(c_i c_j-cs² delta_ij)*du_i/dx_j/cs².
            // An equilibrium-only start would add a kinetic startup transient
            // that is absent from the Navier--Stokes decay oracle.
            let strain = amplitude * k * x.cos() * y.cos();
            for (q, f) in populations.iter_mut().enumerate() {
                let e = D3q27::E[q].map(f64::from);
                *f -= rho * (0.5 + 3.0 * nu) * D3q27::W[q] * strain * (e[0] * e[0] - e[1] * e[1])
                    / D3q27::CS2;
            }
            populations
        },
    )
    .unwrap();
    let gate = CancelGate::new_clock_free();
    let mut w = workspace(dims);
    let initial: f64 = grid
        .kinematics(&gate, &mut w)
        .unwrap()
        .iter()
        .map(|v| v.unwrap().enstrophy_density())
        .sum();
    let steps = 48;
    with_cx(|cx| {
        for _ in 0..steps {
            grid.step_serial(cx).unwrap();
        }
    });
    let final_enstrophy: f64 = grid
        .kinematics(&gate, &mut w)
        .unwrap()
        .iter()
        .map(|v| v.unwrap().enstrophy_density())
        .sum();
    let expected_ratio = (-4.0 * nu * k * k * steps as f64).exp();
    near(
        final_enstrophy / initial,
        expected_ratio,
        0.01 * expected_ratio,
    );
    assert!(final_enstrophy < initial);
}

#[test]
fn forced_grid_curl_uses_the_accepted_spatial_half_force_and_preserves_solids() {
    let dims = [8; 3];
    let grid = Grid3::<D3q19, 19>::try_new(
        dims,
        CollisionModel3::Bgk { tau: 0.8 },
        [0.0; 3],
        [true; 3],
        GridBudget3 {
            maximum_cells: 512,
            maximum_mach: 0.2,
        },
        |p| p == [0; 3],
        |_| equilibrium3(1.0, [0.0; 3]),
    )
    .unwrap();
    let forces: Vec<_> = (0..512)
        .map(|i| {
            if i == 0 {
                [0.0; 3]
            } else {
                [0.0, 0.004 * (TAU * (i % 8) as f64 / 8.0).sin(), 0.0]
            }
        })
        .collect();
    let grid = with_cx(|cx| ForcedGrid3::try_new(grid, cx, forces).unwrap());
    let mut w = workspace(dims);
    let gate = CancelGate::new_clock_free();
    let field = grid.kinematics(&gate, &mut w).unwrap();
    assert!(field[0].is_none());
    for x in 0..8 {
        let curl = field[(3 * 8 + 3) * 8 + x].unwrap().vorticity();
        near(
            curl[2],
            0.002 * (TAU / 8.0).sin() * (TAU * x as f64 / 8.0).cos(),
            1e-16,
        );
    }
    let pointer = field.as_ptr();
    let before = grid.populations().to_vec();
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert_eq!(
        grid.kinematics(&cancelled, &mut w),
        Err(KinematicsError::Cancelled)
    );
    assert_eq!(pointer, grid.kinematics(&gate, &mut w).unwrap().as_ptr());
    assert_eq!(before, grid.populations());
    assert!(grid.kinematics(&gate, &mut workspace([4, 8, 16])).is_err());
    let bytes = KinematicsWorkspace3::required_bytes(dims).unwrap();
    assert!(matches!(
        KinematicsWorkspace3::new(dims, bytes - 1),
        Err(KinematicsError::Budget { .. })
    ));
    assert!(KinematicsWorkspace3::required_bytes([usize::MAX, 2, 2]).is_err());
}

fn nsf_spec(dims: [usize; 3], periodic: [bool; 3]) -> NsfSpec {
    NsfSpec {
        dims,
        periodic,
        gas: GasSpec::dry_air_ussa1976(),
        temperature_k: 293.15,
        pressure_pa: 101325.0,
        dx_m: 5e-5,
        temperature_window_k: [273.15, 333.15],
        maximum_mach: 0.2,
        maximum_temperature_ratio: 0.2,
        bulk_to_shear_ratio: 0.73,
        bulk_ratio_source: "Cramer 2012, Phys. Fluids 24, 066102, N2 at 293 K",
        maximum_cells: dims.into_iter().product(),
        wall: WallThermal::Adiabatic,
        workers: 1,
        transport_points: 9,
    }
}

#[test]
fn nsf_moving_wall_vorticity_and_open_throughflow_keep_physical_boundary_locations() {
    let dims = [4, 8, 4];
    let mut grid = NsfGrid3::<D3q27>::new(
        nsf_spec(dims, [true, false, true]),
        |_| false,
        |p| (1.0, [0.003 * (p[1] as f64 + 0.5), 0.0, 0.0], 293.15),
    )
    .unwrap();
    grid.set_edge_velocity(1, true, [0.024, 0.0, 0.0]).unwrap();
    let gate = CancelGate::new_clock_free();
    let mut w = workspace(dims);
    let state = grid.state_hash();
    for sample in grid.kinematics(&gate, &mut w).unwrap() {
        near(sample.unwrap().vorticity()[2], -0.003, 1e-15);
        near(
            sample
                .unwrap()
                .to_physical(grid.law().dt_s)
                .unwrap()
                .vorticity()[2],
            -0.003 / grid.law().dt_s,
            1e-8,
        );
    }
    assert_eq!(state, grid.state_hash());

    let velocity = [0.003, -0.002, 0.004];
    let through = NsfGrid3::<D3q27>::new(
        nsf_spec([6; 3], [false, true, true]),
        |_| false,
        |_| (1.0, velocity, 293.15),
    )
    .unwrap();
    let reference = velocity.map(|u| u * through.law().dx_m / through.law().dt_s);
    let through = through
        .with_characteristic_axis(NsfCharacteristicAxis {
            axis: 0,
            faces: [NsfCharacteristicBoundary {
                reference_pressure_pa: 101325.0,
                reference_temperature_k: 293.15,
                reference_velocity_m_s: reference,
                incoming_characteristic_m_s: 0.0,
            }; 2],
            maximum_relative_perturbation: 0.05,
        })
        .unwrap();
    for sample in through.kinematics(&gate, &mut workspace([6; 3])).unwrap() {
        near(sample.unwrap().enstrophy_density(), 0.0, 1e-25);
    }
}
