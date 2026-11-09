use super::*;

fn channel_faces(open: usize, normal: usize, inflow: f64, walls: [f64; 2]) -> BoundarySpec3 {
    let mut velocity = [0.0; 3];
    velocity[open] = inflow;
    let mut faces = [FaceBoundary3::Periodic; 6];
    faces[2 * open] = FaceBoundary3::Velocity { velocity };
    faces[2 * open + 1] = FaceBoundary3::Pressure { density: 1.1 };
    for side in 0..2 {
        velocity[open] = walls[side];
        faces[2 * normal + side] = FaceBoundary3::Wall { velocity };
    }
    BoundarySpec3::new(faces)
}

fn channel<L: EquilibriumLattice3<Q>, const Q: usize>(
    cx: &Cx<'_>,
    collision: L::Collision,
    open: usize,
    walls: [f64; 2],
) -> OpenBoundaryGrid3<L, Q> {
    let normal = (open + 1) % 3;
    let mut periodic = [true; 3];
    periodic[open] = false;
    periodic[normal] = false;
    let mut velocity = [0.0; 3];
    velocity[open] = 0.025;
    let grid = Grid3::try_new_equilibrium(
        [6, 5, 4],
        collision,
        [0.0; 3],
        periodic,
        budget(),
        |_| false,
        |_| (1.1, velocity),
    )
    .unwrap();
    OpenBoundaryGrid3::try_new(grid, cx, channel_faces(open, normal, 0.025, walls)).unwrap()
}

fn comotion<L: EquilibriumLattice3<Q>, const Q: usize>(collision: L::Collision) {
    with_cx(|cx| {
        for axis in 0..3 {
            let mut flow = channel::<L, Q>(cx, collision, axis, [0.025; 2]);
            assert!(flow.last_wall_exchange().is_none());
            for _ in 0..6 {
                let open = flow.step_serial(cx, None, [0.5, -0.25, 1.0]).unwrap();
                let wall = flow.last_wall_exchange().unwrap();
                assert!(open.open_links > 0 && wall.measured_links > 0);
                near(wall.fluid_mass_change, 0.0, 2e-14);
                near(wall.wall_work, 0.0, 2e-14);
                for cell in 0..flow.populations().len() {
                    let (rho, velocity) = flow.macroscopic(cell).unwrap();
                    near(rho, 1.1, 2e-13);
                    for a in 0..3 {
                        near(velocity[a], if a == axis { 0.025 } else { 0.0 }, 2e-13);
                    }
                }
            }
        }
    });
}

#[test]
fn uniform_through_flow_comoves_with_walls_in_both_native_stencils() {
    comotion::<D3q19, 19>(CollisionModel3::Bgk { tau: 0.9 });
    comotion::<D3q27, 27>(bgk());
}

fn impulse_balance<L: EquilibriumLattice3<Q>, const Q: usize>(collision: L::Collision) {
    with_cx(|cx| {
        let mut flow = channel::<L, Q>(cx, collision, 0, [0.0; 2]);
        let equilibrium = flow.populations()[0];
        let (mass_before, momentum_before) = totals(&flow);
        let targets = channel_faces(0, 1, 0.025, [0.0, 0.05]);
        let open = flow
            .step_serial_with_open_boundaries(cx, targets, None, [1.25, -0.5, 2.0])
            .unwrap();
        let wall = flow.last_wall_exchange().unwrap();
        let (mass_after, momentum_after) = totals(&flow);
        assert_eq!(flow.boundaries(), targets);
        assert!(wall.wall_work < 0.0);
        near(wall.fluid_mass_change, 0.0, 2e-14);
        near(mass_after - mass_before, open.mass_change, 4e-12);
        for a in 0..3 {
            near(
                momentum_after[a] - momentum_before[a],
                open.impulse[a] + wall.fluid_population_impulse[a],
                4e-12,
            );
            near(
                wall.wall_impulse[a] + wall.fluid_population_impulse[a],
                wall.wall_velocity_mass_impulse[a],
                4e-14,
            );
            near(
                wall.wall_angular_impulse[a] + wall.fluid_population_angular_impulse[a],
                wall.wall_velocity_mass_angular_impulse[a],
                4e-13,
            );
        }
        // Each seam's complete reflected set remains wall owned. With this
        // uniform pre-collision state, every wall row has its analytic value.
        for x in [0, 5] {
            for (y, speed) in [(0, 0.0), (4, 0.05)] {
                let cell = address([x, y, 1], [6, 5, 4]);
                for q in 1..Q {
                    if (y == 0 && L::E[q][1] > 0) || (y == 4 && L::E[q][1] < 0) {
                        let expected = equilibrium[L::OPP[q]]
                            + 2.0 * L::W[q] * 1.1 * f64::from(L::E[q][0]) * speed / L::CS2;
                        near(flow.populations()[cell][q], expected, 3e-16);
                    }
                }
            }
        }
    });
}

#[test]
fn changing_wall_motion_closes_separate_reservoir_and_relative_wall_accounts() {
    impulse_balance::<D3q19, 19>(CollisionModel3::Bgk { tau: 0.9 });
    impulse_balance::<D3q27, 27>(bgk());
}

fn couette<L: EquilibriumLattice3<Q>, const Q: usize>(collision: L::Collision) {
    with_cx(|cx| {
        let make = || {
            Grid3::<L, Q>::try_new_equilibrium(
                [8, 6, 3],
                collision,
                [0.0; 3],
                [true, false, true],
                budget(),
                |_| false,
                |_| (1.0, [0.0; 3]),
            )
            .unwrap()
        };
        let mut faces = [FaceBoundary3::Periodic; 6];
        faces[2] = FaceBoundary3::stationary_wall();
        faces[3] = FaceBoundary3::Wall {
            velocity: [0.04, 0.0, 0.0],
        };
        let mut open = OpenBoundaryGrid3::try_new(make(), cx, BoundarySpec3::new(faces)).unwrap();
        let mut closed = make();
        let probe = closed.wall_momentum_probe(cx, [0.0; 3], |_| true).unwrap();
        let velocities: Vec<_> = probe
            .wall_links()
            .map(|link| {
                if link.exterior_faces[3] {
                    [0.04, 0.0, 0.0]
                } else {
                    [0.0; 3]
                }
            })
            .collect();
        for _ in 0..400 {
            let reservoir = open.step_serial(cx, None, [0.0; 3]).unwrap();
            let wall = closed
                .step_serial_with_moving_wall_momentum(cx, &probe, &velocities, None)
                .unwrap();
            assert_eq!(reservoir.open_links, 0);
            assert_eq!(open.last_wall_exchange(), Some(wall));
            assert_eq!(open.populations(), closed.populations());
        }
        for y in 0..6 {
            let (_, velocity) = open.macroscopic(address([3, y, 1], [8, 6, 3])).unwrap();
            near(velocity[0], 0.04 * (y as f64 + 0.5) / 6.0, 2e-8);
            near(velocity[1], 0.0, 2e-13);
            near(velocity[2], 0.0, 2e-13);
        }
    });
}

#[test]
fn moving_open_owner_matches_native_closed_kernel_and_analytic_couette() {
    couette::<D3q19, 19>(CollisionModel3::Bgk { tau: 0.9 });
    couette::<D3q27, 27>(bgk());
}

#[test]
fn wall_rebinding_is_atomic_and_rejects_normal_motion_and_incompatible_junctions() {
    with_cx(|cx| {
        let mut flow = channel::<D3q27, 27>(cx, bgk(), 0, [0.025; 2]);
        flow.step_serial(cx, None, [0.0; 3]).unwrap();
        let before = bits(&flow);
        let accepted = flow.boundaries();
        let receipt = flow.last_wall_exchange();
        let mut invalid: [FaceBoundary3; 6] =
            core::array::from_fn(|a| accepted.face(Face3::ALL[a]));
        invalid[3] = FaceBoundary3::Wall {
            velocity: [0.025, 1e-20, 0.0],
        };
        assert!(
            flow.rebind_open_boundaries(BoundarySpec3::new(invalid), &CancelGate::new_clock_free())
                .is_err()
        );
        let cancelled = CancelGate::new_clock_free();
        cancelled.request();
        assert!(matches!(
            flow.rebind_open_boundaries(channel_faces(0, 1, 0.025, [-0.02, 0.04]), &cancelled),
            Err(GridError3::Cancelled)
        ));
        assert_eq!(bits(&flow), before);
        assert_eq!(flow.boundaries(), accepted);
        assert_eq!(flow.completed_steps(), 1);
        assert_eq!(flow.last_wall_exchange(), receipt);

        let grid = Grid3::<D3q27, 27>::try_new_equilibrium(
            [4, 4, 4],
            bgk(),
            [0.0; 3],
            [true, false, false],
            budget(),
            |_| false,
            |_| (1.0, [0.0; 3]),
        )
        .unwrap();
        let mut corner = [FaceBoundary3::stationary_wall(); 6];
        corner[0] = FaceBoundary3::Periodic;
        corner[1] = FaceBoundary3::Periodic;
        corner[3] = FaceBoundary3::Wall {
            velocity: [0.02, 0.0, 0.0],
        };
        assert!(OpenBoundaryGrid3::try_new(grid, cx, BoundarySpec3::new(corner)).is_err());
    });
}

#[test]
fn spatial_force_owner_publishes_wall_reaction_with_the_accepted_force_epoch() {
    use fs_lbm::lattice3::open_boundaries::forcing::ForcedOpenBoundaryGrid3;

    with_cx(|cx| {
        let flow = channel::<D3q27, 27>(cx, bgk(), 0, [0.0, 0.05]);
        let cells = flow.populations().len();
        let mut forced =
            ForcedOpenBoundaryGrid3::try_new_preserving_velocity(flow, cx, vec![[0.0; 3]; cells])
                .unwrap();
        let momentum = |flow: &ForcedOpenBoundaryGrid3<D3q27, 27>| {
            (0..cells)
                .map(|cell| {
                    let (rho, velocity) = flow.macroscopic(cell).unwrap();
                    rho * velocity[0]
                })
                .sum::<f64>()
        };
        let before = momentum(&forced);
        let force = [0.00001, 0.0, 0.0];
        let open = forced
            .advance_serial(cx, None, vec![force; cells], [0.0; 3])
            .unwrap();
        let wall = forced.last_wall_exchange().unwrap();
        near(
            momentum(&forced) - before,
            0.5 * cells as f64 * force[0] + open.impulse[0] + wall.fluid_population_impulse[0],
            4e-12,
        );
        let populations = forced.populations().to_vec();
        let mut invalid = vec![force; cells];
        invalid[cells - 1][0] = f64::NAN;
        assert!(forced.advance_serial(cx, None, invalid, [0.0; 3]).is_err());
        assert_eq!(forced.populations(), populations);
        assert_eq!(forced.last_wall_exchange(), Some(wall));
        assert_eq!(forced.completed_steps(), 1);
    });
}

#[test]
fn cancelled_moving_wall_stream_preserves_populations_targets_and_last_reaction() {
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
        let mut flow = channel::<D3q27, 27>(cx, bgk(), 0, [0.025; 2]);
        flow.step_serial(cx, None, [0.0; 3]).unwrap();
        let before = bits(&flow);
        let targets = flow.boundaries();
        let receipt = flow.last_wall_exchange();
        let runner = CancelAfterStream {
            pool: TilePool::for_host(2, SEED),
            calls: AtomicUsize::new(0),
        };
        assert!(matches!(
            flow.step_pooled_with_open_boundaries(
                &runner,
                &CancelGate::new_clock_free(),
                channel_faces(0, 1, 0.025, [0.0, 0.05]),
                None,
                [0.0; 3]
            ),
            Err(GridError3::Cancelled)
        ));
        assert_eq!(bits(&flow), before);
        assert_eq!(flow.boundaries(), targets);
        assert_eq!(flow.last_wall_exchange(), receipt);
        assert_eq!(flow.completed_steps(), 1);
    });
}
