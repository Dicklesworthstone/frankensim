//! G0/G1/G3/G4/G5: evolving law epochs reach native open flow and work.

use super::*;
use fs_lbm::lattice3::open_boundaries::forcing::resistance::ResistanceOpenEndpointExchange3;

fn work_balance(r: ResistanceOpenEndpointExchange3) {
    let r = r.resistance;
    for w in [r.external_work, r.fluid_drag_work, r.skeleton_work] {
        near(
            w.impulse_work - w.endpoint_work,
            w.quadrature_remainder,
            2e-13,
        );
    }
    assert!(r.endpoint_relative_dissipation >= 0.0);
    near(
        r.fluid_drag_work.endpoint_work + r.skeleton_work.endpoint_work,
        -r.endpoint_relative_dissipation,
        2e-13,
    );
    near(
        r.fluid_drag_work.impulse_work + r.skeleton_work.impulse_work,
        -r.endpoint_relative_dissipation + r.relative_quadrature_remainder,
        2e-13,
    );
    assert!(
        r.absolute_relative_quadrature_remainder + 1e-15 >= r.relative_quadrature_remainder.abs()
    );
}

fn manufactured<L: EquilibriumLattice3<Q>, const Q: usize>(collision: L::Collision) {
    with_cx(|cx| {
        for boost in [[0.0; 3], [0.017, -0.009, 0.005]] {
            let n = 45;
            let rho = 1.2;
            let u0 = core::array::from_fn(|a| [0.022, -0.013, 0.018][a] + boost[a]);
            let u1 = core::array::from_fn(|a| [0.018, -0.008, 0.021][a] + boost[a]);
            let mut before = porous_law();
            before.linear.body_force = [0.0003, -0.0002, 0.0001];
            let mut after = DarcyForchheimer3 {
                linear: LinearResistance3 {
                    principal: [0.025, 0.06, 0.04],
                    skeleton_velocity: [-0.003, 0.009, 0.002],
                    ..LinearResistance3::default()
                },
                quadratic: 0.7,
            };
            for a in 0..3 {
                before.linear.skeleton_velocity[a] += boost[a];
                after.linear.skeleton_velocity[a] += boost[a];
            }
            let d0 = drag(before, u0);
            let d1 = drag(after, u1);
            // Manufacture the independent NEXT load from the physical momentum
            // equation; this does not call the native constitutive root solve.
            after.linear.body_force = core::array::from_fn(|a| {
                2.0 * rho * (u1[a] - u0[a]) - before.linear.body_force[a] - d0[a] - d1[a]
            });
            let grid = Grid3::<L, Q>::try_new_equilibrium(
                [5, 3, 3],
                collision,
                [0.0; 3],
                [false, true, true],
                budget(),
                |_| false,
                |_| (rho, u0),
            )
            .unwrap();
            let flow = OpenBoundaryGrid3::try_new(grid, cx, pressure_faces(rho)).unwrap();
            let mut porous = ResistanceOpenGrid3::try_new_preserving_velocity(
                flow,
                cx,
                vec![before; n],
                [0.5, -0.25, 1.0],
            )
            .unwrap();
            let report = porous
                .advance_serial(cx, None, vec![after; n], 1.0)
                .unwrap();
            let r = report.resistance;
            for cell in 0..n {
                let (density, u) = porous.flow().macroscopic(cell).unwrap();
                near(density, rho, 3e-13);
                for a in 0..3 {
                    near(u[a], u1[a], 3e-13);
                    near(
                        porous.flow().forces()[cell][a],
                        after.linear.body_force[a] + d1[a],
                        3e-13,
                    );
                }
            }
            let mean_u = core::array::from_fn(|a| 0.5 * (u0[a] + u1[a]));
            let reaction = core::array::from_fn(|a| -(n as f64) * 0.5 * (d0[a] + d1[a]));
            for a in 0..3 {
                near(
                    r.external_impulse[a],
                    n as f64 * 0.5 * (before.linear.body_force[a] + after.linear.body_force[a]),
                    3e-13,
                );
                near(r.skeleton_impulse[a], reaction[a], 3e-13);
            }
            near(
                r.fluid_drag_work.impulse_work,
                -dot(reaction, mean_u),
                3e-13,
            );
            let slip0 = core::array::from_fn(|a| u0[a] - before.linear.skeleton_velocity[a]);
            let slip1 = core::array::from_fn(|a| u1[a] - after.linear.skeleton_velocity[a]);
            near(
                r.endpoint_relative_dissipation,
                -0.5 * n as f64 * (dot(d0, slip0) + dot(d1, slip1)),
                3e-13,
            );
            near(
                0.5 * n as f64 * rho * (dot(u1, u1) - dot(u0, u0)),
                r.external_work.impulse_work + r.fluid_drag_work.impulse_work,
                3e-13,
            );
            work_balance(report);
            assert!(report.open.open_links > 0);
            assert_eq!(porous.laws(), vec![after; n]);
            assert_eq!(porous.last_endpoint_exchange(), Some(report));
            assert_eq!(porous.last_exchange(), None);
            porous.step_serial(cx, None).unwrap();
            assert_eq!(porous.last_endpoint_exchange(), None);
            assert!(porous.last_exchange().is_some());
            assert_eq!(porous.laws(), vec![after; n]);
        }
    });
}

#[test]
fn changing_skeleton_tensor_and_vector_speed_drag_have_the_authored_endpoint() {
    manufactured::<D3q19, 19>(CollisionModel3::Bgk { tau: 0.9 });
    manufactured::<D3q27, 27>(bgk());
}

fn channel<L: EquilibriumLattice3<Q>, const Q: usize>(
    cx: &Cx<'_>,
    collision: L::Collision,
) -> ResistanceOpenGrid3<L, Q> {
    let mut faces = [FaceBoundary3::Periodic; 6];
    faces[0] = FaceBoundary3::Velocity {
        velocity: [0.015, 0.0, 0.004],
    };
    faces[1] = FaceBoundary3::Pressure { density: 1.1 };
    faces[2] = FaceBoundary3::Wall {
        velocity: [0.0, 0.0, 0.006],
    };
    faces[3] = faces[2];
    let grid = Grid3::<L, Q>::try_new_equilibrium(
        [7, 5, 4],
        collision,
        [0.0; 3],
        [false, false, true],
        budget(),
        |p| p == [3, 2, 1],
        |p| {
            (
                1.1 + 0.0001 * p[0] as f64,
                [0.012, 0.001 * (p[1] as f64 - 2.0), 0.003],
            )
        },
    )
    .unwrap();
    let flow = OpenBoundaryGrid3::try_new(grid, cx, BoundarySpec3::new(faces)).unwrap();
    let laws = (0..140)
        .map(|cell| {
            let (x, y) = (cell % 7, cell / 7 % 5);
            if (x == 0 || x == 6) && (y == 0 || y == 4) {
                DarcyForchheimer3::default()
            } else {
                porous_law()
            }
        })
        .collect();
    ResistanceOpenGrid3::try_new_preserving_velocity(flow, cx, laws, [1.25, -0.75, 0.5]).unwrap()
}

fn next_laws(laws: &[DarcyForchheimer3], tick: usize) -> Vec<DarcyForchheimer3> {
    laws.iter()
        .enumerate()
        .map(|(cell, &law)| {
            let mut next = law;
            // Keep zero-drag seams exactly free while evolving every physical load.
            next.linear.principal = law
                .linear
                .principal
                .map(|d| d * (0.85 + 0.02 * (cell % 7) as f64));
            next.quadratic *= 0.9;
            next.linear.body_force = [0.00002 * tick as f64, -0.00001, 0.00003];
            next.linear.skeleton_velocity = [0.004 * tick as f64, -0.002, 0.003];
            next
        })
        .collect()
}

fn channel_balance<L: EquilibriumLattice3<Q>, const Q: usize>(collision: L::Collision) {
    with_cx(|cx| {
        let mut serial = channel::<L, Q>(cx, collision);
        let mut parallel = channel::<L, Q>(cx, collision);
        let pool = TilePool::for_host(3, SEED);
        for tick in 1..=3 {
            let laws = next_laws(serial.laws(), tick);
            let before = inventory(&serial, Some([3, 2, 1]));
            let report = serial.advance_serial(cx, None, laws.clone(), 1.0).unwrap();
            let pooled = parallel
                .advance_pooled(&pool, &CancelGate::new_clock_free(), None, laws, 1.0)
                .unwrap();
            let after = inventory(&serial, Some([3, 2, 1]));
            assert_eq!(snapshot(&serial), snapshot(&parallel));
            assert_eq!(serial.laws(), parallel.laws());
            assert_eq!(report, pooled);
            let r = report.resistance;
            near(
                after.0 - before.0,
                report.open.mass_change + r.walls.fluid_mass_change,
                2e-11,
            );
            for a in 0..3 {
                near(
                    after.1[a] - before.1[a],
                    r.external_impulse[a] - r.skeleton_impulse[a]
                        + report.open.impulse[a]
                        + r.walls.fluid_population_impulse[a],
                    2e-11,
                );
                near(
                    after.2[a] - before.2[a],
                    r.external_angular_impulse[a] - r.skeleton_angular_impulse[a]
                        + report.open.angular_impulse[a]
                        + r.walls.fluid_population_angular_impulse[a]
                        + r.periodic.angular_momentum_change[a],
                    5e-11,
                );
            }
            for cell in 0..140 {
                if cell == address([3, 2, 1], [7, 5, 4]) {
                    continue;
                }
                let (_, u) = serial.flow().macroscopic(cell).unwrap();
                let law = serial.laws()[cell];
                let d = drag(law, u);
                for a in 0..3 {
                    near(
                        serial.flow().forces()[cell][a],
                        law.linear.body_force[a] + d[a],
                        3e-13,
                    );
                }
            }
            work_balance(report);
            assert!(
                r.walls.measured_links > serial.flow().last_wall_exchange().unwrap().measured_links
            );
            assert!(r.periodic.face_crossings > 0 && report.open.open_links > 0);
        }
    });
}

#[test]
fn evolving_laws_close_reservoir_wall_voxel_and_periodic_balance_on_both_stencils() {
    channel_balance::<D3q19, 19>(CollisionModel3::Bgk { tau: 0.9 });
    channel_balance::<D3q27, 27>(bgk());
}

fn parity<L: EquilibriumLattice3<Q>, const Q: usize>(collision: L::Collision) {
    with_cx(|cx| {
        let mut held = channel::<L, Q>(cx, collision);
        let mut endpoint = channel::<L, Q>(cx, collision);
        for _ in 0..3 {
            let laws = held.laws().to_vec();
            let a = held.step_serial(cx, None).unwrap();
            let b = endpoint.advance_serial(cx, None, laws, 1.0).unwrap();
            assert_eq!(held.flow().populations(), endpoint.flow().populations());
            assert_eq!(held.flow().forces(), endpoint.flow().forces());
            assert_eq!(a.open, b.open);
            assert_eq!(a.walls, b.resistance.walls);
            assert_eq!(a.periodic, b.resistance.periodic);
            assert_eq!(a.external_impulse, b.resistance.external_impulse);
            assert_eq!(a.skeleton_impulse, b.resistance.skeleton_impulse);
            assert_eq!(
                a.external_work.to_bits(),
                b.resistance.external_work.impulse_work.to_bits()
            );
            assert_eq!(
                a.fluid_drag_work.to_bits(),
                b.resistance.fluid_drag_work.impulse_work.to_bits()
            );
            assert_eq!(
                a.skeleton_work.to_bits(),
                b.resistance.skeleton_work.impulse_work.to_bits()
            );
            work_balance(b);
        }
    });
}

#[test]
fn unchanged_laws_preserve_held_native_population_and_impulse_work_bits() {
    parity::<D3q19, 19>(CollisionModel3::Bgk { tau: 0.9 });
    parity::<D3q27, 27>(bgk());
}

#[test]
fn invalid_laws_work_budget_and_late_cancellation_preserve_accepted_history() {
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
        let mut expected = channel::<D3q27, 27>(cx, bgk());
        let laws = next_laws(expected.laws(), 1);
        let receipt = expected
            .advance_serial(cx, None, laws.clone(), 1.0)
            .unwrap();
        assert!(receipt.resistance.absolute_relative_quadrature_remainder > 0.0);
        let mut owner = channel::<D3q27, 27>(cx, bgk());
        let before = snapshot(&owner);
        let old_laws = owner.laws().to_vec();
        let mut invalid = laws.clone();
        invalid[0].linear.principal = [0.1; 3]; // newly active mixed seam
        assert!(matches!(
            owner.advance_serial(cx, None, invalid, 1.0),
            Err(GridError3::MaterialField { cell: 0, .. })
        ));
        let mut invalid = laws.clone();
        invalid[139].quadratic = f64::NAN;
        assert!(matches!(
            owner.advance_serial(cx, None, invalid, 1.0),
            Err(GridError3::MaterialField { cell: 139, .. })
        ));
        assert!(matches!(
            owner.advance_serial(cx, None, laws[..139].to_vec(), 1.0),
            Err(GridError3::Input(_))
        ));
        for limit in [0.0, -1.0, f64::NAN] {
            assert!(matches!(
                owner.advance_serial(cx, None, laws.clone(), limit),
                Err(GridError3::Input(_))
            ));
            assert_eq!(snapshot(&owner), before);
            assert_eq!(owner.laws(), old_laws);
            assert_eq!(owner.last_endpoint_exchange(), None);
        }
        let runner = CancelAfterStream {
            pool: TilePool::for_host(2, SEED),
            calls: AtomicUsize::new(0),
        };
        assert!(matches!(
            owner.advance_pooled(
                &runner,
                &CancelGate::new_clock_free(),
                None,
                laws.clone(),
                1.0
            ),
            Err(GridError3::Cancelled)
        ));
        assert_eq!(snapshot(&owner), before);
        assert_eq!(owner.laws(), old_laws);
        let retry = owner.advance_serial(cx, None, laws, 1.0).unwrap();
        assert_eq!(retry, receipt);
        assert_eq!(snapshot(&owner), snapshot(&expected));
        assert_eq!(owner.last_endpoint_exchange(), Some(receipt));
    });
}
