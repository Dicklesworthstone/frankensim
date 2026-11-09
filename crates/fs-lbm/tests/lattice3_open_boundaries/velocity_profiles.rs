use super::*;
use fs_lbm::d3q19::{BoundaryProfileError3, InletFlux3};
use fs_lbm::lattice3::open_boundaries::forcing::ForcedOpenBoundaryGrid3;

const DIMS: [usize; 3] = [5, 6, 7];

fn fixture<L: EquilibriumLattice3<Q>, const Q: usize>(
    cx: &Cx<'_>,
    collision: L::Collision,
    axis: usize,
    upper: bool,
    force: [f64; 3],
) -> OpenBoundaryGrid3<L, Q> {
    let mut periodic = [true; 3];
    periodic[axis] = false;
    let mut uniform = [0.0; 3];
    uniform[axis] = if upper { -0.02 } else { 0.02 };
    let grid = Grid3::try_new_equilibrium(
        DIMS,
        collision,
        force,
        periodic,
        budget(),
        |_| false,
        |_| (1.1, uniform),
    )
    .unwrap();
    let mut faces = [FaceBoundary3::Periodic; 6];
    faces[2 * axis + usize::from(upper)] = FaceBoundary3::Velocity { velocity: uniform };
    faces[2 * axis + usize::from(!upper)] = FaceBoundary3::Pressure { density: 1.1 };
    OpenBoundaryGrid3::try_new(grid, cx, BoundarySpec3::new(faces)).unwrap()
}

fn trace(p: [usize; 3], axis: usize, upper: bool, epoch: usize) -> [f64; 3] {
    let (a, b) = ((axis + 1) % 3, (axis + 2) % 3);
    let mut value = [0.0; 3];
    value[axis] = (if upper { -1.0 } else { 1.0 })
        * (0.024 + 0.003 * epoch as f64 + 0.0007 * p[a] as f64 + 0.0004 * p[b] as f64);
    value[a] = 0.001 * (p[b] as f64 - 2.0);
    value[b] = -0.0006 * (p[a] as f64 - 3.0);
    value
}

fn resolved<L: EquilibriumLattice3<Q>, const Q: usize>(collision: L::Collision) {
    with_cx(|cx| {
        let force = [0.00013, -0.00009, 0.00007];
        for axis in 0..3 {
            for upper in [false, true] {
                let face = Face3::ALL[2 * axis + usize::from(upper)];
                let sign = if upper { 1.0 } else { -1.0 };
                let mut flow = fixture::<L, Q>(cx, collision, axis, upper, force);
                let mut calls = 0;
                for epoch in 0..3 {
                    let old_state = bits(&flow);
                    let mean = -0.03 - 0.003 * epoch as f64;
                    let report = flow
                        .set_velocity_profile_with_speed_limit(
                            face,
                            InletFlux3::MeanNormal(mean),
                            0.06,
                            &CancelGate::new_clock_free(),
                            |p| {
                                calls += 1;
                                Ok(trace(p, axis, upper, epoch))
                            },
                        )
                        .unwrap();
                    assert_eq!(bits(&flow), old_state);
                    assert_eq!(flow.completed_steps(), epoch);
                    let area = DIMS.iter().product::<usize>() / DIMS[axis];
                    assert_eq!(report.sampled_cells, area);
                    assert_eq!(calls, (epoch + 1) * area);
                    near(sign * report.prescribed_mean_velocity_lu[axis], mean, 2e-17);
                    let profile = flow.velocity_profile(face).unwrap().to_vec();
                    let before = totals(&flow);
                    let receipt = flow.step_serial(cx, None, [0.0; 3]).unwrap();
                    let after = totals(&flow);
                    assert_eq!(calls, (epoch + 1) * area);
                    near(after.0 - before.0, receipt.mass_change, 6e-12);
                    for a in 0..3 {
                        near(
                            after.1[a] - before.1[a],
                            flow.populations().len() as f64 * force[a] + receipt.impulse[a],
                            6e-12,
                        );
                    }
                    let mut index = 0;
                    for cell in 0..flow.populations().len() {
                        let coordinate = cell / DIMS[..axis].iter().product::<usize>() % DIMS[axis];
                        if coordinate != if upper { DIMS[axis] - 1 } else { 0 } {
                            continue;
                        }
                        let actual = flow.macroscopic(cell).unwrap().1;
                        for a in 0..3 {
                            near(actual[a], profile[index][a], 5e-13);
                        }
                        index += 1;
                    }
                    assert_eq!(index, area);
                }
            }
        }
    });
}

#[test]
fn resolved_time_updated_vectors_drive_all_faces_in_native_collision_families() {
    resolved::<D3q19, 19>(CollisionModel3::Bgk { tau: 0.9 });
    resolved::<D3q19, 19>(CollisionModel3::Trt {
        tau: 0.9,
        magic: 3.0 / 16.0,
    });
    resolved::<D3q19, 19>(CollisionModel3::CentralMomentBulk {
        shear_rate: 1.1,
        bulk_rate: 0.8,
        higher_order_rate: 1.2,
    });
    for equilibrium in [Equilibrium27::Quadratic, Equilibrium27::Factorized] {
        resolved::<D3q27, 27>(config(CollisionModel27::Bgk {
            tau: 0.9,
            equilibrium,
        }));
    }
    resolved::<D3q27, 27>(config(CollisionModel27::Cumulant {
        rates: Rates27 { omega: [1.1; 10] },
    }));
}

#[test]
fn uniform_profile_and_clear_keep_the_original_trajectory_bits() {
    with_cx(|cx| {
        let mut ordinary = fixture::<D3q27, 27>(cx, bgk(), 0, false, [0.0001, -0.00003, 0.00002]);
        let mut profiled = fixture::<D3q27, 27>(cx, bgk(), 0, false, ordinary.force_density());
        let report = profiled
            .set_velocity_profile(
                Face3::XMin,
                InletFlux3::MeanNormal(-0.02),
                &CancelGate::new_clock_free(),
                |_| Ok([0.02, 0.0, 0.0]),
            )
            .unwrap();
        assert_eq!(report.normal_velocity_correction_lu, 0.0);
        for _ in 0..3 {
            let a = ordinary.step_serial(cx, None, [0.0; 3]).unwrap();
            let b = profiled.step_serial(cx, None, [0.0; 3]).unwrap();
            assert_eq!(a, b);
            assert_eq!(bits(&ordinary), bits(&profiled));
        }
        profiled.clear_velocity_profile(Face3::XMin);
        assert!(profiled.velocity_profile(Face3::XMin).is_none());
        ordinary.step_serial(cx, None, [0.0; 3]).unwrap();
        profiled.step_serial(cx, None, [0.0; 3]).unwrap();
        assert_eq!(bits(&ordinary), bits(&profiled));
    });
}

#[test]
fn refused_profile_samples_corrections_and_late_cancellation_keep_the_previous_trace() {
    with_cx(|cx| {
        let mut flow = fixture::<D3q27, 27>(cx, bgk(), 0, false, [0.0; 3]);
        flow.set_velocity_profile(
            Face3::XMin,
            InletFlux3::Unconstrained,
            &CancelGate::new_clock_free(),
            |p| Ok(trace(p, 0, false, 0)),
        )
        .unwrap();
        let old = flow.velocity_profile(Face3::XMin).unwrap().to_vec();
        let before = bits(&flow);
        for value in [[0.12, 0.0, 0.0], [f64::NAN, 0.0, 0.0]] {
            assert!(matches!(
                flow.set_velocity_profile(
                    Face3::XMin,
                    InletFlux3::Unconstrained,
                    &CancelGate::new_clock_free(),
                    |_| Ok(value)
                ),
                Err(BoundaryProfileError3::InvalidVelocity { .. })
            ));
        }
        assert!(matches!(
            flow.set_velocity_profile(
                Face3::XMin,
                InletFlux3::MeanNormal(-0.115),
                &CancelGate::new_clock_free(),
                |_| Ok([0.04, 0.02, 0.0])
            ),
            Err(BoundaryProfileError3::InvalidVelocity { .. })
        ));
        assert!(matches!(
            flow.set_velocity_profile_with_speed_limit(
                Face3::XMin,
                InletFlux3::MeanNormal(-0.05),
                0.03,
                &CancelGate::new_clock_free(),
                |_| Ok([0.02, 0.0, 0.0])
            ),
            Err(BoundaryProfileError3::InvalidVelocity { .. })
        ));
        let gate = CancelGate::new_clock_free();
        let mut samples = 0;
        assert!(matches!(
            flow.set_velocity_profile(Face3::XMin, InletFlux3::Unconstrained, &gate, |p| {
                samples += 1;
                if samples == DIMS[1] * DIMS[2] {
                    gate.request();
                }
                Ok(trace(p, 0, false, 1))
            }),
            Err(BoundaryProfileError3::Cancelled)
        ));
        assert_eq!(samples, DIMS[1] * DIMS[2]);
        assert!(matches!(
            flow.set_velocity_profile(
                Face3::XMin,
                InletFlux3::Unconstrained,
                &CancelGate::new_clock_free(),
                |_| Err(BoundaryProfileError3::Source("sample refused".into()))
            ),
            Err(BoundaryProfileError3::Source(_))
        ));
        assert_eq!(flow.velocity_profile(Face3::XMin).unwrap(), old);
        assert_eq!(bits(&flow), before);
        assert_eq!(flow.completed_steps(), 0);
    });
}

#[test]
fn profile_survives_wall_updates_and_keeps_wall_owned_seam_populations() {
    with_cx(|cx| {
        let make = || {
            let grid = Grid3::<D3q27, 27>::try_new_equilibrium(
                [6, 5, 4],
                bgk(),
                [0.0; 3],
                [false, false, true],
                budget(),
                |_| false,
                |_| (1.1, [0.025, 0.0, 0.0]),
            )
            .unwrap();
            OpenBoundaryGrid3::try_new(
                grid,
                cx,
                BoundarySpec3::new([
                    FaceBoundary3::Velocity {
                        velocity: [0.025, 0.0, 0.0],
                    },
                    FaceBoundary3::Pressure { density: 1.1 },
                    FaceBoundary3::Wall {
                        velocity: [0.025, 0.0, 0.0],
                    },
                    FaceBoundary3::Wall {
                        velocity: [0.025, 0.0, 0.0],
                    },
                    FaceBoundary3::Periodic,
                    FaceBoundary3::Periodic,
                ]),
            )
            .unwrap()
        };
        let mut flow = make();
        let mut control = make();
        flow.set_velocity_profile(
            Face3::XMin,
            InletFlux3::Unconstrained,
            &CancelGate::new_clock_free(),
            |p| Ok([0.02 + 0.003 * p[1] as f64, 0.001 * p[2] as f64, 0.0]),
        )
        .unwrap();
        let old = flow.velocity_profile(Face3::XMin).unwrap().to_vec();
        assert!(matches!(
            flow.set_velocity_profile(
                Face3::XMin,
                InletFlux3::MeanNormal(-0.03),
                &CancelGate::new_clock_free(),
                |_| Ok([0.03, 0.0, 0.0])
            ),
            Err(BoundaryProfileError3::MixedWallSeam)
        ));
        let mut faces = core::array::from_fn(|a| flow.boundaries().face(Face3::ALL[a]));
        faces[3] = FaceBoundary3::Wall {
            velocity: [0.05, 0.0, 0.0],
        };
        for grid in [&mut flow, &mut control] {
            grid.rebind_open_boundaries(BoundarySpec3::new(faces), &CancelGate::new_clock_free())
                .unwrap();
        }
        assert_eq!(flow.velocity_profile(Face3::XMin).unwrap(), old);
        flow.step_serial(cx, None, [0.0; 3]).unwrap();
        control.step_serial(cx, None, [0.0; 3]).unwrap();
        assert_eq!(flow.last_wall_exchange(), control.last_wall_exchange());
        for y in [0, 4] {
            let cell = address([0, y, 1], [6, 5, 4]);
            for q in 1..27 {
                if (y == 0 && D3q27::E[q][1] > 0) || (y == 4 && D3q27::E[q][1] < 0) {
                    assert_eq!(
                        flow.populations()[cell][q].to_bits(),
                        control.populations()[cell][q].to_bits()
                    );
                }
            }
        }
        faces[0] = FaceBoundary3::Velocity {
            velocity: [0.03, 0.0, 0.0],
        };
        assert!(
            flow.rebind_open_boundaries(BoundarySpec3::new(faces), &CancelGate::new_clock_free())
                .is_err()
        );
        assert_eq!(flow.velocity_profile(Face3::XMin).unwrap(), old);
        flow.clear_velocity_profile(Face3::XMin);
        flow.rebind_open_boundaries(BoundarySpec3::new(faces), &CancelGate::new_clock_free())
            .unwrap();
    });
}

#[test]
fn spatial_force_history_uses_the_same_resolved_velocity_profile() {
    with_cx(|cx| {
        let original = fixture::<D3q27, 27>(cx, bgk(), 2, true, [0.0; 3]);
        let cells = original.populations().len();
        let mut flow = ForcedOpenBoundaryGrid3::try_new_preserving_velocity(
            original,
            cx,
            (0..cells)
                .map(|cell| [0.00001 * (cell % 5) as f64, -0.00003, 0.00004])
                .collect(),
        )
        .unwrap();
        for epoch in 0..2 {
            flow.set_velocity_profile_with_speed_limit(
                Face3::ZMax,
                InletFlux3::MeanNormal(-0.035),
                0.06,
                &CancelGate::new_clock_free(),
                |p| Ok(trace(p, 2, true, epoch)),
            )
            .unwrap();
            let profile = flow.velocity_profile(Face3::ZMax).unwrap().to_vec();
            flow.advance_serial(
                cx,
                None,
                (0..cells)
                    .map(|cell| [0.00002, -0.00001 * (cell % 7) as f64, -0.00002])
                    .collect(),
                [0.0; 3],
            )
            .unwrap();
            for (index, expected) in profile.iter().enumerate() {
                let cell = (DIMS[2] - 1) * DIMS[0] * DIMS[1] + index;
                let actual = flow.macroscopic(cell).unwrap().1;
                for a in 0..3 {
                    near(actual[a], expected[a], 5e-13);
                }
            }
            assert_eq!(flow.completed_steps(), epoch + 1);
        }
        flow.clear_velocity_profile(Face3::ZMax);
        assert!(flow.velocity_profile(Face3::ZMax).is_none());
    });
}
