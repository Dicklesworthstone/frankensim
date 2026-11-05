use super::*;
use crate::d3q19::CollisionModel3;
use crate::d3q27::{CollisionConfig27, CollisionModel27, D3q27, Equilibrium27, Rates27};
use crate::lattice3::{D3q19, GridBudget3};
use crate::material::plan_scaling_from_material;

fn with_gate<T>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> T) -> T {
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        f(&Cx::new(
            gate,
            arena,
            fs_exec::StreamKey {
                seed: 0x0B03_4E0,
                kernel_id: 3,
                tile: 0,
                iteration: 0,
            },
            fs_exec::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        ))
    })
}
fn with_cx<T>(f: impl FnOnce(&Cx<'_>) -> T) -> T {
    with_gate(&CancelGate::new_clock_free(), f)
}
fn fluid(mu: f64) -> FluidTransportState {
    FluidTransportState::try_new(1000.0, 4000.0, 0.5, mu).unwrap()
}
fn model27(model: CollisionModel27) -> CollisionConfig27 {
    CollisionConfig27 {
        model,
        maximum_mach: 0.2,
    }
}
fn bgk() -> CollisionConfig27 {
    model27(CollisionModel27::Bgk {
        tau: 0.9,
        equilibrium: Equilibrium27::Factorized,
    })
}
#[allow(clippy::too_many_arguments)]
fn make<L: EquilibriumLattice3<Q>, const Q: usize>(
    cx: &Cx<'_>,
    dims: [usize; 3],
    collision: L::Collision,
    force: [f64; 3],
    faces: [FaceBoundary3; 6],
    solid: impl Fn([usize; 3]) -> bool,
    velocity: impl Fn([usize; 3]) -> [f64; 3],
) -> OpenBoundaryGrid3<L, Q> {
    let periodic = core::array::from_fn(|a| matches!(faces[2 * a], FaceBoundary3::Periodic));
    let grid = Grid3::<L, Q>::try_new_equilibrium(
        dims,
        collision,
        force,
        periodic,
        GridBudget3 {
            maximum_cells: dims.iter().product(),
            maximum_mach: 0.2,
        },
        solid,
        |p| (1.3, velocity(p)),
    )
    .unwrap();
    OpenBoundaryGrid3::try_new(grid, cx, BoundarySpec3::new(faces)).unwrap()
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() < tolerance, "{a:.17e} != {b:.17e}");
}
fn bits<L: EquilibriumLattice3<Q>, const Q: usize>(flow: &OpenBoundaryGrid3<L, Q>) -> Vec<u64> {
    flow.populations()
        .iter()
        .flatten()
        .map(|f| f.to_bits())
        .collect()
}

fn uniform_source<L: EquilibriumLattice3<Q> + ViscousLattice3<Q>, const Q: usize>(
    collision: L::Collision,
) {
    with_cx(|cx| {
        let material = fluid(0.01);
        let reference = plan_scaling_from_material(&material, 0.001, 0.01, 0.005).unwrap();
        let velocity = [0.021, -0.017, 0.013];
        for axis in 0..3 {
            let mut faces = [FaceBoundary3::Periodic; 6];
            faces[2 * axis] = FaceBoundary3::Velocity { velocity };
            faces[2 * axis + 1] = FaceBoundary3::Pressure { density: 1.3 };
            let mut flow = make::<L, Q>(
                cx,
                [4, 5, 6],
                collision,
                [0.00013, -0.00009, 0.00007],
                faces,
                |_| false,
                |_| velocity,
            );
            let mut workspace = flow.rheology_workspace().unwrap();
            let mut queries = 0;
            let (report, exchange) = flow
                .step_serial_with_material_rheology_workspace(
                    cx,
                    &reference,
                    &mut workspace,
                    [0.0; 3],
                    |point, _| {
                        queries += 1;
                        near(point.shear_rate_per_s, 0.0, 2e-12);
                        Ok(material)
                    },
                )
                .unwrap();
            assert_eq!(queries, 120);
            assert!(exchange.reconstructed_cells > 0);
            for sample in report.samples.iter().flatten() {
                near(sample.shear_relaxation_time(), 0.8, 2e-15);
            }
        }
    });
}

#[test]
fn uniform_physical_flow_has_no_false_reservoir_shear_in_native_families() {
    for collision in [
        CollisionModel3::Bgk { tau: 0.9 },
        CollisionModel3::Trt {
            tau: 0.9,
            magic: 0.19,
        },
    ] {
        uniform_source::<D3q19, 19>(collision);
    }
    for equilibrium in [Equilibrium27::Quadratic, Equilibrium27::Factorized] {
        uniform_source::<D3q27, 27>(model27(CollisionModel27::Bgk {
            tau: 0.9,
            equilibrium,
        }));
        uniform_source::<D3q27, 27>(model27(CollisionModel27::Trt {
            tau: 0.9,
            magic: 0.19,
            equilibrium,
        }));
    }
    let rates = Rates27 { omega: [1.1; 10] };
    uniform_source::<D3q27, 27>(model27(CollisionModel27::CentralMoment { rates }));
    uniform_source::<D3q27, 27>(model27(CollisionModel27::Cumulant { rates }));
}

#[test]
fn quadratic_open_face_shear_uses_si_rate_and_drives_actual_collision() {
    with_cx(|cx| {
        let material = fluid(0.01);
        let reference = plan_scaling_from_material(&material, 0.001, 0.01, 0.005).unwrap();
        let dims = [5, 6, 7];
        for axis in 0..3 {
            let profile = |p: [usize; 3]| {
                let x = p[axis] as f64;
                let mut u = [0.0; 3];
                u[(axis + 1) % 3] = 0.011 + 0.002 * x + 0.0003 * x * x;
                u
            };
            let mut faces = [FaceBoundary3::Periodic; 6];
            faces[2 * axis] = FaceBoundary3::Velocity {
                velocity: profile([0; 3]),
            };
            let mut high = [0; 3];
            high[axis] = dims[axis] - 1;
            faces[2 * axis + 1] = FaceBoundary3::Velocity {
                velocity: profile(high),
            };
            let mut flow = make::<D3q27, 27>(cx, dims, bgk(), [0.0; 3], faces, |_| false, profile);
            let mut control =
                make::<D3q27, 27>(cx, dims, bgk(), [0.0; 3], faces, |_| false, profile);
            let mut workspace = flow.rheology_workspace().unwrap();
            let (report, _) = flow
                .step_serial_with_material_rheology_workspace(
                    cx,
                    &reference,
                    &mut workspace,
                    [0.0; 3],
                    |point, _| {
                        let expected = (0.002 + 0.0006 * point.coordinates[axis] as f64) / 0.01;
                        near(point.shear_rate_per_s, expected, 2e-12);
                        Ok(fluid(0.01 / (1.0 + point.shear_rate_per_s)))
                    },
                )
                .unwrap();
            let expected: Vec<_> = (0..dims.iter().product())
                .map(|cell| {
                    let x = coordinates(cell, dims)[axis] as f64;
                    let rate = (0.002 + 0.0006 * x) / 0.01;
                    let tau = 0.5 + 0.3 / (1.0 + rate);
                    near(
                        report.samples[cell]
                            .as_ref()
                            .unwrap()
                            .shear_relaxation_time(),
                        tau,
                        1e-13,
                    );
                    D3q27::with_shear_relaxation(bgk(), tau).unwrap()
                })
                .collect();
            control.step_serial(cx, Some(&expected), [0.0; 3]).unwrap();
            for (a, b) in flow
                .populations()
                .iter()
                .flatten()
                .zip(control.populations().iter().flatten())
            {
                near(*a, *b, 5e-14);
            }
        }
    });
}

#[test]
fn pressure_outlets_retain_tangential_halfway_wall_shear_at_seams() {
    with_cx(|cx| {
        let material = fluid(0.01);
        let reference = plan_scaling_from_material(&material, 0.001, 0.01, 0.005).unwrap();
        let faces = [
            FaceBoundary3::Pressure { density: 1.3 },
            FaceBoundary3::Pressure { density: 1.3 },
            FaceBoundary3::Wall { velocity: [0.0; 3] },
            FaceBoundary3::Wall { velocity: [0.0; 3] },
            FaceBoundary3::Periodic,
            FaceBoundary3::Periodic,
        ];
        let mut flow = make::<D3q27, 27>(
            cx,
            [5, 6, 3],
            bgk(),
            [0.0; 3],
            faces,
            |_| false,
            |p| {
                let y = p[1] as f64 + 0.5;
                [0.003 * y * (6.0 - y), 0.0, 0.0]
            },
        );
        let mut workspace = flow.rheology_workspace().unwrap();
        flow.step_serial_with_material_rheology_workspace(
            cx,
            &reference,
            &mut workspace,
            [0.0; 3],
            |point, _| {
                let expected = (0.003 * (5.0 - 2.0 * point.coordinates[1] as f64)).abs() / 0.01;
                near(point.shear_rate_per_s, expected, 2e-12);
                Ok(fluid(0.01 / (1.0 + point.shear_rate_per_s)))
            },
        )
        .unwrap();
    });
}

#[test]
fn blocked_second_donor_falls_back_to_linear_and_pressure_normal_is_zero() {
    with_cx(|cx| {
        let mut faces = [FaceBoundary3::Periodic; 6];
        faces[0] = FaceBoundary3::Velocity { velocity: [0.0; 3] };
        faces[1] = FaceBoundary3::Pressure { density: 1.3 };
        let flow = make::<D3q27, 27>(
            cx,
            [5, 3, 3],
            bgk(),
            [0.0; 3],
            faces,
            |p| p == [2, 1, 1],
            |_| [0.0; 3],
        );
        let velocity: Vec<_> = (0..45)
            .map(|cell| {
                let x = coordinates(cell, [5, 3, 3])[0] as f64;
                [0.0, 0.01 * x * x, 0.0]
            })
            .collect();
        near(flow.rheology_gradient(20, &velocity)[1][0], 0.01, 1e-15);
        near(flow.rheology_gradient(0, &velocity)[1][0], 0.0, 1e-15);
        near(flow.rheology_gradient(24, &velocity)[1][0], 0.0, 1e-15);
    });
}

#[test]
fn pooled_source_epochs_match_serial_and_reuse_report_buffers() {
    with_cx(|cx| {
        let material = fluid(0.01);
        let reference = plan_scaling_from_material(&material, 0.001, 0.01, 0.005).unwrap();
        let mut faces = [FaceBoundary3::Periodic; 6];
        faces[0] = FaceBoundary3::Velocity {
            velocity: [0.02, 0.0, 0.0],
        };
        faces[1] = FaceBoundary3::Pressure { density: 1.3 };
        let factory = || {
            make::<D3q27, 27>(
                cx,
                [8, 5, 4],
                bgk(),
                [2e-5, 0.0, 0.0],
                faces,
                |_| false,
                |p| [0.02, 0.001 * p[0] as f64, 0.0],
            )
        };
        let mut serial = factory();
        let mut pooled = factory();
        let mut serial_workspace = serial.rheology_workspace().unwrap();
        let mut pooled_workspace = pooled.rheology_workspace().unwrap();
        let pool = fs_exec::TilePool::for_host(3, 0x0B03_4E0);
        let gate = CancelGate::new_clock_free();
        let mut first_buffers = None;
        for step in 0..3 {
            let source = |point: RheologyPoint3, _: &FluidTransportState| {
                Ok(fluid(0.01 / (1.0 + point.shear_rate_per_s + point.time_s)))
            };
            let (a, exchange_a) = serial
                .step_serial_with_material_rheology_workspace(
                    cx,
                    &reference,
                    &mut serial_workspace,
                    [0.0; 3],
                    source,
                )
                .unwrap();
            let (b, exchange_b) = pooled
                .step_pooled_with_material_rheology_workspace(
                    &pool,
                    &gate,
                    &reference,
                    &mut pooled_workspace,
                    [0.0; 3],
                    source,
                )
                .unwrap();
            assert_eq!(bits(&serial), bits(&pooled));
            assert_eq!(a.sampled_at_step, step);
            assert_eq!(a.collisions, b.collisions);
            assert_eq!(exchange_a, exchange_b);
            for (a, b) in a.samples.iter().flatten().zip(b.samples.iter().flatten()) {
                assert_eq!(a.point(), b.point());
                assert_eq!(a.transport(), b.transport());
            }
            let buffers = (a.samples.as_ptr(), a.collisions.as_ptr());
            if step == 0 {
                first_buffers = Some(buffers);
            }
            if step == 2 {
                assert_eq!(Some(buffers), first_buffers);
            }
        }
    });
}

#[test]
fn failed_native_admission_and_source_cancellation_preserve_accepted_epoch() {
    let gate = CancelGate::new_clock_free();
    with_gate(&gate, |cx| {
        let material = fluid(0.01);
        let reference = plan_scaling_from_material(&material, 0.001, 0.01, 0.005).unwrap();
        let mut faces = [FaceBoundary3::Periodic; 6];
        faces[0] = FaceBoundary3::Velocity {
            velocity: [0.01, 0.0, 0.0],
        };
        faces[1] = FaceBoundary3::Pressure { density: 1.3 };
        let mut flow = make::<D3q27, 27>(
            cx,
            [4, 3, 3],
            bgk(),
            [0.0; 3],
            faces,
            |_| false,
            |_| [0.01, 0.0, 0.0],
        );
        let mut workspace = flow.rheology_workspace().unwrap();
        flow.step_serial_with_material_rheology_workspace(
            cx,
            &reference,
            &mut workspace,
            [0.0; 3],
            |_, _| Ok(material),
        )
        .unwrap();
        let before = bits(&flow);
        let report_address = workspace.report().unwrap().samples.as_ptr();
        let failure = flow.step_serial_with_material_rheology_workspace(
            cx,
            &reference,
            &mut workspace,
            [f64::NAN, 0.0, 0.0],
            |_, _| Ok(fluid(0.02)),
        );
        assert!(matches!(failure, Err(RheologyError3::Grid(_))));
        assert_eq!(bits(&flow), before);
        assert_eq!(flow.completed_steps(), 1);
        assert_eq!(workspace.report().unwrap().sampled_at_step, 0);
        assert_eq!(workspace.report().unwrap().samples.as_ptr(), report_address);
        let failure = flow.step_serial_with_material_rheology_workspace(
            cx,
            &reference,
            &mut workspace,
            [0.0; 3],
            |point, _| {
                if point.cell == 3 {
                    FluidTransportState::try_new(1000.0, 4000.0, 0.5, -1.0)
                } else {
                    Ok(fluid(0.02))
                }
            },
        );
        assert!(matches!(
            failure,
            Err(RheologyError3::Source { cell: 3, .. })
        ));
        assert_eq!(bits(&flow), before);
        assert_eq!(workspace.report().unwrap().samples.as_ptr(), report_address);
        let failure = flow.step_serial_with_material_rheology_workspace(
            cx,
            &reference,
            &mut workspace,
            [0.0; 3],
            |point, _| {
                if point.cell == 3 {
                    gate.request();
                }
                Ok(fluid(0.02))
            },
        );
        assert!(matches!(
            failure,
            Err(RheologyError3::Grid(GridError3::Cancelled))
        ));
        assert_eq!(bits(&flow), before);
        assert_eq!(workspace.report().unwrap().sampled_at_step, 0);
        assert_eq!(workspace.report().unwrap().samples.as_ptr(), report_address);
        for sample in workspace.report().unwrap().samples.iter().flatten() {
            assert_eq!(sample.transport().dynamic_viscosity_pa_s(), 0.01);
        }
    });
}
