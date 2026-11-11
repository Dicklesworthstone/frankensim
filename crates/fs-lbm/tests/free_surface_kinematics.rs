//! Native sparse liquid curl: actual force frame, free interfaces and solids.

use fs_exec::{CancelGate, Cx};
use fs_lbm::d3q19::freesurface3::sgs::{FreeSurfaceSgsError3, FreeSurfaceSgsWorkspace3};
use fs_lbm::d3q19::freesurface3::{
    ContactModel3, FreeSurface3, FreeSurfaceInitialCell3, FreeSurfaceKinematicsError3,
    FreeSurfaceKinematicsWorkspace3, PlanarWall3, WetGradientStencil3,
};
use fs_lbm::kinematics::{FlowKinematics3, KinematicsError};
use fs_lbm::sgs3::{SgsConfig, SgsModel, SgsUse, eddy_viscosity};

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

fn liquid(
    dims: [usize; 3],
    occupied: impl Fn([usize; 3]) -> bool,
    velocity: impl Fn([usize; 3]) -> [f64; 3],
    solids: &[(u32, u32, u32)],
) -> FreeSurface3 {
    let mut initial = Vec::new();
    for z in 0..dims[2] {
        for y in 0..dims[1] {
            for x in 0..dims[0] {
                let p = [x, y, z];
                if occupied(p) {
                    initial.push(FreeSurfaceInitialCell3 {
                        coordinates: p,
                        fill_fraction: if (x + y + z) % 3 == 0 { 0.6 } else { 1.0 },
                        density: 0.9 + 0.01 * (x + 2 * y + 3 * z) as f64,
                        velocity: velocity(p),
                    });
                }
            }
        }
    }
    with_cx(|cx| {
        FreeSurface3::from_cells(
            dims,
            0.8,
            [0.0006, -0.0002, 0.0004],
            0.0,
            ContactModel3::Neutral,
            solids,
            &initial,
            cx,
        )
        .unwrap()
    })
}

fn workspace(owner: &FreeSurface3) -> FreeSurfaceKinematicsWorkspace3 {
    let tiles = owner.grid().active_tiles();
    let bytes = FreeSurfaceKinematicsWorkspace3::required_bytes(tiles).unwrap();
    FreeSurfaceKinematicsWorkspace3::new(tiles, bytes).unwrap()
}

fn near(actual: f64, expected: f64) {
    assert!((actual - expected).abs() < 5e-14, "{actual} != {expected}");
}

fn same_gradient(actual: FlowKinematics3, expected: [[f64; 3]; 3]) {
    for (a, b) in actual
        .gradient()
        .iter()
        .flatten()
        .zip(expected.iter().flatten())
    {
        near(*a, *b);
    }
    let expected = FlowKinematics3::from_gradient(expected).unwrap();
    for (a, b) in actual.vorticity().into_iter().zip(expected.vorticity()) {
        near(a, b);
    }
    for (a, b) in actual
        .vortex_stretching()
        .into_iter()
        .zip(expected.vortex_stretching())
    {
        near(a, b);
    }
    near(actual.enstrophy_density(), expected.enstrophy_density());
    near(actual.q_criterion(), expected.q_criterion());
}

#[test]
fn native_sparse_full_gradient_keeps_curl_strain_and_unequal_density_force_frame() {
    let gate = CancelGate::new_clock_free();
    let gradient = [
        [-0.001, -0.004, 0.002],
        [0.004, -0.001, -0.003],
        [0.001, 0.003, 0.002],
    ];
    let owner = liquid(
        [8; 3],
        |_| true,
        |p| {
            core::array::from_fn(|a| {
                0.01 + (0..3).map(|b| gradient[a][b] * p[b] as f64).sum::<f64>()
            })
        },
        &[],
    );
    let bits = owner.grid().state_bits();
    let mut cache = workspace(&owner);
    let field = owner.kinematics(&gate, &mut cache).unwrap();
    let mut count = 0;
    field
        .visit_samples(|_, sample| {
            same_gradient(sample.unwrap().kinematics, gradient);
            count += 1;
        })
        .unwrap();
    assert_eq!(count, 512);
    assert_eq!(owner.grid().state_bits(), bits);
    assert!(field.sample_at(usize::MAX, 0, 0).unwrap().is_none());
    let sample = field.sample_at(3, 3, 3).unwrap().unwrap();
    assert_eq!(
        sample.stencils,
        [WetGradientStencil3::CenteredSecondOrder; 3]
    );
    same_gradient(
        sample.kinematics.to_physical(0.25).unwrap(),
        gradient.map(|row| row.map(|g| 4.0 * g)),
    );
}

#[test]
fn free_interfaces_use_contiguous_liquid_quadratics_in_every_orientation() {
    let gate = CancelGate::new_clock_free();
    for axis in 0..3 {
        let flow = (axis + 1) % 3;
        let owner = liquid(
            [8; 3],
            |p| (2..6).contains(&p[axis]),
            |p| {
                let mut u = [0.01; 3];
                u[flow] += 0.0005 * (p[axis] as f64 + 0.5).powi(2);
                u
            },
            &[],
        );
        let mut cache = workspace(&owner);
        let field = owner.kinematics(&gate, &mut cache).unwrap();
        for coordinate in 0..8 {
            let mut p = [3; 3];
            p[axis] = coordinate;
            let actual = field.sample_at(p[0], p[1], p[2]).unwrap();
            if !(2..6).contains(&coordinate) {
                assert!(actual.is_none());
                continue;
            }
            let actual = actual.unwrap();
            let mut gradient = [[0.0; 3]; 3];
            gradient[flow][axis] = 0.001 * (coordinate as f64 + 0.5);
            same_gradient(actual.kinematics, gradient);
            assert_eq!(
                actual.stencils[axis],
                match coordinate {
                    2 => WetGradientStencil3::ForwardSecondOrder,
                    5 => WetGradientStencil3::BackwardSecondOrder,
                    _ => WetGradientStencil3::CenteredSecondOrder,
                }
            );
        }
    }
}

#[test]
fn unresolved_sheet_refuses_normal_derivatives_and_does_not_bridge_gas() {
    let gate = CancelGate::new_clock_free();
    let owner = liquid(
        [8; 3],
        |p| p[1] == 3 || p[1] == 5,
        |p| [0.002 * p[1] as f64, 0.0, 0.0],
        &[],
    );
    let mut cache = workspace(&owner);
    let field = owner.kinematics(&gate, &mut cache).unwrap();
    assert_eq!(
        field.sample_at(3, 3, 3),
        Err(FreeSurfaceKinematicsError3::MissingDerivative {
            cell: [3, 3, 3],
            axis: 1
        })
    );
    let owner = liquid(
        [8; 3],
        |p| p[1] == 3 || p[1] == 4,
        |p| [0.002 * p[1] as f64, 0.0, 0.0],
        &[],
    );
    let field = owner.kinematics(&gate, &mut cache).unwrap();
    for (y, stencil) in [
        (3, WetGradientStencil3::ForwardFirstOrder),
        (4, WetGradientStencil3::BackwardFirstOrder),
    ] {
        let sample = field.sample_at(3, y, 3).unwrap().unwrap();
        same_gradient(sample.kinematics, [[0.0, 0.002, 0.0], [0.0; 3], [0.0; 3]]);
        assert_eq!(sample.stencils[1], stencil);
    }
}

#[test]
fn solid_tiles_and_current_piston_slab_never_enter_the_liquid_stencil() {
    let gate = CancelGate::new_clock_free();
    let owner = liquid(
        [8; 3],
        |p| p.iter().any(|&v| v >= 4),
        |p| [0.01, 0.003 * p[0] as f64, 0.0],
        &[(0, 0, 0)],
    );
    let mut cache = workspace(&owner);
    let field = owner.kinematics(&gate, &mut cache).unwrap();
    assert!(field.sample_at(3, 2, 2).unwrap().is_none());
    let edge = field.sample_at(4, 2, 2).unwrap().unwrap();
    same_gradient(edge.kinematics, [[0.0; 3], [0.003, 0.0, 0.0], [0.0; 3]]);
    assert_eq!(edge.stencils[0], WetGradientStencil3::ForwardSecondOrder);

    let mut owner = FreeSurface3::new(
        8,
        8,
        8,
        0.8,
        [0.0; 3],
        0.0,
        ContactModel3::Neutral,
        |x, _, _| x < 4,
    )
    .unwrap();
    let mut cache = workspace(&owner);
    owner.kinematics(&gate, &mut cache).unwrap();
    owner
        .remap_planar_piston_one_cell(PlanarWall3::XMin, true, 1e-12, &gate)
        .unwrap();
    let field = owner.kinematics(&gate, &mut cache).unwrap();
    assert!(field.sample_at(0, 3, 3).unwrap().is_none());
    assert_eq!(
        field.sample_at(1, 3, 3).unwrap().unwrap().stencils[0],
        WetGradientStencil3::ForwardSecondOrder
    );
}

#[test]
fn sparse_budget_and_cancellation_refuse_without_touching_trajectory_and_retry_rebuilds() {
    let owner = liquid(
        [8; 3],
        |_| true,
        |p| [-0.002 * p[1] as f64, 0.002 * p[0] as f64, 0.0],
        &[],
    );
    let before = owner.grid().state_bits();
    let tiles = owner.grid().active_tiles();
    let required = FreeSurfaceKinematicsWorkspace3::required_bytes(tiles).unwrap();
    assert!(matches!(
        FreeSurfaceKinematicsWorkspace3::new(tiles, required - 1),
        Err(KinematicsError::Budget { .. })
    ));
    assert!(FreeSurfaceKinematicsWorkspace3::required_bytes(usize::MAX).is_err());
    let mut short = FreeSurfaceKinematicsWorkspace3::new(tiles - 1, required).unwrap();
    let gate = CancelGate::new_clock_free();
    assert!(matches!(
        owner.kinematics(&gate, &mut short),
        Err(FreeSurfaceKinematicsError3::Kinematics(
            KinematicsError::Budget { .. }
        ))
    ));
    let mut cache = workspace(&owner);
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(matches!(
        owner.kinematics(&cancelled, &mut cache),
        Err(FreeSurfaceKinematicsError3::Kinematics(
            KinematicsError::Cancelled
        ))
    ));
    let field = owner.kinematics(&gate, &mut cache).unwrap();
    assert!(matches!(
        field.visit_samples(|_, _| gate.request()),
        Err(FreeSurfaceKinematicsError3::Kinematics(
            KinematicsError::Cancelled
        ))
    ));
    let retry = CancelGate::new_clock_free();
    let field = owner.kinematics(&retry, &mut cache).unwrap();
    near(
        field
            .sample_at(3, 3, 3)
            .unwrap()
            .unwrap()
            .kinematics
            .vorticity()[2],
        0.004,
    );
    assert_eq!(owner.grid().state_bits(), before);
}

fn controls(model: SgsModel) -> SgsConfig {
    SgsConfig {
        model,
        filter_width_lu: 1.0,
        maximum_eddy_to_molecular: 100.0,
        usage: SgsUse::Exploration,
    }
}

#[test]
fn native_gradient_sgs_uses_rotation_and_fill_weighted_deviatoric_power() {
    let gate = CancelGate::new_clock_free();
    let g = [
        [-0.002, -0.004, 0.0],
        [0.004, -0.002, 0.0],
        [0.0, 0.0, 0.004],
    ];
    let owner = liquid(
        [4; 3],
        |_| true,
        |p| core::array::from_fn(|a| (0..3).map(|b| g[a][b] * p[b] as f64).sum()),
        &[],
    );
    let mut cache = workspace(&owner);
    let mut sgs = FreeSurfaceSgsWorkspace3::default();
    let mut density_volume = 0.0;
    let mut volume = 0.0;
    for z in 0..4 {
        for y in 0..4 {
            for x in 0..4 {
                let s = owner.sample_at(x, y, z).unwrap();
                density_volume += s.density * s.fill_fraction;
                volume += s.fill_fraction;
            }
        }
    }
    for model in [SgsModel::Wale { cw: 0.4 }, SgsModel::Vreman { cv: 0.1 }] {
        let eddy = eddy_viscosity(model, 1.0, g).unwrap();
        let symmetric = [[-0.002, 0.0, 0.0], [0.0, -0.002, 0.0], [0.0, 0.0, 0.004]];
        assert!((eddy - eddy_viscosity(model, 1.0, symmetric).unwrap()).abs() > 1e-5);
        let report = owner
            .prepare_gradient_sgs(&gate, controls(model), &mut cache, &mut sgs)
            .unwrap();
        near(report.maximum_eddy_viscosity_lu, eddy);
        near(sgs.relaxation_at([0, 0, 0]).unwrap(), 0.8 + 3.0 * eddy);
        near(report.liquid_volume_lu, volume);
        let rate_squared = 2.0 * (0.002_f64.powi(2) * 2.0 + 0.004_f64.powi(2));
        near(
            report.sgs_dissipation_lu,
            density_volume * rate_squared * eddy,
        );
        near(
            report.resolved_dissipation_lu,
            density_volume * rate_squared * 0.1,
        );
        assert_eq!(report.wet_cells, 64);
    }
}

#[test]
fn gradient_sgs_advances_original_bgk_and_trt_with_exact_zero_control() {
    let gate = CancelGate::new_clock_free();
    let initial: Vec<_> = (0..64)
        .map(|index| {
            let p = [index % 4, index / 4 % 4, index / 16];
            FreeSurfaceInitialCell3 {
                coordinates: p,
                density: 1.0,
                fill_fraction: 1.0,
                velocity: [
                    0.003 * p[1] as f64,
                    0.002 * p[2] as f64,
                    0.001 * p[0] as f64,
                ],
            }
        })
        .collect();
    for trt in [false, true] {
        let mut owner = with_cx(|cx| {
            FreeSurface3::from_cells(
                [4; 3],
                0.8,
                [0.0; 3],
                0.0,
                ContactModel3::Neutral,
                &[],
                &initial,
                cx,
            )
            .unwrap()
        });
        if trt {
            owner.set_trt_magic(0.25).unwrap();
        }
        owner.step_with_gate(&gate).unwrap();
        let model = owner.collision_model();
        let mut baseline = owner.clone();
        baseline.step_with_gate(&gate).unwrap();
        let mut cache = workspace(&owner);
        let mut sgs = FreeSurfaceSgsWorkspace3::default();
        for (disabled, enabled) in [
            (SgsModel::Wale { cw: 0.0 }, SgsModel::Wale { cw: 0.5 }),
            (SgsModel::Vreman { cv: 0.0 }, SgsModel::Vreman { cv: 0.1 }),
        ] {
            let mut zero = owner.clone();
            zero.step_with_gradient_sgs(&gate, controls(disabled), &mut cache, &mut sgs)
                .unwrap();
            assert_eq!(zero.grid().state_bits(), baseline.grid().state_bits());
            // The crate-private `mass_snapshot` (the mass of every fluid and interface cell, and the carry) is read through the public accessors: the mass of each of the 64 cells and the carry.
            let masses = |sim: &FreeSurface3| {
                let cells: Vec<f64> = (0..64)
                    .map(|index| sim.region_mass(|x, y, z| (x, y, z) == (index % 4, index / 4 % 4, index / 16)))
                    .collect();
                (cells, sim.carry())
            };
            assert_eq!(masses(&zero), masses(&baseline));
            let mut active = owner.clone();
            let report = active
                .step_with_gradient_sgs(&gate, controls(enabled), &mut cache, &mut sgs)
                .unwrap();
            assert!(report.sgs_dissipation_lu > 0.0);
            assert_eq!(active.collision_model(), model);
            let mut velocity_change: f64 = 0.0;
            for z in 0..4 {
                for y in 0..4 {
                    for x in 0..4 {
                        let a = active.sample_at(x, y, z).unwrap().velocity;
                        let b = baseline.sample_at(x, y, z).unwrap().velocity;
                        for axis in 0..3 {
                            velocity_change = velocity_change.max((a[axis] - b[axis]).abs());
                        }
                    }
                }
            }
            assert!(
                velocity_change > 1e-10,
                "closure must change physical evolution: {velocity_change}"
            );
        }
    }
}

#[test]
fn gradient_sgs_refuses_underresolved_geometry_and_excess_viscosity_atomically() {
    let gate = CancelGate::new_clock_free();
    let mut thin = liquid(
        [8; 3],
        |p| p[1] == 3 || p[1] == 4,
        |p| [0.002 * p[1] as f64, 0.0, 0.0],
        &[],
    );
    let before = thin.grid().state_bits();
    let mut cache = workspace(&thin);
    let mut sgs = FreeSurfaceSgsWorkspace3::default();
    assert!(matches!(
        thin.step_with_gradient_sgs(
            &gate,
            controls(SgsModel::Wale { cw: 0.4 }),
            &mut cache,
            &mut sgs
        ),
        Err(FreeSurfaceSgsError3::Cell { .. })
    ));
    assert_eq!(thin.grid().state_bits(), before);
    let mut owner = liquid(
        [4; 3],
        |_| true,
        |p| {
            [
                0.002 * p[1] as f64,
                0.003 * p[2] as f64,
                0.004 * p[0] as f64,
            ]
        },
        &[],
    );
    let before = owner.grid().state_bits();
    let mut cache = workspace(&owner);
    let excessive = SgsConfig {
        maximum_eddy_to_molecular: 0.0,
        ..controls(SgsModel::Vreman { cv: 0.1 })
    };
    assert!(
        owner
            .step_with_gradient_sgs(&gate, excessive, &mut cache, &mut sgs)
            .is_err()
    );
    assert_eq!(owner.grid().state_bits(), before);
    owner
        .prepare_gradient_sgs(
            &gate,
            controls(SgsModel::Vreman { cv: 0.1 }),
            &mut cache,
            &mut sgs,
        )
        .unwrap();
}
