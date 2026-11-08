use super::*;
use crate::d3q19::{BoundarySpec3, CollisionModel3, FaceBoundary3};
use crate::d3q27::{CollisionConfig27, CollisionModel27, D3q27, Rates27};
use crate::lattice3::{D3q19, EquilibriumLattice3, Grid3, GridBudget3};
use crate::sgs3::{SgsConfig, SgsLattice3, SgsModel, SgsUse, eddy_viscosity};
use crate::thermal::boundary::conduction::{WallCell, WallConductionControls, WallConductionLink};
use crate::thermal::boundary::{FiniteThermalWall, ThermalWall};
use crate::thermal::open::{
    ThermalOpenBalance3, ThermalOpenBoundaries3, ThermalOpenFace3, ThermalOpenField3,
    ThermalOpenLimits3,
};
use fs_exec::{CancelGate, Cx};

const DIMS: [usize; 3] = [6, 5, 4];
const CELLS: usize = 120;
const ORIGIN: [f64; 3] = [0.7, -0.3, 0.2];

fn with_gate<T>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> T) -> T {
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        f(&Cx::new(
            gate,
            arena,
            fs_exec::StreamKey {
                seed: 0x5307_0300,
                kernel_id: 227,
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
fn model27() -> CollisionConfig27 {
    CollisionConfig27 {
        model: CollisionModel27::Cumulant {
            rates: Rates27 {
                omega: [1.1, 1.2, 1.3, 1.1, 1.4, 1.2, 1.3, 1.2, 1.4, 1.1],
            },
        },
        maximum_mach: 0.2,
    }
}
fn models() -> [SgsModel; 3] {
    [
        SgsModel::Smagorinsky { cs: 0.3 },
        SgsModel::Wale { cw: 0.65 },
        SgsModel::Vreman { cv: 0.15 },
    ]
}
fn closure(model: SgsModel) -> SgsConfig {
    SgsConfig {
        model,
        filter_width_lu: 1.3,
        maximum_eddy_to_molecular: 20.0,
        usage: SgsUse::Exploration,
    }
}
fn thermal_config(model: SgsModel) -> OpenThermalSgsConfig3 {
    OpenThermalSgsConfig3 {
        closure: closure(model),
        turbulent_prandtl: 0.73,
        maximum_diffusivity_to_molecular: 20.0,
    }
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}"
    );
}
fn coordinates(cell: usize, dims: [usize; 3]) -> [usize; 3] {
    [
        cell % dims[0],
        (cell / dims[0]) % dims[1],
        cell / (dims[0] * dims[1]),
    ]
}
fn index(p: [usize; 3], dims: [usize; 3]) -> usize {
    p[0] + dims[0] * (p[1] + dims[1] * p[2])
}

fn make<L: EquilibriumLattice3<Q>, const Q: usize>(
    cx: &Cx<'_>,
    dims: [usize; 3],
    collision: L::Collision,
    faces: [FaceBoundary3; 6],
    velocity: impl Fn([usize; 3]) -> [f64; 3],
) -> OpenBoundaryGrid3<L, Q> {
    let grid = Grid3::try_new_equilibrium(
        dims,
        collision,
        [0.00013, -0.00007, 0.00009],
        core::array::from_fn(|axis| matches!(faces[2 * axis], FaceBoundary3::Periodic)),
        GridBudget3 {
            maximum_cells: dims.into_iter().product(),
            maximum_mach: 0.2,
        },
        |_| false,
        |p| (1.2, velocity(p)),
    )
    .unwrap();
    OpenBoundaryGrid3::try_new(grid, cx, BoundarySpec3::new(faces)).unwrap()
}

fn gradient_cases<L: EquilibriumLattice3<Q> + SgsLattice3<Q>, const Q: usize>(
    collision: L::Collision,
) {
    with_cx(|cx| {
        let dims = [5, 6, 7];
        let uniform = [0.027, -0.008, 0.011];
        for axis in 0..3 {
            let mut faces = [FaceBoundary3::Periodic; 6];
            faces[2 * axis] = FaceBoundary3::Velocity { velocity: uniform };
            faces[2 * axis + 1] = FaceBoundary3::Pressure { density: 1.2 };
            let flow = make::<L, Q>(cx, dims, collision, faces, |_| uniform);
            let bytes = flow.required_sgs_workspace_bytes().unwrap();
            assert!(flow.sgs_workspace(bytes - 1).is_err());
            let mut workspace = flow.sgs_workspace(bytes).unwrap();
            assert!(workspace.workspace_bytes() <= bytes);
            for model in models() {
                let report = flow
                    .prepare_sgs(cx, closure(model), None, &mut workspace)
                    .unwrap();
                assert!(
                    report.maximum_eddy_viscosity_lu < 1e-14,
                    "uniform physical flow must not see reservoir wall shear: {report:?}"
                );
            }

            // Accepted initial fields need not equal a later prescribed drive.
            // Nonzero transverse slopes exercise WALE/Vreman's full gradient;
            // probes avoid periodic seams, while the open normal is quadratic.
            let a = axis;
            let b = (axis + 1) % 3;
            let c = (axis + 2) % 3;
            let profile = |p: [usize; 3]| {
                let [x, y, z] = [p[a] as f64, p[b] as f64, p[c] as f64];
                let mut u = [0.0; 3];
                u[a] = 0.025 + 0.001 * x + 0.0001 * x * x + 0.0002 * y;
                u[b] = 0.003 + 0.0013 * x + 0.00008 * x * x - 0.0003 * y + 0.0002 * z;
                u[c] = -0.004 + 0.0007 * x - 0.00005 * x * x + 0.0004 * y - 0.0001 * z;
                u
            };
            for pressure_outlet in [false, true] {
                faces[2 * a] = FaceBoundary3::Velocity {
                    velocity: profile([0; 3]),
                };
                let mut high = [0; 3];
                high[a] = dims[a] - 1;
                faces[2 * a + 1] = if pressure_outlet {
                    FaceBoundary3::Pressure { density: 1.2 }
                } else {
                    FaceBoundary3::Velocity {
                        velocity: profile(high),
                    }
                };
                let flow = make::<L, Q>(cx, dims, collision, faces, profile);
                for model in models() {
                    flow.prepare_sgs(cx, closure(model), None, &mut workspace)
                        .unwrap();
                    for x in [0, dims[a] / 2, dims[a] - 1] {
                        let mut p = [2; 3];
                        p[a] = x;
                        let mut g = [[0.0; 3]; 3];
                        g[a][a] = 0.001 + 0.0002 * x as f64;
                        g[b][a] = 0.0013 + 0.00016 * x as f64;
                        g[c][a] = 0.0007 - 0.0001 * x as f64;
                        g[a][b] = 0.0002;
                        g[b][b] = -0.0003;
                        g[b][c] = 0.0002;
                        g[c][b] = 0.0004;
                        g[c][c] = -0.0001;
                        if pressure_outlet && x == dims[a] - 1 {
                            for row in &mut g {
                                row[a] = 0.0;
                            }
                        }
                        let eddy = eddy_viscosity(model, 1.3, g).unwrap();
                        let cell = index(p, dims);
                        near(workspace.eddy_viscosities()[cell], eddy, 3e-15);
                        near(
                            L::molecular_viscosity(workspace.collisions()[cell]).unwrap(),
                            L::molecular_viscosity(collision).unwrap() + eddy,
                            2e-15,
                        );
                    }
                }
            }
        }
    });
}

#[test]
fn reservoirs_use_physical_full_gradients_in_d3q19_and_d3q27() {
    gradient_cases::<D3q19, 19>(CollisionModel3::Trt {
        tau: 0.9,
        magic: 0.19,
    });
    gradient_cases::<D3q27, 27>(model27());
}

fn duct<L: EquilibriumLattice3<Q>, const Q: usize>(
    cx: &Cx<'_>,
    collision: L::Collision,
) -> OpenBoundaryGrid3<L, Q> {
    let mut faces = [FaceBoundary3::Periodic; 6];
    faces[0] = FaceBoundary3::Velocity {
        velocity: [0.018, 0.0, 0.0],
    };
    faces[1] = FaceBoundary3::Pressure { density: 1.2 };
    faces[2] = FaceBoundary3::Wall { velocity: [0.0; 3] };
    faces[3] = FaceBoundary3::Wall { velocity: [0.0; 3] };
    make::<L, Q>(cx, DIMS, collision, faces, |p| {
        let [x, y, z] =
            core::array::from_fn(|a| core::f64::consts::TAU * (p[a] as f64 + 0.5) / DIMS[a] as f64);
        [
            0.018 + 0.002 * z.cos() + 0.001 * x.sin(),
            0.001 * x.sin(),
            0.001 * y.cos(),
        ]
    })
}
fn thermal<L: EquilibriumLattice3<Q>, const Q: usize>(
    flow: &OpenBoundaryGrid3<L, Q>,
    cx: &Cx<'_>,
) -> ThermalOpenField3 {
    let mut faces = [ThermalOpenFace3::Periodic; 6];
    faces[0] = ThermalOpenFace3::Reservoir { temperature: 0.45 };
    faces[1] = ThermalOpenFace3::Reservoir { temperature: 0.5 };
    for side in 0..2 {
        faces[2 + side] = ThermalOpenFace3::Wall(ThermalWall::FiniteSolid(FiniteThermalWall {
            temperature: if side == 0 { 0.85 } else { 0.65 },
            heat_capacity: 0.2 + 0.1 * side as f64,
            resistance: 0.4 + 0.3 * side as f64,
            heat_per_step: 0.0001,
        }));
    }
    let initial: Vec<_> = (0..CELLS)
        .map(|cell| {
            let [x, y, z] = coordinates(cell, DIMS);
            0.45 + 0.07 * (core::f64::consts::TAU * (z as f64 + 0.5) / DIMS[2] as f64).cos()
                + 0.03 * (core::f64::consts::TAU * x as f64 / DIMS[0] as f64).sin()
                + 0.03 * y as f64 / DIMS[1] as f64
        })
        .collect();
    let mut thermal = flow
        .prepare_thermal_field(
            cx,
            &initial,
            0.85,
            ThermalOpenBoundaries3 {
                faces,
                internal_wall: ThermalWall::Adiabatic,
            },
            ThermalOpenLimits3 {
                maximum_cells: CELLS,
                maximum_storage_bytes: 1024 * 1024,
                temperature_range: [-1.0, 1.2],
                balance_tolerance: 1e-10,
            },
        )
        .unwrap();
    thermal
        .set_wall_conduction(
            1,
            &[
                WallConductionLink {
                    a: WallCell { side: 0, link: 1 },
                    b: WallCell { side: 0, link: 2 },
                    conductance: 0.3,
                },
                WallConductionLink {
                    a: WallCell { side: 0, link: 2 },
                    b: WallCell { side: 1, link: 2 },
                    conductance: 0.15,
                },
            ],
            WallConductionControls {
                max_nodes: 3,
                max_links: 2,
                max_factor_bytes: 1024 * 1024,
                balance_tolerance: 1e-10,
            },
        )
        .unwrap();
    thermal
}
fn molecular<L: SgsLattice3<Q>, const Q: usize>(
    collision: L::Collision,
) -> (Vec<L::Collision>, Vec<f64>) {
    (
        (0..CELLS)
            .map(|cell| {
                L::with_shear_relaxation(collision, 0.86 + 0.005 * (cell % 7) as f64).unwrap()
            })
            .collect(),
        (0..CELLS)
            .map(|cell| 0.78 + 0.012 * (cell % 5) as f64)
            .collect(),
    )
}
fn source() -> Vec<f64> {
    (0..CELLS)
        .map(|cell| 0.0001 * (1 + cell % 4) as f64)
        .collect()
}
#[derive(Debug, PartialEq)]
struct Snapshot {
    flow: Vec<u64>,
    heat: Vec<u64>,
    walls: Vec<Vec<ThermalWall>>,
    flow_step: usize,
    heat_step: usize,
    balance: Option<ThermalOpenBalance3>,
}
fn snapshot<L: EquilibriumLattice3<Q>, const Q: usize>(
    flow: &OpenBoundaryGrid3<L, Q>,
    thermal: &ThermalOpenField3,
) -> Snapshot {
    Snapshot {
        flow: flow
            .populations()
            .iter()
            .flatten()
            .map(|f| f.to_bits())
            .collect(),
        heat: thermal
            .populations()
            .iter()
            .flatten()
            .map(|g| g.to_bits())
            .collect(),
        walls: (0..6)
            .map(|face| thermal.wall_states(face).unwrap().to_vec())
            .collect(),
        flow_step: flow.completed_steps(),
        heat_step: thermal.completed_steps(),
        balance: thermal.last_balance(),
    }
}
fn energy(thermal: &ThermalOpenField3) -> f64 {
    let fluid: f64 = thermal.populations().iter().flatten().sum();
    fluid
        + (0..6)
            .flat_map(|face| thermal.wall_states(face).unwrap())
            .map(|wall| match wall {
                ThermalWall::FiniteSolid(solid) => solid.heat_capacity * solid.temperature,
                _ => 0.0,
            })
            .sum::<f64>()
}
fn wall_source(thermal: &ThermalOpenField3) -> f64 {
    (0..6)
        .flat_map(|face| thermal.wall_states(face).unwrap())
        .map(|wall| match wall {
            ThermalWall::FiniteSolid(solid) => solid.heat_per_step,
            _ => 0.0,
        })
        .sum()
}

fn zero_case<L: EquilibriumLattice3<Q> + SgsLattice3<Q>, const Q: usize>(collision: L::Collision)
where
    L::Collision: core::fmt::Debug + PartialEq,
{
    with_cx(|cx| {
        let (laws, tau) = molecular::<L, Q>(collision);
        let heat = source();
        for model in [
            SgsModel::Smagorinsky { cs: 0.0 },
            SgsModel::Wale { cw: 0.0 },
            SgsModel::Vreman { cv: 0.0 },
        ] {
            let mut flow = duct::<L, Q>(cx, collision);
            let mut reference = duct::<L, Q>(cx, collision);
            let mut field = thermal(&flow, cx);
            let mut expected = thermal(&reference, cx);
            let mut workspace = flow.thermal_sgs_workspace(1024 * 1024).unwrap();
            let pointers = (
                workspace.momentum_collisions().as_ptr(),
                workspace.thermal_relaxation().as_ptr(),
                workspace.eddy_viscosities().as_ptr(),
            );
            for _ in 0..4 {
                let actual = flow
                    .step_serial_with_thermal_sgs(
                        &mut field,
                        cx,
                        thermal_config(model),
                        Some(&laws),
                        &heat,
                        &tau,
                        &mut workspace,
                        ORIGIN,
                    )
                    .unwrap();
                let expected_step = reference
                    .step_serial_with_thermal(
                        &mut expected,
                        cx,
                        Some(&laws),
                        &heat,
                        Some(&tau),
                        ORIGIN,
                    )
                    .unwrap();
                assert_eq!(snapshot(&flow, &field), snapshot(&reference, &expected));
                assert_eq!(actual.transport.momentum, expected_step.momentum);
                assert_eq!(actual.transport.heat, expected_step.heat);
                assert_eq!(actual.sgs.sgs_dissipation_lu, 0.0);
                assert_eq!(actual.maximum_eddy_diffusivity_lu, 0.0);
                assert_eq!(workspace.momentum_collisions(), laws);
                assert_eq!(
                    workspace
                        .thermal_relaxation()
                        .iter()
                        .map(|v| v.to_bits())
                        .collect::<Vec<_>>(),
                    tau.iter().map(|v| v.to_bits()).collect::<Vec<_>>()
                );
                assert_eq!(
                    pointers,
                    (
                        workspace.momentum_collisions().as_ptr(),
                        workspace.thermal_relaxation().as_ptr(),
                        workspace.eddy_viscosities().as_ptr()
                    )
                );
            }
        }
    });
}

#[test]
fn zero_models_preserve_original_coupled_sources_and_conducting_walls_bitwise() {
    zero_case::<D3q19, 19>(CollisionModel3::Trt {
        tau: 0.9,
        magic: 0.19,
    });
    zero_case::<D3q27, 27>(model27());
}

#[test]
fn active_stress_and_heat_match_native_composition_and_close_fluid_solid_energy() {
    with_cx(|cx| {
        let (laws, tau) = molecular::<D3q27, 27>(model27());
        let heat = source();
        for model in models() {
            let config = thermal_config(model);
            let mut flow = duct::<D3q27, 27>(cx, model27());
            let mut reference = duct::<D3q27, 27>(cx, model27());
            let mut molecular_flow = duct::<D3q27, 27>(cx, model27());
            let mut field = thermal(&flow, cx);
            let mut expected = thermal(&reference, cx);
            let mut molecular_field = thermal(&molecular_flow, cx);
            // Equilibrium populations have no viscous defect on their first
            // collision. Warm the actual native carriers before comparison.
            for (flow, field) in [
                (&mut flow, &mut field),
                (&mut reference, &mut expected),
                (&mut molecular_flow, &mut molecular_field),
            ] {
                flow.step_serial_with_thermal(field, cx, Some(&laws), &heat, Some(&tau), ORIGIN)
                    .unwrap();
            }
            let mut workspace = flow.thermal_sgs_workspace(1024 * 1024).unwrap();
            let mut momentum = reference.sgs_workspace(1024 * 1024).unwrap();
            for _ in 0..3 {
                let before = energy(&field);
                let external = heat.iter().sum::<f64>() + wall_source(&field);
                let fit = reference
                    .prepare_sgs(cx, config.closure, Some(&laws), &mut momentum)
                    .unwrap();
                let total_tau: Vec<_> = tau
                    .iter()
                    .zip(momentum.eddy_viscosities())
                    .map(|(&molecular, &eddy)| {
                        if eddy == 0.0 {
                            molecular
                        } else {
                            molecular + 4.0 * (eddy / config.turbulent_prandtl)
                        }
                    })
                    .collect();
                let actual = flow
                    .step_serial_with_thermal_sgs(
                        &mut field,
                        cx,
                        config,
                        Some(&laws),
                        &heat,
                        &tau,
                        &mut workspace,
                        ORIGIN,
                    )
                    .unwrap();
                let expected_step = reference
                    .step_serial_with_thermal(
                        &mut expected,
                        cx,
                        Some(momentum.collisions()),
                        &heat,
                        Some(&total_tau),
                        ORIGIN,
                    )
                    .unwrap();
                molecular_flow
                    .step_serial_with_thermal(
                        &mut molecular_field,
                        cx,
                        Some(&laws),
                        &heat,
                        Some(&tau),
                        ORIGIN,
                    )
                    .unwrap();
                assert_eq!(actual.sgs, fit);
                assert_eq!(snapshot(&flow, &field), snapshot(&reference, &expected));
                assert_eq!(actual.transport.momentum, expected_step.momentum);
                assert_eq!(workspace.thermal_relaxation(), total_tau);
                assert_eq!(workspace.eddy_viscosities(), momentum.eddy_viscosities());
                assert!(actual.sgs.sgs_dissipation_lu > 0.0);
                assert!(actual.maximum_eddy_diffusivity_lu > 0.0);
                near(
                    actual.maximum_eddy_diffusivity_lu,
                    fit.maximum_eddy_viscosity_lu / config.turbulent_prandtl,
                    1e-16,
                );
                let ratio = tau
                    .iter()
                    .zip(momentum.eddy_viscosities())
                    .map(|(&tau, &eddy)| (eddy / config.turbulent_prandtl) / (0.25 * (tau - 0.5)))
                    .fold(0.0_f64, f64::max);
                near(actual.maximum_diffusivity_to_molecular, ratio, 2e-15);
                let reservoir: f64 = actual
                    .transport
                    .heat
                    .reservoirs
                    .iter()
                    .map(|face| face.heat_in)
                    .sum();
                near(energy(&field) - before, external + reservoir, 2e-11);
                assert!(actual.transport.heat.defect.abs() < 1e-11);
                assert!(actual.transport.heat.combined_defect.abs() < 1e-11);
            }
            assert_ne!(
                flow.populations(),
                molecular_flow.populations(),
                "{model:?}"
            );
            assert!(
                field
                    .populations()
                    .iter()
                    .flatten()
                    .zip(molecular_field.populations().iter().flatten())
                    .any(|(a, b)| (a - b).abs() > 1e-14),
                "SGS must change thermal transport for {model:?}"
            );
            assert_ne!(field.wall_states(2), molecular_field.wall_states(2));
        }
    });
}

#[test]
fn warmed_flow_heat_wall_and_receipt_state_survive_refusals_then_retry() {
    with_cx(|cx| {
        let config = thermal_config(SgsModel::Smagorinsky { cs: 0.3 });
        let (laws, tau) = molecular::<D3q19, 19>(CollisionModel3::Trt {
            tau: 0.9,
            magic: 0.19,
        });
        let heat = source();
        let mut flow = duct::<D3q19, 19>(
            cx,
            CollisionModel3::Trt {
                tau: 0.9,
                magic: 0.19,
            },
        );
        let mut reference = duct::<D3q19, 19>(
            cx,
            CollisionModel3::Trt {
                tau: 0.9,
                magic: 0.19,
            },
        );
        let mut field = thermal(&flow, cx);
        let mut expected = thermal(&reference, cx);
        let mut workspace = flow.thermal_sgs_workspace(1024 * 1024).unwrap();
        let mut expected_workspace = reference.thermal_sgs_workspace(1024 * 1024).unwrap();
        let tiny = thermal_config(SgsModel::Smagorinsky { cs: 1e-12 });
        for current in [config, config, tiny] {
            let accepted = flow
                .step_serial_with_thermal_sgs(
                    &mut field,
                    cx,
                    current,
                    Some(&laws),
                    &heat,
                    &tau,
                    &mut workspace,
                    ORIGIN,
                )
                .unwrap();
            reference
                .step_serial_with_thermal_sgs(
                    &mut expected,
                    cx,
                    current,
                    Some(&laws),
                    &heat,
                    &tau,
                    &mut expected_workspace,
                    ORIGIN,
                )
                .unwrap();
            if current == tiny {
                // Positive diffusion below the molecular tau's ulp is normal
                // rounded addition, not an invalid or clipped constitutive law.
                assert!(accepted.maximum_eddy_diffusivity_lu > 0.0);
                assert_eq!(workspace.thermal_relaxation(), tau);
            }
        }
        let before = snapshot(&flow, &field);
        assert!(flow.thermal_sgs_workspace(1).is_err());
        let capped = OpenThermalSgsConfig3 {
            maximum_diffusivity_to_molecular: 1e-12,
            ..config
        };
        assert!(
            flow.step_serial_with_thermal_sgs(
                &mut field,
                cx,
                capped,
                Some(&laws),
                &heat,
                &tau,
                &mut workspace,
                ORIGIN
            )
            .is_err()
        );
        assert_eq!(snapshot(&flow, &field), before);
        let mut bad_tau = tau.clone();
        bad_tau[CELLS - 1] = 0.5;
        assert!(
            flow.step_serial_with_thermal_sgs(
                &mut field,
                cx,
                config,
                Some(&laws),
                &heat,
                &bad_tau,
                &mut workspace,
                ORIGIN
            )
            .is_err()
        );
        assert_eq!(snapshot(&flow, &field), before);
        // The first half-kick stays inside the 1.2 temperature ceiling;
        // the final half-kick refuses after native momentum/wall staging.
        assert!((0..CELLS).all(|cell| field.temperature(cell).unwrap() + 0.45 < 1.2));
        let refused = flow.step_serial_with_thermal_sgs(
            &mut field,
            cx,
            config,
            Some(&laws),
            &vec![0.9; CELLS],
            &tau,
            &mut workspace,
            ORIGIN,
        );
        assert!(refused.is_err());
        assert!(
            format!("{refused:?}").contains("Temperature"),
            "{refused:?}"
        );
        assert_eq!(snapshot(&flow, &field), before);
        // A real executor cancellation is separate from numerical refusal;
        // avoid nondeterministic thread timing to manufacture a late request.
        let cancelled = CancelGate::new_clock_free();
        cancelled.request();
        with_gate(&cancelled, |cx| {
            assert!(
                flow.step_serial_with_thermal_sgs(
                    &mut field,
                    cx,
                    config,
                    Some(&laws),
                    &heat,
                    &tau,
                    &mut workspace,
                    ORIGIN
                )
                .is_err()
            );
        });
        assert_eq!(snapshot(&flow, &field), before);
        let actual = flow
            .step_serial_with_thermal_sgs(
                &mut field,
                cx,
                config,
                Some(&laws),
                &heat,
                &tau,
                &mut workspace,
                ORIGIN,
            )
            .unwrap();
        let clean = reference
            .step_serial_with_thermal_sgs(
                &mut expected,
                cx,
                config,
                Some(&laws),
                &heat,
                &tau,
                &mut expected_workspace,
                ORIGIN,
            )
            .unwrap();
        assert_eq!(snapshot(&flow, &field), snapshot(&reference, &expected));
        assert_eq!(actual.transport.momentum, clean.transport.momentum);
        assert_eq!(actual.transport.heat, clean.transport.heat);
        assert_eq!(actual.sgs, clean.sgs);
    });
}
