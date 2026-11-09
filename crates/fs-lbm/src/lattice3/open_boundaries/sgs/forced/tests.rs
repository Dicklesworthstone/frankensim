use super::*;
use crate::d3q19::{BoundarySpec3, CollisionModel3, FaceBoundary3};
use crate::d3q27::{CollisionConfig27, CollisionModel27, D3q27, Rates27};
use crate::lattice3::open_boundaries::OpenBoundaryGrid3;
use crate::lattice3::{D3q19, Grid3, GridBudget3};
use crate::sgs3::{SgsConfig, SgsModel, SgsUse};
use crate::thermal::boundary::{FiniteThermalWall, ThermalWall};
use crate::thermal::open::{
    ThermalOpenBalance3, ThermalOpenBoundaries3, ThermalOpenFace3, ThermalOpenLimits3,
};
use fs_exec::CancelGate;

const DIMS: [usize; 3] = [6, 4, 3];
const CELLS: usize = 72;
const BYTES: usize = 1 << 20;

fn with_cx<T>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> T) -> T {
    let arena = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    arena.scope(|arena| {
        f(&Cx::new(
            gate,
            arena,
            fs_exec::StreamKey {
                seed: 719,
                kernel_id: 181,
                tile: 0,
                iteration: 0,
            },
            fs_exec::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        ))
    })
}

fn coordinates(cell: usize) -> [usize; 3] {
    [
        cell % DIMS[0],
        cell / DIMS[0] % DIMS[1],
        cell / (DIMS[0] * DIMS[1]),
    ]
}
fn forces() -> Vec<[f64; 3]> {
    (0..CELLS)
        .map(|cell| {
            let [x, y, z] = coordinates(cell);
            [
                0.0004 * (x + 1) as f64,
                -0.0002 * (y + 1) as f64,
                0.0003 * (z + 1) as f64,
            ]
        })
        .collect()
}
fn config(model: SgsModel) -> OpenThermalSgsConfig3 {
    OpenThermalSgsConfig3 {
        closure: SgsConfig {
            model,
            filter_width_lu: 1.2,
            maximum_eddy_to_molecular: 50.0,
            usage: SgsUse::Exploration,
        },
        turbulent_prandtl: 0.7,
        maximum_diffusivity_to_molecular: 50.0,
    }
}
fn models() -> [SgsModel; 3] {
    [
        SgsModel::Smagorinsky { cs: 0.2 },
        SgsModel::Wale { cw: 0.5 },
        SgsModel::Vreman { cv: 0.1 },
    ]
}
fn law27() -> CollisionConfig27 {
    CollisionConfig27 {
        model: CollisionModel27::Cumulant {
            rates: Rates27 {
                omega: [1.1, 1.2, 1.3, 1.1, 1.4, 1.2, 1.3, 1.2, 1.4, 1.1],
            },
        },
        maximum_mach: 0.2,
    }
}
fn base<L: EquilibriumLattice3<Q>, const Q: usize>(
    cx: &Cx<'_>,
    collision: L::Collision,
    uniform: bool,
) -> OpenBoundaryGrid3<L, Q> {
    let mut faces = [FaceBoundary3::Periodic; 6];
    faces[0] = FaceBoundary3::Velocity {
        velocity: [0.02, 0.0, 0.0],
    };
    faces[1] = FaceBoundary3::Pressure { density: 1.2 };
    if !uniform {
        faces[2] = FaceBoundary3::stationary_wall();
        faces[3] = FaceBoundary3::Wall {
            velocity: [0.0, 0.0, 0.01],
        };
    }
    let grid = Grid3::try_new_equilibrium(
        DIMS,
        collision,
        [0.00013, -0.00007, 0.00009],
        [false, uniform, true],
        GridBudget3 {
            maximum_cells: CELLS,
            maximum_mach: 0.2,
        },
        |_| false,
        |p| {
            (
                1.2,
                if uniform {
                    [0.02, 0.0, 0.0]
                } else {
                    [
                        0.02 + 0.001 * p[1] as f64,
                        0.0004 * p[0] as f64,
                        0.0006 * p[0] as f64 + 0.0002 * p[1] as f64,
                    ]
                },
            )
        },
    )
    .unwrap();
    OpenBoundaryGrid3::try_new(grid, cx, BoundarySpec3::new(faces)).unwrap()
}
fn flow<L: EquilibriumLattice3<Q>, const Q: usize>(
    cx: &Cx<'_>,
    collision: L::Collision,
    uniform: bool,
) -> ForcedOpenBoundaryGrid3<L, Q> {
    ForcedOpenBoundaryGrid3::try_new_preserving_velocity(
        base::<L, Q>(cx, collision, uniform),
        cx,
        forces(),
    )
    .unwrap()
}
fn heat<L: EquilibriumLattice3<Q>, const Q: usize>(
    flow: &ForcedOpenBoundaryGrid3<L, Q>,
    cx: &Cx<'_>,
    uniform: bool,
) -> ThermalOpenField3 {
    let mut faces = [ThermalOpenFace3::Periodic; 6];
    faces[0] = ThermalOpenFace3::Reservoir { temperature: 0.6 };
    faces[1] = ThermalOpenFace3::Reservoir { temperature: 0.7 };
    if !uniform {
        faces[2] = ThermalOpenFace3::Wall(ThermalWall::FiniteSolid(FiniteThermalWall {
            temperature: 0.8,
            heat_capacity: 0.2,
            resistance: 0.4,
            heat_per_step: 0.0001,
        }));
        faces[3] = ThermalOpenFace3::Wall(ThermalWall::Adiabatic);
    }
    let initial: Vec<_> = (0..CELLS)
        .map(|cell| 0.5 + 0.01 * coordinates(cell)[1] as f64)
        .collect();
    flow.prepare_thermal_field(
        cx,
        &initial,
        0.8,
        ThermalOpenBoundaries3 {
            faces,
            internal_wall: ThermalWall::Adiabatic,
        },
        ThermalOpenLimits3 {
            maximum_cells: CELLS,
            maximum_storage_bytes: BYTES,
            temperature_range: [0.0, 2.0],
            balance_tolerance: 1e-10,
        },
    )
    .unwrap()
}
fn near(a: f64, b: f64) {
    assert!((a - b).abs() < 3e-14, "{a:.17e} != {b:.17e}");
}
#[derive(Debug, PartialEq)]
struct Snapshot {
    flow: Vec<u64>,
    heat: Vec<u64>,
    forces: Vec<[f64; 3]>,
    walls: Vec<Vec<ThermalWall>>,
    times: [usize; 2],
    balance: Option<ThermalOpenBalance3>,
    wall: Option<crate::lattice3::moving_walls::MovingWallMomentum3>,
}
fn snapshot<L: EquilibriumLattice3<Q>, const Q: usize>(
    flow: &ForcedOpenBoundaryGrid3<L, Q>,
    heat: &ThermalOpenField3,
) -> Snapshot {
    Snapshot {
        flow: flow
            .populations()
            .iter()
            .flatten()
            .map(|v| v.to_bits())
            .collect(),
        heat: heat
            .populations()
            .iter()
            .flatten()
            .map(|v| v.to_bits())
            .collect(),
        forces: flow.forces().to_vec(),
        walls: (0..6)
            .map(|face| heat.wall_states(face).unwrap().to_vec())
            .collect(),
        times: [flow.completed_steps(), heat.completed_steps()],
        balance: heat.last_balance(),
        wall: flow.last_wall_exchange(),
    }
}

fn constant_velocity<L: EquilibriumLattice3<Q> + SgsLattice3<Q>, const Q: usize>(
    law: L::Collision,
) {
    with_cx(&CancelGate::new_clock_free(), |cx| {
        for model in models() {
            let mut actual = flow::<L, Q>(cx, law, true);
            let mut thermal = heat(&actual, cx, true);
            let bytes = actual.required_thermal_sgs_workspace_bytes().unwrap();
            assert!(actual.thermal_sgs_workspace(bytes - 1).is_err());
            let mut workspace = actual.thermal_sgs_workspace(bytes).unwrap();
            let mut wrong = actual.inner.sgs_workspace(BYTES).unwrap();
            let raw = actual
                .inner
                .prepare_sgs(cx, config(model).closure, None, &mut wrong)
                .unwrap();
            assert!(
                raw.maximum_eddy_viscosity_lu > 1e-9,
                "negative control must expose false strain from omitted spatial F_n"
            );
            let report = actual
                .step_serial_with_thermal_sgs(
                    &mut thermal,
                    cx,
                    config(model),
                    None,
                    &vec![0.0001; CELLS],
                    &vec![0.8; CELLS],
                    &mut workspace,
                    [0.0; 3],
                )
                .unwrap();
            assert!(report.sgs.maximum_eddy_viscosity_lu < 1e-14);
            assert!(report.maximum_eddy_diffusivity_lu < 1e-14);
            near(report.transport.heat.source, 0.0001 * CELLS as f64);
            assert!(report.transport.heat.combined_defect.abs() < 1e-10);
        }
    });
}

#[test]
fn constant_physical_velocity_with_spatial_force_does_not_generate_turbulent_diffusion() {
    constant_velocity::<D3q19, 19>(CollisionModel3::Bgk { tau: 0.9 });
    constant_velocity::<D3q27, 27>(law27());
}

fn coupled<L: EquilibriumLattice3<Q> + SgsLattice3<Q>, const Q: usize>(law: L::Collision) {
    with_cx(&CancelGate::new_clock_free(), |cx| {
        for model in models()
            .into_iter()
            .chain([SgsModel::Smagorinsky { cs: 0.0 }])
        {
            for advance in [false, true] {
                let mut actual = flow::<L, Q>(cx, law, false);
                let mut native = flow::<L, Q>(cx, law, false);
                let mut thermal = heat(&actual, cx, false);
                let mut native_heat = heat(&native, cx, false);
                let molecular: Vec<_> = (0..CELLS)
                    .map(|cell| {
                        L::with_shear_relaxation(law, 0.83 + 0.01 * (cell % 5) as f64).unwrap()
                    })
                    .collect();
                let tau: Vec<_> = (0..CELLS)
                    .map(|cell| 0.77 + 0.01 * (cell % 3) as f64)
                    .collect();
                let source = vec![0.0001; CELLS];
                let mut workspace = actual.thermal_sgs_workspace(BYTES).unwrap();
                let pointers = (
                    workspace.momentum_collisions().as_ptr(),
                    workspace.thermal_relaxation().as_ptr(),
                );
                let prepared = actual
                    .inner
                    .prepare_open_thermal_sgs(
                        config(model),
                        Some(&molecular),
                        &tau,
                        &mut workspace,
                        Some(&actual.forces),
                        || cx.checkpoint().map_err(|_| GridError3::Cancelled),
                    )
                    .unwrap();
                let reference = base::<L, Q>(cx, law, false);
                let mut reference_work = reference.sgs_workspace(BYTES).unwrap();
                reference
                    .prepare_sgs(
                        cx,
                        config(model).closure,
                        Some(&molecular),
                        &mut reference_work,
                    )
                    .unwrap();
                for cell in 0..CELLS {
                    let eddy = workspace.eddy_viscosities()[cell];
                    near(eddy, reference_work.eddy_viscosities()[cell]);
                    near(
                        L::molecular_viscosity(workspace.momentum_collisions()[cell]).unwrap(),
                        L::molecular_viscosity(molecular[cell]).unwrap() + eddy,
                    );
                    near(
                        workspace.thermal_relaxation()[cell],
                        tau[cell] + 4.0 * eddy / config(model).turbulent_prandtl,
                    );
                    if matches!(model, SgsModel::Smagorinsky { cs: 0.0 }) {
                        assert_eq!(
                            workspace.thermal_relaxation()[cell].to_bits(),
                            tau[cell].to_bits()
                        );
                    }
                }
                let mut endpoint: Vec<_> = forces()
                    .into_iter()
                    .map(|f| [-0.5 * f[0], 1.5 * f[1], 0.2 * f[2]])
                    .collect();
                let next_force = endpoint.clone();
                let mut native_endpoint = endpoint.clone();
                let expected = if advance {
                    native
                        .advance_serial_with_thermal(
                            &mut native_heat,
                            cx,
                            Some(workspace.momentum_collisions()),
                            &mut native_endpoint,
                            &source,
                            Some(workspace.thermal_relaxation()),
                            [0.0; 3],
                        )
                        .unwrap()
                } else {
                    native
                        .step_serial_with_thermal(
                            &mut native_heat,
                            cx,
                            Some(workspace.momentum_collisions()),
                            &source,
                            Some(workspace.thermal_relaxation()),
                            [0.0; 3],
                        )
                        .unwrap()
                };
                let report = if advance {
                    actual
                        .advance_serial_with_thermal_sgs(
                            &mut thermal,
                            cx,
                            config(model),
                            Some(&molecular),
                            &mut endpoint,
                            &source,
                            &tau,
                            &mut workspace,
                            [0.0; 3],
                        )
                        .unwrap()
                } else {
                    actual
                        .step_serial_with_thermal_sgs(
                            &mut thermal,
                            cx,
                            config(model),
                            Some(&molecular),
                            &source,
                            &tau,
                            &mut workspace,
                            [0.0; 3],
                        )
                        .unwrap()
                };
                assert_eq!(report.sgs, prepared.0);
                assert_eq!(report.transport.heat, expected.heat);
                assert_eq!(snapshot(&actual, &thermal), snapshot(&native, &native_heat));
                assert_eq!(
                    pointers,
                    (
                        workspace.momentum_collisions().as_ptr(),
                        workspace.thermal_relaxation().as_ptr()
                    )
                );
                assert!(report.transport.heat.combined_defect.abs() < 1e-10);
                assert!(actual.last_wall_exchange().unwrap().wall_work.is_finite());
                if advance {
                    assert_eq!(actual.forces(), next_force);
                    assert_eq!(endpoint, forces());
                    assert_eq!(endpoint, native_endpoint);
                }
            }
        }
    });
}

#[test]
fn native_force_heat_and_wall_transaction_consumes_one_incoming_eddy_field() {
    coupled::<D3q19, 19>(CollisionModel3::Trt {
        tau: 0.9,
        magic: 0.19,
    });
    coupled::<D3q27, 27>(law27());
}

#[test]
fn endpoint_heat_closure_and_cancellation_refusals_preserve_all_accepted_epochs() {
    for failure in 0..4 {
        let gate = CancelGate::new_clock_free();
        with_cx(&gate, |cx| {
            let mut actual = flow::<D3q27, 27>(cx, law27(), false);
            let mut thermal = heat(&actual, cx, false);
            let before = snapshot(&actual, &thermal);
            let mut workspace = actual.thermal_sgs_workspace(BYTES).unwrap();
            let mut endpoint = vec![[0.0001, 0.0, -0.0002]; CELLS];
            let mut source = vec![0.0001; CELLS];
            let mut controls = config(models()[0]);
            match failure {
                0 => endpoint[CELLS - 1][0] = f64::NAN,
                1 => source.fill(1.6),
                2 => controls.maximum_diffusivity_to_molecular = 0.0,
                _ => gate.request(),
            }
            let endpoint_before: Vec<_> = endpoint.iter().flatten().map(|v| v.to_bits()).collect();
            assert!(
                actual
                    .advance_serial_with_thermal_sgs(
                        &mut thermal,
                        cx,
                        controls,
                        None,
                        &mut endpoint,
                        &source,
                        &vec![0.8; CELLS],
                        &mut workspace,
                        [0.0; 3],
                    )
                    .is_err()
            );
            assert_eq!(snapshot(&actual, &thermal), before);
            assert_eq!(
                endpoint
                    .iter()
                    .flatten()
                    .map(|v| v.to_bits())
                    .collect::<Vec<_>>(),
                endpoint_before
            );
            if failure < 3 {
                endpoint.fill([0.0001, 0.0, -0.0002]);
                source.fill(0.0001);
                actual
                    .advance_serial_with_thermal_sgs(
                        &mut thermal,
                        cx,
                        config(models()[0]),
                        None,
                        &mut endpoint,
                        &source,
                        &vec![0.8; CELLS],
                        &mut workspace,
                        [0.0; 3],
                    )
                    .unwrap();
                assert_eq!(endpoint, forces());
                assert_eq!(
                    [actual.completed_steps(), thermal.completed_steps()],
                    [1, 1]
                );
            }
        });
    }
}
