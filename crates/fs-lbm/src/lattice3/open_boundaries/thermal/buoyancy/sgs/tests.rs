use super::*;
use crate::d3q19::{BoundarySpec3, CollisionModel3, FaceBoundary3};
use crate::d3q27::{CollisionConfig27, CollisionModel27, D3q27, Rates27};
use crate::lattice3::{D3q19, Grid3, GridBudget3};
use crate::sgs3::{SgsConfig, SgsModel, SgsUse, eddy_viscosity};
use crate::thermal::boundary::FiniteThermalWall;
use crate::thermal::boundary::conduction::WallCell;
use crate::thermal::open::{ThermalOpenBoundaries3, ThermalOpenLimits3};
use fs_exec::CancelGate;

const DIMS: [usize; 3] = [6, 4, 3];
const CELLS: usize = 72;
const ORIGIN: [f64; 3] = [0.2, -0.4, 0.7];

fn with_cx<T>(f: impl FnOnce(&Cx<'_>) -> T) -> T {
    let gate = CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| f(&Cx::new(
        &gate, arena,
        fs_exec::StreamKey { seed: 0x5307_b0a7, kernel_id: 227, tile: 0, iteration: 0 },
        fs_exec::Budget::INFINITE, fs_exec::ExecMode::Deterministic,
    )))
}

fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:.17e} != {b:.17e}");
}

fn model19() -> CollisionModel3 {
    CollisionModel3::Trt { tau: 0.9, magic: 0.19 }
}

fn model27() -> CollisionConfig27 {
    CollisionConfig27 {
        model: CollisionModel27::Cumulant {
            rates: Rates27 { omega: [1.1, 1.2, 1.3, 1.1, 1.4, 1.2, 1.3, 1.2, 1.4, 1.1] },
        },
        maximum_mach: 0.2,
    }
}

fn config(model: SgsModel) -> OpenThermalSgsConfig3 {
    OpenThermalSgsConfig3 {
        closure: SgsConfig {
            model, filter_width_lu: 1.3, maximum_eddy_to_molecular: 20.0,
            usage: SgsUse::Exploration,
        },
        turbulent_prandtl: 0.73,
        maximum_diffusivity_to_molecular: 20.0,
    }
}

fn models() -> [SgsModel; 3] {
    [SgsModel::Smagorinsky { cs: 0.3 }, SgsModel::Wale { cw: 0.65 }, SgsModel::Vreman { cv: 0.15 }]
}

fn coordinates(cell: usize) -> [usize; 3] {
    [cell % DIMS[0], (cell / DIMS[0]) % DIMS[1], cell / (DIMS[0] * DIMS[1])]
}

fn velocity([x, y, z]: [usize; 3]) -> [f64; 3] {
    let [x, y, z] = [x as f64, y as f64, z as f64];
    [
        0.022 + 0.0004 * x + 0.00003 * x * x + 0.0002 * y,
        0.004 + 0.0008 * x + 0.00004 * x * x - 0.0003 * y + 0.0001 * z,
        -0.003 + 0.0005 * x - 0.00002 * x * x + 0.0004 * y - 0.0001 * z,
    ]
}

fn background(epoch: usize) -> Vec<[f64; 3]> {
    (0..CELLS).map(|cell| {
        let [x, y, z] = coordinates(cell);
        let sign = if epoch % 2 == 0 { 1.0 } else { -1.0 };
        [sign * 0.0001 * (1 + x) as f64, -0.00008 * (1 + y) as f64,
            sign * 0.00005 * (1 + z) as f64]
    }).collect()
}

fn flow<L: EquilibriumLattice3<Q>, const Q: usize>(
    cx: &Cx<'_>, collision: L::Collision,
) -> OpenBoundaryGrid3<L, Q> {
    let grid = Grid3::<L, Q>::try_new_equilibrium(
        DIMS, collision, [0.0001, -0.0002, 0.0003], [false, false, true],
        GridBudget3 { maximum_cells: CELLS, maximum_mach: 0.2 },
        |_| false, |p| (1.2, velocity(p)),
    ).unwrap();
    OpenBoundaryGrid3::try_new(grid, cx, BoundarySpec3::new([
        FaceBoundary3::Velocity { velocity: [0.025, 0.0, 0.0] },
        FaceBoundary3::Pressure { density: 1.2 },
        FaceBoundary3::stationary_wall(), FaceBoundary3::stationary_wall(),
        FaceBoundary3::Periodic, FaceBoundary3::Periodic,
    ])).unwrap()
}

fn make<L: EquilibriumLattice3<Q>, const Q: usize>(
    cx: &Cx<'_>, collision: L::Collision,
) -> BoussinesqOpenGrid3<L, Q> {
    let flow = flow::<L, Q>(cx, collision);
    let initial: Vec<_> = (0..CELLS).map(|cell| {
        let [x, y, _] = coordinates(cell);
        0.45 + 0.03 * y as f64 + 0.005 * x as f64
    }).collect();
    let wall = |temperature| ThermalOpenFace3::Wall(ThermalWall::FiniteSolid(FiniteThermalWall {
        temperature, heat_capacity: 0.2, resistance: 0.3, heat_per_step: 0.0,
    }));
    let thermal = flow.prepare_thermal_field(
        cx, &initial, 0.83,
        ThermalOpenBoundaries3 {
            faces: [ThermalOpenFace3::Reservoir { temperature: 0.6 },
                ThermalOpenFace3::Reservoir { temperature: 0.4 }, wall(0.9), wall(0.4),
                ThermalOpenFace3::Periodic, ThermalOpenFace3::Periodic],
            internal_wall: ThermalWall::Adiabatic,
        },
        ThermalOpenLimits3 { maximum_cells: CELLS, maximum_storage_bytes: 2 * 1024 * 1024,
            temperature_range: [-2.0, 3.0], balance_tolerance: 1e-9 },
    ).unwrap();
    let mut coupled = BoussinesqOpenGrid3::try_new(
        flow, thermal, cx,
        BoussinesqLaw3::try_new(1.7, 0.5, 0.02, [0.0, -0.02, 0.0], 0.03).unwrap(),
        background(0),
        BoussinesqOpenLimits3 { maximum_cells: CELLS, maximum_force_storage_bytes: 1024 * 1024,
            force_consistency_tolerance: 1e-12 },
    ).unwrap();
    coupled.set_wall_conduction(1, &[
        WallConductionLink { a: WallCell { side: 0, link: 1 }, b: WallCell { side: 0, link: 2 }, conductance: 0.3 },
        WallConductionLink { a: WallCell { side: 0, link: 2 }, b: WallCell { side: 1, link: 2 }, conductance: 0.15 },
    ], WallConductionControls { max_nodes: 3, max_links: 2, max_factor_bytes: 1024 * 1024,
        balance_tolerance: 1e-10 }).unwrap();
    coupled
}

fn molecular<L: SgsLattice3<Q>, const Q: usize>(collision: L::Collision) -> (Vec<L::Collision>, Vec<f64>) {
    ((0..CELLS).map(|cell| L::with_shear_relaxation(collision, 0.83 + 0.01 * (cell % 4) as f64).unwrap()).collect(),
        (0..CELLS).map(|cell| 0.76 + 0.01 * (cell % 5) as f64).collect())
}

fn sources() -> Vec<f64> {
    (0..CELLS).map(|cell| 0.0001 * (1 + cell % 3) as f64).collect()
}

#[derive(Debug, PartialEq)]
struct Snapshot {
    flow: Vec<u64>, heat: Vec<u64>, forces: Vec<[u64; 3]>, background: Vec<[u64; 3]>,
    walls: [Vec<ThermalWall>; 6], steps: [usize; 2],
    balance: Option<ThermalOpenBalance3>, receipt: Option<BoussinesqOpenStep3>,
}

fn snapshot<L: EquilibriumLattice3<Q>, const Q: usize>(owner: &BoussinesqOpenGrid3<L, Q>) -> Snapshot {
    Snapshot {
        flow: owner.flow().populations().iter().flatten().map(|f| f.to_bits()).collect(),
        heat: owner.thermal().populations().iter().flatten().map(|f| f.to_bits()).collect(),
        forces: owner.flow().forces().iter().map(|f| f.map(f64::to_bits)).collect(),
        background: owner.background_forces().iter().map(|f| f.map(f64::to_bits)).collect(),
        walls: core::array::from_fn(|face| owner.thermal().wall_states(face).unwrap().to_vec()),
        steps: [owner.flow().completed_steps(), owner.thermal().completed_steps()],
        balance: owner.thermal().last_balance(), receipt: owner.last_step(),
    }
}

fn energy(thermal: &ThermalOpenField3) -> f64 {
    thermal.populations().iter().flatten().sum::<f64>() + (0..6)
        .flat_map(|face| thermal.wall_states(face).unwrap())
        .map(|wall| match wall { ThermalWall::FiniteSolid(s) => s.heat_capacity * s.temperature, _ => 0.0 }).sum::<f64>()
}

fn force_case<L: EquilibriumLattice3<Q> + SgsLattice3<Q>, const Q: usize>(collision: L::Collision) {
    with_cx(|cx| {
        for model in models() {
            let cfg = config(model).closure;
            let mut actual = ForcedOpenBoundaryGrid3::try_new_preserving_velocity(
                flow::<L, Q>(cx, collision), cx, background(0),
            ).unwrap();
            let mut reference = ForcedOpenBoundaryGrid3::try_new_preserving_velocity(
                flow::<L, Q>(cx, collision), cx, background(0),
            ).unwrap();
            let bytes = actual.required_sgs_workspace_bytes().unwrap();
            assert!(actual.sgs_workspace(bytes - 1).is_err());
            let mut workspace = actual.sgs_workspace(bytes).unwrap();
            let mut native = reference.sgs_workspace(bytes).unwrap();
            let fit = reference.prepare_sgs(cx, cfg, None, &mut native).unwrap();
            for x in [0, 2, DIMS[0] - 1] {
                let mut gradient = [
                    [0.0004 + 0.00006 * x as f64, 0.0002, 0.0],
                    [0.0008 + 0.00008 * x as f64, -0.0003, 0.0001],
                    [0.0005 - 0.00004 * x as f64, 0.0004, -0.0001],
                ];
                if x == DIMS[0] - 1 { for row in &mut gradient { row[0] = 0.0; } }
                let cell = x + DIMS[0] * (1 + DIMS[1]);
                let expected = eddy_viscosity(model, cfg.filter_width_lu, gradient).unwrap();
                near(native.eddy_viscosities()[cell], expected, 3e-15);
                let (_, velocity) = actual.macroscopic(cell).unwrap();
                for (a, b) in velocity.into_iter().zip(self::velocity([x, 1, 1])) { near(a, b, 2e-16); }
            }
            let expected = reference.advance_serial(cx, Some(native.collisions()), background(1), ORIGIN).unwrap();
            let result = actual.advance_serial_with_sgs(cx, cfg, None, background(1), &mut workspace, ORIGIN).unwrap();
            assert_eq!(result.sgs, fit);
            assert_eq!(result.momentum, expected);
            assert_eq!(actual.populations(), reference.populations());
            assert_eq!(actual.forces(), reference.forces());
            let saved = (actual.populations().to_vec(), actual.forces().to_vec(), actual.completed_steps());
            let mut invalid = background(2);
            invalid[CELLS - 1][0] = f64::NAN;
            assert!(actual.advance_serial_with_sgs(cx, cfg, None, invalid, &mut workspace, ORIGIN).is_err());
            assert_eq!((actual.populations().to_vec(), actual.forces().to_vec(), actual.completed_steps()), saved);
            // Refresh the original molecular baseline at the new physical
            // state; the previously augmented field is never reused as it.
            let fit = reference.prepare_sgs(cx, cfg, None, &mut native).unwrap();
            let expected = reference.step_serial(cx, Some(native.collisions()), ORIGIN).unwrap();
            let result = actual.step_serial_with_sgs(cx, cfg, None, &mut workspace, ORIGIN).unwrap();
            assert_eq!(result.sgs, fit);
            assert_eq!(result.momentum, expected);
            assert_eq!(actual.populations(), reference.populations());
            assert_eq!(actual.forces(), reference.forces());
        }
    });
}

#[test]
fn native_open_sgs_uses_accepted_spatial_force_and_keeps_endpoint_force_separate() {
    force_case::<D3q19, 19>(model19());
    force_case::<D3q27, 27>(model27());
}

fn zero_case<L: EquilibriumLattice3<Q> + SgsLattice3<Q>, const Q: usize>(collision: L::Collision)
where L::Collision: core::fmt::Debug + PartialEq {
    with_cx(|cx| {
        let (laws, tau) = molecular::<L, Q>(collision);
        let heat = sources();
        for model in [SgsModel::Smagorinsky { cs: 0.0 }, SgsModel::Wale { cw: 0.0 }, SgsModel::Vreman { cv: 0.0 }] {
            let mut actual = make::<L, Q>(cx, collision);
            let mut reference = make::<L, Q>(cx, collision);
            let bytes = actual.required_thermal_sgs_workspace_bytes().unwrap();
            assert!(actual.thermal_sgs_workspace(bytes - 1).is_err());
            let mut workspace = actual.thermal_sgs_workspace(bytes).unwrap();
            for epoch in 1..=3 {
                let mut next = background(epoch);
                let mut expected_next = next.clone();
                let result = actual.advance_serial_with_sgs(cx, config(model), Some(&laws), &mut next,
                    &heat, &tau, &mut workspace, ORIGIN).unwrap();
                let expected = reference.advance_serial(cx, Some(&laws), &mut expected_next,
                    &heat, Some(&tau), ORIGIN).unwrap();
                assert_eq!(snapshot(&actual), snapshot(&reference));
                assert_eq!(result.transport, expected);
                assert_eq!(next, expected_next);
                assert_eq!(workspace.momentum_collisions(), laws);
                assert_eq!(workspace.thermal_relaxation(), tau);
                assert_eq!(result.sgs.maximum_eddy_viscosity_lu, 0.0);
                assert_eq!(result.maximum_eddy_diffusivity_lu, 0.0);
            }
            let result = actual.step_serial_with_sgs(cx, config(model), Some(&laws), &heat, &tau, &mut workspace, ORIGIN).unwrap();
            let expected = reference.step_serial(cx, Some(&laws), &heat, Some(&tau), ORIGIN).unwrap();
            assert_eq!(result.transport, expected);
            assert_eq!(snapshot(&actual), snapshot(&reference));
        }
    });
}

#[test]
fn zero_sgs_preserves_joint_source_wall_force_and_thermal_trajectories_bitwise() {
    zero_case::<D3q19, 19>(model19());
    zero_case::<D3q27, 27>(model27());
}

fn active_case<L: EquilibriumLattice3<Q> + SgsLattice3<Q>, const Q: usize>(collision: L::Collision) {
    with_cx(|cx| {
        let (laws, tau) = molecular::<L, Q>(collision);
        let heat = sources();
        for model in models() {
            let cfg = config(model);
            let mut actual = make::<L, Q>(cx, collision);
            let mut reference = make::<L, Q>(cx, collision);
            let mut molecular = make::<L, Q>(cx, collision);
            for owner in [&mut actual, &mut reference, &mut molecular] {
                owner.step_serial(cx, Some(&laws), &heat, Some(&tau), ORIGIN).unwrap();
            }
            let mut workspace = actual.thermal_sgs_workspace(1024 * 1024).unwrap();
            let mut native = reference.flow().sgs_workspace(1024 * 1024).unwrap();
            for epoch in 1..=3 {
                let before_energy = energy(actual.thermal());
                let opening = actual.flow().forces().to_vec();
                let fit = reference.flow().prepare_sgs(cx, cfg.closure, Some(&laws), &mut native).unwrap();
                let total_tau: Vec<_> = tau.iter().zip(native.eddy_viscosities()).map(|(&tau, &eddy)| {
                    if eddy == 0.0 { tau } else { tau + 4.0 * (eddy / cfg.turbulent_prandtl) }
                }).collect();
                let mut next = background(epoch);
                let mut expected_next = next.clone();
                let mut molecular_next = next.clone();
                let result = actual.advance_serial_with_sgs(cx, cfg, Some(&laws), &mut next,
                    &heat, &tau, &mut workspace, ORIGIN).unwrap();
                let expected = reference.advance_serial(cx, Some(native.collisions()), &mut expected_next,
                    &heat, Some(&total_tau), ORIGIN).unwrap();
                molecular.advance_serial(cx, Some(&laws), &mut molecular_next, &heat, Some(&tau), ORIGIN).unwrap();
                assert_eq!(result.transport, expected);
                assert_eq!(result.sgs, fit);
                assert_eq!(snapshot(&actual), snapshot(&reference));
                assert_eq!(next, expected_next);
                assert_eq!(workspace.thermal_relaxation(), total_tau);
                assert!(result.sgs.sgs_dissipation_lu > 0.0);
                near(result.maximum_eddy_diffusivity_lu, fit.maximum_eddy_viscosity_lu / cfg.turbulent_prandtl, 2e-16);
                let reservoir = result.transport.heat.reservoirs.iter().map(|f| f.heat_in).sum::<f64>();
                near(energy(actual.thermal()) - before_energy, heat.iter().sum::<f64>() + reservoir, 2e-10);
                assert!(result.transport.heat.combined_defect.abs() < 1e-9);
                for axis in 0..3 {
                    let impulse = opening.iter().zip(actual.flow().forces()).map(|(a, b)| 0.5 * a[axis] + 0.5 * b[axis]).sum();
                    near(result.transport.body_impulse[axis], impulse, 2e-15);
                    let mechanical = next.iter().zip(actual.background_forces()).map(|(a, b)| 0.5 * a[axis] + 0.5 * b[axis]).sum();
                    near(result.transport.background_impulse[axis], mechanical, 2e-15);
                    near(result.transport.body_impulse[axis], result.transport.background_impulse[axis] + result.transport.buoyancy_impulse[axis], 2e-15);
                }
            }
            assert_ne!(actual.flow().populations(), molecular.flow().populations());
            assert!(actual.thermal().populations().iter().flatten().zip(molecular.thermal().populations().iter().flatten()).any(|(a, b)| (a - b).abs() > 1e-14));
            assert!(actual.flow().forces().iter().flatten().zip(molecular.flow().forces().iter().flatten()).any(|(a, b)| (a - b).abs() > 1e-14));
            assert_ne!(actual.thermal().wall_states(2), molecular.thermal().wall_states(2));
        }
    });
}

#[test]
fn local_sgs_changes_native_buoyant_feedback_and_closes_independent_heat_and_impulse_accounts() {
    active_case::<D3q19, 19>(model19());
    active_case::<D3q27, 27>(model27());
}

#[test]
fn warmed_late_refusals_preserve_all_epochs_and_endpoint_background_then_replay_exactly() {
    with_cx(|cx| {
        let (laws, tau) = molecular::<D3q19, 19>(model19());
        let cfg = config(SgsModel::Wale { cw: 0.65 });
        let heat = sources();
        let mut actual = make::<D3q19, 19>(cx, model19());
        let mut reference = make::<D3q19, 19>(cx, model19());
        let mut workspace = actual.thermal_sgs_workspace(1024 * 1024).unwrap();
        let mut reference_workspace = reference.thermal_sgs_workspace(1024 * 1024).unwrap();
        actual.step_serial_with_sgs(cx, cfg, Some(&laws), &heat, &tau, &mut workspace, ORIGIN).unwrap();
        reference.step_serial_with_sgs(cx, cfg, Some(&laws), &heat, &tau, &mut reference_workspace, ORIGIN).unwrap();
        let before = snapshot(&actual);
        let mut endpoint = background(1);
        let saved_endpoint = endpoint.clone();
        let mut reference_endpoint = endpoint.clone();
        let mut polls = 0;
        let expected = reference.advance_sgs_with_poll(Some(&mut reference_endpoint), cfg, Some(&laws),
            &heat, &tau, &mut reference_workspace, ORIGIN, &mut || { polls += 1; Ok(()) }).unwrap();
        let mut calls = 0;
        let cancelled = actual.advance_sgs_with_poll(Some(&mut endpoint), cfg, Some(&laws),
            &heat, &tau, &mut workspace, ORIGIN, &mut || {
                calls += 1;
                if calls == polls { Err(GridError3::Cancelled) } else { Ok(()) }
            });
        assert!(matches!(cancelled, Err(BoussinesqOpenSgsError3::Transport(BoussinesqOpenError3::Flow(GridError3::Cancelled)))));
        assert_eq!(calls, polls);
        assert_eq!(snapshot(&actual), before);
        assert_eq!(endpoint, saved_endpoint);
        let mut limited = cfg;
        limited.maximum_diffusivity_to_molecular = 0.0;
        assert!(actual.advance_serial_with_sgs(cx, limited, Some(&laws), &mut endpoint,
            &heat, &tau, &mut workspace, ORIGIN).is_err());
        assert_eq!(snapshot(&actual), before);
        assert_eq!(endpoint, saved_endpoint);
        assert!(actual.advance_serial_with_sgs(cx, cfg, Some(&laws), &mut endpoint,
            &vec![2.2; CELLS], &tau, &mut workspace, ORIGIN).is_err());
        assert_eq!(snapshot(&actual), before);
        assert_eq!(endpoint, saved_endpoint);
        let retried = actual.advance_serial_with_sgs(cx, cfg, Some(&laws), &mut endpoint,
            &heat, &tau, &mut workspace, ORIGIN).unwrap();
        assert_eq!(retried.transport, expected.transport);
        assert_eq!(retried.sgs, expected.sgs);
        assert_eq!(snapshot(&actual), snapshot(&reference));
        assert_eq!(endpoint, reference_endpoint);
    });
}
