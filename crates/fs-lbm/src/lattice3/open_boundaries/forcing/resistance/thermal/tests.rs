//! Physical velocity, original heat transport and joint porous rollback.

use super::*;
use crate::d3q19::{BoundarySpec3, CollisionModel3, Face3, FaceBoundary3};
use crate::d3q27::{CollisionConfig27, CollisionModel27, D3q27, Equilibrium27};
use crate::lattice3::open_boundaries::OpenBoundaryGrid3;
use crate::lattice3::open_boundaries::forcing::ForcedOpenBoundaryGrid3;
use crate::lattice3::open_boundaries::forcing::resistance::{
    DarcyForchheimer3, LinearResistance3,
};
use crate::lattice3::{D3q19, Grid3, GridBudget3};
use crate::thermal::boundary::{FiniteThermalWall, ThermalWall};
use crate::thermal::open::ThermalOpenFace3;
use crate::thermal::tensor::DiffusivityTensor;
use fs_exec::CancelGate;

const RHO: f64 = 1.2;

mod endpoint;

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        f(&Cx::new(
            &gate, arena,
            fs_exec::StreamKey { seed: 0xDA2C_7E47, kernel_id: 227, tile: 0, iteration: 0 },
            fs_exec::Budget::INFINITE, fs_exec::ExecMode::Deterministic,
        ))
    })
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!((actual - expected).abs() <= tolerance, "{actual:.17e} != {expected:.17e}");
}

fn bgk27() -> CollisionConfig27 {
    CollisionConfig27 {
        model: CollisionModel27::Bgk { tau: 0.9, equilibrium: Equilibrium27::Factorized },
        maximum_mach: 0.2,
    }
}

fn limits() -> ThermalOpenLimits3 {
    ThermalOpenLimits3 {
        maximum_cells: 512, maximum_storage_bytes: 2 * 1024 * 1024,
        temperature_range: [0.0, 5.0], balance_tolerance: 1e-10,
    }
}

fn periodic_heat() -> ThermalOpenBoundaries3 {
    ThermalOpenBoundaries3 {
        faces: [ThermalOpenFace3::Periodic; 6], internal_wall: ThermalWall::Adiabatic,
    }
}

fn native<L: EquilibriumLattice3<Q>, const Q: usize>(
    cx: &Cx<'_>,
    dims: [usize; 3],
    collision: L::Collision,
    boundaries: BoundarySpec3,
    velocity: impl Fn([usize; 3]) -> [f64; 3],
) -> OpenBoundaryGrid3<L, Q> {
    let periodic = core::array::from_fn(|axis| {
        matches!(boundaries.face(Face3::ALL[2 * axis]), FaceBoundary3::Periodic)
    });
    let grid = Grid3::<L, Q>::try_new_equilibrium(
        dims, collision, [-0.0007, 0.0002, 0.0003], periodic,
        GridBudget3 { maximum_cells: 512, maximum_mach: 0.2 },
        |_| false, |p| (RHO, velocity(p)),
    ).unwrap();
    OpenBoundaryGrid3::try_new(grid, cx, boundaries).unwrap()
}

fn relaxation_and_heat<L: EquilibriumLattice3<Q>, const Q: usize>(collision: L::Collision) {
    with_cx(|cx| {
        let cells = 24;
        let u0 = [0.04, -0.02, 0.03];
        let skeleton = [0.01, 0.005, -0.002];
        let resistance = 0.24;
        let law = DarcyForchheimer3::from(LinearResistance3 {
            principal: [resistance; 3], skeleton_velocity: skeleton,
            ..LinearResistance3::default()
        });
        let mut owner = ResistanceOpenGrid3::try_new_preserving_velocity(
            native::<L, Q>(cx, [4, 3, 2], collision, BoundarySpec3::periodic(), |_| u0),
            cx, vec![law; cells], [0.0; 3],
        ).unwrap();
        let mut heat = owner.prepare_thermal_field(
            cx, &vec![0.7; cells], 0.83, periodic_heat(), limits(),
        ).unwrap();
        let increment = 0.004;
        let tensor = DiffusivityTensor::from_lattice([
            [0.08, 0.015, 0.0], [0.015, 0.1, -0.01], [0.0, -0.01, 0.09],
        ]).unwrap();
        let ratio = (RHO - 0.5 * resistance) / (RHO + 0.5 * resistance);
        for step in 1_i32..=8 {
            let r = if step % 2 == 0 {
                owner.step_serial_with_thermal_tensor_diffusion(
                    &mut heat, cx, None, &vec![increment; cells],
                    DiffusivityField::Uniform(&tensor),
                ).unwrap()
            } else {
                owner.step_serial_with_thermal(
                    &mut heat, cx, None, &vec![increment; cells], None,
                ).unwrap()
            };
            let expected: [f64; 3] = core::array::from_fn(|a| {
                skeleton[a] + ratio.powi(step) * (u0[a] - skeleton[a])
            });
            for cell in 0..cells {
                let (rho, u) = owner.flow().macroscopic(cell).unwrap();
                near(rho, RHO, 3e-13);
                for a in 0..3 {
                    near(u[a], expected[a], 3e-13);
                    near(owner.flow().forces()[cell][a], resistance * (skeleton[a] - expected[a]), 3e-13);
                    if step == 1 {
                        let (plus, minus) = crate::thermal::advection::D3Q7.direction_pair(a);
                        let g = heat.populations()[cell];
                        let expected_flux = (0.7 + 0.5 * increment) * u0[a]
                            + 0.5 * increment * expected[a];
                        near(g[plus] - g[minus], expected_flux, 4e-15);
                    }
                }
                near(heat.temperature(cell).unwrap(), 0.7 + f64::from(step) * increment, 3e-13);
            }
            near(r.heat.source, cells as f64 * increment, 2e-15);
            near(r.heat.stored_change, cells as f64 * increment, 3e-13);
            near(r.heat.defect, 0.0, 3e-13);
            assert!(r.momentum.relative_dissipation > 0.0);
            // Mechanical dissipation is a separate receipt, not extra heat.
            assert_eq!(owner.last_exchange(), Some(r.momentum));
            assert_eq!(owner.flow().completed_steps(), step as usize);
            assert_eq!(heat.completed_steps(), step as usize);
        }
    });
}

#[test]
fn physical_drag_relaxation_and_both_heat_half_kicks_use_their_own_force_epoch() {
    relaxation_and_heat::<D3q19, 19>(CollisionModel3::Bgk { tau: 0.9 });
    relaxation_and_heat::<D3q27, 27>(bgk27());
}

#[test]
fn porous_heat_packet_advection_uses_actual_opening_velocity_in_each_cell() {
    with_cx(|cx| {
        let dims = [9, 4, 2];
        let cells = 72;
        let velocity = |y: usize| 0.01 + 0.008 * y as f64;
        let laws: Vec<_> = (0..cells).map(|cell| DarcyForchheimer3 {
            linear: LinearResistance3 {
                principal: [0.08 * (1 + cell / 9 % 4) as f64; 3],
                skeleton_velocity: [0.003, 0.0, 0.0],
                ..LinearResistance3::default()
            },
            quadratic: 0.8,
        }).collect();
        let mut owner = ResistanceOpenGrid3::try_new_preserving_velocity(
            native::<D3q19, 19>(cx, dims, CollisionModel3::Bgk { tau: 0.9 },
                BoundarySpec3::periodic(), |p| [velocity(p[1]), 0.0, 0.0]),
            cx, laws, [0.0; 3],
        ).unwrap();
        let initial: Vec<_> = (0..cells).map(|cell| {
            if cell % 9 == 4 { (1 + cell / 9 % 4) as f64 } else { 0.0 }
        }).collect();
        let mut heat = owner.prepare_thermal_field(
            cx, &initial, 0.83, periodic_heat(), limits(),
        ).unwrap();
        let expected_shift: f64 = initial.iter().enumerate().map(|(cell, t)| {
            t * velocity(cell / 9 % 4)
        }).sum();
        let before: f64 = initial.iter().enumerate().map(|(cell, t)| (cell % 9) as f64 * t).sum();
        let report = owner.step_serial_with_thermal(
            &mut heat, cx, None, &vec![0.0; cells], None,
        ).unwrap();
        let after: f64 = (0..cells).map(|cell| {
            (cell % 9) as f64 * heat.temperature(cell).unwrap()
        }).sum();
        // The packet does not cross a seam; diffusion has zero first moment.
        near(after - before, expected_shift, 4e-13);
        near(report.heat.stored_change, 0.0, 4e-13);
        near(report.heat.defect, 0.0, 4e-13);
    });
}

fn zero_law_parity<L: EquilibriumLattice3<Q>, const Q: usize>(collision: L::Collision) {
    with_cx(|cx| {
        let cells = 48;
        let boundaries = BoundarySpec3::new([
            FaceBoundary3::Velocity { velocity: [0.025, 0.0, 0.0] },
            FaceBoundary3::Pressure { density: RHO },
            FaceBoundary3::stationary_wall(), FaceBoundary3::stationary_wall(),
            FaceBoundary3::Periodic, FaceBoundary3::Periodic,
        ]);
        let make = || native::<L, Q>(cx, [6, 4, 2], collision, boundaries, |_| [0.018, 0.0, 0.0]);
        let forces: Vec<_> = (0..cells).map(|cell| [1e-5 * (1 + cell % 6) as f64, 0.0, 0.0]).collect();
        let laws = forces.iter().map(|&body_force| {
            DarcyForchheimer3::from(LinearResistance3 { body_force, ..LinearResistance3::default() })
        }).collect();
        let mut porous = ResistanceOpenGrid3::try_new_preserving_velocity(make(), cx, laws, [0.0; 3]).unwrap();
        let mut forced = ForcedOpenBoundaryGrid3::try_new_preserving_velocity(make(), cx, forces).unwrap();
        let thermal_boundaries = ThermalOpenBoundaries3 {
            faces: [
                ThermalOpenFace3::Reservoir { temperature: 0.6 },
                ThermalOpenFace3::Reservoir { temperature: 0.9 },
                ThermalOpenFace3::Wall(ThermalWall::Adiabatic),
                ThermalOpenFace3::Wall(ThermalWall::Adiabatic),
                ThermalOpenFace3::Periodic, ThermalOpenFace3::Periodic,
            ],
            internal_wall: ThermalWall::Adiabatic,
        };
        let mut porous_heat = porous.prepare_thermal_field(cx, &vec![0.5; cells], 0.83, thermal_boundaries, limits()).unwrap();
        let mut forced_heat = forced.prepare_thermal_field(cx, &vec![0.5; cells], 0.83, thermal_boundaries, limits()).unwrap();
        let source: Vec<_> = (0..cells).map(|cell| 0.0001 * (1 + cell % 5) as f64).collect();
        let tensor = DiffusivityTensor::from_lattice([
            [0.08, 0.015, 0.0], [0.015, 0.1, -0.01], [0.0, -0.01, 0.09],
        ]).unwrap();
        let local_tensors = vec![tensor; cells];
        for step in 0..4 {
            let (a, b) = if step % 2 == 0 {
                (
                    porous.step_serial_with_thermal(&mut porous_heat, cx, None, &source, None).unwrap(),
                    forced.step_serial_with_thermal(&mut forced_heat, cx, None, &source, None, [0.0; 3]).unwrap(),
                )
            } else {
                (
                    porous.step_serial_with_thermal_tensor_diffusion(
                        &mut porous_heat, cx, None, &source, DiffusivityField::Uniform(&tensor),
                    ).unwrap(),
                    forced.step_serial_with_thermal_tensor_diffusion(
                        &mut forced_heat, cx, None, &source, DiffusivityField::PerCell(&local_tensors), [0.0; 3],
                    ).unwrap(),
                )
            };
            assert_eq!(porous.flow().populations().iter().flatten().map(|v| v.to_bits()).collect::<Vec<_>>(),
                forced.populations().iter().flatten().map(|v| v.to_bits()).collect::<Vec<_>>());
            assert_eq!(porous.flow().forces(), forced.forces());
            assert_eq!(porous_heat.populations(), forced_heat.populations());
            assert_eq!(a.momentum.open, b.momentum);
            assert_eq!(a.heat, b.heat);
        }
    });
}

#[test]
fn zero_law_preserves_original_open_heat_with_scalar_and_oblique_tensor_diffusion() {
    zero_law_parity::<D3q19, 19>(CollisionModel3::Bgk { tau: 0.9 });
    zero_law_parity::<D3q27, 27>(bgk27());
}

fn channel(cx: &Cx<'_>) -> (ResistanceOpenGrid3<D3q19, 19>, ThermalOpenField3) {
    let cells = 48;
    let boundaries = BoundarySpec3::new([
        FaceBoundary3::Velocity { velocity: [0.025, 0.0, 0.0] },
        FaceBoundary3::Pressure { density: RHO },
        FaceBoundary3::stationary_wall(), FaceBoundary3::stationary_wall(),
        FaceBoundary3::Periodic, FaceBoundary3::Periodic,
    ]);
    let laws = (0..cells).map(|cell| {
        let (x, y) = (cell % 6, cell / 6 % 4);
        let seam = (x == 0 || x == 5) && (y == 0 || y == 3);
        DarcyForchheimer3 {
            linear: LinearResistance3 {
                principal: if seam { [0.0; 3] } else { [0.04; 3] },
                skeleton_velocity: [0.004, 0.0, 0.0],
                body_force: [0.00002, 0.0, 0.0],
                ..LinearResistance3::default()
            },
            quadratic: if seam { 0.0 } else { 0.5 },
        }
    }).collect();
    let mut owner = ResistanceOpenGrid3::try_new_preserving_velocity(
        native::<D3q19, 19>(cx, [6, 4, 2], CollisionModel3::Bgk { tau: 0.9 },
            boundaries, |_| [0.015, 0.0, 0.0]),
        cx, laws, [0.5, -0.25, 1.0],
    ).unwrap();
    let wall = ThermalOpenFace3::Wall(ThermalWall::FiniteSolid(FiniteThermalWall {
        temperature: 0.7, heat_capacity: 1.5, resistance: 0.4, heat_per_step: 0.002,
    }));
    let mut heat = owner.prepare_thermal_field(
        cx, &vec![0.5; cells], 0.83,
        ThermalOpenBoundaries3 {
            faces: [
                ThermalOpenFace3::Reservoir { temperature: 0.5 },
                ThermalOpenFace3::Reservoir { temperature: 0.9 }, wall, wall,
                ThermalOpenFace3::Periodic, ThermalOpenFace3::Periodic,
            ], internal_wall: ThermalWall::Adiabatic,
        },
        ThermalOpenLimits3 { temperature_range: [0.0, 1.0], ..limits() },
    ).unwrap();
    owner.step_serial_with_thermal(&mut heat, cx, None, &vec![0.0001; cells], None).unwrap();
    (owner, heat)
}

#[derive(Debug, PartialEq)]
struct Accepted {
    populations: Vec<[u64; 19]>,
    forces: Vec<[u64; 3]>,
    drag: Vec<[u64; 3]>,
    heat: Vec<[u64; 7]>,
    walls: [Vec<ThermalWall>; 6],
    time: [usize; 2],
    momentum: Option<ResistanceMomentumExchange3>,
    endpoint: Option<ResistanceOpenEndpointExchange3>,
    laws: Vec<DarcyForchheimer3>,
    thermal: Option<ThermalOpenBalance3>,
    boundaries: BoundarySpec3,
}
fn accepted(owner: &ResistanceOpenGrid3<D3q19, 19>, heat: &ThermalOpenField3) -> Accepted {
    Accepted {
        populations: owner.flow().populations().iter().map(|v| v.map(f64::to_bits)).collect(),
        forces: owner.flow().forces().iter().map(|v| v.map(f64::to_bits)).collect(),
        drag: owner.drag.iter().map(|v| v.map(f64::to_bits)).collect(),
        heat: heat.populations().iter().map(|v| v.map(f64::to_bits)).collect(),
        walls: core::array::from_fn(|face| heat.wall_states(face).unwrap().to_vec()),
        time: [owner.flow().completed_steps(), heat.completed_steps()],
        momentum: owner.last_exchange(), endpoint: owner.last_endpoint_exchange(),
        laws: owner.laws().to_vec(), thermal: heat.last_balance(),
        boundaries: owner.flow().boundaries(),
    }
}

#[test]
fn late_heat_refusal_or_cancellation_preserves_drag_finite_walls_receipts_and_retry() {
    with_cx(|cx| {
        let (mut expected, mut expected_heat) = channel(cx);
        let source = vec![0.001; 48];
        let mut polls = 0;
        let expected_receipt = expected.step_with_thermal_poll(
            &mut expected_heat, None, &source, None, None, &mut || { polls += 1; Ok(()) },
        ).unwrap();
        for cancel in [false, true] {
            let (mut owner, mut heat) = channel(cx);
            let before = accepted(&owner, &heat);
            let pending_before = owner.flow.pending_forces.clone();
            let mut seen = 0;
            // Each opening half source remains strictly inside the envelope;
            // the closing half must refuse after porous momentum was staged.
            let excessive_heat: Vec<_> = (0..48).map(|cell| {
                1.8 * (1.0 - heat.temperature(cell).unwrap())
            }).collect();
            let refused = owner.step_with_thermal_poll(
                &mut heat, None, if cancel { &source } else { &excessive_heat },
                None, None, &mut || {
                    seen += 1;
                    if cancel && seen == polls { Err(GridError3::Cancelled) } else { Ok(()) }
                },
            );
            if cancel {
                assert!(matches!(refused, Err(OpenThermalError3::Flow(GridError3::Cancelled))));
            } else {
                assert!(matches!(refused, Err(OpenThermalError3::Thermal(ThermalOpenError3::Temperature { .. }))));
                assert_ne!(owner.flow.pending_forces, pending_before, "native momentum was staged before thermal refusal");
            }
            assert_eq!(accepted(&owner, &heat), before);
            let retry = owner.step_serial_with_thermal(&mut heat, cx, None, &source, None).unwrap();
            assert_eq!(retry, expected_receipt);
            assert_eq!(accepted(&owner, &heat), accepted(&expected, &expected_heat));
            let mut storage = retry.heat.stored_change;
            let mut supplied = retry.heat.source;
            for reservoir in retry.heat.reservoirs { supplied += reservoir.heat_in; }
            for wall in retry.heat.walls.into_iter().flatten() {
                let solid = wall.solid.unwrap();
                storage += solid.stored_change;
                supplied += solid.source;
            }
            near(storage, supplied, 2e-12);
            near(retry.heat.combined_defect, 0.0, 2e-12);
        }
    });
}

#[test]
fn original_time_and_moving_finite_wall_refusals_apply_before_publication() {
    with_cx(|cx| {
        let (mut owner, mut heat) = channel(cx);
        owner.step_serial(cx, None).unwrap();
        let before = accepted(&owner, &heat);
        assert!(matches!(
            owner.step_serial_with_thermal(&mut heat, cx, None, &vec![0.0; 48], None),
            Err(OpenThermalError3::Thermal(ThermalOpenError3::Invalid(_)))
        ));
        assert_eq!(accepted(&owner, &heat), before);

        let (mut owner, mut heat) = channel(cx);
        let mut faces = Face3::ALL.map(|face| owner.flow().boundaries().face(face));
        faces[2] = FaceBoundary3::Wall { velocity: [0.0, 0.0, 0.01] };
        faces[3] = faces[2];
        owner.rebind_open_boundaries(BoundarySpec3::new(faces), &CancelGate::new_clock_free()).unwrap();
        let before = accepted(&owner, &heat);
        assert!(matches!(
            owner.step_serial_with_thermal(&mut heat, cx, None, &vec![0.0; 48], None),
            Err(OpenThermalError3::Thermal(ThermalOpenError3::Invalid(_)))
        ));
        assert_eq!(accepted(&owner, &heat), before);
    });
}
