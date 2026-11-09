use super::*;
use crate::d3q19::{BoundarySpec3, CollisionModel3, FaceBoundary3};
use crate::d3q27::{CollisionConfig27, CollisionModel27, D3q27, Equilibrium27};
use crate::lattice3::{D3q19, Grid3, GridBudget3};
use crate::lattice3::open_boundaries::OpenBoundaryGrid3;
use crate::thermal::boundary::{FiniteThermalWall, ThermalWall};
use crate::thermal::open::ThermalOpenFace3;
use crate::thermal::tensor::DiffusivityTensor;

const RHO: f64 = 1.2;

fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:.17e} != {b:.17e}");
}

fn limits() -> ThermalOpenLimits3 {
    ThermalOpenLimits3 {
        maximum_cells: 512,
        maximum_storage_bytes: 2 * 1024 * 1024,
        temperature_range: [0.0, 5.0],
        balance_tolerance: 1e-10,
    }
}

fn periodic_heat() -> ThermalOpenBoundaries3 {
    ThermalOpenBoundaries3 {
        faces: [ThermalOpenFace3::Periodic; 6],
        internal_wall: ThermalWall::Adiabatic,
    }
}

fn flow<L: EquilibriumLattice3<Q>, const Q: usize>(
    dims: [usize; 3],
    collision: L::Collision,
    boundaries: BoundarySpec3,
    forces: Vec<[f64; 3]>,
    velocity: impl Fn([usize; 3]) -> (f64, [f64; 3]),
) -> ForcedOpenBoundaryGrid3<L, Q> {
    let periodic = core::array::from_fn(|axis| {
        matches!(boundaries.face(crate::d3q19::Face3::ALL[2 * axis]), FaceBoundary3::Periodic)
    });
    let grid = Grid3::<L, Q>::try_new_equilibrium(
        dims, collision, [-0.0007, 0.0002, 0.0003], periodic,
        GridBudget3 { maximum_cells: 512, maximum_mach: 0.2 },
        |_| false, velocity,
    ).unwrap();
    let native = OpenBoundaryGrid3::try_new_with_poll(grid, boundaries, &mut || Ok(())).unwrap();
    ForcedOpenBoundaryGrid3::try_new_preserving_velocity_with_poll(
        native, forces, &mut || Ok(()),
    ).unwrap()
}

fn affine_acceleration<L: EquilibriumLattice3<Q>, const Q: usize>(collision: L::Collision) {
    let cells = 24;
    let u0 = [0.02, -0.01, 0.015];
    let f0 = [0.001, -0.0003, 0.0002];
    let slope = [-0.00007, 0.00002, 0.00001];
    let mut flow = flow::<L, Q>(
        [4, 3, 2], collision, BoundarySpec3::periodic(), vec![f0; cells], |_| (RHO, u0),
    );
    let mut heat = flow.prepare_thermal_field_with_poll(
        &vec![0.7; cells], 0.83, periodic_heat(), limits(), &mut || Ok(()),
    ).unwrap();
    let source = 0.004;
    for step in 1..=12 {
        let time = step as f64;
        let next = core::array::from_fn(|axis| f0[axis] + time * slope[axis]);
        let mut endpoint = vec![next; cells];
        let previous_force = flow.forces().to_vec();
        let report = flow.step_with_thermal_poll(
            &mut heat, Some(&mut endpoint), None, &vec![source; cells], None, None,
            [0.0; 3], &mut || Ok(()),
        ).unwrap();
        assert_eq!(endpoint, previous_force, "accepted buffer is reused only after success");
        assert_eq!(flow.forces(), vec![next; cells]);
        let expected_u: [f64; 3] = core::array::from_fn(|axis| {
            u0[axis] + (time * f0[axis] + 0.5 * time * time * slope[axis]) / RHO
        });
        for cell in 0..cells {
            let (rho, velocity) = flow.macroscopic(cell).unwrap();
            near(rho, RHO, 2e-13);
            for axis in 0..3 {
                near(velocity[axis], expected_u[axis], 2e-13);
                if step == 1 {
                    let pair = crate::thermal::advection::D3Q7.direction_pair(axis);
                    let g = heat.populations()[cell];
                    // Uniform transport changes no local zeroth moment. The
                    // two prescribed heat halves carry opening and endpoint
                    // physical velocities, independently of collision choice.
                    let flux = (0.7 + 0.5 * source) * u0[axis]
                        + 0.5 * source * expected_u[axis];
                    near(g[pair.0] - g[pair.1], flux, 3e-15);
                }
            }
            near(heat.temperature(cell).unwrap(), 0.7 + time * source, 2e-13);
        }
        near(report.heat.source, cells as f64 * source, 2e-15);
        near(report.heat.stored_change, cells as f64 * source, 2e-13);
        near(report.heat.defect, 0.0, 2e-13);
        assert_eq!(flow.completed_steps(), step);
        assert_eq!(heat.completed_steps(), step);
    }
}

#[test]
fn g1_prescribed_acceleration_and_heat_use_both_physical_force_epochs() {
    affine_acceleration::<D3q19, 19>(CollisionModel3::Bgk { tau: 0.9 });
    affine_acceleration::<D3q27, 27>(CollisionConfig27 {
        model: CollisionModel27::Bgk { tau: 0.9, equilibrium: Equilibrium27::Factorized },
        maximum_mach: 0.2,
    });
}

#[test]
fn g1_heat_packet_displacement_uses_each_cells_actual_opening_velocity() {
    let dims = [9, 4, 2];
    let cells = dims.into_iter().product();
    let forces: Vec<_> = (0..cells).map(|cell| {
        [0.002 * (1 + cell / dims[0] % dims[1]) as f64, -0.0003, 0.0001]
    }).collect();
    let velocity = |y: usize| 0.01 + 0.008 * y as f64;
    let mut flow = flow::<D3q19, 19>(
        dims, CollisionModel3::Bgk { tau: 0.9 }, BoundarySpec3::periodic(), forces,
        |point| (RHO, [velocity(point[1]), 0.0, 0.0]),
    );
    let initial: Vec<_> = (0..cells).map(|cell| {
        if cell % dims[0] == 4 { (1 + cell / dims[0] % dims[1]) as f64 } else { 0.0 }
    }).collect();
    let mut heat = flow.prepare_thermal_field_with_poll(
        &initial, 0.83, periodic_heat(), limits(), &mut || Ok(()),
    ).unwrap();
    let expected_shift: f64 = initial.iter().enumerate().map(|(cell, t)| {
        t * velocity(cell / dims[0] % dims[1])
    }).sum();
    let before: f64 = initial.iter().enumerate().map(|(cell, t)| {
        (cell % dims[0]) as f64 * t
    }).sum();
    let report = flow.step_with_thermal_poll(
        &mut heat, None, None, &vec![0.0; cells], None, None, [0.0; 3], &mut || Ok(()),
    ).unwrap();
    // No packet reaches a periodic seam in one step; molecular diffusion has
    // zero global first moment. Only the physical advective displacement remains.
    let after: f64 = (0..cells).map(|cell| {
        (cell % dims[0]) as f64 * heat.temperature(cell).unwrap()
    }).sum();
    near(after - before, expected_shift, 3e-13);
    near(report.heat.stored_change, 0.0, 3e-14);
    near(report.heat.defect, 0.0, 3e-14);
}

fn channel() -> (ForcedOpenBoundaryGrid3<D3q19, 19>, ThermalOpenField3) {
    let cells = 48;
    let u = [0.025, 0.0, 0.0];
    let flow = flow::<D3q19, 19>(
        [6, 4, 2], CollisionModel3::Bgk { tau: 0.9 },
        BoundarySpec3::new([
            FaceBoundary3::Velocity { velocity: u }, FaceBoundary3::Pressure { density: RHO },
            FaceBoundary3::stationary_wall(), FaceBoundary3::stationary_wall(),
            FaceBoundary3::Periodic, FaceBoundary3::Periodic,
        ]),
        (0..cells).map(|cell| [1e-5 * (1 + cell % 6) as f64, 0.0, 0.0]).collect(),
        |_| (RHO, u),
    );
    let wall = ThermalOpenFace3::Wall(ThermalWall::FiniteSolid(FiniteThermalWall {
        temperature: 0.7, heat_capacity: 1.5, resistance: 0.4, heat_per_step: 0.002,
    }));
    let heat = flow.prepare_thermal_field_with_poll(
        &vec![0.5; cells], 0.83,
        ThermalOpenBoundaries3 {
            faces: [
                ThermalOpenFace3::Reservoir { temperature: 0.5 },
                ThermalOpenFace3::Reservoir { temperature: 0.9 }, wall, wall,
                ThermalOpenFace3::Periodic, ThermalOpenFace3::Periodic,
            ],
            internal_wall: ThermalWall::Adiabatic,
        },
        ThermalOpenLimits3 { temperature_range: [0.0, 1.0], ..limits() },
        &mut || Ok(()),
    ).unwrap();
    (flow, heat)
}

#[derive(Debug, PartialEq)]
struct Accepted {
    flow: Vec<[u64; 19]>,
    forces: Vec<[u64; 3]>,
    heat: Vec<[u64; 7]>,
    walls: [Vec<ThermalWall>; 6],
    time: [usize; 2],
    receipt: Option<crate::thermal::open::ThermalOpenBalance3>,
    wall_receipt: Option<Vec<u64>>,
}
fn accepted(flow: &ForcedOpenBoundaryGrid3<D3q19, 19>, heat: &ThermalOpenField3) -> Accepted {
    Accepted {
        flow: flow.populations().iter().map(|v| v.map(f64::to_bits)).collect(),
        forces: flow.forces().iter().map(|v| v.map(f64::to_bits)).collect(),
        heat: heat.populations().iter().map(|v| v.map(f64::to_bits)).collect(),
        walls: core::array::from_fn(|face| heat.wall_states(face).unwrap().to_vec()),
        time: [flow.completed_steps(), heat.completed_steps()],
        receipt: heat.last_balance(),
        wall_receipt: flow.last_wall_exchange().map(|r| {
            let mut bits: Vec<_> = [
                r.wall_impulse, r.fluid_population_impulse, r.wall_velocity_mass_impulse,
                r.wall_angular_impulse, r.fluid_population_angular_impulse,
                r.wall_velocity_mass_angular_impulse,
            ].into_iter().flatten().chain([r.fluid_mass_change, r.wall_work])
                .map(f64::to_bits).collect();
            bits.push(r.measured_links as u64);
            bits
        }),
    }
}

#[test]
fn g3_late_cancel_or_heat_refusal_preserves_forces_heat_and_finite_walls_for_exact_retry() {
    let (mut expected, mut expected_heat) = channel();
    let cells = expected.populations().len();
    let endpoint = vec![[0.0002, 0.0, -0.0001]; cells];
    let source = vec![0.001; cells];
    let mut expected_endpoint = endpoint.clone();
    let mut polls = 0;
    expected.step_with_thermal_poll(
        &mut expected_heat, Some(&mut expected_endpoint), None, &source, None, None,
        [0.0; 3], &mut || { polls += 1; Ok(()) },
    ).unwrap();
    for cancel in [false, true] {
        let (mut actual, mut actual_heat) = channel();
        let mut actual_endpoint = endpoint.clone();
        let before = accepted(&actual, &actual_heat);
        let mut seen = 0;
        let refused_heat = vec![0.96; cells];
        let refusal = actual.step_with_thermal_poll(
            &mut actual_heat, Some(&mut actual_endpoint), None,
            if cancel { &source } else { &refused_heat }, None, None, [0.0; 3],
            &mut || {
                seen += 1;
                if cancel && seen == polls { Err(GridError3::Cancelled) } else { Ok(()) }
            },
        );
        if cancel {
            assert!(matches!(refusal, Err(OpenThermalError3::Flow(GridError3::Cancelled))));
        } else {
            assert!(matches!(refusal, Err(OpenThermalError3::Thermal(ThermalOpenError3::Temperature { .. }))));
        }
        assert_eq!(accepted(&actual, &actual_heat), before);
        assert_eq!(actual_endpoint, endpoint);
        let report = actual.step_with_thermal_poll(
            &mut actual_heat, Some(&mut actual_endpoint), None, &source, None, None,
            [0.0; 3], &mut || Ok(()),
        ).unwrap();
        assert_eq!(accepted(&actual, &actual_heat), accepted(&expected, &expected_heat));
        assert_eq!(actual_endpoint, expected_endpoint);
        near(report.heat.combined_defect, 0.0, 1e-12);
        let mut external = report.heat.source;
        let mut storage = report.heat.stored_change;
        for reservoir in report.heat.reservoirs { external += reservoir.heat_in; }
        for wall in report.heat.walls.into_iter().flatten() {
            let solid = wall.solid.unwrap();
            external += solid.source;
            storage += solid.stored_change;
        }
        near(storage, external, 1e-12);
    }
}

#[test]
fn proposed_reservoirs_publish_only_with_finite_walls_and_closing_heat() {
    let (mut expected, mut expected_heat) = channel();
    let cells = expected.populations().len();
    let mut faces = crate::d3q19::Face3::ALL.map(|f| expected.boundaries().face(f));
    faces[0] = FaceBoundary3::Velocity { velocity: [0.035, 0.0, 0.0] };
    faces[1] = FaceBoundary3::Pressure { density: RHO + 0.001 };
    let targets = BoundarySpec3::new(faces);
    let endpoint = vec![[0.0002, 0.0, -0.0001]; cells];
    let source = vec![0.001; cells];
    let mut expected_endpoint = endpoint.clone();
    let plan = expected.inner.prepare_open_boundaries(targets, &mut || Ok(())).unwrap();
    let mut polls = 0;
    expected.step_with_thermal_plan_poll(&mut expected_heat, plan, Some(&mut expected_endpoint),
        None, &source, None, None, [0.0; 3], &mut || { polls += 1; Ok(()) }).unwrap();
    assert_eq!(expected.boundaries(), targets);
    for cancel in [false, true] {
        let (mut actual, mut heat) = channel();
        let old_targets = actual.boundaries(); let before = accepted(&actual, &heat);
        let plan = actual.inner.prepare_open_boundaries(targets, &mut || Ok(())).unwrap();
        let mut supplied_endpoint = endpoint.clone(); let mut seen = 0;
        let refused_heat = vec![0.96; cells];
        let result = actual.step_with_thermal_plan_poll(&mut heat, plan, Some(&mut supplied_endpoint),
            None, if cancel { &source } else { &refused_heat }, None, None, [0.0; 3], &mut || {
                seen += 1;
                if cancel && seen == polls { Err(GridError3::Cancelled) } else { Ok(()) }
            });
        if cancel { assert!(matches!(result, Err(OpenThermalError3::Flow(GridError3::Cancelled)))); }
        else { assert!(matches!(result, Err(OpenThermalError3::Thermal(ThermalOpenError3::Temperature { .. })))); }
        // Real momentum was staged before the closing heat/cancellation gate.
        assert_ne!(actual.inner.grid.next, actual.inner.grid.f);
        assert_eq!(accepted(&actual, &heat), before); assert_eq!(actual.boundaries(), old_targets);
        assert_eq!(supplied_endpoint, endpoint);
        actual.step_with_thermal_plan_poll(&mut heat, plan, Some(&mut supplied_endpoint), None,
            &source, None, None, [0.0; 3], &mut || Ok(())).unwrap();
        assert_eq!(accepted(&actual, &heat), accepted(&expected, &expected_heat));
        assert_eq!(actual.boundaries(), targets); assert_eq!(supplied_endpoint, expected_endpoint);
    }

    let (mut actual, mut heat) = channel();
    let before = accepted(&actual, &heat); let old_targets = actual.boundaries();
    let mut faces = crate::d3q19::Face3::ALL.map(|f| old_targets.face(f));
    faces[2] = FaceBoundary3::Wall { velocity: [0.0, 0.0, 0.001] };
    faces[3] = faces[2];
    let moving = actual.inner.prepare_open_boundaries(BoundarySpec3::new(faces), &mut || Ok(())).unwrap();
    assert!(matches!(actual.step_with_thermal_plan_poll(&mut heat, moving, None, None, &source,
        None, None, [0.0; 3], &mut || Ok(())),
        Err(OpenThermalError3::Thermal(ThermalOpenError3::Invalid(
            "moving finite-wall storage requires solid heat advection")))));
    assert_eq!(accepted(&actual, &heat), before); assert_eq!(actual.boundaries(), old_targets);
}

#[test]
fn g3_held_and_unchanged_endpoint_forces_agree_with_uniform_or_local_tensor_diffusion() {
    let (mut held, mut held_heat) = channel();
    let (mut advanced, mut advanced_heat) = channel();
    let cells = held.populations().len();
    let tensor = DiffusivityTensor::from_lattice([
        [0.08, 0.015, 0.0], [0.015, 0.1, -0.01], [0.0, -0.01, 0.09],
    ]).unwrap();
    let tensors = vec![tensor; cells];
    let source: Vec<_> = (0..cells).map(|cell| 0.0001 * (1 + cell % 5) as f64).collect();
    let mut endpoint = advanced.forces().to_vec();
    for _ in 0..8 {
        let a = held.step_with_thermal_poll(
            &mut held_heat, None, None, &source, None, Some(DiffusivityField::Uniform(&tensor)),
            [0.0; 3], &mut || Ok(()),
        ).unwrap();
        let b = advanced.step_with_thermal_poll(
            &mut advanced_heat, Some(&mut endpoint), None, &source, None,
            Some(DiffusivityField::PerCell(&tensors)), [0.0; 3], &mut || Ok(()),
        ).unwrap();
        assert_eq!(accepted(&held, &held_heat), accepted(&advanced, &advanced_heat));
        assert_eq!(a.heat, b.heat);
        near(a.heat.combined_defect, 0.0, 1e-12);
    }
}

#[test]
fn g3_moving_wall_receipt_publishes_with_heat_and_survives_late_cancellation() {
    for change_force in [false, true] {
        let dims = [4, 4, 2];
        let cells = dims.into_iter().product();
        let mut flow = flow::<D3q19, 19>(
            dims, CollisionModel3::Bgk { tau: 0.9 },
            BoundarySpec3::new([
                FaceBoundary3::Periodic, FaceBoundary3::Periodic,
                FaceBoundary3::stationary_wall(),
                FaceBoundary3::Wall { velocity: [0.03, 0.0, 0.0] },
                FaceBoundary3::Periodic, FaceBoundary3::Periodic,
            ]),
            vec![[1e-5, 0.0, 0.0]; cells],
            |_| (RHO, [0.005, 0.0, 0.0]),
        );
        let mut boundaries = periodic_heat();
        boundaries.faces[2] = ThermalOpenFace3::Wall(ThermalWall::Adiabatic);
        boundaries.faces[3] = ThermalOpenFace3::Wall(ThermalWall::Adiabatic);
        let mut heat = flow.prepare_thermal_field_with_poll(
            &vec![0.5; cells], 0.83, boundaries, limits(), &mut || Ok(()),
        ).unwrap();
        let source = vec![0.001; cells];
        let mut endpoint = vec![[2e-5, 0.0, 0.0]; cells];
        let mut polls = 0;
        flow.step_with_thermal_poll(
            &mut heat,
            if change_force { Some(&mut endpoint) } else { None },
            None, &source, None, None, [0.0; 3],
            &mut || { polls += 1; Ok(()) },
        ).unwrap();
        let wall = flow.last_wall_exchange().expect("accepted moving-wall receipt");
        assert!(wall.measured_links > 0);
        assert!(wall.wall_work < 0.0, "the faster lid supplies mechanical work");
        // Mechanical wall work has its own receipt, not an extra thermal credit.
        near(heat.last_balance().unwrap().stored_change, cells as f64 * 0.001, 1e-12);
        let before = accepted(&flow, &heat);
        endpoint.fill([3e-5, 0.0, 0.0]);
        let endpoint_before = endpoint.clone();
        let mut seen = 0;
        let refused = flow.step_with_thermal_poll(
            &mut heat,
            if change_force { Some(&mut endpoint) } else { None },
            None, &source, None, None, [0.0; 3],
            &mut || {
                seen += 1;
                if seen == polls { Err(GridError3::Cancelled) } else { Ok(()) }
            },
        );
        assert!(matches!(refused, Err(OpenThermalError3::Flow(GridError3::Cancelled))));
        assert_eq!(accepted(&flow, &heat), before);
        assert_eq!(endpoint, endpoint_before);
        flow.step_with_thermal_poll(
            &mut heat,
            if change_force { Some(&mut endpoint) } else { None },
            None, &source, None, None, [0.0; 3], &mut || Ok(()),
        ).unwrap();
        let after = flow.last_wall_exchange().unwrap();
        assert_eq!(after.measured_links, wall.measured_links);
        assert_ne!(after.wall_work.to_bits(), wall.wall_work.to_bits());
        assert_eq!(flow.completed_steps(), 2);
        assert_eq!(heat.completed_steps(), 2);
    }
}
