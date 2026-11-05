use super::*;
use crate::compressible::tests::{close, gas};

const OPEN: AxisBoundaries = AxisBoundaries::Ends {
    lower: BoundaryCondition::Transmissive,
    upper: BoundaryCondition::Transmissive,
};
const ALL_OPEN: Boundaries3d = Boundaries3d {
    x: OPEN,
    y: OPEN,
    z: OPEN,
};

fn field(grid: CartesianGrid, boundaries: Boundaries3d, velocity: [f64; 3]) -> Euler3d {
    let source = gas(0.4);
    let q = State::from_primitive(&source, 1.0, velocity, 1.0).unwrap();
    Euler3d::new(source, grid, vec![q; grid.len()], boundaries, || false).unwrap()
}

fn face_value(face: usize) -> ReservoirFaceState {
    ReservoirFaceState {
        density: 1.0 + 0.02 * face as f64,
        velocity: [0.03 * face as f64, -0.02 * face as f64, 0.01],
        pressure_pa: 1.2 + 0.04 * face as f64,
    }
}

fn snapshot(flow: &Euler3d) -> Vec<Option<ReservoirFaceState>> {
    BoundarySide::ALL
        .into_iter()
        .flat_map(|side| {
            (0..flow.grid.face_count(side)).map(move |face| flow.reservoir_face(side, face))
        })
        .collect()
}

#[test]
fn g0_uniform_profile_retains_native_arithmetic_and_canonical_face_order() {
    let grid = CartesianGrid::new(3, 2, 2, 0.2, 0.3, 0.4).unwrap();
    let velocity = [0.3, -0.2, 0.1];
    for reconstruction in [Reconstruction::FirstOrder, Reconstruction::Muscl] {
        let source = gas(0.4);
        let q = State::from_primitive(&source, 1.0, velocity, 1.0).unwrap();
        let inputs = BoundarySide::ALL.map(|side| ReservoirBoundary {
            side,
            density: 1.0,
            velocity,
            pressure_pa: 1.0,
        });
        let mut uniform = Euler3d::new_with_reservoirs(
            source,
            grid,
            vec![q; grid.len()],
            ALL_OPEN,
            &inputs,
            || false,
        )
        .unwrap();
        let mut profiled = field(grid, ALL_OPEN, velocity);
        let required = profiled
            .required_reservoir_profile_workspace_bytes(&BoundarySide::ALL)
            .unwrap();
        let mut visited = Vec::new();
        let mut reversed = BoundarySide::ALL;
        reversed.reverse();
        profiled
            .prepare_reservoir_profiles(
                &reversed,
                required,
                |side, face, position| {
                    assert_eq!(position, grid.face_center(side, face).unwrap());
                    visited.push((side, face));
                    Ok(ReservoirFaceState {
                        density: 1.0,
                        velocity,
                        pressure_pa: 1.0,
                    })
                },
                || false,
            )
            .unwrap();
        let expected: Vec<_> = BoundarySide::ALL
            .into_iter()
            .flat_map(|side| (0..grid.face_count(side)).map(move |face| (side, face)))
            .collect();
        assert_eq!(visited, expected);
        assert_eq!(profiled.prepared_reservoir_profile_sides(), [true; 6]);
        assert_eq!(profiled.reservoir_profile_workspace_bytes(), required);
        assert!(profiled.reservoir(BoundarySide::XMin).is_none());
        for side in BoundarySide::ALL {
            assert_eq!(profiled.reservoir_face(side, grid.face_count(side)), None);
        }
        for _ in 0..3 {
            let controls = StepControl {
                reconstruction,
                ..StepControl::default()
            };
            let expected = uniform.step(Time::new(0.002), controls, || false).unwrap();
            let actual = profiled.step(Time::new(0.002), controls, || false).unwrap();
            assert_eq!(actual, expected);
            assert_eq!(profiled.states(), uniform.states());
        }
    }
}

#[test]
fn g1_each_face_uses_its_own_physical_mass_momentum_and_enthalpy_flux() {
    let grid = CartesianGrid::new(3, 2, 2, 0.2, 0.3, 0.4).unwrap();
    for side in BoundarySide::ALL {
        let axis = side.axis();
        let [a, b] = tangents(axis);
        let mut velocity = [0.0; 3];
        velocity[axis] = -4.0 * side.sign();
        let mut pairs = [AxisBoundaries::Periodic; 3];
        pairs[axis] = OPEN;
        let topology = Boundaries3d {
            x: pairs[0],
            y: pairs[1],
            z: pairs[2],
        };
        for reconstruction in [Reconstruction::FirstOrder, Reconstruction::Muscl] {
            let mut flow = field(grid, topology, velocity);
            let required = flow
                .required_reservoir_profile_workspace_bytes(&[side])
                .unwrap();
            flow.prepare_reservoir_profiles(
                &[side],
                required,
                |_, face, point| {
                    let mut u = velocity;
                    u[a] = 0.07 * point[a];
                    u[b] = -0.05 * point[b];
                    Ok(ReservoirFaceState {
                        density: 1.0 + 0.04 * face as f64,
                        velocity: u,
                        pressure_pa: 1.0,
                    })
                },
                || false,
            )
            .unwrap();
            let report = flow
                .step(
                    Time::new(1e-5),
                    StepControl {
                        reconstruction,
                        ..StepControl::default()
                    },
                    || false,
                )
                .unwrap();
            // Both Riemann characteristics point into the box. The upstream
            // physical Euler flux is independent of the evolving interior,
            // giving a closed-form check of both SSP stages on all six sides.
            for (face, impulse) in flow.boundary_impulses(side).iter().enumerate() {
                let p = flow.reservoir_face(side, face).unwrap();
                let q = State::from_primitive(flow.gas(), p.density, p.velocity, p.pressure_pa)
                    .unwrap();
                let mass_flux = p.density * p.velocity[axis];
                let scale = side.sign() * report.elapsed.value() * grid.face_area(side);
                close(impulse.mass, scale * mass_flux, 3e-13);
                for d in 0..3 {
                    let pressure = if d == axis { p.pressure_pa } else { 0.0 };
                    close(
                        impulse.momentum[d],
                        scale * (mass_flux * p.velocity[d] + pressure),
                        3e-13,
                    );
                }
                close(
                    impulse.energy,
                    scale * (q.energy + p.pressure_pa) * p.velocity[axis],
                    3e-13,
                );
            }
            for residual in report.conservation_residual.array() {
                close(residual, 0.0, 6e-12);
            }
        }
    }
}

#[test]
fn g4_profile_budget_source_and_late_cancellation_preserve_every_active_side() {
    let grid = CartesianGrid::new(2, 2, 2, 0.2, 0.3, 0.4).unwrap();
    let sides = [BoundarySide::ZMax, BoundarySide::XMin];
    let mut flow = field(grid, ALL_OPEN, [0.0; 3]);
    let required = flow
        .required_reservoir_profile_workspace_bytes(&sides)
        .unwrap();
    let mut calls = 0;
    assert!(
        flow.prepare_reservoir_profiles(
            &sides,
            required - 1,
            |_, face, _| {
                calls += 1;
                Ok(face_value(face))
            },
            || false
        )
        .is_err()
    );
    assert_eq!(calls, 0);
    assert_eq!(flow.prepared_reservoir_profile_sides(), [false; 6]);
    assert!(
        flow.prepare_reservoir_profiles(
            &[sides[0], sides[0]],
            required,
            |_, _, _| { unreachable!("duplicate geometry is rejected before sampling") },
            || false
        )
        .is_err()
    );
    flow.prepare_reservoir_profiles(
        &sides,
        required,
        |_, face, _| Ok(face_value(face)),
        || false,
    )
    .unwrap();
    flow.step(Time::new(1e-4), StepControl::default(), || false)
        .unwrap();
    let faces = snapshot(&flow);
    let cells = flow.states().to_vec();
    let time = flow.time();
    let report = flow.last_step;
    let impulses = flow.boundary_impulses.clone();
    assert!(
        flow.update_reservoir_profiles(
            |side, face, _| {
                let mut p = face_value(face);
                p.pressure_pa = if side == BoundarySide::ZMax && face == 3 {
                    f64::NAN
                } else {
                    2.0
                };
                Ok(p)
            },
            || false
        )
        .is_err()
    );
    assert_eq!(snapshot(&flow), faces);
    let stop = std::cell::Cell::new(false);
    assert_eq!(
        flow.update_reservoir_profiles(
            |side, face, _| {
                stop.set(side == BoundarySide::ZMax && face == 3);
                Ok(face_value(face + 1))
            },
            || stop.get()
        ),
        Err(Error::Cancelled)
    );
    assert_eq!(snapshot(&flow), faces);
    assert_eq!(flow.states(), cells);
    assert_eq!(flow.time(), time);
    assert_eq!(flow.last_step, report);
    assert_eq!(flow.boundary_impulses, impulses);
    flow.update_reservoir_profiles(|_, face, _| Ok(face_value(face + 1)), || false)
        .unwrap();
    assert_eq!(flow.reservoir_profile_workspace_bytes(), required);
    assert_eq!(
        flow.reservoir_face(BoundarySide::XMin, 0),
        Some(face_value(1))
    );
    assert_eq!(
        flow.reservoir_face(BoundarySide::ZMax, 3),
        Some(face_value(4))
    );
    assert!(flow.last_step.is_none());

    let mut periodic = field(
        grid,
        Boundaries3d {
            x: AxisBoundaries::Periodic,
            ..ALL_OPEN
        },
        [0.0; 3],
    );
    assert!(
        periodic
            .prepare_reservoir_profiles(
                &[BoundarySide::XMin],
                required,
                |_, _, _| unreachable!(),
                || false
            )
            .is_err()
    );
    let source = gas(0.4);
    let q = State::from_primitive(&source, 1.0, [0.0; 3], 1.0).unwrap();
    let mut bound = Euler3d::new_with_open_boundaries(
        source,
        grid,
        vec![q; grid.len()],
        ALL_OPEN,
        &[ReservoirBoundary {
            side: BoundarySide::XMin,
            density: 1.0,
            velocity: [0.0; 3],
            pressure_pa: 1.0,
        }],
        &[PressureOutletBoundary {
            side: BoundarySide::ZMax,
            pressure_pa: 1.0,
        }],
        || false,
    )
    .unwrap();
    for side in sides {
        assert!(
            bound
                .prepare_reservoir_profiles(&[side], required, |_, _, _| unreachable!(), || false)
                .is_err()
        );
    }
}

#[test]
fn g1_hottest_profile_face_bounds_the_whole_unsplit_interval() {
    let source = gas(0.0);
    let grid = CartesianGrid::new(2, 2, 2, 0.1, 0.1, 0.1).unwrap();
    let inside = State::from_primitive(&source, 1.0, [0.0; 3], 0.01).unwrap();
    let walls = AxisBoundaries::Ends {
        lower: BoundaryCondition::Reflecting,
        upper: BoundaryCondition::Reflecting,
    };
    let mut flow = Euler3d::new(
        source,
        grid,
        vec![inside; grid.len()],
        Boundaries3d {
            x: OPEN,
            y: walls,
            z: walls,
        },
        || false,
    )
    .unwrap();
    flow.prepare_reservoir_profiles(
        &[BoundarySide::XMin],
        usize::MAX,
        |_, face, _| {
            Ok(ReservoirFaceState {
                density: 1.0,
                velocity: [0.0; 3],
                pressure_pa: if face == 3 { 32.0 } else { 0.01 },
            })
        },
        || false,
    )
    .unwrap();
    let control = StepControl {
        cfl: 0.2,
        maximum_halvings: 0,
        reconstruction: Reconstruction::FirstOrder,
    };
    let report = flow.step(Time::new(1.0), control, || false).unwrap();
    let bound = control.cfl * 0.1 / (3.0 * (1.4_f64 * 32.0).sqrt());
    assert!(report.elapsed.value() <= bound * (1.0 + 2e-14));
    assert_eq!(report.halvings, 0);
    assert!(report.boundary_change.mass > 0.0);
    assert!(report.maximum_cfl <= control.cfl * (1.0 + 2e-14));
}
