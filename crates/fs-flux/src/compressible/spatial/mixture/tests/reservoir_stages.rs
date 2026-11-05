use super::*;
use crate::compressible::AxisBoundaries;
use crate::compressible::mixture::tests::{close, standard};

#[path = "forced_reservoir_stages.rs"]
mod forced;

const OPEN: AxisBoundaries = AxisBoundaries::Ends {
    lower: BoundaryCondition::Transmissive,
    upper: BoundaryCondition::Transmissive,
};

fn open_axis(axis: usize) -> Boundaries3d {
    let mut axes = [AxisBoundaries::Periodic; 3];
    axes[axis] = OPEN;
    Boundaries3d {
        x: axes[0],
        y: axes[1],
        z: axes[2],
    }
}

fn input(densities: &mut [f64], velocity: [f64; 3], pressure_pa: f64) -> ReservoirFaceState {
    densities.copy_from_slice(&[0.6, 0.4, 0.0]);
    ReservoirFaceState {
        velocity,
        pressure_pa,
    }
}

fn prepared(
    grid: CartesianGrid,
    topology: Boundaries3d,
    sides: &[BoundarySide],
    velocity: [f64; 3],
) -> MixtureEuler3d {
    let gas = standard(0.2);
    let mut states = vec![0.0; grid.len() * 7];
    for row in states.chunks_exact_mut(7) {
        state_into(&gas, &[0.6, 0.4, 0.0], velocity, 1.0, row).unwrap();
    }
    let bytes = MixtureEuler3d::required_workspace_bytes(&gas, &grid).unwrap();
    let mut flow = MixtureEuler3d::new(gas, grid, states, topology, bytes, || false).unwrap();
    let bytes = flow
        .required_reservoir_profile_workspace_bytes(sides)
        .unwrap();
    flow.prepare_reservoir_profiles(
        sides,
        bytes,
        |_, _, _, rho| Ok(input(rho, velocity, 1.0)),
        || false,
    )
    .unwrap();
    flow
}

#[derive(Debug, PartialEq)]
struct Snapshot {
    states: Vec<f64>,
    profiles: Vec<f64>,
    loads: Vec<f64>,
    time: f64,
    report: Option<StepReport>,
    bytes: usize,
}

fn snapshot(flow: &MixtureEuler3d) -> Snapshot {
    let mut profiles = Vec::new();
    let mut loads = Vec::new();
    for side in BoundarySide::ALL {
        for face in 0..flow.grid().face_count(side) {
            if let Some(p) = flow.reservoir_face(side, face) {
                profiles.extend_from_slice(p.partial_densities);
                profiles.extend_from_slice(&p.velocity);
                profiles.push(p.pressure_pa);
            }
        }
        for load in flow.boundary_impulses(side) {
            loads.extend_from_slice(load.species_masses);
            loads.extend_from_slice(&load.momentum);
            loads.push(load.energy);
        }
    }
    Snapshot {
        states: flow.states().to_vec(),
        profiles,
        loads,
        time: flow.time().value(),
        report: flow.last_step().copied(),
        bytes: flow.workspace_bytes(),
    }
}

fn balance(report: &StepReport) {
    for value in report.conservation_residual.species_masses {
        close(value, 0.0, 3e-12);
    }
    for value in report.conservation_residual.momentum {
        close(value, 0.0, 3e-12);
    }
    close(report.conservation_residual.energy, 0.0, 3e-12);
    close(report.absolute_energy_residual_j, 0.0, 3e-10);
}

#[test]
fn g0_constant_stage_profiles_retain_native_transport_and_absolute_face_order() {
    let grid = CartesianGrid::new(2, 3, 2, 0.2, 0.3, 0.4).unwrap();
    let topology = Boundaries3d {
        x: OPEN,
        y: OPEN,
        z: OPEN,
    };
    let velocity = [0.15, -0.1, 0.08];
    let order: Vec<_> = BoundarySide::ALL
        .into_iter()
        .flat_map(|side| (0..grid.face_count(side)).map(move |face| (side, face)))
        .collect();
    for reconstruction in [Reconstruction::FirstOrder, Reconstruction::Muscl] {
        let controls = StepControl {
            reconstruction,
            ..StepControl::default()
        };
        let mut held = prepared(grid, topology, &BoundarySide::ALL, velocity);
        let mut timed = prepared(grid, topology, &BoundarySide::ALL, velocity);
        held.step(Time::new(0.001), controls, || false).unwrap();
        timed.step(Time::new(0.001), controls, || false).unwrap();
        let start = timed.time().value();
        let bytes = timed.workspace_bytes();
        let mut calls = Vec::new();
        let report = timed
            .step_with_reservoir_profiles(
                Time::new(0.002),
                controls,
                2 * order.len(),
                |time, side, face, position, rho| {
                    assert_eq!(position, grid.face_center(side, face).unwrap());
                    assert!(rho.iter().all(|v| v.is_nan()));
                    calls.push((time, side, face));
                    Ok(input(rho, velocity, 1.0))
                },
                || false,
            )
            .unwrap();
        let native = held.step(Time::new(0.002), controls, || false).unwrap();
        assert_eq!(calls.len(), 2 * order.len());
        for (stage, entries) in calls.chunks_exact(order.len()).enumerate() {
            for (&(time, side, face), &(expected_side, expected_face)) in entries.iter().zip(&order)
            {
                assert_eq!((side, face), (expected_side, expected_face));
                assert_eq!(
                    time,
                    start
                        + if stage == 0 {
                            0.0
                        } else {
                            report.elapsed.value()
                        }
                );
            }
        }
        assert_eq!(timed.states(), held.states());
        assert_eq!(snapshot(&timed).loads, snapshot(&held).loads);
        assert_eq!(timed.time(), held.time());
        assert_eq!(report.before, native.before);
        assert_eq!(report.after, native.after);
        assert_eq!(report.boundary_change, native.boundary_change);
        assert_eq!(timed.workspace_bytes(), bytes);
        balance(&report);
    }
}

fn varying(time: f64, side: BoundarySide, face: usize, rho: &mut [f64]) -> ReservoirFaceState {
    let offset = 0.01 * face as f64;
    rho.copy_from_slice(&[0.6 + offset + 0.8 * time, 0.4 + offset - 0.2 * time, 0.0]);
    let axis = side.axis();
    let [a, b] = tangents(axis);
    let mut velocity = [0.0; 3];
    velocity[axis] = -side.sign() * (4.0 + 3.0 * time);
    velocity[a] = 0.15 + 0.4 * time;
    velocity[b] = -0.07 + 0.2 * time;
    ReservoirFaceState {
        velocity,
        pressure_pa: 1.0 + offset + 0.5 * time,
    }
}

#[test]
fn g1_time_varying_species_shear_pressure_and_enthalpy_use_both_original_stages() {
    let grid = CartesianGrid::new(2, 2, 2, 0.3, 0.4, 0.5).unwrap();
    for side in BoundarySide::ALL {
        for reconstruction in [Reconstruction::FirstOrder, Reconstruction::Muscl] {
            let axis = side.axis();
            let mut initial_velocity = [0.0; 3];
            initial_velocity[axis] = -side.sign() * 4.0;
            let mut flow = prepared(grid, open_axis(axis), &[side], initial_velocity);
            let report = flow
                .step_with_reservoir_profiles(
                    Time::new(1e-4),
                    StepControl {
                        reconstruction,
                        ..StepControl::default()
                    },
                    2 * grid.face_count(side),
                    |time, side, face, _, rho| Ok(varying(time, side, face, rho)),
                    || false,
                )
                .unwrap();
            let dt = report.elapsed.value();
            for (face, load) in flow.boundary_impulses(side).enumerate() {
                let mut integrated = [0.0; 7];
                for time in [0.0, dt] {
                    let mut rho = [0.0; 3];
                    let input = varying(time, side, face, &mut rho);
                    let mut q = [0.0; 7];
                    state_into(flow.gas(), &rho, input.velocity, input.pressure_pa, &mut q)
                        .unwrap();
                    let outward = side.sign() * input.velocity[axis];
                    let scale = 0.5 * dt * grid.face_area(side);
                    for d in 0..3 {
                        integrated[d] += scale * outward * rho[d];
                        integrated[d + 3] += scale
                            * (outward * q[d + 3]
                                + if d == axis {
                                    side.sign() * input.pressure_pa
                                } else {
                                    0.0
                                });
                    }
                    integrated[6] += scale * outward * (q[6] + input.pressure_pa);
                }
                for d in 0..3 {
                    close(load.species_masses[d], integrated[d], 3e-12);
                    close(load.momentum[d], integrated[d + 3], 3e-12);
                }
                close(load.energy, integrated[6], 3e-11);
                assert_eq!(load.species_masses[2], 0.0);
                let endpoint = flow.reservoir_face(side, face).unwrap();
                let p = flow.reservoirs.primitive_face(side, face).unwrap();
                let r = endpoint.partial_densities[0] + 0.5 * endpoint.partial_densities[1];
                let temperature = endpoint.pressure_pa / r;
                let cp = 3.5 * endpoint.partial_densities[0]
                    + 0.5 * endpoint.partial_densities[1] * (3.5 + 0.2 * temperature);
                close(p.temperature, temperature, 3e-13);
                close(
                    p.sound_speed.powi(2),
                    cp / (cp - r) * endpoint.pressure_pa / p.density,
                    3e-13,
                );
                let mut rho = [0.0; 3];
                let expected = varying(dt, side, face, &mut rho);
                assert_eq!(endpoint.partial_densities, rho);
                assert_eq!(endpoint.velocity, expected.velocity);
                assert_eq!(endpoint.pressure_pa, expected.pressure_pa);
            }
            balance(&report);
        }
    }
}

#[test]
fn g1_smooth_complete_inflow_refines_second_order_in_physical_time() {
    let grid = CartesianGrid::new(3, 1, 1, 0.5, 0.4, 0.3).unwrap();
    let side = BoundarySide::XMin;
    let end = 0.04_f64;
    // Supersonic incoming species flux is exactly rho_A(t) * u(t). Its
    // independently integrated quadratic isolates boundary temporal error from
    // interior spatial error and source inversion tolerances.
    let exact = grid.face_area(side)
        * (0.6 * 4.0 * end
            + (0.6 * 3.0 + 0.8 * 4.0) * end.powi(2) / 2.0
            + 0.8 * 3.0 * end.powi(3) / 3.0);
    let mut errors = Vec::new();
    for steps in [10, 20, 40] {
        let mut flow = prepared(grid, open_axis(0), &[side], [4.0, 0.0, 0.0]);
        let mut mass = 0.0;
        for step in 0..steps {
            let duration = if step + 1 == steps {
                end - flow.time().value()
            } else {
                end / steps as f64
            };
            let report = flow
                .step_with_reservoir_profiles(
                    Time::new(duration),
                    StepControl::default(),
                    2,
                    |time, side, face, _, rho| Ok(varying(time, side, face, rho)),
                    || false,
                )
                .unwrap();
            assert_eq!(report.halvings, 0);
            mass -= flow.boundary_impulses(side).next().unwrap().species_masses[0];
            balance(&report);
        }
        close(flow.time().value(), end, 2e-14);
        errors.push((mass - exact).abs());
    }
    assert!(errors[0] > 1e-10);
    for pair in errors.windows(2) {
        let ratio = pair[0] / pair[1];
        assert!((3.8..4.2).contains(&ratio), "temporal errors: {errors:?}");
    }
}

fn retry_fixture() -> MixtureEuler3d {
    let grid = CartesianGrid::new(2, 1, 1, 1.0, 1.0, 1.0).unwrap();
    let mut flow = prepared(grid, open_axis(0), &[BoundarySide::XMin], [4.0, 0.0, 0.0]);
    flow.step(Time::new(1e-4), StepControl::default(), || false)
        .unwrap();
    flow
}

#[test]
fn g4_endpoint_cfl_resamples_shortened_time_and_all_source_work_is_bounded() {
    let mut flow = retry_fixture();
    let start = flow.time().value();
    let mut calls = Vec::new();
    let report = flow
        .step_with_reservoir_profiles(
            Time::new(0.02),
            StepControl::default(),
            4,
            |time, _, _, _, rho| {
                calls.push(time);
                Ok(input(rho, [4.0 + 1000.0 * (time - start), 0.0, 0.0], 1.0))
            },
            || false,
        )
        .unwrap();
    assert_eq!(report.halvings, 1);
    assert!(!report.first_order_retry);
    assert_eq!(calls, [start, start + 0.02, start, start + 0.01]);
    assert_eq!(report.elapsed.value(), 0.01);
    close(
        flow.reservoir_face(BoundarySide::XMin, 0).unwrap().velocity[0],
        14.0,
        1e-14,
    );
    assert!(report.maximum_cfl <= StepControl::default().cfl * (1.0 + 2e-14));
    balance(&report);

    // A large prescribed tangential speed bounds the unsplit interval without
    // adding a stronger normal shock. The previous held trace has no shear.
    let mut hot = retry_fixture();
    let controls = StepControl {
        cfl: 0.1,
        maximum_halvings: 0,
        ..StepControl::default()
    };
    let report = hot
        .step_with_reservoir_profiles(
            Time::new(1.0),
            controls,
            2,
            |_, _, _, _, rho| Ok(input(rho, [4.0, 20.0, 0.0], 1.0)),
            || false,
        )
        .unwrap();
    let p = hot
        .reservoirs
        .primitive_face(BoundarySide::XMin, 0)
        .unwrap();
    close(
        report.elapsed.value(),
        controls.cfl / (24.0 + 3.0 * p.sound_speed),
        2e-14,
    );
    assert_eq!(report.halvings, 0);

    for (budget, refusal) in [(0, 0), (3, 0), (20, 1), (20, 2), (20, 3)] {
        let mut rejected = retry_fixture();
        let before = snapshot(&rejected);
        let start = rejected.time().value();
        let mut calls = 0;
        let error = rejected
            .step_with_reservoir_profiles(
                Time::new(0.02),
                StepControl::default(),
                budget,
                |time, _, _, _, rho| {
                    calls += 1;
                    if refusal == 1 && time > start {
                        return Err(Error::Invalid("prescribed source unavailable"));
                    }
                    let mut state = input(rho, [4.0 + 1000.0 * (time - start), 0.0, 0.0], 1.0);
                    if refusal == 2 && time > start {
                        state.pressure_pa = 1e5;
                    }
                    if refusal == 3 {
                        rho[2] = f64::NAN;
                    }
                    Ok(state)
                },
                || false,
            )
            .unwrap_err();
        assert!(!matches!(error, Error::StepRetries { .. }));
        assert_eq!(
            calls,
            if refusal == 3 {
                1
            } else if refusal != 0 {
                2
            } else {
                budget
            }
        );
        assert_eq!(snapshot(&rejected), before);
    }
}

#[test]
fn g4_last_poll_cancellation_preserves_profiles_cells_clock_loads_and_report() {
    let mut reference = retry_fixture();
    let mut polls = 0;
    reference
        .step_with_reservoir_profiles(
            Time::new(0.001),
            StepControl::default(),
            2,
            |time, side, face, _, rho| Ok(varying(time, side, face, rho)),
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    let mut flow = retry_fixture();
    let before = snapshot(&flow);
    let mut seen = 0;
    let error = flow
        .step_with_reservoir_profiles(
            Time::new(0.001),
            StepControl::default(),
            2,
            |time, side, face, _, rho| Ok(varying(time, side, face, rho)),
            || {
                seen += 1;
                seen == polls
            },
        )
        .unwrap_err();
    assert!(matches!(error, Error::Cancelled));
    assert_eq!(snapshot(&flow), before);
    flow.step_with_reservoir_profiles(
        Time::new(0.001),
        StepControl::default(),
        2,
        |time, side, face, _, rho| Ok(varying(time, side, face, rho)),
        || false,
    )
    .unwrap();
    assert_eq!(snapshot(&flow), snapshot(&reference));
}
