use super::*;
use crate::compressible::tests::{close, gas};

fn held(side: BoundarySide, pressure: f64, velocity: [f64; 3]) -> Euler3d {
    let source = gas(0.4);
    let grid = CartesianGrid::new(2, 2, 2, 0.2, 0.2, 0.2).unwrap();
    let initial = State::from_primitive(&source, 1.0, velocity, 1.0).unwrap();
    let mut pairs = [AxisBoundaries::Periodic; 3];
    pairs[side.axis()] = AxisBoundaries::Ends {
        lower: if side.sign() < 0.0 {
            BoundaryCondition::Transmissive
        } else {
            BoundaryCondition::Reflecting
        },
        upper: if side.sign() > 0.0 {
            BoundaryCondition::Transmissive
        } else {
            BoundaryCondition::Reflecting
        },
    };
    Euler3d::new_with_reservoirs(
        source,
        grid,
        vec![initial; grid.len()],
        Boundaries3d {
            x: pairs[0],
            y: pairs[1],
            z: pairs[2],
        },
        &[ReservoirBoundary {
            side,
            density: 1.0,
            velocity: [0.0; 3],
            pressure_pa: pressure,
        }],
        || false,
    )
    .unwrap()
}

fn finite(side: BoundarySide, pressure: f64, volume_m3: f64, velocity: [f64; 3]) -> Euler3d {
    held(side, pressure, velocity)
        .with_finite_reservoirs(&[FiniteReservoirBoundary { side, volume_m3 }], || false)
        .unwrap()
}

fn storage(flow: &Euler3d, side: BoundarySide) -> [f64; 2] {
    let field = flow.totals();
    let tank = flow.finite_reservoir(side).unwrap();
    [
        field.mass + tank.state().density * tank.volume_m3(),
        field.energy + tank.state().energy * tank.volume_m3(),
    ]
}

#[test]
fn g0_finite_tanks_fill_drain_and_conserve_all_axes_and_reconstructions() {
    for side in BoundarySide::ALL {
        for reconstruction in [Reconstruction::FirstOrder, Reconstruction::Muscl] {
            for pressure in [0.5, 2.0] {
                let [a, b] = tangents(side.axis());
                let mut velocity = [0.0; 3];
                if pressure < 1.0 {
                    velocity[a] = 0.2;
                    velocity[b] = -0.3;
                }
                let mut flow = finite(side, pressure, 0.04, velocity);
                let original = storage(&flow, side);
                let tank0 = *flow.finite_reservoir(side).unwrap();
                let temperature0 = flow
                    .gas()
                    .primitive(tank0.state(), || false)
                    .unwrap()
                    .temperature;
                for _ in 0..3 {
                    let report = flow
                        .step(
                            Time::new(0.0002),
                            StepControl {
                                reconstruction,
                                ..StepControl::default()
                            },
                            || false,
                        )
                        .unwrap();
                    assert_eq!(report.halvings, 0);
                    for (got, want) in storage(&flow, side).into_iter().zip(original) {
                        close(got, want, 3e-12);
                    }
                    let tank = flow.finite_reservoir_report(side).unwrap();
                    close(tank.mass_residual_kg, 0.0, 2e-14);
                    close(tank.energy_residual_j, 0.0, 2e-14);
                    close(report.boundary_change.mass, -tank.transfer.mass, 2e-14);
                    close(report.boundary_change.energy, -tank.transfer.energy, 2e-14);
                    assert_eq!(tank.after, *flow.finite_reservoir(side).unwrap());
                    assert!(if pressure < 1.0 {
                        tank.transfer.mass > 0.0
                    } else {
                        tank.transfer.mass < 0.0
                    });
                    if pressure < 1.0 {
                        close(
                            tank.transfer.momentum[a],
                            velocity[a] * tank.transfer.mass,
                            3e-12,
                        );
                        close(
                            tank.transfer.momentum[b],
                            velocity[b] * tank.transfer.mass,
                            3e-12,
                        );
                    }
                    for (dimension, value) in tank.transfer.array().into_iter().enumerate() {
                        close(
                            value,
                            flow.boundary_impulses(side)
                                .iter()
                                .map(|q| q.array()[dimension])
                                .sum(),
                            2e-14,
                        );
                    }
                }
                let tank = flow.finite_reservoir(side).unwrap();
                let p = flow.gas().primitive(tank.state(), || false).unwrap();
                if pressure > 1.0 {
                    assert!(p.pressure < pressure);
                    assert!(p.temperature < temperature0);
                } else {
                    assert!(p.pressure > pressure);
                }
                close(flow.reservoir(side).unwrap().pressure_pa, p.pressure, 2e-14);
            }
        }
    }
}

#[test]
fn g0_finite_tank_equilibrium_and_empty_binding_preserve_native_bits() {
    let side = BoundarySide::YMax;
    let mut tank = finite(side, 1.0, 0.04, [0.0; 3]);
    let mut fixed = held(side, 1.0, [0.0; 3])
        .with_stagnation_reservoirs(&[side], || false)
        .unwrap();
    let original = *tank.finite_reservoir(side).unwrap();
    for _ in 0..4 {
        let a = tank
            .step(Time::new(0.001), StepControl::default(), || false)
            .unwrap();
        let b = fixed
            .step(Time::new(0.001), StepControl::default(), || false)
            .unwrap();
        assert_eq!(a, b);
        assert_eq!(tank.states(), fixed.states());
        assert_eq!(tank.finite_reservoir(side), Some(&original));
        assert_eq!(
            tank.finite_reservoir_report(side).unwrap().transfer.mass,
            0.0
        );
        assert_eq!(
            tank.finite_reservoir_report(side).unwrap().transfer.energy,
            0.0
        );
    }
    let mut ordinary = held(side, 2.0, [0.0; 3]);
    let mut empty = held(side, 2.0, [0.0; 3])
        .with_finite_reservoirs(&[], || false)
        .unwrap();
    assert_eq!(
        ordinary
            .step(Time::new(0.001), StepControl::default(), || false)
            .unwrap(),
        empty
            .step(Time::new(0.001), StepControl::default(), || false)
            .unwrap(),
    );
    assert_eq!(ordinary.states(), empty.states());
}

#[test]
fn g1_tiny_finite_tank_limits_the_initial_step_before_cell_cfl() {
    let side = BoundarySide::XMin;
    let volume = 1e-6;
    let mut flow = finite(side, 2.0, volume, [0.0; 3]);
    let tank = flow.finite_reservoir(side).unwrap().state();
    let flux = flow
        .gas()
        .reservoir_spatial_flux(flow.states()[0].normal(0), tank, -1.0, || false)
        .unwrap();
    let rate =
        0.5 * flow.grid().face_area(side) * flow.grid().face_count(side) as f64 * flux.signal_speed
            / volume;
    let controls = StepControl {
        cfl: 0.1,
        maximum_halvings: 0,
        reconstruction: Reconstruction::FirstOrder,
    };
    let report = flow.step(Time::new(1.0), controls, || false).unwrap();
    assert_eq!(report.halvings, 0);
    assert!(report.elapsed.value() <= controls.cfl / rate * (1.0 + 2e-14));
    assert!(report.maximum_cfl <= controls.cfl * (1.0 + 2e-14));
    assert!(flow.finite_reservoir(side).unwrap().state().density < 1.0);
}

#[test]
fn g1_finite_tank_pressure_feedback_has_second_order_time_convergence() {
    let side = BoundarySide::ZMax;
    let end = 0.02;
    let solve = |steps: usize| {
        let mut flow = finite(side, 2.0, 0.01, [0.0; 3]);
        for _ in 0..steps {
            let dt = end / steps as f64;
            let report = flow
                .step(
                    Time::new(dt),
                    StepControl {
                        reconstruction: Reconstruction::FirstOrder,
                        ..StepControl::default()
                    },
                    || false,
                )
                .unwrap();
            assert_eq!(report.elapsed.value(), dt);
            assert_eq!(report.halvings, 0);
        }
        let tank = flow.finite_reservoir(side).unwrap().state();
        let q = flow.states()[0];
        [
            tank.density,
            tank.energy,
            q.density,
            q.momentum[2],
            q.energy,
        ]
    };
    let reference = solve(512);
    let errors = [16, 32, 64].map(|n| {
        solve(n)
            .into_iter()
            .zip(reference)
            .map(|(a, b)| (a - b).powi(2))
            .sum::<f64>()
            .sqrt()
    });
    assert!(errors[0] / errors[1] > 3.5, "{errors:?}");
    assert!(errors[1] / errors[2] > 3.5, "{errors:?}");
}

#[test]
fn g4_finite_tanks_admit_sources_and_roll_back_the_last_poll() {
    let side = BoundarySide::XMax;
    for volume_m3 in [0.0, -1.0, f64::INFINITY, f64::NAN] {
        assert!(
            held(side, 2.0, [0.0; 3])
                .with_finite_reservoirs(&[FiniteReservoirBoundary { side, volume_m3 }], || false,)
                .is_err()
        );
    }
    let binding = FiniteReservoirBoundary {
        side,
        volume_m3: 0.04,
    };
    assert!(
        held(side, 2.0, [0.0; 3])
            .with_finite_reservoirs(&[binding, binding], || false)
            .is_err()
    );
    assert!(
        held(side, 2.0, [0.0; 3])
            .with_finite_reservoirs(
                &[FiniteReservoirBoundary {
                    side: BoundarySide::YMin,
                    ..binding
                }],
                || false,
            )
            .is_err()
    );
    assert!(
        finite(side, 2.0, 0.04, [0.0; 3])
            .with_stagnation_reservoirs(&[], || false)
            .is_err()
    );
    let mut accepted = finite(side, 2.0, 0.04, [0.0; 3]);
    let mut polls = 0;
    let expected = accepted
        .step(Time::new(0.001), StepControl::default(), || {
            polls += 1;
            false
        })
        .unwrap();
    assert!(polls > 10);
    let mut flow = finite(side, 2.0, 0.04, [0.0; 3]);
    let states = flow.states().to_vec();
    let tank = *flow.finite_reservoir(side).unwrap();
    let exterior = flow.reservoir(side);
    let impulses = flow.boundary_impulses(side).to_vec();
    let mut count = 0;
    assert_eq!(
        flow.step(Time::new(0.001), StepControl::default(), || {
            count += 1;
            count == polls
        }),
        Err(Error::Cancelled)
    );
    assert_eq!(flow.states(), states);
    assert_eq!(flow.finite_reservoir(side), Some(&tank));
    assert_eq!(flow.reservoir(side), exterior);
    assert_eq!(flow.boundary_impulses(side), impulses);
    assert_eq!(flow.time().value(), 0.0);
    assert!(flow.last_step().is_none());
    assert!(flow.finite_reservoir_report(side).is_none());
    assert_eq!(
        flow.step(Time::new(0.001), StepControl::default(), || false)
            .unwrap(),
        expected
    );
    assert_eq!(flow.states(), accepted.states());
    assert_eq!(flow.finite_reservoir(side), accepted.finite_reservoir(side));
    assert_eq!(
        flow.finite_reservoir_report(side),
        accepted.finite_reservoir_report(side)
    );
}

#[test]
fn g0_finite_pressure_storage_composes_with_native_force_work_and_timed_profiles() {
    let side = BoundarySide::XMin;
    let source = gas(0.4);
    let grid = CartesianGrid::new(2, 2, 2, 0.2, 0.2, 0.2).unwrap();
    let state = State::from_primitive(&source, 1.0, [0.0; 3], 1.0).unwrap();
    let mut flow = Euler3d::new_with_reservoirs(
        source,
        grid,
        vec![state; grid.len()],
        Boundaries3d {
            x: AxisBoundaries::Ends {
                lower: BoundaryCondition::Transmissive,
                upper: BoundaryCondition::Transmissive,
            },
            y: AxisBoundaries::Periodic,
            z: AxisBoundaries::Periodic,
        },
        &[ReservoirBoundary {
            side,
            density: 1.0,
            velocity: [0.0; 3],
            pressure_pa: 2.0,
        }],
        || false,
    )
    .unwrap()
    .with_finite_reservoirs(
        &[FiniteReservoirBoundary {
            side,
            volume_m3: 0.04,
        }],
        || false,
    )
    .unwrap();
    let other = BoundarySide::XMax;
    let bytes = flow
        .required_reservoir_profile_workspace_bytes(&[other])
        .unwrap();
    let exterior = ReservoirFaceState {
        density: 1.0,
        velocity: [0.0; 3],
        pressure_pa: 1.0,
    };
    flow.prepare_reservoir_profiles(&[other], bytes, |_, _, _| Ok(exterior), || false)
        .unwrap();
    let initial = storage(&flow, side);
    let mut work = flow
        .forcing_workspace(flow.required_forcing_workspace_bytes().unwrap())
        .unwrap();
    let mut times = Vec::new();
    let report = flow
        .step_with_forcing_and_reservoir_profiles(
            Time::new(0.0001),
            StepControl::default(),
            BodyForce3d::UniformAcceleration([0.2, -0.1, 0.3]),
            ForceControl3d { maximum_mach: 2.0 },
            &mut work,
            1000,
            |t, _, _, _| {
                times.push(t);
                Ok(ReservoirFaceState {
                    pressure_pa: 1.0 + t,
                    ..exterior
                })
            },
            || false,
        )
        .unwrap();
    assert_eq!(report.halvings, 0);
    assert!(times.contains(&0.0));
    assert!(times.contains(&report.elapsed.value()));
    let other_transfer = sum_arrays(
        flow.boundary_impulses(other).iter().map(|q| q.array()),
        &mut || false,
    )
    .unwrap();
    let final_storage = storage(&flow, side);
    close(final_storage[0] - initial[0], -other_transfer[0], 3e-12);
    close(
        final_storage[1] - initial[1],
        report.source_change.work - other_transfer[4],
        3e-12,
    );
    assert!(flow.finite_reservoir(side).unwrap().state().density < 1.0);
}
