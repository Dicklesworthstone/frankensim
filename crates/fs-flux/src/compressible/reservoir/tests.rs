use super::*;
use crate::compressible::tests::{close, gas};

fn field(side: BoundarySide, velocity: [f64; 2]) -> Euler2d {
    let source = gas(0.4);
    let grid = CartesianGrid::new(2, 2, 0.2, 0.3, 0.4).unwrap();
    let initial = State::from_primitive(&source, 1.0, velocity, 1.0).unwrap();
    let mut pairs = [AxisBoundaries::Periodic; 2];
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
    Euler2d::new(
        source,
        grid,
        vec![initial; grid.len()],
        Boundaries2d {
            x: pairs[0],
            y: pairs[1],
        },
        || false,
    )
    .unwrap()
}

fn finite(side: BoundarySide, pressure: f64, volume: f64, velocity: [f64; 2]) -> Euler2d {
    let flow = field(side, velocity);
    let tank = Plenum::new(volume, flow.gas().conserved(1.0, 0.0, pressure).unwrap()).unwrap();
    flow.with_finite_reservoirs(&[FiniteReservoirBoundary { side, tank }], || false)
        .unwrap()
}

fn storage(flow: &Euler2d) -> [f64; 3] {
    let field = flow.totals();
    let mut total = [field.mass, field.energy];
    for side in SIDES {
        if let Some(tank) = flow.finite_reservoir(side) {
            total[0] += tank.state().density * tank.volume_m3();
            total[1] += tank.state().energy * tank.volume_m3();
        }
    }
    [
        total[0],
        total[1],
        total[1] + flow.gas().energy_reference() * total[0],
    ]
}

#[test]
fn g0_finite_planar_tanks_fill_drain_and_conserve_each_orientation_and_order() {
    for side in SIDES {
        for reconstruction in [Reconstruction::FirstOrder, Reconstruction::Muscl] {
            for pressure in [0.5, 2.0] {
                let mut velocity = [0.0; 2];
                if pressure < 1.0 {
                    velocity[1 - side.axis()] = 0.3;
                }
                let mut flow = finite(side, pressure, 0.04, velocity);
                let before = storage(&flow);
                let original = *flow.finite_reservoir(side).unwrap();
                let initial_temperature = flow
                    .gas()
                    .primitive(original.state(), || false)
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
                    for (got, want) in storage(&flow).into_iter().zip(before) {
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
                            tank.transfer.momentum[1 - side.axis()],
                            velocity[1 - side.axis()] * tank.transfer.mass,
                            3e-12,
                        );
                    }
                    for (d, value) in tank.transfer.array().into_iter().enumerate() {
                        close(
                            value,
                            flow.boundary_impulses(side)
                                .iter()
                                .map(|q| q.array()[d])
                                .sum(),
                            2e-14,
                        );
                    }
                }
                let p = flow
                    .gas()
                    .primitive(flow.finite_reservoir(side).unwrap().state(), || false)
                    .unwrap();
                if pressure > 1.0 {
                    assert!(p.pressure < pressure);
                    assert!(p.temperature < initial_temperature);
                } else {
                    assert!(p.pressure > pressure);
                }
            }
        }
    }
}

#[test]
fn g0_equilibrium_and_empty_bindings_preserve_native_bits() {
    let side = BoundarySide::YMax;
    let mut tank = finite(side, 1.0, 0.04, [0.0; 2]);
    let mut ordinary = field(side, [0.0; 2]);
    let original = *tank.finite_reservoir(side).unwrap();
    for _ in 0..3 {
        let a = tank
            .step(Time::new(0.001), StepControl::default(), || false)
            .unwrap();
        let b = ordinary
            .step(Time::new(0.001), StepControl::default(), || false)
            .unwrap();
        assert_eq!(a, b);
        assert_eq!(tank.states(), ordinary.states());
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
    let mut ordinary = field(side, [0.2, 0.1]);
    let mut empty = field(side, [0.2, 0.1])
        .with_finite_reservoirs(&[], || false)
        .unwrap();
    assert_eq!(
        ordinary
            .step(Time::new(0.001), StepControl::default(), || false)
            .unwrap(),
        empty
            .step(Time::new(0.001), StepControl::default(), || false)
            .unwrap()
    );
    assert_eq!(ordinary.states(), empty.states());
}

#[test]
fn g1_tiny_planar_tank_admits_its_aperture_cfl_before_retry() {
    let side = BoundarySide::XMin;
    let volume = 1e-6;
    let mut flow = finite(side, 2.0, volume, [0.0; 2]);
    let tank = flow.finite_reservoir(side).unwrap().state();
    let flux = flow
        .gas()
        .reservoir_normal_flux(flow.states()[0].normal(0), tank, -1.0, || false)
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
fn g1_planar_tank_feedback_uses_both_rk_stages() {
    let side = BoundarySide::XMax;
    let solve = |steps: usize| {
        let mut flow = finite(side, 2.0, 0.01, [0.0; 2]);
        let dt = 0.02 / steps as f64;
        for _ in 0..steps {
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
        let q = flow.states()[1];
        [
            tank.density,
            tank.energy,
            q.density,
            q.momentum[0],
            q.energy,
        ]
    };
    let reference = solve(512);
    let errors = [16, 32, 64].map(|steps| {
        solve(steps)
            .into_iter()
            .zip(reference)
            .map(|(a, b)| (a - b).powi(2))
            .sum::<f64>()
            .sqrt()
    });
    eprintln!("planar finite-tank time errors: {errors:?}");
    assert!(errors[0] / errors[1] > 3.5, "{errors:?}");
    assert!(errors[1] / errors[2] > 3.5, "{errors:?}");
}

#[test]
fn g0_two_independent_tanks_share_field_energy_and_native_body_work() {
    let source = gas(0.4);
    let grid = CartesianGrid::new(3, 2, 0.2, 0.3, 0.4).unwrap();
    let q = State::from_primitive(&source, 1.0, [0.1, 0.2], 1.0).unwrap();
    let tanks =
        [2.0, 0.5].map(|p| Plenum::new(0.04, source.conserved(1.0, 0.0, p).unwrap()).unwrap());
    let mut flow = Euler2d::new(
        source,
        grid,
        vec![q; grid.len()],
        Boundaries2d {
            x: AxisBoundaries::Ends {
                lower: BoundaryCondition::Transmissive,
                upper: BoundaryCondition::Transmissive,
            },
            y: AxisBoundaries::Periodic,
        },
        || false,
    )
    .unwrap()
    .with_finite_reservoirs(
        &[
            FiniteReservoirBoundary {
                side: BoundarySide::XMin,
                tank: tanks[0],
            },
            FiniteReservoirBoundary {
                side: BoundarySide::XMax,
                tank: tanks[1],
            },
        ],
        || false,
    )
    .unwrap();
    let before = storage(&flow);
    let mut work = flow
        .forcing_workspace(flow.required_forcing_workspace_bytes().unwrap())
        .unwrap();
    let report = flow
        .step_with_forcing(
            Time::new(0.0001),
            StepControl::default(),
            BodyForce2d::UniformAcceleration([0.2, -0.1]),
            ForceControl2d { maximum_mach: 2.0 },
            &mut work,
            || false,
        )
        .unwrap();
    let after = storage(&flow);
    close(after[0], before[0], 3e-12);
    close(after[1] - before[1], report.source_change.work, 3e-12);
    close(after[2] - before[2], report.source_change.work, 3e-12);
    assert!(
        flow.finite_reservoir_report(BoundarySide::XMin)
            .unwrap()
            .transfer
            .mass
            < 0.0
    );
    assert!(
        flow.finite_reservoir_report(BoundarySide::XMax)
            .unwrap()
            .transfer
            .mass
            > 0.0
    );
    for residual in report.conservation_residual.array() {
        close(residual, 0.0, 3e-12);
    }
}

#[test]
fn g4_late_cancellation_preserves_accepted_tank_and_retry_feedback() {
    let side = BoundarySide::XMax;
    let prepared = || {
        let mut flow = finite(side, 2.0, 0.04, [0.0; 2]);
        flow.step(Time::new(0.001), StepControl::default(), || false)
            .unwrap();
        flow
    };
    let mut accepted = prepared();
    let mut polls = 0;
    let expected = accepted
        .step(Time::new(0.001), StepControl::default(), || {
            polls += 1;
            false
        })
        .unwrap();
    let mut flow = prepared();
    let states = flow.states().to_vec();
    let tank = *flow.finite_reservoir(side).unwrap();
    let tank_report = *flow.finite_reservoir_report(side).unwrap();
    let report = *flow.last_step().unwrap();
    let impulses = flow.boundary_impulses(side).to_vec();
    let time = flow.time();
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
    assert_eq!(flow.finite_reservoir_report(side), Some(&tank_report));
    assert_eq!(flow.boundary_impulses(side), impulses);
    assert_eq!(flow.time(), time);
    assert_eq!(flow.last_step(), Some(&report));
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
fn g4_binding_admits_source_and_refuses_undeclared_gravitational_ports() {
    use fs_qty::{Acceleration, Length, Pressure, Temperature};
    let side = BoundarySide::XMin;
    let source = gas(0.4);
    let tank = Plenum::new(0.04, source.conserved(1.0, 0.0, 2.0).unwrap()).unwrap();
    let binding = FiniteReservoirBoundary { side, tank };
    assert!(
        field(side, [0.0; 2])
            .with_finite_reservoirs(&[binding, binding], || false)
            .is_err()
    );
    for wrong in [BoundarySide::XMax, BoundarySide::YMin] {
        assert!(
            field(side, [0.0; 2])
                .with_finite_reservoirs(&[FiniteReservoirBoundary { side: wrong, tank }], || false)
                .is_err()
        );
    }
    let unsupported = Plenum::new(
        0.04,
        fs_gas::dynamics::Conserved {
            density: 1.0,
            momentum: 0.0,
            energy: 1e300,
        },
    )
    .unwrap();
    assert!(
        field(side, [0.0; 2])
            .with_finite_reservoirs(
                &[FiniteReservoirBoundary {
                    side,
                    tank: unsupported
                }],
                || false
            )
            .is_err()
    );
    let atmosphere = IsothermalAtmosphere {
        acceleration: Acceleration::new(-0.1),
        temperature: Temperature::new(1.0),
        reference_pressure: Pressure::new(1.0),
        reference_position: Length::new(0.0),
    };
    let mut flow = finite(side, 2.0, 0.04, [0.0; 2]);
    assert!(
        flow.bind_isothermal_atmosphere(0, atmosphere, || false)
            .is_err()
    );
    assert!(flow.hydrostatic_reference().is_none());
    let mut stratified = field(side, [0.0; 2]);
    stratified
        .bind_isothermal_atmosphere(0, atmosphere, || false)
        .unwrap();
    assert!(
        stratified
            .with_finite_reservoirs(&[binding], || false)
            .is_err()
    );
    let mut evolved = field(side, [0.0; 2]);
    evolved
        .step(Time::new(0.001), StepControl::default(), || false)
        .unwrap();
    assert!(
        evolved
            .with_finite_reservoirs(&[binding], || false)
            .is_err()
    );
}
