use super::*;

fn bound(
    side: BoundarySide,
    volume_m3: f64,
    tank_pressure: f64,
    field_pressure: f64,
    velocity: [f64; 3],
) -> MixtureEuler3d {
    let mut flow = supplied(
        standard(0.2),
        CartesianGrid::new(2, 2, 2, 0.1, 0.1, 0.1).unwrap(),
        one_open_side(side),
        ([1.0, 0.0, 0.0], velocity, field_pressure),
        &[ReservoirBoundary {
            side,
            partial_densities: &[0.0, 1.0, 0.0],
            velocity: [0.0; 3],
            pressure_pa: tank_pressure,
        }],
    );
    let bytes = flow
        .required_workspace_bytes_with_finite_reservoirs(1)
        .unwrap();
    flow = flow
        .with_finite_reservoirs(&[FiniteReservoir { side, volume_m3 }], bytes, || false)
        .unwrap();
    assert_eq!(flow.workspace_bytes(), bytes);
    flow
}

fn combined(flow: &MixtureEuler3d, side: BoundarySide) -> ([f64; 3], f64) {
    let totals = flow.totals();
    let tank = flow.finite_reservoir(side).unwrap();
    let q = tank.state();
    let volume = tank.volume_m3();
    let species =
        std::array::from_fn(|s| totals.species_masses[s] + volume * q.partial_densities[s]);
    let tank_absolute = q.energy
        + q.partial_densities
            .iter()
            .zip(flow.gas().energy_references())
            .map(|(rho, reference)| rho * reference)
            .sum::<f64>();
    (
        species,
        totals.absolute_energy(flow.gas()).unwrap() + volume * tank_absolute,
    )
}

#[test]
fn g1_finite_volumetric_tanks_conserve_original_species_and_energy_in_both_flow_directions() {
    for reconstruction in [Reconstruction::FirstOrder, Reconstruction::Muscl] {
        for side in BoundarySide::ALL {
            for fills in [false, true] {
                let pressure = if fills { 0.7 } else { 1.5 };
                let mut flow = bound(side, 0.004, pressure, 1.0, [0.0; 3]);
                let initial = combined(&flow, side);
                let capacity = flow.workspace_bytes();
                for _ in 0..6 {
                    let step = flow
                        .step(
                            Time::new(0.001),
                            StepControl {
                                reconstruction,
                                ..StepControl::default()
                            },
                            || false,
                        )
                        .unwrap();
                    balance(step);
                    let now = combined(&flow, side);
                    for s in 0..3 {
                        close(now.0[s], initial.0[s], 5e-12);
                    }
                    close(now.1, initial.1, 5e-11);
                    let report = flow.finite_reservoir_report(side).unwrap();
                    for residual in report.species_residual_kg {
                        close(residual, 0.0, 2e-13);
                    }
                    close(report.energy_residual_j, 0.0, 2e-13);
                    close(report.absolute_energy_residual_j, 0.0, 3e-12);
                    let mut support = [0.0; 3];
                    for face in flow.boundary_impulses(side) {
                        for axis in 0..3 {
                            support[axis] += face.momentum[axis];
                        }
                    }
                    for axis in 0..3 {
                        close(report.transfer.momentum[axis], support[axis], 2e-13);
                    }
                    assert_eq!(
                        flow.finite_reservoir(side)
                            .unwrap()
                            .state()
                            .partial_densities[2],
                        0.0
                    );
                    for k in 0..flow.cell_count() {
                        assert_eq!(flow.cell(k).unwrap().partial_densities[2], 0.0);
                    }
                    let exterior = flow.reservoir(side).unwrap();
                    let tank = flow.finite_reservoir(side).unwrap();
                    let primitive = flow.gas().primitive(tank.state(), || false).unwrap();
                    assert_eq!(exterior.partial_densities, tank.state().partial_densities);
                    close(exterior.pressure_pa, primitive.pressure, 2e-14);
                    assert_eq!(exterior.velocity, [0.0; 3]);
                    assert_eq!(flow.workspace_bytes(), capacity);
                }
                let tank = flow.finite_reservoir(side).unwrap();
                let p = flow.reservoir(side).unwrap().pressure_pa;
                if fills {
                    assert!(tank.state().partial_densities[0] > 0.0);
                    assert!(p > pressure);
                } else {
                    assert!(tank.state().partial_densities[1] < 1.0);
                    assert!(p < pressure);
                }
            }
        }
    }
}

#[test]
fn g1_finite_volumetric_tank_absorbs_shear_energy_and_native_body_work() {
    for side in BoundarySide::ALL {
        let mut velocity = [0.15, -0.2, 0.25];
        velocity[side.axis()] = 0.1 * side.sign();
        let mut flow = bound(side, 0.01, 0.7, 1.0, velocity);
        let before = combined(&flow, side);
        let mut workspace = flow
            .forcing_workspace(flow.required_forcing_workspace_bytes().unwrap())
            .unwrap();
        let step = flow
            .step_with_forcing(
                Time::new(0.001),
                StepControl::default(),
                BodyForce3d::UniformAcceleration([0.0, 0.25, 0.0]),
                ForceControl3d { maximum_mach: 2.0 },
                &mut workspace,
                || false,
            )
            .unwrap();
        balance(step);
        let after = combined(&flow, side);
        for s in 0..3 {
            close(after.0[s], before.0[s], 3e-12);
        }
        close(after.1 - before.1, step.source_change.work, 3e-11);
        let tank = flow.finite_reservoir_report(side).unwrap();
        assert!(tank.transfer.mass() > 0.0);
        for axis in 0..3 {
            if axis != side.axis() {
                assert!(tank.transfer.momentum[axis] * velocity[axis] > 0.0);
            }
        }
        assert_eq!(flow.finite_reservoir(side).unwrap().state().momentum, 0.0);
    }
}

#[test]
fn g1_finite_volumetric_tank_stage_feedback_refines_at_second_order() {
    let run = |dt: f64| {
        let side = BoundarySide::XMin;
        let mut flow = supplied(
            standard(0.2),
            CartesianGrid::new(1, 1, 1, 0.4, 0.4, 0.4).unwrap(),
            one_open_side(side),
            ([1.0, 0.0, 0.0], [0.0; 3], 1.0),
            &[ReservoirBoundary {
                side,
                partial_densities: &[0.0, 1.0, 0.0],
                velocity: [0.0; 3],
                pressure_pa: 1.3,
            }],
        );
        let bytes = flow
            .required_workspace_bytes_with_finite_reservoirs(1)
            .unwrap();
        flow = flow
            .with_finite_reservoirs(
                &[FiniteReservoir {
                    side,
                    volume_m3: 0.02,
                }],
                bytes,
                || false,
            )
            .unwrap();
        let count = (0.03 / dt).round() as usize;
        for _ in 0..count {
            let step = flow
                .step(Time::new(dt), StepControl::default(), || false)
                .unwrap();
            close(step.elapsed.value(), dt, 2e-14);
            assert_eq!(step.halvings, 0);
        }
        let tank = flow.finite_reservoir(side).unwrap();
        [
            tank.state().partial_densities[1],
            flow.reservoir(side).unwrap().pressure_pa,
        ]
    };
    let reference = run(0.00009375);
    let errors: Vec<_> = [0.006, 0.003, 0.0015]
        .into_iter()
        .map(|dt| {
            let result = run(dt);
            (result[0] - reference[0]).abs() + (result[1] - reference[1]).abs()
        })
        .collect();
    eprintln!("finite 3D mixture tank RK2 errors: {errors:?}");
    assert!(errors[0] > 3.2 * errors[1], "{errors:?}");
    assert!(errors[1] > 3.2 * errors[2], "{errors:?}");
}

#[test]
fn g1_small_finite_volumetric_tank_bounds_the_first_interval_by_its_aperture() {
    let side = BoundarySide::ZMax;
    let mut flow = bound(side, 1e-5, 1.5, 1.0, [0.0; 3]);
    let mut flux = [0.0; 7];
    let tank = flow.finite_reservoir(side).unwrap();
    let receipt = flow
        .gas()
        .reservoir_spatial_flux(
            flow.cell(0).unwrap().normal(side.axis()),
            tank.state(),
            side.sign(),
            &mut flux,
            || false,
        )
        .unwrap();
    let rate =
        0.5 * flow.grid.face_count(side) as f64 * flow.grid.face_area(side) * receipt.signal_speed
            / tank.volume_m3();
    let step = flow
        .step(Time::new(0.01), StepControl::default(), || false)
        .unwrap();
    assert!(step.elapsed.value() < 0.001);
    assert!(step.elapsed.value() * rate <= 0.4 * (1.0 + 3e-14));
    balance(step);
    let flow = flow.with_stagnation_reservoirs(&[], || false).unwrap();
    assert!(flow.uses_stagnation_reservoir(side));
    assert!(flow.finite_reservoir(side).is_some());
}

#[test]
fn g4_finite_volumetric_tank_admission_and_late_cancel_preserve_all_accepted_state() {
    let side = BoundarySide::XMin;
    let make_unbound = || {
        supplied(
            standard(0.2),
            CartesianGrid::new(2, 2, 2, 0.1, 0.1, 0.1).unwrap(),
            one_open_side(side),
            ([1.0, 0.0, 0.0], [0.0; 3], 1.0),
            &[ReservoirBoundary {
                side,
                partial_densities: &[0.0, 1.0, 0.0],
                velocity: [0.0; 3],
                pressure_pa: 1.5,
            }],
        )
    };
    let bytes = make_unbound()
        .required_workspace_bytes_with_finite_reservoirs(1)
        .unwrap();
    let binding = FiniteReservoir {
        side,
        volume_m3: 0.004,
    };
    assert!(
        make_unbound()
            .with_finite_reservoirs(&[binding], bytes - 1, || false)
            .is_err()
    );
    for volume_m3 in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(
            make_unbound()
                .with_finite_reservoirs(
                    &[FiniteReservoir {
                        volume_m3,
                        ..binding
                    }],
                    bytes,
                    || false,
                )
                .is_err()
        );
    }
    assert!(
        make_unbound()
            .with_finite_reservoirs(&[binding, binding], MAX_MIXTURE_WORKSPACE_BYTES, || false,)
            .is_err()
    );
    assert!(
        make_unbound()
            .with_finite_reservoirs(
                &[FiniteReservoir {
                    side: BoundarySide::XMax,
                    ..binding
                }],
                bytes,
                || false,
            )
            .is_err()
    );
    let prepared = || {
        let mut flow = make_unbound()
            .with_finite_reservoirs(&[binding], bytes, || false)
            .unwrap();
        flow.step(Time::new(0.001), StepControl::default(), || false)
            .unwrap();
        flow
    };
    let mut probe = prepared();
    let mut calls = 0;
    probe
        .step(Time::new(0.001), StepControl::default(), || {
            calls += 1;
            false
        })
        .unwrap();
    for cutoff in [1, calls / 2, calls - 1, calls] {
        let mut flow = prepared();
        let cells = flow.states().to_vec();
        let time = flow.time();
        let last = *flow.last_step().unwrap();
        let loads = flow.boundary_impulses.clone();
        let tank = flow.finite_reservoir(side).unwrap();
        let densities = tank.state().partial_densities.to_vec();
        let energy = tank.state().energy;
        let pressure = flow.reservoir(side).unwrap().pressure_pa;
        let tank_report = *flow.finite_reservoir_report(side).unwrap();
        let mut polls = 0;
        assert_eq!(
            flow.step(Time::new(0.001), StepControl::default(), || {
                polls += 1;
                polls == cutoff
            }),
            Err(Error::Cancelled)
        );
        assert_eq!(flow.states(), cells);
        assert_eq!(flow.time(), time);
        assert_eq!(flow.last_step(), Some(&last));
        assert_eq!(flow.boundary_impulses, loads);
        let tank = flow.finite_reservoir(side).unwrap();
        assert_eq!(tank.state().partial_densities, densities);
        assert_eq!(tank.state().energy, energy);
        assert_eq!(flow.reservoir(side).unwrap().pressure_pa, pressure);
        assert_eq!(flow.finite_reservoir_report(side), Some(&tank_report));
        flow.step(Time::new(0.001), StepControl::default(), || false)
            .unwrap();
        assert_eq!(flow.states(), probe.states());
        assert_eq!(
            flow.finite_reservoir_report(side),
            probe.finite_reservoir_report(side)
        );
    }
}
