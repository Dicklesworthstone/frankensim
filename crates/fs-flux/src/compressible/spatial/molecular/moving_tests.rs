use super::*;

fn channel(axis: usize) -> Boundaries3d {
    let mut boundaries = PERIODIC;
    match axis {
        0 => boundaries.x = WALL,
        1 => boundaries.y = WALL,
        2 => boundaries.z = WALL,
        _ => unreachable!(),
    }
    boundaries
}

fn prepare(
    operator: &mut MolecularTransport3d<ConstantTransport>,
    flow: &Euler3d,
    motion: MeshMotion,
) -> WallTransfer {
    let mut power = WallTransfer::default();
    let rate = spatial_with_motion(
        &operator.gas,
        &flow.grid,
        operator.boundaries,
        &operator.walls,
        Some(motion),
        Some(&mut power),
        None,
        &flow.states,
        &mut operator.law,
        &mut operator.work.primitive,
        &mut operator.work.gradient,
        &mut operator.work.viscosity,
        &mut operator.work.flux,
        &mut operator.work.binary,
        #[cfg(feature = "experimental-sgs")]
        operator.sgs.as_mut(),
        &mut || false,
    )
    .unwrap();
    assert!(rate.is_finite() && rate > 0.0);
    power
}

#[test]
fn g1_moving_normal_stress_and_full_wall_work_match_independent_all_side_laws() {
    let spacing = [0.2, 0.3, 0.4];
    let grid = CartesianGrid::new(3, 3, 3, spacing[0], spacing[1], spacing[2]).unwrap();
    let (mu, bulk, conductivity, slip, fraction) = (0.02, 0.009, 0.05, 0.12, 0.35);
    for side in BoundarySide::ALL {
        let axis = side.axis();
        let [a, b] = tangents(axis);
        let velocity = [0.15, -0.1, 0.07];
        let mut wall_velocity = [0.0; 3];
        wall_velocity[a] = 0.4;
        wall_velocity[b] = -0.3;
        let mut motion = MeshMotion::default();
        motion.lower_velocity[axis] = 0.23;
        motion.upper_velocity[axis] = 0.23;
        for mode in 0..3 {
            let model = match mode {
                0 => MolecularWall3d::SlipAdiabatic,
                1 => MolecularWall3d::NoSlip {
                    velocity: wall_velocity,
                    thermal: WallThermalCondition::Isothermal(Temperature::new(3.0)),
                },
                _ => MolecularWall3d::NavierSlip {
                    velocity: wall_velocity,
                    thermal: WallThermalCondition::Isothermal(Temperature::new(3.0)),
                    slip_length_m: slip,
                    fluid_heat_fraction: fraction,
                },
            };
            let flow = field(gas(0.0), grid, channel(axis), |_| {
                [1.0, velocity[0], velocity[1], velocity[2], 2.2]
            });
            let mut operator =
                MolecularTransport3d::new(&flow, law(flow.gas(), mu, bulk, conductivity), || false)
                    .unwrap()
                    .with_wall(side, model)
                    .unwrap();
            let power = prepare(&mut operator, &flow, motion);
            let distance = 0.5 * spacing[axis];
            let mut inward_traction = [0.0; 3];
            inward_traction[axis] =
                (4.0 * mu / 3.0 + bulk) * (motion.lower_velocity[axis] - velocity[axis]) / distance;
            let mut wall = [0.0; 3];
            wall[axis] = motion.lower_velocity[axis];
            let mut dissipation = 0.0;
            let mut heat = 0.0;
            if mode != 0 {
                let length = if mode == 1 { 0.0 } else { slip };
                for tangent in [a, b] {
                    wall[tangent] = wall_velocity[tangent];
                    let relative = wall[tangent] - velocity[tangent];
                    inward_traction[tangent] = mu * relative / (distance + length);
                    dissipation += mu * length * (relative / (distance + length)).powi(2);
                }
                heat = conductivity * (3.0 - 2.2) / distance - (1.0 - fraction) * dissipation;
            }
            let work: f64 = (0..3).map(|d| inward_traction[d] * wall[d]).sum();
            let mut face = [1; 3];
            face[axis] = if side.sign() > 0.0 { 3 } else { 0 };
            let flux = operator.work.flux[axis][grid.cell_face(axis, face)];
            close(
                -side.sign() * flux.normal_momentum,
                inward_traction[axis],
                3e-14,
            );
            for (i, tangent) in [a, b].into_iter().enumerate() {
                close(
                    -side.sign() * flux.transverse_momentum[i],
                    inward_traction[tangent],
                    3e-14,
                );
            }
            close(-side.sign() * flux.energy, work + heat, 4e-14);
            let area = 9.0 * grid.face_area(side);
            let index = 2 * axis + usize::from(side.sign() > 0.0);
            close(power.work[index], area * work, 4e-14);
            close(power.heat[index], area * heat, 4e-14);
        }
    }
}

#[test]
fn g1_uniform_comoving_molecular_box_preserves_all_momenta_for_every_slip_law() {
    let velocity = [0.3, -0.2, 0.15];
    let motion = MeshMotion {
        lower_velocity: velocity,
        upper_velocity: velocity,
    };
    let grid = CartesianGrid::new(3, 2, 2, 0.25, 0.3, 0.4).unwrap();
    for mode in 0..3 {
        let flow = field(gas(0.0), grid, WALLS, |_| {
            [1.0, velocity[0], velocity[1], velocity[2], 2.0]
        });
        let initial = flow.states().to_vec();
        let transport = law(flow.gas(), 0.03, 0.02, 0.05);
        let mut owner = MolecularEuler3d::new(flow, transport, || false).unwrap();
        if mode != 0 {
            for side in BoundarySide::ALL {
                let mut wall = velocity;
                wall[side.axis()] = 0.0;
                owner = owner
                    .with_wall(
                        side,
                        if mode == 1 {
                            MolecularWall3d::NoSlip {
                                velocity: wall,
                                thermal: WallThermalCondition::Adiabatic,
                            }
                        } else {
                            MolecularWall3d::NavierSlip {
                                velocity: wall,
                                thermal: WallThermalCondition::Adiabatic,
                                slip_length_m: 0.1,
                                fluid_heat_fraction: 0.4,
                            }
                        },
                    )
                    .unwrap();
            }
        }
        let report = owner
            .step_moving(
                Time::new(0.002),
                MolecularStepControl::default(),
                motion,
                || false,
            )
            .unwrap();
        for (actual, expected) in owner.flow().states().iter().zip(&initial) {
            for (a, e) in actual.array().into_iter().zip(expected.array()) {
                close(a, e, 4e-13);
            }
        }
        for axis in 0..3 {
            close(
                owner.flow().grid.origin[axis],
                report.transport.elapsed.value() * velocity[axis],
                2e-14,
            );
        }
        assert_eq!(owner.flow().grid.spacing, grid.spacing);
        for value in report
            .molecular_wall_work_j
            .into_iter()
            .chain(report.wall_heat_j)
        {
            close(value, 0.0, 3e-15);
        }
        for value in report.transport.conservation_residual.array() {
            close(value, 0.0, 3e-13);
        }
    }
}

fn heated_box() -> MolecularEuler3d<ConstantTransport> {
    let grid = CartesianGrid::new(3, 3, 2, 0.3, 0.2, 0.4).unwrap();
    let flow = field(gas(0.0), grid, WALLS, |x| {
        [
            1.0 + 0.03 * x[0],
            0.1 + 0.04 * x[0],
            -0.07 + 0.03 * x[1],
            0.05,
            2.0 + 0.1 * x[2],
        ]
    });
    let transport = law(flow.gas(), 0.03, 0.01, 0.07);
    let mut owner = MolecularEuler3d::new(flow, transport, || false).unwrap();
    for side in BoundarySide::ALL {
        let mut velocity = [0.04, -0.03, 0.02];
        velocity[side.axis()] = 0.0;
        owner = owner
            .with_wall(
                side,
                MolecularWall3d::NavierSlip {
                    velocity,
                    thermal: WallThermalCondition::Isothermal(Temperature::new(3.0)),
                    slip_length_m: 0.08,
                    fluid_heat_fraction: 0.6,
                },
            )
            .unwrap();
    }
    owner
}

#[test]
fn g1_zero_motion_keeps_stationary_split_exact_and_retains_nonzero_heat_and_work() {
    let mut stationary = heated_box();
    let mut moving = heated_box();
    let controls = MolecularStepControl::default();
    let expected = stationary
        .step(Time::new(0.001), controls, || false)
        .unwrap();
    let actual = moving
        .step_moving(Time::new(0.001), controls, MeshMotion::default(), || false)
        .unwrap();
    assert_eq!(actual.transport, expected);
    assert_eq!(moving.flow().states(), stationary.flow().states());
    assert_eq!(moving.flow().grid(), stationary.flow().grid());
    for side in BoundarySide::ALL {
        assert_eq!(
            moving.boundary_impulses(side),
            stationary.boundary_impulses(side)
        );
    }
    assert!(actual.wall_heat_j.iter().all(|q| *q > 0.0));
    assert!(actual.molecular_wall_work_j.iter().any(|q| q.abs() > 1e-8));
    close(
        actual.molecular_wall_work_j.into_iter().sum::<f64>()
            + actual.wall_heat_j.into_iter().sum::<f64>(),
        expected.first.boundary_change.energy + expected.second.boundary_change.energy,
        4e-13,
    );
}

#[test]
fn g2_g4_moving_split_uses_new_volume_and_late_cancellation_replays_atomically() {
    let motion = MeshMotion {
        lower_velocity: [-0.05, 0.02, 0.0],
        upper_velocity: [0.12, 0.02, 0.0],
    };
    let controls = MolecularStepControl::default();
    let duration = Time::new(0.001);
    let mut reference = heated_box();
    reference
        .step_moving(duration, controls, motion, || false)
        .unwrap();
    let old_grid = *reference.flow().grid();
    let mut polls = 0;
    let expected = reference
        .step_moving(duration, controls, motion, || {
            polls += 1;
            false
        })
        .unwrap();
    let actual_grid = reference.flow().grid();
    close(
        actual_grid.spacing[0],
        old_grid.spacing[0] + duration.value() * 0.17 / 3.0,
        2e-14,
    );
    close(
        expected.transport.before.mass,
        expected.transport.after.mass,
        3e-13,
    );
    close(
        expected.transport.second.before.mass,
        expected.transport.second.after.mass,
        3e-13,
    );
    for value in expected.transport.conservation_residual.array() {
        close(value, 0.0, 4e-13);
    }
    let molecular_energy = expected.transport.first.boundary_change.energy
        + expected.transport.second.boundary_change.energy;
    close(
        expected.molecular_wall_work_j.into_iter().sum::<f64>()
            + expected.wall_heat_j.into_iter().sum::<f64>(),
        molecular_energy,
        4e-13,
    );
    for cutoff in [polls / 2, polls] {
        let mut owner = heated_box();
        owner
            .step_moving(duration, controls, motion, || false)
            .unwrap();
        let states = owner.flow().states().to_vec();
        let grid = *owner.flow().grid();
        let time = owner.flow().time();
        let report = owner.last_step().copied();
        let euler_report = owner.flow().last_step().copied();
        let loads = BoundarySide::ALL.map(|side| owner.boundary_impulses(side).to_vec());
        let mut calls = 0;
        assert_eq!(
            owner.step_moving(duration, controls, motion, || {
                calls += 1;
                calls == cutoff
            }),
            Err(MolecularError::Cancelled)
        );
        assert_eq!(owner.flow().states(), states);
        assert_eq!(*owner.flow().grid(), grid);
        assert_eq!(owner.flow().time(), time);
        assert_eq!(owner.last_step().copied(), report);
        assert_eq!(owner.flow().last_step().copied(), euler_report);
        for (i, side) in BoundarySide::ALL.into_iter().enumerate() {
            assert_eq!(owner.boundary_impulses(side), loads[i]);
        }
        let retry = owner
            .step_moving(duration, controls, motion, || false)
            .unwrap();
        assert_eq!(retry, expected);
        assert_eq!(owner.flow().states(), reference.flow().states());
    }
    // A refused moving proposal must not leave private geometry attached to
    // the next ordinary stationary split, including after prior box motion.
    let mut cancelled = heated_box();
    let mut stationary = heated_box();
    for owner in [&mut cancelled, &mut stationary] {
        owner
            .step_moving(duration, controls, motion, || false)
            .unwrap();
    }
    let mut calls = 0;
    assert_eq!(
        cancelled.step_moving(duration, controls, motion, || {
            calls += 1;
            calls == polls
        }),
        Err(MolecularError::Cancelled)
    );
    let a = cancelled.step(duration, controls, || false).unwrap();
    let b = stationary.step(duration, controls, || false).unwrap();
    assert_eq!(a, b);
    assert_eq!(cancelled.flow().states(), stationary.flow().states());
}

#[test]
fn g2_moving_split_restarts_both_geometry_halves_after_molecular_work_refusal() {
    let motion = MeshMotion {
        lower_velocity: [-0.05, 0.0, 0.0],
        upper_velocity: [0.12, 0.0, 0.0],
    };
    let controls = MolecularStepControl {
        molecular: MolecularControl {
            maximum_step_s: 0.001,
            maximum_steps: 1,
            ..MolecularControl::default()
        },
        ..MolecularStepControl::default()
    };
    let mut retried = heated_box();
    let report = retried
        .step_moving(Time::new(0.008), controls, motion, || false)
        .unwrap();
    assert!(report.transport.coupling_retries >= 2);
    assert_eq!(
        report.transport.first.elapsed.value() * 2.0,
        report.transport.elapsed.value()
    );
    assert_eq!(
        report.transport.second.elapsed,
        report.transport.first.elapsed
    );
    assert_eq!(report.transport.advection.elapsed, report.transport.elapsed);
    let mut direct = heated_box();
    let expected = direct
        .step_moving(report.transport.elapsed, controls, motion, || false)
        .unwrap();
    assert_eq!(retried.flow().states(), direct.flow().states());
    assert_eq!(retried.flow().grid(), direct.flow().grid());
    assert_eq!(report.transport.first, expected.transport.first);
    assert_eq!(report.transport.advection, expected.transport.advection);
    assert_eq!(report.transport.second, expected.transport.second);
    assert_eq!(report.molecular_wall_work_j, expected.molecular_wall_work_j);
    assert_eq!(report.wall_heat_j, expected.wall_heat_j);
    for side in BoundarySide::ALL {
        assert_eq!(
            retried.boundary_impulses(side),
            direct.boundary_impulses(side)
        );
    }
}

#[test]
fn g1_affine_expansion_recovers_physical_viscous_heating_under_refinement() {
    // R=1, cv=2.5. Uniform rho and T, u=V*x/L, constant molecular viscosities:
    // T(L)=eta*V + (T0-eta*V)*L^(-R/cv), eta=4*mu/3+bulk, L0=rho0=1.
    // This follows the independent homogeneous internal-energy balance;
    // pressure work and viscous heating must both be present in the solver.
    let (mu, bulk, speed, end) = (0.08, 0.025, 0.35, 0.08);
    let eta = 4.0 * mu / 3.0 + bulk;
    let length: f64 = 1.0 + speed * end;
    let expected = eta * speed + (2.0 - eta * speed) * length.powf(-0.4);
    let inviscid = 2.0 * length.powf(-0.4);
    let mut errors = Vec::new();
    for n in [6, 12, 24] {
        let grid = CartesianGrid::new(n, 1, 1, 1.0 / n as f64, 1.0, 1.0).unwrap();
        let flow = field(gas(0.0), grid, channel(0), |x| {
            [1.0, speed * x[0], 0.0, 0.0, 2.0]
        });
        let transport = law(flow.gas(), mu, bulk, 0.03);
        let mut owner = MolecularEuler3d::new(flow, transport, || false).unwrap();
        let motion = MeshMotion {
            lower_velocity: [0.0; 3],
            upper_velocity: [speed, 0.0, 0.0],
        };
        let controls = MolecularStepControl {
            advection: StepControl {
                cfl: 0.25,
                ..StepControl::default()
            },
            ..MolecularStepControl::default()
        };
        while owner.flow().time().value() < end {
            let remaining = end - owner.flow().time().value();
            owner
                .step_moving(Time::new(remaining), controls, motion, || false)
                .unwrap();
        }
        let mean_temperature = owner
            .flow()
            .states()
            .iter()
            .map(|q| {
                q.primitive(owner.flow().gas(), || false)
                    .unwrap()
                    .temperature
            })
            .sum::<f64>()
            / n as f64;
        assert!(mean_temperature > inviscid);
        errors.push((mean_temperature - expected).abs());
        close(owner.flow().totals().mass, 1.0, 2e-12);
    }
    assert!(
        errors[0] > 2.5 * errors[1] && errors[1] > 2.5 * errors[2],
        "errors={errors:?}"
    );
    assert!(errors[2] < 2e-5, "errors={errors:?}");
}
