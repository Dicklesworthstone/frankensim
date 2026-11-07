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

fn prepare(op: &mut MolecularTransport3d<ConstantTransport>, flow: &Euler3d) {
    let rate = spatial(
        &op.gas,
        &op.grid,
        op.boundaries,
        &op.walls,
        None,
        &flow.states,
        &mut op.law,
        &mut op.work.primitive,
        &mut op.work.gradient,
        &mut op.work.viscosity,
        &mut op.work.flux,
        &mut op.work.binary,
        #[cfg(feature = "experimental-sgs")]
        op.sgs.as_mut(),
        &mut || false,
    )
    .unwrap();
    assert!(rate.is_finite() && rate > 0.0);
}

#[test]
fn g1_navier_slip_uses_series_traction_and_explicit_interfacial_heat_on_all_sides() {
    let spacing = [0.2, 0.3, 0.4];
    let grid = CartesianGrid::new(3, 3, 3, spacing[0], spacing[1], spacing[2]).unwrap();
    let (mu, conductivity, slip, fraction) = (0.02, 0.05, 0.12, 0.35);
    for side in BoundarySide::ALL {
        let axis = side.axis();
        let [a, b] = tangents(axis);
        let mut cell_velocity = [0.0; 3];
        cell_velocity[a] = 0.15;
        cell_velocity[b] = -0.1;
        let mut wall_velocity = [0.0; 3];
        wall_velocity[a] = 0.4;
        wall_velocity[b] = -0.3;
        let flow = field(gas(0.0), grid, channel(axis), |_| {
            [
                1.0,
                cell_velocity[0],
                cell_velocity[1],
                cell_velocity[2],
                2.2,
            ]
        });
        let mut op =
            MolecularTransport3d::new(&flow, law(flow.gas(), mu, 0.003, conductivity), || false)
                .unwrap()
                .with_wall(
                    side,
                    MolecularWall3d::NavierSlip {
                        velocity: wall_velocity,
                        thermal: WallThermalCondition::Isothermal(Temperature::new(3.0)),
                        slip_length_m: slip,
                        fluid_heat_fraction: fraction,
                    },
                )
                .unwrap();
        prepare(&mut op, &flow);
        let h = 0.5 * spacing[axis];
        let traction = [a, b].map(|t| mu * (wall_velocity[t] - cell_velocity[t]) / (h + slip));
        let trace = [a, b].map(|t| (h * wall_velocity[t] + slip * cell_velocity[t]) / (h + slip));
        let dissipation: f64 = [a, b]
            .into_iter()
            .map(|t| mu * slip * ((wall_velocity[t] - cell_velocity[t]) / (h + slip)).powi(2))
            .sum();
        assert!(dissipation > 0.0);
        let wall_power = traction[0] * wall_velocity[a] + traction[1] * wall_velocity[b];
        let fluid_power = traction[0] * trace[0] + traction[1] * trace[1];
        close(wall_power - fluid_power, dissipation, 2e-14);
        let inward_heat = conductivity * (3.0 - 2.2) / h;
        let expected_inward_energy = fluid_power + fraction * dissipation + inward_heat;
        let mut face = [1; 3];
        face[axis] = if side.sign() > 0.0 { 3 } else { 0 };
        let flux = op.work.flux[axis][grid.cell_face(axis, face)];
        close(flux.normal_momentum, 0.0, 2e-14);
        for t in 0..2 {
            close(
                -side.sign() * flux.transverse_momentum[t],
                traction[t],
                2e-14,
            );
        }
        close(-side.sign() * flux.energy, expected_inward_energy, 2e-14);
        face[axis] = if side.sign() > 0.0 { 2 } else { 0 };
        for t in [a, b] {
            close(
                op.work.gradient[grid.index(face)][t][axis],
                side.sign() * (wall_velocity[t] - cell_velocity[t]) / (2.0 * (h + slip)),
                2e-14,
            );
        }
    }
}

#[test]
fn g1_navier_fluid_trace_retains_the_physical_tangential_normal_stress() {
    let grid = CartesianGrid::new(3, 3, 3, 0.25, 0.25, 0.25).unwrap();
    let (mu, bulk, slip, slope) = (0.02, 0.003, 0.15, 0.05);
    for side in BoundarySide::ALL {
        let axis = side.axis();
        let t = tangents(axis)[0];
        let flow = field(gas(0.0), grid, channel(axis), |x| {
            let distance = if side.sign() > 0.0 {
                0.75 - x[axis]
            } else {
                x[axis]
            };
            let mut u = [0.0; 3];
            u[t] = slope * x[t] * (distance + slip);
            [1.0, u[0], u[1], u[2], 2.0]
        });
        let mut op = MolecularTransport3d::new(&flow, law(flow.gas(), mu, bulk, 0.0), || false)
            .unwrap()
            .with_wall(
                side,
                MolecularWall3d::NavierSlip {
                    velocity: [0.0; 3],
                    thermal: WallThermalCondition::Adiabatic,
                    slip_length_m: slip,
                    fluid_heat_fraction: 1.0,
                },
            )
            .unwrap();
        prepare(&mut op, &flow);
        let mut face = [1; 3];
        face[axis] = if side.sign() > 0.0 { 3 } else { 0 };
        // u_t at the wall is slope*b*x_t, so d_t u_t = slope*b.
        close(
            op.work.flux[axis][grid.cell_face(axis, face)].normal_momentum,
            -(bulk - 2.0 * mu / 3.0) * slope * slip,
            2e-14,
        );
    }
}

#[test]
fn g1_navier_couette_stepping_preserves_linear_velocity_and_heats_the_fluid() {
    let grid = CartesianGrid::new(3, 3, 6, 0.25, 0.25, 0.2).unwrap();
    let (mu, lower_slip, upper_slip, speed, duration) = (0.02, 0.1, 0.2, 0.3, 0.01);
    let slope = speed / (1.2 + lower_slip + upper_slip);
    let mut flow = field(gas(0.0), grid, channel(2), |x| {
        [1.0, slope * (x[2] + lower_slip), 0.0, 0.0, 2.0]
    });
    let original = flow.states().to_vec();
    let mut op = MolecularTransport3d::new(&flow, law(flow.gas(), mu, 0.0, 0.0), || false).unwrap();
    for (side, velocity, slip) in [
        (BoundarySide::ZMin, 0.0, lower_slip),
        (BoundarySide::ZMax, speed, upper_slip),
    ] {
        op = op
            .with_wall(
                side,
                MolecularWall3d::NavierSlip {
                    velocity: [velocity, 0.0, 0.0],
                    thermal: WallThermalCondition::Adiabatic,
                    slip_length_m: slip,
                    // The wall receives interface friction; only resolved fluid
                    // viscous heating mu*(du/dz)^2 remains in this exact solution.
                    fluid_heat_fraction: 0.0,
                },
            )
            .unwrap();
    }
    let report = op
        .advance(
            &mut flow,
            Time::new(duration),
            MolecularControl::default(),
            || false,
        )
        .unwrap();
    for (state, initial) in flow.states().iter().zip(original) {
        for axis in 0..3 {
            close(state.momentum[axis], initial.momentum[axis], 3e-14);
        }
        close(
            state.energy - initial.energy,
            duration * mu * slope * slope,
            3e-14,
        );
    }
    for residual in report.conservation_residual.array() {
        close(residual, 0.0, 3e-13);
    }
    for side in [BoundarySide::ZMin, BoundarySide::ZMax] {
        for impulse in op.boundary_impulses(side) {
            close(
                impulse.momentum[0],
                -side.sign() * mu * slope * grid.face_area(side) * duration,
                3e-14,
            );
        }
    }
}

#[test]
fn g3_navier_heat_partition_changes_energy_without_changing_viscous_momentum() {
    let grid = CartesianGrid::new(3, 3, 4, 0.25, 0.25, 0.25).unwrap();
    let mut results = Vec::new();
    for fraction in [0.0, 0.4, 1.0] {
        let mut flow = field(gas(0.0), grid, channel(2), |_| [1.0, 0.3, -0.1, 0.0, 2.0]);
        let before = flow.totals();
        let (kinetic_before, internal_before) = partition(&flow);
        let mut op =
            MolecularTransport3d::new(&flow, law(flow.gas(), 0.02, 0.0, 0.0), || false).unwrap();
        for side in [BoundarySide::ZMin, BoundarySide::ZMax] {
            op = op
                .with_wall(
                    side,
                    MolecularWall3d::NavierSlip {
                        velocity: [0.0; 3],
                        thermal: WallThermalCondition::Adiabatic,
                        slip_length_m: 0.2,
                        fluid_heat_fraction: fraction,
                    },
                )
                .unwrap();
        }
        let report = op
            .advance(
                &mut flow,
                Time::new(0.02),
                MolecularControl::default(),
                || false,
            )
            .unwrap();
        let (kinetic_after, internal_after) = partition(&flow);
        assert!(kinetic_after < kinetic_before);
        assert!(internal_after > internal_before);
        let outward: f64 = [BoundarySide::ZMin, BoundarySide::ZMax]
            .into_iter()
            .flat_map(|side| op.boundary_impulses(side))
            .map(|impulse| impulse.energy)
            .sum();
        close(flow.totals().energy - before.energy, -outward, 3e-13);
        if fraction == 1.0 {
            assert_eq!(outward, 0.0);
            close(
                internal_after - internal_before,
                kinetic_before - kinetic_after,
                3e-13,
            );
        } else {
            assert!(outward > 0.0);
        }
        for residual in report.conservation_residual.array() {
            close(residual, 0.0, 3e-13);
        }
        results.push((
            flow.states().iter().map(|q| q.momentum).collect::<Vec<_>>(),
            flow.totals().energy,
        ));
    }
    assert_eq!(results[0].0, results[1].0);
    assert_eq!(results[0].0, results[2].0);
    close(results[1].1, 0.6 * results[0].1 + 0.4 * results[2].1, 3e-13);
}

#[test]
fn g0_zero_navier_slip_matches_no_slip_steps_and_refuses_invalid_parameters() {
    let grid = CartesianGrid::new(3, 3, 4, 0.25, 0.25, 0.25).unwrap();
    let mut results = Vec::new();
    for finite_slip in [false, true] {
        let mut flow = field(gas(0.0), grid, channel(2), |x| {
            [1.0, 0.1 * x[2], -0.03, 0.0, 2.0 + 0.1 * x[2]]
        });
        let mut op =
            MolecularTransport3d::new(&flow, law(flow.gas(), 0.02, 0.003, 0.05), || false).unwrap();
        for (side, speed) in [(BoundarySide::ZMin, -0.1), (BoundarySide::ZMax, 0.3)] {
            let velocity = [speed, 0.0, 0.0];
            let thermal = WallThermalCondition::Isothermal(Temperature::new(2.5));
            let wall = if finite_slip {
                MolecularWall3d::NavierSlip {
                    velocity,
                    thermal,
                    slip_length_m: 0.0,
                    fluid_heat_fraction: 0.0,
                }
            } else {
                MolecularWall3d::NoSlip { velocity, thermal }
            };
            op = op.with_wall(side, wall).unwrap();
        }
        let report = op
            .advance(
                &mut flow,
                Time::new(0.002),
                MolecularControl::default(),
                || false,
            )
            .unwrap();
        results.push((
            flow.states()
                .iter()
                .map(|q| q.array().map(f64::to_bits))
                .collect::<Vec<_>>(),
            report,
            op.boundary_impulses(BoundarySide::ZMax).to_vec(),
        ));
    }
    assert_eq!(results[0], results[1]);
    let flow = field(gas(0.0), grid, channel(2), |_| [1.0, 0.1, 0.0, 0.0, 2.0]);
    for (slip, fraction) in [
        (-0.1, 0.5),
        (f64::NAN, 0.5),
        (f64::INFINITY, 0.5),
        (0.1, -0.1),
        (0.1, 1.1),
        (0.1, f64::NAN),
    ] {
        let op =
            MolecularTransport3d::new(&flow, law(flow.gas(), 0.02, 0.0, 0.0), || false).unwrap();
        assert!(matches!(
            op.with_wall(
                BoundarySide::ZMin,
                MolecularWall3d::NavierSlip {
                    velocity: [0.0; 3],
                    thermal: WallThermalCondition::Adiabatic,
                    slip_length_m: slip,
                    fluid_heat_fraction: fraction,
                }
            ),
            Err(MolecularError::Invalid(_))
        ));
    }
}

#[test]
fn g1_no_slip_wall_uses_half_cell_shear_heat_and_prescribed_wall_work() {
    let h = [0.2, 0.3, 0.4];
    let grid = CartesianGrid::new(3, 3, 3, h[0], h[1], h[2]).unwrap();
    let (mu, conductivity) = (0.02, 0.05);
    for axis in 0..3 {
        let [a, b] = tangents(axis);
        let mut cell_velocity = [0.0; 3];
        cell_velocity[a] = 0.15;
        cell_velocity[b] = -0.1;
        let mut wall_velocity = [0.0; 3];
        wall_velocity[a] = 0.4;
        wall_velocity[b] = -0.3;
        let flow = field(gas(0.0), grid, channel(axis), |_| {
            [
                1.0,
                cell_velocity[0],
                cell_velocity[1],
                cell_velocity[2],
                2.2,
            ]
        });
        for upper in [false, true] {
            let side = BoundarySide::ALL[2 * axis + usize::from(upper)];
            let mut op =
                MolecularTransport3d::new(&flow, law(flow.gas(), mu, 0.003, conductivity), || {
                    false
                })
                .unwrap()
                .with_wall(
                    side,
                    MolecularWall3d::NoSlip {
                        velocity: wall_velocity,
                        thermal: WallThermalCondition::Isothermal(Temperature::new(3.0)),
                    },
                )
                .unwrap();
            prepare(&mut op, &flow);
            let sign = side.sign();
            let distance = 0.5 * h[axis];
            let traction =
                [a, b].map(|t| mu * sign * (wall_velocity[t] - cell_velocity[t]) / distance);
            let heat = -conductivity * sign * (3.0 - 2.2) / distance;
            let mut face = [1; 3];
            face[axis] = if upper { 3 } else { 0 };
            let flux = op.work.flux[axis][op.grid.cell_face(axis, face)];
            close(flux.normal_momentum, 0.0, 2e-14);
            for t in 0..2 {
                close(flux.transverse_momentum[t], -traction[t], 2e-14);
            }
            close(
                flux.energy,
                -traction[0] * wall_velocity[a] - traction[1] * wall_velocity[b] + heat,
                2e-14,
            );
            face[axis] = if upper { 2 } else { 0 };
            let cell = op.grid.index(face);
            for t in [a, b] {
                close(
                    op.work.gradient[cell][t][axis],
                    sign * (wall_velocity[t] - cell_velocity[t]) / h[axis],
                    2e-14,
                );
            }
        }
    }
}

#[test]
fn g1_no_slip_adjoint_normal_traction_tends_to_the_constant_boundary_trace() {
    for n in [3, 6] {
        let h = 0.75 / n as f64;
        let grid = CartesianGrid::new(n, n, n, h, h, h).unwrap();
        for axis in 0..3 {
            let t = tangents(axis)[0];
            for upper in [false, true] {
                let flow = field(gas(0.0), grid, channel(axis), |x| {
                    let distance = if upper { 0.75 - x[axis] } else { x[axis] };
                    let mut velocity = [0.0; 3];
                    velocity[t] = 0.05 * x[t] * distance;
                    [1.0, velocity[0], velocity[1], velocity[2], 2.0]
                });
                let side = BoundarySide::ALL[2 * axis + usize::from(upper)];
                let mut op =
                    MolecularTransport3d::new(&flow, law(flow.gas(), 0.02, 0.003, 0.01), || false)
                        .unwrap()
                        .with_wall(
                            side,
                            MolecularWall3d::NoSlip {
                                velocity: [0.0; 3],
                                thermal: WallThermalCondition::Adiabatic,
                            },
                        )
                        .unwrap();
                prepare(&mut op, &flow);
                let mut face = [1; 3];
                face[axis] = if upper { n } else { 0 };
                // The adjacent centered transverse strain is required by the
                // mechanical adjoint identity. For this smooth zero wall trace,
                // its boundary traction tends to zero linearly with h; this is
                // not a claim of an exact physical trace on a finite grid.
                close(
                    op.work.flux[axis][grid.cell_face(axis, face)].normal_momentum,
                    -(0.003 - 2.0 * 0.02 / 3.0) * 0.05 * (0.5 * h),
                    2e-14,
                );
                face[axis] = if upper { n - 1 } else { 0 };
                close(
                    op.work.gradient[grid.index(face)][t][t],
                    0.05 * (0.5 * h),
                    2e-14,
                );
            }
        }
    }
}

#[test]
fn g0_coupled_co_moving_walls_preserve_uniform_translation_on_every_axis() {
    for axis in 0..3 {
        let [a, b] = tangents(axis);
        let mut velocity = [0.0; 3];
        velocity[a] = 0.25;
        velocity[b] = -0.125;
        let grid = CartesianGrid::new(3, 3, 3, 0.2, 0.3, 0.4).unwrap();
        let flow = field(gas(0.0), grid, channel(axis), |_| {
            [1.0, velocity[0], velocity[1], velocity[2], 2.0]
        });
        let original = flow.states().to_vec();
        let source = law(flow.gas(), 0.02, 0.003, 0.05);
        let mut coupled = MolecularEuler3d::new(flow, source, || false).unwrap();
        for side in &BoundarySide::ALL[2 * axis..2 * axis + 2] {
            coupled = coupled
                .with_wall(
                    *side,
                    MolecularWall3d::NoSlip {
                        velocity,
                        thermal: WallThermalCondition::Adiabatic,
                    },
                )
                .unwrap();
        }
        for _ in 0..3 {
            let report = coupled
                .step(Time::new(0.001), MolecularStepControl::default(), || false)
                .unwrap();
            for (current, initial) in coupled.flow().states().iter().zip(&original) {
                assert_eq!(
                    current.array().map(f64::to_bits),
                    initial.array().map(f64::to_bits)
                );
            }
            for residual in report.conservation_residual.array() {
                close(residual, 0.0, 3e-14);
            }
        }
    }
}

#[test]
fn g3_moving_heated_walls_account_for_actual_total_energy_and_momentum() {
    for heated in [false, true] {
        let grid = CartesianGrid::new(3, 3, 4, 0.25, 0.25, 0.25).unwrap();
        let mut flow = field(gas(0.0), grid, channel(2), |_| [1.0, 0.1, 0.0, 0.0, 2.2]);
        let original = flow.totals();
        let mut op =
            MolecularTransport3d::new(&flow, law(flow.gas(), 0.02, 0.003, 0.05), || false).unwrap();
        for (side, speed, temperature) in [
            (BoundarySide::ZMin, 0.0, 2.0),
            (BoundarySide::ZMax, 0.5, 3.0),
        ] {
            op = op
                .with_wall(
                    side,
                    MolecularWall3d::NoSlip {
                        velocity: [speed, 0.0, 0.0],
                        thermal: if heated {
                            WallThermalCondition::Isothermal(Temperature::new(temperature))
                        } else {
                            WallThermalCondition::Adiabatic
                        },
                    },
                )
                .unwrap();
        }
        let report = op
            .advance(
                &mut flow,
                Time::new(0.002),
                MolecularControl::default(),
                || false,
            )
            .unwrap();
        assert!(flow.totals().energy > original.energy);
        assert!(flow.totals().momentum[0] > original.momentum[0]);
        close(flow.totals().mass, original.mass, 2e-14);
        for residual in report.conservation_residual.array() {
            close(residual, 0.0, 2e-13);
        }
        let mut outward_energy = 0.0;
        for (side, speed) in [(BoundarySide::ZMin, 0.0), (BoundarySide::ZMax, 0.5)] {
            for impulse in op.boundary_impulses(side) {
                outward_energy += impulse.energy;
                if !heated {
                    close(impulse.energy, speed * impulse.momentum[0], 3e-14);
                }
            }
        }
        close(
            flow.totals().energy - original.energy,
            -outward_energy,
            2e-13,
        );
    }
}

#[cfg(feature = "experimental-sgs")]
#[test]
fn g1_sgs_adds_interior_transport_but_preserves_molecular_wall_flux() {
    use crate::compressible::spatial::molecular::sgs::{SgsModel, SgsUse, SubgridConfig};

    let grid = CartesianGrid::new(4, 4, 4, 0.25, 0.25, 0.25).unwrap();
    let flow = field(gas(0.0), grid, channel(2), |x| {
        [1.0, 0.1 * x[2], 0.0, 0.0, 2.0 + 0.2 * x[2]]
    });
    let mut base =
        MolecularTransport3d::new(&flow, law(flow.gas(), 0.002, 0.003, 0.004), || false).unwrap();
    let config = SubgridConfig {
        model: SgsModel::Smagorinsky { cs: 0.17 },
        filter_width_factor: 1.0,
        turbulent_prandtl: 0.8,
        maximum_eddy_to_molecular: 20.0,
        usage: SgsUse::Exploration,
    };
    let mut active =
        MolecularTransport3d::with_sgs(&flow, law(flow.gas(), 0.002, 0.003, 0.004), config, || {
            false
        })
        .unwrap();
    for (side, speed, temperature) in [
        (BoundarySide::ZMin, 0.0, 2.0),
        (BoundarySide::ZMax, 0.1, 2.2),
    ] {
        let wall = MolecularWall3d::NoSlip {
            velocity: [speed, 0.0, 0.0],
            thermal: WallThermalCondition::Isothermal(Temperature::new(temperature)),
        };
        base = base.with_wall(side, wall).unwrap();
        active = active.with_wall(side, wall).unwrap();
    }
    prepare(&mut base, &flow);
    prepare(&mut active, &flow);
    for face in [0, 4] {
        let index = grid.cell_face(2, [1, 1, face]);
        assert_eq!(base.work.flux[2][index], active.work.flux[2][index]);
    }
    let index = grid.cell_face(2, [1, 1, 2]);
    assert!(
        active.work.flux[2][index].transverse_momentum[0]
            < base.work.flux[2][index].transverse_momentum[0]
    );
    assert!(active.work.flux[2][index].energy < base.work.flux[2][index].energy);
}
