use super::*;
use crate::compressible::tests::{close, gas};
use fs_qty::{Acceleration, Length, Pressure, Temperature};

const WALLS: AxisBoundaries = AxisBoundaries::Ends {
    lower: BoundaryCondition::Reflecting,
    upper: BoundaryCondition::Reflecting,
};

fn reference(family: usize, g: f64) -> HydrostaticReference {
    match family {
        0 => HydrostaticReference::Isothermal(IsothermalAtmosphere {
            acceleration: Acceleration::new(g),
            temperature: Temperature::new(2.0),
            reference_pressure: Pressure::new(3.0),
            reference_position: Length::new(0.0),
        }),
        1 => HydrostaticReference::Isentropic(IsentropicAtmosphere {
            acceleration: Acceleration::new(g),
            reference_temperature: Temperature::new(2.0),
            reference_pressure: Pressure::new(3.0),
            reference_position: Length::new(0.0),
        }),
        _ => HydrostaticReference::LapseRate(LapseRateAtmosphere {
            acceleration: Acceleration::new(g),
            reference_temperature: Temperature::new(2.0),
            reference_pressure: Pressure::new(3.0),
            reference_position: Length::new(0.0),
            temperature_gradient_k_m: 0.4 * g.signum(),
        }),
    }
}

fn geometry(axis: usize, ends: AxisBoundaries) -> (CartesianGrid, Boundaries3d) {
    let mut shape = [0; 3];
    shape[axis] = 8;
    shape[(axis + 1) % 3] = 3;
    shape[(axis + 2) % 3] = 4;
    let spacing = shape.map(|n| 1.0 / n as f64);
    let grid = CartesianGrid::new(
        shape[0], shape[1], shape[2], spacing[0], spacing[1], spacing[2],
    )
    .unwrap();
    let mut pairs = [AxisBoundaries::Periodic; 3];
    pairs[axis] = ends;
    (
        grid,
        Boundaries3d {
            x: pairs[0],
            y: pairs[1],
            z: pairs[2],
        },
    )
}

fn atmosphere(axis: usize, reference: HydrostaticReference, ends: AxisBoundaries) -> Euler3d {
    let (grid, boundaries) = geometry(axis, ends);
    Euler3d::new_hydrostatic_atmosphere(
        gas(0.4),
        grid,
        axis,
        reference,
        AtmosphereQuadrature::default(),
        boundaries,
        || false,
    )
    .unwrap()
}

fn balance(report: &StepReport) {
    assert!(report.maximum_cfl <= 0.40000000000001);
    let scale = 1.0 + report.before.mass.abs() + report.before.energy.abs();
    for residual in report.conservation_residual.array() {
        assert!(residual.abs() < 3e-12 * scale, "{report:?}");
    }
}

/// Independent analytic endpoint for the manufactured R=1, Cp=3.5+0.4*T gas.
fn endpoint_pressure(family: usize, g: f64) -> f64 {
    match family {
        0 => 3.0 * (g / 2.0).exp(),
        1 => {
            let cp0 = 4.3_f64;
            let dt = 2.0 * g / (cp0 + (cp0 * cp0 + 0.8 * g).sqrt());
            3.0 * (3.5 * (1.0 + dt / 2.0).ln() + 0.4 * dt).exp()
        }
        _ => {
            let lapse = 0.4 * g.signum();
            3.0 * ((g / lapse) * (1.0 + lapse / 2.0).ln()).exp()
        }
    }
}

#[test]
fn integrated_atmospheres_remain_at_rest_on_every_axis_sign_and_order() {
    for axis in 0..3 {
        for family in 0..3 {
            for sign in [-1.0, 1.0] {
                let g = sign * if family == 0 { 12.0 } else { 3.0 };
                for reconstruction in [Reconstruction::FirstOrder, Reconstruction::Muscl] {
                    let declaration = reference(family, g);
                    let mut flow = atmosphere(axis, declaration, WALLS);
                    assert_eq!(flow.hydrostatic_reference(), Some((axis, declaration)));
                    // Unit transverse area: integral rho dx = (p(1)-p(0))/g.
                    // This tests conserved averages independently of source quadrature.
                    let expected_mass = (endpoint_pressure(family, g) - 3.0) / g;
                    close(flow.totals().mass, expected_mass, 3e-10);
                    let initial = flow.states().to_vec();
                    for _ in 0..8 {
                        let report = flow
                            .step(
                                Time::new(0.02),
                                StepControl {
                                    reconstruction,
                                    ..StepControl::default()
                                },
                                || false,
                            )
                            .unwrap();
                        assert_eq!(report.halvings, 0);
                        assert_eq!(report.source_change.work, 0.0);
                        close(
                            report.source_change.impulse[axis],
                            expected_mass * g * report.elapsed.value(),
                            5e-10,
                        );
                        balance(&report);
                    }
                    assert_eq!(
                        flow.states(),
                        initial,
                        "axis={axis}, family={family}, g={g}"
                    );
                }
            }
        }
    }
}

#[test]
fn three_dimensional_pressure_disturbance_evolves_with_conservative_gravity_work() {
    let reference = reference(0, -8.0);
    let equilibrium = atmosphere(2, reference, WALLS);
    let states = equilibrium
        .states()
        .iter()
        .enumerate()
        .map(|(cell, &q)| {
            let x = equilibrium.grid.cell_center(cell).unwrap();
            if x.into_iter().all(|v| (0.25..0.75).contains(&v)) {
                let p = q.primitive(equilibrium.gas(), || false).unwrap();
                State::from_primitive(
                    equilibrium.gas(),
                    1.01 * q.density,
                    [0.0; 3],
                    1.04 * p.pressure,
                )
                .unwrap()
            } else {
                q
            }
        })
        .collect();
    let mut flow = Euler3d::new(
        equilibrium.gas.clone(),
        equilibrium.grid,
        states,
        equilibrium.boundaries,
        || false,
    )
    .unwrap();
    let initial = flow.states().to_vec();
    flow.bind_hydrostatic_reference(2, reference, AtmosphereQuadrature::default(), || false)
        .unwrap();
    assert_eq!(flow.states(), initial);
    let mass = flow.totals().mass;
    let gas_energy = flow.totals().energy;
    let total_energy = gas_energy + flow.gravitational_potential_energy();
    let mut gravity_work = 0.0;
    for _ in 0..32 {
        let report = flow
            .step(Time::new(0.01), StepControl::default(), || false)
            .unwrap();
        gravity_work += report.source_change.work;
        balance(&report);
    }
    for axis in 0..3 {
        assert!(flow.states().iter().any(|q| q.momentum[axis].abs() > 1e-5));
    }
    assert!(
        flow.states()
            .iter()
            .zip(&initial)
            .any(|(new, old)| (new.density - old.density).abs() > 1e-5)
    );
    assert!(
        flow.states()
            .iter()
            .any(|q| (q.primitive(flow.gas(), || false).unwrap().temperature - 2.0).abs() > 1e-4)
    );
    close(flow.totals().mass, mass, 3e-13);
    close(
        flow.totals().energy + flow.gravitational_potential_energy(),
        total_energy,
        3e-12,
    );
    close(flow.totals().energy - gas_energy, gravity_work, 3e-12);
}

#[test]
fn both_transverse_momenta_and_kinetic_energy_survive_axis_permutation() {
    let reference = reference(1, -3.0);
    let mut canonical = None;
    for axis in 0..3 {
        let equilibrium = atmosphere(axis, reference, WALLS);
        let a = (axis + 1) % 3;
        let b = (axis + 2) % 3;
        let mut velocity = [0.0; 3];
        velocity[a] = 0.12;
        velocity[b] = -0.07;
        let speed_squared: f64 = velocity.iter().map(|u| u * u).sum();
        let states = equilibrium
            .states()
            .iter()
            .map(|&q| State {
                momentum: velocity.map(|u| q.density * u),
                energy: q.energy + 0.5 * q.density * speed_squared,
                ..q
            })
            .collect();
        let mut flow = Euler3d::new(
            equilibrium.gas.clone(),
            equilibrium.grid,
            states,
            equilibrium.boundaries,
            || false,
        )
        .unwrap();
        flow.bind_hydrostatic_reference(axis, reference, AtmosphereQuadrature::default(), || false)
            .unwrap();
        let initial = flow.states().to_vec();
        for _ in 0..8 {
            let report = flow
                .step(Time::new(0.005), StepControl::default(), || false)
                .unwrap();
            balance(&report);
        }
        for (q, before) in flow.states().iter().zip(&initial) {
            close(q.density, before.density, 3e-11);
            close(q.energy, before.energy, 3e-11);
            for component in 0..3 {
                close(
                    q.momentum[component] / q.density,
                    velocity[component],
                    3e-11,
                );
            }
        }
        let totals = flow.totals();
        let permuted = [
            totals.mass,
            totals.energy,
            totals.momentum[a],
            totals.momentum[b],
        ];
        if let Some(expected) = canonical {
            for (actual, expected) in permuted.into_iter().zip(expected) {
                close(actual, expected, 3e-12);
            }
        } else {
            canonical = Some(permuted);
        }
    }
}

#[test]
fn open_faces_transport_potential_with_signed_mass_and_a_nonzero_datum() {
    let (axis, g, datum) = (1, -4.0, 0.3);
    let declaration = HydrostaticReference::Isothermal(IsothermalAtmosphere {
        acceleration: Acceleration::new(g),
        temperature: Temperature::new(2.0),
        reference_pressure: Pressure::new(3.0),
        reference_position: Length::new(datum),
    });
    let equilibrium = atmosphere(
        axis,
        declaration,
        AxisBoundaries::Ends {
            lower: BoundaryCondition::Transmissive,
            upper: BoundaryCondition::Transmissive,
        },
    );
    let velocity = [0.04, 0.12, -0.03];
    let speed_squared: f64 = velocity.iter().map(|u| u * u).sum();
    let states = equilibrium
        .states()
        .iter()
        .map(|&q| State {
            momentum: velocity.map(|u| q.density * u),
            energy: q.energy + 0.5 * q.density * speed_squared,
            ..q
        })
        .collect();
    let mut flow = Euler3d::new(
        equilibrium.gas.clone(),
        equilibrium.grid,
        states,
        equilibrium.boundaries,
        || false,
    )
    .unwrap();
    flow.bind_hydrostatic_reference(axis, declaration, AtmosphereQuadrature::default(), || false)
        .unwrap();
    let mass0 = flow.totals().mass;
    let energy0 = flow.totals().energy + flow.gravitational_potential_energy();
    let mut mass_transfer = 0.0;
    let mut energy_transfer = 0.0;
    let mut outward_mass = [0.0; 2];
    for _ in 0..24 {
        let report = flow
            .step(Time::new(0.01), StepControl::default(), || false)
            .unwrap();
        let mut potential_out = 0.0;
        for end in 0..2 {
            let side = BoundarySide::ALL[2 * axis + end];
            for (face, load) in flow.boundary_impulses(side).iter().enumerate() {
                let x = flow.grid.face_center(side, face).unwrap()[axis];
                potential_out += -g * (x - datum) * load.mass;
                outward_mass[end] += load.mass;
            }
        }
        mass_transfer += report.boundary_change.mass;
        energy_transfer += report.boundary_change.energy - potential_out;
        balance(&report);
    }
    assert!(outward_mass[0] < -1e-4 && outward_mass[1] > 1e-4);
    assert!(mass_transfer.abs() > 1e-4);
    close(flow.totals().mass - mass0, mass_transfer, 3e-12);
    close(
        flow.totals().energy + flow.gravitational_potential_energy() - energy0,
        energy_transfer,
        3e-12,
    );
}

#[test]
fn pressure_outlet_uses_the_hydrostatic_face_instead_of_the_last_cell_average() {
    for axis in 0..3 {
        for g in [-6.0, 6.0] {
            for sign in [-1.0, 1.0] {
                let declaration = reference(0, g);
                let equilibrium = atmosphere(axis, declaration, WALLS);
                let mut c = [0; 3];
                c[axis] = if sign > 0.0 {
                    equilibrium.grid.shape()[axis] - 1
                } else {
                    0
                };
                let last = equilibrium.grid.index(c);
                let profile = equilibrium.hydrostatic.as_ref().unwrap();
                let face_pressure = profile.trace(&equilibrium.grid, last, axis, sign)[4];
                let cell_pressure = profile.cell(&equilibrium.grid, last).1.pressure;
                assert!((cell_pressure - face_pressure).abs() > 1e-3 * face_pressure);
                close(
                    face_pressure,
                    if sign > 0.0 {
                        3.0 * (g / 2.0).exp()
                    } else {
                        3.0
                    },
                    3e-13,
                );
                let ends = AxisBoundaries::Ends {
                    lower: if sign < 0.0 {
                        BoundaryCondition::PressureOutlet {
                            pressure_pa: face_pressure,
                        }
                    } else {
                        BoundaryCondition::Reflecting
                    },
                    upper: if sign > 0.0 {
                        BoundaryCondition::PressureOutlet {
                            pressure_pa: face_pressure,
                        }
                    } else {
                        BoundaryCondition::Reflecting
                    },
                };
                for reconstruction in [Reconstruction::FirstOrder, Reconstruction::Muscl] {
                    let mut flow = atmosphere(axis, declaration, ends);
                    let initial = flow.states().to_vec();
                    for _ in 0..4 {
                        let report = flow
                            .step(
                                Time::new(0.01),
                                StepControl {
                                    reconstruction,
                                    ..StepControl::default()
                                },
                                || false,
                            )
                            .unwrap();
                        balance(&report);
                        close(report.boundary_change.mass, 0.0, 3e-14);
                        close(report.boundary_change.energy, 0.0, 3e-14);
                    }
                    assert_eq!(flow.states(), initial);
                }
            }
        }
    }
}

#[test]
fn cancelled_or_invalid_binding_preserves_the_original_field_and_reference() {
    let declaration = reference(1, -3.0);
    let equilibrium = atmosphere(2, declaration, WALLS);
    let initial = equilibrium.states().to_vec();
    let mut flow = Euler3d::new(
        equilibrium.gas.clone(),
        equilibrium.grid,
        initial.clone(),
        equilibrium.boundaries,
        || false,
    )
    .unwrap();
    // Gravity cannot be periodic, and an invalid axis is refused before indexing.
    for axis in [0, 3] {
        assert!(
            flow.bind_hydrostatic_reference(
                axis,
                declaration,
                AtmosphereQuadrature::default(),
                || false
            )
            .is_err()
        );
        assert_eq!(flow.states(), initial);
        assert_eq!(flow.hydrostatic_reference(), None);
    }
    let mut polls = 0;
    assert_eq!(
        flow.bind_hydrostatic_reference(2, declaration, AtmosphereQuadrature::default(), || {
            polls += 1;
            polls == 20
        }),
        Err(Error::Cancelled)
    );
    assert_eq!(flow.states(), initial);
    assert_eq!(flow.hydrostatic_reference(), None);
    assert_eq!(flow.time(), Time::new(0.0));
    assert_eq!(flow.last_step(), None);
    flow.bind_hydrostatic_reference(2, declaration, AtmosphereQuadrature::default(), || false)
        .unwrap();
    assert!(
        flow.bind_hydrostatic_reference(
            2,
            reference(1, 3.0),
            AtmosphereQuadrature::default(),
            || false
        )
        .is_err()
    );
    assert_eq!(flow.hydrostatic_reference(), Some((2, declaration)));
    assert_eq!(flow.states(), initial);
    let report = flow
        .step(Time::new(0.005), StepControl::default(), || false)
        .unwrap();
    let accepted = flow.states().to_vec();
    assert_eq!(
        flow.step(Time::new(0.005), StepControl::default(), || true),
        Err(Error::Cancelled)
    );
    assert_eq!(flow.last_step(), Some(&report));
    assert_eq!(flow.time(), report.elapsed);
    assert_eq!(flow.states(), accepted);
    assert_eq!(flow.hydrostatic_reference(), Some((2, declaration)));
}

#[test]
fn forced_entry_keeps_gravity_balance_and_both_transverse_accelerations() {
    let mut flow = atmosphere(2, reference(0, -6.0), WALLS);
    let mass = flow.totals().mass;
    let mut workspace = flow.forcing_workspace(usize::MAX).unwrap();
    let report = flow
        .step_with_forcing(
            Time::new(0.005),
            StepControl::default(),
            BodyForce3d::UniformAcceleration([0.7, -0.2, 0.0]),
            ForceControl3d { maximum_mach: 1.0 },
            &mut workspace,
            || false,
        )
        .unwrap();
    let dt = report.elapsed.value();
    for q in flow.states() {
        close(q.momentum[0] / q.density, 0.7 * dt, 3e-10);
        close(q.momentum[1] / q.density, -0.2 * dt, 3e-10);
        close(q.momentum[2] / q.density, 0.0, 3e-10);
    }
    for (axis, acceleration) in [0.7, -0.2, -6.0].into_iter().enumerate() {
        close(
            report.source_change.impulse[axis],
            mass * acceleration * dt,
            3e-12,
        );
    }
    close(
        report.source_change.work,
        0.5 * mass * (0.7_f64.powi(2) + 0.2_f64.powi(2)) * dt * dt,
        3e-12,
    );
    balance(&report);
}

#[test]
fn molecular_trial_preserves_its_original_hydrostatic_reference_after_swaps() {
    use crate::compressible::spatial::molecular::{MolecularEuler3d, MolecularStepControl};
    use fs_gas::multicomponent::dissipation::{
        ConstantTransport, TransportCoefficients, TransportDomain,
    };
    use fs_qty::DynViscosity;
    use fs_thermochem::equilibrium::composition::IdealGasEquilibriumV1;
    use fs_thermochem::transport::kinetic::ThermalConductivityQuantityV1;
    use fs_thermochem::{ChargeVector, ElementId, ElementalMatrix};

    let declaration = reference(0, -6.0);
    let flow = atmosphere(2, declaration, WALLS);
    let initial = flow.states().to_vec();
    let models = flow.gas().source().components().to_vec();
    let ids: Vec<_> = models.iter().map(|m| m.species().clone()).collect();
    let source = IdealGasEquilibriumV1::new(
        ElementalMatrix::new(
            vec![ElementId::new("H").unwrap()],
            ids.clone(),
            vec![vec![1; ids.len()]],
        )
        .unwrap(),
        ChargeVector::new(ids.clone(), vec![0; ids.len()]).unwrap(),
        models,
    )
    .unwrap();
    let law = ConstantTransport::new(
        source,
        "manufactured-spatial-transport",
        1,
        TransportDomain {
            minimum_temperature_k: 0.001,
            maximum_temperature_k: 50.0,
            minimum_pressure_pa: 1e-8,
            maximum_pressure_pa: 1e8,
        },
        TransportCoefficients {
            dynamic_viscosity: DynViscosity::new(0.01),
            bulk_viscosity: DynViscosity::new(0.0),
            thermal_conductivity: ThermalConductivityQuantityV1::new(0.02),
        },
        vec![0.0],
    )
    .unwrap();
    let mut coupled = MolecularEuler3d::new(flow, law, || false).unwrap();
    // Both parity states must retain gravity when accepted/trial owners swap.
    for _ in 0..4 {
        let report = coupled
            .step(Time::new(0.002), MolecularStepControl::default(), || false)
            .unwrap();
        balance(&report.advection);
        assert_eq!(
            coupled.flow().hydrostatic_reference(),
            Some((2, declaration))
        );
        for (q, before) in coupled.flow().states().iter().zip(&initial) {
            close(q.density, before.density, 5e-11);
            close(q.energy, before.energy, 5e-11);
            for momentum in q.momentum {
                close(momentum / q.density, 0.0, 5e-11);
            }
        }
    }
}
