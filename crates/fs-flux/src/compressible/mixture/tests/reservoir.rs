use super::molecular::{law, prepare};
use super::*;
use crate::compressible::mixture::molecular::{
    DissipationControl, DissipationError, MolecularTransport2d, PlanarTransport,
};

fn equal_mass() -> MixtureGas {
    gas(vec![
        source("A", R, 3.5, 0.0, 40.0),
        source("B", R, 4.5, 0.0, -50.0),
        source("absent-C", R, 3.5, 0.0, 2.0),
    ])
}

pub(super) fn bind(flow: MixtureEuler2d, inputs: &[ReservoirBoundary<'_>]) -> MixtureEuler2d {
    let bytes = MixtureEuler2d::required_workspace_bytes_with_reservoirs(
        flow.gas(),
        flow.grid(),
        inputs.len(),
    )
    .unwrap();
    let flow = MixtureEuler2d::new_with_reservoirs(
        flow.gas,
        flow.grid,
        flow.states,
        flow.boundaries,
        inputs,
        bytes,
        || false,
    )
    .unwrap();
    assert_eq!(flow.workspace_bytes(), bytes);
    flow
}

fn sides(axis: usize) -> [BoundarySide; 2] {
    if axis == 0 {
        [BoundarySide::XMin, BoundarySide::XMax]
    } else {
        [BoundarySide::YMin, BoundarySide::YMax]
    }
}

#[test]
fn g1_held_source_reservoirs_supply_true_inflow_and_outflow_in_both_axes() {
    for axis in 0..2 {
        for sign in [-1.0, 1.0] {
            let grid = CartesianGrid::new(1, 1, 0.4, 0.7, 0.2).unwrap();
            let topology = if axis == 0 {
                Boundaries2d {
                    x: OPEN,
                    y: AxisBoundaries::Periodic,
                }
            } else {
                Boundaries2d {
                    x: AxisBoundaries::Periodic,
                    y: OPEN,
                }
            };
            let mut velocity = [0.0; 2];
            velocity[axis] = 3.0 * sign;
            let side = sides(axis);
            let inputs = [
                ReservoirBoundary {
                    side: side[0],
                    partial_densities: &[0.8, 0.2, 0.0],
                    velocity,
                    pressure_pa: 1.0,
                },
                ReservoirBoundary {
                    side: side[1],
                    partial_densities: &[0.1, 0.9, 0.0],
                    velocity,
                    pressure_pa: 1.0,
                },
            ];
            let mut flow = bind(
                field(equal_mass(), grid, topology, |_| {
                    ([0.3, 0.7, 0.0], velocity, 1.0)
                }),
                &inputs,
            );
            let dt = 1e-3;
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
            balance(report);
            // Supersonic same-p/u contact: incoming composition is exactly the
            // upstream bath; the other bath cannot alter outgoing material.
            // One finite cell obeys da/dt=(|u|/dx)*(a_in-a), with SSP RK2.
            let incoming = if sign > 0.0 { 0.8 } else { 0.1 };
            let courant = 3.0 * dt / grid.spacing[axis];
            let expected = 0.3 + (incoming - 0.3) * (courant - 0.5 * courant * courant);
            let q = flow.cell(0).unwrap();
            close(q.partial_densities[0], expected, 3e-12);
            close(q.partial_densities[1], 1.0 - expected, 3e-12);
            assert_eq!(q.partial_densities[2], 0.0);
            let p = q.primitive(flow.gas(), || false).unwrap();
            close(p.temperature, 1.0, 2e-11);
            close(p.pressure, 1.0, 2e-11);
            for d in 0..2 {
                close(p.velocity[d], velocity[d], 2e-12);
            }
            for (index, side) in side.into_iter().enumerate() {
                let exchange = flow.boundary_impulses(side).next().unwrap();
                let area = grid.spacing[1 - axis] * grid.depth;
                let outward = if index == 0 { -1.0 } else { 1.0 };
                close(
                    exchange.species_masses.iter().sum::<f64>(),
                    outward * velocity[axis] * area * dt,
                    3e-12,
                );
                assert_eq!(flow.reservoir(side).unwrap(), inputs[index]);
            }
        }
    }
}

fn diffusion_error(n: usize) -> f64 {
    let dx = 1.0 / n as f64;
    let grid = CartesianGrid::new(n, 3, dx, 1.0 / 3.0, 0.6).unwrap();
    let topology = Boundaries2d {
        x: OPEN,
        y: AxisBoundaries::Periodic,
    };
    let inputs = [
        ReservoirBoundary {
            side: BoundarySide::XMin,
            partial_densities: &[0.2, 0.8, 0.0],
            velocity: [0.0; 2],
            pressure_pa: 1.0,
        },
        ReservoirBoundary {
            side: BoundarySide::XMax,
            partial_densities: &[0.6, 0.4, 0.0],
            velocity: [0.0; 2],
            pressure_pa: 1.0,
        },
    ];
    let z = 0.5 * std::f64::consts::PI * dx;
    let average = z.sin() / z;
    let mut flow = bind(
        field(equal_mass(), grid, topology, |x| {
            let a = 0.2 + 0.4 * x[0] + 0.1 * average * (std::f64::consts::PI * x[0]).sin();
            ([a, 1.0 - a, 0.0], [0.0; 2], 1.0)
        }),
        &inputs,
    );
    let diffusion = 0.03;
    let duration = 0.1;
    let mut operator = prepare(&flow, law(flow.gas(), diffusion, 0.0, 0.0, 0.0));
    let report = operator
        .advance(
            &mut flow,
            Time::new(duration),
            DissipationControl::default(),
            || false,
        )
        .unwrap();
    assert_eq!(report.elapsed.value(), duration);
    assert!(report.boundary_change.species_masses[0] < 0.0);
    close(report.boundary_change.mass(), 0.0, 2e-12);
    close(report.absolute_energy_residual_j, 0.0, 2e-10);
    for residual in report.conservation_residual.species_masses {
        close(residual, 0.0, 2e-12);
    }
    let decay = (-diffusion * std::f64::consts::PI.powi(2) * duration).exp();
    let mut error = 0.0;
    for index in 0..flow.cell_count() {
        let x = grid.cell_center(index).unwrap()[0];
        let expected = 0.2 + 0.4 * x + 0.1 * average * decay * (std::f64::consts::PI * x).sin();
        let q = flow.cell(index).unwrap();
        error += (q.partial_densities[0] - expected).abs();
        assert_eq!(q.partial_densities[2], 0.0);
        let p = q.primitive(flow.gas(), || false).unwrap();
        close(p.temperature, 1.0, 2e-10);
        close(p.pressure, 1.0, 2e-10);
    }
    assert_eq!(flow.time().value(), 0.0);
    error / flow.cell_count() as f64
}

#[test]
fn g1_boundary_species_exchange_carries_original_enthalpy_and_refines_to_dirichlet_solution() {
    let coarse = diffusion_error(8);
    let fine = diffusion_error(16);
    assert!(fine < 0.3 * coarse, "coarse={coarse:e}, fine={fine:e}");
    assert!(fine < 1e-5);
}

#[test]
fn g1_reservoir_heat_and_shear_work_use_actual_boundary_efforts() {
    let topology = Boundaries2d {
        x: OPEN,
        y: AxisBoundaries::Periodic,
    };
    let grid = CartesianGrid::new(1, 1, 0.4, 0.7, 0.2).unwrap();
    let hot = [
        ReservoirBoundary {
            side: BoundarySide::XMin,
            partial_densities: &[1.0, 0.0, 0.0],
            velocity: [0.0; 2],
            pressure_pa: 3.0,
        },
        ReservoirBoundary {
            side: BoundarySide::XMax,
            partial_densities: &[1.0, 0.0, 0.0],
            velocity: [0.0; 2],
            pressure_pa: 3.0,
        },
    ];
    let mut flow = bind(
        field(equal_mass(), grid, topology, |_| {
            ([1.0, 0.0, 0.0], [0.0; 2], 2.0)
        }),
        &hot,
    );
    let conductivity = 0.04;
    let dt = 1e-3;
    // Binary diffusivities remain positive in the original transport law.
    // Pure A has x_A=y_A=1 everywhere, so composition and barodiffusion
    // driving forces vanish exactly even across the thermal pressure jump.
    let mut operator = prepare(&flow, law(flow.gas(), 0.02, 0.0, 0.0, conductivity));
    let report = operator
        .advance(
            &mut flow,
            Time::new(dt),
            DissipationControl::default(),
            || false,
        )
        .unwrap();
    let rate = 4.0 * conductivity / (2.5 * grid.spacing[0].powi(2));
    let expected = 3.0 - (1.0 - rate * dt + 0.5 * (rate * dt).powi(2));
    close(
        flow.cell(0)
            .unwrap()
            .primitive(flow.gas(), || false)
            .unwrap()
            .temperature,
        expected,
        3e-12,
    );
    close(
        report.after.energy - report.before.energy,
        2.5 * grid.cell_volume() * (expected - 2.0),
        3e-12,
    );
    close(report.absolute_energy_residual_j, 0.0, 3e-12);
    for side in sides(0) {
        let impulse = operator.boundary_impulses(side).next().unwrap();
        assert!(impulse.energy < 0.0);
        assert_eq!(impulse.species_masses, [0.0; 3]);
    }

    // Linear Couette field: full stress has constant tau_yx=mu*V/L,
    // so velocities stay fixed. A stationary bath supplies exactly zero
    // traction work; the moving bath supplies mu*V²*A/L, independent of dx.
    for n in [2, 5] {
        let grid = CartesianGrid::new(n, 2, 1.0 / n as f64, 0.5, 0.3).unwrap();
        let speed = 0.4;
        let reservoirs = [
            hot[0],
            ReservoirBoundary {
                velocity: [0.0, speed],
                pressure_pa: 2.0,
                ..hot[1]
            },
        ];
        let reservoirs = [
            ReservoirBoundary {
                pressure_pa: 2.0,
                ..reservoirs[0]
            },
            reservoirs[1],
        ];
        let mut shear = bind(
            field(equal_mass(), grid, topology, |x| {
                ([1.0, 0.0, 0.0], [0.0, speed * x[0]], 2.0)
            }),
            &reservoirs,
        );
        let mu = 0.03;
        let mut operator = prepare(&shear, law(shear.gas(), 0.02, mu, 0.01, 0.0));
        let report = operator
            .advance(
                &mut shear,
                Time::new(dt),
                DissipationControl::default(),
                || false,
            )
            .unwrap();
        close(
            report.after.energy - report.before.energy,
            mu * speed * speed * grid.depth * dt,
            3e-12,
        );
        close(report.absolute_energy_residual_j, 0.0, 3e-12);
        for impulse in operator.boundary_impulses(BoundarySide::XMin) {
            assert_eq!(impulse.energy, 0.0);
        }
        for index in 0..shear.cell_count() {
            let p = shear
                .cell(index)
                .unwrap()
                .primitive(shear.gas(), || false)
                .unwrap();
            close(
                p.velocity[1],
                speed * grid.cell_center(index).unwrap()[0],
                2e-12,
            );
            close(p.temperature, 2.0 + mu * speed * speed * dt / 2.5, 2e-11);
        }
    }
}

#[test]
fn g4_reservoir_admission_budget_and_late_refusal_preserve_fields_and_loads() {
    let grid = CartesianGrid::new(3, 2, 1.0 / 3.0, 0.5, 0.6).unwrap();
    let topology = Boundaries2d {
        x: OPEN,
        y: AxisBoundaries::Periodic,
    };
    let lower = ReservoirBoundary {
        side: BoundarySide::XMin,
        partial_densities: &[0.8, 0.2, 0.0],
        velocity: [0.2, 0.1],
        pressure_pa: 1.2,
    };
    let upper = ReservoirBoundary {
        side: BoundarySide::XMax,
        partial_densities: &[0.2, 0.8, 0.0],
        velocity: [0.1, -0.1],
        pressure_pa: 0.9,
    };
    let raw = || {
        field(equal_mass(), grid, topology, |_| {
            ([0.5, 0.5, 0.0], [0.1, 0.0], 1.0)
        })
    };
    let construct = |inputs: &[ReservoirBoundary<'_>], budget: usize| {
        let f = raw();
        MixtureEuler2d::new_with_reservoirs(
            f.gas,
            f.grid,
            f.states,
            topology,
            inputs,
            budget,
            || false,
        )
    };
    let bytes =
        MixtureEuler2d::required_workspace_bytes_with_reservoirs(raw().gas(), &grid, 2).unwrap();
    assert!(construct(&[lower, upper], bytes - 1).is_err());
    assert!(construct(&[lower, lower], bytes).is_err());
    assert!(
        construct(
            &[
                lower,
                ReservoirBoundary {
                    side: BoundarySide::YMin,
                    ..upper
                }
            ],
            bytes
        )
        .is_err()
    );
    assert!(
        construct(
            &[
                lower,
                ReservoirBoundary {
                    pressure_pa: 100.0,
                    ..upper
                }
            ],
            bytes
        )
        .is_err()
    );
    let incomplete = bind(raw(), &[lower]);
    let source = PlanarTransport::new(
        incomplete.gas().clone(),
        law(incomplete.gas(), 0.02, 0.03, 0.01, 0.02),
        || false,
    )
    .unwrap();
    let budget = MolecularTransport2d::required_workspace_bytes(&incomplete, &source).unwrap();
    assert!(MolecularTransport2d::new(&incomplete, source, budget, || false).is_err());

    let make = || {
        let flow = bind(raw(), &[lower, upper]);
        let op = prepare(&flow, law(flow.gas(), 0.02, 0.03, 0.01, 0.02));
        (flow, op)
    };
    let (mut reference, mut reference_operator) = make();
    let mut polls = 0;
    let expected = reference_operator
        .advance(
            &mut reference,
            Time::new(0.003),
            DissipationControl::default(),
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    let (mut flow, mut op) = make();
    let original = flow.states().to_vec();
    let workspace = flow.workspace_bytes() + op.workspace_bytes();
    let mut seen = 0;
    let failure = op.advance(
        &mut flow,
        Time::new(0.003),
        DissipationControl::default(),
        || {
            seen += 1;
            seen == polls
        },
    );
    assert!(matches!(failure, Err(DissipationError::Cancelled)));
    assert_eq!(flow.states(), original);
    assert_eq!(flow.time().value(), 0.0);
    assert!(op.last_step().is_none());
    assert_eq!(flow.reservoir(BoundarySide::XMin), Some(lower));
    assert_eq!(flow.workspace_bytes() + op.workspace_bytes(), workspace);
    let replay = op
        .advance(
            &mut flow,
            Time::new(0.003),
            DissipationControl::default(),
            || false,
        )
        .unwrap();
    assert_eq!(replay, expected);
    assert_eq!(flow.states(), reference.states());
    assert_eq!(
        super::molecular::impulses(&op),
        super::molecular::impulses(&reference_operator)
    );
}
