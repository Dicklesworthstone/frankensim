use super::*;
use crate::compressible::tests::{close, gas};
use fs_qty::{Acceleration, Pressure, Temperature};

mod isentropic;
mod lapse;

const WALLS: AxisBoundaries = AxisBoundaries::Ends {
    lower: BoundaryCondition::Reflecting,
    upper: BoundaryCondition::Reflecting,
};
const OPEN: AxisBoundaries = AxisBoundaries::Ends {
    lower: BoundaryCondition::Transmissive,
    upper: BoundaryCondition::Transmissive,
};

fn reference(acceleration: f64) -> IsothermalAtmosphere {
    IsothermalAtmosphere {
        acceleration: Acceleration::new(acceleration),
        temperature: Temperature::new(2.0),
        reference_pressure: Pressure::new(3.0),
        reference_position: Length::new(0.0),
    }
}

fn atmosphere(axis: usize, acceleration: f64, ends: AxisBoundaries) -> Euler2d {
    let (grid, boundaries) = if axis == 0 {
        (
            CartesianGrid::new(24, 6, 1.0 / 24.0, 1.0 / 6.0, 0.7).unwrap(),
            Boundaries2d {
                x: ends,
                y: AxisBoundaries::Periodic,
            },
        )
    } else {
        (
            CartesianGrid::new(6, 24, 1.0 / 6.0, 1.0 / 24.0, 0.7).unwrap(),
            Boundaries2d {
                x: AxisBoundaries::Periodic,
                y: ends,
            },
        )
    };
    Euler2d::new_isothermal_atmosphere(
        gas(0.4),
        grid,
        axis,
        reference(acceleration),
        boundaries,
        || false,
    )
    .unwrap()
}

fn perturbed(equilibrium: Euler2d, perturb: impl Fn(usize, State, Primitive) -> State) -> Euler2d {
    let states = equilibrium
        .states
        .iter()
        .enumerate()
        .map(|(cell, q)| {
            let p = q.primitive(equilibrium.gas(), || false).unwrap();
            perturb(cell, *q, p)
        })
        .collect();
    let (axis, reference) = equilibrium.atmosphere().unwrap();
    let mut flow = Euler2d::new(
        equilibrium.gas.clone(),
        equilibrium.grid,
        states,
        equilibrium.boundaries,
        || false,
    )
    .unwrap();
    flow.bind_isothermal_atmosphere(axis, reference, || false)
        .unwrap();
    flow
}

fn balance(report: StepReport) {
    assert!(report.maximum_cfl <= 0.40000000000001);
    let scale = report.before.mass.abs() + report.before.energy.abs();
    for residual in report.conservation_residual.array() {
        assert!(residual.abs() < 2e-12 * (1.0 + scale), "{report:?}");
    }
}

#[test]
fn exact_stratified_rest_uses_integrated_averages_both_axes_and_orders() {
    for axis in 0..2 {
        for g in [-12.0, 12.0] {
            for reconstruction in [Reconstruction::FirstOrder, Reconstruction::Muscl] {
                let mut flow = atmosphere(axis, g, WALLS);
                let initial = flow.states().to_vec();
                // 400:1 density contrast; the exact exponential volume integral
                // distinguishes these averages from centre-point samples.
                let beta = g / 2.0;
                let exact_mass = 0.7 * 1.5 * beta.exp_m1() / beta;
                close(flow.totals().mass, exact_mass, 3e-14);
                let max_density = initial.iter().map(|q| q.density).fold(0.0, f64::max);
                let min_density = initial
                    .iter()
                    .map(|q| q.density)
                    .fold(f64::INFINITY, f64::min);
                assert!(max_density / min_density > 300.0);
                for _ in 0..24 {
                    let report = flow
                        .step(
                            Time::new(0.1),
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
                        exact_mass * g * report.elapsed.value(),
                        5e-14,
                    );
                    balance(report);
                }
                assert_eq!(flow.states(), initial);
            }
        }
    }
}

#[test]
fn two_dimensional_disturbance_evolves_and_conserves_mass_and_total_energy() {
    let equilibrium = atmosphere(1, -8.0, WALLS);
    let source = equilibrium.gas().clone();
    let mut flow = perturbed(equilibrium, |cell, q, p| {
        if (2..4).contains(&(cell % 6)) && (8..14).contains(&(cell / 6)) {
            State::from_primitive(&source, q.density * 1.015, [0.0; 2], p.pressure * 1.04).unwrap()
        } else {
            q
        }
    });
    let initial = flow.states().to_vec();
    let mass = flow.totals().mass;
    let gas_energy = flow.totals().energy;
    let combined_energy = gas_energy + flow.gravitational_potential_energy();
    let mut work = 0.0;
    for _ in 0..64 {
        let report = flow
            .step(Time::new(0.02), StepControl::default(), || false)
            .unwrap();
        work += report.source_change.work;
        balance(report);
    }
    assert!(flow.states().iter().any(|q| q.momentum[0].abs() > 1e-5));
    assert!(flow.states().iter().any(|q| q.momentum[1].abs() > 1e-5));
    assert!(
        flow.states()
            .iter()
            .zip(&initial)
            .any(|(new, old)| (new.density - old.density).abs() > 1e-5)
    );
    assert!(
        flow.states().iter().any(|q| {
            (q.primitive(flow.gas(), || false).unwrap().temperature - 2.0).abs() > 1e-4
        })
    );
    close(flow.totals().mass, mass, 2e-13);
    close(
        flow.totals().energy + flow.gravitational_potential_energy(),
        combined_energy,
        2e-12,
    );
    close(flow.totals().energy - gas_energy, work, 2e-12);
}

#[test]
fn strongly_reduced_density_remains_positive_and_can_evolve() {
    // At |g*dx/(R*T)|=1, an additive perturbation of -90% of the
    // reference cell average makes its lower-density face negative.
    // Relative reconstruction must admit the real positive atmosphere.
    for reconstruction in [Reconstruction::FirstOrder, Reconstruction::Muscl] {
        for disturbed in [false, true] {
            let grid = CartesianGrid::new(3, 6, 1.0, 1.0, 0.7).unwrap();
            let equilibrium = Euler2d::new_isothermal_atmosphere(
                gas(0.4),
                grid,
                1,
                reference(-2.0),
                Boundaries2d {
                    x: AxisBoundaries::Periodic,
                    y: WALLS,
                },
                || false,
            )
            .unwrap();
            let mut flow = perturbed(equilibrium, |cell, q, _| State {
                density: 0.1 * q.density,
                momentum: [0.0; 2],
                energy: 0.1 * q.energy * if disturbed && cell == 7 { 1.03 } else { 1.0 },
            });
            let mass = flow.totals().mass;
            let energy = flow.totals().energy + flow.gravitational_potential_energy();
            for _ in 0..24 {
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
                balance(report);
            }
            assert!(flow.states().iter().all(|q| positive(q.density)));
            let maximum_speed = flow
                .states()
                .iter()
                .map(|q| {
                    q.momentum
                        .iter()
                        .map(|m| (m / q.density).abs())
                        .fold(0.0, f64::max)
                })
                .fold(0.0, f64::max);
            if disturbed {
                assert!(maximum_speed > 1e-5);
            } else {
                assert!(maximum_speed < 1e-8, "{maximum_speed:e}");
            }
            close(flow.totals().mass, mass, 2e-13);
            close(
                flow.totals().energy + flow.gravitational_potential_energy(),
                energy,
                2e-12,
            );
        }
    }
}

#[test]
fn open_boundaries_transport_potential_with_the_same_outward_mass() {
    let equilibrium = atmosphere(1, -4.0, OPEN);
    let source = equilibrium.gas().clone();
    let mut flow = perturbed(equilibrium, |_, q, p| {
        State::from_primitive(&source, q.density, [0.04, 0.12], p.pressure).unwrap()
    });
    let energy0 = flow.totals().energy + flow.gravitational_potential_energy();
    let mass0 = flow.totals().mass;
    let mut energy_transfer = 0.0;
    let mut mass_transfer = 0.0;
    for _ in 0..24 {
        let report = flow
            .step(Time::new(0.02), StepControl::default(), || false)
            .unwrap();
        // phi=4*y, hence zero at YMin and four at YMax. Transverse periodic
        // faces have identical fluxes and cancel in the potential account.
        let potential_out: f64 = flow
            .boundary_impulses(BoundarySide::YMax)
            .iter()
            .map(|load| 4.0 * load.mass)
            .sum();
        energy_transfer += report.boundary_change.energy - potential_out;
        mass_transfer += report.boundary_change.mass;
        balance(report);
    }
    assert!(mass_transfer.abs() > 1e-4);
    close(flow.totals().mass - mass0, mass_transfer, 2e-13);
    close(
        flow.totals().energy + flow.gravitational_potential_energy() - energy0,
        energy_transfer,
        2e-12,
    );
}

#[test]
fn additional_acceleration_uses_each_density_and_retains_gravity_work_account() {
    let mut flow = atmosphere(1, -6.0, WALLS);
    let mass = flow.totals().mass;
    let mut workspace = flow
        .forcing_workspace(flow.required_forcing_workspace_bytes().unwrap())
        .unwrap();
    let report = flow
        .step_with_forcing(
            Time::new(0.005),
            StepControl::default(),
            BodyForce2d::UniformAcceleration([0.7, 0.0]),
            ForceControl2d { maximum_mach: 1.0 },
            &mut workspace,
            || false,
        )
        .unwrap();
    let dv = 0.7 * report.elapsed.value();
    for q in flow.states() {
        close(q.momentum[0] / q.density, dv, 2e-10);
        close(q.momentum[1] / q.density, 0.0, 2e-10);
    }
    close(report.source_change.impulse[0], mass * dv, 2e-12);
    close(
        report.source_change.impulse[1],
        -6.0 * mass * report.elapsed.value(),
        2e-12,
    );
    close(report.source_change.work, 0.5 * mass * dv * dv, 2e-12);
    balance(report);
}

#[test]
fn invalid_reference_rebinding_and_cancelled_trials_publish_nothing() {
    let mut flow = atmosphere(1, -4.0, WALLS);
    let declaration = flow.atmosphere();
    let report = flow
        .step(Time::new(0.005), StepControl::default(), || false)
        .unwrap();
    let initial = flow.states().to_vec();
    let loads = flow.boundary_impulses(BoundarySide::YMin).to_vec();
    assert!(
        flow.bind_isothermal_atmosphere(1, reference(-3.0), || false)
            .is_err()
    );
    flow.bind_isothermal_atmosphere(1, reference(-4.0), || false)
        .unwrap();
    let mut polls = 0;
    assert_eq!(
        flow.step(Time::new(0.005), StepControl::default(), || {
            polls += 1;
            polls == 400
        }),
        Err(Error::Cancelled)
    );
    assert_eq!(flow.atmosphere(), declaration);
    assert_eq!(flow.last_step(), Some(&report));
    assert_eq!(flow.time(), report.elapsed);
    assert_eq!(flow.states(), initial);
    assert_eq!(flow.boundary_impulses(BoundarySide::YMin), loads);

    let mut unbound = Euler2d::new(
        flow.gas.clone(),
        flow.grid,
        initial.clone(),
        flow.boundaries,
        || false,
    )
    .unwrap();
    assert!(
        unbound
            .bind_isothermal_atmosphere(0, reference(-4.0), || false)
            .is_err()
    );
    assert_eq!(
        unbound.bind_isothermal_atmosphere(1, reference(-4.0), || true),
        Err(Error::Cancelled),
    );
    assert_eq!(unbound.atmosphere(), None);
    assert_eq!(unbound.states(), initial);
    unbound
        .step(Time::new(0.001), StepControl::default(), || false)
        .unwrap();
    assert!(
        unbound
            .bind_isothermal_atmosphere(1, reference(-4.0), || false)
            .is_err()
    );
}
