//! G1/G3/G4: the same thermal stratification in full planar Euler dynamics.
use super::*;

fn lapse(gradient: f64) -> LapseRateAtmosphere {
    LapseRateAtmosphere {
        acceleration: Acceleration::new(-3.0),
        reference_temperature: Temperature::new(2.0),
        reference_pressure: Pressure::new(3.0),
        reference_position: Length::new(0.0),
        temperature_gradient_k_m: gradient,
    }
}
fn field(axis: usize, gradient: f64) -> Euler2d {
    let (grid, boundaries) = if axis == 0 {
        (
            CartesianGrid::new(20, 5, 0.05, 0.2, 0.7).unwrap(),
            Boundaries2d {
                x: WALLS,
                y: AxisBoundaries::Periodic,
            },
        )
    } else {
        (
            CartesianGrid::new(5, 20, 0.2, 0.05, 0.7).unwrap(),
            Boundaries2d {
                x: AxisBoundaries::Periodic,
                y: WALLS,
            },
        )
    };
    Euler2d::new_lapse_rate_atmosphere(
        gas(0.4),
        grid,
        axis,
        lapse(gradient),
        AtmosphereQuadrature::default(),
        boundaries,
        || false,
    )
    .unwrap()
}

#[test]
fn g1_stable_and_unstable_rest_share_one_integrated_profile_on_either_axis() {
    for axis in 0..2 {
        for gradient in [-0.4, -1.2] {
            for reconstruction in [Reconstruction::FirstOrder, Reconstruction::Muscl] {
                let mut flow = field(axis, gradient);
                assert_eq!(
                    flow.hydrostatic_reference(),
                    Some((axis, HydrostaticReference::LapseRate(lapse(gradient))))
                );
                let initial = flow.states().to_vec();
                let end_p = 3.0 * (1.0 + gradient / 2.0).powf(-3.0 / gradient);
                let mass = 0.7 * (end_p - 3.0) / -3.0;
                close(flow.totals().mass, mass, 3e-11);
                for _ in 0..16 {
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
                    assert_eq!(report.source_change.work, 0.0);
                    assert_eq!(report.halvings, 0);
                    close(
                        report.source_change.impulse[axis],
                        -3.0 * mass * report.elapsed.value(),
                        3e-11,
                    );
                    balance(report);
                }
                assert_eq!(flow.states(), initial);
            }
        }
    }
}

#[test]
fn g3_lapse_reference_releases_two_dimensional_perturbations_conservatively() {
    let equilibrium = field(1, -0.4);
    let states = equilibrium
        .states()
        .iter()
        .enumerate()
        .map(|(cell, &q)| {
            let p = q.primitive(equilibrium.gas(), || false).unwrap();
            if cell % 5 == 2 && (7..13).contains(&(cell / 5)) {
                State::from_primitive(equilibrium.gas(), p.density, [0.0; 2], 1.04 * p.pressure)
                    .unwrap()
            } else {
                q
            }
        })
        .collect();
    let mut flow = Euler2d::new(
        equilibrium.gas.clone(),
        equilibrium.grid,
        states,
        equilibrium.boundaries,
        || false,
    )
    .unwrap();
    let initial = flow.states().to_vec();
    flow.bind_lapse_rate_atmosphere(1, lapse(-0.4), AtmosphereQuadrature::default(), || false)
        .unwrap();
    assert_eq!(flow.states(), initial);
    let mass = flow.totals().mass;
    let total = flow.totals().energy + flow.gravitational_potential_energy();
    for _ in 0..24 {
        let report = flow
            .step(Time::new(0.005), StepControl::default(), || false)
            .unwrap();
        balance(report);
    }
    assert!(flow.states().iter().any(|q| q.momentum[0].abs() > 1e-5));
    assert!(flow.states().iter().any(|q| q.momentum[1].abs() > 1e-5));
    close(flow.totals().mass, mass, 3e-12);
    close(
        flow.totals().energy + flow.gravitational_potential_energy(),
        total,
        3e-11,
    );
}

#[test]
fn g4_lapse_binding_refuses_replacement_and_publishes_nothing_on_cancellation() {
    let source = field(1, -0.4);
    let states = source.states().to_vec();
    let mut flow = Euler2d::new(
        source.gas.clone(),
        source.grid,
        states.clone(),
        source.boundaries,
        || false,
    )
    .unwrap();
    let mut polls = 0;
    assert_eq!(
        flow.bind_lapse_rate_atmosphere(1, lapse(-0.4), AtmosphereQuadrature::default(), || {
            polls += 1;
            polls == 30
        }),
        Err(Error::Cancelled)
    );
    assert_eq!(flow.states(), states);
    assert_eq!(flow.hydrostatic_reference(), None);
    flow.bind_lapse_rate_atmosphere(1, lapse(-0.4), AtmosphereQuadrature::default(), || false)
        .unwrap();
    assert!(
        flow.bind_lapse_rate_atmosphere(1, lapse(-1.2), AtmosphereQuadrature::default(), || false)
            .is_err()
    );
    assert_eq!(flow.states(), states);
    assert_eq!(
        flow.hydrostatic_reference(),
        Some((1, HydrostaticReference::LapseRate(lapse(-0.4))))
    );
}
