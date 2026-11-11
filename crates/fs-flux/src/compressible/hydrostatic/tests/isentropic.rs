//! G1/G3: independently integrated hydrostatics and evolving pressure waves.

use super::*;

fn isentrope(acceleration: f64) -> IsentropicAtmosphere {
    IsentropicAtmosphere {
        acceleration: Acceleration::new(acceleration),
        reference_temperature: Temperature::new(2.0),
        reference_pressure: Pressure::new(3.0),
        reference_position: Length::new(0.0),
    }
}

fn field(axis: usize, cp_linear: f64, acceleration: f64) -> Euler2d {
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
    Euler2d::new_isentropic_atmosphere(
        gas(cp_linear),
        grid,
        axis,
        isentrope(acceleration),
        AtmosphereQuadrature::default(),
        boundaries,
        || false,
    )
    .unwrap()
}

/// Independent analytic integration for Cp=3.5+a*T and R=1.
fn analytic(cp_linear: f64, acceleration: f64, position: f64) -> (f64, f64) {
    let cp0 = 3.5 + 2.0 * cp_linear;
    let enthalpy = acceleration * position;
    let dt = 2.0 * enthalpy / (cp0 + (cp0 * cp0 + 2.0 * cp_linear * enthalpy).sqrt());
    let temperature = 2.0 + dt;
    let pressure = 3.0 * (3.5 * (temperature / 2.0).ln() + cp_linear * dt).exp();
    (temperature, pressure)
}

#[test]
fn g1_variable_cp_lapse_and_integrated_pressure_balance_remain_at_rest() {
    for axis in 0..2 {
        for a in [0.0, 0.4] {
            for g in [-3.0, 3.0] {
                for reconstruction in [Reconstruction::FirstOrder, Reconstruction::Muscl] {
                    let mut flow = field(axis, a, g);
                    assert_eq!(flow.atmosphere(), None);
                    assert_eq!(
                        flow.hydrostatic_reference(),
                        Some((axis, HydrostaticReference::Isentropic(isentrope(g))))
                    );
                    let profile = flow.hydrostatic.as_ref().unwrap();
                    for (i, (p, _)) in profile.faces.iter().enumerate() {
                        let (temperature, pressure) = analytic(a, g, i as f64 * 0.05);
                        close(p.temperature, temperature, 3e-11);
                        close(p.pressure, pressure, 3e-11);
                    }
                    // Integral rho dx = (p_right-p_left)/g. This tests conserved
                    // cell integration, independently of the profile quadrature.
                    let expected_mass = 0.7 * (analytic(a, g, 1.0).1 - 3.0) / g;
                    close(flow.totals().mass, expected_mass, 3e-10);
                    if a != 0.0 {
                        let frozen_cp_temperature = 2.0 + g / (3.5 + 2.0 * a);
                        assert!((analytic(a, g, 1.0).0 - frozen_cp_temperature).abs() > 0.015);
                    }
                    let initial = flow.states().to_vec();
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
                        balance(report);
                    }
                    assert_eq!(flow.states(), initial);
                }
            }
        }
    }
}

/// A forced step keeps the step its retry loop gave it: the body-force owner kicks half a step before
/// and half a step after the attempt with that step, so the attempt must not shorten it (an unforced
/// attempt on a hydrostatic field does take the step its first stage's face rate admits, as the test
/// above shows). On this atmosphere the first attempt exceeds the CFL budget, because the face rate
/// along the stratified axis exceeds the estimate's cell-centre sum, so the forced step halves as it
/// always did, and the force impulse over the accepted step is the force times its elapsed time. That
/// holds only if the kicks and the stages used one step. The force acts along the periodic axis, where
/// the atmosphere has no gravity, so the whole y impulse is the force's. It is small enough (the
/// velocity it adds, about 1e-16, changes the endpoint's cell rate by far less than the CFL check's
/// 16 eps) that the forcing owner's own endpoint check admits the given step: only the first stage's
/// face rate can shorten it, which is what an attempt that wrongly took the step would show.
#[test]
fn a_forced_atmosphere_step_keeps_its_given_step_and_its_force_impulse_matches_the_elapsed_time() {
    let mut flow = field(0, 0.0, -3.0);
    let mass = flow.totals().mass;
    let acceleration = [0.0, 1e-14];
    let mut scratch = flow
        .forcing_workspace(flow.required_forcing_workspace_bytes().unwrap())
        .unwrap();
    let report = flow
        .step_with_forcing(
            Time::new(0.02),
            StepControl {
                reconstruction: Reconstruction::FirstOrder,
                ..StepControl::default()
            },
            BodyForce2d::UniformAcceleration(acceleration),
            ForceControl2d { maximum_mach: 8.0 },
            &mut scratch,
            || false,
        )
        .unwrap();
    assert_eq!(report.halvings, 1, "the forced first attempt exceeds its CFL budget and halves");
    // The impulse is about 1e-16: compare it relatively, not with `close`'s absolute floor.
    let expected = acceleration[1] * mass * report.elapsed.value();
    let impulse = report.source_change.impulse[1];
    assert!(
        (impulse - expected).abs() <= 1e-9 * expected.abs(),
        "impulse {impulse:e} against force x elapsed {expected:e}"
    );
}

#[test]
fn g3_isentropic_reference_allows_pressure_disturbances_and_conserves_total_energy() {
    let equilibrium = field(1, 0.4, -3.0);
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
    flow.bind_isentropic_atmosphere(1, isentrope(-3.0), AtmosphereQuadrature::default(), || {
        false
    })
    .unwrap();
    assert_eq!(flow.states(), initial);
    let mass = flow.totals().mass;
    let energy = flow.totals().energy + flow.gravitational_potential_energy();
    for _ in 0..48 {
        balance(
            flow.step(Time::new(0.01), StepControl::default(), || false)
                .unwrap(),
        );
    }
    assert!(flow.states().iter().any(|q| q.momentum[0].abs() > 1e-5));
    assert!(flow.states().iter().any(|q| q.momentum[1].abs() > 1e-5));
    assert!(
        flow.states()
            .iter()
            .zip(&initial)
            .any(|(new, old)| (new.density - old.density).abs() > 1e-5)
    );
    close(flow.totals().mass, mass, 2e-13);
    close(
        flow.totals().energy + flow.gravitational_potential_energy(),
        energy,
        2e-12,
    );
}

#[test]
fn g3_source_window_budget_reference_replacement_and_cancellation_are_atomic() {
    let equilibrium = field(1, 0.4, -3.0);
    let initial = equilibrium.states().to_vec();
    let unbound = || {
        Euler2d::new(
            equilibrium.gas.clone(),
            equilibrium.grid,
            initial.clone(),
            equilibrium.boundaries,
            || false,
        )
        .unwrap()
    };
    // Reaching h(T)<h(T_min) is a real thermodynamic refusal, not a clamp.
    let mut flow = unbound();
    assert!(
        flow.bind_isentropic_atmosphere(
            1,
            isentrope(-100.0),
            AtmosphereQuadrature::default(),
            || false
        )
        .is_err()
    );
    assert_eq!(flow.states(), initial);
    assert_eq!(flow.hydrostatic_reference(), None);
    assert!(
        flow.bind_isentropic_atmosphere(
            1,
            isentrope(-3.0),
            AtmosphereQuadrature {
                maximum_panels: 1,
                ..AtmosphereQuadrature::default()
            },
            || false,
        )
        .is_err()
    );
    let mut count = 0;
    flow.bind_isentropic_atmosphere(1, isentrope(-3.0), AtmosphereQuadrature::default(), || {
        count += 1;
        false
    })
    .unwrap();
    for cancel_at in [1, count / 2, count] {
        let mut trial = unbound();
        let mut polls = 0;
        assert_eq!(
            trial.bind_isentropic_atmosphere(
                1,
                isentrope(-3.0),
                AtmosphereQuadrature::default(),
                || {
                    polls += 1;
                    polls == cancel_at
                },
            ),
            Err(Error::Cancelled)
        );
        assert_eq!(trial.states(), initial);
        assert_eq!(trial.hydrostatic_reference(), None);
    }
    assert!(
        flow.bind_isothermal_atmosphere(1, reference(-3.0), || false)
            .is_err()
    );
    assert!(
        flow.bind_isentropic_atmosphere(
            1,
            isentrope(-3.0),
            AtmosphereQuadrature {
                relative_tolerance: 1e-9,
                ..AtmosphereQuadrature::default()
            },
            || false,
        )
        .is_err()
    );
    flow.bind_isentropic_atmosphere(1, isentrope(-3.0), AtmosphereQuadrature::default(), || {
        false
    })
    .unwrap();
    assert_eq!(flow.states(), initial);
    let report = flow
        .step(Time::new(0.005), StepControl::default(), || false)
        .unwrap();
    let loads = flow.boundary_impulses(BoundarySide::YMin).to_vec();
    let mut polls = 0;
    assert_eq!(
        flow.step(Time::new(0.005), StepControl::default(), || {
            polls += 1;
            polls == 400
        }),
        Err(Error::Cancelled)
    );
    assert_eq!(flow.states(), initial);
    assert_eq!(flow.last_step(), Some(&report));
    assert_eq!(flow.time(), report.elapsed);
    assert_eq!(flow.boundary_impulses(BoundarySide::YMin), loads);
}
