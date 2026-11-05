//! Independent prescribed-bath physics, stage timing and complete transactions.

use super::*;
use crate::compressible::spatial::ReservoirFaceState;
use fs_gas::multicomponent::dissipation::{DissipationError, Rkl2Scheme};

const OPEN_AXIS: AxisBoundaries = AxisBoundaries::Ends {
    lower: BoundaryCondition::Transmissive,
    upper: BoundaryCondition::Transmissive,
};
const OPEN: Boundaries3d = Boundaries3d {
    x: OPEN_AXIS,
    y: OPEN_AXIS,
    z: OPEN_AXIS,
};

fn bath(temperature: f64) -> ReservoirFaceState {
    ReservoirFaceState {
        density: 1.0,
        velocity: [0.0; 3],
        pressure_pa: temperature,
    }
}

fn profiled(boundaries: Boundaries3d, sides: &[BoundarySide]) -> Euler3d {
    let mut flow = field(
        gas(0.0),
        CartesianGrid::new(1, 1, 1, 1.0, 1.0, 1.0).unwrap(),
        boundaries,
        |_| [1.0, 0.0, 0.0, 0.0, 2.0],
    );
    let bytes = flow
        .required_reservoir_profile_workspace_bytes(sides)
        .unwrap();
    flow.prepare_reservoir_profiles(sides, bytes, |_, _, _| Ok(bath(2.1)), || false)
        .unwrap();
    flow
}

#[allow(clippy::too_many_arguments)]
fn advance(
    op: &mut MolecularTransport3d<ConstantTransport>,
    flow: &mut Euler3d,
    stabilized: bool,
    start: f64,
    duration: f64,
    maximum_step: f64,
    budget: usize,
    sample: impl FnMut(f64, BoundarySide, usize, [f64; 3]) -> Result<ReservoirFaceState, Error>,
    cancelled: impl FnMut() -> bool,
) -> Result<MolecularReport3d, MolecularError> {
    let controls = MolecularControl {
        maximum_step_s: maximum_step,
        ..MolecularControl::default()
    };
    if stabilized {
        op.advance_stabilized_with_reservoir_profiles(
            flow,
            Time::new(start),
            Time::new(duration),
            controls,
            MolecularStabilization::default(),
            budget,
            sample,
            cancelled,
        )
        .map(|r| r.transport)
    } else {
        op.advance_with_reservoir_profiles(
            flow,
            Time::new(start),
            Time::new(duration),
            controls,
            budget,
            sample,
            cancelled,
        )
    }
}

fn coupled_step(
    flow: &mut MolecularEuler3d<ConstantTransport>,
    stabilized: bool,
    duration: f64,
    budget: usize,
    sample: impl FnMut(f64, BoundarySide, usize, [f64; 3]) -> Result<ReservoirFaceState, Error>,
    cancelled: impl FnMut() -> bool,
) -> Result<MolecularStepReport3d, MolecularError> {
    if stabilized {
        flow.step_stabilized_with_reservoir_profiles(
            Time::new(duration),
            MolecularStepControl::default(),
            MolecularStabilization::default(),
            budget,
            sample,
            cancelled,
        )
        .map(|r| r.transport)
    } else {
        flow.step_with_reservoir_profiles(
            Time::new(duration),
            MolecularStepControl::default(),
            budget,
            sample,
            cancelled,
        )
    }
}

#[test]
fn g0_rkl_times_satisfy_every_native_legendre_recurrence() {
    for stages in 2..=256 {
        let scheme = Rkl2Scheme::new(stages).unwrap();
        assert_eq!(scheme.abscissa(0).unwrap(), 0.0);
        assert_eq!(scheme.abscissa(1).unwrap(), scheme.first_weight());
        assert_eq!(scheme.abscissa(stages).unwrap(), 1.0);
        assert!(scheme.abscissa(stages + 1).is_err());
        for j in 2..=stages {
            let [mu, nu, mu_tilde, gamma_tilde] = scheme.coefficients(j).unwrap();
            let c = scheme.abscissa(j).unwrap();
            close(
                c,
                mu * scheme.abscissa(j - 1).unwrap()
                    + nu * scheme.abscissa(j - 2).unwrap()
                    + mu_tilde
                    + gamma_tilde,
                2e-15,
            );
            assert!(c > scheme.abscissa(j - 1).unwrap());
        }
    }
}

#[test]
fn g0_constant_timed_baths_retain_held_ssp_rkl_and_split_bits() {
    for stabilized in [false, true] {
        let make = || {
            let flow = profiled(OPEN, &BoundarySide::ALL);
            let op = MolecularTransport3d::new(&flow, law(flow.gas(), 0.02, 0.003, 0.3), || false)
                .unwrap()
                .with_reservoir_exchange();
            (flow, op)
        };
        let (mut held, mut held_op) = make();
        let (mut timed, mut timed_op) = make();
        let bytes = timed.reservoir_profile_workspace_bytes();
        let controls = MolecularControl {
            maximum_step_s: 0.01,
            ..MolecularControl::default()
        };
        let expected = if stabilized {
            held_op
                .advance_stabilized(
                    &mut held,
                    Time::new(0.04),
                    controls,
                    MolecularStabilization::default(),
                    || false,
                )
                .unwrap()
                .transport
        } else {
            held_op
                .advance(&mut held, Time::new(0.04), controls, || false)
                .unwrap()
        };
        let actual = advance(
            &mut timed_op,
            &mut timed,
            stabilized,
            1.25,
            0.04,
            0.01,
            1000,
            |_, _, _, _| Ok(bath(2.1)),
            || false,
        )
        .unwrap();
        assert_eq!(actual, expected);
        assert_eq!(timed.states(), held.states());
        assert_eq!(timed.time().value(), 0.0);
        assert_eq!(timed.reservoir_profile_workspace_bytes(), bytes);
        for side in BoundarySide::ALL {
            assert_eq!(
                timed_op.boundary_impulses(side),
                held_op.boundary_impulses(side)
            );
        }
        let make_coupled = || {
            let flow = profiled(OPEN, &BoundarySide::ALL);
            let source = law(flow.gas(), 0.02, 0.003, 0.3);
            MolecularEuler3d::new(flow, source, || false)
                .unwrap()
                .with_reservoir_exchange()
        };
        let mut held = make_coupled();
        let mut timed = make_coupled();
        for _ in 0..2 {
            let expected = if stabilized {
                held.step_stabilized(
                    Time::new(0.01),
                    MolecularStepControl::default(),
                    MolecularStabilization::default(),
                    || false,
                )
                .unwrap()
                .transport
            } else {
                held.step(Time::new(0.01), MolecularStepControl::default(), || false)
                    .unwrap()
            };
            let actual = coupled_step(
                &mut timed,
                stabilized,
                0.01,
                1000,
                |_, _, _, _| Ok(bath(2.1)),
                || false,
            )
            .unwrap();
            assert_eq!(actual, expected);
            assert_eq!(timed.flow().states(), held.flow().states());
            for side in BoundarySide::ALL {
                assert_eq!(timed.boundary_impulses(side), held.boundary_impulses(side));
            }
        }
    }
}

// Exact response of dX/dt = lambda*(amplitude*sin(omega*t) - X), X(0)=0.
fn driven(lambda: f64, amplitude: f64, omega: f64, time: f64) -> f64 {
    amplitude * lambda / (lambda * lambda + omega * omega)
        * (lambda * (omega * time).sin() - omega * (omega * time).cos()
            + omega * (-lambda * time).exp())
}

#[test]
fn g1_oscillating_six_face_heat_bath_has_second_order_absolute_time_response() {
    let start = 0.75;
    let duration = 0.4;
    let conductivity = 0.3;
    let omega = 2.0;
    // Six physical half-cell faces on the unit cube: lambda=12*k/(rho*Cv).
    let lambda = 12.0 * conductivity / 2.5;
    let expected = 2.0 + driven(lambda, 0.3, omega, duration);
    for stabilized in [false, true] {
        let mut errors = Vec::new();
        for step in [0.1, 0.05, 0.025] {
            let mut flow = profiled(OPEN, &BoundarySide::ALL);
            let mut op =
                MolecularTransport3d::new(&flow, law(flow.gas(), 0.0, 0.0, conductivity), || false)
                    .unwrap()
                    .with_reservoir_exchange();
            let report = advance(
                &mut op,
                &mut flow,
                stabilized,
                start,
                duration,
                step,
                10000,
                |time, _, _, _| Ok(bath(2.0 + 0.3 * (omega * (time - start)).sin())),
                || false,
            )
            .unwrap();
            let actual = flow.states()[0].primitive(flow.gas(), || false).unwrap();
            errors.push((actual.temperature - expected).abs());
            assert_eq!(actual.velocity, [0.0; 3]);
            assert_eq!(flow.states()[0].density, 1.0);
            assert_eq!(flow.time().value(), 0.0);
            close(
                report.boundary_change.energy,
                flow.totals().energy - report.before.energy,
                3e-13,
            );
            for residual in report.conservation_residual.array() {
                close(residual, 0.0, 3e-13);
            }
            assert_eq!(
                flow.reservoir_face(BoundarySide::XMin, 0),
                Some(bath(
                    2.0 + 0.3 * (omega * ((start + duration) - start)).sin()
                ))
            );
        }
        assert!(errors[0] < 0.002, "{stabilized}: {errors:?}");
        assert!(errors[1] < 0.30 * errors[0], "{stabilized}: {errors:?}");
        assert!(errors[2] < 0.30 * errors[1], "{stabilized}: {errors:?}");
    }
}

#[test]
fn g1_accelerated_bath_samples_actual_rkl_abscissae_and_endpoint() {
    let mut flow = profiled(OPEN, &BoundarySide::ALL);
    let mut op = MolecularTransport3d::new(&flow, law(flow.gas(), 0.0, 0.0, 20.0), || false)
        .unwrap()
        .with_reservoir_exchange();
    let start = 0.5;
    let dt = 0.1;
    let mut times = Vec::new();
    let report = op
        .advance_stabilized_with_reservoir_profiles(
            &mut flow,
            Time::new(start),
            Time::new(dt),
            MolecularControl::default(),
            MolecularStabilization::default(),
            1000,
            |time, side, face, _| {
                if side == BoundarySide::XMin && face == 0 {
                    times.push(time);
                }
                Ok(bath(2.0 + 0.01 * (time - start)))
            },
            || false,
        )
        .unwrap();
    assert_eq!(report.transport.accepted_steps, 1);
    assert_eq!(report.transport.rejected_steps, 0);
    assert!(report.maximum_stages_used > 2);
    let scheme = Rkl2Scheme::new(report.maximum_stages_used).unwrap();
    let expected: Vec<_> = (0..=scheme.stage_count())
        .map(|j| start + dt * scheme.abscissa(j).unwrap())
        .collect();
    assert_eq!(times, expected);
    assert_eq!(report.operator_evaluations, times.len());
    assert_eq!(
        flow.reservoir_face(BoundarySide::XMin, 0),
        Some(bath(2.0 + 0.01 * dt))
    );
    for residual in report.transport.conservation_residual.array() {
        close(residual, 0.0, 4e-13);
    }
}

#[test]
fn g4_late_source_budget_and_cancellation_refuse_the_complete_molecular_interval() {
    for stabilized in [false, true] {
        let make = || {
            let mut flow = profiled(OPEN, &BoundarySide::ALL);
            let mut op =
                MolecularTransport3d::new(&flow, law(flow.gas(), 0.02, 0.003, 0.3), || false)
                    .unwrap()
                    .with_reservoir_exchange();
            advance(
                &mut op,
                &mut flow,
                stabilized,
                0.0,
                0.01,
                0.01,
                1000,
                |time, _, _, _| Ok(bath(2.1 + time)),
                || false,
            )
            .unwrap();
            (flow, op)
        };
        let (mut flow, mut op) = make();
        let cells = flow.states().to_vec();
        let traces = BoundarySide::ALL.map(|s| flow.reservoir_face(s, 0));
        let loads = BoundarySide::ALL.map(|s| op.boundary_impulses(s).to_vec());
        let bytes = flow.reservoir_profile_workspace_bytes();
        let mut calls = 0;
        let result = advance(
            &mut op,
            &mut flow,
            stabilized,
            0.01,
            0.1,
            0.025,
            18,
            |time, _, _, _| {
                calls += 1;
                Ok(bath(2.1 + time))
            },
            || false,
        );
        assert_eq!(calls, 18);
        assert_eq!(
            result,
            Err(MolecularError::Flow(Error::Invalid(
                "timed reservoir profile face evaluation budget exhausted"
            )))
        );
        let mut calls = 0;
        let result = advance(
            &mut op,
            &mut flow,
            stabilized,
            0.01,
            0.1,
            0.025,
            1000,
            |time, _, _, _| {
                calls += 1;
                if calls == 13 {
                    Err(Error::Invalid("physical trace unavailable"))
                } else {
                    Ok(bath(2.1 + time))
                }
            },
            || false,
        );
        assert_eq!(calls, 13);
        assert_eq!(
            result,
            Err(MolecularError::Flow(Error::Invalid(
                "physical trace unavailable"
            )))
        );
        let (mut probe, mut probe_op) = make();
        let mut total_polls = 0;
        let expected = advance(
            &mut probe_op,
            &mut probe,
            stabilized,
            0.01,
            0.1,
            0.025,
            1000,
            |time, _, _, _| Ok(bath(2.1 + time)),
            || {
                total_polls += 1;
                false
            },
        )
        .unwrap();
        let mut polls = 0;
        assert_eq!(
            advance(
                &mut op,
                &mut flow,
                stabilized,
                0.01,
                0.1,
                0.025,
                1000,
                |time, _, _, _| Ok(bath(2.1 + time)),
                || {
                    polls += 1;
                    polls == total_polls
                }
            ),
            Err(MolecularError::Cancelled)
        );
        assert_eq!(flow.states(), cells);
        assert_eq!(flow.time().value(), 0.0);
        assert_eq!(flow.reservoir_profile_workspace_bytes(), bytes);
        for (i, side) in BoundarySide::ALL.into_iter().enumerate() {
            assert_eq!(flow.reservoir_face(side, 0), traces[i]);
            assert_eq!(op.boundary_impulses(side), loads[i]);
        }
        let actual = advance(
            &mut op,
            &mut flow,
            stabilized,
            0.01,
            0.1,
            0.025,
            1000,
            |time, _, _, _| Ok(bath(2.1 + time)),
            || false,
        )
        .unwrap();
        assert_eq!(actual, expected);
        assert_eq!(flow.states(), probe.states());
    }
}

struct TemperatureConductivity(ConstantTransport);

impl TransportLaw for TemperatureConductivity {
    fn source(&self) -> &IdealGasEquilibriumV1 {
        self.0.source()
    }

    fn evaluate(
        &mut self,
        temperature: Temperature,
        pressure: Pressure,
        fractions: &[f64],
        binary: &mut [f64],
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<TransportCoefficients, DissipationError> {
        let mut value = self
            .0
            .evaluate(temperature, pressure, fractions, binary, cancelled)?;
        value.thermal_conductivity = ThermalConductivityQuantityV1::new(
            value.thermal_conductivity.value() * temperature.value().powi(2),
        );
        Ok(value)
    }
}

#[test]
fn g2_hot_native_law_retries_resample_smaller_ssp_and_rkl_physical_times() {
    for stabilized in [false, true] {
        let mut flow = profiled(OPEN, &BoundarySide::ALL);
        let source = TemperatureConductivity(law(flow.gas(), 0.0, 0.0, 0.01));
        let mut op = MolecularTransport3d::new(&flow, source, || false)
            .unwrap()
            .with_reservoir_exchange();
        let mut times = Vec::new();
        let sample = |time: f64, side, _, _| {
            if side == BoundarySide::XMin {
                times.push(time);
            }
            Ok(bath(2.0 + 380.0 * time))
        };
        let report = if stabilized {
            op.advance_stabilized_with_reservoir_profiles(
                &mut flow,
                Time::new(0.0),
                Time::new(0.1),
                MolecularControl::default(),
                MolecularStabilization::default(),
                10000,
                sample,
                || false,
            )
            .unwrap()
            .transport
        } else {
            op.advance_with_reservoir_profiles(
                &mut flow,
                Time::new(0.0),
                Time::new(0.1),
                MolecularControl::default(),
                10000,
                sample,
                || false,
            )
            .unwrap()
        };
        assert!(report.rejected_steps > 0);
        assert!(times.windows(2).any(|t| t[1] < t[0]));
        assert_eq!(times.last(), Some(&0.1));
        assert_eq!(flow.reservoir_face(BoundarySide::XMin, 0), Some(bath(40.0)));
        for residual in report.conservation_residual.array() {
            close(residual, 0.0, 3e-12);
        }
    }
}

#[test]
fn g3_g4_split_times_budget_and_later_held_calls_share_one_atomic_bath() {
    for stabilized in [false, true] {
        let make = || {
            let flow = profiled(OPEN, &BoundarySide::ALL);
            let source = law(flow.gas(), 0.02, 0.003, 0.3);
            MolecularEuler3d::new(flow, source, || false)
                .unwrap()
                .with_reservoir_exchange()
        };
        let mut flow = make();
        let mut reference = make();
        let dt = 0.01;
        let mut times = Vec::new();
        let report = coupled_step(
            &mut flow,
            stabilized,
            dt,
            1000,
            |time, side, _, _| {
                if side == BoundarySide::XMin {
                    times.push(time);
                }
                Ok(bath(2.1 + time))
            },
            || false,
        )
        .unwrap();
        coupled_step(
            &mut reference,
            stabilized,
            dt,
            1000,
            |time, _, _, _| Ok(bath(2.1 + time)),
            || false,
        )
        .unwrap();
        let stage = 0.5 * dt * Rkl2Scheme::new(2).unwrap().abscissa(1).unwrap();
        let expected = if stabilized {
            vec![
                0.0,
                0.0,
                stage,
                dt / 2.0,
                0.0,
                dt,
                dt / 2.0,
                dt / 2.0 + stage,
                dt,
            ]
        } else {
            vec![0.0, 0.0, dt / 2.0, 0.0, dt, dt / 2.0, dt]
        };
        assert_eq!(times, expected);
        assert_eq!(report.coupling_retries, 0);
        assert_eq!(flow.flow().time().value(), dt);
        assert_eq!(
            flow.flow().reservoir_face(BoundarySide::XMin, 0),
            Some(bath(2.1 + dt))
        );
        let cells = flow.flow().states().to_vec();
        let last = flow.last_step().copied();
        let euler_last = flow.flow().last_step;
        let loads = BoundarySide::ALL.map(|s| flow.boundary_impulses(s).to_vec());
        let mut calls = 0;
        let budget = if stabilized { 36 } else { 30 }; // Exhaust after the central Euler interval.
        assert_eq!(
            coupled_step(
                &mut flow,
                stabilized,
                dt,
                budget,
                |time, _, _, _| {
                    calls += 1;
                    Ok(bath(2.1 + time))
                },
                || false
            ),
            Err(MolecularError::Flow(Error::Invalid(
                "timed reservoir profile face evaluation budget exhausted"
            )))
        );
        assert_eq!(calls, budget);
        let mut calls = 0;
        assert_eq!(
            coupled_step(
                &mut flow,
                stabilized,
                dt,
                1000,
                |time, _, _, _| {
                    calls += 1;
                    if calls == budget + 1 {
                        Err(Error::Invalid("late split trace unavailable"))
                    } else {
                        Ok(bath(2.1 + time))
                    }
                },
                || false
            ),
            Err(MolecularError::Flow(Error::Invalid(
                "late split trace unavailable"
            )))
        );
        assert_eq!(calls, budget + 1);
        let mut probe = make();
        coupled_step(
            &mut probe,
            stabilized,
            dt,
            1000,
            |time, _, _, _| Ok(bath(2.1 + time)),
            || false,
        )
        .unwrap();
        let mut total_polls = 0;
        coupled_step(
            &mut probe,
            stabilized,
            dt,
            1000,
            |time, _, _, _| Ok(bath(2.1 + time)),
            || {
                total_polls += 1;
                false
            },
        )
        .unwrap();
        let mut polls = 0;
        assert_eq!(
            coupled_step(
                &mut flow,
                stabilized,
                dt,
                1000,
                |time, _, _, _| Ok(bath(2.1 + time)),
                || {
                    polls += 1;
                    polls == total_polls
                }
            ),
            Err(MolecularError::Cancelled)
        );
        assert_eq!(flow.flow().states(), cells);
        assert_eq!(flow.flow().time().value(), dt);
        assert_eq!(flow.last_step().copied(), last);
        assert_eq!(flow.flow().last_step, euler_last);
        for (i, side) in BoundarySide::ALL.into_iter().enumerate() {
            assert_eq!(flow.flow().reservoir_face(side, 0), Some(bath(2.1 + dt)));
            assert_eq!(flow.boundary_impulses(side), loads[i]);
        }
        // The failed private split reached a later bath. Two later held swaps
        // must still use exactly the previously accepted endpoint.
        for _ in 0..2 {
            let actual = flow
                .step(Time::new(dt), MolecularStepControl::default(), || false)
                .unwrap();
            let expected = reference
                .step(Time::new(dt), MolecularStepControl::default(), || false)
                .unwrap();
            assert_eq!(actual, expected);
            assert_eq!(flow.flow().states(), reference.flow().states());
        }
    }
}

#[test]
fn g2_central_cfl_restarts_the_whole_timed_split_without_refunding_provider_work() {
    for stabilized in [false, true] {
        let make = || {
            let flow = profiled(OPEN, &BoundarySide::ALL);
            let transport = law(flow.gas(), 0.0, 0.0, 0.002);
            MolecularEuler3d::new(flow, transport, || false)
                .unwrap()
                .with_reservoir_exchange()
        };
        let mut wide = make();
        let mut wide_calls = 0;
        let report = coupled_step(
            &mut wide,
            stabilized,
            0.07,
            10000,
            |time, _, _, _| {
                wide_calls += 1;
                Ok(bath(2.0 + 400.0 * time))
            },
            || false,
        )
        .unwrap();
        assert!(report.coupling_retries > 0);
        let mut direct = make();
        let mut direct_calls = 0;
        let replay = coupled_step(
            &mut direct,
            stabilized,
            report.elapsed.value(),
            10000,
            |time, _, _, _| {
                direct_calls += 1;
                Ok(bath(2.0 + 400.0 * time))
            },
            || false,
        )
        .unwrap();
        assert_eq!(replay.coupling_retries, 0);
        assert_eq!(replay.first, report.first);
        assert_eq!(replay.advection, report.advection);
        assert_eq!(replay.second, report.second);
        assert_eq!(wide.flow().states(), direct.flow().states());
        assert!(wide_calls > direct_calls);
        let mut limited = make();
        let cells = limited.flow().states().to_vec();
        assert_eq!(
            coupled_step(
                &mut limited,
                stabilized,
                0.07,
                direct_calls,
                |time, _, _, _| Ok(bath(2.0 + 400.0 * time)),
                || false
            ),
            Err(MolecularError::Flow(Error::Invalid(
                "timed reservoir profile face evaluation budget exhausted"
            )))
        );
        assert_eq!(limited.flow().states(), cells);
        assert_eq!(limited.flow().time().value(), 0.0);
        assert_eq!(limited.last_step(), None);
        assert_eq!(
            limited.flow().reservoir_face(BoundarySide::XMin, 0),
            Some(bath(2.1))
        );
    }
}

#[test]
fn g0_split_profiles_end_at_the_canonical_outer_clock_without_rounding_outside() {
    for stabilized in [false, true] {
        let initial = profiled(OPEN, &BoundarySide::ALL);
        let transport = law(initial.gas(), 0.02, 0.003, 0.3);
        let mut flow = MolecularEuler3d::new(initial, transport, || false)
            .unwrap()
            .with_reservoir_exchange();
        flow.step(Time::new(0.03), MolecularStepControl::default(), || false)
            .unwrap();
        let start = flow.flow().time().value();
        let end = 0.04;
        let dt = end - start;
        assert!((start + 0.5 * dt) + 0.5 * dt > end);
        let mut last = start;
        let report = coupled_step(
            &mut flow,
            stabilized,
            dt,
            1000,
            |time, _, _, _| {
                if !(start..=end).contains(&time) {
                    return Err(Error::Invalid(
                        "observation outside admitted source interval",
                    ));
                }
                last = time;
                Ok(bath(2.1 + time))
            },
            || false,
        )
        .unwrap();
        assert_eq!(report.elapsed.value(), dt);
        assert_eq!(flow.flow().time().value(), end);
        assert_eq!(last, end);
        assert_eq!(
            flow.flow().reservoir_face(BoundarySide::XMin, 0),
            Some(bath(2.1 + end))
        );
    }
}

#[test]
fn g1_nonauto_strang_resolves_manufactured_shear_and_isothermal_bath_work() {
    let boundaries = Boundaries3d {
        x: OPEN_AXIS,
        ..PERIODIC
    };
    let sides = [BoundarySide::XMin, BoundarySide::XMax];
    let (mu, conductivity, duration) = (0.1, 0.3, 0.5);
    let lambda = 4.0 * mu;
    let exact = |time| driven(lambda, 0.3, 3.0, time);
    let source = |time: f64| {
        let velocity = 0.3 * (3.0 * time).sin();
        // Subtract precisely the viscous conversion into internal heat. With
        // rho=R=1 the native continuum solution keeps T=p=2 and Euler is an
        // exact stationary contact; both operators still execute every stage.
        let temperature = 2.0 - mu / conductivity * (velocity - exact(time)).powi(2);
        ReservoirFaceState {
            density: 2.0 / temperature,
            velocity: [0.0, velocity, 0.0],
            pressure_pa: 2.0,
        }
    };
    for stabilized in [false, true] {
        let mut errors = Vec::new();
        for steps in [8, 16, 32] {
            let mut initial = profiled(boundaries, &sides);
            initial
                .update_reservoir_profiles(|_, _, _| Ok(source(0.0)), || false)
                .unwrap();
            let transport = law(initial.gas(), mu, 0.0, conductivity);
            let mut flow = MolecularEuler3d::new(initial, transport, || false)
                .unwrap()
                .with_reservoir_exchange();
            let dt = duration / steps as f64;
            for _ in 0..steps {
                let report = coupled_step(
                    &mut flow,
                    stabilized,
                    dt,
                    1000,
                    |time, _, _, _| Ok(source(time)),
                    || false,
                )
                .unwrap();
                assert_eq!(report.elapsed.value(), dt);
                assert_eq!(report.coupling_retries, 0);
                for residual in report.conservation_residual.array() {
                    close(residual, 0.0, 4e-13);
                }
            }
            let actual = flow.flow().states()[0]
                .primitive(flow.flow().gas(), || false)
                .unwrap();
            errors.push(
                (actual.velocity[1] - exact(duration)).abs()
                    + (actual.temperature - 2.0).abs()
                    + (flow.flow().states()[0].density - 1.0).abs(),
            );
        }
        assert!(errors[0] < 0.002, "{stabilized}: {errors:?}");
        assert!(errors[1] < 0.35 * errors[0], "{stabilized}: {errors:?}");
        assert!(errors[2] < 0.35 * errors[1], "{stabilized}: {errors:?}");
    }
}
