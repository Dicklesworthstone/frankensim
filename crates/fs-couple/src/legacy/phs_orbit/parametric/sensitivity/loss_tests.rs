//! G1/G3: original nonlinear loss, pump and physical parameter chains together.
use super::*;
use fs_phs::harmonic::equilibrium::LinearizationError;

const ORDER: usize = 12;
#[derive(Clone, Copy, Default)]
struct Loss {
    gamma: f64,
    sigma: f64,
    eta: f64,
}
impl Loss {
    fn source() -> Self {
        Self {
            gamma: 9.0,
            sigma: 0.3,
            eta: 0.25,
        }
    }
    fn slope() -> Self {
        Self {
            gamma: 1.3,
            sigma: -0.07,
            eta: 0.11,
        }
    }
    fn shifted(self, p: Self, h: f64) -> Self {
        Self {
            gamma: self.gamma + h * p.gamma,
            sigma: self.sigma + h * p.sigma,
            eta: self.eta + h * p.eta,
        }
    }
    fn flow(self, x: &[f64], e: &[f64], out: &mut [f64]) -> bool {
        let w = e[1] + self.eta * e[2];
        let d = self.gamma * (1.0 + self.sigma * x[0] * x[0]) * w.powi(3);
        out.copy_from_slice(&[0.0, d, self.eta * d]);
        true
    }
    fn tangent(self, x: &[f64], e: &[f64], dx: &[f64], de: &[f64], out: &mut [f64]) -> bool {
        let w = e[1] + self.eta * e[2];
        let d = self.gamma
            * (2.0 * self.sigma * x[0] * dx[0] * w.powi(3)
                + 3.0 * (1.0 + self.sigma * x[0] * x[0]) * w * w * (de[1] + self.eta * de[2]));
        out.copy_from_slice(&[0.0, d, self.eta * d]);
        true
    }
    fn partial(self, p: Self, x: &[f64], e: &[f64], out: &mut [f64]) -> bool {
        let w = e[1] + self.eta * e[2];
        let geometry = 1.0 + self.sigma * x[0] * x[0];
        let d = self.gamma * geometry * w.powi(3);
        let dp = (p.gamma * geometry + self.gamma * p.sigma * x[0] * x[0]) * w.powi(3)
            + 3.0 * self.gamma * geometry * w * w * p.eta * e[2];
        out.copy_from_slice(&[0.0, dp, self.eta * dp + p.eta * d]);
        true
    }
    // Independent physical source and four accounts; no production D, pump,
    // field, storage gradient, objective or parameter callback is invoked.
    fn physical(self, f: Family, time: f64, x: &[f64; 3]) -> ([f64; 3], [f64; 4]) {
        let theta = f.omega * time;
        let cross = f.cross * (theta + f.cross_phase).cos();
        let e1 = x[1] / 1.2;
        let e2 = f.mu * (x[2] - x[0]) + cross * x[0];
        let w = e1 + self.eta * e2;
        let d = self.gamma * (1.0 + self.sigma * x[0] * x[0]) * w.powi(3);
        let (mut flow, base) = f.physical(time, x);
        flow[1] -= d;
        flow[2] -= self.eta * d;
        (
            flow,
            [
                base[0],
                base[1],
                base[2],
                self.gamma * (1.0 + self.sigma * x[0] * x[0]) * w.powi(4),
            ],
        )
    }
    fn solve(self, f: Family, periods: usize) -> HbOrbit {
        let (system, input, pump) = (f.system(), f.input(), f.pump());
        let flow = |x: &[f64], e: &[f64], out: &mut [f64]| self.flow(x, e, out);
        let tangent = |x: &[f64], e: &[f64], dx: &[f64], de: &[f64], out: &mut [f64]| {
            self.tangent(x, e, dx, de, out)
        };
        let source = ParametricHamiltonian::new(&system, &input, &pump)
            .unwrap()
            .with_dissipation(NonlinearDissipation {
                flow: &flow,
                tangent: &tangent,
                max_tangent_actions: 3,
            })
            .unwrap();
        let view = source.over_periods(periods).unwrap();
        hb(&view, 3, ORDER * periods)
    }
    fn means(self, f: Family, orbit: &HbOrbit) -> [f64; 4] {
        let mut result = [0.0; 4];
        for sample in 0..2048 {
            let theta = TAU * sample as f64 / 2048.0;
            let x = core::array::from_fn(|i| orbit.sample(i, theta));
            let powers = self.physical(f, theta / orbit.omega, &x).1;
            for i in 0..4 {
                result[i] += powers[i] / 2048.0;
            }
        }
        result
    }
}
fn accounts() -> [ParametricPowerObjective; 8] {
    [
        ParametricPowerObjective::MeanPumpPower,
        ParametricPowerObjective::MeanPortPower,
        ParametricPowerObjective::MeanResistancePower,
        ParametricPowerObjective::MeanNonlinearPower,
        ParametricPowerObjective::PumpEnergyPerCycle,
        ParametricPowerObjective::PortEnergyPerCycle,
        ParametricPowerObjective::ResistanceEnergyPerCycle,
        ParametricPowerObjective::NonlinearEnergyPerCycle,
    ]
}

#[test]
fn complete_loss_pump_storage_and_clock_direction_matches_original_neighboring_physics() {
    let (f, p, l, lp) = (
        Family::source(),
        Family::slope(),
        Loss::source(),
        Loss::slope(),
    );
    let (system, input, pump) = (f.system(), f.input(), f.pump());
    let flow = |x: &[f64], e: &[f64], out: &mut [f64]| l.flow(x, e, out);
    let tangent = |x: &[f64], e: &[f64], dx: &[f64], de: &[f64], out: &mut [f64]| {
        l.tangent(x, e, dx, de, out)
    };
    let dp = |x: &[f64], e: &[f64], out: &mut [f64]| l.partial(lp, x, e, out);
    let hessian_calls = Cell::new(0);
    let hessian = |x: &[f64], v: &[f64], out: &mut [f64]| {
        hessian_calls.set(hessian_calls.get() + 1);
        Energy(f).apply(x, v, out)
    };
    let source = ParametricHamiltonian::new_analytic(&system, &input, &pump, &hessian)
        .unwrap()
        .with_dissipation(NonlinearDissipation {
            flow: &flow,
            tangent: &tangent,
            max_tangent_actions: 3,
        })
        .unwrap();
    let orbit = hb(&source, 3, ORDER);
    let mut prepared = source
        .prepare_sensitivity(&orbit, controls(), &scales(3, f.omega), || false)
        .unwrap();
    let (jp, rp, gp) = p.matrices();
    let (up, qp) = (f.input_partial(p), f.pump_partial(p));
    let hp = |x: &[f64], out: &mut [f64]| {
        out.copy_from_slice(&[
            p.k * x[0] + p.beta * x[0].powi(3) + p.mu * (x[0] - x[2]),
            0.0,
            p.mu * (x[2] - x[0]),
        ]);
        true
    };
    let direction = || ParametricParameterDirection {
        physical: PhsForcedParameterDirection {
            structure: PhsStructureParameterDirection {
                interconnection: Some(&jp),
                dissipation: Some(&rp),
                ports: Some(&gp),
            },
            storage: Some(&hp),
            input: Some(&up),
        },
        modulation: Some(&qp),
    };
    let derivative = prepared
        .solve_parameter_with_dissipation(chart(p.omega), direction(), Some(&dp), || false)
        .unwrap();
    let epsilon = 2e-5;
    let (fl, fh) = (f.shifted(p, -epsilon), f.shifted(p, epsilon));
    let (ll, lh) = (l.shifted(lp, -epsilon), l.shifted(lp, epsilon));
    let (lo, hi) = (ll.solve(fl, 1), lh.solve(fh, 1));
    for ((actual, a), b) in derivative
        .coeffs
        .iter()
        .flatten()
        .zip(lo.coeffs.iter().flatten())
        .zip(hi.coeffs.iter().flatten())
    {
        assert!((*actual - (*b - *a).scale(0.5 / epsilon)).abs() < 3e-7);
    }
    // Original independent time-domain equations check the same derivative.
    let time_lo = time_coefficients(fl.omega, 1024, ORDER, |t, x| ll.physical(fl, t, x).0);
    let time_hi = time_coefficients(fh.omega, 1024, ORDER, |t, x| lh.physical(fh, t, x).0);
    for ((actual, a), b) in derivative
        .coeffs
        .iter()
        .flatten()
        .zip(time_lo.iter().flatten())
        .zip(time_hi.iter().flatten())
    {
        assert!((*actual - (*b - *a).scale(0.5 / epsilon)).abs() < 3e-6);
    }
    let (base, low, high) = (l.means(f, &orbit), ll.means(fl, &lo), lh.means(fh, &hi));
    assert!(base[3] > 1e-5 && orbit.coeffs[2][1].abs() > 1e-3);
    close(base[0] + base[1] - base[2] - base[3], 0.0, 3e-11);
    let mut slopes = [0.0; 8];
    for (index, objective) in accounts().into_iter().enumerate() {
        let cycle = index >= 4;
        let mut adjoint = prepared
            .prepare_power_adjoint(objective, objectives(), || false)
            .unwrap();
        close(
            adjoint.value(),
            base[index % 4] * if cycle { TAU / f.omega } else { 1.0 },
            2e-11,
        );
        let calls = hessian_calls.get();
        let result = adjoint
            .gradient_parameter_with_dissipation(chart(p.omega), direction(), Some(&dp), || false)
            .unwrap();
        assert_eq!(
            hessian_calls.get(),
            calls,
            "parameter query must reuse the original tangent/factor"
        );
        let expected = if cycle {
            TAU * (high[index % 4] / fh.omega - low[index % 4] / fl.omega) / (2.0 * epsilon)
        } else {
            (high[index % 4] - low[index % 4]) / (2.0 * epsilon)
        };
        close(result.derivative, expected, 3e-7);
        slopes[index] = result.derivative;
    }
    close(slopes[0] + slopes[1] - slopes[2] - slopes[3], 0.0, 3e-9);
    close(slopes[4] + slopes[5] - slopes[6] - slopes[7], 0.0, 3e-9);
}

#[test]
fn fixed_loss_law_keeps_effort_chain_and_response_period_derivative() {
    let (f, l) = (Family::source(), Loss::source());
    let p = Family {
        k: 0.2,
        pump: 0.11,
        omega: 0.09,
        ..Default::default()
    };
    let (system, input, pump) = (f.system(), f.input(), f.pump());
    let flow = |x: &[f64], e: &[f64], out: &mut [f64]| l.flow(x, e, out);
    let tangent = |x: &[f64], e: &[f64], dx: &[f64], de: &[f64], out: &mut [f64]| {
        l.tangent(x, e, dx, de, out)
    };
    let source = ParametricHamiltonian::new(&system, &input, &pump)
        .unwrap()
        .with_dissipation(NonlinearDissipation {
            flow: &flow,
            tangent: &tangent,
            max_tangent_actions: 3,
        })
        .unwrap();
    let view = source.over_periods(2).unwrap();
    let orbit = hb(&view, 3, 2 * ORDER);
    let mut prepared = view
        .prepare_sensitivity(&orbit, controls(), &scales(3, f.omega / 2.0), || false)
        .unwrap();
    let hp = |x: &[f64], out: &mut [f64]| {
        out.copy_from_slice(&[p.k * x[0], 0.0, 0.0]);
        true
    };
    let qp = f.pump_partial(p);
    let direction = || ParametricParameterDirection {
        physical: PhsForcedParameterDirection {
            storage: Some(&hp),
            ..empty_direction().physical
        },
        modulation: Some(&qp),
    };
    // The old method explicitly fixes the loss law, not the effort seen by it.
    let derivative = prepared
        .solve_parameter(chart(p.omega / 2.0), direction(), || false)
        .unwrap();
    let epsilon = 2e-5;
    let (fl, fh) = (f.shifted(p, -epsilon), f.shifted(p, epsilon));
    let (lo, hi) = (l.solve(fl, 2), l.solve(fh, 2));
    for ((actual, a), b) in derivative
        .coeffs
        .iter()
        .flatten()
        .zip(lo.coeffs.iter().flatten())
        .zip(hi.coeffs.iter().flatten())
    {
        assert!((*actual - (*b - *a).scale(0.5 / epsilon)).abs() < 3e-7);
    }
    let (low, high) = (l.means(fl, &lo), l.means(fh, &hi));
    for (cycle, objective) in [
        (false, ParametricPowerObjective::MeanNonlinearPower),
        (true, ParametricPowerObjective::NonlinearEnergyPerCycle),
    ] {
        let mut adjoint = prepared
            .prepare_power_adjoint(objective, objectives(), || false)
            .unwrap();
        let result = adjoint
            .gradient_parameter(chart(p.omega / 2.0), direction(), || false)
            .unwrap();
        let expected = if cycle {
            2.0 * TAU * (high[3] / fh.omega - low[3] / fl.omega) / (2.0 * epsilon)
        } else {
            (high[3] - low[3]) / (2.0 * epsilon)
        };
        close(result.derivative, expected, 2e-7);
    }
}

#[test]
fn explicit_loss_only_directions_refuse_incomplete_writes_and_reuse_after_cancellation() {
    let (f, l, lp) = (Family::source(), Loss::source(), Loss::slope());
    let (system, input, pump) = (f.system(), f.input(), f.pump());
    let flow = |x: &[f64], e: &[f64], out: &mut [f64]| l.flow(x, e, out);
    let tangent = |x: &[f64], e: &[f64], dx: &[f64], de: &[f64], out: &mut [f64]| {
        l.tangent(x, e, dx, de, out)
    };
    let dp = |x: &[f64], e: &[f64], out: &mut [f64]| l.partial(lp, x, e, out);
    let source = ParametricHamiltonian::new(&system, &input, &pump)
        .unwrap()
        .with_dissipation(NonlinearDissipation {
            flow: &flow,
            tangent: &tangent,
            max_tangent_actions: 3,
        })
        .unwrap();
    let orbit = hb(&source, 3, ORDER);
    let mut prepared = source
        .prepare_sensitivity(&orbit, controls(), &scales(3, f.omega), || false)
        .unwrap();
    let bad = |_: &[f64], _: &[f64], out: &mut [f64]| {
        out[0] = 0.0;
        true
    };
    assert!(matches!(
        prepared.solve_parameter_with_dissipation(
            chart(0.0),
            empty_direction(),
            Some(&bad),
            || false
        ),
        Err(ParametricError::Dissipation(
            LinearizationError::DissipationParameterRefused { .. }
        ))
    ));
    let mut adjoint = prepared
        .prepare_power_adjoint(
            ParametricPowerObjective::NonlinearEnergyPerCycle,
            objectives(),
            || false,
        )
        .unwrap();
    let baseline = adjoint
        .gradient_parameter_with_dissipation(chart(0.0), empty_direction(), Some(&dp), || false)
        .unwrap();
    let epsilon = 2e-5;
    let (ll, lh) = (l.shifted(lp, -epsilon), l.shifted(lp, epsilon));
    let expected = TAU / f.omega
        * (lh.means(f, &lh.solve(f, 1))[3] - ll.means(f, &ll.solve(f, 1))[3])
        / (2.0 * epsilon);
    close(baseline.derivative, expected, 2e-7);
    let mut polls = 0;
    assert!(matches!(
        adjoint.gradient_parameter_with_dissipation(
            chart(0.0),
            empty_direction(),
            Some(&dp),
            || {
                polls += 1;
                polls > 25
            }
        ),
        Err(ParametricError::Sensitivity(
            ForcedSensitivityError::Cancelled
        ))
    ));
    assert_eq!(
        baseline,
        adjoint
            .gradient_parameter_with_dissipation(chart(0.0), empty_direction(), Some(&dp), || false)
            .unwrap()
    );
    assert!(matches!(
        adjoint.gradient_parameter_with_dissipation(
            chart(0.0),
            empty_direction(),
            Some(&bad),
            || false
        ),
        Err(ParametricError::Dissipation(
            LinearizationError::DissipationParameterRefused { .. }
        ))
    ));
    let nondissipative = ParametricHamiltonian::new(&system, &input, &pump).unwrap();
    let original = hb(&nondissipative, 3, ORDER);
    let mut prepared = nondissipative
        .prepare_sensitivity(&original, controls(), &scales(3, f.omega), || false)
        .unwrap();
    assert!(matches!(
        prepared
            .solve_parameter_with_dissipation(chart(0.0), empty_direction(), Some(&dp), || false),
        Err(ParametricError::Invalid(_))
    ));
}
