//! G1/G3: original conservative periods, coupled harmonics and physical units.

use fs_math::{c64::C64, det};
use fs_orbit::sensitivity::{
    AutonomousOrbitSensitivity, AutonomousSensitivityPoint, BackboneSensitivityAnchor,
    BackboneSensitivityDirection, ForcedSensitivityControls, ForcedSensitivityError,
    PreparedBackboneSensitivity,
};
use fs_orbit::{
    HbAnchor, HbBudget, HbOrbit, HbScales, IslandStateJacobian, OrbitProblem,
    solve_hb_backbone_coordinate_scaled_seeded,
};
use std::cell::Cell;

// A conservative Duffing family with a separately relaxing static first state.
// Dense equation multiplication keeps the full frequency port and state tangent
// relevant. No first-coordinate or diagonal-mass shortcut can solve this chart.
const B: [f64; 9] = [2.0, 0.3, 0.1, 0.0, 3.0, 0.2, 0.4, 0.0, 5.0];
struct Model {
    beta: f64,
    tau: f64,
    units: [f64; 3],
    equations: [f64; 3],
    source_calls: Cell<usize>,
    tangent_calls: Cell<usize>,
    exact: bool,
}
impl Model {
    fn new(converted: bool) -> Self {
        Self {
            beta: 8.0,
            tau: 1.3,
            units: if converted {
                [1e4, 1e-9, 2e3]
            } else {
                [1.0; 3]
            },
            equations: if converted {
                [1e-8, 20.0, 2e4]
            } else {
                [1.0; 3]
            },
            source_calls: Cell::new(0),
            tangent_calls: Cell::new(0),
            exact: true,
        }
    }
    fn scales(&self) -> HbScales {
        HbScales {
            states: vec![4.0 * self.units[0], 0.4 * self.units[1], self.units[2]],
            equations: vec![
                2.0 * self.equations[0],
                3.0 * self.equations[1],
                10.0 * self.equations[2],
            ],
            angular_frequency: 2.0 / self.tau,
        }
    }
    fn rows(&self, values: [f64; 3], out: &mut [f64]) {
        for i in 0..3 {
            out[i] = self.equations[i] * (0..3).map(|j| B[3 * i + j] * values[j]).sum::<f64>();
        }
    }
    fn jacobian(&self, _: f64, x: &[f64], out: &mut [f64]) -> bool {
        self.tangent_calls.set(self.tangent_calls.get() + 1);
        let q = x[1] / self.units[1];
        let jac = [
            -7.0,
            0.0,
            0.0,
            0.0,
            0.0,
            1.0,
            0.0,
            -4.0 - 3.0 * self.beta * q * q,
            0.0,
        ];
        for row in 0..3 {
            for col in 0..3 {
                out[3 * row + col] = self.equations[row] / self.units[col]
                    * (0..3)
                        .map(|j| B[3 * row + j] * jac[3 * j + col])
                        .sum::<f64>();
            }
        }
        true
    }
    fn frequency_column(&self, point: AutonomousSensitivityPoint<'_>, out: &mut [f64]) -> bool {
        out.fill(0.0);
        for k in 1..=point.harmonics {
            let re = 3 + (k - 1) * 6;
            let im = re + 3;
            for row in 0..3 {
                for col in 0..3 {
                    let v = k as f64 * self.tau * self.equations[row] * B[3 * row + col]
                        / self.units[col];
                    out[re + row] -= v * point.coefficients[im + col];
                    out[im + row] += v * point.coefficients[re + col];
                }
            }
        }
        true
    }
    fn orbit(&self, amplitude: f64, start: Option<&HbOrbit>) -> HbOrbit {
        let order = 11;
        let anchor = HbAnchor::Backbone {
            amplitude: amplitude * self.units[1],
            omega_guess: 2.1 / self.tau,
        };
        let seed = start.map_or_else(
            || {
                let mut seed = vec![0.0; 3 * (2 * order + 1) + 1];
                seed[0] = (4.0 + 0.02 * self.beta) * self.units[0];
                let phase: f64 = 0.37;
                seed[4] = 0.5 * amplitude * phase.cos() * self.units[1];
                seed[7] = 0.5 * amplitude * phase.sin() * self.units[1];
                seed[5] = -0.5 * 2.1 * amplitude * phase.sin() * self.units[2];
                seed[8] = 0.5 * 2.1 * amplitude * phase.cos() * self.units[2];
                *seed.last_mut().unwrap() = 2.1 / self.tau;
                seed
            },
            |orbit| orbit.packed_seed(&anchor),
        );
        solve_hb_backbone_coordinate_scaled_seeded(
            self,
            1,
            amplitude * self.units[1],
            seed,
            &HbBudget {
                harmonics: order,
                max_newton: 60,
                tolerance: 1e-12,
                deflate_equilibrium: false,
            },
            &self.scales(),
        )
        .unwrap()
        .orbit
    }
    fn prepare(&self, orbit: &HbOrbit) -> PreparedBackboneSensitivity {
        PreparedBackboneSensitivity::new(
            self,
            orbit,
            BackboneSensitivityAnchor {
                coordinate: 1,
                first_harmonic_amplitude: 0.4 * self.units[1],
            },
            &orbit.packed_seed(&HbAnchor::Autonomous {
                omega_guess: orbit.omega,
            }),
            controls(),
            &self.scales(),
            |p, out| self.frequency_column(p, out),
            || false,
        )
        .unwrap()
    }
}
impl OrbitProblem for Model {
    fn dim(&self) -> usize {
        3
    }
    fn autonomous(&self) -> bool {
        true
    }
    fn island(&self, _: f64, state: &[f64], out: &mut [f64]) {
        self.source_calls.set(self.source_calls.get() + 1);
        let x = std::array::from_fn::<_, 3, _>(|i| state[i] / self.units[i]);
        self.rows(
            [
                -7.0 * (x[0] - 4.0 - 0.02 * self.beta),
                x[2],
                -4.0 * x[1] - self.beta * x[1].powi(3),
            ],
            out,
        );
    }
    fn port(&self, s: C64) -> Vec<C64> {
        (0..9)
            .map(|i| s.scale(self.tau * self.equations[i / 3] * B[i] / self.units[i % 3]))
            .collect()
    }
    fn island_state_jacobian(&self) -> Option<IslandStateJacobian<Self>> {
        self.exact.then_some(Self::jacobian)
    }
}
fn controls() -> ForcedSensitivityControls {
    ForcedSensitivityControls {
        samples_per_period: 128,
        maximum_unknowns: 200,
        maximum_payload_bytes: 32_000_000,
        maximum_work: 1_000_000_000,
        maximum_base_residual: 1e-10,
        maximum_tangent_sampling_defect: 1e-9,
        maximum_condition_estimate: 1e9,
        maximum_relative_residual: 1e-11,
        maximum_derivative_residual: 1e-10,
    }
}
fn direction(amplitude: f64) -> BackboneSensitivityDirection {
    BackboneSensitivityDirection {
        parameter_scale: 1.0,
        first_harmonic_amplitude_derivative: amplitude,
    }
}
fn derivative(
    model: &Model,
    prepared: &mut PreparedBackboneSensitivity,
    material: bool,
) -> AutonomousOrbitSensitivity {
    prepared
        .solve_island_direction(
            direction(if material { 0.0 } else { model.units[1] }),
            |_, state, out| {
                let q = state[1] / model.units[1];
                model.rows(
                    if material {
                        [0.14, 0.0, -q.powi(3)]
                    } else {
                        [0.0; 3]
                    },
                    out,
                );
                true
            },
            || false,
        )
        .unwrap()
}

// Match the fixed original full-waveform phase hyperplane, not one arbitrarily
// selected coefficient's phase. Its Parseval factors cancel in this equation.
fn align(reference: &HbOrbit, trial: &mut HbOrbit, scales: &HbScales) {
    let mut phase = 0.0;
    for _ in 0..6 {
        let (mut f, mut df) = (0.0, 0.0);
        for j in 0..3 {
            for k in 1..reference.coeffs[j].len() {
                let a = reference.coeffs[j][k].scale(1.0 / scales.states[j]);
                let b = trial.coeffs[j][k].scale(1.0 / scales.states[j])
                    * C64::new(det::cos(k as f64 * phase), det::sin(k as f64 * phase));
                f += k as f64 * (-a.im * b.re + a.re * b.im);
                df += (k * k) as f64 * (a.re * b.re + a.im * b.im);
            }
        }
        phase -= f / df;
    }
    for row in &mut trial.coeffs {
        for (k, value) in row.iter_mut().enumerate().skip(1) {
            *value = *value * C64::new(det::cos(k as f64 * phase), det::sin(k as f64 * phase));
        }
    }
}
fn sampled_derivative(result: &AutonomousOrbitSensitivity, state: usize) -> f64 {
    result.coeffs[state][0].re
        + 2.0
            * result.coeffs[state]
                .iter()
                .skip(1)
                .map(|z| z.re)
                .sum::<f64>()
}
// Independently differentiate the continuous conservative energy-period
// quadrature, using original displacement/velocity and their physical changes.
fn period_derivative(
    model: &Model,
    base: &HbOrbit,
    result: &AutonomousOrbitSensitivity,
    beta_p: f64,
) -> f64 {
    let q = base.sample(1, 0.0) / model.units[1];
    let v = base.sample(2, 0.0) / model.units[2];
    let qp = sampled_derivative(result, 1) / model.units[1];
    let vp = sampled_derivative(result, 2) / model.units[2];
    let energy = 2.0 * q * q + 0.25 * model.beta * q.powi(4) + 0.5 * v * v;
    let ep = (4.0 * q + model.beta * q.powi(3)) * qp + v * vp + 0.25 * beta_p * q.powi(4);
    let a = (4.0 * energy / (4.0 + (16.0 + 4.0 * model.beta * energy).sqrt())).sqrt();
    let ap = (ep - 0.25 * beta_p * a.powi(4)) / (4.0 * a + model.beta * a.powi(3));
    let panels = 4096;
    let h = std::f64::consts::FRAC_PI_2 / panels as f64;
    let (mut t, mut tp) = (0.0, 0.0);
    for j in 0..=panels {
        let weight = if j == 0 || j == panels {
            1.0
        } else if j % 2 == 0 {
            2.0
        } else {
            4.0
        };
        let shape = 1.0 + (j as f64 * h).sin().powi(2);
        let square = 4.0 + 0.5 * model.beta * a * a * shape;
        t += weight / square.sqrt();
        tp +=
            -0.5 * weight * (0.5 * beta_p * a * a + model.beta * a * ap) * shape / square.powf(1.5);
    }
    let period = 4.0 * h * t / 3.0;
    -std::f64::consts::TAU * (4.0 * h * tp / 3.0) / (model.tau * period * period)
}

#[test]
fn material_and_amplitude_derivatives_match_nonlinear_orbits_and_energy_period() {
    let model = Model::new(false);
    let base = model.orbit(0.4, None);
    let mut prepared = model.prepare(&base);
    assert!(base.coeffs[1][3].abs() > 1e-4);
    for material in [true, false] {
        let result = derivative(&model, &mut prepared, material);
        let h = 2e-5;
        let mut low = Model::new(false);
        let mut high = Model::new(false);
        if material {
            low.beta -= h;
            high.beta += h;
        }
        let mut lo = low.orbit(0.4 - if material { 0.0 } else { h }, Some(&base));
        let mut hi = high.orbit(0.4 + if material { 0.0 } else { h }, Some(&base));
        align(&base, &mut lo, &model.scales());
        align(&base, &mut hi, &model.scales());
        assert!(
            (result.angular_frequency_derivative - (hi.omega - lo.omega) / (2.0 * h)).abs() < 2e-8
        );
        assert!(
            (result.angular_frequency_derivative
                - period_derivative(&model, &base, &result, if material { 1.0 } else { 0.0 }))
            .abs()
                < 2e-8
        );
        for j in 0..3 {
            for k in 0..base.coeffs[j].len() {
                let numerical = (hi.coeffs[j][k] - lo.coeffs[j][k]).scale(0.5 / h);
                assert!(
                    (result.coeffs[j][k] - numerical).abs() < 8e-8,
                    "state {j} harmonic {k}"
                );
            }
        }
        assert!(result.derivative_residual < 1e-10);
    }
}

#[test]
fn complete_mass_time_port_and_coordinate_units_reuse_the_same_factor() {
    let reference = Model::new(false);
    let original = reference.orbit(0.4, None);
    let mut base_prepared = reference.prepare(&original);
    let expected = derivative(&reference, &mut base_prepared, true);
    let converted = Model::new(true);
    let mut base = converted.orbit(0.4, None);
    // Backbone solves have a free phase. Compare physical unit changes in the
    // SAME full-waveform gauge before freezing either sensitivity reference.
    let mut converted_reference = original.clone();
    for j in 0..3 {
        for value in &mut converted_reference.coeffs[j] {
            *value = value.scale(converted.units[j]);
        }
    }
    align(&converted_reference, &mut base, &converted.scales());
    let mut prepared = converted.prepare(&base);
    let calls = (converted.source_calls.get(), converted.tangent_calls.get());
    let parameter_unit = 1e-4;
    let actual = prepared
        .solve_island_direction(
            BackboneSensitivityDirection {
                parameter_scale: parameter_unit,
                first_harmonic_amplitude_derivative: 0.0,
            },
            |_, state, out| {
                converted.rows(
                    [
                        0.14 / parameter_unit,
                        0.0,
                        -(state[1] / converted.units[1]).powi(3) / parameter_unit,
                    ],
                    out,
                );
                true
            },
            || false,
        )
        .unwrap();
    assert!(
        (actual.angular_frequency_derivative * parameter_unit
            - expected.angular_frequency_derivative)
            .abs()
            < 2e-10
    );
    for j in 0..3 {
        for k in 0..base.coeffs[j].len() {
            assert!(
                (actual.coeffs[j][k].scale(parameter_unit / converted.units[j])
                    - expected.coeffs[j][k])
                    .abs()
                    < 2e-10,
                "state {j} harmonic {k}: actual {:?}, expected {:?}, base phases {:?} {:?}",
                actual.coeffs[j][k].scale(parameter_unit / converted.units[j]),
                expected.coeffs[j][k],
                base.coeffs[1][1].scale(1.0 / converted.units[1]),
                original.coeffs[1][1]
            );
        }
    }
    let time = prepared
        .solve(
            direction(0.0),
            |p, out| {
                converted.frequency_column(p, out);
                for value in out {
                    *value *= p.angular_frequency / converted.tau;
                }
                true
            },
            || false,
        )
        .unwrap();
    assert!((time.angular_frequency_derivative + base.omega / converted.tau).abs() < 2e-10);
    for j in 0..3 {
        assert!(
            time.coeffs[j]
                .iter()
                .all(|z| z.abs() / converted.units[j] < 2e-10)
        );
    }
    assert_eq!(
        calls,
        (converted.source_calls.get(), converted.tangent_calls.get())
    );
}

#[test]
fn dissipative_parameter_direction_refuses_and_factor_survives() {
    let model = Model::new(false);
    let base = model.orbit(0.4, None);
    let mut prepared = model.prepare(&base);
    let before = derivative(&model, &mut prepared, true);
    let failure = prepared.solve_island_direction(
        direction(0.0),
        |_, state, out| {
            model.rows([0.0, 0.0, -state[2]], out);
            true
        },
        || false,
    );
    assert!(
        matches!(
            failure,
            Err(ForcedSensitivityError::DerivativeUnresolved { .. })
        ),
        "{failure:?}"
    );
    let after = derivative(&model, &mut prepared, true);
    assert_eq!(
        before.angular_frequency_derivative,
        after.angular_frequency_derivative
    );
    assert_eq!(before.coeffs, after.coeffs);
}

#[test]
fn admission_declared_tangents_and_cancellation_are_transactional() {
    let mut model = Model::new(false);
    let base = model.orbit(0.4, None);
    let phase = base.packed_seed(&HbAnchor::Autonomous {
        omega_guess: base.omega,
    });
    let anchor = BackboneSensitivityAnchor {
        coordinate: 1,
        first_harmonic_amplitude: 0.4,
    };
    let calls = model.source_calls.get();
    let mut small = controls();
    small.maximum_payload_bytes = 1;
    assert!(matches!(
        PreparedBackboneSensitivity::new(
            &model,
            &base,
            anchor,
            &phase,
            small,
            &model.scales(),
            |p, out| model.frequency_column(p, out),
            || false
        ),
        Err(ForcedSensitivityError::Budget { .. })
    ));
    assert_eq!(calls, model.source_calls.get());
    assert!(matches!(
        PreparedBackboneSensitivity::new(
            &model,
            &base,
            anchor,
            &phase,
            controls(),
            &model.scales(),
            |p, out| model.frequency_column(p, out),
            || true
        ),
        Err(ForcedSensitivityError::Cancelled)
    ));
    assert!(matches!(
        PreparedBackboneSensitivity::new(
            &model,
            &base,
            BackboneSensitivityAnchor {
                first_harmonic_amplitude: 0.5,
                ..anchor
            },
            &phase,
            controls(),
            &model.scales(),
            |p, out| model.frequency_column(p, out),
            || false
        ),
        Err(ForcedSensitivityError::BaseUnresolved { .. })
    ));
    assert!(matches!(
        PreparedBackboneSensitivity::new(
            &model,
            &base,
            anchor,
            &phase,
            controls(),
            &model.scales(),
            |_, _| true,
            || false
        ),
        Err(ForcedSensitivityError::FrequencyDirectionRefused { .. })
    ));
    model.exact = false;
    assert!(matches!(
        PreparedBackboneSensitivity::new(
            &model,
            &base,
            anchor,
            &phase,
            controls(),
            &model.scales(),
            |p, out| model.frequency_column(p, out),
            || false
        ),
        Err(ForcedSensitivityError::TangentUnavailable)
    ));
    model.exact = true;
    let mut prepared = model.prepare(&base);
    assert!(matches!(
        prepared.solve(direction(0.0), |_, _| true, || false),
        Err(ForcedSensitivityError::DirectionRefused { .. })
    ));
    assert!(matches!(
        prepared.solve_island_direction(
            direction(1.0),
            |_, _, out| {
                out.fill(0.0);
                true
            },
            || true
        ),
        Err(ForcedSensitivityError::Cancelled)
    ));
    let mut polls = 0;
    assert!(matches!(
        prepared.solve_island_direction(
            direction(1.0),
            |_, _, out| {
                out.fill(0.0);
                true
            },
            || {
                polls += 1;
                polls == 9
            }
        ),
        Err(ForcedSensitivityError::Cancelled)
    ));
    derivative(&model, &mut prepared, false);
}
