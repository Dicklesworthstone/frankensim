use super::*;
use crate::phs_orbit::autonomous_problem;
use fs_math::{c64::C64, det};
use fs_orbit::{HbAnchor, HbBudget, solve_hb_backbone_coordinate_scaled_seeded};
use fs_phs::{ContinuousStorageTangent, Storage};
use std::{cell::Cell, rc::Rc};

// Original momentum/displacement mechanics. In the coupled case the nonlinear
// first mass and second mass exchange energy through a physical elastic spring.
// Constant physical forces move the equilibrium and generate even harmonics.
#[derive(Clone)]
struct Mechanics {
    coupled: bool,
    native: bool,
    beta: f64,
    masses: [f64; 2],
    springs: [f64; 2],
    link: f64,
    j: f64,
    g: [f64; 2],
    bias: [f64; 2],
    gradient_calls: Rc<Cell<usize>>,
    tangent_calls: Rc<Cell<usize>>,
}
impl Mechanics {
    fn new(coupled: bool) -> Self {
        Self {
            coupled,
            native: true,
            beta: 8.0,
            masses: [1.3, 0.8],
            springs: [4.0, 11.0],
            link: if coupled { 0.7 } else { 0.0 },
            j: 0.9,
            g: [0.9, 1.1],
            bias: if coupled { [0.12, -0.09] } else { [0.0; 2] },
            gradient_calls: Rc::default(),
            tangent_calls: Rc::default(),
        }
    }
    fn n(&self) -> usize {
        if self.coupled { 4 } else { 2 }
    }
    fn shifted(&self, p: f64) -> Self {
        let mut changed = self.clone();
        changed.beta += 0.7 * p;
        changed.masses[0] += 0.08 * p;
        changed.masses[1] -= 0.03 * p;
        changed.springs[0] += 0.11 * p;
        changed.springs[1] += 0.17 * p;
        changed.link += 0.04 * p;
        changed.j += 0.06 * p;
        changed.g[0] += 0.05 * p;
        changed.g[1] -= 0.07 * p;
        changed.bias[0] += 0.09 * p;
        changed.bias[1] -= 0.03 * p;
        changed
    }
    fn system(&self) -> PortHamiltonian {
        let (n, m) = (self.n(), self.n() / 2);
        let mut j = vec![0.0; n * n];
        let mut g = vec![0.0; n * m];
        for i in 0..m {
            j[(2 * i) * n + 2 * i + 1] = self.j;
            j[(2 * i + 1) * n + 2 * i] = -self.j;
            g[(2 * i + 1) * m + i] = self.g[i];
        }
        PortHamiltonian::new(n, m, j, vec![0.0; n * n], g, Box::new(self.clone())).unwrap()
    }
    fn scales(&self) -> HbScales {
        HbScales {
            states: [0.4, 1.0, 0.04, 0.2][..self.n()].to_vec(),
            equations: [1.0, 3.0, 0.3, 1.0][..self.n()].to_vec(),
            angular_frequency: 2.0,
        }
    }
    fn orbit(&self, amplitude: f64, prior: Option<&HbOrbit>) -> HbOrbit {
        let (n, h, coordinate) = (self.n(), 11, if self.coupled { 2 } else { 0 });
        let system = self.system();
        let problem = autonomous_problem(&system, &self.bias[..n / 2]).unwrap();
        let anchor = HbAnchor::Backbone {
            amplitude,
            omega_guess: 2.0,
        };
        let seed = prior.map_or_else(
            || {
                let mut seed = vec![0.0; n * (2 * h + 1) + 1];
                let a0 = if self.coupled {
                    amplitude
                        * (self.springs[1] + self.link - self.masses[1] * (2.0 / self.j).powi(2))
                        / self.link
                } else {
                    amplitude
                };
                seed[n] = 0.5 * a0;
                seed[2 * n + 1] = 0.5 * self.masses[0] * 2.0 / self.j * a0;
                if self.coupled {
                    seed[n + 2] = 0.5 * amplitude;
                    seed[2 * n + 3] = 0.5 * self.masses[1] * 2.0 / self.j * amplitude;
                }
                *seed.last_mut().unwrap() = 2.0;
                seed
            },
            |previous| previous.packed_seed(&anchor),
        );
        solve_hb_backbone_coordinate_scaled_seeded(
            &problem,
            coordinate,
            amplitude,
            seed,
            &HbBudget {
                harmonics: h,
                max_newton: 60,
                tolerance: 1e-12,
                deflate_equilibrium: false,
            },
            &self.scales(),
        )
        .unwrap()
        .orbit
    }
}
impl Storage for Mechanics {
    fn hamiltonian(&self, x: &[f64]) -> f64 {
        let mut energy = 0.25 * self.beta * x[0].powi(4);
        for i in 0..self.n() / 2 {
            energy += 0.5 * self.springs[i] * x[2 * i].powi(2)
                + 0.5 * x[2 * i + 1].powi(2) / self.masses[i];
        }
        if self.coupled {
            energy += 0.5 * self.link * (x[0] - x[2]).powi(2);
        }
        energy
    }
    fn gradient(&self, x: &[f64], out: &mut [f64]) {
        self.gradient_calls.set(self.gradient_calls.get() + 1);
        for i in 0..self.n() / 2 {
            out[2 * i] = self.springs[i] * x[2 * i];
            out[2 * i + 1] = x[2 * i + 1] / self.masses[i];
        }
        out[0] += self.beta * x[0].powi(3);
        if self.coupled {
            out[0] += self.link * (x[0] - x[2]);
            out[2] += self.link * (x[2] - x[0]);
        }
    }
    fn continuous_tangent(&self) -> Option<&dyn ContinuousStorageTangent> {
        self.native.then_some(self)
    }
}
impl ContinuousStorageTangent for Mechanics {
    fn apply(&self, x: &[f64], direction: &[f64], out: &mut [f64]) -> bool {
        self.tangent_calls.set(self.tangent_calls.get() + 1);
        for i in 0..self.n() / 2 {
            out[2 * i] = self.springs[i] * direction[2 * i];
            out[2 * i + 1] = direction[2 * i + 1] / self.masses[i];
        }
        out[0] += 3.0 * self.beta * x[0] * x[0] * direction[0];
        if self.coupled {
            out[0] += self.link * (direction[0] - direction[2]);
            out[2] += self.link * (direction[2] - direction[0]);
        }
        true
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
        maximum_condition_estimate: 1e10,
        maximum_relative_residual: 1e-10,
        maximum_derivative_residual: 1e-10,
    }
}
fn direction(amplitude: f64) -> BackboneSensitivityDirection {
    BackboneSensitivityDirection {
        parameter_scale: 1.0,
        first_harmonic_amplitude_derivative: amplitude,
    }
}
fn reference(orbit: &HbOrbit) -> Vec<f64> {
    orbit.packed_seed(&HbAnchor::Autonomous {
        omega_guess: orbit.omega,
    })
}
// Use the same complete Parseval phase plane as the derivative. Independently
// solved backbone phases are free; one-coordinate phase alignment is inadequate.
fn align(base: &HbOrbit, trial: &mut HbOrbit, scales: &HbScales) {
    let mut phase = 0.0;
    for _ in 0..8 {
        let (mut value, mut derivative) = (0.0, 0.0);
        for (state, rows) in base.coeffs.iter().enumerate() {
            for (h, &z) in rows.iter().enumerate().skip(1) {
                let a = z.scale(1.0 / scales.states[state]);
                let b = trial.coeffs[state][h].scale(1.0 / scales.states[state])
                    * C64::new(det::cos(h as f64 * phase), det::sin(h as f64 * phase));
                value += h as f64 * (-a.im * b.re + a.re * b.im);
                derivative += (h * h) as f64 * (a.re * b.re + a.im * b.im);
            }
        }
        phase -= value / derivative;
    }
    for row in &mut trial.coeffs {
        for (h, z) in row.iter_mut().enumerate().skip(1) {
            *z = *z * C64::new(det::cos(h as f64 * phase), det::sin(h as f64 * phase));
        }
    }
}
fn compare(result: &AutonomousOrbitSensitivity, lo: &HbOrbit, hi: &HbOrbit, h: f64) {
    assert!((result.angular_frequency_derivative - (hi.omega - lo.omega) / (2.0 * h)).abs() < 2e-7);
    for state in 0..result.coeffs.len() {
        for harmonic in 0..result.coeffs[state].len() {
            let measured = (hi.coeffs[state][harmonic] - lo.coeffs[state][harmonic]).scale(0.5 / h);
            assert!(
                (result.coeffs[state][harmonic] - measured).abs() < 4e-7,
                "state {state}, harmonic {harmonic}: {:?} versus {measured:?}",
                result.coeffs[state][harmonic]
            );
        }
    }
    assert!(result.derivative_residual < 1e-10);
}

#[test]
fn duffing_material_and_amplitude_slopes_match_original_orbits_and_continuous_period() {
    let model = Mechanics::new(false);
    let system = model.system();
    let base = model.orbit(0.4, None);
    let mut prepared = PreparedPhsBackboneSensitivity::new(
        &system,
        &[0.0],
        &base,
        BackboneSensitivityAnchor {
            coordinate: 0,
            first_harmonic_amplitude: 0.4,
        },
        &reference(&base),
        controls(),
        &model.scales(),
        || false,
    )
    .unwrap();
    assert!(base.coeffs[0][3].abs() > 1e-4);
    for material in [true, false] {
        let tangent_calls = model.tangent_calls.get();
        let partial = |x: &[f64], out: &mut [f64]| {
            out[0] = x[0].powi(3);
            out[1] = 0.0;
            true
        };
        let result = prepared
            .solve_parameter(
                direction(if material { 0.0 } else { 1.0 }),
                PhsBackboneParameterDirection {
                    storage: if material { Some(&partial) } else { None },
                    ..Default::default()
                },
                || false,
            )
            .unwrap();
        // Queries reapply original gradients/parameter columns, not Hessians.
        assert_eq!(model.tangent_calls.get(), tangent_calls);
        let h = 2e-5;
        let (mut low, mut high) = (model.clone(), model.clone());
        if material {
            low.beta -= h;
            high.beta += h;
        }
        let mut lo = low.orbit(0.4 - if material { 0.0 } else { h }, Some(&base));
        let mut hi = high.orbit(0.4 + if material { 0.0 } else { h }, Some(&base));
        align(&base, &mut lo, &model.scales());
        align(&base, &mut hi, &model.scales());
        compare(&result, &lo, &hi, h);

        // Independent continuous energy-period quadrature: no HB residual or
        // Fourier projection is used in this physical period oracle.
        let (q, p) = (base.sample(0, 0.0), base.sample(1, 0.0));
        let sample = |state: usize| {
            result.coeffs[state][0].re
                + 2.0
                    * result.coeffs[state]
                        .iter()
                        .skip(1)
                        .map(|z| z.re)
                        .sum::<f64>()
        };
        let beta_p = if material { 1.0 } else { 0.0 };
        let energy = model.hamiltonian(&[q, p]);
        let ep = (4.0 * q + model.beta * q.powi(3)) * sample(0)
            + p * sample(1) / model.masses[0]
            + 0.25 * beta_p * q.powi(4);
        let a = (4.0 * energy / (4.0 + (16.0 + 4.0 * model.beta * energy).sqrt())).sqrt();
        let ap = (ep - 0.25 * beta_p * a.powi(4)) / (4.0 * a + model.beta * a.powi(3));
        let step = std::f64::consts::FRAC_PI_2 / 4096.0;
        let (mut integral, mut slope) = (0.0, 0.0);
        for j in 0..=4096 {
            let weight = if j == 0 || j == 4096 {
                1.0
            } else if j % 2 == 0 {
                2.0
            } else {
                4.0
            };
            let shape = 1.0 + (j as f64 * step).sin().powi(2);
            let square = 4.0 + 0.5 * model.beta * a * a * shape;
            integral += weight / square.sqrt();
            slope -= 0.5 * weight * (0.5 * beta_p * a * a + model.beta * a * ap) * shape
                / square.powf(1.5);
        }
        let time_scale = model.masses[0].sqrt() / model.j;
        let period = time_scale * 4.0 * step * integral / 3.0;
        let period_p = time_scale * 4.0 * step * slope / 3.0;
        assert!(
            (result.angular_frequency_derivative
                + std::f64::consts::TAU * period_p / (period * period))
                .abs()
                < 2e-8
        );
    }
}

#[test]
fn complete_coupled_material_interconnection_port_and_bias_chain_matches_neighboring_systems() {
    let model = Mechanics::new(true);
    let system = model.system();
    let base = model.orbit(0.035, None);
    let mut prepared = PreparedPhsBackboneSensitivity::new(
        &system,
        &model.bias,
        &base,
        BackboneSensitivityAnchor {
            coordinate: 2,
            first_harmonic_amplitude: 0.035,
        },
        &reference(&base),
        controls(),
        &model.scales(),
        || false,
    )
    .unwrap();
    assert!(
        base.coeffs[0][2].abs() > 1e-5,
        "physical bias must retain even harmonics"
    );
    let mut jp = vec![0.0; 16];
    let mut gp = vec![0.0; 8];
    jp[1] = 0.06;
    jp[4] = -0.06;
    jp[11] = 0.06;
    jp[14] = -0.06;
    gp[2] = 0.05;
    gp[7] = -0.07;
    let storage = |x: &[f64], out: &mut [f64]| {
        out[0] = 0.11 * x[0] + 0.7 * x[0].powi(3) + 0.04 * (x[0] - x[2]);
        out[1] = -0.08 * x[1] / model.masses[0].powi(2);
        out[2] = 0.17 * x[2] + 0.04 * (x[2] - x[0]);
        out[3] = 0.03 * x[3] / model.masses[1].powi(2);
        true
    };
    let result = prepared
        .solve_parameter(
            direction(0.002),
            PhsBackboneParameterDirection {
                structure: PhsStructureParameterDirection {
                    interconnection: Some(&jp),
                    dissipation: None,
                    ports: Some(&gp),
                },
                storage: Some(&storage),
                bias: Some(&[0.09, -0.03]),
            },
            || false,
        )
        .unwrap();
    let h = 2e-5;
    let mut lo = model.shifted(-h).orbit(0.035 - 0.002 * h, Some(&base));
    let mut hi = model.shifted(h).orbit(0.035 + 0.002 * h, Some(&base));
    align(&base, &mut lo, &model.scales());
    align(&base, &mut hi, &model.scales());
    compare(&result, &lo, &hi, h);
}

#[test]
fn exact_tangent_admission_incompatible_loss_and_bad_directions_leave_factor_reusable() {
    let mut model = Mechanics::new(false);
    let base = model.orbit(0.4, None);
    model.native = false;
    let system = model.system();
    let anchor = BackboneSensitivityAnchor {
        coordinate: 0,
        first_harmonic_amplitude: 0.4,
    };
    assert!(matches!(
        PreparedPhsBackboneSensitivity::new(
            &system,
            &[0.0],
            &base,
            anchor,
            &reference(&base),
            controls(),
            &model.scales(),
            || false,
        ),
        Err(PhsBackboneSensitivityError::Sensitivity(
            ForcedSensitivityError::TangentUnavailable
        ))
    ));
    let hessian = |x: &[f64], d: &[f64], out: &mut [f64]| model.apply(x, d, out);
    let before = model.gradient_calls.get();
    let mut tiny = controls();
    tiny.maximum_payload_bytes = 1;
    assert!(matches!(
        PreparedPhsBackboneSensitivity::new_analytic(
            &system,
            &[0.0],
            &base,
            anchor,
            &reference(&base),
            tiny,
            &model.scales(),
            &hessian,
            || false,
        ),
        Err(PhsBackboneSensitivityError::Sensitivity(
            ForcedSensitivityError::Budget { .. }
        ))
    ));
    assert_eq!(before, model.gradient_calls.get());
    let mut prepared = PreparedPhsBackboneSensitivity::new_analytic(
        &system,
        &[0.0],
        &base,
        anchor,
        &reference(&base),
        controls(),
        &model.scales(),
        &hessian,
        || false,
    )
    .unwrap();
    let expected = prepared
        .solve_parameter(direction(1.0), Default::default(), || false)
        .unwrap();
    let damping = [0.0, 0.0, 0.0, 1.0];
    assert!(matches!(
        prepared.solve_parameter(
            direction(0.0),
            PhsBackboneParameterDirection {
                structure: PhsStructureParameterDirection {
                    dissipation: Some(&damping),
                    ..Default::default()
                },
                ..Default::default()
            },
            || false
        ),
        Err(PhsBackboneSensitivityError::Sensitivity(
            ForcedSensitivityError::DerivativeUnresolved { .. }
        ))
    ));
    let incomplete = |_x: &[f64], out: &mut [f64]| {
        out[0] = 0.0;
        true
    };
    assert!(matches!(
        prepared.solve_parameter(
            direction(0.0),
            PhsBackboneParameterDirection {
                storage: Some(&incomplete),
                ..Default::default()
            },
            || false
        ),
        Err(PhsBackboneSensitivityError::Sensitivity(
            ForcedSensitivityError::IslandDirectionRefused { .. }
        ))
    ));
    let bad_j = [0.0, 1.0, 1.0, 0.0];
    assert!(matches!(
        prepared.solve_parameter(
            direction(0.0),
            PhsBackboneParameterDirection {
                structure: PhsStructureParameterDirection {
                    interconnection: Some(&bad_j),
                    ..Default::default()
                },
                ..Default::default()
            },
            || false
        ),
        Err(PhsBackboneSensitivityError::Physical(_))
    ));
    assert!(matches!(
        prepared.solve_parameter(
            direction(0.0),
            PhsBackboneParameterDirection {
                bias: Some(&[]),
                ..Default::default()
            },
            || false
        ),
        Err(PhsBackboneSensitivityError::Source(_))
    ));
    assert!(matches!(
        prepared.solve_parameter(direction(1.0), Default::default(), || true),
        Err(PhsBackboneSensitivityError::Sensitivity(
            ForcedSensitivityError::Cancelled
        ))
    ));
    let recovered = prepared
        .solve_parameter(direction(1.0), Default::default(), || false)
        .unwrap();
    assert_eq!(expected.coeffs, recovered.coeffs);
    assert_eq!(
        expected.angular_frequency_derivative,
        recovered.angular_frequency_derivative
    );
}
