//! G1/G3 full descriptor/material jets against original physical formulas.
//!
//! The small scalar jet below differentiates closed-form physical expressions,
//! independently of the production original-state implicit linear solves.
use fs_math::{c64::C64, det};
use fs_phs::circuit::{Branch, CircuitGraph, assemble_circuit};
use fs_phs::harmonic::{
    HarmonicBudget, HarmonicError, HarmonicFrequencyParameterResponse, HarmonicModel,
    HarmonicParameterDirection as Direction,
};
use fs_phs::{PortHamiltonian, QuadraticStorage, RelaxationBranch};

#[derive(Clone, Copy, Default)]
struct J {
    v: C64,
    w: C64,
    ww: C64,
    t: C64,
    wt: C64,
}
impl J {
    fn real(v: f64, t: f64) -> Self {
        Self {
            v: C64::from_re(v),
            t: C64::from_re(t),
            ..Self::default()
        }
    }
    fn frequency(w: f64) -> Self {
        Self {
            v: C64::from_re(w),
            w: C64::ONE,
            ..Self::default()
        }
    }
    fn complex(v: C64) -> Self {
        Self {
            v,
            ..Self::default()
        }
    }
    fn scale(self, a: f64) -> Self {
        Self {
            v: self.v.scale(a),
            w: self.w.scale(a),
            ww: self.ww.scale(a),
            t: self.t.scale(a),
            wt: self.wt.scale(a),
        }
    }
    fn conj(self) -> Self {
        Self {
            v: self.v.conj(),
            w: self.w.conj(),
            ww: self.ww.conj(),
            t: self.t.conj(),
            wt: self.wt.conj(),
        }
    }
    fn norm(self) -> Self {
        self.conj() * self
    }
    fn reciprocal(self) -> Self {
        let v = C64::ONE / self.v;
        let w = -(self.w * v) * v;
        let t = -(self.t * v) * v;
        let ww = -(self.ww * v) * v + (self.w * self.w * v * v * v).scale(2.0);
        let wt = -(self.wt * v) * v + (self.w * self.t * v * v * v).scale(2.0);
        Self { v, w, ww, t, wt }
    }
    fn exp(self) -> Self {
        let scale = det::exp(self.v.re);
        let v = C64::new(scale * det::cos(self.v.im), scale * det::sin(self.v.im));
        Self {
            v,
            w: v * self.w,
            ww: v * (self.ww + self.w * self.w),
            t: v * self.t,
            wt: v * (self.wt + self.w * self.t),
        }
    }
}
impl core::ops::Add for J {
    type Output = Self;
    fn add(self, b: Self) -> Self {
        Self {
            v: self.v + b.v,
            w: self.w + b.w,
            ww: self.ww + b.ww,
            t: self.t + b.t,
            wt: self.wt + b.wt,
        }
    }
}
impl core::ops::Neg for J {
    type Output = Self;
    fn neg(self) -> Self {
        self.scale(-1.0)
    }
}
impl core::ops::Sub for J {
    type Output = Self;
    fn sub(self, b: Self) -> Self {
        self + (-b)
    }
}
impl core::ops::Mul for J {
    type Output = Self;
    fn mul(self, b: Self) -> Self {
        Self {
            v: self.v * b.v,
            w: self.w * b.v + self.v * b.w,
            ww: self.ww * b.v + (self.w * b.w).scale(2.0) + self.v * b.ww,
            t: self.t * b.v + self.v * b.t,
            wt: self.wt * b.v + self.w * b.t + self.t * b.w + self.v * b.wt,
        }
    }
}
impl core::ops::Div for J {
    type Output = Self;
    fn div(self, b: Self) -> Self {
        self * b.reciprocal()
    }
}

fn budget() -> HarmonicBudget {
    HarmonicBudget {
        max_dimension: 64,
        max_dense_bytes: 2_000_000,
        max_dense_work: 4_000_000,
        relative_residual_tolerance: 1e-12,
    }
}
fn near(a: C64, e: C64, tol: f64) {
    assert!(
        (a - e).abs() <= tol * (1.0 + e.abs()),
        "actual {a:?}, expected {e:?}"
    );
}
fn real(a: f64, e: f64, tol: f64) {
    near(C64::from_re(a), C64::from_re(e), tol);
}
fn check(actual: [C64; 5], expected: J, tol: f64) {
    for (a, e) in
        actual
            .into_iter()
            .zip([expected.v, expected.w, expected.ww, expected.t, expected.wt])
    {
        near(a, e, tol);
    }
}
fn state(r: &HarmonicFrequencyParameterResponse, i: usize, expected: J, tol: f64) {
    let f = &r.frequency_response.frequency_response;
    check(
        [
            f.response.state[i],
            f.state_derivative[i],
            r.frequency_response.state_second_derivative[i],
            r.parameter_response.state_derivative[i],
            r.state_mixed_derivative[i],
        ],
        expected,
        tol,
    );
}
fn effort(r: &HarmonicFrequencyParameterResponse, i: usize, expected: J, tol: f64) {
    let f = &r.frequency_response.frequency_response;
    check(
        [
            f.response.effort[i],
            f.effort_derivative[i],
            r.frequency_response.effort_second_derivative[i],
            r.parameter_response.effort_derivative[i],
            r.effort_mixed_derivative[i],
        ],
        expected,
        tol,
    );
}
fn output(r: &HarmonicFrequencyParameterResponse, i: usize, expected: J, tol: f64) {
    let f = &r.frequency_response.frequency_response;
    check(
        [
            f.response.output[i],
            f.output_derivative[i],
            r.frequency_response.output_second_derivative[i],
            r.parameter_response.output_derivative[i],
            r.output_mixed_derivative[i],
        ],
        expected,
        tol,
    );
}
fn physical(r: &HarmonicFrequencyParameterResponse, energy: J, heat: J, supply: J, tol: f64) {
    let c = &r.frequency_response;
    let f = &c.frequency_response;
    let p = &r.parameter_response;
    for (actual, expected) in [
        (
            [
                f.response.mean_stored_energy_j,
                f.mean_stored_energy_derivative,
                c.mean_stored_energy_second_derivative,
                p.mean_stored_energy_derivative,
                r.mean_stored_energy_mixed_derivative,
            ],
            energy,
        ),
        (
            [
                f.response.dissipated_power_w,
                f.dissipated_power_derivative,
                c.dissipated_power_second_derivative,
                p.dissipated_power_derivative,
                r.dissipated_power_mixed_derivative,
            ],
            heat,
        ),
        (
            [
                f.response.supplied_power_w,
                f.supplied_power_derivative,
                c.supplied_power_second_derivative,
                p.supplied_power_derivative,
                r.supplied_power_mixed_derivative,
            ],
            supply,
        ),
        (
            [
                f.response.port_power_w[0],
                f.port_power_derivative[0],
                c.port_power_second_derivative[0],
                p.port_power_derivative[0],
                r.port_power_mixed_derivative[0],
            ],
            supply,
        ),
    ] {
        for (a, e) in
            actual
                .into_iter()
                .zip([expected.v, expected.w, expected.ww, expected.t, expected.wt])
        {
            real(a, e.re, tol);
        }
    }
    real(
        c.power_residual_second_derivative,
        0.0,
        tol * (1.0 + heat.ww.abs()),
    );
    real(
        r.power_residual_mixed_derivative,
        0.0,
        tol * (1.0 + heat.wt.abs()),
    );
    for residual in [
        f.response.relative_residual,
        f.relative_derivative_residual,
        c.relative_second_derivative_residual,
        p.relative_derivative_residual,
        r.relative_mixed_derivative_residual,
    ] {
        assert!(residual < 1e-12);
    }
}
fn source(omega: f64) -> J {
    let gain = J {
        v: C64::new(0.7, -0.3),
        t: C64::new(0.2, 0.15),
        ..J::default()
    };
    // Positive phase slope is a physical delay in exp(-i omega t).
    gain * (J::complex(C64::new(0.0, 1.0)) * J::frequency(omega) * J::real(0.17, -0.08)).exp()
}
fn solve(
    model: &HarmonicModel,
    omega: f64,
    direction: Direction<'_>,
    input: J,
) -> HarmonicFrequencyParameterResponse {
    model
        .factor(omega, || false)
        .unwrap()
        .response_with_frequency_parameter_derivative(
            &[input.v],
            direction,
            Some(&[input.w]),
            Some(&[input.ww]),
            Some(&[input.t]),
            Some(&[input.wt]),
            || false,
        )
        .unwrap()
}

#[test]
fn oscillator_curvature_and_mixed_jets_keep_mass_loss_ports_and_dispersive_sources() {
    for damping in [0.2, 1e-5] {
        let (mass, stiffness) = (2.0, 8.0);
        let original = fs_phs::mass_spring_damper(mass, stiffness, damping).unwrap();
        let model = HarmonicModel::from_ode(&original, budget(), || false).unwrap();
        for omega in [1.1, 2.0, 3.5] {
            for (dm, dk, dc, dg) in [
                (1.0, 0.0, 0.0, 0.0),
                (0.0, 1.0, 0.0, 0.0),
                (0.0, 0.0, 1.0, 0.0),
                (0.0, 0.0, 0.0, 1.0),
                (0.3, -0.2, 0.1, -0.4),
            ] {
                // At the exact velocity peak several power derivatives vanish
                // by cancellation of terms near 1e15 at this Q. Resolve the
                // complete high-Q work jet one tenth of a physical linewidth
                // away. The moderate-loss case still checks the exact peak;
                // existing resonance tests independently check high-Q peaks.
                let omega = if damping < 1e-4 && omega == 2.0 {
                    omega + 0.1 * damping / mass
                } else {
                    omega
                };
                let input = source(omega);
                let s = J::frequency(omega) * J::complex(C64::new(0.0, -1.0));
                let m = J::real(mass, dm);
                let k = J::real(stiffness, dk);
                let c = J::real(damping, dc);
                let g = J::real(1.0, dg);
                let q = (g * input) / (k + m * s * s + c * s);
                let v = s * q;
                let dq = [dk, 0.0, 0.0, -dm / (mass * mass)];
                let dr = [0.0, 0.0, 0.0, dc];
                let dp = [0.0, dg];
                let direction = Direction {
                    storage: Some(&dq),
                    dissipation: Some(&dr),
                    ports: Some(&dp),
                    ..Default::default()
                };
                let result = solve(&model, omega, direction, input);
                let tol = if damping < 1e-4 { 3e-7 } else { 3e-10 };
                state(&result, 0, q, tol);
                state(&result, 1, m * v, tol);
                effort(&result, 0, k * q, tol);
                effort(&result, 1, v, tol);
                output(&result, 0, g * v, tol);
                physical(
                    &result,
                    (k * q.norm() + m * v.norm()).scale(0.25),
                    (c * v.norm()).scale(0.5),
                    (input.conj() * g * v).scale(0.5),
                    tol,
                );
                let curvature = model
                    .factor(omega, || false)
                    .unwrap()
                    .response_with_frequency_curvature(
                        &[input.v],
                        Some(&[input.w]),
                        Some(&[input.ww]),
                        || false,
                    )
                    .unwrap();
                assert_eq!(curvature, result.frequency_response);
            }
        }
    }
}

#[test]
fn maxwell_mixed_jets_keep_relaxation_state_storage_and_independent_heat() {
    let (mass, equilibrium, damping, branch, tau) = (0.4, 110.0, 0.2, 300.0_f64, 0.025);
    let original = fs_phs::mass_spring_damper(mass, equilibrium, damping)
        .unwrap()
        .with_relaxation_branches(vec![RelaxationBranch {
            projection: vec![1.0, 0.0],
            stiffness: branch,
            relaxation_time_s: tau,
        }])
        .unwrap();
    let model = HarmonicModel::from_ode(&original, budget(), || false).unwrap();
    for omega in [0.1, 10.0, 30.0, 100.0] {
        let input = source(omega);
        let s = J::frequency(omega) * J::complex(C64::new(0.0, -1.0));
        for (dk, dtau) in [(1.0, 0.0), (0.0, 1.0)] {
            let b = J::real(branch, dk);
            let time = J::real(tau, dtau);
            let root = J::real(branch.sqrt(), 0.5 * dk / branch.sqrt());
            let memory = J::real(1.0, 0.0) + s * time;
            let q = input
                / (J::real(equilibrium, 0.0)
                    + s * s * J::real(mass, 0.0)
                    + s * J::real(damping, 0.0)
                    + b * s * time / memory);
            let v = s * q;
            let z = root * q / memory;
            let mismatch = root * q - z;
            let cross = -0.5 * dk / branch.sqrt();
            let dq = [dk, 0.0, cross, 0.0, 0.0, 0.0, cross, 0.0, 0.0];
            let mut dr = [0.0; 9];
            dr[8] = -dtau / (tau * tau);
            let result = solve(
                &model,
                omega,
                Direction {
                    storage: Some(&dq),
                    dissipation: Some(&dr),
                    ..Default::default()
                },
                input,
            );
            state(&result, 0, q, 5e-9);
            state(&result, 1, v.scale(mass), 5e-9);
            state(&result, 2, z, 5e-9);
            effort(&result, 0, q.scale(equilibrium) + root * mismatch, 5e-9);
            effort(&result, 1, v, 5e-9);
            effort(&result, 2, -mismatch, 5e-9);
            output(&result, 0, v, 5e-9);
            physical(
                &result,
                (q.norm().scale(equilibrium) + v.norm().scale(mass) + mismatch.norm()).scale(0.25),
                (v.norm().scale(damping) + mismatch.norm() / time).scale(0.5),
                (input.conj() * v).scale(0.5),
                5e-9,
            );
        }
    }
}

#[test]
fn descriptor_rlc_jets_keep_kirchhoff_reactions_and_all_source_orders() {
    let (r, l, c) = (0.02, 0.5, 0.125);
    let circuit = assemble_circuit(&CircuitGraph {
        node_count: 4,
        branches: vec![
            (1, 0, Branch::VoltageSource { port: 0 }),
            (1, 2, Branch::Inductor { henries: l }),
            (2, 3, Branch::Capacitor { farads: c }),
            (3, 0, Branch::Resistor { ohms: r }),
        ],
        transformers: vec![],
    })
    .unwrap();
    let model = HarmonicModel::from_descriptor(&circuit.system, budget(), || false).unwrap();
    let n = circuit.system.state_dim();
    let flux = circuit.flux_index[0];
    let charge = circuit.charge_index[0];
    let resistor = circuit.node_potential_index[2];
    for omega in [1.0, 4.0, 7.0] {
        let input = source(omega);
        let s = J::frequency(omega) * J::complex(C64::new(0.0, -1.0));
        for (dr, dl, dc) in [(1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)] {
            let resistance = J::real(r, dr);
            let inductance = J::real(l, dl);
            let capacitance = J::real(c, dc);
            let current = input / (resistance + s * inductance + (s * capacitance).reciprocal());
            let stored_charge = current / s;
            let mut dq = vec![0.0; n * n];
            dq[flux * n + flux] = -dl / (l * l);
            dq[charge * n + charge] = -dc / (c * c);
            let mut loss = vec![0.0; n * n];
            loss[resistor * n + resistor] = -dr / (r * r);
            let result = solve(
                &model,
                omega,
                Direction {
                    storage: Some(&dq),
                    dissipation: Some(&loss),
                    ..Default::default()
                },
                input,
            );
            output(&result, 0, current, 3e-8);
            state(&result, flux, inductance * current, 3e-8);
            state(&result, charge, stored_charge, 3e-8);
            effort(&result, flux, current, 3e-8);
            state(&result, circuit.node_potential_index[0], input, 3e-8);
            state(&result, resistor, resistance * current, 3e-8);
            state(
                &result,
                circuit.node_potential_index[1],
                resistance * current + stored_charge / capacitance,
                3e-8,
            );
            physical(
                &result,
                (inductance * current.norm() + stored_charge.norm() / capacitance).scale(0.25),
                (resistance * current.norm()).scale(0.5),
                (input.conj() * current).scale(0.5),
                3e-8,
            );
        }
    }
}

#[test]
fn exact_coherent_notch_has_finite_full_jets_and_nonzero_internal_energy() {
    let mut q = vec![0.0; 16];
    for (i, k) in [2.0, 1.0, 6.0, 1.0].into_iter().enumerate() {
        q[5 * i] = k;
    }
    let mut j = vec![0.0; 16];
    j[1] = 1.0;
    j[4] = -1.0;
    j[11] = 1.0;
    j[14] = -1.0;
    let original = PortHamiltonian::new(
        4,
        1,
        j,
        vec![0.0; 16],
        vec![0.0, 1.0, 0.0, 1.0],
        Box::new(QuadraticStorage::new(q, 4).unwrap()),
    )
    .unwrap();
    let model = HarmonicModel::from_ode(&original, budget(), || false).unwrap();
    let mut dq = [0.0; 16];
    dq[0] = 1.0;
    let input = J::complex(C64::ONE);
    let result = solve(
        &model,
        2.0,
        Direction {
            storage: Some(&dq),
            ..Default::default()
        },
        input,
    );
    let s = J::frequency(2.0) * J::complex(C64::new(0.0, -1.0));
    let k0 = J::real(2.0, 1.0);
    let k1 = J::real(6.0, 0.0);
    let q0 = input / (s * s + k0);
    let q1 = input / (s * s + k1);
    let v0 = s * q0;
    let v1 = s * q1;
    let total = v0 + v1;
    let observed = result.frequency_response.frequency_response.response.output[0];
    // An analytic exact zero can retain a finite solve-rounding residual.
    // No amplitude floor replaces it or assigns a phase to a machine zero.
    near(observed, C64::ZERO, 5e-14);
    if observed == C64::ZERO {
        assert_eq!(
            result
                .parameter_response
                .output_logarithmic_derivative(0)
                .unwrap(),
            None
        );
    }
    output(&result, 0, total, 5e-13);
    state(&result, 0, q0, 5e-13);
    state(&result, 2, q1, 5e-13);
    physical(
        &result,
        (k0 * q0.norm() + v0.norm() + k1 * q1.norm() + v1.norm()).scale(0.25),
        J::default(),
        (input.conj() * total).scale(0.5),
        5e-13,
    );
    assert!(
        result
            .frequency_response
            .frequency_response
            .response
            .mean_stored_energy_j
            > 0.0
    );
}

#[test]
fn curvature_and_mixed_work_retain_weak_remainders_between_strong_states() {
    let identity = vec![1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
    let mut g = vec![0.0; 18];
    for i in 0..3 {
        g[i * 6 + 2 * i] = 1.0;
        g[i * 6 + 2 * i + 1] = 1.0;
    }
    let original = PortHamiltonian::new(
        3,
        6,
        vec![0.0; 9],
        identity.clone(),
        g,
        Box::new(QuadraticStorage::new(identity, 3).unwrap()),
    )
    .unwrap();
    let model = HarmonicModel::from_ode(&original, budget(), || false).unwrap();
    let large = 2.0f64.powi(54);
    let x = [large, 1.0, large];
    let mixed = [1.0, 1.0, -1.0];
    let mut u = [C64::ZERO; 6];
    let mut uw = [C64::ZERO; 6];
    let mut uww = [C64::ZERO; 6];
    let mut uwt = [C64::ZERO; 6];
    for i in 0..3 {
        u[2 * i] = C64::new(x[i], -x[i]);
        uw[2 * i] = C64::new(0.0, -x[i]); // Cancels i*E*x, so x_omega=0.
        uww[2 * i + 1] = C64::new(mixed[i], -mixed[i]);
        uwt[2 * i + 1] = uww[2 * i + 1];
    }
    let result = model
        .factor(1.0, || false)
        .unwrap()
        .response_with_frequency_parameter_derivative(
            &u,
            Direction::default(),
            Some(&uw),
            Some(&uww),
            None,
            Some(&uwt),
            || false,
        )
        .unwrap();
    assert!(
        result
            .frequency_response
            .frequency_response
            .state_derivative
            .iter()
            .all(|&x| x == C64::ZERO)
    );
    for i in 0..3 {
        near(
            result.frequency_response.state_second_derivative[i],
            C64::from_re(mixed[i]),
            3e-14,
        );
        near(
            result.state_mixed_derivative[i],
            C64::from_re(mixed[i]),
            3e-14,
        );
    }
    real(
        result
            .frequency_response
            .mean_stored_energy_second_derivative,
        0.5,
        3e-14,
    );
    real(result.mean_stored_energy_mixed_derivative, 0.5, 3e-14);
    real(
        result.frequency_response.dissipated_power_second_derivative,
        1.0,
        3e-14,
    );
    real(
        result.frequency_response.supplied_power_second_derivative,
        1.0,
        3e-14,
    );
    real(result.dissipated_power_mixed_derivative, 1.0, 3e-14);
    real(result.supplied_power_mixed_derivative, 1.0, 3e-14);
    real(
        result.frequency_response.power_residual_second_derivative,
        0.0,
        3e-14,
    );
    real(result.power_residual_mixed_derivative, 0.0, 3e-14);
}

#[test]
fn complete_jet_admission_cancellation_and_replay_preserve_the_factor() {
    let original = fs_phs::mass_spring_damper(1.0, 4.0, 0.2).unwrap();
    let model = HarmonicModel::from_ode(&original, budget(), || false).unwrap();
    let factor = model.factor(2.0, || false).unwrap();
    let input = [C64::ONE];
    let input_w = [C64::new(0.2, -0.1)];
    let before = factor.response(&input, || false).unwrap();
    let mut polls = 0;
    let expected = factor
        .response_with_frequency_parameter_derivative(
            &input,
            Direction::default(),
            Some(&input_w),
            None,
            None,
            None,
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    for stop in [1, polls / 2, polls] {
        let mut count = 0;
        assert!(matches!(
            factor.response_with_frequency_parameter_derivative(
                &input,
                Direction::default(),
                Some(&input_w),
                None,
                None,
                None,
                || {
                    count += 1;
                    count >= stop
                },
            ),
            Err(HarmonicError::Cancelled)
        ));
    }
    for bad in [vec![], vec![C64::new(f64::NAN, 0.0)]] {
        for slot in 0..4 {
            let mut law = [None; 4];
            law[slot] = Some(bad.as_slice());
            assert!(
                factor
                    .response_with_frequency_parameter_derivative(
                        &input,
                        Direction::default(),
                        law[0],
                        law[1],
                        law[2],
                        law[3],
                        || false,
                    )
                    .is_err()
            );
        }
    }
    assert_eq!(factor.response(&input, || false).unwrap(), before);
    assert_eq!(
        factor
            .response_with_frequency_parameter_derivative(
                &input,
                Direction::default(),
                Some(&input_w),
                None,
                None,
                None,
                || false,
            )
            .unwrap(),
        expected
    );
    // n=2,m=1 has a 100-work ordinary score. Curvature admits exactly three;
    // the full mixed response admits six before computing any state derivative.
    let limited = HarmonicModel::from_ode(
        &original,
        HarmonicBudget {
            max_dense_work: 300,
            ..budget()
        },
        || false,
    )
    .unwrap();
    let limited = limited.factor(2.0, || false).unwrap();
    limited
        .response_with_frequency_curvature(&input, None, None, || false)
        .unwrap();
    assert!(matches!(
        limited.response_with_frequency_parameter_derivative(
            &input,
            Direction::default(),
            None,
            None,
            None,
            None,
            || false,
        ),
        Err(HarmonicError::Invalid { .. })
    ));
}
