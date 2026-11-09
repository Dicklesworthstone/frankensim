//! G1/G3 original mechanical directions against independent physical equations.
use fs_math::{c64::C64, det};
use fs_phs::RelaxationBranch;
use fs_phs::harmonic::{
    HarmonicBudget, HarmonicError, HarmonicFrequencyParameterResponse, HarmonicModel,
};
use fs_phs::mechanical::{
    LinearMechanicalSystem, MechanicalHarmonicDirection, MechanicalHarmonicDirectionControls,
    MechanicalRelaxationDirection,
};

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
fn controls() -> MechanicalHarmonicDirectionControls {
    MechanicalHarmonicDirectionControls {
        max_state_dimension: 64,
        max_dense_bytes: 2_000_000,
        max_multiply_terms: 4_000_000,
    }
}
fn physical(r: &HarmonicFrequencyParameterResponse, energy: J, heat: J, supply: J) {
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
    ] {
        for (a, e) in
            actual
                .into_iter()
                .zip([expected.v, expected.w, expected.ww, expected.t, expected.wt])
        {
            real(a, e.re, 2e-10);
        }
    }
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

// A scalar/2-by-2 physical dynamic stiffness is an independent resolvent:
// D = M*s^2+C*s+K+sum(k*b*b^T*s*tau/(1+s*tau)), with s=-i*omega.
// It does not use the production momentum/memory matrices or their derivatives.
#[allow(clippy::too_many_arguments)]
fn check_physical(
    mass: &[f64],
    damping: &[f64],
    stiffness: &[f64],
    force: &[f64],
    n: usize,
    m: usize,
    branches: Vec<RelaxationBranch>,
    direction: MechanicalHarmonicDirection<'_>,
) {
    let source =
        LinearMechanicalSystem::new(mass, damping, stiffness, force, n, m, branches).unwrap();
    let prepared = source
        .prepare_harmonic_parameter_direction(direction, controls(), || false)
        .unwrap();
    assert_eq!(
        (
            prepared.coordinates(),
            prepared.state_dim(),
            prepared.port_dim()
        ),
        (n, 2 * n + source.branches().len(), m)
    );
    let model = HarmonicModel::from_ode(source.system(), budget(), || false).unwrap();
    let variable = |values: &[f64], change: Option<&[f64]>| -> Vec<J> {
        values
            .iter()
            .enumerate()
            .map(|(i, &v)| J::real(v, change.map_or(0.0, |a| a[i])))
            .collect()
    };
    let mass = variable(mass, direction.mass);
    let damping = variable(damping, direction.damping);
    let stiffness = variable(stiffness, direction.relaxed_stiffness);
    let force = variable(force, direction.force_map);
    for omega in [0.2, 1.6, 2.1, 5.0] {
        let s = J::complex(C64::new(0.0, -1.0)) * J::frequency(omega);
        let input: Vec<_> = (0..m)
            .map(|i| {
                let i = i as f64;
                let gain = J {
                    v: C64::new(0.7 + 0.2 * i, -0.3 + 0.1 * i),
                    t: C64::new(0.2, -0.15 + 0.2 * i),
                    ..J::default()
                };
                gain * (J::complex(C64::new(0.0, 1.0))
                    * J::frequency(omega)
                    * J::real(0.17 + 0.05 * i, -0.08))
                .exp()
            })
            .collect();
        let mut dynamic = vec![J::default(); n * n];
        for i in 0..n * n {
            dynamic[i] = mass[i] * s * s + damping[i] * s + stiffness[i];
        }
        let mut arms = Vec::new();
        for (index, arm) in source.branches().iter().enumerate() {
            let d = direction
                .branches
                .map_or(MechanicalRelaxationDirection::default(), |a| a[index]);
            let projection = variable(&arm.projection, d.projection);
            let k = J::real(arm.stiffness, d.stiffness);
            let root = J::real(
                det::sqrt(arm.stiffness),
                d.stiffness / (2.0 * det::sqrt(arm.stiffness)),
            );
            let tau = J::real(arm.relaxation_time_s, d.relaxation_time_s);
            let memory = J::real(1.0, 0.0) + s * tau;
            for i in 0..n {
                for j in 0..n {
                    dynamic[i * n + j] =
                        dynamic[i * n + j] + k * projection[i] * projection[j] * s * tau / memory;
                }
            }
            arms.push((projection, root, tau, memory));
        }
        let rhs: Vec<_> = (0..n)
            .map(|i| (0..m).fold(J::default(), |sum, p| sum + force[i * m + p] * input[p]))
            .collect();
        let q = if n == 1 {
            vec![rhs[0] / dynamic[0]]
        } else {
            assert_eq!(n, 2);
            let determinant = dynamic[0] * dynamic[3] - dynamic[1] * dynamic[2];
            vec![
                (rhs[0] * dynamic[3] - dynamic[1] * rhs[1]) / determinant,
                (dynamic[0] * rhs[1] - rhs[0] * dynamic[2]) / determinant,
            ]
        };
        let v: Vec<_> = q.iter().map(|&q| s * q).collect();
        let momentum: Vec<_> = (0..n)
            .map(|i| (0..n).fold(J::default(), |sum, j| sum + mass[i * n + j] * v[j]))
            .collect();
        let mut strain_effort: Vec<_> = (0..n)
            .map(|i| (0..n).fold(J::default(), |sum, j| sum + stiffness[i * n + j] * q[j]))
            .collect();
        let mut energy = J::default();
        let mut heat = J::default();
        for i in 0..n {
            for j in 0..n {
                energy = energy
                    + (q[i].conj() * stiffness[i * n + j] * q[j]
                        + v[i].conj() * mass[i * n + j] * v[j])
                        .scale(0.25);
                heat = heat + (v[i].conj() * damping[i * n + j] * v[j]).scale(0.5);
            }
        }
        let mut histories = Vec::new();
        for (projection, root, tau, memory) in arms {
            let strain = (0..n).fold(J::default(), |sum, i| sum + projection[i] * q[i]);
            let z = root * strain / memory;
            let extension = root * strain - z;
            for i in 0..n {
                strain_effort[i] = strain_effort[i] + root * projection[i] * extension;
            }
            energy = energy + extension.norm().scale(0.25);
            heat = heat + (extension.norm() / tau).scale(0.5);
            histories.push((z, -extension));
        }
        let native = model
            .factor(omega, || false)
            .unwrap()
            .response_with_frequency_parameter_derivative(
                &input.iter().map(|j| j.v).collect::<Vec<_>>(),
                prepared.as_harmonic_direction(),
                Some(&input.iter().map(|j| j.w).collect::<Vec<_>>()),
                Some(&input.iter().map(|j| j.ww).collect::<Vec<_>>()),
                Some(&input.iter().map(|j| j.t).collect::<Vec<_>>()),
                Some(&input.iter().map(|j| j.wt).collect::<Vec<_>>()),
                || false,
            )
            .unwrap();
        for i in 0..n {
            state(&native, i, q[i], 2e-10);
            state(&native, n + i, momentum[i], 2e-10);
            effort(&native, i, strain_effort[i], 2e-10);
            effort(&native, n + i, v[i], 2e-10);
        }
        for (arm, (z, e)) in histories.into_iter().enumerate() {
            state(&native, 2 * n + arm, z, 2e-10);
            effort(&native, 2 * n + arm, e, 2e-10);
            if source.branches()[arm].projection.iter().all(|b| *b == 0.0)
                && direction.branches.unwrap()[arm].projection.is_some()
            {
                assert_eq!(z.v, C64::ZERO);
                assert!(z.t.abs() > 0.0, "newly coupled memory must survive");
            }
        }
        let mut supply = J::default();
        for p in 0..m {
            let y = (0..n).fold(J::default(), |sum, i| sum + force[i * m + p] * v[i]);
            output(&native, p, y, 2e-10);
            let work = (input[p].conj() * y).scale(0.5);
            let c = &native.frequency_response;
            let f = &c.frequency_response;
            for (actual, expected) in [
                f.response.port_power_w[p],
                f.port_power_derivative[p],
                c.port_power_second_derivative[p],
                native.parameter_response.port_power_derivative[p],
                native.port_power_mixed_derivative[p],
            ]
            .into_iter()
            .zip([work.v, work.w, work.ww, work.t, work.wt])
            {
                real(actual, expected.re, 2e-10);
            }
            supply = supply + work;
        }
        physical(&native, energy, heat, supply);
    }
}

#[test]
fn original_oscillator_and_maxwell_directions_reach_every_native_field_and_work_order() {
    let branch_direction = [MechanicalRelaxationDirection {
        projection: Some(&[0.2]),
        stiffness: -0.4,
        relaxation_time_s: 0.07,
    }];
    check_physical(
        &[2.0],
        &[0.13],
        &[8.0],
        &[1.3],
        1,
        1,
        vec![RelaxationBranch {
            projection: vec![0.7],
            stiffness: 2.25,
            relaxation_time_s: 0.4,
        }],
        MechanicalHarmonicDirection {
            mass: Some(&[-0.3]),
            damping: Some(&[0.04]),
            relaxed_stiffness: Some(&[0.8]),
            force_map: Some(&[-0.2]),
            branches: Some(&branch_direction),
        },
    );
    check_physical(
        &[2.0],
        &[1e-5],
        &[8.0],
        &[1.3],
        1,
        1,
        vec![],
        MechanicalHarmonicDirection {
            mass: Some(&[0.1]),
            damping: Some(&[-0.2]),
            relaxed_stiffness: Some(&[0.3]),
            force_map: Some(&[0.4]),
            branches: None,
        },
    );
}

#[test]
fn coupled_inertia_gyroscopic_loss_ports_and_all_equal_time_arms_remain_physical() {
    let d = [
        MechanicalRelaxationDirection {
            projection: Some(&[0.1, 0.02]),
            stiffness: 0.4,
            relaxation_time_s: 0.07,
        },
        MechanicalRelaxationDirection {
            projection: Some(&[0.3, -0.2]),
            stiffness: -0.05,
            relaxation_time_s: -0.1,
        },
        MechanicalRelaxationDirection {
            projection: Some(&[-0.05, 0.1]),
            stiffness: -0.2,
            relaxation_time_s: 0.15,
        },
    ];
    check_physical(
        &[2.0, 0.3, 0.3, 1.4],
        &[0.25, 0.19, -0.07, 0.18],
        &[6.0, -1.0, -1.0, 4.0],
        &[1.0, 0.25, -0.4, 1.2],
        2,
        2,
        vec![
            RelaxationBranch {
                projection: vec![0.8, -0.3],
                stiffness: 2.25,
                relaxation_time_s: 0.4,
            },
            RelaxationBranch {
                projection: vec![0.0, 0.0],
                stiffness: 0.5,
                relaxation_time_s: 0.4,
            },
            RelaxationBranch {
                projection: vec![0.2, 1.0],
                stiffness: 0.8,
                relaxation_time_s: 1.2,
            },
        ],
        MechanicalHarmonicDirection {
            mass: Some(&[0.2, -0.05, -0.05, -0.3]),
            damping: Some(&[-0.04, 0.12, -0.08, 0.03]),
            relaxed_stiffness: Some(&[0.5, -0.13, -0.13, -0.2]),
            force_map: Some(&[0.1, -0.04, 0.08, 0.06]),
            branches: Some(&d),
        },
    );
}

#[test]
fn branch_reparameterization_keeps_dynamics_and_returns_the_physical_stress_chain_rule() {
    let source = LinearMechanicalSystem::new(
        &[2.0],
        &[0.2],
        &[3.0],
        &[1.0],
        1,
        1,
        vec![RelaxationBranch {
            projection: vec![0.75],
            stiffness: 4.0,
            relaxation_time_s: 0.5,
        }],
    )
    .unwrap();
    // b(theta)=exp(theta/2)b, k(theta)=exp(-theta)k keeps sqrt(k)b fixed.
    let branch = [MechanicalRelaxationDirection {
        projection: Some(&[0.375]),
        stiffness: -4.0,
        relaxation_time_s: 0.0,
    }];
    let prepared = source
        .prepare_harmonic_parameter_direction(
            MechanicalHarmonicDirection {
                branches: Some(&branch),
                ..Default::default()
            },
            controls(),
            || false,
        )
        .unwrap();
    let d = prepared.as_harmonic_direction();
    assert!(
        d.storage
            .unwrap()
            .iter()
            .chain(d.dissipation.unwrap())
            .chain(d.interconnection.unwrap())
            .chain(d.ports.unwrap())
            .all(|v| *v == 0.0)
    );
    assert_eq!(prepared.branch_root_derivatives(), &[-1.0]);
    let model = HarmonicModel::from_ode(source.system(), budget(), || false).unwrap();
    let r = model
        .factor(1.7, || false)
        .unwrap()
        .response_with_parameter_derivative(&[C64::ONE], d, None, || false)
        .unwrap();
    assert!(r.state_derivative.iter().all(|v| *v == C64::ZERO));
    let e = r.response.effort[2];
    let stress = e.scale(-2.0);
    let stress_derivative = e.scale(-prepared.branch_root_derivatives()[0]);
    assert!(stress_derivative.abs() > 0.0);
    near(
        stress.scale(0.375) + stress_derivative.scale(0.75),
        C64::ZERO,
        1e-15,
    );
}

#[test]
fn physical_mass_and_time_slopes_do_not_square_reciprocals_before_restoring_units() {
    for mass in [1e-200, 1e200] {
        for tau in [1e-200, 1e200] {
            let stiffness = if tau < 1.0 { 1e200 } else { 1e-200 };
            let root = det::sqrt(stiffness);
            let source = LinearMechanicalSystem::new(
                &[mass],
                &[0.0],
                &[1.0],
                &[1.0],
                1,
                1,
                vec![RelaxationBranch {
                    projection: vec![1.0 / root],
                    stiffness,
                    relaxation_time_s: tau,
                }],
            )
            .unwrap();
            let branches = [MechanicalRelaxationDirection {
                projection: None,
                stiffness,
                relaxation_time_s: tau,
            }];
            let prepared = source
                .prepare_harmonic_parameter_direction(
                    MechanicalHarmonicDirection {
                        mass: Some(&[mass]),
                        branches: Some(&branches),
                        ..Default::default()
                    },
                    controls(),
                    || false,
                )
                .unwrap();
            let d = prepared.as_harmonic_direction();
            let expected = -1.0 / mass;
            assert!((d.storage.unwrap()[4] / expected - 1.0).abs() < 2e-14);
            assert!((prepared.branch_root_derivatives()[0] / (0.5 * root) - 1.0).abs() < 2e-14);
            assert!((prepared.branch_rate_derivatives()[0] / (-1.0 / tau) - 1.0).abs() < 2e-14);
            real(d.storage.unwrap()[0], 1.0, 2e-14);
            real(d.storage.unwrap()[2], -0.5, 2e-14);
        }
    }
}

#[test]
fn inverse_mass_chain_rule_retains_the_first_contraction_remainder() {
    // M has the exact binary Cholesky factor [[1,0],[1,1]], hence the
    // cached inverse [[2,-1],[-1,1]] is exact. The p0 diagonal direction
    // is -(4*L-4*(1/4)-4*L)=1. Rounding W*dM before W erases it.
    let source = LinearMechanicalSystem::new(
        &[1.0, 1.0, 1.0, 2.0],
        &[0.0; 4],
        &[1.0, 0.0, 0.0, 1.0],
        &[],
        2,
        0,
        vec![],
    )
    .unwrap();
    let large = 2.0_f64.powi(54);
    let prepared = source
        .prepare_harmonic_parameter_direction(
            MechanicalHarmonicDirection {
                mass: Some(&[large, 0.25, 0.25, -4.0 * large]),
                ..Default::default()
            },
            controls(),
            || false,
        )
        .unwrap();
    assert_eq!(
        prepared.as_harmonic_direction().storage.unwrap()[2 * 4 + 2],
        1.0
    );
}

#[test]
fn bounded_complete_preparation_refuses_bad_directions_and_cancellation_without_mutating_source() {
    let source = LinearMechanicalSystem::new(
        &[2.0],
        &[0.2],
        &[3.0],
        &[1.0],
        1,
        1,
        vec![RelaxationBranch {
            projection: vec![1.0],
            stiffness: 4.0,
            relaxation_time_s: 0.5,
        }],
    )
    .unwrap();
    let exact = MechanicalHarmonicDirectionControls {
        max_state_dimension: 3,
        max_dense_bytes: 328,
        max_multiply_terms: 248,
    };
    let zero = MechanicalHarmonicDirection::default();
    assert!(
        source
            .prepare_harmonic_parameter_direction(zero, exact, || false)
            .is_ok()
    );
    for budget in [
        MechanicalHarmonicDirectionControls {
            max_state_dimension: 2,
            ..exact
        },
        MechanicalHarmonicDirectionControls {
            max_dense_bytes: 327,
            ..exact
        },
        MechanicalHarmonicDirectionControls {
            max_multiply_terms: 247,
            ..exact
        },
    ] {
        assert!(matches!(
            source.prepare_harmonic_parameter_direction(zero, budget, || false),
            Err(HarmonicError::Invalid { .. })
        ));
    }
    for d in [
        MechanicalHarmonicDirection {
            mass: Some(&[]),
            ..zero
        },
        MechanicalHarmonicDirection {
            damping: Some(&[f64::NAN]),
            ..zero
        },
        MechanicalHarmonicDirection {
            branches: Some(&[]),
            ..zero
        },
        MechanicalHarmonicDirection {
            force_map: Some(&[]),
            ..zero
        },
    ] {
        assert!(matches!(
            source.prepare_harmonic_parameter_direction(d, controls(), || false),
            Err(HarmonicError::Invalid { .. })
        ));
    }
    let before = source.state(&[0.7], &[-0.3], &[0.4]).unwrap();
    for stop in [1, 3, 7] {
        let mut polls = 0;
        assert!(matches!(
            source.prepare_harmonic_parameter_direction(zero, controls(), || {
                polls += 1;
                polls == stop
            }),
            Err(HarmonicError::Cancelled)
        ));
    }
    assert_eq!(source.state(&[0.7], &[-0.3], &[0.4]).unwrap(), before);
    assert!(
        source
            .prepare_harmonic_parameter_direction(zero, exact, || false)
            .is_ok()
    );
    let two = LinearMechanicalSystem::new(
        &[2.0, 0.0, 0.0, 1.0],
        &[0.0; 4],
        &[1.0, 0.0, 0.0, 1.0],
        &[],
        2,
        0,
        vec![],
    )
    .unwrap();
    assert!(matches!(
        two.prepare_harmonic_parameter_direction(
            MechanicalHarmonicDirection {
                mass: Some(&[0.0, 1.0, 0.0, 0.0]),
                ..zero
            },
            controls(),
            || false
        ),
        Err(HarmonicError::Invalid { .. })
    ));
}
