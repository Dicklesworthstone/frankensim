use super::*;
use fs_modal::harmonic::resonance::HarmonicPowerCurvatureSample;

// Independent analytic differentiation of the small original physical laws.
#[derive(Clone, Copy)]
struct Jet([C64; 3]);
impl Jet {
    fn constant(x: f64) -> Self {
        Self([C64::from_re(x), C64::ZERO, C64::ZERO])
    }
    fn scale(self, x: f64) -> Self {
        Self(self.0.map(|z| z.scale(x)))
    }
    fn inverse(self) -> Self {
        let a = C64::ONE / self.0[0];
        Self([
            a,
            -(self.0[1] * a * a),
            self.0[1] * self.0[1] * a * a * a.scale(2.0) - self.0[2] * a * a,
        ])
    }
}
impl core::ops::Add for Jet {
    type Output = Self;
    fn add(self, b: Self) -> Self {
        Self(core::array::from_fn(|i| self.0[i] + b.0[i]))
    }
}
impl core::ops::Sub for Jet {
    type Output = Self;
    fn sub(self, b: Self) -> Self {
        self + b.scale(-1.0)
    }
}
impl core::ops::Mul for Jet {
    type Output = Self;
    fn mul(self, b: Self) -> Self {
        Self([
            self.0[0] * b.0[0],
            self.0[1] * b.0[0] + self.0[0] * b.0[1],
            self.0[2] * b.0[0] + (self.0[1] * b.0[1]).scale(2.0) + self.0[0] * b.0[2],
        ])
    }
}
impl core::ops::Div for Jet {
    type Output = Self;
    fn div(self, b: Self) -> Self {
        self * b.inverse()
    }
}
fn work(flow: Jet, effort: Jet) -> [f64; 3] {
    let f = Jet(flow.0.map(C64::conj));
    (f * effort).0.map(|x| 0.5 * x.re)
}
fn power(actual: HarmonicPowerCurvatureSample, expected: [f64; 3]) {
    for (a, b) in [
        actual.sample.average_power_w,
        actual.sample.angular_derivative,
        actual.angular_second_derivative,
    ]
    .into_iter()
    .zip(expected)
    {
        assert!((a - b).abs() < 3e-9 * b.abs().max(1.0), "{a:e} != {b:e}");
    }
}
fn zero_balances(result: &StructuralInteriorRelaxationCurvatureResponse) {
    for b in [
        result.power.structural_balance(),
        result.power.acoustic_balance(),
        result.power.relaxation_balance(),
        result.power.total_balance(),
    ] {
        power(b.unwrap(), [0.0; 3]);
    }
    assert!(result.relative_second_derivative_residual < 1e-9);
}
fn check(actual: [C64; 3], expected: Jet) {
    for (a, b) in actual.into_iter().zip(expected.0) {
        close(a, b, 3e-9);
    }
}
fn solve(
    prepared: &Prepared<'_>,
    omega: f64,
    f: Jet,
    q: Jet,
) -> StructuralInteriorRelaxationCurvatureResponse {
    prepared
        .factor(omega, &mut || false)
        .unwrap()
        .solve_with_frequency_curvature(
            &[f.0[0]],
            &[q.0[0].scale(0.25); 4],
            Some((&[f.0[1]], &[q.0[1].scale(0.25); 4])),
            Some((&[f.0[2]], &[q.0[2].scale(0.25); 4])),
            &mut || false,
        )
        .unwrap()
}
fn check_uniform(
    result: &StructuralInteriorRelaxationCurvatureResponse,
    u: Jet,
    v: Jet,
    p: Jet,
    psi: Jet,
    stress: Jet,
) {
    let r = &result.frequency_response;
    let base = &r.response;
    check(
        [
            base.structural_displacement[0],
            r.structural_displacement_derivative_per_rad_s[0],
            result.structural_displacement_second_derivative[0],
        ],
        u,
    );
    check(
        [
            base.structural_velocity[0],
            r.structural_velocity_derivative_per_rad_s[0],
            result.structural_velocity_second_derivative[0],
        ],
        v,
    );
    check(
        [
            base.fluid_force[0],
            r.fluid_force_derivative_per_rad_s[0],
            result.fluid_force_second_derivative[0],
        ],
        p.scale(0.5),
    );
    for i in 0..4 {
        check(
            [
                base.pressure_pa[i],
                r.pressure_derivative_pa_per_rad_s[i],
                result.pressure_second_derivative[i],
            ],
            p,
        );
        check(
            [
                base.pressure_potential_pa_s[i],
                r.pressure_potential_derivative_pa_s_per_rad_s[i],
                result.pressure_potential_second_derivative[i],
            ],
            psi,
        );
        check(
            [
                base.displaced_flow_m3_s[i],
                r.displaced_flow_derivative_m3_s_per_rad_s[i],
                result.displaced_flow_second_derivative[i],
            ],
            v.scale(0.125),
        );
    }
    check(
        [
            base.memory[0].stress,
            r.relaxation_stress_derivative_per_rad_s[0],
            result.relaxation_stress_second_derivative[0],
        ],
        stress,
    );
}

#[test]
fn full_memory_curvature_matches_original_cubic_with_source_phase_and_power() {
    let cavity = Cavity::new();
    let branches = [branch()];
    let prepared =
        Prepared::new(cavity.source(), &branches, Limits::default(), &mut || false).unwrap();
    let f = Jet([
        C64::new(0.7, -0.2),
        C64::new(0.03, 0.08),
        C64::new(-0.05, 0.02),
    ]);
    let q = Jet([
        C64::new(-0.02, 0.05),
        C64::new(0.01, -0.03),
        C64::new(0.02, 0.01),
    ]);
    for omega in [0.12, 0.5, 1.0, 1.4, 2.0, 8.0] {
        let s = Jet([C64::new(0.0, -omega), C64::new(0.0, -1.0), C64::ZERO]);
        let one = Jet::constant(1.0);
        let u = (f + (q / s).scale(0.6)) / (s * s + one + (s / (one + s)).scale(1.25));
        let v = s * u;
        let p = u.scale(-0.6) + (q / s).scale(1.2);
        let stress = (s * u / (one + s)).scale(1.25);
        let result = solve(&prepared, omega, f, q);
        check_uniform(&result, u, v, p, p / s, stress);
        power(result.power.structural_input, work(v, f));
        power(result.power.acoustic_input, work(p, q));
        power(result.power.fluid_to_structure, work(v, p.scale(0.5)));
        let heat = work(stress, stress).map(|x| x / 1.25);
        power(result.power.relaxation_dissipation, heat);
        power(result.memory[0].dissipation, heat);
        power(result.memory[0].structural_work, work(v, stress));
        power(result.memory[0].acoustic_work, [0.0; 3]);
        let stored = heat.map(|x| 0.5 * x);
        for (a, b) in [
            result.memory[0].mean_stored_energy_j,
            result.memory[0].mean_stored_energy_derivative,
            result.memory[0].mean_stored_energy_second_derivative,
        ]
        .into_iter()
        .zip(stored)
        {
            near(a, b);
        }
        let kinetic = work(v, v).map(|x| 0.5 * x);
        let strain = work(u, u).map(|x| 0.35 * x);
        let compression = work(p, p).map(|x| 5.0 / 12.0 * x);
        near(
            result.mean_stored_energy_derivative,
            kinetic[1] + strain[1] + compression[1] + stored[1],
        );
        near(
            result.mean_stored_energy_second_derivative,
            kinetic[2] + strain[2] + compression[2] + stored[2],
        );
        zero_balances(&result);
        // A common causal source delay rotates all fields but leaves every
        // physical power frequency jet unchanged, including material heat.
        let delay = 1.7;
        let phase = C64::new((omega * delay).cos(), (omega * delay).sin());
        let phase = Jet([
            phase,
            C64::new(0.0, delay) * phase,
            phase.scale(-delay * delay),
        ]);
        let delayed = solve(&prepared, omega, f * phase, q * phase);
        check_uniform(
            &delayed,
            u * phase,
            v * phase,
            p * phase,
            (p / s) * phase,
            stress * phase,
        );
        power(delayed.power.relaxation_dissipation, heat);
        power(delayed.power.structural_input, work(v, f));
        power(delayed.power.acoustic_input, work(p, q));
        zero_balances(&delayed);
    }
}

#[test]
fn wall_relaxation_curvature_keeps_acoustic_heat_and_instantaneous_loss_separate() {
    let mut cavity = Cavity::new();
    cavity.damping[0] = 0.17;
    for i in 1..5 {
        cavity.damping[i * 5 + i] = 0.03;
    }
    let branches = [MaxwellDampingBranch {
        strain: vec![0.0, 0.25, 0.25, 0.25, 0.25],
        stiffness: 0.8,
        relaxation_time_s: 0.4,
    }];
    let prepared =
        Prepared::new(cavity.source(), &branches, Limits::default(), &mut || false).unwrap();
    let f = Jet([
        C64::new(0.6, 0.1),
        C64::new(-0.03, 0.02),
        C64::new(0.01, 0.04),
    ]);
    let q = Jet([
        C64::new(0.07, -0.03),
        C64::new(0.02, 0.01),
        C64::new(-0.01, 0.02),
    ]);
    for omega in [0.2, 0.9, 1.5, 7.0] {
        let s = Jet([C64::new(0.0, -omega), C64::new(0.0, -1.0), C64::ZERO]);
        let admittance = s.scale(5.0 / 6.0)
            + Jet::constant(0.12)
            + Jet::constant(0.8) / (s + Jet::constant(2.5));
        let u = (f * admittance + q.scale(0.5))
            / ((s * s + s.scale(0.17) + Jet::constant(0.7)) * admittance + s.scale(0.25));
        let v = s * u;
        let p = (q - v.scale(0.5)) / admittance;
        let stress = p.scale(0.8) / (s + Jet::constant(2.5));
        let result = solve(&prepared, omega, f, q);
        check_uniform(&result, u, v, p, p / s, stress);
        power(
            result.power.structural_dissipation,
            work(v, v).map(|x| 0.17 * x),
        );
        power(result.power.wall_dissipation, work(p, p).map(|x| 0.12 * x));
        power(
            result.power.relaxation_dissipation,
            work(stress, stress).map(|x| x / (0.8 * 0.4)),
        );
        power(result.power.acoustic_to_relaxation, work(p, stress));
        power(result.power.structural_to_relaxation, [0.0; 3]);
        zero_balances(&result);
    }
}

#[test]
fn mutual_memory_and_gyroscopic_curvature_recover_manufactured_original_rows() {
    let n = 4;
    let mass = [
        1.1, 0.1, 0.0, 0.0, 0.1, 0.8, 0.0, 0.0, 0.0, 0.0, 0.7, 0.1, 0.0, 0.0, 0.1, 0.5,
    ];
    let damping = [
        0.3, 0.8, -0.2, 0.1, -0.6, 0.4, -0.03, -0.07, 0.2, 0.03, 0.02, 0.0, -0.1, 0.07, 0.0, 0.03,
    ];
    let stiffness = [
        4.0, 0.3, 0.0, 0.0, 0.3, 2.7, 0.0, 0.0, 0.0, 0.0, 0.9, -0.9, 0.0, 0.0, -0.9, 0.9,
    ];
    let coupling = [0.2, 0.03, -0.1, 0.07];
    let branches = [
        MaxwellDampingBranch {
            strain: vec![1.0, -0.2, 0.4, -0.1],
            stiffness: 0.8,
            relaxation_time_s: 0.6,
        },
        MaxwellDampingBranch {
            strain: vec![-0.3, 0.6, 0.1, 0.2],
            stiffness: 0.5,
            relaxation_time_s: 0.6,
        },
        MaxwellDampingBranch {
            strain: vec![0.0; n],
            stiffness: 2.0,
            relaxation_time_s: 0.1,
        },
    ];
    let prepared = Prepared::new(
        Source {
            structural: 2,
            coordinates: n,
            mass: &mass,
            damping: &damping,
            stiffness: &stiffness,
            coupling: &coupling,
        },
        &branches,
        Limits::default(),
        &mut || false,
    )
    .unwrap();
    let omega = 2.3;
    let s = Jet([C64::new(0.0, -omega), C64::new(0.0, -1.0), C64::ZERO]);
    let q: Vec<_> = (0..n)
        .map(|i| {
            Jet([
                C64::new(0.2 - 0.05 * i as f64, 0.01 * (i * i) as f64 - 0.08),
                C64::new(0.04 * i as f64 - 0.03, 0.02),
                C64::new(-0.02, 0.01 * i as f64),
            ])
        })
        .collect();
    let stress: Vec<_> = branches
        .iter()
        .map(|b| {
            let strain = q
                .iter()
                .zip(&b.strain)
                .fold(Jet::constant(0.0), |a, (&q, &b)| a + q.scale(b));
            s * strain.scale(b.stiffness) / (s + Jet::constant(b.relaxation_time_s.recip()))
        })
        .collect();
    let input: Vec<_> = (0..n)
        .map(|i| {
            let spatial = (0..n).fold(Jet::constant(0.0), |a, j| {
                a + (Jet::constant(stiffness[i * n + j])
                    + s.scale(damping[i * n + j])
                    + (s * s).scale(mass[i * n + j]))
                    * q[j]
            });
            branches
                .iter()
                .zip(&stress)
                .fold(spatial, |a, (b, &z)| a + z.scale(b.strain[i]))
        })
        .collect();
    let loads: [Vec<C64>; 3] = core::array::from_fn(|d| input.iter().map(|j| j.0[d]).collect());
    let factor = prepared.factor(omega, &mut || false).unwrap();
    let result = factor
        .solve_with_frequency_curvature(
            &loads[0][..2],
            &loads[0][2..],
            Some((&loads[1][..2], &loads[1][2..])),
            Some((&loads[2][..2], &loads[2][2..])),
            &mut || false,
        )
        .unwrap();
    let r = &result.frequency_response;
    for (i, &q) in q.iter().enumerate() {
        let v = s * q;
        if i < 2 {
            check(
                [
                    r.response.structural_displacement[i],
                    r.structural_displacement_derivative_per_rad_s[i],
                    result.structural_displacement_second_derivative[i],
                ],
                q,
            );
            check(
                [
                    r.response.structural_velocity[i],
                    r.structural_velocity_derivative_per_rad_s[i],
                    result.structural_velocity_second_derivative[i],
                ],
                v,
            );
        } else {
            check(
                [
                    r.response.pressure_potential_pa_s[i - 2],
                    r.pressure_potential_derivative_pa_s_per_rad_s[i - 2],
                    result.pressure_potential_second_derivative[i - 2],
                ],
                q,
            );
            check(
                [
                    r.response.pressure_pa[i - 2],
                    r.pressure_derivative_pa_per_rad_s[i - 2],
                    result.pressure_second_derivative[i - 2],
                ],
                v,
            );
        }
    }
    for (i, (&stress, branch)) in stress.iter().zip(&branches).enumerate() {
        check(
            [
                r.response.memory[i].stress,
                r.relaxation_stress_derivative_per_rad_s[i],
                result.relaxation_stress_second_derivative[i],
            ],
            stress,
        );
        power(
            result.memory[i].dissipation,
            work(stress, stress).map(|x| x / (branch.stiffness * branch.relaxation_time_s)),
        );
    }
    zero_balances(&result);
    assert_eq!(result.relaxation_stress_second_derivative[2], C64::ZERO);
    assert!(
        factor
            .solve_with_frequency_curvature(
                &loads[0][..2],
                &loads[0][2..],
                None,
                Some((&[], &[])),
                &mut || false
            )
            .is_err()
    );
    assert!(matches!(
        factor.solve_with_frequency_curvature(
            &loads[0][..2],
            &loads[0][2..],
            None,
            None,
            &mut || true
        ),
        Err(Error::Cancelled)
    ));
}
