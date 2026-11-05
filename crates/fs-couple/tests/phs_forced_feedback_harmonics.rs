//! G1/G3: imposed phase, generated harmonics, original memory and signed work.
use fs_couple::phs_orbit::{
    AutonomousPortFeedback, PeriodicPortInput, PhsOrbitError, audit_forced_feedback_orbit,
    solve_forced_feedback_harmonic_balance_adaptive_scaled_seeded,
    solve_forced_feedback_harmonic_balance_analytic_adaptive_scaled_seeded,
    solve_forced_feedback_harmonic_balance_analytic_scaled_seeded,
    solve_forced_feedback_harmonic_balance_scaled_seeded,
    solve_forced_feedback_shooting_variational,
    solve_forced_feedback_shooting_variational_analytic, solve_forced_harmonic_balance,
};
use fs_math::c64::C64;
use fs_orbit::{
    HbAnchor, HbBudget, HbOrbit, HbRefinementBudget, HbScales, OrbitError, OrbitProblem,
    ShootBudget, TimeDomainRhs, VariationalShootControls, VariationalShootError, solve_shooting,
};
use fs_phs::{PortHamiltonian, QuadraticStorage, RelaxationBranch, Storage};
use std::cell::Cell;

fn budget(harmonics: usize) -> HbBudget {
    HbBudget {
        harmonics,
        max_newton: 30,
        tolerance: 1e-11,
        deflate_equilibrium: false,
    }
}
fn scales(states: &[f64], omega: f64) -> HbScales {
    HbScales {
        states: states.to_vec(),
        equations: states.iter().map(|value| omega * value).collect(),
        angular_frequency: omega,
    }
}
fn refinement(states: &[f64], omega: f64, max_harmonics: usize) -> HbRefinementBudget {
    HbRefinementBudget {
        max_harmonics,
        residual_tolerances: states.iter().map(|v| 1e-8 * omega * v).collect(),
        state_tolerances: states.iter().map(|v| 1e-8 * v).collect(),
        frequency_tolerance: 1e-10,
    }
}
fn shoot_budget(steps_per_period: usize) -> ShootBudget {
    ShootBudget {
        steps_per_period,
        max_newton: 20,
        tolerance: 1e-11,
    }
}
fn shoot_controls(states: &[f64], omega: f64) -> VariationalShootControls {
    VariationalShootControls {
        state_scales: states.to_vec(),
        period_scale: core::f64::consts::TAU / omega,
        minimum_autonomous_excursion: 0.01,
        maximum_dimension: 16,
        maximum_payload_bytes: 2_000_000,
        maximum_work: 10_000_000_000,
        maximum_period_maps: 64,
    }
}
fn waveform(omega: f64, coeffs: Vec<Vec<C64>>) -> HbOrbit {
    HbOrbit {
        omega,
        coeffs,
        residual_trace: vec![],
        residual: 0.0,
        aft_samples: 0,
        sampling_defect: 0.0,
    }
}
fn seed(orbit: &HbOrbit, harmonics: usize) -> Vec<f64> {
    let coeffs = orbit
        .coeffs
        .iter()
        .map(|row| {
            let mut result = vec![C64::ZERO; harmonics + 1];
            for (to, from) in result.iter_mut().zip(row).take(2) {
                *to = from.scale(0.7);
            }
            result
        })
        .collect();
    waveform(orbit.omega, coeffs).packed_seed(&HbAnchor::Forced { omega: orbit.omega })
}
fn phase(angle: f64) -> C64 {
    C64::new(angle.cos(), angle.sin())
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:e} != {b:e} within {tolerance:e}"
    );
}
fn norm_squared(row: &[C64]) -> f64 {
    row[0].re.powi(2)
        + 2.0
            * row[1..]
                .iter()
                .map(|z| z.re * z.re + z.im * z.im)
                .sum::<f64>()
}

struct Cascade {
    calls: Cell<usize>,
}
impl AutonomousPortFeedback for Cascade {
    fn port_count(&self) -> usize {
        3
    }
    fn input(&self, _: &[f64], flow: &[f64], out: &mut [f64]) -> bool {
        self.calls.set(self.calls.get() + 1);
        out[0] = 0.0;
        out[1] = 0.45 * flow[0].powi(2);
        out[2] = 0.30 * flow[0].powi(3) - 0.20 * flow[1];
        true
    }
}
fn cascade_jacobian(states: &[f64; 3], flow: &[f64], out: &mut [f64]) -> bool {
    out.fill(0.0);
    out[3] = 0.90 * flow[0] / states[0];
    out[6] = 0.90 * flow[0].powi(2) / states[0];
    out[7] = -0.20 / states[1];
    true
}
fn cascade(
    states: &[f64; 3],
    omega: f64,
) -> (PortHamiltonian, PeriodicPortInput, HbOrbit, Vec<Vec<C64>>) {
    let rates = [1.3, 0.9, 2.1];
    let mut q = vec![0.0; 9];
    let mut r = vec![0.0; 9];
    let mut g = vec![0.0; 9];
    for i in 0..3 {
        q[4 * i] = 1.0 / states[i].powi(2);
        r[4 * i] = rates[i] * states[i].powi(2);
        g[4 * i] = states[i];
    }
    let system = PortHamiltonian::new(
        3,
        3,
        vec![0.0; 9],
        r,
        g,
        Box::new(QuadraticStorage::new(q, 3).unwrap()),
    )
    .unwrap();
    let (dc, amplitude, phi) = (0.17_f64, 0.72_f64, 0.38_f64);
    let first = phase(phi).scale(amplitude / 2.0);
    let forcing = vec![
        vec![
            C64::from_re(rates[0] * dc),
            C64::new(rates[0], omega) * first,
            C64::ZERO,
            C64::ZERO,
        ],
        vec![
            C64::from_re(-0.13),
            phase(-0.61).scale(0.15),
            C64::ZERO,
            C64::ZERO,
        ],
        vec![
            C64::from_re(0.08),
            C64::new(-0.09, 0.14),
            C64::ZERO,
            C64::ZERO,
        ],
    ];
    // Independent trigonometric identities, including the biased even lines.
    let square = [
        C64::from_re(dc * dc + amplitude * amplitude / 2.0),
        phase(phi).scale(dc * amplitude),
        phase(2.0 * phi).scale(amplitude * amplitude / 4.0),
        C64::ZERO,
    ];
    let cube = [
        C64::from_re(dc.powi(3) + 1.5 * dc * amplitude.powi(2)),
        phase(phi).scale(1.5 * dc * dc * amplitude + 0.375 * amplitude.powi(3)),
        phase(2.0 * phi).scale(0.75 * dc * amplitude.powi(2)),
        phase(3.0 * phi).scale(amplitude.powi(3) / 8.0),
    ];
    let mut normalized = vec![vec![C64::ZERO; 4]; 3];
    normalized[0] = vec![C64::from_re(dc), first, C64::ZERO, C64::ZERO];
    for h in 0..=3 {
        normalized[1][h] =
            (forcing[1][h] + square[h].scale(0.45)) / C64::new(rates[1], h as f64 * omega);
        normalized[2][h] = (forcing[2][h] + cube[h].scale(0.30) - normalized[1][h].scale(0.20))
            / C64::new(rates[2], h as f64 * omega);
    }
    let physical = normalized
        .iter()
        .zip(states)
        .map(|(row, &scale)| row.iter().map(|z| z.scale(scale)).collect())
        .collect();
    (
        system,
        PeriodicPortInput::new(omega, forcing).unwrap(),
        waveform(omega, physical),
        normalized,
    )
}

#[test]
fn phased_multiport_forcing_generates_exact_even_odd_and_memory_harmonics_across_units() {
    let omega = 2.4;
    for unit in [1e-9, 1.0, 1e9] {
        let states = [1e-3 * unit, 0.02 * unit, 5.0 * unit];
        let (system, input, expected, normalized) = cascade(&states, omega);
        let feedback = Cascade {
            calls: Cell::new(0),
        };
        let tangent_calls = Cell::new(0);
        let jacobian = |_: &[f64], flow: &[f64], out: &mut [f64]| {
            tangent_calls.set(tangent_calls.get() + 1);
            cascade_jacobian(&states, flow, out)
        };
        let fixed = solve_forced_feedback_harmonic_balance_scaled_seeded(
            &system,
            &input,
            &feedback,
            &jacobian,
            seed(&expected, 3),
            &budget(3),
            &scales(&states, omega),
        )
        .unwrap();
        let adaptive = solve_forced_feedback_harmonic_balance_adaptive_scaled_seeded(
            &system,
            &input,
            &feedback,
            &jacobian,
            seed(&expected, 1),
            &budget(1),
            &refinement(&states, omega, 7),
            &scales(&states, omega),
        )
        .unwrap();
        assert!(tangent_calls.get() > 0);
        assert!(adaptive.refinement_trace[0].omitted_residual_ratio > 1.0);
        for result in [&fixed.orbit, &adaptive.orbit] {
            assert_eq!(result.omega.to_bits(), omega.to_bits());
            assert_eq!(result.coeffs.len(), 3);
            for i in 0..3 {
                for (h, actual) in result.coeffs[i].iter().enumerate() {
                    let reference = expected.coeffs[i].get(h).copied().unwrap_or(C64::ZERO);
                    assert!(
                        (*actual - reference).abs() / states[i] < 2e-10,
                        "state {i}, harmonic {h}"
                    );
                }
            }
        }
        assert!(adaptive.orbit.coeffs[1][2].abs() / states[1] > 0.01);
        assert!(adaptive.orbit.coeffs[2][3].abs() / states[2] > 0.001);
        let audit =
            audit_forced_feedback_orbit(&system, &input, &feedback, &adaptive.orbit, 128).unwrap();
        let energy = 0.5 * normalized.iter().map(|row| norm_squared(row)).sum::<f64>();
        let loss = normalized
            .iter()
            .zip([1.3, 0.9, 2.1])
            .map(|(row, rate)| rate * norm_squared(row))
            .sum::<f64>();
        close(audit.cycle.mean_stored_energy_j, energy, 2e-10);
        close(audit.cycle.mean_resistance_loss_w, loss, 2e-10);
        close(
            audit.mean_external_drive_power_w + audit.mean_feedback_power_w,
            loss,
            2e-10,
        );
        assert!(audit.cycle.relative_power_residual < 2e-10);
    }
}

struct PassiveDrag {
    gamma: f64,
    offset: f64,
}
impl AutonomousPortFeedback for PassiveDrag {
    fn port_count(&self) -> usize {
        1
    }
    fn input(&self, _: &[f64], flow: &[f64], out: &mut [f64]) -> bool {
        out[0] = -self.gamma * flow[0] * (flow[0] + self.offset).powi(2);
        true
    }
}

struct DrivenDrag<'a> {
    system: &'a PortHamiltonian,
    drag: &'a PassiveDrag,
    omega: f64,
    amplitude: f64,
    phase: f64,
    mass: f64,
}
impl OrbitProblem for DrivenDrag<'_> {
    fn dim(&self) -> usize {
        self.system.state_dim()
    }
    fn autonomous(&self) -> bool {
        false
    }
    fn island(&self, time: f64, state: &[f64], out: &mut [f64]) {
        // Independent direct physical source law, without the Fourier/feedback adapter.
        let velocity = state[1] / self.mass;
        let force = self.amplitude * (self.omega * time + self.phase).cos()
            - self.drag.gamma * velocity * (velocity + self.drag.offset).powi(2);
        self.system
            .vector_field_from_effort_into(&self.system.effort(state), &[force], out)
            .unwrap();
    }
    fn time_domain_rhs(&self) -> Option<TimeDomainRhs<Self>> {
        Some(Self::island)
    }
}

#[test]
fn passive_drag_harmonics_retain_maxwell_memory_and_match_independent_shooting() {
    let (mass, omega, amplitude, phi) = (0.8_f64, 1.1_f64, 0.25_f64, 0.27_f64);
    let system = fs_phs::mass_spring_damper(mass, 1.4, 0.2)
        .unwrap()
        .with_relaxation_branches(vec![RelaxationBranch {
            projection: vec![1.0, 0.0],
            stiffness: 0.35,
            relaxation_time_s: 0.4,
        }])
        .unwrap();
    let input = PeriodicPortInput::new(
        omega,
        vec![vec![C64::ZERO, phase(phi).scale(amplitude / 2.0)]],
    )
    .unwrap();
    let linear = solve_forced_harmonic_balance(&system, &input, 0.3, &budget(1)).unwrap();
    let states = [0.5; 3];
    for offset in [0.0, 0.35] {
        let drag = PassiveDrag { gamma: 0.8, offset };
        let jacobian = |_: &[f64], flow: &[f64], out: &mut [f64]| {
            out.fill(0.0);
            out[1] = -drag.gamma
                * (3.0 * flow[0].powi(2) + 4.0 * offset * flow[0] + offset * offset)
                / mass;
            true
        };
        let result = solve_forced_feedback_harmonic_balance_adaptive_scaled_seeded(
            &system,
            &input,
            &drag,
            &jacobian,
            linear.packed_seed(&HbAnchor::Forced { omega }),
            &budget(1),
            &refinement(&states, omega, 31),
            &scales(&states, omega),
        )
        .unwrap();
        assert!(result.orbit.coeffs[1][3].abs() > 1e-4);
        if offset == 0.0 {
            assert!(result.orbit.coeffs[1][2].abs() < 1e-10);
        } else {
            assert!(result.orbit.coeffs[1][2].abs() > 1e-4);
        }
        assert!(
            result.orbit.coeffs[2][1].abs() > 0.01,
            "the original relaxation state was lost"
        );
        let a = audit_forced_feedback_orbit(&system, &input, &drag, &result.orbit, 256).unwrap();
        let b = audit_forced_feedback_orbit(&system, &input, &drag, &result.orbit, 512).unwrap();
        assert!(a.mean_external_drive_power_w > 0.0 && a.mean_feedback_power_w < 0.0);
        assert!(a.cycle.mean_resistance_loss_w > 0.0 && a.cycle.relative_power_residual < 1e-8);
        close(a.mean_feedback_power_w, b.mean_feedback_power_w, 1e-11);
        let drag_loss = (0..512)
            .map(|j| {
                let velocity = result
                    .orbit
                    .sample(1, core::f64::consts::TAU * j as f64 / 512.0)
                    / mass;
                drag.gamma * velocity.powi(2) * (velocity + offset).powi(2)
            })
            .sum::<f64>()
            / 512.0;
        close(-a.mean_feedback_power_w, drag_loss, 1e-11);
        close(
            a.mean_external_drive_power_w,
            a.cycle.mean_resistance_loss_w + drag_loss,
            1e-9,
        );
        for shift in [0.0, 0.71] {
            let direct = DrivenDrag {
                system: &system,
                drag: &drag,
                omega,
                amplitude,
                phase: phi + shift,
                mass,
            };
            let reference = solve_shooting(
                &direct,
                &[0.0; 3],
                core::f64::consts::TAU / omega,
                &shoot_budget(1024),
            )
            .unwrap();
            let shifted_input = PeriodicPortInput::new(
                omega,
                vec![vec![C64::ZERO, phase(phi + shift).scale(amplitude / 2.0)]],
            )
            .unwrap();
            let shot = solve_forced_feedback_shooting_variational(
                &system,
                &shifted_input,
                &drag,
                &jacobian,
                &[0.0; 3],
                &shoot_budget(1024),
                &shoot_controls(&states, omega),
                || false,
            )
            .unwrap();
            assert_eq!(
                shot.orbit.period.to_bits(),
                (core::f64::consts::TAU / omega).to_bits()
            );
            assert!(shot.period_derivative.is_none());
            assert_eq!(shot.monodromy.len(), 9);
            assert_eq!(shot.orbit.multipliers.len(), 3);
            for (i, &actual) in reference.x0.iter().enumerate() {
                close(result.orbit.sample(i, shift), actual, 2e-8);
                close(shot.orbit.x0[i], actual, 2e-9);
            }
            for &(re, im) in &reference.multipliers {
                assert!(
                    shot.orbit
                        .multipliers
                        .iter()
                        .any(|&(a, b)| (a - re).hypot(b - im) < 2e-6)
                );
            }
        }
    }
}

struct NonlinearStorage {
    beta: f64,
}
impl Storage for NonlinearStorage {
    fn hamiltonian(&self, x: &[f64]) -> f64 {
        0.5 * x[0].powi(2) + 0.25 * self.beta * x[0].powi(4)
    }
    fn gradient(&self, x: &[f64], out: &mut [f64]) {
        out[0] = x[0] + self.beta * x[0].powi(3);
    }
}
struct NonlinearFeedback {
    damping: f64,
}
impl AutonomousPortFeedback for NonlinearFeedback {
    fn port_count(&self) -> usize {
        1
    }
    fn input(&self, state: &[f64], flow: &[f64], out: &mut [f64]) -> bool {
        out[0] = self.damping * (flow[0] - state[0]);
        true
    }
}

#[test]
fn explicit_hessian_and_total_feedback_chain_keep_nonlinear_storage_and_exact_work_split() {
    let (beta, damping, omega, amplitude) = (0.4_f64, 0.7_f64, 1.3_f64, 0.8_f64);
    let system = PortHamiltonian::new(
        1,
        1,
        vec![0.0],
        vec![damping],
        vec![1.0],
        Box::new(NonlinearStorage { beta }),
    )
    .unwrap();
    let feedback = NonlinearFeedback { damping };
    let exact = waveform(
        omega,
        vec![vec![C64::ZERO, phase(0.32).scale(amplitude / 2.0)]],
    );
    let input = PeriodicPortInput::new(
        omega,
        vec![vec![
            C64::ZERO,
            exact.coeffs[0][1] * C64::new(damping, omega),
        ]],
    )
    .unwrap();
    let jacobian = |state: &[f64], _: &[f64], out: &mut [f64]| {
        out[0] = 3.0 * damping * beta * state[0].powi(2);
        true
    };
    let hessian = |state: &[f64], direction: &[f64], out: &mut [f64]| {
        out[0] = (1.0 + 3.0 * beta * state[0].powi(2)) * direction[0];
        true
    };
    assert!(matches!(
        solve_forced_feedback_harmonic_balance_scaled_seeded(
            &system,
            &input,
            &feedback,
            &jacobian,
            seed(&exact, 1),
            &budget(1),
            &scales(&[1.0], omega)
        ),
        Err(PhsOrbitError::Invalid { .. })
    ));
    let fixed = solve_forced_feedback_harmonic_balance_analytic_scaled_seeded(
        &system,
        &input,
        &feedback,
        &jacobian,
        seed(&exact, 1),
        &budget(1),
        &scales(&[1.0], omega),
        &hessian,
    )
    .unwrap();
    let adaptive = solve_forced_feedback_harmonic_balance_analytic_adaptive_scaled_seeded(
        &system,
        &input,
        &feedback,
        &jacobian,
        seed(&exact, 1),
        &budget(1),
        &refinement(&[1.0], omega, 3),
        &scales(&[1.0], omega),
        &hessian,
    )
    .unwrap();
    for result in [&fixed.orbit, &adaptive.orbit] {
        assert!((result.coeffs[0][1] - exact.coeffs[0][1]).abs() < 1e-10);
    }
    let audit =
        audit_forced_feedback_orbit(&system, &input, &feedback, &adaptive.orbit, 64).unwrap();
    let x2 = amplitude.powi(2) / 2.0;
    let x4 = 3.0 * amplitude.powi(4) / 8.0;
    let x6 = 5.0 * amplitude.powi(6) / 16.0;
    close(
        audit.cycle.mean_stored_energy_j,
        0.5 * x2 + 0.25 * beta * x4,
        1e-10,
    );
    close(
        audit.mean_external_drive_power_w,
        damping * (x2 + beta * x4),
        1e-10,
    );
    close(
        audit.mean_feedback_power_w,
        damping * beta * (x4 + beta * x6),
        1e-10,
    );
    close(
        audit.cycle.mean_resistance_loss_w,
        damping * (x2 + 2.0 * beta * x4 + beta * beta * x6),
        1e-10,
    );
    assert!(audit.cycle.relative_power_residual < 1e-10);
    assert!(matches!(
        solve_forced_feedback_shooting_variational(
            &system,
            &input,
            &feedback,
            &jacobian,
            &[0.0],
            &shoot_budget(512),
            &shoot_controls(&[1.0], omega),
            || false
        ),
        Err(PhsOrbitError::Variational(
            VariationalShootError::TangentUnavailable
        ))
    ));
    let shot = solve_forced_feedback_shooting_variational_analytic(
        &system,
        &input,
        &feedback,
        &jacobian,
        &[0.0],
        &shoot_budget(512),
        &shoot_controls(&[1.0], omega),
        &hessian,
        || false,
    )
    .unwrap();
    close(shot.orbit.x0[0], exact.sample(0, 0.0), 5e-10);
    // The total physical tangent reduces exactly to -damping; compare the
    // discrete Floquet multiplier against the independent RK4 polynomial.
    let z = -damping * core::f64::consts::TAU / omega / 512.0;
    let multiplier = (1.0 + z + z * z / 2.0 + z.powi(3) / 6.0 + z.powi(4) / 24.0).powi(512);
    close(shot.monodromy[0], multiplier, 2e-12);
    close(shot.orbit.multipliers[0].0, multiplier, 2e-12);
    close(shot.orbit.multipliers[0].1, 0.0, 1e-15);
    let refuse = |_: &[f64], _: &[f64], _: &mut [f64]| false;
    assert!(matches!(
        solve_forced_feedback_harmonic_balance_analytic_scaled_seeded(
            &system,
            &input,
            &feedback,
            &jacobian,
            seed(&exact, 1),
            &budget(1),
            &scales(&[1.0], omega),
            &refuse
        ),
        Err(PhsOrbitError::Orbit(OrbitError::TangentRefused { .. }))
    ));
    assert!(matches!(
        solve_forced_feedback_shooting_variational_analytic(
            &system,
            &input,
            &feedback,
            &jacobian,
            &[0.0],
            &shoot_budget(512),
            &shoot_controls(&[1.0], omega),
            &refuse,
            || false
        ),
        Err(PhsOrbitError::Variational(
            VariationalShootError::SourceRefused { tangent: true, .. }
        ))
    ));
}

#[test]
fn unresolved_external_lines_and_refusing_sources_never_publish_partial_orbits() {
    let states = [1.0; 3];
    let (system, input, expected, _) = cascade(&states, 2.4);
    let feedback = Cascade {
        calls: Cell::new(0),
    };
    let jacobian_calls = Cell::new(0);
    let jacobian = |_: &[f64], flow: &[f64], out: &mut [f64]| {
        jacobian_calls.set(jacobian_calls.get() + 1);
        cascade_jacobian(&states, flow, out)
    };
    let mut unresolved = vec![vec![C64::ZERO; 10]; 3];
    unresolved[2][9] = C64::new(0.0, f64::from_bits(1));
    let unresolved = PeriodicPortInput::new(2.4, unresolved).unwrap();
    assert!(matches!(
        solve_forced_feedback_harmonic_balance_scaled_seeded(
            &system,
            &unresolved,
            &feedback,
            &jacobian,
            seed(&expected, 1),
            &budget(1),
            &scales(&states, 2.4)
        ),
        Err(PhsOrbitError::Invalid { .. })
    ));
    assert_eq!(feedback.calls.get(), 0);
    assert_eq!(jacobian_calls.get(), 0);
    assert!(matches!(
        solve_forced_feedback_shooting_variational(
            &system,
            &unresolved,
            &feedback,
            &jacobian,
            &[0.0; 3],
            &shoot_budget(16),
            &shoot_controls(&states, 2.4),
            || false
        ),
        Err(PhsOrbitError::Invalid { .. })
    ));
    let mut too_small = shoot_controls(&states, 2.4);
    too_small.maximum_work = 0;
    assert!(matches!(
        solve_forced_feedback_shooting_variational(
            &system,
            &input,
            &feedback,
            &jacobian,
            &[0.0; 3],
            &shoot_budget(64),
            &too_small,
            || false
        ),
        Err(PhsOrbitError::Variational(
            VariationalShootError::Budget { .. }
        ))
    ));
    assert!(matches!(
        solve_forced_feedback_shooting_variational(
            &system,
            &input,
            &feedback,
            &jacobian,
            &[0.0; 3],
            &shoot_budget(64),
            &shoot_controls(&states, 2.4),
            || true
        ),
        Err(PhsOrbitError::Variational(VariationalShootError::Cancelled))
    ));
    assert_eq!(feedback.calls.get(), 0);
    assert_eq!(jacobian_calls.get(), 0);
    let incomplete = |_: &[f64], _: &[f64], out: &mut [f64]| {
        out[0] = 0.0;
        true
    };
    assert!(matches!(
        solve_forced_feedback_harmonic_balance_scaled_seeded(
            &system,
            &input,
            &feedback,
            &incomplete,
            seed(&expected, 3),
            &budget(3),
            &scales(&states, 2.4)
        ),
        Err(PhsOrbitError::Orbit(OrbitError::TangentRefused { .. }))
    ));
    assert!(matches!(
        solve_forced_feedback_shooting_variational(
            &system,
            &input,
            &feedback,
            &incomplete,
            &[0.0; 3],
            &shoot_budget(64),
            &shoot_controls(&states, 2.4),
            || false
        ),
        Err(PhsOrbitError::Variational(
            VariationalShootError::SourceRefused { tangent: true, .. }
        ))
    ));
    struct Refuse;
    impl AutonomousPortFeedback for Refuse {
        fn port_count(&self) -> usize {
            3
        }
        fn input(&self, _: &[f64], _: &[f64], _: &mut [f64]) -> bool {
            false
        }
    }
    assert!(matches!(
        solve_forced_feedback_harmonic_balance_scaled_seeded(
            &system,
            &input,
            &Refuse,
            &jacobian,
            seed(&expected, 3),
            &budget(3),
            &scales(&states, 2.4)
        ),
        Err(PhsOrbitError::SourceRefused)
    ));
    assert!(matches!(
        solve_forced_feedback_shooting_variational(
            &system,
            &input,
            &Refuse,
            &jacobian,
            &[0.0; 3],
            &shoot_budget(64),
            &shoot_controls(&states, 2.4),
            || false
        ),
        Err(PhsOrbitError::SourceRefused)
    ));
    let polls = Cell::new(0);
    assert!(matches!(
        solve_forced_feedback_shooting_variational(
            &system,
            &input,
            &feedback,
            &jacobian,
            &[0.0; 3],
            &shoot_budget(64),
            &shoot_controls(&states, 2.4),
            || {
                polls.set(polls.get() + 1);
                polls.get() > 16
            }
        ),
        Err(PhsOrbitError::Variational(VariationalShootError::Cancelled))
    ));
    assert!(matches!(
        audit_forced_feedback_orbit(&system, &input, &Refuse, &expected, 64),
        Err(PhsOrbitError::SourceRefused)
    ));
    let shifted_clock = input.at_frequency(2.5).unwrap();
    assert!(matches!(
        audit_forced_feedback_orbit(&system, &shifted_clock, &feedback, &expected, 64),
        Err(PhsOrbitError::Invalid { .. })
    ));
}
