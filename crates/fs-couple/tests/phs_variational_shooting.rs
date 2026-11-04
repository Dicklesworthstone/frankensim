//! Complete original pHS sources on the shared discrete variational map.
use fs_couple::phs_orbit::{
    PeriodicPortInput, PhsOrbitError, audit_forced_orbit, solve_forced_shooting_variational,
    solve_forced_shooting_variational_analytic,
};
use fs_math::c64::C64;
use fs_orbit::{HbOrbit, ShootBudget, VariationalShootControls, VariationalShootError};
use fs_phs::{PortHamiltonian, RelaxationBranch, Storage};
use std::{cell::Cell, rc::Rc};

fn budget(steps: usize) -> ShootBudget {
    ShootBudget {
        steps_per_period: steps,
        max_newton: 20,
        tolerance: 1e-11,
    }
}
fn controls(scales: &[f64], period: f64) -> VariationalShootControls {
    VariationalShootControls {
        state_scales: scales.to_vec(),
        period_scale: period,
        minimum_autonomous_excursion: 0.01,
        maximum_dimension: 16,
        maximum_payload_bytes: 2_000_000,
        maximum_work: 10_000_000_000,
        maximum_period_maps: 64,
    }
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e} within {tolerance:e}"
    );
}
fn amplification(pole: C64, period: f64, steps: usize) -> C64 {
    // Independent closed-form RK4 stability polynomial for a constant mode;
    // no integration of the physical source or numerical variational equations.
    let z = pole.scale(period / steps as f64);
    let z2 = z * z;
    let one = C64::ONE + z + z2.scale(0.5) + (z2 * z).scale(1. / 6.) + (z2 * z2).scale(1. / 24.);
    (0..steps).fold(C64::ONE, |value, _| value * one)
}
fn check_spectrum(actual: &[(f64, f64)], expected: &[C64]) {
    assert_eq!(actual.len(), expected.len());
    for &value in expected {
        let error = actual
            .iter()
            .map(|&(re, im)| (C64::new(re, im) - value).abs())
            .fold(f64::INFINITY, f64::min);
        assert!(error < 2e-10, "{actual:?} misses {value:?}: {error:e}");
    }
}
fn waveform(omega: f64, coeffs: Vec<Vec<C64>>) -> HbOrbit {
    HbOrbit {
        omega,
        coeffs,
        residual_trace: vec![],
        residual: 0.,
        aft_samples: 0,
        sampling_defect: 0.,
    }
}

#[test]
fn original_mechanical_forcing_clock_momentum_and_power_survive_physical_scaling() {
    let (m, k, c, w) = (2.0_f64, 200.0_f64, 4.0_f64, 9.0_f64);
    let system = fs_phs::mass_spring_damper(m, k, c).unwrap();
    let period = core::f64::consts::TAU / w;
    let steps = 512;
    let raw = [C64::from_re(1.), C64::new(0.4, -0.1), C64::new(0.05, 0.02)];
    let mut previous: Option<Vec<f64>> = None;
    for scale in [1e-50, 1., 1e50] {
        let force: Vec<_> = raw.iter().map(|z| z.scale(scale)).collect();
        let input = PeriodicPortInput::new(w, vec![force.clone()]).unwrap();
        let q: Vec<_> = force
            .iter()
            .enumerate()
            .map(|(n, &f)| {
                let omega = n as f64 * w;
                f / C64::new(k - m * omega * omega, c * omega)
            })
            .collect();
        let p: Vec<_> = q
            .iter()
            .enumerate()
            .map(|(n, &q)| q * C64::new(0., m * n as f64 * w))
            .collect();
        let exact = waveform(w, vec![q.clone(), p]);
        let result = solve_forced_shooting_variational(
            &system,
            &input,
            &[0., 0.],
            &budget(steps),
            &controls(&[0.02 * scale, 0.4 * scale], period),
            || false,
        )
        .unwrap();
        assert_eq!(result.orbit.period.to_bits(), period.to_bits());
        close(
            result.orbit.x0[0] / scale,
            exact.sample(0, 0.) / scale,
            2e-9,
        );
        close(
            result.orbit.x0[1] / scale,
            exact.sample(1, 0.) / scale,
            3e-8,
        );
        let decay = c / (2. * m);
        let frequency = (k / m - decay * decay).sqrt();
        check_spectrum(
            &result.orbit.multipliers,
            &[
                amplification(C64::new(-decay, frequency), period, steps),
                amplification(C64::new(-decay, -frequency), period, steps),
            ],
        );
        if let Some(previous) = &previous {
            for (&a, &b) in result.monodromy.iter().zip(previous) {
                close(a, b, 3e-12);
            }
        }
        previous = Some(result.monodromy.clone());
        assert_eq!(result.rhs_evaluations, 4 * steps * result.period_maps);
        assert!(result.period_derivative.is_none());
        let x = &result.orbit.x0;
        let energy = system.hamiltonian(x) / (scale * scale);
        close(
            energy,
            0.5 * k * (x[0] / scale).powi(2) + 0.5 / m * (x[1] / scale).powi(2),
            1e-12,
        );
        // Independent source-law audit on the analytic reference waveform.
        let audit = audit_forced_orbit(&system, &input, &exact, 128).unwrap();
        let loss = 2.
            * c
            * q.iter()
                .enumerate()
                .skip(1)
                .map(|(n, q)| (n as f64 * w).powi(2) * (q.scale(1. / scale)).norm_sq())
                .sum::<f64>();
        close(audit.mean_resistance_loss_w / (scale * scale), loss, 1e-12);
        close(audit.mean_supplied_power_w / (scale * scale), loss, 1e-12);
        assert!(audit.relative_power_residual < 1e-12);
    }
}

#[test]
fn complete_maxwell_history_and_all_floquet_modes_use_the_original_material_tangent() {
    let (k, kb, tau, w, force) = (60_000.0_f64, 50_000.0_f64, 1. / 300., 240.0_f64, 2.0_f64);
    let system = fs_phs::mass_spring_damper(1., k, 0.)
        .unwrap()
        .with_relaxation_branches(vec![RelaxationBranch {
            projection: vec![1., 0.],
            stiffness: kb,
            relaxation_time_s: tau,
        }])
        .unwrap();
    let input = PeriodicPortInput::new(w, vec![vec![C64::ZERO, C64::from_re(force / 2.)]]).unwrap();
    let period = core::f64::consts::TAU / w;
    let steps = 1024;
    let result = solve_forced_shooting_variational(
        &system,
        &input,
        &[0.; 3],
        &budget(steps),
        &controls(&[1e-4, 0.03, 0.03], period),
        || false,
    )
    .unwrap();
    assert_eq!(result.orbit.x0.len(), 3);
    assert_eq!(result.monodromy.len(), 9);
    let s = C64::new(0., w);
    let arm = C64::ONE + s.scale(tau);
    let q = C64::from_re(force) / (C64::from_re(k - w * w) + (s.scale(tau) / arm).scale(kb));
    let p = s * q;
    let memory = q.scale(kb.sqrt()) / arm;
    for (&actual, expected) in result.orbit.x0.iter().zip([q, p, memory]) {
        close(actual, expected.re, 2e-10 * expected.abs().max(1e-4));
    }
    let a = (87_500.0_f64).sqrt();
    check_spectrum(
        &result.orbit.multipliers,
        &[
            amplification(C64::new(-50., a), period, steps),
            amplification(C64::new(-50., -a), period, steps),
            amplification(C64::from_re(-200.), period, steps),
        ],
    );
    let x = &result.orbit.x0;
    let strain = kb.sqrt() * x[0] - x[2];
    close(
        system.hamiltonian(x),
        0.5 * (k * x[0] * x[0] + x[1] * x[1] + strain * strain),
        1e-15,
    );
    let exact = waveform(
        w,
        vec![
            vec![C64::ZERO, q.scale(0.5)],
            vec![C64::ZERO, p.scale(0.5)],
            vec![C64::ZERO, memory.scale(0.5)],
        ],
    );
    let audit = audit_forced_orbit(&system, &input, &exact, 128).unwrap();
    let heat = (memory - q.scale(kb.sqrt())).norm_sq() / (2. * tau);
    close(audit.mean_resistance_loss_w, heat, 1e-13);
    close(audit.mean_supplied_power_w, heat, 1e-13);
}

struct Cubic {
    m: f64,
    k: f64,
    alpha: f64,
    calls: Rc<Cell<usize>>,
}
impl Storage for Cubic {
    fn hamiltonian(&self, x: &[f64]) -> f64 {
        0.5 * self.k * x[0] * x[0] + 0.25 * self.alpha * x[0].powi(4) + 0.5 * x[1] * x[1] / self.m
    }
    fn gradient(&self, x: &[f64], out: &mut [f64]) {
        self.calls.set(self.calls.get() + 1);
        out[0] = self.k * x[0] + self.alpha * x[0].powi(3);
        out[1] = x[1] / self.m;
    }
    fn quadratic_part(&self) -> Option<fs_phs::linear_exact::QuadraticPart> {
        Some(fs_phs::linear_exact::QuadraticPart {
            matrix: vec![self.k, 0., 0., 1. / self.m],
            exact: false,
        })
    }
}

#[test]
fn declared_nonlinear_hessian_preserves_original_effort_law_without_trajectory_probes() {
    let (m, k, alpha, c, w, a) = (2.0, 8.0, 5.0, 0.6, 2.0, 0.3);
    let calls = Rc::new(Cell::new(0));
    let system = PortHamiltonian::new(
        2,
        1,
        vec![0., 1., -1., 0.],
        vec![0., 0., 0., c],
        vec![0., 1.],
        Box::new(Cubic {
            m,
            k,
            alpha,
            calls: Rc::clone(&calls),
        }),
    )
    .unwrap();
    let q1 = C64::from_re(a / 2.);
    let p1 = q1 * C64::new(0., m * w);
    let f1 = q1 * C64::new(k - m * w * w + 0.75 * alpha * a * a, c * w);
    let f3 = (q1 * q1 * q1).scale(alpha);
    let input = PeriodicPortInput::new(w, vec![vec![C64::ZERO, f1, C64::ZERO, f3]]).unwrap();
    let period = core::f64::consts::TAU / w;
    let limits = controls(&[a, m * w * a], period);
    calls.set(0);
    assert!(matches!(
        solve_forced_shooting_variational(&system, &input, &[a, 0.], &budget(512), &limits, || {
            false
        }),
        Err(PhsOrbitError::Variational(
            VariationalShootError::TangentUnavailable
        ))
    ));
    assert_eq!(
        calls.get(),
        0,
        "a nonlinear quadratic PART cannot admit a flow tangent"
    );
    let actions = Cell::new(0);
    let hessian = |x: &[f64], v: &[f64], out: &mut [f64]| {
        actions.set(actions.get() + 1);
        out[0] = (k + 3. * alpha * x[0] * x[0]) * v[0];
        out[1] = v[1] / m;
        true
    };
    let result = solve_forced_shooting_variational_analytic(
        &system,
        &input,
        &[0.9 * a, 0.],
        &budget(512),
        &limits,
        &hessian,
        || false,
    )
    .unwrap();
    close(result.orbit.x0[0], a, 2e-8);
    close(result.orbit.x0[1], 0., 2e-8);
    assert_eq!(calls.get(), result.rhs_evaluations);
    assert_eq!(actions.get(), 2 * result.jacobian_evaluations);
    let exact = waveform(
        w,
        vec![
            vec![C64::ZERO, q1, C64::ZERO, C64::ZERO],
            vec![C64::ZERO, p1, C64::ZERO, C64::ZERO],
        ],
    );
    let audit = audit_forced_orbit(&system, &input, &exact, 128).unwrap();
    close(audit.mean_resistance_loss_w, 0.5 * c * w * w * a * a, 1e-13);
    close(audit.mean_supplied_power_w, 0.5 * c * w * w * a * a, 1e-13);
    close(
        audit.mean_stored_energy_j,
        0.25 * (k + m * w * w) * a * a + 3. * alpha * a.powi(4) / 32.,
        1e-13,
    );
    let partial = |_: &[f64], v: &[f64], out: &mut [f64]| {
        out[0] = k * v[0];
        true
    };
    assert!(matches!(
        solve_forced_shooting_variational_analytic(
            &system,
            &input,
            &[a, 0.],
            &budget(512),
            &limits,
            &partial,
            || false
        ),
        Err(PhsOrbitError::Variational(
            VariationalShootError::SourceRefused { tangent: true, .. }
        ))
    ));
}

#[test]
fn physical_bandwidth_extra_resources_and_cancellation_are_admitted_before_publication() {
    let system = fs_phs::mass_spring_damper(2., 200., 4.).unwrap();
    let period = core::f64::consts::TAU / 9.;
    let input = PeriodicPortInput::new(9., vec![vec![C64::ZERO, C64::from_re(0.5)]]).unwrap();
    let base = controls(&[0.02, 0.4], period);
    let wrong = PeriodicPortInput::new(9., vec![]).unwrap();
    assert!(matches!(
        solve_forced_shooting_variational(&system, &wrong, &[0., 0.], &budget(128), &base, || {
            false
        }),
        Err(PhsOrbitError::Invalid { .. })
    ));
    let mut force = vec![C64::ZERO; 6];
    force[5] = C64::ONE;
    let high = PeriodicPortInput::new(9., vec![force]).unwrap();
    assert!(matches!(
        solve_forced_shooting_variational(&system, &high, &[0., 0.], &budget(8), &base, || false),
        Err(PhsOrbitError::Invalid { .. })
    ));
    for (bytes, work) in [(479, 10_000_000_000), (2_000_000, 1)] {
        let limited = VariationalShootControls {
            maximum_payload_bytes: bytes,
            maximum_work: work,
            ..base.clone()
        };
        assert!(matches!(
            solve_forced_shooting_variational(
                &system,
                &input,
                &[0., 0.],
                &budget(128),
                &limited,
                || false
            ),
            Err(PhsOrbitError::Variational(
                VariationalShootError::Budget { .. }
            ))
        ));
    }
    let mut polls = 0;
    assert!(matches!(
        solve_forced_shooting_variational(&system, &input, &[0., 0.], &budget(128), &base, || {
            polls += 1;
            polls > 40
        }),
        Err(PhsOrbitError::Variational(VariationalShootError::Cancelled))
    ));
    let first =
        solve_forced_shooting_variational(&system, &input, &[0., 0.], &budget(128), &base, || {
            false
        })
        .unwrap();
    let second =
        solve_forced_shooting_variational(&system, &input, &[0., 0.], &budget(128), &base, || {
            false
        })
        .unwrap();
    assert_eq!(first.orbit.x0, second.orbit.x0);
    assert_eq!(first.monodromy, second.monodromy);
    assert_eq!(first.orbit.multipliers, second.orbit.multipliers);

    // Reject malformed caller scale storage before cloning/preparing a source
    // or asking an explicitly supplied nonlinear tangent to do any work.
    let calls = Rc::new(Cell::new(0));
    let counted = PortHamiltonian::new(
        2,
        1,
        vec![0., 1., -1., 0.],
        vec![0., 0., 0., 4.],
        vec![0., 1.],
        Box::new(Cubic {
            m: 2.,
            k: 200.,
            alpha: 1.,
            calls: Rc::clone(&calls),
        }),
    )
    .unwrap();
    calls.set(0);
    let actions = Cell::new(0);
    let hessian = |_: &[f64], _: &[f64], _: &mut [f64]| {
        actions.set(actions.get() + 1);
        false
    };
    for state_scales in [
        vec![],
        vec![1.],
        vec![1., 1., 1.],
        vec![0., 1.],
        vec![1., f64::NAN],
        vec![1., f64::INFINITY],
    ] {
        let malformed = VariationalShootControls {
            state_scales,
            ..base.clone()
        };
        assert!(matches!(
            solve_forced_shooting_variational_analytic(
                &counted,
                &input,
                &[0., 0.],
                &budget(128),
                &malformed,
                &hessian,
                || false
            ),
            Err(PhsOrbitError::Invalid { .. })
        ));
    }
    assert_eq!(calls.get(), 0);
    assert_eq!(actions.get(), 0);
}
