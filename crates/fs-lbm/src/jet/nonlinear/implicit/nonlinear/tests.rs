use super::*;
use fs_material::fluid_transport::rheology::{
    CarreauYasuda, ThermalCarreauYasuda, ViscosityShiftMode,
};
use fs_material::visco::temperature::{TemperatureShift, TemperatureShiftModel};

fn liquid(values: [f64; 4], ceiling: f64) -> CarreauSlenderJetParameters {
    let [mu0, mu_inf, time, index] = values;
    CarreauSlenderJetParameters {
        density_kg_m3: 2.0,
        surface_tension_n_m: 0.0,
        temperature_k: 300.0,
        rheology: ThermalCarreauYasuda::new(
            CarreauYasuda::new(mu0, mu_inf, time, index, 2.0, [0.0, ceiling]).unwrap(),
            TemperatureShift::new(
                TemperatureShiftModel::WilliamsLandelFerry {
                    c1: 1.0,
                    c2_k: 100.0,
                },
                300.0,
                [280.0, 360.0],
            )
            .unwrap(),
            ViscosityShiftMode::VerticalAndRate,
        ),
    }
}

fn limits() -> SlenderJetLimits {
    SlenderJetLimits {
        minimum_radius_m: 0.1,
        maximum_abs_slope: 0.75,
        minimum_radius_cells: 2.0,
    }
}
fn budget(dt: f64) -> SlenderJetStepBudget {
    SlenderJetStepBudget {
        maximum_step_s: dt,
        cfl: 0.2,
        maximum_substeps: 20_000,
    }
}
fn implicit(jet: &mut PeriodicSlenderJet, gate: &CancelGate) {
    jet.set_integration_method(SlenderJetIntegrationMethod::ImplicitViscosity, gate)
        .unwrap();
}
fn alternating(p: CarreauSlenderJetParameters, gate: &CancelGate) -> PeriodicSlenderJet {
    let velocity: Vec<_> = (0..16)
        .map(|i| if i % 2 == 0 { 0.125 } else { -0.125 })
        .collect();
    PeriodicSlenderJet::new_carreau(p, 4.0, &[1.0; 16], &velocity, limits(), gate).unwrap()
}

// Direct constitutive algebra at the reference temperature, independent of
// the original material evaluator and all implicit residual/Jacobian code.
fn mu(p: &CarreauSlenderJetParameters, strain: f64) -> f64 {
    let [mu0, mu_inf, time, index, transition] = p.rheology.curve().parameters();
    let rate = 3.0_f64.sqrt() * strain.abs();
    mu_inf
        + (mu0 - mu_inf) * (1.0 + (time * rate).powf(transition)).powf((index - 1.0) / transition)
}
fn manufactured_source(
    p: &CarreauSlenderJetParameters,
    a: &[f64],
    u: &[f64],
    dx: f64,
    dt: f64,
) -> Vec<f64> {
    let n = a.len();
    let stress: Vec<_> = (0..n)
        .map(|i| {
            let e = (u[i] - u[previous(i, n)]) / dx;
            3.0 * mu(p, e) / p.density_kg_m3 * a[i] * e
        })
        .collect();
    (0..n)
        .map(|i| {
            let r = next(i, n);
            average(a[i], a[r]) * u[i] - dt * (stress[r] - stress[i]) / dx
        })
        .collect()
}
fn norm(values: impl IntoIterator<Item = f64>) -> f64 {
    values.into_iter().fold(0.0_f64, f64::hypot)
}

#[test]
fn g1_original_stress_tangent_stage_energy_and_translation_agree() {
    let gate = CancelGate::new_clock_free();
    let n = 17;
    let (dx, dt) = (0.2, 0.03);
    let mut p = liquid([0.4, 0.003, 20.0, 0.2], 1000.0);
    p.density_kg_m3 = 1.7;
    let angle = |i: usize| core::f64::consts::TAU * i as f64 / n as f64;
    let a: Vec<_> = (0..n)
        .map(|i| 1.0 + 0.2 * fs_math::det::cos(angle(i)))
        .collect();
    let exact: Vec<_> = (0..n)
        .map(|i| 0.3 + 0.08 * fs_math::det::sin(angle(i)))
        .collect();
    let direction: Vec<_> = (0..n).map(|i| fs_math::det::cos(2.0 * angle(i))).collect();
    let source = manufactured_source(&p, &a, &exact, dx, dt);
    let mut work = ConstitutiveSolve::new(n, &gate).unwrap();
    let mut linear = Workspace::new(n, &gate).unwrap().solve;
    evaluate(
        &a,
        &exact,
        &source,
        0.0,
        &CarreauProfile::uniform(p),
        dx,
        dt,
        &mut work.flux,
        &mut work.tangent,
        &mut linear.rhs,
        &gate,
    )
    .unwrap();
    let applied: Vec<_> = (0..n)
        .map(|i| {
            let (l, r) = (previous(i, n), next(i, n));
            average(a[i], a[r]) * direction[i] + work.tangent[i] * (direction[i] - direction[l])
                - work.tangent[r] * (direction[r] - direction[i])
        })
        .collect();
    let h = 1e-6;
    let plus: Vec<_> = exact
        .iter()
        .zip(&direction)
        .map(|(u, d)| u + h * d)
        .collect();
    let minus: Vec<_> = exact
        .iter()
        .zip(&direction)
        .map(|(u, d)| u - h * d)
        .collect();
    let rp = manufactured_source(&p, &a, &plus, dx, dt);
    let rm = manufactured_source(&p, &a, &minus, dx, dt);
    let error = norm((0..n).map(|i| (rp[i] - rm[i]) / (2.0 * h) - applied[i]));
    assert!(error < 2e-8 * norm(applied.iter().copied()));
    let quadratic: f64 = direction.iter().zip(&applied).map(|(d, j)| d * j).sum();
    let positive_form: f64 = (0..n)
        .map(|i| {
            average(a[i], a[next(i, n)]) * direction[i].powi(2)
                + work.tangent[i] * (direction[i] - direction[previous(i, n)]).powi(2)
        })
        .sum();
    assert!(quadratic > 0.0 && (quadratic - positive_form).abs() < 2e-13 * quadratic);

    for boost in [0.0, 1000.0] {
        let shifted: Vec<_> = (0..n)
            .map(|i| source[i] + average(a[i], a[next(i, n)]) * boost)
            .collect();
        let mut state = Fields {
            a: a.clone(),
            q: shifted.clone(),
        };
        work.apply(
            &mut linear,
            &mut state,
            &CarreauProfile::uniform(p),
            dx,
            dt,
            SlenderJetImplicitControls::default(),
            &gate,
        )
        .unwrap();
        let actual: Vec<_> = (0..n)
            .map(|i| state.q[i] / average(a[i], a[next(i, n)]))
            .collect();
        assert!(
            actual
                .iter()
                .zip(&exact)
                .all(|(u, expected)| (u - boost - expected).abs() < 2e-10)
        );
        assert!((state.q.iter().sum::<f64>() - shifted.iter().sum::<f64>()).abs() < 1e-10);
        if boost == 0.0 {
            let factor = core::f64::consts::PI * p.density_kg_m3 * dx;
            let before: f64 = (0..n)
                .map(|i| 0.5 * factor * source[i].powi(2) / average(a[i], a[next(i, n)]))
                .sum();
            let after: f64 = (0..n).map(|i| 0.5 * factor * state.q[i] * actual[i]).sum();
            let damping: f64 = (0..n)
                .map(|i| {
                    let m = average(a[i], a[next(i, n)]);
                    0.5 * factor * m * (actual[i] - source[i] / m).powi(2)
                })
                .sum();
            let heat: f64 = (0..n)
                .map(|i| {
                    let e = (actual[i] - actual[previous(i, n)]) / dx;
                    dt * 3.0 * core::f64::consts::PI * dx * a[i] * mu(&p, e) * e * e
                })
                .sum();
            assert!(heat > 0.0 && damping > 0.0 && after < before);
            assert!((after - before + heat + damping).abs() < 3e-12 * before);
        }
    }
}

#[test]
fn g0_algebraic_trials_backtrack_inside_the_original_source_domain() {
    let gate = CancelGate::new_clock_free();
    let n = 16;
    let dx = 0.25;
    let angle = |i: usize| core::f64::consts::TAU * i as f64 / n as f64;
    let a: Vec<_> = (0..n)
        .map(|i| 1.0 + 0.65 * fs_math::det::sin(angle(i) + 0.52))
        .collect();
    let exact: Vec<_> = (0..n)
        .map(|i| {
            0.3 + 0.1 * fs_math::det::sin(angle(i)) + 0.04 * fs_math::det::sin(6.0 * angle(i) + 2.8)
        })
        .collect();
    let maximum = (0..n)
        .map(|i| 3.0_f64.sqrt() * (exact[i] - exact[previous(i, n)]).abs() / dx)
        .fold(0.0_f64, f64::max);
    let mut p = liquid([10.0, 0.001, 30.0, 0.02], 1.00001 * maximum);
    p.density_kg_m3 = 1.7;
    let mut state = Fields {
        a: a.clone(),
        q: manufactured_source(&p, &a, &exact, dx, 0.01),
    };
    let mut work = ConstitutiveSolve::new(n, &gate).unwrap();
    let mut linear = Workspace::new(n, &gate).unwrap().solve;
    let rejected = work
        .apply(
            &mut linear,
            &mut state,
            &CarreauProfile::uniform(p),
            dx,
            0.01,
            SlenderJetImplicitControls::default(),
            &gate,
        )
        .unwrap();
    assert!(
        rejected > 0,
        "the full Newton trial must leave the narrow original support"
    );
    for i in 0..n {
        let actual = state.q[i] / average(a[i], a[next(i, n)]);
        assert!((actual - exact[i]).abs() < 2e-10);
    }
}

#[test]
fn g3_all_newtonian_degeneracies_retain_the_native_implicit_trajectory() {
    let gate = CancelGate::new_clock_free();
    let velocity: Vec<_> = (0..16)
        .map(|i| 0.3 + 0.02 * fs_math::det::sin(core::f64::consts::TAU * i as f64 / 16.0))
        .collect();
    for values in [
        [0.4, 0.01, 20.0, 1.0],
        [0.4, 0.01, 0.0, 0.4],
        [0.4, 0.4, 20.0, 0.4],
    ] {
        for temperature in [300.0, 340.0] {
            let mut p = liquid(values, 1000.0);
            p.temperature_k = temperature;
            let mut actual =
                PeriodicSlenderJet::new_carreau(p, 4.0, &[1.0; 16], &velocity, limits(), &gate)
                    .unwrap();
            let mut native = PeriodicSlenderJet::new(
                actual.parameters(),
                4.0,
                &[1.0; 16],
                &velocity,
                limits(),
                &gate,
            )
            .unwrap();
            implicit(&mut actual, &gate);
            implicit(&mut native, &gate);
            assert_eq!(
                actual.advance(0.005, budget(0.001), &gate),
                native.advance(0.005, budget(0.001), &gate)
            );
            assert_eq!(actual.radius_squared_m2(), native.radius_squared_m2());
            assert_eq!(
                actual.radius_squared_flux_m3_s(),
                native.radius_squared_flux_m3_s()
            );
            assert_eq!(actual.carreau_parameters(), Some(&p));
        }
    }
}

// Alternating faces give q_i=(-1)^i*v, A_i=1-(-1)^i*2*I/dx.
// Their face masses remain exactly one and advective mass flux is zero.
// Thus I'=v, v'=-12*mu(2sqrt(3)|v|/dx)*v/(rho*dx²), with no frozen-area
// approximation. Independent scalar RK4 resolves this nonlinear reference.
fn scalar_reference(p: &CarreauSlenderJetParameters, duration: f64) -> (f64, f64, f64) {
    let dx = 0.25;
    let h = duration / 20_000.0;
    let f = |v| -12.0 * mu(p, 2.0 * v / dx) * v / p.density_kg_m3 / dx / dx;
    let (mut v, mut integral) = (0.125, 0.0);
    for _ in 0..20_000 {
        let k1 = f(v);
        let v2 = v + 0.5 * h * k1;
        let k2 = f(v2);
        let v3 = v + 0.5 * h * k2;
        let k3 = f(v3);
        let v4 = v + h * k3;
        let k4 = f(v4);
        integral += h * (v + 2.0 * v2 + 2.0 * v3 + v4) / 6.0;
        v += h * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
    }
    let heat = 0.5 * core::f64::consts::PI * p.density_kg_m3 * 4.0 * (0.125_f64.powi(2) - v * v);
    (v, integral, heat)
}

#[test]
fn g1_nonlinear_imex_time_and_physical_heat_refine_at_third_order() {
    let gate = CancelGate::new_clock_free();
    let p = liquid([0.8, 0.02, 10.0, 0.5], 1000.0);
    let (expected, integral, heat) = scalar_reference(&p, 0.02);
    let (mut velocity_errors, mut area_errors, mut heat_errors) =
        (Vec::new(), Vec::new(), Vec::new());
    for dt in [0.004, 0.002, 0.001] {
        let mut jet = alternating(p, &gate);
        implicit(&mut jet, &gate);
        let report = jet.advance(0.02, budget(dt), &gate).unwrap();
        assert_eq!(report.outcome, SlenderJetOutcome::Completed);
        velocity_errors.push((jet.velocities_m_s().next().unwrap() - expected).abs());
        area_errors.push((jet.radius_squared_m2()[0] - (1.0 - 8.0 * integral)).abs());
        heat_errors.push((report.viscous_dissipation_j - heat).abs());
        assert!(report.temporal_error.unwrap().maximum_velocity_error_m_s > 0.0);
        assert!(report.viscous_dissipation_j > 0.0 && report.upwind_dissipation_j < 1e-12);
        assert!((report.after.volume_m3 - report.before.volume_m3).abs() < 2e-12);
        assert!(report.after.axial_momentum_kg_m_s.abs() < 2e-12);
    }
    for errors in [&velocity_errors, &area_errors, &heat_errors] {
        assert!(
            errors[1] < 0.22 * errors[0] && errors[2] < 0.22 * errors[1],
            "{errors:?}"
        );
    }
    assert!(velocity_errors[2] < 1e-8 && heat_errors[2] < 1e-7);
}

#[test]
fn g1_stiff_strong_thinning_finishes_beyond_the_explicit_viscous_budget() {
    let gate = CancelGate::new_clock_free();
    let p = liquid([100.0, 0.0, 100.0, 0.02], 1000.0);
    let mut explicit = alternating(p, &gate);
    let mut actual = explicit.clone();
    let initial_mu = actual
        .carreau_sample_at(0)
        .unwrap()
        .unwrap()
        .dynamic_viscosity_pa_s();
    implicit(&mut actual, &gate);
    let short = SlenderJetStepBudget {
        maximum_substeps: 128,
        ..budget(0.002)
    };
    assert_eq!(
        explicit.advance(0.04, short, &gate),
        Err(NonlinearJetError::SubstepBudget)
    );
    assert_eq!(explicit.time_s(), 0.0);
    let solved = actual.advance(0.04, short, &gate).unwrap();
    let reference = explicit.advance(0.04, budget(0.002), &gate).unwrap();
    assert!(solved.accepted_substeps * 20 < reference.accepted_substeps);
    assert!(solved.after.kinetic_energy_j < 1e-12 * solved.before.kinetic_energy_j);
    assert!(solved.viscous_dissipation_j > 0.0);
    assert!(
        (solved.viscous_dissipation_j - solved.before.kinetic_energy_j).abs()
            < 0.005 * solved.before.kinetic_energy_j
    );
    for (a, b) in actual.radii_m().zip(explicit.radii_m()) {
        assert!((a - b).abs() < 2e-6);
    }
    assert!(
        actual
            .carreau_sample_at(0)
            .unwrap()
            .unwrap()
            .dynamic_viscosity_pa_s()
            > 20.0 * initial_mu
    );
    assert_eq!(actual.carreau_parameters(), Some(&p));
}

#[test]
fn g4_nonlinear_work_source_and_cancellation_refusals_preserve_the_full_call() {
    let gate = CancelGate::new_clock_free();
    let mut jet = alternating(liquid([100.0, 0.0, 100.0, 0.02], 1000.0), &gate);
    implicit(&mut jet, &gate);
    let mut fresh = jet.clone();
    let area = jet.radius_squared_m2().to_vec();
    let flux = jet.radius_squared_flux_m3_s().to_vec();
    let diagnostics = jet.diagnostics();
    let short = SlenderJetImplicitControls {
        maximum_iterations: 1,
        ..SlenderJetImplicitControls::default()
    };
    jet.set_implicit_viscosity_controls(short, &gate).unwrap();
    assert!(matches!(
        jet.advance(0.006, budget(0.002), &gate),
        Err(NonlinearJetError::ImplicitViscosity { iterations: 1, .. })
    ));
    assert_eq!(jet.radius_squared_m2(), area);
    assert_eq!(jet.radius_squared_flux_m3_s(), flux);
    assert_eq!(jet.diagnostics(), diagnostics);
    assert_eq!(jet.time_s(), 0.0);
    assert!(
        jet.set_implicit_viscosity_controls(
            SlenderJetImplicitControls {
                maximum_iterations: 0,
                ..short
            },
            &gate
        )
        .is_err()
    );
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert_eq!(
        jet.set_implicit_viscosity_controls(SlenderJetImplicitControls::default(), &cancelled),
        Err(NonlinearJetError::Cancelled)
    );
    assert_eq!(jet.implicit_viscosity_controls(), short);
    jet.set_implicit_viscosity_controls(SlenderJetImplicitControls::default(), &gate)
        .unwrap();
    assert_eq!(
        jet.advance(0.006, budget(0.002), &gate),
        fresh.advance(0.006, budget(0.002), &gate)
    );
    assert_eq!(
        jet.radius_squared_flux_m3_s(),
        fresh.radius_squared_flux_m3_s()
    );

    let mut p = liquid([0.4, 0.005, 20.0, 0.5], 1e-30);
    p.surface_tension_n_m = 1.0;
    let radii: Vec<_> = (0..32)
        .map(|i| 1.0 + 0.05 * fs_math::det::cos(core::f64::consts::TAU * i as f64 / 32.0))
        .collect();
    let mut restricted =
        PeriodicSlenderJet::new_carreau(p, 8.0, &radii, &[0.0; 32], limits(), &gate).unwrap();
    implicit(&mut restricted, &gate);
    let before = restricted.clone();
    let original = p.rheology.evaluate(300.0, 1.0).unwrap_err();
    assert!(matches!(restricted.advance(0.01, budget(0.001), &gate),
        Err(NonlinearJetError::Rheology { cell: Some(_), source }) if source == original));
    assert_eq!(restricted.radius_squared_m2(), before.radius_squared_m2());
    assert_eq!(
        restricted.radius_squared_flux_m3_s(),
        before.radius_squared_flux_m3_s()
    );
    assert_eq!(restricted.diagnostics(), before.diagnostics());
    assert_eq!(restricted.time_s(), 0.0);
    assert_eq!(restricted.carreau_parameters(), Some(&p));
    assert_eq!(
        restricted.advance(0.01, budget(0.001), &cancelled),
        Err(NonlinearJetError::Cancelled)
    );
    assert_eq!(
        restricted.advance(0.0, budget(0.001), &gate).unwrap().after,
        before.diagnostics()
    );
}
