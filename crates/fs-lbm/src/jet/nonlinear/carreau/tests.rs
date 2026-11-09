use super::*;
use fs_material::fluid_transport::rheology::{CarreauYasuda, ViscosityShiftMode};
use fs_material::visco::temperature::{TemperatureShift, TemperatureShiftModel};

fn law(parameters: [f64; 5], support: [f64; 2]) -> ThermalCarreauYasuda {
    let [mu0, mu_inf, time, n, a] = parameters;
    ThermalCarreauYasuda::new(
        CarreauYasuda::new(mu0, mu_inf, time, n, a, support).unwrap(),
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
    )
}

fn parameters(rheology: ThermalCarreauYasuda) -> CarreauSlenderJetParameters {
    CarreauSlenderJetParameters {
        density_kg_m3: 2.0,
        surface_tension_n_m: 0.0,
        temperature_k: 300.0,
        rheology,
    }
}

fn limits() -> SlenderJetLimits {
    SlenderJetLimits {
        minimum_radius_m: 0.1,
        maximum_abs_slope: 0.75,
        minimum_radius_cells: 2.0,
    }
}

fn budget(maximum_step_s: f64) -> SlenderJetStepBudget {
    SlenderJetStepBudget {
        maximum_step_s,
        cfl: 0.2,
        maximum_substeps: 20_000,
    }
}

fn near(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs().max(1e-30),
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn g1_newtonian_limits_reproduce_the_original_trajectory_and_shifted_viscosity() {
    let gate = CancelGate::new_clock_free();
    let radius: Vec<_> = (0..32)
        .map(|i| 1.0 + 0.03 * fs_math::det::cos(core::f64::consts::TAU * i as f64 / 32.0))
        .collect();
    let velocity: Vec<_> = (0..32)
        .map(|i| 0.2 + 0.04 * fs_math::det::sin(core::f64::consts::TAU * (i as f64 + 0.5) / 32.0))
        .collect();
    for curve in [
        [0.4, 0.01, 20.0, 1.0, 2.0],
        [0.4, 0.01, 0.0, 0.5, 2.0],
        [0.4, 0.4, 20.0, 0.5, 0.7],
    ] {
        for temperature_k in [300.0, 340.0] {
            let mut p = parameters(law(curve, [0.0, 1000.0]));
            p.temperature_k = temperature_k;
            p.surface_tension_n_m = 1.0;
            let mut actual =
                PeriodicSlenderJet::new_carreau(p, 8.0, &radius, &velocity, limits(), &gate)
                    .unwrap();
            let shift = p
                .rheology
                .temperature_shift()
                .evaluate(temperature_k)
                .unwrap()
                .factor();
            near(
                actual.parameters().kinematic_viscosity_m2_s,
                0.4 * shift / 2.0,
                1e-15,
            );
            let mut reference = PeriodicSlenderJet::new(
                actual.parameters(),
                8.0,
                &radius,
                &velocity,
                limits(),
                &gate,
            )
            .unwrap();
            let report = actual.advance(0.01, budget(0.001), &gate).unwrap();
            let native = reference.advance(0.01, budget(0.001), &gate).unwrap();
            assert_eq!(report, native);
            assert_eq!(actual.radius_squared_m2(), reference.radius_squared_m2());
            assert_eq!(
                actual.radius_squared_flux_m3_s(),
                reference.radius_squared_flux_m3_s()
            );
            assert_eq!(actual.carreau_parameters(), Some(&p));
            assert!(reference.carreau_sample_at(0).unwrap().is_none());
        }
    }
}

#[test]
fn g1_extensional_tensor_stress_heat_and_cfl_match_independent_carreau_algebra() {
    let gate = CancelGate::new_clock_free();
    let p = parameters(law([0.4, 0.005, 20.0, 0.5, 2.0], [0.0, 1000.0]));
    let n = 16;
    let velocity: Vec<_> = (0..n)
        .map(|i| if i % 2 == 0 { 0.125 } else { -0.125 })
        .collect();
    let mut jet =
        PeriodicSlenderJet::new_carreau(p, 4.0, &[1.0; 16], &velocity, limits(), &gate).unwrap();
    let actual = evaluate_constitutive(
        &jet.accepted,
        jet.parameters,
        jet.carreau.as_ref(),
        jet.limits,
        jet.spacing_m,
        &mut jet.spatial,
        &mut jet.rhs,
        true,
        &gate,
    )
    .unwrap();
    let rhs = jet.rhs.q.clone();
    let mut inviscid_rhs = Fields::zeros(n, &gate).unwrap();
    let inviscid = evaluate(
        &jet.accepted,
        SlenderJetParameters {
            kinematic_viscosity_m2_s: 0.0,
            ..jet.parameters
        },
        jet.limits,
        jet.spacing_m,
        &mut jet.spatial,
        &mut inviscid_rhs,
        &gate,
    )
    .unwrap();
    let dx = jet.spacing_m;
    let mut stresses = vec![0.0; n];
    let mut heat = 0.0;
    let mut largest_mu = 0.0_f64;
    let mut largest_tangent = 0.0_f64;
    for i in 0..n {
        let strain = (velocity[i] - velocity[previous(i, n)]) / dx;
        // Incompressibility fixes two transverse rates, not another shear fit.
        let d = [-0.5 * strain, -0.5 * strain, strain];
        let rate = (2.0 * d.iter().map(|x| x * x).sum::<f64>()).sqrt();
        let x2 = (20.0 * rate) * (20.0 * rate);
        let excess = (0.4 - 0.005) / (1.0 + x2).sqrt().sqrt();
        let mu = 0.005 + excess;
        let tangent = 0.005 + excess * (1.0 + 0.5 * x2) / (1.0 + x2);
        let sample = jet.carreau_sample_at(i).unwrap().unwrap();
        near(sample.shear_rate_per_s(), rate, 2e-15);
        near(sample.dynamic_viscosity_pa_s(), mu, 3e-13);
        near(sample.stress_tangent_pa_s(), tangent, 3e-13);
        stresses[i] = 3.0 * mu / p.density_kg_m3 * strain;
        heat += 3.0 * core::f64::consts::PI * mu * strain * strain * dx;
        largest_mu = largest_mu.max(mu);
        largest_tangent = largest_tangent.max(tangent);
    }
    for i in 0..n {
        near(
            rhs[i] - inviscid_rhs.q[i],
            (stresses[next(i, n)] - stresses[i]) / dx,
            4e-13,
        );
    }
    near(actual.viscous_power_w, heat, 4e-13);
    near(
        actual.rate_per_s - inviscid.rate_per_s,
        12.0 * largest_mu / p.density_kg_m3 / dx / dx,
        4e-13,
    );
    // A tangent-only bound would miss the stronger apparent conductance.
    assert!(largest_mu > 1.5 * largest_tangent);
    assert!(rhs.iter().sum::<f64>().abs() < 1e-12);
}

#[test]
fn g1_g3_strong_thinning_remains_dissipative_and_changes_the_real_evolution() {
    let gate = CancelGate::new_clock_free();
    let p = parameters(law([100.0, 0.0, 100.0, 0.02, 2.0], [0.0, 1000.0]));
    let velocity: Vec<_> = (0..16)
        .map(|i| if i % 2 == 0 { 0.125 } else { -0.125 })
        .collect();
    let mut actual =
        PeriodicSlenderJet::new_carreau(p, 4.0, &[1.0; 16], &velocity, limits(), &gate).unwrap();
    let initial_mu = actual
        .carreau_sample_at(0)
        .unwrap()
        .unwrap()
        .dynamic_viscosity_pa_s();
    let mut frozen = PeriodicSlenderJet::new(
        actual.parameters(),
        4.0,
        &[1.0; 16],
        &velocity,
        limits(),
        &gate,
    )
    .unwrap();
    let report = actual.advance(0.005, budget(1.0), &gate).unwrap();
    let newtonian = frozen.advance(0.005, budget(1.0), &gate).unwrap();
    assert_eq!(report.outcome, SlenderJetOutcome::Completed);
    assert!(report.viscous_dissipation_j > 0.0);
    assert!(report.after.kinetic_energy_j < report.before.kinetic_energy_j);
    assert!(report.viscous_dissipation_j < 0.6 * newtonian.viscous_dissipation_j);
    assert_ne!(
        actual.radius_squared_flux_m3_s(),
        frozen.radius_squared_flux_m3_s()
    );
    near(report.after.volume_m3, report.before.volume_m3, 3e-12);
    assert!(report.after.axial_momentum_kg_m_s.abs() < 1e-12);
    assert!(
        actual
            .carreau_sample_at(0)
            .unwrap()
            .unwrap()
            .dynamic_viscosity_pa_s()
            > initial_mu
    );
}

#[test]
fn g0_g4_source_failure_and_policy_selection_preserve_accepted_state() {
    let gate = CancelGate::new_clock_free();
    let mut p = parameters(law([0.4, 0.005, 20.0, 0.5, 2.0], [0.0, 1e-12]));
    p.surface_tension_n_m = 1.0;
    let radius: Vec<_> = (0..32)
        .map(|i| 1.0 + 0.05 * fs_math::det::cos(core::f64::consts::TAU * i as f64 / 32.0))
        .collect();
    let mut jet =
        PeriodicSlenderJet::new_carreau(p, 8.0, &radius, &[0.0; 32], limits(), &gate).unwrap();
    let area = jet.radius_squared_m2().to_vec();
    let flux = jet.radius_squared_flux_m3_s().to_vec();
    let diagnostics = jet.diagnostics();
    jet.set_integration_method(SlenderJetIntegrationMethod::ImplicitViscosity, &gate)
        .unwrap();
    assert_eq!(
        jet.integration_method(),
        SlenderJetIntegrationMethod::ImplicitViscosity
    );
    jet.set_integration_method(SlenderJetIntegrationMethod::ExplicitSspRk3, &gate)
        .unwrap();
    assert!(matches!(
        jet.advance(0.01, budget(0.001), &gate),
        Err(NonlinearJetError::Rheology { cell: Some(_), .. })
    ));
    assert_eq!(jet.radius_squared_m2(), area);
    assert_eq!(jet.radius_squared_flux_m3_s(), flux);
    assert_eq!(jet.diagnostics(), diagnostics);
    assert_eq!(jet.time_s(), 0.0);
    assert_eq!(jet.carreau_parameters(), Some(&p));
    assert_eq!(
        jet.advance(0.0, budget(0.001), &gate).unwrap().after,
        diagnostics
    );
    assert!(jet.carreau_sample_at(32).is_err());

    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert_eq!(
        jet.advance(0.01, budget(0.001), &cancelled),
        Err(NonlinearJetError::Cancelled)
    );
    p.temperature_k = 400.0;
    assert!(matches!(
        PeriodicSlenderJet::new_carreau(p, 8.0, &radius, &[0.0; 32], limits(), &gate),
        Err(NonlinearJetError::Rheology { cell: None, .. })
    ));
}
