use super::*;
use fs_material::fluid_transport::rheology::{
    CarreauYasuda, ThermalCarreauYasuda, ViscosityShiftMode,
};
use fs_material::visco::temperature::{TemperatureShift, TemperatureShiftModel};

fn material(values: [f64; 4]) -> CarreauSlenderJetParameters {
    let [mu0, mu_inf, time, index] = values;
    CarreauSlenderJetParameters {
        density_kg_m3: 2.0,
        surface_tension_n_m: 0.0,
        temperature_k: 300.0,
        rheology: ThermalCarreauYasuda::new(
            CarreauYasuda::new(mu0, mu_inf, time, index, 2.0, [0.0, 1000.0]).unwrap(),
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
        maximum_substeps: 10000,
    }
}
fn alternating(p: CarreauSlenderJetParameters, gate: &CancelGate) -> PeriodicSlenderJet {
    let u: Vec<_> = (0..16)
        .map(|i| if i % 2 == 0 { 0.125 } else { -0.125 })
        .collect();
    PeriodicSlenderJet::new_carreau(p, 4.0, &[1.0; 16], &u, limits(), gate).unwrap()
}
fn temperatures() -> Vec<f64> {
    (0..16)
        .map(|i| if i % 2 == 0 { 300.0 } else { 330.0 })
        .collect()
}
// Direct WLF/Carreau declaration, independent of preparation/source evaluators.
fn viscosity(p: CarreauSlenderJetParameters, temperature: f64, strain: f64) -> f64 {
    let [mu0, mu_inf, time, index, transition] = p.rheology.curve().parameters();
    let delta = temperature - 300.0;
    let shift = 10.0_f64.powf(-delta / (100.0 + delta));
    let rate = shift * 3.0_f64.sqrt() * strain.abs();
    shift
        * (mu_inf
            + (mu0 - mu_inf)
                * (1.0 + (time * rate).powf(transition)).powf((index - 1.0) / transition))
}
fn near(a: f64, b: f64, relative: f64) {
    assert!(
        (a - b).abs() <= relative * a.abs().max(b.abs()).max(1e-12),
        "{a} vs {b}"
    );
}

#[test]
fn g1_prescribed_cell_temperatures_resolve_original_stress_tangent_and_local_power() {
    let gate = CancelGate::new_clock_free();
    let p = material([0.8, 0.02, 10.0, 0.5]);
    let radius: Vec<_> = (0..16)
        .map(|i| 1.0 + 0.03 * fs_math::det::sin(core::f64::consts::TAU * i as f64 / 16.0))
        .collect();
    let velocity: Vec<_> = (0..16)
        .map(|i| 0.2 + 0.07 * fs_math::det::cos(core::f64::consts::TAU * i as f64 / 16.0))
        .collect();
    let t: Vec<_> = (0..16).map(|i| 290.0 + 3.0 * i as f64).collect();
    let mut jet =
        PeriodicSlenderJet::new_carreau(p, 4.0, &radius, &velocity, limits(), &gate).unwrap();
    jet.set_carreau_cell_temperatures(&t, 100.0, &gate).unwrap();
    jet.enable_cell_viscous_heat(&gate).unwrap();
    let mut full = Fields::zeros(16, &gate).unwrap();
    let mut inviscid = Fields::zeros(16, &gate).unwrap();
    let actual = evaluate_constitutive(
        &jet.accepted,
        jet.parameters,
        jet.carreau.as_ref(),
        jet.limits,
        jet.spacing_m,
        &mut jet.spatial,
        &mut full,
        true,
        &gate,
    )
    .unwrap();
    evaluate_constitutive(
        &jet.accepted,
        jet.parameters,
        jet.carreau.as_ref(),
        jet.limits,
        jet.spacing_m,
        &mut jet.spatial,
        &mut inviscid,
        false,
        &gate,
    )
    .unwrap();
    let u: Vec<_> = jet.velocities_m_s().collect();
    let dx = jet.spacing_m;
    let mut stress = [0.0; 16];
    let mut power = 0.0;
    for i in 0..16 {
        let e = (u[i] - u[previous(i, 16)]) / dx;
        let mu = viscosity(p, t[i], e);
        let sample = jet.carreau_sample_at(i).unwrap().unwrap();
        assert_eq!(sample.temperature().temperature_k(), t[i]);
        assert_eq!(sample.law(), &p.rheology);
        near(sample.dynamic_viscosity_pa_s(), mu, 2e-13);
        stress[i] = 3.0 * mu * jet.accepted.a[i] * e / p.density_kg_m3;
        let expected = 3.0 * core::f64::consts::PI * dx * mu * jet.accepted.a[i] * e * e;
        near(jet.spatial.heat.as_ref().unwrap().power[i], expected, 3e-13);
        power += expected;
        let delta = 1e-6;
        let derivative = (3.0 * viscosity(p, t[i], e + delta) * (e + delta)
            - 3.0 * viscosity(p, t[i], e - delta) * (e - delta))
            / (2.0 * delta);
        near(derivative, 3.0 * sample.stress_tangent_pa_s(), 2e-8);
    }
    for i in 0..16 {
        near(
            full.q[i] - inviscid.q[i],
            (stress[next(i, 16)] - stress[i]) / dx,
            2e-11,
        );
    }
    near(actual.viscous_power_w, power, 3e-13);
}

#[test]
fn g3_uniform_temperature_preparation_preserves_both_native_trajectories() {
    let gate = CancelGate::new_clock_free();
    for values in [
        [0.4, 0.01, 20.0, 1.0],
        [0.4, 0.01, 0.0, 0.4],
        [0.4, 0.4, 20.0, 0.4],
        [0.4, 0.01, 20.0, 0.4],
    ] {
        for method in [
            SlenderJetIntegrationMethod::ExplicitSspRk3,
            SlenderJetIntegrationMethod::ImplicitViscosity,
        ] {
            let mut p = material(values);
            p.temperature_k = 320.0;
            let mut original = alternating(p, &gate);
            let mut prepared = original.clone();
            original.set_integration_method(method, &gate).unwrap();
            prepared.set_integration_method(method, &gate).unwrap();
            prepared
                .set_carreau_cell_temperatures(&[320.0; 16], 100.0, &gate)
                .unwrap();
            prepared.enable_cell_viscous_heat(&gate).unwrap();
            let a = original.advance(0.01, budget(0.001), &gate).unwrap();
            let b = prepared.advance(0.01, budget(0.001), &gate).unwrap();
            assert_eq!(a, b);
            assert_eq!(original.radius_squared_m2(), prepared.radius_squared_m2());
            assert_eq!(
                original.radius_squared_flux_m3_s(),
                prepared.radius_squared_flux_m3_s()
            );
            let heat = prepared.cell_viscous_heat().unwrap();
            assert_eq!(heat.total_heat_j, b.viscous_dissipation_j);
            assert!(heat.cell_heat_j.iter().all(|x| *x >= 0.0));
            assert!(heat.partition_residual_j.abs() < 3e-13 * heat.total_heat_j);
        }
    }
}

// The full conservative alternating mode retains M=1, q_i=(-1)^i*v and
// A_i=1-(-1)^i*2*I/dx. Nonuniform temperatures change the two original cell
// viscosities. Their SUM weighted by the evolving areas governs v; local
// physical heat remains different in the two cell families.
fn reference(p: CarreauSlenderJetParameters, duration: f64) -> [f64; 4] {
    let dx = 0.25;
    let rhs = |y: [f64; 4]| {
        let [v, integral, _, _] = y;
        let area = [1.0 - 2.0 * integral / dx, 1.0 + 2.0 * integral / dx];
        let mu = [
            viscosity(p, 300.0, 2.0 * v / dx),
            viscosity(p, 330.0, 2.0 * v / dx),
        ];
        [
            -6.0 * v * (mu[0] * area[0] + mu[1] * area[1]) / p.density_kg_m3 / dx / dx,
            v,
            12.0 * core::f64::consts::PI * mu[0] * area[0] * v * v / dx,
            12.0 * core::f64::consts::PI * mu[1] * area[1] * v * v / dx,
        ]
    };
    let h = duration / 20000.0;
    let mut y = [0.125, 0.0, 0.0, 0.0];
    for _ in 0..20000 {
        let k1 = rhs(y);
        let k2 = rhs(std::array::from_fn(|i| y[i] + 0.5 * h * k1[i]));
        let k3 = rhs(std::array::from_fn(|i| y[i] + 0.5 * h * k2[i]));
        let k4 = rhs(std::array::from_fn(|i| y[i] + h * k3[i]));
        for i in 0..4 {
            y[i] += h * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0;
        }
    }
    y
}

#[test]
fn g1_nonuniform_thermal_viscosity_and_local_heat_refine_at_third_order() {
    let gate = CancelGate::new_clock_free();
    let p = material([0.8, 0.02, 10.0, 0.5]);
    let expected = reference(p, 0.02);
    let mut errors = Vec::new();
    for dt in [0.004, 0.002, 0.001] {
        let mut jet = alternating(p, &gate);
        jet.set_carreau_cell_temperatures(&temperatures(), 100.0, &gate)
            .unwrap();
        jet.set_integration_method(SlenderJetIntegrationMethod::ImplicitViscosity, &gate)
            .unwrap();
        jet.enable_cell_viscous_heat(&gate).unwrap();
        let report = jet.advance(0.02, budget(dt), &gate).unwrap();
        let heat = jet.cell_viscous_heat().unwrap();
        assert_eq!(report.outcome, SlenderJetOutcome::Completed);
        assert_eq!(heat.total_heat_j, report.viscous_dissipation_j);
        assert!(heat.partition_residual_j.abs() < 3e-13 * heat.total_heat_j);
        let local = heat.cell_heat_j;
        for i in 0..16 {
            near(local[i], local[i % 2], 2e-12);
        }
        assert!(local[0] > local[1]);
        errors.push([
            (jet.velocities_m_s().next().unwrap() - expected[0]).abs(),
            (jet.radius_squared_m2()[0] - (1.0 - 8.0 * expected[1])).abs(),
            (local[0] - expected[2]).abs(),
            (local[1] - expected[3]).abs(),
        ]);
        assert!((report.after.volume_m3 - report.before.volume_m3).abs() < 2e-12);
        assert!(report.after.axial_momentum_kg_m_s.abs() < 2e-12);
    }
    for component in 0..4 {
        assert!(
            errors[1][component] < 0.22 * errors[0][component]
                && errors[2][component] < 0.22 * errors[1][component],
            "{errors:?}"
        );
    }
}

#[test]
fn g3_variable_temperature_newtonian_limits_keep_variable_viscous_stress() {
    let gate = CancelGate::new_clock_free();
    for values in [
        [0.4, 0.01, 20.0, 1.0],
        [0.4, 0.01, 0.0, 0.4],
        [0.4, 0.4, 20.0, 0.4],
    ] {
        let p = material(values);
        let mut jet = alternating(p, &gate);
        let mut frozen = jet.clone();
        // The policy is selected while the law is spatially constant; binding
        // temperatures afterwards must prepare the variable-coefficient solve.
        jet.set_integration_method(SlenderJetIntegrationMethod::ImplicitViscosity, &gate)
            .unwrap();
        frozen
            .set_integration_method(SlenderJetIntegrationMethod::ImplicitViscosity, &gate)
            .unwrap();
        jet.set_carreau_cell_temperatures(&temperatures(), 100.0, &gate)
            .unwrap();
        jet.set_implicit_viscosity_controls(
            SlenderJetImplicitControls {
                maximum_iterations: 1,
                ..SlenderJetImplicitControls::default()
            },
            &gate,
        )
        .unwrap();
        jet.enable_cell_viscous_heat(&gate).unwrap();
        let expected = reference(p, 0.002);
        let report = jet.advance(0.002, budget(0.0001), &gate).unwrap();
        frozen.advance(0.002, budget(0.0001), &gate).unwrap();
        near(jet.velocities_m_s().next().unwrap(), expected[0], 2e-9);
        assert_ne!(
            jet.radius_squared_flux_m3_s(),
            frozen.radius_squared_flux_m3_s()
        );
        let heat = jet.cell_viscous_heat().unwrap();
        near(heat.cell_heat_j[0], expected[2], 2e-8);
        near(heat.cell_heat_j[1], expected[3], 2e-8);
        assert!(heat.partition_residual_j.abs() < 3e-13 * report.viscous_dissipation_j);
    }
}

#[test]
fn g4_temperature_source_and_flow_failures_preserve_accepted_heat_and_history() {
    let gate = CancelGate::new_clock_free();
    let p = material([100.0, 0.0, 100.0, 0.02]);
    let mut jet = alternating(p, &gate);
    jet.set_carreau_cell_temperatures(&temperatures(), 100.0, &gate)
        .unwrap();
    jet.set_integration_method(SlenderJetIntegrationMethod::ImplicitViscosity, &gate)
        .unwrap();
    jet.enable_cell_viscous_heat(&gate).unwrap();
    jet.advance(0.002, budget(0.001), &gate).unwrap();
    let before = jet.clone();
    let mut invalid = temperatures();
    invalid[7] = 400.0;
    assert!(matches!(
        jet.set_carreau_cell_temperatures(&invalid, 100.0, &gate),
        Err(NonlinearJetError::Rheology { cell: Some(7), .. })
    ));
    assert!(
        jet.set_carreau_cell_temperatures(&temperatures(), 1e-30, &gate)
            .is_err()
    );
    assert_eq!(
        jet.prepared_carreau_cells(),
        before.prepared_carreau_cells()
    );
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert_eq!(
        jet.set_carreau_cell_temperatures(&[300.0; 16], 100.0, &cancelled),
        Err(NonlinearJetError::Cancelled)
    );
    assert_eq!(
        jet.advance(0.001, budget(0.001), &cancelled),
        Err(NonlinearJetError::Cancelled)
    );
    let short = SlenderJetStepBudget {
        maximum_substeps: 1,
        ..budget(0.0001)
    };
    assert_eq!(
        jet.advance(0.01, short, &gate),
        Err(NonlinearJetError::SubstepBudget)
    );
    assert_eq!(jet.radius_squared_m2(), before.radius_squared_m2());
    assert_eq!(
        jet.radius_squared_flux_m3_s(),
        before.radius_squared_flux_m3_s()
    );
    assert_eq!(jet.diagnostics(), before.diagnostics());
    assert_eq!(jet.time_s(), before.time_s());
    assert_eq!(
        jet.cell_viscous_heat().unwrap().cell_heat_j,
        before.cell_viscous_heat().unwrap().cell_heat_j
    );
    assert_eq!(
        jet.cell_viscous_heat().unwrap().total_heat_j,
        before.cell_viscous_heat().unwrap().total_heat_j
    );
    let mut fresh = before;
    assert_eq!(
        jet.advance(0.001, budget(0.001), &gate),
        fresh.advance(0.001, budget(0.001), &gate)
    );
    assert_eq!(
        jet.cell_viscous_heat().unwrap().cell_heat_j,
        fresh.cell_viscous_heat().unwrap().cell_heat_j
    );
    jet.advance(0.0, budget(0.001), &gate).unwrap();
    assert!(
        jet.cell_viscous_heat()
            .unwrap()
            .cell_heat_j
            .iter()
            .all(|v| *v == 0.0)
    );
    // An actual stage exceeds the newly prescribed physical rate envelope.
    // Neither its local power nor an algebraic line-search attempt is heat.
    for method in [
        SlenderJetIntegrationMethod::ExplicitSspRk3,
        SlenderJetIntegrationMethod::ImplicitViscosity,
    ] {
        let mut capillary = material([0.4, 0.005, 20.0, 0.5]);
        capillary.surface_tension_n_m = 1.0;
        let radii: Vec<_> = (0..32)
            .map(|i| 1.0 + 0.05 * fs_math::det::cos(core::f64::consts::TAU * i as f64 / 32.0))
            .collect();
        let mut restricted =
            PeriodicSlenderJet::new_carreau(capillary, 8.0, &radii, &[0.0; 32], limits(), &gate)
                .unwrap();
        restricted
            .set_carreau_cell_temperatures(&[300.0; 32], 1e-30, &gate)
            .unwrap();
        restricted.set_integration_method(method, &gate).unwrap();
        restricted.enable_cell_viscous_heat(&gate).unwrap();
        let accepted = restricted.clone();
        assert!(matches!(
            restricted.advance(0.01, budget(0.001), &gate),
            Err(NonlinearJetError::Rheology { cell: Some(_), .. })
        ));
        assert_eq!(restricted.radius_squared_m2(), accepted.radius_squared_m2());
        assert_eq!(
            restricted.radius_squared_flux_m3_s(),
            accepted.radius_squared_flux_m3_s()
        );
        assert_eq!(restricted.time_s(), 0.0);
        assert_eq!(
            restricted.cell_viscous_heat().unwrap().cell_heat_j,
            accepted.cell_viscous_heat().unwrap().cell_heat_j
        );
    }
}

#[test]
fn g0_rest_only_original_source_and_newtonian_heat_account_remain_admitted() {
    let gate = CancelGate::new_clock_free();
    let mut p = material([0.4, 0.01, 20.0, 0.5]);
    p.rheology = ThermalCarreauYasuda::new(
        CarreauYasuda::new(0.4, 0.01, 20.0, 0.5, 2.0, [0.0, 0.0]).unwrap(),
        *p.rheology.temperature_shift(),
        p.rheology.mode(),
    );
    let mut jet =
        PeriodicSlenderJet::new_carreau(p, 4.0, &[1.0; 16], &[0.25; 16], limits(), &gate).unwrap();
    jet.enable_cell_viscous_heat(&gate).unwrap();
    for method in [
        SlenderJetIntegrationMethod::ExplicitSspRk3,
        SlenderJetIntegrationMethod::ImplicitViscosity,
    ] {
        jet.set_integration_method(method, &gate).unwrap();
        assert_eq!(
            jet.advance(0.01, budget(0.001), &gate)
                .unwrap()
                .viscous_dissipation_j,
            0.0
        );
        assert!(
            jet.cell_viscous_heat()
                .unwrap()
                .cell_heat_j
                .iter()
                .all(|v| *v == 0.0)
        );
    }
    let mut native = PeriodicSlenderJet::new(
        SlenderJetParameters {
            density_kg_m3: 2.0,
            surface_tension_n_m: 0.0,
            kinematic_viscosity_m2_s: 0.1,
        },
        4.0,
        &[1.0; 16],
        &(0..16)
            .map(|i| if i % 2 == 0 { 0.125 } else { -0.125 })
            .collect::<Vec<_>>(),
        limits(),
        &gate,
    )
    .unwrap();
    assert!(native.cell_viscous_heat().is_none());
    native.enable_cell_viscous_heat(&gate).unwrap();
    let report = native.advance(0.002, budget(0.001), &gate).unwrap();
    let heat = native.cell_viscous_heat().unwrap();
    assert_eq!(heat.total_heat_j, report.viscous_dissipation_j);
    assert!(heat.partition_residual_j.abs() < 3e-13 * heat.total_heat_j);
}
