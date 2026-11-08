//! Actual source/momentum coupling against continuous thermoviscous evolution.
use super::*;
use fs_lbm::conjugate::enthalpy::advection::feedback::implicit::{
    AdaptiveImplicitFlowControls, AdaptiveImplicitFlowError,
};

fn adaptive(dt: f64, tolerance: f64) -> AdaptiveImplicitFlowControls {
    AdaptiveImplicitFlowControls {
        step: implicit_controls(),
        initial_dt_s: dt,
        minimum_dt_s: 1e-8,
        maximum_dt_s: dt,
        maximum_trials: 1000,
        enthalpy_absolute_tolerance_j_kg: tolerance,
        temperature_absolute_tolerance_k: None,
        relative_tolerance: 0.0,
        energy_tolerance_j: 2e-5,
    }
}

#[test]
fn adaptive_source_viscosity_converges_to_continuous_work_driven_heating() {
    let domain = VoxelDomain::new(2, 2, 1, 1.0).unwrap();
    let config = closed();
    let gate = CancelGate::new();
    let ids: Vec<_> = (0..4).map(|c| Some(CellId(c))).collect();
    let card = source(100.0, true);
    let liquid = liquid(&card);
    let caloric = liquid_caloric(&card);
    let old = initial(&domain, 1.0, &caloric);
    let flow = input(&ids, &liquid, &caloric)
        .feedback(&old, &domain, &config, (), &gate, 32768)
        .unwrap()
        .with_body_force_heating(force(), controls())
        .unwrap();
    // The original source curve has mu=100-9(T-300) on [300,310] K.
    // Native 2x2 circulation gives Tdot=A/mu, A=50/3. Integrating the
    // continuous equation gives 100*dT-4.5*dT²=A*t, not a BE recurrence.
    let duration = 2.0;
    let a = 50.0 / 3.0;
    let exact_h =
        100.0 * (2.0 * a * duration) / (100.0 + (10000.0_f64 - 18.0 * a * duration).sqrt());
    let mut errors = Vec::new();
    for tolerance in [5e-3, 5e-5] {
        let (next, report) = flow
            .advance_implicit_adaptive_to(&old, duration, adaptive(duration, tolerance), || Ok(()))
            .unwrap();
        let h = next.cells()[&CellId(0)].specific_h_j_kg;
        errors.push((h - exact_h).abs());
        for sample in &report.last_step.endpoint.samples {
            near(sample.specific_enthalpy_j_kg, h);
            near(next.cells()[&sample.cell].specific_h_j_kg, h);
            near(sample.temperature_k, 300.0 + h / 100.0);
            near(
                sample.material.dynamic_viscosity_pa_s(),
                100.0 - 9.0 * h / 100.0,
            );
        }
        near(report.thermal.stored_energy_change_j, 4.0 * h);
        near(
            report.deposited_heat_j,
            report.mechanical.supplied_body_work_j,
        );
        near(report.deposited_heat_j, report.thermal.external_heat_j);
        near(report.net_stream_enthalpy_j, 0.0);
        near(report.other_thermal_input_j, 0.0);
        assert!(report.combined_residual_j.abs() <= 2e-5);
        assert!(report.absolute_accepted_energy_residual_j <= 2e-5);
        assert!(report.accuracy_rejections > 0);
        assert!(report.maximum_accepted_normalized_error <= 1.0);
        assert_eq!(report.mechanical.substeps, 2 * report.accepted_intervals);
        assert_eq!(report.thermal.substeps, report.mechanical.substeps);
        assert!(report.momentum_solves > report.last_step.momentum_solves);
        assert_eq!(report.last_step.endpoint.time_s, duration);
    }
    assert!(
        errors[1] < errors[0] / 5.0,
        "continuous time errors {errors:?}"
    );
    let mut floor = adaptive(duration, 1e-5);
    floor.minimum_dt_s = duration;
    assert!(matches!(
        flow.advance_implicit_adaptive_to(&old, duration, floor, || Ok(())),
        Err(AdaptiveImplicitFlowError::TimeAccuracy { .. })
    ));
    let mut bounded = adaptive(duration, 1e-5);
    bounded.maximum_trials = 1;
    assert!(matches!(
        flow.advance_implicit_adaptive_to(&old, duration, bounded, || Ok(())),
        Err(AdaptiveImplicitFlowError::Solve(ImplicitFlowError::Flow(
            MaterialFlowError::Thermal(DepositionError::Budget("adaptive implicit-flow trials"))
        )))
    ));
    assert!(flow.applied().is_none());
    assert_eq!(flow.momentum_solves(), 0);
    assert_eq!(old.time_s(), 0.0);
}

#[test]
fn adaptive_mushy_heating_resolves_latent_storage_independently_of_h_datum() {
    let domain = VoxelDomain::new(2, 2, 1, 1.0).unwrap();
    let config = closed();
    let gate = CancelGate::new();
    let ids: Vec<_> = (0..4).map(|c| Some(CellId(c))).collect();
    let card = source(100.0, false);
    let liquid = liquid(&card);
    let duration = 0.2;
    let primitive = |phi: f64| -0.5 / (phi * phi) + 2.0 / phi + phi.ln();
    // Exact separation of dh/dt=10000/[600+mu/K(phi)] with permeability .5.
    let (mut lo, mut hi) = (50.0, 60.0);
    for _ in 0..64 {
        let h = 0.5 * (lo + hi);
        let elapsed = 0.06 * (h - 50.0) + 2.0 * (primitive(h / 100.0) - primitive(0.5));
        if elapsed < duration {
            lo = h;
        } else {
            hi = h;
        }
    }
    let exact = 0.5 * (lo + hi);
    let mut results = Vec::new();
    for datum in [0.0, 1_048_576.0] {
        let curve = EquilibriumEnthalpyPhaseCurve::try_new(
            card.content_hash(),
            [(0.0, 0.0), (100.0, 1.0)]
                .map(|(h, f)| EnthalpyPhaseKnot {
                    specific_enthalpy_j_kg: datum + h,
                    temperature_k: 500.0,
                    liquid_mass_fraction: f,
                    bulk_density_kg_m3: 1.0,
                })
                .to_vec(),
        )
        .unwrap();
        let old = thermal_state(&domain, &curve, &[datum + 50.0; 4], &[]);
        let mobility = MushyPhaseLaw::try_new(&curve, 1.0, 1.0, 0.5, 0.2).unwrap();
        let flow = MaterialFlowFeedback::try_new(
            &old,
            &domain,
            FlowEnthalpyInput {
                cell_ids: &ids,
                material_id: MAT,
                caloric: &curve,
                inlet_specific_h_j_kg: [Some(datum + 50.0); 6],
                max_cells: 4,
                max_connections: 64,
            },
            &config,
            (),
            &gate,
            32768,
            |_, _: &DepositedCell, t| {
                liquid.at_temperature(Temperature::new(t)).map_err(|e| {
                    fs_material::MaterialError::State {
                        what: e.to_string(),
                    }
                })
            },
        )
        .unwrap()
        .with_body_force_heating(force(), controls())
        .unwrap()
        .with_mushy_mobility(&mobility)
        .unwrap();
        let mut controls = adaptive(duration, 1e-5);
        controls.temperature_absolute_tolerance_k = Some(1e-10);
        let (next, report) = flow
            .advance_implicit_adaptive_to(&old, duration, controls, || Ok(()))
            .unwrap();
        let h = next.cells()[&CellId(0)].specific_h_j_kg - datum;
        assert!(
            (h - exact).abs() < 0.002,
            "latent h {h} != continuous {exact}"
        );
        assert!(report.accuracy_rejections > 0);
        assert_eq!(
            report.last_error.maximum_temperature_difference_k,
            Some(0.0)
        );
        assert!(report.maximum_accepted_normalized_error <= 1.0);
        assert!(report.mechanical.darcy_heat_j > 0.0);
        assert!(report.mechanical.molecular_heat_j > 0.0);
        near(report.deposited_heat_j, report.mechanical.generated_heat_j);
        near(report.thermal.stored_energy_change_j, 4.0 * (h - 50.0));
        for sample in &report.last_step.endpoint.samples {
            near(sample.temperature_k, 500.0);
        }
        results.push(h);
        assert_eq!(flow.momentum_solves(), 0);
    }
    assert!(
        (results[0] - results[1]).abs() < 2e-6,
        "datum response {results:?}"
    );
}

#[test]
fn all_comparison_solves_share_one_budget_and_native_step_counts_stay_local() {
    let domain = VoxelDomain::new(2, 2, 1, 1.0).unwrap();
    let config = closed();
    let gate = CancelGate::new();
    let ids: Vec<_> = (0..4).map(|c| Some(CellId(c))).collect();
    let card = source(100.0, false);
    let liquid = liquid(&card);
    let caloric = liquid_caloric(&card);
    let old = initial(&domain, 1.0, &caloric);
    let make = |budget| {
        input(&ids, &liquid, &caloric)
            .feedback(&old, &domain, &config, (), &gate, budget)
            .unwrap()
            .with_body_force_heating(force(), controls())
            .unwrap()
    };
    let flow = make(4096);
    let (_, coarse) = flow
        .step_implicit(&old, 0.5, implicit_controls(), &[0.0; 4], || Ok(()))
        .unwrap();
    let (half, first) = flow
        .step_implicit(&old, 0.25, implicit_controls(), &[0.0; 4], || Ok(()))
        .unwrap();
    let (_, second) = flow
        .step_implicit(&half, 0.5, implicit_controls(), &[0.0; 4], || Ok(()))
        .unwrap();
    let total = coarse.momentum_solves + first.momentum_solves + second.momentum_solves;
    let (_, report) = flow
        .advance_implicit_adaptive_to(&old, 0.5, adaptive(0.5, 1e-4), || Ok(()))
        .unwrap();
    assert_eq!(report.trials, 1);
    assert_eq!(report.momentum_solves, total);
    assert_eq!(report.last_step.momentum_solves, second.momentum_solves);
    near(report.mechanical.duration_s, 0.5);
    near(report.deposited_heat_j, 100.0 / 3.0);
    assert!(matches!(
        make(total - 1).advance_implicit_adaptive_to(&old, 0.5, adaptive(0.5, 1e-4), || Ok(())),
        Err(AdaptiveImplicitFlowError::Solve(ImplicitFlowError::Flow(
            MaterialFlowError::Thermal(DepositionError::Budget("live material momentum solves"))
        )))
    ));
    let mut polls = 0;
    assert!(matches!(
        flow.advance_implicit_adaptive_to(&old, 0.5, adaptive(0.5, 1e-4), || {
            polls += 1;
            if polls == 20 {
                gate.request();
            }
            Ok(())
        }),
        Err(AdaptiveImplicitFlowError::Solve(ImplicitFlowError::Flow(
            MaterialFlowError::Thermal(DepositionError::Cancelled)
        )))
    ));
    assert!(flow.applied().is_none());
    assert_eq!(flow.momentum_solves(), 0);
    assert_eq!(old.time_s(), 0.0);
}
