//! G1/G3: original reflected heat, physical reservoir entropy and native trajectories.
use super::open::{reference, spatial};
use super::*;
use fs_conduction::transient::enthalpy::adaptive::{
    AdaptiveEnthalpyConfig, AdaptiveEnthalpyStatus,
};
use fs_conduction::transient::enthalpy::heterogeneous::EnthalpyRadiation;
use fs_conduction::transient::enthalpy::sdirk::{EnthalpySdirk2Input, SDIRK2_GAMMA};

#[test]
fn closed_enclosure_entropy_preserves_native_be_and_adaptive_heat_exchange() {
    let fixture = Fixture::new(false);
    let boundary = boundary(&fixture, false);
    let source = ScalarField::Uniform(0.0);
    let old = [
        2000.0, 2000.0, 2000.0, 2000.0, 1000.0, 1000.0, 1000.0, 1000.0,
    ];
    let capacity = RHO * CP * 2.0_f64.sqrt() / 12.0;
    let area = 3.0_f64.sqrt();
    let exchange = area * STEFAN_BOLTZMANN_W_M2_K4 / (1.0 / EPS[0] + 1.0 / EPS[1] - 1.0);
    // Independent two-isothermal-body BE balance. The total temperature is
    // exactly 700 K and radiation is internal to the closed accounting domain.
    let (mut lo, mut hi) = (0.0, 100.0);
    for _ in 0..80 {
        let difference = 0.5 * (lo + hi);
        let hot = 350.0_f64 + 0.5 * difference;
        let cold = 350.0_f64 - 0.5 * difference;
        let residual =
            capacity * (difference - 100.0) + 2.0 * DT * exchange * (hot.powi(4) - cold.powi(4));
        if residual > 0.0 {
            hi = difference;
        } else {
            lo = difference;
        }
    }
    let expected = [350.0 + 0.25 * (lo + hi), 350.0 - 0.25 * (lo + hi)];
    let heat = exchange * (expected[0].powi(4) - expected[1].powi(4));
    with_cx(|cx| {
        let radiation = PointwiseEnclosureRadiation::new(
            cx,
            &fixture.mesh,
            &boundary,
            &fixture.enclosure,
            PointwiseEnclosureTolerance::default(),
        )
        .unwrap();
        let p = problem(&fixture, &boundary, &source);
        let engine = fixture.engine(cx, false);
        let native = engine
            .advance_with_pointwise_enclosure_radiation(
                cx,
                p,
                None,
                &old,
                DT,
                step_config(),
                &radiation,
            )
            .unwrap();
        let (step, entropy) = engine
            .advance_with_radiation_and_entropy(
                cx,
                p,
                None,
                &old,
                DT,
                step_config(),
                EnthalpyRadiation::PointwiseEnclosure(&radiation),
                1e-9,
            )
            .unwrap();
        assert_eq!(step.specific_enthalpy_j_kg, native.specific_enthalpy_j_kg);
        assert_eq!(
            step.pointwise_enclosure_radiation,
            native.pointwise_enclosure_radiation
        );
        for (body, temperature) in expected.into_iter().enumerate() {
            for node in 0..4 {
                close(step.temperature[4 * body + node], temperature, 2e-10);
            }
        }
        let report = step.pointwise_enclosure_radiation.as_ref().unwrap();
        assert!(report.ambient.is_none());
        close(report.outward_heat_w, 0.0, 1e-9);
        close(step.stored_energy_change_j, 0.0, 1e-8);
        assert_eq!(entropy.radiation_entropy_input_j_k, 0.0);
        close(
            entropy.radiation_generation_j_k,
            DT * heat * (1.0 / expected[1] - 1.0 / expected[0]),
            2e-10,
        );
        close(
            entropy.thermal_storage_change_j_k,
            capacity * ((expected[0] / 400.0).ln() + (expected[1] / 300.0).ln()),
            2e-10,
        );
        assert!(entropy.radiation_generation_j_k > 0.0);
        assert!(entropy.temporal_defect_j_k > 0.0);
        close(entropy.balance_residual_j_k, 0.0, 1e-9);
        assert!(!entropy.has_unclassified_inputs);

        let policy = AdaptiveEnthalpyConfig {
            enthalpy_absolute_tolerance_j_kg: 100.0,
            temperature_absolute_tolerance_k: 10.0,
            liquid_fraction_absolute_tolerance: 1e-6,
            relative_tolerance: 0.0,
            minimum_step_s: 1e-6,
            maximum_step_s: DT,
            endpoint: step_config(),
        };
        let mut ordinary = engine.adaptive_state(cx, &old, 0.0, DT, policy).unwrap();
        let mut audited = engine
            .adaptive_state_with_entropy(cx, &old, 0.0, DT, policy, 1e-8)
            .unwrap();
        let native_history = ordinary
            .advance_to_with_pointwise_enclosure_radiation(cx, p, None, &radiation, DT, 10)
            .unwrap();
        let history = audited
            .advance_to_with_radiation(
                cx,
                p,
                None,
                EnthalpyRadiation::PointwiseEnclosure(&radiation),
                DT,
                10,
            )
            .unwrap();
        assert_eq!(history.thermal, native_history);
        assert_eq!(history.thermal.status, AdaptiveEnthalpyStatus::ReachedEnd);
        assert_eq!(history.thermal.accepted, 1);
        assert_eq!(
            audited.history().specific_enthalpy_j_kg(),
            ordinary.specific_enthalpy_j_kg()
        );
        assert_eq!(history.entropy.radiation_entropy_input_j_k, 0.0);
        let state_entropy: f64 = audited
            .history()
            .temperature()
            .iter()
            .enumerate()
            .map(|(i, &t)| capacity / 4.0 * (t / if i < 4 { 400.0 } else { 300.0 }).ln())
            .sum();
        close(
            history.entropy.thermal_storage_change_j_k,
            state_entropy,
            1e-9,
        );
        close(history.entropy.balance_residual_j_k, 0.0, 1e-8);
    });
    assert_eq!(
        old,
        [
            2000.0, 2000.0, 2000.0, 2000.0, 1000.0, 1000.0, 1000.0, 1000.0
        ]
    );
}

#[test]
fn open_enclosure_sdirk_retains_each_actual_stage_bath_and_reflected_nodal_heat() {
    let fixture = Fixture::new(false);
    let boundary = boundary(&fixture, false);
    let source = ScalarField::Uniform(0.0);
    let initial = [380.0, 380.0, 380.0, 380.0, 330.0, 330.0, 330.0, 330.0];
    let old = initial.map(|t| CP * (t - 200.0));
    with_cx(|cx| {
        let first = spatial(&fixture, &boundary, cx, EPS, 280.0);
        let second = first.with_ambient_temperature_k(500.0).unwrap();
        let p = problem(&fixture, &boundary, &source);
        let inputs = [&first, &second].map(|radiation| EnthalpySdirk2Input {
            problem: p,
            interfaces: None,
            radiation: EnthalpyRadiation::PointwiseEnclosure(radiation),
            contacts: None,
        });
        let engine = fixture.engine(cx, false);
        let native = engine
            .advance_sdirk2(cx, inputs, &old, DT, step_config())
            .unwrap();
        let (step, entropy) = engine
            .advance_sdirk2_with_entropy(cx, inputs, &old, DT, step_config(), 1e-9)
            .unwrap();
        assert_eq!(
            step.first_stage.specific_enthalpy_j_kg,
            native.first_stage.specific_enthalpy_j_kg
        );
        assert_eq!(
            step.second_stage.specific_enthalpy_j_kg,
            native.second_stage.specific_enthalpy_j_kg
        );
        // Original independent simplex moments and explicit two-patch radiosity
        // equations: the bath has its own reciprocal row and no solid node.
        let a = reference(&fixture, &step.first_stage.temperature, EPS, 280.0);
        let b = reference(&fixture, &step.second_stage.temperature, EPS, 500.0);
        assert!(a.2 < 0.0 && b.2 > 0.0);
        let surface_rate =
            |load: &[f64], t: &[f64]| -> f64 { load.iter().zip(t).map(|(q, t)| q / t).sum() };
        let input = (1.0 - SDIRK2_GAMMA) * a.2 / 280.0 + SDIRK2_GAMMA * b.2 / 500.0;
        let generation = -(1.0 - SDIRK2_GAMMA)
            * (surface_rate(&a.0, &step.first_stage.temperature) + a.2 / 280.0)
            - SDIRK2_GAMMA * (surface_rate(&b.0, &step.second_stage.temperature) + b.2 / 500.0);
        close(entropy.radiation_entropy_input_j_k, DT * input, 2e-10);
        close(entropy.radiation_generation_j_k, DT * generation, 2e-10);
        let nodal_capacity = RHO * CP * 2.0_f64.sqrt() / 48.0;
        let state_entropy: f64 = step
            .second_stage
            .temperature
            .iter()
            .zip(initial)
            .map(|(&t, old_t)| nodal_capacity * (t / old_t).ln())
            .sum();
        close(entropy.thermal_storage_change_j_k, state_entropy, 2e-10);
        close(
            step.stored_energy_change_j + step.radiation_out_j,
            0.0,
            1e-8,
        );
        close(entropy.balance_residual_j_k, 0.0, 1e-9);
        assert!(entropy.radiation_generation_j_k > 0.0);
        assert_eq!(entropy.source_at_node_j_k, 0.0);
        assert_eq!(entropy.dirichlet_entropy_input_j_k, 0.0);
        assert_eq!(entropy.robin_entropy_input_j_k, 0.0);
    });
}
