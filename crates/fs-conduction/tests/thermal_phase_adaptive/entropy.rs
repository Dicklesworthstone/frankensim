//! Native accepted-pair entropy, schedule knots and transactional refusals.
use super::*;
use fs_conduction::bc::history::march::enthalpy::adaptive::EntropyAdaptiveHistoryError;
use fs_conduction::transient::enthalpy::entropy::{EnthalpyEntropyBalance, EnthalpyEntropyError};
use fs_conduction::transient::enthalpy::heterogeneous::EnthalpyRadiation;
use fs_conduction::transient::enthalpy::sdirk::{EnthalpySdirk2Input, SDIRK2_GAMMA};

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:.17e}, expected={expected:.17e}, allowance={tolerance:.3e}"
    );
}

fn plateau_account(account: EnthalpyEntropyBalance, heat_j: f64) {
    close(account.thermal_storage_change_j_k, heat_j / 310.0, 1e-10);
    close(account.source_at_node_j_k, heat_j / 310.0, 1e-10);
    close(account.conduction_generation_j_k, 0.0, 1e-10);
    close(account.temporal_defect_j_k, 0.0, 1e-10);
    close(account.algebraic_residual_j_k, 0.0, 1e-10);
    close(account.balance_residual_j_k, 0.0, 1e-10);
    assert_eq!(account.dirichlet_entropy_input_j_k, 0.0);
    assert_eq!(account.robin_entropy_input_j_k, 0.0);
    assert_eq!(account.radiation_entropy_input_j_k, 0.0);
    assert_eq!(account.robin_generation_j_k, 0.0);
    assert_eq!(account.radiation_generation_j_k, 0.0);
    account.admit(1e-9).unwrap();
}

fn initial_plateau(model: &EnthalpyHistoryRun<'_, '_, '_>) {
    assert_eq!(model.specific_enthalpy_j_kg(), &[150.0; 4]);
    assert_eq!(model.temperature(), &[310.0; 4]);
    assert_eq!(model.liquid_mass_fraction(), &[0.125; 4]);
    assert_eq!(model.time_s(), 0.0);
    assert_eq!(model.accepted_steps(), 0);
}

#[test]
fn g1_entropy_excludes_rejected_latent_trials_and_preserves_original_trajectory() {
    let mesh = mesh();
    let material = ConductivityModel::isotropic_declared(2.0).unwrap();
    with_cx(|cx| {
        let history = boundary(cx, &mesh);
        for offset in [0.0, 10000.0] {
            let chart = curve(offset);
            let source =
                || HeatSourceHistory::Uniform(ramp(0.0, 40.0, 1.0, VOLUMETRIC_SOURCE_DIMS));
            let mut reference = run(cx, &history, &material, &chart, offset, source());
            let expected = reference
                .advance_adaptive_step(cx, 1.0, 1.0, controls(), adaptive(), phase())
                .unwrap();
            let mut audited = run(cx, &history, &material, &chart, offset, source());
            let got = audited
                .advance_adaptive_step_with_entropy(
                    cx,
                    1.0,
                    1.0,
                    controls(),
                    adaptive(),
                    phase(),
                    1e-9,
                )
                .unwrap();

            assert_eq!(got.physical.end_time_s, 0.125);
            assert_eq!(got.physical.rejected_attempts, 3);
            assert_eq!(got.physical.endpoint_solves, 12);
            assert_eq!(got.physical.end_time_s, expected.end_time_s);
            assert_eq!(got.physical.normalized_error, expected.normalized_error);
            assert_eq!(
                got.physical.first_half.specific_enthalpy_j_kg,
                expected.first_half.specific_enthalpy_j_kg
            );
            assert_eq!(
                got.physical.second_half.specific_enthalpy_j_kg,
                expected.second_half.specific_enthalpy_j_kg
            );
            assert_eq!(
                audited.specific_enthalpy_j_kg(),
                reference.specific_enthalpy_j_kg()
            );
            assert_eq!(audited.accepted_steps(), 2);

            // V=1/6, rho=1, h'=40*t: only the accepted two endpoints supply
            // delta H=5*dt^2. The plateau stores latent entropy delta H/310.
            let dt = got.physical.end_time_s;
            let heat = 5.0 * dt * dt;
            plateau_account(got.half_entropy[0], heat / 3.0);
            plateau_account(got.half_entropy[1], 2.0 * heat / 3.0);
            plateau_account(got.entropy, heat);
            assert!(got.entropy.has_unclassified_inputs);
            close(
                got.physical.first_half.stored_energy_change_j
                    + got.physical.second_half.stored_energy_change_j,
                heat,
                1e-10,
            );
            assert_eq!(
                got.entropy.thermal_storage_change_j_k,
                audited
                    .storage()
                    .thermal_entropy_change_j_k(
                        cx,
                        &[150.0 + offset; 4],
                        audited.specific_enthalpy_j_kg()
                    )
                    .unwrap()
            );
        }
    });
}

#[test]
fn g1_sdirk_entropy_uses_one_fine_window_allowance_and_ignores_coarse_entropy() {
    let mesh = mesh();
    let material = ConductivityModel::isotropic_declared(2.0).unwrap();
    let chart = curve(0.0);
    with_cx(|cx| {
        let history = boundary(cx, &mesh);
        let base = EnthalpyHistoryRun::new(
            cx,
            problem(
                &history,
                &material,
                HeatSourceHistory::Fixed(ScalarField::Uniform(90.0)),
            ),
            EnthalpyHistoryMaterials {
                table: &[ReferenceEnthalpyMaterial {
                    curve: &chart,
                    reference_density_kg_m3: 1.0,
                }],
                of_element: &[0],
                phase_conductivity: None,
            },
            &[0.0; 4],
            EnthalpyHistoryBudget {
                spatial: EnthalpyBudget {
                    max_vertices: 4,
                    max_elements: 1,
                },
                max_steps: 2,
            },
        )
        .unwrap();
        // Cp=10, total mass=1/6, q=90: T(t)=300+9*t and total heat=15 J.
        // The actual SDIRK source entropy uses its two original stage temperatures.
        let capacity = 10.0 / 6.0;
        let source_entropy = |temperature: f64, dt: f64| {
            capacity
                * 9.0
                * dt
                * ((1.0 - SDIRK2_GAMMA) / (temperature + SDIRK2_GAMMA * 9.0 * dt)
                    + SDIRK2_GAMMA / (temperature + 9.0 * dt))
        };
        let first_source = source_entropy(300.0, 0.5);
        let second_source = source_entropy(304.5, 0.5);
        let first_defect = capacity * (304.5_f64 / 300.0).ln() - first_source;
        let second_defect = capacity * (309.0_f64 / 304.5).ln() - second_source;
        let storage = capacity * (309.0_f64 / 300.0).ln();
        let total_defect = storage - first_source - second_source;
        let strict = 1.1 * (-first_defect).max(-second_defect);
        let lenient = 1.2 * -total_defect;
        assert!(strict < -total_defect && total_defect < 0.0);

        // This rejected coarse account must not block the admissible fine path.
        let boundary = history.at(cx, 0.0).unwrap();
        let source = ScalarField::Uniform(90.0);
        let input = EnthalpySdirk2Input {
            problem: ConductionProblem {
                mesh: &mesh,
                boundary: &boundary,
                material: &material,
                source: &source,
                element_materials: None,
            },
            interfaces: None,
            contacts: None,
            radiation: EnthalpyRadiation::None,
        };
        assert!(matches!(
            base.storage().advance_sdirk2_with_entropy(
                cx,
                [input; 2],
                &[0.0; 4],
                1.0,
                controls(),
                lenient
            ),
            Err(EnthalpyEntropyError::NegativeProduction { .. })
        ));

        let mut audited = base.into_sdirk2();
        let refused = audited.advance_adaptive_step_with_entropy(
            cx,
            1.0,
            1.0,
            controls(),
            adaptive(),
            phase(),
            strict,
        );
        assert!(matches!(
            refused,
            Err(EntropyAdaptiveHistoryError::Entropy(
                EnthalpyEntropyError::NegativeProduction {
                    account: "thermal storage minus inputs",
                    ..
                }
            ))
        ));
        assert_eq!(audited.conduction().specific_enthalpy_j_kg(), &[0.0; 4]);
        assert_eq!(audited.conduction().temperature(), &[300.0; 4]);
        assert_eq!(audited.conduction().time_s(), 0.0);
        assert_eq!(audited.conduction().accepted_steps(), 0);

        let accepted = audited
            .advance_adaptive_step_with_entropy(
                cx,
                1.0,
                1.0,
                controls(),
                adaptive(),
                phase(),
                lenient,
            )
            .unwrap();
        for half in accepted.half_entropy {
            half.admit(strict).unwrap();
        }
        assert!(accepted.entropy.admit(strict).is_err());
        accepted.entropy.admit(lenient).unwrap();
        assert_eq!(accepted.physical.rejected_attempts, 0);
        assert_eq!(accepted.physical.endpoint_solves, 3);
        assert_eq!(audited.conduction().accepted_steps(), 2);
        assert!(
            audited
                .conduction()
                .temperature()
                .iter()
                .all(|t| (*t - 309.0).abs() < 1e-9)
        );
        close(
            accepted.physical.first_half.source_j + accepted.physical.second_half.source_j,
            15.0,
            1e-10,
        );
        close(accepted.entropy.thermal_storage_change_j_k, storage, 1e-11);
        close(
            accepted.entropy.source_at_node_j_k,
            first_source + second_source,
            1e-11,
        );
        close(accepted.entropy.temporal_defect_j_k, total_defect, 1e-11);
        close(
            accepted.half_entropy[0].temporal_defect_j_k,
            first_defect,
            1e-11,
        );
        close(
            accepted.half_entropy[1].temporal_defect_j_k,
            second_defect,
            1e-11,
        );
        close(accepted.entropy.conduction_generation_j_k, 0.0, 1e-11);
        close(accepted.entropy.balance_residual_j_k, 0.0, 1e-11);
        assert!(accepted.entropy.has_unclassified_inputs);
    });
}

#[test]
fn g1_entropy_adaptive_histories_keep_held_pulse_left_limits_and_resume_budget() {
    let mesh = mesh();
    let material = ConductivityModel::isotropic_declared(2.0).unwrap();
    let chart = curve(0.0);
    with_cx(|cx| {
        let history = boundary(cx, &mesh);
        let source = || {
            HeatSourceHistory::Uniform(TimeSignal::Table {
                times: vec![0.0, 0.25, 0.5],
                values: vec![120.0, 0.0, 90.0],
                dims: VOLUMETRIC_SOURCE_DIMS,
                interp: Interp::Hold,
            })
        };
        let mut be = run(cx, &history, &material, &chart, 0.0, source());
        for (end, heat) in [(0.25, 5.0), (0.5, 0.0)] {
            let got = be
                .advance_adaptive_step_with_entropy(
                    cx,
                    1.0,
                    1.0,
                    controls(),
                    adaptive(),
                    phase(),
                    1e-9,
                )
                .unwrap();
            assert_eq!(got.physical.end_time_s, end);
            assert_eq!(got.physical.rejected_attempts, 0);
            let h1 = got.physical.midpoint_time_s - got.physical.start_time_s;
            let h2 = got.physical.end_time_s - got.physical.midpoint_time_s;
            close(
                h1 * got.physical.first_half.source_w + h2 * got.physical.second_half.source_w,
                heat,
                1e-10,
            );
            plateau_account(got.half_entropy[0], heat / 2.0);
            plateau_account(got.half_entropy[1], heat / 2.0);
            plateau_account(got.entropy, heat);
        }
        assert_eq!(be.time_s(), 0.5);
        assert_eq!(be.accepted_steps(), 4);
        let retained = be.specific_enthalpy_j_kg().to_vec();
        assert!(retained.iter().all(|h| (*h - 180.0).abs() < 1e-9));
        assert!(
            be.advance_adaptive_step_with_entropy(
                cx,
                1.0,
                1.0,
                controls(),
                adaptive(),
                phase(),
                1e-9
            )
            .is_err()
        );
        assert_eq!(be.specific_enthalpy_j_kg(), retained.as_slice());
        assert_eq!(be.time_s(), 0.5);

        let mut sdirk = run(cx, &history, &material, &chart, 0.0, source()).into_sdirk2();
        for (end, heat) in [(0.25, 5.0), (0.5, 0.0)] {
            let got = sdirk
                .advance_adaptive_step_with_entropy(
                    cx,
                    1.0,
                    1.0,
                    controls(),
                    adaptive(),
                    phase(),
                    1e-9,
                )
                .unwrap();
            assert_eq!(got.physical.end_time_s, end);
            assert_eq!(got.physical.rejected_attempts, 0);
            close(
                got.physical.first_half.source_j + got.physical.second_half.source_j,
                heat,
                1e-10,
            );
            plateau_account(got.half_entropy[0], heat / 2.0);
            plateau_account(got.half_entropy[1], heat / 2.0);
            plateau_account(got.entropy, heat);
        }
        assert_eq!(sdirk.conduction().time_s(), 0.5);
        assert_eq!(sdirk.conduction().accepted_steps(), 4);
        let retained = sdirk.conduction().specific_enthalpy_j_kg().to_vec();
        assert!(retained.iter().all(|h| (*h - 180.0).abs() < 1e-9));
        assert!(
            sdirk
                .advance_adaptive_step_with_entropy(
                    cx,
                    1.0,
                    1.0,
                    controls(),
                    adaptive(),
                    phase(),
                    1e-9
                )
                .is_err()
        );
        assert_eq!(
            sdirk.conduction().specific_enthalpy_j_kg(),
            retained.as_slice()
        );
        assert_eq!(sdirk.conduction().time_s(), 0.5);
    });
}

#[test]
fn g0_entropy_pair_refusal_after_a_valid_first_half_is_atomic_and_cancel_safe() {
    let mesh = mesh();
    let material = ConductivityModel::isotropic_declared(2.0).unwrap();
    let chart = curve(0.0);
    let history = with_cx(|cx| boundary(cx, &mesh));
    let mut model = with_cx(|cx| {
        run(
            cx,
            &history,
            &material,
            &chart,
            0.0,
            HeatSourceHistory::Uniform(ramp(2400.0, 800.0, 1.0, VOLUMETRIC_SOURCE_DIMS)),
        )
    });
    with_cancelled_cx(|cx| {
        assert!(
            model
                .advance_adaptive_step_with_entropy(
                    cx,
                    1.0,
                    1.0,
                    controls(),
                    adaptive(),
                    phase(),
                    1e-9
                )
                .is_err()
        )
    });
    initial_plateau(&model);
    with_cx(|cx| {
        let boundary = history.at(cx, 0.0).unwrap();
        let source = ScalarField::Uniform(1600.0);
        let problem = ConductionProblem {
            mesh: &mesh,
            boundary: &boundary,
            material: &material,
            source: &source,
            element_materials: None,
        };
        // Coarse h=950 and first-half h=950 are inside the original chart.
        // The second half requires h=1350, beyond its upper support of 1200.
        let (first, _) = model
            .storage()
            .advance_with_entropy(cx, problem, None, &[150.0; 4], 0.5, controls(), 1e-9)
            .unwrap();
        assert!(
            first
                .specific_enthalpy_j_kg
                .iter()
                .all(|h| (*h - 950.0).abs() < 1e-8)
        );
        assert!(
            model
                .advance_adaptive_step_with_entropy(
                    cx,
                    1.0,
                    1.0,
                    controls(),
                    adaptive(),
                    phase(),
                    1e-9
                )
                .is_err()
        );
        initial_plateau(&model);
    });
}
