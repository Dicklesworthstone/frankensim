//! Adaptive contact/radiation with an independent continuous-time P1 reference.
use super::*;
use fs_conduction::transient::enthalpy::adaptive::contact::AdaptiveContactEnthalpyConfig;
use fs_conduction::transient::enthalpy::adaptive::{
    AdaptiveEnthalpyConfig, AdaptiveEnthalpyStatus,
};
use fs_conduction::{
    AmbientRadiationPatch, PointwiseAmbientRadiation, PointwiseEnclosureRadiation,
    PointwiseEnclosureTolerance,
};
use std::ops::ControlFlow;

fn time_config(h_tolerance: f64) -> AdaptiveEnthalpyConfig {
    AdaptiveEnthalpyConfig {
        enthalpy_absolute_tolerance_j_kg: h_tolerance,
        temperature_absolute_tolerance_k: h_tolerance / 10.,
        liquid_fraction_absolute_tolerance: h_tolerance / 100.,
        relative_tolerance: 0.,
        minimum_step_s: 1e-8,
        maximum_step_s: 0.5,
        endpoint: coupled_newton().step,
    }
}

fn contacts_config() -> AdaptiveContactEnthalpyConfig {
    let original = coupled_newton();
    AdaptiveContactEnthalpyConfig {
        max_residual_evaluations_per_endpoint: original.max_residual_evaluations,
        max_face_pairs: original.max_face_pairs,
        residual_tolerance_j: original.residual_tolerance_j,
    }
}

fn bits(values: &[f64]) -> Vec<u64> {
    values.iter().map(|value| value.to_bits()).collect()
}

// Independent time-continuous eight-node P1 balance, integrated by RK4 with
// 100,000 equal steps; halving dt changed every h by less than 1.5e-9 J/kg.
// The reference uses explicit tetrahedral k=2 stiffness, diagonal masses
// rho*V/4, consistent triangle contact blocks (1+delta)/24R, R(mean T), and
// the two piecewise-linear phase charts. No production endpoint supplies it.
const CONTINUOUS_H_AT_ONE_SECOND: [f64; 8] = [
    310.5863459823301,
    331.1459045221456,
    331.1459045221456,
    121.51465628348586,
    368.4550936466141,
    354.06048422385055,
    354.06048422385055,
    456.22753225066884,
];

#[test]
fn g1_adaptive_contact_resolves_spatial_melting_against_independent_time_continuous_balance() {
    context(&CancelGate::new_clock_free(), |cx| {
        let mesh = mesh();
        let boundary = ThermalBoundaryBuilder::new(&mesh)
            .adiabatic_remainder()
            .finish()
            .unwrap();
        let card = card();
        let contacts =
            TemperatureContacts::new(cx, &mesh, &boundary, vec![], vec![law(&mesh, &card)], 1)
                .unwrap();
        let cold = chart(false);
        let hot = chart(true);
        let engine = HeterogeneousEnthalpyBackwardEuler::new(
            cx,
            &mesh,
            &[
                ReferenceEnthalpyMaterial {
                    curve: &cold,
                    reference_density_kg_m3: 1.,
                },
                ReferenceEnthalpyMaterial {
                    curve: &hot,
                    reference_density_kg_m3: 2.,
                },
            ],
            &[0, 1],
            EnthalpyBudget {
                max_vertices: 8,
                max_elements: 2,
            },
        )
        .unwrap();
        let material = ConductivityModel::isotropic_declared(2.).unwrap();
        let source = ScalarField::Nodal(vec![0., 0., 0., 0., 60., 60., 60., 60.]);
        let problem = ConductionProblem {
            mesh: &mesh,
            boundary: &boundary,
            material: &material,
            element_materials: None,
            source: &source,
        };
        let initial = [80., 80., 80., 80., 450., 450., 450., 450.];
        let mut errors = Vec::new();
        let mut accepted = Vec::new();
        for tolerance in [0.2, 0.002] {
            let mut state = engine
                .adaptive_state(cx, &initial, 0., 0.5, time_config(tolerance))
                .unwrap();
            let report = state
                .advance_to_with_temperature_contacts(
                    cx,
                    problem,
                    &contacts,
                    contacts_config(),
                    1.,
                    8192,
                )
                .unwrap();
            assert_eq!(report.status, AdaptiveEnthalpyStatus::ReachedEnd);
            assert_eq!(state.time_s().to_bits(), 1_f64.to_bits());
            assert!(report.rejected > 0);
            assert!(report.endpoint_solves >= 3 * report.accepted + report.rejected);
            assert!(report.endpoint_solves <= 3 * (report.accepted + report.rejected));
            assert_eq!(state.endpoint_solves(), report.endpoint_solves);
            assert!(report.last_error.unwrap().normalized <= 1.);
            close(report.integrated_net_input_j, 10., 1e-9);
            close(report.stored_energy_change_j, 10., 2e-7);
            close(report.energy_residual_j, 0., 2e-7);
            // The far cold node is still melting while the contact nodes have
            // left the plateau. Temperature alone cannot control this error.
            close(state.temperature()[3], 310., 1e-10);
            assert!((0.1..0.4).contains(&state.liquid_mass_fraction()[3]));
            assert_eq!(state.liquid_mass_fraction()[0], 1.);
            let error = state
                .specific_enthalpy_j_kg()
                .iter()
                .zip(CONTINUOUS_H_AT_ONE_SECOND)
                .map(|(&h, reference)| (h - reference).abs())
                .fold(0_f64, f64::max);
            errors.push(error);
            accepted.push(report.accepted);
        }
        assert!(errors[1] < 0.5, "fine time error = {} J/kg", errors[1]);
        assert!(
            errors[1] < 0.25 * errors[0],
            "refinement errors: {errors:?}"
        );
        assert!(accepted[1] > accepted[0]);
        assert_eq!(initial, [80., 80., 80., 80., 450., 450., 450., 450.]);
    });
}

#[test]
fn g3_adaptive_radiative_contact_observer_resume_preserves_accepted_heat_and_exact_history() {
    context(&CancelGate::new_clock_free(), |cx| {
        let mesh = mesh();
        let boundary = pointwise::ambient_boundary(&mesh);
        let card = card();
        let contacts =
            TemperatureContacts::new(cx, &mesh, &boundary, vec![], vec![law(&mesh, &card)], 1)
                .unwrap();
        let cold = chart(false);
        let hot = chart(true);
        let engine = HeterogeneousEnthalpyBackwardEuler::new(
            cx,
            &mesh,
            &[
                ReferenceEnthalpyMaterial {
                    curve: &cold,
                    reference_density_kg_m3: 1.,
                },
                ReferenceEnthalpyMaterial {
                    curve: &hot,
                    reference_density_kg_m3: 2.,
                },
            ],
            &[0, 1],
            EnthalpyBudget {
                max_vertices: 8,
                max_elements: 2,
            },
        )
        .unwrap();
        let material = ConductivityModel::isotropic_declared(2.).unwrap();
        let source = ScalarField::Uniform(0.);
        let problem = ConductionProblem {
            mesh: &mesh,
            boundary: &boundary,
            material: &material,
            element_materials: None,
            source: &source,
        };
        let patches =
            [
                AmbientRadiationPatch::new(
                    "exterior",
                    pointwise::emissivity("exterior", 0.2),
                    340.,
                )
                .unwrap(),
            ];
        let radiation = PointwiseAmbientRadiation::new(cx, &mesh, &boundary, &patches).unwrap();
        let initial = [80., 80., 80., 80., 450., 450., 450., 450.];
        let mut whole = engine
            .adaptive_state(cx, &initial, 0., 0.2, time_config(0.01))
            .unwrap();
        let mut resumed = whole.clone();
        let expected = whole
            .advance_to_with_temperature_contacts_and_pointwise_ambient_radiation(
                cx,
                problem,
                &contacts,
                contacts_config(),
                &radiation,
                0.2,
                4096,
            )
            .unwrap();
        let mut visits = 0;
        let mut integrated = 0.;
        let mut previous_time = 0.0_f64;
        let mut observe =
            |trial: fs_conduction::transient::enthalpy::adaptive::EnthalpyAcceptedTrial<'_>| {
                visits += 1;
                assert_eq!(trial.start_time_s.to_bits(), previous_time.to_bits());
                assert!(trial.error.normalized <= 1.);
                for (half, dt) in [
                    (trial.first_half, trial.midpoint_time_s - trial.start_time_s),
                    (trial.second_half, trial.end_time_s - trial.midpoint_time_s),
                ] {
                    let radiation = half.pointwise_radiation.as_ref().unwrap();
                    assert!(radiation.outward_heat_w < 0.);
                    assert!(half.robin_out_w > 0.);
                    integrated -= dt * (half.robin_out_w + radiation.outward_heat_w);
                    close(
                        half.energy_residual_j,
                        0.,
                        coupled_newton().step.energy_tolerance_j,
                    );
                    assert_eq!(half.contact_fluxes.len(), 1);
                }
                previous_time = trial.end_time_s;
                if visits == 1 {
                    ControlFlow::Break(())
                } else {
                    ControlFlow::Continue(())
                }
            };
        let first = resumed
            .advance_to_with_temperature_contacts_and_pointwise_ambient_radiation_observed(
                cx,
                problem,
                &contacts,
                contacts_config(),
                &radiation,
                0.2,
                4096,
                &mut observe,
            )
            .unwrap();
        assert_eq!(first.status, AdaptiveEnthalpyStatus::ObserverStop);
        assert_eq!(first.accepted, 1);
        let rest = resumed
            .advance_to_with_temperature_contacts_and_pointwise_ambient_radiation_observed(
                cx,
                problem,
                &contacts,
                contacts_config(),
                &radiation,
                0.2,
                4096,
                &mut observe,
            )
            .unwrap();
        assert_eq!(rest.status, AdaptiveEnthalpyStatus::ReachedEnd);
        assert_eq!(expected.status, AdaptiveEnthalpyStatus::ReachedEnd);
        assert_eq!(visits, resumed.accepted_steps());
        assert_eq!(previous_time.to_bits(), 0.2_f64.to_bits());
        assert_eq!(
            bits(resumed.specific_enthalpy_j_kg()),
            bits(whole.specific_enthalpy_j_kg())
        );
        assert_eq!(
            bits(resumed.liquid_mass_fraction()),
            bits(whole.liquid_mass_fraction())
        );
        assert_eq!(
            resumed.next_step_s().to_bits(),
            whole.next_step_s().to_bits()
        );
        assert_eq!(
            first.endpoint_solves + rest.endpoint_solves,
            expected.endpoint_solves
        );
        close(rest.integrated_net_input_j, integrated, 1e-11);
        close(rest.stored_energy_change_j, integrated, 1e-7);
        assert!(rest.stored_energy_change_j > 0.);
    });
}

#[test]
fn g1_adaptive_closed_enclosure_and_contact_exchange_do_not_create_external_heat() {
    context(&CancelGate::new_clock_free(), |cx| {
        let mesh = mesh();
        let boundary = pointwise::enclosure_boundary(&mesh);
        let card = card();
        let contacts =
            TemperatureContacts::new(cx, &mesh, &boundary, vec![], vec![law(&mesh, &card)], 1)
                .unwrap();
        let cold = chart(false);
        let hot = chart(true);
        let engine = HeterogeneousEnthalpyBackwardEuler::new(
            cx,
            &mesh,
            &[
                ReferenceEnthalpyMaterial {
                    curve: &cold,
                    reference_density_kg_m3: 1.,
                },
                ReferenceEnthalpyMaterial {
                    curve: &hot,
                    reference_density_kg_m3: 2.,
                },
            ],
            &[0, 1],
            EnthalpyBudget {
                max_vertices: 8,
                max_elements: 2,
            },
        )
        .unwrap();
        let material = ConductivityModel::isotropic_declared(2.).unwrap();
        let source = ScalarField::Uniform(0.);
        let problem = ConductionProblem {
            mesh: &mesh,
            boundary: &boundary,
            material: &material,
            element_materials: None,
            source: &source,
        };
        let enclosure = pointwise::enclosure(&mesh, [0.2, 0.7]);
        let radiation = PointwiseEnclosureRadiation::new(
            cx,
            &mesh,
            &boundary,
            &enclosure,
            PointwiseEnclosureTolerance::default(),
        )
        .unwrap();
        let initial = [80., 80., 80., 80., 450., 450., 450., 450.];
        let mut state = engine
            .adaptive_state(cx, &initial, 0., 0.2, time_config(0.02))
            .unwrap();
        let report = state
            .advance_to_with_temperature_contacts_and_pointwise_enclosure_radiation(
                cx,
                problem,
                &contacts,
                contacts_config(),
                &radiation,
                0.2,
                4096,
            )
            .unwrap();
        assert_eq!(report.status, AdaptiveEnthalpyStatus::ReachedEnd);
        close(report.integrated_net_input_j, 0., 1e-9);
        close(report.stored_energy_change_j, 0., 1e-7);
        close(report.energy_residual_j, 0., 1e-7);
        assert!(state.specific_enthalpy_j_kg()[0] > initial[0]);
        assert!(state.specific_enthalpy_j_kg()[4] < initial[4]);
        let radiation = state
            .last_endpoint()
            .unwrap()
            .pointwise_enclosure_radiation
            .as_ref()
            .unwrap();
        assert!(radiation.patches[0].outward_heat_w < -1.);
        assert!(radiation.patches[1].outward_heat_w > 1.);
        close(radiation.outward_heat_w, 0., 1e-9);
    });
}

#[test]
fn g4_adaptive_contact_bounds_failed_work_and_retries_nonlinear_refusals_without_partial_history() {
    context(&CancelGate::new_clock_free(), |cx| {
        let mesh = mesh();
        let boundary = ThermalBoundaryBuilder::new(&mesh)
            .adiabatic_remainder()
            .finish()
            .unwrap();
        let card = card();
        let contacts =
            TemperatureContacts::new(cx, &mesh, &boundary, vec![], vec![law(&mesh, &card)], 1)
                .unwrap();
        let cold = chart(false);
        let hot = chart(true);
        let engine = HeterogeneousEnthalpyBackwardEuler::new(
            cx,
            &mesh,
            &[
                ReferenceEnthalpyMaterial {
                    curve: &cold,
                    reference_density_kg_m3: 1.,
                },
                ReferenceEnthalpyMaterial {
                    curve: &hot,
                    reference_density_kg_m3: 2.,
                },
            ],
            &[0, 1],
            EnthalpyBudget {
                max_vertices: 8,
                max_elements: 2,
            },
        )
        .unwrap();
        let material = ConductivityModel::isotropic_declared(2.).unwrap();
        let source = ScalarField::Uniform(0.);
        let problem = ConductionProblem {
            mesh: &mesh,
            boundary: &boundary,
            material: &material,
            element_materials: None,
            source: &source,
        };
        let initial = [80., 80., 80., 80., 450., 450., 450., 450.];
        let mut state = engine
            .adaptive_state(cx, &initial, 0., 0.2, time_config(0.02))
            .unwrap();
        let mut controls = contacts_config();
        controls.max_residual_evaluations_per_endpoint = 1;
        assert!(matches!(
            state.advance_to_with_temperature_contacts(cx, problem, &contacts, controls, 0.2, 8),
            Err(EnthalpyError::InvalidInput(
                "contact enthalpy residual allowance exhausted"
            ))
        ));
        assert_eq!(state.endpoint_solves(), 1);
        assert_eq!(state.accepted_steps(), 0);
        assert_eq!(state.time_s(), 0.);
        assert_eq!(state.specific_enthalpy_j_kg(), initial);
        assert!(state.last_endpoint().is_none());
        assert!(matches!(
            state.advance_to_with_temperature_contacts(
                cx,
                problem,
                &contacts,
                contacts_config(),
                0.2,
                usize::MAX
            ),
            Err(EnthalpyError::InvalidInput(
                "adaptive contact work allowance overflows"
            ))
        ));
        assert_eq!(state.endpoint_solves(), 1);
        let mut limited = time_config(0.02);
        limited.endpoint.max_newton_iterations = 1;
        limited.minimum_step_s = 0.1;
        let mut refused = engine
            .adaptive_state(cx, &initial, 0., 0.2, limited)
            .unwrap();
        let report = refused
            .advance_to_with_temperature_contacts(cx, problem, &contacts, contacts_config(), 0.2, 8)
            .unwrap();
        assert_eq!(report.status, AdaptiveEnthalpyStatus::MinimumStep);
        assert_eq!(
            (report.accepted, report.rejected, report.endpoint_solves),
            (0, 2, 2)
        );
        assert!(matches!(
            report.last_endpoint_refusal,
            Some(EnthalpyError::NotConverged(_))
        ));
        assert_eq!(refused.specific_enthalpy_j_kg(), initial);
        assert_eq!(refused.time_s(), 0.);
        close(report.integrated_net_input_j, 0., 0.);
        let recovered = state
            .advance_to_with_temperature_contacts(
                cx,
                problem,
                &contacts,
                contacts_config(),
                0.05,
                4096,
            )
            .unwrap();
        assert_eq!(recovered.status, AdaptiveEnthalpyStatus::ReachedEnd);
        assert_eq!(state.endpoint_solves(), recovered.endpoint_solves + 1);
        close(recovered.energy_residual_j, 0., 1e-7);
    });
}
