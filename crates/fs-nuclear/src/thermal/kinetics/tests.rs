// A child of thermal::tests reuses its original caloric, cell and isotope fixtures.
use super::*;
use crate::kinetics::{
    DelayedNeutronGroup, KineticsBudget, KineticsError, PointKineticsModel, PointKineticsState,
    TemperatureFeedback, TemperatureFeedbackLaw,
};

fn kinetic_budget() -> KineticHeatingBudget {
    KineticHeatingBudget {
        irradiation: irradiation_budget(),
        kinetics: KineticsBudget {
            max_delayed_groups: 8,
            max_squarings: 32,
            max_taylor_terms: 24,
            population_tolerance: 1e-8,
            integral_tolerance_population_s: 1e-8,
        },
        max_feedback_terms: 4,
    }
}

fn kinetic_source(reactivity: f64, delayed_heat: bool) -> KineticCellSource {
    let model = Arc::new(
        PointKineticsModel::new("manufactured-normalized-region-v1".into(), 1.0, vec![]).unwrap(),
    );
    KineticCellSource {
        irradiation: irradiated_source(1000.0, delayed_heat),
        kinetics: PointKineticsState::new(model, 1.0, vec![]).unwrap(),
        controls: KineticsControlSchedule::new(
            KineticsControl {
                reactivity,
                external_source_per_s: 0.0,
            },
            vec![],
        )
        .unwrap(),
        feedback: vec![],
    }
}

fn feedback(coefficient: f64, maximum_k: f64) -> TemperatureFeedback {
    TemperatureFeedback::new(
        "manufactured-cell-feedback-v1".into(),
        300.0,
        290.0,
        maximum_k,
        TemperatureFeedbackLaw::Linear {
            coefficient_per_kelvin: coefficient,
        },
    )
    .unwrap()
}

#[test]
fn g1_neutron_population_integral_drives_actual_reaction_counts_and_heat() {
    let initial = thermal(
        vec![cell(1, 1, 2.0, 1200.0, 3.0)],
        &[(MaterialId(1), &HOT)],
        vec![],
    );
    let heating =
        KineticHeating::new(&initial, vec![kinetic_source(0.5, false)], kinetic_budget()).unwrap();
    let (heated, source, audit) = heating
        .advance_to(&initial, 1.0, controls(10.0), || Ok(()))
        .unwrap();
    // beta=0, Lambda=1 and rho=0.5 give n=exp(0.5t). The real
    // parent survival is exp(-0.2 integral(n dt)), not exp(-0.2 n_end dt).
    let integrated = 0.5_f64.exp_m1() / 0.5;
    let loss = -(-0.2 * integrated).exp_m1();
    let expected_fissions = 95.0 * loss;
    let expected_captures = 85.0 * loss;
    let deposited = 8.0 * expected_fissions + 3.0 * expected_captures;
    let report = source.last_step(CellId(1)).unwrap();
    close(report.kinetics.integrated_population_s, integrated, 2e-13);
    close(report.mean_neutron_population, integrated, 2e-13);
    close(
        source
            .source(CellId(1))
            .unwrap()
            .kinetics
            .neutron_population(),
        0.5_f64.exp(),
        2e-13,
    );
    close(source.audit().fissions, expected_fissions, 2e-11);
    close(source.audit().captures, expected_captures, 2e-11);
    close(audit.external_heat_j, deposited + 3.0, 2e-10);
    close(
        heated.temperature_k(CellId(1)).unwrap(),
        300.0 + (deposited + 3.0) / 8.0,
        2e-11,
    );
    assert_eq!(heated.cells()[&CellId(1)].boundary.source_w, 3.0);
    assert_eq!(heated.cells()[&CellId(1)].mass_kg, 2.0);
    assert_eq!(audit.substeps, 1);
}

#[test]
fn g3_temperature_feedback_closes_the_loop_and_converges_to_an_independent_invariant() {
    let initial = thermal(
        vec![cell(1, 1, 2.0, 1200.0, 0.0)],
        &[(MaterialId(1), &HOT)],
        vec![],
    );
    let reference_rho = 0.05;
    let coefficient = -0.02;
    let mut coupled = kinetic_source(reference_rho, false);
    coupled.feedback.push(feedback(coefficient, 1000.0));
    let heating = KineticHeating::new(&initial, vec![coupled], kinetic_budget()).unwrap();
    let unregulated = KineticHeating::new(
        &initial,
        vec![kinetic_source(reference_rho, false)],
        kinetic_budget(),
    )
    .unwrap();
    let (_, unregulated, _) = unregulated
        .advance_to(&initial, 1.0, controls(0.025), || Ok(()))
        .unwrap();
    let mut errors = Vec::new();
    for dt in [0.1, 0.05, 0.025] {
        let (heated, state, _) = heating
            .advance_to(&initial, 1.0, controls(dt), || Ok(()))
            .unwrap();
        let source = state.source(CellId(1)).unwrap();
        let neutron_population = source.kinetics.neutron_population();
        let integrated =
            -(source.irradiation.inventories[0].remaining_parent_nuclei() / 100.0).ln() / 0.2;
        // For these two parents the deposited heat is 1015*(1-exp(-0.2I)).
        // Heat capacity is 8 J/K, so eliminating t from n'=rho(T)n gives
        // n(I)=1+rho0*I+alpha*(1015/8)*(I+expm1(-0.2I)/0.2).
        let invariant_n = 1.0
            + reference_rho * integrated
            + coefficient * (1015.0 / 8.0) * (integrated + (-0.2 * integrated).exp_m1() / 0.2);
        errors.push((neutron_population - invariant_n).abs());
        assert!(heated.temperature_k(CellId(1)).unwrap() > 300.0);
        assert!(neutron_population < 1.0);
        assert!(
            neutron_population
                < unregulated
                    .source(CellId(1))
                    .unwrap()
                    .kinetics
                    .neutron_population()
        );
        assert!(state.audit().fissions < unregulated.audit().fissions);
        assert!(state.last_step(CellId(1)).unwrap().applied_reactivity < reference_rho);
        close(
            heated.energy_j() - initial.energy_j(),
            state.cumulative_released().unwrap().deposited_j,
            1e-8,
        );
    }
    for pair in errors.windows(2) {
        assert!(
            (1.8..2.2).contains(&(pair[0] / pair[1])),
            "feedback refinement errors {errors:?}"
        );
    }
    assert!(errors[2] < 0.006, "feedback refinement errors {errors:?}");
}

#[test]
fn g3_controls_reference_flux_shutdown_and_checkpoint_restart_share_one_timeline() {
    let initial = thermal(
        vec![cell(1, 1, 2.0, 1200.0, 0.0)],
        &[(MaterialId(1), &HOT)],
        vec![],
    );
    let mut source = kinetic_source(0.0, true);
    let model = Arc::new(
        PointKineticsModel::new(
            "manufactured-precursor-region-v1".into(),
            0.1,
            vec![DelayedNeutronGroup {
                delayed_fraction: 0.1,
                decay_constant_s_inv: 0.2,
            }],
        )
        .unwrap(),
    );
    source.kinetics = PointKineticsState::new(
        Arc::clone(&model),
        1.0,
        model.equilibrium_precursors(1.0).unwrap(),
    )
    .unwrap();
    source.controls = KineticsControlSchedule::new(
        KineticsControl {
            reactivity: 0.0,
            external_source_per_s: 0.0,
        },
        vec![KineticsControlChange {
            time_s: 1.0,
            control: KineticsControl {
                reactivity: -0.1,
                external_source_per_s: 0.0,
            },
        }],
    )
    .unwrap();
    let heating = KineticHeating::new(&initial, vec![source], kinetic_budget()).unwrap();
    let (whole_t, whole_s, audit) = heating
        .advance_to(&initial, 5.0, controls(10.0), || Ok(()))
        .unwrap();
    assert_eq!(
        audit.substeps, 3,
        "control change at1, reference-flux shutdown at2, target at5"
    );
    let (at1, s1, _) = heating
        .advance_to(&initial, 1.0, controls(10.0), || Ok(()))
        .unwrap();
    let (at2, s2, _) = s1.advance_to(&at1, 2.0, controls(10.0), || Ok(())).unwrap();
    let (split_t, split_s, after) = s2.advance_to(&at2, 5.0, controls(10.0), || Ok(())).unwrap();
    let whole = whole_s.source(CellId(1)).unwrap();
    let split = split_s.source(CellId(1)).unwrap();
    assert_eq!(whole_t.cells(), split_t.cells());
    assert_eq!(whole.kinetics, split.kinetics);
    assert_eq!(whole.irradiation.inventories, split.irradiation.inventories);
    assert_eq!(split_s.audit().fissions, 0.0);
    assert_eq!(split_s.audit().captures, 0.0);
    assert!(
        after.external_heat_j > 0.0,
        "effective decay heat persists after incident flux is zero"
    );
    assert!(split.kinetics.neutron_population() > 0.0);
    assert!(split.kinetics.precursor_populations()[0] > 0.0);
}

#[test]
fn g4_feedback_endpoint_and_nested_checkpoint_refusal_publish_no_coupled_history() {
    let initial = thermal(
        vec![cell(1, 1, 2.0, 1200.0, 0.0)],
        &[(MaterialId(1), &HOT)],
        vec![],
    );
    let mut source = kinetic_source(0.0, false);
    source.feedback.push(feedback(-0.001, 300.1));
    let original_kinetics = source.kinetics.clone();
    let original_inventory = source.irradiation.inventories.clone();
    let heating = KineticHeating::new(&initial, vec![source], kinetic_budget()).unwrap();
    assert!(matches!(
        heating.advance_to(&initial, 0.5, controls(1.0), || Ok(())),
        Err(FissionThermalError::Kinetics(KineticsError::InvalidInput(
            _
        )))
    ));
    assert_eq!(initial.time_s(), 0.0);
    assert_eq!(initial.temperature_k(CellId(1)).unwrap(), 300.0);
    assert_eq!(heating.time_s(), 0.0);
    assert_eq!(
        heating.source(CellId(1)).unwrap().kinetics,
        original_kinetics
    );
    assert_eq!(
        heating.source(CellId(1)).unwrap().irradiation.inventories,
        original_inventory
    );

    let source = kinetic_source(0.1, false);
    let mut polls = 0;
    let failure = source.trial_source(0.0, 1.0, 300.0, kinetic_budget(), &mut || {
        polls += 1;
        if polls == 5 {
            Err(FissionError::InvalidInput("checkpoint marker"))
        } else {
            Ok(())
        }
    });
    assert!(matches!(
        failure,
        Err(FissionThermalError::Fission(FissionError::InvalidInput(
            "checkpoint marker"
        )))
    ));
    assert_eq!(source.kinetics.neutron_population(), 1.0);
    assert_eq!(source.irradiation.inventories[0].cumulative_fissions(), 0.0);
}

#[test]
fn g0_control_admission_feedback_budget_and_idle_source_have_explicit_semantics() {
    let initial = KineticsControl {
        reactivity: 0.1,
        external_source_per_s: 0.0,
    };
    let changed = KineticsControl {
        reactivity: -0.2,
        external_source_per_s: 0.4,
    };
    let schedule = KineticsControlSchedule::new(
        initial,
        vec![KineticsControlChange {
            time_s: 1.0,
            control: changed,
        }],
    )
    .unwrap();
    assert_eq!(schedule.control_at(0.999).unwrap(), initial);
    assert_eq!(schedule.control_at(1.0).unwrap(), changed);
    assert!(schedule.control_at(f64::NAN).is_err());
    assert!(
        KineticsControlSchedule::new(
            initial,
            vec![
                KineticsControlChange {
                    time_s: 1.0,
                    control: changed
                },
                KineticsControlChange {
                    time_s: 1.0,
                    control: initial
                },
            ]
        )
        .is_err()
    );
    assert!(
        KineticsControlSchedule::new(
            KineticsControl {
                reactivity: 0.0,
                external_source_per_s: -1.0
            },
            vec![]
        )
        .is_err()
    );
    let mut source = kinetic_source(0.1, false);
    let (idle, report) = source
        .trial_source(0.0, 0.0, 300.0, kinetic_budget(), &mut || Ok(()))
        .unwrap();
    assert_eq!(idle.kinetics, source.kinetics);
    assert_eq!(idle.irradiation.inventories, source.irradiation.inventories);
    assert_eq!(report.mean_neutron_population, 1.0);
    assert_eq!(report.kinetics.integrated_population_s, 0.0);
    assert_eq!(report.irradiation.fissions, 0.0);
    assert_eq!(report.mean_population_roundoff_s, 0.0);
    source.feedback.push(feedback(-0.001, 1000.0));
    let too_few = KineticHeatingBudget {
        max_feedback_terms: 0,
        ..kinetic_budget()
    };
    assert!(matches!(
        source.admit(0.0, 300.0, too_few),
        Err(FissionThermalError::Kinetics(KineticsError::Budget(_)))
    ));
    // Both authored inputs are finite; their physical flux mapping is not.
    let mut overflowing = source.clone();
    overflowing.kinetics =
        PointKineticsState::new(Arc::new(source.kinetics.model().clone()), 1e300, vec![]).unwrap();
    assert!(matches!(
        overflowing.admit(0.0, 300.0, kinetic_budget()),
        Err(FissionThermalError::Arithmetic)
    ));
}
