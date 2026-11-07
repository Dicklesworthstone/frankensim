//! G1/G3/G4: physical switching time, thermal holds and atomic bounded refinement.
use super::*;
use fs_couple::circuit_enthalpy::dielectric::hysteretic::trajectory::HystereticThermalPathControls;

fn controls() -> HystereticThermalPathControls {
    HystereticThermalPathControls {
        thermal: heat_controls(),
        path: HystereticPathControls {
            field: field_controls(),
            relative_tolerance: 1e-5,
            absolute_polarization_tolerance_c_m2: 1e-8,
            absolute_heat_tolerance_j: 1e-8,
            absolute_work_tolerance_j: 0.02,
            energy_tolerance_j: 2e-8,
            maximum_field_solves: 8192,
            maximum_operator_applications: 2_000_000,
            maximum_refinement_depth: 20,
        },
        absolute_thermal_tolerance_j: 0.002,
        relative_thermal_tolerance: 1e-4,
        maximum_thermal_evaluations: 100_000,
    }
}

fn cooling(curve: &ResolvedHeatCapacityCurve, temperature: f64) -> DepositionSimulation<'_> {
    let seed = thermal(curve, 1, false, temperature);
    let mut cell = seed.cells()[&CellId(0)].clone();
    cell.boundary.ambient_k = 300.0;
    cell.boundary.htc_w_m2_k = CAPACITY; // area=1 m²: time constant C/G=1 second.
    DepositionSimulation::new(
        vec![TimedEvent {
            time_s: 0.0,
            event: Event::Birth(cell),
        }],
        &[(MaterialId(0), curve.curve() as &dyn ThermalLaw)],
        seed.budget(),
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 1.0,
            energy_tolerance_j: 1e-9,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}

#[test]
fn switching_heat_is_released_at_the_switching_time_during_cooling() {
    let card = card(true);
    let curve = caloric(&card);
    let context = point(300.0);
    let bindings = cells(&card, &curve, &context, 1, None);
    let network = MaterialNetwork {
        cells: &bindings,
        resistors: &[],
        links: &[],
    };
    let thermal = cooling(&curve, 300.0);
    let coupled = ThermalHystereticElectrostatics::bind(
        field(&card, true, 300.0, [0.6, 0.8]),
        &thermal,
        &network,
        &regions(),
        1e-12,
        || false,
    )
    .unwrap();
    let smeared = coupled
        .step_candidate(
            &thermal,
            1.0,
            &[0.0; 4],
            &voltage(12.0),
            heat_controls(),
            controls().path,
            &[0.0],
            || false,
        )
        .unwrap();
    let mut loose_controls = controls();
    loose_controls.absolute_thermal_tolerance_j *= 10.0;
    loose_controls.path.absolute_work_tolerance_j *= 10.0;
    let loose = coupled
        .advance_affine_candidate(
            &thermal,
            1.0,
            &[0.0; 4],
            &voltage(12.0),
            loose_controls,
            || false,
        )
        .unwrap();
    let fine = coupled
        .advance_affine_candidate(&thermal, 1.0, &[0.0; 4], &voltage(12.0), controls(), || {
            false
        })
        .unwrap();
    // E(t)=12t. Independent scalar play releases constant heat only while
    // Ec<E<Ec+Ps/chi. Solve C*theta'=Qdot-G*theta over that active interval,
    // then cool with Qdot=0 to t=1. The total physical Q is still Ec*volume*Ps.
    let onset = EC / 12.0;
    let saturation = (EC + PS / CHI) / 12.0;
    let power = DEPTH * EC * CHI * 12.0;
    let exact = 300.0
        + power / CAPACITY * (1.0 - (-(saturation - onset)).exp()) * (-(1.0 - saturation)).exp();
    let temperature = fine.thermal().temperature_k(CellId(0)).unwrap();
    let loose_error = (loose.thermal().temperature_k(CellId(0)).unwrap() - exact).abs();
    let error = (temperature - exact).abs();
    assert!(
        error < 0.4 * loose_error,
        "fine={error:e}, coarse={loose_error:e}"
    );
    assert!(error < 1e-3, "temperature error {error:e}");
    assert!((smeared.thermal().temperature_k(CellId(0)).unwrap() - exact).abs() > 5.0 * error);
    near(fine.report().deposited_joule_j, DEPTH * EC * PS);
    near(fine.report().electrical.switching_heat_j, DEPTH * EC * PS);
    near(
        fine.report().thermal.stored_energy_change_j,
        CAPACITY * (temperature - 300.0),
    );
    near(fine.report().total_energy_residual_j, 0.0);
    assert!(fine.report().non_electrical_heat_j < 0.0);
    assert!(fine.report().electrical.numerical_dissipation_j > 0.0);
    assert!(fine.summary().accepted_steps > loose.summary().accepted_steps);
    assert!(fine.summary().rejected_intervals > 0);
    assert_eq!(
        fine.report().circuit_trials,
        fine.summary().thermal_evaluations
    );
    assert_eq!(
        fine.report().thermal.substeps,
        fine.summary().accepted_steps
    );
    assert_eq!(fine.field().time_s(), 1.0);
    assert_eq!(fine.thermal().time_s(), 1.0);
    let DielectricMaterial::Hysteretic(source) = &fine.field().field().problem().materials()[0]
    else {
        panic!("lost original source")
    };
    near(source.source().material().query().axes()["T"], temperature);
    let summary = *fine.summary();
    let (_, _, report, retained) = fine.into_parts();
    assert_eq!(summary, retained);
    near(report.deposited_joule_j, DEPTH * EC * PS);
}

#[test]
fn a_zero_field_hold_still_refines_the_original_heat_assembly_timeline() {
    let card = card(true);
    let curve = caloric(&card);
    let context = point(320.0);
    let bindings = cells(&card, &curve, &context, 1, None);
    let network = MaterialNetwork {
        cells: &bindings,
        resistors: &[],
        links: &[],
    };
    let thermal = cooling(&curve, 320.0);
    let coupled = ThermalHystereticElectrostatics::bind(
        field(&card, true, 320.0, [1.0, 0.0]),
        &thermal,
        &network,
        &regions(),
        1e-12,
        || false,
    )
    .unwrap();
    let electrical = ElectricalContactNetwork {
        bulk: &network,
        contacts: &[],
        electrical_bases: &[],
    };
    let heat = MaterialHeatNetwork {
        electrical: &electrical,
        thermal_contacts: &[],
        radiation: None,
    };
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    let (coarse, fine) = pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            fs_exec::StreamKey {
                seed: 17,
                kernel_id: 1,
                tile: 0,
                iteration: 0,
            },
            fs_exec::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        );
        let prepared = heat.prepare(&cx, &thermal).unwrap();
        let mut coarse_controls = controls();
        coarse_controls.absolute_thermal_tolerance_j = 0.08;
        coarse_controls.relative_thermal_tolerance = 0.01;
        let coarse = coupled
            .advance_affine_candidate_with_heat_assembly(
                &cx,
                &prepared,
                &thermal,
                1.0,
                &[0.0; 4],
                &voltage(0.0),
                coarse_controls,
            )
            .unwrap();
        let mut fine_controls = coarse_controls;
        fine_controls.absolute_thermal_tolerance_j = 0.02;
        fine_controls.relative_thermal_tolerance = 0.001;
        let fine = coupled
            .advance_affine_candidate_with_heat_assembly(
                &cx,
                &prepared,
                &thermal,
                1.0,
                &[0.0; 4],
                &voltage(0.0),
                fine_controls,
            )
            .unwrap();
        (coarse, fine)
    });
    let exact = 300.0 + 20.0 * (-1.0f64).exp();
    let coarse_error = (coarse.thermal().temperature_k(CellId(0)).unwrap() - exact).abs();
    let error = (fine.thermal().temperature_k(CellId(0)).unwrap() - exact).abs();
    assert!(
        error < 0.25 * coarse_error,
        "fine={error:e}, coarse={coarse_error:e}"
    );
    assert!(fine.summary().accepted_steps > coarse.summary().accepted_steps);
    assert!(fine.summary().estimated_thermal_energy_error_j > 0.0);
    assert_eq!(fine.summary().estimated_polarization_error_c_m2, 0.0);
    assert_eq!(fine.summary().estimated_switching_heat_error_j, 0.0);
    assert_eq!(fine.summary().estimated_source_work_error_j, 0.0);
    assert_eq!(fine.report().deposited_joule_j, 0.0);
    assert_eq!(fine.report().electrical.numerical_dissipation_j, 0.0);
    near(
        fine.report().thermal.stored_energy_change_j,
        fine.report().non_electrical_heat_j,
    );
    near(fine.report().total_energy_residual_j, 0.0);
}

#[test]
fn chronology_work_caps_and_late_cancel_preserve_both_histories_for_exact_retry() {
    let card = card(true);
    let curve = caloric(&card);
    let context = point(300.0);
    let bindings = cells(&card, &curve, &context, 1, None);
    let network = MaterialNetwork {
        cells: &bindings,
        resistors: &[],
        links: &[],
    };
    let thermal = cooling(&curve, 300.0);
    let coupled = ThermalHystereticElectrostatics::bind(
        field(&card, true, 300.0, [1.0, 0.0]),
        &thermal,
        &network,
        &regions(),
        1e-12,
        || false,
    )
    .unwrap();
    let mut settings = controls();
    settings.path.absolute_work_tolerance_j = 0.2;
    settings.absolute_thermal_tolerance_j = 0.02;
    let mut polls = 0;
    let expected = coupled
        .advance_affine_candidate(&thermal, 1.0, &[0.0; 4], &voltage(12.0), settings, || {
            polls += 1;
            false
        })
        .unwrap();
    let mut stopped = 0;
    assert!(matches!(
        coupled.advance_affine_candidate(
            &thermal,
            1.0,
            &[0.0; 4],
            &voltage(12.0),
            settings,
            || {
                stopped += 1;
                stopped >= polls - 1
            }
        ),
        Err(CircuitEnthalpyError::Cancelled)
    ));
    let mut limited = settings;
    limited.path.maximum_field_solves = 2;
    assert!(matches!(
        coupled
            .advance_affine_candidate(&thermal, 1.0, &[0.0; 4], &voltage(12.0), limited, || false),
        Err(CircuitEnthalpyError::Source(
            HystereticDielectricThermalError::Path(HystereticPathError::FieldSolveBudget)
        ))
    ));
    limited = settings;
    limited.maximum_thermal_evaluations = 1;
    assert!(matches!(
        coupled
            .advance_affine_candidate(&thermal, 1.0, &[0.0; 4], &voltage(12.0), limited, || false),
        Err(CircuitEnthalpyError::CircuitBudget)
    ));
    limited = settings;
    limited.path.maximum_operator_applications = 1;
    assert!(matches!(
        coupled
            .advance_affine_candidate(&thermal, 1.0, &[0.0; 4], &voltage(12.0), limited, || false),
        Err(CircuitEnthalpyError::Source(
            HystereticDielectricThermalError::Path(HystereticPathError::OperatorBudget)
        ))
    ));
    limited = settings;
    limited.path.maximum_refinement_depth = 0;
    assert!(matches!(
        coupled
            .advance_affine_candidate(&thermal, 1.0, &[0.0; 4], &voltage(12.0), limited, || false),
        Err(CircuitEnthalpyError::Source(
            HystereticDielectricThermalError::Path(HystereticPathError::RefinementLimit { .. })
        ))
    ));
    let retried = coupled
        .advance_affine_candidate(&thermal, 1.0, &[0.0; 4], &voltage(12.0), settings, || false)
        .unwrap();
    assert_eq!(retried.thermal().cells(), expected.thermal().cells());
    assert_eq!(
        retried.report().electrical.cells,
        expected.report().electrical.cells
    );
    assert_eq!(retried.summary(), expected.summary());
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(thermal.temperature_k(CellId(0)).unwrap(), 300.0);
    assert_eq!(coupled.time_s(), 0.0);
    assert_eq!(
        coupled
            .field()
            .cell_history(0)
            .unwrap()
            .polarizations_c_m2(),
        &[[0.0; 3]]
    );
    assert!(retried.summary().attempted_steps <= settings.path.maximum_field_solves);
    assert!(
        retried.summary().field_operator_applications
            <= settings.path.maximum_operator_applications
    );
    assert!(retried.summary().thermal_evaluations <= settings.maximum_thermal_evaluations);
}
