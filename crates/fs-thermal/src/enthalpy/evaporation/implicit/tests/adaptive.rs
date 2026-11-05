//! Local h/T/mass accuracy and conservative publication of the accepted fine path.
use super::*;

fn binding(id: u64, law: &Interface, film: GasFilm) -> EvaporatingCellBinding<'_, Interface> {
    EvaporatingCellBinding {
        id: CellId(id),
        law,
        film,
        pure_liquid_enthalpy_range_j_kg: range(0.0),
    }
}

fn no_sources(
    _: &DepositionSimulation<'_>,
    _: &BTreeMap<CellId, LiquidSurfaceState>,
) -> Result<BoundaryExchangeStepCoefficients<()>, ThermalEvaporationError> {
    Ok(coefficients(()))
}

fn policy(span: f64, temperature_tolerance: f64) -> AdaptiveEvaporationControls {
    let mut solve = EvaporationNetworkImplicitControls::default();
    solve.maximum_iterations = 10_000;
    solve.maximum_evaluations = 100_000;
    let mut controls = AdaptiveEvaporationControls::new(solve);
    controls.initial_dt_s = span;
    controls.minimum_dt_s = 1e-4;
    controls.maximum_dt_s = span;
    controls.maximum_comparisons = 512;
    controls.maximum_reductions = 16;
    controls.enthalpy_absolute_tolerance_j_kg = CP * temperature_tolerance;
    controls.temperature_absolute_tolerance_k = temperature_tolerance;
    controls.mass_absolute_tolerance_kg = 1e-12;
    controls.relative_tolerance = 0.0;
    controls
}

fn heated() -> DepositionSimulation<'static> {
    let mut liquid = cell(1, 0.1, 300.0, 0.0, 0.0);
    liquid.boundary.ambient_k = 310.0;
    liquid.boundary.htc_w_m2_k = 400.0;
    simulation(&CALORIC, vec![liquid], None, vec![])
}

fn unchanged(actual: &DepositionSimulation<'_>, expected: &DepositionSimulation<'_>) {
    assert_eq!(actual.time_s, expected.time_s);
    assert_eq!(actual.cursor, expected.cursor);
    assert_eq!(actual.cells, expected.cells);
    assert_eq!(actual.links, expected.links);
    assert_eq!(actual.covered_areas, expected.covered_areas);
    assert_eq!(actual.schedule, expected.schedule);
    assert_eq!(actual.budget, expected.budget);
}

#[test]
fn tighter_time_tolerances_improve_the_independent_variable_mass_solution() {
    let old = heated();
    let bindings = [binding(1, &CONSTANT, film(0.0))];
    let duration = 200.0;
    let rate = transfer(&CONSTANT, 300.0, film(0.0), 0.01, 1.0);
    let equilibrium = 310.0 - rate * LATENT / 4.0;
    // Exact m(t) Cp T' = G(Tbath-T) - r L, with m(t) = m0-r t.
    let exact =
        equilibrium + (300.0 - equilibrium) * (1.0 - rate * duration / 0.1).powf(4.0 / (CP * rate));
    let run = |tolerance| {
        let controls = policy(duration, tolerance);
        let (next, _, report) = old
            .advance_adaptive_evaporating_network_to(
                &bindings,
                duration,
                controls,
                no_sources,
                || Ok(()),
            )
            .unwrap();
        close(next.mass_kg(), 0.1 - rate * duration, 1e-14);
        close(report.physical.energy_residual_j, 0.0, 1e-6);
        assert!(report.absolute_energy_residual_j <= controls.solve.total_energy_tolerance_j);
        assert!(report.maximum_accepted_normalized_error <= 1.0);
        assert_eq!(report.physical.substeps, 2 * report.accepted_comparisons);
        (
            (next.temperature_k(CellId(1)).unwrap() - exact).abs(),
            report,
        )
    };
    let (loose_error, loose) = run(0.1);
    let (tight_error, tight) = run(0.001);
    assert!(loose_error > 0.01);
    assert!(
        tight_error < 0.25 * loose_error,
        "{loose_error:e}, {tight_error:e}"
    );
    assert!(tight.accepted_comparisons > loose.accepted_comparisons);
    assert!(tight.accuracy_rejections > 0);
    assert_eq!(old.time_s(), 0.0);
    assert_eq!(old.cells()[&CellId(1)].mass_kg, 0.1);
}

#[test]
fn equal_temperatures_cannot_hide_liquid_mass_or_plateau_enthalpy_error() {
    let interface = Interface {
        pressure_slope: 100.0,
        latent_slope: 0.0,
    };
    let gas = GasFilm {
        mass_transfer_coefficient_m_s: 10.0,
        ..film(0.0)
    };
    let r0 = transfer(&interface, 300.0, gas, 0.01, 1.0);
    let r1 = transfer(&interface, 301.0, gas, 0.01, 1.0);
    let feedback = (0.1 - r1) * CP / (r1 - r0);
    let pulse = 0.1 * CP + (r1 - r0) * (LATENT + feedback);
    let old = simulation(&CALORIC, vec![cell(1, 0.1, 300.0, 0.0, 0.0)], None, vec![]);
    let before = old.clone();
    let bindings = [binding(1, &interface, gas)];
    let mut controls = policy(2.0, 1e-6);
    controls.minimum_dt_s = 2.0;
    controls.enthalpy_absolute_tolerance_j_kg = 1e-3;
    controls.mass_absolute_tolerance_kg = 1e-8;
    controls.relative_tolerance = 0.01;
    // The coarse endpoint is 300 K. Fine endpoints are 301 K then 300 K,
    // but the hotter midpoint exports more liquid. Scaling by retained mass
    // instead of transferred mass would incorrectly admit this comparison.
    let result = old.advance_adaptive_evaporating_network_to(
        &bindings,
        2.0,
        controls,
        |thermal, liquids| {
            let time = thermal.time_s();
            let mass = liquids[&CellId(1)].mass_kg;
            let power =
                r0 * LATENT + pulse * time * (2.0 - time) + feedback * (mass - (0.1 - r0 * time));
            let mut result = coefficients(());
            result
                .transport
                .heat
                .bulk
                .sources_w
                .insert(CellId(1), power);
            Ok(result)
        },
        || Ok(()),
    );
    let error = match result {
        Err(AdaptiveEvaporationError::TimeAccuracy {
            time_s,
            span_s,
            error,
        }) => {
            assert_eq!(time_s, 0.0);
            assert_eq!(span_s, 2.0);
            error
        }
        Err(error) => panic!("wrong refusal: {error:?}"),
        Ok(_) => panic!("equal temperature hid a resolved mass discrepancy"),
    };
    assert!(error.maximum_temperature_difference_k < controls.temperature_absolute_tolerance_k);
    assert!(error.maximum_enthalpy_difference_j_kg < controls.enthalpy_absolute_tolerance_j_kg);
    close(error.maximum_mass_difference_kg, r1 - r0, 1e-11);
    let expected_ratio =
        (r1 - r0) / (controls.mass_absolute_tolerance_kg + controls.relative_tolerance * (r1 + r0));
    close(error.normalized, expected_ratio, 1e-6);
    assert!(error.normalized > 1.0);
    unchanged(&old, &before);

    // A separate heat-only latent plateau has no liquid/vapor binding. Its
    // original h must still be controlled even though T and mass stay fixed.
    struct Plateau;
    impl ThermalLaw for Plateau {
        fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
            if (-100.0..=200.0).contains(&h) {
                Ok(300.0)
            } else {
                Err(DepositionError::MaterialDomain)
            }
        }
        fn maximum_temperature_slope(&self) -> f64 {
            0.0
        }
    }
    let plateau = Plateau;
    let old = simulation(&plateau, vec![cell(1, 1.0, 300.0, 0.0, 0.0)], None, vec![]);
    let empty_bindings: [EvaporatingCellBinding<'_, Interface>; 0] = [];
    let mut controls = policy(1.0, 1e-6);
    controls.minimum_dt_s = 1.0;
    controls.enthalpy_absolute_tolerance_j_kg = 1.0;
    let result = old.advance_adaptive_evaporating_network_to(
        &empty_bindings,
        1.0,
        controls,
        |thermal, liquids| {
            assert!(liquids.is_empty());
            let mut result = coefficients(());
            result.transport.heat.bulk.sources_w.insert(
                CellId(1),
                10.0 * (100.0 - thermal.cells()[&CellId(1)].specific_h_j_kg),
            );
            Ok(result)
        },
        || Ok(()),
    );
    match result {
        Err(AdaptiveEvaporationError::TimeAccuracy { error, .. }) => {
            assert_eq!(error.maximum_temperature_difference_k, 0.0);
            assert_eq!(error.maximum_mass_difference_kg, 0.0);
            assert!(error.maximum_enthalpy_difference_j_kg > 6.0);
            assert!(error.normalized > 1.0);
        }
        Err(error) => panic!("wrong plateau refusal: {error:?}"),
        Ok(_) => panic!("a temperature-only controller admitted latent inventory error"),
    }
    assert_eq!(old.cells()[&CellId(1)].specific_h_j_kg, 0.0);
}

#[test]
fn rejected_comparisons_never_enter_surface_totals_or_endpoint_event_payloads() {
    let initial = cell(1, 0.1, 300.0, 0.0, 0.0);
    let mut covered = initial.boundary;
    covered.source_w = 4.0;
    covered.reservoir = Some(Reservoir {
        area_m2: 0.01,
        conductance_w_k: 0.0,
        temperature_k: 300.0,
    });
    let mut final_boundary = covered;
    final_boundary.source_w = 99.0;
    let old = simulation(
        &CALORIC,
        vec![initial],
        None,
        vec![
            TimedEvent {
                time_s: 500.0,
                event: Event::Boundary(CellId(1), covered),
            },
            TimedEvent {
                time_s: 1000.0,
                event: Event::Boundary(CellId(1), final_boundary),
            },
        ],
    );
    let bindings = [binding(1, &CONSTANT, film(0.0))];
    let mut calls = 0;
    let (next, payload, report) = old
        .advance_adaptive_evaporating_network_to(
            &bindings,
            1000.0,
            policy(500.0, 0.001),
            |thermal, liquids| {
                calls += 1;
                Ok(coefficients((
                    thermal.time_s(),
                    liquids[&CellId(1)].area_m2,
                    thermal.cells()[&CellId(1)].boundary.source_w,
                )))
            },
            || Ok(()),
        )
        .unwrap();
    let expected_mass = transfer(&CONSTANT, 300.0, film(0.0), 0.01, 500.0);
    assert!(report.accuracy_rejections > 0);
    assert_eq!(
        report.comparisons,
        report.accepted_comparisons + report.accuracy_rejections
    );
    assert_eq!(report.physical.substeps, 2 * report.accepted_comparisons);
    assert_eq!(calls, report.residual_evaluations);
    close(report.physical.transferred_mass_kg, expected_mass, 1e-14);
    close(
        report.physical.latent_energy_j,
        expected_mass * LATENT,
        1e-7,
    );
    close(report.physical.external_heat_j, 2000.0, 1e-7);
    close(next.mass_kg(), 0.1 - expected_mass, 1e-14);
    close(
        next.energy_j() - old.energy_j() + report.physical.vapor_energy_j,
        2000.0,
        1e-6,
    );
    let total = &report.surface_totals[&CellId(1)];
    close(total.transferred_mass_kg, expected_mass, 1e-14);
    close(total.latent_energy_j, expected_mass * LATENT, 1e-7);
    close(total.vapor_energy_j, report.physical.vapor_energy_j, 1e-8);
    close(total.recession_m, expected_mass / (0.01 * 1000.0), 1e-15);
    assert_eq!(total.intervals, report.physical.substeps);
    assert_eq!(report.physical.events, 2);
    assert_eq!(payload.unwrap(), (1000.0, 0.0, 4.0));
    assert_eq!(next.cells()[&CellId(1)].boundary, final_boundary);
    assert_eq!(report.last_surfaces[&CellId(1)].transferred_mass_kg, 0.0);
    assert_eq!(report.last_surfaces[&CellId(1)].liquid.area_m2, 0.0);
    assert_eq!(old.time_s(), 0.0);
    assert_eq!(old.cells()[&CellId(1)].mass_kg, 0.1);
}

#[test]
fn rejected_work_shares_all_budgets_and_final_cancellation_is_atomic() {
    let old = heated();
    let before = old.clone();
    let bindings = [binding(1, &CONSTANT, film(0.0))];
    let (_, _, comparison) = old
        .advance_adaptive_evaporating_network_to(
            &bindings,
            200.0,
            policy(200.0, 10.0),
            no_sources,
            || Ok(()),
        )
        .unwrap();
    assert_eq!(comparison.comparisons, 1);
    assert_eq!(comparison.accepted_comparisons, 1);
    assert!(comparison.newton_iterations > 0);
    for limit in 0..3 {
        let mut controls = policy(200.0, 0.001);
        match limit {
            0 => controls.solve.maximum_evaluations = comparison.residual_evaluations,
            1 => controls.solve.maximum_iterations = comparison.newton_iterations,
            _ => controls.maximum_comparisons = 1,
        }
        let mut calls = 0;
        let result = old.advance_adaptive_evaporating_network_to(
            &bindings,
            200.0,
            controls,
            |_, _| {
                calls += 1;
                Ok(coefficients(()))
            },
            || Ok(()),
        );
        if limit == 2 {
            assert!(matches!(
                result,
                Err(AdaptiveEvaporationError::ComparisonBudget)
            ));
        } else {
            assert!(matches!(
                result,
                Err(AdaptiveEvaporationError::Implicit(
                    EvaporationImplicitError::Physical(ThermalEvaporationError::Thermal(
                        DepositionError::Budget(_)
                    ))
                ))
            ));
        }
        assert_eq!(calls, comparison.residual_evaluations);
        unchanged(&old, &before);
    }
    let controls = policy(200.0, 0.05);
    let mut checkpoints = 0;
    let (_, _, completed) = old
        .advance_adaptive_evaporating_network_to(&bindings, 200.0, controls, no_sources, || {
            checkpoints += 1;
            Ok(())
        })
        .unwrap();
    assert!(completed.accuracy_rejections > 0);
    assert!(completed.accepted_comparisons > 1);
    let mut visited = 0;
    let cancelled =
        old.advance_adaptive_evaporating_network_to(&bindings, 200.0, controls, no_sources, || {
            visited += 1;
            if visited == checkpoints {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        });
    assert!(matches!(
        cancelled,
        Err(AdaptiveEvaporationError::Implicit(
            EvaporationImplicitError::Physical(ThermalEvaporationError::Thermal(
                DepositionError::Cancelled
            ))
        ))
    ));
    assert_eq!(visited, checkpoints);
    unchanged(&old, &before);
}

#[test]
fn evaporating_and_condensing_surfaces_keep_individual_totals_in_canonical_order() {
    let dry = film(0.0);
    // Stefan logarithms have equal magnitude and opposite sign: 25/24 and 24/25.
    let humid = film(7840.0);
    let exposed = 0.008;
    let rate_a = transfer(&CONSTANT, 301.0, dry, exposed, 1.0);
    let rate_b = transfer(&CONSTANT, 299.0, humid, exposed, 1.0);
    let old = simulation(
        &CALORIC,
        vec![
            cell(1, 0.1, 301.0, 0.0, rate_a * LATENT + 4.0),
            cell(2, 0.1, 299.0, 0.0, rate_b * LATENT - 4.0),
        ],
        Some(2.0),
        vec![],
    );
    let forward = [binding(1, &CONSTANT, dry), binding(2, &CONSTANT, humid)];
    let reverse = [binding(2, &CONSTANT, humid), binding(1, &CONSTANT, dry)];
    let controls = policy(100.0, 1e-3);
    let (next, _, report) = old
        .advance_adaptive_evaporating_network_to(&forward, 100.0, controls, no_sources, || Ok(()))
        .unwrap();
    let (reordered, _, reordered_report) = old
        .advance_adaptive_evaporating_network_to(&reverse, 100.0, controls, no_sources, || Ok(()))
        .unwrap();
    assert_eq!(next.cells(), reordered.cells());
    assert_eq!(report, reordered_report);
    assert_eq!(report.accepted_comparisons, 1);
    assert_eq!(report.physical.substeps, 2);
    for (id, rate, temperature) in [(CellId(1), rate_a, 301.0), (CellId(2), rate_b, 299.0)] {
        let total = &report.surface_totals[&id];
        close(next.temperature_k(id).unwrap(), temperature, 1e-10);
        close(next.cells()[&id].mass_kg, 0.1 - rate * 100.0, 1e-15);
        close(total.transferred_mass_kg, rate * 100.0, 1e-15);
        close(total.latent_energy_j, rate * 100.0 * LATENT, 1e-8);
        close(
            total.vapor_energy_j,
            rate * 100.0 * (LATENT + CP * (temperature - 300.0)),
            1e-8,
        );
        close(total.recession_m, rate * 100.0 / (exposed * 1000.0), 1e-15);
        assert_eq!(total.intervals, 2);
        assert_eq!(
            report.last_surfaces[&id].liquid.mass_kg,
            next.cells()[&id].mass_kg
        );
    }
    assert!(report.surface_totals[&CellId(1)].transferred_mass_kg > 1e-4);
    assert!(report.surface_totals[&CellId(2)].transferred_mass_kg < -1e-4);
    close(report.physical.transferred_mass_kg, 0.0, 1e-14);
    close(report.physical.external_heat_j, 0.0, 1e-8);
    close(
        next.energy_j() - old.energy_j() + report.physical.vapor_energy_j,
        0.0,
        1e-7,
    );
    close(report.physical.energy_residual_j, 0.0, 1e-7);
}
