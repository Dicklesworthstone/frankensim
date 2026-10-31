use super::*;
use fs_render::volumes::heating::exchange::refinement::{
    RadiativeExchangeRefinement, RadiativeExchangeRefinementError,
};

fn controls() -> RadiativeExchangeRefinement {
    RadiativeExchangeRefinement {
        initial_max_dt_s: 0.5,
        minimum_max_dt_s: 1e-5,
        maximum_halvings: 16,
        maximum_total_substeps: 20_000,
        maximum_total_source_evaluations: 20_000,
        enthalpy_absolute_tolerance_j_kg: 0.4,
        temperature_absolute_tolerance_k: 0.04,
        radiative_energy_absolute_tolerance_j: 0.4,
        relative_tolerance: 0.0,
        conservation_tolerance_j: 1e-8,
    }
}

// These success studies request two resolved .04 K/.4 J comparisons. Their
// explicit native/source budgets must admit the fine trajectories as well as
// the study's separate aggregate cap. Ordinary core fixtures keep their 1000.
fn refinement_thermal(
    cells: Vec<DepositedCell>,
    law: &'static dyn ThermalLaw,
    events: Vec<TimedEvent>,
) -> DepositionSimulation<'static> {
    let mut schedule: Vec<_> = cells
        .into_iter()
        .map(|cell| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(cell),
        })
        .collect();
    schedule.extend(events);
    DepositionSimulation::new(
        schedule,
        &[(MaterialId(7), law)],
        DepositionBudget {
            max_cells: 8,
            max_links: 8,
            max_events: 16,
            max_substeps: 8192,
        },
    )
    .unwrap()
    .advance_to(0.0, config(0.1), || Ok(()))
    .unwrap()
    .0
}

fn refinement_observer(
    op: RadiativeExchangeOperator<1>,
    weight_nm: f64,
    bath: RadiativeBath,
    bindings: &[RadiativeThermalCell],
    thermal: &DepositionSimulation<'_>,
) -> VolumeRadiativeExchange<1> {
    VolumeRadiativeExchange::try_new(
        op,
        [weight_nm],
        bath,
        bindings,
        thermal,
        RadiativeExchangeBudget {
            max_source_evaluations: 8192,
            ..budget()
        },
        tolerance(),
        || Ok(()),
    )
    .unwrap()
}

#[test]
fn nonlinear_cooling_refines_against_independent_fine_rk4_planck_solution() {
    with_cx(|cx| {
        let original = refinement_thermal(vec![warm_cell(11, 500.0)], &RADIATIVE_SENSIBLE, vec![]);
        let op = operator(cx, &[1.0], 0.7, 1.0, photon_config());
        let weight = 3000.0;
        let scale = op.cells[0].escape.mean_m2_sr[0] * weight;
        let exchange = refinement_observer(
            op,
            weight,
            RadiativeBath::Vacuum,
            &[binding(0, 11)],
            &original,
        );
        // Independent ODE owner for this test only: constant 10 J/(kg K),
        // 1 kg, direct SI Planck, and the prescribed sampled optical operator.
        // This verifies time integration conditional on that operator, not its
        // Monte Carlo coefficient uncertainty or spectral quadrature accuracy.
        let target = 1.0;
        let dt = target / 20_000.0;
        let rhs = |temperature| -scale * planck(temperature) / 10.0;
        let mut reference = 500.0;
        for _ in 0..20_000 {
            let a = rhs(reference);
            let b = rhs(reference + 0.5 * dt * a);
            let c = rhs(reference + 0.5 * dt * b);
            let d = rhs(reference + dt * c);
            reference += dt * (a + 2.0 * b + 2.0 * c + d) / 6.0;
        }
        let (cooled, accepted, report) = exchange
            .advance_refined(&original, target, controls(), || Ok(()))
            .unwrap();
        assert!((cooled.temperature_k(CellId(11)).unwrap() - reference).abs() < 0.04);
        assert!(report.halvings >= 2);
        assert!(report.normalized_difference <= 1.0);
        assert!(report.previous_normalized_difference <= 1.0);
        assert!(report.thermal.substeps > report.penultimate_substeps);
        assert!(report.total_substeps > report.thermal.substeps);
        assert_eq!(report.total_substeps, report.total_source_evaluations);
        assert_eq!(accepted.source_evaluations(), report.thermal.substeps);
        assert!(!report.roundoff_limited);
        assert!(report.observed_order.is_some_and(|p| p > 0.8));
        relative_close(
            report.thermal.external_heat_j,
            accepted.energy().thermal_external_j,
        );
        assert_eq!(original.time_s(), 0.0);
        assert_eq!(exchange.time_s(), 0.0);
    });
}

#[test]
fn unchanged_native_cfl_grids_cannot_claim_time_refinement() {
    with_cx(|cx| {
        let original = thermal(vec![warm_cell(11, 500.0)], &RADIATIVE_SENSIBLE, vec![]);
        let op = operator(cx, &[1.0], 0.7, 1.0, photon_config());
        let exchange = observer(
            op,
            3000.0,
            RadiativeBath::Vacuum,
            &[binding(0, 11)],
            &original,
        );
        let mut controls = controls();
        controls.initial_max_dt_s = 1e30;
        controls.minimum_max_dt_s = 1e20;
        controls.maximum_halvings = 3;
        controls.enthalpy_absolute_tolerance_j_kg = 1e9;
        controls.temperature_absolute_tolerance_k = 1e9;
        controls.radiative_energy_absolute_tolerance_j = 1e9;
        // All these large maxima collapse to the same radiative CFL steps.
        // Zero endpoint disagreement between identical executions is not a
        // time-resolution comparison, even with arbitrarily loose tolerances.
        assert!(matches!(
            exchange.advance_refined(&original, 20.0, controls, || Ok(())),
            Err(RadiativeExchangeRefinementError::Convergence { .. })
        ));
        assert_eq!(original.time_s(), 0.0);
        assert_eq!(exchange.energy().thermal_external_j, 0.0);
    });
}

#[test]
fn refinement_charges_all_trials_and_rolls_back_budgets_and_cancellation() {
    with_cx(|cx| {
        let original = thermal(vec![warm_cell(11, 500.0)], &RADIATIVE_SENSIBLE, vec![]);
        let op = operator(cx, &[1.0], 0.7, 1.0, photon_config());
        let exchange = observer(
            op,
            3000.0,
            RadiativeBath::Vacuum,
            &[binding(0, 11)],
            &original,
        );
        let original_cells = original.cells().clone();
        let mut steps = controls();
        steps.maximum_total_substeps = 1;
        assert!(
            exchange
                .advance_refined(&original, 1.0, steps, || Ok(()))
                .is_err()
        );
        let mut sources = controls();
        sources.maximum_total_source_evaluations = 1;
        assert!(
            exchange
                .advance_refined(&original, 1.0, sources, || Ok(()))
                .is_err()
        );
        let mut calls = 0;
        let cancelled = exchange.advance_refined(&original, 1.0, controls(), || {
            calls += 1;
            if calls == 20 {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        });
        assert!(matches!(
            cancelled,
            Err(RadiativeExchangeRefinementError::Exchange(
                RadiativeExchangeError::Thermal(DepositionError::Cancelled)
            ))
        ));
        assert_eq!(calls, 20);
        assert_eq!(original.cells(), &original_cells);
        assert_eq!(original.time_s(), 0.0);
        assert_eq!(exchange.time_s(), 0.0);
        assert_eq!(exchange.source_evaluations(), 0);
        assert_eq!(exchange.energy().thermal_external_j, 0.0);
    });
}

#[test]
fn refined_two_cell_reabsorption_and_bath_energy_keep_native_event_accounting() {
    with_cx(|cx| {
        let original = refinement_thermal(
            vec![warm_cell(11, 450.0), warm_cell(12, 300.0)],
            &RADIATIVE_SENSIBLE,
            vec![TimedEvent {
                time_s: 0.5,
                event: Event::Boundary(CellId(12), boundary(10.0)),
            }],
        );
        let op = operator(cx, &[1.0; 2], 0.7, 1.0, photon_config());
        let exchange = refinement_observer(
            op,
            3000.0,
            blackbody(300.0),
            &[binding(0, 11), binding(1, 12)],
            &original,
        );
        let (accepted, radiation, report) = exchange
            .advance_refined(&original, 1.0, controls(), || Ok(()))
            .unwrap();
        assert_eq!(accepted.time_s(), 1.0);
        assert_eq!(report.thermal.events, 1);
        assert!(accepted.temperature_k(CellId(11)).unwrap() < 450.0);
        assert!(accepted.temperature_k(CellId(12)).unwrap() > 300.0);
        assert!(radiation.energy().internal_transfer_magnitude_j[0] > 0.0);
        assert!(radiation.energy().thermal_external_j < 0.0);
        relative_close(
            report.thermal.external_heat_j,
            radiation.energy().thermal_external_j + 5.0,
        );
        let storage: f64 = accepted
            .cells()
            .iter()
            .map(|(id, c)| c.mass_kg * (c.specific_h_j_kg - original.cells()[id].specific_h_j_kg))
            .sum();
        relative_close(storage, report.thermal.external_heat_j);
        assert!(report.thermal.substeps > report.penultimate_substeps);
        assert!(
            report.normalized_difference <= 1.0 && report.previous_normalized_difference <= 1.0
        );
        // Refinement work caps do not permanently reduce a returned observer's
        // ordinary native-advance source allowance.
        assert_eq!(
            radiation.budget().max_source_evaluations,
            exchange.budget().max_source_evaluations
        );
    });
}

struct FormationDatum;
impl ThermalLaw for FormationDatum {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        let relative = h - 1e12;
        if !relative.is_finite() || !(-2000.0..=5000.0).contains(&relative) {
            return Err(DepositionError::MaterialDomain);
        }
        Ok(300.0 + relative / 10.0)
    }
    fn maximum_temperature_slope(&self) -> f64 {
        0.1
    }
}
static FORMATION: FormationDatum = FormationDatum;

#[test]
fn exact_equal_temperatures_do_not_certify_sub_ulp_formation_energy_accuracy() {
    with_cx(|cx| {
        let mut initial = cell(11, 1.0, 1.0, 0.0);
        initial.specific_h_j_kg = 1e12;
        let original = thermal(vec![initial], &FORMATION, vec![]);
        let op = operator(cx, &[1.0], 0.7, 1.0, photon_config());
        let exchange = observer(op, 3000.0, blackbody(300.0), &[binding(0, 11)], &original);
        let mut controls = controls();
        controls.maximum_halvings = 2;
        controls.minimum_max_dt_s = controls.initial_max_dt_s * 0.25;
        controls.enthalpy_absolute_tolerance_j_kg = 1e-6;
        controls.temperature_absolute_tolerance_k = 1e-6;
        controls.radiative_energy_absolute_tolerance_j = 1e-6;
        // Heat is exactly zero, but native h=1e12 has a resolution much larger
        // than 1e-6 J/kg. Neither zero discrepancy nor a large datum may turn
        // this into a false claim of resolved sub-ULP accuracy.
        assert!(matches!(
            exchange.advance_refined(&original, 1.0, controls, || Ok(())),
            Err(RadiativeExchangeRefinementError::Convergence { .. })
        ));
        assert_eq!(original.cells()[&CellId(11)].specific_h_j_kg, 1e12);
        assert_eq!(exchange.energy().thermal_external_j, 0.0);
    });
}
