//! Local constitutive feedback uses the same reciprocal/native thermal owners.
#[path = "feedback/refinement.rs"]
mod refinement;

#[path = "feedback/combined.rs"]
mod combined;

use super::*;
use fs_render::volumes::heating::exchange::feedback::{
    RadiativeFeedbackBudget, RadiativeFeedbackError, TemperatureRadiativeExchange,
};
use fs_render::volumes::heating::feedback::TemperatureAbsorptionTable;
use fs_render::volumes::heating::temperature_field::{
    TemperatureOpticalField, TemperatureOpticalFieldBudget, TemperatureScatteringTable,
    ThermalScatteringCell,
};
use fs_thermal::enthalpy::DepositionObserver;
use std::sync::Arc;

fn absorption(rows: &[(f64, f64)]) -> TemperatureAbsorptionTable {
    let values: Vec<_> = rows
        .iter()
        .map(|&(t, value)| {
            (
                t,
                MassAbsorptionSpectrum::source_table(
                    &[(8000.0, value), (18000.0, value)],
                    ContentHash([91; 32]),
                )
                .unwrap(),
            )
        })
        .collect();
    TemperatureAbsorptionTable::try_new(ContentHash([92; 32]), &values, 16).unwrap()
}
fn scattering(rows: &[(f64, f64)]) -> TemperatureScatteringTable {
    let values: Vec<_> = rows
        .iter()
        .map(|&(t, value)| {
            (
                t,
                MassScatteringSpectrum::source_table(
                    &[(8000.0, value), (18000.0, value)],
                    ContentHash([93; 32]),
                )
                .unwrap(),
            )
        })
        .collect();
    TemperatureScatteringTable::try_new(ContentHash([94; 32]), &values, 16).unwrap()
}
fn field<'g>(
    cx: &Cx<'_>,
    density: &'g [f64],
    thermal: &DepositionSimulation<'_>,
    absorption: &[TemperatureAbsorptionTable],
    scattering: &[TemperatureScatteringTable],
    ids: &[u64],
) -> Arc<TemperatureOpticalField<'g, 1>> {
    assert_eq!(density.len(), ids.len());
    let bindings: Vec<_> = ids
        .iter()
        .enumerate()
        .map(|(index, &id)| ThermalScatteringCell {
            grid_cell: index,
            thermal_cell: CellId(id),
            absorption_table: index.min(absorption.len() - 1),
            scattering_table: index.min(scattering.len() - 1),
            phase: 0,
        })
        .collect();
    Arc::new(
        TemperatureOpticalField::try_new(
            grid(cx, density, MassAbsorptionSpectrum::constant(0.0).unwrap()),
            [10000.0],
            absorption,
            scattering,
            &[PhaseFunction::Isotropic],
            &bindings,
            thermal,
            TemperatureOpticalFieldBudget {
                max_grid_cells: 8,
                max_tables: 16,
                max_total_temperature_knots: 128,
            },
            || Ok(()),
        )
        .unwrap(),
    )
}
fn budget() -> RadiativeFeedbackBudget {
    RadiativeFeedbackBudget {
        max_source_evaluations: 1000,
        max_traced_segments: 20_000_000,
        max_thermal_pairs: 16,
    }
}
fn photon_config() -> RadiativeExchangeConfig {
    RadiativeExchangeConfig {
        histories_per_emitter_lane: 128,
        ..super::photon_config()
    }
}
fn observer<'g>(
    field: Arc<TemperatureOpticalField<'g, 1>>,
    thermal: &DepositionSimulation<'_>,
    weight: f64,
    bath: RadiativeBath,
) -> TemperatureRadiativeExchange<'g, 1> {
    TemperatureRadiativeExchange::try_new(
        field,
        [weight],
        bath,
        thermal,
        budget(),
        photon_config(),
        tolerance(),
        || Ok(()),
    )
    .unwrap()
}

#[test]
fn constant_optical_packet_reuses_transport_and_matches_original_fixed_owner() {
    with_cx(|cx| {
        let original = thermal(vec![warm_cell(11, 500.0)], &RADIATIVE_SENSIBLE, vec![]);
        let field = field(
            cx,
            &[1.0],
            &original,
            &[absorption(&[(100.0, 0.7), (1000.0, 0.7)])],
            &[scattering(&[(100.0, 1.0), (1000.0, 1.0)])],
            &[11],
        );
        let exchange = observer(field, &original, 2000.0, RadiativeBath::Vacuum);
        assert!(exchange.operator().is_none());
        assert_eq!(exchange.traced_segments(), 0);
        assert_eq!(exchange.source_evaluations(), 0);
        let (cooled, dynamic, audit) = original
            .advance_observed(&exchange, 0.4, config(0.05), || Ok(()))
            .unwrap();
        assert!(cooled.temperature_k(CellId(11)).unwrap() < 500.0);
        assert_eq!(dynamic.operator_builds(), 1);
        assert_eq!(dynamic.source_evaluations(), audit.substeps);
        assert_eq!(
            dynamic.traced_segments(),
            dynamic.operator().unwrap().traced_segments
        );
        let frozen = super::observer(
            dynamic.operator().unwrap().clone(),
            2000.0,
            RadiativeBath::Vacuum,
            &[super::binding(0, 11)],
            &original,
        );
        let (reference, fixed, fixed_audit) = original
            .advance_observed(&frozen, 0.4, config(0.05), || Ok(()))
            .unwrap();
        assert_eq!(cooled.cells(), reference.cells());
        assert_eq!(dynamic.energy(), fixed.energy());
        assert_eq!(audit, fixed_audit);
        assert!(dynamic.power().unwrap().thermal_external_w > frozen.power().thermal_external_w);
    });
}

#[test]
fn temperature_laws_rebuild_actual_emission_capacity_and_conservative_heat() {
    with_cx(|cx| {
        let original = thermal(
            vec![warm_cell(11, 500.0), warm_cell(12, 300.0)],
            &RADIATIVE_SENSIBLE,
            vec![],
        );
        let laws = [
            absorption(&[(200.0, 0.2), (800.0, 2.0)]),
            absorption(&[(200.0, 1.0), (800.0, 0.4)]),
        ];
        let field = field(
            cx,
            &[1.0, 0.5],
            &original,
            &laws,
            &[scattering(&[(200.0, 0.1), (800.0, 1.2)])],
            &[11, 12],
        );
        let exchange = observer(field, &original, 500.0, blackbody(300.0));
        let mut sampled = exchange.clone();
        sampled.sources(&original, &mut || Ok(())).unwrap();
        let op = sampled.operator().unwrap();
        relative_close(
            op.cells[0].emission_capacity_m2_sr[0],
            4.0 * std::f64::consts::PI * 1.1,
        );
        relative_close(
            op.cells[1].emission_capacity_m2_sr[0],
            4.0 * std::f64::consts::PI * 0.5 * 0.9,
        );
        let (advanced, changed, audit) = original
            .advance_observed(&exchange, 0.1, config(0.025), || Ok(()))
            .unwrap();
        assert!(changed.operator_builds() > 1);
        assert_eq!(changed.operator_builds(), changed.source_evaluations());
        assert!(advanced.temperature_k(CellId(11)).unwrap() < 500.0);
        assert!(advanced.temperature_k(CellId(12)).unwrap() > 300.0);
        assert!(
            changed.operator().unwrap().cells[0].emission_capacity_m2_sr[0]
                < op.cells[0].emission_capacity_m2_sr[0]
        );
        assert_eq!(
            changed.sampling_config().seed,
            sampled.sampling_config().seed
        );
        assert_eq!(
            changed.sampling_config().source_id,
            sampled.sampling_config().source_id
        );
        relative_close(audit.external_heat_j, changed.energy().thermal_external_j);
        relative_close(
            audit.stored_energy_change_j,
            changed.energy().thermal_external_j,
        );
        assert!(changed.energy().internal_transfer_magnitude_j[0] > 0.0);
    });
}

#[test]
fn differing_temperature_laws_still_have_exact_samplewise_bath_equilibrium() {
    with_cx(|cx| {
        let original = thermal(
            vec![warm_cell(11, 400.0), warm_cell(12, 400.0)],
            &RADIATIVE_SENSIBLE,
            vec![],
        );
        let field = field(
            cx,
            &[1.0, 0.5],
            &original,
            &[
                absorption(&[(200.0, 0.1), (800.0, 3.0)]),
                absorption(&[(200.0, 2.0), (800.0, 0.3)]),
            ],
            &[scattering(&[(200.0, 0.1), (800.0, 2.0)])],
            &[11, 12],
        );
        let exchange = TemperatureRadiativeExchange::try_new(
            field,
            [2000.0],
            blackbody(400.0),
            &original,
            budget(),
            photon_config(),
            RadiativeExchangeTolerance {
                max_cell_standard_error_w: 0.0,
                max_total_standard_error_w: 0.0,
                max_unresolved_exchange_w: 0.0,
                max_unresolved_standard_error_w: 0.0,
                ..tolerance()
            },
            || Ok(()),
        )
        .unwrap();
        let (same, accepted, audit) = original
            .advance_observed(&exchange, 0.3, config(0.1), || Ok(()))
            .unwrap();
        assert_eq!(same.cells(), original.cells());
        assert_eq!(accepted.operator_builds(), 1);
        assert_eq!(accepted.energy().thermal_external_j, 0.0);
        assert_eq!(accepted.energy().internal_transfer_magnitude_j, [0.0]);
        assert_eq!(accepted.power().unwrap().total_standard_error_w, 0.0);
        assert_eq!(audit.external_heat_j, 0.0);
    });
}

#[test]
fn initially_transparent_material_can_become_an_emitter_after_native_heating() {
    with_cx(|cx| {
        let mut heated = warm_cell(11, 300.0);
        heated.boundary.source_w = 1000.0;
        let original = thermal(vec![heated], &RADIATIVE_SENSIBLE, vec![]);
        let field = field(
            cx,
            &[1.0],
            &original,
            &[absorption(&[(300.0, 0.0), (500.0, 2.0)])],
            &[scattering(&[(300.0, 0.0), (500.0, 0.0)])],
            &[11],
        );
        let exchange = observer(field, &original, 1000.0, RadiativeBath::Vacuum);
        let mut initial = exchange.clone();
        assert!(
            initial
                .sources(&original, &mut || Ok(()))
                .unwrap()
                .is_empty()
        );
        assert!(initial.operator().unwrap().cells.is_empty());
        assert_eq!(initial.traced_segments(), 0);
        let (warm, emitted, audit) = original
            .advance_observed(&exchange, 0.2, config(0.1), || Ok(()))
            .unwrap();
        // The radiative CFL can introduce an additional endpoint remainder;
        // every actual changed-temperature source needs its own operator.
        assert!(emitted.operator_builds() >= 2);
        assert_eq!(emitted.operator_builds(), audit.substeps);
        assert!(emitted.operator().unwrap().cells[0].emission_capacity_m2_sr[0] > 0.0);
        assert!(emitted.traced_segments() > 0);
        assert!(emitted.energy().thermal_external_j < 0.0);
        assert!(warm.temperature_k(CellId(11)).unwrap() > 300.0);
        relative_close(
            audit.external_heat_j,
            200.0 + emitted.energy().thermal_external_j,
        );
    });
}

#[test]
fn cached_zero_work_continuation_and_changed_replica_preserve_past_energy() {
    with_cx(|cx| {
        // A native event gives whole and split runs the same physical cut.
        // Arbitrarily splitting an interval would alter the explicit CFL grid.
        let original = thermal(
            vec![warm_cell(11, 500.0)],
            &RADIATIVE_SENSIBLE,
            vec![TimedEvent {
                time_s: 0.1,
                event: Event::Boundary(CellId(11), boundary(0.0)),
            }],
        );
        let field = field(
            cx,
            &[1.0],
            &original,
            &[absorption(&[(100.0, 0.7), (1000.0, 0.7)])],
            &[scattering(&[(100.0, 0.0), (1000.0, 0.0)])],
            &[11],
        );
        let exchange = observer(field, &original, 1000.0, RadiativeBath::Vacuum);
        let (whole, whole_observer, _) = original
            .advance_observed(&exchange, 0.2, config(0.05), || Ok(()))
            .unwrap();
        let (middle, first, _) = original
            .advance_observed(&exchange, 0.1, config(0.05), || Ok(()))
            .unwrap();
        let (split, second, _) = middle
            .advance_observed(
                &first.with_work_limits(100, 0),
                0.2,
                config(0.05),
                || Ok(()),
            )
            .unwrap();
        assert_eq!(whole.cells(), split.cells());
        assert_eq!(whole_observer.energy(), second.energy());
        assert_eq!(second.operator_builds(), 0);
        assert_eq!(second.traced_segments(), 0);
        assert!(second.source_evaluations() > 0);
        let different = first.with_replica_seed(77);
        assert_eq!(different.energy(), first.energy());
        assert_ne!(
            different.sampling_config().seed,
            first.sampling_config().seed
        );
        let doubled = first.with_histories_per_emitter_lane(256).unwrap();
        assert_eq!(doubled.sampling_config().seed, first.sampling_config().seed);
        assert_eq!(
            doubled.sampling_config().source_id,
            first.sampling_config().source_id
        );
        let (_, replicated, _) = middle
            .advance_observed(&different, 0.2, config(0.05), || Ok(()))
            .unwrap();
        assert_eq!(replicated.operator_builds(), 1);
        assert!(replicated.energy().thermal_external_j < first.energy().thermal_external_j);
        let (_, refined, _) = middle
            .advance_observed(&doubled, 0.2, config(0.05), || Ok(()))
            .unwrap();
        assert_eq!(refined.operator_builds(), 1);
        assert_eq!(
            refined.operator().unwrap().cells[0].sampled_histories,
            [256]
        );
    });
}

#[test]
fn aggregate_source_transport_budgets_and_cancellation_roll_back_feedback() {
    with_cx(|cx| {
        let original = thermal(vec![warm_cell(11, 500.0)], &RADIATIVE_SENSIBLE, vec![]);
        let field = field(
            cx,
            &[1.0],
            &original,
            &[absorption(&[(100.0, 0.1), (1000.0, 1.5)])],
            &[scattering(&[(100.0, 0.1), (1000.0, 0.5)])],
            &[11],
        );
        let exchange = observer(field, &original, 1000.0, RadiativeBath::Vacuum);
        let (_, reference, audit) = original
            .advance_observed(&exchange, 0.2, config(0.05), || Ok(()))
            .unwrap();
        assert!(audit.substeps >= 4);
        assert!(reference.operator_builds() >= 4);
        assert!(matches!(
            original.advance_observed(
                &exchange.with_work_limits(2, u64::MAX),
                0.2,
                config(0.05),
                || Ok(())
            ),
            Err(RadiativeFeedbackError::Budget)
        ));
        assert!(
            original
                .advance_observed(
                    &exchange.with_work_limits(1000, reference.traced_segments() - 1),
                    0.2,
                    config(0.05),
                    || Ok(())
                )
                .is_err()
        );
        let mut polls = 0usize;
        let cancelled = original.advance_observed(&exchange, 0.2, config(0.05), || {
            polls += 1;
            if polls >= 500 {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        });
        assert!(matches!(
            cancelled,
            Err(RadiativeFeedbackError::Exchange(
                RadiativeExchangeError::Thermal(DepositionError::Cancelled)
            ))
        ));
        assert_eq!(original.time_s(), 0.0);
        assert_eq!(original.cells()[&CellId(11)].specific_h_j_kg, 5000.0);
        assert!(exchange.operator().is_none());
        assert_eq!(exchange.energy().thermal_external_j, 0.0);
        assert_eq!(exchange.source_evaluations(), 0);
    });
}

#[test]
fn latent_storage_and_scheduled_native_sources_retain_split_history() {
    with_cx(|cx| {
        let mut initial = cell(11, 1.0, 1.0, 0.0);
        initial.specific_h_j_kg = 120.0;
        let original = thermal(
            vec![initial],
            &LATENT,
            vec![TimedEvent {
                time_s: 1.0,
                event: Event::Boundary(CellId(11), boundary(4.0)),
            }],
        );
        let field = field(
            cx,
            &[1.0],
            &original,
            &[absorption(&[(290.0, 0.5), (350.0, 1.5)])],
            &[scattering(&[(290.0, 0.2), (350.0, 1.0)])],
            &[11],
        );
        let exchange = observer(field, &original, 3000.0, blackbody(300.0));
        let (whole, complete, audit) = original
            .advance_observed(&exchange, 2.0, config(0.125), || Ok(()))
            .unwrap();
        let (middle, first, _) = original
            .advance_observed(&exchange, 1.0, config(0.125), || Ok(()))
            .unwrap();
        let (split, second, _) = middle
            .advance_observed(&first, 2.0, config(0.125), || Ok(()))
            .unwrap();
        assert_eq!(whole.cells(), split.cells());
        assert_eq!(complete.energy(), second.energy());
        assert_eq!(whole.temperature_k(CellId(11)).unwrap(), 310.0);
        assert!((10.0..110.0).contains(&whole.cells()[&CellId(11)].specific_h_j_kg));
        assert_eq!(audit.events, 1);
        relative_close(
            audit.external_heat_j,
            complete.energy().thermal_external_j + 4.0,
        );
    });
}

#[test]
fn no_emission_still_checks_material_temperature_domain_with_zero_photon_budget() {
    with_cx(|cx| {
        let mut initial = warm_cell(11, 300.0);
        initial.boundary.source_w = 1000.0;
        let original = thermal(vec![initial], &RADIATIVE_SENSIBLE, vec![]);
        let field = field(
            cx,
            &[1.0],
            &original,
            &[absorption(&[(290.0, 0.0), (305.0, 0.0)])],
            &[scattering(&[(290.0, 1.0), (305.0, 2.0)])],
            &[11],
        );
        let exchange =
            observer(field, &original, 1000.0, RadiativeBath::Vacuum).with_work_limits(10, 0);
        let (valid, accepted, _) = original
            .advance_observed(&exchange, 0.025, config(0.025), || Ok(()))
            .unwrap();
        assert_eq!(accepted.traced_segments(), 0);
        assert_eq!(accepted.energy().thermal_external_j, 0.0);
        assert_eq!(valid.temperature_k(CellId(11)).unwrap(), 302.5);
        assert!(
            original
                .advance_observed(&exchange, 0.1, config(0.1), || Ok(()))
                .is_err()
        );
        assert_eq!(original.temperature_k(CellId(11)).unwrap(), 300.0);
        assert!(exchange.operator().is_none());
    });
}
