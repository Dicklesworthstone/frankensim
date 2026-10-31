//! Whole-trajectory uncertainty, resolution and work of local radiative feedback.
use super::*;
use fs_render::volumes::heating::exchange::feedback::refinement::{
    RadiativeEnsembleTolerance, RadiativeEnsembleWork, RadiativeFeedbackRefinement,
    RadiativeFeedbackRefinementError,
};

fn accuracy(value: f64) -> RadiativeEnsembleTolerance {
    RadiativeEnsembleTolerance {
        enthalpy_absolute_j_kg: value,
        temperature_absolute_k: value,
        radiative_energy_absolute_j: value,
        relative: 0.0,
    }
}
fn controls() -> RadiativeFeedbackRefinement {
    RadiativeFeedbackRefinement {
        initial_max_dt_s: 0.25,
        minimum_max_dt_s: 1e-5,
        maximum_halvings: 12,
        replicas: 3,
        ensemble_seed: 267,
        maximum_photon_doublings: 1,
        standard_error_multiplier: 3.0,
        temporal: accuracy(100.0),
        sampling: accuracy(100.0),
        photon_sensitivity: accuracy(100.0),
        conservation_tolerance_j: 1e-8,
        maximum_total_substeps: 20_000,
        maximum_total_source_evaluations: 20_000,
        maximum_total_traced_segments: 2_000_000,
        maximum_retained_trajectories: 12,
        maximum_observable_cells: 8,
        maximum_retained_operator_pair_slots: 1024,
        maximum_retained_operator_cell_slots: 2048,
    }
}
fn sampled<'g>(
    optical_field: Arc<TemperatureOpticalField<'g, 1>>,
    original: &DepositionSimulation<'_>,
    weight: f64,
    bath: RadiativeBath,
) -> TemperatureRadiativeExchange<'g, 1> {
    TemperatureRadiativeExchange::try_new(
        optical_field,
        [weight],
        bath,
        original,
        budget(),
        RadiativeExchangeConfig {
            histories_per_emitter_lane: 8,
            ..photon_config()
        },
        tolerance(),
        || Ok(()),
    )
    .unwrap()
}

#[test]
fn pure_absorption_feedback_uses_h2h_and_matches_independent_thin_limit_cooling() {
    with_cx(|cx| {
        let original = thermal(vec![warm_cell(11, 500.0)], &RADIATIVE_SENSIBLE, vec![]);
        let source = sampled(
            field(
                cx,
                &[1.0],
                &original,
                &[absorption(&[(300.0, 1e-6), (800.0, 3e-6)])],
                &[scattering(&[(300.0, 0.0), (800.0, 0.0)])],
                &[11],
            ),
            &original,
            3000.0,
            RadiativeBath::Vacuum,
        );
        let mut c = controls();
        c.initial_max_dt_s = 100.0;
        c.minimum_max_dt_s = 0.1;
        c.temporal = RadiativeEnsembleTolerance {
            enthalpy_absolute_j_kg: 1e-4,
            temperature_absolute_k: 1e-5,
            radiative_energy_absolute_j: 1e-4,
            relative: 0.0,
        };
        c.sampling = c.temporal;
        c.photon_sensitivity = c.temporal;
        let result = source
            .advance_ensemble_refined(&original, 1000.0, c, || Ok(()))
            .unwrap();
        // Independent direct SI Planck/linear κa(T) RK4 owner in the thin limit.
        // For the unit cube, every outgoing chord is at most sqrt(3), so omitted
        // reabsorption changes this endpoint by <2e-6 K at these coefficients.
        // This is a physical limiting reference, separate from the production
        // collision sampler, table interpolation and native explicit stepper.
        let rhs = |temperature: f64| {
            let absorption = 1e-6 + (temperature - 300.0) * (2e-6 / 500.0);
            -4.0 * std::f64::consts::PI * absorption * 3000.0 * planck(temperature) / 10.0
        };
        let mut expected = 500.0;
        let dt = 1000.0 / 20_000.0;
        for _ in 0..20_000 {
            let a = rhs(expected);
            let b = rhs(expected + 0.5 * dt * a);
            let d = rhs(expected + 0.5 * dt * b);
            let e = rhs(expected + dt * d);
            expected += dt * (a + 2.0 * b + 2.0 * d + e) / 6.0;
        }
        assert!((result.cells[&CellId(11)].temperature_k.mean - expected).abs() < 1.5e-5);
        assert!(result.diagnostics.temporal_refinement_observed);
        assert!(
            result
                .diagnostics
                .temporal
                .unwrap()
                .maximum_normalized_width
                <= 1.0
        );
        assert!(
            result
                .diagnostics
                .sampling
                .unwrap()
                .maximum_normalized_width
                <= 1.0
        );
        assert!(
            result
                .diagnostics
                .photon_sensitivity
                .unwrap()
                .maximum_normalized_width
                <= 1.0
        );
        assert_eq!(result.diagnostics.photon_doublings, 1);
        assert_eq!(result.diagnostics.histories_per_emitter_lane, 16);
        assert!(result.diagnostics.work.operator_builds > c.replicas as usize * 4);
        for candidate in result.trajectories {
            assert_eq!(candidate.optical.histories_per_emitter_lane(), 16);
            assert!(candidate.thermal.temperature_k(CellId(11)).unwrap() < 500.0);
            relative_close(
                candidate.audit.external_heat_j,
                candidate.optical.energy().thermal_external_j,
            );
        }
        assert_eq!(original.time_s(), 0.0);
        assert!(source.operator().is_none());
    });
}

struct Curved;
impl ThermalLaw for Curved {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        if !h.is_finite() || !(0.0..=1000.0).contains(&h) {
            return Err(DepositionError::MaterialDomain);
        }
        Ok(250.0 + 0.1 * h + 0.001 * h * h)
    }
    fn maximum_temperature_slope(&self) -> f64 {
        2.1
    }
}
static CURVED: Curved = Curved;

#[test]
fn independent_operator_trajectories_replay_and_average_actual_caloric_temperature() {
    with_cx(|cx| {
        let mut initial_cell = warm_cell(11, 300.0);
        initial_cell.specific_h_j_kg = 300.0;
        let original = thermal(vec![initial_cell], &CURVED, vec![]);
        let source = sampled(
            field(
                cx,
                &[1.0],
                &original,
                &[absorption(&[(200.0, 0.3), (600.0, 1.3)])],
                &[scattering(&[(200.0, 0.8), (600.0, 0.2)])],
                &[11],
            ),
            &original,
            500.0,
            RadiativeBath::Vacuum,
        );
        let mut c = controls();
        c.replicas = 6;
        let result = source
            .advance_ensemble_refined(&original, 1.0, c, || Ok(()))
            .unwrap();
        let replay = source
            .advance_ensemble_refined(&original, 1.0, c, || Ok(()))
            .unwrap();
        assert_eq!(result.diagnostics, replay.diagnostics);
        assert_eq!(result.cells, replay.cells);
        assert_eq!(result.radiative_energy, replay.radiative_energy);
        let hs: Vec<_> = result
            .trajectories
            .iter()
            .map(|r| r.thermal.cells()[&CellId(11)].specific_h_j_kg)
            .collect();
        let ts: Vec<_> = result
            .trajectories
            .iter()
            .map(|r| r.thermal.temperature_k(CellId(11)).unwrap())
            .collect();
        let count = hs.len() as f64;
        let h_mean = hs.iter().sum::<f64>() / count;
        let t_mean = ts.iter().sum::<f64>() / count;
        let h_variance = hs.iter().map(|h| (h - h_mean).powi(2)).sum::<f64>() / (count - 1.0);
        let t_variance = ts.iter().map(|t| (t - t_mean).powi(2)).sum::<f64>() / (count - 1.0);
        let estimate = result.cells[&CellId(11)];
        relative_close(estimate.specific_h_j_kg.mean, h_mean);
        relative_close(estimate.temperature_k.mean, t_mean);
        relative_close(
            estimate.specific_h_j_kg.standard_error,
            (h_variance / count).sqrt(),
        );
        relative_close(
            estimate.temperature_k.standard_error,
            (t_variance / count).sqrt(),
        );
        assert!(estimate.temperature_k.standard_error > 0.0);
        // Independent quadratic-caloric identity: E[T(h)]-T(E[h])=b Var(h).
        let mean_state_temperature = CURVED.temperature_k(h_mean).unwrap();
        assert!(t_mean - mean_state_temperature > 1e-7);
        relative_close(
            t_mean - mean_state_temperature,
            0.001 * h_variance * (count - 1.0) / count,
        );
        let mut bath_j = 0.0;
        for (candidate, again) in result.trajectories.iter().zip(&replay.trajectories) {
            assert_eq!(candidate.thermal.cells(), again.thermal.cells());
            assert_eq!(candidate.optical.energy(), again.optical.energy());
            assert_eq!(
                candidate.optical.budget().max_source_evaluations,
                source.budget().max_source_evaluations
            );
            assert_eq!(
                candidate.optical.budget().max_traced_segments,
                source.budget().max_traced_segments
            );
            relative_close(
                candidate.audit.external_heat_j,
                candidate.optical.energy().thermal_external_j,
            );
            relative_close(
                candidate.audit.stored_energy_change_j,
                candidate.audit.external_heat_j,
            );
            bath_j += candidate.optical.energy().thermal_external_j;
        }
        relative_close(
            result.radiative_energy.thermal_external_j.mean,
            bath_j / count,
        );
        assert!(result.radiative_energy.thermal_external_j.standard_error > 0.0);
        assert!(
            result.diagnostics.work.traced_segments
                > result
                    .trajectories
                    .iter()
                    .map(|t| t.optical.traced_segments())
                    .sum::<u64>()
        );
    });
}

#[test]
fn sampling_photon_and_real_time_gates_are_independent() {
    with_cx(|cx| {
        let original = thermal(vec![warm_cell(11, 500.0)], &RADIATIVE_SENSIBLE, vec![]);
        let optical_field = field(
            cx,
            &[1.0],
            &original,
            &[absorption(&[(100.0, 0.6), (1000.0, 0.6)])],
            &[scattering(&[(100.0, 0.0), (1000.0, 0.0)])],
            &[11],
        );
        let source = sampled(
            optical_field.clone(),
            &original,
            1000.0,
            RadiativeBath::Vacuum,
        );
        let mut c = controls();
        c.initial_max_dt_s = 0.025;
        c.sampling = accuracy(1e-8);
        let Err(RadiativeFeedbackRefinementError::Unresolved(noise)) =
            source.advance_ensemble_refined(&original, 0.1, c, || Ok(()))
        else {
            panic!("finite independent-state sampling noise must refuse")
        };
        assert!(noise.temporal_refinement_observed);
        assert!(noise.temporal.unwrap().maximum_normalized_width <= 1.0);
        assert!(noise.sampling.unwrap().maximum_normalized_width > 1.0);
        assert!(noise.photon_sensitivity.is_none());
        c.sampling = accuracy(100.0);
        c.photon_sensitivity = accuracy(1e-8);
        let Err(RadiativeFeedbackRefinementError::Unresolved(photon)) =
            source.advance_ensemble_refined(&original, 0.1, c, || Ok(()))
        else {
            panic!("finite H/2H sensitivity must refuse")
        };
        assert!(photon.sampling.unwrap().maximum_normalized_width <= 1.0);
        assert!(
            photon
                .photon_sensitivity
                .unwrap()
                .maximum_radiative_energy_width_j
                > 1e-8
        );
        assert!(photon.photon_sensitivity.unwrap().maximum_normalized_width > 1.0);

        // Exactly isothermal heat is zero, but a high conductance still sets
        // the native grid. Lowering an inactive maximum dt proves nothing.
        let equilibrium = sampled(optical_field, &original, 1e6, blackbody(500.0));
        let mut c = controls();
        c.initial_max_dt_s = 100.0;
        c.minimum_max_dt_s = 1.0;
        c.maximum_halvings = 2;
        let Err(RadiativeFeedbackRefinementError::Unresolved(time)) =
            equilibrium.advance_ensemble_refined(&original, 0.02, c, || Ok(()))
        else {
            panic!("unchanged native CFL grids cannot demonstrate refinement")
        };
        assert!(!time.temporal_refinement_observed);
        assert_eq!(time.temporal.unwrap().maximum_enthalpy_width_j_kg, 0.0);
        assert!(time.sampling.is_none());
        assert_eq!(source.energy().thermal_external_j, 0.0);
    });
}

#[test]
fn all_discarded_work_and_sparse_operator_storage_are_bounded_transactionally() {
    with_cx(|cx| {
        let original = thermal(vec![warm_cell(11, 500.0)], &RADIATIVE_SENSIBLE, vec![]);
        let source = sampled(
            field(
                cx,
                &[1.0],
                &original,
                &[absorption(&[(100.0, 0.2), (1000.0, 0.2)])],
                &[scattering(&[(100.0, 0.0), (1000.0, 0.0)])],
                &[11],
            ),
            &original,
            20.0,
            RadiativeBath::Vacuum,
        );
        let c = controls();
        let result = source
            .advance_ensemble_refined(&original, 1.0, c, || Ok(()))
            .unwrap();
        assert_eq!(result.diagnostics.halvings, 2);
        // Positive radiative conductance harmonically shortens the proposed
        // maximum dt, including when that maximum dominates the native bound.
        // Sum independently invoked native trials, not target/max_dt guesses.
        let high = 2 * source.histories_per_emitter_lane();
        let mut expected = RadiativeEnsembleWork::default();
        for (max_dt_s, histories) in [
            (c.initial_max_dt_s, high),
            (c.initial_max_dt_s * 0.5, high),
            (c.initial_max_dt_s * 0.25, high),
            (c.initial_max_dt_s * 0.25, high / 2),
        ] {
            for replica in 0..c.replicas {
                let trial = source
                    .with_replica_seed(c.ensemble_seed.wrapping_add(u64::from(replica)))
                    .with_histories_per_emitter_lane(histories)
                    .unwrap();
                let (_, observed, audit) = original
                    .advance_observed(
                        &trial,
                        1.0,
                        StepConfig {
                            max_dt_s,
                            energy_tolerance_j: c.conservation_tolerance_j,
                        },
                        || Ok(()),
                    )
                    .unwrap();
                expected.substeps += audit.substeps;
                expected.source_evaluations += observed.source_evaluations();
                expected.traced_segments += observed.traced_segments();
                expected.operator_builds += observed.operator_builds();
            }
        }
        assert_eq!(result.diagnostics.work, expected);
        // Constant packets build once per trial; each pure-absorption source
        // history takes one extinction segment, regardless of its endpoint.
        assert_eq!(result.diagnostics.work.operator_builds, 3 * 4);
        assert_eq!(
            result.diagnostics.work.traced_segments,
            3 * (16 + 16 + 16 + 8)
        );
        let mut cases = [controls(); 5];
        cases[0].maximum_total_substeps = result.diagnostics.work.substeps - 1;
        cases[1].maximum_total_source_evaluations = result.diagnostics.work.source_evaluations - 1;
        cases[2].maximum_total_traced_segments = result.diagnostics.work.traced_segments - 1;
        cases[3].maximum_retained_operator_pair_slots =
            result.diagnostics.admitted_peak_operator_pair_slots - 1;
        cases[4].maximum_retained_operator_cell_slots =
            result.diagnostics.admitted_peak_operator_cell_slots - 1;
        for c in cases {
            assert!(matches!(
                source.advance_ensemble_refined(&original, 1.0, c, || Ok(())),
                Err(RadiativeFeedbackRefinementError::Budget(_))
                    | Err(RadiativeFeedbackRefinementError::Physics(_))
            ));
        }
        let mut checkpoints = 0;
        let cancelled = source.advance_ensemble_refined(&original, 1.0, controls(), || {
            checkpoints += 1;
            if checkpoints == 120 {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        });
        assert!(matches!(
            cancelled,
            Err(RadiativeFeedbackRefinementError::Physics(
                RadiativeFeedbackError::Exchange(RadiativeExchangeError::Thermal(
                    DepositionError::Cancelled
                ))
            )) | Err(RadiativeFeedbackRefinementError::Physics(
                RadiativeFeedbackError::Coefficients(
                    fs_render::volumes::heating::feedback::OpticalFeedbackError::Thermal(
                        DepositionError::Cancelled
                    )
                )
            ))
        ));
        assert_eq!(original.time_s(), 0.0);
        assert_eq!(original.temperature_k(CellId(11)).unwrap(), 500.0);
        assert!(source.operator().is_none());
        assert_eq!(source.energy().thermal_external_j, 0.0);
    });
}

#[test]
fn accepted_cached_trajectories_restore_original_source_limits_for_continuation() {
    with_cx(|cx| {
        let original = thermal(vec![warm_cell(11, 500.0)], &RADIATIVE_SENSIBLE, vec![]);
        let source = sampled(
            field(
                cx,
                &[1.0],
                &original,
                &[absorption(&[(100.0, 0.2), (1000.0, 0.2)])],
                &[scattering(&[(100.0, 0.0), (1000.0, 0.0)])],
                &[11],
            ),
            &original,
            20.0,
            RadiativeBath::Vacuum,
        );
        let mut c = controls();
        c.maximum_total_source_evaluations = 160;
        let result = source
            .advance_ensemble_refined(&original, 1.0, c, || Ok(()))
            .unwrap();
        assert!(result.diagnostics.work.source_evaluations <= 160);
        let accepted = result.trajectories.last().unwrap();
        assert_eq!(
            accepted.optical.budget().max_source_evaluations,
            source.budget().max_source_evaluations
        );
        let before = accepted.optical.energy().thermal_external_j;
        let (continued, observed, audit) = accepted
            .thermal
            .advance_observed(
                &accepted.optical,
                3.0,
                StepConfig {
                    max_dt_s: 0.01,
                    energy_tolerance_j: c.conservation_tolerance_j,
                },
                || Ok(()),
            )
            .unwrap();
        // This exceeds even the WHOLE study allowance, and therefore every
        // tightened remaining allowance inherited by its accepted trials.
        // A cached fixed owner must use the restored per-advance source limit.
        assert!(observed.source_evaluations() > c.maximum_total_source_evaluations);
        assert!(observed.source_evaluations() <= source.budget().max_source_evaluations);
        assert_eq!(observed.source_evaluations(), audit.substeps);
        assert_eq!(observed.operator_builds(), 0);
        assert_eq!(observed.traced_segments(), 0);
        assert_eq!(continued.time_s(), 3.0);
        assert_eq!(accepted.thermal.time_s(), 1.0);
        assert!(
            continued.temperature_k(CellId(11)).unwrap()
                < accepted.thermal.temperature_k(CellId(11)).unwrap()
        );
        relative_close(
            audit.external_heat_j,
            observed.energy().thermal_external_j - before,
        );
        relative_close(audit.stored_energy_change_j, audit.external_heat_j);
    });
}

#[test]
fn vacuum_grid_packets_count_in_memory_and_zero_emission_needs_zero_photon_work() {
    with_cx(|cx| {
        let mut active = warm_cell(11, 500.0);
        active.boundary.source_w = 1.0;
        let original = thermal(vec![active], &RADIATIVE_SENSIBLE, vec![]);
        let optical_field = Arc::new(
            TemperatureOpticalField::try_new(
                grid(
                    cx,
                    &[1.0, 0.0, 0.0, 0.0],
                    MassAbsorptionSpectrum::constant(0.0).unwrap(),
                ),
                [10_000.0],
                &[absorption(&[(100.0, 0.0), (1000.0, 0.0)])],
                &[scattering(&[(100.0, 1.0), (1000.0, 1.0)])],
                &[PhaseFunction::Isotropic],
                &[ThermalScatteringCell {
                    grid_cell: 0,
                    thermal_cell: CellId(11),
                    absorption_table: 0,
                    scattering_table: 0,
                    phase: 0,
                }],
                &original,
                TemperatureOpticalFieldBudget {
                    max_grid_cells: 4,
                    max_tables: 2,
                    max_total_temperature_knots: 4,
                },
                || Ok(()),
            )
            .unwrap(),
        );
        let source = TemperatureRadiativeExchange::try_new(
            optical_field,
            [1000.0],
            RadiativeBath::Vacuum,
            &original,
            RadiativeFeedbackBudget {
                max_traced_segments: 0,
                ..budget()
            },
            RadiativeExchangeConfig {
                max_material_cells: 1,
                max_emitter_cells: 0,
                max_pairs: 0,
                max_total_segments: 0,
                histories_per_emitter_lane: 8,
                ..photon_config()
            },
            tolerance(),
            || Ok(()),
        )
        .unwrap();
        let mut c = controls();
        c.maximum_total_traced_segments = 0;
        c.maximum_retained_operator_pair_slots = 0;
        // Material/emitter metadata alone would need 10 slots. Mutable packets
        // also span the four original geometry cells, including vacuum.
        c.maximum_retained_operator_cell_slots = 10;
        assert!(matches!(
            source.advance_ensemble_refined(&original, 1.0, c, || Ok(())),
            Err(RadiativeFeedbackRefinementError::Budget(_))
        ));
        c.maximum_retained_operator_cell_slots = 50;
        let result = source
            .advance_ensemble_refined(&original, 1.0, c, || Ok(()))
            .unwrap();
        assert_eq!(result.diagnostics.admitted_peak_operator_cell_slots, 50);
        assert_eq!(result.diagnostics.work.traced_segments, 0);
        assert!(result.diagnostics.photon_sensitivity.is_some());
        assert_eq!(result.radiative_energy.thermal_external_j.mean, 0.0);
        assert_eq!(
            result.radiative_energy.thermal_external_j.standard_error,
            0.0
        );
        assert!(result.cells[&CellId(11)].temperature_k.mean > 500.0);
    });
}

struct FormationDatum;
impl ThermalLaw for FormationDatum {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        if !h.is_finite() {
            return Err(DepositionError::MaterialDomain);
        }
        Ok(300.0)
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0
    }
}
static FORMATION: FormationDatum = FormationDatum;

#[test]
fn identical_rounded_enthalpies_cannot_claim_sub_ulp_accuracy() {
    with_cx(|cx| {
        let mut active = warm_cell(11, 300.0);
        active.specific_h_j_kg = 1e20;
        let original = thermal(vec![active], &FORMATION, vec![]);
        let source = sampled(
            field(
                cx,
                &[1.0],
                &original,
                &[absorption(&[(100.0, 0.2), (1000.0, 0.2)])],
                &[scattering(&[(100.0, 0.0), (1000.0, 0.0)])],
                &[11],
            ),
            &original,
            20.0,
            blackbody(300.0),
        );
        let mut c = controls();
        c.maximum_halvings = 2;
        c.temporal.enthalpy_absolute_j_kg = 1e-3;
        let Err(RadiativeFeedbackRefinementError::Unresolved(report)) =
            source.advance_ensemble_refined(&original, 1.0, c, || Ok(()))
        else {
            panic!("equal rounded enthalpy is not sub-ULP convergence evidence")
        };
        assert!(report.temporal_refinement_observed);
        let time = report.temporal.unwrap();
        assert_eq!(time.maximum_enthalpy_width_j_kg, 0.0);
        assert!(time.maximum_coordinate_resolution_ratio > 1.0);
        assert!(time.maximum_normalized_width > 1.0);
        assert!(report.sampling.is_none());
    });
}
