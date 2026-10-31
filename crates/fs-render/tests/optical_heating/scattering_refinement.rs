//! Whole-trajectory numerical assessment of the coupled native observer.
use super::*;
use fs_blake3::ContentHash;
use fs_render::volumes::heating::feedback::{OpticalFeedbackError, TemperatureAbsorptionTable};
use fs_render::volumes::heating::scattering::ScatteringHeatingTolerance;
use fs_render::volumes::heating::scattering_feedback::{
    IncidentScatteringBeam, ScatteringEnsembleTolerance, ScatteringFeedbackBudget,
    ScatteringFeedbackRefinement, ScatteringRefinementError, TemperatureScatteringHeating,
    TemperatureScatteringTable, ThermalScatteringCell,
};
use fs_render::volumes::transport::scattering::deposition::ScatteringDepositionConfig;

fn observer<'g>(
    cx: &Cx<'_>,
    density: &'g [f64],
    original: &DepositionSimulation<'_>,
    power: f64,
    absorption: [f64; 2],
    scattering: [f64; 2],
    upper_k: f64,
) -> TemperatureScatteringHeating<'g, 1> {
    let a = TemperatureAbsorptionTable::try_new(
        ContentHash([51; 32]),
        &[
            (
                300.0,
                MassAbsorptionSpectrum::constant(absorption[0]).unwrap(),
            ),
            (
                upper_k,
                MassAbsorptionSpectrum::constant(absorption[1]).unwrap(),
            ),
        ],
        2,
    )
    .unwrap();
    let s = TemperatureScatteringTable::try_new(
        ContentHash([52; 32]),
        &[
            (
                300.0,
                MassScatteringSpectrum::constant(scattering[0]).unwrap(),
            ),
            (
                upper_k,
                MassScatteringSpectrum::constant(scattering[1]).unwrap(),
            ),
        ],
        2,
    )
    .unwrap();
    TemperatureScatteringHeating::try_new(
        grid(cx, density, MassAbsorptionSpectrum::constant(1.0).unwrap()),
        IncidentScatteringBeam {
            origin: Point3::new(-1.0, 0.5, 0.5),
            unit_direction: Vec3::new(1.0, 0.0, 0.0),
            wavelengths_nm: [550.0],
            incident_power_w: [power],
        },
        &[a],
        &[s],
        &[PhaseFunction::Isotropic],
        &[ThermalScatteringCell {
            grid_cell: 0,
            thermal_cell: CellId(11),
            absorption_table: 0,
            scattering_table: 0,
            phase: 0,
        }],
        original,
        ScatteringFeedbackBudget {
            max_grid_cells: 1,
            max_tables: 2,
            max_total_temperature_knots: 4,
            max_source_evaluations: 10_000,
            max_traced_segments: 1_000_000,
        },
        ScatteringDepositionConfig {
            histories_per_lane: 8,
            max_collisions_per_history: 32,
            max_total_segments: 50_000,
            max_deposition_cells: 1,
            seed: 125,
            source_id: 7,
        },
        ScatteringHeatingTolerance {
            power_roundoff_w: 1e-10,
            max_unresolved_w: 100.0,
            max_total_standard_error_w: 100.0,
            max_cell_standard_error_w: 100.0,
        },
        || Ok(()),
    )
    .unwrap()
}

fn tolerance(value: f64) -> ScatteringEnsembleTolerance {
    ScatteringEnsembleTolerance {
        enthalpy_absolute_j_kg: value,
        temperature_absolute_k: value,
        optical_energy_absolute_j: value,
        relative: 0.0,
    }
}
fn controls() -> ScatteringFeedbackRefinement {
    ScatteringFeedbackRefinement {
        initial_max_dt_s: 0.5,
        minimum_max_dt_s: 0.001,
        maximum_halvings: 8,
        replicas: 8,
        ensemble_seed: 173,
        maximum_photon_doublings: 1,
        standard_error_multiplier: 3.0,
        temporal: tolerance(100.0),
        sampling: tolerance(100.0),
        photon_sensitivity: tolerance(100.0),
        conservation_tolerance_j: 1e-9,
        maximum_total_substeps: 20_000,
        maximum_total_source_evaluations: 20_000,
        maximum_total_traced_segments: 2_000_000,
        maximum_retained_trajectories: 16,
        maximum_observable_cells: 1,
    }
}

#[test]
fn deterministic_limit_refines_the_native_nonlinear_beer_solution_without_sampling_noise() {
    with_cx(|cx| {
        let original = thermal(vec![cell(11, 1.0, 1.0, 0.0)], &SENSIBLE, vec![]);
        let optical = observer(cx, &[1.0], &original, 4.0, [0.2, 10.2], [0.0; 2], 400.0);
        let mut c = controls();
        c.replicas = 3;
        c.maximum_retained_trajectories = 6;
        c.maximum_photon_doublings = 0;
        c.temporal = ScatteringEnsembleTolerance {
            enthalpy_absolute_j_kg: 0.004,
            temperature_absolute_k: 0.001,
            optical_energy_absolute_j: 0.004,
            relative: 0.0,
        };
        c.sampling = tolerance(1e-12);
        let result = optical
            .advance_ensemble_refined(&original, 5.0, c, || Ok(()))
            .unwrap();
        // Independent exact solution: y=exp[a+b(T-T0)], y'=b P/(m cp)(y-1).
        let exact_t = 300.0 + ((1.0 + 0.2_f64.exp_m1() * 0.5_f64.exp()).ln() - 0.2) / 0.1;
        let estimate = result.cells[&CellId(11)];
        assert!((estimate.temperature_k.mean - exact_t).abs() < 0.0011);
        assert_eq!(estimate.temperature_k.standard_error, 0.0);
        assert_eq!(estimate.specific_h_j_kg.standard_error, 0.0);
        assert_eq!(result.optical_energy.absorbed_j[0].standard_error, 0.0);
        assert_eq!(result.diagnostics.photon_sensitivity, None);
        assert_eq!(result.diagnostics.photon_doublings, 0);
        assert!(
            result
                .diagnostics
                .temporal
                .unwrap()
                .maximum_normalized_width
                <= 1.0
        );
        assert!(result.diagnostics.halvings >= 2);
        assert_eq!(result.trajectories.len(), 3);
        for candidate in &result.trajectories {
            close(
                candidate.thermal.temperature_k(CellId(11)).unwrap(),
                estimate.temperature_k.mean,
            );
            close(
                candidate.audit.external_heat_j,
                candidate.optical.energy().thermal_input_j,
            );
            close(
                candidate.optical.energy().absorbed_j[0],
                candidate.audit.external_heat_j,
            );
        }
        assert_eq!(original.time_s(), 0.0);
        assert_eq!(optical.energy().thermal_input_j, 0.0);
    });
}

struct CurvedCaloric;
impl ThermalLaw for CurvedCaloric {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        if !h.is_finite() || !(0.0..=1000.0).contains(&h) {
            return Err(DepositionError::MaterialDomain);
        }
        Ok(300.0 + 0.1 * h + 0.001 * h * h)
    }
    fn maximum_temperature_slope(&self) -> f64 {
        2.1
    }
}
static CURVED: CurvedCaloric = CurvedCaloric;

#[test]
fn scattering_ensembles_replay_and_average_each_actual_caloric_temperature() {
    with_cx(|cx| {
        let original = thermal(vec![cell(11, 1.0, 1.0, 0.0)], &CURVED, vec![]);
        let optical = observer(cx, &[1.0], &original, 30.0, [0.4, 1.4], [2.0, 4.0], 500.0);
        let a = optical
            .advance_ensemble_refined(&original, 0.5, controls(), || Ok(()))
            .unwrap();
        let b = optical
            .advance_ensemble_refined(&original, 0.5, controls(), || Ok(()))
            .unwrap();
        assert_eq!(a.cells, b.cells);
        assert_eq!(a.optical_energy, b.optical_energy);
        assert_eq!(a.diagnostics, b.diagnostics);
        assert_eq!(a.diagnostics.halvings, 2);
        assert_eq!(a.diagnostics.photon_doublings, 1);
        assert_eq!(a.diagnostics.histories_per_lane, 16);
        // Three time levels (1+2+4 steps) plus paired low photons at four steps.
        assert_eq!(a.diagnostics.work.substeps, 8 * (1 + 2 + 4 + 4));
        assert_eq!(
            a.diagnostics.work.substeps,
            a.diagnostics.work.source_evaluations
        );
        assert!(a.diagnostics.work.traced_segments > 0);
        let hs: Vec<_> = a
            .trajectories
            .iter()
            .map(|r| r.thermal.cells()[&CellId(11)].specific_h_j_kg)
            .collect();
        let ts: Vec<_> = a
            .trajectories
            .iter()
            .map(|r| r.thermal.temperature_k(CellId(11)).unwrap())
            .collect();
        let mean_h = hs.iter().sum::<f64>() / hs.len() as f64;
        let mean_t = ts.iter().sum::<f64>() / ts.len() as f64;
        let variance_h = hs.iter().map(|h| (h - mean_h).powi(2)).sum::<f64>() / hs.len() as f64;
        let estimate = a.cells[&CellId(11)];
        close(estimate.specific_h_j_kg.mean, mean_h);
        close(estimate.temperature_k.mean, mean_t);
        close(
            estimate.specific_h_j_kg.standard_error,
            (variance_h / (hs.len() - 1) as f64).sqrt(),
        );
        assert!(estimate.temperature_k.standard_error > 0.0);
        // For this convex manufactured caloric law, E[T(h)]-T(E[h])=0.001 Var(h).
        let gap = estimate.temperature_k.mean - CURVED.temperature_k(mean_h).unwrap();
        close(gap, 0.001 * variance_h);
        assert!(gap > 1e-7);
        assert!(a.diagnostics.sampling.unwrap().maximum_normalized_width <= 1.0);
        assert!(
            a.diagnostics
                .photon_sensitivity
                .unwrap()
                .maximum_normalized_width
                <= 1.0
        );
        for (left, right) in a.trajectories.iter().zip(&b.trajectories) {
            assert_eq!(
                left.thermal.cells()[&CellId(11)].specific_h_j_kg,
                right.thermal.cells()[&CellId(11)].specific_h_j_kg
            );
            close(
                left.audit.external_heat_j,
                left.optical.energy().thermal_input_j,
            );
            let e = left.optical.energy();
            close(
                e.incident_j[0],
                e.ballistic_escape_j[0]
                    + e.diffuse_escape_j[0]
                    + e.absorbed_j[0]
                    + e.unresolved_j[0]
                    + e.residual_j[0],
            );
        }
        assert_eq!(original.cells()[&CellId(11)].specific_h_j_kg, 0.0);
        assert_eq!(optical.time_s(), 0.0);
        assert_eq!(optical.energy().thermal_input_j, 0.0);
    });
}

#[test]
fn paired_time_agreement_does_not_replace_final_sampling_or_photon_resolution_admission() {
    with_cx(|cx| {
        let original = thermal(vec![cell(11, 1.0, 1.0, 0.0)], &CURVED, vec![]);
        let optical = observer(cx, &[1.0], &original, 30.0, [0.4, 1.4], [2.0, 4.0], 500.0);
        let mut c = controls();
        c.sampling = tolerance(1e-16);
        let result = optical.advance_ensemble_refined(&original, 0.5, c, || Ok(()));
        let Err(ScatteringRefinementError::Unresolved(report)) = result else {
            panic!("sampling gate must refuse");
        };
        assert!(report.temporal.unwrap().maximum_normalized_width <= 1.0);
        assert!(report.sampling.unwrap().maximum_normalized_width > 1.0);
        assert_eq!(report.photon_sensitivity, None);
        c = controls();
        c.photon_sensitivity = tolerance(1e-16);
        let result = optical.advance_ensemble_refined(&original, 0.5, c, || Ok(()));
        let Err(ScatteringRefinementError::Unresolved(report)) = result else {
            panic!("photon gate must refuse");
        };
        assert!(report.temporal.unwrap().maximum_normalized_width <= 1.0);
        assert!(report.sampling.unwrap().maximum_normalized_width <= 1.0);
        assert!(report.photon_sensitivity.unwrap().maximum_normalized_width > 1.0);
        assert_eq!(optical.source_evaluations(), 0);
        assert_eq!(original.time_s(), 0.0);
    });
}

#[test]
fn whole_study_work_storage_cancellation_and_source_domains_are_transactional() {
    with_cx(|cx| {
        let original = thermal(vec![cell(11, 1.0, 1.0, 0.0)], &CURVED, vec![]);
        let optical = observer(cx, &[1.0], &original, 30.0, [0.4, 1.4], [2.0, 4.0], 500.0);
        for c in [
            ScatteringFeedbackRefinement {
                maximum_total_substeps: 24,
                ..controls()
            },
            ScatteringFeedbackRefinement {
                maximum_total_source_evaluations: 24,
                ..controls()
            },
            ScatteringFeedbackRefinement {
                maximum_total_traced_segments: 1,
                ..controls()
            },
            ScatteringFeedbackRefinement {
                maximum_retained_trajectories: 15,
                ..controls()
            },
            ScatteringFeedbackRefinement {
                maximum_observable_cells: 0,
                ..controls()
            },
        ] {
            assert!(
                optical
                    .advance_ensemble_refined(&original, 0.5, c, || Ok(()))
                    .is_err()
            );
        }
        let mut calls = 0;
        let cancelled = optical.advance_ensemble_refined(&original, 0.5, controls(), || {
            calls += 1;
            if calls > 25 {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        });
        assert!(matches!(
            cancelled,
            Err(ScatteringRefinementError::Physics(
                OpticalFeedbackError::Thermal(DepositionError::Cancelled)
            ))
        ));
        let narrow = observer(cx, &[1.0], &original, 30.0, [0.4, 0.4], [2.0; 2], 300.001);
        assert!(matches!(
            narrow.advance_ensemble_refined(&original, 0.5, controls(), || Ok(())),
            Err(ScatteringRefinementError::Physics(
                OpticalFeedbackError::TemperatureDomain { .. }
            ))
        ));
        assert_eq!(original.time_s(), 0.0);
        assert_eq!(original.cells()[&CellId(11)].specific_h_j_kg, 0.0);
        assert_eq!(optical.time_s(), 0.0);
        assert_eq!(optical.energy().thermal_input_j, 0.0);
        assert_eq!(optical.source_evaluations(), 0);
        assert_eq!(optical.traced_segments(), 0);
    });
}

#[test]
fn inactive_timestep_ceilings_and_sub_resolution_accuracy_cannot_masquerade_as_convergence() {
    with_cx(|cx| {
        let original = thermal(vec![cell(11, 1.0, 1.0, 0.0)], &SENSIBLE, vec![]);
        let optical = observer(cx, &[1.0], &original, 4.0, [0.0; 2], [0.0; 2], 400.0);
        let mut c = controls();
        c.initial_max_dt_s = 1.0;
        c.maximum_halvings = 2;
        // All three requested maxima exceed the physical interval, so they
        // produce the same single native step and no time-resolution evidence.
        let result = optical.advance_ensemble_refined(&original, 0.05, c, || Ok(()));
        let Err(ScatteringRefinementError::Unresolved(report)) = result else {
            panic!("clipped identical time grids do not establish refinement");
        };
        assert!(!report.temporal_refinement_observed);
        assert!(report.temporal.unwrap().maximum_normalized_width <= 1.0);
        c.maximum_halvings = 3;
        let admitted = optical
            .advance_ensemble_refined(&original, 0.5, c, || Ok(()))
            .unwrap();
        assert_eq!(admitted.diagnostics.halvings, 3);
        assert!(admitted.diagnostics.temporal_refinement_observed);
        assert_eq!(
            admitted.cells[&CellId(11)].temperature_k.standard_error,
            0.0
        );
        // Temporal differences can be exactly zero in rounded native states.
        // This does not provide temperature accuracy below a binary64 ULP.
        c = controls();
        c.maximum_halvings = 2;
        c.temporal = tolerance(1e-16);
        let result = optical.advance_ensemble_refined(&original, 0.5, c, || Ok(()));
        let Err(ScatteringRefinementError::Unresolved(report)) = result else {
            panic!("sub-resolution temperature tolerance must be refused");
        };
        assert!(report.temporal_refinement_observed);
        assert_eq!(report.temporal.unwrap().maximum_temperature_width_k, 0.0);
        assert!(report.temporal.unwrap().maximum_coordinate_resolution_ratio > 1.0);
        assert_eq!(original.time_s(), 0.0);
        assert_eq!(optical.time_s(), 0.0);
    });
}
