//! Joint visible absorption and infrared cooling use complete native replicas.
use super::*;
use fs_render::volumes::heating::joint::refinement::{
    JointOpticalRefinement, JointOpticalRefinementError, JointOpticalTolerance,
};

fn allowance(value: f64) -> JointOpticalTolerance {
    JointOpticalTolerance {
        enthalpy_absolute_j_kg: value,
        temperature_absolute_k: value,
        illumination_energy_absolute_j: value,
        radiation_energy_absolute_j: value,
        net_optical_heat_absolute_j: value,
        relative: 0.0,
    }
}
fn study() -> JointOpticalRefinement {
    JointOpticalRefinement {
        initial_max_dt_s: 0.25,
        minimum_max_dt_s: 1e-6,
        maximum_halvings: 2,
        replicas: 3,
        ensemble_seed: 2781,
        maximum_photon_doublings: 1,
        standard_error_multiplier: 3.0,
        temporal: allowance(100.0),
        sampling: allowance(100.0),
        illumination_photon_sensitivity: allowance(100.0),
        radiation_photon_sensitivity: allowance(100.0),
        conservation_tolerance_j: 1e-8,
        maximum_total_substeps: 20_000,
        maximum_total_illumination_source_evaluations: 20_000,
        maximum_total_radiation_source_evaluations: 20_000,
        maximum_total_illumination_traced_segments: 5_000_000,
        maximum_total_radiation_traced_segments: 5_000_000,
        maximum_retained_trajectories: 12,
        maximum_observable_cells: 8,
        maximum_retained_operator_pair_slots: 4096,
        maximum_retained_operator_cell_slots: 4096,
        maximum_retained_illumination_report_slots: 4096,
    }
}
fn near(a: f64, b: f64) {
    assert!(
        (a - b).abs() <= 2e-10 * a.abs().max(b.abs()).max(1.0),
        "{a:.17e} != {b:.17e}"
    );
}
fn mean_se(values: &[f64]) -> (f64, f64) {
    let count = values.len() as f64;
    let mean = values.iter().sum::<f64>() / count;
    let variance = values.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / (count - 1.0);
    (mean, (variance / count).sqrt())
}

#[test]
fn independent_visible_and_infrared_replicas_are_actual_joint_native_trajectories() {
    with_cx(|cx| {
        let original = thermal(vec![joint_cell(11, 500.0)], &JOINT_CALORIC, vec![]);
        let optical = joint_fixture(cx, &[1.0], &original, 25.0, 0.6, true);
        let c = study();
        let result = optical
            .advance_ensemble_refined(&original, 1.0, c, || Ok(()))
            .unwrap();
        assert_eq!(result.trajectories.len(), c.replicas as usize);
        assert!(result.diagnostics.temporal_refinement_observed);
        for gate in [
            result.diagnostics.temporal,
            result.diagnostics.sampling,
            result.diagnostics.illumination_photon_sensitivity,
            result.diagnostics.radiation_photon_sensitivity,
        ] {
            assert!(gate.unwrap().maximum_normalized_width <= 1.0);
        }
        assert_eq!(result.diagnostics.illumination_histories_per_lane, 32);
        assert_eq!(result.diagnostics.radiation_histories_per_emitter_lane, 32);
        let mut final_temperatures = Vec::new();
        for (replica, actual) in result.trajectories.iter().enumerate() {
            let seed = c.ensemble_seed.wrapping_add(replica as u64);
            let incident = optical
                .illumination()
                .with_replica_seed(seed)
                .with_histories_per_lane(result.diagnostics.illumination_histories_per_lane)
                .unwrap();
            let radiation = optical
                .radiation()
                .with_replica_seed(seed)
                .with_histories_per_emitter_lane(
                    result.diagnostics.radiation_histories_per_emitter_lane,
                )
                .unwrap();
            let pair = fs_thermal::observers::CoupledObservers::new(incident, radiation);
            let (expected, owners, audit) = original
                .advance_observed(
                    &pair,
                    1.0,
                    StepConfig {
                        max_dt_s: result.diagnostics.maximum_dt_s,
                        energy_tolerance_j: c.conservation_tolerance_j,
                    },
                    || Ok(()),
                )
                .unwrap();
            assert_eq!(actual.thermal.cells(), expected.cells());
            assert_eq!(actual.audit, audit);
            assert_eq!(
                actual.optical.illumination().energy(),
                owners.first.energy()
            );
            assert_eq!(actual.optical.radiation().energy(), owners.second.energy());
            assert_eq!(
                actual
                    .optical
                    .illumination()
                    .optical_field()
                    .wavelengths_nm(),
                &[450.0, 650.0]
            );
            assert_eq!(
                actual.optical.radiation().field().wavelengths_nm(),
                &[10_000.0]
            );
            near(
                actual.audit.external_heat_j,
                owners.first.energy().thermal_input_j + owners.second.energy().thermal_external_j,
            );
            assert_eq!(
                actual.optical.illumination().budget(),
                optical.illumination().budget()
            );
            assert_eq!(
                actual.optical.radiation().budget().max_source_evaluations,
                optical.radiation().budget().max_source_evaluations
            );
            final_temperatures.push(actual.thermal.temperature_k(CellId(11)).unwrap());
        }
        let (mean, se) = mean_se(&final_temperatures);
        near(result.cells[&CellId(11)].temperature_k.mean, mean);
        near(result.cells[&CellId(11)].temperature_k.standard_error, se);
        assert!(se > 0.0);
        assert_eq!(original.time_s(), 0.0);
        assert_eq!(optical.time_s(), 0.0);
        assert_eq!(optical.illumination().source_evaluations(), 0);
        assert!(optical.radiation().operator().is_none());
    });
}

struct CurvedJointCaloric;
impl ThermalLaw for CurvedJointCaloric {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        if h.is_finite() && (0.0..=800.0).contains(&h) {
            Ok(250.0 + 0.1 * h + 0.001 * h * h)
        } else {
            Err(DepositionError::MaterialDomain)
        }
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.7
    }
}
static CURVED_JOINT: CurvedJointCaloric = CurvedJointCaloric;

#[test]
fn actual_caloric_temperatures_and_joint_heat_covariance_survive_averaging() {
    with_cx(|cx| {
        let mut initial = joint_cell(11, 500.0);
        initial.specific_h_j_kg = 300.0;
        let original = thermal(vec![initial], &CURVED_JOINT, vec![]);
        let optical = joint_fixture(cx, &[1.0], &original, 80.0, 1.2, true);
        let mut c = study();
        c.replicas = 6;
        let result = optical
            .advance_ensemble_refined(&original, 0.5, c, || Ok(()))
            .unwrap();
        let replay = optical
            .advance_ensemble_refined(&original, 0.5, c, || Ok(()))
            .unwrap();
        assert_eq!(result.diagnostics, replay.diagnostics);
        assert_eq!(result.cells, replay.cells);
        assert_eq!(result.energy, replay.energy);
        let hs: Vec<_> = result
            .trajectories
            .iter()
            .map(|x| x.thermal.cells()[&CellId(11)].specific_h_j_kg)
            .collect();
        let ts: Vec<_> = result
            .trajectories
            .iter()
            .map(|x| x.thermal.temperature_k(CellId(11)).unwrap())
            .collect();
        let incident: Vec<_> = result
            .trajectories
            .iter()
            .map(|x| x.optical.illumination().energy().thermal_input_j)
            .collect();
        let radiation: Vec<_> = result
            .trajectories
            .iter()
            .map(|x| x.optical.radiation().energy().thermal_external_j)
            .collect();
        let net: Vec<_> = incident
            .iter()
            .zip(&radiation)
            .map(|(a, b)| a + b)
            .collect();
        let (h_mean, h_se) = mean_se(&hs);
        let (t_mean, t_se) = mean_se(&ts);
        let (net_mean, net_se) = mean_se(&net);
        near(result.cells[&CellId(11)].specific_h_j_kg.mean, h_mean);
        near(
            result.cells[&CellId(11)].specific_h_j_kg.standard_error,
            h_se,
        );
        near(result.cells[&CellId(11)].temperature_k.mean, t_mean);
        near(result.cells[&CellId(11)].temperature_k.standard_error, t_se);
        assert!(t_mean - CURVED_JOINT.temperature_k(h_mean).unwrap() > 1e-7);
        near(result.energy.net_optical_heat_j.mean, net_mean);
        near(result.energy.net_optical_heat_j.standard_error, net_se);
        let independent_se = (mean_se(&incident).1.powi(2) + mean_se(&radiation).1.powi(2)).sqrt();
        assert!((net_se - independent_se).abs() > 1e-9);
    });
}

#[test]
fn sampling_and_each_photon_source_have_independent_accuracy_gates() {
    with_cx(|cx| {
        let original = thermal(vec![joint_cell(11, 500.0)], &JOINT_CALORIC, vec![]);
        let optical = joint_fixture(cx, &[1.0], &original, 40.0, 1.5, false);
        let c = study();
        let base = optical
            .advance_ensemble_refined(&original, 0.5, c, || Ok(()))
            .unwrap();
        let incident_width = base
            .diagnostics
            .illumination_photon_sensitivity
            .unwrap()
            .maximum_illumination_energy_width_j;
        let radiation_width = base
            .diagnostics
            .radiation_photon_sensitivity
            .unwrap()
            .maximum_radiation_energy_width_j;
        let sampling_width = base
            .diagnostics
            .sampling
            .unwrap()
            .maximum_illumination_energy_width_j;
        assert!(incident_width > 0.0 && radiation_width > 0.0 && sampling_width > 0.0);
        let mut narrow = c;
        narrow
            .illumination_photon_sensitivity
            .illumination_energy_absolute_j = 0.25 * incident_width;
        match optical.advance_ensemble_refined(&original, 0.5, narrow, || Ok(())) {
            Err(JointOpticalRefinementError::Unresolved(d)) => {
                assert!(d.temporal.unwrap().maximum_normalized_width <= 1.0);
                assert!(d.sampling.unwrap().maximum_normalized_width <= 1.0);
                assert!(
                    d.illumination_photon_sensitivity
                        .unwrap()
                        .maximum_normalized_width
                        > 1.0
                );
                assert!(d.radiation_photon_sensitivity.is_none());
            }
            _ => panic!("tight illumination-only photon gate must refuse"),
        }
        let mut narrow = c;
        narrow
            .radiation_photon_sensitivity
            .radiation_energy_absolute_j = 0.25 * radiation_width;
        match optical.advance_ensemble_refined(&original, 0.5, narrow, || Ok(())) {
            Err(JointOpticalRefinementError::Unresolved(d)) => {
                assert!(
                    d.illumination_photon_sensitivity
                        .unwrap()
                        .maximum_normalized_width
                        <= 1.0
                );
                assert!(
                    d.radiation_photon_sensitivity
                        .unwrap()
                        .maximum_normalized_width
                        > 1.0
                );
            }
            _ => panic!("tight radiation-only photon gate must refuse"),
        }
        let mut narrow = c;
        narrow.sampling.illumination_energy_absolute_j = 0.25 * sampling_width;
        match optical.advance_ensemble_refined(&original, 0.5, narrow, || Ok(())) {
            Err(JointOpticalRefinementError::Unresolved(d)) => {
                assert!(d.sampling.unwrap().maximum_normalized_width > 1.0);
                assert!(d.illumination_photon_sensitivity.is_none());
                assert!(d.radiation_photon_sensitivity.is_none());
            }
            _ => panic!("tight sampling-only gate must refuse"),
        }
        assert_eq!(optical.time_s(), 0.0);
    });
}

#[test]
fn both_owners_charge_all_discarded_trials_and_conservative_storage() {
    with_cx(|cx| {
        let original = thermal(vec![joint_cell(11, 500.0)], &JOINT_CALORIC, vec![]);
        let optical = joint_fixture(cx, &[1.0], &original, 40.0, 0.8, false);
        let c = study();
        let base = optical
            .advance_ensemble_refined(&original, 0.5, c, || Ok(()))
            .unwrap();
        let work = base.diagnostics.work;
        // Coarse/fine/finer plus isolated low-illumination and low-radiation
        // ensembles. Constant coefficients need one build per complete trial.
        assert_eq!(work.illumination_transport_builds, 5 * c.replicas as usize);
        assert_eq!(work.radiation_operator_builds, 5 * c.replicas as usize);
        assert_eq!(work.illumination_source_evaluations, work.substeps);
        assert_eq!(work.radiation_source_evaluations, work.substeps);
        assert!(work.illumination_traced_segments > 0 && work.radiation_traced_segments > 0);
        for lane in 0..8 {
            let mut limited = c;
            match lane {
                0 => limited.maximum_total_substeps = work.substeps - 1,
                1 => {
                    limited.maximum_total_illumination_source_evaluations =
                        work.illumination_source_evaluations - 1
                }
                2 => {
                    limited.maximum_total_radiation_source_evaluations =
                        work.radiation_source_evaluations - 1
                }
                3 => {
                    limited.maximum_total_illumination_traced_segments =
                        work.illumination_traced_segments - 1
                }
                4 => {
                    limited.maximum_total_radiation_traced_segments =
                        work.radiation_traced_segments - 1
                }
                5 => {
                    limited.maximum_retained_operator_pair_slots =
                        base.diagnostics.admitted_peak_operator_pair_slots - 1
                }
                6 => {
                    limited.maximum_retained_operator_cell_slots =
                        base.diagnostics.admitted_peak_operator_cell_slots - 1
                }
                _ => {
                    limited.maximum_retained_illumination_report_slots =
                        base.diagnostics.admitted_peak_illumination_report_slots - 1
                }
            }
            assert!(
                optical
                    .advance_ensemble_refined(&original, 0.5, limited, || Ok(()))
                    .is_err()
            );
        }
        let mut polls = 0;
        assert!(
            optical
                .advance_ensemble_refined(&original, 0.5, c, || {
                    polls += 1;
                    if polls == 23 {
                        Err(DepositionError::Cancelled)
                    } else {
                        Ok(())
                    }
                })
                .is_err()
        );
        assert_eq!(original.time_s(), 0.0);
        assert_eq!(optical.time_s(), 0.0);
        assert_eq!(optical.illumination().source_evaluations(), 0);
        assert_eq!(optical.radiation().source_evaluations(), 0);
        assert_eq!(optical.illumination().energy().thermal_input_j, 0.0);
        assert_eq!(optical.radiation().energy().thermal_external_j, 0.0);
    });
}

#[test]
fn accepted_pure_absorption_joint_trajectories_continue_with_restored_cached_budgets() {
    with_cx(|cx| {
        let original = thermal(vec![joint_cell(11, 500.0)], &JOINT_CALORIC, vec![]);
        let optical = joint_fixture(cx, &[1.0], &original, 30.0, 0.0, false);
        let mut c = study();
        c.maximum_total_illumination_source_evaluations = 160;
        c.maximum_total_radiation_source_evaluations = 160;
        let result = optical
            .advance_ensemble_refined(&original, 1.0, c, || Ok(()))
            .unwrap();
        assert!(result.diagnostics.illumination_photon_sensitivity.is_none());
        assert!(result.diagnostics.radiation_photon_sensitivity.is_some());
        assert_eq!(result.diagnostics.illumination_histories_per_lane, 16);
        assert_eq!(result.diagnostics.radiation_histories_per_emitter_lane, 32);
        let accepted = result.trajectories.last().unwrap();
        let (light, radiation) = accepted.optical.clone().into_parts();
        // Cached response fields need no photons, but both native source clocks
        // must be free to exceed the former WHOLE-study source allowance.
        let cached = JointOpticalHeating::try_new(
            light.with_work_limits(5000, 0),
            radiation.with_work_limits(5000, 0),
            &accepted.thermal,
            || Ok(()),
        )
        .unwrap();
        let (continued, observed, audit) = accepted
            .thermal
            .advance_observed(
                &cached,
                3.0,
                StepConfig {
                    max_dt_s: 0.0051,
                    energy_tolerance_j: c.conservation_tolerance_j,
                },
                || Ok(()),
            )
            .unwrap();
        assert!(observed.illumination().source_evaluations() > 160);
        assert!(observed.radiation().source_evaluations() > 160);
        assert_eq!(observed.illumination().transport_builds(), 0);
        assert_eq!(observed.radiation().operator_builds(), 0);
        assert_eq!(observed.illumination().traced_segments(), 0);
        assert_eq!(observed.radiation().traced_segments(), 0);
        assert_eq!(continued.time_s(), 3.0);
        near(
            audit.external_heat_j,
            observed.illumination().energy().thermal_input_j
                - accepted.optical.illumination().energy().thermal_input_j
                + observed.radiation().energy().thermal_external_j
                - accepted.optical.radiation().energy().thermal_external_j,
        );
        near(audit.stored_energy_change_j, audit.external_heat_j);
    });
}
