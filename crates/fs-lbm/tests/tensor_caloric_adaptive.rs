#![cfg(feature = "enthalpy-coupling")]
//! Embedded original-energy accuracy and accepted physical heat in native tensors.

use fs_blake3::ContentHash;
use fs_exec::CancelGate;
use fs_lbm::conjugate::{
    ChtError, ConductivityTensor, EnergyConfig, SolidMaterial, TensorCaloricAdaptiveConfig,
    TensorCaloricAdaptiveError, TensorCaloricForcing, TensorCaloricLaw, TensorCaloricStorage,
    TensorNonlinearConfig, ThermalFace, ThermalSetup, Voxel, VoxelDomain,
    advance_caloric_tensor_conduction_adaptive,
};
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve};

fn config(absolute: f64) -> TensorCaloricAdaptiveConfig {
    TensorCaloricAdaptiveConfig {
        step: TensorNonlinearConfig {
            linear: EnergyConfig {
                tolerance: 1e-10,
                ..EnergyConfig::default()
            },
            ..TensorNonlinearConfig::default()
        },
        initial_dt_s: 1.0,
        minimum_dt_s: 1e-8,
        maximum_dt_s: 1.0,
        maximum_attempts: 4096,
        energy_absolute_tolerance_j_kg: absolute,
        relative_tolerance: 0.0,
        energy_tolerance_j: 1e-6,
    }
}

fn tensor(k: f64) -> ConductivityTensor {
    ConductivityTensor::new([[k, 0.0, 0.0], [0.0, k, 0.0], [0.0, 0.0, k]]).unwrap()
}

fn sensible() -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_new(
        ContentHash([0xB1; 32]),
        [(-100.0, 200.0), (100.0, 400.0)]
            .map(|(e, t)| EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: e,
                temperature_k: t,
                liquid_mass_fraction: 0.0,
                bulk_density_kg_m3: 1.0,
            })
            .to_vec(),
    )
    .unwrap()
}

fn latent(datum: f64) -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_new(
        ContentHash([0xB2; 32]),
        [
            (-100.0, 250.0, 0.0),
            (0.0, 300.0, 0.0),
            (100.0, 300.0, 1.0),
            (300.0, 400.0, 1.0),
        ]
        .map(|(e, t, fraction)| EnthalpyPhaseKnot {
            specific_enthalpy_j_kg: datum + e,
            temperature_k: t,
            liquid_mass_fraction: fraction,
            bulk_density_kg_m3: 1.0,
        })
        .to_vec(),
    )
    .unwrap()
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} != {expected}"
    );
}

#[test]
fn tighter_energy_accuracy_refines_a_native_diffusion_mode_with_two_stages_per_attempt() {
    // Six unit faces, each with half-cell resistance (1/2)/k:
    // k=1/12, unit mass and Cp=1 give e'=-e exactly.
    let domain = VoxelDomain::from_fn(1, 1, 1, 1.0, |_| Voxel::Solid(0)).unwrap();
    let solids = [SolidMaterial::new("unit decay mode", 1.0 / 12.0)];
    let conductivity = tensor(1.0 / 12.0);
    let law = sensible();
    let laws: [&dyn TensorCaloricLaw; 1] = [&law];
    let setup = ThermalSetup::new([ThermalFace::Temperature(300.0); 6]);
    let initial = [1.0];
    let mut errors = Vec::new();
    let mut accepted = Vec::new();
    for tolerance in [1e-3, 1e-5] {
        let mut calls = 0;
        let result = advance_caloric_tensor_conduction_adaptive(
            &domain,
            &solids,
            &[&conductivity],
            TensorCaloricStorage {
                laws: &laws,
                reference_mass_kg: &[1.0],
                previous_specific_energy_j_kg: &initial,
            },
            0.0,
            1.0,
            &[],
            &config(tolerance),
            |sample| {
                calls += 1;
                assert!(sample.time_s > sample.interval_s[0]);
                assert!(sample.time_s <= sample.interval_s[1]);
                assert!(!sample.left_limit);
                Ok(TensorCaloricForcing::borrowed(&setup, &[]))
            },
            &CancelGate::new(),
        )
        .unwrap();
        assert_eq!(calls, 2 * result.report.attempts);
        assert_eq!(result.report.accepted_steps, result.records.len());
        assert!(result.report.accuracy_rejections > 0);
        assert_eq!(result.report.solver_rejections, 0);
        assert!(result.report.maximum_accepted_normalized_error <= 1.0);
        assert_eq!(result.record.time_s, 1.0);
        assert!(result.records.windows(2).all(|w| w[0].time_s < w[1].time_s));
        let energy = result.endpoint.specific_energy_j_kg[0];
        close(
            result.endpoint.solution.solution.temperature[0],
            300.0 + energy,
            1e-10,
        );
        close(result.record.stored_energy_change_j, energy - 1.0, 1e-12);
        close(result.record.source_j, 0.0, 0.0);
        close(result.record.boundary_outflow_j, 1.0 - energy, 1e-8);
        close(result.record.closure_j, 0.0, 1e-8);
        assert!(result.report.absolute_accepted_closure_j <= 1e-6);
        close(
            result.records.iter().map(|r| r.source_j).sum(),
            result.record.source_j,
            1e-12,
        );
        errors.push((energy - (-1.0_f64).exp()).abs());
        accepted.push(result.report.accepted_steps);
    }
    assert!(errors[1] < errors[0] / 5.0, "{errors:?}");
    assert!(accepted[1] > accepted[0]);
    assert_eq!(initial, [1.0]);
}

#[test]
fn latent_energy_controls_time_accuracy_even_when_every_temperature_is_constant() {
    let domain = VoxelDomain::from_fn(1, 1, 1, 1.0, |_| Voxel::Solid(0)).unwrap();
    let solids = [SolidMaterial::new("isothermal latent inventory", 1.0)];
    let conductivity = tensor(1.0);
    let mut changes = Vec::new();
    for datum in [0.0, 1e9] {
        let law = latent(datum);
        let laws: [&dyn TensorCaloricLaw; 1] = [&law];
        let initial = [datum + 20.0];
        let mut controls = config(1e-2);
        // The 1e9 datum quantizes returned energies. The spatial tolerance
        // admits that representable heat while remaining well below the
        // requested 1e-2 J/kg temporal accuracy; no datum allowance is added.
        controls.step.linear.tolerance = 1e-4;
        controls.energy_tolerance_j = 1e-4;
        let mut calls = 0;
        let result = advance_caloric_tensor_conduction_adaptive(
            &domain,
            &solids,
            &[&conductivity],
            TensorCaloricStorage {
                laws: &laws,
                reference_mass_kg: &[1.0],
                previous_specific_energy_j_kg: &initial,
            },
            0.0,
            1.0,
            &[],
            &controls,
            |sample| {
                calls += 1;
                let mut setup = ThermalSetup::new([ThermalFace::Adiabatic; 6]);
                setup.power_w = vec![20.0 * sample.time_s];
                Ok(TensorCaloricForcing::owned(setup, vec![]))
            },
            &CancelGate::new(),
        )
        .unwrap();
        assert_eq!(calls, 2 * result.report.attempts);
        assert!(result.report.accuracy_rejections > 0);
        assert!(result.report.accepted_steps > 1);
        assert!(result.records.iter().all(|r| r.max_temperature_k == 300.0));
        assert_eq!(result.endpoint.solution.solution.temperature, vec![300.0]);
        close(result.record.source_j, 10.0, 1e-10); // integral_0^1 20t dt
        close(result.record.boundary_outflow_j, 0.0, 1e-10);
        close(result.endpoint.specific_energy_j_kg[0], datum + 30.0, 2e-5);
        close(
            law.values_at_specific_enthalpy(result.endpoint.specific_energy_j_kg[0])
                .unwrap()
                .liquid_mass_fraction,
            0.3,
            2e-7,
        );
        changes.push(result.endpoint.specific_energy_j_kg[0] - initial[0]);
        assert_eq!(initial, [datum + 20.0]);
    }
    close(changes[0], changes[1], 2e-5);
}

#[test]
fn declared_transitions_use_left_limits_and_allow_a_new_prescribed_mask() {
    let domain = VoxelDomain::from_fn(1, 1, 1, 1.0, |_| Voxel::Solid(0)).unwrap();
    let solids = [SolidMaterial::new("piecewise forcing", 1.0)];
    let conductivity = tensor(1.0);
    let law = sensible();
    let laws: [&dyn TensorCaloricLaw; 1] = [&law];
    let initial = [0.0];
    for prescribed_after_event in [false, true] {
        let mut left_limits = 0;
        let result = advance_caloric_tensor_conduction_adaptive(
            &domain,
            &solids,
            &[&conductivity],
            TensorCaloricStorage {
                laws: &laws,
                reference_mass_kg: &[1.0],
                previous_specific_energy_j_kg: &initial,
            },
            0.0,
            1.0,
            &[0.5],
            &config(1e-8),
            |sample| {
                let left_branch =
                    sample.time_s < 0.5 || (sample.time_s == 0.5 && sample.left_limit);
                if sample.left_limit {
                    left_limits += 1;
                    assert_eq!(sample.time_s, 0.5);
                    assert_eq!(sample.interval_s, [0.0, 0.5]);
                }
                assert!(sample.interval_s[1] <= 0.5 || sample.interval_s[0] >= 0.5);
                let mut setup = ThermalSetup::new([ThermalFace::Adiabatic; 6]);
                setup.power_w = vec![if left_branch { 2.0 } else { 5.0 }];
                let prescribed = if prescribed_after_event && !left_branch {
                    setup.fixed_temperature = vec![(0, 307.0)];
                    vec![(0, 7.0)]
                } else {
                    vec![]
                };
                Ok(TensorCaloricForcing::owned(setup, prescribed))
            },
            &CancelGate::new(),
        )
        .unwrap();
        assert_eq!(left_limits, 1);
        assert_eq!(result.records.len(), 2);
        assert_eq!(result.records[0].time_s, 0.5);
        assert_eq!(result.records[1].time_s, 1.0);
        let final_energy = if prescribed_after_event { 7.0 } else { 3.5 };
        close(result.endpoint.specific_energy_j_kg[0], final_energy, 1e-10);
        close(result.record.stored_energy_change_j, final_energy, 1e-10);
        close(result.record.closure_j, 0.0, 1e-10);
        if prescribed_after_event {
            close(result.record.source_j, 1.0, 1e-10);
            close(result.record.boundary_outflow_j, -6.0, 1e-10);
        } else {
            close(result.record.source_j, 3.5, 1e-10);
            close(result.record.boundary_outflow_j, 0.0, 1e-10);
        }
    }
    assert_eq!(initial, [0.0]);
}

#[test]
fn attempt_and_accuracy_limits_refuse_without_publishing_a_thermal_prefix() {
    let domain = VoxelDomain::from_fn(1, 1, 1, 1.0, |_| Voxel::Solid(0)).unwrap();
    let solids = [SolidMaterial::new("bounded decay", 1.0 / 12.0)];
    let conductivity = tensor(1.0 / 12.0);
    let law = sensible();
    let laws: [&dyn TensorCaloricLaw; 1] = [&law];
    let setup = ThermalSetup::new([ThermalFace::Temperature(300.0); 6]);
    let initial = [1.0];
    let solve = |controls: TensorCaloricAdaptiveConfig| {
        advance_caloric_tensor_conduction_adaptive(
            &domain,
            &solids,
            &[&conductivity],
            TensorCaloricStorage {
                laws: &laws,
                reference_mass_kg: &[1.0],
                previous_specific_energy_j_kg: &initial,
            },
            0.0,
            1.0,
            &[],
            &controls,
            |_| Ok(TensorCaloricForcing::borrowed(&setup, &[])),
            &CancelGate::new(),
        )
    };
    let mut controls = config(1e-8);
    controls.maximum_attempts = 1;
    assert!(matches!(
        solve(controls),
        Err(TensorCaloricAdaptiveError::AttemptBudget {
            attempts: 1,
            time_s: 0.0
        })
    ));
    controls.minimum_dt_s = 1.0;
    assert!(matches!(
        solve(controls),
        Err(TensorCaloricAdaptiveError::TimeAccuracy { dt_s: 1.0, .. })
    ));
    assert_eq!(initial, [1.0]);
}

#[test]
fn undeclared_mask_changes_source_refusal_and_late_cancellation_are_fatal() {
    let domain = VoxelDomain::from_fn(1, 1, 1, 1.0, |_| Voxel::Solid(0)).unwrap();
    let solids = [SolidMaterial::new("atomic schedule", 1.0)];
    let conductivity = tensor(1.0);
    let law = sensible();
    let laws: [&dyn TensorCaloricLaw; 1] = [&law];
    let initial = [1.0];
    let storage = TensorCaloricStorage {
        laws: &laws,
        reference_mass_kg: &[1.0],
        previous_specific_energy_j_kg: &initial,
    };
    for between_attempts in [false, true] {
        let mut controls = config(1e-4);
        controls.initial_dt_s = 0.25;
        controls.maximum_dt_s = 0.25;
        let bad_mask = advance_caloric_tensor_conduction_adaptive(
            &domain,
            &solids,
            &[&conductivity],
            storage,
            0.0,
            1.0,
            &[],
            &controls,
            |sample| {
                let mut setup = ThermalSetup::new([ThermalFace::Adiabatic; 6]);
                let prescribed_now = if between_attempts {
                    sample.time_s > 0.25
                } else {
                    sample.time_s == sample.interval_s[1]
                };
                let prescribed = if prescribed_now {
                    setup.fixed_temperature = vec![(0, 307.0)];
                    vec![(0, 7.0)]
                } else {
                    vec![]
                };
                Ok(TensorCaloricForcing::owned(setup, prescribed))
            },
            &CancelGate::new(),
        );
        assert!(matches!(
            bad_mask,
            Err(TensorCaloricAdaptiveError::Invalid(_))
        ));
    }
    for cancel in [false, true] {
        let mut controls = config(1e-4);
        controls.initial_dt_s = 0.25;
        controls.maximum_dt_s = 0.25;
        let gate = CancelGate::new();
        let mut calls = 0;
        let result = advance_caloric_tensor_conduction_adaptive(
            &domain,
            &solids,
            &[&conductivity],
            storage,
            0.0,
            1.0,
            &[],
            &controls,
            |sample| {
                calls += 1;
                if sample.interval_s[0] >= 0.25 {
                    if cancel {
                        gate.request();
                    } else {
                        return Err(ChtError::InvalidInput {
                            field: "test.original_source",
                            reason: "source support ended".into(),
                        });
                    }
                }
                let mut setup = ThermalSetup::new([ThermalFace::Adiabatic; 6]);
                setup.power_w = vec![1.0];
                Ok(TensorCaloricForcing::owned(setup, vec![]))
            },
            &gate,
        );
        if cancel {
            assert!(matches!(
                result,
                Err(TensorCaloricAdaptiveError::Solve(ChtError::Cancelled))
            ));
        } else {
            assert!(matches!(
                result,
                Err(TensorCaloricAdaptiveError::Solve(
                    ChtError::InvalidInput { .. }
                ))
            ));
        }
        // The first interval was privately accepted. The next source query is
        // fatal and is never turned into a timestep retry or a partial result.
        assert_eq!(calls, 3);
        assert_eq!(initial, [1.0]);
    }
}
