#![cfg(feature = "enthalpy-coupling")]
//! G1/G3/G4: original adaptive energy, scheduled current and accepted-only heat.

use fs_blake3::ContentHash;
use fs_exec::CancelGate;
use fs_lbm::conjugate::{
    CaloricElectrothermalAdaptiveConfig, CaloricElectrothermalConfig, CaloricElectrothermalForcing,
    ChtError, ConductivityTensor, ElectricSetup, Electrode, ElectrodeDrive, EnergyConfig,
    RadiationConfig, STEFAN_BOLTZMANN, SolidMaterial, TensorCaloricAdaptiveError,
    TensorCaloricForcingSample, TensorCaloricLaw, TensorCaloricStorage, TensorConductivityLaw,
    TensorNonlinearConfig, ThermalFace, ThermalSetup, Voxel, VoxelDomain,
    advance_caloric_tensor_electrothermal_adaptive,
};
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, SolidLiquidPhase};
use std::{borrow::Cow, cell::Cell};

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}"
    );
}

fn phase(datum: f64, latent: bool) -> EquilibriumEnthalpyPhaseCurve {
    let points = if latent {
        vec![
            (-200.0, 200.0, 0.0),
            (0.0, 300.0, 0.0),
            (100.0, 300.0, 1.0),
            (500.0, 500.0, 1.0),
        ]
    } else {
        vec![(-200.0, 200.0, 0.0), (400.0, 500.0, 0.0)]
    };
    let knots = points
        .into_iter()
        .map(
            |(energy, temperature_k, liquid_mass_fraction)| EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: datum + energy,
                temperature_k,
                liquid_mass_fraction,
                bulk_density_kg_m3: 2.0,
            },
        )
        .collect();
    if latent {
        EquilibriumEnthalpyPhaseCurve::try_new(ContentHash([0xEA; 32]), knots)
    } else {
        EquilibriumEnthalpyPhaseCurve::try_single_phase(
            ContentHash([0xEA; 32]),
            SolidLiquidPhase::Solid,
            knots,
        )
    }
    .unwrap()
}

fn bar(
    exposed: bool,
) -> (
    VoxelDomain,
    [SolidMaterial; 1],
    ElectricSetup,
    ConductivityTensor,
) {
    let domain = VoxelDomain::from_fn(3, if exposed { 2 } else { 1 }, 1, 1.0, |point| {
        if point[1] < 1.0 {
            Voxel::Solid(0)
        } else {
            Voxel::Fluid
        }
    })
    .unwrap();
    let electric = ElectricSetup {
        conductors: vec![],
        electrodes: vec![
            Electrode {
                lo: [0, 0, 0],
                hi: [1, 1, 1],
                drive: ElectrodeDrive::Current(4.0),
            },
            Electrode {
                lo: [2, 0, 0],
                hi: [3, 1, 1],
                drive: ElectrodeDrive::Voltage(0.0),
            },
        ],
    };
    let tensor =
        ConductivityTensor::new([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]).unwrap();
    (
        domain,
        [SolidMaterial::new("manufactured caloric conductor", 1.0)],
        electric,
        tensor,
    )
}

fn controls(energy_tolerance: f64) -> CaloricElectrothermalAdaptiveConfig {
    CaloricElectrothermalAdaptiveConfig {
        step: CaloricElectrothermalConfig {
            thermal: TensorNonlinearConfig {
                linear: EnergyConfig {
                    tolerance: 1e-11,
                    ..EnergyConfig::default()
                },
                max_newton_iterations: 256,
                max_backtracks: 256,
            },
            maximum_passes: 64,
            specific_energy_tolerance_j_kg: 1e-8,
            relative_tolerance: 0.0,
            power_tolerance_w: 1e-8,
            energy_tolerance_j: 1e-6,
        },
        initial_dt_s: 1.0,
        minimum_dt_s: 1e-6,
        maximum_dt_s: 1.0,
        maximum_attempts: 8192,
        energy_absolute_tolerance_j_kg: energy_tolerance,
        source_absolute_tolerance_j: 2.0 * energy_tolerance,
        relative_tolerance: 0.0,
        energy_tolerance_j: 1e-4,
    }
}

#[test]
fn adaptive_latent_feedback_matches_exponential_and_retains_only_accepted_heat() {
    let (domain, solids, electric, k) = bar(false);
    let setup = ThermalSetup::new([ThermalFace::Adiabatic; 6]);
    let expected = 120.0 * (0.16_f64).exp() - 100.0;
    let mut answers = Vec::new();
    for (datum, tolerance) in [(0.0, 1e-2), (0.0, 1e-3), (1e6, 1e-3)] {
        let law = phase(datum, true);
        let laws: [&dyn TensorCaloricLaw; 3] = [&law; 3];
        let previous = [datum + 20.0; 3];
        let config = controls(tolerance);
        let result = advance_caloric_tensor_electrothermal_adaptive(
            &domain,
            &solids,
            &[&k as &dyn TensorConductivityLaw; 3],
            TensorCaloricStorage {
                laws: &laws,
                reference_mass_kg: &[2.0; 3],
                previous_specific_energy_j_kg: &previous,
            },
            0.0,
            2.0,
            &[],
            &config,
            None,
            |_| {
                Ok(CaloricElectrothermalForcing::borrowed(
                    &setup,
                    &[],
                    &electric,
                    &[],
                ))
            },
            |_, _, _, energy, temperature| {
                assert_eq!(temperature, 300.0);
                Ok(Some(1.0 + 0.01 * (energy - datum)))
            },
            &CancelGate::new(),
        )
        .unwrap();
        let endpoint = result.thermal.endpoint.specific_energy_j_kg[1] - datum;
        answers.push(endpoint);
        assert_eq!(result.thermal.record.time_s, 2.0);
        assert_eq!(
            result.thermal.endpoint.solution.solution.temperature,
            vec![300.0; 3]
        );
        assert!(result.thermal.report.accuracy_rejections > 0);
        assert!(result.thermal.report.maximum_accepted_normalized_error <= 1.0);
        assert_eq!(result.coupling.len(), result.thermal.records.len());
        assert!(result.radiation.is_empty());
        assert_eq!(result.radiated_energy_j, 0.0);
        let accepted_passes: usize = result.coupling.iter().map(|r| r.passes).sum();
        assert!(result.report.completed_trial_passes > accepted_passes);
        assert_eq!(
            result.report.completed_trial_current_solves,
            2 * (result.report.completed_trial_passes
                + result.thermal.report.accepted_steps
                + result.thermal.report.accuracy_rejections)
        );
        let accepted_native_iterations: usize =
            result.thermal.records.iter().map(|r| r.iterations).sum();
        // The exact latent plateau has a zero condensed face RHS, so native
        // Krylov work may be zero even while source passes are rejected.
        assert!(result.thermal.report.completed_trial_iterations >= accepted_native_iterations);
        near(
            result.report.applied_electrode_energy_j,
            2.0 * (endpoint - 20.0),
            config.energy_tolerance_j,
        );
        near(
            result.report.applied_joule_energy_j,
            result.thermal.record.source_j,
            1e-8,
        );
        near(result.report.thermostat_joule_energy_j, 0.0, 1e-12);
        near(
            result.endpoint_electric[1].joule_w[1],
            16.0 * (1.0 + 0.01 * endpoint),
            1e-8,
        );
        assert!(result.report.absolute_energy_defect_j <= config.energy_tolerance_j);
        assert_eq!(previous, [datum + 20.0; 3]);
    }
    assert!((answers[1] - expected).abs() < 0.4 * (answers[0] - expected).abs());
    near(answers[1], expected, 2e-3);
    near(answers[2], answers[1], 5e-6);

    let law = phase(0.0, true);
    let laws: [&dyn TensorCaloricLaw; 3] = [&law; 3];
    let mut limited = controls(1e-3);
    limited.maximum_attempts = 1;
    let refused = advance_caloric_tensor_electrothermal_adaptive(
        &domain,
        &solids,
        &[&k as &dyn TensorConductivityLaw; 3],
        TensorCaloricStorage {
            laws: &laws,
            reference_mass_kg: &[2.0; 3],
            previous_specific_energy_j_kg: &[20.0; 3],
        },
        0.0,
        2.0,
        &[],
        &limited,
        None,
        |_| {
            Ok(CaloricElectrothermalForcing::borrowed(
                &setup,
                &[],
                &electric,
                &[],
            ))
        },
        |_, _, _, energy, _| Ok(Some(1.0 + 0.01 * energy)),
        &CancelGate::new(),
    );
    assert!(matches!(
        refused,
        Err(TensorCaloricAdaptiveError::AttemptBudget { attempts: 1, .. })
    ));
}

#[test]
fn thermostat_heat_has_its_own_temporal_accuracy_when_storage_is_constant() {
    let (domain, solids, electric, k) = bar(false);
    let law = phase(0.0, true);
    let laws: [&dyn TensorCaloricLaw; 3] = [&law; 3];
    let previous = [20.0; 3];
    let mut setup = ThermalSetup::new([ThermalFace::Adiabatic; 6]);
    setup.fixed_temperature = vec![(0, 300.0), (1, 300.0), (2, 300.0)];
    let prescribed = [(0, 20.0), (1, 20.0), (2, 20.0)];
    let mut errors = Vec::new();
    for tolerance in [0.02, 0.002] {
        let mut config = controls(1.0);
        config.source_absolute_tolerance_j = tolerance;
        let result = advance_caloric_tensor_electrothermal_adaptive(
            &domain,
            &solids,
            &[&k as &dyn TensorConductivityLaw; 3],
            TensorCaloricStorage {
                laws: &laws,
                reference_mass_kg: &[2.0; 3],
                previous_specific_energy_j_kg: &previous,
            },
            0.0,
            1.0,
            &[],
            &config,
            None,
            |sample| {
                let mut drive = electric.clone();
                drive.electrodes[0].drive = ElectrodeDrive::Current(2.0 + sample.time_s);
                let mut forcing =
                    CaloricElectrothermalForcing::borrowed(&setup, &prescribed, &electric, &[]);
                forcing.electric = Cow::Owned(drive);
                Ok(forcing)
            },
            |_, _, _, _, _| Ok(Some(1.0)),
            &CancelGate::new(),
        )
        .unwrap();
        assert_eq!(result.thermal.endpoint.specific_energy_j_kg, previous);
        near(result.thermal.record.source_j, 0.0, 1e-12);
        near(result.thermal.record.stored_energy_change_j, 0.0, 1e-12);
        near(result.thermal.record.boundary_outflow_j, 0.0, 1e-12);
        assert!(result.thermal.report.accuracy_rejections > 0);
        assert!(result.thermal.records.len() > 2);
        near(
            result.report.thermostat_joule_energy_j,
            result.report.applied_joule_energy_j,
            1e-10,
        );
        near(
            result.report.applied_electrode_energy_j,
            result.report.applied_joule_energy_j,
            1e-10,
        );
        // Analytic integral of the actual prescribed I(t)^2 R=(2+t)^2 W.
        errors.push((result.report.applied_joule_energy_j - 19.0 / 3.0).abs());
    }
    assert!(errors[1] < 0.4 * errors[0]);
    assert!(errors[1] < 2e-4);
}

#[test]
fn scheduled_current_and_resistivity_share_left_limits_and_atomic_refusals() {
    let (domain, solids, electric, k) = bar(false);
    let setup = ThermalSetup::new([ThermalFace::Adiabatic; 6]);
    let law = phase(0.0, true);
    let laws: [&dyn TensorCaloricLaw; 3] = [&law; 3];
    let previous = [20.0; 3];
    let before = |sample: TensorCaloricForcingSample| {
        sample.time_s < 0.5 || (sample.time_s == 0.5 && sample.left_limit)
    };
    for mode in 0..3 {
        let force_left = Cell::new(false);
        let resistance_left = Cell::new(false);
        let saw_later = Cell::new(false);
        let gate = CancelGate::new();
        let result = advance_caloric_tensor_electrothermal_adaptive(
            &domain,
            &solids,
            &[&k as &dyn TensorConductivityLaw; 3],
            TensorCaloricStorage {
                laws: &laws,
                reference_mass_kg: &[2.0; 3],
                previous_specific_energy_j_kg: &previous,
            },
            0.0,
            1.0,
            &[0.5],
            &controls(1e-3),
            None,
            |sample| {
                if sample.time_s == 0.5 {
                    force_left.set(sample.left_limit);
                }
                let mut drive = electric.clone();
                drive.electrodes[0].drive =
                    ElectrodeDrive::Current(if before(sample) { 4.0 } else { 2.0 });
                let mut forcing =
                    CaloricElectrothermalForcing::borrowed(&setup, &[], &electric, &[]);
                forcing.electric = Cow::Owned(drive);
                Ok(forcing)
            },
            |sample, _, _, _, _| {
                if sample.time_s == 0.5 {
                    resistance_left.set(sample.left_limit);
                }
                if sample.interval_s[0] == 0.5 {
                    saw_later.set(true);
                    if mode == 1 {
                        return Err(ChtError::InvalidInput {
                            field: "test.adaptive_electrical_source",
                            reason: "later source unavailable".into(),
                        });
                    }
                    if mode == 2 {
                        gate.request();
                    }
                }
                Ok(Some(if before(sample) { 1.0 } else { 2.0 }))
            },
            &gate,
        );
        assert!(force_left.get() && resistance_left.get() && saw_later.get());
        match mode {
            0 => {
                let result = result.unwrap();
                assert_eq!(
                    result
                        .thermal
                        .records
                        .iter()
                        .map(|r| r.time_s)
                        .collect::<Vec<_>>(),
                    vec![0.5, 1.0]
                );
                near(result.thermal.endpoint.specific_energy_j_kg[1], 26.0, 1e-8);
                near(result.report.applied_electrode_energy_j, 12.0, 1e-8);
                near(result.thermal.record.stored_energy_change_j, 12.0, 1e-8);
            }
            1 => assert!(matches!(
                result,
                Err(TensorCaloricAdaptiveError::Solve(ChtError::InvalidInput {
                    field: "test.adaptive_electrical_source",
                    ..
                }))
            )),
            _ => assert!(matches!(
                result,
                Err(TensorCaloricAdaptiveError::Solve(ChtError::Cancelled))
            )),
        }
        assert_eq!(previous, [20.0; 3]);
    }
}

#[test]
fn adaptive_radiating_current_matches_independent_nonlinear_cooling() {
    let (domain, solids, electric, k) = bar(true);
    let law = phase(0.0, false);
    let laws: [&dyn TensorCaloricLaw; 6] = [&law; 6];
    let previous = [0.0, 200.0, 0.0, 0.0, 0.0, 0.0];
    let mut setup = ThermalSetup::new([ThermalFace::Adiabatic; 6]);
    setup.fixed_temperature = [(0, 300.0), (2, 300.0), (3, 300.0), (4, 300.0), (5, 300.0)].to_vec();
    let prescribed = [(0, 0.0), (2, 0.0), (3, 0.0), (4, 0.0), (5, 0.0)];
    let mut radiation = RadiationConfig::new(vec![0.04], [Some(250.0); 6]);
    radiation.rays_per_face = 16;
    radiation.patch_size = 1;
    radiation.max_iterations = 64;
    radiation.tolerance = 1e-10;
    // Independent scalar ODE: three unit conduction paths and one emitting
    // face meet the free cell. Its M*Cp is 4 J/K; no native solver is reused.
    let rate = |t: f64| {
        (16.0 * (1.0 + 0.002 * (t - 300.0))
            - 3.0 * (t - 300.0)
            - 0.04 * STEFAN_BOLTZMANN * (t.powi(4) - 250.0_f64.powi(4)))
            / 4.0
    };
    let mut reference = 400.0;
    let dt = 0.5 / 8192.0;
    for _ in 0..8192 {
        let k1 = rate(reference);
        let k2 = rate(reference + 0.5 * dt * k1);
        let k3 = rate(reference + 0.5 * dt * k2);
        let k4 = rate(reference + dt * k3);
        reference += dt * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
    }
    let mut errors = Vec::new();
    for tolerance in [0.02, 0.002] {
        let config = controls(tolerance);
        let result = advance_caloric_tensor_electrothermal_adaptive(
            &domain,
            &solids,
            &[&k as &dyn TensorConductivityLaw; 6],
            TensorCaloricStorage {
                laws: &laws,
                reference_mass_kg: &[2.0; 6],
                previous_specific_energy_j_kg: &previous,
            },
            0.0,
            0.5,
            &[],
            &config,
            Some(&radiation),
            |_| {
                Ok(CaloricElectrothermalForcing::borrowed(
                    &setup,
                    &prescribed,
                    &electric,
                    &[],
                ))
            },
            |_, _, occupancy, energy, _| {
                Ok(match occupancy {
                    Voxel::Fluid => None,
                    Voxel::Solid(_) => Some(1.0 + 0.001 * energy),
                })
            },
            &CancelGate::new(),
        )
        .unwrap();
        errors.push((result.thermal.endpoint.solution.solution.temperature[1] - reference).abs());
        assert_eq!(result.thermal.record.time_s, 0.5);
        assert_eq!(result.radiation.len(), result.thermal.records.len());
        assert!(result.thermal.report.accuracy_rejections > 0);
        let accepted_solves: usize = result
            .radiation
            .iter()
            .flatten()
            .map(|r| r.iterations)
            .sum();
        assert!(result.completed_trial_radiation_solves > accepted_solves);
        assert!(result.radiated_energy_j > 0.0);
        assert!(result.radiation.iter().flatten().all(|r| r.rays == 48));
        near(
            result.report.applied_electrode_energy_j,
            result.thermal.record.source_j,
            config.energy_tolerance_j,
        );
        near(
            result.thermal.record.stored_energy_change_j,
            result.thermal.record.source_j - result.thermal.record.boundary_outflow_j,
            config.energy_tolerance_j,
        );
        assert_eq!(previous, [0.0, 200.0, 0.0, 0.0, 0.0, 0.0]);
    }
    assert!(errors[1] < 0.4 * errors[0]);
    assert!(errors[1] < 0.005);
}
