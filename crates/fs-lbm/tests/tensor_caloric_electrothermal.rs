#![cfg(feature = "enthalpy-coupling")]
//! G1/G3/G4: native current heats original latent energy on the tensor grid.
use fs_blake3::ContentHash;
use fs_exec::CancelGate;
use fs_lbm::conjugate::{
    CaloricElectrothermalConfig, ChtError, ConductivityTensor, ElectricSetup, ElectricalContact,
    Electrode, ElectrodeDrive, EnergyConfig, RadiationConfig, STEFAN_BOLTZMANN, SolidMaterial,
    TensorCaloricLaw, TensorCaloricStorage, TensorNonlinearConfig, ThermalFace, ThermalSetup,
    Voxel, VoxelDomain, advance_caloric_tensor_electrothermal,
    advance_caloric_tensor_electrothermal_radiating,
};
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve};

const CURRENT: f64 = 4.0;
const MASS: f64 = 2.0;
const DURATION: f64 = 1.0;
const OPENING_H: f64 = 20.0;
const RESISTIVITY_SLOPE: f64 = 0.01;

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}"
    );
}

fn phase(datum: f64) -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_new(
        ContentHash([if datum == 0.0 { 0xE1 } else { 0xE2 }; 32]),
        [
            (-200.0, 200.0, 0.0),
            (0.0, 300.0, 0.0),
            (100.0, 300.0, 1.0),
            (500.0, 500.0, 1.0),
        ]
        .map(
            |(h, temperature_k, liquid_mass_fraction)| EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: datum + h,
                temperature_k,
                liquid_mass_fraction,
                bulk_density_kg_m3: MASS,
            },
        )
        .to_vec(),
    )
    .unwrap()
}

fn bar() -> (
    VoxelDomain,
    [SolidMaterial; 1],
    ElectricSetup,
    ConductivityTensor,
) {
    let domain = VoxelDomain::from_fn(3, 1, 1, 1.0, |_| Voxel::Solid(0)).unwrap();
    let electric = ElectricSetup {
        conductors: vec![],
        electrodes: vec![
            Electrode {
                lo: [0, 0, 0],
                hi: [1, 1, 1],
                drive: ElectrodeDrive::Current(CURRENT),
            },
            Electrode {
                lo: [2, 0, 0],
                hi: [3, 1, 1],
                drive: ElectrodeDrive::Voltage(0.0),
            },
        ],
    };
    let k = ConductivityTensor::new([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]).unwrap();
    (
        domain,
        [SolidMaterial::new("manufactured latent conductor", 1.0)],
        electric,
        k,
    )
}

fn controls() -> CaloricElectrothermalConfig {
    CaloricElectrothermalConfig {
        thermal: TensorNonlinearConfig {
            linear: EnergyConfig {
                tolerance: 1e-10,
                ..EnergyConfig::default()
            },
            max_newton_iterations: 128,
            max_backtracks: 128,
        },
        maximum_passes: 64,
        specific_energy_tolerance_j_kg: 1e-8,
        relative_tolerance: 0.0,
        power_tolerance_w: 1e-8,
        energy_tolerance_j: 1e-6,
    }
}

/// Unit bulk length and cross-section give R=rho. A single latent cell has
/// M(h-h0)/dt = I^2*(1 + alpha*h) + independently authored/partitioned heat.
fn exact_endpoint(extra_power_w: f64) -> f64 {
    let a = DURATION * CURRENT * CURRENT / MASS;
    (OPENING_H + a + DURATION * extra_power_w / MASS) / (1.0 - a * RESISTIVITY_SLOPE)
}

fn exposed_bar() -> (
    VoxelDomain,
    [SolidMaterial; 1],
    ElectricSetup,
    ConductivityTensor,
    ThermalSetup,
    RadiationConfig,
) {
    let (_, solids, electric, conductivity) = bar();
    // Every top face of the three-cell solid bar borders transparent fluid.
    // Its upward hemisphere contains no other solid, so every ray escapes:
    // F=1 exactly, independently of the ray count or radiosity mode.
    let domain = VoxelDomain::from_fn(3, 2, 1, 1.0, |position| {
        if position[1] < 1.0 {
            Voxel::Solid(0)
        } else {
            Voxel::Fluid
        }
    })
    .unwrap();
    let mut thermal = ThermalSetup::new([ThermalFace::Adiabatic; 6]);
    thermal.fixed_temperature = vec![(3, 300.0), (4, 300.0), (5, 300.0)];
    let mut radiation = RadiationConfig::new(vec![0.04], [Some(250.0); 6]);
    radiation.rays_per_face = 32;
    radiation.patch_size = 1;
    radiation.max_iterations = 64;
    radiation.tolerance = 1e-10;
    (domain, solids, electric, conductivity, thermal, radiation)
}

fn exposed_bar_radiation_w(temperature_k: f64) -> f64 {
    0.04 * STEFAN_BOLTZMANN * (temperature_k.powi(4) - 250.0_f64.powi(4))
}

#[test]
fn current_and_phase_energy_converge_even_when_every_temperature_is_constant() {
    let (domain, solids, electric, k) = bar();
    let setup = ThermalSetup::new([ThermalFace::Adiabatic; 6]);
    let expected_h = exact_endpoint(0.0);
    for datum in [0.0, 1e6] {
        let law = phase(datum);
        let laws: [&dyn TensorCaloricLaw; 3] = [&law; 3];
        let previous = [datum + OPENING_H; 3];
        let mut source_states = Vec::new();
        let result = advance_caloric_tensor_electrothermal(
            &domain,
            &solids,
            &electric,
            &[],
            &setup,
            &[&k; 3],
            TensorCaloricStorage {
                laws: &laws,
                reference_mass_kg: &[MASS; 3],
                previous_specific_energy_j_kg: &previous,
            },
            &[],
            0.0,
            DURATION,
            &controls(),
            |cell, occupancy, energy, temperature| {
                assert_eq!(cell, 1);
                assert_eq!(occupancy, Voxel::Solid(0));
                assert_eq!(temperature, 300.0);
                source_states.push(energy - datum);
                Ok(Some(1.0 + RESISTIVITY_SLOPE * (energy - datum)))
            },
            &CancelGate::new(),
        )
        .unwrap();
        assert!(result.report.passes > 2);
        assert_eq!(result.report.current_solves, result.report.passes + 1);
        assert_eq!(source_states.len(), result.report.current_solves);
        near(source_states[0], OPENING_H, 1e-12);
        assert!(source_states.last().unwrap() > &source_states[0]);
        assert_eq!(
            result.thermal.endpoint.solution.solution.temperature,
            vec![300.0; 3]
        );
        let endpoint = &result.thermal.endpoint.specific_energy_j_kg;
        near(endpoint[1] - datum, expected_h, 3e-8);
        near(endpoint[0] - datum, OPENING_H, 1e-9);
        near(endpoint[2] - datum, OPENING_H, 1e-9);
        near(
            law.values_at_specific_enthalpy(endpoint[1])
                .unwrap()
                .liquid_mass_fraction,
            expected_h / 100.0,
            1e-9,
        );
        let stored = MASS * (endpoint[1] - previous[1]);
        near(result.thermal.record.stored_energy_change_j, stored, 1e-8);
        near(
            result.thermal.record.source_j,
            result.report.applied_joule_energy_j,
            1e-8,
        );
        near(result.thermal.record.boundary_outflow_j, 0.0, 1e-8);
        near(result.report.applied_electrode_energy_j, stored, 1e-6);
        near(result.report.thermostat_joule_energy_j, 0.0, 1e-12);
        near(
            result.endpoint_electric.total_joule_w,
            CURRENT * CURRENT * (1.0 + RESISTIVITY_SLOPE * expected_h),
            1e-7,
        );
        assert!(result.report.maximum_normalized_energy_change <= 1.0);
        assert!(result.report.source_mismatch_j <= controls().energy_tolerance_j);
        assert!(result.report.absolute_energy_defect_j <= controls().energy_tolerance_j);
        assert_eq!(previous, [datum + OPENING_H; 3]);
    }
}

#[test]
fn electrical_contact_heat_and_a_fixed_electrode_have_one_physical_destination() {
    let (domain, solids, electric, k) = bar();
    let law = phase(0.0);
    let laws: [&dyn TensorCaloricLaw; 3] = [&law; 3];
    let previous = [OPENING_H; 3];
    let contact = [ElectricalContact {
        cells: [0, 1],
        resistance_ohm_m2: 0.25,
        heat_fraction_to_first: 0.25,
    }];
    let mut setup = ThermalSetup::new([ThermalFace::Adiabatic; 6]);
    setup.fixed_temperature = vec![(0, 300.0)];
    setup.power_w = vec![0.0, 5.0, 0.0];
    let saved_setup = setup.clone();
    let contact_power = CURRENT * CURRENT * contact[0].resistance_ohm_m2;
    let thermostat_heat = DURATION * contact_power * contact[0].heat_fraction_to_first;
    let expected_h = exact_endpoint(5.0 + contact_power * 0.75);
    let result = advance_caloric_tensor_electrothermal(
        &domain,
        &solids,
        &electric,
        &contact,
        &setup,
        &[&k; 3],
        TensorCaloricStorage {
            laws: &laws,
            reference_mass_kg: &[MASS; 3],
            previous_specific_energy_j_kg: &previous,
        },
        &[(0, OPENING_H)],
        0.0,
        DURATION,
        &controls(),
        |_, _, energy, _| Ok(Some(1.0 + RESISTIVITY_SLOPE * energy)),
        &CancelGate::new(),
    )
    .unwrap();
    near(
        result.thermal.endpoint.specific_energy_j_kg[1],
        expected_h,
        3e-8,
    );
    assert_eq!(result.thermal.endpoint.specific_energy_j_kg[0], OPENING_H);
    assert_eq!(
        result.thermal.endpoint.solution.solution.temperature,
        vec![300.0; 3]
    );
    near(
        result.applied_electric.joule_w[0] * DURATION,
        thermostat_heat,
        1e-9,
    );
    near(
        result.report.thermostat_joule_energy_j,
        thermostat_heat,
        1e-9,
    );
    // Native prescribed cells omit applied source power; the coupling reports
    // their equal Joule source and immediate thermostat removal explicitly.
    near(
        result.report.applied_joule_energy_j + 5.0 * DURATION,
        result.thermal.record.source_j + result.report.thermostat_joule_energy_j,
        1e-8,
    );
    near(
        result.report.applied_electrode_energy_j + 5.0 * DURATION,
        result.thermal.record.stored_energy_change_j
            + result.thermal.record.boundary_outflow_j
            + result.report.thermostat_joule_energy_j,
        controls().energy_tolerance_j,
    );
    near(result.thermal.record.boundary_outflow_j, 0.0, 1e-8);
    assert!(result.report.absolute_energy_defect_j <= controls().energy_tolerance_j);
    assert_eq!(setup, saved_setup);
    assert_eq!(previous, [OPENING_H; 3]);
}

#[test]
fn exhausted_coupling_and_late_electrical_refusals_publish_no_heated_state() {
    let (domain, solids, electric, k) = bar();
    let law = phase(0.0);
    let laws: [&dyn TensorCaloricLaw; 3] = [&law; 3];
    let previous = [OPENING_H; 3];
    let setup = ThermalSetup::new([ThermalFace::Adiabatic; 6]);
    let saved_setup = setup.clone();
    let saved_electric = electric.clone();
    for mode in 0..3 {
        let gate = CancelGate::new();
        let mut config = controls();
        if mode == 0 {
            config.maximum_passes = 1;
        }
        let mut calls = 0;
        let result = advance_caloric_tensor_electrothermal(
            &domain,
            &solids,
            &electric,
            &[],
            &setup,
            &[&k; 3],
            TensorCaloricStorage {
                laws: &laws,
                reference_mass_kg: &[MASS; 3],
                previous_specific_energy_j_kg: &previous,
            },
            &[],
            0.0,
            DURATION,
            &config,
            |_, _, energy, _| {
                calls += 1;
                if calls == 2 {
                    if mode == 1 {
                        gate.request();
                    }
                    if mode == 2 {
                        return Err(ChtError::InvalidInput {
                            field: "test.electrical_source",
                            reason: "original phase response unavailable".into(),
                        });
                    }
                }
                Ok(Some(1.0 + RESISTIVITY_SLOPE * energy))
            },
            &gate,
        );
        assert_eq!(calls, 2);
        match mode {
            0 => assert!(matches!(result, Err(ChtError::SolverNotConverged { .. }))),
            1 => assert_eq!(result.unwrap_err(), ChtError::Cancelled),
            _ => assert!(matches!(
                result,
                Err(ChtError::InvalidInput {
                    field: "test.electrical_source",
                    ..
                })
            )),
        }
        assert_eq!(setup, saved_setup);
        assert_eq!(electric, saved_electric);
        assert_eq!(previous, [OPENING_H; 3]);
    }
}

#[test]
fn current_feedback_and_spatial_conduction_match_a_three_cell_energy_system() {
    let (domain, solids, electric, k) = bar();
    let law = EquilibriumEnthalpyPhaseCurve::try_new(
        ContentHash([0xE3; 32]),
        [(-200.0, 200.0), (400.0, 500.0)]
            .map(
                |(specific_enthalpy_j_kg, temperature_k)| EnthalpyPhaseKnot {
                    specific_enthalpy_j_kg,
                    temperature_k,
                    liquid_mass_fraction: 0.0,
                    bulk_density_kg_m3: MASS,
                },
            )
            .to_vec(),
    )
    .unwrap();
    let laws: [&dyn TensorCaloricLaw; 3] = [&law; 3];
    let previous = [0.0; 3];
    let setup = ThermalSetup::new([ThermalFace::Adiabatic; 6]);
    let result = advance_caloric_tensor_electrothermal(
        &domain,
        &solids,
        &electric,
        &[],
        &setup,
        &[&k; 3],
        TensorCaloricStorage {
            laws: &laws,
            reference_mass_kg: &[MASS; 3],
            previous_specific_energy_j_kg: &previous,
        },
        &[],
        0.0,
        DURATION,
        &controls(),
        |_, _, energy, _| Ok(Some(1.0 + RESISTIVITY_SLOPE * energy)),
        &CancelGate::new(),
    )
    .unwrap();

    // Cp=2 J/(kg K), G=1 W/K between unit cells, and M/dt=2 kg/s.
    // Each end row gives e_end=e_mid/5. The middle row then gives
    // (2 + .5*(2-2/5) - I^2*alpha)*e_mid = I^2.
    let capacity = 2.0;
    let conductance = 1.0;
    let storage_rate = MASS / DURATION;
    let endpoint_ratio = (conductance / capacity) / (storage_rate + conductance / capacity);
    let middle = CURRENT * CURRENT
        / (storage_rate + 2.0 * conductance / capacity * (1.0 - endpoint_ratio)
            - CURRENT * CURRENT * RESISTIVITY_SLOPE);
    let expected = [endpoint_ratio * middle, middle, endpoint_ratio * middle];
    for (cell, energy) in expected.iter().enumerate() {
        near(
            result.thermal.endpoint.specific_energy_j_kg[cell],
            *energy,
            3e-8,
        );
        near(
            result.thermal.endpoint.solution.solution.temperature[cell],
            300.0 + energy / capacity,
            2e-8,
        );
    }
    assert!(result.thermal.endpoint.specific_energy_j_kg[0] > 0.0);
    assert!(result.thermal.endpoint.specific_energy_j_kg[2] > 0.0);
    let stored = MASS * expected.iter().sum::<f64>();
    near(result.thermal.record.stored_energy_change_j, stored, 1e-7);
    near(result.report.applied_electrode_energy_j, stored, 1e-6);
    near(result.thermal.record.boundary_outflow_j, 0.0, 1e-8);
    near(
        result.endpoint_electric.total_joule_w,
        CURRENT * CURRENT * (1.0 + RESISTIVITY_SLOPE * middle),
        1e-7,
    );
    assert!(result.report.absolute_energy_defect_j <= controls().energy_tolerance_j);
    assert_eq!(previous, [0.0; 3]);
}

#[test]
fn radiating_current_feedback_closes_latent_contact_and_thermostat_energy() {
    let (domain, solids, electric, k, mut setup, mut radiation) = exposed_bar();
    setup.fixed_temperature.push((0, 300.0));
    setup.power_w = vec![0.0, 5.0, 0.0, 0.0, 0.0, 0.0];
    let saved = setup.clone();
    let contact = [ElectricalContact {
        cells: [0, 1],
        resistance_ohm_m2: 0.25,
        heat_fraction_to_first: 0.25,
    }];
    let law = phase(0.0);
    let laws: [&dyn TensorCaloricLaw; 6] = [&law; 6];
    let previous = [OPENING_H; 6];
    let prescribed = [
        (0, OPENING_H),
        (3, OPENING_H),
        (4, OPENING_H),
        (5, OPENING_H),
    ];
    let radiation_w = exposed_bar_radiation_w(300.0);
    // Contact heat is 4 W: 1 W goes to the fixed electrode and 3 W to
    // the free bulk cell. Radiation cools both free solid cells; the fixed
    // emitter's thermostat supplies its radiation without changing storage.
    let expected_middle = exact_endpoint(5.0 + 3.0 - radiation_w);
    let expected_end = OPENING_H - DURATION * radiation_w / MASS;
    for surface_exchange in [false, true] {
        radiation.surface_exchange = surface_exchange;
        let (result, report) = advance_caloric_tensor_electrothermal_radiating(
            &domain,
            &solids,
            &electric,
            &contact,
            &setup,
            &[&k; 6],
            TensorCaloricStorage {
                laws: &laws,
                reference_mass_kg: &[MASS; 6],
                previous_specific_energy_j_kg: &previous,
            },
            &prescribed,
            0.0,
            DURATION,
            &controls(),
            &radiation,
            |_, occupancy, energy, _| {
                Ok(match occupancy {
                    Voxel::Fluid => None,
                    Voxel::Solid(_) => Some(1.0 + RESISTIVITY_SLOPE * energy),
                })
            },
            &CancelGate::new(),
        )
        .unwrap();
        let endpoint = &result.thermal.endpoint;
        assert_eq!(endpoint.solution.solution.temperature, vec![300.0; 6]);
        near(endpoint.specific_energy_j_kg[1], expected_middle, 3e-8);
        near(endpoint.specific_energy_j_kg[2], expected_end, 3e-8);
        near(
            law.values_at_specific_enthalpy(endpoint.specific_energy_j_kg[1])
                .unwrap()
                .liquid_mass_fraction,
            expected_middle / 100.0,
            1e-9,
        );
        let stored = MASS * (expected_middle + expected_end - 2.0 * OPENING_H);
        near(result.thermal.record.stored_energy_change_j, stored, 1e-7);
        near(
            result.thermal.record.boundary_outflow_j,
            2.0 * radiation_w * DURATION,
            1e-7,
        );
        near(result.thermal.record.radiated_w, 3.0 * radiation_w, 1e-8);
        near(report.radiated_w, 3.0 * radiation_w, 1e-8);
        near(result.report.thermostat_joule_energy_j, DURATION, 1e-9);
        near(
            result.report.applied_joule_energy_j + 5.0 * DURATION,
            result.thermal.record.source_j + result.report.thermostat_joule_energy_j,
            1e-8,
        );
        near(
            result.report.applied_electrode_energy_j + 5.0 * DURATION,
            stored + 2.0 * radiation_w * DURATION + result.report.thermostat_joule_energy_j,
            controls().energy_tolerance_j,
        );
        assert_eq!(report.exposed_faces, 3);
        assert_eq!(report.rays, 3 * radiation.rays_per_face);
        assert_eq!(report.patches, if surface_exchange { 3 } else { 0 });
        assert!(report.relative_power_defect <= radiation.tolerance);
        assert!(result.report.passes > 2);
        assert_eq!(result.report.current_solves, result.report.passes + 1);
        assert!(result.report.absolute_energy_defect_j <= controls().energy_tolerance_j);
        assert_eq!(previous, [OPENING_H; 6]);
    }
    assert_eq!(setup, saved);
}

#[test]
fn radiating_electrothermal_endpoint_matches_independent_quartic_balance() {
    let (domain, solids, electric, k, mut setup, radiation) = exposed_bar();
    setup.fixed_temperature.extend([(0, 300.0), (2, 300.0)]);
    let law = EquilibriumEnthalpyPhaseCurve::try_new(
        ContentHash([0xE4; 32]),
        [(-200.0, 200.0), (400.0, 500.0)]
            .map(
                |(specific_enthalpy_j_kg, temperature_k)| EnthalpyPhaseKnot {
                    specific_enthalpy_j_kg,
                    temperature_k,
                    liquid_mass_fraction: 0.0,
                    bulk_density_kg_m3: MASS,
                },
            )
            .to_vec(),
    )
    .unwrap();
    let laws: [&dyn TensorCaloricLaw; 6] = [&law; 6];
    let previous = [0.0, 200.0, 0.0, 0.0, 0.0, 0.0];
    let dt = 0.25;
    let cp = 2.0;
    let alpha = 0.001;
    // Three unit thermal conductances meet the free cell: two electrode
    // thermostats and the fluid above. The independent scalar residual uses
    // the physical T^4 loss and energy-dependent current at the same endpoint.
    let outflow = |t: f64| 3.0 * (t - 300.0) + exposed_bar_radiation_w(t);
    let residual = |t: f64| {
        MASS * cp * (t - 400.0)
            + dt * (outflow(t) - CURRENT * CURRENT * (1.0 + alpha * cp * (t - 300.0)))
    };
    let (mut low, mut high) = (300.0, 400.0);
    assert!(residual(low) < 0.0 && residual(high) > 0.0);
    for _ in 0..80 {
        let middle = low + 0.5 * (high - low);
        if residual(middle) > 0.0 {
            high = middle;
        } else {
            low = middle;
        }
    }
    let expected = low + 0.5 * (high - low);
    let (result, report) = advance_caloric_tensor_electrothermal_radiating(
        &domain,
        &solids,
        &electric,
        &[],
        &setup,
        &[&k; 6],
        TensorCaloricStorage {
            laws: &laws,
            reference_mass_kg: &[MASS; 6],
            previous_specific_energy_j_kg: &previous,
        },
        &[(0, 0.0), (2, 0.0), (3, 0.0), (4, 0.0), (5, 0.0)],
        0.0,
        dt,
        &controls(),
        &radiation,
        |_, occupancy, energy, _| {
            Ok(match occupancy {
                Voxel::Fluid => None,
                Voxel::Solid(_) => Some(1.0 + alpha * energy),
            })
        },
        &CancelGate::new(),
    )
    .unwrap();
    near(
        result.thermal.endpoint.solution.solution.temperature[1],
        expected,
        1e-7,
    );
    near(
        result.thermal.endpoint.specific_energy_j_kg[1],
        cp * (expected - 300.0),
        2e-7,
    );
    near(
        result.thermal.record.stored_energy_change_j,
        MASS * cp * (expected - 400.0),
        5e-7,
    );
    near(
        result.thermal.record.boundary_outflow_j,
        dt * outflow(expected),
        5e-7,
    );
    near(
        result.report.applied_electrode_energy_j,
        MASS * cp * (expected - 400.0) + dt * outflow(expected),
        controls().energy_tolerance_j,
    );
    near(
        report.radiated_w,
        exposed_bar_radiation_w(expected) + 2.0 * exposed_bar_radiation_w(300.0),
        1e-7,
    );
    assert!(report.iterations > 1);
    assert!(result.report.absolute_energy_defect_j <= controls().energy_tolerance_j);
    assert_eq!(previous, [0.0, 200.0, 0.0, 0.0, 0.0, 0.0]);
}
