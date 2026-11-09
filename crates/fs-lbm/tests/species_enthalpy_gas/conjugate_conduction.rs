//! Mixed original gas enthalpy and solid calorics share one native face heat.
use super::*;
use fs_blake3::ContentHash;
use fs_lbm::conjugate::{
    ConductivityTensor, FaceContactResistance, SolidMaterial, TENSOR_SDIRK2_GAMMA,
    TensorCaloricLaw, TensorNonlinearConfig, ThermalFace, ThermalSetup,
    ThermochemicalConductionAdvance, ThermochemicalConductionControls,
    ThermochemicalConductionScheme, ThermochemicalSolidConductionCell,
};
use fs_material::phase::{
    EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, HeatCapacityKnot, SolidLiquidPhase,
};
use fs_material::thermal::transport::storage::{
    HeatStorageAmount, HeatStorageComponentQuery, HeatStorageThermalLaw, PreparedHeatStorage,
};

const HEAT_BUDGET_J: f64 = 1e-4;

fn controls(scheme: ThermochemicalConductionScheme) -> ThermochemicalConductionControls {
    let mut thermal = TensorNonlinearConfig::default();
    thermal.linear.tolerance = 1e-11;
    ThermochemicalConductionControls {
        thermal,
        scheme,
        conductivity: fs_lbm::conjugate::ThermochemicalConductivityModel::HeldAtOpening,
        maximum_caloric_queries: 100_000,
        maximum_enthalpy_defect_j: HEAT_BUDGET_J,
    }
}

fn gas_h(bank: &ThermochemicalGasModel, temperature: f64) -> f64 {
    bank.evaluate(&source::point(temperature, 100_000.0))
        .unwrap()
        .specific_enthalpy_j_kg()
}

fn isotropic(k: f64) -> ConductivityTensor {
    ConductivityTensor::new([[k, 0.0, 0.0], [0.0, k, 0.0], [0.0, 0.0, k]]).unwrap()
}

fn contact_setup(resistance: f64) -> ThermalSetup {
    let mut setup = ThermalSetup::new([ThermalFace::Adiabatic; 6]);
    setup.face_contacts.push(FaceContactResistance {
        cells: [0, 1],
        resistance_m2_k_w: resistance,
    });
    setup
}

fn melting(density: f64) -> EquilibriumEnthalpyPhaseCurve {
    // Cp=2 J/(kg K), with an actual 100 J/kg isothermal latent interval.
    EquilibriumEnthalpyPhaseCurve::try_new(
        ContentHash([0xC1; 32]),
        [
            (-100.0, 250.0, 0.0),
            (0.0, 300.0, 0.0),
            (100.0, 300.0, 1.0),
            (300.0, 400.0, 1.0),
        ]
        .map(|(h, t, liquid)| EnthalpyPhaseKnot {
            specific_enthalpy_j_kg: h,
            temperature_k: t,
            liquid_mass_fraction: liquid,
            bulk_density_kg_m3: density,
        })
        .to_vec(),
    )
    .unwrap()
}

fn constant_cp(density: f64) -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_from_heat_capacity_exact(
        ContentHash([0xC2; 32]),
        SolidLiquidPhase::Solid,
        -100.0,
        &[250.0, 450.0].map(|temperature_k| HeatCapacityKnot {
            temperature_k,
            specific_heat_capacity_j_kg_k: 2.0,
            bulk_density_kg_m3: density,
        }),
    )
    .unwrap()
}

#[test]
fn conjugate_graded_capacities_keep_exact_be_sdirk_modes_and_contact_heat() {
    let mut domain = VoxelDomain::graded([vec![0.5, 1.5], vec![1.0], vec![1.0]]).unwrap();
    domain.set(1, 0, 0, Voxel::Solid(0));
    let mass = [0.02, 0.0];
    let carrier = carrier(&domain, &mass);
    let before = carrier.clone();
    let species = [[mass[0]], [0.0]];
    let bank = source::constant_cp_model(-4000.0, 1.0);
    let names = [SpeciesId::new("manufactured").unwrap()];
    let context = source::point(400.0, 100_000.0);
    let d = [0.0; 2];
    let gate = CancelGate::new_clock_free();
    let gas = GasSpeciesEnthalpySource::new(&bank, &names, &context, &d, energy(), &gate).unwrap();
    let initial_t = [420.0, 280.0];
    let h = [gas_h(&bank, initial_t[0]), 0.0];
    let opening = gas.prepare(&domain, &carrier, &species, &h, &gate).unwrap();
    let gas_k = opening.transport(0).unwrap().thermal_conductivity_w_m_k();
    let solid_k = 3.0;
    let conductivity = isotropic(solid_k);
    let solid_mass = 2.0;
    let caloric = constant_cp(solid_mass / 1.5);
    let initial_solid_h = 2.0 * (initial_t[1] - 300.0);
    let solids = [SolidMaterial::new("graded original solid", solid_k)];
    let bindings = [ThermochemicalSolidConductionCell {
        cell: 1,
        reference_mass_kg: solid_mass,
        specific_energy_j_kg: initial_solid_h,
        caloric: &caloric,
        conductivity: &conductivity,
        prescribed_specific_energy_j_kg: None,
    }];
    let resistance = 0.75;
    let setup = contact_setup(resistance);
    let cp = 3.0 * R / MOLAR_MASS[0];
    let capacity = [mass[0] * cp, solid_mass * 2.0];
    // Geometry is graded: half widths .25 and .75, with unit interface area.
    let conductance = 1.0 / (0.25 / gas_k + 0.75 / solid_k + resistance);
    let z = 0.8;
    let duration = z / (conductance * (1.0 / capacity[0] + 1.0 / capacity[1]));
    let total_capacity = capacity[0] + capacity[1];
    let mean = (capacity[0] * initial_t[0] + capacity[1] * initial_t[1]) / total_capacity;
    for scheme in [
        ThermochemicalConductionScheme::BackwardEuler,
        ThermochemicalConductionScheme::Sdirk2,
    ] {
        let result = carrier
            .conduct_thermochemical_species_with_solids(
                &domain,
                &species,
                &h,
                &solids,
                &bindings,
                &setup,
                2.0,
                duration,
                &gas,
                controls(scheme),
                &gate,
            )
            .unwrap();
        // Independent two-capacity mode, not a call back into a caloric inverse.
        let factor = match scheme {
            ThermochemicalConductionScheme::BackwardEuler => 1.0 / (1.0 + z),
            ThermochemicalConductionScheme::Sdirk2 => {
                let gamma = TENSOR_SDIRK2_GAMMA;
                (1.0 + (2.0 * gamma - 1.0) * z) / (1.0 + gamma * z).powi(2)
            }
        };
        let difference = (initial_t[0] - initial_t[1]) * factor;
        let expected = [
            mean + capacity[1] * difference / total_capacity,
            mean - capacity[0] * difference / total_capacity,
        ];
        let endpoint = result.thermal.endpoint();
        let heat = &endpoint.solution.solution;
        for cell in 0..2 {
            near(heat.temperature[cell], expected[cell], 2e-7);
        }
        near(result.state.temperatures_k()[0], expected[0], 2e-7);
        near(
            result.gas_specific_enthalpy_j_kg[0] - h[0],
            cp * (expected[0] - initial_t[0]),
            3e-4,
        );
        near(
            endpoint.specific_energy_j_kg[1],
            2.0 * (expected[1] - 300.0),
            1e-6,
        );
        let q = conductance * (expected[0] - expected[1]);
        near(heat.outward_heat_w[0][1], q, 1e-7);
        near(heat.outward_heat_w[1][0], -q, 1e-7);
        near(
            heat.face_temperature_k[0][1],
            expected[0] - q * 0.25 / gas_k,
            2e-7,
        );
        near(
            heat.face_temperature_k[1][0],
            expected[1] + q * 0.75 / solid_k,
            2e-7,
        );
        near(
            heat.face_temperature_k[0][1] - heat.face_temperature_k[1][0],
            q * resistance,
            2e-7,
        );
        near(
            mass[0] * (endpoint.specific_energy_j_kg[0] - h[0])
                + solid_mass * (endpoint.specific_energy_j_kg[1] - initial_solid_h),
            0.0,
            HEAT_BUDGET_J,
        );
        near(
            result.thermal.record().stored_energy_change_j,
            0.0,
            HEAT_BUDGET_J,
        );
        near(
            result.thermal.record().boundary_outflow_j,
            0.0,
            HEAT_BUDGET_J,
        );
        near(result.thermal.record().closure_j, 0.0, HEAT_BUDGET_J);
        near(result.thermal.record().time_s, 2.0 + duration, 1e-10);
        // Complete mixed thermal state and reusable gas state have different masks.
        assert_eq!(result.gas_specific_enthalpy_j_kg[1], 0.0);
        assert_eq!(result.state.temperatures_k()[1], 0.0);
        assert!(result.state.transport(1).is_none());
        assert_eq!(result.state.diffusivity_m2_s(), d);
        assert!(heat.temperature[1] > initial_t[1] + 1.0);
    }
    assert_eq!(carrier, before);
    assert_eq!(species, [[mass[0]], [0.0]]);
    assert_eq!(h, [gas_h(&bank, initial_t[0]), 0.0]);
}

/// Integrated heat Q enters a 2 kg solid initially at h=20 J/kg and leaves
/// a constant-Cp gas initially at 360 K. The original piecewise law gives
/// Ts=300 on the plateau, then Ts=260+Q/4 after Q exceeds 160 J.
/// "carry" is the actual prior-stage heat, not a temperature increment.
fn latent_heat_oracle(a: f64, carry: f64, gas_capacity: f64) -> f64 {
    let plateau = (carry + a * 60.0) / (1.0 + a / gas_capacity);
    if plateau <= 160.0 {
        plateau
    } else {
        (carry + a * 100.0) / (1.0 + a * (1.0 / gas_capacity + 0.25))
    }
}

#[test]
fn conjugate_original_latent_plateau_and_crossing_keep_conservative_staged_storage() {
    let mut domain = VoxelDomain::new(2, 1, 1, 1.0).unwrap();
    domain.set(1, 0, 0, Voxel::Solid(0));
    let mass = [0.01, 0.0];
    let carrier = carrier(&domain, &mass);
    let species = [[mass[0]], [0.0]];
    let bank = source::constant_cp_model(-4000.0, 1.0);
    let names = [SpeciesId::new("manufactured").unwrap()];
    let context = source::point(400.0, 100_000.0);
    let d = [0.0; 2];
    let gate = CancelGate::new_clock_free();
    let gas = GasSpeciesEnthalpySource::new(&bank, &names, &context, &d, energy(), &gate).unwrap();
    let h = [gas_h(&bank, 360.0), 0.0];
    let opening = gas.prepare(&domain, &carrier, &species, &h, &gate).unwrap();
    let gas_k = opening.transport(0).unwrap().thermal_conductivity_w_m_k();
    let conductivity = isotropic(1.0);
    let caloric = melting(2.0);
    let solids = [SolidMaterial::new("original latent solid", 1.0)];
    let bindings = [ThermochemicalSolidConductionCell {
        cell: 1,
        reference_mass_kg: 2.0,
        specific_energy_j_kg: 20.0,
        caloric: &caloric,
        conductivity: &conductivity,
        prescribed_specific_energy_j_kg: None,
    }];
    let setup = contact_setup(2.0);
    let conductance = 1.0 / (0.5 / gas_k + 0.5 + 2.0);
    let gas_capacity = mass[0] * 3.0 * R / MOLAR_MASS[0];
    for a in [1.0, 5.0] {
        for scheme in [
            ThermochemicalConductionScheme::BackwardEuler,
            ThermochemicalConductionScheme::Sdirk2,
        ] {
            let duration = a / conductance;
            let gamma = TENSOR_SDIRK2_GAMMA;
            let stage1_q = latent_heat_oracle(gamma * a, 0.0, gas_capacity);
            let q = match scheme {
                ThermochemicalConductionScheme::BackwardEuler => {
                    latent_heat_oracle(a, 0.0, gas_capacity)
                }
                ThermochemicalConductionScheme::Sdirk2 => {
                    latent_heat_oracle(gamma * a, (1.0 - gamma) * stage1_q / gamma, gas_capacity)
                }
            };
            let result = carrier
                .conduct_thermochemical_species_with_solids(
                    &domain,
                    &species,
                    &h,
                    &solids,
                    &bindings,
                    &setup,
                    0.0,
                    duration,
                    &gas,
                    controls(scheme),
                    &gate,
                )
                .unwrap();
            let endpoint = result.thermal.endpoint();
            let heat = &endpoint.solution.solution;
            let expected_solid_h = 20.0 + q / 2.0;
            let expected_gas_t = 360.0 - q / gas_capacity;
            let expected_solid_t = if expected_solid_h <= 100.0 {
                300.0
            } else {
                300.0 + 0.5 * (expected_solid_h - 100.0)
            };
            near(heat.temperature[0], expected_gas_t, 2e-7);
            near(heat.temperature[1], expected_solid_t, 2e-7);
            near(endpoint.specific_energy_j_kg[1], expected_solid_h, 1e-6);
            near(
                mass[0] * (endpoint.specific_energy_j_kg[0] - h[0]),
                -q,
                HEAT_BUDGET_J,
            );
            near(
                2.0 * (endpoint.specific_energy_j_kg[1] - 20.0),
                q,
                HEAT_BUDGET_J,
            );
            let phase = caloric
                .values_at_specific_enthalpy(endpoint.specific_energy_j_kg[1])
                .unwrap();
            if a == 1.0 {
                assert_eq!(heat.temperature[1], 300.0);
                assert_eq!(
                    caloric
                        .temperature_slope_kg_k_j(endpoint.specific_energy_j_kg[1])
                        .unwrap(),
                    0.0
                );
                near(phase.liquid_mass_fraction, expected_solid_h / 100.0, 1e-8);
                assert!(phase.liquid_mass_fraction > 0.2 && phase.liquid_mass_fraction < 1.0);
            } else {
                assert!(expected_solid_h > 100.0);
                assert_eq!(phase.liquid_mass_fraction, 1.0);
                assert!(heat.temperature[1] > 300.0);
            }
            if let ThermochemicalConductionAdvance::Sdirk2(staged) = &result.thermal {
                near(
                    staged.stages[0].specific_energy_j_kg[1],
                    20.0 + stage1_q / 2.0,
                    1e-6,
                );
                near(
                    staged.stages[0].solution.solution.temperature[0],
                    360.0 - stage1_q / gas_capacity,
                    2e-7,
                );
            }
            let endpoint_q = conductance * (expected_gas_t - expected_solid_t);
            near(heat.outward_heat_w[0][1], endpoint_q, 1e-7);
            near(heat.outward_heat_w[1][0], -endpoint_q, 1e-7);
            near(
                result.thermal.record().stored_energy_change_j,
                0.0,
                HEAT_BUDGET_J,
            );
            near(result.thermal.record().closure_j, 0.0, HEAT_BUDGET_J);
            assert_eq!(result.gas_specific_enthalpy_j_kg[1], 0.0);
            assert_eq!(result.state.temperatures_k()[1], 0.0);
        }
    }
    assert_eq!(carrier.mass_kg, mass);
    assert_eq!(bindings[0].specific_energy_j_kg, 20.0);
}

#[test]
fn conjugate_prescribed_plateau_keeps_explicit_latent_energy_and_thermostat_heat() {
    let mut domain = VoxelDomain::new(2, 1, 1, 1.0).unwrap();
    domain.set(1, 0, 0, Voxel::Solid(0));
    let mass = [0.01, 0.0];
    let carrier = carrier(&domain, &mass);
    let species = [[mass[0]], [0.0]];
    let bank = source::constant_cp_model(-4000.0, 1.0);
    let names = [SpeciesId::new("manufactured").unwrap()];
    let context = source::point(400.0, 100_000.0);
    let d = [0.0; 2];
    let gate = CancelGate::new_clock_free();
    let gas = GasSpeciesEnthalpySource::new(&bank, &names, &context, &d, energy(), &gate).unwrap();
    let h = [gas_h(&bank, 360.0), 0.0];
    let opening = gas.prepare(&domain, &carrier, &species, &h, &gate).unwrap();
    let gas_k = opening.transport(0).unwrap().thermal_conductivity_w_m_k();
    let caloric = melting(2.0);
    let conductivity = isotropic(1.0);
    let solids = [SolidMaterial::new("explicit latent thermostat", 1.0)];
    let base = ThermochemicalSolidConductionCell {
        cell: 1,
        reference_mass_kg: 2.0,
        specific_energy_j_kg: 20.0,
        caloric: &caloric,
        conductivity: &conductivity,
        prescribed_specific_energy_j_kg: None,
    };
    let mut setup = contact_setup(2.0);
    setup.fixed_temperature.push((1, 300.0));
    let saved = setup.clone();
    let conductance = 1.0 / (0.5 / gas_k + 0.5 + 2.0);
    let duration = 1.0 / conductance;
    let gas_capacity = mass[0] * 3.0 * R / MOLAR_MASS[0];
    let expected_gas_t = (gas_capacity * 360.0 + 300.0) / (gas_capacity + 1.0);
    let gas_heat = gas_capacity * (360.0 - expected_gas_t);
    let solve = |binding: ThermochemicalSolidConductionCell<'_>, setup: &ThermalSetup| {
        carrier.conduct_thermochemical_species_with_solids(
            &domain,
            &species,
            &h,
            &solids,
            &[binding],
            setup,
            0.0,
            duration,
            &gas,
            controls(ThermochemicalConductionScheme::BackwardEuler),
            &gate,
        )
    };
    for prescribed in [40.0, 80.0] {
        let result = solve(
            ThermochemicalSolidConductionCell {
                prescribed_specific_energy_j_kg: Some(prescribed),
                ..base
            },
            &setup,
        )
        .unwrap();
        let endpoint = result.thermal.endpoint();
        near(endpoint.specific_energy_j_kg[1], prescribed, 1e-12);
        assert_eq!(endpoint.solution.solution.temperature[1], 300.0);
        near(result.state.temperatures_k()[0], expected_gas_t, 2e-7);
        near(
            endpoint.solution.solution.outward_heat_w[0][1],
            gas_heat / duration,
            1e-7,
        );
        let solid_heat = 2.0 * (prescribed - 20.0);
        near(
            result.thermal.record().stored_energy_change_j,
            solid_heat - gas_heat,
            HEAT_BUDGET_J,
        );
        near(
            result.thermal.record().boundary_outflow_j,
            gas_heat - solid_heat,
            HEAT_BUDGET_J,
        );
        near(result.thermal.record().source_j, 0.0, HEAT_BUDGET_J);
        near(result.thermal.record().closure_j, 0.0, HEAT_BUDGET_J);
    }
    // A temperature on the plateau cannot select either physical energy.
    assert!(solve(base, &setup).is_err());
    assert!(
        solve(
            ThermochemicalSolidConductionCell {
                prescribed_specific_energy_j_kg: Some(160.0),
                ..base
            },
            &setup,
        )
        .is_err()
    );
    assert!(
        solve(
            ThermochemicalSolidConductionCell {
                prescribed_specific_energy_j_kg: Some(40.0),
                ..base
            },
            &contact_setup(2.0),
        )
        .is_err()
    );
    assert_eq!(setup, saved);
    assert_eq!(base.specific_energy_j_kg, 20.0);
}

/// A real amount-bound original storage source, not a substitute ThermalLaw.
fn capacity_assembly() -> PreparedHeatStorage {
    use fs_evidence::ValidityDomain;
    use fs_matdb::{
        ClaimSelection, ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId,
        PropertyClaim, PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy,
        UncertaintyModel,
    };
    use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
    use fs_qty::{Dims, QuantitySpec, Temperature};
    let temperature = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ));
    let dims = Dims([2, 0, -2, -1, 0, 0]);
    let capacity = PropertyKey::new("conjugate-reference-capacity", dims);
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: capacity.clone(),
            value: PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![(300.0, 2.0), (320.0, 2.0)],
                dims,
            },
            validity: ValidityDomain::unconstrained().with_quantity("T", temperature, 300.0, 320.0),
            interpolation: InterpolationPolicy::LinearInside,
            uncertainty: UncertaintyModel::Unstated,
            observations: vec![],
            provenance: Provenance {
                source: "authored fixed-inventory conjugate test".into(),
                license: "test-only".into(),
                artifact: None,
            },
        })
        .unwrap();
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic solid".into(),
            phase: "solid".into(),
            process: "fixed amount".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let lower = QueryPoint::new()
        .with_quantity("T", temperature, 300.0)
        .unwrap();
    let upper = QueryPoint::new()
        .with_quantity("T", temperature, 320.0)
        .unwrap();
    PreparedHeatStorage::try_new(&[HeatStorageComponentQuery {
        id: 1,
        card: &card,
        capacity: &capacity,
        amount: HeatStorageAmount::MassKg(2.0),
        temperature_axis: "T",
        lower: &lower,
        upper: &upper,
        selection: ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
    }])
    .unwrap()
}

#[test]
fn conjugate_bindings_preserve_source_mass_constraints_and_refuse_atomically() {
    let mut domain = VoxelDomain::new(2, 1, 1, 1.0).unwrap();
    domain.set(1, 0, 0, Voxel::Solid(0));
    let carrier = carrier(&domain, &[0.01, 0.0]);
    let before = carrier.clone();
    let species = [[0.01], [0.0]];
    let bank = source::constant_cp_model(0.0, 1.0);
    let names = [SpeciesId::new("manufactured").unwrap()];
    let context = source::point(400.0, 100_000.0);
    let d = [0.0; 2];
    let gate = CancelGate::new_clock_free();
    let gas = GasSpeciesEnthalpySource::new(&bank, &names, &context, &d, energy(), &gate).unwrap();
    let h = [gas_h(&bank, 310.0), 0.0];
    let storage = capacity_assembly();
    let caloric = HeatStorageThermalLaw::try_new(&storage, 2.0, 300.0).unwrap();
    let conductivity = isotropic(1.0);
    let solids = [SolidMaterial::new("original fixed-mass assembly", 1.0)];
    let binding = ThermochemicalSolidConductionCell {
        cell: 1,
        reference_mass_kg: 2.0,
        specific_energy_j_kg: 20.0, // Original source T=310 K.
        caloric: &caloric,
        conductivity: &conductivity,
        prescribed_specific_energy_j_kg: None,
    };
    let setup = contact_setup(2.0);
    let solve = |bindings: &[ThermochemicalSolidConductionCell<'_>], gate: &CancelGate| {
        carrier.conduct_thermochemical_species_with_solids(
            &domain,
            &species,
            &h,
            &solids,
            bindings,
            &setup,
            0.0,
            1.0,
            &gas,
            controls(ThermochemicalConductionScheme::BackwardEuler),
            gate,
        )
    };
    let unchanged = solve(&[binding], &gate).unwrap();
    near(
        unchanged.thermal.endpoint().solution.solution.temperature[1],
        310.0,
        1e-9,
    );
    near(
        unchanged.thermal.record().stored_energy_change_j,
        0.0,
        HEAT_BUDGET_J,
    );
    for invalid in [
        vec![],
        vec![binding, binding],
        vec![ThermochemicalSolidConductionCell { cell: 0, ..binding }],
        vec![ThermochemicalSolidConductionCell { cell: 2, ..binding }],
        vec![ThermochemicalSolidConductionCell {
            reference_mass_kg: 0.0,
            ..binding
        }],
        vec![ThermochemicalSolidConductionCell {
            reference_mass_kg: 1.0,
            ..binding
        }],
        vec![ThermochemicalSolidConductionCell {
            reference_mass_kg: 4.0,
            ..binding
        }],
        vec![ThermochemicalSolidConductionCell {
            specific_energy_j_kg: 41.0,
            ..binding
        }],
    ] {
        assert!(solve(&invalid, &gate).is_err());
    }
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(matches!(
        solve(&[binding], &cancelled),
        Err(ChtError::Cancelled)
    ));
    assert_eq!(carrier, before);
    assert_eq!(species, [[0.01], [0.0]]);
    assert_eq!(h, [gas_h(&bank, 310.0), 0.0]);
    assert_eq!(binding.reference_mass_kg, 2.0);
    assert_eq!(binding.specific_energy_j_kg, 20.0);
}

#[test]
fn radiating_gas_solid_interval_cools_latent_inventory_and_reuses_original_gas_sources() {
    use fs_lbm::conjugate::{RadiationConfig, ThermochemicalConductivityModel};
    use fs_material::thermochemical_gas::conductivity::GasConductivityQuadrature;

    let mut domain = VoxelDomain::new(2, 1, 1, 1.0).unwrap();
    domain.set(1, 0, 0, Voxel::Solid(0));
    let mass = [0.02, 0.0];
    let carrier = carrier(&domain, &mass);
    let before = carrier.clone();
    let species = [[mass[0]], [0.0]];
    let bank = source::constant_cp_model(-4000.0, 1.0);
    let names = [SpeciesId::new("manufactured").unwrap()];
    let context = source::point(350.0, 100_000.0);
    let d = [0.0; 2];
    let gate = CancelGate::new_clock_free();
    let source =
        GasSpeciesEnthalpySource::new(&bank, &names, &context, &d, energy(), &gate).unwrap();
    let h = [gas_h(&bank, 350.0), 0.0];
    let opening = source
        .prepare(&domain, &carrier, &species, &h, &gate)
        .unwrap();
    let gas_k = opening.transport(0).unwrap().thermal_conductivity_w_m_k();
    let solid_k = 3.0;
    let conductivity = isotropic(solid_k);
    let caloric = melting(2.0);
    let bindings = [ThermochemicalSolidConductionCell {
        cell: 1,
        reference_mass_kg: 2.0,
        specific_energy_j_kg: 60.0,
        caloric: &caloric,
        conductivity: &conductivity,
        prescribed_specific_energy_j_kg: None,
    }];
    let solids = [SolidMaterial::new("radiating latent solid", solid_k)];
    let setup = contact_setup(0.75);
    let mut radiation = RadiationConfig::new(vec![0.8], [Some(250.0); 6]);
    radiation.surface_exchange = false;
    radiation.rays_per_face = 64;
    radiation.tolerance = 1e-10;
    let duration = 0.1;
    let cp = 3.0 * R / MOLAR_MASS[0];
    let capacity = mass[0] * cp;
    let conductance = 1.0 / (0.5 / gas_k + 0.5 / solid_k + 0.75);
    // One unit-area solid face looks through transparent gas into a uniform
    // enclosure. Every ray sees the same temperature, so no fitted view factor
    // is used by this independent physical-energy calculation.
    let radiated_w = 0.8 * 5.670_374_419e-8 * (300.0_f64.powi(4) - 250.0_f64.powi(4));
    for scheme in [
        ThermochemicalConductionScheme::BackwardEuler,
        ThermochemicalConductionScheme::Sdirk2,
    ] {
        let mut control = controls(scheme);
        control.conductivity = ThermochemicalConductivityModel::SourceIntegral {
            quadrature: GasConductivityQuadrature {
                absolute_tolerance_w_m_k: 1e-12,
                relative_tolerance: 1e-10,
                maximum_subdivisions: 256,
                maximum_depth: 16,
            },
            maximum_evaluations: 2_000_000,
        };
        let result = carrier
            .conduct_thermochemical_species_with_solids_radiating(
                &domain, &species, &h, &solids, &bindings, &setup, 0.0, duration, &source, control,
                &radiation, &gate,
            )
            .unwrap();
        let z = duration * conductance / capacity;
        let decay = match scheme {
            ThermochemicalConductionScheme::BackwardEuler => 1.0 / (1.0 + z),
            ThermochemicalConductionScheme::Sdirk2 => {
                let gamma = TENSOR_SDIRK2_GAMMA;
                (1.0 + (2.0 * gamma - 1.0) * z) / (1.0 + gamma * z).powi(2)
            }
        };
        let expected_gas_t = 300.0 + 50.0 * decay;
        let gas_change = capacity * (expected_gas_t - 350.0);
        let expected_solid_h = 60.0 + (-gas_change - duration * radiated_w) / 2.0;
        assert!(expected_solid_h > 0.0 && expected_solid_h < 100.0);
        let endpoint = result.thermal.endpoint();
        near(
            endpoint.solution.solution.temperature[0],
            expected_gas_t,
            2e-7,
        );
        near(endpoint.solution.solution.temperature[1], 300.0, 1e-10);
        near(endpoint.specific_energy_j_kg[1], expected_solid_h, 2e-6);
        near(
            result.gas_specific_enthalpy_j_kg[0] - h[0],
            gas_change / mass[0],
            2e-4,
        );
        assert_eq!(result.gas_specific_enthalpy_j_kg[1], 0.0);
        near(result.thermal.record().radiated_w, radiated_w, 1e-8);
        near(
            result.thermal.record().boundary_outflow_j,
            duration * radiated_w,
            2e-7,
        );
        assert!(result.thermal.record().closure_j.abs() < HEAT_BUDGET_J);
        assert_eq!(
            result.thermal.radiation_reports().len(),
            match scheme {
                ThermochemicalConductionScheme::BackwardEuler => 1,
                ThermochemicalConductionScheme::Sdirk2 => 2,
            }
        );
        assert!(result.report.conductivity_evaluations > 0);

        let exhausted = ThermochemicalConductionControls {
            maximum_caloric_queries: 1,
            ..control
        };
        assert!(
            carrier
                .conduct_thermochemical_species_with_solids_radiating(
                    &domain, &species, &h, &solids, &bindings, &setup, 0.0, duration, &source,
                    exhausted, &radiation, &gate,
                )
                .is_err()
        );
        let cancelled = CancelGate::new_clock_free();
        cancelled.request();
        assert!(matches!(
            carrier.conduct_thermochemical_species_with_solids_radiating(
                &domain, &species, &h, &solids, &bindings, &setup, 0.0, duration, &source, control,
                &radiation, &cancelled,
            ),
            Err(ChtError::Cancelled)
        ));
        assert_eq!(carrier, before);
        assert_eq!(bindings[0].specific_energy_j_kg, 60.0);
        assert_eq!(species, [[0.02], [0.0]]);
    }
}
