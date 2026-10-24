//! G1/G3/G4: finite source gas inside the ORIGINAL hereditary heat feedback.
use super::*;
use fs_couple::cavity_phs::piston::thermal::HeatExchangeSolve;
use fs_couple::finite_gas_thermal::FiniteGasThermalEnvironment;
use fs_material::state_point::MaterialPropertySelection;
use fs_material::thermochemical_gas::{
    SpeciesTransportMaterial, ThermochemicalGasModel,
    energy::{GasEnergySolve, IsochoricGasParcel},
};

fn gas(temperature_k: f64, volume_m3: f64) -> IsochoricGasParcel {
    use fs_matdb::{ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter};
    use fs_qty::{DynViscosity, MolarMass};
    use fs_thermochem::mixture::FrozenIdealGasMixtureModelV1;
    use fs_thermochem::{
        Composition, CompositionBasis, ElementalReferenceIdV1, NASA9_LAW_ID_V1,
        NASA9_LAW_VERSION_V1, NASA9_STATE_SCHEMA_VERSION_V1, Nasa9RegionV1,
        Nasa9StandardStateModelV1, ReferenceEquationOfStateV1, SpeciesId,
        StandardStateConventionV1, StandardStatePhaseV1, required_nasa9_parameter_dims,
    };
    let typed = |kind| QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static));
    let provenance = Provenance {
        source: "manufactured finite thermal gas, not measured data".into(),
        license: "CC0-1.0".into(),
        artifact: None,
    };
    let mut parameters = std::collections::BTreeMap::new();
    for ((key, dims), value) in required_nasa9_parameter_dims()
        .into_iter().take(9).zip([0.0, 0.0, 3.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    {
        parameters.insert(key.to_owned(), LawParameter { value, dims });
    }
    parameters.insert("reference_pressure".into(), LawParameter {
        value: 100000.0, dims: Pressure::DIMS,
    });
    let nasa = ConstitutiveModelCard {
        law: LawId(NASA9_LAW_ID_V1.into()), law_version: NASA9_LAW_VERSION_V1,
        parameters, state_schema_version: NASA9_STATE_SCHEMA_VERSION_V1,
        initial_state: InitialStatePolicy::ZeroInternalState,
        validity: ValidityDomain::unconstrained().with("T", 200.0, 900.0),
        sources: vec![], provenance: provenance.clone(),
    };
    let species = SpeciesId::new("manufactured-thermal-gas").unwrap();
    let thermo = Nasa9StandardStateModelV1::new(
        species.clone(), MolarMass::new(0.028),
        StandardStateConventionV1::new(
            StandardStatePhaseV1::Gas, ReferenceEquationOfStateV1::IdealGas,
            Pressure::new(100000.0), ElementalReferenceIdV1::new("manufactured-reference").unwrap(),
        ).unwrap(),
        vec![Nasa9RegionV1::from_card(nasa).unwrap()],
    ).unwrap();
    let mut claims = ClaimSet::new();
    for (name, dims, value) in [
        ("dynamic_viscosity", DynViscosity::DIMS, 1.7e-5),
        ("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS, 0.025),
    ] {
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(name, dims),
            value: PropertyValue::Scalar { value, dims },
            validity: ValidityDomain::unconstrained()
                .with_quantity("temperature", typed(QuantityKind::AbsoluteTemperature), 200.0, 900.0)
                .with_quantity("pressure", typed(QuantityKind::Pressure), 1000.0, 1e7),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::LinearInside,
            observations: vec![], provenance: provenance.clone(),
        }).unwrap();
    }
    let transport = SpeciesTransportMaterial {
        species,
        material: MaterialCard::assemble(MaterialStateId {
            chemistry: "manufactured gas".into(), phase: "gas".into(),
            process: "numerical fixture".into(), revision: 0,
        }, claims, vec![]).unwrap(),
        selection: MaterialPropertySelection::SingleClaimOnly,
    };
    let model = ThermochemicalGasModel::new(
        FrozenIdealGasMixtureModelV1::new(
            vec![thermo], Composition::new(CompositionBasis::MoleFraction, vec![1.0]).unwrap(),
        ).unwrap(),
        vec![transport], None,
    ).unwrap();
    let point = QueryPoint::new()
        .with_quantity("temperature", typed(QuantityKind::AbsoluteTemperature), temperature_k).unwrap()
        .with_quantity("pressure", typed(QuantityKind::Pressure), 100000.0).unwrap()
        .with("material-lot", 7.0).unwrap();
    IsochoricGasParcel::new(model, point, volume_m3).unwrap()
}

fn gas_solve() -> GasEnergySolve {
    GasEnergySolve {
        minimum_temperature_k: 200.0, maximum_temperature_k: 800.0,
        energy_tolerance_j_kg: 1e-7, maximum_iterations: 128,
    }
}
fn exchange(bracket: f64) -> HeatExchangeSolve {
    HeatExchangeSolve {
        heat_bracket_j: [-bracket, bracket], source_tolerance_j: 1e-12,
        // An inner thermal report cannot exceed its enclosing material gate.
        energy_tolerance_j: config().feedback.energy_tolerance_j, maximum_iterations: 64,
    }
}
fn boundary(source: &Sources, film: f64) -> LumpedThermalEnvironment {
    LumpedThermalEnvironment {
        temperature_k: 300.0, radiation_temperature_k: 300.0,
        convection_w_per_m2_k: film, transport: source.transport.clone(),
        enthalpy_tolerance_j_kg: 1e-12, maximum_thermal_residual_j: 1e-11,
    }
}
fn scope(f: impl FnOnce(&Cx<'_>, &CancelGate)) {
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(&gate, arena,
            StreamKey { seed: 0, kernel_id: 0x464754, tile: 0, iteration: 0 },
            Budget::INFINITE, ExecMode::Deterministic);
        f(&cx, &gate);
    });
}

#[test]
fn g1_finite_source_gas_matches_two_capacity_solution_and_opposite_heat() {
    scope(|cx, _| {
        let source = sources(2.0, 1000.0, false, false);
        for (body_k, gas_k) in [(300.0, 320.0), (320.0, 300.0)] {
            let old = Oscillator::try_new(&source.mechanical, &source.caloric,
                geometry(0.0), body_k, 0.0, 0.0).unwrap();
            let reservoir = gas(gas_k, 1e-6);
            let saved = reservoir.clone();
            let mut body = boundary(&source, 100.0);
            body.temperature_k = 12345.0; // Replaced by actual trial gas, never used as a bath.
            let environment = FiniteGasThermalEnvironment {
                gas: &reservoir, body: &body, gas_heat_j: 2e-4,
                gas_solve: gas_solve(), exchange: exchange(0.02),
                maximum_gas_temperature_increment_k: 100.0,
            };
            let input = UniformEnthalpyStepInput {
                curve: source.caloric.curve(), initial: old.point().caloric_state(),
                mass_kg: old.point().mass_kg(), volume_m3: old.point().volume_m3(),
                surface_area_m2: old.geometry().thermal_surface_m2,
                internal_heat_j: 4e-4, duration_s: 0.01,
            };
            let result = environment.advance(cx, input).unwrap();
            let cs = input.mass_kg * 2.0;
            let cv = 2.5 * fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / 0.028;
            let cg = reservoir.mass_kg() * cv;
            let a = input.duration_s * body.convection_w_per_m2_k * input.surface_area_m2;
            let source_body_k = body_k + input.internal_heat_j / cs;
            let source_gas_k = gas_k + environment.gas_heat_j / cg;
            let expected_q = a * (source_body_k - source_gas_k)
                / (1.0 + a * (1.0 / cs + 1.0 / cg));
            let gas_tolerance = gas_solve().energy_tolerance_j_kg / cv;
            let q_tolerance = exchange(0.02).source_tolerance_j
                + a * (body.maximum_thermal_residual_j / cs + gas_tolerance)
                + 64.0 * f64::EPSILON;
            close(result.report.solid_to_gas_heat_j, expected_q, 2.0 * q_tolerance);
            close(result.state.temperature_k(), source_body_k - expected_q / cs,
                2.0 * q_tolerance / cs + body.maximum_thermal_residual_j / cs + 1e-10);
            close(result.report.next_gas.state().state().temperature,
                source_gas_k + expected_q / cg,
                2.0 * q_tolerance / cg + gas_tolerance + 1e-10);
            assert_eq!(result.external_heat_j, -result.report.solid_to_gas_heat_j);
            assert_eq!(result.report.gas.supplied_heat_j,
                environment.gas_heat_j + result.report.solid_to_gas_heat_j);
            let storage = input.mass_kg * (result.state.specific_enthalpy_j_kg()
                - input.initial.specific_enthalpy_j_kg())
                + reservoir.mass_kg() * (result.report.next_gas.state().specific_internal_energy_j_kg()
                    - reservoir.state().specific_internal_energy_j_kg());
            close(storage, input.internal_heat_j + environment.gas_heat_j,
                environment.exchange.energy_tolerance_j);
            assert!(result.report.source_residual_j.abs() <= environment.exchange.source_tolerance_j);
            assert!(result.report.solid_energy_residual_j.abs() <= environment.exchange.energy_tolerance_j);
            assert_eq!(result.report.next_gas.point().axes().get("material-lot"), Some(&7.0));
            assert_eq!(reservoir, saved, "proposal cannot debit the accepted reservoir");
        }
    });
}

#[test]
fn g3_finite_inventory_changes_actual_hereditary_vibration_without_heat_leaks() {
    scope(|cx, gate| {
        let source = sources_with_modulus_slope(0.02, 1000.0, true, true, None, -0.005);
        let initial = Oscillator::try_new(&source.mechanical, &source.caloric,
            geometry(0.0), 300.0, 0.0, 0.5).unwrap();
        let mut states = [initial.clone(), initial.clone()];
        let mut gases = [gas(330.0, 1e-7), gas(330.0, 1e-4)];
        let original_gases = gases.clone();
        let body = boundary(&source, 100.0);
        let mut numerical = [0.0; 2];
        let mut reversible = [0.0; 2];
        let mut max_velocity_difference = 0.0_f64;
        for _ in 0..32 {
            for index in 0..2 {
                let old = &states[index];
                let environment = FiniteGasThermalEnvironment {
                    gas: &gases[index], body: &body, gas_heat_j: 0.0,
                    gas_solve: gas_solve(), exchange: exchange(8e-4),
                    maximum_gas_temperature_increment_k: 100.0,
                };
                let step = old.step(0.0, 0.00005, config(),
                    |input| environment.advance(cx, input), gate).unwrap();
                let expected = old.point().mechanical().advance(
                    step.next.point().caloric_state().temperature_k(),
                    step.next.point().mechanical().total_strain(),
                    0.00005, config().feedback.energy_tolerance_j / old.point().volume_m3(),
                ).unwrap();
                assert_eq!(step.next.point().mechanical(), &expected.point);
                assert!(step.material_heat_j >= 0.0);
                close(step.applied_material_heat_j,
                    step.material_heat_j + step.reversible_material_heat_j,
                    config().feedback.energy_tolerance_j);
                close(step.reversible_material_heat_j
                    + old.point().volume_m3() * step.material_energy.modulus_work
                    + step.mechanical_entropy_energy_change_j, 0.0, 1e-12);
                close(step.next.point().mechanical_internal_energy_j()
                    - old.point().mechanical_internal_energy_j(),
                    step.material_internal_energy_change_j, 1e-12);
                assert_eq!(step.transport.external_heat_j,
                    -step.transport.report.solid_to_gas_heat_j);
                assert!(step.transport.report.combined_energy_residual_j.abs() <= exchange(8e-4).energy_tolerance_j);
                numerical[index] += step.numerical_loss_j;
                reversible[index] += step.reversible_material_heat_j.abs();
                gases[index] = step.transport.report.next_gas;
                states[index] = step.next;
            }
            max_velocity_difference = max_velocity_difference
                .max((states[0].velocity_m_s() - states[1].velocity_m_s()).abs());
        }
        assert!(states[1].point().caloric_state().temperature_k()
            > states[0].point().caloric_state().temperature_k());
        assert!(max_velocity_difference > 1e-7);
        for index in 0..2 {
            let state = &states[index];
            let body_storage = state.point().mass_kg()
                * (state.point().caloric_state().specific_enthalpy_j_kg()
                    - initial.point().caloric_state().specific_enthalpy_j_kg());
            let gas_storage = gases[index].mass_kg()
                * (gases[index].state().specific_internal_energy_j_kg()
                    - original_gases[index].state().specific_internal_energy_j_kg());
            let stored = state.kinetic_energy_j()
                + state.point().mechanical_internal_energy_j()
                + body_storage + gas_storage;
            close(stored + numerical[index], initial.kinetic_energy_j(), 5e-8);
            assert!(numerical[index] > 0.0, "BE loss is an explicit nonthermal account");
            assert!(reversible[index] > 0.0, "stiffness and strain changes exchange reversible heat");
        }
    });
}

#[test]
fn g3_insulation_preserves_original_motion_and_radiation_stays_external() {
    scope(|cx, gate| {
        let source = sources(2.0, 1000.0, true, false);
        let state = Oscillator::try_new(&source.mechanical, &source.caloric,
            geometry(0.0), 300.0, 0.0, 0.5).unwrap();
        let reservoir = gas(300.0, 1e-6);
        let body = boundary(&source, 0.0);
        let environment = FiniteGasThermalEnvironment {
            gas: &reservoir, body: &body, gas_heat_j: 0.0,
            gas_solve: gas_solve(), exchange: exchange(0.02),
            maximum_gas_temperature_increment_k: 100.0,
        };
        let finite = state.step(0.0, 0.00005, config(),
            |input| environment.advance(cx, input), gate).unwrap();
        let original = state.step(0.0, 0.00005, config(),
            |input| body.advance(cx, input), gate).unwrap();
        assert_eq!(finite.next.point(), original.next.point());
        assert_eq!(finite.next.velocity_m_s(), original.next.velocity_m_s());
        assert_eq!(finite.material_heat_j, original.material_heat_j);
        assert_eq!(finite.numerical_loss_j, original.numerical_loss_j);
        assert_eq!(finite.external_work_j, original.external_work_j);
        assert_eq!(finite.transport.report.next_gas, reservoir);
        assert_eq!(finite.transport.report.trials, 1);
        let radiating = LumpedThermalEnvironment {
            radiation_temperature_k: 400.0,
            transport: LumpedThermalTransport::try_declared(10.0, 0.4).unwrap(),
            ..body.clone()
        };
        let with_radiation = FiniteGasThermalEnvironment { body: &radiating, ..environment };
        let result = with_radiation.advance(cx, UniformEnthalpyStepInput {
            curve: source.caloric.curve(), initial: state.point().caloric_state(),
            mass_kg: state.point().mass_kg(), volume_m3: state.point().volume_m3(),
            surface_area_m2: state.geometry().thermal_surface_m2,
            internal_heat_j: 0.0, duration_s: 0.001,
        }).unwrap();
        assert!(result.report.radiation_heat_j > 0.0);
        assert_eq!(result.external_heat_j, result.report.radiation_heat_j);
        assert_eq!(result.report.solid_to_gas_heat_j, 0.0);
        assert_eq!(result.report.gas.supplied_heat_j, 0.0);
        assert_eq!(result.report.next_gas, reservoir);
    });
}

#[test]
fn g4_failed_finite_gas_trials_and_late_cancellation_preserve_both_owners() {
    scope(|cx, gate| {
        let source = sources_with_modulus_slope(0.02, 1000.0, true, true, None, -0.005);
        let initial = Oscillator::try_new(&source.mechanical, &source.caloric,
            geometry(0.0), 300.0, 0.0, 0.5).unwrap();
        let old = initial.step(0.0, 0.00005, config(), adiabatic, gate).unwrap().next;
        let saved = old.clone();
        let reservoir = gas(330.0, 1e-7);
        let saved_gas = reservoir.clone();
        let body = boundary(&source, 100.0);
        let environment = FiniteGasThermalEnvironment {
            gas: &reservoir, body: &body, gas_heat_j: 0.0,
            gas_solve: gas_solve(), exchange: exchange(8e-4),
            maximum_gas_temperature_increment_k: 100.0,
        };
        for failure in 0..4 {
            let mut invalid_body = body.clone();
            let mut invalid = environment.clone();
            match failure {
                0 => invalid.exchange.heat_bracket_j = [-1e-9, 1e-9],
                1 => invalid.gas_heat_j = 1e6,
                2 => invalid_body.convection_w_per_m2_k = 1e9,
                _ => invalid.maximum_gas_temperature_increment_k = 1e-9,
            }
            invalid.body = &invalid_body;
            assert!(old.step(0.0, 0.00005, config(),
                |input| invalid.advance(cx, input), gate).is_err());
            assert_eq!(old.point(), saved.point());
            assert_eq!(old.velocity_m_s(), saved.velocity_m_s());
            assert_eq!(reservoir, saved_gas);
        }
        let mut reached_candidate = false;
        let result = old.step(0.0, 0.00005, config(), |input| {
            let candidate = environment.advance(cx, input)?;
            reached_candidate = true;
            gate.request();
            Ok::<_, fs_couple::finite_gas_thermal::Error>(candidate)
        }, gate);
        assert!(reached_candidate && result.is_err());
        assert_eq!(old.point(), saved.point());
        assert_eq!(old.velocity_m_s(), saved.velocity_m_s());
        assert_eq!(reservoir, saved_gas);
    });
}
