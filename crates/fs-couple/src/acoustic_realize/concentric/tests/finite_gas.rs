//! G1/G3/G4: finite source gas closes the layered string's energy path.
//! Synthetic NASA-9, transport and solid cards are manufactured fixtures,
//! not measured instrument or atmosphere data. No reference audio is used.
use super::*;
use crate::acoustic_realize::concentric::finite_gas::{
    ConcentricFiniteGasSolve, ConcentricGasContact,
};
use crate::acoustic_realize::finite_gas::FiniteGasStringSolve;
use crate::cavity_phs::piston::thermal::{ConductanceSample, HeatExchangeSolve};
use fs_material::thermochemical_gas::energy::{GasEnergySolve, IsochoricGasParcel};
use fs_material::thermochemical_gas::{SpeciesTransportMaterial, ThermochemicalGasModel};
use fs_qty::{MolarMass, QuantitySpec, semantic::{QuantityKind, SemanticType, ValueForm}};
use std::cell::Cell;

fn parcel(light_fraction: f64, pressure_pa: f64, temperature_k: f64) -> IsochoricGasParcel {
    use fs_matdb::{ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter};
    use fs_thermochem::mixture::FrozenIdealGasMixtureModelV1;
    use fs_thermochem::{
        Composition, CompositionBasis, ElementalReferenceIdV1, NASA9_LAW_ID_V1,
        NASA9_LAW_VERSION_V1, NASA9_STATE_SCHEMA_VERSION_V1, Nasa9RegionV1,
        Nasa9StandardStateModelV1, ReferenceEquationOfStateV1, SpeciesId,
        StandardStateConventionV1, StandardStatePhaseV1, required_nasa9_parameter_dims,
    };
    let typed = |kind| QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static));
    let provenance = Provenance {
        source: "synthetic two-species layered-string gas".into(),
        license: "CC0-1.0".into(), artifact: None,
    };
    let mut species_models = Vec::new();
    let mut transport = Vec::new();
    for (name, molar_mass, cp_over_r) in [
        ("synthetic-heavy", 0.028, 3.5),
        ("synthetic-light", 0.004, 2.5),
    ] {
        let species = SpeciesId::new(name).unwrap();
        let mut parameters = BTreeMap::new();
        for ((key, dims), value) in required_nasa9_parameter_dims().into_iter().take(9)
            .zip([0.0, 0.0, cp_over_r, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
        {
            parameters.insert(key.into(), LawParameter { value, dims });
        }
        parameters.insert("reference_pressure".into(),
            LawParameter { value: 100000.0, dims: Pressure::DIMS });
        let nasa = ConstitutiveModelCard {
            law: LawId(NASA9_LAW_ID_V1.into()), law_version: NASA9_LAW_VERSION_V1,
            parameters, state_schema_version: NASA9_STATE_SCHEMA_VERSION_V1,
            initial_state: InitialStatePolicy::ZeroInternalState,
            validity: ValidityDomain::unconstrained().with("T", 200.0, 900.0),
            sources: vec![], provenance: provenance.clone(),
        };
        species_models.push(Nasa9StandardStateModelV1::new(
            species.clone(), MolarMass::new(molar_mass),
            StandardStateConventionV1::new(StandardStatePhaseV1::Gas,
                ReferenceEquationOfStateV1::IdealGas, Pressure::new(100000.0),
                ElementalReferenceIdV1::new("synthetic-reference").unwrap()).unwrap(),
            vec![Nasa9RegionV1::from_card(nasa).unwrap()],
        ).unwrap());
        let mut claims = ClaimSet::new();
        for (name, dims, value) in [
            ("dynamic_viscosity", DynViscosity::DIMS, 1.7e-5),
            ("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS, 0.025),
        ] {
            claims.insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Curve {
                    abscissa: "temperature".into(), abscissa_dims: Temperature::DIMS,
                    knots: vec![(200.0, value), (900.0, 2.0*value)], dims,
                },
                validity: ValidityDomain::unconstrained()
                    .with_quantity("temperature", typed(QuantityKind::AbsoluteTemperature), 200.0, 900.0)
                    .with_quantity("pressure", typed(QuantityKind::Pressure), 1000.0, 1e7),
                interpolation: InterpolationPolicy::LinearInside,
                uncertainty: UncertaintyModel::Unstated, observations: vec![],
                provenance: provenance.clone(),
            }).unwrap();
        }
        transport.push(SpeciesTransportMaterial {
            species,
            material: MaterialCard::assemble(MaterialStateId {
                chemistry: name.into(), phase: "gas".into(), process: "synthetic".into(), revision: 0,
            }, claims, vec![]).unwrap(),
            selection: MaterialPropertySelection::SingleClaimOnly,
        });
    }
    let model = ThermochemicalGasModel::new(
        FrozenIdealGasMixtureModelV1::new(species_models,
            Composition::new(CompositionBasis::MoleFraction,
                vec![1.0-light_fraction, light_fraction]).unwrap()).unwrap(),
        transport, None,
    ).unwrap();
    let point = QueryPoint::new()
        .with_quantity("temperature", typed(QuantityKind::AbsoluteTemperature), temperature_k).unwrap()
        .with_quantity("pressure", typed(QuantityKind::Pressure), pressure_pa).unwrap()
        .with("material-lot", 7.0).unwrap();
    IsochoricGasParcel::new(model, point, 0.001).unwrap()
}

fn limits() -> ConcentricFiniteGasSolve {
    ConcentricFiniteGasSolve {
        thermal: controls(),
        gas: FiniteGasStringSolve {
            gas: GasEnergySolve {
                minimum_temperature_k: 200.0, maximum_temperature_k: 800.0,
                energy_tolerance_j_kg: 1e-7, maximum_iterations: 128,
            },
            maximum_gas_temperature_increment_k: 100.0,
        },
        exchange: HeatExchangeSolve {
            heat_bracket_j: [-1e-3, 1e-3], source_tolerance_j: 1e-11,
            energy_tolerance_j: 2e-9, maximum_iterations: 64,
        },
    }
}

fn contact() -> ConcentricGasContact { ConcentricGasContact { area_m2: 0.001 } }

fn runtime<'a>(cx: &Cx<'_>, sources: &'a Sources, gas: &GasState)
    -> (ThermalConcentricStringRuntime<'a>, DepositionSimulation<'a>)
{
    // Reuse the parent's exact material cards, geometry and thermal inventory.
    let (_, thermal) = sources.runtime(cx);
    let cards: Vec<_> = sources.cards.iter().collect();
    let law = ConcentricThermoelasticLaw::try_new(ConcentricSectionTemperatureLaw::try_new(
        sources.reference.clone(), &cards, &["T", "T"]).unwrap()).unwrap();
    let bindings = sources.caloric.iter().enumerate().map(|(i, caloric)|
        ConcentricStringThermalCell { cell: CellId(i as u64), caloric }).collect();
    let live = ThermalConcentricStringRuntime::try_new(cx, template(), law, bindings,
        &thermal, REFERENCE, 0.02, Some(Pluck { station_frac: 0.4, height_m: 0.0005 }),
        *gas, 1.0, RATE).unwrap();
    (live, thermal)
}

fn no_heat(_: &DepositionSimulation<'_>) -> Result<StepCoefficients<()>, &'static str> {
    Ok(StepCoefficients { sources_w: BTreeMap::new(), conductances_w_k: BTreeMap::new(), payload: () })
}

fn supplied_heat(_: &DepositionSimulation<'_>) -> Result<StepCoefficients<()>, &'static str> {
    Ok(StepCoefficients { sources_w: BTreeMap::from([(CellId(0), 16.0), (CellId(1), 48.0)]),
        conductances_w_k: BTreeMap::new(), payload: () })
}

#[test]
fn g1_ohmic_layered_string_and_finite_gas_close_energy_without_double_depositing_drag_or_q() {
    with_cx(|cx, _| {
        let sources = Sources::new(4e-5);
        let light = 0.25;
        let mut gas = parcel(light, 100000.0, 330.0);
        let (mut live, mut thermal) = runtime(cx, &sources, gas.state().state());
        let cv = fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K
            * ((1.0-light)*2.5+light*1.5) / ((1.0-light)*0.028+light*0.004);
        let gas_heat_j = 0.002;
        for sample in 0..8 {
            let old_u = independent_internal(&live, 4e-5);
            let old_h = thermal.energy_j();
            let old_gas = gas.clone();
            let frame = live.step_with_isochoric_gas_exchange(cx, &mut thermal, &mut gas,
                limits(), &[0.002, 0.006], &[0.0; 2], gas_heat_j,
                Some(OhmicDrive::Voltage(2.0)), contact(), supplied_heat, |_, _|
                    Ok(ConductanceSample { conductance_w_k: 0.1, source: () })).unwrap();
            let electrical = frame.string.source.electrical.as_ref().unwrap();
            let supplied_solid = 64.0/f64::from(RATE)+electrical.joule_heat_j;
            near(electrical.joule_heat_j,
                2.0*electrical.electrical.current_a/f64::from(RATE), 2e-12);
            for i in 0..2 {
                let inside = if i == 0 { 0.0_f64 } else { RADII[i-1] };
                let area = core::f64::consts::PI*(RADII[i].powi(2)-inside.powi(2));
                let temperature = thermal.temperature_k(CellId(i as u64)).unwrap();
                let resistivity = (i+1) as f64*1e-6*(1.0+0.002*(temperature-300.0));
                near(electrical.electrical.layers[i].current_a,
                    2.0*area/(resistivity*LENGTH), 2e-12);
            }
            assert!(frame.solid_to_gas_heat_j < 0.0);
            assert!(frame.drag_heat_j > 0.0);
            assert_eq!(frame.drag_heat_j, frame.string.vibration.dissipation.air_loss_j);
            near(frame.external_solid_heat_j, supplied_solid, 2e-12);
            assert_eq!(frame.external_gas_heat_j, gas_heat_j);
            near(frame.string.external_heat_j, supplied_solid-frame.solid_to_gas_heat_j, 2e-12);
            near(frame.string.source.gas.supplied_heat_j,
                gas_heat_j+frame.drag_heat_j+frame.solid_to_gas_heat_j, 2e-12);
            let expected_q = 0.1/f64::from(RATE)
                *(thermal.temperature_k(CellId(1)).unwrap()-gas.state().state().temperature);
            assert!((frame.solid_to_gas_heat_j-expected_q).abs()
                <= limits().exchange.source_tolerance_j);
            let expected_t = old_gas.state().state().temperature
                +(gas_heat_j+frame.drag_heat_j+frame.solid_to_gas_heat_j)/(gas.mass_kg()*cv);
            assert!((gas.state().state().temperature-expected_t).abs()
                <= limits().gas.gas.energy_tolerance_j_kg/cv+64.0*f64::EPSILON*expected_t.abs());
            let gas_storage = gas.mass_kg()*(gas.state().specific_internal_energy_j_kg()
                -old_gas.state().specific_internal_energy_j_kg());
            let loss = &frame.string.vibration.dissipation;
            let defect = frame.string.vibration.acoustic.input_work_j+supplied_solid+gas_heat_j
                -(independent_internal(&live, 4e-5)-old_u+thermal.energy_j()-old_h
                    +gas_storage+loss.authored_loss_j+loss.roundoff_residual_j);
            assert!(defect.abs() <= frame.energy_tolerance_j);
            assert!(frame.energy_balance_residual_j.abs() <= frame.energy_tolerance_j);
            near(frame.string.layers.iter().map(|l| l.material_heat_j).sum(),
                loss.material_heat_j, 2e-12);
            assert_eq!(live.accepted_samples(), sample+1);
            assert_eq!(live.gas(), gas.state().state());
            assert_eq!(&frame.string.source.next_gas, &gas);
            assert_eq!(gas.mass_kg(), old_gas.mass_kg());
            assert_eq!(gas.volume_m3(), old_gas.volume_m3());
            assert_eq!(gas.point().axes().get("material-lot"), Some(&7.0));
        }
    });
}

#[test]
fn g3_composition_and_pressure_data_change_live_drag_motion_and_sound_on_identical_layers() {
    with_cx(|cx, _| {
        let sources = Sources::new(4e-5);
        let run = |fraction, pressure| {
            let mut gas = parcel(fraction, pressure, 300.0);
            let (mut live, mut thermal) = runtime(cx, &sources, gas.state().state());
            let initial_gas = *live.gas();
            let initial_mode = live.modes()[0];
            let initial_motion = live.states().to_vec();
            let geometry = [live.string().length_m, live.string().width_m,
                live.string().lin_density_kg_m, live.string().bending_stiffness_n_m2];
            let mut pressures = Vec::new();
            let mut drag_heat = 0.0;
            for _ in 0..32 {
                let frame = live.step_with_isochoric_gas_exchange(cx, &mut thermal, &mut gas,
                    limits(), &[0.0; 2], &[0.0; 2], 0.0, None, contact(), no_heat, |_, _|
                        Ok(ConductanceSample { conductance_w_k: 0.0, source: () })).unwrap();
                assert_eq!(frame.solid_to_gas_heat_j, 0.0);
                assert_eq!(frame.string.source.gas.supplied_heat_j, frame.drag_heat_j);
                assert_eq!(live.gas(), gas.state().state());
                pressures.push(frame.string.vibration.acoustic.observer_pressure_pa);
                drag_heat += frame.drag_heat_j;
            }
            (initial_gas, initial_mode, initial_motion, geometry, live.states().to_vec(),
                pressures, drag_heat)
        };
        let ordinary = run(0.25, 100000.0);
        let lighter = run(0.75, 100000.0);
        let compressed = run(0.25, 200000.0);
        near(compressed.0.density, 2.0*ordinary.0.density, 2e-13);
        assert!(lighter.0.sound_speed > ordinary.0.sound_speed);
        assert!(compressed.1.damping_ratio > ordinary.1.damping_ratio);
        assert_ne!(lighter.1.damping_ratio, ordinary.1.damping_ratio);
        for altered in [&lighter, &compressed] {
            assert_eq!(altered.1.angular_frequency_rad_s, ordinary.1.angular_frequency_rad_s);
            assert_eq!(altered.2, ordinary.2);
            assert_eq!(altered.3, ordinary.3);
            assert_ne!(altered.4, ordinary.4);
            assert_ne!(altered.5, ordinary.5);
            assert!(altered.5.iter().any(|p| *p != 0.0));
            assert!(altered.6 > 0.0);
        }
        assert!(ordinary.5.iter().any(|p| *p != 0.0));
    });
}

#[test]
fn g4_late_source_refusal_and_cancellation_preserve_three_owners_and_exact_retry() {
    with_cx(|cx, gate| {
        let sources = Sources::new(4e-5);
        let mut accepted_gas = parcel(0.25, 100000.0, 330.0);
        let (mut accepted, mut accepted_thermal) = runtime(cx, &sources, accepted_gas.state().state());
        for _ in 0..32 {
            accepted.step_with_isochoric_gas_exchange(cx, &mut accepted_thermal, &mut accepted_gas,
                limits(), &[0.0; 2], &[0.0; 2], 0.0, None, contact(), no_heat, |_, _|
                    Ok(ConductanceSample { conductance_w_k: 0.1, source: () })).unwrap();
        }
        for cancel in [false, true] {
            let mut live = accepted.clone();
            let mut thermal = accepted_thermal.clone();
            let mut gas = accepted_gas.clone();
            let reached_source = Cell::new(false);
            let refused = live.step_with_isochoric_gas_exchange(cx, &mut thermal, &mut gas,
                limits(), &[0.002, 0.006], &[0.0; 2], 0.002,
                Some(OhmicDrive::Voltage(2.0)), contact(), supplied_heat, |_, _| {
                    reached_source.set(true);
                    if cancel {
                        gate.request();
                        Ok(ConductanceSample { conductance_w_k: 0.1, source: () })
                    } else {
                        Err("synthetic transport unavailable at candidate state")
                    }
                });
            assert!(reached_source.get() && refused.is_err());
            assert_eq!(gas, accepted_gas);
            assert_eq!(thermal.cells(), accepted_thermal.cells());
            assert_eq!(thermal.time_s(), accepted_thermal.time_s());
            assert_eq!(live.states(), accepted.states());
            assert_eq!(live.source(), accepted.source());
            assert_eq!(live.gas(), accepted.gas());
            assert_eq!(live.accepted_samples(), accepted.accepted_samples());
            // A fresh cancellation scope retries the exact same accepted histories.
            with_cx(|retry_cx, _| {
                let mut control = accepted.clone();
                let mut control_thermal = accepted_thermal.clone();
                let mut control_gas = accepted_gas.clone();
                let actual = live.step_with_isochoric_gas_exchange(retry_cx, &mut thermal, &mut gas,
                    limits(), &[0.002, 0.006], &[0.0; 2], 0.002,
                    Some(OhmicDrive::Voltage(2.0)), contact(), supplied_heat, |_, _|
                        Ok(ConductanceSample { conductance_w_k: 0.1, source: () })).unwrap();
                let expected = control.step_with_isochoric_gas_exchange(retry_cx,
                    &mut control_thermal, &mut control_gas, limits(), &[0.002, 0.006], &[0.0; 2],
                    0.002, Some(OhmicDrive::Voltage(2.0)), contact(), supplied_heat, |_, _|
                        Ok(ConductanceSample { conductance_w_k: 0.1, source: () })).unwrap();
                assert_eq!(actual.string.vibration, expected.string.vibration);
                assert_eq!(actual.string.layers, expected.string.layers);
                assert_eq!(actual.solid_to_gas_heat_j, expected.solid_to_gas_heat_j);
                assert_eq!(actual.drag_heat_j, expected.drag_heat_j);
                assert_eq!(actual.energy_balance_residual_j, expected.energy_balance_residual_j);
                assert_eq!(gas, control_gas);
                assert_eq!(thermal.cells(), control_thermal.cells());
                assert_eq!(live.states(), control.states());
                assert_eq!(live.source(), control.source());
                assert_eq!(live.gas(), control.gas());
                assert_eq!(live.accepted_samples(), control.accepted_samples());
            });
        }
    });
}

mod flow;
