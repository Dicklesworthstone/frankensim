//! G1/G3/G4: source gas, live acoustic state and finite gas energy.
//! Synthetic NASA-9 and transport data check coupling, not measured gas fidelity.

use super::*;
use fs_couple::acoustic_realize::finite_gas::{FiniteGasStringHeat, FiniteGasStringSolve};
use fs_material::thermochemical_gas::energy::{GasEnergySolve, IsochoricGasParcel};
use fs_material::thermochemical_gas::{SpeciesTransportMaterial, ThermochemicalGasModel};

fn parcel(molar_mass: f64, pressure: f64) -> IsochoricGasParcel {
    parcel_named("synthetic-string-gas", molar_mass, pressure)
}

fn parcel_named(name: &str, molar_mass: f64, pressure: f64) -> IsochoricGasParcel {
    use fs_matdb::{ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter};
    use fs_qty::{MolarMass, Temperature};
    use fs_thermochem::mixture::FrozenIdealGasMixtureModelV1;
    use fs_thermochem::{
        Composition, CompositionBasis, ElementalReferenceIdV1, NASA9_LAW_ID_V1,
        NASA9_LAW_VERSION_V1, NASA9_STATE_SCHEMA_VERSION_V1, Nasa9RegionV1,
        Nasa9StandardStateModelV1, ReferenceEquationOfStateV1, SpeciesId,
        StandardStateConventionV1, StandardStatePhaseV1, required_nasa9_parameter_dims,
    };
    let typed = |kind| QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static));
    let provenance = Provenance {
        source: "synthetic finite-gas string fixture".into(),
        license: "CC0-1.0".into(),
        artifact: None,
    };
    let mut parameters = std::collections::BTreeMap::new();
    let coefficients = [0.0, 0.0, 3.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
    for ((key, dims), value) in required_nasa9_parameter_dims()
        .into_iter()
        .take(9)
        .zip(coefficients)
    {
        parameters.insert(key.to_owned(), LawParameter { value, dims });
    }
    parameters.insert(
        "reference_pressure".into(),
        LawParameter {
            value: 100000.0,
            dims: Pressure::DIMS,
        },
    );
    let nasa = ConstitutiveModelCard {
        law: LawId(NASA9_LAW_ID_V1.into()),
        law_version: NASA9_LAW_VERSION_V1,
        parameters,
        state_schema_version: NASA9_STATE_SCHEMA_VERSION_V1,
        initial_state: InitialStatePolicy::ZeroInternalState,
        validity: ValidityDomain::unconstrained().with("T", 200.0, 900.0),
        sources: vec![],
        provenance: provenance.clone(),
    };
    let species = SpeciesId::new(name).unwrap();
    let thermo = Nasa9StandardStateModelV1::new(
        species.clone(),
        MolarMass::new(molar_mass),
        StandardStateConventionV1::new(
            StandardStatePhaseV1::Gas,
            ReferenceEquationOfStateV1::IdealGas,
            Pressure::new(100000.0),
            ElementalReferenceIdV1::new("synthetic-reference").unwrap(),
        )
        .unwrap(),
        vec![Nasa9RegionV1::from_card(nasa).unwrap()],
    )
    .unwrap();
    let mut claims = ClaimSet::new();
    for (name, dims, value) in [
        ("dynamic_viscosity", DynViscosity::DIMS, 1.7e-5),
        (
            "thermal_conductivity",
            fs_material::state_point::THERMAL_CONDUCTIVITY_DIMS,
            0.025,
        ),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Curve {
                    abscissa: "temperature".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![(200.0, value), (900.0, 2.0 * value)],
                    dims,
                },
                validity: ValidityDomain::unconstrained()
                    .with_quantity(
                        "temperature",
                        typed(QuantityKind::AbsoluteTemperature),
                        200.0,
                        900.0,
                    )
                    .with_quantity("pressure", typed(QuantityKind::Pressure), 1000.0, 1e7),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: vec![],
                provenance: provenance.clone(),
            })
            .unwrap();
    }
    let transport = SpeciesTransportMaterial {
        species,
        material: MaterialCard::assemble(
            MaterialStateId {
                chemistry: name.into(),
                phase: "gas".into(),
                process: "synthetic".into(),
                revision: 0,
            },
            claims,
            vec![],
        )
        .unwrap(),
        selection: MaterialPropertySelection::SingleClaimOnly,
    };
    let model = ThermochemicalGasModel::new(
        FrozenIdealGasMixtureModelV1::new(
            vec![thermo],
            Composition::new(CompositionBasis::MoleFraction, vec![1.0]).unwrap(),
        )
        .unwrap(),
        vec![transport],
        None,
    )
    .unwrap();
    let point = QueryPoint::new()
        .with_quantity(
            "temperature",
            typed(QuantityKind::AbsoluteTemperature),
            300.0,
        )
        .unwrap()
        .with_quantity("pressure", typed(QuantityKind::Pressure), pressure)
        .unwrap()
        .with("material-lot", 7.0)
        .unwrap();
    IsochoricGasParcel::new(model, point, 0.001).unwrap()
}

fn solve() -> FiniteGasStringSolve {
    FiniteGasStringSolve {
        gas: GasEnergySolve {
            minimum_temperature_k: 200.0,
            maximum_temperature_k: 800.0,
            energy_tolerance_j_kg: 1e-7,
            maximum_iterations: 128,
        },
        maximum_gas_temperature_increment_k: 100.0,
    }
}

fn runtime(cx: &fs_exec::Cx<'_>, gas: &fs_material::gas::GasState) -> ThermalMaterialStringRuntime {
    let (mechanical, card, curve) =
        heated_string_parts(cx, 48000, false, InterpolationPolicy::LinearInside);
    let bound = mechanical.specimen();
    let specimen = with_uniform_circular_material_and_constraints(
        PrestressedString {
            n_modes: 2,
            ..bound.string()
        },
        bound.geometry_constraint(),
        bound.material(),
        bound.prestress(),
    )
    .unwrap()
    .with_kelvin_voigt_bending_loss()
    .unwrap();
    let mechanical = LinearMaterialStringRuntime::try_new_with_gas(
        cx,
        specimen,
        Some(Pluck {
            station_frac: 0.4,
            height_m: 1e-4,
        }),
        gas,
        1.0,
        48000,
    )
    .unwrap();
    ThermalMaterialStringRuntime::try_new(cx, mechanical, card, curve, 0.0, 5.0).unwrap()
}

fn energy(runtime: &ThermalMaterialStringRuntime) -> f64 {
    runtime
        .mechanical()
        .modes()
        .iter()
        .zip(runtime.mechanical().states())
        .map(|(m, s)| {
            0.5 * (s.velocity_m_sqrt_kg_per_s.powi(2)
                + (m.angular_frequency_rad_s * s.displacement_m_sqrt_kg).powi(2))
        })
        .sum()
}

fn accelerations(runtime: &ThermalMaterialStringRuntime) -> Vec<f64> {
    runtime
        .mechanical()
        .modes()
        .iter()
        .zip(runtime.mechanical().states())
        .map(|(m, s)| {
            -m.angular_frequency_rad_s.powi(2) * s.displacement_m_sqrt_kg
                - 2.0 * m.damping_ratio * m.angular_frequency_rad_s * s.velocity_m_sqrt_kg_per_s
        })
        .collect()
}

#[test]
fn g3_live_source_gas_rebind_preserves_motion_and_changes_physical_drag() {
    with_string_cx(|cx, _| {
        let cold = parcel(0.028, 100000.0);
        let other = parcel(0.004, 200000.0);
        let mut live = runtime(cx, cold.state().state());
        live.step(cx, 0, &[0.0; 2], 0.0).unwrap();
        let before = live.clone();
        let update = live.rebind_gas(cx, 1, other.state().state()).unwrap();
        assert_eq!(update.parameter_work_j, 0.0);
        assert_eq!(live.mechanical().states(), before.mechanical().states());
        assert_eq!(live.mechanical().specimen(), before.mechanical().specimen());
        assert_eq!(live.thermal(), before.thermal());
        assert_eq!(live.mechanical().accepted_samples(), 1);
        let string = live.mechanical().specimen().string();
        for (old, new) in before
            .mechanical()
            .modes()
            .iter()
            .zip(live.mechanical().modes())
        {
            assert_eq!(old.angular_frequency_rad_s, new.angular_frequency_rad_s);
            let drag = |gas: &fs_material::gas::GasState| {
                let pi = core::f64::consts::PI;
                2.0 * pi * gas.dynamic_viscosity
                    + 2.0
                        * pi
                        * string.width_m
                        * (gas.dynamic_viscosity * gas.density * new.angular_frequency_rad_s / 2.0)
                            .sqrt()
            };
            let expected = (drag(other.state().state()) - drag(cold.state().state()))
                / (2.0 * string.lin_density_kg_m * new.angular_frequency_rad_s);
            assert!(
                ((new.damping_ratio - old.damping_ratio) - expected).abs()
                    <= 2e-11 * expected.abs()
            );
        }
        let mut control = live.clone();
        let mut incoherent = *other.state().state();
        incoherent.pressure *= 2.0;
        assert!(live.rebind_gas(cx, 2, &incoherent).is_err());
        assert!(live.rebind_gas(cx, 1, cold.state().state()).is_err());
        assert_eq!(
            live.step(cx, 2, &[0.0; 2], 0.0).unwrap(),
            control.step(cx, 2, &[0.0; 2], 0.0).unwrap()
        );
    });
}

#[test]
fn g1_finite_gas_energy_and_endpoint_pressure_use_the_same_accepted_medium() {
    with_string_cx(|cx, _| {
        let mut gas = parcel(0.028, 100000.0);
        let mut live = runtime(cx, gas.state().state());
        for epoch in 0..64 {
            live.step_with_isochoric_gas(
                cx,
                epoch,
                &[0.0; 2],
                &mut gas,
                FiniteGasStringHeat::default(),
                solve(),
            )
            .unwrap();
        }
        let before = live.clone();
        let old_gas = gas.clone();
        let previous = accelerations(&before);
        let cv = 2.5 * fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / 0.028;
        let heat = FiniteGasStringHeat {
            string_heat_j: 5e-8,
            gas_heat_j: gas.mass_kg() * cv * 25.0,
        };
        let frame = live
            .step_with_isochoric_gas(cx, 64, &[0.0; 2], &mut gas, heat, solve())
            .unwrap();
        assert_eq!(live.mechanical().epoch(), 65);
        assert_eq!(live.mechanical().accepted_samples(), 65);
        assert_eq!(live.mechanical().ambient(), gas.state().state());
        assert_eq!(gas.mass_kg(), old_gas.mass_kg());
        assert_eq!(gas.volume_m3(), old_gas.volume_m3());
        assert_eq!(gas.point().axes().get("material-lot"), Some(&7.0));
        assert_eq!(
            gas.state().components()[0].card_identity(),
            old_gas.state().components()[0].card_identity()
        );
        assert!(frame.drag_heat_j > 0.0);
        assert_eq!(
            frame.drag_heat_j,
            frame.string.vibration.dissipation.air_loss_j
        );
        assert_eq!(
            frame.gas.supplied_heat_j,
            heat.gas_heat_j + frame.drag_heat_j
        );
        assert_eq!(frame.string.external_heat_j, heat.string_heat_j);
        let expected_temperature = old_gas.state().state().temperature
            + (heat.gas_heat_j + frame.drag_heat_j) / (gas.mass_kg() * cv);
        let temperature_tolerance = solve().gas.energy_tolerance_j_kg / cv
            + 32.0 * f64::EPSILON * expected_temperature.abs();
        assert!(
            (gas.state().state().temperature - expected_temperature).abs() <= temperature_tolerance
        );
        close(
            gas.state().state().pressure / old_gas.state().state().pressure,
            gas.state().state().temperature / old_gas.state().state().temperature,
        );

        let solid = live.mechanical().specimen().mass_kg()
            * (live.thermal().specific_enthalpy_j_kg() - before.thermal().specific_enthalpy_j_kg());
        let gas_storage = gas.mass_kg()
            * (gas.state().specific_internal_energy_j_kg()
                - old_gas.state().specific_internal_energy_j_kg());
        let defect = heat.string_heat_j + heat.gas_heat_j
            - ((energy(&live) - energy(&before))
                + solid
                + gas_storage
                + frame.string.vibration.dissipation.authored_loss_j
                + frame.string.vibration.dissipation.roundoff_residual_j);
        assert!(defect.abs() <= frame.energy_tolerance_j);
        assert!(frame.energy_balance_residual_j.abs() <= frame.energy_tolerance_j);
        let mut frozen = before.clone();
        frozen.step(cx, 64, &[0.0; 2], heat.string_heat_j).unwrap();
        assert_eq!(live.mechanical().states(), frozen.mechanical().states());

        let string = live.mechanical().specimen().string();
        let pi = core::f64::consts::PI;
        let acceleration = accelerations(&live);
        let pressure = |medium: &fs_material::gas::GasState| {
            let factor = medium.density * string.width_m
                / (4.0 * pi * (string.lin_density_kg_m * string.length_m / 2.0).sqrt());
            factor * 2.0 * string.length_m / pi * acceleration[0]
                - factor * string.length_m.powi(2) / (2.0 * pi * medium.sound_speed)
                    * (acceleration[1] - previous[1])
                    * 48000.0
        };
        close(
            frame.string.vibration.acoustic.observer_pressure_pa,
            pressure(gas.state().state()),
        );
        assert!((pressure(gas.state().state()) - pressure(old_gas.state().state())).abs() > 1e-12);
    });
}

#[test]
fn g4_finite_gas_late_refusal_epoch_and_cancellation_preserve_both_owners() {
    with_string_cx(|cx, gate| {
        let gas = parcel(0.028, 100000.0);
        let mut pristine = runtime(cx, gas.state().state());
        for epoch in 0..64 {
            pristine.step(cx, epoch, &[0.0; 2], 0.0).unwrap();
        }
        for (heat, limits) in [
            (
                FiniteGasStringHeat {
                    string_heat_j: 0.0,
                    gas_heat_j: 512.0,
                },
                solve(),
            ),
            (
                FiniteGasStringHeat {
                    string_heat_j: 0.0,
                    gas_heat_j: 1.0,
                },
                FiniteGasStringSolve {
                    maximum_gas_temperature_increment_k: 1e-6,
                    ..solve()
                },
            ),
            (
                FiniteGasStringHeat::default(),
                FiniteGasStringSolve {
                    gas: GasEnergySolve {
                        maximum_iterations: 0,
                        ..solve().gas
                    },
                    ..solve()
                },
            ),
        ] {
            let mut live = pristine.clone();
            let mut candidate_gas = gas.clone();
            assert!(
                live.step_with_isochoric_gas(cx, 64, &[0.0; 2], &mut candidate_gas, heat, limits)
                    .is_err()
            );
            assert_eq!(candidate_gas, gas);
            assert_eq!(live.mechanical().states(), pristine.mechanical().states());
            assert_eq!(live.mechanical().ambient(), pristine.mechanical().ambient());
            assert_eq!(live.thermal(), pristine.thermal());
            assert_eq!(live.mechanical().epoch(), 64);
            let mut control = pristine.clone();
            assert_eq!(
                live.step(cx, 64, &[0.0; 2], 0.0).unwrap(),
                control.step(cx, 64, &[0.0; 2], 0.0).unwrap()
            );
        }
        let mut live = pristine.clone();
        let mut other_gas = parcel(0.028, 200000.0);
        let other_saved = other_gas.clone();
        assert!(
            live.step_with_isochoric_gas(
                cx,
                64,
                &[0.0; 2],
                &mut other_gas,
                FiniteGasStringHeat::default(),
                solve()
            )
            .is_err()
        );
        assert_eq!(other_gas, other_saved);
        let mut same_gas = gas.clone();
        assert!(
            live.step_with_isochoric_gas(
                cx,
                65,
                &[0.0; 2],
                &mut same_gas,
                FiniteGasStringHeat::default(),
                solve()
            )
            .is_err()
        );
        gate.request();
        assert!(
            live.step_with_isochoric_gas(
                cx,
                64,
                &[0.0; 2],
                &mut same_gas,
                FiniteGasStringHeat::default(),
                solve()
            )
            .is_err()
        );
        assert_eq!(same_gas, gas);
        assert_eq!(live.mechanical().states(), pristine.mechanical().states());
        assert_eq!(live.thermal(), pristine.thermal());
        assert_eq!(live.mechanical().epoch(), 64);
    });
}

#[path = "finite_gas/electrical.rs"]
mod electrical;
#[path = "finite_gas/exchange.rs"]
mod exchange;
#[path = "finite_gas/flow.rs"]
mod flow;
#[path = "finite_gas/reactive.rs"]
mod reactive;
