//! Manufactured isomers verify integration, not empirical combustion accuracy.
use super::*;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, ConstitutiveModelCard, InitialStatePolicy, InterpolationPolicy, LawId, LawParameter,
    MaterialCard, MaterialStateId, PropertyClaim, PropertyKey, PropertyValue, Provenance,
    QueryPoint, UncertaintyModel,
};
use fs_material::state_point::MaterialPropertySelection;
use fs_material::thermochemical_gas::SpeciesTransportMaterial;
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{DynViscosity, MolarMass, Pressure, QuantitySpec};
use fs_thermochem::equilibrium::ReactionProgressRateQuantityV1;
use fs_thermochem::equilibrium::composition::IdealGasEquilibriumV1;
use fs_thermochem::kinetics::{ArrheniusRateV1, ElementaryReactionV1};
use fs_thermochem::mixture::FrozenIdealGasMixtureModelV1;
use fs_thermochem::{
    ChargeVector, Composition, CompositionBasis, ElementId, ElementalMatrix,
    ElementalReferenceIdV1, MolarEnergyQuantityV1, NASA9_LAW_ID_V1, NASA9_LAW_VERSION_V1,
    NASA9_STATE_SCHEMA_VERSION_V1, Nasa9RegionV1, Nasa9StandardStateModelV1, ReactionId,
    ReferenceEquationOfStateV1, StandardStateConventionV1, StandardStatePhaseV1,
    UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R, required_nasa9_parameter_dims,
};

fn id(name: &str) -> SpeciesId {
    SpeciesId::new(name).unwrap()
}
fn quantity(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}
fn provenance() -> Provenance {
    Provenance {
        source: "manufactured reversible isomer fixture".into(),
        license: "CC0-1.0".into(),
        artifact: None,
    }
}
fn species(name: &str, formation_j_mol: f64) -> Nasa9StandardStateModelV1 {
    let mut coefficients = [0.0; 9];
    coefficients[2] = 4.0;
    coefficients[7] = formation_j_mol / R;
    let mut parameters = std::collections::BTreeMap::new();
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
    let card = ConstitutiveModelCard {
        law: LawId(NASA9_LAW_ID_V1.into()),
        law_version: NASA9_LAW_VERSION_V1,
        parameters,
        state_schema_version: NASA9_STATE_SCHEMA_VERSION_V1,
        initial_state: InitialStatePolicy::ZeroInternalState,
        validity: ValidityDomain::unconstrained().with("T", 200.0, 2000.0),
        sources: vec![],
        provenance: provenance(),
    };
    Nasa9StandardStateModelV1::new(
        id(name),
        MolarMass::new(0.01),
        StandardStateConventionV1::new(
            StandardStatePhaseV1::Gas,
            ReferenceEquationOfStateV1::IdealGas,
            Pressure::new(100000.0),
            ElementalReferenceIdV1::new("manufactured-isomers").unwrap(),
        )
        .unwrap(),
        vec![Nasa9RegionV1::from_card(card).unwrap()],
    )
    .unwrap()
}
fn transport(name: &str, scale: f64, maximum_temperature: f64) -> SpeciesTransportMaterial {
    let mut claims = ClaimSet::new();
    for (property, dims, value) in [
        ("dynamic_viscosity", DynViscosity::DIMS, 1.7e-5 * scale),
        (
            "thermal_conductivity",
            fs_material::state_point::THERMAL_CONDUCTIVITY_DIMS,
            0.025 * scale,
        ),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(property, dims),
                value: PropertyValue::Curve {
                    abscissa: "temperature".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![(200.0, value), (maximum_temperature, 2.0 * value)],
                    dims,
                },
                validity: ValidityDomain::unconstrained()
                    .with_quantity(
                        "temperature",
                        quantity(QuantityKind::AbsoluteTemperature),
                        200.0,
                        maximum_temperature,
                    )
                    .with_quantity("pressure", quantity(QuantityKind::Pressure), 100.0, 1e8),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: vec![],
                provenance: provenance(),
            })
            .unwrap();
    }
    SpeciesTransportMaterial {
        species: id(name),
        material: MaterialCard::assemble(
            MaterialStateId {
                chemistry: name.into(),
                phase: "gas".into(),
                process: "manufactured".into(),
                revision: 0,
            },
            claims,
            vec![],
        )
        .unwrap(),
        selection: MaterialPropertySelection::SingleClaimOnly,
    }
}
fn fixture(
    formation: f64,
    maximum_transport_temperature: f64,
) -> (ThermochemicalGasModel, IdealGasKineticsV1, ClosedGasParcel) {
    let models = vec![species("A", 0.0), species("B", formation)];
    let ids = vec![id("A"), id("B")];
    let chemistry = IdealGasEquilibriumV1::new(
        ElementalMatrix::new(
            vec![ElementId::new("H").unwrap()],
            ids.clone(),
            vec![vec![1, 1]],
        )
        .unwrap(),
        ChargeVector::new(ids, vec![0, 0]).unwrap(),
        models.clone(),
    )
    .unwrap();
    let rate = ArrheniusRateV1::new(
        Temperature::new(400.0),
        ReactionProgressRateQuantityV1::new(4.0),
        0.0,
        MolarEnergyQuantityV1::new(0.0),
    )
    .unwrap();
    let reaction = ElementaryReactionV1::new(
        ReactionId::new("isomerization").unwrap(),
        vec![(id("A"), 1)],
        vec![(id("B"), 1)],
        rate,
    )
    .unwrap();
    let mechanism = IdealGasKineticsV1::new(
        chemistry,
        "manufactured-isomer-rate",
        1,
        vec![reaction],
        0.0,
    )
    .unwrap();
    let bank = ThermochemicalGasModel::new(
        FrozenIdealGasMixtureModelV1::new(
            models,
            Composition::new(CompositionBasis::MoleFraction, vec![0.5, 0.5]).unwrap(),
        )
        .unwrap(),
        vec![
            transport("A", 1.0, maximum_transport_temperature),
            transport("B", 5.0, maximum_transport_temperature),
        ],
        None,
    )
    .unwrap();
    let initial = bank
        .with_mole_amounts(&[Amount::new(1.0), Amount::new(0.0)])
        .unwrap();
    let point = QueryPoint::new()
        .with_quantity(
            "temperature",
            quantity(QuantityKind::AbsoluteTemperature),
            400.0,
        )
        .unwrap()
        .with_quantity("pressure", quantity(QuantityKind::Pressure), R * 400.0)
        .unwrap()
        .with("sample_id", 7.0)
        .unwrap();
    let parcel = ClosedGasParcel::new(initial, point, 0.01).unwrap();
    (bank, mechanism, parcel)
}
fn controls() -> MaterialKineticsSolve {
    MaterialKineticsSolve {
        chemistry: IsochoricKineticsSolve {
            minimum_temperature_k: 200.0,
            maximum_temperature_k: 1900.0,
            energy_tolerance_j: 1e-8,
            relative_tolerance: 2e-6,
            relative_conservation_tolerance: 1e-10,
            absolute_amount_tolerance_mol: 1e-13,
            entropy_tolerance_j_k: 1e-9,
            initial_step_s: 1e-4,
            maximum_step_s: 1e-2,
            minimum_step_s: 1e-16,
            maximum_steps: 50000,
            maximum_newton_iterations: 24,
            maximum_backtracks: 24,
            maximum_temperature_iterations: 64,
        },
        material_energy: GasEnergySolve {
            minimum_temperature_k: 200.0,
            maximum_temperature_k: 1900.0,
            energy_tolerance_j_kg: 1e-5,
            maximum_iterations: 96,
        },
        mass_tolerance_kg: 1e-13,
        relative_state_tolerance: 1e-9,
    }
}

#[test]
fn g1_actual_elapsed_reaction_matches_analytic_isomer_and_preserves_storage() {
    let (bank, mechanism, mut parcel) = fixture(0.0, 2000.0);
    let before = parcel.clone();
    let report = advance_closed_parcel_adiabatic(
        &mut parcel,
        &bank,
        &mechanism,
        Time::new(0.1),
        controls(),
        || false,
    )
    .unwrap();
    // With equal standard-state energies, kf=kr=4/s and xB=(1-exp(-8t))/2.
    let expected = 0.5 * (1.0 - (-0.8_f64).exp());
    let fractions = parcel
        .model()
        .thermodynamics()
        .mole_composition()
        .fractions();
    assert!((fractions[1] - expected).abs() < 3e-4, "{fractions:?}");
    assert!(fractions[1] < 0.4); // Finite time does not jump to equilibrium0.5.
    assert_eq!(report.chemistry.elapsed_time(), Time::new(0.1));
    assert_eq!(parcel.mass_kg().to_bits(), before.mass_kg().to_bits());
    assert_eq!(parcel.volume_m3().to_bits(), before.volume_m3().to_bits());
    assert_eq!(parcel.point().axes()["sample_id"], 7.0);
    assert!((parcel.state().state().temperature - 400.0).abs() < 1e-6);
    assert!(
        report.balance.energy_residual_j.abs()
            <= parcel.mass_kg() * controls().material_energy.energy_tolerance_j_kg
    );
    assert_eq!(
        parcel.model().thermodynamics().components(),
        bank.thermodynamics().components()
    );
    assert_eq!(
        parcel.model().transport_materials(),
        bank.transport_materials()
    );
    assert!(
        parcel.state().state().thermal_conductivity > before.state().state().thermal_conductivity
    );
}

#[test]
fn g0_formation_energy_drives_heating_pressure_and_common_material_properties() {
    let (bank, mechanism, mut parcel) = fixture(-24000.0, 2000.0);
    let before = parcel.clone();
    let report = advance_closed_parcel_adiabatic(
        &mut parcel,
        &bank,
        &mechanism,
        Time::new(0.05),
        controls(),
        || false,
    )
    .unwrap();
    let gas = parcel.state().state();
    let fraction = parcel
        .model()
        .thermodynamics()
        .mole_composition()
        .fractions()[1];
    let independently_expected = 400.0 + 24000.0 * fraction / (3.0 * R);
    assert!((gas.temperature - independently_expected).abs() < 1e-6);
    assert!(gas.temperature > 500.0);
    assert!(
        (gas.pressure / before.state().state().pressure - gas.temperature / 400.0).abs() < 1e-10
    );
    gas.validate_coherence().unwrap();
    assert!(gas.sound_speed > before.state().state().sound_speed);
    assert!(gas.dynamic_viscosity > before.state().state().dynamic_viscosity);
    assert!(report.balance.freeze_energy_residual_j.abs() < 2e-8);
    assert_eq!(parcel.mass_kg().to_bits(), before.mass_kg().to_bits());
    assert_eq!(parcel.volume_m3().to_bits(), before.volume_m3().to_bits());
}

#[test]
fn g0_zero_time_is_exact_and_source_identity_refusals_are_atomic() {
    let (bank, mechanism, mut parcel) = fixture(0.0, 2000.0);
    let before = parcel.clone();
    let report = advance_closed_parcel_adiabatic(
        &mut parcel,
        &bank,
        &mechanism,
        Time::new(0.0),
        controls(),
        || false,
    )
    .unwrap();
    assert_eq!(parcel, before);
    assert_eq!(report.chemistry.elapsed_time(), Time::new(0.0));
    let wrong_water = ThermochemicalGasModel::new(
        bank.thermodynamics().clone(),
        bank.transport_materials().to_vec(),
        Some(&id("A")),
    )
    .unwrap();
    assert!(matches!(
        advance_closed_parcel_adiabatic(
            &mut parcel,
            &wrong_water,
            &mechanism,
            Time::new(0.01),
            controls(),
            || false
        ),
        Err(MaterialKineticsError::Invalid(_))
    ));
    assert_eq!(parcel, before);
    let wrong_transport = ThermochemicalGasModel::new(
        bank.thermodynamics().clone(),
        vec![transport("A", 2.0, 2000.0), transport("B", 5.0, 2000.0)],
        None,
    )
    .unwrap();
    assert!(matches!(
        advance_closed_parcel_adiabatic(
            &mut parcel,
            &wrong_transport,
            &mechanism,
            Time::new(0.01),
            controls(),
            || false
        ),
        Err(MaterialKineticsError::SourceMismatch(_))
    ));
    assert_eq!(parcel, before);
}

#[test]
fn g4_early_late_cancellation_and_exhausted_work_leave_parcel_unchanged() {
    let (bank, mechanism, original) = fixture(0.0, 2000.0);
    let mut parcel = original.clone();
    assert_eq!(
        advance_closed_parcel_adiabatic(
            &mut parcel,
            &bank,
            &mechanism,
            Time::new(0.03),
            controls(),
            || true
        ),
        Err(MaterialKineticsError::Cancelled)
    );
    assert_eq!(parcel, original);
    let mut polls = 0;
    advance_closed_parcel_adiabatic(
        &mut parcel,
        &bank,
        &mechanism,
        Time::new(0.03),
        controls(),
        || {
            polls += 1;
            false
        },
    )
    .unwrap();
    assert!(polls > 4);
    for stop in [polls / 2, polls] {
        let mut candidate = original.clone();
        let mut count = 0;
        assert_eq!(
            advance_closed_parcel_adiabatic(
                &mut candidate,
                &bank,
                &mechanism,
                Time::new(0.03),
                controls(),
                || {
                    count += 1;
                    count >= stop
                }
            ),
            Err(MaterialKineticsError::Cancelled)
        );
        assert_eq!(candidate, original);
    }
    let mut candidate = original.clone();
    let mut limited = controls();
    limited.chemistry.maximum_steps = 1;
    assert!(
        advance_closed_parcel_adiabatic(
            &mut candidate,
            &bank,
            &mechanism,
            Time::new(0.03),
            limited,
            || false
        )
        .is_err()
    );
    assert_eq!(candidate, original);
}

#[test]
fn g0_final_transport_domain_failure_refuses_the_entire_reaction_update() {
    let (bank, mechanism, mut parcel) = fixture(-24000.0, 450.0);
    let before = parcel.clone();
    assert!(
        advance_closed_parcel_adiabatic(
            &mut parcel,
            &bank,
            &mechanism,
            Time::new(0.05),
            controls(),
            || false
        )
        .is_err()
    );
    assert_eq!(parcel, before);
}
