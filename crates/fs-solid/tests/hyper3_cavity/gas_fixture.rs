//! Original synthetic NASA-9 and transport sources; no measured-gas claim.
use fs_evidence::ValidityDomain;
use fs_material::state_point::{MaterialPropertySelection, THERMAL_CONDUCTIVITY_DIMS};
use fs_material::thermochemical_gas::{SpeciesTransportMaterial, ThermochemicalGasModel};
use fs_material::thermochemical_gas::energy::parcel::ClosedGasParcel;
use fs_matdb::{
    ClaimSet, ConstitutiveModelCard, InitialStatePolicy, InterpolationPolicy, LawId,
    LawParameter, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey, PropertyValue,
    Provenance, QueryPoint, UncertaintyModel,
};
use fs_qty::{DynViscosity, MolarMass, Pressure, QuantitySpec, Temperature};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_thermochem::{
    Composition, CompositionBasis, ElementalReferenceIdV1, Nasa9RegionV1,
    Nasa9StandardStateModelV1, ReferenceEquationOfStateV1, SpeciesId,
    StandardStateConventionV1, StandardStatePhaseV1, NASA9_LAW_ID_V1,
    NASA9_LAW_VERSION_V1, NASA9_STATE_SCHEMA_VERSION_V1,
    UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R, required_nasa9_parameter_dims,
};
use fs_thermochem::mixture::FrozenIdealGasMixtureModelV1;

fn typed(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}

pub(super) fn gas(
    volume_m3: f64, pressure_pa: f64, cp_slope: f64, energy_datum_j_kg: f64,
) -> ClosedGasParcel {
    let provenance = Provenance {
        source: "synthetic sealed-cavity source".into(), license: "CC0-1.0".into(), artifact: None,
    };
    let coefficients = [0.0, 0.0, 3.5, cp_slope, 0.0, 0.0, 0.0,
        0.028 * energy_datum_j_kg / R, 0.0];
    let mut parameters: std::collections::BTreeMap<_, _> = required_nasa9_parameter_dims()
        .into_iter().take(9).zip(coefficients)
        .map(|((name, dims), value)| (name.to_owned(), LawParameter { value, dims })).collect();
    parameters.insert("reference_pressure".into(), LawParameter { value: 100000.0, dims: Pressure::DIMS });
    let species = SpeciesId::new("synthetic-cavity-gas").unwrap();
    let standard = Nasa9StandardStateModelV1::new(
        species.clone(), MolarMass::new(0.028),
        StandardStateConventionV1::new(StandardStatePhaseV1::Gas, ReferenceEquationOfStateV1::IdealGas,
            Pressure::new(100000.0), ElementalReferenceIdV1::new("synthetic-cavity-reference").unwrap()).unwrap(),
        vec![Nasa9RegionV1::from_card(ConstitutiveModelCard {
            law: LawId(NASA9_LAW_ID_V1.into()), law_version: NASA9_LAW_VERSION_V1,
            parameters, state_schema_version: NASA9_STATE_SCHEMA_VERSION_V1,
            initial_state: InitialStatePolicy::ZeroInternalState,
            validity: ValidityDomain::unconstrained().with("T", 150.0, 2000.0),
            sources: vec![], provenance: provenance.clone(),
        }).unwrap()],
    ).unwrap();
    let mut claims = ClaimSet::new();
    for (name, dims, value) in [
        ("dynamic_viscosity", DynViscosity::DIMS, 1.7e-5),
        ("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS, 0.025),
    ] {
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(name, dims),
            value: PropertyValue::Curve { abscissa: "temperature".into(),
                abscissa_dims: Temperature::DIMS, knots: vec![(150.0, value), (2000.0, 2.0 * value)], dims },
            validity: ValidityDomain::unconstrained()
                .with_quantity("temperature", typed(QuantityKind::AbsoluteTemperature), 150.0, 2000.0)
                .with_quantity("pressure", typed(QuantityKind::Pressure), 1e-2, 1e12),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::LinearInside,
            observations: vec![], provenance: provenance.clone(),
        }).unwrap();
    }
    let transport = SpeciesTransportMaterial {
        species,
        material: MaterialCard::assemble(MaterialStateId {
            chemistry: "synthetic-cavity-gas".into(), phase: "gas".into(),
            process: "synthetic".into(), revision: 0,
        }, claims, vec![]).unwrap(),
        selection: MaterialPropertySelection::SingleClaimOnly,
    };
    let model = ThermochemicalGasModel::new(
        FrozenIdealGasMixtureModelV1::new(vec![standard],
            Composition::new(CompositionBasis::MoleFraction, vec![1.0]).unwrap()).unwrap(),
        vec![transport], None,
    ).unwrap();
    ClosedGasParcel::new(model, QueryPoint::new()
        .with_quantity("temperature", typed(QuantityKind::AbsoluteTemperature), 300.0).unwrap()
        .with_quantity("pressure", typed(QuantityKind::Pressure), pressure_pa).unwrap()
        .with("material-lot", 47.0).unwrap(), volume_m3).unwrap()
}
