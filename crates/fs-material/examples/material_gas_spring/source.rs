//! Example-only synthetic input, lowered into EXISTING source-bearing cards.
//! Not a production material schema, empirical gas database or accuracy claim.
use super::Result;
use fs_evidence::ValidityDomain;
use fs_material::state_point::MaterialPropertySelection;
use fs_material::thermochemical_gas::{SpeciesTransportMaterial, ThermochemicalGasModel};
use fs_matdb::{ClaimSet, ConstitutiveModelCard, InitialStatePolicy, InterpolationPolicy,
    LawId, LawParameter, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Dims, DynViscosity, MolarMass, Pressure, QuantitySpec, Temperature};
use fs_thermochem::mixture::FrozenIdealGasMixtureModelV1;
use fs_thermochem::{Composition, CompositionBasis, ElementalReferenceIdV1,
    NASA9_LAW_ID_V1, NASA9_LAW_VERSION_V1, NASA9_STATE_SCHEMA_VERSION_V1,
    Nasa9RegionV1, Nasa9StandardStateModelV1, ReferenceEquationOfStateV1,
    SpeciesId, StandardStateConventionV1, StandardStatePhaseV1, required_nasa9_parameter_dims};

fn typed(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}
pub fn point(t: f64, p: f64) -> Result<QueryPoint> {
    Ok(QueryPoint::new().with_quantity("temperature", typed(QuantityKind::AbsoluteTemperature), t)?
        .with_quantity("pressure", typed(QuantityKind::Pressure), p)?)
}
pub fn load(text: &str, source: &str) -> Result<ThermochemicalGasModel> {
    if text.len() > 65536 { return Err("synthetic mixture file exceeds 64 KiB".into()); }
    let provenance = Provenance { source: format!("SYNTHETIC example: {source}"),
        license: "CC0-1.0".into(), artifact: None };
    let (mut species, mut transport, mut fractions) = (Vec::new(), Vec::new(), Vec::new());
    for (line, row) in text.lines().enumerate() {
        let row = row.split('#').next().unwrap_or("").trim();
        if row.is_empty() { continue; }
        let fields: Vec<_> = row.split_whitespace().collect();
        if fields.len() != 7 || species.len() >= 32 {
            return Err(format!("line {}: seven fields required, at most 32 species", line+1).into());
        }
        let v: Vec<f64> = fields[1..].iter().map(|x| x.parse()).collect::<std::result::Result<_, _>>()?;
        if v.iter().any(|x| !x.is_finite()) || v[0] <= 0.0 || v[1] <= 0.0
            || v[2] < 2.5 || v[3] < 0.0 || v[4] <= 0.0 || v[5] <= 0.0 {
            return Err(format!("line {}: invalid synthetic caloric/transport data", line+1).into());
        }
        let id = SpeciesId::new(fields[0])?;
        let mut coefficients = [0.0; 9]; coefficients[2] = v[2]; coefficients[3] = v[3];
        let mut parameters = std::collections::BTreeMap::new();
        for ((key, dims), value) in required_nasa9_parameter_dims().into_iter().take(9).zip(coefficients) {
            parameters.insert(key.to_owned(), LawParameter { value, dims });
        }
        parameters.insert("reference_pressure".into(), LawParameter { value: 100000.0, dims: Pressure::DIMS });
        let card = ConstitutiveModelCard { law: LawId(NASA9_LAW_ID_V1.into()),
            law_version: NASA9_LAW_VERSION_V1, parameters, state_schema_version: NASA9_STATE_SCHEMA_VERSION_V1,
            initial_state: InitialStatePolicy::ZeroInternalState,
            validity: ValidityDomain::unconstrained().with("T", 200.0, 2000.0),
            sources: vec![], provenance: provenance.clone() };
        let convention = StandardStateConventionV1::new(StandardStatePhaseV1::Gas,
            ReferenceEquationOfStateV1::IdealGas, Pressure::new(100000.0),
            ElementalReferenceIdV1::new("synthetic-gas-spring")?)?;
        species.push(Nasa9StandardStateModelV1::new(id.clone(), MolarMass::new(v[1]), convention,
            vec![Nasa9RegionV1::from_card(card)?])?);
        fractions.push(v[0]);
        let mut claims = ClaimSet::new();
        for (name, dims, value) in [("dynamic_viscosity", DynViscosity::DIMS, v[4]),
            ("thermal_conductivity", Dims([1, 1, -3, -1, 0, 0]), v[5])] {
            claims.insert_claim(PropertyClaim { key: PropertyKey::new(name, dims),
                value: PropertyValue::Curve { abscissa: "temperature".into(), abscissa_dims: Temperature::DIMS,
                    knots: vec![(200.0, value), (2000.0, 2.0*value)], dims },
                validity: ValidityDomain::unconstrained()
                    .with_quantity("temperature", typed(QuantityKind::AbsoluteTemperature), 200.0, 2000.0)
                    .with_quantity("pressure", typed(QuantityKind::Pressure), 1000.0, 1e7),
                uncertainty: UncertaintyModel::Unstated, interpolation: InterpolationPolicy::LinearInside,
                observations: vec![], provenance: provenance.clone() })?;
        }
        transport.push(SpeciesTransportMaterial { species: id,
            material: MaterialCard::assemble(MaterialStateId { chemistry: fields[0].into(),
                phase: "gas".into(), process: "synthetic-gas-spring".into(), revision: 0 }, claims, vec![])?,
            selection: MaterialPropertySelection::SingleClaimOnly });
    }
    if species.is_empty() { return Err("explicit nonempty synthetic mixture required".into()); }
    Ok(ThermochemicalGasModel::new(FrozenIdealGasMixtureModelV1::new(species,
        Composition::new(CompositionBasis::MoleFraction, fractions)?)?, transport, None)?)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn malformed_source_data_refuses_without_an_air_fallback() {
        for text in ["", "# only comment", "A 1 0.028 3.5 0 1e-5",
            "A 1 0.028 3.5 0 NaN 0.02", "A 1 0.028 0.5 0 1e-5 0.02",
            "A 0.5 0.028 3.5 0 1e-5 0.02\nA 0.5 0.028 3.5 0 1e-5 0.02"] {
            assert!(load(text, "synthetic refusal test").is_err());
        }
    }
}
