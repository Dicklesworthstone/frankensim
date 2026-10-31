//! Manufactured canonical source packs, not empirical gas or phase data.
//! Exercise the public pack ingress without a new thermochemistry dependency.
use fs_blake3::{ContentHash, hash_domain};
use fs_evidence::ValidityDomain;
use fs_matdb::*;
use fs_material::state_point::{MaterialPropertySelection, THERMAL_CONDUCTIVITY_DIMS};
use fs_material::thermochemical_gas::{ThermochemicalGasModel, packs::{SpeciesPackRef, from_packs}};
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve};
use fs_qty::{Dims, DynViscosity, Pressure, Temperature, QuantitySpec};
use fs_qty::chemistry::SpeciesId;
use fs_qty::semantic::{Composition, CompositionBasis, QuantityKind, SemanticType, ValueForm};

fn hash(text: &str) -> ContentHash {
    hash_domain("org.frankensim.tests.finite-gas-convection.v1", text.as_bytes())
}
pub fn typed(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}
struct Packs {
    species: NormalizedSpeciesPack,
    nasa9: NormalizedModelPack,
    transport: NormalizedMaterialCardPack,
}
impl Packs {
    fn binding(&self) -> SpeciesPackRef<'_> {
        SpeciesPackRef { species: &self.species, nasa9: Some(&self.nasa9),
            transport: Some((&self.transport, MaterialPropertySelection::SingleClaimOnly)) }
    }
}
fn packs(name: &str, mass: f64, cp_r: f64, mu: f64, k: f64,
    variable_transport: bool, energy_reference_k: f64) -> Packs {
    let artifact = hash(name);
    let provenance = Provenance { source: format!("synthetic {name}; not measured"),
        license: "CC0-1.0".into(), artifact: Some(artifact) };
    let association = SpeciesAssociation::new(SpeciesId::new(name).unwrap(), mass,
        "gas", "ideal-gas", 100000.0, "synthetic-common-reference", vec![artifact], provenance.clone()).unwrap();
    let species = NormalizedSpeciesPack::new(name, "synthetic-convection-v1", artifact, "CC0-1.0",
        association, vec![
            SpeciesNormalizationReceipt::new(SpeciesNormalizationTarget::MolarMass, hash("mass"),
                SPECIES_MOLAR_MASS_DIMS, 1.0, 0.0, "kg/mol", SPECIES_PACK_TARGET_BASIS),
            SpeciesNormalizationReceipt::new(SpeciesNormalizationTarget::ReferencePressure, hash("p0"),
                SPECIES_REFERENCE_PRESSURE_DIMS, 1.0, 0.0, "Pa", SPECIES_PACK_TARGET_BASIS),
        ]).unwrap();
    // Explicit version-one NASA-9 source format, not a replacement evaluator.
    let mut parameters = std::collections::BTreeMap::new();
    for (i, power) in [2, 1, 0, -1, -2, -3, -4, 1, 0].into_iter().enumerate() {
        parameters.insert(format!("a{i}"), LawParameter {
            value: if i == 2 { cp_r } else if i == 7 { energy_reference_k } else { 0.0 },
            dims: Dims([0, 0, 0, power, 0, 0]) });
    }
    parameters.insert("reference_pressure".into(), LawParameter { value: 100000.0, dims: Pressure::DIMS });
    let model = ConstitutiveModelCard { law: LawId("nasa9-standard-state".into()), law_version: 1,
        parameters, state_schema_version: 0, initial_state: InitialStatePolicy::ZeroInternalState,
        validity: ValidityDomain::unconstrained().with("T", 200.0, 1000.0),
        sources: vec![artifact], provenance: provenance.clone() };
    let identity = model.content_hash();
    let mut normalized = model.parameters.iter().map(|(name, p)| ModelNormalizationReceipt::new(
        ModelNormalizationTarget::Parameter { model: identity, parameter: name.clone() },
        hash(&format!("{name} {}", p.value)), p.dims, 1.0, 0.0, "authored SI",
        MODEL_PACK_TARGET_BASIS, None, None)).collect::<Vec<_>>();
    for (side, value) in [(ValidityBoundSide::Lower, 200.0), (ValidityBoundSide::Upper, 1000.0)] {
        normalized.push(ModelNormalizationReceipt::new(
            ModelNormalizationTarget::ValidityBound { model: identity, axis: "T".into(), side },
            hash(&format!("{value} K")), Temperature::DIMS, 1.0, 0.0, "K",
            MODEL_PACK_TARGET_BASIS, None, None));
    }
    let nasa9 = NormalizedModelPack::new(name, "synthetic-convection-v1", artifact,
        "CC0-1.0", vec![model], normalized).unwrap();
    let mut claims = ClaimSet::new();
    let observation = claims.register_observation(ObservationDataset { specimen: name.into(),
        method: "manufactured fixture".into(), artifact, caveats: "not empirical".into(),
        provenance: provenance.clone() }).unwrap();
    for (name, dims, value) in [("dynamic_viscosity", DynViscosity::DIMS, mu),
        ("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS, k)] {
        claims.insert_claim(PropertyClaim { key: PropertyKey::new(name, dims),
            value: if variable_transport {
                PropertyValue::Curve { abscissa: "temperature".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![(200.0, 0.5 * value), (1000.0, 2.5 * value)], dims }
            } else { PropertyValue::Scalar { value, dims } },
            validity: ValidityDomain::unconstrained()
                .with_quantity("temperature", typed(QuantityKind::AbsoluteTemperature), 200.0, 1000.0)
                .with_quantity("pressure", typed(QuantityKind::Pressure), 1000.0, 1e7),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: if variable_transport { InterpolationPolicy::LinearInside }
                else { InterpolationPolicy::ConstantWithinValidity },
            observations: vec![observation.clone()], provenance: provenance.clone() }).unwrap();
    }
    let pack = NormalizedPack::new(name, "synthetic-convection-v1", artifact, "CC0-1.0",
        claims, vec![], vec![]).unwrap();
    let transport = NormalizedMaterialCardPack::new(MaterialStateId { chemistry: name.into(),
        phase: "gas".into(), process: "manufactured".into(), revision: 0 }, pack).unwrap();
    Packs { species, nasa9, transport }
}
pub fn model(fraction: f64) -> ThermochemicalGasModel {
    model_transport(fraction, false)
}
pub fn model_transport(fraction: f64, variable_transport: bool) -> ThermochemicalGasModel {
    model_transport_with_reference(fraction, variable_transport, 0.0)
}
pub fn model_transport_with_reference(
    fraction: f64, variable_transport: bool, energy_reference_k: f64,
) -> ThermochemicalGasModel {
    let a = packs("A", 0.028, 3.5, 1.8e-5, 0.026, variable_transport, energy_reference_k);
    let b = packs("B", 0.004, 2.5, 2.0e-5, 0.08, variable_transport, energy_reference_k);
    from_packs(&[a.binding(), b.binding()], Composition::new(
        CompositionBasis::MoleFraction, vec![fraction, 1.0-fraction]).unwrap(), None).unwrap()
}
pub fn point(t: f64, p: f64) -> QueryPoint {
    QueryPoint::new().with_quantity("temperature", typed(QuantityKind::AbsoluteTemperature), t).unwrap()
        .with_quantity("pressure", typed(QuantityKind::Pressure), p).unwrap()
}
pub fn body_curve() -> EquilibriumEnthalpyPhaseCurve {
    let card = MaterialCard::assemble(MaterialStateId { chemistry: "synthetic body".into(),
        phase: "equilibrium phase".into(), process: "manufactured".into(), revision: 0 },
        ClaimSet::new(), vec![]).unwrap();
    // cp=1000 J/(kg K), latent heat=60000 J/kg at 450 K; fixed density.
    EquilibriumEnthalpyPhaseCurve::try_new(card.content_hash(), [
        (0.0, 300.0, 0.0), (150000.0, 450.0, 0.0),
        (210000.0, 450.0, 1.0), (360000.0, 600.0, 1.0),
    ].map(|(h, t, x)| EnthalpyPhaseKnot { specific_enthalpy_j_kg: h,
        temperature_k: t, liquid_mass_fraction: x, bulk_density_kg_m3: 1000.0 }).to_vec()).unwrap()
}
