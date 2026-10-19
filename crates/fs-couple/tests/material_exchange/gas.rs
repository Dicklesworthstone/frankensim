//! Synthetic pure gases through the original NASA-9/transport/pressure-work owners.
use fs_couple::cavity_phs::piston::{GasPiston,PistonInput,PistonSolve};
use fs_couple::cavity_phs::piston::thermal::HeatExchangeSolve;
use fs_couple::cavity_phs::piston::thermal::boundary::ports::network::{Chamber,Contact};
use fs_couple::cavity_phs::piston::thermal::boundary::ports::network::material::MaterialContact;
use fs_couple::cavity_phs::piston::thermal::boundary::network::material::GasConductionGeometry;
use fs_evidence::ValidityDomain;
use fs_material::state_point::{MaterialPropertySelection,THERMAL_CONDUCTIVITY_DIMS};
use fs_material::thermochemical_gas::{SpeciesTransportMaterial,ThermochemicalGasModel};
use fs_material::thermochemical_gas::energy::{GasEnergySolve,parcel::ClosedGasParcel};
use fs_matdb::*;
use fs_qty::{DynViscosity,MolarMass,Pressure,QuantitySpec};
use fs_qty::semantic::{QuantityKind,SemanticType,ValueForm};
use fs_thermochem::*;
use fs_thermochem::mixture::FrozenIdealGasMixtureModelV1;
use fs_thermal::enthalpy::CellId;
fn quantity(k:QuantityKind)->QuantitySpec {QuantitySpec::semantic(SemanticType::new(k,ValueForm::Static))}
fn gas(mass:f64)->ClosedGasParcel {
    let id=SpeciesId::new("A").unwrap();
    let mut coefficients=[0.0;9];coefficients[2]=3.5;
    let mut parameters=required_nasa9_parameter_dims().into_iter().take(9).zip(coefficients)
        .map(|((key,dims),value)|(key.to_owned(),LawParameter {value,dims})).collect::<std::collections::BTreeMap<_,_>>();
    parameters.insert("reference_pressure".into(),LawParameter {value:100000.0,dims:Pressure::DIMS});
    let species=Nasa9StandardStateModelV1::new(id.clone(),MolarMass::new(mass),
        StandardStateConventionV1::new(StandardStatePhaseV1::Gas,ReferenceEquationOfStateV1::IdealGas,Pressure::new(100000.0),
            ElementalReferenceIdV1::new("synthetic").unwrap()).unwrap(),
        vec![Nasa9RegionV1::from_card(ConstitutiveModelCard {law:LawId(NASA9_LAW_ID_V1.into()),law_version:NASA9_LAW_VERSION_V1,
            parameters,state_schema_version:NASA9_STATE_SCHEMA_VERSION_V1,initial_state:InitialStatePolicy::ZeroInternalState,
            validity:ValidityDomain::unconstrained().with("T",200.0,2000.0),sources:vec![],
            provenance:Provenance {source:"synthetic caloric gas".into(),license:"CC0-1.0".into(),artifact:None}}).unwrap()]).unwrap();
    let mut claims=ClaimSet::new();
    for (key,dims,value) in [("dynamic_viscosity",DynViscosity::DIMS,2e-5),("thermal_conductivity",THERMAL_CONDUCTIVITY_DIMS,0.03)] {
        claims.insert_claim(PropertyClaim {key:PropertyKey::new(key,dims),value:PropertyValue::Scalar {value,dims},
            validity:ValidityDomain::unconstrained().with_quantity("temperature",quantity(QuantityKind::AbsoluteTemperature),200.0,2000.0)
                .with_quantity("pressure",quantity(QuantityKind::Pressure),1000.0,1e7),
            uncertainty:UncertaintyModel::Unstated,interpolation:InterpolationPolicy::ConstantWithinValidity,observations:vec![],
            provenance:Provenance {source:"synthetic transport".into(),license:"CC0-1.0".into(),artifact:None}}).unwrap();
    }
    let transport=SpeciesTransportMaterial {species:id,material:MaterialCard::assemble(MaterialStateId {
        chemistry:"synthetic gas".into(),phase:"gas".into(),process:"authored".into(),revision:0},claims,vec![]).unwrap(),
        selection:MaterialPropertySelection::SingleClaimOnly};
    let model=ThermochemicalGasModel::new(FrozenIdealGasMixtureModelV1::new(vec![species],
        Composition::new(CompositionBasis::MoleFraction,vec![1.0]).unwrap()).unwrap(),vec![transport],None).unwrap();
    let point=QueryPoint::new().with_quantity("temperature",quantity(QuantityKind::AbsoluteTemperature),400.0).unwrap()
        .with_quantity("pressure",quantity(QuantityKind::Pressure),100000.0).unwrap();
    ClosedGasParcel::new(model,point,0.001).unwrap()
}
pub fn pistons()->[GasPiston;2] {[0.028,0.020].map(|m|GasPiston::new(gas(m),1.0,0.01,0.0).unwrap())}
pub fn chamber(id:u64,piston:&GasPiston,dt:f64)->Chamber<'_> {
    Chamber {id,piston,input:PistonInput {duration_s:dt,heat_j:0.0,exterior_pressure_pa:100000.0,external_force_n:0.0},
        piston_solve:PistonSolve {volume_bracket_m3:[0.0007,0.0015],gas:GasEnergySolve {
            minimum_temperature_k:200.0,maximum_temperature_k:1800.0,energy_tolerance_j_kg:1e-7,maximum_iterations:128},
            momentum_tolerance_n_s:1e-11,energy_tolerance_j:1e-7,maximum_iterations:128},
        exchange:HeatExchangeSolve {heat_bracket_j:[-1.0,1.0],source_tolerance_j:1e-9,energy_tolerance_j:1e-7,maximum_iterations:128}}
}
pub fn patches()->[MaterialContact;2] {
    [(101,10,CellId(30)),(202,20,CellId(10))].map(|(id,chamber,cell)|MaterialContact {
        contact:Contact {id,chamber,cell},geometry:GasConductionGeometry {area_m2:0.1,gas_distance_m:0.001,cell_distance_m:0.02}})
}
