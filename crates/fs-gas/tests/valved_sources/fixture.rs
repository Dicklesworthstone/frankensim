//! Manufactured original NASA/kinetics sources; no empirical material claims.
use super::*;
use fs_evidence::ValidityDomain;
use fs_matdb::{ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, Provenance};
use fs_qty::{DynViscosity, MolarMass, Pressure};
use fs_thermochem::equilibrium::composition::IdealGasEquilibriumV1;
use fs_thermochem::kinetics::{ArrheniusRateV1, ElementaryReactionV1, ReactionProgressRateQuantityV1};
use fs_thermochem::{
    ChargeVector, ElementId, ElementalMatrix, ElementalReferenceIdV1, MolarEnergyQuantityV1,
    Nasa9RegionV1, Nasa9StandardStateModelV1, ReactionId, ReferenceEquationOfStateV1,
    SpeciesId, StandardStateConventionV1, StandardStatePhaseV1,
    UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R, required_nasa9_parameter_dims,
};
use fs_thermochem::transport::kinetic::ThermalConductivityQuantityV1;
use fs_gas::multicomponent::dissipation::{TransportCoefficients, TransportDomain};

fn source(name: &str, formation: f64) -> Nasa9StandardStateModelV1 {
    let mut parameters = std::collections::BTreeMap::new();
    let a = [0.0, 0.0, 3.5, 0.0, 0.0, 0.0, 0.0, formation, 0.0];
    for ((name, dims), value) in required_nasa9_parameter_dims().into_iter().take(9).zip(a) {
        parameters.insert(name.to_owned(), LawParameter { value, dims });
    }
    parameters.insert("reference_pressure".into(), LawParameter { value: 1.0, dims: Pressure::DIMS });
    let card = ConstitutiveModelCard {
        law: LawId("nasa9-standard-state".into()), law_version: 1, parameters,
        state_schema_version: 0, initial_state: InitialStatePolicy::ZeroInternalState,
        validity: ValidityDomain::unconstrained().with("T", 0.001, 50.0), sources: vec![],
        provenance: Provenance { source: "manufactured coupled-valve caloric source".into(),
            license: "fixture".into(), artifact: None },
    };
    Nasa9StandardStateModelV1::new(SpeciesId::new(name).unwrap(), MolarMass::new(R),
        StandardStateConventionV1::new(StandardStatePhaseV1::Gas, ReferenceEquationOfStateV1::IdealGas,
            Pressure::new(1.0), ElementalReferenceIdV1::new("valved-source-fixture").unwrap()).unwrap(),
        vec![Nasa9RegionV1::from_card(card).unwrap()]).unwrap()
}
pub(super) fn mechanism(heat: f64, rate: f64) -> IdealGasKineticsV1 {
    let models = vec![source("A", 2.0), source("B", 2.0 - heat)];
    let ids: Vec<_> = models.iter().map(|m| m.species().clone()).collect();
    let bank = IdealGasEquilibriumV1::new(
        ElementalMatrix::new(vec![ElementId::new("H").unwrap()], ids.clone(), vec![vec![1, 1]]).unwrap(),
        ChargeVector::new(ids, vec![0, 0]).unwrap(), models).unwrap();
    let reaction = ElementaryReactionV1::new(ReactionId::new("valved-isomerization").unwrap(),
        vec![(SpeciesId::new("A").unwrap(), 1)], vec![(SpeciesId::new("B").unwrap(), 1)],
        ArrheniusRateV1::new(Temperature::new(1.0), ReactionProgressRateQuantityV1::new(rate),
            0.0, MolarEnergyQuantityV1::new(0.0)).unwrap()).unwrap();
    IdealGasKineticsV1::new(bank, "valved-source-isomers", 1, vec![reaction], 1e-12).unwrap()
}
pub(super) fn gas(mechanism: &IdealGasKineticsV1) -> MixtureGas {
    MixtureGas::new(mechanism.thermodynamics().clone(),
        CaloricDomain::new(Temperature::new(0.001), Temperature::new(50.0)), || false).unwrap()
}
pub(super) fn chemistry() -> IsochoricKineticsSolve {
    IsochoricKineticsSolve {
        minimum_temperature_k: 0.001, maximum_temperature_k: 50.0,
        energy_tolerance_j: 1e-10, relative_tolerance: 2e-7,
        absolute_amount_tolerance_mol: 1e-14, entropy_tolerance_j_k: 1e-11,
        initial_step_s: 1e-4, maximum_step_s: 1.0, maximum_steps: 30_000,
        ..IsochoricKineticsSolve::default()
    }
}
pub(super) fn law(g: &MixtureGas, strength: f64) -> ConstantTransport {
    ConstantTransport::new(g.source().clone(), "coupled-valve-constant-transport-fixture", 1,
        TransportDomain { minimum_temperature_k: 0.001, maximum_temperature_k: 50.0,
            minimum_pressure_pa: 1e-8, maximum_pressure_pa: 1e8 },
        TransportCoefficients { dynamic_viscosity: DynViscosity::new(0.003 * strength),
            bulk_viscosity: DynViscosity::new(0.001 * strength),
            thermal_conductivity: ThermalConductivityQuantityV1::new(0.02 * strength) },
        vec![0.0, 0.01 * strength, 0.01 * strength, 0.0]).unwrap()
}
pub(super) fn raw(mechanism: &IdealGasKineticsV1) -> MixtureNetwork1d {
    let g = gas(mechanism);
    let branches = [([0.9, 0.1], 2.0), ([0.2, 0.8], 1.0)].into_iter().map(|(d, p)| {
        let states = (0..4).map(|i| g.conserved(&d, 0.0, if i < 2 { 0.9 * p } else { p }).unwrap()).collect();
        MixtureEuler1d::new(g.clone(), Length::new(0.1), 0.1, states,
            Boundaries::Ends { left: Boundary::Reflecting, right: Boundary::Reflecting }, || false).unwrap()
    }).collect();
    MixtureNetwork1d::new(branches, vec![
        Plenum::new(0.03, g.conserved(&[0.85, 0.15], 0.0, 2.1).unwrap()).unwrap(),
        Plenum::new(0.025, g.conserved(&[0.25, 0.75], 0.0, 0.9).unwrap()).unwrap(),
    ], vec![PlenumConnection { branch: 0, end: BranchEnd::Right, plenum: 0 },
        PlenumConnection { branch: 1, end: BranchEnd::Left, plenum: 1 }], || false).unwrap()
}
pub(super) fn orifice(area: f64) -> Orifice {
    Orifice::new(Area::new(area), 0.8, 1.0, ValveDirection::Bidirectional).unwrap()
}
pub(super) fn sources(n: &MixtureNetwork1d) -> (NetworkWallOperator<'static>, NetworkValveOperator) {
    let contact = WallContact { drag: WallDrag::Linear { rate_per_s: 0.4 },
        thermal: WallThermal::Lumped { temperature_k: 3.0, heat_capacity_j_per_k: 0.2, conductance_w_per_k: 0.03 } };
    let walls = NetworkWallOperator::new(n,
        n.branches().iter().map(|b| vec![contact; b.cell_count()]).collect(),
        vec![WallContact { drag: WallDrag::None, ..contact }; n.plenums().len()], || false).unwrap();
    let valves = NetworkValveOperator::new(n,
        vec![ValveLink { first: 0, second: 1, orifice: orifice(0.01) }], || false).unwrap();
    (walls, valves)
}
pub(super) fn molecular(n: &MixtureNetwork1d) -> Diffusion {
    NetworkDissipationOperator::new(n, law(n.gas(), 1.0), vec![Length::new(0.1); n.connections().len()], || false).unwrap()
}
pub(super) fn reacting(heat: f64, rate: f64) -> (ReactingNetwork1d, Diffusion, NetworkWallOperator<'static>, NetworkValveOperator) {
    let m = mechanism(heat, rate);
    let n = raw(&m);
    let (wall, valve) = sources(&n);
    let diffusion = molecular(&n);
    (ReactingNetwork1d::new(n, m, chemistry(), || false).unwrap(), diffusion, wall, valve)
}
