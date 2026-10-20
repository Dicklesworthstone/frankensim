//! Existing complete assembly fixtures; actual source-backed gas transport.
#[path = "flow/snapshot.rs"]
mod snapshot_tests;

use super::*;
use fs_couple::cavity_phs::piston::thermal::boundary::ports::network::material::state::flow::*;
use fs_material::thermochemical_gas::energy::parcel::transfer::{GasMassTransferSolve,
    nozzle::{GasNozzle, GasNozzleStepSolve, network::{self, GasNetworkLink, GasNetworkSolve}}};
use fs_ports::ModuleWindowAudit;

fn source_gas(name: &str, molar_mass: f64, pressure: f64) -> GasPiston {
    let (species, transport) = gas_fixture::component(name, molar_mass, 0.0);
    let model = ThermochemicalGasModel::new(FrozenIdealGasMixtureModelV1::new(vec![species],
        Composition::new(CompositionBasis::MoleFraction, vec![1.0]).unwrap()).unwrap(), vec![transport], None).unwrap();
    parcel_piston(model, pressure)
}
fn parcel_piston(model: ThermochemicalGasModel, pressure: f64) -> GasPiston {
    let q = |kind| QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static));
    let point = QueryPoint::new().with_quantity("temperature", q(QuantityKind::AbsoluteTemperature), 400.0).unwrap()
        .with_quantity("pressure", q(QuantityKind::Pressure), pressure).unwrap();
    GasPiston::new(ClosedGasParcel::new(model, point, 0.001).unwrap(), 1.0, 0.01, 0.0).unwrap()
}
fn gases(molar_mass: f64) -> BTreeMap<u64, GasPiston> {
    BTreeMap::from([(10, source_gas("A", 0.028, 100000.0)), (20, source_gas("B", molar_mass, 180000.0))])
}
fn drives(gases: &BTreeMap<u64, GasPiston>) -> Vec<ChamberDrive> {
    gases.iter().map(|(&id, gas)| ChamberDrive { exterior_pressure_pa: gas.gas().state().state().pressure,
        ..loads(id, gas) }).collect()
}
fn flow_controls() -> GasNetworkSolve {
    GasNetworkSolve { maximum_parcels: 2, maximum_links: 2, maximum_withdrawn_fraction: 0.05,
        mass_tolerance_kg: 1e-13, energy_tolerance_j: 1e-7 }
}
fn links() -> [GasNetworkLink; 1] {
    let energy = chemical().gas;
    [GasNetworkLink { id: 901, parcels: [20,10], nozzle: GasNozzle { area_m2: 2e-6, discharge_coefficient: 0.8 },
        solve: GasNozzleStepSolve { nozzle: energy, transfer: GasMassTransferSolve {
            donor: energy, receiver: energy, mass_tolerance_kg: 1e-14, energy_tolerance_j: 1e-8 }, maximum_donor_fraction: 0.05 } }]
}
#[test]
fn no_links_reduce_exactly_and_open_flow_matches_the_original_staged_owners() {
    fixture([0.2,0.2],150.0,|assembly,thermal,circuit|with_cx(|cx,_| {
        let prepared = assembly.prepare(cx,&thermal).unwrap();
        let state = MaterialExchangeState::from_parts(thermal,circuit,gases(0.004));
        let drives = drives(state.chambers());
        let (plain, _) = state.step_candidate(&prepared,cx,0.00005,&[],&drives,&[],&gas_fixture::patches(),controls(),&[0.0;2]).unwrap();
        let (idle, _) = state.flow_step_candidate(&prepared,cx,0.00005,&[],&drives,&[],&gas_fixture::patches(),&[],flow_controls(),controls(),&[0.0;2]).unwrap();
        assert_eq!(plain.chambers(),idle.chambers());assert_eq!(plain.circuit(),idle.circuit());assert_eq!(plain.thermal().cells(),idle.thermal().cells());
        let parcels = state.chambers().iter().map(|(&id,g)| (id,g.gas().clone())).collect();
        let (transferred, _) = network::step_candidate(&parcels,&links(),0.00005,flow_controls(),||false).unwrap().into_parts();
        let staged = MaterialExchangeState::from_parts(state.thermal().clone(),state.circuit().clone(),
            transferred.into_iter().map(|(id,g)| (id,GasPiston::new(g,1.0,0.01,0.0).unwrap())).collect());
        let (expected, _) = staged.step_candidate(&prepared,cx,0.00005,&[],&drives,&[],&gas_fixture::patches(),controls(),&[0.0;2]).unwrap();
        let (actual, report) = state.flow_step_candidate(&prepared,cx,0.00005,&[],&drives,&[],&gas_fixture::patches(),&links(),flow_controls(),controls(),&[0.0;2]).unwrap();
        assert_eq!(actual.chambers(),expected.chambers());assert_eq!(actual.circuit(),expected.circuit());assert_eq!(actual.thermal().cells(),expected.thermal().cells());
        assert!(report.flow.links[0].transfer.transferred_mass_kg>0.0);
        for gas in actual.chambers().values() {
            assert_eq!(gas.last_window_duration_s(),Some(0.00005));
            assert_eq!(gas.last_window_transfers().unwrap().len(),2);
            assert!(matches!(gas.audit_last_window().unwrap(),ModuleWindowAudit::Unavailable {..}));
        }
        assert_eq!(state.time_s(),0.0);
    }));
}
#[test]
fn arriving_species_data_changes_heat_exchange_electrical_history_and_motion() {
    let run = |molar_mass| fixture([0.8,0.8],150.0,|assembly,thermal,circuit|with_cx(|cx,_| {
        let prepared = assembly.prepare(cx,&thermal).unwrap();
        let mut state = MaterialExchangeState::from_parts(thermal,circuit,gases(molar_mass));
        let drives = drives(state.chambers()); let initial = state.chambers().clone();
        let energy0 = state.thermal().energy_j()+50.0+initial.values().map(|g|g.gas().mass_kg()*g.gas().state().specific_internal_energy_j_kg()).sum::<f64>();
        let mass0: f64 = initial.values().map(|g|g.gas().mass_kg()).sum();let mut work = 0.0;
        for _ in 0..40 {
            let target=state.time_s()+0.00005;
            let report=state.advance_with_flow(&prepared,cx,target,&[],&drives,&[],&gas_fixture::patches(),&links(),flow_controls(),controls(),&[0.0;2]).unwrap();
            assert!(report.total_energy_residual_j.abs()<=controls().total_energy_tolerance_j);
            assert!(report.total_mass_residual_kg.abs()<=flow_controls().mass_tolerance_kg);
            for c in &report.exchange.report.source.chambers {work+=c.boundary.piston.exterior_pressure_work_j;}
            let storage:f64=state.chambers().values().map(|g|g.gas().mass_kg()*g.gas().state().specific_internal_energy_j_kg()+0.5*g.mass_kg()*g.velocity_m_s().powi(2)).sum();
            assert!((state.thermal().energy_j()+report.exchange.report.electrical.stored_after_j+storage+work-energy0).abs()<4e-5);
        }
        assert!((state.chambers().values().map(|g|g.gas().mass_kg()).sum::<f64>()-mass0).abs()<1e-12);
        assert_eq!(state.chambers()[&10].gas().model().thermodynamics().components().len(),2);
        (state.thermal().temperature_k(CellId(30)).unwrap(),state.circuit().charges_c()[0],
            state.chambers()[&10].position_m(),state.chambers()[&10].gas().state().state().sound_speed)
    }));
    let light=run(0.004);let heavy=run(0.044);
    assert!((light.0-heavy.0).abs()>1e-8);assert!((light.1-heavy.1).abs()>1e-12);
    assert!((light.2-heavy.2).abs()>1e-10);assert!((light.3-heavy.3).abs()>1e-3);
}
#[test]
fn reactions_see_the_transferred_composition_exactly_once() {
    fixture([0.2,0.2],150.0,|assembly,thermal,circuit|with_cx(|cx,_| {
        let law=reaction(50.0);
        let product=law.catalog().with_composition(Composition::new(CompositionBasis::MoleFraction,vec![0.0,1.0]).unwrap()).unwrap();
        let prepared=assembly.prepare(cx,&thermal).unwrap();
        let state=MaterialExchangeState::from_parts(thermal,circuit,BTreeMap::from([(10,piston(&law)),(20,parcel_piston(product,180000.0))]));
        let drives=drives(state.chambers());
        let parcels=state.chambers().iter().map(|(&id,g)|(id,g.gas().clone())).collect();
        let staged=network::step_candidate(&parcels,&links(),0.00005,flow_controls(),||false).unwrap();
        let expected=law.parcel_candidate(&staged.parcels()[&10],0.00005,chemical(),||false).unwrap();
        let (_,result)=state.flow_step_candidate(&prepared,cx,0.00005,&[],&drives,
            &[ChamberReaction {chamber:10,reaction:&law,controls:chemical()}],&gas_fixture::patches(),&links(),flow_controls(),controls(),&[0.0;2]).unwrap();
        assert_eq!(result.exchange.report.source.reactions[0].reaction,*expected.step());
        assert!(result.total_energy_residual_j.abs()<=controls().total_energy_tolerance_j);
        assert_eq!(state.time_s(),0.0);
    }));
}
#[test]
fn failure_after_transport_preserves_every_live_owner_and_retry_advances_once() {
    fixture([0.2,0.2],150.0,|assembly,thermal,circuit|with_cx(|cx,gate| {
        let prepared=assembly.prepare(cx,&thermal).unwrap();
        let mut state=MaterialExchangeState::from_parts(thermal,circuit,gases(0.004));let drives=drives(state.chambers());
        let cells=state.thermal().cells().clone();let circuit=state.circuit().clone();let gases=state.chambers().clone();
        let mut limited=controls();limited.network.maximum_circuit_trials=1;
        assert!(state.advance_with_flow(&prepared,cx,0.00005,&[],&drives,&[],&gas_fixture::patches(),&links(),flow_controls(),limited,&[0.0;2]).is_err());
        assert_eq!(state.thermal().cells(),&cells);assert_eq!(state.circuit(),&circuit);assert_eq!(state.chambers(),&gases);assert_eq!(state.time_s(),0.0);
        let (preview,_)=state.flow_step_candidate(&prepared,cx,0.00005,&[],&drives,&[],&gas_fixture::patches(),&links(),flow_controls(),controls(),&[0.0;2]).unwrap();
        state.advance_with_flow(&prepared,cx,0.00005,&[],&drives,&[],&gas_fixture::patches(),&links(),flow_controls(),controls(),&[0.0;2]).unwrap();
        assert_eq!(state.chambers(),preview.chambers());assert_eq!(state.circuit(),preview.circuit());
        let saved=state.chambers().clone();gate.request();
        assert!(matches!(state.advance_with_flow(&prepared,cx,0.0001,&[],&drives,&[],&gas_fixture::patches(),&links(),flow_controls(),controls(),&[0.0;2]),Err(MaterialFlowError::Cancelled)));
        assert_eq!(state.chambers(),&saved);assert_eq!(state.time_s(),0.00005);
        for (&id,g) in state.chambers() {assert_eq!(g.last_window_transfers(),saved[&id].last_window_transfers());}
    }));
}
