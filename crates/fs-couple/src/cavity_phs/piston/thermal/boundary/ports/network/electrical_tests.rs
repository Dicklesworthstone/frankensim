//! Synthetic source cards through the original circuit, thermal and gas solvers.
use super::*;
use crate::cavity_phs::piston::tests::{controls as piston_controls, drive, gas};
use crate::cavity_phs::piston::thermal::boundary::network::material::{GasConductionGeometry,resolve_contact};
use crate::circuit_enthalpy::contacts::material::*;
use crate::circuit_enthalpy::material::*;
use fs_evidence::ValidityDomain;
use fs_material::conductor::{ElectricalPropertyBasis,ELECTRICAL_RESISTIVITY_DIMS};
use fs_material::electrical_interface::{ELECTRICAL_CONTACT_RESISTANCE_PROPERTY,ELECTRICAL_CONTACT_RESISTANCE_DIMS};
use fs_material::phase::{EnthalpyPhaseKnot,EquilibriumEnthalpyPhaseCurve,SolidLiquidPhase};
use fs_material::state_point::{MaterialPropertySelection,ResolvedMaterialStatePoint,THERMAL_CONDUCTIVITY_DIMS};
use fs_material::thermal::{source::ThermalMaterialPoint,interface::ThermalInterfacePoint};
use fs_matdb::*;
use fs_qty::Temperature;
use fs_thermal::enthalpy::*;

fn fixture<R>(contact_ohm:f64, run:impl FnOnce(&ElectricalContactNetwork<'_>,DepositionSimulation<'_>,MaterialCircuit)->R)->R {
    let cards=[2.0,4.0].map(|r| {
        let mut claims=ClaimSet::new();
        for (name,dims,value) in [("electrical_resistivity",ELECTRICAL_RESISTIVITY_DIMS,r*0.0025),
            ("thermal_conductivity",THERMAL_CONDUCTIVITY_DIMS,0.1*r)] {
            claims.insert_claim(PropertyClaim {key:PropertyKey::new(name,dims),
                value:PropertyValue::Scalar {value,dims},
                validity:ValidityDomain::unconstrained().with("T",200.0,800.0).with("lot",7.0,7.0),
                uncertainty:UncertaintyModel::Unstated,interpolation:InterpolationPolicy::ConstantWithinValidity,
                observations:vec![],provenance:Provenance {source:"synthetic contact/chamber fixture".into(),license:"CC0-1.0".into(),artifact:None}}).unwrap();
        }
        MaterialCard::assemble(MaterialStateId {chemistry:format!("fixture-{r}"),phase:"solid".into(),process:"authored".into(),revision:0},claims,vec![]).unwrap()
    });
    let curves=cards.each_ref().map(|card|EquilibriumEnthalpyPhaseCurve::try_single_phase(card.content_hash(),SolidLiquidPhase::Solid,
        [(0.0,200.0),(300000.0,800.0)].map(|(h,t)|EnthalpyPhaseKnot {specific_enthalpy_j_kg:h,temperature_k:t,
            liquid_mass_fraction:0.0,bulk_density_kg_m3:4000.0}).to_vec()).unwrap());
    let ids=[CellId(30),CellId(10)];
    let mut events=ids.iter().enumerate().map(|(i,&id)|TimedEvent {time_s:0.0,event:Event::Birth(DepositedCell {
        id,segment_id:i as u64,body_id:i as u64,material_id:MaterialId(i as u64),mass_kg:0.001,specific_h_j_kg:100000.0,
        geometry:CellGeometry {start_m:[i as f64,0.0,0.0],end_m:[i as f64+0.01,0.0,0.0],volume_m3:2.5e-7,surface_area_m2:0.04},
        boundary:Boundary {ambient_k:400.0,htc_w_m2_k:0.0,reservoir:None,source_w:0.0}})}).collect::<Vec<_>>();
    events.push(TimedEvent {time_s:0.0,event:Event::Contact(ThermalLink {id:7,a:ids[0],b:ids[1],area_m2:0.002,conductance_w_k:999.0})});
    let thermal=DepositionSimulation::new(events,&[(MaterialId(0),&curves[0] as &dyn ThermalLaw),(MaterialId(1),&curves[1] as &dyn ThermalLaw)],
        DepositionBudget {max_cells:2,max_links:1,max_events:3,max_substeps:1000}).unwrap()
        .advance_to(0.0,StepConfig {max_dt_s:0.001,energy_tolerance_j:1e-7},||Ok(())).unwrap().0;
    let point=QueryPoint::new().with("T",400.0).unwrap().with("lot",7.0).unwrap();
    let key=PropertyKey::new("thermal_conductivity",THERMAL_CONDUCTIVITY_DIMS);
    let cells=[0,1].map(|i|MaterialCell {cell:ids[i],material_id:MaterialId(i as u64),curve:&curves[i],
        point:ThermalMaterialPoint {material:&cards[i],context:&point,temperature_axis:"T",liquid_fraction_axis:None},
        conductivity:Some((&key,MaterialPropertySelection::SingleClaimOnly))});
    let resistors=[0,1].map(|i|CellResistor {destination:circuit_enthalpy::ResistorCell {branch:i+1,cell:ids[i]},
        selection:MaterialPropertySelection::SingleClaimOnly,length_m:0.01,area_m2:0.000025});
    let links=[ConductionLink {link:7,distance_a_m:0.001,distance_b_m:0.001}];
    let bulk=MaterialNetwork {cells:&cells,resistors:&resistors,links:&links};
    let mut claims=ClaimSet::new();
    claims.insert_claim(PropertyClaim {key:PropertyKey::new(ELECTRICAL_CONTACT_RESISTANCE_PROPERTY,ELECTRICAL_CONTACT_RESISTANCE_DIMS),
        value:PropertyValue::Curve {abscissa:"mean".into(),abscissa_dims:Temperature::DIMS,
            knots:vec![(200.0,0.8*contact_ohm*0.002),(800.0,1.4*contact_ohm*0.002)],dims:ELECTRICAL_CONTACT_RESISTANCE_DIMS},
        validity:ValidityDomain::unconstrained().with("TA",200.0,800.0).with("TB",200.0,800.0).with("mean",200.0,800.0).with("lot",7.0,7.0),
        uncertainty:UncertaintyModel::Unstated,interpolation:InterpolationPolicy::LinearInside,observations:vec![],
        provenance:Provenance {source:"synthetic electrical interface".into(),license:"CC0-1.0".into(),artifact:None}}).unwrap();
    let interface=InterfaceSystemCard::assemble(SurfaceSpec {material:cards[0].id().clone(),texture_frame:"A".into()},
        SurfaceSpec {material:cards[1].id().clone(),texture_frame:"B".into()},
        SystemContext {medium:"dry".into(),third_body:None,environment:"authored".into(),history:"initial".into()},claims,vec![]).unwrap();
    let context=QueryPoint::new().with("TA",400.0).unwrap().with("TB",400.0).unwrap().with("mean",400.0).unwrap().with("lot",7.0).unwrap();
    let bindings=[MaterialElectricalContact {link:7,branch:1,point:ThermalInterfacePoint {
        interface:&interface,materials:[&cards[0],&cards[1]],context:&context,temperature_axes:["TA","TB"],
        liquid_fraction_axes:[None,None],mean_temperature_axis:Some("mean")},
        selection:MaterialPropertySelection::SingleClaimOnly,heat_fraction_to_first:0.25}];
    let bases=[(1,ElectricalPropertyBasis::Resistivity),(2,ElectricalPropertyBasis::Resistivity)];
    let network=ElectricalContactNetwork {bulk:&bulk,contacts:&bindings,electrical_bases:&bases};
    let initial=network.resolve(&thermal,||false).unwrap();
    let circuit=MaterialCircuit::new_with_contacts(fs_phs::circuit::CircuitGraph {node_count:3,branches:vec![
        (1,0,fs_phs::circuit::Branch::Capacitor {farads:0.01}),(1,2,fs_phs::circuit::Branch::Resistor {ohms:999.0}),
        (2,0,fs_phs::circuit::Branch::Resistor {ohms:999.0})],transformers:vec![]},
        initial.bulk.resistors,initial.contacts,vec![],vec![1.0],&[]).unwrap();
    run(&network,thermal,circuit)
}
fn solve()->NetworkSolve {
    NetworkSolve {network:CircuitEnthalpySolve {maximum_cells:2,maximum_circuit_trials:300,maximum_newton_iterations:32,
        difference_step_j:1e-4,cell_energy_tolerance_j:1e-7,thermal_energy_tolerance_j:1e-7,total_energy_tolerance_j:1e-7,
        circuit:crate::material_circuit::MaterialCircuitSolve {energy_tolerance_j:1e-7,power_tolerance_w:1e-7},
        newton:fs_solver::NewtonKrylovConfig {absolute_tolerance:1e-8,relative_tolerance:1e-14,forcing_minimum:1e-4,linear_restart:4,
            ..fs_solver::NewtonKrylovConfig::default()}},maximum_chambers:2,maximum_contacts:4,maximum_boundary_solves:600,total_energy_tolerance_j:1e-6}
}
fn chamber(id:u64,piston:&GasPiston,dt:f64)->Chamber<'_> {
    Chamber {id,piston,input:PistonInput {duration_s:dt,..drive(0.0)},piston_solve:piston_controls(),
        exchange:HeatExchangeSolve {heat_bracket_j:[-1.0,1.0],source_tolerance_j:1e-9,energy_tolerance_j:1e-7,maximum_iterations:128}}
}
fn routes()->[ContactHeatCells;1] { [ContactHeatCells {id:7,cells:[CellId(30),CellId(10)],heat_fraction_to_first:0.25}] }
fn patches()->[Contact;2] { [Contact {id:101,chamber:10,cell:CellId(30)},Contact {id:202,chamber:20,cell:CellId(10)}] }
fn geometry()->GasConductionGeometry { GasConductionGeometry {area_m2:0.008,gas_distance_m:0.001,cell_distance_m:0.02} }
fn pistons()->[GasPiston;2] { [0.75,0.25].map(|x|GasPiston::new(gas(x,0.001),1.0,0.01,0.0).unwrap()) }
type Candidate<'m> = NetworkCandidate<'m,ContactThermalSample<MaterialNetworkSample>,ResolvedMaterialStatePoint>;
type Error = CircuitEnthalpyError<NetworkError<ElectricalContactNetworkError,MaterialNetworkError>>;
fn candidate<'m>(network:&ElectricalContactNetwork<'_>,thermal:&DepositionSimulation<'m>,circuit:&MaterialCircuit,
    gases:&[GasPiston;2],patches:&[Contact],cancelled:impl FnMut()->bool)->Result<Candidate<'m>,Error> {
    let target=thermal.time_s()+0.00005;
    let dt=target-thermal.time_s();
    let chambers=[chamber(10,&gases[0],dt),chamber(20,&gases[1],dt)];
    let destinations=network.bulk.resistors.iter().map(|r|r.destination).collect::<Vec<_>>();
    step_with_contacts_candidate(thermal,target,CircuitEnthalpyInput {circuit,ports:&[],destinations:&destinations},
        &routes(),&chambers,patches,solve(),&[0.0;2],|t|network.resolve(t,||false),
        |p,t,s,g|resolve_contact(network.bulk,p.cell,t,s,g,geometry()),cancelled)
}

#[test]
fn isolated_gases_retain_original_electrical_contact_and_thermal_endpoints() {
    fixture(0.2,|network,thermal,circuit| {
        let destinations=network.bulk.resistors.iter().map(|r|r.destination).collect::<Vec<_>>();
        let direct=electrical_contacts::step_candidate(&thermal,0.00005,
            CircuitEnthalpyInput {circuit:&circuit,ports:&[],destinations:&destinations},&routes(),solve().network,
            &[0.0;2],|t|network.resolve(t,||false),||false).unwrap();
        let gases=pistons();
        let combined=candidate(network,&thermal,&circuit,&gases,&[],||false).unwrap();
        assert_eq!(combined.thermal().cells(),direct.thermal().cells());
        assert_eq!(combined.circuit(),direct.circuit());
        assert_eq!(combined.material().destinations,routes().to_vec());
        for (i,c) in combined.chambers().iter().enumerate() {
            let mut expected=gases[i].clone();let input=chamber(c.id,&gases[i],0.00005);
            expected.step(input.input,input.piston_solve,||false).unwrap();
            assert_eq!(c.boundary.piston(),&expected);
        }
        let (_,_,report,step)=combined.into_parts();
        assert_eq!(report.electrical.contacts,direct.report().electrical.contacts);
        assert!(report.electrical.contacts[0].heat_j>0.0);
        assert!(step.external_network_heat_j.abs()<1e-7);
    });
}

#[test]
fn interface_data_changes_both_material_bodies_and_both_gas_motions() {
    let run=|r|fixture(r,|network,mut thermal,mut circuit| {
        let mut gases=pistons();let old=gases.clone();let energy0=thermal.energy_j()+50.0;let mut work=0.0;
        for _ in 0..80 {
            let result=candidate(network,&thermal,&circuit,&gases,&patches(),||false).unwrap();
            assert!(result.step().total_energy_residual_j.abs()<=solve().total_energy_tolerance_j);
            assert!(result.step().external_network_heat_j.abs()<1e-7);
            for c in result.chambers() {for p in c.boundary.source() {
                let cell=result.material().material.cells.iter().find(|s|s.cell==p.source.0.cell).unwrap();
                assert_eq!(cell.conductivity.as_ref(),Some(&p.source.1));
            }}
            let (t,c,report,_)=result.into_parts();
            assert!(report.electrical.contacts[0].heat_j>0.0);
            for (i,g) in report.source.chambers.into_iter().enumerate() {
                work+=g.boundary.step().piston.exterior_pressure_work_j;
                gases[i]=g.boundary.into_parts().0;
            }
            let gas_energy:f64=gases.iter().zip(&old).map(|(g,o)|g.gas().mass_kg()
                *(g.gas().state().specific_internal_energy_j_kg()-o.gas().state().specific_internal_energy_j_kg())
                +0.5*g.mass_kg()*g.velocity_m_s().powi(2)).sum();
            assert!((t.energy_j()+report.electrical.stored_after_j+gas_energy+work-energy0).abs()<2e-5);
            thermal=t;circuit=c;
        }
        (thermal.temperature_k(CellId(30)).unwrap(),thermal.temperature_k(CellId(10)).unwrap(),circuit.charges_c()[0],gases)
    });
    let a=run(0.2);let b=run(2.0);
    assert!((a.0-b.0).abs()>1e-4 && (a.1-b.1).abs()>1e-4 && (a.2-b.2).abs()>1e-6);
    for i in 0..2 {
        assert!((a.3[i].position_m()-b.3[i].position_m()).abs()>1e-12);
        assert!((a.3[i].gas().state().state().sound_speed-b.3[i].gas().state().state().sound_speed).abs()>1e-7);
    }
}

#[test]
fn late_cancellation_and_missing_routes_publish_no_partial_history() {
    fixture(0.2,|network,thermal,circuit| {
        let gases=pistons();let before=gases.clone();let cells=thermal.cells().clone();let original=circuit.clone();
        let mut polls=0;
        candidate(network,&thermal,&circuit,&gases,&patches(),||{polls+=1;false}).unwrap();
        let mut n=0;
        assert!(matches!(candidate(network,&thermal,&circuit,&gases,&patches(),||{n+=1;n==polls}),Err(CircuitEnthalpyError::Cancelled)));
        let destinations=network.bulk.resistors.iter().map(|r|r.destination).collect::<Vec<_>>();
        let mut calls=0;
        let result=step_with_contacts_candidate(&thermal,0.00005,CircuitEnthalpyInput {circuit:&circuit,ports:&[],destinations:&destinations},
            &[],&[chamber(10,&gases[0],0.00005)],&patches()[..1],solve(),&[0.0;2],|t|network.resolve(t,||false),
            |p,t,s,g|resolve_contact(network.bulk,p.cell,t,s,g,geometry()),||{calls+=1;false});
        assert!(matches!(result,Err(CircuitEnthalpyError::Invalid(_))));
        assert!(calls<=2); // Before nonlinear/material/physical trial work.
        assert_eq!(gases,before);assert_eq!(thermal.cells(),&cells);assert_eq!(circuit,original);
        assert!(candidate(network,&thermal,&circuit,&gases,&patches(),||false).is_ok());
    });
}

#[test]
fn contact_free_entry_point_and_empty_contact_route_are_identical() {
    fixture(0.2,|network,thermal,circuit| {
        let circuit=circuit.rebind_contacts_candidate(vec![],&[]).unwrap();
        let destinations=network.bulk.resistors.iter().map(|r|r.destination).collect::<Vec<_>>();
        let input=CircuitEnthalpyInput {circuit:&circuit,ports:&[],destinations:&destinations};
        let gases=pistons();let chambers=[chamber(10,&gases[0],0.00005),chamber(20,&gases[1],0.00005)];
        let original=step_candidate(&thermal,0.00005,input,&chambers,&patches(),solve(),&[0.0;2],
            |t|network.bulk.resolve(t,||false),|p,t,s,g|resolve_contact(network.bulk,p.cell,t,s,g,geometry()),||false).unwrap();
        let contacts=step_with_contacts_candidate(&thermal,0.00005,input,&[],&chambers,&patches(),solve(),&[0.0;2],
            |t|network.bulk.resolve(t,||false).map(|bulk|ContactThermalCoefficients {bulk,contacts:vec![]}),
            |p,t,s,g|resolve_contact(network.bulk,p.cell,t,s,g,geometry()),||false).unwrap();
        assert_eq!(contacts.thermal().cells(),original.thermal().cells());
        assert_eq!(contacts.circuit(),original.circuit());assert_eq!(contacts.step(),original.step());
        assert_eq!(&contacts.material().material,original.material());
        assert!(contacts.material().destinations.is_empty());
        for (a,b) in contacts.chambers().iter().zip(original.chambers()) {
            assert_eq!(a.boundary.piston(),b.boundary.piston());
            assert_eq!(a.boundary.source(),b.boundary.source());
        }
    });
}
