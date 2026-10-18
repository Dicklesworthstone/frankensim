//! Contact-only data changes through the original circuit and two thermal bodies.
#[path="circuit_enthalpy/common.rs"] mod common;
#[path="electrical_contacts/source.rs"] mod source;
use common::*;
use fs_couple::circuit_enthalpy::*;
use fs_couple::circuit_enthalpy::material::*;
use fs_couple::circuit_enthalpy::contacts::{self,ContactHeatCells};
use fs_couple::circuit_enthalpy::contacts::material::*;
use fs_material::conductor::ElectricalPropertyBasis;
use fs_material::phase::*;
use fs_material::state_point::{MaterialPropertySelection,THERMAL_CONDUCTIVITY_DIMS};
use fs_material::thermal::{interface::ThermalInterfacePoint,source::ThermalMaterialPoint};
use fs_matdb::*;
use fs_thermal::enthalpy::*;

fn cells<'a>(cards:&'a [MaterialCard;2],curves:[&'a EquilibriumEnthalpyPhaseCurve;2],point:&'a QueryPoint,key:&'a PropertyKey)
    ->[MaterialCell<'a>;2] {
    [0,1].map(|i|MaterialCell {cell:DESTINATIONS[i].cell,material_id:MaterialId(i as u64),curve:curves[i],
        point:ThermalMaterialPoint {material:&cards[i],context:point,temperature_axis:"T",liquid_fraction_axis:None},
        conductivity:Some((key,MaterialPropertySelection::SingleClaimOnly))})
}
fn resistors()->[CellResistor;2] {[0,1].map(|i|CellResistor {destination:DESTINATIONS[i],
    selection:MaterialPropertySelection::SingleClaimOnly,length_m:1.0,area_m2:1.0})}
fn links()->[ConductionLink;1] {[ConductionLink {link:7,distance_a_m:0.5,distance_b_m:0.5}]}
fn bases()->[(usize,ElectricalPropertyBasis);2] {[(1,ElectricalPropertyBasis::Resistivity),(2,ElectricalPropertyBasis::Resistivity)]}
fn contact<'a>(cards:&'a [MaterialCard;2],card:&'a InterfaceSystemCard,point:&'a QueryPoint,fraction:f64)->MaterialElectricalContact<'a> {
    MaterialElectricalContact {link:7,branch:1,point:ThermalInterfacePoint {interface:card,materials:[&cards[0],&cards[1]],
        context:point,temperature_axes:["Ta","Tb"],liquid_fraction_axes:[Some("melt_a"),Some("melt_b")],mean_temperature_axis:None},
        selection:MaterialPropertySelection::SingleClaimOnly,heat_fraction_to_first:fraction}
}
fn routes(fraction:f64)->[ContactHeatCells;1] {[ContactHeatCells {id:7,cells:[CellId(30),CellId(10)],heat_fraction_to_first:fraction}]}

#[test]
fn contact_data_drives_both_body_temperatures_and_finite_electrical_history() {
    let run=|contact_scale| {
        let (a,ca)=material(1.0,1.0,0.005);let (b,cb)=material(3.0,2.0,0.0);
        let cards=[a,b];let curves=[ca,cb];let point=point(300.0);let interface_point=source::point();
        let key=PropertyKey::new("thermal_conductivity",THERMAL_CONDUCTIVITY_DIMS);
        let cells=cells(&cards,[curves[0].curve(),curves[1].curve()],&point,&key);let resistors=resistors();let links=links();let bases=bases();
        let bulk=MaterialNetwork {cells:&cells,resistors:&resistors,links:&links};
        let card=source::card(&cards[0],&cards[1],contact_scale,0.01);
        let bindings=[contact(&cards,&card,&interface_point,0.25)];
        let network=ElectricalContactNetwork {bulk:&bulk,contacts:&bindings,electrical_bases:&bases};
        let mut thermal=initial(&curves,999.0,None);
        let mut electrical=circuit(&thermal,&cards,&curves).rebind_contacts_candidate(network.resolve(&thermal,||false).unwrap().contacts,&[]).unwrap();
        let total0=thermal.energy_j()+50.0;
        for n in 1..=40 {
            let target=f64::from(n)*0.05;let dt=target-thermal.time_s();
            let candidate=network.step_candidate(&thermal,target,&electrical,&[],config(),&[0.0;2],||false).unwrap();
            let report=candidate.report();let source=&report.source.material;
            let a=source.cells.iter().find(|c|c.cell==CellId(30)).unwrap().phase;
            let b=source.cells.iter().find(|c|c.cell==CellId(10)).unwrap().phase;
            let contact_heat=report.electrical.contacts[0].heat_j;
            near(report.electrical.contacts[0].contact.resistance_ohm(),contact_scale*(1.0+0.01*(a.temperature_k()-300.0)));
            let exchanged=dt*source.links[0].conductance_w_k*(b.temperature_k()-a.temperature_k());
            let da=candidate.thermal().cells()[&CellId(30)].specific_h_j_kg-thermal.cells()[&CellId(30)].specific_h_j_kg;
            let db=candidate.thermal().cells()[&CellId(10)].specific_h_j_kg-thermal.cells()[&CellId(10)].specific_h_j_kg;
            near(da,report.electrical.resistors[0].heat_j+0.25*contact_heat+exchanged);
            near(db,report.electrical.resistors[1].heat_j+0.75*contact_heat-exchanged);
            near(candidate.thermal().energy_j()+report.electrical.stored_after_j,total0);
            near(report.deposited_joule_j,report.electrical.heat_j);
            near(report.non_electrical_heat_j,0.0);
            assert_eq!(report.source.destinations,routes(0.25));
            let (t,c,_)=candidate.into_parts();thermal=t;electrical=c;
        }
        (thermal.temperature_k(CellId(30)).unwrap(),thermal.temperature_k(CellId(10)).unwrap(),electrical.charges_c()[0])
    };
    let lower=run(0.2);let higher=run(2.0);
    assert!(higher.2>lower.2);
    assert!((higher.0-lower.0).abs()>0.1 && (higher.1-lower.1).abs()>0.1);
}

#[test]
fn explicit_heat_partition_matches_analytic_isolated_two_body_limit() {
    for fraction in [0.0,0.25,1.0] {
        let (a,ca)=material_with_conductivity(1.0,1.0,0.0,0.0);let (b,cb)=material_with_conductivity(3.0,2.0,0.0,0.0);
        let cards=[a,b];let curves=[ca,cb];let point=point(300.0);let interface_point=source::point();
        let key=PropertyKey::new("thermal_conductivity",THERMAL_CONDUCTIVITY_DIMS);
        let cells=cells(&cards,[curves[0].curve(),curves[1].curve()],&point,&key);let resistors=resistors();let links=links();let bases=bases();
        let bulk=MaterialNetwork {cells:&cells,resistors:&resistors,links:&links};let card=source::card(&cards[0],&cards[1],2.0,0.0);
        let bindings=[contact(&cards,&card,&interface_point,fraction)];
        let network=ElectricalContactNetwork {bulk:&bulk,contacts:&bindings,electrical_bases:&bases};
        let thermal=initial(&curves,999.0,None);
        let electrical=circuit(&thermal,&cards,&curves).rebind_contacts_candidate(network.resolve(&thermal,||false).unwrap().contacts,&[]).unwrap();
        let result=network.step_candidate(&thermal,0.25,&electrical,&[],config(),&[0.0;2],||false).unwrap();
        let alpha=0.25/12.0;let q1=10.0*(1.0-alpha)/(1.0+alpha);let current=0.5*(10.0+q1)/6.0;
        let unit_heat=0.25*current*current;
        near(result.circuit().charges_c()[0],q1);
        near(result.thermal().temperature_k(CellId(30)).unwrap(),300.0+unit_heat*(1.0+2.0*fraction));
        near(result.thermal().temperature_k(CellId(10)).unwrap(),300.0+unit_heat*(3.0+2.0*(1.0-fraction))/2.0);
    }
}

#[test]
fn cancellation_missing_routes_and_source_switching_publish_no_partial_history() {
    let (a,ca)=material(1.0,1.0,0.005);let (b,cb)=material(3.0,2.0,0.0);
    let cards=[a,b];let curves=[ca,cb];let point=point(300.0);let interface_point=source::point();
    let key=PropertyKey::new("thermal_conductivity",THERMAL_CONDUCTIVITY_DIMS);
    let cells=cells(&cards,[curves[0].curve(),curves[1].curve()],&point,&key);let resistors=resistors();let links=links();let bases=bases();
    let bulk=MaterialNetwork {cells:&cells,resistors:&resistors,links:&links};let card=source::card(&cards[0],&cards[1],2.0,0.01);
    let mut bindings=[contact(&cards,&card,&interface_point,0.25)];
    let thermal=initial(&curves,999.0,None);
    let electrical=circuit(&thermal,&cards,&curves).rebind_contacts_candidate(
        ElectricalContactNetwork {bulk:&bulk,contacts:&bindings,electrical_bases:&bases}.resolve(&thermal,||false).unwrap().contacts,&[]).unwrap();
    let saved_cells=thermal.cells().clone();let saved_circuit=electrical.clone();
    for fraction in [-0.1,1.1,f64::NAN] {
        bindings[0].heat_fraction_to_first=fraction;
        assert!(ElectricalContactNetwork {bulk:&bulk,contacts:&bindings,electrical_bases:&bases}
            .step_candidate(&thermal,0.1,&electrical,&[],config(),&[0.0;2],||false).is_err());
    }
    bindings[0].heat_fraction_to_first=0.25;
    let network=ElectricalContactNetwork {bulk:&bulk,contacts:&bindings,electrical_bases:&bases};
    assert!(contacts::step_candidate(&thermal,0.1,input(&electrical),&[],config(),&[0.0;2],
        |s|network.resolve(s,||false),||false).is_err());
    let other=source::card(&cards[0],&cards[1],4.0,0.01);
    assert!(contacts::step_candidate(&thermal,0.1,input(&electrical),&routes(0.25),config(),&[0.0;2],|s| {
        let mut c=network.resolve(s,||false)?;c.contacts=vec![source::contact(7,1,&other,1.0)];
        Ok::<_,ElectricalContactNetworkError>(c)
    },||false).is_err());
    let mut polls=0;
    network.step_candidate(&thermal,0.1,&electrical,&[],config(),&[0.0;2],||{polls+=1;false}).unwrap();
    let mut actual=0;
    assert!(matches!(network.step_candidate(&thermal,0.1,&electrical,&[],config(),&[0.0;2],||{actual+=1;actual==polls}),
        Err(CircuitEnthalpyError::Cancelled)));
    assert!(network.step_candidate(&thermal,0.1,&electrical,&[],CircuitEnthalpySolve {maximum_circuit_trials:1,..config()},&[0.0;2],||false).is_err());
    let driven=contacts::step_candidate(&thermal,0.1,input(&electrical),&routes(0.25),config(),&[0.0;2],|s| {
        let mut c=network.resolve(s,||false)?;c.bulk.sources_w.insert(CellId(30),2.0);Ok::<_,ElectricalContactNetworkError>(c)
    },||false).unwrap();
    near(driven.report().non_electrical_heat_j,0.2);
    near(driven.report().electrical.stored_after_j-50.0+driven.thermal().energy_j()-thermal.energy_j(),0.2);
    assert_eq!(thermal.cells(),&saved_cells);assert_eq!(electrical,saved_circuit);
}

#[test]
fn latent_fraction_changes_contact_current_without_an_invented_temperature_rise() {
    use fs_evidence::ValidityDomain;
    use fs_material::electrical_interface::{ELECTRICAL_CONTACT_RESISTANCE_PROPERTY,ELECTRICAL_CONTACT_RESISTANCE_DIMS};
    use fs_qty::Dims;
    let (a,_)=material_with_conductivity(1.0,1.0,0.0,0.0);let (b,cb)=material_with_conductivity(3.0,2.0,0.0,0.0);
    let cards=[a,b];
    let knot=|h,t,x|EnthalpyPhaseKnot {specific_enthalpy_j_kg:h,temperature_k:t,liquid_mass_fraction:x,bulk_density_kg_m3:1.0};
    let first=EquilibriumEnthalpyPhaseCurve::try_new(cards[0].content_hash(),vec![
        knot(-100.0,250.0,0.0),knot(0.0,300.0,0.0),knot(20.0,300.0,1.0),knot(1020.0,800.0,1.0)]).unwrap();
    let point=point(300.0);let interface_point=source::point();
    let key=PropertyKey::new("thermal_conductivity",THERMAL_CONDUCTIVITY_DIMS);
    let cells=cells(&cards,[&first,cb.curve()],&point,&key);let resistors=resistors();let links=links();let bases=bases();
    let bulk=MaterialNetwork {cells:&cells,resistors:&resistors,links:&links};
    let mut claims=ClaimSet::new();
    claims.insert_claim(PropertyClaim {key:PropertyKey::new(ELECTRICAL_CONTACT_RESISTANCE_PROPERTY,ELECTRICAL_CONTACT_RESISTANCE_DIMS),
        value:PropertyValue::Curve {abscissa:"melt_a".into(),abscissa_dims:Dims::NONE,
            knots:vec![(0.0,2.0),(1.0,12.0)],dims:ELECTRICAL_CONTACT_RESISTANCE_DIMS},
        validity:ValidityDomain::unconstrained().with("Ta",200.0,2000.0).with("Tb",200.0,2000.0)
            .with("melt_a",0.0,1.0).with("pressure",100000.0,100000.0),
        uncertainty:UncertaintyModel::Unstated,interpolation:InterpolationPolicy::LinearInside,observations:vec![],
        provenance:Provenance {source:"authored phase contact, not measured qualification".into(),license:"CC0-1.0".into(),artifact:None}}).unwrap();
    let card=InterfaceSystemCard::assemble(SurfaceSpec {material:cards[0].id().clone(),texture_frame:"first finish".into()},
        SurfaceSpec {material:cards[1].id().clone(),texture_frame:"second finish".into()},
        SystemContext {medium:"authored contact".into(),third_body:None,environment:"test".into(),history:"fixed".into()},claims,vec![]).unwrap();
    let bindings=[contact(&cards,&card,&interface_point,0.25)];
    let network=ElectricalContactNetwork {bulk:&bulk,contacts:&bindings,electrical_bases:&bases};
    let mut events=[0,1].map(|i|TimedEvent {time_s:0.0,event:Event::Birth(DepositedCell {
        id:DESTINATIONS[i].cell,segment_id:i as u64,body_id:i as u64,material_id:MaterialId(i as u64),mass_kg:1.0,
        specific_h_j_kg:if i==0 {5.0} else {200.0},geometry:CellGeometry {start_m:[i as f64,0.0,0.0],
            end_m:[i as f64+1.0,0.0,0.0],volume_m3:1.0,surface_area_m2:6.0},boundary:boundary()})}).to_vec();
    events.push(TimedEvent {time_s:0.0,event:Event::Contact(ThermalLink {id:7,a:CellId(30),b:CellId(10),area_m2:1.0,conductance_w_k:0.0})});
    let thermal=DepositionSimulation::new(events,&[(MaterialId(0),&first as &dyn ThermalLaw),(MaterialId(1),cb.curve() as &dyn ThermalLaw)],
        DepositionBudget {max_cells:2,max_links:1,max_events:3,max_substeps:100}).unwrap()
        .advance_to(0.0,StepConfig {max_dt_s:0.1,energy_tolerance_j:1e-7},||Ok(())).unwrap().0;
    let initial=network.resolve(&thermal,||false).unwrap();
    let electrical=fs_couple::material_circuit::MaterialCircuit::new_with_contacts(fs_phs::circuit::CircuitGraph {
        node_count:3,branches:vec![(1,0,fs_phs::circuit::Branch::Capacitor {farads:1.0}),
            (1,2,fs_phs::circuit::Branch::Resistor {ohms:0.0}),(2,0,fs_phs::circuit::Branch::Resistor {ohms:0.0})],transformers:vec![]},
        initial.bulk.resistors,initial.contacts,vec![],vec![10.0],&[]).unwrap();
    let frozen=electrical.step_with_contacts_candidate(electrical.resistors().to_vec(),electrical.contacts().to_vec(),&[],0.25,config().circuit,||false).unwrap();
    let result=network.step_candidate(&thermal,0.25,&electrical,&[],config(),&[0.0;2],||false).unwrap();
    let phase=result.report().source.material.cells.iter().find(|c|c.cell==CellId(30)).unwrap().phase;
    assert_eq!(phase.temperature_k(),300.0);
    assert!(phase.liquid_mass_fraction()>0.25 && phase.liquid_mass_fraction()<1.0);
    near(result.report().electrical.contacts[0].contact.resistance_ohm(),2.0+10.0*phase.liquid_mass_fraction());
    assert!(result.circuit().charges_c()[0]>frozen.circuit().charges_c()[0]);
    near(result.thermal().energy_j()-thermal.energy_j()+result.report().electrical.stored_after_j-50.0,0.0);
}
