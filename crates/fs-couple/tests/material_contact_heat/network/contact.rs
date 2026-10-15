//! Source-selected contact area also controls thermal links and exposed cooling.
use super::*;
use fs_conduction::interface::{AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY as PROPERTY,
    AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS as DIMS};
use fs_couple::render::plate::impact::material::MaterialContact;
use fs_couple::render::plate::impact::material::heating::network::contact::ContactNetworkLink;
use fs_couple::render::plate::impact::material::heating::feedback::resolve_contact_at_temperatures;
use fs_contact::normal_patch::source::state::PointContactMaterialCards;
use fs_contact::interface_binding::NormalContactModelSelection;
use fs_material::state_point::MaterialPropertySelection;
use fs_evidence::ValidityDomain;
use fs_matdb::{MaterialCard,InterfaceSystemCard,SelectionPolicy,PropertyClaim,PropertyKey,
    PropertyValue,Provenance,UncertaintyModel,InterpolationPolicy,SystemContext};
use fs_tribo::InputAuthority;

fn interface(a:&MaterialCard,b:&MaterialCard,resistance:Option<f64>)->InterfaceSystemCard {
    let source=crate::common::bind(a,b,Some(0.4),0.0);
    let old=source.normal_state().elastic().interface().card();let mut claims=old.claims().clone();
    if let Some(value)=resistance {
        claims.insert_claim(PropertyClaim {key:PropertyKey::new(PROPERTY,DIMS),
            value:PropertyValue::Scalar {value,dims:DIMS},validity:ValidityDomain::unconstrained().with("T",280.0,320.0),
            uncertainty:UncertaintyModel::Unstated,interpolation:InterpolationPolicy::ConstantWithinValidity,
            observations:vec![],provenance:Provenance {source:"synthetic thermal-interface coefficient".into(),
                license:"CC0-1.0".into(),artifact:None}}).unwrap();
    }
    InterfaceSystemCard::assemble(old.surface_a().clone(),old.surface_b().clone(),
        SystemContext {medium:"dry".into(),third_body:None,environment:"test".into(),history:"initial".into()},
        claims,old.models().to_vec()).unwrap()
}
fn caloric(card:&MaterialCard,rho:f64)->EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_single_phase(card.content_hash(),SolidLiquidPhase::Solid,
        vec![EnthalpyPhaseKnot {specific_enthalpy_j_kg:0.0,temperature_k:300.0,liquid_mass_fraction:0.0,bulk_density_kg_m3:rho},
            EnthalpyPhaseKnot {specific_enthalpy_j_kg:0.4,temperature_k:320.0,liquid_mass_fraction:0.0,bulk_density_kg_m3:rho}]).unwrap()
}
fn cards<'a>(a:&'a MaterialCard,b:&'a MaterialCard,i:&'a InterfaceSystemCard)->PointContactMaterialCards<'a> {
    PointContactMaterialCards {surface_a:a,surface_b:b,interface:i,
        bulk_selection:[MaterialPropertySelection::SingleClaimOnly,MaterialPropertySelection::SingleClaimOnly],
        interface_selection:MaterialPropertySelection::SingleClaimOnly,
        model_selection:NormalContactModelSelection::SingleSupported,authority:InputAuthority::SyntheticFixture}
}
fn setup<'a>(cards:&PointContactMaterialCards<'_>,curves:[&'a EquilibriumEnthalpyPhaseCurve;2])
    ->(HeatedImpactSystem<'a>,DepositionSimulation<'a>) {
    let source=cards.resolve(&query(),PointContactGeometry::SpherePlane {effective_radius_m:0.1},envelope()).unwrap();
    let (a,wa)=ImpactBody::free_mass(0.001,-2e-5,0.04).unwrap();
    let (b,wb)=ImpactBody::free_mass(0.002,0.0,0.0).unwrap();
    let contacts=[("touch",0.0),("other",1e-6)].into_iter().map(|(key,gap)|
        MaterialContact::new(key.into(),source.clone(),vec![wa,-wb],gap).unwrap()).collect();
    let mechanics=ImpactSystem::new_with_material_contacts(vec![a,b],contacts,vec![],vec![],config()).unwrap();
    let spec=[("left",9,0,0.2,1e-6),("right",3,1,0.04,0.5e-6),("support",77,0,0.0,1e-6)];
    let mut stores=Vec::new();let mut events=Vec::new();
    for (key,id,material,h,volume) in spec {
        stores.push(CaloricReservoir::new(key.into(),curves[material],curves[material].state_at_specific_enthalpy(h).unwrap(),
            0.001,volume,0.0006).unwrap());
        events.push(TimedEvent {time_s:0.0,event:Event::Birth(DepositedCell {id:CellId(id),segment_id:id,body_id:id,
            material_id:MaterialId(material as u64),mass_kg:0.001,specific_h_j_kg:h,
            geometry:CellGeometry {start_m:[id as f64,0.0,0.0],end_m:[id as f64+0.01,0.0,0.0],volume_m3:volume,surface_area_m2:0.0006},
            boundary:Boundary {ambient_k:300.0,htc_w_m2_k:0.5,reservoir:None,source_w:0.0}})});
    }
    events.push(TimedEvent {time_s:0.0,event:Event::Contact(ThermalLink {id:102,a:CellId(3),b:CellId(77),area_m2:1e-6,conductance_w_k:0.01})});
    let network=DepositionSimulation::new(events,&[(MaterialId(0),curves[0] as &dyn ThermalLaw),(MaterialId(1),curves[1] as &dyn ThermalLaw)],
        DepositionBudget {max_cells:3,max_links:3,max_events:16,max_substeps:16}).unwrap()
        .advance_to(0.0,thermal_config(5e-7),||Ok(())).unwrap().0;
    (HeatedImpactSystem::new(mechanics.prepare_analytic().unwrap(),stores,vec![route("touch"),route("other")],1e-10).unwrap(),network)
}
fn links()->Vec<ContactNetworkLink> {
    [("other",201),("touch",101)].into_iter().map(|(contact,link_id)|ContactNetworkLink {
        contact:contact.into(),link_id,cells:[CellId(9),CellId(3)],selection:SelectionPolicy::SingleClaimOnly}).collect()
}

#[test]
fn interface_data_alone_changes_network_heat_while_real_footprints_open_and_close() {
    let a=solid("A",1000.0,1e7);let b=solid("B",2000.0,2e7);let ca=caloric(&a,1000.0);let cb=caloric(&b,2000.0);
    let ia=interface(&a,&b,Some(1e-4));let ib=interface(&a,&b,Some(4e-4));let ma=cards(&a,&b,&ia);let mb=cards(&a,&b,&ib);
    let (mut fast,mut na)=setup(&ma,[&ca,&cb]);let (mut slow,mut nb)=setup(&mb,[&ca,&cb]);let gate=CancelGate::new_clock_free();
    let initial=fast.mechanics().stored_energy_j();let mut boundary=0.0;let mut active=0;let mut separated=0;
    for _ in 0..1600 {
        let fa=fast.step_with_contact_network(&mut na,&mut (),&bindings(),&links(),thermal_config(5e-7),&[0.0;2],&gate,
            |c,s|resolve_contact_at_temperatures(c,&ma,&query(),[&s[0],&s[1]],s[0].state().temperature_k()).map(Some)).unwrap();
        slow.step_with_contact_network(&mut nb,&mut (),&bindings(),&links(),thermal_config(5e-7),&[0.0;2],&gate,
            |c,s|resolve_contact_at_temperatures(c,&mb,&query(),[&s[0],&s[1]],s[0].state().temperature_k()).map(Some)).unwrap();
        assert_eq!(fast.mechanics().state(),slow.mechanics().state());
        assert_eq!(fa.contacts.iter().map(|r|r.link_id).collect::<Vec<_>>(),[101,201]);
        let mut area=0.0;
        for report in &fa.contacts {
            assert_eq!(report.interface.card_identity(),ia.content_hash());area+=report.area_m2;
            if report.area_m2==0.0 {separated+=1;assert!(!na.links().contains_key(&report.link_id));}
            else {active+=1;near(na.links()[&report.link_id].conductance_w_k,report.area_m2/1e-4,1e-15);}
        }
        near(na.exposed_area_m2(CellId(9)).unwrap(),0.0006-area,1e-18);
        assert!(na.links().contains_key(&102));boundary+=fa.network.non_contact_heat_j;
    }
    assert!(active>0 && separated>0);assert!(!na.links().contains_key(&101));
    assert!(na.temperature_k(CellId(9)).unwrap()<nb.temperature_k(CellId(9)).unwrap());
    near(fast.mechanics().stored_energy_j()+fast.reservoirs().iter().map(|r|r.storage_change_j()).sum::<f64>()-boundary,initial,1e-9);
}

#[derive(Clone,Default)]
struct Refuse {steps:usize}
impl DepositionObserver for Refuse {
    type Error=DepositionError;
    fn event(&mut self,_:&Event,_:&DepositionSimulation<'_>)->Result<(),Self::Error>{Ok(())}
    fn step(&mut self,_:f64,_:f64,_:&BTreeMap<CellId,f64>,_:&DepositionSimulation<'_>,
        _:&mut dyn FnMut()->Result<(),DepositionError>)->Result<(),Self::Error>{
        self.steps+=1;if self.steps==2{Err(DepositionError::MaterialDomain)}else{Ok(())}
    }
}
#[test]
fn late_refusal_keeps_physical_link_topology_and_sources_at_the_accepted_time() {
    let a=solid("A",1000.0,1e7);let b=solid("B",2000.0,2e7);let ca=caloric(&a,1000.0);let cb=caloric(&b,2000.0);
    let i=interface(&a,&b,Some(1e-4));let cards=cards(&a,&b,&i);let (mut s,mut n)=setup(&cards,[&ca,&cb]);let gate=CancelGate::new_clock_free();
    for _ in 0..350 {s.step_with_contact_network(&mut n,&mut (),&bindings(),&links(),thermal_config(5e-7),&[0.0;2],&gate,unchanged).unwrap();}
    assert!(n.links().contains_key(&101));let old_links=n.links().clone();let old_cells=n.cells().clone();let time=n.time_s();
    let motion=s.mechanics().state().to_vec();let stores=s.reservoirs().to_vec();let mut observer=Refuse::default();
    assert!(s.step_with_contact_network(&mut n,&mut observer,&bindings(),&links(),thermal_config(5e-7),&[0.0;2],&gate,unchanged).is_err());
    assert_eq!(observer.steps,0);assert_eq!(n.links(),&old_links);assert_eq!(n.cells(),&old_cells);assert_eq!(n.time_s(),time);
    assert_eq!(s.mechanics().state(),motion);assert_eq!(s.reservoirs(),stores);
    // A separate observer cannot silently overwrite the physical source law.
    assert!(s.step_with_contact_network(&mut n,&mut MaterialObserver::default(),&bindings(),&links(),thermal_config(5e-7),&[0.0;2],&gate,unchanged).is_err());
    assert_eq!(n.links(),&old_links);assert_eq!(s.mechanics().state(),motion);
}

#[test]
fn missing_interface_law_reversed_material_roles_and_link_aliases_refuse() {
    let a=solid("A",1000.0,1e7);let b=solid("B",2000.0,2e7);let ca=caloric(&a,1000.0);let cb=caloric(&b,2000.0);
    let absent=interface(&a,&b,None);let present=interface(&a,&b,Some(1e-4));let gate=CancelGate::new_clock_free();
    for i in [&absent,&present] {
        let cards=cards(&a,&b,i);let (mut s,mut n)=setup(&cards,[&ca,&cb]);let mut requested=links();
        if i.content_hash()==present.content_hash(){requested[0].cells.reverse();}
        assert!(s.step_with_contact_network(&mut n,&mut (),&bindings(),&requested,thermal_config(5e-7),&[0.0;2],&gate,unchanged).is_err());
        assert_eq!(s.mechanics().samples(),0);assert_eq!(n.time_s(),0.0);assert_eq!(n.links().len(),1);
    }
    let cards=cards(&a,&b,&present);let (mut s,mut n)=setup(&cards,[&ca,&cb]);let mut aliased=links();aliased[0].link_id=102;
    assert!(s.step_with_contact_network(&mut n,&mut (),&bindings(),&aliased,thermal_config(5e-7),&[0.0;2],&gate,unchanged).is_err());
    assert_eq!(n.links().len(),1);assert_eq!(s.mechanics().samples(),0);
}
