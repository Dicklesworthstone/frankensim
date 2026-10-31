//! Independent circuit/interface heat checks on continuing material memories.
use super::*;
use crate::circuit_enthalpy::material::ConductionLink;
use crate::circuit_enthalpy::contacts::material::{ElectricalContactNetwork, MaterialElectricalContact};
use crate::material_circuit::contacts::SeriesElectricalContact;
use fs_material::electrical_interface::{ELECTRICAL_CONTACT_RESISTANCE_DIMS, ELECTRICAL_CONTACT_RESISTANCE_PROPERTY};
use fs_material::thermal::interface::ThermalInterfacePoint;
use fs_material::state_point::THERMAL_CONDUCTIVITY_DIMS;
use fs_matdb::{InterfaceSystemCard, SurfaceSpec, SystemContext};

const LINK:u64=7;
const CONTACT_AREA:f64=0.01;
#[derive(Clone,Copy)]
struct Contact { resistance:f64, fraction:f64, upper_k:f64 }
fn spec()->Contact {Contact {resistance:15.0,fraction:0.35,upper_k:340.0}}
fn sources(layer:usize)->(MaterialCard,ResolvedHeatCapacityCurve) {
    let (card,_)=powered_source(layer,1.0,false);
    let mut claims=ClaimSet::new();
    for (_,claim) in card.claims().claims_ordered(){claims.insert_claim(claim.clone()).unwrap();}
    // Explicit insulating thermal law, independently of finite electrical contact.
    claims.insert_claim(PropertyClaim {key:PropertyKey::new("thermal_conductivity",THERMAL_CONDUCTIVITY_DIMS),
        value:PropertyValue::Scalar {value:0.0,dims:THERMAL_CONDUCTIVITY_DIMS},
        validity:ValidityDomain::unconstrained().with_quantity("T",QuantitySpec::dimensional(Temperature::DIMS),290.0,340.0),
        interpolation:InterpolationPolicy::ConstantWithinValidity,uncertainty:UncertaintyModel::Unstated,
        observations:vec![],provenance:Provenance {source:"authored insulation fixture".into(),license:"test-only".into(),artifact:None}}).unwrap();
    let card=MaterialCard::assemble(card.id().clone(),claims,vec![]).unwrap();
    let curve=ResolvedHeatCapacityCurve::try_new(&card,HeatCapacitySourceQuery {
        heat_capacity:&PropertyKey::new("specific_heat_capacity",SPECIFIC_HEAT_CAPACITY_DIMS),
        density:&PropertyKey::new("density",Density::DIMS),temperature_axis:"T",
        lower:&point(290.0),upper:&point(340.0),selection:SelectionPolicy::SingleClaimOnly,
    },SolidLiquidPhase::Solid,0.0,1e-6).unwrap();
    (card,curve)
}
fn interface(cards:&[MaterialCard;2],s:Contact)->InterfaceSystemCard {
    let mut claims=ClaimSet::new();let value=CONTACT_AREA*s.resistance;
    claims.insert_claim(PropertyClaim {
        key:PropertyKey::new(ELECTRICAL_CONTACT_RESISTANCE_PROPERTY,ELECTRICAL_CONTACT_RESISTANCE_DIMS),
        value:PropertyValue::Curve {abscissa:"Ta".into(),abscissa_dims:Temperature::DIMS,
            knots:vec![(290.0,0.8*value),(s.upper_k,value*(1.0+0.02*(s.upper_k-300.0)))],
            dims:ELECTRICAL_CONTACT_RESISTANCE_DIMS},
        validity:ValidityDomain::unconstrained()
            .with_quantity("Ta",QuantitySpec::dimensional(Temperature::DIMS),290.0,s.upper_k)
            .with_quantity("Tb",QuantitySpec::dimensional(Temperature::DIMS),290.0,340.0),
        interpolation:InterpolationPolicy::LinearInside,uncertainty:UncertaintyModel::Unstated,
        observations:vec![],provenance:Provenance {source:"manufactured electrical layer interface".into(),license:"test-only".into(),artifact:None},
    }).unwrap();
    InterfaceSystemCard::assemble(
        SurfaceSpec {material:cards[0].id().clone(),texture_frame:"first".into()},
        SurfaceSpec {material:cards[1].id().clone(),texture_frame:"second".into()},
        SystemContext {medium:"declared contact".into(),third_body:None,
            environment:"synthetic".into(),history:"fixed geometry".into()},claims,vec![]).unwrap()
}
fn fixture<R>(s:Contact,run:impl FnOnce(&MaterialFluidFrame<'_,FluidTransportState>,
    &MaterialNetwork<'_>,&ElectricalContactNetwork<'_>,DepositionSimulation<'_>,MaterialCircuit)->R)->R {
    let (a,ac)=sources(0);let (b,bc)=sources(1);let cards=[a,b];let curves=[ac,bc];
    let f=frame(&cards);let m=model(&f);let query=point(300.0);
    let key=PropertyKey::new("thermal_conductivity",THERMAL_CONDUCTIVITY_DIMS);
    let cells:[_;2]=core::array::from_fn(|i|MaterialCell {cell:CellId(i as u64),material_id:MaterialId(i as u64),curve:curves[i].curve(),
        point:ThermalMaterialPoint {material:&cards[i],context:&query,temperature_axis:"T",liquid_fraction_axis:None},
        conductivity:Some((&key,MaterialPropertySelection::SingleClaimOnly))});
    let mut events=Vec::new();
    for i in 0..2 {events.push(TimedEvent {time_s:0.0,event:Event::Birth(DepositedCell {
        id:CellId(i as u64),segment_id:0,body_id:0,material_id:MaterialId(i as u64),
        mass_kg:[3.0,5.0][i]*area(i),specific_h_j_kg:10.0*CP,
        geometry:CellGeometry {start_m:[0.0;3],end_m:[1.0,0.0,0.0],volume_m3:area(i),surface_area_m2:1.0},
        boundary:Boundary {ambient_k:300.0,htc_w_m2_k:0.0,reservoir:None,source_w:0.0},
    })});}
    events.push(TimedEvent {time_s:0.0,event:Event::Contact(ThermalLink {
        id:LINK,a:CellId(0),b:CellId(1),area_m2:CONTACT_AREA,conductance_w_k:999.0,
    })});
    let heat=DepositionSimulation::new(events,&[(MaterialId(0),curves[0].curve() as &dyn ThermalLaw),
        (MaterialId(1),curves[1].curve() as &dyn ThermalLaw)],
        DepositionBudget {max_cells:2,max_links:1,max_events:3,max_substeps:128}).unwrap()
        .advance_to(0.0,StepConfig {max_dt_s:1.0/256.0,energy_tolerance_j:1e-11},||Ok(())).unwrap().0;
    let links=[ConductionLink {link:LINK,distance_a_m:0.5,distance_b_m:0.5}];
    let plain=MaterialNetwork {cells:&cells,resistors:&[],links:&links};
    let resistors:[_;2]=core::array::from_fn(|i|CellResistor {destination:ResistorCell {branch:i+1,cell:CellId(i as u64)},
        selection:MaterialPropertySelection::SingleClaimOnly,length_m:1.0,area_m2:area(i)});
    let bulk=MaterialNetwork {cells:plain.cells,resistors:&resistors,links:plain.links};
    let bases=[(1,ElectricalPropertyBasis::Resistivity),(2,ElectricalPropertyBasis::Resistivity)];
    let card=interface(&cards,s);
    let q=QueryPoint::new().with_quantity("Ta",QuantitySpec::dimensional(Temperature::DIMS),300.0).unwrap()
        .with_quantity("Tb",QuantitySpec::dimensional(Temperature::DIMS),300.0).unwrap();
    let contacts=[MaterialElectricalContact {link:LINK,branch:1,
        point:ThermalInterfacePoint {interface:&card,materials:[&cards[0],&cards[1]],context:&q,
            temperature_axes:["Ta","Tb"],liquid_fraction_axes:[None,None],mean_temperature_axis:None},
        selection:MaterialPropertySelection::SingleClaimOnly,heat_fraction_to_first:s.fraction}];
    let source=ElectricalContactNetwork {bulk:&bulk,contacts:&contacts,electrical_bases:&bases};
    let initial=source.prepare(&heat,||false).unwrap().resolve(&heat,||false).unwrap();
    let circuit=MaterialCircuit::new_with_contacts(CircuitGraph {node_count:3,branches:vec![
        (1,0,Branch::Capacitor {farads:CAPACITANCE}),
        (1,2,Branch::Resistor {ohms:999.0}),(2,0,Branch::Resistor {ohms:999.0})],transformers:vec![]},
        initial.bulk.resistors,initial.contacts,vec![],vec![INITIAL_CHARGE],&[]).unwrap();
    run(&m,&plain,&source,heat,circuit)
}

#[test]
fn interface_heat_reaches_both_original_layers_and_closes_rc_and_total_energy() {
    for fraction in [0.0,0.35,1.0] {
        let s=Contact {fraction,..spec()};
        fixture(s,|m,plain,source,mut heat,circuit| {
            let gate=CancelGate::new_clock_free();let (_,material,native)=native(m);
            let initial=native.total_energy_j().unwrap()+caloric(&heat)+0.5*INITIAL_CHARGE.powi(2)/CAPACITANCE;
            let mut owner=FrameRelaxationThermalSystem::new(material,native,&heat,plain,&DEST,budget(),&gate).unwrap();
            owner.bind_material_contact_circuit(&heat,source,circuit,&gate).unwrap();
            let routes=owner.electrical_contact_heat_destinations();
            assert_eq!(routes.len(),1);assert_eq!(routes[0].cells,[CellId(0),CellId(1)]);
            assert_eq!(routes[0].heat_fraction_to_first,fraction);
            let force=vec![0.0;owner.mechanical().mode_count()];let mut outside=0.0;
            for _ in 0..8 {
                let before=heat.cells().clone();let q0=owner.electrical().unwrap().charges_c()[0];
                let report=owner.step_electrical(&mut heat,&force,RelaxationElectricalInput {ports:&[],solve:electrical_controls()},
                    controls(),&[0.0;2],&gate).unwrap();
                let e=report.electrical.as_ref().unwrap();assert_eq!(e.contacts.len(),1);
                let rc=e.contacts[0].contact.resistance_ohm();
                assert!((rc-s.resistance*(1.0+0.02*(heat.temperature_k(CellId(0)).unwrap()-300.0))).abs()<1e-9);
                let r=e.resistors.iter().map(|r|r.material.conductor.resistance_ohm()).sum::<f64>()+rc;
                let a=1.0/(256.0*2.0*CAPACITANCE*r);
                let expected=q0*(1.0-a)/(1.0+a);let q1=owner.electrical().unwrap().charges_c()[0];
                assert!((q1-expected).abs()<1e-12);
                let current=(q0+q1)/(2.0*CAPACITANCE*r);
                let contact_heat=current*current*rc/256.0;
                assert!((e.contacts[0].heat_j-contact_heat).abs()<1e-12);
                for i in 0..2 {
                    let cell=CellId(i as u64);
                    let h_bulk=e.resistors.iter().find(|r|r.material.branch==i+1).unwrap().heat_j;
                    let h_relax=report.layers.iter().find(|r|r.destination.cell==cell).unwrap().heat_j;
                    let h_contact=contact_heat*if i==0 {fraction}else{1.0-fraction};
                    let delta=heat.cells()[&cell].mass_kg*(heat.cells()[&cell].specific_h_j_kg-before[&cell].specific_h_j_kg);
                    assert!((delta-h_bulk-h_relax-h_contact).abs()<1e-10);
                }
                assert_eq!(report.source.links[0].conductance_w_k,0.0);
                near(report.deposited_joule_heat_j,e.resistors.iter().map(|r|r.heat_j).sum::<f64>()+e.contacts[0].heat_j);
                outside+=report.external_mechanical_loss_j+report.mechanical.numerical_dissipation_j;
                near(owner.mechanical().total_energy_j().unwrap()+caloric(&heat)+e.stored_after_j-initial+outside,0.0);
                near(report.total_energy_residual_j,0.0);
            }
        });
    }
}

#[test]
fn changing_only_interface_source_changes_charge_both_temperatures_and_continuing_motion() {
    let mut outcomes=Vec::new();
    for resistance in [2.0,40.0] {
        fixture(Contact {resistance,..spec()},|m,plain,source,mut heat,circuit| {
            let gate=CancelGate::new_clock_free();let (placement,material,native)=native(m);
            let mut owner=FrameRelaxationThermalSystem::new(material,native,&heat,plain,&DEST,budget(),&gate).unwrap();
            owner.bind_material_contact_circuit(&heat,source,circuit,&gate).unwrap();
            let mut force=placement.project_nodal_force(&[0.0,0.0,0.0,0.0,0.002,0.0],&gate).unwrap();
            force.resize(owner.mechanical().mode_count(),0.0);
            for _ in 0..12 {owner.step_electrical(&mut heat,&force,RelaxationElectricalInput {ports:&[],solve:electrical_controls()},
                controls(),&[0.0;2],&gate).unwrap();}
            outcomes.push([placement.physical_state(owner.mechanical().components(),&gate).unwrap()[0],
                owner.electrical().unwrap().charges_c()[0],heat.temperature_k(CellId(0)).unwrap(),heat.temperature_k(CellId(1)).unwrap()]);
        });
    }
    for i in 0..4 {assert!((outcomes[1][i]-outcomes[0][i]).abs()>1e-12);}
    assert!(outcomes[1][1]>outcomes[0][1]);
}

#[test]
fn initial_area_and_late_interface_domain_refusal_preserve_complete_owners() {
    for narrow in [false,true] {
        fixture(Contact {upper_k:if narrow {300.000001}else{340.0},..spec()},|m,plain,source,mut heat,circuit| {
            let gate=CancelGate::new_clock_free();let (placement,material,native)=native(m);
            let mut owner=FrameRelaxationThermalSystem::new(material,native,&heat,plain,&DEST,budget(),&gate).unwrap();
            let original=&circuit.contacts()[0];
            let wrong=SeriesElectricalContact::new(original.id(),original.branch(),original.source().clone(),
                2.0*original.area_m2()).unwrap();
            let bad=MaterialCircuit::new_with_contacts(circuit.graph().clone(),circuit.resistors().to_vec(),vec![wrong],
                circuit.fluxes_wb().to_vec(),circuit.charges_c().to_vec(),&[]).unwrap();
            assert!(owner.bind_material_contact_circuit(&heat,source,bad,&gate).is_err());
            assert!(owner.electrical().is_none());
            owner.bind_material_contact_circuit(&heat,source,circuit.clone(),&gate).unwrap();
            let cells=heat.cells().clone();let memory=owner.mechanical().material_memory().to_vec();
            let physical=placement.physical_state(owner.mechanical().components(),&gate).unwrap();
            let material=owner.material().members()[0].source.clone();
            let force=vec![0.0;owner.mechanical().mode_count()];
            let mut exhausted=controls();exhausted.maximum_material_trials=1;
            assert!(owner.step_electrical(&mut heat,&force,RelaxationElectricalInput {ports:&[],solve:electrical_controls()},
                exhausted,&[0.0;2],&gate).is_err());
            let cancelled=CancelGate::new_clock_free();cancelled.request();
            assert!(owner.step_electrical(&mut heat,&force,RelaxationElectricalInput {ports:&[],solve:electrical_controls()},
                controls(),&[0.0;2],&cancelled).is_err());
            if narrow {assert!(owner.step_electrical(&mut heat,&force,RelaxationElectricalInput {ports:&[],solve:electrical_controls()},
                controls(),&[0.0;2],&gate).is_err());}
            assert_eq!(owner.electrical().unwrap(),&circuit);assert_eq!(heat.cells(),&cells);
            assert_eq!(owner.mechanical().material_memory(),memory);
            assert_eq!(owner.material().members()[0].source,material);
            assert_eq!(placement.physical_state(owner.mechanical().components(),&gate).unwrap(),physical);
            assert_eq!(owner.mechanical().samples_rendered(),0);assert_eq!(heat.time_s(),0.0);
            if !narrow {owner.step_electrical(&mut heat,&force,RelaxationElectricalInput {ports:&[],solve:electrical_controls()},
                controls(),&[0.0;2],&gate).unwrap();assert_eq!(owner.mechanical().samples_rendered(),1);}
        });
    }
}

#[test]
fn complete_heat_retains_original_electrical_interfaces_and_refuses_changed_partitions() {
    use crate::circuit_enthalpy::contacts::heat::MaterialHeatNetwork;
    use fs_exec::Cx;
    fixture(spec(),|m,plain,source,heat,circuit| {
        let gate=CancelGate::new_clock_free();let pool=fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
        pool.scope(|arena| {
            let cx=Cx::new(&gate,arena,fs_exec::StreamKey {seed:11,kernel_id:52,tile:0,iteration:0},
                fs_exec::Budget::INFINITE,fs_exec::ExecMode::Deterministic);
            let network=MaterialHeatNetwork {electrical:source,thermal_contacts:&[],radiation:None};
            let prepared=network.prepare(&cx,&heat).unwrap();
            let altered:Vec<_>=source.contacts.iter().map(|c|MaterialElectricalContact {link:c.link,branch:c.branch,
                point:c.point,selection:c.selection.clone(),heat_fraction_to_first:1.0-c.heat_fraction_to_first}).collect();
            let foreign=ElectricalContactNetwork {bulk:source.bulk,contacts:&altered,electrical_bases:source.electrical_bases};
            let wrong=MaterialHeatNetwork {electrical:&foreign,thermal_contacts:&[],radiation:None};
            let wrong=wrong.prepare(&cx,&heat).unwrap();
            let mut outcomes=Vec::new();
            for full in [false,true] {
                let mut heat=heat.clone();let (placement,material,native)=native(m);
                let mut owner=FrameRelaxationThermalSystem::new(material,native,&heat,plain,&DEST,budget(),&gate).unwrap();
                owner.bind_material_contact_circuit(&heat,source,circuit.clone(),&gate).unwrap();
                let force=vec![0.0;owner.mechanical().mode_count()];
                let memory=owner.mechanical().material_memory().to_vec();let cells=heat.cells().clone();
                assert!(owner.step_heat_network(&mut heat,&force,RelaxationHeatInput {network:&wrong,cx:&cx},
                    Some(RelaxationElectricalInput {ports:&[],solve:electrical_controls()}),controls(),&[0.0;2],&gate).is_err());
                assert_eq!(owner.electrical().unwrap(),&circuit);assert_eq!(heat.cells(),&cells);
                assert_eq!(owner.mechanical().material_memory(),memory);assert_eq!(owner.mechanical().samples_rendered(),0);
                for _ in 0..3 {
                    let input=RelaxationElectricalInput {ports:&[],solve:electrical_controls()};
                    if full {owner.step_heat_network(&mut heat,&force,RelaxationHeatInput {network:&prepared,cx:&cx},
                        Some(input),controls(),&[0.0;2],&gate).unwrap();}
                    else {owner.step_electrical(&mut heat,&force,input,controls(),&[0.0;2],&gate).unwrap();}
                }
                outcomes.push((heat.cells().clone(),owner.electrical().unwrap().clone(),
                    owner.mechanical().material_memory().to_vec(),
                    placement.physical_state(owner.mechanical().components(),&gate).unwrap()));
            }
            assert_eq!(outcomes[0],outcomes[1]);
        })
    });
}
