//! Source-preserving electrical interfaces through real material deposition.
use super::*;
use fs_evidence::ValidityDomain;
use fs_material::electrical_interface::{
    ELECTRICAL_CONTACT_RESISTANCE_DIMS as RESISTANCE_DIMS,
    ELECTRICAL_CONTACT_RESISTANCE_PROPERTY as RESISTANCE,
};
use fs_matdb::{ClaimSet, InterfaceSystemCard, InterpolationPolicy, PropertyClaim,
    PropertyValue, Provenance, SurfaceSpec, SystemContext, UncertaintyModel};
use fs_qty::Temperature;
use fs_couple::circuit_enthalpy::contacts::material::ElectricalContactNetworkError;

fn interface(m: &[Material;2], expired: bool) -> InterfaceSystemCard {
    let mut claims=ClaimSet::new();
    for (lo,hi,r) in [(280.0,if expired {300.02}else{320.0},1.0),(300.01,320.0,2.0)] {
        claims.insert_claim(PropertyClaim {
            key:PropertyKey::new(RESISTANCE,RESISTANCE_DIMS),
            value:PropertyValue::Curve {abscissa:"Ta".into(),abscissa_dims:Temperature::DIMS,
                knots:vec![(lo,r*(1.0+0.005*(lo-300.0))),(hi,r*(1.0+0.005*(hi-300.0)))],
                dims:RESISTANCE_DIMS},
            validity:ValidityDomain::unconstrained().with("Ta",lo,hi).with("Tb",280.0,320.0)
                .with("pressure",100000.0,100000.0),
            uncertainty:UncertaintyModel::Unstated,interpolation:InterpolationPolicy::LinearInside,
            observations:vec![],provenance:Provenance {source:"synthetic competing contact sources".into(),
                license:"CC0-1.0".into(),artifact:None},
        }).unwrap();
    }
    InterfaceSystemCard::assemble(
        SurfaceSpec {material:m[0].card.id().clone(),texture_frame:"first finish".into()},
        SurfaceSpec {material:m[1].card.id().clone(),texture_frame:"second finish".into()},
        SystemContext {medium:"authored contact".into(),third_body:None,
            environment:"test".into(),history:"fixed".into()},claims,vec![]).unwrap()
}

fn with_successor<R>(m:&[Material;2],card:&InterfaceSystemCard,grown:bool,
    run:impl FnOnce(&ElectricalContactNetwork<'_>)->R)->R {
    let context=point(300.0); let key=key();
    let ids=if grown {vec![(99,1),(55,1),(90,0)]}else{vec![(55,1),(90,0)]};
    let bound=cells(m,&context,&key,&ids); let resistors=resistors();
    let links=[ConductionLink {link:7,distance_a_m:0.0005,distance_b_m:0.0005}];
    let bulk=MaterialNetwork {cells:&bound,resistors:&resistors,links:&links};
    let context=contact_source::point();
    let contacts=[MaterialElectricalContact {link:7,branch:1,
        point:ThermalInterfacePoint {interface:card,materials:[&m[0].card,&m[1].card],
            context:&context,temperature_axes:["Ta","Tb"],
            liquid_fraction_axes:[Some("melt_a"),Some("melt_b")],mean_temperature_axis:None},
        selection:MaterialPropertySelection::SingleClaimOnly,heat_fraction_to_first:0.25}];
    let bases=[(1,ElectricalPropertyBasis::Resistivity)];
    run(&ElectricalContactNetwork {bulk:&bulk,contacts:&contacts,electrical_bases:&bases})
}

#[test]
fn retained_interface_claim_survives_birth_competition_and_further_heating_and_sag() {
    with_selected_contact(false,|m|interface(m,false),|m,old,circuit,p,card| {
        let original=p.resolve(old.thermal(),||false).unwrap();
        let claim=original.contacts[0].source().resolved().property(RESISTANCE).unwrap().answer().receipt.selected;
        let (born,charged,r)=old.step_electrical_with_contacts_candidate(circuit,p,&[],0.2,
            1.0,&[0.02],controls(),&[0.0;2],||false).unwrap();
        assert!(born.thermal().temperature_k(CellId(90)).unwrap()>300.01);
        with_successor(m,card,true,|network| {
            assert!(network.prepare(born.thermal(),||false).is_err(),"competing contact fits are now eligible");
            let continued=p.extend_fixed_circuit(network,born.thermal(),||false).unwrap();
            let sample=continued.resolve(born.thermal(),||false).unwrap();
            assert_eq!(sample.contacts[0].source().resolved().property(RESISTANCE).unwrap().answer().receipt.selected,claim);
            near(sample.contacts[0].resistance_ohm(),1000.0*(1.0+0.005*
                (born.thermal().temperature_k(CellId(90)).unwrap()-300.0)),1e-8);
            let (after,final_circuit,r2)=born.step_electrical_with_contacts_candidate(&charged,&continued,&[],
                0.4,1.0,&[0.02;2],controls(),&[0.0;3],||false).unwrap();
            assert!(final_circuit.charges_c()[0]<charged.charges_c()[0]);
            assert_eq!(final_circuit.contacts().len(),1);
            assert!(r2.interval.layer_energy_j[&(3,LAYERS[0])].dissipated>0.0);
            assert!(!r2.birth_equilibrium_pending);
            near(after.thermal().energy_j()+after.bridge().stored_energy_j()+r2.interval.electrical.stored_after_j
                +r.interval.mechanical.energy_j.numerical_dissipation+r2.interval.mechanical.energy_j.numerical_dissipation
                -r.interval.mechanical.load_work_j-r2.interval.mechanical.load_work_j,
                Q0*Q0/(2.0*CAP),2.0*config().total_energy_tolerance_j);
        });
        assert_eq!(old.thermal().time_s(),0.0); assert_eq!(circuit.charges_c(),&[Q0]);
    });
}

#[test]
fn expired_interface_claim_refuses_instead_of_reselecting_a_hot_fit() {
    with_selected_contact(false,|m|interface(m,true),|m,old,circuit,p,card| {
        // Use the existing independent heat source to leave electrical validity.
        let (warm,_)=old.step_candidate(0.1,0.0,&[0.0],controls().bridge,&[0.0;2],|t| {
            let mut c=coefficients(t)?; c.sources_w.insert(CellId(90),0.01);
            Ok::<_,core::convert::Infallible>(c)
        },||false).unwrap();
        assert!(warm.thermal().temperature_k(CellId(90)).unwrap()>300.02);
        with_successor(m,card,false,|network| {
            assert!(network.prepare(warm.thermal(),||false).is_ok(),"only the replacement fit is eligible");
            assert!(p.extend_fixed_circuit(network,warm.thermal(),||false).is_err());
        });
        assert!(p.resolve(old.thermal(),||false).is_ok());
        assert_eq!(circuit.charges_c(),&[Q0]); assert_eq!(old.thermal().time_s(),0.0);
    });
}

#[test]
fn changed_contact_meaning_is_not_a_source_preserving_extension() {
    with_selected_contact(false,|m|interface(m,false),|m,old,_,p,card| {
        with_successor(m,card,false,|network| {
            let other=contact_source::card(&m[0].card,&m[1].card,3.0,0.0);
            let pressure=contact_source::point().with("pressure",200000.0).unwrap();
            for mode in 0..6 {
                let b=&network.contacts[0];
                let mut changed=MaterialElectricalContact {link:b.link,branch:b.branch,point:b.point,
                    selection:b.selection.clone(),heat_fraction_to_first:b.heat_fraction_to_first};
                match mode {
                    0=>changed.link=8,
                    1=>changed.branch=2,
                    2=>changed.heat_fraction_to_first=0.75,
                    3=>changed.point.interface=&other,
                    4=>changed.point.context=&pressure,
                    _=>changed.point.materials.swap(0,1),
                }
                let contacts=[changed];
                let next=ElectricalContactNetwork {bulk:network.bulk,contacts:&contacts,
                    electrical_bases:network.electrical_bases};
                assert!(p.extend_fixed_circuit(&next,old.thermal(),||false).is_err());
            }
            let removed=ElectricalContactNetwork {bulk:network.bulk,contacts:&[],electrical_bases:network.electrical_bases};
            assert!(p.extend_fixed_circuit(&removed,old.thermal(),||false).is_err());
            let bases=[(1,ElectricalPropertyBasis::Conductivity)];
            let changed=ElectricalContactNetwork {bulk:network.bulk,contacts:network.contacts,electrical_bases:&bases};
            assert!(p.extend_fixed_circuit(&changed,old.thermal(),||false).is_err());
            let bases=[(1,ElectricalPropertyBasis::Resistivity);2];
            let duplicate=ElectricalContactNetwork {bulk:network.bulk,contacts:network.contacts,electrical_bases:&bases};
            assert!(p.extend_fixed_circuit(&duplicate,old.thermal(),||false).is_err());
        });
    });
}

#[test]
fn stale_contact_area_refuses_while_reordering_and_final_cancellation_are_transactional() {
    with_selected_contact(false,|m|interface(m,false),|m,old,_,p,card| {
        let initial=p.resolve(old.thermal(),||false).unwrap();
        let mut events=old.thermal().cells().values().cloned().map(|cell|
            TimedEvent {time_s:0.0,event:Event::Birth(cell)}).collect::<Vec<_>>();
        let mut link=old.thermal().links()[&7].clone(); link.area_m2*=0.5;
        events.push(TimedEvent {time_s:0.0,event:Event::Contact(link)});
        let changed=DepositionSimulation::new(events,&[
            (MaterialId(0),m[0].caloric.curve() as &dyn ThermalLaw),
            (MaterialId(1),m[1].caloric.curve() as &dyn ThermalLaw)],
            DepositionBudget {max_cells:8,max_links:8,max_events:16,max_substeps:100}).unwrap()
            .advance_to(0.0,StepConfig {max_dt_s:0.1,energy_tolerance_j:1e-9},||Ok(())).unwrap().0;
        assert!(p.resolve(&changed,||false).is_err());
        with_successor(m,card,false,|network| {
            assert!(p.extend_fixed_circuit(network,&changed,||false).is_err());
            let mut polls=0;
            let next=p.extend_fixed_circuit(network,old.thermal(),||{polls+=1;false}).unwrap();
            let actual=next.resolve(old.thermal(),||false).unwrap();
            assert_eq!(actual.bulk.source,initial.bulk.source);
            assert_eq!(actual.bulk.resistors,initial.bulk.resistors);
            assert_eq!(actual.contacts,initial.contacts);
            let mut count=0;
            assert!(matches!(p.extend_fixed_circuit(network,old.thermal(),||{count+=1;count==polls}),
                Err(ElectricalContactNetworkError::Thermal(DepositionError::Cancelled))));
            let replay=p.extend_fixed_circuit(network,old.thermal(),||false).unwrap();
            assert_eq!(replay.resolve(old.thermal(),||false).unwrap().contacts,initial.contacts);
        });
        assert_eq!(p.resolve(old.thermal(),||false).unwrap().contacts,initial.contacts);
    });
}
