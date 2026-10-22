//! G0/G3: retain actual source claims, not stale coefficients, across heat-network changes.
use super::*;
use fs_conduction::SURFACE_EMISSIVITY_PROPERTY;
use fs_conduction::transient::deposition::radiation::{MaterialRadiation, RadiatingCell};
use fs_couple::circuit_enthalpy::contacts::heat::{HeatNetworkError, MaterialHeatNetwork};
use fs_material::electrical_interface::ELECTRICAL_CONTACT_RESISTANCE_PROPERTY;
use fs_material::state_point::MaterialPropertySelection;
use fs_matdb::ClaimSelection;

fn remap<'a>(old: &RadiatingCell<'a>, surface: usize, selection: MaterialPropertySelection) -> RadiatingCell<'a> {
    RadiatingCell { surface, cell: old.cell, material_id: old.material_id,
        curve: old.curve, material: old.material, selection }
}

#[test]
fn source_pins_survive_surface_reordering_and_reach_the_next_physical_interval() {
    with_case(0.7,0.5,true,false,vec![],|m,cx,_,network,t,b,c| {
        let original=network.prepare(cx,t).unwrap();
        let before=original.resolve(cx,t).unwrap();
        // A foreign-property claim would be invalid for both heat mechanisms.
        // Continuation must reuse the old pins, not try successor policies first.
        let foreign=before.contacts[0].source().resolved()
            .property(ELECTRICAL_CONTACT_RESISTANCE_PROPERTY).unwrap().answer().receipt.selected;
        let contacts=network.thermal_contacts.iter().map(|old| {
            let mut next=*old; next.selection=ClaimSelection::Pinned(foreign); next
        }).collect::<Vec<_>>();
        let old_rad=network.radiation.unwrap();
        let surfaces=old_rad.cells.iter().map(|old|remap(old,1-old.surface,
            MaterialPropertySelection::PinnedByProperty(vec![(SURFACE_EMISSIVITY_PROPERTY.into(),foreign)])))
            .collect::<Vec<_>>();
        let radiation=MaterialRadiation { enclosure:old_rad.enclosure,cells:&surfaces,
            power_tolerance_w:old_rad.power_tolerance_w,radiosity_tolerance_w_m2:old_rad.radiosity_tolerance_w_m2 };
        let successor=MaterialHeatNetwork { electrical:network.electrical,
            thermal_contacts:&contacts,radiation:Some(&radiation) };
        assert!(successor.prepare(cx,t).is_err());
        let retained=original.extend_fixed_circuit(cx,&successor,t,&[Some(1),Some(0)]).unwrap();
        let after=retained.resolve(cx,t).unwrap();
        assert_eq!(before.bulk.source.thermal_contacts,after.bulk.source.thermal_contacts);
        let left=before.bulk.source.radiation.as_ref().unwrap();
        let right=after.bulk.source.radiation.as_ref().unwrap();
        for surface in left.surfaces() {
            let mapped=&right.surfaces()[1-surface.surface];
            assert_eq!(surface.cell,mapped.cell);
            assert_eq!(surface.material,mapped.material);
            near(left.sources_w()[&surface.cell],right.sources_w()[&surface.cell],1e-12);
        }
        let reference=step_heat_network_candidate(cx,t,b,&map(m,false),c,&original,&[],DT,
            1.0,&[0.02],controls(),&[0.0;2]).unwrap();
        let actual=step_heat_network_candidate(cx,t,b,&map(m,false),c,&retained,&[],DT,
            1.0,&[0.02],controls(),&[0.0;2]).unwrap();
        for id in CELLS { near(actual.thermal().temperature_k(id).unwrap(),reference.thermal().temperature_k(id).unwrap(),1e-7); }
        near(actual.circuit().charges_c()[0],reference.circuit().charges_c()[0],1e-12);
        near(actual.bridge().deformed_position_m(LENGTH,0.0).unwrap()[1],
            reference.bridge().deformed_position_m(LENGTH,0.0).unwrap()[1],1e-10);
        assert_eq!(t.time_s(),0.0); assert_eq!(c.charges_c(),&[CHARGE]);
    });
}

#[test]
fn maps_must_be_complete_injective_and_preserve_physical_receivers() {
    with_case(0.7,0.5,true,false,vec![],|_,cx,_,network,t,_,_| {
        let original=network.prepare(cx,t).unwrap();
        for mapping in [vec![],vec![Some(0)],vec![Some(0),Some(0)],vec![Some(9),Some(1)],vec![Some(1),Some(0)]] {
            assert!(original.extend_fixed_circuit(cx,network,t,&mapping).is_err());
        }
        let missing=MaterialHeatNetwork { electrical:network.electrical,thermal_contacts:&[],radiation:network.radiation };
        assert!(original.extend_fixed_circuit(cx,&missing,t,&[Some(0),Some(1)]).is_err(),
            "an active interface cannot silently lose series resistance");
        let retired=MaterialHeatNetwork { electrical:network.electrical,
            thermal_contacts:network.thermal_contacts,radiation:None };
        assert!(original.extend_fixed_circuit(cx,&retired,t,&[Some(0),Some(1)]).is_err());
        let next=original.extend_fixed_circuit(cx,&retired,t,&[None,None]).unwrap();
        assert!(next.resolve(cx,t).unwrap().bulk.source.radiation.is_none());
        assert!(original.resolve(cx,t).unwrap().bulk.source.radiation.is_some());
        assert_eq!(t.time_s(),0.0);
    });
}

#[test]
fn new_surfaces_and_interfaces_use_authored_policies_not_unrelated_old_pins() {
    with_case(0.7,0.5,true,false,vec![],|_,cx,_,network,t,_,_| {
        let bare=MaterialHeatNetwork { electrical:network.electrical,thermal_contacts:&[],radiation:None };
        let original=bare.prepare(cx,t).unwrap();
        let added=original.extend_fixed_circuit(cx,network,t,&[]).unwrap();
        let expected=network.prepare(cx,t).unwrap().resolve(cx,t).unwrap();
        let actual=added.resolve(cx,t).unwrap();
        assert_eq!(actual.bulk.source.thermal_contacts,expected.bulk.source.thermal_contacts);
        assert_eq!(actual.bulk.source.radiation,expected.bulk.source.radiation);
        let wrong=actual.bulk.source.thermal_contacts[0].interface.receipt().selected;
        let old_rad=network.radiation.unwrap();
        let surfaces=old_rad.cells.iter().map(|s|remap(s,s.surface,
            MaterialPropertySelection::PinnedByProperty(vec![(SURFACE_EMISSIVITY_PROPERTY.into(),wrong)])))
            .collect::<Vec<_>>();
        let radiation=MaterialRadiation { enclosure:old_rad.enclosure,cells:&surfaces,
            power_tolerance_w:old_rad.power_tolerance_w,radiosity_tolerance_w_m2:old_rad.radiosity_tolerance_w_m2 };
        let invalid=MaterialHeatNetwork { electrical:network.electrical,
            thermal_contacts:network.thermal_contacts,radiation:Some(&radiation) };
        assert!(original.extend_fixed_circuit(cx,&invalid,t,&[]).is_err());
        assert!(original.resolve(cx,t).unwrap().bulk.source.radiation.is_none());
    });
}

#[test]
fn changed_fixed_conditions_and_scope_cancellation_leave_old_sources_usable() {
    with_case(0.7,0.5,true,false,vec![],|_,cx,gate,network,t,b,c| {
        let original=network.prepare(cx,t).unwrap();
        let old_rad=network.radiation.unwrap();
        let changed=old_rad.cells[0].material.context.clone().with("lot",9.0).unwrap();
        let mut surfaces=old_rad.cells.iter().map(|s|remap(s,s.surface,s.selection.clone())).collect::<Vec<_>>();
        surfaces[0].material.context=&changed;
        let radiation=MaterialRadiation { enclosure:old_rad.enclosure,cells:&surfaces,
            power_tolerance_w:old_rad.power_tolerance_w,radiosity_tolerance_w_m2:old_rad.radiosity_tolerance_w_m2 };
        let invalid=MaterialHeatNetwork { electrical:network.electrical,
            thermal_contacts:network.thermal_contacts,radiation:Some(&radiation) };
        assert!(original.extend_fixed_circuit(cx,&invalid,t,&[Some(0),Some(1)]).is_err());
        let intact=original.resolve(cx,t).unwrap();
        assert_eq!(intact.bulk.source.thermal_contacts.len(),1);
        gate.request();
        assert!(matches!(original.extend_fixed_circuit(cx,network,t,&[Some(0),Some(1)]),Err(HeatNetworkError::Cancelled)));
        assert_eq!(t.time_s(),0.0); assert_eq!(b.time_s(),0.0); assert_eq!(c.charges_c(),&[CHARGE]);
    });
}
