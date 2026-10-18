//! Original resistance law with independent thermal states, not measured contact qualification.
#[path = "contact_material_state/common.rs"]
mod common;
use common::*;
use fs_conduction::InterfaceResistance;
use fs_matdb::{ClaimSet, SelectionPolicy};
use fs_qty::{Dims, QuantitySpec, Temperature};

#[test]
fn phase_changes_resistance_at_fixed_kelvin_and_retains_complete_source_point() {
    let cards=[bulk("A"),bulk("B")]; let curves=[curve(&cards[0]),curve(&cards[1])];
    let (source,id)=source(&cards,1.0); let context=point(); let binding=binding(&cards,&source,&context);
    for (h,expected) in [(25.0,1.25),(75.0,1.75)] {
        let states=[curves[0].state_at_specific_enthalpy(h).unwrap(),
            curves[1].state_at_specific_enthalpy(-50.0).unwrap()];
        let point=binding.at(states).unwrap();
        assert_eq!(point.axes()["aT"],400.0); assert_eq!(point.axes()["bT"],350.0);
        assert_eq!(point.axes()["T"],375.0); assert_eq!(point.axes()["pressure"],200000.0);
        assert_eq!(point.axes()["lot"],7.0); assert_eq!(point.axes()["fa"],h/100.0);
        assert_eq!(point.axis_quantities()["aT"],QuantitySpec::dimensional(Temperature::DIMS));
        let resistance=InterfaceResistance::from_card_pinned("contact",&source,&point,id).unwrap();
        assert_eq!(resistance.value_m2_k_per_w(),expected);
        assert_eq!(resistance,InterfaceResistance::from_card_pinned("contact",&source,&point,id).unwrap());
        assert_eq!(resistance.card_identity(),source.content_hash());
    }
    assert_eq!(context.axes()["fa"],0.0);
}

#[test]
fn ordered_surface_and_chart_identity_and_axis_aliasing_are_not_silently_repaired() {
    let cards=[bulk("A"),bulk("B")]; let curves=[curve(&cards[0]),curve(&cards[1])];
    let (source,_)=source(&cards,1.0); let context=point(); let mut binding=binding(&cards,&source,&context);
    let states=[curves[0].state_at_specific_enthalpy(25.0).unwrap(),curves[1].state_at_specific_enthalpy(-50.0).unwrap()];
    assert!(binding.at([states[1],states[0]]).is_err());
    binding.materials.swap(0,1); assert!(binding.at([states[1],states[0]]).is_err());
    binding.materials.swap(0,1); binding.temperature_axes=["aT","aT"];
    assert!(binding.at(states).is_err()); binding.temperature_axes=["aT","bT"];
    binding.mean_temperature_axis=Some("aT"); assert!(binding.at(states).is_err());
    binding.mean_temperature_axis=Some("T"); binding.liquid_fraction_axes=[Some("fa"),Some("fa")];
    assert!(binding.at(states).is_err());
    let bad=fs_matdb::QueryPoint::new().with_quantity("aT",QuantitySpec::dimensional(Dims::NONE),400.0).unwrap();
    binding.context=&bad; binding.liquid_fraction_axes=[Some("fa"),Some("fb")];
    assert!(binding.at(states).is_err());
}

#[test]
fn pins_domains_and_static_source_support_remain_load_bearing() {
    let cards=[bulk("A"),bulk("B")]; let curves=[curve(&cards[0]),curve(&cards[1])];
    let states=[curves[0].state_at_specific_enthalpy(25.0).unwrap(),curves[1].state_at_specific_enthalpy(-50.0).unwrap()];
    let mut claims=ClaimSet::new(); let first=claims.insert_claim(claim("fa",Dims::NONE,1.0)).unwrap();
    claims.insert_claim(claim("fa",Dims::NONE,2.0)).unwrap();
    let source=interface(&cards,claims); let context=point(); let b=binding(&cards,&source,&context);
    let query=b.at(states).unwrap();
    assert!(InterfaceResistance::from_card("contact",&source,&query,SelectionPolicy::SingleClaimOnly).is_err());
    assert_eq!(InterfaceResistance::from_card_pinned("contact",&source,&query,first)
        .unwrap().value_m2_k_per_w(),1.25);
    let wrong_pressure=context.clone().with("pressure",100000.0).unwrap();
    let wrong_query=binding(&cards,&source,&wrong_pressure).at(states).unwrap();
    assert!(InterfaceResistance::from_card_pinned("contact",&source,&wrong_query,first).is_err());
    assert!(b.is_thermal_axis("fa",Dims::NONE));
    assert!(b.is_thermal_axis("aT",Temperature::DIMS));
    assert!(b.is_thermal_axis("T",Temperature::DIMS));
    assert!(!b.is_thermal_axis("frequency",Dims::NONE));
    assert!(!b.is_thermal_axis("aT",Dims::NONE));
}
