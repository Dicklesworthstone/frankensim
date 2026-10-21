//! Canonical gas packets preserve the original material evaluation and receipts.
use fs_material::gas_species::input::{GasRecord,ThermochemicalPhase,ThermochemicalSources,read_gas_record};
use fs_material::gas_species::input::packs::{SourcePackBytes,read_source_packs,read_source_packet};
#[path = "../src/render/schedule/reed/field_input/tests/support.rs"]
mod support;
use support::{Packs,packs};

fn encoded(source:&Packs)->[Vec<u8>;3] {
    [source.species.to_bytes(),source.nasa9.to_bytes(),source.transport.to_bytes()]
}
fn sources(bytes:&[[Vec<u8>;3]])->Vec<SourcePackBytes<'_>> {
    bytes.iter().map(|b|SourcePackBytes{species:&b[0],nasa9:&b[1],transport:&b[2]}).collect()
}
fn packet(bytes:&[[Vec<u8>;3]])->(Vec<u8>,Vec<u32>) {
    (bytes.iter().flatten().flatten().copied().collect(),
        bytes.iter().flatten().map(|b|u32::try_from(b.len()).unwrap()).collect())
}
fn read(text:&str,bytes:&[[Vec<u8>;3]])->GasRecord {
    read_source_packs(text,&sources(bytes),None,ThermochemicalPhase::FrozenGas).unwrap()
}
fn receipts(record:&GasRecord)->&[fs_material::thermochemical_gas::packs::ThermochemicalPackReceipt] {
    let GasRecord::Thermochemical(gas)=record else {panic!("original thermochemical source required")};
    gas.source_packs()
}

#[test]
fn arbitrary_species_packets_match_the_borrowed_registry_and_outlive_input_bytes() {
    for basis in ["mole","mass"] {
        let (actual,reference)= {
            let a=packs("synthetic-alpha",0.028,3.5,200.0);
            let b=packs("synthetic-beta",0.044,4.0,200.0);
            let text=format!("thermochemical 300 100000 {basis} 2 synthetic-beta 0.3 synthetic-alpha 0.7");
            let bindings=[a.binding(),b.binding()];
            let registry=ThermochemicalSources{species:&bindings,water_species:None,
                phase:ThermochemicalPhase::FrozenGas};
            let reference=read_gas_record(&mut text.split_ascii_whitespace(),None,Some(&registry)).unwrap();
            let bytes=[encoded(&a),encoded(&b)];
            assert_eq!(read(&text,&bytes),reference);
            let (payload,lengths)=packet(&bytes);
            let actual=read_source_packet(&text,&payload,&lengths,None,ThermochemicalPhase::FrozenGas).unwrap();
            assert_eq!(receipts(&actual)[0].species_pack,a.species.content_hash());
            assert_eq!(receipts(&actual)[0].nasa9_pack,a.nasa9.content_hash());
            assert_eq!(receipts(&actual)[0].transport_pack,a.transport.content_hash());
            (actual,reference)
        }; // Every source pack, packet byte and borrowed registry has been dropped.
        assert_eq!(actual,reference);
        actual.state().validate_coherence().unwrap();
        assert_eq!(receipts(&actual).len(),2);
    }
}

#[test]
fn source_properties_and_composition_change_the_same_physical_gas() {
    let text="thermochemical 300 100000 mole 1 synthetic-alpha 1";
    let baseline=read(text,&[encoded(&packs("synthetic-alpha",0.028,3.5,200.0))]);
    let heavier=read(text,&[encoded(&packs("synthetic-alpha",0.056,3.5,200.0))]);
    let higher_cp=read(text,&[encoded(&packs("synthetic-alpha",0.028,4.5,200.0))]);
    assert!((heavier.state().density/baseline.state().density-2.0).abs()<1e-13);
    assert!(heavier.state().sound_speed<baseline.state().sound_speed);
    assert_eq!(higher_cp.state().density,baseline.state().density);
    assert!(higher_cp.state().specific_heat_cp>baseline.state().specific_heat_cp);
    assert!(higher_cp.state().sound_speed<baseline.state().sound_speed);
    assert_eq!(higher_cp.state().dynamic_viscosity,baseline.state().dynamic_viscosity);
    assert_ne!(receipts(&heavier)[0].species_pack,receipts(&baseline)[0].species_pack);
    assert_ne!(receipts(&higher_cp)[0].nasa9_pack,receipts(&baseline)[0].nasa9_pack);
    assert_eq!(receipts(&higher_cp)[0].transport_pack,receipts(&baseline)[0].transport_pack);
    let bytes=[encoded(&packs("synthetic-alpha",0.028,3.5,200.0)),
        encoded(&packs("synthetic-beta",0.044,4.0,200.0))];
    let a=read("thermochemical 300 100000 mole 2 synthetic-alpha 0.8 synthetic-beta 0.2",&bytes);
    let b=read("thermochemical 300 100000 mole 2 synthetic-alpha 0.2 synthetic-beta 0.8",&bytes);
    assert_eq!(receipts(&a),receipts(&b));
    assert!(b.state().density>a.state().density);
    assert!(b.state().sound_speed<a.state().sound_speed);
}

#[test]
fn framing_missing_or_ambiguous_sources_and_source_role_mismatch_refuse() {
    let bytes=[encoded(&packs("synthetic-alpha",0.028,3.5,200.0))];
    let text="thermochemical 300 100000 mole 1 synthetic-alpha 1";
    let (payload,lengths)=packet(&bytes);
    let admit=|data:&[u8],sizes:&[u32]|read_source_packet(text,data,sizes,None,ThermochemicalPhase::FrozenGas);
    let mut trailing=payload.clone();trailing.push(0);
    let mut too_long=lengths.clone();too_long[2]+=1;
    let mut too_short=lengths.clone();too_short[2]-=1;
    let mut shifted=lengths.clone();shifted[0]+=1;shifted[1]-=1;
    assert!(admit(&[],&[]).is_err());
    assert!(admit(&payload,&lengths[..2]).is_err());
    assert!(admit(&payload,&[0,lengths[1],lengths[2]]).is_err());
    assert!(admit(&payload[..payload.len()-1],&lengths).is_err());
    assert!(admit(&trailing,&lengths).is_err());
    for wrong in [&too_long,&too_short,&shifted] {assert!(admit(&payload,wrong).is_err());}
    let duplicate=[sources(&bytes)[0];2];
    assert!(read_source_packs(text,&duplicate,None,ThermochemicalPhase::FrozenGas).is_err());
    let mut wrong=sources(&bytes);wrong[0].species=wrong[0].nasa9;
    assert!(read_source_packs(text,&wrong,None,ThermochemicalPhase::FrozenGas).is_err());
    let other=encoded(&packs("synthetic-beta",0.044,4.0,200.0));
    let mut mismatched=sources(&bytes);mismatched[0].nasa9=&other[1];
    assert!(read_source_packs(text,&mismatched,None,ThermochemicalPhase::FrozenGas).is_err());
    for record in [format!("{text} trailing"),text.replace("synthetic-alpha","missing-species"),
        "293.15 100000 0".into()] {
        assert!(read_source_packet(&record,&payload,&lengths,None,ThermochemicalPhase::FrozenGas).is_err());
    }
}

#[test]
fn packet_admission_preserves_explicit_phase_water_identity_and_source_domains() {
    let bytes=[encoded(&packs("synthetic-alpha",0.028,3.5,200.0)),
        encoded(&packs("synthetic-water",0.018,4.0,200.0))];
    let (payload,lengths)=packet(&bytes);
    let wet="thermochemical 300 100000 mole 2 synthetic-alpha 0.5 synthetic-water 0.5";
    let frozen=read_source_packet(wet,&payload,&lengths,Some("synthetic-water"),ThermochemicalPhase::FrozenGas).unwrap();
    assert_eq!(frozen.state().water_mole_fraction,0.5);
    for (record,water) in [(wet,Some("synthetic-water")),(wet,None),
        ("thermochemical 199 100000 mole 1 synthetic-alpha 1",Some("synthetic-water"))] {
        assert!(read_source_packet(record,&payload,&lengths,water,ThermochemicalPhase::LiquidWaterLimit).is_err());
    }
    let dry="thermochemical 250 100000 mole 1 synthetic-alpha 1";
    let checked=read_source_packet(dry,&payload,&lengths,Some("synthetic-water"),ThermochemicalPhase::LiquidWaterLimit).unwrap();
    assert_eq!(checked,read_source_packet(dry,&payload,&lengths,Some("synthetic-water"),ThermochemicalPhase::FrozenGas).unwrap());
    assert_eq!(checked.state().water_mole_fraction,0.0);
    assert_eq!(receipts(&checked).len(),1);
}
