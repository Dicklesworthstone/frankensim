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

/// G1/G3: the finite input retains source calorics, not only initial coefficients.
#[test]
fn finite_parcel_ingress_retains_sources_for_nonlinear_heating_after_bytes_are_dropped() {
    use fs_material::gas_species::input::read_isochoric_parcel;
    use fs_material::gas_species::input::packs::read_source_parcel;
    use fs_material::thermochemical_gas::energy::GasEnergySolve;
    use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;
    for basis in ["mole", "mass"] {
        let (mut parcel, expected_receipts) = {
            let a = packs("synthetic-alpha", 0.028, 3.5, 200.0);
            let b = packs("synthetic-beta", 0.044, 4.0, 200.0);
            let text = format!(
                "thermochemical 300 100000 {basis} 2 synthetic-beta 0.3 synthetic-alpha 0.7"
            );
            let bindings = [a.binding(), b.binding()];
            let registry = ThermochemicalSources {
                species: &bindings,
                water_species: None,
                phase: ThermochemicalPhase::FrozenGas,
            };
            let borrowed = read_isochoric_parcel(
                &mut text.split_ascii_whitespace(), &registry, 0.002,
            ).unwrap();
            let bytes = [encoded(&a), encoded(&b)];
            let owned = read_source_parcel(
                &text, &sources(&bytes), None, ThermochemicalPhase::FrozenGas, 0.002,
            ).unwrap();
            assert_eq!(owned, borrowed);
            let GasRecord::Thermochemical(snapshot) = read(&text, &bytes)
                else { panic!("original thermochemical state required"); };
            assert_eq!(owned.state(), &snapshot);
            (owned, snapshot.source_packs().to_vec())
        }; // All packs, bytes, borrowed bindings and the snapshot have gone away.
        let xa = if basis == "mole" { 0.7 }
            else { (0.7 / 0.028) / (0.7 / 0.028 + 0.3 / 0.044) };
        let xb = 1.0 - xa;
        let molar_mass = xa * 0.028 + xb * 0.044;
        let cp_intercept = xa * 3.5 + xb * 4.0;
        let expected_mass = 100000.0 * 0.002 * molar_mass / (R * 300.0);
        assert!((parcel.mass_kg() / expected_mass - 1.0).abs() < 1e-13);
        let old = parcel.clone();
        let target = 420.0_f64;
        // Fixture NASA-9 has cp/R = cp_intercept + 0.001*T.
        // Integrate cv(T), not the initial cp or a frozen-gamma approximation.
        let delta_u = R / molar_mass * ((cp_intercept - 1.0) * (target - 300.0)
            + 0.0005 * (target.powi(2) - 300.0_f64.powi(2)));
        let heat = expected_mass * delta_u;
        let controls = GasEnergySolve {
            minimum_temperature_k: 200.0,
            maximum_temperature_k: 800.0,
            energy_tolerance_j_kg: 1e-7,
            maximum_iterations: 128,
        };
        let step = parcel.add_heat(heat, controls, || false).unwrap();
        assert!((parcel.state().state().temperature - target).abs() < 1e-8);
        assert!((parcel.state().state().pressure / 140000.0 - 1.0).abs() < 1e-10);
        assert!((parcel.state().state().density / old.state().state().density - 1.0).abs()
            <= 64.0 * f64::EPSILON);
        assert_eq!(parcel.mass_kg(), old.mass_kg());
        assert_eq!(parcel.volume_m3(), old.volume_m3());
        assert_eq!(parcel.state().source_packs(), expected_receipts.as_slice());
        assert!(parcel.state().state().specific_heat_cp > old.state().state().specific_heat_cp);
        assert!(parcel.state().state().gamma < old.state().state().gamma);
        assert_eq!(step.supplied_heat_j, heat);
        assert!(step.residual_j.abs() <= parcel.mass_kg() * controls.energy_tolerance_j_kg);
        let measured = parcel.mass_kg() * (parcel.state().specific_internal_energy_j_kg()
            - old.state().specific_internal_energy_j_kg());
        assert!((measured - heat).abs() <= parcel.mass_kg() * controls.energy_tolerance_j_kg);
        // The same live source owner still refuses a request beyond its bracket.
        let accepted = parcel.clone();
        assert!(parcel.add_heat(1e9, controls, || false).is_err());
        assert_eq!(parcel, accepted);
    }
}

/// G4: even a late volume/source refusal must leave the caller's record intact.
#[test]
fn finite_parcel_prefix_and_refusals_are_transactional() {
    use fs_material::gas_species::input::{ParcelInputError, read_isochoric_parcel};
    let a = packs("synthetic-alpha", 0.028, 3.5, 200.0);
    let bindings = [a.binding()];
    let registry = ThermochemicalSources {
        species: &bindings, water_species: None, phase: ThermochemicalPhase::FrozenGas,
    };
    let valid = "thermochemical 300 100000 mole 1 synthetic-alpha 1";
    let text = format!("{valid} caller-owned 7");
    let mut fields = text.split_ascii_whitespace();
    read_isochoric_parcel(&mut fields, &registry, 0.001).unwrap();
    assert_eq!(fields.collect::<Vec<_>>(), ["caller-owned", "7"]);
    for (text, volume) in [
        (valid, 0.0), (valid, -1.0), (valid, f64::NAN), (valid, f64::INFINITY),
        ("300 100000 0", 0.001),
        ("mixture 300 100000 mole 1 N2 1", 0.001),
        ("thermochemical 300 100000 mole 1 missing 1", 0.001),
        ("thermochemical 199 100000 mole 1 synthetic-alpha 1", 0.001),
        ("thermochemical 300 500 mole 1 synthetic-alpha 1", 0.001),
        ("thermochemical 300 100000 mole 129", 0.001),
    ] {
        let mut fields = text.split_ascii_whitespace();
        assert!(read_isochoric_parcel(&mut fields, &registry, volume).is_err());
        assert_eq!(fields.collect::<Vec<_>>(), text.split_ascii_whitespace().collect::<Vec<_>>());
    }
    let ambiguous = [a.binding(), a.binding()];
    let absent_model = [fs_material::thermochemical_gas::packs::SpeciesPackRef {
        nasa9: None, ..a.binding()
    }];
    for sources in [&ambiguous[..], &absent_model[..]] {
        let bad = ThermochemicalSources {
            species: sources, water_species: None, phase: ThermochemicalPhase::FrozenGas,
        };
        let mut fields = valid.split_ascii_whitespace();
        assert!(matches!(read_isochoric_parcel(&mut fields, &bad, 0.001),
            Err(ParcelInputError::Input(_))));
        assert_eq!(fields.collect::<Vec<_>>(), valid.split_ascii_whitespace().collect::<Vec<_>>());
    }
}

/// G3/G4: snapshot phase checking never becomes an evolving condensation claim.
#[test]
fn finite_parcel_requires_frozen_phase_and_owned_sources_reject_suffixes() {
    use fs_material::gas_species::input::{ParcelInputError, read_isochoric_parcel};
    use fs_material::gas_species::input::packs::read_source_parcel;
    let a = packs("synthetic-alpha", 0.028, 3.5, 200.0);
    let bindings = [a.binding()];
    let checked = ThermochemicalSources {
        species: &bindings, water_species: None, phase: ThermochemicalPhase::LiquidWaterLimit,
    };
    let text = "thermochemical 300 100000 mole 1 synthetic-alpha 1";
    let mut fields = text.split_ascii_whitespace();
    assert!(matches!(read_isochoric_parcel(&mut fields, &checked, 0.001),
        Err(ParcelInputError::Input(_))));
    assert_eq!(fields.collect::<Vec<_>>(), text.split_ascii_whitespace().collect::<Vec<_>>());
    let bytes = [encoded(&a)];
    let source_bytes = sources(&bytes);
    for (record, phase) in [
        (text.to_owned(), ThermochemicalPhase::LiquidWaterLimit),
        (format!("{text} trailing"), ThermochemicalPhase::FrozenGas),
        ("300 100000 0".to_owned(), ThermochemicalPhase::FrozenGas),
    ] {
        assert!(read_source_parcel(&record, &source_bytes, Some("synthetic-water"), phase, 0.001).is_err());
    }
    let duplicate = [source_bytes[0]; 2];
    assert!(read_source_parcel(text, &duplicate, None, ThermochemicalPhase::FrozenGas, 0.001).is_err());
    let malformed = [SourcePackBytes { nasa9: &[], ..source_bytes[0] }];
    assert!(read_source_parcel(text, &malformed, None, ThermochemicalPhase::FrozenGas, 0.001).is_err());
    // The existing checked snapshot remains admitted for explicitly dry input.
    let snapshot = read_source_packs(
        text, &source_bytes, Some("synthetic-water"), ThermochemicalPhase::LiquidWaterLimit,
    ).unwrap();
    assert_eq!(snapshot, read(text, &bytes));
}
