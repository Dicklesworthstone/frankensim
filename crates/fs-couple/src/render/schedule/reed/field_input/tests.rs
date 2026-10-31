//! Source-to-field-to-pressure wiring checks with explicit synthetic packs.
use super::*;
use super::super::{
    ObservationPolicy, PressureRenderer, ReedBoreVoice, ReedPerformance, RenderVoice,
};
use fs_material::gas_species::input::ThermochemicalPhase;
use fs_material::thermochemical_gas::packs::{SpeciesPackRef, from_packs};
use fs_qty::semantic::{Composition, CompositionBasis};
mod support;
use support::packs;

const BLOCK: &str = "gas_field 8 3 mole 2 A B\ngas_station 0 300 1 0\ngas_station 0.285 320 0.75 0.25\ngas_station 0.57 310 0.5 0.5\n";
const INPUT_FIXTURE: &str = include_str!("../../../../../examples/reed-flanged-baffled.performance");
/// The fixture's physics without its sidecar declaration: these library tests
/// add their own observation records, and the library never refuses an
/// undeclared law (ruling 6532 puts that in the sidecar writer).
static INPUT: std::sync::LazyLock<String> =
    std::sync::LazyLock::new(|| INPUT_FIXTURE.replace("observer_law stokes-kirchhoff\n", ""));

fn parse(text: &str, sources: &ThermochemicalSources<'_>) -> (AxialGasField, Vec<GasRecord>) {
    let mut reader = Reader { lines: text.lines(), line: 0 };
    let mut records = Vec::new();
    let (field, resolution) = read(&mut reader, Some(sources), 100000.0, &mut records).unwrap();
    assert_eq!(resolution, GasFieldResolution::Midpoint(8));
    assert!(reader.lines.next().is_none());
    (field, records)
}
fn base() -> String {
    INPUT.replace("audio 48000 12000 2", "audio 48000 4801 2")
        .replace("ambient mixture 293.15 100000 mole 2 N2 0.7 O2 0.3", "ambient 300 100000 0")
        .replace("reed 0.0004 0.013 6000 0 0 0.35", "reed 0.0004 0.013 6000 0.00001 100 0.35")
        .replace("segments 2\ncylinder 0.0015 0.2\ncylinder 0.0015 0.3 gas mixture 293.15 100000 mole 3 N2 0.4 O2 0.2 CO2 0.4",
            "segments 1\ncylinder 0.0015 0.57")
}
fn with_block(block: &str) -> String {
    base().replace("compile_limits", &format!("{block}compile_limits"))
}
fn load(text: &str, sources: &ThermochemicalSources<'_>, block: usize) -> ReedPerformance {
    ReedPerformance::from_bytes_with_gas_sources(
        text.as_bytes(), block, ObservationPolicy::Diagnostic, Some(sources),
    ).unwrap()
}
fn voice(performance: &ReedPerformance) -> &ReedBoreVoice {
    let RenderVoice::ReedBore(voice) = &performance.renderer.context().voices[0] else {
        panic!("original reed was lost")
    };
    voice
}
fn motion(performance: &ReedPerformance) -> [u64;4] {
    let v = voice(performance);
    [v.reed_y, v.reed_v, v.p_plus_prev, v.line.incoming()].map(f64::to_bits)
}
fn render(mut performance: ReedPerformance, block: usize) -> (Vec<f64>, [u64;4]) {
    let mut pressure = vec![0.0;performance.info().samples as usize];
    for chunk in pressure.chunks_mut(block) { performance.pressure_block(chunk).unwrap(); }
    assert_eq!(performance.remaining_samples(), 0);
    assert!(pressure.iter().all(|p| p.is_finite()));
    (pressure, motion(&performance))
}

#[test]
fn authored_zero_inlet_species_and_column_permutations_preserve_original_sources() {
    let (a,b) = (packs("A",0.028,3.5,200.0), packs("B",0.044,4.5,200.0));
    let bindings = [a.binding(),b.binding()];
    let sources = ThermochemicalSources { species: &bindings, water_species: None, phase: ThermochemicalPhase::FrozenGas };
    let (field,records) = parse(BLOCK,&sources);
    assert_eq!(field.stations_m(), &[0.0,0.285,0.57]);
    assert_eq!(field.primitive_states()[0].composition.fractions(), &[1.0,0.0]);
    assert_eq!(field.model().source_packs().len(),2);
    let inlet = field.state_at(0.0).unwrap();
    assert_eq!(inlet.components().len(),1);
    assert_eq!(inlet.source_packs().len(),1);
    assert_eq!(inlet.source_packs()[0].species_pack,a.species.content_hash());
    assert_eq!(inlet.thermodynamics().properties().molar_mass().value(),0.028);
    let middle = field.state_at(0.1425).unwrap();
    assert_eq!(middle.components().len(),2);
    assert_eq!(middle.state().temperature,310.0);
    assert!((middle.thermodynamics().properties().molar_mass().value()-0.030).abs()<1e-16);
    assert_eq!(records.len(),3);
    for (record,position) in records.iter().zip(field.stations_m()) {
        assert_eq!(record,&GasRecord::Thermochemical(field.state_at(*position).unwrap()));
    }
    let permuted = BLOCK.replace("mole 2 A B","mole 2 B A")
        .replace("300 1 0","300 0 1").replace("320 0.75 0.25","320 0.25 0.75");
    let (other,other_records) = parse(&permuted,&sources);
    assert_eq!(field,other);
    assert_eq!(records,other_records);
}

#[test]
fn catalogue_seed_never_reaches_physical_samples_or_the_actual_pressure_path() {
    let (a,b) = (packs("A",0.028,3.5,200.0), packs("B",0.044,4.5,200.0));
    let bindings = [a.binding(),b.binding()];
    let sources = ThermochemicalSources { species: &bindings, water_species: None, phase: ThermochemicalPhase::FrozenGas };
    let (field,_) = parse(BLOCK,&sources);
    let different_seed = from_packs(&bindings,
        Composition::new(CompositionBasis::MoleFraction,vec![0.875,0.125]).unwrap(),None).unwrap();
    assert_ne!(field.model().thermodynamics().mole_composition(),
        different_seed.thermodynamics().mole_composition());
    let direct = AxialGasField::new(field.stations_m().to_vec(),100000.0,
        field.primitive_states().to_vec(),different_seed,QueryPoint::new(),
        ThermochemicalPhase::FrozenGas).unwrap();
    for x in [0.0,0.07125,0.1425,0.285,0.4275,0.57] {
        assert_eq!(field.state_at(x).unwrap(),direct.state_at(x).unwrap());
    }
    for (mode,resolution,segments) in [
        ("8",GasFieldResolution::Midpoint(8),8),
        ("midpoint 8",GasFieldResolution::Midpoint(8),8),
        ("continuous 16",GasFieldResolution::Continuous(fs_duct::RiccatiSteps(16)),1),
    ] {
        let text = with_block(&BLOCK.replace("gas_field 8 ",&format!("gas_field {mode} ")));
        let file = load(&text,&sources,37);
        let inlet = field.state_at(0.0).unwrap();
        assert_eq!(file.gas_field_resolution(),Some(resolution));
        assert_eq!(voice(&file).rho,inlet.state().density);
        assert_eq!(voice(&file).inlet_gas.sound_speed,inlet.state().sound_speed);
        assert_eq!(file.segment_gases().len(),segments);
        assert_eq!(file.outlet_receiver().unwrap().gas,*field.state_at(0.57).unwrap().state());
        assert_eq!(file.gas_records().len(),4,"ambient plus three authored field knots");
        let injected = ReedPerformance::from_bytes_with_field_resolution(
            base().as_bytes(),37,ObservationPolicy::Diagnostic,Some(&sources),direct.clone(),resolution).unwrap();
        assert_eq!(render(file,37),render(injected,37));
        assert_eq!(render(load(&text,&sources,37),37),
            render(load(&text,&sources,113),113));
    }
}

#[test]
fn turning_added_gas_off_is_a_data_change_and_downstream_gas_changes_reed_motion() {
    let (a,b) = (packs("A",0.028,3.5,200.0), packs("B",0.044,4.5,200.0));
    let bindings = [a.binding(),b.binding()];
    let sources = ThermochemicalSources { species: &bindings, water_species: None, phase: ThermochemicalPhase::FrozenGas };
    let off = BLOCK.replace("320 0.75 0.25","300 1 0").replace("310 0.5 0.5","300 1 0");
    let pure = off.replace("mole 2 A B","mole 1 A").replace("300 1 0","300 1");
    let fixed_exterior = |block| with_block(block).replace("termination flanged","exterior 300 100000 0\ntermination flanged");
    let changed_text = fixed_exterior(BLOCK);
    let unchanged_text = fixed_exterior(&off);
    let mut changed = load(&changed_text,&sources,37);
    let mut unchanged = load(&unchanged_text,&sources,37);
    assert_eq!(voice(&changed).rho,voice(&unchanged).rho);
    assert_eq!(voice(&changed).zc,voice(&unchanged).zc);
    assert_eq!(changed.outlet_receiver().unwrap().gas,unchanged.outlet_receiver().unwrap().gas);
    let (mut energy,mut difference) = (0.0,0.0);
    let mut changed_motion = false;
    for _ in 0..4801 {
        let (mut p,mut q) = ([0.0],[0.0]);
        changed.pressure_block(&mut p).unwrap(); unchanged.pressure_block(&mut q).unwrap();
        energy += q[0]*q[0]; difference += (p[0]-q[0])*(p[0]-q[0]);
        changed_motion |= voice(&changed).reed_y != voice(&unchanged).reed_y
            || voice(&changed).reed_v != voice(&unchanged).reed_v;
    }
    assert!(changed_motion);
    assert!(energy>0.0 && difference>1e-8*energy);
    assert_eq!(render(load(&unchanged_text,&sources,37),37),
        render(load(&fixed_exterior(&pure),&sources,37),37));
}

#[test]
fn inactive_species_does_not_query_its_out_of_domain_source_or_invent_water() {
    let (a,b) = (packs("A",0.028,3.5,200.0), packs("B",0.044,4.5,400.0));
    let bindings = [a.binding(),b.binding()];
    let sources = ThermochemicalSources { species: &bindings, water_species: None, phase: ThermochemicalPhase::FrozenGas };
    let off = BLOCK.replace("320 0.75 0.25","320 1 0").replace("310 0.5 0.5","310 1 0");
    let (field,_) = parse(&off,&sources);
    for x in [0.0,0.1425,0.285,0.57] {
        assert_eq!(field.state_at(x).unwrap().components().len(),1);
    }
    assert_eq!(field.model().source_packs().len(),2);
    let (a,b) = (packs("A",0.028,3.5,200.0), packs("B",0.018,4.0,200.0));
    let bindings = [a.binding(),b.binding()];
    let water = b.species.association().species();
    let sources = ThermochemicalSources { species: &bindings, water_species: Some(water),
        phase: ThermochemicalPhase::LiquidWaterLimit };
    let dry = "gas_field 8 2 mole 2 A B\ngas_station 0 300 1 0\ngas_station 0.57 300 1 0\n";
    let (field,_) = parse(dry,&sources);
    assert_eq!(field.state_at(0.0).unwrap().state().water_mole_fraction,0.0);
    let wet = dry.replace("0.57 300 1 0","0.57 300 0.5 0.5");
    let mut reader = Reader { lines: wet.lines(), line: 0 };
    assert!(read(&mut reader,Some(&sources),100000.0,&mut Vec::new()).is_err());
}

#[test]
fn field_admission_refuses_conflicting_sources_domains_dimensions_and_work() {
    let (a,b) = (packs("A",0.028,3.5,200.0), packs("B",0.044,4.5,200.0));
    let bindings = [a.binding(),b.binding()];
    let sources = ThermochemicalSources { species: &bindings, water_species: None, phase: ThermochemicalPhase::FrozenGas };
    for text in [
        BLOCK.replace("8 3 mole","0 3 mole"), BLOCK.replace("8 3 mole","65 3 mole"),
        BLOCK.replace("8 3 mole","midpoint 0 3 mole"),
        BLOCK.replace("8 3 mole","midpoint 65 3 mole"),
        BLOCK.replace("8 3 mole","continuous 0 3 mole"),
        BLOCK.replace("8 3 mole","continuous 65537 3 mole"),
        BLOCK.replace("8 3 mole","unknown 8 3 mole"),
        BLOCK.replace("8 3 mole","8 1 mole"), BLOCK.replace("8 3 mole","8 4097 mole"),
        BLOCK.replace("mole 2","mass 2"), BLOCK.replace("mole 2","mole 0"),
        BLOCK.replace("mole 2","mole 129"), BLOCK.replace("8 3 mole 2","8 4096 mole 128"),
        BLOCK.replace("A B","A A"), BLOCK.replace("A B","A missing"),
        BLOCK.replace("gas_station 0 300","gas_station 0.1 300"),
        BLOCK.replace("gas_station 0.285","gas_station 0"),
        BLOCK.replace("300 1 0","300 1 0.1"),
        BLOCK.replace("300 1 0","300 1.1 -0.1"),
        BLOCK.replace("300 1 0","NaN 1 0"),
        BLOCK.replace("320 0.75 0.25","199 0.75 0.25"),
    ] {
        let mut reader = Reader { lines: text.lines(), line: 0 };
        assert!(read(&mut reader,Some(&sources),100000.0,&mut Vec::new()).is_err(),"{text}");
    }
    let mut reader = Reader { lines: BLOCK.lines(), line: 0 };
    assert!(read(&mut reader,None,100000.0,&mut Vec::new()).is_err());
    let duplicate = [a.binding(),a.binding(),b.binding()];
    let ambiguous = ThermochemicalSources { species: &duplicate, water_species: None, phase: ThermochemicalPhase::FrozenGas };
    let mut reader = Reader { lines: BLOCK.lines(), line: 0 };
    assert!(read(&mut reader,Some(&ambiguous),100000.0,&mut Vec::new()).is_err());
    let missing_transport = [a.binding(),SpeciesPackRef { transport: None, ..b.binding() }];
    let incomplete = ThermochemicalSources { species: &missing_transport, water_species: None, phase: ThermochemicalPhase::FrozenGas };
    let mut reader = Reader { lines: BLOCK.lines(), line: 0 };
    assert!(read(&mut reader,Some(&incomplete),100000.0,&mut Vec::new()).is_err());
    for text in [
        with_block(BLOCK).replace("cylinder 0.0015 0.57","cylinder 0.0015 0.57 gas 300 0"),
        with_block(BLOCK).replace("cylinder 0.0015 0.57","cylinder 0.0015 0.56"),
        with_block(BLOCK).replace("compile_limits", &format!("{BLOCK}compile_limits")),
    ] {
        assert!(ReedPerformance::from_bytes_with_gas_sources(
            text.as_bytes(),37,ObservationPolicy::Diagnostic,Some(&sources)).is_err());
    }
    let (field,_) = parse(BLOCK,&sources);
    assert!(ReedPerformance::from_bytes_with_gas_field(
        with_block(BLOCK).as_bytes(),37,ObservationPolicy::Diagnostic,Some(&sources),field,8).is_err());
}

#[test]
fn absent_field_keeps_legacy_parser_and_pressure_bytes_unchanged() {
    let a = packs("A",0.028,3.5,200.0);
    let bindings = [a.binding()];
    let sources = ThermochemicalSources { species: &bindings, water_species: None, phase: ThermochemicalPhase::FrozenGas };
    let original = ReedPerformance::from_bytes(base().as_bytes(),37).unwrap();
    let admitted = load(&base(),&sources,113);
    assert!(original.gas_field().is_none());
    assert!(admitted.gas_field().is_none());
    assert_eq!(original.gas_records(),admitted.gas_records());
    assert_eq!(render(original,37),render(admitted,113));
}

#[test]
fn checked_in_field_example_reaches_the_shared_pressure_owner() {
    // These are explicit synthetic packs, not source claims for the named gases.
    let (a,b) = (packs("N2",0.028,3.5,200.0), packs("CO2",0.044,4.5,200.0));
    let bindings = [a.binding(),b.binding()];
    let sources = ThermochemicalSources { species: &bindings, water_species: None, phase: ThermochemicalPhase::FrozenGas };
    let input = include_str!("../../../../../examples/reed-axial-gas-field.performance");
    let performance = load(input,&sources,113);
    assert_eq!(performance.gas_field_resolution(),Some(
        GasFieldResolution::Continuous(fs_duct::RiccatiSteps(16))));
    let (field,steps) = performance.gas_field().unwrap();
    assert_eq!(steps,16);
    assert_eq!(field.state_at(0.0).unwrap().components().len(),1);
    assert_eq!(field.state_at(0.57).unwrap().components().len(),2);
    assert_eq!(field.model().source_packs().len(),2);
    let (pressure,_) = render(performance,113);
    assert!(pressure.iter().any(|p| *p!=0.0));
}


#[test]
fn primitive_field_drives_actual_plate_dg_inlet_slit_and_causal_pressure() {
    use crate::bernoulli_aperture::performance::{CoupledAperture, file::PlateValvePerformance};
    use fs_exec::CancelGate;
    let (a, b) = (packs("A", 0.028, 3.5, 200.0), packs("B", 0.044, 4.5, 200.0));
    let bindings = [a.binding(), b.binding()];
    let sources = ThermochemicalSources {
        species: &bindings, water_species: None, phase: ThermochemicalPhase::FrozenGas,
    };
    let example = include_str!("../../../../../examples/plate-valve-dg-bore.performance");
    let header = example.split_once("\nschedule\n").unwrap().0
        .replace("audio 96000 9602", "audio 96000 256")
        .replace("1e-11 224 4", "1e-11 224 4 gas_field 3 mole 2 A B")
        .replace("observation inlet", "gas_station 0 300 1 0\ngas_station 0.025 320 0.5 0.5\ngas_station 0.05 310 0 1\nobservation inlet");
    let input = format!("{header}\nschedule\nfrankensim-gesture-schedule-v1\ncontrol_rate_hz\t700\ntracks\t1\ntrack\tmouth\tblowing-pressure\ninitial\tpressure-pa\t5e0\nevents\t0\n");
    let render = |text: &str, block| {
        let source = PlateValvePerformance::from_bytes_with_gas_sources(
            text.as_bytes(), block, &CancelGate::new_clock_free(), Some(&sources),
        ).unwrap();
        let CoupledAperture::Dg(dg) = source.renderer().system() else { panic!("original DG") };
        let duct = dg.load().duct();
        let field = duct.gas_field().unwrap();
        assert_eq!(field.primitive_states()[0].composition.fractions(), &[1.0, 0.0]);
        assert_eq!(field.primitive_states()[2].composition.fractions(), &[0.0, 1.0]);
        assert_eq!(duct.inlet_gas().temperature, 300.0);
        assert_eq!(duct.outlet_gas().temperature, 310.0);
        assert_ne!(duct.inlet_gas(), source.ambient_gas());
        assert_eq!(dg.aperture().viscous_slit().unwrap().dynamic_viscosity_pa_s(),
            duct.inlet_gas().dynamic_viscosity);
        assert!(duct.spatial_viscothermal_admission().is_some());
        let mut renderer = source.into_renderer();
        let mut pressure = vec![0.0; 256];
        for chunk in pressure.chunks_mut(block) { renderer.block(chunk).unwrap(); }
        pressure
    };
    let reference = render(&input, 37);
    assert_eq!(render(&input, 1), reference);
    assert!(reference.iter().any(|p| p.abs() > 1e-6));
    // Inlet, outlet and solid data stay fixed: the interior primitive change
    // must reach the original acoustic/wall operator and actual coupled trace.
    let changed = input.replace("gas_station 0.025 320 0.5 0.5",
        "gas_station 0.025 350 0.9 0.1");
    assert!(render(&changed, 37).iter().zip(&reference)
        .any(|(a, b)| (a - b).abs() > 1e-8));
}
