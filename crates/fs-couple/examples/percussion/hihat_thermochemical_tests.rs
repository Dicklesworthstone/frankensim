//! Actual species/NASA-9/transport sources reach native percussion mechanics.
use super::*;
use crate::{ambient::Ambient, drum_spec, thermochemical};
use fs_material::gas_species::input::{GasRecord, ThermochemicalPhase, ThermochemicalSources};
use fs_material::thermochemical_gas::packs::SpeciesPackRef;

#[path = "../../src/render/schedule/reed/field_input/tests/support.rs"]
mod support;

const RECORD: &str = "thermochemical 293.15 101325 mole 1 percussion-atmosphere 1";
const HUMID: &str = "293.15 101325 0.65";

fn registry<'a>(species: &'a [SpeciesPackRef<'a>]) -> ThermochemicalSources<'a> {
    ThermochemicalSources { species, water_species: None, phase: ThermochemicalPhase::FrozenGas }
}

fn film() -> String {
    format!("frankensim-squeeze-film-v1\nannulus,0.04,0.085,2,4\n\
        limits,0.000001,0.04,200000\nboundary,sealed,sealed\ngas,{RECORD}\n")
}

#[test]
fn native_source_loading_reaches_drum_profile_mesh_nested_pair_and_override() {
    let packs = support::packs("percussion-atmosphere", 0.028, 3.5, 200.0);
    let args = ["hihat", "pair.fshh", "64", "--gas-component", "species.pack", "nasa9.pack",
        "transport.pack", "--gas-water", "none", "--gas-phase", "frozen"].map(String::from);
    let (remaining, options) = thermochemical::Options::extract(&args).unwrap();
    assert_eq!(remaining, ["hihat", "pair.fshh", "64"]);
    let loaded = options.load_with(|path| match path.to_str().unwrap() {
        "species.pack" => Ok(packs.species.to_bytes()),
        "nasa9.pack" => Ok(packs.nasa9.to_bytes()),
        "transport.pack" => Ok(packs.transport.to_bytes()),
        _ => Err("unexpected source path".into()),
    }).unwrap();
    let bindings = loaded.bindings(); let sources = loaded.registry(&bindings).unwrap();
    let declared = Ambient::described_with_gas_sources(RECORD, Some(&sources)).unwrap();
    let GasRecord::Thermochemical(source) = declared.record() else { panic!("lost original gas source"); };
    assert_eq!(source.source_packs()[0].species_pack, packs.species.content_hash());
    assert_eq!(source.source_packs()[0].nasa9_pack, packs.nasa9.content_hash());
    assert_eq!(source.source_packs()[0].transport_pack, packs.transport.content_hash());
    assert!(declared.receipt().contains("thermochemical_gases"));

    let dir = std::env::temp_dir().join(format!("fs-percussion-thermochemical-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let profile = include_str!("estimated-splash.profile").replace(HUMID, RECORD);
    let mut shell = specimen::Specimen::parse_with_gas_sources(&profile, Some(&sources)).unwrap();
    shell.azimuths = 8;
    let mesh = format!("{}ambient,{RECORD}\n", shell.mesh_text().unwrap());
    std::fs::write(dir.join("upper.profile"), &profile).unwrap();
    std::fs::write(dir.join("lower.mesh"), &mesh).unwrap();
    std::fs::write(dir.join("ambient.gas"), RECORD).unwrap();
    let drum = include_str!("estimated_drum.fsd").replace(HUMID, RECORD);
    std::fs::write(dir.join("drum.fsd"), drum).unwrap();
    let pair = include_str!("estimated-hihat.fshh").replace(HUMID, RECORD)
        .replace("shell,upper,profile,estimated-splash.profile", "shell,upper,profile,upper.profile")
        .replace("shell,lower,profile,estimated-splash.profile", "shell,lower,mesh,lower.mesh");
    std::fs::write(dir.join("pair.fshh"), pair).unwrap();
    let drum = drum_spec::Spec::load_with_gas_sources(dir.join("drum.fsd").to_str().unwrap(), Some(&sources)).unwrap();
    assert_eq!(drum.ambient, declared);
    for (profile, mesh) in [(Some(dir.join("upper.profile")), None), (None, Some(dir.join("lower.mesh")))] {
        let shell = specimen::load_selection_with_gas_sources(profile.as_deref(), mesh.as_deref(), Some(&sources)).unwrap().unwrap();
        assert_eq!(shell.ambient.as_ref(), Some(&declared));
    }
    let pair = input::Spec::load_with_gas_sources(&dir.join("pair.fshh"), Some(&sources)).unwrap();
    assert_eq!(pair.ambient, declared);
    for shell in pair.shells_with_gas_sources(Some(&sources)).unwrap() {
        assert_eq!(shell.ambient.as_ref(), Some(&declared));
    }
    assert!(pair.shells().is_err(), "nested sources must never become reference air");
    let override_gas = crate::ambient::load_override(Some(&dir.join("ambient.gas")), Some(&sources)).unwrap().unwrap();
    assert_eq!(override_gas.record(), declared.record());
    assert_eq!(override_gas.source, crate::ambient::Source::CommandLine);
    assert_eq!(crate::ambient::resolve(&Ambient::reference(), Some(&override_gas)), override_gas);
}

#[test]
fn nasa9_source_changes_real_drum_pressure_feedback_without_retuning_heads() {
    let original = support::packs("percussion-atmosphere", 0.028, 3.5, 200.0);
    let changed = support::packs("percussion-atmosphere", 0.028, 4.5, 200.0);
    let bindings = [original.binding()]; let replacement = [changed.binding()];
    let a = Ambient::described_with_gas_sources(RECORD, Some(&registry(&bindings))).unwrap();
    let b = Ambient::described_with_gas_sources(RECORD, Some(&registry(&replacement))).unwrap();
    assert_eq!(a.gas().density, b.gas().density);
    assert_ne!(a.gas().sound_speed, b.gas().sound_speed);
    let build = |ambient| crate::drum_with_spec(128, 2e-6, false, false, None, false,
        Stroke { speed_m_s: 4.0, position_m: Some([0.06, 0.01]) }, false, None,
        Some(drum_spec::Spec { ambient, radial_intervals: 2, azimuths: 16, ..drum_spec::Spec::reference() })).unwrap();
    let mut first = build(a.clone()); let mut second = build(b); let mut replay = build(a.clone());
    assert_eq!(first.ambient.record(), a.record());
    assert_eq!(first.observer_a, second.observer_a);
    assert_eq!(first.pressure.as_ref().unwrap().areas, second.pressure.as_ref().unwrap().areas);
    assert_eq!(first.system.state(), second.system.state());
    // Independent NASA-9 polynomial in the supplied synthetic pack: cp/R = a2 + a3*T.
    let cp_r = 3.5 + 0.001 * 293.15;
    let expected_bulk = 101325.0 * cp_r / (cp_r - 1.0);
    assert!((first.pressure.as_ref().unwrap().bulk_modulus_pa / expected_bulk - 1.0).abs() < 1e-12);
    let gate = CancelGate::new_clock_free(); let mut changed_motion = 0.0_f64;
    for _ in 0..128 {
        for e in [&mut first, &mut second, &mut replay] {
            let frame = e.system.step(&e.force, &gate).unwrap();
            assert!(frame.balance_residual_j.abs() < 2e-7);
        }
        assert_eq!(first.system.state(), replay.system.state());
        for (x, y) in first.system.state().iter().zip(second.system.state()) {
            changed_motion = changed_motion.max((x - y).abs());
        }
    }
    assert!(changed_motion > 1e-16, "source heat capacity must change accepted mechanical feedback");
}

#[test]
fn thermochemical_film_keeps_original_sources_and_species_mass_on_real_skins() {
    let light = support::packs("percussion-atmosphere", 0.028, 3.5, 200.0);
    let heavy = support::packs("percussion-atmosphere", 0.044, 3.5, 200.0);
    let mut shell = specimen::Specimen::reference(); shell.azimuths = 8;
    let mut initial_masses = Vec::new();
    for packs in [&light, &heavy] {
        let bindings = [packs.binding()]; let sources = registry(&bindings);
        let config = squeeze::Config::parse_with_gas_sources(&film(), Some(&sources)).unwrap();
        let source = config.material_gas().unwrap();
        let spec = input::Spec::parse(include_str!("estimated-hihat.fshh")).unwrap();
        let pair = build_with_squeeze(&spec, &shell, &shell,
            Stroke { speed_m_s: 0.0, position_m: None }, None, 128, 2e-6, false, Some(&config)).unwrap();
        let coordinate = pair.upper_modes.start; let mut e = pair.experiment;
        let Mechanics::Reference(system) = &e.system else { panic!("original impact owner required"); };
        assert_eq!(system.gas_film_source(), Some(source));
        let initial = system.gas_film_observation().unwrap().unwrap();
        assert_eq!(initial.minimum_pressure_pa, source.state().pressure);
        initial_masses.push(initial.mass_kg);
        e.system = e.system.into_analytic_nonlinear().unwrap().with_impact_substeps(
            fs_couple::render::plate::impact::ImpactSubstepConfig { max_depth: 8, max_attempts: 511 }).unwrap();
        let mut force = vec![0.0; e.force.len()]; force[coordinate] = -5.0;
        let gate = CancelGate::new_clock_free(); let mut peak = initial.maximum_pressure_pa;
        for _ in 0..128 {
            e.system.step(&force, &gate).unwrap();
            let gas = play::gas_observation(&e.system).unwrap().unwrap();
            assert!((gas.mass_kg / initial.mass_kg - 1.0).abs() < 1e-8);
            peak = peak.max(gas.maximum_pressure_pa);
        }
        assert!(peak > initial.maximum_pressure_pa + 0.001, "actual skin closure must compress source-backed gas");
    }
    assert!((initial_masses[1] / initial_masses[0] - 0.044 / 0.028).abs() < 1e-12,
        "same physical pocket at fixed p,T must contain the source species mass pV/(RT)");
}

#[test]
fn missing_ambiguous_or_out_of_domain_sources_refuse_and_legacy_air_is_unchanged() {
    let packs = support::packs("percussion-atmosphere", 0.028, 3.5, 200.0);
    let late = support::packs("percussion-atmosphere", 0.028, 3.5, 350.0);
    let mut no_nasa9 = packs.binding(); no_nasa9.nasa9 = None;
    let mut no_transport = packs.binding(); no_transport.transport = None;
    let drum = include_str!("estimated_drum.fsd").replace(HUMID, RECORD);
    let profile = include_str!("estimated-splash.profile").replace(HUMID, RECORD);
    let pair = include_str!("estimated-hihat.fshh").replace(HUMID, RECORD);
    assert!(Ambient::described(RECORD).is_err());
    for bindings in [vec![], vec![no_nasa9], vec![no_transport], vec![late.binding()],
        vec![packs.binding(), packs.binding()]] {
        let sources = registry(&bindings);
        assert!(drum_spec::Spec::read_with_gas_sources(&drum, Some(&sources)).is_err());
        assert!(specimen::Specimen::parse_with_gas_sources(&profile, Some(&sources)).is_err());
        assert!(input::Spec::parse_with_gas_sources(&pair, Some(&sources)).is_err());
        assert!(squeeze::Config::parse_with_gas_sources(&film(), Some(&sources)).is_err());
    }
    let bindings = [packs.binding()]; let sources = registry(&bindings);
    for record in [HUMID, "mixture 293.15 101325 mole 2 N2 0.8 O2 0.2"] {
        assert_eq!(Ambient::described(record).unwrap(),
            Ambient::described_with_gas_sources(record, Some(&sources)).unwrap());
    }
    for bad in [format!("{RECORD} extra"), RECORD.replace("293.15", "1100"), RECORD.replace(" 1 percussion", " 0 percussion")] {
        assert!(Ambient::described_with_gas_sources(&bad, Some(&sources)).is_err());
    }
}
