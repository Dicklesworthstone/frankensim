//! Synthetic source-backed membranes exercise the shared owners, not measured instruments.
use super::*;
use crate::driving_point::regional::RegionalDuct;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialStateId, NormalizedMaterialCardPack, NormalizedPack,
    ObservationDataset, PropertyClaim, PropertyKey, PropertyValue, Provenance, UncertaintyModel,
};
use fs_qty::{Density, Dims, DynViscosity, Pressure};

const FIXTURE: &str = include_str!("../../../../examples/reed-baffled-mixture.performance");
const WALL: &str =
    "bore_wall isotropic wall.matpack 300 1000 0.00001 stress-free 0.000001 5000 kelvin-voigt\n";
const RATE: u32 = 8000;
const SAMPLES: usize = 4801;

fn input(wall: &str) -> String {
    FIXTURE
        .replace("audio 48000 4801 2", "audio 8000 4801 2")
        .replace("baffled_outlet 0 0 0.5 4 16 1000", "baffled_outlet 0 0 0.5 4 16 800")
        .replace("observer_law stokes-kirchhoff", "observer_law none")
        .replace("reed 0.0004 0.013 6000 0 0 0.35", "reed 0.0004 0.013 6000 0.00001 100 0.35")
        .replace("compile_limits", &format!("{wall}compile_limits"))
}

fn wall_pack(young: f64, omit: Option<&str>) -> NormalizedMaterialCardPack {
    let artifact = hash_domain("org.frankensim.test.scheduled-reed-wall", &young.to_le_bytes());
    let provenance = Provenance {
        source: "synthetic local membrane; not measured".into(),
        license: "CC0-1.0".into(),
        artifact: Some(artifact),
    };
    let mut claims = ClaimSet::new();
    let observation = claims.register_observation(ObservationDataset {
        specimen: "synthetic wall".into(),
        method: "manufactured integration regression".into(),
        artifact,
        caveats: "not a measured material or instrument".into(),
        provenance: provenance.clone(),
    }).unwrap();
    for (name, dims, value) in [
        ("density", Density::DIMS, 1000.0),
        ("young_modulus", Pressure::DIMS, young),
        ("poisson_ratio", Dims::NONE, 0.3),
        ("kelvin_voigt_young_viscosity", DynViscosity::DIMS, 10.0),
    ] {
        if omit == Some(name) { continue; }
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(name, dims),
            value: PropertyValue::Scalar { value, dims },
            validity: ValidityDomain::unconstrained()
                .with("T", 200.0, 800.0).with("omega", 0.0, 1_000_000.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![observation.clone()],
            provenance: provenance.clone(),
        }).unwrap();
    }
    NormalizedMaterialCardPack::new(
        MaterialStateId {
            chemistry: "synthetic".into(), phase: "solid".into(),
            process: "manufactured regression".into(), revision: 0,
        },
        NormalizedPack::new("synthetic-wall", "authored-v1", artifact, "CC0-1.0",
            claims, vec![], vec![]).unwrap(),
    ).unwrap()
}

fn load(text: &str, pack: &NormalizedMaterialCardPack, block: usize)
    -> Result<ReedPerformance, ReedPerformanceError>
{
    ReedPerformance::from_bytes_with_sources(text.as_bytes(), block,
        ObservationPolicy::Diagnostic, None, |path| {
            assert_eq!(path, "wall.matpack");
            Ok(pack.clone())
        })
}

fn voice(performance: &ReedPerformance) -> &ReedBoreVoice {
    let RenderVoice::ReedBore(voice) = &performance.renderer.context().voices[0] else {
        panic!("expected original reed voice")
    };
    voice
}

fn motion(voice: &ReedBoreVoice) -> [u64; 5] {
    [voice.reed_y, voice.reed_v, voice.p_plus_prev, voice.f_jet_prev, voice.line.incoming()]
        .map(f64::to_bits)
}

fn rendered(mut performance: ReedPerformance, block: usize) -> (Vec<u64>, [u64; 5]) {
    let mut pressure = vec![0.0; performance.info().samples as usize];
    for chunk in pressure.chunks_mut(block) { performance.pressure_block(chunk).unwrap(); }
    assert!(pressure.iter().all(|p| p.is_finite()));
    assert_eq!(performance.remaining_samples(), 0);
    (pressure.into_iter().map(f64::to_bits).collect(), motion(voice(&performance)))
}

fn reed() -> BeatingReed {
    BeatingReed {
        rest_opening_m: 0.0004, width_m: 0.013, closing_pressure_pa: 6000.0,
        mass_kg: 0.00001, stiffness_n_m: 100.0, damping_ratio: 0.35,
        blowing_pressure_pa: 0.0, attack_s: 0.0,
    }
}

fn geometry() -> Duct {
    Duct { segments: vec![
        Segment::Cylinder { radius: 0.0015, length: 0.2 },
        Segment::Cylinder { radius: 0.0015, length: 0.3 },
    ] }
}

#[test]
fn omitted_and_explicit_rigid_walls_preserve_pressure_and_reed_history() {
    for legacy in [false, true] {
        let make = |wall| {
            let text = input(wall);
            if legacy {
                text.replace("baffled_outlet 0 0 0.5 4 16 800", "listener_m 1")
            } else { text }
        };
        let plain = ReedPerformance::from_bytes(make("").as_bytes(), 37).unwrap();
        let rigid = ReedPerformance::from_bytes_with_sources(
            make("bore_wall rigid\n").as_bytes(), 113, ObservationPolicy::Diagnostic, None,
            |_| panic!("rigid wall must not load a material source"),
        ).unwrap();
        assert!(plain.bore_wall().is_none() && rigid.bore_wall().is_none());
        assert_eq!(rendered(plain, 37), rendered(rigid, 113));
    }
}

#[test]
fn one_source_drives_scheduled_feedback_and_the_matched_physical_microphone() {
    let pack = wall_pack(6e7, None);
    let text = input(WALL);
    let mut actual = load(&text, &pack, 37).unwrap();
    let record = actual.bore_wall().unwrap().clone();
    assert_eq!(record.path(), "wall.matpack");
    assert_eq!(record.pack(), &pack);
    assert_eq!(record.wall().material().state_coordinate("T"), Some(300.0));
    assert_ne!(actual.ambient_gas().temperature, 300.0);
    let info = actual.outlet_receiver().unwrap();
    let regional = RegionalDuct::new(&geometry(), actual.segment_gases()).unwrap()
        .with_exterior_gas(info.gas).unwrap();
    let wall = fs_duct::WallModel::Law(record.wall());
    let mut direct = ReedBoreVoice::new_with_regional_wall(&regional, reed(),
        Termination::FlangedOpen, PlateBank::default(), 0.5, RATE, SAMPLES, wall).unwrap();
    let mut microphone = regional.baffled_outlet_with_wall(
        RATE, info.fft_samples, info.receiver, wall).unwrap();
    microphone.step(direct.p_plus_prev).unwrap();
    let schedule = GestureSchedule::from_canonical_bytes(
        text.split_once("\nschedule\n").unwrap().1.as_bytes()).unwrap();
    let mut changed = load(&text, &wall_pack(2e8, None), 37).unwrap();
    assert_eq!(actual.segment_gases(), changed.segment_gases());
    assert_eq!(voice(&actual).zc, voice(&changed).zc);
    let (mut mechanics_changed, mut pressure_changed, mut audible) = (false, false, false);
    for n in 0..SAMPLES {
        direct.set_blowing_pressure(schedule.sample("blow", n as u64 * 700 / u64::from(RATE)).unwrap());
        direct.step_block(&mut [0.0]).unwrap();
        let (mut sample, mut other) = ([0.0], [0.0]);
        actual.pressure_block(&mut sample).unwrap();
        changed.pressure_block(&mut other).unwrap();
        assert_eq!(motion(voice(&actual)), motion(&direct));
        assert_eq!(sample[0].to_bits(), microphone.step(direct.p_plus_prev).unwrap().to_bits());
        assert!(sample[0].is_finite() && other[0].is_finite());
        mechanics_changed |= voice(&actual).reed_y.to_bits() != voice(&changed).reed_y.to_bits();
        pressure_changed |= sample[0].to_bits() != other[0].to_bits();
        audible |= sample[0].abs() > 1e-9;
    }
    assert!(mechanics_changed && pressure_changed && audible);
    assert_eq!(rendered(load(&text, &pack, 37).unwrap(), 37),
        rendered(load(&text, &pack, 113).unwrap(), 113));
}

#[test]
fn material_wall_reaches_the_register_shunt_and_each_physical_mouth() {
    let pack = wall_pack(6e7, None);
    let text = input(WALL).replace("segments 2\ncylinder 0.0015 0.2",
        "segments 4\ncylinder 0.0015 0.1\nhole 0.0005 0.002 0.0015 0.3\ncylinder 0.0015 0.1");
    let mut actual = load(&text, &pack, 37).unwrap();
    let record = actual.bore_wall().unwrap().clone();
    let duct = Duct { segments: vec![
        Segment::Cylinder { radius: 0.0015, length: 0.1 },
        Segment::ToneHole { hole_radius: 0.0005, chimney_height: 0.002,
            bore_radius: 0.0015, state: HoleState::Vent(0.3) },
        Segment::Cylinder { radius: 0.0015, length: 0.1 },
        Segment::Cylinder { radius: 0.0015, length: 0.3 },
    ] };
    let info = actual.outlet_receiver().unwrap();
    let regional = RegionalDuct::new(&duct, actual.segment_gases()).unwrap()
        .with_exterior_gas(info.gas).unwrap();
    let wall = fs_duct::WallModel::Law(record.wall());
    let mut direct = ReedBoreVoice::new_with_regional_wall(&regional, reed(),
        Termination::FlangedOpen, PlateBank::default(), 0.5, RATE, SAMPLES, wall).unwrap();
    let grid = direct.line.history().len();
    let (mut bell, mut holes, single) = regional
        .baffled_outlet_and_open_holes_and_shunt_with_delay_speeds(
            RATE, grid, info.receiver, wall, info.gas.sound_speed,
            &[info.gas.sound_speed], 1.0).unwrap();
    assert_eq!(holes.len(), 1);
    assert_eq!(actual.open_hole_observers().count(), 1);
    let mut shunt = regional.select_register_hole_shunt(single, RATE, grid, wall)
        .unwrap().expect("the vented chimney has an original nonlinear shunt");
    shunt.prime(direct.p_plus_prev).unwrap();
    direct.set_register_shunt(shunt);
    bell.step(direct.p_plus_prev).unwrap();
    holes[0].step(direct.p_plus_prev).unwrap();
    let schedule = GestureSchedule::from_canonical_bytes(
        text.split_once("\nschedule\n").unwrap().1.as_bytes()).unwrap();
    let mut hole_audible = false;
    for n in 0..SAMPLES {
        direct.set_blowing_pressure(schedule.sample("blow", n as u64 * 700 / u64::from(RATE)).unwrap());
        direct.step_block(&mut [0.0]).unwrap();
        let mut sample = [0.0];
        actual.pressure_block(&mut sample).unwrap();
        assert_eq!(motion(voice(&actual)), motion(&direct));
        let bell_pa = bell.step(direct.p_plus_prev).unwrap();
        let hole_pa = holes[0].step(direct.p_plus_prev).unwrap();
        assert_eq!(sample[0].to_bits(), (bell_pa + hole_pa).to_bits());
        hole_audible |= hole_pa.abs() > 1e-9;
    }
    assert!(hole_audible);
    assert!(actual.last_outlet_error().is_none());
}

#[test]
fn a_continuous_gas_field_and_held_wall_share_the_same_feedback_and_outlet() {
    use crate::driving_point::regional::field::tests::field;
    let pack = wall_pack(6e7, None);
    let text = input(WALL)
        .replace("ambient mixture 293.15 100000 mole 2 N2 0.7 O2 0.3", "ambient 300 101325 0")
        .replace("baffled_outlet 0 0 0.5 4 16 800",
            "baffled_outlet 0 0 0.5 4 16 800\nexterior 300 101325 0")
        .replace("segments 2\ncylinder 0.0015 0.2\ncylinder 0.0015 0.3 gas mixture 293.15 100000 mole 3 N2 0.4 O2 0.2 CO2 0.4",
            "segments 1\ncylinder 0.0015 0.57");
    let source = field(280.0);
    let steps = fs_duct::RiccatiSteps(16);
    let mut actual = ReedPerformance::from_bytes_with_material_environment(
        text.as_bytes(), 37, ObservationPolicy::Diagnostic, None,
        Some((source.clone(), GasFieldResolution::Continuous(steps))), None,
        &mut |_| Ok(pack.clone()),
    ).unwrap();
    let record = actual.bore_wall().unwrap().clone();
    assert_eq!(record.wall().material().state_coordinate("T"), Some(300.0));
    assert_eq!(actual.segment_gases().len(), 1);
    let duct = Duct { segments: vec![Segment::Cylinder { radius: 0.0015, length: 0.57 }] };
    let info = actual.outlet_receiver().unwrap();
    let regional = RegionalDuct::from_continuous_gas_field(&duct, &source, steps).unwrap()
        .with_exterior_gas(info.gas).unwrap();
    let wall = fs_duct::WallModel::Law(record.wall());
    let mut direct = ReedBoreVoice::new_with_regional_wall(&regional, reed(),
        Termination::FlangedOpen, PlateBank::default(), 0.5, RATE, SAMPLES, wall).unwrap();
    let mut microphone = regional.baffled_outlet_with_wall(
        RATE, info.fft_samples, info.receiver, wall).unwrap();
    microphone.step(direct.p_plus_prev).unwrap();
    let schedule = GestureSchedule::from_canonical_bytes(
        text.split_once("\nschedule\n").unwrap().1.as_bytes()).unwrap();
    for n in 0..SAMPLES {
        direct.set_blowing_pressure(schedule.sample("blow", n as u64 * 700 / u64::from(RATE)).unwrap());
        direct.step_block(&mut [0.0]).unwrap();
        let mut sample = [0.0];
        actual.pressure_block(&mut sample).unwrap();
        assert_eq!(motion(voice(&actual)), motion(&direct));
        assert_eq!(sample[0].to_bits(), microphone.step(direct.p_plus_prev).unwrap().to_bits());
    }
}

#[test]
fn missing_conflicting_and_out_of_domain_wall_sources_refuse() {
    let text = input(WALL);
    let line = text.lines().position(|row| row.starts_with("bore_wall ")).unwrap() + 1;
    assert!(matches!(ReedPerformance::from_bytes(text.as_bytes(), 37),
        Err(ReedPerformanceError::MaterialSource { line: actual, .. }) if actual == line));
    let error = load(&text, &wall_pack(6e7, Some("young_modulus")), 37)
        .err().expect("a missing modulus must not become a rigid wall");
    assert!(matches!(error, ReedPerformanceError::MaterialWall(_)));
    assert!(std::error::Error::source(&error).is_some());
    let pack = wall_pack(6e7, None);
    for invalid in [
        text.replace("0.000001 5000", "0.000001 500"),
        text.replace("wall.matpack", "../wall.matpack"),
        text.replace("300 1000 0.00001", "1000 1000 0.00001"),
        input("bore_wall rigid extra\n"),
    ] {
        assert!(load(&invalid, &pack, 37).is_err(), "accepted {invalid}");
    }
    // Every positive audio bin is covered, but the shunt's original
    // 1e-3 rad/s probe is not. Its failed query must not become zero DC bins.
    let unsupported_dc = text
        .replace("segments 2\ncylinder 0.0015 0.2",
            "segments 4\ncylinder 0.0015 0.1\nhole 0.0005 0.002 0.0015 0.3\ncylinder 0.0015 0.1")
        .replace("0.000001 5000", "0.01 5000");
    let error = load(&unsupported_dc, &pack, 37).err()
        .expect("a shunted material wall must support its near-DC query");
    assert!(error.to_string().contains("source band"), "{error}");
    let molecular = include_str!("../../../../examples/reed-molecular-synthetic.performance")
        .replace("\nmolecular ", &format!("\n{WALL}molecular "));
    assert!(matches!(load(&molecular, &pack, 37),
        Err(ReedPerformanceError::Input { what, .. }) if what.contains("molecular_wall")));
}
