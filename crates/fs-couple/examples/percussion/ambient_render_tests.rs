//! Only the declared air changes; the drum and cymbal respond physically.
use super::*;
use ambient::Ambient;

const AIR: &str = "293.15 101325 0";
// The shipped descriptions' declared environment (ruling 6205).
const HUMID: &str = "293.15 101325 0.65";
const THIN: &str = "293.15 70927.5 0";
const CO2: &str = "mixture 293.15 101325 mole 1 CO2 1";
const HELIOX: &str = "mixture 293.15 101325 mole 2 He 0.8 O2 0.2";

fn described(record: &str) -> Ambient { Ambient::described(record).unwrap() }

#[test]
fn the_drum_compresses_exactly_its_declared_air() {
    let reference = drum_with_spec(4,2e-6,false,false,None,false,Stroke::default(),false,None,None).unwrap();
    assert_eq!(reference.ambient, Ambient::reference());
    for record in [AIR, THIN, CO2, HELIOX] {
        let declared = described(record);
        let spec = drum_spec::Spec {ambient: declared.clone(), ..drum_spec::Spec::reference()};
        let e = drum_with_spec(4,2e-6,false,false,None,false,Stroke::default(),false,None,Some(spec)).unwrap();
        assert_eq!(e.ambient, declared);
        // The enclosed volume spring is the adiabatic gamma p of the declared gas.
        let bulk = e.pressure.as_ref().unwrap().bulk_modulus_pa;
        assert!((bulk / (declared.gas().gamma * declared.gas().pressure) - 1.0).abs() < 1e-12, "{record}: {bulk}");
    }
    let thin = described(THIN).bulk_modulus_pa() / described(AIR).bulk_modulus_pa();
    assert!((thin - 0.7).abs() < 1e-12, "{thin}");
}

fn splash_pcm(record: &str, frames: usize) -> Vec<i16> {
    let mut shell = specimen::Specimen::reference();
    shell.azimuths = 8;
    shell.ambient = Some(described(record));
    let steps = (frames * acoustics::SUBSTEPS) as u64;
    let mut e = splash_with_mount(steps,acoustics::MECHANICAL_DT,true,Stroke::default(),Some(shell),&[],None,
        None,&mallets::Selection::default(),&shaft_playing::Selection::default(),None,None,None).unwrap();
    assert_eq!(e.ambient.gas(), described(record).gas());
    let wav = acoustics::render(&mut e,frames,2.0,acoustics::Receiver::FarField([1.5,0.7,1.5])).unwrap();
    wav[44..].chunks_exact(2).map(|b| i16::from_le_bytes([b[0],b[1]])).collect()
}

#[test]
fn the_cymbal_radiates_with_its_declared_density_and_sound_speed() {
    // No enclosed air couples the shell, so its motion ignores the gas. At equal c
    // the far field therefore scales with rho; a slower gas delays the wavefront.
    let frames = 1200;
    let air = splash_pcm(AIR, frames);
    let thin = splash_pcm(THIN, frames);
    let co2 = splash_pcm(CO2, frames);
    let heliox = splash_pcm(HELIOX, frames);
    let peak = air.iter().map(|s| s.unsigned_abs()).max().unwrap();
    assert!(peak > 100 && peak < i16::MAX.unsigned_abs(), "unclipped audible reference peak {peak} LSB");
    let ratio = described(THIN).gas().density / described(AIR).gas().density;
    for (n, (a, t)) in air.iter().zip(&thin).enumerate() {
        assert!((f64::from(*t) - ratio * f64::from(*a)).abs() <= 1.0, "sample {n}: {a} {t}");
    }
    // The CO2 wavefront lags by the path length times (1/c_co2 - 1/c_air). The
    // pure delay covers (range - R) for the enclosing radius R = 0.1016 m; the
    // fitted near-source transfer adds between 0 and 2R of path. Measure the lag
    // by cross-correlation, which a density change cannot bias.
    let (c_air, c_co2) = (described(AIR).gas().sound_speed, described(CO2).gas().sound_speed);
    let per_metre = (1.0 / c_co2 - 1.0 / c_air) * f64::from(acoustics::OUTPUT_RATE);
    let (range, radius) = ((1.5_f64 * 1.5 + 0.7 * 0.7 + 1.5 * 1.5).sqrt(), 0.1016);
    let correlation = |lag: usize| air.iter().zip(&co2[lag..])
        .map(|(a, c)| f64::from(*a) * f64::from(*c)).sum::<f64>();
    let lag = (0..400).max_by(|a, b| correlation(*a).total_cmp(&correlation(*b))).unwrap() as f64;
    let (low, high) = ((range - radius) * per_metre, (range + radius) * per_metre);
    assert!(low > 50.0 && lag >= low - 2.0 && lag <= high + 2.0, "lag {lag} outside [{low}, {high}]");
    // 80/20 He-O2 is faster: air now lags the heliox wavefront.
    let c_heliox = described(HELIOX).gas().sound_speed;
    let ahead = (1.0 / c_air - 1.0 / c_heliox) * f64::from(acoustics::OUTPUT_RATE);
    let correlation = |lag: usize| heliox.iter().zip(&air[lag..])
        .map(|(h, a)| f64::from(*h) * f64::from(*a)).sum::<f64>();
    let lag = (0..400).max_by(|a, b| correlation(*a).total_cmp(&correlation(*b))).unwrap() as f64;
    let (low, high) = ((range - radius) * ahead, (range + radius) * ahead);
    assert!(low > 50.0 && lag >= low - 2.0 && lag <= high + 2.0, "heliox lead {lag} outside [{low}, {high}]");
}

#[test]
fn imported_descriptions_must_declare_air_and_bad_records_refuse() {
    let text = include_str!("estimated_drum.fsd");
    assert!(drum_spec::Spec::read(text).is_ok());
    let missing: String = text.lines().filter(|l| !l.starts_with("ambient,")).map(|l| format!("{l}\n")).collect();
    let error = drum_spec::Spec::read(&missing).unwrap_err().to_string();
    assert!(error.contains("ambient"), "{error}");
    for bad in ["ambient,293.15 101325", "ambient,mixture 293.15 101325 mole 1 Xe 1", "ambient,293.15 0 0"] {
        assert!(drum_spec::Spec::read(&text.replace("ambient,293.15 101325 0.65", bad)).is_err(), "{bad}");
    }
    let mut shell = specimen::Specimen::parse(include_str!("estimated-splash.profile")).unwrap();
    assert_eq!(shell.ambient.as_ref().unwrap().gas(), described(HUMID).gas());
    shell.ambient = None;
    let refused = splash_with_mount(4,2e-6,false,Stroke::default(),Some(shell),&[],None,None,
        &mallets::Selection::default(),&shaft_playing::Selection::default(),None,None,None);
    assert!(refused.is_err_and(|e| e.to_string().contains("ambient")));
}
