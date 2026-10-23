//! Hammers from the felt pack (`Config::felt_pack_hammers`), which the browser plays:
//! Stulov's 88 measured Abel hammers instead of RT-0425's per-key K_H and R_H.
use fs_live_piano::{Config, Description, LivePiano};

const RATE: u32 = 48_000;

fn config(felt_pack: bool) -> Config {
    let mut cfg = Config::live(RATE);
    if felt_pack {
        cfg.felt_pack_hammers = true;
    } else {
        cfg.rt0425_hammer_stiffness = true;
        cfg.rt0425_hammer_dissipation = true;
    }
    cfg
}

/// Pressure and per-sample hammer force for one strike.
fn strike(felt_pack: bool, midi: u8, velocity: f64, seconds: f64) -> (Vec<f64>, Vec<f64>) {
    let mut p = LivePiano::new(&Description::embedded().unwrap(), config(felt_pack)).unwrap();
    p.note_on(midi, velocity).unwrap();
    (0..(seconds * f64::from(RATE)) as usize)
        .map(|_| (p.step(), p.hammer_force(midi).unwrap_or(0.0)))
        .unzip()
}

/// Samples whose force exceeds 0.1 % of the stroke's peak.
fn contact_samples(force: &[f64]) -> usize {
    let peak = force.iter().copied().fold(0.0, f64::max);
    force.iter().filter(|&&f| f > 1.0e-3 * peak).count()
}

#[test]
fn felt_pack_hammers_book_their_loss_and_close_the_ledger() {
    for velocity in [0.5, 4.5] {
        let mut p = LivePiano::new(&Description::embedded().unwrap(), config(true)).unwrap();
        p.note_on(60, velocity).unwrap();
        for _ in 0..(0.1 * f64::from(RATE)) as usize {
            p.step();
        }
        let l = p.close_ledger();
        println!(
            "C4 {velocity} m/s, felt pack: input {:.4e} J, felt loss {:.4e} J, residual {:.2e} J",
            l.input_j, l.felt_loss_j, l.residual_j
        );
        assert!(
            l.felt_loss_j > 1.0e-2 * l.input_j,
            "felt loss {:e} of {:e} J",
            l.felt_loss_j,
            l.input_j
        );
        assert!(
            l.residual_j.abs() < 1.0e-9 * l.input_j,
            "residual {:e} of {:e} J",
            l.residual_j,
            l.input_j
        );
    }
}

/// The measured hammers are not RT-0425's: the strike changes. The contact
/// times and peak forces are reported, not tuned toward either card.
#[test]
fn felt_pack_hammers_change_the_strike_from_rt0425() {
    for midi in [36, 60, 84] {
        let (pressure_pack, force_pack) = strike(true, midi, 2.5, 0.03);
        let (pressure_rt, force_rt) = strike(false, midi, 2.5, 0.03);
        let peak = |f: &[f64]| f.iter().copied().fold(0.0, f64::max);
        println!(
            "key {midi} at 2.5 m/s: contact RT-0425 {} / felt pack {} samples, peak {:.1} / {:.1} N",
            contact_samples(&force_rt),
            contact_samples(&force_pack),
            peak(&force_rt),
            peak(&force_pack)
        );
        assert!(contact_samples(&force_pack) > 0 && peak(&force_pack) > 0.0);
        assert!(
            pressure_pack.iter().zip(&pressure_rt).any(|(a, b)| a != b),
            "key {midi}: the felt pack changed nothing"
        );
    }
}
