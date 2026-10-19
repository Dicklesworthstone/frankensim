//! Hammer felt stiffness (`Config::rt0425_hammer_stiffness`). The default stays the reference
//! preset's legacy card, bit for bit. RT-0425's K_H on each unison string, which the browser
//! plays, is the unison count times stiffer per string, so the felt lets go of the strings sooner.
use fs_live_piano::{Config, Description, LivePiano};

const RATE: u32 = 48_000;

/// Pressure and per-sample hammer force for one C4 strike at 2.5 m/s.
fn strike(rt0425: Option<bool>, seconds: f64) -> (Vec<f64>, Vec<f64>) {
    let mut cfg = Config::live(RATE);
    if let Some(on) = rt0425 {
        cfg.rt0425_hammer_stiffness = on;
    }
    let mut p = LivePiano::new(&Description::embedded().unwrap(), cfg).unwrap();
    p.note_on(60, 2.5).unwrap();
    (0..(seconds * f64::from(RATE)) as usize)
        .map(|_| (p.step(), p.hammer_force(60).unwrap_or(0.0)))
        .unzip()
}

/// Samples whose force exceeds 0.1 % of the stroke's peak.
fn contact_samples(force: &[f64]) -> usize {
    let peak = force.iter().copied().fold(0.0, f64::max);
    force.iter().filter(|&&f| f > 1.0e-3 * peak).count()
}

#[test]
fn the_legacy_card_remains_the_default_bit_for_bit() {
    let (default, _) = strike(None, 0.1);
    let (legacy, _) = strike(Some(false), 0.1);
    assert!(default.iter().zip(&legacy).all(|(a, b)| a.to_bits() == b.to_bits()));
    let (published, _) = strike(Some(true), 0.1);
    assert!(default.iter().zip(&published).any(|(a, b)| a != b), "the published stiffness changed nothing");
}

#[test]
fn per_string_stiffness_shortens_the_contact_and_raises_the_force() {
    let (_, legacy) = strike(Some(false), 0.03);
    let (_, published) = strike(Some(true), 0.03);
    let (n_legacy, n_published) = (contact_samples(&legacy), contact_samples(&published));
    let peak = |f: &[f64]| f.iter().copied().fold(0.0, f64::max);
    println!("C4 2.5 m/s contact: legacy {n_legacy} samples, peak {:.1} N; RT-0425 {n_published} samples, peak {:.1} N",
        peak(&legacy), peak(&published));
    assert!(n_published > 0 && n_published < n_legacy, "contact {n_published} vs legacy {n_legacy} samples");
    assert!(peak(&published) > peak(&legacy), "peak force must rise with the stiffer felt");
}

#[test]
fn published_dissipation_needs_the_published_stiffness() {
    let mut cfg = Config::live(RATE);
    cfg.rt0425_hammer_dissipation = true;
    assert!(LivePiano::new(&Description::embedded().unwrap(), cfg).is_err());
}

#[test]
fn published_dissipation_books_its_loss_and_closes_the_ledger() {
    // R_H replaces the estimated crush and Prony terms, and the K_H card is conservative, so
    // every joule of felt loss is R_H's (measured 3.5 % of the input at 0.5 m/s, 6.5 % at 4.5).
    // The ledger explains the input to roundoff, soft and hard.
    for velocity in [0.5, 4.5] {
        let mut cfg = Config::live(RATE);
        cfg.rt0425_hammer_stiffness = true;
        cfg.rt0425_hammer_dissipation = true;
        let mut p = LivePiano::new(&Description::embedded().unwrap(), cfg).unwrap();
        p.note_on(60, velocity).unwrap();
        for _ in 0..(0.1 * f64::from(RATE)) as usize {
            p.step();
        }
        let l = p.close_ledger();
        println!("C4 {velocity} m/s with R_H: input {:.4e} J, felt loss {:.4e} J, residual {:.2e} J", l.input_j, l.felt_loss_j, l.residual_j);
        assert!(l.felt_loss_j > 1.0e-2 * l.input_j, "felt loss {:e} of {:e} J", l.felt_loss_j, l.input_j);
        assert!(l.residual_j.abs() < 1.0e-9 * l.input_j, "residual {:e} of {:e} J", l.residual_j, l.input_j);
    }
}
