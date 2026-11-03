//! The card row `jet_response forced-bickley SPACING WINDOW [exit kutta|measured]` (MB 14979):
//! no tail is the Kutta exit, `exit measured` integrates the measured in-phase receptivity
//! `eta(0, t) = (h / U0) v` with `R = h / b` (the D4 flute's lip jet: b = 0.75 h, R = 4/3), and
//! any other law is refused. Parsed on the committed D4 flute card at the coarse row
//! (`jet_response forced-bickley 0.25 40`, `forced_jet::COARSE_SPEC`). `b = 0.75 h` (uniform exit)
//! and `b = 0.4 h` (Poiseuille exit) are Verge's matching as Segoufin states it (thesis, Univ.
//! Paris 6, 2000, sec. 5.2 p. 69), the constants of `forced_jet`; `R = h / b` is Terrien et al.'s
//! Eq. 1 (`eta_0 = (h / U_j) v_ac`) in `b / U0` units with `U_j = U0`.

use fs_aeroac::forced_jet::ExitReceptivity;
use fs_couple::flue_loop::forced_jet::{JetExit, JetExitLaw};
use fs_couple::flue_loop::network::performance::file::FlueFilePerformance;
use fs_exec::CancelGate;

const CARD: &str =
    include_str!("../../../data/instruments/descriptions/flute/flute-d4.performance");
const JET: &str = "jet 90 1e-6 3 0\n";

fn parse(tail: &str) -> Result<FlueFilePerformance, String> {
    assert_eq!(CARD.matches(JET).count(), 1, "the D4 jet row changed; update this test");
    let text = CARD.replace(JET, &format!("{JET}jet_response forced-bickley 0.25 40{tail}\n"));
    FlueFilePerformance::from_bytes(text.as_bytes(), 256, &CancelGate::new_clock_free())
        .map_err(|e| format!("{e:?}"))
}

#[test]
fn no_tail_is_kutta() {
    let card = parse("").expect("parses");
    let jet = card.renderer().voice().forced_jet().expect("forced jet attached");
    assert_eq!(jet.exit_law(), JetExitLaw::Kutta);
    assert_eq!(jet.step_response().exit(), ExitReceptivity::Kutta);
    let named = parse(" exit kutta").expect("parses");
    let named = named.renderer().voice().forced_jet().expect("forced jet attached");
    assert_eq!(named.exit_law(), JetExitLaw::Kutta);
    for (a, b) in jet.step_response().values().iter().zip(named.step_response().values()) {
        assert_eq!(a.to_bits(), b.to_bits());
    }
}

#[test]
fn exit_measured_is_h_over_b_of_the_lip_jet() {
    let card = parse(" exit measured").expect("parses");
    let jet = card.renderer().voice().forced_jet().expect("forced jet attached");
    assert_eq!(jet.exit_law(), JetExitLaw::Measured);
    let r = 1.0 / JetExit::UNIFORM.half_width_over_height;
    assert_eq!(
        jet.step_response().exit(),
        ExitReceptivity::Measured { displacement_per_velocity: r }
    );
    assert!((r - 4.0 / 3.0).abs() < 1e-15);
    // The Poiseuille slit's law is h / (0.4 h).
    assert_eq!(
        JetExitLaw::Measured.receptivity(JetExit::POISEUILLE.half_width_over_height),
        ExitReceptivity::Measured { displacement_per_velocity: 2.5 }
    );
    println!("D4 measured exit: R {r}");
}

#[test]
fn other_exit_laws_are_refused() {
    for tail in [" exit fletcher", " exit", " measured", " exit measured extra"] {
        let err = parse(tail).err().unwrap_or_else(|| panic!("'{tail}' must be refused"));
        println!("'{tail}' refused: {err}");
    }
}
