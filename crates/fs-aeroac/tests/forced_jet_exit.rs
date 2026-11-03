//! The forced jet's opt-in measured exit receptivity (fs_aeroac::forced_jet::ExitReceptivity,
//! MB 14979): the jet issues displaced by `eta(0, t) = R v_ac(t)`, in phase, with no carried term
//! (Terrien, Vergez, de la Cuadra and Fabre, arXiv 1402.3143, Eq. 1, after de la Cuadra; thesis,
//! Stanford 2005, sec. 2.3.5, Fig. 2.28: V_ac and eta_0 in phase). The Kutta default is the same
//! arithmetic as before; the measured exit answers R in phase at the exit at low frequency; a
//! wrong-direction plant (the same response pi out of phase) and the Kutta seed both fail that
//! check.

use fs_aeroac::forced_jet::{
    ExitReceptivity, ForcedJetResolution, bickley_step_response, bickley_step_response_with_exit,
};
use fs_math::c64::C64;
use fs_math::det;

const WORKERS: usize = 4;
/// The lip jet's `h / b` (b = 0.75 h): the measured receptivity with `U_j = U0`. `b = 0.75 h` is
/// the uniform exit's Bickley half-width from the momentum conserved between exit and jet (Verge's
/// matching as Segoufin states it, thesis, Univ. Paris 6, 2000, sec. 5.2 p. 69; fs-couple
/// `flue_loop::forced_jet::UNIFORM_HALF_WIDTH_OVER_HEIGHT`).
const LIP_R: f64 = 4.0 / 3.0;
/// The Poiseuille slit's `h / b` (b = 0.4 h, the same matching for the Poiseuille exit, Segoufin
/// 2000 sec. 5.2 p. 69; fs-couple `flue_loop::forced_jet::POISEUILLE_HALF_WIDTH_OVER_HEIGHT`).
const SLIT_R: f64 = 2.5;
/// Two coarse cells downstream of the exit: the displacement there is the exit's. A numerical
/// station of this check, not physics.
const NEAR_EXIT: f64 = 0.5;
/// The low frequency of the check [w = omega b / U0]: the in-phase law is the quasi-static limit
/// (de la Cuadra 2005 Fig. 2.28 measures V_ac and eta_0 in phase), so the check sits at w << 1.
/// A numerical choice, not physics.
const LOW_W: f64 = 0.01;

/// Numerics, not physics: the coarse row (`spacing 0.25`, `window 40`, fs-couple
/// `flue_loop::forced_jet::COARSE_SPEC`), the band half-height 8 b and Courant 0.8 (fs-couple
/// `BAND_HALF_HEIGHT`, `COURANT`), Re 600 as in `tests/forced_jet.rs`, and 256 cells (64 b) to
/// hold the window's front plus `FRONT_MARGIN`; the solver refuses a shorter box.
fn coarse() -> ForcedJetResolution {
    ForcedJetResolution {
        spacing: 0.25,
        streamwise_cells: 256,
        half_height: 8.0,
        window: 40.0,
        courant: 0.8,
        reynolds: 600.0,
        workers: WORKERS,
    }
}

fn measured(r: f64) -> ExitReceptivity {
    ExitReceptivity::Measured { displacement_per_velocity: r }
}

/// The cited law at the exit: the near-exit transfer at low frequency is `R` [b per v/U0] in
/// phase: `|G / R - 1| < 1 %` and `|arg G| < 0.01 rad` at `w = 0.01` (test tolerances, not
/// physics: the coarse grid's error, far inside the pi a wrong direction would give).
fn check_in_phase_at_exit(g: C64, r: f64) -> Result<(), String> {
    let (mag, arg) = (g.abs(), det::atan2(g.im, g.re));
    println!("near-exit G({LOW_W}) = {mag:.6} at {arg:+.6} rad against R {r:.6}");
    if (mag / r - 1.0).abs() >= 0.01 {
        return Err(format!("|G| {mag:.6} is not R {r:.6}"));
    }
    if arg.abs() >= 0.01 {
        return Err(format!("arg G {arg:+.6} rad is not in phase"));
    }
    Ok(())
}

#[test]
fn kutta_is_the_default_and_unchanged() {
    assert_eq!(ExitReceptivity::default(), ExitReceptivity::Kutta);
    let res = coarse();
    // Station 11 b: a numerical choice inside the box; the bit identity holds at any station.
    let old = bickley_step_response(11.0, &res).expect("kutta");
    let named = bickley_step_response_with_exit(11.0, &res, ExitReceptivity::Kutta).expect("named");
    assert_eq!(old.exit(), ExitReceptivity::Kutta);
    assert_eq!(old.values().len(), named.values().len());
    for (a, b) in old.values().iter().zip(named.values()) {
        assert_eq!(a.to_bits(), b.to_bits(), "the named Kutta exit must be the same bits");
    }
}

#[test]
fn measured_exit_answers_r_in_phase_at_the_exit() {
    let res = coarse();
    let lip = bickley_step_response_with_exit(NEAR_EXIT, &res, measured(LIP_R)).expect("lip");
    check_in_phase_at_exit(lip.transfer(LOW_W).expect("transfer"), LIP_R)
        .expect("the lip jet's measured exit is h/b in phase");
    // Linear in R: the slit's receptivity scales every value by 2.5 / (4/3) (to 1e-12 relative,
    // a rounding tolerance).
    let slit = bickley_step_response_with_exit(NEAR_EXIT, &res, measured(SLIT_R)).expect("slit");
    check_in_phase_at_exit(slit.transfer(LOW_W).expect("transfer"), SLIT_R)
        .expect("the slit jet's measured exit is h/b in phase");
    for (a, b) in lip.values().iter().zip(slit.values()) {
        assert!((b - a * (SLIT_R / LIP_R)).abs() <= 1e-12 * (1.0 + b.abs()));
    }
    // Nothing is carried: the held exit displacement is R, not the Kutta ramp (1 %, the same
    // test tolerance as the exit check).
    assert!((lip.held_gain() / LIP_R - 1.0).abs() < 0.01, "held {}", lip.held_gain());
}

#[test]
fn planted_wrong_direction_and_kutta_fail_the_exit_check() {
    let res = coarse();
    let lip = bickley_step_response_with_exit(NEAR_EXIT, &res, measured(LIP_R)).expect("lip");
    let g = lip.transfer(LOW_W).expect("transfer");
    check_in_phase_at_exit(g, LIP_R).expect("the real response passes");
    // Plant 1: the same response pi out of phase (the wrong direction) must fail.
    let flipped = C64::new(-g.re, -g.im);
    let err = check_in_phase_at_exit(flipped, LIP_R).expect_err("pi out of phase must fail");
    println!("flipped plant refused: {err}");
    // Plant 2: the Kutta seed (Fletcher's 1/omega, nothing at the exit) must fail.
    let kutta = bickley_step_response_with_exit(NEAR_EXIT, &res, ExitReceptivity::Kutta)
        .expect("kutta");
    let err = check_in_phase_at_exit(kutta.transfer(LOW_W).expect("transfer"), LIP_R)
        .expect_err("the Kutta seed is not the measured exit");
    println!("kutta plant refused: {err}");
}

#[test]
fn measured_exit_refuses_non_positive_receptivity() {
    let res = coarse();
    for r in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(bickley_step_response_with_exit(11.0, &res, measured(r)).is_err(), "R {r}");
    }
}
