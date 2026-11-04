use super::*;
use crate::compose::{ComposeOptions, tone_hole_scattering};
use crate::{HoleState, LossModel, Segment, WallModel};
use fs_material::gas::{GasSpec, GasState};
use fs_math::c64::C64;

const TAU: f64 = core::f64::consts::TAU;
const PI: f64 = core::f64::consts::PI;
// A clarinet-like hole: radius 4 mm on a 7.5 mm bore (d = 8/15), chimney 6 mm.
const B: f64 = 4.0e-3;
const A: f64 = 7.5e-3;
const T: f64 = 6.0e-3;

fn air() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).expect("air")
}

fn near(a: f64, b: f64, tolerance: f64) -> bool {
    (a - b).abs() <= tolerance * a.abs().max(b.abs())
}

fn close(a: C64, b: C64, tolerance: f64) -> bool {
    (a - b).abs() <= tolerance * a.abs().max(b.abs())
}

/// The block of one hole as the chain builds it (lossless, rigid walls), `[A, B, C, D]`.
fn block(lengths: HoleLengths, state: HoleState, omega: f64) -> [C64; 4] {
    let gas = air();
    let segment = Segment::ToneHole { hole_radius: B, chimney_height: T, bore_radius: A, state };
    let reference = gas.density * gas.sound_speed / (PI * A * A);
    tone_hole_scattering(
        &segment, &gas, &gas, omega, LossModel::Lossless, WallModel::Rigid, reference, 0.0, lengths,
        ComposeOptions::default().maximum_condition,
    )
    .unwrap()
    .abcd()
    .unwrap()
}

#[test]
fn the_engine_terms_are_the_original_functions_bit_for_bit() {
    assert_eq!(HoleLengths::from_compact(true), HoleLengths::Engine);
    assert_eq!(HoleLengths::from_compact(false), HoleLengths::Uncorrected);
    assert_eq!(HoleLengths::default(), HoleLengths::Engine);
    for (b, a, t) in [(B, A, T), (2.5e-3, 7.35e-3, 6.5e-3), (7.0e-3, 7.5e-3, 3.0e-3), (8.0e-3, 7.5e-3, 3.0e-3)] {
        let engine = HoleLengths::Engine;
        assert_eq!(engine.open_extra_m(b, t, a, 40.0), fs_phs::side_hole_inner_length(b, a));
        assert_eq!(engine.sealed_extra_m(b, a), 0.0);
        let series = fs_phs::side_hole_series_length(b, a);
        assert_eq!(engine.series_m(b, t, a, 1.0, 0.0), series * 1.0 + 0.0);
        assert_eq!(engine.series_m(b, t, a, 0.3, 1.0e-5), series * 0.3 + 1.0e-5);
        // A closed hole has no series term in the engine, whatever mutual length is offered.
        assert_eq!(engine.series_m(b, t, a, 0.0, 1.0e-5), 0.0);
        let none = HoleLengths::Uncorrected;
        assert_eq!((none.open_extra_m(b, t, a, 40.0), none.sealed_extra_m(b, a), none.series_m(b, t, a, 1.0, 1.0e-5)), (0.0, 0.0, 0.0));
    }
}

#[test]
fn the_printed_lengths_follow_the_printed_equations_for_a_hand_computed_hole() {
    // d = 8/15, ka = 0.3 (k = 40 /m). Each value below was evaluated by hand from the Lefebvre and Scavone
    // (2012) equation named beside it, in plain power form, independently of this module's Horner form:
    //   Eq. 31 t_i/b = 0.539207619423868, Eq. 32 G = 1.058532386666667,
    //   Eq. 8 t_m = 2.75040711111111e-4, Eq. 11 t_r = 2.57754131520697e-3,
    //   Eq. 33 t_a(o) = -3.29996986190055e-4, Eq. 34 t_a(c) = -3.29666958445790e-4.
    let printed = HoleLengths::Printed;
    let open = 2.28307491319039609e-3 + 2.75040711111111127e-4 + (2.57754131520696690e-3 - 0.8216 * B);
    assert!(near(printed.open_extra_m(B, T, A, 40.0), open, 1.0e-12), "{}", printed.open_extra_m(B, T, A, 40.0));
    assert!(near(printed.sealed_extra_m(B, A), 2.75040711111111127e-4, 1.0e-12));
    let (ta_open, ta_closed) = (-3.29996986190055498e-4, -3.29666958445790450e-4);
    assert!(near(printed.series_m(B, T, A, 1.0, 0.0), ta_open, 1.0e-12));
    assert!(near(printed.series_m(B, T, A, 0.0, 0.0), ta_closed, 1.0e-12));
    // The mutual length belongs to the open part only.
    assert!(near(printed.series_m(B, T, A, 1.0, 2.0e-5), ta_open + 2.0e-5, 1.0e-12));
    assert!(near(printed.series_m(B, T, A, 0.0, 2.0e-5), ta_closed, 1.0e-12));
    assert!(near(printed.series_m(B, T, A, 0.25, 0.0), 0.25 * ta_open + 0.75 * ta_closed, 1.0e-12));
    assert!(near(printed.series_m(B, T, A, 0.25, 0.0), -3.29749465381856753e-4, 1.0e-12));
}

#[test]
fn a_hole_outside_the_bore_keeps_the_engine_terms_and_the_frequency_factor_stops_at_ka_one() {
    for b in [A, 1.2 * A] {
        assert_eq!(HoleLengths::Printed.open_extra_m(b, T, A, 40.0), HoleLengths::Engine.open_extra_m(b, T, A, 40.0));
        assert_eq!(HoleLengths::Printed.sealed_extra_m(b, A), 0.0);
        assert_eq!(HoleLengths::Printed.series_m(b, T, A, 1.0, 0.0), HoleLengths::Engine.series_m(b, T, A, 1.0, 0.0));
    }
    // Eq. 32 is printed for ka <= 1: ka = 3 holds the ka = 1 value instead of the quartic's excursion.
    let at_one = HoleLengths::Printed.open_extra_m(B, T, A, 1.0 / A);
    let beyond = HoleLengths::Printed.open_extra_m(B, T, A, 3.0 / A);
    assert!(near(at_one, beyond, 1.0e-12), "{at_one} vs {beyond}");
    assert!(at_one.is_finite() && at_one > 0.0);
    // And the factor grows from ka = 0 for this d (H(d) > 0).
    assert!(HoleLengths::Printed.open_extra_m(B, T, A, 0.0) < at_one);
}

#[test]
fn a_closed_hole_is_the_printed_symmetric_t_section_and_the_engine_closed_hole_is_a_pure_shunt() {
    // Eq. 2 with Eq. 16 (k t_i dropped, see the module documentation) and Eq. 6, for exp(-i w t):
    //   Zs = i Z0h cot k(t + t_m),  Za = -i k Z0 t_a(c),
    //   T = [[1 + Za/2Zs, Za (1 + Za/4Zs)], [1/Zs, 1 + Za/2Zs]].
    let gas = air();
    let omega = TAU * 1500.0;
    let k = omega / gas.sound_speed;
    let z0h = gas.density * gas.sound_speed / (PI * B * B);
    let z0 = gas.density * gas.sound_speed / (PI * A * A);
    let d = B / A;
    let tm = B * d * (1.0 + 0.207 * d.powi(3)) / 8.0;
    let zs = C64::new(0.0, z0h / (k * (T + tm)).tan());
    let ta = (-0.12 - 0.17 * (2.4 * T / B).tanh()) * B * d * d;
    let za = C64::new(0.0, -k * z0 * ta);
    let half = C64::ONE + za / zs.scale(2.0);
    let expected = [half, za * (C64::ONE + za / zs.scale(4.0)), zs.recip(), half];
    let got = block(HoleLengths::Printed, HoleState::Closed, omega);
    for (g, e) in got.iter().zip(&expected) {
        assert!(close(*g, *e, 1.0e-9), "{g:?} vs {e:?}");
    }
    // The engine's closed hole: Zs = i Z0h cot(k t) and nothing in series.
    let plain = block(HoleLengths::Engine, HoleState::Closed, omega);
    let zs_engine = C64::new(0.0, z0h / (k * T).tan());
    assert!(close(plain[2], zs_engine.recip(), 1.0e-9), "{:?}", plain[2]);
    assert!(plain[1].abs() < 1.0e-9 * z0, "{:?}", plain[1]);
    assert!(close(plain[0], C64::ONE, 1.0e-9));
}

#[test]
fn an_open_hole_is_the_printed_open_shunt_and_series_at_low_frequency() {
    // Eq. 7 reactance, -Z0h [k t_i G + tan k(t + t_m + t_r)] for exp(-i w t), with the series Za = -i k Z0 t_a(o) of
    // Eq. 4, at 100 Hz where the flanged mouth is its 0.8216 b mass to better than 1e-4. The engine puts the inner
    // length into the line, not beside the tangent: the difference is second order in k t_i.
    let gas = air();
    let omega = TAU * 100.0;
    let k = omega / gas.sound_speed;
    let z0h = gas.density * gas.sound_speed / (PI * B * B);
    let z0 = gas.density * gas.sound_speed / (PI * A * A);
    let d = B / A;
    let ka = k * A;
    let ti = (0.822 - 0.095 * d - 1.566 * d.powi(2) + 2.138 * d.powi(3) - 1.640 * d.powi(4) + 0.502 * d.powi(5)) * B
        * (1.0 + (1.0 - 4.56 * d + 6.55 * d * d) * (0.17 * ka + 0.92 * ka.powi(2) + 0.16 * ka.powi(3) - 0.29 * ka.powi(4)));
    let tm = B * d * (1.0 + 0.207 * d.powi(3)) / 8.0;
    let tr = (0.822 - 0.47 * (B / (A + T)).powf(0.8)) * B;
    let reactance = -z0h * (k * ti + (k * (T + tm + tr)).tan());
    let ta = (-0.35 + 0.06 * (2.7 * T / B).tanh()) * B * d * d;
    let za = C64::new(0.0, -k * z0 * ta);
    let got = block(HoleLengths::Printed, HoleState::Open, omega);
    let zs = got[2].recip();
    assert!(near(zs.im, reactance, 1.0e-3), "{} vs {reactance}", zs.im);
    // The series element of the block: B = Za (1 + Za/4Zs) and A = 1 + Za/2Zs for the engine's own Zs.
    let half = C64::ONE + za / zs.scale(2.0);
    assert!(close(got[0], half, 1.0e-9) && close(got[3], half, 1.0e-9));
    assert!(close(got[1], za * (C64::ONE + za / zs.scale(4.0)), 1.0e-9));
}
