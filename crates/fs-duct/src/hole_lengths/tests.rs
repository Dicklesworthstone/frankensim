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
fn a_hole_outside_the_bore_keeps_the_engine_terms() {
    for b in [A, 1.2 * A] {
        assert_eq!(HoleLengths::Printed.open_extra_m(b, T, A, 40.0), HoleLengths::Engine.open_extra_m(b, T, A, 40.0));
        assert_eq!(HoleLengths::Printed.sealed_extra_m(b, A), 0.0);
        assert_eq!(HoleLengths::Printed.series_m(b, T, A, 1.0, 0.0), HoleLengths::Engine.series_m(b, T, A, 1.0, 0.0));
    }
}

/// `G` of Eq. 32 exactly as this module computed it before the continuation above `ka = 1`: the printed quartic in
/// Horner form on `ka` clamped to `0..=1`.
fn printed_factor_as_it_was(d: f64, ka: f64) -> f64 {
    let x = ka.clamp(0.0, 1.0);
    let h = 1.0 + d * (-4.56 + d * 6.55);
    let i = x * (0.17 + x * (0.92 + x * (0.16 - 0.29 * x)));
    1.0 + h * i
}

/// Planted test 1 (MB 18434): the printed quartic is untouched for `ka <= 1`, bit for bit.
#[test]
fn the_frequency_factor_is_bit_identical_to_the_printed_quartic_up_to_ka_one() {
    let mut grid = vec![-0.5, -1.0e-12, 0.0, 1.0e-9, 0.1, 0.25, 0.5, 0.75, 0.999, 1.0 - f64::EPSILON, 1.0];
    grid.extend((0..=200).map(|i| f64::from(i) / 200.0));
    for d in [0.0, 0.1, 8.0 / 15.0, 0.9, 1.0] {
        for &ka in &grid {
            assert_eq!(
                frequency_factor(d, ka).to_bits(),
                printed_factor_as_it_was(d, ka).to_bits(),
                "d = {d}, ka = {ka}"
            );
        }
    }
}

/// Planted test 2 (MB 18434): above `ka = 1` the declared continuation `I(1) + I'(1) tanh(ka - 1)` takes the printed
/// value and slope at the edge (`I(1) = 0.96`, `I'(1) = 0.17 + 2 * 0.92 + 3 * 0.16 - 4 * 0.29 = 1.33`, so `G - 1 = H(d) I`
/// has slope `1.33 H(d)`) and saturates at `I(1) + I'(1) = 2.29`. The old hold had slope 0 to the right of `ka = 1`, so
/// the slope assertion FAILS on the old file: the kink that gave the reflectance its pre-arrival skirt.
#[test]
fn above_ka_one_the_declared_continuation_has_the_printed_value_and_slope_and_saturates() {
    for d in [0.2_f64, 8.0 / 15.0, 0.8] {
        let h = 1.0 + d * (-4.56 + d * 6.55);
        assert!(h.abs() > 0.3, "d = {d}: H = {h} keeps the slope test sensitive");
        let g = |ka: f64| frequency_factor(d, ka);
        // Value at the edge, continuous across it.
        assert!(near(g(1.0), 1.0 + h * 0.96, 1.0e-13), "d = {d}: {}", g(1.0));
        assert!((g(1.0 + 1.0e-9) - g(1.0 - 1.0e-9)).abs() < 1.0e-8, "d = {d}");
        // Slope at the edge, the same from both sides: the printed quartic's 1.33 H(d).
        let step = 1.0e-7;
        let (left, right) = ((g(1.0) - g(1.0 - step)) / step, (g(1.0 + step) - g(1.0)) / step);
        assert!(near(left, 1.33 * h, 1.0e-6), "d = {d}: left slope {left} vs {}", 1.33 * h);
        assert!(near(right, 1.33 * h, 1.0e-6), "d = {d}: right slope {right} vs {}", 1.33 * h);
        // Saturation at I(1) + I'(1) = 2.29, approached from below without overshoot (H > 0 here; mirrored if H < 0).
        let ceiling = 1.0 + h * 2.29;
        assert!(near(g(40.0), ceiling, 1.0e-12), "d = {d}: {} vs {ceiling}", g(40.0));
        let mut last = g(1.0);
        for i in 1..=400 {
            let value = g(1.0 + f64::from(i) * 0.1);
            assert!((value - last) * h >= 0.0, "d = {d}: G turned back at ka = {}", 1.0 + f64::from(i) * 0.1);
            assert!((ceiling - value) * h >= -1.0e-12, "d = {d}: G passed its ceiling at ka = {}", 1.0 + f64::from(i) * 0.1);
            last = value;
        }
    }
}

#[test]
fn the_open_hole_inner_length_past_ka_one_follows_the_declared_continuation() {
    // ka = 3 is no longer the ka = 1 value: only t_i depends on ka, so the open length differs by t_i/b * b * H * 1.33 tanh(2).
    let d = B / A;
    let h = 1.0 + d * (-4.56 + d * 6.55);
    let at_one = HoleLengths::Printed.open_extra_m(B, T, A, 1.0 / A);
    let beyond = HoleLengths::Printed.open_extra_m(B, T, A, 3.0 / A);
    assert!(at_one.is_finite() && at_one > 0.0 && beyond > at_one, "{at_one} vs {beyond}");
    let expected = inner_over_b(d) * B * h * 1.33 * 2.0_f64.tanh();
    assert!(near(beyond - at_one, expected, 1.0e-9), "{} vs {expected}", beyond - at_one);
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
