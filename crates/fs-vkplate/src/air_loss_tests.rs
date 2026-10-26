//! Air loss laws (.372): the Bessel quadrature against an independent series,
//! Stuart's cylinder drag against the paper's own limits, and the plate
//! viscous layer's closed form and refusals.

use crate::GasState;
use crate::air_loss::{
    Complex, oscillating_cylinder, plate_viscous_layer, scaled_bessel_k01, sliding_cylinder,
};
use std::f64::consts::{FRAC_PI_4, PI};

const EULER_GAMMA: f64 = 0.577_215_664_901_532_9;

fn humid_air() -> GasState {
    GasState::try_new_moist_air(293.15, 101_325.0, 0.65).expect("humid default air")
}

fn c_mul(a: Complex, b: Complex) -> Complex {
    Complex::new(a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re)
}

fn c_div(a: Complex, b: Complex) -> Complex {
    let d = b.re * b.re + b.im * b.im;
    Complex::new(
        (a.re * b.re + a.im * b.im) / d,
        (a.im * b.re - a.re * b.im) / d,
    )
}

fn c_exp(z: Complex) -> Complex {
    let m = z.re.exp();
    Complex::new(m * z.im.cos(), m * z.im.sin())
}

fn c_ln(z: Complex) -> Complex {
    Complex::new(z.abs().ln(), z.im.atan2(z.re))
}

fn rel(a: Complex, b: Complex) -> f64 {
    Complex::new(a.re - b.re, a.im - b.im).abs() / b.abs()
}

/// K0 and K1 from the ascending series, A&S 9.6.10, 9.6.11 (n = 1) and 9.6.13.
/// Independent of the quadrature; accurate in double precision for |z| <= 2.
fn series_k01(z: Complex) -> (Complex, Complex) {
    let q = c_mul(z, z);
    let q = Complex::new(q.re / 4.0, q.im / 4.0);
    let mut term = Complex::new(1.0, 0.0);
    let mut harmonic = 0.0;
    let (mut i0, mut k0_sum, mut i1, mut k1_sum) = (
        Complex::new(0.0, 0.0),
        Complex::new(0.0, 0.0),
        Complex::new(0.0, 0.0),
        Complex::new(0.0, 0.0),
    );
    for k in 0..60 {
        if k > 0 {
            term = c_mul(term, q);
            term = Complex::new(term.re / f64::from(k * k), term.im / f64::from(k * k));
            harmonic += 1.0 / f64::from(k);
        }
        let over = f64::from(k + 1);
        let psi_sum = 2.0 * (harmonic - EULER_GAMMA) + 1.0 / over;
        i0 = Complex::new(i0.re + term.re, i0.im + term.im);
        k0_sum = Complex::new(
            k0_sum.re + term.re * harmonic,
            k0_sum.im + term.im * harmonic,
        );
        i1 = Complex::new(i1.re + term.re / over, i1.im + term.im / over);
        k1_sum = Complex::new(
            k1_sum.re + psi_sum * term.re / over,
            k1_sum.im + psi_sum * term.im / over,
        );
    }
    let half = Complex::new(z.re / 2.0, z.im / 2.0);
    let ln_half = c_ln(half);
    let lead = Complex::new(ln_half.re + EULER_GAMMA, ln_half.im);
    let k0 = Complex::new(k0_sum.re, k0_sum.im);
    let k0 = Complex::new(k0.re - c_mul(lead, i0).re, k0.im - c_mul(lead, i0).im);
    let i1 = c_mul(half, i1);
    let quarter = Complex::new(z.re / 4.0, z.im / 4.0);
    let inv = c_div(Complex::new(1.0, 0.0), z);
    let a = c_mul(ln_half, i1);
    let b = c_mul(quarter, k1_sum);
    let k1 = Complex::new(inv.re + a.re - b.re, inv.im + a.im - b.im);
    (k0, k1)
}

#[test]
fn bessel_k_quadrature_matches_the_ascending_series_on_the_kelvin_ray() {
    // K0(x e^{i pi/4}) = ker x + i kei x (A&S 9.9.2): the ray Stuart's law uses.
    let mut worst = 0.0_f64;
    for x in [0.001, 0.01, 0.1, 0.5, 1.0, 2.0] {
        let z = Complex::new(x * FRAC_PI_4.cos(), x * FRAC_PI_4.sin());
        let (s0, s1) = scaled_bessel_k01(z).expect("Re z > 0");
        let undo = c_exp(Complex::new(-z.re, -z.im));
        let (q0, q1) = (c_mul(s0, undo), c_mul(s1, undo));
        let (r0, r1) = series_k01(z);
        let e = rel(q0, r0).max(rel(q1, r1));
        worst = worst.max(e);
        println!(
            "{{\"test\":\"bessel_k\",\"x\":{x},\"ker\":{:.12},\"kei\":{:.12},\"rel_diff\":{e:e}}}",
            q0.re, q0.im
        );
        assert!(e <= 1e-13, "x {x}: quadrature and series differ by {e:e}");
    }
    println!("{{\"test\":\"bessel_k\",\"worst_rel_diff\":{worst:e}}}");
}

#[test]
fn bessel_k_refuses_outside_the_right_half_plane() {
    for z in [
        Complex::new(0.0, 1.0),
        Complex::new(-1.0, 0.5),
        Complex::new(f64::NAN, 1.0),
    ] {
        assert!(scaled_bessel_k01(z).is_err(), "{z:?} must refuse");
    }
}

/// a / delta for radius `a` at `omega` in `gas`.
fn radius_for(r: f64, omega: f64, gas: &GasState) -> f64 {
    r * (2.0 * gas.dynamic_viscosity / (gas.density * omega)).sqrt()
}

#[test]
fn stuart_drag_meets_the_papers_low_and_high_frequency_limits() {
    // Dolfo, Vigue & Lhuillier eqs. (50)-(52) and (55)-(56), with the accuracies
    // the paper states: k' within 1 % and k within 7 % for a/delta < 0.1 and a/delta > 1.
    // The low-frequency denominator is L^2 + (pi/4)^2, as the small-argument forms
    // K0 ~ -ln(alpha/2) - gamma, K1 ~ 1/alpha give from eq. (60).
    let gas = humid_air();
    let omega = 2.0 * PI * 500.0;
    let low = |r: f64| {
        let l = -std::f64::consts::LN_2 / 2.0 + EULER_GAMMA + r.ln();
        let d = l * l + PI * PI / 16.0;
        (1.0 + (PI / 2.0) / d / (r * r), (-2.0 * l / d) / (r * r))
    };
    let high = |r: f64| {
        let e = 1.0 / r;
        (
            1.0 + 2.0 * e + e.powi(3) / 8.0,
            2.0 * e + e * e - e.powi(3) / 8.0,
        )
    };
    for (r, limit, k_band) in [
        (0.001, low(0.001), true),
        (0.01, low(0.01), true),
        (0.03, low(0.03), true),
        (0.1, low(0.1), false),
        (1.0, high(1.0), true),
        (3.0, high(3.0), true),
        (10.0, high(10.0), true),
        (100.0, high(100.0), true),
    ] {
        let drag = oscillating_cylinder(radius_for(r, omega, &gas), omega, &gas).expect("drag");
        let (k, kp) = (drag.added_mass_coefficient, drag.damping_coefficient);
        let (ek, ekp) = ((k - limit.0).abs() / k, (kp - limit.1).abs() / kp);
        println!(
            "{{\"test\":\"stuart_limits\",\"a_over_delta\":{r},\"k\":{k:.9},\"k_prime\":{kp:.9},\"limit_k\":{:.9},\"limit_k_prime\":{:.9},\"rel_k\":{ek:e},\"rel_k_prime\":{ekp:e}}}",
            limit.0, limit.1
        );
        assert!((drag.radius_over_penetration() - r).abs() <= 1e-12 * r);
        assert!(
            ekp <= 0.01,
            "a/delta {r}: k' {kp} is {ekp:e} from its limit {}",
            limit.1
        );
        if k_band {
            assert!(
                ek <= 0.07,
                "a/delta {r}: k {k} is {ek:e} from its limit {}",
                limit.0
            );
        }
    }
}

#[test]
fn stuart_drag_per_length_scales_with_the_gas() {
    // Far above the layer (a/delta ~ 1e3) the resistance tends to
    // 2 pi a sqrt(2 rho mu omega): doubling mu at fixed rho multiplies it by sqrt 2.
    let gas = humid_air();
    let omega = 2.0 * PI * 2000.0;
    let a = 0.05;
    let base = oscillating_cylinder(a, omega, &gas).expect("drag");
    let mut viscous = gas;
    viscous.dynamic_viscosity *= 2.0;
    let doubled = oscillating_cylinder(a, omega, &viscous).expect("drag");
    let ratio = doubled.resistance_n_s_per_m2 / base.resistance_n_s_per_m2;
    let asymptote = 2.0 * PI * a * (2.0 * gas.density * gas.dynamic_viscosity * omega).sqrt();
    println!(
        "{{\"test\":\"stuart_scaling\",\"a_over_delta\":{:.1},\"ratio\":{ratio:.6},\"resistance\":{:e},\"asymptote\":{asymptote:e}}}",
        base.radius_over_penetration(),
        base.resistance_n_s_per_m2
    );
    assert!(base.radius_over_penetration() > 500.0);
    assert!((ratio - 2.0_f64.sqrt()).abs() <= 2e-3, "ratio {ratio}");
    assert!((base.resistance_n_s_per_m2 / asymptote - 1.0).abs() <= 2e-3);
    // The added mass tends to the displaced gas, pi a^2 rho.
    assert!((base.added_mass_kg_per_m / (PI * a * a * gas.density) - 1.0).abs() <= 3e-3);
}

#[test]
fn sliding_cylinder_meets_its_asymptotes_and_half_the_normal_drag() {
    // R = 2 pi mu Re(z K1(z) / K0(z)), z = (1 + i) a / delta. Large z (A&S 9.7.2):
    // z K1/K0 = z + 1/2 - 1/(8 z) + O(z^-2), so Re = r + 1/2 - 1/(16 r). Small z
    // (A&S 9.6.13, the K0 ~ -ln(z/2) - gamma and K1 ~ 1/z of the test above):
    // Re = -l / (l^2 + (pi/4)^2), l = ln(r / sqrt 2) + gamma.
    let gas = humid_air();
    let omega = 2.0 * PI * 1000.0;
    let mu = gas.dynamic_viscosity;
    for (r, expected, tol) in [
        (
            0.001,
            -(0.001_f64.ln() - std::f64::consts::LN_2 / 2.0 + EULER_GAMMA),
            0.01,
        ),
        (
            0.01,
            -(0.01_f64.ln() - std::f64::consts::LN_2 / 2.0 + EULER_GAMMA),
            0.01,
        ),
        (3.0, 3.5 - 1.0 / 48.0, 2e-3),
        (10.0, 10.5 - 1.0 / 160.0, 1e-4),
        (100.0, 100.5 - 1.0 / 1600.0, 1e-6),
    ] {
        let re = if r < 1.0 {
            expected / (expected * expected + FRAC_PI_4 * FRAC_PI_4)
        } else {
            expected
        };
        let resistance =
            sliding_cylinder(radius_for(r, omega, &gas), omega, &gas).expect("sliding drag");
        let got = resistance / (2.0 * PI * mu);
        println!(
            "{{\"test\":\"sliding_limits\",\"a_over_delta\":{r},\"re\":{got:.9},\"limit\":{re:.9}}}"
        );
        assert!(
            (got / re - 1.0).abs() <= tol,
            "a/delta {r}: Re(z K1/K0) {got} against {re}"
        );
    }
    // Far above the layer the normal (Stuart) resistance is twice the sliding one.
    let a = radius_for(100.0, omega, &gas);
    let normal = oscillating_cylinder(a, omega, &gas)
        .expect("drag")
        .resistance_n_s_per_m2;
    let ratio = normal / sliding_cylinder(a, omega, &gas).expect("sliding drag");
    assert!((ratio - 2.0).abs() <= 0.02, "normal over sliding {ratio}");
    assert!(sliding_cylinder(0.0, omega, &gas).is_err());
    assert!(sliding_cylinder(1e-3, f64::NAN, &gas).is_err());
}

#[test]
fn oscillating_cylinder_refuses_nonphysical_inputs() {
    let gas = humid_air();
    assert!(oscillating_cylinder(0.0, 1000.0, &gas).is_err());
    assert!(oscillating_cylinder(1e-3, -1.0, &gas).is_err());
    assert!(oscillating_cylinder(f64::NAN, 1000.0, &gas).is_err());
    let mut no_viscosity = gas;
    no_viscosity.dynamic_viscosity = 0.0;
    let refusal = oscillating_cylinder(1e-3, 1000.0, &no_viscosity).unwrap_err();
    assert!(refusal.contains("dynamic viscosity"), "{refusal}");
}

#[test]
fn plate_viscous_layer_is_rho_delta_over_rho_s_with_the_near_field_slip() {
    let gas = humid_air();
    // The default crash bow: 0.95 mm of the kit's B20 at 8600 kg/m^3, 2 kHz.
    let areal = 8600.0 * 0.95e-3;
    let omega = 2.0 * PI * 2000.0;
    let delta = (2.0 * gas.dynamic_viscosity / (gas.density * omega)).sqrt();
    let far_below = plate_viscous_layer(areal, omega, 0.0, 2, &gas).expect("layer");
    assert!((far_below.loss_factor / (gas.density * delta / areal) - 1.0).abs() <= 1e-14);
    assert!((far_below.slip_over_normal - 1.0).abs() <= 1e-15);
    let one_face = plate_viscous_layer(areal, omega, 0.0, 1, &gas).expect("layer");
    assert!((2.0 * one_face.loss_factor / far_below.loss_factor - 1.0).abs() <= 1e-15);
    // f / f_c = 2000 / 18561 for the crash in this air: slip^2 = 1 / (1 - f/f_c).
    let crash = plate_viscous_layer(areal, omega, 2000.0 / 18_561.0, 2, &gas).expect("layer");
    let expected = gas.density * delta / areal / (1.0 - 2000.0 / 18_561.0);
    assert!((crash.loss_factor / expected - 1.0).abs() <= 1e-14);
    println!(
        "{{\"test\":\"plate_viscous_layer\",\"case\":\"crash bow 2 kHz humid air\",\"delta_m\":{:e},\"eta\":{:e},\"t60_s\":{:.1}}}",
        crash.penetration_depth_m,
        crash.loss_factor,
        2.0 * 6.907_755_278_982_137 / (crash.loss_factor * omega)
    );
    for f_over in [1.0, 1.5, -0.1, f64::NAN] {
        assert!(
            plate_viscous_layer(areal, omega, f_over, 2, &gas).is_err(),
            "f/f_c {f_over}"
        );
    }
    for faces in [0, 3] {
        assert!(plate_viscous_layer(areal, omega, 0.1, faces, &gas).is_err());
    }
}
