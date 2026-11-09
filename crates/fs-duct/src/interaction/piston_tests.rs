//! The baffled-piston form of the open-hole mutual impedance (`OpenHoleInteraction::with_piston`): limits, the measured error of the large-distance product against the exact uniform-piston
//! pair resistance, and the passivity of the radiation load on the five saxophone page geometries with the engine's own self load.

use super::*;
use crate::{FLANGED_MODES, flanged_pipe_impedance};
use fs_material::gas::GasSpec;

const TAU: f64 = core::f64::consts::TAU;

fn air() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).expect("air")
}

fn mono() -> OpenHoleInteraction {
    OpenHoleInteraction::new(2.0).unwrap()
}

fn piston() -> OpenHoleInteraction {
    mono().with_piston().unwrap()
}

/// The five page fingerings of the WIAT Selmer Super Action Series II alto lattice (cards `Sax.alto.{C4,E3,F3,Ab4,E5}.lattice.performance`): per fingering the open holes in
/// axial order as (axial station [m], mouth radius [m]), then the bell as (bore length [m], bell radius [m]).
const GEOMETRIES: [(&str, &[(f64, f64)], (f64, f64)); 5] = [
    ("C4", &[(0.4494807268933825, 0.01175), (0.5279807268933824, 0.01175), (0.5615307268933825, 0.0118), (0.6115807268933825, 0.01485), (0.6524307268933826, 0.0177), (0.7294807268933826, 0.01425), (0.8672807268933825, 0.02005), (0.9234307268933825, 0.0192)], (1.0617307268933824, 0.0595)),
    ("E3", &[(0.8300807268933826, 0.01485), (0.8672807268933825, 0.02005), (0.9234307268933825, 0.0192)], (1.0617307268933824, 0.0595)),
    ("F3", &[(0.7294807268933826, 0.01425), (0.8672807268933825, 0.02005), (0.9234307268933825, 0.0192)], (1.0617307268933824, 0.0595)),
    ("Ab4", &[(0.2608057268933825, 0.001075), (0.6115807268933825, 0.01485), (0.6524307268933826, 0.0177), (0.7294807268933826, 0.01425), (0.8672807268933825, 0.02005), (0.9234307268933825, 0.0192)], (1.0617307268933824, 0.0595)),
    ("E5", &[(0.11888072689338253, 0.00085), (0.3235307268933826, 0.0048), (0.34843072689338256, 0.0097), (0.38178072689338255, 0.00755), (0.41663072689338254, 0.0099), (0.4494807268933825, 0.01175), (0.5279807268933824, 0.01175), (0.5615307268933825, 0.0118), (0.6115807268933825, 0.01485), (0.6524307268933826, 0.0177), (0.7294807268933826, 0.01425), (0.8672807268933825, 0.02005), (0.9234307268933825, 0.0192)], (1.0617307268933824, 0.0595)),
];

/// The exact resistance of two uniform circular pistons in an infinite rigid baffle, from the plane-wave spectrum (propagating part):
/// `R_nm = (rho c k^2/2 pi) Int_0^{pi/2} sin t D_n(k sin t) D_m(k sin t) J0(k d sin t) dt`, fixed 600-node Gauss-Legendre in `t`. Reference only.
fn exact_pair_resistance(gas: &GasState, omega: f64, a_n: f64, a_m: f64, d: f64) -> f64 {
    let k = omega / gas.sound_speed;
    let j0 = |x: f64| crate::modal::bessel_j0_certified(x);
    let dir = |x: f64| if x < 1.0e-9 { 1.0 } else { 2.0 * fs_math::special::bessel_j1(x) / x };
    let mut sum = 0.0;
    for (x, w) in crate::levine_schwinger::gauss_legendre(600) {
        let t = 0.5 * core::f64::consts::PI * x;
        let s = t.sin();
        sum += 0.5 * core::f64::consts::PI * w * s * dir(k * a_n * s) * dir(k * a_m * s) * j0(k * d * s);
    }
    gas.density * gas.sound_speed * k * k / (2.0 * core::f64::consts::PI) * sum
}

#[test]
fn the_piston_token_needs_the_half_space_and_changes_nothing_else() {
    assert!(!mono().is_piston());
    assert!(piston().is_piston());
    assert_ne!(mono(), piston());
    assert_eq!(piston().epsilon(), 2.0);
    for epsilon in [1.0, 3.0, 4.0, 2.0000001] {
        assert!(OpenHoleInteraction::new(epsilon).unwrap().with_piston().is_err(), "{epsilon}");
    }
    assert!(!OpenHoleInteraction::independent().is_piston());
}

#[test]
fn small_mouths_give_the_monopole_and_the_form_is_symmetric() {
    let gas = air();
    let (omega, d) = (TAU * 1500.0, 0.04);
    let monopole = mutual_radiation_impedance(&gas, omega, d, mono()).unwrap();
    let tiny = mutual_piston_impedance(&gas, omega, d, 1.0e-9, 1.0e-9, piston()).unwrap();
    assert!((tiny - monopole).abs() <= 1.0e-12 * monopole.abs(), "{tiny:?} vs {monopole:?}");
    let a = mutual_piston_impedance(&gas, omega, d, 0.01, 0.02, piston()).unwrap();
    let b = mutual_piston_impedance(&gas, omega, d, 0.02, 0.01, piston()).unwrap();
    assert_eq!(a, b);
    assert!(mutual_piston_impedance(&gas, omega, d, 0.0, 0.01, piston()).is_err());
    assert!(mutual_piston_impedance(&gas, omega, d, 0.01, f64::NAN, piston()).is_err());
}

#[test]
fn the_exact_pair_integral_reproduces_the_rayleigh_piston_and_the_monopole() {
    // The reference used below is itself checked: d = 0, equal radii is the Rayleigh self resistance, tiny radii the monopole of Eq. 16.
    let gas = air();
    for f in [1000.0, 8000.0, 20000.0] {
        let omega = TAU * f;
        let a = 0.01;
        let (r, _) = fs_phs::baffled_piston_impedance(gas.density, gas.sound_speed, a, omega).unwrap();
        let exact = exact_pair_resistance(&gas, omega, a, a, 0.0);
        assert!((exact - r).abs() <= 1.0e-8 * r, "{f} Hz: {exact} vs {r}");
        let monopole = mutual_radiation_impedance(&gas, omega, 0.04, mono()).unwrap().re;
        let small = exact_pair_resistance(&gas, omega, 1.0e-6, 1.0e-6, 0.04);
        assert!((small - monopole).abs() <= 1.0e-7 * monopole.abs().max(1.0), "{f} Hz: {small} vs {monopole}");
    }
}

#[test]
fn the_large_distance_form_is_the_far_limit_of_the_exact_pair_resistance() {
    let gas = air();
    for f in [1000.0, 20000.0] {
        let omega = TAU * f;
        let product = mutual_piston_impedance(&gas, omega, 3.0, 0.01, 0.02, piston()).unwrap().re;
        let exact = exact_pair_resistance(&gas, omega, 0.01, 0.02, 3.0);
        assert!((product - exact).abs() <= 2.0e-3 * exact.abs().max(product.abs()).max(1.0e-3), "{f} Hz: {product} vs {exact}");
    }
}

/// Resistance matrix of a page geometry at one frequency: the engine's own self load (Zorumski, `FLANGED_MODES` modes) on the diagonal, the chosen mutual form off it, distances by the engine's rule.
fn resistance_matrix(gas: &GasState, omega: f64, holes: &[(f64, f64)], bell: (f64, f64), use_piston: bool) -> Vec<f64> {
    let mut ports: Vec<CoupledPort> = holes
        .iter()
        .map(|&(x, a)| CoupledPort { position_m: x, radius_m: a, chimney: [C64::ONE, C64::ZERO, C64::ZERO, C64::ONE], self_radiation: C64::ZERO, mouth_radius_m: a })
        .collect();
    ports.push(CoupledPort::open_end(bell.0, C64::ZERO).with_mouth_radius(bell.1));
    let q = ports.len();
    let mut r = vec![0.0; q * q];
    for i in 0..q {
        let radius = if i + 1 == q { bell.1 } else { holes[i].1 };
        r[i * q + i] = flanged_pipe_impedance(gas, radius, omega, FLANGED_MODES).unwrap().re;
        for j in 0..i {
            let d = pair_distance(&ports[i], &ports[j]);
            let z = if use_piston {
                mutual_piston_impedance(gas, omega, d, ports[i].mouth_radius_m, ports[j].mouth_radius_m, piston()).unwrap()
            } else {
                mutual_radiation_impedance(gas, omega, d, mono()).unwrap()
            };
            r[i * q + j] = z.re;
            r[j * q + i] = z.re;
        }
    }
    r
}

/// Test (1) of the plan: is the radiation load passive? The smallest eigenvalue of the resistance matrix over the mean diagonal, at 941 frequencies (500 to 24000 Hz, 25 Hz), the engine's own
/// self load on the diagonal. Prints the counts of failing points for both forms.
#[test]
fn the_piston_form_radiation_resistance_is_psd_with_the_engine_self_load_on_the_five_page_geometries() {
    let gas = air();
    for (name, holes, bell) in GEOMETRIES {
        let q = holes.len() + 1;
        let (mut fail_piston, mut fail_mono, mut worst_piston, mut worst_mono, mut points) = (0usize, 0usize, 0.0f64, 0.0f64, 0usize);
        for step in 0..941 {
            let omega = TAU * (500.0 + 25.0 * step as f64);
            for use_piston in [true, false] {
                let r = resistance_matrix(&gas, omega, holes, bell, use_piston);
                let mean_diagonal = (0..q).map(|i| r[i * q + i]).sum::<f64>() / q as f64;
                let (eigenvalues, _) = fs_la::eigen::jacobi_eigh(&r, q);
                let lowest = eigenvalues[0] / mean_diagonal;
                if use_piston {
                    worst_piston = worst_piston.min(lowest);
                    fail_piston += usize::from(lowest < -1.0e-9);
                } else {
                    worst_mono = worst_mono.min(lowest);
                    fail_mono += usize::from(lowest < -1.0e-9);
                }
            }
            points += 1;
        }
        eprintln!("PSD {name}: {points} points; piston form fails {fail_piston} (worst {worst_piston:.3e}); monopole form fails {fail_mono} (worst {worst_mono:.3e})");
        assert_eq!(fail_piston, 0, "{name}: the piston form's load is not passive at {fail_piston} of {points} points (worst {worst_piston:e})");
        assert!(fail_mono > 0, "{name}: the monopole form is expected to be non-passive somewhere on the grid");
    }
}

/// The error of the large-distance form against the exact pair resistance on the page geometries, as an error of the coupling coefficient |R_form - R_exact| / sqrt(R_nn R_mm) with the exact
/// Rayleigh self resistance, over pairs at the engine's distances and 500 to 24000 Hz (50 Hz). Measured on the C4 geometry with the same 600 nodes in a script: p50 1.1e-3, max 4.5e-2.
#[test]
fn the_large_distance_form_stays_within_a_measured_coupling_error_of_the_exact_pair_resistance() {
    let gas = air();
    let (_, holes, bell) = GEOMETRIES[0];
    let mut ports: Vec<CoupledPort> = holes
        .iter()
        .map(|&(x, a)| CoupledPort { position_m: x, radius_m: a, chimney: [C64::ONE, C64::ZERO, C64::ZERO, C64::ONE], self_radiation: C64::ZERO, mouth_radius_m: a })
        .collect();
    ports.push(CoupledPort::open_end(bell.0, C64::ZERO).with_mouth_radius(bell.1));
    let mut worst = 0.0f64;
    let (mut sum, mut n) = (0.0f64, 0usize);
    for step in 0..471 {
        let omega = TAU * (500.0 + 50.0 * step as f64);
        let k = omega / gas.sound_speed;
        let self_r = |a: f64| {
            let (r, _) = fs_phs::baffled_piston_impedance(gas.density, gas.sound_speed, a, omega).unwrap();
            r
        };
        for i in 0..ports.len() {
            for j in 0..i {
                let d = pair_distance(&ports[i], &ports[j]);
                let (a, b) = (ports[i].mouth_radius_m, ports[j].mouth_radius_m);
                let exact = exact_pair_resistance(&gas, omega, a, b, d);
                let form = mutual_piston_impedance(&gas, omega, d, a, b, piston()).unwrap().re;
                let error = (form - exact).abs() / (self_r(a) * self_r(b)).sqrt();
                worst = worst.max(error);
                sum += error;
                n += 1;
                assert!(k > 0.0);
            }
        }
    }
    eprintln!("piston form coupling-coefficient error vs exact pair resistance: mean {:.3e}, max {worst:.3e} over {n} pair-points", sum / n as f64);
    assert!(worst < 0.06, "max coupling error {worst}");
}

#[test]
fn without_the_token_the_mouth_radius_changes_no_bit() {
    // `mouth_radius_m` is read by the piston form only: a card without the token couples exactly as before.
    let gas = air();
    let omega = TAU * 800.0;
    let plain = [CoupledPort::open_end(0.0, C64::new(2.0e4, -9.0e4)), CoupledPort::open_end(0.05, C64::new(3.0e4, -5.0e4))];
    let with_radius = [plain[0].with_mouth_radius(0.02), plain[1].with_mouth_radius(0.03)];
    let a = couple(&plain, &gas, omega, mono()).unwrap();
    let b = couple(&with_radius, &gas, omega, mono()).unwrap();
    for (x, y) in a.admittance.iter().zip(&b.admittance) {
        assert_eq!(x.re.to_bits(), y.re.to_bits());
        assert_eq!(x.im.to_bits(), y.im.to_bits());
    }
}

#[test]
fn the_piston_form_needs_every_mouth_radius_and_reduces_the_coupling_of_large_mouths() {
    let gas = air();
    let omega = TAU * 12000.0;
    let plain = [CoupledPort::open_end(0.0, C64::new(2.0e4, -9.0e4)), CoupledPort::open_end(0.05, C64::new(3.0e4, -5.0e4))];
    assert!(couple(&plain, &gas, omega, piston()).is_err(), "no mouth radius given: the piston form refuses");
    let ports = [plain[0].with_mouth_radius(0.02), plain[1].with_mouth_radius(0.0595)];
    assert!(couple(&ports, &gas, omega, piston()).is_ok());
    let monopole = mutual_radiation_impedance(&gas, omega, 0.05, mono()).unwrap();
    let piston_z = mutual_piston_impedance(&gas, omega, 0.05, 0.02, 0.0595, piston()).unwrap();
    assert!(piston_z.abs() < 0.2 * monopole.abs(), "{} vs {}", piston_z.abs(), monopole.abs());
}
