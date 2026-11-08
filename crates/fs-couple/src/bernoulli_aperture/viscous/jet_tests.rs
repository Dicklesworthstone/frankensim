//! Tests of the separated-jet slit law (van Zon et al. 1990) as a [`ViscousSlit`]: the public entry points follow it, the spatial and
//! condensed paths refuse it, and every other constructor is still the fully developed law. The pure pressure law, its constants
//! and the Fig. 2 theory curve are tested in `separated_jet`; here the reference law is written out again without its series,
//! Newton solve or tangent formulas (bisection on the printed closed form).

use super::*;
use fs_material::gas::GasSpec;

const WIDTH: f64 = 0.013;

fn gas() -> GasState { GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap() }

fn jet(length: f64) -> ViscousSlit { ViscousSlit::separated_jet_from_gas(&gas(), length).unwrap() }

fn close(actual: f64, expected: f64, relative: f64) {
    assert!((actual - expected).abs() <= relative * expected.abs().max(1e-250), "actual={actual:e}, expected={expected:e}");
}

/// van Zon 1990 written out: eq. 3 below L/h = 3; eqs. 4-5 by bisection while the boundary layers are growing; eq. 1 per length
/// past l_c after the Bernoulli core (eq. 5 at d_c), with eq. 6 for l_c and S = 2 h.
fn reference_drop(speed: f64, h: f64, length: f64, rho: f64, mu: f64) -> f64 {
    let alpha = 0.45;
    if length < 3.0 * h { return rho * speed * speed / (2.0 * alpha * alpha); }
    let x = length / h - 2.0;
    let reynolds = speed * h / (mu / rho);
    let d_c = (4.0 / 9.0) * (1.0 - (5.0_f64 / 32.0).sqrt());
    let braces = |d: f64| 4.0 * d + 9.0 * (1.0 - d).ln() + 5.0 * d / (1.0 - d);
    let c = braces(d_c) / 6.0;
    if c * reynolds >= x {
        let y = 6.0 * x / reynolds;
        let (mut lo, mut hi) = (0.0, d_c);
        for _ in 0..200 {
            let mid = 0.5 * (lo + hi);
            if braces(mid) < y { lo = mid; } else { hi = mid; }
        }
        let d = 0.5 * (lo + hi);
        rho * speed * speed / (2.0 * (1.0 - d) * (1.0 - d))
    } else {
        let l_c = 2.0 * h + c * reynolds * h;
        rho * speed * speed / (2.0 * (1.0 - d_c) * (1.0 - d_c)) + 12.0 * mu * (length - l_c) * speed / (h * h)
    }
}

// Geometries: contracted (L/h < 3), entrance at card speeds, and long channels, at the rest opening 0.4 mm of the reed cards.
const GEOMETRIES: [(f64, f64); 6] = [(6e-4, 1.4e-3), (4e-4, 1.4e-3), (2e-4, 1.4e-3), (4e-5, 1.4e-3), (3e-4, 25e-3), (4e-5, 25e-3)];

#[test]
fn only_the_new_constructor_builds_the_jet_law_and_it_survives_a_gas_rebind() {
    let g = gas();
    assert!(!ViscousSlit::new(1.81e-5, 1.4e-3).unwrap().is_separated_jet());
    assert!(!ViscousSlit::from_gas(&g, 1.4e-3).unwrap().is_separated_jet());
    let law = jet(1.4e-3);
    assert!(law.is_separated_jet());
    assert_eq!(law.gas(), Some(&g));
    assert_eq!(law.dynamic_viscosity_pa_s(), g.dynamic_viscosity);
    assert_eq!(law.length_m(), 1.4e-3);
    assert_ne!(law, ViscousSlit::from_gas(&g, 1.4e-3).unwrap());
    let hot = GasState::try_new(&GasSpec::dry_air_ussa1976(), 313.15, 101325.0).unwrap();
    let rebound = law.rebind_source_gas(&hot).unwrap();
    assert!(rebound.is_separated_jet());
    assert_eq!(rebound.gas(), Some(&hot));
    assert!(!ViscousSlit::from_gas(&g, 1.4e-3).unwrap().rebind_source_gas(&hot).unwrap().is_separated_jet());
    // a density from another gas is refused, as for the fully developed law
    assert!(law.evaluate(WIDTH, 2e-4, hot.density, 100.0).is_err());
}

#[test]
fn evaluate_solves_the_printed_pressure_law_signed_with_a_closing_power_balance() {
    let g = gas();
    let (rho, mu) = (g.density, g.dynamic_viscosity);
    for (h, length) in GEOMETRIES {
        let law = jet(length);
        for dp in [1e-6, 0.3, 40.0, 900.0, 6000.0, 15000.0] {
            let forward = law.evaluate(WIDTH, h, rho, dp).unwrap();
            let reverse = law.evaluate(WIDTH, h, rho, -dp).unwrap();
            assert!(forward.flow_m3_s > 0.0);
            assert_eq!(reverse.flow_m3_s, -forward.flow_m3_s);
            let speed = forward.flow_m3_s / (WIDTH * h);
            close(reference_drop(speed, h, length, rho, mu), dp, 2e-8);
            close(forward.dissipated_power_w(), dp * forward.flow_m3_s, 1e-12);
            assert!(forward.viscous_dissipation_w >= 0.0 && forward.separation_dissipation_w > 0.0);
            close(forward.separation_dissipation_w, 0.5 * rho * speed * speed * forward.flow_m3_s, 1e-12);
            assert_eq!(reverse.viscous_dissipation_w, forward.viscous_dissipation_w);
            assert!(forward.flow_per_pressure_m3_s_pa > 0.0);
        }
    }
}

#[test]
fn the_pressure_and_gap_tangents_match_independent_differences_on_every_branch() {
    let g = gas();
    let rho = g.density;
    for (h, length) in GEOMETRIES {
        let law = jet(length);
        for dp in [0.5, 60.0, 2500.0, 9000.0] {
            let value = law.evaluate(WIDTH, h, rho, dp).unwrap();
            let (ep, eh) = (dp * 1e-6, h * 1e-6);
            let flow = |gap: f64, pressure: f64| law.evaluate(WIDTH, gap, rho, pressure).unwrap().flow_m3_s;
            close(value.flow_per_pressure_m3_s_pa, (flow(h, dp + ep) - flow(h, dp - ep)) / (2.0 * ep), 2e-5);
            close(value.flow_per_opening_m2_s, (flow(h + eh, dp) - flow(h - eh, dp)) / (2.0 * eh), 1e-4);
        }
    }
}

#[test]
fn below_l_over_h_3_the_flow_is_alpha_times_the_bernoulli_flow_and_the_rest_tangent_is_infinite() {
    let g = gas();
    let (rho, h, length) = (g.density, 6e-4, 1.4e-3);
    let law = jet(length);
    for dp in [0.01, 100.0, 5000.0] {
        let flow = law.evaluate(WIDTH, h, rho, dp).unwrap();
        close(flow.flow_m3_s, 0.45 * WIDTH * h * (2.0 * dp / rho).sqrt(), 1e-12);
    }
    let rest = law.evaluate(WIDTH, h, rho, 0.0).unwrap();
    assert_eq!((rest.flow_m3_s, rest.flow_per_opening_m2_s), (0.0, 0.0));
    assert!(rest.flow_per_pressure_m3_s_pa.is_infinite());
    // reattached at the same length: finite Poiseuille tangent at rest
    let attached = law.evaluate(WIDTH, 4e-4, rho, 0.0).unwrap();
    close(attached.flow_per_pressure_m3_s_pa, WIDTH * 4e-4 * 4e-4 * 4e-4 / (12.0 * g.dynamic_viscosity * (1.4e-3 - 8e-4)), 1e-12);
}

#[test]
fn at_the_card_length_the_jet_law_carries_the_entrance_loss_the_fully_developed_law_lacks() {
    let g = gas();
    let (rho, mu) = (g.density, g.dynamic_viscosity);
    let old = ViscousSlit::from_gas(&g, 25e-3).unwrap();
    let new = jet(25e-3);
    // Low drive: Poiseuille over L - S = L - 2 h instead of L, so MORE flow, by L / (L - 2 h).
    let (h, dp) = (3e-4, 1e-3);
    let ratio = new.evaluate(WIDTH, h, rho, dp).unwrap().flow_m3_s / old.evaluate(WIDTH, h, rho, dp).unwrap().flow_m3_s;
    close(ratio, 25e-3 / (25e-3 - 2.0 * h), 2e-6);
    // Card drive (5 kPa across 0.3 mm): the quadratic coefficient 0.7438 against 0.5 takes flow away.
    let (jet_flow, old_flow) = (new.evaluate(WIDTH, h, rho, 5000.0).unwrap().flow_m3_s, old.evaluate(WIDTH, h, rho, 5000.0).unwrap().flow_m3_s);
    assert!(jet_flow < old_flow, "{jet_flow} {old_flow}");
    // and the drop is the printed closed form on the Poiseuille branch: 0.7438 rho v^2 + 12 mu (L - 2 h) v / h^2
    let v = jet_flow / (WIDTH * h);
    close(0.743_847_343_454_391_3 * rho * v * v + 12.0 * mu * (25e-3 - 2.0 * h) * v / (h * h), 5000.0, 1e-9);
}

#[test]
fn the_characteristic_junction_recovers_manufactured_bidirectional_flow() {
    let g = gas();
    let (rho, mu) = (g.density, g.dynamic_viscosity);
    let impedance = 1e6;
    for (h, length) in GEOMETRIES {
        let law = jet(length);
        for velocity in [-30.0_f64, -2.0, -1e-4, 1e-4, 2.0, 30.0] {
            let dp = velocity.signum() * reference_drop(velocity.abs(), h, length, rho, mu);
            for incoming in [-200.0, 0.0, 200.0] {
                for moving in [-1e-6, 0.0, 1e-6] {
                    let expected = incoming + impedance * (WIDTH * h * velocity + moving);
                    let upstream = dp + expected + incoming;
                    let actual = law.characteristic_pressure(WIDTH, h, rho, impedance, incoming, upstream, moving).unwrap();
                    assert!((actual - expected).abs() <= 3e-9 * (1.0 + expected.abs()), "h={h} L={length} v={velocity}: {actual} vs {expected}");
                }
            }
        }
        // a shut slit and a zero drive return the no-jet wave, as for the fully developed law
        assert_eq!(law.characteristic_pressure(WIDTH, 0.0, rho, impedance, 75.0, 1000.0, 2e-7).unwrap(), 75.2);
        // (the drive is upstream - 2 incoming - Z moving = 20 - 2 * 10 - 0)
        assert_eq!(law.characteristic_pressure(WIDTH, h, rho, impedance, 10.0, 20.0, 0.0).unwrap(), 10.0);
    }
}

#[test]
fn staple_recovery_closes_on_the_jet_law_and_a_vanishing_recovery_is_the_plain_junction() {
    let g = gas();
    let (rho, mu) = (g.density, g.dynamic_viscosity);
    let recovery = StapleRecovery::new(0.8, 4e-6).unwrap();
    let (width, impedance) = (0.0056, 6e7);
    let (mut forward, mut reverse) = (0, 0);
    for (h, length) in [(3e-4, 1.4e-3), (1.5e-4, 1.4e-3), (6e-4, 1.4e-3), (4e-5, 25e-3)] {
        let law = jet(length);
        let area = width * h;
        if recovery.forward_factor(area).is_err() { continue; }
        for incoming in [-300.0, 0.0, 1000.0] {
            for mouth in [-4000.0, -50.0, 50.0, 800.0, 7400.0, 30000.0] {
                for moving in [-8e-5, 0.0, 6e-5] {
                    let j = law.characteristic_junction_with_recovery(recovery, width, h, rho, impedance, incoming, mouth, moving).unwrap();
                    let drive = mouth - 2.0 * incoming - impedance * moving;
                    let velocity = j.flow_m3_s / area;
                    // The junction equation: slit drop - recovered head + Z q = drive, with the slit's own drop from the printed law.
                    let reference = velocity.signum() * reference_drop(velocity.abs(), h, length, rho, mu);
                    if j.flow_m3_s != 0.0 {
                        close(j.constriction_drop_pa, reference, 2e-8);
                    }
                    let balance = (j.constriction_drop_pa - j.recovered_pa) + impedance * j.flow_m3_s;
                    assert!((balance - drive).abs() <= 3e-9 * (drive.abs() + balance.abs() + 1.0), "{balance} vs {drive}");
                    if j.flow_m3_s > 0.0 {
                        forward += 1;
                        assert!(j.recovered_pa > 0.0 && j.recovered_pa < j.constriction_drop_pa);
                    } else if j.flow_m3_s < 0.0 {
                        reverse += 1;
                        assert_eq!(j.recovered_pa, 0.0);
                    }
                }
            }
        }
    }
    assert!(forward > 40 && reverse > 40, "forward {forward}, reverse {reverse}");
    // C_P -> 0 is the plain jet-law junction
    let tiny = StapleRecovery::new(1e-14, 4e-6).unwrap();
    for (h, length, mouth, incoming, moving) in [(3e-4, 1.4e-3, 7400.0, 75.0, 1e-7), (1e-4, 1.4e-3, 800.0, -30.0, -5e-8), (6e-4, 1.4e-3, -2000.0, 10.0, 0.0)] {
        let law = jet(length);
        let j = law.characteristic_junction_with_recovery(tiny, width, h, rho, impedance, incoming, mouth, moving).unwrap();
        let plain = law.characteristic_pressure(width, h, rho, impedance, incoming, mouth, moving).unwrap();
        assert!((j.outgoing_pa - plain).abs() <= 1e-9 * (1.0 + plain.abs()), "{} vs {plain}", j.outgoing_pa);
    }
    // a shut slit hands the blade the whole difference, as for the fully developed law
    let shut = jet(1.4e-3).characteristic_junction_with_recovery(recovery, 0.0056, 0.0, rho, 6e7, 75.0, 7400.0, 1e-7).unwrap();
    assert_eq!((shut.flow_m3_s, shut.recovered_pa), (0.0, 0.0));
}

#[test]
fn the_spatial_and_condensed_paths_refuse_the_jet_law_and_keep_the_fully_developed_one() {
    let g = gas();
    let rho = g.density;
    let law = jet(1.4e-3);
    assert!(law.linear_gap_flow(WIDTH, [1e-4, 3e-4], rho, 100.0).is_err());
    assert!(law.fixed_opening_pressure_drop(WIDTH, 2e-4, rho, 1e6, 100.0).is_err());
    let old = ViscousSlit::from_gas(&g, 1.4e-3).unwrap();
    assert!(old.linear_gap_flow(WIDTH, [1e-4, 3e-4], rho, 100.0).is_ok());
    assert!(old.fixed_opening_pressure_drop(WIDTH, 2e-4, rho, 1e6, 100.0).is_ok());
}

#[test]
fn invalid_and_unrepresentable_inputs_are_refused_and_a_shut_slit_is_still_shut() {
    let g = gas();
    let rho = g.density;
    let law = jet(1.4e-3);
    assert_eq!(law.evaluate(WIDTH, 0.0, rho, 500.0).unwrap(), ViscousSlitFlow::default());
    for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(law.evaluate(bad, 2e-4, rho, 100.0).is_err());
        assert!(law.evaluate(WIDTH, 2e-4, bad, 100.0).is_err());
        assert!(law.characteristic_pressure(WIDTH, 2e-4, rho, bad, 0.0, 100.0, 0.0).is_err());
        assert!(ViscousSlit::separated_jet_from_gas(&g, bad).is_err());
    }
    assert!(law.evaluate(WIDTH, 2e-4, rho, f64::NAN).is_err());
    assert!(law.evaluate(WIDTH, -1e-4, rho, 100.0).is_err());
    for dp in [1e-20, 1e-100] {
        let forward = law.evaluate(WIDTH, 2e-4, rho, dp).unwrap();
        let reverse = law.evaluate(WIDTH, 2e-4, rho, -dp).unwrap();
        assert!(forward.flow_m3_s > 0.0);
        assert_eq!(reverse.flow_m3_s, -forward.flow_m3_s);
    }
}
