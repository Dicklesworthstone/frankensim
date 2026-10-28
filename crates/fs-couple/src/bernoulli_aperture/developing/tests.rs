use super::*;
use crate::bernoulli_aperture::viscous::ViscousSlit;
use fs_material::gas::GasSpec;

fn air() -> GasState { GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap() }

#[test]
fn g1_poiseuille_limit_is_the_fully_developed_slit_plus_shahs_incremental_drop() {
    // Gate 1 (ruling 8736), ladder-free: as x+ grows the held relation must become the
    // fully developed owner's 12 mu L U/h^2 + rho U^2/2 plus Shah's 0.674 rho U^2/2. The
    // correlation's blend leaves f_app Re - 24 - 0.674/(4x+) ~ -6.96e-4/x+^2, so the speed
    // error is below 4e-5/x+^2; leaving out the 0.674 term would cost 7e-4 at x+ = 10.
    let gas = air();
    let (width, h, speed) = (0.01, 1e-3, 1.0);
    let nu = gas.dynamic_viscosity / gas.density;
    let reynolds = speed * 2.0 * h / nu;
    for x_plus in [10.0, 100.0, 1000.0] {
        let length = x_plus * 2.0 * h * reynolds;
        let head = 0.5 * gas.density * speed * speed;
        let developed = 12.0 * gas.dynamic_viscosity * length * speed / (h * h) + head;
        let pressure = developed + SHAH_K_INFINITY * head;
        let fully = ViscousSlit::from_gas(&gas, length).unwrap()
            .evaluate(width, h, gas.density, developed).unwrap().flow_m3_s / (width * h);
        assert!((fully - speed).abs() <= 1e-12, "the owners' fully developed relations must agree: {fully}");
        let held = DevelopingSlit::from_gas(&gas, length).unwrap().held(width, h, pressure).unwrap();
        let error = (held.bulk_speed_m_s - speed).abs();
        assert!(error <= 4e-5 / (x_plus * x_plus), "x+ {x_plus}: error {error:e}");
        let without = ViscousSlit::from_gas(&gas, length).unwrap()
            .evaluate(width, h, gas.density, pressure).unwrap().flow_m3_s / (width * h);
        assert!((without - speed).abs() > 100.0 * 4e-5 / (x_plus * x_plus),
            "x+ {x_plus}: the incremental drop must be resolved, not hidden in the tolerance");
        assert!((held.x_plus / x_plus - 1.0).abs() < 1e-6 && held.developed);
    }
}

#[test]
fn g1_correlation_reaches_both_published_asymptotes() {
    // Short duct: 3.44/sqrt(x+) (the developing-layer limit); long duct: 24 + 0.674/(4x+).
    // The blend's remainder there is about 0.674 x+/(4 * 0.000029), i.e. 1690 x+^1.5 of the
    // leading term.
    for x_plus in [1e-8, 1e-7, 1e-6] {
        let short = 3.44 / f64::sqrt(x_plus);
        assert!((shah_apparent_friction_re(x_plus) / short - 1.0).abs() < 2000.0 * x_plus.powf(1.5));
    }
    for x_plus in [1.0, 10.0, 100.0] {
        let long = 24.0 + 0.674 / (4.0 * x_plus);
        assert!((shah_apparent_friction_re(x_plus) - long).abs() <= 7e-4 / (x_plus * x_plus));
    }
}

#[test]
fn g2_held_speed_solves_the_relation_below_the_inviscid_core() {
    let gas = air();
    let law = DevelopingSlit::from_gas(&gas, 0.02455).unwrap();
    let mut previous = 0.0;
    for pressure in [1e-3, 1.0, 67.15, 757.5, 5000.0] {
        let held = law.held(0.008, 0.00148, pressure).unwrap();
        let residual = law.pressure_at(0.00148, held.bulk_speed_m_s) - pressure;
        assert!(residual.abs() <= 1e-12 * pressure, "{pressure} Pa: residual {residual:e}");
        assert!(held.bulk_speed_m_s > previous && held.bulk_speed_m_s < held.core_speed_m_s);
        assert!((held.flow_m3_s - held.bulk_speed_m_s * 0.008 * 0.00148).abs() <= 1e-15 * held.flow_m3_s);
        previous = held.bulk_speed_m_s;
    }
    // The prediction mail's whistle D5 card (8740): 67.15 Pa holds 7.198 m/s, 17.5 % below the
    // fully developed law's 8.728 m/s, with the core still inviscid at the exit.
    let held = law.held(0.008, 0.00148, 67.15).unwrap();
    assert!((held.bulk_speed_m_s - 7.198).abs() < 2e-3 && !held.developed, "{held:?}");
    assert_eq!(law.held(0.008, 0.00148, 0.0).unwrap().bulk_speed_m_s, 0.0);
}

#[test]
fn g3_gas_and_geometry_are_explicit_and_refused_when_unphysical() {
    let gas = air();
    let hot = GasState::try_new(&GasSpec::dry_air_ussa1976(), 393.15, 101325.0).unwrap();
    let speed = |gas: &GasState| DevelopingSlit::from_gas(gas, 0.03).unwrap()
        .held(0.008, 0.0015, 100.0).unwrap().bulk_speed_m_s;
    assert_ne!(speed(&gas), speed(&hot), "density and viscosity must come from the gas");
    let mut corrupt = gas; corrupt.dynamic_viscosity *= 2.0;
    assert!(DevelopingSlit::from_gas(&corrupt, 0.03).is_err());
    for length in [0.0, -1.0, f64::NAN] { assert!(DevelopingSlit::from_gas(&gas, length).is_err()); }
    let law = DevelopingSlit::from_gas(&gas, 0.03).unwrap();
    assert!(law.held(0.0, 0.0015, 100.0).is_err());
    assert!(law.held(0.008, -0.0015, 100.0).is_err());
    assert!(law.held(0.008, 0.0015, -1.0).is_err());
    assert!(law.held(0.008, 0.0015, f64::INFINITY).is_err());
}

#[test]
fn g4_validity_refuses_short_ducts_and_low_reynolds_numbers() {
    let gas = air();
    // Whistle D5 under W1 (93.5 Pa holds the measured 6.2 l/min): inside every bound.
    let whistle = DevelopingSlit::from_gas(&gas, 0.02455).unwrap();
    let held = whistle.held(0.008, 0.00148, 93.5).unwrap();
    assert!(held.validity().is_ok(), "{held:?}");
    // A 1 mm slit at 5 kPa: x+ below 4e-4.
    let short = DevelopingSlit::from_gas(&gas, 0.001).unwrap().held(0.008, 0.00148, 5000.0).unwrap();
    assert!(short.x_plus < MINIMUM_X_PLUS);
    assert!(short.validity().unwrap_err().contains("x+"));
    // A 1 mm square-section slit at 0.1 Pa: Re about 30, where the inlet head is a third of
    // the drop and Chen's 38/Re moves it by about 38 %. (In a long slit at low Re the viscous
    // drop dominates and the inlet term cannot move it by 2.4 %.)
    let slow = DevelopingSlit::from_gas(&gas, 0.001).unwrap().held(0.008, 0.001, 0.1).unwrap();
    assert!(slow.x_plus >= MINIMUM_X_PLUS && slow.chen_inlet_shift > STATED_ACCURACY, "{slow:?}");
    assert!(slow.validity().unwrap_err().contains("Reynolds number too low"));
}

#[test]
fn g4_validity_refuses_unstable_layers_and_channels_that_can_sustain_turbulence() {
    let gas = air();
    // The whistle windway at 3 kPa keeps its core with Re_delta* near 505; at 4 kPa the
    // layer passes 520 (near 548) while x+ and Chen's term stay inside their bounds.
    let whistle = DevelopingSlit::from_gas(&gas, 0.02455).unwrap();
    let below = whistle.held(0.008, 0.00148, 3000.0).unwrap();
    assert!(!below.developed && below.reynolds_displacement < LAYER_CRITICAL_REYNOLDS_DISPLACEMENT);
    assert!(below.validity().is_ok(), "{below:?}");
    let above = whistle.held(0.008, 0.00148, 4000.0).unwrap();
    assert!(!above.developed && above.x_plus >= MINIMUM_X_PLUS && above.chen_inlet_shift <= STATED_ACCURACY);
    assert!(above.validity().unwrap_err().contains("Re_delta* 520"), "{above:?}");
    // A 60 mm, 1 mm channel is developed at the exit: 300 Pa holds Re_Dh near 1600 and
    // is admitted; 500 Pa holds about 2245, above 1760, and is refused even though its
    // Re_delta* estimate (near 395) would pass the layer bound.
    let channel = DevelopingSlit::from_gas(&gas, 0.06).unwrap();
    let laminar = channel.held(0.008, 0.001, 300.0).unwrap();
    assert!(laminar.developed && laminar.validity().is_ok(), "{laminar:?}");
    let sustained = channel.held(0.008, 0.001, 500.0).unwrap();
    assert!(sustained.developed && sustained.reynolds_hydraulic > CHANNEL_CRITICAL_REYNOLDS_HYDRAULIC);
    assert!(sustained.reynolds_displacement < LAYER_CRITICAL_REYNOLDS_DISPLACEMENT);
    assert!(sustained.validity().unwrap_err().contains("Re_Dh 1760"), "{sustained:?}");
    assert_eq!(CHANNEL_CRITICAL_REYNOLDS_HYDRAULIC, 1760.0);
}
