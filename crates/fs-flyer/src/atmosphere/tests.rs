//! G0/G3: hydrostatic pressure, common gas ownership, wind frames and real loads.

use super::*;
use crate::aircraft::wright_openloop_v1;
use fs_atmo::{DEC17_AIR, FlatSiteLogLaw, TurbulenceField};
use fs_material::gas::{ConductivityModel, GasSpec};

fn wind(speed: f64, turbulence: f64) -> Atmosphere {
    Atmosphere {
        mean: FlatSiteLogLaw {
            scenario_effective_z0_m: 0.003,
            displacement_height_m: 0.0,
            reference_height_m: 10.0,
            reference_speed_mps: speed,
        },
        turbulence: TurbulenceField::build(19, 8, turbulence, 20.0, speed).unwrap(),
        air: DEC17_AIR,
    }
}

fn close(a: f64, b: f64, relative: f64) {
    assert!((a - b).abs() <= relative * b.abs().max(1e-15), "{a} != {b}");
}

fn flight() -> AtmosphericFlightState {
    AtmosphericFlightState {
        position_m: [0.0, 0.0, 10.0],
        ground_velocity_mps: [13.0, 0.0],
        pitch_rad: 0.06,
        canard_rad: 0.1,
        propeller_rad_s: 45.0,
        pitch_rate_rad_s: 0.0,
    }
}

#[test]
fn profile_pressure_density_and_transport_are_the_original_material_state() {
    let profile = MaterialAirProfile::standard1976(1500.0).unwrap();
    let source = StandardAtmosphere1976::new();
    let wind = wind(6.0, 0.3);
    for h in [0.0, 10.0, 3500.0, 20000.0] {
        let original = source.at_geometric_altitude(1500.0 + h).unwrap();
        let sample = profile.sample(&wind, [13.0, -2.0, h], 240).unwrap();
        assert_eq!(sample.gas, original.gas);
        sample.gas.validate_coherence().unwrap();
        close(
            sample.gas.pressure,
            sample.gas.density * sample.gas.specific_gas_constant * sample.gas.temperature,
            3e-15,
        );
        let p_above = profile.gas_at_height(h + 0.1).unwrap().pressure;
        let p_below = source
            .at_geometric_altitude(1500.0 + h - 0.1)
            .unwrap()
            .gas
            .pressure;
        close(
            (p_below - p_above) / 0.2,
            sample.gas.density * original.reference_gravity_m_s2,
            2e-9,
        );
    }
    let low = profile.gas_at_height(0.0).unwrap();
    let high = profile.gas_at_height(3500.0).unwrap();
    assert!(high.pressure < low.pressure && high.density < low.density);
    assert!(high.sound_speed < low.sound_speed);
    assert!(high.dynamic_viscosity < low.dynamic_viscosity);
}

#[test]
fn both_wind_clocks_preserve_vectors_gradients_and_read_only_sampling() {
    let profile = MaterialAirProfile::standard1976(0.0).unwrap();
    let mut frozen = wind(6.0, 0.3);
    let [x, y, h] = [13.0, -2.0, 10.0];
    let original = frozen.sample_air_state(-x, y, h, 240).unwrap();
    let actual = profile.sample(&frozen, [x, y, h], 240).unwrap();
    for i in 0..3 {
        let sign = [-1.0, 1.0, 1.0];
        assert_eq!(actual.wind.u[i], sign[i] * original.velocity_mps[i]);
        for j in 0..3 {
            assert_eq!(
                actual.wind.grad[i][j],
                sign[i] * sign[j] * original.grad[i][j]
            );
        }
    }
    // Legacy constant air cannot leak into this source-bound pressure/density.
    frozen.air.pressure_pa = 1.0;
    frozen.air.rho_kg_m3 = 1e6;
    assert_eq!(actual, profile.sample(&frozen, [x, y, h], 240).unwrap());
    let mut evolving = frozen.into_evolving(71, vec![1.8; 8]).unwrap();
    evolving.turbulence.advance_to(240).unwrap();
    let before = evolving.clone();
    let a = profile.sample_evolving(&evolving, [x, y, h]).unwrap();
    let b = profile.sample_evolving(&evolving, [x, y, h]).unwrap();
    assert_eq!(a, b);
    assert_eq!(evolving, before);
    let native = evolving.sample_air_state(-x, y, h).unwrap();
    assert_eq!(
        a.wind.u,
        [
            -native.velocity_mps[0],
            native.velocity_mps[1],
            native.velocity_mps[2]
        ]
    );
    assert_eq!(a.gas, actual.gas);
}

#[test]
fn arbitrary_uniform_gas_preserves_composition_transport_and_pressure_scaling() {
    let spec = GasSpec {
        molar_mass: 0.039948,
        gamma: 5.0 / 3.0,
        conductivity: ConductivityModel::Eucken,
        ..GasSpec::dry_air_ussa1976()
    };
    let low = GasState::try_new(&spec, 300.0, 60000.0).unwrap();
    let high = GasState::try_new(&spec, 300.0, 120000.0).unwrap();
    let wind = wind(4.0, 0.0);
    let a = MaterialAirProfile::uniform(low)
        .unwrap()
        .sample(&wind, [0.0, 0.0, 10.0], 0)
        .unwrap();
    let b = MaterialAirProfile::uniform(high)
        .unwrap()
        .sample(&wind, [0.0, 0.0, 10.0], 0)
        .unwrap();
    assert_eq!(a.gas, low);
    assert_eq!(b.gas, high);
    close(
        b.dynamic_pressure_pa([13.0, 0.0, 0.0]).unwrap()
            / a.dynamic_pressure_pa([13.0, 0.0, 0.0]).unwrap(),
        2.0,
        1e-15,
    );
    close(
        b.reynolds(1.5, [13.0, 0.0, 0.0]).unwrap() / a.reynolds(1.5, [13.0, 0.0, 0.0]).unwrap(),
        2.0,
        1e-15,
    );
    assert_eq!(low.sound_speed, high.sound_speed);
    assert_eq!(low.dynamic_viscosity, high.dynamic_viscosity);
}

#[test]
fn hydrostatic_density_drives_actual_wing_and_propeller_loads() {
    let design = wright_openloop_v1();
    let profile = MaterialAirProfile::standard1976(0.0).unwrap();
    let still = wind(0.0, 0.0);
    let low_state = flight();
    let high_state = AtmosphericFlightState {
        position_m: [0.0, 0.0, 5000.0],
        ..low_state
    };
    let low = design
        .force_buildup_in_atmosphere(&profile, &still, low_state, 0)
        .unwrap();
    let high = design
        .force_buildup_in_atmosphere(&profile, &still, high_state, 0)
        .unwrap();
    for sample in [&low, &high] {
        assert_eq!(
            sample.loads,
            design
                .force_buildup(
                    sample.airspeed_mps,
                    sample.angle_of_attack_rad,
                    low_state.canard_rad,
                    low_state.propeller_rad_s,
                    low_state.pitch_rate_rad_s,
                    sample.air.gas.density,
                )
                .unwrap()
        );
    }
    let ratio = high.air.gas.density / low.air.gas.density;
    assert!(ratio < 0.7);
    for (a, b) in [
        (low.loads.lift_n, high.loads.lift_n),
        (low.loads.thrust_n[0], high.loads.thrust_n[0]),
        (low.loads.thrust_n[1], high.loads.thrust_n[1]),
    ] {
        close(b / a, ratio, 2e-8);
    }
    assert_ne!(
        low.loads.torque_imbalance_nm,
        high.loads.torque_imbalance_nm
    );
}

#[test]
fn force_path_uses_relative_wind_once_and_evolving_queries_do_not_advance() {
    let design = wright_openloop_v1();
    let profile = MaterialAirProfile::standard1976(0.0).unwrap();
    let state = flight();
    let frozen = wind(4.0, 0.0);
    let loads = design
        .force_buildup_in_atmosphere(&profile, &frozen, state, 240)
        .unwrap();
    close(loads.airspeed_mps, 17.0, 2e-15);
    close(loads.angle_of_attack_rad, state.pitch_rad, 2e-15);
    let mut evolving = frozen.into_evolving(71, vec![1.8; 8]).unwrap();
    evolving.turbulence.advance_to(240).unwrap();
    let before = evolving.clone();
    let dynamic = design
        .force_buildup_in_evolving_atmosphere(&profile, &evolving, state)
        .unwrap();
    assert_eq!(loads.air, dynamic.air);
    assert_eq!(loads.loads, dynamic.loads);
    assert_eq!(before, evolving);
}

#[test]
fn pressure_profile_and_wind_domain_refusals_are_not_replaced_by_constants() {
    for altitude in [-1.0, 80000.1, f64::INFINITY, f64::NAN] {
        assert!(MaterialAirProfile::standard1976(altitude).is_err());
    }
    let profile = MaterialAirProfile::standard1976(1000.0).unwrap();
    assert!(profile.gas_at_height(79000.0).is_ok());
    for height in [-1.0, 79000.1, f64::MAX, f64::INFINITY, f64::NAN] {
        assert!(profile.gas_at_height(height).is_err());
    }
    let mut bad = profile.gas_at_height(0.0).unwrap();
    bad.pressure *= 2.0;
    assert!(MaterialAirProfile::uniform(bad).is_err());
    let wind = wind(6.0, 0.3);
    assert!(profile.sample(&wind, [f64::NAN, 0.0, 10.0], 0).is_err());
    let design = wright_openloop_v1();
    let invalid_state = AtmosphericFlightState {
        ground_velocity_mps: [f64::INFINITY, 0.0],
        ..flight()
    };
    assert!(
        design
            .force_buildup_in_atmosphere(&profile, &wind, invalid_state, 0)
            .is_err()
    );
}
