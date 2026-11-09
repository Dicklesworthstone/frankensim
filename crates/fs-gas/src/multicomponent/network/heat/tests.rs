use super::*;
use crate::multicomponent::tests::{gas, source};
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;

fn model(linear: f64, formation: f64) -> MixtureGas {
    gas(vec![
        source("A", R, 3.5, linear, formation),
        source("B", 2.0 * R, 4.5, 0.5 * linear, -2.0 * formation),
        source("absent-C", R, 3.5, 0.0, 7.0 * formation),
    ])
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * (1.0 + expected.abs()),
        "{actual:e} != {expected:e}"
    );
}
fn tank(gas: &MixtureGas) -> Plenum {
    Plenum::new(0.5, gas.conserved(&[0.4, 0.6, 0.0], 0.0, 1.4).unwrap()).unwrap()
}

#[test]
fn g1_plenum_heat_matches_original_two_body_implicit_limit_without_extra_work() {
    let gas = model(0.0, 30.0);
    let before = tank(&gas);
    for (wall_temperature, conductance) in [(5.0, 7.0), (0.5, 7.0), (5.0, 1e4)] {
        let dt = 0.4;
        let capacity = 3.0;
        let thermal = WallThermal::Lumped {
            temperature_k: wall_temperature,
            heat_capacity_j_per_k: capacity,
            conductance_w_per_k: conductance,
        };
        let mut after = before.clone();
        let pointer = after.state().partial_densities.as_ptr();
        let update = before
            .heat_candidate_into(
                &gas,
                Time::new(dt),
                thermal,
                WallControl::default(),
                &mut after,
                || false,
            )
            .unwrap();
        let gas_capacity = 0.5 * (0.4 * 2.5 + 0.6 * 0.5 * 3.5);
        let integrated = conductance * dt;
        let effective = integrated / (1.0 + integrated / capacity);
        let expected_t =
            (gas_capacity * 2.0 + effective * wall_temperature) / (gas_capacity + effective);
        let p = gas.primitive(after.state(), || false).unwrap();
        close(p.temperature, expected_t, 3e-9);
        close(p.pressure, 0.7 * expected_t, 3e-9);
        let expected_heat = gas_capacity * (expected_t - 2.0);
        close(update.report.heat_to_gas_j, expected_heat, 3e-9);
        close(update.report.wall_energy_change_j, -expected_heat, 3e-9);
        assert_eq!(update.report.external_heat_to_gas_j, 0.0);
        assert_eq!(update.report.mechanical_work_to_gas_j, 0.0);
        assert_eq!(update.report.wall_impulse_ns, 0.0);
        assert_eq!(
            after.state().partial_densities,
            before.state().partial_densities
        );
        assert_eq!(after.state().partial_densities.as_ptr(), pointer);
        assert_eq!(after.volume_m3(), before.volume_m3());
        close(
            (after.state().energy - before.state().energy) * 0.5
                + update.report.wall_energy_change_j,
            0.0,
            3e-9,
        );
    }
}

#[test]
fn g1_plenum_heat_retains_variable_cp_reference_gauge_and_prescribed_bath_account() {
    let mut reference = None;
    for formation in [0.0, -1e18, 1e18] {
        let gas = model(0.2, formation);
        let before = tank(&gas);
        let mut after = before.clone();
        let dt = 0.3;
        let conductance = 2.5;
        let wall_t = 5.0;
        let update = before
            .heat_candidate_into(
                &gas,
                Time::new(dt),
                WallThermal::Reservoir {
                    temperature_k: wall_t,
                    conductance_w_per_k: conductance,
                },
                WallControl::default(),
                &mut after,
                || false,
            )
            .unwrap();
        // Independent NASA-9 Cv integral with the two original molecular masses.
        let a = 0.5 * (0.4 * 2.5 + 0.6 * 0.5 * 3.5);
        let b = 0.5 * 0.5 * (0.4 * 0.2 + 0.6 * 0.5 * 0.1);
        let h = conductance * dt;
        let c = a * 2.0 + b * 4.0 + h * wall_t;
        let t = 2.0 * c / (a + h + ((a + h) * (a + h) + 4.0 * b * c).sqrt());
        let primitive = gas.primitive(after.state(), || false).unwrap();
        close(primitive.temperature, t, 2e-9);
        close(update.report.heat_to_gas_j, h * (wall_t - t), 2e-9);
        assert_eq!(update.report.wall_energy_change_j, 0.0);
        assert_eq!(
            update.report.external_heat_to_gas_j,
            update.report.heat_to_gas_j
        );
        assert_eq!(after.state().partial_densities[2], 0.0);
        let actual = [
            primitive.temperature,
            primitive.pressure,
            update.report.heat_to_gas_j,
        ];
        if let Some(expected) = reference {
            for (a, b) in actual.into_iter().zip(expected) {
                close(a, b, 3e-9);
            }
        } else {
            reference = Some(actual);
        }
    }
}

#[test]
fn g0_plenum_heat_refusal_and_last_poll_preserve_prepared_output() {
    let gas = model(0.2, 10.0);
    let before = tank(&gas);
    let saved = Plenum::new(0.8, gas.conserved(&[0.2, 0.8, 0.0], 0.0, 1.0).unwrap()).unwrap();
    let thermal = WallThermal::Reservoir {
        temperature_k: 5.0,
        conductance_w_per_k: 2.0,
    };
    let mut output = saved.clone();
    let mut count = 0;
    before
        .heat_candidate_into(
            &gas,
            Time::new(0.1),
            thermal,
            WallControl::default(),
            &mut output,
            || {
                count += 1;
                false
            },
        )
        .unwrap();
    output = saved.clone();
    let pointer = output.state().partial_densities.as_ptr();
    let mut seen = 0;
    let failure = before.heat_candidate_into(
        &gas,
        Time::new(0.1),
        thermal,
        WallControl::default(),
        &mut output,
        || {
            seen += 1;
            seen == count
        },
    );
    assert_eq!(failure, Err(WallExchangeError::Cancelled));
    assert_eq!(output, saved);
    assert_eq!(output.state().partial_densities.as_ptr(), pointer);
    assert!(
        before
            .heat_candidate_into(
                &gas,
                Time::new(0.1),
                WallThermal::Caloric {
                    mass_kg: 1.0,
                    specific_energy_j_kg: 1.0,
                    conductance_w_per_k: 2.0
                },
                WallControl::default(),
                &mut output,
                || false
            )
            .is_err()
    );
    assert_eq!(output, saved);
    // At extreme G*dt the represented endpoint temperatures cannot satisfy a
    // sub-nanojoule flux residual. Keep that original precision refusal.
    let stiff_gas = model(0.0, 30.0);
    let stiff_before = tank(&stiff_gas);
    let mut stiff_output = stiff_before.clone();
    assert!(matches!(
        stiff_before.heat_candidate_into(
            &stiff_gas,
            Time::new(0.4),
            WallThermal::Lumped {
                temperature_k: 5.0,
                heat_capacity_j_per_k: 3.0,
                conductance_w_per_k: 1e8,
            },
            WallControl::default(),
            &mut stiff_output,
            || false,
        ),
        Err(WallExchangeError::Balance { .. }),
    ));
    assert_eq!(stiff_output, stiff_before);
    assert_eq!(
        crate::multicomponent::Error::from(WallExchangeError::Cancelled),
        crate::multicomponent::Error::Cancelled
    );
    assert!(matches!(
        crate::multicomponent::Error::from(WallExchangeError::Invalid("heat")),
        crate::multicomponent::Error::Thermal(_)
    ));
}
