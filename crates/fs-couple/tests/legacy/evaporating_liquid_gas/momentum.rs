//! Actual native/liquid and original gas consumers of retained mean drift.

use super::*;
use fs_couple::liquid_gas::evaporation::momentum::{
    MovingGasEvaporationCandidate, MovingGasParcel, advance_liquid_with_moving_gas,
};
use fs_material::evaporation::gas::GasEvaporationError;
use std::sync::atomic::{AtomicUsize, Ordering};

fn liquid_velocity(liquid: &MaterialFreeSurfaceEnthalpy3<'_, '_>) -> [f64; 3] {
    liquid
        .flow()
        .sample_at_si(INTERFACE[0], INTERFACE[1], INTERFACE[2])
        .unwrap()
        .unwrap()
        .velocity_m_s
}

fn check_source_and_thermal(
    liquid: &MaterialFreeSurfaceEnthalpy3<'_, '_>,
    gas: &MovingGasParcel,
    binding: &GasLiquidBinding,
    next: &MovingGasEvaporationCandidate<'_, '_, '_>,
) {
    let report = next.report();
    let native = &report.native;
    let motion = &report.momentum;
    for d in 0..3 {
        near(
            native.mass_exchange.momentum_change_kg_m_s[d] + next.gas().momentum_kg_m_s()[d]
                - gas.momentum_kg_m_s()[d],
            0.0,
            controls().native.momentum_tolerance_kg_m_s,
        );
    }
    // Independent actual caloric inventory plus SOURCE kinetic changes. The
    // following VOF step's mechanical truncation is not relabeled as heat.
    let liquid_mass =
        next.liquid().flow().ledger_mass_kg().unwrap() - liquid.flow().ledger_mass_kg().unwrap();
    let enthalpy = next.liquid().total_enthalpy_j().unwrap() - liquid.total_enthalpy_j().unwrap();
    let before = gas.parcel();
    let after = next.gas().parcel();
    let internal = after.mass_kg() * after.state().specific_internal_energy_j_kg()
        - before.mass_kg() * before.state().specific_internal_energy_j_kg();
    let gas_kinetic = next.gas().kinetic_energy_j() - gas.kinetic_energy_j();
    near(gas_kinetic, motion.gas_kinetic_energy_change_j, 1e-17);
    near(
        enthalpy
            + binding.liquid_enthalpy_offset_j_kg * liquid_mass
            + internal
            + motion.liquid_kinetic_energy_change_j
            + gas_kinetic,
        native.flow.external_heat_j + native.flow.viscous_heat_j + native.flow.body_heat_j,
        controls().energy_tolerance_j,
    );
    near(
        liquid_mass + after.mass_kg() - before.mass_kg(),
        0.0,
        controls().mass_tolerance_kg,
    );
    near(native.mass_residual_kg, 0.0, controls().mass_tolerance_kg);
    near(
        native.volatile_mass_residual_kg,
        0.0,
        controls().mass_tolerance_kg,
    );
    near(native.energy_residual_j, 0.0, controls().energy_tolerance_j);
    near(motion.energy_residual_j, 0.0, controls().energy_tolerance_j);
    assert_eq!(native.duration_s, DT);
    assert_eq!(motion.duration_s, DT);
    assert_eq!(
        native.mass_exchange.physical_time_s,
        liquid.physical_time_s()
    );
    assert_eq!(next.liquid().accepted_steps(), liquid.accepted_steps() + 1);
    assert_eq!(
        next.liquid().physical_time_s(),
        liquid.physical_time_s() + DT
    );
    assert_eq!(after.volume_m3(), before.volume_m3());
    assert_eq!(after.point().axes().get("specimen-lot"), Some(&27.0));
}

#[test]
fn evaporation_carries_native_momentum_and_heats_the_real_gas_receiver() {
    with_liquid(0.0, 0.4, |liquid, gate| {
        let source = source(-13.4e6);
        let gas = MovingGasParcel::new(gas(&source, 300.0, 0.0), [30.0, -10.0, 5.0]).unwrap();
        let binding = binding(&source, 0.0);
        let before = gas.clone();
        let bits = liquid.flow().lattice().grid().state_bits();
        let donor = liquid_velocity(liquid);
        let next = advance_liquid_with_moving_gas(
            liquid,
            &gas,
            &source,
            surface(&Phase, &binding),
            controls(),
            None,
            None,
            gate,
        )
        .unwrap();
        let native = &next.report().native;
        let motion = &next.report().momentum;
        let dm = native.phase.transferred_mass_kg;
        assert!(dm > 0.0);
        assert_eq!(motion.donor_velocity_m_s, donor);
        for d in 0..3 {
            near(
                next.gas().momentum_kg_m_s()[d],
                gas.momentum_kg_m_s()[d] + dm * donor[d],
                1e-22,
            );
        }
        let difference = (0..3)
            .map(|d| (donor[d] - gas.velocity_m_s()[d]).powi(2))
            .sum::<f64>();
        let expected_heat =
            0.5 * dm * gas.parcel().mass_kg() / (gas.parcel().mass_kg() + dm) * difference;
        near(motion.gas_mixing_heat_j, expected_heat, 1e-22);
        assert!(expected_heat > 0.0);
        assert_eq!(motion.liquid_mixing_heat_j, 0.0);
        let heat = motion.gas_heat.unwrap();
        near(heat.supplied_heat_j, expected_heat, 1e-22);
        // Both explicitly authored ideal-gas species have Cv/R=2.5.
        let capacity = next.gas().parcel().mass_kg()
            * 2.5
            * next.gas().parcel().state().state().specific_gas_constant;
        near(
            next.gas().parcel().state().state().temperature,
            native.phase.gas_temperature_k[1] + heat.stored_energy_change_j / capacity,
            1e-10,
        );
        assert!(
            next.gas().parcel().state().state().temperature > native.phase.gas_temperature_k[1]
        );
        near(
            next.gas().parcel().species_mass_kg(&carrier()).unwrap(),
            gas.parcel().species_mass_kg(&carrier()).unwrap(),
            1e-22,
        );
        check_source_and_thermal(liquid, &gas, &binding, &next);
        assert_eq!(gas, before);
        assert_eq!(liquid.flow().lattice().grid().state_bits(), bits);
    });
}

#[test]
fn condensation_uses_retained_upstream_drift_and_real_native_mixing() {
    with_liquid(0.0, 0.4, |liquid, gate| {
        let source = source(-13.4e6);
        let gas =
            MovingGasParcel::new(gas(&source, 320.0, 6000.0), [0.025, -0.015, 0.005]).unwrap();
        let binding = binding(&source, 0.0);
        let face = surface(&Phase, &binding);
        let legacy =
            advance_liquid_with_gas(liquid, gas.parcel(), &source, face, controls(), gate).unwrap();
        let next = advance_liquid_with_moving_gas(
            liquid,
            &gas,
            &source,
            face,
            controls(),
            None,
            None,
            gate,
        )
        .unwrap();
        let native = &next.report().native;
        let motion = &next.report().momentum;
        let dm = -native.phase.transferred_mass_kg;
        assert!(dm > 0.0);
        assert_eq!(
            native.phase.transferred_mass_kg,
            legacy.report().phase.transferred_mass_kg
        );
        assert_eq!(motion.donor_velocity_m_s, gas.velocity_m_s());
        for d in 0..3 {
            near(
                native.mass_exchange.carried_momentum_kg_m_s[d],
                dm * gas.velocity_m_s()[d],
                1e-24,
            );
            near(next.gas().velocity_m_s()[d], gas.velocity_m_s()[d], 1e-15);
        }
        let initial_mass = liquid
            .project_interface_mass_exchange(INTERFACE, 0.0)
            .unwrap()
            .initial_mass_kg;
        let projected = liquid
            .project_interface_mass_exchange(INTERFACE, dm)
            .unwrap();
        let u = liquid_velocity(liquid);
        let squared = (0..3)
            .map(|d| (u[d] - gas.velocity_m_s()[d]).powi(2))
            .sum::<f64>();
        near(
            motion.liquid_mixing_heat_j,
            0.5 * dm * (initial_mass / projected.final_mass_kg) * squared,
            1e-26,
        );
        assert!(motion.liquid_mixing_heat_j > 0.0);
        assert_eq!(motion.gas_mixing_heat_j, 0.0);
        assert!(motion.gas_heat.is_none());
        assert_ne!(
            next.liquid().flow().lattice().grid().state_bits(),
            legacy.liquid().flow().lattice().grid().state_bits()
        );
        check_source_and_thermal(liquid, &gas, &binding, &next);
    });
}

#[test]
fn continuation_retains_changed_gas_momentum_instead_of_resetting_an_atmosphere() {
    with_liquid(0.0, 0.4, |liquid, gate| {
        let source = source(-13.4e6);
        let gas = MovingGasParcel::new(gas(&source, 300.0, 0.0), [0.03, -0.01, 0.0]).unwrap();
        let binding = binding(&source, 0.0);
        let face = surface(&Phase, &binding);
        let first = advance_liquid_with_moving_gas(
            liquid,
            &gas,
            &source,
            face,
            controls(),
            None,
            None,
            gate,
        )
        .unwrap();
        assert_ne!(first.gas().velocity_m_s(), gas.velocity_m_s());
        let next = advance_liquid_with_moving_gas(
            first.liquid(),
            first.gas(),
            first.catalog(),
            face,
            controls(),
            None,
            None,
            gate,
        )
        .unwrap();
        assert_eq!(
            next.report().momentum.gas_momentum_before_kg_m_s,
            first.gas().momentum_kg_m_s()
        );
        let donor = liquid_velocity(first.liquid());
        let dm = next.report().native.phase.transferred_mass_kg;
        assert!(dm > 0.0);
        for d in 0..3 {
            near(
                next.gas().momentum_kg_m_s()[d],
                first.gas().momentum_kg_m_s()[d] + dm * donor[d],
                1e-24,
            );
        }
        check_source_and_thermal(first.liquid(), first.gas(), &binding, &next);
    });
}

#[test]
fn source_boost_preserves_both_receiver_losses_and_original_thermal_transfer() {
    for vapor in [0.0, 6000.0] {
        let run = |boost: [f64; 3]| {
            with_liquid_velocity(0.0, 0.4, boost, |liquid, gate| {
                let source = source(-13.4e6);
                let gas_velocity = core::array::from_fn(|d| [0.025, -0.015, 0.005][d] + boost[d]);
                let gas = MovingGasParcel::new(gas(&source, 300.0, vapor), gas_velocity).unwrap();
                let binding = binding(&source, 0.0);
                let next = advance_liquid_with_moving_gas(
                    liquid,
                    &gas,
                    &source,
                    surface(&Phase, &binding),
                    controls(),
                    None,
                    None,
                    gate,
                )
                .unwrap();
                check_source_and_thermal(liquid, &gas, &binding, &next);
                (
                    next.report().native.phase.transferred_mass_kg,
                    next.report().momentum.liquid_mixing_heat_j,
                    next.report().momentum.gas_mixing_heat_j,
                    next.report().native.phase.surface_temperature_k[1],
                )
            })
        };
        let a = run([0.0; 3]);
        let b = run([0.01, 0.015, -0.01]);
        near(a.0, b.0, 1e-21);
        near(a.1, b.1, 1e-24);
        near(a.2, b.2, 1e-24);
        near(a.3, b.3, 1e-10);
    }
}

#[test]
fn equal_mean_drift_preserves_original_phase_endpoints_without_extra_heat() {
    for vapor in [0.0, 6000.0] {
        with_liquid(0.0, 0.4, |liquid, gate| {
            let source = source(-13.4e6);
            let gas =
                MovingGasParcel::new(gas(&source, 300.0, vapor), liquid_velocity(liquid)).unwrap();
            let binding = binding(&source, 0.0);
            let face = surface(&Phase, &binding);
            let original =
                advance_liquid_with_gas(liquid, gas.parcel(), &source, face, controls(), gate)
                    .unwrap();
            let next = advance_liquid_with_moving_gas(
                liquid,
                &gas,
                &source,
                face,
                controls(),
                None,
                None,
                gate,
            )
            .unwrap();
            assert_eq!(
                next.report().native.phase.transferred_mass_kg,
                original.report().phase.transferred_mass_kg
            );
            near(next.report().momentum.gas_mixing_heat_j, 0.0, 1e-30);
            near(next.report().momentum.liquid_mixing_heat_j, 0.0, 1e-30);
            let next_gas = next.gas().parcel();
            let original_gas = original.gas();
            assert_eq!(next_gas.mass_kg(), original_gas.mass_kg());
            assert_eq!(next_gas.volume_m3(), original_gas.volume_m3());
            assert_eq!(next_gas.model(), original_gas.model());
            assert_eq!(
                next_gas.point().axes().get("specimen-lot"),
                original_gas.point().axes().get("specimen-lot")
            );
            for species in [water(), carrier()] {
                assert_eq!(
                    next_gas.species_mass_kg(&species).unwrap(),
                    original_gas.species_mass_kg(&species).unwrap()
                );
            }
            near(
                next_gas.state().state().temperature,
                original_gas.state().state().temperature,
                1e-8,
            );
            near(
                next_gas.state().state().pressure,
                original_gas.state().state().pressure,
                1e-5,
            );
            near(
                next.liquid().total_enthalpy_j().unwrap(),
                original.liquid().total_enthalpy_j().unwrap(),
                controls().energy_tolerance_j,
            );
            let next_velocity = liquid_velocity(next.liquid());
            let original_velocity = liquid_velocity(original.liquid());
            for d in 0..3 {
                near(
                    next.report().native.mass_exchange.momentum_change_kg_m_s[d],
                    original.report().mass_exchange.momentum_change_kg_m_s[d],
                    controls().native.momentum_tolerance_kg_m_s,
                );
                near(next_velocity[d], original_velocity[d], 1e-10);
            }
            check_source_and_thermal(liquid, &gas, &binding, &next);
        });
    }
}

#[test]
fn zero_phase_flux_keeps_drift_while_reciprocal_sensible_heat_advances_once() {
    with_liquid(0.0, 0.4, |liquid, gate| {
        let source = source(-13.4e6);
        let gas = MovingGasParcel::new(gas(&source, 320.0, 6000.0), [0.025, 0.01, 0.0]).unwrap();
        let binding = binding(&source, 0.0);
        let mut face = surface(&Phase, &binding);
        face.mass_transfer_coefficient_m_s = 0.0;
        let next = advance_liquid_with_moving_gas(
            liquid,
            &gas,
            &source,
            face,
            controls(),
            Some((50.0, heat_controls())),
            None,
            gate,
        )
        .unwrap();
        assert_eq!(next.report().native.phase.transferred_mass_kg, 0.0);
        assert_eq!(next.gas().momentum_kg_m_s(), gas.momentum_kg_m_s());
        assert_eq!(next.report().momentum.gas_mixing_heat_j, 0.0);
        assert_eq!(next.report().momentum.liquid_mixing_heat_j, 0.0);
        assert!(next.report().native.heat.unwrap().heat_into_liquid_j > 0.0);
        check_source_and_thermal(liquid, &gas, &binding, &next);
    });
}

#[test]
fn accumulated_sensible_and_receiver_heat_obey_the_original_single_interface_cap() {
    with_liquid_interfaces(0.0, &[(INTERFACE, 0.4, 305.0)], |liquid, gate| {
        let source = source(-13.4e6);
        let gas = MovingGasParcel::new(gas(&source, 300.0, 0.0), [300.0, -100.0, 50.0]).unwrap();
        let original = gas.clone();
        let bits = liquid.flow().lattice().grid().state_bits();
        let binding = binding(&source, 0.0);
        let face = surface(&Phase, &binding);
        let mut solve = controls();
        solve.phase.maximum_temperature_increment_k = 10.0;
        let baseline = advance_liquid_with_moving_gas(
            liquid,
            &gas,
            &source,
            face,
            solve,
            Some((50.0, heat_controls())),
            None,
            gate,
        )
        .unwrap();
        let native = &baseline.report().native;
        assert!(native.heat.is_some() && native.phase.gas_temperature_k[0] > 300.0);
        assert!(baseline.report().momentum.gas_mixing_heat_j > 0.0);
        let final_gas = baseline.gas().parcel().state().state().temperature;
        let original_increment = final_gas - gas.parcel().state().state().temperature;
        let phase_increment = (native.exchange.gas_temperature_k[1] - 300.0).abs();
        let reset_increment = (final_gas - native.phase.gas_temperature_k[0]).abs();
        let liquid_increment = (native.exchange.surface_temperature_k[1] - 305.0).abs();
        let before_mixing_limit = phase_increment.max(reset_increment).max(liquid_increment);
        assert!(original_increment > before_mixing_limit + 1e-6);
        solve.phase.maximum_temperature_increment_k =
            0.5 * (before_mixing_limit + original_increment);
        // The original thermal/phase owner fits the same complete-interval cap.
        advance_liquid_with_gas_and_heat(
            liquid,
            gas.parcel(),
            &source,
            face,
            50.0,
            solve,
            heat_controls(),
            gate,
        )
        .unwrap();
        let error = advance_liquid_with_moving_gas(
            liquid,
            &gas,
            &source,
            face,
            solve,
            Some((50.0, heat_controls())),
            None,
            gate,
        )
        .err()
        .expect("combined original-to-final temperature increment must refuse");
        assert!(matches!(
            error,
            EvaporatingLiquidGasError::Phase(GasEvaporationError::StepSizeLimit {
                reason: "phase plus receiver mixing exceeds the original temperature increment bound",
                ..
            })
        ));
        assert_eq!(gas, original);
        assert_eq!(liquid.flow().lattice().grid().state_bits(), bits);
        assert_eq!(liquid.physical_time_s(), 0.0);
    });
}

#[test]
fn optional_pressure_uses_the_actual_mixing_heated_gas_eos() {
    with_liquid_pressure(0.0, 0.4, Some(100000.0), |liquid, gate| {
        let source = source(-13.4e6);
        let gas = MovingGasParcel::new(gas(&source, 300.0, 0.0), [30.0, -10.0, 5.0]).unwrap();
        let binding = binding(&source, 0.0);
        let next = advance_liquid_with_moving_gas(
            liquid,
            &gas,
            &source,
            surface(&Phase, &binding),
            controls(),
            None,
            Some(0.1),
            gate,
        )
        .unwrap();
        let pressure = next.report().native.pressure.unwrap();
        assert_eq!(
            pressure.applied_pressure_pa,
            next.gas().parcel().state().state().pressure
        );
        assert!(next.report().momentum.gas_heat.is_some());
        check_source_and_thermal(liquid, &gas, &binding, &next);
    });
}

#[test]
fn unsupported_drift_temperature_and_late_cancellation_preserve_all_owners() {
    with_liquid(0.0, 0.4, |liquid, gate| {
        let source = source(-13.4e6);
        let binding = binding(&source, 0.0);
        let bits = liquid.flow().lattice().grid().state_bits();
        let fast = MovingGasParcel::new(gas(&source, 320.0, 6000.0), [1.0, 0.0, 0.0]).unwrap();
        let original = fast.clone();
        assert!(
            advance_liquid_with_moving_gas(
                liquid,
                &fast,
                &source,
                surface(&Phase, &binding),
                controls(),
                None,
                None,
                gate,
            )
            .is_err()
        );
        assert_eq!(fast, original);
        // The gas remains source-resolved at entry; excessive drift-derived
        // receiver heat must refuse the source interval, not be clipped away.
        let hot = MovingGasParcel::new(gas(&source, 300.0, 0.0), [1e8, 0.0, 0.0]).unwrap();
        let before = hot.clone();
        assert!(
            advance_liquid_with_moving_gas(
                liquid,
                &hot,
                &source,
                surface(&Phase, &binding),
                controls(),
                None,
                None,
                gate,
            )
            .is_err()
        );
        assert_eq!(hot, before);
        assert!(MovingGasParcel::new(gas(&source, 300.0, 0.0), [f64::NAN, 0.0, 0.0]).is_err());

        struct CancelOnCall<'a> {
            calls: AtomicUsize,
            stop: usize,
            gate: &'a CancelGate,
        }
        impl EvaporationLaw for CancelOnCall<'_> {
            fn temperature_range_k(&self) -> [f64; 2] {
                Phase.temperature_range_k()
            }
            fn properties(
                &self,
                t: f64,
            ) -> Result<fs_thermal::evaporation::EvaporationProperties, EvaporationError>
            {
                let properties = Phase.properties(t)?;
                if self.calls.fetch_add(1, Ordering::Relaxed) + 1 == self.stop {
                    self.gate.request();
                }
                Ok(properties)
            }
        }
        let count = CancelOnCall {
            calls: AtomicUsize::new(0),
            stop: usize::MAX,
            gate,
        };
        let gas = MovingGasParcel::new(gas(&source, 300.0, 0.0), [30.0, -10.0, 5.0]).unwrap();
        advance_liquid_with_moving_gas(
            liquid,
            &gas,
            &source,
            surface(&count, &binding),
            controls(),
            None,
            None,
            gate,
        )
        .unwrap();
        // The last source query is the actual post-mixing endpoint admission,
        // after private native mass/momentum and original gas heat have updated.
        let cancel = CancelOnCall {
            calls: AtomicUsize::new(0),
            stop: count.calls.load(Ordering::Relaxed),
            gate,
        };
        let original = gas.clone();
        assert!(matches!(
            advance_liquid_with_moving_gas(
                liquid,
                &gas,
                &source,
                surface(&cancel, &binding),
                controls(),
                None,
                None,
                gate,
            ),
            Err(EvaporatingLiquidGasError::Cancelled)
        ));
        assert_eq!(gas, original);
        assert_eq!(liquid.flow().lattice().grid().state_bits(), bits);
        assert_eq!(liquid.physical_time_s(), 0.0);
    });
}
