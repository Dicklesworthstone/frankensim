//! G0/G1/G3/G4: finite material heat and bidirectional molecular mechanics.
use super::*;
use crate::rocket_blowdown::effusion::{
    exchange::{EffusionAmbientPort, EffusionReservoirPort},
    thermal::exchange::EffusionWallExchangeCandidate,
};
use fs_material::thermochemical_gas::energy::parcel::transfer::species::{
    effusion::drift::EffusiveDriftInterface,
    network::{SpeciesNetworkError, transport::GasTransportError},
};

fn bath(before: &GasRocketState, t: f64, p: f64) -> ClosedGasParcel {
    ClosedGasParcel::new(
        before.reservoir().model().clone(), source::point(t, p).unwrap(), VOLUME,
    ).unwrap()
}
fn supply<'a>(gas: &'a ClosedGasParcel, law: &'a EffusiveMomentumInterface) -> EffusionReservoirPort<'a> {
    EffusionReservoirPort {
        id: 9, source: gas, interface: law, area_m2: AREA,
        outward_normal_world: Vec3::new(0.0, 1.0, 0.0),
    }
}
fn drifting(law: EffusiveMomentumInterface) -> EffusiveDriftInterface {
    EffusiveDriftInterface::try_new(
        law, 1000.0, "synthetic dilute nondeflecting molecular exposure over the declared drift ball",
    ).unwrap()
}
fn ambient<'a>(gas: &'a ClosedGasParcel, law: &'a EffusiveDriftInterface, velocity: Vec3) -> EffusionAmbientPort<'a> {
    EffusionAmbientPort {
        id: 3, source: gas, interface: law, velocity_world_m_s: velocity,
        inertial_frame: FRAME, area_m2: AREA, outward_normal_world: Vec3::new(1.0, 0.0, 0.0),
    }
}
fn accounts<S>(before: &GasRocketState, input: UniformEnthalpyStepInput<'_>, result: &EffusionWallExchangeCandidate<S>) {
    let after = result.endpoint();
    let report = result.step();
    let wall_energy = input.mass_kg * (result.wall().specific_enthalpy_j_kg() - input.initial.specific_enthalpy_j_kg());
    near(after.reservoir().mass_kg() - before.reservoir().mass_kg(), report.exchange.net_mass_into_gas_kg, 2e-15);
    vector_near(after.vehicle().velocity_m_s.scale(mass(after))
        .sub(before.vehicle().velocity_m_s.scale(mass(before)))
        .sub(report.exchange.net_momentum_into_vessel_world_kg_m_s)
        .sub(report.exchange.mechanics.external_impulse_n_s), Vec3::ZERO, 3e-12);
    near(internal(after) - internal(before) + wall_energy + kinetic(after) - kinetic(before)
        - report.exchange.net_energy_into_vessel_world_j - report.exchange.mechanics.external_work_j
        - report.external_thermal_heat_j, 0.0, 5e-7);
    assert!(report.absolute_energy_residual_sum_j + report.energy_uncertainty_bound_j <= solve().combined_energy_tolerance_j);
    assert_eq!(after.vehicle().dry_mass_kg, before.vehicle().dry_mass_kg);
    assert_eq!(after.reservoir().volume_m3(), before.reservoir().volume_m3());
    assert_eq!(report.duration_s, input.duration_s);
}

#[test]
fn finite_wall_changes_comoving_supply_and_discharge_through_one_energy_and_momentum_account() {
    scope(|cx, _| {
        let before = initial(PURE, Vec3::ZERO);
        let saved = before.clone();
        let source = bath(&before, 300.0, 0.9 * P0);
        let saved_source = source.clone();
        let law = interface(before.reservoir().model());
        let chart = curve([1000.0; 2]);
        let input = wall(&chart, 600.0);
        let result = before.advance_effusion_exchange_with_wall(
            cx, free(), &[outlet(&law, Vec3::new(1.0, 0.0, 0.0))], &[supply(&source, &law)], &[],
            input, &body(100.0), 0.0, solve(),
        ).unwrap();

        // Independent constant-Cv implicit heat balance followed by gross
        // molecular beams. No native sampler supplies these expected values.
        let n0 = P0 * VOLUME / (R * T0);
        let gas_capacity = 1.5 * R * n0;
        let wall_capacity = WALL_MASS * WALL_CP;
        let conductance_time = DT * 100.0 * WALL_AREA;
        let q = conductance_time * (600.0 - T0) /
            (1.0 + conductance_time * (1.0 / gas_capacity + 1.0 / wall_capacity));
        let heated = T0 + q / gas_capacity;
        let nout = escaped_moles(0.028, 1.0) * (heated / T0).sqrt();
        let nin = DT * AREA * 0.9 * P0 / (2.0 * core::f64::consts::PI * 0.028 * R * 300.0).sqrt();
        let n = n0 + nin - 2.0 * nout;
        let crossing_energy = 2.0 * R * (nin * 300.0 - 2.0 * nout * heated);
        let final_t = (1.5 * R * n0 * T0 + q + crossing_energy) / (1.5 * R * n);
        near(result.step().thermal.report.solid_to_gas_heat_j, q, 3e-8);
        near(result.wall().temperature_k(), 600.0 - q / wall_capacity, 2e-8);
        near(result.endpoint().reservoir().mass_kg(), n * 0.028, 2e-15);
        near(result.endpoint().reservoir().state().state().temperature, final_t, 5e-7);
        near(result.step().exchange.incoming.mass_kg, nin * 0.028, 2e-15);
        near(result.step().exchange.outgoing.mass_kg, 2.0 * nout * 0.028, 2e-15);
        near(result.step().exchange.net_relative_energy_into_gas_j, crossing_energy, 3e-8);
        let impulse = DT * AREA * P0 * heated / (2.0 * T0);
        vector_near(result.step().exchange.outgoing.relative_momentum_kg_m_s,
            Vec3::new(impulse, impulse, 0.0), 3e-12);
        vector_near(result.step().exchange.incoming.relative_momentum_kg_m_s,
            Vec3::new(0.0, -DT * AREA * 0.9 * P0 / 2.0, 0.0), 3e-12);
        accounts(&before, input, &result);
        assert_eq!(before, saved);
        assert_eq!(source, saved_source);
        assert_eq!(input.initial.temperature_k(), 600.0);
    });
}

#[test]
fn mixed_ambient_ports_preserve_tangential_particle_energy_and_common_inertial_boosts() {
    scope(|cx, _| {
        let before = initial(MIXTURE, Vec3::new(2.0, -1.0, 3.0));
        let boost = Vec3::new(31.0, -7.0, 13.0);
        let boosted = initial(MIXTURE, before.vehicle().velocity_m_s.add(boost));
        let held = bath(&before, 300.0, 0.8 * P0);
        let ambient_gas = bath(&before, 700.0, 1.3 * P0);
        let saved_ambient = ambient_gas.clone();
        let law = drifting(interface(before.reservoir().model()));
        let chart = curve([1000.0; 2]);
        let mut input = wall(&chart, 600.0);
        input.internal_heat_j = 0.2;
        let slip = Vec3::new(-20.0, 35.0, -12.0);
        let run = |state: &GasRocketState| {
            state.advance_effusion_exchange_with_wall(
                cx, free(), &[outlet(law.momentum_interface(), Vec3::new(0.0, 0.0, 1.0))],
                &[supply(&held, law.momentum_interface())],
                &[ambient(&ambient_gas, &law, state.vehicle().velocity_m_s.add(slip))],
                input, &body(100.0), 0.3, solve(),
            ).unwrap()
        };
        let a = run(&before);
        let b = run(&boosted);
        assert_eq!(a.step().thermal, b.step().thermal);
        assert_eq!(a.endpoint().reservoir(), b.endpoint().reservoir());
        assert_eq!(a.wall(), b.wall());
        let port = a.step().exchange.ports.iter().find(|p| p.id == 3).unwrap();
        assert!(port.incoming.mass_kg > 0.0);
        near(port.incoming.relative_momentum_kg_m_s.y, port.incoming.mass_kg * slip.y, 2e-14);
        near(a.step().external_thermal_heat_j, 0.5, 1e-14);
        vector_near(b.endpoint().vehicle().velocity_m_s.sub(a.endpoint().vehicle().velocity_m_s), boost, 1e-12);
        vector_near(b.endpoint().vehicle().position_m.sub(a.endpoint().vehicle().position_m), boost.scale(DT), 1e-12);
        let x = &a.step().exchange;
        let y = &b.step().exchange;
        vector_near(y.net_momentum_into_vessel_world_kg_m_s,
            x.net_momentum_into_vessel_world_kg_m_s.add(boost.scale(x.net_mass_into_gas_kg)), 3e-12);
        near(y.net_energy_into_vessel_world_j, x.net_energy_into_vessel_world_j
            + boost.dot(x.net_momentum_into_vessel_world_kg_m_s)
            + 0.5 * boost.norm_squared() * x.net_mass_into_gas_kg, 2e-9);
        accounts(&before, input, &a);
        accounts(&boosted, input, &b);
        assert_eq!(ambient_gas, saved_ambient);
    });
}

#[test]
fn insulated_wall_is_exact_original_bidirectional_owner_and_initial_inventory_spends_budget() {
    scope(|cx, _| {
        let before = initial(MIXTURE, Vec3::new(3.0, -2.0, 1.0));
        let source = bath(&before, 600.0, 1.2 * P0);
        let law = drifting(interface(before.reservoir().model()));
        let outlets = [outlet(law.momentum_interface(), Vec3::new(0.0, 0.0, 1.0))];
        let supplies = [supply(&source, law.momentum_interface())];
        let ambients = [ambient(&source, &law, Vec3::new(-30.0, 10.0, 0.0))];
        let chart = curve([1000.0; 2]);
        let input = wall(&chart, 600.0);
        let native = before.advance_effusion_exchange_with_ambient(
            free(), &outlets, &supplies, &ambients, DT, endpoint(), controls(), || false,
        ).unwrap();
        let combined = before.advance_effusion_exchange_with_wall(
            cx, free(), &outlets, &supplies, &ambients, input, &body(0.0), 0.0, solve(),
        ).unwrap();
        assert_eq!(combined.endpoint(), native.endpoint());
        assert_eq!(&combined.step().exchange, native.step());
        assert_eq!(*combined.wall(), input.initial);
        accounts(&before, input, &combined);

        let mut rounded = before.clone();
        rounded.vehicle.propellant_mass_kg = rounded.vehicle.propellant_mass_kg.next_up();
        let saved = rounded.clone();
        let mut exact = solve();
        exact.combined_energy_tolerance_j = 0.0;
        assert!(matches!(rounded.advance_effusion_exchange_with_wall(
            cx, free(), &[], &[], &[], input, &body(0.0), 0.0, exact,
        ), Err(EffusionWallError::EnergyBalance { uncertainty_bound_j, .. }) if uncertainty_bound_j > 0.0));
        assert_eq!(rounded, saved);
    });
}

#[test]
fn original_source_gates_and_late_withdrawal_callback_and_cancellation_refusals_are_atomic() {
    scope(|cx, gate| {
        let before = initial(PURE, Vec3::ZERO);
        let saved = before.clone();
        let source = bath(&before, 600.0, 2.0 * P0);
        let saved_source = source.clone();
        let law = interface(before.reservoir().model());
        let chart = curve([1000.0; 2]);
        let input = wall(&chart, 600.0);
        let boundary = body(100.0);
        let callback = |gas: &fs_material::thermochemical_gas::ResolvedThermochemicalGasState| {
            Ok::<_, &'static str>(ConvectionCoefficient {
                coefficient_w_per_m2_k: gas.state().temperature / 4.0,
                source: gas.state().temperature,
            })
        };
        let candidate = before.advance_effusion_exchange_with_wall_convection(
            cx, free(), &[], &[supply(&source, &law)], &[], input, &boundary, 0.0, solve(), callback,
        ).unwrap();
        assert_eq!(candidate.step().thermal.report.convection.source,
            candidate.step().thermal.report.next_gas.state().state().temperature);
        accounts(&before, input, &candidate);

        let narrow = EffusiveMomentumInterface::try_new(EffusiveInterface::try_new(
            law.interface().gas_catalog().clone(), law.interface().transmissions().to_vec(),
            EffusiveSupport { temperature_k: (410.0, 1800.0), ..law.interface().support() },
            law.interface().regime().clone(),
        ).unwrap(), law.provenance()).unwrap();
        let mut calls = 0;
        assert!(matches!(before.advance_effusion_exchange_with_wall_convection(
            cx, free(), &[], &[supply(&source, &narrow)], &[], input, &boundary, 0.0, solve(),
            |gas| { calls += 1; callback(gas) },
        ), Err(EffusionWallError::Effusion(EffusionRocketError::Source { .. }))));
        assert_eq!(calls, 0, "heating cannot hide the unqualified original gas");
        let closed = before.advance_effusion_exchange_with_wall(
            cx, free(), &[], &[EffusionReservoirPort { area_m2: 0.0, ..supply(&source, &narrow) }], &[],
            input, &boundary, 0.0, solve(),
        ).unwrap();
        assert_eq!(closed.step().exchange.incoming.mass_kg, 0.0);
        assert_eq!(closed.step().exchange.outgoing.mass_kg, 0.0);
        assert!(closed.step().thermal.report.solid_to_gas_heat_j > 0.0);

        let mut limited = solve();
        limited.recoil.transport.exchange.maximum_species_fraction = 1e-12;
        let failed = before.advance_effusion_exchange_with_wall_convection(
            cx, free(), &[], &[supply(&source, &law)], &[], input, &boundary, 0.0, limited,
            |gas| { calls += 1; callback(gas) },
        );
        assert!(calls > 1);
        assert!(matches!(failed, Err(EffusionWallError::Effusion(EffusionRocketError::Gas(
            GasTransportError::Exchange(SpeciesNetworkError::Withdrawal { .. })
        )))));
        assert!(matches!(before.advance_effusion_exchange_with_wall_convection(
            cx, free(), &[], &[supply(&source, &law)], &[], input, &boundary, 0.0, solve(),
            |_| Err::<ConvectionCoefficient<()>, _>("unsupported heat source"),
        ), Err(EffusionWallError::Thermal(ConvectionError::Convection("unsupported heat source")))));
        assert!(before.advance_effusion_exchange_with_wall_convection(
            cx, free(), &[], &[supply(&source, &law)], &[], input, &boundary, 0.0, solve(),
            |gas| { gate.request(); callback(gas) },
        ).is_err());
        assert_eq!(before, saved);
        assert_eq!(source, saved_source);
        assert_eq!(input.initial.temperature_k(), 600.0);
    });
}
