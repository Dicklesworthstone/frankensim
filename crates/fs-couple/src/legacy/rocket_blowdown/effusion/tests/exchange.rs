//! G0/G1/G3/G4: original finite gas, simultaneous molecular ports and full ledgers.
mod ambient;
use super::*;
use crate::rocket_blowdown::effusion::exchange::{
    EffusionExchangeCandidate, EffusionReservoirPort,
};
use fs_material::thermochemical_gas::energy::parcel::transfer::species::network::SpeciesNetworkError;
use fs_thermochem::{Composition, CompositionBasis, SpeciesId};

fn environment() -> EffusionRocketEnvironment {
    EffusionRocketEnvironment {
        inertial_frame: FRAME,
        gravity_world_m_s2: Vec3::ZERO,
        external_force_world_n: Vec3::ZERO,
    }
}

fn bath(before: &GasRocketState, temperature: f64, pressure: f64) -> ClosedGasParcel {
    ClosedGasParcel::new(
        before.reservoir().model().clone(),
        source::point(temperature, pressure).unwrap(),
        VOLUME,
    )
    .unwrap()
}

fn port<'a>(
    source: &'a ClosedGasParcel,
    interface: &'a EffusiveMomentumInterface,
) -> EffusionReservoirPort<'a> {
    EffusionReservoirPort {
        id: 9,
        source,
        interface,
        area_m2: AREA,
        outward_normal_world: Vec3::new(1.0, 0.0, 0.0),
    }
}

fn advance(
    before: &GasRocketState,
    environment: EffusionRocketEnvironment,
    outlets: &[EffusionOutlet<'_>],
    ports: &[EffusionReservoirPort<'_>],
) -> EffusionExchangeCandidate {
    before
        .advance_effusion_exchange(
            environment,
            outlets,
            ports,
            DT,
            endpoint(),
            controls(),
            || false,
        )
        .unwrap()
}

fn accounts(before: &GasRocketState, result: &EffusionExchangeCandidate) {
    let after = result.endpoint();
    let step = result.step();
    let dp = after
        .vehicle()
        .velocity_m_s
        .scale(mass(after))
        .sub(before.vehicle().velocity_m_s.scale(mass(before)));
    vector_near(
        dp.sub(step.net_momentum_into_vessel_world_kg_m_s)
            .sub(step.mechanics.external_impulse_n_s),
        Vec3::ZERO,
        2e-12,
    );
    near(
        internal(after) - internal(before) + kinetic(after)
            - kinetic(before)
            - step.net_energy_into_vessel_world_j
            - step.mechanics.external_work_j,
        0.0,
        3e-8,
    );
    near(
        after.reservoir().mass_kg() - before.reservoir().mass_kg(),
        step.net_mass_into_gas_kg,
        2e-15,
    );
    near(
        internal(after) - internal(before),
        step.net_relative_energy_into_gas_j,
        2e-8,
    );
    let (net_mass, relative_energy, net_momentum, energy) =
        step.ports
            .iter()
            .fold((0.0, 0.0, Vec3::ZERO, 0.0), |(m, r, p, e), port| {
                (
                    m + port.net_mass_into_gas_kg,
                    r + port.net_relative_energy_into_gas_j,
                    p.add(port.net_momentum_into_vessel_world_kg_m_s),
                    e + port.net_energy_into_vessel_world_j,
                )
            });
    near(net_mass, step.net_mass_into_gas_kg, 2e-15);
    near(relative_energy, step.net_relative_energy_into_gas_j, 2e-12);
    vector_near(
        net_momentum,
        step.net_momentum_into_vessel_world_kg_m_s,
        2e-12,
    );
    near(energy, step.net_energy_into_vessel_world_j, 2e-12);
    vector_near(step.momentum_residual_kg_m_s, Vec3::ZERO, 1e-12);
    near(
        step.energy_residual_j,
        0.0,
        controls().transport.exchange.energy_tolerance_j,
    );
    assert_eq!(
        after.vehicle().propellant_mass_kg,
        after.reservoir().mass_kg()
    );
    assert_eq!(after.vehicle().dry_mass_kg, before.vehicle().dry_mass_kg);
    assert_eq!(
        after.reservoir().volume_m3(),
        before.reservoir().volume_m3()
    );
}

#[test]
fn hot_comoving_source_fills_and_heats_the_actual_gas_with_both_beam_moments() {
    let before = initial(PURE, Vec3::ZERO);
    let saved = before.clone();
    let source = ClosedGasParcel::new(
        before.reservoir().model().clone(),
        source::point(2.0 * T0, 3.0 * P0).unwrap(),
        1e-12,
    )
    .unwrap();
    let saved_source = source.clone();
    let law = interface(before.reservoir().model());
    let result = advance(&before, environment(), &[], &[port(&source, &law)]);
    let no = escaped_moles(0.028, 1.0);
    let ni = no * 3.0 / 2.0_f64.sqrt();
    let n0 = P0 * VOLUME / (R * T0);
    let dm = (ni - no) * 0.028;
    let impulse = AREA * P0 * DT / 2.0;
    let step = result.step();
    near(step.incoming.mass_kg, ni * 0.028, 2e-18);
    assert!(step.incoming.mass_kg > source.mass_kg());
    near(step.outgoing.mass_kg, no * 0.028, 2e-18);
    vector_near(
        step.incoming.relative_momentum_kg_m_s,
        Vec3::new(-3.0 * impulse, 0.0, 0.0),
        2e-16,
    );
    vector_near(
        step.outgoing.relative_momentum_kg_m_s,
        Vec3::new(impulse, 0.0, 0.0),
        2e-16,
    );
    near(
        step.incoming.relative_energy_j,
        2.0 * R * ni * 2.0 * T0,
        2e-12,
    );
    near(step.outgoing.relative_energy_j, 2.0 * R * no * T0, 2e-12);
    near(
        result.endpoint().reservoir().mass_kg(),
        before.reservoir().mass_kg() + dm,
        2e-15,
    );
    near(
        result.endpoint().reservoir().state().state().temperature,
        (n0 * T0 + (4.0 / 3.0) * (ni * 2.0 * T0 - no * T0)) / (n0 + ni - no),
        2e-7,
    );
    let dv = -4.0 * impulse * (dm / mass(&before)).ln_1p() / dm;
    vector_near(
        result.endpoint().vehicle().velocity_m_s,
        Vec3::new(dv, 0.0, 0.0),
        2e-11,
    );
    accounts(&before, &result);
    assert_eq!(before, saved);
    assert_eq!(source, saved_source); // This is a prescribed source, not a depleted finite receiver.
}

#[test]
fn hot_counterflow_and_vacuum_share_one_reconstruction_at_zero_net_mass() {
    let before = initial(PURE, Vec3::ZERO);
    let law = interface(before.reservoir().model());
    // Source mass flux is twice the vessel beam: 4p/sqrt(4T).
    let source = bath(&before, 4.0 * T0, 4.0 * P0);
    let vacuum = outlet(&law, Vec3::new(1.0, 0.0, 0.0));
    let result = advance(&before, environment(), &[vacuum], &[port(&source, &law)]);
    let no = escaped_moles(0.028, 1.0);
    near(result.step().incoming.mass_kg, 2.0 * no * 0.028, 2e-18);
    near(
        result.step().outgoing.mass_kg,
        result.step().incoming.mass_kg,
        2e-18,
    );
    near(
        result.endpoint().reservoir().mass_kg(),
        before.reservoir().mass_kg(),
        2e-15,
    );
    let n0 = P0 * VOLUME / (R * T0);
    near(
        result.endpoint().reservoir().state().state().temperature,
        T0 + (4.0 / 3.0) * (2.0 * no * (4.0 * T0 - T0)) / n0,
        2e-7,
    );
    // All incoming and outgoing normal moments add; zero net mass is not zero stress.
    let dv = -3.0 * AREA * P0 * DT / mass(&before);
    vector_near(
        result.endpoint().vehicle().velocity_m_s,
        Vec3::new(dv, 0.0, 0.0),
        2e-11,
    );
    vector_near(
        result.endpoint().vehicle().position_m,
        Vec3::new(0.5 * dv * DT, 0.0, 0.0),
        2e-12,
    );
    assert_eq!(
        result
            .step()
            .ports
            .iter()
            .map(|row| row.id)
            .collect::<Vec<_>>(),
        vec![4, 9]
    );
    let dm = no * 0.028;
    let energy = 2.0 * R * no * T0;
    let impulse = AREA * P0 * DT / 2.0;
    let mean = result.step().mechanics.mean_velocity_world_m_s;
    let mean_squared = result.step().mechanics.mean_speed_squared_m2_s2;
    for (receipt, (m, e, j)) in result
        .step()
        .ports
        .iter()
        .zip([(-dm, -energy, -impulse), (dm, 7.0 * energy, -5.0 * impulse)])
    {
        near(receipt.net_mass_into_gas_kg, m, 2e-18);
        near(receipt.net_relative_energy_into_gas_j, e, 2e-12);
        vector_near(
            receipt.net_momentum_into_vessel_world_kg_m_s,
            mean.scale(m).add(Vec3::new(j, 0.0, 0.0)),
            2e-12,
        );
        near(
            receipt.net_energy_into_vessel_world_j,
            e + j * mean.x + 0.5 * m * mean_squared,
            2e-12,
        );
    }
    accounts(&before, &result);
}

#[test]
fn selected_new_species_enters_without_outflow_using_the_retained_full_catalog() {
    let catalog = source::load(MIXTURE, "effusive recoil analytic fixture").unwrap();
    let pure = |fractions| {
        catalog
            .with_composition(Composition::new(CompositionBasis::MoleFraction, fractions).unwrap())
            .unwrap()
    };
    let gas =
        ClosedGasParcel::new(pure(vec![1.0, 0.0]), source::point(T0, P0).unwrap(), VOLUME).unwrap();
    let before = GasRocketState::new(gas, 0.01, Vec3::ZERO, Vec3::ZERO, FRAME).unwrap();
    let source =
        ClosedGasParcel::new(pure(vec![0.0, 1.0]), source::point(T0, P0).unwrap(), VOLUME).unwrap();
    let original = interface(&catalog);
    let b = SpeciesId::new("B").unwrap();
    let transmissions = catalog
        .thermodynamics()
        .components()
        .iter()
        .map(|species| {
            EffusiveSpeciesTransmission::try_new(
                species.species().clone(),
                if species.species() == &b { 1.0 } else { 0.0 },
                "synthetic B-selective reciprocal aperture",
            )
            .unwrap()
        })
        .collect();
    let law = EffusiveMomentumInterface::try_new(
        EffusiveInterface::try_new(
            catalog.clone(),
            transmissions,
            original.interface().support(),
            original.interface().regime().clone(),
        )
        .unwrap(),
        "synthetic uncollimated angle-independent B transmission",
    )
    .unwrap();
    let result = advance(&before, environment(), &[], &[port(&source, &law)]);
    let ni = escaped_moles(0.004, 1.0);
    let n0 = P0 * VOLUME / (R * T0);
    assert_eq!(result.step().outgoing.mass_kg, 0.0);
    near(result.step().incoming.mass_kg, ni * 0.004, 2e-18);
    assert_eq!(before.reservoir().species_mass_kg(&b).unwrap(), 0.0);
    near(
        result.endpoint().reservoir().species_mass_kg(&b).unwrap(),
        ni * 0.004,
        2e-15,
    );
    near(
        result
            .endpoint()
            .reservoir()
            .species_mass_kg(&SpeciesId::new("A").unwrap())
            .unwrap(),
        before.reservoir().mass_kg(),
        2e-15,
    );
    near(
        result.endpoint().reservoir().state().state().temperature,
        T0 * (n0 + (4.0 / 3.0) * ni) / (n0 + ni),
        2e-7,
    );
    accounts(&before, &result);
}

#[test]
fn incoming_and_outgoing_inertial_accounts_transform_with_boost_and_external_work() {
    let velocity = Vec3::new(2.0, -1.0, 3.0);
    let boost = Vec3::new(19.0, -7.0, 11.0);
    let before = initial(MIXTURE, velocity);
    let boosted = initial(MIXTURE, velocity.add(boost));
    let law = interface(before.reservoir().model());
    let source = bath(&before, 600.0, 2.4 * P0);
    let source_port = EffusionReservoirPort {
        outward_normal_world: Vec3::new(0.6, 0.8, 0.0),
        ..port(&source, &law)
    };
    let vacuum = outlet(&law, Vec3::new(0.0, 0.0, -1.0));
    let loads = EffusionRocketEnvironment {
        gravity_world_m_s2: Vec3::new(0.3, -0.5, 0.2),
        external_force_world_n: Vec3::new(1e-3, 2e-4, -3e-4),
        ..environment()
    };
    let a = advance(&before, loads, &[vacuum], &[source_port]);
    let b = advance(&boosted, loads, &[vacuum], &[source_port]);
    vector_near(
        b.endpoint().vehicle().velocity_m_s,
        a.endpoint().vehicle().velocity_m_s.add(boost),
        3e-11,
    );
    vector_near(
        b.endpoint().vehicle().position_m,
        a.endpoint().vehicle().position_m.add(boost.scale(DT)),
        3e-12,
    );
    for (original, shifted) in [
        (&a.step().incoming, &b.step().incoming),
        (&a.step().outgoing, &b.step().outgoing),
    ] {
        near(shifted.mass_kg, original.mass_kg, 2e-18);
        vector_near(
            shifted.momentum_world_kg_m_s,
            original
                .momentum_world_kg_m_s
                .add(boost.scale(original.mass_kg)),
            2e-12,
        );
        near(
            shifted.energy_world_j,
            original.energy_world_j
                + boost.dot(original.momentum_world_kg_m_s)
                + 0.5 * original.mass_kg * boost.norm_squared(),
            3e-8,
        );
    }
    vector_near(
        b.step().mechanics.external_impulse_n_s,
        a.step().mechanics.external_impulse_n_s,
        2e-12,
    );
    near(
        b.step().mechanics.external_work_j,
        a.step().mechanics.external_work_j + boost.dot(a.step().mechanics.external_impulse_n_s),
        3e-8,
    );
    accounts(&before, &a);
    accounts(&boosted, &b);
}

#[test]
fn incoming_supply_cannot_fund_shared_gross_original_species_withdrawal() {
    let before = initial(PURE, Vec3::ZERO);
    let saved = before.clone();
    let source = bath(&before, T0, 2.0 * P0);
    let law = interface(before.reservoir().model());
    let vacuum = outlet(&law, Vec3::new(1.0, 0.0, 0.0));
    let source_port = port(&source, &law);
    let mut bounded = controls();
    bounded.transport.exchange.maximum_species_fraction =
        1.5 * escaped_moles(0.028, 1.0) * 0.028 / before.reservoir().mass_kg();
    let attempt = |outlets: &[EffusionOutlet<'_>], ports: &[EffusionReservoirPort<'_>]| {
        before.advance_effusion_exchange(
            environment(),
            outlets,
            ports,
            DT,
            endpoint(),
            bounded,
            || false,
        )
    };
    assert!(attempt(&[vacuum], &[]).is_ok());
    assert!(attempt(&[], &[source_port]).is_ok());
    assert!(matches!(
        attempt(&[vacuum], &[source_port]),
        Err(EffusionRocketError::Gas(GasTransportError::Exchange(
            SpeciesNetworkError::Withdrawal { id: 0, .. }
        )))
    ));
    assert_eq!(before, saved);
}

#[test]
fn late_source_and_cancel_refusals_are_atomic_and_closed_foreign_ports_coast() {
    let before = initial(PURE, Vec3::new(4.0, -3.0, 2.0));
    let saved = before.clone();
    let source = bath(&before, T0, 2.0 * P0);
    let law = interface(before.reservoir().model());
    let path = port(&source, &law);
    let mut polls = 0;
    let accepted = before
        .advance_effusion_exchange(
            environment(),
            &[],
            &[path],
            DT,
            endpoint(),
            controls(),
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    let final_poll = polls;
    polls = 0;
    assert!(matches!(
        before.advance_effusion_exchange(
            environment(),
            &[],
            &[path],
            DT,
            endpoint(),
            controls(),
            || {
                polls += 1;
                polls == final_poll
            }
        ),
        Err(EffusionRocketError::Cancelled)
    ));
    // Initial gas/source are inside this support; influx heats the final gas beyond it.
    let narrow = EffusiveMomentumInterface::try_new(
        EffusiveInterface::try_new(
            before.reservoir().model().clone(),
            law.interface().transmissions().to_vec(),
            EffusiveSupport {
                temperature_k: (200.0, T0 + 1e-6),
                ..law.interface().support()
            },
            law.interface().regime().clone(),
        )
        .unwrap(),
        "synthetic narrow original support",
    )
    .unwrap();
    assert!(
        before
            .advance_effusion_exchange(
                environment(),
                &[],
                &[EffusionReservoirPort {
                    interface: &narrow,
                    ..path
                }],
                DT,
                endpoint(),
                controls(),
                || false
            )
            .is_err()
    );
    let foreign = initial("A 1 0.028 3.5 0 1.7e-5 0.025", Vec3::ZERO);
    let foreign_path = EffusionReservoirPort {
        source: foreign.reservoir(),
        ..path
    };
    assert!(
        before
            .advance_effusion_exchange(
                environment(),
                &[],
                &[foreign_path],
                DT,
                endpoint(),
                controls(),
                || false
            )
            .is_err()
    );
    assert!(
        before
            .advance_effusion_exchange(
                EffusionRocketEnvironment {
                    inertial_frame: FRAME + 1,
                    ..environment()
                },
                &[],
                &[path],
                DT,
                endpoint(),
                controls(),
                || false
            )
            .is_err()
    );
    assert!(
        before
            .advance_effusion_exchange(
                environment(),
                &[],
                &[EffusionReservoirPort {
                    outward_normal_world: Vec3::ZERO,
                    ..path
                }],
                DT,
                endpoint(),
                controls(),
                || false
            )
            .is_err()
    );
    let coast = advance(
        &before,
        environment(),
        &[],
        &[EffusionReservoirPort {
            area_m2: 0.0,
            ..foreign_path
        }],
    );
    assert_eq!(coast.endpoint().reservoir(), before.reservoir());
    assert_eq!(
        coast.endpoint().vehicle().velocity_m_s,
        before.vehicle().velocity_m_s
    );
    vector_near(
        coast.endpoint().vehicle().position_m,
        before
            .vehicle()
            .position_m
            .add(before.vehicle().velocity_m_s.scale(DT)),
        2e-14,
    );
    assert_eq!(coast.step().incoming.mass_kg, 0.0);
    assert_eq!(coast.step().outgoing.mass_kg, 0.0);
    accounts(&before, &coast);
    assert_eq!(advance(&before, environment(), &[], &[path]), accepted);
    assert_eq!(before, saved);
}
