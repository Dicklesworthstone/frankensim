//! G1/G3/G4: shared-controller parity, new species, accepted receipts and refusal.
use super::*;
use crate::rocket_blowdown::effusion::{
    exchange::{EffusionExchangeStep, EffusionReservoirPort},
    trajectory::exchange::EffusionExchangeTrajectoryCandidate,
};
use fs_thermochem::{
    Composition, CompositionBasis, SpeciesId, mixture::FrozenIdealGasMixtureModelV1,
};

fn port<'a>(
    source: &'a ClosedGasParcel,
    law: &'a EffusiveMomentumInterface,
) -> EffusionReservoirPort<'a> {
    EffusionReservoirPort {
        id: 9,
        source,
        interface: law,
        area_m2: AREA,
        outward_normal_world: Vec3::new(1.0, 0.0, 0.0),
    }
}
fn exchange(
    before: &GasRocketState,
    outlets: &[EffusionOutlet<'_>],
    ports: &[EffusionReservoirPort<'_>],
    duration: f64,
    accuracy: EffusionAdaptiveSolve,
) -> EffusionExchangeTrajectoryCandidate {
    before
        .advance_effusion_exchange_adaptive(
            quiet(),
            outlets,
            ports,
            duration,
            endpoint(),
            controls(),
            accuracy,
            || false,
        )
        .unwrap()
}

#[test]
fn no_reservoirs_preserve_the_original_vacuum_adaptive_trajectory() {
    let before = initial(MIXTURE, Vec3::new(2.0, -1.0, 0.5));
    let law = interface(before.reservoir().model());
    let outlets = [outlet(&law, Vec3::new(0.6, 0.8, 0.0))];
    let original = evolve(&before, &outlets, 1.0, adaptive(0.1));
    let result = exchange(&before, &outlets, &[], 1.0, adaptive(0.1));
    let (a, b) = (original.endpoint(), result.endpoint());
    near(
        a.reservoir().state().state().temperature,
        b.reservoir().state().state().temperature,
        1e-8,
    );
    near(
        a.reservoir().state().state().pressure,
        b.reservoir().state().state().pressure,
        1e-5,
    );
    near(a.reservoir().mass_kg(), b.reservoir().mass_kg(), 1e-13);
    vector_near(a.vehicle().position_m, b.vehicle().position_m, 1e-10);
    vector_near(a.vehicle().velocity_m_s, b.vehicle().velocity_m_s, 1e-9);
    near(
        original.step().total_escaped_mass_kg,
        result.step().total_outgoing.mass_kg,
        1e-13,
    );
    vector_near(
        original.step().total_escaped_momentum_world_kg_m_s,
        result.step().total_outgoing.momentum_world_kg_m_s,
        1e-10,
    );
    near(
        original.step().total_escaped_energy_world_j,
        result.step().total_outgoing.energy_world_j,
        1e-7,
    );
    assert_eq!(result.step().total_incoming.mass_kg, 0.0);
    assert_eq!(result.step().duration_s, 1.0);
    assert!(result.step().maximum_accepted_error <= 1.0);
}

#[test]
fn species_absent_from_original_catalog_controls_refinement_from_zero_mass() {
    let catalog = source::load(MIXTURE, "effusive recoil analytic fixture").unwrap();
    let pure = |fractions| {
        catalog
            .with_composition(Composition::new(CompositionBasis::MoleFraction, fractions).unwrap())
            .unwrap()
    };
    let gas =
        ClosedGasParcel::new(pure(vec![1.0, 0.0]), source::point(T0, P0).unwrap(), VOLUME).unwrap();
    let before = GasRocketState::new(gas, 0.01, Vec3::ZERO, Vec3::ZERO, FRAME).unwrap();
    let saved = before.clone();
    let source =
        ClosedGasParcel::new(pure(vec![0.0, 1.0]), source::point(T0, P0).unwrap(), VOLUME).unwrap();
    let original = interface(&catalog);
    let b = SpeciesId::new("B").unwrap();
    assert_eq!(
        before
            .reservoir()
            .model()
            .thermodynamics()
            .components()
            .len(),
        1
    );
    assert!(
        before
            .reservoir()
            .model()
            .thermodynamics()
            .components()
            .iter()
            .all(|species| species.species() != &b)
    );
    let transmissions = catalog
        .thermodynamics()
        .components()
        .iter()
        .map(|species| {
            EffusiveSpeciesTransmission::try_new(
                species.species().clone(),
                if species.species() == &b { 1.0 } else { 0.0 },
                "synthetic B-only reciprocal aperture",
            )
            .unwrap()
        })
        .collect();
    let law = EffusiveMomentumInterface::try_new(
        EffusiveInterface::try_new(
            catalog,
            transmissions,
            original.interface().support(),
            original.interface().regime().clone(),
        )
        .unwrap(),
        "synthetic uncollimated B transmission",
    )
    .unwrap();
    let ports = [port(&source, &law)];
    // Only the entering B inventory can cause temporal rejection here. A is
    // retained exactly, while all other dimensional local tolerances are loose.
    let accuracy = |mass_tolerance| EffusionAdaptiveSolve {
        relative_tolerance: 0.0,
        temperature_tolerance_k: 1e6,
        pressure_tolerance_pa: 1e9,
        position_tolerance_m: 1e6,
        velocity_tolerance_m_s: 1e6,
        species_mass_tolerance_kg: mass_tolerance,
        ..adaptive(1.0)
    };
    let loose = exchange(&before, &[], &ports, 1.0, accuracy(2e-8));
    let tight = exchange(&before, &[], &ports, 1.0, accuracy(2e-10));
    let mut reference = before.clone();
    for _ in 0..512 {
        reference = reference
            .advance_effusion_exchange(
                quiet(),
                &[],
                &ports,
                1.0 / 512.0,
                endpoint(),
                controls(),
                || false,
            )
            .unwrap()
            .into_parts()
            .0;
    }
    let mass_b = |state: &GasRocketState| state.reservoir().species_mass_kg(&b).unwrap();
    assert_eq!(mass_b(&before), 0.0);
    assert!(mass_b(tight.endpoint()) > 0.0);
    assert!(tight.step().accepted_intervals > loose.step().accepted_intervals);
    assert!(tight.step().rejected_intervals > 0);
    assert!(
        (mass_b(tight.endpoint()) - mass_b(&reference)).abs()
            < 0.25 * (mass_b(loose.endpoint()) - mass_b(&reference)).abs()
    );
    let exhausted = EffusionAdaptiveSolve {
        maximum_attempts: 1,
        ..accuracy(2e-10)
    };
    assert!(matches!(
        before.advance_effusion_exchange_adaptive(
            quiet(),
            &[],
            &ports,
            1.0,
            endpoint(),
            controls(),
            exhausted,
            || false
        ),
        Err(EffusionTrajectoryError::WorkBudget { .. })
    ));
    let fixed = EffusionAdaptiveSolve {
        minimum_step_s: 1.0,
        ..accuracy(2e-10)
    };
    assert!(matches!(
        before.advance_effusion_exchange_adaptive(
            quiet(),
            &[],
            &ports,
            1.0,
            endpoint(),
            controls(),
            fixed,
            || false
        ),
        Err(EffusionTrajectoryError::MinimumStep { .. })
    ));
    assert_eq!(before, saved);
}

#[test]
fn accepted_exchange_history_owns_all_gross_net_load_and_conservation_totals() {
    let before = initial(PURE, Vec3::new(2.0, -1.0, 0.5));
    let saved = before.clone();
    let source = ClosedGasParcel::new(
        before.reservoir().model().clone(),
        source::point(600.0, 3.0 * P0).unwrap(),
        VOLUME,
    )
    .unwrap();
    let law = interface(before.reservoir().model());
    let outlets = [outlet(&law, Vec3::new(0.0, 1.0, 0.0))];
    let ports = [port(&source, &law)];
    let loads = EffusionRocketEnvironment {
        gravity_world_m_s2: Vec3::new(0.0, -0.2, 0.0),
        external_force_world_n: Vec3::new(0.001, 0.0, 0.0),
        ..quiet()
    };
    let accuracy = EffusionAdaptiveSolve {
        initial_step_s: 0.1,
        maximum_step_s: 0.1,
        relative_tolerance: 0.0,
        temperature_tolerance_k: 1e-5,
        pressure_tolerance_pa: 0.01,
        species_mass_tolerance_kg: 1e-11,
        position_tolerance_m: 1e-7,
        velocity_tolerance_m_s: 1e-6,
        ..adaptive(1.0)
    };
    let mut polls = 0;
    let result = before
        .advance_effusion_exchange_adaptive(
            loads,
            &outlets,
            &ports,
            0.1,
            endpoint(),
            controls(),
            accuracy,
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    let step = result.step();
    assert_eq!(step.duration_s, 0.1);
    assert!(step.rejected_intervals > 0);
    assert_eq!(step.accepted_steps.len(), 2 * step.accepted_intervals);
    assert_eq!(
        step.attempted_intervals,
        step.accepted_intervals + step.rejected_intervals
    );
    assert!(step.maximum_accepted_error <= 1.0);
    let sum = |field: fn(&EffusionExchangeStep) -> f64| {
        step.accepted_steps.iter().map(field).sum::<f64>()
    };
    let vector = |field: fn(&EffusionExchangeStep) -> Vec3| {
        step.accepted_steps
            .iter()
            .fold(Vec3::ZERO, |sum, row| sum.add(field(row)))
    };
    let scalars: [(f64, fn(&EffusionExchangeStep) -> f64); 13] = [
        (step.total_incoming.mass_kg, |s| s.incoming.mass_kg),
        (step.total_outgoing.mass_kg, |s| s.outgoing.mass_kg),
        (step.total_incoming.relative_energy_j, |s| {
            s.incoming.relative_energy_j
        }),
        (step.total_outgoing.relative_energy_j, |s| {
            s.outgoing.relative_energy_j
        }),
        (step.total_incoming.energy_world_j, |s| {
            s.incoming.energy_world_j
        }),
        (step.total_outgoing.energy_world_j, |s| {
            s.outgoing.energy_world_j
        }),
        (step.total_net_mass_into_gas_kg, |s| s.net_mass_into_gas_kg),
        (step.total_net_relative_energy_into_gas_j, |s| {
            s.net_relative_energy_into_gas_j
        }),
        (step.total_net_energy_into_vessel_world_j, |s| {
            s.net_energy_into_vessel_world_j
        }),
        (step.total_external_work_j, |s| s.mechanics.external_work_j),
        (step.total_mass_residual_kg, |s| {
            s.gas.exchange.total_mass_residual_kg
        }),
        (step.total_energy_residual_j, |s| s.energy_residual_j),
        (step.total_gravity_work_uncertainty_j, |s| {
            s.gravity_work_uncertainty_j
        }),
    ];
    for (total, field) in scalars {
        near(total, sum(field), 2e-12);
    }
    let vectors: [(Vec3, fn(&EffusionExchangeStep) -> Vec3); 8] = [
        (step.total_incoming.relative_momentum_kg_m_s, |s| {
            s.incoming.relative_momentum_kg_m_s
        }),
        (step.total_outgoing.relative_momentum_kg_m_s, |s| {
            s.outgoing.relative_momentum_kg_m_s
        }),
        (step.total_incoming.momentum_world_kg_m_s, |s| {
            s.incoming.momentum_world_kg_m_s
        }),
        (step.total_outgoing.momentum_world_kg_m_s, |s| {
            s.outgoing.momentum_world_kg_m_s
        }),
        (step.total_net_momentum_into_vessel_world_kg_m_s, |s| {
            s.net_momentum_into_vessel_world_kg_m_s
        }),
        (step.total_external_impulse_n_s, |s| {
            s.mechanics.external_impulse_n_s
        }),
        (step.total_momentum_residual_kg_m_s, |s| {
            s.momentum_residual_kg_m_s
        }),
        (step.total_gravity_impulse_uncertainty_n_s, |s| {
            s.gravity_impulse_uncertainty_n_s
        }),
    ];
    for (total, field) in vectors {
        vector_near(total, vector(field), 2e-12);
    }
    let mut replay = before.clone();
    for accepted in &step.accepted_steps {
        let next = replay
            .advance_effusion_exchange(
                loads,
                &outlets,
                &ports,
                accepted.duration_s,
                endpoint(),
                controls(),
                || false,
            )
            .unwrap();
        assert_eq!(next.step(), accepted);
        replay = next.into_parts().0;
    }
    assert_eq!(&replay, result.endpoint());
    near(
        replay.reservoir().mass_kg() - before.reservoir().mass_kg(),
        step.total_net_mass_into_gas_kg,
        1e-13,
    );
    vector_near(
        replay
            .vehicle()
            .velocity_m_s
            .scale(mass(&replay))
            .sub(before.vehicle().velocity_m_s.scale(mass(&before)))
            .sub(step.total_net_momentum_into_vessel_world_kg_m_s)
            .sub(step.total_external_impulse_n_s),
        Vec3::ZERO,
        1e-10,
    );
    near(
        internal(&replay) - internal(&before) + kinetic(&replay)
            - kinetic(&before)
            - step.total_net_energy_into_vessel_world_j
            - step.total_external_work_j,
        0.0,
        2e-7,
    );
    let mut count = 0;
    assert!(matches!(
        before.advance_effusion_exchange_adaptive(
            loads,
            &outlets,
            &ports,
            0.1,
            endpoint(),
            controls(),
            accuracy,
            || {
                count += 1;
                count == polls
            }
        ),
        Err(EffusionTrajectoryError::Cancelled)
    ));
    let foreign = initial("A 1 0.028 3.5 0 1.7e-5 0.025", Vec3::ZERO);
    assert!(matches!(
        before.advance_effusion_exchange_adaptive(
            loads,
            &outlets,
            &[port(foreign.reservoir(), &law)],
            0.1,
            endpoint(),
            controls(),
            accuracy,
            || false
        ),
        Err(EffusionTrajectoryError::Invalid(_))
    ));
    assert_eq!(before, saved);
}

#[test]
fn closed_33_species_gas_preserves_thermochemical_limit_and_nonbinary_final_time() {
    let (mut species, mut transport) = (Vec::new(), Vec::new());
    for i in 0..33 {
        let model = source::load(
            &format!("A{i:02} 1 0.028 2.5 0 1.7e-5 0.025"),
            "original 33-species closed-gas fixture",
        )
        .unwrap();
        species.push(model.thermodynamics().components()[0].clone());
        transport.push(model.transport_materials()[0].clone());
    }
    let mut fractions = vec![1.0 / 64.0; 33];
    fractions[32] = 0.5;
    let model = ThermochemicalGasModel::new(
        FrozenIdealGasMixtureModelV1::new(
            species,
            Composition::new(CompositionBasis::MoleFraction, fractions).unwrap(),
        )
        .unwrap(),
        transport,
        None,
    )
    .unwrap();
    let gas = ClosedGasParcel::new(model, source::point(T0, P0).unwrap(), VOLUME).unwrap();
    let before =
        GasRocketState::new(gas, 0.01, Vec3::ZERO, Vec3::new(3.0, -2.0, 1.0), FRAME).unwrap();
    let foreign = initial(PURE, Vec3::ZERO);
    let law = interface(foreign.reservoir().model());
    let ports = [EffusionReservoirPort {
        area_m2: 0.0,
        ..port(foreign.reservoir(), &law)
    }];
    let accuracy = EffusionAdaptiveSolve {
        initial_step_s: 0.0125,
        maximum_step_s: 0.0125,
        ..adaptive(1.0)
    };
    let old = evolve(&before, &[], 0.1, accuracy);
    let result = exchange(&before, &[], &ports, 0.1, accuracy);
    assert_eq!(result.step().duration_s, 0.1);
    assert_eq!(old.endpoint().reservoir(), before.reservoir());
    assert_eq!(result.endpoint().reservoir(), before.reservoir());
    vector_near(
        result.endpoint().vehicle().position_m,
        before.vehicle().velocity_m_s.scale(0.1),
        2e-14,
    );
    vector_near(
        old.endpoint().vehicle().position_m,
        result.endpoint().vehicle().position_m,
        2e-14,
    );
    assert_eq!(result.step().total_incoming.mass_kg, 0.0);
    assert_eq!(result.step().total_outgoing.mass_kg, 0.0);
    assert_eq!(result.step().rejected_intervals, 0);
}
