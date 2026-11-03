//! G0/G1/G3/G4: refine the original aperture operator, retain accepted accounts.
use super::*;
use crate::molecular_exchange::trajectory::{
    MolecularExchangeAdaptiveControls, MolecularExchangeTrajectoryCandidate,
    MolecularExchangeTrajectoryError, advance_candidate,
};

fn accuracy(scale: f64) -> MolecularExchangeAdaptiveControls {
    MolecularExchangeAdaptiveControls {
        initial_step_s: 0.5,
        minimum_step_s: 1e-6,
        maximum_step_s: 0.5,
        maximum_attempts: 2048,
        maximum_accepted_intervals: 1024,
        relative_tolerance: 1e-3 * scale,
        temperature_tolerance_k: 0.1 * scale,
        pressure_tolerance_pa: 10.0 * scale,
        species_mass_tolerance_kg: 1e-9 * scale,
        velocity_tolerance_m_s: 0.01 * scale,
    }
}

fn evolve(
    before: &[MovingGasParcel; 2],
    law: &EffusiveDriftInterface,
    aperture: MolecularExchangeSurface,
    duration: f64,
    adaptive: MolecularExchangeAdaptiveControls,
) -> MolecularExchangeTrajectoryCandidate {
    advance_candidate(
        [&before[0], &before[1]],
        law,
        aperture,
        duration,
        controls(),
        adaptive,
        || false,
    )
    .unwrap()
}

#[test]
fn tighter_controls_refine_both_finite_gases_and_only_accepted_native_accounts_survive() {
    let model = catalog(PURE, 0.0);
    let law = interface(&model, 1000.0);
    let before = [
        gas(model.clone(), 600.0, 2.0 * P, [120.0, -40.0, 15.0]),
        gas(model, 300.0, 0.8 * P, [-70.0, 80.0, -10.0]),
    ];
    let saved = before.clone();
    let aperture = MolecularExchangeSurface {
        area_m2: 4.0 * AREA,
        ..surface()
    };
    let duration = 0.5;
    let loose = evolve(&before, &law, aperture, duration, accuracy(1.0));
    let tight = evolve(&before, &law, aperture, duration, accuracy(0.01));
    let mut reference = before.clone();
    for _ in 0..512 {
        reference = step_candidate(
            [&reference[0], &reference[1]],
            &law,
            aperture,
            duration / 512.0,
            controls(),
            || false,
        )
        .unwrap()
        .into_parts()
        .0;
    }
    // This native-step reference refines the interface operator only. It does
    // not supply remaining wall tractions or assert isolated-vessel entropy.
    let distance = |states: &[MovingGasParcel; 2]| {
        (0..2)
            .map(|side| {
                let mut value = (temperature(&states[side]) - temperature(&reference[side])).abs()
                    / T
                    + (states[side].parcel().mass_kg() - reference[side].parcel().mass_kg()).abs()
                        / before[side].parcel().mass_kg();
                for d in 0..3 {
                    value += (states[side].velocity_m_s()[d] - reference[side].velocity_m_s()[d])
                        .abs()
                        / (R * T / M).sqrt();
                }
                value
            })
            .sum::<f64>()
    };
    assert!(tight.step().accepted_intervals > loose.step().accepted_intervals);
    assert!(tight.step().rejected_intervals > 0);
    assert!(distance(tight.gases()) < 0.5 * distance(loose.gases()));
    let receipt = tight.step();
    assert_eq!(receipt.duration_s, duration);
    assert_eq!(receipt.accepted_steps.len(), 2 * receipt.accepted_intervals);
    assert_eq!(
        receipt.attempted_intervals,
        receipt.accepted_intervals + receipt.rejected_intervals
    );
    assert!(receipt.maximum_accepted_error <= 1.0);

    let mut masses = [&before[0], &before[1]].map(|g| g.parcel().mass_kg());
    let mut momenta = [&before[0], &before[1]].map(MovingGasParcel::momentum_kg_m_s);
    for native in &receipt.accepted_steps {
        for side in 0..2 {
            let source = if side == 0 {
                &native.exchange.left_to_right
            } else {
                &native.exchange.right_to_left
            };
            vector_near(
                source.relative_bulk_velocity_m_s,
                core::array::from_fn(|d| {
                    momenta[side][d] / masses[side] - aperture.velocity_world_m_s[d]
                }),
                1e-8,
            );
            for d in 0..3 {
                momenta[side][d] += native.gas[side].actual_impulse_n_s[d];
            }
            masses[side] += native.gas[side].caloric.mass_into_gas_kg
                + native.gas[side].caloric.mass_residual_kg;
        }
    }
    let first = &receipt
        .accepted_steps
        .first()
        .unwrap()
        .exchange
        .left_to_right;
    let last = &receipt
        .accepted_steps
        .last()
        .unwrap()
        .exchange
        .left_to_right;
    assert!((first.relative_bulk_velocity_m_s[0] - last.relative_bulk_velocity_m_s[0]).abs() > 1.0);
    for side in 0..2 {
        let gross = receipt
            .accepted_steps
            .iter()
            .flat_map(|s| &s.exchange.species)
            .map(|s| s.mass_kg[side])
            .sum::<f64>();
        near(receipt.totals.gross_mass_kg[side], gross, 1e-17);
        let relative_energy = receipt
            .accepted_steps
            .iter()
            .flat_map(|s| &s.exchange.species)
            .map(|s| s.relative_energy_into_gas_j[side])
            .sum::<f64>();
        near(
            receipt.totals.relative_energy_into_gas_j[side],
            relative_energy,
            1e-11,
        );
    }
    near(
        receipt.totals.net_mass_left_to_right_kg,
        tight.gases()[1].parcel().mass_kg() - before[1].parcel().mass_kg(),
        1e-13,
    );
    vector_near(
        receipt.totals.net_momentum_left_to_right_n_s,
        core::array::from_fn(|d| {
            tight.gases()[1].momentum_kg_m_s()[d] - before[1].momentum_kg_m_s()[d]
        }),
        1e-11,
    );
    near(
        tight
            .gases()
            .iter()
            .map(|g| internal(g) + g.kinetic_energy_j())
            .sum::<f64>()
            - before
                .iter()
                .map(|g| internal(g) + g.kinetic_energy_j())
                .sum::<f64>(),
        0.0,
        2e-6,
    );
    near(
        receipt.totals.energy_residual_j,
        receipt
            .accepted_steps
            .iter()
            .map(|s| s.energy_residual_j)
            .sum(),
        1e-18,
    );
    assert_eq!(before, saved);
}

#[test]
fn both_original_catalog_species_enter_error_control_from_zero_inventory() {
    let model = catalog(MIXTURE, 0.0);
    let law = interface(&model, 1000.0);
    let pure = |fractions| {
        model
            .with_composition(Composition::new(CompositionBasis::MoleFraction, fractions).unwrap())
            .unwrap()
    };
    let before = [
        gas(pure(vec![1.0, 0.0]), T, P, [0.0; 3]),
        gas(pure(vec![0.0, 1.0]), 500.0, 0.8 * P, [0.0, 40.0, -20.0]),
    ];
    let saved = before.clone();
    let a = SpeciesId::new("A").unwrap();
    let b = SpeciesId::new("B").unwrap();
    assert_eq!(before[0].parcel().species_mass_kg(&b).unwrap(), 0.0);
    assert_eq!(before[1].parcel().species_mass_kg(&a).unwrap(), 0.0);
    let adaptive = |mass_tolerance| MolecularExchangeAdaptiveControls {
        relative_tolerance: 0.0,
        temperature_tolerance_k: 1e6,
        pressure_tolerance_pa: 1e9,
        velocity_tolerance_m_s: 1e6,
        species_mass_tolerance_kg: mass_tolerance,
        ..accuracy(1.0)
    };
    let loose = evolve(&before, &law, surface(), 0.5, adaptive(2e-8));
    let tight = evolve(&before, &law, surface(), 0.5, adaptive(2e-10));
    assert!(tight.step().accepted_intervals > loose.step().accepted_intervals);
    assert!(tight.step().rejected_intervals > 0);
    assert!(tight.gases()[0].parcel().species_mass_kg(&b).unwrap() > 0.0);
    assert!(tight.gases()[1].parcel().species_mass_kg(&a).unwrap() > 0.0);
    for id in [&a, &b] {
        near(
            tight
                .gases()
                .iter()
                .map(|g| g.parcel().species_mass_kg(id).unwrap())
                .sum(),
            before
                .iter()
                .map(|g| g.parcel().species_mass_kg(id).unwrap())
                .sum(),
            1e-13,
        );
    }
    assert_eq!(before, saved);
}

#[test]
fn only_duration_withdrawal_and_local_error_retry_and_whole_windows_refuse_atomically() {
    let model = catalog(PURE, 0.0);
    let law = interface(&model, 1000.0);
    let before = [
        gas(model.clone(), T, P, [0.0; 3]),
        gas(model.clone(), T, P, [0.0; 3]),
    ];
    let saved = before.clone();
    let adaptive = MolecularExchangeAdaptiveControls {
        initial_step_s: 0.04,
        maximum_step_s: 0.04,
        ..accuracy(1.0)
    };
    let physical = MolecularExchangeControls {
        maximum_duration_s: 0.015,
        maximum_species_fraction: 1e-4,
        ..controls()
    };
    let duration = 0.04;
    let mut polls = 0;
    let accepted = advance_candidate(
        [&before[0], &before[1]],
        &law,
        surface(),
        duration,
        physical,
        adaptive,
        || {
            polls += 1;
            false
        },
    )
    .unwrap();
    assert!(accepted.step().rejected_intervals > 0);
    assert!(
        accepted
            .step()
            .accepted_steps
            .iter()
            .all(|s| s.duration_s <= 0.0075)
    );
    assert!(matches!(
        advance_candidate(
            [&before[0], &before[1]],
            &law,
            surface(),
            duration,
            physical,
            MolecularExchangeAdaptiveControls {
                maximum_attempts: 1,
                ..adaptive
            },
            || false
        ),
        Err(MolecularExchangeTrajectoryError::WorkBudget {
            elapsed_s: 0.0,
            attempted_intervals: 1,
            ..
        })
    ));
    assert!(
        matches!(advance_candidate([&before[0], &before[1]], &law, surface(), duration, physical,
        MolecularExchangeAdaptiveControls { maximum_accepted_intervals: 1, ..adaptive }, || false),
        Err(MolecularExchangeTrajectoryError::WorkBudget { elapsed_s, accepted_intervals: 1, .. }) if elapsed_s > 0.0)
    );
    let strict = MolecularExchangeAdaptiveControls {
        initial_step_s: 0.25,
        minimum_step_s: 0.25,
        maximum_step_s: 0.25,
        relative_tolerance: 0.0,
        temperature_tolerance_k: 1e-12,
        pressure_tolerance_pa: 1e-8,
        species_mass_tolerance_kg: 1e-20,
        velocity_tolerance_m_s: 1e-14,
        ..adaptive
    };
    assert!(matches!(
        advance_candidate(
            [&before[0], &before[1]],
            &law,
            surface(),
            0.25,
            controls(),
            strict,
            || false
        ),
        Err(MolecularExchangeTrajectoryError::MinimumStep { elapsed_s: 0.0, .. })
    ));
    let bounded = interface(&model, 0.0);
    assert!(
        matches!(advance_candidate([&before[0], &before[1]], &bounded, surface(), duration, controls(), adaptive, || false),
        Err(MolecularExchangeTrajectoryError::Step { elapsed_s: 0.0, duration_s, error: MolecularExchangeError::Material(_) }) if duration_s == duration)
    );
    let mut count = 0;
    assert!(
        matches!(advance_candidate([&before[0], &before[1]], &law, surface(), duration, physical, adaptive, || {
        count += 1; count == polls
    }), Err(MolecularExchangeTrajectoryError::Cancelled { elapsed_s }) if elapsed_s == duration)
    );

    let foreign = interface(&catalog(PURE, 1.0), 0.0);
    let closed = evolve(
        &before,
        &foreign,
        MolecularExchangeSurface {
            area_m2: 0.0,
            normal_left_to_right: [0.0; 3],
            velocity_world_m_s: [f64::NAN; 3],
        },
        duration,
        adaptive,
    );
    assert_eq!(closed.gases(), &before);
    assert_eq!(closed.step().totals.gross_mass_kg, [0.0; 2]);
    assert_eq!(before, saved);
}
