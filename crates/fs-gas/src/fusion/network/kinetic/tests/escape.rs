//! G0/G1/G3/G4 selective alpha escape, physical-energy sampling and full accounts.

use super::*;
use crate::fusion::kinetic::AlphaEscapeTable;
use crate::fusion::volume::{FusionVolumeHistory, FusionVolumeKnot};
use crate::fusion::{FusionEnergyConfinement, FusionParticleExchange};
use fs_qty::{Frequency, Volume};

fn constant(rate: f64) -> AlphaEscapeTable {
    AlphaEscapeTable::new(
        "manufactured constant alpha escape".into(),
        vec![Energy::new(0.001 * KEV_J), Energy::new(10_000.0 * KEV_J)],
        vec![Frequency::new(rate); 2],
        || false,
    )
    .unwrap()
}

#[test]
fn g1_escape_exact_energy_maps_retain_weak_losses_and_deep_survivors() {
    let energies = [1.0, 5.5, 10.0, 100.0].map(|e| Energy::new(e * KEV_J));
    let grid = AlphaEnergyGrid::new(energies.to_vec(), || false).unwrap();
    let table = AlphaEscapeTable::new(
        "manufactured selective loss".into(),
        [1.0, 10.0, 100.0, 1000.0]
            .map(|e| Energy::new(e * KEV_J))
            .to_vec(),
        [1e-16, 2.0, 120.0, 120.0].map(Frequency::new).to_vec(),
        || false,
    )
    .unwrap();
    close(
        table.hazard_at(energies[1]).unwrap().value(),
        1.0,
        3e-15,
        0.0,
    );
    assert_eq!(table.data_id(), "manufactured selective loss");
    assert_eq!(table.energies().len(), table.hazards().len());
    let mut before = [0.0; MAX_ALPHA_BINS];
    before[..4].copy_from_slice(&[1e-5, 2e-5, 3e-5, 4e-5]);
    let dt = 0.3;
    let step = table
        .prepare(&grid, 1.0, &mut || false)
        .unwrap()
        .advance(&grid, 1.0, &before, dt, &mut || false)
        .unwrap();
    let mut lost = 0.0;
    let mut exported = 0.0;
    for (i, hazard) in [1e-16_f64, 1.0, 2.0, 120.0].into_iter().enumerate() {
        let expected = before[i] * (-hazard * dt).exp();
        let removed = before[i] * -(-hazard * dt).exp_m1();
        close(step.amounts[i], expected, 4e-14, 0.0);
        lost += removed;
        exported += removed * AVOGADRO * energies[i].value();
    }
    assert!(step.amounts[3] > 0.0); // exp(-36), not subtraction of nearly all stock
    close(step.removed_mol, lost, 4e-14, 0.0);
    close(step.removed_energy_j, exported, 4e-14, 1e-8);
    close(step.retained_mol, step.amounts.iter().sum(), 3e-15, 0.0);
    // A tiny removal remains a nonzero direct account even below one ULP of n.
    before[1..].fill(0.0);
    let weak = table
        .prepare(&grid, 1.0, &mut || false)
        .unwrap()
        .advance(&grid, 1.0, &before, dt, &mut || false)
        .unwrap();
    close(weak.removed_mol, before[0] * 3e-17, 3e-15, 0.0);
    assert_eq!(weak.amounts[0], before[0]);
    assert!(weak.removed_energy_j > 0.0);
}

#[test]
fn g1_escape_constant_hazard_preserves_exact_number_with_collisional_slowing() {
    let (initial, closed) = scenario(1e16, &grid(5000.0));
    let alpha = closed.clone().with_escape(constant(50.0));
    let duration = 0.01;
    let result = evolve(
        &initial,
        &alpha,
        FusionNetworkChannels::NONE,
        duration,
        controls(),
    )
    .unwrap();
    let before = initial.alpha_amount().value();
    let expected = before * (-50.0_f64 * duration).exp();
    let lost = before * -(-50.0_f64 * duration).exp_m1();
    close(
        result.state().alpha_amount().value(),
        expected,
        5e-12,
        1e-22,
    );
    close(
        result.alpha_balance().escaped_alpha_amount_mol,
        lost,
        5e-12,
        1e-22,
    );
    close(
        result.alpha_balance().escaped_electrons_mol,
        2.0 * lost,
        5e-12,
        1e-22,
    );
    let electron = initial.plasma().electron_index();
    close(
        result.state().amounts()[electron].value(),
        initial.amounts()[electron].value() - 2.0 * lost,
        3e-13,
        1e-20,
    );
    assert!(result.alpha_balance().escaped_alpha_energy_j > 0.0);
    assert!(result.alpha_balance().escaped_electron_thermal_energy_j > 0.0);
    assert!(result.alpha_balance().deposited_electron_energy_j > 0.0);
    assert_eq!(
        result.balance().particles,
        crate::fusion::FusionParticleBalance::default()
    );
    assert_eq!(result.alpha_transport().escape(), alpha.escape());
    check_accounts(&result);

    let no_loss = evolve(
        &initial,
        &closed,
        FusionNetworkChannels::NONE,
        duration,
        controls(),
    )
    .unwrap();
    assert!(
        result.alpha_balance().deposited_electron_energy_j
            < no_loss.alpha_balance().deposited_electron_energy_j
    );
    let zero = evolve(
        &initial,
        &closed.with_escape(constant(0.0)),
        FusionNetworkChannels::NONE,
        duration,
        controls(),
    )
    .unwrap();
    assert_eq!(zero.state(), no_loss.state());
    assert_eq!(zero.balance(), no_loss.balance());
    assert_eq!(zero.alpha_balance(), no_loss.alpha_balance());
}

#[test]
fn g1_escape_samples_physical_energy_under_compression_against_analytic_ode() {
    // A low-density plasma makes collisional transfer negligible, with that
    // error checked below. The independent source ODE has E=E0 exp(bt) and
    // lambda(E)=lambda0 E/E0, hence n=n0 exp[-lambda0 expm1(bt)/b].
    let (plasma, reference) = fixture::plasma([0.0, 5e7, 5e7, 0.0, 1e4], 1.0, [10.0, 5.0]);
    let grid = grid(5000.0);
    let index = 80;
    let initial_energy = grid.energies()[index].value();
    let initial_amount = reference.amounts()[plasma.indices()[4]];
    let mut populations = vec![Amount::new(0.0); grid.len()];
    populations[index] = initial_amount;
    let initial = KineticFusionState::new(
        &reference,
        &plasma,
        &grid,
        &populations,
        controls().plasma_limits,
        || false,
    )
    .unwrap();
    let lambda0 = 20.0;
    let support = [0.001 * KEV_J, 10_000.0 * KEV_J];
    let law = AlphaEscapeTable::new(
        "manufactured lambda proportional to physical energy".into(),
        support.map(Energy::new).to_vec(),
        support
            .map(|e| Frequency::new(lambda0 * e / initial_energy))
            .to_vec(),
        || false,
    )
    .unwrap();
    let alpha = alpha_transport(collision_source(&reference, plasma.indices()[4])).with_escape(law);
    let duration = 0.02;
    let ratio = 0.5;
    let geometry = FusionVolumeHistory::new(
        vec![
            FusionVolumeKnot {
                time: Time::new(0.0),
                volume: Volume::new(1.0),
            },
            FusionVolumeKnot {
                time: Time::new(duration),
                volume: Volume::new(ratio),
            },
        ],
        || false,
    )
    .unwrap();
    let mut c = controls();
    c.relative_tolerance = 2e-8;
    c.absolute_amount_tolerance_mol = 1e-32;
    c.energy_tolerance_j = 1e-18;
    let result = advance_kinetic_fusion_with_volume(
        &initial,
        FusionNetworkChannels::NONE,
        heating(),
        &FusionNetworkTransport::default(),
        &alpha,
        &geometry,
        c,
        || false,
    )
    .unwrap();
    let b = -(2.0 / 3.0) * ratio.ln() / duration;
    let scale = (b * duration).exp();
    let exposure = lambda0 * (b * duration).exp_m1() / b;
    let survivor = (-exposure).exp();
    let loss = -(-exposure).exp_m1();
    let n0 = initial_amount.value();
    let u0 = initial.alpha_energy().value();
    let expected_work = u0 * (b / lambda0) * loss;
    let expected_export = u0 * ((1.0 + b / lambda0) * loss - (scale - 1.0) * survivor);
    let a = result.alpha_balance();
    close(
        result.state().alpha_amount().value(),
        n0 * survivor,
        4e-7,
        1e-32,
    );
    close(
        result.state().alpha_energy().value(),
        u0 * scale * survivor,
        4e-7,
        1e-19,
    );
    close(a.escaped_alpha_amount_mol, n0 * loss, 4e-7, 1e-32);
    close(a.escaped_alpha_energy_j, expected_export, 8e-7, 1e-19);
    close(a.boundary_work_on_alphas_j, expected_work, 8e-7, 1e-19);
    let electron_energy_per_mol =
        1.5 * BOLTZMANN * AVOGADRO * initial.electron_temperature().value();
    close(
        a.escaped_electron_thermal_energy_j,
        expected_export * (2.0 * electron_energy_per_mol / (AVOGADRO * initial_energy)),
        8e-7,
        1e-21,
    );
    for (actual, original) in [
        (
            result.state().background_heavy_temperature(),
            initial.background_heavy_temperature(),
        ),
        (
            result.state().electron_temperature(),
            initial.electron_temperature(),
        ),
    ] {
        close(actual.value(), original.value() * scale, 2e-10, 1e-5);
    }
    assert!(a.deposited_heavy_energy_j.abs() + a.deposited_electron_energy_j.abs() < 1e-9 * u0);
    assert!(result.balance().total_energy_residual_j.abs() < 1e-18);
    check_accounts(&result);
}

#[test]
fn g3_escape_competes_with_births_fueling_and_thermal_exhaust_and_restarts() {
    let (plasma, reference) = fixture::plasma([2e19, 3e20, 7e19, 4e19, 1e19], 0.025, [15.0, 4.0]);
    let initial = KineticFusionState::from_thermal(
        &reference,
        &plasma,
        &grid(5000.0),
        controls().plasma_limits,
        || false,
    )
    .unwrap();
    let alpha = alpha_transport(collision_source(&reference, plasma.indices()[4]))
        .with_escape(constant(80.0));
    let mut thermal = thermal_transport(&alpha, 0.2);
    thermal.particles = Some(
        FusionParticleExchange::new(
            [0.0, 2e-5, 1e-5, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 10.0].map(Frequency::new),
            Temperature::new(2.0 * KELVIN_PER_KEV),
            Temperature::new(KELVIN_PER_KEV),
        )
        .unwrap(),
    );
    thermal.energy_confinement = Some(
        FusionEnergyConfinement::new(
            Time::new(0.2),
            Time::new(0.13),
            Temperature::new(2.0 * KELVIN_PER_KEV),
            Temperature::new(KELVIN_PER_KEV),
        )
        .unwrap(),
    );
    let run = |state, duration| {
        advance_kinetic_fusion(
            state,
            FusionNetworkChannels::ALL,
            heating(),
            &thermal,
            &alpha,
            Time::new(duration),
            controls(),
            || false,
        )
        .unwrap()
    };
    let result = run(&initial, 0.006);
    let half = run(&initial, 0.003);
    let resumed = run(half.state(), 0.003);
    let b = result.balance();
    let a = result.alpha_balance();
    assert!(b.reaction_extents_mol.iter().all(|x| *x > 0.0));
    assert!(a.born_alpha_amount_mol > a.escaped_alpha_amount_mol);
    assert!(a.escaped_alpha_amount_mol > 0.0 && a.escaped_alpha_energy_j > 0.0);
    assert!(b.particles.exhausted_ions_mol[4] > 0.0);
    assert!(b.photon_energy_escaped_j > 0.0 && b.electron_to_heavy_coulomb_energy_j != 0.0);
    assert!(b.heavy_transport_energy_out_j > 0.0 && b.electron_transport_energy_out_j > 0.0);
    for solved in [&result, &half, &resumed] {
        check_accounts(solved);
    }
    for (actual, expected) in resumed
        .state()
        .amounts()
        .iter()
        .zip(result.state().amounts())
    {
        close(actual.value(), expected.value(), 3e-6, 1e-19);
    }
    close(
        resumed.state().alpha_amount().value(),
        result.state().alpha_amount().value(),
        3e-6,
        1e-20,
    );
    close(
        resumed.state().alpha_energy().value(),
        result.state().alpha_energy().value(),
        4e-6,
        1e-6,
    );
    close(
        resumed.state().electron_temperature().value(),
        result.state().electron_temperature().value(),
        3e-6,
        1e-3,
    );
    close(
        half.alpha_balance().escaped_alpha_energy_j
            + resumed.alpha_balance().escaped_alpha_energy_j,
        a.escaped_alpha_energy_j,
        5e-6,
        1e-6,
    );
}

#[test]
fn g0_escape_invalid_support_work_exhaustion_and_late_cancellation_are_atomic() {
    for (energies, hazards) in [
        (vec![Energy::new(1.0)], vec![Frequency::new(1.0)]),
        (
            vec![Energy::new(2.0), Energy::new(1.0)],
            vec![Frequency::new(1.0); 2],
        ),
        (vec![Energy::new(1.0); 2], vec![Frequency::new(1.0); 2]),
        (
            vec![Energy::new(1.0), Energy::new(f64::NAN)],
            vec![Frequency::new(1.0); 2],
        ),
        (
            vec![Energy::new(1.0), Energy::new(2.0)],
            vec![Frequency::new(-1.0); 2],
        ),
        (
            vec![Energy::new(1.0), Energy::new(2.0)],
            vec![Frequency::new(f64::INFINITY); 2],
        ),
    ] {
        assert!(AlphaEscapeTable::new("invalid".into(), energies, hazards, || false).is_err());
    }
    assert!(constant(1.0).hazard_at(Energy::new(0.0)).is_err());
    assert!(
        constant(1.0)
            .hazard_at(Energy::new(20_000.0 * KEV_J))
            .is_err()
    );
    let (initial, closed) = scenario(1e16, &grid(5000.0));
    let original = initial.clone();
    let unsupported = AlphaEscapeTable::new(
        "unsupported".into(),
        vec![Energy::new(1.0 * KEV_J), Energy::new(1000.0 * KEV_J)],
        vec![Frequency::new(1.0); 2],
        || false,
    )
    .unwrap();
    assert_eq!(
        evolve(
            &initial,
            &closed.clone().with_escape(unsupported),
            DT,
            0.0,
            controls()
        ),
        Err(FusionNetworkError::OutsideDomain),
    );
    let alpha = closed.with_escape(constant(50.0));
    let mut c = controls();
    c.maximum_steps = 1;
    assert!(evolve(&initial, &alpha, DT, 0.01, c).is_err());
    let mut calls = 0;
    advance_kinetic_fusion(
        &initial,
        DT,
        heating(),
        &FusionNetworkTransport::default(),
        &alpha,
        Time::new(0.001),
        controls(),
        || {
            calls += 1;
            false
        },
    )
    .unwrap();
    let mut cancelled_calls = 0;
    assert_eq!(
        advance_kinetic_fusion(
            &initial,
            DT,
            heating(),
            &FusionNetworkTransport::default(),
            &alpha,
            Time::new(0.001),
            controls(),
            || {
                cancelled_calls += 1;
                cancelled_calls == calls
            },
        ),
        Err(FusionNetworkError::Cancelled),
    );
    assert_eq!(initial, original);
}

#[test]
fn g1_escape_stiff_births_refuse_saturated_trials_under_a_declared_work_cap() {
    // Collision/burn feedback is negligible at this density. A continuous
    // constant birth B and loss lambda have n=B/lambda*(1-exp(-lambda*t)),
    // even when separate midpoint injection/removal maps all underflow to zero.
    let (plasma, reference) = fixture::plasma([0.0, 5e7, 5e7, 0.0, 0.0], 1.0, [10.0, 5.0]);
    let initial = KineticFusionState::from_thermal(
        &reference,
        &plasma,
        &grid(5000.0),
        controls().plasma_limits,
        || false,
    )
    .unwrap();
    let original = initial.clone();
    let lambda = 1e9;
    let duration = 0.001;
    let alpha = alpha_transport(collision_source(&reference, plasma.indices()[4]))
        .with_escape(constant(lambda));
    let mut c = controls();
    c.initial_step_s = duration;
    c.maximum_step_s = duration;
    c.maximum_steps = 1;
    c.absolute_amount_tolerance_mol = 1e-50;
    let [_, d, t, _, _] = plasma.indices();
    let birth_rate = AVOGADRO
        * reference.amounts()[d].value()
        * reference.amounts()[t].value()
        * FusionReactionV1::DeuteriumTritium
            .reactivity(reference.heavy_temperature())
            .unwrap()
            .value_m3_s()
        / reference.volume().value();
    let expected_stock = birth_rate / lambda * -(-lambda * duration).exp_m1();
    assert!(expected_stock > 1e6 * c.absolute_amount_tolerance_mol);
    assert_eq!((-lambda * duration / 4.0).exp(), 0.0);
    assert_eq!(
        evolve(&initial, &alpha, DT, duration, c),
        Err(FusionNetworkError::StepBudget),
    );
    // There is no escape/collision/birth competition in an empty, source-off
    // population; that exact identity remains usable with the same work cap.
    let empty = evolve(&initial, &alpha, FusionNetworkChannels::NONE, duration, c).unwrap();
    assert_eq!(empty.balance().accepted_steps, 1);
    assert_eq!(empty.balance().rejected_steps, 0);
    assert_eq!(empty.state().alpha_amount().value(), 0.0);
    assert_eq!(empty.alpha_balance(), KineticAlphaBalance::default());
    assert_eq!(initial, original);
}
