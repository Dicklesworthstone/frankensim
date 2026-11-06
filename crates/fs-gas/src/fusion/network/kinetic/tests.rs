//! G0/G1/G3/G4 finite-time alpha heating, conservation and source-domain checks.

use super::*;
use fs_qty::Temperature;
use fs_thermochem::coulomb::{
    CoulombLimitsV1, CoulombModelV1, CoulombScreeningV1, CoulombShortRangeCutoffV1,
};
use fs_thermochem::fusion::KELVIN_PER_KEV;
use fs_thermochem::fusion::plasma::DtFusionPlasmaLimitsV1;
use fs_thermochem::plasma_radiation::{
    BremsstrahlungControlsV1, BremsstrahlungErrorV1, BremsstrahlungGauntFactorV1,
    ElectronIonBremsstrahlungV1,
};

#[path = "../tests/fixture.rs"]
mod fixture;

mod compression;
mod confinement;
mod history;

const KEV_J: f64 = 1.602_176_634e-16;
const BOLTZMANN: f64 = 1.380_649e-23;
const DT: FusionNetworkChannels = FusionNetworkChannels::new([true, false, false, false]);

fn close(actual: f64, expected: f64, relative: f64, absolute: f64) {
    assert!(
        (actual - expected).abs() <= absolute + relative * expected.abs(),
        "actual={actual:.17e}, expected={expected:.17e}"
    );
}

fn controls() -> FusionNetworkSolve {
    let mut c = FusionNetworkSolve::new(
        Temperature::new(0.5 * KELVIN_PER_KEV),
        Temperature::new(100.0 * KELVIN_PER_KEV),
        Temperature::new(1e4),
        Temperature::new(20.0 * KELVIN_PER_KEV),
        DtFusionPlasmaLimitsV1 {
            maximum_coulomb_coupling: 0.01,
            maximum_thermal_phase_space_density: 0.01,
            maximum_thermal_energy_rest_energy_ratio: 0.05,
        },
    );
    c.relative_tolerance = 1e-6;
    c.absolute_amount_tolerance_mol = 1e-19;
    c.energy_tolerance_j = 1e-7;
    c.initial_step_s = 1e-3;
    c.maximum_step_s = 0.02;
    c.maximum_steps = 10_000;
    c
}

fn grid(maximum_kev: f64) -> AlphaEnergyGrid {
    AlphaEnergyGrid::uniform_speed(
        Energy::new(0.05 * KEV_J),
        Energy::new(maximum_kev * KEV_J),
        128,
        || false,
    )
    .unwrap()
}

fn heating() -> FusionNetworkHeating {
    FusionNetworkHeating {
        product_electron_fractions: [[0.0, 0.0], [0.65, 0.0], [0.4, 0.9], [0.0, 0.9]],
        electron_power: Power::new(0.0),
        heavy_power: Power::new(0.0),
    }
}

fn collision_source(
    reference: &FrozenIonizedGasStateV1,
    projectile: usize,
) -> TestParticleCoulombV1 {
    let background = CoulombModelV1::new_with_cutoff(
        reference,
        CoulombScreeningV1::ElectronDebye,
        CoulombShortRangeCutoffV1::QuantumDiffraction,
        CoulombLimitsV1 {
            minimum_coulomb_logarithm: 4.0,
            minimum_particles_per_debye_sphere: 10.0,
            maximum_coulomb_coupling: 0.01,
            maximum_thermal_phase_space_density: 0.01,
            maximum_thermal_energy_rest_energy_ratio: 0.05,
            maximum_quantum_to_classical_ratio: 0.1,
            maximum_relative_drift_thermal_ratio: 0.1,
        },
        || false,
    )
    .unwrap();
    TestParticleCoulombV1::new(background, projectile, 0.01).unwrap()
}

fn alpha_transport(source: TestParticleCoulombV1) -> KineticAlphaTransport {
    KineticAlphaTransport::new(source, 0.05, 0.25, 1e-4, 100_000).unwrap()
}

fn thermal_transport(
    alpha: &KineticAlphaTransport,
    ion_variance_limit: f64,
) -> FusionNetworkTransport {
    FusionNetworkTransport {
        energy_confinement: None,
        coulomb: Some(alpha.collisions().background().clone()),
        bremsstrahlung: Some(
            ElectronIonBremsstrahlungV1::new(
                alpha
                    .collisions()
                    .background()
                    .initial_state()
                    .thermodynamics()
                    .clone(),
                BremsstrahlungControlsV1 {
                    frequency_integrated_gaunt_factor: BremsstrahlungGauntFactorV1::Constant(1.2),
                    minimum_electron_temperature_k: 1e6,
                    maximum_electron_temperature_k: 20.0 * KELVIN_PER_KEV,
                    maximum_electron_thermal_rest_energy_ratio: 0.05,
                    maximum_ion_to_electron_thermal_variance_ratio: ion_variance_limit,
                    maximum_electron_thermal_phase_space_density: 0.01,
                },
            )
            .unwrap(),
        ),
    }
}

fn scenario(
    alpha_density: f64,
    grid: &AlphaEnergyGrid,
) -> (KineticFusionState, KineticAlphaTransport) {
    let (plasma, reference) =
        fixture::plasma([0.0, 5e19, 5e19, 0.0, alpha_density], 1.0, [10.0, 5.0]);
    let transport = alpha_transport(collision_source(&reference, plasma.indices()[4]));
    let state = if alpha_density == 0.0 {
        KineticFusionState::from_thermal(
            &reference,
            &plasma,
            grid,
            controls().plasma_limits,
            || false,
        )
        .unwrap()
    } else {
        // Author the initial two-cell population independently of the birth
        // injector: exactly the supplied total He4, at mean energy 3.5 MeV.
        let birth_energy = 3500.0 * KEV_J;
        let high = grid
            .energies()
            .partition_point(|e| e.value() < birth_energy);
        let low = high - 1;
        let total = reference.amounts()[plasma.indices()[4]].value();
        let upper = total * (birth_energy - grid.energies()[low].value())
            / (grid.energies()[high].value() - grid.energies()[low].value());
        let mut amounts = vec![Amount::new(0.0); grid.len()];
        amounts[low] = Amount::new(total - upper);
        amounts[high] = Amount::new(upper);
        KineticFusionState::new(
            &reference,
            &plasma,
            grid,
            &amounts,
            controls().plasma_limits,
            || false,
        )
        .unwrap()
    };
    (state, transport)
}

fn evolve(
    initial: &KineticFusionState,
    alpha: &KineticAlphaTransport,
    channels: FusionNetworkChannels,
    duration: f64,
    c: FusionNetworkSolve,
) -> Result<KineticFusionResult, FusionNetworkError> {
    advance_kinetic_fusion(
        initial,
        channels,
        heating(),
        &FusionNetworkTransport::default(),
        alpha,
        Time::new(duration),
        c,
        || false,
    )
}

fn distribution_moments(state: &KineticFusionState) -> (f64, f64, f64) {
    let amount: f64 = state.alpha_populations_mol().iter().sum();
    let energy: f64 = state
        .alpha_populations_mol()
        .iter()
        .zip(state.energy_grid().energies())
        .map(|(n, e)| n * AVOGADRO * e.value())
        .sum();
    if amount == 0.0 {
        return (0.0, energy, 0.0);
    }
    let mean_kev = energy / (amount * AVOGADRO * KEV_J);
    let variance_kev2: f64 = state
        .alpha_populations_mol()
        .iter()
        .zip(state.energy_grid().energies())
        .map(|(n, e)| n * (e.value() / KEV_J - mean_kev).powi(2))
        .sum::<f64>()
        / amount;
    (amount, energy, variance_kev2)
}

fn check_accounts(result: &KineticFusionResult) {
    let (initial, state) = (result.initial_state(), result.state());
    let b = result.balance();
    let a = result.alpha_balance();
    let emission = result.neutron_emission().unwrap();
    assert_eq!(emission.duration(), result.elapsed_time());
    for (i, channel) in emission.channels().iter().enumerate() {
        assert_eq!(channel.reaction(), FusionReactionV1::ALL[i]);
        assert_eq!(channel.amount().value(), b.reaction_extents_mol[i]);
        assert_eq!(
            channel.kinetic_energy().value(),
            b.product_kinetic_energies_j[i][1],
        );
    }
    close(
        emission.total_amount().value(),
        b.neutron_amount_mol,
        2e-13,
        0.0,
    );
    close(
        emission.total_kinetic_energy().value(),
        b.neutron_energy_escaped_j,
        2e-13,
        0.0,
    );
    assert!(
        state
            .alpha_populations_mol()
            .iter()
            .all(|n| n.is_finite() && *n >= 0.0)
    );
    let (amount, alpha_energy, _) = distribution_moments(state);
    close(amount, state.alpha_amount().value(), 3e-12, 1e-22);
    close(alpha_energy, state.alpha_energy().value(), 3e-12, 1e-8);

    let thermal = state.thermal_amounts();
    let ion_amount: f64 = state
        .plasma()
        .indices()
        .iter()
        .map(|&i| thermal[i].value())
        .sum();
    let electron_amount = state.amounts()[state.plasma().electron_index()].value();
    // Independent monatomic fixture calorics: kinetic He4 contributes no Cv.
    let heavy_energy =
        1.5 * BOLTZMANN * AVOGADRO * ion_amount * state.background_heavy_temperature().value();
    let electron_energy =
        1.5 * BOLTZMANN * AVOGADRO * electron_amount * state.electron_temperature().value();
    close(
        state.thermal_heavy_energy().value(),
        heavy_energy,
        3e-12,
        1e-7,
    );
    close(
        state.thermal_electron_energy().value(),
        electron_energy,
        3e-12,
        1e-7,
    );
    close(
        state.pressure().value(),
        2.0 * (heavy_energy + electron_energy + alpha_energy) / (3.0 * state.volume().value()),
        3e-12,
        1e-7,
    );

    close(
        state.alpha_amount().value() - initial.alpha_amount().value(),
        a.born_alpha_amount_mol,
        3e-12,
        1e-22,
    );
    close(
        alpha_energy - initial.alpha_energy().value()
            + a.deposited_heavy_energy_j
            + a.deposited_electron_energy_j,
        a.born_alpha_energy_j + a.boundary_work_on_alphas_j,
        3e-12,
        2e-6,
    );
    close(
        state.thermal_heavy_energy().value() - initial.thermal_heavy_energy().value(),
        a.deposited_heavy_energy_j + b.charged_product_heavy_energy_j
            - b.reacting_ion_energy_removed_j
            + b.electron_to_heavy_coulomb_energy_j
            + b.external_heavy_energy_j
            + b.boundary_work_on_heavy_j
            - b.heavy_transport_energy_out_j,
        3e-12,
        2e-6,
    );
    close(
        state.thermal_electron_energy().value() - initial.thermal_electron_energy().value(),
        a.deposited_electron_energy_j + b.charged_product_electron_energy_j
            - b.electron_to_heavy_coulomb_energy_j
            - b.photon_energy_escaped_j
            + b.external_electron_energy_j
            + b.boundary_work_on_electrons_j
            - b.electron_transport_energy_out_j,
        3e-12,
        2e-6,
    );
    close(
        state.internal_energy().value() - initial.internal_energy().value()
            + b.neutron_energy_escaped_j
            + b.photon_energy_escaped_j
            + b.heavy_transport_energy_out_j
            + b.electron_transport_energy_out_j,
        b.nuclear_energy_released_j
            + b.external_heavy_energy_j
            + b.external_electron_energy_j
            + b.boundary_work_on_heavy_j
            + b.boundary_work_on_electrons_j
            + a.boundary_work_on_alphas_j,
        3e-12,
        2e-6,
    );

    for weights in [[1.0, 1.0, 1.0, 2.0, 2.0], [1.0, 2.0, 3.0, 3.0, 4.0]] {
        let inventory = |s: &KineticFusionState| -> f64 {
            s.plasma()
                .indices()
                .into_iter()
                .zip(weights)
                .map(|(i, w)| w * s.amounts()[i].value())
                .sum()
        };
        let escaped = if weights[4] == 4.0 {
            b.neutron_amount_mol
        } else {
            0.0
        };
        close(inventory(state) + escaped, inventory(initial), 3e-12, 1e-20);
    }
    assert!(b.maximum_normalized_local_error <= 1.0);
    assert_eq!(state.plasma(), initial.plasma());
    if result.volume_history().is_none() {
        assert_eq!(state.energy_grid(), initial.energy_grid());
    } else {
        let factor = (initial.volume().value() / state.volume().value()).powf(2.0 / 3.0);
        assert_eq!(state.energy_grid().len(), initial.energy_grid().len());
        for (actual, original) in state
            .energy_grid()
            .energies()
            .iter()
            .zip(initial.energy_grid().energies())
        {
            close(actual.value(), original.value() * factor, 3e-14, 0.0);
        }
    }
}

// Independent fixed-grid reference: SciPy 1.17 Radau, rtol=2e-12, directly
// integrating alpha populations, both finite thermal reservoirs and reaction
// extents/product energies. Coulomb Psi uses scipy.special.gammainc(3/2,x).
// There are no production reaction maps, uniformization or split substeps.
// Doubling temporal accuracy changes the heat accounts by <2e-8 J.
// The same midpoint-face SG model is intentional: these are time/physics
// coupling checks, not a claim that 128 cells resolve the continuum shape.
#[test]
fn g1_seeded_alphas_slow_and_heat_both_finite_reservoirs() {
    let (initial, alpha) = scenario(1e16, &grid(5000.0));
    let unchanged = initial.clone();
    close(initial.alpha_energy().value(), 5607.618219, 2e-14, 1e-9);
    assert_eq!(
        initial.thermal_amounts()[initial.plasma().indices()[4]].value(),
        0.0
    );
    let result = evolve(
        &initial,
        &alpha,
        FusionNetworkChannels::NONE,
        0.1,
        controls(),
    )
    .unwrap();
    let state = result.state();
    let a = result.alpha_balance();
    close(
        state.background_heavy_temperature().value() / KELVIN_PER_KEV,
        10.006042136213585,
        2e-6,
        1e-9,
    );
    close(
        state.electron_temperature().value() / KELVIN_PER_KEV,
        5.1743821568428165,
        2e-6,
        1e-9,
    );
    close(state.alpha_energy().value(), 1270.706247835396, 3e-5, 1e-5);
    close(a.deposited_heavy_energy_j, 145.2085419124674, 3e-5, 1e-5);
    close(a.deposited_electron_energy_j, 4191.703429252431, 3e-5, 1e-5);
    close(
        distribution_moments(state).2,
        29728.942705058027,
        2e-4,
        1e-3,
    );
    assert_eq!(state.amounts(), initial.amounts());
    assert_eq!(result.balance().reaction_extents_mol, [0.0; 4]);
    assert_eq!(a.born_alpha_amount_mol, 0.0);
    assert_eq!(a.born_alpha_energy_j, 0.0);
    assert!(a.deposited_heavy_energy_j > 0.0 && a.deposited_electron_energy_j > 0.0);
    assert!(
        state
            .alpha_populations_mol()
            .iter()
            .filter(|n| **n > 0.0)
            .count()
            > 2
    );
    assert_eq!(result.alpha_transport(), &alpha);
    assert_eq!(initial, unchanged);
    check_accounts(&result);
}

#[test]
fn g1_g3_dt_birth_retains_alpha_storage_and_restarts_without_thermalizing_it() {
    let (initial, alpha) = scenario(0.0, &grid(5000.0));
    assert!(initial.alpha_populations_mol().iter().all(|n| *n == 0.0));
    let whole = evolve(&initial, &alpha, DT, 0.1, controls()).unwrap();
    let state = whole.state();
    let b = whole.balance();
    let a = whole.alpha_balance();
    close(
        state.background_heavy_temperature().value() / KELVIN_PER_KEV,
        10.00265658781198,
        2e-6,
        1e-9,
    );
    close(
        state.electron_temperature().value() / KELVIN_PER_KEV,
        5.306561481163247,
        2e-6,
        1e-9,
    );
    close(b.reaction_extents_mol[0], 4.713862290577633e-8, 3e-6, 1e-20);
    close(state.alpha_energy().value(), 8651.279136037416, 3e-5, 1e-5);
    close(a.deposited_heavy_energy_j, 178.19523518957058, 3e-5, 1e-5);
    close(a.deposited_electron_energy_j, 7367.484630062573, 3e-5, 1e-5);
    close(
        b.product_kinetic_energies_j[0][0],
        16196.959001289579,
        3e-6,
        1e-5,
    );
    close(b.neutron_energy_escaped_j, 64053.03869618796, 3e-6, 1e-5);
    close(
        b.reacting_ion_energy_removed_j,
        250.83220998177535,
        3e-6,
        1e-5,
    );
    close(state.pressure().value(), 250957.24629420502, 2e-6, 1e-5);
    close(distribution_moments(state).2, 660067.6625228622, 2e-4, 1e-3);
    assert_eq!(b.charged_product_heavy_energy_j, 0.0);
    assert_eq!(b.charged_product_electron_energy_j, 0.0);
    assert_eq!(
        state.thermal_amounts()[state.plasma().indices()[4]].value(),
        0.0
    );
    close(
        state.amounts()[state.plasma().indices()[4]].value(),
        state.alpha_amount().value(),
        3e-12,
        1e-22,
    );
    for &i in &state.plasma().indices()[1..=2] {
        close(
            state.amounts()[i].value(),
            8.297981473578671e-5,
            3e-9,
            1e-20,
        );
    }
    assert_eq!(
        state.amounts()[state.plasma().electron_index()],
        initial.amounts()[initial.plasma().electron_index()]
    );
    check_accounts(&whole);

    let first = evolve(&initial, &alpha, DT, 0.05, controls()).unwrap();
    let second = evolve(first.state(), first.alpha_transport(), DT, 0.05, controls()).unwrap();
    close(second.state().elapsed_time().value(), 0.1, 0.0, 2e-16);
    assert_eq!(second.elapsed_time(), Time::new(0.05));
    close(
        second.state().alpha_energy().value(),
        state.alpha_energy().value(),
        3e-5,
        1e-5,
    );
    close(
        second.state().electron_temperature().value(),
        state.electron_temperature().value(),
        2e-6,
        1e-3,
    );
    close(
        first.alpha_balance().deposited_electron_energy_j
            + second.alpha_balance().deposited_electron_energy_j,
        a.deposited_electron_energy_j,
        3e-5,
        1e-5,
    );
    close(
        first.balance().reaction_extents_mol[0] + second.balance().reaction_extents_mol[0],
        b.reaction_extents_mol[0],
        3e-6,
        1e-20,
    );
    let l1: f64 = second
        .state()
        .alpha_populations_mol()
        .iter()
        .zip(state.alpha_populations_mol())
        .map(|(a, b)| (a - b).abs())
        .sum();
    assert!(l1 / state.alpha_amount().value() < 2e-4);
    assert_eq!(second.alpha_transport(), &alpha);
    check_accounts(&first);
    check_accounts(&second);
}

#[test]
fn g1_fixed_grid_time_refinement_approaches_independent_birth_oracle() {
    let (initial, alpha) = scenario(0.0, &grid(5000.0));
    let run = |step| {
        let mut c = controls();
        c.relative_tolerance = 0.05;
        c.absolute_amount_tolerance_mol = 1e-8;
        c.absolute_temperature_tolerance_k = 1e5;
        c.energy_tolerance_j = 1000.0;
        c.initial_step_s = step;
        c.maximum_step_s = step;
        evolve(&initial, &alpha, DT, 0.1, c).unwrap()
    };
    let coarse = run(0.01);
    let fine = run(0.005);
    assert_eq!(coarse.balance().rejected_steps, 0);
    assert_eq!(fine.balance().rejected_steps, 0);
    let error =
        |r: &KineticFusionResult| (r.state().alpha_energy().value() - 8651.279136037416).abs();
    assert!(error(&coarse) > 1e-5);
    let ratio = error(&fine) / error(&coarse);
    assert!(
        (0.1..0.4).contains(&ratio),
        "time-refinement error ratio={ratio}"
    );
    check_accounts(&coarse);
    check_accounts(&fine);
}

#[test]
fn g4_hybrid_refuses_source_substitution_dilution_and_invalid_birth_domains() {
    let (initial, alpha) = scenario(1e16, &grid(5000.0));
    let unchanged = initial.clone();
    let (other_plasma, other) =
        fixture::plasma_with_helium3_limit([0.0, 5e19, 5e19, 0.0, 1e16], 1.0, [10.0, 5.0], 1e9);
    let substituted = alpha_transport(collision_source(&other, other_plasma.indices()[4]));
    let wrong_projectile = alpha_transport(collision_source(
        alpha.collisions().background().initial_state(),
        initial.plasma().indices()[1],
    ));
    for source in [&substituted, &wrong_projectile] {
        assert!(matches!(
            evolve(
                &initial,
                source,
                FusionNetworkChannels::NONE,
                0.0,
                controls()
            ),
            Err(FusionNetworkError::Invalid(_))
        ));
    }
    for (number_limit, energy_limit) in [(5e-5, 0.25), (0.05, 0.01)] {
        let limited = KineticAlphaTransport::new(
            alpha.collisions().clone(),
            number_limit,
            energy_limit,
            1e-4,
            100_000,
        )
        .unwrap();
        assert_eq!(
            evolve(
                &initial,
                &limited,
                FusionNetworkChannels::NONE,
                0.0,
                controls()
            )
            .unwrap_err(),
            FusionNetworkError::OutsideDomain
        );
    }
    for channel in [0, 3] {
        let mut prompt = heating();
        prompt.product_electron_fractions[channel][0] = 0.5;
        assert!(matches!(
            advance_kinetic_fusion(
                &initial,
                DT,
                prompt,
                &FusionNetworkTransport::default(),
                &alpha,
                Time::new(0.0),
                controls(),
                || false
            ),
            Err(FusionNetworkError::Invalid(_))
        ));
    }

    let mut edge_populations = vec![Amount::new(0.0); initial.energy_grid().len()];
    *edge_populations.last_mut().unwrap() = initial.alpha_amount();
    let edge = KineticFusionState::new(
        alpha.collisions().background().initial_state(),
        initial.plasma(),
        initial.energy_grid(),
        &edge_populations,
        controls().plasma_limits,
        || false,
    )
    .unwrap();
    assert_eq!(
        evolve(&edge, &alpha, FusionNetworkChannels::NONE, 0.0, controls()).unwrap_err(),
        FusionNetworkError::OutsideDomain
    );

    let (empty, empty_alpha) = scenario(0.0, &grid(3000.0));
    let mut bounded = controls();
    bounded.minimum_step_s = 1e-4;
    assert_eq!(
        evolve(&empty, &empty_alpha, DT, 0.01, bounded).unwrap_err(),
        FusionNetworkError::OutsideDomain
    );
    let (too_energetic, too_energetic_source) = scenario(0.0, &grid(50_000.0));
    assert!(matches!(
        evolve(&too_energetic, &too_energetic_source, DT, 0.0, controls()),
        Err(FusionNetworkError::Invalid(_))
    ));

    let mut polls = 0;
    let cancellation = advance_kinetic_fusion(
        &initial,
        FusionNetworkChannels::NONE,
        heating(),
        &FusionNetworkTransport::default(),
        &alpha,
        Time::new(0.1),
        controls(),
        || {
            polls += 1;
            polls >= 2000
        },
    );
    assert_eq!(cancellation.unwrap_err(), FusionNetworkError::Cancelled);
    assert_eq!(polls, 2000);
    assert_eq!(initial, unchanged);
}

#[test]
fn g0_no_alpha_population_matches_thermal_owner_under_external_heating() {
    let (initial, alpha) = scenario(0.0, &grid(5000.0));
    let mut held = heating();
    held.heavy_power = Power::new(200.0);
    held.electron_power = Power::new(-100.0);
    let result = advance_kinetic_fusion(
        &initial,
        FusionNetworkChannels::NONE,
        held,
        &FusionNetworkTransport::default(),
        &alpha,
        Time::new(0.02),
        controls(),
        || false,
    )
    .unwrap();
    let thermal = advance_fusion_network(
        alpha.collisions().background().initial_state(),
        initial.plasma(),
        FusionNetworkChannels::NONE,
        held,
        &FusionNetworkTransport::default(),
        Time::new(0.02),
        controls(),
        || false,
    )
    .unwrap();
    assert_eq!(result.state().alpha_amount().value(), 0.0);
    assert_eq!(result.state().alpha_energy().value(), 0.0);
    assert!(
        result
            .state()
            .alpha_populations_mol()
            .iter()
            .all(|n| *n == 0.0)
    );
    assert_eq!(result.state().amounts(), thermal.state().amounts());
    close(
        result.state().background_heavy_temperature().value(),
        thermal.state().heavy_temperature().value(),
        2e-13,
        1e-5,
    );
    close(
        result.state().electron_temperature().value(),
        thermal.state().electron_temperature().value(),
        2e-13,
        1e-5,
    );
    close(
        result.state().pressure().value(),
        thermal.state().pressure().value(),
        2e-13,
        1e-7,
    );
    check_accounts(&result);
}

#[test]
fn g4_collision_work_exhaustion_retries_shorter_steps_within_the_network_budget() {
    let (initial, unrestricted) = scenario(1e16, &grid(5000.0));
    let unchanged = initial.clone();
    let limited =
        KineticAlphaTransport::new(unrestricted.collisions().clone(), 0.05, 0.25, 1e-4, 30)
            .unwrap();
    let mut c = controls();
    // Loose time tolerances isolate the per-collision work cap: the same
    // initial interval is accepted in one attempt when collision work is ample.
    c.relative_tolerance = 0.01;
    c.absolute_amount_tolerance_mol = 1e-12;
    c.initial_step_s = 0.004;
    c.maximum_step_s = 0.004;
    c.minimum_step_s = 1e-8;
    c.maximum_steps = 500;
    let mut one_attempt = c;
    one_attempt.maximum_steps = 1;
    let baseline = evolve(
        &initial,
        &unrestricted,
        FusionNetworkChannels::NONE,
        0.004,
        one_attempt,
    )
    .unwrap();
    assert_eq!(baseline.balance().accepted_steps, 1);
    assert_eq!(baseline.balance().rejected_steps, 0);
    assert_eq!(
        evolve(
            &initial,
            &limited,
            FusionNetworkChannels::NONE,
            0.004,
            one_attempt
        )
        .unwrap_err(),
        FusionNetworkError::StepBudget
    );

    let recovered = evolve(&initial, &limited, FusionNetworkChannels::NONE, 0.004, c).unwrap();
    assert!(recovered.balance().rejected_steps > 0);
    assert!(recovered.balance().accepted_steps > 1);
    assert_eq!(
        recovered.alpha_transport().maximum_matrix_applications(),
        30
    );
    close(recovered.state().elapsed_time().value(), 0.004, 0.0, 1e-17);
    close(
        recovered.state().alpha_energy().value(),
        baseline.state().alpha_energy().value(),
        3e-5,
        1e-5,
    );
    close(
        recovered.alpha_balance().deposited_electron_energy_j,
        baseline.alpha_balance().deposited_electron_energy_j,
        3e-5,
        1e-5,
    );
    check_accounts(&baseline);
    check_accounts(&recovered);

    // One application cannot resolve the Poisson tail before the declared
    // minimum time step is reached; recovery must remain bounded and refuse.
    let exhausted =
        KineticAlphaTransport::new(unrestricted.collisions().clone(), 0.05, 0.25, 1e-4, 1).unwrap();
    assert_eq!(
        evolve(&initial, &exhausted, FusionNetworkChannels::NONE, 0.004, c).unwrap_err(),
        FusionNetworkError::StepBudget
    );
    assert_eq!(initial, unchanged);
}

#[test]
fn g1_all_channels_share_alpha_storage_coulomb_exchange_and_radiation_accounts() {
    let (plasma, reference) = fixture::plasma([2e19, 5e19, 5e19, 2e19, 0.0], 1.0, [10.0, 5.0]);
    let original_reference = reference.clone();
    let alpha = alpha_transport(collision_source(&reference, plasma.indices()[4]));
    let thermal = thermal_transport(&alpha, 0.1);
    let initial = KineticFusionState::from_thermal(
        &reference,
        &plasma,
        &grid(5000.0),
        controls().plasma_limits,
        || false,
    )
    .unwrap();
    let result = advance_kinetic_fusion(
        &initial,
        FusionNetworkChannels::ALL,
        heating(),
        &thermal,
        &alpha,
        Time::new(0.002),
        controls(),
        || false,
    )
    .unwrap();
    let b = result.balance();
    let a = result.alpha_balance();
    assert!(b.reaction_extents_mol.iter().all(|extent| *extent > 0.0));
    assert!(
        b.product_kinetic_energies_j
            .iter()
            .flatten()
            .all(|energy| *energy > 0.0)
    );
    close(
        a.born_alpha_amount_mol,
        b.reaction_extents_mol[0] + b.reaction_extents_mol[3],
        3e-12,
        1e-22,
    );
    close(
        result.state().alpha_amount().value(),
        a.born_alpha_amount_mol,
        3e-12,
        1e-22,
    );
    close(
        a.born_alpha_energy_j,
        b.product_kinetic_energies_j[0][0] + b.product_kinetic_energies_j[3][0],
        3e-12,
        1e-8,
    );
    assert_eq!(
        result.state().thermal_amounts()[plasma.indices()[4]].value(),
        0.0
    );

    // Only DD He3/T/protons and the DHe3 proton deposit promptly. Both alpha
    // products enter the common kinetic distribution with their full energy.
    let products = b.product_kinetic_energies_j;
    let prompt_electron =
        0.65 * products[1][0] + 0.4 * products[2][0] + 0.9 * products[2][1] + 0.9 * products[3][1];
    let prompt_heavy =
        0.35 * products[1][0] + 0.6 * products[2][0] + 0.1 * products[2][1] + 0.1 * products[3][1];
    close(
        b.charged_product_electron_energy_j,
        prompt_electron,
        3e-12,
        1e-8,
    );
    close(b.charged_product_heavy_energy_j, prompt_heavy, 3e-12, 1e-8);
    assert!(a.deposited_heavy_energy_j > 0.0 && a.deposited_electron_energy_j > 0.0);
    assert!(b.electron_to_heavy_coulomb_energy_j < 0.0);
    assert!(b.photon_energy_escaped_j > 0.0);
    assert_eq!(result.thermal_transport(), &thermal);
    assert_eq!(result.alpha_transport(), &alpha);
    assert_eq!(result.state().plasma(), &plasma);
    assert_eq!(reference, original_reference);
    check_accounts(&result);

    let (seeded, seed_alpha) = scenario(1e16, &grid(5000.0));
    let strict = thermal_transport(&seed_alpha, 0.05);
    // The thermal radiation model admits these temperatures. The hybrid owner
    // additionally rejects the fastest represented alpha against this limit.
    strict
        .bremsstrahlung
        .as_ref()
        .unwrap()
        .evaluate(seed_alpha.collisions().background().initial_state(), || {
            false
        })
        .unwrap();
    assert!(matches!(
        advance_kinetic_fusion(
            &seeded,
            FusionNetworkChannels::NONE,
            heating(),
            &strict,
            &seed_alpha,
            Time::new(0.0),
            controls(),
            || false,
        ),
        Err(FusionNetworkError::Radiation(
            BremsstrahlungErrorV1::Domain(_)
        ))
    ));
}
