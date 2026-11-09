//! G0/G1/G3/G4 coupled reaction-network and original-source regressions.

use super::*;
use fs_qty::Temperature;
use fs_thermochem::coulomb::{
    CoulombLimitsV1, CoulombModelV1, CoulombScreeningV1, CoulombShortRangeCutoffV1,
};
use fs_thermochem::fusion::plasma::DtFusionPlasmaLimitsV1;
use fs_thermochem::fusion::{KELVIN_PER_KEV, SPEED_OF_LIGHT_M_S};
use fs_thermochem::plasma_radiation::{
    BremsstrahlungControlsV1, BremsstrahlungGauntFactorV1, ElectronIonBremsstrahlungV1,
};

mod compression;
mod confinement;
mod fixture;
mod history;
mod particles;

const MIXED_DENSITY: [f64; 5] = [2e19, 3e20, 7e19, 4e19, 1e19];

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
    c.relative_tolerance = 1e-9;
    c.absolute_amount_tolerance_mol = 1e-20;
    c.energy_tolerance_j = 1e-8;
    c.initial_step_s = 1e-3;
    c.maximum_step_s = 0.1;
    c
}

fn heating() -> FusionNetworkHeating {
    FusionNetworkHeating::uniform_charged_fraction(0.8, Power::new(0.0), Power::new(0.0))
}

fn mixed_heating() -> FusionNetworkHeating {
    FusionNetworkHeating {
        product_electron_fractions: [[0.8, 0.0], [0.65, 0.0], [0.4, 0.9], [0.75, 0.9]],
        ..heating()
    }
}

fn transport(state: &FrozenIonizedGasStateV1) -> FusionNetworkTransport {
    FusionNetworkTransport {
        particles: None,
        energy_confinement: None,
        coulomb: Some(
            CoulombModelV1::new_with_cutoff(
                state,
                CoulombScreeningV1::StaticAllChargedDebye,
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
            .unwrap(),
        ),
        bremsstrahlung: Some(
            ElectronIonBremsstrahlungV1::new(
                state.thermodynamics().clone(),
                BremsstrahlungControlsV1 {
                    frequency_integrated_gaunt_factor: BremsstrahlungGauntFactorV1::Constant(1.2),
                    minimum_electron_temperature_k: 1e6,
                    maximum_electron_temperature_k: 20.0 * KELVIN_PER_KEV,
                    maximum_electron_thermal_rest_energy_ratio: 0.05,
                    maximum_ion_to_electron_thermal_variance_ratio: 0.05,
                    maximum_electron_thermal_phase_space_density: 0.01,
                },
            )
            .unwrap(),
        ),
    }
}

fn pure_deuterium(c: FusionNetworkSolve) -> FusionNetworkResult {
    let (plasma, state) = fixture::plasma([0.0, 2e21, 0.0, 0.0, 0.0], 1.0, [20.0, 8.0]);
    advance_fusion_network(
        &state,
        &plasma,
        FusionNetworkChannels::ALL,
        heating(),
        &FusionNetworkTransport::default(),
        Time::new(1.0),
        c,
        || false,
    )
    .unwrap()
}

fn mixed(volume: f64, duration: f64, c: FusionNetworkSolve) -> FusionNetworkResult {
    let (plasma, state) = fixture::plasma(MIXED_DENSITY, volume, [15.0, 4.0]);
    advance_fusion_network(
        &state,
        &plasma,
        FusionNetworkChannels::ALL,
        mixed_heating(),
        &transport(&state),
        Time::new(duration),
        c,
        || false,
    )
    .unwrap()
}

struct Reference {
    amounts: [f64; 5],
    extents: [f64; 4],
    temperature: [f64; 2],
    products: [[f64; 2]; 4],
    // Q, reacting-ion debit, neutrons, charged->electrons, charged->ions,
    // signed electron->ion Coulomb exchange, photons.
    energy: [f64; 7],
}

// Independent direct-SI fixed-step RK4 with compensated increments, platform
// exp/cbrt and separately integrated ion/electron energies and product accounts.
// Bosch-Hale (1992), Eq.14/Table VII: https://doi.org/10.1088/0029-5515/32/4/I07
// The independent reference uses neither production reaction maps nor the
// production thermal/source routines. A: 1024/2048 and B: 16384/32768 refinement
// agree below the assertions' tolerances. These test the declared Maxwellian,
// instantaneous-thermalization closures, not measured plasma accuracy.
const A: Reference = Reference {
    amounts: [
        8.518235228612281e-6,
        0.0032827303761589315,
        5.582018783049605e-6,
        9.187535642987263e-6,
        2.9362164455626753e-6,
    ],
    extents: [
        2.89819486822852e-6,
        9.225557220321419e-6,
        8.480213651278123e-6,
        3.802157733415542e-8,
    ],
    temperature: [253441709.7371331, 193567205.26762873],
    products: [
        [1000597.8780947154, 3943215.5388482716],
        [765830.9372046244, 2233686.9908631616],
        [859366.4371133078, 2521651.962646281],
        [13711.799843990371, 54124.89755506886],
    ],
    energy: [
        11195169.818225566,
        197016.6239438558,
        6176902.529711433,
        4172227.1299663903,
        1043056.7824915973,
        0.0,
        0.0,
    ],
};
const B: Reference = Reference {
    amounts: [
        8.305228965099501e-7,
        1.2432121766416968e-5,
        2.8852102129232982e-6,
        1.6607530641672823e-6,
        4.36121284347422e-7,
    ],
    extents: [
        2.0961832252485085e-8,
        2.386822949106136e-10,
        2.2867762155147293e-10,
        2.4685301475233243e-11,
    ],
    temperature: [116215365.1549411, 116966146.86008859],
    products: [
        [7203.893427029916, 28485.050447162568],
        [19.37252586418339, 57.189742618841535],
        [22.75572602554566, 67.43756930394288],
        [8.851296592914387, 35.049400637422195],
    ],
    energy: [
        35782.38746658981,
        117.21266864552088,
        28542.24018978141,
        5883.685919237785,
        1473.6740262161395,
        -14566.377558932865,
        2552.422396387565,
    ],
};

fn check_reference(result: &FusionNetworkResult, reference: &Reference) {
    let b = result.balance();
    for (i, expected) in result.plasma().indices().into_iter().zip(reference.amounts) {
        close(result.state().amounts()[i].value(), expected, 3e-7, 1e-19);
    }
    for (actual, expected) in b.reaction_extents_mol.into_iter().zip(reference.extents) {
        close(actual, expected, 3e-7, 1e-20);
    }
    for (actual, expected) in [
        result.state().heavy_temperature().value(),
        result.state().electron_temperature().value(),
    ]
    .into_iter()
    .zip(reference.temperature)
    {
        close(actual, expected, 3e-8, 1e-3);
    }
    for (actual, expected) in b
        .product_kinetic_energies_j
        .into_iter()
        .flatten()
        .zip(reference.products.into_iter().flatten())
    {
        close(actual, expected, 3e-7, 1e-7);
    }
    for (actual, expected) in [
        b.nuclear_energy_released_j,
        b.reacting_ion_energy_removed_j,
        b.neutron_energy_escaped_j,
        b.charged_product_electron_energy_j,
        b.charged_product_heavy_energy_j,
        b.electron_to_heavy_coulomb_energy_j,
        b.photon_energy_escaped_j,
    ]
    .into_iter()
    .zip(reference.energy)
    {
        close(actual, expected, 3e-7, 1e-7);
    }
}

fn check_accounts(result: &FusionNetworkResult) {
    let b = result.balance();
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
    let nuclei = [
        NuclearSpeciesV1::Proton,
        NuclearSpeciesV1::Deuterium,
        NuclearSpeciesV1::Tritium,
        NuclearSpeciesV1::Helium3,
        NuclearSpeciesV1::Helium4,
    ];
    let before = result
        .plasma()
        .indices()
        .map(|i| result.initial_state().amounts()[i].value());
    let after = result
        .plasma()
        .indices()
        .map(|i| result.state().amounts()[i].value());
    let mut baryons = [0.0; 2];
    let mut charges = [0.0; 2];
    let mut rest_mass = [0.0; 2];
    for (s, nucleus) in nuclei.into_iter().enumerate() {
        let expected = before[s]
            + FusionReactionV1::ALL
                .into_iter()
                .zip(b.reaction_extents_mol)
                .map(|(r, x)| f64::from(r.stoichiometric_coefficient(nucleus)) * x)
                .sum::<f64>();
        close(after[s], expected, 2e-12, 1e-18);
        for (pool, amount) in [before[s], after[s]].into_iter().enumerate() {
            baryons[pool] += f64::from(nucleus.baryon_number()) * amount;
            charges[pool] += f64::from(nucleus.charge_number()) * amount;
            rest_mass[pool] += nucleus.mass_kg() * (AVOGADRO * amount);
        }
    }
    baryons[1] += b.neutron_amount_mol;
    rest_mass[1] += NuclearSpeciesV1::Neutron.mass_kg() * (AVOGADRO * b.neutron_amount_mol);
    close(baryons[1], baryons[0], 2e-13, 1e-18);
    close(charges[1], charges[0], 2e-13, 1e-18);
    close(
        b.nuclear_energy_released_j,
        (rest_mass[0] - rest_mass[1]) * SPEED_OF_LIGHT_M_S.powi(2),
        3e-9,
        1e-5,
    );
    assert_eq!(
        result.state().amounts()[result.plasma().electron_index()],
        result.initial_state().amounts()[result.plasma().electron_index()]
    );
    assert_eq!(
        result.state().thermodynamics(),
        result.initial_state().thermodynamics()
    );
    let mut products = 0.0;
    let mut charged = [0.0; 2];
    let mut neutrons = 0.0;
    for (r, reaction) in FusionReactionV1::ALL.into_iter().enumerate() {
        for (p, product) in reaction.products().into_iter().enumerate() {
            let energy = b.product_kinetic_energies_j[r][p];
            assert!(energy.is_finite() && energy >= 0.0);
            products += energy;
            if product == NuclearSpeciesV1::Neutron {
                neutrons += energy;
            } else {
                let f = result.heating().product_electron_fractions[r][p];
                charged[0] += f * energy;
                charged[1] += (1.0 - f) * energy;
            }
        }
    }
    close(
        products,
        b.nuclear_energy_released_j + b.reacting_ion_energy_removed_j,
        3e-14,
        1e-7,
    );
    close(neutrons, b.neutron_energy_escaped_j, 3e-14, 1e-7);
    close(charged[0], b.charged_product_electron_energy_j, 3e-14, 1e-7);
    close(charged[1], b.charged_product_heavy_energy_j, 3e-14, 1e-7);
    close(
        b.neutron_amount_mol,
        b.reaction_extents_mol[0] + b.reaction_extents_mol[1],
        3e-14,
        1e-20,
    );
    assert!(b.total_energy_residual_j.abs() < 1e-6);
    assert!(b.heavy_energy_residual_j.abs() < 1e-6);
    assert!(b.electron_energy_residual_j.abs() < 1e-6);
    assert!(b.relative_baryon_residual.abs() < 1e-12);
    assert!(b.relative_charge_residual.abs() < 1e-12);
    assert!(b.maximum_stoichiometry_residual_mol < 1e-16);
    assert!(b.maximum_normalized_local_error <= 1.0);
}

#[test]
fn g1_identical_pair_depletion_consumes_two_ions_per_event_exactly() {
    let (amount, rate, volume, duration) = (2e-4, 1.1e-22, 0.01, 5.0);
    let expected = amount / (1.0 + AVOGADRO * rate * amount * duration / volume);
    let (remaining, extent) = solve::dd_depletion(amount, rate, volume, duration).unwrap();
    close(remaining, expected, 2e-14, 0.0);
    close(remaining + 2.0 * extent, amount, 2e-14, 0.0);
    let (half, x1) = solve::dd_depletion(amount, rate, volume, 0.5 * duration).unwrap();
    let (twice, x2) = solve::dd_depletion(half, rate, volume, 0.5 * duration).unwrap();
    close(twice, remaining, 2e-14, 0.0);
    close(x1 + x2, extent, 2e-14, 0.0);
    assert_eq!(
        solve::dd_depletion(0.0, rate, volume, duration).unwrap(),
        (0.0, 0.0)
    );
    assert_eq!(
        solve::dd_depletion(amount, rate, volume, 0.0).unwrap(),
        (amount, 0.0)
    );
    let (trace, events) = solve::dd_depletion(amount, rate, volume, 1e18).unwrap();
    assert!(trace > 0.0 && trace < 1e-16 * amount);
    assert!(events <= 0.5 * amount);
}

#[test]
fn g1_pure_deuterium_produces_secondary_fuels_and_burns_them_without_seeds() {
    let result = pure_deuterium(controls());
    let [p, _, t, he3, he4] = result.plasma().indices();
    for i in [p, t, he3, he4] {
        assert_eq!(result.initial_state().amounts()[i].value(), 0.0);
        assert!(result.state().amounts()[i].value() > 0.0);
    }
    assert!(
        result
            .balance()
            .reaction_extents_mol
            .iter()
            .all(|&x| x > 0.0)
    );
    check_reference(&result, &A);
    check_accounts(&result);
}

#[test]
fn g1_mixed_burn_couples_each_product_split_to_original_coulomb_and_radiation_sources() {
    let result = mixed(0.025, 0.2, controls());
    check_reference(&result, &B);
    check_accounts(&result);
    assert!(result.balance().electron_to_heavy_coulomb_energy_j < 0.0);
    assert!(result.balance().photon_energy_escaped_j > 0.0);
    assert_eq!(
        result.heating().product_electron_fractions,
        mixed_heating().product_electron_fractions
    );
}

#[test]
fn g1_competing_network_has_second_order_fixed_step_refinement() {
    let mut errors = Vec::new();
    for intervals in [16, 32, 64] {
        let mut c = controls();
        c.relative_tolerance = 0.1;
        c.absolute_amount_tolerance_mol = 1e-3;
        c.absolute_temperature_tolerance_k = 1e6;
        c.energy_tolerance_j = 1e6;
        c.initial_step_s = 1.0 / f64::from(intervals);
        c.minimum_step_s = c.initial_step_s;
        c.maximum_step_s = c.initial_step_s;
        let result = pure_deuterium(c);
        assert_eq!(result.balance().rejected_steps, 0);
        errors.push(
            (result.state().heavy_temperature().value() - A.temperature[0]).abs()
                + (result.state().electron_temperature().value() - A.temperature[1]).abs(),
        );
    }
    for pair in errors.windows(2) {
        let ratio = pair[0] / pair[1];
        assert!(
            (3.4..4.7).contains(&ratio),
            "errors={errors:?}, ratio={ratio}"
        );
    }
}

#[test]
fn g0_channel_selection_preserves_branch_identity_and_checks_absent_fuel_domains() {
    let (plasma, state) = fixture::plasma([0.0, 2e21, 0.0, 0.0, 0.0], 1.0, [20.0, 8.0]);
    let transport = FusionNetworkTransport::default();
    for selected in [1, 2] {
        let mut enabled = [false; 4];
        enabled[selected] = true;
        let result = advance_fusion_network(
            &state,
            &plasma,
            FusionNetworkChannels::new(enabled),
            heating(),
            &transport,
            Time::new(0.1),
            controls(),
            || false,
        )
        .unwrap();
        for (index, extent) in result
            .balance()
            .reaction_extents_mol
            .into_iter()
            .enumerate()
        {
            if index == selected {
                assert!(extent > 0.0);
            } else {
                assert_eq!(extent, 0.0);
            }
        }
        let [p, _, t, he3, he4] = plasma.indices();
        assert_eq!(result.state().amounts()[he4].value(), 0.0);
        if selected == 1 {
            assert!(result.state().amounts()[he3].value() > 0.0);
            assert_eq!(result.state().amounts()[p].value(), 0.0);
            assert_eq!(result.state().amounts()[t].value(), 0.0);
        } else {
            assert!(result.state().amounts()[p].value() > 0.0);
            assert!(result.state().amounts()[t].value() > 0.0);
            assert_eq!(result.state().amounts()[he3].value(), 0.0);
            assert_eq!(result.neutron_amount().value(), 0.0);
        }
        check_accounts(&result);
    }
    let mut low = controls();
    low.minimum_heavy_temperature_k = 0.2 * KELVIN_PER_KEV;
    assert!(matches!(
        advance_fusion_network(
            &state,
            &plasma,
            FusionNetworkChannels::ALL,
            heating(),
            &transport,
            Time::new(0.0),
            low,
            || false
        ),
        Err(FusionNetworkError::Reaction(_))
    ));
    // Disabling D-He3 legitimately restores the other sources' 0.2-keV lower bound.
    advance_fusion_network(
        &state,
        &plasma,
        FusionNetworkChannels::new([true, true, true, false]),
        heating(),
        &transport,
        Time::new(0.0),
        low,
        || false,
    )
    .unwrap();
    let (cold_plasma, cold) = fixture::plasma([0.0, 2e21, 0.0, 0.0, 0.0], 1.0, [0.1, 0.2]);
    low.minimum_heavy_temperature_k = 0.05 * KELVIN_PER_KEV;
    low.maximum_heavy_temperature_k = 0.15 * KELVIN_PER_KEV;
    let mut held = heating();
    held.heavy_power = Power::new(1000.0);
    held.electron_power = Power::new(-200.0);
    let result = advance_fusion_network(
        &cold,
        &cold_plasma,
        FusionNetworkChannels::NONE,
        held,
        &transport,
        Time::new(0.1),
        low,
        || false,
    )
    .unwrap();
    assert_eq!(result.state().amounts(), cold.amounts());
    assert_eq!(result.balance().reaction_extents_mol, [0.0; 4]);
    close(
        result.state().heavy_internal_energy().value() - cold.heavy_internal_energy().value(),
        100.0,
        1e-12,
        1e-8,
    );
    close(
        result.state().electron_internal_energy().value() - cold.electron_internal_energy().value(),
        -20.0,
        1e-12,
        1e-8,
    );
}

#[test]
fn g3_coupled_network_scales_with_volume_and_restarts_on_retained_sources() {
    let whole = mixed(0.025, 0.02, controls());
    let mut larger_controls = controls();
    larger_controls.absolute_amount_tolerance_mol *= 7.0;
    larger_controls.energy_tolerance_j *= 7.0;
    let large = mixed(0.175, 0.02, larger_controls);
    close(
        large.state().heavy_temperature().value(),
        whole.state().heavy_temperature().value(),
        3e-8,
        1e-3,
    );
    close(
        large.state().electron_temperature().value(),
        whole.state().electron_temperature().value(),
        3e-8,
        1e-3,
    );
    close(
        large.state().pressure().value(),
        whole.state().pressure().value(),
        3e-8,
        1e-6,
    );
    for (large_amount, small_amount) in large.state().amounts().iter().zip(whole.state().amounts())
    {
        close(
            large_amount.value(),
            7.0 * small_amount.value(),
            3e-7,
            1e-18,
        );
    }
    for (a, b) in large
        .balance()
        .reaction_extents_mol
        .into_iter()
        .zip(whole.balance().reaction_extents_mol)
    {
        close(a, 7.0 * b, 3e-7, 1e-20);
    }
    for (a, b) in [
        large.balance().neutron_energy_escaped_j,
        large.balance().photon_energy_escaped_j,
        large.balance().electron_to_heavy_coulomb_energy_j,
    ]
    .into_iter()
    .zip([
        whole.balance().neutron_energy_escaped_j,
        whole.balance().photon_energy_escaped_j,
        whole.balance().electron_to_heavy_coulomb_energy_j,
    ]) {
        close(a, 7.0 * b, 3e-7, 1e-7);
    }
    let first = mixed(0.025, 0.01, controls());
    let second = advance_fusion_network(
        first.state(),
        first.plasma(),
        first.channels(),
        first.heating(),
        first.transport(),
        Time::new(0.01),
        first.controls(),
        || false,
    )
    .unwrap();
    close(
        second.state().heavy_temperature().value(),
        whole.state().heavy_temperature().value(),
        3e-8,
        1e-3,
    );
    close(
        second.state().electron_temperature().value(),
        whole.state().electron_temperature().value(),
        3e-8,
        1e-3,
    );
    for (part, total) in second.state().amounts().iter().zip(whole.state().amounts()) {
        close(part.value(), total.value(), 3e-7, 1e-18);
    }
    for i in 0..4 {
        close(
            first.balance().reaction_extents_mol[i] + second.balance().reaction_extents_mol[i],
            whole.balance().reaction_extents_mol[i],
            3e-7,
            1e-20,
        );
    }
    close(
        first.balance().photon_energy_escaped_j + second.balance().photon_energy_escaped_j,
        whole.balance().photon_energy_escaped_j,
        3e-7,
        1e-7,
    );
    assert_eq!(second.plasma(), first.plasma());
    assert_eq!(second.transport(), first.transport());
    check_accounts(&second);
}

#[test]
fn g4_network_cancellation_work_and_source_domain_refuse_without_partial_state() {
    let (plasma, state) = fixture::plasma(MIXED_DENSITY, 0.025, [15.0, 4.0]);
    let initial = state.clone();
    let transport = transport(&state);
    let cancelled = advance_fusion_network(
        &state,
        &plasma,
        FusionNetworkChannels::ALL,
        mixed_heating(),
        &transport,
        Time::new(0.2),
        controls(),
        || true,
    );
    assert_eq!(cancelled.unwrap_err(), FusionNetworkError::Cancelled);
    let mut polls = 0;
    let cancelled = advance_fusion_network(
        &state,
        &plasma,
        FusionNetworkChannels::ALL,
        mixed_heating(),
        &transport,
        Time::new(0.2),
        controls(),
        || {
            polls += 1;
            polls >= 500
        },
    );
    assert_eq!(cancelled.unwrap_err(), FusionNetworkError::Cancelled);
    let mut small_budget = controls();
    small_budget.maximum_steps = 1;
    assert_eq!(
        advance_fusion_network(
            &state,
            &plasma,
            FusionNetworkChannels::ALL,
            mixed_heating(),
            &transport,
            Time::new(0.2),
            small_budget,
            || false
        )
        .unwrap_err(),
        FusionNetworkError::StepBudget
    );
    let mut invalid_heating = mixed_heating();
    invalid_heating.product_electron_fractions[0][1] = 0.5;
    assert!(matches!(
        advance_fusion_network(
            &state,
            &plasma,
            FusionNetworkChannels::ALL,
            invalid_heating,
            &transport,
            Time::new(0.2),
            controls(),
            || false
        ),
        Err(FusionNetworkError::Invalid(_))
    ));
    let mut narrow = controls();
    narrow.maximum_heavy_temperature_k = 14.0 * KELVIN_PER_KEV;
    assert!(
        advance_fusion_network(
            &state,
            &plasma,
            FusionNetworkChannels::ALL,
            mixed_heating(),
            &transport,
            Time::new(0.2),
            narrow,
            || false
        )
        .is_err()
    );
    let (limited_plasma, limited) =
        fixture::plasma_with_helium3_limit([0.0, 2e21, 0.0, 0.0, 0.0], 1.0, [20.0, 8.0], 1e7);
    assert!(matches!(
        advance_fusion_network(
            &limited,
            &limited_plasma,
            FusionNetworkChannels::ALL,
            heating(),
            &FusionNetworkTransport::default(),
            Time::new(0.2),
            controls(),
            || false
        ),
        Err(FusionNetworkError::Plasma(_))
    ));
    assert_eq!(state, initial);
}

#[test]
fn g4_network_stiff_differential_coulomb_heating_requires_resolution() {
    let (plasma, initial) = fixture::plasma([0.0, 0.0, 0.0, 0.0, 5e19], 1.0, [5.0; 2]);
    let original = initial.clone();
    let mut source = transport(&initial);
    source.bremsstrahlung = None;
    let ch = initial.heavy_heat_capacity_j_k();
    let ce = initial.electron_heat_capacity_j_k();
    let conductance = source
        .coulomb
        .as_ref()
        .unwrap()
        .evaluate_at_state(&initial, || false)
        .unwrap()
        .electron_heavy_thermal_conductance_w_m3_k()
        * initial.volume().value();
    let rate = conductance / ch + conductance / ce;
    let duration = 100.0 / rate;
    let power = (ch * ce / (ch + ce)) * (0.001 * initial.heavy_temperature().value()) * rate;
    let mut heat = heating();
    heat.heavy_power = Power::new(power);
    heat.electron_power = Power::new(-power);
    let mut c = controls();
    c.initial_step_s = duration;
    c.minimum_step_s = duration;
    c.maximum_step_s = duration;
    c.maximum_steps = 1;
    assert_eq!(
        advance_fusion_network(
            &initial,
            &plasma,
            FusionNetworkChannels::NONE,
            heat,
            &source,
            Time::new(duration),
            c,
            || false,
        ),
        Err(FusionNetworkError::StepBudget),
    );
    assert_eq!(initial, original);
}
