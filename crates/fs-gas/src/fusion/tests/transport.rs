use super::*;
use fs_thermochem::coulomb::{CoulombLimitsV1, CoulombScreeningV1, CoulombShortRangeCutoffV1};
use fs_thermochem::plasma_radiation::{BremsstrahlungControlsV1, BremsstrahlungGauntFactorV1};

mod confinement;

fn thermal_sources(state: &FrozenIonizedGasStateV1) -> DtFusionTransport {
    DtFusionTransport {
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

fn coupled(
    density: [f64; 3],
    volume: f64,
    temperatures: [f64; 2],
    fraction: f64,
    duration: f64,
    c: DtFusionSolve,
) -> DtFusionResult {
    let (model, state) = source::plasma(density, volume, temperatures);
    advance_dt_fusion_with_transport(
        &state,
        &model,
        heating(fraction),
        &thermal_sources(&state),
        Time::new(duration),
        c,
        || false,
    )
    .unwrap()
}

#[test]
fn g1_free_free_cooling_matches_exact_monatomic_solution_without_fuel() {
    let (model, state) = source::plasma([0.0, 0.0, 5e19], 1.0, [10.0, 5.0]);
    let mut transport = thermal_sources(&state);
    transport.coulomb = None;
    let coefficient = transport
        .bremsstrahlung
        .as_ref()
        .unwrap()
        .evaluate(&state, || false)
        .unwrap()
        .total_power_over_sqrt_temperature_w_k_sqrt();
    let time = 0.3;
    let ce = state.electron_heat_capacity_j_k();
    let root = state.electron_temperature().value().sqrt() - coefficient * time / (2.0 * ce);
    let expected = root * root;
    let result = advance_dt_fusion_with_transport(
        &state,
        &model,
        heating(0.8),
        &transport,
        Time::new(time),
        controls(),
        || false,
    )
    .unwrap();
    close(
        result.state().electron_temperature().value(),
        expected,
        1e-12,
        1e-4,
    );
    close(
        result.state().heavy_temperature().value(),
        state.heavy_temperature().value(),
        1e-14,
        1e-5,
    );
    close(
        result.balance().photon_energy_escaped_j,
        ce * (state.electron_temperature().value() - expected),
        1e-12,
        1e-6,
    );
    assert_eq!(result.state().amounts(), state.amounts());
    assert_eq!(result.balance().reaction_extent_mol, 0.0);
    assert_eq!(result.balance().electron_to_heavy_coulomb_energy_j, 0.0);
    assert!(result.balance().total_energy_residual_j.abs() < 1e-6);
}

#[test]
fn g0_coulomb_exchange_relaxes_both_pools_without_creating_energy() {
    let (model, state) = source::plasma([0.0, 5e19, 0.0], 1.0, [10.0, 2.0]);
    let mut transport = thermal_sources(&state);
    transport.bremsstrahlung = None;
    let result = advance_dt_fusion_with_transport(
        &state,
        &model,
        heating(0.8),
        &transport,
        Time::new(0.2),
        controls(),
        || false,
    )
    .unwrap();
    let th = result.state().heavy_temperature().value();
    let te = result.state().electron_temperature().value();
    let equilibrium = 6.0 * KELVIN_PER_KEV;
    assert!(th < state.heavy_temperature().value() && th >= equilibrium);
    assert!(te > state.electron_temperature().value() && te <= equilibrium);
    assert!(result.balance().electron_to_heavy_coulomb_energy_j < 0.0);
    close(
        result.state().heavy_internal_energy().value() - state.heavy_internal_energy().value(),
        result.balance().electron_to_heavy_coulomb_energy_j,
        1e-12,
        1e-7,
    );
    assert_eq!(result.balance().reaction_extent_mol, 0.0);
    assert_eq!(result.balance().photon_energy_escaped_j, 0.0);
    assert!(result.balance().total_energy_residual_j.abs() < 1e-6);
}

#[test]
fn g1_coupled_fuel_heat_exchange_and_photons_match_independent_rk4_oracles() {
    // Independently integrated pool energies and outgoing accounts, using
    // Python math and direct SI pair formulas. 16/32/64 steps show fourth-order
    // RK4 convergence; 1024/2048 differences are below 5e-5 K and 4e-8 J.
    // The references test the declared Gaunt/cutoff/thermalization closure.
    let cases = [
        (
            [5e19, 5e19, 0.0],
            1.0,
            [10.0, 5.0],
            0.8,
            0.1,
            [
                3.732168289935527e-8,
                99949359.62306,
                79544401.33905877,
                50706.82653053553,
                12818.793815298794,
                186.8271760512673,
                -35757.474454975105,
                1441.409810658896,
            ],
        ),
        (
            [3e20, 7e20, 0.0],
            0.01,
            [10.0, 3.0],
            0.65,
            0.02,
            [
                4.2529643908190474e-9,
                83373224.37124175,
                73160703.5958682,
                5777.16385124651,
                1459.9810032886762,
                19.42017060370427,
                -7262.271117375493,
                269.66578367335904,
            ],
        ),
    ];
    for (density, volume, temperatures, fraction, duration, expected) in cases {
        let result = coupled(
            density,
            volume,
            temperatures,
            fraction,
            duration,
            controls(),
        );
        let b = result.balance();
        let actual = [
            b.reaction_extent_mol,
            result.state().heavy_temperature().value(),
            result.state().electron_temperature().value(),
            b.neutron_energy_escaped_j,
            b.alpha_energy_deposited_j,
            b.reacting_ion_energy_removed_j,
            b.electron_to_heavy_coulomb_energy_j,
            b.photon_energy_escaped_j,
        ];
        for (index, (a, e)) in actual.into_iter().zip(expected).enumerate() {
            close(
                a,
                e,
                if index == 1 || index == 2 { 3e-8 } else { 3e-7 },
                if index == 0 { 1e-20 } else { 1e-7 },
            );
        }
        assert!(b.total_energy_residual_j.abs() < 1e-6);
        assert!(b.electron_energy_residual_j.abs() < 1e-6);
        assert!(b.heavy_energy_residual_j.abs() < 1e-6);
        assert!(b.relative_baryon_residual.abs() < 1e-12);
        assert!(b.relative_charge_residual.abs() < 1e-12);
        assert!(b.maximum_normalized_local_error <= 1.0);
    }
}

#[test]
fn g1_coupled_midpoint_has_second_order_fixed_step_convergence() {
    let reference = [99949359.62306, 79544401.33905877];
    let mut errors = Vec::new();
    for intervals in [4, 8, 16] {
        let mut c = controls();
        c.relative_tolerance = 0.1;
        c.absolute_amount_tolerance_mol = 1e-6;
        c.absolute_temperature_tolerance_k = 1e6;
        c.energy_tolerance_j = 1e3;
        c.initial_step_s = 0.1 / f64::from(intervals);
        c.minimum_step_s = c.initial_step_s;
        c.maximum_step_s = c.initial_step_s;
        let result = coupled([5e19, 5e19, 0.0], 1.0, [10.0, 5.0], 0.8, 0.1, c);
        assert_eq!(result.balance().rejected_steps, 0);
        errors.push(
            (result.state().heavy_temperature().value() - reference[0]).abs()
                + (result.state().electron_temperature().value() - reference[1]).abs(),
        );
    }
    for pair in errors.windows(2) {
        let ratio = pair[0] / pair[1];
        assert!(
            (3.0..5.5).contains(&ratio),
            "errors={errors:?}, ratio={ratio}"
        );
    }
}

#[test]
fn g1_temperature_and_charge_dependent_gaunt_fit_matches_independent_burn() {
    let (model, state) = source::plasma([5e19, 5e19, 0.0], 1.0, [10.0, 5.0]);
    let mut transport = thermal_sources(&state);
    let mut radiation = transport.bremsstrahlung.as_ref().unwrap().controls();
    radiation.frequency_integrated_gaunt_factor = BremsstrahlungGauntFactorV1::ItohNonrelativistic;
    transport.bremsstrahlung =
        Some(ElectronIonBremsstrahlungV1::new(state.thermodynamics().clone(), radiation).unwrap());
    let result = advance_dt_fusion_with_transport(
        &state,
        &model,
        heating(0.8),
        &transport,
        Time::new(0.1),
        controls(),
        || false,
    )
    .unwrap();
    // Independent direct-power polynomial evaluation for each ion charge and
    // 2048-step RK4. Halving from 1024 steps changes T by <6e-8 K and En by <1e-9 J.
    let b = result.balance();
    let actual = [
        b.reaction_extent_mol,
        result.state().heavy_temperature().value(),
        result.state().electron_temperature().value(),
        b.neutron_energy_escaped_j,
        b.alpha_energy_deposited_j,
        b.reacting_ion_energy_removed_j,
        b.electron_to_heavy_coulomb_energy_j,
        b.photon_energy_escaped_j,
    ];
    let expected = [
        3.732378518218533e-8,
        99955214.72411978,
        79564004.12157919,
        50709.684203134755,
        12819.5169105779,
        186.84015043764742,
        -35745.48568658078,
        1389.4026753618346,
    ];
    for (index, (a, e)) in actual.into_iter().zip(expected).enumerate() {
        close(
            a,
            e,
            if index == 1 || index == 2 { 3e-8 } else { 3e-7 },
            if index == 0 { 1e-20 } else { 1e-7 },
        );
    }
    assert!(b.total_energy_residual_j.abs() < 1e-6);
    assert!(b.photon_energy_escaped_j < 1441.409810658896);

    // Hydrogen's published gamma-squared domain ends near 13.6 keV. The fitted
    // operator refuses hotter electrons instead of silently changing its law.
    let (model, too_hot) = source::plasma([5e19, 5e19, 0.0], 1.0, [10.0, 15.0]);
    assert!(matches!(
        advance_dt_fusion_with_transport(
            &too_hot,
            &model,
            heating(0.8),
            &transport,
            Time::new(0.1),
            controls(),
            || false
        ),
        Err(DtFusionError::Radiation(BremsstrahlungErrorV1::Domain(_)))
    ));
}

#[test]
fn g3_coupled_volume_scaling_and_restart_preserve_physical_accounts() {
    let a = coupled([5e19, 5e19, 0.0], 1.0, [10.0, 5.0], 0.8, 0.1, controls());
    let b = coupled([5e19, 5e19, 0.0], 0.01, [10.0, 5.0], 0.8, 0.1, controls());
    close(
        a.state().heavy_temperature().value(),
        b.state().heavy_temperature().value(),
        1e-8,
        1e-3,
    );
    close(
        a.state().electron_temperature().value(),
        b.state().electron_temperature().value(),
        1e-8,
        1e-3,
    );
    close(
        a.balance().photon_energy_escaped_j,
        100.0 * b.balance().photon_energy_escaped_j,
        1e-8,
        1e-6,
    );
    let first = coupled([5e19, 5e19, 0.0], 1.0, [10.0, 5.0], 0.8, 0.04, controls());
    // Retain the same original operator models; refresh their populations from
    // the accepted changed state, including previously absent helium ash.
    let second = advance_dt_fusion_with_transport(
        first.state(),
        first.plasma(),
        first.heating(),
        first.transport(),
        Time::new(0.06),
        controls(),
        || false,
    )
    .unwrap();
    close(
        a.state().heavy_temperature().value(),
        second.state().heavy_temperature().value(),
        1e-8,
        1e-3,
    );
    close(
        a.state().electron_temperature().value(),
        second.state().electron_temperature().value(),
        1e-8,
        1e-3,
    );
    close(
        a.balance().reaction_extent_mol,
        first.balance().reaction_extent_mol + second.balance().reaction_extent_mol,
        1e-7,
        1e-18,
    );
    close(
        a.balance().photon_energy_escaped_j,
        first.balance().photon_energy_escaped_j + second.balance().photon_energy_escaped_j,
        1e-8,
        1e-6,
    );
}

#[test]
fn g4_coupled_cancellation_and_source_domain_refusal_publish_no_state() {
    let (model, state) = source::plasma([5e19, 5e19, 0.0], 1.0, [10.0, 5.0]);
    let before = state.clone();
    let transport = thermal_sources(&state);
    let mut polls = 0;
    assert!(matches!(
        advance_dt_fusion_with_transport(
            &state,
            &model,
            heating(0.8),
            &transport,
            Time::new(0.1),
            controls(),
            || {
                polls += 1;
                polls == 50
            }
        ),
        Err(DtFusionError::Cancelled)
    ));
    assert_eq!(state, before);
    let (model, cold) = source::plasma([5e19, 5e19, 0.0], 1.0, [10.0, 0.05]);
    assert!(matches!(
        advance_dt_fusion_with_transport(
            &cold,
            &model,
            heating(0.8),
            &transport,
            Time::new(0.1),
            controls(),
            || false
        ),
        Err(DtFusionError::Radiation(BremsstrahlungErrorV1::Domain(_)))
    ));
}
