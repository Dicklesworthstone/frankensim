use super::*;
use fs_qty::Temperature;
use fs_thermochem::coulomb::{CoulombLimitsV1, CoulombModelV1, CoulombShortRangeCutoffV1};
use fs_thermochem::equilibrium::composition::IdealGasEquilibriumV1;
use fs_thermochem::fusion::{KELVIN_PER_KEV, NuclearSpeciesV1};

#[path = "../../tests/fixture.rs"]
mod fixture;

const KEV_J: f64 = 1.602_176_634e-16;

fn limits() -> DtFusionPlasmaLimitsV1 {
    DtFusionPlasmaLimitsV1 {
        maximum_coulomb_coupling: 0.01,
        maximum_thermal_phase_space_density: 0.01,
        maximum_thermal_energy_rest_energy_ratio: 0.05,
    }
}

fn close(actual: f64, expected: f64, relative: f64, absolute: f64) {
    assert!(
        (actual - expected).abs() <= absolute + relative * expected.abs(),
        "actual={actual:.17e}, expected={expected:.17e}"
    );
}

fn transport(
    reference: &FrozenIonizedGasStateV1,
    projectile: usize,
    energy_rest_limit: f64,
) -> KineticAlphaTransport {
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
    KineticAlphaTransport::new(
        TestParticleCoulombV1::new(background, projectile, energy_rest_limit).unwrap(),
        0.05,
        0.25,
        1e-4,
        100_000,
    )
    .unwrap()
}

fn grid() -> AlphaEnergyGrid {
    AlphaEnergyGrid::new(
        [0.05, 1000.0, 3500.0, 5000.0]
            .map(|e| Energy::new(e * KEV_J))
            .into(),
        || false,
    )
    .unwrap()
}

fn seeded(
    reference: &FrozenIonizedGasStateV1,
    plasma: &LightIonFusionPlasmaV1,
) -> KineticFusionState {
    let fast = 0.5 * reference.amounts()[plasma.indices()[4]].value();
    let lower = 0.25 * fast;
    let populations = [0.0, lower, fast - lower, 0.0].map(Amount::new);
    KineticFusionState::new(reference, plasma, &grid(), &populations, limits(), || false).unwrap()
}

fn scenario() -> (KineticFusionState, KineticAlphaTransport) {
    let (plasma, reference) = fixture::plasma([0.0, 5e19, 5e19, 0.0, 2e16], 1.0, [8.0, 2.0]);
    let alpha = transport(&reference, plasma.indices()[4], 0.01);
    let state = seeded(&reference, &plasma);
    (state, alpha)
}

fn stores(state: &KineticFusionState) -> [f64; 3] {
    [
        state.thermal_heavy_energy().value(),
        state.thermal_electron_energy().value(),
        state.alpha_energy().value(),
    ]
}

fn works(result: &KineticFusionVolumeResult) -> [f64; 3] {
    [
        result.boundary_work_on_heavy().value(),
        result.boundary_work_on_electrons().value(),
        result.boundary_work_on_alphas().value(),
    ]
}

#[test]
fn g1_kinetic_volume_motion_scales_three_stores_and_pressure_without_thermalization() {
    let (initial, alpha) = scenario();
    // An already elapsed hybrid state remains at the same physical time.
    let initial = KineticFusionState::from_parts(
        &initial.reference,
        initial.plasma(),
        initial.energy_grid(),
        initial.alpha_amounts,
        initial.thermal_helium4_amount,
        initial.alpha_amount().value(),
        Time::new(3.25),
        limits(),
        &mut || false,
    )
    .unwrap();
    let unchanged = initial.clone();
    for ratio in [0.125_f64, 8.0] {
        let result = advance_adiabatic_kinetic_fusion_volume(
            &initial,
            &alpha,
            Volume::new(ratio),
            limits(),
            || false,
        )
        .unwrap();
        let final_state = result.state();
        let factor = ratio.powf(-2.0 / 3.0);
        for ((before, after), work) in stores(&initial)
            .into_iter()
            .zip(stores(final_state))
            .zip(works(&result))
        {
            close(after, factor * before, 4e-14, 0.0);
            close(work, (factor - 1.0) * before, 4e-14, 0.0);
            close(after - before, work, 4e-14, 0.0);
        }
        close(
            final_state.pressure().value(),
            initial.pressure().value() * factor / ratio,
            4e-14,
            0.0,
        );
        for (before, after) in initial
            .energy_grid()
            .energies()
            .iter()
            .zip(final_state.energy_grid().energies())
        {
            close(after.value(), factor * before.value(), 4e-14, 0.0);
        }
        for (before, after) in initial
            .energy_grid()
            .widths_j()
            .iter()
            .zip(final_state.energy_grid().widths_j())
        {
            close(*after, factor * before, 4e-14, 0.0);
        }
        assert_eq!(final_state.elapsed_time(), Time::new(3.25));
        assert_eq!(
            final_state.alpha_populations_mol(),
            initial.alpha_populations_mol()
        );
        assert_eq!(final_state.alpha_amount(), initial.alpha_amount());
        assert_eq!(final_state.thermal_amounts(), initial.thermal_amounts());
        assert_eq!(final_state.amounts(), initial.amounts());
        assert_eq!(final_state.plasma(), initial.plasma());
        close(
            final_state.internal_energy().value() - initial.internal_energy().value(),
            result.total_boundary_work().value(),
            4e-14,
            0.0,
        );
        close(
            result.energy_residual_j(),
            0.0,
            0.0,
            3e-13 * initial.internal_energy().value(),
        );
    }
    assert_eq!(initial, unchanged);
}

#[test]
fn g0_kinetic_volume_reversal_noop_and_zero_alpha_match_existing_thermal_motion() {
    let (initial, alpha) = scenario();
    let forward = advance_adiabatic_kinetic_fusion_volume(
        &initial,
        &alpha,
        Volume::new(0.37),
        limits(),
        || false,
    )
    .unwrap();
    let reverse = advance_adiabatic_kinetic_fusion_volume(
        forward.state(),
        &alpha,
        initial.volume(),
        limits(),
        || false,
    )
    .unwrap();
    assert_eq!(forward.initial_state(), &initial);
    assert_eq!(
        reverse.state().alpha_populations_mol(),
        initial.alpha_populations_mol()
    );
    assert_eq!(reverse.state().amounts(), initial.amounts());
    for (actual, expected) in stores(reverse.state()).into_iter().zip(stores(&initial)) {
        close(actual, expected, 4e-14, 0.0);
    }
    for (actual, expected) in reverse
        .state()
        .energy_grid()
        .energies()
        .iter()
        .zip(initial.energy_grid().energies())
    {
        close(actual.value(), expected.value(), 4e-14, 0.0);
    }
    close(
        forward.total_boundary_work().value() + reverse.total_boundary_work().value(),
        0.0,
        0.0,
        3e-13 * initial.internal_energy().value(),
    );
    let noop = advance_adiabatic_kinetic_fusion_volume(
        &initial,
        &alpha,
        initial.volume(),
        limits(),
        || false,
    )
    .unwrap();
    assert_eq!(noop.state(), &initial);
    assert_eq!(
        noop.state().energy_grid().energies().as_ptr(),
        initial.energy_grid().energies().as_ptr()
    );
    assert_eq!(works(&noop), [0.0; 3]);
    assert_eq!(noop.energy_residual_j(), 0.0);

    let thermal = KineticFusionState::from_thermal(
        &initial.reference,
        initial.plasma(),
        initial.energy_grid(),
        limits(),
        || false,
    )
    .unwrap();
    let hybrid_result = advance_adiabatic_kinetic_fusion_volume(
        &thermal,
        &alpha,
        Volume::new(1.7),
        limits(),
        || false,
    )
    .unwrap();
    let thermal_result = advance_adiabatic_fusion_volume(
        &initial.reference,
        initial.plasma().into(),
        Volume::new(1.7),
        limits(),
        || false,
    )
    .unwrap();
    assert_eq!(hybrid_result.state().alpha_energy().value(), 0.0);
    assert_eq!(hybrid_result.boundary_work_on_alphas().value(), 0.0);
    close(
        hybrid_result.state().pressure().value(),
        thermal_result.state().pressure().value(),
        4e-14,
        0.0,
    );
    close(
        hybrid_result.total_boundary_work().value(),
        thermal_result.total_boundary_work().value(),
        4e-14,
        0.0,
    );
    close(
        hybrid_result.state().internal_energy().value(),
        thermal_result.state().internal_energy().value(),
        4e-14,
        0.0,
    );
}

#[test]
fn g1_adjacent_volumes_retain_nonzero_work_on_kinetic_alphas() {
    let (initial, alpha) = scenario();
    for bits in [1.0_f64.to_bits() - 1, 1.0_f64.to_bits() + 1] {
        let volume = f64::from_bits(bits);
        let result = advance_adiabatic_kinetic_fusion_volume(
            &initial,
            &alpha,
            Volume::new(volume),
            limits(),
            || false,
        )
        .unwrap();
        let fraction = ((-2.0 / 3.0) * (volume - 1.0).ln_1p()).exp_m1();
        for (before, work) in stores(&initial).into_iter().zip(works(&result)) {
            assert_ne!(work, 0.0);
            close(work, before * fraction, 4e-15, 0.0);
            assert_eq!(work.is_sign_positive(), volume < 1.0);
        }
        assert_eq!(
            result.state().alpha_populations_mol(),
            initial.alpha_populations_mol()
        );
    }
}

#[test]
fn g4_kinetic_volume_source_particle_domains_and_late_cancellation_are_atomic() {
    let (initial, alpha) = scenario();
    let unchanged = initial.clone();
    let (other_plasma, other_reference) =
        fixture::plasma_with_helium3_limit([0.0, 5e19, 5e19, 0.0, 2e16], 1.0, [8.0, 2.0], 1e9);
    let substituted = transport(&other_reference, other_plasma.indices()[4], 0.01);
    let wrong_projectile = transport(&initial.reference, initial.plasma().indices()[1], 0.01);
    for source in [&substituted, &wrong_projectile] {
        assert!(matches!(
            advance_adiabatic_kinetic_fusion_volume(
                &initial,
                source,
                initial.volume(),
                limits(),
                || false
            ),
            Err(FusionNetworkError::Invalid(_))
        ));
    }
    let limited_energy = transport(&initial.reference, initial.plasma().indices()[4], 0.002);
    assert_eq!(
        advance_adiabatic_kinetic_fusion_volume(
            &initial,
            &limited_energy,
            Volume::new(0.125),
            limits(),
            || false
        )
        .unwrap_err(),
        FusionNetworkError::OutsideDomain
    );
    for (number, energy) in [(1e-5, 0.25), (0.05, 0.001)] {
        let limited =
            KineticAlphaTransport::new(alpha.collisions().clone(), number, energy, 1e-4, 100_000)
                .unwrap();
        assert_eq!(
            advance_adiabatic_kinetic_fusion_volume(
                &initial,
                &limited,
                initial.volume(),
                limits(),
                || false
            )
            .unwrap_err(),
            FusionNetworkError::OutsideDomain
        );
    }
    let mut edge_populations = [Amount::new(0.0); 4];
    edge_populations[3] = initial.alpha_amount();
    let edge = KineticFusionState::new(
        &initial.reference,
        initial.plasma(),
        initial.energy_grid(),
        &edge_populations,
        limits(),
        || false,
    )
    .unwrap();
    assert_eq!(
        advance_adiabatic_kinetic_fusion_volume(&edge, &alpha, edge.volume(), limits(), || false)
            .unwrap_err(),
        FusionNetworkError::OutsideDomain
    );
    assert!(matches!(
        advance_adiabatic_kinetic_fusion_volume(
            &initial,
            &alpha,
            Volume::new(0.0),
            limits(),
            || false
        ),
        Err(FusionNetworkError::Invalid(_))
    ));
    let mut calls = 0;
    advance_adiabatic_kinetic_fusion_volume(&initial, &alpha, Volume::new(0.5), limits(), || {
        calls += 1;
        false
    })
    .unwrap();
    let mut cancelled_calls = 0;
    assert_eq!(
        advance_adiabatic_kinetic_fusion_volume(
            &initial,
            &alpha,
            Volume::new(0.5),
            limits(),
            || {
                cancelled_calls += 1;
                cancelled_calls == calls
            }
        ),
        Err(FusionNetworkError::Cancelled)
    );
    assert_eq!(initial, unchanged);
}

#[test]
fn g3_kinetic_volume_motion_rejects_gaps_in_absent_product_sources() {
    use fs_evidence::ValidityDomain;
    use fs_thermochem::{Nasa9RegionV1, Nasa9StandardStateModelV1};

    let (original, _) = scenario();
    let bank = original.plasma().thermodynamics();
    let ids = original
        .plasma()
        .indices()
        .map(|i| bank.models()[i].species().clone());
    let mut models = bank.models().to_vec();
    let absent_index = original.plasma().indices()[3];
    assert_eq!(original.amounts()[absent_index].value(), 0.0);
    let source = &models[absent_index];
    let regions =
        [(1e4, 9.0 * KELVIN_PER_KEV), (11.0 * KELVIN_PER_KEV, 2e9)].map(|(lower, upper)| {
            let mut card = source.regions()[0].card().clone();
            card.validity = ValidityDomain::unconstrained().with("T", lower, upper);
            Nasa9RegionV1::from_card(card).unwrap()
        });
    models[absent_index] = Nasa9StandardStateModelV1::new(
        source.species().clone(),
        source.molar_mass(),
        source.convention().clone(),
        regions.into(),
    )
    .unwrap();
    let plasma = LightIonFusionPlasmaV1::new(
        IdealGasEquilibriumV1::new(
            bank.elemental_matrix().clone(),
            bank.charge_vector().clone(),
            models,
        )
        .unwrap(),
        ids[0].clone(),
        ids[1].clone(),
        ids[2].clone(),
        ids[3].clone(),
        ids[4].clone(),
    )
    .unwrap();
    let reference = FrozenIonizedGasStateV1::new(
        plasma.thermodynamics(),
        original.amounts(),
        original.volume(),
        original.background_heavy_temperature(),
        original.electron_temperature(),
        64.0 * f64::EPSILON,
        || false,
    )
    .unwrap();
    let initial = seeded(&reference, &plasma);
    let alpha = transport(&reference, plasma.indices()[4], 0.01);
    // Both endpoint temperatures are supported, but the absent He3 source
    // has no law along part of the monotone compression path.
    let endpoint = FrozenIonizedGasStateV1::new(
        plasma.thermodynamics(),
        original.amounts(),
        Volume::new(0.125),
        Temperature::new(32.0 * KELVIN_PER_KEV),
        Temperature::new(8.0 * KELVIN_PER_KEV),
        64.0 * f64::EPSILON,
        || false,
    )
    .unwrap();
    plasma
        .admit_state_with_limits(&endpoint, limits(), || false)
        .unwrap();
    assert!(matches!(
        advance_adiabatic_kinetic_fusion_volume(&initial, &alpha, endpoint.volume(), limits(), || false),
        Err(FusionNetworkError::Plasma(DtFusionPlasmaErrorV1::Invalid(message))) if message.contains("gap")
    ));
}
