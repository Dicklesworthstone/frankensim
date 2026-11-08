//! G0/G1/G3/G4: native reaction energy and source refresh after Maxwell work.
use super::*;
use fs_couple::ionized_wave_heating::IonizedGasReactingHeatControls;
use fs_gas::chemistry::collisional_radiative::CollisionalRadiativeSolve;
use fs_gas::chemistry::ionization_network::IonizationNetworkSolve;
use fs_material::ionized_gas::{
    CarrierMomentumSource, IonizedGasMomentumSource, IonizedGasTransportBank,
};
use fs_qty::{Area, Energy};
use fs_thermochem::electron_collision::{
    ElectronCollisionKindV1, ElectronCrossSectionPointV1, ElectronCrossSectionV1,
};
use fs_thermochem::kinetics::ElectronImpactRateV1;
use fs_thermochem::kinetics::atomic_ionization::photoionization::{
    AtomicPhotoionizationV1, PhotoionizationCrossSectionPointV1,
};
use fs_thermochem::kinetics::atomic_ionization::radiative_recombination::AtomicRadiativeRecombinationV1;
use fs_thermochem::kinetics::atomic_ionization::{AtomicIonizationV1, AtomicStatisticalWeightsV1};
use fs_thermochem::{ReactionId, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R};

fn source(heating: &MagnetizedIonizedGasHeating<'_, '_>) -> AtomicRadiativeRecombinationV1 {
    let bank = heating.source().state().thermodynamics();
    let gap = bank.models()[1].regions()[0].coefficients()[7] * R / fixture::NA;
    let provenance = Provenance {
        source: "manufactured ramp-impact and constant bound-free cross sections on original H levels; not measured data".into(),
        license: "analytic test fixture".into(),
        artifact: None,
    };
    let impact = ElectronCrossSectionV1::tabulated(
        fixture::id("H"),
        ElectronCollisionKindV1::Ionization,
        provenance.clone(),
        vec![
            ElectronCrossSectionPointV1 {
                energy: Energy::new(gap),
                cross_section: Area::new(0.0),
            },
            ElectronCrossSectionPointV1 {
                energy: Energy::new(gap + 1000.0 * fixture::E),
                cross_section: Area::new(1e-18),
            },
        ],
        Area::new(1e-18),
    )
    .unwrap();
    let atomic = AtomicIonizationV1::new(
        bank.clone(),
        ReactionId::new("wave-heated-H-ionization").unwrap(),
        fixture::id("H+"),
        ElectronImpactRateV1::new(fixture::id("e-"), impact, controls().rate_controls, 0.1).unwrap(),
        AtomicStatisticalWeightsV1 { neutral: 2, ion: 1 },
        1e-12,
    )
    .unwrap();
    let photon = AtomicPhotoionizationV1::new(
        atomic,
        provenance.clone(),
        vec![
            PhotoionizationCrossSectionPointV1 {
                photon_energy: Energy::new(gap),
                cross_section: Area::new(1e-24),
            },
            PhotoionizationCrossSectionPointV1 {
                photon_energy: Energy::new(gap + 1000.0 * fixture::E),
                cross_section: Area::new(1e-24),
            },
        ],
    )
    .unwrap();
    AtomicRadiativeRecombinationV1::new(photon, Area::new(1e-24), provenance).unwrap()
}

fn transport_bank(
    heating: &MagnetizedIonizedGasHeating<'_, '_>,
    ion_card: &MaterialCard,
    electron_cross_section_scale: f64,
) -> IonizedGasTransportBank {
    let original = heating.source();
    let mut sources = Vec::new();
    for carrier in original.carriers() {
        let source = if carrier.is_electron() {
            let CarrierMomentumSource::ElectronNeutral {
                contributions,
                maximum_target_thermal_speed_ratio,
                ..
            } = carrier.relaxation().source()
            else {
                unreachable!()
            };
            IonizedGasMomentumSource::ElectronNeutral {
                cross_sections: vec![fixture::momentum_cross_section(
                    "H",
                    2e-20 * electron_cross_section_scale,
                )],
                controls: contributions[0].evaluation().controls(),
                maximum_target_thermal_speed_ratio: *maximum_target_thermal_speed_ratio,
            }
        } else {
            IonizedGasMomentumSource::Material {
                template: carrier.relaxation().clone(),
                original_card: ion_card.clone(),
            }
        };
        sources.push(source);
    }
    IonizedGasTransportBank::new(
        original.state(),
        sources,
        original.coulomb_model().cloned(),
        original.limits(),
        || false,
    )
    .unwrap()
}

fn reaction_controls() -> IonizedGasReactingHeatControls {
    let mut local =
        CollisionalRadiativeSolve::new(Temperature::new(2000.0), Temperature::new(20000.0));
    local.relative_tolerance = 1e-6;
    local.absolute_amount_tolerance_mol = 1e-16;
    local.temperature_tolerance_k = 1e-7;
    local.energy_tolerance_j = 1e-7;
    local.initial_step_s = 1e-9;
    local.maximum_step_s = 1e-7;
    IonizedGasReactingHeatControls {
        network: IonizationNetworkSolve::new(local),
        minimum_heavy_temperature_k: 2000.0,
        maximum_heavy_temperature_k: 20000.0,
        maximum_heavy_caloric_iterations: 96,
        maximum_temperature_increment_k: 200.0,
        maximum_absolute_population_increment_mol: 0.0,
        maximum_relative_population_increment: 0.1,
    }
}

#[test]
fn g1_wave_work_changes_atomic_storage_and_next_maxwell_response() {
    with_heating(
        false,
        8000.0,
        8000.0,
        true,
        |heating, ion_card, background_card| {
            let bank = transport_bank(heating, ion_card, 1.0);
            let duration = Time::new(1e-6);
            let step = heating
                .advance_reacting(
                    duration,
                    &[source(heating)],
                    &bank,
                    reaction_controls(),
                    background_card,
                    || false,
                )
                .unwrap();
            let result = step.chemistry();
            let next = result.state();
            let initial = heating.source().state();
            let balance = result.balance();
            assert_eq!(result.initial_state(), initial);
            assert_eq!(next.thermodynamics(), initial.thermodynamics());
            assert_eq!(next.volume(), initial.volume());
            assert_eq!(result.elapsed_time(), duration);
            assert_eq!(result.photon_field(), &[]);
            assert_eq!(balance.absorbed_radiation_energy_j, 0.0);
            assert!(balance.emitted_radiation_energy_j > 0.0);
            assert!(balance.chemical_storage_change_j > 0.0);
            assert!(next.amounts()[2].value() > initial.amounts()[2].value());
            assert_eq!(next.amounts()[1], next.amounts()[2]);
            near(next.amounts()[0].value() + next.amounts()[1].value(), 1.0, 2e-14);
            near(next.mass().value(), initial.mass().value(), 2e-14);
            near(
                balance.electron_deposited_energy_j,
                heating.absorption().electron_power_w() * duration.value(),
                3e-13,
            );
            near(
                balance.heavy_deposited_energy_j,
                heating.absorption().heavy_power_w() * duration.value(),
                3e-13,
            );
            // Independently reconstruct original H-level and monatomic storage.
            let extent = next.amounts()[1].value() - initial.amounts()[1].value();
            let chemical = extent * 13.598_434_599_702 * fixture::E * fixture::NA;
            let electron = 1.5
                * R
                * (next.amounts()[2].value() * next.electron_temperature().value()
                    - initial.amounts()[2].value() * initial.electron_temperature().value());
            let heavy =
                1.5 * R * (next.heavy_temperature().value() - initial.heavy_temperature().value());
            assert!((chemical - balance.chemical_storage_change_j).abs() < 2e-9);
            assert!((electron - balance.electron_sensible_energy_change_j).abs() < 2e-9);
            assert!((heavy - balance.heavy_sensible_energy_change_j).abs() < 2e-8);
            assert!(
                (chemical + electron + heavy + balance.emitted_radiation_energy_j
                    - step.absorbed_energy_j())
                .abs()
                    < 3e-7
            );
            let frozen = heating
                .advance(
                    duration,
                    controls(),
                    200.0,
                    &[ion_card],
                    background_card,
                    || false,
                )
                .unwrap();
            assert!(
                next.electron_temperature().value()
                    < frozen.energy().state().electron_temperature().value()
            );
            let refreshed = step.refreshed_transport().unwrap();
            assert_eq!(refreshed.state(), next);
            assert_eq!(step.transport_bank(), &bank);
            assert_ne!(refreshed.conductivity_s_m(), heating.source().conductivity_s_m());
            assert_eq!(
                step.refreshed_background().card_identity(),
                background_card.content_hash()
            );
            assert_eq!(
                step.refreshed_background().state_coordinate("p"),
                Some(next.pressure().value())
            );
            let previous = heating.absorption().wave();
            let harmonic = refreshed.prepare_harmonic(|| false).unwrap();
            let transport = harmonic
                .prepare_magnetized(
                    previous.interface().transport().omega_rad_s(),
                    previous.interface().transport().magnetic_field_t(),
                    magnetic_limits(),
                    || false,
                )
                .unwrap();
            let interface = MagnetizedIonizedGasInterface::try_new(
                previous.interface().incident_response(),
                &transport,
                step.refreshed_background(),
                previous.interface().outgoing_normal(),
                || false,
            )
            .unwrap();
            let next_wave = interface
                .prepare_wave(previous.incident_electric_v_m(), || false)
                .unwrap();
            assert_ne!(
                next_wave.transmitted_power_w_m2(),
                previous.transmitted_power_w_m2()
            );
            assert!(next_wave.power_balance_relative_residual() < 1e-10);
        },
    );
}

#[test]
fn g3_reacting_coulomb_exchange_retains_the_signed_wave_work_once() {
    with_heating(
        true,
        6000.0,
        12000.0,
        true,
        |heating, ion_card, background_card| {
            let bank = transport_bank(heating, ion_card, 1.0);
            let duration = Time::new(1e-10);
            let step = heating
                .advance_reacting(
                    duration,
                    &[source(heating)],
                    &bank,
                    reaction_controls(),
                    background_card,
                    || false,
                )
                .unwrap();
            let balance = step.chemistry().balance();
            let thermal = step.chemistry().thermal_inputs().unwrap();
            assert_eq!(
                thermal.collisions.as_ref().unwrap().coulomb.as_ref(),
                heating.source().coulomb_model()
            );
            assert_eq!(
                thermal.collisions.as_ref().unwrap().momentum_sources,
                heating.momentum_sources()
            );
            assert_ne!(heating.absorption().exchange_correction_to_electrons_w(), 0.0);
            near(
                balance.electron_deposited_energy_j,
                duration.value()
                    * (heating.absorption().electron_friction_power_w()
                        + heating.absorption().exchange_correction_to_electrons_w()),
                3e-13,
            );
            near(
                balance.heavy_deposited_energy_j,
                duration.value()
                    * (heating.absorption().heavy_friction_power_w()
                        - heating.absorption().exchange_correction_to_electrons_w()),
                3e-13,
            );
            near(
                balance.electron_deposited_energy_j + balance.heavy_deposited_energy_j,
                step.absorbed_energy_j(),
                3e-11,
            );
            assert!(balance.electron_to_heavy_coulomb_energy_j > 0.0);
            assert!(balance.electron_to_heavy_neutral_energy_j > 0.0);
            near(
                balance.electron_to_heavy_energy_j,
                balance.electron_to_heavy_neutral_energy_j
                    + balance.electron_to_heavy_coulomb_energy_j,
                2e-13,
            );
            assert!(balance.energy_residual_j.abs() < 1e-7);
            assert!(
                step.refreshed_transport()
                    .unwrap()
                    .coulomb_model()
                    .unwrap()
                    .is_bound_to_populations(step.chemistry().state())
            );
        },
    );
}

#[test]
fn g3_zero_wave_allows_chemical_evolution_without_adding_field_energy() {
    with_heating(
        false,
        8000.0,
        8000.0,
        false,
        |heating, ion_card, background_card| {
            let bank = transport_bank(heating, ion_card, 1.0);
            let step = heating
                .advance_reacting(
                    Time::new(1e-6),
                    &[source(heating)],
                    &bank,
                    reaction_controls(),
                    background_card,
                    || false,
                )
                .unwrap();
            let balance = step.chemistry().balance();
            assert_eq!(step.absorbed_energy_j(), 0.0);
            assert_eq!(balance.electron_deposited_energy_j, 0.0);
            assert_eq!(balance.heavy_deposited_energy_j, 0.0);
            assert_ne!(step.chemistry().state().amounts(), heating.source().state().amounts());
            assert!(
                (balance.chemical_storage_change_j
                    + balance.electron_sensible_energy_change_j
                    + balance.heavy_sensible_energy_change_j
                    + balance.emitted_radiation_energy_j)
                .abs()
                    < 1e-7
            );
        },
    );
}

#[test]
fn g4_reacting_increment_source_and_cancellation_refusals_are_atomic() {
    with_heating(
        false,
        8000.0,
        8000.0,
        true,
        |heating, ion_card, background_card| {
            let bank = transport_bank(heating, ion_card, 1.0);
            let retained_bank = bank.clone();
            let retained_source = heating.source().clone();
            let retained_power = heating.absorption().absorbed_power_w();
            let sources = [source(heating)];
            let duration = Time::new(1e-6);
            let mut restricted = reaction_controls();
            restricted.maximum_relative_population_increment = 1e-12;
            assert!(matches!(
                heating.advance_reacting(
                    duration,
                    &sources,
                    &bank,
                    restricted,
                    background_card,
                    || false
                ),
                Err(IonizedWaveHeatingError::PopulationIncrement { .. })
            ));
            restricted = reaction_controls();
            restricted.maximum_temperature_increment_k = 1e-12;
            assert!(matches!(
                heating.advance_reacting(
                    duration,
                    &sources,
                    &bank,
                    restricted,
                    background_card,
                    || false
                ),
                Err(IonizedWaveHeatingError::TemperatureIncrement { .. })
            ));
            let altered = transport_bank(heating, ion_card, 2.0);
            assert!(matches!(
                heating.advance_reacting(
                    duration,
                    &sources,
                    &altered,
                    reaction_controls(),
                    background_card,
                    || false
                ),
                Err(IonizedWaveHeatingError::Invalid(_))
            ));
            assert!(matches!(
                heating.advance_reacting(
                    duration,
                    &sources,
                    &bank,
                    reaction_controls(),
                    ion_card,
                    || false
                ),
                Err(IonizedWaveHeatingError::Invalid(_))
            ));
            assert!(matches!(
                heating.advance_reacting(
                    duration,
                    &sources,
                    &bank,
                    reaction_controls(),
                    background_card,
                    || true
                ),
                Err(IonizedWaveHeatingError::Cancelled)
            ));
            let mut polls = 0;
            assert!(matches!(
                heating.advance_reacting(
                    duration,
                    &sources,
                    &bank,
                    reaction_controls(),
                    background_card,
                    || {
                        polls += 1;
                        polls > 100
                    }
                ),
                Err(IonizedWaveHeatingError::Cancelled)
            ));
            assert_eq!(bank, retained_bank);
            assert_eq!(heating.source(), &retained_source);
            assert_eq!(heating.absorption().absorbed_power_w(), retained_power);
            assert!(
                heating
                    .advance_reacting(
                        duration,
                        &sources,
                        &bank,
                        reaction_controls(),
                        background_card,
                        || false
                    )
                    .is_ok()
            );
        },
    );
}
