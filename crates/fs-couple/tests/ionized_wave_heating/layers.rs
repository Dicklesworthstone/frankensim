//! Finite coherent Maxwell films, source-bound gas radiation and carrier birth.
use super::*;
use fs_couple::ionized_wave_heating::IonizedGasReactingHeatControls;
use fs_couple::ionized_wave_heating::layers::{
    IonizedGasLayerHeatControls, IonizedGasLayerInput, LayeredIonizedGasHeating,
    LayeredIonizedGasHeatingError as Error,
};
use fs_em::plane_wave::{
    Polarization,
    layers::{DielectricLayer, IonizedGasLayerExcitation, PreparedIonizedGasLayerIllumination},
    plasma::IonizedGasWaveResponse,
};
use fs_gas::chemistry::{
    collisional_radiative::CollisionalRadiativeSolve,
    diffuse_recombination::DiffuseRecombinationSolve,
    ionization_network::{
        IonizationNetworkSolve,
        transport::{NetworkRadiationCell, diffuse::advance_network_diffuse_recombination},
    },
};
use fs_material::ionized_gas::{IonizedGasMomentumSource, IonizedGasTransportBank};
use fs_qty::{Area, Energy};
use fs_thermochem::electron_collision::{
    ElectronCollisionKindV1, ElectronCrossSectionPointV1, ElectronCrossSectionV1,
};
use fs_thermochem::equilibrium::composition::FrozenIonizedGasStateV1;
use fs_thermochem::kinetics::{
    ElectronImpactRateV1,
    atomic_ionization::{
        AtomicIonizationV1, AtomicStatisticalWeightsV1,
        photoionization::{AtomicPhotoionizationV1, PhotoionizationCrossSectionPointV1},
        radiative_recombination::AtomicRadiativeRecombinationV1,
    },
};
use fs_thermochem::plasma_radiation::{
    BremsstrahlungControlsV1, BremsstrahlungGauntFactorV1, ElectronIonBremsstrahlungV1,
};
use fs_thermochem::{ReactionId, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R};

fn bank(state: &FrozenIonizedGasStateV1, momentum_scale: f64) -> IonizedGasTransportBank {
    let (template, original_card) = fixture::collision_source_and_card(
        "H+",
        &[
            ("T_h", state.heavy_temperature().value()),
            ("T_e", state.electron_temperature().value()),
        ],
        state.pressure().value(),
        1e9,
    );
    IonizedGasTransportBank::new(
        state,
        vec![
            IonizedGasMomentumSource::Material {
                template,
                original_card,
            },
            IonizedGasMomentumSource::ElectronNeutral {
                cross_sections: vec![fixture::momentum_cross_section("H", momentum_scale * 2e-20)],
                controls: super::controls().rate_controls,
                maximum_target_thermal_speed_ratio: 0.1,
            },
        ],
        None,
        fixture::limits(),
        || false,
    )
    .unwrap()
}

fn atomic(state: &FrozenIonizedGasStateV1) -> AtomicRadiativeRecombinationV1 {
    let gap = 13.598_434_599_702 * fixture::E;
    let provenance = Provenance {
        source: "manufactured original H impact ramp and constant bound-free film fixture".into(),
        license: "analytic test only".into(),
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
        state.thermodynamics().clone(),
        ReactionId::new("coherent-film-H").unwrap(),
        fixture::id("H+"),
        ElectronImpactRateV1::new(
            fixture::id("e-"),
            impact,
            super::controls().rate_controls,
            0.1,
        )
        .unwrap(),
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

struct Parcel {
    state: FrozenIonizedGasStateV1,
    bank: IonizedGasTransportBank,
    background: ResolvedMaterialStatePoint,
    card: MaterialCard,
    atomic: [AtomicRadiativeRecombinationV1; 1],
    cooling: ElectronIonBremsstrahlungV1,
}
fn parcel(fraction: f64) -> Parcel {
    let state = fixture::frozen_state(fraction, 8000.0, 8000.0, 1.0);
    let bank = bank(&state, 1.0);
    // The same original background can represent exact neutrality. Its actual
    // pressure must be queried separately, not inherited from a charged gas.
    let populated = bank
        .evaluate(&fixture::frozen_state(0.001, 8000.0, 8000.0, 1.0), || false)
        .unwrap()
        .unwrap();
    let (template, card) = super::background(&populated);
    let query = template
        .query()
        .with("p", state.pressure().value())
        .unwrap();
    let background = template.at_query_point(&card, &query).unwrap();
    let atomic = [atomic(&state)];
    let cooling = ElectronIonBremsstrahlungV1::new(
        state.thermodynamics().clone(),
        BremsstrahlungControlsV1 {
            frequency_integrated_gaunt_factor: BremsstrahlungGauntFactorV1::Constant(1.2),
            minimum_electron_temperature_k: 2000.0,
            maximum_electron_temperature_k: 20000.0,
            maximum_electron_thermal_rest_energy_ratio: 0.01,
            maximum_ion_to_electron_thermal_variance_ratio: 0.01,
            maximum_electron_thermal_phase_space_density: 0.01,
        },
    )
    .unwrap();
    Parcel {
        state,
        bank,
        background,
        card,
        atomic,
        cooling,
    }
}

fn solve_controls() -> IonizedGasLayerHeatControls {
    let mut local =
        CollisionalRadiativeSolve::new(Temperature::new(2000.0), Temperature::new(20000.0));
    local.relative_tolerance = 1e-5;
    local.absolute_amount_tolerance_mol = 1e-15;
    local.temperature_tolerance_k = 1e-5;
    local.energy_tolerance_j = 1e-6;
    local.initial_step_s = 1e-9;
    local.maximum_step_s = 1e-9;
    IonizedGasLayerHeatControls::new(
        IonizedGasReactingHeatControls {
            network: IonizationNetworkSolve::new(local),
            minimum_heavy_temperature_k: 2000.0,
            maximum_heavy_temperature_k: 20000.0,
            maximum_heavy_caloric_iterations: 96,
            maximum_temperature_increment_k: 200.0,
            maximum_absolute_population_increment_mol: 1e-6,
            maximum_relative_population_increment: 0.1,
        },
        DiffuseRecombinationSolve {
            relative_tolerance: 1e-5,
            ..Default::default()
        },
    )
}

fn with_column(
    neutral_receiver: bool,
    excited: bool,
    free_free: bool,
    test: impl FnOnce(&LayeredIonizedGasHeating<'_, '_, '_>),
) {
    let parcels = [
        parcel(0.001),
        parcel(if neutral_receiver { 0.0 } else { 0.002 }),
    ];
    let transports: Vec<_> = parcels
        .iter()
        .map(|p| p.bank.evaluate(&p.state, || false).unwrap())
        .collect();
    let harmonics: Vec<_> = transports
        .iter()
        .map(|t| t.as_ref().map(|t| t.prepare_harmonic(|| false).unwrap()))
        .collect();
    let omega = 1e13;
    let films: Vec<_> = harmonics
        .iter()
        .zip(&parcels)
        .map(|(harmonic, p)| DielectricLayer {
            response: match harmonic {
                Some(h) => DielectricResponse::IonizedGas(
                    IonizedGasWaveResponse::try_new(h, &p.background, omega, || false).unwrap(),
                ),
                None => DielectricResponse::Lossless(&p.background),
            },
            thickness_m: 0.1,
        })
        .collect();
    let amplitude = if excited { 1e4 } else { 0.0 };
    let wave = PreparedIonizedGasLayerIllumination::try_new(
        DielectricResponse::Lossless(&parcels[0].background),
        DielectricResponse::Lossless(&parcels[1].background),
        omega,
        &films,
        IonizedGasLayerExcitation {
            incidence_cosine: 0.8,
            polarization: Polarization::Te,
            incoming_tangential_e_v_m: [
                C64::new(amplitude, 0.2 * amplitude),
                C64::new(-0.3 * amplitude, 0.4 * amplitude),
            ],
        },
        || false,
    )
    .unwrap();
    let inputs: Vec<_> = parcels
        .iter()
        .map(|p| {
            let input = IonizedGasLayerInput::new(&p.state, &p.atomic, &p.bank, &p.card);
            if free_free {
                input.with_optically_thin_bremsstrahlung(&p.cooling)
            } else {
                input
            }
        })
        .collect();
    let heating =
        LayeredIonizedGasHeating::try_new(&wave, &inputs, Area::new(10.0), || false).unwrap();
    test(&heating);
}

#[test]
fn coherent_film_work_and_both_radiation_losses_close_before_next_maxwell_solve() {
    with_column(false, true, true, |heating| {
        let dt = Time::new(1e-9);
        let step = heating.advance(dt, solve_controls(), || false).unwrap();
        let mut storage = 0.0;
        for (i, cell) in step.chemistry().cells().iter().enumerate() {
            let old = cell.initial_state();
            let next = cell.state();
            let powers = heating.wave().layer_carrier_heating(i).unwrap();
            let b = cell.balance();
            near(
                b.electron_deposited_energy_j,
                powers.electron_w_m2 * heating.area().value() * dt.value(),
                5e-13,
            );
            near(
                b.heavy_deposited_energy_j,
                powers.heavy_w_m2 * heating.area().value() * dt.value(),
                5e-13,
            );
            assert_eq!(next.amounts()[1], next.amounts()[2]);
            near(
                next.amounts()[0].value() + next.amounts()[1].value(),
                1.0,
                2e-14,
            );
            near(next.mass().value(), old.mass().value(), 2e-14);
            storage += (next.amounts()[1].value() - old.amounts()[1].value())
                * 13.598_434_599_702
                * fixture::E
                * fixture::NA
                + 1.5
                    * R
                    * (next.amounts()[2].value() * next.electron_temperature().value()
                        - old.amounts()[2].value() * old.electron_temperature().value())
                + 1.5 * R * (next.heavy_temperature().value() - old.heavy_temperature().value());
            let refreshed = &step.refreshed_layers()[i];
            assert_eq!(refreshed.transport().unwrap().state(), next);
            assert_eq!(
                refreshed.background().card_identity(),
                heating.inputs()[i].original_background_card.content_hash()
            );
            assert_eq!(
                refreshed.background().state_coordinate("p"),
                Some(next.pressure().value())
            );
            assert!(b.free_free_escaped_energy_j > 0.0);
            assert_eq!(
                cell.free_free_cooling_source(),
                heating.inputs()[i].optically_thin_bremsstrahlung
            );
        }
        let r = step.chemistry().diffuse_radiation().unwrap().balance();
        assert!(r.escaped_energy_j > 0.0);
        assert!(step.absorbed_energy_j() > 0.0);
        assert!(
            (storage + r.escaped_energy_j + r.free_free_escaped_energy_j
                - step.absorbed_energy_j())
            .abs()
                < 1e-6
        );
        assert!(step.energy_residual_j().abs() < 1e-6);
        let harmonics: Vec<_> = step
            .refreshed_layers()
            .iter()
            .map(|p| p.transport().unwrap().prepare_harmonic(|| false).unwrap())
            .collect();
        let films: Vec<_> = harmonics
            .iter()
            .zip(step.refreshed_layers())
            .map(|(h, p)| DielectricLayer {
                response: DielectricResponse::IonizedGas(
                    IonizedGasWaveResponse::try_new(
                        h,
                        p.background(),
                        heating.wave().omega_rad_s(),
                        || false,
                    )
                    .unwrap(),
                ),
                thickness_m: 0.1,
            })
            .collect();
        let next_wave = PreparedIonizedGasLayerIllumination::try_new(
            heating.wave().incident_response(),
            heating.wave().exit_response(),
            heating.wave().omega_rad_s(),
            &films,
            heating.wave().excitation(),
            || false,
        )
        .unwrap();
        assert_ne!(
            next_wave.layer_absorbed_power_w_m2(0),
            heating.wave().layer_absorbed_power_w_m2(0)
        );
    });
}

#[test]
fn neutral_receiver_is_born_from_diffuse_photons_and_dark_drive_matches_native_column() {
    with_column(true, false, false, |heating| {
        let dt = Time::new(1e-9);
        let step = heating.advance(dt, solve_controls(), || false).unwrap();
        let cells = step.chemistry().cells();
        assert_eq!(cells[1].initial_state().amounts()[2].value(), 0.0);
        assert!(cells[1].state().amounts()[2].value() > 0.0);
        assert!(
            cells[1]
                .diffuse_radiation()
                .unwrap()
                .nonlocal_absorbed_photons
                > 0.0
        );
        assert!(step.refreshed_layers()[1].transport().is_some());
        assert_eq!(step.absorbed_energy_j(), 0.0);
        let original: Vec<_> = heating
            .inputs()
            .iter()
            .zip(cells)
            .map(|(input, cell)| {
                NetworkRadiationCell::new(input.initial, input.sources)
                    .with_thermal(cell.thermal_inputs().unwrap())
            })
            .collect();
        let native = advance_network_diffuse_recombination(
            &original,
            heating.area(),
            &[],
            dt,
            step.chemistry().diffuse_radiation().unwrap().controls(),
            || false,
        )
        .unwrap();
        assert_eq!(&native, step.chemistry());
        for cell in cells {
            assert_eq!(cell.balance().electron_deposited_energy_j, 0.0);
            assert_eq!(cell.balance().heavy_deposited_energy_j, 0.0);
        }
    });
}

#[test]
fn independent_source_geometry_and_interval_refusals_leave_original_column_reusable() {
    let neutral = parcel(0.0);
    let threshold = neutral.atomic[0].source().threshold().value();
    let film = DielectricLayer {
        response: DielectricResponse::Lossless(&neutral.background),
        thickness_m: 0.1,
    };
    let ionizing_drive = PreparedIonizedGasLayerIllumination::try_new(
        film.response,
        film.response,
        1.01 * threshold * core::f64::consts::TAU / 6.626_070_15e-34,
        &[film],
        IonizedGasLayerExcitation {
            incidence_cosine: 1.0,
            polarization: Polarization::Te,
            incoming_tangential_e_v_m: [C64::ONE, C64::ZERO],
        },
        || false,
    )
    .unwrap();
    assert!(matches!(
        LayeredIonizedGasHeating::try_new(
            &ionizing_drive,
            &[IonizedGasLayerInput::new(
                &neutral.state,
                &neutral.atomic,
                &neutral.bank,
                &neutral.card
            )],
            Area::new(10.0),
            || false
        ),
        Err(Error::Invalid {
            layer: Some(0),
            reason: "coherent drive must lie below every original photoionization threshold"
        })
    ));
    with_column(false, true, false, |heating| {
        let inputs = heating.inputs();
        assert!(matches!(
            LayeredIonizedGasHeating::try_new(heating.wave(), inputs, Area::new(11.0), || false),
            Err(Error::Invalid { layer: Some(0), .. })
        ));
        let changed = fixture::frozen_state(0.003, 8000.0, 8000.0, 1.0);
        let mut wrong = inputs.to_vec();
        wrong[0].initial = &changed;
        assert!(matches!(
            LayeredIonizedGasHeating::try_new(heating.wave(), &wrong, heating.area(), || false),
            Err(Error::Invalid { layer: Some(0), .. })
        ));
        let changed_bank = bank(inputs[0].initial, 1.1);
        wrong = inputs.to_vec();
        wrong[0].transport_bank = &changed_bank;
        assert!(matches!(
            LayeredIonizedGasHeating::try_new(heating.wave(), &wrong, heating.area(), || false),
            Err(Error::Invalid { layer: Some(0), .. })
        ));
        let wrong_card = inputs[0]
            .transport_bank
            .sources()
            .find_map(|(_, law)| match law {
                IonizedGasMomentumSource::Material { original_card, .. } => Some(original_card),
                _ => None,
            })
            .unwrap();
        wrong = inputs.to_vec();
        wrong[0].original_background_card = wrong_card;
        assert!(matches!(
            LayeredIonizedGasHeating::try_new(heating.wave(), &wrong, heating.area(), || false),
            Err(Error::Invalid { layer: Some(0), .. })
        ));
        let dt = Time::new(1e-9);
        assert!(matches!(
            heating.advance(dt, solve_controls(), || true),
            Err(Error::Cancelled)
        ));
        let mut limits = solve_controls();
        limits.reacting.maximum_temperature_increment_k = 1e-30;
        assert!(matches!(
            heating.advance(dt, limits, || false),
            Err(Error::TemperatureIncrement { .. })
        ));
        limits = solve_controls();
        limits.reacting.maximum_absolute_population_increment_mol = 0.0;
        limits.reacting.maximum_relative_population_increment = 1e-30;
        assert!(matches!(
            heating.advance(dt, limits, || false),
            Err(Error::PopulationIncrement { .. })
        ));
        let accepted = heating.advance(dt, solve_controls(), || false).unwrap();
        for (cell, input) in accepted.chemistry().cells().iter().zip(inputs) {
            assert_eq!(cell.initial_state(), input.initial);
        }
        let replay = heating.advance(dt, solve_controls(), || false).unwrap();
        assert_eq!(accepted.chemistry(), replay.chemistry());
    });
}
