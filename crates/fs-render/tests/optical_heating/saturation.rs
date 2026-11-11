//! G1/G3/G4: finite emitter inventories, actual pump heat and accepted afterglow.

use super::*;
use fs_blake3::ContentHash;
use fs_render::volumes::heating::fluorescence::FluorescenceSource;
use fs_render::volumes::heating::luminescence::saturation::{
    PrescribedPhotonPump, SaturableEmitterCell, SaturableLuminescentHeating,
};
use fs_render::volumes::heating::luminescence::{
    LuminescenceSource, SingletTripletLimits, SingletTripletSource,
};
use fs_render::volumes::transport::luminescence::LuminescenceGrid;
use fs_thermal::enthalpy::DepositionObserver;
use fs_thermal::luminescence::ArrheniusQuenching;

const HC: f64 = 6.626_070_15e-34 * 299_792_458.0 * 1e9;
const COUNT: f64 = 1e19;

#[path = "saturation/beam.rs"]
mod beam;

#[path = "saturation/beams.rs"]
mod beams;

#[path = "saturation/annihilation.rs"]
mod annihilation;

#[path = "saturation/recycling.rs"]
mod recycling;

include!("saturation_environment.rs");

fn source(level_nm: f64, lifetime_s: f64, weights: [f64; 2], id: u8) -> LuminescenceSource<2> {
    LuminescenceSource::try_new(
        FluorescenceSource::try_new(
            ContentHash([id; 32]),
            (380.0, 500.0),
            (290.0, 800.0),
            1.0,
            [650.0, 800.0],
            weights,
        )
        .unwrap(),
        level_nm,
        lifetime_s,
    )
    .unwrap()
}

fn setup(
    cx: &Cx<'_>,
    rates: [f64; 2],
    populations: [f64; 3],
    crossing: f64,
) -> (
    DepositionSimulation<'static>,
    SaturableLuminescentHeating<2>,
) {
    let density = [1.0];
    let medium = grid(cx, &density, transparent_host());
    let thermal = thermal(vec![cell(1, 1.0, 1.0, 0.0)], &SENSIBLE, vec![]);
    let pump = PrescribedPhotonPump::try_new(
        ContentHash([231; 32]),
        [400.0, 450.0],
        [1e-20; 2],
        [rates[0] * 1e20, rates[1] * 1e20],
    )
    .unwrap();
    let singlet = source(500.0, 1.0, [1.0, 0.0], 232);
    let triplet = SingletTripletSource::try_new(
        source(600.0, 5.0, [0.0, 1.0], 233),
        crossing,
        ArrheniusQuenching::try_new(0.0, 0.0).unwrap(),
        SingletTripletLimits::default(),
    )
    .unwrap();
    let heating = SaturableLuminescentHeating::try_new(
        &medium,
        singlet,
        triplet,
        &[SaturableEmitterCell {
            grid_cell: 0,
            thermal_cell: CellId(1),
            initial_populations: populations,
            pump,
        }],
        &thermal,
        1,
        1e-7,
        None,
        || Ok(()),
    )
    .unwrap();
    (thermal, heating)
}

fn relative(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-10 * actual.abs().max(expected.abs()).max(1e-30),
        "{actual:.17e} != {expected:.17e}"
    );
}

fn balance(heating: &SaturableLuminescentHeating<2>) {
    let e = heating.energy();
    relative(
        e.initial_excitation_j + heating.absorbed_pump_energy_j().iter().sum::<f64>(),
        e.stored_excitation_j
            + e.escaped_j.iter().sum::<f64>()
            + e.unresolved_radiation_j.iter().sum::<f64>()
            + heating.thermal_input_j(),
    );
    assert!(e.residual_j.abs() < 1e-7);
    assert!(e.absolute_roundoff_j < 1e-7);
    relative(heating.populations(0).unwrap().iter().sum(), COUNT);
}

#[test]
fn g1_saturation_limits_population_and_debits_actual_absorbed_photons() {
    with_cx(|cx| {
        let (thermal, heating) = setup(cx, [9.0, 0.0], [COUNT, 0.0, 0.0], 0.0);
        let (next, accepted, audit) = thermal
            .advance_observed(&heating, 2.0, config(2.0), || Ok(()))
            .unwrap();
        let excited = COUNT * 0.9 * (1.0 - (-20.0_f64).exp());
        let integrated_singlets = COUNT * 0.9 * (2.0 - (1.0 - (-20.0_f64).exp()) / 10.0);
        let pump_events = 9.0 * (COUNT * 2.0 - integrated_singlets);
        let counts = accepted.populations(0).unwrap();
        relative(counts[1], excited);
        relative(counts[0], COUNT - excited);
        assert_eq!(counts[2], 0.0);
        relative(
            accepted.absorbed_pump_energy_j()[0],
            pump_events * HC / 400.0,
        );
        relative(
            accepted.energy().escaped_j[0],
            integrated_singlets * HC / 650.0,
        );
        assert!(accepted.absorbed_pump_energy_j()[0] < 0.2 * 9.0 * COUNT * 2.0 * HC / 400.0);
        relative(next.energy_j(), accepted.thermal_input_j());
        relative(audit.external_heat_j, accepted.thermal_input_j());
        assert_eq!(audit.substeps, 1);
        balance(&accepted);
    });
}

#[test]
fn g1_long_lived_triplets_shelve_emitters_and_reduce_the_pump_rate() {
    with_cx(|cx| {
        let (thermal, heating) = setup(cx, [10.0, 0.0], [COUNT, 0.0, 0.0], 3.0);
        let (_, accepted, _) = thermal
            .advance_observed(&heating, 100.0, config(100.0), || Ok(()))
            .unwrap();
        // Steady state: 4 S = 10 G, T = 3 S / 0.2, G+S+T=N.
        let singlet = COUNT / 16.4;
        relative(accepted.populations(0).unwrap()[1], singlet);
        relative(accepted.populations(0).unwrap()[2], 15.0 * singlet);
        assert!(accepted.populations(0).unwrap()[0] < 0.03 * COUNT);
        balance(&accepted);
    });
}

#[test]
fn g3_pump_edges_split_the_native_clock_and_afterglow_reaches_the_original_observer() {
    with_cx(|cx| {
        let (thermal, heating) = setup(cx, [4.0, 0.0], [COUNT, 0.0, 0.0], 3.0);
        let heating = heating.with_pump_window(0.25, 1.0).unwrap();
        let (_, whole, audit) = thermal
            .advance_observed(&heating, 3.0, config(3.0), || Ok(()))
            .unwrap();
        let (middle, charged, _) = thermal
            .advance_observed(&heating, 1.0, config(3.0), || Ok(()))
            .unwrap();
        let (_, split, _) = middle
            .advance_observed(&charged, 3.0, config(3.0), || Ok(()))
            .unwrap();
        assert_eq!(audit.substeps, 3);
        assert_eq!(whole.populations(0), split.populations(0));
        assert_eq!(whole.energy(), split.energy());
        assert_eq!(
            whole.absorbed_pump_energy_j(),
            charged.absorbed_pump_energy_j()
        );
        let density = [1.0];
        let medium = grid(cx, &density, transparent_host());
        let view = LuminescenceGrid::try_new(
            cx,
            &medium,
            whole.emission_source_identity(),
            [650.0, 800.0],
            whole.cell_emission_w(),
            1,
        )
        .unwrap();
        let seen = view
            .transport(
                cx,
                Point3::new(-1.0, 0.5, 0.5),
                Vec3::new(1.0, 0.0, 0.0),
                f64::INFINITY,
            )
            .unwrap();
        for lane in 0..2 {
            relative(
                seen.emitted_radiance_w_m2_sr[lane],
                whole.cell_emission_w()[0].1[lane] / (4.0 * core::f64::consts::PI),
            );
        }
        assert!(whole.last_mean_cell_emission_w()[0].1[0] > whole.cell_emission_w()[0].1[0]);
        balance(&whole);
    });
}

#[test]
fn g3_precharged_emitters_decay_without_inventing_pump_energy() {
    with_cx(|cx| {
        let (thermal, heating) =
            setup(cx, [0.0; 2], [0.25 * COUNT, 0.25 * COUNT, 0.5 * COUNT], 0.0);
        let (_, accepted, _) = thermal
            .advance_observed(&heating, 5.0, config(5.0), || Ok(()))
            .unwrap();
        relative(
            accepted.energy().initial_excitation_j,
            COUNT * HC * (0.25 / 500.0 + 0.5 / 600.0),
        );
        assert_eq!(accepted.absorbed_pump_energy_j(), [0.0; 2]);
        relative(
            accepted.populations(0).unwrap()[1],
            0.25 * COUNT * (-5.0_f64).exp(),
        );
        relative(
            accepted.populations(0).unwrap()[2],
            0.5 * COUNT * (-1.0_f64).exp(),
        );
        balance(&accepted);
    });
}

#[test]
fn g4_retry_cancel_and_caloric_failure_publish_no_partial_population() {
    with_cx(|cx| {
        let (thermal, heating) = setup(cx, [9.0, 0.0], [COUNT, 0.0, 0.0], 0.0);
        let mut candidate = heating.clone();
        candidate
            .step_sources(&thermal, 0.0, 2.0, 2.0, &mut || Ok(()))
            .unwrap();
        let retry = candidate
            .step_sources(&thermal, 0.0, 0.5, 0.5, &mut || Ok(()))
            .unwrap();
        let fresh = heating
            .clone()
            .step_sources(&thermal, 0.0, 0.5, 0.5, &mut || Ok(()))
            .unwrap();
        assert_eq!(retry, fresh);
        assert_eq!(candidate.populations(0), heating.populations(0));
        assert_eq!(candidate.energy(), heating.energy());
        assert!(matches!(
            thermal.advance_observed(&heating, 1.0, config(1.0), || Err(
                DepositionError::Cancelled
            )),
            Err(OpticalHeatingError::Thermal(DepositionError::Cancelled))
        ));
        assert!(
            thermal
                .advance_observed(&heating, 1e5, config(1e5), || Ok(()))
                .is_err()
        );
        assert_eq!(heating.populations(0), Some([COUNT, 0.0, 0.0]));
        assert_eq!(heating.time_s(), 0.0);
        assert_eq!(thermal.energy_j(), 0.0);
    });
}

#[test]
fn g0_explicit_pump_inputs_and_medium_support_are_admitted() {
    assert!(
        PrescribedPhotonPump::try_new(ContentHash([0; 32]), [400.0; 2], [1e-20; 2], [1e20; 2])
            .is_err()
    );
    assert!(
        PrescribedPhotonPump::try_new(ContentHash([1; 32]), [400.0; 2], [-1.0, 0.0], [1.0; 2])
            .is_err()
    );
    with_cx(|cx| {
        let density = [1.0];
        let medium = grid(cx, &density, MassAbsorptionSpectrum::constant(1.0).unwrap());
        let thermal = thermal(vec![cell(1, 1.0, 1.0, 0.0)], &SENSIBLE, vec![]);
        let triplet = SingletTripletSource::try_new(
            source(600.0, 1.0, [0.0, 1.0], 2),
            0.0,
            ArrheniusQuenching::try_new(0.0, 0.0).unwrap(),
            SingletTripletLimits::default(),
        )
        .unwrap();
        assert!(
            SaturableLuminescentHeating::try_new(
                &medium,
                source(500.0, 1.0, [1.0, 0.0], 1),
                triplet,
                &[],
                &thermal,
                1,
                1e-7,
                None,
                || Ok(())
            )
            .is_err()
        );
    });
}
