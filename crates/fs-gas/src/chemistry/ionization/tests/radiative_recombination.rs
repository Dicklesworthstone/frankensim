use super::{EV, NA, R, initial};
use crate::chemistry::radiative_recombination::*;
use fs_matdb::Provenance;
use fs_qty::{Area, Energy, Temperature, Time, Volume};
use fs_thermochem::equilibrium::composition::FrozenIonizedGasStateV1;
use fs_thermochem::kinetics::atomic_ionization::photoionization::{
    AtomicPhotoionizationV1, PhotoionizationCrossSectionPointV1,
};
use fs_thermochem::kinetics::atomic_ionization::radiative_recombination::AtomicRadiativeRecombinationV1;

fn fixture(
    x: f64,
    te: f64,
    elemental: f64,
    charge: f64,
) -> (FrozenIonizedGasStateV1, AtomicRadiativeRecombinationV1) {
    let (original, channel) = initial(x, te, 6000.0, elemental, charge);
    let state = FrozenIonizedGasStateV1::new(
        original.thermodynamics(),
        original.amounts(),
        Volume::new(0.2),
        original.heavy_temperature(),
        original.electron_temperature(),
        1e-12,
        || false,
    )
    .unwrap();
    let gap = channel.ionization_energy_j_mol() / NA;
    let provenance = Provenance {
        source: "manufactured constant bound-free cross section and global constant envelope"
            .into(),
        license: "analytic test fixture".into(),
        artifact: None,
    };
    let photon = AtomicPhotoionizationV1::new(
        channel,
        provenance.clone(),
        vec![
            PhotoionizationCrossSectionPointV1 {
                photon_energy: Energy::new(gap),
                cross_section: Area::new(1e-22),
            },
            PhotoionizationCrossSectionPointV1 {
                photon_energy: Energy::new(100.0 * EV),
                cross_section: Area::new(1e-22),
            },
        ],
    )
    .unwrap();
    let source = AtomicRadiativeRecombinationV1::new(photon, Area::new(1e-22), provenance).unwrap();
    (state, source)
}
fn controls(tolerance: f64) -> RadiativeRecombinationSolve {
    let mut controls =
        RadiativeRecombinationSolve::new(Temperature::new(1000.0), Temperature::new(50000.0));
    controls.relative_tolerance = tolerance;
    controls
}
fn evolve(
    state: &FrozenIonizedGasStateV1,
    source: &AtomicRadiativeRecombinationV1,
    seconds: f64,
    tolerance: f64,
) -> RadiativeRecombinationResult {
    advance_radiative_recombination(
        state,
        source,
        Time::new(seconds),
        controls(tolerance),
        || false,
    )
    .unwrap()
}

#[test]
fn g1_radiative_capture_matches_independent_coupled_population_and_temperature_oracles() {
    let (initial, source) = fixture(0.1, 10000.0, 0.0, 0.0);
    // Independent DOP853 integration of extent, temperature and both escaping
    // energy rates, with closed-form constant-cross-section Maxwellian moments.
    // The omitted 100 eV tail is <2e-37 of either moment in these cases.
    for (time, extent, te, kinetic, radiation) in [
        (
            2e-5,
            0.015031597039235413,
            10411.80209805202,
            1438.308126718066,
            21160.50718621601,
        ),
        (
            4e-5,
            0.026003850514376876,
            10769.88899624853,
            2532.622943375727,
            36650.96145519997,
        ),
    ] {
        let result = evolve(&initial, &source, time, 1e-9);
        let state = result.state();
        let balance = result.balance();
        assert!((result.recombination_extent().value() - extent).abs() < 7e-10);
        assert!((state.electron_temperature().value() - te).abs() < 7e-4);
        assert!((balance.captured_electron_energy_j - kinetic).abs() < 7e-4);
        assert!((balance.emitted_radiation_energy_j - radiation).abs() < 7e-4);
        assert_eq!(state.amounts()[1], state.amounts()[2]);
        assert!((state.amounts()[0].value() + state.amounts()[1].value() - 1.0).abs() < 2e-14);
        assert!((state.mass().value() / initial.mass().value() - 1.0).abs() < 2e-14);
        assert_eq!(state.volume(), initial.volume());
        assert_eq!(state.heavy_temperature(), initial.heavy_temperature());
        assert_eq!(state.thermodynamics(), initial.thermodynamics());
        assert_eq!(result.elapsed_time(), Time::new(time));
        assert!(state.electron_temperature().value() > initial.electron_temperature().value());
        assert!(balance.electron_sensible_energy_change_j < 0.0);
        assert!(balance.chemical_storage_change_j < 0.0);
        assert!(
            (balance.emitted_photons / NA - result.recombination_extent().value()).abs() < 1e-16
        );
        assert!(balance.energy_residual_j.abs() < 1e-7);
        assert!(
            (state.internal_energy().value() - initial.internal_energy().value()
                + balance.emitted_radiation_energy_j)
                .abs()
                < 2e-7
        );
        assert!(balance.maximum_normalized_local_error <= 1.0);
    }
}

#[test]
fn g1_refinement_reduces_coupled_error_without_clipping_the_temperature() {
    let (initial, source) = fixture(0.1, 10000.0, 0.0, 0.0);
    let coarse = evolve(&initial, &source, 4e-5, 1e-5);
    let fine = evolve(&initial, &source, 4e-5, 1e-7);
    let expected = 10769.88899624853;
    let coarse_error = (coarse.state().electron_temperature().value() - expected).abs();
    let fine_error = (fine.state().electron_temperature().value() - expected).abs();
    assert!(fine_error < 0.15 * coarse_error);
    assert!(fine.balance().accepted_steps > coarse.balance().accepted_steps);
}

#[test]
fn g3_escaping_energy_is_invariant_to_conserved_formation_references() {
    let (initial, source) = fixture(0.1, 10000.0, 0.0, 0.0);
    let base = evolve(&initial, &source, 2e-5, 1e-7);
    let (shifted, shifted_source) = fixture(0.1, 10000.0, 1e12, 1e8);
    let changed = evolve(&shifted, &shifted_source, 2e-5, 1e-7);
    assert!(
        (base.recombination_extent().value() - changed.recombination_extent().value()).abs() < 1e-9
    );
    assert!(
        (base.state().electron_temperature().value()
            - changed.state().electron_temperature().value())
        .abs()
            < 1e-5
    );
    assert!(
        (base.balance().emitted_radiation_energy_j - changed.balance().emitted_radiation_energy_j)
            .abs()
            < 1e-4
    );
    assert!(changed.balance().energy_residual_j.abs() < 1e-7);
    let end = changed.state();
    let sensible = 1.5
        * R
        * (end.amounts()[2].value() * end.electron_temperature().value()
            - shifted.amounts()[2].value() * shifted.electron_temperature().value());
    assert!((sensible + changed.balance().captured_electron_energy_j).abs() < 1e-7);
}

#[test]
fn g0_empty_reactants_and_zero_duration_are_exact_identities() {
    let (neutral, source) = fixture(0.0, 10000.0, 0.0, 0.0);
    let result = evolve(&neutral, &source, 1.0, 1e-7);
    assert_eq!(result.state(), &neutral);
    assert_eq!(result.recombination_extent().value(), 0.0);
    assert_eq!(result.balance().emitted_radiation_energy_j, 0.0);
    let (charged, source) = fixture(0.1, 10000.0, 0.0, 0.0);
    assert_eq!(evolve(&charged, &source, 0.0, 1e-7).state(), &charged);
}

#[test]
fn g4_invalid_inputs_work_exhaustion_source_binding_and_cancellation_refuse() {
    let (initial, source) = fixture(0.1, 10000.0, 0.0, 0.0);
    for seconds in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(
            advance_radiative_recombination(
                &initial,
                &source,
                Time::new(seconds),
                controls(1e-7),
                || false,
            )
            .is_err()
        );
    }
    let mut bounded = controls(1e-7);
    bounded.maximum_steps = 1;
    assert!(matches!(
        advance_radiative_recombination(&initial, &source, Time::new(4e-5), bounded, || false,),
        Err(RadiativeRecombinationError::StepBudget)
    ));
    let mut bounded = controls(1e-7);
    bounded.maximum_electron_temperature_k = 10000.0;
    bounded.minimum_step_s = 1e-10;
    assert!(matches!(
        advance_radiative_recombination(&initial, &source, Time::new(4e-5), bounded, || false,),
        Err(RadiativeRecombinationError::OutsideDomain)
    ));
    let (different, _) = fixture(0.1, 10000.0, 1e6, 0.0);
    assert!(matches!(
        advance_radiative_recombination(
            &different,
            &source,
            Time::new(1e-5),
            controls(1e-7),
            || false,
        ),
        Err(RadiativeRecombinationError::Invalid(_))
    ));
    assert_eq!(
        advance_radiative_recombination(&initial, &source, Time::new(1e-5), controls(1e-7), || {
            true
        },),
        Err(RadiativeRecombinationError::Cancelled)
    );
    let mut polls = 0;
    let cancelled =
        advance_radiative_recombination(&initial, &source, Time::new(1e-5), controls(1e-7), || {
            polls += 1;
            polls >= 20
        });
    assert_eq!(cancelled, Err(RadiativeRecombinationError::Cancelled));
    assert_eq!(polls, 20);
}
