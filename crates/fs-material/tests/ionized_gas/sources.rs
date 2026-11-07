//! G0/G1/G3: source-complete transport across exact population boundaries.
use super::*;
use fs_material::ionized_gas::{
    CarrierMomentumSource, IonizedGasMomentumSource, IonizedGasRefreshError,
    IonizedGasTransportBank,
};
use fs_math::c64::C64;
use fs_qty::{Amount, Temperature, Volume};
use fs_thermochem::electron_collision::MaxwellianRateControlsV1;
use fs_thermochem::equilibrium::composition::{FrozenIonizedGasStateV1, IdealGasEquilibriumV1};
use fs_thermochem::{ChargeVector, ElementId, ElementalMatrix};

fn electron_source(targets: &[(&str, f64)]) -> IonizedGasMomentumSource {
    IonizedGasMomentumSource::ElectronNeutral {
        cross_sections: targets
            .iter()
            .map(|&(name, sigma)| momentum_cross_section(name, sigma))
            .collect(),
        controls: MaxwellianRateControlsV1 {
            relative_tolerance: 1e-10,
            maximum_subdivisions: 4096,
        },
        maximum_target_thermal_speed_ratio: 0.1,
    }
}

fn ion_source(state: &FrozenIonizedGasStateV1) -> IonizedGasMomentumSource {
    let (template, original_card) = collision_source_and_card(
        "H+",
        &[
            ("T_h", state.heavy_temperature().value()),
            ("T_e", state.electron_temperature().value()),
        ],
        state.pressure().value(),
        1e9,
    );
    IonizedGasMomentumSource::Material {
        template,
        original_card,
    }
}

fn bank(state: &FrozenIonizedGasStateV1, coulomb: bool) -> IonizedGasTransportBank {
    IonizedGasTransportBank::new(
        state,
        vec![electron_source(&[("H", 2e-20)]), ion_source(state)],
        coulomb.then(|| coulomb_model(state)),
        limits(),
        || false,
    )
    .unwrap()
}

#[test]
fn neutral_birth_extinction_and_trace_carriers_reach_native_electrical_response() {
    let neutral = frozen_state(0.0, 6000.0, 12000.0, 1.0);
    let sources = bank(&neutral, false);
    let retained = sources.clone();
    assert_eq!(sources.evaluate(&neutral, || false).unwrap(), None);
    assert_eq!(sources.sources().len(), 2);
    for x in [0.001, 0.003, 1e-18] {
        let state = frozen_state(x, 6500.0, 15000.0, 1.0);
        let response = sources.evaluate(&state, || false).unwrap().unwrap();
        assert_eq!(response.state(), &state);
        assert!(response.equilibrium().is_none());
        assert_eq!(response.carriers().len(), 2);
        // Independent shifted-Maxwellian momentum moment for constant sigma:
        // (2/3)<sigma v E/kT> = (4/3)*sigma*sqrt(8 kT/(pi me)).
        let momentum =
            (4.0 / 3.0) * 2e-20 * (8.0 * KB * 15000.0 / (core::f64::consts::PI * ME)).sqrt();
        let electron_nu = (1.0 - x) * NA * momentum;
        let electron_sigma = x * NA * E * E / (ME * electron_nu);
        let ion_sigma = x * NA * E * E / ((MH - ME) * 1e9);
        relative(
            response.conductivity_s_m(),
            electron_sigma + ion_sigma,
            2e-9,
        );
        let efield = response.maximum_electric_field_v_m() * 0.01;
        let current = response
            .ohmic_conductivity()
            .current_density_a_m2([efield, 0.0]);
        let heat = response.joule_heating_w_m3([efield, 0.0]).unwrap();
        relative(heat.total_power_density_w_m3(), efield * current[0], 2e-14);
        // The result is the native carrier owner, including settled inertia.
        let omega = 3e8;
        let harmonic = response.prepare_harmonic(|| false).unwrap();
        let actual = harmonic
            .response(omega, [C64::from_re(efield), C64::ZERO, C64::ZERO], || {
                false
            })
            .unwrap();
        let expected = C64::from_re(electron_sigma * electron_nu) / C64::new(electron_nu, omega)
            + C64::from_re(ion_sigma * 1e9) / C64::new(1e9, omega);
        relative(actual.conductivity_s_m.re, expected.re, 2e-9);
        relative(actual.conductivity_s_m.im, expected.im, 2e-9);
    }
    assert_eq!(sources.evaluate(&neutral, || false).unwrap(), None);
    assert_eq!(sources, retained);
}

fn buffered_state(neutral_h: f64) -> FrozenIonizedGasStateV1 {
    let mut models = chemistry().models().to_vec();
    // Manufactured inert monatomic buffer; H nuclei and He nuclei remain
    // independently conserved throughout this test's recombination endpoints.
    models.insert(2, standard_state("He", 4.0 * MH, 1.0, 0.0));
    let names: Vec<_> = models.iter().map(|m| m.species().clone()).collect();
    let elements = ElementalMatrix::new(
        vec![ElementId::new("H").unwrap(), ElementId::new("He").unwrap()],
        names.clone(),
        vec![vec![1, 1, 0, 0], vec![0, 0, 1, 0]],
    )
    .unwrap();
    let chemistry = IdealGasEquilibriumV1::new(
        elements,
        ChargeVector::new(names, vec![0, 1, 0, -1]).unwrap(),
        models,
    )
    .unwrap();
    let ion = 0.001 - neutral_h;
    FrozenIonizedGasStateV1::new(
        &chemistry,
        &[
            Amount::new(neutral_h),
            Amount::new(ion),
            Amount::new(1.0),
            Amount::new(ion),
        ],
        Volume::new(1.0),
        Temperature::new(6000.0),
        Temperature::new(12000.0),
        1e-12,
        || false,
    )
    .unwrap()
}

#[test]
fn recombination_activates_new_neutral_targets_from_the_original_complete_bank() {
    let initial = buffered_state(0.0);
    let sources = IonizedGasTransportBank::new(
        &initial,
        vec![
            ion_source(&initial),
            electron_source(&[("He", 3e-20), ("H", 2e-20)]),
        ],
        None,
        limits(),
        || false,
    )
    .unwrap();
    for neutral_h in [0.0, 0.0005, 0.0] {
        let state = buffered_state(neutral_h);
        let response = sources.evaluate(&state, || false).unwrap().unwrap();
        let electron = response
            .carriers()
            .iter()
            .find(|c| c.is_electron())
            .unwrap();
        let coefficient = (4.0 / 3.0) * (8.0 * KB * 12000.0 / (core::f64::consts::PI * ME)).sqrt();
        relative(
            electron.relaxation().rate_per_s(),
            NA * (neutral_h * 2e-20 + 3e-20) * coefficient,
            2e-9,
        );
        let CarrierMomentumSource::ElectronNeutral { contributions, .. } =
            electron.relaxation().source()
        else {
            panic!("original cross-section law must remain retained");
        };
        assert_eq!(contributions.len(), if neutral_h == 0.0 { 1 } else { 2 });
        assert!(
            contributions
                .iter()
                .all(|c| c.number_density_per_m3() > 0.0)
        );
    }
    assert_eq!(
        sources.evaluate(&buffered_state(0.001), || false).unwrap(),
        None
    );
}

#[test]
fn initially_neutral_coulomb_prescription_rebinds_actual_created_carriers() {
    let neutral = frozen_state(0.0, 6000.0, 12000.0, 1.0);
    let sources = bank(&neutral, true);
    assert_eq!(sources.evaluate(&neutral, || false).unwrap(), None);
    for x in [0.001, 0.003] {
        let state = frozen_state(x, 6500.0, 15000.0, 1.0);
        let result = sources.evaluate(&state, || false).unwrap().unwrap();
        let independently_rebuilt = IonizedGasConductivity::from_state(
            state.clone(),
            vec![
                source_two_temperature("H+", 6500.0, 15000.0, state.pressure().value(), 1e9),
                electron_relaxation(&state, 2e-20),
            ],
            limits(),
            || false,
        )
        .unwrap()
        .with_coulomb_drag(coulomb_model(&state), || false)
        .unwrap();
        assert!(
            result
                .coulomb_model()
                .unwrap()
                .is_bound_to_populations(&state)
        );
        assert_eq!(result.coulomb_pairs().len(), 1);
        relative(
            result.conductivity_s_m(),
            independently_rebuilt.conductivity_s_m(),
            2e-13,
        );
        let field = [result.maximum_electric_field_v_m() * 0.01, 0.0];
        assert_eq!(
            result.joule_heating_w_m3(field),
            independently_rebuilt.joule_heating_w_m3(field)
        );
    }
    assert_eq!(sources.evaluate(&neutral, || false).unwrap(), None);
}

#[test]
fn incomplete_foreign_duplicate_and_unsupported_banks_refuse_without_mutation() {
    let neutral = frozen_state(0.0, 6000.0, 12000.0, 1.0);
    let good_electron = electron_source(&[("H", 2e-20)]);
    for sources in [
        vec![ion_source(&neutral)],
        vec![ion_source(&neutral), ion_source(&neutral)],
        vec![ion_source(&neutral), electron_source(&[])],
        vec![ion_source(&neutral), electron_source(&[("H+", 2e-20)])],
        vec![
            ion_source(&neutral),
            electron_source(&[("H", 2e-20), ("H", 2e-20)]),
        ],
    ] {
        assert!(IonizedGasTransportBank::new(&neutral, sources, None, limits(), || false).is_err());
    }
    let mut wrong_card = ion_source(&neutral);
    if let IonizedGasMomentumSource::Material { original_card, .. } = &mut wrong_card {
        *original_card = collision_source_and_card(
            "H+",
            &[("T_h", 6000.0), ("T_e", 12000.0)],
            neutral.pressure().value(),
            2e9,
        )
        .1;
    }
    assert!(
        IonizedGasTransportBank::new(
            &neutral,
            vec![wrong_card, good_electron],
            None,
            limits(),
            || false,
        )
        .is_err()
    );
    let mut domain = limits();
    domain.maximum_charged_mole_fraction = 0.0;
    assert!(
        IonizedGasTransportBank::new(
            &neutral,
            vec![ion_source(&neutral), electron_source(&[("H", 2e-20)])],
            None,
            domain,
            || false,
        )
        .is_err()
    );
    let sources = bank(&neutral, false);
    let retained = sources.clone();
    assert!(
        sources
            .evaluate(&frozen_state(0.001, 6000.0, 12000.0, 2.0), || false)
            .is_err()
    );
    assert!(
        sources
            .evaluate(&frozen_state(0.2, 6000.0, 12000.0, 1.0), || false)
            .is_err()
    );
    let next = frozen_state(0.001, 6500.0, 15000.0, 1.0);
    assert_eq!(
        sources.evaluate(&next, || true),
        Err(IonizedGasRefreshError::Cancelled)
    );
    let mut polls = 0;
    assert_eq!(
        sources.evaluate(&next, || {
            polls += 1;
            polls > 4
        }),
        Err(IonizedGasRefreshError::Cancelled)
    );
    assert_eq!(sources, retained);
}
