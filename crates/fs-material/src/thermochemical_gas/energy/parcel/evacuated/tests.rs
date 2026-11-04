//! G0/G1/G3/G4: exact-vacuum storage transition through original gas laws.
//! Sources are synthetic; these are caloric/conservation checks, not gas data.

use super::*;
use crate::thermochemical_gas::tests::{point, species_with_energy_reference, transport_in_domain};
use fs_thermochem::{
    Composition, CompositionBasis, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R,
    mixture::FrozenIdealGasMixtureModelV1,
};

const VOLUME: f64 = 0.001;
const MASS_TOLERANCE: f64 = 1e-13;
const ENERGY_TOLERANCE: f64 = 1e-8;

fn catalog(
    rows: &[(&str, f64, f64, f64)],
    fractions: Vec<f64>,
    datum: f64,
) -> ThermochemicalGasModel {
    ThermochemicalGasModel::new(
        FrozenIdealGasMixtureModelV1::new(
            rows.iter()
                .enumerate()
                .map(|(i, &(id, m, cp, slope))| {
                    species_with_energy_reference(id, m, cp, slope, datum * (i + 1) as f64)
                })
                .collect(),
            Composition::new(CompositionBasis::MoleFraction, fractions).unwrap(),
        )
        .unwrap(),
        rows.iter()
            .map(|&(id, _, _, _)| {
                transport_in_domain(id, 1.7e-5, 0.025, [200.0, 2000.0], [1.0, 1e7])
            })
            .collect(),
        None,
    )
    .unwrap()
}

fn controls() -> GasEnergySolve {
    GasEnergySolve {
        minimum_temperature_k: 200.0,
        maximum_temperature_k: 1800.0,
        energy_tolerance_j_kg: 1e-7,
        maximum_iterations: 128,
    }
}

fn receiver() -> EvacuatedGasVolume {
    EvacuatedGasVolume::new(
        VOLUME,
        QueryPoint::new().with("receiver-lot", 17.0).unwrap(),
    )
    .unwrap()
}

fn masses(rows: &[(&str, f64)]) -> BTreeMap<SpeciesId, f64> {
    rows.iter()
        .map(|&(id, value)| (SpeciesId::new(id).unwrap(), value))
        .collect()
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}"
    );
}

#[test]
fn g1_first_fill_matches_molecular_energy_for_monatomic_and_diatomic_gases() {
    let empty = receiver();
    let received = 2e-5;
    let molar_mass = 0.028;
    let source_temperature = 400.0;
    for cp_r in [2.5, 3.5] {
        let model = catalog(&[("A", molar_mass, cp_r, 0.0)], vec![1.0], 0.0);
        // An effusive molecule carries u(T_source) + R*T_source/2 per mole.
        // Complete thermalization at conserved mass therefore heats the receiver
        // to T_source * (Cv + R/2)/Cv. This does not use the inversion as oracle.
        let energy = received * R * source_temperature / (2.0 * molar_mass);
        let result = empty
            .fill_candidate(
                &model,
                &masses(&[("A", received)]),
                source_temperature,
                energy,
                controls(),
                MASS_TOLERANCE,
                ENERGY_TOLERANCE,
                || false,
            )
            .unwrap();
        let expected_t = source_temperature * (cp_r - 0.5) / (cp_r - 1.0);
        let expected_p = received / molar_mass * R * expected_t / VOLUME;
        let gas = result.gas();
        near(gas.state().state().temperature, expected_t, 1e-8);
        near(gas.state().state().pressure, expected_p, 1e-6);
        assert_eq!(gas.mass_kg(), received);
        assert_eq!(gas.volume_m3(), VOLUME);
        assert_eq!(gas.model(), &model);
        assert_eq!(gas.point().axes().get("receiver-lot"), Some(&17.0));
        assert!(result.step().energy_residual_j.abs() <= ENERGY_TOLERANCE);
        assert!(result.step().density_mass_residual_kg.abs() <= MASS_TOLERANCE);
        assert!(result.step().maximum_species_mass_residual_kg <= MASS_TOLERANCE);
        near(
            result.step().stored_relative_energy_j,
            energy,
            ENERGY_TOLERANCE,
        );
    }
    assert!(!empty.context().axes().contains_key("temperature"));
    assert!(!empty.context().axes().contains_key("pressure"));
}

#[test]
fn g3_received_mixture_uses_original_variable_calorics_independent_of_energy_reference() {
    let rows = [
        ("A", 0.028, 3.5, 0.001),
        ("B", 0.004, 2.5, 0.0),
        ("C", 0.032, 3.5, 0.0),
    ];
    let received = masses(&[("B", 4e-6), ("C", 0.0), ("A", 3e-5)]);
    let empty = receiver();
    let target_t: f64 = 650.0;
    // Independent integral of each specified Cv(T), using RECEIVED mass rather
    // than the source catalog's unrelated 1/4, 1/4, 1/2 mole fractions.
    let energy = |reference: f64| {
        3e-5 * R / 0.028
            * (2.5 * (target_t - reference) + 0.0005 * (target_t.powi(2) - reference.powi(2)))
            + 4e-6 * R / 0.004 * 1.5 * (target_t - reference)
    };
    for datum in [0.0, 1e12] {
        let model = catalog(&rows, vec![0.25, 0.25, 0.5], datum);
        let original = model.clone();
        for reference in [300.0, 400.0] {
            let result = empty
                .fill_candidate(
                    &model,
                    &received,
                    reference,
                    energy(reference),
                    controls(),
                    MASS_TOLERANCE,
                    ENERGY_TOLERANCE,
                    || false,
                )
                .unwrap();
            let gas = result.gas();
            near(gas.state().state().temperature, target_t, 1e-8);
            near(
                gas.state().state().pressure,
                (3e-5 / 0.028 + 4e-6 / 0.004) * R * target_t / VOLUME,
                1e-6,
            );
            let actual = gas.species_inventory_kg().unwrap();
            near(actual[&SpeciesId::new("A").unwrap()], 3e-5, MASS_TOLERANCE);
            near(actual[&SpeciesId::new("B").unwrap()], 4e-6, MASS_TOLERANCE);
            assert!(!actual.contains_key(&SpeciesId::new("C").unwrap()));
            assert_eq!(
                gas.model().thermodynamics().components(),
                &model.thermodynamics().components()[..2]
            );
            assert_eq!(
                gas.model().transport_materials(),
                &model.transport_materials()[..2]
            );
            assert_eq!(gas.point().axes().get("receiver-lot"), Some(&17.0));
            assert!(result.step().energy_residual_j.abs() <= ENERGY_TOLERANCE);
        }
        assert_eq!(model, original);
    }
}

#[test]
fn g0_zero_relative_energy_uses_the_declared_chart_without_a_seed_gas() {
    let empty = receiver();
    let model = catalog(&[("A", 0.028, 3.5, 0.001)], vec![1.0], 1e12);
    let result = empty
        .fill_candidate(
            &model,
            &masses(&[("A", 2e-5)]),
            400.0,
            0.0,
            controls(),
            MASS_TOLERANCE,
            ENERGY_TOLERANCE,
            || false,
        )
        .unwrap();
    assert_eq!(result.gas().state().state().temperature, 400.0);
    assert_eq!(result.step().stored_relative_energy_j, 0.0);
    assert_eq!(result.step().energy_residual_j, 0.0);
    assert_eq!(result.step().iterations, 0);
    assert_eq!(empty, receiver());
}

#[test]
fn g4_invalid_inventory_source_support_and_late_cancellation_publish_no_gas() {
    for volume in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(EvacuatedGasVolume::new(volume, QueryPoint::new()).is_err());
    }
    for context in [
        QueryPoint::new().with("temperature", 0.0).unwrap(),
        QueryPoint::new().with("pressure", 0.0).unwrap(),
        point(400.0, 100000.0),
    ] {
        assert!(EvacuatedGasVolume::new(VOLUME, context).is_err());
    }
    let model = catalog(&[("A", 0.028, 3.5, 0.001)], vec![1.0], 0.0);
    let saved_model = model.clone();
    let empty = receiver();
    let saved_empty = empty.clone();
    let valid = masses(&[("A", 2e-5)]);
    for invalid in [
        BTreeMap::new(),
        masses(&[("A", 0.0)]),
        masses(&[("A", -1.0)]),
        masses(&[("A", f64::NAN)]),
        masses(&[("A", f64::INFINITY)]),
        masses(&[("unknown", 2e-5)]),
    ] {
        assert!(
            empty
                .fill_candidate(
                    &model,
                    &invalid,
                    400.0,
                    0.1,
                    controls(),
                    MASS_TOLERANCE,
                    ENERGY_TOLERANCE,
                    || false
                )
                .is_err()
        );
    }
    assert!(
        empty
            .fill_candidate(
                &model,
                &valid,
                190.0,
                0.1,
                controls(),
                MASS_TOLERANCE,
                ENERGY_TOLERANCE,
                || false
            )
            .is_err()
    );
    for energy in [f64::NAN, f64::INFINITY, 1e20] {
        assert!(
            empty
                .fill_candidate(
                    &model,
                    &valid,
                    400.0,
                    energy,
                    controls(),
                    MASS_TOLERANCE,
                    ENERGY_TOLERANCE,
                    || false
                )
                .is_err()
        );
    }
    assert!(
        empty
            .fill_candidate(
                &model,
                &valid,
                400.0,
                0.1,
                controls(),
                0.0,
                ENERGY_TOLERANCE,
                || false
            )
            .is_err()
    );
    // Thermodynamics can solve this state, but its actual positive pressure is
    // below the retained transport sources. A chart pressure cannot admit it.
    let unsupported = EvacuatedGasVolume::new(1e9, QueryPoint::new()).unwrap();
    assert!(matches!(
        unsupported.fill_candidate(
            &model,
            &valid,
            400.0,
            0.1,
            controls(),
            MASS_TOLERANCE,
            ENERGY_TOLERANCE,
            || false
        ),
        Err(GasEnergyError::Material(_))
    ));
    let mut total_polls = 0;
    empty
        .fill_candidate(
            &model,
            &valid,
            400.0,
            0.1,
            controls(),
            MASS_TOLERANCE,
            ENERGY_TOLERANCE,
            || {
                total_polls += 1;
                false
            },
        )
        .unwrap();
    let mut polls = 0;
    let result = empty.fill_candidate(
        &model,
        &valid,
        400.0,
        0.1,
        controls(),
        MASS_TOLERANCE,
        ENERGY_TOLERANCE,
        || {
            polls += 1;
            polls == total_polls
        },
    );
    assert!(matches!(result, Err(GasEnergyError::Cancelled)));
    assert_eq!(polls, total_polls);
    assert_eq!(empty, saved_empty);
    assert_eq!(model, saved_model);
}
