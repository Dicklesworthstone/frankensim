//! Manufactured analytic checks; no empirical accuracy claim for named gases.
use super::*;
use crate::thermochemical_gas::energy::{
    GasEnergySolve,
    parcel::transfer::species::network::{
        SpeciesNetworkSolve,
        transport::{GasTransportPaths, GasTransportSolve},
    },
};
use crate::thermochemical_gas::tests::{model, point, species_with_energy_reference, transport};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialStateId, PropertyClaim, PropertyKey, PropertyValue,
    Provenance, UncertaintyModel,
};
use fs_thermochem::{Composition, CompositionBasis, mixture::FrozenIdealGasMixtureModelV1};
use std::collections::BTreeMap;

fn scalar(name: &str, dims: Dims, value: f64, maximum_t: f64) -> MaterialCard {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(name, dims),
            value: PropertyValue::Curve {
                abscissa: "temperature".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![(200.0, value), (maximum_t, value)],
                dims,
            },
            validity: ValidityDomain::unconstrained()
                .with_quantity(
                    "temperature",
                    QuantitySpec::semantic(SemanticType::new(
                        QuantityKind::AbsoluteTemperature,
                        ValueForm::Static,
                    )),
                    200.0,
                    maximum_t,
                )
                .with_quantity(
                    "pressure",
                    QuantitySpec::semantic(SemanticType::new(
                        QuantityKind::Pressure,
                        ValueForm::Static,
                    )),
                    1000.0,
                    1e7,
                ),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::LinearInside,
            observations: vec![],
            provenance: Provenance {
                source: "manufactured thermal diffusion".into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        })
        .unwrap();
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic".into(),
            phase: "gas".into(),
            process: name.into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn cards(heats: [f64; 2], conductivity: f64) -> [MaterialCard; 4] {
    [
        scalar(
            BINARY_DIFFUSIVITY_PROPERTY,
            BINARY_DIFFUSIVITY_DIMS,
            2e-5,
            1800.0,
        ),
        scalar(
            THERMAL_CONDUCTIVITY_PROPERTY,
            THERMAL_CONDUCTIVITY_DIMS,
            conductivity,
            1800.0,
        ),
        scalar(
            HEAT_OF_TRANSPORT_PROPERTY,
            HEAT_OF_TRANSPORT_DIMS,
            heats[0],
            1800.0,
        ),
        scalar(
            HEAT_OF_TRANSPORT_PROPERTY,
            HEAT_OF_TRANSPORT_DIMS,
            heats[1],
            1800.0,
        ),
    ]
}
fn prepared<'a>(
    bank: &ThermochemicalGasModel,
    cards: &'a [MaterialCard; 4],
    context: &'a QueryPoint,
) -> GasThermodiffusionModel<'a> {
    let source = |i| GasThermodiffusionMaterial {
        material: &cards[i],
        context,
        selection: MaterialPropertySelection::SingleClaimOnly,
    };
    GasThermodiffusionModel::new(
        bank,
        &[BinaryDiffusionMaterial {
            species: [SpeciesId::new("B").unwrap(), SpeciesId::new("A").unwrap()],
            material: &cards[0],
            context,
            selection: MaterialPropertySelection::SingleClaimOnly,
        }],
        source(1),
        &[
            SpeciesHeatOfTransportMaterial {
                species: SpeciesId::new("B").unwrap(),
                source: source(3),
            },
            SpeciesHeatOfTransportMaterial {
                species: SpeciesId::new("A").unwrap(),
                source: source(2),
            },
        ],
        || false,
    )
    .unwrap()
}
fn gas(bank: &ThermochemicalGasModel, x: f64, t: f64, p: f64) -> ClosedGasParcel {
    let composed = bank
        .with_composition(
            Composition::new(CompositionBasis::MoleFraction, vec![x, 1.0 - x]).unwrap(),
        )
        .unwrap();
    ClosedGasParcel::new(composed, point(t, p), 0.001).unwrap()
}
fn geometry() -> GasDiffusionGeometry {
    GasDiffusionGeometry {
        area_m2: 0.01,
        distance_m: 0.01,
    }
}
fn solve() -> GasThermodiffusionSolve {
    GasThermodiffusionSolve {
        maximum_temperature_difference_k: 1000.0,
        maximum_pressure_difference_pa: 1e6,
        diffusion: MaxwellStefanSolveV1::default(),
    }
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} != {expected} (tolerance {tolerance})"
    );
}

#[test]
fn binary_soret_matches_independent_mass_frame_solution_and_reciprocal_heat() {
    let bank = model(vec![0.4, 0.6], 0.001);
    let cards = cards([100000.0, 0.0], 0.02);
    let context = point(350.0, 90000.0);
    let law = prepared(&bank, &cards, &context);
    let mut work = law.workspace().unwrap();
    let (left, right) = (
        gas(&bank, 0.4, 300.0, 100000.0),
        gas(&bank, 0.4, 500.0, 100000.0),
    );
    let flux = law
        .face_flux(&left, &right, geometry(), solve(), &mut work, || false)
        .unwrap();
    let mean_mass = 0.4 * 0.028 + 0.6 * 0.004;
    let h_a = 100000.0 * (1.0 - 0.4 * 0.028 / mean_mass);
    let d_a = -0.4 * 0.028 * h_a * (1.0 / 500.0 - 1.0 / 300.0) / (R * geometry().distance_m);
    let z_a = -(100000.0 / R / 400.0) * 2e-5 * d_a * 0.004 / mean_mass;
    near(flux.species[0].molar_flux_mol_m2_s, z_a, z_a.abs() * 2e-12);
    near(
        flux.species[1].molar_flux_mol_m2_s,
        -z_a * 0.028 / 0.004,
        z_a.abs() * 2e-12,
    );
    near(
        flux.species
            .iter()
            .map(|s| s.mass_rate_left_to_right_kg_s)
            .sum(),
        0.0,
        1e-20,
    );
    near(
        flux.reciprocal_heat_flux_w_m2,
        100000.0 * z_a * 0.028,
        2e-11,
    );
    near(flux.fourier_heat_flux_w_m2, -400.0, 1e-12);
    // Independent integration for Cp/R=a+bT: Delta(h/T-s)=-R[a ln(Tr/Tl)+b DeltaT/2].
    let beta = 1.0 / 500.0 - 1.0 / 300.0;
    let log_t = (500.0_f64 / 300.0).ln();
    let effort = [-R * (3.5 * log_t + 0.5 * 0.001 * 200.0), -R * 2.5 * log_t];
    let upstream_h = [
        R * (3.5 * 500.0 + 0.5 * 0.001 * 500.0 * 500.0),
        R * 2.5 * 300.0,
    ];
    let expected_correction: f64 = (0..2)
        .map(|i| flux.species[i].molar_flux_mol_m2_s * (effort[i] / beta - upstream_h[i]))
        .sum();
    near(
        flux.enthalpy_correction_flux_w_m2,
        expected_correction,
        2e-11,
    );
    let carried: f64 = (0..2)
        .map(|i| flux.species[i].molar_flux_mol_m2_s * upstream_h[i])
        .sum();
    let chemical: f64 = (0..2)
        .map(|i| flux.species[i].molar_flux_mol_m2_s * effort[i])
        .sum();
    let paired_entropy =
        ((carried + flux.additional_heat_rate_left_to_right_w / geometry().area_m2) * beta
            - chemical)
            * geometry().area_m2;
    near(flux.entropy_production_w_k.unwrap(), paired_entropy, 2e-13);
    assert!(flux.entropy_production_w_k.unwrap() > 0.0);
    assert!(flux.maximum_gradient_residual_per_m < 1e-10);
    assert_eq!(
        flux.pairs[0].sample.face.state_coordinate("temperature"),
        Some(400.0)
    );
    assert_eq!(
        flux.pairs[0].sample.endpoints[0].state_coordinate("temperature"),
        Some(300.0)
    );
    assert_eq!(
        flux.pairs[0].sample.endpoints[1].state_coordinate("temperature"),
        Some(500.0)
    );
    assert_eq!(
        flux.pairs[0]
            .sample
            .face
            .property(BINARY_DIFFUSIVITY_PROPERTY)
            .unwrap()
            .answer()
            .receipt
            .selected,
        flux.pairs[0].sample.endpoints[0]
            .property(BINARY_DIFFUSIVITY_PROPERTY)
            .unwrap()
            .answer()
            .receipt
            .selected
    );
    let reverse = law
        .face_flux(&right, &left, geometry(), solve(), &mut work, || false)
        .unwrap();
    near(
        reverse.species[0].mass_rate_left_to_right_kg_s,
        -flux.species[0].mass_rate_left_to_right_kg_s,
        1e-20,
    );
    near(
        reverse.additional_heat_rate_left_to_right_w,
        -flux.additional_heat_rate_left_to_right_w,
        2e-12,
    );
}

#[test]
fn original_caloric_contrast_does_not_invent_soret_and_common_heat_offsets_cancel() {
    let bank = model(vec![0.4, 0.6], 0.001);
    let context = point(400.0, 100000.0);
    let uncoupled = cards([1e12, 1e12], 0.02);
    let law = prepared(&bank, &uncoupled, &context);
    let mut work = law.workspace().unwrap();
    let (left, right) = (
        gas(&bank, 0.4, 300.0, 100000.0),
        gas(&bank, 0.4, 700.0, 100000.0),
    );
    let zero = law
        .face_flux(&left, &right, geometry(), solve(), &mut work, || false)
        .unwrap();
    assert!(
        zero.species
            .iter()
            .all(|s| s.mass_rate_left_to_right_kg_s.abs() < 1e-25)
    );
    near(zero.additional_heat_rate_left_to_right_w, -8.0, 2e-12);
    let ordinary = cards([100000.0, -30000.0], 0.02);
    let shifted = cards([1e12 + 100000.0, 1e12 - 30000.0], 0.02);
    let a = prepared(&bank, &ordinary, &context)
        .face_flux(&left, &right, geometry(), solve(), &mut work, || false)
        .unwrap();
    let b = prepared(&bank, &shifted, &context)
        .face_flux(&left, &right, geometry(), solve(), &mut work, || false)
        .unwrap();
    assert_eq!(a.species, b.species);
    assert_eq!(
        a.additional_heat_rate_left_to_right_w,
        b.additional_heat_rate_left_to_right_w
    );
}

fn shifted_bank(offsets: [f64; 2]) -> ThermochemicalGasModel {
    ThermochemicalGasModel::new(
        FrozenIdealGasMixtureModelV1::new(
            vec![
                species_with_energy_reference("A", 0.028, 3.5, 0.001, offsets[0]),
                species_with_energy_reference("B", 0.004, 2.5, 0.0, offsets[1]),
            ],
            Composition::new(CompositionBasis::MoleFraction, vec![0.5, 0.5]).unwrap(),
        )
        .unwrap(),
        vec![transport("A", 1.7e-5, 0.025), transport("B", 2e-5, 0.15)],
        None,
    )
    .unwrap()
}

#[test]
fn formation_references_cancel_before_heat_and_force_arithmetic() {
    let cards = cards([100000.0, -30000.0], 0.02);
    let context = point(400.0, 100000.0);
    let mut results = Vec::new();
    for bank in [shifted_bank([0.0, 0.0]), shifted_bank([1e24, -3e24])] {
        let law = prepared(&bank, &cards, &context);
        let result = law
            .face_flux(
                &gas(&bank, 0.6, 300.0, 110000.0),
                &gas(&bank, 0.2, 700.0, 90000.0),
                geometry(),
                solve(),
                &mut law.workspace().unwrap(),
                || false,
            )
            .unwrap();
        results.push(result);
    }
    assert_eq!(results[0], results[1]);
}

#[test]
fn common_zero_stays_absent_positive_traces_survive_and_one_sided_zero_refuses() {
    let bank = model(vec![0.5, 0.5], 0.001);
    let cards = cards([100000.0, 0.0], 0.0);
    let context = point(400.0, 100000.0);
    let law = prepared(&bank, &cards, &context);
    let mut work = law.workspace().unwrap();
    let pure = gas(&bank, 0.0, 300.0, 100000.0);
    let pure_hot = gas(&bank, 0.0, 500.0, 100000.0);
    let result = law
        .face_flux(&pure, &pure_hot, geometry(), solve(), &mut work, || false)
        .unwrap();
    assert!(
        result
            .species
            .iter()
            .all(|s| s.mass_rate_left_to_right_kg_s == 0.0)
    );
    let trace = gas(&bank, 1e-100, 300.0, 100000.0);
    let trace_hot = gas(&bank, 1e-100, 500.0, 100000.0);
    let result = law
        .face_flux(&trace, &trace_hot, geometry(), solve(), &mut work, || false)
        .unwrap();
    assert!(result.species[0].mass_rate_left_to_right_kg_s < 0.0);
    assert!(
        law.face_flux(&pure, &trace_hot, geometry(), solve(), &mut work, || false)
            .is_err()
    );
    let mut calls = 0;
    assert_eq!(
        law.face_flux(&trace, &trace_hot, geometry(), solve(), &mut work, || {
            calls += 1;
            calls == 4
        }),
        Err(GasEnergyError::Cancelled)
    );
    assert_eq!(
        law.face_flux(&trace, &trace_hot, geometry(), solve(), &mut work, || false)
            .unwrap(),
        result
    );
}

#[test]
fn source_admission_checks_original_endpoints_not_only_the_average_face() {
    let bank = model(vec![0.5, 0.5], 0.001);
    let mut cards = cards([100000.0, 0.0], 0.02);
    cards[2] = scalar(
        HEAT_OF_TRANSPORT_PROPERTY,
        HEAT_OF_TRANSPORT_DIMS,
        100000.0,
        450.0,
    );
    let context = point(400.0, 100000.0);
    let law = prepared(&bank, &cards, &context);
    let result = law.face_flux(
        &gas(&bank, 0.5, 300.0, 100000.0),
        &gas(&bank, 0.5, 500.0, 100000.0),
        geometry(),
        solve(),
        &mut law.workspace().unwrap(),
        || false,
    );
    assert!(result.is_err()); // The 400 K face is supported; the 500 K original trace is not.
    let replacement = shifted_bank([10.0, 0.0]);
    assert!(
        law.admit_endpoints(
            &gas(&replacement, 0.5, 350.0, 100000.0),
            &gas(&replacement, 0.5, 400.0, 100000.0),
            solve(),
            || false
        )
        .is_err()
    );
}

fn energy_solve() -> GasEnergySolve {
    GasEnergySolve {
        minimum_temperature_k: 200.0,
        maximum_temperature_k: 1800.0,
        energy_tolerance_j_kg: 1e-7,
        maximum_iterations: 128,
    }
}
fn network_solve() -> GasTransportSolve {
    GasTransportSolve {
        exchange: SpeciesNetworkSolve {
            maximum_species_fraction: 0.2,
            mass_tolerance_kg: 1e-14,
            energy_tolerance_j: 1e-8,
        },
        entropy_tolerance_j_k: 1e-10,
    }
}

#[test]
fn finite_network_combines_soret_reciprocal_heat_and_original_enthalpy_once() {
    let bank = model(vec![0.4, 0.6], 0.001);
    let cards = cards([100000.0, 0.0], 0.02);
    let context = point(400.0, 100000.0);
    let law = prepared(&bank, &cards, &context);
    let original = BTreeMap::from([
        (1, gas(&bank, 0.4, 300.0, 100000.0)),
        (2, gas(&bank, 0.4, 500.0, 100000.0)),
        (3, gas(&bank, 0.3, 450.0, 100000.0)),
    ]);
    let saved = original.clone();
    let links = [
        network::GasThermodiffusionLink {
            id: 2,
            parcels: [2, 3],
            geometry: geometry(),
            model: &law,
            solve: solve(),
        },
        network::GasThermodiffusionLink {
            id: 1,
            parcels: [1, 2],
            geometry: geometry(),
            model: &law,
            solve: solve(),
        },
    ];
    let endpoints = BTreeMap::from([
        (1, energy_solve()),
        (2, energy_solve()),
        (3, energy_solve()),
    ]);
    let mut scratch = [law.workspace().unwrap(), law.workspace().unwrap()];
    let candidate = network::step_candidate(
        &original,
        &links,
        &mut scratch,
        GasTransportPaths::default(),
        0.001,
        &endpoints,
        network_solve(),
        || false,
    )
    .unwrap();
    assert!(candidate.transport().step().entropy_change_j_k > 0.0);
    assert_eq!(
        candidate
            .links()
            .iter()
            .map(|row| row.id)
            .collect::<Vec<_>>(),
        vec![1, 2]
    );
    assert!(candidate.parcels()[&1].state().state().temperature > 300.0);
    assert!(candidate.parcels()[&2].state().state().temperature < 500.0);
    assert_ne!(
        candidate.parcels()[&1]
            .model()
            .thermodynamics()
            .mole_composition(),
        original[&1].model().thermodynamics().mole_composition()
    );
    for row in &candidate.transport().step().exchange.parcels {
        assert!(row.energy_residual_j.abs() < 1e-8);
        assert!(row.maximum_species_mass_residual_kg < 1e-14);
    }
    // The actual heat record equals the correction; donor enthalpy is owned elsewhere.
    let expected = candidate
        .links()
        .iter()
        .find(|row| row.id == 1)
        .unwrap()
        .flux
        .additional_heat_rate_left_to_right_w
        * 0.001;
    let first = candidate
        .transport()
        .step()
        .exchange
        .parcels
        .iter()
        .find(|row| row.id == 1)
        .unwrap();
    near(first.net_additional_heat_j, -expected, 1e-14);
    assert!(
        network::step_candidate(
            &original,
            &links,
            &mut scratch,
            GasTransportPaths::default(),
            1e8,
            &endpoints,
            network_solve(),
            || false
        )
        .is_err()
    );
    assert_eq!(
        network::step_candidate(
            &original,
            &links,
            &mut scratch,
            GasTransportPaths::default(),
            0.001,
            &endpoints,
            network_solve(),
            || true
        ),
        Err(network::GasThermodiffusionNetworkError::Cancelled)
    );
    assert_eq!(original, saved);
}
