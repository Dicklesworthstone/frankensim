//! Original native capacities, finite-body heat balance and source admission.
use super::*;
use crate::thermal::transport::conduction::pair::{
    self as finite, FiniteCaloricBody, FiniteCaloricExchangeControls,
    FiniteCaloricExchangeError,
};
use crate::phase::exchange::pair::FinitePhaseExchangeError;

fn body<'l, 'a, 'm>(
    layer: &'l CaloricConductionLayer<'a, 'm>,
    initial_specific_energy_j_kg: f64,
    mass_kg: f64,
) -> FiniteCaloricBody<'l, 'a, 'm> {
    FiniteCaloricBody { layer, initial_specific_energy_j_kg, mass_kg }
}

// Independent two-body integral when BOTH total capacities are K*(T+400).
// Energy conservation gives (Th+400)^2+(Tc+400)^2=S; integrate using their sum.
fn affine_pair_exposure(coefficient: f64) -> f64 {
    let invariant: f64 = 1_665_000.0;
    let a = (invariant / 2.0).sqrt();
    let initial_sum = 1800.0;
    let final_sum = 940.0 + (invariant - 940.0_f64.powi(2)).sqrt();
    coefficient / 2.0 * (
        a * ((final_sum / (2.0 * a)).atanh() - (initial_sum / (2.0 * a)).atanh())
            - (final_sum - initial_sum)
    )
}

#[test]
fn g1_native_cv_cp_pair_matches_independent_nonlinear_oracle_and_reciprocity() {
    let (cv_source, cp_source) = (Fixture::new(true), Fixture::new(false));
    let (cv, cp) = (cv_source.cv(), cp_source.cp());
    let cv_k = IsochoricConductionLaw::try_new(
        &cv_source.card, &cv, &cv_source.conductivity, SELECTION,
    ).unwrap();
    let cp_k = IsobaricConductionLaw::try_new(
        &cp_source.card, &cp, &cp_source.conductivity, SELECTION,
    ).unwrap();
    let cv_layer = CaloricConductionLayer::from_isochoric(&cv_k).unwrap();
    let cp_layer = CaloricConductionLayer::from_isobaric(&cp_k).unwrap();
    let donor = body(&cv_layer, energy(0, 650.0), 0.8);
    let receiver = body(&cp_layer, energy(1, 350.0), 0.4);
    let controls = FiniteCaloricExchangeControls::default();
    let x = affine_pair_exposure(0.4);
    let result = finite::exchange(donor, receiver, x, controls).unwrap();
    near(result.heat_first_to_second_j, 43_780.0, 2e-5);
    near(result.temperatures_k[0], 540.0, 8e-8);
    near(result.temperatures_k[1], (1_665_000.0 - 940.0_f64.powi(2)).sqrt() - 400.0, 8e-8);
    assert_eq!(result.storage_bases, [
        ThermalStorageBasis::InternalEnergy, ThermalStorageBasis::Enthalpy,
    ]);
    assert_eq!(result.caloric_identities, [cv.caloric_identity(), cp.caloric_identity()]);
    assert!(result.root_iterations > 0 && result.quadrature_panels > 0);
    assert!(result.heat_error_bound_j <= controls.absolute_heat_tolerance_j
        + controls.relative_heat_tolerance * result.heat_first_to_second_j);
    for side in 0..2 {
        let supplied = if side == 0 { -result.heat_first_to_second_j } else { result.heat_first_to_second_j };
        assert_eq!(result.energy_residual_j[side], result.represented_energy_change_j[side] - supplied);
        assert_eq!([&cv_layer, &cp_layer][side].state_at_specific_energy(
            result.final_specific_energies_j_kg[side],
        ).unwrap().temperature_k, result.temperatures_k[side]);
    }
    let reversed = finite::exchange(receiver, donor, x, controls).unwrap();
    near(reversed.heat_first_to_second_j, -result.heat_first_to_second_j, 2e-8);
    near(reversed.temperatures_k[0], result.temperatures_k[1], 2e-10);
    let first = finite::exchange(donor, receiver, x / 2.0, controls).unwrap();
    let second = finite::exchange(
        body(&cv_layer, first.final_specific_energies_j_kg[0], 0.8),
        body(&cp_layer, first.final_specific_energies_j_kg[1], 0.4),
        x / 2.0, controls,
    ).unwrap();
    near(first.heat_first_to_second_j + second.heat_first_to_second_j,
        result.heat_first_to_second_j, 4e-5);
    near(second.temperatures_k[0], result.temperatures_k[0], 1e-7);
}

#[test]
fn g1_native_pair_scales_original_volume_amount_and_free_expansion_once() {
    let (source, keys) = additional_sources();
    let cv_source = Fixture::new(true);
    let (cp, cv) = (source.cp(), cv_source.cv());
    let cp_k = IsobaricConductionLaw::try_new(
        &source.card, &cp, &source.conductivity, SELECTION,
    ).unwrap();
    let cv_k = IsochoricConductionLaw::try_new(
        &cv_source.card, &cv, &cv_source.conductivity, SELECTION,
    ).unwrap();
    let cp_layer = CaloricConductionLayer::from_isobaric(&cp_k).unwrap();
    let cv_layer = CaloricConductionLayer::from_isochoric(&cv_k).unwrap();
    let volume = PreparedVolumetricThermalTransport::try_new(
        &source.card, VolumetricThermalPathQuery {
            volumetric_heat_capacity: &keys[0], conductivity: &source.conductivity,
            temperature_axis: "T", lower: &source.lower, upper: &source.upper,
            selections: [SELECTION; 2],
        },
    ).unwrap();
    let volume_law = FixedVolumeThermalLaw::try_new(&volume, 750.0, 300.0).unwrap();
    let volume_layer = CaloricConductionLayer::from_fixed_volume(&volume_law).unwrap();
    let storage = PreparedHeatStorage::try_new(&[
        HeatStorageComponentQuery {
            id: 1, card: &source.card, capacity: &source.capacity,
            amount: HeatStorageAmount::MassKg(1.5), temperature_axis: "T",
            lower: &source.lower, upper: &source.upper, selection: SELECTION,
        },
        HeatStorageComponentQuery {
            id: 2, card: &source.card, capacity: &keys[0],
            amount: HeatStorageAmount::ReferenceVolumeM3(0.002), temperature_axis: "T",
            lower: &source.lower, upper: &source.upper, selection: SELECTION,
        },
    ]).unwrap();
    let assembled = HeatStorageThermalLaw::try_new(&storage, 3.0, 300.0).unwrap();
    let assembled_k = HeatStorageConductionLaw::try_new(
        &source.card, &assembled, &source.conductivity, "T", &source.lower, &source.upper, SELECTION,
    ).unwrap();
    let assembled_layer = CaloricConductionLayer::from_heat_storage(&assembled_k).unwrap();
    let density = PreparedFreeExpansionDensity::try_new(
        &source.card, FreeExpansionSourceQuery {
            density: &keys[1], density_selection: SELECTION, reference: &source.lower,
            strains: FreeExpansionStrainSources::Isotropic { strain: &keys[2], selection: SELECTION },
            temperature_axis: "T", lower: &source.lower, upper: &source.upper,
        },
    ).unwrap();
    let expansion = IsobaricFreeExpansion::try_new(&source.card, &cp, &density).unwrap();
    let expanding_k = IsobaricExpansionConductionLaw::try_new(&expansion, &cp_k).unwrap();
    let expanding_layer = CaloricConductionLayer::from_expanding_isobaric(&expanding_k).unwrap();
    let donor = body(&cv_layer, energy(0, 650.0), 6.0);
    let controls = FiniteCaloricExchangeControls::default();
    for (layer, basis) in [
        (&cp_layer, ThermalStorageBasis::Enthalpy),
        (&volume_layer, ThermalStorageBasis::DeclaredSensibleStorage),
        (&assembled_layer, ThermalStorageBasis::DeclaredSensibleStorage),
        (&expanding_layer, ThermalStorageBasis::Enthalpy),
    ] {
        let result = finite::exchange(
            donor, body(layer, energy(1, 350.0), 3.0), affine_pair_exposure(3.0), controls,
        ).unwrap();
        near(result.heat_first_to_second_j, 328_350.0, 1e-4);
        near(result.temperatures_k[0], 540.0, 8e-8);
        assert_eq!(result.storage_bases[1], basis);
        assert!(result.energy_residual_j.iter().all(|e| e.abs() < 1e-7));
    }
    // An amount-bound assembly cannot silently become a different inventory.
    assert!(finite::exchange(donor, body(&assembled_layer, energy(1, 350.0), 1.5),
        0.0, controls).is_err());
    assert_eq!(assembled.specific_entropy_change_j_kg_k(
        energy(1, 350.0), energy(1, 500.0),
    ).unwrap(), None);
}

#[test]
fn g4_native_phase_pair_retains_original_inverse_and_narrow_conductivity_support() {
    use crate::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, SolidLiquidPhase};
    use crate::thermal::{conduction::source::ResolvedConductionLayer, source::ThermalMaterialPoint};
    use crate::phase::exchange::pair::{self as phase_pair, FinitePhaseBody};
    let source = Fixture::new(false);
    // Explicit authored constant-C h/T source; the material card supplies k.
    let curve = EquilibriumEnthalpyPhaseCurve::try_single_phase(
        source.card.content_hash(),
        SolidLiquidPhase::Solid,
        vec![
            EnthalpyPhaseKnot { specific_enthalpy_j_kg: 0.0, temperature_k: 300.0,
                liquid_mass_fraction: 0.0, bulk_density_kg_m3: 1000.0 },
            EnthalpyPhaseKnot { specific_enthalpy_j_kg: 400_000.0, temperature_k: 700.0,
                liquid_mass_fraction: 0.0, bulk_density_kg_m3: 1000.0 },
        ],
    ).unwrap();
    let point = ThermalMaterialPoint {
        material: &source.card, context: &source.lower, temperature_axis: "T",
        liquid_fraction_axis: None,
    };
    let full = ResolvedConductionLayer::new(
        &curve, point, &source.conductivity, SELECTION, SolidLiquidPhase::Solid,
        [300.0, 700.0], || Ok(()),
    ).unwrap();
    let narrow = ResolvedConductionLayer::new(
        &curve, point, &source.conductivity, SELECTION, SolidLiquidPhase::Solid,
        [325.0, 425.0], || Ok(()),
    ).unwrap();
    let (full_layer, narrow_layer) = (
        CaloricConductionLayer::from_equilibrium_phase(&full),
        CaloricConductionLayer::from_equilibrium_phase(&narrow),
    );
    let donor = body(&full_layer, 350_000.0, 1.0);
    let receiver = body(&narrow_layer, 50_000.0, 1.0);
    let controls = FiniteCaloricExchangeControls::default();
    let short = finite::exchange(donor, receiver, 1.0, controls).unwrap();
    near(short.heat_first_to_second_j, 150_000.0 * -(-0.002_f64).exp_m1(), 2e-8);
    assert_eq!(short.root_iterations, 0);
    assert_eq!(short.quadrature_panels, 0);
    let endpoint = finite::exchange(donor, receiver, 500.0 * 2.0_f64.ln(), controls).unwrap();
    near(endpoint.temperatures_k[1], 425.0, 2e-10);
    near(endpoint.final_specific_energies_j_kg[1], narrow.enthalpy_range_j_kg()[1], 2e-7);
    assert!(matches!(finite::exchange(donor, receiver, 600.0, controls),
        Err(FiniteCaloricExchangeError::Integration(FinitePhaseExchangeError::OutsideSourceSupport))));
    // The unchanged full-chart API can reach the same native state beyond
    // the selected layer's conductivity support and retains original entropy.
    let unrestricted = phase_pair::exchange(
        FinitePhaseBody { curve: &curve, initial: curve.state_at_specific_enthalpy(350_000.0).unwrap(), mass_kg: 1.0 },
        FinitePhaseBody { curve: &curve, initial: curve.state_at_specific_enthalpy(50_000.0).unwrap(), mass_kg: 1.0 },
        600.0, controls,
    ).unwrap();
    assert!(unrestricted.states[1].temperature_k() > 425.0);
    assert!(unrestricted.total_entropy_change_j_k > 0.0);
}

#[test]
fn g4_native_pair_preserves_weak_heat_residuals_and_whole_call_refusals() {
    let (cv_source, cp_source) = (Fixture::new(true), Fixture::new(false));
    let (cv, cp) = (cv_source.cv(), cp_source.cp());
    let cv_k = IsochoricConductionLaw::try_new(
        &cv_source.card, &cv, &cv_source.conductivity, SELECTION,
    ).unwrap();
    let cp_k = IsobaricConductionLaw::try_new(
        &cp_source.card, &cp, &cp_source.conductivity, SELECTION,
    ).unwrap();
    let cv_layer = CaloricConductionLayer::from_isochoric(&cv_k).unwrap();
    let cp_layer = CaloricConductionLayer::from_isobaric(&cp_k).unwrap();
    let donor = body(&cv_layer, energy(0, 650.0), 0.8);
    let receiver = body(&cp_layer, 0.0, 0.4);
    let controls = FiniteCaloricExchangeControls::default();
    let weak = finite::exchange(donor, receiver, 1e-20, controls).unwrap();
    near(weak.heat_first_to_second_j, 350.0e-20, 2e-31);
    assert_eq!(weak.temperatures_k, [650.0, 300.0]);
    assert_eq!(weak.represented_energy_change_j[0], 0.0);
    assert_eq!(weak.energy_residual_j[0], weak.heat_first_to_second_j);
    assert!(weak.final_specific_energies_j_kg[1] > 0.0);
    near(weak.energy_residual_j[1], 0.0, 2e-32);
    assert!(matches!(finite::exchange(donor, receiver, 1.0,
        FiniteCaloricExchangeControls { maximum_chart_intervals: 1, ..controls }),
        Err(FiniteCaloricExchangeError::Integration(FinitePhaseExchangeError::WorkLimit))));
    for x in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(finite::exchange(donor, receiver, x, controls).is_err());
    }
    let mut calls = 0;
    let expected = finite::exchange_with_checkpoint(donor, receiver, 1.0, controls, || {
        calls += 1; false
    }).unwrap();
    assert!(calls > 5);
    let mut replay = 0;
    assert!(matches!(
        finite::exchange_with_checkpoint(donor, receiver, 1.0, controls, || {
            replay += 1; replay == calls
        }),
        Err(FiniteCaloricExchangeError::Integration(FinitePhaseExchangeError::Cancelled))
    ));
    assert_eq!(finite::exchange(donor, receiver, 1.0, controls).unwrap(), expected);
}

