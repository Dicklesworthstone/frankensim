//! G0/G1/G3: analytic chemistry, conserved storage and source-preserving consumers.
use super::*;
use crate::thermochemical_gas::energy::equilibrium::tests::isomers;
use crate::thermochemical_gas::tests::{point, transport};
use fs_thermochem::equilibrium::composition::volume::EquilibriumVolumeSolveV1;
use fs_thermochem::{Composition, CompositionBasis, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R};

fn controls() -> GasIsochoricEquilibriumSolve {
    GasIsochoricEquilibriumSolve {
        equilibrium: EquilibriumInternalEnergySolveV1 {
            volume: EquilibriumVolumeSolveV1::default(),
            minimum_temperature_k: 200.0,
            maximum_temperature_k: 1500.0,
            energy_tolerance_j: 1e-6,
            maximum_temperature_iterations: 128,
        },
        mass_tolerance_kg: 1e-12,
    }
}

fn feed() -> Vec<(SpeciesId, Amount)> {
    vec![(SpeciesId::new("A").unwrap(), Amount::new(2.0))]
}

fn original_parcel(bank: &ThermochemicalGasModel) -> ClosedGasParcel {
    let pure = bank
        .with_composition(Composition::new(CompositionBasis::MoleFraction, vec![1.0, 0.0]).unwrap())
        .unwrap();
    ClosedGasParcel::new(
        pure,
        point(300.0, 3e5).with("sample_id", 7.0).unwrap(),
        2.0 * R * 300.0 / 3e5,
    )
    .unwrap()
}

#[test]
fn g1_uv_isomer_temperature_pressure_and_transport_share_one_conserved_state() {
    let (bank, chemistry) = isomers();
    let volume = Volume::new(2.0 * R * 300.0 / 3e5);
    let result = bank
        .equilibrium_chamber_uv(
            &chemistry,
            &point(250.0, 1e5).with("sample_id", 7.0).unwrap(),
            &feed(),
            Temperature::new(300.0),
            volume,
            controls(),
            || false,
        )
        .unwrap();
    // Independent equal-cv isomer balance. The B formation reference is
    // -24 kJ/mol and x_B=1/(1+exp(-24000/(R T))). UV uses cv=3R, not cp=4R.
    let (mut lo, mut hi) = (300.0, 1500.0);
    for _ in 0..80 {
        let t = 0.5 * (lo + hi);
        let x_b = 1.0 / (1.0 + (-24000.0 / (R * t)).exp());
        if 3.0 * R * (t - 300.0) - 24000.0 * x_b < 0.0 {
            lo = t;
        } else {
            hi = t;
        }
    }
    let parcel = result.parcel();
    let gas = parcel.state().state();
    let expected_t = 0.5 * (lo + hi);
    assert!((gas.temperature - expected_t).abs() < 1e-7);
    assert!((gas.pressure / (3e5 * expected_t / 300.0) - 1.0).abs() < 1e-10);
    assert_eq!(parcel.volume_m3(), volume.value());
    assert_eq!(parcel.point().axes()["sample_id"], 7.0);
    assert!((parcel.mass_kg() - 0.02).abs() < 1e-15);
    assert!((gas.density * parcel.volume_m3() - parcel.mass_kg()).abs() < 1e-15);
    assert!(result.balance().energy_residual_j.abs() < 1e-6);
    assert!(result.balance().freeze_energy_residual_j.abs() < 1e-6);
    assert!((result.balance().parcel_internal_energy_j - 2.0 * 3.0 * R * 300.0).abs() < 1e-6);
    assert!(gas.dynamic_viscosity > 1.7e-5);
    assert!(gas.thermal_conductivity > 0.025);
    assert!(
        (gas.sound_speed.powi(2) / (gas.gamma * gas.specific_gas_constant * gas.temperature) - 1.0)
            .abs()
            < 1e-14
    );
    assert_eq!(
        parcel.model().transport_materials(),
        bank.transport_materials()
    );
    let frozen_cv = parcel.mass_kg() * (gas.specific_heat_cp - gas.specific_gas_constant);
    assert!(
        result
            .equilibrium()
            .thermodynamic()
            .response()
            .heat_capacity_v_j_k
            > frozen_cv
    );
}

#[test]
fn g0_closed_parcel_equilibration_preserves_exact_mass_volume_and_original_sources() {
    let (bank, chemistry) = isomers();
    let mut parcel = original_parcel(&bank);
    let initial = parcel.clone();
    let energy = parcel.mass_kg() * parcel.state().specific_internal_energy_j_kg();
    let (equilibrium, balance) = parcel
        .equilibrate_adiabatic(&bank, &chemistry, controls(), || false)
        .unwrap();
    assert_eq!(parcel.mass_kg().to_bits(), initial.mass_kg().to_bits());
    assert_eq!(parcel.volume_m3().to_bits(), initial.volume_m3().to_bits());
    assert_eq!(parcel.point().axes()["sample_id"], 7.0);
    assert!(parcel.state().state().temperature > initial.state().state().temperature);
    assert!(parcel.state().state().pressure > initial.state().state().pressure);
    assert!(balance.energy_residual_j.abs() < 1e-6);
    assert!(
        (parcel.mass_kg() * parcel.state().specific_internal_energy_j_kg() - energy).abs() < 1e-6
    );
    assert_eq!(
        equilibrium.target_internal_energy(),
        Some(Energy::new(energy))
    );
    assert_eq!(
        parcel.model().thermodynamics().components(),
        chemistry.models()
    );
    let first = parcel.clone();
    parcel
        .equilibrate_adiabatic(&bank, &chemistry, controls(), || false)
        .unwrap();
    assert!((parcel.state().state().temperature - first.state().state().temperature).abs() < 1e-7);
    assert_eq!(parcel.mass_kg().to_bits(), initial.mass_kg().to_bits());
}

#[test]
fn g3_uv_extensive_scaling_and_catalog_fractions_do_not_change_intensive_physics() {
    let (bank, chemistry) = isomers();
    let mut reference: Option<(f64, f64)> = None;
    for scale in [0.25, 1.0, 8.0] {
        let feed = vec![(SpeciesId::new("A").unwrap(), Amount::new(2.0 * scale))];
        let other = bank
            .with_composition(
                Composition::new(CompositionBasis::MoleFraction, vec![0.9, 0.1]).unwrap(),
            )
            .unwrap();
        let result = other
            .equilibrium_chamber_uv(
                &chemistry,
                &point(300.0, 3e5),
                &feed,
                Temperature::new(300.0),
                Volume::new(scale * 2.0 * R * 300.0 / 3e5),
                controls(),
                || false,
            )
            .unwrap();
        let gas = result.parcel().state().state();
        if let Some((t, p)) = reference {
            assert!((gas.temperature - t).abs() < 1e-7);
            assert!((gas.pressure / p - 1.0).abs() < 1e-9);
        } else {
            reference = Some((gas.temperature, gas.pressure));
        }
    }
}

#[test]
fn g0_source_energy_and_late_cancellation_refusals_leave_the_whole_parcel_unchanged() {
    let (bank, chemistry) = isomers();
    let initial = original_parcel(&bank);
    let mut parcel = initial.clone();
    let mut wrong_transport = bank.transport_materials().to_vec();
    wrong_transport[0] = transport("A", 1.9e-5, 0.025);
    let wrong =
        ThermochemicalGasModel::new(bank.thermodynamics().clone(), wrong_transport, None).unwrap();
    assert!(matches!(
        parcel.equilibrate_adiabatic(&wrong, &chemistry, controls(), || false),
        Err(GasIsochoricEquilibriumError::Material(
            GasEquilibriumChamberError::SourceMismatch(_)
        ))
    ));
    assert_eq!(parcel, initial);
    let mut outside = controls();
    outside.equilibrium.maximum_temperature_k = 500.0;
    assert!(
        parcel
            .equilibrate_adiabatic(&bank, &chemistry, outside, || false)
            .is_err()
    );
    assert_eq!(parcel, initial);
    let mut polls = 0;
    let mut completed = initial.clone();
    completed
        .equilibrate_adiabatic(&bank, &chemistry, controls(), || {
            polls += 1;
            false
        })
        .unwrap();
    for stop in [1, polls / 2, polls] {
        let mut counter = 0;
        assert!(matches!(
            parcel.equilibrate_adiabatic(&bank, &chemistry, controls(), || {
                counter += 1;
                counter >= stop
            }),
            Err(GasIsochoricEquilibriumError::Cancelled)
        ));
        assert_eq!(parcel, initial);
    }
}

#[test]
fn g0_water_identity_is_retained_for_active_feed_and_initially_absent_products() {
    // These manufactured species test the explicit label's semantics only,
    // not the physical thermochemistry or condensation of water.
    let (bank, chemistry) = isomers();
    let water_a = SpeciesId::new("A").unwrap();
    let labeled_feed = ThermochemicalGasModel::new(
        bank.thermodynamics().clone(),
        bank.transport_materials().to_vec(),
        Some(&water_a),
    )
    .unwrap();
    let mut parcel = original_parcel(&labeled_feed);
    let before = parcel.clone();
    assert!(matches!(
        parcel.equilibrate_adiabatic(&bank, &chemistry, controls(), || false),
        Err(GasIsochoricEquilibriumError::Invalid(_))
    ));
    assert_eq!(parcel, before);
    let water_b = SpeciesId::new("B").unwrap();
    let labeled_product = ThermochemicalGasModel::new(
        bank.thermodynamics().clone(),
        bank.transport_materials().to_vec(),
        Some(&water_b),
    )
    .unwrap();
    let mut parcel = original_parcel(&labeled_product);
    assert_eq!(parcel.model().water_species(), None);
    parcel
        .equilibrate_adiabatic(&labeled_product, &chemistry, controls(), || false)
        .unwrap();
    assert_eq!(parcel.model().water_species(), Some(&water_b));
    assert!(parcel.state().state().water_mole_fraction > 0.5);
}
