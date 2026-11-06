//! G0/G1/G3/G4: reacting isentropes, physical material closure and atomic refusal.
use super::*;
use crate::thermochemical_gas::energy::equilibrium::tests::isomers;
use crate::thermochemical_gas::tests::{point, transport_in_domain};
use fs_qty::Temperature;
use fs_thermochem::equilibrium::composition::volume::EquilibriumVolumeSolveV1;
use fs_thermochem::{
    Composition, CompositionBasis, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R,
};

fn controls() -> GasIsentropicVolumeSolve {
    let mut volume = EquilibriumVolumeSolveV1::default();
    volume.equilibrium.relative_balance_tolerance = 1e-12;
    volume.equilibrium.dimensionless_potential_tolerance = 1e-12;
    volume.relative_volume_tolerance = 1e-12;
    GasIsentropicVolumeSolve {
        equilibrium: EquilibriumVolumeEntropySolveV1 {
            volume,
            minimum_temperature_k: 200.0,
            maximum_temperature_k: 1500.0,
            entropy_tolerance_j_k: 1e-8,
            maximum_temperature_iterations: 128,
        },
        mass_tolerance_kg: 1e-12,
        freeze_energy_tolerance_j: 1e-6,
    }
}

fn feed(amount: f64) -> Vec<(SpeciesId, Amount)> {
    vec![(SpeciesId::new("A").unwrap(), Amount::new(amount))]
}

fn mixing_entropy(t: f64) -> f64 {
    let y = 1.0 / (1.0 + (-24000.0 / (R * t)).exp());
    -y * y.ln() - (1.0 - y) * (1.0 - y).ln()
}

#[test]
fn g1_reacting_compression_and_expansion_retain_original_material_sources() {
    let (bank, chemistry) = isomers();
    let original = chemistry
        .equilibrate_tv(
            Temperature::new(900.0),
            Volume::new(0.02),
            &feed(2.0),
            controls().equilibrium.volume,
            || false,
        )
        .unwrap();
    let entropy = original.thermodynamic().entropy_j_k();
    let context = point(300.0, 1e5).with("sample_id", 7.0).unwrap();
    for ratio in [0.5, 1.0, 2.0] {
        let result = bank
            .equilibrium_chamber_sv(
                &chemistry,
                &context,
                &feed(2.0),
                entropy,
                Volume::new(0.02 * ratio),
                controls(),
                || false,
            )
            .unwrap();
        // Independent scalar isomer entropy root, not a second call to SV.
        // Cp=4R for both species; N is unchanged and Cv_frozen=3R.
        let (mut low, mut high) = (200.0_f64, 1500.0_f64);
        for _ in 0..80 {
            let t = 0.5 * (low + high);
            let residual = 3.0 * (t / 900.0).ln() + ratio.ln()
                + mixing_entropy(t) - mixing_entropy(900.0);
            if residual < 0.0 {
                low = t;
            } else {
                high = t;
            }
        }
        let expected = 0.5 * (low + high);
        let parcel = result.parcel();
        let gas = parcel.state().state();
        assert!((gas.temperature - expected).abs() < 2e-6);
        assert!((gas.pressure / (2.0 * R * expected / (0.02 * ratio)) - 1.0).abs() < 1e-8);
        assert_eq!(parcel.volume_m3(), 0.02 * ratio);
        assert!((parcel.mass_kg() - 0.02).abs() < 1e-14);
        assert!((gas.density * parcel.volume_m3() - parcel.mass_kg()).abs() < 1e-14);
        assert_eq!(parcel.point().axes()["sample_id"], 7.0);
        assert_eq!(parcel.model().thermodynamics().components(), chemistry.models());
        assert_eq!(parcel.model().transport_materials(), bank.transport_materials());
        assert!(gas.dynamic_viscosity > 0.0 && gas.thermal_conductivity > 0.0);
        assert!(result.balance().entropy_residual_j_k.abs() <= 1e-8);
        assert!(result.balance().freeze_entropy_residual_j_k.abs() <= 1e-8);
        assert!(result.balance().freeze.freeze_energy_residual_j.abs() <= 1e-6);
        assert_eq!(result.equilibrium().target_internal_energy(), None);
        if ratio < 1.0 {
            assert!(gas.temperature < 900.0 * ratio.powf(-1.0 / 3.0));
            assert!(
                result.equilibrium().internal_energy().value()
                    > original.internal_energy().value()
            );
        }
    }
}

#[test]
fn g3_extensive_scaling_and_catalog_fractions_preserve_the_sv_material_state() {
    let (bank, chemistry) = isomers();
    let alternate = bank
        .with_composition(
            Composition::new(CompositionBasis::MoleFraction, vec![0.9, 0.1]).unwrap(),
        )
        .unwrap();
    let initial = chemistry
        .equilibrate_tv(
            Temperature::new(900.0),
            Volume::new(0.02),
            &feed(2.0),
            controls().equilibrium.volume,
            || false,
        )
        .unwrap();
    let entropy = initial.thermodynamic().entropy_j_k();
    let mut reference: Option<(f64, f64)> = None;
    for scale in [0.25, 1.0, 8.0] {
        let mut limits = controls();
        limits.equilibrium.entropy_tolerance_j_k *= scale;
        limits.mass_tolerance_kg *= scale;
        limits.freeze_energy_tolerance_j *= scale;
        let result = alternate
            .equilibrium_chamber_sv(
                &chemistry,
                &point(300.0, 1e5),
                &feed(2.0 * scale),
                entropy * scale,
                Volume::new(0.01 * scale),
                limits,
                || false,
            )
            .unwrap();
        let gas = result.parcel().state().state();
        if let Some((t, p)) = reference {
            assert!((gas.temperature - t).abs() < 2e-6);
            assert!((gas.pressure / p - 1.0).abs() < 1e-8);
        } else {
            reference = Some((gas.temperature, gas.pressure));
        }
        assert!((result.parcel().mass_kg() - 0.02 * scale).abs() < limits.mass_tolerance_kg);
    }
}

#[test]
fn g0_g4_source_refusal_and_cancellation_publish_no_sv_material_result() {
    let (bank, chemistry) = isomers();
    let before = bank.clone();
    let initial = chemistry
        .equilibrate_tv(
            Temperature::new(900.0),
            Volume::new(0.02),
            &feed(2.0),
            controls().equilibrium.volume,
            || false,
        )
        .unwrap();
    let entropy = initial.thermodynamic().entropy_j_k();
    let context = point(300.0, 1e5);
    let mut polls = 0;
    bank.equilibrium_chamber_sv(
        &chemistry, &context, &feed(2.0), entropy, Volume::new(0.01), controls(),
        || { polls += 1; false },
    ).unwrap();
    for stop in [1, polls / 2, polls] {
        let mut count = 0;
        assert!(matches!(
            bank.equilibrium_chamber_sv(
                &chemistry, &context, &feed(2.0), entropy, Volume::new(0.01), controls(),
                || { count += 1; count >= stop },
            ),
            Err(GasIsentropicVolumeError::Cancelled)
        ));
    }
    for tolerance in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        let mut invalid = controls();
        invalid.freeze_energy_tolerance_j = tolerance;
        assert!(matches!(
            bank.equilibrium_chamber_sv(
                &chemistry, &context, &feed(2.0), entropy, Volume::new(0.01), invalid,
                || false,
            ),
            Err(GasIsentropicVolumeError::Invalid(_))
        ));
    }
    // Thermochemistry is supported, but the original transport source is not
    // valid at the solved temperature. Do not substitute a convenient material.
    let limited = ThermochemicalGasModel::new(
        bank.thermodynamics().clone(),
        vec![
            transport_in_domain("A", 1.7e-5, 0.025, [200.0, 350.0], [1000.0, 1e7]),
            transport_in_domain("B", 2.0e-5, 0.15, [200.0, 350.0], [1000.0, 1e7]),
        ],
        None,
    ).unwrap();
    assert!(matches!(
        limited.equilibrium_chamber_sv(
            &chemistry, &context, &feed(2.0), entropy, Volume::new(0.01), controls(),
            || false,
        ),
        Err(GasIsentropicVolumeError::Volume(
            GasIsochoricEquilibriumError::Material(GasEquilibriumChamberError::State(_))
        ))
    ));
    let incomplete = bank
        .with_composition(
            Composition::new(CompositionBasis::MoleFraction, vec![1.0, 0.0]).unwrap(),
        ).unwrap();
    assert!(incomplete.equilibrium_chamber_sv(
        &chemistry, &context, &feed(2.0), entropy, Volume::new(0.01), controls(),
        || false,
    ).is_err());
    assert_eq!(bank, before);
}
