//! Original source-owner synthetic fixtures; no measured-gas accuracy claim.
use super::*;
use crate::gas_field::{AxialGasField, PrimitiveGasState};
use crate::gas_scattering::AVOGADRO_PER_MOL;
use crate::gas_species::input::ThermochemicalPhase;
use crate::thermochemical_gas::ThermochemicalGasModel;
use crate::thermochemical_gas::tests::{model, point, species_in_domain, transport};
use fs_blake3::ContentHash;
use fs_matdb::QueryPoint;
use fs_thermochem::mixture::FrozenIdealGasMixtureModelV1;
use fs_thermochem::{Composition, CompositionBasis, Nasa9StandardStateModelV1};

fn optical(species: &str, refractivity: f64, delta: f64, source: u8) -> SpeciesMolecularRayleigh {
    SpeciesMolecularRayleigh {
        species: SpeciesId::new(species).unwrap(),
        optical: MolecularRayleigh::try_new(
            refractivity,
            2.5e25,
            550.0,
            delta,
            ContentHash([source; 32]),
        )
        .unwrap(),
    }
}

fn optics() -> [SpeciesMolecularRayleigh; 2] {
    [optical("A", 0.00015, 0.0, 1), optical("B", 0.0006, 0.3, 2)]
}

fn near(got: f64, expected: f64) {
    assert!(
        (got - expected).abs() < 5.0e-13 * expected.abs().max(1.0e-300),
        "{got:.17e} != {expected:.17e}"
    );
}

#[test]
fn g1_attachment_uses_original_nasa_state_molar_mass_and_number_fractions() {
    let source = model(vec![0.75, 0.25], 0.001);
    let accepted = source.evaluate(&point(700.0, 150_000.0)).unwrap();
    let expected = accepted.clone();
    let original_components = accepted.components().as_ptr();
    let [a, b] = optics();
    let attached = ThermochemicalRayleighState::try_new(accepted, &[b.clone(), a.clone()]).unwrap();
    assert_eq!(attached.resolved(), &expected);
    assert_eq!(
        attached.resolved().components().as_ptr(),
        original_components
    );
    assert_eq!(attached.rayleigh().gas(), expected.state());
    assert_eq!(attached.optical_sources(), &[a.clone(), b.clone()]);
    assert_ne!(
        expected.state().gamma,
        source
            .evaluate(&point(300.0, 150_000.0))
            .unwrap()
            .state()
            .gamma
    );
    let mass = expected.thermodynamics().properties().molar_mass().value();
    assert_eq!(attached.rayleigh().molar_mass_kg_mol(), mass);
    let volume =
        0.75 * a.optical.lorentz_lorenz_volume_m3() + 0.25 * b.optical.lorentz_lorenz_volume_m3();
    let mass_refractivity = volume / mass * AVOGADRO_PER_MOL;
    near(
        attached.rayleigh().mass_refractivity_m3_kg(),
        mass_refractivity,
    );
    let ll = expected.state().density * mass_refractivity;
    for wavelength in [380.0, 417.3, 550.0, 709.7, 780.0] {
        near(
            attached
                .rayleigh()
                .refractive_index(wavelength)
                .unwrap()
                .phase_index,
            ((1.0 + 2.0 * ll) / (1.0 - ll)).sqrt(),
        );
        let sigma = 0.75 * a.optical.cross_section_m2(wavelength).unwrap()
            + 0.25 * b.optical.cross_section_m2(wavelength).unwrap();
        near(
            attached
                .rayleigh()
                .mass_scattering_m2_kg(wavelength)
                .unwrap(),
            sigma / mass * AVOGADRO_PER_MOL,
        );
        near(
            attached.rayleigh().scattering_per_m(wavelength).unwrap(),
            expected.state().density * sigma / mass * AVOGADRO_PER_MOL,
        );
        for cosine in [-1.0, -0.4, 0.0, 0.7, 1.0] {
            let differential: f64 = [(0.75, &a), (0.25, &b)]
                .iter()
                .map(|(fraction, row)| {
                    let delta = row.optical.depolarization_ratio();
                    let phase = 3.0 * ((1.0 + delta) + (1.0 - delta) * cosine * cosine)
                        / (8.0 * core::f64::consts::PI * (2.0 + delta));
                    fraction * row.optical.cross_section_m2(wavelength).unwrap() * phase
                })
                .sum();
            near(
                sigma * attached.rayleigh().phase_per_sr(cosine).unwrap(),
                differential,
            );
        }
    }
}

fn field(left: f64, right: f64) -> AxialGasField {
    let composition =
        |a| Composition::new(CompositionBasis::MoleFraction, vec![a, 1.0 - a]).unwrap();
    AxialGasField::new(
        vec![0.0, 2.0],
        101_325.0,
        vec![
            PrimitiveGasState {
                temperature_k: 300.0,
                composition: composition(left),
            },
            PrimitiveGasState {
                temperature_k: 900.0,
                composition: composition(right),
            },
        ],
        model(vec![0.75, 0.25], 0.001),
        QueryPoint::new(),
        ThermochemicalPhase::FrozenGas,
    )
    .unwrap()
}

#[test]
fn g3_axial_field_samples_attach_directly_and_retain_composition_changes() {
    let field = field(0.75, 0.75);
    let optics = optics();
    let states: Vec<_> = [0.0, 1.0, 2.0]
        .map(|x| {
            let accepted = field.state_at(x).unwrap();
            let expected = accepted.clone();
            let attached = ThermochemicalRayleighState::try_new(accepted, &optics).unwrap();
            assert_eq!(attached.resolved(), &expected);
            attached
        })
        .into_iter()
        .collect();
    drop(field);
    for state in &states {
        assert_eq!(
            state.rayleigh().model_identity(),
            states[0].rayleigh().model_identity()
        );
        assert_eq!(
            state.rayleigh().mass_scattering_reference_m2_kg(),
            states[0].rayleigh().mass_scattering_reference_m2_kg()
        );
        near(
            state.rayleigh().scattering_per_m(550.0).unwrap(),
            states[0].rayleigh().scattering_per_m(550.0).unwrap() * 300.0
                / state.resolved().state().temperature,
        );
    }
    assert_ne!(
        states[0].rayleigh().content_identity(),
        states[1].rayleigh().content_identity()
    );
    let changing = self::field(0.75, 0.25);
    let changed =
        ThermochemicalRayleighState::try_new(changing.state_at(1.0).unwrap(), &optics).unwrap();
    assert_ne!(
        changed.rayleigh().model_identity(),
        states[1].rayleigh().model_identity()
    );
    assert_ne!(
        changed.rayleigh().mass_scattering_reference_m2_kg(),
        states[1].rayleigh().mass_scattering_reference_m2_kg()
    );
}

#[test]
fn g3_temperature_region_changes_keep_fixed_optics_and_complete_original_sources() {
    let low = species_in_domain("A", 0.028, 3.5, 0.001, 0.0, [200.0, 500.0]);
    let high = species_in_domain("A", 0.028, 4.0, 0.001, 50_000.0, [500.0, 2000.0]);
    let species = Nasa9StandardStateModelV1::new(
        low.species().clone(),
        low.molar_mass(),
        low.convention().clone(),
        vec![low.regions()[0].clone(), high.regions()[0].clone()],
    )
    .unwrap();
    let source = ThermochemicalGasModel::new(
        FrozenIdealGasMixtureModelV1::new(
            vec![species],
            Composition::new(CompositionBasis::MoleFraction, vec![1.0]).unwrap(),
        )
        .unwrap(),
        vec![transport("A", 1.7e-5, 0.025)],
        None,
    )
    .unwrap();
    let optical = [optical("A", 0.0003, 0.03, 1)];
    let cold = ThermochemicalRayleighState::try_new(
        source.evaluate(&point(300.0, 100_000.0)).unwrap(),
        &optical,
    )
    .unwrap();
    let hot = ThermochemicalRayleighState::try_new(
        source.evaluate(&point(800.0, 100_000.0)).unwrap(),
        &optical,
    )
    .unwrap();
    assert_eq!(
        cold.rayleigh().model_identity(),
        hot.rayleigh().model_identity()
    );
    assert_ne!(
        cold.rayleigh().content_identity(),
        hot.rayleigh().content_identity()
    );
    let cold_receipt = &cold
        .resolved()
        .thermodynamics()
        .receipt()
        .component_receipts()[0];
    let hot_receipt = &hot
        .resolved()
        .thermodynamics()
        .receipt()
        .component_receipts()[0];
    assert_eq!(cold_receipt.region_index(), 0);
    assert_eq!(hot_receipt.region_index(), 1);
    assert_ne!(cold_receipt.card_identity(), hot_receipt.card_identity());

    // A changed caloric reference can leave the ENTIRE GasState and optics
    // identical. Its original owner/receipt must still change the full identity.
    let shifted = crate::thermochemical_gas::tests::species_with_energy_reference(
        "A", 0.028, 3.5, 0.001, 25_000.0,
    );
    let shifted_source = ThermochemicalGasModel::new(
        FrozenIdealGasMixtureModelV1::new(
            vec![shifted],
            Composition::new(CompositionBasis::MoleFraction, vec![1.0]).unwrap(),
        )
        .unwrap(),
        vec![transport("A", 1.7e-5, 0.025)],
        None,
    )
    .unwrap();
    let shifted = ThermochemicalRayleighState::try_new(
        shifted_source.evaluate(&point(300.0, 100_000.0)).unwrap(),
        &optical,
    )
    .unwrap();
    assert_eq!(cold.rayleigh().gas(), shifted.rayleigh().gas());
    assert_eq!(
        cold.rayleigh().model_identity(),
        shifted.rayleigh().model_identity()
    );
    assert_ne!(
        cold.resolved().specific_enthalpy_j_kg(),
        shifted.resolved().specific_enthalpy_j_kg()
    );
    assert_ne!(
        cold.rayleigh().content_identity(),
        shifted.rayleigh().content_identity()
    );
}

#[test]
fn g0_species_join_refuses_missing_duplicate_and_foreign_optical_rows() {
    let accepted = model(vec![0.75, 0.25], 0.001)
        .evaluate(&point(600.0, 100_000.0))
        .unwrap();
    let [a, b] = optics();
    for rows in [
        vec![],
        vec![a.clone()],
        vec![a.clone(), a.clone()],
        vec![a.clone(), optical("foreign", 0.0006, 0.3, 2)],
        vec![a.clone(), b.clone(), b.clone()],
    ] {
        assert!(matches!(
            ThermochemicalRayleighState::try_new(accepted.clone(), &rows),
            Err(GasScatteringError::Components)
        ));
    }
    let canonical =
        ThermochemicalRayleighState::try_new(accepted.clone(), &[a.clone(), b.clone()]).unwrap();
    let permuted = ThermochemicalRayleighState::try_new(accepted.clone(), &[b, a.clone()]).unwrap();
    assert_eq!(canonical, permuted);
    let changed_source =
        ThermochemicalRayleighState::try_new(accepted, &[a, optical("B", 0.0006, 0.3, 3)]).unwrap();
    assert_eq!(
        canonical.rayleigh().mass_scattering_reference_m2_kg(),
        changed_source.rayleigh().mass_scattering_reference_m2_kg()
    );
    assert_ne!(
        canonical.rayleigh().model_identity(),
        changed_source.rayleigh().model_identity()
    );
}
