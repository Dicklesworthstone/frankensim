//! Synthetic original-source fixtures; no measured-gas accuracy claim.
use super::*;
use crate::gas_scattering::thermochemical::SpeciesMolecularRayleigh;
use fs_qty::semantic::CompositionBasis;
use fs_thermochem::SpeciesId;

fn optics() -> [MolecularRayleigh; 2] {
    [0.00015, 0.0006]
        .map(|r| MolecularRayleigh::try_new(r, 2.5e25, 550.0, 0.03, ContentHash([1; 32])).unwrap())
}

fn fixture() -> (
    MolecularRayleighMixture,
    RayleighGasState,
    MolecularRefractionModel,
) {
    let mut specs = [GasSpec::dry_air_ussa1976(); 2];
    specs[1].molar_mass *= 2.0;
    let original = MolecularRayleighMixture::try_new(
        &specs,
        Composition::new(CompositionBasis::MoleFraction, vec![0.75, 0.25]).unwrap(),
        &optics(),
        None,
    )
    .unwrap();
    let state = original.evaluate(300.0, 100_000.0).unwrap();
    let model = MolecularRefractionModel::from_mixture(&original, state.gas()).unwrap();
    (original, state, model)
}

fn near(a: f64, b: f64) {
    assert!(
        (a - b).abs() <= 5.0e-13 * b.abs().max(1.0e-300),
        "{a:.17e} != {b:.17e}"
    );
}

#[test]
fn g3_composition_model_retains_original_masses_and_partial_density_units() {
    let (original, state, model) = fixture();
    let n_total = state.gas().density / state.molar_mass_kg_mol() * AVOGADRO_PER_MOL;
    let numbers = [0.75 * n_total, 0.25 * n_total];
    let masses: [f64; 2] = core::array::from_fn(|i| {
        let component = &model.components()[i];
        assert_eq!(component.optical(), optics()[i]);
        assert_eq!(component.species(), None);
        near(component.molar_mass_kg_mol(), original.specs[i].molar_mass);
        numbers[i] / AVOGADRO_PER_MOL * component.molar_mass_kg_mol()
    });
    near(masses.iter().sum(), state.gas().density);
    for wavelength in [380.0, 550.0, 780.0] {
        let a = model
            .refractive_index(&numbers, PartialDensityBasis::MoleculesPerM3, wavelength)
            .unwrap();
        let b = model
            .refractive_index(&masses, PartialDensityBasis::MassKgPerM3, wavelength)
            .unwrap();
        let expected = state.refractive_index(wavelength).unwrap();
        near(a.phase_index, expected.phase_index);
        near(a.refractivity, expected.refractivity);
        near(b.refractivity, expected.refractivity);
        near(a.lorentz_lorenz, b.lorentz_lorenz);
    }
    let changed = GasState::try_new(&GasSpec::water_vapor_nist(), 300.0, 100_000.0).unwrap();
    assert_eq!(
        MolecularRefractionModel::from_mixture(&original, &changed),
        Err(GasScatteringError::StateMismatch)
    );
    let swapped = MolecularRayleighMixture::try_new(
        &original.specs,
        original.composition.clone(),
        &[optics()[1], optics()[0]],
        None,
    )
    .unwrap();
    let changed = MolecularRefractionModel::from_mixture(
        &swapped,
        swapped.evaluate(300.0, 100_000.0).unwrap().gas(),
    )
    .unwrap();
    assert_ne!(changed.content_identity(), model.content_identity());
}

#[test]
fn g3_thermochemical_model_retains_canonical_species_and_original_reference() {
    let source = crate::thermochemical_gas::tests::model(vec![0.75, 0.25], 0.001);
    let accepted = source
        .evaluate(&crate::thermochemical_gas::tests::point(700.0, 150_000.0))
        .unwrap();
    let rows: Vec<_> = ["B", "A"]
        .into_iter()
        .zip([optics()[1], optics()[0]])
        .map(|(name, optical)| SpeciesMolecularRayleigh {
            species: SpeciesId::new(name).unwrap(),
            optical,
        })
        .collect();
    let attached = ThermochemicalRayleighState::try_new(accepted, &rows).unwrap();
    let model = MolecularRefractionModel::from_thermochemical(&attached).unwrap();
    assert_eq!(model.components()[0].species(), Some("A"));
    assert_eq!(model.components()[1].species(), Some("B"));
    let n_total = attached.rayleigh().gas().density / attached.rayleigh().molar_mass_kg_mol()
        * AVOGADRO_PER_MOL;
    let result = model
        .refractive_index(
            &[0.75 * n_total, 0.25 * n_total],
            PartialDensityBasis::MoleculesPerM3,
            550.0,
        )
        .unwrap();
    near(
        result.refractivity,
        attached
            .rayleigh()
            .refractive_index(550.0)
            .unwrap()
            .refractivity,
    );
    for (component, receipt) in model.components().iter().zip(
        attached
            .resolved()
            .thermodynamics()
            .receipt()
            .component_receipts(),
    ) {
        assert_eq!(
            component.molar_mass_kg_mol().to_bits(),
            receipt.molar_mass_bits()
        );
    }
}

#[test]
fn g1_composition_changes_phase_at_constant_total_mass_and_retains_nonlinear_curvature() {
    let (_, _, model) = fixture();
    let basis = PartialDensityBasis::MassKgPerM3;
    let a = model.refractive_index(&[0.8, 0.2], basis, 550.0).unwrap();
    let b = model.refractive_index(&[0.2, 0.8], basis, 550.0).unwrap();
    assert!(b.phase_index > a.phase_index);
    for index in [a, b] {
        let l = index.lorentz_lorenz;
        let d = 1.0 - l;
        let n = ((1.0 + 2.0 * l) / d).sqrt();
        near(index.phase_index, n);
        near(index.derivative_lorentz_lorenz, 3.0 / (2.0 * n * d * d));
        near(
            index.second_derivative_lorentz_lorenz,
            3.0 * (1.0 + 8.0 * l) / (4.0 * n.powi(3) * d.powi(4)),
        );
    }
}

#[test]
fn g0_partial_density_vacuum_and_invalid_inputs_are_explicit() {
    let (_, _, model) = fixture();
    for basis in [
        PartialDensityBasis::MassKgPerM3,
        PartialDensityBasis::MoleculesPerM3,
    ] {
        let vacuum = model.refractive_index(&[0.0, 0.0], basis, 550.0).unwrap();
        assert_eq!(vacuum.phase_index, 1.0);
        assert_eq!(vacuum.refractivity, 0.0);
        for missing in [&[][..], &[1.0][..], &[1.0, 0.0, 0.0][..]] {
            assert_eq!(
                model.refractive_index(missing, basis, 550.0),
                Err(GasScatteringError::Components)
            );
        }
        for bad in [-1.0, f64::NAN, f64::INFINITY, f64::MAX] {
            assert_eq!(
                model.refractive_index(&[bad, 0.0], basis, 550.0),
                Err(GasScatteringError::RefractivityDomain)
            );
        }
        assert_eq!(
            model.refractive_index(&[1.0, 0.0], basis, 379.0),
            Err(GasScatteringError::Wavelength)
        );
        assert_eq!(
            model.refractive_index(&[f64::from_bits(1), 0.0], basis, 550.0),
            Err(GasScatteringError::Unrepresentable)
        );
    }
}
