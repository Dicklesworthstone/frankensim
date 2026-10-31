use super::*;
use fs_qty::semantic::CompositionBasis;

fn near(got: f64, expected: f64) {
    assert!(
        (got - expected).abs() <= 5.0e-13 * expected.abs().max(1.0e-300),
        "{got:.17e} != {expected:.17e}"
    );
}

fn optical(refractivity: f64, delta: f64, source: u8) -> MolecularRayleigh {
    MolecularRayleigh::try_new(
        refractivity,
        2.5e25,
        550.0,
        delta,
        ContentHash([source; 32]),
    )
    .unwrap()
}

fn model(fractions: &[f64], optics: &[MolecularRayleigh]) -> MolecularRayleighMixture {
    let mut specs = vec![GasSpec::dry_air_ussa1976(); fractions.len()];
    for (i, spec) in specs.iter_mut().enumerate() {
        spec.molar_mass *= (i + 1) as f64;
    }
    MolecularRayleighMixture::try_new(
        &specs,
        Composition::new(CompositionBasis::MoleFraction, fractions.to_vec()).unwrap(),
        optics,
        None,
    )
    .unwrap()
}

#[test]
fn g1_reference_cross_section_and_inverse_fourth_law() {
    let molecule = optical(0.0003, 0.03, 1);
    // Independent 60-decimal evaluation via electric dipole polarizability:
    // alpha_volume=3(n²-1)/(4pi N(n²+2)), sigma=(8pi/3)k⁴ alpha_volume² F.
    // This is a synthetic reference, not a measured dry-air data assertion.
    let expected = 5.473_755_023_921_817_5e-31;
    near(molecule.cross_section_m2(550.0).unwrap(), expected);
    for wavelength in [380.0, 400.0, 550.0, 700.0, 780.0] {
        near(
            molecule.cross_section_m2(wavelength).unwrap()
                / molecule.cross_section_m2(550.0).unwrap(),
            (550.0 / wavelength).powi(4),
        );
    }
    let zero_delta = optical(0.0003, 0.0, 1);
    near(
        molecule.cross_section_m2(550.0).unwrap() / zero_delta.cross_section_m2(550.0).unwrap(),
        (6.0 + 3.0 * 0.03) / (6.0 - 7.0 * 0.03),
    );
}

#[test]
fn g3_pressure_temperature_follow_original_gas_density() {
    let model = model(&[1.0], &[optical(0.0003, 0.03, 1)]);
    let a = model.evaluate(300.0, 100_000.0).unwrap();
    let pressure = model.evaluate(300.0, 200_000.0).unwrap();
    let temperature = model.evaluate(600.0, 100_000.0).unwrap();
    for wavelength in [380.0, 550.0, 780.0] {
        near(
            pressure.scattering_per_m(wavelength).unwrap(),
            2.0 * a.scattering_per_m(wavelength).unwrap(),
        );
        near(
            temperature.scattering_per_m(wavelength).unwrap(),
            0.5 * a.scattering_per_m(wavelength).unwrap(),
        );
        near(
            a.scattering_per_m(wavelength).unwrap(),
            a.gas().density * a.mass_scattering_m2_kg(wavelength).unwrap(),
        );
    }
    assert_eq!(
        a.mass_scattering_reference_m2_kg(),
        pressure.mass_scattering_reference_m2_kg()
    );
    assert_ne!(a.content_identity(), pressure.content_identity());
    near(
        a.scattering_per_m(550.0).unwrap(),
        a.gas().density / a.molar_mass_kg_mol()
            * AVOGADRO_PER_MOL
            * a.molecular_cross_section_reference_m2(),
    );
}

#[test]
fn g1_mixture_adds_number_weighted_differential_cross_sections() {
    let optics = [optical(0.00015, 0.0, 1), optical(0.0006, 0.3, 2)];
    let mixture = model(&[0.75, 0.25], &optics)
        .evaluate(300.0, 100_000.0)
        .unwrap();
    let pure_a = model(&[1.0], &[optics[0]])
        .evaluate(300.0, 100_000.0)
        .unwrap();
    let pure_b = model(&[1.0], &[optics[1]])
        .evaluate(300.0, 100_000.0)
        .unwrap();
    let sigma_a = optics[0].cross_section_m2(550.0).unwrap();
    let sigma_b = optics[1].cross_section_m2(550.0).unwrap();
    near(
        mixture.molecular_cross_section_reference_m2(),
        0.75 * sigma_a + 0.25 * sigma_b,
    );
    for cosine in [-1.0, -0.6, 0.0, 0.4, 1.0] {
        let want = 0.75 * sigma_a * pure_a.phase_per_sr(cosine).unwrap()
            + 0.25 * sigma_b * pure_b.phase_per_sr(cosine).unwrap();
        near(
            mixture.molecular_cross_section_reference_m2() * mixture.phase_per_sr(cosine).unwrap(),
            want,
        );
    }
    assert!(
        (mixture.depolarization_ratio() - 0.25 * 0.3).abs() > 0.1,
        "strongly scattering species must dominate the angular law"
    );
}

#[test]
fn g0_molecular_phase_is_normalized_and_reciprocal() {
    for delta in [0.0, 0.03, 0.4, 0.85] {
        let state = model(&[1.0], &[optical(0.0003, delta, 1)])
            .evaluate(300.0, 100_000.0)
            .unwrap();
        // Simpson integrates this quadratic phase exactly.
        let total = core::f64::consts::TAU / 3.0
            * (state.phase_per_sr(-1.0).unwrap()
                + 4.0 * state.phase_per_sr(0.0).unwrap()
                + state.phase_per_sr(1.0).unwrap());
        near(total, 1.0);
        for cosine in [0.0, 0.2, 0.9, 1.0] {
            near(
                state.phase_per_sr(cosine).unwrap(),
                state.phase_per_sr(-cosine).unwrap(),
            );
        }
    }
}

#[test]
fn g3_mass_mole_basis_equivalence_and_original_state_binding() {
    let optics = [optical(0.0002, 0.0, 1), optical(0.0005, 0.17, 2)];
    let mole_model = model(&[0.65, 0.35], &optics);
    let masses: Vec<_> = mole_model
        .specs
        .iter()
        .map(|spec| MolarMass::new(spec.molar_mass))
        .collect();
    let mass = mole_model.composition.to_mass_fractions(&masses).unwrap();
    let mass_model =
        MolecularRayleighMixture::try_new(&mole_model.specs, mass, &optics, None).unwrap();
    let a = mole_model.evaluate(300.0, 100_000.0).unwrap();
    let b = mass_model.evaluate(300.0, 100_000.0).unwrap();
    near(
        a.scattering_per_m(550.0).unwrap(),
        b.scattering_per_m(550.0).unwrap(),
    );
    near(a.depolarization_ratio(), b.depolarization_ratio());
    let original = GasState::try_new_mixture(
        &mole_model.specs,
        &mole_model.composition,
        300.0,
        100_000.0,
        None,
    )
    .unwrap();
    assert_eq!(a.gas(), &original);
    assert_eq!(mole_model.bind_gas_state(&original).unwrap(), a);
    // A different but internally coherent gas cannot borrow these optics.
    let changed = GasState::try_new(&GasSpec::water_vapor_nist(), 300.0, 100_000.0).unwrap();
    assert_eq!(
        mole_model.bind_gas_state(&changed),
        Err(GasScatteringError::StateMismatch)
    );
}

#[test]
fn g0_transparent_component_and_optical_source_identity_are_retained() {
    let optics = [optical(0.0, 0.0, 1), optical(0.0003, 0.03, 2)];
    let mixed = model(&[0.5, 0.5], &optics)
        .evaluate(300.0, 100_000.0)
        .unwrap();
    near(
        mixed.molecular_cross_section_reference_m2(),
        0.5 * optics[1].cross_section_m2(550.0).unwrap(),
    );
    near(mixed.depolarization_ratio(), 0.03);
    let transparent = model(&[1.0], &optics[..1])
        .evaluate(300.0, 100_000.0)
        .unwrap();
    assert_eq!(transparent.scattering_per_m(380.0).unwrap(), 0.0);
    let replaced_source = model(&[1.0], &[optical(0.0, 0.0, 3)])
        .evaluate(300.0, 100_000.0)
        .unwrap();
    assert_ne!(
        transparent.content_identity(),
        replaced_source.content_identity()
    );
}

#[test]
fn g0_invalid_sources_ranges_and_composition_refuse() {
    for delta in [-0.01, 6.0 / 7.0, 1.0, f64::NAN] {
        assert!(
            MolecularRayleigh::try_new(0.0003, 2.5e25, 550.0, delta, ContentHash([1; 32])).is_err()
        );
    }
    for density in [0.0, -1.0, f64::NAN, f64::INFINITY, 1.0e-300, 1.0e300] {
        assert!(
            MolecularRayleigh::try_new(0.0003, density, 550.0, 0.0, ContentHash([1; 32])).is_err()
        );
    }
    assert!(MolecularRayleigh::try_new(0.0003, 2.5e25, 550.0, 0.0, ContentHash([0; 32])).is_err());
    let molecule = optical(0.0003, 0.0, 1);
    for wavelength in [0.0, 379.0, 781.0, f64::NAN] {
        assert_eq!(
            molecule.cross_section_m2(wavelength),
            Err(GasScatteringError::Wavelength)
        );
    }
    let volume_composition = Composition::new(CompositionBasis::VolumeFraction, vec![1.0]).unwrap();
    let invalid = MolecularRayleighMixture::try_new(
        &[GasSpec::dry_air_ussa1976()],
        volume_composition,
        &[molecule],
        None,
    )
    .unwrap();
    assert!(invalid.evaluate(300.0, 100_000.0).is_err());
    let valid = model(&[1.0], &[molecule]);
    assert!(valid.evaluate(f64::NAN, 100_000.0).is_err());
    assert!(valid.evaluate(300.0, 0.0).is_err());
}
