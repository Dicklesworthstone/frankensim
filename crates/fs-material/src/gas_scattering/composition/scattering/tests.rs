//! Synthetic optical sources exercise the law, not a measured-gas error bound.
use super::*;
use fs_qty::semantic::CompositionBasis;

fn fixture(
    refractivities: [f64; 2],
) -> (
    MolecularRefractionModel,
    RayleighGasState,
    [MolecularRayleigh; 2],
) {
    let mut specs = [GasSpec::dry_air_ussa1976(); 2];
    specs[1].molar_mass *= 2.0;
    let optical = core::array::from_fn(|i| {
        MolecularRayleigh::try_new(
            refractivities[i],
            2.5e25,
            550.0,
            [0.0, 0.3][i],
            ContentHash([i as u8 + 1; 32]),
        )
        .unwrap()
    });
    let original = MolecularRayleighMixture::try_new(
        &specs,
        Composition::new(CompositionBasis::MoleFraction, vec![0.75, 0.25]).unwrap(),
        &optical,
        None,
    )
    .unwrap();
    let state = original.evaluate(300.0, 100_000.0).unwrap();
    let model = MolecularRefractionModel::from_mixture(&original, state.gas()).unwrap();
    (model, state, optical)
}

fn dispersed(
    model: &MolecularRefractionModel,
    optical: [MolecularRayleigh; 2],
) -> MolecularRefractionModel {
    let shapes: [MolecularCauchyDispersion; 2] = core::array::from_fn(|i| {
        MolecularCauchyDispersion::try_new(
            optical[i],
            [0.0, 2.0e-5][i],
            [0.0, 1.0e-6][i],
            [380.0, 780.0],
            ContentHash([i as u8 + 3; 32]),
        )
        .unwrap()
    });
    model.with_dispersion(&shapes).unwrap()
}

fn near(a: f64, b: f64) {
    assert!(
        (a - b).abs() <= 5.0e-13 * a.abs().max(b.abs()),
        "{a:.17e} != {b:.17e}"
    );
}

// Independent electric-dipole form: sigma = (8 pi/3) k^4 alpha_volume^2 F.
fn dipole_cross_section(r: f64, delta: f64, wavelength_nm: f64) -> f64 {
    let pi = core::f64::consts::PI;
    let alpha = 3.0 / (4.0 * pi * 2.5e25) * (r * (2.0 + r) / (3.0 + r * (2.0 + r)));
    let k = 2.0 * pi / (wavelength_nm * 1.0e-9);
    (8.0 * pi / 3.0) * k.powi(4) * alpha * alpha * (6.0 + 3.0 * delta) / (6.0 - 7.0 * delta)
}

fn phase(delta: f64, cosine: f64) -> f64 {
    3.0 * ((1.0 + delta) + (1.0 - delta) * cosine * cosine)
        / (8.0 * core::f64::consts::PI * (2.0 + delta))
}

#[test]
fn g1_dispersive_cross_sections_and_full_angular_mixture_match_independent_molecules() {
    let (reference, _, optical) = fixture([0.00015, 0.0006]);
    let model = dispersed(&reference, optical);
    let numbers = [2.0e25, 1.0e24];
    let masses: [f64; 2] = core::array::from_fn(|i| {
        numbers[i] / AVOGADRO_PER_MOL * model.components()[i].molar_mass_kg_mol()
    });
    let mut deltas = [0.0; 3];
    for (i, lambda) in [380.0, 550.0, 780.0].into_iter().enumerate() {
        let prepared = model
            .at_wavelength(lambda, PartialDensityBasis::MoleculesPerM3)
            .unwrap();
        let sample = prepared.scattering(&numbers).unwrap();
        let x = 1.0 / (lambda * 1.0e-3).powi(2);
        let xr = 1.0 / 0.55_f64.powi(2);
        let r = 0.0006 + 2.0e-5 * (x - xr) + 1.0e-6 * (x * x - xr * xr);
        let sigma = [
            dipole_cross_section(0.00015, 0.0, lambda),
            dipole_cross_section(r, 0.3, lambda),
        ];
        let beta = [numbers[0] * sigma[0], numbers[1] * sigma[1]];
        near(sample.scattering_per_m, beta[0] + beta[1]);
        // Pure species and a mixture share the SAME prepared spectral owner.
        near(
            prepared
                .scattering(&[0.0, numbers[1]])
                .unwrap()
                .scattering_per_m,
            beta[1],
        );
        // Check the entire angular law, independently of the w-to-delta reduction.
        for cosine in [-1.0, -0.37, 0.0, 0.42, 1.0] {
            near(
                sample.scattering_per_m * phase(sample.depolarization_ratio, cosine),
                beta[0] * phase(0.0, cosine) + beta[1] * phase(0.3, cosine),
            );
        }
        let mass_sample = model
            .at_wavelength(lambda, PartialDensityBasis::MassKgPerM3)
            .unwrap()
            .scattering(&masses)
            .unwrap();
        near(mass_sample.scattering_per_m, sample.scattering_per_m);
        near(
            mass_sample.depolarization_ratio,
            sample.depolarization_ratio,
        );
        near(mass_sample.mass_density_kg_m3, sample.mass_density_kg_m3);
        assert_eq!(mass_sample.mass_density_kg_m3, masses[0] + masses[1]);
        deltas[i] = sample.depolarization_ratio;
    }
    // Relative polarizabilities change: a gray mixture delta would be wrong.
    assert!(deltas[0] > deltas[1] && deltas[1] > deltas[2]);
    assert!(deltas[0] - deltas[2] > 0.03);
}

#[test]
fn g1_constant_and_zero_shape_paths_retain_original_molecular_and_mixture_laws() {
    let (reference, state, optical) = fixture([0.00015, 0.0006]);
    let zero = optical.map(|row| {
        MolecularCauchyDispersion::try_new(row, 0.0, 0.0, [380.0, 780.0], ContentHash([4; 32]))
            .unwrap()
    });
    let zero = reference.with_dispersion(&zero).unwrap();
    let total_number = state.gas().density / state.molar_mass_kg_mol() * AVOGADRO_PER_MOL;
    let numbers = [0.75 * total_number, 0.25 * total_number];
    let masses: [f64; 2] = core::array::from_fn(|i| {
        numbers[i] / AVOGADRO_PER_MOL * reference.components()[i].molar_mass_kg_mol()
    });
    for lambda in [380.0, 550.0, 780.0] {
        for (basis, densities) in [
            (PartialDensityBasis::MoleculesPerM3, numbers),
            (PartialDensityBasis::MassKgPerM3, masses),
        ] {
            let original = reference
                .at_wavelength(lambda, basis)
                .unwrap()
                .scattering(&densities)
                .unwrap();
            let constant = zero
                .at_wavelength(lambda, basis)
                .unwrap()
                .scattering(&densities)
                .unwrap();
            assert_eq!(constant, original);
            near(
                original.scattering_per_m,
                state.scattering_per_m(lambda).unwrap(),
            );
            near(original.mass_density_kg_m3, state.gas().density);
            near(original.depolarization_ratio, state.depolarization_ratio());
        }
        let prepared = reference
            .at_wavelength(lambda, PartialDensityBasis::MoleculesPerM3)
            .unwrap();
        for i in 0..2 {
            let mut pure = [0.0; 2];
            pure[i] = 1.0;
            assert_eq!(
                prepared
                    .scattering(&pure)
                    .unwrap()
                    .scattering_per_m
                    .to_bits(),
                optical[i].cross_section_m2(lambda).unwrap().to_bits()
            );
        }
    }
}

#[test]
fn g2_scattering_shares_full_admission_and_retains_vacuum_and_inert_mass() {
    let (reference, _, _) = fixture([0.00015, 0.0006]);
    let mass = reference
        .at_wavelength(550.0, PartialDensityBasis::MassKgPerM3)
        .unwrap();
    assert_eq!(
        mass.scattering(&[0.0; 2]).unwrap(),
        CompositionRayleighScattering {
            mass_density_kg_m3: 0.0,
            scattering_per_m: 0.0,
            depolarization_ratio: 0.0,
        }
    );
    for invalid in [&[][..], &[1.0][..], &[0.0; 3][..]] {
        assert_eq!(
            mass.scattering(invalid),
            Err(GasScatteringError::Components)
        );
    }
    for bad in [-1.0, f64::NAN, f64::INFINITY, 1_000.0] {
        assert_eq!(
            mass.scattering(&[bad, 0.0]),
            Err(GasScatteringError::RefractivityDomain)
        );
    }
    assert_eq!(
        mass.scattering(&[f64::from_bits(1), 0.0]),
        Err(GasScatteringError::Unrepresentable)
    );
    let number = reference
        .at_wavelength(550.0, PartialDensityBasis::MoleculesPerM3)
        .unwrap();
    // Refraction and physical mass remain representable, but a positive photon
    // collision rate cannot silently underflow into transparent vacuum.
    assert!(number.refractive_index(&[1.0e-293, 0.0]).is_ok());
    assert_eq!(
        number.scattering(&[1.0e-293, 0.0]),
        Err(GasScatteringError::Unrepresentable)
    );
    let (inert, _, _) = fixture([0.0; 2]);
    let sample = inert
        .at_wavelength(550.0, PartialDensityBasis::MassKgPerM3)
        .unwrap()
        .scattering(&[0.4, 0.7])
        .unwrap();
    assert_eq!(sample.mass_density_kg_m3, 0.4 + 0.7);
    assert_eq!(sample.scattering_per_m, 0.0);
    assert_eq!(sample.depolarization_ratio, 0.0);
}

#[test]
fn g2_componentwise_density_envelopes_bound_every_spectral_collision_rate() {
    let (reference, _, optical) = fixture([0.00015, 0.0006]);
    let model = dispersed(&reference, optical);
    for lambda in [380.0, 550.0, 780.0] {
        let prepared = model
            .at_wavelength(lambda, PartialDensityBasis::MassKgPerM3)
            .unwrap();
        let bound = prepared.scattering(&[2.0, 2.0]).unwrap();
        for local in [[0.0, 0.0], [2.0, 0.0], [0.0, 2.0], [0.6, 1.7], [2.0, 2.0]] {
            assert!(
                prepared.scattering(&local).unwrap().scattering_per_m <= bound.scattering_per_m
            );
        }
        let half = prepared.scattering(&[1.0, 1.0]).unwrap();
        near(half.scattering_per_m, 0.5 * bound.scattering_per_m);
        assert_eq!(half.depolarization_ratio, bound.depolarization_ratio);
    }
}
