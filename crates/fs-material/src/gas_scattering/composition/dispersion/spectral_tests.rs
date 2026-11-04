//! Explicit synthetic broadband sources; no measured-gas accuracy is asserted.
use super::*;
use fs_qty::semantic::CompositionBasis;

fn fixture() -> (MolecularRefractionModel, [MolecularRayleigh; 2]) {
    let mut specs = [GasSpec::dry_air_ussa1976(); 2];
    specs[1].molar_mass *= 2.0;
    let optics = core::array::from_fn(|i| {
        MolecularRayleigh::try_new(
            [0.00015, 0.0006][i],
            2.5e25,
            550.0,
            [0.0, 0.08][i],
            ContentHash([i as u8 + 1; 32]),
        )
        .unwrap()
    });
    let original = MolecularRayleighMixture::try_new(
        &specs,
        Composition::new(CompositionBasis::MoleFraction, vec![0.75, 0.25]).unwrap(),
        &optics,
        None,
    )
    .unwrap();
    let model = MolecularRefractionModel::from_mixture(
        &original,
        original.evaluate(300.0, 100_000.0).unwrap().gas(),
    )
    .unwrap();
    (model, optics)
}

fn shapes(optics: [MolecularRayleigh; 2], zero: bool) -> [MolecularCauchyDispersion; 2] {
    core::array::from_fn(|i| {
        MolecularCauchyDispersion::try_new_in_band(
            optics[i],
            if zero { 0.0 } else { [2e-6, 8e-6][i] },
            if zero { 0.0 } else { [1e-7, 3e-7][i] },
            [300.0, 2500.0],
            ContentHash([i as u8 + 3; 32]),
        )
        .unwrap()
    })
}

fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:.17e} != {b:.17e}");
}

#[test]
fn broadband_pure_source_and_mixture_group_delay_retain_the_original_equations() {
    let (original, optics) = fixture();
    let model = original.with_dispersion(&shapes(optics, false)).unwrap();
    assert_eq!(model.wavelength_range_nm(), [300.0, 2500.0]);
    let basis = PartialDensityBasis::MoleculesPerM3;
    for lambda in [300.0_f64, 550.0, 1550.0, 2500.0] {
        let x = 1.0 / (lambda * 1e-3).powi(2);
        let xr = 1.0 / 0.55_f64.powi(2);
        let r = 0.00015 + 2e-6 * (x - xr) + 1e-7 * (x * x - xr * xr);
        let pure = model
            .refractive_index(&[2.5e25, 0.0], basis, lambda)
            .unwrap();
        close(pure.refractivity, r, 3e-18);
        close(pure.group_refractivity, r + 4e-6 * x + 4e-7 * x * x, 3e-18);
        let numbers = [1.7e25, 0.6e25];
        let mixture = model.refractive_index(&numbers, basis, lambda).unwrap();
        let masses = core::array::from_fn::<_, 2, _>(|i| {
            numbers[i] * model.components()[i].molar_mass_kg_mol() / AVOGADRO_PER_MOL
        });
        let mass = model
            .refractive_index(&masses, PartialDensityBasis::MassKgPerM3, lambda)
            .unwrap();
        close(mixture.group_refractivity, mass.group_refractivity, 3e-18);
        if lambda > 300.0 && lambda < 2500.0 {
            let h = lambda * 1e-4;
            let lo = model.refractive_index(&numbers, basis, lambda - h).unwrap();
            let hi = model.refractive_index(&numbers, basis, lambda + h).unwrap();
            close(
                mixture.group_refractivity,
                mixture.refractivity - lambda * (hi.refractivity - lo.refractivity) / (2.0 * h),
                3e-12,
            );
        }
        assert!(mixture.group_index > mixture.phase_index);
    }
}

#[test]
fn broadband_molecular_scattering_uses_each_species_spectral_volume() {
    let (original, optics) = fixture();
    let numbers = [1.7e25, 0.6e25];
    for zero in [false, true] {
        let laws = shapes(optics, zero);
        let model = original.with_dispersion(&laws).unwrap();
        for lambda in [300.0_f64, 550.0, 1550.0, 2500.0] {
            let sample = model
                .at_wavelength(lambda, PartialDensityBasis::MoleculesPerM3)
                .unwrap()
                .scattering(&numbers)
                .unwrap();
            let mut total = 0.0;
            let mut isotropic = 0.0;
            for i in 0..2 {
                let [b, c] = laws[i].coefficients();
                let x = 1.0 / (lambda * 1e-3).powi(2);
                let xr = 1.0 / 0.55_f64.powi(2);
                let r = [0.00015, 0.0006][i] + b * (x - xr) + c * (x * x - xr * xr);
                let n2 = (1.0 + r).powi(2);
                let volume = ((n2 - 1.0) / (n2 + 2.0)) / 2.5e25;
                let delta = optics[i].depolarization_ratio();
                let sigma = 24.0
                    * core::f64::consts::PI.powi(3)
                    * ((6.0 + 3.0 * delta) / (6.0 - 7.0 * delta))
                    * volume.powi(2)
                    / (lambda * 1e-9).powi(4);
                total += numbers[i] * sigma;
                isotropic += numbers[i] * sigma * 3.0 * delta / (2.0 + delta);
            }
            // The independent n² subtraction loses dilute precision; retain a
            // strict relative tolerance consistent with that reference's units.
            close(sample.scattering_per_m, total, total * 3e-12);
            let w = isotropic / total;
            close(sample.depolarization_ratio, 2.0 * w / (3.0 - w), 3e-13);
            if zero && (380.0..=780.0).contains(&lambda) {
                assert_eq!(
                    sample,
                    original
                        .at_wavelength(lambda, PartialDensityBasis::MoleculesPerM3)
                        .unwrap()
                        .scattering(&numbers)
                        .unwrap()
                );
            }
        }
    }
}

#[test]
fn explicit_bands_preserve_legacy_sources_and_refuse_extrapolation_or_invalid_fits() {
    let (original, optics) = fixture();
    let broad = shapes(optics, false);
    let model = original.with_dispersion(&broad).unwrap();
    assert!(
        original
            .at_wavelength(1550.0, PartialDensityBasis::MassKgPerM3)
            .is_err()
    );
    assert!(
        MolecularCauchyDispersion::try_new(
            optics[0],
            2e-6,
            1e-7,
            [300.0, 2500.0],
            ContentHash([3; 32])
        )
        .is_err()
    );
    for bad in [
        300.0_f64.next_down(),
        2500.0_f64.next_up(),
        0.0,
        f64::NAN,
        f64::INFINITY,
    ] {
        assert!(
            model
                .at_wavelength(bad, PartialDensityBasis::MassKgPerM3)
                .is_err()
        );
    }
    for band in [
        [0.0, 2500.0],
        [300.0, f64::INFINITY],
        [600.0, 2500.0],
        [550.0, 550.0],
    ] {
        assert!(
            MolecularCauchyDispersion::try_new_in_band(
                optics[0],
                2e-6,
                1e-7,
                band,
                ContentHash([3; 32])
            )
            .is_err()
        );
    }
    // A distant source endpoint with negative refractivity must refuse even
    // though the anchor itself is admissible.
    assert!(
        MolecularCauchyDispersion::try_new_in_band(
            optics[0],
            1e-3,
            0.0,
            [550.0, 2500.0],
            ContentHash([3; 32])
        )
        .is_err()
    );
    let mut mixed = broad;
    mixed[1] = MolecularCauchyDispersion::try_new(
        optics[1],
        8e-6,
        3e-7,
        [400.0, 780.0],
        ContentHash([4; 32]),
    )
    .unwrap();
    let narrowed = original.with_dispersion(&mixed).unwrap();
    assert_eq!(narrowed.wavelength_range_nm(), [400.0, 780.0]);
    assert!(
        narrowed
            .at_wavelength(1550.0, PartialDensityBasis::MassKgPerM3)
            .is_err()
    );
    assert_ne!(narrowed.content_identity(), model.content_identity());
    let visible = MolecularCauchyDispersion::try_new(
        optics[0],
        2e-6,
        1e-7,
        [380.0, 780.0],
        ContentHash([3; 32]),
    )
    .unwrap();
    assert_eq!(
        visible,
        MolecularCauchyDispersion::try_new_in_band(
            optics[0],
            2e-6,
            1e-7,
            [380.0, 780.0],
            ContentHash([3; 32])
        )
        .unwrap()
    );
}
