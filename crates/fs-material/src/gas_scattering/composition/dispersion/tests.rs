//! Synthetic source fits, not measured constants or a gas accuracy claim.
use super::*;
use fs_qty::semantic::CompositionBasis;

fn fixture() -> (MolecularRefractionModel, [MolecularCauchyDispersion; 2]) {
    let mut specs = [GasSpec::dry_air_ussa1976(); 2];
    specs[1].molar_mass *= 2.0;
    let optical = [0.00015, 0.0006]
        .map(|r| MolecularRayleigh::try_new(r, 2.5e25, 550.0, 0.03, ContentHash([1; 32])).unwrap());
    let original = MolecularRayleighMixture::try_new(
        &specs,
        Composition::new(CompositionBasis::MoleFraction, vec![0.75, 0.25]).unwrap(),
        &optical,
        None,
    )
    .unwrap();
    let model = MolecularRefractionModel::from_mixture(
        &original,
        original.evaluate(300.0, 100_000.0).unwrap().gas(),
    )
    .unwrap();
    let shapes = core::array::from_fn(|i| {
        MolecularCauchyDispersion::try_new(
            optical[i],
            [2.0e-6, 8.0e-6][i],
            [1.0e-7, 3.0e-7][i],
            [380.0, 780.0],
            ContentHash([i as u8 + 2; 32]),
        )
        .unwrap()
    });
    (model, shapes)
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:.17e} != {b:.17e}");
}

#[test]
fn g1_zero_dispersion_retains_frozen_bits_and_nonzero_dispersion_retains_reference() {
    let (reference, shapes) = fixture();
    let zero = shapes.map(|shape| {
        MolecularCauchyDispersion::try_new(
            shape.reference(),
            0.0,
            0.0,
            [380.0, 780.0],
            shape.source_identity(),
        )
        .unwrap()
    });
    let constant = reference.with_dispersion(&zero).unwrap();
    let spectral = reference.with_dispersion(&shapes).unwrap();
    for basis in [
        PartialDensityBasis::MassKgPerM3,
        PartialDensityBasis::MoleculesPerM3,
    ] {
        let densities = if basis == PartialDensityBasis::MassKgPerM3 {
            [1.3, 0.4]
        } else {
            [2.0e25, 1.0e25]
        };
        for lambda in [380.0, 550.0, 780.0] {
            assert_eq!(
                constant
                    .refractive_index(&densities, basis, lambda)
                    .unwrap(),
                reference
                    .refractive_index(&densities, basis, lambda)
                    .unwrap()
            );
        }
        let original = reference
            .refractive_index(&densities, basis, 550.0)
            .unwrap();
        let anchored = spectral.refractive_index(&densities, basis, 550.0).unwrap();
        assert_eq!(
            original.phase_index.to_bits(),
            anchored.phase_index.to_bits()
        );
        assert_eq!(
            original.refractivity.to_bits(),
            anchored.refractivity.to_bits()
        );
        assert!(anchored.group_index > anchored.phase_index);
    }
}

#[test]
fn g1_spectral_group_derivative_obeys_the_full_nonlinear_mixture_and_both_density_bases() {
    let (reference, shapes) = fixture();
    let model = reference.with_dispersion(&shapes).unwrap();
    let masses = [1.3, 0.4];
    let numbers: [f64; 2] = core::array::from_fn(|i| {
        masses[i] / model.components()[i].molar_mass_kg_mol() * AVOGADRO_PER_MOL
    });
    for lambda in [400.0, 550.0, 760.0] {
        let mass = model
            .at_wavelength(lambda, PartialDensityBasis::MassKgPerM3)
            .unwrap();
        assert_eq!(mass.coefficients().len(), 2);
        let index = mass.refractive_index(&masses).unwrap();
        let number = model
            .refractive_index(&numbers, PartialDensityBasis::MoleculesPerM3, lambda)
            .unwrap();
        near(index.refractivity, number.refractivity, 2.0e-18);
        near(index.group_refractivity, number.group_refractivity, 2.0e-18);
        // Independently differentiate the full index, not constituent groups.
        let h = 0.005;
        let lo = model
            .refractive_index(&masses, PartialDensityBasis::MassKgPerM3, lambda - h)
            .unwrap();
        let hi = model
            .refractive_index(&masses, PartialDensityBasis::MassKgPerM3, lambda + h)
            .unwrap();
        let expected =
            index.refractivity - lambda * (hi.refractivity - lo.refractivity) / (2.0 * h);
        near(index.group_refractivity, expected, 2.0e-13);
        near(index.group_index, 1.0 + expected, 2.0e-13);
        assert!(index.group_index > index.phase_index);
        // At the original pure number density, LL conversion inverts exactly.
        let pure = model
            .refractive_index(&[2.5e25, 0.0], PartialDensityBasis::MoleculesPerM3, lambda)
            .unwrap();
        let x = 1.0 / (lambda * 1.0e-3).powi(2);
        let xr = 1.0 / 0.55_f64.powi(2);
        let r = 0.00015 + 2.0e-6 * (x - xr) + 1.0e-7 * (x * x - xr * xr);
        near(pure.refractivity, r, 2.0e-18);
        near(
            pure.group_refractivity,
            r + 4.0e-6 * x + 4.0e-7 * x * x,
            2.0e-18,
        );
    }
    let vacuum = model
        .refractive_index(&[0.0; 2], PartialDensityBasis::MassKgPerM3, 550.0)
        .unwrap();
    assert_eq!(vacuum.phase_index, 1.0);
    assert_eq!(vacuum.group_index, 1.0);
    assert_eq!(vacuum.group_refractivity, 0.0);
}

#[test]
fn g1_spectral_curvature_recovers_cauchy_and_differentiates_the_original_mixture() {
    let (reference, shapes) = fixture();
    let model = reference.with_dispersion(&shapes).unwrap();
    let numbers = [7.0e25, 2.0e25];
    let masses: [f64; 2] = core::array::from_fn(|i| {
        numbers[i] * model.components()[i].molar_mass_kg_mol() / AVOGADRO_PER_MOL
    });
    for lambda in [400.0, 550.0, 760.0] {
        let spectral = model
            .at_wavelength(lambda, PartialDensityBasis::MoleculesPerM3)
            .unwrap();
        // At its original pure number density, Lorentz--Lorenz must invert
        // exactly to the source Cauchy law, including its ordinary curvature.
        let pure = spectral.refractive_index_jet(&[2.5e25, 0.0]).unwrap();
        let x = (lambda * 1e-3).powi(-2);
        near(
            pure.scaled_wavelength_derivatives[0],
            -4e-6 * x - 4e-7 * x * x,
            2e-18,
        );
        near(
            pure.scaled_wavelength_derivatives[1],
            12e-6 * x + 20e-7 * x * x,
            2e-18,
        );
        let number = spectral.refractive_index_jet(&numbers).unwrap();
        let mass = model
            .at_wavelength(lambda, PartialDensityBasis::MassKgPerM3)
            .unwrap()
            .refractive_index_jet(&masses)
            .unwrap();
        for i in 0..2 {
            near(
                number.scaled_wavelength_derivatives[i],
                mass.scaled_wavelength_derivatives[i],
                3e-18,
            );
        }
        // Independent five-point differentiation of phase REFRACTIVITY keeps
        // subtraction of the unit vacuum index out of the curvature oracle.
        let h = 0.2;
        let r = |offset: f64| {
            model
                .refractive_index(
                    &numbers,
                    PartialDensityBasis::MoleculesPerM3,
                    lambda + offset * h,
                )
                .unwrap()
                .refractivity
        };
        let slope = (r(-2.0) - 8.0 * r(-1.0) + 8.0 * r(1.0) - r(2.0)) / (12.0 * h);
        let curvature =
            (-r(-2.0) + 16.0 * r(-1.0) - 30.0 * r(0.0) + 16.0 * r(1.0) - r(2.0)) / (12.0 * h * h);
        near(
            number.scaled_wavelength_derivatives[0],
            lambda * slope,
            5e-14,
        );
        near(
            number.scaled_wavelength_derivatives[1],
            lambda * lambda * curvature,
            3e-11,
        );
        assert_eq!(number.index, spectral.refractive_index(&numbers).unwrap());
    }
}

#[test]
fn g0_spectral_curvature_retains_vacuum_frozen_and_trace_composition_limits() {
    let (reference, shapes) = fixture();
    let model = reference.with_dispersion(&shapes).unwrap();
    let constant = reference
        .at_wavelength(550.0, PartialDensityBasis::MoleculesPerM3)
        .unwrap();
    assert_eq!(
        constant
            .refractive_index_jet(&[2e25, 1e25])
            .unwrap()
            .scaled_wavelength_derivatives,
        [0.0; 2]
    );
    let spectral = model
        .at_wavelength(550.0, PartialDensityBasis::MoleculesPerM3)
        .unwrap();
    let vacuum = spectral.refractive_index_jet(&[0.0; 2]).unwrap();
    assert_eq!(vacuum.index.phase_index, 1.0);
    assert_eq!(vacuum.scaled_wavelength_derivatives, [0.0; 2]);
    let trace = spectral.refractive_index_jet(&[1e-20, 1e-20]).unwrap();
    assert_eq!(trace.index.phase_index, 1.0);
    assert_eq!(trace.index.group_index, 1.0);
    assert!(trace.scaled_wavelength_derivatives[0] < 0.0);
    assert!(trace.scaled_wavelength_derivatives[1] > 0.0);
    for densities in [
        &[][..],
        &[-1.0, 1.0],
        &[f64::NAN, 1.0],
        &[f64::INFINITY, 1.0],
    ] {
        assert!(spectral.refractive_index_jet(densities).is_err());
    }
}

#[test]
fn g0_spectral_source_support_order_and_identity_are_explicit() {
    let (reference, shapes) = fixture();
    let original = shapes[0].reference();
    for b in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(
            MolecularCauchyDispersion::try_new(
                original,
                b,
                0.0,
                [380.0, 780.0],
                ContentHash([2; 32])
            )
            .is_err()
        );
    }
    for range in [
        [379.0, 780.0],
        [380.0, f64::NAN],
        [550.0, 550.0],
        [600.0, 700.0],
    ] {
        assert!(
            MolecularCauchyDispersion::try_new(original, 0.0, 0.0, range, ContentHash([2; 32]))
                .is_err()
        );
    }
    assert!(
        MolecularCauchyDispersion::try_new(
            original,
            1.0,
            0.0,
            [380.0, 780.0],
            ContentHash([2; 32])
        )
        .is_err()
    );
    assert!(
        MolecularCauchyDispersion::try_new(
            original,
            0.0,
            0.0,
            [380.0, 780.0],
            ContentHash([0; 32])
        )
        .is_err()
    );
    assert_eq!(
        reference.with_dispersion(&[]),
        Err(GasScatteringError::Components)
    );
    assert_eq!(
        reference.with_dispersion(&[shapes[1], shapes[0]]),
        Err(GasScatteringError::StateMismatch)
    );
    let narrow = [
        MolecularCauchyDispersion::try_new(
            shapes[0].reference(),
            2.0e-6,
            1.0e-7,
            [400.0, 700.0],
            ContentHash([2; 32]),
        )
        .unwrap(),
        MolecularCauchyDispersion::try_new(
            shapes[1].reference(),
            8.0e-6,
            3.0e-7,
            [450.0, 750.0],
            ContentHash([3; 32]),
        )
        .unwrap(),
    ];
    let model = reference.with_dispersion(&narrow).unwrap();
    assert_eq!(model.wavelength_range_nm(), [450.0, 700.0]);
    assert_eq!(
        model.refractive_index(&[1.0, 1.0], PartialDensityBasis::MassKgPerM3, 440.0),
        Err(GasScatteringError::RefractivityDomain)
    );
    assert_eq!(
        model.with_dispersion(&shapes).unwrap().content_identity(),
        reference
            .with_dispersion(&shapes)
            .unwrap()
            .content_identity()
    );
    let mut sourced = shapes;
    sourced[0] = MolecularCauchyDispersion::try_new(
        original,
        2.0e-6,
        1.0e-7,
        [380.0, 780.0],
        ContentHash([9; 32]),
    )
    .unwrap();
    assert_ne!(
        reference
            .with_dispersion(&sourced)
            .unwrap()
            .content_identity(),
        reference
            .with_dispersion(&shapes)
            .unwrap()
            .content_identity()
    );
}
