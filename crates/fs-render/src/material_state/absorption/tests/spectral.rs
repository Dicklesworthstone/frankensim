//! Finite-band absorption follows the material owner's original alpha law.

use super::*;

#[derive(Clone, Copy)]
enum Family {
    Tabulated,
    Sellmeier,
    Conductor,
}

fn fixture(
    family: Family,
    lossless: bool,
) -> (
    MaterialCard,
    ResolvedMaterialStatePoint,
    SpectralOpticalStatePoint,
) {
    let length = QuantitySpec::dimensional(Dims([1, 0, 0, 0, 0, 0]));
    let dimensionless = QuantitySpec::dimensional(Dims::NONE);
    let mut claims = vec![curve(density_key(), 1000.0, 800.0)];
    let wavelengths_m = [
        8e-6, 8.25e-6, 9e-6, 9.5e-6, 10e-6, 10.25e-6, 11e-6, 12e-6, 14e-6,
    ];
    let alpha = [2.0, 3.0, 5.0, 7.0, 11.0, 13.0, 17.0, 19.0, 23.0];
    for (index, wavelength_m) in wavelengths_m.into_iter().enumerate() {
        claims.push(curve(
            PropertyKey::with_quantity(SPECTRAL_OPTICAL_WAVELENGTH_M_PROPERTIES[index], length),
            wavelength_m,
            if index == 8 { 16e-6 } else { wavelength_m },
        ));
        let (name, quantity, loss) = if matches!(family, Family::Conductor) {
            claims.push(curve(
                PropertyKey::with_quantity(SPECTRAL_CONDUCTOR_ETA_PROPERTIES[index], dimensionless),
                1.5,
                1.5,
            ));
            (
                SPECTRAL_CONDUCTOR_K_PROPERTIES[index],
                SPECTRAL_CONDUCTOR_K_QUANTITY,
                0.01 * (index + 1) as f64,
            )
        } else {
            if matches!(family, Family::Tabulated) {
                claims.push(curve(
                    PropertyKey::with_quantity(
                        SPECTRAL_DIELECTRIC_IOR_PROPERTIES[index],
                        dimensionless,
                    ),
                    1.5,
                    1.5,
                ));
            }
            (
                SPECTRAL_DIELECTRIC_ABSORPTION_PER_M_PROPERTIES[index],
                SPECTRAL_DIELECTRIC_ABSORPTION_QUANTITY,
                alpha[index],
            )
        };
        let loss = if lossless { 0.0 } else { loss };
        claims.push(curve(
            PropertyKey::with_quantity(name, quantity),
            loss,
            2.0 * loss,
        ));
    }
    if matches!(family, Family::Sellmeier) {
        for index in 0..3 {
            claims.push(curve(
                PropertyKey::with_quantity(
                    VISIBLE_DIELECTRIC_SELLMEIER_B_PROPERTIES[index],
                    dimensionless,
                ),
                [1.1, 0.2, 0.8][index],
                [1.1, 0.2, 0.8][index],
            ));
            let resonance = [0.01, 0.04, 0.1][index] * 1e-12;
            claims.push(curve(
                PropertyKey::new(
                    VISIBLE_DIELECTRIC_SELLMEIER_C_M2_PROPERTIES[index],
                    Dims([2, 0, 0, 0, 0, 0]),
                ),
                resonance,
                resonance,
            ));
        }
        for (name, bound) in SPECTRAL_DIELECTRIC_SELLMEIER_WAVELENGTH_BOUNDS_M_PROPERTIES
            .into_iter()
            .zip([9e-6, 12e-6])
        {
            claims.push(curve(
                PropertyKey::with_quantity(name, length),
                bound,
                bound,
            ));
        }
    }
    let requirements = claims
        .iter()
        .map(|claim| {
            ScalarPropertyRequirement::try_with_key(&claim.key, ScalarAdmissibility::Finite)
                .unwrap()
        })
        .collect::<Vec<_>>();
    let card = card(claims);
    let accepted = resolve_material_state_point(
        &card,
        &point(300.0),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let optical = match family {
        Family::Tabulated => SpectralOpticalStatePoint::TabulatedDielectric(
            SpectralTabulatedDielectricStatePoint::try_from_material(accepted.clone()).unwrap(),
        ),
        Family::Sellmeier => SpectralOpticalStatePoint::DielectricSellmeier(
            SpectralSellmeierStatePoint::try_from_material(accepted.clone()).unwrap(),
        ),
        Family::Conductor => SpectralOpticalStatePoint::Conductor(
            SpectralConductorStatePoint::try_from_material(accepted.clone()).unwrap(),
        ),
    };
    (card, accepted, optical)
}

#[test]
fn g1_infrared_tabulated_alpha_drives_original_density_opacity() {
    let (_, accepted, optical) = fixture(Family::Tabulated, false);
    let SpectralOpticalStatePoint::TabulatedDielectric(source) = &optical else {
        unreachable!()
    };
    let binding = MaterialAbsorptionBinding::try_spectral_tabulated_dielectric(
        &accepted,
        source,
        &density_key(),
    )
    .unwrap();
    assert_eq!(
        binding,
        MaterialAbsorptionBinding::try_spectral_optical(&accepted, &optical, &density_key())
            .unwrap()
    );
    assert_eq!(binding.source_identity(), source.identity());
    assert_eq!(binding.accepted_temperature_k(), Some(300.0));
    assert!(core::ptr::eq(binding.optical_state(), source.resolved()));
    for sample in source.samples() {
        close(
            binding
                .mass_absorption_at(sample.wavelength_m * 1e9)
                .unwrap(),
            sample.absorption_per_m / 1000.0,
        );
    }
    let bounds = binding.wavelength_bounds_nm();
    let wavelengths = [bounds[0], 9750.0, bounds[1]];
    let samples = binding.sample(&wavelengths).unwrap();
    let density = [1000.0];
    let grid = VolumeGrid::new([1; 3], &density, [0.0; 3], [1.0, 1.0, 0.25]);
    let image =
        render_density_spectral_box(&grid, settings(), samples.as_density_spectrum(1.0, 1.0))
            .unwrap();
    for (observed, alpha) in image.samples().iter().zip([2.0, 9.0, 23.0]) {
        close(observed.optical_depth(), alpha * 0.25);
        close(observed.transmittance(), (-alpha * 0.25_f64).exp());
    }
}

#[test]
fn g1_sellmeier_absorption_keeps_its_independent_domain_and_native_alpha() {
    let (_, accepted, optical) = fixture(Family::Sellmeier, false);
    let SpectralOpticalStatePoint::DielectricSellmeier(source) = &optical else {
        unreachable!()
    };
    let binding =
        MaterialAbsorptionBinding::try_spectral_sellmeier(&accepted, source, &density_key())
            .unwrap();
    assert_eq!(
        binding,
        MaterialAbsorptionBinding::try_spectral_optical(&accepted, &optical, &density_key())
            .unwrap()
    );
    assert_eq!(source.dispersion_wavelength_bounds_m(), (9e-6, 12e-6));
    let bounds = binding.wavelength_bounds_nm();
    assert!(bounds[0] < 9000.0 && bounds[1] > 12000.0);
    let wavelengths = [bounds[0], 13000.0, bounds[1]];
    let samples = binding.sample(&wavelengths).unwrap();
    for (coefficient, alpha) in samples
        .mass_absorption_m2_kg()
        .iter()
        .zip([2.0, 21.0, 23.0])
    {
        close(*coefficient, alpha / 1000.0);
    }
    let density = [1000.0];
    let grid = VolumeGrid::new([1; 3], &density, [0.0; 3], [1.0, 1.0, 0.25]);
    let image =
        render_density_spectral_box(&grid, settings(), samples.as_density_spectrum(1.0, 1.0))
            .unwrap();
    close(image.samples()[1].optical_depth(), 21.0 * 0.25);
    assert_eq!(samples.source_identity(), source.identity());
}

#[test]
fn g3_finite_band_refresh_updates_alpha_density_coordinates_and_refuses_stale_sources() {
    for family in [Family::Tabulated, Family::Sellmeier] {
        let (card, accepted, optical) = fixture(family, false);
        let before =
            MaterialAbsorptionBinding::try_spectral_optical(&accepted, &optical, &density_key())
                .unwrap();
        let hot = accepted.at_query_point(&card, &point(400.0)).unwrap();
        let hot_optical = optical.at_query_point(&card, &point(400.0)).unwrap();
        let after =
            MaterialAbsorptionBinding::try_spectral_optical(&hot, &hot_optical, &density_key())
                .unwrap();
        let wavelength = before.wavelength_bounds_nm()[0];
        close(
            after.mass_absorption_at(wavelength).unwrap()
                / before.mass_absorption_at(wavelength).unwrap(),
            2.5,
        );
        assert_eq!(after.material_density_kg_m3(), 800.0);
        assert_eq!(after.accepted_temperature_k(), Some(400.0));
        close(after.wavelength_bounds_nm()[1], 16000.0);
        assert!(before.mass_absorption_at(15000.0).is_err());
        assert!(after.mass_absorption_at(15000.0).is_ok());
        assert_ne!(after.identity(), before.identity());
        assert_ne!(after.source_identity(), before.source_identity());
        assert!(
            MaterialAbsorptionBinding::try_spectral_optical(&hot, &optical, &density_key())
                .is_err()
        );
        assert!(
            MaterialAbsorptionBinding::try_spectral_optical(
                &accepted,
                &hot_optical,
                &density_key()
            )
            .is_err()
        );
        card.claims()
            .verify_receipt(&after.density_property().answer().receipt)
            .unwrap();
        assert!(core::ptr::eq(
            after.density_property(),
            hot.property_by_key(&density_key()).unwrap()
        ));
    }
}

#[test]
fn g0_lossless_finite_band_support_stays_closed_and_conductor_dispatch_preserves_identity() {
    for family in [Family::Tabulated, Family::Sellmeier, Family::Conductor] {
        let (_, accepted, optical) = fixture(family, true);
        let binding =
            MaterialAbsorptionBinding::try_spectral_optical(&accepted, &optical, &density_key())
                .unwrap();
        let bounds = binding.wavelength_bounds_nm();
        assert_eq!(
            binding.sample(&bounds).unwrap().mass_absorption_m2_kg(),
            &[0.0, 0.0]
        );
        for wavelength in [bounds[0].next_down(), bounds[1].next_up(), 550.0, f64::NAN] {
            assert!(matches!(
                binding.mass_absorption_at(wavelength),
                Err(MaterialAbsorptionError::WavelengthOutsideSource { .. })
            ));
        }
        if let SpectralOpticalStatePoint::Conductor(source) = &optical {
            let original = MaterialAbsorptionBinding::try_spectral(
                &accepted,
                source.complex_index(),
                &density_key(),
            )
            .unwrap();
            assert_eq!(binding, original);
            assert_eq!(binding.identity(), original.identity());
        }
    }
    assert_eq!(
        validate_si_absorption_band((1.0, f64::MAX)),
        Err(MaterialAbsorptionError::UnrepresentableWavelengthBand)
    );
}
