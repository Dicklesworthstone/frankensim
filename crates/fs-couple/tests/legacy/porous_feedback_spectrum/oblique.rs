use super::*;
use fs_couple::porous_heating::feedback::covariance::{
    PorousCovarianceBin, PorousCrossSpectrumHeating,
};
use fs_duct::porous::PressureCovariance;

fn covariance(bins: &[PorousHarmonicDrive]) -> Vec<PorousCovarianceBin> {
    bins.iter()
        .map(|b| PorousCovarianceBin {
            angular_frequency_rad_s: b.angular_frequency_rad_s,
            area_m2: b.area_m2,
            reference_impedance_pa_s_m3: b.reference_impedance_pa_s_m3,
            covariance: PressureCovariance::from_peak_pressure(b.incident_pressure_pa).unwrap(),
        })
        .collect()
}

#[test]
fn g0_zero_tangential_constructors_preserve_original_native_arithmetic() {
    let materials = materials();
    let routes = routes(&materials);
    let thermal = thermal(events([0.0, 8000.0], [0.0; 2]), [&SENSIBLE; 2]);
    let bins = bins();
    let plain_tone = TemperaturePorousHeating::try_new(
        &thermal,
        &routes,
        gas(),
        bins[0],
        ComposeOptions::default(),
        || Ok(()),
    )
    .unwrap();
    let zero_tone = TemperaturePorousHeating::try_new_oblique(
        &thermal,
        &routes,
        gas(),
        bins[0],
        0.0,
        ComposeOptions::default(),
        || Ok(()),
    )
    .unwrap();
    assert_eq!(plain_tone.dissipation(), zero_tone.dissipation());
    let plain =
        PorousSpectrumHeating::try_new(&thermal, &routes, gas(), &bins, options(), || Ok(()))
            .unwrap();
    let zero = PorousSpectrumHeating::try_new_oblique(
        &thermal,
        &routes,
        gas(),
        &bins,
        &[0.0; 3],
        options(),
        || Ok(()),
    )
    .unwrap();
    assert_eq!(plain.power(), zero.power());
    let plain = plain
        .advance_to(&thermal, 0.02, config(0.01), || Ok(()))
        .unwrap();
    let zero = zero
        .advance_to(&thermal, 0.02, config(0.01), || Ok(()))
        .unwrap();
    assert_eq!(plain.step.thermal.energy_j(), zero.step.thermal.energy_j());
    assert_eq!(plain.work, zero.work);
    let covariance = covariance(&bins);
    let plain = PorousCrossSpectrumHeating::try_new(
        &thermal,
        &routes,
        gas(),
        &covariance,
        options(),
        || Ok(()),
    )
    .unwrap();
    let zero = PorousCrossSpectrumHeating::try_new_oblique(
        &thermal,
        &routes,
        gas(),
        &covariance,
        &[0.0; 3],
        options(),
        || Ok(()),
    )
    .unwrap();
    assert_eq!(plain.power(), zero.power());
    assert_eq!(zero.tangential_wavenumber_rad_m(3), None);
}

#[test]
fn g3_canonical_oblique_bins_preserve_modes_partial_coherence_and_shared_refresh() {
    let materials = materials();
    let routes = routes(&materials);
    let thermal = thermal(events([0.0, 8000.0], [0.0; 2]), [&SENSIBLE; 2]);
    let ambient = GasState::try_new(&SPEC, 300.0, 101325.0).unwrap();
    let mut bins = bins();
    let cosines = [0.7_f64, 0.4, 0.85];
    let tangential = std::array::from_fn::<_, 3, _>(|i| {
        bins[i].reference_impedance_pa_s_m3 /= cosines[i];
        bins[i].angular_frequency_rad_s / ambient.sound_speed
            * (1.0 - cosines[i] * cosines[i]).sqrt()
    });
    let coherent = PorousSpectrumHeating::try_new_oblique(
        &thermal,
        &routes,
        gas(),
        &bins,
        &tangential,
        options(),
        || Ok(()),
    )
    .unwrap();
    for (canonical, original) in [1, 2, 0].into_iter().enumerate() {
        assert_eq!(
            coherent.bins()[canonical].drive().angular_frequency_rad_s,
            bins[original].angular_frequency_rad_s
        );
        assert_eq!(
            coherent.bins()[canonical].tangential_wavenumber_rad_m(),
            tangential[original]
        );
    }
    let covariance = covariance(&bins);
    let cross = PorousCrossSpectrumHeating::try_new_oblique(
        &thermal,
        &routes,
        gas(),
        &covariance,
        &tangential,
        options(),
        || Ok(()),
    )
    .unwrap();
    assert_eq!(coherent.power(), cross.power());
    for (canonical, original) in [1, 2, 0].into_iter().enumerate() {
        assert_eq!(
            cross.tangential_wavenumber_rad_m(canonical),
            Some(tangential[original])
        );
    }
    // Uncorrelated left and right RMS waves are an incoherent sum, not
    // one phasor pair. Both share the SAME oblique mode per frequency.
    let independent: Vec<_> = covariance
        .iter()
        .map(|b| PorousCovarianceBin {
            covariance: PressureCovariance::try_new([
                C64::from_re(20_000.0),
                C64::ZERO,
                C64::ZERO,
                C64::from_re(5000.0),
            ])
            .unwrap(),
            ..*b
        })
        .collect();
    let mixed = PorousCrossSpectrumHeating::try_new_oblique(
        &thermal,
        &routes,
        gas(),
        &independent,
        &tangential,
        options(),
        || Ok(()),
    )
    .unwrap();
    let initial_power = mixed.power();
    let mut expected = 0.0;
    for pressure in [
        [C64::from_re(200.0), C64::ZERO],
        [C64::ZERO, C64::from_re(100.0)],
    ] {
        let drives: Vec<_> = bins
            .iter()
            .map(|b| PorousHarmonicDrive {
                incident_pressure_pa: pressure,
                ..*b
            })
            .collect();
        expected += PorousSpectrumHeating::try_new_oblique(
            &thermal,
            &routes,
            gas(),
            &drives,
            &tangential,
            options(),
            || Ok(()),
        )
        .unwrap()
        .power()
        .dissipation
        .absorbed_power_w;
    }
    close(mixed.power().dissipation.absorbed_power_w, expected, 3e-12);
    let evolved = mixed
        .advance_to(&thermal, 0.02, config(0.01), || Ok(()))
        .unwrap();
    assert_eq!(evolved.work.gas_evaluations, 4);
    assert_eq!(evolved.work.stack_solves, 6);
    for (canonical, original) in [1, 2, 0].into_iter().enumerate() {
        assert_eq!(
            evolved.heating.tangential_wavenumber_rad_m(canonical),
            Some(tangential[original])
        );
        assert_eq!(
            evolved.heating.bins()[canonical].reference_impedance_pa_s_m3,
            bins[original].reference_impedance_pa_s_m3
        );
    }
    assert!(evolved.step.thermal.energy_j() > thermal.energy_j());
    for invalid in [&tangential[..2], &[0.0, f64::NAN, 0.0], &[0.0, -1.0, 0.0]] {
        assert!(
            PorousSpectrumHeating::try_new_oblique(
                &thermal,
                &routes,
                gas(),
                &bins,
                invalid,
                options(),
                || Ok(())
            )
            .is_err()
        );
        assert!(
            PorousCrossSpectrumHeating::try_new_oblique(
                &thermal,
                &routes,
                gas(),
                &independent,
                invalid,
                options(),
                || Ok(())
            )
            .is_err()
        );
    }
    assert!(
        mixed
            .advance_to(&thermal, 0.02, config(0.01), || Err(
                DepositionError::Cancelled
            ))
            .is_err()
    );
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(mixed.power(), initial_power);
}
