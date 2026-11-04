use super::*;
use fs_blake3::ContentHash;
use fs_render::volumes::heating::fluorescence::{FluorescenceSource, FluorescentOpticalHeating};

fn source<const M: usize>(
    quantum_yield: f64,
    excitation: (f64, f64),
    wavelengths: [f64; M],
    weights: [f64; M],
) -> FluorescenceSource<M> {
    FluorescenceSource::try_new(
        ContentHash([93; 32]),
        excitation,
        (290.0, 800.0),
        quantum_yield,
        wavelengths,
        weights,
    )
    .unwrap()
}

fn spectrum(absorption: f64) -> MassAbsorptionSpectrum {
    MassAbsorptionSpectrum::source_table(
        &[
            (380.0, absorption),
            (500.0, absorption),
            (600.0, 0.0),
            (800.0, 0.0),
        ],
        ContentHash([94; 32]),
    )
    .unwrap()
}

#[test]
fn photon_weighted_stokes_emission_and_native_heat_conserve_the_original_pump() {
    with_cx(|cx| {
        let density = [1.0, 1.0];
        let medium = grid(cx, &density, spectrum(core::f64::consts::LN_2));
        let deposition = medium
            .deposit_beam_power(
                cx,
                Point3::new(-1.0, 0.5, 0.5),
                Vec3::new(1.0, 0.0, 0.0),
                4.0,
                &[400.0, 500.0],
                &[10.0, 6.0],
                2,
            )
            .unwrap();
        let original = deposition.clone();
        let thermal = thermal(
            vec![cell(1, 2.0, 1.0, 0.0), cell(2, 1.0, 1.0, 0.0)],
            &SENSIBLE,
            vec![],
        );
        let heating = FluorescentOpticalHeating::try_new(
            &medium,
            source(0.8, (380.0, 500.0), [600.0, 800.0], [1.0, 3.0]),
            deposition,
            &[(0, CellId(1)), (1, CellId(2))],
            &thermal,
            2,
            1e-12,
            || Ok(()),
        )
        .unwrap();
        // Independent photon accounting: total absorbed 7.5 W at 400 nm and
        // 4.5 W at 500 nm; photons/s = P lambda/(h c), with SI lambda.
        const HC_J_NM: f64 = 6.626_070_15e-34 * 299_792_458.0 * 1e9;
        let absorbed_photons = (7.5 * 400.0 + 4.5 * 500.0) / HC_J_NM;
        let emitted_photons = heating.power().escaped_w[0] * 600.0 / HC_J_NM
            + heating.power().escaped_w[1] * 800.0 / HC_J_NM;
        assert!((emitted_photons / absorbed_photons - 0.8).abs() < 2e-15);
        close(heating.power().escaped_w[0], 0.8 * 0.25 * 5250.0 / 600.0);
        close(heating.power().escaped_w[1], 0.8 * 0.75 * 5250.0 / 800.0);
        close(heating.power().nonradiative_heat_w, 12.0 * 0.2);
        close(
            heating.power().stokes_heat_w,
            12.0 * 0.8 - heating.power().escaped_w.iter().sum::<f64>(),
        );
        let heat = 12.0 - heating.power().escaped_w.iter().sum::<f64>();
        close(heating.pump_power().thermal_input_w, heat);
        assert_eq!(heating.pump_power().absorbed_w, original.absorbed_power_w);
        close(heating.sources_w()[&CellId(1)], heat * 2.0 / 3.0);
        close(heating.sources_w()[&CellId(2)], heat / 3.0);
        let (heated, heating, audit) = thermal
            .advance_observed(&heating, 2.0, config(0.125), || Ok(()))
            .unwrap();
        close(heated.temperature_k(CellId(1)).unwrap(), 300.0 + heat / 6.0);
        close(heated.temperature_k(CellId(2)).unwrap(), 300.0 + heat / 6.0);
        close(audit.external_heat_j, 2.0 * heat);
        close(heating.pump_energy().thermal_input_j, audit.external_heat_j);
        close(
            heating.energy().escaped_j.iter().sum::<f64>() + audit.external_heat_j,
            heating.pump_energy().absorbed_j.iter().sum::<f64>(),
        );
        close(heating.energy().duration_s, 2.0);
    });
}

#[test]
fn thin_opaque_and_zero_yield_sources_keep_absorption_and_escape_separate() {
    with_cx(|cx| {
        let density = [1.0];
        for absorption in [0.0, 1e-20, 1e4] {
            for yield_ in [0.0, 1.0] {
                let medium = grid(cx, &density, spectrum(absorption));
                let deposition = medium
                    .deposit_beam_power(
                        cx,
                        Point3::new(0.0, 0.5, 0.5),
                        Vec3::new(1.0, 0.0, 0.0),
                        1.0,
                        &[400.0],
                        &[10.0],
                        1,
                    )
                    .unwrap();
                let absorbed = 10.0 * -(-absorption).exp_m1();
                let bindings = if absorbed == 0.0 {
                    vec![]
                } else {
                    vec![(0, CellId(1))]
                };
                let thermal = thermal(vec![cell(1, 1.0, 1.0, 0.0)], &SENSIBLE, vec![]);
                let heating = FluorescentOpticalHeating::try_new(
                    &medium,
                    source(yield_, (380.0, 500.0), [800.0], [1.0]),
                    deposition,
                    &bindings,
                    &thermal,
                    1,
                    1e-13,
                    || Ok(()),
                )
                .unwrap();
                let relative = |actual: f64, expected: f64| {
                    assert!(
                        (actual - expected).abs() <= 3e-14 * expected.abs().max(1e-30),
                        "{actual} != {expected}"
                    );
                };
                relative(heating.pump_power().absorbed_w[0], absorbed);
                relative(heating.power().escaped_w[0], absorbed * yield_ * 0.5);
                relative(
                    heating.pump_power().thermal_input_w,
                    absorbed * (1.0 - 0.5 * yield_),
                );
            }
        }
    });
}

#[test]
fn near_zero_stokes_heat_remains_a_positive_local_loss() {
    with_cx(|cx| {
        let pump: f64 = 512.0;
        let emitted = f64::from_bits(pump.to_bits() + 1);
        let absorption = MassAbsorptionSpectrum::source_table(
            &[(400.0, 1e4), (pump, 1e4), (emitted, 0.0), (800.0, 0.0)],
            ContentHash([95; 32]),
        )
        .unwrap();
        let density = [1.0];
        let medium = grid(cx, &density, absorption);
        let deposition = medium
            .deposit_beam_power(
                cx,
                Point3::new(0.0, 0.5, 0.5),
                Vec3::new(1.0, 0.0, 0.0),
                1.0,
                &[pump],
                &[1.0],
                1,
            )
            .unwrap();
        let thermal = thermal(vec![cell(1, 1.0, 1.0, 0.0)], &SENSIBLE, vec![]);
        let heating = FluorescentOpticalHeating::try_new(
            &medium,
            source(1.0, (400.0, pump), [emitted], [1.0]),
            deposition,
            &[(0, CellId(1))],
            &thermal,
            1,
            1e-14,
            || Ok(()),
        )
        .unwrap();
        let expected = (emitted - pump) / emitted;
        assert!(heating.power().stokes_heat_w > 0.0);
        assert_eq!(heating.power().stokes_heat_w, expected);
        assert_eq!(heating.pump_power().thermal_input_w, expected);
        assert_eq!(heating.power().nonradiative_heat_w, 0.0);
    });
}

#[test]
fn source_and_medium_domains_refuse_unmodeled_fluorescence() {
    assert!(
        FluorescenceSource::try_new(
            ContentHash([1; 32]),
            (400.0, 600.0),
            (290.0, 800.0),
            0.8,
            [500.0],
            [1.0]
        )
        .is_err()
    );
    assert!(
        FluorescenceSource::try_new(
            ContentHash([1; 32]),
            (400.0, 500.0),
            (290.0, 800.0),
            1.01,
            [600.0],
            [1.0]
        )
        .is_err()
    );
    with_cx(|cx| {
        let density = [1.0, 1.0];
        let clear_emission = grid(cx, &density, spectrum(1.0));
        let opaque_emission = grid(
            cx,
            &density,
            MassAbsorptionSpectrum::source_table(
                &[(380.0, 1.0), (800.0, 1.0)],
                ContentHash([96; 32]),
            )
            .unwrap(),
        );
        let thermal = thermal(vec![cell(1, 1.0, 1.0, 0.0)], &SENSIBLE, vec![]);
        let deposit = |medium: &AbsorbingGrid<'_>| {
            medium
                .deposit_beam_power(
                    cx,
                    Point3::new(0.0, 0.5, 0.5),
                    Vec3::new(1.0, 0.0, 0.0),
                    1.0,
                    &[400.0],
                    &[1.0],
                    1,
                )
                .unwrap()
        };
        let model = source(0.8, (380.0, 500.0), [600.0], [1.0]);
        assert!(
            FluorescentOpticalHeating::try_new(
                &opaque_emission,
                model,
                deposit(&opaque_emission),
                &[(0, CellId(1))],
                &thermal,
                2,
                1e-12,
                || Ok(())
            )
            .is_err()
        );
        assert!(
            FluorescentOpticalHeating::try_new(
                &clear_emission,
                model,
                deposit(&opaque_emission),
                &[(0, CellId(1))],
                &thermal,
                2,
                1e-12,
                || Ok(())
            )
            .is_err()
        );
        assert!(
            FluorescentOpticalHeating::try_new(
                &clear_emission,
                source(0.8, (380.0, 500.0), [900.0], [1.0]),
                deposit(&clear_emission),
                &[(0, CellId(1))],
                &thermal,
                2,
                1e-12,
                || Ok(())
            )
            .is_err()
        );
        assert!(
            FluorescentOpticalHeating::try_new(
                &clear_emission,
                model,
                deposit(&clear_emission),
                &[(0, CellId(1))],
                &thermal,
                1,
                1e-12,
                || Ok(())
            )
            .is_err()
        );
        assert!(
            FluorescentOpticalHeating::try_new(
                &clear_emission,
                model,
                deposit(&clear_emission),
                &[(0, CellId(1))],
                &thermal,
                2,
                1e-12,
                || Err(DepositionError::Cancelled)
            )
            .is_err()
        );
    });
}

#[test]
fn source_temperature_domain_failure_rolls_back_native_heat_and_radiation() {
    with_cx(|cx| {
        let density = [1.0];
        let medium = grid(cx, &density, spectrum(1e4));
        let deposition = medium
            .deposit_beam_power(
                cx,
                Point3::new(0.0, 0.5, 0.5),
                Vec3::new(1.0, 0.0, 0.0),
                1.0,
                &[400.0],
                &[10.0],
                1,
            )
            .unwrap();
        let thermal = thermal(vec![cell(1, 1.0, 1.0, 0.0)], &SENSIBLE, vec![]);
        let model = FluorescenceSource::try_new(
            ContentHash([97; 32]),
            (380.0, 500.0),
            (299.0, 300.1),
            0.8,
            [800.0],
            [1.0],
        )
        .unwrap();
        let heating = FluorescentOpticalHeating::try_new(
            &medium,
            model,
            deposition,
            &[(0, CellId(1))],
            &thermal,
            1,
            1e-12,
            || Ok(()),
        )
        .unwrap();
        assert!(
            thermal
                .advance_observed(&heating, 1.0, config(0.1), || Ok(()))
                .is_err()
        );
        assert_eq!(thermal.time_s(), 0.0);
        assert_eq!(thermal.temperature_k(CellId(1)).unwrap(), 300.0);
        assert_eq!(heating.pump_energy().thermal_input_j, 0.0);
        assert_eq!(heating.energy().escaped_j, [0.0]);
    });
}
