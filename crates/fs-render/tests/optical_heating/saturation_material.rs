// Included inside phosphorescence::material_rates to reuse its original
// material-card fixtures without exposing or duplicating their implementation.

#[test]
fn g1_saturable_material_rates_refresh_on_the_actual_thermal_clock() {
    use fs_render::volumes::heating::luminescence::saturation::{
        PrescribedPhotonPump, SaturableEmitterCell, SaturableLuminescentHeating,
    };

    with_cx(|cx| {
        let mut data = claims();
        // Source singlet radiative rate: 1/s at 300 K, increasing by .01/s/K.
        // Remove crossing to give an independently soluble saturation oracle.
        data[2].value = PropertyValue::Curve {
            abscissa: "T".into(),
            abscissa_dims: Temperature::DIMS,
            knots: vec![(300.0, 1.0), (400.0, 2.0)],
            dims: PHOSPHORESCENCE_RATE_DIMS,
        };
        data[2].validity = ValidityDomain::unconstrained().with("T", 300.0, 400.0);
        data[2].interpolation = InterpolationPolicy::LinearInside;
        data[6].value = PropertyValue::Scalar {
            value: 0.0,
            dims: PHOSPHORESCENCE_RATE_DIMS,
        };
        let (singlet, coupling) = native_sources(card(data), &point()).unwrap();
        // The fixture's legacy spectral quantum yield is zero. The original
        // material's radiative rate must still produce nonzero singlet light.
        let density = [1.0];
        let medium = grid(cx, &density, MassAbsorptionSpectrum::constant(0.0).unwrap());
        let thermal = thermal(vec![cell(1, 1.0, 1.0, 0.0)], &SENSIBLE, vec![]);
        let pump = PrescribedPhotonPump::try_new(
            ContentHash([244; 32]),
            [400.0; 2],
            [1e-20, 0.0],
            [9e20, 0.0],
        )
        .unwrap();
        let count = 1e19;
        let heating = SaturableLuminescentHeating::try_new(
            &medium,
            singlet,
            coupling,
            &[SaturableEmitterCell {
                grid_cell: 0,
                thermal_cell: CellId(1),
                initial_populations: [count, 0.0, 0.0],
                pump,
            }],
            &thermal,
            1,
            1e-8,
            None,
            || Ok(()),
        )
        .unwrap()
        .with_pump_window(0.0, 1.0)
        .unwrap();
        let (warm, charged, audit) = thermal
            .advance_observed(&heating, 1.0, config(1.0), || Ok(()))
            .unwrap();
        let excited = count * 0.9 * (1.0 - (-10.0_f64).exp());
        let singlet_residence = count * 0.9 * (1.0 - (1.0 - (-10.0_f64).exp()) / 10.0);
        close(charged.populations(0).unwrap()[1] / count, excited / count);
        close(
            charged.energy().escaped_j[0],
            singlet_residence * HC_J_NM / 650.0,
        );
        close(warm.energy_j(), charged.thermal_input_j());
        close(audit.external_heat_j, charged.thermal_input_j());
        let kr = 1.0 + 0.01 * (warm.temperature_k(CellId(1)).unwrap() - 300.0);
        assert!(kr > 1.0);
        close(
            charged.cell_emission_w()[0].1[0],
            excited * kr * HC_J_NM / 650.0,
        );
        let (_, dark, _) = warm
            .advance_observed(&charged, 1.5, config(0.5), || Ok(()))
            .unwrap();
        close(
            dark.populations(0).unwrap()[1] / count,
            excited * (-kr * 0.5).exp() / count,
        );
        assert_eq!(
            dark.absorbed_pump_energy_j(),
            charged.absorbed_pump_energy_j()
        );
        assert!(dark.cell_emission_w()[0].1[0] > 0.0);
        close(
            dark.populations(0).unwrap().iter().sum::<f64>() / count,
            1.0,
        );
        close(
            dark.absorbed_pump_energy_j().iter().sum(),
            dark.energy().stored_excitation_j
                + dark.energy().escaped_j.iter().sum::<f64>()
                + dark.thermal_input_j(),
        );
    });
}
