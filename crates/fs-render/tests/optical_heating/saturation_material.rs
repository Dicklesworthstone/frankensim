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

mod depleting_material {
    use super::*;
    use fs_material::phosphorescence::pump::{
        GROUND_STATE_PUMP_CROSS_SECTION_DIMS, GROUND_STATE_PUMP_CROSS_SECTION_PROPERTIES,
        GroundStatePumpStatePoint,
    };
    use fs_render::volumes::heating::feedback::IncidentBeam;
    use fs_render::volumes::heating::luminescence::saturation::beam::{
        BeamEmitterCrossSections, PumpColumnLimits,
    };
    use fs_render::volumes::heating::luminescence::saturation::{
        PrescribedPhotonPump, SaturableEmitterCell, SaturableLuminescentHeating,
    };

    const COUNT: f64 = 1e19;

    // Independently authored affine molecular law. The material database only
    // receives four corner samples, not this evaluation function.
    fn areas(wavelength_nm: f64, temperature_k: f64) -> [f64; 2] {
        let wavelength = (wavelength_nm - 400.0) / 50.0;
        let temperature = temperature_k - 300.0;
        [
            (1.0 + wavelength + 0.2 * temperature) * 1e-19,
            (0.3 + 0.4 * wavelength + 0.1 * temperature) * 1e-19,
        ]
    }

    fn molecular_source(
        maximum_temperature_k: f64,
    ) -> (Arc<MaterialCard>, GroundStatePumpStatePoint) {
        let wavelength_quantity = QuantitySpec::dimensional(fs_qty::Dims([1, 0, 0, 0, 0, 0]));
        let temperature_quantity = QuantitySpec::semantic(SemanticType::new(
            QuantityKind::AbsoluteTemperature,
            ValueForm::Static,
        ));
        let mut input = claims();
        for (target, name) in GROUND_STATE_PUMP_CROSS_SECTION_PROPERTIES
            .into_iter()
            .enumerate()
        {
            let mut claim = input[0].clone();
            claim.key = PropertyKey::new(name, GROUND_STATE_PUMP_CROSS_SECTION_DIMS);
            claim.value = PropertyValue::Surface2D {
                axes: ["lambda".into(), "T".into()],
                axis_dims: [wavelength_quantity.dims(), Temperature::DIMS],
                samples: [
                    [380.0, 300.0],
                    [500.0, 300.0],
                    [380.0, maximum_temperature_k],
                    [500.0, maximum_temperature_k],
                ]
                .into_iter()
                .map(|[wavelength, temperature]| {
                    (
                        [wavelength * 1e-9, temperature],
                        areas(wavelength, temperature)[target],
                    )
                })
                .collect(),
                triangles: vec![[0, 1, 2], [1, 3, 2]],
                dims: GROUND_STATE_PUMP_CROSS_SECTION_DIMS,
            };
            claim.validity = ValidityDomain::unconstrained()
                .with_quantity("lambda", wavelength_quantity, 380.0 * 1e-9, 500.0 * 1e-9)
                .with_quantity("T", temperature_quantity, 300.0, maximum_temperature_k);
            claim.interpolation = InterpolationPolicy::LinearInside;
            claim.provenance.source =
                format!("manufactured temperature-dependent molecular pump: {name}");
            input.push(claim);
        }
        let card = card(input);
        let source = GroundStatePumpStatePoint::resolve(
            &card,
            &point()
                .with_quantity("lambda", wavelength_quantity, 400.0 * 1e-9)
                .unwrap(),
            MaterialPropertySelection::SingleClaimOnly,
            "lambda",
            Some("T"),
        )
        .unwrap();
        (card, source)
    }

    fn prepare(
        cx: &Cx<'_>,
        medium: &AbsorbingGrid<'_>,
        thermal: &DepositionSimulation<'_>,
        card: Arc<MaterialCard>,
        populations: [f64; 3],
        sections: BeamEmitterCrossSections,
    ) -> SaturableLuminescentHeating<2> {
        let (singlet, coupling) = native_sources(card, &point()).unwrap();
        let pump = PrescribedPhotonPump::try_new(
            ContentHash([170; 32]),
            [400.0, 450.0],
            [0.0; 2],
            [0.0; 2],
        )
        .unwrap();
        SaturableLuminescentHeating::try_new(
            medium,
            singlet,
            coupling,
            &[SaturableEmitterCell {
                grid_cell: 0,
                thermal_cell: CellId(1),
                initial_populations: populations,
                pump,
            }],
            thermal,
            1,
            1e-8,
            None,
            || Ok(()),
        )
        .unwrap()
        .with_depleting_beam(
            cx,
            medium,
            ContentHash([171; 32]),
            IncidentBeam {
                origin: Point3::new(-1.0, 0.5, 0.5),
                unit_direction: Vec3::new(1.0, 0.0, 0.0),
                max_distance_m: 3.0,
                wavelengths_nm: [400.0, 450.0],
                incident_power_w: [0.8, 1.2],
            },
            &[sections],
            thermal,
            1,
            PumpColumnLimits::default(),
            || Ok(()),
        )
        .unwrap()
    }

    #[test]
    fn g3_material_depleting_beam_refreshes_both_wavelengths_after_accepted_heating() {
        with_cx(|cx| {
            let (card, source) = molecular_source(310.0);
            let density = [1.0];
            let medium = grid(cx, &density, MassAbsorptionSpectrum::constant(0.0).unwrap());
            // An original independent boundary source ensures a measurable
            // thermal change even at this modest optical pump power.
            let thermal = thermal(vec![cell(1, 1.0, 1.0, 40.0)], &SENSIBLE, vec![]);
            let material = prepare(
                cx,
                &medium,
                &thermal,
                card.clone(),
                [COUNT, 0.0, 0.0],
                BeamEmitterCrossSections::from_material(0, card.clone(), source.clone()).unwrap(),
            );
            let authored = prepare(
                cx,
                &medium,
                &thermal,
                card.clone(),
                [COUNT, 0.0, 0.0],
                BeamEmitterCrossSections::authored(
                    0,
                    ContentHash([172; 32]),
                    [400.0, 450.0].map(|w| areas(w, 300.0)),
                )
                .unwrap(),
            );
            let (warm, charged, audit) = thermal
                .advance_observed(&material, 0.1, config(0.1), || Ok(()))
                .unwrap();
            let (control_warm, first_control, _) = thermal
                .advance_observed(&authored, 0.1, config(0.1), || Ok(()))
                .unwrap();
            assert_eq!(audit.substeps, 1);
            for (actual, expected) in charged
                .populations(0)
                .unwrap()
                .into_iter()
                .zip(first_control.populations(0).unwrap())
            {
                close(actual / COUNT, expected / COUNT);
            }
            close(warm.energy_j(), control_warm.energy_j());
            for lane in 0..2 {
                close(
                    charged.depleting_pump_energy().unwrap().absorbed_j[lane],
                    first_control.depleting_pump_energy().unwrap().absorbed_j[lane],
                );
                close(
                    charged.depleting_pump_energy().unwrap().transmitted_j[lane],
                    first_control.depleting_pump_energy().unwrap().transmitted_j[lane],
                );
            }

            let temperature = warm.temperature_k(CellId(1)).unwrap();
            assert!(temperature > 300.9);
            let refreshed = [400.0, 450.0].map(|w| areas(w, temperature));
            for (lane, wavelength) in [400.0, 450.0].into_iter().enumerate() {
                let native = source
                    .at_wavelength_temperature(&card, wavelength * 1e-9, Some(temperature))
                    .unwrap();
                for target in 0..2 {
                    close(
                        native.cross_sections_m2()[target] * COUNT,
                        refreshed[lane][target] * COUNT,
                    );
                }
            }
            let control = prepare(
                cx,
                &medium,
                &warm,
                card,
                charged.populations(0).unwrap(),
                BeamEmitterCrossSections::authored(0, ContentHash([173; 32]), refreshed).unwrap(),
            );
            let (next, accepted, _) = warm
                .advance_observed(&charged, 0.2, config(0.1), || Ok(()))
                .unwrap();
            let (expected_thermal, expected, _) = warm
                .advance_observed(&control, 0.2, config(0.1), || Ok(()))
                .unwrap();
            for (actual, expected) in accepted
                .populations(0)
                .unwrap()
                .into_iter()
                .zip(expected.populations(0).unwrap())
            {
                close(actual / COUNT, expected / COUNT);
            }
            close(next.energy_j(), expected_thermal.energy_j());
            close(
                accepted.thermal_input_j() - charged.thermal_input_j(),
                expected.thermal_input_j(),
            );
            for lane in 0..2 {
                close(
                    accepted.depleting_pump_energy().unwrap().absorbed_j[lane]
                        - charged.depleting_pump_energy().unwrap().absorbed_j[lane],
                    expected.depleting_pump_energy().unwrap().absorbed_j[lane],
                );
                close(
                    accepted.depleting_pump_energy().unwrap().transmitted_j[lane]
                        - charged.depleting_pump_energy().unwrap().transmitted_j[lane],
                    expected.depleting_pump_energy().unwrap().transmitted_j[lane],
                );
            }
        });
    }

    #[test]
    fn g4_material_depleting_beam_endpoint_source_refusal_preserves_cached_accounts() {
        with_cx(|cx| {
            let (card, source) = molecular_source(300.01);
            let density = [1.0];
            let medium = grid(cx, &density, MassAbsorptionSpectrum::constant(0.0).unwrap());
            let thermal = thermal(vec![cell(1, 1.0, 1.0, 40.0)], &SENSIBLE, vec![]);
            let heating = prepare(
                cx,
                &medium,
                &thermal,
                card.clone(),
                [COUNT, 0.0, 0.0],
                BeamEmitterCrossSections::from_material(0, card, source).unwrap(),
            );
            let before = heating.clone();
            // The thermal law and kinetic spectra support this temperature;
            // only the original pump cross-section source ends at 300.01 K.
            assert!(
                thermal
                    .advance_observed(&heating, 0.1, config(0.1), || Ok(()))
                    .is_err()
            );
            assert_eq!(heating.populations(0), before.populations(0));
            assert_eq!(heating.energy(), before.energy());
            assert_eq!(
                heating.depleting_pump_energy(),
                before.depleting_pump_energy()
            );
            assert_eq!(heating.time_s(), 0.0);
            assert_eq!(thermal.energy_j(), 0.0);
            let (warm, retry, _) = thermal
                .advance_observed(&heating, 1e-4, config(1e-4), || Ok(()))
                .unwrap();
            let (_, fresh, _) = thermal
                .advance_observed(&before, 1e-4, config(1e-4), || Ok(()))
                .unwrap();
            assert!(warm.temperature_k(CellId(1)).unwrap() < 300.01);
            assert_eq!(retry.populations(0), fresh.populations(0));
            assert_eq!(retry.energy(), fresh.energy());
            assert_eq!(retry.depleting_pump_energy(), fresh.depleting_pump_energy());
            assert!(
                retry
                    .depleting_pump_energy()
                    .unwrap()
                    .absorbed_j
                    .iter()
                    .sum::<f64>()
                    > 0.0
            );
        });
    }
}
