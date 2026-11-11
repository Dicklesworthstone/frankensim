// Included inside the saturation integration module to reuse its source,
// thermal, and conservation helpers. Authored transport fractions isolate
// this new thermal adapter; the transport owner's tests exercise ray sampling.

#[test]
fn g1_finite_afterglow_routes_both_spectra_into_passive_cells_once() {
    use fs_render::volumes::heating::exchange::RadiativeThermalCell;
    use fs_render::volumes::heating::luminescence::{
        AbsorbingLuminescenceEnvironment, LuminescenceAbsorptionTolerance,
    };
    use fs_render::volumes::transport::scattering::luminescence::{
        FractionEstimate, LuminescenceAbsorption, LuminescenceTransferOperator,
        LuminescenceTransferRow, LuminescenceTransportConfig,
    };

    with_cx(|cx| {
        let density = [1.0, 1.0];
        // A scattering grid takes absorption and scattering over one domain;
        // both are authored over the scattering table's 380 to 900 nm.
        let optical = grid(
            cx,
            &density,
            MassAbsorptionSpectrum::source_table(
                &[(380.0, 1.0), (900.0, 1.0)],
                ContentHash([243; 32]),
            )
            .unwrap(),
        );
        let medium = ScatteringGrid::try_new(
            cx,
            optical,
            MassScatteringSpectrum::source_table(
                &[(380.0, 0.5), (900.0, 0.5)],
                ContentHash([245; 32]),
            )
            .unwrap(),
            PhaseFunction::Isotropic,
        )
        .unwrap();
        let thermal = thermal(
            vec![cell(1, 1.0, 1.0, 0.0), cell(2, 1.0, 1.0, 0.0)],
            &SENSIBLE,
            vec![],
        );
        let fraction = |mean| FractionEstimate {
            mean,
            standard_error: [0.0; 2],
        };
        let operator = LuminescenceTransferOperator {
            medium_identity: medium.content_identity(),
            absorption_identity: medium.absorption_grid().content_identity(),
            wavelength_bounds_nm: medium.wavelength_bounds_nm(),
            wavelengths_nm: [650.0, 800.0],
            grid_cell_count: 2,
            material_cells: vec![0, 1],
            rows: vec![LuminescenceTransferRow {
                emitter_cell: 0,
                absorption: vec![
                    LuminescenceAbsorption {
                        cell_index: 0,
                        fraction: fraction([0.1, 0.2]),
                    },
                    LuminescenceAbsorption {
                        cell_index: 1,
                        fraction: fraction([0.5, 0.3]),
                    },
                ],
                absorbed: fraction([0.6, 0.5]),
                escape: fraction([0.3, 0.5]),
                unresolved: fraction([0.1, 0.0]),
                precision_limited_histories: [0; 2],
            }],
            traced_segments: 4,
            config: LuminescenceTransportConfig {
                histories_per_emitter_lane: 2,
                max_collisions_per_history: 2,
                max_total_segments: 16,
                max_material_cells: 2,
                max_emitter_cells: 1,
                max_absorption_links: 2,
                seed: 71,
                source_id: 9,
            },
        };
        let bindings = [0, 1].map(|index| RadiativeThermalCell {
            grid_cell: index,
            thermal_cell: CellId(index as u64 + 1),
            minimum_temperature_k: 290.0,
            maximum_temperature_k: 800.0,
        });
        let environment = AbsorbingLuminescenceEnvironment::try_new(
            &medium,
            operator,
            &bindings,
            &thermal,
            2,
            LuminescenceAbsorptionTolerance {
                fraction_roundoff: 1e-13,
                max_fraction_standard_error: 0.1,
                max_unresolved_energy_j: 1.0,
            },
            || Ok(()),
        )
        .unwrap();
        let triplet = SingletTripletSource::try_new(
            source(600.0, 5.0, [0.0, 1.0], 246),
            0.0,
            ArrheniusQuenching::try_new(0.0, 0.0).unwrap(),
            SingletTripletLimits::default(),
        )
        .unwrap();
        let pump =
            PrescribedPhotonPump::try_new(ContentHash([247; 32]), [400.0; 2], [1e-20; 2], [0.0; 2])
                .unwrap();
        let heating = SaturableLuminescentHeating::try_new(
            medium.absorption_grid(),
            source(500.0, 1.0, [1.0, 0.0], 248),
            triplet,
            &[SaturableEmitterCell {
                grid_cell: 0,
                thermal_cell: CellId(1),
                initial_populations: [0.0, 0.25 * COUNT, 0.75 * COUNT],
                pump,
            }],
            &thermal,
            2,
            1e-8,
            Some(environment),
            || Ok(()),
        )
        .unwrap();
        let (next, accepted, audit) = thermal
            .advance_observed(&heating, 2.0, config(2.0), || Ok(()))
            .unwrap();
        let events = [
            0.25 * COUNT * (1.0 - (-2.0_f64).exp()),
            0.75 * COUNT * (1.0 - (-0.4_f64).exp()),
        ];
        let emitted = [events[0] * HC / 650.0, events[1] * HC / 800.0];
        let stokes = events[0] * HC * (1.0 / 500.0 - 1.0 / 650.0)
            + events[1] * HC * (1.0 / 600.0 - 1.0 / 800.0);
        for (actual, expected) in accepted.energy().emitted_j.into_iter().zip(emitted) {
            relative(actual, expected);
        }
        relative(accepted.energy().escaped_j[0], 0.3 * emitted[0]);
        relative(accepted.energy().escaped_j[1], 0.5 * emitted[1]);
        relative(
            accepted.energy().unresolved_radiation_j[0],
            0.1 * emitted[0],
        );
        assert_eq!(accepted.energy().unresolved_radiation_j[1], 0.0);
        let emitter_heat = stokes + 0.1 * emitted[0] + 0.2 * emitted[1];
        let receiver_heat = 0.5 * emitted[0] + 0.3 * emitted[1];
        relative(
            4.0 * (next.temperature_k(CellId(1)).unwrap() - 300.0),
            emitter_heat,
        );
        relative(
            4.0 * (next.temperature_k(CellId(2)).unwrap() - 300.0),
            receiver_heat,
        );
        relative(
            accepted.energy().secondary_absorption_heat_j,
            0.6 * emitted[0] + 0.5 * emitted[1],
        );
        relative(audit.external_heat_j, emitter_heat + receiver_heat);
        relative(accepted.thermal_input_j(), emitter_heat + receiver_heat);
        assert_eq!(accepted.absorbed_pump_energy_j(), [0.0; 2]);
        assert_eq!(accepted.populations(1), None);
        relative(
            accepted.populations(0).unwrap()[1],
            0.25 * COUNT * (-2.0_f64).exp(),
        );
        relative(
            accepted.populations(0).unwrap()[2],
            0.75 * COUNT * (-0.4_f64).exp(),
        );
        balance(&accepted);
    });
}

#[test]
fn g1_uniform_emitter_cells_keep_their_own_inventory_and_pump_photon_energy() {
    with_cx(|cx| {
        let density = [1.0; 3];
        let medium = grid(cx, &density, transparent_host());
        let thermal = thermal(
            (1..=3).map(|id| cell(id, 1.0, 1.0, 0.0)).collect(),
            &SENSIBLE,
            vec![],
        );
        let inventories = [COUNT, 0.5 * COUNT, 2.0 * COUNT];
        let wavelengths = [400.0, 400.0, 450.0];
        let bindings: Vec<_> = (0..3)
            .map(|index| SaturableEmitterCell {
                grid_cell: index,
                thermal_cell: CellId(index as u64 + 1),
                initial_populations: [inventories[index], 0.0, 0.0],
                pump: PrescribedPhotonPump::try_new(
                    ContentHash([249; 32]),
                    [wavelengths[index]; 2],
                    [1e-20, 0.0],
                    [9e20, 0.0],
                )
                .unwrap(),
            })
            .collect();
        let heating = SaturableLuminescentHeating::try_new(
            &medium,
            source(500.0, 1.0, [1.0, 0.0], 250),
            SingletTripletSource::try_new(
                source(600.0, 5.0, [0.0, 1.0], 251),
                0.0,
                ArrheniusQuenching::try_new(0.0, 0.0).unwrap(),
                SingletTripletLimits::default(),
            )
            .unwrap(),
            &bindings,
            &thermal,
            3,
            1e-7,
            None,
            || Ok(()),
        )
        .unwrap();
        let (next, accepted, audit) = thermal
            .advance_observed(&heating, 1.0, config(1.0), || Ok(()))
            .unwrap();
        let excited_fraction = 0.9 * (1.0 - (-10.0_f64).exp());
        let residence_per_emitter = 0.9 * (1.0 - (1.0 - (-10.0_f64).exp()) / 10.0);
        let pump_events_per_emitter = 9.0 * (1.0 - residence_per_emitter);
        let mut absorbed = 0.0;
        let mut native_heat = 0.0;
        for index in 0..3 {
            let count = inventories[index];
            let populations = accepted.populations(index).unwrap();
            relative(populations.iter().sum(), count);
            relative(populations[1], count * excited_fraction);
            let pump_energy = count * pump_events_per_emitter * HC / wavelengths[index];
            let emitted = count * residence_per_emitter * HC / 650.0;
            let stored = count * excited_fraction * HC / 500.0;
            let heat = pump_energy - emitted - stored;
            relative(
                4.0 * (next.temperature_k(CellId(index as u64 + 1)).unwrap() - 300.0),
                heat,
            );
            absorbed += pump_energy;
            native_heat += heat;
        }
        relative(accepted.absorbed_pump_energy_j()[0], absorbed);
        relative(accepted.thermal_input_j(), native_heat);
        relative(audit.external_heat_j, native_heat);
        relative(
            absorbed,
            accepted.energy().stored_excitation_j
                + accepted.energy().escaped_j.iter().sum::<f64>()
                + native_heat,
        );
    });
}
