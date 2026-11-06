use super::*;
use fs_blake3::ContentHash;
use fs_render::volumes::heating::exchange::RadiativeThermalCell;
use fs_render::volumes::heating::fluorescence::FluorescenceSource;
use fs_render::volumes::heating::illumination::schedule::{
    SpectralPowerSchedule, SpectralPowerSegment,
};
use fs_render::volumes::heating::luminescence::{
    AbsorbingLuminescenceEnvironment, LuminescenceAbsorptionTolerance, LuminescenceSource,
    LuminescentOpticalHeating,
};
use fs_render::volumes::transport::scattering::luminescence::{
    FractionEstimate, LuminescenceAbsorption, LuminescenceTransferOperator,
    LuminescenceTransferRow, LuminescenceTransportConfig,
};

fn source() -> LuminescenceSource<1> {
    LuminescenceSource::try_new(
        FluorescenceSource::try_new(
            ContentHash([122; 32]),
            (380.0, 500.0),
            (290.0, 800.0),
            0.8,
            [800.0],
            [1.0],
        )
        .unwrap(),
        500.0,
        2.0,
    )
    .unwrap()
}

fn medium<'a>(cx: &Cx<'_>, density: &'a [f64]) -> ScatteringGrid<'a> {
    ScatteringGrid::try_new(
        cx,
        grid(
            cx,
            density,
            MassAbsorptionSpectrum::source_table(
                &[(380.0, 1e4), (900.0, 1e4)],
                ContentHash([123; 32]),
            )
            .unwrap(),
        ),
        MassScatteringSpectrum::source_table(&[(380.0, 0.0), (900.0, 0.0)], ContentHash([124; 32]))
            .unwrap(),
        PhaseFunction::Isotropic,
    )
    .unwrap()
}

fn sampling() -> LuminescenceTransportConfig {
    LuminescenceTransportConfig {
        histories_per_emitter_lane: 2,
        max_collisions_per_history: 2,
        max_total_segments: 16,
        max_material_cells: 2,
        max_emitter_cells: 1,
        max_absorption_links: 2,
        seed: 71,
        source_id: 9,
    }
}

fn fraction(value: f64) -> FractionEstimate<1> {
    FractionEstimate {
        mean: [value],
        standard_error: [0.0],
    }
}

/// Authored probability fixture isolates the thermal adapter from ray sampling.
/// The separate transport tests exercise actual grid-derived probabilities.
fn operator(medium: &ScatteringGrid<'_>, tail: f64) -> LuminescenceTransferOperator<1> {
    LuminescenceTransferOperator {
        medium_identity: medium.content_identity(),
        absorption_identity: medium.absorption_grid().content_identity(),
        wavelength_bounds_nm: medium.wavelength_bounds_nm(),
        wavelengths_nm: [800.0],
        grid_cell_count: 2,
        material_cells: vec![0, 1],
        rows: vec![LuminescenceTransferRow {
            emitter_cell: 0,
            absorption: vec![
                LuminescenceAbsorption {
                    cell_index: 0,
                    fraction: fraction(0.25),
                },
                LuminescenceAbsorption {
                    cell_index: 1,
                    fraction: fraction(0.5 - tail),
                },
            ],
            absorbed: fraction(0.75 - tail),
            escape: fraction(0.25),
            unresolved: fraction(tail),
            precision_limited_histories: [0],
        }],
        traced_segments: 2,
        config: sampling(),
    }
}

fn bindings(receiver_max_k: f64) -> [RadiativeThermalCell; 2] {
    [
        RadiativeThermalCell {
            grid_cell: 0,
            thermal_cell: CellId(1),
            minimum_temperature_k: 290.0,
            maximum_temperature_k: 800.0,
        },
        RadiativeThermalCell {
            grid_cell: 1,
            thermal_cell: CellId(2),
            minimum_temperature_k: 290.0,
            maximum_temperature_k: receiver_max_k,
        },
    ]
}

fn tolerance(tail_j: f64) -> LuminescenceAbsorptionTolerance {
    LuminescenceAbsorptionTolerance {
        fraction_roundoff: 1e-13,
        max_fraction_standard_error: 0.1,
        max_unresolved_energy_j: tail_j,
    }
}

fn setup(
    cx: &Cx<'_>,
    tail: f64,
    allowed_tail_j: f64,
    receiver_max_k: f64,
) -> (
    DepositionSimulation<'static>,
    LuminescentOpticalHeating<1, 1>,
) {
    let density = [1.0, 1.0];
    let medium = medium(cx, &density);
    let optical = medium.absorption_grid();
    let deposition = optical
        .deposit_beam_power(
            cx,
            Point3::new(0.0, 0.5, 0.5),
            Vec3::new(1.0, 0.0, 0.0),
            1.0,
            &[400.0],
            &[10.0],
            2,
        )
        .unwrap();
    let thermal = thermal(
        vec![cell(1, 1.0, 1.0, 0.0), cell(2, 1.0, 1.0, 0.0)],
        &SENSIBLE,
        vec![],
    );
    let environment = AbsorbingLuminescenceEnvironment::try_new(
        &medium,
        operator(&medium, tail),
        &bindings(receiver_max_k),
        &thermal,
        2,
        tolerance(allowed_tail_j),
        || Ok(()),
    )
    .unwrap();
    let heating = LuminescentOpticalHeating::try_new_with_emission_transport(
        &medium,
        source(),
        deposition,
        &[(0, CellId(1))],
        &thermal,
        2,
        1e-12,
        1e-9,
        environment,
        || Ok(()),
    )
    .unwrap();
    (thermal, heating)
}

fn partition(heating: &LuminescentOpticalHeating<1, 1>) {
    let energy = heating.energy();
    let heat = energy.prompt_heat_j
        + energy.nonradiative_heat_j
        + energy.relaxation_heat_j
        + energy.secondary_absorption_heat_j;
    close(heat, heating.pump_energy().thermal_input_j);
    close(
        heating.pump_energy().absorbed_j[0],
        heat + energy.stored_excitation_j + energy.escaped_j[0] + energy.unresolved_radiation_j[0],
    );
    close(
        energy.emitted_j[0],
        energy.secondary_absorption_heat_j + energy.escaped_j[0] + energy.unresolved_radiation_j[0],
    );
    assert!(energy.residual_j.abs() < 1e-11);
}

#[test]
fn g1_self_return_and_an_unpumped_receiver_receive_integrated_emitted_energy_once() {
    with_cx(|cx| {
        let (thermal, heating) = setup(cx, 0.0, 0.0, 800.0);
        let (next, accepted, audit) = thermal
            .advance_observed(&heating, 1.0, config(1.0), || Ok(()))
            .unwrap();
        let stored = 16.0 * (1.0 - (-0.5_f64).exp());
        let released = 8.0 - stored;
        let emitted = 0.5 * released;
        close(accepted.energy().stored_excitation_j, stored);
        close(accepted.energy().emitted_j[0], emitted);
        close(accepted.energy().escaped_j[0], 0.25 * emitted);
        close(
            accepted.energy().secondary_absorption_heat_j,
            0.75 * emitted,
        );
        close(
            next.temperature_k(CellId(1)).unwrap(),
            300.0 + (2.0 + 0.625 * released) / 4.0,
        );
        close(
            next.temperature_k(CellId(2)).unwrap(),
            300.0 + 0.25 * released / 4.0,
        );
        close(
            audit.external_heat_j,
            accepted.pump_energy().thermal_input_j,
        );
        close(accepted.last_mean_cell_emission_w()[0].1[0], emitted);
        // Endpoint power is larger while pumping; it must not drive this heat.
        assert!(accepted.cell_emission_w()[0].1[0] > emitted);
        partition(&accepted);
        assert_eq!(heating.energy().emitted_j, [0.0]);
        assert_eq!(thermal.temperature_k(CellId(2)).unwrap(), 300.0);
    });
}

#[test]
fn g1_dark_afterglow_keeps_heating_a_remote_absorber_and_split_steps_agree() {
    with_cx(|cx| {
        let (thermal, heating) = setup(cx, 0.0, 0.0, 800.0);
        let heating = heating
            .with_pump_schedule(
                SpectralPowerSchedule::try_new(
                    &[SpectralPowerSegment {
                        start_s: 0.0,
                        end_s: 1.0,
                        start_power_w: [10.0],
                        end_power_w: [10.0],
                    }],
                    1,
                )
                .unwrap(),
            )
            .unwrap();
        let (whole, all, audit) = thermal
            .advance_observed(&heating, 8.0, config(8.0), || Ok(()))
            .unwrap();
        let (at_one, first, _) = thermal
            .advance_observed(&heating, 1.0, config(1.0), || Ok(()))
            .unwrap();
        let (split, second, _) = at_one
            .advance_observed(&first, 8.0, config(7.0), || Ok(()))
            .unwrap();
        assert_eq!(audit.substeps, 2);
        assert!(whole.temperature_k(CellId(2)).unwrap() > at_one.temperature_k(CellId(2)).unwrap());
        close(all.energy().emitted_j[0], second.energy().emitted_j[0]);
        close(all.energy().escaped_j[0], second.energy().escaped_j[0]);
        close(
            all.energy().stored_excitation_j,
            second.energy().stored_excitation_j,
        );
        close(
            whole.temperature_k(CellId(2)).unwrap(),
            split.temperature_k(CellId(2)).unwrap(),
        );
        close(all.pump_energy().incident_j[0], 10.0);
        partition(&all);
    });
}

#[test]
fn g1_unresolved_radiation_is_neither_heat_nor_escape_and_its_limit_is_transactional() {
    with_cx(|cx| {
        let (thermal, allowed) = setup(cx, 0.25, 10.0, 800.0);
        let (_, accepted, _) = thermal
            .advance_observed(&allowed, 1.0, config(1.0), || Ok(()))
            .unwrap();
        close(
            accepted.energy().unresolved_radiation_j[0],
            0.25 * accepted.energy().emitted_j[0],
        );
        partition(&accepted);
        let (thermal, strict) = setup(cx, 0.25, 0.1, 800.0);
        let (short, initial, _) = thermal
            .advance_observed(&strict, 0.1, config(0.1), || Ok(()))
            .unwrap();
        assert!(initial.energy().unresolved_radiation_j[0] > 0.0);
        assert!(
            short
                .advance_observed(&initial, 1.0, config(1.0), || Ok(()))
                .is_err()
        );
        assert_eq!(initial.time_s(), 0.1);
        assert_eq!(short.time_s(), 0.1);
        assert!(initial.energy().unresolved_radiation_j[0] < 0.1);
    });
}

#[test]
fn g0_passive_receiver_temperature_domain_is_checked_at_the_accepted_endpoint() {
    with_cx(|cx| {
        let (thermal, heating) = setup(cx, 0.0, 0.0, 300.01);
        assert!(
            thermal
                .advance_observed(&heating, 1.0, config(1.0), || Ok(()))
                .is_err()
        );
        assert_eq!(thermal.time_s(), 0.0);
        assert_eq!(heating.energy().secondary_absorption_heat_j, 0.0);
        assert_eq!(heating.stored_energy_j(0), Some(0.0));
    });
}

#[test]
fn g1_actual_opaque_cell_transport_recovers_luminescent_energy_as_native_heat() {
    with_cx(|cx| {
        let density = [1.0];
        let medium = medium(cx, &density);
        let optical = medium.absorption_grid();
        let transfer = medium
            .estimate_luminescence_transfer(cx, &[800.0], &[0], sampling())
            .unwrap();
        assert!(transfer.rows[0].absorbed.mean[0] > 0.0);
        let absorbed_fraction = transfer.rows[0].absorbed.mean[0];
        let deposition = optical
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
        assert!(
            LuminescentOpticalHeating::try_new(
                optical,
                source(),
                deposition.clone(),
                &[(0, CellId(1))],
                &thermal,
                1,
                1e-12,
                1e-9,
                || Ok(())
            )
            .is_err()
        );
        let environment = AbsorbingLuminescenceEnvironment::try_new(
            &medium,
            transfer,
            &bindings(800.0)[..1],
            &thermal,
            1,
            tolerance(0.0),
            || Ok(()),
        )
        .unwrap();
        let heating = LuminescentOpticalHeating::try_new_with_emission_transport(
            &medium,
            source(),
            deposition,
            &[(0, CellId(1))],
            &thermal,
            1,
            1e-12,
            1e-9,
            environment,
            || Ok(()),
        )
        .unwrap();
        let (_, accepted, audit) = thermal
            .advance_observed(&heating, 1.0, config(1.0), || Ok(()))
            .unwrap();
        close(
            accepted.energy().secondary_absorption_heat_j,
            absorbed_fraction * accepted.energy().emitted_j[0],
        );
        close(
            audit.external_heat_j,
            10.0 - accepted.energy().stored_excitation_j - accepted.energy().escaped_j[0],
        );
        partition(&accepted);
    });
}

#[test]
fn g0_incomplete_bindings_mutated_rows_and_unsupported_sampling_refuse() {
    with_cx(|cx| {
        let density = [1.0, 1.0];
        let medium = medium(cx, &density);
        let thermal = thermal(
            vec![cell(1, 1.0, 1.0, 0.0), cell(2, 1.0, 1.0, 0.0)],
            &SENSIBLE,
            vec![],
        );
        let validate = |report, bound: &[RadiativeThermalCell]| {
            AbsorbingLuminescenceEnvironment::try_new(
                &medium,
                report,
                bound,
                &thermal,
                2,
                tolerance(0.0),
                || Ok(()),
            )
        };
        assert!(validate(operator(&medium, 0.0), &bindings(800.0)[..1]).is_err());
        let mut mutated = operator(&medium, 0.0);
        mutated.rows[0].escape.mean[0] = 0.9;
        assert!(validate(mutated, &bindings(800.0)).is_err());
        let mut imprecise = operator(&medium, 0.0);
        imprecise.rows[0].escape.standard_error[0] = 0.5;
        assert!(validate(imprecise, &bindings(800.0)).is_err());
        let mut mismatched = operator(&medium, 0.0);
        mismatched.medium_identity = ContentHash([233; 32]);
        assert!(validate(mismatched, &bindings(800.0)).is_err());
        let mut omitted = operator(&medium, 0.0);
        omitted.material_cells.pop();
        assert!(validate(omitted, &bindings(800.0)).is_err());
    });
}

#[test]
fn g0_transparent_luminescence_rejects_an_out_of_grid_deposition_with_a_valid_thermal_binding() {
    with_cx(|cx| {
        let density = [1.0];
        let optical = grid(
            cx,
            &density,
            MassAbsorptionSpectrum::source_table(
                &[(380.0, 1e4), (500.0, 1e4), (800.0, 0.0), (900.0, 0.0)],
                ContentHash([125; 32]),
            )
            .unwrap(),
        );
        let deposition = optical
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
        assert!(
            LuminescentOpticalHeating::try_new(
                &optical,
                source(),
                deposition.clone(),
                &[(0, CellId(1))],
                &thermal,
                1,
                1e-12,
                1e-9,
                || Ok(()),
            )
            .is_ok()
        );

        // Keep the original identity and conserved powers, and supply an
        // active thermal receiver: only the optical cell index is invalid.
        let mut outside = deposition;
        assert_eq!(outside.cells.len(), 1);
        outside.cells[0].cell_index = optical.cell_count();
        assert!(matches!(
            LuminescentOpticalHeating::try_new(
                &optical,
                source(),
                outside,
                &[(optical.cell_count(), CellId(1))],
                &thermal,
                1,
                1e-12,
                1e-9,
                || Ok(()),
            ),
            Err(OpticalHeatingError::Invalid(_))
        ));
    });
}

#[test]
fn g0_emission_scattering_is_allowed_but_hidden_pump_scattering_refuses() {
    with_cx(|cx| {
        let density = [1.0, 1.0];
        let absorption = MassAbsorptionSpectrum::source_table(
            &[(380.0, 1e4), (900.0, 1e4)],
            ContentHash([126; 32]),
        )
        .unwrap();
        let optical = grid(cx, &density, absorption);
        let deposition = optical
            .deposit_beam_power(
                cx,
                Point3::new(0.0, 0.5, 0.5),
                Vec3::new(1.0, 0.0, 0.0),
                1.0,
                &[400.0],
                &[10.0],
                2,
            )
            .unwrap();
        let thermal = thermal(
            vec![cell(1, 1.0, 1.0, 0.0), cell(2, 1.0, 1.0, 0.0)],
            &SENSIBLE,
            vec![],
        );
        assert_eq!(deposition.scattered_power_w, [0.0]);

        for pump_scattering in [0.0, 1.0] {
            let medium = ScatteringGrid::try_new(
                cx,
                grid(cx, &density, absorption),
                MassScatteringSpectrum::source_table(
                    &[
                        (380.0, pump_scattering),
                        (500.0, pump_scattering),
                        (800.0, 2.0),
                        (900.0, 2.0),
                    ],
                    ContentHash([127; 32]),
                )
                .unwrap(),
                PhaseFunction::Isotropic,
            )
            .unwrap();
            assert_eq!(
                medium.absorption_grid().content_identity(),
                deposition.medium_identity,
            );
            assert_eq!(
                medium.scattering_at_cell(0).unwrap().at(800.0).unwrap(),
                2.0
            );
            // Each environment matches its full snapshot. Rejection of the
            // second case must not rely on a mismatched operator identity.
            let environment = AbsorbingLuminescenceEnvironment::try_new(
                &medium,
                operator(&medium, 0.0),
                &bindings(800.0),
                &thermal,
                2,
                tolerance(0.0),
                || Ok(()),
            )
            .unwrap();
            let result = LuminescentOpticalHeating::try_new_with_emission_transport(
                &medium,
                source(),
                deposition.clone(),
                &[(0, CellId(1))],
                &thermal,
                2,
                1e-12,
                1e-9,
                environment,
                || Ok(()),
            );
            if pump_scattering == 0.0 {
                assert!(result.is_ok(), "emission-band scattering is supported");
            } else {
                assert!(matches!(result, Err(OpticalHeatingError::Invalid(_))));
            }
        }
    });
}

#[test]
fn g3_preparation_preserves_native_cancellation_and_budget_refusals() {
    with_cx(|cx| {
        let density = [1.0, 1.0];
        let medium = medium(cx, &density);
        let thermal = thermal(
            vec![cell(1, 1.0, 1.0, 0.0), cell(2, 1.0, 1.0, 0.0)],
            &SENSIBLE,
            vec![],
        );
        let mut calls = 0;
        let error = AbsorbingLuminescenceEnvironment::try_new(
            &medium,
            operator(&medium, 0.0),
            &bindings(800.0),
            &thermal,
            2,
            tolerance(0.0),
            || {
                calls += 1;
                if calls == 3 {
                    Err(DepositionError::Cancelled)
                } else {
                    Ok(())
                }
            },
        )
        .unwrap_err();
        assert!(matches!(
            error,
            OpticalHeatingError::Thermal(DepositionError::Cancelled)
        ));
        assert!(
            AbsorbingLuminescenceEnvironment::try_new(
                &medium,
                operator(&medium, 0.0),
                &bindings(800.0),
                &thermal,
                1,
                tolerance(0.0),
                || Ok(())
            )
            .is_err()
        );
    });
}
