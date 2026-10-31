use super::*;
use fs_blake3::ContentHash;
use fs_render::volumes::heating::exchange::{
    RadiativeBath, RadiativeExchangeBudget, RadiativeExchangeError, RadiativeExchangeTolerance,
    RadiativeThermalCell, VolumeRadiativeExchange,
};
use fs_render::volumes::transport::scattering::exchange::{
    RadiativeExchangeConfig, RadiativeExchangeOperator,
};

struct Sensible;
impl ThermalLaw for Sensible {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        if !h.is_finite() || !(10.0..=20_000.0).contains(&h) {
            return Err(DepositionError::MaterialDomain);
        }
        Ok(h / 10.0)
    }
    fn maximum_temperature_slope(&self) -> f64 {
        0.1
    }
}
static RADIATIVE_SENSIBLE: Sensible = Sensible;

fn warm_cell(id: u64, temperature_k: f64) -> DepositedCell {
    let mut value = cell(id, 1.0, 1.0, 0.0);
    value.specific_h_j_kg = 10.0 * temperature_k;
    value
}
fn photon_config() -> RadiativeExchangeConfig {
    RadiativeExchangeConfig {
        histories_per_emitter_lane: 1024,
        max_collisions_per_history: 96,
        max_total_segments: 1_000_000,
        max_emitter_cells: 8,
        max_material_cells: 8,
        max_pairs: 16,
        seed: 48,
        source_id: 518,
    }
}
fn operator(
    cx: &Cx<'_>,
    density: &[f64],
    absorption: f64,
    scattering: f64,
    config: RadiativeExchangeConfig,
) -> RadiativeExchangeOperator<1> {
    let absorption = MassAbsorptionSpectrum::source_table(
        &[(8000.0, absorption), (18000.0, absorption)],
        ContentHash([81; 32]),
    )
    .unwrap();
    let scattering = MassScatteringSpectrum::source_table(
        &[(8000.0, scattering), (18000.0, scattering)],
        ContentHash([82; 32]),
    )
    .unwrap();
    ScatteringGrid::try_new(
        cx,
        grid(cx, density, absorption),
        scattering,
        PhaseFunction::Isotropic,
    )
    .unwrap()
    .estimate_thermal_exchange(cx, &[10_000.0], config)
    .unwrap()
}
fn binding(index: usize, id: u64) -> RadiativeThermalCell {
    RadiativeThermalCell {
        grid_cell: index,
        thermal_cell: CellId(id),
        minimum_temperature_k: 1.0,
        maximum_temperature_k: 1000.0,
    }
}
fn budget() -> RadiativeExchangeBudget {
    RadiativeExchangeBudget {
        max_bindings: 8,
        max_thermal_pairs: 8,
        max_source_evaluations: 1000,
    }
}
fn tolerance() -> RadiativeExchangeTolerance {
    RadiativeExchangeTolerance {
        coefficient_roundoff_m2_sr: 1e-10,
        power_roundoff_w: 1e-9,
        max_cell_standard_error_w: 1e4,
        max_total_standard_error_w: 1e4,
        max_unresolved_exchange_w: 1.0,
        max_unresolved_standard_error_w: 1.0,
    }
}
fn blackbody(temperature_k: f64) -> RadiativeBath {
    RadiativeBath::Blackbody { temperature_k }
}
fn observer(
    op: RadiativeExchangeOperator<1>,
    weight_nm: f64,
    bath: RadiativeBath,
    bindings: &[RadiativeThermalCell],
    thermal: &DepositionSimulation<'_>,
) -> VolumeRadiativeExchange<1> {
    VolumeRadiativeExchange::try_new(
        op,
        [weight_nm],
        bath,
        bindings,
        thermal,
        budget(),
        tolerance(),
        || Ok(()),
    )
    .unwrap()
}
// Independent direct SI Planck expression, used only at ordinary finite states.
fn planck(temperature_k: f64) -> f64 {
    let wavelength_m: f64 = 1e-5;
    let h = 6.626_070_15e-34;
    let c = 299_792_458.0;
    let k = 1.380_649e-23;
    2.0 * h * c * c / wavelength_m.powi(5) / (h * c / (wavelength_m * k * temperature_k)).exp_m1()
        * 1e-9
}
fn relative_close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 1e-10 * expected.abs().max(1.0),
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn sampled_reciprocal_field_has_exact_uniform_thermal_equilibrium() {
    with_cx(|cx| {
        let original = thermal(
            vec![warm_cell(11, 300.0), warm_cell(12, 300.0)],
            &RADIATIVE_SENSIBLE,
            vec![],
        );
        let op = operator(cx, &[1.0, 0.3], 0.7, 1.0, photon_config());
        assert!(!op.pairs.is_empty());
        assert!(
            op.cells
                .iter()
                .any(|c| c.escape.standard_error_m2_sr[0] > 0.0)
        );
        let heating = VolumeRadiativeExchange::try_new(
            op,
            [3000.0],
            blackbody(300.0),
            &[binding(0, 11), binding(1, 12)],
            &original,
            budget(),
            RadiativeExchangeTolerance {
                max_cell_standard_error_w: 0.0,
                max_total_standard_error_w: 0.0,
                max_unresolved_exchange_w: 0.0,
                max_unresolved_standard_error_w: 0.0,
                ..tolerance()
            },
            || Ok(()),
        )
        .unwrap();
        assert!(heating.sources_w().is_empty());
        assert!(heating.transfers().is_empty());
        assert_eq!(heating.power().net_bath_to_material_w, [0.0]);
        assert_eq!(heating.power().total_standard_error_w, 0.0);
        assert!(heating.conductance_bounds_w_k().values().all(|g| *g > 0.0));
        let (accepted, radiative, audit) = original
            .advance_observed(&heating, 2.0, config(1.0), || Ok(()))
            .unwrap();
        assert_eq!(accepted.cells(), original.cells());
        assert_eq!(audit.external_heat_j, 0.0);
        assert_eq!(audit.stored_energy_change_j, 0.0);
        assert_eq!(radiative.energy().thermal_external_j, 0.0);
        assert_eq!(radiative.energy().internal_transfer_magnitude_j, [0.0]);
    });
}

#[test]
fn hotter_volume_cools_and_reabsorbed_heat_reaches_original_cold_cell() {
    with_cx(|cx| {
        let original = thermal(
            vec![warm_cell(11, 500.0), warm_cell(12, 300.0)],
            &RADIATIVE_SENSIBLE,
            vec![],
        );
        let op = operator(cx, &[1.0; 2], 0.7, 1.0, photon_config());
        let coefficient = op
            .pairs
            .iter()
            .find(|p| (p.cell_a, p.cell_b) == (0, 1))
            .unwrap()
            .exchange
            .mean_m2_sr[0];
        let escaped_hot = op.cells[0].escape.mean_m2_sr[0];
        let weight = 3000.0;
        let expected_pair_w = coefficient * weight * (planck(500.0) - planck(300.0));
        let expected_bath_w = -escaped_hot * weight * (planck(500.0) - planck(300.0));
        let exchange = observer(
            op,
            weight,
            blackbody(300.0),
            &[binding(0, 11), binding(1, 12)],
            &original,
        );
        assert_eq!(exchange.transfers().len(), 1);
        relative_close(exchange.transfers()[0].power_w, expected_pair_w);
        relative_close(exchange.sources_w()[&CellId(11)], expected_bath_w);
        assert!(expected_pair_w > 0.0 && expected_bath_w < 0.0);
        let duration = 1e-4;
        let (heated, accepted, audit) = original
            .advance_observed(&exchange, duration, config(1.0), || Ok(()))
            .unwrap();
        assert_eq!(audit.substeps, 1);
        assert!(heated.temperature_k(CellId(11)).unwrap() < 500.0);
        assert!(heated.temperature_k(CellId(12)).unwrap() > 300.0);
        relative_close(
            heated.cells()[&CellId(11)].specific_h_j_kg,
            5000.0 + duration * (expected_bath_w - expected_pair_w),
        );
        relative_close(
            heated.cells()[&CellId(12)].specific_h_j_kg,
            3000.0 + duration * expected_pair_w,
        );
        relative_close(audit.external_heat_j, duration * expected_bath_w);
        relative_close(audit.external_heat_j, accepted.energy().thermal_external_j);
        relative_close(audit.stored_energy_change_j, audit.external_heat_j);
        relative_close(
            accepted.energy().internal_transfer_magnitude_j[0],
            duration * expected_pair_w,
        );
        assert!(accepted.power().total_standard_error_w > 0.0);
    });
}

#[test]
fn native_radiative_limiter_prevents_stiff_vacuum_cooling_overshoot() {
    with_cx(|cx| {
        let original = thermal(vec![warm_cell(11, 500.0)], &RADIATIVE_SENSIBLE, vec![]);
        let op = operator(cx, &[1.0], 2.0, 0.5, photon_config());
        let exchange = observer(
            op,
            10_000.0,
            RadiativeBath::Vacuum,
            &[binding(0, 11)],
            &original,
        );
        let initial_cooling_w = -exchange.power().thermal_external_w;
        assert!(initial_cooling_w > 0.0);
        // This interval would remove more than the entire initial thermal store
        // if the original sampled watts were applied without the native limiter.
        let duration = 1.5 * 5000.0 / initial_cooling_w;
        let (cooled, accepted, audit) = original
            .advance_observed(
                &exchange,
                duration,
                StepConfig {
                    max_dt_s: duration * 10.0,
                    energy_tolerance_j: 1e-8,
                },
                || Ok(()),
            )
            .unwrap();
        assert!(audit.substeps > 1);
        let final_t = cooled.temperature_k(CellId(11)).unwrap();
        assert!(final_t > 1.0 && final_t < 500.0);
        assert!(accepted.energy().thermal_external_j < 0.0);
        relative_close(audit.external_heat_j, accepted.energy().thermal_external_j);
        relative_close(audit.stored_energy_change_j, audit.external_heat_j);
    });
}

#[test]
fn black_bath_heats_without_crossing_its_temperature() {
    with_cx(|cx| {
        let original = thermal(vec![warm_cell(11, 100.0)], &RADIATIVE_SENSIBLE, vec![]);
        let op = operator(cx, &[1.0], 2.0, 0.5, photon_config());
        let exchange = observer(op, 10_000.0, blackbody(500.0), &[binding(0, 11)], &original);
        let duration = 2.0 * (5000.0 - 1000.0) / exchange.power().thermal_external_w;
        let (heated, accepted, audit) = original
            .advance_observed(
                &exchange,
                duration,
                StepConfig {
                    max_dt_s: duration * 10.0,
                    energy_tolerance_j: 1e-8,
                },
                || Ok(()),
            )
            .unwrap();
        assert!(audit.substeps > 1);
        let final_t = heated.temperature_k(CellId(11)).unwrap();
        assert!(final_t > 100.0 && final_t <= 500.0);
        assert!(accepted.energy().thermal_external_j > 0.0);
        relative_close(audit.external_heat_j, accepted.energy().thermal_external_j);
    });
}

#[test]
fn mapping_intersects_optical_domains_and_cancels_same_receiver_exchange() {
    with_cx(|cx| {
        let original = thermal(vec![warm_cell(11, 400.0)], &RADIATIVE_SENSIBLE, vec![]);
        let op = operator(cx, &[1.0; 2], 0.7, 1.0, photon_config());
        let mut bindings = [binding(0, 11), binding(1, 11)];
        bindings[0].minimum_temperature_k = 300.0;
        bindings[0].maximum_temperature_k = 600.0;
        bindings[1].minimum_temperature_k = 350.0;
        bindings[1].maximum_temperature_k = 500.0;
        let exchange = observer(op.clone(), 3000.0, blackbody(300.0), &bindings, &original);
        assert!(exchange.transfers().is_empty());
        let expected = op.cells.iter().map(|c| c.escape.mean_m2_sr[0]).sum::<f64>()
            * 3000.0
            * (planck(300.0) - planck(400.0));
        relative_close(exchange.power().thermal_external_w, expected);
        bindings[1].minimum_temperature_k = 650.0;
        bindings[1].maximum_temperature_k = 700.0;
        assert!(matches!(
            VolumeRadiativeExchange::try_new(
                op,
                [3000.0],
                blackbody(300.0),
                &bindings,
                &original,
                budget(),
                tolerance(),
                || Ok(())
            ),
            Err(RadiativeExchangeError::Invalid(_))
        ));
    });
}

#[test]
fn radiative_cooling_preserves_native_latent_events_and_split_history() {
    with_cx(|cx| {
        let mut initial = cell(11, 1.0, 1.0, 0.0);
        initial.specific_h_j_kg = 120.0;
        let original = thermal(
            vec![initial],
            &LATENT,
            vec![TimedEvent {
                time_s: 1.0,
                event: Event::Boundary(CellId(11), boundary(4.0)),
            }],
        );
        let op = operator(cx, &[1.0], 0.7, 1.0, photon_config());
        let exchange = observer(op, 3000.0, blackbody(300.0), &[binding(0, 11)], &original);
        let (whole, whole_exchange, audit) = original
            .advance_observed(&exchange, 2.0, config(0.125), || Ok(()))
            .unwrap();
        let (middle, middle_exchange, _) = original
            .advance_observed(&exchange, 1.0, config(0.125), || Ok(()))
            .unwrap();
        let (split, split_exchange, _) = middle
            .advance_observed(&middle_exchange, 2.0, config(0.125), || Ok(()))
            .unwrap();
        assert_eq!(whole.temperature_k(CellId(11)).unwrap(), 310.0);
        assert!((10.0..110.0).contains(&whole.cells()[&CellId(11)].specific_h_j_kg));
        assert_eq!(whole.cells(), split.cells());
        assert_eq!(whole_exchange.energy(), split_exchange.energy());
        assert_eq!(audit.events, 1);
        relative_close(
            audit.external_heat_j,
            whole_exchange.energy().thermal_external_j + 4.0,
        );
    });
}

#[test]
fn domains_sampling_and_unresolved_radiation_refuse_before_partial_heat() {
    with_cx(|cx| {
        let mut hot = warm_cell(11, 400.0);
        hot.boundary.source_w = 800.0;
        let original = thermal(vec![hot], &RADIATIVE_SENSIBLE, vec![]);
        let op = operator(cx, &[1.0], 0.7, 1.0, photon_config());
        let binding = RadiativeThermalCell {
            minimum_temperature_k: 395.0,
            maximum_temperature_k: 405.0,
            ..binding(0, 11)
        };
        let exchange = observer(op.clone(), 3000.0, blackbody(400.0), &[binding], &original);
        assert!(matches!(
            original.advance_observed(&exchange, 0.1, config(1.0), || Ok(())),
            Err(RadiativeExchangeError::TemperatureDomain {
                thermal_cell: CellId(11),
                ..
            })
        ));
        for strict in [
            RadiativeExchangeTolerance {
                max_cell_standard_error_w: 0.0,
                ..tolerance()
            },
            RadiativeExchangeTolerance {
                max_total_standard_error_w: 0.0,
                ..tolerance()
            },
        ] {
            assert!(matches!(
                VolumeRadiativeExchange::try_new(
                    op.clone(),
                    [3000.0],
                    blackbody(300.0),
                    &[binding],
                    &original,
                    budget(),
                    strict,
                    || Ok(())
                ),
                Err(RadiativeExchangeError::Invalid(_))
            ));
        }
        let truncated = operator(
            cx,
            &[1.0],
            0.05,
            10.0,
            RadiativeExchangeConfig {
                max_collisions_per_history: 1,
                ..photon_config()
            },
        );
        assert!(truncated.cells[0].unresolved.mean_m2_sr[0] > 0.0);
        assert!(matches!(
            VolumeRadiativeExchange::try_new(
                truncated,
                [3000.0],
                blackbody(300.0),
                &[binding],
                &original,
                budget(),
                RadiativeExchangeTolerance {
                    max_unresolved_exchange_w: 0.0,
                    ..tolerance()
                },
                || Ok(())
            ),
            Err(RadiativeExchangeError::Invalid(_))
        ));
        assert_eq!(original.time_s(), 0.0);
        assert_eq!(original.temperature_k(CellId(11)).unwrap(), 400.0);
        assert_eq!(exchange.energy().thermal_external_j, 0.0);
    });
}

#[test]
fn exchange_source_budget_cancellation_and_report_admission_are_transactional() {
    with_cx(|cx| {
        let empty = thermal(vec![], &RADIATIVE_SENSIBLE, vec![]);
        let vacuum = operator(cx, &[0.0], 0.7, 1.0, photon_config());
        assert!(matches!(
            VolumeRadiativeExchange::try_new(
                vacuum.clone(),
                [3000.0],
                blackbody(300.0),
                &[],
                &empty,
                budget(),
                tolerance(),
                || Ok(())
            ),
            Err(RadiativeExchangeError::Invalid(_))
        ));
        let original = thermal(vec![warm_cell(11, 400.0)], &RADIATIVE_SENSIBLE, vec![]);
        assert!(
            VolumeRadiativeExchange::try_new(
                vacuum,
                [3000.0],
                blackbody(300.0),
                &[],
                &original,
                budget(),
                tolerance(),
                || Ok(())
            )
            .is_ok()
        );
        let op = operator(cx, &[1.0], 0.7, 1.0, photon_config());
        let exchange = observer(
            op.clone(),
            3000.0,
            blackbody(300.0),
            &[binding(0, 11)],
            &original,
        );
        let limited = exchange.with_work_limit(1);
        assert!(matches!(
            original.advance_observed(&limited, 1.0, config(0.125), || Ok(())),
            Err(RadiativeExchangeError::Budget)
        ));
        let mut calls = 0;
        let result = original.advance_observed(&exchange, 1.0, config(0.125), || {
            calls += 1;
            if calls == 20 {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        });
        assert!(matches!(
            result,
            Err(RadiativeExchangeError::Thermal(DepositionError::Cancelled))
        ));
        let mut forged = op;
        forged.cells[0].escape.mean_m2_sr[0] += 1.0;
        assert!(matches!(
            VolumeRadiativeExchange::try_new(
                forged,
                [3000.0],
                blackbody(300.0),
                &[binding(0, 11)],
                &original,
                budget(),
                tolerance(),
                || Ok(())
            ),
            Err(RadiativeExchangeError::Transport(_))
        ));
        assert_eq!(original.time_s(), 0.0);
        assert_eq!(exchange.time_s(), 0.0);
        assert_eq!(exchange.energy().thermal_external_j, 0.0);
    });
}

#[test]
fn separately_heated_transparent_and_scattering_materials_keep_optical_domains() {
    use fs_render::volumes::transport::scattering::{CellPhaseFunction, CellScatteringGrid};
    with_cx(|cx| {
        for scattering_coefficient in [0.0, 2.0] {
            let density = [1.0; 2];
            let spectrum = MassAbsorptionSpectrum::source_table(
                &[(8000.0, 0.0), (18000.0, 0.0)],
                ContentHash([91; 32]),
            )
            .unwrap();
            let geometry = grid(cx, &density, spectrum);
            let absorption = [[0.7], [0.0]];
            let scattering = [[0.5], [scattering_coefficient]];
            let phase = PhaseFunction::Isotropic;
            let field = CellScatteringGrid::try_new(
                cx,
                &geometry,
                [10_000.0],
                &absorption,
                &scattering,
                CellPhaseFunction::Uniform(&phase),
                (8000.0, 18000.0),
                ContentHash([92; 32]),
            )
            .unwrap();
            let op = field
                .estimate_thermal_exchange(cx, photon_config())
                .unwrap();
            assert_eq!(op.material_cells, vec![0, 1]);
            assert_eq!(op.cells.len(), 1);
            assert_eq!(op.cells[0].cell_index, 0);
            let mut illuminated_scatterer = warm_cell(12, 400.0);
            illuminated_scatterer.boundary.source_w = 800.0;
            let original = thermal(
                vec![warm_cell(11, 400.0), illuminated_scatterer],
                &RADIATIVE_SENSIBLE,
                vec![],
            );
            assert!(matches!(
                VolumeRadiativeExchange::try_new(
                    op.clone(),
                    [3000.0],
                    blackbody(400.0),
                    &[binding(0, 11)],
                    &original,
                    budget(),
                    tolerance(),
                    || Ok(())
                ),
                Err(RadiativeExchangeError::Heating(
                    OpticalHeatingError::UnmappedCell(1)
                ))
            ));
            let bindings = [
                binding(0, 11),
                RadiativeThermalCell {
                    minimum_temperature_k: 395.0,
                    maximum_temperature_k: 405.0,
                    ..binding(1, 12)
                },
            ];
            let exchange = observer(op, 3000.0, blackbody(400.0), &bindings, &original);
            assert!(exchange.sources_w().is_empty());
            assert!(exchange.transfers().is_empty());
            assert!(matches!(
                original.advance_observed(&exchange, 0.1, config(1.0), || Ok(())),
                Err(RadiativeExchangeError::TemperatureDomain {
                    thermal_cell: CellId(12),
                    ..
                })
            ));
            assert_eq!(original.temperature_k(CellId(12)).unwrap(), 400.0);
            assert_eq!(original.time_s(), 0.0);
            assert_eq!(exchange.energy().thermal_external_j, 0.0);
        }
    });
}
