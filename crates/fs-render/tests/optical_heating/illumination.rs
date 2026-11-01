use super::*;
use fs_blake3::ContentHash;
use fs_render::volumes::heating::feedback::{OpticalFeedbackError, TemperatureAbsorptionTable};
use fs_render::volumes::heating::illumination::{
    IlluminationBudget, IncidentIlluminationSource, TemperatureIlluminationHeating,
};
use fs_render::volumes::heating::scattering::ScatteringHeatingTolerance;
use fs_render::volumes::heating::scattering_feedback::IncidentScatteringBeam;
use fs_render::volumes::heating::temperature_field::{
    TemperatureOpticalField, TemperatureOpticalFieldBudget, TemperatureScatteringTable,
    ThermalScatteringCell,
};
use fs_render::volumes::transport::scattering::deposition::ScatteringDepositionConfig;
use fs_thermal::enthalpy::DepositionObserver;
use std::sync::{
    Arc,
    atomic::{AtomicUsize, Ordering},
};

fn limits() -> IlluminationBudget {
    IlluminationBudget {
        max_sources: 8,
        max_source_evaluations: 1000,
        max_traced_segments: 1_000_000,
        max_retained_deposition_cells: 64,
    }
}
fn photons() -> ScatteringDepositionConfig {
    ScatteringDepositionConfig {
        histories_per_lane: 64,
        max_collisions_per_history: 64,
        max_total_segments: 1_000_000,
        max_deposition_cells: 8,
        seed: 173,
        source_id: 91,
    }
}
fn tolerance() -> ScatteringHeatingTolerance {
    ScatteringHeatingTolerance {
        power_roundoff_w: 1e-10,
        max_unresolved_w: 100.0,
        max_total_standard_error_w: 100.0,
        max_cell_standard_error_w: 100.0,
    }
}
fn illumination_field<'g>(
    cx: &Cx<'_>,
    density: &'g [f64],
    thermal: &DepositionSimulation<'_>,
    absorption: &[(f64, f64)],
    scattering: f64,
) -> Arc<TemperatureOpticalField<'g, 1>> {
    let rows: Vec<_> = absorption
        .iter()
        .map(|&(t, value)| (t, MassAbsorptionSpectrum::constant(value).unwrap()))
        .collect();
    let a = TemperatureAbsorptionTable::try_new(ContentHash([71; 32]), &rows, 16).unwrap();
    let s = TemperatureScatteringTable::try_new(
        ContentHash([72; 32]),
        &[
            (
                absorption[0].0,
                MassScatteringSpectrum::constant(scattering).unwrap(),
            ),
            (
                absorption.last().unwrap().0,
                MassScatteringSpectrum::constant(scattering).unwrap(),
            ),
        ],
        16,
    )
    .unwrap();
    let bindings: Vec<_> = (0..density.len())
        .map(|i| ThermalScatteringCell {
            grid_cell: i,
            thermal_cell: CellId(11 + i as u64),
            absorption_table: 0,
            scattering_table: 0,
            phase: 0,
        })
        .collect();
    Arc::new(
        TemperatureOpticalField::try_new(
            grid(cx, density, MassAbsorptionSpectrum::constant(0.0).unwrap()),
            [550.0],
            &[a],
            &[s],
            &[PhaseFunction::Isotropic],
            &bindings,
            thermal,
            TemperatureOpticalFieldBudget {
                max_grid_cells: 8,
                max_tables: 2,
                max_total_temperature_knots: 32,
            },
            || Ok(()),
        )
        .unwrap(),
    )
}
fn source(id: u64, power: f64, right: bool, cells: usize) -> IncidentIlluminationSource<1> {
    IncidentIlluminationSource {
        source_id: id,
        beam: IncidentScatteringBeam {
            origin: Point3::new(if right { cells as f64 + 1.0 } else { -1.0 }, 0.5, 0.5),
            unit_direction: Vec3::new(if right { -1.0 } else { 1.0 }, 0.0, 0.0),
            wavelengths_nm: [550.0],
            incident_power_w: [power],
        },
    }
}
fn observer<'g>(
    field: Arc<TemperatureOpticalField<'g, 1>>,
    sources: &[IncidentIlluminationSource<1>],
    thermal: &DepositionSimulation<'_>,
) -> TemperatureIlluminationHeating<'g, 1> {
    TemperatureIlluminationHeating::from_field(
        field,
        sources,
        thermal,
        limits(),
        photons(),
        tolerance(),
        || Ok(()),
    )
    .unwrap()
}
fn close(a: f64, b: f64) {
    assert!(
        (a - b).abs() < 1e-10 * a.abs().max(b.abs()).max(1.0),
        "{a} != {b}"
    );
}

#[test]
fn opposing_sources_match_independent_beer_balance_and_count_time_once() {
    with_cx(|cx| {
        let original = thermal(
            vec![cell(11, 1.0, 1.0, 0.0), cell(12, 2.0, 1.0, 0.0)],
            &SENSIBLE,
            vec![],
        );
        let field = illumination_field(
            cx,
            &[1.0, 2.0],
            &original,
            &[(290.0, 0.4), (600.0, 0.4)],
            0.0,
        );
        let initial = observer(
            field,
            &[source(20, 7.0, true, 2), source(10, 12.0, false, 2)],
            &original,
        );
        assert_eq!(initial.traced_segments(), 0);
        assert!(initial.deposition(10).is_none());
        let (hot, accepted, audit) = original
            .advance_observed(&initial, 0.5, config(0.125), || Ok(()))
            .unwrap();
        let t0 = (-0.4_f64).exp();
        let t1 = (-0.8_f64).exp();
        let q0 = (12.0 + 7.0 * t1) * (1.0 - t0);
        let q1 = (7.0 + 12.0 * t0) * (1.0 - t1);
        close(hot.cells()[&CellId(11)].specific_h_j_kg, 0.5 * q0);
        close(hot.cells()[&CellId(12)].specific_h_j_kg, 0.5 * q1 / 2.0);
        close(accepted.energy().duration_s, 0.5);
        close(accepted.energy().incident_j[0], 19.0 * 0.5);
        close(accepted.energy().absorbed_j[0], (q0 + q1) * 0.5);
        close(
            accepted.energy().ballistic_escape_j[0],
            19.0 * t0 * t1 * 0.5,
        );
        assert_eq!(accepted.energy().diffuse_escape_j, [0.0]);
        assert_eq!(accepted.energy().unresolved_j, [0.0]);
        close(audit.external_heat_j, accepted.energy().thermal_input_j);
        close(audit.stored_energy_change_j, audit.external_heat_j);
        assert_eq!(accepted.source_evaluations(), audit.substeps);
        assert_eq!(accepted.transport_builds(), 1);
        assert_eq!(accepted.traced_segments(), 2);
        assert_eq!(accepted.retained_deposition_cells(), 4);

        // Fixed coefficients reuse all sources even with no new photon budget.
        let cached = accepted.with_work_limits(1000, 0);
        let (_, continued, _) = hot
            .advance_observed(&cached, 1.0, config(0.125), || Ok(()))
            .unwrap();
        assert_eq!(continued.traced_segments(), 0);
        assert_eq!(continued.transport_builds(), 0);
        close(continued.energy().duration_s, 1.0);
        close(
            continued.energy().thermal_input_j,
            2.0 * accepted.energy().thermal_input_j,
        );
        assert!(
            hot.advance_observed(&cached.with_replica_seed(1), 1.0, config(0.125), || Ok(()))
                .is_err()
        );
        assert!(
            hot.advance_observed(
                &cached.with_histories_per_lane(128).unwrap(),
                1.0,
                config(0.125),
                || Ok(())
            )
            .is_err()
        );
    });
}

#[test]
fn canonical_source_ids_preserve_streams_under_reordering_and_inactive_sources() {
    with_cx(|cx| {
        let original = thermal(vec![cell(11, 1.0, 1.0, 0.0)], &SENSIBLE, vec![]);
        let field = illumination_field(cx, &[1.0], &original, &[(290.0, 0.4), (600.0, 0.4)], 1.5);
        let beams = [source(30, 12.0, false, 1), source(10, 7.0, true, 1)];
        let forward = observer(Arc::clone(&field), &beams, &original);
        let reverse = observer(Arc::clone(&field), &[beams[1], beams[0]], &original);
        assert_eq!(forward.model_identity(), reverse.model_identity());
        let (hot_a, a, audit_a) = original
            .advance_observed(&forward, 0.25, config(0.125), || Ok(()))
            .unwrap();
        let (hot_b, b, audit_b) = original
            .advance_observed(&reverse, 0.25, config(0.125), || Ok(()))
            .unwrap();
        assert_eq!(hot_a.cells(), hot_b.cells());
        assert_eq!(a.energy(), b.energy());
        assert_eq!(audit_a, audit_b);
        for id in [10, 30] {
            assert_eq!(a.deposition(id), b.deposition(id));
        }
        let stream_a = a.deposition(10).unwrap().config;
        let stream_b = a.deposition(30).unwrap().config;
        assert_ne!(
            (stream_a.seed, stream_a.source_id),
            (stream_b.seed, stream_b.source_id)
        );
        let extra = observer(
            Arc::clone(&field),
            &[source(5, 0.0, false, 1), beams[0], beams[1]],
            &original,
        );
        let (_, extra, _) = original
            .advance_observed(&extra, 0.25, config(0.125), || Ok(()))
            .unwrap();
        assert_eq!(a.energy(), extra.energy());
        for id in [10, 30] {
            assert_eq!(a.deposition(id), extra.deposition(id));
        }
        assert_eq!(extra.deposition(5).unwrap().traced_segments, 0);
        assert!(matches!(
            TemperatureIlluminationHeating::from_field(
                field,
                &[beams[0], beams[0]],
                &original,
                limits(),
                photons(),
                tolerance(),
                || Ok(())
            ),
            Err(OpticalFeedbackError::Invalid(_))
        ));
    });
}

#[test]
fn temperature_feedback_rebuilds_all_sources_and_converges_to_the_absorption_ode() {
    with_cx(|cx| {
        let original = thermal(vec![cell(11, 1.0, 1.0, 0.0)], &SENSIBLE, vec![]);
        // κ(T)=0.1+0.005*(T-300), T=300+h/4; h'=19*(1-exp(-κ)).
        let field = illumination_field(cx, &[1.0], &original, &[(290.0, 0.05), (600.0, 1.6)], 0.0);
        let illumination = observer(
            field,
            &[source(10, 12.0, false, 1), source(20, 7.0, true, 1)],
            &original,
        );
        let (coarse, _, _) = original
            .advance_observed(&illumination, 0.5, config(0.1), || Ok(()))
            .unwrap();
        let (fine, source, audit) = original
            .advance_observed(&illumination, 0.5, config(0.025), || Ok(()))
            .unwrap();
        let c = 0.005 / 4.0;
        let exact = ((0.1_f64.exp_m1() * (c * 19.0 * 0.5_f64).exp()).ln_1p() - 0.1) / c;
        let coarse_h = coarse.cells()[&CellId(11)].specific_h_j_kg;
        let fine_h = fine.cells()[&CellId(11)].specific_h_j_kg;
        assert!(coarse_h < fine_h && fine_h < exact);
        assert!((exact - fine_h) < 0.3 * (exact - coarse_h));
        close(fine_h, audit.external_heat_j);
        close(audit.external_heat_j, source.energy().thermal_input_j);
        assert_eq!(source.transport_builds(), audit.substeps);
        assert_eq!(source.traced_segments(), 2 * audit.substeps as u64);
    });
}

#[test]
fn whole_source_work_storage_and_native_transactions_are_bounded() {
    with_cx(|cx| {
        let original = thermal(vec![cell(11, 1.0, 1.0, 0.0)], &SENSIBLE, vec![]);
        let field = illumination_field(cx, &[1.0], &original, &[(290.0, 0.4), (600.0, 0.4)], 0.0);
        let beams = [source(10, 12.0, false, 1), source(20, 7.0, true, 1)];
        for budget in [
            IlluminationBudget {
                max_traced_segments: 1,
                ..limits()
            },
            IlluminationBudget {
                max_retained_deposition_cells: 1,
                ..limits()
            },
            IlluminationBudget {
                max_source_evaluations: 1,
                ..limits()
            },
        ] {
            let input = TemperatureIlluminationHeating::from_field(
                Arc::clone(&field),
                &beams,
                &original,
                budget,
                photons(),
                tolerance(),
                || Ok(()),
            )
            .unwrap();
            assert!(
                original
                    .advance_observed(&input, 0.5, config(0.125), || Ok(()))
                    .is_err()
            );
            assert_eq!(input.time_s(), 0.0);
            assert_eq!(input.energy().incident_j, [0.0]);
            assert_eq!(input.source_evaluations(), 0);
            assert_eq!(original.cells()[&CellId(11)].specific_h_j_kg, 0.0);
        }
        assert!(
            TemperatureIlluminationHeating::from_field(
                Arc::clone(&field),
                &beams,
                &original,
                IlluminationBudget {
                    max_sources: 1,
                    ..limits()
                },
                photons(),
                tolerance(),
                || Ok(())
            )
            .is_err()
        );
        let input = observer(Arc::clone(&field), &beams, &original);
        let mut checks = 0;
        assert!(
            original
                .advance_observed(&input, 0.5, config(0.125), || {
                    checks += 1;
                    if checks == 40 {
                        Err(DepositionError::Cancelled)
                    } else {
                        Ok(())
                    }
                })
                .is_err()
        );
        assert_eq!(input.time_s(), 0.0);
        assert_eq!(input.energy().thermal_input_j, 0.0);
        let empty = TemperatureIlluminationHeating::from_field(
            field,
            &[],
            &original,
            IlluminationBudget {
                max_sources: 0,
                max_traced_segments: 0,
                max_retained_deposition_cells: 0,
                ..limits()
            },
            photons(),
            tolerance(),
            || Ok(()),
        )
        .unwrap();
        let (unchanged, empty, audit) = original
            .advance_observed(&empty, 0.5, config(0.125), || Ok(()))
            .unwrap();
        assert_eq!(unchanged.cells(), original.cells());
        assert_eq!(audit.external_heat_j, 0.0);
        assert_eq!(empty.energy().duration_s, 0.5);
        assert_eq!(empty.traced_segments(), 0);
    });
}

#[test]
fn source_error_allowances_apply_to_the_whole_illumination() {
    with_cx(|cx| {
        let original = thermal(vec![cell(11, 1.0, 1.0, 0.0)], &SENSIBLE, vec![]);
        let field = illumination_field(cx, &[1.0], &original, &[(290.0, 0.4), (600.0, 0.4)], 2.0);
        let beams = [source(10, 12.0, false, 1), source(20, 7.0, true, 1)];
        let photon = ScatteringDepositionConfig {
            max_collisions_per_history: 1,
            ..photons()
        };
        let mut observed = TemperatureIlluminationHeating::from_field(
            Arc::clone(&field),
            &beams,
            &original,
            limits(),
            photon,
            tolerance(),
            || Ok(()),
        )
        .unwrap();
        observed.sources(&original, &mut || Ok(())).unwrap();
        let a = observed.deposition(10).unwrap();
        let b = observed.deposition(20).unwrap();
        let error = |r: &fs_render::volumes::transport::scattering::deposition::ScatteringPowerDeposition<1>| r.absorbed.standard_error_w[0].max(r.unresolved.standard_error_w[0]);
        let error_cap = 0.5 * (error(a).max(error(b)) + error(a) + error(b));
        let cell_cap = 0.5
            * (a.absorbed.standard_error_w[0].max(b.absorbed.standard_error_w[0])
                + a.absorbed.standard_error_w[0]
                + b.absorbed.standard_error_w[0]);
        let unresolved_cap = 0.5
            * (a.unresolved.mean_w[0].max(b.unresolved.mean_w[0])
                + a.unresolved.mean_w[0]
                + b.unresolved.mean_w[0]);
        assert!(error(a) > 0.0 && error(b) > 0.0);
        for tolerance in [
            ScatteringHeatingTolerance {
                max_total_standard_error_w: error_cap,
                ..tolerance()
            },
            ScatteringHeatingTolerance {
                max_cell_standard_error_w: cell_cap,
                ..tolerance()
            },
            ScatteringHeatingTolerance {
                max_unresolved_w: unresolved_cap,
                ..tolerance()
            },
        ] {
            // Each physical source fits individually; their combined demand does not.
            for beam in &beams {
                let mut single = TemperatureIlluminationHeating::from_field(
                    Arc::clone(&field),
                    &[*beam],
                    &original,
                    limits(),
                    photon,
                    tolerance,
                    || Ok(()),
                )
                .unwrap();
                single.sources(&original, &mut || Ok(())).unwrap();
            }
            let mut together = TemperatureIlluminationHeating::from_field(
                Arc::clone(&field),
                &beams,
                &original,
                limits(),
                photon,
                tolerance,
                || Ok(()),
            )
            .unwrap();
            assert!(together.sources(&original, &mut || Ok(())).is_err());
            assert!(together.sources_w().is_none());
        }
    });
}

#[test]
fn material_temperature_is_prepared_once_for_all_source_rays() {
    struct Counted {
        reads: AtomicUsize,
    }
    impl ThermalLaw for Counted {
        fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
            self.reads.fetch_add(1, Ordering::Relaxed);
            SENSIBLE.temperature_k(h)
        }
        fn maximum_temperature_slope(&self) -> f64 {
            SENSIBLE.maximum_temperature_slope()
        }
    }
    with_cx(|cx| {
        let law = Counted {
            reads: AtomicUsize::new(0),
        };
        let original = DepositionSimulation::new(
            vec![TimedEvent {
                time_s: 0.0,
                event: Event::Birth(cell(11, 1.0, 1.0, 0.0)),
            }],
            &[(MaterialId(7), &law)],
            DepositionBudget {
                max_cells: 8,
                max_links: 8,
                max_events: 16,
                max_substeps: 1000,
            },
        )
        .unwrap()
        .advance_to(0.0, config(0.1), || Ok(()))
        .unwrap()
        .0;
        let field = illumination_field(cx, &[1.0], &original, &[(290.0, 0.4), (600.0, 0.4)], 0.0);
        let mut illumination = observer(
            field,
            &[
                source(10, 12.0, false, 1),
                source(20, 7.0, true, 1),
                source(30, 9.0, false, 1),
            ],
            &original,
        );
        law.reads.store(0, Ordering::Relaxed);
        illumination.sources(&original, &mut || Ok(())).unwrap();
        assert_eq!(law.reads.load(Ordering::Relaxed), 1);
        assert_eq!(illumination.transport_builds(), 1);
        illumination.sources(&original, &mut || Ok(())).unwrap();
        assert_eq!(law.reads.load(Ordering::Relaxed), 2);
        assert_eq!(illumination.transport_builds(), 1);
        assert_eq!(illumination.traced_segments(), 3);
    });
}
