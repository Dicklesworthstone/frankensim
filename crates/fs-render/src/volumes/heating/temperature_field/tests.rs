use super::*;
use crate::volumes::VolumeGrid;
use crate::volumes::heating::scattering::ScatteringHeatingTolerance;
use crate::volumes::heating::scattering_feedback::{
    IncidentScatteringBeam, ScatteringFeedbackBudget, TemperatureScatteringHeating,
};
use crate::volumes::transport::scattering::deposition::ScatteringDepositionConfig;
use crate::volumes::transport::{AbsorbingGridBudget, MassAbsorptionSpectrum};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_geom::{Point3, Vec3};
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, DepositedCell, DepositionBudget, Event, MaterialId, StepConfig,
    ThermalLaw, TimedEvent,
};
use std::sync::{
    Arc,
    atomic::{AtomicUsize, Ordering},
};

struct Caloric {
    calls: AtomicUsize,
}
impl ThermalLaw for Caloric {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        self.calls.fetch_add(1, Ordering::Relaxed);
        if h.is_finite() && (1.0..=1000.0).contains(&h) {
            Ok(h)
        } else {
            Err(DepositionError::MaterialDomain)
        }
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0
    }
}
fn caloric() -> Caloric {
    Caloric {
        calls: AtomicUsize::new(0),
    }
}
fn thermal<'a>(
    law: &'a dyn ThermalLaw,
    temperatures: &[f64],
    source_w: f64,
) -> DepositionSimulation<'a> {
    let events = temperatures
        .iter()
        .enumerate()
        .map(|(index, &temperature)| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(index as u64 + 1),
                segment_id: index as u64 + 1,
                body_id: 1,
                material_id: MaterialId(7),
                mass_kg: 1.0,
                specific_h_j_kg: temperature,
                geometry: CellGeometry {
                    start_m: [0.0; 3],
                    end_m: [1.0, 0.0, 0.0],
                    volume_m3: 1.0,
                    surface_area_m2: 6.0,
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w,
                },
            }),
        })
        .collect();
    DepositionSimulation::new(
        events,
        &[(MaterialId(7), law)],
        DepositionBudget {
            max_cells: 8,
            max_links: 8,
            max_events: 16,
            max_substeps: 1000,
        },
    )
    .unwrap()
    .advance_to(0.0, step(1.0), || Ok(()))
    .unwrap()
    .0
}
fn step(dt: f64) -> StepConfig {
    StepConfig {
        max_dt_s: dt,
        energy_tolerance_j: 1e-10,
    }
}
fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        f(&Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 37,
                kernel_id: 1,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}
fn geometry<'g>(cx: &Cx<'_>, density: &'g [f64]) -> AbsorbingGrid<'g> {
    AbsorbingGrid::try_new(
        cx,
        VolumeGrid::new([density.len(), 1, 1], density, [0.0; 3], [1.0; 3]),
        MassAbsorptionSpectrum::constant(0.0).unwrap(),
        None,
        AbsorbingGridBudget {
            max_grid_cells: 8,
            max_steps_per_segment: 32,
        },
    )
    .unwrap()
}
fn absorption(id: u8, lo: f64, hi: f64, coefficients: [f64; 2]) -> TemperatureAbsorptionTable {
    TemperatureAbsorptionTable::try_new(
        ContentHash([id; 32]),
        &[
            (
                lo,
                MassAbsorptionSpectrum::source_table(
                    &[(8000.0, coefficients[0]), (20000.0, coefficients[0])],
                    ContentHash([id; 32]),
                )
                .unwrap(),
            ),
            (
                hi,
                MassAbsorptionSpectrum::source_table(
                    &[(8000.0, coefficients[1]), (20000.0, coefficients[1])],
                    ContentHash([id; 32]),
                )
                .unwrap(),
            ),
        ],
        2,
    )
    .unwrap()
}
fn scattering(lo: f64, hi: f64, coefficients: [f64; 2]) -> TemperatureScatteringTable {
    TemperatureScatteringTable::try_new(
        ContentHash([8; 32]),
        &[
            (
                lo,
                MassScatteringSpectrum::source_table(
                    &[(8000.0, coefficients[0]), (20000.0, coefficients[0])],
                    ContentHash([8; 32]),
                )
                .unwrap(),
            ),
            (
                hi,
                MassScatteringSpectrum::source_table(
                    &[(8000.0, coefficients[1]), (20000.0, coefficients[1])],
                    ContentHash([8; 32]),
                )
                .unwrap(),
            ),
        ],
        2,
    )
    .unwrap()
}
fn binding(cell: usize, receiver: u64, table: usize) -> ThermalScatteringCell {
    ThermalScatteringCell {
        grid_cell: cell,
        thermal_cell: CellId(receiver),
        absorption_table: table,
        scattering_table: 0,
        phase: 0,
    }
}
fn budget() -> TemperatureOpticalFieldBudget {
    TemperatureOpticalFieldBudget {
        max_grid_cells: 8,
        max_tables: 8,
        max_total_temperature_knots: 16,
    }
}

#[test]
fn shared_tables_read_each_native_temperature_once_and_preserve_scalar_interpolation() {
    with_cx(|cx| {
        let law = caloric();
        let state = thermal(&law, &[350.0, 450.0], 0.0);
        let density = [1.0; 4];
        let field = TemperatureOpticalField::try_new(
            geometry(cx, &density),
            [10000.0, 15000.0],
            &[
                absorption(3, 300.0, 500.0, [1.0, 5.0]),
                absorption(4, 300.0, 500.0, [2.0, 6.0]),
            ],
            &[scattering(300.0, 500.0, [0.0, 2.0])],
            &[PhaseFunction::Isotropic],
            &[
                binding(0, 1, 0),
                binding(1, 1, 0),
                binding(2, 1, 1),
                binding(3, 2, 0),
            ],
            &state,
            budget(),
            || Ok(()),
        )
        .unwrap();
        assert_eq!(field.unique_temperature_count(), 2);
        assert_eq!(field.unique_material_count(), 3);
        law.calls.store(0, Ordering::Relaxed);
        let packet = field.initial_state(&state, || Ok(())).unwrap();
        assert_eq!(law.calls.load(Ordering::Relaxed), 2);
        // Independent two-knot arithmetic at 1/4 and 3/4 of the interval.
        assert_eq!(
            packet.absorption_m2_kg(),
            &[[2.0; 2], [2.0; 2], [3.0; 2], [4.0; 2]]
        );
        assert_eq!(
            packet.scattering_m2_kg(),
            &[[0.5; 2], [0.5; 2], [0.5; 2], [1.5; 2]]
        );
        assert_eq!(field.material_bindings().len(), 4);
        assert!(field.has_scattering());
    });
}

#[test]
fn coefficient_change_detection_reuses_flat_laws_and_rejects_out_of_domain_state() {
    with_cx(|cx| {
        let law = caloric();
        let original = thermal(&law, &[350.0], 10.0);
        let warmed = original.advance_to(1.0, step(1.0), || Ok(())).unwrap().0;
        let density = [1.0];
        let flat = TemperatureOpticalField::try_new(
            geometry(cx, &density),
            [10000.0],
            &[absorption(3, 300.0, 500.0, [1.0, 1.0])],
            &[scattering(300.0, 500.0, [0.0, 0.0])],
            &[PhaseFunction::Isotropic],
            &[binding(0, 1, 0)],
            &original,
            budget(),
            || Ok(()),
        )
        .unwrap();
        let mut flat_packet = flat.initial_state(&original, || Ok(())).unwrap();
        assert!(
            !flat
                .prepare(&warmed, &mut flat_packet, &mut || Ok(()))
                .unwrap()
        );
        assert_eq!(flat_packet.absorption_m2_kg(), &[[1.0]]);
        let varying = TemperatureOpticalField::try_new(
            geometry(cx, &density),
            [10000.0],
            &[absorption(3, 300.0, 500.0, [1.0, 3.0])],
            &[scattering(300.0, 500.0, [0.0, 2.0])],
            &[PhaseFunction::Isotropic],
            &[binding(0, 1, 0)],
            &original,
            budget(),
            || Ok(()),
        )
        .unwrap();
        let mut varying_packet = varying.initial_state(&original, || Ok(())).unwrap();
        assert!(
            varying
                .prepare(&warmed, &mut varying_packet, &mut || Ok(()))
                .unwrap()
        );
        assert!((varying_packet.absorption_m2_kg()[0][0] - 1.6).abs() < 2e-15);
        assert!((varying_packet.scattering_m2_kg()[0][0] - 0.6).abs() < 2e-15);
        assert!(
            !varying
                .prepare(&warmed, &mut varying_packet, &mut || Ok(()))
                .unwrap()
        );
        assert!(matches!(
            flat.field(&varying_packet, &mut || Ok(())),
            Err(VolumeTransportError::SourceIdentity)
        ));
        let outside = original.advance_to(16.0, step(1.0), || Ok(())).unwrap().0;
        assert!(matches!(
            flat.prepare(&outside, &mut flat_packet, &mut || Ok(())),
            Err(OpticalFeedbackError::TemperatureDomain { .. })
        ));
        assert_eq!(flat_packet.absorption_m2_kg(), &[[1.0]]);
    });
}

#[test]
fn nonemitting_material_domains_and_native_error_identity_remain_authoritative() {
    with_cx(|cx| {
        let law = caloric();
        let state = thermal(&law, &[400.0], 0.0);
        let density = [1.0, 1.0, 0.0];
        let absorption_tables = [absorption(4, 300.0, 500.0, [0.0, 0.0])];
        let scattering_tables = [scattering(350.0, 450.0, [0.0, 2.0])];
        let field = TemperatureOpticalField::try_new(
            geometry(cx, &density),
            [10000.0],
            &absorption_tables,
            &scattering_tables,
            &[PhaseFunction::Isotropic],
            &[binding(0, 1, 0), binding(1, 1, 0)],
            &state,
            budget(),
            || Ok(()),
        )
        .unwrap();
        assert_eq!(field.material_bindings().len(), 2);
        for material in field.material_bindings() {
            assert_eq!(
                (
                    material.minimum_temperature_k,
                    material.maximum_temperature_k
                ),
                (350.0, 450.0)
            );
        }
        let missing = TemperatureOpticalField::try_new(
            geometry(cx, &density),
            [10000.0],
            &absorption_tables,
            &scattering_tables,
            &[PhaseFunction::Isotropic],
            &[binding(0, 1, 0)],
            &state,
            budget(),
            || Ok(()),
        );
        assert!(matches!(
            missing,
            Err(OpticalFeedbackError::Heating(
                OpticalHeatingError::UnmappedCell(1)
            ))
        ));
        let error = field
            .initial_state(&state, || Err(DepositionError::Budget("native checkpoint")))
            .unwrap_err();
        assert_eq!(
            error,
            OpticalFeedbackError::Thermal(DepositionError::Budget("native checkpoint"))
        );
        let other_absorption = [absorption(5, 401.0, 500.0, [1.0, 1.0])];
        assert!(matches!(
            TemperatureOpticalField::try_new(
                geometry(cx, &density),
                [10000.0],
                &other_absorption,
                &scattering_tables,
                &[PhaseFunction::Isotropic],
                &[binding(0, 1, 0), binding(1, 1, 0)],
                &state,
                budget(),
                || Ok(())
            ),
            Err(OpticalFeedbackError::TemperatureDomain { .. })
        ));
    });
}

#[test]
fn incident_observer_reuses_the_exact_compiled_field_and_original_photon_streams() {
    with_cx(|cx| {
        let law = caloric();
        let state = thermal(&law, &[350.0], 0.0);
        let density = [1.0];
        let beam = IncidentScatteringBeam {
            origin: Point3::new(-1.0, 0.5, 0.5),
            unit_direction: Vec3::new(1.0, 0.0, 0.0),
            wavelengths_nm: [10000.0],
            incident_power_w: [2.0],
        };
        let budget = ScatteringFeedbackBudget {
            max_grid_cells: 8,
            max_tables: 8,
            max_total_temperature_knots: 16,
            max_source_evaluations: 16,
            max_traced_segments: 100_000,
        };
        let config = ScatteringDepositionConfig {
            histories_per_lane: 64,
            max_collisions_per_history: 64,
            max_total_segments: 100_000,
            max_deposition_cells: 8,
            seed: 137,
            source_id: 48,
        };
        let tolerance = ScatteringHeatingTolerance {
            power_roundoff_w: 1e-10,
            max_unresolved_w: 1.0,
            max_total_standard_error_w: 1.0,
            max_cell_standard_error_w: 1.0,
        };
        let original = TemperatureScatteringHeating::try_new(
            geometry(cx, &density),
            beam,
            &[absorption(4, 300.0, 500.0, [1.0, 2.0])],
            &[scattering(300.0, 500.0, [1.0, 1.0])],
            &[PhaseFunction::Isotropic],
            &[binding(0, 1, 0)],
            &state,
            budget,
            config,
            tolerance,
            || Ok(()),
        )
        .unwrap();
        let shared = TemperatureScatteringHeating::from_field(
            original.optical_field().clone(),
            beam,
            &state,
            budget,
            config,
            tolerance,
            || Ok(()),
        )
        .unwrap();
        assert!(Arc::ptr_eq(
            original.optical_field(),
            shared.optical_field()
        ));
        let (first, a, _) = state
            .advance_observed(&original, 0.1, step(0.1), || Ok(()))
            .unwrap();
        let (second, b, _) = state
            .advance_observed(&shared, 0.1, step(0.1), || Ok(()))
            .unwrap();
        assert_eq!(first.cells(), second.cells());
        assert_eq!(a.deposition(), b.deposition());
        assert_eq!(a.energy(), b.energy());
        let mut wrong = beam;
        wrong.wavelengths_nm = [12000.0];
        assert!(matches!(
            TemperatureScatteringHeating::from_field(
                original.optical_field().clone(),
                wrong,
                &state,
                budget,
                config,
                tolerance,
                || Ok(())
            ),
            Err(OpticalFeedbackError::Invalid(_))
        ));
    });
}
