//! One material field drives incident absorption, reabsorption and cooling.
//!
//! Run: cargo run -p fs-render --example temperature_radiative_balance \
//!   --no-default-features --features optical-heating
//!
//! The IR tables, heat capacity and contact conductance below are explicit
//! synthetic demonstration inputs, not measured material presets. Three
//! wavelength bins cover 8–16 micrometres only; this is NOT bolometric cooling.
//! Both radiative observers see the same native temperature and same compiled
//! coefficient laws at every substep. CSV heat accounts are cumulative joules;
//! source/work counts refer to the most recent two-second interval. This
//! prescribed-step demonstration does not certify final temperature accuracy.

use fs_blake3::ContentHash;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_geom::{Point3, Vec3};
use fs_render::volumes::VolumeGrid;
use fs_render::volumes::heating::exchange::feedback::{
    RadiativeFeedbackBudget, TemperatureRadiativeExchange,
};
use fs_render::volumes::heating::exchange::{RadiativeBath, RadiativeExchangeTolerance};
use fs_render::volumes::heating::feedback::TemperatureAbsorptionTable;
use fs_render::volumes::heating::scattering::ScatteringHeatingTolerance;
use fs_render::volumes::heating::scattering_feedback::{
    IncidentScatteringBeam, ScatteringFeedbackBudget, TemperatureScatteringHeating,
};
use fs_render::volumes::heating::temperature_field::{
    TemperatureOpticalField, TemperatureOpticalFieldBudget, TemperatureScatteringTable,
    ThermalScatteringCell,
};
use fs_render::volumes::transport::scattering::deposition::ScatteringDepositionConfig;
use fs_render::volumes::transport::scattering::exchange::RadiativeExchangeConfig;
use fs_render::volumes::transport::scattering::{MassScatteringSpectrum, PhaseFunction};
use fs_render::volumes::transport::{AbsorbingGrid, AbsorbingGridBudget, MassAbsorptionSpectrum};
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionError,
    DepositionSimulation, Event, MaterialId, StepConfig, ThermalLaw, ThermalLink, TimedEvent,
};
use fs_thermal::observers::CoupledObservers;
use std::error::Error;
use std::sync::Arc;

const WAVELENGTHS_NM: [f64; 3] = [9000.0, 11_000.0, 14_000.0];
const BIN_WIDTHS_NM: [f64; 3] = [2000.0, 2000.0, 4000.0];
const HEAT_CAPACITY_J_KG_K: f64 = 500.0;

struct DemonstrationCaloricLaw;
impl ThermalLaw for DemonstrationCaloricLaw {
    fn temperature_k(&self, specific_h_j_kg: f64) -> Result<f64, DepositionError> {
        let t = specific_h_j_kg / HEAT_CAPACITY_J_KG_K;
        if t.is_finite() && (250.0..=800.0).contains(&t) {
            Ok(t)
        } else {
            Err(DepositionError::MaterialDomain)
        }
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0 / HEAT_CAPACITY_J_KG_K
    }
}
static CALORIC: DemonstrationCaloricLaw = DemonstrationCaloricLaw;

fn cell(index: u64, temperature_k: f64) -> DepositedCell {
    DepositedCell {
        id: CellId(index),
        segment_id: index,
        body_id: 1,
        material_id: MaterialId(1),
        mass_kg: 1.0,
        specific_h_j_kg: HEAT_CAPACITY_J_KG_K * temperature_k,
        geometry: CellGeometry {
            start_m: [(index - 1) as f64 * 0.1, 0.05, 0.05],
            end_m: [index as f64 * 0.1, 0.05, 0.05],
            volume_m3: 0.001,
            surface_area_m2: 0.06,
        },
        boundary: Boundary {
            ambient_k: 300.0,
            htc_w_m2_k: 0.0,
            reservoir: None,
            source_w: 0.0,
        },
    }
}

fn run(cx: &Cx<'_>) -> Result<(), Box<dyn Error>> {
    let schedule = vec![
        TimedEvent {
            time_s: 0.0,
            event: Event::Birth(cell(1, 500.0)),
        },
        TimedEvent {
            time_s: 0.0,
            event: Event::Birth(cell(2, 350.0)),
        },
        TimedEvent {
            time_s: 0.0,
            event: Event::Contact(ThermalLink {
                id: 1,
                a: CellId(1),
                b: CellId(2),
                area_m2: 0.01,
                conductance_w_k: 0.5,
            }),
        },
    ];
    let step = StepConfig {
        max_dt_s: 0.25,
        energy_tolerance_j: 1e-6,
    };
    let mut thermal = DepositionSimulation::new(
        schedule,
        &[(MaterialId(1), &CALORIC)],
        DepositionBudget {
            max_cells: 2,
            max_links: 1,
            max_events: 3,
            max_substeps: 1000,
        },
    )?
    .advance_to(0.0, step, || Ok(()))?
    .0;
    let initial_energy_j = thermal.energy_j();
    let absorption = TemperatureAbsorptionTable::try_new(
        ContentHash([61; 32]),
        &[
            (
                250.0,
                MassAbsorptionSpectrum::source_table(
                    &[(8000.0, 0.004), (16_000.0, 0.003)],
                    ContentHash([62; 32]),
                )?,
            ),
            (
                500.0,
                MassAbsorptionSpectrum::source_table(
                    &[(8000.0, 0.008), (16_000.0, 0.006)],
                    ContentHash([63; 32]),
                )?,
            ),
            (
                800.0,
                MassAbsorptionSpectrum::source_table(
                    &[(8000.0, 0.002), (16_000.0, 0.004)],
                    ContentHash([64; 32]),
                )?,
            ),
        ],
        3,
    )?;
    let scattering = TemperatureScatteringTable::try_new(
        ContentHash([65; 32]),
        &[
            (
                250.0,
                MassScatteringSpectrum::source_table(
                    &[(8000.0, 0.002), (16_000.0, 0.001)],
                    ContentHash([66; 32]),
                )?,
            ),
            (
                800.0,
                MassScatteringSpectrum::source_table(
                    &[(8000.0, 0.004), (16_000.0, 0.002)],
                    ContentHash([67; 32]),
                )?,
            ),
        ],
        2,
    )?;
    let density = [1000.0; 2];
    let geometry = AbsorbingGrid::try_new(
        cx,
        VolumeGrid::new([2, 1, 1], &density, [0.0; 3], [0.1; 3]),
        MassAbsorptionSpectrum::source_table(
            &[(8000.0, 0.004), (16_000.0, 0.003)],
            ContentHash([62; 32]),
        )?,
        None,
        AbsorbingGridBudget {
            max_grid_cells: 2,
            max_steps_per_segment: 16,
        },
    )?;
    let bindings = [0, 1].map(|i| ThermalScatteringCell {
        grid_cell: i,
        thermal_cell: CellId(i as u64 + 1),
        absorption_table: 0,
        scattering_table: 0,
        phase: 0,
    });
    let field = Arc::new(TemperatureOpticalField::try_new(
        geometry,
        WAVELENGTHS_NM,
        &[absorption],
        &[scattering],
        &[PhaseFunction::henyey_greenstein(0.4)?],
        &bindings,
        &thermal,
        TemperatureOpticalFieldBudget {
            max_grid_cells: 2,
            max_tables: 2,
            max_total_temperature_knots: 5,
        },
        || Ok(()),
    )?);
    let incident = TemperatureScatteringHeating::from_field(
        Arc::clone(&field),
        IncidentScatteringBeam {
            origin: Point3::new(-0.1, 0.05, 0.05),
            unit_direction: Vec3::new(1.0, 0.0, 0.0),
            wavelengths_nm: WAVELENGTHS_NM,
            incident_power_w: [500.0; 3],
        },
        &thermal,
        ScatteringFeedbackBudget {
            max_grid_cells: 2,
            max_tables: 2,
            max_total_temperature_knots: 5,
            max_source_evaluations: 1000,
            max_traced_segments: 8_000_000,
        },
        ScatteringDepositionConfig {
            histories_per_lane: 512,
            max_collisions_per_history: 128,
            max_total_segments: 200_000,
            max_deposition_cells: 2,
            seed: 73,
            source_id: 1,
        },
        ScatteringHeatingTolerance {
            power_roundoff_w: 1e-8,
            max_unresolved_w: 1e-6,
            max_total_standard_error_w: 100.0,
            max_cell_standard_error_w: 100.0,
        },
        || Ok(()),
    )?;
    let radiative = TemperatureRadiativeExchange::try_new(
        field,
        BIN_WIDTHS_NM,
        RadiativeBath::Blackbody {
            temperature_k: 300.0,
        },
        &thermal,
        RadiativeFeedbackBudget {
            max_source_evaluations: 1000,
            max_traced_segments: 8_000_000,
            max_thermal_pairs: 1,
        },
        RadiativeExchangeConfig {
            histories_per_emitter_lane: 512,
            max_collisions_per_history: 128,
            max_total_segments: 400_000,
            max_material_cells: 2,
            max_emitter_cells: 2,
            max_pairs: 1,
            seed: 79,
            source_id: 2,
        },
        RadiativeExchangeTolerance {
            coefficient_roundoff_m2_sr: 1e-12,
            power_roundoff_w: 1e-8,
            max_cell_standard_error_w: 100.0,
            max_total_standard_error_w: 200.0,
            max_unresolved_exchange_w: 1e-6,
            max_unresolved_standard_error_w: 1e-6,
        },
        || Ok(()),
    )?;
    assert!(Arc::ptr_eq(incident.optical_field(), radiative.field()));
    let mut coupled = CoupledObservers::new(incident, radiative);
    eprintln!(
        "Synthetic 8–16 micrometre material example; shared coefficients, native conduction and enthalpy; no final-state uncertainty certificate."
    );
    println!(
        "time_s,left_temperature_k,right_temperature_k,incident_absorbed_j,bath_to_material_j,internal_radiation_j,energy_residual_j,source_evaluations,operator_builds"
    );
    for output in 1..=10 {
        let (next, observers, audit) =
            thermal.advance_observed(&coupled, 2.0 * f64::from(output), step, || Ok(()))?;
        thermal = next;
        coupled = observers;
        let incident_j = coupled.first.energy().thermal_input_j;
        let bath_j = coupled.second.energy().thermal_external_j;
        let internal_j: f64 = coupled
            .second
            .energy()
            .internal_transfer_magnitude_j
            .iter()
            .sum();
        let residual = thermal.energy_j() - initial_energy_j - incident_j - bath_j;
        if residual.abs() > 1e-6 || audit.residual_j.abs() > 1e-6 {
            return Err("combined native energy account exceeded its declared tolerance".into());
        }
        println!(
            "{:.2},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9e},{},{}",
            thermal.time_s(),
            thermal.temperature_k(CellId(1))?,
            thermal.temperature_k(CellId(2))?,
            incident_j,
            bath_j,
            internal_j,
            residual,
            coupled.second.source_evaluations(),
            coupled.second.operator_builds()
        );
    }
    Ok(())
}

fn main() -> Result<(), Box<dyn Error>> {
    let gate = CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        run(&Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 73,
                kernel_id: 1,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}
