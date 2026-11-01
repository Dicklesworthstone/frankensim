//! Absorption evolves through melting and freezing on the original phase chart.
//!
//! Run: cargo run -p fs-render --example phase_change_optical_heating \
//!   --no-default-features --features optical-heating,material-state
//!
//! A synthetic 1 g specimen receives 40 W of visible light. Its independently
//! supplied visible and infrared coefficients change with native enthalpy while
//! the original material chart holds temperature at 350 K during melting. A
//! native convection event at 12 s cools and refreezes it under continued light;
//! neither the thermal state nor the optical energy histories are reset.
//!
//! This is a uniform-cell, fixed-density integration example. Optical tables
//! explicitly prescribe an effective coefficient versus enthalpy; their linear
//! interpolation is not a universal solid/liquid mixing rule. Data are synthetic,
//! infrared quadrature covers only 6–18 micrometres, and the prescribed time,
//! spatial and photon resolutions are not a material or accuracy certification.

use fs_blake3::ContentHash;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_geom::{Point3, Vec3};
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve};
use fs_render::volumes::VolumeGrid;
use fs_render::volumes::heating::exchange::feedback::{
    RadiativeFeedbackBudget, TemperatureRadiativeExchange,
};
use fs_render::volumes::heating::exchange::{RadiativeBath, RadiativeExchangeTolerance};
use fs_render::volumes::heating::illumination::{
    IlluminationBudget, IncidentIlluminationSource, TemperatureIlluminationHeating,
};
use fs_render::volumes::heating::joint::JointOpticalHeating;
use fs_render::volumes::heating::scattering::ScatteringHeatingTolerance;
use fs_render::volumes::heating::scattering_feedback::IncidentScatteringBeam;
use fs_render::volumes::heating::temperature_field::{
    EnthalpyAbsorptionTable, EnthalpyOpticalFieldBudget, EnthalpyScatteringCell,
    EnthalpyScatteringTable, TemperatureOpticalField,
};
use fs_render::volumes::transport::scattering::deposition::ScatteringDepositionConfig;
use fs_render::volumes::transport::scattering::exchange::RadiativeExchangeConfig;
use fs_render::volumes::transport::scattering::{MassScatteringSpectrum, PhaseFunction};
use fs_render::volumes::transport::{AbsorbingGrid, AbsorbingGridBudget, MassAbsorptionSpectrum};
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionSimulation, Event,
    MaterialId, StepConfig, TimedEvent,
};
use std::error::Error;
use std::sync::Arc;

const VISIBLE_NM: [f64; 2] = [450.0, 650.0];
const INFRARED_NM: [f64; 3] = [8000.0, 12_000.0, 16_000.0];
const MATERIAL: MaterialId = MaterialId(1);
const CELL: CellId = CellId(1);
const QUENCH_TIME_S: f64 = 12.0;

#[allow(clippy::too_many_arguments)]
fn field<'g, const N: usize>(
    cx: &Cx<'_>,
    density: &'g [f64],
    thermal: &DepositionSimulation<'_>,
    caloric: &EquilibriumEnthalpyPhaseCurve,
    wavelengths: [f64; N],
    support: (f64, f64),
    coefficients: [f64; 4],
    identity: u8,
) -> Result<Arc<TemperatureOpticalField<'g, N>>, Box<dyn Error>> {
    let rows: Vec<_> = caloric
        .knots()
        .iter()
        .zip(coefficients)
        .map(|(knot, kappa)| {
            Ok((
                knot.specific_enthalpy_j_kg,
                MassAbsorptionSpectrum::source_table(
                    &[(support.0, kappa), (support.1, kappa)],
                    ContentHash([identity; 32]),
                )?,
            ))
        })
        .collect::<Result<_, Box<dyn Error>>>()?;
    let zero_scattering = MassScatteringSpectrum::source_table(
        &[(support.0, 0.0), (support.1, 0.0)],
        ContentHash([identity + 1; 32]),
    )?;
    Ok(Arc::new(TemperatureOpticalField::try_new_enthalpy(
        AbsorbingGrid::try_new(
            cx,
            VolumeGrid::new([1; 3], density, [0.0; 3], [0.01; 3]),
            MassAbsorptionSpectrum::constant(0.0)?,
            None,
            AbsorbingGridBudget {
                max_grid_cells: 1,
                max_steps_per_segment: 8,
            },
        )?,
        wavelengths,
        &[EnthalpyAbsorptionTable::try_new(
            ContentHash([identity + 2; 32]),
            caloric.identity(),
            (300.0, 600.0),
            &rows,
            4,
        )?],
        &[EnthalpyScatteringTable::try_new(
            ContentHash([identity + 3; 32]),
            caloric.identity(),
            (300.0, 600.0),
            &[(0.0, zero_scattering), (400_000.0, zero_scattering)],
            2,
        )?],
        &[PhaseFunction::Isotropic],
        &[EnthalpyScatteringCell {
            grid_cell: 0,
            thermal_cell: CELL,
            material_id: MATERIAL,
            absorption_table: 0,
            scattering_table: 0,
            phase: 0,
        }],
        thermal,
        EnthalpyOpticalFieldBudget {
            max_grid_cells: 1,
            max_tables: 2,
            max_total_enthalpy_knots: 6,
        },
        || Ok(()),
    )?))
}

fn run(cx: &Cx<'_>) -> Result<(), Box<dyn Error>> {
    // The SAME immutable chart owns the native thermal law and observed phase.
    let caloric = EquilibriumEnthalpyPhaseCurve::try_new(
        ContentHash([201; 32]),
        [
            (0.0, 300.0, 0.0),
            (50_000.0, 350.0, 0.0),
            (150_000.0, 350.0, 1.0),
            (400_000.0, 600.0, 1.0),
        ]
        .map(
            |(specific_enthalpy_j_kg, temperature_k, liquid_mass_fraction)| EnthalpyPhaseKnot {
                specific_enthalpy_j_kg,
                temperature_k,
                liquid_mass_fraction,
                bulk_density_kg_m3: 1000.0,
            },
        )
        .to_vec(),
    )?;
    let boundary = |htc_w_m2_k| Boundary {
        ambient_k: 300.0,
        htc_w_m2_k,
        reservoir: None,
        source_w: 0.0,
    };
    let schedule = vec![
        TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CELL,
                segment_id: 1,
                body_id: 1,
                material_id: MATERIAL,
                mass_kg: 0.001,
                specific_h_j_kg: 0.0,
                geometry: CellGeometry {
                    start_m: [0.0, 0.005, 0.005],
                    end_m: [0.01, 0.005, 0.005],
                    volume_m3: 1e-6,
                    surface_area_m2: 6e-4,
                },
                boundary: boundary(0.0),
            }),
        },
        TimedEvent {
            time_s: QUENCH_TIME_S,
            event: Event::Boundary(CELL, boundary(2000.0)),
        },
    ];
    let step = StepConfig {
        max_dt_s: 0.05,
        energy_tolerance_j: 1e-8,
    };
    let mut thermal = DepositionSimulation::new(
        schedule,
        &[(MATERIAL, &caloric)],
        DepositionBudget {
            max_cells: 1,
            max_links: 0,
            max_events: 2,
            max_substeps: 1000,
        },
    )?
    .advance_to(0.0, step, || Ok(()))?
    .0;
    let initial_energy_j = thermal.energy_j();
    let density = [1000.0];
    let visible = field(
        cx,
        &density,
        &thermal,
        &caloric,
        VISIBLE_NM,
        (400.0, 700.0),
        [0.04, 0.04, 0.12, 0.12],
        211,
    )?;
    let infrared = field(
        cx,
        &density,
        &thermal,
        &caloric,
        INFRARED_NM,
        (6000.0, 18_000.0),
        [0.05, 0.05, 0.15, 0.15],
        221,
    )?;
    let illumination = TemperatureIlluminationHeating::from_field(
        visible,
        &[IncidentIlluminationSource {
            source_id: 1,
            beam: IncidentScatteringBeam {
                origin: Point3::new(-0.01, 0.005, 0.005),
                unit_direction: Vec3::new(1.0, 0.0, 0.0),
                wavelengths_nm: VISIBLE_NM,
                incident_power_w: [20.0; 2],
            },
        }],
        &thermal,
        IlluminationBudget {
            max_sources: 1,
            max_source_evaluations: 1000,
            max_traced_segments: 1_000_000,
            max_retained_deposition_cells: 1,
        },
        ScatteringDepositionConfig {
            histories_per_lane: 64,
            max_collisions_per_history: 32,
            max_total_segments: 100_000,
            max_deposition_cells: 1,
            seed: 431,
            source_id: 1,
        },
        ScatteringHeatingTolerance {
            power_roundoff_w: 1e-10,
            max_unresolved_w: 1e-8,
            max_total_standard_error_w: 0.0,
            max_cell_standard_error_w: 0.0,
        },
        || Ok(()),
    )?;
    let radiation = TemperatureRadiativeExchange::try_new(
        infrared,
        [4000.0; 3],
        RadiativeBath::Blackbody {
            temperature_k: 300.0,
        },
        &thermal,
        RadiativeFeedbackBudget {
            max_source_evaluations: 1000,
            max_traced_segments: 1_000_000,
            max_thermal_pairs: 0,
        },
        RadiativeExchangeConfig {
            histories_per_emitter_lane: 128,
            max_collisions_per_history: 32,
            max_total_segments: 100_000,
            max_material_cells: 1,
            max_emitter_cells: 1,
            max_pairs: 0,
            seed: 432,
            source_id: 2,
        },
        RadiativeExchangeTolerance {
            coefficient_roundoff_m2_sr: 1e-14,
            power_roundoff_w: 1e-10,
            max_cell_standard_error_w: 0.1,
            max_total_standard_error_w: 0.1,
            max_unresolved_exchange_w: 1e-8,
            max_unresolved_standard_error_w: 1e-8,
        },
        || Ok(()),
    )?;
    let mut joint = JointOpticalHeating::try_new(illumination, radiation, &thermal, || Ok(()))?;
    let mut boundary_heat_j = 0.0;
    let mut melted = false;
    let mut heating_plateau = false;
    let mut freezing_plateau = false;
    eprintln!(
        "Synthetic 1 g specimen; 350 K melting plateau; 40 W visible illumination; convection quench at 12 s; finite 6–18 micrometre radiation."
    );
    println!(
        "time_s,temperature_k,liquid_mass_fraction,specific_h_j_kg,visible_kappa_m2_kg,absorbed_j,bath_to_material_j,convection_to_material_j,energy_residual_j"
    );
    for output in 1..=30 {
        let prior_optical_j = joint.illumination().energy().thermal_input_j
            + joint.radiation().energy().thermal_external_j;
        let (next, accepted, audit) =
            thermal.advance_observed(&joint, f64::from(output), step, || Ok(()))?;
        thermal = next;
        joint = accepted;
        let absorbed_j = joint.illumination().energy().thermal_input_j;
        let bath_j = joint.radiation().energy().thermal_external_j;
        // Native external heat includes all three owners. Separate the original
        // optical accounts to retain the convection term over every interval.
        boundary_heat_j += audit.external_heat_j - (absorbed_j + bath_j - prior_optical_j);
        let residual_j =
            thermal.energy_j() - initial_energy_j - absorbed_j - bath_j - boundary_heat_j;
        if residual_j.abs() > 1e-8 || audit.residual_j.abs() > 1e-8 {
            return Err("phase-change energy balance exceeded its declared tolerance".into());
        }
        let h = thermal.cells()[&CELL].specific_h_j_kg;
        let phase = caloric.state_at_specific_enthalpy(h)?;
        let fraction = phase.liquid_mass_fraction();
        melted |= fraction == 1.0;
        if fraction > 0.0 && fraction < 1.0 {
            heating_plateau |= thermal.time_s() < QUENCH_TIME_S;
            freezing_plateau |= thermal.time_s() > QUENCH_TIME_S;
        }
        let packet = joint
            .illumination()
            .optical_field()
            .initial_state(&thermal, || Ok(()))?;
        println!(
            "{:.2},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9e}",
            thermal.time_s(),
            phase.temperature_k(),
            fraction,
            h,
            packet.absorption_m2_kg()[0][0],
            absorbed_j,
            bath_j,
            boundary_heat_j,
            residual_j
        );
    }
    let final_phase = caloric.state_at_specific_enthalpy(thermal.cells()[&CELL].specific_h_j_kg)?;
    if !heating_plateau || !melted || !freezing_plateau || final_phase.liquid_mass_fraction() != 0.0
    {
        return Err("example did not complete both the melting and freezing trajectory".into());
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
                seed: 431,
                kernel_id: 1,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}
