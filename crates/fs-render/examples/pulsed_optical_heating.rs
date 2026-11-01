//! Finite spectral pulses melt a specimen; native cooling then refreezes it.
//!
//! Run: cargo run -p fs-render --example pulsed_optical_heating \
//!   --no-default-features --features optical-heating,material-state
//!
//! A synthetic 1 g, 1 cm cube receives a 320 J blue box pulse and a 200 J red
//! triangular pulse. Their different absorption coefficients evolve through a
//! native 350 K fusion plateau. Convection increases when the source switches
//! off at 8.5 s; infrared exchange and the thermal history remain continuous.
//!
//! Geometry and density stay fixed. Effective optical coefficients are supplied
//! versus the same material chart's enthalpy, not inferred from a phase mixing
//! rule. These synthetic data and the finite 6–18 micrometre radiation packet
//! demonstrate coupling, not measured-material or bolometric accuracy. Native
//! time and spatial/photon resolutions remain explicit modeling choices.

use fs_blake3::ContentHash;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_geom::{Point3, Vec3};
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve};
use fs_render::volumes::VolumeGrid;
use fs_render::volumes::heating::exchange::feedback::{
    RadiativeFeedbackBudget, TemperatureRadiativeExchange,
};
use fs_render::volumes::heating::exchange::{RadiativeBath, RadiativeExchangeTolerance};
use fs_render::volumes::heating::illumination::schedule::{
    SpectralPowerSchedule, SpectralPowerSegment,
};
use fs_render::volumes::heating::illumination::{
    IlluminationBudget, ScheduledIlluminationSource, TemperatureIlluminationHeating,
};
use fs_render::volumes::heating::joint::JointOpticalHeating;
use fs_render::volumes::heating::scattering::ScatteringHeatingTolerance;
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

const MATERIAL: MaterialId = MaterialId(1);
const CELL: CellId = CellId(1);
const SWITCH_OFF_S: f64 = 8.5;

#[allow(clippy::too_many_arguments)]
fn field<'g, const N: usize>(
    cx: &Cx<'_>,
    density: &'g [f64],
    thermal: &DepositionSimulation<'_>,
    caloric: &EquilibriumEnthalpyPhaseCurve,
    wavelengths: [f64; N],
    support: (f64, f64),
    coefficients: [[f64; N]; 4],
    identity: u8,
) -> Result<Arc<TemperatureOpticalField<'g, N>>, Box<dyn Error>> {
    let rows = caloric
        .knots()
        .iter()
        .zip(coefficients)
        .map(|(knot, kappa)| {
            let mut samples = vec![(support.0, kappa[0])];
            samples.extend(wavelengths.into_iter().zip(kappa));
            samples.push((support.1, kappa[N - 1]));
            Ok((
                knot.specific_enthalpy_j_kg,
                MassAbsorptionSpectrum::source_table(&samples, ContentHash([identity; 32]))?,
            ))
        })
        .collect::<Result<Vec<_>, Box<dyn Error>>>()?;
    let zero = MassScatteringSpectrum::source_table(
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
            &[(0.0, zero), (400_000.0, zero)],
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
    let caloric = EquilibriumEnthalpyPhaseCurve::try_new(
        ContentHash([241; 32]),
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
    let step = StepConfig {
        max_dt_s: 0.05,
        energy_tolerance_j: 1e-8,
    };
    let mut thermal = DepositionSimulation::new(
        vec![
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
                    boundary: boundary(10.0),
                }),
            },
            TimedEvent {
                time_s: SWITCH_OFF_S,
                event: Event::Boundary(CELL, boundary(2000.0)),
            },
        ],
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
        [450.0, 650.0],
        (400.0, 700.0),
        [[0.04, 0.03], [0.04, 0.03], [0.12, 0.10], [0.12, 0.10]],
        211,
    )?;
    let infrared = field(
        cx,
        &density,
        &thermal,
        &caloric,
        [8000.0, 12_000.0, 16_000.0],
        (6000.0, 18_000.0),
        [[0.05; 3], [0.05; 3], [0.15; 3], [0.15; 3]],
        221,
    )?;
    let power = SpectralPowerSchedule::try_new(
        &[
            SpectralPowerSegment {
                start_s: 0.5,
                end_s: 4.5,
                start_power_w: [80.0, 0.0],
                end_power_w: [80.0, 0.0],
            },
            SpectralPowerSegment {
                start_s: 4.5,
                end_s: 6.5,
                start_power_w: [0.0; 2],
                end_power_w: [0.0, 100.0],
            },
            SpectralPowerSegment {
                start_s: 6.5,
                end_s: SWITCH_OFF_S,
                start_power_w: [0.0, 100.0],
                end_power_w: [0.0; 2],
            },
        ],
        3,
    )?;
    let illumination = TemperatureIlluminationHeating::from_scheduled_field(
        visible,
        &[ScheduledIlluminationSource {
            source_id: 1,
            origin: Point3::new(-0.01, 0.005, 0.005),
            unit_direction: Vec3::new(1.0, 0.0, 0.0),
            power,
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
        3,
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
    let mut convection_j = 0.0;
    let mut melting = false;
    let mut melted = false;
    let mut freezing = false;
    eprintln!(
        "Synthetic 1 g cube; 320 J blue box plus 200 J red triangle; source off and convection quench at 8.5 s."
    );
    println!(
        "time_s,temperature_k,liquid_fraction,specific_h_j_kg,blue_incident_j,red_incident_j,absorbed_j,escaped_j,bath_to_material_j,convection_to_material_j,energy_residual_j"
    );
    for output in 1..=40 {
        let target_s = 0.5 * f64::from(output);
        let prior_optical_j = joint.illumination().energy().thermal_input_j
            + joint.radiation().energy().thermal_external_j;
        let (next, accepted, audit) =
            thermal.advance_observed(&joint, target_s, step, || Ok(()))?;
        thermal = next;
        joint = accepted;
        let energy = joint.illumination().energy();
        let bath_j = joint.radiation().energy().thermal_external_j;
        convection_j += audit.external_heat_j - (energy.thermal_input_j + bath_j - prior_optical_j);
        let residual_j =
            thermal.energy_j() - initial_energy_j - energy.thermal_input_j - bath_j - convection_j;
        if residual_j.abs() > 1e-8 || audit.residual_j.abs() > 1e-8 {
            return Err("pulsed thermal balance exceeded its declared energy tolerance".into());
        }
        let h = thermal.cells()[&CELL].specific_h_j_kg;
        let phase = caloric.state_at_specific_enthalpy(h)?;
        let fraction = phase.liquid_mass_fraction();
        melted |= fraction == 1.0;
        if fraction > 0.0 && fraction < 1.0 {
            melting |= target_s < SWITCH_OFF_S;
            freezing |= target_s > SWITCH_OFF_S;
        }
        if target_s > SWITCH_OFF_S && joint.illumination().traced_segments() != 0 {
            return Err("a wholly dark interval traced incident photons".into());
        }
        let escaped_j: f64 = energy
            .ballistic_escape_j
            .iter()
            .chain(&energy.diffuse_escape_j)
            .sum();
        println!(
            "{:.2},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9e}",
            target_s,
            phase.temperature_k(),
            fraction,
            h,
            energy.incident_j[0],
            energy.incident_j[1],
            energy.thermal_input_j,
            escaped_j,
            bath_j,
            convection_j,
            residual_j
        );
    }
    let final_phase = caloric.state_at_specific_enthalpy(thermal.cells()[&CELL].specific_h_j_kg)?;
    let incident = joint.illumination().energy().incident_j;
    if (incident[0] - 320.0).abs() > 1e-9 || (incident[1] - 200.0).abs() > 1e-9 {
        return Err("spectral pulse fluence disagreed with the authored box/triangle areas".into());
    }
    if !melting || !melted || !freezing || final_phase.liquid_mass_fraction() != 0.0 {
        return Err(
            "pulsed example did not complete melting and subsequent dark refreezing".into(),
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
