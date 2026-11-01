//! Distributed visible illumination and infrared cooling of one thermal body.
//!
//! Run: cargo run -p fs-render --example distributed_optical_heating \
//!   --no-default-features --features optical-heating
//!
//! Six prescribed weighted rays illuminate four physical cells. Visible source
//! powers and three infrared quadrature bins have distinct wavelength packets;
//! both optical fields read the same native enthalpy state. Conduction and net
//! radiation use the original native owners. Source tables are synthetic model
//! inputs. This finite-band example is an integration/conservation demonstration,
//! without measured-material, source-quadrature or final-accuracy certification.

use fs_blake3::ContentHash;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_geom::{Point3, Vec3};
use fs_render::volumes::VolumeGrid;
use fs_render::volumes::heating::exchange::feedback::{
    RadiativeFeedbackBudget, TemperatureRadiativeExchange,
};
use fs_render::volumes::heating::exchange::{RadiativeBath, RadiativeExchangeTolerance};
use fs_render::volumes::heating::feedback::TemperatureAbsorptionTable;
use fs_render::volumes::heating::illumination::{
    IlluminationBudget, IncidentIlluminationSource, TemperatureIlluminationHeating,
};
use fs_render::volumes::heating::joint::JointOpticalHeating;
use fs_render::volumes::heating::scattering::ScatteringHeatingTolerance;
use fs_render::volumes::heating::scattering_feedback::IncidentScatteringBeam;
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
use std::error::Error;
use std::sync::Arc;

const VISIBLE_NM: [f64; 2] = [450.0, 650.0];
const INFRARED_NM: [f64; 3] = [9000.0, 11_000.0, 14_000.0];
const INFRARED_WIDTHS_NM: [f64; 3] = [2000.0, 2000.0, 4000.0];
const CP: f64 = 500.0;

struct Caloric;
impl ThermalLaw for Caloric {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        let t = h / CP;
        if t.is_finite() && (250.0..=900.0).contains(&t) {
            Ok(t)
        } else {
            Err(DepositionError::MaterialDomain)
        }
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0 / CP
    }
}
static CALORIC: Caloric = Caloric;

#[allow(clippy::too_many_arguments)]
fn field<'g, const N: usize>(
    cx: &Cx<'_>,
    density: &'g [f64],
    thermal: &DepositionSimulation<'_>,
    wavelengths: [f64; N],
    support: (f64, f64),
    absorption: [f64; 3],
    scattering: f64,
    identity: u8,
) -> Result<Arc<TemperatureOpticalField<'g, N>>, Box<dyn Error>> {
    let absorption_rows: Vec<_> = [250.0, 500.0, 900.0]
        .into_iter()
        .zip(absorption)
        .map(|(t, k)| {
            Ok((
                t,
                MassAbsorptionSpectrum::source_table(
                    &[(support.0, k), (support.1, k)],
                    ContentHash([identity; 32]),
                )?,
            ))
        })
        .collect::<Result<_, Box<dyn Error>>>()?;
    let scattering_spectrum = MassScatteringSpectrum::source_table(
        &[(support.0, scattering), (support.1, scattering)],
        ContentHash([identity + 1; 32]),
    )?;
    let bindings: Vec<_> = (0..density.len())
        .map(|i| ThermalScatteringCell {
            grid_cell: i,
            thermal_cell: CellId(i as u64 + 1),
            absorption_table: 0,
            scattering_table: 0,
            phase: 0,
        })
        .collect();
    Ok(Arc::new(TemperatureOpticalField::try_new(
        AbsorbingGrid::try_new(
            cx,
            VolumeGrid::new([2, 2, 1], density, [0.0; 3], [0.1; 3]),
            MassAbsorptionSpectrum::constant(0.0)?,
            None,
            AbsorbingGridBudget {
                max_grid_cells: 4,
                max_steps_per_segment: 32,
            },
        )?,
        wavelengths,
        &[TemperatureAbsorptionTable::try_new(
            ContentHash([identity + 2; 32]),
            &absorption_rows,
            3,
        )?],
        &[TemperatureScatteringTable::try_new(
            ContentHash([identity + 3; 32]),
            &[(250.0, scattering_spectrum), (900.0, scattering_spectrum)],
            2,
        )?],
        &[PhaseFunction::henyey_greenstein(0.35)?],
        &bindings,
        thermal,
        TemperatureOpticalFieldBudget {
            max_grid_cells: 4,
            max_tables: 2,
            max_total_temperature_knots: 5,
        },
        || Ok(()),
    )?))
}

fn run(cx: &Cx<'_>) -> Result<(), Box<dyn Error>> {
    let mut schedule = Vec::new();
    for (i, t) in [450.0, 400.0, 350.0, 300.0].into_iter().enumerate() {
        let x = (i % 2) as f64;
        let y = (i / 2) as f64;
        schedule.push(TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(i as u64 + 1),
                segment_id: i as u64 + 1,
                body_id: 1,
                material_id: MaterialId(1),
                mass_kg: 1.0,
                specific_h_j_kg: CP * t,
                geometry: CellGeometry {
                    start_m: [0.1 * x, 0.1 * (y + 0.5), 0.05],
                    end_m: [0.1 * (x + 1.0), 0.1 * (y + 0.5), 0.05],
                    volume_m3: 0.001,
                    surface_area_m2: 0.06,
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w: 0.0,
                },
            }),
        });
    }
    for (id, a, b, conductance) in [
        (1, 1, 2, 0.5),
        (2, 3, 4, 0.5),
        (3, 1, 3, 0.25),
        (4, 2, 4, 0.25),
    ] {
        schedule.push(TimedEvent {
            time_s: 0.0,
            event: Event::Contact(ThermalLink {
                id,
                a: CellId(a),
                b: CellId(b),
                area_m2: 0.01,
                conductance_w_k: conductance,
            }),
        });
    }
    let step = StepConfig {
        max_dt_s: 0.25,
        energy_tolerance_j: 1e-6,
    };
    let mut thermal = DepositionSimulation::new(
        schedule,
        &[(MaterialId(1), &CALORIC)],
        DepositionBudget {
            max_cells: 4,
            max_links: 8,
            max_events: 8,
            max_substeps: 1000,
        },
    )?
    .advance_to(0.0, step, || Ok(()))?
    .0;
    let initial_energy = thermal.energy_j();
    let density = [1000.0; 4];
    let visible = field(
        cx,
        &density,
        &thermal,
        VISIBLE_NM,
        (400.0, 700.0),
        [0.004, 0.008, 0.005],
        0.002,
        101,
    )?;
    let infrared = field(
        cx,
        &density,
        &thermal,
        INFRARED_NM,
        (8000.0, 16_000.0),
        [0.002, 0.006, 0.004],
        0.001,
        111,
    )?;

    let mut sources = Vec::new();
    // A 1000 W prescribed spatial quadrature, with 25%/75% lower/upper power.
    for (i, (y, weight)) in [(0.025, 0.1), (0.075, 0.15), (0.125, 0.3), (0.175, 0.45)]
        .into_iter()
        .enumerate()
    {
        sources.push(IncidentIlluminationSource {
            source_id: i as u64 + 1,
            beam: IncidentScatteringBeam {
                origin: Point3::new(-0.1, y, 0.05),
                unit_direction: Vec3::new(1.0, 0.0, 0.0),
                wavelengths_nm: VISIBLE_NM,
                incident_power_w: [600.0 * weight, 400.0 * weight],
            },
        });
    }
    // Two opposed 150 W fill rays use a different source spectrum.
    for (id, y) in [(11, 0.05), (12, 0.15)] {
        sources.push(IncidentIlluminationSource {
            source_id: id,
            beam: IncidentScatteringBeam {
                origin: Point3::new(0.3, y, 0.05),
                unit_direction: Vec3::new(-1.0, 0.0, 0.0),
                wavelengths_nm: VISIBLE_NM,
                incident_power_w: [45.0, 105.0],
            },
        });
    }
    let illumination = TemperatureIlluminationHeating::from_field(
        visible,
        &sources,
        &thermal,
        IlluminationBudget {
            max_sources: 6,
            max_source_evaluations: 1000,
            max_traced_segments: 10_000_000,
            max_retained_deposition_cells: 24,
        },
        ScatteringDepositionConfig {
            histories_per_lane: 256,
            max_collisions_per_history: 128,
            max_total_segments: 200_000,
            max_deposition_cells: 4,
            seed: 421,
            source_id: 1,
        },
        ScatteringHeatingTolerance {
            power_roundoff_w: 1e-8,
            max_unresolved_w: 1e-6,
            max_total_standard_error_w: 200.0,
            max_cell_standard_error_w: 200.0,
        },
        || Ok(()),
    )?;
    let radiation = TemperatureRadiativeExchange::try_new(
        infrared,
        INFRARED_WIDTHS_NM,
        RadiativeBath::Blackbody {
            temperature_k: 300.0,
        },
        &thermal,
        RadiativeFeedbackBudget {
            max_source_evaluations: 1000,
            max_traced_segments: 10_000_000,
            max_thermal_pairs: 6,
        },
        RadiativeExchangeConfig {
            histories_per_emitter_lane: 256,
            max_collisions_per_history: 128,
            max_total_segments: 500_000,
            max_material_cells: 4,
            max_emitter_cells: 4,
            max_pairs: 6,
            seed: 422,
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
    let mut joint = JointOpticalHeating::try_new(illumination, radiation, &thermal, || Ok(()))?;
    eprintln!(
        "Synthetic four-cell specimen; six visible rays totaling 1300 W; finite 8–16 micrometre cooling; prescribed spatial/spectral quadrature and time steps."
    );
    println!(
        "time_s,t1_k,t2_k,t3_k,t4_k,absorbed_j,bath_to_material_j,internal_radiation_j,energy_residual_j,illumination_builds,radiation_builds"
    );
    for output in 1..=10 {
        let (next, accepted, audit) =
            thermal.advance_observed(&joint, 2.0 * f64::from(output), step, || Ok(()))?;
        thermal = next;
        joint = accepted;
        let absorbed = joint.illumination().energy().thermal_input_j;
        let bath = joint.radiation().energy().thermal_external_j;
        let internal: f64 = joint
            .radiation()
            .energy()
            .internal_transfer_magnitude_j
            .iter()
            .sum();
        let residual = thermal.energy_j() - initial_energy - absorbed - bath;
        if residual.abs() > 1e-6 || audit.residual_j.abs() > 1e-6 {
            return Err("joint physical energy balance exceeded its declared tolerance".into());
        }
        println!(
            "{:.2},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9e},{},{}",
            thermal.time_s(),
            thermal.temperature_k(CellId(1))?,
            thermal.temperature_k(CellId(2))?,
            thermal.temperature_k(CellId(3))?,
            thermal.temperature_k(CellId(4))?,
            absorbed,
            bath,
            internal,
            residual,
            joint.illumination().transport_builds(),
            joint.radiation().operator_builds()
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
                seed: 421,
                kernel_id: 1,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}
