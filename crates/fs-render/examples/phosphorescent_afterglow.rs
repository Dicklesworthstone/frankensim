//! Pump a singlet, transfer population into a triplet, and observe dark afterglow.
//!
//! Run: cargo run -p fs-render --example phosphorescent_afterglow \
//!   --no-default-features --features optical-heating
//!
//! All coefficients are SYNTHETIC: this demonstrates spin-state kinetics,
//! thermal conservation and volume observation, not a measured phosphor.
//! The dilute model has no traps, saturation or oxygen chemistry. A 0.5 s
//! optical pulse excites a 1 cm cube; the singlet has a 100 ns intrinsic
//! lifetime while the triplet continues to emit for seconds.

use fs_blake3::ContentHash;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_geom::{Point3, Vec3};
use fs_render::volumes::VolumeGrid;
use fs_render::volumes::heating::fluorescence::FluorescenceSource;
use fs_render::volumes::heating::illumination::schedule::{
    SpectralPowerSchedule, SpectralPowerSegment,
};
use fs_render::volumes::heating::luminescence::{
    LuminescenceSource, LuminescentOpticalHeating, SingletTripletLimits, SingletTripletSource,
};
use fs_render::volumes::transport::luminescence::LuminescenceGrid;
use fs_render::volumes::transport::{AbsorbingGrid, AbsorbingGridBudget, MassAbsorptionSpectrum};
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionError,
    DepositionSimulation, Event, MaterialId, StepConfig, ThermalLaw, TimedEvent,
};
use fs_thermal::luminescence::ArrheniusQuenching;
use std::error::Error;

struct Caloric;
impl ThermalLaw for Caloric {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        if !h.is_finite() || !(0.0..=100_000.0).contains(&h) {
            return Err(DepositionError::MaterialDomain);
        }
        Ok(300.0 + h / 1000.0)
    }
    fn maximum_temperature_slope(&self) -> f64 {
        0.001
    }
}

fn state(
    level_nm: f64,
    lifetime_s: f64,
    yield_: f64,
    weights: [f64; 2],
    id: u8,
) -> Result<LuminescenceSource<2>, Box<dyn Error>> {
    Ok(LuminescenceSource::try_new(
        FluorescenceSource::try_new(
            ContentHash([id; 32]),
            (380.0, 420.0),
            (295.0, 400.0),
            yield_,
            [600.0, 750.0],
            weights,
        )?,
        level_nm,
        lifetime_s,
    )?)
}

fn run(cx: &Cx<'_>) -> Result<(), Box<dyn Error>> {
    let density = [1000.0];
    let medium = AbsorbingGrid::try_new(
        cx,
        VolumeGrid::new([1; 3], &density, [0.0; 3], [0.01; 3]),
        MassAbsorptionSpectrum::source_table(
            &[(380.0, 10.0), (420.0, 10.0), (600.0, 0.0), (750.0, 0.0)],
            ContentHash([191; 32]),
        )?,
        None,
        AbsorbingGridBudget {
            max_grid_cells: 1,
            max_steps_per_segment: 8,
        },
    )?;
    let deposited = medium.deposit_beam_power(
        cx,
        Point3::new(-0.01, 0.005, 0.005),
        Vec3::new(1.0, 0.0, 0.0),
        0.03,
        &[400.0],
        &[0.001],
        1,
    )?;
    let caloric = Caloric;
    let config = StepConfig {
        max_dt_s: 0.25,
        energy_tolerance_j: 1e-12,
    };
    let empty = DepositionSimulation::new(
        vec![TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(1),
                segment_id: 1,
                body_id: 1,
                material_id: MaterialId(1),
                mass_kg: 0.001,
                specific_h_j_kg: 0.0,
                geometry: CellGeometry {
                    start_m: [0.0, 0.005, 0.005],
                    end_m: [0.01, 0.005, 0.005],
                    volume_m3: 1e-6,
                    surface_area_m2: 6e-4,
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w: 0.0,
                },
            }),
        }],
        &[(MaterialId(1), &caloric)],
        DepositionBudget {
            max_cells: 1,
            max_links: 0,
            max_events: 1,
            max_substeps: 1000,
        },
    )?;
    let mut thermal = empty.advance_to(0.0, config, || Ok(()))?.0;
    let singlet = state(450.0, 1e-7, 0.6, [1.0, 0.0], 192)?;
    let triplet = state(550.0, 2.0, 0.8, [0.0, 1.0], 193)?;
    let heating = LuminescentOpticalHeating::try_new(
        &medium,
        singlet,
        deposited,
        &[(0, CellId(1))],
        &thermal,
        1,
        1e-12,
        1e-10,
        || Ok(()),
    )?
    .with_singlet_triplet(
        SingletTripletSource::try_new(
            triplet,
            1e7,
            ArrheniusQuenching::try_new(0.0, 0.0)?,
            SingletTripletLimits::default(),
        )?,
        &thermal,
        || Ok(()),
    )?;
    let mut heating = heating.with_pump_schedule(SpectralPowerSchedule::try_new(
        &[SpectralPowerSegment {
            start_s: 0.0,
            end_s: 0.5,
            start_power_w: [0.001],
            end_power_w: [0.001],
        }],
        1,
    )?)?;
    println!("Synthetic dilute singlet/triplet example; pump switches off at 0.5 s.");
    println!(
        "time_s,singlet_W,triplet_W,stored_J,native_heat_J,temperature_K,triplet_radiance_W_m2_sr,closure_J"
    );
    for target in [0.25, 0.5, 1.0, 2.0, 5.0, 10.0] {
        let (next, accepted, _) = thermal.advance_observed(&heating, target, config, || Ok(()))?;
        thermal = next;
        heating = accepted;
        let emission = heating.cell_emission_w()[0].1;
        let view = LuminescenceGrid::try_new(
            cx,
            &medium,
            heating.emission_source_identity(),
            [600.0, 750.0],
            heating.cell_emission_w(),
            1,
        )?;
        let observed = view.transport(
            cx,
            Point3::new(-0.01, 0.005, 0.005),
            Vec3::new(1.0, 0.0, 0.0),
            0.03,
        )?;
        println!(
            "{target:.2},{:.9e},{:.9e},{:.9e},{:.9e},{:.9},{:.9e},{:.3e}",
            emission[0],
            emission[1],
            heating.energy().stored_excitation_j,
            heating.pump_energy().thermal_input_j,
            thermal.temperature_k(CellId(1))?,
            observed.emitted_radiance_w_m2_sr[1],
            heating.energy().residual_j,
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
                seed: 0,
                kernel_id: 194,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}
