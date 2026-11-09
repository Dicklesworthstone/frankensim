//! Finite-emitter pumping, triplet shelving and observed afterglow.
//!
//! Run: cargo run -p fs-render --example saturable_afterglow \
//!   --no-default-features --features optical-heating
//!
//! All coefficients are SYNTHETIC, not measurements of a particular phosphor.
//! A prescribed local photon flux pumps 1e13 emitters in a 1 cm cube for 1 us.
//! The maximum ground-state optical depth is N*sigma/area = 1e-3; pump-field
//! attenuation is neglected explicitly. Only actual absorbed photons supply
//! energy. The 100 ns singlet feeds a 20 us triplet; neither lifetime sets the
//! thermal timestep. Emitted light escapes this transparent observation medium.

use fs_blake3::ContentHash;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_geom::{Point3, Vec3};
use fs_render::volumes::VolumeGrid;
use fs_render::volumes::heating::fluorescence::FluorescenceSource;
use fs_render::volumes::heating::luminescence::{
    LuminescenceSource, SingletTripletLimits, SingletTripletSource,
    saturation::{PrescribedPhotonPump, SaturableEmitterCell, SaturableLuminescentHeating},
};
use fs_render::volumes::transport::luminescence::LuminescenceGrid;
use fs_render::volumes::transport::{AbsorbingGrid, AbsorbingGridBudget, MassAbsorptionSpectrum};
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionError,
    DepositionSimulation, Event, MaterialId, StepConfig, ThermalLaw, TimedEvent,
};
use fs_thermal::luminescence::ArrheniusQuenching;
use std::error::Error;

const EMITTERS: f64 = 1e13;

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

fn source(
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
        MassAbsorptionSpectrum::constant(0.0)?,
        None,
        AbsorbingGridBudget {
            max_grid_cells: 1,
            max_steps_per_segment: 8,
        },
    )?;
    let caloric = Caloric;
    let config = StepConfig {
        max_dt_s: 5e-6,
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
    let singlet = source(450.0, 1e-7, 0.6, [1.0, 0.0], 240)?;
    let triplet = SingletTripletSource::try_new(
        source(550.0, 2e-5, 0.8, [0.0, 1.0], 241)?,
        1e7,
        ArrheniusQuenching::try_new(0.0, 0.0)?,
        SingletTripletLimits::default(),
    )?;
    let pump = PrescribedPhotonPump::try_new(
        ContentHash([242; 32]),
        [400.0; 2],
        [1e-20, 0.0],
        [5e26, 0.0],
    )?;
    let mut heating = SaturableLuminescentHeating::try_new(
        &medium,
        singlet,
        triplet,
        &[SaturableEmitterCell {
            grid_cell: 0,
            thermal_cell: CellId(1),
            initial_populations: [EMITTERS, 0.0, 0.0],
            pump,
        }],
        &thermal,
        1,
        1e-12,
        None,
        || Ok(()),
    )?
    .with_pump_window(0.0, 1e-6)?;
    println!("Synthetic optically thin prescribed-flux pulse; pump switches off at 1 us.");
    println!(
        "time_s,ground_fraction,singlet_fraction,triplet_fraction,absorbed_pump_J,singlet_W,triplet_W,stored_J,native_heat_J,temperature_K,singlet_radiance_W_m2_sr,triplet_radiance_W_m2_sr,closure_J"
    );
    for target in [1e-7, 1e-6, 2e-6, 1e-5, 5e-5, 1e-4] {
        let previous_heat = heating.thermal_input_j();
        let (next, accepted, audit) =
            thermal.advance_observed(&heating, target, config, || Ok(()))?;
        thermal = next;
        heating = accepted;
        if (audit.external_heat_j - (heating.thermal_input_j() - previous_heat)).abs() > 1e-12 {
            return Err("native heat audit disagrees with luminescence heat".into());
        }
        let counts = heating
            .populations(0)
            .ok_or("missing accepted emitter cell")?;
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
            "{target:.6e},{:.6e},{:.6e},{:.6e},{:.9e},{:.9e},{:.9e},{:.9e},{:.9e},{:.12},{:.9e},{:.9e},{:.3e}",
            counts[0] / EMITTERS,
            counts[1] / EMITTERS,
            counts[2] / EMITTERS,
            heating.absorbed_pump_energy_j().iter().sum::<f64>(),
            emission[0],
            emission[1],
            heating.energy().stored_excitation_j,
            heating.thermal_input_j(),
            thermal.temperature_k(CellId(1))?,
            observed.emitted_radiance_w_m2_sr[0],
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
                kernel_id: 243,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}
