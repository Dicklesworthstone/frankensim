//! Finite-emitter pumping, triplet shelving and observed afterglow.
//!
//! Run: cargo run -p fs-render --example saturable_afterglow \
//!   --no-default-features --features optical-heating
//! Add `-- --depleting` to transport and deplete the incident pump itself.
//!
//! All coefficients are SYNTHETIC, not measurements of a particular phosphor.
//! A prescribed local photon flux pumps 1e13 emitters in a 1 cm cube for 1 us.
//! The maximum ground-state optical depth is N*sigma/area = 1e-3; pump-field
//! attenuation is neglected explicitly. Only actual absorbed photons supply
//! energy. The 100 ns singlet feeds a 20 us triplet; neither lifetime sets the
//! thermal timestep. Emitted light escapes this transparent observation medium.
//!
//! Depleting mode instead puts 1e17 emitters in that cube (initial optical
//! depth 10), illuminates its full 1 cm² face, and pumps for 10 us. Incident
//! flux is 5e26 photons/(m² s), or 5e22 photons/s at 400 nm. The original
//! ground inventory sets absorption as the sample bleaches and recovers.
//! Populations are homogeneous within this single cell; a resolved bleaching
//! front needs spatial refinement. Passive host absorption is explicitly zero.

use fs_blake3::ContentHash;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_geom::{Point3, Vec3};
use fs_render::volumes::VolumeGrid;
use fs_render::volumes::heating::feedback::IncidentBeam;
use fs_render::volumes::heating::fluorescence::FluorescenceSource;
use fs_render::volumes::heating::luminescence::{
    LuminescenceSource, SingletTripletLimits, SingletTripletSource,
    saturation::{
        PrescribedPhotonPump, SaturableEmitterCell, SaturableLuminescentHeating,
        beam::{BeamEmitterCrossSections, PumpColumnLimits},
    },
};
use fs_render::volumes::transport::luminescence::LuminescenceGrid;
use fs_render::volumes::transport::{AbsorbingGrid, AbsorbingGridBudget, MassAbsorptionSpectrum};
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionError,
    DepositionSimulation, Event, MaterialId, StepConfig, ThermalLaw, TimedEvent,
};
use fs_thermal::luminescence::ArrheniusQuenching;
use std::error::Error;

const HC_J_NM: f64 = 6.626_070_15e-34 * 299_792_458.0 * 1e9;
const CROSS_SECTION_M2: f64 = 1e-20;

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

fn run(cx: &Cx<'_>, depleting: bool) -> Result<(), Box<dyn Error>> {
    let emitters = if depleting { 1e17 } else { 1e13 };
    let pulse_end_s = if depleting { 1e-5 } else { 1e-6 };
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
        [CROSS_SECTION_M2, 0.0],
        if depleting { [0.0; 2] } else { [5e26, 0.0] },
    )?;
    let mut heating = SaturableLuminescentHeating::try_new(
        &medium,
        singlet,
        triplet,
        &[SaturableEmitterCell {
            grid_cell: 0,
            thermal_cell: CellId(1),
            initial_populations: [emitters, 0.0, 0.0],
            pump,
        }],
        &thermal,
        1,
        1e-12,
        None,
        || Ok(()),
    )?;
    if depleting {
        heating = heating.with_depleting_beam(
            cx,
            &medium,
            ContentHash([244; 32]),
            IncidentBeam {
                origin: Point3::new(-0.001, 0.005, 0.005),
                unit_direction: Vec3::new(1.0, 0.0, 0.0),
                max_distance_m: 0.012,
                wavelengths_nm: [400.0; 2],
                // Authored full-face area 1e-4 m² times photon flux 5e26.
                incident_power_w: [5e22 * HC_J_NM / 400.0, 0.0],
            },
            &[BeamEmitterCrossSections::authored(
                0,
                ContentHash([245; 32]),
                [[CROSS_SECTION_M2, 0.0], [0.0; 2]],
            )?],
            &thermal,
            1,
            PumpColumnLimits {
                relative_tolerance: 1e-4,
                absolute_energy_tolerance_j: 1e-12,
                maximum_substeps: 4096,
                ..PumpColumnLimits::default()
            },
            || Ok(()),
        )?;
        println!(
            "Synthetic depleting beam: 1e17 emitters, sigma=1e-20 m², initial optical depth=10."
        );
        println!(
            "Full-face 400 nm pump: 5e22 photons/s for 10 us; uniform-cell closure, zero passive host absorption."
        );
    } else {
        println!("Synthetic optically thin prescribed-flux pulse; pump switches off at 1 us.");
    }
    heating = heating.with_pump_window(0.0, pulse_end_s)?;
    println!(
        "time_s,ground_fraction,singlet_fraction,triplet_fraction,absorbed_pump_J,singlet_W,triplet_W,stored_J,native_heat_J,temperature_K,singlet_radiance_W_m2_sr,triplet_radiance_W_m2_sr,closure_J{}",
        if depleting {
            ",incident_J,transmitted_J,molecular_absorption_J,passive_absorption_J,pump_closure_J,cumulative_transmission_fraction,ground_optical_depth,kinetic_substeps"
        } else {
            ""
        }
    );
    let report_times = if depleting {
        &[1e-7, 1e-6, 2e-6, 5e-6, 1e-5, 2e-5, 5e-5, 1e-4][..]
    } else {
        &[1e-7, 1e-6, 2e-6, 1e-5, 5e-5, 1e-4][..]
    };
    for &target in report_times {
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
        let beam_accounts = if let Some(pump) = heating.depleting_pump_energy() {
            let incident = pump.incident_j.iter().sum::<f64>();
            let transmitted = pump.transmitted_j.iter().sum::<f64>();
            let molecular = pump.absorbed_j.iter().sum::<f64>();
            let passive = pump.passive_absorbed_j.iter().sum::<f64>();
            let residual = pump.residual_j.iter().sum::<f64>();
            format!(
                ",{incident:.9e},{transmitted:.9e},{molecular:.9e},{passive:.9e},{residual:.3e},{:.6e},{:.6e},{}",
                if incident > 0.0 {
                    transmitted / incident
                } else {
                    0.0
                },
                counts[0] * CROSS_SECTION_M2 * 0.01 / 1e-6,
                pump.accepted_kinetic_substeps,
            )
        } else {
            String::new()
        };
        println!(
            "{target:.6e},{:.6e},{:.6e},{:.6e},{:.9e},{:.9e},{:.9e},{:.9e},{:.9e},{:.12},{:.9e},{:.9e},{:.3e}{beam_accounts}",
            counts[0] / emitters,
            counts[1] / emitters,
            counts[2] / emitters,
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
    let arguments: Vec<_> = std::env::args().skip(1).collect();
    let depleting = match arguments.as_slice() {
        [] => false,
        [argument] if argument == "--depleting" => true,
        _ => return Err("usage: saturable_afterglow [--depleting]".into()),
    };
    let gate = CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        run(
            &Cx::new(
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
            ),
            depleting,
        )
    })
}
