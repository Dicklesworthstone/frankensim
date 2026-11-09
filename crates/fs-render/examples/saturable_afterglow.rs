//! Finite-emitter pumping, triplet shelving and observed afterglow.
//!
//! Run: cargo run -p fs-render --example saturable_afterglow \
//!   --no-default-features --features optical-heating
//! Add `-- --depleting` to transport and deplete the incident pump itself.
//! Add `-- --counterpropagating` for two differently colored, staggered pumps
//! sharing the same finite inventory in two original thermal/optical cells.
//! Add `--annihilation` to any mode to include explicitly authored triplet-pair
//! quenching and delayed singlet fluorescence; compare the late singlet column.
//! Add `--recycling` to any mode for finite-ground photon reabsorption. Used
//! alone, it retains the existing explicitly thin/maintained local pump and
//! enables photon return without a depleted beam. Recycling explicitly changes
//! the synthetic triplet level to 650 nm:
//! 600 nm singlet light can excite it, while 750 nm triplet light escapes.
//! The endpoint view includes the same accepted molecular absorption.
//! Add `--scattered-recycling` to include original elastic scattering in both
//! photon return and the endpoint view; this implies `--recycling`. Its explicit
//! synthetic spectrum is zero at both pump colors, preserving the direct pump,
//! and 0.05/0.1 m²/kg at 600/750 nm with isotropic scattering. Recycling uses
//! 32 complete walks per cell and wavelength; each endpoint view averages 128
//! native histories with a fixed seed. Refine both sample counts and seeds to
//! assess sampling sensitivity separately from kinetic/thermal refinement.
//! Add `--background-pump` to a depleted-pump mode to also supply an explicitly
//! maintained 400 nm local field: 5e24 photons/(m² s), or a 5e4/s excitation
//! hazard per available ground-state emitter, under the same global pulse.
//! This separate source competes with the beams for the original populations;
//! its photons are not also included in a beam's declared incident power.
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
//!
//! Counterpropagating mode divides the same cube and inventory into two cells.
//! A 400 nm beam enters from the left for 0–8 us, and a 410 nm beam enters
//! from the right for 3–10 us, each supplying 2.5e22 photons/s. Their incoherent
//! intensities compete for each cell's shared ground state during the overlap.
//! The cross sections remain explicitly synthetic and constant at both colors.
//! Output retains separate beam energy accounts and each cell's ground fraction;
//! the reported temperature is the equal-mass mean of the two native cells.
//! The optional pair coefficient is 1e-18 m³/s per event, with a 0.5 singlet
//! yield and each cell's original active volume. Pair events restore ground
//! availability during pumping and continue producing singlets after switch-off.

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
        annihilation::{AnnihilationLimits, PhosphorescenceAnnihilation, TripletAnnihilation},
        beam::{BeamEmitterCrossSections, DepletingPumpBeam, PumpBeamsLimits, PumpColumnLimits},
        recycling::{
            PhotonRecyclingGeometry, RecyclingEmitterCrossSections, SaturableRecyclingLimits,
        },
    },
};
use fs_render::volumes::transport::luminescence::{
    LuminescenceGrid, LuminescenceScatteringBudget, ScatteringLuminescenceGrid,
};
use fs_render::volumes::transport::scattering::{
    MassScatteringSpectrum, PhaseFunction, ScatteringGrid,
};
use fs_render::volumes::transport::{AbsorbingGrid, AbsorbingGridBudget, MassAbsorptionSpectrum};
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionError,
    DepositionSimulation, Event, MaterialId, StepConfig, ThermalLaw, TimedEvent,
};
use fs_thermal::luminescence::ArrheniusQuenching;
use std::error::Error;

const HC_J_NM: f64 = 6.626_070_15e-34 * 299_792_458.0 * 1e9;
const CROSS_SECTION_M2: f64 = 1e-20;
const SCATTERED_VIEW_SAMPLES: u32 = 128;
const SCATTERED_VIEW_SEED: u64 = 251;

fn observe_scattered<'g>(
    cx: &Cx<'_>,
    field: &LuminescenceGrid<'_, 'g, 2>,
    medium: &ScatteringGrid<'g>,
) -> Result<[f64; 2], Box<dyn Error>> {
    let view = ScatteringLuminescenceGrid::try_new(cx, field, medium)?;
    let mut radiance = [0.0; 2];
    let mut remaining_visits = 32_768_u32;
    for path in 0..SCATTERED_VIEW_SAMPLES {
        let sample = view.sample_ray(
            cx,
            Point3::new(-0.01, 0.005, 0.005),
            Vec3::new(1.0, 0.0, 0.0),
            LuminescenceScatteringBudget {
                max_scattering_events: 128,
                max_cell_visits: remaining_visits.min(512),
            },
            SCATTERED_VIEW_SEED,
            u64::from(path),
        )?;
        remaining_visits = remaining_visits
            .checked_sub(sample.cell_visits)
            .ok_or("scattered endpoint observation exhausted its aggregate visit budget")?;
        for lane in 0..2 {
            radiance[lane] +=
                (sample.radiance_w_m2_sr[lane] - radiance[lane]) / f64::from(path + 1);
        }
    }
    if radiance.iter().any(|value| !value.is_finite()) {
        return Err("scattered endpoint radiance is not representable".into());
    }
    Ok(radiance)
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum PumpMode {
    Prescribed,
    Depleting,
    Counterpropagating,
}

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
    recycling: bool,
) -> Result<LuminescenceSource<2>, Box<dyn Error>> {
    Ok(LuminescenceSource::try_new(
        FluorescenceSource::try_new(
            ContentHash([id; 32]),
            (
                380.0,
                if recycling {
                    level_nm.min(600.0)
                } else {
                    420.0
                },
            ),
            (295.0, 400.0),
            yield_,
            [600.0, 750.0],
            weights,
        )?,
        level_nm,
        lifetime_s,
    )?)
}

fn run(
    cx: &Cx<'_>,
    mode: PumpMode,
    annihilation: bool,
    recycling: bool,
    scattered_recycling: bool,
    background_pump: bool,
) -> Result<(), Box<dyn Error>> {
    let depleting = mode != PumpMode::Prescribed;
    let counterpropagating = mode == PumpMode::Counterpropagating;
    let emitters = if depleting { 1e17 } else { 1e13 };
    let pulse_end_s = if depleting { 1e-5 } else { 1e-6 };
    let cell_count = if counterpropagating { 2 } else { 1 };
    let cell_width_m = 0.01 / cell_count as f64;
    let density = vec![1000.0; cell_count];
    let medium = AbsorbingGrid::try_new(
        cx,
        VolumeGrid::new(
            [cell_count, 1, 1],
            &density,
            [0.0; 3],
            [cell_width_m, 0.01, 0.01],
        ),
        MassAbsorptionSpectrum::constant(0.0)?,
        None,
        AbsorbingGridBudget {
            max_grid_cells: cell_count,
            max_steps_per_segment: 8,
        },
    )?;
    let (absorption, scattering) = if scattered_recycling {
        (
            None,
            Some(ScatteringGrid::try_new(
                cx,
                medium,
                MassScatteringSpectrum::source_table(
                    &[
                        (380.0, 0.0),
                        (400.0, 0.0),
                        (410.0, 0.0),
                        (600.0, 0.05),
                        (750.0, 0.1),
                        (780.0, 0.1),
                    ],
                    ContentHash([250; 32]),
                )?,
                PhaseFunction::Isotropic,
            )?),
        )
    } else {
        (Some(medium), None)
    };
    let medium = match &scattering {
        Some(medium) => medium.absorption_grid(),
        None => absorption.as_ref().expect("original absorption snapshot"),
    };
    let caloric = Caloric;
    let config = StepConfig {
        max_dt_s: 5e-6,
        energy_tolerance_j: 1e-12,
    };
    let empty = DepositionSimulation::new(
        (0..cell_count)
            .map(|index| TimedEvent {
                time_s: 0.0,
                event: Event::Birth(DepositedCell {
                    id: CellId(index as u64 + 1),
                    segment_id: index as u64 + 1,
                    body_id: 1,
                    material_id: MaterialId(1),
                    mass_kg: 0.001 / cell_count as f64,
                    specific_h_j_kg: 0.0,
                    geometry: CellGeometry {
                        start_m: [index as f64 * cell_width_m, 0.005, 0.005],
                        end_m: [(index + 1) as f64 * cell_width_m, 0.005, 0.005],
                        volume_m3: 1e-6 / cell_count as f64,
                        surface_area_m2: 2.0 * (1e-4 + 0.02 * cell_width_m),
                    },
                    boundary: Boundary {
                        ambient_k: 300.0,
                        htc_w_m2_k: 0.0,
                        reservoir: None,
                        source_w: 0.0,
                    },
                }),
            })
            .collect(),
        &[(MaterialId(1), &caloric)],
        DepositionBudget {
            max_cells: cell_count,
            max_links: 0,
            max_events: cell_count,
            max_substeps: 1000,
        },
    )?;
    let mut thermal = empty.advance_to(0.0, config, || Ok(()))?.0;
    let singlet = source(450.0, 1e-7, 0.6, [1.0, 0.0], 240, recycling)?;
    let triplet = SingletTripletSource::try_new(
        source(
            if recycling { 650.0 } else { 550.0 },
            2e-5,
            0.8,
            [0.0, 1.0],
            241,
            recycling,
        )?,
        1e7,
        ArrheniusQuenching::try_new(0.0, 0.0)?,
        SingletTripletLimits::default(),
    )?;
    let pump = PrescribedPhotonPump::try_new(
        ContentHash([242; 32]),
        [400.0; 2],
        [CROSS_SECTION_M2, 0.0],
        if background_pump {
            [5e24, 0.0]
        } else if depleting {
            [0.0; 2]
        } else {
            [5e26, 0.0]
        },
    )?;
    let bindings: Vec<_> = (0..cell_count)
        .map(|grid_cell| SaturableEmitterCell {
            grid_cell,
            thermal_cell: CellId(grid_cell as u64 + 1),
            initial_populations: [emitters / cell_count as f64, 0.0, 0.0],
            pump,
        })
        .collect();
    let mut heating = SaturableLuminescentHeating::try_new(
        medium,
        singlet,
        triplet,
        &bindings,
        &thermal,
        cell_count,
        1e-12,
        None,
        || Ok(()),
    )?;
    if background_pump {
        println!(
            "Synthetic maintained background: 400 nm, 5e24 photons/(m² s), shared ground-state hazard=5e4/s, gated by the global pump window."
        );
    }
    if annihilation {
        let pairs: Vec<_> = (0..cell_count)
            .map(|grid_cell| {
                Ok(PhosphorescenceAnnihilation {
                    grid_cell,
                    law: TripletAnnihilation::try_new(1e-18, 1e-6 / cell_count as f64, 0.5)?,
                })
            })
            .collect::<Result<_, fs_thermal::luminescence::transfer::SingletTripletError>>()?;
        heating = heating.with_triplet_annihilation(
            &pairs,
            AnnihilationLimits {
                absolute_energy_tolerance_j: 1e-12,
                relative_tolerance: 1e-4,
                maximum_substeps: 4096,
            },
        )?;
        println!(
            "Synthetic triplet pairs: k=1e-18 m³/s per event, productive singlet probability=0.5; two triplets consumed per event."
        );
    }
    if counterpropagating {
        let cross_sections: Vec<_> = (0..cell_count)
            .map(|grid_cell| {
                BeamEmitterCrossSections::authored(
                    grid_cell,
                    ContentHash([245; 32]),
                    [[CROSS_SECTION_M2, 0.0], [0.0; 2]],
                )
            })
            .collect::<Result<_, _>>()?;
        let beams = [
            DepletingPumpBeam {
                source_identity: ContentHash([244; 32]),
                beam: IncidentBeam {
                    origin: Point3::new(-0.001, 0.005, 0.005),
                    unit_direction: Vec3::new(1.0, 0.0, 0.0),
                    max_distance_m: 0.012,
                    wavelengths_nm: [400.0; 2],
                    incident_power_w: [2.5e22 * HC_J_NM / 400.0, 0.0],
                },
                cross_sections: &cross_sections,
                pump_window_s: Some((0.0, 8e-6)),
            },
            DepletingPumpBeam {
                source_identity: ContentHash([246; 32]),
                beam: IncidentBeam {
                    origin: Point3::new(0.011, 0.005, 0.005),
                    unit_direction: Vec3::new(-1.0, 0.0, 0.0),
                    max_distance_m: 0.012,
                    wavelengths_nm: [410.0; 2],
                    incident_power_w: [2.5e22 * HC_J_NM / 410.0, 0.0],
                },
                cross_sections: &cross_sections,
                pump_window_s: Some((3e-6, 1e-5)),
            },
        ];
        heating = heating.with_depleting_beams(
            cx,
            medium,
            &beams,
            &thermal,
            2,
            2 * cell_count,
            PumpBeamsLimits {
                column: PumpColumnLimits {
                    relative_tolerance: 1e-4,
                    absolute_energy_tolerance_j: 1e-12,
                    maximum_substeps: 4096,
                    ..PumpColumnLimits::default()
                },
                maximum_coupling_iterations: 128,
            },
            || Ok(()),
        )?;
        println!(
            "Synthetic counterpropagating pumps: 1e17 shared emitters in two cells, initial total optical depth=10."
        );
        println!(
            "Left 400 nm: 0–8 us; right 410 nm: 3–10 us; each 2.5e22 photons/s, full-face, incoherent, zero passive absorption."
        );
    } else if depleting {
        heating = heating.with_depleting_beam(
            cx,
            medium,
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
    if recycling {
        let sections: Vec<_> = (0..cell_count)
            .map(|grid_cell| {
                RecyclingEmitterCrossSections::authored(
                    grid_cell,
                    ContentHash([247; 32]),
                    [[0.0, 1e-20], [0.0; 2]],
                )
            })
            .collect::<Result<_, _>>()?;
        let geometry = PhotonRecyclingGeometry {
            histories_per_emitter: 32,
            // Two states per wavelength, plus any external pump paths.
            maximum_paths: if scattered_recycling { 512 } else { 256 },
            maximum_total_segments: if scattered_recycling { 16_384 } else { 1024 },
            seed: 248,
            source_id: 249,
            minimum_temperature_k: 295.0,
            maximum_temperature_k: 400.0,
            fraction_roundoff: 1e-12,
        };
        let limits = SaturableRecyclingLimits {
            pump: PumpBeamsLimits {
                column: PumpColumnLimits {
                    relative_tolerance: 1e-4,
                    absolute_energy_tolerance_j: 1e-12,
                    maximum_substeps: 4096,
                    maximum_cell_evaluations: 1_048_576,
                    ..PumpColumnLimits::default()
                },
                maximum_coupling_iterations: 128,
            },
            annihilation: AnnihilationLimits {
                relative_tolerance: 1e-4,
                absolute_energy_tolerance_j: 1e-12,
                maximum_substeps: 4096,
            },
            maximum_feedback_iterations: 128,
            maximum_secant_history: 8,
        };
        heating = if let Some(scattering) = &scattering {
            heating.with_scattered_photon_recycling(
                cx, scattering, &sections, &thermal, geometry, 128, limits, || Ok(()),
            )?
        } else {
            heating.with_photon_recycling(
                cx, medium, &sections, &thermal, geometry, limits, || Ok(()),
            )?
        };
        if scattered_recycling {
            println!(
                "Synthetic elastic recycling: 32 uniform-cell isotropic walks/cell/wavelength, sigma(600nm,G→T)=1e-20 m², triplet level=650nm; no molecular absorption at 750nm."
            );
            println!(
                "Original isotropic scattering: kappa_s(400/410nm)=0, kappa_s(600nm)=0.05 m²/kg, kappa_s(750nm)=0.1 m²/kg. Direct pumps remain unscattered."
            );
            println!(
                "Endpoint radiance uses the same snapshot and {SCATTERED_VIEW_SAMPLES} native histories, seed={SCATTERED_VIEW_SEED}; complete-walk sampling sensitivity is separate from kinetic and thermal refinement."
            );
        } else {
            println!(
                "Synthetic active recycling: 32 uniform-cell isotropic rays/cell, sigma(600nm,G→T)=1e-20 m², triplet level=650nm; 750nm remains transparent."
            );
            println!(
                "Ballistic open-grid quadrature, quasistatic flight; ray-count/seed refinement is separate from kinetic and thermal refinement."
            );
        }
    }
    heating = heating.with_pump_window(0.0, pulse_end_s)?;
    let absorbed_label = if background_pump {
        println!(
            "Total external absorption includes the maintained background and beams; each beam's energy columns contain only that beam."
        );
        "total_external_absorbed_J"
    } else {
        "absorbed_pump_J"
    };
    println!(
        "time_s,ground_fraction,singlet_fraction,triplet_fraction,{absorbed_label},singlet_W,triplet_W,stored_J,native_heat_J,mean_temperature_K,singlet_radiance_W_m2_sr,triplet_radiance_W_m2_sr,closure_J{}{}{}",
        if annihilation { ",pair_heat_J" } else { "" },
        if recycling {
            ",recycled_excitation_J,secondary_heat_J,escaped_emission_J"
        } else {
            ""
        },
        if counterpropagating {
            ",left_400nm_incident_J,left_transmitted_J,left_molecular_absorption_J,left_passive_absorption_J,left_pump_closure_J,right_410nm_incident_J,right_transmitted_J,right_molecular_absorption_J,right_passive_absorption_J,right_pump_closure_J,left_ground_fraction,right_ground_fraction,kinetic_substeps,coupling_iterations"
        } else if depleting {
            ",incident_J,transmitted_J,molecular_absorption_J,passive_absorption_J,pump_closure_J,cumulative_transmission_fraction,ground_optical_depth,kinetic_substeps"
        } else {
            ""
        }
    );
    let report_times = if counterpropagating {
        &[1e-7, 1e-6, 3e-6, 5e-6, 8e-6, 1e-5, 2e-5, 5e-5, 1e-4][..]
    } else if depleting {
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
        let mut counts = [0.0; 3];
        let mut mean_temperature_k = 0.0;
        for cell in 0..cell_count {
            let population = heating
                .populations(cell)
                .ok_or("missing accepted emitter cell")?;
            for level in 0..3 {
                counts[level] += population[level];
            }
            mean_temperature_k +=
                thermal.temperature_k(CellId(cell as u64 + 1))? / cell_count as f64;
        }
        let emission = heating
            .cell_emission_w()
            .iter()
            .fold([0.0; 2], |sum, (_, power)| {
                [sum[0] + power[0], sum[1] + power[1]]
            });
        let view = if recycling {
            heating.recycling_emission_grid(cx, medium, cell_count)?
        } else {
            LuminescenceGrid::try_new(
                cx,
                medium,
                heating.emission_source_identity(),
                [600.0, 750.0],
                heating.cell_emission_w(),
                cell_count,
            )?
        };
        let observed = if let Some(scattering) = &scattering {
            observe_scattered(cx, &view, scattering)?
        } else {
            view.transport(
                cx,
                Point3::new(-0.01, 0.005, 0.005),
                Vec3::new(1.0, 0.0, 0.0),
                0.03,
            )?
            .emitted_radiance_w_m2_sr
        };
        let beam_accounts = if counterpropagating {
            let report = heating
                .depleting_pump_report()
                .ok_or("missing accepted simultaneous-pump accounts")?;
            let mut columns = String::new();
            // Each beam owns its wavelengths: do not combine ordinal lanes
            // from distinct colors into one spectral-energy packet.
            for pump in &report.beams {
                columns.push_str(&format!(
                    ",{:.9e},{:.9e},{:.9e},{:.9e},{:.3e}",
                    pump.incident_j.iter().sum::<f64>(),
                    pump.transmitted_j.iter().sum::<f64>(),
                    pump.absorbed_j.iter().sum::<f64>(),
                    pump.passive_absorbed_j.iter().sum::<f64>(),
                    pump.residual_j.iter().sum::<f64>(),
                ));
            }
            columns.push_str(&format!(
                ",{:.6e},{:.6e},{},{}",
                heating.populations(0).ok_or("missing left cell")?[0]
                    / (emitters / cell_count as f64),
                heating.populations(1).ok_or("missing right cell")?[0]
                    / (emitters / cell_count as f64),
                report.accepted_kinetic_substeps,
                report.coupling_iterations,
            ));
            columns
        } else if let Some(pump) = heating.depleting_pump_energy() {
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
        let pair_accounts = if annihilation {
            format!(",{:.9e}", heating.energy().annihilation_heat_j)
        } else {
            String::new()
        };
        let recycling_accounts = if recycling {
            format!(
                ",{:.9e},{:.9e},{:.9e}",
                heating.energy().recycled_excitation_j,
                heating.energy().secondary_absorption_heat_j,
                heating.energy().escaped_j.iter().sum::<f64>()
            )
        } else {
            String::new()
        };
        println!(
            "{target:.6e},{:.6e},{:.6e},{:.6e},{:.9e},{:.9e},{:.9e},{:.9e},{:.9e},{:.12},{:.9e},{:.9e},{:.3e}{pair_accounts}{recycling_accounts}{beam_accounts}",
            counts[0] / emitters,
            counts[1] / emitters,
            counts[2] / emitters,
            heating.absorbed_pump_energy_j().iter().sum::<f64>(),
            emission[0],
            emission[1],
            heating.energy().stored_excitation_j,
            heating.thermal_input_j(),
            mean_temperature_k,
            observed[0],
            observed[1],
            heating.energy().residual_j,
        );
    }
    Ok(())
}

fn main() -> Result<(), Box<dyn Error>> {
    let arguments: Vec<_> = std::env::args().skip(1).collect();
    let mut mode = PumpMode::Prescribed;
    let mut annihilation = false;
    let mut recycling = false;
    let mut scattered_recycling = false;
    let mut background_pump = false;
    for argument in &arguments {
        match argument.as_str() {
            "--depleting" if mode == PumpMode::Prescribed => mode = PumpMode::Depleting,
            "--counterpropagating" if mode == PumpMode::Prescribed => {
                mode = PumpMode::Counterpropagating
            }
            "--annihilation" if !annihilation => annihilation = true,
            "--recycling" if !recycling => recycling = true,
            "--scattered-recycling" if !scattered_recycling => scattered_recycling = true,
            "--background-pump" if !background_pump => background_pump = true,
            _ => return Err(
                "usage: saturable_afterglow [--depleting | --counterpropagating] [--annihilation] [--recycling] [--scattered-recycling] [--background-pump]"
                    .into(),
            ),
        }
    }
    if background_pump && mode == PumpMode::Prescribed {
        return Err("--background-pump requires --depleting or --counterpropagating".into());
    }
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
            mode,
            annihilation,
            recycling || scattered_recycling,
            scattered_recycling,
            background_pump,
        )
    })
}
