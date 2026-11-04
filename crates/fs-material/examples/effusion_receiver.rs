//! Adaptive molecular filling of an exactly evacuated finite chamber.
//!
//! Run with GAS.tsv, exposed area [m2], aperture diameter [m], wall thickness
//! [m], qualified minimum mean free path [m], authored synthetic minimum
//! pressure [Pa], optional reporting steps and reporting dt [s]. The original
//! source loader's seven TSV columns are unchanged: species, mole fraction,
//! molar mass [kg/mol], Cp/R, Cp/R slope [1/K], viscosity [Pa s], conductivity
//! [W/(m K)] at 200 K. These are synthetic sources, not measured gas data.
//!
//! Both rigid volumes are 0.001 m3. The donor starts at 400 K and 100000 Pa;
//! the receiver starts with exactly zero gas. The first reporting interval
//! adaptively resolves capture and receiver backflow; subsequent intervals
//! use ordinary reciprocal finite-gas transport. No pressure floor is applied.
//!
//! The authored collision bound must cover both evolving mixtures over
//! 200..1800 K and minimum_pressure..1e7 Pa. Unit transmission, independent
//! thin apertures, complete retained-beam thermalization, zero bulk momentum,
//! fixed walls and no external heat/work are explicit assumptions. Wall inertia,
//! molecular flight delays, nonequilibrium startup and accommodation are absent.
//! Adaptivity estimates local error; it is not a certified global error bound.

#[allow(dead_code)] // Other examples use the shared loader's default entry point.
#[path = "material_gas_spring/source.rs"]
mod source;

use fs_matdb::QueryPoint;
use fs_material::thermochemical_gas::energy::{
    GasEnergySolve,
    parcel::{
        ClosedGasParcel,
        evacuated::EvacuatedGasVolume,
        transfer::species::{
            effusion::{
                EffusiveInterface, EffusiveRegime, EffusiveSpeciesTransmission, EffusiveSupport,
            },
            network::{
                SpeciesNetworkSolve,
                transport::{
                    GasTransportPaths, GasTransportSolve,
                    effusion::{
                        EffusionNetworkLink,
                        filling::{
                            self, EffusionFillBalance, EffusionFillLink, EffusionReceiver,
                            trajectory as filling_trajectory,
                        },
                        trajectory::{self, EffusionAdvanceControls},
                    },
                },
            },
        },
    },
};
use fs_thermochem::SpeciesId;
use std::{
    collections::BTreeMap,
    io::{Read, Write},
};

type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;
const VOLUME_M3: f64 = 0.001;
const USAGE: &str = "usage: effusion_receiver GAS.tsv area_m2 aperture_diameter_m wall_thickness_m qualified_minimum_mean_free_path_m synthetic_minimum_pressure_pa [steps=100] [dt_s=0.01]";

#[derive(Default)]
struct Work {
    accepted_substeps: usize,
    accepted_intervals: usize,
    attempted_steps: usize,
    rejected_intervals: usize,
    maximum_accepted_error: f64,
}

fn gas_columns(output: &mut impl Write, gas: Option<&ClosedGasParcel>) -> Result<()> {
    if let Some(gas) = gas {
        let state = gas.state().state();
        write!(
            output,
            ",{:.12e},{:.12e},{:.12e},{:.12e}",
            state.pressure,
            state.temperature,
            gas.mass_kg(),
            state.sound_speed
        )?;
    } else {
        // Exact vacuum has zero pressure/mass, but no gas temperature or sound
        // speed. Empty fields preserve that distinction in ordinary CSV readers.
        write!(output, ",0,,0,")?;
    }
    Ok(())
}

fn composition_columns(
    output: &mut impl Write,
    gas: Option<&ClosedGasParcel>,
    species: &[SpeciesId],
) -> Result<()> {
    for id in species {
        if let Some(gas) = gas {
            let thermo = gas.model().thermodynamics();
            let fraction = thermo
                .components()
                .binary_search_by(|s| s.species().cmp(id))
                .map_or(0.0, |i| thermo.mole_composition().fractions()[i]);
            write!(output, ",{fraction:.12e}")?;
        } else {
            write!(output, ",")?;
        }
    }
    Ok(())
}

fn main() -> Result<()> {
    let args: Vec<_> = std::env::args().skip(1).take(9).collect();
    if args == ["--help"] {
        println!("{USAGE}");
        return Ok(());
    }
    if !(6..=8).contains(&args.len()) {
        return Err(USAGE.into());
    }
    let area: f64 = args[1].parse()?;
    let diameter: f64 = args[2].parse()?;
    let thickness: f64 = args[3].parse()?;
    let mean_free_path: f64 = args[4].parse()?;
    let minimum_pressure: f64 = args[5].parse()?;
    let steps: u32 = args.get(6).map_or(Ok(100), |v| v.parse())?;
    let dt: f64 = args.get(7).map_or(Ok(0.01), |v| v.parse())?;
    let minimum_step = dt * 2.0_f64.powi(-24);
    if !area.is_finite()
        || area < 0.0
        || !minimum_pressure.is_finite()
        || minimum_pressure <= 0.0
        || minimum_pressure > 100000.0
        || !(1..=10000).contains(&steps)
        || !dt.is_finite()
        || dt <= 0.0
        || minimum_step <= 0.0
        || dt * f64::from(steps) > 3600.0
    {
        return Err("nonnegative finite area, authored 0 < minimum pressure <= 100000 Pa, 1..=10000 reporting steps, positive representable dt and at most 3600 s required".into());
    }
    let mut text = String::new();
    std::fs::File::open(&args[0])?
        .take(65537)
        .read_to_string(&mut text)?;
    let model = source::load_with_pressure_support(&text, &args[0], (minimum_pressure, 1e7))?;
    let species: Vec<_> = model
        .thermodynamics()
        .components()
        .iter()
        .map(|s| s.species().clone())
        .collect();
    let transmissions = species
        .iter()
        .map(|id| {
            EffusiveSpeciesTransmission::try_new(
                id.clone(),
                1.0,
                "authored open thin apertures; reciprocal and independent of molecular energy",
            )
        })
        .collect::<std::result::Result<Vec<_>, _>>()?;
    let provenance = format!(
        "authored CLI bound lambda_min={mean_free_path:.17e} m for both evolving mixtures at 200..1800 K and {minimum_pressure:.17e}..1e7 Pa; no empirical validation"
    );
    let interface = EffusiveInterface::try_new(
        model.clone(),
        transmissions,
        EffusiveSupport {
            temperature_k: (200.0, 1800.0),
            pressure_pa: (minimum_pressure, 1e7),
        },
        EffusiveRegime::try_new(diameter, thickness, mean_free_path, provenance)?,
    )?;
    let original = ClosedGasParcel::new(model, source::point(400.0, 100000.0)?, VOLUME_M3)?;
    let mut gases = BTreeMap::from([(0, original.clone())]);
    let mut empty = Some(EvacuatedGasVolume::new(VOLUME_M3, QueryPoint::new())?);
    let solve = GasEnergySolve {
        minimum_temperature_k: 200.0,
        maximum_temperature_k: 1800.0,
        energy_tolerance_j_kg: 1e-7,
        maximum_iterations: 128,
    };
    let endpoints = BTreeMap::from([(0, solve), (1, solve)]);
    let physical = GasTransportSolve {
        exchange: SpeciesNetworkSolve {
            maximum_species_fraction: 0.1,
            mass_tolerance_kg: 1e-13,
            energy_tolerance_j: 1e-8,
        },
        entropy_tolerance_j_k: 1e-10,
    };
    let adaptive = EffusionAdvanceControls {
        initial_step_s: dt,
        minimum_step_s: minimum_step,
        maximum_step_s: dt,
        maximum_attempts: 10000,
        maximum_accepted_intervals: 2048,
        relative_tolerance: 1e-5,
        temperature_tolerance_k: 1e-4,
        pressure_tolerance_pa: 0.025,
        species_mass_tolerance_kg: 1e-11,
        stream_energy_tolerance_j: 1e-7,
    };
    let fill_link = EffusionFillLink {
        id: 0,
        donor: 0,
        receiver: 1,
        interface: &interface,
        area_m2: area,
    };
    let links = [EffusionNetworkLink {
        id: 0,
        parcels: [0, 1],
        interface: &interface,
        area_m2: area,
    }];
    let paths = GasTransportPaths {
        effusion: &links,
        ..GasTransportPaths::default()
    };
    eprintln!(
        "SYNTHETIC source: {}; pressure-independent curves authored over {minimum_pressure}..1e7 Pa. Unsupported positive receiver pressures refuse; no floor is applied.",
        args[0]
    );
    eprintln!(
        "Two rigid 0.001 m3 chambers; donor 400 K/100000 Pa, receiver exact vacuum; complete molecular-energy thermalization, zero bulk momentum and no external heat/work are authored assumptions."
    );
    eprintln!("Qualification: {}", interface.regime().provenance());
    eprintln!(
        "Reporting dt: {dt} s; local error controls: {adaptive:?}. Pair residuals use the original donor plus vacuum over the complete elapsed time; entropy is physical change, not a zero-target residual."
    );
    let stdout = std::io::stdout();
    let mut output = std::io::BufWriter::new(stdout.lock());
    write!(
        output,
        "time_s,donor_pressure_pa,donor_temperature_k,donor_mass_kg,donor_sound_speed_m_s,receiver_state,receiver_pressure_pa,receiver_temperature_k,receiver_mass_kg,receiver_sound_speed_m_s,pair_mass_residual_kg,pair_maximum_species_mass_residual_kg,pair_energy_residual_j,pair_entropy_change_j_k,accepted_substeps,accepted_intervals,attempted_steps,rejected_intervals,maximum_accepted_error"
    )?;
    for side in ["donor", "receiver"] {
        for id in &species {
            write!(
                output,
                ",\"{side}_mole_fraction_{}\"",
                id.as_str().replace('"', "\"\"")
            )?;
        }
    }
    writeln!(output)?;
    let mut balance = EffusionFillBalance {
        total_mass_residual_kg: 0.0,
        maximum_species_mass_residual_kg: 0.0,
        energy_residual_j: 0.0,
        entropy_change_j_k: 0.0,
    };
    for i in 0..=steps {
        let mut work = Work::default();
        if i > 0 {
            let context = |error: &dyn std::fmt::Display| {
                format!(
                    "receiver interval {i} refused at t={} s: {error}",
                    f64::from(i - 1) * dt
                )
            };
            if let Some(receiver) = empty.as_ref() {
                let candidate = filling_trajectory::advance_candidate(
                    &gases[&0],
                    receiver,
                    fill_link,
                    dt,
                    solve,
                    solve,
                    physical,
                    adaptive,
                    || false,
                )
                .map_err(|e| context(&e))?;
                let (donor, receiver, report) = candidate.into_parts();
                work = Work {
                    accepted_substeps: 1
                        + usize::from(report.first_reciprocal.is_some())
                        + report
                            .continuation
                            .as_ref()
                            .map_or(0, |r| r.accepted_steps.len()),
                    accepted_intervals: report.accepted_intervals,
                    attempted_steps: report.attempted_steps,
                    rejected_intervals: report.rejected_intervals,
                    maximum_accepted_error: report.maximum_accepted_error,
                };
                gases.insert(0, donor);
                empty = match receiver {
                    EffusionReceiver::Evacuated(volume) => Some(volume),
                    EffusionReceiver::Filled(gas) => {
                        gases.insert(1, gas);
                        None
                    }
                };
            } else {
                let candidate = trajectory::advance_candidate(
                    &gases,
                    paths,
                    &[],
                    &[],
                    &[],
                    dt,
                    &endpoints,
                    physical,
                    adaptive,
                    || false,
                )
                .map_err(|e| context(&e))?;
                let (next, report) = candidate.into_parts();
                work = Work {
                    accepted_substeps: report.accepted_steps.len(),
                    accepted_intervals: report.accepted_intervals,
                    attempted_steps: report.attempted_steps,
                    rejected_intervals: report.rejected_intervals,
                    maximum_accepted_error: report.maximum_accepted_error,
                };
                gases = next;
            }
            balance = filling::audit_pair(
                &original,
                (0, &gases[&0]),
                gases.get(&1).map(|gas| (1, gas)),
                physical,
                &mut || false,
            )
            .map_err(|e| context(&e))?;
        }
        write!(output, "{:.12e}", f64::from(i) * dt)?;
        gas_columns(&mut output, Some(&gases[&0]))?;
        write!(
            output,
            ",{}",
            if empty.is_some() {
                "evacuated"
            } else {
                "filled"
            }
        )?;
        gas_columns(&mut output, gases.get(&1))?;
        write!(
            output,
            ",{:.12e},{:.12e},{:.12e},{:.12e},{},{},{},{},{:.12e}",
            balance.total_mass_residual_kg,
            balance.maximum_species_mass_residual_kg,
            balance.energy_residual_j,
            balance.entropy_change_j_k,
            work.accepted_substeps,
            work.accepted_intervals,
            work.attempted_steps,
            work.rejected_intervals,
            work.maximum_accepted_error
        )?;
        composition_columns(&mut output, Some(&gases[&0]), &species)?;
        composition_columns(&mut output, gases.get(&1), &species)?;
        writeln!(output)?;
    }
    output.flush()?;
    Ok(())
}
