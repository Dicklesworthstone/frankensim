//! Species-selective vacuum blowdown using the original material gas owner.
//!
//! Run with GAS.tsv, exposed area [m2], aperture diameter [m], wall thickness
//! [m], qualified minimum mean free path [m], optional steps, dt [s] and mode
//! (`fixed`, the default, or `adaptive`). In adaptive mode dt is the reporting
//! interval; the original finite-gas owner refines its internal timesteps.
//! TSV columns: species, mole fraction, molar mass [kg/mol], Cp/R constant,
//! Cp/R slope [1/K], viscosity [Pa s], conductivity [W/(m K)] at 200 K.
//! Inputs use the existing SYNTHETIC example format, not measured gas data.
//! The CLI mean-free-path bound must cover every retained mixture over the
//! entire 200–1800 K / 1000–1e7 Pa support; this program cannot establish it.
//! Fixed volume, stationary Maxwellian gas, independent thin apertures and
//! unit transmission are authored assumptions. No vacuum temperature or pump
//! model is inferred. Adaptivity controls a first-order local error estimate,
//! not a certified global error. Rejected steps contribute no output accounts.

#[path = "material_gas_spring/source.rs"]
mod source;
use fs_material::thermochemical_gas::energy::parcel::transfer::species::{
    effusion::{EffusiveInterface, EffusiveRegime, EffusiveSpeciesTransmission, EffusiveSupport},
    network::{
        SpeciesNetworkSolve, SpeciesNetworkStep,
        transport::{
            GasTransportPaths, GasTransportSolve,
            effusion::{
                trajectory::{EffusionAdvanceControls, advance_candidate},
                vacuum::{VacuumEffusionLink, step_candidate},
            },
        },
    },
};
use fs_material::thermochemical_gas::energy::{GasEnergySolve, parcel::ClosedGasParcel};
use fs_math::dd::Dd;
use std::{
    collections::BTreeMap,
    io::{Read, Write},
};
type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;

const USAGE: &str = "usage: effusion_blowdown GAS.tsv area_m2 aperture_diameter_m wall_thickness_m qualified_minimum_mean_free_path_m [steps=100] [dt_s=0.01] [mode=fixed|adaptive]";

// Aggregate only original accepted receipts, retaining small residuals when
// a reporting interval contains many fine substeps.
fn accumulate(accounts: &mut [Dd; 4], exchange: &SpeciesNetworkStep) -> Result<()> {
    for boundary in &exchange.boundaries {
        accounts[0] = accounts[0] + Dd::from_f64(boundary.transfer.outgoing_mass_kg);
        accounts[1] = accounts[1] - Dd::from_f64(boundary.energy_into_gas_j);
    }
    accounts[2] = accounts[2] + Dd::from_f64(exchange.total_mass_residual_kg);
    accounts[3] = accounts[3] + Dd::from_f64(exchange.total_energy_residual_j);
    if accounts.iter().any(|v| !v.to_f64().is_finite()) {
        return Err("accepted blowdown accounts exceed the finite output range".into());
    }
    Ok(())
}

fn main() -> Result<()> {
    let args: Vec<_> = std::env::args().skip(1).take(9).collect();
    if args == ["--help"] {
        println!("{USAGE}");
        return Ok(());
    }
    if !(5..=8).contains(&args.len()) {
        return Err(USAGE.into());
    }
    let area: f64 = args[1].parse()?;
    let diameter: f64 = args[2].parse()?;
    let thickness: f64 = args[3].parse()?;
    let mean_free_path: f64 = args[4].parse()?;
    let steps: u32 = args.get(5).map_or(Ok(100), |v| v.parse())?;
    let dt: f64 = args.get(6).map_or(Ok(0.01), |v| v.parse())?;
    let adaptive_mode = match args.get(7).map_or("fixed", String::as_str) {
        "fixed" => false,
        "adaptive" => true,
        _ => return Err("mode must be fixed or adaptive".into()),
    };
    if !area.is_finite()
        || area < 0.0
        || !(1..=10_000).contains(&steps)
        || !dt.is_finite()
        || dt <= 0.0
        || dt * f64::from(steps) > 3600.0
    {
        return Err(
            "nonnegative finite area, 1..=10000 steps, positive dt and at most 3600 s required"
                .into(),
        );
    }
    let mut text = String::new();
    std::fs::File::open(&args[0])?
        .take(65_537)
        .read_to_string(&mut text)?;
    let model = source::load(&text, &args[0])?;
    let transmissions = model
        .thermodynamics()
        .components()
        .iter()
        .map(|species| {
            EffusiveSpeciesTransmission::try_new(
                species.species().clone(),
                1.0,
                "authored open thin apertures; reciprocal and independent of molecular energy",
            )
        })
        .collect::<std::result::Result<Vec<_>, _>>()?;
    let provenance = format!(
        "authored CLI bound lambda_min={mean_free_path:.17e} m for all source mixtures at 200..1800 K, 1000..1e7 Pa; no empirical validation"
    );
    let interface = EffusiveInterface::try_new(
        model.clone(),
        transmissions,
        EffusiveSupport {
            temperature_k: (200.0, 1800.0),
            pressure_pa: (1000.0, 1e7),
        },
        EffusiveRegime::try_new(diameter, thickness, mean_free_path, provenance)?,
    )?;
    let gas = ClosedGasParcel::new(model, source::point(400.0, 100_000.0)?, 0.001)?;
    let mut gases = BTreeMap::from([(0, gas)]);
    let endpoints = BTreeMap::from([(
        0,
        GasEnergySolve {
            minimum_temperature_k: 200.0,
            maximum_temperature_k: 1800.0,
            energy_tolerance_j_kg: 1e-7,
            maximum_iterations: 128,
        },
    )]);
    let controls = GasTransportSolve {
        exchange: SpeciesNetworkSolve {
            maximum_species_fraction: 0.1,
            mass_tolerance_kg: 1e-13,
            energy_tolerance_j: 1e-8,
        },
        entropy_tolerance_j_k: 1e-10,
    };
    let adaptive = EffusionAdvanceControls {
        initial_step_s: dt,
        minimum_step_s: dt * 2.0_f64.powi(-24),
        maximum_step_s: dt,
        maximum_attempts: 10_000,
        maximum_accepted_intervals: 2048,
        relative_tolerance: 1e-5,
        temperature_tolerance_k: 1e-4,
        pressure_tolerance_pa: 0.025,
        // Temporal stream accuracy is separate from the stricter 1e-13 kg
        // conservation budget above; neither allowance scales with substeps.
        species_mass_tolerance_kg: 1e-11,
        stream_energy_tolerance_j: 1e-7,
    };
    let links = [VacuumEffusionLink {
        id: 0,
        parcel: 0,
        interface: &interface,
        area_m2: area,
    }];
    eprintln!(
        "SYNTHETIC source: {}; initial gas 0.001 m3, 400 K, 100000 Pa; geometry and collision bound are authored, not empirical.",
        args[0]
    );
    eprintln!("Qualification: {}", interface.regime().provenance());
    if adaptive_mode {
        eprintln!("Reporting interval: {dt} s; first-order local error controls: {adaptive:?}");
    }
    let stdout = std::io::stdout();
    let mut output = std::io::BufWriter::new(stdout.lock());
    write!(
        output,
        "time_s,pressure_pa,temperature_k,mass_kg,sound_speed_m_s,step_escaped_mass_kg,step_escaped_energy_j,total_escaped_mass_kg,total_escaped_energy_j,step_mass_residual_kg,step_energy_residual_j"
    )?;
    for species in gases[&0].model().thermodynamics().components() {
        write!(
            output,
            ",\"mole_fraction_{}\"",
            species.species().as_str().replace('"', "\"\"")
        )?;
    }
    writeln!(
        output,
        ",accepted_substeps,attempted_steps,rejected_intervals,maximum_accepted_error"
    )?;
    let (mut total_mass, mut total_energy) = (Dd::ZERO, Dd::ZERO);
    for i in 0..=steps {
        let mut accounts = [Dd::ZERO; 4];
        let mut counts = [0_usize; 3];
        let mut maximum_error = 0.0;
        if i > 0 {
            let context = |error: &dyn std::fmt::Display| {
                format!(
                    "blowdown interval {i} refused at t={} s: {error}",
                    f64::from(i - 1) * dt
                )
            };
            let next = if adaptive_mode {
                let proposal = advance_candidate(
                    &gases,
                    GasTransportPaths::default(),
                    &links,
                    &[],
                    &[],
                    dt,
                    &endpoints,
                    controls,
                    adaptive,
                    || false,
                )
                .map_err(|error| context(&error))?;
                let (next, report) = proposal.into_parts();
                counts = [
                    report.accepted_steps.len(),
                    report.attempted_steps,
                    report.rejected_intervals,
                ];
                maximum_error = report.maximum_accepted_error;
                for step in &report.accepted_steps {
                    accumulate(&mut accounts, &step.exchange)?;
                }
                next
            } else {
                let proposal = step_candidate(
                    &gases,
                    GasTransportPaths::default(),
                    &links,
                    dt,
                    &endpoints,
                    controls,
                    || false,
                )
                .map_err(|error| context(&error))?;
                let (next, step) = proposal.into_parts();
                counts = [1, 1, 0];
                accumulate(&mut accounts, &step.exchange)?;
                next
            };
            gases = next;
            total_mass = total_mass + accounts[0];
            total_energy = total_energy + accounts[1];
            if !total_mass.to_f64().is_finite() || !total_energy.to_f64().is_finite() {
                return Err("cumulative blowdown accounts exceed the finite output range".into());
            }
        }
        let [escaped_mass, escaped_energy, mass_residual, energy_residual] =
            accounts.map(Dd::to_f64);
        let gas = &gases[&0];
        let state = gas.state().state();
        write!(
            output,
            "{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e}",
            f64::from(i) * dt,
            state.pressure,
            state.temperature,
            gas.mass_kg(),
            state.sound_speed,
            escaped_mass,
            escaped_energy,
            total_mass.to_f64(),
            total_energy.to_f64(),
            mass_residual,
            energy_residual
        )?;
        for fraction in gas.model().thermodynamics().mole_composition().fractions() {
            write!(output, ",{fraction:.12e}")?;
        }
        writeln!(
            output,
            ",{},{},{},{maximum_error:.12e}",
            counts[0], counts[1], counts[2]
        )?;
    }
    output.flush()?;
    Ok(())
}
