//! Species-selective vacuum blowdown using the original material gas owner.
//!
//! Run with GAS.tsv, exposed area [m2], aperture diameter [m], wall thickness
//! [m], qualified minimum mean free path [m], optional steps and dt [s].
//! TSV columns: species, mole fraction, molar mass [kg/mol], Cp/R constant,
//! Cp/R slope [1/K], viscosity [Pa s], conductivity [W/(m K)] at 200 K.
//! Inputs use the existing SYNTHETIC example format, not measured gas data.
//! The CLI mean-free-path bound must cover every retained mixture over the
//! entire 200–1800 K / 1000–1e7 Pa support; this program cannot establish it.
//! Fixed volume, stationary Maxwellian gas, independent thin apertures and
//! unit transmission are authored assumptions. No vacuum temperature or pump
//! model is inferred. Timesteps require refinement; rejected steps emit no row.

#[path = "material_gas_spring/source.rs"]
mod source;
use fs_material::thermochemical_gas::energy::parcel::transfer::species::{
    effusion::{EffusiveInterface, EffusiveRegime, EffusiveSpeciesTransmission, EffusiveSupport},
    network::{
        SpeciesNetworkSolve,
        transport::{
            GasTransportPaths, GasTransportSolve,
            effusion::vacuum::{VacuumEffusionLink, step_candidate},
        },
    },
};
use fs_material::thermochemical_gas::energy::{GasEnergySolve, parcel::ClosedGasParcel};
use std::{
    collections::BTreeMap,
    io::{Read, Write},
};
type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;

fn main() -> Result<()> {
    let args: Vec<_> = std::env::args().skip(1).take(8).collect();
    if !(5..=7).contains(&args.len()) {
        return Err("usage: effusion_blowdown GAS.tsv area_m2 aperture_diameter_m wall_thickness_m qualified_minimum_mean_free_path_m [steps=100] [dt_s=0.01]".into());
    }
    let area: f64 = args[1].parse()?;
    let diameter: f64 = args[2].parse()?;
    let thickness: f64 = args[3].parse()?;
    let mean_free_path: f64 = args[4].parse()?;
    let steps: u32 = args.get(5).map_or(Ok(100), |v| v.parse())?;
    let dt: f64 = args.get(6).map_or(Ok(0.01), |v| v.parse())?;
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
    writeln!(output)?;
    let (mut total_mass, mut total_energy) = (0.0, 0.0);
    for i in 0..=steps {
        let (mut escaped_mass, mut escaped_energy, mut mass_residual, mut energy_residual) =
            (0.0, 0.0, 0.0, 0.0);
        if i > 0 {
            let proposal = step_candidate(
                &gases,
                GasTransportPaths::default(),
                &links,
                dt,
                &endpoints,
                controls,
                || false,
            )?;
            let (next, step) = proposal.into_parts();
            for boundary in &step.exchange.boundaries {
                escaped_mass += boundary.transfer.outgoing_mass_kg;
                escaped_energy -= boundary.energy_into_gas_j;
            }
            mass_residual = step.exchange.total_mass_residual_kg;
            energy_residual = step.exchange.total_energy_residual_j;
            gases = next;
            total_mass += escaped_mass;
            total_energy += escaped_energy;
        }
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
            total_mass,
            total_energy,
            mass_residual,
            energy_residual
        )?;
        for fraction in gas.model().thermodynamics().mole_composition().fractions() {
            write!(output, ",{fraction:.12e}")?;
        }
        writeln!(output)?;
    }
    output.flush()?;
    Ok(())
}
