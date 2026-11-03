//! Finite source gas -> molecular vacuum escape -> vessel recoil and exhaust.
//!
//! Run from the repository root:
//! `cargo run -p fs-couple --example effusion_recoil --
//! crates/fs-material/examples/material_gas_spring/mixture-a.tsv
//! 0.01 1e-7 1e-9 0 1e-7 100 0.01`
//!
//! GAS.tsv uses the existing SYNTHETIC seven-column example format: species,
//! mole fraction, molar mass [kg/mol], Cp/R constant, Cp/R slope [1/K],
//! viscosity [Pa s], conductivity [W/(m K)] at 200 K. No measured accuracy is
//! claimed. Initial gas is 0.001 m3 at 400 K and 100000 Pa; the vessel starts at
//! rest, with a +x outlet and no gravity or external force. Recoil points -x.
//! The authored minimum mean free path must cover every evolving mixture over
//! 200..1800 K and 1000..1e7 Pa. Unit, angle-independent, nondeflecting
//! transmission is assumed. Time refinement is required; refused steps emit
//! no endpoint row. Exported momentum and energy use the vessel's inertial frame.

#[path = "../../fs-material/examples/material_gas_spring/source.rs"]
mod source;
use fs_couple::rocket_blowdown::{
    GasRocketState,
    effusion::{EffusionOutlet, EffusionRocketSolve},
};
use fs_material::thermochemical_gas::energy::parcel::transfer::species::{
    effusion::{
        EffusiveInterface, EffusiveRegime, EffusiveSpeciesTransmission, EffusiveSupport,
        momentum::EffusiveMomentumInterface,
    },
    network::{SpeciesNetworkSolve, transport::GasTransportSolve},
};
use fs_material::thermochemical_gas::energy::{GasEnergySolve, parcel::ClosedGasParcel};
use fs_mbd::Vec3;
use std::io::{Read, Write};
type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;
const USAGE: &str = "usage: effusion_recoil GAS.tsv dry_mass_kg exposed_area_m2 aperture_diameter_m wall_thickness_m qualified_minimum_mean_free_path_m [steps=100] [dt_s=0.01]";

fn main() -> Result<()> {
    let args: Vec<_> = std::env::args().skip(1).take(9).collect();
    if args == ["--help"] {
        println!("{USAGE}");
        return Ok(());
    }
    if !(6..=8).contains(&args.len()) {
        return Err(USAGE.into());
    }
    let dry_mass: f64 = args[1].parse()?;
    let area: f64 = args[2].parse()?;
    let diameter: f64 = args[3].parse()?;
    let thickness: f64 = args[4].parse()?;
    let mean_free_path: f64 = args[5].parse()?;
    let steps: u32 = args.get(6).map_or(Ok(100), |v| v.parse())?;
    let dt: f64 = args.get(7).map_or(Ok(0.01), |v| v.parse())?;
    if !dry_mass.is_finite()
        || dry_mass <= 0.0
        || !area.is_finite()
        || area < 0.0
        || !(1..=10_000).contains(&steps)
        || !dt.is_finite()
        || dt <= 0.0
        || dt * f64::from(steps) > 3600.0
    {
        return Err("positive finite dry mass, nonnegative finite area, 1..=10000 steps, positive finite dt and at most 3600 s required".into());
    }
    let mut text = String::new();
    std::fs::File::open(&args[0])?
        .take(65_537)
        .read_to_string(&mut text)?;
    let model = source::load(&text, &args[0])?; // Original loader enforces 64 KiB / 32 species.
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
    let qualification = format!(
        "authored CLI bound lambda_min={mean_free_path:.17e} m over all source mixtures at 200..1800 K, 1000..1e7 Pa; no empirical validation"
    );
    let interface = EffusiveMomentumInterface::try_new(
        EffusiveInterface::try_new(
            model.clone(),
            transmissions,
            EffusiveSupport {
                temperature_k: (200.0, 1800.0),
                pressure_pa: (1000.0, 1e7),
            },
            EffusiveRegime::try_new(diameter, thickness, mean_free_path, qualification)?,
        )?,
        "authored uncollimated angle-independent nondeflecting transmission; no wall scattering or measured validation",
    )?;
    let gas = ClosedGasParcel::new(model, source::point(400.0, 100_000.0)?, 0.001)?;
    let mut vessel = GasRocketState::new(gas, dry_mass, Vec3::ZERO, Vec3::ZERO, 0)?;
    let outlets = [EffusionOutlet {
        id: 0,
        interface: &interface,
        area_m2: area,
        outward_normal_world: Vec3::new(1.0, 0.0, 0.0),
    }];
    let endpoint = GasEnergySolve {
        minimum_temperature_k: 200.0,
        maximum_temperature_k: 1800.0,
        energy_tolerance_j_kg: 1e-7,
        maximum_iterations: 128,
    };
    let controls = EffusionRocketSolve {
        transport: GasTransportSolve {
            exchange: SpeciesNetworkSolve {
                maximum_species_fraction: 0.1,
                mass_tolerance_kg: 1e-13,
                energy_tolerance_j: 1e-8,
            },
            entropy_tolerance_j_k: 1e-10,
        },
        momentum_tolerance_n_s: 1e-12,
    };
    eprintln!(
        "SYNTHETIC source: {}; 0.001 m3, 400 K, 100000 Pa; +x vacuum emission, free -x recoil, dry mass {dry_mass} kg.",
        args[0]
    );
    eprintln!(
        "Qualification: {}",
        interface.interface().regime().provenance()
    );
    eprintln!("Angular qualification: {}", interface.provenance());
    let stdout = std::io::stdout();
    let mut output = std::io::BufWriter::new(stdout.lock());
    write!(
        output,
        "time_s,temperature_k,pressure_pa,gas_mass_kg,position_x_m,velocity_x_m_s,step_escaped_mass_kg,total_escaped_mass_kg,step_escaped_momentum_x_kg_m_s,total_escaped_momentum_x_kg_m_s,step_escaped_energy_j,total_escaped_energy_j,step_mass_residual_kg,step_momentum_residual_x_kg_m_s,step_energy_residual_j"
    )?;
    for species in vessel.reservoir().model().thermodynamics().components() {
        write!(
            output,
            ",\"mole_fraction_{}\"",
            species.species().as_str().replace('"', "\"\"")
        )?;
    }
    writeln!(output)?;
    let (mut totals, mut compensation) = ([0.0; 3], [0.0; 3]);
    for i in 0..=steps {
        let (mut escaped, mut residuals) = ([0.0; 3], [0.0; 3]);
        if i > 0 {
            let candidate = vessel
                .advance_effusion(0, &outlets, dt, endpoint, controls, || false)
                .map_err(|error| {
                    format!(
                        "joint gas/vessel step {i} refused at t={} s: {error}",
                        f64::from(i - 1) * dt
                    )
                })?;
            let (next, step) = candidate.into_parts();
            escaped = [
                step.expelled_mass_kg,
                step.exhaust_momentum_world_kg_m_s.x,
                step.exhaust_energy_world_j,
            ];
            residuals = [
                step.gas.exchange.total_mass_residual_kg,
                step.momentum_residual_kg_m_s.x,
                step.energy_residual_j,
            ];
            for k in 0..3 {
                let increment = escaped[k] - compensation[k];
                let total = totals[k] + increment;
                compensation[k] = (total - totals[k]) - increment;
                totals[k] = total;
            }
            vessel = next;
        }
        let gas = vessel.reservoir();
        let state = gas.state().state();
        let row = [
            f64::from(i) * dt,
            state.temperature,
            state.pressure,
            gas.mass_kg(),
            vessel.vehicle().position_m.x,
            vessel.vehicle().velocity_m_s.x,
            escaped[0],
            totals[0],
            escaped[1],
            totals[1],
            escaped[2],
            totals[2],
            residuals[0],
            residuals[1],
            residuals[2],
        ];
        write!(output, "{:.12e}", row[0])?;
        for value in &row[1..] {
            write!(output, ",{value:.12e}")?;
        }
        for fraction in gas.model().thermodynamics().mole_composition().fractions() {
            write!(output, ",{fraction:.12e}")?;
        }
        writeln!(output)?;
    }
    output.flush()?;
    Ok(())
}
