//! Finite-reservoir rocket performance from the existing source-backed gas APIs.
//!
//! From the repository root:
//! `cargo run -p fs-material --example rocket_blowdown -- --synthetic
//! crates/fs-material/examples/material_gas_spring/mixture-a.tsv 500 300000`
//!
//! The supplied example data are manufactured laws, not measured propellant
//! properties. `--packs` consumes the existing canonical gas-record/pack API.
//! See README.md for source limits, arguments and CSV time conventions.

#[path = "../material_gas_spring/source.rs"]
mod source;

use fs_material::gas_species::input::{
    ThermochemicalPhase,
    packs::{MAX_RECORD_BYTES, MAX_SOURCE_PACK_BYTES, SourcePackBytes, read_source_parcel},
};
use fs_material::thermochemical_gas::energy::{
    GasEnergySolve,
    entropy::GasEntropySolve,
    parcel::{
        ClosedGasParcel,
        transfer::nozzle::{
            blowdown::{GasBlowdownSolve, RocketBlowdownSolve},
            rocket::{RocketNozzle, RocketNozzleSolve},
        },
    },
};
use fs_qty::Temperature;
use std::{
    error::Error,
    io::{Read, Write},
    path::Path,
};

type Result<T> = std::result::Result<T, Box<dyn Error>>;

const USAGE: &str = "usage:
  rocket_blowdown --synthetic MIXTURE.tsv TEMPERATURE_K PRESSURE_PA [OPTIONS]
  rocket_blowdown --packs STATE.gas SPECIES.pack NASA9.pack TRANSPORT.pack [OPTIONS]
options (SI):
  --steps 100 --dt 0.02 --volume 0.002 --throat-area 0.000001
  --area-ratio 1.5 --ambient 0 --temperature-min 200 --temperature-max 1800
The synthetic route reuses the existing manufactured gas-spring source format.
The pack route accepts one canonical species/NASA-9/material-card triple and
an existing thermochemical gas record. All nozzle and tank states must remain
within the supplied source domains; overexpanded flow is refused.";

#[derive(Clone, Copy)]
struct Config {
    steps: usize,
    dt_s: f64,
    volume_m3: f64,
    throat_area_m2: f64,
    area_ratio: f64,
    ambient_pressure_pa: f64,
    minimum_temperature_k: f64,
    maximum_temperature_k: f64,
}
impl Default for Config {
    fn default() -> Self {
        Self {
            steps: 100,
            dt_s: 0.02,
            volume_m3: 0.002,
            throat_area_m2: 1e-6,
            area_ratio: 1.5,
            ambient_pressure_pa: 0.0,
            minimum_temperature_k: 200.0,
            maximum_temperature_k: 1800.0,
        }
    }
}
impl Config {
    fn parse(words: &[String]) -> Result<Self> {
        let mut config = Self::default();
        if words.len() % 2 != 0 {
            return Err("every option needs a value".into());
        }
        for pair in words.chunks_exact(2) {
            let value = &pair[1];
            match pair[0].as_str() {
                "--steps" => config.steps = value.parse()?,
                "--dt" => config.dt_s = value.parse()?,
                "--volume" => config.volume_m3 = value.parse()?,
                "--throat-area" => config.throat_area_m2 = value.parse()?,
                "--area-ratio" => config.area_ratio = value.parse()?,
                "--ambient" => config.ambient_pressure_pa = value.parse()?,
                "--temperature-min" => config.minimum_temperature_k = value.parse()?,
                "--temperature-max" => config.maximum_temperature_k = value.parse()?,
                other => return Err(format!("unknown option {other:?}\n{USAGE}").into()),
            }
        }
        config.validate()?;
        Ok(config)
    }

    fn validate(self) -> Result<()> {
        let positive = |x: f64| x.is_finite() && x > 0.0;
        if !(1..=10000).contains(&self.steps)
            || !positive(self.dt_s)
            || !(self.dt_s * self.steps as f64).is_finite()
            || !positive(self.volume_m3)
            || !positive(self.throat_area_m2)
            || !self.area_ratio.is_finite()
            || self.area_ratio < 1.0
            || !positive(self.throat_area_m2 * self.area_ratio)
            || !self.ambient_pressure_pa.is_finite()
            || self.ambient_pressure_pa < 0.0
            || !positive(self.minimum_temperature_k)
            || !positive(self.maximum_temperature_k)
            || self.minimum_temperature_k >= self.maximum_temperature_k
        {
            return Err("require 1..10000 steps, positive finite dt/volume/areas, area ratio >= 1, nonnegative ambient pressure and increasing positive temperature bounds".into());
        }
        Ok(())
    }

    fn nozzle(self) -> RocketNozzle {
        RocketNozzle {
            throat_area_m2: self.throat_area_m2,
            exit_area_m2: self.throat_area_m2 * self.area_ratio,
        }
    }

    fn controls(self) -> RocketBlowdownSolve {
        RocketBlowdownSolve {
            nozzle: RocketNozzleSolve {
                nozzle: GasEnergySolve {
                    minimum_temperature_k: self.minimum_temperature_k,
                    maximum_temperature_k: self.maximum_temperature_k,
                    energy_tolerance_j_kg: 1e-7,
                    maximum_iterations: 128,
                },
                relative_mass_flow_tolerance: 1e-9,
            },
            reservoir: GasBlowdownSolve {
                entropy: GasEntropySolve {
                    minimum_temperature_k: self.minimum_temperature_k,
                    maximum_temperature_k: self.maximum_temperature_k,
                    entropy_tolerance_j_kg_k: 1e-10,
                    maximum_iterations: 128,
                },
                mass_tolerance_kg: 1e-13,
            },
            maximum_donor_fraction: 0.05,
            maximum_iterations: 32,
        }
    }
}

fn read_bounded(path: &Path, limit: usize) -> Result<Vec<u8>> {
    let mut bytes = Vec::new();
    std::fs::File::open(path)?
        .take(limit as u64 + 1)
        .read_to_end(&mut bytes)?;
    if bytes.is_empty() || bytes.len() > limit {
        return Err(format!("{} must contain 1..{limit} bytes", path.display()).into());
    }
    Ok(bytes)
}

fn load(words: &[String]) -> Result<(ClosedGasParcel, Config)> {
    let mode = words.first().map(String::as_str).ok_or(USAGE)?;
    match mode {
        "--synthetic" if words.len() >= 4 => {
            let config = Config::parse(&words[4..])?;
            let text = String::from_utf8(read_bounded(Path::new(&words[1]), 65536)?)?;
            let model = source::load(&text, &words[1])?;
            let point = source::point(words[2].parse()?, words[3].parse()?)?;
            eprintln!(
                "Manufactured source laws: numerical demonstration, not calibrated propellant properties."
            );
            Ok((
                ClosedGasParcel::new(model, point, config.volume_m3)?,
                config,
            ))
        }
        "--packs" if words.len() >= 5 => {
            let config = Config::parse(&words[5..])?;
            let record = String::from_utf8(read_bounded(Path::new(&words[1]), MAX_RECORD_BYTES)?)?;
            let species = read_bounded(Path::new(&words[2]), MAX_SOURCE_PACK_BYTES)?;
            let nasa9 = read_bounded(Path::new(&words[3]), MAX_SOURCE_PACK_BYTES)?;
            let transport = read_bounded(Path::new(&words[4]), MAX_SOURCE_PACK_BYTES)?;
            let parcel = read_source_parcel(
                &record,
                &[SourcePackBytes {
                    species: &species,
                    nasa9: &nasa9,
                    transport: &transport,
                }],
                None,
                ThermochemicalPhase::FrozenGas,
                config.volume_m3,
            )?;
            Ok((ClosedGasParcel::from_isochoric(parcel), config))
        }
        _ => Err(USAGE.into()),
    }
}

#[derive(Default)]
struct Sum {
    value: f64,
    correction: f64,
}
impl Sum {
    fn add(&mut self, increment: f64) {
        let adjusted = increment - self.correction;
        let next = self.value + adjusted;
        self.correction = (next - self.value) - adjusted;
        self.value = next;
    }
}

struct Summary {
    reservoir: ClosedGasParcel,
    total_impulse_n_s: f64,
    mass_balance_residual_kg: f64,
    energy_balance_residual_j: f64,
}

fn stored_energy_change(initial: &ClosedGasParcel, current: &ClosedGasParcel) -> Result<f64> {
    let t0 = initial.state().state().temperature;
    let t1 = current.state().state().temperature;
    let dh = initial
        .model()
        .thermodynamics()
        .specific_enthalpy_change(Temperature::new(t0), Temperature::new(t1))?
        .value();
    let du = dh - initial.state().state().specific_gas_constant * (t1 - t0);
    Ok(current.mass_kg() * du
        + (current.mass_kg() - initial.mass_kg()) * initial.state().specific_internal_energy_j_kg())
}

fn simulate(initial: ClosedGasParcel, config: Config, mut output: impl Write) -> Result<Summary> {
    config.validate()?;
    let controls = config.controls();
    let nozzle = config.nozzle();
    let mut reservoir = initial.clone();
    let (mut expelled, mut exhausted, mut impulse) =
        (Sum::default(), Sum::default(), Sum::default());
    let mut flow = reservoir.rocket_nozzle_flow(config.ambient_pressure_pa, nozzle, controls.nozzle, || false)
        .map_err(|error| format!("initial nozzle refused: {error}; the source temperature and pressure domains must cover the tank, sonic throat and supersonic exit"))?;
    let (mut mass_residual, mut energy_residual) = (0.0, 0.0);
    writeln!(
        output,
        "time_s,flow_time_s,reservoir_mass_kg,reservoir_temperature_k,reservoir_pressure_pa,mass_flow_kg_s,exit_mach_number,exit_pressure_pa,thrust_n,specific_impulse_s,total_impulse_n_s,expelled_mass_kg,exhausted_enthalpy_j,mass_balance_residual_kg,energy_balance_residual_j,nozzle_relative_mass_flow_residual,nozzle_energy_residual_j_kg,step_mass_quadrature_residual_kg"
    )?;
    for i in 0..=config.steps {
        let mut quadrature = 0.0;
        if i > 0 {
            let candidate = reservoir.rocket_blowdown(config.ambient_pressure_pa, nozzle,
                config.dt_s, controls, || false).map_err(|error| format!(
                    "step {i} at t={} s refused: {error}; reduce dt for a step-bound refusal, or use source data and nozzle conditions that admit the full path",
                    (i - 1) as f64 * config.dt_s))?;
            quadrature = candidate.mass_quadrature_residual_kg;
            impulse.add(candidate.axial_impulse_n_s);
            let (next, step, midpoint_flow) = candidate.into_parts();
            expelled.add(step.expelled_mass_kg);
            exhausted.add(step.exhausted_enthalpy_j);
            reservoir = next;
            flow = midpoint_flow;
            mass_residual = (reservoir.mass_kg() - initial.mass_kg()) + expelled.value;
            energy_residual = stored_energy_change(&initial, &reservoir)? + exhausted.value;
        }
        let state = reservoir.state().state();
        let time = i as f64 * config.dt_s;
        let flow_time = if i == 0 {
            0.0
        } else {
            (i as f64 - 0.5) * config.dt_s
        };
        write_row(
            &mut output,
            [
                time,
                flow_time,
                reservoir.mass_kg(),
                state.temperature,
                state.pressure,
                flow.mass_flow_kg_s,
                flow.exit_mach_number,
                flow.exit().state().pressure,
                flow.thrust_n,
                flow.specific_impulse_s,
                impulse.value,
                expelled.value,
                exhausted.value,
                mass_residual,
                energy_residual,
                flow.relative_mass_flow_residual,
                flow.energy_residual_j_kg,
                quadrature,
            ],
        )?;
    }
    Ok(Summary {
        reservoir,
        total_impulse_n_s: impulse.value,
        mass_balance_residual_kg: mass_residual,
        energy_balance_residual_j: energy_residual,
    })
}

fn write_row(output: &mut impl Write, values: [f64; 18]) -> Result<()> {
    for (i, value) in values.into_iter().enumerate() {
        if i > 0 {
            write!(output, ",")?;
        }
        write!(output, "{value:.12e}")?;
    }
    writeln!(output)?;
    Ok(())
}

fn main() -> Result<()> {
    let words: Vec<_> = std::env::args().skip(1).take(65).collect();
    if words.len() > 64 {
        return Err("at most 64 command-line arguments are supported".into());
    }
    if words == ["--help"] {
        println!("{USAGE}");
        return Ok(());
    }
    let (initial, config) = load(&words)?;
    let stdout = std::io::stdout();
    let mut output = std::io::BufWriter::new(stdout.lock());
    let summary = simulate(initial, config, &mut output)?;
    output.flush()?;
    eprintln!(
        "Completed {:.6} s: remaining mass {:.9e} kg, temperature {:.6} K, pressure {:.6} Pa, impulse {:.9e} N s; mass residual {:.3e} kg, energy residual {:.3e} J.",
        config.steps as f64 * config.dt_s,
        summary.reservoir.mass_kg(),
        summary.reservoir.state().state().temperature,
        summary.reservoir.state().state().pressure,
        summary.total_impulse_n_s,
        summary.mass_balance_residual_kg,
        summary.energy_balance_residual_j
    );
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    const MIXTURE: &str = include_str!("../material_gas_spring/mixture-a.tsv");

    #[test]
    fn supplied_manufactured_data_runs_finite_blowdown_and_closes_inventory() {
        let config = Config::default();
        let initial = ClosedGasParcel::new(
            source::load(MIXTURE, "existing synthetic fixture").unwrap(),
            source::point(500.0, 300000.0).unwrap(),
            config.volume_m3,
        )
        .unwrap();
        let initial_mass = initial.mass_kg();
        let initial_flow = initial
            .rocket_nozzle_flow(0.0, config.nozzle(), config.controls().nozzle, || false)
            .unwrap();
        let mut csv = Vec::new();
        let summary = simulate(initial, config, &mut csv).unwrap();
        assert!(summary.reservoir.mass_kg() < initial_mass);
        assert!(summary.reservoir.state().state().temperature < 500.0);
        assert!(summary.reservoir.state().state().pressure < 300000.0);
        assert!(summary.total_impulse_n_s > 0.0);
        assert!(
            summary.total_impulse_n_s < initial_flow.thrust_n * config.dt_s * config.steps as f64
        );
        assert!(summary.mass_balance_residual_kg.abs() < 1e-13);
        assert!(summary.energy_balance_residual_j.abs() < 1e-6);
        assert_eq!(
            String::from_utf8(csv).unwrap().lines().count(),
            config.steps + 2
        );
    }
}
