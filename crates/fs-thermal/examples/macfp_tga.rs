//! Prescribed-temperature kinetic comparison against retained MaCFP FDS output.
//!
//! Run from the repository root:
//! `cargo run -p fs-thermal --example macfp_tga -- [CSV_PATH] [MAX_DT_S]`
//! Defaults are the retained NIST/MaCFP CSV and a maximum step of 0.1 s.
//! The CSV's material-temperature history is linearly interpolated, and the
//! original pyrolysis Arrhenius owner advances its first-order inventory with
//! backward Euler at each substep endpoint. No solid heat equation is solved.
//!
//! This is model-to-model kinetic verification, not an experimental validation
//! or an independent PMMA property fit. A=2.85e13 /s and E=191000 J/mol are the
//! NIST 2021 MaCFP submission's single first-order reaction. The retained source
//! describes calibration from TGA; using it here is not validation of that fit.
//! Source repository: https://github.com/MaCFP/matl-db
//! Retained revision: 9ae21bebd5dbc7478b8637fbc4df1c81291ec6ea.
//!
//! Despite its `[-]` mass heading, the source starts at 5. Both histories are
//! therefore normalized to their initial inventories; no physical mass unit is
//! guessed. The comparison stops at the first zero source mass because its
//! material-temperature history freezes after depletion. Our finite-rate
//! remaining inventory is reported as computed, without a burnout threshold.

use fs_thermal::enthalpy::pyrolysis::Arrhenius;
use std::error::Error;
use std::io::{self, BufWriter, Write};

const DEFAULT_INPUT: &str = "data/reference/compartment_fire/pyrolysis/macfp_nist_tga_fds.csv";
const MAX_FILE_BYTES: u64 = 1_048_576;
const MAX_ROWS: usize = 100_000;
const MAX_STEPS: usize = 10_000_000;

#[derive(Clone, Copy)]
struct Sample {
    time_s: f64,
    temperature_k: f64,
    mass: f64,
}

#[derive(Default)]
struct Errors {
    samples: usize,
    squared_sum: f64,
    maximum: f64,
}
impl Errors {
    fn add(&mut self, reference: f64, predicted: f64) {
        let error = (predicted - reference).abs();
        self.samples += 1;
        self.squared_sum += error * error;
        self.maximum = self.maximum.max(error);
    }
    fn print(&self, population: &str) {
        eprintln!(
            "{population}: samples={}, normalized_rmse={:.12e}, normalized_max_error={:.12e}",
            self.samples,
            (self.squared_sum / self.samples as f64).sqrt(),
            self.maximum,
        );
    }
}

fn read_samples(path: &str) -> Result<Vec<Sample>, Box<dyn Error>> {
    if std::fs::metadata(path)?.len() > MAX_FILE_BYTES {
        return Err("TGA input exceeds the 1 MiB example budget".into());
    }
    let text = std::fs::read_to_string(path)?;
    if text.len() as u64 > MAX_FILE_BYTES {
        return Err("TGA input exceeds the 1 MiB example budget".into());
    }
    let mut lines = text.lines();
    if lines.next().map(str::trim) != Some("t,T,m")
        || lines.next().map(str::trim) != Some("[s],[K],[-]")
    {
        return Err("expected the retained MaCFP t,T,m and [s],[K],[-] headings".into());
    }
    let mut samples: Vec<Sample> = Vec::new();
    for (index, line) in lines.enumerate() {
        if line.trim().is_empty() {
            continue;
        }
        if samples.len() == MAX_ROWS {
            return Err("TGA input exceeds the 100000-row example budget".into());
        }
        let fields: Vec<_> = line.split(',').map(str::trim).collect();
        if fields.len() != 3 {
            return Err(format!("TGA row {} must contain three numeric fields", index + 3).into());
        }
        let row = Sample {
            time_s: fields[0].parse()?,
            temperature_k: fields[1].parse()?,
            mass: fields[2].parse()?,
        };
        if !row.time_s.is_finite()
            || row.time_s < 0.0
            || !row.temperature_k.is_finite()
            || row.temperature_k <= 0.0
            || !row.mass.is_finite()
            || row.mass < 0.0
            || samples.last().is_some_and(|last| row.time_s <= last.time_s)
        {
            return Err(format!(
                "invalid temperature, mass or increasing time at TGA row {}",
                index + 3
            )
            .into());
        }
        samples.push(row);
    }
    if samples.len() < 2 || samples[0].mass <= 0.0 {
        return Err("need at least two TGA rows and a positive initial inventory".into());
    }
    Ok(samples)
}

fn run() -> Result<(), Box<dyn Error>> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    if args.len() > 2 {
        return Err("usage: macfp_tga [CSV_PATH] [MAX_DT_S]".into());
    }
    let input = args.first().map(String::as_str).unwrap_or(DEFAULT_INPUT);
    let maximum_dt: f64 = args.get(1).map(|s| s.parse()).transpose()?.unwrap_or(0.1);
    if !maximum_dt.is_finite() || maximum_dt <= 0.0 {
        return Err("MAX_DT_S must be positive and finite".into());
    }
    let samples = read_samples(input)?;
    let first_zero = samples.iter().position(|s| s.mass == 0.0);
    let stop_index = first_zero.unwrap_or(samples.len() - 1);
    let law = Arrhenius {
        prefactor_s_inv: 2.85e13,
        activation_energy_j_mol: 191_000.0,
    };
    let mut predicted = 1.0;
    let mut steps = 0usize;
    let mut positive_errors = Errors::default();
    let mut core_errors = Errors::default();
    let mut out = BufWriter::new(io::stdout().lock());
    writeln!(
        out,
        "time_s,temperature_k,reference_mass_fraction,predicted_mass_fraction,absolute_error"
    )?;
    for i in 0..=stop_index {
        let row = samples[i];
        if i > 0 {
            let previous = samples[i - 1];
            let duration = row.time_s - previous.time_s;
            let count = (duration / maximum_dt).ceil();
            if !count.is_finite() || count > (MAX_STEPS - steps) as f64 {
                return Err("requested TGA integration exceeds the 10000000-step budget".into());
            }
            let count = count as usize;
            let mut time = previous.time_s;
            for j in 1..=count {
                let fraction = j as f64 / count as f64;
                let endpoint = if j == count {
                    row.time_s
                } else {
                    previous.time_s + fraction * duration
                };
                let temperature = if j == count {
                    row.temperature_k
                } else {
                    previous.temperature_k + fraction * (row.temperature_k - previous.temperature_k)
                };
                predicted = law
                    .isothermal_step(predicted, temperature, endpoint - time)?
                    .virgin_mass_kg;
                time = endpoint;
            }
            steps += count;
        }
        let reference = row.mass / samples[0].mass;
        if reference > 0.0 {
            positive_errors.add(reference, predicted);
            if reference >= 0.01 {
                core_errors.add(reference, predicted);
            }
        }
        writeln!(
            out,
            "{:.9},{:.9},{:.15e},{:.15e},{:.15e}",
            row.time_s,
            row.temperature_k,
            reference,
            predicted,
            (predicted - reference).abs()
        )?;
    }
    out.flush()?;
    eprintln!(
        "scope=prescribed-temperature model-to-model kinetic verification; not experimental validation"
    );
    eprintln!(
        "source_revision=9ae21bebd5dbc7478b8637fbc4df1c81291ec6ea; A=2.85e13/s; E=191000 J/mol; order=1"
    );
    eprintln!(
        "reference_initial_value={}; normalized_initial_inventory=1; maximum_dt_s={maximum_dt}; actual_steps={steps}",
        samples[0].mass
    );
    positive_errors.print("strictly_positive_reference_mass");
    core_errors.print("reference_mass_fraction_at_least_0.01");
    if let Some(index) = first_zero {
        eprintln!(
            "first_reference_zero_time_s={}; predicted_remaining_fraction_at_first_zero={predicted:.15e}; no_burnout_clipping=true",
            samples[index].time_s
        );
        eprintln!(
            "comparison_stopped_at_first_zero=true; later_source_rows_excluded={}",
            samples.len() - index - 1
        );
    } else {
        eprintln!(
            "first_reference_zero=absent; final_time_s={}; predicted_remaining_fraction={predicted:.15e}",
            samples[stop_index].time_s
        );
    }
    Ok(())
}

fn main() {
    if let Err(error) = run() {
        eprintln!("macfp_tga: {error}");
        std::process::exit(1);
    }
}
