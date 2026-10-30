//! The same accepted pressure/thermal history for CLI and attached experiments.
//! File output is bounded before either artifact is opened; stdout keeps the
//! original streaming CSV path. No resampling or level normalization occurs.

use super::{History, Result, flow};
use std::{
    ffi::OsString,
    io::{self, Write},
    path::{Path, PathBuf},
};

pub(super) const MAX_HISTORY_BYTES: usize = 16 * 1024 * 1024;

#[derive(Default)]
pub(super) struct Options {
    pub(super) history: Option<PathBuf>,
    required_sample_rate_hz: Option<u32>,
}

impl Options {
    pub(super) fn extract(args: &[String]) -> Result<(Vec<OsString>, Self)> {
        let mut options = Self::default();
        let mut retained = args.iter().take(3).map(OsString::from).collect::<Vec<_>>();
        let mut rest = args.iter().skip(3);
        while let Some(arg) = rest.next() {
            match arg.as_str() {
                "--history-output" => {
                    if options.history.is_some() {
                        return Err("--history-output cannot be repeated".into());
                    }
                    let path = rest
                        .next()
                        .filter(|s| !s.is_empty() && !s.starts_with("--"))
                        .ok_or("--history-output needs a new CSV path")?;
                    options.history = Some(PathBuf::from(path));
                }
                "--require-sample-rate" => {
                    if options.required_sample_rate_hz.is_some() {
                        return Err("--require-sample-rate cannot be repeated".into());
                    }
                    let rate = rest
                        .next()
                        .ok_or("--require-sample-rate needs a rate in Hz")?
                        .parse::<u32>()?;
                    if !(8_000..=192_000).contains(&rate) {
                        return Err("required sample rate must be 8000..192000 Hz".into());
                    }
                    options.required_sample_rate_hz = Some(rate);
                }
                _ => retained.push(OsString::from(arg)),
            }
        }
        Ok((retained, options))
    }

    /// Assert an output capability without changing the authored physical clock.
    pub(super) fn admit(&self, sample_rate_hz: u32, wav: &Path) -> Result<()> {
        if let Some(required) = self.required_sample_rate_hz {
            if required != sample_rate_hz {
                return Err(format!("scenario sample rate {sample_rate_hz} Hz differs from --require-sample-rate {required} Hz; no resampling is performed").into());
            }
        }
        if self.history.as_deref() == Some(wav) {
            return Err("WAV and history output paths must differ".into());
        }
        for path in std::iter::once(wav).chain(self.history.as_deref()) {
            match std::fs::symlink_metadata(path) {
                Ok(_) => return Err(format!("output {} already exists", path.display()).into()),
                Err(error) if error.kind() == io::ErrorKind::NotFound => {}
                Err(error) => return Err(error.into()),
            }
        }
        Ok(())
    }
}

struct LimitedBytes {
    bytes: Vec<u8>,
    limit: usize,
}
impl Write for LimitedBytes {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        if self
            .bytes
            .len()
            .checked_add(bytes.len())
            .is_none_or(|size| size > self.limit)
        {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "trajectory CSV exceeds the history byte limit; shorten the declared run",
            ));
        }
        self.bytes.extend_from_slice(bytes);
        Ok(bytes.len())
    }

    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}

pub(super) fn csv_bytes(history: &History, limit: usize) -> Result<Vec<u8>> {
    let mut output = LimitedBytes {
        bytes: Vec::new(),
        limit,
    };
    write_csv(history, &mut output)?;
    Ok(output.bytes)
}

/// Every column observes the same accepted samples as the original WAV.
pub(super) fn write_csv(history: &History, mut output: impl Write) -> Result<()> {
    write!(
        output,
        "time_s,pressure_pa,displacement_m,velocity_m_s,solid_temperature_k,gas_temperature_k,gas_pressure_pa,material_heat_j,drag_heat_j,solid_to_gas_heat_j,joule_heat_j,prescribed_solid_heat_j,prescribed_gas_heat_j,mechanical_work_j,storage_change_j,energy_defect_j,energy_tolerance_j,current_a,endpoint_current_a,resistance_ohm,endpoint_resistance_ohm,first_mode_damping_ratio,tension_n,first_mode_frequency_hz,reversible_heat_j,axial_free_energy_change_j,mechanical_entropy_energy_change_j"
    )?;
    if let Some(temperatures) = history.layer_temperatures.first() {
        for layer in 0..temperatures.len() {
            write!(output, ",layer_{layer}_temperature_k")?;
        }
    }
    if let Some(gas_flow) = &history.gas_flow {
        for name in flow::COLUMNS {
            write!(output, ",{name}")?;
        }
        for species in &gas_flow.species {
            let escaped = species.replace('"', "\"\"");
            write!(
                output,
                ",\"gas_mole_fraction_{escaped}\",\"reservoir_mole_fraction_{escaped}\""
            )?;
        }
    }
    if history.electrical_feedback.is_some() {
        write!(
            output,
            ",cumulative_electrical_source_work_j,cumulative_electrical_source_residual_j,electrical_feedback_trials"
        )?;
    }
    writeln!(output)?;
    for (sample, record) in history.records.iter().enumerate() {
        for (i, value) in record.iter().enumerate() {
            if i != 0 {
                write!(output, ",")?;
            }
            write!(output, "{value:.17e}")?;
        }
        if let Some(temperatures) = history.layer_temperatures.get(sample) {
            for temperature in temperatures {
                write!(output, ",{temperature:.17e}")?;
            }
        }
        if let Some(gas_flow) = &history.gas_flow {
            for value in &gas_flow.records[sample] {
                write!(output, ",{value:.17e}")?;
            }
        }
        if let Some(feedback) = &history.electrical_feedback {
            for value in feedback[sample] {
                write!(output, ",{value:.17e}")?;
            }
        }
        writeln!(output)?;
    }
    output.flush()?;
    Ok(())
}
