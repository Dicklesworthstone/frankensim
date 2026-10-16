//! Opt-in computed broadband receivers over the existing v1 modal-file parser
//! and physical-force compiler. No second force scheduler or instrument model.
use super::{
    ForceInitialization, MAX_MODAL_PERFORMANCE_BYTES, MODAL_PERFORMANCE_SCHEMA,
    ModalPerformanceInfo, ParsedPerformance, Reader, input,
};
use crate::broadband_radiation::modal::admit_modal_samples;
use crate::modal_acoustic_time::broadband::{
    BroadbandModalControls, BroadbandModalObserver, BroadbandModalReport, ModalTransferSamples,
};
use crate::pcm_wav::observation::PressureRenderer;
use crate::render::{
    RenderError,
    schedule::force::program::{PhysicalForceInput, ProjectedForceEvent, ProjectedForceProgram},
};
use fs_math::c64::C64;

/// Existing mechanical performance with one declared receiver per source voice.
/// The appended block's exact bytes, including computed transfer samples, are
/// bound into the input hash. All sources retain EstimateOnly authority.
pub struct BroadbandModalPerformance {
    info: ModalPerformanceInfo,
    observers: Vec<BroadbandModalObserver>,
    reports: Vec<BroadbandModalReport>,
    forces: Vec<Vec<f64>>,
    events: Vec<ProjectedForceEvent>,
    next_event: usize,
    completed: u64,
    max_block: usize,
    poisoned: bool,
}
impl BroadbandModalPerformance {
    /// Detect only the explicit opt-in record. Invalid UTF-8 still refuses in
    /// the ordinary parser; a missing block retains the legacy renderer.
    pub fn has_section(bytes: &[u8]) -> bool {
        std::str::from_utf8(bytes).is_ok_and(|s| s.contains("\nbroadband "))
    }
    /// Read a v1 modal performance followed by `broadband VOICE_COUNT` and,
    /// for each voice, `transfer PRODUCER DELAY_S PRODUCER_REL_ERROR ORDER
    /// ITERATIONS DELAY_HALF MAX_REL_ERROR TRAINING_COUNT HELD_OUT_COUNT`, then
    /// mode-complete `frequency OMEGA H_RE H_IM ...` and `withheld ...` rows.
    /// Counts, full-file bytes and operator budgets are checked before fitting.
    /// Coupled/nonlinear performances refuse: their within-step trajectories
    /// cannot be replaced by this independent-mode exponential.
    pub fn from_bytes(bytes: &[u8], max_block: usize) -> Result<Self, String> {
        if bytes.len() > MAX_MODAL_PERFORMANCE_BYTES {
            return Err("broadband input exceeds 4 MiB".into());
        }
        let text = std::str::from_utf8(bytes).map_err(|_| "broadband input must be UTF-8")?;
        let offset = text
            .find("\nbroadband ")
            .ok_or("explicit broadband block required")?;
        let parsed = ParsedPerformance::from_bytes(&bytes[..offset], max_block)
            .map_err(|e| e.to_string())?;
        if parsed.info.schema != MODAL_PERFORMANCE_SCHEMA || parsed.coupled.is_some() {
            return Err("broadband receiver requires independent exact-ZOH v1 mechanics; coupled/nonlinear trajectories need their own within-step observer".into());
        }
        let mut reader = Reader {
            lines: text[offset + 1..].lines(),
            line: text[..offset].lines().count(),
        };
        let count: usize = reader.one("broadband").map_err(|e| e.to_string())?;
        if count != parsed.voices.len() {
            return Err("one complete transfer block is required per mechanical voice".into());
        }
        let mut bundles = Vec::with_capacity(count);
        let mut settings = Vec::with_capacity(count);
        for voice in &parsed.voices {
            let modes = voice.model.modes().len();
            let mut row = reader.row("transfer").map_err(|e| e.to_string())?;
            let source_id = row.word().map_err(|e| e.to_string())?.to_owned();
            let propagation_delay_s = row.scalar().map_err(|e| e.to_string())?;
            let producer_relative_error_estimate = row.scalar().map_err(|e| e.to_string())?;
            let controls = BroadbandModalControls {
                fit_order: row.parse().map_err(|e| e.to_string())?,
                fit_iterations: row.parse().map_err(|e| e.to_string())?,
                delay_kernel_half_len: row.parse().map_err(|e| e.to_string())?,
                max_held_out_rel_error: row.scalar().map_err(|e| e.to_string())?,
            };
            let training: usize = row.parse().map_err(|e| e.to_string())?;
            let held: usize = row.parse().map_err(|e| e.to_string())?;
            row.finish().map_err(|e| e.to_string())?;
            admit_modal_samples(modes, training, held).map_err(|e| e.to_string())?;
            let mut omega_rad_s = Vec::with_capacity(training);
            let mut per_mode = vec![Vec::with_capacity(training); modes];
            for _ in 0..training {
                let (omega, values) = read_frequency(&mut reader, "frequency", modes)?;
                omega_rad_s.push(omega);
                for (column, value) in per_mode.iter_mut().zip(values) {
                    column.push(value);
                }
            }
            let mut held_out = Vec::with_capacity(held);
            for _ in 0..held {
                held_out.push(read_frequency(&mut reader, "withheld", modes)?);
            }
            bundles.push(ModalTransferSamples {
                source_id,
                propagation_delay_s,
                producer_relative_error_estimate,
                omega_rad_s,
                per_mode,
                held_out,
            });
            settings.push(controls);
        }
        if reader.lines.next().is_some() {
            return Err(input(reader.line + 1, "unexpected broadband trailing record").to_string());
        }
        let mut inputs = Vec::with_capacity(count);
        let mut models = Vec::with_capacity(count);
        for voice in parsed.voices {
            inputs.push(
                PhysicalForceInput::new(voice.port_shapes, voice.initial_force_n)
                    .map_err(|e| e.to_string())?,
            );
            models.push((voice.model, voice.initialization));
        }
        let (forces, events) =
            ProjectedForceProgram::compile(inputs, parsed.events, parsed.force_config)
                .map_err(|e| e.to_string())?
                .into_parts();
        let mut observers = Vec::with_capacity(count);
        let mut reports = Vec::with_capacity(count);
        for (((mut model, init), force), (bundle, controls)) in models
            .into_iter()
            .zip(&forces)
            .zip(bundles.iter().zip(settings))
        {
            if init == ForceInitialization::StaticPreload {
                model
                    .initialize_static_equilibrium(force)
                    .map_err(|e| e.to_string())?;
            }
            let (observer, report) =
                BroadbandModalObserver::build(&model, bundle, controls, parsed.info.sample_rate_hz)
                    .map_err(|e| e.to_string())?;
            observers.push(observer);
            reports.push(report);
        }
        let mut info = parsed.info;
        info.input_hash =
            fs_blake3::hash_domain("org.frankensim.modal-performance.broadband.v1", bytes);
        Ok(Self {
            info,
            observers,
            reports,
            forces,
            events,
            next_event: 0,
            completed: 0,
            max_block,
            poisoned: false,
        })
    }
    pub const fn info(&self) -> ModalPerformanceInfo {
        self.info
    }
    pub fn reports(&self) -> &[BroadbandModalReport] {
        &self.reports
    }
    pub fn compiled_controls(&self) -> usize {
        self.events.len()
    }
}
fn read_frequency(
    reader: &mut Reader<'_>,
    kind: &str,
    modes: usize,
) -> Result<(f64, Vec<C64>), String> {
    let mut row = reader.row(kind).map_err(|e| e.to_string())?;
    let omega = row.scalar().map_err(|e| e.to_string())?;
    let mut values = Vec::with_capacity(modes);
    for _ in 0..modes {
        values.push(C64::new(
            row.scalar().map_err(|e| e.to_string())?,
            row.scalar().map_err(|e| e.to_string())?,
        ));
    }
    row.finish().map_err(|e| e.to_string())?;
    Ok((omega, values))
}
fn refused() -> RenderError {
    RenderError::Control {
        what: "broadband receiver failed its admitted model/operator/window budget",
    }
}
impl PressureRenderer for BroadbandModalPerformance {
    fn samples_rendered(&self) -> u64 {
        self.completed
    }
    fn max_block_len(&self) -> usize {
        self.max_block
    }
    fn validate_sample_rate(&self, rate: u32) -> Result<(), RenderError> {
        if self.poisoned || rate != self.info.sample_rate_hz {
            return Err(refused());
        }
        Ok(())
    }
    fn validate_sample_count(&self, samples: u64) -> Result<(), RenderError> {
        if self.poisoned || samples > self.info.samples - self.completed {
            return Err(refused());
        }
        Ok(())
    }
    fn block(&mut self, output: &mut [f64]) -> Result<(), RenderError> {
        self.validate_sample_count(output.len() as u64)?;
        if output.is_empty() || output.len() > self.max_block {
            return Err(refused());
        }
        for slot in output {
            while self.next_event < self.events.len()
                && self.events[self.next_event].sample == self.completed
            {
                let event = self.events[self.next_event];
                self.forces[event.voice][event.mode] = event.force_n_per_sqrt_kg;
                self.next_event += 1;
            }
            let mut pressure = 0.0;
            for (observer, force) in self.observers.iter_mut().zip(&self.forces) {
                match observer.step(force) {
                    Ok(p) => pressure += p,
                    Err(_) => {
                        self.poisoned = true;
                        return Err(refused());
                    }
                }
            }
            if !pressure.is_finite() {
                self.poisoned = true;
                return Err(refused());
            }
            *slot = pressure;
            self.completed += 1;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    const BASE: &str = "frankensim-modal-performance-v1\nsample_rate_hz 48000\nsamples 41\nfull_scale_pa 1\nlimits 0.9 1 100000 1000000000 10000000\ncompile_limits 4 32\nvoices 1\nvoice retain-state 1 1\nmode 3141.592653589793 0.01 0 0 0 0\nport 0 1\nevents 2\nforce 0 0 0 1\nforce 13 0 0 0\n";
    fn constant_input() -> String {
        format!(
            "{BASE}broadband 1\ntransfer test:computed-constant 0 0.03 1 4 1 1e-6 2 1\nfrequency 100 2 -0.3\nfrequency 300 2 -0.3\nwithheld 200 2 -0.3\n"
        )
    }
    #[test]
    fn explicit_block_is_complete_and_callback_partition_is_bitwise() {
        assert!(!BroadbandModalPerformance::has_section(BASE.as_bytes()));
        let text = constant_input();
        assert!(BroadbandModalPerformance::has_section(text.as_bytes()));
        let mut a = BroadbandModalPerformance::from_bytes(text.as_bytes(), 41).unwrap();
        let mut b = BroadbandModalPerformance::from_bytes(text.as_bytes(), 17).unwrap();
        assert_eq!(a.reports()[0].producer_relative_error_estimate, 0.03);
        let mut whole = [0.0; 41];
        a.block(&mut whole).unwrap();
        let mut pieces = [0.0; 41];
        for block in pieces.chunks_mut(17) {
            b.block(block).unwrap();
        }
        for (x, y) in whole.into_iter().zip(pieces) {
            assert_eq!(x.to_bits(), y.to_bits());
        }
        assert!(b.block(&mut [0.0]).is_err());
    }
    #[test]
    fn malformed_or_missing_operator_never_silently_uses_legacy() {
        assert!(BroadbandModalPerformance::from_bytes(BASE.as_bytes(), 17).is_err());
        let good = constant_input();
        for bad in [
            good.replace("broadband 1", "broadband 2"),
            good.replace("withheld 200 2 -0.3", "withheld 100 2 -0.3"),
            good.replace(" 2 1\nfrequency", " 513 1\nfrequency"),
            good.replace(
                "transfer test:computed-constant 0 0.03",
                "transfer test:computed-constant 0 -1",
            ),
            format!("{good}trailing ignored-data\n"),
        ] {
            assert!(BroadbandModalPerformance::from_bytes(bad.as_bytes(), 17).is_err());
        }
    }
}
