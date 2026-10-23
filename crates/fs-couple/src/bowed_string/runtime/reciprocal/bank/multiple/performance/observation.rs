//! Finite-distance observation of the original closed material shell.
//!
//! This is an output owner only: the accepted source-bound air reaction remains
//! inside the existing mechanics. Boundary traces, gas and geometry policy are
//! shared with that load. The existing observation producer supplies centroid
//! receiver p/grad rows and the existing FIR retains their causal history.
//! Finite-band reconstruction convergence does not certify spatial accuracy or
//! molecular gas attenuation.

use super::{MultiBowPerformance, render_invalid};
use crate::bowed_string::runtime::reciprocal::bank::radiation::bem::ShellBemError;
use crate::render::RenderError;
use crate::render::schedule::force::coupled::contact::friction::multiple::MultiFrictionModalSystem;
use fs_bem::helmholtz::{Medium, PanelRule};
use fs_bem::observation_transfers::{
    BemObservation, ObservationTransferBake, ProducerUncertainty, ReceiverSpec,
    SourceId, TransferDeclaration, TransferError, MAX_OBSERVATION_RUNGS,
};
use fs_blake3::ContentHash;
use fs_exec::CancelGate;
use fs_material::gas::GasState;
use fs_vfit::observation::{MAX_OBSERVATION_ENTRIES, ObservationFir};

/// Complete receiver and realization declaration. No microphone distance,
/// reconstruction delay, tolerance or work allowance is inferred from a name.
#[derive(Clone, Debug, PartialEq)]
pub struct ShellObserverOptions {
    pub receiver: ReceiverSpec,
    pub minimum_clearance_m: f64,
    pub pass_band_upper_hz: f64,
    pub taper_stop_hz: f64,
    pub filter_latency_samples: usize,
    pub fft_lengths: Vec<usize>,
    pub tail_energy_tolerance: f64,
    pub dtft_relative_tolerance: f64,
    pub maximum_frequency_entries: usize,
    /// Conservative frequency-count * dense shell algebraic work bound, using
    /// the same panel/mode expression as the original loading owner.
    pub maximum_sweep_terms: u64,
    /// Maximum retained FIR multiply-add terms PER accepted physical sample,
    /// including all four pressure/gradient outputs.
    pub maximum_replay_terms: u64,
}

#[derive(Debug)]
pub enum ShellObserverError {
    Invalid(&'static str),
    Shell(ShellBemError),
    Transfer(TransferError),
    Cancelled,
}
impl core::fmt::Display for ShellObserverError {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(formatter, "closed-shell observer: {self:?}")
    }
}
impl std::error::Error for ShellObserverError {}
impl From<TransferError> for ShellObserverError {
    fn from(error: TransferError) -> Self { Self::Transfer(error) }
}
fn poll(gate: &CancelGate) -> Result<(), ShellObserverError> {
    if gate.is_requested() { Err(ShellObserverError::Cancelled) } else { Ok(()) }
}
fn invalid(what: &'static str) -> ShellObserverError { ShellObserverError::Invalid(what) }

impl ShellObserverOptions {
    /// Reuse the common finite-band declaration admission before shell or pack
    /// work. The temporary one-source shape admits numerical options only;
    /// the full original basis and physical identities are admitted at binding.
    pub fn validate(&self, rate: u32) -> Result<(), ShellObserverError> {
        if rate == 0 || self.receiver.id.is_empty() || self.receiver.id.len() > 1024
            || self.receiver.position_m.iter().any(|x| !x.is_finite())
            || self.fft_lengths.len() > MAX_OBSERVATION_RUNGS
            || !self.minimum_clearance_m.is_finite() || self.minimum_clearance_m <= 0.0
            || !(1..=100_000_000_000).contains(&self.maximum_sweep_terms)
            || self.maximum_replay_terms == 0
            || self.maximum_replay_terms > MAX_OBSERVATION_ENTRIES as u64
        {
            return Err(invalid("invalid receiver clearance or observation work allowance"));
        }
        let declaration = self.declaration(
            rate, vec![SourceId::StructuralMode("numerical admission only".into())],
            vec![0.0], ContentHash([1; 32]), ContentHash([1; 32]),
        );
        declaration.admit_work()?;
        let finest = self.fft_lengths.last().copied()
            .ok_or_else(|| invalid("missing observer FFT ladder"))?;
        if self.pass_band_upper_hz < f64::from(rate) / (finest as f64) {
            return Err(invalid("observer pass band must contain a resolved positive FFT bin"));
        }
        Ok(())
    }

    fn declaration(
        &self, rate: u32, sources: Vec<SourceId>, geometric_delays_s: Vec<f64>,
        producer_hash: ContentHash, loading_identity: ContentHash,
    ) -> TransferDeclaration {
        TransferDeclaration {
            producer_hash, sources, receivers: vec![self.receiver.clone()],
            fs_hz: f64::from(rate), geometric_delays_s,
            pass_band_upper_hz: self.pass_band_upper_hz,
            taper_stop_hz: self.taper_stop_hz,
            filter_latency_samples: self.filter_latency_samples,
            fft_lengths: self.fft_lengths.clone(),
            tail_energy_tolerance: self.tail_energy_tolerance,
            dtft_relative_tolerance: self.dtft_relative_tolerance,
            producer_uncertainty: ProducerUncertainty::Unresolved { source_hash: loading_identity },
            maximum_frequency_entries: self.maximum_frequency_entries,
        }
    }
}

/// Retained finite-band producer data and the exact original load/source binding.
/// This reports realization estimates, not continuum or measured-source proof.
#[derive(Clone, Debug)]
pub struct ShellObserverReport {
    pub options: ShellObserverOptions,
    pub basis_identity: ContentHash,
    pub loading_identity: ContentHash,
    pub gas: GasState,
    pub body_modes: usize,
    pub surface_panels: usize,
    pub sweep_terms: u64,
    /// Complete retained convolution terms per accepted sample.
    pub replay_terms: u64,
    transfer: ObservationTransferBake,
}
impl ShellObserverReport {
    pub const fn transfer(&self) -> &ObservationTransferBake { &self.transfer }
}

pub(super) struct BoundShellObserver {
    report: ShellObserverReport,
    observer: ObservationFir,
    body_component: usize,
    velocity: Vec<f64>,
    maximum_abs_pressure_pa: f64,
    samples: u64,
}
impl BoundShellObserver {
    /// Called once, after the original mechanics accepts the complete sample.
    /// A failure poisons the outer performance; it never pretends to roll back
    /// accepted mechanics or reconstruct an unknown observation prehistory.
    pub(super) fn pressure(
        &mut self, system: &MultiFrictionModalSystem,
    ) -> Result<f64, RenderError> {
        if self.samples.checked_add(1) != Some(system.samples_rendered()) {
            return Err(render_invalid("closed-shell observer lost its accepted mechanical clock"));
        }
        let body = system.components().get(self.body_component)
            .ok_or_else(|| render_invalid("closed-shell observer lost its original body"))?;
        if body.states().len() != self.velocity.len() {
            return Err(render_invalid("closed-shell observer lost its complete original modal basis"));
        }
        for (velocity, state) in self.velocity.iter_mut().zip(body.states()) {
            *velocity = state.velocity_m_sqrt_kg_per_s;
        }
        let output = self.observer.step(&self.velocity)
            .map_err(|_| render_invalid("closed-shell observation failed on accepted modal velocities"))?;
        let pressure = output[0];
        if !pressure.is_finite() || pressure.abs() > self.maximum_abs_pressure_pa {
            return Err(render_invalid("closed-shell receiver exceeded the original pressure allowance"));
        }
        self.samples = system.samples_rendered();
        Ok(pressure)
    }
}

impl MultiBowPerformance {
    /// Replace compact pressure with one finite-distance closed-shell receiver.
    /// The source-bound BEM load, accepted mechanics, force schedule and static
    /// preload are retained. The FIR begins with explicit zero prehistory, so
    /// a bank that has accepted any physical sample is refused.
    ///
    /// A raw supplied radiation model has no geometry-route evidence and cannot
    /// substitute for the original material-shell BEM binding. Cavity/open-body
    /// coupling and molecular gas attenuation are outside this observation path.
    pub fn with_material_shell_observer(
        mut self, options: ShellObserverOptions, gate: &CancelGate,
    ) -> Result<Self, ShellObserverError> {
        poll(gate)?;
        if self.poisoned || self.shell_observer.is_some() || self.base_sample != 0
            || self.bank.system.samples_rendered() != 0
        {
            return Err(invalid("shell observation requires one unstarted physical history"));
        }
        let source = self.bank.source.bank();
        let rate = source.options().sample_rate_hz;
        options.validate(rate)?;
        if source.material_cavity().is_some() || source.prepared_linear_budget().is_some() {
            return Err(invalid("closed-shell observation requires the exact loaded owner without a cavity"));
        }
        let binding = source.body_radiation()
            .ok_or_else(|| invalid("shell observation requires the matching original BEM air load"))?;
        let loading = binding.material_shell_bem()
            .ok_or_else(|| invalid("raw radiation models have no material-shell observation binding"))?;
        if binding.basis_identity() != loading.basis_identity || binding.gas() != source.gas()
            || &loading.gas != source.gas()
        {
            return Err(invalid("shell observation no longer matches its original body basis and gas"));
        }
        source.gas().validate_coherence()
            .map_err(|_| invalid("incoherent original shell observation gas"))?;
        let finest = options.fft_lengths.last().copied()
            .ok_or_else(|| invalid("missing observer FFT ladder"))?;
        if options.taper_stop_hz > loading.options.maximum_hz
            || f64::from(rate) / (finest as f64) < loading.options.minimum_hz
        {
            return Err(invalid("every positive observation frequency must lie inside the original load band"));
        }
        let body_component = source.strings().len();
        let body_modes = self.bank.system.components().get(body_component)
            .ok_or_else(|| invalid("missing original shell body component"))?.modes().len();
        if body_modes != loading.body_modes {
            return Err(invalid("shell observation must retain every original body mode in order"));
        }
        let replay_terms = (body_modes as u64).checked_mul(4)
            .and_then(|n| n.checked_mul((finest / 2) as u64))
            .ok_or_else(|| invalid("complete observation replay extent overflow"))?;
        if replay_terms > options.maximum_replay_terms {
            return Err(invalid("complete shell observer exceeds its per-sample replay allowance"));
        }
        // Nested FFT grids share a cache. This bounds every positive frequency
        // below the stop band, conservatively including its last boundary bin.
        let frequency_samples =
            (options.taper_stop_hz * finest as f64 / f64::from(rate)).ceil() as usize;
        let (surface, traces, sweep_terms) = source.shell_surface_traces(
            body_modes, loading.options.maximum_surface_panels, options.maximum_sweep_terms,
            frequency_samples, gate,
        ).map_err(ShellObserverError::Shell)?;
        let surface_panels = surface.areas().len();
        if surface_panels != loading.surface_panels {
            return Err(invalid("shell observation geometry differs from the original radiation load"));
        }
        let sources = (0..body_modes).map(|index| SourceId::StructuralMode(
            format!("material-shell-mode:{index}:{}", loading.basis_identity),
        )).collect::<Vec<_>>();
        let gas = *source.gas();
        let owner = BemObservation::new_routed(
            &surface, Medium { density: gas.density, sound_speed: gas.sound_speed },
            &sources, &traces, core::slice::from_ref(&options.receiver), options.minimum_clearance_m,
            PanelRule::Triangle { max_pairs: surface_panels.saturating_pow(2) },
        )?;
        let declaration = options.declaration(
            rate, sources, owner.geometric_delays_s(), owner.operator_hash(), loading.loading_identity,
        );
        declaration.admit_work()?;
        let transfer = owner.bake(&declaration, || gate.is_requested())?;
        let observer = transfer.observer()?;
        poll(gate)?;
        let report = ShellObserverReport {
            options, basis_identity: loading.basis_identity, loading_identity: loading.loading_identity,
            gas, body_modes, surface_panels, sweep_terms, replay_terms, transfer,
        };
        self.shell_observer = Some(BoundShellObserver {
            report, observer, body_component, velocity: vec![0.0; body_modes],
            maximum_abs_pressure_pa: source.options().network.maximum_abs_pressure_pa,
            samples: 0,
        });
        Ok(self)
    }

    /// Original bound receiver declaration and completed finite-band realization.
    pub fn material_shell_observer(&self) -> Option<&ShellObserverReport> {
        self.shell_observer.as_ref().map(|observer| &observer.report)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn options() -> ShellObserverOptions {
        ShellObserverOptions {
            receiver: ReceiverSpec { id: "front".into(), position_m: [0.2, 0.2, 1.0] },
            minimum_clearance_m: 0.01, pass_band_upper_hz: 120.0, taper_stop_hz: 400.0,
            filter_latency_samples: 256, fft_lengths: vec![1024, 2048, 4096, 8192],
            tail_energy_tolerance: 1e-10, dtft_relative_tolerance: 1e-4,
            maximum_frequency_entries: 65_536, maximum_sweep_terms: 10_000_000_000,
            maximum_replay_terms: 65_536,
        }
    }

    #[test]
    fn g0_receiver_declaration_keeps_common_fft_tolerances_and_explicit_work_bounds() {
        options().validate(8000).unwrap();
        let mut changed = options();
        changed.maximum_replay_terms = 0;
        assert!(changed.validate(8000).is_err());
        changed = options();
        changed.maximum_sweep_terms = 0;
        assert!(changed.validate(8000).is_err());
        changed = options();
        changed.fft_lengths = vec![1024, 2048, 4096];
        assert!(changed.validate(8000).is_err());
        changed = options();
        changed.pass_band_upper_hz = 0.1;
        assert!(changed.validate(8000).is_err(), "no unresolved all-zero frequency grid");
        changed = options();
        changed.tail_energy_tolerance = 1e-3;
        assert!(changed.validate(8000).is_err());
        changed = options();
        changed.receiver.position_m[0] = f64::NAN;
        assert!(changed.validate(8000).is_err());
        changed = options();
        changed.minimum_clearance_m = 0.0;
        assert!(changed.validate(8000).is_err());
    }
}
