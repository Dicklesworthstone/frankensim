//! Causal pressure observation of the complete coupled source trajectory.
use super::{CoupledModalSystem, CoupledModalTones, CoupledToneBudget, ModalCouplingError, poll};
pub use crate::broadband_radiation::modal::{
    BroadbandModalControls, BroadbandModalError, BroadbandModalReport, ModalTransferSamples,
};
use crate::broadband_radiation::modal::{
    MAX_MODAL_DELAY_KERNEL_HALF_LEN, MAX_MODAL_TRANSFER_ORDER, admit_modal_samples,
};
use crate::modal_acoustic_time::broadband::{
    Delay, MAX_MODAL_OBSERVER_DELAY_SAMPLES, prepare_modal_transfer_observation,
};
use fs_exec::CancelGate;
use fs_math::c64::C64;
use fs_phs::linear_exact::PreparedToneObserver;
use fs_vfit::{PoleTerm, RationalModel};

/// Original source or receiver realization/observation refusal.
#[derive(Debug)]
pub enum BroadbandCoupledModalError {
    Modal(ModalCouplingError),
    Observation(BroadbandModalError),
}
impl From<ModalCouplingError> for BroadbandCoupledModalError {
    fn from(error: ModalCouplingError) -> Self {
        Self::Modal(error)
    }
}
impl From<BroadbandModalError> for BroadbandCoupledModalError {
    fn from(error: BroadbandModalError) -> Self {
        Self::Observation(error)
    }
}
impl core::fmt::Display for BroadbandCoupledModalError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Modal(error) => error.fmt(f),
            Self::Observation(error) => error.fmt(f),
        }
    }
}
impl std::error::Error for BroadbandCoupledModalError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        Some(match self {
            Self::Modal(error) => error,
            Self::Observation(error) => error,
        })
    }
}

/// One fitted receiver driven by the ORIGINAL continuous coupled velocities.
///
/// Every signed transfer row addresses the same flattened mechanical mode as
/// its generalized force. Maxwell and reciprocal-air states remain in the
/// source operator and change those velocities; they are not extra radiators.
/// Continuous rational filters share an exact source/filter exponential. The
/// source's actual work and heat remain entirely with its original owner.
///
/// Build at the source's initial clock. Filter and delay histories start at
/// zero. Existing component pressure maps are explicitly cleared on accepted
/// and candidate models; this changes no physical coefficient, kernel identity,
/// state, source force, heat, radiation load or mechanical limit. This receiver
/// owns the network's original pressure ceiling. Its pressure is returned
/// separately from the original mechanical frame.
pub struct BroadbandCoupledModalTones {
    source: CoupledModalTones,
    observation: PreparedToneObserver,
    models: Vec<RationalModel>,
    state: Vec<f64>,
    candidate: Vec<f64>,
    output: [f64; 1],
    delay: Delay,
    report: BroadbandModalReport,
    last_pressure: Option<f64>,
}

impl BroadbandCoupledModalTones {
    /// Fit the solver-neutral computed transfers and admit the actual delayed
    /// realization against the same independent held-out samples as the
    /// existing modal receiver. The declared sample clock must match bitwise.
    ///
    /// The additional receiver byte allowance covers conservative fit QR/eigen
    /// workspace, retained rational models, dense A/B/C/D assembly, filter
    /// histories and the delay ring. Remaining bytes and the explicit dense
    /// multiplication allowance go to the shared exact observer preparation.
    /// Fitting keeps its own original order/iteration controls. No mode, pole,
    /// tone or physical source state is removed to fit these allowances.
    pub fn build(
        mut source: CoupledModalTones,
        samples: &ModalTransferSamples,
        controls: BroadbandModalControls,
        sample_rate_hz: u32,
        observer_budget: CoupledToneBudget,
        gate: &CancelGate,
    ) -> Result<(Self, BroadbandModalReport), BroadbandCoupledModalError> {
        poll(Some(gate))?;
        let network = source.network();
        let modes = network.mode_count();
        if sample_rate_hz == 0
            || network.samples_rendered() != 0
            || network.sample_period_s().to_bits() != (1.0 / f64::from(sample_rate_hz)).to_bits()
            || samples.per_mode.len() != modes
        {
            return Err(BroadbandModalError::Invalid(
                "coupled receiver requires one cold source, its exact clock and complete original modal basis",
            ).into());
        }
        let reserve = receiver_reservation(samples, controls, sample_rate_hz, observer_budget)?;
        let (bank, delay) = prepare_modal_transfer_observation(samples, controls, sample_rate_hz)?;
        poll(Some(gate))?;
        let models = if let Some(constant) = bank.constant {
            if constant.iter().any(|h| h.im != 0.0) {
                return Err(BroadbandModalError::Invalid(
                    "a complex constant quadrature is not a causal real broadband transfer",
                )
                .into());
            }
            constant
                .into_iter()
                .map(|h| RationalModel {
                    terms: Vec::new(),
                    d: h.re,
                    e: 0.0,
                })
                .collect()
        } else {
            bank.models
        };
        let states = models
            .iter()
            .try_fold(0usize, |n, model| n.checked_add(model.order()))
            .ok_or(BroadbandModalError::Invalid(
                "coupled receiver filter extent overflow",
            ))?;
        let physical =
            2 * modes + network.material_memory().len() + network.radiation_state().len();
        if states > fs_phs::linear_exact::MAX_LINEAR_EXPONENTIAL_BLOCK_DIM
            || states > observer_budget.max_dimension
        {
            return Err(BroadbandModalError::Cap {
                what: "complete coupled receiver filter states",
                found: states,
                cap: observer_budget
                    .max_dimension
                    .min(fs_phs::linear_exact::MAX_LINEAR_EXPONENTIAL_BLOCK_DIM),
            }
            .into());
        }
        let matrices = states
            .checked_mul(states)
            .and_then(|v| v.checked_add(states.checked_mul(physical)?))
            .and_then(|v| v.checked_add(3 * states + physical + 1))
            .and_then(|v| v.checked_mul(core::mem::size_of::<f64>()))
            .and_then(|v| v.checked_add(reserve))
            .ok_or(BroadbandModalError::Invalid(
                "coupled receiver matrix byte overflow",
            ))?;
        admit_bytes(matrices, observer_budget.max_dense_bytes)?;
        let mut a = vec![0.0; states * states];
        let mut b = vec![0.0; states * physical];
        let mut c = vec![0.0; states];
        let mut d = vec![0.0; physical];
        let mut offset = 0;
        for (mode, model) in models.iter().enumerate() {
            poll(Some(gate))?;
            let ss = model.state_space();
            for i in 0..ss.n {
                a[(offset + i) * states + offset..(offset + i) * states + offset + ss.n]
                    .copy_from_slice(&ss.a[i * ss.n..(i + 1) * ss.n]);
                b[(offset + i) * physical + 2 * mode + 1] = ss.b[i];
                c[offset + i] = ss.c[i];
            }
            d[2 * mode + 1] = ss.d;
            offset += ss.n;
        }
        let observation = source.prepared.prepare_observer(
            states,
            &a,
            &b,
            &c,
            &d,
            CoupledToneBudget {
                max_dense_bytes: observer_budget.max_dense_bytes - matrices,
                ..observer_budget
            },
            gate,
        )?;
        poll(Some(gate))?;
        for model in source
            .network
            .models
            .iter_mut()
            .chain(&mut source.network.candidates)
        {
            *model = model.without_pressure_transfers();
        }
        let mut report = bank.report;
        // The older independent-mode operation-order shortcut is not used by
        // this complete coupled source/filter exponential, even for real D.
        report.constant_narrow_band_compatibility = false;
        let result = Self {
            source,
            observation,
            models,
            state: vec![0.0; states],
            candidate: vec![0.0; states],
            output: [0.0],
            delay,
            report: report.clone(),
            last_pressure: None,
        };
        Ok((result, report))
    }

    #[must_use]
    pub const fn source(&self) -> &CoupledModalTones {
        &self.source
    }
    #[must_use]
    pub fn network(&self) -> &CoupledModalSystem {
        self.source.network()
    }
    /// Recover the accepted source with its observation-only pressure maps
    /// cleared. Mechanical coefficients, histories and ledgers are unchanged.
    #[must_use]
    pub fn into_source(self) -> CoupledModalTones {
        self.source
    }
    #[must_use]
    pub const fn report(&self) -> &BroadbandModalReport {
        &self.report
    }
    #[must_use]
    pub const fn last_pressure_pa(&self) -> Option<f64> {
        self.last_pressure
    }

    pub fn advance_force_phase(
        &self,
        tone: usize,
        cosine: f64,
        sine: f64,
    ) -> Result<(f64, f64), ModalCouplingError> {
        self.source.advance_force_phase(tone, cosine, sine)
    }

    /// Fitted modal-velocity/pressure operator, including the actual sampled
    /// delay, in exp(-i*omega*t). Its sampled holdout is an estimate, not a
    /// continuum error bound or a force-to-pressure transfer.
    pub fn transfer_at(&self, mode: usize, omega: f64) -> Result<C64, BroadbandCoupledModalError> {
        if mode >= self.models.len()
            || !omega.is_finite()
            || omega <= 0.0
            || omega * self.network().sample_period_s() >= core::f64::consts::PI
        {
            return Err(BroadbandModalError::Invalid(
                "operator query requires an original mode and positive sub-Nyquist frequency",
            )
            .into());
        }
        Ok(self.models[mode].eval_iw(omega).conj()
            * self.delay.response(omega, self.network().sample_period_s()))
    }

    /// Advance the actual continuous source and return accepted receiver Pa.
    /// Source, filter, delay, clock and last pressure remain unchanged on error.
    pub fn step(
        &mut self,
        cosine: &[f64],
        sine: &[f64],
    ) -> Result<f64, BroadbandCoupledModalError> {
        self.step_inner(cosine, sine, None)
    }
    pub fn step_under_gate(
        &mut self,
        cosine: &[f64],
        sine: &[f64],
        gate: &CancelGate,
    ) -> Result<f64, BroadbandCoupledModalError> {
        self.step_inner(cosine, sine, Some(gate))
    }
    fn step_inner(
        &mut self,
        cosine: &[f64],
        sine: &[f64],
        gate: Option<&CancelGate>,
    ) -> Result<f64, BroadbandCoupledModalError> {
        self.source.stage(cosine, sine, gate)?;
        self.finish_staged(gate)
    }

    /// Observe an affine physical tone envelope along the same complete
    /// continuous trajectory. Both endpoint pairs use the interval-start
    /// carrier basis and the source must have been prepared with `new_affine`.
    /// A release retains source, filter and delay histories for ringdown.
    pub fn step_affine(
        &mut self,
        cosine_start: &[f64],
        sine_start: &[f64],
        cosine_end: &[f64],
        sine_end: &[f64],
    ) -> Result<f64, BroadbandCoupledModalError> {
        self.step_affine_inner(cosine_start, sine_start, cosine_end, sine_end, None)
    }

    /// Cancellable affine source/filter/delay transaction.
    pub fn step_affine_under_gate(
        &mut self,
        cosine_start: &[f64],
        sine_start: &[f64],
        cosine_end: &[f64],
        sine_end: &[f64],
        gate: &CancelGate,
    ) -> Result<f64, BroadbandCoupledModalError> {
        self.step_affine_inner(cosine_start, sine_start, cosine_end, sine_end, Some(gate))
    }

    fn step_affine_inner(
        &mut self,
        cosine_start: &[f64],
        sine_start: &[f64],
        cosine_end: &[f64],
        sine_end: &[f64],
        gate: Option<&CancelGate>,
    ) -> Result<f64, BroadbandCoupledModalError> {
        self.source
            .stage_affine(cosine_start, sine_start, cosine_end, sine_end, gate)?;
        self.finish_staged(gate)
    }

    fn finish_staged(
        &mut self,
        gate: Option<&CancelGate>,
    ) -> Result<f64, BroadbandCoupledModalError> {
        self.source.prepared.preview_observer(
            &mut self.observation,
            &self.state,
            &mut self.candidate,
            &mut self.output,
            gate,
        )?;
        let raw = self.output[0];
        let pressure = self.delay.preview(raw)?;
        let maximum = self.source.network.config.maximum_abs_pressure_pa;
        if pressure.abs() > maximum {
            return Err(ModalCouplingError::Budget {
                what: "broadband receiver pressure",
                value: pressure.abs(),
                limit: maximum,
            }
            .into());
        }
        poll(gate)?;
        self.source.network.publish_staged();
        core::mem::swap(&mut self.state, &mut self.candidate);
        self.delay.commit(raw);
        self.last_pressure = Some(pressure);
        Ok(pressure)
    }
}

fn admit_bytes(found: usize, cap: usize) -> Result<(), BroadbandModalError> {
    if found >= cap {
        return Err(BroadbandModalError::Cap {
            what: "additional coupled receiver bytes",
            found,
            cap,
        });
    }
    Ok(())
}

/// Conservative additional peak storage, admitted before fitting or delay
/// allocation. The original fitter controls bound its numerical iterations;
/// its QR stores two tall matrices and the relocation eigenproblem is small.
fn receiver_reservation(
    samples: &ModalTransferSamples,
    controls: BroadbandModalControls,
    rate: u32,
    budget: CoupledToneBudget,
) -> Result<usize, BroadbandModalError> {
    admit_modal_samples(
        samples.per_mode.len(),
        samples.omega_rad_s.len(),
        samples.held_out.len(),
    )?;
    if controls.fit_order == 0
        || controls.fit_order > MAX_MODAL_TRANSFER_ORDER
        || controls.fit_iterations == 0
        || controls.fit_iterations > 128
        || controls.delay_kernel_half_len == 0
        || controls.delay_kernel_half_len > MAX_MODAL_DELAY_KERNEL_HALF_LEN
    {
        return Err(BroadbandModalError::Invalid(
            "bounded explicit fit and sampled-delay controls required",
        ));
    }
    let delay_samples = samples.propagation_delay_s * f64::from(rate);
    if !delay_samples.is_finite()
        || delay_samples < 0.0
        || delay_samples > MAX_MODAL_OBSERVER_DELAY_SAMPLES as f64
    {
        return Err(BroadbandModalError::Invalid(
            "receiver delay exceeds its existing history capacity",
        ));
    }
    let order = controls.fit_order;
    let rows = 2 * samples.omega_rad_s.len() + 1;
    let columns = 2 * order + 2;
    let half = controls.delay_kernel_half_len;
    let fit_entries = 16 * rows * columns + 32 * order * order + 64 * (rows + order);
    let fit_bytes = fit_entries.checked_mul(core::mem::size_of::<f64>()).ok_or(
        BroadbandModalError::Invalid("receiver fit workspace byte overflow"),
    )?;
    let retained_models = samples
        .per_mode
        .len()
        .checked_mul(
            order * core::mem::size_of::<PoleTerm>() + core::mem::size_of::<RationalModel>(),
        )
        .ok_or(BroadbandModalError::Invalid(
            "receiver rational model byte overflow",
        ))?;
    let delay_bytes = (delay_samples.ceil() as usize + half + 1)
        .checked_mul(core::mem::size_of::<f64>())
        // The shared tap owner grows a Vec by pushing; reserve its possible
        // spare capacity as well as the logically retained taps.
        .and_then(|v| v.checked_add((2 * (2 * half + 1) + 4) * core::mem::size_of::<(usize, f64)>()))
        .ok_or(BroadbandModalError::Invalid("receiver delay byte overflow"))?;
    let pressure_map_copy = 2
        * samples.per_mode.len()
        * (core::mem::size_of::<crate::modal_acoustic_time::ModalAcousticMode>()
            + core::mem::size_of::<crate::modal_acoustic_time::ModalAcousticState>());
    let bytes = fit_bytes
        .checked_add(retained_models)
        .and_then(|v| v.checked_add(delay_bytes))
        .and_then(|v| v.checked_add(pressure_map_copy))
        .and_then(|v| v.checked_add(samples.source_id.len().checked_mul(2)?))
        .and_then(|v| v.checked_add(core::mem::size_of::<BroadbandCoupledModalTones>()))
        .ok_or(BroadbandModalError::Invalid(
            "receiver retained byte overflow",
        ))?;
    admit_bytes(bytes, budget.max_dense_bytes)?;
    Ok(bytes)
}

#[cfg(test)]
mod tests;
