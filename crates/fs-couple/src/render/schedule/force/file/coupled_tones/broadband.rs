//! A complete coupled source and one explicitly sampled causal receiver.
use super::super::{Reader, broadband::read_transfer};
use super::{MAX_PREPARATION_WORK, MAX_STORAGE_BYTES};
use crate::modal_acoustic_time::broadband::{
    BroadbandModalControls, BroadbandModalReport, ModalTransferSamples,
};
use crate::render::{
    RenderError,
    schedule::force::coupled::{
        CoupledModalSystem, ModalCouplingError,
        tones::{
            CoupledModalTones, CoupledToneBudget,
            broadband::{BroadbandCoupledModalError, BroadbandCoupledModalTones},
        },
    },
};
use fs_exec::CancelGate;

pub(super) struct ReceiverSpec {
    samples: ModalTransferSamples,
    controls: BroadbandModalControls,
    budget: CoupledToneBudget,
}

impl ReceiverSpec {
    /// Read the optional suffix after the original physical tone records. Every
    /// frequency row contains the complete flattened source basis, in the same
    /// voice-major order as the original mechanical file.
    pub(super) fn read_optional(
        reader: &mut Reader<'_>,
        modes: usize,
    ) -> Result<Option<Self>, String> {
        let Some(next) = reader.lines.clone().next() else {
            return Ok(None);
        };
        if next.split_ascii_whitespace().next() != Some("coupled_broadband") {
            return Err("unexpected trailing coupled-tone record".into());
        }
        let mut row = reader.row("coupled_broadband").map_err(|e| e.to_string())?;
        let budget = CoupledToneBudget {
            max_dimension: row
                .count(fs_phs::linear_exact::MAX_LINEAR_EXPONENTIAL_BLOCK_DIM)
                .map_err(|e| e.to_string())?,
            max_dense_bytes: row.count(MAX_STORAGE_BYTES).map_err(|e| e.to_string())?,
            max_multiply_terms: row.parse().map_err(|e| e.to_string())?,
        };
        row.finish().map_err(|e| e.to_string())?;
        if budget.max_dimension == 0
            || budget.max_dense_bytes == 0
            || budget.max_multiply_terms == 0
            || budget.max_multiply_terms > MAX_PREPARATION_WORK
        {
            return Err(
                "coupled receiver needs positive bounded dimension, storage and work".into(),
            );
        }
        let (samples, controls) = read_transfer(reader, modes)?;
        if reader.lines.next().is_some() {
            return Err("unexpected trailing coupled receiver record".into());
        }
        Ok(Some(Self {
            samples,
            controls,
            budget,
        }))
    }
}

pub(super) enum Source {
    Instantaneous(CoupledModalTones),
    Broadband(Box<BroadbandCoupledModalTones>),
}

impl Source {
    pub(super) fn new(
        source: CoupledModalTones,
        receiver: Option<ReceiverSpec>,
        sample_rate_hz: u32,
        gate: &CancelGate,
    ) -> Result<Self, String> {
        match receiver {
            None => Ok(Self::Instantaneous(source)),
            Some(spec) => {
                let (receiver, _) = BroadbandCoupledModalTones::build(
                    source,
                    &spec.samples,
                    spec.controls,
                    sample_rate_hz,
                    spec.budget,
                    gate,
                )
                .map_err(|e| e.to_string())?;
                Ok(Self::Broadband(Box::new(receiver)))
            }
        }
    }

    fn physical(&self) -> &CoupledModalTones {
        match self {
            Self::Instantaneous(source) => source,
            Self::Broadband(receiver) => receiver.source(),
        }
    }

    pub(super) fn network(&self) -> &CoupledModalSystem {
        self.physical().network()
    }

    pub(super) fn tone_count(&self) -> usize {
        self.physical().tone_count()
    }

    pub(super) fn advance_force_phase(
        &self,
        tone: usize,
        cosine: f64,
        sine: f64,
    ) -> Result<(f64, f64), ModalCouplingError> {
        self.physical().advance_force_phase(tone, cosine, sine)
    }

    pub(super) fn report(&self) -> Option<&BroadbandModalReport> {
        match self {
            Self::Instantaneous(_) => None,
            Self::Broadband(receiver) => Some(receiver.report()),
        }
    }

    pub(super) fn step(&mut self, cosine: &[f64], sine: &[f64]) -> Result<f64, RenderError> {
        match self {
            Self::Instantaneous(source) => source
                .step(cosine, sine)
                .map(|frame| frame.observer_pressure_pa)
                .map_err(RenderError::Coupled),
            Self::Broadband(receiver) => receiver.step(cosine, sine).map_err(|error| match error {
                BroadbandCoupledModalError::Modal(error) => RenderError::Coupled(error),
                BroadbandCoupledModalError::Observation(_) => RenderError::Control {
                    what: "causal coupled receiver refused the candidate sample",
                },
            }),
        }
    }
}
