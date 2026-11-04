//! Molecular time-error control in the original complete network transaction.

use super::super::dissipation::AdaptiveNetworkDissipationControl;
use super::*;

/// Acoustic controls and independent physical tolerances for molecular time error.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AdaptiveDissipatingStepControl {
    pub transport: StepControl,
    /// Full/two-half RKL2 probes share one RHS budget across the complete step.
    pub molecular: AdaptiveNetworkDissipationControl,
    pub maximum_coupling_retries: usize,
}

impl Default for AdaptiveDissipatingStepControl {
    fn default() -> Self {
        Self {
            transport: StepControl::default(),
            molecular: AdaptiveNetworkDissipationControl::default(),
            maximum_coupling_retries: 20,
        }
    }
}

/// Accepted physical balance and molecular temporal diagnostics.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AdaptiveDissipatingNetworkReport {
    pub coupled: DissipatingNetworkReport,
    pub dissipation: SplitAdaptiveNetworkDissipationReport,
    /// Includes all estimator probes and all discarded complete compositions.
    pub operator_evaluations: usize,
}

impl<L: TransportLaw> DissipatingNetwork1d<L> {
    /// Advance D/2, A, D/2 with automatically controlled molecular time error.
    ///
    /// The acoustic timestep still determines the complete composition. Each
    /// diffusion half internally selects its macrosteps using the original
    /// source-backed transport and its finite-plenum and port-transfer errors.
    /// The reported estimator covers molecular integration only; it does not
    /// estimate advection, source integration, spatial or splitting error.
    pub fn step_adaptive(
        &mut self,
        maximum_duration: Time,
        controls: AdaptiveDissipatingStepControl,
        cancelled: impl FnMut() -> bool,
    ) -> Result<AdaptiveDissipatingNetworkReport, DissipatingNetworkError> {
        self.step_adaptive_with_sources(
            maximum_duration,
            controls,
            None,
            NetworkSources::default(),
            cancelled,
        )
    }

    /// Compose error-controlled diffusion with prescribed baths, walls and valves.
    ///
    /// Both molecular halves use the same tolerances and complete-step RHS
    /// allowance. Temporal probes and coupling retries never reset that budget.
    /// Only the accepted two-half diffusion paths contribute to physical port
    /// exchanges. Final state, source storage and reports publish atomically.
    pub fn step_adaptive_with_sources(
        &mut self,
        maximum_duration: Time,
        controls: AdaptiveDissipatingStepControl,
        reservoirs: Option<&[[Option<Cell<'_>>; 2]]>,
        sources: NetworkSources<'_, '_>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AdaptiveDissipatingNetworkReport, DissipatingNetworkError> {
        poll(&mut cancelled)?;
        controls.molecular.validate()?;
        let mut evaluations = 0;
        let (coupled, diagnostics) = self.step_method(
            maximum_duration,
            DissipatingStepControl {
                transport: controls.transport,
                molecular: controls.molecular.stabilized.transport,
                maximum_coupling_retries: controls.maximum_coupling_retries,
            },
            sources,
            reservoirs,
            MolecularMethod::Adaptive {
                controls: controls.molecular,
                evaluations: &mut evaluations,
            },
            &mut cancelled,
        )?;
        Ok(AdaptiveDissipatingNetworkReport {
            coupled,
            dissipation: diagnostics
                .adaptive
                .expect("both molecular halves were adaptive"),
            operator_evaluations: evaluations,
        })
    }
}
