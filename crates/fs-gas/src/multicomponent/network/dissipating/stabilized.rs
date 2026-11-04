//! Stabilized transport through the original molecular/advection transaction.

use super::*;
use crate::multicomponent::dissipation::StabilizedDissipationControl;

/// Acoustic and complete-composition controls with explicitly selected RKL2.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StabilizedDissipatingStepControl {
    pub transport: StepControl,
    /// Its RHS allowance covers both halves and every discarded coupled attempt.
    pub molecular: StabilizedDissipationControl,
    pub maximum_coupling_retries: usize,
}
impl Default for StabilizedDissipatingStepControl {
    fn default() -> Self {
        Self {
            transport: StepControl::default(),
            molecular: StabilizedDissipationControl::default(),
            maximum_coupling_retries: 20,
        }
    }
}

/// Complete accepted D/A/D balance and the actual work of its RKL2 stages.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StabilizedDissipatingNetworkReport {
    pub coupled: DissipatingNetworkReport,
    pub dissipation: SplitStabilizedNetworkDissipationReport,
    /// Includes accepted halves and all work discarded by coupled retries.
    pub operator_evaluations: usize,
}

impl<L: TransportLaw> DissipatingNetwork1d<L> {
    /// Advance D/2, A, D/2 using the shared network RKL2 operator.
    ///
    /// Finite-plenum inventories, constitutive reevaluation and conservative
    /// port transfers participate in every molecular stage. The original
    /// acoustic timestep still controls the complete composition. RKL2's
    /// extended stability interval is not a temporal accuracy or positivity
    /// certificate; callers retain an accuracy-appropriate maximum macrostep.
    pub fn step_stabilized(
        &mut self,
        maximum_duration: Time,
        controls: StabilizedDissipatingStepControl,
        cancelled: impl FnMut() -> bool,
    ) -> Result<StabilizedDissipatingNetworkReport, DissipatingNetworkError> {
        self.step_stabilized_with_sources(
            maximum_duration,
            controls,
            None,
            NetworkSources::default(),
            cancelled,
        )
    }

    /// Compose RKL2 with optional prescribed baths, finite walls and valves.
    ///
    /// Baths drive both molecular halves and characteristic advection. Every
    /// source and molecular half replays the same accepted duration on retry.
    /// The actual RHS allowance never resets between halves or retries; its
    /// exhaustion is terminal and publishes no state or source report.
    pub fn step_stabilized_with_sources(
        &mut self,
        maximum_duration: Time,
        controls: StabilizedDissipatingStepControl,
        reservoirs: Option<&[[Option<Cell<'_>>; 2]]>,
        sources: NetworkSources<'_, '_>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<StabilizedDissipatingNetworkReport, DissipatingNetworkError> {
        poll(&mut cancelled)?;
        controls.molecular.validate()?;
        let mut evaluations = 0;
        let (coupled, dissipation) = self.step_method(
            maximum_duration,
            DissipatingStepControl {
                transport: controls.transport,
                molecular: controls.molecular.transport,
                maximum_coupling_retries: controls.maximum_coupling_retries,
            },
            sources,
            reservoirs,
            MolecularMethod::Stabilized {
                controls: controls.molecular,
                evaluations: &mut evaluations,
            },
            &mut cancelled,
        )?;
        Ok(StabilizedDissipatingNetworkReport {
            coupled,
            dissipation: dissipation.expect("both molecular halves used RKL2"),
            operator_evaluations: evaluations,
        })
    }
}
