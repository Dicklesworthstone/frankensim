//! Physical endpoint entropy for the original finite-network RKL2 owner.

use super::*;
use crate::multicomponent::dissipation::StabilizedDissipationControl;
use crate::multicomponent::network::dissipation::StabilizedNetworkDissipationReport;

impl<L: TransportLaw> NetworkDissipationOperator<L> {
    /// Advance a closed molecular network with source-conjugate RKL2 transport.
    ///
    /// Both branch interiors and finite-plenum ports retain the original
    /// Maxwell--Stefan, divided enthalpy, Fourier, reciprocal cross-effect and
    /// stress-work laws. The existing Legendre recurrence advances every finite
    /// inventory and its signed port transfers with the same stage weights.
    ///
    /// Each physical macrostep endpoint and the direct original-to-final
    /// interval must satisfy the declared entropy and species-continuity bounds.
    /// Internal Legendre stages retain native constitutive, caloric and rate
    /// admission; their signed recurrence is not interpreted as a sequence of
    /// independently entropy-producing physical substeps. Refused endpoints
    /// restart from the unchanged anchor under the original halving controls.
    ///
    /// Native endpoint temperatures supply entropy without another EOS inverse,
    /// flattened volume array, stage allocation or extra constitutive RHS.
    /// Stage/attempt/RHS budgets and the caller's physical macrostep accuracy cap
    /// remain authoritative. The native extended linear stability interval is
    /// not a nonlinear positivity, entropy or temporal-error theorem.
    ///
    /// Geometry is stationary and constant-area; exterior molecular species
    /// and heat fluxes are zero. Every active face requires common species
    /// support. This entry excludes reservoirs, prescribed motion and sources.
    /// Cancellation or refusal preserves all accepted fields, clocks and reports.
    pub fn advance_stabilized_entropy_consistent_for(
        &mut self,
        network: &mut MixtureNetwork1d,
        duration: Time,
        controls: StabilizedDissipationControl,
        entropy: MolecularEntropyControl,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(StabilizedNetworkDissipationReport, MolecularEntropyReport), DissipationError>
    {
        let report = self.advance_stabilized_entropy_for_trial(
            network,
            duration,
            controls,
            entropy,
            &mut 0,
            &mut cancelled,
        )?;
        self.publish_trial_report();
        Ok(report)
    }

    /// A composing owner retains one RHS budget across both halves and retries,
    /// and publishes these receipts only after complete physical admission.
    #[allow(clippy::too_many_arguments)]
    pub(in crate::multicomponent::network) fn advance_stabilized_entropy_for_trial(
        &mut self,
        network: &mut MixtureNetwork1d,
        duration: Time,
        controls: StabilizedDissipationControl,
        entropy: MolecularEntropyControl,
        evaluations: &mut usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(StabilizedNetworkDissipationReport, MolecularEntropyReport), DissipationError>
    {
        controls.validate()?;
        self.admit_plenum_volumes(network, cancelled)?;
        let mut account = initial_account(network, entropy, cancelled)?;
        self.trial_report = None;
        let report = super::super::stabilized::advance_entropy_consistent(
            self,
            network,
            duration,
            controls,
            evaluations,
            &mut account,
            cancelled,
        )?;
        self.trial_report = Some(report.transport);
        Ok((report, account.report()))
    }
}
