//! Source-conjugate RKL2 inside the original complete reacting transaction.

use super::*;
use crate::multicomponent::dissipation::StabilizedDissipationControl;

/// Complete chemistry/transport entropy with both actual RKL2 half intervals.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StabilizedReactiveNetworkEntropyReport {
    /// Original coupled balances, component entropy and complete source entropy.
    pub coupled: ReactiveNetworkEntropyReport,
    /// Actual accepted molecular halves, including their native stability work.
    pub dissipation: SplitStabilizedNetworkDissipationReport,
    /// All actual molecular RHS calls, including discarded coupled attempts.
    pub operator_evaluations: usize,
}

impl ReactingNetwork1d {
    /// Compose source-conjugate RKL2, chemistry and advection with full entropy.
    ///
    /// The original D/2-C/2-A-C/2-D/2 retry owner advances all branch and plenum
    /// inventories together. Each molecular macrostep and complete half interval
    /// receives its native endpoint entropy gate. One immutable pre-step source
    /// account also checks the actual complete original-to-final state, retaining
    /// physical chemical entropy constants and removing only the reference
    /// contribution of bounded numerical source-continuity defects.
    ///
    /// Acoustic shortening or entropy refusal replays every half from accepted
    /// state. One caller-bounded actual RHS count spans both molecular intervals
    /// and every discarded composition; exhaustion is terminal. Both molecular
    /// receipts stay provisional until complete temperature, balance, entropy
    /// and final cancellation checks succeed. No reaction or viscous heat is
    /// added twice, and internal finite ports carry no exterior entropy supply.
    ///
    /// Closed fixed constant-area networks and common transport species support
    /// are required. The caller still chooses the molecular accuracy cap and
    /// chemistry controls. Neither signed Legendre stages nor symmetric splitting
    /// establish a global order, nonlinear positivity or entropy-stability theorem.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_stabilized_dissipation_and_entropy<L: TransportLaw>(
        &mut self,
        operator: &mut NetworkDissipationOperator<L>,
        maximum_duration: Time,
        controls: ReactiveStepControl,
        dissipation_controls: StabilizedDissipationControl,
        entropy: MixtureEntropyControl,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<StabilizedReactiveNetworkEntropyReport, ReactiveNetworkError> {
        poll(&mut cancelled)?;
        dissipation_controls.validate()?;
        let mut account = CoupledEntropy::new(&self.accepted, entropy, &mut cancelled)?;
        let mut evaluations = 0;
        let (coupled, diagnostics) = self.step_method(
            maximum_duration,
            controls,
            Some(operator),
            MolecularMethod::EntropyStabilized {
                controls: dissipation_controls,
                entropy,
                evaluations: &mut evaluations,
            },
            NetworkSources::default(),
            None,
            None,
            Some(&mut account),
            &mut cancelled,
        )?;
        Ok(StabilizedReactiveNetworkEntropyReport {
            coupled: account.complete(coupled),
            dissipation: diagnostics
                .stabilized
                .expect("both molecular halves used RKL2"),
            operator_evaluations: evaluations,
        })
    }
}
