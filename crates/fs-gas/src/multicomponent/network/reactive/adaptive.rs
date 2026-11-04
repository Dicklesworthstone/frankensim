//! Adaptive molecular intervals composed with existing chemistry and sources.

use super::super::dissipation::AdaptiveNetworkDissipationControl;
use super::*;

/// Complete accepted network balance and both molecular temporal estimates.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AdaptiveReactiveNetworkReport {
    pub coupled: ReactiveNetworkReport,
    pub dissipation: SplitAdaptiveNetworkDissipationReport,
    /// Actual molecular RHS calls, including all probes and discarded compositions.
    pub operator_evaluations: usize,
}

impl ReactingNetwork1d {
    /// Advance D/2, C/2, A, C/2, D/2 with adaptive RKL2 molecular intervals.
    ///
    /// Chemistry and wave evolution retain their original owners and controls.
    /// Heat-release retries restart the complete composition and preserve all
    /// molecular work already charged to the complete-step allowance.
    pub fn step_with_adaptive_dissipation<L: TransportLaw>(
        &mut self,
        operator: &mut NetworkDissipationOperator<L>,
        maximum_duration: Time,
        controls: ReactiveStepControl,
        dissipation_controls: AdaptiveNetworkDissipationControl,
        cancelled: impl FnMut() -> bool,
    ) -> Result<AdaptiveReactiveNetworkReport, ReactiveNetworkError> {
        self.step_with_adaptive_dissipation_and_sources(
            operator,
            maximum_duration,
            controls,
            dissipation_controls,
            NetworkSources::default(),
            None,
            None,
            cancelled,
        )
    }

    /// Compose molecular error control, chemistry, baths, sources and plenum motion.
    ///
    /// Each molecular interval evaluates its error at the actual current plenum
    /// volumes. Advection owns their subsequent motion and pressure work; the
    /// molecular estimator does not integrate a second copy of that work.
    /// Failed source, temporal, physical or cancellation checks publish no partial
    /// field, geometry, finite-wall storage or accepted diagnostic.
    ///
    /// The local estimator controls molecular time error. It does not certify
    /// the complete reaction/source/advection splitting error, spatial error,
    /// global accuracy, or entropy stability of the coupled method.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_adaptive_dissipation_and_sources<L: TransportLaw>(
        &mut self,
        operator: &mut NetworkDissipationOperator<L>,
        maximum_duration: Time,
        controls: ReactiveStepControl,
        dissipation_controls: AdaptiveNetworkDissipationControl,
        sources: NetworkSources<'_, '_>,
        reservoirs: Option<&[[Option<Cell<'_>>; 2]]>,
        volume_rates_m3_s: Option<&[f64]>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AdaptiveReactiveNetworkReport, ReactiveNetworkError> {
        poll(&mut cancelled)?;
        dissipation_controls.validate()?;
        let mut evaluations = 0;
        let (coupled, diagnostics) = self.step_method(
            maximum_duration,
            controls,
            Some(operator),
            MolecularMethod::Adaptive {
                controls: dissipation_controls,
                evaluations: &mut evaluations,
            },
            sources,
            volume_rates_m3_s,
            reservoirs,
            &mut cancelled,
        )?;
        Ok(AdaptiveReactiveNetworkReport {
            coupled,
            dissipation: diagnostics
                .adaptive
                .expect("both molecular halves were adaptive"),
            operator_evaluations: evaluations,
        })
    }
}
