//! RKL2 molecular transport in the original reacting network transaction.

use super::*;
use crate::multicomponent::dissipation::StabilizedDissipationControl;

/// Accepted chemical/wave/source balances and both actual RKL2 intervals.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StabilizedReactiveNetworkReport {
    pub coupled: ReactiveNetworkReport,
    pub dissipation: SplitStabilizedNetworkDissipationReport,
    /// All actual molecular RHS calls, including discarded coupled attempts.
    pub operator_evaluations: usize,
}

impl ReactingNetwork1d {
    /// Advance D/2, C/2, A, C/2, D/2 with RKL2 on both molecular intervals.
    ///
    /// The existing reactor and acoustic timestep own chemistry and wave
    /// evolution. Heat-release retries replay the complete composition from
    /// accepted fields, including the original finite-plenum inventories.
    /// The RKL2 RHS allowance is shared by both halves and all retries.
    pub fn step_with_stabilized_dissipation<L: TransportLaw>(
        &mut self,
        operator: &mut NetworkDissipationOperator<L>,
        maximum_duration: Time,
        controls: ReactiveStepControl,
        dissipation_controls: StabilizedDissipationControl,
        cancelled: impl FnMut() -> bool,
    ) -> Result<StabilizedReactiveNetworkReport, ReactiveNetworkError> {
        self.step_with_stabilized_dissipation_and_sources(
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

    /// Compose stabilized diffusion with optional walls, valves, baths and
    /// prescribed plenum volume rates in one atomic accepted step.
    ///
    /// Reservoir layout and volume rates have the same meaning as the existing
    /// explicit network APIs. Each diffusion half uses its actual stage volume;
    /// moving-volume advection retains pressure work exactly once. Prescribed
    /// baths reach both molecular halves and characteristic advection. Source
    /// geometry, finite walls and reports publish only after complete success.
    ///
    /// RKL2 extends the estimated diffusion stability interval, not the acoustic
    /// timestep or chemistry validity. Its maximum macrostep remains an explicit
    /// accuracy choice. No global temporal-order or entropy claim is inferred
    /// for the complete source/chemistry/advection composition.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_stabilized_dissipation_and_sources<L: TransportLaw>(
        &mut self,
        operator: &mut NetworkDissipationOperator<L>,
        maximum_duration: Time,
        controls: ReactiveStepControl,
        dissipation_controls: StabilizedDissipationControl,
        sources: NetworkSources<'_, '_>,
        reservoirs: Option<&[[Option<Cell<'_>>; 2]]>,
        volume_rates_m3_s: Option<&[f64]>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<StabilizedReactiveNetworkReport, ReactiveNetworkError> {
        poll(&mut cancelled)?;
        dissipation_controls.validate()?;
        let mut evaluations = 0;
        let (coupled, dissipation) = self.step_method(
            maximum_duration,
            controls,
            Some(operator),
            MolecularMethod::Stabilized {
                controls: dissipation_controls,
                evaluations: &mut evaluations,
            },
            sources,
            volume_rates_m3_s,
            reservoirs,
            &mut cancelled,
        )?;
        Ok(StabilizedReactiveNetworkReport {
            coupled,
            dissipation: dissipation.expect("both molecular halves used RKL2"),
            operator_evaluations: evaluations,
        })
    }
}
