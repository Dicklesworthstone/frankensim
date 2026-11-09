//! Body-force histories inside the original complete reacting transaction.

use super::*;
use crate::multicomponent::dissipation::StabilizedDissipationControl;

impl ReactingNetwork1d {
    /// Advance chemistry and distributed branch acceleration on one accepted clock.
    ///
    /// The original C/2,A,C/2 owner encloses the native forced transport map.
    /// Chemical source states and both mechanical half intervals are rebuilt
    /// from accepted storage whenever acoustic, force or chemistry admission
    /// shortens the duration. Source work remains in transport.body_force_work_j;
    /// formation-energy change remains the original chemical account.
    pub fn step_with_body_forces(
        &mut self,
        maximum_duration: Time,
        controls: ReactiveStepControl,
        forces: NetworkBodyForces<'_>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<ReactiveNetworkReport, ReactiveNetworkError> {
        self.step_with_body_forces_and_sources(
            maximum_duration, controls, forces, NetworkSources::default(), None, None, cancelled,
        )
    }

    /// Couple acceleration to native chemistry, finite walls/valves and boundaries.
    ///
    /// Optional reservoirs and plenum volume rates retain the original source
    /// ownership. No chemistry, wall storage, valve report, gas field or clock
    /// publishes before the complete accepted composition. All-zero forces
    /// retain the corresponding original composition arithmetic.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_body_forces_and_sources(
        &mut self,
        maximum_duration: Time,
        controls: ReactiveStepControl,
        forces: NetworkBodyForces<'_>,
        sources: NetworkSources<'_, '_>,
        reservoirs: Option<&[[Option<Cell<'_>>; 2]]>,
        volume_rates_m3_s: Option<&[f64]>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<ReactiveNetworkReport, ReactiveNetworkError> {
        self.step_forcing_method::<ConstantTransport>(
            maximum_duration, controls, None,
            MolecularMethod::Explicit(DissipationControl::default()), sources,
            volume_rates_m3_s, reservoirs, Some(forces), None, cancelled,
        ).map(|(report, _)| report)
    }

    /// Add native molecular diffusion, heat and viscosity to the forced reacting step.
    ///
    /// Each D/2,C/2,forced-A,C/2,D/2 candidate uses actual composition, volume
    /// and temperature at that stage. Molecular port enthalpy and viscous
    /// heating remain internal transfers; distributed body work is counted
    /// once in the complete physical energy balance.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_dissipation_and_body_forces<L: TransportLaw>(
        &mut self,
        operator: &mut NetworkDissipationOperator<L>,
        maximum_duration: Time,
        controls: ReactiveStepControl,
        dissipation_controls: DissipationControl,
        forces: NetworkBodyForces<'_>,
        sources: NetworkSources<'_, '_>,
        reservoirs: Option<&[[Option<Cell<'_>>; 2]]>,
        volume_rates_m3_s: Option<&[f64]>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<ReactiveNetworkReport, ReactiveNetworkError> {
        self.step_forcing_method(
            maximum_duration, controls, Some(operator),
            MolecularMethod::Explicit(dissipation_controls), sources,
            volume_rates_m3_s, reservoirs, Some(forces), None, cancelled,
        ).map(|(report, _)| report)
    }

    /// Use the original RKL2 molecular stages with the same acceleration history.
    ///
    /// The molecular RHS allowance spans both halves and every discarded
    /// complete attempt. RKL2 extends diffusion stability only; actual source
    /// speed extrema still constrain the common acoustic duration. Linear
    /// histories use that accepted physical interval, not an RKL stage index.
    /// No global entropy, hydrostatic balance or temporal-order theorem is inferred.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_stabilized_dissipation_and_body_forces<L: TransportLaw>(
        &mut self,
        operator: &mut NetworkDissipationOperator<L>,
        maximum_duration: Time,
        controls: ReactiveStepControl,
        dissipation_controls: StabilizedDissipationControl,
        forces: NetworkBodyForces<'_>,
        sources: NetworkSources<'_, '_>,
        reservoirs: Option<&[[Option<Cell<'_>>; 2]]>,
        volume_rates_m3_s: Option<&[f64]>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<StabilizedReactiveNetworkReport, ReactiveNetworkError> {
        poll(&mut cancelled)?;
        dissipation_controls.validate()?;
        let mut evaluations = 0;
        let (coupled, diagnostics) = self.step_forcing_method(
            maximum_duration, controls, Some(operator),
            MolecularMethod::Stabilized {
                controls: dissipation_controls, evaluations: &mut evaluations,
            },
            sources, volume_rates_m3_s, reservoirs, Some(forces), None, &mut cancelled,
        )?;
        Ok(StabilizedReactiveNetworkReport {
            coupled,
            dissipation: diagnostics.stabilized.expect("both molecular halves used RKL2"),
            operator_evaluations: evaluations,
        })
    }
}
