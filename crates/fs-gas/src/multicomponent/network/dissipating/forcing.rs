//! Native acceleration in complete molecular/wave/source transactions.

use super::*;

impl<L: TransportLaw> DissipatingNetwork1d<L> {
    /// Enclose accelerated finite-port transport in the original molecular halves.
    ///
    /// Actual source-stage mass, force-history time, heat, viscosity and species
    /// diffusion remain with their native owners. Any shortened transport trial
    /// restarts the entire composition, including finite-plenum inventories.
    pub fn step_with_body_forces(
        &mut self,
        maximum_duration: Time,
        controls: DissipatingStepControl,
        forces: NetworkBodyForces<'_>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<DissipatingNetworkReport, DissipatingNetworkError> {
        self.step_with_body_forces_and_sources(
            maximum_duration, controls, forces, None, NetworkSources::default(), cancelled,
        )
    }

    /// Compose branch acceleration, molecular transport, reservoirs, walls and valves.
    ///
    /// Source work is retained in transport.body_force_work_j and removed once
    /// from the complete conservation residual. Molecular dissipation changes
    /// internal energy; it is not an additional exterior energy supply. All
    /// fields, finite-wall storage, source reports and clocks publish together.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_body_forces_and_sources(
        &mut self,
        maximum_duration: Time,
        controls: DissipatingStepControl,
        forces: NetworkBodyForces<'_>,
        reservoirs: Option<&[[Option<Cell<'_>>; 2]]>,
        sources: NetworkSources<'_, '_>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<DissipatingNetworkReport, DissipatingNetworkError> {
        self.step_forcing_method(
            maximum_duration, controls, sources, reservoirs, Some(forces),
            MolecularMethod::Explicit(controls.molecular), cancelled,
        ).map(|(report, _)| report)
    }

    /// Native RKL2 molecular transport with the same distributed body-force law.
    ///
    /// The original molecular RHS allowance spans accepted halves and all
    /// retries; exhaustion is terminal. Acceleration extrema still constrain
    /// acoustic CFL. Sources are evaluated on accepted physical time, never
    /// on an artificial molecular stage clock. All-zero forcing reduces to the
    /// corresponding original stabilized composition.
    #[allow(clippy::too_many_arguments)]
    pub fn step_stabilized_with_body_forces(
        &mut self,
        maximum_duration: Time,
        controls: StabilizedDissipatingStepControl,
        forces: NetworkBodyForces<'_>,
        reservoirs: Option<&[[Option<Cell<'_>>; 2]]>,
        sources: NetworkSources<'_, '_>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<StabilizedDissipatingNetworkReport, DissipatingNetworkError> {
        poll(&mut cancelled)?;
        controls.molecular.validate()?;
        let mut evaluations = 0;
        let (coupled, diagnostics) = self.step_forcing_method(
            maximum_duration,
            DissipatingStepControl {
                transport: controls.transport,
                molecular: controls.molecular.transport,
                maximum_coupling_retries: controls.maximum_coupling_retries,
            },
            sources, reservoirs, Some(forces),
            MolecularMethod::Stabilized {
                controls: controls.molecular, evaluations: &mut evaluations,
            },
            &mut cancelled,
        )?;
        Ok(StabilizedDissipatingNetworkReport {
            coupled,
            dissipation: diagnostics.stabilized.expect("both molecular halves used RKL2"),
            operator_evaluations: evaluations,
        })
    }
}
