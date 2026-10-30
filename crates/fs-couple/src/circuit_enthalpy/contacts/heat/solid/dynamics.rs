//! Complete original heat transport on a retained inertial J2 trajectory.

use super::*;
use crate::solid_enthalpy::plastic::dynamics::{
    ThermoplasticDynamics, ThermoplasticDynamicsCandidate,
};

impl<'m> ThermoplasticDynamics<'m> {
    /// Advance source-bound inertia, plastic history and the complete material
    /// heat network as one candidate.
    ///
    /// Conduction, thermal interfaces, integrated layers and enclosure radiation
    /// use the existing prepared network at the actual caloric substeps. Their
    /// original cards, caloric objects and fixed material conditions must agree
    /// with the moving solid before any proposal is made. Additional finite
    /// thermal bodies remain in the same thermal timeline.
    ///
    /// This calls the native inertial thermoplastic solve, retaining full nodal
    /// velocity, displacement, J2 history and the physical clock. Physical
    /// yield/overstress loss and reversible heat enter that one thermal solve.
    /// Consistent reference-mass kinetic energy, hardening storage, constitutive
    /// and inertial numerical losses keep their original separate accounts.
    ///
    /// The additional closed audit excludes internal radiative roundoff from
    /// external supply and checks both corrected residuals against the same
    /// energy tolerance. No second mechanical or thermal advance is performed.
    /// Electrical bindings refuse because this entry point owns no circuit work.
    ///
    /// On continuation, pass the returned observer's `second` history as
    /// `inner`; the prepared network supplies exactly one heat observer each
    /// interval. Accept the returned thermal state and dynamics together.
    /// Source refusal, observer failure, cancellation and failed energy gates
    /// leave all inputs unchanged.
    ///
    /// Native small-strain/rotation, fully solid, fixed reference geometry and
    /// source-range restrictions still apply. Thermal substeps do not refine
    /// inertia or constitutive loading: refine physical intervals separately.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub fn advance_heat_network_candidate<'t, 'p, 'n, 'a, 'scope, O: DepositionObserver>(
        &self,
        thermal: &DepositionSimulation<'t>,
        prepared: &'p PreparedHeatNetwork<'n, 'a>,
        inner: &O,
        target_s: f64,
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        settings: IncrementSettings,
        thermal_config: StepConfig,
        config: ThermoplasticConfig,
        cx: &'p Cx<'scope>,
    ) -> Result<
        (
            ThermoplasticDynamicsCandidate<
                'm,
                't,
                CoupledObservers<MaterialHeatObserver<'p, 'n, 'a, 'scope>, O>,
            >,
            SolidHeatAudit,
        ),
        SolidHeatError<O::Error>,
    > {
        cx.checkpoint().map_err(|_| SolidHeatError::Cancelled)?;
        self.solid().admit_heat_network(thermal, prepared, cx)?;
        let observer = prepared
            .observe(cx, inner.clone(), config.energy_tolerance_j)
            .map_err(SolidHeatError::Heat)?;
        let candidate = self
            .advance(
                thermal,
                &observer,
                target_s,
                prescribed_m,
                nodal_forces_n,
                settings,
                thermal_config,
                config,
                cx,
            )
            .map_err(SolidHeatError::Coupled)?;
        let report = &candidate.report().material;
        let audit = closed_account(
            report.external_heat_j,
            report.thermal_residual_j,
            report.combined_residual_j,
            candidate.observer().first.closed_radiation_net_j(),
            config.energy_tolerance_j,
        )
        .map_err(SolidHeatError::Invalid)?;
        cx.checkpoint().map_err(|_| SolidHeatError::Cancelled)?;
        Ok((candidate, audit))
    }
}
