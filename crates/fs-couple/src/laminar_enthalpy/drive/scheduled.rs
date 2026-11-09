//! Exact authored mechanical events on the existing joint adaptive controller.

use super::super::adaptive::{
    LaminarAdaptiveCandidate, LaminarAdaptiveControls, LaminarAdaptiveError,
};
use super::*;

/// Final native state plus accepted mechanical-source and thermal receipts.
pub type LaminarScheduledCandidate<'m, D, G, T> =
    LaminarAdaptiveCandidate<'m, G, LaminarDrivenExchange<D, T>>;

impl<'s, 'mesh, 'source, F: FluidTransport + ?Sized> LaminarEnthalpy<'s, 'mesh, 'source, F> {
    /// Advance prescribed and temperature-dependent mechanics with exact
    /// authored event stops. Events are finite, nonnegative, strictly ordered
    /// absolute times; past/future entries may remain for deterministic resume.
    /// A query at a closing event retains its preceding segment. The following
    /// interval uses the next segment without shifting either physical clock.
    ///
    /// Every coarse, first-fine and last-fine thermal candidate resolves its
    /// own actual mechanical snapshot. Equal times with different candidate
    /// thermal states must not reuse a lagged drive. Within each segment the
    /// caller authors smooth data and fixed boundary kinds; discontinuities
    /// and kind changes require explicit event times. Finite nonlinear probes
    /// must stay in the original source's declared domain.
    ///
    /// The original adaptive error/work/heat gates remain authoritative.
    /// Body, native prescribed-wall/supplied-traction/backflow work and signed
    /// carried-kinetic exchanges enter separate transfer-error channels with
    /// signed and resolved-magnitude defects. Numerical upwind loss stays
    /// separate; none of these observations adds a thermal source.
    /// No thermal inventory event is activated and no accepted prefix escapes
    /// source, native, accuracy, work or late cancellation failure.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub fn advance_scheduled<'m, D: LaminarMechanicalDrive, G: ThermalFluidTransport, T, E>(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'m>,
        target_s: f64,
        history: LaminarEnthalpyHistory<'_>,
        event_times_s: &[f64],
        controls: LaminarAdaptiveControls,
        resolve_drive: impl Fn(&Cx<'_>, LaminarDriveQuery, &DepositionSimulation<'m>) -> Result<D, E>,
        resolve_material: impl Fn(MaterialFluxRheologyPoint, f64, &F, &D) -> Result<G, MaterialError>,
        exchange: impl Fn(
            &Cx<'_>,
            &DepositionSimulation<'m>,
            &MaterialNavierStokesRheologyStep<LaminarThermalMaterial<G>>,
            &D,
        ) -> Result<LaminarHeatExchange<T>, E>,
    ) -> Result<LaminarScheduledCandidate<'m, D, G, T>, LaminarAdaptiveError<E>> {
        self.advance_adaptive_native(
            cx,
            thermal,
            target_s,
            history,
            event_times_s,
            controls,
            |old, query, accepted, step, guess| {
                self.advance_driven(
                    cx,
                    old,
                    query,
                    accepted,
                    step,
                    guess,
                    &resolve_drive,
                    &resolve_material,
                    &exchange,
                )
            },
        )
    }

    /// The same scheduled mechanics with original density-driven natural
    /// convection at every joint endpoint. Buoyancy is resolved separately
    /// from the drive's additional force and never duplicated as source heat.
    /// Hydrostatic reference inventory and reduced pressures retain the
    /// original advance_boussinesq_driven conventions.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub fn advance_boussinesq_scheduled<
        'm,
        D: LaminarMechanicalDrive,
        G: ThermalFluidTransport,
        T,
        E,
    >(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'m>,
        target_s: f64,
        history: LaminarEnthalpyHistory<'_>,
        event_times_s: &[f64],
        config: boussinesq::BoussinesqConfig<2>,
        controls: LaminarAdaptiveControls,
        resolve_drive: impl Fn(&Cx<'_>, LaminarDriveQuery, &DepositionSimulation<'m>) -> Result<D, E>,
        resolve_material: impl Fn(MaterialFluxRheologyPoint, f64, &F, &D) -> Result<G, MaterialError>,
        exchange: impl Fn(
            &Cx<'_>,
            &DepositionSimulation<'m>,
            &MaterialNavierStokesRheologyStep<
                LaminarThermalMaterial<boussinesq::LaminarBoussinesqProjection<G>>,
            >,
            &D,
        ) -> Result<LaminarHeatExchange<T>, E>,
    ) -> Result<
        LaminarScheduledCandidate<'m, D, boussinesq::LaminarBoussinesqProjection<G>, T>,
        LaminarAdaptiveError<E>,
    > {
        self.advance_adaptive_native(
            cx,
            thermal,
            target_s,
            history,
            event_times_s,
            controls,
            |old, query, accepted, step, guess| {
                self.advance_boussinesq_driven(
                    cx,
                    old,
                    query,
                    accepted,
                    config,
                    step,
                    guess,
                    &resolve_drive,
                    &resolve_material,
                    &exchange,
                )
            },
        )
    }
}
