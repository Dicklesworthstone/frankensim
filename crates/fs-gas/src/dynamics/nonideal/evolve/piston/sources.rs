//! Distributed mechanical work on the original finite-piston transaction.

use super::*;

impl RealGasEuler1d {
    /// Advance distributed body-force work, the real-gas field and a finite piston.
    ///
    /// Original mechanical half-kicks surround the pressure/inertia solve, so
    /// the first kick changes the pressure waves and piston response during
    /// this interval. Each kick supplies its actual impulse times midpoint gas
    /// velocity as total-energy work, without changing internal energy or
    /// adding that work again as heat. The second kick uses transported density
    /// and the actual moved cell volumes. `gas.source_change` retains the
    /// extensive impulse, signed work and energy roundoff from both kicks.
    ///
    /// Acceleration and force density are authored external forcing. No
    /// gravitational potential inventory is inferred. Force arrays follow
    /// mesh cell order; each value stays fixed over the attempted interval.
    /// Zero forcing follows [`Self::step_with_piston`] without source scratch.
    ///
    /// # Errors
    /// Invalid force controls/shape, an inadmissible source endpoint, pressure
    /// solve or full-interval CFL refusal, exhausted retries and cancellation
    /// leave accepted gas, thermodynamic cache, geometry, clock and piston
    /// unchanged. A refused source endpoint shortens the complete interval.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_piston_and_body_force(
        &mut self,
        piston: &mut Piston,
        maximum_duration: Time,
        gas_controls: StepControl,
        load: PistonLoad,
        controls: PistonControl,
        force: BodyForce<'_>,
        force_controls: ForceControl,
        cancelled: impl FnMut() -> bool,
    ) -> Result<PistonStep, PistonError> {
        self.step_with_piston_sources(
            piston,
            maximum_duration,
            gas_controls,
            load,
            None,
            controls,
            force,
            force_controls,
            None,
            cancelled,
        )
        .map(|step| step.piston)
    }

    /// Compose body work, optional thermal walls and optional piston suspension.
    ///
    /// The order is body-force half, wall-source half, original pressure/inertia
    /// solve, wall-source half, body-force half. The first source endpoint is
    /// computed once per duration and reused by all midpoint-speed proposals.
    /// Native wall owners retain finite/shared wall storage, external reservoir
    /// heat, sliding-wall work and slip heating. The spring stores energy and
    /// piston damping remains a separate loss. None of these exchanges replaces
    /// pressure work or the independently reported distributed body work.
    ///
    /// Gas, piston, mesh, wall inventories and reports publish together. Mach
    /// admission applies at each body-force endpoint, and the complete final
    /// state must satisfy the full-interval relative-face CFL. Source buffers
    /// allocate once and are reused. Backward-Euler wall heat retains its first
    /// order accuracy despite the symmetric source ordering.
    ///
    /// The original constant-area, impermeable-end, frozen-composition real-gas
    /// restrictions remain. Chemistry, phase flash and evolving force fields
    /// are not inferred. See [`Self::step_with_piston_and_body_force`] for the
    /// authored forcing convention and [`Self::step_with_piston_and_walls`] for
    /// thermal wall conventions.
    ///
    /// # Errors
    /// All original source, wall binding, pressure, suspension and work-balance
    /// admissions apply. Refusal or cancellation publishes none of the owners.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_piston_sources(
        &mut self,
        piston: &mut Piston,
        maximum_duration: Time,
        gas_controls: StepControl,
        load: PistonLoad,
        spring: Option<PistonSuspension>,
        controls: PistonControl,
        force: BodyForce<'_>,
        force_controls: ForceControl,
        walls: Option<(&mut RealGasWallOperator<'_>, WallControl)>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<SuspendedPistonStep, PistonError> {
        self.step_with_piston_suspension(
            piston,
            maximum_duration,
            gas_controls,
            load,
            spring,
            controls,
            walls,
            Some((force, force_controls)),
            cancelled,
        )
    }
}
