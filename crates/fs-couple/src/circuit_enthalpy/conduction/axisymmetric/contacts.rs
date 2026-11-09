//! Explicit source-bound surface resistance and two-sided Joule heat.
//!
//! Uses the original spatial contact owner for source pinning, actual thermal
//! and phase coordinates, receiving-cell admission and authored heat fractions.
//! Bulk and surface heat enter one backward-Euler finite-enthalpy residual.
//! Swept contact measures and radial basis moments remain native field data.

use super::*;

impl ConductionThermalNetwork<'_> {
    /// Advance axisymmetric bulk current, finite contact voltage jumps and
    /// both receiving bodies' enthalpies to one accepted thermal endpoint.
    ///
    /// Every native contact requires its original `ConductionContactThermal`
    /// binding and explicit fraction of heat entering surface A. The remainder
    /// enters B; receivers come from the original adjacent volume routes.
    /// The native field's surface heat is additional to its bulk cell powers
    /// and is deposited exactly once, with no inferred split or interface store.
    ///
    /// Interface temperature and optional liquid mass fractions are refreshed
    /// from the same original caloric charts as the bulk material. Original
    /// surface order, card/source choices, fixed pressure/history coordinates
    /// and declared phase-axis ownership remain pinned. Contact indices follow
    /// the native owner's canonical `contacts()` order.
    ///
    /// All four `ConductionDrive` modes and original field/enthalpy limits
    /// apply. Matching stationary traces and constant volume/density retain
    /// the native owners' limits; no contact motion, spreading or arcing law
    /// is supplied. Original Fourier links remain independently authored.
    /// The returned thermal and electrical owners must be published together;
    /// any source, balance, solve or cancellation failure changes neither input.
    #[allow(clippy::too_many_arguments)]
    pub fn step_axisymmetric_contact_candidate<'t, C: PlanarConductorMaterial>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        problem: &AxisymmetricConduction<C>,
        drive: ConductionDrive<'_>,
        contacts: &[ConductionContactThermal<'_>],
        assembly: WeightedAssemblyLimits,
        field_controls: ConductionControls,
        controls: CircuitEnthalpySolve,
        initial_energy_increments_j: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<AxisymmetricConductionEnthalpyCandidate<'t, C>, Error> {
        let cancel = RefCell::new(cancelled);
        self.step_axisymmetric_with_contacts_candidate(
            thermal,
            target_s,
            problem,
            drive,
            contacts,
            assembly,
            field_controls,
            controls,
            initial_energy_increments_j,
            &cancel,
        )
    }
}
