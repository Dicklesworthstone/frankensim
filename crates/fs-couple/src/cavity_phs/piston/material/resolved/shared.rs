//! Two resolved chambers acting on the original material-derived wall inertia.

use super::MaterialResolvedPiston;
use fs_gas::dynamics::{
    StepControl,
    piston::{PistonControl, PistonError, PistonStep, PistonSuspension, SuspendedPistonStep},
};
use fs_gas::multicomponent::{
    MixtureEuler1d, SharedPistonGasReport, SharedPistonThermal, SharedPistonThermalError,
    SharedPistonThermalStepReport, wall_exchange::WallControl,
};
use fs_qty::Time;

impl MaterialResolvedPiston {
    /// Advance two resolved gas chambers against this one source-bound piston.
    ///
    /// The original shared-wall solver admits matching interface positions,
    /// areas and clocks, and closes both pressure works against the SAME body
    /// inertia. Positive velocity expands `left` and compresses `right`.
    /// `external_force_n` excludes both resolved gas pressure forces.
    ///
    /// Both gas fields, meshes, clocks and wall velocity publish together.
    /// Original species axes, caloric sources and density receipts are retained.
    /// This is the native negligible-thickness, rigid axial partition model;
    /// the material specimen's finite geometry supplies inertia, not a new
    /// moving-mesh geometry or deformable wall.
    #[allow(clippy::too_many_arguments)]
    pub fn step_shared_mixture(
        &mut self,
        left: &mut MixtureEuler1d,
        right: &mut MixtureEuler1d,
        maximum_duration: Time,
        gas_controls: StepControl,
        external_force_n: f64,
        controls: PistonControl,
        cancelled: impl FnMut() -> bool,
    ) -> Result<PistonStep<SharedPistonGasReport>, PistonError> {
        left.step_with_shared_piston(
            right, &mut self.piston, maximum_duration, gas_controls,
            external_force_n, controls, cancelled,
        )
    }

    /// Resolve both gas pressures and a stationary spring/damper in one step.
    ///
    /// The native midpoint solve includes the support force in every trial.
    /// Spring storage and physical damping loss retain their original accounts;
    /// this entry sends no damper heat into either gas. A refused step or
    /// cancellation preserves both chambers and the material-bound piston.
    #[allow(clippy::too_many_arguments)]
    pub fn step_shared_mixture_suspended(
        &mut self,
        left: &mut MixtureEuler1d,
        right: &mut MixtureEuler1d,
        maximum_duration: Time,
        gas_controls: StepControl,
        external_force_n: f64,
        suspension: PistonSuspension,
        controls: PistonControl,
        cancelled: impl FnMut() -> bool,
    ) -> Result<SuspendedPistonStep<SharedPistonGasReport>, PistonError> {
        left.step_with_shared_suspended_piston(
            right, &mut self.piston, maximum_duration, gas_controls,
            external_force_n, suspension, controls, cancelled,
        )
    }

    /// Couple both resolved chambers to one finite thermal store and this body.
    ///
    /// The original heat/pressure/heat composition owns gas calorics, moved
    /// volumes, pressure work and wall storage. Both gas fields, the finite
    /// store and this piston publish together after their original admissions.
    /// The supplied thermal store is explicit; it does not re-resolve the
    /// rigid specimen's frozen density or infer expansion. When it represents
    /// the same specimen, the caller must bind its compatible caloric inventory.
    ///
    /// The native split-step accuracy and rigid partition scope apply. This
    /// adds no second pressure, heat, geometry or material evolution solver.
    #[allow(clippy::too_many_arguments)]
    pub fn step_shared_mixture_with_heat(
        &mut self,
        left: &mut MixtureEuler1d,
        right: &mut MixtureEuler1d,
        thermal: &mut SharedPistonThermal<'_>,
        maximum_duration: Time,
        gas_controls: StepControl,
        heat_controls: WallControl,
        external_force_n: f64,
        controls: PistonControl,
        cancelled: impl FnMut() -> bool,
    ) -> Result<PistonStep<SharedPistonThermalStepReport>, SharedPistonThermalError> {
        left.step_with_shared_piston_and_heat(
            right, thermal, &mut self.piston, maximum_duration, gas_controls,
            heat_controls, external_force_n, controls, cancelled,
        )
    }

    /// Return an explicit share of suspension loss to the original finite wall.
    ///
    /// The native suspended heat composition uses the same sourced inertia for
    /// both pressure fields, kinetic storage and spring/damper feedback. The
    /// authored fraction in [0, 1] deposits accepted physical damper loss once;
    /// its remainder leaves the represented thermal system. Heat then changes
    /// actual gas state and subsequent motion through the original wall owner.
    /// A refused trial never publishes gas, thermal or mechanical endpoints.
    ///
    /// The separate thermal-inventory and fixed-material assumptions of
    /// [`Self::step_shared_mixture_with_heat`] still apply.
    #[allow(clippy::too_many_arguments)]
    pub fn step_shared_mixture_suspended_with_heat(
        &mut self,
        left: &mut MixtureEuler1d,
        right: &mut MixtureEuler1d,
        thermal: &mut SharedPistonThermal<'_>,
        maximum_duration: Time,
        gas_controls: StepControl,
        heat_controls: WallControl,
        external_force_n: f64,
        suspension: PistonSuspension,
        damping_heat_fraction_to_wall: f64,
        controls: PistonControl,
        cancelled: impl FnMut() -> bool,
    ) -> Result<SuspendedPistonStep<SharedPistonThermalStepReport>, SharedPistonThermalError> {
        left.step_with_shared_suspended_piston_and_heat(
            right, thermal, &mut self.piston, maximum_duration, gas_controls,
            heat_controls, external_force_n, suspension, damping_heat_fraction_to_wall,
            controls, cancelled,
        )
    }
}
