//! Finite wall recoil on the original source-driven thermal liquid.

use super::{
    PlanarWallCandidate, PlanarWallControls, PlanarWallError, PlanarWallLiquid, PlanarWallMotion,
    Vec3, sealed,
};
use fs_exec::CancelGate;
use fs_lbm::material::free_surface::MaterialFreeSurface3;
use fs_lbm::material::free_surface::enthalpy::{
    FreeSurfaceCellHeat, FreeSurfaceEnthalpyReport, MaterialFreeSurfaceEnthalpy3,
    viscous::FreeSurfaceViscousCell,
};

/// Actual native thermal liquid and a guided finite body with two free tangent
/// translations and axial spin. All constructor and motion bounds are the
/// original planar mechanical owner's bounds; the full inertia is retained.
pub type FreeSurfaceThermalPlanarWall<'s, 'm> =
    PlanarWallMotion<MaterialFreeSurfaceEnthalpy3<'s, 'm>>;

/// Thermal receipt from the final native replay, never from the discarded pilot.
#[derive(Clone, Debug, PartialEq)]
pub struct ThermalPlanarWallReport {
    pub enthalpy: FreeSurfaceEnthalpyReport,
    /// Explicitly selected original initial-state molecular heat samples.
    /// `None` means this source was not selected, rather than inferred zero loss.
    pub molecular_heat: Option<Vec<FreeSurfaceViscousCell>>,
}

/// Original complete mechanical audit plus the final thermal source/transport
/// report. Continue from `motion`; its thermal clock has advanced exactly once.
pub type ThermalPlanarWallStep<'s, 'm> =
    PlanarWallCandidate<MaterialFreeSurfaceEnthalpy3<'s, 'm>, ThermalPlanarWallReport>;

impl sealed::NativeLiquid for MaterialFreeSurfaceEnthalpy3<'_, '_> {}
impl PlanarWallLiquid for MaterialFreeSurfaceEnthalpy3<'_, '_> {
    fn wall_flow(&self) -> &MaterialFreeSurface3<'_> {
        self.flow()
    }
}

impl<'s, 'm> PlanarWallMotion<MaterialFreeSurfaceEnthalpy3<'s, 'm>> {
    /// Advance held tangent force/couple and authored cell heat with a finite wall.
    ///
    /// One original rigid-link pilot measures all three momentum responses.
    /// The existing mass-whitened joint midpoint solve includes finite body
    /// recoil and translation/spin cross coupling. The final native thermal
    /// replay supplies the actual impulse, torque and work to the original MBD
    /// integrator; its endpoint velocity is never replaced by the solved guess.
    ///
    /// Both trials start from this unchanged source/enthalpy/population state.
    /// Heat, local viscosity/rheology and optional buoyancy act before the same
    /// native collision. Actual mass routes carry extensive enthalpy. Only the
    /// final body and fully admitted thermal candidate escape, with one tick
    /// on their accepted clocks. No prepared matrix/core is cloned.
    ///
    /// Wall work, applied/guide work, kinetic change and authored heat remain
    /// distinct. This is fixed flooded-plane tangential motion, with no normal
    /// coverage, finite moving edge, tilt or newly invented phase/heat boundary.
    /// Existing thermal source splitting and flow convergence limits remain.
    pub fn step_candidate(
        &self,
        external_force_n: Vec3,
        external_torque_nm: f64,
        heat: &[FreeSurfaceCellHeat],
        controls: PlanarWallControls,
        gate: &CancelGate,
    ) -> Result<ThermalPlanarWallStep<'s, 'm>, PlanarWallError> {
        self.step_candidate_with(
            external_force_n,
            external_torque_nm,
            controls,
            gate,
            |motion, mass_tolerance, gate| {
                let (candidate, report) = self.liquid.planar_rigid_wall_motion_candidate_si(
                    motion,
                    mass_tolerance,
                    heat,
                    gate,
                )?;
                Ok((
                    candidate,
                    report.wall,
                    ThermalPlanarWallReport {
                        enthalpy: report.flow,
                        molecular_heat: None,
                    },
                ))
            },
        )
    }

    /// Select the original held initial-state molecular stress-power heat too.
    /// The discarded pilot does not deposit a second copy of either heat source.
    /// Molecular heat is not reconstructed from actuator work or numerical loss.
    /// The original explicit source split requires timestep refinement when
    /// viscosity and temperature strongly feed back on each other.
    pub fn step_candidate_with_viscous_heating(
        &self,
        external_force_n: Vec3,
        external_torque_nm: f64,
        heat: &[FreeSurfaceCellHeat],
        controls: PlanarWallControls,
        gate: &CancelGate,
    ) -> Result<ThermalPlanarWallStep<'s, 'm>, PlanarWallError> {
        self.step_candidate_with(
            external_force_n,
            external_torque_nm,
            controls,
            gate,
            |motion, mass_tolerance, gate| {
                let (candidate, report) = self
                    .liquid
                    .planar_rigid_wall_motion_candidate_with_viscous_heating_si(
                        motion,
                        mass_tolerance,
                        heat,
                        gate,
                    )?;
                Ok((
                    candidate,
                    report.flow.wall,
                    ThermalPlanarWallReport {
                        enthalpy: report.flow.flow,
                        molecular_heat: Some(report.applied),
                    },
                ))
            },
        )
    }
}
