//! Both finite wall mechanisms share one original thermal liquid trajectory.

use super::{
    PlanarWallControls, PlanarWallPairCandidate, PlanarWallPairError, PlanarWallPairMotion,
};
use crate::free_surface_wall::planar::thermal::ThermalPlanarWallReport;
use fs_exec::CancelGate;
use fs_lbm::material::free_surface::enthalpy::{FreeSurfaceCellHeat, MaterialFreeSurfaceEnthalpy3};
use fs_mbd::Vec3;

/// Two finite bodies and the one original source-bound thermal liquid.
/// Each body retains two tangent translations and its principal normal spin.
pub type FreeSurfaceThermalPlanarWallPair<'s, 'm> =
    PlanarWallPairMotion<MaterialFreeSurfaceEnthalpy3<'s, 'm>>;

/// Original complete body/wall accounts and the final thermal source report.
/// The retained liquid advances one tick, regardless of its time at binding.
pub type ThermalPlanarWallPairStep<'s, 'm> =
    PlanarWallPairCandidate<MaterialFreeSurfaceEnthalpy3<'s, 'm>, ThermalPlanarWallReport>;

impl<'s, 'm> PlanarWallPairMotion<MaterialFreeSurfaceEnthalpy3<'s, 'm>> {
    /// Advance both finite walls, authored heat and native thermal flow together.
    ///
    /// One common pilot measures both original three-mode wall responses after
    /// the same heat/material preparation. The existing mass-whitened midpoint
    /// solves and one final pair replay supply both actual body endpoints. The
    /// two bodies retain their original impulse, spin/orbital, guide-reaction
    /// and work audits. The final thermal candidate and both bodies publish
    /// together; discarded pilots advance neither accepted heat nor time.
    ///
    /// Cell heat refreshes the original temperature-dependent viscosity,
    /// rheology and optional buoyancy before collision. The one native mass
    /// transport carries extensive enthalpy. Body and boundary mechanical work
    /// remain separate from caloric sources; this does not infer heat from
    /// actuator work or numerical loss. The original explicit source split and
    /// fixed flooded-plane geometry retain their existing limits.
    pub fn step_candidate(
        &self,
        external_forces_n: [Vec3; 2],
        external_torques_nm: [f64; 2],
        heat: &[FreeSurfaceCellHeat],
        controls: PlanarWallControls,
        gate: &CancelGate,
    ) -> Result<ThermalPlanarWallPairStep<'s, 'm>, PlanarWallPairError> {
        self.step_candidate_with(
            external_forces_n,
            external_torques_nm,
            controls,
            gate,
            |motions, mass_tolerance, gate| {
                let (candidate, report) = self.liquid.planar_rigid_wall_pair_candidate_si(
                    motions,
                    mass_tolerance,
                    heat,
                    gate,
                )?;
                Ok((
                    candidate,
                    report.walls,
                    ThermalPlanarWallReport {
                        enthalpy: report.flow,
                        molecular_heat: None,
                    },
                ))
            },
        )
    }

    /// Select the original held initial-state molecular stress-power heat.
    /// One common liquid source is admitted once in the accepted interval;
    /// neither the second body nor the discarded pilot deposits another copy.
    /// Strong viscosity/temperature feedback retains the original source-split
    /// timestep-refinement requirement.
    pub fn step_candidate_with_viscous_heating(
        &self,
        external_forces_n: [Vec3; 2],
        external_torques_nm: [f64; 2],
        heat: &[FreeSurfaceCellHeat],
        controls: PlanarWallControls,
        gate: &CancelGate,
    ) -> Result<ThermalPlanarWallPairStep<'s, 'm>, PlanarWallPairError> {
        self.step_candidate_with(
            external_forces_n,
            external_torques_nm,
            controls,
            gate,
            |motions, mass_tolerance, gate| {
                let (candidate, report) = self
                    .liquid
                    .planar_rigid_wall_pair_candidate_with_viscous_heating_si(
                        motions,
                        mass_tolerance,
                        heat,
                        gate,
                    )?;
                Ok((
                    candidate,
                    report.flow.walls,
                    ThermalPlanarWallReport {
                        enthalpy: report.flow.flow,
                        molecular_heat: Some(report.applied),
                    },
                ))
            },
        )
    }
}
