//! Source-bound SI inputs and receipts for the original wall/WALE/Vreman step.

use fs_exec::CancelGate;

use super::{
    MaterialFreeSurface3, MaterialFreeSurfaceScaling, MaterialPlanarWallFieldReport3,
    MaterialPlanarWallVelocityField3,
};
use crate::d3q19::freesurface3::{
    FreeSurfaceError3, FreeSurfaceKinematicsWorkspace3, PlanarWallFieldControls3,
    PlanarWallVelocityField3,
};
use crate::d3q19::freesurface3::sgs::{
    FreeSurfaceSgsError3, FreeSurfaceSgsReport3, FreeSurfaceSgsWorkspace3,
};
use crate::material::{MaterialLbmError, finite, positive};
use crate::sgs3::{SgsConfig, SgsModel, SgsUse};

/// Explicit exploration controls. No coefficient or physical filter is inferred.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaterialFreeSurfaceSgsConfig3 {
    /// WALE or Vreman; the native stress-only Smagorinsky path is separate.
    pub model: SgsModel,
    /// Isotropic spatial filter width [m], converted through the original dx.
    pub filter_width_m: f64,
    /// Dimensionless refusal boundary, never an eddy-viscosity clamp.
    pub maximum_eddy_to_molecular: f64,
    pub usage: SgsUse,
}

/// Instantaneous pre-step liquid diagnostics in the original physical units.
/// Powers are fill/density-weighted deviatoric shear dissipation, not heat
/// deposited into this isothermal owner or a whole-fluid energy balance.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaterialFreeSurfaceSgsReport3 {
    /// Includes zero-fill interface reconstruction carriers.
    pub wet_cells: usize,
    pub liquid_volume_m3: f64,
    pub resolved_dissipation_w: f64,
    pub sgs_dissipation_w: f64,
    pub maximum_eddy_viscosity_m2_s: f64,
    pub maximum_eddy_to_molecular: f64,
}

#[derive(Debug)]
pub enum MaterialFreeSurfaceSgsError3 {
    Material(MaterialLbmError),
    Sgs(FreeSurfaceSgsError3),
}

impl core::fmt::Display for MaterialFreeSurfaceSgsError3 {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "FS-LBM-MATERIAL-FREE-SURFACE-SGS3: {self:?}")
    }
}
impl std::error::Error for MaterialFreeSurfaceSgsError3 {}
impl From<MaterialLbmError> for MaterialFreeSurfaceSgsError3 {
    fn from(error: MaterialLbmError) -> Self {
        Self::Material(error)
    }
}
impl From<FreeSurfaceSgsError3> for MaterialFreeSurfaceSgsError3 {
    fn from(error: FreeSurfaceSgsError3) -> Self {
        Self::Sgs(error)
    }
}

impl MaterialFreeSurfaceSgsConfig3 {
    fn lattice(self, dx: f64) -> Result<SgsConfig, MaterialFreeSurfaceSgsError3> {
        positive(self.filter_width_m, "SGS filter width must be positive metres")?;
        let config = SgsConfig {
            model: self.model,
            filter_width_lu: positive(
                self.filter_width_m / dx,
                "SGS filter width is unrepresentable in lattice units",
            )?,
            maximum_eddy_to_molecular: self.maximum_eddy_to_molecular,
            usage: self.usage,
        };
        config.validate().map_err(FreeSurfaceSgsError3::from)?;
        Ok(config)
    }
}

struct SgsUnits {
    volume: f64,
    viscosity: f64,
    power: f64,
}

impl SgsUnits {
    fn new(scale: &MaterialFreeSurfaceScaling<'_>) -> Result<Self, MaterialLbmError> {
        let dx = scale.cell_size_m;
        Ok(Self {
            volume: positive(dx * dx * dx, "SGS volume unit is unrepresentable")?,
            viscosity: positive(
                dx * scale.velocity_scale,
                "SGS viscosity unit is unrepresentable",
            )?,
            // rho0*dx^5/dt^3: native power already includes density and fill.
            power: positive(
                (scale.mass_unit_kg * scale.velocity_scale)
                    * (scale.velocity_scale / scale.timestep_s),
                "SGS power unit is unrepresentable",
            )?,
        })
    }

    fn report(
        self,
        raw: FreeSurfaceSgsReport3,
    ) -> Result<MaterialFreeSurfaceSgsReport3, MaterialLbmError> {
        Ok(MaterialFreeSurfaceSgsReport3 {
            wet_cells: raw.wet_cells,
            liquid_volume_m3: scaled(
                raw.liquid_volume_lu, self.volume, "SGS liquid volume is unrepresentable",
            )?,
            resolved_dissipation_w: scaled(
                raw.resolved_dissipation_lu, self.power, "resolved liquid power is unrepresentable",
            )?,
            sgs_dissipation_w: scaled(
                raw.sgs_dissipation_lu, self.power, "modeled liquid power is unrepresentable",
            )?,
            maximum_eddy_viscosity_m2_s: scaled(
                raw.maximum_eddy_viscosity_lu, self.viscosity, "SGS viscosity is unrepresentable",
            )?,
            maximum_eddy_to_molecular: raw.maximum_eddy_to_molecular,
        })
    }
}

fn scaled(value: f64, unit: f64, reason: &'static str) -> Result<f64, MaterialLbmError> {
    let result = finite(value * unit, reason)?;
    if value != 0.0 && result == 0.0 {
        return Err(MaterialLbmError::Invalid(reason));
    }
    Ok(result)
}

fn poll(gate: &CancelGate) -> Result<(), MaterialFreeSurfaceSgsError3> {
    if gate.is_requested() {
        Err(MaterialLbmError::from(FreeSurfaceError3::Cancelled).into())
    } else {
        Ok(())
    }
}

impl MaterialFreeSurface3<'_> {
    /// Propose one source-bound moving-wall WALE/Vreman step in physical units.
    ///
    /// The original molecular material, force frame, capillarity and scaling
    /// remain authoritative. Accepted interior liquid velocities provide the
    /// gradient; no wall-law derivative or unresolved wall turbulence is added.
    /// The original fixed, fully flooded, tangential wall-field rules apply.
    ///
    /// Both wall exchange and SGS SI conversions, source-scenario speed limits,
    /// budgets and cancellation finish before returning the SAME native candidate.
    /// No thermal mass-route history or second fluid clone is allocated.
    /// The tuple contains candidate flow, actual wall exchange, and pre-step SGS.
    pub fn planar_wall_velocity_field_gradient_sgs_candidate_si(
        &self,
        motion: MaterialPlanarWallVelocityField3<'_>,
        controls: PlanarWallFieldControls3,
        gate: &CancelGate,
        config: MaterialFreeSurfaceSgsConfig3,
        kinematics: &mut FreeSurfaceKinematicsWorkspace3,
        workspace: &mut FreeSurfaceSgsWorkspace3,
    ) -> Result<
        (Self, MaterialPlanarWallFieldReport3, MaterialFreeSurfaceSgsReport3),
        MaterialFreeSurfaceSgsError3,
    > {
        poll(gate)?;
        let config = config.lattice(self.scaling.cell_size_m)?;
        let units = SgsUnits::new(&self.scaling)?;
        let (velocities, origin) =
            self.prepare_planar_wall_velocity_field_si(motion, controls, gate)?;
        let (candidate, wall, sgs) =
            self.lattice.planar_wall_velocity_field_gradient_sgs_candidate(
                PlanarWallVelocityField3 {
                    wall: motion.wall,
                    moment_origin: origin,
                    tangential_velocity: &velocities,
                },
                controls,
                gate,
                config,
                kinematics,
                workspace,
            )?;
        let wall = self.planar_wall_field_report_si(wall, motion.moment_origin_m, gate)?;
        let sgs = units.report(sgs)?;
        poll(gate)?;
        Ok((
            Self { lattice: candidate, scaling: self.scaling },
            wall,
            sgs,
        ))
    }

    /// Publish one combined physical wall/SGS step after every native and SI gate.
    /// On admitted units/support, zero coefficient preserves the original wall
    /// trajectory and wall receipts.
    pub fn step_with_planar_wall_velocity_field_and_gradient_sgs_si(
        &mut self,
        motion: MaterialPlanarWallVelocityField3<'_>,
        controls: PlanarWallFieldControls3,
        gate: &CancelGate,
        config: MaterialFreeSurfaceSgsConfig3,
        kinematics: &mut FreeSurfaceKinematicsWorkspace3,
        workspace: &mut FreeSurfaceSgsWorkspace3,
    ) -> Result<
        (MaterialPlanarWallFieldReport3, MaterialFreeSurfaceSgsReport3),
        MaterialFreeSurfaceSgsError3,
    > {
        let (candidate, wall, sgs) =
            self.planar_wall_velocity_field_gradient_sgs_candidate_si(
                motion,
                controls,
                gate,
                config,
                kinematics,
                workspace,
            )?;
        poll(gate)?;
        *self = candidate;
        Ok((wall, sgs))
    }
}
